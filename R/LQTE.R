###################################################################################
# IPW estimator for the Local Quantile Treatment Effect
#' IPW estimator for the Local Quantile Treatment Effect
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcome of interest.
#' @param z An \eqn{n} x \eqn{1} vector of binary instruments.
#' @param d An \eqn{n} x \eqn{1} vector of binary treatment adoption indicators.
#' @param x An \eqn{n} x \eqn{k}  matrix of covariates used in the propensity score estimation
#' @param ps An \eqn{n} x \eqn{1} vector of fitted propensity scores.
#' @param beta.lin.rep An \eqn{n} x \eqn{k}  matrix of estimates of the asymptotic linear representaion of the propensity score parameters (used to compute std. errors).
#' @param tau An \eqn{l} x \eqn{1} vector of quantile to compute the LQTE at. If NULL, then we set tau = 0.5.
#' @param bw  Bandwidth choice to compute densities (used to compute std. errors). Options are "ucv","nrd", "nrd0", "bcv", "SJ" - see bw.nrd for additional details.
#' Default choice is "nrd0".
#' @param trim Logical argument to whether one should trim propensity scores. Deafault is FALSE.
#' @param trim.at Only used if trim=TRUE. If a scalar, trim all propensity score below trim.at and above 1 - trim.at.
#'If a  \eqn{2} x \eqn{1} vector, trim all propensity scores below trim.at[1] and all propensity scores above trim.at[2].
#'If NULL, trim.at is set to 1e-10.
#' @param whs An optional \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' 
#' @return A list containing the following components:
#' \item{lqte}{The estimated LQTE}
#' \item{lqte.se}{Estimated (pointwise) std. error of the LQTE.}
#' \item{lqte.inf}{Estimated influence function of LQTE estimator.}
#' \item{tau}{The evaluation points of LQTE.}
#'
#' @references
#'       Sant'Anna, Pedro H. C, Song, Xiaojun, and Xu, Qi (2019), \emph{Covariate Distribution Balance via Propensity Scores},
#'       Working Paper <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551>.
#' @export

LQTE <- function(y, z, d, x, ps, beta.lin.rep, tau = 0.5, bw = "nrd0",
                 trim = FALSE, trim.at = NULL,
                 whs = NULL){
  # Define some underlying variables
  z <-  base::as.vector(z)
  d <-  base::as.vector(d)
  x <- base::as.matrix(x)
  ps <- base::as.vector(ps)
  beta.lin.rep <- base::as.matrix(beta.lin.rep)
  n <- base::dim(x)[1]
  k <- base::dim(x)[2]
  st <- length(tau)
  tau <- as.matrix(tau)
  if(is.null(whs)) whs <- rep(1, n)
  if(!is.numeric(whs)) base::stop("weights must be a NULL or a numeric vector")
  #-----------------------------------------------------------------------------
  # If triming = TRUE, delete observations below threshold and above threshold
  if(is.null(trim.at)) trim.at <- 1e-10
  len.trim <- length(trim.at)
  if(len.trim==1){
    ps.min <- trim.at
    ps.max <- 1 - trim.at
  }
  if(len.trim==2){
    ps.min <- trim.at[1]
    ps.max <- trim.at[2]
  }
  if(len.trim>2){
    base::stop("trim.at must be a scalar or a vector with two elements")
  }
  ps.keep <- base::as.vector((ps>ps.min) * (ps<ps.max))
  
  if(trim){
    # Trimming message
    if(base::any(ps < ps.min)) {
      base::warning(paste0("Fitted propensity scores smaller than ", ps.min," were provided. We trimmed them", sep=" "))
    }
    if(base::any(ps > ps.max)) {
      base::warning(paste0("Fitted propensity scores bigger than ", ps.max," were provided. We trimmed them", sep=" "))
    }
  }
  
  #-----------------------------------------------------------------------------
  ps.d1 <- sum(ps.keep[d==1])
  ps.d0 <- sum(ps.keep[d==0])
  if((ps.d1<20) || (ps.d0<20)){
    base::warning(paste0("Less than 20 observations with pscore between ", ps.min," and ", ps.max,
                         " in either treated or comparison group. Proceed with caution!" , sep=" "))
  }
  
  ps.d1.int <- sum( (ps[d==1]>0.01)*(ps[d==1]<0.99))
  ps.d0.int <- sum( (ps[d==0]>0.01)*(ps[d==0]<0.99))
  
  ps.d1.f = min(ps.d1.int, ps.d1 )
  ps.d0.f = min(ps.d0.int, ps.d0 )
  
  if((ps.d1.f>0) && (ps.d0.f>0)){
    # Compute instrument pscore weights
    ps <- base::as.vector(ps)
    
    # First subindex is for d, second for z
    w11.ps <-  base::as.vector(whs * d * (z/ps))
    w10.ps <-  base::as.vector(whs * d * ((1 - z)/(1 - ps)))
    w01.ps <-  base::as.vector(whs * (1 - d) * (z/ps))
    w00.ps <-  base::as.vector(whs * (1 - d) * ((1 - z) / (1 - ps)))
    
    if(trim){
      w11.ps <- w11.ps * ps.keep
      w10.ps <- w10.ps * ps.keep
      w01.ps <- w01.ps * ps.keep
      w00.ps <- w00.ps * ps.keep
    }
    
    # Complier weights
    w1c <- w11.ps - w10.ps
    w0c <- (w01.ps - w00.ps)
    kappa1 <- base::mean(w1c)
    kappa0 <- base::mean(w0c)
    
    # Normalized complier weights
    w1c <- w1c/kappa1
    w0c <- w0c/kappa0
    #-----------------------------------------------------------------------------
    # Estimate QTE
    
    # First compute weighted empirical cdfs
    y.sorted <-  base::as.vector(base::sort(y))
    y.order <-  base::as.vector(base::order(y))
    F1.c.hat <- w.ecdf(y.sorted,  base::as.vector(w1c[y.order] / n))
    F0.c.hat <- w.ecdf(y.sorted,  base::as.vector(w0c[y.order] / n))
    
    
    #Rearrange CDF Y(1) for compliers (guarantee it is non-decreasing)
    F1.c.hat <- Rearrangement::rearrangement(data.frame(y.sorted), F1.c.hat(y.sorted))
    F1.c.hat[F1.c.hat > 1] <- 1
    F1.c.hat[F1.c.hat < 0] <- 0
    F1.c.hat <- r.ecdf(y.sorted, F1.c.hat)
    
    #rearrange CDF Y(0) for compliers
    F0.c.hat <- Rearrangement::rearrangement(data.frame(y.sorted), F0.c.hat(y.sorted))
    F0.c.hat[F0.c.hat > 1] <- 1
    F0.c.hat[F0.c.hat < 0] <- 0
    F0.c.hat <- r.ecdf(y.sorted, F0.c.hat)
    
    # Now compute quantile
    lq1.hat <- base::as.numeric(cdfinv(F1.c.hat, probs = as.numeric(tau)))
    lq0.hat <- base::as.numeric(cdfinv(F0.c.hat, probs = as.numeric(tau)))
    
    # LQTE
    lqte.hat <- lq1.hat - lq0.hat
    #-----------------------------------------------------------------------------
    # Compute influence function of LQTE
    # First, get the score of lqreg - n x length(tau)
    lq.summand.Y1 <- base::matrix(w1c * outer(y, lq1.hat, "<="), nrow = n)
    lq.summand.Y0 <- base::matrix(w0c * outer(y, lq0.hat, "<="), nrow = n)
    
    #Numerator of influence function -  without estimation effect - n x length(tau)
    lq.Y1.inf1 <- lq.summand.Y1 - base::matrix(outer(w1c, tau, "*"), nrow = n)
    lq.Y0.inf1 <- lq.summand.Y0 - base::matrix(outer(w0c, tau, "*"), nrow = n)
    
    # Estimation effects
    #ps derivative - n by k
    ps.dot.prime <- (ps * (1 - ps)) * x
    
    
    # estimate the expectations of derivatives wrt pscore parameters - 
    g1.c <-  (w11.ps/ps + w10.ps/(1 - ps))/kappa1
    g0.c <-  (w01.ps/ps + w00.ps/(1 - ps))/kappa0
    
    g1.c <- base::matrix(g1.c * outer(y, lq1.hat, "<="), nrow = n) - 
      base::matrix(outer(g1.c, tau, "*"), nrow = n)
    g0.c <- base::matrix(g0.c * outer(y, lq0.hat, "<="), nrow = n) - 
      base::matrix(outer(g0.c, tau, "*"), nrow = n)
    
    G1.beta <- base::crossprod(ps.dot.prime,
                               g1.c) / n   # k x n.tau
    G0.beta <-  base::crossprod(ps.dot.prime,
                                g0.c) / n  # k x n.tau
    
    # Kernel density estimates influence function for each tau
    
    # initialize to speed up
    lq.Y1.inf <- matrix(0, nrow = n, ncol = st)
    lq.Y0.inf <- matrix(0, nrow = n, ncol = st)
    
    # Sanity check for the bandwidth choice
    if( (bw %in% c("ucv","nrd", "nrd0", "bcv", "SJ")) == FALSE) {
      stop ("bw must be one of c('ucv','nrd', 'nrd0', 'bcv', 'SJ)")
    }
    

    # Compute kernel density
    f1.c <- complier.kde.gaussian(y, lq1.hat, whs = w1c, gaussian = T)
    f1.c <- base::pmax(f1.c, .1/n)
    f0.c <-  complier.kde.gaussian(y, lq0.hat, whs = w0c, gaussian = T)
    f0.c <- base::pmax(f0.c, .1/n)

    # Compute the influence function at each tau
    for (s in 1:st){
      # Compute the influence function at each tau
      lq.Y1.inf[,s] <- -(lq.Y1.inf1[,s] - beta.lin.rep %*% G1.beta[,s])/f1.c[s]
      lq.Y0.inf[,s] <- -(lq.Y0.inf1[,s] - beta.lin.rep %*% G0.beta[,s])/f0.c[s]
    }
    
    # Influence function of the qte
    lqte.inf <- lq.Y1.inf - lq.Y0.inf   # n x st
    #-----------------------------------------------------------------------------
    # Compute standard error
    lqte.var <- base::diag(stats::cov(lqte.inf))
    lqte.se <- base::sqrt(lqte.var/n)
  } else {
    lqte.hat <- NA
    lqte.se <- NA
    lqte.inf <- NA
    base::warning(paste0("No observations with pscore between ",
                         max(ps.min, 0.01)," and ",
                         min(ps.max, 0.99),
                         " in either treated or comparison group." , sep=" ")) }
  #-----------------------------------------------------------------------------
  # Return qte.hat, qte.se, qte.inf, and tau
  out <- list(lqte = lqte.hat,
              lqte.se = lqte.se,
              lqte.inf = lqte.inf,
              tau = tau)
  return(out)
}


