###################################################################################
# IPW estimator for the Distributional Treatment Effect
#' IPW estimator for the Distributional Treatment Effect
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcome of interest.
#' @param z An \eqn{n} x \eqn{1} vector of binary instruments.
#' @param d An \eqn{n} x \eqn{1} vector of binary treatment adoption indicators.
#' @param x An \eqn{n} x \eqn{k}  matrix of covariates used in the propensity score estimation
#' @param ps An \eqn{n} x \eqn{1} vector of fitted propensity scores.
#' @param beta.lin.rep An \eqn{n} x \eqn{k}  matrix of estimates of the asymptotic linear representaion of the propensity score parameters (used to compute std. errors).
#' @param ysup An \eqn{l} x \eqn{1} vector of points in the support of y to compute the LDTE at.
#' If NULL, then we set ysup to be all unique points in the support of y.
#' @param trim Logical argument to whether one should trim propensity scores. Deafault is FALSE.
#' @param trim.at Only used if trim=TRUE. If a scalar, trim all propensity score below trim.at and above 1 - trim.at.
#'If a  \eqn{2} x \eqn{1} vector, trim all propensity scores below trim.at[1] and all propensity scores above trim.at[2].
#'If NULL, trim.at is set to 1e-10.
#' @param whs An optional \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' 
#' @return A list containing the following components:
#' \item{ldte}{The estimated LDTE}
#' \item{ldte.se}{Estimated (pointwise) std. error of the LDTE.}
#' \item{ldte.inf}{Estimated influence function of LDTE estimator.}
#' \item{ysup}{The evaluation points of LDTE.}
#'
#' @references
#'       Sant'Anna, Pedro H. C, Song, Xiaojun, and Xu, Qi (2019), \emph{Covariate Distribution Balance via Propensity Scores},
#'       Working Paper <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551>.
#' @export


LDTE <- function(y, z, d, x, ps, beta.lin.rep, ysup = NULL,
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
  
  if(is.null(whs)) whs <- rep(1, n)
  if(!is.numeric(whs)) base::stop("weights must be a NULL or a numeric vector")
  
  #-----------------------------------------------------------------------------
  # set up all unique points in the data (use it as ysup, if NULL)
  yy <- sort(unique(y))
  if (is.null(ysup) == TRUE){
    ysup <- yy
  }
  st <- length(ysup)
  ysup <- base::as.vector(ysup)
  
  #-----------------------------------------------------------------------------
  # If triming = TRUE, delete observations below threshold and above threshold
  if(is.null(trim.at)) trim.at<- 1e-10
  
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
  ps.keep <- base::as.vector((ps>ps.min)*(ps<ps.max))
  
  # Trimming message
  if(trim){
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
    w11.ps <- base::as.vector(whs * d * (z/ps))
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
    w0c <- w01.ps - w00.ps
    kappa1 <- base::mean(w1c)
    kappa0 <- base::mean(w0c)
    
    # Normalized complier weights
    w1c <- w1c/kappa1
    w0c <- w0c/kappa0
    #-----------------------------------------------------------------------------
    # Estimate DTE
    
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
    
    # Now compute rearranged CDF
    F1.hat.est <- F1.c.hat(ysup)
    F0.hat.est <- F0.c.hat(ysup)
    
    # LQTE
    ldte.hat <- F1.hat.est - F0.hat.est
    #-----------------------------------------------------------------------------
    # Compute influence function of LQTE
    # First, get the score of 1(Y(1)<= y) - n x length(FY(j)(y))
    ld.summand.Y1 <- base::matrix(w1c * outer(y, ysup, "<="), nrow = n)
    ld.summand.Y0 <- base::matrix(w0c * outer(y, ysup, "<="), nrow = n)
    
    #Numerator of influence function -  without estimation effect - n x length(ysup)
    ld.Y1.inf1 <- ld.summand.Y1 - base::matrix(outer(w1c, F1.hat.est, "*"), nrow = n)
    ld.Y0.inf1 <- ld.summand.Y0 - base::matrix(outer(w0c, F0.hat.est, "*"), nrow = n)
    
    # Estimation effects
    #ps derivative - n by k
    ps.dot.prime <- (ps * (1 - ps)) * x
    
    
    # estimate the expectations of derivatives wrt pscore parameters - 
    g1.c <- (w11.ps/ps + w10.ps/(1 - ps))/kappa1
    g0.c <- (w01.ps/ps + w00.ps/(1 - ps))/kappa0
    
    g1.c <- base::matrix(g1.c * outer(y, ysup, "<="), nrow = n) - 
      base::matrix(outer(g1.c, F1.hat.est, "*"), nrow = n)
    g0.c <- base::matrix(g0.c * outer(y, ysup, "<="), nrow = n) - 
      base::matrix(outer(g0.c, F0.hat.est, "*"), nrow = n)
    
    G1.beta <- base::crossprod(ps.dot.prime,
                               g1.c) / n   # k x n.tau
    G0.beta <-  base::crossprod(ps.dot.prime,
                                g0.c) / n  # k x n.tau
    
    
    
    # Compute the influence function at each ysup
    ld.Y1.inf <- ld.Y1.inf1 - beta.lin.rep %*% G1.beta
    ld.Y0.inf <- ld.Y0.inf1 - beta.lin.rep %*% G0.beta
    
    
    # Influence function of the dte
    ldte.inf <- ld.Y1.inf - ld.Y0.inf   # n x st
    #-----------------------------------------------------------------------------
    # Compute standard error
    ldte.var <- base::diag(stats::cov(ldte.inf))
    ldte.se <- base::sqrt(ldte.var/n)
  } else {
    ldte.hat <- NA
    ldte.se <- NA
    ldte.inf <- NA
    base::warning(paste0("No observations with pscore between ",
                         max(ps.min, 0.01)," and ",
                         min(ps.max, 0.99),
                         " in either treated or comparison group." , sep=" "))}
  #-----------------------------------------------------------------------------
  # Return dte.hat, dte.se, dte.inf, and ysup
  out <- list(ldte = ldte.hat,
              ldte.se = ldte.se,
              ldte.inf = ldte.inf,
              ysup = ysup)
  return(out)
}

