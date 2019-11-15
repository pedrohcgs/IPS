###################################################################################
# IPW estimator for the Quantile Treatment Effect
#' IPW estimator for the Quantile Treatment Effect
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcome of interest.
#' @param d An \eqn{n} x \eqn{1} vector of binary treatment adoption indicators.
#' @param x An \eqn{n} x \eqn{k}  matrix of covariates used in the propensity score estimation
#' @param ps An \eqn{n} x \eqn{1} vector of fitted propensity scores.
#' @param beta.lin.rep An \eqn{n} x \eqn{k}  matrix of estimates of the asymptotic linear representaion of the propensity score parameters (used to compute std. errors).
#' @param tau An \eqn{l} x \eqn{1} vector of quantile to compute the QTE at. If NULL, then we set tau = 0.5.
#' @param bw  Bandwidth choice to compute densities (used to compute std. errors). Options are "ucv","nrd", "nrd0", "bcv", "SJ" - see bw.nrd for additional details.
#' Default choice is "nrd0".
#' @param trim Logical argument to whether one should trim propensity scores. Deafault is FALSE.
#' @param trim.at Only used if trim=TRUE. If a scalar, trim all propensity score below trim.at and above 1 - trim.at.
#'If a  \eqn{2} x \eqn{1} vector, trim all propensity scores below trim.at[1] and all propensity scores above trim.at[2].
#'If NULL, trim.at is set to 1e-10.
#' @param whs An optional \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' 
#' @return A list containing the following components:
#' \item{qte}{The estimated QTE}
#' \item{qte.se}{Estimated (pointwise) std. error of the QTE.}
#' \item{qte.inf}{Estimated influence function of QTE estimator.}
#' \item{tau}{The evaluation points of QTE.}
#'
#' @references
#'       Sant'Anna, Pedro H. C, Song, Xiaojun, and Xu, Qi (2019), \emph{Covariate Distribution Balance via Propensity Scores},
#'       Working Paper <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551>.
#' @export

QTE <- function(y, d, x, ps, beta.lin.rep, tau = 0.5, bw = "nrd0",
                trim = FALSE, trim.at = NULL,
                whs = NULL){
  # Define some underlying variables
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
  
  ps.d1.int <- sum((ps[d==1]>0.01)*(ps[d==1]<0.99))
  ps.d0.int <- sum((ps[d==0]>0.01)*(ps[d==0]<0.99))
  
  ps.d1.f = min(ps.d1.int, ps.d1 )
  ps.d0.f = min(ps.d0.int, ps.d0 )
  
  if((ps.d1.f>0) && (ps.d0.f>0)){
    # Compute pscore weights
    ps <- base::as.vector(ps)
    w1.ps <- whs * (d/ps)
    w0.ps <- whs * ((1-d) / (1-ps))
    
    if(trim){
      w1.ps <- w1.ps * ps.keep
      w0.ps <- w0.ps * ps.keep
    }
    # Normalize pscore weights
    w1.ps <- w1.ps/mean(w1.ps)
    w0.ps <- w0.ps/mean(w0.ps)
    #-----------------------------------------------------------------------------
    # Estimate QTE
    
    # First compute weighted empirical cdfs
    F1.hat <- w.ecdf(base::sort(y), w1.ps[base::order(y)] / n)
    F0.hat <- w.ecdf(base::sort(y), w0.ps[base::order(y)] / n)
    
    # Now compute quantile
    q1.hat <- base::as.numeric(cdfinv(F1.hat, probs = as.numeric(tau)))
    q0.hat <- base::as.numeric(cdfinv(F0.hat, probs = as.numeric(tau)))
    
    # QTE
    qte.hat <- q1.hat - q0.hat
    #-----------------------------------------------------------------------------
    # Compute influence function of QTE
    # First, get the score of qreg
    q.summand.Y1 <- base::matrix(w1.ps * outer(y, q1.hat, "<="), nrow = n)
    q.summand.Y0 <- base::matrix(w0.ps * outer(y, q0.hat, "<="), nrow = n)
    
    #Numerator of influence function -  without estimation effect
    q.Y1.inf1 <- q.summand.Y1 - matrix(outer(w1.ps, tau, "*"), nrow = n)
    q.Y0.inf1 <- q.summand.Y0 - matrix(outer(w0.ps, tau, "*"), nrow = n)
    
    # Estimation effects
    #ps derivative
    ps.dot.prime <- (ps * (1 - ps)) * x
    # estimate the expectations of derivatives wrt pscore parameters
    G1.beta <- base::crossprod(ps.dot.prime/ps,
                               q.Y1.inf1) / n           # k x n.tau
    G0.beta <-  base::crossprod(ps.dot.prime/(1 - ps),
                                q.Y0.inf1) / n          # k x n.tau
    
    # Kernel density estimates influence function for each tau
    
    # initialize to speed up
    q.Y1.inf <- matrix(0, nrow = n, ncol = st)
    q.Y0.inf <- matrix(0, nrow = n, ncol = st)
   
    
     # Sanity check for the bandwidth choice
    if( (bw %in% c("ucv","nrd", "nrd0", "bcv", "SJ")) == FALSE) {
      stop ("bw must be one of c('ucv','nrd', 'nrd0', 'bcv', 'SJ)")
    }
    
    # Compute kernel density
    f1 <- stats::density(y, kernel = "gaussian", bw = bw, weights = w1.ps/n)
    f0 <- stats::density(y, kernel = "gaussian", bw = bw, weights = w0.ps/n)
    
    # Evaluate Kernel density at the exact points
    f1 <- stats::approx(f1$x, f1$y, xout = q1.hat)$y
    f0 <- stats::approx(f0$x, f0$y, xout = q0.hat)$y
    
   
    for (s in 1:st){
      # Compute the influence function at each tau
      q.Y1.inf[,s] <- -(q.Y1.inf1[,s] - beta.lin.rep %*% G1.beta[,s])/f1[s]
      q.Y0.inf[,s] <- -(q.Y0.inf1[,s] + beta.lin.rep %*% G0.beta[,s])/f0[s]
    }
    
    # Influence function of the qte
    qte.inf <- q.Y1.inf - q.Y0.inf   # n x st
    #-----------------------------------------------------------------------------
    # Compute standard error
    qte.var <- base::diag(stats::cov(qte.inf))
    qte.se <- base::sqrt(qte.var/n)
  } else {
    qte.hat <- NA
    qte.se <- NA
    qte.inf <- NA
    base::warning(paste0("No observations with pscore between ",
                         max(ps.min, 0.01)," and ",
                         min(ps.max, 0.99),
                         " in either treated or comparison group." , sep=" "))
  }
  #-----------------------------------------------------------------------------
  # Return qte.hat, qte.se, qte.inf, and tau
  out <- list(qte = qte.hat,
              qte.se = qte.se,
              qte.inf = qte.inf,
              tau = tau)
  return(out)
}


