###################################################################################
# IPW estimator for the Distributional Treatment Effect on the Treated
#' IPW estimator for the Distributional Treatment Effect on the Treated
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcome of interest.
#' @param d An \eqn{n} x \eqn{1} vector of binary treatment adoption indicators.
#' @param x An \eqn{n} x \eqn{k}  matrix of covariates used in the propensity score estimation
#' @param ps An \eqn{n} x \eqn{1} vector of fitted propensity scores.
#' @param beta.lin.rep An \eqn{n} x \eqn{k}  matrix of estimates of the asymptotic linear representaion of the propensity score parameters (used to compute std. errors).
#' @param ysup An \eqn{l} x \eqn{1} vector of points in the support of y to compute the DTT at.
#' If NULL, then we set ysup to be all unique points in the support of y.
#' @param trim Logical argument to whether one should trim propensity scores. Deafault is FALSE.
#' @param trim.at Only used if trim=TRUE. If a scalar, trim all propensity score below trim.at and above 1 - trim.at.
#'If a  \eqn{2} x \eqn{1} vector, trim all propensity scores below trim.at[1] and all propensity scores above trim.at[2].
#'If NULL, trim.at is set to 1e-10.
#' @param whs An optional \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' 
#' @return A list containing the following components:
#' \item{dtt}{The estimated DTT}
#' \item{dtt.se}{Estimated (pointwise) std. error of the DTT.}
#' \item{dtt.inf}{Estimated influence function of DTT estimator.}
#' \item{ysup}{The evaluation points of DTT.}
#'
#' @references
#'       Sant'Anna, Pedro H. C, Song, Xiaojun, and Xu, Qi (2019), \emph{Covariate Distribution Balance via Propensity Scores},
#'       Working Paper <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551>.
#' @export

DTT <- function(y, d, x, ps, beta.lin.rep, ysup = NULL,
                trim = FALSE, trim.at = NULL,
                whs = NULL){
  # Define some underlying variables
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
  if(is.null(trim.at)) trim.at <- c(0, 1 - 1e-10)
  
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
  
  ps.d1.int <- sum((ps[d==1]>0)*(ps[d==1]<0.99))
  ps.d0.int <- sum((ps[d==0]>0)*(ps[d==0]<0.99))
  
  ps.d1.f = min(ps.d1.int, ps.d1 )
  ps.d0.f = min(ps.d0.int, ps.d0 )
  
  if((ps.d1.f>0) && (ps.d0.f>0)){
    # Compute pscore weights
    ps <- base::as.vector(ps)
    w1.ps <- whs * d
    w0.ps <- whs * ((1-d) * ps / (1-ps))
    
    if(trim){
      w1.ps <- w1.ps * ps.keep
      w0.ps <- w0.ps * ps.keep
    }
    # Normalize pscore weights
    w1.ps <- w1.ps/mean(w1.ps)
    w0.ps <- w0.ps/mean(w0.ps)
    #-----------------------------------------------------------------------------
    # Estimate DTT
    # First compute weighted empirical cdfs
    F1.hat <- w.ecdf(base::sort(y), w1.ps[base::order(y)] / n)
    F0.hat <- w.ecdf(base::sort(y), w0.ps[base::order(y)] / n)
    
    # DTT
    F1.hat.est <- F1.hat(ysup)
    F0.hat.est <- F0.hat(ysup)
    
    dtt.hat <- F1.hat.est - F0.hat.est
    #-----------------------------------------------------------------------------
    # Compute influence function of DTT
    # First, get the score of CDF
    cdf.summand.Y1 <- base::matrix(w1.ps * outer(y, ysup, "<="), nrow = n)
    cdf.summand.Y0 <- base::matrix(w0.ps * outer(y, ysup, "<="), nrow = n)
    
    # Estimation effects
    cdf.Y1.inf1 <- cdf.summand.Y1 - w1.ps %*% base::t(F1.hat.est)
    cdf.Y0.inf1 <- cdf.summand.Y0 - w0.ps %*% base::t(F0.hat.est)
    
    # Estimation effects
    #ps derivative
    ps.dot.prime <- (ps * (1 - ps)) * x
    # estimate the expectations of derivatives wrt pscore parameters
    
    G0.beta <-  base::crossprod(ps.dot.prime/(ps * (1 - ps)),
                                cdf.Y0.inf1) / n          # k x n.ysup
    
    # Estimations effect themselves
    cdf.Y0.est.eff <- beta.lin.rep %*% G0.beta
    
    # dtt Influence function
    cdf.Y1.inf <- cdf.Y1.inf1
    cdf.Y0.inf <- cdf.Y0.inf1 + cdf.Y0.est.eff
    dtt.inf <- cdf.Y1.inf - cdf.Y0.inf
    #-----------------------------------------------------------------------------
    # Compute standard error
    dtt.var <- base::diag(stats::cov(dtt.inf))
    dtt.se <- base::sqrt(dtt.var/n)
  } else {
    dte.hat <- NA
    dte.se <- NA
    dte.inf <- NA
    base::warning(paste0("No observations with pscore between ",
                         max(ps.min, 0)," and ",
                         min(ps.max, 0.99),
                         " in either treated or comparison group." , sep=" "))
  }
  #-----------------------------------------------------------------------------
  # Return dtt.hat, dtt.se, dtt.inf, and ysup
  out <- list(dtt = dtt.hat,
              dtt.se = dtt.se,
              dtt.inf = dtt.inf,
              ysup = ysup)
  return(out)
}


