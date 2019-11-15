###################################################################################
# IPW estimator for the Average Treatment Effect
#' IPW estimator for the Average Treatment Effect
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcome of interest.
#' @param d An \eqn{n} x \eqn{1} vector of binary treatment adoption indicators.
#' @param x An \eqn{n} x \eqn{k}  matrix of covariates used in the propensity score estimation
#' @param ps An \eqn{n} x \eqn{1} vector of fitted propensity scores.
#' @param beta.lin.rep An \eqn{n} x \eqn{k}  matrix of estimates of the asymptotic linear representaion of the propensity score parameters (used to compute std. errors)
#' @param trim Logical argument to whether one should trim propensity scores. Deafault is FALSE.
#' @param trim.at Only used if trim=TRUE. If a scalar, trim all propensity score below trim.at and above 1 - trim.at.
#'If a  \eqn{2} x \eqn{1} vector, trim all propensity scores below trim.at[1] and all propensity scores above trim.at[2].
#'If NULL, trim.at is set to 1e-10.
#' @param whs An optional \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' 
#' @return A list containing the following components:
#' \item{ate}{The estimated ATE}
#' \item{ate.se}{Estimated std. error of the ATE.}
#' \item{ate.inf}{Estimated influence function of ATE estimator.}
#'
#' @references
#'       Sant'Anna, Pedro H. C, Song, Xiaojun, and Xu, Qi (2019), \emph{Covariate Distribution Balance via Propensity Scores},
#'       Working Paper <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551>.
#' @export


# This function computes the ATE using logistic propensity score with normalization
ATE <- function(y, d, x, ps, beta.lin.rep,
                trim = FALSE, trim.at = NULL,
                whs = NULL){
  #-----------------------------------------------------------------------------
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
  #-----------------------------------------------------------------------------
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
    #-----------------------------------------------------------------------------
    # Estimate ATE
    mu.summand.Y1 <- w1.ps * y
    mu.summand.Y0 <- w0.ps * y
    
    mu.Y1 <- base::mean(mu.summand.Y1)
    mu.Y0 <- base::mean(mu.summand.Y0)
    
    ate.hat <- mu.Y1 - mu.Y0
    #-----------------------------------------------------------------------------
    # Compute influence function of ATE
    # Estimation effects
    mu.Y1.inf1 <- mu.summand.Y1 - w1.ps %*% base::t(mu.Y1)
    mu.Y0.inf1 <- mu.summand.Y0 - w0.ps %*% base::t(mu.Y0)
    
    #ps derivative
    ps.dot.prime <- (ps * (1 - ps)) * x
    # estimate the expectations of derivatives wrt pscore parameters
    G1.beta <- base::crossprod(ps.dot.prime/ps,
                               mu.Y1.inf1) / n
    G0.beta <-  base::crossprod(ps.dot.prime/(1 - ps),
                                mu.Y0.inf1) / n
    # Estimations effect themselves
    mu.Y1.est.eff <- beta.lin.rep %*% G1.beta
    mu.Y0.est.eff <- beta.lin.rep %*% G0.beta
    
    # ate Influence function
    mu.Y1.inf <- mu.Y1.inf1 - mu.Y1.est.eff
    mu.Y0.inf <- mu.Y0.inf1 + mu.Y0.est.eff
    ate.inf <- mu.Y1.inf - mu.Y0.inf
    #-----------------------------------------------------------------------------
    # Compute standard error
    ate.var <- stats::var(ate.inf)
    ate.se <- base::sqrt(ate.var/n)
  } else {
    ate.hat <- NA
    ate.se <- NA
    ate.inf <- NA
    base::warning(paste0("No observations with pscore between ",
                         max(ps.min, 0.01)," and ",
                         min(ps.max, 0.99),
                         " in either treated or comparison group." , sep=" ")) 
    }
  #-----------------------------------------------------------------------------
  # Return ate.hat, ate.se, and ate.inf
  out <- list(ate = ate.hat,
              ate.se = ate.se,
              ate.inf = ate.inf)
  return(out)
}


