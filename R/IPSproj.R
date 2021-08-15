###################################################################################
# Integrated Propensity Score estimator based on projection weighting function
#' Integrated Propensity Score estimator based on projection weighting function
#'
#' @param d An \eqn{n} x \eqn{1} vector of binary treatment adoption indicators.
#' @param x An \eqn{n} x \eqn{k}  matrix of covariates to be used in the propensity score. First element must be a vector of 1's.
#' @param xbal An \eqn{n} x \eqn{l}, \eqn{l\leq k}, matrix of ``raw'' covariares to be balanced (does not need to include interaction terms). Default is \code{NULL}, which will use the same as x. 
#' @param Treated Default is FALSE, which aims to achieve covariate distribution balance among treated, untreated and overall subpopulations.
#' If TRUE, then the estimator aims to achieve covariate distribution balance for the treated subpopulation.
#' @param beta.initial An optional \eqn{k} x \eqn{1} vector of initial values for the parameters to be optimized over.
#' @param lin.rep Logical argument to whether an estimator for the asymptotic linear representation of the IPS
#' parameters should be provided. Deafault is TRUE.
#' @param whs An optional \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' @param maxit The maximum number of iterations. Defaults to 50000
#' @param allRows Default is FALSE, which attempts to fist check if all rows of the matrix x are unique. If there are draws, it tries to optimize the code, but requires 'x' to be sorted such that all unique rows are together. 
#' @param x_keep Default is FALSE. If TRUE, we return covariate matrix in the output.
#' 
#' @return A list containing the following components:
#' \item{coefficients}{The estimated IPS_proj coefficients}
#' \item{fitted.values}{The IPS_proj fitted probabilities}
#' \item{linear.predictors}{The IPS_proj estimated index (X'beta)}
#' \item{lin.rep}{An estimator of the IPS_proj coefficients' asymptotic linear representation}
#' \item{converged}{An integer code. 0 indicates successful completion}
#' \item{x}{The model matrix (i.e. the matrix of covariates used to estimate the IPS_proj parameters). Only returned if \code{x_keep = TRUE}.}
#'
#'
#' @references
#'       Sant'Anna, Pedro H. C, Song, Xiaojun, and Xu, Qi (2019), \emph{Covariate Distribution Balance via Propensity Scores},
#'       Working Paper <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551>.
#'       
#' @export

#-------------------------------------------------------------------------------
IPS_proj = function(d, x, xbal = NULL,Treated = FALSE,
                    beta.initial = NULL, lin.rep = TRUE,
                    whs = NULL,  x_keep = FALSE,
                    maxit = 50000,
                    allRows = FALSE) {
  #-----------------------------------------------------------------------------
  # Define some underlying variables
  d <-  base::as.matrix(d)
  x <- base::as.matrix(x)
  n <- base::dim(x)[1]
  k <- base::dim(x)[2]
  treated.flag <- base::as.numeric(base::isTRUE(Treated))
  if(is.null(whs)) whs <- rep(1, n)
  if(!is.numeric(whs)) base::stop("weights must be a NULL or a numeric vector")
  #-----------------------------------------------------------------------------
  # FIRST ELEMENT OF X MUST BE A CONSTANT
  if(all.equal(as.numeric(x[,1]), rep(1,n)) == FALSE) {
    stop(" first element of x must be a vector of 1's")
  }
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  #Weight function based on prjection weights
  # data_ips <- cbind(d,x)
  
  if(is.null(xbal)) {
    xbal <- x
  } else {
    xbal <- base::as.matrix(xbal)
  }
  
  # Test if all observations are unique, as this allow us to speed up the codes
  n.unique <- dplyr::n_distinct(xbal) 
  
  if( ((n - n.unique) > 500) &&  allRows == FALSE) {
    #   use code that avoid number double calculations
    x1 <- data.table::data.table(xbal)
    x1 <- data.table::data.table(x1, key = colnames(x1))
    if(max(abs(xbal-x1))>0) {
      stop("Matrix 'x' must be sorted such that all unique rows are together./n Otherwise set 'uniqueRows = T', though is usually slower.")
    }
    x_unique <- as.matrix(plyr::count(x1)[,-(dim(x1)[2]+1)])
    vec_rep  <- as.vector(plyr::count(x1)[,(dim(x1)[2]+1)])
    w.proj <- weightIPSproj_uniq(x_unique, vec_rep) 
    
  }
  else {
    w.proj <- weightIPSproj(xbal)
  }
  
  #-----------------------------------------------------------------------------
  # initial parameter value for IPS
  if (is.null(beta.initial)==TRUE){
    beta.initial <- base::suppressWarnings(CBPS::CBPS(d ~ x[,-1],
                                                      ATT = 0)$coefficients)
  }
  #-----------------------------------------------------------------------------
  # Define the Objective function for exponential weights
  #-----------------------------------------------------------------------------
  # Define the gradient of the objective function
  #-----------------------------------------------------------------------------
  # Now we are ready to estimate the pscore parameters
  ips.est.proj <- stats::optim(par = beta.initial,
                               fn = objIPS,
                               gr = gradIPS,
                               method = "BFGS",
                               control =  list(maxit = maxit, abstol = 1e-8, reltol=1e-8),
                               d = d,
                               X = x,
                               w = w.proj,
                               treated_flag = treated.flag,
                               whs = whs)
  
  beta.hat.ips <- ips.est.proj$par
  converged <- ips.est.proj$convergence
  linear.predictors <- x %*% beta.hat.ips
  ps.hat <- as.numeric(1/(1 + exp(-linear.predictors)))
  probs.min <- 1e-10
  if(base::any(ps.hat<probs.min)) {
    base::message("IPS.proj: fitted probabilities smaller than 1e-10 occurred. We truncate these.")
  }
  if(base::any(ps.hat>(1-probs.min))) {
    base::message("IPS.proj: fitted probabilities bigger than 1 - 1e-10 occurred. We truncate these.")
  }
  ps.hat <- base::pmin(1 - probs.min, ps.hat)
  ps.hat <- base::pmax(probs.min, ps.hat)
  
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Next, we compute an estimate of the asymptotic linear representation of
  # beta.hat - beta
  if (lin.rep == TRUE){
    lin.rep.hat <- linIPS(beta.hat.ips, d, ps.hat, x, w.proj, treated.flag, whs)
    covSing <- (Matrix::rankMatrix(base::crossprod(lin.rep.hat))[1] == base::dim(lin.rep.hat)[2])
    if(covSing==FALSE) base::message("IPS.proj: The variance-Covariance matrix is close to singular. Used Generalized-Inverse to compute std. errors.")
    
  }
  if(converged!=0) base::warning("IPS.proj: IPS optmization did not converge.")
  
  if(x_keep != TRUE){
    x = NULL
  }
  
  
  out <- list(coefficients = beta.hat.ips,
              fitted.values = ps.hat,
              linear.predictors = linear.predictors,
              lin.rep = lin.rep.hat,
              converged = converged,
              x = x,
              treated.flag = Treated)
  return(out)
}