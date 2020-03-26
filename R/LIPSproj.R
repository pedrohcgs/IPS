###################################################################################
# 'Local' Integrated Propensity Score estimator based on projection weighting function
#' 'Local' Integrated Propensity Score estimator based on projection weighting function
#'
#' @param z An \eqn{n} x \eqn{1} vector of binary instruments.
#' @param d An \eqn{n} x \eqn{1} vector of binary treatment adoption indicators.
#' @param x An \eqn{n} x \eqn{k}  matrix of covariates to be used in the propensity score. First element must be a vector of 1's.
#' @param beta.initial An optional \eqn{k} x \eqn{1} vector of initial values for the parameters to be optimized over.
#' @param lin.rep Logical argument to whether an estimator for the asymptotic linear representation of the LIPS
#' parameters should be provided. Deafault is TRUE.
#' @param whs An optional \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' @param maxit The maximum number of iterations. Defaults to 50000.  = FALSE). Deafault is 999 if boot = TRUE
#' @param allRows Default is FALSE, which attempts to fist check if all rows of the matrix x are unique. If there are draws, it tries to optimize the code, but requires 'x' to be sorted such that all unique rows are together. 
#' 
#' @return A list containing the following components:
#' \item{coefficients}{The estimated LIPS_proj coefficients}
#' \item{fitted.values}{The LIPS_proj fitted probabilities}
#' \item{linear.predictors}{The LIPS_proj estimated index (X'beta)}
#' \item{lin.rep}{An estimator of the LIPS_proj coefficients' asymptotic linear representation}
#' \item{converged}{An integer code. 0 indicates successful completion}
#' \item{x}{The model matrix (i.e. the matrix of covariates used to estimate the LIPS_proj parameters)}
#'
#'
#' @references
#'       Sant'Anna, Pedro H. C, Song, Xiaojun, and Xu, Qi (2019), \emph{Covariate Distribution Balance via Propensity Scores},
#'       Working Paper <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551>.
#'       
#' @export

#-------------------------------------------------------------------------------
LIPS_proj = function(z, d, x,
                    beta.initial = NULL, lin.rep = TRUE,
                    whs = NULL,
                    maxit = 50000,
                    allRows = F) {
  #-----------------------------------------------------------------------------
  # Define some underlying variables
  d <-  base::as.matrix(d)
  x <- base::as.matrix(x)
  z <- base::as.matrix(z)
  n <- base::dim(x)[1]
  k <- base::dim(x)[2]

  if(is.null(whs)) whs <- rep(1, n)
  if(!is.numeric(whs)) base::stop("weights must be a NULL or a numeric vector")
  #-----------------------------------------------------------------------------
  # FIRST ELEMENT OF X MUST BE A CONSTANT
  if(all.equal(x[,1], rep(1,n)) == F) {
    stop(" first element of x must be a vector of 1's")
  }
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  #Weight function based on prjection weights
  data_ips <- cbind(d,z,x)
  
  # Test if all observations are unique, as this allow us to speed up the codes
  
  # Test if all observations are unique, as this allow us to speed up the codes
  n.unique <- dplyr::n_distinct(x) 
  if( (n.unique != n) &&  allRows == F) {
    #   use code that avoid number double calculations
    x1 <- data.table::data.table(x)
    x1 <- data.table::data.table(x1, key = colnames(x1))
    if(max(abs(x-x1))>0) {
      stop("Matrix 'x' must be sorted such that all unique rows are together./n Otherwise set 'uniqueRows = T', though is usually slower.")
    }
    x_unique <- as.matrix(plyr::count(x1)[,-(dim(x1)[2]+1)])
    vec_rep  <- as.vector(plyr::count(x1)[,(dim(x1)[2]+1)])
    w.proj <- weightIPSproj_uniq(x_unique, vec_rep) 
    
  }
  else {
    w.proj <- weightIPSproj(x)
  }
  
  
  
  
  #-----------------------------------------------------------------------------
  # initial parameter value for LIPS
  if (is.null(beta.initial)==TRUE){
    beta.initial <- base::suppressWarnings(CBPS::CBPS(z ~ x[,-1],
                                                      ATT = 0)$coefficients)
   }
  #-----------------------------------------------------------------------------
  # Define the Objective function for exponential weights
  #-----------------------------------------------------------------------------
  # Define the gradient of the objective function
  #-----------------------------------------------------------------------------
  # Now we are ready to estimate the pscore parameters
  ips.est.proj <- stats::optim(par = beta.initial,
                               fn = objLIPS,
                               gr = gradLIPS,
                               method = "BFGS",
                               control =  list(maxit = maxit, abstol = 1e-8, reltol=1e-8),
                               d = d,
                               z = z,
                               X = x,
                               w = w.proj, 
                               whs = whs)
  
  beta.hat.ips <- ips.est.proj$par
  converged <- ips.est.proj$convergence
  linear.predictors <- x %*% beta.hat.ips
  ips.hat <- as.numeric(1/(1 + exp(-linear.predictors)))
  probs.min <- 1e-10
  if(base::any(ips.hat<probs.min)) {
    base::message("LIPS.proj: fitted probabilities smaller than 1e-10 occurred. We truncate these.")
  }
  if(base::any(ips.hat>(1-probs.min))) {
    base::message("LIPS.proj: fitted probabilities bigger than 1 - 1e-10 occurred. We truncate these.")
  }
  ips.hat <- base::pmin(1 - probs.min, ips.hat)
  ips.hat <- base::pmax(probs.min, ips.hat)
  
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Next, we compute an estimate of the asymptotic linear representation of
  # beta.hat - beta
  if (lin.rep == TRUE){
    lin.rep.hat <- linLIPS(beta.hat.ips, d, z, ips.hat, x, w.proj, whs)
    covSing <- (Matrix::rankMatrix(base::crossprod(lin.rep.hat))[1] == base::dim(lin.rep.hat)[2])
    if(covSing==F) base::message("LIPS.proj: The variance-Covariance matrix is close to singular. Used Generalized-Inverse to compute std. errors.")
  }
  
  if(converged!=0) base::warning("LIPS.proj: LIPS optmization did not converge.")
  
  
  
  out <- list(coefficients = beta.hat.ips,
              fitted.values = ips.hat,
              linear.predictors = linear.predictors,
              lin.rep = lin.rep.hat,
              converged = converged,
              x = x
              )
  return(out)
}