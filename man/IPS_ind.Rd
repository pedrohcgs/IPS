% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IPSind.R
\name{IPS_ind}
\alias{IPS_ind}
\title{Integrated Propensity Score estimator based on indicator weighting function}
\usage{
IPS_ind(
  d,
  x,
  Treated = FALSE,
  beta.initial = NULL,
  lin.rep = TRUE,
  whs = NULL,
  maxit = 50000
)
}
\arguments{
\item{d}{An \eqn{n} x \eqn{1} vector of binary treatment adoption indicators.}

\item{x}{An \eqn{n} x \eqn{k}  matrix of covariates to be used in the propensity score. First element must be a vector of 1's.}

\item{Treated}{Default is FALSE, which aims to achieve covariate distribution balance among treated, untreated and overall subpopulations.
If TRUE, then the estimator aims to achieve covariate distribution balance for the treated subpopulation.}

\item{beta.initial}{An optional \eqn{k} x \eqn{1} vector of initial values for the parameters to be optimized over.}

\item{lin.rep}{Logical argument to whether an estimator for the asymptotic linear representation of the IPS
parameters should be provided. Deafault is TRUE.}

\item{whs}{An optional \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.}

\item{maxit}{The maximum number of iterations. Defaults to 50000.  = FALSE). Deafault is 999 if boot = TRUE}
}
\value{
A list containing the following components:
\item{coefficients}{The estimated IPS_ind coefficients}
\item{fitted.values}{The IPS_ind fitted probabilities}
\item{linear.predictors}{The IPS_ind estimated index (X'beta)}
\item{lin.rep}{An estimator of the IPS_ind coefficients' asymptotic linear representation}
\item{converged}{An integer code. 0 indicates successful completion}
\item{x}{The model matrix (i.e. the matrix of covariates used to estimate the IPS_ind parameters)}
}
\description{
Integrated Propensity Score estimator based on indicator weighting function
}
\references{
Sant'Anna, Pedro H. C, Song, Xiaojun, and Xu, Qi (2019), \emph{Covariate Distribution Balance via Propensity Scores},
      Working Paper <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551>.
}
