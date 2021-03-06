% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DTE.R
\name{DTE}
\alias{DTE}
\title{IPW estimator for the Distributional Treatment Effect}
\usage{
DTE(
  y,
  d,
  x,
  ps,
  beta.lin.rep,
  ysup = NULL,
  trim = FALSE,
  trim.at = NULL,
  whs = NULL
)
}
\arguments{
\item{y}{An \eqn{n} x \eqn{1} vector of outcome of interest.}

\item{d}{An \eqn{n} x \eqn{1} vector of binary treatment adoption indicators.}

\item{x}{An \eqn{n} x \eqn{k}  matrix of covariates used in the propensity score estimation}

\item{ps}{An \eqn{n} x \eqn{1} vector of fitted propensity scores.}

\item{beta.lin.rep}{An \eqn{n} x \eqn{k}  matrix of estimates of the asymptotic linear representaion of the propensity score parameters (used to compute std. errors).}

\item{ysup}{An \eqn{l} x \eqn{1} vector of points in the support of y to compute the DTE at.
If NULL, then we set ysup to be all unique points in the support of y.}

\item{trim}{Logical argument to whether one should trim propensity scores. Deafault is FALSE.}

\item{trim.at}{Only used if trim=TRUE. If a scalar, trim all propensity score below trim.at and above 1 - trim.at.
If a  \eqn{2} x \eqn{1} vector, trim all propensity scores below trim.at[1] and all propensity scores above trim.at[2].
If NULL, trim.at is set to 1e-10.}

\item{whs}{An optional \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.}
}
\value{
A list containing the following components:
\item{dte}{The estimated DTE}
\item{dte.se}{Estimated (pointwise) std. error of the DTE.}
\item{dte.inf}{Estimated influence function of DTE estimator.}
\item{ysup}{The evaluation points of DTE.}
}
\description{
IPW estimator for the Distributional Treatment Effect
}
\references{
Sant'Anna, Pedro H. C, Song, Xiaojun, and Xu, Qi (2019), \emph{Covariate Distribution Balance via Propensity Scores},
      Working Paper <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551>.
}
