\name{IPS-package}
\alias{IPS-package}
\alias{IPS}
\docType{package}
\title{
  IPS: Covariate Distribution Balance via Integrated Propensity Scores
}
\description{
 This package implements the different integrated propensity score (IPS) estimators proposed in Sant'Anna, Song and Xu (2019), Covariate Distribution Balance via Propensity Scores<https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551>, and also the inverse probabily weigthed (IPW) estimators for the average, quantile and distributional treatment effects that build on these IPS estimators.

The IPS is estimated by fully exploiting the covariate balancing of the propensity score, i.e., by maximing the entire covariate distribution balance between the treated, untreated, and combined groups. The IPS estimators are data-driven, do not rely on tuning parameters such as bandwidths, and admit an asymptotic linear representation, which, in turn, facilitates the statistical analysis of IPW estimators for the average, quantile and distributional treatment effects.

We emphasize that the IPS can be used under different "research designs", including not only the unconfounded treatment assignment setup, but also the "local treatment effect" setup, where selection into treatment is possibly endogenous but a binary instrumental variable is available, see Sant'Anna, Song and Xu (2019)<https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551> for further details.

}

\author{
Pedro H. C. Sant'Anna, Vanderbilt University.

Xiaojun Song, Peking University.

Qi Xu, Vanderbilt University

Maintainer: Pedro H. C. Sant'Anna <pedro.h.santanna@vanderbilt.edu>
}
\references{
  Sant'Anna, Pedro H., Song, Xiaojun, and Xu, Qi. (2019), Covariate Distribution Balance via Propensity Scores <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551>, Working paper.
}
