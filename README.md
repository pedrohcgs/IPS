# IPS: Covariate Distribution Balance via Integrated Propensity Scores

## Overview 


This `R` package implements the different integrated propensity score (IPS) estimators proposed in Sant'Anna, Song and Xu (2019), [Covariate Distribution Balance via Propensity Scores](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551), and also the inverse probabily weigthed (IPW) estimators for the average, quantile and distributional treatment effects that build on these IPS estimators.

The IPS is estimated by fully exploiting the covariate balancing of the propensity score, i.e., by maximing the entire covariate distribution balance between the treated, untreated, and combined groups. The IPS estimators are data-driven, do not rely on tuning parameters such as bandwidths, and admit an asymptotic linear representation, which, in turn, facilitates the statistical analysis of IPW estimators for the average, quantile and distributional treatment effects.

We emphasize that the IPS can be used under different "research designs", including not only the unconfounded treatment assignment setup, but also the "local treatment effect" setup, where selection into treatment is possibly endogenous but a binary instrumental variable is available, see [Sant'Anna, Song and Xu (2019)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551) for further details.


At the moment, The `IPS` package implements three IPS estimators and three local IPS (LIPS) estimators, the latter aiming to balancing covariate distribution among compliers: 

**IPS ESTIMATORS**
        
* `IPS_exp` - This implements the IPS estimator with the exponential weigthing function.

* `IPS_proj` - This implements the IPS estimator with the projection weigthing function.

* `IPS_ind` - This implements the IPS estimator with the indicator weigthing function --- we do not recommend using this estimator when the number of covariates is moderate or large.

**LOCAL IPS ESTIMATORS (suitable for setups with treatment noncompliance)**

      
* `LIPS_exp` - This implements the local IPS estimator with the exponential weigthing function.

* `LIPS_proj` - This implements the local IPS estimator with the projection weigthing function.

* `LIPS_ind` - This implements the local IPS estimator with the indicator weigthing function --- we do not recommend using this estimator, but included it here for completeness and transparency.


On top of the aforementioned propensity score estimators, the `IPS` package also implements IPW estimators for the average, distributional and quantile treatment effects: Check out the commands `ATE`, `ATT`, `QTE`, `QTT`, `DTE`, `DTT` for treatment effect measures under unconfoundedness, and `LATE`, `LQTE`, and `LDTE` for treatment effect measures under the local treatment effect setup.


For further details, please see the paper Sant'Anna, Song and Xu (2019), Sant'Anna, Song and Xu (2019), [Covariate Distribution Balance via Propensity Scores](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551). ***This is still a work in progress, so in case you have any comments and/or questions, please contact Pedro Sant'Anna (see email below)***.

## Installing IPS
This github website hosts the source code, and it always has the most updated version of the package.

To install the most recent version of the `IPS` package from GitHub (this is what we recommend):

        library(devtools)
        devtools::install_github("pedrohcgs/IPS")

 If you are a macOS user and are facing issues installing our package, make sure you have Xcode installed in your machine. [Here is a detailed guidelines on how to compile Rcpp codes in macOS](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/).
## Authors 

Pedro H. C. Sant'Anna, Vanderbilt University, Nashville, TN. E-mail: pedro.h.santanna [at] vanderbilt [dot] edu.

Xiaojun Song, Peking University, Beijing, China. E-mail: sxj [at] gsm [dot] pku [dot] edu [dot] cn.

Qi Xu, Vanderbilt University, Nashville, TN. E-mail: qi.xu.1 [at] vanderbilt [dot] edu.


## References

* Sant'Anna, Pedro H., Song, Xiaojun, and Xu, Qi. (2019), [Covariate Distribution Balance via Propensity Scores](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3258551), Working paper.

