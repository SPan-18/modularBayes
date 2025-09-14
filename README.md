# modularBayes: Modular Bayesian inference for a linear model

This R package implements modular Bayesian inference for linear regression when 
some predictors are themselves outputs of an earlier Bayesian analysis. Instead 
of plugging in point estimates, this approach borrows posterior samples of the 
predictors and, for each such sample, fits the current response model by 
sampling the regression coefficients from the conditional posterior 
distribution. This ensures proper propagation of uncertainty from the previous 
analysis to the current one, while avoiding unwanted feedback between the two 
models. See [Plummer, 2014](https://www.doi.org/10.1007/s11222-014-9503-z)),
[Bayarri et al. (2009)](https://www.doi.org/10.1214/09-BA404) and
[Jacob et al. (2017)](https://www.doi.org/10.48550/arXiv.1708.08719) for
references.

## Installation
For a quick installation of the development version, run the following command 
in R.

```r
# install.packages("devtools")
devtools::install_github("SPan-18/modularBayes")
```

## Usage
Once successfully installed, load the library in R.
```r
library(spStack)
```
View example code from the vignette by running `vignette("modularBayes")`.
