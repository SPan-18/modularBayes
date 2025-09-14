# modularBayes: Modular Bayesian inference for a linear model

Implements modular Bayesian inference for linear regression when some predictors
are themselves outputs of an earlier Bayesian analysis. Instead of plugging in
point estimates, this approach borrows posterior samples of the predictors and,
for each such sample, fits the current response model by sampling the regression
coefficients from the conditional posterior distribution. This ensures proper
propagation of uncertainty from the previous analysis to the current one, while
avoiding unwanted feedback between the two models. The method can be viewed as a
cut model ([Plummer, 2014](https://www.doi.org/10.1007/s11222-014-9503-z)),
which severs feedback in the underlying graphical model. Such modular strategies
are useful for post-hoc analyses or combining evidence from multiple data
sources; see [Bayarri et al. (2009)](https://www.doi.org/10.1214/09-BA404) and
[Jacob et al. (2017)](https://www.doi.org/10.48550/arXiv.1708.08719).