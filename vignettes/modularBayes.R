## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(modularBayes)

## -----------------------------------------------------------------------------
set.seed(1729)

n <- 100
beta0 <- c(5, -2, 1, 4)
sd0 <- 1
wts0 <- rep(1, n)

## -----------------------------------------------------------------------------
p <- length(beta0)
X <- cbind(rep(1, n), sapply(1:(p - 1), function(x) rnorm(n)))
mu <- X %*% beta0
y <- rnorm(n, mu, sd0 / sqrt(wts0))
dat <- cbind(y, X[, -1])
dat <- as.data.frame(dat)
names(dat) <- c("y", paste("x", 1:(p - 1), sep = ""))

## -----------------------------------------------------------------------------
head(dat)

## -----------------------------------------------------------------------------
nS <- 1000
modular.x3 <- matrix(rnorm(length(dat$x3) * nS, mean = rep(dat$x3, nS), sd = 0.1),
                     nrow = length(dat$x3), ncol = nS)

## -----------------------------------------------------------------------------
mod1 <- modularLM(y ~ x1 + x2, post.var = list(x3 = modular.x3),
                  data = dat)

## -----------------------------------------------------------------------------
mod1.summary <- credibleInterval(mod1)
print(mod1.summary)

## -----------------------------------------------------------------------------
mod2 <- modularLM(y ~ x1 * x2,
                  post.var = list(x3 = modular.x3,
                                  `x1:x3` = modular.x3 * dat$x1),
                  data = dat)

mod2.summary <- credibleInterval(mod2)
print(mod2.summary)

