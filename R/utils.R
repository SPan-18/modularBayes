#' @importFrom stats is.empty.model model.matrix model.response terms
parseFormula <- function(formula, data, intercept = TRUE, justX = FALSE) {

    # extract Y, X, and variable names for model formula and frame
    mt <- terms(formula, data = data)
    if (missing(data))
        data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$intercept <- mf$justX <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    if (!intercept) {
        attributes(mt)$intercept <- 0
    }

    # null model support
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf)
    X <- as.matrix(X)  # X matrix
    xvars <- dimnames(X)[[2L]]  # X variable names
    xobs <- dimnames(X)[[1L]]  # X observation names
    if (justX) {
        Y <- NULL
    } else {
        Y <- as.matrix(model.response(mf, "numeric"))  # Y matrix
    }

    return(list(Y, X, xvars, xobs))

}

# internal function: sample from multivariate normal distribution
#' @importFrom stats rnorm
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V), p))))
    stop("Dimension problem!")
  D <- chol(V)
  matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p))
}

# internal function: sample from normal inverse gamma distribution
#' @importFrom stats rgamma
normalIGammaSampler <- function(n.samples, mu, V, a, b){

  sigmasq <-  1/rgamma(n.samples, a, b)

  p <- length(mu)
  beta <- matrix(0, nrow = n.samples, ncol = p)

  for (i in 1:n.samples) {
    beta[i, ] <- rmvn(1, mu, sigmasq[i] * V)
  }

  Output <- as.matrix(cbind(beta, sigmasq))
  colnames(Output) <- c(paste0("beta", 1:ncol(beta)), "sigmaSq")

  return(Output)

}

# internal function: log-sum-exp function
logSumExp <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  if (length(x) == 0) {
    return(NA)
  }
  max_x <- max(x)
  if (is.infinite(max_x)) {
    return(max_x)
  }
  sum_exp <- sum(exp(x - max_x))
  return(max_x + log(sum_exp))
}