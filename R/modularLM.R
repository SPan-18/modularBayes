#' Fit a modular Bayesian linear model
#' @param formula a symbolic description of the non-modular component of the
#' regression model to be fit.
#' @param data An optional data frame containing the variables in the model.
#' @param post.var A list tagged with modular predictor names, each containing
#' a matrix of posterior samples for that modular predictor. In case of
#' interaction terms involving a modular predictor, the interaction term should
#' also be included as a separate entry in this list.
#' @param wts A numeric vector of weights for the observations.
#' @param priors A list containing prior specifications for model parameters.
#' @param ... Additional arguments (not currently used).
#' @return An object of class \code{modularLM} containing the fitted model and
#' posterior samples.
#' @examples
#' set.seed(1729)
#' n <- 100
#' beta0 <- c(5, -2, 1, 4)
#' sd0 <- 1
#' wts0 <- rep(1, n)
#'
#' p <- length(beta0)
#' X <- cbind(rep(1, n), sapply(1:(p - 1), function(x) rnorm(n)))
#' mu <- X %*% beta0
#' y <- rnorm(n, mu, sd0 / sqrt(wts0))
#' dat <- cbind(y, X[, -1])
#' dat <- as.data.frame(dat)
#' names(dat) <- c("y", paste("x", 1:(p - 1), sep = ""))
#' nS <- 100
#' modular.x3 <- matrix(rnorm(length(dat$x3) * nS, mean = rep(dat$x3, nS), sd = 0.01),
#'                      nrow = length(dat$x3), ncol = nS)
#'
#' prior.list <- list(beta.Norm = list(rep(0.0, p), diag(1E4, p)),
#'                    sigma.sq.IG = c(0.01, 0.01))
#' mod1 <- modularLM(y ~ x1 + x2, post.var = list(x3 = modular.x3),
#'                   data = dat, priors = prior.list)
#' mod1.summary <- credibleInterval(mod1)
#' @export
modularLM <- function(formula, data, post.var, wts, priors, ...){

    ##### check for unused args #####
    formal.args <- names(formals(sys.function(sys.parent())))
    elip.args <- names(list(...))
    for(i in elip.args){
        if (!i %in% formal.args)
        warning("'", i, "' is not an argument")
    }

    ##### formula #####
    if(missing(formula)){
        stop("error: formula must be specified!")
    }

    if(inherits(formula, "formula")){
        holder <- parseFormula(formula, data)
        y <- holder[[1]]
        X <- as.matrix(holder[[2]])
        X.names <- holder[[3]]
    }else{
        stop("error: formula is misspecified")
    }

    p <- ncol(X)
    n <- nrow(X)

    if(missing(post.var)){
        stop("error: post.var must be specified!")
    }else{
        if(!is.list(post.var)){
            stop("error: post.var must be a list!")
        }
        if(length(post.var) < 1){
            stop("error: post.var must contain at least one element!")
        }
        if(any(sapply(post.var, function(x) !is.matrix(x)))){
            stop("error: each element of post.var must be a matrix!")
        }
    }
    q <- length(post.var)
    p1 <- p + q
    modular.names <- names(post.var)
    n.samples <- ncol(post.var[[1]])

    for(i in 1:q){
        if(ncol(post.var[[i]]) != n.samples){
            stop("error: each element of post.var must have the same number of columns!")
        }
        if(nrow(post.var[[i]]) != n){
            stop("error: each element of post.var must have the same number of rows as X!")
        }
    }

    all.names <- c(X.names, modular.names)

    if(missing(wts)){
        wts <- rep(1, n)
    }else{
        if(length(wts) != n){
            stop("The length of wts must be equal to the length of y.")
        }
    }

    # Check priors
    if(missing(priors)){
        beta.Norm <- list(rep(0.0, p1), diag(1E4, p1))
        sigma.sq.IG <- c(0.01, 0.01)
    }else{
        if(length(priors$beta.Norm) != 2){
            stop("The length of beta.Norm must be equal to two")
        }
        if(length(priors$beta.Norm[[1]]) != p1){
            stop(paste("beta.Norm[[1]] must have length equal to", p1, ".",
            sep = ""))
        }
        if(length(priors$beta.Norm[[2]]) != p1 * p1){
            stop(paste("beta.Norm[[2]] must be a "), p1, " by ", p1, " matrix.",
            sep = "")
        }
        if (length(priors$sigma.sq.IG) != 2) {
            stop("The length of sigma.sq.IG must be equal to 2.")
        }
        beta.Norm <- priors$beta.Norm
        sigma.sq.IG <- priors$sigma.sq.IG
    }

    beta.prior.mean <- beta.Norm[[1]]
    beta.prior.cov <- as.matrix(beta.Norm[[2]])
    prior.shape <- sigma.sq.IG[1]
    prior.rate <- sigma.sq.IG[2]

    V.beta.inv <- solve(beta.prior.cov)
    mu <- beta.prior.mean

    # Compute the necessary quantities
    tX1X1 <- crossprod(X / sqrt(wts))
    tX1y <- crossprod(X / wts, y)
    yty <- sum(y * y / wts)
    muTVBetaInvmu <- t(mu) %*% V.beta.inv %*% mu
    posterior.shape <- prior.shape + n/2

    # storage for modular predictors
    Z <- matrix(0, nrow = n, ncol = q)

    # storage for posterior samples
    posterior.samples <- array(0.0, dim = c(n.samples, p1 + 1))
    colnames(posterior.samples) <- c(paste0("beta", 1:p1), "sigmaSq")

    ptm <- proc.time()

    for(i in 1:n.samples){

        for(j in 1:q){
            Z[, j] <- post.var[[j]][, i]
        }
        tX1Z <- crossprod(X / wts, Z)
        tZZ <- crossprod(Z / sqrt(wts))
        tXX <- rbind(cbind(tX1X1, tX1Z), cbind(t(tX1Z), tZZ))
        tZy <- crossprod(Z / wts, y)
        tXy <- c(tX1y, tZy)

        # compute the posterior parameters
        posterior.precision <- V.beta.inv + tXX
        posterior.variance <- chol2inv(chol(posterior.precision))
        posterior.mean <- posterior.variance%*%(V.beta.inv%*%mu + tXy)
        posterior.rate <- prior.rate + 0.5*(muTVBetaInvmu + yty - t(posterior.mean) %*% posterior.precision %*% posterior.mean)

        posterior.samples[i, ] <- as.matrix(normalIGammaSampler(1, posterior.mean, posterior.variance, posterior.shape, posterior.rate))

    }

    run.time <- proc.time() - ptm

    # Extract the samples
    out <- list()
    out$y <- y
    out$X <- X
    out$modularX <- post.var
    out$X.names <- all.names
    out$wts <- wts
    out$priors <- list(beta.Norm = beta.Norm, sigma.sq.IG = sigma.sq.IG)
    out$samples <- list(beta = posterior.samples[, grep("beta", colnames(posterior.samples))],
                        sigmaSq = posterior.samples[, grep("sigmaSq", colnames(posterior.samples))])
    out$run.time <- run.time
    class(out) <- "modularLM"
    return(out)

}