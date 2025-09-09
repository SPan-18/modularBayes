#' Fit a Bayesian linear model
#' @param formula a symbolic description of the regression model to be fit.
#' @param data An optional data frame containing the variables in the model.
#' @param wts A numeric vector of weights for the observations.
#' @param priors A list containing prior specifications for model parameters.
#' @param n.samples Number of posterior samples to draw.
#' @param ... Additional arguments (not currently used).
#' @return An object of class \code{fitLM} containing the fitted model and
#' posterior samples.
#' @export
simpleLM <- function(formula, data, wts, priors, n.samples = 1000, ...){

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

    if(missing(wts)){
        wts <- rep(1, n)
    }else{
        if(length(wts) != n){
            stop("The length of wts must be equal to the length of y.")
        }
    }

    # Check priors
    if(missing(priors)){
        beta.Norm <- list(rep(0.0, p), diag(1000.0, p))
        sigma.sq.IG <- c(0.01, 0.01)
    }else{
        if(length(priors$beta.Norm) != 2){
            stop("The length of beta.Norm must be equal to two")
        }
        if(length(priors$beta.Norm[[1]]) != p){
            stop(paste("beta.Norm[[1]] must have length equal to", p, ".",
            sep = ""))
        }
        if(length(priors$beta.Norm[[2]]) != p * p){
            stop(paste("beta.Norm[[2]] must be a "), p, " by ", p, " matrix.",
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

    tXX <- crossprod(X / sqrt(wts))
    tXy <- crossprod(X / wts, y)

    yty <- sum(y * y / wts)

    muTVBetaInvmu <- t(mu) %*% V.beta.inv %*% mu

    posterior.precision <- V.beta.inv + tXX
    posterior.variance <- chol2inv(chol(posterior.precision))
    posterior.mean <- posterior.variance%*%(V.beta.inv%*%mu + tXy)

    posterior.shape <- prior.shape + n/2
    posterior.rate <- prior.rate + 0.5*(muTVBetaInvmu + yty - t(posterior.mean) %*% posterior.precision %*% posterior.mean)

    ptm <- proc.time()

    posterior.samples <- as.matrix(normalIGammaSampler(n.samples, posterior.mean, posterior.variance, posterior.shape, posterior.rate))

    run.time <- proc.time() - ptm

    # Extract the samples
    out <- list()
    out$y <- y
    out$X <- X
    out$X.names <- X.names
    out$wts <- wts
    out$priors <- list(beta.Norm = beta.Norm, sigma.sq.IG = sigma.sq.IG)
    out$samples <- list(beta = posterior.samples[, grep("beta", colnames(posterior.samples))],
                        sigmaSq = posterior.samples[, grep("sigmaSq", colnames(posterior.samples))])
    out$run.time <- run.time
    class(out) <- "simpleLM"
    return(out)

}
