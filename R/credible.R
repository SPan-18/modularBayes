#' Generate credible intervals for model parameters
#' @param mod.out model output from running `simpleLM` or `modularLM`.
#' @param param (optional) A character vector of variable names for which
#' credible intervals are desired. If not specified, credible intervals for all
#' parameters are returned.
#' @param prob Specifies the high probability density region of the posterior
#' distribution. Default is 0.95, which corresponds to a 95% credible interval.
#' @return An object of class `credibleInterval` containing the lower and upper
#' bounds of the credible intervals for the specified parameters, along with
#' posterior medians.
#' @importFrom stats quantile
#' @export
credibleInterval <- function(mod.out, param, prob = 0.95){


    if(!inherits(mod.out, c("simpleLM", "modularLM"))){
        stop("mod.out must be an object of class 'simpleLM' or 'modularLM'")
    }

    if(missing(param)){
        param <- names(mod.out$samples)
    }else{
        if(!all(param %in% names(mod.out$samples))){
            stop("param must be either missing or one of 'beta' or 'sigmaSq'")
        }
    }

    alpha <- 1 - prob
    probs <- c(0.5, alpha / 2, 1 - alpha / 2)

    out_list <- list()

    for(p in param){

        samp <- mod.out$samples[[p]]

        if(p == "beta"){

            if(!is.matrix(samp)){
                stop("'beta' must be stored as a matrix with cols = coefficients and rows = draws")
            }
            if(ncol(samp) != length(mod.out$X.names)){
                stop("Number of columns in 'beta' samples must match number of predictors in the model")
            }

            ci <- t(apply(samp, 2, quantile, probs = probs))

            colnames(ci) <- paste0(c(50, round(100 * alpha / 2, 2), round(100 * (1 - alpha / 2), 2)), "%")
            rownames(ci) <- mod.out$X.names
            out_list[[p]] <- ci

        }else{

            if(!is.numeric(samp)){
                stop("Parameter ", p, " must be a numeric vector of draws")
            }
            ci <- quantile(samp, probs = probs)
            out_list[[p]] <- matrix(ci, nrow = 1, dimnames = list(paste0(p), names(ci)))
        }

    }

    out <- do.call("rbind", out_list)
    return(out)

}