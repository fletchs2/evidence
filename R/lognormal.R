#' Profile likelihood for log-normal expected value
#'
#' For data (y), generates likelihood for mean under
#' iid log-normal model.
#' 
#' If x = ln(y) has N(meanlog, varlog) then y has
#' log-normal distribution with mean exp(meanlog + 0.5*varlog).
#' 
#' @param y Observations.
#' @param varlog True variance of ln(y), if known.
#' @param lo Lower parameter bound to likelihood calculation (optional).
#' @param hi Upper parameter bound to likelihood calculation (optional).
#' @param estimated Calculate estimated LF (defaults to FALSE)
#' @param normal Calculate profile LF under normal model (defaults to FALSE)
#' @param lpoints The number of evenly-spaced points in the interval over
#' which to calculate profile likelihood (defaults to 1000).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE.
#' @return likelihood object.
#' @keywords likelihood lognormal
#' @export

llnorm <- function(y, varlog=NA, lo=NA, hi=NA, estimated=FALSE, 
    normal=FALSE, lpoints=1000, scale=TRUE) {
    
    n <- length(y)
    x <- log(y)
    
    # MLE for mean
    mle <- exp(mean(x) + var(x)/2)
    
    if (is.na(lo)) {
        lo <- max(0.01, mle * (1 - 5 * sqrt((var(x) + var((x - mean(x))^2)/4)/n)))
    }
    if (is.na(hi)) {
        hi <- mle * (1 + 5 * sqrt((var(x) + var((x - mean(x))^2)/4)/n))
    }
    
    z <- seq(lo, hi, length=lpoints)
    
    if (normal==T) {
        
        like <- -((n - 1)/2)*log(mean(y^2) - 2 * z*mean(y) + z^2)
    }
    else {
        
        if (!is.na(varlog)) {
            
            if (estimated==T) {
                print("Cannot specify both varlog and estimated LF. Calculating LF with varlog.")
            }
            sigsq <- varlog
        }
        else if (estimated==T) {
            # MLE of variance
            sigsq <- (var(x)*(n - 1))/n
        }
        else {
            c1 <- mean(x^2) - 2*mean(x)*log(z) + (log(z))^2
            sigsq <- 2*(sqrt(1 + c1) - 1)
        }
        
        like <- -(n/2)*log(sigsq) - (sum(x^2) - 2*sum(x)*(log(z) - sigsq/2) + n*(log(z) - sigsq/2)^2)/(2 * sigsq)
            
    }

    if (scale==T) {
        like <- like - max(like)
    }
    
    # Return likelihood object
    likelihood <- list(x=z, lx=exp(like))
    class(likelihood) <- "likelihood"
    likelihood$name <- "Mean"
    likelihood
}