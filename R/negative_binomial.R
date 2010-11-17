#' Negative binomial likelihood
#' 
#' Original code by Richard Royall Oct 18 2000
#' For vector of counts (y), generates likelihood for E[Y_i] under
#' iid Poisson model.
#' 
#' @param y Vector of counts.
#' @param r Fixed scale parameter (optional)
#' @param lo Lower parameter bound to likelihood calculation (optional).
#' @param hi Upper parameter bound to likelihood calculation (optional).
#' @param robust Flag for calculating robust profile likelihood function. Can be used with fixed or unknown r (defaults to FALSE).
#' @param model Alternative model specification ("normal" or "poisson").
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#' @return likelihood object.
#' @keywords likelihood
#' @export
lnbinom <- function(y, r=NA, lo=NA, hi=NA, robust=F, model="", scale=T)
{
    n <- length(y)
    # Determine upper/lower bounds
    if(is.na(hi)) hi <- mean(y) + 4 * sqrt(var(y)/(n - 1))
    if(is.na(lo)) lo <- max(0.00001, mean(y) - 4 * sqrt(var(y)/(n - 1)))
    # Values for abscissa (negative binomial means)
    z <- seq(lo, hi, (hi - lo)/1000)
    
    # Normal model (Kalbfleisch and Sprott)
    if (tolower(model)=="normal") {
        llike <-  - ((n - 1)/2) * log(mean(y^2) - 2 * z * mean(y) + z^2)
    }
    # Poisson model
    else if (tolower(model)=="poisson") {
        llike <- sum(y) * log(z) - n * z
    }
    # Invalid model
    else if (model!="") {
        stop("Invalid model speficied in lnbinom")
    }
    else {
        llike <- numeric(length(z))
        for (i in 1:length(z)) {
            if (is.na(r)) {
                r <- rmax(y, z[i])
            }
            llike[i] <- sum(log(dnbinom(y, r, r/(r + z[i]))))
            
            # Robust profile LF
            if (robust==T) {
                if (!is.na(r)) {
                    r <- rmax(y, mean(y))
                } 
                llike <- llike * ((sum(y) * (r + mean(y)))/(r * var(y) * (n - 1)))
            }
        }
        
    }

    # Scale to max=1
    if (scale==T) {
        like <- exp(llike - max(llike))        
    }
    else {
        like <- exp(llike)
    }
    
    # Instantiate likelihood object
    likelihood <- list(x=z, lx=like)
    class(likelihood) <- "likelihood"
    likelihood$dist <- "Negative binomial"
    likelihood
}