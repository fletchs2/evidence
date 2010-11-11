#' Poisson likelihood
#' 
#' Original code by Richard Royall Oct 17 2000
#' For vector of counts (y), generates likelihood for E[Y_i] under
#' iid Poisson model.
#' 
#' @param y Vector of counts.
#' @param lo Lower parameter bound to likelihood calculation (optional).
#' @param hi Upper parameter bound to likelihood calculation (optional).
#' @param robust Flag for calculating robust correction (defaults to FALSE).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE.
#' @return likelihood object.
#' @keywords likelihood
#' @export
lpois <- function(y, lo=NA, hi=NA, robust=F, scale=T)
{

    # Set default bounds if none provided
    if(is.na(lo)) lo <- mean(y) * max(1 - 4/sqrt(sum(y)), 1/4)
    if(is.na(hi)) hi <- mean(y) * (1 + 4/sqrt(sum(y)))
    # Range of parameter
    z <- seq(lo, hi, (hi - lo)/1000)
    # Calculate likelihood over range of z
    like <- sum(y) * log(z) - length(y) * z
    # Correct for model failure
    if(robust == T) {
      like <- (like * (mean(y) * length(y)))/
          (var(y) * (length(y) - 1))
    }
    if (scale==T) {
        like <- like - max(like)
    }
    
    # Return likelihood object
    likelihood <- list(x=z, lx=exp(like))
    class(likelihood) <- "likelihood"
    likelihood$dist <- "Poisson"
    likelihood
}