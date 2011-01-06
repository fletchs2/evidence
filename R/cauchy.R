#' Cauchy likelihood
#' 
#' Original code by Richard Royall March 29 1996
#' For vector of counts (y), profile likelihood function for 
#  location parameter in Cauchy model (unknown scale).
#' 
#' @param y Vector of counts.
#' @param lo Lower parameter bound to likelihood calculation (optional).
#' @param hi Upper parameter bound to likelihood calculation (optional).
#' @param lpoints The number of evenly-spaced points in the interval over which to calculate profile likelihood (defaults to 1000).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE.
#' @return likelihood object.
#' @keywords likelihood
#' @export
lcauchy <- function(y, lo=NA, hi=NA, lpoints=1000, scale=T)
{

    # Set default bounds if none provided
    halfwidth <- max(abs(mean(y) - y)) + sqrt(var(y))
    if(is.na(lo)) lo <- mean(y) - halfwidth
    if(is.na(hi)) hi <- mean(y) + halfwidth
    # Range of parameter
    z <- seq(lo, hi, (hi - lo)/lpoints)
    
    llike <- numeric(length(z))
	for(i in 1:length(z)) {
	    
		scale <- quantile(abs(y - z[i]), 0.5)
		step <- scale/10
		
		llike.i <- - sum(log(1 + ((y - z[i])/scale)^2) + log(scale))
		done <- FALSE
		while(done==FALSE) {
			newscale <- scale + step
			new.like <-  - sum(log(1 + ((y - z[i])/newscale)^2) + log(newscale))
			if (new.like < llike.i) done <- TRUE
			llike.i <- new.like
			scale <- newscale
		}
		done <- FALSE
		while(done==FALSE) {
			newscale <- scale - step
			new.like <-  - sum(log(1 + ((y - z[i])/newscale)^2) + log(newscale))
			if (new.like < llike.i) done <- TRUE
			llike.i <- new.like
			scale <- newscale
		}
		step <- (step * 2)/100
		done <- FALSE
		while(done==FALSE) {
			newscale <- scale + step
			new.like <-  - sum(log(1 + ((y - z[i])/newscale)^2) + log(newscale))
			if (new.like > llike.i) {
    			llike.i <- new.like
    			scale <- newscale
            }
            else done <- TRUE
		}
		llike[i] <- llike.i
	}
	
	if (scale==TRUE) {
	    like <- exp(llike - max(llike))
	}
	else {
	    like <- exp(llike)
	}
    
    # Return likelihood object
    likelihood <- list(x=z, lx=like)
    class(likelihood) <- "likelihood"
    likelihood$dist <- "Cauchy"
    likelihood
}

# TODO Add function for calculating probability of misleading or weak evidence