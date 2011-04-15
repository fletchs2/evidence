#' Gamma likelihood
#' 
#' Original code by Richard Royall Oct 26 2000
#' For vector of counts (y), generates likelihood for E[Y_i] under
#' iid gamma model.
#' 
#' @param y Vector of data.
#' @param beta Fixed scale parameter (optional)
#' @param lo Lower parameter bound to likelihood calculation (optional).
#' @param hi Upper parameter bound to likelihood calculation (optional).
#' @param robust Flag for calculating robust profile likelihood function. Can be used with fixed or unknown beta (defaults to FALSE).
#' @param estimated Estimated likelihood function; beta estimated by MLE (defaults to FALSE)
#' @param normal Normal model specification (defaults to FALSE).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#' @return likelihood object.
#' @keywords likelihood
#' @export

gamma_like <- function(y, beta=NA, lo=NA, hi=NA, robust=FALSE, estimated=FALSE, normal=FALSE, lpoints=1000, scale=TRUE) {
    
    n <- length(y)
    x <- log(y)
    
    # Setting limits using log-normal model
    
    mle <- exp(mean(x) + var(x)/2)
    
    if (is.na(lo)) {
        lo <- max(0.01, mle * (1 - 5 * sqrt((var(x) + var((x - mean(x))^2)/4)/n)))
    }
    if (is.na(hi)) {
        hi <- mle * (1 + 5 * sqrt((var(x) + var((x - mean(x))^2)/4)/n))
    }
    
    z <- seq(lo, hi, length=lpoints)
    
    # Normal model
    if (normal==T) {
        like <- -((n - 1)/2) * log(mean(y^2) - 2*z*mean(y) + z^2)
    }
    else {
        # Estimated LF
        if (estimated==T) {
            beta <- gamma_scale_max(y, mean(y))
        }
        # Profile LF
        else if (is.na(beta)) {

            # Standard profile likelihood
            beta <- sapply(z, function(v) gamma_scale_max(y, v))
        }
        
        # Calculate log-likelihood
        like <- n*beta*(log(beta/z) - lgamma(beta)/beta +
            ((beta - 1)*mean(log(y)))/beta - mean(y)/z)     
        
        # Robust LF
        if (robust==T) {
            if (length(beta)>1) {
                beta <- gamma_scale_max(y, mean(y))
            }
            like <- like^(mean(y)^2/(beta*var(y)))
        }       
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
   
#' Gamma scale maximizer
#'
#' For a sample y, finds value of scale parameter (beta) in gamma model with
#' fixed mean (expval) that maximizes the likelihood.
#' For use in profile likelihood in lgamma.
#' 
#' @param y Vector of data. 
#' @param expval Mean of gamma distribution.
#' @return scalar value beta
#' @keywords maximizer, support function
#' @export
gamma_scale_max <- function(y, expval) {
   
   # Initial beta value
   beta <- max(floor(expval^2/mean((y - expval)^2)), 1)
	
	# Loop over increment values
   for (i in c(-1, 1, -0.1, -0.01, 0.01)) {
      
      # Propose new beta value
      new_beta <- max(beta + i, -i)
	   repeat {
   	   log_ratio <- log(((new_beta/beta) * 
   	      dgamma((y * new_beta)/expval, new_beta))/
   	      dgamma((y * beta)/expval, beta))
   	   # Accept if log ratio is positive
   		if(sum(log_ratio, na.rm = T) > 0) {
   			beta <- new_beta
   			new_beta <- max(beta + i, -i)
   		}
   		else break
   	}
      
   }
	beta   
}