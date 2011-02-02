#' Normal likelihood for mean
#'
#' For data (y), generates likelihood for mean under
#' iid Normal model.
#' 
#' @param y Observations.
#' @param sigma Standard deviations of normal distribution. If numeric 
#' value is given, sd is taken to be fixed; if TRUE, likelihood is to be
#' calculated for sd; if FALSE, no likelihood is calculated for sd 
#' (defaults to TRUE).
#' @param lo Lower parameter bound to likelihood calculation (optional).
#' @param hi Upper parameter bound to likelihood calculation (optional).
#' @param lpoints The number of evenly-spaced points in the interval over
#' which to calculate profile likelihood (defaults to 1000).
#' @param profile Plot Student's t profile LF (defaults to FALSE)
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE.
#' @return likelihood object.
#' @keywords likelihood normal
#' @export
lnorm.mean <- function(y, sigma=NA, lo=NA, hi=NA, lpoints=1000, profile=FALSE, 
    estimated=F, scale=T) {
    
    n <- length(y)
    if(n >= 5)	k <- 3
    else if (n == 4) k <- 4
    else if (n == 3) k <- 7
	else if (n == 2) k <- 18
    else stop("Must have more than one value of y")
		
	if(is.na(lo))
		lo <- mean(y) - k * sqrt(var(y)/n)
	if(is.na(hi))
		hi <- mean(y) + k * sqrt(var(y)/n)
		
	mu <- seq(lo, hi, length = lpoints)
	
    # Estimated likelihood with MLE for variance
    if (estimated==T) {
        like <- -(n/2) * ((mean(y) - mu)^2)/((var(y) * (n - 1))/n)
    }
    # Likelihood assuming sigma is known
    else if (!is.na(sigma)) {
        like <- -(n * (mean(y) - mu)^2)/(2 * sigma^2)
    }
    else {
        like <- -((n - 1)/2) * log(mean(y^2) - 2 * mu * mean(y) + mu^2)

        # Student t profile
        if (profile==T) {
           like <- (like * n)/(n - 1)
        }
    }
    
    if (scale==T) {
        like <- like - max(like)
    }
    
    # Return likelihood object
    likelihood <- list(x=mu, lx=exp(like))
    class(likelihood) <- "likelihood"
    likelihood$dist <- "Normal"
    likelihood
}

#' Normal likelihood for variance
#'
#' For data (y), generates likelihood for variance under
#' iid Normal model, based on marginal chi-square likelihood function
#' for sample variance.
#' 
#' @param y Observations.
#' @param mu Mean of normal distribution. If numeric value is given, mean
#' is taken to be fixed; if TRUE, likelihood is to be calculated for mean; 
#' if False, no likelihood is calculated for mean (defaults to TRUE)
#' @param lo Lower parameter bound to likelihood calculation (optional).
#' @param hi Upper parameter bound to likelihood calculation (optional).
#' @param lpoints The number of evenly-spaced points in the interval over
#' which to calculate profile likelihood (defaults to 1000).
#' @param estimated Calculate estimated likelihood by replacing
#' true mean by sample mean in full likelihood.
#' @profile Calculate profile likelihood.
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE.
#' @return likelihood object.
#' @keywords likelihood normal
#' @export

lnorm.var <- function(y, mu=NA, lo=NA, hi=NA, lpoints=1000, estimated=FALSE, 
    profile=FALSE, scale=TRUE) {
    
    n <- length(y)
    ss <- sum((y - mean(y))^2)
    if(missing(lo))
		lo <- ss/qchisq(0.999, n-1)
	if(lo == 0)
		lo <- 0.001
	if(missing(hi))
		hi <- ss/qchisq(0.001, n-1)
	sig <- seq(lo, hi, length=lpoints)
	
	# Estimated likelihood
	if (estimated==TRUE) {
	    like <- log((ss/(n*sig))^(n/2)) - (ss/(n*sig) - 1) * n/2
	}
	# Profile likelihood
	else if (profile==TRUE) {
	    like <- log((ss/((n - 1)*sig))^(n/2)) - ss/(2*sig)
	}
	# Likelihood function
	else {
	    # Assuming mu is known
    	if (!is.na(mu)) {
    	    ss <- sum((y - mu)^2)
    	    like <- (1 + log(ss/(n*sig)) - ss/(n*sig))*n/2
    	}
    	like <- (1 + log(ss/((n - 1)*sig)) - ss/((n - 1)*sig))*(n-1)/2
	}
	
	if (scale==T) {
        like <- like - max(like)
    }
	
	# Return likelihood object
    likelihood <- list(x=sig, lx=exp(like))
    class(likelihood) <- "likelihood"
    likelihood$dist <- "Normal"
    likelihood
}

# TODO Implement probability of finding weak or misleading evidence for normal models

# TODO Add likelihood for normal CV (NORM_CV.R)

# TODO Implement likelihood for difference between normal means

# TODO Implement likelihood for ratio of normal means