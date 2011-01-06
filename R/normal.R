#' Normal likelihood for mean
#'
#' For data (y), generates likelihood for mean and/or variance under
#' iid Normal model.
#' 
#' @param y Observations.
#' @param mu Mean of normal distribution. If numeric value is given, mean
#' is taken to be fixed; if TRUE, likelihood is to be calculated for mean; 
#' if False, no likelihood is calculated for mean (defaults to TRUE)
#' @param sigma Standard deviations of normal distribution. If numeric 
#' value is given, sd is taken to be fixed; if TRUE, likelihood is to be
#' calculated for sd; if FALSE, no likelihood is calculated for sd (defaults to TRUE)
#' @param lo Lower parameter bound to likelihood calculation (optional).
#' @param hi Upper parameter bound to likelihood calculation (optional).
#' @param lpoints The number of evenly-spaced points in the interval over
#' which to calculate profile likelihood (defaults to 1000).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE.
#' @return likelihood object.
#' @keywords likelihood
#' @export
lnorm.mean <- function(y, sigma=NA, lo=NA, hi=NA, lpoints=1000, t.prof=F, estimated=F, scale=T) {
    
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
		
	mu <- seq(lo, hi, length = m)
	
    # Estimated likelihood with MLE for variance
    if (estimated==T) {
        like <- -(n/2) * ((mean(y) - mu)^2)/((var(y) * (n - 1))/n)
    }
    else if (!is.na(sigma)) {
        like <- -(n * (mean(y) - mu)^2)/(2 * sigma^2)
    }
    else {
        like <- -((n - 1)/2) * log(mean(y^2) - 2 * mu * mean(y) + mu^2)

        # Student t profile
        if (t.prof==T) {
           like <- (like * n)/(n - 1)
        }
    }
    
    if (scale==T) {
        like <- like - max(like)
    }
    
    # Return likelihood object
    likelihood <- list(x=z, lx=exp(like))
    class(likelihood) <- "likelihood"
    likelihood$dist <- "Normal"
    likelihood
}

# TODO Implement normal variance likelihood

lnorm.var <- function(y, mu=F, lo=NA, hi=NA, lpoints=1000, model="", scale=T) {
    
    
    
}

# TODO Implement probability of finding weak or misleading evidence for normal models

# TODO Add likelihood for normal CV (NORM_CV.R)

# TODO Implement likelihood for difference between normal means

# TODO Implement likelihood for ratio of normal means