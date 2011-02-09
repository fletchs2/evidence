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
    likelihood$name <- "Mean"
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
	var <- seq(lo, hi, length=lpoints)
	
	# Estimated likelihood
	if (estimated==TRUE) {
	    like <- log((ss/(n*var))^(n/2)) - (ss/(n*var) - 1) * n/2
	}
	# Profile likelihood
	else if (profile==TRUE) {
	    like <- log((ss/((n - 1)*var))^(n/2)) - ss/(2*var)
	}
	# Likelihood function
	else {
	    # Assuming mu is known
    	if (!is.na(mu)) {
    	    ss <- sum((y - mu)^2)
    	    like <- (1 + log(ss/(n*var)) - ss/(n*var))*n/2
    	}
    	like <- (1 + log(ss/((n - 1)*var)) - ss/((n - 1)*var))*(n-1)/2
	}
	
	if (scale==T) {
        like <- like - max(like)
    }
	
	# Return likelihood object
    likelihood <- list(x=var, lx=exp(like))
    class(likelihood) <- "likelihood" 
    likelihood$name <- "Variance"
    likelihood
}


#' Probability of weak or misleading evidence for 
#' mean of normal distribution
#' 
#' Misleading: Pr(L(mean+delta)/L(mean) >k ; mean)
#' Weak: Pr(1/k< L(mean+delta)/L(mean < k ; mean)
#' 
#' @param mu Normal mean.
#' @param sigma Normal standard deviation.
#' @param hi Upper parameter bound to likelihood calculation.
#' @param weak Calculates probability of failing to find strong evidence 
#'      supporting mu.
#' @keywords error normal
#' @export

enorm.mean <- function(mu, sigma, k=8, hi=NA, lpoints=1000, weak=FALSE) {
    
    if(is.na(hi)) {
	    hi <- 10 * log(k)
	}
	n <- seq(0.001, hi, length = lpoints)

    if (weak==T) {
        bad <- pnorm(sqrt(n)/2 + log(k)/sqrt(n)) - 
                pnorm(sqrt(n)/2 - log(k)/sqrt(n))
    }
    else {
	    bad <- pnorm(-sqrt(n)/2 - log(k)/sqrt(n))
    }
    
	# Return error object
    error <- list(x=n, px=bad)
    class(error) <- "error"
    error$name <- "Sample size"
    error
	
}


#' Normal likelihood for coefficient of variation
#'
#' For data (y), calculates profile likelihood for coefficient under
#' iid Normal model.
#' 
#' @param y Observations.
#' @param lo Lower parameter bound to likelihood calculation (Defaults to -10).
#' @param hi Upper parameter bound to likelihood calculation (Defaults to 10).
#' @param lpoints The number of evenly-spaced points in the interval over
#' which to calculate profile likelihood (defaults to 1000).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE.
#' @return likelihood object.
#' @keywords likelihood normal cv
#' @export

lnorm.cv <- function(y, lo=-10, hi=10, lpoints=1000, scale=TRUE) {
    
    z <- seq(lo, hi, length = lpoints)
	n <- length(y)
	s.hat <- (-mean(y)/z + sqrt((mean(y)/z)^2 + 4*mean(y^2)))/2
	like <- -n*(log(s.hat) + 0.5*(mean(y^2)/(s.hat)^2 
	    - (2*mean(y))/(z*s.hat) + (1/z)^2))
	    
	if (scale==T) {
	    like <- like - max(like)
	}

	# Return likelihood object
	likelihood <- list(x=z, lx=exp(like))
    class(likelihood) <- "likelihood" 
    likelihood$name <- "Coefficient of variation"
    likelihood
}

#' Normal profile likelihood for differences between means
#'
#' Calculates profile likelihood for the difference between two
#' iid Normal means (x,y). Variances can be assumed equal or unequal. Uses
#' exponent (n-1)/2 instead of n/2 in profile likelihood for
#' normal mean, as suggested by Kalbfleisch and Sprott(1971).
#' 
#' @param x First set of observations.
#' @param y Second set of observations.
#' @param lo Lower parameter bound to likelihood calculation.
#' @param hi Upper parameter bound to likelihood calculation.
#' @param lpoints The number of evenly-spaced points in the interval over
#' which to calculate profile likelihood (defaults to 1000).
#' @param equal.var Flag for assuming variances of variables are equal.
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE.
#' @return likelihood object.
#' @keywords likelihood normal cv
#' @export

lnorm.diff <- function(x, y, lo=NA, hi=NA, equal.var=TRUE, lpoints=1000, scale=TRUE) {
    
    m <- length(x)
	n <- length(y)
	
	diff <- mean(x) - mean(y)
	
	# Determine range of x-axis
	k <- 3
	if(min(m, n) == 2)
		k <- 13
	if(min(m, n) == 3)
		k <- 5
	if(min(m, n) == 4)
		k <- 4.5
		
	se <- sqrt(var(x)/m + var(y)/n)
	
	if(is.na(lo))
		lo <- diff - k * se
	if(is.na(hi))
		hi <- diff + k * se
    
    z <- seq(lo, hi, length=lpoints)
    
    z.0 <- c(0, z)
    
    like <- numeric(length(z.0))
        
    like.mat <- matrix(1, length(z.0), 100)
    
    for (i in 1:length(z.0)) {
        
        if (equal.var==T) {
            tmp <- c(x - z.0[i], y)
            like[i] <- 1/((var(tmp)^((m+n-2)/2)))            
        } 
        else {
            mu.y <- seq(min(mean(x) - z.0[i], mean(y)), 
                max(mean(x) - z.0[i], mean(y)), length = ncol(like.mat))
            
            for (j in 1:ncol(like.mat)) {
                like.mat[i,j] <- 1/((mean((x - z.0[i] - mu.y[j])^2)^((m - 1)/2)) 
                    * (mean((y - mu.y[j])^2)^((n - 1)/2)))
            }

            like[i] <- max(like.mat[i,])            
        }
    }
    
    LR <- max(like)/like[1]
    like <- like[2:length(like)]
        
    if (scale==T) {
        like <- like/max(like)
    }
    
    # Return likelihood object
    likelihood <- list(x=z, lx=like)
    class(likelihood) <- "likelihood" 
    likelihood$name <- "Difference of means"
    likelihood
}

#' Normal profile likelihood for ratio of two means
#'
#' Calculates profile likelihood for the ratio of two
#' iid Normal means (x,y). Variances can be assumed equal or unequal. Uses
#' exponent (n-1)/2 instead of n/2 in profile likelihood for
#' normal mean, as suggested by Kalbfleisch and Sprott(1971).
#' 
#' @param x First set of observations.
#' @param y Second set of observations.
#' @param lo Lower parameter bound to likelihood calculation.
#' @param hi Upper parameter bound to likelihood calculation.
#' @param lpoints The number of evenly-spaced points in the interval over
#' which to calculate profile likelihood (defaults to 1000).
#' @param equal.var Flag for assuming variances of variables are equal.
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE.
#' @return likelihood object.
#' @keywords likelihood normal cv
#' @export

lnorm.ratio <- function(x, y, lo=0, hi=10, equal.var=TRUE, lpoints=1000, scale=TRUE) {
    
    m <- length(x)
	n <- length(y)
	
	z <- seq(lo, hi, length=lpoints)
		
	if (equal.var==T) {
	    like <-  - ((m + n - 1)/2) * log(((m - 1) * var(x) + (n - 1) * var(y))/(m + n) + 
	        (((m * n)/(m + n)) * ((z * mean(x) - mean(y))^2))/(m + n * z^2))
	}
	else {
	    z.1 <- c(1, z)

    	logfactor.1 <-  - ((m - 1)/2) * log((var(x) * (
    		m - 1))/m + ((n * z.1 * (z.1 * mean(x) - 
    		mean(y)))/(m + n * z.1^2))^2)
    	logfactor.2 <-  - ((n - 1)/2) * log((var(y) * (
    		n - 1))/n + ((m * (z.1 * mean(x) - mean(
    		y)))/(m + n * z.1^2))^2)
    		
    	like <- (logfactor.1 + logfactor.2 - max(logfactor.1 + logfactor.2))[2:(lpoints+1)]
	}
	
	if (scale==T) {
	    like <- like - max(like)
	}
    
    # Return likelihood object
    likelihood <- list(x=z, lx=exp(like))
    class(likelihood) <- "likelihood" 
    likelihood$name <- "Ratio of means"
    likelihood
}
