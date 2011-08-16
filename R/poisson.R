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
poisson_like <- function(y, lo=NA, hi=NA, robust=F, scale=T)
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
    likelihood
}


#' Probability of misleading evidence for Poisson distribution
#' 
#' L(p)/L(trueprob) >= k
#' 
#' @param n Binomial sample size.
#' @param truemean True mean.
#' @param truevar True variance for normal approximation only (optional).
#' @param lo Lower parameter bound to likelihood calculation (optional).
#' @param hi Upper parameter bound to likelihood calculation (optional).
#' @param robust Flag for calculating robust correction (defaults to FALSE)
#' @param p.pmis Flag for calculating probability of failing to find strong  
#'      evidence supporting trueprob (defaults to FALSE).
poisson_pmis <- function(n, truemean, truevar=F, r=NA, lo=NA, hi=NA, 
        k=8, robust=F, approx=F, weak=F) {
    
    # Set default bounds if none provided
    if(is.na(lo)) lo <- max(truemean - 10 * sqrt(truemean/n), 0.001)
    if(is.na(hi)) hi <- truemean + 10 * sqrt(truemean/n)
    
    # Less than true mean
    s.lo <- (truemean - lo)/500
    mean.lo <- seq(lo, truemean - s.lo, s.lo)
    crit.lo <- pmax((log(k) + n*(mean.lo - truemean))/log(mean.lo/truemean),  0)

    # Greather than true mean
    s.hi <- (truemean - hi)/500
    mean.hi <- seq(hi, truemean-s.hi, s.hi)
    crit.hi <- (log(k) + n*(mean.hi - truemean))/log(mean.hi/truemean)
    
    mean <- c(mean.lo, truemean, mean.hi)
    
    # Normal approximation
    if (approx==T) {
        mislead <- 1 - pnorm(((sqrt(n) * abs(mean - truemean))/2 + (truemean * log(8))/(sqrt(n) * abs(mean - truemean)))/sqrt(truevar))
    }
    # Negative binomial model
    else if (!is.na(r)) {
        
        # a/b= negbin mean / negbin var = p
        if (robust==T) {
            k1 <- k^(1/(r/(r + truemean)))
			mislead.lo <- pnbinom(floor((n * (mean.lo - truemean) + log(k1))/log(mean.lo/truemean)), n * r, r/(r + truemean))
			mislead.hi <- 1 - pnbinom(floor((n * (mean.hi - truemean) + log(k1))/log(mean.hi/truemean)), n * r, r/(r + truemean))
    		mislead <- c(mislead.lo, 0, mislead.hi)
        }
        else {
            mislead.lo <- pnbinom(floor((n * (mean.lo - truemean) + log(k))/log(mean.lo/truemean)), n * r, r/(r + truemean))
            mislead.hi <- 1 - pnbinom(floor((n * (mean.hi - truemean) + log(k))/log(mean.hi/truemean)), n * r, r/(r + truemean))
    		mislead <- c(mislead.lo, 0, mislead.hi)
        }
    }
    else {
        mislead.lo <- ppois(floor(crit.lo), n*truemean)
        mislead.hi <- 1 - ppois(ceiling(crit.hi) - 1, n*truemean)
        mislead <- c(mislead.lo, 0, mislead.hi)
        
    }
    
    pmis <- list(x=mean, px=mislead)
    class(pmis) <- "pmis"
    pmis
}


#' Poisson profile likelihood for ratio of two means
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
#' @param robust Flag for returning robust adjusted likelihood.
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE.
#' @return likelihood object.
#' @keywords likelihood poisson mean
#' @export
poisson_ratio_like <- function(x, y, lo=NA, hi=NA, lpoints=1000, robust=FALSE, scale=TRUE) {
    
    # Set default bounds if none provided
    if(is.na(lo)) lo <- 0
    if(is.na(hi)) hi <- (2 * mean(y))/mean(x)
    
    # Range of parameter
    z <- seq(lo, hi, length=1000)
    
    m <- length(x)
	n <- length(y)
	b <- m/(m + n * z)
	
	like <- sum(x) * log(b) + sum(y) * log(1 - b) 
	    - max(sum(x) * log(b) + sum(y) * log(1 - b))
    
    if (robust==T) {
        ratio <- mean(y)/mean(x)
		adj <- (m * mean(y) + n * ratio^2 * mean(x)) / 
		    ((m * var(y) * (n - 1))/n + (n * ratio^2 * var(x) * (m - 1))/m)
    }
    else {
        adj <- 1
    }
    
    if (scale==T) {
        like <- like - max(like)
    }
    
    # Return likelihood object
    likelihood <- list(x=z, lx=exp(like)^adj)
    class(likelihood) <- "likelihood" 
    likelihood$name <- "Ratio of means"
    likelihood
    
}