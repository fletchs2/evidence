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

# TODO Add probability of failing to find strong evidence to Poisson.

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
#' @param p.fail Flag for calculating probability of failing to find strong  
#'      evidence supporting trueprob (defaults to FALSE).
epois <- function(n, truemean, truevar=F, r=NA, lo=NA, hi=NA, k=8, robust=F, approx=F, p.fail=F) {
    
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
    
    error <- list(mislead=list(x=mean, px=mislead), dist="Poisson")
    class(error) <- "error"
    error
}

# TODO Add likelihood for ratio of Poisson means