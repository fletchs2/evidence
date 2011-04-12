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
nbinomial_like <- function(y, r=NA, lo=NA, hi=NA, robust=F, model="", scale=T)
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
    likelihood
}

#' Probability of weak or misleading evidence for negative binomial distribution
#' 
#' L(p)/L(trueprob) >= k
#' 
#' @param r Negative binomial r parameter.
#' @param trueprob True probability of success.
#' @param lo Lower parameter bound to likelihood calculation (defaults to 0).
#' @param hi Upper parameter bound to likelihood calculation (defaults to 1).
#' @param p.weak Calculates probability of failing to find strong evidence 
#'      supporting trueprob.
#' @keywords likelihood
#' @export
enbinom <- function(r, trueprob, lo=0, hi=1, k=8, points=1000, weak=F) {
    
    if (!r > 0) stop("Parameter r must be positive in enbionom")
    
    z1 <- seq(lo+0.001, trueprob-0.001, by=1/points)
    z2 <- seq(trueprob+0.001, hi-0.001, by=1/points)
    z <- c(z1, trueprob, z2)
    
    if (weak==T) {
        fail.hi <- pnbinom(floor(( - log(k) - r * log(z2/trueprob))/log((1 - z2)/(1 - trueprob))) - 0.001, r, trueprob) -pnbinom((log(k) - r * log(z2/trueprob))/log((1 - z2)/(1 - trueprob)), r, trueprob)
        fail.lo <- pnbinom(floor((log(k) - r * log(z1/trueprob))/log((1 - z1)/(1 - trueprob))) - 0.001, r, trueprob) - pnbinom(( - log(k) - r * log(z1/trueprob))/log((1 - z1)/(1 - trueprob)), r, trueprob)
        fail <- c(fail.lo, 1, fail.hi)
        fail <- list(x=z, px=fail)

    } else {
        mislead.hi <- pnbinom((log(k) - r * log(z2/trueprob))/log((1 - z2)/(1 - trueprob)), r, trueprob)
        mislead.lo <- 1 - pnbinom(floor((log(k) - r * log(z1/trueprob))/log((1 - z1)/(1 - trueprob)) - 0.001), r, trueprob)
        mislead <- c(mislead.lo, 0, mislead.hi)
        fail <- list(x=z, px=mislead)
    }
    
    class(fail) <- "fail"
    fail
}


#' Finds integer value for negative binomial parameter r that maximizes
#' the likelihood for fixed mean expval. Used in profile likelihood.
#'
#' Original code by Richard Royall May 1996.
# 
# @param y Observed data.
# @param expval Expected value for negative binomial.
rmax <- function(y, expval)
{
	z <- 10000000
	newz <- 1000000
	llikz <- sum(log(dnbinom(y, z, z/(z + expval)))
		)
	while(z > 1) {
		lliknewz <- sum(log(dnbinom(y, newz, newz/(newz + expval))))
		if(llikz <= lliknewz) {
			z <- newz
			llikz <- lliknewz
			newz <- newz/10
		}
		else break
	}
	if(z == 10000000)
		z
	else if(z <= 10) {
		z <- 1
		repeat {
			if(sum(log(dnbinom(y, z, z/(z + 
				expval)))) < sum(log(
				dnbinom(y, z + 1, (z + 
				1)/((z + 1) + expval)))
				))
				z <- z + 1
			else break
		}
	}
	else {
###(now z/10 < max < z*10 )###
		c1 <- 10^0.25
		z <- z/10
		llikz <- sum(log(dnbinom(y, z, z/(z + 
			expval))))
		repeat {
			newz <- ceiling(c1 * z)
			lliknewz <- sum(log(dnbinom(y, newz, newz/(newz + 
				expval))))
			if(llikz <= lliknewz) {
				z <- newz
				llikz <- lliknewz
			}
			else {
				z <- newz
				llikz <- lliknewz
				break	
	###(NOW MAX TO LEFT OF Z, BETWEEN Z/(c1)^2 AND Z)###
			}
		}
		f1 <- 1/(c1)^0.1
		repeat {
			newz <- ceiling(f1 * z)
			lliknewz <- sum(log(dnbinom(y, newz, newz/(newz + 
				expval))))
			if(llikz <= lliknewz) {
				z <- newz
				llikz <- lliknewz
			}
			else break
		}
	}
	z
}
