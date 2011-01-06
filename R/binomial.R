#' Binomial likelihood
#'
#' For binomial count (y), generates likelihood for p under
#' iid Poisson model.
#' 
#' @param n Binomial sample size.
#' @param y Binomial count.
#' @param lo Lower parameter bound to likelihood calculation (defaults to 0).
#' @param hi Upper parameter bound to likelihood calculation (defaults to 1).
#' @param points Number of points in [0,1] to calculate likelihood 
#'      (defaults to 1000).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#' @return likelihood object.
#' @keywords likelihood
#' @export
lbinom <- function(n, y, lo=0, hi=1, points=1000, scale=T) {
    
    # Range of parameter
    p <- seq(lo, hi, length=points)
    # Calculate likelihood over range of p
    like <- exp(y*log(p) + (n-y) * log(1 - p))
    if (scale==T) {
        like <- like/max(like)
    }

    # Instantiate likelihood object
    likelihood <- list(x=p, lx=like)
    class(likelihood) <- "likelihood"
    likelihood$dist <- "Binomial"
    likelihood
}

# TODO Add probability of failing to find strong evidence to binomial model

#' Probability of misleading evidence for binomial distribution
#' 
#' L(p)/L(trueprob) >= k
#' 
#' @param n Binomial sample size.
#' @param trueprob True probability of success.
#' @param lo Lower parameter bound to likelihood calculation (defaults to 0).
#' @param hi Upper parameter bound to likelihood calculation (defaults to 1).
#' @param p.fail Calculates probability of failing to find strong evidence 
#'      supporting trueprob.
#' @keywords likelihood
#' @export
ebinom <- function(n, trueprob, lo=0, hi=1, k=8, points=1000, p.fail=F) {
    
    z <- seq(lo+0.001, hi-0.001, by=1/points)
    x <- 0:n
    
    mislead <- numeric(length(z))
    if (p.fail==T) {
        fail <- numeric(length(z))
    } else fail <- NULL
    
    for (i in 1:length(z)) {
        
        x.mis <- x[dbinom(x, n, z[i]) >= k*dbinom(x, n, trueprob)]
        
        mislead[i] <- ifelse(length(x.mis)==0, 0, sum(dbinom(x.mis, n, trueprob)))
        
        # Probability of failing to find strong evidence supporting trueprob
        if (p.fail==T) {
            
            x.fail <- x[dbinom(x, n, z[i]) > dbinom(x, n, trueprob)/k]
            
            fail[i] <- ifelse(length(x.fail)==0, 0, sum(dbinom(x.fail, n, trueprob)))
        }
    }
    
    error <- list(mislead=list(x=z, px=mislead), fail=list(x=z, px=fail),
                true=trueprob, dist="Binomial")
    class(error) <- "error"
    error
}