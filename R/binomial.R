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
binomial_like <- function(n, y, lo=0, hi=1, points=1000, scale=T) {
    
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
    likelihood
}


#' Probability of misleading evidence for binomial distribution
#' 
#' L(p)/L(trueprob) >= k
#' 
#' @param n Binomial sample size.
#' @param trueprob True probability of success.
#' @param lo Lower parameter bound to likelihood calculation (defaults to 0).
#' @param hi Upper parameter bound to likelihood calculation (defaults to 1).
#' @param weak Calculates probability of failing to find strong evidence 
#'      supporting trueprob.
#' @keywords likelihood
#' @export
binomial_error <- function(n, trueprob, lo=0, hi=1, k=8, points=1000, weak=F) {
    
    z <- seq(lo+0.001, hi-0.001, by=1/points)
    x <- 0:n
    
    bad <- numeric(length(z))
        
    for (i in 1:length(z)) {
        
        
        # Probability of failing to find strong evidence supporting trueprob
        if (weak==T) {
            x.e <- x[dbinom(x, n, z[i]) > dbinom(x, n, trueprob)/k]
        }
        else {
            x.e <- x[dbinom(x, n, z[i]) >= k*dbinom(x, n, trueprob)]
        }
        bad[i] <- ifelse(length(x.e)==0, 0, sum(dbinom(x.e, n, trueprob)))
    }
    
    error <- list(x=z, px=bad)
    class(error) <- "error"
    error
}