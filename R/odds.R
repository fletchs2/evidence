#' Odds ratio likelihood for binomial model
#'
#' For binomial counts x and y, with sample sizes m and n, respectively,
#' calculates conditional likelihood.
#' 
#' @param y1 Binomial count, sample 1.
#' @param n1 Binomial sample size, sample 1.
#' @param y2 Binomial count, sample 2.
#' @param n2 Binomial sample size, sample 2.
#' @param lo Lower parameter bound to likelihood calculation (defaults to 0).
#' @param hi Upper parameter bound to likelihood calculation (defaults to 1).
#' @param points Number of points in [0,1] to calculate likelihood 
#'      (defaults to 1000).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#' @return likelihood object.
#' @keywords likelihood odds ratio binomial
#' @export
odds_like <- function(y1, n1, y2, n2, lo=NA, hi=NA, points=1000, scale=T) {
    
    k <- y1 + y2
    s <- max(k - n2, 0):min(n1, k)
    lcoeff <- lgamma(y1 + 1) + lgamma(n1 - y1 + 1) + lgamma(k - y1 + 1) + lgamma(n2 - k + y1 + 1) - lgamma(s + 1) - lgamma(n1 - s + 1) - lgamma(k - s + 1) - lgamma(n2 - k + s + 1)
        
    psi <- (y1 * (n2 - y2))/(y2 * (n1 - y1))  
    std <- sqrt(1/y1 + 1/y2 + 1/(n1 - y1) + 1/(n2 - y2))
	ll <- exp(log(psi) - 3 * std)
	ul <- exp(log(psi) + 2.7 * std)
    
    if(is.na(lo)) lo <- ll
    if(is.na(hi)) hi <- ul
    
    z <- seq(lo, hi, length=points)
    
    Tz <- matrix(z, nrow = length(z), ncol = length(
		s))
	Tz <- (t(Tz))^(s - y1)
	like <- 1/(t(Tz)) %*% exp(lcoeff)
    
    psihat <- z[like==max(like)]
    aa <- matrix(1, nrow=1, ncol=length(s))
	LR <- max(like) * (aa %*% exp(lcoeff))
    
    if (scale==T) {
        like <- like/max(like)
    }

    # Instantiate likelihood object
    likelihood <- list(x=z, lx=like, LR=LR)
    class(likelihood) <- "likelihood" 
    likelihood
}