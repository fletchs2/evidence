#' Odds ratio likelihood for binomial model
#'
#' For binomial counts x and y, with sample sizes m and n, respectively,
#' calculates conditional likelihood.
#' 
#' @param x Binomial count, sample 1.
#' @param m Binomial sample size, sample 1.
#' @param y Binomial count, sample 2.
#' @param n Binomial sample size, sample 2.
#' @param lo Lower parameter bound to likelihood calculation (defaults to 0).
#' @param hi Upper parameter bound to likelihood calculation (defaults to 1).
#' @param points Number of points in [0,1] to calculate likelihood 
#'      (defaults to 1000).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#' @return likelihood object.
#' @keywords likelihood
#' @export
odds_like <- function(x, m, y, n, lo=NA, hi=NA, points=1000, scale=T) {
    
    k <- x + y
    s <- max(k - n, 0):min(m, k)
    lcoeff <- lgamma(x + 1) + lgamma(m - x + 1) + lgamma(k - x + 1) + lgamma(n - k + x + 1) - lgamma(s + 1) - lgamma(m - s + 1) - lgamma(k - s + 1) - lgamma(n - k + s + 1)
        
    psi <- (x * (n - y))/(y * (m - x))  
    std <- sqrt(1/x + 1/y + 1/(m - x) + 1/(n - y))
	ll <- exp(log(psi) - 3 * std)
	ul <- exp(log(psi) + 2.7 * std)
    
    if(is.na(lo)) lo <- ll
    if(is.na(hi)) hi <- ul
    
    z <- seq(lo, hi, length=points)
    
    Tz <- matrix(z, nrow = length(z), ncol = length(
		s))
	Tz <- (t(Tz))^(s - x)
	like <- 1/(t(Tz)) %*% exp(lcoeff)
    
    psihat <- z[like==max(like)]
    aa <- matrix(1, nrow=1, ncol=length(s))
	LR <- max(like) * (aa %*% exp(lcoeff))
    
    if (scale==T) {
        like <- like/max(like)
    }

    # Instantiate likelihood object
    likelihood <- list(x=z, lx=like)
    class(likelihood) <- "likelihood" 
    likelihood
}