#' Negative binomial likelihood
#' 
#' Original code by Richard Royall Oct 18 2000
#' For vector of counts (y), generates likelihood for E[Y_i] under
#' iid Poisson model.
#' 
#' @param y Vector of counts.
#' @param r Scale parameter (optional)
#' @param lo Lower parameter bound to likelihood calculation (optional).
#' @param hi Upper parameter bound to likelihood calculation (optional).
#' @param robust Flag for calculating robust correction (defaults to FALSE).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#' @return likelihood object.
#' @keywords likelihood
#' @export
lnbinom <- function(y, r=NA, lo = NA, hi = NA, robust = F, scale= T)
{
    n <- length(y)
    if(is.na(hi)) hi <- mean(y) + 4 * sqrt(var(y)/(n - 1))
    if(is.na(lo)) lo <- max(0.00001, mean(y) - 4 * sqrt(var(y)/(n - 1)))
    z <- seq(lo, hi, (hi - lo)/1000)
}