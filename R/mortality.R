#' Standardized mortality ratio likelihood
#'
#' Calculates likelihood function for expected value
#' of the standardized mortality ratio (SMR) 
#' i.e., E(observed deaths/expected deaths),  under Poisson
#' model for observed deaths, taking expected deaths to be fixed.
#' 
#' @param obs_deaths Observed number of deaths.
#' @param exp_deaths Expected number of deaths.
#' @param lo Lower parameter bound to likelihood calculation (defaults to 0.1).
#' @param hi Upper parameter bound to likelihood calculation (defaults to 5).
#' @param points Number of points in [0,1] to calculate likelihood 
#'      (defaults to 1000).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#' @return likelihood object.
#' @keywords likelihood mortality ratio poisson
#' @export
odds_like <- function(obs_deaths, exp_deaths, lo=0.1, hi=5, points=1000, scale=T) {
    
    # Range of ratio
    z <- seq(lo, hi, length=points)
    
    # Calculate likelihood
	like <- obs_deaths * log(z*exp_deaths) - z*exp_deaths - lgamma(obs_deaths + 1)
    
    if (scale==T) {
        like <- like - max(like)
    }

    # Instantiate likelihood object
    likelihood <- list(x=z, lx=exp(like))
    class(likelihood) <- "likelihood" 
    likelihood
}


#' Standardized ratio of mortality rates likelihood
#'
#' Calculates likelihood function for ratio of
#' SMR from two causes. 
#' i.e., E(d1/e1)/E(d2/e2)= [E(d1)/E(d2)]/[e1/e2], under independent
#' Poisson model for d1 and d2.
#' 
#' @param obs1 Observed number of deaths, cause 1.
#' @param exp1 Expected number of deaths, cause 1.
#' @param obs2 Observed number of deaths, cause 2.
#' @param exp2 Expected number of deaths, cause 2.
#' @param lo Lower parameter bound to likelihood calculation (defaults to 0.1).
#' @param hi Upper parameter bound to likelihood calculation (defaults to 5).
#' @param points Number of points in [0,1] to calculate likelihood 
#'      (defaults to 1000).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#' @return likelihood object.
#' @keywords likelihood mortality ratio poisson
#' @export
odds_like <- function(obs1, exp1, obs2, exp2, lo=0.1, hi=5, points=1000, scale=T) {
    
    # Range of ratio
    z <- seq(lo, hi, length=points)
    
    # Calculate likelihood
    c1 <- exp2/exp1
	like <- obs1 * log(z/(z + c1)) + obs2 * log(c1/(z + c1))
    
    if (scale==T) {
        like <- like - max(like)
    }

    # Instantiate likelihood object
    likelihood <- list(x=z, lx=exp(like))
    class(likelihood) <- "likelihood" 
    likelihood
}