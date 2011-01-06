# TODO Add gamma likelihood and probability of finding weak or misleading evidence
# TODO Add robust likelihood for gamma (GAMMAROB.R)

#' Gamma likelihood
#' 
#' Original code by Richard Royall Oct 26 2000
#' For vector of counts (y), generates likelihood for E[Y_i] under
#' iid gamma model.
#' 
#' @param y Vector of counts.
#' @param r Fixed scale parameter (optional)
#' @param lo Lower parameter bound to likelihood calculation (optional).
#' @param hi Upper parameter bound to likelihood calculation (optional).
#' @param robust Flag for calculating robust profile likelihood function. Can be used with fixed or unknown r (defaults to FALSE).
#' @param normal Normal model specification (defaults to FALSE).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#' @return likelihood object.
#' @keywords likelihood
#' @export