#' Relative risk likelihood for conditional binomial model
#'
#' For binomial counts x and y, with sample sizes m and n, respectively,
#' calculates conditional likelihood.
#' 
#' @param y1 Binomial count, sample 1.
#' @param n1 Binomial sample size, sample 1.
#' @param y2 Binomial count, sample 2.
#' @param n2 Binomial sample size, sample 2.
#' @param lo Lower parameter bound to likelihood calculation.
#' @param hi Upper parameter bound to likelihood calculation.
#' @param points Number of points in [0,1] to calculate likelihood 
#'      (defaults to 1000).
#' @param scale Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#' @return likelihood object.
#' @keywords likelihood relative risk binomial
#' @export
rr_like <- function(y1, n1, y2, n2, lo=0, hi=6, points=1000, scale=T) {

if(n1>1000 | n2>1000) stop("Sample sizes must be less than 1001")

if((n1 - y1) * (n2 - y2) == 0)
     stop("Conditional model requires at least one failure in each group.") else 
	 {
     z <- seq(lo, hi, length=points)
     
    like <- matrix(1,nrow=length(z), ncol=1) 
    j <- seq(n1-y1,n1+y2,1)	## summation index
    Mx <- matrix(z, nrow = length(j), ncol = length(z), byrow=T)
    Mx <- Mx^(j - n1)
     
    co <- exp(lgamma(n1 + n2 - j - 1 + 1)-lgamma(n2 - y2 - 1 + 1)-lgamma((n1 + n2 - j - 1)-(n2 - y2 - 1) + 1)-max(lgamma(n1 + n2 - j - 1 + 1)-lgamma(n2 - y2 - 1 + 1)-lgamma((n1 + n2 - j - 1)-(n2 - y2 - 1) + 1)))*exp(lgamma(j - 1 + 1)-lgamma(n1 - y1 - 1 + 1)-lgamma((j - 1)-(n1 - y1 - 1) + 1)-max(lgamma(j - 1 + 1)-lgamma(n1 - y1 - 1 + 1)-lgamma((j - 1)-(n1 - y1 - 1) + 1)))
    
    like <- 1/(t(Mx) %*% co)

    if(y2== 0) like <- like * exp(co[length(j)])   
    if(scale==T)  like <- like/max(like[like!=Inf],na.rm=TRUE)
    
	rrhat <- z[like == max(like[like!=Inf],na.rm=TRUE)]
	rrhat <- max(rrhat,na.rm=TRUE)
	rrhat <- round(rrhat, digits = 3)

	if(like[length(z)] == max(like[like!=Inf],na.rm=TRUE)) rrhat <- NA
	if(is.na(like[1])==FALSE & like[1] == max(like[like!=Inf],na.rm=TRUE)) rrhat <- NA
	
	LR <- round(min(max(like[like!=Inf],na.rm=TRUE)/like[abs(z - 1) == min(abs(z - 1))],1000000), digits = 1)
	
# Instantiate likelihood object
    likelihood <- list(x=z, lx=like, LR=LR)
    class(likelihood) <- "likelihood" 
    likelihood
    
    }
    
}
