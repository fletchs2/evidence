#' Evidential sample size calculation 
#'
#' Evidential sample size required to ensure the specified strong evidence bound for each 
#' alternatively hypothesized delta.
#' 
#' @param k A vector of likelihood ratios of k-strength (defaults to c(8,20,32)).
#' @param S A vector of specified strong evidence bounds (defaults to c(0.80,0.85,0.90)).
#' @param delta A vector of effect sizes for each alternative hypotehsis.
#' @return results data frame.
#' @keywords likelihood sample size 
#' @export

ss.mean <- function(k=c(8,20,32),S=c(0.80,0.85,0.90), delta=c(deltas)){

k=k
if(k<1) stop("Frame your likelihood ratio such that k is greater than 1")

delta=delta
S=S

results <- data.frame(k = NA, S = NA, delta = NA, ss = NA, stringsAsFactors = FALSE)
i=1
for(x in 1:length(k)) {
	for(y in 1:length(delta)) {
		for(z in 1:length(S)) {
		
		ss <- (qnorm(S[z]) + sqrt(qnorm(S[z])^2 + 2 * log(k[x])))^2/(delta[y])^2
		
		results[i, "k"] <- k[x]
   		results[i, "S"] <- S[z]
  		results[i, "delta"] <- delta[y]
   		results[i, "ss"] <- ss
   		
   		i<-i+1

		}
	}
}
names(results) = c("k", "Strong.evidence", "Effect.size","Sample.size")

return(results)
 
}








