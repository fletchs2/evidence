# =====================
# = Utility functions =
# =====================

rmaxer.2 <- function(y, expval)
{
#	FUNCTION: "rmaxer.2"       May 20/96
#  For sample y finds integer value of parameter "r"
#  in negative binomial model with fixed mean
#  (expval) that maximizes the likelihood, L(expval,r)
#  (For use in profile likelihood) 
#  Used in "negbinmean.lik" and "multipoisson.lik"
#
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
