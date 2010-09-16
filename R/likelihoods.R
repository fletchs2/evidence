library(ggplot2)

setClass("likelihood", representation(x="numeric", lx="numeric"))

setMethod("plot",
    c("likelihood"),
    function(x, y, intervals = c(8, 32), xlabel="", ylabel="Likelihood") {
        dframe <- data.frame(x=x@x, lx=x@lx)
        pl <- ggplot(dframe, aes(x)) + geom_line(aes(y=lx), size=0.6)
        pl <- pl + opts(axis.line=theme_segment(), 
        panel.background=theme_blank(), 
        panel.grid.major=theme_blank(), 
        panel.grid.minor=theme_blank(),
        axis.text.x=theme_text(colour = "black", vjust = 1),
        axis.text.y=theme_text(colour = "black", hjust = 1)) + ylab(ylabel) + xlab(xlabel)
        for (i in 1:length(intervals)) {
            xx <- dframe$x[dframe$lx >= 1/intervals[i]]
            yy <- rep(1/intervals[i], length(xx))
            interval <- data.frame(x=xx, y=yy)
            pl <- pl + geom_line(data=interval, aes(x, y), size=0.3) 
            labels <- data.frame(x=x@x[x@lx==1], y=yy[1]+0.005, int=intervals[i])
            pl <- pl + geom_text(aes(x = x, y = y, label = paste("1/",int)), data=labels, size = 3, hjust = 0, vjust = 0)
        }
        pl
    }
)

# Poisson likelihood
# 
# Original code by Richard Royall Oct 17 2000
# For vector of counts (y), generates likelihood for E[Y_i] under
# iid Poisson model.

lpois <- function(y, lo = NA, hi = NA, robust = F)
{

    # Set default bounds if none provided
	if(is.na(lo)) lo <- mean(y) * max(1 - 4/sqrt(sum(y)), 1/4)
	if(is.na(hi)) hi <- mean(y) * (1 + 4/sqrt(sum(y)))
	z <- seq(lo, hi, (hi - lo)/1000)
	llik <- sum(y) * log(z) - length(y) * z
	lik <- exp(llik - max(llik))
	
	# Correct for model failure
	if(robust == T) {
      llik <- (llik * (mean(y) * length(y)))/
          (var(y) * (length(y) - 1))
      lik <- exp(llik - max(llik))
    }
    new("likelihood", x=z, lx=lik)
}

a <- lpois(rpois(100, 3))
plot(a)
showClass("likelihood")