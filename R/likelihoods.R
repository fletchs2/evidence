# ==============================================
# = Statistical likelihood objects and methods =
# ==============================================
library(ggplot2)

# Define likelihood class
setClass("likelihood", representation(x="numeric", lx="numeric"))

# Method for plotting likelihood
setMethod("plot",
    signature=c("likelihood"),
    function(x, y, intervals=c(8, 32), xlabel="", ylabel="Likelihood", main="") {
        # Create data frame
        dframe <- data.frame(x=x@x, lx=x@lx)
        # Plot likelihood curve
        pl <- ggplot(dframe, aes(x)) + geom_line(aes(y=lx), size=0.6)
        pl <- pl + opts(axis.line=theme_blank(), 
        panel.background=theme_rect(), 
        panel.grid.major=theme_blank(), 
        panel.grid.minor=theme_blank(),
        axis.text.x=theme_text(colour = "black", vjust = 1),
        axis.text.y=theme_text(colour = "black", hjust = 1),
        title=main) + ylab(ylabel) + xlab(xlabel) 
        maxlike <- max(x@lx)
        for (i in 1:length(intervals)) {
            # Calculate intervals
            xx <- dframe$x[dframe$lx >= (maxlike/intervals[i])]
            yy <- rep(maxlike/intervals[i], length(xx))
            # Plot intervals
            interval <- data.frame(x=xx, y=yy)
            pl <- pl + geom_line(data=interval, aes(x, y), size=0.3) 
            # Plot labels
            labels <- data.frame(x=x@x[x@lx==maxlike], y=yy[1], int=intervals[i])
            pl <- pl + geom_text(aes(x = x, y = y, label = paste("1/",int,sep="")), data=labels, size = 3, hjust = 0, vjust = -1)
        }
        pl
    }
)

# Poisson likelihood
# 
# Original code by Richard Royall Oct 17 2000
# For vector of counts (y), generates likelihood for E[Y_i] under
# iid Poisson model.
# 
# Args:
#     y: Vector of counts.
#     lo: Lower parameter bound to likelihood calculation (optional).
#     hi: Upper parameter bound to likelihood calculation (optional).
#     robust: Flag for calculating robust correction (defaults to FALSE).
#     scale: Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#     
# Returns:
#     likelihood object.

lpois <- function(y, lo=NA, hi=NA, robust=F, scale=T)
{

    # Set default bounds if none provided
	if(is.na(lo)) lo <- mean(y) * max(1 - 4/sqrt(sum(y)), 1/4)
	if(is.na(hi)) hi <- mean(y) * (1 + 4/sqrt(sum(y)))
	# Range of parameter
	z <- seq(lo, hi, (hi - lo)/1000)
	# Calculate likelihood over range of z
	like <- sum(y) * log(z) - length(y) * z
	# Correct for model failure
	if(robust == T) {
      like <- (like * (mean(y) * length(y)))/
          (var(y) * (length(y) - 1))
    }
	if (scale==T) {
	    like <- like - max(like)
	}
	
    # Instantiate likelihood object
    new("likelihood", x=z, lx=exp(like))
}

# Binomial likelihood
#
# For binomial count (y), generates likelihood for p under
# iid Poisson model.
# 
# Args:
#     n: Binomial sample size.
#     y: Binomial count.
#     lo: Lower parameter bound to likelihood calculation (defaults to 0).
#     hi: Upper parameter bound to likelihood calculation (defaults to 1).
#     points: Number of points in [0,1] to calculate likelihood 
#         (defaults to 1000).
#     scale: Flag for scaling maximum likelihood to 1 (defaults to TRUE).
#     
# Returns:
#     likelihood object.

lbinom <- function(n, y, lo=0, hi=1, points=1000, scale=T) {
    
    # Range of parameter
    p <- seq(lo, hi, length=points)
    # Calculate likelihood over range of p
    like <- exp(y*log(p) + (n-y) * log(1 - p))
    if (scale==T) {
        like <- like/max(like)
    }

    # Instantiate likelihood object
    new("likelihood", x=p, lx=like)
}

x <- c(2,3,3,4,5,5)
a <- lpois(x, scale=F)
plot(a)