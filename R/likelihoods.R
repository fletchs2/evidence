# ==============================================
# = Statistical likelihood objects and methods =
# ==============================================
library(ggplot2)

#' Method for calculating likelihood intervals
#' 
interval <- function(x, ...) { UseMethod("interval") }

interval.likelihood <- function(x, interval) {
    maxlike <- max(x$lx)
    interval.points <- x$x[x$lx >= (maxlike/interval)]
    list(endpoints=c(min(interval.points), max(interval.points)), like=maxlike/interval)
}

interpolate.likelihood <- function(like_obj) {
	interpolate_points <- (approx(like_obj$x, like_obj$lx[,1], n=10000))
	z <- interpolate_points$x
	like <- interpolate_points$y
	LR <- like_obj$LR
    likelihood <- list(x=z, lx=like, LR=LR)
    class(likelihood) <- "likelihood" 
    likelihood
}

#' Method for plotting likelihood
#' 
#' @param like_obj likelihood A likelihood object
#' @param y NA Not used
plot.likelihood <- function(like_obj, y=NA, intervals=c(8, 32), 
    xlabel="", ylabel="Likelihood", main="", color="black", linetype=1, 
    int.color="black", int.linetype=1) {
    # Create data frame
    dframe <- data.frame(x=like_obj$x, lx=like_obj$lx)
    

    # Plot likelihood curve
    pl <- ggplot(dframe, aes(x)) + geom_line(aes(y=lx), color=color, size=0.6)
    pl <- pl + opts(axis.line=theme_blank(), 
        panel.background=theme_rect(),
        panel.grid.major=theme_blank(),
        panel.grid.minor=theme_blank(),
        axis.text.x=theme_text(colour = "black", vjust = 1),
        axis.text.y=theme_text(colour = "black", hjust = 1),
        title=main) + ylab(ylabel) + xlab(xlabel)
    
    
    for (i in 1:length(intervals)) {
        # Calculate intervals
        if(length(like_obj$x) < 1000) int <- interval(interpolate.likelihood(like_obj), intervals[i])
        else int <- interval(like_obj, intervals[i]) 
        # Plot intervals
        pl <- pl + geom_line(data=data.frame(int), 
            aes(endpoints, like), color=int.color, linetype=int.linetype, size=0.3) 
        # Plot labels
        labs <- data.frame(x=like_obj$x[like_obj$lx==max(like_obj$lx)], y=int$like, 
            int=intervals[i])
        pl <- pl + geom_text(aes(x = x, y = y, 
            label = paste("1/",int,sep="")), 
            data=labs, 
            size = 3, 
            hjust = 0, 
            vjust = -1)
    }
    pl
}

#' Method for plotting probabilites of misleading evidence
#' 
#' @param pmis_obj pmis An pmis object
#' @param y NA Not used
plot.pmis <- function(pmis_obj, y=NA, xlabel="", ylabel="Probability", main="") {
    # Create data frame
    dframe <- data.frame(x=pmis_obj$x, px=pmis_obj$px)
    #Plot likelihood curve
    pl <- ggplot(dframe, aes(x)) + geom_line(aes(y=px), size=0.6)
    pl <- pl + opts(axis.line=theme_blank(), 
        panel.background=theme_rect(),
        panel.grid.major=theme_blank(),
        panel.grid.minor=theme_blank(),
        axis.text.x=theme_text(colour = "black", vjust = 1),
        axis.text.y=theme_text(colour = "black", hjust = 1),
        title=main) + ylab(ylabel) + xlab(xlabel)
    if (!is.null(pmis_obj$pmis$px)) {
        dframe <- data.frame(pmis_obj$pmis)
        pl <- pl + geom_line(data=dframe, aes(y=px), size=0.6, linetype="dashed")
    }
    pl
}

add_plot <- function(x, plot, ...) { UseMethod("add_plot") }

# Method for adding likelihood to existing plot
add_plot.likelihood <- function(likelihood, existing.plot, size=0.6, linetype=2) {
    dframe <- data.frame(x=likelihood$x, lx=likelihood$lx)
    existing.plot <- existing.plot + geom_line(data=dframe, aes(y=lx), 
        size=size, linetype=linetype)
    existing.plot
}

add_plot.pmis <- function(pmis, existing.plot, size=0.3) {
    dframe <- data.frame(pmis$mislead)
    existing.plot <- existing.plot + geom_line(data=dframe, aes(y=px), size=size)
    
    if (!is.null(pmis$pmis$px)) {
        dframe <- data.frame(pmis$pmis)
        existing.plot <- existing.plot 
            + geom_line(data=dframe, aes(y=px), size=size, linetype="dashed")
    }
    existing.plot
}




