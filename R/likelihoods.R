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

#' Method for plotting likelihood
#' 
#' @param like.obj likelihood A likelihood object
#' @param y NA Not used
plot.likelihood <- function(like.obj, y=NA, intervals=c(8, 32), xlabel="", ylabel="Likelihood", main="") {
    # Create data frame
    dframe <- data.frame(like.obj$x, like.obj$lx)
    # Plot likelihood curve
    pl <- ggplot(dframe, aes(x)) + geom_line(aes(y=lx), size=0.6)
    pl <- pl + opts(axis.line=theme_blank(), 
        panel.background=theme_rect(),
        panel.grid.major=theme_blank(),
        panel.grid.minor=theme_blank(),
        axis.text.x=theme_text(colour = "black", vjust = 1),
        axis.text.y=theme_text(colour = "black", hjust = 1),
        title=main) + ylab(ylabel) + xlab(xlabel)
    for (i in 1:length(intervals)) {
        # Calculate intervals
        int <- interval(like.obj, intervals[i])
        # Plot intervals
        pl <- pl + geom_line(data=data.frame(int), aes(endpoints, like), size=0.3) 
        # Plot labels
        labels <- data.frame(x=like.obj$x[like.obj$lx==max(like.obj$lx)], y=int$like, int=intervals[i])
        pl <- pl + geom_text(aes(x = x, y = y, 
            label = paste("1/",int,sep="")), 
            data=labels, 
            size = 3, 
            hjust = 0, 
            vjust = -1)
    }
    pl
}

#' Method for plotting probabilites of misleading evidence
#' 
#' @param error.obj error An error object
#' @param y NA Not used
plot.error <- function(error.obj, y=NA, xlabel="", ylabel="Probability", main="") {
    # Create data frame
    dframe <- data.frame(error.obj$mislead)
    # Plot likelihood curve
    pl <- ggplot(dframe, aes(x)) + geom_line(aes(y=px), size=0.6)
    pl <- pl + opts(axis.line=theme_blank(), 
        panel.background=theme_rect(),
        panel.grid.major=theme_blank(),
        panel.grid.minor=theme_blank(),
        axis.text.x=theme_text(colour = "black", vjust = 1),
        axis.text.y=theme_text(colour = "black", hjust = 1),
        title=main) + ylab(ylabel) + xlab(xlabel)
    if (!is.null(error.obj$fail$px)) {
        dframe <- data.frame(error.obj$fail)
        pl <- pl + geom_line(data=dframe, aes(y=px), size=0.3, linetype="dashed")
    }
    pl
}

add.plot <- function(x, plot, ...) { UseMethod("add.plot") }

# Method for adding likelihood to existing plot
add.plot.likelihood <- function(likelihood, existing.plot, size=0.6, linetype=2) {
    dframe <- data.frame(x=likelihood$x, lx=likelihood$lx)
    existing.plot <- existing.plot + geom_line(data=dframe, aes(y=lx), size=size, linetype=linetype)
    existing.plot
}

add.plot.error <- function(error, existing.plot, size=0.6) {
    dframe <- data.frame(error$mislead)
    existing.plot <- existing.plot + geom_line(data=dframe, aes(y=px), size=size)
    
    if (!is.null(error$fail$px)) {
        dframe <- data.frame(error$fail)
        existing.plot <- existing.plot + geom_line(data=dframe, aes(y=px), size=size, linetype="dashed")
    }
}




