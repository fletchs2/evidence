# TODO Add unit tests for all functions

library(testthat)

test_that("Poisson likelihood validation",
    {
        x <- rpois(5, 3)
        a <- lpois(x)
        pl <- plot(a)
        b <- lpois(x, robust=T)
        plot(b, pl)
        
        expect_that(a$x[a$lx==1], equals(mean(x)))
    }
)