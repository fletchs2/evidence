# TODO Add unit tests for all functions

library(testthat)

test_that("Poisson likelihood validation",
    {
        x <- rpois(30, 3)
        a <- poisson_like(x)
        pl <- plot(a)
        b <- poisson_like(x, robust=T)
        plot(b, pl)
        
        expect_that(a$x[a$lx==1], equals(mean(x)))
    }
)

test_that("Binomial likelihood validation",
    {
        x <- rbinom(30, 10, 0.2)
        a <- binomial_like(10, x)
        pl <- plot(a)
        
        expect_that(round(a$x[a$lx==1], 3), equals(round(mean(x/10), 3)))
    }
)