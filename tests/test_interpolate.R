library(testthat)


test_that("Interpolate likelihood interval validation",
    {
        a <- rr_like(12,20,15,20,points=10)
        plot(a)
        
    }
)

test_that("Interpolate likelihood interval validation",
    {
        a <- rr_like(12,20,15,20,points=5)
        plot(a)
        
    }
)

test_that("Interpolate likelihood interval validation",
    {
        a <- rr_like(12,20,15,20,points=2)
        plot(a)
        
    }
)
