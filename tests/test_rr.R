library(testthat)

test_that("Relative risk likelihood validation",
    {
        a <- rr_like(10,20,12,20,points=10000)
        pl <- plot(a)
        
        expect_that(round(a$x[a$lx==1],2), equals(round((10/20)/(12/20),2)))
    }
)

test_that("Relative risk likelihood validation",
    {
        a <- rr_like(99,1000,10,1000,lo=0,hi=50,points=10000)
        pl <- plot(a)
        expect_that(round(max(a$x[a$lx==1],na.rm=TRUE),1), equals(round((99/1000)/(10/1000),1)))
    }
)


test_that("Relative risk likelihood validation",
    {     
       expect_that(rr_like(10,10000,12,20), throws_error())
    }
)


test_that("Relative risk likelihood validation",
    {     
       expect_that(rr_like(16,16,20,20), throws_error())
    }
)

