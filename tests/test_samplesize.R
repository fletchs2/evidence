library(testthat)

test_that("Sample size likelihood validation",
    {
    	a <- ss.mean(k=8,S=c(0.80,0.85,.90,0.95), delta=c(.2,.5,.7,1.2))
		ss <- (qnorm(.80) + sqrt(qnorm(.8)^2 + 2 * log(8)))^2/(.2)^2
        
        expect_that(a$Sample.size[a$Strong.evidence==0.80&a$Effect.size==.2], equals(ss))
    }
)

test_that("Sample size likelihood validation",
    { 		
 		expect_that(ss.mean(k=1/8,S=c(0.80,0.85,.90,0.95), delta=c(.2,.5,.7)), throws_error())
    }
)


test_that("Sample size likelihood validation",
    { 		
 		expect_that(ss.mean(k=8,S=c(0.80,0.85,.90,0.95), delta=c(.2,.5,.7,"1.2")), throws_error())
    }
)
