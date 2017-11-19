library(gwpcR)
context('Parameter Estimation')

EPS.ABS <- 5e-2
EPS.REL <- 5e-2

test_mom <- function(efficiency, lambda0, n.obs, molecules, threshold) {
  test_that(paste0('method of moments (E#', efficiency, ', D=lambda0',
                   ', n.obs=', n.obs, ', mol=', molecules, ', th=', threshold,
                   ')'), {
    x <- rgwpcrpois(n = n.obs, efficiency=efficiency, lambda0=lambda0,
                    molecules=molecules, threshold=threshold)
    m <- gwpcrpois.est(x, method='mom', molecules=molecules, threshold=threshold)
    
    expect_lt(abs(efficiency - m$efficiency), EPS.ABS)
    expect_lt(abs(lambda0 - m$lambda0)/lambda0, EPS.REL)
    expect_equal(m$n.obs, n.obs)
    expect_equal(m$n.umis, n.obs)
    expect_equal(m$molecules, molecules)
    expect_equal(m$threshold, threshold)
  })
}

test_mom(0.01, 18, 1e4, 1, 1)
test_mom(0.05, 17, 1e4, 1, 1)
test_mom(0.1 , 16, 1e4, 1, 1)
test_mom(0.2 , 15, 1e4, 1, 1)
test_mom(0.3 , 14, 1e4, 1, 1)
test_mom(0.4 , 13, 1e4, 1, 1)
test_mom(0.5 , 12, 1e4, 1, 1)
test_mom(0.6 , 11, 1e4, 1, 1)
test_mom(0.7 , 10, 1e4, 1, 1)
test_mom(0.8 , 9, 1e4, 1, 1)
test_mom(0.9 , 8, 1e4, 1, 1)
test_mom(0.95, 7, 1e4, 1, 1)
test_mom(0.99, 6, 1e4, 1, 1)
test_mom(1.0 , 5, 1e4, 1, 1)
