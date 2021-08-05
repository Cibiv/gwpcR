library(gwpcR)
context('Galton-Watson PCR-Poisson (gwpcrpois) parameter estimation')

N <- 1e5
EPS.ABS <- 5e-2
EPS.REL <- 5e-2

test_that('Precompute Distribution', {
  suppressMessages(gwpcr.molecules.precompute(3))
  expect_true(TRUE)
})

test_est <- function(efficiency, lambda0, n.obs, molecules, threshold, method) {
  skip_on_cran()
  
  p0 <- pgwpcrpois(threshold - 1, efficiency=efficiency, lambda0=lambda0,
                   molecules=molecules, threshold=0)
  test_that(paste0('method of moments (E#', efficiency, ', D=', lambda0,
                   ', n.obs=', n.obs, ', mol=', molecules, ', th=', threshold,
                   ')'), {
    x <- rgwpcrpois(n = n.obs, efficiency=efficiency, lambda0=lambda0,
                    molecules=molecules, threshold=threshold)
    m <- gwpcrpois.est(x, method=method, molecules=molecules, threshold=threshold)
    
    expect_equal(m$convergence, 0)
    expect_lt(abs(efficiency - m$efficiency), EPS.ABS)
    expect_lt(abs(lambda0 - m$lambda0)/lambda0, EPS.REL)
    expect_lt(abs(p0 - m$p0), EPS.ABS)
    expect_equal(m$loss, m$p0)
    expect_equal(m$n.tot, m$n.umis / (1 - m$loss))
    expect_equal(m$n.obs, n.obs)
    expect_equal(m$n.umis, n.obs)
    expect_equal(m$threshold, threshold)
    expect_equal(m$molecules, molecules)
  })
}

# ***************************************************************************************
# Method of Moments
# ***************************************************************************************
# Molecules: 1, Threshold: 1
test_est(0.01, 18, N, 1, 1, 'mom')
test_est(0.05, 17, N, 1, 1, 'mom')
test_est(0.1 , 16, N, 1, 1, 'mom')
test_est(0.2 , 15, N, 1, 1, 'mom')
test_est(0.3 , 14, N, 1, 1, 'mom')
test_est(0.4 , 13, N, 1, 1, 'mom')
test_est(0.5 , 12, N, 1, 1, 'mom')
test_est(0.6 , 11, N, 1, 1, 'mom')
test_est(0.7 , 10, N, 1, 1, 'mom')
test_est(0.8 , 9, N, 1, 1, 'mom')
test_est(0.9 , 8, N, 1, 1, 'mom')
test_est(0.95, 7, N, 1, 1, 'mom')
test_est(0.99, 6, N, 1, 1, 'mom')
test_est(1.0 , 5, N, 1, 1, 'mom')

# Molecules: 1, Threshold > 1
test_est(0.3 , 14, N, 1, 5, 'mom')
test_est(0.5 , 12, N, 1, 5, 'mom')
test_est(0.9 , 8, N, 1, 4, 'mom')
test_est(1.0 , 5, N, 1, 4, 'mom')

# Molecules: 2, Threshold > 1
test_est(0.5 , 12, N, 2, 5, 'mom')
test_est(1.0 , 5, N, 2, 4, 'mom')

# Molecules: 3, Threshold > 1
test_est(0.5 , 12, N, 3, 5, 'mom')
test_est(1.0 , 5, N, 3, 4, 'mom')

# ***************************************************************************************
# Maximum Likelihood
# ***************************************************************************************
test_est(0.01, 18, N, 1, 1, 'mle')
test_est(0.05, 17, N, 1, 1, 'mle')
test_est(0.1 , 16, N, 1, 1, 'mle')
test_est(0.2 , 15, N, 1, 1, 'mle')
test_est(0.3 , 14, N, 1, 1, 'mle')
test_est(0.4 , 13, N, 1, 1, 'mle')
test_est(0.5 , 12, N, 1, 1, 'mle')
test_est(0.6 , 11, N, 1, 1, 'mle')
test_est(0.7 , 10, N, 1, 1, 'mle')
test_est(0.8 , 9, N, 1, 1, 'mle')
test_est(0.9 , 8, N, 1, 1, 'mle')
test_est(0.95, 7, N, 1, 1, 'mle')
test_est(0.99, 6, N, 1, 1, 'mle')
test_est(1.0 , 5, N, 1, 1, 'mle')

# Molecules: 1, Threshold > 1
test_est(0.3 , 14, N, 1, 5, 'mle')
test_est(0.5 , 12, N, 1, 5, 'mle')
test_est(0.9 , 8, N, 1, 4, 'mle')
test_est(1.0 , 5, N, 1, 4, 'mle')

# Molecules: 2, Threshold > 1
test_est(0.5 , 12, N, 2, 5, 'mle')
test_est(1.0 , 5, N, 2, 2, 'mle')

# Molecules: 3, Threshold > 1
test_est(0.5 , 12, N, 3, 5, 'mle')
test_est(1.0 , 5, N, 3, 2, 'mle')
