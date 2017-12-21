library(gwpcR)
library(data.table)
context('Galton-Watson PCR (gwpcr) distribution functions')

P.TH <- 0.1
TRIALS <- 6

test_gwpcr <- function(efficiency, molecules, N) {
  test_that(paste0("rgwpcr, pgwpcr (efficiency=", 100*efficiency, "%)"), {
    i <- 0
    p <- 0
    while((p <= P.TH) && (i < TRIALS)) {
      i <- i + 1
      x <- rgwpcr(n=N, efficiency=efficiency, molecules=molecules)
      p <- ks.test(x=x, y="pgwpcr", efficiency=efficiency, molecules=molecules)$p.value
    }
    expect_gt(p, P.TH)
  })
}

test_gwpcr(0.00, 1, 1e4)
test_gwpcr(0.01, 1, 1e4)
test_gwpcr(0.05, 1, 1e4)
test_gwpcr(0.10, 1, 1e4)
test_gwpcr(0.20, 1, 1e4)
test_gwpcr(0.40, 1, 1e4)
test_gwpcr(0.80, 1, 1e4)
test_gwpcr(0.95, 1, 1e4)
test_gwpcr(0.99, 1, 1e4)

test_gwpcr(0.00, 2, 1e4)
test_gwpcr(0.01, 2, 1e4)
test_gwpcr(0.05, 2, 1e4)
test_gwpcr(0.10, 2, 1e4)
test_gwpcr(0.20, 2, 1e4)
test_gwpcr(0.40, 2, 1e4)
test_gwpcr(0.80, 2, 1e4)
test_gwpcr(0.95, 2, 1e4)
test_gwpcr(0.99, 2, 1e4)

test_gwpcr(0.00, 5, 1e4)
test_gwpcr(0.01, 5, 1e4)
test_gwpcr(0.05, 5, 1e4)
test_gwpcr(0.10, 5, 1e4)
test_gwpcr(0.20, 5, 1e4)
test_gwpcr(0.40, 5, 1e4)
test_gwpcr(0.80, 5, 1e4)
test_gwpcr(0.95, 5, 1e4)
test_gwpcr(0.99, 5, 1e4)

test_that("gwpcr (efficiency=0%, finite cycle count)", {
  expect_identical(rgwpcr(n=5, efficiency=0, molecules=1, cycles=0), c(1,1,1,1,1))
  expect_identical(rgwpcr(n=5, efficiency=0, molecules=1, cycles=1), c(1,1,1,1,1))
  expect_identical(rgwpcr(n=5, efficiency=0, molecules=1, cycles=2), c(1,1,1,1,1))
})

test_that("gwpcr (efficiency=1e-6, finite cycle count)", {
  # Will fail with probability =~ 5e-6.
  expect_identical(rgwpcr(n=5, efficiency=1e-6, molecules=1, cycles=0), c(1,1,1,1,1))
  expect_identical(rgwpcr(n=5, efficiency=1e-6, molecules=1, cycles=1), c(1,1,1,1,1)/(1+1e-6))
  expect_identical(rgwpcr(n=5, efficiency=1e-6, molecules=1, cycles=2), c(1,1,1,1,1)/(1+1e-6)^2)
})

test_that("pgwpcr", {
  # XXX: Shouldn't this raise an error instead?
  expect_identical(pgwpcr(-1.0, efficiency=-1), as.numeric(NA))
})

test_that("pgwpcr.fun", {
  f <- dgwpcr.fun(efficiency=0.5, molecules=1)
  expect_equal(f(-0.5), 0.0)
  expect_equal(f(0), 0.0)
  expect_equal(f(max(gwpcR:::GWPCR$lambda)), 0.0)
})
