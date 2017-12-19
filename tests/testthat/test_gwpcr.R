library(gwpcR)
library(data.table)
context('Galton-Watson PCR Distribution (gwpcR)')

P.TH <- 0.1
TRIALS <- 6

test_gwpcr <- function(efficiency, molecules, N) {
  test_that(paste0("gwpcr(efficiency=", 100*efficiency, "%)"), {
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
