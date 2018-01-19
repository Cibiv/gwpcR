library(gwpcR)
library(data.table)
context('Galton-Watson PCR (gwpcrpois) distribution functions')

P.TH <- 0.1
TRIALS <- 6

test_gwpcrpois <- function(efficiency, lambda0, threshold, molecules, N) {
  test_that(paste0("rgwpcrpois, pgwpcrpois (efficiency=", 100*efficiency,
                   "%. lambda0=", lambda0 , ")"), {
    i <- 0
    p <- 0
    while((p <= P.TH) && (i < TRIALS)) {
      i <- i + 1
      x <- rgwpcrpois(n=N, efficiency=efficiency, lambda0=lambda0, molecules=molecules, threshold=threshold)
      x.uniq <- sort(unique(x))
      cdf <- stepfun(x=c(x.uniq, x.uniq[length(x.uniq)]+1),
                     y=c(0, pgwpcrpois(x.uniq, efficiency=efficiency, lambda0=lambda0,
                                       threshold=threshold, molecules=molecules), 1.0))
      # Must use K-S test from dgof, because the distribution is discrete!
      p <- dgof::ks.test(x, cdf, exact=FALSE)$p.value
    }
    expect_gt(p, P.TH)
  })
}

test_gwpcrpois(0.00, 5  , 1, 1, 1e4)
test_gwpcrpois(0.01, 50 , 0, 1, 1e4)
test_gwpcrpois(0.05, 500, 1, 1, 1e4)
test_gwpcrpois(0.10, 4  , 0, 1, 1e4)
test_gwpcrpois(0.20, 40 , 1, 1, 1e4)
test_gwpcrpois(0.40, 400, 0, 1, 1e4)
test_gwpcrpois(0.80, 3  , 1, 1, 1e4)
test_gwpcrpois(0.95, 30 , 0, 1, 1e4)
test_gwpcrpois(0.99, 300, 1, 1, 1e4)

test_gwpcrpois(0.00, 5  , 0, 2, 1e4)
test_gwpcrpois(0.01, 50 , 1, 2, 1e4)
test_gwpcrpois(0.05, 500, 0, 2, 1e4)
test_gwpcrpois(0.10, 4  , 1, 2, 1e4)
test_gwpcrpois(0.20, 40 , 0, 2, 1e4)
test_gwpcrpois(0.40, 400, 1, 2, 1e4)
test_gwpcrpois(0.80, 3  , 0, 2, 1e4)
test_gwpcrpois(0.95, 30 , 1, 2, 1e4)
test_gwpcrpois(0.99, 300, 0, 2, 1e4)

test_gwpcrpois(0.00, 5  , 2, 5, 1e4)
test_gwpcrpois(0.01, 50 , 1, 5, 1e4)
test_gwpcrpois(0.05, 500, 2, 5, 1e4)
test_gwpcrpois(0.10, 4  , 1, 5, 1e4)
test_gwpcrpois(0.20, 40 , 2, 5, 1e4)
test_gwpcrpois(0.40, 400, 1, 5, 1e4)
test_gwpcrpois(0.80, 3  , 2, 5, 1e4)
test_gwpcrpois(0.95, 30 , 1, 5, 1e4)
test_gwpcrpois(0.99, 300, 2, 5, 1e4)

test_gwpcrpois_1cycle <- function(efficiency, lambda0, N) {
  test_that(paste0("rgwpcrpois, pgwpcrpois (efficiency=", 100*efficiency,
                   "%, lambda0=", lambda0 , ", cycles=1)"), {
                     i <- 0
                     p <- 0
                     while((p <= P.TH) && (i < TRIALS)) {
                       i <- i + 1
                       x <- rgwpcrpois(n=N, efficiency=efficiency, lambda0=lambda0, molecules=1, threshold=0, cycles=1)
                       x.uniq <- sort(unique(x))
                       cdf <- stepfun(x=c(x.uniq, x.uniq[length(x.uniq)]+1),
                                      y=c(0, (ppois(x.uniq, lambda=lambda0     / (1+efficiency)) * (1 - efficiency) +
                                              ppois(x.uniq, lambda=lambda0 * 2 / (1+efficiency)) *      efficiency  ), 1.0))
                       # Must use K-S test from dgof, because the distribution is discrete!
                       p <- dgof::ks.test(x, cdf, exact=FALSE)$p.value
                     }
                     expect_gt(p, P.TH)
                   })
}

test_gwpcrpois_1cycle(0.01, 5, 1e4)
test_gwpcrpois_1cycle(0.05, 5, 1e4)
test_gwpcrpois_1cycle(0.5 , 5, 1e4)
test_gwpcrpois_1cycle(0.99, 5, 1e4)
