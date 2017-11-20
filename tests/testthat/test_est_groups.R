library(gwpcR)
library(data.table)
context('Groupwise Parameter Estimation')

EPS.REL <- 2e-1

test_that("between-group variance estimates", {
  n.obs <- rep(c(5, 15, 50, 150, 500), times=c(50, 50, 50, 50, 50))
  efficiency <- 0.5
  efficiency.sd <- 0.1
  lambda0 <- 5
  lambda0.sd <- 1

  g.in <- data.table(g=1:length(n.obs), n.obs=n.obs)
  g.in[, efficiency.true := rbeta(n=.N, shape1=efficiency*((efficiency*(1-efficiency))/efficiency.sd^2 - 1),
                                  shape2=(1-efficiency)*((efficiency*(1-efficiency))/efficiency.sd^2 - 1))]
  g.in[, lambda0.true := rnorm(n=.N, mean=lambda0, sd=lambda0.sd)]

  g.obs <- g.in[, list(count=rgwpcrpois(n=n.obs, threshold=0, molecules=1,
                                        efficiency=efficiency.true,
                                        lambda0=lambda0.true))
                , by=g]
  g.est <- gwpcrpois.est.groups(count ~ g, g.obs, threshold=0, molecules=1)

  expect_lt(abs(sqrt(g.est$efficiency.grp.var[1]) - efficiency.sd) / efficiency.sd, EPS.REL)
  expect_lt(abs(sqrt(g.est$lambda0.grp.var[1]) - lambda0.sd) / lambda0.sd, EPS.REL)
})
