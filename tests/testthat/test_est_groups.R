library(gwpcR)
library(data.table)
context('Groupwise Parameter Estimation')

EPS.ABS <- 5e-2
EPS.REL <- 5e-2

test_that("test", {
  n.obs <- rep(c(5, 50, 500, 5000), times=c(2000, 200, 200  ,200))
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
  print(head(g.obs))
  g.est <- gwpcrpois.est.groups(count ~ g, g.obs, threshold=0, molecules=1)
  View(g.est)
})