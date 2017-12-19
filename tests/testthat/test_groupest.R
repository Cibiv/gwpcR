library(gwpcR)
library(data.table)
context('Groupwise Parameter Estimation')

EPS.REL.VAR <- 25e-2
ESP.REL.EFFICIENCY <- 10e-2
ESP.REL.LAMBDA0 <- 5e-2

test_est_grp <- function(efficiency, efficiency.sd, lambda0, lambda0.sd) {
  test_that(paste0("group estimation (",
                   "efficiency: ", 100*efficiency, '% +/- ', 100*efficiency.sd, "%",
                   "lambda0: ", lambda0, "% +/- ", lambda0.sd), {
    n.obs <- rep(c(5, 15, 50, 150, 500), times=c(200, 200, 200, 200, 200))
    efficiency <- 0.5
    efficiency.sd <- 0.1
    lambda0 <- 5
    lambda0.sd <- 1

    g.in <- data.table(g=1:length(n.obs), n.obs=n.obs)
    g.in[, efficiency.true := rbeta(n=.N, shape1=efficiency*((efficiency*(1-efficiency))/efficiency.sd^2 - 1),
                                    shape2=(1-efficiency)*((efficiency*(1-efficiency))/efficiency.sd^2 - 1))]
    g.in[, lambda0.true := rnorm(n=.N, mean=lambda0, sd=lambda0.sd)]
    setkey(g.in, g)

    g.obs <- g.in[, list(count=rgwpcrpois(n=n.obs, threshold=0, molecules=1,
                                          efficiency=efficiency.true,
                                          lambda0=lambda0.true))
                  , by=g]
    g.est <- gwpcrpois.groupest(count ~ g, g.obs, threshold=0, molecules=1)
    g.est[, efficiency.true := g.in[g.est, efficiency.true]]
    g.est[, efficiency.true := g.in[g.est, efficiency.true]]
    g.est[, lambda0.true := g.in[g.est, lambda0.true]]

    # Check that the between-group standard deviations of efficiency and lambda0 are approximately recovered
    expect_lt(abs(sqrt(g.est$efficiency.grp.var[1]) - efficiency.sd) / efficiency.sd, EPS.REL.VAR)
    expect_lt(abs(sqrt(g.est$lambda0.grp.var[1]) - lambda0.sd) / lambda0.sd, EPS.REL.VAR)

    expect_lt(g.est[n.obs==500, mean(abs(efficiency - efficiency.true) / efficiency.true)], ESP.REL.EFFICIENCY)
    expect_lt(g.est[n.obs==500, mean(abs(lambda0 - lambda0.true) / lambda0)], ESP.REL.LAMBDA0)
  })
}

test_est_grp(0.2, 0.05,  5, 1)
test_est_grp(0.2, 0.05, 15, 3)
test_est_grp(0.5, 0.1 ,  5, 1)
test_est_grp(0.5, 0.1 , 15, 3)
test_est_grp(0.8, 0.1 ,  5, 1)
test_est_grp(0.8, 0.1 , 15, 3)
