# **********************************************************************
# Mixture of Poisson Distributions with Parameter l * l0, where
# l is distributed according to the Galton-Watson PCR Distribution 
# **********************************************************************

rgwpcrpois <- function(samples, efficiency, lambda0, threshold=1, molecules=1) {
  # Determine probability of a sample being accepted (i.e., of being >= threshold)
  p.th <- if (threshold > 0)
    1.0 - pgwpcrpois(threshold-1, efficiency=efficiency, lambda0=lambda0,
                     threshold=0, molecules=molecules)
  else
    1.0

  # Sample until we have enough samples
  n <- 0
  r <- rep(as.integer(NA), samples)
  while(n < samples) {
    # Compute number of samples generate, taking the average acception rate p.th
    # into account, as well as the fluctuations around that rate, which are binomial.
    k <- ceiling((samples - n) / p.th)
    k <- k + 2*sqrt(k*p.th*(1-p.th))
    # Generate lambda values for samples by sampling from the PCR distribution
    s.l <- rgwpcr(samples=k + 2*sqrt(k), efficiency=efficiency, molecules=molecules)
    # Generate counts by sampling from a (different) poissin distribution per sample
    s.c <- sapply(s.l, function(l) { stats::rpois(n=1, lambda=l*lambda0) })
    # Remove samples below threshold, and determine how many to use
    s.c <- s.c[s.c >= threshold]
    c.n <- min(length(s.c), samples - n)
    # Append to output
    stopifnot(c.n > 0)
    r[(n+1):(n+c.n)] <- head(s.c, c.n)
    n <- n + c.n
  }

  r
}

dgwpcrpois <- function(c, efficiency, lambda0, threshold=1, molecules=1, clamp.efficiency=TRUE) {
  if (clamp.efficiency)
    efficiency <- pmin(pmax(GWPCR$efficiency[1], efficiency), tail(GWPCR$efficiency,1))

  handle.parameters(list(c=c, efficiency=efficiency, lambda0=lambda0,
                         threshold=threshold, molecules=molecules),
                    by=c('efficiency', 'lambda0', 'threshold', 'molecules'), {
                      # Compute P[X < Th]
                      p <- if (threshold >= 1)
                        gwpcr.mixture(threshold-1, stats::ppois, efficiency=efficiency,
                                      lambda0=lambda0, molecules=molecules)
                      else
                        0
                      # Compute probabilities for those X which are >= TH
                      c.v <- (c >= threshold)
                      d <- rep(0, length(c))
                      d[c.v] <- gwpcr.mixture(c[c.v], stats::dpois, efficiency=efficiency,
                                              lambda0=lambda0, molecules=molecules)
                      # And compute P(X = c | X >= Th)
                      d / (1.0 - p)
                    })
}

pgwpcrpois <- function(c, efficiency, lambda0, threshold=1, molecules=1, clamp.efficiency=TRUE) {
  if (clamp.efficiency)
    efficiency <- pmin(pmax(GWPCR$efficiency[1], efficiency), tail(GWPCR$efficiency,1))

    handle.parameters(list(c=c, efficiency=efficiency, lambda0=lambda0,
                         threshold=threshold, molecules=molecules),
                    by=c('efficiency', 'lambda0', 'threshold', 'molecules'), {
                      # Compute P[X < Th]
                      p <- if (threshold > 0)
                        gwpcr.mixture(threshold-1, stats::ppois, efficiency=efficiency,
                                      lambda0=lambda0, molecules=molecules)
                      else
                        0
                      # Compute probabilities for those X which are >= TH
                      c.v <- (c >= threshold)
                      d <- rep(0, length(c))
                      d[c.v] <- gwpcr.mixture(c[c.v], stats::ppois, efficiency=efficiency,
                                              lambda0=lambda0, molecules=molecules)
                      # And compute P(X <= c | X >= Th)
                      (d - p) / (1.0 - p)
                    })
}
