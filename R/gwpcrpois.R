#' @title Poissonian Sampling Distribution of the PCR Product Distribution
#'
#' @description Test
#'
#' @param n
#'
#' @param lambda
#'
#' @param efficiency
#'
#' @param molecules
#'
#' @name gwpcrpois
#'
#' @seealso gwpcr
NULL

#' @rdname gwpcrpois
#' @useDynLib gwpcR gwpcrpois_simulate_c
#' @export
rgwpcrpois <- function(samples, efficiency, lambda0, threshold=1, molecules=1, cycles=NA) {
  if (!is.numeric(efficiency) || (length(efficiency) != 1) || (efficiency <= 0) || (efficiency > 1))
    stop('efficiency must be a numeric scalar within (0,1]')
  if (!is.numeric(lambda0) || (length(lambda0) != 1) || (lambda0 <= 0))
    stop('lambda0 must be a numeric scalar within (0,Infinity)')
  if (!is.numeric(threshold) || (length(threshold) != 1) || (threshold != floor(threshold)) || (threshold < 0))
    stop('threshold must be a non-negative integral scalar')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  # Determine how many cycles are necessary on average to produce 1e6
  # molecules. After that point, we assume that the additional variability
  # is negligible.
  if (is.na(cycles))
    cycles <- ceiling(log(1e6 / molecules) / log(1+efficiency))

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
    # We add 3 standard deviations to make it relatively unlikely that we'll need
    # more than one iteration.
    k <- ceiling((samples - n) / p.th)
    k <- ceiling(k + 3*sqrt(k*p.th*(1-p.th)))
    # Generate lambda values for samples by sampling from the PCR distribution
    s.c <- .C(gwpcrpois_simulate_c,
              nsamples=as.integer(k),
              samples=double(k),
              efficiency=as.double(efficiency),
              lambda0=as.double(lambda0),
              molecules=as.double(molecules),
              ncycles=as.integer(cycles),
              NAOK=TRUE)$samples
    # Remove samples below threshold, and determine how many to use
    s.c <- s.c[s.c >= threshold]
    if (length(s.c) == 0)
      next
    c.n <- min(length(s.c), samples - n)
    # Append to output
    stopifnot(c.n > 0)
    r[(n+1):(n+c.n)] <- head(s.c, c.n)
    n <- n + c.n
  }

  r
}

#' @rdname gwpcrpois
#' @export
dgwpcrpois <- function(c, efficiency, lambda0, threshold=1, molecules=1) {
  handle.parameters(list(c=c, efficiency=efficiency, lambda0=lambda0,
                         threshold=threshold, molecules=molecules),
                    by=c('efficiency', 'lambda0', 'threshold', 'molecules'), {
                      if (is.finite(lambda0) && (lambda0 > 0)) {
                        # To ensure reasonable accuracy, we want the effective lambda
                        # grid we use to compute the mixture (i.e. AFTER scaling with lambda0)
                        # to be much finer than the variance of the poisson distribution.
                        # We thus demand that:
                        #    lambda0 * width <= (1/10) * Sqrt(l * lambda0),
                        # i.e. that
                        #    width <= Sqrt(l / lambda0) / 10
                        # At l=0, we instead chose width such that we get 10 samples of l*lambda0
                        # within [0, 1].
                        grid.width.fun <- function(l) { ifelse(l > 0,
                                                               sqrt(l / lambda0) / 10,
                                                               1 / (10 * lambda0)) }

                        # Compute P[X < Th]
                        p <- if (threshold >= 1)
                          gwpcr.mixture(threshold-1, function(x,l) { stats::ppois(x,l*lambda0) },
                                        efficiency=efficiency, molecules=molecules,
                                        grid.width.fun=grid.width.fun)
                        else
                          0
                        # Compute probabilities for those X which are >= TH
                        c.v <- (c >= threshold)
                        d <- rep(0, length(c))
                        if (sum(c.v) > 0)
                          d[c.v] <- gwpcr.mixture(c[c.v], function(x,l) { stats::dpois(x,l*lambda0) },
                                                  efficiency=efficiency, molecules=molecules,
                                                  grid.width.fun=grid.width.fun)
                        # And compute P(X = c | X >= Th)
                        d / (1.0 - p)
                      } else if ((lambda0 == 0.0) && (threshold == 0)) {
                        ifelse(c == 0, 1.0, 0.0)
                      } else {
                        rep(as.numeric(NA), length(c))
                      }
                    })
}

#' @rdname gwpcrpois
#' @export
pgwpcrpois <- function(c, efficiency, lambda0, threshold=1, molecules=1) {
  handle.parameters(list(c=c, efficiency=efficiency, lambda0=lambda0,
                         threshold=threshold, molecules=molecules),
                    by=c('efficiency', 'lambda0', 'threshold', 'molecules'), {
                      if (is.finite(lambda0) && (lambda0 > 0)) {
                        # To ensure reasonable accuracy, we want the effective lambda
                        # grid we use to compute the mixture (i.e. AFTER scaling with lambda0)
                        # to be much finer than the variance of the poisson distribution.
                        # We thus demand that:
                        #    lambda0 * width <= (1/10) * Sqrt(l * lambda0),
                        # i.e. that
                        #    width <= Sqrt(l / lambda0) / 10
                        # At l=0, we instead chose width such that we get 10 samples of l*lambda0
                        # within [0, 1].
                        grid.width.fun <- function(l) { ifelse(l > 0,
                                                               sqrt(l / lambda0) / 10,
                                                               1 / (10 * lambda0)) }

                        # Compute P[X < Th]
                        p <- if (threshold > 0)
                          gwpcr.mixture(threshold-1, function(x,l) { stats::ppois(x,l*lambda0) },
                                        efficiency=efficiency, molecules=molecules,
                                        grid.width.fun=grid.width.fun)
                        else
                          0
                        # Compute probabilities for those X which are >= TH
                        c.v <- (c >= threshold)
                        d <- rep(0, length(c))
                        if (sum(c.v) > 0)
                          d[c.v] <- gwpcr.mixture(c[c.v], function(x,l) { stats::ppois(x,l*lambda0) },
                                                  efficiency=efficiency, molecules=molecules,
                                                  grid.width.fun=grid.width.fun)
                        # And compute P(X <= c | X >= Th)
                        (d - p) / (1.0 - p)
                      } else if ((lambda0 == 0.0) && (threshold == 0)) {
                        ifelse(c == 0, 1.0, 0.0)
                      } else {
                        rep(as.numeric(NA), length(c))
                      }
                    })
}
