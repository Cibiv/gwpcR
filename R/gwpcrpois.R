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
#' @useDynLib gwpcR gwpcrpois_simulate
#' @export
rgwpcrpois <- function(samples, efficiency, lambda0, threshold=1, molecules=1, cycles=NA) {
  # Determine how many cycles are necessary on average to produce 1e6
  # molecules. After that point, we assume that the additional variability
  # is negligible.
  if (is.na(cycles))
    cycles <- ceiling(log(1e6 / molecules)/log(1+efficiency))

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
    k <- ceiling(k + 3*sqrt(k*p.th*(1-p.th)))
    # Generate lambda values for samples by sampling from the PCR distribution
    s.c <- .C(gwpcrpois_simulate,
              nsamples=as.integer(k),
              samples=double(k),
              efficiency=as.double(efficiency),
              lambda0=as.double(lambda0),
              molecules=as.double(molecules),
              cycles=cycles,
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

#' @rdname gwpcrpois
#' @export
pgwpcrpois <- function(c, efficiency, lambda0, threshold=1, molecules=1) {
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
