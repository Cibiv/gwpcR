#' @title Poissonian Sampling Distribution of the PCR Product Distribution
#'
#' @description Write Me
#'
#' @param n number of random samples to generate
#'
#' @param c number of observations of a particular molecular family
#'
#' @param lambda0 average number of observations per molecular family
#'
#' @param efficiency efficiency of amplification
#'
#' @param molecules initial copy number
#'
#' @param threshold minimal number of observations a molecular family must have
#'   to count as \dfn{unambiguously detected}. Setting this to a value \eqn{v >=
#'   0} \emph{conditions} the distribution on \eqn{c >= c}, i.e every value of
#'   \eqn{c} less than that gets assigned probability zero.
#'
#' @param cycles number of amplification cycles used for simulation. By default,
#'   a large enough value is used to make the results virtually idistinguishable
#'   from the limit for \eqn{cycles \to \infty}
#'
#' @name gwpcrpois
#'
#' @seealso \code{\link{gwpcr}}
NULL

#' @rdname gwpcrpois
#' @useDynLib gwpcR gwpcrpois_simulate_c
#' @export
rgwpcrpois <- function(n, efficiency, lambda0, threshold=1, molecules=1, cycles=Inf) {
  if (!is.numeric(efficiency) || (length(efficiency) != 1) || (efficiency < 0) || (efficiency > 1))
    stop('efficiency must be a numeric scalar within [0,1]')
  if (!is.numeric(lambda0) || (length(lambda0) != 1) || (lambda0 <= 0))
    stop('lambda0 must be a numeric scalar within (0,Infinity)')
  if (!is.numeric(threshold) || (length(threshold) != 1) || (threshold != floor(threshold)) || (threshold < 0))
    stop('threshold must be a non-negative integral scalar')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')
  if (!is.numeric(cycles) || (length(cycles) != 1) || (cycles != floor(cycles)) || (cycles < 0))
    stop('cycles must be a positive integral scalar or +Infinity')
  
  # Determine a suitable cycle count if set to "Infinity", which we take
  # to mean "as many as necessary so that the results are virtually
  # indistinguishable from the limit case
  method <- if (is.infinite(cycles) && (efficiency >= E.MIN)) {
    # Determine how many cycles are necessary on average to produce 1e6
    # molecules. After that point, we assume that the additional variability
    # is negligible.
    cycles <- ceiling(log(1e6 / molecules) / log(1+efficiency))
    "simulate"
  } else if (is.infinite(cycles) && (efficiency < E.MIN)) {
    # Instead of simulating the PCR process, generate samples using the
    # gamma approximation of the PCR distribution
    "gamma"
  } else if (is.finite(cycles)) {
    # For finitely many cycles, always simulate
    "simulate"
  }
  else {
    # Should never happen
    NA
  }

  # Determine probability of a sample being accepted (i.e., of being >= threshold)
  p.th <- if (threshold > 0)
    1.0 - pgwpcrpois(threshold-1, efficiency=efficiency, lambda0=lambda0,
                     threshold=0, molecules=molecules)
  else
    1.0

  # Sample until we have enough samples
  j <- 0
  r <- rep(as.double(NA), j)
  while(j < n) {
    # Compute number of samples generate, taking the average acception rate p.th
    # into account, as well as the fluctuations around that rate, which are binomial.
    # We add 3 standard deviations to make it relatively unlikely that we'll need
    # more than one iteration.
    k <- ceiling((n - j) / p.th)
    k <- ceiling(k + 3*sqrt(k*p.th*(1-p.th)))
    
    s.c <- if (method == "simulate") {
      # Generate read counts by sampling from the PCR-Poisson distribution
      .C(gwpcrpois_simulate_c,
         nsamples=as.integer(k),
         samples=double(k),
         efficiency=as.double(efficiency),
         lambda0=as.double(lambda0),
         molecules=as.double(molecules),
         ncycles=as.integer(cycles),
         NAOK=TRUE)$samples
    } else if (method == "gamma") {
      # Generate read counts by replacing the PCR distribution with its gamma
      # approximation. Note that mixing Poisson distribution with Gamma rates
      # yields a Negative Binomial distribution with parameters (see also gwpcr.sd)
      nb.size <- 1/gwpcr.sd(efficiency, molecules)^2
      nb.prob <- nb.size / (nb.size + lambda0)
      rnbinom(n, size=nb.size, prob=nb.prob)
    } else
      stop(paste0("Unknown method ", method))

    # Remove samples below threshold, and determine how many to use
    s.c <- s.c[s.c >= threshold]
    if (length(s.c) == 0)
      next
    c <- min(length(s.c), n - j)
    # Append to output
    stopifnot(c > 0)
    r[(j+1):(j+c)] <- head(s.c, c)
    j <- j + c
  }

  r
}

#' @rdname gwpcrpois
#' @export
dgwpcrpois <- function(c, efficiency, lambda0, threshold=1, molecules=1) {
  handle.parameters(list(c=c, efficiency=efficiency, lambda0=lambda0,
                         threshold=threshold, molecules=molecules),
                    by=c('efficiency', 'lambda0', 'threshold', 'molecules'), {
                      r <- rep(as.numeric(NA), length(c))

                      r <- if (!is.finite(efficiency) || (efficiency < 0) || (efficiency > 1))
                        rep(as.numeric(NA), length(c))
                      else if (!is.finite(lambda0) || (lambda0 < 0))
                        rep(as.numeric(NA), length(c))
                      else if (!is.finite(molecules) || (molecules != floor(molecules)) || (molecules < 1))
                        rep(as.numeric(NA), length(c))
                      else if (!is.finite(threshold) || (threshold != floor(threshold)) || (threshold < 0))
                        rep(as.numeric(NA), length(c))
                      else if ((lambda0 == 0) && (threshold == 0))
                        ifelse(c == 0, 1.0, 0.0)
                      else if ((lambda0 == 0) && (threshold > 0))
                        rep(as.numeric(NA), length(c))
                      else
                        rep(0.0, length(c))

                      if ((lambda0 > 0) && (efficiency >= E.MIN)) {
                        # Within precomputed lambda range
                        #
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
                        # And compute P(X = c | X >= Th) and scale appropriately,
                        # (in case we're within the gamma cutover range).
                        r <- r + (1 - gamma.factor(efficiency)) * d / (1.0 - p)
                      }

                      if ((lambda0 > 0) && (efficiency <= E.GAMMA.TH)) {
                        # Within gamma approximation range
                        #
                        # Mixing Poissons with Gamma rates yields a Negative Binomial,
                        # in our case the parameters are (see also gwpcr.sd()):
                        nb.size <- 1/gwpcr.sd(efficiency, molecules)^2
                        nb.prob <- nb.size / (nb.size + lambda0)

                        # Compute P[X < Th]
                        p <- if (threshold >= 1)
                          pnbinom(threshold-1, size=nb.size, prob=nb.prob)
                        else
                          0
                        # Compute probabilities for those X which are >= TH
                        c.v <- (c >= threshold)
                        d <- rep(0, length(c))
                        if (sum(c.v) > 0)
                          d[c.v] <- dnbinom(c[c.v], size=nb.size, prob=nb.prob)
                        # And compute P(X = c | X >= Th) and scale appropriately,
                        # (in case we're within the gamma cutover range).
                        r <- r + gamma.factor(efficiency) * d / (1.0 - p)
                      }

                      r
                    })
}

#' @rdname gwpcrpois
#' @export
pgwpcrpois <- function(c, efficiency, lambda0, threshold=1, molecules=1) {
  handle.parameters(list(c=c, efficiency=efficiency, lambda0=lambda0,
                         threshold=threshold, molecules=molecules),
                    by=c('efficiency', 'lambda0', 'threshold', 'molecules'), {
                      r <- rep(as.numeric(NA), length(c))

                      r <- if (!is.finite(efficiency) || (efficiency < 0) || (efficiency > 1))
                        rep(as.numeric(NA), length(c))
                      else if (!is.finite(lambda0) || (lambda0 < 0))
                        rep(as.numeric(NA), length(c))
                      else if (!is.finite(molecules) || (molecules != floor(molecules)) || (molecules < 1))
                        rep(as.numeric(NA), length(c))
                      else if (!is.finite(threshold) || (threshold != floor(threshold)) || (threshold < 0))
                        rep(as.numeric(NA), length(c))
                      else if ((lambda0 == 0) && (threshold == 0))
                        ifelse(c >= 0, 1.0, 0.0)
                      else if ((lambda0 == 0) && (threshold > 0))
                        rep(as.numeric(NA), length(c))
                      else
                        rep(0.0, length(c))

                      if ((lambda0 > 0) && (efficiency >= E.MIN)) {
                        # Within precomputed lambda range
                        #
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
                          d[c.v] <- gwpcr.mixture(c[c.v], function(x,l) { stats::ppois(x,l*lambda0) },
                                                  efficiency=efficiency, molecules=molecules,
                                                  grid.width.fun=grid.width.fun) - p
                        # And compute P(X = c | X >= Th) and scale appropriately,
                        # (in case we're within the gamma cutover range).
                        r <- r + (1 - gamma.factor(efficiency)) * d / (1.0 - p)
                      }

                      if ((lambda0 > 0) && (efficiency <= E.GAMMA.TH)) {
                        # Within gamma approximation range
                        #
                        # Mixing Poissons with Gamma rates yields a Negative Binomial,
                        # in our case the parameters are (see also gwpcr.sd()):
                        nb.size <- 1/gwpcr.sd(efficiency, molecules)^2
                        nb.prob <- nb.size / (nb.size + lambda0)

                        # Compute P[X < Th]
                        p <- if (threshold >= 1)
                          pnbinom(threshold-1, size=nb.size, prob=nb.prob)
                        else
                          0
                        # Compute probabilities for those X which are >= TH
                        c.v <- (c >= threshold)
                        d <- rep(0, length(c))
                        if (sum(c.v) > 0)
                          d[c.v] <- pnbinom(c[c.v], size=nb.size, prob=nb.prob) - p
                        # And compute P(X = c | X >= Th) and scale appropriately,
                        # (in case we're within the gamma cutover range).
                        r <- r + gamma.factor(efficiency) * d / (1.0 - p)
                      }

                      r
                    })
}
