# **********************************************************************
#    PCR Distribution as produced by a simple Galton-Watson process
# **********************************************************************

rgwpcr <- function(samples, efficiency, molecules=1) {
  if (!is.numeric(efficiency) || (length(efficiency) != 1) || (efficiency < 0) || (efficiency > 1))
    stop('efficiency must be a numeric scalar within [0,1]')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  # Determine how many cycles are necessary on average to produce 1e6
  # molecules. After that point, we assume that the additional variability
  # is negligible.
  cycles <- ceiling(log(1e6 / molecules)/log(1+efficiency))

  # Start with initial copy number
  s <- rep(as.double(molecules), samples)
  if ((efficiency > 0) && (efficiency < 1.0)) {
    # Flip a coin for each molecule in each sample to determine whether
    # its copied. That can be done efficiently by sampling from a binomial
    # distribution for each sample.
    # Note: Since this requires looping over all samples, the code was
    # translated to C. See gwpcr_simulate in simulate.c.
    if (TRUE)
      s <- .C(gwpcr_simulate,
              nsamples=as.integer(length(s)),
              samples=as.double(s),
              efficiency=as.double(efficiency),
              cycles=as.integer(cycles),
              NAOK=TRUE)$samples
    else
      for(i in 1:cycles) {
        sp <- sapply(s, FUN=function(z) { rbinom(n=1, size=z, prob=efficiency) })
        s <- s + sp
      }

    # Scale with expected number of molecules
    s <- s / (molecules * (1+efficiency)**cycles)
  }

  # Return samples
  s
}

dgwpcr <- function(lambda, efficiency, molecules=1) {
  # Process each number of molecules separately
  handle.parameters(list(lambda=lambda, efficiency=efficiency, molecules=molecules), by="molecules", {
    # Moleculec ount must be >= 1 and integral
    if (is.finite(molecules) && (molecules == floor(molecules)) && (molecules > 0)) {
      # Determine molecule count filter (p.m), efficiency and lambda vectors
      # (e, l), their validity filter (p.e.v, p.l.v) and overall filter (p.v)
      e <- efficiency
      l <- lambda
      p.e.v <- (e >= head(GWPCR$efficiency, 1)) & (e <= tail(GWPCR$efficiency, 1))
      p.l.v <- (l >= 0) & (l <= tail(GWPCR$lambda, 1))
      p.v <- p.e.v & p.l.v

      # Interpolate density at requested efficiencies and lambda values, but
      # only were those values lie within the data matrix's range (i.e., no
      # extrapolation!)
      d <- rep(as.numeric(NA), length(e))
      d[p.v] <- density.interpolate(x0=e[p.v], y0=l[p.v], m=molecules)
      # And set density to zero if lambda is outside the data matrix's range.
      d[p.e.v & !p.l.v] <- 0
      # Result
      d
    } else
      NA
  })
}

dgwpcr.fun <- function(lambda, efficiency, molecules=1) {
  if (!is.numeric(efficiency) || (length(efficiency) != 1) || (efficiency < 0) || (efficiency > 1))
    stop('efficiency must be a numeric scalar within [0, 1]')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  # Create piece-wise monotone cubic spline interpolant
  d <- dgwpcr(GWPCR$lambda, efficiency=efficiency, molecules=molecules)
  splinefun(c(sum(head(GWPCR$lambda, 2) * c(3, -2)),
              sum(head(GWPCR$lambda, 2) * c(2, -1)),
              GWPCR$lambda,
              sum(tail(GWPCR$lambda, 2) * c(-1, 2)),
              sum(tail(GWPCR$lambda, 2) * c(-2, 3))),
            c(0, 0, d, 0, 0), method='monoH.FC')
}

pgwpcr.fun <- function(efficiency, molecules=1) {
  if (!is.numeric(efficiency) || (length(efficiency) != 1) || (efficiency < 0) || (efficiency > 1))
    stop('efficiency must be a numeric scalar within [0, 1]')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  # Handle corner-cases
  if (efficiency %in% c(0.0, 1.0))
    return(function(lambda) { ifelse(lambda < 1.0, 0, 1) })

  # Compute density at desired efficiency midpoints
  d <- dgwpcr(lambda=GWPCR$lambda, efficiency=efficiency, molecules=molecules)
  # Compute partial riemann sums up to interval endpoint
  p.midpoints <- cumsum(d * GWPCR$lambda.weights)
  p.midpoints <- p.midpoints / tail(p.midpoints, 1)

  # Use monotone interpolation to evaluate CDF at requested points.
  if (any(is.na(p.midpoints)))
    return(rep(as.numeric(NA), length(lambda)))
  splinefun(c(2*head(GWPCR$lambda, 1) - GWPCR$lambda.midpoints[1],
              head(GWPCR$lambda, 1),
              GWPCR$lambda.midpoints,
              tail(GWPCR$lambda, 1),
              2*tail(GWPCR$lambda.midpoints, 1) - tail(GWPCR$lambda, 1)),
            c(0, 0, p.midpoints, 1),
            method='monoH.FC')
}

pgwpcr <- function(lambda, efficiency, molecules=1){
  # Process each combination of efficiency and molecule count separately
  handle.parameters(list(lambda=lambda, efficiency=efficiency, molecules=molecules), by=c("efficiency", "molecules"), {
    tryCatch({ pgwpcr.fun(efficiency=efficiency, molecules=molecules)(lambda) },
             error=function(e) { rep(as.numeric(NA), length(lambda)) })
  })
}

gwpcr.sd <- function(efficiency, molecules=1) {
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  if (!is.list(GWPCR$sd.fun))
    GWPCR$sd.fun <- list()

  if ((molecules > length(GWPCR$sd.fun)) || is.null(GWPCR$sd.fun[[molecules]])) {
    # Ensure that we have a data matrix for the requested molecule count
    if ((molecules > length(GWPCR$data)) || is.null(GWPCR$data[[molecules]]))
      gwpcr.molecules.precompute(molecules=molecules)

    # Compute std. dev. of PCR distribution for pre-computed efficiencies
    s <- unlist(lapply(1:length(GWPCR$efficiency), FUN=function(i) {
      # Compute std.dev. of lambda for efficiency with index i
      sqrt(sum((GWPCR$lambda - 1.0)^2 * GWPCR$data[[molecules]][i,] * GWPCR$lambda.weights))
    }))

    # Estimate y-intercept of curve, i.e. maximal sd, assuming that
    # the std.dev. is reasonably linear within [0, min(efficiency)]
    s0 <- s[1] - GWPCR$efficiency[1] * (s[2]-s[1]) / (GWPCR$efficiency[2]-GWPCR$efficiency[1])

    # Create interpolating (monotone) spline.
    GWPCR$sd.fun[[molecules]] <-
      splinefun(c(0, GWPCR$efficiency, 1), c(s0, s, 0), method='monoH.FC')
  }

  # Evaluate spline at requested points
  r <- rep(as.numeric(NA), length(efficiency))
  e.v <- (efficiency >= 0) & (efficiency <= 1)
  r[e.v] <- GWPCR$sd.fun[[molecules]](efficiency[e.v])
  r
}

gwpcr.sd.inv <- function(sd, molecules=1) {
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  if (!is.list(GWPCR$sd.inv.fun))
    GWPCR$sd.inv.fun <- list()

  if ((molecules > length(GWPCR$sd.inv.fun)) || is.null(GWPCR$sd.inv.fun[[molecules]])) {
    y <- c(0, GWPCR$efficiency, 1)
    x <- gwpcr.sd(y, molecules=molecules)
    GWPCR$sd.inv.fun[[molecules]] <-
      splinefun(c(sum(head(x,2) * c(2, -1)), x, sum(tail(x,2) * c(-1, 2))),
                c(y[1],                      y, tail(y,1)                ),
                method='monoH.FC')
  }

  # Evaluate spline at requested points
  # Evaluate spline at requested points
  r <- rep(as.numeric(NA), length(sd))
  sd.v <- (sd >= 0)
  r[sd.v] <- GWPCR$sd.inv.fun[[molecules]](sd[sd.v])
  r
}

gwpcr.mixture <- function(x, FUN, efficiency, lambda0, molecules=1) {
  if (!is.function(FUN))
    stop("FUN must be a function with signature FUN(x, lambda)")

  # Process each combination of efficiency, lambda0 and molecule count separately
  handle.parameters(list(x=x, efficiency=efficiency, lambda0=lambda0, molecules=molecules), by=c("efficiency", "lambda0", "molecules"), {
    # Get probabilities for lambda lying within the intervals defined by
    # GWPCR$lambda.midpoints by multiplying the densities at the points lambda
    # by the interval lengths. Then set the probability to zero whereever it
    # is smaller than 1e-6 / <number of points> -- note that the the total
    # probability of intervals ignored that way is <= 1e-6.
    w <- dgwpcr(lambda=GWPCR$lambda, efficiency=efficiency, molecules=molecules) * GWPCR$lambda.weights
    w[w <= (1e-6 / length(w))] <- 0.0
    w <- w / sum(w)

    # For each lambda, evaluate the function at X. Note that each function call
    # can yield a vector -- we turn those into row vectors, and concatenate them
    # vertically to produce matrix d whose rows correspond to different lambdas.
    l <- GWPCR$lambda
    d <- do.call(rbind, lapply(1:length(l), function(i) {
      if (!is.na(w[i]) && (w[i] > 0.0)) FUN(as.vector(x), lambda=l[i]*lambda0) else rep(0, length(x))
    }))

    # Average the function's results over the range of lambda values, weighting
    # the value with the probability of the corresponding lambda interval from
    # above (vector w).
    apply(d, MARGIN=2, FUN=function(c) { sum(c * w) })
  })
}
