# **********************************************************************
#    PCR Distribution as produced by a simple Galton-Watson process
# **********************************************************************

rgwpcr <- function(samples, efficiency, molecules=1) {
  if (!is.numeric(efficiency) || (length(efficiency) != 1))
    stop('efficiency must be a numeric scalar')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  # Determine how many cycles are necessary on average to produce 1e6
  # molecules. After that point, we assume that the additional variability
  # is negligible.
  cycles <- ceiling(log(1e6 / molecules)/log(1+efficiency))

  # Start with initial copy number
  s <- rep(molecules, samples)
  if (efficiency == 1.0)

  if (efficiency < 1.0) {
    # Flip a coin for each molecule in each sample to determine whether
    # its copied. That can be done efficiently by sampling from a binomial
    # distribution for each sample.
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
  if (!is.numeric(lambda))
    stop('lambda must be a numeric vector')
  if (!is.numeric(efficiency))
    stop('efficiency must be a numeric vector')
  if (!is.numeric(molecules) || any(floor(molecules) != molecules) || any(molecules < 1))
    stop('molecules must be a positive integral vector')

  # Bring efficiency, lambda and molecules to the same length.
  n <- max(length(efficiency), length(lambda), length(molecules))
  if ((n %% length(efficiency) != 0) || (n %% length(lambda) != 0) ||
      (n %% length(molecules) != 0))
    warning("longest object length is not a multiple of all shorter object lengths")
  efficiency <- rep(efficiency, length.out=n)
  lambda <- rep(lambda, length.out=n)
  molecules <- rep(molecules, length.out=n)

  # Group parameters by initial molecule count
  d <- rep(as.numeric(NA), n)
  for(m in unique(molecules)) {
    # Moleculec ount must be >= 1 and integral
    if ((m != as.integer(m)) || (m <= 0))
      next

    # Ensure that we have a data matrix for the requested molecule count
    if ((m > length(GWPCR$data)) || is.null(GWPCR$data[[m]]))
      gwpcr.molecules.precompute(molecules=m)

    # Determine molecule count filter (p.m), efficiency and lambda vectors
    # (e, l), their validity filter (p.e.v, p.l.v) and overall filter (p.v)
    p.m <- (molecules == m)
    e <- efficiency[p.m]
    l <- lambda[p.m]
    p.e.v <- (e >= head(GWPCR$efficiency, 1)) & (e <= tail(GWPCR$efficiency, 1))
    p.l.v <- (l >= 0) & (l <= tail(GWPCR$lambda, 1))
    p.v <- p.m & p.e.v & p.l.v

    # Interpolate density at requested efficiencies and lambda values, but
    # only were those values lie within the data matrix's range (i.e., no
    # extrapolation!)
    r <- akima::bicubic(x=GWPCR$efficiency, y=GWPCR$lambda, z=GWPCR$data[[m]],
                        x0=efficiency[p.v],
                        y0=lambda[p.v])

    # Store interpolated densities are correct output positions
    d[p.v] <- pmax(r$z, 0)
    # And set density to zero if lambda is outside the data matrix's range.
    d[p.m & p.e.v & !p.l.v] <- 0
  }

  # Return result
  return(d)
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
  if (!is.numeric(lambda))
    stop('lambda must be a numeric vector')

  pgwpcr.fun(efficiency=efficiency, molecules=molecules)(lambda)
}
