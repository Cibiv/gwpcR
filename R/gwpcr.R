# **********************************************************************
#    PCR Distribution as produced by a simple Galton-Watson process
# **********************************************************************

handle.parameters <- function(parameters, by, expr) {
  n <- 0
  for(p in names(parameters)) {
    if (!is.numeric(parameters[[p]]) && (length(parameters[[p]]) == 0))
      stop(p, ' must be a non-empty numeric vector')
    n <- max(n, length(parameters[[p]]))
  }
  for(p in names(parameters)) {
    if (n %% length(parameters[[p]]) != 0)
      warning("longest parameter length ", n, " s not a multiple of length of ", p)
  }

  # Create data.table containing the parameter values as columns
  t <- suppressWarnings(do.call(data.table::data.table,
                                c(list(`__result__` = as.numeric(NA)),
                                  parameters)))

  # Add result column by evaluting the provided expression within each by group
  # It's a bit tricky to get the expression to be evaluated with the correct
  # stack of nested frames so that both parameters (i.e. columns of t) *and*
  # arbitrary stuff defined in our caller's frame are accessible. First, we
  # define a function r(...), which evaluates the given expression in an
  # environment which contains its parameters, and as the enclosing env our
  # parent's frame.
  r <- function(...) {
    eval(expr, envir=list(...), enclos=r.enclos)
  }
  # Then we make sure the function can access those things regardless of how
  # its called, and that is contains the actual expression and enclosing env.
  environment(r) <- list2env(list(expr=substitute(expr),
                                  r.enclos=parent.frame()),
                             parent=baseenv())
  # We now create an expression e which reads
  #   as.numeric(r(parameter1=parameter2, parameter2=parameter2, ...))
  # for all the parameters in our parameter list.
  p <- lapply(names(parameters), as.symbol)
  names(p) <- names(parameters)
  e <- bquote(as.numeric(.(e)), list(e=as.call(c(list(as.symbol('r')), p))))
  # Finally, we evaluate the expression for each group, and store the result
  t[, `__result__` := eval(e), by=by]

  # Return result
  return(t$`__result__`)
}

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

dgwpcrpois <- function(c, efficiency, lambda0, threshold=0, molecules=1) {
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

pgwpcrpois <- function(c, efficiency, lambda0, threshold=0, molecules=1) {
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
                      d[c.v] <- gwpcr.mixture(c[c.v], stats::ppois, efficiency=efficiency,
                                              lambda0=lambda0, molecules=molecules)
                      # And compute P(X = c | X >= Th)
                      d / (1.0 - p)
                    })
}
