#' Binomial Galton-Watson Model for the Polymerase Chain Reaction (PCR)
#'
#' @description Limit distribution of amount of PCR product relative to the
#'   expected value after many cycles as prediced by a binomial Galton-Watson
#'   model of the polymerase chain reaction (PCR), Parameters are
#'   \var{efficiency} and \var{molecules} -- the former is the probability that
#'   a molecule is duplicated during a particular PCR cycle, the latter the
#'   initial number of molecules in the reaction mix.
#'
#'   Note that a \emph{molecule} refers to single strand here, and complementary
#'   strands are not distinguished. Setting \var{molecules=2} thus models a PCR
#'   reaction starting from a single piece of double-stranded DNA.
#'
#'   \code{dgwpcr} (resp. \code{pgwpcr}) evaluate the density (resp. the CDF) at
#'   the given relative amount of PCR product \code{l} for parameters
#'   \var{efficiency} and \var{molecules}. \code{dgwpcr.fun} (resp.
#'   \code{pgwpcr.fun}) return a unary function which represents the density
#'   (resp. CDF) for the given parameters. \code{rgwpcr} draws random samples by
#'   simulation. If the number of PCR cycles to use is not specified, the
#'   simulation is stopped once the expected absolute amount of PCR product
#'   reaches one million molecules, at which point the distribution is
#'   considered to be close to the limit distribution for infinitly many cycles.
#'
#' @param n number of random samples to generate
#'
#' @param l molecular family size relative to the average
#'
#' @param efficiency efficiency of amplification
#'
#' @param molecules initial copy number
#'
#' @param cycles number of amplification cycles used for simulation. By default,
#'   a large enough value is used to make the results virtually idistinguishable
#'   from the limit for \eqn{cycles \to \infty}
#'
#' @details The binomial Galton-Watson PCR model treats PCR as a branching
#'   process. At time 0, the absolute number of molecules \eqn{c_n} is the
#'   initial copy number \var{molecules}. Each time step from \eqn{c_n} to
#'   \eqn{c_{n+1}}{c_(n+1)} corresponds to a PCR cycle and duplicates each of
#'   the \eqn{c_n} molecules with the probability specified in parameter
#'   \var{efficiency} (\var{E}). Thus,
#'
#'   \deqn{c_{n+1} = c_n + \textrm{Binomial}(c_n, E).}{c_(c+1) = c_n +
#'   Binomial(c_n, E).}
#'
#'   Each cycle thus increases the expected molecule count by a factor of
#'   \eqn{(1+E)}. The \emph{relative} amount of PCR product after \eqn{n} cycles
#'   is therefore
#'
#'   \deqn{l_n = c_n \cdot (1+E)^{-n}.}{l_n = c_n (1+E)^-n.}
#'
#'   \code{dgwpcr} (resp. \code{rgwpcr}) is the density (resp. CDF) of the a.c.
#'   limit in distribution of \eqn{l_n} for \eqn{n \to \infty}{n to \infty}.
#'   Essentially, the reason that \eqn{l_n} converges in distribution is that
#'   the larger the number of molecules \eqn{c_n}, the smaller the additional
#'   variability introduced into \eqn{c_{n+1}}{c_(n+1)} by the term
#'   \eqn{\textrm{Binomial}(c_n, E).}{Binomial(c_n, E).}. For reasonably large
#'   efficiencies, that becomes quickly negligible compared to \eqn{c_n}. Assume
#'   e.g. an efficiency of 50\%, and that \eqn{c_n=10^4}{c_n=2500}. The standard
#'   deviation of the binomial term is then already only 25, i.e. 1\% of
#'   \eqn{c_n}, and this state can be expected to be reached after about 20
#'   cycles.
#'
#'   For this reason, using the limit distribution still yields a model of the
#'   polymerase chain reaction that is sufficiently accurate for most purposes.
#'
#' @seealso \code{\link{gwpcr.sd}}
#' @seealso \code{\link{gwpcr.mixture}}
#' @seealso \code{\link{gwpcrpois}}
#'
#' @name gwpcr
NULL

#' @rdname gwpcr
#' @useDynLib gwpcR gwpcr_simulate_c
#' @export
rgwpcr <- function(n, efficiency, molecules=1, cycles=NA) {
  if (!is.numeric(efficiency) || (length(efficiency) != 1) || (efficiency <= 0) || (efficiency > 1))
    stop('efficiency must be a numeric scalar within (0,1]')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  # Determine how many cycles are necessary on average to produce 1e6 molecules.
  # After that point, we assume that the additional variability is negligible.
  if (is.na(cycles))
    cycles <- ceiling(log(1e6 / molecules)/log(1+efficiency))

  # Start with initial copy number
  if (efficiency < 1.0) {
    # Flip a coin for each molecule in each sample to determine whether its
    # copied. That can be done efficiently by sampling from a binomial
    # distribution for each sample.gwpcrpois_simulate_c Note: Since this
    # requires looping over all samples, the code was translated to C. See
    # gwpcr_simulate in simulate.c.
    if (TRUE)
      .C(gwpcr_simulate_c,
         nsamples=as.integer(n),
         samples=double(n),
         efficiency=as.double(efficiency),
         molecules=as.double(molecules),
         cycles=as.integer(cycles),
         NAOK=TRUE)$samples
    else {
      s <- rep(1.0, n)
      for(i in 1:cycles) {
        sp <- sapply(s, FUN=function(z) { rbinom(n=1, size=z, prob=efficiency) })
        s <- s + sp
      }
      s / (molecules * (1+efficiency)**cycles)
    }
  }
  else
    rep(1.0, n)
}

#' @rdname gwpcr
#' @export
dgwpcr <- function(l, efficiency, molecules=1) {
  # Process each number of molecules separately
  handle.parameters(list(l=l, efficiency=efficiency, molecules=molecules), by="molecules", {
    # Moleculec ount must be >= 1 and integral
    if (is.finite(molecules) && (molecules == floor(molecules)) && (molecules > 0)) {
      # Determine molecule count filter (p.m), efficiency and lambda vectors
      # (e, l), their validity filter (p.e.i, p.e.g, p.l.v) and overall filters
      # (p.i, p.g). *.i and *.g stand for "interpolate" and "gamma", and select
      # whether the density is interpolated, or approximated using a gamma dist.
      e <- efficiency
      p.e.i <- (e > E.MIN) & (e <= E.MAX)
      p.e.g <- (e >= 0) & (e < E.GAMMA.TH)
      p.l.v <- (l >= 0) & (l <= tail(GWPCR$lambda, 1))
      p.i <- p.e.i & p.l.v
      p.g <- p.e.g & p.l.v

      # Prepare result vector.
      d <- rep(as.numeric(NA), length(e))
      d[p.i | p.g] <- 0

      # In the cut-over region between precomputed and gamma density, determine
      # weight factor f for gamma density.
      f <- gamma.factor(e)
      stopifnot(all(f[p.e.i] < 1))
      stopifnot(all(f[p.e.g] > 0))
      stopifnot(all((0 <= f) & (f <= 1)))

      # Interpolate density at requested efficiencies and lambda values, but
      # only were those values lie within the data matrix's range (i.e., no
      # extrapolation!)
      if (sum(p.i) > 0)
        d[p.i] <- d[p.i] + (1 - f[p.i]) * density.interpolate(x0=e[p.i], y0=l[p.i], m=molecules)

      # Use gamma for efficiencies below the smallest precomputed efficiency
      # The gamma parameters are chosen to produce a distribution with mean 1 and
      #               1 - E     1
      #   Var( L ) = ------- * ---,
      #               1 + E     m
      # which matches the variance of the true distribution (see gwpcr.sd).
      g.par <- molecules * (1+e[p.g]) / (1-e[p.g])
      if (length(g.par) > 0)
        d[p.g] <- d[p.g] + f[p.g] * dgamma(l[p.g], shape=g.par, rate=g.par)

      # And set density to zero if lambda is outside the data matrix's range.
      d[(p.e.i | p.e.g) & !p.l.v] <- 0

      # Result
      d
    } else
      NA
  })
}

#' @rdname gwpcr
#' @export
dgwpcr.fun <- function(efficiency, molecules=1) {
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

#' @rdname gwpcr
#' @export
pgwpcr.fun <- function(efficiency, molecules=1) {
  if (!is.numeric(efficiency) || (length(efficiency) != 1) || (efficiency < 0) || (efficiency > 1))
    stop('efficiency must be a numeric scalar within [0, 1]')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  # "Round" efficiencies either down to the maximal simulated efficiency, or up to 1,
  # depending on which is closer.
  efficiency <- if ((efficiency < 0) || (efficiency > 1.0))
    NA
  else if (efficiency <= 0.5 + 0.5*E.MAX)
    pmin(efficiency, E.MAX)
  else
    1.0

  # In the cut-over region between precomputed and gamma density, determine
  # weight factor f for gamma density.
  f <- gamma.factor(efficiency)

  # CDF computed from pre-computed densities by numeric integration
  cdf.i <- if ((efficiency > E.MIN) && (efficiency <= E.MAX)) {
    # Determine integration knots and weights for integration on [0, 1]
    gq <- statmod::gauss.quad(2, kind = "legendre")
    q <- (1 + gq$nodes) / 2
    w <- gq$weights / 2

    # We integrate over each of the intervals between points, using
    # (scaled versions of) the knots from above. The densities at
    # the knots are determined by interpolation (see density.interpolate),
    # and the value of the definitive integral at (original) point i
    # is then determined by summing over the definitive integrals of
    # the precending intervals. This yield the value of Int(0, l) for
    # the points l from GWPCR$lambda.
    x0 <- head(GWPCR$lambda, -1)
    dx <- diff(GWPCR$lambda)
    p <- matrix(NA, nrow=2, ncol=length(x0))
    p[1,] <- x0 + dx * q[1]
    p[2,] <- x0 + dx * q[2]
    v <- density.interpolate(rep(efficiency, length(p)), as.vector(p), molecules)
    dim(v) <- c(2, length(v)/2)
    y <- cumsum(c(0, (v[1,] * w[1] + v[2,] * w[2]) * dx))
    y <- y / tail(y, 1)

    # Interpolate between the CDF at points monotonically.
    splinefun(GWPCR$lambda, y, method="monoH.FC")
  } else {
    # Determine dummy function if outside of precomputed efficiency range
    function(l) { rep(0, length(l)) }
  }

  # CDF computed from the gamma approximation
  cdf.g <- if (efficiency < E.GAMMA.TH) {
    # Gamma approximation of the true CDF, starts being reasonable below E,GAMMA.TH
    g.par <- molecules * (1+efficiency) / (1-efficiency)
    function(l) {
      pgamma(l, shape=g.par, rate=g.par)
    }
  } else {
    # Determine dummy function if outside of gamma-approximated range
    function(l) { rep(0, length(l)) }
  }

  # Return a function approximating the CDF.
  function(lambda) {
    # Create result value and determine valid lambda values
    r <- rep(NA, length(lambda))
    p.v <- (lambda >= 0) & (lambda <= tail(GWPCR$lambda, 1))

    # Compute results, use 0 resp. 1 outside the valid lambda range.
    r[lambda < 0] <- 0
    r[p.v] <- f * cdf.g(lambda[p.v]) + (1 - f) * cdf.i(lambda[p.v])
    r[lambda > tail(GWPCR$lambda, 1)] <- 1

    # Return result
    r
  }
}

#' @rdname gwpcr
#' @export
pgwpcr <- function(l, efficiency, molecules=1){
  # Process each combination of efficiency and molecule count separately
  handle.parameters(list(l=l, efficiency=efficiency, molecules=molecules), by=c("efficiency", "molecules"), {
    tryCatch({ pgwpcr.fun(efficiency=efficiency, molecules=molecules)(l) },
             error=function(e) { rep(as.numeric(NA), length(l)) })
  })
}

#' PCR Product Distribution Standard Deviation
#'
#' @description For efficiency \var{E} and initial number of molecules \var{m},
#'   the PCR product distribution has \emph{variance} (not standard deviation!)
#'   \eqn{\frac{1-E}{1+E}\cdot\frac{1}{m}}{(1+E) / ((1-E) * m)}.
#'
#'   Function \code{gwpcr.sd} uses this to compute the standard deviation (i.e.
#'   the square of the above) for a given efficiency and initial copy number,
#'   and \code{gwpcr.sd.inv} does the reverse and computes the efficiency given
#'   standard deviation and copy number.
#'
#' @inheritParams gwpcr
#'
#' @param sd Standard deviation of the PCR product distribution
#'
#' @seealso \code{\link{gwpcr}}
#'
#' @export
gwpcr.sd <- function(efficiency, molecules=1) {
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  # The formula
  #               1 - E     1
  #   Var( L ) = ------- * ---
  #               1 + E     m
  # follows from the theory of general Galton-Watson processes by re-scaling
  # with the expectation to find the convergent process. For m initial molecules,
  # the variance is the variance of the m-fold average of the single molecule
  # case, i.e. the std. dev. is one sqrt(m)-th of that single-molecule std. dev.
  r <- rep(as.numeric(NA), length(efficiency))
  e.v <- (efficiency >= 0) & (efficiency <= 1)
  r[e.v] <- sqrt( (1-efficiency[e.v]) / (molecules * (1 + efficiency[e.v])) )
  r
}

#' @rdname gwpcr.sd
#' @export
gwpcr.sd.inv <- function(sd, molecules=1) {
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  # The formula
  #        1 - m * Var(L)
  #   E = ----------------
  #        1 + m * Var(L)
  # is found by solving for E in the formula in gwpcr.sd()
  r <- rep(as.numeric(NA), length(sd))
  e.v <- (sd >= 0)
  v <- pmin(molecules * sd[e.v] ^ 2, 1)
  r[e.v] <- (1 - v) / (1 + v)
  r
}

#' Mixtures of Distributions with PCR-distributed Weights
#'
#' Computes mixtures (i.e. convex linear combinations) of arbitrary functions
#' with PCR-distributed weights.
#'
#' @details Numerically approximates the integral
#'
#'   \deqn{\int_0^\infty F(x,\lambda) \cdot \textrm{dgwpcr}(\lambda)
#'   \,d\lambda}{Int F(x,\lambda) dgwpcr(\lambda) d\lambda over [0, Infinity)}
#'
#'   where \eqn{\textrm{dgwpcr}}{dgwpcr} is the density function of the PCR
#'   product distribution for the specified efficiency and initial number of
#'   molecules.
#'
#'   Function \eqn{F} is usually the pdf (probability density function) or cdf
#'   (cumulative density function) of a probability distribution, in which case
#'   \code{gwpcr.mixture} computes the pdf (resp. cdf) of mixture of F's with
#'   PCR-distributed weights.
#'
#'   \code{grid.width.fun} can be used to control the coarseness of the
#'   integration grid. This should be a unary function that takes a value of
#'   \eqn{\lambda} and returns the maximum acceptable distance between grid
#'   points in the vicinity of \eqn{\lambda} value. If \eqn{F(.,\lambda)} is a
#'   pdf or cdf with \eqn{\lambda}-dependent variance, \code{grid.width.fun} can
#'   be used to enforce a finer grid in regions where the variance is small.
#'
#' @inheritParams gwpcr
#'
#' @param x value to evaluate the mixture at
#'
#' @param FUN function to mix
#'
#' @param grid.width.fun functions which returns the maximum grid size (i.e.
#'   distance between points) depending on \eqn{\lambda}. If the variance of
#'   \eqn{F} depends strongly on \eqn{\lambda}, this can be used to ensure that
#'   \eqn{F} is evaluated on finer grid for values of \eqn{lambda} where the
#'   variance is small.
#'
#' @seealso \code{\link{gwpcr}}
#'
#' @export
gwpcr.mixture <- function(x, FUN, efficiency, molecules=1, grid.width.fun = NULL) {
  if (is.null(grid.width.fun))
    grid.width.fun <- function(x) { Inf }
  if (!is.function(FUN))
    stop("FUN must be a function with signature FUN(x, lambda)")
  if (!is.numeric(efficiency) || (length(efficiency) != 1) || (efficiency < 0) || (efficiency > 1))
    stop('efficiency must be a numeric scalar within [0, 1]')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')
  if (!is.function(grid.width.fun))
    stop("grid.width.fun must be a function with signature grid.width.fun(lambda)")

  # "Round" efficiencies either down to the maximal simulated efficiency, or up to 1,
  # depending on which is closer.
  efficiency <- if ((efficiency < 0) || (efficiency > 1.0))
    NA
  else if (efficiency <= 0.5 + 0.5*E.MAX)
    pmin(efficiency, E.MAX)
  else
    1.0

  if (is.na(efficiency))
    rep(as.numeric(NA), length(x))
  else if (efficiency < 1.0) {
    # Refine the lambda grid to ensure that the grid width is at least
    # grid.scale times the standard deviation, at all values of lambda we
    # integrate over.
    r <- refine(GWPCR$lambda, grid.width.fun(GWPCR$lambda))
    l <- r$points
    l.w <- r$weights

    # Get probabilities for lambda lying within the intervals defined by
    # GWPCR$lambda.midpoints by multiplying the densities at the points lambda
    # by the interval lengths. Then set the probability to zero whereever it
    # is smaller than 1e-6 / <number of points> -- note that the the total
    # probability of intervals ignored that way is <= 1e-6.
    w <- dgwpcr(l, efficiency=efficiency, molecules=molecules) * l.w
    w[w <= (1e-6 / length(w))] <- 0.0
    w <- w / sum(w)

    # For each lambda, evaluate the function at X. Note that each function call
    # can yield a vector -- we turn those into row vectors, and concatenate them
    # vertically to produce matrix d whose rows correspond to different lambdas.
    d <- do.call(rbind, lapply(1:length(l), function(i) {
      if (!is.na(w[i]) && (w[i] > 0.0)) FUN(as.vector(x), l[i]) else rep(0, length(x))
    }))

    # Average the function's results over the range of lambda values, weighting
    # the value with the probability of the corresponding lambda interval from
    # above (vector w).
    apply(d, MARGIN=2, FUN=function(c) { sum(c * w) })
  } else {
    # For efficiency 1.0, the PCR distribution has a single mass at lambda=1
    FUN(as.vector(x), 1)
  }
}
