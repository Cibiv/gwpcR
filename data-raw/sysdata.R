#!/usr/bin/env Rscript
library(parallel)
library(Rcpp)
library(inline)

CORES <- 48

# GWPCR precomputed data
GWPCR <- new.env()
GWPCR$samples <- 1e9
GWPCR$until.molecules <- 1e7
GWPCR$density.threshold <- 1e-6
GWPCR$bandwidth.min <- 0.0025
GWPCR$bandwidth.max <- 0.01
GWPCR$lambda.def <- list(list(to=0.025,by=0.0025),
                         list(to=0.05, by=0.005),
                         list(to=0.95, by=0.01),
                         list(to=1.0,  by=0.005),
                         list(to=1.05, by=0.0025),
                         list(to=1.2,  by=0.005),
                         list(to=1.5,  by=0.01),
                         list(to=2,    by=0.02),
                         list(to=10,   by=0.1),
                         list(to=15,   by=0.5),
                         list(to=50,   by=5))
GWPCR$efficiency.def <- list(list(to=0.90, by=0.01),
                             list(to=0.94, by=0.005),
                             list(to=0.99, by=0.002))


make.grid <- function(def) {
  c(0, unlist(local({
    p <- list(to=0);
    lapply(def, function(e) {
      r <- tail(seq(from=p$to, to=e$to, by=e$by), -1)
      p <<- e
      r
    })
  })))
}

# Generate efficiency grid from grid definition
GWPCR$efficiency <- tail(make.grid(GWPCR$efficiency.def), -1)
message('Evaluating density for ', length(GWPCR$efficiency), ' efficiencies')

GWPCR$lambda <- make.grid(GWPCR$lambda.def)
message('Evaluating density for ', length(GWPCR$lambda), ' values of lambda')

# Generate lambda integration weights. We use the length an interval around
# each knot, places such that the break points between intervals lie exactly
# half-way between knots. Things look like this:
#     Knot i-1     Knot i         Knot i+1
#          .    |     .        |      .
#               |<- Weight i ->|
GWPCR$lambda.weights <- c((GWPCR$lambda[2] - GWPCR$lambda[1])/2,
                          (tail(GWPCR$lambda, -2) - head(GWPCR$lambda, -2))/2,
                          (GWPCR$lambda[length(GWPCR$lambda)] - GWPCR$lambda[length(GWPCR$lambda)-1])/2)

# Generate lambda knot midpoints.
GWPCR$lambda.midpoints <- (head(GWPCR$lambda, -1) + tail(GWPCR$lambda, -1))/2

# Density matrix
GWPCR$data <- list(matrix(nrow=length(GWPCR$efficiency),
                          ncol=length(GWPCR$lambda),
                          data=as.numeric(NA)))

# Density estimation bandwidth
GWPCR$bandwidth <- rep(as.numeric(NA), length(GWPCR$efficiency))

# Density estimation on non-uniform grid with non-uniform bandwidths
density.noniform <- function(r, efficiency) {
  # Since stats::density can only do bandwith estimation with a uniform
  # bandwidth, we use transform t before bandwidth estimation to decrease
  # the effective bandwidth within [0, 0.2] at low efficiencies (coincidentally,
  # also <= 0.2). This accounts for the large derivative of the density close
  # to zero for such low efficiencies (the density is approximately gamma!).
  # Transform t was chosen to have derivative 4 at 0, and approximately 1
  # on [0.2, inf].

  # Determine transformation parameters. Alpha determines where the
  # derivative of t reaches (almost) 1. beta is the derivative at zero
  # minus 1. If beta is zero, the transform is the identity.
  alpha <- 3/0.2
  beta <- 3 * min(1, max(1 - (efficiency - 0.1)/0.2, 0))
  t  <- function(x) { x - beta * exp(-x*alpha) / alpha }
  tp <- function(x) { 1 + beta * exp(-x*alpha) }

  # Transform samples r, and transform points x at which density is to be estimated.
  rp <- t(r)
  lp <- t(GWPCR$lambda)

  # Estimate density of transformed data. Bandwidth is estimates using
  # the original data, and clipped to [bandwidth.min, bandwith.max].
  # 10 estimates per bandwidth (i.e. kernel std.dev.) are produced.
  lp.mm <- range(lp)
  message('Efficiency=', efficiency, ': Estimating bandwidth')
  bw.raw <- bw.nrd0(r)
  bw <- min(max(GWPCR$bandwidth.min, bw.raw), GWPCR$bandwidth.max)

  message('Efficiency=', efficiency, ': Estimating density using bandwidth ', bw, ' (bw.nrd0 estimate was ', bw.raw, ')')
  d <- stats::density(rp, bw=bw, from=lp.mm[1], to=lp.mm[2],
                      n=ceiling((lp.mm[2]-lp.mm[1]) / (0.1*bw)))

  message('Efficiency=', efficiency, ': Computing density estimates on lambda grid')
  # Use monotone interpolation to find density estimates at requested points
  f <- splinefun(d$x, c(0, tail(d$y, -1)), method="monoH.FC")
  return(list(d=f(lp) * tp(GWPCR$lambda), bw=bw))
}


simulate.cycle <- cxxfunction(signature(samples_="numeric", efficiency_="numeric"), body='
  using namespace Rcpp;
  NumericVector result;
  RNGScope rngScope;

  // Arguments, and create result vector
  const NumericVector samples = as<NumericVector>(samples_);
  const double efficiency = as<double>(efficiency_);
  const std::size_t n = samples.size();
  result = NumericVector(n);

  // Simulate cycle
  for(std::size_t i=0; i < n; ++i)
    result[i] = samples[i] + R::rbinom(samples[i], efficiency);

  // Return result
  return(result);
', plugin='Rcpp', convention='.Call', language='C++')

simulate <- function(efficiency, cycles, samples=1) {
  # Produce one CORE-th of the total number of required samples
  n.c <- ceiling(samples / CORES)

  r <- mclapply(1:CORES, FUN=function(c) {
    # Start with one copy
    s.c <- rep(1, n.c)
    for(i in 1:cycles)
      # Simulate a cycle
      s.c <- simulate.cycle(s.c, efficiency)
    s.c
  }, mc.preschedule=TRUE, mc.cores=CORES, mc.silent=FALSE)

  # Combine results from all cores, throw away unnecessary samples, and scale
  # to have expectation 1
  f <- ((1+efficiency)**cycles)
  s <- rep(as.numeric(NA), samples)
  for(i in 1:length(r)) {
    l <- n.c * (i - 1) + 1
    u <- min(l + n.c, samples)
    if (l <= u)
      s[l:u] <- r[[i]][1:(u-l+1)] / f
  }

  return(s)
}

# Setup
RNGkind("L'Ecuyer-CMRG")

# Compute PCR distribution
for(e in 1:length(GWPCR$efficiency)) {
  efficiency <- GWPCR$efficiency[e]
  cycles <- ceiling(log(GWPCR$until.molecules)/log(1+efficiency))

  # Generate samples
  message('Efficiency=', efficiency, ': Generating ', GWPCR$samples, ' samples using ', cycles, ' cycles on ', CORES, ' cores')
  s <- simulate(efficiency, cycles, GWPCR$samples)

  # Estimate density. density() only supported outputting the esitmated density on a uniform grid!
  L.MAX <- max(GWPCR$lambda)
  L.INC <- min(diff(GWPCR$lambda))
  if (any(s < 0) || any(s > L.MAX) || any(is.na(s))) {
    cat('Indices containing negative samples:\n')
    print(which(s < 0))
    print(s[s < 0])
    cat('Indices containing samples beyond ', L.MAX, ':\n')
    print(which(s > L.MAX))
    print(s[!is.na(s) & (s > L.MAX)])
    cat('Indices containing infinite samples:\n')
    print(which(s == Inf))
    print(s[!is.na(s) & (s == Inf)])
    cat('Indices containing NA samples:\n')
    print(which(is.na(s)))
    print(s[is.na(s)])
    stop('invalid samples generated')
  }

  res <- density.noniform(s, efficiency=efficiency)
  GWPCR$bandwidth[e] <- res$bw
  d <- res$d

  # Pick a set of intervals (start with the ones where the density is smallest),
  # which all combined have probability less than one in a million, and set
  # the probability to zero there.
  p <- d * GWPCR$lambda.weights
  pc <- c(0, cumsum(sort(p)))
  z <- (p <= pc[which.max(pc >= GWPCR$density.threshold)-1])
  p[z] <- 0
  d[z] <- 0

  # Normalize so that riemann sum is exactly one
  d <- d / sum(p)

  # Store density
  GWPCR$data[[1]][e,] <- d

  # Save file
  message('Efficiency=', efficiency, ': Saving')
  save(GWPCR, file="R/sysdata.new.rda", compress='bzip2', compression_level=9)
}
