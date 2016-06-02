#!/usr/bin/env Rscript
library(parallel)

CORES <- 40

# GWPCR precomputed data
GWPCR <- new.env()
GWPCR$samples <- 1e8
GWPCR$until.molecules <- 1e6
GWPCR$lambda.def <- list(list(to=2, by=0.01),
                         list(to=10, by=0.1),
                         list(to=50, by=1))
GWPCR$efficiency.def <- list(list(to=0.05, by=0.05),
                             list(to=0.99, by=0.01))

# Generate efficiency grid from grid definition
GWPCR$efficiency <- unlist(local({
  p <- list(to=0);
  lapply(GWPCR$efficiency.def, function(e) {
    r <- tail(seq(from=p$to, to=e$to, by=e$by), -1)
    p <<- e
    r
  })
}))

# Generate lambda grid from grid definition
GWPCR$lambda <- c(0, unlist(local({
  p <- list(to=0);
  lapply(GWPCR$lambda.def, function(e) {
    r <- tail(seq(from=p$to, to=e$to, by=e$by), -1)
    p <<- e
    r
  })
})))

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

simulate <- function(efficiency, cycles, samples=1) {
  r <- mclapply(1:CORES, FUN=function(c) {
    # Produce one CORE-th of the total number of required samples
    samples.c <- ceiling(samples / CORES)
    # Start with one copy
    s.c <- rep(1, samples.c)
    for(i in 1:cycles)
      # Simulate a cycle
      s.c <- s.c + sapply(s.c, FUN=function(z) { rbinom(n=1, size=z, prob=efficiency) })
    s.c
  }, mc.preschedule=TRUE, mc.cores=CORES, mc.silent=FALSE)

  # Combine results from all cores, throw away unnecessary samples
  s <- unlist(r)[1:samples]

  # Scale to have expectation 1
  s / ((1+efficiency)**cycles)
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
  stopifnot(s >= 0)
  stopifnot(s <= L.MAX)
  message('Efficiency=', efficiency, ': Estimating bandwidth')
  GWPCR$bandwidth[e] <- bw.nrd0(s)
  message('Efficiency=', efficiency, ': Estimating density using bandwidth ', GWPCR$bandwidth[e])
  l <- density(s, kernel='gaussian', adjust=2, bw=GWPCR$bandwidth[e], from=0, to=L.MAX, n=L.MAX/L.INC)

  # Find indices into uniformly spaced density output that correspond to our non-uniform lambda grid
  i <- c(1, unlist(local({
    p <- list(to=0); m <- 1;
    lapply(GWPCR$lambda.def, function(e) {
      r <- tail(seq(from=m, to=m+(e$to-p$to)/L.INC, by=e$by/L.INC), -1)
      p <<- e
      m <<- max(r)
      r
    })
  })))
  stopifnot(abs(l$x[i] - GWPCR$lambda) <= 1e-12)
  d <- c(0, l$y[tail(i,-1)])

  # Normalize so that riemann sum is exactly one
  d <- d / sum(d * GWPCR$lambda.weights)

  # Store density
  GWPCR$data[[1]][e,] <- d

  # Save file
  message('Efficiency=', efficiency, ': Saving')
  save(GWPCR, file="sysdata.rda", compress='bzip2', compression_level=9)
}
