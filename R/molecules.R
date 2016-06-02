# ******************************************************************************************
#    Extend Precomputed Galton-Watson PCR Distribution to Different Nr. of initial Molecules
# ******************************************************************************************

gwpcr.molecules.precompute <- function(molecules) {
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  # To compute a data matrix for M molecules, we need to convolve each
  # row of the data matrx for 1 molecules M times with itself. This can
  # be done efficiently by fourier transforming, taking the M-th power,
  # and reverse transforming. Since we use the discrete transform, we need
  # to be carefully, however, to first pad the vector to a sufficient length
  # to be able to hold the resulting distribution, which (as a sum of M
  # i.i.d. random variables) has a M-fold larger domain than the original.
  # The grid on which the original density is defined also needs to be uniform,
  # which the raw data row violate (GWPCR$lambda is non-uniform!). After
  # computing the convolution, we transform back to the original non-uniform
  # grid.
  data <- matrix(NA, nrow=nrow(GWPCR$data[[1]]), ncol=ncol(GWPCR$data[[1]]))
  # Determine smallest grid size used within rows of the data matrix, and
  # the largest x-value for which data is provided. Then synthesize uniform
  # grid l, which is what will be used for the convolution. Note that is
  # extend up to M times the original largest x-value!
  w <- min(diff(GWPCR$lambda))
  l.max <- max(GWPCR$lambda)
  l <- seq(from=0, to=l.max*molecules, by=w)
  lp <- l / molecules

  # Process each row (i.e. efficiency) separately.
  for(e in 1:length(GWPCR$efficiency)) {
    # Transform original density to grid l by spline interpolation. Be carefull
    # to interpolate only within the original data range, and zero-extend beyond
    # that.
    d <- c(pmax(splinefun(GWPCR$lambda, GWPCR$data[[1]][e,])(l[l <= l.max]), 0),
           rep(0, sum(l > l.max)))
    # Compute the convolution, and transform back to the original grid, again
    # by spline interpolation
    dp <- pmax(splinefun(lp, Re(fft(fft(d)**molecules, inverse = TRUE)))(GWPCR$lambda), 0)
    # Rescale the transformed density to have riemann sum 1
    dp <- dp / sum(dp * GWPCR$lambda.weights)
    # Store into output data matrix
    data[e,] <- dp
  }

  GWPCR$data[[molecules]] <- data
}
