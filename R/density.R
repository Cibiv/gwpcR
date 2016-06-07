# ***************************************************************************************
# Density Interpolation, calles Fotran Code from Package "akima"
# ***************************************************************************************

density.interpolate <- function(x0, y0, m) {
  # Akima's Fortran code allows some pre-computation steps (finding derivatives)
  # to be done only once for fixed knots and their z-values, instead of
  # re-doing them for each list of points to interpolate at (i.e. each call
  # of the interpolation function). However, the function bicubic() in package
  # akima does not expose that functionality - it always passes '1' for
  # parameter md. We therefore call akima's procedure rgbi3p directly here,
  # instead of going through the akima package. We still required the package,
  # though, because it includes the necessary shared library containing the
  # compiled Fortran code.
  requireNamespace('akima', quietly=TRUE)

  # Ensure that we have a data matrix for the requested molecule count
  stopifnot(m >= 1)
  if ((m > length(GWPCR$data)) || is.null(GWPCR$data[[m]]))
    gwpcr.molecules.precompute(molecules=m)

  # Pre-process data by calling Akima's procedure with MD=1.
  if (is.null(GWPCR$data.akima.wk))
    GWPCR$data.akima.wk <- list()
  if ((m > length(GWPCR$data.akima.wk)) || is.null(GWPCR$data.akima.wk[[m]])) {
    # These checks cannot fail, unless the data in sysdata.rda is tampared with
    #stopifnot(is.double(GWPCR$efficiency) && is.double(GWPCR$lambda))
    #stopifnot(is.double(GWPCR$data[[m]]))
    #stopifnot(length(GWPCR$efficiency) == nrow(GWPCR$data[[m]]))
    #stopifnot(length(GWPCR$lambda) == ncol(GWPCR$data[[m]]))
    r <- .Fortran("rgbi3p",
                  md=as.integer(1),
                  nxd=length(GWPCR$efficiency), nyd=length(GWPCR$lambda),
                  xd=GWPCR$efficiency, yd=GWPCR$lambda, zd=GWPCR$data[[m]],
                  nip=as.integer(1),
                  xi=GWPCR$efficiency[1], yi=GWPCR$lambda[1], zi=double(1),
                  ier=integer(1),
                  wk=double(3*length(GWPCR$efficiency)*length(GWPCR$lambda)),
                  PACKAGE="akima", NAOK=TRUE)
    stopifnot(r$ier == 0)
    GWPCR$data.akima.wk[[m]] <- r$wk
  }

  # These checks cannot fail, unless the data in sysdata.rda is tampared with
  #stopifnot(is.double(GWPCR$efficiency) && is.double(GWPCR$lambda))
  #stopifnot(is.double(GWPCR$data[[m]]))
  #stopifnot(length(GWPCR$efficiency) == nrow(GWPCR$data[[m]]))
  #stopifnot(length(GWPCR$lambda) == ncol(GWPCR$data[[m]]))
  #stopifnot(is.double(GWPCR$data.akima.wk[[m]]))
  #stopifnot(length(GWPCR$data.akima.wk[[m]]) == 3*length(GWPCR$efficiency)*length(GWPCR$lambda))
  stopifnot(is.double(x0) && is.double(y0))
  stopifnot(length(x0) == length(y0))
  if (length(x0) == 0)
    return(double(0))
  r <- .Fortran("rgbi3p",
                md=as.integer(2),
                nxd=length(GWPCR$efficiency), nyd=length(GWPCR$lambda),
                xd=GWPCR$efficiency, yd=GWPCR$lambda, zd=GWPCR$data[[m]],
                nip=length(x0),
                xi=x0, yi=y0, zi=double(length(x0)),
                ier=integer(1),
                wk=GWPCR$data.akima.wk[[m]],
                PACKAGE="akima", NAOK=TRUE)
  stopifnot(r$ier == 0)
  return(pmax(r$zi, 0))
}


# ***************************************************************************************
# Extend Precomputed Galton-Watson PCR Distribution to Different Nr. of initial Molecules
# ***************************************************************************************

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

