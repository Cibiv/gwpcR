# density.R, Copyright 2016,2017 Florian G. Pflug
#
# This file is part of gwpcR
#
# gwpcR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gwpcR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with gwpcR.  If not, see <http://www.gnu.org/licenses/>.

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
    r <- .Fortran('rgbi3p',
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
  r <- .Fortran('rgbi3p',
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


#' Precompute PCR Product Distribution for a larger Number of initial Molecules
#'
#' @description This function is called automatically whenever a value of
#'   \code{molecules} larger than one is passed to any of the package's other
#'   functions for the first time, so this function usually isn't needed and
#'   shouldn't be called.
#'
#'   However, when using the \code{parallel} package to spread computions over
#'   multiple CPUs, it can be advantageous to do the necessary computations
#'   once, instead within each of the worker processes launched by the
#'   \code{parallel} package. In this case, call
#'   \code{gwpcr.molecules.precompute} before lauching any workers.
#'
#' @inheritParams gwpcrpois
#'
#' @export
gwpcr.molecules.precompute <- function(molecules) {
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')

  # Data for a single molecule is computed by sysdata.R and loaded from sysdata.rda
  # Thus, for molecules==1 we always exit immediately
  if ((molecules <= length(GWPCR$data)) && !is.null(GWPCR$data[[molecules]]))
    return(GWPCR$data[[molecules]])

  # The distribution for <molecules> initial molecules is computed by
  # convolving the densities for m1 and m_2 initial molecules, where m = m1+m2.
  m1 <- floor(molecules / 2)
  m2 <- ceiling(molecules / 2)
  if ((m1 > length(GWPCR$data)) || is.null(GWPCR$data[[m1]]))
    gwpcr.molecules.precompute(m1)
  if ((m2 > length(GWPCR$data)) || is.null(GWPCR$data[[m2]]))
    gwpcr.molecules.precompute(m2)
  message('gwpcr.molecules.precompute: Computing PCR distribution for ',
          molecules, ' initial molecules from those for ', m1, ' and ', m2)

  # Determine smallest grid size used within rows of the data matrix, and
  # the largest x-value for which data is provided. Then synthesize uniform
  # grid l, which is what will be used for the convolution. Note that since
  # it needs to hold the *unscaled* sum, it must extend to twice l.max!
  w <- min(diff(GWPCR$lambda))
  l.max <- max(GWPCR$lambda)
  l <- seq(from=0, to=2*l.max, by=w)

  # If m1 == m2, then the resulting distribution for m1 + m2 == molecules is
  # simply the convolution of the distribution for m1 (resp. m2) with itself,
  # and then re-scaled in x-direction to again have expectation 1. But if
  # m1 == m2 - 1, then what we actually compute is the weighted average
  #   X_m = (m1 / m) X_m1 + (m2 / m) * X_m2
  #       = (m2 / m) [ (m1 / m2) X_m1 + X_m2 ].
  # We thus need to re-scale X_m1 in x-direction to have mean m1 / m2 < 1,
  # then X_m2, and re-scale the sum from mean (m1 / m2 + 1) back to 1.
  l1 <- if (m1 != m2) GWPCR$lambda * (m1 / m2) else GWPCR$lambda
  lp <- if (m1 != m2) l * (m2 / molecules) else l / 2

  # Get data matrices for the both summands.
  data <- matrix(NA, nrow=nrow(GWPCR$data[[1]]), ncol=ncol(GWPCR$data[[1]]))
  data1 <- GWPCR$data[[m1]]
  data2 <- GWPCR$data[[m2]]

  # Process each row (i.e. efficiency) separately.
  for(e in 1:length(GWPCR$efficiency)) {
    # Apply x-asis scaling to distribution X_m1 for m1 molecules as explained
    # above, and evaluate on the uniform grid l. Beyond the original domain
    # of X_m1's density don't extrapolate with splinefun(), but instead set zero!
    f1 <- fft(c(pmax(splinefun(l1, data1[e,])(l[l <= l.max]), 0),
                rep(0, sum(l > l.max))))

    # Evaluated istribution X_m2 for m2 molecules on the uniform grid l. Beyond
    # the original domain of X_m1's density don't extrapolate with splinefun(),
    # but instead set zero!
    f2 <- if (m1 != m2)
      fft(c(pmax(splinefun(GWPCR$lambda, data2[e,])(l[l <= l.max]), 0),
            rep(0, sum(l > l.max))))
    else
      f1

    # Convolute to compute (weighted) sum of X_m1 and X_m2, and evaluate result
    # on the non-uniform grid used to store distributions (GWPCR$lambda).
    d <- pmax(splinefun(lp, Re(fft(f1 * f2, inverse = TRUE)))(GWPCR$lambda), 0)

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
    stopifnot(abs(sum(GWPCR$lambda * d * GWPCR$lambda.weights) - 1) < 1e-2)

    # Store into output data matrix
    data[e,] <- d
  }

  GWPCR$data[[molecules]] <- data
}
