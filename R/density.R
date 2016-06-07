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
  if ((m > length(GWPCR$data)) || is.null(GWPCR$data[[m]]))
    gwpcr.molecules.precompute(molecules=m)

  # Pre-process data by calling Akima's procedure with MD=1.
  if (is.null(GWPCR$data.akima.wk))
    GWPCR$data.akima.wk <- list()
  if ((m > length(GWPCR$data.akima.wk)) || is.null(GWPCR$data.akima.wk[[m]])) {
    #stopifnot(is.double(GWPCR$efficiency) && is.double(GWPCR$lambda))
    #stopifnot(is.double(GWPCR$data[[m]]))
    #stopifnot(length(GWPCR$efficiency) == nrow(GWPCR$data[[m]]))
    #stopifnot(length(GWPCR$lambda) == ncol(GWPCR$data[[m]]))
    r <- .Fortran("rgbi3p",
                  md=as.integer(1),
                  nxd=length(GWPCR$efficiency), nyd=length(GWPCR$lambda),
                  xd=GWPCR$efficiency, yd=GWPCR$lambda, zd=GWPCR$data[[m]],
                  nip=as.integer(0),
                  xi=double(), yi=double(), zi=double(),
                  ier=integer(1),
                  wk=double(3*length(GWPCR$efficiency)*length(GWPCR$lambda)),
                  PACKAGE="akima", NAOK=TRUE)
    GWPCR$data.akima.wk[[m]] <- r$wk
  }

  #stopifnot(is.double(GWPCR$efficiency) && is.double(GWPCR$lambda))
  #stopifnot(is.double(GWPCR$data[[m]]))
  #stopifnot(length(GWPCR$efficiency) == nrow(GWPCR$data[[m]]))
  #stopifnot(length(GWPCR$lambda) == ncol(GWPCR$data[[m]]))
  #stopifnot(is.double(GWPCR$data.akima.wk[[m]]))
  #stopifnot(length(GWPCR$data.akima.wk[[m]]) == 3*length(GWPCR$efficiency)*length(GWPCR$lambda))
  #stopifnot(is.double(x0) && is.double(y0))
  #stopifnot(length(x0) == length(y0))
  r <- .Fortran("rgbi3p",
                md=as.integer(2),
                nxd=length(GWPCR$efficiency), nyd=length(GWPCR$lambda),
                xd=GWPCR$efficiency, yd=GWPCR$lambda, zd=GWPCR$data[[m]],
                nip=length(x0),
                xi=x0, yi=y0, zi=double(length(x0)),
                ier=integer(1),
                wk=double(3*length(GWPCR$efficiency)*length(GWPCR$lambda)),
                PACKAGE="akima", NAOK=TRUE)
  return(pmax(r$zi, 0))
}
