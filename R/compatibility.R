# estimate.R, Copyright 2016,2017 Florian G. Pflug
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

#' Compatibility wrapper of \code{\link{gwpcrpois.est}}
#'
#' @seealso \code{\link{gwpcrpois.est}}
#'
#' @export
gwpcrpois.mom <- function(mean, var, threshold=1, molecules=1,
                          ctrl=list(), nonconvergence.is.error=FALSE)
{
  gwpcrpois.est(mean=mean, var=var, threshold=threshold, molecules=molecules,
                must.converge=nonconvergence.is.error, ctrl=ctrl)
}

#' Compatibility wrapper of \code{\link{gwpcrpois.est}}
#'
#' @seealso \code{\link{gwpcrpois.est}}
#'
#' @export
gwpcrpois.mle <- function(c, threshold=1, molecules=1) {
  gwpcrpois.mle(c, threshold=threshold, molecules=molecules)
}

#' Compatibility wrapper of \code{\link{gwpcrpois.groupest}}
#'
#' Uses the \emph{method of moments} and sets the \var{ctrl} parameters to
#' their previous defaults, i.e. \code{include.mean.var=TRUE},
#' \code{obs.min.ingroup=2} and \code{use.nonconv.groupest=TRUE}.
#'
#' @seealso \code{\link{gwpcrpois.groupest}}
#'
#' @export
gwpcrpois.mom.groupwise <- function(formula, data, threshold=1, molecules=1,
                                    loss=expression(p0), ctrl=list()) {
  # Backwards-compatible defaults
  if (!('include.mean.var' %in% names(ctrl)))
    ctrl$include.mean.var <- TRUE
  if (!('obs.min.ingroup' %in% names(ctrl)))
    ctrl$obs.min.ingroup <- 2
  if (!('use.nonconv.groupest' %in% names(ctrl)))
    ctrl$use.nonconv.groupest <- TRUE

  gwpcrpois.groupest(formula, data, method='mom', threshold=threshold,
                     molecules=molecules, ctrl=ctrl)
}
