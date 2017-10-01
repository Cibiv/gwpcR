# utility.R, Copyright 2016,2017 Florian G. Pflug
#
# This file is part of gwpcR
#
# Foobar is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
# Foobar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

#' @importFrom data.table data.table :=
handle.parameters <- function(parameters, by, expr) {
  n <- 0
  for(p in names(parameters)) {
    if (!is.numeric(parameters[[p]]) || (length(parameters[[p]]) == 0))
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

#' @useDynLib gwpcR gwpcr_refine_c
refine <- function(points, width) .Call(gwpcr_refine_c, points, width)

# Experimentally determined by KS-testing rgwpcr against pgwpcr. For
# efficiency 0.01, the KS-Test found essentially no difference at 1e6 samples.
# At efficiency 0.02, both the interpolated simulation results and the gamma
# approximation yield reduced p-values, but the simulation results are still
# much better. This was thus chosen as as the cutoff point.
E.GAMMA.TH <- 0.02

# These are just for convenience and readability.
delayedAssign('E.MIN', head(GWPCR$efficiency, 1))
delayedAssign('E.MAX', tail(GWPCR$efficiency, 1))

# This the weight factor use for slow cutover to the gamma approximation of
# the gwpcr distribution. It is supposed to be C^1-smooth, and takes
# value 1 at E.MIN and below, and 0 at E.GAMMA.TH and above, and increases
# monotonically inbetween.
gamma.factor <- function(efficiency) {
  c <- pmin(pmax(0, (efficiency - E.MIN)/(E.GAMMA.TH - E.MIN)), 1)
  f <- ifelse(c <= 0.5, 1 - 2*c^2, 2*(1 - c)^2)
}
