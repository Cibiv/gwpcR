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

E.GAMMA.TH <- 1e-1

delayedAssign('E.MIN', head(GWPCR$efficiency, 1))

delayedAssign('E.MAX', tail(GWPCR$efficiency, 1))

gamma.factor <- function(efficiency) {
  c <- pmin(pmax(0, (efficiency - E.MIN)/(E.GAMMA.TH - E.MIN)), 1)
  f <- ifelse(c <= 0.5, 1 - 2*c^2, 2*(1 - c)^2)
}
