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

#' Parameter Estimation for PCR-Poisson Mixture
#'
#' Estimates the parameters \var{efficiency} and \var{lambda0} from a vector
#' of observed read counts per molecular family, or (depending on the estimation
#' method) the mean and variance of these observations. Supports arbitrary
#' detection thresholds and initial molecule counts, but estimation may be
#' considerably faster in the (unrealistic) case \var{threshold=0} than in the
#' general one.
#'
#' @param x a numeric vector containing the observed read counts per molecular
#'   family, after removal of families below the detection threshold. The
#'   vector must thus contain only whole numbers not smaller than \var{threshold}.
#'   This parameter and the combination of parameters \var{mean} and \var{var}
#'   are mutually exclusive.
#'
#' @param mean average number of observations per molecular family
#'   \emph{computed over the unambiguously detected famililies}, i.e. over those
#'   families which were observed at least \var{threshold} times. This parameter
#'   and specifying an observation vector through parameter \var{x} are mutually
#'   exclusive, and specifying \var{mean} and \var{var} instead of the full
#'   observation vector is only possible for estimation method 'mom'.
#'
#' @param var standard deviations of number of observations per molecular
#'   family, also \emph{computed over the unambiguously detected famililies},
#'   i.e. over those families which were observed at least \var{threshold} times.
#'   This parameter and specifying an observation vector through parameter \var{x}
#'   are mutually exclusive, and specifying \var{mean} and \var{var} instead of
#'   the full observation vector is only possible for estimation method 'mom'.
#'
#' @param n.umis the number of observed \acronym{UMI}s, used in the estimation
#'    of \code{n.tot} (i.e. the total number of molecules/\acronym{UMI}s in the
#'    original sample). See also the dicussion of the parameter \code{loss}, and
#'    the result value \code{n.tot}.
#'
#' @param method the estimation method to use, either 'mle' for \emph{maximum
#'   likelihood estimation} or 'mom' for \emph{method of moments}. See Details.
#'
#' @param must.converge if set to \code{TRUE}, an error is reported of the
#'   parameter estimation fails to converge. If \code{FALSE}, a warning is
#'   reported instead.
#'
#' @param ctrl a list of settings controlling the estimation procedure.
#'   Difference estimation methods recognize different possible \var{ctrl}
#'   settings, unrecognized settings are ignored without warning. See Details
#'   for the settings relevant to each estimation method.
#'
#' @inheritParams gwpcrpois
#'
#' @param loss an expression specifying how the loss, i.e. the percentage of
#'   all molecules (or UMIs) that was not observed, or removed by the read
#'   count threshold. In the simple case of each read count observation
#'   representing a separate molecules, the default value \var{p0} is correct
#'   -- the lost molecules are then simply those whose read count lies below
#'   the specified \var{threshold}. In more complex scenarios, e.g. if a
#'   single molecule produces separate read count for each strand, which are
#'   then either both rejected or both accepted, the additional rejection cases
#'   must be considered by a custom loss expression
#'
#' @param ctrl a list of settings influecing the numerical method(s) used for
#'   estimating the parameters. The behaviour and supported \var{ctrl} settings
#'   depend on the estimation method used.
#'
#' @return A list containing the values
#'
#'   \item{convergence}{flag indicating whether the estimation converged.
#'   \var{0} indicates convergence.}
#'
#'   \item{efficiency}{parameter estimate for \var{efficiency} (see
#'   \link{gwpcrpois})}
#'
#'   \item{lambda0}{parameter estimate for \var{lambda0} (see \link{gwpcrpois})}
#'
#'   \item{p0}{probability of observing a read count less than the specified
#'   \var{threshold}}
#'
#'   \item{loss}{the estimation loss according to the specified loss expression,
#'    i.e. percentage of molecules not observed or filtered out}
#'
#'   \item{n.tot}{the estimated total number of molecules in the sample, i.e.
#'   \code{n.ums / (1 - loss)}}
#'
#'   \item{n.obs}{The length of the observation vector \var{x} used for parameter
#'   estimation, or \code{NA} if \var{mean} and \var{var} were specified
#'   directly.}
#'
#'   \item{n.umis}{The number of observed molecules/\acronym{UMI}s specified in
#'   the call to \code{gwpcrpois.est}}
#'
#'   \item{threshold}{detection threshold specified in the call to
#'   \code{gwpcrpois.est}}
#'
#'   \item{molecules}{initial molecule count specified in the call to
#'   \code{gwpcrpois.est}}
#'
#' @details \describe{
#'   \item{\code{mom}}{
#'   For the (unrealistic) uncensored case, i.e. \var{threshold=0}, the
#'   specified mean is the method-of-moments estimate for \var{lambda0}, and a
#'   closed formula is used to compute \var{efficiency} from \var{mean} and
#'   \var{var}. The \var{ctrl} argument is not used in this case.
#'
#'   In case of a censored distribution, the sample mean is not a consistent
#'   estimator for \var{lambda0} because the expectation of the censored
#'   distribution is in general larger than \var{lambda0}. The sample variance
#'   simiarly deviates from the variance of the uncensored distribution.
#'
#'   An interative approach is used to find method-of-moment estimates in this
#'   case. Initial estimates are computed as if \var{threshold} were zero. From
#'   these the probability \var{pdetect} of detecting a particular family is
#'   found, and used to correct for the biases in the sample mean and variance.
#'   Then the parameter estimates are updated. This process continues until it
#'   either converges or reaches the maximum allowed number of iterations. Both
#'   termination criteria can be controlled via the \var{ctrl} parameter, which
#'   is a list that can contain the following components: \describe{
#'
#'   \item{maxit}{Maximum number of iterations. Defaults to 150.}
#'
#'   \item{rel.tol}{Relative convergence tolerance. Applied to \var{efficiency},
#'   \var{lambda0} and the detection probability \var{pdetect}. Defaults to
#'   \code{1e-4}.}
#'
#'   \item{rel.tol}{Absolute convergence tolerance. Only used as the tolerance
#'   around zero, where the relative tolerance becomes meaningless. Defaults to
#'   \code{1e-4}.}
#'
#'   \item{trace}{Output estimates after each round}
#'   }
#'   }
#'   \item{\code{mle}}{
#'   The parameters are estimated by maximizing the log-likelihood using
#'   \code{\link{optim}}. The \var{ctrl} settings are passed through to
#'   \code{\link{optim}}, except for \var{fnscale} and \var{parscale} which
#'   are overwritten. The method of moments (method \code{'mom'})) estimate
#'   is used as the starting point during liklihood optimization.
#'   }
#'   }
#'
#' @seealso \code{\link{gwpcrpois}}
#'
#' @export
gwpcrpois.est <- function(x=NULL, mean=NULL, var=NULL, n.umis=NULL, method="mom",
                          must.converge=TRUE, threshold=1, molecules=1, loss=expression(p0),
                          ctrl=list())
{
  if (!is.null(x) && !is.numeric(x))
    stop('if specified, observation vector x must be a numeric vector')
  if (!is.null(mean) && (!is.numeric(mean) || (length(mean) != 1) || (mean < 0)))
    stop('if specified, mean must be a non-negative scalar')
  if (!is.null(var) && (!is.numeric(var) || (length(var) != 1) || (var < 0)))
    stop('if specified, var must be a non-negative scalar')
  method <- match.arg(method, c("mom", "mle"))
  if (!(method %in% c("mom", "mle")))
    stop("method must be 'mom' or 'mle'")
  if (!is.logical(must.converge) || (length(must.converge) != 1) || is.na(must.converge))
    stop('nonconvergence.is.error must be TRUE or FALSE')
  if (!is.numeric(threshold) || (length(threshold) != 1) || (threshold != floor(threshold)) || (threshold < 0))
    stop('threshold must be a non-negative integral scalar')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a strictly positive integral scalar')
  if (!is.expression(loss) || (length(loss) != 1))
    stop('loss must be an expression')
  if (!is.list(ctrl))
    stop('ctrl must be a list')

  # Call either gwpcrpois.estimator.mom or gwpcrpois.estimator.mle
  r <- if (is.numeric(x)) {
    if (!is.null(mean) || !is.null(var))
      stop("either mean and var OR a vector or observations, but not both, must be specified")
    # x is a vector of per-UMI read counts (after threshold application)
    if (any(x != floor(x)) || any(is.na(x)) || any(!is.finite(x)) || any(x < threshold))
      stop("observations in vector x must be whole numbers >= threshold")
    if (method == "mom")
      gwpcrpois.estimator.mom(mean=mean(x), var=var(x), threshold=threshold,
                              molecules=molecules, ctrl=ctrl)
    else if (method == "mle")
      gwpcrpois.estimator.mle(x, threshold=threshold, molecules=molecules, ctrl=ctrl)
  } else if (is.null(x) && is.numeric(mean) && is.numeric(var)) {
    if (method != "mom")
      stop("if mean and var instead of the full observation vector is specified, only method 'mom' is supported")
    gwpcrpois.estimator.mom(mean=mean, var=var, threshold=threshold,
                            molecules=molecules, ctrl=ctrl)
  } else {
    stop("either mean and var or a vector of observations must be specified")
  }

  # Handle nonconvergence.is.error
  if (r$convergence != 0) {
    if (nonconvergence.is.error)
      stop("gwpcrpois.mom did not converge")
    else
      warning("gwpcrpois.mom did not converge, returning best estimate so far")
  }

  # Evaluate loss expression
  r$loss <- tryCatch(eval(loss, envir=r, enclos=parent.frame()),
                     error=function(e) { warning(conditionMessage(e)) ; NA})

  # Set n.obs, n.umis and n.tot = n.umis / (1 - loss)
  r$n.obs <- if (!is.null(x)) length(x) else NA
  r$n.umis <- if (!is.null(n.umis)) n.umis else if (!is.null(x)) length(x) else NA
  r$n.tot <- r$n.umis / (1 - r$loss)

  return(r)
}

#' Group-wise Parameter Estimation for PCR-Poisson Mixture
#'
#' @seealso \code{\link{gwpcrpois}}
#'
#' @export
gwpcrpois.est.groups <- function(formula, data, method="mom", threshold=1, molecules=1,
                                 loss=expression(p0), ctrl=list())
{
  if (!class(formula) == "formula")
    stop("formula must be a 'formula', see help(formula)")
  data <- as.data.table(data)
  method <- match.arg(method, c("mom", "mle"))
  if (!(method %in% c("mom", "mle")))
    stop("method must be 'mom' or 'mle'")
  if (!is.numeric(threshold) || (length(threshold) != 1) || (threshold != floor(threshold)) || (threshold < 0))
    stop('threshold must be a non-negative integral scalar')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a strictly positive integral scalar')
  if (!is.expression(loss) || (length(loss) != 1))
    stop('loss must be an expression')
  if (!is.list(ctrl))
    stop('ctrl must be a list')

  # Convert to "frame", i.e. a data.table with group key columns, "n.umis" and "count".
  frame <- as.gwpcrpois.frame(formula, data)
  if (nrow(frame[!is.finite(count) | (count < threshold)]) > 0)
    stop("data contains non-finite counts or counts below the threshold")

  # Run shrinkage estimator
  r <- gwpcrpois.estimator.groups.shrink(frame, method, threshold=threshold,
                                         molecules=molecules, loss=loss, ctrl=ctrl)

  # Correct for unobserved UMIs, i.e. compute
  #            n
  #   n   = ----------
  #    tot   1 - loss
  if (verbose)
    message("Computing n.tot")
  r[, n.tot := ifelse(is.finite(loss), n.umis / (1 - loss), NA) ]

  return(r)
}

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

# ***************************************************************************************
# Method of Moments (MoM) Estimator
# ***************************************************************************************

gwpcrpois.estimator.mom <- function(mean, var, threshold, molecules, ctrl) {
  # If the threshold is non-zero, we resort to a numerical solution
  exact.solution <- (threshold == 0)
  
  ctrl.get <- function(key, default) {
    r <- ctrl[[key]]
    if (!is.null(r))
      r
    else
      default
  }
  rel.tol <- as.numeric(ctrl.get("rel.tol", 1e-4))
  abs.tol <- as.numeric(ctrl.get("abs.tol", 1e-4))
  maxit <- as.integer(ctrl.get("maxit", 150))
  trace <- as.logical(ctrl.get("trace", FALSE))
  initial <- as.list(ctrl.get("initial", NULL))
  within.tol <- function(old, new) {
    if (abs(new) < abs.tol)
      return(TRUE)
    else if (abs((new - old) / new) < rel.tol)
      return(TRUE)
    else
      return(FALSE)
  }

  # A random variable C distributed according to the PCR-Poisson mixture
  # with parameters E (efficiency) and lambda0 has mean
  #   E( C | E, lambda0 ) = lambda0,
  # and variance
  #   V( C | E, lambda0 ) = lambda0 + lambda0^2 * V( L | E ),
  # where L is distributed according to the PCR distribution with efficiency E.
  # A method-of-moments estimator for the parameters of C is thus to set
  #   lambda0 = Avg(C),
  # and to set
  #   v' = ( Var(C) - lambda0 ) / lambda0^2
  # and find E such that
  #   V( L | E ) = v',
  # which is what gwpcr.sd.inv does
  if (exact.solution || (length(initial) == 0)) {
    pdetect <- 1
    lambda0 <- mean
    efficiency <- gwpcr.sd.inv(sqrt(max((var - lambda0) / ( lambda0^2 ), 0)),
                               molecules=molecules)
  } else {
    # If the user specified an initial value, use that
    pdetect <- 1 - initial$p0
    lambda0 <- initial$lambda0
    efficiency <- initial$efficiency
  }
    
  # If the data is censored, i.e. if only counts >= threshold > 0 are
  # observed, the sample mean will over-estimate lambda0, and the sample
  # variance also won't reflect the uncensored distribution's variance.
  #
  # Let c_m be the expected sampling mean, i.e.
  #   c_m = E ( C | C >= TH, E, lambda0 ),
  # and v_m the expected sampling variance, i.e.
  #   v_m = V ( C | C >= TH, E, lambda0 ),
  # and set:
  #   m_0 = Sum c   * P( C=c | E, lambda0 ) for 0 <= c < TH,
  #   v_0 = Sum c^2 * P( C=c | E, lambda0 ) for 0 <= c < TH,
  #   p_0 =       Sum P( C=c | E, lambda0 ) for 0 <= c < TH,
  # then the following relationship holds:
  #   (A)  E( C | E, lambda0 ) = lambda_0 = m_0 + (1 - p_0) * c_m,
  #   (B)  V( C | E, lambda0 ) = v = v_0 + (1 - p_0) * ( v_m + c_m^2 ) - lambda_0^2.
  # Together with the relationship between v and E from above, i.e.
  #   (C)  V( L | E ) = ( v - lambda0 ) / lambda0^2
  # this yields a system of equations in E and lambda with parameters
  # c_m, v_m and TH. The code below solves this iteratively by starting with
  # the estimates for lambda0 and E for TH=0 from abovel. (A) is then used
  # to update lambda0, (B) to compute a new v, and (C) to find the updated (E).
  #
  # Experiments showed that it is beneficial to use the updated lambda0 when
  # evaluting (B) and (C). However, for performance reasons, m_0, v_0 and p_0
  # are not re-evaluated immediately after updating lambda0, but instead the
  # old values to used to update the efficiency. The values ARE then re-computed
  # during the next round, however,
  if (!exact.solution) {
    stop <- FALSE
    i <- 0
    if (trace)
      print(cbind(iteration=i, efficiency=efficiency,
                  lambda0=lambda0, pdetect=pdetect))
    while (!stop) {
      x <- seq(from=0, to=threshold-1, by=1)
      d <- dgwpcrpois(x, efficiency=efficiency, lambda0=lambda0,
                      threshold=0, molecules=molecules)
      pdetect.p <- (1 - sum(d))
      lambda0.p <- sum(x * d)   + pdetect.p * mean
      v         <- sum(x^2 * d) + pdetect.p * (var + mean^2) - lambda0.p^2
      efficiency.p <- gwpcr.sd.inv(sqrt(max((v - lambda0.p) / ( lambda0.p^2 ), 0)),
                                   molecules=molecules)
      i <- i + 1

      if (within.tol(lambda0, lambda0.p) && within.tol(efficiency, efficiency.p) &&
          within.tol(pdetect, pdetect.p)) {
        stop <- TRUE
        convergence <- 0
      } else if (i >= maxit) {
        stop <- TRUE
        convergence <- 1
      }

      lambda0 <- lambda0.p
      efficiency <- efficiency.p
      pdetect <- pdetect.p

      if (trace)
        print(cbind(iteration=i, efficiency=efficiency,
                    lambda0=lambda0, pdetect=pdetect))
    }
  } else {
    pdetect <- 0
    convergence <- 0
  }

  if (trace)
    cat(paste0("Convergence: ", if (convergence == 0) "Yes" else "No"))

  # Return results
  return(list(lambda0=lambda0, efficiency=efficiency, p0=1-pdetect,
              convergence=convergence, threshold=threshold, molecules=molecules))
}

# ***************************************************************************************
# Maximum Likelihood (ML) Estimator
# ***************************************************************************************

gwpcrpois.estimator.mle <- function(c, threshold, molecules, ctrl) {
  # Since evaluating the PCR-Poisson mixture is slow, we optimize by evaluating it only
  # once for each unique observed count, and multiplying the log-likelihood with the
  # number of times that count was observed.
  v <- rle(sort(c))

  # Use method of moments estimates as initial parameters. Since we do the parameter
  # search with clamp.efficiency set to FALSE, we must take care to clamp it to the
  # range of efficiencies found in GWPCR here
  mom <- gwpcrpois.estimator.mom(mean=mean(c), var=var(c), must.converge=FALSE,
                                 molecules=molecules, threshold=threshold)
  mom$efficiency <- pmin(pmax(GWPCR$efficiency[1], mom$efficiency), tail(GWPCR$efficiency,1))

  # Optimize likelihood
  # We set parscale such that on the par/parscale scale that optim uses, a
  # step of minus one away from the initial parameters corresponds to a decrease
  # of 10%.  fnscale is set such that on the fn/fnscale scale of optim, such a step
  # corresponds approximately to a unit change of the goal function. Remeber that
  # optim minimizes, so to maximize fnscale contains an additional factor "-1".
  logl <- function(p) {
    e <- p['efficiency']
    l <- p['lambda0']
    # If parameters are valid, compute log-likelihood, otherwise
    # return NA.
    if ((e >= 0) && (e <= 1) && (l > 0))
      sum(log(dgwpcrpois(c=v$values, efficiency=e, lambda0=l,
                         threshold=threshold, molecules=molecules)) * v$lengths)
    else
      as.numeric(NA)
  }
  p0 <- c(efficiency=mom$efficiency, lambda0=mom$lambda0)
  ctrl$fnscale <- -abs(logl(p0) - logl(p0 * c(0.9, 0.9)))
  ctrl$parscale <- c(efficiency=p0['efficiency']/10, lambda0=p0['lambda0']/10)
  r <- optim(par=p0, fn=logl, method="Nelder-Mead", control=ctrl)

  # Compute p0
  p0 <- if (threshold > 0)
    pgwpcrpois(threshold-1, efficiency=r$par['efficiency'], lambda0=r$par['lambda0'],
               threshold=0, molecules=molecules)
  else
    0

  # Return result
  list(lambda0=as.vector(r$par['lambda0']),
       efficiency=as.vector(r$par['efficiency']),
       p0=p0,
       convergence=r$convergence,
       threshold=threshold,
       molecules=molecules)
}

# ***************************************************************************************
# Data Grouping
# ***************************************************************************************
as.gwpcrpois.frame <- function(formula, data) {
  # Subset data to have a single "counts" columns, follows by the key columns defining
  # the groups (e.g. the gene name, but there can be more than one, e.g. sample and gene).
  # If the LHS of the formula contains an expressionlike c(col1, col2), the resulting
  # table will have more rows than "data". This allows counts in multiple columns to be
  # treated together (e.g. read counts for the plus and minus strand).
  formula.t <- terms(formula)
  group.key <- labels(formula.t)
  counts.expr <- attr(formula.t, "variables")[[1 + attr(formula.t, "response")]]
  list(group.key=group.key, frame=data[, list(n.umis=.N, count=eval(counts.expr)), keyby=group.key])
}

gwpcrpois.grouped.frame <- function(frame) {
  frame[, list(), keyby=key(frame)]
}

# ***************************************************************************************
# Raw Groupwise Estimator
# ***************************************************************************************
gwpcrpois.estimator.groups.raw <- function(frame, frame.grp=NULL, method, threshold,
                                           molecules, loss=expression(p0), ctrl=list())
{
  # Get control parameters
  cores <- as.numeric(ctrl.get("cores", 1))
  my.lapply <- if (cores > 1) {
    function(...) { mclapply(..., mc.cores=cores, mc.preschedule=TRUE, mc.cleanup=TRUE) }
  } else {
    lapply
  }
  use.nonconv.groupest <- as.logical(ctrl.get("use.nonconv.groupest", FALSE))
  include.mean.var <- as.logical(ctrl.get("include.mean.var", FALSE))
  verbose <- as.logical(ctrl.get("verbose", FALSE))

  # Generate frame with one raw per group (unless one was provided)
  if (is.null(frame.grp))
    frame.grp <- gwpcrpois.grouped.frame(frame)

  # Compute raw (group-wise) parameter estimates
  if (verbose)
    message("Computing group-wise parameter estimates on ", cores, " cores")
  group.key <- key(frame.grp)
  groups <- do.call(mapply, c(list(list), frame.grp[, ..group.key], list(SIMPLIFY=FALSE)))
  model.na <- list(lambda0=NA, efficiency=NA, p0=NA, loss=NA)
  frame.grp <- rbindlist(my.lapply(groups, function(k) {
    # Fetch observations for group
    r <- frame.grp[k,]

    # Fit model
    obs <- frame[k, count]
    m <- if (length(obs) >= obs.min.ingroup)
      tryCatch(gwpcrpois.est(obs, must.converge=!use.nonconv.groupest, loss=loss,
                             threshold=threshold, molecules=molecules, ctrl=list()),
               error=function(e) {
                 cat(paste0("Failed to estimate parameters for group (",
                            paste0(lapply(k, as.character), collapse=","),
                            "): ", conditionMessage(e), "\n"), file = stderr())
                 model.na })
    else
      model.na

    # Backwards-compatibility hack: Insert raw group mean and variance for method of moments
    if (include.mean.var)
      r[, c("mean.raw", "var.raw") := list(mean(obs), var(obs)) ]

    # Generate row
    r[, c("n.umis", "n.obs", "efficiency.raw", "lambda0.raw", "loss.raw") :=
        list(frame[k, n.umis], length(obs), m$efficiency, m$lambda0, m$loss)]
  }))
  setkeyv(frame.grp, group.key)

  # Return grouped frame
  return(frame.grp)
}

# ***************************************************************************************
# Shrinkage Estimator
# ***************************************************************************************
gwpcrpois.estimator.groups.shrink <- function(frame, method, threshold, molecules,
                                              loss, ctrl=list())
{
  verbose <- as.logical(ctrl.get("verbose", FALSE))
  var.est.distfree <- as.logical(ctrl.get("var.est.distfree", TRUE))
  include.mean.var <- as.logical(ctrl.get("include.mean.var", FALSE))

  # Compute global parameter estimates
  if (verbose)
    message("Computing overall parameter estimates")
  mean.all = frame[, mean(count)]
  var.all = frame[, var(count)]
  m.all <- gwpcrpois.est(frame$count, method=method, must.converge=!use.nonconv.globalest,
                         loss=loss, threshold=threshold, molecules=molecules, ctrl=ctrl)

  # Group data, add global estimate as columns for convenience
  frame.grp <- gwpcrpois.grouped.frame(frame)
  if (include.mean.var)
    frame.grp[, c("mean.all", "var.all") := list(mean(frame$count), var(frame$count)) ]
  frame.grp[, c("efficiency.all", "lambda0.all", "loss.all") :=
              list(m.all$efficiency, m.all$lambda0, m.all$loss) ]

  # Estimate raw group-wise parameters efficiency.raw, lambda0.raw, loss.raw
  frame.grp <- gwpcrpois.estimator.groups.raw(frame, frame.grp, method, threshold=threshold,
                                              molecules=molecules, loss=loss, ctrl)

  # Estimate variance of loss^raw (and also lambda0^raw, efficiency^raw) as a function
  # of the number of observed UMIs n.
  #
  # This is the distribution-free version of the estimator. We assume the mean loss is
  # independent from n^obs. For a fixed n^obs, the expected value of (loss^raw - m)^2
  # is then v_g + v_e / n^obs. We find v_g and v_e by mininizing the squared deviation
  # of the (single-point) variance estimate (loss^raw - m)^2 and the predictor. Since
  # variances are generally bigger for small n^obs, minimizing an unweighted sum of
  # square deviations would allow loss estimates for small n^obs to influence the
  # estimation unduely strongly. We correct for this by weigthing the squared deivations
  # with n^obs (which matches the expected drop in scale), and thus minimize the sum
  #
  #                                               2
  #    /                                        \
  #   |  (loss^raw - m)^2  -  v_g - v_e / n^obs |   * w(n^obs), where w(n) = n / (1 + n/W)
  #    \                                        /
  #
  # over all genes. m is the sample mean over all genes. We start the numerical
  # search with v_g = v_e = v/2, where v is the sampling variance of the loss.
  variance.estimates.distfree <- function(x.name) {
    # Parameter W controls the grow behaviour of weights. For small n.obs, the weights
    # equal n.obs, but as n.obs grows, weights eventually converge to W.
    W <- 100
    if (verbose)
      message("Estimating variances of ", x.name, " using the distribution-free algorithm")
    x <- parse(text=paste0(x.name, ".raw"))
    # Find least-squared estimates for v_g, v_e.
    r <- frame.grp[is.finite(eval(x)) & (n.obs > 0), {
      y <- (eval(x) - mean(eval(x), na.rm=TRUE))^2
      x <- cbind(v_g=rep(1, .N), v_e=1/n.obs)
      w <- n.obs / (1 + n.obs / W)
      lsfit(x, y, wt=w, intercept=FALSE)$coefficients
    } ]
    if (all(r < 0)) {
      stop("both variance estimators are negative")
    } else if (r["v_g"] < 0) {
      warning("variance of ", x.name, " between groups estimated to be negative, setting to zero")
      # Between-gene variance estimate is negative. Assume between-gene variance
      # is negligible, and estimate only the estimator variance.
      r <- c(v_g=0, frame.grp[is.finite(eval(x)) & (n.obs > 0), {
        y <- (eval(x) - mean(eval(x), na.rm=TRUE))^2
        x <- cbind(v_e=1/n.obs)
        w <- n.obs
        lsfit(x, y, wt=w, intercept=FALSE)$coefficients["v_e"]
      } ])
    } else if (r["v_e"] < 0) {
      warning("estimator variance of ", x.name, " within groups estimated to be negative, setting to zero")
      # Estimate of estimator variance is negative. Assume estimator variance is
      # negligible, and estimate only the between-gene variance.
      r <- c(v_g=frame.grp[is.finite(eval(x)) & (n.obs > 0), var(eval(x), na.rm=TRUE)], v_e=0)
    }
    if (any(r < 0))
      stop("some variance estimators are negative")
    # Add columns
    frame.grp[, paste0(x.name, ".grp.var") := r["v_g"]]
    frame.grp[, paste0(x.name, ".raw.var") := r["v_e"] / n.obs]
  }

  # Estimate variance of p_0^raw (and also lambda0^raw, efficiency^raw) as a function
  # of the number of observed UMIs n.
  #
  # We assume that p_0^raw | n ~ Beta( a(n), b(n) ), with a fixed expectation m and
  # n-dependent variance comprising a estimation error that decreases with 1/n and a
  # constant part that reflects the variance of the loss between groups,
  #
  #        raw
  #    V( p    | n) = v(n) = v + v  /  n.
  #        0                  g   e
  #
  # From the mean (assumed constant) and variance of the Beta distribution it follows that
  #
  #                                                      m * (1 - m)
  #   a(n) = m * f(n), b(n) = (1-m) * f(n) where f(n) = ------------- - 1.
  #                                                           v
  #                                                            e
  #                                                       v + ---
  #                                                        g   n
  #
  # For m we used the sample mean over all genes, and for v_g and v_e ML estimates
  # found using numerical optimzation. We start the numerical search with
  # v_g = v_e = v/2, where v is the sampling variance of the loss over all genes (limited
  # to the possible range (0, m * (1 - m) ).
  #
  # efficiency^raw is handled similarly. For lambda0^raw, the beta distribution
  # is replaced by a normal distribution.
  variance.estimates.beta <- function(x.name) {
    if (verbose)
      message("Estimating variances of ", x.name, " assuming a beta distribution")
    x <- parse(text=paste0(x.name, ".raw"))
    # For the mean use the sampling mean, since we assume it's independent of n
    m <- frame.grp[, mean(eval(x), na.rm=TRUE) ]
    # Compute quantity to minimize, i.e. negative log-liklihood
    # We contract "x" a bit towards 0.5 to ensure that all values are strictly
    # greater than 0 and less than 1 -- otherwise, the log-likelihood becomes -Inf.
    # M is the contraction factor, i.e. we move X to the interval [M/2, 1-M/2]
    M <- 1e-6
    logl <- function(p) {
      # Reject invalid parameters. Note that the beta variance is always <= m * (1-m).
      if (any(p <= 0) || (p["v_g"] >= m * (1-m)))
        return(NA)
      # Evaluate likelihoods. We strict shape1,2 to >= 1 to ensure monomodality
      -sum(frame.grp[is.finite(eval(x)) & (n.obs > 0), {
        f <- pmax(m * (1-m) / (p["v_g"] + p["v_e"] / n.obs) - 1, 1/m, 1/(1-m))
        dbeta(M/2 + eval(x)*(1-M), shape1=m*f, shape2=(1-m)*f, log=TRUE)
      }])
    }
    # Optimize v_g and v_e. Initially, we split the total observed variance evenly
    # between v_g and v_e / n.obs for the smallest positive group size n.obs.
    v <- frame.grp[, min(var(eval(x), na.rm=TRUE), m * (1-m) ) ]
    n.min <- frame.grp[n.obs > 0, min(n.obs) ]
    r <- optim(fn=logl, par=c(v_g=v/2, v_e=n.min*v/2), method="Nelder-Mead")
    # Correct variances for the previous contraction
    v_g <- r$par["v_g"] / (1-M)^2
    v_e <- r$par["v_e"] / (1-M)^2
    # Add columns
    frame.grp[, paste0(x.name, ".raw.var") := v_e / n.obs]
    frame.grp[, paste0(x.name, ".grp.var") := v_g]
  }
  variance.estimates.normal <- function(x.name) {
    if (verbose)
      message("Estimating variances of ", x.name, " assuming a normal distribution")
    x <- parse(text=paste0(x.name, ".raw"))
    m <- frame.grp[, mean(eval(x), na.rm=TRUE) ]
    v <- frame.grp[, var(eval(x), na.rm=TRUE) ]
    # Compute quantity to minimize, i.e. negative log-liklihood
    logl <- function(p) {
      # Reject invalid parameters
      if (any(p < 0))
        return(NA)
      # Evaluate likelihoods
      -sum(frame.grp[is.finite(eval(x)) & (n.obs > 0),
                     dnorm(eval(x), mean=m, sd=sqrt(p["v_g"] + p["v_e"] / n.obs), log=TRUE)
                     ])
    }
    r <- optim(fn=logl, par=c(v_g=v/2, v_e=v/2), method="Nelder-Mead")
    # Add columns
    frame.grp[, paste0(x.name, ".raw.var") := r$par["v_e"] / n.obs ]
    frame.grp[, paste0(x.name, ".grp.var") := r$par["v_g"] ]
  }
  if (var.est.distfree) {
    variance.estimates.distfree("loss")
    variance.estimates.distfree("efficiency")
    variance.estimates.distfree("lambda0")
  } else {
    variance.estimates.beta("loss")
    variance.estimates.beta("efficiency")
    variance.estimates.normal("lambda0")
  }

  # Compute shrunken estimates of loss, efficiency and lambda0 i.e.
  #
  #         all   v_e      raw
  #        p    * ---  +  p    * v_g
  #         0      n       0
  #   p = -----------------------------------.
  #    0          v_e
  #               ---  +  v_g
  #                n
  #
  # and similarly for the other two quantities
  if (verbose)
    message("Computing final parameter estimates")
  frame.grp[, loss := (loss.all * loss.raw.var + loss.raw * loss.grp.var) / (loss.raw.var + loss.grp.var) ]
  frame.grp[, efficiency := (efficiency.all * efficiency.raw.var + efficiency.raw * efficiency.grp.var) / (efficiency.raw.var + efficiency.grp.var) ]
  frame.grp[, lambda0 := (lambda0.all * lambda0.raw.var + lambda0.raw * lambda0.grp.var) / (lambda0.raw.var + lambda0.grp.var) ]
  # If the local estimate is NA, use the global one
  frame.grp[!is.finite(loss), loss := loss.all ]
  frame.grp[!is.finite(efficiency), efficiency := efficiency.all ]
  frame.grp[!is.finite(lambda0), lambda0 := lambda0.all ]

  return(frame.grp)
}

# ***************************************************************************************
# Original Monolithic Shrinkage Estimator Implementation
# ***************************************************************************************
gwpcrpois.mom.groupwise <- function(formula, data, threshold=1, molecules=1, loss=expression(p0), ctrl=list()) {
  # Get control parameters
  ctrl.get <- function(key, default) {
    r <- ctrl[[key]]
    if (!is.null(r))
      r
    else
      default
  }
  cores <- as.numeric(ctrl.get("cores", 1))
  my.lapply <- if (cores > 1) {
    function(...) { mclapply(..., mc.cores=cores, mc.preschedule=TRUE, mc.cleanup=TRUE) }
  } else {
    lapply
  }
  verbose <- as.logical(ctrl.get("verbose", FALSE))
  var.est.distfree <- as.logical(ctrl.get("var.est.distfree", TRUE))
  obs.min.ingroup <- as.numeric(ctrl.get("obs.min.ingroup", 6))
  use.nonconv.groupest <- as.logical(ctrl.get("use.nonconv.groupest", FALSE))
  
  # Subset data to have a single "counts" columns, follows by the key columns defining
  # the groups (e.g. the gene name, but there can be more than one, e.g. sample and gene).
  # If the LHS of the formula contains an expressionlike c(col1, col2), the resulting
  # table will have more rows than "data". This allows counts in multiple columns to be
  # treated together (e.g. read counts for the plus and minus strand).
  formula.t <- terms(formula)
  group.key <- labels(formula.t)
  counts.expr <- attr(formula.t, "variables")[[1 + attr(formula.t, "response")]]
  data <- data[, list(n.umis=.N, count=eval(counts.expr)), keyby=group.key]
  if (nrow(data[!is.finite(count) | (count < threshold)]) > 0)
    stop("data contains non-finite counts or counts below the threshold")

  # Compute global parameter estimates
  if (verbose)
    message("Computing overall parameter estimates")
  mean.all = data[, mean(count)]
  var.all = data[, var(count)]
  m.all <- gwpcrpois.mom(mean.all, var.all,
                         threshold=threshold, molecules=molecules,
                         ctrl=ctrl)
  eval.loss.env <- new.env(parent=parent.frame())
  eval.loss <- function(efficiency, lambda0, p0) {
    eval.loss.env$efficiency <- efficiency
    eval.loss.env$lambda0 <- lambda0
    eval.loss.env$p0 <- p0
    tryCatch(eval(loss, envir=eval.loss.env),
             error=function(e) { warning(conditionMessage(e)) ; NA})
  }

  # Use the global estimates are starting point for further group-wise parameter estimation
  ctrl$initial <- m.all
  
  # Compute raw (group-wise) mean and variance estimates, add global estimates are columns
  if (verbose)
    message("Aggregating data per group")
  data.gen <- data[, list(mean.all=mean.all,
                          var.all=var.all,
                          efficiency.all=m.all$efficiency,
                          lambda0.all=m.all$lambda0,
                          loss.all=eval.loss(m.all$efficiency, m.all$lambda0, m.all$p0),
                          mean.raw=mean(count),
                          var.raw=var(count),
                          n.umis=n.umis[1],
                          n.obs=.N),
                   keyby=group.key]
  
  # Compute raw (group-wise) parameter estimates
  if (verbose)
    message("Computing group-wise parameter estimates on ", cores, " cores")
  data.gen.keys <- do.call(mapply, c(list(list), data.gen[, ..group.key], list(SIMPLIFY=FALSE)))
  data.gen <- rbindlist(my.lapply(data.gen.keys, function(k) {
    data.gen[k,][, c('efficiency.raw', 'lambda0.raw', 'loss.raw') := {
      tryCatch({
        if (is.finite(mean.raw) && is.finite(var.raw) && (n.obs >= obs.min.ingroup)) {
          m <- gwpcrpois.mom(mean.raw, var.raw, threshold=threshold, molecules=molecules,
                             ctrl=ctrl, nonconvergence.is.error=!use.nonconv.groupest)
          list(m$efficiency, m$lambda0, eval.loss(m$efficiency, m$lambda0, m$p0))
        } else
          list(NA, NA, NA)
      }, error=function(e) {
        cat(paste0("Failed to solve parameters for group (",
                   paste0(lapply(k, as.character), collapse=","),
                   "): ", conditionMessage(e), "\n"), file = stderr())
        list(NA, NA, NA)
      })
    }]
  }))
  setkeyv(data.gen, group.key)
  
  # Estimate variance of loss^raw (and also lambda0^raw, efficiency^raw) as a function
  # of the number of observed UMIs n.
  #
  # This is the distribution-free version of the estimator. We assume the mean loss is
  # independent from n^obs. For a fixed n^obs, the expected value of (loss^raw - m)^2
  # is then v_g + v_e / n^obs. We find v_g and v_e by mininizing the squared deviation
  # of the (single-point) variance estimate (loss^raw - m)^2 and the predictor. Since
  # variances are generally bigger for small n^obs, minimizing an unweighted sum of
  # square deviations would allow loss estimates for small n^obs to influence the
  # estimation unduely strongly. We correct for this by weigthing the squared deivations
  # with n^obs (which matches the expected drop in scale), and thus minimize the sum
  #
  #                                               2
  #    /                                        \
  #   |  (loss^raw - m)^2  -  v_g - v_e / n^obs |   * w(n^obs), where w(n) = n / (1 + n/W)
  #    \                                        /
  #
  # over all genes. m is the sample mean over all genes. We start the numerical
  # search with v_g = v_e = v/2, where v is the sampling variance of the loss.
  variance.estimates.distfree <- function(x.name) {
    # Parameter W controls the grow behaviour of weights. For small n.obs, the weights
    # equal n.obs, but as n.obs grows, weights eventually converge to W.
    W <- 100
    if (verbose)
      message("Estimating variances of ", x.name, " using the distribution-free algorithm")
    x <- parse(text=paste0(x.name, ".raw"))
    # Find least-squared estimates for v_g, v_e.
    r <- data.gen[is.finite(eval(x)) & (n.obs > 0), {
      y <- (eval(x) - mean(eval(x), na.rm=TRUE))^2
      x <- cbind(v_g=rep(1, .N), v_e=1/n.obs)
      w <- n.obs / (1 + n.obs / W)
      lsfit(x, y, wt=w, intercept=FALSE)$coefficients
    } ]
    if (all(r < 0)) {
      stop("both variance estimators are negative")
    } else if (r["v_g"] < 0) {
      warning("variance of ", x.name, " between groups estimated to be negative, setting to zero")
      # Between-gene variance estimate is negative. Assume between-gene variance
      # is negligible, and estimate only the estimator variance.
      r <- c(v_g=0, data.gen[is.finite(eval(x)) & (n.obs > 0), {
        y <- (eval(x) - mean(eval(x), na.rm=TRUE))^2
        x <- cbind(v_e=1/n.obs)
        w <- n.obs
        lsfit(x, y, wt=w, intercept=FALSE)$coefficients["v_e"]
      } ])
    } else if (r["v_e"] < 0) {
      warning("estimator variance of ", x.name, " within groups estimated to be negative, setting to zero")
      # Estimate of estimator variance is negative. Assume estimator variance is
      # negligible, and estimate only the between-gene variance.
      r <- c(v_g=data.gen[is.finite(eval(x)) & (n.obs > 0), var(eval(x), na.rm=TRUE)], v_e=0)
    }
    if (any(r < 0))
      stop("some variance estimators are negative")
    # Add columns
    data.gen[, paste0(x.name, ".grp.var") := r["v_g"]]
    data.gen[, paste0(x.name, ".raw.var") := r["v_e"] / n.obs]
  }

  # Estimate variance of p_0^raw (and also lambda0^raw, efficiency^raw) as a function
  # of the number of observed UMIs n.
  #
  # We assume that p_0^raw | n ~ Beta( a(n), b(n) ), with a fixed expectation m and
  # n-dependent variance comprising a estimation error that decreases with 1/n and a
  # constant part that reflects the variance of the loss between groups,
  #
  #        raw
  #    V( p    | n) = v(n) = v + v  /  n. 
  #        0                  g   e
  #
  # From the mean (assumed constant) and variance of the Beta distribution it follows that
  #
  #                                                      m * (1 - m)
  #   a(n) = m * f(n), b(n) = (1-m) * f(n) where f(n) = ------------- - 1.
  #                                                           v
  #                                                            e
  #                                                       v + ---
  #                                                        g   n
  #
  # For m we used the sample mean over all genes, and for v_g and v_e ML estimates
  # found using numerical optimzation. We start the numerical search with
  # v_g = v_e = v/2, where v is the sampling variance of the loss over all genes (limited
  # to the possible range (0, m * (1 - m) ).
  #
  # efficiency^raw is handled similarly. For lambda0^raw, the beta distribution
  # is replaced by a normal distribution.
  variance.estimates.beta <- function(x.name) {
    if (verbose)
      message("Estimating variances of ", x.name, " assuming a beta distribution")
    x <- parse(text=paste0(x.name, ".raw"))
    # For the mean use the sampling mean, since we assume it's independent of n
    m <- data.gen[, mean(eval(x), na.rm=TRUE) ]
    # Compute quantity to minimize, i.e. negative log-liklihood
    # We contract "x" a bit towards 0.5 to ensure that all values are strictly
    # greater than 0 and less than 1 -- otherwise, the log-likelihood becomes -Inf.
    # M is the contraction factor, i.e. we move X to the interval [M/2, 1-M/2]
    M <- 1e-6
    logl <- function(p) {
      # Reject invalid parameters. Note that the beta variance is always <= m * (1-m).
      if (any(p <= 0) || (p["v_g"] >= m * (1-m)))
        return(NA)
      # Evaluate likelihoods. We strict shape1,2 to >= 1 to ensure monomodality
      -sum(data.gen[is.finite(eval(x)) & (n.obs > 0), {
        f <- pmax(m * (1-m) / (p["v_g"] + p["v_e"] / n.obs) - 1, 1/m, 1/(1-m))
        dbeta(M/2 + eval(x)*(1-M), shape1=m*f, shape2=(1-m)*f, log=TRUE)
      }])
    }
    # Optimize v_g and v_e. Initially, we split the total observed variance evenly
    # between v_g and v_e / n.obs for the smallest positive group size n.obs.
    v <- data.gen[, min(var(eval(x), na.rm=TRUE), m * (1-m) ) ]
    n.min <- data.gen[n.obs > 0, min(n.obs) ]
    r <- optim(fn=logl, par=c(v_g=v/2, v_e=n.min*v/2), method="Nelder-Mead")
    # Correct variances for the previous contraction
    v_g <- r$par["v_g"] / (1-M)^2
    v_e <- r$par["v_e"] / (1-M)^2
    # Add columns
    data.gen[, paste0(x.name, ".raw.var") := v_e / n.obs]
    data.gen[, paste0(x.name, ".grp.var") := v_g]
  }
  variance.estimates.normal <- function(x.name) {
    if (verbose)
      message("Estimating variances of ", x.name, " assuming a normal distribution")
    x <- parse(text=paste0(x.name, ".raw"))
    m <- data.gen[, mean(eval(x), na.rm=TRUE) ]
    v <- data.gen[, var(eval(x), na.rm=TRUE) ]
    # Compute quantity to minimize, i.e. negative log-liklihood
    logl <- function(p) {
      # Reject invalid parameters
      if (any(p < 0))
        return(NA)
      # Evaluate likelihoods
      -sum(data.gen[is.finite(eval(x)) & (n.obs > 0),
        dnorm(eval(x), mean=m, sd=sqrt(p["v_g"] + p["v_e"] / n.obs), log=TRUE)
      ])
    }
    r <- optim(fn=logl, par=c(v_g=v/2, v_e=v/2), method="Nelder-Mead")
    # Add columns
    data.gen[, paste0(x.name, ".raw.var") := r$par["v_e"] / n.obs ]
    data.gen[, paste0(x.name, ".grp.var") := r$par["v_g"] ]
  }
  if (var.est.distfree) {
    variance.estimates.distfree("loss")
    variance.estimates.distfree("efficiency")
    variance.estimates.distfree("lambda0")
  } else {
    variance.estimates.beta("loss")
    variance.estimates.beta("efficiency")
    variance.estimates.normal("lambda0")
  }
  
  # Compute shrunken estimates of loss, efficiency and lambda0 i.e.
  # 
  #         all   v_e      raw
  #        p    * ---  +  p    * v_g 
  #         0      n       0
  #   p = -----------------------------------.
  #    0          v_e
  #               ---  +  v_g
  #                n
  #
  # and similarly for the other two quantities
  if (verbose)
    message("Computing final parameter estimates")
  data.gen[, loss := (loss.all * loss.raw.var + loss.raw * loss.grp.var) / (loss.raw.var + loss.grp.var) ]
  data.gen[, efficiency := (efficiency.all * efficiency.raw.var + efficiency.raw * efficiency.grp.var) / (efficiency.raw.var + efficiency.grp.var) ]
  data.gen[, lambda0 := (lambda0.all * lambda0.raw.var + lambda0.raw * lambda0.grp.var) / (lambda0.raw.var + lambda0.grp.var) ]
  # If the local estimate is NA, use the global one
  data.gen[!is.finite(loss), loss := loss.all ]
  data.gen[!is.finite(efficiency), efficiency := efficiency.all ]
  data.gen[!is.finite(lambda0), lambda0 := lambda0.all ]
  
  # Correct for unobserved UMIs, i.e. compute
  #                    n
  #   n   = -----------------------
  #    tot   (1 - p_0^)^cliquesize
  #
  # cliquesize is the number of molecules that are subjected to thresholding
  # together, i.e. ALL of them are dropped if ANY has a read count below the threshold.
  if (verbose)
    message("Computing n.tot")
  data.gen[, n.tot := ifelse(is.finite(loss),
                             n.umis / (1 - loss),
                             NA) ]

  # Return data
  return(data.gen)
}
