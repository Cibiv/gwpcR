# estimate.R, Copyright 2016,2017 Florian G. Pflug
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

#' Method-of-Moments Parameter Estimation for PCR-Poisson Mixture
#'
#' Estimates the parameters \var{efficiency} and \var{lambda0} from the sample
#' mean and variance. Supports arbitrary detection thresholds and initial
#' molecule counts, but estimation is considerably faster in the (unrealistic)
#' case \var{threshold=0} than in the general one.
#'
#' @inheritParams gwpcrpois
#'
#' @param mean average number of observations per molecular family
#'   \emph{computed over the unambiguously detected famililies}, i.e. over those
#'   families which were observed at least \var{threshold} times.
#'
#' @param var standard deviations of number of observations per molecular
#'   family, also \emph{computed over the unambiguously detected famililies},
#'   i.e. over those families which were observed at least \var{threshold}
#'   times.
#'
#' @return A list containing the values
#'
#'   \item{convergence}{flag indicating whether the estimation converged.
#'   \code{0} indicates convergence.}
#'
#'   \item{efficiency}{parameter estimate for \var{efficiency} (see
#'   \link{gwpcrpois})}
#'
#'   \item{lambda0}{parameter estimate for \var{lambda0} (see \link{gwpcrpois})}
#'
#'   \item{pdetect}{probability of unambiguosly detecting a particular molecular
#'   family, i.e of it having at least \var{threshold} observations}
#'   
#'   \item{p0}{probability of not unambiguosly detecting a particular molecular
#'   family, i.e of it having fewer than \var{threshold} observations. Always
#'   equal to \code{1-pdetect}}
#'
#'   \item{threshold}{detection threshold specified in the call to
#'   \code{gwpcrpois.mom}}
#'
#'   \item{molecules}{initial molecule count specified in the call to
#'   \code{gwpcrpois.mom}}
#'
#' @details For the (unrealistic) uncensored case, i.e. \var{threshold=0}, the
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
#'
#'   }
#'
#' @seealso \code{\link{gwpcrpois}}
#'
#' @export
gwpcrpois.mom <- function(mean, var, threshold=1, molecules=1, ctrl=list(), nonconvergence.is.error=FALSE) {
  if (!is.numeric(mean) || (length(mean) != 1) || (mean <= 0) || (mean >= Inf))
    stop('mean must be a strictly positive numeric scalar')
  if (!is.numeric(var) || (length(var) != 1) || (var <= 0) || (var >= Inf))
    stop('var must be a strictly positive numeric scalar')
  if (!is.numeric(threshold) || (length(threshold) != 1) || (threshold != floor(threshold)) || (threshold < 0))
    stop('threshold must be a non-negative integral scalar')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a strictly positive integral scalar')

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

  if (convergence != 0) {
    if (nonconvergence.is.error)
      stop("gwpcrpois.mom did not converge")
    else
      warning("gwpcrpois.mom did not converge, returning best estimate so far")
  }

  # Return results
  return(list(lambda0=lambda0, efficiency=efficiency, pdetect=pdetect, p0=1-pdetect,
              convergence=convergence, threshold=threshold, molecules=molecules))
}

#' Maximum-Likelihood Parameter Estimation for PCR-Poisson Mixture
#'
#' @inheritParams gwpcrpois
#'
#' @export
gwpcrpois.mle <- function(c, threshold=1, molecules=1) {
  # Since evaluating the PCR-Poisson mixture is slow, we optimize by evaluating it only
  # once for each unique observed count, and multiplying the log-likelihood with the
  # number of times that count was observed.
  v <- rle(sort(c))

  # Use method of moments estimates as initial parameters. Since we do the parameter
  # search with clamp.efficiency set to FALSE, we must take care to clamp it to the
  # range of efficiencies found in GWPCR here
  mom <- gwpcrpois.mom(mean=mean(c), var=var(c), molecules=molecules, threshold=threshold)
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
  r <- optim(par=p0, fn=logl, method="Nelder-Mead",
             control=list(fnscale=-abs(logl(p0) - logl(p0 * c(0.9, 0.9))),
                          parscale=c(efficiency=p0['efficiency']/10,
                                     lambda0=p0['lambda0']/10)))

  # Return result
  # XXX: Also return p0 and pdetect
  list(convergence=r$convergence,
       efficiency=as.vector(r$par['efficiency']),
       lambda0=as.vector(r$par['lambda0']),
       threshold=threshold,
       molecules=molecules)
}

#' Shrinkage Estimator for Individual Parameters of PCR-Poisson Mixture
#' 
#' Low-noise estimates of \var{efficiency} and \var{lambda0} for subsets of the
#' data (e.g.' for individual genes of an RNA-Seq experiment), even if those
#' subsets contain only few UMIs.
#' 
#' XXX
#'
#' @inheritParams gwpcrpois
#'
#' @export
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
                             ctrl=ctrl, nonconvergence.is.error=TRUE)
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
