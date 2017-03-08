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
gwpcrpois.mom <- function(mean, var, threshold=1, molecules=1, ctrl=list()) {
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

  if (convergence != 0)
    warning("gwpcrpois.mom did not converge, returning best estimate so far")

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

#' Shrinkage Estimator for Individual Parameters of PCR-Poisson Mixture (OLD VERSION)
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
gwpcrpois.mom.groupwise.old <- function(formula, data, threshold=1, molecules=1, clique.size=1, ctrl=list()) {
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
  shrink.only <- as.logical(ctrl.get("shrink.only", FALSE))
  include.unshrunken.estimates <- as.logical(ctrl.get("include.unshrunken.estimates", FALSE))
  
  # Subset data to have the "count" column (left side of formula) as the first column,
  # followed by the key columns defining the groups (e.g. the gene name, but can be more
  # than one).
  data.f <- data.table(model.frame(formula, data))
  group.key <- labels(terms(formula))
  colnames(data.f)[1] <- "count"
  data.th <- data.f[count >= threshold, ]
  setkeyv(data.th, group.key)
  
  # Global mean and variance
  mean.global = data.th[, mean(count)]
  var.global = data.th[, var(count)]
  
  # Use method of moments to find global parameter estimates, and use those
  # as initial values for further parameter searches
  if (!shrink.only) {
    m.global <- gwpcrpois.mom(mean.global, var.global,
                              threshold=threshold, molecules=molecules,
                              ctrl=ctrl)
    ctrl$initial <- m.global
  }
  
  # Groupwise means and variances
  data.gen <- data.th[data.f[, list(dummy=1), by=group.key]
                      , c(list(mean.global=mean.global,
                               var.global=var.global),
                          if (!shrink.only)
                            list(efficiency.global=m.global$efficiency,
                                 lambda0.lgobal=m.global$lambda0,
                                 p0.global=m.global$p0)
                          else list(),
                          list(mean.withingroup=mean(count),
                               var.withingroup=var(count),
                               n=.N))
                      , on=group.key, by=.EACHI]
  setkeyv(data.gen, group.key)
  
  if (include.unshrunken.estimates) {
    # Estimates using raw (unshrunken) groupwise mean and variance
    data.gen.keys <- do.call(mapply, c(list(list), data.gen[, ..group.key], list(SIMPLIFY=FALSE)))
    data.gen <- rbindlist(my.lapply(data.gen.keys, function(k) {
      data.gen[k,][, c('efficiency.unshrunken', 'lambda0.unshrunken', 'p0.unshrunken') := {
        m <- tryCatch({
          gwpcrpois.mom(mean.withingroup, var.withingroup, threshold=threshold, molecules=molecules, ctrl=ctrl)
        }, error=function(e) {
          cat(paste0("Failed to solve for unshrunken parameters for group (",
                     paste0(lapply(k, as.character), collapse=","),
                     "): ", conditionMessage(e), "\n"), file = stderr())
          list(efficiency=NA, lambda0=NA, p0=NA)
        })
        list(m$efficiency, m$lambda0, m$p0)
      }]
    }))
    setkeyv(data.gen, group.key)
    
    # Correct for unobserved UMIs
    data.gen[, n.tot.unshrunken := ifelse(is.finite(p0.unshrunken),
                                          n / ((1 - p0.unshrunken)^clique.size),
                                          NA) ]
  }
  
  # Between-groups variance of groupwise mean/variance XXX 1:500
  mean.intergroup.var <- data.gen[order(n, decreasing=TRUE), var(mean.withingroup, na.rm=TRUE)]
  var.intergroup.var <- data.gen[order(n, decreasing=TRUE), var(var.withingroup, na.rm=TRUE)]
  
  # Shrink genewise means towards global mean according to estimated standard error
  # of the genewise mean. We assume that the variance is the same between groups when
  # computing the standard error
  data.gen[, mean.intergroup.var := mean.intergroup.var ]
  data.gen[, mean.withingroup.stderr := var.global / n ]
  data.gen[, mean := ifelse(is.finite(mean.withingroup.stderr),
                            (mean.withingroup * mean.intergroup.var + mean.global * mean.withingroup.stderr) /
                              (mean.intergroup.var + mean.withingroup.stderr),
                            mean.global) ]
  
  # Shrink genewise variances towards global variance according to estimated standard error
  # of the genwise variance. We assume that the 4th moment (m4.glob) is the same between
  # groups when computing the standard error.
  m4.glob <- data.th[, mean((count - mean.global)^4) ]
  data.gen[, var.intergroup.var := var.intergroup.var ]
  data.gen[, var.withingroup.stderr := 
             ifelse(n >= 3,
                    pmax(m4.glob - (var.global ^ 2) * (n-3)/(n-1), 0) / n,
                    Inf) ]
  data.gen[, var :=
             ifelse(is.finite(var.withingroup.stderr),
                    (var.withingroup * var.intergroup.var + var.global * var.withingroup.stderr) /
                      (var.intergroup.var + var.withingroup.stderr),
                    var.global) ]
  
  # If we're only supposed to compute shrunken estimates of mean and variance,
  # we're done here.
  if (shrink.only)
    return(data.gen)
  
  # Use method of moments to find model parameters.
  data.gen.keys <- do.call(mapply, c(list(list), data.gen[, ..group.key], list(SIMPLIFY=FALSE)))
  data.gen <- rbindlist(my.lapply(data.gen.keys, function(k) {
    data.gen[k,][, c('efficiency', 'lambda0', 'p0') := {
      m <- tryCatch({
        gwpcrpois.mom(mean, var, threshold=threshold, molecules=molecules, ctrl=ctrl)
      }, error=function(e) {
        cat(paste0("Failed to solve for parameters for group (",
                   paste0(lapply(k, as.character), collapse=","),
                   "): ", conditionMessage(e), "\n"), file = stderr())
        list(efficiency=NA, lambda0=NA, p0=NA)
      })
      list(m$efficiency, m$lambda0, m$p0)
    }]
  }))
  setkeyv(data.gen, group.key)
  
  # Correct for unobserved UMIs
  data.gen[, n.tot := ifelse(is.finite(p0),
                             n / ((1 - p0)^clique.size),
                             NA) ]
  
  # Return data
  return(data.gen)
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
gwpcrpois.mom.groupwise <- function(formula, data, threshold=1, molecules=1, clique.size=1, ctrl=list()) {
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
  
  # Subset data to have the "count" column (left side of formula) as the first column,
  # followed by the key columns defining the groups (e.g. the gene name, but can be more
  # than one).
  data.f <- data.table(model.frame(formula, data))
  group.key <- labels(terms(formula))
  colnames(data.f)[1] <- "count"
  data.th <- data.f[count >= threshold, ]
  setkeyv(data.th, group.key)
  
  # Compute global parameter estimates
  if (verbose)
    message("Computing overall parameter estimates")
  mean.all = data.th[, mean(count)]
  var.all = data.th[, var(count)]
  m.all <- gwpcrpois.mom(mean.all, var.all,
                            threshold=threshold, molecules=molecules,
                            ctrl=ctrl)
  
  # Use the global estimates are starting point for further group-wise parameter estimation
  ctrl$initial <- m.all
  
  # Compute raw (group-wise) mean and variance estimates, add global estimates are columns
  if (verbose)
    message("Aggregating data per group")
  data.gen <- data.th[data.f[, list(dummy=1), by=group.key]
                      , list(mean.all=mean.all,
                             var.all=var.all,
                             efficiency.all=m.all$efficiency,
                             lambda0.all=m.all$lambda0,
                             p0.all=m.all$p0,
                             mean.raw=mean(count),
                             var.raw=var(count),
                             n=.N)
                      , on=group.key, by=.EACHI]
  setkeyv(data.gen, group.key)
  
  # Compute raw (group-wise) parameter estimates
  if (verbose)
    message("Computing group-wise parameter estimates on ", cores, " cores")
  data.gen.keys <- do.call(mapply, c(list(list), data.gen[, ..group.key], list(SIMPLIFY=FALSE)))
  data.gen <- rbindlist(my.lapply(data.gen.keys, function(k) {
    data.gen[k,][, c('efficiency.raw', 'lambda0.raw', 'p0.raw') := {
      m <- tryCatch({
        if (is.finite(mean.raw) && is.finite(var.raw))
          gwpcrpois.mom(mean.raw, var.raw, threshold=threshold, molecules=molecules, ctrl=ctrl)
        else
          list(efficiency=NA, lambda0=NA, p0=NA)
      }, error=function(e) {
        cat(paste0("Failed to solve for unshrunken parameters for group (",
                   paste0(lapply(k, as.character), collapse=","),
                   "): ", conditionMessage(e), "\n"), file = stderr())
        list(efficiency=NA, lambda0=NA, p0=NA)
      })
      list(m$efficiency, m$lambda0, m$p0)
    }]
  }))
  setkeyv(data.gen, group.key)
  
  # Estimate variance of p_0^raw as a function of the number of observed UMIs n.
  # We assume that p_0^raw | n ~ Beta( a(n), b(n) ), with a fixed expectation m and
  # n-dependent variance comprising a estimation error that decreases with 1/n and a
  # constant part that reflects the variance of p0 between groups,
  #
  #        raw
  #    V( p    | n) = v(n) = v + v  /  n. 
  #        0                  g   e
  #
  # From the mean and variance of the Beta distribution it follows that
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
  # v_g = v_e = v/2, where v is the sampling variance of p0 over all genes (limited
  # to the possible range (0, m * (1 - m) ).
  if (verbose)
    message("Estimating variances of p0.raw")
  m <- data.gen[, mean(p0.raw, na.rm=TRUE) ]
  v <- data.gen[, min(var(p0.raw, na.rm=TRUE), m * (1-m) * 0.9)]
  logl <- function(p) {
    # Reject invalid parameters ( minimal n = 1, hence v_g + v_e < m * (1-m) ).
    if (any(p <=0) || (sum(p) >= m * (1-m)))
      return(NA)
    # Compute quantity to minimize, i.e. negative log-liklihood
    data.gen[is.finite(p0.raw), {
      f <- m * (1-m) / (p["v_g"] + p["v_e"] / n) - 1
      -sum(dbeta(p0.raw, shape1=m*f, shape2=(1-m)*f, log=TRUE))
    } ]
  }
  r <- optim(fn=logl, par=c(v_g=v/2, v_e=v/2), method="Nelder-Mead")
  v_g <- r$par["v_g"]
  v_e <- r$par["v_e"]
  data.gen[, p0.raw.var := v_e / n]
  data.gen[, p0.grp.var := v_g]
  
  # Compute shrunken estimates of p0, i.e.
  # 
  #         all   v_e      raw
  #        p    * ---  +  p    * v_g 
  #         0      n       0
  #   p = -----------------------------------.
  #    0          v_e
  #               ---  +  v_g
  #                n
  #
  if (verbose)
    message("Computing final (shrinkage) estimates for p0")
  data.gen[, p0 := (p0.all * p0.raw.var + p0.raw * p0.grp.var) / (p0.raw.var + p0.grp.var) ]
  
  # Correct for unobserved UMIs, i.e. compute
  #                    n
  #   n   = -----------------------
  #    tot   (1 - p_0^)^cliquesize
  #
  # cliquesize is the number of molecules that are subjected to thresholding
  # together, i.e. ALL of them are dropped if ANY has a read count below the threshold.
  if (verbose)
    message("Computing n.tot")
  data.gen[, n.tot := ifelse(is.finite(p0),
                             n / ((1 - p0)^clique.size),
                             NA) ]

  # Return data
  return(data.gen)
}
