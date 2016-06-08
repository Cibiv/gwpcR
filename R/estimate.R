#' @export
gwpcrpois.mom <- function(c, threshold=1, molecules=1) {
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
  # which is what gwpcr.sd.inv does using a pre-computed piecewise polynomial
  # approximation.
  #
  # XXX: In principle, the relationship between variance and efficiency is
  # analytically tractable, so gwpcr.sd.inv could be replaced by the exact
  # formula.
  c.mean <- mean(c)
  c.var <- var(c)
  lambda0 <- c.mean
  efficiency <- gwpcr.sd.inv(sqrt(max((c.var - lambda0) / ( lambda0^2 ), 0)),
                             molecules=molecules)

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
  #   (B)  V( C | E, lambda0 ) = v = v_1 + (1 - p_0) * ( v_m + c_m^2 ) - lambda_0^2.
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
  if (threshold > 0) {
    converged <- FALSE
    rel.tol <- 1e-3
    while (!converged) {
      x <- seq(from=0, to=threshold-1, by=1)
      d <- dgwpcrpois(x, efficiency=efficiency, lambda0=lambda0,
                      threshold=0, molecules=molecules)
      p <- (1 - sum(d))
      lambda0.p <- sum(x * d)   + p * c.mean
      v         <- sum(x^2 * d) + p * (c.var + c.mean^2) - lambda0.p^2
      efficiency.p <- gwpcr.sd.inv(sqrt(max((v - lambda0.p) / ( lambda0.p^2 ), 0)),
                                   molecules=molecules)

      if ((abs(lambda0 - lambda0.p) / lambda0 <= rel.tol) &&
          (abs(efficiency - efficiency.p) / efficiency <= rel.tol))
        converged <- TRUE

      lambda0 <- lambda0.p
      efficiency <- efficiency.p
    }
  }

  # Return results
  return(list(lambda0=lambda0, efficiency=efficiency,
              threshold=threshold, molecules=molecules))
}

#' @export
gwpcrpois.mle <- function(c, threshold=1, molecules=1) {
  # Since evaluating the PCR-Poisson mixture is slow, we optimize by evaluating it only
  # once for each unique observed count, and multiplying the log-likelihood with the
  # number of times that count was observed.
  v <- rle(sort(c))

  # Use method of moments estimates as initial parameters. Since we do the parameter
  # search with clamp.efficiency set to FALSE, we must take care to clamp it to the
  # range of efficiencies found in GWPCR here
  mom <- gwpcrpois.mom(c, molecules=molecules)
  mom$efficiency <- pmin(pmax(GWPCR$efficiency[1], mom$efficiency), tail(GWPCR$efficiency,1))

  # Optimize likelihood
  r <- optim(par=c(efficiency.inv=1/mom$efficiency, lambda0=mom$lambda0),
             fn=function(p) {
               e <- 1 / p['efficiency.inv']
               l <- p['lambda0']
               # If parameters are valid, compute log-likelihood, otherwise
               # return NA.
               if ((e >= 0) && (e <= 1) && (l > 0))
                 sum(log(dgwpcrpois(c=v$values, efficiency=e, lambda0=l,
                                    threshold=threshold, molecules=molecules,
                                    clamp.efficiency=FALSE)) * v$lengths)
               else
                 as.numeric(NA)
             }, method="Nelder-Mead",
             control=list(fnscale=-1, parscale=c(efficiency.inv=mom$efficiency, lambda0=1/mom$lambda0)))

  # Return result
  list(efficiency=as.vector(1/r$par['efficiency.inv']),
       lambda0=as.vector(r$par['lambda0']),
       threshold=threshold,
       molecules=molecules)
}
