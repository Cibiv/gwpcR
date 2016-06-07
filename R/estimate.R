gwpcrpois.mom <- function(c, molecules=1) {
  # A random variable C distributed according to the PCR-Poisson mixture
  # with parameters E (efficiency) and lambda0 has mean
  #   E( C | E, lambda0 ) = lambda0,
  # and variance
  #   V( C | E, lambda0 ) = lambda0 + lambda0^2 * V( L | E ),
  # where L is distributed according to the PCR distribution with efficiency E.
  # A method-of-moments estimator for the parameters of C is thus to set
  #   lambda0 = Avg(C),
  # and to set 
  #   v = ( Var(C) - lambda0 ) / lambda0^2
  # and find E such that
  #   V( L | E ) = v,
  # which is what gwpcr.sd.inv does using a pre-computed piecewise polynomial
  # approximation.
  #
  # NOTE: Setting a threshold, i.e. using a conditional version of the
  # Poisson distribution is not supported yet. M-o-M estimates will thus
  # be skewed unless the probability of non-observation is negigible.
  lambda0 <- mean(c)
  efficiency <- gwpcr.sd.inv(sqrt(max((var(c) - lambda0) / ( lambda0^2 ), 0)),
                             molecules=molecules)
  return(list(lambda0=lambda0, efficiency=efficiency, molecules=molecules))
}

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
  r <- optim(par=c(efficiency=mom$efficiency, lambda0=mom$lambda0),
             fn=function(p) {
               e <- p['efficiency']
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
             control=list(fnscale=-1, parscale=c(efficiency=1/mom$efficiency, lambda0=1/mom$lambda0)))
  
  # Return result
  list(efficiency=as.vector(r$par['efficiency']),
       lambda0=as.vector(r$par['lambda0']),
       threshold=threshold,
       molecules=molecules)
}
