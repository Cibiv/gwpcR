# gwpcrpois.R, Copyright 2016,2017,2020c Florian G. Pflug
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

#' @title Simulate sequencing from a pool of DNA molecules with different abundances
#'
#' @description Simulated amplification and sequencing of a pool of (DNA or RNA)
#'   fragments with the specified \var{abundances} (i.e. copy numbers). The actual
#'   pre-amplification abundances are multiplied with \var{molecules}, which would
#'   typically be set to "2" if fragments are double-stranded and both strands are
#'   possible templates for the DNA polymerase (although leaving \var{molecules} at
#'   one and doubling the \var{abundances} would have the same effect). By default,
#'   the PCR simulation uses an approximation that assumes many (actually, infinitly)
#'   many cycles, but a finite number of \var{cycles} can be specified instead.
#'   
#'   After amplification, reads are sampled from the (random) post-amplification abundance
#'   distribution, such that the total library size is \var{reads.target} \emph{on
#'   average}. The actual number of reads will fluctuate randomly around that average,
#'   which reflects reality more closely than always exactly hitting the targetted
#'   number of sequencing reads.
#'   
#'   The actual read count simulation is performed by \code{rgwpcrpois}
#'
#' @param abundances vector of (integral) abundances of different (DNA or RNA) fragments
#'
#' @param reads.target the target for the total number of sequencing reads
#'
#' @param efficiency efficiency of amplification
#'
#' @param molecules initial copy number of fragments
#' 
#' @param method the method used to draw from the PCR distribution. "simulate"
#'   simulates a Galton-Watson branching process modeling PCR, "gamma" uses
#'   approximates the PCR distribution with a Gamma distribution. By default,
#'   the Gamma approximation is used for small efficiencies, where it is quite
#'   good and where simulations are computationally expensive.
#'
#' @param cycles number of amplification cycles used for simulation. By default,
#'   a large enough value is used to make the results virtually indistinguishable
#'   from the limit for \eqn{cycles \to \infty}
#'
#' @return a vector with the same length as \var{abundances} containing a read count
#'   for every fragment
#'
#' @seealso \code{\link{gwpcrpois}}
#'
#' @export
seqsim <- function(abundances, reads.target, efficiency, molecules=1, method=NULL, cycles=Inf) {
  # Validate parameters
  if (!is.numeric(abundances) || any(abundances != floor(abundances)) || any(!is.finite(abundances)) || any(abundances < 0))
    stop("abundances must be non-negative integers")
  if (!is.numeric(reads.target) || (length(reads.target) != 1) || (reads.target != floor(reads.target)) || (reads.target < 0))
    stop('reads.target must be a non-negative integral scalar')
  if (!is.numeric(efficiency) || (length(efficiency) != 1) || (efficiency < 0) || (efficiency > 1))
    stop('efficiency must be a numeric scalar within [0,1]')
  if (!is.numeric(molecules) || (length(molecules) != 1) || (molecules != floor(molecules)) || (molecules < 1))
    stop('molecules must be a positive integral scalar')
  if (!is.null(method) && (!is.character(method) || (length(method) != 1)))
    stop('method must be a single character value')
  if (!is.numeric(cycles) || (length(cycles) != 1) || (cycles != floor(cycles)) || (cycles < 0))
    stop('cycles must be a positive integral scalar or +Infinity')

  # Find the unique values appearing in the abundances and how often they appear
  a.i <- order(abundances)
  a.rle <- rle(abundances[a.i])

  # Simulate PCR and sequencing. We have to call rgwpcrpois separately for each
  # unique abundance value because it excepts a scalar for the number of molecules,
  # and then distribute the result appropriately
  r <- numeric(length(abundances))
  total.a <- sum(abundances)
  r.rle <- list(lengths=a.rle$lengths)
  r[a.i] <- unlist(lapply(1:length(a.rle$values), function(j) {
    if (a.rle$values[j] == 0) return(numeric(a.rle$lengths[j]))
    rgwpcrpois(a.rle$lengths[j], efficiency=efficiency,
               lambda0=as.numeric(a.rle$values[j])*reads.target/total.a,
               threshold=0, molecules=molecules*a.rle$values[j],
               method=method, cycles=cycles)
  }))
  return(r)
}
