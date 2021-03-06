% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqsim.R
\name{seqsim}
\alias{seqsim}
\title{Simulate sequencing from a pool of DNA molecules with different abundances}
\usage{
seqsim(abundances, reads.target, efficiency, molecules = 1, cycles = Inf)
}
\arguments{
\item{abundances}{vector of (integral) abundances of different (DNA or RNA) fragments}

\item{efficiency}{efficiency of amplification}

\item{cycles}{number of amplification cycles used for simulation. By default,
a large enough value is used to make the results virtually indistinguishable
from the limit for \eqn{cycles \to \infty}}

\item{total.reads.target}{the target for the total number of sequencing reads}
}
\value{
a vector with the same length as \var{abundances} containing a read count
  for every fragment
}
\description{
Simulated amplification and sequencing of a pool of (DNA or RNA)
  fragments with the specified \var{abundances} (i.e. copy numbers). The actual
  pre-amplification abundances are multiplied with \var{molecules}, which would
  typically be set to "2" if fragments are double-stranded and both strands are
  possible templates for the DNA polymerase (although leaving \var{molecules} at
  one and doubling the \var{abundances} would have the same effect). By default,
  the PCR simulation uses an approximation that assumes many (actually, infinitly)
  many cycles, but a finite number of \var{cycles} can be specified instead.
  
  After amplification, reads are sampled from the (random) post-amplification abundance
  distribution, such that the total library size is \var{reads.target} \emph{on
  average}. The actual number of reads will fluctuate randomly around that average,
  which reflects reality more closely than always exactly hitting the targetted
  number of sequencing reads.
  
  The actual read count simulation is performed by \code{rgwpcrpois}
}
\seealso{
\code{\link{gwpcrpois}}
}
