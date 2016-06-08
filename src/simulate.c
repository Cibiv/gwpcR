#include <R.h>
#include <Rmath.h>

void gwpcr_simulate(int *nsamples_, double *samples, double *efficiency_, double* molecules_, int *ncycles_)
{
  GetRNGstate();

  const int nsamples = *nsamples_;
  const double molecules = *molecules_;
  for(int j=0; j < nsamples; ++j)
    samples[j] = molecules;

  const double efficiency = *efficiency_;
  const int ncycles = *ncycles_;
  for(int i=0; i < ncycles; ++i)
    for(int j=0; j < nsamples; ++j)
      samples[j] += rbinom(samples[j], efficiency);

  const double scale = 1.0 / (molecules * pow(1+efficiency, ncycles));
  for(int j=0; j < nsamples; ++j)
    samples[j] *= scale;

  PutRNGstate();
}

void gwpcrpois_simulate(int *nsamples_, double *samples, double *efficiency_, double *lambda0_, double* molecules_, int *ncycles_)
{
  GetRNGstate();

  gwpcr_simulate(nsamples_, samples, efficiency_, molecules_, ncycles_);

  const int nsamples = *nsamples_;
  const double lambda0 = *lambda0_;
  for(int j=0; j < nsamples; ++j)
    samples[j] = rpois(samples[j] * lambda0);

  PutRNGstate();
}
