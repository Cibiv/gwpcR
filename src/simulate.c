#include <R.h>
#include <Rmath.h>

void gwpcr_simulate(int *nsamples_, double *samples, double *efficiency_, int *ncycles_)
{
  const int nsamples = *nsamples_;
  const double efficiency = *efficiency_;
  const int ncycles = *ncycles_;
  for(int i=0; i < ncycles; ++i)
    for(int j=0; j < nsamples; ++j)
      samples[j] += rbinom(samples[j], efficiency);
}
