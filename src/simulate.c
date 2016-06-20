#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

static const int BATCH_MASK = (1 << 18) - 1;

void gwpcr_simulate_c(int *nsamples_, double *samples, double *efficiency_, double* molecules_, int *ncycles_)
{
  GetRNGstate();

  const int nsamples = *nsamples_;
  const double molecules = *molecules_;
  for(int j=0; j < nsamples; ++j)
    samples[j] = molecules;

  const double efficiency = *efficiency_;
  const int ncycles = *ncycles_;
  for(int i=0; i < ncycles; ++i) {
    for(int j=0; j < nsamples; ++j) {
      if (!(j & BATCH_MASK))
        R_CheckUserInterrupt();

      samples[j] += rbinom(samples[j], efficiency);
    }
  }

  const double scale = 1.0 / (molecules * pow(1+efficiency, ncycles));
  for(int j=0; j < nsamples; ++j)
    samples[j] *= scale;

  PutRNGstate();
}

void gwpcrpois_simulate_c(int *nsamples_, double *samples, double *efficiency_, double *lambda0_, double* molecules_, int *ncycles_)
{
  GetRNGstate();

  gwpcr_simulate_c(nsamples_, samples, efficiency_, molecules_, ncycles_);

  const int nsamples = *nsamples_;
  const double lambda0 = *lambda0_;
  for(int j=0; j < nsamples; ++j) {
    if (!(j & BATCH_MASK))
      R_CheckUserInterrupt();

    samples[j] = rpois(samples[j] * lambda0);
  }

  PutRNGstate();
}
