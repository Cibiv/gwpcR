/*
 simulate.c, Copyright 2016,2017 Florian G. Pflug

 This file is part of gwpcR

 Foobar is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

 Foobar is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

static const int BATCH_MASK = (1 << 18) - 1;

void gwpcr_simulate_c(int *nsamples_, double *samples, double *samples_tmp, double *efficiency_, double* molecules_, int *mincycles_, int *maxcycles_)
{
  GetRNGstate();

  const int nsamples = *nsamples_;
  const double molecules = *molecules_;
  for(int j=0; j < nsamples; ++j)
    samples[j] = molecules;

  /* Simulate for <mincycles> cycles */
  const double efficiency = *efficiency_;
  const int mincycles = *mincycles_;
  int cycles=0;
  for(; cycles < mincycles; ++cycles) {
    /* Simulate cycle */
    for(int j=0; j < nsamples; ++j) {
      if (!(j & BATCH_MASK))
        R_CheckUserInterrupt();

      samples[j] += rbinom(samples[j], efficiency);
    }
  }

  /* Simulate until no ties remain, but at most until maxcycles are reached */
  const int maxcycles = *maxcycles_;
  for(; cycles <= maxcycles; ++cycles) {
    /* Stop if there are no ties */
    int ties=0;
    memcpy(samples_tmp, samples, nsamples * sizeof(*samples_tmp));
    R_rsort(samples_tmp, nsamples);
    for(int j=1; j < nsamples; ++j)
      if (samples_tmp[j-1] == samples_tmp[j]) {
        ties=1;
        break;
      }
    if (!ties)
      break;

    /* Simulate cycle */
    for(int j=0; j < nsamples; ++j) {
      if (!(j & BATCH_MASK))
        R_CheckUserInterrupt();

      samples[j] += rbinom(samples[j], efficiency);
    }
  }

  /* Scale to expectation 1.0 */
  const double scale = 1.0 / (molecules * pow(1+efficiency, cycles));
  for(int j=0; j < nsamples; ++j)
    samples[j] *= scale;

  PutRNGstate();
}

void gwpcrpois_simulate_c(int *nsamples_, double *samples, double *samples_tmp, double *efficiency_, double *lambda0_, double* molecules_, int *mincycles_, int *maxcycles_)
{
  GetRNGstate();

  gwpcr_simulate_c(nsamples_, samples, samples_tmp, efficiency_, molecules_, mincycles_, maxcycles_);

  const int nsamples = *nsamples_;
  const double lambda0 = *lambda0_;
  for(int j=0; j < nsamples; ++j) {
    if (!(j & BATCH_MASK))
      R_CheckUserInterrupt();

    samples[j] = rpois(samples[j] * lambda0);
  }

  PutRNGstate();
}
