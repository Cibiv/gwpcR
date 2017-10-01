/*
 refine.c, Copyright 2016,2017 Florian G. Pflug

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
#include <R.h>
#include <Rinternals.h>

SEXP gwpcr_refine_c(SEXP points_, SEXP width_) {
  /* Check parameters */
  const int npoints = length(points_);
  if (npoints < 2)
    error("need at least two points");
  const double* points = REAL(points_);
  const int nwidths = length(width_);
  if (nwidths < 1)
    error("need at least one width");
  const double* widths = REAL(width_);

  /* Compute for each pair p_i, p_{i+1} the number of intermediate points
   * n_i and their distance w_i, and store these in vectors ns and ws.
   */
  SEXP ns_ = PROTECT(allocVector(REALSXP, npoints-1));
  double* const ns = REAL(ns_);
  SEXP ws_ = PROTECT(allocVector(REALSXP, npoints-1));
  double* const ws = REAL(ws_);

  int nrefined = 1;
  for(int i=0; i < (npoints-1); ++i) {
    const double d = points[i+1] - points[i];
    if (d <= 0) {
      UNPROTECT(2);
      error("points must be strictly increasing");
    }

    const double width = widths[i % nwidths];
    if (width == R_PosInf) {
      ns[i] = 1;
      ws[i] = d;
    }
    else {
      if (!(width > 0))
        error("widths must be positive");

      const double n = ceil(d / width);
      const double w = d / n;
      ns[i] = n;
      ws[i] = w;
    }

    nrefined += ns[i];
  }

  /* Compute refined grid (refined), and the associated weights, i.e.
   * the lengths of centered intervals around each point in the refined
   * grid. At the edges, the intervals are one-sided instead of centered.
   */
  SEXP refined_ = PROTECT(allocVector(REALSXP, nrefined));
  double* const refined = REAL(refined_);
  SEXP weights_ = PROTECT(allocVector(REALSXP, nrefined));
  double* const weights = REAL(weights_);
  int j = 0;
  for(int i=0; i < (npoints-1); ++i) {
    const int n = ns[i];
    for(int k = 0; k < n; ++k, ++j) {
      refined[j] = points[i] + k * ws[i];
      if (j >= 2)
        weights[j-1] = (refined[j] - refined[j-2]) * 0.5;
    }
  }
  /* The loop above never reaches the right edge point */
  refined[nrefined-1] = points[npoints-1];
  /* and thus also fails to compute the weight of the last intermediate point */
  if (nrefined >= 3)
    weights[nrefined-2] = (refined[nrefined-1] - refined[nrefined-3]) * 0.5;
  /* The weights of the edge points are computed differently, since the
   * intervals are one-sided there.
   */
  weights[0] = (refined[1] - refined[0]) * 0.5;
  weights[nrefined-1] = (refined[nrefined-1] - refined[nrefined-2]) * 0.5;

  /* Create result list and add refined points and weights */
  SEXP result_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result_, 0, refined_);
  SET_VECTOR_ELT(result_, 1, weights_);

  /* Set names of result list elements */
  SEXP result_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(result_names_, 0, mkChar("points"));
  SET_STRING_ELT(result_names_, 1, mkChar("weights"));
  setAttrib(result_, R_NamesSymbol, result_names_);

  /* Result (1), its names (1), its elements (2), and temporaries ws and ns (2) */
  UNPROTECT(6);

  return result_;
}
