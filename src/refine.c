#include <R.h>
#include <Rinternals.h>

SEXP gwpcr_refine_c(SEXP points_, SEXP width_) {
  const int npoints = length(points_);
  if (length(points_) < 2)
    error("need at least two points");
  const double* points = REAL(points_);
  if (length(width_) != 1)
    error("width must be scalar");
  const double width = REAL(width_)[0];
  if (!R_FINITE(width) || (width <= 0))
    error("width must be positive and finite");

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
    const double n = ceil(d / width);
    const double w = d / n;
    ns[i] = n;
    ws[i] = w;
    nrefined += n;
  }

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
  refined[nrefined-1] = points[npoints-1];
  weights[0] = (refined[1] - refined[0]) * 0.5;
  if (nrefined >= 3)
    weights[nrefined-2] = (refined[nrefined-1] - refined[nrefined-3]) * 0.5;
  weights[nrefined-1] = (refined[nrefined-1] - refined[nrefined-2]) * 0.5;

  SEXP result_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result_, 0, refined_);
  SET_VECTOR_ELT(result_, 1, weights_);

  SEXP result_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(result_names_, 0, mkChar("points"));
  SET_STRING_ELT(result_names_, 1, mkChar("weights"));
  setAttrib(result_, R_NamesSymbol, result_names_);

  UNPROTECT(6);

  return result_;
}
