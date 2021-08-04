#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void gwpcr_simulate_c(int *nsamples_, double *samples, double *samples_tmp, double *efficiency_, double* molecules_, int *mincycles_, int *maxcycles_);
extern void gwpcrpois_simulate_c(int *nsamples_, double *samples, double *samples_tmp, double *efficiency_, double *lambda0_, double* molecules_, int *mincycles_, int *maxcycles_);

/* .Call calls */
extern SEXP gwpcr_refine_c(SEXP points_, SEXP width_);

static const R_CMethodDef CEntries[] = {
    {"gwpcr_simulate_c",     (DL_FUNC) &gwpcr_simulate_c,     7},
    {"gwpcrpois_simulate_c", (DL_FUNC) &gwpcrpois_simulate_c, 8},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"gwpcr_refine_c", (DL_FUNC) &gwpcr_refine_c, 2},
    {NULL, NULL, 0}
};

void R_init_gwpcR(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
