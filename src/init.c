#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern void F77_NAME(pp_etas_trigger)(int *tflag, int *n,
                                      double *mu, double *k, double *c,
                                      double *p, double *g, double *d,
                                      double *q,
                                      double *x, double *y, double *t,
                                      double *m, double *predictor,
                                      double *l);

static const R_FortranMethodDef FortranEntries[] = {
  {"pp_etas_trigger", (DL_FUNC) &F77_NAME(pp_etas_trigger), 15},
  {NULL, NULL, 0}
};

void R_init_ppnexus(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
