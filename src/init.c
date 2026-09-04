#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP C_gsim_hapnest_founders(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                     SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_hapnest_segment_endpoint(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_hapnest_copy_segment(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_hapnest_pair(SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
    {"C_gsim_hapnest_founders", (DL_FUNC) &C_gsim_hapnest_founders, 14},
    {"C_gsim_hapnest_segment_endpoint", (DL_FUNC) &C_gsim_hapnest_segment_endpoint, 4},
    {"C_gsim_hapnest_copy_segment", (DL_FUNC) &C_gsim_hapnest_copy_segment, 5},
    {"C_gsim_hapnest_pair", (DL_FUNC) &C_gsim_hapnest_pair, 2},
    {NULL, NULL, 0}
};

void attribute_visible R_init_gsim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
