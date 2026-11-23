#include "simulation_workflow.h"
#include <stdexcept>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
extern "C" SEXP C_gsim_simulation_workflow(SEXP wire,SEXP path){
 try {if(TYPEOF(wire)!=STRSXP||XLENGTH(wire)!=1||STRING_ELT(wire,0)==NA_STRING||TYPEOF(path)!=STRSXP||XLENGTH(path)!=1||STRING_ELT(path,0)==NA_STRING)throw std::runtime_error("recipe/output must be nonmissing strings");
 gsim::native::simulate_recipe(Rf_translateCharUTF8(STRING_ELT(wire,0)),Rf_translateCharUTF8(STRING_ELT(path,0)));return Rf_ScalarLogical(1);
 }catch(const std::exception&e){Rf_error("simulation workflow: %s",e.what());}return R_NilValue;}
