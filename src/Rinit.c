// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


void freelsgen(int* sampsize, double* xdata, double* ydata, int* ord,
      int* numknot, int* initseed, int* initstream, double* optknot, 
           double* trace_hat, double* GCV, double* GSJS_value);


static R_NativePrimitiveArgType freelsgen_t[] = {
    INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

void freelsgold(int* sampsize, double* xdata, double* ydata, int* ord,
      int* numknot, int* initseed, int* initstream, double* optknot, 
      double* trace_hat, double* GCV, double* GSJS_value);

static R_NativePrimitiveArgType freelsgold_t[] = {
    INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP
};


void freepsgen(int* sampsize, double* xdata, double* ydata, int* ord,
      int* numknot, int* initseed, int* initstream, double* lambda,
      double* optknot, double* trace_hat, double* GCV, double* GSJS_value);

static R_NativePrimitiveArgType freepsgen_t[] = {
    INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

void freepsgold(int* sampsize, double* xdata, double* ydata, int* ord,
      int* numknot, int* initseed, int* initstream, double* lambda,
      double* optknot, double* trace_hat, double* GCV, double* GSJS_value);

static R_NativePrimitiveArgType freepsgold_t[] = {
    INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static const R_CMethodDef cMethods[] = {
   {"freelsgen", (DL_FUNC) &freelsgen, 11, freelsgen_t},
   {"freelsgold", (DL_FUNC) &freelsgold, 11, freelsgold_t},
   {"freepsgen", (DL_FUNC) &freepsgen, 12, freepsgen_t},
   {"freepsgold", (DL_FUNC) &freepsgold, 12, freepsgold_t},
   {NULL, NULL, 0, NULL}
};

void R_init_freeknotsplines(DllInfo* info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
