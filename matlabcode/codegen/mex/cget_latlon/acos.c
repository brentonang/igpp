/*
 * acos.c
 *
 * Code generation for function 'acos'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cget_latlon.h"
#include "acos.h"
#include "eml_error.h"

/* Variable Definitions */
static emlrtRSInfo j_emlrtRSI = { 15, "acos",
  "C:\\Program Files\\MATLAB\\MATLAB Production Server\\R2015a\\toolbox\\eml\\lib\\matlab\\elfun\\acos.m"
};

/* Function Definitions */
void b_acos(const emlrtStack *sp, real_T *x)
{
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  if ((*x < -1.0) || (1.0 < *x)) {
    st.site = &j_emlrtRSI;
    b_eml_error(&st);
  }

  *x = muDoubleScalarAcos(*x);
}

/* End of code generation (acos.c) */
