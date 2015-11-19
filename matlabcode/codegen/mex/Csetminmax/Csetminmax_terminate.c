/*
 * Csetminmax_terminate.c
 *
 * Code generation for function 'Csetminmax_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Csetminmax.h"
#include "Csetminmax_terminate.h"
#include "Csetminmax_data.h"

/* Function Definitions */
void Csetminmax_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void Csetminmax_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (Csetminmax_terminate.c) */
