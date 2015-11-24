/*
 * cget_latlon_terminate.c
 *
 * Code generation for function 'cget_latlon_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cget_latlon.h"
#include "cget_latlon_terminate.h"
#include "cget_latlon_data.h"

/* Function Definitions */
void cget_latlon_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void cget_latlon_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (cget_latlon_terminate.c) */
