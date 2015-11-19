/*
 * cadjlon_terminate.c
 *
 * Code generation for function 'cadjlon_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cadjlon.h"
#include "cadjlon_terminate.h"
#include "cadjlon_data.h"

/* Function Definitions */
void cadjlon_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void cadjlon_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (cadjlon_terminate.c) */
