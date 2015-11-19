/*
 * Csetminmax_initialize.c
 *
 * Code generation for function 'Csetminmax_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Csetminmax.h"
#include "Csetminmax_initialize.h"
#include "Csetminmax_data.h"

/* Function Definitions */
void Csetminmax_initialize(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (Csetminmax_initialize.c) */
