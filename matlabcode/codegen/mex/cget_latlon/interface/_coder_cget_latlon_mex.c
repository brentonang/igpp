/*
 * _coder_cget_latlon_mex.c
 *
 * Code generation for function '_coder_cget_latlon_mex'
 *
 */

/* Include files */
#include "cget_latlon.h"
#include "_coder_cget_latlon_mex.h"
#include "cget_latlon_terminate.h"
#include "_coder_cget_latlon_api.h"
#include "cget_latlon_initialize.h"
#include "cget_latlon_data.h"

/* Function Declarations */
static void cget_latlon_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T nrhs,
  const mxArray *prhs[6]);

/* Function Definitions */
static void cget_latlon_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T nrhs,
  const mxArray *prhs[6])
{
  int32_T n;
  const mxArray *inputs[6];
  const mxArray *outputs[2];
  int32_T b_nlhs;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 6, 4,
                        11, "cget_latlon");
  }

  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 11,
                        "cget_latlon");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }

  /* Call the function. */
  cget_latlon_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  cget_latlon_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(cget_latlon_atexit);

  /* Module initialization. */
  cget_latlon_initialize();

  /* Dispatch the entry-point. */
  cget_latlon_mexFunction(nlhs, plhs, nrhs, prhs);
}

/* End of code generation (_coder_cget_latlon_mex.c) */
