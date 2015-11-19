/*
 * _coder_cadjlon_mex.c
 *
 * Code generation for function '_coder_cadjlon_mex'
 *
 */

/* Include files */
#include "cadjlon.h"
#include "_coder_cadjlon_mex.h"
#include "cadjlon_terminate.h"
#include "_coder_cadjlon_api.h"
#include "cadjlon_initialize.h"
#include "cadjlon_data.h"

/* Function Declarations */
static void cadjlon_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T nrhs,
  const mxArray *prhs[1]);

/* Function Definitions */
static void cadjlon_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T nrhs,
  const mxArray *prhs[1])
{
  int32_T n;
  const mxArray *inputs[1];
  const mxArray *outputs[1];
  int32_T b_nlhs;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 1, 4, 7,
                        "cadjlon");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 7,
                        "cadjlon");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }

  /* Call the function. */
  cadjlon_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  cadjlon_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(cadjlon_atexit);

  /* Module initialization. */
  cadjlon_initialize();

  /* Dispatch the entry-point. */
  cadjlon_mexFunction(nlhs, plhs, nrhs, prhs);
}

/* End of code generation (_coder_cadjlon_mex.c) */
