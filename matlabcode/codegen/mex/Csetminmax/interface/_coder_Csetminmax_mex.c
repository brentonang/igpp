/*
 * _coder_Csetminmax_mex.c
 *
 * Code generation for function '_coder_Csetminmax_mex'
 *
 */

/* Include files */
#include "Csetminmax.h"
#include "_coder_Csetminmax_mex.h"
#include "Csetminmax_terminate.h"
#include "_coder_Csetminmax_api.h"
#include "Csetminmax_initialize.h"
#include "Csetminmax_data.h"

/* Function Declarations */
static void Csetminmax_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T nrhs,
  const mxArray *prhs[3]);

/* Function Definitions */
static void Csetminmax_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T nrhs,
  const mxArray *prhs[3])
{
  int32_T n;
  const mxArray *inputs[3];
  const mxArray *outputs[1];
  int32_T b_nlhs;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 3, 4,
                        10, "Csetminmax");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 10,
                        "Csetminmax");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }

  /* Call the function. */
  Csetminmax_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  Csetminmax_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(Csetminmax_atexit);

  /* Module initialization. */
  Csetminmax_initialize();

  /* Dispatch the entry-point. */
  Csetminmax_mexFunction(nlhs, plhs, nrhs, prhs);
}

/* End of code generation (_coder_Csetminmax_mex.c) */
