/*
 * File: _coder_Csetminmax_api.h
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 19-Nov-2015 11:10:19
 */

#ifndef ___CODER_CSETMINMAX_API_H__
#define ___CODER_CSETMINMAX_API_H__

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_Csetminmax_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void Csetminmax(real_T value[1302], real_T lowlim, real_T uplim);
extern void Csetminmax_api(const mxArray *prhs[3], const mxArray *plhs[1]);
extern void Csetminmax_atexit(void);
extern void Csetminmax_initialize(void);
extern void Csetminmax_terminate(void);
extern void Csetminmax_xil_terminate(void);

#endif

/*
 * File trailer for _coder_Csetminmax_api.h
 *
 * [EOF]
 */
