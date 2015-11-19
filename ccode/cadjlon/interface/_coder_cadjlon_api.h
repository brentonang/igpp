/*
 * File: _coder_cadjlon_api.h
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 19-Nov-2015 10:37:00
 */

#ifndef ___CODER_CADJLON_API_H__
#define ___CODER_CADJLON_API_H__

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_cadjlon_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void cadjlon(real32_T *theta);
extern void cadjlon_api(const mxArray * const prhs[1], const mxArray *plhs[1]);
extern void cadjlon_atexit(void);
extern void cadjlon_initialize(void);
extern void cadjlon_terminate(void);
extern void cadjlon_xil_terminate(void);

#endif

/*
 * File trailer for _coder_cadjlon_api.h
 *
 * [EOF]
 */
