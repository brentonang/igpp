/*
 * File: _coder_cget_latlon_api.h
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 24-Nov-2015 01:13:59
 */

#ifndef ___CODER_CGET_LATLON_API_H__
#define ___CODER_CGET_LATLON_API_H__

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_cget_latlon_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void cget_latlon(real_T blat, real_T blon, real_T ranges, real_T azm,
  real_T major_axis, real_T esquared, real_T *plat, real_T *plon);
extern void cget_latlon_api(const mxArray * const prhs[6], const mxArray *plhs[2]);
extern void cget_latlon_atexit(void);
extern void cget_latlon_initialize(void);
extern void cget_latlon_terminate(void);
extern void cget_latlon_xil_terminate(void);

#endif

/*
 * File trailer for _coder_cget_latlon_api.h
 *
 * [EOF]
 */
