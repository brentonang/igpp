/*
 * cget_latlon.h
 *
 * Code generation for function 'cget_latlon'
 *
 */

#ifndef __CGET_LATLON_H__
#define __CGET_LATLON_H__

/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "cget_latlon_types.h"

/* Function Declarations */
extern void cget_latlon(const emlrtStack *sp, real_T blat, real_T blon, real_T
  ranges, real_T azm, real_T major_axis, real_T esquared, real_T *plat, real_T
  *plon);

#endif

/* End of code generation (cget_latlon.h) */
