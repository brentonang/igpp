/*
 * rdivide.h
 *
 * Code generation for function 'rdivide'
 *
 */

#ifndef __RDIVIDE_H__
#define __RDIVIDE_H__

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
extern void b_rdivide(const emlrtStack *sp, const real_T x_data[], const int32_T
                      x_size[2], const real_T y_data[], const int32_T y_size[2],
                      real_T z_data[], int32_T z_size[2]);
extern real_T rdivide(real_T x, real_T y);

#ifdef __WATCOMC__

#pragma aux rdivide value [8087];

#endif
#endif

/* End of code generation (rdivide.h) */
