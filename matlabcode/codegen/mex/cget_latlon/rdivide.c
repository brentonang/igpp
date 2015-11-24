/*
 * rdivide.c
 *
 * Code generation for function 'rdivide'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cget_latlon.h"
#include "rdivide.h"

/* Variable Definitions */
static emlrtRTEInfo b_emlrtRTEI = { 13, 15, "rdivide",
  "C:\\Program Files\\MATLAB\\MATLAB Production Server\\R2015a\\toolbox\\eml\\lib\\matlab\\ops\\rdivide.m"
};

/* Function Definitions */
void b_rdivide(const emlrtStack *sp, const real_T x_data[], const int32_T
               x_size[2], const real_T y_data[], const int32_T y_size[2], real_T
               z_data[], int32_T z_size[2])
{
  int8_T varargin_1[2];
  int8_T varargin_2[2];
  int32_T k;
  boolean_T p;
  boolean_T b_p;
  boolean_T exitg1;
  int32_T loop_ub;
  for (k = 0; k < 2; k++) {
    varargin_1[k] = (int8_T)x_size[k];
    varargin_2[k] = (int8_T)y_size[k];
  }

  p = false;
  b_p = true;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 2)) {
    if (!(varargin_1[k] == varargin_2[k])) {
      b_p = false;
      exitg1 = true;
    } else {
      k++;
    }
  }

  if (!b_p) {
  } else {
    p = true;
  }

  if (p) {
  } else {
    emlrtErrorWithMessageIdR2012b(sp, &b_emlrtRTEI, "MATLAB:dimagree", 0);
  }

  z_size[0] = x_size[0];
  z_size[1] = x_size[1];
  loop_ub = x_size[0] * x_size[1];
  for (k = 0; k < loop_ub; k++) {
    z_data[k] = x_data[k] / y_data[k];
  }
}

real_T rdivide(real_T x, real_T y)
{
  return x / y;
}

/* End of code generation (rdivide.c) */
