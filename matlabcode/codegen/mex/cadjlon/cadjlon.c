/*
 * cadjlon.c
 *
 * Code generation for function 'cadjlon'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cadjlon.h"
#include "cadjlon_data.h"

/* Function Definitions */
void cadjlon(const emlrtStack *sp, real32_T *theta)
{
  int32_T ii_size_idx_0;
  int32_T ii_size_idx_1;
  int32_T ii_data[1];
  int32_T loop_ub;
  int32_T i0;
  int32_T ind_data[1];
  int32_T tmp_data[1];
  real32_T fv0[1];
  real32_T fv1[1];

  /* CADJLON        reduces argument to range from -pi to pi for single value,  */
  /*  use Csetminmax instead  */
  /*  */
  /*        function [theta]=cadjlon(theta); */
  /*  */
  if (*theta > 3.1415926535897931) {
    ii_size_idx_0 = 1;
    ii_size_idx_1 = 1;
    ii_data[0] = 1;
  } else {
    ii_size_idx_0 = 0;
    ii_size_idx_1 = 0;
  }

  loop_ub = ii_size_idx_0 * ii_size_idx_1;
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    ind_data[0] = 1;
    i0 = 1;
  }

  while (!((ii_size_idx_0 == 0) || (ii_size_idx_1 == 0))) {
    i0 = 0;
    while (i0 <= 0) {
      ii_data[0] = ind_data[0] - 1;
      i0 = 1;
    }

    i0 = 0;
    while (i0 <= 0) {
      tmp_data[0] = ind_data[0];
      i0 = 1;
    }

    fv0[0] = *theta;
    i0 = 0;
    while (i0 <= 0) {
      fv0[tmp_data[0] - 1] = *theta - 6.28318548F;
      i0 = 1;
    }

    *theta = fv0[0];
    if (fv0[0] > 3.1415926535897931) {
      ii_size_idx_0 = 1;
      ii_size_idx_1 = 1;
      ii_data[0] = 1;
    } else {
      ii_size_idx_0 = 0;
      ii_size_idx_1 = 0;
    }

    loop_ub = ii_size_idx_0 * ii_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      ind_data[0] = ii_data[0];
      i0 = 1;
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  if (*theta < -3.1415926535897931) {
    ii_size_idx_0 = 1;
    ii_size_idx_1 = 1;
    ii_data[0] = 1;
  } else {
    ii_size_idx_0 = 0;
    ii_size_idx_1 = 0;
  }

  loop_ub = ii_size_idx_0 * ii_size_idx_1;
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    ind_data[0] = ii_data[0];
    i0 = 1;
  }

  while (!((ii_size_idx_0 == 0) || (ii_size_idx_1 == 0))) {
    i0 = 0;
    while (i0 <= 0) {
      ii_data[0] = ind_data[0] - 1;
      i0 = 1;
    }

    i0 = 0;
    while (i0 <= 0) {
      tmp_data[0] = ind_data[0];
      i0 = 1;
    }

    fv1[0] = *theta;
    i0 = 0;
    while (i0 <= 0) {
      fv1[tmp_data[0] - 1] = *theta + 6.28318548F;
      i0 = 1;
    }

    *theta = fv1[0];
    if (fv1[0] < -3.1415926535897931) {
      ii_size_idx_0 = 1;
      ii_size_idx_1 = 1;
      ii_data[0] = 1;
    } else {
      ii_size_idx_0 = 0;
      ii_size_idx_1 = 0;
    }

    loop_ub = ii_size_idx_0 * ii_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      ind_data[0] = ii_data[0];
      i0 = 1;
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }
}

/* End of code generation (cadjlon.c) */
