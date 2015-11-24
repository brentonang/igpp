/*
 * cadjlon.c
 *
 * Code generation for function 'cadjlon'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cget_latlon.h"
#include "cadjlon.h"
#include "cget_latlon_data.h"

/* Function Definitions */
void cadjlon(const emlrtStack *sp, real_T *theta)
{
  int32_T ii_size_idx_0;
  int32_T ii_size_idx_1;
  int32_T ii_data[1];
  int32_T loop_ub;
  int32_T i1;
  int32_T ind_data[1];
  int32_T tmp_data[1];
  real_T dv9[1];
  real_T dv10[1];

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
  i1 = 0;
  while (i1 <= loop_ub - 1) {
    ind_data[0] = 1;
    i1 = 1;
  }

  while (!((ii_size_idx_0 == 0) || (ii_size_idx_1 == 0))) {
    i1 = 0;
    while (i1 <= 0) {
      ii_data[0] = ind_data[0] - 1;
      i1 = 1;
    }

    i1 = 0;
    while (i1 <= 0) {
      tmp_data[0] = ind_data[0];
      i1 = 1;
    }

    dv9[0] = *theta;
    i1 = 0;
    while (i1 <= 0) {
      dv9[tmp_data[0] - 1] = *theta - 6.2831853071795862;
      i1 = 1;
    }

    *theta = dv9[0];
    if (dv9[0] > 3.1415926535897931) {
      ii_size_idx_0 = 1;
      ii_size_idx_1 = 1;
      ii_data[0] = 1;
    } else {
      ii_size_idx_0 = 0;
      ii_size_idx_1 = 0;
    }

    loop_ub = ii_size_idx_0 * ii_size_idx_1;
    i1 = 0;
    while (i1 <= loop_ub - 1) {
      ind_data[0] = ii_data[0];
      i1 = 1;
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
  i1 = 0;
  while (i1 <= loop_ub - 1) {
    ind_data[0] = ii_data[0];
    i1 = 1;
  }

  while (!((ii_size_idx_0 == 0) || (ii_size_idx_1 == 0))) {
    i1 = 0;
    while (i1 <= 0) {
      ii_data[0] = ind_data[0] - 1;
      i1 = 1;
    }

    i1 = 0;
    while (i1 <= 0) {
      tmp_data[0] = ind_data[0];
      i1 = 1;
    }

    dv10[0] = *theta;
    i1 = 0;
    while (i1 <= 0) {
      dv10[tmp_data[0] - 1] = *theta + 6.2831853071795862;
      i1 = 1;
    }

    *theta = dv10[0];
    if (dv10[0] < -3.1415926535897931) {
      ii_size_idx_0 = 1;
      ii_size_idx_1 = 1;
      ii_data[0] = 1;
    } else {
      ii_size_idx_0 = 0;
      ii_size_idx_1 = 0;
    }

    loop_ub = ii_size_idx_0 * ii_size_idx_1;
    i1 = 0;
    while (i1 <= loop_ub - 1) {
      ind_data[0] = ii_data[0];
      i1 = 1;
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }
}

/* End of code generation (cadjlon.c) */
