/*
 * File: cadjlon.c
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 24-Nov-2015 01:15:08
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "cget_latlon.h"
#include "cadjlon.h"

/* Function Definitions */

/*
 * Arguments    : double *theta
 * Return Type  : void
 */
void cadjlon(double *theta)
{
  int ii_size_idx_0;
  int ii_size_idx_1;
  int ii_data[1];
  int loop_ub;
  int i2;
  int ind_data[1];
  double dv9[1];
  double dv10[1];

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
  i2 = 0;
  while (i2 <= loop_ub - 1) {
    ind_data[0] = 1;
    i2 = 1;
  }

  while (!((ii_size_idx_0 == 0) || (ii_size_idx_1 == 0))) {
    i2 = 0;
    while (i2 <= 0) {
      ii_data[0] = ind_data[0] - 1;
      i2 = 1;
    }

    dv9[0] = *theta;
    i2 = 0;
    while (i2 <= 0) {
      dv9[ind_data[0] - 1] = *theta - 6.2831853071795862;
      i2 = 1;
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
    i2 = 0;
    while (i2 <= loop_ub - 1) {
      ind_data[0] = ii_data[0];
      i2 = 1;
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
  i2 = 0;
  while (i2 <= loop_ub - 1) {
    ind_data[0] = ii_data[0];
    i2 = 1;
  }

  while (!((ii_size_idx_0 == 0) || (ii_size_idx_1 == 0))) {
    i2 = 0;
    while (i2 <= 0) {
      ii_data[0] = ind_data[0] - 1;
      i2 = 1;
    }

    dv10[0] = *theta;
    i2 = 0;
    while (i2 <= 0) {
      dv10[ind_data[0] - 1] = *theta + 6.2831853071795862;
      i2 = 1;
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
    i2 = 0;
    while (i2 <= loop_ub - 1) {
      ind_data[0] = ii_data[0];
      i2 = 1;
    }
  }
}

/*
 * File trailer for cadjlon.c
 *
 * [EOF]
 */
