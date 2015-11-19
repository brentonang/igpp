//
// File: cadjlon.cpp
//
// MATLAB Coder version            : 2.8
// C/C++ source code generated on  : 19-Nov-2015 10:37:00
//

// Include Files
#include "rt_nonfinite.h"
#include "cadjlon.h"

// Function Definitions

//
// Arguments    : float *theta
// Return Type  : void
//
void cadjlon(float *theta)
{
  int ii_size_idx_0;
  int ii_size_idx_1;
  int ii_data[1];
  int loop_ub;
  int i0;
  int ind_data[1];
  float fv0[1];
  float fv1[1];

  // CADJLON        reduces argument to range from -pi to pi for single value,
  //  use Csetminmax instead
  //
  //        function [theta]=cadjlon(theta);
  //
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

    fv0[0] = *theta;
    i0 = 0;
    while (i0 <= 0) {
      fv0[ind_data[0] - 1] = *theta - 6.28318548F;
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

    fv1[0] = *theta;
    i0 = 0;
    while (i0 <= 0) {
      fv1[ind_data[0] - 1] = *theta + 6.28318548F;
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
  }
}

//
// File trailer for cadjlon.cpp
//
// [EOF]
//
