/*
 * File: Csetminmax.c
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 19-Nov-2015 11:10:19
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Csetminmax.h"

/* Function Definitions */

/*
 * Arguments    : double value[1302]
 *                double lowlim
 *                double uplim
 * Return Type  : void
 */
void Csetminmax(double value[1302], double lowlim, double uplim)
{
  double diff;
  int idx;
  short ii_data[1302];
  short ii_size[2];
  int i0;
  int ii;
  boolean_T exitg4;
  boolean_T guard4 = false;
  short ind_data[1302];
  double value_data[1302];
  boolean_T exitg3;
  boolean_T guard3 = false;
  boolean_T exitg2;
  boolean_T guard2 = false;
  boolean_T exitg1;
  boolean_T guard1 = false;

  /*  Function CSETMINMAX	resets value to within given limits */
  /* 	useful for ensuring that longitude is in correct range */
  /*       [value] = Csetminmax(value,lowlim,uplim) */
  /* 	 */
  /* 	input: */
  /* 		value - value to be limited */
  /* 		lolim = lower limit */
  /* 		uplimit = upper limit */
  /* 	output: */
  /* 		value - new value (between limits) */
  /*  */
  /* 	EXAMPLE: */
  /* 	a) force longitude 0<=longitude<360 */
  /* 		longitude = -180; */
  /* 		[longitude] = Csetminmax(longitude,0,360) */
  /* 	b) force longitude -pi<=longitude<pi */
  /* 		longitude = 2*pi; */
  /* 		[longitude] = Csetminmax(longitude,-pi,pi) */
  diff = uplim - lowlim;
  if (diff <= 0.0) {
  } else {
    idx = 0;
    for (i0 = 0; i0 < 2; i0++) {
      ii_size[i0] = (short)(1 + 1301 * i0);
    }

    ii = 1;
    exitg4 = false;
    while ((!exitg4) && (ii < 1303)) {
      guard4 = false;
      if (value[ii - 1] >= uplim) {
        idx++;
        ii_data[idx - 1] = (short)ii;
        if (idx >= 1302) {
          exitg4 = true;
        } else {
          guard4 = true;
        }
      } else {
        guard4 = true;
      }

      if (guard4) {
        ii++;
      }
    }

    if (1 > idx) {
      ii = 0;
    } else {
      ii = idx;
    }

    idx = ii_size[0] * ii;
    for (i0 = 0; i0 < idx; i0++) {
      ind_data[i0] = ii_data[i0];
    }

    while (!(ii == 0)) {
      for (i0 = 0; i0 < ii; i0++) {
        value_data[i0] = value[ind_data[i0] - 1] - diff;
      }

      for (i0 = 0; i0 < ii; i0++) {
        value[ind_data[i0] - 1] = value_data[i0];
      }

      idx = 0;
      for (i0 = 0; i0 < 2; i0++) {
        ii_size[i0] = (short)(1 + 1301 * i0);
      }

      ii = 1;
      exitg3 = false;
      while ((!exitg3) && (ii < 1303)) {
        guard3 = false;
        if (value[ii - 1] >= uplim) {
          idx++;
          ii_data[idx - 1] = (short)ii;
          if (idx >= 1302) {
            exitg3 = true;
          } else {
            guard3 = true;
          }
        } else {
          guard3 = true;
        }

        if (guard3) {
          ii++;
        }
      }

      if (1 > idx) {
        ii = 0;
      } else {
        ii = idx;
      }

      idx = ii_size[0] * ii;
      for (i0 = 0; i0 < idx; i0++) {
        ind_data[i0] = ii_data[i0];
      }
    }

    idx = 0;
    for (i0 = 0; i0 < 2; i0++) {
      ii_size[i0] = (short)(1 + 1301 * i0);
    }

    ii = 1;
    exitg2 = false;
    while ((!exitg2) && (ii < 1303)) {
      guard2 = false;
      if (value[ii - 1] < lowlim) {
        idx++;
        ii_data[idx - 1] = (short)ii;
        if (idx >= 1302) {
          exitg2 = true;
        } else {
          guard2 = true;
        }
      } else {
        guard2 = true;
      }

      if (guard2) {
        ii++;
      }
    }

    if (1 > idx) {
      ii = 0;
    } else {
      ii = idx;
    }

    idx = ii_size[0] * ii;
    for (i0 = 0; i0 < idx; i0++) {
      ind_data[i0] = ii_data[i0];
    }

    while (!(ii == 0)) {
      for (i0 = 0; i0 < ii; i0++) {
        value_data[i0] = value[ind_data[i0] - 1] + diff;
      }

      for (i0 = 0; i0 < ii; i0++) {
        value[ind_data[i0] - 1] = value_data[i0];
      }

      idx = 0;
      for (i0 = 0; i0 < 2; i0++) {
        ii_size[i0] = (short)(1 + 1301 * i0);
      }

      ii = 1;
      exitg1 = false;
      while ((!exitg1) && (ii < 1303)) {
        guard1 = false;
        if (value[ii - 1] < lowlim) {
          idx++;
          ii_data[idx - 1] = (short)ii;
          if (idx >= 1302) {
            exitg1 = true;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          ii++;
        }
      }

      if (1 > idx) {
        ii = 0;
      } else {
        ii = idx;
      }

      idx = ii_size[0] * ii;
      for (i0 = 0; i0 < idx; i0++) {
        ind_data[i0] = ii_data[i0];
      }
    }
  }
}

/*
 * File trailer for Csetminmax.c
 *
 * [EOF]
 */
