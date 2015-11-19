/*
 * Csetminmax.c
 *
 * Code generation for function 'Csetminmax'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Csetminmax.h"
#include "Csetminmax_data.h"

/* Variable Definitions */
static emlrtMCInfo emlrtMCI = { 24, 2, "Csetminmax",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\Csetminmax.m"
};

static emlrtRSInfo emlrtRSI = { 24, "Csetminmax",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\Csetminmax.m"
};

/* Function Declarations */
static void disp(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location);
static const mxArray *emlrt_marshallOut(const emlrtStack *sp, const char_T u[50]);

/* Function Definitions */
static void disp(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(sp, 0, NULL, 1, &pArray, "disp", true, location);
}

static const mxArray *emlrt_marshallOut(const emlrtStack *sp, const char_T u[50])
{
  const mxArray *y;
  static const int32_T iv0[2] = { 1, 50 };

  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateCharArray(2, iv0);
  emlrtInitCharArrayR2013a(sp, 50, m0, &u[0]);
  emlrtAssign(&y, m0);
  return y;
}

void Csetminmax(const emlrtStack *sp, real_T value[1302], real_T lowlim, real_T
                uplim)
{
  real_T diff;
  static const char_T cv0[50] = { 'E', 'R', 'R', 'O', 'R', ':', ' ', 'u', 'p',
    'p', 'e', 'r', ' ', 'l', 'i', 'm', 'i', 't', ' ', 'm', 'u', 's', 't', ' ',
    'b', 'e', ' ', 'h', 'i', 'g', 'h', 'e', 'r', ' ', 't', 'h', 'a', 'n', ' ',
    'l', 'o', 'w', 'e', 'r', ' ', 'l', 'i', 'm', 'i', 't' };

  int32_T idx;
  int16_T ii_data[1302];
  int16_T ii_size[2];
  int32_T i1;
  int32_T ii;
  boolean_T exitg4;
  boolean_T guard4 = false;
  int16_T ind_data[1302];
  real_T value_data[1302];
  boolean_T exitg3;
  boolean_T guard3 = false;
  boolean_T exitg2;
  boolean_T guard2 = false;
  boolean_T exitg1;
  boolean_T guard1 = false;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;

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
    st.site = &emlrtRSI;
    disp(&st, emlrt_marshallOut(&st, cv0), &emlrtMCI);
  } else {
    idx = 0;
    for (i1 = 0; i1 < 2; i1++) {
      ii_size[i1] = (int16_T)(1 + 1301 * i1);
    }

    ii = 1;
    exitg4 = false;
    while ((!exitg4) && (ii < 1303)) {
      guard4 = false;
      if (value[ii - 1] >= uplim) {
        idx++;
        ii_data[idx - 1] = (int16_T)ii;
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
    for (i1 = 0; i1 < idx; i1++) {
      ind_data[i1] = ii_data[i1];
    }

    while (!(ii == 0)) {
      for (i1 = 0; i1 < ii; i1++) {
        ii_data[i1] = ind_data[i1];
      }

      for (i1 = 0; i1 < ii; i1++) {
        value_data[i1] = value[ind_data[i1] - 1] - diff;
      }

      for (i1 = 0; i1 < ii; i1++) {
        value[ii_data[i1] - 1] = value_data[i1];
      }

      idx = 0;
      for (i1 = 0; i1 < 2; i1++) {
        ii_size[i1] = (int16_T)(1 + 1301 * i1);
      }

      ii = 1;
      exitg3 = false;
      while ((!exitg3) && (ii < 1303)) {
        guard3 = false;
        if (value[ii - 1] >= uplim) {
          idx++;
          ii_data[idx - 1] = (int16_T)ii;
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
      for (i1 = 0; i1 < idx; i1++) {
        ind_data[i1] = ii_data[i1];
      }

      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }

    idx = 0;
    for (i1 = 0; i1 < 2; i1++) {
      ii_size[i1] = (int16_T)(1 + 1301 * i1);
    }

    ii = 1;
    exitg2 = false;
    while ((!exitg2) && (ii < 1303)) {
      guard2 = false;
      if (value[ii - 1] < lowlim) {
        idx++;
        ii_data[idx - 1] = (int16_T)ii;
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
    for (i1 = 0; i1 < idx; i1++) {
      ind_data[i1] = ii_data[i1];
    }

    while (!(ii == 0)) {
      for (i1 = 0; i1 < ii; i1++) {
        ii_data[i1] = ind_data[i1];
      }

      for (i1 = 0; i1 < ii; i1++) {
        value_data[i1] = value[ind_data[i1] - 1] + diff;
      }

      for (i1 = 0; i1 < ii; i1++) {
        value[ii_data[i1] - 1] = value_data[i1];
      }

      idx = 0;
      for (i1 = 0; i1 < 2; i1++) {
        ii_size[i1] = (int16_T)(1 + 1301 * i1);
      }

      ii = 1;
      exitg1 = false;
      while ((!exitg1) && (ii < 1303)) {
        guard1 = false;
        if (value[ii - 1] < lowlim) {
          idx++;
          ii_data[idx - 1] = (int16_T)ii;
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
      for (i1 = 0; i1 < idx; i1++) {
        ind_data[i1] = ii_data[i1];
      }

      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }
  }
}

/* End of code generation (Csetminmax.c) */
