/*
 * _coder_Csetminmax_api.c
 *
 * Code generation for function '_coder_Csetminmax_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "Csetminmax.h"
#include "_coder_Csetminmax_api.h"
#include "Csetminmax_data.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[1302];
static void b_emlrt_marshallOut(const real_T u[1302], const mxArray *y);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *lowlim,
  const char_T *identifier);
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[1302];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *value,
  const char_T *identifier))[1302];
static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[1302]
{
  real_T (*y)[1302];
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static void b_emlrt_marshallOut(const real_T u[1302], const mxArray *y)
{
  static const int32_T iv1[2] = { 1, 1302 };

  mxSetData((mxArray *)y, (void *)u);
  emlrtSetDimensions((mxArray *)y, iv1, 2);
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *lowlim,
  const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = d_emlrt_marshallIn(sp, emlrtAlias(lowlim), &thisId);
  emlrtDestroyArray(&lowlim);
  return y;
}

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = f_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[1302]
{
  real_T (*ret)[1302];
  int32_T iv2[2];
  int32_T i0;
  for (i0 = 0; i0 < 2; i0++) {
    iv2[i0] = 1 + 1301 * i0;
  }

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, iv2);
  ret = (real_T (*)[1302])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *value,
  const char_T *identifier))[1302]
{
  real_T (*y)[1302];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(sp, emlrtAlias(value), &thisId);
  emlrtDestroyArray(&value);
  return y;
}

static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, 0);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void Csetminmax_api(const mxArray *prhs[3], const mxArray *plhs[1])
{
  real_T (*value)[1302];
  real_T lowlim;
  real_T uplim;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, true, -1);

  /* Marshall function inputs */
  value = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "value");
  lowlim = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "lowlim");
  uplim = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "uplim");

  /* Invoke the target function */
  Csetminmax(&st, *value, lowlim, uplim);

  /* Marshall function outputs */
  b_emlrt_marshallOut(*value, prhs[0]);
  plhs[0] = prhs[0];
}

/* End of code generation (_coder_Csetminmax_api.c) */
