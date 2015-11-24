/*
 * _coder_cget_latlon_api.c
 *
 * Code generation for function '_coder_cget_latlon_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cget_latlon.h"
#include "_coder_cget_latlon_api.h"
#include "cget_latlon_data.h"

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *blat, const
  char_T *identifier);
static const mxArray *emlrt_marshallOut(const real_T u);

/* Function Definitions */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, 0);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *blat, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(sp, emlrtAlias(blat), &thisId);
  emlrtDestroyArray(&blat);
  return y;
}

static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m0);
  return y;
}

void cget_latlon_api(const mxArray * const prhs[6], const mxArray *plhs[2])
{
  real_T blat;
  real_T blon;
  real_T ranges;
  real_T azm;
  real_T major_axis;
  real_T esquared;
  real_T plon;
  real_T plat;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  blat = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "blat");
  blon = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "blon");
  ranges = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "ranges");
  azm = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "azm");
  major_axis = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "major_axis");
  esquared = emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "esquared");

  /* Invoke the target function */
  cget_latlon(&st, blat, blon, ranges, azm, major_axis, esquared, &plat, &plon);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(plat);
  plhs[1] = emlrt_marshallOut(plon);
}

/* End of code generation (_coder_cget_latlon_api.c) */
