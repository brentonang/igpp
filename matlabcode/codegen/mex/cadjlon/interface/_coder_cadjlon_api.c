/*
 * _coder_cadjlon_api.c
 *
 * Code generation for function '_coder_cadjlon_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cadjlon.h"
#include "_coder_cadjlon_api.h"
#include "cadjlon_data.h"

/* Function Declarations */
static real32_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real32_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);
static real32_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *theta,
  const char_T *identifier);
static const mxArray *emlrt_marshallOut(const real32_T u);

/* Function Definitions */
static real32_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real32_T y;
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real32_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  real32_T ret;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "single", false, 0U, 0);
  ret = *(real32_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real32_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *theta,
  const char_T *identifier)
{
  real32_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(sp, emlrtAlias(theta), &thisId);
  emlrtDestroyArray(&theta);
  return y;
}

static const mxArray *emlrt_marshallOut(const real32_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
  *(real32_T *)mxGetData(m0) = u;
  emlrtAssign(&y, m0);
  return y;
}

void cadjlon_api(const mxArray * const prhs[1], const mxArray *plhs[1])
{
  real32_T theta;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  theta = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "theta");

  /* Invoke the target function */
  cadjlon(&st, &theta);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(theta);
}

/* End of code generation (_coder_cadjlon_api.c) */
