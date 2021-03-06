/*
 * eml_error.c
 *
 * Code generation for function 'eml_error'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cget_latlon.h"
#include "eml_error.h"

/* Variable Definitions */
static emlrtRTEInfo emlrtRTEI = { 20, 5, "eml_error",
  "C:\\Program Files\\MATLAB\\MATLAB Production Server\\R2015a\\toolbox\\eml\\lib\\matlab\\eml\\eml_error.m"
};

/* Function Definitions */
void b_eml_error(const emlrtStack *sp)
{
  static const char_T varargin_1[4] = { 'a', 'c', 'o', 's' };

  emlrtErrorWithMessageIdR2012b(sp, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
    3, 4, 4, varargin_1);
}

void eml_error(const emlrtStack *sp)
{
  static const char_T varargin_1[4] = { 's', 'q', 'r', 't' };

  emlrtErrorWithMessageIdR2012b(sp, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
    3, 4, 4, varargin_1);
}

/* End of code generation (eml_error.c) */
