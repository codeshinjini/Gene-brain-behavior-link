/*
 * sum.c
 *
 * Code generation for function 'sum'
 *
 * C source code generated on: Sun Aug 23 01:10:53 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GPExpand.h"
#include "sum.h"

/* Variable Definitions */
static emlrtRSInfo b_emlrtRSI = { 17, "sum",
  "/afs/ece.cmu.edu/support/matlab/2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtRSInfo c_emlrtRSI = { 20, "sum",
  "/afs/ece.cmu.edu/support/matlab/2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtMCInfo emlrtMCI = { 18, 9, "sum",
  "/afs/ece.cmu.edu/support/matlab/2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtMCInfo b_emlrtMCI = { 17, 19, "sum",
  "/afs/ece.cmu.edu/support/matlab/2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtMCInfo c_emlrtMCI = { 23, 9, "sum",
  "/afs/ece.cmu.edu/support/matlab/2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtMCInfo d_emlrtMCI = { 20, 19, "sum",
  "/afs/ece.cmu.edu/support/matlab/2013a/toolbox/eml/lib/matlab/datafun/sum.m" };

/* Function Declarations */
static void error(const mxArray *b, emlrtMCInfo *location);
static const mxArray *message(const mxArray *b, emlrtMCInfo *location);

/* Function Definitions */
static void error(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 0, NULL, 1, &pArray, "error", TRUE,
                        location);
}

static const mxArray *message(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  const mxArray *m3;
  pArray = b;
  return emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 1, &m3, 1, &pArray, "message",
    TRUE, location);
}

real_T sum(const real_T x_data[125], const int32_T x_size[1])
{
  real_T y;
  boolean_T p;
  boolean_T b_p;
  int32_T i;
  int32_T exitg1;
  int32_T b_i;
  const mxArray *b_y;
  static const int32_T iv7[2] = { 1, 30 };

  const mxArray *m0;
  char_T cv0[30];
  static const char_T cv1[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 's', 'u', 'm', '_', 's', 'p', 'e', 'c', 'i', 'a',
    'l', 'E', 'm', 'p', 't', 'y' };

  const mxArray *c_y;
  static const int32_T iv8[2] = { 1, 36 };

  char_T cv2[36];
  static const char_T cv3[36] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'a', 'u', 't', 'o', 'D', 'i', 'm', 'I', 'n', 'c',
    'o', 'm', 'p', 'a', 't', 'i', 'b', 'i', 'l', 'i', 't', 'y' };

  p = FALSE;
  b_p = FALSE;
  i = 0;
  do {
    exitg1 = 0;
    if (i < 2) {
      if (i + 1 <= 1) {
        b_i = x_size[0];
      } else {
        b_i = 1;
      }

      if (b_i != 0) {
        exitg1 = 1;
      } else {
        i++;
      }
    } else {
      b_p = TRUE;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (!b_p) {
  } else {
    p = TRUE;
  }

  if (!p) {
  } else {
    emlrtPushRtStackR2012b(&b_emlrtRSI, emlrtRootTLSGlobal);
    b_y = NULL;
    m0 = mxCreateCharArray(2, iv7);
    for (i = 0; i < 30; i++) {
      cv0[i] = cv1[i];
    }

    emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 30, m0, cv0);
    emlrtAssign(&b_y, m0);
    error(message(b_y, &emlrtMCI), &b_emlrtMCI);
    emlrtPopRtStackR2012b(&b_emlrtRSI, emlrtRootTLSGlobal);
  }

  if ((x_size[0] == 1) || (x_size[0] != 1)) {
    p = TRUE;
  } else {
    p = FALSE;
  }

  if (p) {
  } else {
    emlrtPushRtStackR2012b(&c_emlrtRSI, emlrtRootTLSGlobal);
    c_y = NULL;
    m0 = mxCreateCharArray(2, iv8);
    for (i = 0; i < 36; i++) {
      cv2[i] = cv3[i];
    }

    emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 36, m0, cv2);
    emlrtAssign(&c_y, m0);
    error(message(c_y, &c_emlrtMCI), &d_emlrtMCI);
    emlrtPopRtStackR2012b(&c_emlrtRSI, emlrtRootTLSGlobal);
  }

  if (x_size[0] == 0) {
    y = 0.0;
  } else {
    y = x_data[0];
    for (i = 2; i <= x_size[0]; i++) {
      y += x_data[i - 1];
    }
  }

  return y;
}

/* End of code generation (sum.c) */
