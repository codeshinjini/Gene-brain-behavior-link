/*
 * GPExpand_terminate.c
 *
 * Code generation for function 'GPExpand_terminate'
 *
 * C source code generated on: Sun Aug 23 01:10:53 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GPExpand.h"
#include "GPExpand_terminate.h"

/* Function Definitions */
void GPExpand_atexit(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void GPExpand_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (GPExpand_terminate.c) */
