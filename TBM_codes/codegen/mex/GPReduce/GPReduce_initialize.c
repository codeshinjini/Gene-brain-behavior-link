/*
 * GPReduce_initialize.c
 *
 * Code generation for function 'GPReduce_initialize'
 *
 * C source code generated on: Sun Aug 23 01:10:26 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GPReduce.h"
#include "GPReduce_initialize.h"
#include "GPReduce_data.h"

/* Function Definitions */
void GPReduce_initialize(emlrtContext *aContext)
{
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, FALSE, 0U, 0);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (GPReduce_initialize.c) */
