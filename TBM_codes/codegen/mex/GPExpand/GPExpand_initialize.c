/*
 * GPExpand_initialize.c
 *
 * Code generation for function 'GPExpand_initialize'
 *
 * C source code generated on: Sun Aug 23 01:10:53 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GPExpand.h"
#include "GPExpand_initialize.h"
#include "GPExpand_data.h"

/* Function Definitions */
void GPExpand_initialize(emlrtContext *aContext)
{
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, FALSE, 0U, 0);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (GPExpand_initialize.c) */
