/*
 * ceil.c
 *
 * Code generation for function 'ceil'
 *
 * C source code generated on: Sun Aug 23 01:10:26 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GPReduce.h"
#include "ceil.h"

/* Function Definitions */
void b_ceil(real_T x[3])
{
  int32_T k;
  for (k = 0; k < 3; k++) {
    x[k] = muDoubleScalarCeil(x[k]);
  }
}

/* End of code generation (ceil.c) */
