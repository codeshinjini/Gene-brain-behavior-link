/*
 * sum.c
 *
 * Code generation for function 'sum'
 *
 * C source code generated on: Sun Aug 23 01:10:26 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GPReduce.h"
#include "sum.h"

/* Function Definitions */
real_T sum(const real_T x[125])
{
  real_T y;
  int32_T k;
  y = x[0];
  for (k = 0; k < 124; k++) {
    y += x[k + 1];
  }

  return y;
}

/* End of code generation (sum.c) */
