/*
 * floor.c
 *
 * Code generation for function 'floor'
 *
 * C source code generated on: Sun Aug 23 01:10:53 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GPExpand.h"
#include "floor.h"

/* Function Definitions */
void b_floor(real_T x[5])
{
  int32_T k;
  for (k = 0; k < 5; k++) {
    x[k] = muDoubleScalarFloor(x[k]);
  }
}

/* End of code generation (floor.c) */
