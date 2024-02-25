/*
 * GPExpand.c
 *
 * Code generation for function 'GPExpand'
 *
 * C source code generated on: Sun Aug 23 01:10:53 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GPExpand.h"
#include "sum.h"
#include "floor.h"
#include "GPExpand_emxutil.h"
#include "GPExpand_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 68, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtRTEInfo emlrtRTEI = { 4, 20, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtRTEInfo c_emlrtRTEI = { 54, 3, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtBCInfo emlrtBCI = { -1, -1, 55, 6, "I2", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 55, 17, "I2", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 55, 28, "I2", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtECInfo emlrtECI = { -1, 55, 3, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtECInfo b_emlrtECI = { -1, 56, 3, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtBCInfo d_emlrtBCI = { -1, -1, 56, 38, "I2", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 56, 26, "I2", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtECInfo c_emlrtECI = { -1, 56, 23, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtECInfo d_emlrtECI = { -1, 57, 3, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtBCInfo f_emlrtBCI = { -1, -1, 57, 40, "I2", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 57, 28, "I2", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtECInfo e_emlrtECI = { -1, 57, 23, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtECInfo f_emlrtECI = { -1, 58, 3, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtBCInfo h_emlrtBCI = { -1, -1, 58, 42, "I2", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtBCInfo i_emlrtBCI = { -1, -1, 58, 30, "I2", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtECInfo g_emlrtECI = { -1, 58, 23, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtRTEInfo e_emlrtRTEI = { 61, 3, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtRTEInfo f_emlrtRTEI = { 62, 4, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtRTEInfo g_emlrtRTEI = { 63, 5, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtECInfo h_emlrtECI = { 3, 67, 10, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m"
};

static emlrtDCInfo emlrtDCI = { 19, 17, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  1 };

static emlrtDCInfo b_emlrtDCI = { 19, 17, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  4 };

static emlrtBCInfo j_emlrtBCI = { -1, -1, 67, 12, "I", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtDCInfo c_emlrtDCI = { 67, 12, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  1 };

static emlrtBCInfo k_emlrtBCI = { -1, -1, 67, 25, "I", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtDCInfo d_emlrtDCI = { 67, 25, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  1 };

static emlrtBCInfo l_emlrtBCI = { -1, -1, 67, 38, "I", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtDCInfo e_emlrtDCI = { 67, 38, "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  1 };

static emlrtBCInfo m_emlrtBCI = { -1, -1, 68, 14, "IResult", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtBCInfo n_emlrtBCI = { -1, -1, 68, 21, "IResult", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

static emlrtBCInfo o_emlrtBCI = { -1, -1, 68, 27, "IResult", "GPExpand",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPExpand.m",
  0 };

/* Function Definitions */
void GPExpand(const emxArray_real_T *I, const real_T newdim[3], emxArray_real_T *
              IResult)
{
  int16_T dim[3];
  int32_T i0;
  real_T b;
  int32_T loop_ub;
  real_T Wt3[125];
  int32_T i;
  int32_T i1;
  static const real_T dv0[5] = { 0.05, 0.25, 0.4, 0.25, 0.05 };

  int16_T iv0[3];
  emxArray_real_T *I2;
  int32_T i2;
  int32_T tmp_size_idx_0;
  int32_T b_loop_ub;
  int32_T i3;
  int32_T idx;
  int32_T b_tmp_size_idx_0;
  int32_T tmp_data[257];
  int32_T b_tmp_data[257];
  int32_T c_tmp_data[257];
  int32_T b_Wt3[3];
  int32_T d_tmp_data[258];
  int32_T e_tmp_data[258];
  emxArray_real_T *b_I2;
  int32_T iv1[3];
  int32_T c_I2[3];
  emxArray_real_T *d_I2;
  emxArray_real_T *e_I2;
  int32_T iv2[3];
  emxArray_real_T *f_I2;
  emxArray_real_T *g_I2;
  int32_T iv3[3];
  emxArray_real_T *h_I2;
  emxArray_real_T *i_I2;
  int32_T iv4[3];
  emxArray_real_T *j_I2;
  emxArray_real_T *k_I2;
  int32_T iv5[2];
  int32_T l_I2[2];
  emxArray_real_T *m_I2;
  emxArray_real_T *n_I2;
  int32_T iv6[2];
  emxArray_real_T *o_I2;
  int32_T j;
  int32_T k;
  real_T pixeli[5];
  int8_T ii_data[5];
  boolean_T exitg3;
  real_T dv1[5];
  boolean_T guard3 = FALSE;
  int32_T f_tmp_data[5];
  int8_T b_ii_data[5];
  int8_T idxi_data[5];
  real_T pixelj[5];
  boolean_T exitg2;
  boolean_T guard2 = FALSE;
  int8_T idxj_data[5];
  real_T pixelk[5];
  boolean_T exitg1;
  boolean_T guard1 = FALSE;
  int8_T idxk_data[5];
  int32_T I2_size[3];
  int32_T Wt3_size[3];
  int8_T b_idxk_data[5];
  int8_T b_idxi_data[5];
  int8_T b_idxj_data[5];
  int8_T c_idxk_data[5];
  real_T A_data[125];
  int32_T A[1];
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  Expand an image as per the Gaussian Pyramid. */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  for (i0 = 0; i0 < 3; i0++) {
    dim[i0] = (int16_T)I->size[i0];
  }

  /*  newdim = zeros(1,3);  */
  /*  newdim(1) = dim(1); newdim(2) = dim(2); newdim(3) = dim(3);   */
  /*  for i = 1:numel(dim) */
  /*      if mod(dim(i),2)==0 */
  /*          newdim(i) = dim(i)*2; */
  /*      else  */
  /*          newdim(i) = dim(i)*2-1;  */
  /*      end */
  /*  end */
  for (i0 = 0; i0 < 3; i0++) {
    b = emlrtNonNegativeCheckFastR2012b(newdim[i0], &b_emlrtDCI,
      emlrtRootTLSGlobal);
    emlrtIntegerCheckFastR2012b(b, &emlrtDCI, emlrtRootTLSGlobal);
  }

  i0 = IResult->size[0] * IResult->size[1] * IResult->size[2];
  IResult->size[0] = (int32_T)newdim[0];
  emxEnsureCapacity((emxArray__common *)IResult, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  i0 = IResult->size[0] * IResult->size[1] * IResult->size[2];
  IResult->size[1] = (int32_T)newdim[1];
  emxEnsureCapacity((emxArray__common *)IResult, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  i0 = IResult->size[0] * IResult->size[1] * IResult->size[2];
  IResult->size[2] = (int32_T)newdim[2];
  emxEnsureCapacity((emxArray__common *)IResult, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = (int32_T)newdim[0] * (int32_T)newdim[1] * (int32_T)newdim[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    IResult->data[i0] = 0.0;
  }

  /*  Initialize the array in the beginning .. */
  for (i0 = 0; i0 < 125; i0++) {
    Wt3[i0] = 1.0;
  }

  for (i = 0; i < 5; i++) {
    for (i0 = 0; i0 < 5; i0++) {
      for (i1 = 0; i1 < 5; i1++) {
        Wt3[(i + 5 * i1) + 25 * i0] *= dv0[i];
      }

      for (i1 = 0; i1 < 5; i1++) {
        Wt3[(i1 + 5 * i) + 25 * i0] *= dv0[i];
      }
    }

    for (i0 = 0; i0 < 5; i0++) {
      for (i1 = 0; i1 < 5; i1++) {
        Wt3[(i1 + 5 * i0) + 25 * i] *= dv0[i];
      }
    }

    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
  }

  /* 		%% Pad the boundaries */
  for (i0 = 0; i0 < 3; i0++) {
    iv0[i0] = (int16_T)(dim[i0] + 2);
  }

  emxInit_real_T(&I2, 3, &c_emlrtRTEI, TRUE);
  i0 = I2->size[0] * I2->size[1] * I2->size[2];
  I2->size[0] = iv0[0];
  emxEnsureCapacity((emxArray__common *)I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  i0 = I2->size[0] * I2->size[1] * I2->size[2];
  I2->size[1] = iv0[1];
  emxEnsureCapacity((emxArray__common *)I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  i0 = I2->size[0] * I2->size[1] * I2->size[2];
  I2->size[2] = iv0[2];
  emxEnsureCapacity((emxArray__common *)I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = iv0[0] * iv0[1] * iv0[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    I2->data[i0] = 0.0;
  }

  if (2 > 1 + dim[0]) {
    i0 = 0;
    i1 = 0;
  } else {
    i0 = 1;
    i1 = iv0[0];
    i2 = 1 + dim[0];
    i1 = emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &emlrtBCI,
      emlrtRootTLSGlobal);
  }

  if (2 > 1 + dim[1]) {
    i2 = 0;
    tmp_size_idx_0 = 0;
  } else {
    i2 = 1;
    tmp_size_idx_0 = iv0[1];
    b_loop_ub = 1 + dim[1];
    tmp_size_idx_0 = emlrtDynamicBoundsCheckFastR2012b(b_loop_ub, 1,
      tmp_size_idx_0, &b_emlrtBCI, emlrtRootTLSGlobal);
  }

  if (2 > 1 + dim[2]) {
    b_loop_ub = 0;
    i3 = 0;
  } else {
    b_loop_ub = 1;
    i3 = iv0[2];
    idx = 1 + dim[2];
    i3 = emlrtDynamicBoundsCheckFastR2012b(idx, 1, i3, &c_emlrtBCI,
      emlrtRootTLSGlobal);
  }

  b_tmp_size_idx_0 = i1 - i0;
  loop_ub = i1 - i0;
  for (i1 = 0; i1 < loop_ub; i1++) {
    tmp_data[i1] = i0 + i1;
  }

  i = tmp_size_idx_0 - i2;
  loop_ub = tmp_size_idx_0 - i2;
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_tmp_data[i0] = i2 + i0;
  }

  tmp_size_idx_0 = i3 - b_loop_ub;
  loop_ub = i3 - b_loop_ub;
  for (i0 = 0; i0 < loop_ub; i0++) {
    c_tmp_data[i0] = b_loop_ub + i0;
  }

  b_Wt3[0] = b_tmp_size_idx_0;
  b_Wt3[1] = i;
  b_Wt3[2] = tmp_size_idx_0;
  emlrtSubAssignSizeCheckR2012b(b_Wt3, 3, *(int32_T (*)[3])I->size, 3, &emlrtECI,
    emlrtRootTLSGlobal);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < i; i1++) {
      for (i2 = 0; i2 < b_tmp_size_idx_0; i2++) {
        I2->data[(tmp_data[i2] + I2->size[0] * b_tmp_data[i1]) + I2->size[0] *
          I2->size[1] * c_tmp_data[i0]] = I->data[(i2 + b_tmp_size_idx_0 * i1) +
          b_tmp_size_idx_0 * i * i0];
      }
    }
  }

  loop_ub = I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  b_loop_ub = I2->size[2];
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&b_I2, 3, &emlrtRTEI, TRUE);
  iv1[0] = 1;
  iv1[1] = loop_ub;
  iv1[2] = b_loop_ub;
  loop_ub = I2->size[1];
  b_loop_ub = I2->size[2];
  i0 = b_I2->size[0] * b_I2->size[1] * b_I2->size[2];
  b_I2->size[0] = 1;
  b_I2->size[1] = loop_ub;
  b_I2->size[2] = b_loop_ub;
  emxEnsureCapacity((emxArray__common *)b_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_I2->data[b_I2->size[0] * i1 + b_I2->size[0] * b_I2->size[1] * i0] =
        I2->data[(I2->size[0] * i1 + I2->size[0] * I2->size[1] * i0) + 1];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    c_I2[i0] = b_I2->size[i0];
  }

  emxFree_real_T(&b_I2);
  emxInit_real_T(&d_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv1, 3, c_I2, 3, &b_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[1];
  idx = I2->size[2];
  i0 = d_I2->size[0] * d_I2->size[1] * d_I2->size[2];
  d_I2->size[0] = 1;
  d_I2->size[1] = i;
  d_I2->size[2] = idx;
  emxEnsureCapacity((emxArray__common *)d_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < idx; i0++) {
    for (i1 = 0; i1 < i; i1++) {
      d_I2->data[d_I2->size[0] * i1 + d_I2->size[0] * d_I2->size[1] * i0] =
        I2->data[(I2->size[0] * i1 + I2->size[0] * I2->size[1] * i0) + 1];
    }
  }

  loop_ub = d_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_loop_ub = d_I2->size[1];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      I2->data[I2->size[0] * d_tmp_data[i1] + I2->size[0] * I2->size[1] *
        e_tmp_data[i0]] = d_I2->data[d_I2->size[0] * i1 + d_I2->size[0] *
        d_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&d_I2);
  i0 = I2->size[0];
  i1 = I2->size[0];
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &e_emlrtBCI, emlrtRootTLSGlobal);
  loop_ub = I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  b_loop_ub = I2->size[2];
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&e_I2, 3, &emlrtRTEI, TRUE);
  i0 = I2->size[0];
  i1 = (int32_T)((real_T)I2->size[0] - 1.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &d_emlrtBCI, emlrtRootTLSGlobal);
  iv2[0] = 1;
  iv2[1] = loop_ub;
  iv2[2] = b_loop_ub;
  loop_ub = I2->size[1];
  b_loop_ub = I2->size[2];
  i = I2->size[0];
  i0 = e_I2->size[0] * e_I2->size[1] * e_I2->size[2];
  e_I2->size[0] = 1;
  e_I2->size[1] = loop_ub;
  e_I2->size[2] = b_loop_ub;
  emxEnsureCapacity((emxArray__common *)e_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      e_I2->data[e_I2->size[0] * i1 + e_I2->size[0] * e_I2->size[1] * i0] =
        I2->data[((i + I2->size[0] * i1) + I2->size[0] * I2->size[1] * i0) - 2];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    c_I2[i0] = e_I2->size[i0];
  }

  emxFree_real_T(&e_I2);
  emxInit_real_T(&f_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv2, 3, c_I2, 3, &c_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[0] - 1;
  idx = I2->size[0] - 2;
  tmp_size_idx_0 = I2->size[1];
  b_tmp_size_idx_0 = I2->size[2];
  i0 = f_I2->size[0] * f_I2->size[1] * f_I2->size[2];
  f_I2->size[0] = 1;
  f_I2->size[1] = tmp_size_idx_0;
  f_I2->size[2] = b_tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)f_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      f_I2->data[f_I2->size[0] * i1 + f_I2->size[0] * f_I2->size[1] * i0] =
        I2->data[(idx + I2->size[0] * i1) + I2->size[0] * I2->size[1] * i0];
    }
  }

  loop_ub = f_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_loop_ub = f_I2->size[1];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      I2->data[(i + I2->size[0] * d_tmp_data[i1]) + I2->size[0] * I2->size[1] *
        e_tmp_data[i0]] = f_I2->data[f_I2->size[0] * i1 + f_I2->size[0] *
        f_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&f_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  b_loop_ub = I2->size[2];
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&g_I2, 3, &emlrtRTEI, TRUE);
  iv3[0] = loop_ub;
  iv3[1] = 1;
  iv3[2] = b_loop_ub;
  loop_ub = I2->size[0];
  b_loop_ub = I2->size[2];
  i0 = g_I2->size[0] * g_I2->size[1] * g_I2->size[2];
  g_I2->size[0] = loop_ub;
  g_I2->size[1] = 1;
  g_I2->size[2] = b_loop_ub;
  emxEnsureCapacity((emxArray__common *)g_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      g_I2->data[i1 + g_I2->size[0] * g_I2->size[1] * i0] = I2->data[(i1 +
        I2->size[0]) + I2->size[0] * I2->size[1] * i0];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    c_I2[i0] = g_I2->size[i0];
  }

  emxFree_real_T(&g_I2);
  emxInit_real_T(&h_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv3, 3, c_I2, 3, &d_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[0];
  idx = I2->size[2];
  i0 = h_I2->size[0] * h_I2->size[1] * h_I2->size[2];
  h_I2->size[0] = i;
  h_I2->size[1] = 1;
  h_I2->size[2] = idx;
  emxEnsureCapacity((emxArray__common *)h_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < idx; i0++) {
    for (i1 = 0; i1 < i; i1++) {
      h_I2->data[i1 + h_I2->size[0] * h_I2->size[1] * i0] = I2->data[(i1 +
        I2->size[0]) + I2->size[0] * I2->size[1] * i0];
    }
  }

  loop_ub = h_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_loop_ub = h_I2->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      I2->data[d_tmp_data[i1] + I2->size[0] * I2->size[1] * e_tmp_data[i0]] =
        h_I2->data[i1 + h_I2->size[0] * h_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&h_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  i0 = I2->size[1];
  i1 = I2->size[1];
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &g_emlrtBCI, emlrtRootTLSGlobal);
  b_loop_ub = I2->size[2];
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&i_I2, 3, &emlrtRTEI, TRUE);
  i0 = I2->size[1];
  i1 = (int32_T)((real_T)I2->size[1] - 1.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &f_emlrtBCI, emlrtRootTLSGlobal);
  iv4[0] = loop_ub;
  iv4[1] = 1;
  iv4[2] = b_loop_ub;
  loop_ub = I2->size[0];
  b_loop_ub = I2->size[2];
  i = I2->size[1];
  i0 = i_I2->size[0] * i_I2->size[1] * i_I2->size[2];
  i_I2->size[0] = loop_ub;
  i_I2->size[1] = 1;
  i_I2->size[2] = b_loop_ub;
  emxEnsureCapacity((emxArray__common *)i_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      i_I2->data[i1 + i_I2->size[0] * i_I2->size[1] * i0] = I2->data[(i1 +
        I2->size[0] * (i - 2)) + I2->size[0] * I2->size[1] * i0];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    c_I2[i0] = i_I2->size[i0];
  }

  emxFree_real_T(&i_I2);
  emxInit_real_T(&j_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv4, 3, c_I2, 3, &e_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[1] - 1;
  idx = I2->size[0];
  tmp_size_idx_0 = I2->size[1] - 2;
  b_tmp_size_idx_0 = I2->size[2];
  i0 = j_I2->size[0] * j_I2->size[1] * j_I2->size[2];
  j_I2->size[0] = idx;
  j_I2->size[1] = 1;
  j_I2->size[2] = b_tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)j_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < idx; i1++) {
      j_I2->data[i1 + j_I2->size[0] * j_I2->size[1] * i0] = I2->data[(i1 +
        I2->size[0] * tmp_size_idx_0) + I2->size[0] * I2->size[1] * i0];
    }
  }

  loop_ub = j_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_loop_ub = j_I2->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      I2->data[(d_tmp_data[i1] + I2->size[0] * i) + I2->size[0] * I2->size[1] *
        e_tmp_data[i0]] = j_I2->data[i1 + j_I2->size[0] * j_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&j_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  b_loop_ub = I2->size[1];
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    e_tmp_data[i0] = i0;
  }

  b_emxInit_real_T(&k_I2, 2, &emlrtRTEI, TRUE);
  iv5[0] = loop_ub;
  iv5[1] = b_loop_ub;
  loop_ub = I2->size[0];
  b_loop_ub = I2->size[1];
  i0 = k_I2->size[0] * k_I2->size[1];
  k_I2->size[0] = loop_ub;
  k_I2->size[1] = b_loop_ub;
  emxEnsureCapacity((emxArray__common *)k_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      k_I2->data[i1 + k_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        I2->size[0] * I2->size[1]];
    }
  }

  for (i0 = 0; i0 < 2; i0++) {
    l_I2[i0] = k_I2->size[i0];
  }

  emxFree_real_T(&k_I2);
  b_emxInit_real_T(&m_I2, 2, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv5, 2, l_I2, 2, &f_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[0];
  idx = I2->size[1];
  i0 = m_I2->size[0] * m_I2->size[1];
  m_I2->size[0] = i;
  m_I2->size[1] = idx;
  emxEnsureCapacity((emxArray__common *)m_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < idx; i0++) {
    for (i1 = 0; i1 < i; i1++) {
      m_I2->data[i1 + m_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        I2->size[0] * I2->size[1]];
    }
  }

  loop_ub = m_I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_loop_ub = m_I2->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      I2->data[d_tmp_data[i1] + I2->size[0] * e_tmp_data[i0]] = m_I2->data[i1 +
        m_I2->size[0] * i0];
    }
  }

  emxFree_real_T(&m_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  b_loop_ub = I2->size[1];
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    e_tmp_data[i0] = i0;
  }

  b_emxInit_real_T(&n_I2, 2, &emlrtRTEI, TRUE);
  i0 = I2->size[2];
  i1 = I2->size[2];
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &i_emlrtBCI, emlrtRootTLSGlobal);
  i0 = I2->size[2];
  i1 = (int32_T)((real_T)I2->size[2] - 1.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &h_emlrtBCI, emlrtRootTLSGlobal);
  iv6[0] = loop_ub;
  iv6[1] = b_loop_ub;
  loop_ub = I2->size[0];
  b_loop_ub = I2->size[1];
  i = I2->size[2];
  i0 = n_I2->size[0] * n_I2->size[1];
  n_I2->size[0] = loop_ub;
  n_I2->size[1] = b_loop_ub;
  emxEnsureCapacity((emxArray__common *)n_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      n_I2->data[i1 + n_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        I2->size[0] * I2->size[1] * (i - 2)];
    }
  }

  for (i0 = 0; i0 < 2; i0++) {
    l_I2[i0] = n_I2->size[i0];
  }

  emxFree_real_T(&n_I2);
  b_emxInit_real_T(&o_I2, 2, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv6, 2, l_I2, 2, &g_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[2] - 1;
  idx = I2->size[0];
  tmp_size_idx_0 = I2->size[1];
  b_tmp_size_idx_0 = I2->size[2] - 2;
  i0 = o_I2->size[0] * o_I2->size[1];
  o_I2->size[0] = idx;
  o_I2->size[1] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)o_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < idx; i1++) {
      o_I2->data[i1 + o_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        I2->size[0] * I2->size[1] * b_tmp_size_idx_0];
    }
  }

  loop_ub = o_I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_loop_ub = o_I2->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      I2->data[(d_tmp_data[i1] + I2->size[0] * e_tmp_data[i0]) + I2->size[0] *
        I2->size[1] * i] = o_I2->data[i1 + o_I2->size[0] * i0];
    }
  }

  emxFree_real_T(&o_I2);

  /* clear I2; */
  emlrtForLoopVectorCheckR2012b(0.0, 1.0, newdim[0] - 1.0, mxDOUBLE_CLASS,
    (int32_T)((newdim[0] - 1.0) + 1.0), &e_emlrtRTEI, emlrtRootTLSGlobal);
  i = 0;
  while (i <= (int32_T)((newdim[0] - 1.0) + 1.0) - 1) {
    emlrtForLoopVectorCheckR2012b(0.0, 1.0, newdim[1] - 1.0, mxDOUBLE_CLASS,
      (int32_T)((newdim[1] - 1.0) + 1.0), &f_emlrtRTEI, emlrtRootTLSGlobal);
    j = 0;
    while (j <= (int32_T)((newdim[1] - 1.0) + 1.0) - 1) {
      emlrtForLoopVectorCheckR2012b(0.0, 1.0, newdim[2] - 1.0, mxDOUBLE_CLASS,
        (int32_T)((newdim[2] - 1.0) + 1.0), &g_emlrtRTEI, emlrtRootTLSGlobal);
      k = 0;
      while (k <= (int32_T)((newdim[2] - 1.0) + 1.0) - 1) {
        for (i0 = 0; i0 < 5; i0++) {
          pixeli[i0] = ((real_T)i - (-2.0 + (real_T)i0)) / 2.0 + 2.0;
        }

        idx = 0;
        tmp_size_idx_0 = 1;
        exitg3 = FALSE;
        while ((exitg3 == FALSE) && (tmp_size_idx_0 < 6)) {
          for (i0 = 0; i0 < 5; i0++) {
            dv1[i0] = pixeli[i0];
          }

          b_floor(dv1);
          guard3 = FALSE;
          if (dv1[tmp_size_idx_0 - 1] == pixeli[tmp_size_idx_0 - 1]) {
            idx++;
            ii_data[idx - 1] = (int8_T)tmp_size_idx_0;
            if (idx >= 5) {
              exitg3 = TRUE;
            } else {
              guard3 = TRUE;
            }
          } else {
            guard3 = TRUE;
          }

          if (guard3 == TRUE) {
            tmp_size_idx_0++;
          }
        }

        if (1 > idx) {
          loop_ub = 0;
        } else {
          loop_ub = idx;
        }

        for (i0 = 0; i0 < loop_ub; i0++) {
          f_tmp_data[i0] = 1 + i0;
        }

        for (i0 = 0; i0 < loop_ub; i0++) {
          i1 = 0;
          while (i1 <= 0) {
            b_ii_data[i0] = ii_data[f_tmp_data[i0] - 1];
            i1 = 1;
          }
        }

        for (i0 = 0; i0 < loop_ub; i0++) {
          ii_data[i0] = b_ii_data[i0];
        }

        for (i0 = 0; i0 < loop_ub; i0++) {
          idxi_data[i0] = ii_data[i0];
        }

        for (i0 = 0; i0 < 5; i0++) {
          pixelj[i0] = ((real_T)j - (-2.0 + (real_T)i0)) / 2.0 + 2.0;
        }

        idx = 0;
        tmp_size_idx_0 = 1;
        exitg2 = FALSE;
        while ((exitg2 == FALSE) && (tmp_size_idx_0 < 6)) {
          for (i0 = 0; i0 < 5; i0++) {
            dv1[i0] = pixelj[i0];
          }

          b_floor(dv1);
          guard2 = FALSE;
          if (dv1[tmp_size_idx_0 - 1] == pixelj[tmp_size_idx_0 - 1]) {
            idx++;
            ii_data[idx - 1] = (int8_T)tmp_size_idx_0;
            if (idx >= 5) {
              exitg2 = TRUE;
            } else {
              guard2 = TRUE;
            }
          } else {
            guard2 = TRUE;
          }

          if (guard2 == TRUE) {
            tmp_size_idx_0++;
          }
        }

        if (1 > idx) {
          b_loop_ub = 0;
        } else {
          b_loop_ub = idx;
        }

        for (i0 = 0; i0 < b_loop_ub; i0++) {
          f_tmp_data[i0] = 1 + i0;
        }

        for (i0 = 0; i0 < b_loop_ub; i0++) {
          i1 = 0;
          while (i1 <= 0) {
            b_ii_data[i0] = ii_data[f_tmp_data[i0] - 1];
            i1 = 1;
          }
        }

        for (i0 = 0; i0 < b_loop_ub; i0++) {
          ii_data[i0] = b_ii_data[i0];
        }

        for (i0 = 0; i0 < b_loop_ub; i0++) {
          idxj_data[i0] = ii_data[i0];
        }

        for (i0 = 0; i0 < 5; i0++) {
          pixelk[i0] = ((real_T)k - (-2.0 + (real_T)i0)) / 2.0 + 2.0;
        }

        idx = 0;
        tmp_size_idx_0 = 1;
        exitg1 = FALSE;
        while ((exitg1 == FALSE) && (tmp_size_idx_0 < 6)) {
          for (i0 = 0; i0 < 5; i0++) {
            dv1[i0] = pixelk[i0];
          }

          b_floor(dv1);
          guard1 = FALSE;
          if (dv1[tmp_size_idx_0 - 1] == pixelk[tmp_size_idx_0 - 1]) {
            idx++;
            ii_data[idx - 1] = (int8_T)tmp_size_idx_0;
            if (idx >= 5) {
              exitg1 = TRUE;
            } else {
              guard1 = TRUE;
            }
          } else {
            guard1 = TRUE;
          }

          if (guard1 == TRUE) {
            tmp_size_idx_0++;
          }
        }

        if (1 > idx) {
          b_tmp_size_idx_0 = 0;
        } else {
          b_tmp_size_idx_0 = idx;
        }

        for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
          f_tmp_data[i0] = 1 + i0;
        }

        for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
          i1 = 0;
          while (i1 <= 0) {
            b_ii_data[i0] = ii_data[f_tmp_data[i0] - 1];
            i1 = 1;
          }
        }

        for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
          ii_data[i0] = b_ii_data[i0];
        }

        for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
          idxk_data[i0] = ii_data[i0];
        }

        for (i0 = 0; i0 < loop_ub; i0++) {
          i1 = I2->size[0];
          b = pixeli[idxi_data[i0] - 1];
          i2 = (int32_T)emlrtIntegerCheckFastR2012b(b, &c_emlrtDCI,
            emlrtRootTLSGlobal);
          emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &j_emlrtBCI,
            emlrtRootTLSGlobal);
        }

        for (i0 = 0; i0 < b_loop_ub; i0++) {
          i1 = I2->size[1];
          b = pixelj[idxj_data[i0] - 1];
          i2 = (int32_T)emlrtIntegerCheckFastR2012b(b, &d_emlrtDCI,
            emlrtRootTLSGlobal);
          emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &k_emlrtBCI,
            emlrtRootTLSGlobal);
        }

        for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
          i1 = I2->size[2];
          b = pixelk[idxk_data[i0] - 1];
          i2 = (int32_T)emlrtIntegerCheckFastR2012b(b, &e_emlrtDCI,
            emlrtRootTLSGlobal);
          emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &l_emlrtBCI,
            emlrtRootTLSGlobal);
        }

        I2_size[0] = loop_ub;
        I2_size[1] = b_loop_ub;
        I2_size[2] = b_tmp_size_idx_0;
        Wt3_size[0] = loop_ub;
        Wt3_size[1] = b_loop_ub;
        Wt3_size[2] = b_tmp_size_idx_0;
        for (i0 = 0; i0 < 3; i0++) {
          c_I2[i0] = I2_size[i0];
          b_Wt3[i0] = Wt3_size[i0];
        }

        emlrtSizeEqCheckNDR2012b(c_I2, b_Wt3, &h_emlrtECI, emlrtRootTLSGlobal);
        for (i0 = 0; i0 < loop_ub; i0++) {
          ii_data[i0] = idxi_data[i0];
        }

        for (i0 = 0; i0 < b_loop_ub; i0++) {
          b_ii_data[i0] = idxj_data[i0];
        }

        for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
          b_idxk_data[i0] = idxk_data[i0];
        }

        for (i0 = 0; i0 < loop_ub; i0++) {
          b_idxi_data[i0] = idxi_data[i0];
        }

        for (i0 = 0; i0 < b_loop_ub; i0++) {
          b_idxj_data[i0] = idxj_data[i0];
        }

        for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
          c_idxk_data[i0] = idxk_data[i0];
        }

        for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
          for (i1 = 0; i1 < b_loop_ub; i1++) {
            for (i2 = 0; i2 < loop_ub; i2++) {
              A_data[(i2 + loop_ub * i1) + loop_ub * b_loop_ub * i0] = I2->data
                [(((int32_T)pixeli[ii_data[i2] - 1] + I2->size[0] * ((int32_T)
                    pixelj[b_ii_data[i1] - 1] - 1)) + I2->size[0] * I2->size[1] *
                  ((int32_T)pixelk[b_idxk_data[i0] - 1] - 1)) - 1] * Wt3
                [((b_idxi_data[i2] + 5 * (b_idxj_data[i1] - 1)) + 25 *
                  (c_idxk_data[i0] - 1)) - 1];
            }
          }
        }

        emlrtPushRtStackR2012b(&emlrtRSI, emlrtRootTLSGlobal);
        A[0] = loop_ub * b_loop_ub * b_tmp_size_idx_0;
        b = sum(A_data, A);
        i0 = IResult->size[0];
        i1 = (int32_T)((real_T)i + 1.0);
        i2 = IResult->size[1];
        tmp_size_idx_0 = (int32_T)((real_T)j + 1.0);
        b_loop_ub = IResult->size[2];
        i3 = (int32_T)((real_T)k + 1.0);
        IResult->data[((emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &m_emlrtBCI,
          emlrtRootTLSGlobal) + IResult->size[0] *
                        (emlrtDynamicBoundsCheckFastR2012b(tmp_size_idx_0, 1, i2,
          &n_emlrtBCI, emlrtRootTLSGlobal) - 1)) + IResult->size[0] *
                       IResult->size[1] * (emlrtDynamicBoundsCheckFastR2012b(i3,
          1, b_loop_ub, &o_emlrtBCI, emlrtRootTLSGlobal) - 1)) - 1] = 8.0 * b;
        emlrtPopRtStackR2012b(&emlrtRSI, emlrtRootTLSGlobal);
        k++;
        emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar,
          emlrtRootTLSGlobal);
      }

      j++;
      emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
    }

    i++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
  }

  emxFree_real_T(&I2);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (GPExpand.c) */
