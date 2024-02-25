/*
 * GPReduce.c
 *
 * Code generation for function 'GPReduce'
 *
 * C source code generated on: Sun Aug 23 01:10:26 2015
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GPReduce.h"
#include "sum.h"
#include "GPReduce_emxutil.h"
#include "ceil.h"
#include "GPReduce_data.h"

/* Variable Definitions */
static emlrtRTEInfo emlrtRTEI = { 4, 20, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtRTEInfo b_emlrtRTEI = { 44, 3, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtBCInfo emlrtBCI = { -1, -1, 45, 6, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 45, 17, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 45, 28, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtECInfo emlrtECI = { -1, 45, 3, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtECInfo b_emlrtECI = { -1, 46, 3, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtECInfo c_emlrtECI = { -1, 46, 23, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtBCInfo d_emlrtBCI = { -1, -1, 46, 58, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 46, 46, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtECInfo d_emlrtECI = { -1, 46, 43, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtBCInfo f_emlrtBCI = { -1, -1, 46, 86, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 46, 72, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtECInfo e_emlrtECI = { -1, 46, 69, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtECInfo f_emlrtECI = { -1, 47, 3, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtECInfo g_emlrtECI = { -1, 47, 23, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtBCInfo h_emlrtBCI = { -1, -1, 47, 60, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo i_emlrtBCI = { -1, -1, 47, 48, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtECInfo h_emlrtECI = { -1, 47, 43, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtBCInfo j_emlrtBCI = { -1, -1, 47, 88, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo k_emlrtBCI = { -1, -1, 47, 74, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtECInfo i_emlrtECI = { -1, 47, 69, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtECInfo j_emlrtECI = { -1, 48, 3, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtECInfo k_emlrtECI = { -1, 48, 23, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtBCInfo l_emlrtBCI = { -1, -1, 48, 62, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo m_emlrtBCI = { -1, -1, 48, 50, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtECInfo l_emlrtECI = { -1, 48, 43, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtBCInfo n_emlrtBCI = { -1, -1, 48, 90, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo o_emlrtBCI = { -1, -1, 48, 76, "I2", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtECInfo m_emlrtECI = { -1, 48, 69, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m"
};

static emlrtDCInfo emlrtDCI = { 10, 17, "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  1 };

static emlrtBCInfo p_emlrtBCI = { -1, -1, 54, 12, "I", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo q_emlrtBCI = { -1, -1, 54, 20, "I", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo r_emlrtBCI = { -1, -1, 54, 28, "I", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo s_emlrtBCI = { -1, -1, 55, 14, "IResult", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo t_emlrtBCI = { -1, -1, 55, 18, "IResult", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

static emlrtBCInfo u_emlrtBCI = { -1, -1, 55, 22, "IResult", "GPReduce",
  "/afs/.ece.cmu.edu/project/cbi/users/gustavor/Shinjini/Autism/TBM_features/GPReduce.m",
  0 };

/* Function Definitions */
void GPReduce(const emxArray_real_T *I, emxArray_real_T *IResult)
{
  int16_T dim[3];
  int32_T i0;
  real_T newdim[3];
  int32_T loop_ub;
  real_T Wt3[125];
  int32_T i;
  int32_T i1;
  static const real_T dv0[5] = { 0.05, 0.25, 0.4, 0.25, 0.05 };

  int16_T iv0[3];
  emxArray_real_T *I2;
  int32_T i2;
  int32_T i3;
  int32_T i4;
  int32_T i5;
  int32_T tmp_size_idx_0;
  int32_T b_tmp_size_idx_0;
  int32_T tmp_data[258];
  int32_T c_tmp_size_idx_0;
  int32_T b_tmp_data[258];
  int32_T c_tmp_data[258];
  int32_T b_I2[3];
  int32_T d_tmp_data[260];
  int32_T e_tmp_data[260];
  emxArray_real_T *c_I2;
  int32_T iv1[3];
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
  int32_T iv5[3];
  emxArray_real_T *l_I2;
  emxArray_real_T *m_I2;
  int32_T iv6[3];
  emxArray_real_T *n_I2;
  emxArray_real_T *o_I2;
  int32_T iv7[3];
  emxArray_real_T *p_I2;
  emxArray_real_T *q_I2;
  int32_T iv8[3];
  emxArray_real_T *r_I2;
  emxArray_real_T *s_I2;
  int32_T iv9[2];
  int32_T t_I2[2];
  emxArray_real_T *u_I2;
  emxArray_real_T *v_I2;
  int32_T iv10[2];
  emxArray_real_T *w_I2;
  emxArray_real_T *x_I2;
  int32_T iv11[2];
  emxArray_real_T *y_I2;
  emxArray_real_T *ab_I2;
  int32_T iv12[2];
  emxArray_real_T *bb_I2;
  real_T y;
  real_T b_y;
  real_T c_y;
  real_T A[125];
  int32_T i6;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  Reduce an image applying Gaussian Pyramid. */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  for (i0 = 0; i0 < 3; i0++) {
    dim[i0] = (int16_T)I->size[i0];
  }

  for (i0 = 0; i0 < 3; i0++) {
    newdim[i0] = (real_T)dim[i0] * 0.5;
  }

  b_ceil(newdim);
  for (i0 = 0; i0 < 3; i0++) {
    emlrtIntegerCheckFastR2012b(newdim[i0], &emlrtDCI, emlrtRootTLSGlobal);
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

  /* 		%% Pad the boundaries. */
  for (i0 = 0; i0 < 3; i0++) {
    iv0[i0] = (int16_T)(dim[i0] + 4);
  }

  emxInit_real_T(&I2, 3, &b_emlrtRTEI, TRUE);
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

  if (3 > 2 + dim[0]) {
    i0 = 0;
    i1 = 0;
  } else {
    i0 = 2;
    i1 = iv0[0];
    i2 = 2 + dim[0];
    i1 = emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &emlrtBCI,
      emlrtRootTLSGlobal);
  }

  if (3 > 2 + dim[1]) {
    i2 = 0;
    i3 = 0;
  } else {
    i2 = 2;
    i3 = iv0[1];
    i4 = 2 + dim[1];
    i3 = emlrtDynamicBoundsCheckFastR2012b(i4, 1, i3, &b_emlrtBCI,
      emlrtRootTLSGlobal);
  }

  if (3 > 2 + dim[2]) {
    i4 = 0;
    i5 = 0;
  } else {
    i4 = 2;
    i5 = iv0[2];
    tmp_size_idx_0 = 2 + dim[2];
    i5 = emlrtDynamicBoundsCheckFastR2012b(tmp_size_idx_0, 1, i5, &c_emlrtBCI,
      emlrtRootTLSGlobal);
  }

  b_tmp_size_idx_0 = i1 - i0;
  loop_ub = i1 - i0;
  for (i1 = 0; i1 < loop_ub; i1++) {
    tmp_data[i1] = i0 + i1;
  }

  c_tmp_size_idx_0 = i3 - i2;
  loop_ub = i3 - i2;
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_tmp_data[i0] = i2 + i0;
  }

  tmp_size_idx_0 = i5 - i4;
  loop_ub = i5 - i4;
  for (i0 = 0; i0 < loop_ub; i0++) {
    c_tmp_data[i0] = i4 + i0;
  }

  b_I2[0] = b_tmp_size_idx_0;
  b_I2[1] = c_tmp_size_idx_0;
  b_I2[2] = tmp_size_idx_0;
  emlrtSubAssignSizeCheckR2012b(b_I2, 3, *(int32_T (*)[3])I->size, 3, &emlrtECI,
    emlrtRootTLSGlobal);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < c_tmp_size_idx_0; i1++) {
      for (i2 = 0; i2 < b_tmp_size_idx_0; i2++) {
        I2->data[(tmp_data[i2] + I2->size[0] * b_tmp_data[i1]) + I2->size[0] *
          I2->size[1] * c_tmp_data[i0]] = I->data[(i2 + b_tmp_size_idx_0 * i1) +
          b_tmp_size_idx_0 * c_tmp_size_idx_0 * i0];
      }
    }
  }

  loop_ub = I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  tmp_size_idx_0 = I2->size[2];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&c_I2, 3, &emlrtRTEI, TRUE);
  iv1[0] = 1;
  iv1[1] = loop_ub;
  iv1[2] = tmp_size_idx_0;
  loop_ub = I2->size[1];
  tmp_size_idx_0 = I2->size[2];
  i0 = c_I2->size[0] * c_I2->size[1] * c_I2->size[2];
  c_I2->size[0] = 1;
  c_I2->size[1] = loop_ub;
  c_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)c_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      c_I2->data[c_I2->size[0] * i1 + c_I2->size[0] * c_I2->size[1] * i0] =
        I2->data[(I2->size[0] * i1 + I2->size[0] * I2->size[1] * i0) + 2];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    b_I2[i0] = c_I2->size[i0];
  }

  emxFree_real_T(&c_I2);
  emxInit_real_T(&d_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv1, 3, b_I2, 3, &b_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[1];
  tmp_size_idx_0 = I2->size[2];
  i0 = d_I2->size[0] * d_I2->size[1] * d_I2->size[2];
  d_I2->size[0] = 1;
  d_I2->size[1] = i;
  d_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)d_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < i; i1++) {
      d_I2->data[d_I2->size[0] * i1 + d_I2->size[0] * d_I2->size[1] * i0] =
        I2->data[(I2->size[0] * i1 + I2->size[0] * I2->size[1] * i0) + 2];
    }
  }

  loop_ub = d_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = d_I2->size[1];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[I2->size[0] * d_tmp_data[i1] + I2->size[0] * I2->size[1] *
        e_tmp_data[i0]] = d_I2->data[d_I2->size[0] * i1 + d_I2->size[0] *
        d_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&d_I2);
  loop_ub = I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  tmp_size_idx_0 = I2->size[2];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&e_I2, 3, &emlrtRTEI, TRUE);
  iv2[0] = 1;
  iv2[1] = loop_ub;
  iv2[2] = tmp_size_idx_0;
  loop_ub = I2->size[1];
  tmp_size_idx_0 = I2->size[2];
  i0 = e_I2->size[0] * e_I2->size[1] * e_I2->size[2];
  e_I2->size[0] = 1;
  e_I2->size[1] = loop_ub;
  e_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)e_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      e_I2->data[e_I2->size[0] * i1 + e_I2->size[0] * e_I2->size[1] * i0] =
        I2->data[(I2->size[0] * i1 + I2->size[0] * I2->size[1] * i0) + 2];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    b_I2[i0] = e_I2->size[i0];
  }

  emxFree_real_T(&e_I2);
  emxInit_real_T(&f_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv2, 3, b_I2, 3, &c_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[1];
  tmp_size_idx_0 = I2->size[2];
  i0 = f_I2->size[0] * f_I2->size[1] * f_I2->size[2];
  f_I2->size[0] = 1;
  f_I2->size[1] = i;
  f_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)f_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < i; i1++) {
      f_I2->data[f_I2->size[0] * i1 + f_I2->size[0] * f_I2->size[1] * i0] =
        I2->data[(I2->size[0] * i1 + I2->size[0] * I2->size[1] * i0) + 2];
    }
  }

  loop_ub = f_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = f_I2->size[1];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[I2->size[0] * d_tmp_data[i1] + I2->size[0] * I2->size[1] *
        e_tmp_data[i0]] = f_I2->data[f_I2->size[0] * i1 + f_I2->size[0] *
        f_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&f_I2);
  i0 = I2->size[0];
  i1 = I2->size[0];
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &e_emlrtBCI, emlrtRootTLSGlobal);
  loop_ub = I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  tmp_size_idx_0 = I2->size[2];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&g_I2, 3, &emlrtRTEI, TRUE);
  i0 = I2->size[0];
  i1 = (int32_T)((real_T)I2->size[0] - 2.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &d_emlrtBCI, emlrtRootTLSGlobal);
  iv3[0] = 1;
  iv3[1] = loop_ub;
  iv3[2] = tmp_size_idx_0;
  loop_ub = I2->size[1];
  tmp_size_idx_0 = I2->size[2];
  i = I2->size[0];
  i0 = g_I2->size[0] * g_I2->size[1] * g_I2->size[2];
  g_I2->size[0] = 1;
  g_I2->size[1] = loop_ub;
  g_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)g_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      g_I2->data[g_I2->size[0] * i1 + g_I2->size[0] * g_I2->size[1] * i0] =
        I2->data[((i + I2->size[0] * i1) + I2->size[0] * I2->size[1] * i0) - 3];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    b_I2[i0] = g_I2->size[i0];
  }

  emxFree_real_T(&g_I2);
  emxInit_real_T(&h_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv3, 3, b_I2, 3, &d_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[0] - 1;
  tmp_size_idx_0 = I2->size[0] - 3;
  b_tmp_size_idx_0 = I2->size[1];
  c_tmp_size_idx_0 = I2->size[2];
  i0 = h_I2->size[0] * h_I2->size[1] * h_I2->size[2];
  h_I2->size[0] = 1;
  h_I2->size[1] = b_tmp_size_idx_0;
  h_I2->size[2] = c_tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)h_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < c_tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < b_tmp_size_idx_0; i1++) {
      h_I2->data[h_I2->size[0] * i1 + h_I2->size[0] * h_I2->size[1] * i0] =
        I2->data[(tmp_size_idx_0 + I2->size[0] * i1) + I2->size[0] * I2->size[1]
        * i0];
    }
  }

  loop_ub = h_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = h_I2->size[1];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[(i + I2->size[0] * d_tmp_data[i1]) + I2->size[0] * I2->size[1] *
        e_tmp_data[i0]] = h_I2->data[h_I2->size[0] * i1 + h_I2->size[0] *
        h_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&h_I2);
  i0 = I2->size[0];
  i1 = (int32_T)((real_T)I2->size[0] - 1.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &g_emlrtBCI, emlrtRootTLSGlobal);
  loop_ub = I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  tmp_size_idx_0 = I2->size[2];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&i_I2, 3, &emlrtRTEI, TRUE);
  i0 = I2->size[0];
  i1 = (int32_T)((real_T)I2->size[0] - 2.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &f_emlrtBCI, emlrtRootTLSGlobal);
  iv4[0] = 1;
  iv4[1] = loop_ub;
  iv4[2] = tmp_size_idx_0;
  loop_ub = I2->size[1];
  tmp_size_idx_0 = I2->size[2];
  i = I2->size[0];
  i0 = i_I2->size[0] * i_I2->size[1] * i_I2->size[2];
  i_I2->size[0] = 1;
  i_I2->size[1] = loop_ub;
  i_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)i_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      i_I2->data[i_I2->size[0] * i1 + i_I2->size[0] * i_I2->size[1] * i0] =
        I2->data[((i + I2->size[0] * i1) + I2->size[0] * I2->size[1] * i0) - 3];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    b_I2[i0] = i_I2->size[i0];
  }

  emxFree_real_T(&i_I2);
  emxInit_real_T(&j_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv4, 3, b_I2, 3, &e_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[0] - 2;
  tmp_size_idx_0 = I2->size[0] - 3;
  b_tmp_size_idx_0 = I2->size[1];
  c_tmp_size_idx_0 = I2->size[2];
  i0 = j_I2->size[0] * j_I2->size[1] * j_I2->size[2];
  j_I2->size[0] = 1;
  j_I2->size[1] = b_tmp_size_idx_0;
  j_I2->size[2] = c_tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)j_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < c_tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < b_tmp_size_idx_0; i1++) {
      j_I2->data[j_I2->size[0] * i1 + j_I2->size[0] * j_I2->size[1] * i0] =
        I2->data[(tmp_size_idx_0 + I2->size[0] * i1) + I2->size[0] * I2->size[1]
        * i0];
    }
  }

  loop_ub = j_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = j_I2->size[1];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[(i + I2->size[0] * d_tmp_data[i1]) + I2->size[0] * I2->size[1] *
        e_tmp_data[i0]] = j_I2->data[j_I2->size[0] * i1 + j_I2->size[0] *
        j_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&j_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  tmp_size_idx_0 = I2->size[2];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&k_I2, 3, &emlrtRTEI, TRUE);
  iv5[0] = loop_ub;
  iv5[1] = 1;
  iv5[2] = tmp_size_idx_0;
  loop_ub = I2->size[0];
  tmp_size_idx_0 = I2->size[2];
  i0 = k_I2->size[0] * k_I2->size[1] * k_I2->size[2];
  k_I2->size[0] = loop_ub;
  k_I2->size[1] = 1;
  k_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)k_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      k_I2->data[i1 + k_I2->size[0] * k_I2->size[1] * i0] = I2->data[(i1 +
        (I2->size[0] << 1)) + I2->size[0] * I2->size[1] * i0];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    b_I2[i0] = k_I2->size[i0];
  }

  emxFree_real_T(&k_I2);
  emxInit_real_T(&l_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv5, 3, b_I2, 3, &f_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[0];
  tmp_size_idx_0 = I2->size[2];
  i0 = l_I2->size[0] * l_I2->size[1] * l_I2->size[2];
  l_I2->size[0] = i;
  l_I2->size[1] = 1;
  l_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)l_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < i; i1++) {
      l_I2->data[i1 + l_I2->size[0] * l_I2->size[1] * i0] = I2->data[(i1 +
        (I2->size[0] << 1)) + I2->size[0] * I2->size[1] * i0];
    }
  }

  loop_ub = l_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = l_I2->size[0];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[d_tmp_data[i1] + I2->size[0] * I2->size[1] * e_tmp_data[i0]] =
        l_I2->data[i1 + l_I2->size[0] * l_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&l_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  tmp_size_idx_0 = I2->size[2];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&m_I2, 3, &emlrtRTEI, TRUE);
  iv6[0] = loop_ub;
  iv6[1] = 1;
  iv6[2] = tmp_size_idx_0;
  loop_ub = I2->size[0];
  tmp_size_idx_0 = I2->size[2];
  i0 = m_I2->size[0] * m_I2->size[1] * m_I2->size[2];
  m_I2->size[0] = loop_ub;
  m_I2->size[1] = 1;
  m_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)m_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      m_I2->data[i1 + m_I2->size[0] * m_I2->size[1] * i0] = I2->data[(i1 +
        (I2->size[0] << 1)) + I2->size[0] * I2->size[1] * i0];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    b_I2[i0] = m_I2->size[i0];
  }

  emxFree_real_T(&m_I2);
  emxInit_real_T(&n_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv6, 3, b_I2, 3, &g_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[0];
  tmp_size_idx_0 = I2->size[2];
  i0 = n_I2->size[0] * n_I2->size[1] * n_I2->size[2];
  n_I2->size[0] = i;
  n_I2->size[1] = 1;
  n_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)n_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < i; i1++) {
      n_I2->data[i1 + n_I2->size[0] * n_I2->size[1] * i0] = I2->data[(i1 +
        (I2->size[0] << 1)) + I2->size[0] * I2->size[1] * i0];
    }
  }

  loop_ub = n_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = n_I2->size[0];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[(d_tmp_data[i1] + I2->size[0]) + I2->size[0] * I2->size[1] *
        e_tmp_data[i0]] = n_I2->data[i1 + n_I2->size[0] * n_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&n_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  i0 = I2->size[1];
  i1 = I2->size[1];
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &i_emlrtBCI, emlrtRootTLSGlobal);
  tmp_size_idx_0 = I2->size[2];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&o_I2, 3, &emlrtRTEI, TRUE);
  i0 = I2->size[1];
  i1 = (int32_T)((real_T)I2->size[1] - 2.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &h_emlrtBCI, emlrtRootTLSGlobal);
  iv7[0] = loop_ub;
  iv7[1] = 1;
  iv7[2] = tmp_size_idx_0;
  loop_ub = I2->size[0];
  tmp_size_idx_0 = I2->size[2];
  i = I2->size[1];
  i0 = o_I2->size[0] * o_I2->size[1] * o_I2->size[2];
  o_I2->size[0] = loop_ub;
  o_I2->size[1] = 1;
  o_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)o_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      o_I2->data[i1 + o_I2->size[0] * o_I2->size[1] * i0] = I2->data[(i1 +
        I2->size[0] * (i - 3)) + I2->size[0] * I2->size[1] * i0];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    b_I2[i0] = o_I2->size[i0];
  }

  emxFree_real_T(&o_I2);
  emxInit_real_T(&p_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv7, 3, b_I2, 3, &h_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[1] - 1;
  tmp_size_idx_0 = I2->size[0];
  b_tmp_size_idx_0 = I2->size[1] - 3;
  c_tmp_size_idx_0 = I2->size[2];
  i0 = p_I2->size[0] * p_I2->size[1] * p_I2->size[2];
  p_I2->size[0] = tmp_size_idx_0;
  p_I2->size[1] = 1;
  p_I2->size[2] = c_tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)p_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < c_tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      p_I2->data[i1 + p_I2->size[0] * p_I2->size[1] * i0] = I2->data[(i1 +
        I2->size[0] * b_tmp_size_idx_0) + I2->size[0] * I2->size[1] * i0];
    }
  }

  loop_ub = p_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = p_I2->size[0];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[(d_tmp_data[i1] + I2->size[0] * i) + I2->size[0] * I2->size[1] *
        e_tmp_data[i0]] = p_I2->data[i1 + p_I2->size[0] * p_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&p_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  i0 = I2->size[1];
  i1 = (int32_T)((real_T)I2->size[1] - 1.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &k_emlrtBCI, emlrtRootTLSGlobal);
  tmp_size_idx_0 = I2->size[2];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  emxInit_real_T(&q_I2, 3, &emlrtRTEI, TRUE);
  i0 = I2->size[1];
  i1 = (int32_T)((real_T)I2->size[1] - 2.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &j_emlrtBCI, emlrtRootTLSGlobal);
  iv8[0] = loop_ub;
  iv8[1] = 1;
  iv8[2] = tmp_size_idx_0;
  loop_ub = I2->size[0];
  tmp_size_idx_0 = I2->size[2];
  i = I2->size[1];
  i0 = q_I2->size[0] * q_I2->size[1] * q_I2->size[2];
  q_I2->size[0] = loop_ub;
  q_I2->size[1] = 1;
  q_I2->size[2] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)q_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      q_I2->data[i1 + q_I2->size[0] * q_I2->size[1] * i0] = I2->data[(i1 +
        I2->size[0] * (i - 3)) + I2->size[0] * I2->size[1] * i0];
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    b_I2[i0] = q_I2->size[i0];
  }

  emxFree_real_T(&q_I2);
  emxInit_real_T(&r_I2, 3, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv8, 3, b_I2, 3, &i_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[1] - 2;
  tmp_size_idx_0 = I2->size[0];
  b_tmp_size_idx_0 = I2->size[1] - 3;
  c_tmp_size_idx_0 = I2->size[2];
  i0 = r_I2->size[0] * r_I2->size[1] * r_I2->size[2];
  r_I2->size[0] = tmp_size_idx_0;
  r_I2->size[1] = 1;
  r_I2->size[2] = c_tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)r_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < c_tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      r_I2->data[i1 + r_I2->size[0] * r_I2->size[1] * i0] = I2->data[(i1 +
        I2->size[0] * b_tmp_size_idx_0) + I2->size[0] * I2->size[1] * i0];
    }
  }

  loop_ub = r_I2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = r_I2->size[0];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[(d_tmp_data[i1] + I2->size[0] * i) + I2->size[0] * I2->size[1] *
        e_tmp_data[i0]] = r_I2->data[i1 + r_I2->size[0] * r_I2->size[1] * i0];
    }
  }

  emxFree_real_T(&r_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  tmp_size_idx_0 = I2->size[1];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  b_emxInit_real_T(&s_I2, 2, &emlrtRTEI, TRUE);
  iv9[0] = loop_ub;
  iv9[1] = tmp_size_idx_0;
  loop_ub = I2->size[0];
  tmp_size_idx_0 = I2->size[1];
  i0 = s_I2->size[0] * s_I2->size[1];
  s_I2->size[0] = loop_ub;
  s_I2->size[1] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)s_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      s_I2->data[i1 + s_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        (I2->size[0] * I2->size[1] << 1)];
    }
  }

  for (i0 = 0; i0 < 2; i0++) {
    t_I2[i0] = s_I2->size[i0];
  }

  emxFree_real_T(&s_I2);
  b_emxInit_real_T(&u_I2, 2, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv9, 2, t_I2, 2, &j_emlrtECI, emlrtRootTLSGlobal);
  i = I2->size[0];
  tmp_size_idx_0 = I2->size[1];
  i0 = u_I2->size[0] * u_I2->size[1];
  u_I2->size[0] = i;
  u_I2->size[1] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)u_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < i; i1++) {
      u_I2->data[i1 + u_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        (I2->size[0] * I2->size[1] << 1)];
    }
  }

  loop_ub = u_I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = u_I2->size[0];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[d_tmp_data[i1] + I2->size[0] * e_tmp_data[i0]] = u_I2->data[i1 +
        u_I2->size[0] * i0];
    }
  }

  emxFree_real_T(&u_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  tmp_size_idx_0 = I2->size[1];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  b_emxInit_real_T(&v_I2, 2, &emlrtRTEI, TRUE);
  iv10[0] = loop_ub;
  iv10[1] = tmp_size_idx_0;
  loop_ub = I2->size[0];
  tmp_size_idx_0 = I2->size[1];
  i0 = v_I2->size[0] * v_I2->size[1];
  v_I2->size[0] = loop_ub;
  v_I2->size[1] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)v_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      v_I2->data[i1 + v_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        (I2->size[0] * I2->size[1] << 1)];
    }
  }

  for (i0 = 0; i0 < 2; i0++) {
    t_I2[i0] = v_I2->size[i0];
  }

  emxFree_real_T(&v_I2);
  b_emxInit_real_T(&w_I2, 2, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv10, 2, t_I2, 2, &k_emlrtECI,
    emlrtRootTLSGlobal);
  i = I2->size[0];
  tmp_size_idx_0 = I2->size[1];
  i0 = w_I2->size[0] * w_I2->size[1];
  w_I2->size[0] = i;
  w_I2->size[1] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)w_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < i; i1++) {
      w_I2->data[i1 + w_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        (I2->size[0] * I2->size[1] << 1)];
    }
  }

  loop_ub = w_I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = w_I2->size[0];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[(d_tmp_data[i1] + I2->size[0] * e_tmp_data[i0]) + I2->size[0] *
        I2->size[1]] = w_I2->data[i1 + w_I2->size[0] * i0];
    }
  }

  emxFree_real_T(&w_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  tmp_size_idx_0 = I2->size[1];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  b_emxInit_real_T(&x_I2, 2, &emlrtRTEI, TRUE);
  i0 = I2->size[2];
  i1 = I2->size[2];
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &m_emlrtBCI, emlrtRootTLSGlobal);
  i0 = I2->size[2];
  i1 = (int32_T)((real_T)I2->size[2] - 2.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &l_emlrtBCI, emlrtRootTLSGlobal);
  iv11[0] = loop_ub;
  iv11[1] = tmp_size_idx_0;
  loop_ub = I2->size[0];
  tmp_size_idx_0 = I2->size[1];
  i = I2->size[2];
  i0 = x_I2->size[0] * x_I2->size[1];
  x_I2->size[0] = loop_ub;
  x_I2->size[1] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)x_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      x_I2->data[i1 + x_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        I2->size[0] * I2->size[1] * (i - 3)];
    }
  }

  for (i0 = 0; i0 < 2; i0++) {
    t_I2[i0] = x_I2->size[i0];
  }

  emxFree_real_T(&x_I2);
  b_emxInit_real_T(&y_I2, 2, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv11, 2, t_I2, 2, &l_emlrtECI,
    emlrtRootTLSGlobal);
  i = I2->size[2] - 1;
  tmp_size_idx_0 = I2->size[0];
  b_tmp_size_idx_0 = I2->size[1];
  c_tmp_size_idx_0 = I2->size[2] - 3;
  i0 = y_I2->size[0] * y_I2->size[1];
  y_I2->size[0] = tmp_size_idx_0;
  y_I2->size[1] = b_tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)y_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      y_I2->data[i1 + y_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        I2->size[0] * I2->size[1] * c_tmp_size_idx_0];
    }
  }

  loop_ub = y_I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = y_I2->size[0];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[(d_tmp_data[i1] + I2->size[0] * e_tmp_data[i0]) + I2->size[0] *
        I2->size[1] * i] = y_I2->data[i1 + y_I2->size[0] * i0];
    }
  }

  emxFree_real_T(&y_I2);
  loop_ub = I2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_tmp_data[i0] = i0;
  }

  tmp_size_idx_0 = I2->size[1];
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    e_tmp_data[i0] = i0;
  }

  b_emxInit_real_T(&ab_I2, 2, &emlrtRTEI, TRUE);
  i0 = I2->size[2];
  i1 = (int32_T)((real_T)I2->size[2] - 1.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &o_emlrtBCI, emlrtRootTLSGlobal);
  i0 = I2->size[2];
  i1 = (int32_T)((real_T)I2->size[2] - 2.0);
  emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &n_emlrtBCI, emlrtRootTLSGlobal);
  iv12[0] = loop_ub;
  iv12[1] = tmp_size_idx_0;
  loop_ub = I2->size[0];
  tmp_size_idx_0 = I2->size[1];
  i = I2->size[2];
  i0 = ab_I2->size[0] * ab_I2->size[1];
  ab_I2->size[0] = loop_ub;
  ab_I2->size[1] = tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)ab_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      ab_I2->data[i1 + ab_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        I2->size[0] * I2->size[1] * (i - 3)];
    }
  }

  for (i0 = 0; i0 < 2; i0++) {
    t_I2[i0] = ab_I2->size[i0];
  }

  emxFree_real_T(&ab_I2);
  b_emxInit_real_T(&bb_I2, 2, &emlrtRTEI, TRUE);
  emlrtSubAssignSizeCheckR2012b(iv12, 2, t_I2, 2, &m_emlrtECI,
    emlrtRootTLSGlobal);
  i = I2->size[2] - 2;
  tmp_size_idx_0 = I2->size[0];
  b_tmp_size_idx_0 = I2->size[1];
  c_tmp_size_idx_0 = I2->size[2] - 3;
  i0 = bb_I2->size[0] * bb_I2->size[1];
  bb_I2->size[0] = tmp_size_idx_0;
  bb_I2->size[1] = b_tmp_size_idx_0;
  emxEnsureCapacity((emxArray__common *)bb_I2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < b_tmp_size_idx_0; i0++) {
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      bb_I2->data[i1 + bb_I2->size[0] * i0] = I2->data[(i1 + I2->size[0] * i0) +
        I2->size[0] * I2->size[1] * c_tmp_size_idx_0];
    }
  }

  loop_ub = bb_I2->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    tmp_size_idx_0 = bb_I2->size[0];
    for (i1 = 0; i1 < tmp_size_idx_0; i1++) {
      I2->data[(d_tmp_data[i1] + I2->size[0] * e_tmp_data[i0]) + I2->size[0] *
        I2->size[1] * i] = bb_I2->data[i1 + bb_I2->size[0] * i0];
    }
  }

  emxFree_real_T(&bb_I2);

  /* clear I2; */
  i = 0;
  while (i <= (int32_T)((newdim[0] - 1.0) + 1.0) - 1) {
    b_tmp_size_idx_0 = 0;
    while (b_tmp_size_idx_0 <= (int32_T)((newdim[1] - 1.0) + 1.0) - 1) {
      c_tmp_size_idx_0 = 0;
      while (c_tmp_size_idx_0 <= (int32_T)((newdim[2] - 1.0) + 1.0) - 1) {
        y = 2.0 * (real_T)i;
        b_y = 2.0 * (real_T)b_tmp_size_idx_0;
        c_y = 2.0 * (real_T)c_tmp_size_idx_0;
        for (i0 = 0; i0 < 5; i0++) {
          for (i1 = 0; i1 < 5; i1++) {
            for (i2 = 0; i2 < 5; i2++) {
              i3 = I2->size[0];
              i4 = (int32_T)((y + (-2.0 + (real_T)i2)) + 3.0);
              i5 = I2->size[1];
              tmp_size_idx_0 = (int32_T)((b_y + (-2.0 + (real_T)i1)) + 3.0);
              loop_ub = I2->size[2];
              i6 = (int32_T)((c_y + (-2.0 + (real_T)i0)) + 3.0);
              A[(i2 + 5 * i1) + 25 * i0] = I2->data
                [((emlrtDynamicBoundsCheckFastR2012b(i4, 1, i3, &p_emlrtBCI,
                    emlrtRootTLSGlobal) + I2->size[0] *
                   (emlrtDynamicBoundsCheckFastR2012b(tmp_size_idx_0, 1, i5,
                     &q_emlrtBCI, emlrtRootTLSGlobal) - 1)) + I2->size[0] *
                  I2->size[1] * (emlrtDynamicBoundsCheckFastR2012b(i6, 1,
                    loop_ub, &r_emlrtBCI, emlrtRootTLSGlobal) - 1)) - 1] * Wt3
                [(i2 + 5 * i1) + 25 * i0];
            }
          }
        }

        i0 = IResult->size[0];
        i1 = (int32_T)((real_T)i + 1.0);
        i2 = IResult->size[1];
        i3 = (int32_T)((real_T)b_tmp_size_idx_0 + 1.0);
        i4 = IResult->size[2];
        i5 = (int32_T)((real_T)c_tmp_size_idx_0 + 1.0);
        IResult->data[((emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &s_emlrtBCI,
          emlrtRootTLSGlobal) + IResult->size[0] *
                        (emlrtDynamicBoundsCheckFastR2012b(i3, 1, i2,
          &t_emlrtBCI, emlrtRootTLSGlobal) - 1)) + IResult->size[0] *
                       IResult->size[1] * (emlrtDynamicBoundsCheckFastR2012b(i5,
          1, i4, &u_emlrtBCI, emlrtRootTLSGlobal) - 1)) - 1] = sum(A);
        c_tmp_size_idx_0++;
        emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar,
          emlrtRootTLSGlobal);
      }

      b_tmp_size_idx_0++;
      emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
    }

    i++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, emlrtRootTLSGlobal);
  }

  emxFree_real_T(&I2);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (GPReduce.c) */
