dnl **********************************************************************
dnl * Generates alpha, AB, x, beta, and y, where AB is a symmetric       *
dnl * packed matrix; and computes r_true.                                *
dnl **********************************************************************
dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
define(`GBMV_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) $1_array
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) $2_array
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) $3_array
 *
 * beta         (input/output) $1_array
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) $1_array
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */')dnl
dnl
dnl
define(`GBMV_TESTGEN_NAME', 
  `BLAS_$1gbmv`'ifelse(`$2&&$3', `$1&&$1', `', `_$2_$3')_testgen')dnl
dnl
dnl
define(`GBMV_TESTGEN_HEAD', 
  `void GBMV_TESTGEN_NAME($1, $2, $3)(int norm, enum blas_order_type order, dnl
       enum blas_trans_type trans, int m, int n, int kl, int ku, dnl
       $1_array alpha, int alpha_flag, $2_array AB, int lda, $3_array x, dnl
       $1_array beta, int beta_flag, $1_array y, int *seed, dnl
       double *r_true_l, double *r_true_t)')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: GBMV_TESTGEN(aby_typeltr, AB_typeltr, x_typeltr)
dnl        produce gbmv_prepare signature
dnl ---------------------------------------------------------------------
define(`GBMV_TESTGEN', 
  `GBMV_TESTGEN_HEAD($1, $2, $3)
   GBMV_TESTGEN_COMMENT($1, $2, $3)
   GBMV_TESTGEN_BODY($1, $2, $3)')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: GBMV_TESTGEN_BODY(aby_typeltr, AB_typeltr, x_typeltr)
dnl        produce gbmv_prepare signature
dnl ---------------------------------------------------------------------
define(`GBMV_TESTGEN_BODY',
`{
  PTR_CAST(x, $3_type)
  PTR_CAST(y, $1_type)
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  IF_COMPLEX($2_type, `int j;')
  DECLARE_VECTOR(temp, $2_type)
  int m_i, n_i;
  int max_mn;
  int incy, incAB, incx;
  DECLARE(y_elem, $1_type)

  max_mn = MAX(m, n);
  incx = incy  = incAB = 1;
  INC_ADJUST(incy, $1_type)
  INC_ADJUST(incAB, $2_type)
  INC_ADJUST(incx, $3_type)

  if (trans==blas_no_trans){
    m_i=m; n_i=n;
  }
  else{
    m_i=n; n_i=m;
  }

  MALLOC_VECTOR(temp, $2_type, max_mn)

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */  
  for(i=0; i<m_i; i++){
    /* copy AB to temp */
    $2gbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i, 
                   &n_fix2, &n_mix, &ysize);

    if (i == 1){
      /* from now on, fix alpha and beta */
      alpha_flag = 1; beta_flag = 1;
    }
  
    DOT_TESTGEN_NAME($1, $3, $2)(ysize, n_fix2, n_mix, norm, dnl
        blas_no_conj, alpha, alpha_flag, beta, beta_flag, dnl
        x, temp, seed, PASS_BY_REF(y_elem, $1_type), 
         &r_true_l[i*incy], &r_true_t[i*incy]);
    SET_VECTOR_ELEMENT(y_i, i*incy, y_elem, $1_type)

    IF_COMPLEX($2_type, 
      `if (trans == blas_conj_trans){
         for(j=0; j<n_i*incAB; j += 2){       
           temp[j+1] = -temp[j+1];
         }
       }')

    /* copy temp to AB */
    $2gbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  /* Zero out trailing part of x */
  for (i = ysize; i < n_i; i++) {
    SET_ZERO_VECTOR_ELEMENT(x_i, i*incx, $3_type)
  }

  FREE_VECTOR(temp, $2_type)
}')dnl
dnl
dnl
define(`PROTOTYPES', `FOREACH(`PREC_ARGS', `
GBMV_TESTGEN_HEAD(arg);')')dnl
dnl
dnl
define(`SOURCE', `dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

FOREACH(`PREC_ARGS', `
GBMV_TESTGEN(arg)')
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
