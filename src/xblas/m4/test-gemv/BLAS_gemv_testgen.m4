dnl **********************************************************************
dnl * Generates alpha, A, x, beta, and y, where A is a symmetric matrix  *
dnl * and computes r_true.                                               *
dnl **********************************************************************
dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
define(`GEMV_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) $1_array
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A           (output) $2_array
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
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
 */')
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: GEMV_TESTGEN(aby_typeltr, A_typeltr, x_typeltr)
dnl        produce gemv_prepare signature
dnl ---------------------------------------------------------------------
define(`GEMV_TESTGEN_NAME', 
  `BLAS_$1gemv`'ifelse(`$2&&$3', `$1&&$1', `', `_$2_$3')_testgen')dnl
dnl
dnl
define(`GEMV_TESTGEN_HEAD', 
  `void GEMV_TESTGEN_NAME($1, $2, $3)(int norm, enum blas_order_type order, dnl
      enum blas_trans_type trans, int m, int n, $1_array alpha, dnl
      int alpha_flag, $2_array A, int lda, $3_array x, $1_array beta, dnl
      int beta_flag, $1_array y, int *seed, double *r_true_l, dnl
      double *r_true_t)')dnl
dnl
dnl
define(`GEMV_TESTGEN', 
  `GEMV_TESTGEN_HEAD($1, $2, $3)
   GEMV_TESTGEN_COMMENT($1, $2, $3)
   GEMV_TESTGEN_BODY($1, $2, $3)')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: GEMV_TESTGEN_BODY(aby_typeltr, A_typeltr, x_typeltr)
dnl        produce gemv_prepare signature
dnl ---------------------------------------------------------------------
define(`GEMV_TESTGEN_BODY',
`{
  PTR_CAST(y, $1_type)
  int n_fix2;
  int n_mix;
  int i;
  DECLARE_VECTOR(temp, $2_type)
  int m_i, n_i;
  int max_mn;
  int incy, incA;
  DECLARE(y_elem, $1_type)

  incy = incA = 1;
  INC_ADJUST(incy, $1_type)
  INC_ADJUST(incA, $2_type)

  max_mn = MAX(m, n);
  
  if (trans==blas_no_trans) {
    m_i=m; n_i=n;
  } else {
    m_i=n; n_i=m;
  }

  MALLOC_VECTOR(temp, $2_type, max_mn)

  /* calling dot_testgen n time. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */  
  n_fix2 = n_mix = 0;
  for(i=0; i<m_i; i++) {

    if (i == 0) {
      n_fix2 = 0;
      n_mix  = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix  = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1; beta_flag = 1;
    }
        
    DOT_TESTGEN_NAME($1, $3, $2)(n_i, n_fix2, n_mix, norm, blas_no_conj, dnl
        alpha, alpha_flag, beta, beta_flag, x, temp, seed, dnl
        PASS_BY_REF(y_elem, $1_type), &r_true_l[i*incy], &r_true_t[i*incy]);
    SET_VECTOR_ELEMENT(y_i, i*incy, y_elem, $1_type)

    /* copy temp to A */
    $2ge_commit_row(order, trans, m_i, n_i, A, lda, temp, i);
  }  

  FREE_VECTOR(temp, $2_type)
}')dnl
dnl
dnl
define(`PROTOTYPES', `FOREACH(`PREC_ARGS', `
GEMV_TESTGEN_HEAD(arg);')')dnl
dnl
dnl
define(`SOURCE', `dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

FOREACH(`PREC_ARGS', `
GEMV_TESTGEN(arg)')
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
