dnl **********************************************************************
dnl * Generates alpha, A, x, beta, and y, where x two parts: 
dnl * (head_x, tail_x) and computes r_true.                              *
dnl **********************************************************************
dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
define(`SYMV2_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * uplo         (input) enum blas_uplo_type
 *              Which half of the symmetric matrix A is to be stored.
 *
 * n            (input) int
 *              Size of matrix A, length of vectors x and y.
 *
 * alpha        (input/output) $1_array
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) $2_array
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) $3_array
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
dnl
dnl ---------------------------------------------------------------------
dnl Usage: SYMV2_TESTGEN(aby_typeltr, A_typeltr, x_typeltr)
dnl        produce symv2_prepare signature
dnl ---------------------------------------------------------------------
define(`SYMV2_TESTGEN_NAME', 
  `BLAS_$1symv2`'ifelse(`$2&&$3', `$1&&$1', `', `_$2_$3')_testgen')dnl
dnl
dnl
define(`SYMV2_TESTGEN_HEAD', 
  `void SYMV2_TESTGEN_NAME($1, $2, $3)(int norm, dnl
       enum blas_order_type order, enum blas_uplo_type uplo, int n, dnl
       $1_array alpha, int alpha_flag, $2_array A, dnl
       int lda, $3_array head_x, $3_array tail_x, $1_array beta, dnl
       int beta_flag, $1_array y, int *seed, dnl
       double *r_true_l, double *r_true_t)')dnl
dnl
dnl
define(`SYMV2_TESTGEN', 
  `SYMV2_TESTGEN_HEAD($1, $2, $3)
   SYMV2_TESTGEN_COMMENT($1, $2, $3)
   SYMV2_TESTGEN_BODY($1, $2, $3)')dnl
dnl
dnl
dnl ---------------------------------------------------------------------
dnl Usage: SYMV2_TESTGEN_BODY(aby_typeltr, A_typeltr, x_typeltr)
dnl        produce symv2_prepare signature
dnl ---------------------------------------------------------------------
define(`SYMV2_TESTGEN_BODY',
`{
  PTR_CAST(y, $1_type)
  int n_fix2;
  int n_mix;
  int i;
  DECLARE_VECTOR(temp, $2_type)
  int incy, incA;
  DECLARE(y_elem, $1_type)

  incy = incA = 1;
  INC_ADJUST(incy, $1_type)
  INC_ADJUST(incA, $2_type)
  MALLOC_VECTOR(temp, $2_type, n*incA)

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */  
  n_mix = 0;
  for (i = 0; i < n; i++) {
    n_fix2 = i;
    if (i >= 1) {
      n_mix = n - n_fix2;
      alpha_flag = 1; beta_flag = 1;
    }

    $2sy_copy_row(order, uplo, n, A, lda, temp, i);
    DOT2_TESTGEN_NAME($1, $3, $2)(n, n_fix2, n_mix, norm, blas_no_conj, dnl
        alpha, alpha_flag, beta, beta_flag, head_x, tail_x, temp, seed, dnl
        PASS_BY_REF(y_elem, $1_type), &r_true_l[i*incy], &r_true_t[i*incy]);
    SET_VECTOR_ELEMENT(y_i, i*incy, y_elem, $1_type)

    /* copy temp to A */
    $2sy_commit_row(order, uplo, n, A, lda, temp, i);
  }
  FREE_VECTOR(temp, $2_type)
}')dnl
dnl
dnl
define(`PROTOTYPES', `FOREACH(`PREC_ARGS', `SYMV2_TESTGEN_HEAD(arg);
')')dnl
dnl
dnl
define(`SOURCE', `
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

FOREACH(`PREC_ARGS', `SYMV2_TESTGEN(arg)
')')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
