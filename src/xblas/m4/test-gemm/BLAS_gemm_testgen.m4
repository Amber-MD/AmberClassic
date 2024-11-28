dnl
dnl BLAS_gemm_testgen.m4
dnl
dnl Test case generator for gemm routines.
dnl Generates test cases for alpha, A, beta, B, and C, and
dnl computes r_true in double-double precision.
dnl
dnl
dnl C Interface has the form
dnl
dnl    BLAS_xgemm{_a_b}{_x}_testgen(int norm,
dnl                              enum blas_order_type order,
dnl                              enum blas_trans_type transa,
dnl                              enum blas_trans_type transb,
dnl                              int m, int n, int k, 
dnl                              int randomize, 
dnl                              SCALAR *alpha, int alpha_flag, 
dnl                              ARRAY a, int lda, 
dnl                              SCALAR *beta, int beta_flag, 
dnl                              ARRAY b, int ldb, 
dnl                              ARRAY c, int ldc, 
dnl                              int seed, 
dnl                              double *HEAD(r_true), 
dnl                              double *TAIL(r_true))
dnl                              
dnl Arguments
dnl   norm    (in) blas_norm_type
dnl
dnl   order   (in) blas_order_type
dnl                  determines the storage format for matrices.
dnl
dnl   transa  (in) blas_trans_type
dnl                  whether matrix A is to be transposed.
dnl   transb  (in) blas_trans_type
dnl                  whether matrix B is to be transposed.
dnl
dnl   m, n, k (in) int
dnl                  the dimensions of matrices
dnl                  matrix C is m-by-n
dnl                  matrix A is m-by-k (after transpose, if necessary)
dnl                  matrix B is k-by-n (after transpose, if necessary)
dnl
dnl   randomize (in) int
dnl                  if 0, entries in matrices A, B will be chosen for
dnl                        maximum cancellation, but with less randomness.
dnl                  if 1, every entry in the matrix A, B will be 
dnl                        random.
dnl
dnl   alpha   (in/out) SCALAR
dnl                  if alpha_flag = 1, alpha is input
dnl                  if alpha_flag = 0, alpha is output
dnl   alpha_flag (in) int
dnl                  see above
dnl
dnl   a, b, c (out)  matrix a, b, c.
dnl   lda, ldb, ldc (in) leading dimensions of matrices above.
dnl
dnl   beta    (in/out) SCALAR
dnl                  if beta_flag = 1, beta is input
dnl                  if beta_flag = 0, beta is output
dnl   beta_flag  (in) int
dnl                  see above
dnl
dnl   seed    (in/out) int *
dnl
dnl   HEAD(r_true)   (out) double*
dnl                 the leading part of the truth double-double
dnl   TAIL(r_true)   (out) double*
dnl                 the trailing part of the truth double-double
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
dnl  Usage: GEMM_TESTGEN_NAME($1, $2, $3)
dnl         GEMM_TESTGEN_PARAMS($1, $2, $3)
dnl         GEMM_TESTGEN_HEAD($1, $2, $3)
dnl    $1 -- type for matrix C, scalars alpha, beta
dnl    $2 -- type for matrix A
dnl    $3 -- type for matrix B
dnl
define(`GEMM_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1gemm_testgen', `BLAS_$1gemm_$2_$3_testgen')')dnl
dnl
define(`GEMM_TESTGEN_PARAMS',
       `int norm, enum blas_order_type order,
        enum blas_trans_type transa, enum blas_trans_type transb, 
        int m, int n, int k, int randomize, 
        $1_array alpha, int alpha_flag, $2_array a, int lda, 
        $1_array beta, int beta_flag, $3_array b, int ldb, 
        $1_array c, int ldc, int *seed, double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
define(`GEMM_TESTGEN_HEAD',
  `void GEMM_TESTGEN_NAME($1, $2, $3)(GEMM_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
dnl  Usage: GEMM_TESTGEN_BODY($1, $2, $3)
dnl    $1 -- type for matrix C, scalars alpha, beta
dnl    $2 -- type for matrix A.
dnl    $3 -- type for matrix B.
dnl
define(`GEMM_TESTGEN_BODY', `{
  pushdef(`__IS_MIXED_PREC', `IS_MIXED_PREC($1_type, $2_type, $3_type)')dnl
  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  DECLARE(c_elem, $1_type)
  DECLARE(a_elem, $2_type)
  DECLARE(b_elem, $3_type)

  DECLARE(r_true_elem, EXTRA_TYPE($1_type))

  PTR_CAST(c, $1_type)
  PTR_CAST(b, $3_type)
  PTR_CAST(a, $2_type)
  PTR_CAST(alpha, $1_type)
  PTR_CAST(beta, $1_type)

  /* Temporary storage space for vectors of length k */
  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(b_vec, $3_type)

  inca = incb = 1;
  INC_ADJUST(inca, $2_type)  
  INC_ADJUST(incb, $3_type)  

  MALLOC_VECTOR(a_vec, $2_type, k)
  for (i = 0; i < k*inca; i += inca) {
        SET_ZERO_VECTOR_ELEMENT(a_vec, i, $2_type)
  }
  MALLOC_VECTOR(b_vec, $3_type, k)
  for (i = 0; i < k*incb; i += incb) {
        SET_ZERO_VECTOR_ELEMENT(b_vec, i, $3_type)
  }

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  INC_ADJUST(incci, $1_type)
  INC_ADJUST(inccij, $1_type)

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    DOT_TESTGEN_NAME($1, $3, $2)(k, 0, 0, norm, blas_no_conj, 
                                 alpha, alpha_flag, beta, beta_flag, 
                                 b_vec, a_vec, seed, PASS_BY_REF(c_elem, $1_type), 
                                 PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));
    SET_VECTOR_ELEMENT(c_i, cij, c_elem, $1_type)
    SET_VECTOR_ELEMENT(r_true, cij, r_true_elem, EXTRA_TYPE($1_type))

    /* Copy a_vec to the first row of A */
    $2ge_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      $3ge_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      DOT_TESTGEN_NAME($1, $3, $2)(k, 0, k, norm, blas_no_conj, alpha, 1, 
                                   beta, 1, b_vec, a_vec, seed, 
                                   PASS_BY_REF(c_elem, $1_type), 
                                   PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));

      SET_VECTOR_ELEMENT(c_i, cij, c_elem, $1_type)
      SET_VECTOR_ELEMENT(r_true, cij, r_true_elem, EXTRA_TYPE($1_type));
      $2ge_commit_row(order, transa, m, k, a, lda, a_vec, i);    
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
        GET_VECTOR_ELEMENT(c_elem, c_i, ci, $1_type)
        SET_VECTOR_ELEMENT(c_i, cij, c_elem, $1_type)
        GET_VECTOR_ELEMENT(r_true_elem, r_true, ci, EXTRA_TYPE($1_type))
        SET_VECTOR_ELEMENT(r_true, cij, r_true_elem, EXTRA_TYPE($1_type))
      }
    }
  } else {

    ifelse(IS_MIXED($1, $2), `t', 
      `DECLARE_VECTOR(aa_vec, COMPLEX_TYPE($2_type))')
    ifelse(IS_MIXED($1, $3), `t', 
      `DECLARE_VECTOR(bb_vec, COMPLEX_TYPE($3_type))')

    ifelse(IS_MIXED($1, $2), `t', 
      `MALLOC_VECTOR(aa_vec, COMPLEX_TYPE($2_type), k)')
    ifelse(IS_MIXED($1, $3), `t', 
      `MALLOC_VECTOR(bb_vec, COMPLEX_TYPE($3_type), k)')


    if (alpha_flag == 0) {
      RANDOM(c_elem, $1_type, __IS_MIXED_PREC)
      SET_VECTOR_ELEMENT(alpha_i, 0, c_elem, $1_type)
    }
    if (beta_flag == 0) {
      RANDOM(c_elem, $1_type, __IS_MIXED_PREC)
      SET_VECTOR_ELEMENT(beta_i, 0, c_elem, $1_type)
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
        (order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
        (order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }
    
    INC_ADJUST(incbi, $3_type)
    INC_ADJUST(incbij, $3_type)
    INC_ADJUST(incai, $2_type)
    INC_ADJUST(incaij, $2_type)

    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
        RANDOM(b_elem, $3_type, __IS_MIXED_PREC)
        SET_VECTOR_ELEMENT(b_i, bij, b_elem, $3_type)
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
        RANDOM(a_elem, $2_type, __IS_MIXED_PREC)
        SET_VECTOR_ELEMENT(a_i, aij, a_elem, $2_type)
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      $2ge_copy_row(order, transa, m, k, a, lda, a_vec, i);

      ifelse(IS_MIXED($1, $2), `t', 
        `{
         int r;
         for (r = 0; r < k; r++) {
           aa_vec[2*r] = a_vec[r];
           aa_vec[2*r+1] = 0.0;
         }
        }')

      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
        $3ge_copy_col(order, transb, k, n, b, ldb, b_vec, j);

        ifelse(IS_MIXED($1, $3), `t', 
          `{
           int r;
           for (r = 0; r < k; r++) {
             bb_vec[2*r] = b_vec[r];
             bb_vec[2*r+1] = 0.0;
           }
          }')

        ifelse(IS_MIXED($1, $2, $3), `t', 
        `DOT_TESTGEN_NAME($1, $1, $1)(k, k, 0, norm, blas_no_conj, alpha, 1, 
                         beta, 1, 
                         ifelse(IS_MIXED($1, $3), `t', `bb_vec', `b_vec'), 
                         ifelse(IS_MIXED($1, $2), `t', `aa_vec', `a_vec'), 
                         seed, 
                         PASS_BY_REF(c_elem, $1_type), 
                         PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));',
        `DOT_TESTGEN_NAME($1, $3, $2)(k, k, 0, norm, blas_no_conj, alpha, 1, 
                                     beta, 1, b_vec, a_vec, seed, 
                                     PASS_BY_REF(c_elem, $1_type), 
                                     PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));')

        SET_VECTOR_ELEMENT(c_i, cij, c_elem, $1_type)
        SET_VECTOR_ELEMENT(r_true, cij, r_true_elem, EXTRA_TYPE($1_type))
      }
    }

    ifelse(IS_MIXED($1, $2), `t', 
      `FREE_VECTOR(aa_vec, COMPLEX_TYPE($2_type))')
    ifelse(IS_MIXED($1, $3), `t', 
      `FREE_VECTOR(bb_vec, COMPLEX_TYPE($3_type))')
  }

  FREE_VECTOR(a_vec, $2_type)
  FREE_VECTOR(b_vec, $3_type)
  popdef(`__IS_MIXED_PREC')dnl
}')dnl
dnl
dnl
dnl
dnl  Usage: GEMM_TESTGEN_NAME($1, $2, $3)
dnl         GEMM_TESTGEN_PARAMS($1, $2, $3)
dnl         GEMM_TESTGEN_HEAD($1, $2, $3)
dnl    $1 -- type for matrix C, scalars alpha, beta
dnl    $2 -- type for matrix A
dnl    $3 -- type for matrix B
dnl
define(`GEMM_TESTGEN', 
  `GEMM_TESTGEN_HEAD($1, $2, $3)
   GEMM_TESTGEN_BODY($1, $2, $3)')dnl
dnl
dnl
define(`PROTOTYPES', `FOREACH(`PREC_ARGS', `
GEMM_TESTGEN_HEAD(arg);')')dnl
dnl
dnl
define(`SOURCE', `
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

FOREACH(`PREC_ARGS', `GEMM_TESTGEN(arg)
')')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
