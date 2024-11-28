dnl
dnl symm_testgen.m4
dnl
dnl Test case generator for symm routines.
dnl Generates test cases for alpha, A, beta, B and C and
dnl computes r_true in double-double precision.
dnl
dnl
dnl C Interface has the form 
dnl
dnl   BLAS_xsymm{_a_b}{_x}_testgen(int norm, 
dnl                             enum blas_order_type order, 
dnl                             enum blas_uplo_type uplo, 
dnl                             enum blas_side_type side, 
dnl                             int m, int n, int randomize, 
dnl                             SCALAR *alpha, int alpha_flag, 
dnl                             SCALAR *beta, int beta_flag, 
dnl                             ARRAY a, int lda, 
dnl                             ARRAY b, int ldb, 
dnl                             ARRAY c, int ldc, 
dnl                             int *seed, 
dnl                             double *HEAD(r_true), double *TAIL(r_true))
dnl
dnl Arguments
dnl   norm      (in) int
dnl
dnl   order     (in) blas_order_type
dnl                    determines the storage format for matrices
dnl   
dnl   uplo      (in) blas_uplo_type
dnl                    determines whether the upper triangular portion
dnl                    or the lower triangular portion of the symmetric
dnl                    matrix A is used.
dnl 
dnl   side      (in) blas_side_type
dnl                    determines on which side the symmetric matrix A
dnl                    is multiplied.
dnl
dnl   m, n      (in) int
dnl                    the dimensions of matrices
dnl                    matrix B and C are m-by-n
dnl                    matrix A is m-by-m if side = left, 
dnl                                n-by-n if side = right.
dnl   
dnl   randomize (in) int
dnl                  if 0, entries in matrices A, B will be chosen for
dnl                        maximum cancellation, but with less randomness.
dnl                  if 1, every entry in the matrix A, B will be 
dnl                        random.
dnl
dnl   alpha     (in/out) SCALAR
dnl                    if alpha_flag = 1, alpha is input
dnl                    if alpha_flag = 0, alpha is output
dnl
dnl   alpha_flag (in) int
dnl                    see above
dnl
dnl   beta      (in/out) SCALAR
dnl                    if beta_flag = 1, beta is input
dnl                    if beta_flag = 0, beta is output
dnl
dnl   beta_flag (in) int
dnl                    see above
dnl
dnl   a, b, c (out) matrices a, b, c.
dnl   lda, ldb, ldc (in) leading dimensions of matrices above
dnl
dnl   seed      (in/out) int
dnl
dnl   HEAD(r_true)  (out) double *
dnl   TAIL(r_true)  (out) double *
dnl               the leading/trailing part of the true in double-double
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
define(`SYMM_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1symm$4', 
  `BLAS_$1symm_$2_$3$4')')dnl
dnl
dnl
dnl  SYMM_TESTGEN
dnl     |
dnl     |-- SYMM_TESTGEN_HEAD
dnl     |      |
dnl     |      |-- SYMM_TESTGEN_NAME
dnl     |      |
dnl     |      |-- SYMM_TESTGEN_PARAMS
dnl     |      
dnl     |-- SYMM_TESTGEN_COMMENT
dnl     |
dnl     |-- SYMM_TESTGEN_BODY
dnl
dnl  Usage:
dnl    SYMM_TESTGEN        ($1, $2, $3)
dnl    SYMM_TESTGEN_HEAD   ($1, $2, $3)
dnl    SYMM_TESTGEN_NAME   ($1, $2, $3)
dnl    SYMM_TESTGEN_PARAMS ($1, $2, $3)
dnl    SYMM_TESTGEN_COMMENT($1, $2, $3)
dnl    SYMM_TESTGEN_BODY   ($1, $2, $3)
dnl
dnl    $1 -- type of alpha, beta, c
dnl    $2 -- type of a
dnl    $3 -- type of b
dnl
define(`SYMM_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1symm_testgen', `BLAS_$1symm_$2_$3_testgen')')dnl
dnl
dnl
define(`SYMM_TESTGEN_PARAMS',
  `int norm, enum blas_order_type order, dnl
   enum blas_uplo_type uplo, enum blas_side_type side, dnl
   int m, int n,  int randomize, dnl
   $1_array alpha, int alpha_flag, $1_array beta, int beta_flag, dnl
   $2_array a, int lda, $3_array b, int ldb, $1_array c, int ldc, dnl
   int *seed, double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
dnl
define(`SYMM_TESTGEN_HEAD', 
  `void SYMM_TESTGEN_NAME($1, $2, $3)(SYMM_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`SYMM_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to SYMM_NAME($1, $2, $3, `'){_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_side_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
 *
 * side    (input) enum blas_side_type
 *           which side of matrix b matrix a is to be multiplied.
 *
 * m n     (input) int
 *           sizes of matrices a, b, c:
 *              matrix a is m-by-m for left multiplication
 *                          n-by-n otherwise, 
 *              matrices b, c are m-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, B will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, B will be 
 *              random.
 *
 * alpha   (input/output) $1_array
 *           if alpha_flag = 1, alpha is input.
 *           if alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *           = 0: alpha is free, and is output.
 *           = 1: alpha is fixed on input.
 *
 * beta    (input/output) $1_array
 *           if beta_flag = 1, beta is input.
 *           if beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *           = 0: beta is free, and is output.
 *           = 1: beta is fixed on input.
 * 
 * a       (input/output) $2_array
 * 
 * lda     (input) lda
 *         leading dimension of matrix A.
 *
 * b       (input/output) $3_array
 *
 * ldb     (input) int
 *         leading dimension of matrix B.
 * 
 * c       (input/output) $1_array
 *         generated matrix C that will be used as an input to SYMM.
 * 
 * ldc     (input) int
 *         leading dimension of matrix C.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * double  (output) *HEAD(r_true)
 *         the leading part of the truth in double-double.
 *
 * double  (output) *TAIL(r_true)
 *         the trailing part of the truth in double-double
 *
 */')dnl
dnl
dnl
define(`SYMM_TESTGEN_BODY', `{

  int i, j;
  int cij, ci;
  int bij, bi;
  int aij, ai;
  int inccij, incci;
  int incbij, incbi;
  int incaij, incai;
  int inca, incb;
  int m_i, n_i;

  DECLARE(c_elem, $1_type)
  DECLARE(a_elem, $2_type)
  DECLARE(b_elem, $3_type)
  DECLARE(r_true_elem, EXTRA_TYPE($1_type))

  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(b_vec, $3_type)

  PTR_CAST(c, $1_type)
  PTR_CAST(alpha, $1_type)
  PTR_CAST(beta, $1_type)
  PTR_CAST(a, $2_type)
  PTR_CAST(b, $3_type)

  if (side == blas_left_side) {
    m_i = m; n_i = n;
  } else {
    m_i = n; n_i = m;
  }

  inca = incb = 1;
  INC_ADJUST(inca, $2_type)
  INC_ADJUST(incb, $3_type)
  MALLOC_VECTOR(a_vec, $2_type, m_i)
  for (i = 0; i < m_i*inca; i += inca) {
        SET_ZERO_VECTOR_ELEMENT(a_vec, i, $2_type)
  }
  MALLOC_VECTOR(b_vec, $3_type, m_i)
  for (i = 0; i < m_i*incb; i += incb) {
        SET_ZERO_VECTOR_ELEMENT(b_vec, i, $3_type)
  }

  if ((order == blas_colmajor && side == blas_left_side) ||
      (order == blas_rowmajor && side == blas_right_side)) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  INC_ADJUST(incci, $1_type)
  INC_ADJUST(inccij, $1_type)


  if (randomize == 0) {
    /* First fill in the first row of A and the first column/row of B */ 

    DOT_TESTGEN_NAME($1, $3, $2)(m_i, 0, 0, norm, blas_no_conj, 
                                 alpha, alpha_flag, beta, beta_flag, 
                                 b_vec, a_vec, seed, PASS_BY_REF(c_elem, $1_type), 
                                 PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));

    cij = 0;
    SET_VECTOR_ELEMENT(c_i, cij, c_elem, $1_type)
    SET_VECTOR_ELEMENT(r_true, cij, r_true_elem, EXTRA_TYPE($1_type))
 
    /* Copy a_vec to first row of A */
    $2sy_commit_row(order, uplo, m_i, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n_i; j++) {
      if (side == blas_left_side)
        $3ge_commit_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
      else
        $3ge_commit_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
    }

    /* Fill in rest of matrix A */
    cij = incci;
    for (i = 1; i < m_i; i++, cij += incci) {
      $2sy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      DOT_TESTGEN_NAME($1, $3, $2)(m_i, i, m_i-i, norm, 
                                   blas_no_conj, alpha, 1, 
                                   beta, 1, b_vec, a_vec, seed, 
                                   PASS_BY_REF(c_elem, $1_type), 
                                   PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));

      $2sy_commit_row(order, uplo, m_i, a, lda, a_vec, i);

      SET_VECTOR_ELEMENT(c_i, cij, c_elem, $1_type)
      SET_VECTOR_ELEMENT(r_true, cij, r_true_elem, EXTRA_TYPE($1_type))
    }

    /* Now fill in c and r_true */  
    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n_i; j++, cij += inccij) {
        GET_VECTOR_ELEMENT(c_elem, c_i, ci, $1_type)
        SET_VECTOR_ELEMENT(c_i, cij, c_elem, $1_type)
        COPY_VECTOR_ELEMENT(r_true, cij, r_true, ci, EXTRA_TYPE($1_type))
      } 
    }
  } else {

    ifelse(IS_MIXED($1, $2), `t', 
      `DECLARE_VECTOR(aa_vec, COMPLEX_TYPE($2_type))')
    ifelse(IS_MIXED($1, $3), `t', 
      `DECLARE_VECTOR(bb_vec, COMPLEX_TYPE($3_type))')

    ifelse(IS_MIXED($1, $2), `t', 
      `MALLOC_VECTOR(aa_vec, COMPLEX_TYPE($2_type), m_i)')
    ifelse(IS_MIXED($1, $3), `t', 
      `MALLOC_VECTOR(bb_vec, COMPLEX_TYPE($3_type), m_i)')

    if (alpha_flag == 0) {
      RANDOM(c_elem, $1_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
      SET_VECTOR_ELEMENT(alpha_i, 0, c_elem, $1_type)
    }
    if (beta_flag == 0) {
      RANDOM(c_elem, $1_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
      SET_VECTOR_ELEMENT(beta_i, 0, c_elem, $1_type)      
    }

    if ((order == blas_colmajor && side == blas_left_side) ||
        (order == blas_rowmajor && side == blas_right_side)) {
      incai = incbi = 1;
      incbij = ldb;
      incaij = lda;
    } else {
      incai = lda;
      incbi = ldb;
      incaij = incbij = 1;
    }

    INC_ADJUST(incbi, $3_type)
    INC_ADJUST(incbij, $3_type)
    INC_ADJUST(incai, $2_type)
    INC_ADJUST(incaij, $2_type)

    for (i = 0, ai = 0; i < m_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < m_i; j++, aij += incaij) {
        RANDOM(a_elem, $2_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
        SET_VECTOR_ELEMENT(a_i, aij, a_elem, $2_type) 
      }
    }

    for (i = 0, bi = 0; i < m_i; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n_i; j++, bij += incbij) {
        RANDOM(b_elem, $3_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
        SET_VECTOR_ELEMENT(b_i, bij, b_elem, $3_type) 
      }
    }

    for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
      $2sy_copy_row(order, uplo, m_i, a, lda, a_vec, i);
      ifelse(IS_MIXED($1, $2), `t', 
        `{
         int r;
         for (r = 0; r < m_i; r++) {
           aa_vec[2*r] = a_vec[r];
           aa_vec[2*r+1] = 0.0;
         }
        }')

      for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {

        if (side == blas_left_side)
          $3ge_copy_col(order, blas_no_trans, m, n, b, ldb, b_vec, j);
        else
          $3ge_copy_row(order, blas_no_trans, m, n, b, ldb, b_vec, j);
        ifelse(IS_MIXED($1, $3), `t', 
          `{
           int r;
           for (r = 0; r < m_i; r++) {
             bb_vec[2*r] = b_vec[r];
             bb_vec[2*r+1] = 0.0;
           }
          }')


        ifelse(IS_MIXED($1, $2, $3), `t', 
        `DOT_TESTGEN_NAME($1, $1, $1)(m_i, m_i, 0, norm, blas_no_conj, alpha, 1, 
                         beta, 1, 
                         ifelse(IS_MIXED($1, $3), `t', `bb_vec', `b_vec'), 
                         ifelse(IS_MIXED($1, $2), `t', `aa_vec', `a_vec'), 
                         seed, 
                         PASS_BY_REF(c_elem, $1_type), 
                         PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));', 
        `DOT_TESTGEN_NAME($1, $3, $2)(m_i, m_i, 0, norm, blas_no_conj, alpha, 1, 
                                     beta, 1, b_vec, a_vec, seed, 
                                     PASS_BY_REF(c_elem, $1_type), 
                                     PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));')

        SET_VECTOR_ELEMENT(c_i, cij, c_elem, $1_type)
        SET_VECTOR_ELEMENT(r_true, cij, r_true_elem, EXTRA_TYPE($1_type))
      }
    }

    ifelse(IS_MIXED($1, $2), `t', `FREE_VECTOR(aa_vec, COMPLEX_TYPE($2_type))')
    ifelse(IS_MIXED($1, $3), `t', `FREE_VECTOR(bb_vec, COMPLEX_TYPE($3_type))')
  }

  FREE_VECTOR(a_vec, $2_type)
  FREE_VECTOR(b_vec, $3_type)
}')dnl
dnl
dnl
dnl
define(`SYMM_TESTGEN', 
  `SYMM_TESTGEN_HEAD($1, $2, $3)
   SYMM_TESTGEN_COMMENT($1, $2, $3)
   SYMM_TESTGEN_BODY($1, $2, $3)')dnl
dnl
dnl
define(`PROTOTYPES', `FOREACH(`PREC_ARGS', `SYMM_TESTGEN_HEAD(arg);
')')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

FOREACH(`PREC_ARGS', `SYMM_TESTGEN(arg)
')
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
