dnl
dnl symv_testgen.m4
dnl
dnl Test case generator for symv routines.
dnl Generates test cases for alpha, A, beta, x and y and
dnl computes r_true in double-double precision.
dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
define(`SYMV_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1symv$4', 
  `BLAS_$1symv_$2_$3$4')')dnl
dnl
dnl
dnl  SYMV_TESTGEN
dnl     |
dnl     |-- SYMV_TESTGEN_HEAD
dnl     |      |
dnl     |      |-- SYMV_TESTGEN_NAME
dnl     |      |
dnl     |      |-- SYMV_TESTGEN_PARAMS
dnl     |      
dnl     |-- SYMV_TESTGEN_COMMENT
dnl     |
dnl     |-- SYMV_TESTGEN_BODY
dnl
dnl  Usage:
dnl    SYMV_TESTGEN        ($1, $2, $3)
dnl    SYMV_TESTGEN_HEAD   ($1, $2, $3)
dnl    SYMV_TESTGEN_NAME   ($1, $2, $3)
dnl    SYMV_TESTGEN_PARAMS ($1, $2, $3)
dnl    SYMV_TESTGEN_COMMENT($1, $2, $3)
dnl    SYMV_TESTGEN_BODY   ($1, $2, $3)
dnl
dnl    $1 -- type of alpha, beta, y
dnl    $2 -- type of a
dnl    $3 -- type of x
dnl
define(`SYMV_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1symv_testgen', `BLAS_$1symv_$2_$3_testgen')')dnl
dnl
dnl
define(`SYMV_TESTGEN_PARAMS',
  `int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n,  int randomize, 
   $1_array alpha, int alpha_flag, $1_array beta, int beta_flag, 
   $2_array a, int lda, $3_array x, int incx, $1_array y, int incy, 
   int *seed, double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
dnl
define(`SYMV_TESTGEN_HEAD', 
  `void SYMV_TESTGEN_NAME($1, $2, $3)(SYMV_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`SYMV_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to SYMV_NAME($1, $2, $3, `'){_x}
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
 *           which half of the symvetric matrix a is to be stored.
 *
 * n       (input) int
 *           sizes of symmetrical matrix a, size of vectors x, y:
 *              matrix a is n-by-n.
 * 
 * randomize (input) int
 *           if 0, entries in matrices A, x will be chosen for
 *              maximum cancellation, but with less randomness.
 *           if 1, every entry in the matrix A, x will be 
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
 * x       (input/output) $3_array
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) $1_array
 *         generated vector y that will be used as an input to SYMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
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
define(`SYMV_TESTGEN_BODY', `{

  int i, j;
  int yi;
  int xi;
  int aij, ai, ri;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int incaij, incai;
  int inca;
  int n_i;

  DECLARE(y_elem, $1_type)
  DECLARE(a_elem, $2_type)
  DECLARE(x_elem, $3_type)
  DECLARE(r_true_elem, EXTRA_TYPE($1_type))

  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(x_vec, $3_type)

  PTR_CAST(y, $1_type)
  PTR_CAST(alpha, $1_type)
  PTR_CAST(beta, $1_type)
  PTR_CAST(a, $2_type)
  PTR_CAST(x, $3_type)

  n_i = n;

  /*x_vec, a_vec must have stride of 1*/
  inca = 1;
  INC_ADJUST(inca, $2_type)

  MALLOC_VECTOR(a_vec, $2_type, n_i)
  for (i = 0; i < n_i; i += inca) {
        SET_ZERO_VECTOR_ELEMENT(a_vec, i, $2_type)
  }

  incyi = incy;
  INC_ADJUST(incyi, $1_type)
  if (incyi < 0) {
        y_starti = (-n+1) * incyi;
  }
  else {
        y_starti = 0;
  }



  incri = 1;
  INC_ADJUST(incri, EXTRA_TYPE($1_type))

  incxi = incx;
  incx_veci = 1;
  INC_ADJUST(incx_veci, $3_type)
  INC_ADJUST(incxi, $3_type)

  if (incxi < 0) {
        x_starti = (-n+1) * incxi;
  }
  else {
        x_starti = 0;
  }

  MALLOC_VECTOR(x_vec, $3_type, n_i);
  for (i = 0; i < n_i; i += incx_veci) {
        SET_ZERO_VECTOR_ELEMENT(x_vec, i, $3_type)
  }

  if (randomize == 0) {
    /* First fill in the first row of A and all of x */ 

    DOT_TESTGEN_NAME($1, $3, $2)(n_i, 0, 0, norm, blas_no_conj, 
                                 alpha, alpha_flag, beta, beta_flag, 
                                 x_vec, a_vec, seed, PASS_BY_REF(y_elem, $1_type), 
                                 PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));

    yi = y_starti;
    ri = 0;
    SET_VECTOR_ELEMENT(y_i, yi, y_elem, $1_type)
    SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))
 
    /* Copy a_vec to first row of A */
    $2sy_commit_row(order, uplo, n_i, a, lda, a_vec, 0);

    /* commit x_vec to x */
    $3copy_vector(x_vec, n_i, 1, x_i, incx);

    /* Fill in rest of matrix A */
    for(i = 1, yi += incyi, ri = incri; i < n_i; i++, ri += incri, 
                yi += incyi) {
      $2sy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      DOT_TESTGEN_NAME($1, $3, $2)(n_i, i, n_i-i, norm, 
                                   blas_no_conj, alpha, 1, 
                                   beta, 1, x_vec, a_vec, seed, 
                                   PASS_BY_REF(y_elem, $1_type), 
                                   PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));

      $2sy_commit_row(order, uplo, n_i, a, lda, a_vec, i);

        /*commits an element to the generated y */
      SET_VECTOR_ELEMENT(y_i, yi, y_elem, $1_type)
      SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))
    }

  } else {

    dnl NOTE: There are several issues with the dot test generator that
    dnl       affects the code in this section:
    dnl
    dnl         - in mixed precisions (i.e., single mixed with double)
    dnl           cases, alpha, beta, and the vector with higher precision
    dnl           are rounded to single precision.
    dnl
    dnl         - in mixed complex / real cases, alpha, beta, and complex
    dnl           vector are not truly random: they are multiples of (1+i)
    dnl           or (1-i), or something similar.
    dnl
    dnl       This is why this section is not as clean since we need to do
    dnl       certain things to comply with the assumption of the dot test
    dnl       generator (which should be eventually fixed).
    ifelse(IS_MIXED($1, $2), `t', 
      `DECLARE_VECTOR(aa_vec, COMPLEX_TYPE($2_type))')
    ifelse(IS_MIXED($1, $3), `t', 
      `DECLARE_VECTOR(xx_vec, COMPLEX_TYPE($3_type))')

    ifelse(IS_MIXED($1, $2), `t', 
      `MALLOC_VECTOR(aa_vec, COMPLEX_TYPE($2_type), n_i)')
    ifelse(IS_MIXED($1, $3), `t', 
      `MALLOC_VECTOR(xx_vec, COMPLEX_TYPE($3_type), n_i)')

        /* randomly select alpha, beta */
    if (alpha_flag == 0) {
      RANDOM(y_elem, $1_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
      SET_VECTOR_ELEMENT(alpha_i, 0, y_elem, $1_type)
    }
    if (beta_flag == 0) {
      RANDOM(y_elem, $1_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
      SET_VECTOR_ELEMENT(beta_i, 0, y_elem, $1_type)      
    }


    /*set a, x randomly */

    if(order == blas_colmajor)
    {
        incai = 1;
        incaij = lda;
    } else {
        incai = lda;
        incaij = 1;
    }   

    incxi = incx;
    INC_ADJUST(incxi, $3_type)

    if (incxi < 0) {
          x_starti = (-n+1) * incxi;
    }
    else {
          x_starti = 0;
    }

    INC_ADJUST(incai, $2_type)
    INC_ADJUST(incaij, $2_type)

    for (i = 0, ai = 0; i < n_i; i++, ai += incai) {
      for (j = 0, aij = ai; j < n_i; j++, aij += incaij) {
        RANDOM(a_elem, $2_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
        SET_VECTOR_ELEMENT(a_i, aij, a_elem, $2_type) 
      }
    }

    for (i = 0, xi = x_starti; i < n_i; i++, xi += incxi) {
      RANDOM(x_elem, $3_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
      SET_VECTOR_ELEMENT(x_i, xi, x_elem, $3_type)
    }

    /* now compute appropriate y vector */
      
      /* get x */
      $3copy_vector(x_i, n_i, incx, x_vec, 1);
      ifelse(IS_MIXED($1, $3), `t', 
        `{
          /* promote to complex */
         int r;
         for (r = 0; r < n_i; r++) {
           xx_vec[2*r] = x_vec[r];
           xx_vec[2*r+1] = 0.0;
         }
        }')

    for (i = 0, yi = y_starti, ri = 0; 
        i < n_i; i++, yi += incyi, ri += incri) {
      $2sy_copy_row(order, uplo, n_i, a, lda, a_vec, i);
      ifelse(IS_MIXED($1, $2), `t', 
        `{
                /* promote to complex */
         int r;
         for (r = 0; r < n_i; r++) {
           aa_vec[2*r] = a_vec[r];
           aa_vec[2*r+1] = 0.0;
         }
        }')

      ifelse(IS_MIXED($1, $2, $3), `t', 
      `DOT_TESTGEN_NAME($1, $1, $1)(n_i, n_i, 0, norm, blas_no_conj, alpha, 1, 
                       beta, 1, 
                       ifelse(IS_MIXED($1, $3), `t', `xx_vec', `x_vec'), 
                       ifelse(IS_MIXED($1, $2), `t', `aa_vec', `a_vec'), 
                       seed, 
                       PASS_BY_REF(y_elem, $1_type), 
                       PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));', 
      `DOT_TESTGEN_NAME($1, $3, $2)(n_i, n_i, 0, norm, blas_no_conj, alpha, 1, 
                       beta, 1, x_vec, a_vec, seed, 
                       PASS_BY_REF(y_elem, $1_type), 
                       PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));')

      SET_VECTOR_ELEMENT(y_i, yi, y_elem, $1_type)
      SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))    
    }

    ifelse(IS_MIXED($1, $2), `t', `FREE_VECTOR(aa_vec, COMPLEX_TYPE($2_type))')
    ifelse(IS_MIXED($1, $3), `t', `FREE_VECTOR(xx_vec, COMPLEX_TYPE($3_type))')
  }

  FREE_VECTOR(a_vec, $2_type)
  FREE_VECTOR(x_vec, $3_type)
}')dnl
dnl
dnl
define(`SYMV_TESTGEN', 
  `SYMV_TESTGEN_HEAD($1, $2, $3)
   SYMV_TESTGEN_COMMENT($1, $2, $3)
   SYMV_TESTGEN_BODY($1, $2, $3)')dnl
dnl
dnl
define(`PROTOTYPES', `FOREACH(`PREC_ARGS', `SYMV_TESTGEN_HEAD(arg);
')')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

FOREACH(`PREC_ARGS', `SYMV_TESTGEN(arg)
')
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
