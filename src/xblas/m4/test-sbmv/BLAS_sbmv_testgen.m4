dnl
dnl sbmv_testgen.m4
dnl
dnl Test case generator for sbmv routines.
dnl Generates test cases for alpha, A, beta, x and y and
dnl computes r_true in double-double precision.
dnl
dnl
dnl C Interface has the form 
dnl
dnl   BLAS_xsbmv{_a_b}{_x}_testgen(int norm, 
dnl                             enum blas_order_type order, 
dnl                             enum blas_uplo_type uplo, 
dnl                             int n, int k,int randomize, 
dnl                             SCALAR *alpha, int alpha_flag, 
dnl                             SCALAR *beta, int beta_flag, 
dnl                             ARRAY a, int lda, 
dnl                             ARRAY x, int incx, 
dnl                             ARRAY y, int incy, 
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
dnl                    or the lower triangular portion of the sbmvetric
dnl                    matrix A is used.
dnl
dnl   n         (in) int
dnl                    the size of the
dnl                    vectors x and y is n
dnl                    matrix A is n-by-n.
dnl
dnl   k         (in) int
dnl                     the number of superdiagonals in
dnl                     the symmetric matrix A.
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
dnl   a         (out) matrix a.
dnl   lda       (in) leading dimensions of matrix a.
dnl
dnl   x, y      (out) vectors, x, y.
dnl   incx, incy (in) strides for vectors x, y.
dnl
dnl   seed      (in/out) int
dnl
dnl   HEAD(r_true)  (out) double *  (these are vectors of size n)
dnl   TAIL(r_true)  (out) double *
dnl               the leading/trailing part of the true in double-double
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
define(`SBMV_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1sbmv$4', 
  `BLAS_$1sbmv_$2_$3$4')')
dnl
dnl
dnl  SBMV_TESTGEN
dnl     |
dnl     |-- SBMV_TESTGEN_HEAD
dnl     |      |
dnl     |      |-- SBMV_TESTGEN_NAME
dnl     |      |
dnl     |      |-- SBMV_TESTGEN_PARAMS
dnl     |      
dnl     |-- SBMV_TESTGEN_COMMENT
dnl     |
dnl     |-- SBMV_TESTGEN_BODY
dnl
dnl  Usage:
dnl    SBMV_TESTGEN        ($1, $2, $3)
dnl    SBMV_TESTGEN_HEAD   ($1, $2, $3)
dnl    SBMV_TESTGEN_NAME   ($1, $2, $3)
dnl    SBMV_TESTGEN_PARAMS ($1, $2, $3)
dnl    SBMV_TESTGEN_COMMENT($1, $2, $3)
dnl    SBMV_TESTGEN_BODY   ($1, $2, $3)
dnl
dnl    $1 -- type of alpha, beta, y
dnl    $2 -- type of a
dnl    $3 -- type of x
dnl
define(`SBMV_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1sbmv_testgen', `BLAS_$1sbmv_$2_$3_testgen')')dnl
dnl
dnl
define(`SBMV_TESTGEN_PARAMS',
  `int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n, int randomize, 
   $1_array alpha, int alpha_flag, $1_array beta, int beta_flag, 
   $2_array a, int k, int lda, $3_array x, int incx, $1_array y, int incy, 
   int *seed, double *HEAD(r_true), double *TAIL(r_true)')dnl
dnl
dnl
define(`SBMV_TESTGEN_HEAD', 
  `void SBMV_TESTGEN_NAME($1, $2, $3)(SBMV_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`SBMV_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to SBMV_NAME($1, $2, $3, `'){_x}
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
 *           which half of the sbmvetric matrix a is to be stored.
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
 *         generated vector y that will be used as an input to SBMV.
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
define(`SBMV_TESTGEN_BODY', `{

  int i, j;
  int yi;
  int xi;
  int ri;
  int alpha_fixed, beta_fixed;
  int incyi, incri;
  int incxi, incx_veci, x_starti, y_starti;
  int inca, inca_vec, a_veci;
  int n_i;
  int n_elements, x_fixed;

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
  alpha_fixed = alpha_flag;
  beta_fixed = beta_flag;

  /*x_vec, a_vec must have stride of 1*/
  inca = 1;
  INC_ADJUST(inca, $2_type)

  MALLOC_VECTOR(a_vec, $2_type, 2*n_i)
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

  MALLOC_VECTOR(x_vec, $3_type, 2*n_i);
  for (xi = 0, i = 0; i < n_i; i++, xi += incx_veci) {
        SET_ZERO_VECTOR_ELEMENT(x_vec, xi, $3_type)
  }

  if (randomize == 0) {

        /* fill in Matrix A */
    for(i = 0, yi = y_starti, ri = 0; i < n_i; i++, ri += incri, 
                yi += incyi) {
      $2sbmv_copy_row(order, uplo, n_i, a_i, k, lda, a_vec, i);

        /* n_elements is the total number of x values we will
        *       have after this iteration of the loop.
        *  x_fixed is the number of fixed x values at
        *     the beginning of the loop
        */
        n_elements = MIN(i + k + 1, n_i);
        x_fixed = (i != 0) ? MIN(i + k, n_i) : 0;

      DOT_TESTGEN_NAME($1, $3, $2)(n_elements, i, 
                                   x_fixed-i, norm, 
                                   blas_no_conj, alpha, alpha_fixed, 
                                   beta, beta_fixed, x_vec, a_vec, seed, 
                                   PASS_BY_REF(y_elem, $1_type), 
                                   PASS_BY_REF(r_true_elem, EXTRA_TYPE($1_type)));
                                
        beta_fixed = alpha_fixed = 1;

        /* ignores portion that should be zero */
      $2sbmv_commit_row(order, uplo, n_i, a_i, k, lda, a_vec, i);

        /*commits an element to the generated y */
      SET_VECTOR_ELEMENT(y_i, yi, y_elem, $1_type)
      SET_VECTOR_ELEMENT(r_true, ri, r_true_elem, EXTRA_TYPE($1_type))
    }
        /* copy x_vec to output vector x */
    $3copy_vector(x_vec, n_i, 1, x_i, incx);

  } else {

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

    incxi = incx;
    INC_ADJUST(incxi, $3_type)

    if (incxi < 0) {
          x_starti = (-n+1) * incxi;
    }
    else {
          x_starti = 0;
    }
    inca_vec = 1;
    INC_ADJUST(inca_vec, $2_type)

    for (i = 0; i < n_i; i++) {
      for (j = 0, a_veci = 0; j < n_i; j++, a_veci += inca_vec) {
        RANDOM(a_elem, $2_type, IS_MIXED_PREC($1_type, $2_type, $3_type))
        SET_VECTOR_ELEMENT(a_vec, a_veci, a_elem, $2_type) 
      }
      $2sbmv_commit_row(order, uplo, n_i, a_i, k, lda, a_vec, i);
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
      $2sbmv_copy_row(order, uplo, n_i, a_i, k, lda, a_vec, i);
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
define(`SBMV_TESTGEN', 
  `SBMV_TESTGEN_HEAD($1, $2, $3)
   SBMV_TESTGEN_COMMENT($1, $2, $3)
   SBMV_TESTGEN_BODY($1, $2, $3)')dnl
dnl
dnl
define(`PROTOTYPES', `FOREACH(`PREC_ARGS', `SBMV_TESTGEN_HEAD(arg);
')')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

FOREACH(`PREC_ARGS', `SBMV_TESTGEN(arg)
')
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
