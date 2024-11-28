dnl
dnl spmv_testgen.m4
dnl
dnl Test case generator for spmv routines.
dnl Generates test cases for alpha, A, beta, x and y and
dnl computes r_true in double-double precision.
dnl
dnl
dnl C Interface has the form 
dnl
dnl   BLAS_xspmv{_a_b}{_x}_testgen(int norm, 
dnl                             enum blas_order_type order, 
dnl                             enum blas_uplo_type uplo, 
dnl                             int n, int randomize, 
dnl                             SCALAR *alpha, int alpha_flag, 
dnl                             SCALAR *beta, int beta_flag, 
dnl                             ARRAY a, 
dnl                             ARRAY x, int incx, 
dnl                             ARRAY y, int incy, 
dnl                             int *seed, 
dnl                             double *r_true_l, double *r_true_t)
dnl
dnl Arguments
dnl   norm      (in) int
dnl
dnl   order     (in) blas_order_type
dnl                    determines the storage format for the A matrix
dnl   
dnl   uplo      (in) blas_uplo_type
dnl                    determines whether the upper triangular portion
dnl                    or the lower triangular portion of the symmetric
dnl                    matrix A is used.
dnl
dnl   n         (in) int
dnl                    the size of the
dnl                    vectors x and y is n
dnl                    matrix A is n-by-n.
dnl   
dnl   randomize (in) int
dnl                  if 0, entries in matrices A, x will be chosen for
dnl                        maximum cancellation, but with less randomness.
dnl                  if 1, every entry in the matrix A, x will be 
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
dnl
dnl   x, y      (out) vectors, x, y.
dnl   incx, incy (in) strides for vectors x, y.
dnl
dnl   seed      (in/out) int
dnl
dnl   r_true_l  (out) double *  (these are vectors of size n)
dnl   r_true_t  (out) double *
dnl               the leading/trailing part of the true in double-double
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl
dnl
define(`SYMV_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1symv_testgen', `BLAS_$1symv_$2_$3_testgen')')dnl
dnl
define(`SPMV_PACK_MATRIX_NAME', `$1spmv_pack_matrix')dnl
dnl
dnl
dnl
define(`SPMV_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1spmv$4', 
  `BLAS_$1spmv_$2_$3$4')')dnl
dnl
dnl
dnl  SPMV_TESTGEN
dnl     |
dnl     |-- SPMV_TESTGEN_HEAD
dnl     |      |
dnl     |      |-- SPMV_TESTGEN_NAME
dnl     |      |
dnl     |      |-- SPMV_TESTGEN_PARAMS
dnl     |      
dnl     |-- SPMV_TESTGEN_COMMENT
dnl     |
dnl     |-- SPMV_TESTGEN_BODY
dnl
dnl  Usage:
dnl    SPMV_TESTGEN        ($1, $2, $3)
dnl    SPMV_TESTGEN_HEAD   ($1, $2, $3)
dnl    SPMV_TESTGEN_NAME   ($1, $2, $3)
dnl    SPMV_TESTGEN_PARAMS ($1, $2, $3)
dnl    SPMV_TESTGEN_COMMENT($1, $2, $3)
dnl    SPMV_TESTGEN_BODY   ($1, $2, $3)
dnl
dnl    $1 -- type of alpha, beta, y
dnl    $2 -- type of a
dnl    $3 -- type of x
dnl
define(`SPMV_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1spmv_testgen', `BLAS_$1spmv_$2_$3_testgen')')dnl
dnl
dnl
define(`SPMV_TESTGEN_PARAMS',
  `int norm, enum blas_order_type order, dnl
   enum blas_uplo_type uplo, dnl
   int n,  int randomize, dnl
   $1_array alpha, int alpha_flag, $1_array beta, int beta_flag, dnl
   $2_array a, $3_array x, int incx, $1_array y, int incy, dnl
   int *seed, double *r_true_l, double *r_true_t')dnl
dnl
dnl
define(`SPMV_TESTGEN_HEAD', 
  `void SPMV_TESTGEN_NAME($1, $2, $3)(SPMV_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`SPMV_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to SPMV_NAME($1, $2, $3, `'){_x}
 *
 * Arguments
 * =========
 *
 * norm    (input) int
 *           = -1: the vectors are scaled with norms near underflow.
 *           =  0: the vectors have norms of order 1.
 *           =  1: the vectors are scaled with norms near overflow.
 *
 * order   (input) enum blas_order_type
 *           storage format of the matrices
 * 
 * uplo    (input) enum blas_uplo_type
 *           which half of the symmetric matrix a is to be stored.
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
 * a       (input/output) $2_array The Packed Symmetric Matrix
 * 
 *
 * x       (input/output) $3_array
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) $1_array
 *         generated vector y that will be used as an input to SPMV.
 * 
 * incy    (input) int
 *         leading dimension of vector y.
 *
 * seed    (input/output) int *
 *         seed for the random number generator.
 *
 * r_true_l  (output) double *
 *         the leading part of the truth in double-double.
 *
 * r_true_t  (output) double *
 *         the trailing part of the truth in double-double
 *
 */')dnl
dnl
dnl
dnl
dnl
define(`SPMV_TESTGEN_BODY', `{

  /* Strategy:  

                Use the SYMV generator, then simply pack the 
                        output.
  */

  DECLARE_VECTOR(a_full, $2_type)
  MALLOC_VECTOR(a_full, $2_type, n * n)

  SYMV_TESTGEN_NAME($1, $2, $3)(norm, order, uplo,
        n, randomize, alpha, alpha_flag, beta, beta_flag,
        a_full, n /* lda */, x, incx, y, incy,
        seed, r_true_l, r_true_t);
  SPMV_PACK_MATRIX_NAME($2)(order, uplo, n, a, a_full, n /*lda*/);

  FREE_VECTOR(a_full, $2_type)
}')dnl
dnl
dnl
dnl
dnl
define(`SPMV_TESTGEN', 
  `SPMV_TESTGEN_HEAD($1, $2, $3)
   SPMV_TESTGEN_COMMENT($1, $2, $3)
{
   SPMV_TESTGEN_BODY($1, $2, $3)
} /* end SPMV_TESTGEN_NAME($1, $2, $3) */')dnl
dnl
dnl
define(`PROTOTYPES', `FOREACH(`PREC_ARGS', `SPMV_TESTGEN_HEAD(arg);
')')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

FOREACH(`PREC_ARGS', `SPMV_TESTGEN(arg)
')')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
