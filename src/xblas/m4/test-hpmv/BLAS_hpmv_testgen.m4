dnl
dnl hpmv_testgen.m4
dnl
dnl Test case generator for hpmv routines.
dnl Generates test cases for alpha, A, beta, x and y and
dnl computes r_true in double-double precision.
dnl
dnl
dnl C Interface has the form 
dnl
dnl   BLAS_xhpmv{_a_b}{_x}_testgen(int norm, 
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
dnl                    or the lower triangular portion of the hermitian
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
define(`HEMV_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1hemv_testgen', `BLAS_$1hemv_$2_$3_testgen')')dnl
dnl
define(`HPMV_PACK_MATRIX_NAME', `$1hpmv_pack_matrix')dnl
dnl
dnl
dnl
define(`HPMV_NAME', `ifelse(
  `$2&&$3', `$1&&$1', `BLAS_$1hpmv$4', 
  `BLAS_$1hpmv_$2_$3$4')')dnl
dnl
dnl
dnl  HPMV_TESTGEN
dnl     |
dnl     |-- HPMV_TESTGEN_HEAD
dnl     |      |
dnl     |      |-- HPMV_TESTGEN_NAME
dnl     |      |
dnl     |      |-- HPMV_TESTGEN_PARAMS
dnl     |      
dnl     |-- HPMV_TESTGEN_COMMENT
dnl     |
dnl     |-- HPMV_TESTGEN_BODY
dnl
dnl  Usage:
dnl    HPMV_TESTGEN        ($1, $2, $3)
dnl    HPMV_TESTGEN_HEAD   ($1, $2, $3)
dnl    HPMV_TESTGEN_NAME   ($1, $2, $3)
dnl    HPMV_TESTGEN_PARAMS ($1, $2, $3)
dnl    HPMV_TESTGEN_COMMENT($1, $2, $3)
dnl    HPMV_TESTGEN_BODY   ($1, $2, $3)
dnl
dnl    $1 -- type of alpha, beta, y
dnl    $2 -- type of a
dnl    $3 -- type of x
dnl
define(`HPMV_TESTGEN_NAME', `ifelse(`$2&&$3', `$1&&$1', 
  `BLAS_$1hpmv_testgen', `BLAS_$1hpmv_$2_$3_testgen')')dnl
dnl
dnl
define(`HPMV_TESTGEN_PARAMS',
  `int norm, enum blas_order_type order, 
   enum blas_uplo_type uplo,
   int n,  int randomize, 
   $1_array alpha, int alpha_flag, $1_array beta, int beta_flag, 
   $2_array a, $3_array x, int incx, $1_array y, int incy, 
   int *seed, double *r_true_l, double *r_true_t')dnl
dnl
dnl
define(`HPMV_TESTGEN_HEAD', 
  `void HPMV_TESTGEN_NAME($1, $2, $3)(HPMV_TESTGEN_PARAMS($1, $2, $3))')dnl
dnl
dnl
define(`HPMV_TESTGEN_COMMENT', `
/*
 * Purpose
 * =======
 *
 *   Generates the test inputs to HPMV_NAME($1, $2, $3, `'){_x}
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
 *           which half of the hermitian matrix a is to be stored.
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
 * a       (input/output) $2_array The Packed Hermitian Matrix
 * 
 *
 * x       (input/output) $3_array
 *
 * incx    (input) int
 *         stride of vector x.
 * 
 * y       (input/output) $1_array
 *         generated vector y that will be used as an input to HPMV.
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
define(`HPMV_TESTGEN_BODY', `{

  /* Strategy:  
         R1 = alpha * A1 * x + beta * y1
         R2 = alpha * A2 * x + beta * y2
     where all the matrices and vectors are real.  Then let R = R1 + i R2, 
     A = A1 + i A2, y = y1 + i y2.  To make A hermitian, A1 is 
     symmetric, and A2 is a skew matrix (trans(A2) = -A2).

        hemv_testgen_body does all that, so just call it, then
        finally, at the end, pack it all properly.
  */

  DECLARE_VECTOR(a_full, $2_type)
  MALLOC_VECTOR(a_full, $2_type, n * n)

  HEMV_TESTGEN_NAME($1, $2, $3)(norm, order, uplo,
        n, randomize, alpha, alpha_flag, beta, beta_flag,
        a_full, n /* lda */, x, incx, y, incy,
        seed, r_true_l, r_true_t);
  HPMV_PACK_MATRIX_NAME($2)(order, uplo, n, a, a_full, n /*lda*/);

  FREE_VECTOR(a_full, $2_type)
}
')dnl
dnl
dnl
dnl
dnl
define(`HPMV_TESTGEN', 
  `HPMV_TESTGEN_HEAD($1, $2, $3)
   HPMV_TESTGEN_COMMENT($1, $2, $3)
{
dnl   char *routine_name = "HPMV_TESTGEN_NAME($1, $2, $3)";
   HPMV_TESTGEN_BODY($1, $2, $3)
} /* end HPMV_TESTGEN_NAME($1, $2, $3) */')dnl
dnl
dnl
define(`PROTOTYPES', `
HPMV_TESTGEN_HEAD(c, c, c);
HPMV_TESTGEN_HEAD(z, z, z);
HPMV_TESTGEN_HEAD(z, c, z);
HPMV_TESTGEN_HEAD(z, z, c);
HPMV_TESTGEN_HEAD(z, c, c);
HPMV_TESTGEN_HEAD(z, z, d);
HPMV_TESTGEN_HEAD(c, c, s);
')dnl
dnl
dnl
define(`SOURCE', `
#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

HPMV_TESTGEN(c, c, c)
HPMV_TESTGEN(z, z, z)
HPMV_TESTGEN(z, c, z)
HPMV_TESTGEN(z, z, c)
HPMV_TESTGEN(z, c, c)
HPMV_TESTGEN(z, z, d)
HPMV_TESTGEN(c, c, s)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
