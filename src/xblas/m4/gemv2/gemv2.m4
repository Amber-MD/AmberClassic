dnl ----------------------------------------------------------------------
dnl GEMV2 -- General matrix 2-vector multiply
dnl   y <--- alpha * op(A) * head_x + alpha * op(A) * tail_x + beta * y
dnl
dnl   where op can be no-op, tranpose, or conjugate transpose
dnl ----------------------------------------------------------------------
dnl
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(gemv2-common.m4)dnl
dnl
dnl
dnl  Usage:
dnl    GEMV2        ($1, $2, $3, $4)
dnl    GEMV2_HEAD   ($1, $2, $3, $4)
dnl    GEMV2_NAME   ($1, $2, $3, $4)
dnl    GEMV2_PARAMS ($1, $2, $3, $4)
dnl    GEMV2_COMMENT($1, $2, $3, $4)
dnl
dnl      $1 -- type of alpha, beta, y.
dnl      $2 -- type of A
dnl      $3 -- type of head_x, tail_x
dnl      $4 -- Set to `_x' for _x routines.  Otherwise set to `'.
dnl
dnl
define(`GEMV2_COMMENT', `
/*
 * Purpose
 * =======
 *
 * Computes y = alpha * op(A) * head_x + alpha * op(A) * tail_x + beta * y,
 * where A is a general matrix.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans        (input) blas_trans_type
 *              Transpose of A: no trans, trans, or conjugate trans
 *
 * m            (input) int
 *              Dimension of A
 *
 * n            (input) int
 *              Dimension of A and the length of vector x and z
 *
 * alpha        (input) $1_scalar
 *              
 * A            (input) const $2_array
 *
 * lda          (input) int 
 *              Leading dimension of A
 *
 * head_x
 * tail_x       (input) const $3_array
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) $1_scalar
 *
 * y            (input) const $1_array
 * 
 * incy         (input) int
 *              The stride for vector y.
 * 
PREC_COMMENT($4)dnl
 */')dnl
dnl
dnl
dnl  Usage: GEMV2_BODY($1, $2, $3, $4, $5, $6)
dnl  Generates the main body of the product code.
dnl    $1 - type of alpha, beta, y
dnl    $2 - type of A
dnl    $3 - type of x
dnl    $4 - type of sum/prod
dnl    $5 - type of temp
dnl    $6 - [optional] String `FPU' is passed if FPU fix 
dnl         is needed.  Empty string is passed otherwise.
dnl
define(`GEMV2_BODY', `
  int i, j;
  int iy, jx, kx, ky;
  int lenx, leny;
  int ai, aij;
  int incai, incaij;

  PTR_CAST(a, $2, `const')
  PTR_CAST(head_x, $3, `const')
  PTR_CAST(tail_x, $3, `const')
  PTR_CAST(y, $1)
  SCALAR_CAST(alpha, $1)
  SCALAR_CAST(beta, $1)
  DECLARE(a_elem, $2)
  DECLARE(x_elem, $3)
  DECLARE(y_elem, $1)
  DECLARE(prod, $4)
  DECLARE(sum, $4)
  DECLARE(sum2, $4)
  DECLARE(tmp1, $5)
  DECLARE(tmp2, $5)
  ifelse(`$6', `FPU', `FPU_FIX_DECL;')

  /* all error calls */
  if (m < 0)
    BLAS_error(routine_name, -3, m, 0);
  else if (n <= 0)
    BLAS_error(routine_name, -4, n, 0);
  else if (incx == 0)
    BLAS_error(routine_name, -10, incx, 0);
  else if (incy == 0)
    BLAS_error(routine_name, -13, incy, 0);

  if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
    lenx = n;
    leny = m;
    incai = lda;
    incaij = 1;
  } else if ((order == blas_rowmajor) && (trans != blas_no_trans)) {
    lenx = m;
    leny = n;
    incai = 1;
    incaij = lda;
  } else if ((order == blas_colmajor) && (trans == blas_no_trans)) {
    lenx = n;
    leny = m;
    incai = 1;
    incaij = lda;
  } else {    /* colmajor and blas_trans */
    lenx = m;
    leny = n;
    incai = lda;
    incaij = 1;
  }

  if (lda < leny)
    BLAS_error(routine_name,  -7,  lda, NULL);

  ifelse(`$6', `FPU', `FPU_FIX_START;')

  INC_ADJUST(incx, $3)
  INC_ADJUST(incy, $1)
  INC_ADJUST(incai, $2)
  INC_ADJUST(incaij, $2)

  if (incx > 0) kx = 0;
  else kx = (1 - lenx)*incx;
  if (incy > 0) ky = 0;
  else ky = (1 - leny)*incy;

  /* No extra-precision needed for alpha = 0 */
  if ( TEST_0(alpha_i, $1) ) {
      if ( TEST_0(beta_i, $1) ) {
          iy = ky;
          for (i = 0; i < leny; i++) {
              SET_ZERO_VECTOR_ELEMENT(y_i, iy, $1)
              iy += incy;
          }
      } else if ( !(TEST_0(beta_i, $1)) ) {
          iy = ky;
          for (i = 0; i < leny; i++) {
              GET_VECTOR_ELEMENT(y_elem, y_i, iy, $1)
              MUL(tmp1, $5, y_elem, $1, beta_i, $1)
              SET_ROUND_VECTOR_ELEMENT(y_i, iy, tmp1, $5)
              iy += incy;
          }
      }
  } else { /* alpha != 0 */
    IF_COMPLEX($2, 
      `if (trans == blas_conj_trans){
        GEMV2_LOOP($1, $2, $3, $4, $5, blas_conj)
      } else{
        GEMV2_LOOP($1, $2, $3, $4, $5, blas_no_conj)
      }', 
      `GEMV2_LOOP($1, $2, $3, $4, $5, blas_no_conj)')
  }

  ifelse(`$6', `FPU', `FPU_FIX_STOP;')
')dnl
dnl
dnl
define(`GEMV2_LOOP',`
  /* if beta = 0, we can save m multiplies:
     y = alpha*A*head_x + alpha*A*tail_x  */
  if ( TEST_0(beta_i, $1) ) {
      if ( TEST_1(alpha_i, $1) ) {
          /* save m more multiplies if alpha = 1 */
          ai = 0;
          iy = ky;
          for (i = 0; i < leny; i++) {
              ZERO(sum, $4)
              ZERO(sum2, $4)
              aij = ai;
              jx = kx;
              for (j = 0; j < lenx; j++) {
                  GET_VECTOR_ELEMENT(a_elem, a_i, aij, $2)
                  CONJ(a_elem, $2, $6)
                  GET_VECTOR_ELEMENT(x_elem, head_x_i, jx, $3)
                  MUL(prod, $4, a_elem, $2, x_elem, $3)
                  ADD(sum, $4, sum, $4, prod, $4)
                  GET_VECTOR_ELEMENT(x_elem, tail_x_i, jx, $3)
                  MUL(prod, $4, a_elem, $2, x_elem, $3)
                  ADD(sum2, $4, sum2, $4, prod, $4)
                  aij += incaij;
                  jx += incx;
              }
              ADD(sum, $4, sum, $4, sum2, $4)
              ifelse(`$5', `$4', 
              `SET_ROUND_VECTOR_ELEMENT(y_i, iy, sum, $4)',
              `ASSIGN(tmp1, $5, sum, $4)
              SET_ROUND_VECTOR_ELEMENT(y_i, iy, tmp1, $5)')
              ai += incai;
              iy += incy;
          } /* end for */
      } else { /* alpha != 1 */
          ai = 0;
          iy = ky;
          for (i = 0; i < leny; i++) {
              ZERO(sum, $4)
              ZERO(sum2, $4)
              aij = ai;
              jx = kx;
              for (j = 0; j < lenx; j++) {
                  GET_VECTOR_ELEMENT(a_elem, a_i, aij, $2)
                  CONJ(a_elem, $2, $6)
                  GET_VECTOR_ELEMENT(x_elem, head_x_i, jx, $3)
                  MUL(prod, $4, a_elem, $2, x_elem, $3)
                  ADD(sum, $4, sum, $4, prod, $4)
                  GET_VECTOR_ELEMENT(x_elem, tail_x_i, jx, $3)
                  MUL(prod, $4, a_elem, $2, x_elem, $3)
                  ADD(sum2, $4, sum2, $4, prod, $4)
                  aij += incaij;
                  jx += incx;
              }
              MUL(tmp1, $5, sum, $4, alpha_i, $1)
              MUL(tmp2, $5, sum2, $4, alpha_i, $1)
              ADD(tmp1, $5, tmp1, $5, tmp2, $5)
              SET_ROUND_VECTOR_ELEMENT(y_i, iy, tmp1, $5)
              ai += incai;
              iy += incy;
          }
      }
  } else { /* beta != 0 */
      if ( TEST_1(alpha_i, $1) ) {
          /* save m multiplies if alpha = 1 */
          ai = 0;
          iy = ky;
          for (i = 0; i < leny; i++) {
              ZERO(sum, $4);
              ZERO(sum2, $4);
              aij = ai;
              jx = kx;
              for (j = 0; j < lenx; j++) {
                  GET_VECTOR_ELEMENT(a_elem, a_i, aij, $2)
                  CONJ(a_elem, $2, $6)
                  GET_VECTOR_ELEMENT(x_elem, head_x_i, jx, $3)
                  MUL(prod, $4, a_elem, $2, x_elem, $3)
                  ADD(sum, $4, sum, $4, prod, $4)
                  GET_VECTOR_ELEMENT(x_elem, tail_x_i, jx, $3)
                  MUL(prod, $4, a_elem, $2, x_elem, $3)
                  ADD(sum2, $4, sum2, $4, prod, $4)
                  aij += incaij;
                  jx += incx;
              }
              ADD(sum, $4, sum, $4, sum2, $4)
              GET_VECTOR_ELEMENT(y_elem, y_i, iy, $1)
              MUL(tmp1, $5, y_elem, $1, beta_i, $1)
              ifelse(`$5', `$4', 
              `ADD(tmp2, $5, sum, $4, tmp1, $5)',
              `ASSIGN(tmp2, $5, sum, $4)
               ADD(tmp2, $5, tmp2, $5, tmp1, $5)')
              SET_ROUND_VECTOR_ELEMENT(y_i, iy, tmp2, $5)
              ai += incai;
              iy += incy;
          }
      } else { /* alpha != 1, the most general form:
                  y = alpha*A*head_x + alpha*A*tail_x + beta*y */
          ai = 0;
          iy = ky;
          for (i = 0; i < leny; i++) {
              ZERO(sum, $4);
              ZERO(sum2, $4);
              aij = ai;
              jx = kx;
              for (j = 0; j < lenx; j++) {
                  GET_VECTOR_ELEMENT(a_elem, a_i, aij, $2)
                  CONJ(a_elem, $2, $6)
                  GET_VECTOR_ELEMENT(x_elem, head_x_i, jx, $3)
                  MUL(prod, $4, a_elem, $2, x_elem, $3)
                  ADD(sum, $4, sum, $4, prod, $4)
                  GET_VECTOR_ELEMENT(x_elem, tail_x_i, jx, $3)
                  MUL(prod, $4, a_elem, $2, x_elem, $3)
                  ADD(sum2, $4, sum2, $4, prod, $4)
                  aij += incaij;
                  jx += incx;
              }
              MUL(tmp1, $5, sum, $4, alpha_i, $1)
              MUL(tmp2, $5, sum2, $4, alpha_i, $1)
              ADD(tmp1, $5, tmp1, $5, tmp2, $5)
              GET_VECTOR_ELEMENT(y_elem, y_i, iy, $1)
              MUL(tmp2, $5, y_elem, $1, beta_i, $1)
              ADD(tmp1, $5, tmp1, $5, tmp2, $5)
              SET_ROUND_VECTOR_ELEMENT(y_i, iy, tmp1, $5)
              ai += incai;
              iy += incy;
          }
      }
  }
')dnl
dnl
dnl
dnl  Usage: SWITCH_prec($1, $2, $3, $4, $5, $6, $7, $8, $9)
dnl  Generates a 3-way switch statement based on prec.
dnl    $1      --  type of alpha, beta, r
dnl    $2      --  type of A
dnl    $3      --  type of x, y
dnl    $4, $5  --  type of `sum' and `tmp' in single case
dnl    $6, $7  --  type of `sum' and `tmp' in double/indigenous case
dnl    $8, $9  --  type of `sum' and `tmp' in extra case
dnl
define(`SWITCH_prec', 
       `switch ( prec ) {
        case blas_prec_single: ifelse(`$4&&$5', `$6&&$7', `', `{
          GEMV2_BODY($1, $2, $3, $4, $5)
          break;
        }
        ')dnl
        case blas_prec_double:
        case blas_prec_indigenous: {
          GEMV2_BODY($1, $2, $3, $6, $7)
          break;
        }
        case blas_prec_extra: { 
          GEMV2_BODY($1, $2, $3, $8, $9, FPU) }
          break;
        }')dnl
dnl
dnl
dnl
dnl  Usage: GEMV2_X_BODY($1, $2, $3)
dnl  Generates the main body of the extended version of dot code.
dnl    $1 -- type of alpha, beta, y
dnl    $2 -- type of A
dnl    $3 -- type of x
dnl
define(`GEMV2_X_BODY', 
  `SWITCH_prec($1, $2, $3, 
    SUM_TYPE_X($1, $2, $3, S), TMP_TYPE_X($1, S), 
    SUM_TYPE_X($1, $2, $3, D), TMP_TYPE_X($1, D), 
    SUM_TYPE_X($1, $2, $3, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
define(`GEMV2', 
  `GEMV2_HEAD($1, $2, $3, $4)
   GEMV2_COMMENT($1, $2, $3, $4)
   {
     static const char routine_name[] = "GEMV2_NAME($1, $2, $3, $4)";
     ifelse($4, _x, `GEMV2_X_BODY($1_type, $2_type, $3_type)', 
     `GEMV2_BODY($1_type, $2_type, $3_type, 
       SUM_TYPE($1_type, $2_type, $3_type), $1_type)')
   }')dnl
dnl
dnl
