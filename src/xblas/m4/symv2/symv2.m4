dnl -------------------------------------------------------------------
dnl  SYMV2 ---- Symmetric Matrix Vector Multiply 2
dnl
dnl    y  <---  alpha * A * (x_head + x_tail) + beta * y
dnl
dnl  where A is a symmetric matrix.
dnl -------------------------------------------------------------------
dnl
#include <blas_extended.h>
#include <blas_extended_private.h>
#include <blas_fpu.h>
include(cblas.m4)dnl
include(symv2-common.m4)dnl
dnl
dnl The call trees of the macros are as follows:
dnl
dnl     SYMV2  |-- SYMV2_HEAD
dnl            |      |
dnl            |      |-- SYMV2_NAME
dnl            |      |
dnl            |      |-- SYMV2_PARAMS
dnl            |
dnl            |-- SYMV2_COMMENT
dnl            |
dnl            |-- SYMV2_BODY
dnl
dnl          or (for _x routines)
dnl
dnl     SYMV2  |-- SYMV2_HEAD
dnl            |      |
dnl            |      |-- SYMV2_NAME
dnl            |      |
dnl            |      |-- SYMV2_PARAMS
dnl            |
dnl            |-- SYMV2_COMMENT
dnl            |
dnl            |-- SYMV2_X_BODY
dnl                   |
dnl                   |-- SWITCH_prec
dnl                          |
dnl                          |-- SYMV2_BODY
dnl
dnl
dnl
dnl
dnl  Usage:
dnl    SYMV2        ($1, $2, $3, $4)
dnl    SYMV2_HEAD   ($1, $2, $3, $4)
dnl    SYMV2_NAME   ($1, $2, $3, $4)
dnl    SYMV2_PARAMS ($1, $2, $3, $4)
dnl    SYMV2_COMMENT($1, $2, $3, $4)
dnl
dnl    $1 -- type of alpha, beta, y.
dnl    $2 -- type of A matrix
dnl    $3 -- type of x
dnl    $4 -- set to `_x' for _x routines.
dnl          Otherwise, `'.
dnl
define(`SYMV2_COMMENT', `
/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * (x_head + x_tail) + beta * y
 * 
 * where A is a symmetric matrix.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input symmetric matrix A.
 * 
 * uplo    (input) enum blas_uplo_type
 *         Determines which half of matrix A (upper or lower triangle)
 *         is accessed.
 *
 * n       (input) int
 *         Dimension of A and size of vectors x, y.
 *
 * alpha   (input) $1_scalar
 * 
 * a       (input) $2_array
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * x_head  (input) $3_array
 *         Vector x_head
 *
 * x_tail  (input) $3_array
 *         Vector x_tail
 *   
 * incx    (input) int
 *         Stride for vector x.
 *
 * beta    (input) $1_scalar
 * 
 * y       (input) $3_array
 *         Vector y.
 *
 * incy    (input) int
 *         Stride for vector y.
 *
PREC_COMMENT($4)dnl
 */')dnl
dnl
dnl
dnl ------------------------------------------------------------
dnl Usage:  SYMV2_BODY($1, $2, $3, $4, $5, $6)
dnl     $1 - type of alpha, beta, and vector y
dnl     $2 - type of matrix a
dnl     $3 - type of vector x
dnl     $4 - type of sum/prod
dnl     $5 - type of temp
dnl     $6 - [optional] String `FPU' is passed if FPU fix 
dnl          is needed.  Empty string is passed otherwise.
dnl ------------------------------------------------------------
dnl
define(`SYMV2_BODY', `
  int i, j;
  int xi, yi, xi0, yi0;
  int aij, ai;
  int incai;
  int incaij, incaij2;

  PTR_CAST(a, $2, `const')
  PTR_CAST(x_head, $3, `const')
  PTR_CAST(x_tail, $3, `const')
  PTR_CAST(y, $1, `')
  SCALAR_CAST(alpha, $1)
  SCALAR_CAST(beta,  $1)
  DECLARE(a_elem, $2)
  DECLARE(x_elem, $3)
  DECLARE(y_elem, $1)
  DECLARE(prod1, $4)
  DECLARE(prod2, $4)
  DECLARE(sum1, $4)
  DECLARE(sum2, $4)
  DECLARE(tmp1, $5)
  DECLARE(tmp2, $5)
  DECLARE(tmp3, $5)

  ifelse(`$6', `FPU', `FPU_FIX_DECL;')

  /* Test for no-op */
  if (n <= 0) {
    return;
  }
  if (TEST_0(alpha_i, $1) && TEST_1(beta_i, $1)) {
    return;
  }

  /* Check for error conditions. */
  if (n < 0) {
    BLAS_error(routine_name,  -3,  n, NULL);
  }
  if (lda < n) {
    BLAS_error(routine_name,  -6,  n, NULL);
  }
  if (incx == 0) {
    BLAS_error(routine_name,  -9,  incx, NULL);
  }
  if (incy == 0) {
    BLAS_error(routine_name,  -12,  incy, NULL);
  }

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai = lda;
    incaij = 1;
    incaij2 = lda;
  } else {
    incai = 1;
    incaij = lda;
    incaij2 = 1;
  }

  INC_ADJUST(incx, $3)
  INC_ADJUST(incy, $1)
  INC_ADJUST(incai, $2)
  INC_ADJUST(incaij, $2)
  INC_ADJUST(incaij2, $2)
  xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
  yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);

  ifelse(`$6', `FPU', `FPU_FIX_START;')

  /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
  for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
    ZERO(sum1, $4)
    ZERO(sum2, $4)
      
    for (j = 0, aij = ai, xi = xi0; j < i; j++, aij += incaij, xi += incx) {
      GET_VECTOR_ELEMENT(a_elem, a_i, aij, $2)
      GET_VECTOR_ELEMENT(x_elem, x_head_i, xi, $3)
      MUL(prod1, $4, a_elem, $2, x_elem, $3)
      ADD(sum1, $4, sum1, $4, prod1, $4)
      GET_VECTOR_ELEMENT(x_elem, x_tail_i, xi, $3)
      MUL(prod2, $4, a_elem, $2, x_elem, $3)
      ADD(sum2, $4, sum2, $4, prod2, $4)
    }
    for (; j < n; j++, aij += incaij2, xi += incx) {
      GET_VECTOR_ELEMENT(a_elem, a_i, aij, $2)
      GET_VECTOR_ELEMENT(x_elem, x_head_i, xi, $3)
      MUL(prod1, $4, a_elem, $2, x_elem, $3)
      ADD(sum1, $4, sum1, $4, prod1, $4)
      GET_VECTOR_ELEMENT(x_elem, x_tail_i, xi, $3)
      MUL(prod2, $4, a_elem, $2, x_elem, $3)
      ADD(sum2, $4, sum2, $4, prod2, $4)
    }
    ADD(sum1, $4, sum1, $4, sum2, $4)
    MUL(tmp1, $5, sum1, $4, alpha_i, $1)
    GET_VECTOR_ELEMENT(y_elem, y_i, yi, $1)
    MUL(tmp2, $5, y_elem, $1, beta_i, $1)
    ADD(tmp3, $5, tmp1, $5, tmp2, $5)
    SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp3, $5)
  }

  ifelse(`$6', `FPU', `FPU_FIX_STOP;')
')dnl
dnl
dnl Usage: SWITCH_prec($1, $2, $3, $4, $5, $6, $7, $8, $9)
dnl Generates a 4-way switch statement based on prec.
dnl Arguments
dnl    $1      --  type of vector y
dnl    $2      --  type of matrix A
dnl    $3      --  type of vector x
dnl    $4, $5  --  type of `sum' and `tmp' in single case
dnl    $6, $7  --  type of `sum' and `tmp' in double/indigenous case
dnl    $8, $9  --  type of `sum' and `tmp' in extra case
dnl
define(`SWITCH_prec', 
 `switch (prec) {

    case blas_prec_single: {
      SYMV2_BODY($1, $2, $3, $4, $5)
      break;
    }

    case blas_prec_double:
    case blas_prec_indigenous: {
      SYMV2_BODY($1, $2, $3, $6, $7)
      break;
    }

    case blas_prec_extra: { 
      SYMV2_BODY($1, $2, $3, $8, $9, FPU) 
      break;
    }
  }')dnl
dnl
dnl 
dnl ------------------------------------------------------------ 
dnl  Usage : SYMV2_X_BODY(y type, A type, x type)
dnl ------------------------------------------------------------
dnl 
define(`SYMV2_X_BODY', 
  `SWITCH_prec($1, $2, $3, 
    SUM_TYPE_X($1, $2, $3, S), TMP_TYPE_X($1, S), 
    SUM_TYPE_X($1, $2, $3, D), TMP_TYPE_X($1, D), 
    SUM_TYPE_X($1, $2, $3, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
dnl ------------------------------------------------------------
dnl  Usage: SYMV2(y type, A type, x type, extra)
dnl   where
dnl    y,A,x types are one of { s,d,c,z }
dnl    extra can be either `_x' or ommited.
dnl
dnl    The types represent single, double, complex, or double complex
dnl     precision and type for the respective arguments.
dnl     Note that y type also specifies the type of alpha, beta
dnl    if extra is set to `_x', then the precision parameter, prec,
dnl     is added to the parameter list of the function. 
dnl ------------------------------------------------------------
dnl
define(`SYMV2', 
  `SYMV2_HEAD($1, $2, $3, $4) 
   SYMV2_COMMENT($1, $2, $3, $4)
  {
  /* Routine name */
  const char routine_name[] = "SYMV2_NAME($1, $2, $3, $4)";
    ifelse($4, `_x', `SYMV2_X_BODY($1_type, $2_type, $3_type)', 
      `SYMV2_BODY($1_type, $2_type, $3_type, 
        SUM_TYPE($1_type, $2_type, $3_type), $1_type)')
  } /* end SYMV2_NAME($1, $2, $3, $4) */')dnl
dnl
dnl
