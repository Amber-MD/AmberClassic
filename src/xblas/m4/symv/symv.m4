dnl ---------------------------------------------------------
dnl  SYMV ---- Symmetric Matrix Vector Multiply
dnl
dnl    y  <---   alpha * A * x + beta * y
dnl
dnl    where matrix A is symmetric.
dnl ---------------------------------------------------------
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(symv-common.m4)dnl
dnl
dnl
dnl  Usage:
dnl    SYMV        ($1, $2, $3, $4)
dnl    SYMV_HEAD   ($1, $2, $3, $4)
dnl    SYMV_NAME   ($1, $2, $3, $4)
dnl    SYMV_PARAMS ($1, $2, $3, $4)
dnl    SYMV_COMMENT($1, $2, $3, $4)
dnl
dnl    $1 -- type of alpha, beta, y.
dnl    $2 -- type of A matrix
dnl    $3 -- type of x
dnl    $4 -- set to `_x' for _x routines.
dnl          Otherwise, `'.
dnl
define(`SYMV_COMMENT', `
/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * x  +  beta * y
 * 
 * where A is a Symmetric matrix.
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
 * x       (input) $3_array
 *         Vector x.
 *   
 * incx    (input) int
 *         Stride for vector x.
 *
 * beta    (input) $1_scalar
 * 
 * y       (input/output) $1_array
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
dnl Usage:  SYMV_BODY($1, $2, $3, $4, $5, $6)
dnl     $1 - type of matrix a
dnl     $2 - type of vector x
dnl     $3 - type of alpha, beta, and vector y
dnl     $4 - type of sum/prod
dnl     $5 - type of temp
dnl     $6 - [optional] String `FPU' is passed if FPU fix 
dnl          is needed.  Empty string is passed otherwise.
dnl ------------------------------------------------------------
dnl
define(`SYMV_BODY', `
  /* Integer Index Variables */
  int i, k;

  int xi, yi;
  int aik, astarti, x_starti, y_starti;

  int incai;
  int incaik, incaik2;

  int n_i;

  /* Input Matrices */
  PTR_CAST(a, $1, `const')
  PTR_CAST(x, $2, `const')

  /* Output Vector */
  PTR_CAST(y, $3)

  /* Input Scalars */
  SCALAR_CAST(alpha, $3)
  SCALAR_CAST(beta, $3)

  /* Temporary Floating-Point Variables */
  DECLARE(a_elem, $1)
  DECLARE(x_elem, $2)
  DECLARE(y_elem, $3)
  DECLARE(prod, $4)
  DECLARE(sum, $4)
  DECLARE(tmp1, $5)
  DECLARE(tmp2, $5)

  ifelse(`$6', `FPU', `FPU_FIX_DECL;')

  /* Test for no-op */
  if (n <= 0) {
    return;
  }
  if (TEST_0(alpha_i, $3) && TEST_1(beta_i, $3)) {
    return;
  }

  /* Check for error conditions. */
  if (lda < n) {
    BLAS_error(routine_name,  -3,  n, NULL);
  }
  if (incx == 0) {
    BLAS_error(routine_name,  -8,  incx, NULL);
  }
  if (incy == 0) {
    BLAS_error(routine_name,  -11,  incy, NULL);
  }


  /* Set Index Parameters */
  n_i = n;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai = lda;
    incaik = 1;
    incaik2 = lda;
  } else {
    incai = 1;
    incaik = lda;
    incaik2 = 1;
  }

  /* Adjustment to increments (if any) */
  INC_ADJUST(incx, $2)
  INC_ADJUST(incy, $3)
  INC_ADJUST(incai, $1)
  INC_ADJUST(incaik, $1)
  INC_ADJUST(incaik2, $1)
  if (incx < 0) {
    x_starti = (-n + 1) * incx;
  }
  else {        
    x_starti = 0;
  }
  if (incy < 0) {
    y_starti = (-n + 1) * incy;
  } else {
    y_starti = 0; 
  } 

  ifelse(`$6', `FPU', `FPU_FIX_START;')

  /* alpha = 0.  In this case, just return beta * y */
  if (TEST_0(alpha_i, $3)) {
    for (i = 0, yi = y_starti; 
        i < n_i; i++, yi += incy) {
      GET_VECTOR_ELEMENT(y_elem, y_i, yi, $3)
      MUL(tmp1, $5, y_elem, $3, beta_i, $3)
      SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)
    }
  } else if (TEST_1(alpha_i, $3)) {
    
    /* Case alpha == 1. */
    
    if (TEST_0(beta_i, $3)) {
      /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
      for (i = 0, yi = y_starti, astarti = 0; 
        i < n_i; i++, yi += incy, astarti += incai) {
        ZERO(sum, $4)
          
        for (k = 0, aik = astarti, xi = x_starti; k < i; 
                k++, aik += incaik, xi += incx) {
          GET_VECTOR_ELEMENT(a_elem, a_i, aik, $1)
          GET_VECTOR_ELEMENT(x_elem, x_i, xi, $2)
          MUL(prod, $4, a_elem, $1, x_elem, $2)
          ADD(sum, $4, sum, $4, prod, $4)
        }
        for (; k < n_i; k++, aik += incaik2, xi += incx) {
          GET_VECTOR_ELEMENT(a_elem, a_i, aik, $1)
          GET_VECTOR_ELEMENT(x_elem, x_i, xi, $2)
          MUL(prod, $4, a_elem, $1, x_elem, $2)
          ADD(sum, $4, sum, $4, prod, $4)
        }
        ifelse(`$5', `$4', 
          `SET_ROUND_VECTOR_ELEMENT(y_i, yi, sum, $4)',
          `ASSIGN(tmp1, $5, sum, $4)
           SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)')
      }
    } else {
      /* Case alpha = 1, but beta != 0. 
         We compute  y  <--- A * x + beta * y */
      for (i = 0, yi = y_starti, astarti = 0; 
        i < n_i; i++, yi += incy, astarti += incai) {
        ZERO(sum, $4)
          
        for (k = 0, aik = astarti, xi = x_starti; 
                k < i; k++, aik += incaik, xi += incx) {
          GET_VECTOR_ELEMENT(a_elem, a_i, aik, $1)
          GET_VECTOR_ELEMENT(x_elem, x_i, xi, $2)
          MUL(prod, $4, a_elem, $1, x_elem, $2)
          ADD(sum, $4, sum, $4, prod, $4)
        }
        for (; k < n_i; k++, aik += incaik2, xi += incx) {
          GET_VECTOR_ELEMENT(a_elem, a_i, aik, $1)
          GET_VECTOR_ELEMENT(x_elem, x_i, xi, $2)
          MUL(prod, $4, a_elem, $1, x_elem, $2)
          ADD(sum, $4, sum, $4, prod, $4)
        }
        GET_VECTOR_ELEMENT(y_elem, y_i, yi, $3)
        MUL(tmp2, $5, y_elem, $3, beta_i, $3)
        ASSIGN(tmp1, $5, sum, $4)
        ADD(tmp1, $5, tmp2, $5, tmp1, $5)
        SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)
      }
    }
    
  } else {
    /* The most general form,   y <--- alpha * A * x + beta * y */
      for (i = 0, yi = y_starti, astarti = 0; 
        i < n_i; i++, yi += incy, astarti += incai) {
        ZERO(sum, $4)
          
        for (k = 0, aik = astarti, xi = x_starti; 
                k < i; k++, aik += incaik, xi += incx) {
          GET_VECTOR_ELEMENT(a_elem, a_i, aik, $1)
          GET_VECTOR_ELEMENT(x_elem, x_i, xi, $2)
          MUL(prod, $4, a_elem, $1, x_elem, $2)
          ADD(sum, $4, sum, $4, prod, $4)
        }
        for (; k < n_i; k++, aik += incaik2, xi += incx) {
          GET_VECTOR_ELEMENT(a_elem, a_i, aik, $1)
          GET_VECTOR_ELEMENT(x_elem, x_i, xi, $2)
          MUL(prod, $4, a_elem, $1, x_elem, $2)
          ADD(sum, $4, sum, $4, prod, $4)
        }
        GET_VECTOR_ELEMENT(y_elem, y_i, yi, $3)
        MUL(tmp2, $5, y_elem, $3, beta_i, $3)
        MUL(tmp1, $5, sum, $4, alpha_i, $3)
        ADD(tmp1, $5, tmp2, $5, tmp1, $5)
        SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)
      }
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

    case blas_prec_single: ifelse(`$4&&$5', `$6&&$7', `', `{
      SYMV_BODY($2, $3, $1, $4, $5)
      break;
    }
    ')dnl
    case blas_prec_double:
    case blas_prec_indigenous: {
      SYMV_BODY($2, $3, $1, $6, $7)
      break;
    }

    case blas_prec_extra: { 
      SYMV_BODY($2, $3, $1, $8, $9, FPU) 
      break;
    }
  }')dnl
dnl
dnl 
dnl ------------------------------------------------------------ 
dnl  Usage : SYMV_X_BODY(y type, A type, x type)
dnl ------------------------------------------------------------
dnl 
define(`SYMV_X_BODY', 
  `SWITCH_prec($1, $2, $3, 
    SUM_TYPE_X($1, $2, $3, S), TMP_TYPE_X($1, S), 
    SUM_TYPE_X($1, $2, $3, D), TMP_TYPE_X($1, D), 
    SUM_TYPE_X($1, $2, $3, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
dnl ------------------------------------------------------------
dnl  Usage: SYMV(y type, A type, x type, extra)
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
define(`SYMV', 
  `SYMV_HEAD($1, $2, $3, $4) 
   SYMV_COMMENT($1, $2, $3, $4)
  {
  /* Routine name */
  static const char routine_name[] = "SYMV_NAME($1, $2, $3, $4)";
    ifelse($4, `_x', `SYMV_X_BODY($1_type, $2_type, $3_type)', 
      `SYMV_BODY($2_type, $3_type, $1_type, 
        SUM_TYPE($1_type, $2_type, $3_type), $1_type)')
  } /* end SYMV_NAME($1, $2, $3, $4) */')dnl
dnl
dnl
