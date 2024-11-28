dnl ---------------------------------------------------------
dnl  HEMV ---- Hermitian Matrix Vector Multiply
dnl
dnl    y  <---   alpha * A * x + beta * y
dnl
dnl    where matrix A is hermitian.
dnl ---------------------------------------------------------
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(hemv-common.m4)dnl
dnl
dnl
dnl  Usage:
dnl    HEMV        ($1, $2, $3, $4)
dnl    HEMV_HEAD   ($1, $2, $3, $4)
dnl    HEMV_NAME   ($1, $2, $3, $4)
dnl    HEMV_PARAMS ($1, $2, $3, $4)
dnl    HEMV_COMMENT($1, $2, $3, $4)
dnl
dnl    $1 -- type of alpha, beta, y.
dnl    $2 -- type of A matrix
dnl    $3 -- type of x
dnl    $4 -- set to `_x' for _x routines.
dnl          Otherwise, `'.
dnl
dnl
define(`HEMV_COMMENT', `
/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * x  +  beta * y
 * 
 * where A is a Hermitian matrix.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input hermitian matrix A.
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
dnl  Special for loops for use in main body loops :
dnl   Usage :
dnl     $1 - type of matrix a
dnl     $2 - type of vector x
dnl     $3 - type of alpha, beta, and vector y
dnl     $4 - type of sum/prod
dnl     $5 - conj on first or second loop.
dnl ------------------------------------------------------------
define(`HEMV_FIRST_FOR_LOOP',
 `for (k = 0, aik = astarti, xi = x_starti; 
       k < i; k++, aik += incaik, xi += incx) {
    GET_VECTOR_ELEMENT(a_elem, a_i, aik, $1)
    ifelse($5, `conj_first_loop', `CONJ_AUX(a_elem, $1)', `')
    GET_VECTOR_ELEMENT(x_elem, x_i, xi, $2)
    MUL(prod, $4, a_elem, $1, x_elem, $2)
    ADD(sum, $4, sum, $4, prod, $4)
  }')dnl
dnl
dnl
define(`HEMV_DO_DIAGONAL_ELEMENT',
  `GET_VECTOR_ELEMENT(a_elem[0], a_i, aik, REAL_TYPE($1))
   GET_VECTOR_ELEMENT(x_elem, x_i, xi, $2)
   MUL(prod, $4, x_elem, $2, a_elem[0], REAL_TYPE($1))
   ADD(sum, $4, sum, $4, prod, $4)
   k++;
   aik += incaik2;
   xi += incx;')dnl
dnl
dnl
define(`HEMV_SECOND_FOR_LOOP',
  `for (; k < n_i; k++, aik += incaik2, xi += incx) {
     GET_VECTOR_ELEMENT(a_elem, a_i, aik, $1)
     ifelse($5, `conj_second_loop', `CONJ_AUX(a_elem, $1)', `')
     GET_VECTOR_ELEMENT(x_elem, x_i, xi, $2)
     MUL(prod, $4, a_elem, $1, x_elem, $2)
     ADD(sum, $4, sum, $4, prod, $4)
   }')dnl
dnl
dnl
dnl ------------------------------------------------------------
dnl     HEMV_MAIN_LOOPS
dnl  Usage :
dnl     $1 - type of matrix a
dnl     $2 - type of vector x
dnl     $3 - type of alpha, beta, and vector y
dnl     $4 - type of sum/prod
dnl     $5 - type of temp
dnl     $6 - string, either `conj_first_loop' or `conj_second_loop'
dnl ------------------------------------------------------------
define(`HEMV_MAIN_LOOPS', `
    /* Case alpha == 1. */
    if (TEST_1(alpha_i, $3)) {    
    
      if (TEST_0(beta_i, $3)) {
        /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
        for (i = 0, yi = y_starti, astarti = 0; 
             i < n_i; i++, yi += incy, astarti += incaik2) {
          ZERO(sum, $4)
          HEMV_FIRST_FOR_LOOP($1, $2, $3, $4, $6)
          HEMV_DO_DIAGONAL_ELEMENT($1, $2, $3, $4)
          HEMV_SECOND_FOR_LOOP($1, $2, $3, $4, $6)
          SET_ROUND_VECTOR_ELEMENT(y_i, yi, sum, $4)
        }
      } else {
        /* Case alpha = 1, but beta != 0. 
           We compute  y  <--- A * x + beta * y */
        for (i = 0, yi = y_starti, astarti = 0; 
            i < n_i; i++, yi += incy, astarti += incaik2) {
          ZERO(sum, $4)

          HEMV_FIRST_FOR_LOOP($1, $2, $3, $4, $6)
          HEMV_DO_DIAGONAL_ELEMENT($1, $2, $3, $4)
          HEMV_SECOND_FOR_LOOP($1, $2, $3, $4, $6)
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
                  i < n_i; i++, yi += incy, astarti += incaik2) {
        ZERO(sum, $4)
            
        HEMV_FIRST_FOR_LOOP($1, $2, $3, $4, $6)
        HEMV_DO_DIAGONAL_ELEMENT($1, $2, $3, $4)
        HEMV_SECOND_FOR_LOOP($1, $2, $3, $4, $6)
        GET_VECTOR_ELEMENT(y_elem, y_i, yi, $3)
        MUL(tmp2, $5, y_elem, $3, beta_i, $3)
        MUL(tmp1, $5, sum, $4, alpha_i, $3)
        ADD(tmp1, $5, tmp2, $5, tmp1, $5)
        SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)
      }
    }')dnl ---- end main loop
dnl
dnl
dnl ------------------------------------------------------------
dnl Usage:  HEMV_BODY($1, $2, $3, $4, $5, $6)
dnl     $1 - type of matrix a
dnl     $2 - type of vector x
dnl     $3 - type of alpha, beta, and vector y
dnl     $4 - type of sum/prod
dnl     $5 - type of temp
dnl     $6 - [optional] String `FPU' is passed if FPU fix 
dnl          is needed. 
dnl ------------------------------------------------------------
dnl
define(`HEMV_BODY', `
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
  if (order != blas_colmajor && order != blas_rowmajor) {
    BLAS_error(routine_name,  -1,  order, NULL);
  }
  if (uplo != blas_upper && uplo != blas_lower) {
    BLAS_error(routine_name,  -2,  uplo, NULL);
  }
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
  INC_ADJUST(incx, $2) dnl
  INC_ADJUST(incy, $3) dnl
  INC_ADJUST(incai, $1) dnl
  INC_ADJUST(incaik, $1) dnl
  INC_ADJUST(incaik2, $1) dnl
  if (incx < 0) {
    x_starti = (-n_i + 1) * incx;
  } else {      
    x_starti = 0;
  }
  if (incy < 0) {
    y_starti = (-n_i + 1) * incy;
  } else {
    y_starti = 0; 
  }
 
  ifelse(`$6', `FPU', `FPU_FIX_START;')

  /* alpha = 0.  In this case, just return beta * y */
  if (TEST_0(alpha_i, $3)) {
    for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
      GET_VECTOR_ELEMENT(y_elem, y_i, yi, $3)
      MUL(tmp1, $5, y_elem, $3, beta_i, $3)
      SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)
    }
  } else {
    /*  determine whether we conjugate in first loop or second loop */

    if (uplo == blas_lower) {
      /*  conjugate second*/
      HEMV_MAIN_LOOPS($1, $2, $3, $4, $5, conj_second_loop)
    } else {
      /*  conjugate first loop */
      HEMV_MAIN_LOOPS($1, $2, $3, $4, $5, conj_first_loop)
    }
  }
  ifelse(`$6', `FPU', `FPU_FIX_STOP;')
')dnl
dnl
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
define(`SWITCH_prec', `switch (prec) {

    case blas_prec_single: ifelse(`$4&&$5', `$6&&$7', `', `{
      HEMV_BODY($2, $3, $1, $4, $5)
      break;
    }
    ')dnl
    case blas_prec_indigenous: 
    case blas_prec_double: {
      HEMV_BODY($2, $3, $1, $6, $7)
      break;
    }

    case blas_prec_extra: { 
      HEMV_BODY($2, $3, $1, $8, $9, FPU) 
      break;
    }

    default:
      BLAS_error(routine_name,  -12,  prec, NULL);
      break;
  }')dnl
dnl
dnl 
dnl 
dnl ------------------------------------------------------------ 
dnl  Usage : HEMV_X_BODY(y type, A type, x type)
dnl ------------------------------------------------------------
dnl 
define(`HEMV_X_BODY', 
  `SWITCH_prec($1, $2, $3, 
    SUM_TYPE_X($1, $2, $2, S), TMP_TYPE_X($1, S), 
    SUM_TYPE_X($1, $2, $2, D), TMP_TYPE_X($1, D), 
    SUM_TYPE_X($1, $2, $2, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
define(`HEMV', 
  `HEMV_HEAD($1, $2, $3, $4) 
   HEMV_COMMENT($1, $2, $3, $4)
  {
  /* Routine name */
  static const char routine_name[] = "HEMV_NAME($1, $2, $3)";
    ifelse($4, `_x', `HEMV_X_BODY($1_type, $2_type, $3_type)', 
      `HEMV_BODY($2_type, $3_type, $1_type, 
       SUM_TYPE($1_type, $2_type, $3_type), $1_type)')
  } /* end HEMV_NAME($1, $2, $3) */')dnl
dnl
dnl
