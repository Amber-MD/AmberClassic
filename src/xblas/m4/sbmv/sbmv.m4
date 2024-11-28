dnl ---------------------------------------------------------
dnl  SBMV ---- Symmetric Band Matrix Vector Multiply
dnl
dnl    y  <---   alpha * A * x + beta * y
dnl
dnl    where matrix A is a symmetric band matrix.
dnl ---------------------------------------------------------
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(sbmv-common.m4)dnl
dnl
dnl
dnl  Usage:
dnl    SBMV        ($1, $2, $3, $4)
dnl    SBMV_HEAD   ($1, $2, $3, $4)
dnl    SBMV_NAME   ($1, $2, $3, $4)
dnl    SBMV_PARAMS ($1, $2, $3, $4)
dnl    SBMV_COMMENT($1, $2, $3, $4)
dnl
dnl    $1 -- type of alpha, beta, y.
dnl    $2 -- type of A matrix
dnl    $3 -- type of x
dnl    $4 -- set to `_x' for _x routines.
dnl          Otherwise, `'.
dnl
dnl
define(`SBMV_COMMENT', `
/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * x  +  beta * y
 * 
 * where A is a symmetric band matrix.
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
 * k       (input) int
 *         Number of subdiagonals ( = number of superdiagonals)
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
 *
 *  Notes on storing a symmetric band matrix:
 * 
 *      Integers in the below arrays represent values of
 *              type ifelse(`$2',`s',`float',
                        `$2', `d', double,
                        `$2', `c', `complex float',
                        `$2', `z', `complex double').
 *
 *    if we have a symettric matrix:
 *
 *      1  2  3  0  0
 *      2  4  5  6  0
 *      3  5  7  8  9
 *      0  6  8  10 11
 *      0  0  9  11 12
 *
 *     This matrix has n == 5, and k == 2. It can be stored in the
 *      following ways:
 *
 *      Notes for the examples:
 *      Each column below represents a contigous vector.
 *      Columns are strided by lda.
 *      An asterisk (*) represents a position in the 
 *       matrix that is not used.
 *      Note that the minimum lda (size of column) is 3 (k+1).
 *       lda may be arbitrarily large; an lda > 3 would mean
 *       there would be unused data at the bottom of the below
 *       columns.        
 *
 *    blas_colmajor and blas_upper:
 *      *  *  3  6  9  
 *      *  2  5  8  11 
 *      1  4  7  10 12
 *
 *
 *    blas_colmajor and blas_lower
 *      1  4  7  10  12
 *      2  5  8  11  *
 *      3  6  9  *   *
 *
 *
 *    blas_rowmajor and blas_upper 
 *      Columns here also represent contigous arrays.
 *      1  4  7  10  12
 *      2  5  8  11  *
 *      3  6  9  *   *
 *
 *
 *    blas_rowmajor and blas_lower
 *      Columns here also represent contigous arrays.
 *      *  *  3  6   9
 *      *  2  5  8   11
 *      1  4  7  10  12
 *
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
dnl ------------------------------------------------------------
define(`SBMV_FIRST_FOR_LOOP',
 `for (j = 0, aij = astarti, xi = x_starti; 
       j < maxj_first; j++, aij += incaij, xi += incx) {
    GET_VECTOR_ELEMENT(a_elem, a_i, aij, $1)
    GET_VECTOR_ELEMENT(x_elem, x_i, xi, $2)
    MUL(prod, $4, a_elem, $1, x_elem, $2)
    ADD(sum, $4, sum, $4, prod, $4)
  }')dnl
dnl
dnl
define(`SBMV_SECOND_FOR_LOOP',
  `for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
     GET_VECTOR_ELEMENT(a_elem, a_i, aij, $1)
     GET_VECTOR_ELEMENT(x_elem, x_i, xi, $2)
     MUL(prod, $4, a_elem, $1, x_elem, $2)
     ADD(sum, $4, sum, $4, prod, $4)
   }')dnl
dnl
dnl
define(`SBMV_SET_NEW_COUNTERS',
  `if (i+1 >= (n_i - k))
     maxj_second--;
   if (i >= k) {
     astarti += (incaij + incaij2);
     x_starti += incx;
   } else {
     maxj_first++;
     astarti += incaij2;
   }')dnl
dnl
dnl
dnl ------------------------------------------------------------
dnl     SBMV_MAIN_LOOPS
dnl  Usage :
dnl     $1 - type of matrix a
dnl     $2 - type of vector x
dnl     $3 - type of alpha, beta, and vector y
dnl     $4 - type of sum/prod
dnl     $5 - type of temp
dnl ------------------------------------------------------------
define(`SBMV_MAIN_LOOPS', `
    /* Case alpha == 1. */
    if (TEST_1(alpha_i, $3)) {    
    
      if (TEST_0(beta_i, $3)) {
        /* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
        for (i = 0, yi = y_starti; 
             i < n_i; i++, yi += incy) {
          ZERO(sum, $4)
          SBMV_FIRST_FOR_LOOP($1, $2, $3, $4)
          SBMV_SECOND_FOR_LOOP($1, $2, $3, $4)
          ifelse(`$5', `$4', 
            `SET_ROUND_VECTOR_ELEMENT(y_i, yi, sum, $4)',
            `ASSIGN(tmp1, $5, sum, $4)
             SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)')
          SBMV_SET_NEW_COUNTERS()
        }
      } else {
        /* Case alpha = 1, but beta != 0. 
           We compute  y  <--- A * x + beta * y */
        for (i = 0, yi = y_starti; 
            i < n_i; i++, yi += incy) {
          ZERO(sum, $4)

          SBMV_FIRST_FOR_LOOP($1, $2, $3, $4)
          SBMV_SECOND_FOR_LOOP($1, $2, $3, $4)
          GET_VECTOR_ELEMENT(y_elem, y_i, yi, $3)
          MUL(tmp2, $5, y_elem, $3, beta_i, $3)
          ASSIGN(tmp1, $5, sum, $4)
          ADD(tmp1, $5, tmp2, $5, tmp1, $5)
          SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)
          SBMV_SET_NEW_COUNTERS()
        }
      }
    } else {
      if (TEST_0(beta_i, $3)) {
        /* Case alpha != 1, but beta == 0. 
           We compute  y  <--- A * x * a */
        for (i = 0, yi = y_starti; 
            i < n_i; i++, yi += incy) {
          ZERO(sum, $4)

          SBMV_FIRST_FOR_LOOP($1, $2, $3, $4)
          SBMV_SECOND_FOR_LOOP($1, $2, $3, $4)
          GET_VECTOR_ELEMENT(y_elem, y_i, yi, $3)
          MUL(tmp1, $5, sum, $4, alpha_i, $3)
          SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)
          SBMV_SET_NEW_COUNTERS()
        }
      } else {
      /* The most general form,   y <--- alpha * A * x + beta * y */
        for (i = 0, yi = y_starti; 
                  i < n_i; i++, yi += incy) {
          ZERO(sum, $4)
            
          SBMV_FIRST_FOR_LOOP($1, $2, $3, $4)
          SBMV_SECOND_FOR_LOOP($1, $2, $3, $4)
          GET_VECTOR_ELEMENT(y_elem, y_i, yi, $3)
          MUL(tmp2, $5, y_elem, $3, beta_i, $3)
          MUL(tmp1, $5, sum, $4, alpha_i, $3)
          ADD(tmp1, $5, tmp2, $5, tmp1, $5)
          SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)
          SBMV_SET_NEW_COUNTERS()
        }
      }
    }')dnl ---- end main loop
dnl ------------------------------------------------------------
dnl Usage:  SBMV_BODY($1, $2, $3, $4, $5, $6)
dnl     $1 - type of matrix a
dnl     $2 - type of vector x
dnl     $3 - type of alpha, beta, and vector y
dnl     $4 - type of sum/prod
dnl     $5 - type of temp
dnl     $6 - [optional] String `FPU' is passed if FPU fix 
dnl          is needed.  
dnl ------------------------------------------------------------
dnl
define(`SBMV_BODY', `
  /* Integer Index Variables */
  int i, j;
  int xi, yi;
  int aij, astarti, x_starti, y_starti;
  int incaij, incaij2;
  int n_i;
  int maxj_first, maxj_second;
 
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
    BLAS_error(routine_name, -1, order, 0);
  }
  if (uplo != blas_upper && uplo != blas_lower) {
    BLAS_error(routine_name, -2, uplo, 0);
  }
  if (n < 0) {
    BLAS_error(routine_name, -3, n, 0);
  }
  if (k < 0 || k > n) {
    BLAS_error(routine_name, -4, k, 0);
  }
  if ((lda < k + 1) || (lda < 1)) {
    BLAS_error(routine_name, -7, lda, 0);
  }
  if (incx == 0) {
    BLAS_error(routine_name, -9, incx, 0);
  }
  if (incy == 0) {
    BLAS_error(routine_name, -12, incy, 0);
  }

  /* Set Index Parameters */
  n_i = n;

  if ( ((uplo == blas_upper) && (order == blas_colmajor)) ||
        ((uplo == blas_lower) && (order == blas_rowmajor)) ) {
    incaij = 1;             /* increment in first loop */
    incaij2 = lda - 1;      /* increment in second loop */
    astarti = k;    /* does not start on zero element */
  } else {
    incaij = lda - 1;
    incaij2 = 1;
    astarti = 0; /* start on first element of array */
  }
  /* Adjustment to increments (if any) */
  INC_ADJUST(incx, $2) dnl
  INC_ADJUST(incy, $3) dnl
  INC_ADJUST(astarti, $1) dnl
  INC_ADJUST(incaij, $1) dnl
  INC_ADJUST(incaij2, $1) dnl
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
    for (i = 0, yi = y_starti; 
         i < n_i; i++, yi += incy) {
      GET_VECTOR_ELEMENT(y_elem, y_i, yi, $3)
      MUL(tmp1, $5, y_elem, $3, beta_i, $3)
      SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)
    }
  } else {
    /*  determine the loop interation counts */
    /* number of elements done in first loop 
       (this will increase by one over each column up to a limit) */
    maxj_first = 0;
    /* number of elements done in second loop the first time*/
    maxj_second = MIN(k+1, n_i); 
    SBMV_MAIN_LOOPS($1, $2, $3, $4, $5)
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
define(`SWITCH_prec', `switch (prec) {

    case blas_prec_single: ifelse(`$4&&$5', `$6&&$7', `', `{
      SBMV_BODY($2, $3, $1, $4, $5)
      break;
    }
    ')dnl
    case blas_prec_indigenous: 
    case blas_prec_double: {
      SBMV_BODY($2, $3, $1, $6, $7)
      break;
    }

    case blas_prec_extra: { 
      SBMV_BODY($2, $3, $1, $8, $9, FPU) 
      break;
    }

    default:
      BLAS_error(routine_name, -13, prec, 0);
      break;
  }')dnl
dnl
dnl 
define(`SBMV_X_BODY', 
  `SWITCH_prec($1, $2, $3, 
    SUM_TYPE_X($1, $2, $3, S), TMP_TYPE_X($1, S), 
    SUM_TYPE_X($1, $2, $3, D), TMP_TYPE_X($1, D), 
    SUM_TYPE_X($1, $2, $3, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
define(`SBMV', 
  `SBMV_HEAD($1, $2, $3, $4) 
   SBMV_COMMENT($1, $2, $3, $4)
   {
     static const char routine_name[] = "SBMV_NAME($1, $2, $3, $4)";
     ifelse($4, `_x', `SBMV_X_BODY($1_type, $2_type, $3_type)', 
       `SBMV_BODY($2_type, $3_type, $1_type, 
         SUM_TYPE($1_type, $2_type, $3_type), $1_type)')
   } /* end SBMV_NAME($1, $2, $3) */')dnl
dnl
dnl
