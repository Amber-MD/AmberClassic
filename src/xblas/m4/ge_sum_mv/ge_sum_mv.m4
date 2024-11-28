dnl ----------------------------------------------------------------------
dnl GE_SUM_MV -- Summed Matrix Vector Multiplies
dnl   y <--- alpha * A * x + beta * B * x
dnl ----------------------------------------------------------------------
dnl
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(ge_sum_mv-common.m4)dnl
dnl
dnl 
dnl  Usage:
dnl    GE_SUM_MV        ($1, $2, $3, $4)
dnl    GE_SUM_MV_HEAD   ($1, $2, $3, $4)
dnl    GE_SUM_MV_NAME   ($1, $2, $3, $4)
dnl    GE_SUM_MV_PARAMS ($1, $2, $3, $4)
dnl    GE_SUM_MV_COMMENT($1, $2, $3, $4)
dnl
dnl    $1 -- type of alpha, beta, y.
dnl    $2 -- type of A matrix
dnl    $3 -- type of x
dnl    $4 -- set to `_x' for _x routines.  Otherwise set to `'.
dnl
dnl
define(`GE_SUM_MV_COMMENT', `
/*
 * Purpose
 * =======
 *
 * Computes y = alpha * A * x + beta * B * y, 
 *     where A, B are general matricies.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * m            (input) int
 *              Row Dimension of A, B, length of output vector y
 *
 * n            (input) int
 *              Column Dimension of A, B and the length of vector x
 *
 * alpha        (input) $1_scalar
 *              
 * A            (input) const $2_array
 *
 * lda          (input) int 
 *              Leading dimension of A
 *
 * x            (input) const $3_array
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) $1_scalar
 *
 * b            (input) const $2_array
 *
 * ldb          (input) int 
 *              Leading dimension of B
 *
 * y            (input/output) $1_array
 *
 * incy         (input) int
 *              The stride for vector y.
 * 
PREC_COMMENT($4)dnl
 */')dnl
dnl
dnl
define(`GE_SUM_MV_BODY',
 `int i, j;
  int xi, yi;
  int x_starti, y_starti, incxi, incyi;
  int lda_min;
  int ai;
  int incai;
  int aij;
  int incaij;
  int bi;
  int incbi;
  int bij;
  int incbij;

  PTR_CAST(a, $2, `const')
  PTR_CAST(b, $2, `const')
  PTR_CAST(x, $3, `const')
  PTR_CAST(y, $1)
  SCALAR_CAST(alpha, $1)
  SCALAR_CAST(beta, $1)
  DECLARE(a_elem, $2)
  DECLARE(b_elem, $2)
  DECLARE(x_elem, $3)
  DECLARE(prod, $4)
  DECLARE(sumA, $4)
  DECLARE(sumB, $4)
  DECLARE(tmp1, $5)
  DECLARE(tmp2, $5)

  ifelse(`$6', `FPU', `FPU_FIX_DECL;')

  /* m is number of rows */
  /* n is number of columns */

  if (m == 0 || n == 0)
    return; 


  /* all error calls */
  if (order == blas_rowmajor) {
    lda_min = n;
    incai = lda; /* row stride */
    incbi = ldb; 
    incbij = incaij = 1;  /* column stride */
  } else if (order == blas_colmajor) {
    lda_min = m;
    incai = incbi = 1; /*row stride */
    incaij = lda;  /* column stride */
    incbij = ldb;
  } else {
    /* error, order not blas_colmajor not blas_rowmajor */
    BLAS_error(routine_name, -1, order, 0);
    return;
  }

  if (m < 0)
    BLAS_error(routine_name, -2, m, 0);
  else if (n < 0)
    BLAS_error(routine_name, -3, n, 0);
  if (lda < lda_min)
    BLAS_error(routine_name, -6, lda, 0);
  else if (ldb < lda_min)
    BLAS_error(routine_name, -11, ldb, 0);
  else if (incx == 0)
    BLAS_error(routine_name, -8, incx, 0);
  else if (incy == 0)
    BLAS_error(routine_name, -13, incy, 0);

  incxi = incx;
  incyi = incy;
  INC_ADJUST(incxi, $3)
  INC_ADJUST(incyi, $1)
  INC_ADJUST(incai, $2)
  INC_ADJUST(incaij, $2)
  INC_ADJUST(incbi, $2)
  INC_ADJUST(incbij, $2)

  if (incxi > 0)
    x_starti = 0;
  else
    x_starti = (1 - n)*incxi;

  if (incyi > 0)
    y_starti = 0;
  else
    y_starti = (1 - m)*incyi;

  ifelse(`$6', `FPU', `FPU_FIX_START;')

  if (TEST_0(alpha_i, $1)) {
    if (TEST_0(beta_i, $1)) {
      /* alpha, beta are 0.0 */
      for ( i=0, yi = y_starti; i < m; i++, yi += incyi ) {
        SET_ZERO_VECTOR_ELEMENT(y_i, yi, $1)
      } 
    } else if (TEST_1(beta_i, $1)) {
      /* alpha is 0.0, beta is 1.0 */
      GE_SUM_MV_LOOPS($1, $2, $3, $4, $5, `alpha_zero', `beta_one')
    } else {
      /* alpha is 0.0, beta not 1.0 nor 0.0 */
      GE_SUM_MV_LOOPS($1, $2, $3, $4, $5, `alpha_zero', `')
    }
  } else if(TEST_1(alpha_i, $1)) {
    if (TEST_0(beta_i, $1)) {
      /* alpha is 1.0, beta is 0.0 */
      GE_SUM_MV_LOOPS($1, $2, $3, $4, $5, `alpha_one', `beta_zero')
    } else if (TEST_1(beta_i, $1)) {
      /* alpha is 1.0, beta is 1.0 */
      GE_SUM_MV_LOOPS($1, $2, $3, $4, $5, `alpha_one', `beta_one')
    } else {
      /* alpha is 1.0, beta is other */
      GE_SUM_MV_LOOPS($1, $2, $3, $4, $5, `alpha_one', `')
    }
  } else {
    if (TEST_0(beta_i, $1)) {
      /* alpha is other, beta is 0.0 */
      GE_SUM_MV_LOOPS($1, $2, $3, $4, $5, `', `beta_zero')
    } else if (TEST_1(beta_i, $1)) {
      /* alpha is other, beta is 1.0 */
      GE_SUM_MV_LOOPS($1, $2, $3, $4, $5, `', `beta_one')
    } else {
      /* most general form, alpha, beta are other */
      GE_SUM_MV_LOOPS($1, $2, $3, $4, $5, `', `')
    }
  }
  ifelse(`$6', `FPU', `FPU_FIX_STOP;')
')dnl
dnl
dnl
define(`GE_SUM_MV_LOOPS', `
  ifelse($6, `alpha_zero', `', `ai = 0;')
  ifelse($7, `beta_zero',  `', `bi = 0;')
  for (i = 0, yi = y_starti; i < m; i++, yi += incyi) {
    ifelse($6, `alpha_zero', `', 
      `ZERO(sumA, $4)
      aij = ai;')
    ifelse($7, `beta_zero', `', 
      `ZERO(sumB, $4)
      bij = bi;')
    for (j = 0, xi = x_starti; j < n; j++, xi += incxi) {
      GET_VECTOR_ELEMENT(x_elem, x_i, xi, $3)
      ifelse($6, `alpha_zero', `',
        `GET_VECTOR_ELEMENT(a_elem, a_i, aij, $2)
          MUL(prod, $4, a_elem, $2, x_elem, $3)
         ADD(sumA, $4, sumA, $4, prod, $4)
        aij += incaij;')
      ifelse($7, `beta_zero', `',
        `GET_VECTOR_ELEMENT(b_elem, b_i, bij, $2)
          MUL(prod, $4, b_elem, $2, x_elem, $3)
         ADD(sumB, $4, sumB, $4, prod, $4)
        bij += incbij;')
    }
    /* now put the result into y_i */
    ifelse($6, `alpha_zero',
      `ifelse($7, `beta_one',
       `ifelse($5, $4,
          `SET_ROUND_VECTOR_ELEMENT(y_i, yi, sumB, $4)',
          `ASSIGN(tmp1, $5, sumB, $4)
            SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)')',
       `MUL(tmp1, $5, sumB, $4, beta_i, $1)
        SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)')',
      `$7', `beta_zero',
      `ifelse($6, `alpha_one',
       `ifelse($5, $4,
          `SET_ROUND_VECTOR_ELEMENT(y_i, yi, sumA, $4)',
          `ASSIGN(tmp1, $5, sumA, $4)
            SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)')',
       `MUL(tmp1, $5, sumA, $4, alpha_i, $1)
        SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)')',
      `    dnl most general case ---
      ifelse($6, `alpha_one', 
        `ASSIGN(tmp1, $5, sumA, $4)',
        `MUL(tmp1, $5, sumA, $4, alpha_i, $1)')
      ifelse($7, `beta_one', 
        `ASSIGN(tmp2, $5, sumB, $4)',
        `MUL(tmp2, $5, sumB, $4, beta_i, $1)')
      ADD(tmp1, $5, tmp1, $5, tmp2, $5)
      SET_ROUND_VECTOR_ELEMENT(y_i, yi, tmp1, $5)')
    ifelse($6, `alpha_zero', `', `ai += incai;')
    ifelse($7, `beta_zero', `',  `bi += incbi;')
  }')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SWITCH_prec($1, $2, $3, $4, $5, $6, $7, $8, $9) ... generate
dnl        a 3-way switch statement based on prec.
dnl        $4 and $5 are the types of 'sum' and 'tmp' in single case.
dnl        $6 and $7 are the types of 'sum' and 'tmp' in double/indigenous case.
dnl        $7 and $8 are the types of 'sum' and 'tmp' in extra case.
dnl ----------------------------------------------------------------------
define(`SWITCH_prec', 
 `switch ( prec ) {
    case blas_prec_single: ifelse(`$4&&$5', `$6&&$7', `', `{
      GE_SUM_MV_BODY($1, $2, $3, $4, $5)
      break;
    }
    ')dnl
    case blas_prec_indigenous:
    case blas_prec_double: 
    { GE_SUM_MV_BODY($1, $2, $3, $6, $7) }
    break;

    case blas_prec_extra: 
    { GE_SUM_MV_BODY($1, $2, $3, $8, $9, FPU) }
    break;

    default: 
    { BLAS_error(routine_name, -14, prec, 0); }
    break;
  }
')dnl
dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: GE_SUM_MV_X_BODY(aby_type, A_type, x_type) ... dispatches
dnl        GE_SUM_MV with appropriate type and precision info of
dnl        the specified internal prec.       
dnl Each type specifier can be one of
dnl        s   ... real and single
dnl     d   ... real and double
dnl     c  ... complex and single
dnl     z  ... complex and double
dnl ----------------------------------------------------------------------
define(`GE_SUM_MV_X_BODY', 
  `SWITCH_prec($1, $2, $3, 
    SUM_TYPE_X($1, $2, $3, S), TMP_TYPE_X($1, S), 
    SUM_TYPE_X($1, $2, $3, D), TMP_TYPE_X($1, D), 
    SUM_TYPE_X($1, $2, $3, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
define(`GE_SUM_MV', 
  `GE_SUM_MV_HEAD($1, $2, $3, $4) 
   GE_SUM_MV_COMMENT($1, $2, $3, $4)
  {
    /* Routine name */
    static const char routine_name[] = "GE_SUM_MV_NAME($1, $2, $3)";
    ifelse($4, `_x', `GE_SUM_MV_X_BODY($1_type, $2_type, $3_type)', 
      `GE_SUM_MV_BODY($1_type, $2_type, $3_type, 
        SUM_TYPE($1_type, $2_type, $3_type), $1_type)')
  } /* end GE_SUM_MV_NAME($1, $2, $3) */')dnl
dnl
dnl
