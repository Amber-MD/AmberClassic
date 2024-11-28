dnl ----------------------------------------------------------------------
dnl AXPBY --- Scaled Vector Accumulation
dnl   y <--- alpha * x + beta * y
dnl ----------------------------------------------------------------------
dnl
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(axpby-common.m4)dnl
dnl
dnl
dnl  Usage:
dnl    AXPBY        ($1, $2, $3)
dnl    AXPBY_HEAD   ($1, $2, $3)
dnl    AXPBY_NAME   ($1, $2, $3)
dnl    AXPBY_PARAMS ($1, $2, $3)
dnl    AXPBY_COMMENT($1, $2, $3)
dnl
dnl      $1 -- type of alpha, beta, y
dnl      $2 -- type of x
dnl      $3 -- Set to `_x' for _x routines.  Otherwise set to `'.
dnl
dnl
define(`AXPBY_COMMENT', `
/*
 * Purpose
 * =======
 *
 * This routine computes:
 *
 *      y <- alpha * x + beta * y.
 *
 * Arguments
 * =========
 * 
 * n         (input) int
 *           The length of vectors x and y.
 * 
 * alpha     (input) $1_scalar
 *
 * x         (input) const $2_array
 *           Array of length n.
 *
 * incx      (input) int
 *           The stride used to access components x[i].
 * 
 * beta      (input) $1_scalar
 *
 * y         (input) $1_array
 *           Array of length n.
 * 
 * incy      (input) int
 *           The stride used to access components y[i].
 *
PREC_COMMENT($3)dnl
 */')dnl
dnl
dnl
dnl  Usage: AXPBY_BODY($1, $2, $3, $4)
dnl  Generates the main body of the product code.
dnl    $1 - type of alpha, beta, y
dnl    $2 - type of x
dnl    $3 - type of temp
dnl    $4 - [optional] String `FPU' is passed if FPU fix 
dnl         is needed.  Empty string is passed otherwise.
dnl
define(`AXPBY_BODY', `
  int i, ix = 0, iy = 0;
  PTR_CAST(x, $2, `const')
  PTR_CAST(y, $1)
  SCALAR_CAST(alpha, $1)
  SCALAR_CAST(beta, $1)
  DECLARE(x_ii, $2)
  DECLARE(y_ii, $1) 
  DECLARE(tmpx, $3)
  DECLARE(tmpy, $3)
  ifelse(`$4', `FPU', `FPU_FIX_DECL;')

  /* Test the input parameters. */
  if (incx == 0) 
    BLAS_error(routine_name, -4, incx, NULL);
  else if ( incy == 0 )
    BLAS_error(routine_name, -7, incy, NULL);

  /* Immediate return */
  if (n <= 0 || (TEST_0(alpha_i, $1) && TEST_1(beta_i, $1)))
    return;

  ifelse(`$4', `FPU', `FPU_FIX_START;')
  
  INC_ADJUST(incx, $2)
  INC_ADJUST(incy, $1)
  if (incx < 0) ix = (-n+1)*incx;
  if (incy < 0) iy = (-n+1)*incy;
  
  for (i = 0; i < n; ++i) {         
    GET_VECTOR_ELEMENT(x_ii, x_i, ix, $2)
    GET_VECTOR_ELEMENT(y_ii, y_i, iy, $1)
    MUL(tmpx, $3, alpha_i, $1, x_ii, $2) /* tmpx  = alpha * x[ix] */
    MUL(tmpy, $3, beta_i, $1, y_ii, $1) /* tmpy = beta * y[iy] */
    ADD(tmpy, $3, tmpy, $3, tmpx, $3)
    SET_ROUND_VECTOR_ELEMENT(y_i, iy, tmpy, $3)
    ix += incx;
    iy += incy; 
  } /* endfor */

  ifelse(`$4', `FPU', `FPU_FIX_STOP;')
')dnl
dnl
dnl 
dnl  Usage: SWITCH_prec($1, $2, $3, $4, $5)
dnl  Generates a 3-way switch statement based on prec.
dnl    $1  --  type of alpha, beta, y
dnl    $2  --  type of x
dnl    $3  --  type of tmp in single case
dnl    $4  --  type of tmp in double/indigenous case
dnl    $5  --  type of tmp in extra case
dnl
define(`SWITCH_prec', `
  switch (prec) {
  case blas_prec_single: ifelse(`$3', `$4', `', `{
    AXPBY_BODY($1, $2, $3)
    break;
  }
  ')dnl
  case blas_prec_double:
  case blas_prec_indigenous:
    { AXPBY_BODY($1, $2, $4) }
    break;
  case blas_prec_extra:
    { AXPBY_BODY($1, $2, $5, FPU) }
    break;
}')dnl
dnl
dnl
dnl  Usage: AXPBY_X_BODY($1, $2)
dnl  Generates the main body of the extended version of axpby code.
dnl    $1 -- type of alpha, beta, y
dnl    $2 -- type of x
dnl
define(`AXPBY_X_BODY', 
  `SWITCH_prec($1, $2, $1, 
    REAL_COMPLEX($1)`_'HIGHER_PREC(PREC($1), D), 
    REAL_COMPLEX($1)`_'HIGHER_PREC(PREC($1), E))')dnl
dnl
dnl
define(`AXPBY',
  `AXPBY_HEAD($1, $2, $3) 
   AXPBY_COMMENT($1, $2, $3)
   {
     static const char routine_name[] = "AXPBY_NAME($1, $2, $3)";
     ifelse($3, _x, `AXPBY_X_BODY($1_type, $2_type)', 
       `AXPBY_BODY($1_type, $2_type, $1_type)')
   } /* end AXPBY_NAME($1, $2, $3) */
')dnl
