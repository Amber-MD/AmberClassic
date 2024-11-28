dnl **********************************************************************
dnl * Perform  sum <- SUM_{i=0, n-1} x[i]                                *
dnl **********************************************************************
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(sum-common.m4)dnl
dnl
dnl
define(`SUM_COMMENT',`
/*
 * Purpose
 * =======
 * 
 * This routine computes the summation:
 * 
 *     sum <- SUM_{i=0, n-1} x[i].
 * 
 * Arguments
 * =========
 *
 * n      (input) int
 *        The length of vector x.
 * 
 * x      (input) const $1_array
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * sum    (output) $1_array
 * 
PREC_COMMENT($2)dnl
 */')dnl
dnl
dnl
define(`SUM', 
  `SUM_HEAD($1, $2)
   SUM_COMMENT($1, $2) {
     static const char routine_name[] = "SUM_NAME($1, $2)";
     ifelse(`$2', `_x', `SUM_X_BODY($1_type, $2_type)',
       `SUM_BODY($1_type, $1_type)')
   }')dnl
dnl
dnl
dnl Usage:  SUM($1, $2, [$3])
dnl
dnl    $1 -- type of x (input vector).
dnl    $2 -- type of accumulator tmp
dnl    $3 -- [optional] String `FPU' passed if FPU fix
dnl          is needed.  Otherwise, it is an empty string.
dnl
define(`SUM_BODY', `
  int i, xi;
  PTR_CAST(sum, $1)
  PTR_CAST(x, $1, `const')
  DECLARE(x_elem, $1)
  DECLARE(tmp, $2)
  ifelse(`$3', `FPU', `FPU_FIX_DECL;')

  /* Test the input parameters. */
  if ( n < 0 )
    BLAS_error(routine_name,  -1,  n, NULL);
  if ( incx == 0 )
    BLAS_error(routine_name,  -3,  incx, NULL);

  /* Immediate return. */
  if ( n <= 0 ) {
    ZERO_OUT(sum_i, $1)
    return;
  }

  ifelse(`$3', `FPU', `FPU_FIX_START;')

  ZERO(tmp, $2)

  INC_ADJUST(incx, $1)
  if (incx < 0) 
    xi = -(n-1)*incx;
  else
    xi = 0;

  for (i = 0; i < n; i++, xi += incx) {
    GET_VECTOR_ELEMENT(x_elem, x_i, xi, $1)
    ADD(tmp, $2, tmp, $2, x_elem, $1) 
  }
  ROUND(sum, $1, tmp, $2)

  ifelse(`$3', `FPU', `FPU_FIX_STOP;')
')dnl
dnl
dnl
dnl Usage: SWITCH_prec($1, $2, $3, $4)
dnl   $1 -- type of sum / input vector
dnl   $2 -- type of temp in single case
dnl   $3 -- type of temp in double/indigenous case
dnl   $4 -- type of temp in extra case
define(`SWITCH_prec', `switch (prec) {
    case blas_prec_single: ifelse(`$4&&$5', `$6&&$7', `{
      SUM_BODY($1, $2)
      break;
    }
    ')dnl
    case blas_prec_double:
    case blas_prec_indigenous: {
      SUM_BODY($1, $3)
      break;
    }

    case blas_prec_extra:
      { SUM_BODY($1, $4, FPU) }
      break;
  }')dnl
dnl
dnl
define(`SUM_X_BODY', `ifelse(
  `$1', `real_S',
  `SWITCH_prec($1, real_S, real_D, real_E)', 

  `$1', `real_D', 
  `SWITCH_prec($1, real_D, real_D, real_E)', 

  `$1', `complex_S', 
  `SWITCH_prec($1, complex_S, complex_D, complex_E)', 

  `$1', `complex_D', 
  `SWITCH_prec($1, complex_D, complex_D, complex_E)'
)')dnl
dnl
dnl
