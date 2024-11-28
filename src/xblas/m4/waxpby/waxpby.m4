dnl **********************************************************************
dnl * Perform w = alpha*x + beta*y                                       *
dnl **********************************************************************
dnl
#include "blas_extended.h"
#include "blas_extended_private.h"
include(cblas.m4)dnl
include(waxpby-common.m4)dnl
dnl
dnl
define(`WAXPBY_COMMENT',`
/*
 * Purpose
 * =======
 *
 * This routine computes:
 *
 *     w <- alpha * x + beta * y
 * 
 * Arguments
 * =========
 *
 * n     (input) int
 *       The length of vectors x, y, and w.
 * 
 * alpha (input) $1_scalar
 *
 * x     (input) const $2_array
 *       Array of length n.
 * 
 * incx  (input) int
 *       The stride used to access components x[i].
 *
 * beta  (input) $1_scalar
 *
 * y     (input) $3_array
 *       Array of length n.
 *
 * incy  (input) int
 *       The stride used to access components y[i].
 *
 * w     (output) $1_array
 *       Array of length n.
 *
 * incw  (input) int
 *       The stride used to write components w[i].
 *
PREC_COMMENT($4)dnl
 */')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: WAXPBY(a_b_w_type, x_type, y_type, temp_type)
dnl        ... w = a * x + b * y
dnl        w, x, y: vectors      a, b: scalars
dnl
dnl  $1    a_b_w_type : the type and precision of alpha, beta and w
dnl  $2    x_type   : the type and precision of x
dnl  $3    y_type   : the type and precision of y
dnl  $4    temp_type: the internal prec
dnl
dnl
dnl
dnl
dnl Each type and precision specifier can be one of
dnl        real_S       ... real and single
dnl        real_D       ... real and double
dnl        real_I       ... real and indigenous
dnl        real_E       ... real and extra
dnl        complex_S    ... complex and single
dnl        complex_D    ... complex and double
dnl        complex_I    ... complex and indigeneous
dnl        complex_E    ... complex and extra
dnl ----------------------------------------------------------------------
dnl
define(`WAXPBY_BODY', `
  int i, ix = 0, iy = 0, iw = 0;
  PTR_CAST(w, $1)
  PTR_CAST(x, $2, `const')
  PTR_CAST(y, $3, `const')
  SCALAR_CAST(alpha, $1)
  SCALAR_CAST(beta, $1)
  DECLARE(x_ii, $2)
  DECLARE(y_ii, $3) 
  DECLARE(tmpx, $4)
  DECLARE(tmpy, $4)

  ifelse(`$5', `FPU', `FPU_FIX_DECL;')

  /* Test the input parameters. */
  if ( incx == 0 )
    BLAS_error(routine_name,  -4,  incx, NULL);
  else if ( incy == 0 )
    BLAS_error(routine_name,  -7,  incy, NULL);
  else if ( incw == 0 )
    BLAS_error(routine_name,  -9,  incw, NULL);


  /* Immediate return */
  if ( n <= 0 ) {
    return;
  }

  ifelse(`$5', `FPU', `FPU_FIX_START;')
  
  INC_ADJUST(incx, $2)
  INC_ADJUST(incy, $3)
  INC_ADJUST(incw, $1)
  if ( incx < 0 ) ix = (-n+1)*incx;
  if ( incy < 0 ) iy = (-n+1)*incy;
  if ( incw < 0 ) iw = (-n+1)*incw;

  for (i = 0; i < n; ++i) {         
    GET_VECTOR_ELEMENT(x_ii, x_i, ix, $2)
    GET_VECTOR_ELEMENT(y_ii, y_i, iy, $3)
    MUL(tmpx, $4, alpha_i, $1, x_ii, $2) /* tmpx  = alpha * x[ix] */
    MUL(tmpy, $4, beta_i, $1, y_ii, $3) /* tmpy = beta * y[iy] */
    ADD(tmpy, $4, tmpy, $4, tmpx, $4)
    SET_ROUND_VECTOR_ELEMENT(w_i, iw, tmpy, $4)
    ix += incx;
    iy += incy; 
    iw += incw;
  } /* endfor */

  ifelse(`$5', `FPU', `FPU_FIX_STOP;')
')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SWITCH_prec($1, $2, $3, $4, $5, $6) ... generate
dnl        $4 is the type of 'tmpx' and 'tmpy' in single case.
dnl        $5 is the type of 'tmpx' and 'tmpy' in double/indigenous case.
dnl        $6 is the type of 'tmpx' and 'tmpy' in extra case.
dnl ----------------------------------------------------------------------
dnl 
define(`SWITCH_prec',
 `switch ( prec ) {
    case blas_prec_single: ifelse(`$4', `$5', `', `{
      WAXPBY_BODY($1, $2, $3, $4)
      break;
    }
    ')dnl
    case blas_prec_double:
    case blas_prec_indigenous: {
      WAXPBY_BODY($1, $2, $3, $5)
      break;
    }

    case blas_prec_extra: { 
      WAXPBY_BODY($1, $2, $3, $6, FPU)
      break;
    }
  }')dnl
dnl
dnl
define(`WAXPBY_X_BODY', 
  `SWITCH_prec($1, $2, $3, 
    TMP_TYPE_X($1, S), TMP_TYPE_X($1, D), TMP_TYPE_X($1, E))')dnl
dnl
dnl
define(`WAXPBY', `
  WAXPBY_HEAD($1, $2, $3, $4)
  WAXPBY_COMMENT($1, $2, $3, $4)
  {
    char *routine_name = "WAXPBY_NAME($1, $2, $3, $4)";
    ifelse($4, _x, `WAXPBY_X_BODY($1_type, $2_type, $3_type)', 
      `WAXPBY_BODY($1_type, $2_type, $3_type, $1_type)')
  }')dnl
dnl
dnl
