dnl ----------------------------------------------------------------------
dnl DOT2 --- Dot Product
dnl   r <--- beta * r + alpha * x^T*(head_y + tail_y)
dnl   r <--- beta * r + alpha * x^T*(head_y + tail_y)
dnl ----------------------------------------------------------------------
dnl
dnl
include(cblas.m4)dnl
dnl 
dnl
dnl  Usage:
dnl    DOT2        ($1, $2, $3, $4)
dnl    DOT2_HEAD   ($1, $2, $3, $4)
dnl    DOT2_NAME   ($1, $2, $3, $4)
dnl    DOT2_PARAMS ($1, $2, $3, $4)
dnl    DOT2_COMMENT($1, $2, $3, $4)
dnl
dnl      $1 -- type of alpha, beta, r.
dnl      $2 -- type of a
dnl      $3 -- type of b
dnl      $4 -- Set to `_x' for _x routines.  Otherwise set to `'.
dnl
dnl
define(`DOT2_NAME', `ifelse(
        `$2&&$3', `$1&&$1',`BLAS_$1dot2$4',
        `BLAS_$1dot2_$2_$3$4')')dnl
dnl
dnl
define(`DOT2_PARAMS', 
  `enum blas_conj_type conj, int n, $1_scalar alpha, 
   const $2_array x, int incx, $1_scalar beta,
   const $3_array head_y, const $3_array tail_y, int incy, 
   $1_array r`'ifelse(`$4', `_x', `, enum blas_prec_type prec')')dnl
dnl
dnl
define(`DOT2_HEAD', 
  `void DOT2_NAME($1, $2, $3, $4)(DOT2_PARAMS($1, $2, $3, $4))')dnl
dnl
dnl
define(`DOT2_COMMENT',`
/*
 * Purpose
 * =======
 * 
 * This routine computes the inner product:
 * 
 *   r <- beta * r + alpha * SUM_{i=0, n-1} x[i] * (head_y[i] + tail_y[i]).
 * 
 * Arguments
 * =========
 *  
 * conj   (input) enum blas_conj_type
 *        When x and y are complex vectors, specifies whether vector
 *        components x[i] are used unconjugated or conjugated. 
 * 
 * n      (input) int
 *        The length of vectors x and y.
 * 
 * alpha  (input) $1_scalar
 * 
 * x      (input) const $2_array
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * beta   (input) $1_scalar
 *
 * head_y 
 * tail_y (input) const $3_array
 *        Array of length n.
 *        head_y is the leading part of Y, tail_y is the trailing part.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input/output) $1_array
 * 
PREC_COMMENT($4)dnl
 */')dnl
dnl
dnl
dnl  Usage: DOT2_BODY($1, $2, $3, $4, $5, $6)
dnl  Generates the main body of the product code.
dnl    $1 - type of alpha, beta, r
dnl    $2 - type of x
dnl    $3 - type of (head_y, tail_y)
dnl    $4 - type of sum/prod
dnl    $5 - type of temp
dnl    $6 - [optional] String `FPU' is passed if FPU fix 
dnl         is needed.  Empty string is passed otherwise.
dnl
define(`DOT2_BODY', `
  int i, ix = 0, iy = 0;
  PTR_CAST(r, $1)
  PTR_CAST(x, $2, `const')
  PTR_CAST(head_y, $3, `const')
  PTR_CAST(tail_y, $3, `const')
  SCALAR_CAST(alpha, $1)
  SCALAR_CAST(beta, $1)
  DECLARE(x_ii, $2)
  DECLARE(y_ii, $3)
  DECLARE(r_v, $1)
  DECLARE(prod, $4)
  DECLARE(sum, $4)
  DECLARE(tmp1, $5)
  DECLARE(tmp2, $5)
  ifelse(`$6', `FPU', `FPU_FIX_DECL;')

  /* Test the input parameters. */
  if (n < 0)
    BLAS_error(routine_name, -2, n, NULL);
  else if (incx == 0)
    BLAS_error(routine_name, -5, incx, NULL);
  else if (incy == 0)
    BLAS_error(routine_name, -8, incy, NULL);

  /* Immediate return. */
  if ((TEST_1(beta_i, $1)) && (n == 0 || (TEST_0(alpha_i, $1))))
    return;

  ifelse(`$6', `FPU', `FPU_FIX_START;')

  GET_VECTOR_ELEMENT(r_v, r_i, 0, $1)
  ZERO(sum, $4)
  INC_ADJUST(incx, $2)
  INC_ADJUST(incy, $3)
  if (incx < 0) ix = (-n+1)*incx;
  if (incy < 0) iy = (-n+1)*incy;
dnl *** Only bother to check the value of conj if x is complex.
IF_COMPLEX($2, `
  if (conj == blas_conj) {
    for (i = 0; i < n; ++i) {
      GET_VECTOR_ELEMENT(x_ii, x_i, ix, $2)
      GET_VECTOR_ELEMENT(y_ii, head_y_i, iy, $3)
      CONJ(x_ii, $2, blas_conj)
      MUL(prod, $4, x_ii, $2, y_ii, $3) /* prod = x[i] * head_y[i] */
      ADD(sum, $4, sum, $4, prod, $4) /* sum = sum+prod */
      GET_VECTOR_ELEMENT(y_ii, tail_y_i, iy, $3)
      MUL(prod, $4, x_ii, $2, y_ii, $3) /* prod = x[i] * tail_y[i] */
      ADD(sum, $4, sum, $4, prod, $4) /* sum = sum+prod */
      ix += incx;
      iy += incy;
    } /* endfor */
  } else { 
    /* do not conjugate */
')
  for (i = 0; i < n; ++i) {
    GET_VECTOR_ELEMENT(x_ii, x_i, ix, $2)
    GET_VECTOR_ELEMENT(y_ii, head_y_i, iy, $3)
    CONJ(x_ii, $2, blas_no_conj)
    MUL(prod, $4, x_ii, $2, y_ii, $3) /* prod = x[i] * head_y[i] */
    ADD(sum, $4, sum, $4, prod, $4) /* sum = sum+prod */
    GET_VECTOR_ELEMENT(y_ii, tail_y_i, iy, $3)
    MUL(prod, $4, x_ii, $2, y_ii, $3) /* prod = x[i] * tail_y[i] */
    ADD(sum, $4, sum, $4, prod, $4) /* sum = sum+prod */
    ix += incx;
    iy += incy;
  } /* endfor */
dnl *** Close the outer if, if x is complex
  IF_COMPLEX($2, `}')

  MUL(tmp1, $5, sum, $4, alpha_i, $1) /* tmp1 = sum*alpha */
  MUL(tmp2, $5, r_v, $1, beta_i, $1) /* tmp2 = r*beta */
  ADD(tmp1, $5, tmp1, $5, tmp2, $5) /* tmp1 = tmp1+tmp2 */
  ROUND(r, $1, tmp1, $5) /* r = tmp1 */

  ifelse(`$6', `FPU', `FPU_FIX_STOP;')
')dnl
dnl
dnl
dnl  Usage: SWITCH_prec($1, $2, $3, $4, $5, $6, $7, $8, $9)
dnl  Generates a 3-way switch statement based on prec.
dnl    $1      --  type of alpha, beta, r
dnl    $2      --  type of a
dnl    $3      --  type of b
dnl    $4, $5  --  type of `sum' and `tmp' in single case
dnl    $6, $7  --  type of `sum' and `tmp' in double/indigenous case
dnl    $8, $9  --  type of `sum' and `tmp' in extra case
dnl
define(`SWITCH_prec', `
  switch ( prec ) {
  case blas_prec_single:
    { DOT2_BODY($1, $2, $3, $4, $5) }
    break;
  case blas_prec_double:
  case blas_prec_indigenous:
    { DOT2_BODY($1, $2, $3, $6, $7) }
    break;
  case blas_prec_extra:
    { DOT2_BODY($1, $2, $3, $8, $9, FPU) }
    break;
 }')dnl
dnl
dnl
dnl  Usage: DOT2_X_BODY($1, $2, $3)
dnl  Generates the main body of the extended version of dot code.
dnl    $1 -- type of alpha, beta, r
dnl    $2 -- type of a
dnl    $3 -- type of b
dnl
define(`DOT2_X_BODY', 
  `SWITCH_prec($1, $2, $3, 
    SUM_TYPE_X($1, $2, $3, S), TMP_TYPE_X($1, S), 
    SUM_TYPE_X($1, $2, $3, D), TMP_TYPE_X($1, D), 
    SUM_TYPE_X($1, $2, $3, E), TMP_TYPE_X($1, E))')dnl
dnl
dnl
define(`DOT2', 
  `DOT2_HEAD($1, $2, $3, $4)
   DOT2_COMMENT($1, $2, $3, $4)
   {
     static const char routine_name[] = "DOT2_NAME($1, $2, $3, $4)";
     ifelse($4, _x, `DOT2_X_BODY($1_type, $2_type, $3_type)', 
       `DOT2_BODY($1_type, $2_type, $3_type, 
         SUM_TYPE($1_type, $2_type, $3_type), $1_type)')
   }')dnl
dnl
dnl
define(`PROTOTYPES', `
DOT2_HEAD(s, s, s, _x);
DOT2_HEAD(d, d, d, _x);
DOT2_HEAD(c, c, c, _x);
DOT2_HEAD(z, z, z, _x);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include "blas_extended.h"
#include "blas_extended_private.h"

DOT2(s, s, s, _x)
DOT2(d, d, d, _x)
DOT2(c, c, c, _x)
DOT2(z, z, z, _x)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
