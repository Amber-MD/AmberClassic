dnl **********************************************************************
dnl * Computes the ratio of computed error over error bound              *
dnl **********************************************************************
dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: TEST_DOT2_COMMENT(abr_typeltr, x_typeltr, y_typeltr, extended)
dnl
dnl        abr_typeltr : the type and precision of alpha, beta and r
dnl        x_typeltr   : the type and precision of x
dnl        y_typeltr   : the type and precision of y
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ----------------------------------------------------------------------
define(`TEST_DOT2_COMMENT',`
/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from DOT2 over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vectors X and Y.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input) $1_scalar
 *
 * beta    (input) $1_scalar
 *
 * rin     (input) $1_scalar
 *
 * rout    (input) $1_scalar
 *         This result was computed by some other routine, and will be
 *         tested by this routine by comparing it with the truth.
 *
 * r_true_l (input) R_TRUE_TYPE($1)
 *         The leading part of the truth.
 *
 * r_true_t (input) R_TRUE_TYPE($1)
 *         The trailing part of the truth.
 *
 * x       (input) $2_array
 *
 * head_y
 * tail_y  (input) $3_array
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) float*
 *         The ratio of computed error for r over the error bound.
 */')dnl
dnl
dnl
define(`TEST_DOT2_NAME', `test_BLAS_$1dot2`'ifelse(`$2&&$3', `$1&&$1', `', `_$2_$3')')dnl
dnl
dnl
define(`TEST_DOT2_HEAD', 
  `void TEST_DOT2_NAME($1, $2, $3)(int n, enum blas_conj_type conj, dnl
       $1_scalar alpha, $1_scalar beta, $1_scalar rin, $1_scalar rout, dnl
       R_TRUE_TYPE($1) r_true_l, R_TRUE_TYPE($1) r_true_t, $2_array x, dnl
       int incx, $3_array head_y, $3_array tail_y, int incy, dnl
       double eps_int, double un_int, double *test_ratio)')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: TEST_DOT2(abr_typeltr, x_typeltr, y_typeltr, extended)
dnl        ...generate the top level of test_dot2
dnl
dnl        abr_typeltr : the type and precision of alpha, beta and r
dnl        x_typeltr   : the type and precision of x
dnl        y_typeltr   : the type and precision of y
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ----------------------------------------------------------------------
define(`TEST_DOT2', 
 `TEST_DOT2_HEAD($1, $2, $3)
  TEST_DOT2_COMMENT($1, $2, $3)
  IF_REAL($1_type, `TEST_DOT2_BODY_REAL($1_type)', 
    `TEST_DOT2_BODY_COMPLEX($1_type, $2_type, $3_type)')dnl
  /* end of test_BLAS_$1dot$4 */')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: TEST_DOT2_BODY_REAL(abr_typeltr)
dnl        ...generate the body of test_dot2 for real cases
dnl
dnl        abr_type : the type and precision of alpha, beta and r
dnl Each type and precision specifier can be one of
dnl        real_S    ... real and single
dnl        real_D    ... real and double
dnl ----------------------------------------------------------------------
define(`TEST_DOT2_BODY_REAL',`{
    int i, ix, iy;
    double eps_accurate, eps_out, tmp1, S, S1, S2, U;
    double un_d, un_accurate, un_out;

    /* Set the starting position */
    ix = 0; iy = 0;
    if (incx < 0) ix = -(n-1)*incx;
    if (incy < 0) iy = -(n-1)*incy;

    /* computing S */
    S = S1 = S2 = 0.;
    for (i = 0; i < n; ++i) {
      S += fabs(x[ix] * head_y[iy]) + fabs(x[ix] * tail_y[iy]);
      S1 += fabs(x[ix]);
      S2 += fabs(head_y[iy]) + fabs(tail_y[iy]);
      ix += incx; iy += incy;
    }
    S *= fabs(alpha);
    S += fabs(beta * rin);

    SET_UN(real_D, un_d)
    S = MAX(S, un_d);

    SET_EPS(real_E, eps_accurate)
    SET_UN(real_E, un_accurate)
    SET_EPS($1, eps_out)
    SET_UN($1, un_out)
    tmp1 = fabs((rout - r_true_l) - r_true_t);

    /* underflow */
    U = 2*fabs(alpha)*n + 3;
    U = MAX(U, S1 + 2*n + 1);
    U = MAX(U, S2 + 2*n + 1) * (un_int + un_accurate) + un_out;

    *test_ratio = tmp1 /( (n+2)*(eps_int + eps_accurate)*S
                          + eps_out*fabs(r_true_l) + U);
}')dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: TEST_DOT2_BODY_COMPLEX(abr_type, x_type, y_type)
dnl        ...generate the body of test_dot2 for complex cases
dnl
dnl        abr_type : the type and precision of alpha, beta and r
dnl        x_type   : the type and precision of x
dnl        y_type   : the type and precision of y
dnl Each type and precision specifier can be one of
dnl        real_S       ... real and single
dnl        real_D       ... real and double
dnl        complex_S    ... complex and single
dnl        complex_D    ... complex and double
dnl ----------------------------------------------------------------------
define(`TEST_DOT2_BODY_COMPLEX',`{
    int i, ix, iy;
    double eps_accurate, eps_out, tmp1, S, S1, S2, U, prod[2], tmp[2];
    double un_d, un_accurate, un_out;
    PTR_CAST(x, $2)
    PTR_CAST(head_y, $3)
    PTR_CAST(tail_y, $3)
    SCALAR_CAST(alpha, $1)
    SCALAR_CAST(beta, $1)
    SCALAR_CAST(rin, $1)
    SCALAR_CAST(rout, $1)
    DECLARE(x_ii, $2)
    DECLARE(y_ii, $3)

    /* Set the starting position */
    INC_ADJUST(incx, $2)
    INC_ADJUST(incy, $3)
    ix = 0; iy = 0;
    if (incx < 0) ix = -(n-1)*incx;
    if (incy < 0) iy = -(n-1)*incy;

    /* computing S */
    S = S1 = S2 = 0.;
    for (i = 0; i < n; ++i) {
      GET_VECTOR_ELEMENT(x_ii, x_i, ix, $2)
      GET_VECTOR_ELEMENT(y_ii, head_y_i, iy, $3)
dnl *** Only bother to check the value of conj if x is complex.
      IF_REAL($2, 
        `S1 += fabs(x_ii);', 
        `if (conj == blas_conj) CONJ(x_ii, $2, blas_conj)
         S1 += CABS(x_ii);')
      IF_REAL($3, 
        `S2 += fabs(y_ii);', 
        `S2 += CABS(y_ii);')

      MUL(prod, $1, x_ii, $2, y_ii, $3)
      S += CABS(prod);

      /* Now get the tail of y */
      GET_VECTOR_ELEMENT(y_ii, tail_y_i, iy, $3)
      MUL(prod, $1, x_ii, $2, y_ii, $3) /* prod = x[i]*y[i] */
      S += CABS(prod);
        
      IF_REAL($3, 
        `S2 += fabs(y_ii);', 
        `S2 += CABS(y_ii);')
      ix += incx; iy += incy;
    }
    S *= CABS(alpha_i);
    MUL(prod, $1, beta_i, $1, rin_i, $1)
    S += CABS(prod);

    SET_UN(real_D, un_d)
    S = MAX(S, un_d);

    SET_EPS(real_E, eps_accurate)
    SET_UN(real_E, un_accurate)
    SET_EPS($1, eps_out)
    SET_UN($1, un_out)
    tmp[0] = (rout_i[0] - r_true_l[0]) - r_true_t[0];
    tmp[1] = (rout_i[1] - r_true_l[1]) - r_true_t[1];
    tmp1 = CABS(tmp);

    /* underflow */
    U = 2*CABS(alpha_i)*n + 3;
    U = MAX(U, S1 + 2*n + 1);
    U = MAX(U, S2 + 2*n + 1) * (un_int + un_accurate) + un_out;
    U *= 2*sqrt(2.);
dnl U = 2*sqrt(2.)*((CABS(alpha_i)*n+2)*(un_int+un_accurate) + un_out);

    *test_ratio = tmp1 /( 2*sqrt(2.)*(n+2)*(eps_int + eps_accurate)*S
                          + sqrt(2.)*eps_out*CABS(r_true_l) + U);
}
')dnl
dnl
dnl
define(`PROTOTYPES', `FOREACH(`PREC_ARGS', `
TEST_DOT2_HEAD(arg);')')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdlib.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

FOREACH(`PREC_ARGS', `
TEST_DOT2(arg)')')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
