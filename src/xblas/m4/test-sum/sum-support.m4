dnl **********************************************************************
dnl * Computes the ratio of computed error over error bound              *
dnl **********************************************************************
dnl
dnl
include(cblas.m4)dnl
dnl
dnl
define(`TEST_SUM_COMMENT',`
/* Purpose
 * ======= 
 *
 * Computes ratio of the computed error from sum over the expected
 * error bound.
 * 
 * Arguments
 * =========
 *
 * n       (input) int
 *         The length of the vector X.
 *
 * sum_comp    (input) $1_scalar
 *             This result was computed by some other routine, and will be
 *             tested by this routine by comparing it with the truth.
 *
 * sum_true_l  (input) R_TRUE_TYPE($1)
 *             leading part of true sum value
 *
 * sum_true_t  (input) R_TRUE_TYPE($1) 
 *             tailing part of true sum value
 *
 * x       (input) $1_array
 *
 * incx    (input) int
 *
 * eps_int (input) double
 *         The internal machine precision.
 *
 * un_int  (input) double
 *         The internal underflow threshold.
 *
 * test_ratio (output) double*
 *         The ratio of computed error for r over the error bound.
 */')dnl
dnl
dnl
define(`TEST_SUM', 
  `TEST_SUM_HEAD($1)
   TEST_SUM_COMMENT($1)
   {
     TEST_SUM_BODY($1)
   }')dnl
dnl
dnl
define(`TEST_SUM_NAME', `test_BLAS_$1sum')dnl
dnl
dnl
define(`TEST_SUM_PARAMS', 
  `int n, $1_scalar sum_comp, R_TRUE_TYPE($1) sum_true_l, 
   R_TRUE_TYPE($1) sum_true_t, $1_array x, int incx, 
   double eps_int, double un_int, double *test_ratio')dnl
dnl
dnl
define(`TEST_SUM_HEAD', 
  `void TEST_SUM_NAME($1)(TEST_SUM_PARAMS($1))')dnl
dnl
dnl
define(`TEST_SUM_BODY',
`{
        int i;
        int xi;
        int incxi = 1;
        DECLARE(y_elem, $1_type)
        DECLARE(alpha, $1_type)
        DECLARE(beta, $1_type)
        DECLARE(r, $1_type)

        DECLARE_VECTOR(y, $1_type)
        MALLOC_VECTOR(y, $1_type, n)

        INC_ADJUST(incxi, $1_type)

        /* set each element of y to be 1.0 */
        ONE(y_elem, $1_type)
        for (i = 0, xi = 0; i < n; i++, xi += incxi) {
          SET_VECTOR_ELEMENT(y, xi, y_elem, $1_type)
        }

        ONE(alpha, $1_type)
        ZERO(beta, $1_type)
        ZERO(r, $1_type)

        test_BLAS_$1dot(n, blas_no_conj, alpha, beta, r, sum_comp,
                     sum_true_l, sum_true_t, x, incx, y, 1, 
                     eps_int, un_int, test_ratio);
        FREE_VECTOR(y, $1_type)
}')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
TEST_SUM_HEAD(s);
TEST_SUM_HEAD(d);
TEST_SUM_HEAD(c);
TEST_SUM_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

TEST_SUM(s)
TEST_SUM(d)
TEST_SUM(c)
TEST_SUM(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
