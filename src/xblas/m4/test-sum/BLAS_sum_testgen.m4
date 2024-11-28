dnl **********************************************************************
dnl Generate routines to test DOT.
dnl **********************************************************************
dnl
include(cblas.m4)dnl
dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SUM_TESTGEN(sum_type, x_type, int_typeltr)
dnl        sum_type    : the type and precision of sum
dnl        x_type      : the type and precision of x
dnl        int_typeltr : the internal precision used to generate test vectors;
dnl                      it can be one of s, d, c, or z
dnl ----------------------------------------------------------------------
define(`SUM_TESTGEN_BODY',
`{  
        int i;
        int xi, incxi = 1;
        DECLARE(alpha, $1_type)
        DECLARE(beta, $1_type)
        DECLARE(r, $1_type)
        DECLARE(x_elem, $1_type)
        DECLARE_VECTOR(tmp, $1_type)

        PTR_CAST(x, $1_type)

        ONE(alpha, $1_type)
        ZERO(beta, $1_type)
        
        MALLOC_VECTOR(tmp, $1_type, n)

        INC_ADJUST(incxi, $1_type)

        ONE(x_elem, $1_type)
        for (i = 0, xi = 0; i < n; i++, xi += incxi) {
            SET_VECTOR_ELEMENT(tmp, xi, x_elem, $1_type)
        }

        /* Call generator now. */
        testgen_BLAS_$1dot(n, 0, n, norm, blas_conj,    
                        PASS_BY_REF(alpha, $1_type), 1, PASS_BY_REF(beta, $1_type), 1,     
                        tmp, x_i, seed, &r, sum_true_l,
                        sum_true_t);

        FREE_VECTOR(tmp, $1_type)
}')dnl
dnl
dnl ----------------------------------------------------------------------
dnl Usage: SUM_TESTGEN_COMMENT(sum_typeltr, x_typeltr)
dnl        ... generate the leading comment for the summation
dnl            testing routines.
dnl Each _typeltr specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ----------------------------------------------------------------------
dnl
define(`SUM_TESTGEN_COMMENT',`
/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_$1sum{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vector X.
 *
 * norm    (input) int
 *         = -1 : the vector is scaled with norms near underflow.
 *         = 0  : the vector has norms of order 1.
 *         = 1  : the vector is scaled with norms near overflow.
 *
 * x       (output) $1_array contains generated test values
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 *
 * sum_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * sum_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */')dnl
dnl
dnl
define(`SUM_TESTGEN_NAME', `BLAS_$1sum_testgen')dnl
dnl
dnl
define(`SUM_TESTGEN_PARAMS', 
  `int n, int norm, $1_array x, int *seed, 
   double *sum_true_l, double *sum_true_t')dnl
dnl
dnl
define(`SUM_TESTGEN_HEAD', 
  `void SUM_TESTGEN_NAME($1)(SUM_TESTGEN_PARAMS($1))')dnl
dnl
dnl
define(`SUM_TESTGEN', 
  `SUM_TESTGEN_HEAD($1)
   SUM_TESTGEN_COMMENT($1)
   {
     SUM_TESTGEN_BODY($1)
   }')dnl
dnl
dnl
define(`PROTOTYPES', `dnl
SUM_TESTGEN_HEAD(s);
SUM_TESTGEN_HEAD(d);
SUM_TESTGEN_HEAD(c);
SUM_TESTGEN_HEAD(z);
')dnl
dnl
dnl
define(`SOURCE', `dnl
#include <stdlib.h>
#include "blas_extended.h"
#include "blas_extended_test.h"

SUM_TESTGEN(s)
SUM_TESTGEN(d)
SUM_TESTGEN(c)
SUM_TESTGEN(z)
')dnl
dnl
dnl
ifdef(`prototypes_only', `PROTOTYPES()', `SOURCE()')dnl
dnl
dnl
