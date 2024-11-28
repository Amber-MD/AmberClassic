dnl Generates test code for sum
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"
dnl
dnl
#define NORM_START -1
#define NORM_END    1
#define INC_START  -2
#define INC_END     2
#define PREC_START  0
#define PREC_END    2
dnl
dnl
include(cblas.m4)dnl
include(test-common.m4)dnl
include(sum/sum-common.m4)dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_SUM_COMMENT(extended)
dnl        if extended, then print info about prec loop in the code 
dnl        structure
dnl -----------------------------------------------------------------------
define(`DO_TEST_SUM_COMMENT',`
/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on sum  
 *
 * Arguments
 * =========
 *
 * n         (input) int
 *           The size of vector being tested
 *
 * ntests    (input) int
 *           The number of tests to run for each set of attributes.
 *
 * seed      (input/output) int         
 *           The seed for the random number generator used in testgen().
 *
 * thresh    (input) double
 *           When the ratio returned from test() exceeds the specified
 *           threshold, the current size, r_true, r_comp, and ratio will be
 *           printed.  (Since ratio is supposed to be O(1), we can set thresh
 *           to ~10.)
 *
 * debug     (input) int
 *           If debug=3, print summary 
 *           If debug=2, print summary only if the number of bad ratios > 0
 *           If debug=1, print complete info if tests fail
 *           If debug=0, return max ratio
 *
 * min_ratio (output) double
 *           The minimum ratio
 * 
 * num_bad_ratio (output) int
 *               The number of tests fail; they are above the threshold.
 *
 * num_tests (output) int
 *           The number of tests is being performed.
 *
 */')dnl
dnl
dnl
define(`DO_TEST_SUM_NAME', `do_test_$1sum$2')dnl
dnl
dnl
define(`DO_TEST_SUM_PARAMS', 
  `int n, int ntests, int *seed, double thresh, 
   int debug, float test_prob, double *min_ratio, double *max_ratio, 
   int *num_bad_ratio, int *num_tests')dnl
dnl
dnl
define(`DO_TEST_SUM_HEAD', 
  `void DO_TEST_SUM_NAME($1, $2)(DO_TEST_SUM_PARAMS($1, $2))')dnl
dnl
dnl
define(`DO_TEST_SUM', 
  `DO_TEST_SUM_HEAD($1, $2)
   DO_TEST_SUM_COMMENT($1, $2)
   {
     DO_TEST_SUM_BODY($1, $2)
   }')dnl
dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: DO_TEST_SUM_BODY(sum_typeltr, x_typeltr, extended)
dnl
dnl        sum_typeltr : the type and precision of sum
dnl        x_typeltr   : the type and precision of x
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ---------------------------------------------------------------------
define(`DO_TEST_SUM_BODY',
`
  /* function name */
  const char fname[] = "BLAS_$1sum$2";

  int i, j;
  int norm;
  int xi, incx_val, incx;
  double ratio_max, ratio_min;
  double ratio;
  int saved_seed;
  double eps_int;
  double un_int;

  int test_count;
  int bad_ratio_count;

  DECLARE_VECTOR(x, $1_type)
  DECLARE_VECTOR(x_gen, $1_type)
  DECLARE(x_elem, $1_type)

  DECLARE(sum_true, EXTRA_TYPE($1_type))

  DECLARE(sum, $1_type) 

  ifelse(`$2', `_x', `int prec_val;')
  enum blas_prec_type prec;
  int x_gen_i, incx_gen;

  FPU_FIX_DECL;


  /* test for bad arguments */
  if (n < 0 || ntests < 0) 
    BLAS_error(fname, 0, 0, NULL);   

  /* if there is nothing to test, return all zero */
  if (n == 0 || ntests == 0){
    *min_ratio = 0.0;
    *max_ratio = 0.0;
    *num_bad_ratio = 0;
    *num_tests = 0;
    return;
  }

  FPU_FIX_START;

  incx_gen = 1;
  INC_ADJUST(incx_gen, $1_type)

  /* get space for calculation */
  MALLOC_VECTOR(x, $1_type, n*2)
  MALLOC_VECTOR(x_gen, $1_type, n)

  /* initialization */
  ratio_min = 1e308;
  ratio_max = 0.0;
  test_count = 0;
  bad_ratio_count = 0;

  ifelse($2, _x, `
  /* varying extra precs */
  for (prec_val = PREC_START; prec_val <= PREC_END; prec_val++) {')
    SET_INTERNAL_PARAMS($1_type, $2)

       for (norm = NORM_START; norm <= NORM_END; norm++) {
           
         for (incx_val = INC_START; incx_val <= INC_END; incx_val++) {

            if (incx_val == 0)
              continue;
 
            for (i = 0; i < ntests; i++){               

                /* For the sake of speed, we throw out this case at random */
                if ( xrand(seed) >= test_prob ) continue;

                saved_seed = *seed;

                BLAS_$1sum_testgen(n, norm, x_gen, seed, 
                                PASS_BY_REF(sum_true, EXTRA_TYPE($1_type)));

                test_count++;

                incx = incx_val;
                INC_ADJUST(incx, $1_type)               
 
                xi = 0;
                if (incx < 0)
                  xi = - (n-1) * incx;

                /* copy x_gen to x */
                incx_gen = 1;
                INC_ADJUST(incx_gen, $1_type)
                for (j = 0, x_gen_i = 0; j < n; 
                     j++, x_gen_i += incx_gen, xi += incx) {
                  GET_VECTOR_ELEMENT(x_elem, x_gen, x_gen_i, $1_type)
                  SET_VECTOR_ELEMENT(x, xi, x_elem, $1_type)
                }

                FPU_FIX_STOP;
                BLAS_$1sum$2(n, x, incx_val, PASS_BY_REF(sum, $1_type)
                          ifelse(`$2', `_x', `, prec'));
                FPU_FIX_START;

                test_BLAS_$1sum(n, sum, HEAD(sum_true), TAIL(sum_true), x, incx_val, 
                             eps_int, un_int, &ratio);

               /* The !<= below causes NaN error to be detected.
                  Note that (NaN > thresh) is always false */
                if ( !(ratio <= thresh) ) {
                  bad_ratio_count++;

                  if (debug == 3) {
                    printf("Seed = %d\n", saved_seed);
                    printf("n = %d\n", n);
                    printf("norm = %d\n", norm);
                    ifelse(`$2', `_x', `PRINT_PREC(prec)')
 
                    printf("incx = %d\n", incx);
                    
                    /* print out the vector */
                    xi = (incx < 0) ? -(n-1) * incx : 0;
                    printf(" [ ");
                    for (j = 0; j < n; j++, xi += incx) {
                      PRINT_ARRAY_ELEM(x, xi, $1_type)
                    }
                    printf("]\n");

                    PRINT_VAR(sum, $1_type) 
                    printf("\n");
                    printf("ratio = %.4e\n", ratio);
                    PRINT_VAR(sum_true, EXTRA_TYPE($1_type))

                  }  /* end of if (debug == 3) */
                } /* end of if (ratio > thresh) */

                if (ratio > ratio_max)
                  ratio_max = ratio;
                if (ratio != 0.0 && ratio < ratio_min)
                  ratio_min = ratio;

              } /* end of ntests loop */
            } /* end of incx loop */
          } /* end of norm loop */

ifelse(`$2', `_x', `} /* end of prec loop */')

   FPU_FIX_STOP;

   *max_ratio = ratio_max;
   *min_ratio = ratio_min;
   *num_tests = test_count;
   *num_bad_ratio = bad_ratio_count;

   FREE_VECTOR(x, $1_type)
   FREE_VECTOR(x_gen, $1_type)

' )dnl
dnl
dnl
define(`CALL_DO_TEST_SUM', 
  `fname = "SUM_NAME($1, $2)";
   printf("Testing %s...\n", fname);
   total_tests = 0;
   total_bad_ratios = 0;
   total_min_ratio = 1e308;
   total_max_ratio = 0.0;

   for (n = 0; n <= nsizes; n++) {
     DO_TEST_SUM_NAME($1, $2)(n, ntests, &seed, thresh, debug, 
           test_prob, &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);
     if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
       printf("  n = %d:  ", n);
       printf("bad/total = %d/%d, min_ratio = %g, max_ratio = %g\n", 
              num_bad_ratio, num_tests, min_ratio, max_ratio);
     }
     total_tests += num_tests;
     total_bad_ratios += num_bad_ratio;
     if (total_min_ratio > min_ratio)
       total_min_ratio = min_ratio;
     if (total_max_ratio < max_ratio)
       total_max_ratio = max_ratio;
   }

   nr_routines++;
   if (total_bad_ratios == 0)
     printf("PASS> ");
   else {
     printf("FAIL> ");
     nr_failed_routines++;
   }
  printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n",
      fname, total_bad_ratios, total_tests, max_ratio);
   printf("\n");
  ')

FOREACH(`SUM_ARGS', `
DO_TEST_SUM(arg)')dnl

MAIN(`', `

FOREACH(`SUM_ARGS', `
CALL_DO_TEST_SUM(arg)')')dnl

