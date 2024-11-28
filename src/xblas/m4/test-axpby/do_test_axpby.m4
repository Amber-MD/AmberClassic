dnl Generates test code for axpby
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"
include(cblas.m4)dnl
include(test-common.m4)dnl
include(axpby/axpby-common.m4)dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_AXPBY_COMMENT(extended)
dnl        if extended, then print info about prec loop in the code 
dnl        structure
dnl -----------------------------------------------------------------------
define(`DO_TEST_AXPBY_COMMENT',`
/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on axpby  
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
 *           threshold, the current size, y_true, y_comp, and ratio will be
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
 * Return value
 * ============
 *
 * The maximum ratio if run successfully, otherwise return -1 
 *
 * Code structure
 * ==============
 * 
 *  debug loop  -- if debug is one, the first loop computes the max ratio
 *              -- and the last(second) loop outputs debugging information,
 *              -- if the test fail and its ratio > 0.5 * max ratio.
 *              -- if debug is zero, the loop is executed once
 *    alpha loop  -- varying alpha: 0, 1, or random
 *      beta loop   -- varying beta: 0, 1, or random
ifelse(`$1', `_x',` *        prec loop   -- varying internal prec: single, double, or extra', `')
 *          norm loop   -- varying norm: near undeflow, near one, or 
 *                        -- near overflow
 *            numtest loop  -- how many times the test is perform with 
 *                            -- above set of attributes
 *                incx loop     -- varying incx: -2, -1, 1, 2
 *                  incy loop     -- varying incy: -2, -1, 1, 2
 */')dnl
dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_AXPBY(aby_typeltr, x_typeltr, extended)
dnl
dnl        aby_typeltr : the type and precision of alpha, beta and y
dnl        x_typeltr   : the type and precision of x
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ----------------------------------------------------------------------
define(`DO_TEST_AXPBY_NAME', `do_test_$1axpby`'ifelse(`$1', `$2', `$3', `_$2$3')')dnl
dnl
dnl
define(`DO_TEST_AXPBY', 
`double DO_TEST_AXPBY_NAME($1, $2, $3)(int n, int ntests, int *seed, double thresh, dnl
    int debug, double *min_ratio, int *num_bad_ratio, int *num_tests)
DO_TEST_AXPBY_COMMENT($3)
DO_TEST_AXPBY_BODY($1, $2, $3) /* end of DO_TEST_AXPBY_NAME($1, $2, $3) */')dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: DO_TEST_AXPBY_BODY(aby_typeltr, x_typeltr, extended)
dnl
dnl        aby_typeltr : the type and precision of alpha, beta and y
dnl        x_typeltr   : the type and precision of x
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ---------------------------------------------------------------------
define(`DO_TEST_AXPBY_BODY',
`{                                        
  /* function name */
  const char fname[] = "AXPBY_NAME($1, $2, $3)";

  /* max number of debug lines to print */
  const int max_print = 32;

  /* Variables in the "x_val" form are loop vars for corresponding
     variables */
  int i;            /* iterate through the repeating tests */
  int ix, iy;       /* use to index x and y respectively */
  int incx_val, incy_val, /* for testing different inc values */
      incx, incy, incx_gen, incy_gen,  
      ygen_val, xgen_val, test_val;   
  int d_count;      /* counter for debug */
  int find_max_ratio; /* find_max_ratio = 1 only if debug = 3 */
  int p_count;      /* counter for the number of debug lines printed*/
  int tot_tests;    /* total number of tests to be done */
  int norm;         /* input values of near underflow/one/overflow */
  double ratio_max; /* the current maximum ratio */
  double ratio_min; /* the current minimum ratio */
  double ratio;     /* the per-use test ratio from test() */
  double new_ratio;
  int bad_ratios;   /* the number of ratios over the threshold */
  double eps_int;   /* the internal epsilon expected--2^(-24) for float */
  double un_int;    /* the internal underflow threshold */
  DECLARE(alpha, $1_type)
  DECLARE(beta, $1_type)
  DECLARE_VECTOR(x, $2_type)
  DECLARE(y_fix, $2_type)
  DECLARE_VECTOR(y_ori, $1_type)
  DECLARE_VECTOR(y_comp, $1_type)  /* the y computed  by AXPBY_NAME($1, $2, $3) */

  int ixmax, iymax;             
  int ixmax_max, iymax_max;      
                
  /* x_gen and y_gen are used to store vectors generated by testgen.
     they eventually are copied back to x and y */
  DECLARE_VECTOR(x_gen, $2_type)
  DECLARE_VECTOR(y_gen, $1_type) 
 
  /* the true y calculated by testgen(), in double-double */
  DECLARE_VECTOR(y_true, EXTRA_TYPE($1_type))
  int alpha_val;
  int alpha_flag;   /* input flag for DOT_TESTGEN_NAME($1, $2, $2) */
  int beta_val;
  int beta_flag;    /* input flag for DOT_TESTGEN_NAME($1, $2, $2) */
  ifelse(`$3', `_x', `int prec_val;')
  enum blas_prec_type prec;
  int saved_seed;   /* for saving the original seed */
  int count, old_count;  /* use for counting the number of testgen calls * 2 */ 

  FPU_FIX_DECL;

  /* test for bad arguments */
  if (n < 0 || ntests < 0) 
    BLAS_error(fname, 0, 0, NULL);   

  /* if there is nothing to test, return all zero */
  if (n == 0 || ntests == 0){
    *min_ratio = 0.0;
    *num_bad_ratio = 0;
    *num_tests = 0;
    return 0.0;
  }

  FPU_FIX_START;

  /* get space for calculation */
  MALLOC_VECTOR(x, $2_type, n*2)
  MALLOC_VECTOR(y_ori, $1_type, n*2)
  MALLOC_VECTOR(y_comp, $1_type, n*2)
  MALLOC_VECTOR(y_true, EXTRA_TYPE($1_type), n)
  MALLOC_VECTOR(x_gen, $2_type, n)
  MALLOC_VECTOR(y_gen, $1_type, n)

  /* initialization */
  saved_seed = *seed;
  ratio_min = 1e308;
  ratio_max = 0.0;
  tot_tests = 0;
  p_count = 0;
  count = 0;
  bad_ratios = 0;
  ixmax_max = 0;
  iymax_max = 0;
  alpha_flag = 0;
  beta_flag = 0;
  old_count = 0;
  
  find_max_ratio = 0;
  if (debug == 3)
    find_max_ratio = 1;
  ONE(y_fix, $2_type)

  /* initialize incx_gen and incy_gen */
  incx_gen = 1;
  incy_gen = 1;
  INC_ADJUST(incx_gen, $2_type)
  INC_ADJUST(incy_gen, $1_type)

  /* The debug iteration:
     If debug=1, then will execute the iteration twice. First, compute the
     max ratio. Second, print info if ratio > (50% * ratio_max). */
  for (d_count=0; d_count<= find_max_ratio; d_count++) {
    bad_ratios = 0; /* set to zero */ 
  
    if ((debug == 3) && (d_count == find_max_ratio))
      *seed = saved_seed; /* restore the original seed */

    /* varying alpha */
    for (alpha_val=0; alpha_val<3; alpha_val++) { 
      SET_ALPHA($1_type)
        
      /* varying beta */
      for (beta_val=0; beta_val<3; beta_val++) {
        SET_BETA($1_type)

        ifelse($3, _x, `
        /* varying extra precs */
        for (prec_val = 0; prec_val <= 2; prec_val++) {')
          SET_INTERNAL_PARAMS($1_type, $3)
           
          /* values near underflow, 1, or overflow */
          for (norm = -1; norm <= 1; norm++) {

            /* number of tests */
            for (i=0; i<ntests; i++){
        
              /* generate test inputs */                
              DOT_TESTGEN_NAME($1, $2, $2)(1, 0, 1, norm, blas_no_conj, dnl
                PASS_BY_REF(alpha, $1_type), alpha_flag, dnl
                PASS_BY_REF(beta, $1_type), beta_flag, PASS_BY_REF(y_fix, $2_type), dnl
                x_gen, seed, y_gen, HEAD(y_true), TAIL(y_true));
              xgen_val = incx_gen;
              for ( ygen_val = incy_gen; ygen_val < n*incy_gen; ygen_val+=incy_gen ) {
                DOT_TESTGEN_NAME($1, $2, $2)(1, 0, 1, norm, blas_no_conj, dnl
                  PASS_BY_REF(alpha, $1_type), 1, PASS_BY_REF(beta, $1_type), 1, dnl
                  PASS_BY_REF(y_fix, $2_type), &x_gen[xgen_val], seed, &y_gen[ygen_val], dnl
                  &HEAD(y_true)[ygen_val], &TAIL(y_true)[ygen_val]);
                xgen_val += incx_gen; 
              }
         
              count++;
              
              /* varying incx */
              for (incx_val = -2; incx_val <= 2; incx_val++){
                if (incx_val == 0) continue;           

                /* setting incx */
                incx = incx_val;
                INC_ADJUST(incx, $2_type)

                $2copy_vector(x_gen, n, 1, x, incx_val);

                /* varying incy */
                for (incy_val = -2; incy_val <= 2; incy_val++){
                  if (incy_val == 0) continue;

                    /* setting incy */
                    incy = incy_val;
                    INC_ADJUST(incy, $1_type)           

                    $1copy_vector(y_gen, n, 1, y_ori, incy_val);
                    $1copy_vector(y_gen, n, 1, y_comp, incy_val);

                    /* call AXPBY_NAME($1, $2, $3) to get y_comp */
                    FPU_FIX_STOP;
                    AXPBY_NAME($1, $2, $3)(n, alpha, x, incx_val, beta,
                                           y_comp, incy_val ifelse(`$3', `_x', `, prec'));
                    FPU_FIX_START;

                        
                    /* computing the ratio */
                    ix = 0;
                    if (incx < 0) ix = -(n-1)*incx;
                    iy = 0;
                    if (incy < 0) iy = -(n-1)*incy;
                    ratio = 0.0;

                    ixmax = -1;
                    iymax = -1;
                    for ( test_val = 0; test_val < n*incy_gen; test_val+=incy_gen ) {
                      TEST_DOT_NAME($1, $2, $2, $3)(1, blas_no_conj, alpha, beta, dnl
                        VECTOR_ELEMENT(y_ori, iy, $1_type), VECTOR_ELEMENT(y_comp, iy, $1_type), dnl
                        VECTOR_ELEMENT(y_true, test_val, EXTRA_TYPE($1_type)), dnl
                        PASS_BY_REF(y_fix, $2_type), incx, &x[ix], incx, eps_int, dnl
                        un_int, &new_ratio);
                                                
                      ix += incx;
                      iy += incy;
                      if (MAX(ratio, new_ratio) == new_ratio) {
                        ixmax = ix;
                        iymax = iy;
                      }
                      ratio = MAX(ratio, new_ratio);
                    }

                    /* increase the number of bad ratio, if the ratio
                       is bigger than the threshold.
                       The !<= below causes NaN error to be detected.
                       Note that (NaN > thresh) is always false. */
                    if ( !(ratio <= thresh) ) {
                      bad_ratios++;
                        
                        
                      if ((debug == 3) &&         /* print only when debug is on */
                         (count != old_count) &&  /* print if old vector is different */
                                                  /* from the current one */ 
                         (d_count == find_max_ratio) &&  
                         (p_count <= max_print) &&
                         (ratio > 0.5*ratio_max))
                      {
                        old_count = count;      
                     
                        printf("FAIL> %s: n = %d, ntests = %d, threshold = %4.2f,\n", 
                                fname, n, ntests, thresh);
                        printf("seed = %d\n", *seed);
                        printf("norm = %d\n", norm);


                        /* Print test info */
                        PRINT_PREC(prec)
                        PRINT_NORM(norm)

                        printf("incx=%d, incy=%d:\n", incx, incy);
                        
                        $2print_vector(x, n, incx_val, "x");
                        $1print_vector(y_ori, n, incy_val, "y_ori");
                        $1print_vector(y_comp, n, incy_val, "y_comp");

                        printf("      "); PRINT_VAR(alpha, $1_type) printf("; "); 
                        PRINT_VAR(beta, $1_type) printf("\n");
                        printf("      ratio=%.4e\n", ratio);
                        printf("iymax = %d\n", iymax);
                        printf("ixmax = %d\n", ixmax);
                        p_count++;     
                      }
                    }
                    if (d_count == 0) {
                      if (ratio > ratio_max)
                        ratio_max = ratio;
                        iymax_max = iymax;
                        ixmax_max = ixmax;
                      if (ratio != 0.0 && ratio < ratio_min)
                        ratio_min = ratio;

                      tot_tests++;
                    }
                  } /* incy */
                } /* incx */
            } /* tests */
          } /* norm */
ifelse(`$3', `_x', `         } /* prec */')
      } /* beta */
    } /* alpha */
  } /* debug */

  if ((debug == 2) ||
     ((debug == 1) && (bad_ratios > 0))){
    printf("      %s:  n = %d, ntests = %d, thresh = %4.2f\n", 
            fname, n, ntests, thresh);
    if ( ratio_min == 1.0e+308 )
        ratio_min = 0.0;
    printf("      bad/total = %d/%d=%3.2f, min_ratio = %.4e, max_ratio = %.4e\n\n",
            bad_ratios, tot_tests,
           ((double)bad_ratios)/((double)tot_tests), ratio_min, ratio_max);
    printf("iymax_max = %d, ixmax_max = %d\n", iymax_max, ixmax_max);
  }

  FREE_VECTOR(x, $2_type)
  FREE_VECTOR(y_ori, $1_type)
  FREE_VECTOR(y_comp, $1_type)
  FREE_VECTOR(y_true, EXTRA_TYPE($1_type))
  FREE_VECTOR(x_gen, $2_type)
  FREE_VECTOR(y_gen, $1_type)

  *min_ratio = ratio_min;
  *num_bad_ratio = bad_ratios;
  *num_tests = tot_tests;

  FPU_FIX_STOP;
  return ratio_max;
}' )dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: DO_TEST_AXPBY_NAME(aby_typeltr, x_typeltr, extended)
dnl        create do_test_axpby name
dnl -------------------------------------------------------------------
define(`DO_TEST_AXPBY_NAME', `ifelse(
        `$1', `$2',
        `do_test_$1axpby$3',
        `do_test_$1axpby_$2$3')') dnl
dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: CALL_DO_TEST_AXPBY(aby_typeltr, x_typeltr, extended)
dnl
dnl        aby_type : the type and precision of alpha, beta and y
dnl        x_type   : the type and precision of x
dnl        extended : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl -------------------------------------------------------------------
define(`CALL_DO_TEST_AXPBY',           
  `min_ratio = 1e308; max_ratio = 0.0;
   total_bad_ratios = 0; total_tests = 0;
   fname = "AXPBY_NAME($1, $2, $3)";
   printf("Testing %s...\n", fname);
   for(n=0; n<=nsizes; n++){
     total_max_ratio = DO_TEST_AXPBY_NAME($1, $2, $3)(n, ntests, dnl
         &seed, thresh, debug, &total_min_ratio, &num_bad_ratio, dnl
         &num_tests);
     if (total_max_ratio > max_ratio)
       max_ratio = total_max_ratio;

     if (total_min_ratio < min_ratio)
       min_ratio = total_min_ratio;

     total_bad_ratios += num_bad_ratio;
     total_tests += num_tests;
   }

   if (total_bad_ratios == 0)
     printf("PASS> ");
   else{
     printf("FAIL> ");
     nr_failed_routines++;
   }
   nr_routines++;

   printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n", dnl
          fname, total_bad_ratios, total_tests, max_ratio);')
dnl
dnl
dnl
FOREACH(`AXPBY_ARGS', `
DO_TEST_AXPBY(arg)
')dnl

MAIN(`', `FOREACH(`AXPBY_ARGS', `
CALL_DO_TEST_AXPBY(arg)
')')dnl

