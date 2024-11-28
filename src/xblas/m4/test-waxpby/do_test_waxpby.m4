dnl Generates test code for waxpby
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"
include(cblas.m4)dnl
include(test-common.m4)dnl
include(waxpby/waxpby-common.m4)dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_WAXPBY_COMMENT(extended)
dnl        if extended, then print info about prec loop in the code 
dnl        structure
dnl -----------------------------------------------------------------------
define(`DO_TEST_WAXPBY_COMMENT',`
/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on waxpby  
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
 *           threshold, the current size, w_true, w, and ratio will be
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
dnl Usage: DO_TEST_WAXPBY(abw_typeltr, x_typeltr, y_typeltr, extended)
dnl
dnl        abw_typeltr : the type and precision of alpha, beta and w
dnl        x_typeltr   : the type and precision of x
dnl        y_typeltr   : the type and precision of y
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ----------------------------------------------------------------------
define(`DO_TEST_WAXPBY', `ifelse(
        `$1&&$1', `$2&&$3',`
double do_test_$1waxpby$4(int n,
                     int ntests,
                     int *seed,
                     double thresh, 
                     int debug, float test_prob, 
                     double *min_ratio,
                     int *num_bad_ratio,
                     int *num_tests)
DO_TEST_WAXPBY_COMMENT($4)
DO_TEST_WAXPBY_BODY($1, $2, $3, $4) /* end of do_test_$1waxpby$4 */',`
double do_test_$1waxpby_$2_$3$4 (int n, 
                          int ntests,
                          int *seed,
                          double thresh, 
                          int debug, float test_prob, 
                          double *min_ratio,
                          int *num_bad_ratio,
                          int *num_tests)
DO_TEST_WAXPBY_COMMENT($4)
DO_TEST_WAXPBY_BODY($1, $2, $3, $4) /* end of do_test_$1waxpby_$2_$3$4 */')') dnl
dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: DO_TEST_WAXPBY_BODY(abw_typeltr, x_typeltr, y_typeltr, extended)
dnl
dnl        abw_typeltr : the type and precision of alpha, beta and w
dnl        x_typeltr   : the type and precision of x
dnl        y_typeltr   : the type and precision of y
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ---------------------------------------------------------------------
define(`DO_TEST_WAXPBY_BODY',
`{                                        
  /* function name */
  const char fname[] = "WAXPBY_NAME($1, $2, $3, $4)";

  /* max number of debug lines to print */
  const int max_print = 32;

  /* Variables in the "x_val" form are loop vars for corresponding
     variables */
  int i;            /* iterate through the repeating tests */
  int j;            /* multipurpose counter */
  int ix, iy, iw;       /* use to index x, y, w respectively */
  int incx_val, incy_val, incw_val, /* for testing different inc values */
      incx, incy, incw, gen_val, test_val;   
  int d_count;      /* counter for debug */
  int find_max_ratio; /* find_max_ratio = 1 only if debug = 3 */
  int p_count;      /* counter for the number of debug lines printed*/
  int tot_tests;    /* total number of tests to be done */
  int norm;         /* input values of near underflow/one/overflow */
  int X_int;
  double X;
  double ratio_max; /* the current maximum ratio */
  double ratio_min; /* the current minimum ratio */
  double ratio;     /* the per-use test ratio from test() */
  double new_ratio;
  int bad_ratios;   /* the number of ratios over the threshold */
  double eps_int;   /* the internal epsilon expected--2^(-24) for float */
  double un_int;    /* the internal underflow threshold */
  DECLARE(x_i, $2_type)
  DECLARE(y_i, $3_type)
  DECLARE(alpha, $1_type)
  DECLARE(beta, $1_type)
  DECLARE_VECTOR(x, $2_type)
  DECLARE_VECTOR(y, $3_type)
  DECLARE_VECTOR(w, $1_type)  /* the w computed  by WAXPBY_NAME($1, $2, $3, $4) */
  DECLARE(x_fix1, $2_type)
  DECLARE(x_fix2, $3_type)
  DECLARE(zero, $1_type)
  DECLARE(one, $1_type)
  DECLARE(dummy, $1_type)

  /* x_gen and y_gen are used to store vectors generated by testgen.
     they eventually are copied back to x and y */
  DECLARE_VECTOR(x_gen, $2_type)
  DECLARE_VECTOR(y_gen, $3_type) 
  DECLARE_VECTOR(temp_ab, $1_type)
  DECLARE_VECTOR(temp_xy, $2_type)  


  /* added by DY */
  DECLARE(x_genj, $2_type)
  DECLARE(y_genj, $3_type)
  int incy_gen, incx_gen, incw_gen;
  int xgen_val, ygen_val, wgen_val;  
  int iymax, ixmax;
  float xtemp;
  float ytemp;
  float atemp;
  float btemp;
  double wltemp;
  double wttemp;
  float x_fix1_temp;

  /* the true w calculated by testgen(), in double-double */
  DECLARE_VECTOR(w_true, EXTRA_TYPE($1_type))
  ifelse(`$4', `_x', `int prec_val;')
  enum blas_prec_type prec;
  int saved_seed;   /* for saving the original seed */
  int count, old_count;  /* use for counting the number of testgen calls * 2 */ 

  FPU_FIX_DECL;

  /* There are there to get rid of compiler warnings.
     Should modify M4 code to not even produce these variables when not
     needed. */
  xtemp = ytemp = atemp = btemp = 0.0;
  wltemp = wttemp = x_fix1_temp = 0.0;
  ZERO(x_i, $2_type)
  ZERO(y_i, $3_type)
  X = 0.0;
  X_int = 0;  
  gen_val = 0;

  /* test for bad arguments */
  if (n < 0 )
    BLAS_error(fname,  -1,  n, NULL);
  if (ntests < 0) 
    BLAS_error(fname,  -2,  ntests, NULL);

  /* if there is nothing to test, return all zero */
  if (n == 0 || ntests == 0){
    *min_ratio = 0.0;
    *num_bad_ratio = 0;
    *num_tests = 0;
    return 0.0;
  }

  FPU_FIX_START;

  incw_gen = 1;
  incx_gen = 1;
  incy_gen = 1;
  INC_ADJUST(incw_gen, $1_type)
  INC_ADJUST(incx_gen, $2_type)
  INC_ADJUST(incy_gen, $3_type)
 
  /* get space for calculation */
  MALLOC_VECTOR(x, $2_type, n*2)
  MALLOC_VECTOR(y, $3_type, n*2)
  MALLOC_VECTOR(w, $1_type, n*2)
  MALLOC_VECTOR(w_true, EXTRA_TYPE($1_type), n)
  MALLOC_VECTOR(x_gen, $2_type, n)
  MALLOC_VECTOR(y_gen, $3_type, n)
  MALLOC_VECTOR(temp_ab, $1_type, 2)
  MALLOC_VECTOR(temp_xy, $2_type, 2)

  /* initialization */
  saved_seed = *seed;
  ratio_min = 1e308;
  ratio_max = 0.0;
  tot_tests = 0;
  p_count = 0;
  count = 0;
  old_count = 0;
  bad_ratios = 0;

  find_max_ratio = 0;
  if (debug == 3)
    find_max_ratio = 1;
  ONE(x_fix1, $2_type)
  ONE(x_fix2, $3_type)
  ZERO(zero, $1_type)
  ONE(one, $1_type)
  ZERO(dummy, $1_type);
  

  /* The debug iteration:
     If debug=1, then will execute the iteration twice. First, compute the
     max ratio. Second, print info if ratio > (50% * ratio_max). */
  for (d_count=0; d_count<= find_max_ratio; d_count++) {
    bad_ratios = 0; /* set to zero */ 

    if ((debug == 3) && (d_count == find_max_ratio))
      *seed = saved_seed; /* restore the original seed */
              
      ifelse($4, _x, `
      /* varying extra precs */
      for (prec_val = 0; prec_val <= 2; prec_val++) {')
      SET_INTERNAL_PARAMS($1_type, $4)
           
          /* values near underflow, 1, or overflow */
          for (norm = -1; norm <= 1; norm++) {

            /* number of tests */
            for (i=0; i<ntests; i++){

              /* generate test inputs */                
             ifelse(`$1&&$1', `$2&&$3', `TESTGEN_CASE1($1, $2, $3)',
                     `$1', `c', `TESTGEN_CASE4($1, $2, $3)',
                     `$1', `z', `TESTGEN_CASE3($1, $2, $3)',
                     `TESTGEN_CASE2($1, $2, $3)')
        
                count++;
        
        
                /* varying incx */
                for (incx_val = -2; incx_val <= 2; incx_val++){
                  if (incx_val == 0) continue;         

                  /* setting incx */
                  incx = incx_val;
                  INC_ADJUST(incx, $2_type)               

                  /* set x starting index */
                  ix=0;
                  if (incx < 0) ix = -(n-1)*incx;
                
                  /* copy x_gen to x */
                  for(j=0 ; j<n*incx_gen; j+=incx_gen){
                    GET_ARRAY_ELEMENT(x_genj, $2_type, x_gen, $2_type, j)
                    SET_ARRAY_ELEMENT(x_genj, $2_type, x, $2_type, ix)
                    ix += incx;
                  } 

                  /* varying incy */
                  for (incy_val = -2; incy_val <= 2; incy_val++){
                    if (incy_val == 0) continue;

                    /* setting incy */
                    incy = incy_val;
                    INC_ADJUST(incy, $3_type)           

                    /* set y starting index */
                    iy=0;
                    if (incy < 0) iy = -(n-1)*incy;

                    /* copy y_gen to y */
                    for(j=0 ; j<n*incy_gen; j+=incy_gen){
                      GET_ARRAY_ELEMENT(y_genj, $3_type, y_gen, $3_type, j)
                      SET_ARRAY_ELEMENT(y_genj, $3_type, y, $3_type, iy)
                      iy += incy;
                    } 
                  
                    /* varying incw */
                    for (incw_val = -2; incw_val <= 2; incw_val++){
                      if (incw_val == 0) continue;

                      /* setting incw */
                      incw = incw_val;
                      INC_ADJUST(incw, $1_type)         

                      /* For the sake of speed, we throw out this case at random */
                      if ( xrand(seed) >= test_prob ) continue;

                      /* call WAXPBY_NAME($1, $2, $3, $4) to get w */
                      FPU_FIX_STOP;
                      WAXPBY_NAME($1, $2, $3, $4)(n, alpha, x, incx_val, beta, y, incy_val,
                                                  w, incw_val ifelse(`$4', `_x', `, prec'));
                      FPU_FIX_START;

                      /* computing the ratio */
                      ifelse(`$1', `$3', `TEST_CASE1($1, $2, $3)',
                             `ifelse(`$1', `$2', `TEST_CASE2($1, $2, $3)',
                                                 `TEST_CASE3($1, $2, $3)')')

                    /* Increase the number of bad ratio, if the ratio
                       is bigger than the threshold.
                       The !<= below causes NaN error to be detected.
                       Note that (NaN > thresh) is always false. */
                      if ( !(ratio <= thresh) ) {
                        bad_ratios++;
                  
                        if ((debug == 3) &&        /* print only when debug is on */
                           (count != old_count) && /* print if old vector is different 
                                                      from the current one */ 
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

                          printf("incx=%d, incy=%d, incw=%d:\n", incx, incy, incw);
                        
                          ix=0; iy=0; iw=0;
                          if (incx < 0) ix = -(n-1)*incx;
                          if (incy < 0) iy = -(n-1)*incy;
                          if (incw < 0) iw = -(n-1)*incw;
  
                          for (j=0; j<n; j++){
                            printf("      "); PRINT_ARRAY_ELEM(x, ix, $2_type) printf("; ");
                            PRINT_ARRAY_ELEM(y, iy, $3_type) printf("; ");
                            PRINT_ARRAY_ELEM(w, iw, $1_type) printf("; ");
                            ix += incx; iy += incy; iw += incw;
                          }

                          printf("      "); PRINT_VAR(alpha, $1_type) printf("; "); 
                          PRINT_VAR(beta, $1_type) printf("\n");
                          printf("      ratio=%.4e\n", ratio);
                          p_count++;
                        }
                      }
                      if (d_count == 0) {

                        if (ratio > ratio_max)
                          ratio_max = ratio;

                        if (ratio != 0.0 && ratio < ratio_min)
                          ratio_min = ratio;

                        tot_tests++;
                      }
                    } /* incw */
                  } /* incy */
                } /* incx */
            } /* tests */
          } /* norm */
ifelse(`$4', `_x', `         } /* prec */')
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
  }

  FREE_VECTOR(x, $2_type)
  FREE_VECTOR(y, $3_type)
  FREE_VECTOR(w, $1_type)
  FREE_VECTOR(w_true, EXTRA_TYPE($1_type))
  FREE_VECTOR(x_gen, $2_type)
  FREE_VECTOR(y_gen, $3_type)
  FREE_VECTOR(temp_ab, $1_type)
  FREE_VECTOR(temp_xy, $2_type)

  *min_ratio = ratio_min;
  *num_bad_ratio = bad_ratios;
  *num_tests = tot_tests;
  FPU_FIX_STOP;
  return ratio_max;
}' )dnl
dnl
dnl
dnl--------------------------------------------------------------------------
dnl Usage TESTGEN_CASE1(abw_typeltr, x_typeltr, y_typeltr)
dnl typeltr can be:
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl--------------------------------------------------------------------------
define(`TESTGEN_CASE1', 
`       DOT_TESTGEN_NAME($1, $2, $2)(1, 0, 1, norm, blas_no_conj,
                       &alpha, 0, &beta, 0,
                       &x_fix1, &x_gen[0], seed,
                                    &y_gen[0], &HEAD(w_true)[0], &TAIL(w_true)[0]);
        
        xgen_val = incx_gen;
        ygen_val = incy_gen;
        for ( wgen_val = incw_gen; wgen_val < n*incw_gen; wgen_val+=incw_gen ) {
           DOT_TESTGEN_NAME($1, $2, $2)(1, 0, 1, norm, blas_no_conj,
                                        &alpha, 1, &beta, 1,
                                        &x_fix1, &x_gen[xgen_val], seed,
                                        &y_gen[ygen_val], &HEAD(w_true)[wgen_val], &TAIL(w_true)[wgen_val]);
           xgen_val+=incx_gen;
           ygen_val+=incy_gen;
        }')dnl
dnl
dnl
dnl
dnl--------------------------------------------------------------------------
dnl Usage TESTGEN_CASE2(abw_typeltr, x_typeltr, y_typeltr)
dnl typeltr can be:
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl--------------------------------------------------------------------------
define(`TESTGEN_CASE2', 
`       X = xrand(seed);      
        X_int = X * (power(2,12)-1);
        X = X_int;

        alpha = X*X*X*X / power(2,48);
        beta = (X*X+X+1)*(X*X-X+1) / power(2,48);

        x_i = X*X / power(2,24);
        y_i = -(X*X-1) / power(2,24);
        
        xgen_val = 0;
        ygen_val = 0;
        for ( wgen_val = 0; wgen_val < n*incw_gen; wgen_val+=incw_gen ) {
           x_gen[xgen_val] = x_i;
           y_gen[ygen_val] = y_i;
           HEAD(w_true)[wgen_val] = 1.0 / power(2,72);
           TAIL(w_true)[wgen_val] = 0.0;
           xgen_val += incx_gen;
           ygen_val += incy_gen;
        }')dnl
dnl
dnl
dnl
dnl
dnl--------------------------------------------------------------------------
dnl Usage TESTGEN_CASE3(abw_typeltr, x_typeltr, y_typeltr)
dnl typeltr can be:
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl--------------------------------------------------------------------------
define(`TESTGEN_CASE3', 
`       X = xrand(seed);
        X_int = X * (power(2,12)-1);
        X = X_int;

        ifelse(`$2', `c', 
                        alpha[0] =  X*X*X*X / power(2,48);
                        alpha[1] =  X*X*X*X / power(2,48);
                        `x_i[0] = 0.0;
                         x_i[1] = X*X / power(2,24);',
               `$2', `z',
                        alpha[0] =  X*X*X*X / power(2,48);
                        alpha[1] =  X*X*X*X / power(2,48);
                        `x_i[0] = 0.0;
                         x_i[1] = X*X / power(2,24);',
                  
                        `alpha[0] =  - X*X*X*X / power(2,48);
                         alpha[1] =  X*X*X*X / power(2,48);
                         x_i = X*X / power(2,24);')

        ifelse(`$3', `c', 
                        `beta[0] = (X*X+X+1)*(X*X-X+1) / power(2,48);
                         beta[1] = (X*X+X+1)*(X*X-X+1) / power(2,48);
                         y_i[0] = 0.0;
                         y_i[1] = -(X*X-1) / power(2,24);',
               `$3', `z', 
                        `beta[0] = (X*X+X+1)*(X*X-X+1) / power(2,48);
                         beta[1] = (X*X+X+1)*(X*X-X+1) / power(2,48);
                         y_i[0] = 0.0;
                         y_i[1] = -(X*X-1) / power(2,24);',

                        `beta[0] = - (X*X+X+1)*(X*X-X+1) / power(2,48);
                         beta[1] = (X*X+X+1)*(X*X-X+1) / power(2,48);
                         y_i = -(X*X-1) / power(2,24);')

       
        xgen_val = 0;
        ygen_val = 0;
        for ( wgen_val = 0; wgen_val < n*incw_gen; wgen_val+=incw_gen ) {
           SET_ARRAY_ELEMENT(x_i, $2_type, x_gen, $2_type, xgen_val)
           SET_ARRAY_ELEMENT(y_i, $3_type, y_gen, $3_type, ygen_val)  
           HEAD(w_true)[wgen_val] = - 1.0 / power(2,72);
           HEAD(w_true)[wgen_val+1] = 1.0 / power(2,72);
           TAIL(w_true)[wgen_val] = 0.0;
           TAIL(w_true)[wgen_val+1] = 0.0;
           xgen_val += incx_gen;
           ygen_val += incy_gen;
        }')dnl
dnl
dnl
dnl
dnl--------------------------------------------------------------------------
dnl Usage TESTGEN_CASE4(abw_typeltr, x_typeltr, y_typeltr)
dnl typeltr can be:
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl--------------------------------------------------------------------------
define(`TESTGEN_CASE4',
`       x_fix1_temp = 1.0;
        BLAS_sdot_testgen(1, 0, 1, norm, blas_no_conj,
                       &atemp, 0, &btemp, 0,
                       &x_fix1_temp, &xtemp, seed,
                       &ytemp, &wltemp, &wttemp);
        ifelse(`$2',`c',`x_gen[0] = 0.0;
                         x_gen[1] = xtemp;
                         alpha[0] = atemp;
                         alpha[1] = atemp;'
               ,`x_gen[0] = xtemp;
                 alpha[0] = -atemp;
                 alpha[1] = atemp;')
        
        ifelse(`$3',`c',`y_gen[0] = 0.0;
                         y_gen[1] = ytemp;
                         beta[0] = btemp;
                         beta[1] = btemp;'
               ,`y_gen[0] = ytemp;
                 beta[0] = -btemp;
                 beta[1] = btemp;')
        
        HEAD(w_true)[0] = -wltemp;
        HEAD(w_true)[1] = wltemp;
        TAIL(w_true)[0] = 0.0;
        TAIL(w_true)[1] = 0.0;      
        
        xgen_val = incx_gen;
        ygen_val = incy_gen;
        for ( wgen_val = incw_gen; wgen_val < n*incw_gen; wgen_val+=incw_gen ) {
           BLAS_sdot_testgen(1, 0, 1, norm, blas_no_conj,
                          &atemp, 1, &btemp, 1,
                          &x_fix1_temp, &xtemp, seed,
                          &ytemp, &wltemp, &wttemp);
           
           ifelse(`$2',`c',`x_gen[xgen_val] = 0;
                            x_gen[xgen_val+1] = xtemp;'
                  ,`x_gen[xgen_val] = xtemp;')
        
           ifelse(`$3',`c',`y_gen[ygen_val] = 0;
                            y_gen[ygen_val+1] = ytemp;'
                            
                  ,`y_gen[ygen_val] = ytemp;')
           
        
           HEAD(w_true)[wgen_val] = -wltemp;
           HEAD(w_true)[wgen_val+1] = wltemp;
           TAIL(w_true)[wgen_val] = 0.0;
           TAIL(w_true)[wgen_val+1] = 0.0;
           xgen_val+=incx_gen;
           ygen_val+=incy_gen;
        }')dnl
dnl
dnl
dnl     
dnl--------------------------------------------------------------------------
dnl Usage TEST_CASE1(abw_typeltr, x_typeltr, y_typeltr)
dnl typeltr can be:
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl--------------------------------------------------------------------------
define(`TEST_CASE1',
`       ix = 0;
        if (incx < 0) ix = -(n-1)*incx;
        iy = 0;
        if (incy < 0) iy = -(n-1)*incy;
        iw = 0;
        if (incw < 0) iw = -(n-1)*incw;
        ratio = 0.0;

        for ( test_val = 0; test_val < n*incw_gen; test_val+=incw_gen ) {
           TEST_DOT_NAME($1, $2, $2, $4) TEST_CASE1_ARGS($1, $2, $2, $4)
           ix += incx;
           iy += incy;
           iw += incw;
           if (MAX(ratio, new_ratio) == new_ratio) {
                iymax = iy - incy;
                ixmax = ix - incx;
           }
           ratio = MAX(ratio, new_ratio);
        }')dnl
dnl
dnl
dnl
dnl
dnl--------------------------------------------------------------------------
dnl Usage TEST_CASE1_ARGS(abw_typeltr, x_typeltr, y_typeltr)
dnl typeltr can be:
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl--------------------------------------------------------------------------
define(`TEST_CASE1_ARGS', `ifelse(

        `$1', `c', `(1, blas_no_conj, alpha, beta, 
                    &y[iy], &w[iw], 
                    &HEAD(w_true)[test_val], &TAIL(w_true)[test_val],
                    &x_fix1, incx, &x[ix], incx,
                    eps_int, un_int, &new_ratio);',

        `$1', `z', `(1, blas_no_conj, alpha, beta, 
                    &y[iy], &w[iw], 
                    &HEAD(w_true)[test_val], &TAIL(w_true)[test_val],
                    &x_fix1, incx, &x[ix], incx,
                    eps_int, un_int, &new_ratio);',
        
        `(1, blas_no_conj, alpha, beta, 
        y[iy], w[iw],
        HEAD(w_true)[test_val], TAIL(w_true)[test_val], 
        &x_fix1, incx, &x[ix], incx,
        eps_int, un_int, &new_ratio);')')dnl
dnl
dnl
dnl
dnl--------------------------------------------------------------------------
dnl Usage TEST_CASE2(abw_typeltr, x_typeltr, y_typeltr)
dnl typeltr can be:
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl--------------------------------------------------------------------------
define(`TEST_CASE2',
`       ix = 0;
        if (incx < 0) ix = -(n-1)*incx;
        iy = 0;
        if (incy < 0) iy = -(n-1)*incy;
        iw = 0;
        if (incw < 0) iw = -(n-1)*incw;
        ratio = 0.0;

        for ( test_val = 0; test_val < n*incw_gen; test_val+=incw_gen ) {
           TEST_DOT_NAME($1, $3, $3, $4) TEST_CASE2_ARGS($1, $3, $3, $4)
           ix += incx;
           iy += incy;
           iw += incw;
           if (MAX(ratio, new_ratio) == new_ratio) {
                iymax = iy - incy;
                ixmax = ix - incx;
           }
           ratio = MAX(ratio, new_ratio);
        }')dnl
dnl
dnl
dnl
dnl--------------------------------------------------------------------------
dnl Usage TEST_CASE2_ARGS(abw_typeltr, x_typeltr, y_typeltr)
dnl typeltr can be:
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl--------------------------------------------------------------------------
define(`TEST_CASE2_ARGS', `ifelse(

        `$1', `c', `(1, blas_no_conj, beta, alpha, 
                                         &x[ix], &w[iw], 
                                         &HEAD(w_true)[test_val], &TAIL(w_true)[test_val], 
                                         &x_fix2, incy, &y[iy], incy, 
                                         eps_int, un_int, &new_ratio);',
        
        `$1', `z', `(1, blas_no_conj, beta, alpha, 
                                         &x[ix], &w[iw], 
                                         &HEAD(w_true)[test_val], &TAIL(w_true)[test_val], 
                                         &x_fix2, incy, &y[iy], incy, 
                                         eps_int, un_int, &new_ratio);',
        
        `(1, blas_no_conj, beta, alpha, 
        x[ix], w[iw],
        HEAD(w_true)[test_val], TAIL(w_true)[test_val], 
        &x_fix2, incy, &y[iy], incy, 
        eps_int, un_int, &new_ratio);')')dnl
dnl
dnl
dnl
dnl
dnl--------------------------------------------------------------------------
dnl Usage TEST_CASE3(abw_typeltr, x_typeltr, y_typeltr)
dnl typeltr can be:
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl--------------------------------------------------------------------------
define(`TEST_CASE3',
`       ix = 0;
        if (incx < 0) ix = -(n-1)*incx;
        iy = 0;
        if (incy < 0) iy = -(n-1)*incy;
        iw = 0;
        if (incw < 0) iw = -(n-1)*incw;
        ratio = 0.0;
                
        SET_VECTOR_ELEMENT(temp_ab, 0, alpha, $1_type)
        SET_VECTOR_ELEMENT(temp_ab, incw_gen, beta, $1_type)    
 
        for ( test_val = 0; test_val < n*incw_gen; test_val+=incw_gen ) {
           GET_ARRAY_ELEMENT(x_genj, $2_type, x, $2_type, ix)
           SET_ARRAY_ELEMENT(x_genj, $2_type, temp_xy, $2_type, 0)

           GET_ARRAY_ELEMENT(y_genj, $3_type, y, $3_type, iy)
           SET_ARRAY_ELEMENT(y_genj, $3_type, temp_xy, $3_type, incy_gen)

           TEST_DOT_NAME($1, $1, $2, $4) TEST_CASE3_ARGS($1, $1, $2, $4)
           if (MAX(ratio, new_ratio) == new_ratio) {
                iymax = iy;
                ixmax = ix;
           }
           ratio = MAX(ratio, new_ratio);

           ix += incx;
           iy += incy;
           iw += incw;
        }')dnl
dnl
dnl
dnl
dnl
dnl--------------------------------------------------------------------------
dnl Usage TEST_CASE3_ARGS(abw_typeltr, x_typeltr, y_typeltr)
dnl typeltr can be:
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl--------------------------------------------------------------------------
define(`TEST_CASE3_ARGS', `ifelse(
        `$1', `c', `(2, blas_no_conj, one, zero, 
                                         dummy, &w[iw], 
                                         &HEAD(w_true)[test_val], &TAIL(w_true)[test_val], 
                                         &temp_ab[0], 1, &temp_xy[0], 1, 
                                         eps_int, un_int, &new_ratio);',

        `$1', `z', `(2, blas_no_conj, one, zero, 
                                         dummy, &w[iw], 
                                         &HEAD(w_true)[test_val], &TAIL(w_true)[test_val], 
                                         &temp_ab[0], 1, &temp_xy[0], 1, 
                                         eps_int, un_int, &new_ratio);',

        `(2, blas_no_conj, one, zero, 
                                         dummy, w[iw], 
                                         HEAD(w_true)[test_val], TAIL(w_true)[test_val], 
                                         temp_ab, 1, temp_xy, 1, 
                                         eps_int, un_int, &new_ratio);')')dnl
dnl
dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: DO_TEST_WAXPBY_NAME(abw_typeltr, x_typeltr, y_typeltr, extended)
dnl        create do_test_waxpby name
dnl -------------------------------------------------------------------
define(`DO_TEST_WAXPBY_NAME', `ifelse(
        `$1&&$1', `$2&&$3',
        `do_test_$1waxpby$4',
        `do_test_$1waxpby_$2_$3$4')') dnl
dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: CALL_DO_TEST_WAXPBY(abw_typeltr, x_typeltr, y_typeltr, extended)
dnl
dnl        abw_type : the type and precision of alpha, beta and w
dnl        x_type   : the type and precision of x
dnl        y_type   : the type and precision of y
dnl        extended : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl -------------------------------------------------------------------
define(`CALL_DO_TEST_WAXPBY',           
        `  fname = "WAXPBY_NAME($1, $2, $3, $4)";
           printf("Testing %s...\n", fname);
           min_ratio = 1e308; max_ratio = 0.0;
           total_bad_ratios = 0; total_tests = 0;
           for(n=0; n<=nsizes; n++){
            
              total_max_ratio = DO_TEST_WAXPBY_NAME($1, $2, $3, $4)(n, ntests, &seed, thresh, debug, test_prob, 
                                            &total_min_ratio, &num_bad_ratio, &num_tests);
              if (total_max_ratio > max_ratio)
                max_ratio = total_max_ratio;

              if (total_min_ratio < min_ratio)
                min_ratio = total_min_ratio;

              total_bad_ratios += num_bad_ratio;
              total_tests += num_tests; 
           }

           nr_routines++;
           if (total_bad_ratios == 0)
             printf("PASS> ");
           else{
             printf("FAIL> ");
             nr_failed_routines++;
           }

           printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n",
               fname, total_bad_ratios, total_tests, max_ratio);
')
dnl
dnl
dnl
FOREACH(`WAXPBY_ARGS', `
DO_TEST_WAXPBY(arg)')dnl

MAIN(`', `

FOREACH(`WAXPBY_ARGS', `
CALL_DO_TEST_WAXPBY(arg)')')dnl

