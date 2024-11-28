dnl Generates test code for tpmv
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"
include(cblas.m4)dnl
include(test-common.m4)dnl
include(tpmv/tpmv-common.m4)dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_TPMV_COMMENT(extended)
dnl        if extended, then print info about prec loop in the code 
dnl        structure
dnl -----------------------------------------------------------------------
define(`DO_TEST_TPMV_COMMENT',`
/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on TPMV.
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
ifelse(`$1', `_x',` *        prec loop   -- varying internal prec: single, double, or extra', ` *')
 *        norm loop   -- varying norm: near undeflow, near one, or 
 *                    -- near overflow
 *          order loop  -- varying order type: rowmajor or colmajor
 *            uplo loop   -- varying uplo type: upper or lower
 *              trans loop  -- varying trans type: no trans, trans, and conj trans  
 *                diag loop   -- varying diag type: non unit, and unit
 *                  numtest loop  -- how many times the test is perform with 
 *                                -- above set of attributes
 *                        incx loop     -- varying incx: -2, -1, 1, 2
 */')dnl
dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_TPMV(ax_typeltr, tp_typeltr, extended)
dnl
dnl        ax_typeltr : the type and precision of alpha, and x
dnl        tp_typeltr  : the type and precision of the matrix T
dnl        extended   : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ----------------------------------------------------------------------
define(`DO_TEST_TPMV', `ifelse(
        `$1&&$2', `$1&&$1',`
double do_test_$1tpmv$3(int n,
                     int ntests,
                     int *seed,
                     double thresh, 
                     int debug,
                     float test_prob,
                     double *min_ratio,
                     int *num_bad_ratio,
                     int *num_tests)
DO_TEST_TPMV_COMMENT($3)
DO_TEST_TPMV_BODY($1, $2, $3) /* end of do_test_$1tpmv$3 */',`
double do_test_$1tpmv_$2$3 (int n, 
                          int ntests,
                          int *seed,
                          double thresh, 
                          int debug,
                          float test_prob,
                          double *min_ratio,
                          int *num_bad_ratio,
                          int *num_tests)
DO_TEST_TPMV_COMMENT($3)
DO_TEST_TPMV_BODY($1, $2, $3) /* end of do_test_$1tpmv_$2$3 */')') dnl
dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: DO_TEST_TPMV_BODY(ax_typeltr, tp_typeltr, extended)
dnl
dnl        ax_typeltr  : the type and precision of alpha, and x
dnl        tp_typeltr  : the type and precision of the matrix T
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ---------------------------------------------------------------------
define(`DO_TEST_TPMV_BODY',
`{                                        
  /* function name */
  const char fname[] = "TPMV_NAME($1, $2, $3)";

  /* max number of debug lines to print */
  const int max_print = 3;

  /* Variables in the "x_val" form are loop vars for corresponding
     variables */
  int i;             /* iterate through the repeating tests */
  int j;             /* multipurpose counters */
  int ix;            /* use to index x  */
  int incx_val, incx;/* for testing different inc values */
  int inctp, inc_index;/* 1 if real, 2 if complex */   
  int d_count;       /* counter for debug */
  int find_max_ratio;/* find_max_ratio = 1 only if debug = 3 */
  int p_count;       /* counter for the number of debug lines printed*/
  int tot_tests;     /* total number of tests to be done */
  int norm;          /* input values of near underflow/one/overflow */
  double ratio_max;  /* the current maximum ratio */
  double ratio_min;  /* the current minimum ratio */
  double *ratios;    /* a temporary variable for calculating ratio */
  double ratio = 0.0;      /* the per-use test ratio from test() */
  int bad_ratios = 0;    /* the number of ratios over the threshold */
  double eps_int;    /* the internal epsilon expected--2^(-24) for float */
  double un_int;     /* the internal underflow threshold */
  DECLARE(alpha, $1_type)
  DECLARE(beta, $1_type)
  DECLARE_VECTOR(tp, $2_type)
  DECLARE_VECTOR(x, $1_type)
  DECLARE_VECTOR(x_gen, $1_type)/* used to store vectors generated by testgen; eventually copied to x */
  DECLARE_VECTOR(temp, $2_type) /* use for calculating ratio */
  
  /* the true r calculated by testgen(), in double-double */
  DECLARE_VECTOR(r_true, EXTRA_TYPE($1_type))

  int alpha_val;
  int alpha_flag = 0, beta_flag = 0;
  int order_val;
  enum blas_order_type order_type = 0;
  ifelse(`$3', `_x', `int prec_val;')
  enum blas_prec_type prec = 0;
  int uplo_val;
  enum blas_uplo_type uplo_type = 0;
  int trans_val;
  enum blas_trans_type trans_type = 0;
  int diag_val;
  enum blas_diag_type diag_type = 0;

  DECLARE(beta_zero_fake, $1_type)
  DECLARE(rin_zero_fake, $1_type)

  int saved_seed;   /* for saving the original seed */
  int count, old_count = -1;  
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



  incx = inctp = inc_index = 1;
  INC_ADJUST(incx, $1_type)
  INC_ADJUST(inctp, $2_type)
  INC_ADJUST(inc_index, $1_type)

  /* get space for calculation */
  MALLOC_VECTOR(x, $1_type, 3*n*2)
  MALLOC_VECTOR(x_gen, $1_type, 3*n*2)
  MALLOC_VECTOR(tp, $2_type, 2*n*n*n*inctp)
  MALLOC_VECTOR(temp, $2_type, 3*n*2)
  MALLOC_VECTOR(r_true, EXTRA_TYPE($1_type), n)
  MALLOC_VECTOR(ratios, real_D, n)


  /* initialization */
  saved_seed = *seed;
  ratio_min = 1e308;
  ratio_max = 0.0;
  tot_tests = 0;
  p_count = 0;
  count = 0; 
  find_max_ratio = 0;
  ZERO(beta, $1_type)
  beta_flag=1;

  ZERO(beta_zero_fake, $1_type)
  ZERO(rin_zero_fake, $1_type)

  if (debug == 3)
    find_max_ratio = 1;

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

      ifelse($3, _x, `
      /* varying extra precs */
      for (prec_val = 0; prec_val <= 2; prec_val++) {')
      SET_INTERNAL_PARAMS($1_type, $3)
           
      /* values near underflow, 1, or overflow */
      for (norm = -1; norm <= 1; norm++) {

        /* row or col major */
        for (order_val=0; order_val<2; order_val++) {
          switch(order_val){
          case 0: order_type = blas_rowmajor; break;
          case 1: order_type = blas_colmajor; break;
        }

          /* upper or lower */
          for (uplo_val=0; uplo_val<2; uplo_val++) {
            switch(uplo_val){
            case 0: uplo_type = blas_upper; break;
            case 1: uplo_type = blas_lower; break;
            }

            /* no trans, trans, or conj trans */
            for (trans_val=0; trans_val<3; trans_val++){
              switch(trans_val){
              case 0: trans_type = blas_no_trans; break;
              case 1: trans_type = blas_trans; break;
              case 2: trans_type = blas_conj_trans; break;
              }

              /* non_unit_diag, unit_diag */
              for (diag_val=0; diag_val<2; diag_val++){
                switch(diag_val){
                case 0: diag_type = blas_non_unit_diag; break;
                case 1: diag_type = blas_unit_diag; break;
                }             

                /* number of tests */
                for (i=0; i<ntests; i++){               
                  /* For the sake of speed, we throw out this case at random */
                  if ( xrand(seed) >= test_prob ) continue;
                  
                        TESTGEN_TPMV_NAME($1, $2)(norm, order_type, uplo_type, 
                                    trans_type, diag_type, n,
                                    &alpha, alpha_flag, tp, x_gen, 
                                    seed, HEAD(r_true), TAIL(r_true));

                  count++;

                  /* varying incx */
                  for (incx_val = -2; incx_val <= 2; incx_val++){
                    if (incx_val == 0) continue;               
                    
                    /* setting incx */
                    incx = incx_val;
                    INC_ADJUST(incx, $1_type)             

                    /* set x starting index */
                    ix = 0;
                    if (incx < 0) ix = -(n-1)*incx;
                
                    /* copy x_gen to x */
                    for(j = 0 ; j < n; j++){
                      COPY_VECTOR_ELEMENT(x, ix, x_gen, j*inc_index, $1_type)
                      ix += incx;
                    } 

                    /* call TPMV_NAME($1, $2, $3) */
                    FPU_FIX_STOP;
                    TPMV_NAME($1, $2, $3)(order_type, uplo_type, trans_type, diag_type,
                              n, alpha, tp,
                              x, incx_val ifelse(`$3', `_x', `, prec'));
                    FPU_FIX_START;

                    /* set x starting index */
                    ix = 0;
                    if (incx < 0) ix = -(n-1)*incx;

                    /* computing the ratio */
                    for(j=0; j<n; j++){
                      /* copy row j of tp to temp */
                      $2tpmv_copy_row(order_type, uplo_type, trans_type, n, tp, temp, j);

                      TEST_DOT_NAME($1, $2, $1)(n, blas_no_conj,
                                    alpha, beta_zero_fake, rin_zero_fake, 
                                    VECTOR_ELEMENT(x, ix, $1_type), 
                                    VECTOR_ELEMENT(r_true, j*inc_index, EXTRA_TYPE($1_type)), 
                                    temp, 1, x_gen, 1, eps_int, un_int,
                                    &ratios[j]);

                      /* take the max ratio */
                      if (j==0)
                        ratio = ratios[0];
                      else if (ratios[j] > ratio)
                        ratio = ratios[j];
                                
                      ix += incx;                       
                    }

                    /* increase the number of bad ratio, if the ratio is 
                       bigger than the threshold */
                    if (ratio > thresh) {
                      bad_ratios++;
                     
                      if ((debug==3) &&          /* print only when debug is on */
                          (count != old_count) && /* print if old vector is different 
                                                     from the current one */ 
                          (d_count == find_max_ratio) &&  
                          (p_count <= max_print) &&
                          (ratio > 0.5*ratio_max))
                      {
                        p_count++;
                        old_count = count;      

                        printf("FAIL> %s: n = %d, ntests = %d, threshold = %4.2f,\n", 
                                fname, n, ntests, thresh);
        
                        /* Print test info */
                        PRINT_PREC(prec)
                        PRINT_NORM(norm)
                        PRINT_ORDER(order_type)
                        PRINT_UPLO(uplo_type)
                        PRINT_TRANS(trans_type)
                        PRINT_DIAG(diag_type)

                        printf(" incx=%d:\n", incx);
                        ix=0;
                        if (incx < 0) ix = -(n-1)*incx;
                        
                        printf("      TP=");     
                        for (j=0; j<n; j++){
                          /* copy row j of tp to temp */
                          $2tpmv_copy_row(order_type, uplo_type, trans_type, n, tp, temp, j);
                         
                          if (j>0) printf("        ");
                          $2print_vector(temp, n, 1, NULL);
                        }

                        ix=0;
                        if (incx < 0) ix = -(n-1)*incx;
                        
                        for (j=0; j<n; j++){
                          printf("      "); 
                          PRINT_ARRAY_ELEM(x_gen, j*inc_index, $1_type, x) printf("; ");
                          PRINT_ARRAY_ELEM(x, ix, $1_type, x_final) printf("\n");
                          ix += incx;
                        }

                        printf("      "); PRINT_VAR(alpha, $1_type) printf("; "); printf("\n");
                        for (j=0; j<n; j++){
                          printf("      ");

                          PRINT_ARRAY_ELEM(r_true, j*inc_index, EXTRA_TYPE($1_type))
                          printf(", ratio[%d]=%.4e\n", j, ratios[j]); 
                        }
                        printf("      ratio=%.4e\n", ratio);
                      }
                    }
                    if (d_count==0) {
                      if (ratio > ratio_max)
                        ratio_max = ratio;
                      if (ratio != 0.0 && ratio < ratio_min)
                        ratio_min = ratio;
                      tot_tests++;
                    }
                  } /* incx */

                } /* numtests */
              } /* diag */
            } /* trans */
          } /* uplo */
        } /* order */
      } /* norm */
ifelse(`$3', `_x', `       } /* prec */')
    } /* alpha */
  } /* debug */

  if ((debug == 2) ||
     ((debug == 1) && bad_ratios > 0)){
    printf("      %s:  n = %d, ntests = %d, thresh = %4.2f\n", 
            fname, n, ntests, thresh);
    printf("      bad/total = %d/%d=%3.2f, min_ratio = %.4e, max_ratio = %.4e\n\n",
            bad_ratios, tot_tests,
           ((double)bad_ratios)/((double)tot_tests), ratio_min, ratio_max);
  }

  FPU_FIX_STOP;

  FREE_VECTOR(x, $1_type)
  FREE_VECTOR(x_gen, $1_type)
  FREE_VECTOR(temp, $2_type)
  FREE_VECTOR(tp, $2_type)
  FREE_VECTOR(r_true, EXTRA_TYPE($1_type))
  FREE_VECTOR(ratios, real_D)

  *min_ratio = ratio_min;
  *num_bad_ratio = bad_ratios;
  *num_tests = tot_tests;
  return ratio_max;
}' )dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: TESTGEN_TPMV_NAME(ax_typeltr, tp_typeltr) 
dnl        create a testgen_trmv name 
dnl --------------------------------------------------------------------
define(`TESTGEN_TPMV_NAME', `ifelse(
        `$1&&$2', `$1&&$1', `BLAS_$1tpmv_testgen',
                            `BLAS_$1tpmv_$2_testgen')')dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: DO_TEST_TPMV_NAME(abx_typeltr, tp_typeltr, extended)
dnl        create do_test_dot name
dnl -------------------------------------------------------------------
define(`DO_TEST_TPMV_NAME', `ifelse(
        `$1&&$2', `$1&&$1',
        `do_test_$1tpmv$3',
        `do_test_$1tpmv_$2$3')') dnl
dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: CALL_DO_TEST_TPMV(ax_typeltr, tp_typeltr, extended)
dnl
dnl        ax_type  : the type and precision of alpha and x
dnl        tp_type   : the type and precision of tp
dnl        extended : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl -------------------------------------------------------------------
define(`CALL_DO_TEST_TPMV',           
        `  fname = "TPMV_NAME($1, $2, $3, $4)";
           printf("Testing %s...\n", fname);
           min_ratio = 1e308; max_ratio = 0.0;
           total_bad_ratios = 0; total_tests = 0;
           for(i=0; i<nsizes; i++){
            
              n = sizes[i];
              total_max_ratio = DO_TEST_TPMV_NAME($1, $2, $3)(n, ntests, &seed, thresh, debug,
                                            test_prob, 
                                            &total_min_ratio, &num_bad_ratio, &num_tests);
              if (total_max_ratio > max_ratio)
                max_ratio = total_max_ratio;

              if (total_min_ratio != 0.0 && total_min_ratio < min_ratio)
                min_ratio = total_min_ratio;

              total_bad_ratios += num_bad_ratio;
              total_tests += num_tests; 
           }

           nr_routines++;
           if (total_bad_ratios == 0)
             printf("PASS> ");
           else{
             nr_failed_routines++;
             printf("FAIL> ");
           }

           if (min_ratio == 1e308) min_ratio = 0.0;
           printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n",
               fname, total_bad_ratios, total_tests, max_ratio);
')
dnl
FOREACH(`TPMV_ARGS', `
DO_TEST_TPMV(arg)')dnl

#define NSIZES 12
MAIN(`
  int i;
  int sizes[NSIZES]={0, 1, 2, 3, 4, 5, 6, 10, 17, 28, 37, 58};', `

FOREACH(`TPMV_ARGS', `
CALL_DO_TEST_TPMV(arg)')')dnl

