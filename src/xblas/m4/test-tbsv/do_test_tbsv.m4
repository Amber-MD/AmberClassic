dnl Generates test code for tbsv
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"
include(cblas.m4)dnl
include(test-common.m4)dnl
include(tbsv/tbsv-common.m4)dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_TBSV_COMMENT(extended)
dnl        if extended, then print info about prec loop in the code 
dnl        structure
dnl -----------------------------------------------------------------------
define(`DO_TEST_TBSV_COMMENT',`
/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on TBSV.
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
 * test_prob (input) float
 *           The specified test will be performed only if the generated 
 *           random exceeds this threshold.
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
 *                    ldt loop      -- varying ldt: n, n+1, and 2n
 *                      incx loop     -- varying incx: -2, -1, 1, 2
 */')dnl
dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_TBSV(ax_typeltr, TP_typeltr, extended)
dnl
dnl        ax_typeltr : the type and precision of alpha, and x
dnl        TP_typeltr  : the type and precision of the matrix A
dnl        extended   : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ----------------------------------------------------------------------
define(`DO_TEST_TBSV', `ifelse(
        `$1&&$2', `$1&&$1',`
double do_test_$1tbsv$3(int n,
                     int ntests,
                     int *seed,
                     double thresh, 
                     int debug,
                     float test_prob,
                     double *min_ratio,
                     int *num_bad_ratio,
                     int *num_tests)
DO_TEST_TBSV_COMMENT($3)
DO_TEST_TBSV_BODY($1, $2, $3) /* end of do_test_$1tbsv$3 */',`
double do_test_$1tbsv_$2$3 (int n, 
                          int ntests,
                          int *seed,
                          double thresh, 
                          int debug,
                          float test_prob,
                          double *min_ratio,
                          int *num_bad_ratio,
                          int *num_tests)
DO_TEST_TBSV_COMMENT($3)
DO_TEST_TBSV_BODY($1, $2, $3) /* end of do_test_$1tbsv_$2$3 */')') dnl
dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: DO_TEST_TBSV_BODY(ax_typeltr, TP_typeltr, extended)
dnl
dnl        ax_typeltr  : the type and precision of alpha, and x
dnl        TP_typeltr  : the type and precision of the matrix A
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl ---------------------------------------------------------------------
define(`DO_TEST_TBSV_BODY',
`{                                        
  /* function name */
  const char fname[] = "TBSV_NAME($1, $2, $3)";

  /* max number of debug lines to print */
  /*const int max_print = 3;*/

  /* Variables in the "x_val" form are loop vars for corresponding
     variables */
  int i;            /* iterate through the repeating tests */
  int j;         /* multipurpose counters */
  int ix;       /* use to `index' x respectively */
  int ldt_val, ldt; /* for testing different values for ldt */
  int k_val, k;
  int incx_val, incx, incx_unadj; /* for testing different inc values */
  int incx_gen=1, incx_gen_unadj=1, incT=1; /* 1 if real, 2 if complex */   
  int d_count;      /* counter for debug */
  int find_max_ratio; /* find_max_ratio = 1 only if debug = 3 */
  int p_count;      /* counter for the number of debug lines printed*/
  int tot_tests;    /* total number of tests to be done */
  int norm;         /* input values of near underflow/one/overflow */
  double ratio_max; /* the current maximum ratio */
  double ratio_min; /* the current minimum ratio */
  double *ratios;   /* a temporary variable for calculating ratio */
  double ratio;     /* the per-use test ratio from test() */
  int bad_ratios=0;   /* the number of ratios over the threshold */
  double eps_int;   /* the internal epsilon expected--2^(-24) for float */
  double un_int;    /* the internal underflow threshold */
  DECLARE(alpha, $1_type)
  DECLARE(beta, $1_type)
  DECLARE_VECTOR(T, $2_type)
  DECLARE_VECTOR(x, $1_type)
  DECLARE_VECTOR(x_gen, $1_type)
  DECLARE_VECTOR(temp, $2_type) /* use for calculating ratio */
  
  /* x_gen and y_gen are used to store vectors generated by testgen.
     they eventually are copied back to x and y */
  DECLARE_VECTOR(x_gen, $3_type)

  /* the true r calculated by testgen(), in double-double */
  DECLARE_VECTOR(r_true, EXTRA_TYPE($1_type))
  int alpha_val;
  int alpha_flag, beta_flag;/* input flag for TESTGEN_TBSV_NAME($1, $2, $3) */
  int order_val;
  enum blas_order_type order_type;
  ifelse($3, `_x', 
        `int prec_val;')
  enum blas_prec_type prec;
  int uplo_val;
  enum blas_uplo_type uplo_type;
  int trans_val;
  enum blas_trans_type trans_type;
  int diag_val;
  enum blas_diag_type diag_type;
  int row; 
  int saved_seed;   /* for saving the original seed */
  int test_count = 0, count, old_count= -1;  /* use for counting the number of testgen calls * 2 */ 
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

  INC_ADJUST(incx_gen, $1_type)
  INC_ADJUST(incT, $2_type)

  /* get space for calculation */
  MALLOC_VECTOR(x, $1_type, n*4)
  MALLOC_VECTOR(x_gen, $1_type, n*2)
  MALLOC_VECTOR(T, $2_type, 16*n*n)
  MALLOC_VECTOR(temp, $2_type, n)
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
          default:
          case 0: order_type = blas_rowmajor; break;
          case 1: order_type = blas_colmajor; break;
        }

          /* upper or lower */
          for (uplo_val=0; uplo_val<2; uplo_val++) {
            switch(uplo_val){
            default:
            case 0: uplo_type = blas_upper; break;
            case 1: uplo_type = blas_lower; break;
            }

            /* no trans, trans, or conj trans */
            for (trans_val=0; trans_val<3; trans_val++){
              switch(trans_val){
              default:
              case 0: trans_type = blas_no_trans; break;
              case 1: trans_type = blas_trans; break;
              case 2: trans_type = blas_conj_trans; break;
              }

              /* non_unit_diag, unit_diag */
              for (diag_val=0; diag_val<2; diag_val++){
                switch(diag_val){
                default:
                case 0: diag_type = blas_non_unit_diag; break;
                case 1: diag_type = blas_unit_diag; break;
                }             

                /* number of tests */
                for (i=0; i<ntests; i++){               

                for (ldt_val=0; ldt_val<3; ldt_val++){  
                  switch(ldt_val){
                  case 0: ldt = n; break;
                  case 1: ldt = n+1; break;
                  case 2: ldt = 2*n; break;
                  default: ldt = n; break;
                  }

                
                /* For the sake of speed, we throw out this case at random */
                if ( xrand(seed) >= test_prob ) continue;

                for(k_val=0; k_val<ldt-1 && k_val<n; k_val++) {

                        /* For the sake of speed, we throw out this case
                           at random.  We put this in a second time
                           to make the k testing go much faster when prob < 1 and n is large*/
                        if ( xrand(seed) >= test_prob + (2.0/n) ) continue;
        
                        k = k_val;
                        
                        for (row=0; row<n; row++){
                          TESTGEN_TBSV_NAME($1, $2)(norm, order_type, uplo_type, 
                                            trans_type, diag_type, n, k, 0 /*randomize */,
                                            &alpha, alpha_flag, T, ldt, x_gen, 
                                            seed, HEAD(r_true), TAIL(r_true), row, prec);
        
                          count++;
        
                          /* varying incx */
                          for (incx_val = -2; incx_val <= 2; incx_val++){
                            if (incx_val == 0) continue;               
        
                            
                            /* setting incx */
                            incx = incx_unadj = incx_val;
                            INC_ADJUST(incx, $1_type)             
        
                            /* set x starting index */
                            ix=0;
                            if (incx < 0) ix= -(n-1)*incx;
                        
                            /* copy x_gen to x */
                            $1copy_vector(x_gen, n, incx_gen_unadj, x, incx_unadj); 

                            /* call TBSV_NAME($1, $2, $3, $4) */
                            FPU_FIX_STOP;
                            TBSV_NAME($1, $2, $3)(order_type, uplo_type, trans_type, diag_type,
                                      n, k, alpha, T, ldt,
                                      x, incx_unadj ifelse(`$3', `_x', `, prec'));
                            FPU_FIX_START;
                                test_count++;                       
        
                            ix = 0;
                            if (incx < 0) ix= -(n-1)*incx;
        
                            /* computing the ratio */
                            ratio = 0.0;
                            for(j=0; j<n; j++){
                              if (j == row){
                                ifelse(`$1', `c', `float minus_one[2] = {-1.0, 0.0};',
                                        `$1', `z', `double minus_one[2] = {-1.0, 0.0};',
                                        `int minus_one = -1.0;')
                                /* copy row j of T to temp */
                                $2tbsv_copy(order_type, uplo_type, trans_type, n, k, T, ldt, temp, j);
        
                                        /* zero out the row element, not included in the dot product */
                                SET_ZERO_VECTOR_ELEMENT(temp, `(row * incT)', $2_type)
                                        /* in doing this we rely that temp[row*incT] == 1.0 */
                                TEST_DOT_NAME($1, $1, $2)(n, blas_no_conj, minus_one, alpha, 
                                  VECTOR_ELEMENT(x_gen, j*incx_gen, $1_type), VECTOR_ELEMENT(x, ix, $1_type), 
                                  VECTOR_ELEMENT(r_true, j*incx_gen, EXTRA_TYPE($1_type)), 
                                  x, incx_unadj, temp, 1, eps_int, un_int, &ratios[j]);
                              }
                              else{
                        ifelse(`$1', `s',
        `                       double  eps_out = power(2,-BITS_S);',
                               `$1', `c',
        `                       double  eps_out = power(2,-BITS_S);',
                               `$1', `d',
        `                       double  eps_out = power(2,-BITS_D);',
                               `$1', `z',
        `                       double  eps_out = power(2,-BITS_D);')
        
                        ifelse(`$1&&$2', `c&&c', `COMPUTE_PURE_COMPLEX_RATIO',
                               `$1&&$2', `z&&z', `COMPUTE_PURE_COMPLEX_RATIO',
                               `$1&&$2', `z&&c', `COMPUTE_PURE_COMPLEX_RATIO',
                               `$1&&$2', `c&&s', `COMPUTE_MIX_COMPLEX_RATIO', 
                               `$1&&$2', `z&&d', `COMPUTE_MIX_COMPLEX_RATIO',
                               `double tmp = fabs((x[ix] - HEAD(r_true)[j]) - TAIL(r_true)[j]);
        
                                if (HEAD(r_true)[j]==0.0)
                                  ratios[j]=0.0;
                                else
                                  ratios[j] = tmp/(HEAD(r_true)[j]*2.0*eps_out);')
                              }
        
                              /* take the max ratio */
                              if (j==0) {
                                ratio = ratios[0];
                              } else if ( !(ratios[j] <= ratio) ) {
                              /* The !<= test causes NaN error to be detected.
                                 Note that (NaN > thresh) is always false. */
                                ratio = ratios[j];
                              }
                              ix += incx;                       
                            }
        
                            /* Increment the bad ratio count, if the ratio
                               is bigger than the threshold.
                               The !<= below causes NaN error to be detected.
                               Note that (NaN > thresh) is always false. */
                            if ( !(ratio <= thresh) ) {
                
                              bad_ratios++;
                             
                              if ((debug==3) &&          /* print only when debug is on */
                                  (1 || count != old_count) )
                              {
                                p_count++;
                                old_count = count;      
                                printf("\ntest %d\n", test_count);      

                                printf("FAIL: %s: n = %d, ntests = %d, threshold = %4.2f,\n", 
                                        fname, n, ntests, thresh);
                
                                /* Print test info */
                                PRINT_PREC(prec)
                                PRINT_NORM(norm)
                                PRINT_ORDER(order_type)
                                PRINT_UPLO(uplo_type)
                                PRINT_TRANS(trans_type)
                                PRINT_DIAG(diag_type)
                                
                                printf("row=%d, ldt=%d, incx=%d:\n", row, ldt, incx_unadj);
                                printf("k = %d\n", k);                          
                                ix=0;
                                if (incx < 0) ix = -(n-1)*incx;
                                
                                printf("      T=");   
                                $2print_tbsv_matrix(T, n, k, ldt, 
                                        order_type,
                                        uplo_type, 
                                        trans_type);  
                                $2print_vector(T, n*ldt, 1, "T");
        
                                ix=0;
                                if (incx < 0) ix = -(n-1)*incx;
                                $1print_vector(x_gen, n, incx_gen_unadj, "x");
                                for (j=0; j<n; j++){
                                  printf("      "); 
                                  PRINT_ARRAY_ELEM(x_gen, j*incx_gen, $1_type, x) printf("; ");
                                  PRINT_ARRAY_ELEM(x, ix, $1_type, x_final) printf("\n");
                                  ix += incx;
                                }
        
                                printf("      "); PRINT_VAR(alpha, $1_type) printf("; "); printf("\n");
                                for (j=0; j<n; j++){
                                  if (j==row) printf("    =>");
                                  else        printf("      ");
                                  PRINT_ARRAY_ELEM(r_true, j*incx_gen, EXTRA_TYPE($1_type))
                                  printf(", ratio[%d]=%.4e\n", j, ratios[j]); 
                                }
                                printf("      ratio=%.4e\n", ratio);
                              }
                              if (bad_ratios >= MAX_BAD_TESTS) {
                                printf("\ntoo many failures, exiting....");
                                printf("\nTesting and compilation");
                                printf(" are incomplete\n\n");
                                goto end;
                              }
                              if (0 && !(ratio <= TOTAL_FAILURE_THRESHOLD) ) {
                                printf("\nFlagrant ratio %e, exiting...", ratio);
                                printf("\nTesting and compilation");
                                printf(" are incomplete\n\n");
                                goto end;
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
                        } /* row */
                    } /* k */
                  } /* ldt */
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

end:
  FPU_FIX_STOP;

  FREE_VECTOR(x, $1_type)
  FREE_VECTOR(x_gen, $1_type)
  FREE_VECTOR(temp, $2_type)
  FREE_VECTOR(T, $2_type)
  FREE_VECTOR(r_true, EXTRA_TYPE($1_type))
  FREE_VECTOR(ratios, real_D)

  *min_ratio = ratio_min;
  *num_bad_ratio = bad_ratios;
  *num_tests = tot_tests;
  return ratio_max;
}' )dnl
dnl
dnl
dnl
define(`COMPUTE_PURE_COMPLEX_RATIO',
`                       double tmp = sqrt(((x[ix] - HEAD(r_true)[j*incx_gen]) - TAIL(r_true)[j*incx_gen])*
                                          ((x[ix] - HEAD(r_true)[j*incx_gen]) - TAIL(r_true)[j*incx_gen])+
                                          ((x[ix+1] - HEAD(r_true)[j*incx_gen+1]) - TAIL(r_true)[j*incx_gen+1])*
                                          ((x[ix+1] - HEAD(r_true)[j*incx_gen+1]) - TAIL(r_true)[j*incx_gen+1]));

                        if (HEAD(r_true)[j*incx_gen] == 0.0 && HEAD(r_true)[j*incx_gen+1] == 0.0)
                                ratios[j] = 0.0;
                        else    
                                ratios[j] = tmp/(sqrt(HEAD(r_true)[j*incx_gen]*HEAD(r_true)[j*incx_gen]+
                                              HEAD(r_true)[j*incx_gen+1]*HEAD(r_true)[j*incx_gen+1])* 
                                          8.0*sqrt(2.0)*eps_out);')
dnl
dnl
dnl
define(`COMPUTE_MIX_COMPLEX_RATIO',
`                       double tmp = sqrt(((x[ix] - HEAD(r_true)[j*incx_gen]) - TAIL(r_true)[j*incx_gen])*
                                          ((x[ix] - HEAD(r_true)[j*incx_gen]) - TAIL(r_true)[j*incx_gen])+
                                          ((x[ix+1] - HEAD(r_true)[j*incx_gen+1]) - TAIL(r_true)[j*incx_gen+1])*
                                          ((x[ix+1] - HEAD(r_true)[j*incx_gen+1]) - TAIL(r_true)[j*incx_gen+1]));
                        if (HEAD(r_true)[j*incx_gen] == 0.0 && HEAD(r_true)[j*incx_gen+1] == 0.0)
                                ratios[j] = 0.0;
                        else    
                                ratios[j] = tmp/(sqrt(HEAD(r_true)[j*incx_gen]*HEAD(r_true)[j*incx_gen]+
                                              HEAD(r_true)[j*incx_gen+1]*HEAD(r_true)[j*incx_gen+1])* 
                                          2.0*sqrt(2.0)*eps_out);')
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: TESTGEN_TBSV_NAME(ax_typeltr, T_typeltr) 
dnl        create a testgen_tbsv name 
dnl --------------------------------------------------------------------
define(`TESTGEN_TBSV_NAME', `ifelse(
        `$1&&$2', `$1&&$1', `BLAS_$1tbsv_testgen',
                            `BLAS_$1tbsv_$2_testgen')')dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: DO_TEST_TBSV_NAME(abx_typeltr, T_typeltr, extended)
dnl        create do_test_dot name
dnl -------------------------------------------------------------------
define(`DO_TEST_TBSV_NAME', `ifelse(
        `$1&&$2', `$1&&$1',
        `do_test_$1tbsv$3',
        `do_test_$1tbsv_$2$3')') dnl
dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: CALL_DO_TEST_TBSV(ax_typeltr, T_typeltr, extended)
dnl
dnl        ax_type  : the type and precision of alpha and x
dnl        T_type   : the type and precision of T
dnl        extended : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s    ... real and single
dnl        d    ... real and double
dnl        c    ... complex and single
dnl        z    ... complex and double
dnl -------------------------------------------------------------------
define(`CALL_DO_TEST_TBSV',           
        `  fname = "TBSV_NAME($1, $2, $3, $4)";
           printf("Testing %s...\n", fname);
           min_ratio = 1e308; max_ratio = 0.0;
           total_bad_ratios = 0; total_tests = 0;
           for(i=0; i<nsizes; i++){
            
              n = sizes[i];
              total_max_ratio = DO_TEST_TBSV_NAME($1, $2, $3)(n, ntests, &seed, thresh, debug, 
                   test_prob, &total_min_ratio, &num_bad_ratio, &num_tests);
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
dnl
dnl
FOREACH(`TBSV_ARGS', `
DO_TEST_TBSV(arg)')dnl

#define NSIZES 12
MAIN(`
  int i;
  int sizes[NSIZES]={0, 1, 2, 3, 4, 5, 6, 10, 15};', `

FOREACH(`TBSV_ARGS', `
CALL_DO_TEST_TBSV(arg)')')dnl

