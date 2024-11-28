dnl Generates test code for trmv
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"
include(cblas.m4)dnl
include(test-common.m4)dnl
include(trmv/trmv-common.m4)dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_TRMV_COMMENT(extended)
dnl        if extended, then print info about prec loop in the code 
dnl        structure
dnl -----------------------------------------------------------------------
define(`DO_TEST_TRMV_COMMENT',`
/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on TRMV.
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
 *                    lda loop      -- varying lda: n, n+1, and 2n
 *                        incx loop     -- varying incx: -2, -1, 1, 2
 */')dnl
dnl
dnl
define(`DO_TEST_TRMV_NAME', `ifelse(
  `$1', `$2', `do_test_$1trmv$3', `do_test_$1trmv_$2$3')')dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_TRMV(ax_typeltr, T_typeltr, extended)
dnl
dnl        ax_type  : the type and precision of alpha and x
dnl        T_type  : the type and precision of the matrix T
dnl        extended   : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl     s   ... real and single
dnl     d   ... real and double
dnl     c   ... complex and single
dnl     z   ... complex and double
dnl ----------------------------------------------------------------------
define(`DO_TEST_TRMV', `
double DO_TEST_TRMV_NAME($1, $2, $3)(int n, int ntests, int *seed, 
  double thresh, int debug, float test_prob, double *min_ratio, 
  int *num_bad_ratio, int *num_tests)
  DO_TEST_TRMV_COMMENT($3)
  DO_TEST_TRMV_BODY($1, $2, $3) /* end of DO_TEST_TRMV_NAME($1, $2, $3) */
')dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: DO_TEST_TRMV_BODY(a_typeltr, T_typeltr, x_typeltr, extended)
dnl
dnl      ax_typeltr : the type and precision of alpha, and x
dnl      T_typeltr  : the type and precision of the matrix T
dnl      extended   : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl     s   ... real and single
dnl     d   ... real and double
dnl     c   ... complex and single
dnl     z   ... complex and double
dnl ---------------------------------------------------------------------
define(`DO_TEST_TRMV_BODY',
`{                                        
  /* function name */
  const char fname[] = "TRMV_NAME($1, $2, $3)";

  /* Variables in the "x_val" form are loop vars for corresponding
     variables */
  int i;             /* iterate through the repeating tests */
  int j;             /* multipurpose counters */
  int ix;        /* use to index x and y respectively */
  int lda_val, lda;  /* for testing different values for lda */
  int incx_val, incx;/* for testing different inc values */
  int d_count;       /* counter for debug */
  int find_max_ratio;/* find_max_ratio = 1 only if debug = 3 */
  int p_count;       /* counter for the number of debug lines printed*/
  int tot_tests;     /* total number of tests to be done */
  int norm;          /* input values of near underflow/one/overflow */
  double ratio_max;  /* the current maximum ratio */
  double ratio_min;  /* the current minimum ratio */
  double *ratios;    /* a temporary variable for calculating ratio */
  double ratio;      /* the per-use test ratio from test() */
  int bad_ratios;    /* the number of ratios over the threshold */
  double eps_int;    /* the internal epsilon expected--2^(-24) for float */
  double un_int;     /* the internal underflow threshold */
  DECLARE(alpha, $1_type)
  DECLARE(beta, $1_type)
  DECLARE(rin, $1_type)
  DECLARE_VECTOR(T, $2_type)
  DECLARE_VECTOR(x, $1_type)

  /* used to store vectors generated by testgen; eventually copied to x */
  DECLARE_VECTOR(x_gen, $1_type)

  DECLARE_VECTOR(temp, $2_type) /* use for calculating ratio */
  
  /* the true r calculated by testgen(), in double-double */
  DECLARE_VECTOR(r_true, EXTRA_TYPE($1_type))

  int alpha_val;
  int alpha_flag, beta_flag; /* input flag for TESTGEN_TRMV_NAME($1, $2, $3) */
  int order_val;
  enum blas_order_type order_type;
  ifelse(`$3', `_x', `int prec_val;')
  enum blas_prec_type prec;
  int uplo_val;
  enum blas_uplo_type uplo_type;
  int trans_val;
  enum blas_trans_type trans_type;
  int diag_val;
  enum blas_diag_type diag_type;
  int row; 
  int inc_xgen = 1;
  int ixgen;

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

  bad_ratios = 0;
  incx = 1;
  INC_ADJUST(incx, $1_type)
  INC_ADJUST(inc_xgen, $1_type)

  /* get space for calculation */
  MALLOC_VECTOR(x, $1_type, n*2)
  MALLOC_VECTOR(x_gen, $1_type, n)
  MALLOC_VECTOR(T, $2_type, 2*n*n)
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
  ZERO(rin, $1_type)
  beta_flag=1;

  if (debug == 3)
    find_max_ratio = 1;

  /* The debug iteration:
     If debug=1, then will execute the iteration twice. First, compute the
     max ratio. Second, print info if ratio > (50% * ratio_max). */
  for (d_count = 0; d_count <= find_max_ratio; d_count++) {
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
          default:
          case 1: order_type = blas_colmajor; break;
        }

          /* upper or lower */
          for (uplo_val=0; uplo_val<2; uplo_val++) {
            switch(uplo_val){
              case 0: uplo_type = blas_upper; break;
              default:
              case 1: uplo_type = blas_lower; break;
            }

            /* no trans, trans, or conj trans */
            for (trans_val=0; trans_val<3; trans_val++) {
              switch(trans_val){
                case 0: trans_type = blas_no_trans; break;
                case 1: trans_type = blas_trans; break;
                default:
                case 2: trans_type = blas_conj_trans; break;
              }

              /* non_unit_diag, unit_diag */
              for (diag_val=0; diag_val<2; diag_val++) {
                switch(diag_val){
                  case 0: diag_type = blas_non_unit_diag; break;
                  default:
                  case 1: diag_type = blas_unit_diag; break;
                }        

                /* number of tests */
                for (i=0; i<ntests; i++) {    

                  for (lda_val=0; lda_val<3; lda_val++) {
                    switch(lda_val) {
                      case 0: lda = n; break;
                      case 1: lda = n+1; break;
                      default:
                      case 2: lda = 2*n; break;
                    }

                    /* For the sake of speed, we throw out this case at random */
                    if ( xrand(seed) >= test_prob ) continue;

                    /* generate test case */
                    TESTGEN_TRMV_NAME($1, $2)(norm, order_type, uplo_type, 
                      trans_type, diag_type, n, PASS_BY_REF(alpha, $1_type), 
                      alpha_flag, T, lda, x_gen, seed, HEAD(r_true), TAIL(r_true));

                    count++;

                    /* varying incx */
                    for (incx_val = -2; incx_val <= 2; incx_val++) {
                      if (incx_val == 0) continue;
        
                      /* setting incx */
                      incx = incx_val;
                      INC_ADJUST(incx, $1_type)      

                      /* set x starting index */
                      ix=0;
                      if (incx < 0) ix = -(n-1)*incx;
    
                      /* copy x_gen to x */
                      ixgen = 0;
                      for(j=0 ; j<n; j++) {
                        COPY_VECTOR_ELEMENT(x, ix, x_gen, ixgen, $1_type)
                        ix += incx;
                        ixgen += inc_xgen;
                      } 

                      /* call TRMV_NAME($1, $2, $3) */
                      FPU_FIX_STOP;
                      TRMV_NAME($1, $2, $3)(order_type, uplo_type, trans_type, 
                          diag_type, n, alpha, T, lda, x, incx_val 
                          ifelse(`$3', `_x', `, prec'));
                      FPU_FIX_START;

                      /* set x starting index */
                      ix = 0;
                      if (incx < 0) ix = -(n-1)*incx;

                      /* computing the ratio */
                      ratio = 0.0;
                      row = 0;
                      for(j=0; j<n; j++) {
                        /* copy row j of T to temp */
                        $2tr_copy_row(order_type, uplo_type, trans_type, n, T, 
                            lda, temp, j);

                        TEST_DOT_NAME($1, $2, $1, $3)(n, blas_no_conj,
                            alpha, beta, rin, 
                            VECTOR_ELEMENT(x, ix, $1_type), 
                            VECTOR_ELEMENT(r_true, j*inc_xgen, EXTRA_TYPE($1_type)),
                            temp, 1, x_gen, 1, eps_int, un_int, &ratios[j]);

                        /* take the max ratio */
                        if (j==0)
                          ratio = ratios[0];
                        else if (ratios[j] > ratio) {
                          ratio = ratios[j];
                          row = j;
                        }
          
                          ix += incx;      
                      }

                      /* increase the number of bad ratio, if the ratio is 
                         bigger than the threshold */
                      if (ratio > thresh) {
                        bad_ratios++;
                     
                        if (debug==3) {
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

                          printf("row=%d, lda=%d, incx=%d:\n", row, lda, incx);
      
                          ix=0;
                          if (incx < 0) ix = -(n-1)*incx;
        
                          printf("      T=");     
                          for (j=0; j<n; j++) {
                            /* copy row j of T to temp */
                            $2tr_copy_row(order_type, uplo_type, trans_type, n, 
                                T, lda, temp, j);
       
                            if (j>0) printf("        ");
                            $2print_vector(temp, n, 1, NULL);
                          }
 
                          ix=0;
                          if (incx < 0) ix = -(n-1)*incx;
         
                          for (j=0; j<n; j++){
                            printf("      "); 
                            PRINT_ARRAY_ELEM(x_gen, j, $1_type, x) printf("; ");
                            PRINT_ARRAY_ELEM(x, ix, $1_type, x_final) printf("\n");
                            ix += incx;
                          }

                          printf("      "); 
                          PRINT_VAR(alpha, $1_type) printf("; "); printf("\n");
                          for (j=0; j<n; j++) {
                            if (j==row) printf("    =>");
                            else        printf("      ");
                            PRINT_ARRAY_ELEM(r_true, j, EXTRA_TYPE($1_type))
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
                  } /* lda */
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
dnl --------------------------------------------------------------------
dnl Usage: TESTGEN_TRMV_NAME(a_typeltr, T_typeltr, x_typeltr) 
dnl        create a testgen_trmv name 
dnl --------------------------------------------------------------------
define(`TESTGEN_TRMV_NAME',
  `ifelse(`$1', `$2', `BLAS_$1trmv_testgen', `BLAS_$1trmv_$2_testgen')')dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: DO_TEST_TRMV_NAME(a_typeltr, T_typeltr, x_typeltr, extended)
dnl        create do_test_dot name
dnl -------------------------------------------------------------------
define(`DO_TEST_TRMV_NAME', `ifelse(
        `$1', `$2',
        `do_test_$1trmv$3',
        `do_test_$1trmv_$2$3')') dnl
dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: CALL_DO_TEST_TRMV(a_typeltr, T_typeltr, x_typeltr, extended)
dnl
dnl        a_type   : the type and precision of alpha
dnl      T_type   : the type and precision of T
dnl        x_type   : the type and precision of x
dnl        extended : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s   ... real and single
dnl     d   ... real and double
dnl     c  ... complex and single
dnl     z  ... complex and double
dnl -------------------------------------------------------------------
define(`CALL_DO_TEST_TRMV',           
  `  fname = "TRMV_NAME($1, $2, $3)";
     printf("Testing %s...\n", fname);
     min_ratio = 1e308; max_ratio = 0.0;
     total_bad_ratios = 0; total_tests = 0;
     for(i=0; i<nsizes; i++){
      
     n = sizes[i];
     total_max_ratio = DO_TEST_TRMV_NAME($1, $2, $3)(n, ntests, &seed, 
     thresh, debug, test_prob, &total_min_ratio, &num_bad_ratio, &num_tests);
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
     else {
       nr_failed_routines++;
       printf("FAIL> ");
     }

     if (min_ratio == 1e308) min_ratio = 0.0;

     printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n",
         fname, total_bad_ratios, total_tests, max_ratio);
')dnl
dnl
dnl
FOREACH(`TRMV_ARGS', `
DO_TEST_TRMV(arg)')dnl

#define NSIZES  12
MAIN(`
  int i;
  int sizes[NSIZES]={0, 1, 2, 3, 4, 5, 6, 10, 17, 28, 37, 58};', `

FOREACH(`TRMV_ARGS', `
CALL_DO_TEST_TRMV(arg)')')dnl

