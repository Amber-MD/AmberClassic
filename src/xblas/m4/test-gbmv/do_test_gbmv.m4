dnl Generates test code for gbmv
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"
include(cblas.m4)dnl
include(test-common.m4)dnl
include(gbmv/gbmv-common.m4)dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_GBMV_COMMENT(extended)
dnl        if extended, then print info about prec loop in the code 
dnl        structure
dnl -----------------------------------------------------------------------
define(`DO_TEST_GBMV_COMMENT',`
/*
 * Purpose  
 * =======
 *
 * Runs a series of tests on GBMV.
 *
 * Arguments
 * =========
 *  
 * m         (input) int
 *           The number of rows
 *
 * n         (input) int
 *           The number of columns
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
 *      beta loop   -- varying beta: 0, 1, or random
ifelse(`$1', `_x',`   *        prec loop   -- varying internal prec: single, double, or extra', `')
 *          norm loop   -- varying norm: near undeflow, near one, or 
 *                        -- near overflow
 *            numtest loop  -- how many times the test is perform with 
 *                            -- above set of attributes
 *              order loop   -- varying order type: rowmajor or colmajor
 *                trans loop    -- varying trans type: no trans, trans, or conj trans
 *                  ku loop       -- varying ku: 0 to n-1
 *                    kl loop       -- varying kl: 0 to m-1
 *                      lda loop      -- varying lda: ku+kl+1, ku+kl+2, 2*(ku+kl+1) 
 *                        incx loop     -- varying incx: -2, -1, 1, 2
 *                          incy loop     -- varying incy: -2, -1, 1, 2
 */')dnl
dnl
dnl
dnl
dnl -----------------------------------------------------------------------
dnl Usage: DO_TEST_GBMV(abr_typeltr, x_typeltr, y_typeltr, extended)
dnl
dnl        aby_typeltr : the type and precision of alpha, beta and r
dnl      AB_typeltr   : the type and precision of x
dnl      x_typeltr   : the type and precision of y
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s   ... real and single
dnl     d   ... real and double
dnl     c  ... complex and single
dnl     z  ... complex and double
dnl ----------------------------------------------------------------------
define(`DO_TEST_GBMV', 
  `double DO_TEST_GBMV_NAME($1, $2, $3, $4)(int m, int n, int ntests, dnl
       int *seed, double thresh, int debug, float test_prob, dnl
       double *min_ratio, int *num_bad_ratio, int *num_tests)
DO_TEST_GBMV_COMMENT($4)
DO_TEST_GBMV_BODY($1, $2, $3, $4)')dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: DO_TEST_GBMV_BODY(aby_typeltr, AB_typeltr, x_typeltr, extended)
dnl
dnl        abr_typeltr : the type and precision of alpha, beta
dnl      x_typeltr   : the type and precision of x
dnl      y_typeltr   : the type and precision of y
dnl        extended    : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s   ... real and single
dnl     d   ... real and double
dnl     c  ... complex and single
dnl     z  ... complex and double
dnl ---------------------------------------------------------------------
define(`DO_TEST_GBMV_BODY',
`{                                        
  /* function name */
  const char fname[] = "GBMV_NAME($1, $2, $3, $4)";

  /* max number of debug lines to print */
  const int max_print = 8;

  /* Variables in the "x_val" form are loop vars for corresponding
     variables */
  int i;            /* iterate through the repeating tests */
  int j;
  int k;         /* multipurpose counters or variables */
  int ix, iy;       /* use to index x and y respectively */
  int incx_val, incy_val, /* for testing different inc values */
      incx, incy;   
  int incx_gen, incy_gen; /* for complex case inc=2, for real case inc=1 */
  int d_count;      /* counter for debug */
  int find_max_ratio; /* find_max_ratio = 1 only if debug = 3 */
  int p_count;      /* counter for the number of debug lines printed*/
  int tot_tests;    /* total number of tests to be done */
  int norm;         /* input values of near underflow/one/overflow */
  double ratio_max; /* the current maximum ratio */
  double ratio_min; /* the current minimum ratio */
  double *ratios;   /* a temporary variable for calculating ratio */
  double ratio;     /* the per-use test ratio from test() */
  int bad_ratios = 0;   /* the number of ratios over the threshold */
  double eps_int;   /* the internal epsilon expected--2^(-24) for float */
  double un_int;    /* the internal underflow threshold */
  DECLARE(alpha, $1_type)
  DECLARE(beta, $1_type)
  DECLARE_VECTOR(AB, $2_type)
  DECLARE_VECTOR(x, $3_type)
  DECLARE_VECTOR(y, $1_type)
  DECLARE_VECTOR(temp, $2_type) /* use for calculating ratio */
  
  /* x_gen and y_gen are used to store vectors generated by testgen.
     they eventually are copied back to x and y */
  DECLARE_VECTOR(x_gen, $3_type)
  DECLARE_VECTOR(y_gen, $1_type) 

  /* the true r calculated by testgen(), in double-double */
  DECLARE_VECTOR(r_true, EXTRA_TYPE($1_type))
  int alpha_val;
  int alpha_flag = 0;   /* input flag for TESTGEN_GBMV_NAME($1, $2, $3, $4) */
  int beta_val;
  int beta_flag = 0;    /* input flag for TESTGEN_GBMV_NAME($1, $2, $3, $4) */
  int order_val;
  enum blas_order_type order_type = 0;
  ifelse(`$4', `_x', `int prec_val;')
  enum blas_prec_type prec = 0;
  int trans_val;
  enum blas_trans_type trans_type = 0;
  int m_i = 0;
  int n_i = 0;
  int max_mn; /* the max of m and n */
  int ku;
  int kl;
  int lda_val;
  int lda = 0;
  int saved_seed;   /* for saving the original seed */

  /* use for counting the number of testgen calls * 2 */  
  int count, old_count = -1;  

  FPU_FIX_DECL; 

  /* test for bad arguments */
  if (n < 0 || m < 0 || ntests < 0) 
    BLAS_error(fname, 0, 0, NULL);   

  /* initialization */
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  saved_seed = *seed;
  ratio_min = 1e308;
  ratio_max = 0.0;
  ratio = 0.0;
  tot_tests = 0;
  p_count = 0;
  count = 0; 
  find_max_ratio = 0;
  if (debug == 3)
    find_max_ratio = 1;
  max_mn = MAX(m, n);

  if (m == 0 || n == 0) {
    return 0.0;
  }

  FPU_FIX_START;

  incx_gen = incy_gen = 1;
  INC_ADJUST(incx_gen, $3_type)
  INC_ADJUST(incy_gen, $1_type)

  /* get space for calculation */
  MALLOC_VECTOR(x, $3_type, max_mn*2)
  MALLOC_VECTOR(y, $1_type, max_mn*2)
  MALLOC_VECTOR(x_gen, $3_type, max_mn)
  MALLOC_VECTOR(y_gen, $1_type, max_mn)
  MALLOC_VECTOR(temp, $2_type, max_mn)
  MALLOC_VECTOR(r_true, EXTRA_TYPE($1_type), max_mn)
  MALLOC_VECTOR(ratios, real_D, max_mn)

  MALLOC_VECTOR(AB, $2_type, (m-1+n-1+1)*max_mn*2)

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

        ifelse($4, _x, `
        /* varying extra precs */
        for (prec_val = 0; prec_val <= 2; prec_val++) {')
          SET_INTERNAL_PARAMS($1_type, $4)
           
          /* values near underflow, 1, or overflow */
          for (norm = -1; norm <= 1; norm++) {

            /* number of tests */
            for (i=0; i<ntests; i++){    

              /* row or col major */
              for (order_val=0; order_val<2; order_val++) {
                switch(order_val){
                  case 0: order_type = blas_rowmajor; break;
                  case 1: order_type = blas_colmajor; break;
                }

                /* no_trans, trans, or conj_trans */
                for (trans_val=0; trans_val<3; trans_val++) {
                  switch(trans_val){
                    case 0: trans_type = blas_no_trans; m_i=m; n_i=n; break;
                    case 1: trans_type = blas_trans; m_i=n; n_i=m; break;
                    case 2: trans_type = blas_conj_trans; m_i=n; n_i=m; break;
                  }

                  /* ku from 0 to n-1 */
                  for (ku=0; ku<n; ku++) {  
                    if (ku == n && ku != 0) continue; /* the purpose of doing this is 
                                                         to test for ku=0 */
    
                    /* kl from 0 to m-1 */
                    for (kl=0; kl<m; kl++) {
                      if (kl == n && kl != 0) continue; /* the purpose of doing this is 
                                                           to test for kl=0 */

                      /* lda=ku+kl+1, ku+kl+2, 2*(ku+kl+1) */
                      for (lda_val=0; lda_val<3; lda_val++){
                        switch(lda_val){
                          case 0: lda=ku+kl+1; break;
                          case 1: lda=ku+kl+2; break;
                          case 2: lda=2*(ku+kl+1); break;
                        }
      
                        if ((order_type == blas_rowmajor && lda < n) ||
                           (order_type == blas_colmajor && lda < m))
                          continue;

                          /* For the sake of speed, we throw out this case at random */
                          if ( xrand(seed) >= test_prob ) continue;

                          /* in the trivial cases, no need to run testgen */
                          if (m > 0 && n > 0) 
                            TESTGEN_GBMV_NAME($1, $2, $3)(norm, order_type, dnl
                                trans_type, m, n, kl, ku, &alpha, alpha_flag, dnl
                                AB, lda, x_gen, &beta, beta_flag, y_gen, dnl
                                seed, HEAD(r_true), TAIL(r_true));
                            count++;

                            /* varying incx */
                            for (incx_val = -2; incx_val <= 2; incx_val++){
                              if (incx_val == 0) continue;         

                            /* setting incx */
                            incx = incx_val;
                            INC_ADJUST(incx, $3_type)      

                            $3copy_vector(x_gen, n_i, 1, x, incx_val);

                            /* varying incy */
                            for (incy_val = -2; incy_val <= 2; incy_val++){
                              if (incy_val == 0) continue;

                              /* setting incy */
                              incy = incy_val;
                              INC_ADJUST(incy, $1_type)    

                              $1copy_vector(y_gen, m_i, 1, y, incy_val);
    
                              /* call GBMV_NAME($1, $2, $3, $4) */
                              FPU_FIX_STOP;
                              GBMV_NAME($1, $2, $3, $4)(order_type, dnl
                                  trans_type, m, n, kl, ku, alpha, AB, dnl
                                  lda, x, incx_val, beta, y, incy_val dnl
                                  ifelse(`$4', `_x', `, prec'));
                              FPU_FIX_START;

                              /* set y starting index */
                              iy=0;
                              if (incy < 0) iy = -(m_i-1)*incy;

                              /* computing the ratio */
                              for(j=0; j<m_i; j++){
                                /* copy row j of AB to temp */
                                $2gbmv_copy(order_type, trans_type, m, n, dnl
                                    kl, ku, AB, lda, temp, j);

                                TEST_DOT_NAME($1, $2, $3, $4)(n_i, blas_no_conj, alpha, beta, 
                                    VECTOR_ELEMENT(y_gen, j*incy_gen, $1_type), 
                                    VECTOR_ELEMENT(y, iy, $1_type), 
                                    VECTOR_ELEMENT(r_true, j*incy_gen, EXTRA_TYPE($1_type)),
                                    temp, 1, x, incx_val, eps_int, un_int,
                                    &ratios[j]);

                                /* take the max ratio */
                                if (j==0) {
                                  ratio = ratios[0];
                                  /* The !<= below causes NaN error to be detected.
                                     Note that (NaN > thresh) is always false. */
                                } else if ( !(ratios[j] <= ratio) ) {
                                  ratio = ratios[j];
                                }

                                iy += incy;      
                              }

                              /* Increase the number of bad ratio, if the ratio
                                 is bigger than the threshold.
                                 The !<= below causes NaN error to be detected.
                                 Note that (NaN > thresh) is always false. */
                              if ( !(ratio <= thresh) ) {
                                bad_ratios++;
                  
                                if ((debug==3) &&          /* print only when debug is on */
                                    (count != old_count) && /* print if old vector is different 
                                                             from the current one */ 
                                    (d_count == find_max_ratio) &&  
                                    (p_count <= max_print) &&
                                    (ratio > 0.5*ratio_max)) {
                                  old_count = count;   
                                  printf("FAIL> %s: m = %d, n = %d, ntests = %d, threshold = %4.2f,\n", 
                                      fname, m, n, ntests, thresh);
  
                                  /* Print test info */
                                  PRINT_PREC(prec)
                                  PRINT_NORM(norm)
                                  PRINT_ORDER(order_type)
                                  PRINT_TRANS(trans_type)

                                  printf("ku=%d, kl=%d, lda=%d, incx=%d, incy=%d:\n", ku, kl, lda, incx, incy);
      
                                  ix=0; iy=0;
                                  if (incx < 0) ix = -(n_i-1)*incx;
                                  if (incy < 0) iy = -(m_i-1)*incy;

                                  printf("      A=");     
                                  for (j=0; j<m_i; j++) {
                                    /* copy row j of A to temp */
                                    $2gbmv_copy(order_type, trans_type, m, n, kl, ku, AB, lda, temp, j);
         
                                    if (j>0) printf("        ");
                                    $2print_vector(temp, n_i, 1, NULL);
                                  }

                                  for (j=0, k=0; j<n_i || k<m_i; j++, k++) {
                                    if (j<n_i) {
                                      printf("      "); PRINT_ARRAY_ELEM(x, ix, $3_type) printf("\n");
                                    }
                                    if (k<m_i) {
                                      printf("      "); PRINT_ARRAY_ELEM(y_gen, k*incy_gen, $1_type) printf("\n");
                                      printf("      "); PRINT_ARRAY_ELEM(y, iy, $1_type, y_final) printf("\n");            
                                    }
                                    ix += incx; iy += incy;
                                  }

                                  printf("      "); PRINT_VAR(alpha, $1_type)
                                  printf("\n      "); PRINT_VAR(beta, $1_type) printf("\n");
                                  for (j=0; j<m_i; j++) {
                                    printf("      ");
                                    PRINT_ARRAY_ELEM(r_true, j*incy_gen, EXTRA_TYPE($1_type))
                                    printf(", ratio[%d]=%.4e\n", j, ratios[j]); 
                                  }
 
                                  printf("      ratio=%.4e\n", ratio);
                                  p_count++;
                                }
                                if (bad_ratios >= MAX_BAD_TESTS) {
                                  printf("\ntoo many failures, exiting....");
                                  printf("\nTesting and compilation");
                                  printf(" are incomplete\n\n");
                                  goto end;
                                }
                                if ( !(ratio <= TOTAL_FAILURE_THRESHOLD) ) {
                                  printf("\nFlagrant ratio error, exiting...");
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
                          }/* incy */
                        }/* incx */
                      }/* lda */
                    } /* kl */
                  } /* ku */
                } /* trans */
              } /* order */
            } /* tests */
          } /* norm */
ifelse(`$4', `_x', `         } /* prec */')
      } /* beta */
    } /* alpha */
  } /* debug */

  if ((debug == 2) ||
     ((debug == 1) && bad_ratios > 0)){
    printf("      %s:  m = %d, n = %d, ntests = %d, thresh = %4.2f\n", dnl
        fname, m, n, ntests, thresh);
    printf("      bad/total = %d/%d=%3.2f, min_ratio = %.4e, max_ratio = %.4e\n\n",
        bad_ratios, tot_tests, dnl
        ((double)bad_ratios)/((double)tot_tests), ratio_min, ratio_max);
  }

end:
  FREE_VECTOR(x, $3_type)
  FREE_VECTOR(y, $1_type)
  FREE_VECTOR(x_gen, $3_type)
  FREE_VECTOR(y_gen, $1_type)
  FREE_VECTOR(temp, $2_type)
  FREE_VECTOR(AB, $2_type)
  FREE_VECTOR(r_true, EXTRA_TYPE($1_type))
  FREE_VECTOR(ratios, real_D)

  FPU_FIX_STOP;

  *min_ratio = ratio_min;
  *num_bad_ratio = bad_ratios;
  *num_tests = tot_tests;
  return ratio_max;
}')dnl
dnl
dnl
dnl --------------------------------------------------------------------
dnl Usage: TESTGEN_GBMV_NAME(aby_typeltr, AB_typeltr, x_typeltr) 
dnl        create a testgen_gbmv name 
dnl --------------------------------------------------------------------
define(`TESTGEN_GBMV_NAME', `ifelse(
        `$2&&$3', `$1&&$1', `BLAS_$1gbmv_testgen',
                            `BLAS_$1gbmv_$2_$3_testgen')')dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: DO_TEST_GBMV_NAME(abr_typeltr, x_typeltr, y_typeltr, extended)
dnl        create do_test_gbmv name
dnl -------------------------------------------------------------------
define(`DO_TEST_GBMV_NAME', `ifelse(
        `$2&&$3', `$1&&$1',
        `do_test_$1gbmv$4',
        `do_test_$1gbmv_$2_$3$4')') dnl
dnl
dnl
dnl
dnl -------------------------------------------------------------------
dnl Usage: CALL_DO_TEST_GBMV(abr_typeltr, x_typeltr, y_typeltr, extended)
dnl
dnl        abr_type : the type and precision of alpha, beta and r
dnl      x_type   : the type and precision of x
dnl      y_type   : the type and precision of y
dnl        extended : `' if no extended, or `_x' if extended
dnl Each type and precision specifier can be one of
dnl        s   ... real and single
dnl     d   ... real and double
dnl     c  ... complex and single
dnl     z  ... complex and double
dnl -------------------------------------------------------------------
define(`CALL_DO_TEST_GBMV',           
 `min_ratio = 1e308; max_ratio = 0.0;
  total_bad_ratios = 0; total_tests = 0;
  fname = "GBMV_NAME($1, $2, $3, $4)";
  printf("Testing %s...\n", fname);
  for(i=0; i<nsizes; i++){
    m = mn_pairs[i][0];
    n = mn_pairs[i][1];
    total_max_ratio = DO_TEST_GBMV_NAME($1, $2, $3, $4)(m, n, 1, dnl
        &seed, thresh, debug, test_prob, &total_min_ratio, dnl
        &num_bad_ratio, &num_tests);
   if (total_max_ratio > max_ratio)
      max_ratio = total_max_ratio;

   if (total_min_ratio != 0.0 && total_min_ratio < min_ratio)
           min_ratio = total_min_ratio;

   total_bad_ratios += num_bad_ratio;
   total_tests += num_tests;   
  }

  if (total_bad_ratios == 0)
    printf("PASS> ");
  else {
    nr_failed_routines++;
    printf("FAIL> ");
  }
  nr_routines++;

  if (min_ratio == 1e308) min_ratio = 0.0;

  printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n", dnl
      fname, total_bad_ratios, total_tests, max_ratio);')
dnl
dnl
dnl
FOREACH(`GBMV_ARGS', `
DO_TEST_GBMV(arg)')dnl

#define NUMPAIRS 12
MAIN(`
  int m, i;
  int mn_pairs[NUMPAIRS][2]={{0, 0}, {1, 0}, {0, 1}, {1, 2}, {2, 1}, {1, 3},
      {3, 1}, {2, 3}, {4, 3}, {2, 4}, {6, 6}, {10, 8}};', 

`FOREACH(`GBMV_ARGS', `
CALL_DO_TEST_GBMV(arg)')')dnl

