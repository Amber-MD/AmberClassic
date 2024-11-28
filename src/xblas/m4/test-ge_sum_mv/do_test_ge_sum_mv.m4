dnl Generates test code for ge_sum_mv
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

/* 0 -- 1 */
#define ORDER_START  0
#define ORDER_END    1

/* 0 -- 2 */
#define ALPHA_START  0
#define ALPHA_END    2

/* 0 -- 2 */
#define BETA_START   0
#define BETA_END     2

/* -1 -- 1 */
#define NORM_START   -1
#define NORM_END     1

/* 0 -- 2 */
#define LDA_START    0
#define LDA_END      2

/* 0 -- 2 */
#define PREC_START   0
#define PREC_END     2

/* 0 -- 1 */
#define RANDOMIZE_START 0
#define RANDOMIZE_END   1

/* -2 -- 2 (Stride) */
#define INCX_START -2
#define INCX_END 2

/* -2 -- 2 (Stride) */
#define INCY_START -2
#define INCY_END 2

#define NUM_DATA 7

include(cblas.m4)dnl
include(test-common.m4)dnl
include(ge_sum_mv/ge_sum_mv-common.m4)dnl
dnl
dnl
define(`DO_TEST_GE_SUM_MV_NAME', 
  `ifelse(`$1&&$1', `$2&&$3', `do_test_$1ge_sum_mv$4', `do_test_$1ge_sum_mv_$2_$3$4')')
dnl
dnl
define(`DO_TEST_GE_SUM_MV_PARAMS', 
  `int m, int n,  
   int ntests, int *seed, double thresh, int debug, float test_prob, 
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests')
dnl
dnl
define(`DO_TEST_GE_SUM_MV_HEAD', 
  `void DO_TEST_GE_SUM_MV_NAME($1, $2, $3, $4) 
           (DO_TEST_GE_SUM_MV_PARAMS($1, $2, $3, $4))')
dnl
dnl
define(`DO_TEST_GE_SUM_MV', 
  `DO_TEST_GE_SUM_MV_HEAD($1, $2, $3, $4) {
   DO_TEST_GE_SUM_MV_BODY($1, $2, $3, $4)
  }')
dnl
dnl
dnl
define(`TESTGEN_GE_SUM_MV_NAME', `ifelse(
       `$2&&$3', `$1&&$1', `BLAS_$1ge_sum_mv_testgen', `BLAS_$1ge_sum_mv_$2_$3_testgen')')
dnl
dnl
define(`DO_TEST_GE_SUM_MV_BODY', 
  `
  /* Function name */
  const char fname[] = "GE_SUM_MV_NAME($1, $2, $3, $4)";

  int i;
  int yi;
  int incyi, y_starti, incx_veci;
  int test_count;
  int bad_ratio_count;

  int ri;
  int incri;
  int inca, incx, incy;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;  /* internal machine epsilon     */
  double un_int;   /* internal underflow threshold */

  DECLARE(rin, $1_type)
  DECLARE(rout, $1_type)
  DECLARE(r_true_elem, EXTRA_TYPE($1_type))

  enum blas_order_type order_type;
  enum blas_prec_type prec;
 
  int order_val;
  int lda_val, incx_val, incy_val;
  int ldb_val;
  int alpha_val, beta_val;
  int randomize_val;
  
  ifelse(`$4', `_x', `int prec_val;', `')
 
  int lda, ldb;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  int n_i, m_i;
  int inca_veci;

  DECLARE(alpha, $1_type)
  DECLARE(beta, $1_type)
  DECLARE(beta_zero_fake, $1_type)
  DECLARE(alpha_use, $1_type)
  DECLARE_VECTOR(a, $2_type)
  DECLARE_VECTOR(a_use, $2_type)
  DECLARE_VECTOR(B, $2_type)
  DECLARE_VECTOR(B_use, $2_type)
  DECLARE_VECTOR(x, $3_type)
  DECLARE_VECTOR(y, $1_type)
  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(x_vec, $3_type)


  DECLARE_VECTOR(ratios, real_D)

  /* true result calculated by testgen, in double-double */
  DECLARE_VECTOR(r_true, EXTRA_TYPE($1_type))

  FPU_FIX_DECL;
  
  ZERO(beta_zero_fake, $1_type) 

  if (n < 0 || ntests < 0)
    BLAS_error(fname,  -3,  n, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0)
    return;

  FPU_FIX_START;
 
  n_i = n;
  m_i = m;

  inca = incx = incy = 1;
  INC_ADJUST(inca, $2_type)
  INC_ADJUST(incx, $3_type)
  INC_ADJUST(incy, $1_type)

  /* allocate memory for arrays */
  MALLOC_VECTOR(y, $1_type, 4 * m_i)
  MALLOC_VECTOR(a, $2_type, 2 * n_i * m_i * m_i * n_i)
  MALLOC_VECTOR(a_use, $2_type, 2*m_i * n_i *m_i * n_i)
  MALLOC_VECTOR(B, $2_type, 2*n_i * n_i * m_i * m_i)
  MALLOC_VECTOR(B_use, $2_type, 2*n_i * n_i * m_i * m_i)
  MALLOC_VECTOR(x, $3_type, 4 * n_i)

  inca_veci = 1;
  INC_ADJUST(inca_veci, $2_type)
  MALLOC_VECTOR(a_vec, $2_type, 2*n_i)
  MALLOC_VECTOR(x_vec, $3_type, 2*n_i)
  MALLOC_VECTOR(r_true, EXTRA_TYPE($1_type), m_i)
  MALLOC_VECTOR(ratios, real_D, m_i)

  test_count = 0;
  bad_ratio_count = 0;

    /* vary alpha */
    for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {
      
      SET_ALPHA($1_type)

      /* vary beta */
      for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
        SET_BETA($1_type)

        ifelse($4, _x, `
        /* varying extra precs */
        for (prec_val = PREC_START; prec_val <= PREC_END; prec_val++) {')
        SET_INTERNAL_PARAMS($1_type, $4)

        /* vary norm -- underflow, approx 1, overflow */
        for (norm = NORM_START; norm <= NORM_END; norm++) {
         
          /* number of tests */
          for (test_no = 0; test_no < ntests; test_no++) {


            /* vary storage format */
            for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

              order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

                /* vary lda = n_i, n_i+1, 2*n_i */
                for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {
                  
                        if(order_type == blas_rowmajor) {
                            lda = (lda_val == 0) ? n_i : 
                              (lda_val == 1) ? n_i+1 : n_i*n_i;
                        } else {
                            lda = (lda_val == 0) ? m_i : 
                              (lda_val == 1) ? m_i+1 : m_i*m_i;
                        }

                  /* vary ldb = n_i, n_i+1, 2*n_i */
                  for (ldb_val = LDA_START; ldb_val <= LDA_END; ldb_val++) {
                  
                        if(order_type == blas_rowmajor) {
                            ldb = (ldb_val == 0) ? n_i : 
                              (ldb_val == 1) ? n_i+1 : n_i*n_i;
                        } else {      
                            ldb = (ldb_val == 0) ? m_i : 
                              (ldb_val == 1) ? m_i+1 : m_i*m_i;
                        }

                      for (randomize_val = RANDOMIZE_START; 
                           randomize_val <= RANDOMIZE_END; randomize_val++) {

                        /* For the sake of speed, we throw out this case at random */
                        if ( xrand(seed) >= test_prob ) continue;

                        /* finally we are here to generate the test case */
                        /* alpha_use, a_use, B_use are the generated alpha, a, B
                        *  before any scaling.  
                        *  That is, in the generator, alpha == beta == alpha_use 
                        *  before scaling. */

                        saved_seed = *seed;
                        TESTGEN_GE_SUM_MV_NAME($1, $2, $3) (norm, order_type, 
                              m, n, randomize_val, &alpha, 
                              alpha_flag, &beta, beta_flag, a, lda, B, ldb, x_vec, 1, 
                              &alpha_use, a_use, B_use,
                              seed, HEAD(r_true), TAIL(r_true));

                      /* vary incx = 1, 2*/
                      for (incx_val = INCX_START; incx_val <= INCX_END; incx_val++) {
                      
                        incx = incx_val;
                        if (0 == incx)
                          continue;     
                      
                        $3copy_vector(x_vec, n_i, 1, x, incx);

                        /* vary incy = 1, 2 */
                        for (incy_val = INCY_START; incy_val <= INCY_END; incy_val++) {
                        
                          incy = incy_val;
                          if (0 == incy)
                                continue;
                                
                            test_count++;
        
                            /* call ge_sum_mv routines to be tested */ 
                            FPU_FIX_STOP;
                            GE_SUM_MV_NAME($1, $2, $3, $4)(order_type,
                              m, n, alpha, a, lda, x, incx, beta, B, ldb, 
                              y, incy ifelse(`$4', `_x', `, prec'));
                            FPU_FIX_START;

                            /* now compute the ratio using test_BLAS_xdot */
                            /* copy a row from A, use x, run 
                               dot test */

                            incyi = incy;

                            incri = 1;
                            incx_veci = 1;
                            INC_ADJUST(incx_veci, $3_type)
                            INC_ADJUST(incyi, $1_type)
                            INC_ADJUST(incri, EXTRA_TYPE($1_type))
                            if (incy < 0) {
                                y_starti = (-m_i + 1) * incyi;
                            } else {
                                y_starti = 0;
                            }
                                /* make two copies of x into x_vec. redundant*/
                            $3copy_vector(x, n_i, incx, x_vec, 1);
                            $3copy_vector(x, n_i, incx, (x_vec + (n_i*incx_veci)), 1); 
                            for (i = 0, yi = y_starti, ri = 0; 
                                    i < m_i; i++, yi += incyi, ri += incri) {
                              $2ge_copy_row(order_type, blas_no_trans, dnl
                                  m_i, n_i, a_use, lda, a_vec, i);
                              $2ge_copy_row(order_type, blas_no_trans, dnl
                                  m_i, n_i, B_use, ldb, (a_vec + inca_veci*n_i), i);    
                                
                              ZERO(rin, $1_type)
                              GET_VECTOR_ELEMENT(rout, y, yi, $1_type)
                              GET_VECTOR_ELEMENT(r_true_elem, r_true, 
                                  ri, EXTRA_TYPE($1_type))

                              TEST_DOT_NAME($1, $2, $3, $4)(2*n_i, 
                                blas_no_conj, 
                                alpha_use, beta_zero_fake, rin, rout, 
                                HEAD(r_true_elem), TAIL(r_true_elem), 
                                a_vec, 1, x_vec, 1, eps_int, un_int, 
                                &ratios[i]);

                                /* take the max ratio */
                                if (i == 0) {
                                  ratio = ratios[0];
                                /* The !<= below causes NaN errors
                                 *  to be included.
                                 * Note that (NaN > 0) is false */
                                } else if (!(ratios[i] <= ratio)) {
                                  ratio = ratios[i];
                                }
                            }  /* end of dot-test loop */

                            /* The !<= below causes NaN errors
                             *  to be included.
                             * Note that (NaN > 0) is false */
                            if (  !(ratio <= thresh) ) {

                              if (debug == 3) {
                                printf("\n\t\tTest # %d\n", test_count);
                                printf("y type : $1, a type : $2, x type : $3\n");
                                printf("Seed = %d\t", saved_seed);
                                printf("n %d, m %d\n", 
                                        n, m);
                                printf("LDA %d  LDB %d, INCX %d  INCY %d\n", lda, ldb, incx, incx);

                                if (order_type == blas_rowmajor)
                                  printf("row ");
                                else
                                  printf("col ");

                                printf("NORM %d, ALPHA %d, BETA %d\n", 
                                       norm, alpha_val, beta_val);
                                printf("randomize %d\n",
                                        randomize_val);

                                /* print out info */
                                PRINT_VAR(alpha, $1_type);
                                printf("   ");
                                PRINT_VAR(beta, $1_type);
                                printf("\n");
                                PRINT_VAR(alpha_use, $1_type);
                                printf("\n");

                                $2ge_print_matrix(a, m_i, n_i, lda, order_type, "A");
                                $2ge_print_matrix(B, m_i, n_i, ldb, order_type, "B");
                                $3print_vector(x, n_i, incx, "x");

                                $1print_vector(y, m_i, incy, "y");

                                DOUBLE_ABBREV($1)print_vector(HEAD(r_true), m_i, 1, "HEAD(r_true)");

                                $2ge_print_matrix(a_use, m_i, n_i, lda, order_type, "A_use");
                                $2ge_print_matrix(B_use, m_i, n_i, ldb, order_type, "B_use");

                                dprint_vector(ratios, m_i, 1, "ratios");
                                printf("ratio = %g\n", ratio);
                                fflush(stdout);
                              }
                              bad_ratio_count++;
                              if (bad_ratio_count >= MAX_BAD_TESTS) {
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
                              }                                             }

                            if (!(ratio <= ratio_max))
                              ratio_max = ratio;

                            if (ratio != 0.0 && !(ratio >= ratio_min))
                              ratio_min = ratio;

                        } /* end of incy loop */

                     } /* end of incx loop */

                  } /* end of randmize loop */

                } /* end of ldb loop */     
 
              } /* end of lda loop */     
                
            } /* end of order loop */

          } /* end of nr test loop */
          
        } /* end of norm loop */

dnl ---------------------------------------------------------
ifelse(`$4', `_x', `
          } /* end of prec loop */  ')
dnl ---------------------------------------------------------

      } /* end of beta loop */

    } /* end of alpha loop */

  FPU_FIX_STOP;

end:
  FREE_VECTOR(y, $1_type)
  FREE_VECTOR(a, $2_type)
  FREE_VECTOR(a_use, $2_type)
  FREE_VECTOR(B, $2_type)
  FREE_VECTOR(B_use, $2_type)
  FREE_VECTOR(x, $3_type)
  FREE_VECTOR(r_true, EXTRA_TYPE($1_type))
  FREE_VECTOR(ratios, real_D)
  FREE_VECTOR(a_vec, $2_type)
  FREE_VECTOR(x_vec, $3_type)

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;
') 
dnl
dnl
dnl
define(`CALL_DO_TEST_GE_SUM_MV', 
 `fname = "GE_SUM_MV_NAME($1, $2, $3, $4)";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = n_data[i][0];
    n = n_data[i][1];

    DO_TEST_GE_SUM_MV_NAME($1, $2, $3, $4)(m, n, 
             ntests, &seed, thresh, debug,
             test_prob,   
             &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
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
') 

FOREACH(`GE_SUM_MV_ARGS', `
DO_TEST_GE_SUM_MV(arg)')dnl

MAIN(`
  int m, i;
  int n_data[NUM_DATA][2] = 
      {{1, 1}, {1, 2}, {3, 2}, {8, 6}, {9, 10}, {4, 4}, {7, 7}};', `

FOREACH(`GE_SUM_MV_ARGS', `
CALL_DO_TEST_GE_SUM_MV(arg)')')dnl
