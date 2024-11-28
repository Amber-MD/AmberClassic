dnl Generates test code for symm
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

/* 0 -- 1 */
#define UPLO_START 0
#define UPLO_END   1

/* 0 -- 1 */
#define SIDE_START 0
#define SIDE_END   1

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
#define LDB_START    0
#define LDB_END      2

/* 0 -- 2 */
#define LDC_START    0
#define LDC_END      2

/* 0 -- 2 */
#define PREC_START   0
#define PREC_END     2

/* 0 -- 1 */
#define RANDOMIZE_START 0
#define RANDOMIZE_END   1

#define NUM_DATA 7

include(cblas.m4)dnl
include(test-common.m4)dnl
include(symm/symm-common.m4)dnl
dnl
dnl
define(`DO_TEST_SYMM_NAME', 
  `ifelse(`$1&&$1', `$2&&$3', `do_test_$1symm$4', `do_test_$1symm_$2_$3$4')')
dnl
dnl
define(`DO_TEST_SYMM_PARAMS', 
  `int m, int n,  
   int ntests, int *seed, double thresh, int debug, float test_prob,
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests')
dnl
dnl
define(`DO_TEST_SYMM_HEAD', 
  `void DO_TEST_SYMM_NAME($1, $2, $3, $4) 
           (DO_TEST_SYMM_PARAMS($1, $2, $3, $4))')
dnl
dnl
define(`DO_TEST_SYMM', 
  `DO_TEST_SYMM_HEAD($1, $2, $3, $4) {
   DO_TEST_SYMM_BODY($1, $2, $3, $4)
  }')
dnl
dnl
dnl
define(`TESTGEN_SYMM_NAME', `ifelse(
       `$2&&$3', `$1&&$1', `BLAS_$1symm_testgen',
       `BLAS_$1symm_$2_$3_testgen')')
dnl
dnl
define(`DO_TEST_SYMM_BODY',`
  /* Function name */
  const char fname[] = "SYMM_NAME($1, $2, $3, $4)";

  int i, j;
  int ci, cij;
  int incci, inccij;
  int test_count;
  int bad_ratio_count;

  int ri, rij;
  int incri, incrij;
  int inca, incb, incc;

  double ratio;

  double ratio_min, ratio_max;

  double eps_int;  /* internal machine epsilon     */
  double un_int;   /* internal underflow threshold */

  DECLARE(rin, $1_type)
  DECLARE(rout, $1_type)
  DECLARE(r_true_elem, EXTRA_TYPE($1_type))

  enum blas_order_type order_type;
  enum blas_side_type side_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;
 
  int order_val, uplo_val, side_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;
  
  ifelse(`$4', `_x', `int prec_val;', `')
 
  int lda, ldb, ldc;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;
  int max_mn;

  int m_i, n_i;

  DECLARE(alpha, $1_type)
  DECLARE(beta, $1_type)
  DECLARE_VECTOR(a, $2_type)
  DECLARE_VECTOR(b, $3_type)
  DECLARE_VECTOR(c, $1_type)
  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(b_vec, $3_type)

  /* generated test values for c */
  DECLARE_VECTOR(c_gen, $1_type)

  DECLARE_VECTOR(ratios, real_D)

  /* true result calculated by testgen, in double-double */
  DECLARE_VECTOR(r_true, EXTRA_TYPE($1_type))

  FPU_FIX_DECL;
  
  if (n < 0 || m < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (m == 0 || n == 0)
    return;

  FPU_FIX_START;
 
  max_mn = (m > n) ? m : n;
 
  inca = incb = incc = 1;
  INC_ADJUST(inca, $2_type)
  INC_ADJUST(incb, $3_type)
  INC_ADJUST(incc, $1_type)

  /* allocate memory for arrays */
  MALLOC_VECTOR(c, $1_type, 2 * m * n)
  MALLOC_VECTOR(c_gen, $1_type, 2 * m * n)
  MALLOC_VECTOR(a, $2_type, 2 * max_mn * max_mn)
  MALLOC_VECTOR(b, $3_type, 2 * m * n)
  MALLOC_VECTOR(a_vec, $2_type, max_mn)
  MALLOC_VECTOR(b_vec, $3_type, max_mn)
  MALLOC_VECTOR(r_true, EXTRA_TYPE($1_type), 2 * m * n)
  MALLOC_VECTOR(ratios, real_D, 2 * m * n)

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

              /* vary upper / lower variation */
              for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

                uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;
                
                /* vary left / right multiplication */
                for (side_val = SIDE_START; side_val <= SIDE_END; side_val++) {

                  side_type = (side_val == 0) ? 
                    blas_left_side : blas_right_side;
                        
                  if (side_type == blas_left_side) {
                    m_i = m; n_i = n;
                  } else {
                    m_i = n; n_i = m;
                  }

                  /* vary lda = m_i, m_i+1, 2*m_i */
                  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {
                  
                    lda = (lda_val == 0) ? m_i : 
                      (lda_val == 1) ? m_i+1 : 2*m_i;
                      
                    /* vary ldb = n_i, n_i+1, 2*n_i */
                    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {
                      
                      if (order_type == blas_colmajor)
                        ldb = (ldb_val == 0) ? m : 
                          (ldb_val == 1) ? m+1 : 2*m;
                      else
                        ldb = (ldb_val == 0) ? n :
                          (ldb_val == 1) ? n+1 : 2*n;
                      
                      /* vary ldc = k, k+1, 2*k */
                      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {
                        
                        if (order_type == blas_colmajor)
                          ldc = (ldc_val == 0) ? m : 
                            (ldc_val == 1) ? m+1 : 2*m;
                        else
                          ldc = (ldc_val == 0) ? n :
                            (ldc_val == 1) ? n+1 : 2*n;

                        for (randomize_val = RANDOMIZE_START; 
                             randomize_val <= RANDOMIZE_END; randomize_val++) {

                          /* For the sake of speed, we throw out this case at random */
                          if ( xrand(seed) >= test_prob ) continue;

                            saved_seed = *seed;

                            /* finally we are here to generate the test case */
                            TESTGEN_SYMM_NAME($1, $2, $3) (norm, order_type, 
                              uplo_type, side_type, m, n, randomize_val, &alpha, 
                              alpha_flag, &beta, beta_flag, a, lda, b, ldb, 
                              c, ldc, seed, HEAD(r_true), TAIL(r_true));
                            test_count++;

                            /* copy generated C matrix since this will be
                               over written */
                            $1ge_copy_matrix(order_type, m, n, c_gen, ldc, c, ldc);

                            /* call symm routines to be tested */ 
                            FPU_FIX_STOP;
                            SYMM_NAME($1, $2, $3, $4)(order_type, side_type, 
                              uplo_type, m, n, alpha, a, lda, b, ldb, beta, 
                              c, ldc ifelse(`$4', `_x', `, prec'));
                            FPU_FIX_START;

                            /* now compute the ratio using test_c_xdot */
                            /* copy a row from A, a column from B, run 
                               dot test */

                            if ((order_type == blas_colmajor && 
                                 side_type == blas_left_side) ||
                                (order_type == blas_rowmajor &&
                                 side_type == blas_right_side)) {
                              incci = 1;
                              inccij = ldc;
                            } else {
                              incci = ldc;
                              inccij = 1;
                            }

                            incri = incci;
                            incrij = inccij;
                            INC_ADJUST(incci, $1_type)
                            INC_ADJUST(inccij, $1_type)

                            for (i = 0, ci = 0, ri = 0; 
                                    i < m_i; i++, ci += incci, ri += incri) {
                              $2sy_copy_row(order_type, uplo_type, m_i, a, lda, a_vec, i);
                              for (j = 0, cij = ci, rij = ri; 
                                    j < n_i; j++, cij += inccij, rij += incrij) {
                                /* copy i-th row of A and j-th col of B */
                                if (side_type == blas_left_side)
                                  $3ge_copy_col(order_type, blas_no_trans, 
                                                  m, n, b, ldb, b_vec, j);
                                else
                                  $3ge_copy_row(order_type, blas_no_trans, 
                                                  m, n, b, ldb, b_vec, j);
                                GET_VECTOR_ELEMENT(rin, c_gen, cij, $1_type)
                                GET_VECTOR_ELEMENT(rout, c, cij, $1_type)
                                GET_VECTOR_ELEMENT(r_true_elem, r_true, cij, EXTRA_TYPE($1_type))

                                TEST_DOT_NAME($1, $2, $3, $4)(m_i, 
                                  blas_no_conj, 
                                  alpha, beta, rin, rout, 
                                  HEAD(r_true)_elem, TAIL(r_true)_elem, 
                                  a_vec, 1, b_vec, 1, eps_int, un_int, 
                                  &ratios[rij]);

                                /* take the max ratio */
                                if (rij == 0) {
                                  ratio = ratios[0];
                                /* The !<= below causes NaN error to be detected.
                                   Note that (NaN > thresh) is always false. */
                                } else if ( !(ratios[rij] <= ratio) ) {
                                  ratio = ratios[rij];
                                }

                              }
                            }  /* end of dot-test loop */

                           /* Increase the number of bad ratio, if the ratio
                              is bigger than the threshold.
                              The !<= below causes NaN error to be detected.
                              Note that (NaN > thresh) is always false. */
                           if ( !(ratio <= thresh) ) {

                              if (debug == 3) {
                                printf("Seed = %d\n", saved_seed);
                                printf("m %d   n %d\n", m, n);
                                printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

                                if (order_type == blas_rowmajor)
                                  printf("row ");
                                else
                                  printf("col ");

                                if (uplo_type == blas_upper)
                                  printf("upper ");
                                else 
                                  printf("lower");

                                if (side_type == blas_left_side)
                                  printf(" left\n");
                                else
                                  printf(" right\n");

                                printf("NORM %d, ALPHA %d, BETA %d\n", 
                                       norm, alpha_val, beta_val);

                                /* print out info */
                                PRINT_VAR(alpha, $1_type);
                                printf("   ");
                                PRINT_VAR(beta, $1_type);
                                printf("\n");

                                printf("a\n");
                                $2sy_print_matrix(a, m_i, lda, order_type, uplo_type);
                                $3ge_print_matrix(b, m, n, ldb, order_type, "B");
                                $1ge_print_matrix(c_gen, m, n, ldc, order_type, "C_gen");
                                $1ge_print_matrix(c, m, n, ldc, order_type, "C");

                                printf("ratio = %g\n", ratio);
                              }
                              bad_ratio_count++;
                              if (bad_ratio_count >= MAX_BAD_TESTS) {
                                printf("\ntoo many failures, exiting....");
                                printf("\nTesting and compilation");
                                printf(" are incomplete\n\n");
                                goto end;
                              }
                              if ( !(ratio <= TOTAL_FAILURE_THRESHOLD) ) {
                                printf("\nFlagrant ratio %e, exiting...", ratio);
                                printf("\nTesting and compilation");
                                printf(" are incomplete\n\n");
                                goto end;
                              }                         
                            }

                            if (ratio > ratio_max)
                              ratio_max = ratio;

                            if (ratio != 0.0 && ratio < ratio_min)
                              ratio_min = ratio;

                        } /* end of randmize loop */

                      } /* end of ldc loop */

                    } /* end of ldb loop */

                  } /* end of lda loop */     

                } /* end of side loop */
                
              } /* end of uplo loop */

            } /* end of order loop */

          } /* end of nr test loop */
          
        } /* end of norm loop */

dnl ---------------------------------------------------------
ifelse(`$4', `_x', `
          } /* end of prec loop */  ')
dnl ---------------------------------------------------------

      } /* end of beta loop */

    } /* end of alpha loop */

end:
  FPU_FIX_STOP;

  FREE_VECTOR(c, $1_type)
  FREE_VECTOR(a, $2_type)
  FREE_VECTOR(b, $3_type)
  FREE_VECTOR(c_gen, $1_type)
  FREE_VECTOR(r_true, EXTRA_TYPE($1_type))
  FREE_VECTOR(ratios, real_D)
  FREE_VECTOR(a_vec, $2_type)
  FREE_VECTOR(b_vec, $3_type)

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;
') 
dnl
dnl
dnl
define(`CALL_DO_TEST_SYMM', 
 `fname = "SYMM_NAME($1, $2, $3, $4)";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mn_data[i][0];
    n = mn_data[i][1];

    DO_TEST_SYMM_NAME($1, $2, $3, $4)(m, n, ntests, &seed, thresh, debug, 
        test_prob, &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", m, n);
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
   
FOREACH(`SYMM_ARGS', `
DO_TEST_SYMM(arg)')dnl

MAIN(`
  int m, i;
  int mn_data[NUM_DATA][2] = 
      {{4, 4}, {2, 3}, {4, 2}, {8, 8}, {10, 1}, {1, 1}, {1, 7}};', `

FOREACH(`SYMM_ARGS', `
CALL_DO_TEST_SYMM(arg)')')dnl

