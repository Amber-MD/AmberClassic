dnl Generates test code for gemm
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"
include(cblas.m4)dnl
include(test-common.m4)dnl
include(gemm/gemm-common.m4)dnl
dnl
dnl
dnl -------------------------------------------------
dnl Usage: DO_TEST_GEMM($1, $2, $3, $4)
dnl   $1 -- type of alpha, beta, and matrix C
dnl   $2 -- type of matrix A
dnl   $3 -- type of matrix B
dnl   $4 -- `_x' for extended precision, `' otherwise
dnl -------------------------------------------------
define(`DO_TEST_GEMM_NAME', 
  `ifelse(`$1&&$1', `$2&&$3', `do_test_$1gemm$4', `do_test_$1gemm_$2_$3$4')')
dnl
dnl
define(`DO_TEST_GEMM_PARAMS', 
  `int m, int n, int k, dnl
   int ntests, int *seed, double thresh, int debug, float test_prob, dnl
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests')
dnl
dnl
define(`DO_TEST_GEMM_HEAD', 
  `void DO_TEST_GEMM_NAME($1, $2, $3, $4)(DO_TEST_GEMM_PARAMS($1, $2, $3, $4))')
dnl
dnl
define(`DO_TEST_GEMM', 
  `DO_TEST_GEMM_HEAD($1, $2, $3, $4) {
   DO_TEST_GEMM_BODY($1, $2, $3, $4)
  }')
dnl
dnl
define(`TESTGEN_GEMM_NAME', `ifelse(
       `$2&&$3', `$1&&$1', `BLAS_$1gemm_testgen',
       `BLAS_$1gemm_$2_$3_testgen')') dnl
dnl
dnl
define(`DO_TEST_GEMM_BODY', 
  `
  /* Function name */
  const char fname[] = "GEMM_NAME($1, $2, $3, $4)";

  int i, j;
  int ci, cij;
  int ri, rij;
  int incci, inccij;
  int incri, incrij;
  int inca, incb, incc;

  int test_count;      /* number of tests done so far   */
  int bad_ratio_count; /* number of failed tests so far */

  double ratio;
  double ratio_min, ratio_max;

  double eps_int;  /* internal machine epsilon     */
  double un_int;   /* internal underflow threshold */

  enum blas_order_type order;
  enum blas_trans_type transa;
  enum blas_trans_type transb;
  enum blas_prec_type  prec;

  int order_val, transa_val, transb_val;
  int lda_val, ldb_val, ldc_val;
  int alpha_val, beta_val;
  int randomize_val;
  
  ifelse(`$4', `_x', `int prec_val;', `')
  
  int lda, ldb, ldc;
  
  int lda_0, ldb_0;
  int tda_0, tdb_0;
  int lda_1, ldb_1, ldc_1;
  
  int alpha_flag, beta_flag;

  int saved_seed;

  int norm;

  int test_no;

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

  if (n < 0 || m < 0 || k < 0 || ntests < 0)
    BLAS_error(fname, 0, 0, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *max_ratio = 0.0;
  *min_ratio = 0.0;
  *num_bad_ratio = 0;
  *num_tests = 0;

  if (n == 0 || m == 0 || k == 0)
    return;  

  FPU_FIX_START;

  inca = incb = incc = 1;
  INC_ADJUST(inca, $2_type)
  INC_ADJUST(incb, $3_type)
  INC_ADJUST(incc, $1_type)

  /* allocate memory for arrays */
  MALLOC_VECTOR(c, $1_type, 2 * m * n)
  MALLOC_VECTOR(c_gen, $1_type, 2 * m * n)
  MALLOC_VECTOR(a, $2_type, 2 * m * k)
  MALLOC_VECTOR(b, $3_type, 2 * k * n)
  MALLOC_VECTOR(a_vec, $2_type, k)
  MALLOC_VECTOR(b_vec, $3_type, k)
  MALLOC_VECTOR(r_true, EXTRA_TYPE($1_type), 2*m*n)
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

              if (order_val == 0) {
                /* row major storage */
                lda_0 = k; ldb_0 = n; ldc_1 = n;
                tda_0 = m; tdb_0 = k;
                order = blas_rowmajor;
              } else {
                /* column major storage */
                lda_0 = m; ldb_0 = k; ldc_1 = m;
                tda_0 = k; tdb_0 = n;
                order = blas_colmajor;
              }

              /* vary transpositions of A */
              for (transa_val = TRANSA_START; transa_val <= TRANSA_END; 
                   transa_val++) {

                transa = (transa_val == 0) ? blas_no_trans : 
                  (transa_val == 1) ? blas_trans : blas_conj_trans;

                if (transa == blas_no_trans) {
                  lda_1 = lda_0;
                } else {
                  lda_1 = tda_0;
                }
                  
                /* vary transpositions of B */
                for (transb_val = TRANSB_START; transb_val <= TRANSB_END; 
                     transb_val++) {

                  transb = (transb_val == 0) ? blas_no_trans :
                    (transb_val == 1) ? blas_trans : blas_conj_trans;

                  if (transb == blas_no_trans) {
                    ldb_1 = ldb_0;
                  } else {
                    ldb_1 = tdb_0;
                  }

                  /* vary lda = k, k+1, 2*k */
                  for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {
                   
                    lda = (lda_val == 0) ? lda_1 :
                      (lda_val == 1) ? lda_1+1 : 2*lda_1;
                      
                    /* vary ldb = n, n+1, 2*n */
                    for (ldb_val = LDB_START; ldb_val <= LDB_END; ldb_val++) {
                      
                      ldb = (ldb_val == 0) ? ldb_1 :
                        (ldb_val == 1) ? ldb_1+1 : 2*ldb_1;

                      /* vary ldc = k, k+1, 2*k */
                      for (ldc_val = LDC_START; ldc_val <= LDC_END; ldc_val++) {
                        
                        ldc = (ldc_val == 0) ? ldc_1 :
                          (ldc_val == 1) ? ldc_1+1 : 2*ldc_1;
                        
                        for (randomize_val = RANDOMIZE_START; 
                             randomize_val <= RANDOMIZE_END; randomize_val++) {

                           /* For the sake of speed, we throw out this case at random */
                           if ( xrand(seed) >= test_prob ) continue;
                           
                           /* finally we are here to generate the test case */
                           TESTGEN_GEMM_NAME($1, $2, $3) (norm, order, 
                             transa, transb, m, n, k, randomize_val, &alpha, 
                             alpha_flag, a, lda, &beta, beta_flag, b, ldb, 
                             c, ldc, seed, HEAD(r_true), TAIL(r_true));
                           test_count++;

                           /* copy generated C matrix since this will be
                              over written */
                           $1ge_copy_matrix(order, m, n, c_gen, ldc, c, ldc);

                           /* call GEMM routines to be tested */
                           FPU_FIX_STOP;
                           GEMM_NAME($1, $2, $3, $4)(order, transa, 
                             transb, m, n, k, alpha, a, lda, b, ldb, beta, 
                             c, ldc ifelse(`$4', `_x', `, prec'));
                           FPU_FIX_START;

                           /* now compute the ratio using test_c_xdot */
                           /* copy a row from A, a column from B, run 
                              dot test */

                           if (order == blas_colmajor) {
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

                           for (i = 0, ci = 0, ri = 0; i < m; 
                                i++, ci += incci, ri += incri) {
                             $2ge_copy_row(order, transa, 
                                             m, k, a, lda, a_vec, i);
                             for (j = 0, cij = ci, rij = ri; j < n; 
                                  j++, cij += inccij, rij += incrij) {
                               /* copy i-th row of A and j-th col of B */
                               $3ge_copy_col(order, transb, 
                                               k, n, b, ldb, b_vec, j);

                               TEST_DOT_NAME($1, $2, $3, $4)(k, blas_no_conj, 
                                 alpha, beta, VECTOR_ELEMENT(c_gen, cij, $1_type), 
                                 VECTOR_ELEMENT(c, cij, $1_type), 
                                 VECTOR_ELEMENT(r_true, cij, EXTRA_TYPE($1_type)), 
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

                               printf("\nm %d   n %d   k %d\n", m, n, k);
                               printf("LDA %d  LDB %d  LDC %d\n", lda, ldb, ldc);

                               PRINT_TRANS(transa, A)
                               PRINT_TRANS(transb, B)

                               printf("NORM %d, ALPHA %d, BETA %d\n", 
                                      norm, alpha_val, beta_val);

                               PRINT_ORDER(order)
                               ifelse(`$4', `_x', `PRINT_PREC(prec)')

                               if (randomize_val == 0)
                                 printf("Not randomized\n");
                               else 
                                 printf("Randomized\n");

                               /* print out info */
                               PRINT_VAR(alpha, $1_type);
                               printf("   ");
                               PRINT_VAR(beta, $1_type);
                               printf("\n");

                               $2ge_print_matrix(a, m, k, lda, order, "A");
                               $3ge_print_matrix(b, k, n, ldb, order, "B");
                               $1ge_print_matrix(c_gen, m, n, ldc, order, "C_gen");
                               $1ge_print_matrix(c, m, n, ldc, order, "C");
                               DOUBLE_ABBREV($1)ge_print_matrix(HEAD(r_true), m, n, ldc, order, "truth");

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

                        } /* end of randomize loop */

                      } /* end of ldc loop */

                    } /* end of ldb loop */

                  } /* end of lda loop */     

                } /* end of transb loop */
                
              } /* end of transa loop */

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
define(`CALL_DO_TEST_GEMM', 
 `fname = "GEMM_NAME($1, $2, $3, $4)";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    m = mnk_data[i][0];
    n = mnk_data[i][1];
    k = mnk_data[i][2];

    DO_TEST_GEMM_NAME($1, $2, $3, $4)(m, n, k, ntests, &seed, thresh, debug, 
      test_prob, &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d %d]: ", m, n, k);
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
dnl   
dnl

/* 0 -- 2 */
#define TRANSA_START    0
#define TRANSA_END      2

/* 0 -- 2 */
#define TRANSB_START    0
#define TRANSB_END      2

/* 0 -- 1 */
#define ORDER_START     0
#define ORDER_END       1

/* 0 -- 2 */
#define ALPHA_START     0
#define ALPHA_END       2

/* 0 -- 2 */
#define BETA_START      0
#define BETA_END        2

/* -1 -- 1 */
#define NORM_START     -1
#define NORM_END        1

/* 0 -- 2 */
#define LDA_START       0
#define LDA_END         2

/* 0 -- 2 */
#define LDB_START       0 
#define LDB_END         2

/* 0 -- 2 */
#define LDC_START       0
#define LDC_END         2

/* 0 -- 2 */
#define PREC_START      0
#define PREC_END        2

/* 0 -- 1 */
#define RANDOMIZE_START 0
#define RANDOMIZE_END   1

#define NUM_DATA 9

FOREACH(`GEMM_ARGS', `
DO_TEST_GEMM(arg)')dnl

MAIN(`
  int m, k, i;
  int mnk_data[NUM_DATA][3] = {{2, 3, 4},  {3, 6, 4}, {5, 1, 7}, 
  {6, 15, 4}, {5, 2, 3}, {8, 4, 1}, 
  {1, 3, 1}, {8, 8, 8}, {1, 1, 1}};', `

FOREACH(`GEMM_ARGS', `
CALL_DO_TEST_GEMM(arg)')')dnl

