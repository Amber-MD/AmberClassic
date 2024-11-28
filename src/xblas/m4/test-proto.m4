#ifndef BLAS_EXTENDED_TEST_H
#define BLAS_EXTENDED_TEST_H
define(`prototypes_only')dnl

dnl #define TEST_PROBABILITY_DEFAULT 0.001
#define MAX_BAD_TESTS 100
#define TOTAL_FAILURE_THRESHOLD 1000

double power(int i1, int i2);
double xrand(int *is);
int FixedBits(double r_true_l, double r_true_t);

void ddmuld(double dda_l, double dda_t, double db, double *ddc_l, dnl
    double *ddc_t);
void ddadd(double dda_l, double dda_t, double ddb_l, double ddb_t, dnl
    double *ddc_l, double *ddc_t);
void dddiv(double dda_l, double dda_t, double ddb_l, double ddb_t, dnl
    double *ddc_l, double *ddc_t);
void z_dddivd(double *dda_l, double *dda_t, double *db, dnl
    double *ddc_l, double *ddc_t);

void testgen_BLAS_sdot(int n, int n_fix2, int n_mix, int norm, dnl
    enum blas_conj_type conj, float *alpha, int alpha_flag, dnl
    float *beta, int beta_flag, float *x, float *y, int *seed, dnl
    float *r, double *r_true_l, double *r_true_t);

void testgen_BLAS_cdot(int n, int n_fix2, int n_mix, int norm, dnl
    enum blas_conj_type conj, void *alpha, int alpha_flag, dnl
    void *beta, int beta_flag, void *x, void *y, int *seed, dnl
    void *r, double r_true_l[], double r_true_t[]);

void testgen_BLAS_ddot(int n, int n_fix2, int n_mix, int norm, dnl
    enum blas_conj_type conj, double *alpha, int alpha_flag, dnl
    double *beta, int beta_flag, double *x, double *y, int *seed, dnl
    double *r, double *r_true_l, double *r_true_t);

void testgen_BLAS_zdot(int n, int n_fix2, int n_mix, int norm, dnl
    enum blas_conj_type conj, void *alpha, int alpha_flag, dnl
    void *beta, int beta_flag, void *x, void *y, int *seed, dnl
    void *r, double r_true_l[], double r_true_t[]);

void testgen_BLAS_sdot2(int n, int n_fix2, int n_mix, int norm, dnl
   enum blas_conj_type conj, float *alpha, int alpha_flag, dnl
   float *beta, int beta_flag, float *head_x, float *tail_x, dnl
   float *y, int *seed, float *r, double *r_true_l, double *r_true_t);

void testgen_BLAS_cdot2(int n, int n_fix2, int n_mix, int norm, dnl
   enum blas_conj_type conj, void *alpha, int alpha_flag, dnl
   void *beta, int beta_flag, void *head_x, void *tail_x, dnl
   void *y, int *seed, void *r, double r_true_l[], double r_true_t[]);

void testgen_BLAS_ddot2(int n, int n_fix2, int n_mix, int norm, dnl
   enum blas_conj_type conj, double *alpha, int alpha_flag, dnl
   double *beta, int beta_flag, double *head_x, double *tail_x, dnl
   double *y, int *seed, double *r, double *r_true_l, double *r_true_t);

void testgen_BLAS_zdot2(int n, int n_fix2, int n_mix, int norm, dnl
   enum blas_conj_type conj, void *alpha, int alpha_flag, dnl
   void *beta, int beta_flag, void *head_x, void *tail_x, dnl
   void *y, int *seed, void *r, double r_true_l[], double r_true_t[]);

include(test-dot/BLAS_dot_testgen.m4)dnl
include(test-dot/copy_vector.m4)dnl
include(test-dot/print_vector.m4)dnl
include(test-dot/test_dot.m4)dnl

include(test-dot2/BLAS_dot2_testgen.m4)dnl
include(test-dot2/dot2.m4)dnl
include(test-dot2/r_truth2.m4)dnl
include(test-dot2/test_dot2.m4)dnl

include(test-gemv/BLAS_gemv_testgen.m4)dnl
include(test-gemv/gemv-support.m4)dnl

include(test-gbmv/BLAS_gbmv_testgen.m4)dnl
include(test-gbmv/gbmv-support.m4)dnl

include(test-gemm/BLAS_gemm_testgen.m4)dnl

include(test-gemv2/BLAS_gemv2_testgen.m4)dnl

include(test-symv2/BLAS_symv2_testgen.m4)dnl

include(test-hemv2/BLAS_hemv2_testgen.m4)dnl

include(test-gbmv2/BLAS_gbmv2_testgen.m4)dnl

include(test-ge_sum_mv/BLAS_ge_sum_mv_testgen.m4)dnl

include(test-symv/BLAS_symv_testgen.m4)dnl
include(test-symv/symv-support.m4)dnl

include(test-sbmv/BLAS_sbmv_testgen.m4)dnl
include(test-sbmv/sbmv-support.m4)dnl

include(test-hemv/BLAS_hemv_testgen.m4)dnl
include(test-hemv/hemv-support.m4)dnl

include(test-hbmv/BLAS_hbmv_testgen.m4)dnl
include(test-hbmv/hbmv-support.m4)dnl

include(test-symm/BLAS_symm_testgen.m4)dnl

include(test-hemm/BLAS_hemm_testgen.m4)dnl

include(test-hpmv/BLAS_hpmv_testgen.m4)dnl
include(test-hpmv/hpmv-support.m4)dnl

include(test-spmv/BLAS_spmv_testgen.m4)dnl
include(test-spmv/spmv-support.m4)dnl

include(test-sum/BLAS_sum_testgen.m4)dnl
include(test-sum/sum-support.m4)dnl

void BLAS_sdot_x_testgen(int n, int n_fix2, int n_mix, int norm, dnl
    enum blas_conj_type conj, float *alpha, int alpha_flag, dnl
    float *beta, int beta_flag, double *x_l, double *x_t, dnl
    float *y, int *seed, float *r, double *r_true_l, double *r_true_t);

void BLAS_ddot_x_testgen(int n, int n_fix2, int n_mix, int norm, dnl
    enum blas_conj_type conj, double *alpha, int alpha_flag, dnl
    double *beta, int beta_flag, double *x_l, double *x_t, dnl
    double *y, int *seed, double *r, double *r_true_l, double *r_true_t);

void BLAS_ddot_s_x_testgen(int n, int n_fix2, int n_mix, int norm, dnl
    enum blas_conj_type conj, double *alpha, int alpha_flag, dnl
    double *beta, int beta_flag, double *x_l, double *x_t, dnl
    float *y, int *seed, double *r, double *r_true_l, double *r_true_t);

void testgen_BLAS_sdot_x(int n, int n_fix2, int n_mix, int norm, dnl
    enum blas_conj_type conj, float *alpha, int alpha_flag, dnl
    float *beta, int beta_flag, double *x_l, double *x_t, float *y, dnl
    int *seed, float *r, double *r_true_l, double *r_true_t);

include(test-trsv/BLAS_trsv_testgen.m4)dnl
include(test-trsv/trsv-support.m4)dnl

include(test-tbsv/BLAS_tbsv_testgen.m4)dnl
include(test-tbsv/tbsv-support.m4)dnl

include(test-tpmv/BLAS_tpmv_testgen.m4)dnl
include(test-tpmv/tpmv-support.m4)dnl

include(test-trmv/BLAS_trmv_testgen.m4)dnl
include(test-trmv/trmv-support.m4)dnl

#endif /* BLAS_EXTENDED_TEST_H */
