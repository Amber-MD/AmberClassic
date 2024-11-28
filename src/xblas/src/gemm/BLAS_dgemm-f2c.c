
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dgemm(enum blas_order_type order, enum blas_trans_type transa,
   enum blas_trans_type transb, int m, int n, int k,
   double alpha, const double* a, int lda, const double* b, int ldb,
   double beta, double* c, int ldc);


extern void FC_FUNC_(blas_dgemm,BLAS_DGEMM)
  ( int*  transa, int*  transb, int* m, int* n, int* k, double* alpha, const double* a, int* lda, const double* b, int* ldb, double* beta, double* c, int* ldc )
{
	 BLAS_dgemm(  blas_colmajor,  (enum blas_trans_type)* transa,
    (enum blas_trans_type)* transb, * m, * n, * k,
   * alpha,   a, * lda,   b, * ldb,
   * beta,  c, * ldc);
}
