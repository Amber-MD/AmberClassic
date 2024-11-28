
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dgemv(enum blas_order_type order, enum blas_trans_type trans,
   int m, int n, double alpha, const double* a, int lda,
   const double* x, int incx, double beta, double* y, 
   int incy);


extern void FC_FUNC_(blas_dgemv,BLAS_DGEMV)
  ( int*  trans, int* m, int* n, double* alpha, const double* a, int* lda, const double* x, int* incx, double* beta, double* y, int* incy )
{
	 BLAS_dgemv(  blas_colmajor,  (enum blas_trans_type)* trans,
   * m, * n, * alpha,   a, * lda,
     x, * incx, * beta,  y, 
   * incy);
}
