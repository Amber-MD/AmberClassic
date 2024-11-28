
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dgemv2(enum blas_order_type order, enum blas_trans_type trans,
   int m, int n, double alpha, const double* a, int lda,
   const double* head_x, const double* tail_x, int incx, 
   double beta, double* y, int incy);


extern void FC_FUNC_(blas_dgemv2,BLAS_DGEMV2)
  ( int*  trans, int* m, int* n, double* alpha, const double* a, int* lda, const double* head_x, const double* tail_x, int* incx, double* beta, double* y, int* incy )
{
	 BLAS_dgemv2(  blas_colmajor,  (enum blas_trans_type)* trans,
   * m, * n, * alpha,   a, * lda,
     head_x,   tail_x, * incx, 
   * beta,  y, * incy);
}
