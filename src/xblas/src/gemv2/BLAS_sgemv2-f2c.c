
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_sgemv2(enum blas_order_type order, enum blas_trans_type trans,
		 int m, int n, float alpha, const float *a, int lda,
		 const float *head_x, const float *tail_x, int incx,
		 float beta, float *y, int incy);


extern void FC_FUNC_(blas_sgemv2, BLAS_SGEMV2)
 
  (int *trans, int *m, int *n, float *alpha, const float *a, int *lda,
   const float *head_x, const float *tail_x, int *incx, float *beta, float *y,
   int *incy) {
  BLAS_sgemv2(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *alpha, a,
	      *lda, head_x, tail_x, *incx, *beta, y, *incy);
}
