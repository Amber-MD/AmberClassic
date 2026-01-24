
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_sgemv(enum blas_order_type order, enum blas_trans_type trans,
		int m, int n, float alpha, const float *a, int lda,
		const float *x, int incx, float beta, float *y, int incy);


extern void FC_FUNC_(blas_sgemv, BLAS_SGEMV)
 
  (int *trans, int *m, int *n, float *alpha, const float *a, int *lda,
   const float *x, int *incx, float *beta, float *y, int *incy) {
  BLAS_sgemv(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *alpha, a,
	     *lda, x, *incx, *beta, y, *incy);
}
