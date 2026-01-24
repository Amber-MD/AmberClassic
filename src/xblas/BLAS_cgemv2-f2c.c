
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cgemv2(enum blas_order_type order, enum blas_trans_type trans,
		 int m, int n, const void *alpha, const void *a, int lda,
		 const void *head_x, const void *tail_x, int incx,
		 const void *beta, void *y, int incy);


extern void FC_FUNC_(blas_cgemv2, BLAS_CGEMV2)
 
  (int *trans, int *m, int *n, const void *alpha, const void *a, int *lda,
   const void *head_x, const void *tail_x, int *incx, const void *beta,
   void *y, int *incy) {
  BLAS_cgemv2(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, alpha, a,
	      *lda, head_x, tail_x, *incx, beta, y, *incy);
}
