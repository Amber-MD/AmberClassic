
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cgemm(enum blas_order_type order, enum blas_trans_type transa,
		enum blas_trans_type transb, int m, int n, int k,
		const void *alpha, const void *a, int lda, const void *b,
		int ldb, const void *beta, void *c, int ldc);


extern void FC_FUNC_(blas_cgemm, BLAS_CGEMM)
 
  (int *transa, int *transb, int *m, int *n, int *k, const void *alpha,
   const void *a, int *lda, const void *b, int *ldb, const void *beta,
   void *c, int *ldc) {
  BLAS_cgemm(blas_colmajor, (enum blas_trans_type) *transa,
	     (enum blas_trans_type) *transb, *m, *n, *k, alpha, a, *lda, b,
	     *ldb, beta, c, *ldc);
}
