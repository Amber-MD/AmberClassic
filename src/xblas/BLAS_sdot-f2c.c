
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_sdot(enum blas_conj_type conj, int n, float alpha,
	       const float *x, int incx, float beta,
	       const float *y, int incy, float *r);


extern void FC_FUNC_(blas_sdot, BLAS_SDOT)
 
  (int *conj, int *n, float *alpha, const float *x, int *incx, float *beta,
   const float *y, int *incy, float *r) {
  BLAS_sdot((enum blas_conj_type) *conj, *n, *alpha, x, *incx, *beta, y,
	    *incy, r);
}
