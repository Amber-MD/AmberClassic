
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zdot(enum blas_conj_type conj, int n, const void *alpha,
	       const void *x, int incx, const void *beta,
	       const void *y, int incy, void *r);


extern void FC_FUNC_(blas_zdot, BLAS_ZDOT)
 
  (int *conj, int *n, const void *alpha, const void *x, int *incx,
   const void *beta, const void *y, int *incy, void *r) {
  BLAS_zdot((enum blas_conj_type) *conj, *n, alpha, x, *incx, beta, y, *incy,
	    r);
}
