
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cdot(enum blas_conj_type conj, int n, const void *alpha,
	       const void *x, int incx, const void *beta,
	       const void *y, int incy, void *r);


extern void FC_FUNC_(blas_cdot, BLAS_CDOT)
 
  (int *conj, int *n, const void *alpha, const void *x, int *incx,
   const void *beta, const void *y, int *incy, void *r) {
  BLAS_cdot((enum blas_conj_type) *conj, *n, alpha, x, *incx, beta, y, *incy,
	    r);
}
