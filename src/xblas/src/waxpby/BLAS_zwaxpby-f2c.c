
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zwaxpby(int n, const void *alpha, const void *x, int incx,
		  const void *beta, const void *y, int incy, void *w,
		  int incw);


extern void FC_FUNC_(blas_zwaxpby, BLAS_ZWAXPBY)
 
  (int *n, const void *alpha, const void *x, int *incx, const void *beta,
   const void *y, int *incy, void *w, int *incw) {
  BLAS_zwaxpby(*n, alpha, x, *incx, beta, y, *incy, w, *incw);
}
