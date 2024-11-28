
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cwaxpby(int n, const void *alpha, const void *x, int incx,
		  const void *beta, const void *y, int incy, void *w,
		  int incw);


extern void FC_FUNC_(blas_cwaxpby, BLAS_CWAXPBY)
 
  (int *n, const void *alpha, const void *x, int *incx, const void *beta,
   const void *y, int *incy, void *w, int *incw) {
  BLAS_cwaxpby(*n, alpha, x, *incx, beta, y, *incy, w, *incw);
}
