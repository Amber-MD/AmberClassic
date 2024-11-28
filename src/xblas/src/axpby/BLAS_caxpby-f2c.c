
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_caxpby(int n, const void *alpha, const void *x, int incx,
		 const void *beta, void *y, int incy);


extern void FC_FUNC_(blas_caxpby, BLAS_CAXPBY)
 
  (int *n, const void *alpha, const void *x, int *incx, const void *beta,
   void *y, int *incy) {
  BLAS_caxpby(*n, alpha, x, *incx, beta, y, *incy);
}
