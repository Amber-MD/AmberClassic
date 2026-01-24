
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_swaxpby(int n, float alpha, const float *x, int incx,
		  float beta, const float *y, int incy, float *w, int incw);


extern void FC_FUNC_(blas_swaxpby, BLAS_SWAXPBY)
 
  (int *n, float *alpha, const float *x, int *incx, float *beta,
   const float *y, int *incy, float *w, int *incw) {
  BLAS_swaxpby(*n, *alpha, x, *incx, *beta, y, *incy, w, *incw);
}
