
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_saxpby(int n, float alpha, const float *x, int incx,
		 float beta, float *y, int incy);


extern void FC_FUNC_(blas_saxpby, BLAS_SAXPBY)
 
  (int *n, float *alpha, const float *x, int *incx, float *beta, float *y,
   int *incy) {
  BLAS_saxpby(*n, *alpha, x, *incx, *beta, y, *incy);
}
