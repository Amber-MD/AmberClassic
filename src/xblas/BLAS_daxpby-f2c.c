
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_daxpby(int n, double alpha, const double* x, int incx, 
   double beta, double* y, 
   int incy);


extern void FC_FUNC_(blas_daxpby,BLAS_DAXPBY)
  ( int* n, double* alpha, const double* x, int* incx, double* beta, double* y, int* incy )
{
	 BLAS_daxpby(* n, * alpha,   x, * incx, 
   * beta,  y, 
   * incy);
}
