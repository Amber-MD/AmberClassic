
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dwaxpby(int n, double alpha, const double* x, int incx,
   double beta, const double* y, int incy, double* w, 
   int incw);


extern void FC_FUNC_(blas_dwaxpby,BLAS_DWAXPBY)
  ( int* n, double* alpha, const double* x, int* incx, double* beta, const double* y, int* incy, double* w, int* incw )
{
	 BLAS_dwaxpby(* n, * alpha,   x, * incx,
   * beta,   y, * incy,  w, 
   * incw);
}
