
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ddot(enum blas_conj_type conj, int n, double alpha, 
   const double* x, int incx, double beta,
   const double* y, int incy, 
   double* r);


extern void FC_FUNC_(blas_ddot,BLAS_DDOT)
  ( int*  conj, int* n, double* alpha, const double* x, int* incx, double* beta, const double* y, int* incy, double* r )
{
	 BLAS_ddot( (enum blas_conj_type)* conj, * n, * alpha, 
     x, * incx, * beta,
     y, * incy, 
    r);
}
