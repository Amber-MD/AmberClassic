#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_swaxpby(int n, float alpha, const float *x, int incx,
		  float beta, const float *y, int incy, float *w, int incw)

/*
 * Purpose
 * =======
 *
 * This routine computes:
 *
 *     w <- alpha * x + beta * y
 * 
 * Arguments
 * =========
 *
 * n     (input) int
 *       The length of vectors x, y, and w.
 * 
 * alpha (input) float
 *
 * x     (input) const float*
 *       Array of length n.
 * 
 * incx  (input) int
 *       The stride used to access components x[i].
 *
 * beta  (input) float
 *
 * y     (input) float*
 *       Array of length n.
 *
 * incy  (input) int
 *       The stride used to access components y[i].
 *
 * w     (output) float*
 *       Array of length n.
 *
 * incw  (input) int
 *       The stride used to write components w[i].
 *
 */
{
  char *routine_name = "BLAS_swaxpby";

  int i, ix = 0, iy = 0, iw = 0;
  float *w_i = w;
  const float *x_i = x;
  const float *y_i = y;
  float alpha_i = alpha;
  float beta_i = beta;
  float x_ii;
  float y_ii;
  float tmpx;
  float tmpy;



  /* Test the input parameters. */
  if (incx == 0)
    BLAS_error(routine_name, -4, incx, NULL);
  else if (incy == 0)
    BLAS_error(routine_name, -7, incy, NULL);
  else if (incw == 0)
    BLAS_error(routine_name, -9, incw, NULL);


  /* Immediate return */
  if (n <= 0) {
    return;
  }






  if (incx < 0)
    ix = (-n + 1) * incx;
  if (incy < 0)
    iy = (-n + 1) * incy;
  if (incw < 0)
    iw = (-n + 1) * incw;

  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii = y_i[iy];
    tmpx = alpha_i * x_ii;	/* tmpx  = alpha * x[ix] */
    tmpy = beta_i * y_ii;	/* tmpy = beta * y[iy] */
    tmpy = tmpy + tmpx;
    w_i[iw] = tmpy;
    ix += incx;
    iy += incy;
    iw += incw;
  }				/* endfor */



}
