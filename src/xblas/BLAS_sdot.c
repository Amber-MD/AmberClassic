#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_sdot(enum blas_conj_type conj, int n, float alpha,
	       const float *x, int incx, float beta,
	       const float *y, int incy, float *r)

/*
 * Purpose
 * =======
 * 
 * This routine computes the inner product:
 * 
 *     r <- beta * r + alpha * SUM_{i=0, n-1} x[i] * y[i].
 * 
 * Arguments
 * =========
 *  
 * conj   (input) enum blas_conj_type
 *        When x and y are complex vectors, specifies whether vector
 *        components x[i] are used unconjugated or conjugated. 
 * 
 * n      (input) int
 *        The length of vectors x and y.
 * 
 * alpha  (input) float
 * 
 * x      (input) const float*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * beta   (input) float
 *
 * y      (input) const float*
 *        Array of length n.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input/output) float*
 * 
 */
{
  static const char routine_name[] = "BLAS_sdot";

  int i, ix = 0, iy = 0;
  float *r_i = r;
  const float *x_i = x;
  const float *y_i = y;
  float alpha_i = alpha;
  float beta_i = beta;
  float x_ii;
  float y_ii;
  float r_v;
  float prod;
  float sum;
  float tmp1;
  float tmp2;


  /* Test the input parameters. */
  if (n < 0)
    BLAS_error(routine_name, -2, n, NULL);
  else if (incx == 0)
    BLAS_error(routine_name, -5, incx, NULL);
  else if (incy == 0)
    BLAS_error(routine_name, -8, incy, NULL);

  /* Immediate return. */
  if ((beta_i == 1.0) && (n == 0 || (alpha_i == 0.0)))
    return;



  r_v = r_i[0];
  sum = 0.0;


  if (incx < 0)
    ix = (-n + 1) * incx;
  if (incy < 0)
    iy = (-n + 1) * incy;

  for (i = 0; i < n; ++i) {
    x_ii = x_i[ix];
    y_ii = y_i[iy];

    prod = x_ii * y_ii;		/* prod = x[i]*y[i] */
    sum = sum + prod;		/* sum = sum+prod */
    ix += incx;
    iy += incy;
  }				/* endfor */


  tmp1 = sum * alpha_i;		/* tmp1 = sum*alpha */
  tmp2 = r_v * beta_i;		/* tmp2 = r*beta */
  tmp1 = tmp1 + tmp2;		/* tmp1 = tmp1+tmp2 */
  *r = tmp1;			/* r = tmp1 */



}
