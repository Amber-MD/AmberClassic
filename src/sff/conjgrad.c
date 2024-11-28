#include <stdio.h>
#include "sff.h"
#include "memutil.h"
#include "timer.h"

REAL_T seconds();


/***********************************************************************
                            CONJGRAD()
************************************************************************/

/*
**  ---routine carry out conjugate gradient optimization
**
**     x(n)     contains the variables to be updated.
**     n        is the number of variables.
**     f        will contain the final value of the objective function.
**     func     is the name of the function that computes the objective
**                   function and its gradient.
**     rsmgrad  return when rms gradient is less than rmsgrad.
**     dfpred   expected decrease in the function on the first iteration.
**
**     return codes:
**             >0    converged, final iteration number
**             -1    bad line search, probably an error in the relation
**                     of the funtion to its gradient (perhaps from
**                     round-off if you push too hard on the minimization).
**             -2    search direction was uphill
**             -3    exceeded the maximum number of iterations
**             -4    could not further reduce function value
**
*/

int conjgrad( REAL_T x[], int *n, REAL_T *f,
	      REAL_T ( *func )( REAL_T*, REAL_T*, int* ),
	      REAL_T *rmsgrad, REAL_T *dfpred, int *maxiter ) 
{
  REAL_T *w, *g;
  int maxlin = 10;
  int mxfcon = 4;
    
  int ret_val;
  REAL_T r__1, r__2, r__3;
  REAL_T dgrad;
    
  REAL_T gama, beta, fmin, gmin, dfpr, gnew;
  REAL_T step, work, finit, ginit, gsqrd, gspln;
  REAL_T stmin, gamden, ddspln, stepch, sbound;
  REAL_T fch, sum;
  int niter, niterm1, i, nfbeg, iterc, irsdg, igopt;
  int nfopt, irsdx, ixopt, iginit;
  int iterfm, iterrs, iretry, ier;
  int ncopy;

  /* All arrays, even the x array, will be indexed from 1 to *n. */
    
  --x;

  ncopy = *n;
  w = vector( 1, 6*ncopy );
  g = vector( 1, ncopy );
  dgrad = (*rmsgrad) * (*rmsgrad) * ncopy;
  ier = 0;

  irsdx = *n;
  irsdg = irsdx + *n;
  iginit = irsdg + *n;
  ixopt = iginit + *n;
  igopt = ixopt + *n;

  iterc = 0;
  niter = 0;
  iterfm = iterc;

 L10:
  ++niter;

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

  /*
  **   Compute the function value in "f" and its
  **   gradient (with respect to x) in "g".
  */

  niterm1 = niter - 1;

  *f = ( *func )( &x[ 1 ], &g[ 1 ], &niterm1 );

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

  if (niter >= 2) {
    goto L30;
  }
 L20:
  for (i = 1; i <= *n; ++i) w[i] = -g[i];
  iterrs = 0;
  if (iterc > 0) goto L40;

 L30:
  gnew = (REAL_T)0.;
  sum = (REAL_T)0.;
  for (i = 1; i <= *n; ++i) {
    gnew += w[i] * g[i];
    r__1 = g[i];
    sum += r__1 * r__1;
  }
  if (niter == 1) {
    fmin = *f;
    gsqrd = sum;
    nfopt = niter;
    for (i = 1; i <= *n; ++i) {
      w[ixopt + i] = x[i];
      w[igopt + i] = g[i];
    }
    if (sum <= dgrad) goto L100;
  } else {
    fch = *f - fmin;

    if (fch < (REAL_T)0. || (fch == (REAL_T)0. && gnew / gmin >= (REAL_T)-1.)) {
      fmin = *f;
      gsqrd = sum;
      nfopt = niter;
      for (i = 1; i <= *n; ++i) {
	w[ixopt + i] = x[i];
	w[igopt + i] = g[i];
      }
      if (sum <= dgrad) goto L100;
    }
  }
  if (niter == *maxiter + 1) { ier = -3; goto L100; }
  if (niter > 1) goto L60;
  dfpr = *dfpred;
  stmin = *dfpred / gsqrd;

 L40:
  ++iterc;

  finit = *f;
  ginit = (REAL_T)0.;
  for (i = 1; i <= *n; ++i) {
    w[iginit + i] = g[i];
    ginit += w[i] * g[i];
  }
  if (ginit >= (REAL_T)0.) { ier = -2; goto L90; }
  gmin = ginit;
  sbound = (REAL_T)-1.;
  nfbeg = niter;
  iretry = -1;
  r__2 = stmin, r__3 = (r__1 = dfpr / ginit, ((r__1) >= 0 ? (r__1) : -(r__1)));
  stepch =  (r__2) <= (r__3) ? (r__2) : (r__3);
  stmin = (REAL_T)0.;

 L50:
  step = stmin + stepch;
  work = (REAL_T)0.;
  for (i = 1; i <= *n; ++i) {
    x[i] = w[ixopt + i] + stepch * w[i];
    r__2 = work, r__3 = (r__1 = x[i] - w[ixopt + i], ((r__1) >= 0 ? (r__1) : -(r__1)));
    work = ((r__2) >= (r__3) ? (r__2) : (r__3));
  }
  if (work > (REAL_T)0.) goto L10;
  if (niter > nfbeg + 1) ier = -1;
  if ((r__1 = gmin / ginit, ((r__1) >= 0 ? (r__1) : -(r__1))) > (REAL_T).2)
    ier = -1;
  goto L90;

 L60:
  work = (fch + fch) / stepch - gnew - gmin;
  ddspln = (gnew - gmin) / stepch;
  if (niter > nfopt) {
    sbound = step;
  }
  if (niter <= nfopt) {
    if (gmin * gnew <= (REAL_T)0.) sbound = stmin;
    stmin = step;
    gmin = gnew;
    stepch = -stepch;
  }
  if (fch != (REAL_T)0.) ddspln += (work + work) / stepch;

  if (gmin == (REAL_T)0.) goto L90;
  if (niter > nfbeg + 1) {
    if ((r__1 = gmin / ginit, ((r__1) >= 0 ? (r__1) :
			       -(r__1))) <= (REAL_T).2) {
      goto L90;
    }
    if (niter >= nfopt + maxlin) { ier = -1; goto L90; }
  }

 L70:
  stepch = (sbound - stmin) * (REAL_T).5;
  if (sbound < (REAL_T)-.5) {
    stepch = stmin * (REAL_T)9.;
  }
  gspln = gmin + stepch * ddspln;
  if (gmin * gspln < (REAL_T)0.) stepch = stepch * gmin / (gmin - gspln);
  goto L50;

 L80:
  sum = (REAL_T)0.;
  for (i = 1; i <= *n; ++i) sum += g[i] * w[iginit + i];
  beta = (gsqrd - sum) / (gmin - ginit);

  if ((r__1 = beta * gmin, ((r__1) >= 0 ? (r__1) : -(r__1))) > gsqrd *
      (REAL_T).2) {
    ++iretry;
    if (iretry <= 0) {
      if (niter >= nfopt + maxlin) { ier = -1; goto L90; }
      else goto L70;
    }
  }

  if (*f < finit) iterfm = iterc;
  if (iterc >= iterfm + mxfcon) { ier = -4; goto L100; }
  dfpr = stmin * ginit;

  if (iretry > 0) {
    goto L20;
  }
  if (iterrs != 0 && iterc - iterrs < *n && (sum >= 0 ? sum : -sum) <
      gsqrd * (REAL_T).2) {

    gama = (REAL_T)0.;
    sum = (REAL_T)0.;
    for (i = 1; i <= *n; ++i) {
      gama += g[i] * w[irsdg + i];
      sum += g[i] * w[irsdx + i];
    }
    gama /= gamden;
    if ((r__1 = beta * gmin + gama * sum, ((r__1) >= 0 ? 
					   (r__1) : -(r__1))) < gsqrd * (REAL_T).2) {
      for (i = 1; i <= *n; ++i)
	w[i] = -g[i] + beta * w[i] + gama * w[irsdx + i];
      goto L40;
    }
  }

  gamden = gmin - ginit;
  for (i = 1; i <= *n; ++i) {
    w[irsdx + i] = w[i];
    w[irsdg + i] = g[i] - w[iginit + i];
    w[i] = -g[i] + beta * w[i];
  }
  iterrs = iterc;
  goto L40;

 L90:
  if (niter != nfopt) {
    *f = fmin;
    for (i = 1; i <= *n; ++i) {
      x[i] = w[ixopt + i];
      g[i] = w[igopt + i];
    }
  }
  if (ier == 0) goto L80;
 L100:

  free_vector( w, 1, 6*ncopy );  free_vector( g, 1, ncopy );
  if( ier == 0 ) ret_val = niter;
  else ret_val = ier;

  return ret_val;
} 

