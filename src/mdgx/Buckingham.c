#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Macros.h"
#include "Buckingham.h"

//-----------------------------------------------------------------------------
// BuckPotential: computes the Buckingham potential given the three    
//                coefficients eps, gam, and sig.  The complexity of the
//                potential is the reason for having a separate function just
//                to calculate it.                        
//                                                                      
//  U(r) = (eps / (1-6/gam)) * [ (6/gam)*exp(gam*(1-r/sig))-(sig/r)^6 ] 
//                                                                      
// Arguments:                                                           
//   r:        the distance at which to compute the potential           
//-----------------------------------------------------------------------------
static double BuckPotential(double eps, double gam, double sig, double r)
{
  double u;

  u = eps/(1-6.0/gam)*( (6.0/gam)*exp(gam*(1.0-r/sig)) - pow(sig/r, 6.0));

  return u;
}

//-----------------------------------------------------------------------------
// BuckDerivative: computes the Buckingham potential derivative given the
//                 coefficients eps, gam, and sig.                  
//                                                                      
// Arguments:                                                           
//   r:        the distance at which to compute the potential           
//-----------------------------------------------------------------------------
static double BuckDerivative(double eps, double gam, double sig, double r)
{
  double du;

  du = eps/(1-6.0/gam)*( -(6.0/sig)*exp(gam*(1.0-r/sig)) +
			 6.0*pow(sig/r, 6.0)/r);

  printf("du = %16.8lf at r = %16.8lf\n", du, r);

  return du;
}

//-----------------------------------------------------------------------------
// BuckDerivative2: computes the Buckingham potential second derivative given
//                  the coefficients eps, gam, and sig.           
//                                                                      
// Arguments:                                                           
//   r:        the distance at which to compute the potential           
//-----------------------------------------------------------------------------
static double BuckDerivative2(double eps, double gam, double sig, double r)
{
  double ddu;

  ddu = eps/(1-6.0/gam)*( (6.0*gam/(sig*sig))*exp(gam*(1.0-r/sig)) -
			  42.0*pow(sig/r, 6.0)/(r*r));

  return ddu;
}

//-----------------------------------------------------------------------------
// BuckDerivative3: computes the Buckingham potential second derivative given
//                  the coefficients eps, gam, and sig.           
//                                                                      
// Arguments:                                                           
//   r:        the distance at which to compute the potential           
//-----------------------------------------------------------------------------
static double BuckDerivative3(double eps, double gam, double sig, double r)
{
  double ddu;

  ddu = eps/(1-6.0/gam)*( (-6.0*gam*gam/(sig*sig*sig))*exp(gam*(1.0-r/sig)) -
			  336.0*pow(sig/r, 6.0)/(r*r*r));

  return ddu;
}

//-----------------------------------------------------------------------------
// DetermineDCoef: this routine computes the final coefficient for the modified
//                 Buckingham potential of the form:           
//                                                                      
//  U(r) = (eps / (1-6/gam)) * [ (6/gam)*exp(gam*(1-r/sig))-(sig/r)^6 ] 
//                                                                      
// Arguments:                                                           
//   [abc]val: the A, B, and C coefficients above                       
//-----------------------------------------------------------------------------
double DetermineDCoef(double eps, double gam, double sig)
{
  int i, j;
  double r, dr, rmin, rmax, d2Ur, d2Udr, dval;

  // Simply return zero if there are indeed no vdW parameters
  if (eps < 1.0e-8 && sig < 1.0e-8) {
    return 0.0;
  }

  // Numerically determine the point at which the unmodified
  // Buckingham potential would reach an inflection point.  
  // First, we find bounds on the range to search by        
  // determining where the unmodified potential reaches a   
  // maximum and later crosses zero.                        
  rmax = -1.0;
  r = 1.0e-4;
  while (rmax < 0.0) {
    if (BuckDerivative(eps, gam, sig, r) < 0.0) {
      rmax = r;
      rmin = MAX(0.0, rmax-0.5);
    }
    r += 0.5;
  }
  for (i = 0; i < 12; i++) {
    dr = 0.1*(rmax-rmin);
    for (j = 0; j < 10; j++) {
      r = rmin + (j+1)*dr;
      if (BuckDerivative(eps, gam, sig, r) < 0.0) {
        rmax = r;
        rmin = MAX(0.0, r-dr);
        dr = 0.1*(rmax-rmin);
	j = 10;
      }
    }
  }
  r = 0.5*(rmax+rmin);
  printf("U(r) reaches a maximum of %12.8lf at %12.8lf\n",
	 BuckPotential(eps, gam, sig, r), r);
  printf("Check: %12.8lf vs. %12.8lf\n", BuckPotential(eps, gam, sig, r-0.01),
	 BuckPotential(eps, gam, sig, r+0.01));

  rmax = -1.0;
  r = rmin;
  while (rmax < 0.0) {
    if (BuckPotential(eps, gam, sig, r) < 0.0) {
      rmax = r;
    }
    r += 0.01;
  }
  printf("Range is %12.8lf -> %12.8lf\n", rmin, rmax);

  for (i = 0; i < 12; i++) {
    dr = 0.1*(rmax-rmin);
    for (j = 0; j < 10; j++) {
      r = rmin + j*dr;
      d2Ur = BuckDerivative2(eps, gam, sig, r);
      d2Udr = BuckDerivative2(eps, gam, sig, r+dr);
      if ((d2Ur <= 0.0 && d2Udr >= 0.0) || (d2Ur >= 0.0 && d2Udr <= 0.0)) {
	rmax = r+dr;
	rmin = r;
	dr = 0.1*(rmax-rmin);
	j = 10;
      }
    }
  }
  r = 0.5*(rmin+rmax);
  printf("Inflection point found at %12.8lf\n", 0.5*(rmin+rmax));
  printf("Check: %12.8lf vs. %12.8lf\n",
	 BuckDerivative2(eps, gam, sig, r-0.01),
         BuckDerivative2(eps, gam, sig, r+0.01));

  printf("Third derivative of the Buckingham potential at r = %16.8lf\n", r);
  printf("is %16.8lf\n", BuckDerivative3(eps, gam, sig, r));

  dval = 0.5*eps*pow(sig, 6.0)/(1.0-(6.0/gam))*pow(r, 6.0);
  printf("D coefficient: %16.8lf\n", dval);

  return dval;
}
