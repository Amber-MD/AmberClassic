#include <math.h>
#include <stdio.h>
#include "Constants.h"
#include "SpecialMath.h"

#include "CompFrcDS.h"

//-----------------------------------------------------------------------------
// perfc: this erfc() function was taken from sander source code _erfcfun.f.
//        The algorithm is based on the observation that, for the portion of
//        the domain where erfc() changes rapidly, it appears very much like a
//        product of a Gaussian and a rational function.  Computing these, and
//        optimizing the rational function with seventh-order polynomials,
//        produces the desired result to machine precision.
//
// Arguments:                                                           
//   x:  the argument of the complimentary error function               
//-----------------------------------------------------------------------------
double perfc(double x)
{
  double absx, c, p, q, nonexperfc, tserf, tserfc;

  absx = fabs(x);
  if (x > 26.0) {
    tserfc = 0.0;
  }
  else if (x < -5.5) {
    tserfc = 2.0;
  }
  else if (absx <= 0.5) {
    c = x*x;
    p = (((-0.356098437018154e-1)*c+6.99638348861914)*c + 21.9792616182942)*c +
      242.667955230532;
    q=((c+15.0827976304078)*c+91.1649054045149)*c + 215.058875869861;
    tserf = x*p/q;
    tserfc = 1.0-tserf;
  }
  else if (absx < 4.0) {
    c = absx;
    p = (((((((-0.136864857382717e-6)*c + 0.564195517478974)*c +
             7.21175825088309)*c + 43.1622272220567)*c +
           152.989285046940)*c + 339.320816734344)*c +
         451.918953711873)*c + 300.459261020162;
    q = ((((((c+12.7827273196294)*c + 77.0001529352295)*c +
            277.585444743988)*c + 638.980264465631)*c +
          931.354094850610)*c + 790.950925327898)*c + 300.459260956983;
    nonexperfc = (x > 0.0) ? p/q : 2.0*exp(x*x) - p/q;
    tserfc = exp(-absx*absx)*nonexperfc;
    if (x < 0.0) {
      tserfc = 2.0 - tserfc;
    }
  }
  else {
    c = 1.0/(x*x);
    p=(((0.0223192459734185*c + 0.278661308609648)*c +
        0.226956593539687)*c + 0.0494730910623251)*c + 0.00299610707703542;
    q=(((c + 1.98733201817135)*c + 1.05167510706793)*c + 0.191308926107830)*c +
      0.0106209230528468;
    c = (-c*p/q + 0.564189583547756)/absx;
    nonexperfc = (x > 0.0) ? c : 2.0*exp(x*x) - c;
    tserfc = exp(-absx*absx)*nonexperfc;
    if (x < 0.0) {
      tserfc = 2.0 - tserfc;
    }
  }

  return tserfc;
}

//-----------------------------------------------------------------------------
// PerfcFast: this is a faster way to compute the complimentary error function.
//            It assumes that x is positive, which is usually the case in the
//            sorts of calculations done for molecular dynamics simulations.
//                                                                      
// Arguments:                                                           
//   x:  the argument of the complimentary error function               
//-----------------------------------------------------------------------------
double PerfcFast(double x)
{
  double c, p, q, tserf, tserfc, invx;

  // Most likely, x will be between 0.5 and 4.0
  if (x > 0.5 && x < 4.0) {
    p = (((((((-0.136864857382717e-6)*x + 0.564195517478974)*x +
             7.21175825088309)*x + 43.1622272220567)*x +
           152.989285046940)*x + 339.320816734344)*x +
         451.918953711873)*x + 300.459261020162;
    q = ((((((x+12.7827273196294)*x + 77.0001529352295)*x +
            277.585444743988)*x + 638.980264465631)*x +
          931.354094850610)*x + 790.950925327898)*x + 300.459260956983;
    tserfc = exp(-x*x)*p/q;
    return tserfc;
  }

  // In a minority of cases, x will be less than 0.5
  if (x <= 0.5) {
    c = x*x;
    p = (((-0.356098437018154e-1)*c+6.99638348861914)*c + 21.9792616182942)*c +
      242.667955230532;
    q=((c+15.0827976304078)*c+91.1649054045149)*c + 215.058875869861;
    tserf = x*p/q;
    tserfc = 1.0-tserf;
    return tserfc;
  }

  // Otherwise, we have one more condition to evaluate about x
  if (x >= 4.0 && x < 5.65) {
    invx = 1.0/x;
    c = invx*invx;
    p=(((0.0223192459734185*c + 0.278661308609648)*c +
        0.226956593539687)*c + 0.0494730910623251)*c + 0.00299610707703542;
    q=(((c + 1.98733201817135)*c + 1.05167510706793)*c + 0.191308926107830)*c +
      0.0106209230528468;
    c = (-c*p/q + 0.564189583547756)*invx;
    tserfc = exp(-x*x)*c;
    return tserfc;
  }
  
  // If x is bigger than 5.65, the answer is easy
  return 0.0;
}

//-----------------------------------------------------------------------------
// SinFast: rapidly compute the sine function, to a precision of 3.0e-8 or
//          less.                                                    
//                                                                      
// Arguments:                                                           
//   x:  the argument of the sine function                              
//   stab: the spline interpolation table                               
//-----------------------------------------------------------------------------
double SinFast(double x, CSpln* stab)
{
  int ix;
  double y;

  x = x - TWOPI*floor(x*INV_TWOPI);
  ix = INV_SIN_DSC*x;
  y = ((stab[ix].A*x + stab[ix].B)*x + stab[ix].C)*x + stab[ix].D;

  return y;
}

//-----------------------------------------------------------------------------
// CosFast: rapidly compute the cosine function, to a precision of 3.0e-8 or
//          less.                                             
//                                                                      
// Arguments:                                                           
//   x:    the argument of the cosine function                          
//   ctab: the spline interpolation table                               
//-----------------------------------------------------------------------------
double CosFast(double x, CSpln* ctab)
{
  int ix;
  double y;

  x = x - TWOPI*floor(x*INV_TWOPI);
  ix = INV_SIN_DSC*x;
  y = ((ctab[ix].A*x + ctab[ix].B)*x + ctab[ix].C)*x + ctab[ix].D;

  return y;
}

//-----------------------------------------------------------------------------
// AsinFast: rapidly computes the arcsine function, to an accuracy of roughly
//           8.0e-10 if the argument x <= 0.76 or if x >= 0.76.  For other
//           values, the table lookup has an accuracy of 7.0e-8 or better.
//                                                                      
// Arguments:                                                           
//   x:     the argument of the arcsine function                        
//   astab: the spline interpolation table                              
//-----------------------------------------------------------------------------
double AsinFast(double x, CSpln* astab)
{
  int ix;
  double y;

  x = x - 2.0*floor(0.5*x+0.5);
  if (x >= ACOS_TAB_BUF) {
    y = ((((3.3858269499514e-3)*x - 2.1276912580328e-2)*x +
          7.00470014781691e-2)*x - 2.07657184085058e-1)*x + 1.5697148349467;
    y = HALFPI - y*sqrt(1.0-x);
    return y;
  }
  else if (x <= -ACOS_TAB_BUF) {
    x = -x;
    y = ((((3.3858269499514e-3)*x - 2.1276912580328e-2)*x +
          7.00470014781691e-2)*x - 2.07657184085058e-1)*x + 1.5697148349467;
    y = y*sqrt(1.0-x) - HALFPI;
    return y;
  }

  // If we're still here, it's OK to use the table
  x += 1.0;
  ix = INV_ASIN_DSC*x;
  y = ((astab[ix].A*x + astab[ix].B)*x + astab[ix].C)*x + astab[ix].D;

  return y;
}

//-----------------------------------------------------------------------------
// AcosFast: rapidly computes the arccosine function, to an accuracy of 
//           roughly 7.0e-8 or better.  For the most commonly needed acos
//           evaluations in molecular dynamics (near 1 or -1), this function
//           will give an accuracy of 8.0e-10 or less.         
//                                                                      
// Arguments:                                                           
//   x:     the argument of the arccosine function                      
//   actab: the spline interpolation table                              
//-----------------------------------------------------------------------------
double AcosFast(double x, CSpln* actab)
{
  int ix;
  double y;

  x = x - 2.0*floor(0.5*x+0.5);
  if (x >= ACOS_TAB_BUF) {
    y = ((((3.3858269499514e-3)*x - 2.1276912580328e-2)*x + 
	  7.00470014781691e-2)*x - 2.07657184085058e-1)*x + 1.5697148349467;
    y *= sqrt(1.0-x);
    return y;
  }
  else if (x <= -ACOS_TAB_BUF) {
    x = -x;
    y = ((((3.3858269499514e-3)*x - 2.1276912580328e-2)*x +
          7.00470014781691e-2)*x - 2.07657184085058e-1)*x + 1.5697148349467;
    y = PI-y*sqrt(1.0-x);
    return y;
  }

  // If we're still here, it's OK to use the table
  x += 1.0;
  ix = INV_ASIN_DSC*x;
  y = ((actab[ix].A*x + actab[ix].B)*x + actab[ix].C)*x + actab[ix].D;

  return y;
}

//-----------------------------------------------------------------------------
// SafeSqrt: a function for doing a "safe" square root operation, which first
//           checks that the argument to the square root function is positive.
//           Function returns zero in any other case.     
//                                                                      
// Arguments:                                                           
//   x:     the argument of the square root function
//-----------------------------------------------------------------------------
double SafeSqrt(double x)
{
  return (x > 1.0e-12) ? sqrt(x) : 0.0;
}
