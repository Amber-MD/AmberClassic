#include <math.h>
#include <stdio.h>
#include "pmeDirect.h"
#include "SpecialMath.h"

//-----------------------------------------------------------------------------
// EwaldCoefficient: returns the Ewald coefficient.                     
//                                                                      
// Arguments:                                                           
//   Lcut:  the direct space (electrostatic) cutoff                     
//   Dtol:  the direct sum tolerance                                    
//-----------------------------------------------------------------------------
double EwaldCoefficient(double Lcut, double Dtol)
{
  int h, i;
  double calpha, dela, maxerr, amin, amax, testa, da;

  calpha = 1.0;
  maxerr = fabs(perfc(calpha*Lcut)/Lcut - Dtol);
  amin = 0.0;
  amax = 100.0;
  for (h = 0; h < 15; h++) {
    dela = 0.005*(amax-amin);
    for (i = 0; i <= 200; i++) {
      da = amin + dela*i;
      testa = fabs(perfc(da*Lcut)/Lcut - Dtol);
      if (testa < maxerr) {
        maxerr = testa;
	    calpha = da;
      }
    }
    amin = calpha - 5.0*dela;;
    amin = (amin < 0.0) ? 0.0 : amin;
    amax = calpha + 5.0*dela;
  }

  return calpha;
}
