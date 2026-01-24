#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "Constants.h"
#include "Random.h"

#include "TrajectoryDS.h"

//-----------------------------------------------------------------------------
// InitPRNG: initialize the pseudo-random number generator.             
//                                                                      
// Arguments:                                                           
//   tj:      trajectory control information, contains user-specified random
//            seed and process rank                              
//-----------------------------------------------------------------------------
void InitPRNG(trajcon *tj)
{
  int i;

  tj->rndcon = -1;
  ran2(&tj->rndcon);
  tj->rndcon = SetPRNG(tj);

  // Run forward 100 numbers to ensure that
  // separate threads become de-correlated 
  if (tj->nthreads > 1 || tj->SyncRNG == 1) {
    for (i = 0; i < 100; i++) {
      ran2(&tj->rndcon);
    }
  }
}

//-----------------------------------------------------------------------------
// SetPRNG: this sets the pseudo-random number generator based on a seed.  If
//          the seed is less than zero, then the seed is taken from the wall
//          time and returned.                            
//                                                                      
// Arguments:                                                           
//   tj:      trajectory control information, contains user-specified random
//            seed and process rank                              
//-----------------------------------------------------------------------------
long int SetPRNG(trajcon *tj)
{
  struct timeval trng;

  if (tj->igseed < 0) {
    gettimeofday(&trng, NULL);
    tj->igseed = 1.0e6*(trng.tv_sec % 1000) + trng.tv_usec;
  }
  if (tj->SyncRNG == 0) {
    tj->igseed += tj->tid;
  }

  return tj->igseed;
}

//-----------------------------------------------------------------------------
// RAN2: function for returning a single random number from a uniform   
//       distribution in the range (0, 1), exclusive of the endpoints.  
//                                                                      
// Arguments:                                                           
//   nwin:  pointer to a long unsigned integer, the window into the 3-digit
//          random number state vector                         
//-----------------------------------------------------------------------------
double ran2(long *nwin) 
{ 
  int j;
  long k, nwinlocal;
  static long nwin2=123456789; 
  static long iy=0;
  static long iv[NTAB]; 
  double temp;

  // Local value for the modifiable input variable
  nwinlocal = *nwin;

  if (nwinlocal <= 0) { 
    nwinlocal = (-nwinlocal < 1) ? 1 : -nwinlocal;
    nwin2 = nwinlocal;
    for (j = NTAB+7; j >= 0 ; j--) {
      k = nwinlocal/IQ1;
      nwinlocal = IA1*(nwinlocal - k*IQ1) - k*IR1; 
      if (nwinlocal < 0) {
	nwinlocal += IM1;
      }
      if (j < NTAB) {
	iv[j] = nwinlocal;
      }
    }
    iy = iv[0];
  }
  k = nwinlocal/IQ1;
  nwinlocal = IA1*(nwinlocal - k*IQ1) - k*IR1; 
  if (nwinlocal < 0) {
    nwinlocal += IM1;
  }
  k = nwin2/IQ2;
  nwin2 = IA2*(nwin2 - k*IQ2) - k*IR2; 
  if (nwin2 < 0) {
    nwin2 += IM2;
  }
  j = iy/NDIV;
  iy = iv[j] - nwin2; 
  iv[j] = nwinlocal;

  // Copy local value back to modifiable input variable
  *nwin = nwinlocal;

  if (iy < 1) iy += IMM1; 
  temp = AM*iy;
  if (temp > RNMX) {
    return RNMX;
  }
  else {
    return temp;
  }
}

//-----------------------------------------------------------------------------
// GaussBoxMuller: obtain a Gaussian distribution of random numbers from a
//                 uniform one.  The mean of the Gaussian is 0.0 and its width 
//                 is 1.0.  This function returns only one number; although it
//                 can generate two they would be correlated with each other.
//
// Arguments:                                                           
//   counter:  pointer to a long unsigned (typically 32-bit) integer    
//-----------------------------------------------------------------------------
double GaussBoxMuller(long *counter)
{
  double x1, x2, y;

  x1 = sqrt(-2.0*log(ran2(counter)));
  x2 = sin(TWOPI*ran2(counter));
  y = x1 * x2;

  return y;
}
