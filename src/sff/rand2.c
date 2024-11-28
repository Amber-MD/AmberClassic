#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "sff.h"

#define	IM1	2147483563
#define	IM2	2147483399
#define	AM	( 1.0 / IM1 )
#define	IMM1	( IM1 - 1 )
#define	IA1	40014
#define	IA2	40692
#define	IQ1	53668
#define	IQ2	52774
#define	IR1	12211
#define	IR2	3791
#define	NTAB	32
#define	NDIV	( 1 + IMM1 / NTAB )
#define	EPS	1.2e-13
#define	RNMX	( 1.0 - EPS )

static	int	seed2 = 0;
static  int seed3 = -1;   /* surrogate for seed  */
static	int	iy = 0;
static	int	iv[ NTAB ];

/*
     Get a pseudo-random number in the range 0 -> 1.  If "seed" is
     negative, reset the seeds so that the same sequence can be
     re-created.

     This is based on the MRG32k3a RNG from Pierre L'Ecuyer:
     Random Number Generation,
     In: Handbook on Simulation,
     Jerry Banks, ed. Wiley, New York, 1997

     A reference implementation (courtesy of M.L. Dodson) is appended at
     the end of this file.

     Note that the "driver" routine, rand2(), which is provided at the 
     bottom of this file, passes the static variable "seed3" to rand2a().
     This means that the state of the system is stored internally.

     If the argument is negative, the sequence is reset to a new starting
     position, and the first pseudo-random number of that sequence is 
     returned.

*/

static
REAL_T	rand2a( int *seed )
{
	int		j, k;
	REAL_T	temp;

	if( *seed <= 0 ){
		if( -*seed < 1 )
			*seed = 1;
		else
			*seed = -*seed;
		seed2 = *seed;
		for( j = NTAB + 7; j >= 0; j-- ){
			k = *seed / IQ1;
			*seed = IA1 * ( *seed - k * IQ1 ) - k * IR1;
			if( *seed < 0 )
				*seed += IM1;
			if( j < NTAB )
				iv[ j ] = *seed;
		}
		iy = iv[ 0 ];
	}
	k = *seed / IQ1;
	*seed = IA1 * ( *seed - k * IQ1 ) - k * IR1;
	if( *seed < 0 )
		*seed += IM1;
	k = seed2 / IQ2;
	seed2 = IA2 * ( seed2 - k * IQ2 ) - k * IR2;
	if( seed2 < 0 )
		seed2 += IM2; 
	j = iy / NDIV;
	iy = iv[ j ] - seed2;
	iv[ j ] = *seed;
	if( iy < 1 )
		iy += IMM1;
	if( ( temp = AM * iy ) > RNMX )
		return( RNMX );
	else
		return( temp );
}

/*
   Generate a pseudo-random sequence of numbers with a given mean
   and standard deviation.  Use the Box & Mueller method, but only
   use the first value, since the two values are correlated.
*/

static
REAL_T gaussa( REAL_T *mean, REAL_T *sd, int *seed )
{
	REAL_T fac,gdev1,rsq,s1,s2;

		do {
			s1 = 2.*rand2a(seed) - 1.;
			s2 = 2.*rand2a(seed) - 1.;
			rsq = s1*s1 + s2*s2;
		} while ( rsq >= 1. || rsq == 0.0 );
		fac = sqrt(-2.*log(rsq)/rsq);
		gdev1 = s1*fac;

		return( *sd*gdev1 + *mean );

}

/*   
    Driver routine for randa(), so that state information is kept within
    this routine.  Same for gaussa().
*/

REAL_T  rand2( void )
{

	return rand2a( &seed3 );
}

REAL_T gauss2( REAL_T *mean, REAL_T *sd )
{
	return gaussa( mean, sd, &seed3 );
}

/*
    Set the seed with a given value
*/

int   setseed( int *seed4 )
{
	if( *seed4 >= 0 ){
		fprintf( stderr, "argument to setseed must be negative!\n" );
		return 1;
	} else {
		seed3 = *seed4;
		return 0;
	}
}

/*
    Randomize the seed using gettimeofday(); return the new seed, so that the run
    could be repeated if desired.
*/

int   rseed()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	seed3 = -((int)(tv.tv_sec ^ tv.tv_usec));
	return seed3;
}

#if 0
/***********************************************************************

   uniform [0,1] random number generator
   developed by Pierre Lecuyer based on a clever
   and tested combination of two linear congruential
   sequences

   s1 and s2 are the seeds (positive integers)

   Pierre Lecuyer,
     Random Number Generation,
     In: Handbook on Simulation,
     Jerry Banks, ed. Wiley, New York, 1997

***********************************************************************/

static long s1 = 55555;
static long s2 = 99999;


int lecuyer()
{
        int k;

	/**  s1 = 40014 * s1 % 2147483563  **/
        k  = s1 / 53668;
        s1 = 40014*(s1%53668)-k*12211;
        if (s1 < 0) s1 += 2147483563;

	/**  s2 = 40692 * s2 % 2147483399  **/
        k  = s2 / 52774;
        s2 =40692*(s2%52774)-k*3791;
        if (s2 < 0) s2 += 2147483399;

        /**  combine s1 and s2  **/
        k = (s1 - 2147483563) + s2;
        if (k < 1) k += 2147483562;

        return k;
}

double uni()
{
        const double factor = 1.0/2147483399.0;
        int k;

	/**  s1 = 40014 * s1 % 2147483563  **/
        k  = s1 / 53668;
        s1 = 40014*(s1%53668)-k*12211;
        if (s1 < 0) s1 += 2147483563;

	/**  s2 = 40692 * s2 % 2147483399  **/
        k  = s2 / 52774;
        s2 =40692*(s2%52774)-k*3791;
        if (s2 < 0) s2 += 2147483399;

        /**  combine s1 and s2  **/
        k = (s1 - 2147483563) + s2;
        if (k < 1) k += 2147483562;

        return(((double)(k))*factor);
}

#endif
