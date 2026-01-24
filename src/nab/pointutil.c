#include <stdio.h>

#include "nab.h"

POINT_T	*NAB_ptcpy( POINT_T, POINT_T );
POINT_T	*NAB_ptadd( POINT_T, POINT_T, POINT_T );
POINT_T	*NAB_ptsub( POINT_T, POINT_T, POINT_T );
POINT_T	*NAB_ptscl( POINT_T, REAL_T, POINT_T );
POINT_T	*NAB_ptcrs( POINT_T, POINT_T, POINT_T ); 
REAL_T	NAB_ptdot( POINT_T, POINT_T );

POINT_T	*NAB_ptcpy( POINT_T ps, POINT_T pt )
{

	ps[ 0 ] = pt[ 0 ];
	ps[ 1 ] = pt[ 1 ];
	ps[ 2 ] = pt[ 2 ];
	return( ( POINT_T * )ps );
}

POINT_T	*NAB_ptadd( POINT_T ps, POINT_T p1, POINT_T p2 )
{

	ps[ 0 ] = p1[ 0 ] + p2[ 0 ];
	ps[ 1 ] = p1[ 1 ] + p2[ 1 ];
	ps[ 2 ] = p1[ 2 ] + p2[ 2 ];
	return( ( POINT_T * )ps );
}

POINT_T	*NAB_ptsub( POINT_T pd, POINT_T p1, POINT_T p2 )
{

	pd[ 0 ] = p1[ 0 ] - p2[ 0 ];
	pd[ 1 ] = p1[ 1 ] - p2[ 1 ];
	pd[ 2 ] = p1[ 2 ] - p2[ 2 ];
	return( ( POINT_T * )pd );
}

POINT_T	*NAB_ptscl( POINT_T ps, REAL_T a, POINT_T p1 )
{

	ps[ 0 ] = a * p1[ 0 ];
	ps[ 1 ] = a * p1[ 1 ];
	ps[ 2 ] = a * p1[ 2 ];
	return( ( POINT_T * )ps );
}

POINT_T	*NAB_ptcrs( POINT_T pc, POINT_T p1, POINT_T p2 ) 
{

	pc[ 0 ] = p1[ 1 ] * p2[ 2 ] - p2[ 1 ] * p1[ 2 ];
	pc[ 1 ] = p2[ 0 ] * p1[ 2 ] - p1[ 0 ] * p2[ 2 ];
	pc[ 2 ] = p1[ 0 ] * p2[ 1 ] - p2[ 0 ] * p1[ 1 ];
	return( ( POINT_T * )pc );
}

REAL_T	NAB_ptdot( POINT_T p1, POINT_T p2 )
{

	return( p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2] );
}
