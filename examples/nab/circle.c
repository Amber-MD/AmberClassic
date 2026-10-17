#include <stdio.h>
#include <math.h>

#define	PI	3.14159265358979323844
#define	EPS	1e-4

float	vmag();
float	vangle();
float	vdot();

float	**matrix();
float	*vector();
int	*ivector();

float	**mat = NULL;
float	*b = NULL;
int	*index = NULL;

int	circle( p1, p2, p3, pc )
float	p1[];
float	p2[];
float	p3[];
float	pc[];
{
	float	v12[ 3 ], v13[ 3 ], v23[ 3 ], vn[ 3 ], d;
	float	p12m[ 3 ], p13m[ 3 ];
	float	a1, a2, a3;

	/* allocate stuff for nr		*/
	if( !mat ){
		mat = matrix( 1, 3, 1, 3 );
		b = vector( 1, 3 );
		index = ivector( 1, 3 );
	}

	/* check if any 2 pts are coincident 	*/ 

	mat[ 1 ][ 1 ] = v12[ 0 ] = p2[ 0 ] - p1[ 0 ];
	mat[ 1 ][ 2 ] = v12[ 1 ] = p2[ 1 ] - p1[ 1 ];
	mat[ 1 ][ 3 ] = v12[ 2 ] = p2[ 2 ] - p1[ 2 ];
	if( vmag( v12, 3 ) < EPS ){
		fprintf( stderr, "circle: p1, p2 are too close\n" );
		return( 1 );
	}
	mat[ 2 ][ 1 ] = v13[ 0 ] = p3[ 0 ] - p1[ 0 ];
	mat[ 2 ][ 2 ] = v13[ 1 ] = p3[ 1 ] - p1[ 1 ];
	mat[ 2 ][ 3 ] = v13[ 2 ] = p3[ 2 ] - p1[ 2 ];
	if( vmag( v13, 3 ) < EPS ){
		fprintf( stderr, "circle: p1, p3 are too close\n" );
		return( 1 );
	}
	v23[ 0 ] = p3[ 0 ] - p2[ 0 ];
	v23[ 1 ] = p3[ 1 ] - p2[ 1 ];
	v23[ 2 ] = p3[ 2 ] - p2[ 2 ];
	if( vmag( v23, 3 ) < EPS ){
		fprintf( stderr, "circle: p2, p3 are too close\n" );
		return( 1 );
	}

	/* check for colinearity		*/
	if( ( a1 = fabs( vangle( v12, v13, 3 ) ) ) < EPS || a1 > PI - EPS ){
		fprintf( stderr, "circle: a1: p1, p2, p3 are colinear\n" );
		return( 1 );
	}
	if( ( a2 = fabs( vangle( v12, v23, 3 ) ) ) < EPS || a2 > PI - EPS ){
		fprintf( stderr, "circle: a2: p1, p2, p3 are colinear\n" );
		return( 1 );
	}
	if( ( a3 = fabs( vangle( v13, v23, 3 ) ) ) < EPS || a3 > PI - EPS ){
		fprintf( stderr, "circle: a3: p1, p2, p3 are colinear\n" );
		return( 1 );
	}

	p12m[ 0 ] = p1[ 0 ] + 0.5 * mat[ 1 ][ 1 ];
	p12m[ 1 ] = p1[ 1 ] + 0.5 * mat[ 1 ][ 2 ];
	p12m[ 2 ] = p1[ 2 ] + 0.5 * mat[ 1 ][ 3 ];
	p13m[ 0 ] = p1[ 0 ] + 0.5 * mat[ 2 ][ 1 ];
	p13m[ 1 ] = p1[ 1 ] + 0.5 * mat[ 2 ][ 2 ];
	p13m[ 2 ] = p1[ 2 ] + 0.5 * mat[ 2 ][ 3 ];
	vcross( v12, v13, vn );
	vnorm( vn, 3 );
	mat[ 3 ][ 1 ] = vn[ 0 ];
	mat[ 3 ][ 2 ] = vn[ 1 ];
	mat[ 3 ][ 3 ] = vn[ 2 ];

	b[ 1 ] = vdot( v12, p12m, 3 );
	b[ 2 ] = vdot( v13,  p13m, 3 );
	b[ 3 ] = vdot( vn, p1, 3 );

	ludcmp( mat, 3, index, &d ); 
	lubksb( mat, 3, index, b );

	pc[ 0 ] = b[ 1 ];
	pc[ 1 ] = b[ 2 ];
	pc[ 2 ] = b[ 3 ];
}
