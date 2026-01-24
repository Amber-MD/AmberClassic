#include <stdio.h>
#include <math.h>

#include "nab.h"

#define	MKAXIS(ax,x,y,z)	((ax)[0]=(x),(ax)[1]=(y),(ax)[2]=(z))

static	void	vnorm( REAL_T v[], int n );
static	void	vcross( REAL_T v1[], REAL_T v2[], REAL_T vc[] );

int	axis2frame( MOLECULE_T *m, POINT_T p1, POINT_T p2 )
{
	REAL_T	dx, dy, dz;
	REAL_T	xax[ 3 ], yax[ 3 ], zax[ 3 ];

	dx = p2[ 0 ] - p1[ 0 ];
	dy = p2[ 1 ] - p1[ 1 ];
	dz = p2[ 2 ] - p1[ 2 ];
	MKAXIS( zax, dx, dy, dz );
	vnorm( zax, 3 );
	if( dz ){
		MKAXIS( yax, 0.0, 1.0, -dy/dz );
		vnorm( yax, 3 );
		vcross( yax, zax, xax );
		vnorm( xax, 3 );
	}else if( dy ){		/* in XY plane */
		MKAXIS( yax, 1.0, -dx/dy, 0.0 );
		vnorm( yax, 3 );
		vcross( yax, zax, xax );
		vnorm( xax, 3 );
	}else if( dx ){
		MKAXIS( yax, 0.0, 1.0, 0.0 );
		MKAXIS( xax, 0.0, 0.0, -1.0 );
	}else{
		MKAXIS( xax, 1.0, 0.0, 0.0 );
		MKAXIS( yax, 0.0, 1.0, 0.0 );
		MKAXIS( zax, 0.0, 0.0, 1.0 );
	}
	m->m_frame[ 0 ][ 0 ] = p1[ 0 ];
	m->m_frame[ 0 ][ 1 ] = p1[ 1 ];
	m->m_frame[ 0 ][ 2 ] = p1[ 2 ];

	m->m_frame[ 1 ][ 0 ] = xax[ 0 ] + p1[ 0 ];
	m->m_frame[ 1 ][ 1 ] = xax[ 1 ] + p1[ 1 ];
	m->m_frame[ 1 ][ 2 ] = xax[ 2 ] + p1[ 2 ];

	m->m_frame[ 2 ][ 0 ] = yax[ 0 ] + p1[ 0 ];
	m->m_frame[ 2 ][ 1 ] = yax[ 1 ] + p1[ 1 ];
	m->m_frame[ 2 ][ 2 ] = yax[ 2 ] + p1[ 2 ];

	m->m_frame[ 3 ][ 0 ] = zax[ 0 ] + p1[ 0 ];
	m->m_frame[ 3 ][ 1 ] = zax[ 1 ] + p1[ 1 ];
	m->m_frame[ 3 ][ 2 ] = zax[ 2 ] + p1[ 2 ];

	return( 0 );
}

static	void	vnorm( REAL_T v[], int n )
{
	int	i;
	REAL_T	vn;

	for( vn = 0.0, i = 0; i < n; i++ )
		vn += v[ i ] * v[ i ];
	if( vn ){
		vn = sqrt( vn );
		for( i = 0; i < n; i++ )
			v[ i ] /= vn;
	}
}

static	void	vcross( REAL_T v1[], REAL_T v2[], REAL_T vc[] )
{

	vc[ 0 ] = v1[ 1 ] * v2[ 2 ] - v2[ 1 ] * v1[ 2 ];
	vc[ 1 ] = -( v1[ 0 ] * v2[ 2 ] - v2[ 0 ] * v1[ 2 ] );
	vc[ 2 ] = v1[ 0 ] * v2[ 1 ] - v2[ 0 ] * v1[ 1 ];
}
