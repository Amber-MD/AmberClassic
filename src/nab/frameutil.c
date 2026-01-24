#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "nab.h"
#include "matop.h"
#include "memutil.h"
#include "molutil.h"
#include "select_atoms.h"

static	POINT_T	org;
static	POINT_T	xtail;
static	POINT_T	xhead;
static	POINT_T	ytail;
static	POINT_T	yhead;

static	FRAME_T	gl_frame = {
	{ 0.0, 0.0, 0.0 },
	{ 1.0, 0.0, 0.0 },
	{ 0.0, 1.0, 0.0 },
	{ 0.0, 0.0, 1.0 } };

static	int	count_atoms( MOLECULE_T * );
static	int	count_sel_atoms( MOLECULE_T * );
static	void	get_apos( MOLECULE_T *, POINT_T [] );
static	void	get_sel_apos( MOLECULE_T *, POINT_T [] );
static	void	put_apos( MOLECULE_T *, POINT_T [] );
static	void	makeframe( int, MOLECULE_T * );
static	void	crossp( REAL_T *, REAL_T *, REAL_T * );
static	REAL_T	dotp( REAL_T *, REAL_T *, int );
static	void	normalize( REAL_T [], int );

static	void	overlay( int, int [], REAL_T [][3], REAL_T [][3], MATRIX_T );
static	void	copypts( void *, void *, int );
static	void	cenmas( int, int [], REAL_T [][3],
	REAL_T *, REAL_T *, REAL_T *, REAL_T * );
static	void	trans( int, REAL_T [][3], REAL_T, REAL_T, REAL_T );
static	void	rotate( int, REAL_T [][3], REAL_T [][3] );
static	void	lsqrot( int, int [], REAL_T [][3], REAL_T [][3], REAL_T [][3] );

int	jacobi( REAL_T **, int, REAL_T [], REAL_T **, int * );
void	eigsrt( REAL_T [], REAL_T **, int );

int	setframe( int usevec, MOLECULE_T *mol, char origin[],
	char x_tail[], char x_head[], char y_tail[], char y_head[] )
{
	int	rval = 0;

	if( (rval = setpoint( mol, origin, org )) )
		goto clean_up;
	if( (rval = setpoint( mol, x_tail, xtail )) )
		goto clean_up;
	if( (rval = setpoint( mol, x_head, xhead )) )
		goto clean_up;
	if( (rval = setpoint( mol, y_tail, ytail )) )
		goto clean_up;
	if( (rval = setpoint( mol, y_head, yhead )) )
		goto clean_up;

	makeframe( usevec, mol );

clean_up : ;

	return( rval );
}

int	setframep( int usevec, MOLECULE_T *mol, POINT_T origin,
	POINT_T x_tail, POINT_T x_head, POINT_T y_tail, POINT_T y_head )
{
	copypts( org, origin, 1 );
	copypts( xtail, x_tail, 1 );
	copypts( xhead, x_head, 1 );
	copypts( ytail, y_tail, 1 );
	copypts( yhead, y_head, 1 );
	makeframe( usevec, mol );
	return( 0 );
}

int	alignframe( MOLECULE_T *mol, MOLECULE_T *mtemplate )
{
	int	r, a, na, p, i;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap;
	POINT_T		*g_pos;
	POINT_T		*a_pos;
	int		*a_mask;
	static	MATRIX_T	mat;

	for( na = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++ )
			na += sp->s_residues[ r ]->r_natoms;
	}

	na += 4;	/* add space for the frame	*/
	if( ( g_pos = ( POINT_T * )malloc( ( na ) * sizeof( POINT_T ) ) )
		== NULL ){
		fprintf( stderr, "can't allocate new a_pos array\n" );
		exit( 1 );
	}
	if( ( a_pos = ( POINT_T * )malloc( ( na ) * sizeof( POINT_T ) ) )
		== NULL ){
		fprintf( stderr, "can't allocate new a_pos array\n" );
		exit( 1 );
	}
	if( ( a_mask = ( int * )malloc( ( na ) * sizeof( int ) ) ) == NULL ){
		fprintf( stderr, "can't allocate new a_mask array\n" );
		exit( 1 );
	}

	for( i = 0; i < 4; i++ )
		a_mask[ i ] = 1;
	for( i = 4; i < na; i++ )
		a_mask[ i ] = 0;

	copypts( a_pos, mol->m_frame, 4 );
	for( p = 4, sp = mol->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			for( a = 0; a < res->r_natoms; a++ ){
				ap = &res->r_atoms[ a ];
				for( i = 0; i < 3; i++ ){
					a_pos[ p ][ i ] = ap->a_pos[ i ];
					g_pos[ p ][ i ] = 0.0;
				}
				p++;
			}
		}
	}

	if( !mtemplate ){
		/* align molecule to global XYZ	*/
		copypts( g_pos, gl_frame, 4 );
		overlay( na, a_mask, a_pos, g_pos, mat );
	}else{
		/* align molecule to template 	*/
		copypts( g_pos, mtemplate->m_frame, 4 );
		overlay( na, a_mask, a_pos, g_pos, mat );
	}

	for( p = 4, sp = mol->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			for( a = 0; a < res->r_natoms; a++ ){
				ap = &res->r_atoms[ a ];
				for( i = 0; i < 3; i++ )
					ap->a_pos[ i ] = a_pos[ p ][ i ];
				p++;
			}
		}
	}
	
	free( a_mask );
	free( a_pos );
	free( g_pos );

	return( 0 );
}

REF_MATRIX_T	superimpose( MOLECULE_T *mol, char aexpr[],
	MOLECULE_T *r_mol, char r_aexpr[] )
{
	int	n_atoms, s_atoms;
	int	s_ratoms;
	int	a, na;
	int	*a_mask;
	POINT_T	*a_pos;
	POINT_T	*ra_pos;
	static	MATRIX_T	mat;
	int	i,j;

	select_atoms( mol, aexpr );
	n_atoms = count_atoms( mol );
	s_atoms = count_sel_atoms( mol );

	select_atoms( r_mol, r_aexpr );
	s_ratoms = count_sel_atoms( r_mol );

	if( s_atoms != s_ratoms ){
		fprintf( stderr,
			"superimpose: atom mismatch: mol(%d) & r_mol(%d)\n", 
			s_atoms, s_ratoms );
		for( i = 0; i < 4; i++ ){
			for( j = 0; j < 4; j++ )
				mat[i][j] = 0.;
		}
		return( mat );
	}

	na = s_atoms + n_atoms;
	if( ( a_mask = ( int * )malloc( na * sizeof( int ) ) )
		== NULL ){
		fprintf( stderr, "can't allocate new a_mask array\n" );
		exit( 1 );
	}
	if( ( a_pos = ( POINT_T * )malloc( na * sizeof( POINT_T ) ) ) ==
		NULL ){
		fprintf( stderr, "can't allocate new a_pos array\n" );
		exit( 1 );
	}
	if( ( ra_pos = ( POINT_T * )malloc( na * sizeof( POINT_T ) ) ) ==
		NULL ){
		fprintf( stderr, "can't allocate new ra_pos array\n" );
		exit( 1 );
	}

	for( a = 0; a < s_atoms; a++ )
		a_mask[ a ] = 1;
	for( a = s_atoms; a < na; a++ )
		a_mask[ a ] = 0;

	get_sel_apos( mol, a_pos );
	get_apos( mol, &a_pos[ s_atoms ] );
	get_sel_apos( r_mol, ra_pos );
	for( a = s_atoms; a < na; a++ ){
		ra_pos[ a ][ 0 ] = 0.0;
		ra_pos[ a ][ 1 ] = 0.0;
		ra_pos[ a ][ 2 ] = 0.0;
	}
	overlay( na, a_mask, a_pos, ra_pos, mat );
	put_apos( mol, &a_pos[ s_atoms ] );

	free( ra_pos );
	free( a_pos );
	free( a_mask );
	return( mat );
}

int	rmsd( MOLECULE_T **m1, char **aex1,
	MOLECULE_T **m2, char **aex2, REAL_T *r )
{
	int	rval = 0;
	int	a, na1, na2;
	REAL_T	dx, dy, dz, d2;
	POINT_T	*a1pos = NULL;
	POINT_T	*a2pos = NULL;

	*r = 0.0;
	select_atoms( *m1, aex1 ? *aex1 : NULL );
	if( !( na1 = count_sel_atoms( *m1 ) ) ){
		fprintf( stderr, "rms: no atoms selected in mol 1\n" );
		return( 1 );
	}
	select_atoms( *m2, aex2 ? *aex2 : NULL );
	if( !( na2 = count_sel_atoms( *m2 ) ) ){
		fprintf( stderr, "rms: no atoms selected in mol 2\n" );
		return( 1 );
	}
	if( na1 != na2 ){
		fprintf( stderr, "rms: atom mismatch m1(%d) & m2(%d)\n",
			na1, na2 );
		return( 1 );
	}

	if( ( a1pos = ( POINT_T * )malloc( na1 * sizeof( POINT_T ) ) ) ==
		NULL ){
		fprintf( stderr, "rms: can't allocate new a1pos array\n" );
		rval = 1;
		goto clean_up;
	}

	if( ( a2pos = ( POINT_T * )malloc( na2 * sizeof( POINT_T ) ) ) ==
		NULL ){
		fprintf( stderr, "rms: can't allocate new a2pos array\n" );
		rval = 1;
		goto clean_up;
	}

	get_sel_apos( *m1, a1pos );
	get_sel_apos( *m2, a2pos );

	for( d2 = 0.0, a = 0; a < na1; a++ ){
		dx = a1pos[ a ][ 0 ] - a2pos[ a ][ 0 ];
		dy = a1pos[ a ][ 1 ] - a2pos[ a ][ 1 ];
		dz = a1pos[ a ][ 2 ] - a2pos[ a ][ 2 ];
		d2 += dx * dx + dy * dy + dz * dz;
	}
	*r = sqrt( d2 / na1 );

clean_up : ;
	if( a1pos )
		free( a1pos );
	if( a2pos )
		free( a2pos );

	return( rval );
}

static	int	count_atoms( MOLECULE_T *m )
{
	STRAND_T	*sp;
	int		ta, r;

	for( ta = 0, sp = m->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++ )
			ta += sp->s_residues[ r ]->r_natoms;
	}
	return( ta );
}

static	int	count_sel_atoms( MOLECULE_T *m )
{
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap;
	int		ta, a, r;

	for( ta = 0, sp = m->m_strands; sp; sp = sp->s_next ){
		if( AT_SELECT & sp->s_attr ){
			for( r = 0; r < sp->s_nresidues; r++ ){
				res = sp->s_residues[ r ];
				if( AT_SELECT & res->r_attr ){
					for( a = 0; a < res->r_natoms; a++ ){
						ap = &res->r_atoms[ a ];
						if( AT_SELECT & ap->a_attr )
							ta++;
					}
				}
			}
		}
	}
	return( ta );
}

static	void	get_apos( MOLECULE_T *m, POINT_T a_pos[] )
{
	int		r, a, i;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap;

	for( i = 0, sp = m->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			for( a = 0; a < res->r_natoms; a++ ){
				ap = &res->r_atoms[ a ];
				a_pos[ i ][ 0 ] = ap->a_pos[ 0 ];
				a_pos[ i ][ 1 ] = ap->a_pos[ 1 ];
				a_pos[ i ][ 2 ] = ap->a_pos[ 2 ];
				i++;
			}
		}
	}
}

static	void	get_sel_apos( MOLECULE_T *m, POINT_T a_pos[] )
{
	int		r, a, i;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap;

	for( i = 0, sp = m->m_strands; sp; sp = sp->s_next ){
		if( AT_SELECT & sp->s_attr ){
			for( r = 0; r < sp->s_nresidues; r++ ){
				res = sp->s_residues[ r ];
				if( AT_SELECT & res->r_attr ){
					for( a = 0; a < res->r_natoms; a++ ){
						ap = &res->r_atoms[ a ];
						if( AT_SELECT & ap->a_attr ){
							a_pos[i][0]
								= ap->a_pos[0]; 
							a_pos[i][1]
								= ap->a_pos[1]; 
							a_pos[i][2]
								= ap->a_pos[2]; 
							i++;
						}
					}
				}
			}
		}
	}
}

static	void	put_apos( MOLECULE_T *m, POINT_T a_pos[] )
{
	int		r, a, i;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap;

	for( i = 0, sp = m->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			for( a = 0; a < res->r_natoms; a++ ){
				ap = &res->r_atoms[ a ];
				ap->a_pos[ 0 ] = a_pos[ i ][ 0 ];
				ap->a_pos[ 1 ] = a_pos[ i ][ 1 ];
				ap->a_pos[ 2 ] = a_pos[ i ][ 2 ];
				i++;
			}
		}
	}
}

static	void	makeframe( int usevec, MOLECULE_T *mol )
{
	int	i;
	REAL_T	r_xvec[ 3 ];
	REAL_T	r_yvec[ 3 ];
	REAL_T	xvec[ 3 ];
	REAL_T	yvec[ 3 ];
	REAL_T	zvec[ 3 ];

	for( i = 0; i < 3; i++ ){
		r_xvec[ i ] = xhead[ i ] - xtail[ i ];
		r_yvec[ i ] = yhead[ i ] - ytail[ i ];
	}
	crossp( r_xvec, r_yvec, zvec );
	if( usevec == 1 ){	/* use X */
		for( i = 0; i < 3; i++ )
			xvec[ i ] = r_xvec[ i ];
		crossp( zvec, xvec, yvec );
	}else{			/* use Y */
		for( i = 0; i < 3; i++ )
			yvec[ i ] = r_yvec[ i ];
		crossp( yvec, zvec, xvec );
	}
	normalize( xvec, 3 );
	normalize( yvec, 3 );
	normalize( zvec, 3 );

	for( i = 0; i < 3; i++ )
		mol->m_frame[ 0 ][ i ] = org[ i ];
	for( i = 0; i < 3; i++ )
		mol->m_frame[ 1 ][ i ] = org[ i ] + xvec[ i ];
	for( i = 0; i < 3; i++ )
		mol->m_frame[ 2 ][ i ] = org[ i ] + yvec[ i ];
	for( i = 0; i < 3; i++ )
		mol->m_frame[ 3 ][ i ] = org[ i ] + zvec[ i ];
}

static	void	crossp( REAL_T *v1, REAL_T *v2, REAL_T *vx )
{

	vx[ 0 ] =    v1[ 1 ] * v2[ 2 ] - v2[ 1 ] * v1[ 2 ];
	vx[ 1 ] = -( v1[ 0 ] * v2[ 2 ] - v2[ 0 ] * v1[ 2 ] );
	vx[ 2 ] =    v1[ 0 ] * v2[ 1 ] - v2[ 0 ] * v1[ 1 ];
}

static	REAL_T	dotp( REAL_T *v1, REAL_T *v2, int vlen )
{
	int	i;
	REAL_T	dp;

	for( dp = 0.0, i = 0; i < vlen; i++ )
		dp += v1[ i ] * v2[ i ];
	return( dp );
}

static	void	normalize( REAL_T v[], int vlen )
{
	REAL_T	vmag;
	int	i;

	vmag = sqrt( dotp( v, v, vlen ) );
	if( vmag == 0 )
		return;
	for( i = 0; i < vlen; i++ )
		v[ i ] /= vmag;
}

static	void	overlay( int n_atoms, int a_mask[],
	REAL_T a_pos[][ 3 ], REAL_T a_pos0[][ 3 ], MATRIX_T mat )
{
	REAL_T	mtot0, mtot;
	REAL_T	xcen0, ycen0, zcen0;
	REAL_T	xcen, ycen, zcen;
	REAL_T	rotmat[ 3 ][ 3 ];
	int	i, j;
	MATRIX_T	tmat1, rmat, tmat2, tmp1;

	cenmas( n_atoms, a_mask, a_pos, &xcen, &ycen, &zcen, &mtot );
	trans( n_atoms, a_pos, xcen, ycen, zcen );
	cenmas( n_atoms, a_mask, a_pos0, &xcen0, &ycen0, &zcen0, &mtot0 );
	trans( n_atoms, a_pos0, xcen0, ycen0, zcen0 );
	NAB_matcpy( tmat1, newtransform( -xcen, -ycen, -zcen, 0., 0., 0. ) );  

	lsqrot( n_atoms, a_mask, a_pos, a_pos0, rotmat );
	rotate( n_atoms, a_pos, rotmat );

	trans( n_atoms, a_pos, -xcen0, -ycen0, -zcen0 );
	trans( n_atoms, a_pos0, -xcen0, -ycen0, -zcen0 );

	for( i = 0; i < 3; i++ ){
		for( j = 0; j < 3; j++ )
			rmat[j][i] = rotmat[i][j];
		rmat[i][3] = 0.0;
	}
	rmat[3][0] = 0.0;
	rmat[3][1] = 0.0;
	rmat[3][2] = 0.0;
	rmat[3][3] = 1.0;

	NAB_matcpy( tmat2,
		newtransform( xcen0, ycen0, zcen0, 0., 0., 0. ) );  

	NAB_matcpy( tmp1, MAT_concat( tmat1, rmat ) );
	NAB_matcpy( mat, MAT_concat( tmp1, tmat2 ) );
}

static	void	copypts( void *dstpts, void *srcpts, int n )
{

	memcpy( dstpts, srcpts, n * sizeof( POINT_T ) );
}

static	void	cenmas( int n_atom, int a_mask[], REAL_T a_pos[][3],
	REAL_T *xcen, REAL_T *ycen, REAL_T *zcen, REAL_T *mtot )
{
	int	i;

	*mtot = 0;
	*xcen = *ycen = *zcen = 0.0;
	for( i = 0; i < n_atom; i++ ){
		*mtot += a_mask[ i ];
		*xcen += a_mask[ i ] * a_pos[ i ][ 0 ];
		*ycen += a_mask[ i ] * a_pos[ i ][ 1 ];
		*zcen += a_mask[ i ] * a_pos[ i ][ 2 ];
	}
	if( *mtot != 0 ){
		*xcen /= *mtot;
		*ycen /= *mtot;
		*zcen /= *mtot;
	}
}

static	void	trans( int n_atom, REAL_T a_pos[][3],
	REAL_T xcen, REAL_T ycen, REAL_T zcen )
{
	int	i;

	for( i = 0; i < n_atom; i++ ){
		a_pos[ i ][ 0 ] -= xcen;
		a_pos[ i ][ 1 ] -= ycen;
		a_pos[ i ][ 2 ] -= zcen;
	}
}

static	void	rotate( int n_atom, REAL_T a_pos[][3], REAL_T rotmat[][3] )
{
	int	i, j;
	REAL_T	xrot, yrot, zrot;

	for( i = 0; i < n_atom; i++ ){
		xrot = yrot = zrot = 0.0;
		for( j = 0; j < 3; j++ ){
			xrot += rotmat[ 0 ][ j ] * a_pos[ i ][ j ];
			yrot += rotmat[ 1 ][ j ] * a_pos[ i ][ j ];
			zrot += rotmat[ 2 ][ j ] * a_pos[ i ][ j ];
		}
		a_pos[ i ][ 0 ] = xrot;
		a_pos[ i ][ 1 ] = yrot;
		a_pos[ i ][ 2 ] = zrot;
	}
}

static	void	lsqrot( int n_atom, int a_mask[],
	REAL_T a_pos[][3], REAL_T a_pos0[][3], REAL_T rotmat[][3] )
{
	int	i, j, k;
	int	nrot;
	REAL_T	trm;
	REAL_T	sum;
	REAL_T	m[ 3 ][ 3 ];
	REAL_T	p[ 10 ];
	REAL_T	rhomat[ 4 ][ 4 ] ;
	REAL_T	rhovec[ 4 ];
	static	REAL_T	**a = NULL;
	static	REAL_T	**v = NULL;
	static	REAL_T	*d = NULL;

	if( a == NULL ){
		a = matrix( 1, 4, 1, 4 );
		v = matrix( 1, 4, 1, 4 );
		d = vector( 1, 4 );
	}

	for( i = 0; i < 3; i ++ ){
		for( j = 0; j < 3; j++ )
			m[ i ][ j ] = 0.0;
	}

	for( i = 0; i < n_atom; i++ ){
		for( j = 0; j < 3; j++ ){
			for( k = 0; k < 3; k++ ){
				m[j][k] += a_mask[i]*a_pos[i][j]*a_pos0[i][k];
			}
		}
	}

	trm =  m[0][0] + m[1][1] + m[2][2];
	p[0] = m[0][0] + m[0][0] - 2.0 * trm;
	p[1] = m[0][1] + m[1][0];
	p[2] = m[1][1] + m[1][1] - 2.0 * trm;
	p[3] = m[0][2] + m[2][0];
	p[4] = m[1][2] + m[2][1];
	p[5] = m[2][2] + m[2][2] - 2.0 * trm;

	p[6] = m[1][2] - m[2][1];
	p[7] = m[2][0] - m[0][2];
	p[8] = m[0][1] - m[1][0];

	p[9] = 0.0;

	/* pack p[] into a[][]*/             /****************/
	a[ 1 ][ 1 ] =               p[ 0 ];  /* 1  2  4  7	*/
	a[ 1 ][ 2 ] = a[ 2 ][ 1 ] = p[ 1 ];  /*    3  5  8	*/
	a[ 2 ][ 2 ] = a[ 2 ][ 2 ] = p[ 2 ];  /*       6  9   */
	a[ 1 ][ 3 ] = a[ 3 ][ 1 ] = p[ 3 ];  /*         10	*/
	a[ 2 ][ 3 ] = a[ 3 ][ 2 ] = p[ 4 ];  /****************/
	a[ 3 ][ 3 ] = a[ 3 ][ 3 ] = p[ 5 ];
	a[ 1 ][ 4 ] = a[ 4 ][ 1 ] = p[ 6 ];
	a[ 2 ][ 4 ] = a[ 4 ][ 2 ] = p[ 7 ];
	a[ 3 ][ 4 ] = a[ 4 ][ 3 ] = p[ 8 ];
	a[ 4 ][ 4 ] = a[ 4 ][ 4 ] = p[ 9 ];

/*
	for( i = 1; i <= 4; i++ ){
		fprintf( stderr, "a[%d][1:4] ", i );
		for( j = 1; j <= 4; j++ )
			fprintf( stderr, "%12.4f", a[i][j] );
		putc( '\n', stderr );
	}
	putc( '\n', stderr );
*/

	jacobi( a, 4, d, v, &nrot );
	eigsrt( d, v, 4 );

/*
	for( i = 1; i <= 4; i++ ){
		fprintf( stderr, "v[%d][1:4] ", i );
		for( j = 1; j <= 4; j++ )
			fprintf( stderr, "%12.4f", v[i][j] );
		putc( '\n', stderr );
	}
	putc( '\n', stderr );
	for( i = 1; i <= 4; i++ ){
		fprintf( stderr, "d[%d]    = %12.4f\n", i, d[i] );
	}
*/

	for( i = 0; i < 4; i++ )
		rhovec[ i ] = v[ i + 1 ][ 1 ];

	for( i = 0; i < 3; i++ ){
		j = ( i + 1 ) % 3;
		k = ( j + 1 ) % 3;
		rhomat[ i ][ i ] =  rhovec[ 3 ];
		rhomat[ i ][ j ] = -rhovec[ k ];
		rhomat[ i ][ k ] =  rhovec[ j ];
		rhomat[ i ][ 3 ] =  rhovec[ i ];
		rhomat[ 3 ][ i ] = -rhovec[ i ];
	}
	rhomat[ 3 ][ 3 ] = rhovec[ 3 ];
	for( i = 0; i < 3; i++ ){
		for( j = 0; j < 3; j++ ){
			sum = -rhomat[ i ][ 3 ] * rhomat[ 3 ][ j ];
			for( k = 0; k < 3; k++ ){
				sum += rhomat[ i ][ k ] * rhomat[ k ][ j ];
			}	
			rotmat[ i ][ j ] = sum;
		}
	}
}
