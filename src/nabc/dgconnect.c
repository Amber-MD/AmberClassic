#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nabc.h"
#include "errormsg.h"
#include "nrutil.h"

#define	CLOSEST		2.0
#define	FURTHEST	10000.0

static	int	mk_boundsmat( char [], MOLECULE_T *, int );
static	int	set_bounds_ij( RESIDUE_T *, RESIDUE_T *, int , int *, int * );
static	void	mk_intrares_dist( float **, int, RESIDUE_T * );
static	void	mk_interres_dist( float **, int, int,
	RESIDUE_T *, int, RESIDUE_T *, int );
static	int	is_12_connected( RESIDUE_T *, int, int );
static	int	is_13_connected( RESIDUE_T *, int, int, int * );
static	int	is_14_connected( RESIDUE_T *, int, int, int *, int * );
static	float	dist( RESIDUE_T *, int, int );
static	void	get_14_bounds( RESIDUE_T *, int, int, int, int,
	float *, float * );
static	void	get_measured_bounds( RESIDUE_T *, int, int, float *, float * );
static	int	writematrix( char [], float **, int );

int	dgconnect( char bfname[], MOLECULE_T *mol, int mlen )
{

fprintf( stderr, "dgconnect: mlen = %d\n", mlen );

	if( mk_boundsmat( bfname, mol, mlen ) )
		return( 1 );
	return( 0 );
}

static	int	mk_boundsmat( char bfname[], MOLECULE_T *mol, int mlen )
{
	int	i, j, ri, rj, rjlim;
	int	ta, tai, taj, tres, tresi, tresj; 
	int	adj, xyz2b;
	STRAND_T	*sp;
	RESIDUE_T	*resi, *resj;
	RESIDUE_T	**rindex;
	float		**bmat;
	int	*roff;

	for( ta = 0, tres = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		tres += sp->s_nresidues;
		for( ri = 0; ri < sp->s_nresidues; ri++ ){
			resi = sp->s_residues[ ri ];
			ta += resi->r_natoms;
		}
	}

	if( ( roff = ( int * )malloc( ta * sizeof( int ) ) ) == NULL ){
		rt_errormsg_s( TRUE, E_NOMEM_FOR_S, "roff[]" );
		return( 1 );
	}

	if( ( rindex = ( RESIDUE_T ** )malloc( tres * sizeof( RESIDUE_T * ) ) )
		== NULL ){
		rt_errormsg_s( TRUE, E_NOMEM_FOR_S, "rindex[]" );
		return( 1 );
	}

	for( ta = 0, tres = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		for( ri = 0; ri < sp->s_nresidues; ri++, tres++ ){
			roff[ tres ] = ta;
			resi = sp->s_residues[ ri ];
			rindex[ tres ] = resi;
			ta += resi->r_natoms;
fprintf( stderr, "roff[%2d] (%s) = %3d\n", tres + 1, resi->r_name, ta );
		}
	}

fprintf( stderr, "dgc: mkm, %d total atoms %d total residues\n", ta, tres );

	if( ( bmat = matrix( 0, ta - 1, 0, ta - 1 ) ) == NULL )
		return( 1 );

	for( i = 0; i < ta; i++ ){
		for( j = i + 1; j < ta; j++ ){
			bmat[ j ][ i ] = CLOSEST;
			bmat[ i ][ j ] = FURTHEST;
		}
	}
	for( i = 0; i < ta; i++ )
		bmat[ i ][ i ] = 0.0;

	for( ri = 0; ri < tres; ri++ ){
		resi = rindex[ ri ];
		mk_intrares_dist( bmat, roff[ ri ], resi );
		for( rj = ri + 1; rj < tres; rj++ ){
			resj = rindex[ rj ];
			if( set_bounds_ij( resi, resj, mlen, &xyz2b, &adj ) ){

fprintf( stderr, "ri = %d, roff[i] = %d, rj = %d, roff[j] = %d, %s, %s\n",
	ri, roff[ ri ], rj, roff[ rj ],
	xyz2b ? "xyz2b" : "NOTxyz2b", adj ? "Adj" : "NOTAdj" );

			mk_interres_dist( bmat, xyz2b, adj,
				resi, roff[ ri ], resj, roff[ rj ] );
			}
		}
	}

	writematrix( bfname, bmat, ta );

	free_matrix( bmat, 0, ta - 1, 0, ta - 1 );

	free( rindex );

	free( roff );

	return( 0 );
}

static	int	set_bounds_ij( RESIDUE_T *resi, RESIDUE_T *resj,
	int mlen, int *xyz2b, int *adj )
{
	int	ri, rj, nrj;
	STRAND_T 	*spi, *spj;

	*xyz2b = *adj = FALSE;
	spi = resi->r_strand;
	for( ri = 0; ri < spi->s_nresidues; ri++ )
		if( spi->s_residues[ ri ] == resi )
			break;
	spj = resj->r_strand;
	nrj = spj->s_nresidues;
	for( rj = 0; rj < spj->s_nresidues; rj++ )
		if( spj->s_residues[ rj ] == resj )
			break;
	if( spi == spj ){
		if( rj == ri + 1 )
			*adj = TRUE;
		else if( rj <= ri + mlen )
			*xyz2b = TRUE;
	}else if( ri == nrj - 1 - rj  )
		*xyz2b = TRUE;	/* paired, KLUDGE!!!	*/

	return( *xyz2b || *adj );
}

static	void	mk_intrares_dist( float **bmat, int roff, RESIDUE_T *res )
{
	int	ai, aj, a2, a3;
	int	aki, akj;
	ATOM_T	*api, *apj;
	float	ub, lb;

fprintf( stderr, "res %s, roff = %d\n", res->r_name, roff );
	for( ai = 0; ai < res->r_natoms - 1; ai++ ){
		api = &res->r_atoms[ ai ];
		aki = api->a_kind;
		for( aj = ai + 1; aj < res->r_natoms; aj++ ){
			apj = &res->r_atoms[ aj ];
			akj = apj->a_kind;
			if( is_12_connected( res, ai, aj ) ){
				lb = ub = dist( res, ai, aj );
			}else if( is_13_connected( res, ai, aj, &a2 ) ){
				lb = ub = dist( res, ai, aj );
			}else if( is_14_connected( res, ai, aj, &a2, &a3 ) ){
				get_14_bounds( res, ai, a2, a3, aj, &lb, &ub );
			}else{
				get_measured_bounds( res, ai, aj, &lb, &ub );
			}
		bmat[ ai + roff ][ aj + roff ] = ub;
		bmat[ aj + roff ][ ai + roff ] = lb;
/*
		if( ub < FURTHEST )
			fprintf( stderr, "d(%4s:%d,%4s:%d) = %10.3f\n",
				api->a_name, ai, apj->a_name, aj, ub );
*/
		}
	} 
}

static	void	mk_interres_dist( float **bmat, int xyz2b, int adj,
	RESIDUE_T *res1, int r1off, RESIDUE_T *res2, int r2off )
{
	int	a1, a2;
	ATOM_T	*ap1, *ap2;
	float	d, dx, dy, dz;
	float	lb, ub;

fprintf( stderr, "res1 = %s, roff = %d\n", res1->r_name, r1off );
fprintf( stderr, "res2 = %s, roff = %d\n", res2->r_name, r2off );
	for( a1 = 0; a1 < res1->r_natoms; a1++ ){
		ap1 = &res1->r_atoms[ a1 ];
		for( a2 = 0; a2 < res2->r_natoms; a2++ ){
			ap2 = &res2->r_atoms[ a2 ];
			if( xyz2b &&
				ap1->a_kind == AK_BASE &&
				ap2->a_kind == AK_BASE ){
				dx = ap1->a_pos[ 0 ] - ap2->a_pos[ 0 ];
				dy = ap1->a_pos[ 1 ] - ap2->a_pos[ 1 ];
				dz = ap1->a_pos[ 2 ] - ap2->a_pos[ 2 ];
				d = sqrt( dx * dx + dy * dy + dz * dz );
				lb = 0.95 * d;
				ub = 1.05 * d;
				bmat[ r2off + a2 ][ r1off + a1 ] = lb;
				bmat[ r1off + a1 ][ r2off + a2 ] = ub;
/*
				fprintf( stderr,
					"d(%4s:%d,%4s:%d) = %10.3f:%10.3f\n",
					ap1->a_name, a1, ap2->a_name, a2,
					lb, ub );
*/
			}else if( adj ){
				get_ij_bounds( res1, a1, res2, a2, &lb, &ub );
				bmat[ r2off + a2 ][ r1off + a1 ] =  lb;
				bmat[ r1off + a1 ][ r2off + a2 ] =  ub;
/*
				fprintf( stderr,
					"d(%4s:%d,%4s:%d) = %10.3f:%10.3f\n",
					ap1->a_name, a1, ap2->a_name, a2,
					lb, ub );
*/
			}
		}
	}
}

static	int	is_12_connected( RESIDUE_T *res, int ai, int aj )
{
	int	a;
	ATOM_T	*api;

	api = &res->r_atoms[ ai ];
	for( a = 0; a < api->a_nconnect; a++ ){
		if( api->a_connect[ a ] == aj )
			return( TRUE );
	}
	return( FALSE );
}

static	int	is_13_connected( RESIDUE_T *res, int ai, int aj, int *a2 )
{
	int	a;
	ATOM_T	*api;

	api = &res->r_atoms[ ai ];
	for( a = 0; a < api->a_nconnect; a++ ){
		if( is_12_connected( res, api->a_connect[ a ], aj ) ){
			*a2 = api->a_connect[ a ];
			return( TRUE );
		}
	}
	return( FALSE );
}

static	int	is_14_connected( RESIDUE_T *res, int ai, int aj,
	int *a2, int *a3 )
{
	int	a;
	ATOM_T	*api;

	api = &res->r_atoms[ ai ];
	for( a = 0; a < api->a_nconnect; a++ ){
		if( is_13_connected( res, api->a_connect[ a ], aj, a3 ) ){
			*a2 = api->a_connect[ a ];
			return( TRUE );
		}
	}
	return( FALSE );
}

static	float	dist( RESIDUE_T *res, int ai, int aj )
{
	float	dx, dy, dz;
	ATOM_T	*api, *apj;

	api = &res->r_atoms[ ai ];
	apj = &res->r_atoms[ aj ];
	dx = api->a_pos[ 0 ] - apj->a_pos[ 0 ];
	dy = api->a_pos[ 1 ] - apj->a_pos[ 1 ];
	dz = api->a_pos[ 2 ] - apj->a_pos[ 2 ];
	return( sqrt( dx * dx + dy * dy + dz * dz ) );
}

static	void	get_14_bounds( RESIDUE_T *res, int a1, int a2, int a3, int a4,
	float *lb, float *ub )
{
	float	d12, d23, d13, d34, d24, d14;
	float	x1, y1, x2, y2, x3, y3, x4, y4;
	float	dx, dy, dz;

	d12 = dist( res, a1, a2 );
	d23 = dist( res, a2, a3 );
	d13 = dist( res, a1, a3 );
	d34 = dist( res, a3, a4 );
	d24 = dist( res, a2, a4 );

	x2 = y2 = y3 = 0.0;
	x3 = d23;
	x1 = - ( ( d13 * d13 - d23 * d23 - d12 * d12 ) / ( 2.0 * d23 ) );
	y1 = sqrt( d12 * d12 - x1 * x1 );

	x4 = ( d24 * d24 - d23 * d23 - d34 * d34 ) / ( 2.0 * d23 );
	y4 = sqrt( d34 * d34 - x4 * x4 );
	x4 += d23;

	dx = x1 - x4;
	dy = y1 - y4;
	*lb = sqrt( dx * dx + dy * dy );
	dy = y1 + y4;
	*ub = sqrt( dx * dx + dy * dy );
}

static	void	get_measured_bounds( RESIDUE_T *res, int a1, int a2,
	float *lb, float *ub )
{

	*lb = res->r_bounds[ a2 ][ a1 ]; 
	if( *lb == 100.0 )
		*lb = 2.0;
	*ub = res->r_bounds[ a1 ][ a2 ]; 
	if( *ub == 0.0 )
		*ub = 12.0;
}

static	int	writematrix( char fname[], float **mat, int n )
{
	int	i, j;
	FILE	*fp;

	fp = fopen( fname, "w" );
	for( i = 0; i < n; i++ ){
		for( j = i + 1; j < n; j++ )
			fprintf( fp, "%4d %4d %10.4f %10.4f\n",
				i, j, mat[ j ][ i ], mat[ i ][ j ] );
	}
	fclose( fp );
}
