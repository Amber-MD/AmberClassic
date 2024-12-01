#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "nab.h"
#include "errormsg.h"
#include "traceback.h"
#include "memutil.h"
#include "molutil.h"
#include "select_atoms.h"

#define	B_NAME_SIZE	32
typedef	struct	bdbent	{
	char	b_name[ B_NAME_SIZE ];
	REAL_T	b_mean;
	REAL_T	b_stddev;
	REAL_T	b_scale;
	REAL_T	b_min;
	REAL_T	b_max;
} BDBENT;

#define	IBSIZE		40
#define	JBSIZE		10

#define	CLOSEST		0.5
#define	FURTHEST	10000.0

#define	EPSILON		1e-6

	/* used to add "slop" when using "model distances"	*/
#define	LB_MULT		0.99
#define	UB_MULT		1.01

	/* bounds ops	*/
#define	B_SET		0
#define	B_OR		1
#define	B_AND		2
#define	B_FROM_MODEL	3
#define	B_SHOW		4

	/* radii, from dgeom 95:	*/

#define	R_H	0.10
#define	R_N	1.55
#define	R_O	1.40
#define	R_C	1.65
#define	R_S	1.80
#define	R_P	1.90
#define	R_Fe	1.40
#define	R_X	1.80

#define	MAXRAND		2147483647

static	char	*b_opt;

	/* measured 1-4 bounds from NDB		*/
#define	BDB_SIZE	200
static	int	n_bdb;
static	BDBENT	bdb[ BDB_SIZE ];

	/* routines used to create bounds matrix:	*/

REAL_T	distp();

BOUNDS_T	*newbounds( MOLECULE_T *, char [] );
static	BOUNDS_T	*bsetup( MOLECULE_T * );
static	int	mk_bmat( MOLECULE_T *, BOUNDS_T * );
static	void	resolve_intbonds( BOUNDS_T *, int, int );
static	void	resolve_extbonds( BOUNDS_T *, STRAND_T *, int );
static	int	is_12_connected( BOUNDS_T *, int, int, REAL_T *, REAL_T * );
static	int	is_13_connected( BOUNDS_T *, int, int *, int,
	REAL_T *, REAL_T * );
static	int	is_14_connected( BOUNDS_T *, int, int *, int *, int,
	REAL_T *, REAL_T * );
#if 0
static	int	set_14_cis( BOUNDS_T *, int, int, int, int );
static	int	set_14_trans( BOUNDS_T *, int, int, int, int );
#endif
static	int	setbound( BOUNDS_T *, int, int, REAL_T, REAL_T );
static	int	showbound( BOUNDS_T *, int, int );
static	int	orbound( BOUNDS_T *, int, int, REAL_T, REAL_T );
static	int	andbound( BOUNDS_T *, int, int, REAL_T, REAL_T );
REAL_T	vdw_radius( BOUNDS_T *, int );
static	void	read_bdb( int *, BDBENT [] );

#if 0
	/* these 2 funcs are not used	*/
static	int	lookup_b14( BOUNDS_T *, int, int, int, int,
	REAL_T *lb, REAL_T *ub );
#endif

	/* OK back to used funcs	*/
void	mkdgname( BOUNDS_T *, int, char [] );
int	dumpbounds( FILE *, BOUNDS_T *, int );
#if 0
static	void	details( FILE *, BOUNDS_T *, int, int, int,
	REAL_T **, REAL_T ** );
#endif
static	int	rd_bmat( char [], BOUNDS_T *, int );
static	int	boundsop( int, BOUNDS_T *, MOLECULE_T *, char [], char [],
	REAL_T, REAL_T );
int	setbounds( BOUNDS_T *, MOLECULE_T *, char [], char [], REAL_T, REAL_T );
int	andbounds( BOUNDS_T *, MOLECULE_T *, char [], char [], REAL_T, REAL_T );
int	orbounds( BOUNDS_T *, MOLECULE_T *, char [], char [], REAL_T, REAL_T );
int	usemodeldist( BOUNDS_T *, MOLECULE_T *, char [], char [] );
int	showbounds( BOUNDS_T *, MOLECULE_T *, char [], char [] );
int	setchivol( BOUNDS_T *, MOLECULE_T *, char [], char [], char [], char [],
	REAL_T );
int	useboundsfrom( BOUNDS_T *, MOLECULE_T *, char [], MOLECULE_T *, char [],
	REAL_T );
REAL_T	dumpboundsviolations( FILE *, BOUNDS_T *, REAL_T );
REAL_T	dumpchiviolations( FILE *, BOUNDS_T *, REAL_T );

char	*getenv();

BOUNDS_T	*newbounds( MOLECULE_T *mol, char options[] )
{
	BOUNDS_T	*bp;
	int			binary;
	
	b_opt = options ? options : "";

	if( !n_bdb )
		read_bdb( &n_bdb, bdb );

	if( !( bp = bsetup( mol ) ) ){
		fprintf( stderr, "Failed in bsetup!\n" );
		exit( 1 );
	}

	if( strstr( b_opt, "-nocov" ) ) return( bp );

	if( strstr( b_opt, "-binary") ) binary = 1;
	else binary = 0;

	if( !strstr( b_opt, "-rbm" ) ){
		mk_bmat( mol, bp );
	}else
		rd_bmat( b_opt, bp, binary );

	return( bp );
}

static	BOUNDS_T	*bsetup( MOLECULE_T *mol )
{
	int	r, tr;
	int	i, j, c, tc;
	int	ai, a1, ac;
	int	natoms, nres, nchi, schi;
	int	user_chi;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	BOUNDS_T	*bp;
	int	*a_off;
	ATOM_T	**a_index, *ap1;
	int	*r_nconnect;
	int	**r_connect;
	char	*op,*bfnp;
	REAL_T	**bmat;
	CHIRAL_T	*chi;

	if( ( bp = ( BOUNDS_T * )malloc( sizeof( BOUNDS_T ) ) ) == NULL ){
		rt_errormsg_s( TRUE, E_NOMEM_FOR_S, "bp" );
		return( NULL );
	}
	bp->b_mol = mol;

	for( nres = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		nres += sp->s_nresidues;
	}
	bp->b_nres = nres;

	if( ( a_off = ( int * )malloc( nres * sizeof( int ) ) ) == NULL ){
		rt_errormsg_s( TRUE, E_NOMEM_FOR_S, "a_off" );
		return( NULL );
	}
	bp->b_aoff = a_off;

	nchi = nres = natoms = 0;
	for( sp = mol->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++ ){
			a_off[ nres++ ] = natoms;
			res = sp->s_residues[ r ];
			natoms += res->r_natoms;
			nchi += res->r_nchiral;
		}
	}
	user_chi = 4;
	if( (op = strstr( b_opt, "-nchi" )) ){
		if( ( bfnp = strchr( op, '=' ) ) ) {
			bfnp++;
			if( *bfnp )
				user_chi = atoi(bfnp);	
		}
	}
	bp->b_natoms = natoms;
	bp->b_nchiral = nchi;
	schi = nchi + user_chi * nres; /* allow users to add chirality */
	bp->b_schiral = schi;

	if( ( a_index = ( ATOM_T ** )ipvector( 0, natoms - 1 ) ) == NULL ){
		rt_errormsg_s( TRUE, E_NOMEM_FOR_S, "a_index" );
		return( NULL );
	}
	bp->b_aindex = a_index;
	if( ( chi = ( CHIRAL_T * )malloc( schi * sizeof( CHIRAL_T ) ) )==NULL ){
		rt_errormsg_s( TRUE, E_NOMEM_FOR_S, "chi" );
		return( NULL );
	}
	bp->b_chiral = chi;
	for( tc = 0, tr = 0, ai = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++, tr++ ){
			res = sp->s_residues[ r ];
			for( ap1 = res->r_atoms, a1 = 0; a1 < res->r_natoms;
				a1++, ap1++ ){
				bp->b_aindex[ ai ] = ap1;
				ai++;
			}
			for( c = 0; c < res->r_nchiral; c++, tc++ ){
				for( ac = 0; ac < 4; ac++ ){
					chi[ tc ].c_anum[ ac ] =
						res->r_chiral[ c ].c_anum[ ac ]
						+ a_off[ tr ];
				}
				chi[ tc ].c_dist = res->r_chiral[ c ].c_dist;
			}
		}
	}
	if( ( r_nconnect = ( int * )ivector( 0, natoms - 1 ) ) == NULL ){
		rt_errormsg_s( TRUE, E_NOMEM_FOR_S, "r_nconnect" );
		return( NULL );
	} 
	bp->b_nconnect = r_nconnect;
	if( ( r_connect =
		( int ** )imatrix( 0, natoms - 1, 0, A_CONNECT_SIZE - 1 ) ) ==
		NULL ){
		rt_errormsg_s( TRUE, E_NOMEM_FOR_S, "r_connect" );
		return( NULL );
	}
	bp->b_connect = r_connect;

	if( ( bmat = matrix( 0, natoms - 1, 0, natoms - 1 ) ) == NULL ){
		rt_errormsg_s( TRUE, E_NOMEM_FOR_S, "bmat" );
		return( NULL );
	}
	for( i = 0; i < natoms; i++ ){
		bmat[ i ][ i ] = 0.0;
		for( j = i + 1; j < natoms; j++ ){
			bmat[ i ][ j ] = FURTHEST;
			bmat[ j ][ i ] = CLOSEST;
		}
	}
	bp->b_bmat = bmat;

#ifdef PRINT_RESCHI
{
int	a, c, ca;
ATOM_T	*ap;
CHIRAL_T	*cp;

	for( c = 0; c < nchi; c++ ){
		cp = &chi[ c ];
		ap = a_index[ cp->c_anum[ 0 ] ];	
		fprintf( stderr, "%-4s:", ap->a_residue->r_resname );
		for( a = 0; a < 4; a++ ){
			ap = a_index[ cp->c_anum[ a ] ];	
			fprintf(stderr,"%4s/%-4d",ap->a_atomname,cp->c_anum[a]);
		}
		fprintf( stderr, " %8.3f\n", cp->c_dist );
	}
	
}
	fprintf(stderr, "Set up %d chirality constraints\n", nchi );
#endif
	return( bp );
}

static	int	mk_bmat( MOLECULE_T *mol, BOUNDS_T *bp )
{
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap1;
	REAL_T		lb, ub;
	int		tr0, tr, r, a1, a2, ai, aj;

	for( tr0 = tr = 0, ai = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++, tr++ ){
			res = sp->s_residues[ r ];
			for( ap1 = res->r_atoms, a1 = 0; a1 < res->r_natoms;
				a1++, ap1++ ){
				resolve_intbonds( bp, ai, tr );
				ai++;
			}
		}
		resolve_extbonds( bp, sp, tr0 );
		tr0 = tr;
	}

	for( a1 = 0; a1 < bp->b_natoms; a1++ )
		bp->b_bmat[ a1 ][ a1 ] = 0.0;
	for( a1 = 0; a1 < bp->b_natoms - 1; a1++ ){
		for( a2 = a1 + 1; a2 < bp->b_natoms; a2++ ){
			if( is_12_connected( bp, a1, a2, &lb, &ub ) ){
				setbound( bp, a1, a2, lb, ub );
			}else if( is_13_connected( bp,a1,&ai,a2,&lb,&ub ) ){
				setbound( bp, a1, a2, lb, ub );
			}else if(
				is_14_connected( bp,a1,&ai,&aj,a2,&lb,&ub ) ){
				andbound( bp, a1, a2, lb, ub );
			}else{
				lb = vdw_radius( bp,a1 ) + vdw_radius( bp,a2 );
				ub = FURTHEST;
				andbound( bp, a1, a2, lb, ub );
			}
		}
	}

	return( 0 );
}

#if 0
static  void	dumpBoundList( BOUNDS_T *bp )
{
        int i, j, natoms;

        natoms = bp->b_natoms;
        fprintf( stderr, "dumpBoundList:\n" );
        for( i = 0; i < natoms; i++ ){
                fprintf( stderr, "%d: ", i );
                for( j = 0; j < natoms; j++ ){
                        if ( bp->b_boundList[ i ][ j ] == -1 )
                                continue;
                        fprintf(stderr, "%d ", bp->b_boundList[ i ][ j ] );
                }
                fprintf( stderr, "\n" );
        }
}
#endif

static	void	resolve_intbonds( BOUNDS_T *bp, int ai, int r )
{
	ATOM_T	*ap;
	int	c;

	ap = bp->b_aindex[ ai ];
	bp->b_nconnect[ ai ] = ap->a_nconnect;
	for( c = 0; c < ap->a_nconnect; c++ ){
		bp->b_connect[ ai ][ c ] = ap->a_connect[ c ] + bp->b_aoff[ r ];
	}
}

static	void	resolve_extbonds( BOUNDS_T *bp, STRAND_T *sp, int tr )
{
	int		r, r1;
	int		a, a1;
	int		c;
	RESIDUE_T	*res;
	EXTBOND_T	*ebp;

	for( r = 0; r < sp->s_nresidues; r++ ){
		res = sp->s_residues[ r ];
		for( ebp = res->r_extbonds; ebp; ebp = ebp->eb_next ){
			a = ebp->eb_anum - 1;
			r1 = ebp->eb_rnum - 1;
			a1 = ebp->eb_ranum - 1;
			a += bp->b_aoff[ tr + r ];
			a1 += bp->b_aoff[ tr + r1 ];
			c = bp->b_nconnect[ a ];
			bp->b_connect[ a ][ c ] = a1;
			bp->b_nconnect[ a ] = c + 1;
		}
	}
}

static	int	is_12_connected( BOUNDS_T *bp, int a1, int a2,
	REAL_T *lb, REAL_T *ub )
{
	int	c;
/*	REAL_T	llb, lub;   */
	ATOM_T	*ap1, *ap2;
	REAL_T	dx, dy, dz;

	for( c = 0; c < bp->b_nconnect[ a1 ]; c++ ){
		if( bp->b_connect[ a1 ][ c ] == a2 ){
			ap1 = bp->b_aindex[ a1 ];
			ap2 = bp->b_aindex[ a2 ];
			dx = ap1->a_pos[ 0 ] - ap2->a_pos[ 0 ];
			dy = ap1->a_pos[ 1 ] - ap2->a_pos[ 1 ];
			dz = ap1->a_pos[ 2 ] - ap2->a_pos[ 2 ];
			*lb = *ub = sqrt( dx * dx + dy * dy + dz * dz );
	/*		if( lookup_b12( bp, a1, a2, &llb, &lub ) ){
				*lb = llb;
				*ub = lub;
			}   */
			return( TRUE );
		}
	}
	return( FALSE );
}

static	int	is_13_connected( BOUNDS_T *bp, int a1, int *ai, int a2,
	REAL_T *lb, REAL_T *ub )
{
	int	c;
	REAL_T	llb, lub;
	ATOM_T	*ap1, *ap2;
	REAL_T	dx, dy, dz;

	for( c = 0; c < bp->b_nconnect[ a1 ]; c++ ){
		*ai = bp->b_connect[ a1 ][ c ];
		if( is_12_connected( bp, *ai, a2, &llb, &lub ) ){
			ap1 = bp->b_aindex[ a1 ];
			ap2 = bp->b_aindex[ a2 ];
			dx = ap1->a_pos[ 0 ] - ap2->a_pos[ 0 ];
			dy = ap1->a_pos[ 1 ] - ap2->a_pos[ 1 ];
			dz = ap1->a_pos[ 2 ] - ap2->a_pos[ 2 ];
			*lb = *ub = sqrt( dx * dx + dy * dy + dz * dz );
/*			if( lookup_b13( bp, a1, *ai, a2, &llb, &lub ) ){
				*lb = llb;
				*ub = lub;
			}  */
			return( TRUE );
		}
	}
	return( FALSE );
}

static	int	is_14_connected( BOUNDS_T *bp,
	int a1, int *ai, int *aj, int a2, REAL_T *lb, REAL_T *ub )
{
	int	c;
	REAL_T	llb, lub, clb;
	ATOM_T	*ap1, *ap2, *api, *apj;
	REAL_T	d12, d23, d13, d34, d24;
	REAL_T	dx, dy, dz;
	REAL_T	x1, y1, x2, y2;
	REAL_T	cos13, cos34;

	for( c = 0; c < bp->b_nconnect[ a1 ]; c++ ){
		*ai = bp->b_connect[ a1 ][ c ];
		if( is_13_connected( bp, *ai, aj, a2, &llb, &lub ) ){
			ap1 = bp->b_aindex[ a1 ];
			api = bp->b_aindex[ *ai ];
			apj = bp->b_aindex[ *aj ];
			ap2 = bp->b_aindex[ a2 ];

			dx = ap1->a_pos[ 0 ] - api->a_pos[ 0 ];
			dy = ap1->a_pos[ 1 ] - api->a_pos[ 1 ];
			dz = ap1->a_pos[ 2 ] - api->a_pos[ 2 ];
			d12 = sqrt( dx * dx + dy * dy + dz * dz );
			dx = api->a_pos[ 0 ] - apj->a_pos[ 0 ];
			dy = api->a_pos[ 1 ] - apj->a_pos[ 1 ];
			dz = api->a_pos[ 2 ] - apj->a_pos[ 2 ];
			d23 = sqrt( dx * dx + dy * dy + dz * dz );
			dx = ap1->a_pos[ 0 ] - apj->a_pos[ 0 ];
			dy = ap1->a_pos[ 1 ] - apj->a_pos[ 1 ];
			dz = ap1->a_pos[ 2 ] - apj->a_pos[ 2 ];
			d13 = sqrt( dx * dx + dy * dy + dz * dz );

			dx = apj->a_pos[ 0 ] - ap2->a_pos[ 0 ];
			dy = apj->a_pos[ 1 ] - ap2->a_pos[ 1 ];
			dz = apj->a_pos[ 2 ] - ap2->a_pos[ 2 ];
			d34 = sqrt( dx * dx + dy * dy + dz * dz );
			dx = api->a_pos[ 0 ] - ap2->a_pos[ 0 ];
			dy = api->a_pos[ 1 ] - ap2->a_pos[ 1 ];
			dz = api->a_pos[ 2 ] - ap2->a_pos[ 2 ];
			d24 = sqrt( dx * dx + dy * dy + dz * dz );

			cos13 = (d12*d12+d23*d23-d13*d13)/(2.*d12*d23);
			cos34 = (d23*d23+d24*d24-d34*d34)/(2.*d23*d24);

			x1 = d12 * cos13;
			y1 = sqrt( d12 * d12 - x1 * x1 );

			x2 = d24 * cos34;
			y2 = sqrt( d24 * d24 - x2 * x2 );

			dx = x1 - x2;
			dy = y1 - y2;
			clb = *lb = sqrt( dx * dx + dy * dy );

			if( ( llb = vdw_radius( bp, a1 ) + 
				vdw_radius( bp, a2 ) ) > *lb )
				*lb = llb;

			y2 = -y2;
			dy = y1 - y2;
			*ub = sqrt( dx * dx + dy * dy );

			if( *ub < *lb )
				*ub = *lb;

/*			if( set_14_cis( bp, a1, *ai, *aj, a2 ) )  */
				*lb = clb;

/*			if( set_14_trans( bp, a1, *ai, *aj, a2 ) )
				*lb = *ub;


			if( lookup_b14( bp, a1, *ai, *aj, a2, &llb, &lub ) ){
				*lb = llb;
				*ub = lub;
			}
*/
			return( TRUE );
		}
	}
	return( FALSE );
}

#if 0
static	int	set_14_cis( BOUNDS_T *bp, int a1, int a2, int a3, int a4 )
{
	ATOM_T	*a1p, *a2p, *a3p, *a4p;
	RESIDUE_T	*r1p, *r3p;
	int	set;

	a1p = bp->b_aindex[ a1 ];
	a2p = bp->b_aindex[ a2 ];
	a3p = bp->b_aindex[ a3 ];
	a4p = bp->b_aindex[ a4 ];

	if(!strcmp( a1p->a_atomname, "CA" ) 	&&
		!strcmp( a2p->a_atomname, "C" )	&&
		!strcmp( a3p->a_atomname, "N" ) 	&&
		!strcmp( a4p->a_atomname, "HN" ) )
		set = TRUE;
	else if(!strcmp( a1p->a_atomname, "O" ) 	&&
		!strcmp( a2p->a_atomname, "C" )	&&
		!strcmp( a3p->a_atomname, "N" ) 	&&
		!strcmp( a4p->a_atomname, "CA" ) )
		set = TRUE;
	else
		return( FALSE );

	r1p = a1p->a_residue;
	r3p = a3p->a_residue;

	if( !strcmp( r3p->r_resname, "PRO" ) )
		return( FALSE );

	return( TRUE );
}

static	int	set_14_trans( BOUNDS_T *bp, int a1, int a2, int a3, int a4 )
{
	ATOM_T	*a1p, *a2p, *a3p, *a4p;
	RESIDUE_T	*r1p, *r3p;
	int	set;

	a1p = bp->b_aindex[ a1 ];
	a2p = bp->b_aindex[ a2 ];
	a3p = bp->b_aindex[ a3 ];
	a4p = bp->b_aindex[ a4 ];

	if( 	!strcmp( a1p->a_atomname, "CA" )	&&
		!strcmp( a2p->a_atomname, "C" )	&&
		!strcmp( a3p->a_atomname, "N" ) 	&&
		!strcmp( a4p->a_atomname, "CA" ) )
		set = TRUE;
	else if(!strcmp( a1p->a_atomname, "O" ) 	&&
		!strcmp( a2p->a_atomname, "C" )	&&
		!strcmp( a3p->a_atomname, "N" ) 	&&
		!strcmp( a4p->a_atomname, "HN" ) )
		set = TRUE;
	else
		return( FALSE );

	r1p = a1p->a_residue;
	r3p = a3p->a_residue;

	if( !strcmp( r3p->r_resname, "PRO" ) )
		return( FALSE );

	return( TRUE );
}
#endif

static	int	setbound( BOUNDS_T *bp, int a1, int a2, REAL_T lb, REAL_T ub )
{
	char	a1nm[ 40 ], a2nm[ 40 ];
	REAL_T beps = 1.e-9;
	
	/* the following is necessary since somtimes apparently equal */
	/* lb and ub are evaluated as lb >ub, leading to an error     */
	if ((lb > ub) && (lb - ub < beps)) lb = ub;
	
	if( lb > ub ){
		mkdgname( bp, a1, a1nm );
		mkdgname( bp, a2, a2nm );
		fprintf( stderr, "setbound: ERROR: lb > ub: %s -> %s %f %f\n",
			a1nm, a2nm, lb, ub );
		return( 0 );
	}
	bp->b_bmat[ a1 ][ a2 ] = ub;
	bp->b_bmat[ a2 ][ a1 ] = lb;
	return( 1 );
}

static	int	showbound( BOUNDS_T *bp, int a1, int a2 )
{
	char	a1nm[ 40 ], a2nm[ 40 ];
	mkdgname( bp, a1, a1nm );
	mkdgname( bp, a2, a2nm );
	
	if (a1 >= a2) {
		fprintf( stderr, "showbounds:  %s-%s,\t%f - %f\n",
		a1nm, a2nm, bp->b_bmat[ a1 ][ a2 ], bp->b_bmat[ a2 ][ a1 ]); 
	}
	else {
		fprintf( stderr, "showbounds:  %s-%s,\t%f - %f\n",
		a1nm, a2nm, bp->b_bmat[ a2 ][ a1 ], bp->b_bmat[ a1 ][ a2 ]); 
	}
	return( 1 );
}

static	int	orbound( BOUNDS_T *bp, int a1, int a2, REAL_T lb, REAL_T ub )
{
	char	a1nm[ 40 ], a2nm[ 40 ];
	int	rval;

	if( lb > ub ){
		mkdgname( bp, a1, a1nm );
		mkdgname( bp, a2, a2nm );
		fprintf( stderr, "orbound: ERROR: lb > ub: %s -> %s %f %f\n",
			a1nm, a2nm, lb, ub );
		return( 0 );
	}
	rval = 0;
	if( ub > bp->b_bmat[ a1 ][ a2 ] ){
		bp->b_bmat[ a1 ][ a2 ] = ub;
		rval = 1;
	}
	if( lb < bp->b_bmat[ a2 ][ a1 ] ){
		bp->b_bmat[ a2 ][ a1 ] = lb;
		rval = 1;
	}
	return( rval );
}

static	int	andbound( BOUNDS_T *bp, int a1, int a2, REAL_T lb, REAL_T ub )
{
	char	a1nm[ 40 ], a2nm[ 40 ];
	int	rval;

	if( lb > ub ){
		mkdgname( bp, a1, a1nm );
		mkdgname( bp, a2, a2nm );
		fprintf( stderr, "andbound: ERROR: lb > ub: %s -> %s %f %f\n",
			a1nm, a2nm, lb, ub );
		return( 0 );
	}
	rval = 0;
	if( ub < bp->b_bmat[ a1 ][ a2 ] ){
		if( ub >= bp->b_bmat[ a2 ][ a1 ] ){
			bp->b_bmat[ a1 ][ a2 ] = ub;
			rval = 1;
		}
	}
	if( lb > bp->b_bmat[ a2 ][ a1 ] ){
		if( lb <= bp->b_bmat[ a1 ][ a2 ] ){
			bp->b_bmat[ a2 ][ a1 ] = lb;
			rval = 1;
		}
	}
	return( rval );
}

REAL_T	vdw_radius( BOUNDS_T *bp, int a )
{
	ATOM_T	*ap;
	REAL_T	r, delta = .00001;

	ap = bp->b_aindex[ a ];

/*  Use value in a_radius if it is not zero, otherwise use default values */

	if( ap->a_radius > delta ) return( ap->a_radius );
		
	switch( *ap->a_atomname ){
	case '1' :
	case '2' :
	case 'H' :
	case 'h' :
		r = R_H;
		break;
	case 'N' :
	case 'n' :
		r = R_N;
		break;
	case 'O' :
	case 'o' :
		r = R_O;
		break;
	case 'C' :
	case 'c' :
		r = R_C;
		break;
	case 'S' :
	case 's' :
		r = R_S;
		break;
	case 'P' :
	case 'p' :
		r = R_P;
		break;
	case 'F' :
	case 'f' :
		if( ap->a_atomname[ 1 ] == 'E' || ap->a_atomname[ 1 ] == 'e' )
			r = R_Fe;
		else
			r = R_X;
		break;
	default :
		r = R_X;
	}
	return( r );
}

static	void	read_bdb( int *n_bdb, BDBENT bdb[] )
{
	FILE	*fp;
	char	line[ 256 ];
	BDBENT	*bp;
	char	bdfname[ 256 ];

	sprintf( bdfname, "%s/dat/reslib/bounds.db", getenv("AMBERHOME") );

	if( ( fp = fopen( bdfname, "r" ) ) == NULL ){
		fprintf( stderr, "read_b14tab: can't open db %s\n", bdfname );
		return;
	}

	for( *n_bdb = 0; fgets( line, sizeof( line ), fp ); ){
		if( *line == '#' )
			continue;
		bp = &bdb[ *n_bdb ];
#ifdef NAB_DOUBLE_PRECISION
		sscanf( line, "%s %lf %lf %lf %lf %lf",
#else
		sscanf( line, "%s %f %f %f %f %f",
#endif
			bp->b_name, &bp->b_mean, &bp->b_stddev,
			&bp->b_scale, &bp->b_min, &bp->b_max );
		( *n_bdb )++;
	}

	fclose( fp );
}

#if 0
static	int	lookup_b12( BOUNDS_T *bp, int a1, int a2,
	REAL_T *lb, REAL_T *ub )
{
	ATOM_T	*ap1, *ap2;
	char	bname[ B_NAME_SIZE ];
	int	b;

	ap1 = bp->b_aindex[ a1 ];
	ap2 = bp->b_aindex[ a2 ];
	if( strcmp( ap1->a_atomname, ap2->a_atomname ) < 0 )
		sprintf( bname, "%s-%s",
			ap1->a_atomname, ap2->a_atomname );
	else
		sprintf( bname, "%s-%s",
			ap2->a_atomname, ap1->a_atomname );
	for( b = 0; b < n_bdb; b++ ){
		if( !strcmp( bdb[ b ].b_name, bname ) ){
			*lb = bdb[ b ].b_mean;
			*ub = bdb[ b ].b_mean;
			return( 1 );
		}
	}
	return( 0 );
}

static	int	lookup_b13( BOUNDS_T *bp, int a1, int a2, int a3,
	REAL_T *lb, REAL_T *ub )
{
	ATOM_T	*ap1, *ap2, *ap3;
	char	bname[ B_NAME_SIZE ];
	int	b;

	ap1 = bp->b_aindex[ a1 ];
	ap2 = bp->b_aindex[ a2 ];
	ap3 = bp->b_aindex[ a3 ];
	if( strcmp( ap1->a_atomname, ap3->a_atomname ) < 0 )
		sprintf( bname, "%s-%s-%s",
			ap1->a_atomname, ap2->a_atomname, ap3->a_atomname );
	else
		sprintf( bname, "%s-%s-%s",
			ap3->a_atomname, ap2->a_atomname, ap1->a_atomname );
	for( b = 0; b < n_bdb; b++ ){
		if( !strcmp( bdb[ b ].b_name, bname ) ){
			*lb = bdb[ b ].b_mean;
			*ub = bdb[ b ].b_mean;
			return( 1 );
		}
	}
	return( 0 );
}

static	int	lookup_b14( BOUNDS_T *bp, int a1, int a2, int a3, int a4,
	REAL_T *lb, REAL_T *ub )
{
	REAL_T	d12, d23, d34, d13, d24;
	REAL_T	junk;
	REAL_T	cos13, cos34;
	int	x1, y1, x2, y2;
	REAL_T	dx, dy;
	REAL_T	llb;

	if( !lookup_b12( bp, a1, a2, &d12, &junk ) )
		return( 0 );
	if( !lookup_b12( bp, a2, a3, &d23, &junk ) )
		return( 0 );
	if( !lookup_b12( bp, a3, a4, &d34, &junk ) )
		return( 0 );
	if( !lookup_b13( bp, a1, a2, a3, &d13, &junk ) )
		return( 0 );
	if( !lookup_b13( bp, a2, a3, a4, &d24, &junk ) )
		return( 0 );

	cos13 = (d12*d12+d23*d23-d13*d13)/(2.*d12*d23);
	cos34 = (d23*d23+d24*d24-d34*d34)/(2.*d23*d24);

	x1 = d12 * cos13;
	y1 = sqrt( d12 * d12 - x1 * x1 );

	x2 = d24 * cos34;
	y2 = sqrt( d24 * d24 - x2 * x2 );

	dx = x1 - x2;
	dy = y1 - y2;
	*lb = sqrt( dx * dx + dy * dy );

	if( ( llb = vdw_radius( bp, a1 ) + 
		vdw_radius( bp, a2 ) ) > *lb )
		*lb = llb;

	y2 = -y2;
	dy = y1 - y2;
	*ub = sqrt( dx * dx + dy * dy );

	if( *ub < *lb )
		*ub = *lb;

	return( 1 );
}
#endif

void	mkdgname( BOUNDS_T *bp, int i, char nm[] )
{
	ATOM_T	*ap;
	RESIDUE_T	*res;
	STRAND_T	*sp;
	MOLECULE_T	*mp;
	int		r, anl;
	char		anm[ 40 ];

	ap = bp->b_aindex[ i ];
	res = ap->a_residue;
	sp = res->r_strand;
	mp = sp->s_molecule;
	if( !mp->m_nvalid )
		upd_molnumbers( mp );
	for( r = 0; r < sp->s_nresidues; r++ ){
		if( sp->s_residues[ r ] == res ){ break; }
	}
	if( isdigit( *ap->a_atomname ) ){
		anl = strlen( ap->a_atomname );
		strcpy( anm, &ap->a_atomname[ 1 ] );
		anm[ anl - 1 ] = *ap->a_atomname;
		anm[ anl ] = '\0';
	}else
		strcpy( anm, ap->a_atomname );
	sprintf( nm, "%s:%d_%s:%s",
		sp->s_strandname, res->r_resnum, res->r_resname, anm );
}

int	dumpbounds( FILE *fp, BOUNDS_T *bp, int binary )
{
	int	a1, a2, cnt, nchiral;
/*	char	a1nm[ 20 ], a2nm[ 20 ];   */

	if( binary ){
		for( a1 = 0; a1 < bp->b_natoms; a1++ ){
			fwrite( bp->b_bmat[ a1 ], sizeof( REAL_T ), bp->b_natoms, fp );
		}
		cnt = 0;
		nchiral = bp->b_nchiral;
		fwrite( &nchiral, sizeof( int ), 1, fp );
		fwrite( bp->b_chiral, sizeof( CHIRAL_T ), nchiral, fp );
		fprintf( stderr, "wrote %d chiral volumes\n", nchiral );
	} else {
		for( cnt = 0, a1 = 0; a1 < bp->b_natoms - 1; a1++ ){
/*			mkdgname( bp, a1, a1nm );    */
			for( a2 = a1 + 1; a2 < bp->b_natoms; a2++, cnt++ ){
/*
**				mkdgname( bp, a2, a2nm );
**				fprintf( fp, "%-12s %-12s %d %d %f %f\n",
**					a1nm, a2nm, a1, a2,
**					bp->b_bmat[ a2 ][ a1 ],bp->b_bmat[ a1 ][ a2 ] );
*/
				if( bp->b_bmat[ a2 ][ a1 ] == CLOSEST &&
			    	bp->b_bmat[ a1 ][ a2 ] == FURTHEST ) continue;
				fprintf( fp, "%4d %4d %7.2f %7.2f\n",
					a1, a2, bp->b_bmat[ a2 ][ a1 ],bp->b_bmat[ a1 ][ a2 ] );
			}
		}
	}
	return( cnt );
}

#if 0
static	void	details( FILE *fp, BOUNDS_T *bp,
	int i, int j, int k, REAL_T **l_bmat, REAL_T **u_bmat )
{
	char	inm[ 40 ], jnm[ 40 ], knm[ 40 ];

	mkdgname( bp, i, inm );
	mkdgname( bp, j, jnm );
	mkdgname( bp, k, knm );
	fprintf( fp, "i - j - k: %s - %s - %s\n", inm, jnm, knm );
	fprintf( stderr, "ub[ %s,%s ] = %9.3f\n", inm, jnm, u_bmat[ i ][ j ] );
	fprintf( stderr, "ub[ %s,%s ] = %9.3f\n", jnm, knm, u_bmat[ j ][ k ] );
	fprintf( stderr, "ub[ %s,%s ] = %9.3f\n", inm, knm, u_bmat[ i ][ k ] );
	fprintf( stderr, "lb[ %s,%s ] = %9.3f\n", inm, jnm, l_bmat[ i ][ j ] );
	fprintf( stderr, "lb[ %s,%s ] = %9.3f\n", jnm, knm, l_bmat[ j ][ k ] );
	fprintf( stderr, "lb[ %s,%s ] = %9.3f\n", inm, knm, l_bmat[ i ][ k ] );
}
#endif

static	int	rd_bmat( char op[], BOUNDS_T *bp, int binary )
{
	char	*bmfname;
	FILE	*bmfp;
	int		i, ai, aj;
	char	line[ 256 ];
/*	char	ainame[ 50 ], ajname[ 50 ];   */
	REAL_T	l, u;
	int		nchiral;

	if( !( bmfname = strchr( op, '=' ) ) ){
		fprintf( stderr, "rd_bmat: syntax: %s\n", op );
		exit( 1 );
	}
	bmfname++;
	if( !*bmfname ){
		fprintf( stderr, "rd_bmat: syntax: %s\n", op );
		exit( 1 );
	}
	if( ( bmfp = fopen( bmfname, "r" ) ) == NULL ){
		fprintf( stderr, "rd_bmat: can't open bounds mat file %s\n",
			bmfname );
		exit( 1 );
	}
	for( i = 0; i < bp->b_natoms; i++ )
		bp->b_bmat[ i ][ i ] = 0.0;

	if( binary ){
		for( ai=0; ai<bp->b_natoms; ai++ ){
			fread( bp->b_bmat[ ai ], sizeof( REAL_T ), bp->b_natoms, bmfp );
		}
		fread( &nchiral, sizeof( int ), 1, bmfp );
		bp->b_nchiral = nchiral;
		if( nchiral > bp->b_schiral ){
			fprintf( stderr, "no room for reading chiral restraints: %d %d\n",
				nchiral, bp->b_schiral );
			exit(1);
		}
		fread( bp->b_chiral, sizeof( CHIRAL_T ), nchiral, bmfp );
		fprintf( stderr, "read %d chiral volumes\n", nchiral );
	} else {
		while( fgets( line, sizeof( line ), bmfp ) ){
#ifdef NAB_DOUBLE_PRECISION
			sscanf( line , "%d %d %lf %lf", &ai, &aj, &l, &u );
#else
			sscanf( line , "%d %d %f %f", &ai, &aj, &l, &u );
#endif
			bp->b_bmat[ ai ][ aj ] = u;
			bp->b_bmat[ aj ][ ai ] = l;
		}
	}
	fclose( bmfp );
	return( 0 );
}

static	int	boundsop( int bop, BOUNDS_T *bp, MOLECULE_T *mol,
	char bex1[], char bex2[], REAL_T lb, REAL_T ub )
{
	int	ai, aj;
	int	cnt;
	ATOM_T	*api, *apj;
	REAL_T	dx, dy, dz, d;

	clear_attr( mol, AT_SELECTED );
	select_atoms( mol, bex1 );
	set_attr_if( mol, AT_SELECTED, AT_SELECT );
	select_atoms( mol, bex2 );
	for( cnt = 0, ai = 0; ai < bp->b_natoms; ai++ ){
		api = bp->b_aindex[ ai ];
		if( !( api->a_attr & AT_SELECTED ) )
			continue;
		for( aj = 0; aj < bp->b_natoms; aj++ ){
			if( aj == ai )
				continue;
			apj = bp->b_aindex[ aj ];
			if( !( apj->a_attr & AT_SELECT ) )
				continue;

		   if (ai > aj) {
			switch( bop ){
			case B_SET :
				cnt += setbound( bp, aj, ai, lb, ub );
				break;
			case B_AND :
				cnt += andbound( bp, aj, ai, lb, ub );
				break;
			case B_OR :
				cnt += orbound( bp, aj, ai, lb, ub );
				break;
			case B_FROM_MODEL :
				dx = apj->a_pos[ 0 ] - api->a_pos[ 0 ];
				dy = apj->a_pos[ 1 ] - api->a_pos[ 1 ];
				dz = apj->a_pos[ 2 ] - api->a_pos[ 2 ];
				d = sqrt( dx * dx + dy * dy + dz * dz );
				cnt += setbound( bp, aj, ai, d, d );
				break;
			case B_SHOW :
				cnt += showbound( bp, ai, aj );
				break;
			}
		    }
		    else {
			switch( bop ){
			case B_SET :
				cnt += setbound( bp, ai, aj, lb, ub );
				break;
			case B_AND :
				cnt += andbound( bp, ai, aj, lb, ub );
				break;
			case B_OR :
				cnt += orbound( bp, ai, aj, lb, ub );
				break;
			case B_FROM_MODEL :
				dx = api->a_pos[ 0 ] - apj->a_pos[ 0 ];
				dy = api->a_pos[ 1 ] - apj->a_pos[ 1 ];
				dz = api->a_pos[ 2 ] - apj->a_pos[ 2 ];
				d = sqrt( dx * dx + dy * dy + dz * dz );
				cnt += setbound( bp, ai, aj, d, d );
				break;
			case B_SHOW :
				cnt += showbound( bp, ai, aj );
				break;
			}
		    }
		}
	}
	return( cnt );
}

int	setbounds( BOUNDS_T *bp, MOLECULE_T *mol,
	char bex1[], char bex2[], REAL_T lb, REAL_T ub )
{
	int	cnt;

	cnt = boundsop( B_SET, bp, mol, bex1, bex2, lb, ub );
	return( cnt );
}

int	andbounds( BOUNDS_T *bp, MOLECULE_T *mol,
	char bex1[], char bex2[], REAL_T lb, REAL_T ub )
{
	int	cnt;

	cnt = boundsop( B_AND, bp, mol, bex1, bex2, lb, ub );
	return( cnt );
}

int	orbounds( BOUNDS_T *bp, MOLECULE_T *mol,
	char bex1[], char bex2[], REAL_T lb, REAL_T ub )
{
	int	cnt;

	cnt = boundsop( B_OR, bp, mol, bex1, bex2, lb, ub );
	return( cnt );
}

int	usemodeldist( BOUNDS_T *bp, MOLECULE_T *mol,
	char bex1[], char bex2[] )
{
	int	cnt;

	cnt = boundsop( B_FROM_MODEL, bp, mol, bex1, bex2, CLOSEST, FURTHEST );
	return( cnt );
}

int	showbounds( BOUNDS_T *bp, MOLECULE_T *mol,
	char bex1[], char bex2[] )
{
	int	cnt;

	cnt = boundsop( B_SHOW, bp, mol, bex1, bex2, CLOSEST, FURTHEST );
	return( cnt );
}

int	setchivol( BOUNDS_T *bp, MOLECULE_T *mol,
	char aex1[], char aex2[], char aex3[], char aex4[], REAL_T vol )
{
	ATOM_T	*api;
	int	ai, c;
	int	anums[ 4 ];
	CHIRAL_T	*chi;

	if( bp->b_nchiral >= bp->b_schiral ){
		fprintf( stderr, "setchivol: chiral table overflow\n" );
		exit( 1 );
	}

	if ( countmolatoms( mol, aex1 ) != 1 ){
		fprintf( stderr, "setchivol got bad selection: %s\n", aex1 );
		exit( 1 );
	}
	for( ai = 0; ai < bp->b_natoms; ai++ ){
		api = bp->b_aindex[ ai ];
		if(  api->a_attr & AT_SELECT ) { anums[ 0 ] = ai; break; }
	}

	if ( countmolatoms( mol, aex2 ) != 1 ){
		fprintf( stderr, "setchivol got bad selection: %s\n", aex2 );
		exit( 1 );
	}
	for( ai = 0; ai < bp->b_natoms; ai++ ){
		api = bp->b_aindex[ ai ];
		if(  api->a_attr & AT_SELECT ) { anums[ 1 ] = ai; break; }
	}

	if ( countmolatoms( mol, aex3 ) != 1 ){
		fprintf( stderr, "setchivol got bad selection: %s\n", aex3 );
		exit( 1 );
	}
	for( ai = 0; ai < bp->b_natoms; ai++ ){
		api = bp->b_aindex[ ai ];
		if(  api->a_attr & AT_SELECT ) { anums[ 2 ] = ai; break; }
	}

	if ( countmolatoms( mol, aex4 ) != 1 ){
		fprintf( stderr, "setchivol got bad selection: %s\n", aex4 );
		exit( 1 );
	}
	for( ai = 0; ai < bp->b_natoms; ai++ ){
		api = bp->b_aindex[ ai ];
		if(  api->a_attr & AT_SELECT ) { anums[ 3 ] = ai; break; }
	}

	chi = &bp->b_chiral[ bp->b_nchiral++ ];
	for( c = 0; c < 4; c++ )
		chi->c_anum[ c ] = anums[ c ];
	chi->c_dist = vol;

#ifdef PRINT_CHIRAL_VOLUMES
	{
	int	a, ca;
	ATOM_T	*ap;

		fprintf( stderr, "scv: new chiral center:\n" );
		for( a = 0; a < 4; a++ ){
			ap = bp->b_aindex[ chi->c_anum[ a ] ];	
			fprintf( stderr, "%-4s:", ap->a_residue->r_resname );
			fprintf( stderr,"%4s/%-4d",ap->a_atomname,chi->c_anum[a] );
		}
		fprintf( stderr, " %8.3f\n", chi->c_dist );
	}
#endif
	return( 0 );
}

int	useboundsfrom( BOUNDS_T *b, MOLECULE_T *m1, char ae1[],
	MOLECULE_T *m2, char ae2[], REAL_T slack )
{
	int	n_idx1, n_idx2;
	int	*idx1 = NULL;
	ATOM_T	**idx2 = NULL;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap;
	int	r, a, ta, sa;
	int	i, j;
	int	bi, bj;
	ATOM_T	*aip, *ajp;
	REAL_T	lb, ub;

	n_idx1 = countmolatoms( m1, ae1 );
	if( ( idx1 = ( int * )malloc( n_idx1 * sizeof( int ) ) ) == NULL ){
		fprintf( stderr, "useboundsfrom: can't allocate idx1\n" );
		exit( 1 );
	}
	for( sa = ta = 0, sp = m1->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			ap = res->r_atoms;
			for( a = 0; a < res->r_natoms; a++, ap++ ){
				if( ap->a_attr & AT_SELECT ){
					idx1[ sa ] = ta;
					sa++;
				}
				ta++;
			}
		}
	}

	n_idx2 = countmolatoms( m2, ae2 );
	if( n_idx1 != n_idx2 ){
		free( idx1 );
		fprintf( stderr, "useboundsfrom: atom mismatch: m1(%d) & m2(%d)\n",
			n_idx1, n_idx2 );
		return( 1 );
	}
	if( ( idx2=( ATOM_T ** )malloc( n_idx2*sizeof( ATOM_T * ) ) ) == NULL ){
		free( idx1 );
		fprintf( stderr, "useboundsfrom: can't allocate idx2\n" );
		exit( 1 );
	}
	for( sa = 0, sp = m2->m_strands; sp; sp = sp->s_next ){
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			ap = res->r_atoms;
			for( a = 0; a < res->r_natoms; a++, ap++ ){
				if( ap->a_attr & AT_SELECT ){
					idx2[ sa ] = ap;
					sa++;
				} 
			}
		}
	}

	for( i = 0; i < n_idx1 - 1; i++ ){
		bi = idx1[ i ];
		aip = idx2[ i ];
		for( j = i + 1; j < n_idx1; j++ ){
			bj = idx1[ j ];
			ajp = idx2[ j ];
			lb = ub = distp( aip->a_pos, ajp->a_pos );
			if( ( lb - slack*lb ) < 0.0 )
				b->b_bmat[ bj ][ bi ] = 0.0;
			else
				b->b_bmat[ bj ][ bi ] = lb - slack*lb;
			b->b_bmat[ bi ][ bj ] = ub + slack*ub;
		}
	}

	free( idx1 );
	free( idx2 );
	fprintf( stderr, "useboundsfrom set bounds for %d atoms\n", n_idx1 );
	return( 0 );
}

REAL_T	dumpboundsviolations( FILE *fp, BOUNDS_T *bp, REAL_T cutoff )
{
	int	i, j;
	ATOM_T	*ai, *aj;
	RESIDUE_T	*resi, *resj;
	REAL_T	dx, dy, dz, d;
	REAL_T	viol, tviol;

	NAB_ari( bp->b_aindex[ 0 ], "resnum" );
	for( tviol = 0.0, i = 0; i < bp->b_natoms - 1; i++ ){
		ai = bp->b_aindex[ i ];
		resi = ai->a_residue;
		for( j = i + 1; j < bp->b_natoms; j++ ){
			aj = bp->b_aindex[ j ];
			resj = aj->a_residue;
			dx = ai->a_pos[ 0 ] - aj->a_pos[ 0 ];
			dy = ai->a_pos[ 1 ] - aj->a_pos[ 1 ];
			dz = ai->a_pos[ 2 ] - aj->a_pos[ 2 ];
			d = sqrt( dx * dx + dy * dy + dz * dz );
			if( ( viol = bp->b_bmat[ j ][ i ] - d ) > cutoff ){
				tviol += viol;
				fprintf( fp,
				"%4d %-4s %-4s %4d %-4s %-4s %8.3f %8.3f %8.3f %8.3f\n",
				resi->r_tresnum,resi->r_resname,ai->a_atomname,
				resj->r_tresnum,resj->r_resname,aj->a_atomname,
				bp->b_bmat[j][i],
				bp->b_bmat[i][j],
				d, -viol );
			}else if( ( viol = d - bp->b_bmat[ i ][ j ]) > cutoff ){
				tviol += viol;
				fprintf( fp,
				"%4d %-4s %-4s %4d %-4s %-4s %8.3f %8.3f %8.3f %8.3f\n",
				resi->r_tresnum,resi->r_resname,ai->a_atomname,
				resj->r_tresnum,resj->r_resname,aj->a_atomname,
				bp->b_bmat[j][i],
				bp->b_bmat[i][j],
				d, viol );
			}
		}
	}
	return( tviol );
}

REAL_T	dumpchiviolations( FILE *fp, BOUNDS_T *bp, REAL_T cutoff )
{
	int	i, nchi;
	CHIRAL_T    *chi, *cp;   /* chirality tetrads    */
	ATOM_T  *ai, *aj, *ak, *al;
	REAL_T	viol, tviol;
	REAL_T	x0, y0, z0;
	REAL_T	x1, y1, z1;
	REAL_T	x2, y2, z2;
	REAL_T	x3, y3, z3;
	REAL_T	a1,a2,a3, b1,b2,b3, c1,c2,c3;
	REAL_T	gq1, gq2, gq3, vol;

	nchi = bp->b_nchiral;
	chi = bp->b_chiral;

	/*   chirality violations:    */

	NAB_ari( bp->b_aindex[ 0 ], "resnum" );
	for( tviol=0.0, i = 0; i < nchi; i++ ){

		cp = &chi[ i ];
		ai = bp->b_aindex[ cp->c_anum[0] ];
		aj = bp->b_aindex[ cp->c_anum[1] ];
		ak = bp->b_aindex[ cp->c_anum[2] ];
		al = bp->b_aindex[ cp->c_anum[3] ];

		x0 = ai->a_pos[ 0 ];
		y0 = ai->a_pos[ 1 ];
		z0 = ai->a_pos[ 2 ];
		x1 = aj->a_pos[ 0 ];
		y1 = aj->a_pos[ 1 ];
		z1 = aj->a_pos[ 2 ];
		x2 = ak->a_pos[ 0 ];
		y2 = ak->a_pos[ 1 ];
		z2 = ak->a_pos[ 2 ];
		x3 = al->a_pos[ 0 ];
		y3 = al->a_pos[ 1 ];
		z3 = al->a_pos[ 2 ];

		a1 = x1 - x0; a2 = y1 - y0; a3 = z1 - z0;
		b1 = x2 - x0; b2 = y2 - y0; b3 = z2 - z0;
		c1 = x3 - x0; c2 = y3 - y0; c3 = z3 - z0;

		gq1 = b2*c3 - b3*c2;
		gq2 = b3*c1 - b1*c3;
		gq3 = b1*c2 - b2*c1;
		vol = a1*gq1 + a2*gq2 + a3*gq3;
#define SIXTH 0.1666666667
		vol *= SIXTH;

		if( cp->c_dist > 0.0 && (viol = vol - cp->c_dist) > cutoff) {
			tviol += viol;
			fprintf( fp, "\tchirality: %5d %5d %5d %5d %5d %9.2f %9.2f %9.2f\n",
				i+1, cp->c_anum[0]+1, cp->c_anum[1]+1, cp->c_anum[2]+1,
				cp->c_anum[3]+1, vol, cp->c_dist, viol );
		}else if( cp->c_dist < 0.0 && (viol = cp->c_dist - vol) > cutoff) {
			tviol += viol;
			fprintf( fp, "\tchirality: %5d %5d %5d %5d %5d %9.2f %9.2f %9.2f\n",
				i+1, cp->c_anum[0]+1, cp->c_anum[1]+1, cp->c_anum[2]+1,
				cp->c_anum[3]+1, vol, cp->c_dist, viol );
		}
	}

	return( tviol );
}
