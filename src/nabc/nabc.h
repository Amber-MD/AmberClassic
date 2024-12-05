#ifndef	NABC_H
#define	NABC_H

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <sff.h>

	/* Fundamental nab types, that are not also in sff.h:	*/

	/* geometric types: Frames are not directly accessible 	*/
	/* at the nab level.					*/

typedef	REAL_T	POINT_T[ 3 ];		/* 0 == x, 1 == y, 2 == z */
typedef	REAL_T	FRAME_T[ 4 ][ 3 ];	/* org, x, y, z axes	*/
typedef	REAL_T	MATRIX_T[ 4 ][ 4 ];	
/* If MATRIX_T m; then REF_MATRIX_T is the type of the array name m. */
/* Functions returning, in principle, MATRIX_T should have a return value */
/* type of REF_MATRIX_T.  Note that MATRIX_T* is not REF_MATRIX_T. */
typedef	REAL_T	( *REF_MATRIX_T )[ 4 ];

	/* atom, residue & molecule types, includes types	*/
	/* INTBOND_T, EXTBOND_T, CHIRAL_T & STRAND_T that have	*/
	/* no nab equivalent					*/

#define	A_CONNECT_SIZE	8
typedef	struct	atom_t	{
	STRING_T *a_atomname;
	STRING_T *a_atomtype;
	INT_T	a_attr;
	INT_T	a_nconnect;
	INT_T	a_connect[ A_CONNECT_SIZE ];
	struct	residue_t	*a_residue;
	REAL_T	a_charge;
	REAL_T	a_radius;
	REAL_T	a_bfact;
	REAL_T	a_occ;
	STRING_T *a_element;	/* from SD (mol) files		*/
	INT_T	a_int1;		/* user property		*/
	REAL_T	a_float1;	/* user property		*/
	REAL_T	a_float2;	/* user property		*/
	INT_T	a_tatomnum;
	INT_T	a_atomnum;
	STRING_T *a_fullname;
	POINT_T	a_pos;
	REAL_T	a_w;		/* 4th dimension		*/
} ATOM_T;

typedef	struct	extbond_t	{
	struct	extbond_t	*eb_next;
	INT_T	eb_anum;	/* atom in current residue	*/
	INT_T	eb_rnum;	/* other residue number		*/
	INT_T	eb_ranum;	/* atom in other residue	*/
} EXTBOND_T;

typedef	INT_T	INTBOND_T[ 2 ];

typedef	struct	residue_t	{
	struct	residue_t	*r_next;
	INT_T	r_num;
	INT_T	r_tresnum;	/* set by NAB_rri( a, "tresnum" )  */
	INT_T	r_resnum;	/* set by NAB_rri( a, "resnum" )  */
	STRING_T *r_resname;	/* set by NAB_rrc( a, "resname" ) */
	STRING_T *r_resid;	/* set by NAB_rrc( a, "resid" ) */
	INT_T	r_attr;
	INT_T	r_kind;
	INT_T	r_atomkind;
	struct	strand_t	*r_strand;
	EXTBOND_T	*r_extbonds;	
	INT_T	r_nintbonds;	/* INTERNAL bonds		*/
	INTBOND_T	*r_intbonds;
	INT_T		r_nchiral;
	CHIRAL_T	*r_chiral;
	INT_T	r_natoms;
	INT_T	*r_aindex;
	ATOM_T	*r_atoms;
} RESIDUE_T;

typedef	struct	strand_t	{	/* not visible in nab 	*/
	STRING_T *s_strandname;
	INT_T	s_strandnum;
	INT_T	s_attr;
	struct	molecule_t	*s_molecule;
	struct	strand_t	*s_next;
	INT_T	s_nresidues;
	INT_T	s_res_size;
	RESIDUE_T	**s_residues;
} STRAND_T;

typedef	struct	molecule_t	{
	FRAME_T	m_frame;
	INT_T	m_nstrands;
	STRAND_T	*m_strands;
	INT_T	m_nresidues;
	INT_T	m_natoms;
	INT_T	m_nvalid;	/* "numbers" valid	*/
	PARMSTRUCT_T	*m_prm;	/* points forcefield stuff */
} MOLECULE_T;

	/* attributes of atoms, residues and molecules.	*/
	/* AT_SELECT is set by select_atoms() when it 	*/
	/* evaluates a regular expression.  AT_SELECTED	*/
	/* is used to store the previous select value	*/
	/* and is used in the construction of bounds	*/
	/* and dist matrices				*/
	
#define	AT_SELECT	0001
#define	AT_SELECTED	0002
#define	AT_WORK		0200

	/* Residue & Atom "kinds"	*/

#define	RT_UNDEF	0
#define	RT_DNA		1
#define	RT_RNA		2
#define	RT_AA		3

#define	RAT_UNDEF	0
#define	RAT_UNITED	1
#define	RAT_ALLATOM	2

	/* nab builtins (but no libc or libm calls):	*/

INT_T		addresidue( MOLECULE_T*, STRING_T*, RESIDUE_T* );
INT_T		addstrand( MOLECULE_T*, STRING_T* );
INT_T		alignframe( MOLECULE_T*, MOLECULE_T* );
INT_T		andbounds( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T*, REAL_T, REAL_T );
REAL_T		angle( MOLECULE_T*, STRING_T*, STRING_T*, STRING_T* );
REAL_T		anglep( POINT_T, POINT_T, POINT_T );
INT_T		axis2frame( MOLECULE_T*, POINT_T, POINT_T );
MOLECULE_T	*bdna( STRING_T** );
INT_T		bonded_atoms( ATOM_T*, ATOM_T** );
INT_T		cap( MOLECULE_T*, STRING_T*, INT_T, INT_T );
INT_T		circle( REAL_T*, REAL_T*, REAL_T*, REAL_T* );
INT_T		connectres( MOLECULE_T*, STRING_T*, INT_T, STRING_T*, INT_T, STRING_T* );
MOLECULE_T	*copymolecule( MOLECULE_T* );
INT_T		countmolatoms( MOLECULE_T*, STRING_T* );
INT_T		countmolres( MOLECULE_T*, STRING_T* );
INT_T		countmolstrands( MOLECULE_T*, STRING_T* );
INT_T		countstrandresidues( MOLECULE_T*, INT_T );
STRING_T	*date();
REAL_T		db_viol( REAL_T*, REAL_T*, INT_T* );
REAL_T		db_viol3( REAL_T*, REAL_T*, INT_T* );
MOLECULE_T	*dg_helix( STRING_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T**, 
		REAL_T*, REAL_T*, REAL_T*, REAL_T*, STRING_T** );
INT_T		dg_options( BOUNDS_T*, STRING_T* );
/* REAL_T		dist( MOLECULE_T*, STRING_T*, STRING_T* );  --duplicated in molio.c  */
REAL_T		distp( POINT_T, POINT_T );
BOUNDS_T	*dt_to_bmat( MOLECULE_T**, STRING_T**, STRING_T** );
INT_T		dt_to_prmtop( STRING_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T** );
INT_T		dumpbounds( FILE_T*, BOUNDS_T*, INT_T );
REAL_T		dumpboundsviolations( FILE_T*, BOUNDS_T*, REAL_T );
REAL_T		dumpchiviolations( FILE_T*, BOUNDS_T*, REAL_T);
INT_T		dumpmatrix( FILE_T*, MATRIX_T );
INT_T		dumpmolecule( FILE_T*, MOLECULE_T*, INT_T, INT_T, INT_T);
INT_T		embed( BOUNDS_T*, REAL_T* );
MOLECULE_T     *fd_helix( STRING_T*, STRING_T* );
INT_T		freemolecule( MOLECULE_T* );
INT_T		freeparm( MOLECULE_T* );
INT_T		freeresidue ( RESIDUE_T* );
#ifndef WIN32
STRING_T	*ftime( STRING_T* );
#endif
REAL_T      gauss( REAL_T*, REAL_T* );
INT_T		geodesics( BOUNDS_T* );
REAL_T		getchivol( MOLECULE_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T** );
REAL_T		getchivolp( POINT_T, POINT_T, POINT_T, POINT_T );
MOLECULE_T	*getcif( STRING_T*, STRING_T* );
MOLECULE_T	*getcompcif( STRING_T*, STRING_T* );
STRING_T	*NAB_getline( FILE_T* );
REF_MATRIX_T	getmatrix( STRING_T* );
INT_T		getseq_from_pdb( STRING_T**, INT_T*, STRING_T**, STRING_T**, STRING_T** );
MOLECULE_T	*getpdb( STRING_T*, STRING_T* );
MOLECULE_T	*getpdb_prm( STRING_T**, STRING_T**, STRING_T**, INT_T* );
RESIDUE_T	*getresidue( STRING_T* );
STRING_T	*getreslibkind( STRING_T* );
STRING_T	*getresname( RESIDUE_T* );
INT_T		getxyz_from_pdb( STRING_T**, MOLECULE_T**, STRING_T**, INT_T* );
INT_T		helixanal( MOLECULE_T** );
INT_T		length( STRING_T*** );
MOLECULE_T 	*linkprot( STRING_T**, STRING_T**, STRING_T** );
MOLECULE_T 	*link_na( STRING_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T** );
            
INT_T		mergestr( MOLECULE_T*, STRING_T*, STRING_T*, MOLECULE_T*, STRING_T*, STRING_T* );
INT_T		metrize( BOUNDS_T*, INT_T );
INT_T           mme2_timer();
INT_T		mme_init( MOLECULE_T*, STRING_T*, STRING_T*, REAL_T*, STRING_T* );
INT_T		mme_timer();
REAL_T		molsurf( MOLECULE_T**, STRING_T**, REAL_T* );
INT_T           mpierror( INT_T );
INT_T           mpifinalize( void );
INT_T           mpiinit( INT_T*, STRING_T**, INT_T*, INT_T* );
BOUNDS_T	*newbounds( MOLECULE_T*, STRING_T* );
MOLECULE_T	*newmolecule();
REF_MATRIX_T	newtransform( REAL_T, REAL_T, REAL_T, REAL_T, REAL_T, REAL_T );
INT_T           nm_timer();
INT_T		orbounds( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T*, REAL_T, REAL_T );
REAL_T		pair_ener( STRING_T*, INT_T );
INT_T		plane( MOLECULE_T**, STRING_T**, REAL_T*, REAL_T*, REAL_T* );
INT_T		putarc( STRING_T**, MOLECULE_T** );
INT_T		putbnd( STRING_T*, MOLECULE_T* );
INT_T		putcif( STRING_T*, STRING_T*, MOLECULE_T* );
INT_T		putdist( STRING_T*, MOLECULE_T*, STRING_T*, STRING_T* );
INT_T		putmatrix( STRING_T*, MATRIX_T );
INT_T		putpdb( STRING_T*, MOLECULE_T*, STRING_T* );
INT_T		putx( STRING_T**, MOLECULE_T** );
REAL_T		rand2( void );
INT_T		readparm( MOLECULE_T*, STRING_T* );
INT_T		rmsd( MOLECULE_T**, STRING_T**, MOLECULE_T**, STRING_T**, REAL_T*);
REF_MATRIX_T	rot4( MOLECULE_T*, STRING_T*, STRING_T*, REAL_T );
REF_MATRIX_T	rot4p( POINT_T, POINT_T, REAL_T );
INT_T       rseed( void );
FILE_T		*safe_fopen( STRING_T*, STRING_T* );
INT_T		setbounds( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T*, REAL_T, REAL_T );
INT_T		setboundsfromdb( BOUNDS_T**, MOLECULE_T**, STRING_T**, STRING_T**, 
		STRING_T**, REAL_T* );
INT_T		setchiplane( BOUNDS_T**, MOLECULE_T**, STRING_T** );
INT_T		setchivol( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T*, STRING_T*, STRING_T*, 
		REAL_T );
INT_T		setframe( INT_T, MOLECULE_T*, STRING_T*, STRING_T*, STRING_T*, STRING_T*, STRING_T* );
INT_T		setframep( INT_T, MOLECULE_T*, POINT_T, POINT_T, POINT_T, POINT_T, POINT_T );
INT_T		setmol_from_xyz( MOLECULE_T**, STRING_T**, REAL_T* );
INT_T		setmol_from_xyzw( MOLECULE_T**, STRING_T**, REAL_T* );
INT_T		setpoint( MOLECULE_T*, STRING_T*, POINT_T );
INT_T		setreskind( MOLECULE_T*, STRING_T*, STRING_T* );
INT_T		setreslibkind( STRING_T* );
INT_T       setseed( INT_T* );
INT_T		setxyz_from_mol( MOLECULE_T**, STRING_T**, POINT_T* );
INT_T		setxyzw_from_mol( MOLECULE_T**, STRING_T**, REAL_T* );
INT_T		showbounds( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T* );
INT_T		split( STRING_T*, STRING_T**, STRING_T* );
REAL_T		step_ener( STRING_T*, INT_T );
STRING_T	*substr( STRING_T*, INT_T, INT_T );
INT_T		sugarpuckeranal( MOLECULE_T**, INT_T*, INT_T*, INT_T* );
REF_MATRIX_T		superimpose( MOLECULE_T*, STRING_T*, MOLECULE_T*, STRING_T* );
STRING_T	*timeofday();
REF_MATRIX_T	trans4( MOLECULE_T*, STRING_T*, STRING_T*, REAL_T );
REF_MATRIX_T	trans4p( POINT_T, POINT_T, REAL_T );
REAL_T		torsion( MOLECULE_T*, STRING_T*, STRING_T*, STRING_T*, STRING_T* );
REAL_T		torsionp( POINT_T, POINT_T, POINT_T, POINT_T );
INT_T		transformmol( MATRIX_T, MOLECULE_T*, STRING_T* );
INT_T		transformpts( MATRIX_T, POINT_T*, INT_T );
RESIDUE_T	*transformres( MATRIX_T, RESIDUE_T*, STRING_T* );
INT_T		tsmooth( BOUNDS_T*, REAL_T );
INT_T		useboundsfrom( BOUNDS_T*, MOLECULE_T*, STRING_T*, MOLECULE_T*, STRING_T*, REAL_T );
INT_T		usemodeldist( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T* );
REF_MATRIX_T	updtransform( MATRIX_T, MATRIX_T );
MOLECULE_T	*wc_basepair( RESIDUE_T*, RESIDUE_T* );
STRING_T	*wc_complement( STRING_T*, STRING_T* );
MOLECULE_T	*wc_helix( STRING_T*, STRING_T*, STRING_T*, STRING_T*, 
                       REAL_T, REAL_T, REAL_T, REAL_T, STRING_T* );
INT_T		writeparm( MOLECULE_T*, STRING_T* );

void  copy_mat( MATRIX_T, MATRIX_T );

/*  from stringutil.c:  */
char    *substr( char [], int, int );
char    *NAB_getline( FILE * );
int split( char [], char *[], char * );
int split_n( char [], int, char *[], char * );
int NAB_index( char [], char [] );
char    *NAB_strcat( char [], char [] );
char    *NAB_strcpy( char **, char [] );
int NAB_strcmp( char *, char * );
char    *NAB_readstring( char ** );
int NAB_newstring( char **, int * );

REF_MATRIX_T    NAB_matcpy( REF_MATRIX_T, REF_MATRIX_T );

/* Functions defined in matop.c */

int	MAT_fprint( FILE *, int, MATRIX_T [] );
int	MAT_sprint( char [], int, MATRIX_T [] );
int	MAT_fscan( FILE *, int, MATRIX_T [] );
int	MAT_sscan( char [], int, MATRIX_T [] );
REF_MATRIX_T	MAT_concat( MATRIX_T, MATRIX_T );
int	MAT_count( char [] );
char	*MAT_getsyminfo( void );
int	MAT_istrue( MATRIX_T );

	/* functions for accessing parts of atoms & residues */

INT_T		*NAB_ari( ATOM_T *ap, char *key );
REAL_T		*NAB_arf( ATOM_T *ap, char *key );
POINT_T		*NAB_arp( ATOM_T *ap, char key[] );
INT_T		*NAB_rri( RESIDUE_T *res, char *key );
STRING_T	**NAB_rrc( RESIDUE_T *res, char *key );
INT_T		*NAB_mri( MOLECULE_T *mol, char *key );

	/* functions for for( a in m ) etc	*/

ATOM_T		*NAB_mnext (MOLECULE_T *mol, ATOM_T *cap);
ATOM_T		*NAB_anext( RESIDUE_T *res, ATOM_T *cap );
RESIDUE_T	*NAB_rnext( MOLECULE_T *mol, RESIDUE_T *crp );

    /* trig functions in degrees:   */

#define R2D 57.29577951308232090712
#define D2R  0.01745329251994329576
#define ACOS(c) (R2D*acos(c))
#define ASIN(s) (R2D*asin(s))
#define ATAN(t) (R2D*atan(t))
#define ATAN2(y,x)  (R2D*atan2((y),(x)))
#define COS(a)  cos(D2R*(a))
#define SIN(a)  sin(D2R*(a))
#define TAN(a)  tan(D2R*(a))

#define	NAB_RSBUF_SIZE	10000
#endif
