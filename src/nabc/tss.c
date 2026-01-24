#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "nabcode.h"

	/* Compile phase stuff: */

#define	T_INT		0
#define	T_FLT		1
#define	T_STR		2


#define	S_PARMNAME_SIZE	20
#define	S_PARMVALUE_SIZE	60
typedef	struct	symparm_t	{
	struct	symparm_t	*s_next;
	int		s_addr;
	int		s_type;
	char		s_name[ S_PARMNAME_SIZE ];
	char		s_value[ S_PARMVALUE_SIZE ];
} SYMPARM_T;

#define	S_OP_SIZE	16
#define	S_OPNAME_SIZE	32
#define S_MAX_INNER	4
typedef	struct	symop_t	{
	char	s_op[ S_OP_SIZE ];
	int	s_id;
	char	s_name[ S_OPNAME_SIZE ];
	int	s_addr;
	int	s_active;
	int	s_isloop;
	SYMPARM_T	*s_plist;
	int	s_n_inner;
	struct	symop_t *s_inner[ S_MAX_INNER ];
} SYMOP_T;

#define	MAXSYMOPS 30
static	SYMOP_T	*symops[ MAXSYMOPS ];
static	int	n_symop;

static	int	parse();
static	int	compile();
static	int	set_vars();
static	void	dump_vars();
static	void	dump_prog();
static	int	findvar();
static	int	getop();

	/* Compile and runtime stuff:	*/

#define	OP_CYC_2	0
#define	OP_CYC_3	1
#define	OP_CYC_4	2
#define	OP_CYC_5	3
#define	OP_CYC_6	4
#define	OP_CYC_N	5
#define	OP_CUBE		6
#define	OP_DIHED	7
#define	OP_HELIX	8
#define	OP_ICO		9
#define	OP_OCTA		10
#define	OP_ORIENT	11
#define	OP_ROTATE	12
#define	OP_TETRA	13
#define	OP_TRANS	14
#define	OP_MERGE	15
#define	OP_MULT		16
#define	OP_READ		17

#define	I_PARM_SIZE	10
typedef	struct	inst_t	{
	int	i_op;
	int	i_parms[ I_PARM_SIZE ];
} INST_T;

INST_T	SG_prog[ MAXSYMOPS ];
int	SG_s_prog;
int	SG_pc;

	/* Run time stuff: */

#define W_CUR 0
#define W_LO  1
#define W_HI  2
#define W_INC 3

typedef	struct	mat_t	{
	int	nmats;
	MATRIX_T	*mats;
} MAT_T;

static	MAT_T	*SG_mats = NULL;
static	int	SG_n_mats;

static int	SG_ovfl = 0;
static int	SG_nact = 0;
static int	*SG_act = NULL;
static int	SG_nvar = 0;
static double	*SG_val = NULL;
static double	*SG_lo = NULL;
static double	*SG_hi = NULL;
static double	*SG_inc = NULL;

static int	SG_nstr = 0;
static char	**SG_str = NULL;

static	int	SG_cyclic_2();
static	int 	SG_cyclic_3();
static	int 	SG_cyclic_4();
static	int 	SG_cyclic_5();
static	int 	SG_cyclic_6();
static	int 	SG_cyclic_N();
static	int 	SG_cube();
static	int 	SG_dihedral();
static	int 	SG_dodeca();
static	int 	SG_helix();
static	int 	SG_ico();
static	int 	SG_octa();
static	int 	SG_orient();
static	int	SG_polyhedron();
static	int 	SG_rotate();
static	int 	SG_tetra();
static	int 	SG_translate();

static	int 	SG_merge();
static	int 	SG_multiply();
static	int 	SG_read();

#define	MAT_ALLOC(n)	( ( MATRIX_T * )malloc( ( n )*sizeof( MATRIX_T ) ) )
#define	BPN	14

static	int	mk_omats();

int	tss_init( matf, rangef )
char	matf[];
char	rangef[];
{
	FILE	*mfp, *rfp;
	int	rval;

	if( ( mfp = fopen( matf, "r" ) ) == NULL ){
		fprintf( stderr, "syminit: can't open mat-file %s\n",
			matf );
		return( 1 );
	}
	rval = parse( mfp, stdin );
	fclose( mfp );
	if( rval != 0 )
		return( 1 );

	if( ( rfp = fopen( rangef, "r" ) ) == NULL ){
		fclose( mfp );
		fprintf( stderr, "syminit: can't open range-file %s\n",
			rangef );
		return( 1 );
	}
	rval = compile( rfp, stdout );
	fclose( rfp );

	return( rval );
}

static	int	parse( fp )
FILE	*fp;
{
	char	line[ 256 ];
	char	*lp;
	SYMOP_T	*symp, *symp1;
	int	s, s1, si, nlp, id, dup;
	SYMOP_T	*s_stk[ 20 ];
	int	s_stkp;
	SYMPARM_T	*sparm, *sp, *sp1;
	char	jstr[ 256 ];
	int	junk;

	/* read the mat file, collect parms:	*/

	s_stkp = -1;
	symp = NULL;
	for( nlp = 0, n_symop = 0; fgets( line, sizeof( line ), fp ); ){
		if( *line != '#' )
			continue;
		else if( !strncmp( line, "#S{", 3 ) ){
			s_stkp++;
			s_stk[ s_stkp ] = symp;
			symp = ( SYMOP_T * )malloc( sizeof( SYMOP_T ) );
			sscanf( &line[ 4 ], "%s %d", symp->s_op, &symp->s_id );
			if( strcmp( symp->s_op, "merge" ) &&
				strcmp( symp->s_op, "multiply" ) ){
				sscanf( &line[ 4 ], "%s %d %s",
					jstr, &junk, symp->s_name );
			}else
				*symp->s_name = '\0';
			symp->s_active = 1;
			symp->s_isloop = 0;
			symp->s_plist = NULL;
			if( !strcmp( symp->s_op, "merge" ) ||
				!strcmp( symp->s_op, "multiply" ) ||
				!strcmp( symp->s_op, "read" ) )
				symp->s_isloop = 0;
			else{
				symp->s_isloop = 1;
				nlp++;
			}
			symp->s_n_inner = 0;
			for( s = 0; s < 4; s++ )
				symp->s_inner[ s ] = NULL;
			symops[ n_symop ] = symp;
			if( s_stk[ s_stkp ] ){
				s_stk[s_stkp]->
					s_inner[s_stk[s_stkp]->s_n_inner]=symp;
				s_stk[s_stkp]->s_n_inner++;
			}
			n_symop++;
		}else if( !strncmp( line, "#S+", 3 ) ){
			sparm = ( SYMPARM_T * )malloc( sizeof( SYMPARM_T ) );
			sparm->s_next = NULL;
			lp = &line[ 3 ];
			sscanf( &lp[ 3 ], "%s", sparm->s_name );
			lp += strspn( lp, " \t" ); 		
			lp = strpbrk( lp, " \t" );
			lp += strspn( lp, " \t" ); 		
			sscanf( lp, "%[^\n]", sparm->s_value );
			sp1 = NULL;
			for( sp = symp->s_plist; sp; sp = sp->s_next )
				sp1 = sp;
			if( sp1 == NULL )
				symp->s_plist = sparm;
			else
				sp1->s_next = sparm;
		}else if( !strncmp( line, "#S}", 3 ) ){
			symp = s_stk[ s_stkp ];
			s_stkp--;
		}else{
			fprintf( stderr,
				"parse: illegal symmetry directive:\n" );
			fprintf( stderr, "\t%s", line );
			return( 1 );
		}
	}

	/* remove duplicated code: */

	for( s = 0; s < n_symop - 1; s++ ){
		id = symops[ s ]->s_id;
		for( s1 = s + 1; s1 < n_symop; s1++ ){
			if( symops[ s1 ]->s_id == id ){
				symops[ s ]->s_active = 0;
				break;
			}
		}
	}
	for( s = 0; s < n_symop - 1; s++ ){
		if( !symops[ s ]->s_active )
			continue;
		for( si = 0; si < symops[ s ]->s_n_inner; si++ ){
			symp = symops[ s ]->s_inner[ si ];
			id = symp->s_id;
			for( s1 = s + 1; s1 < n_symop; s1++ ){
				if( !symops[ s1 ]->s_active )
					continue;
				if( symops[ s1 ]->s_id == id ){
					symops[s]->s_inner[si] = symops[s1];
					break;
				}
			}
		}
	}

	/* check for duplicate symop names, error if any found: */

	for( dup = 0, s = n_symop - 1; s > 0; s-- ){
		symp = symops[ s ];
		if( !*symp->s_name )
			continue;
		for( s1 = s - 1; s1 >= 0; s1-- ){
			symp1 = symops[ s1 ];
			if( !*symp1->s_name )
				continue;
			if( !strcmp( symp->s_name, symp1->s_name ) )
				dup = 1;
		}
	}
	if( dup ){
		fprintf( stderr,
		"parse: Matrix file contains duplicate symop names.\n" );
		for( s = n_symop - 1; s >= 0; s-- ){
			symp = symops[ s ];
			if( !*symp1->s_name )
				continue;
			fprintf( stderr, "symop[%2d] = %s\n", s, symp->s_name );
		}
		return( 1 );
	}

	return( 0 );
}

static	int	compile( rfp, ofp )
FILE	*rfp;
FILE	*ofp;
{
	int	s, si, pa;
	SYMOP_T	*symp, *i_symp;
	SYMPARM_T	*sp;

	/* free previously assigned matrices: */

	if( SG_mats != NULL ){
		for( s = 0; s < SG_n_mats; s++ )
			free( SG_mats[ s ].mats );
		free( SG_mats );
		SG_mats = NULL;
	}

	/* free previously assigned state vars: */

	if( SG_act != NULL ){
		free( SG_act );
		SG_nact = 0;
	}
	if( SG_val != NULL ){
		free( SG_val );
		free( SG_lo );
		free( SG_hi );
		free( SG_inc );
		free( SG_str );
	}

	/* assign addresses, compile prog: */

	SG_n_mats = 1; 	/* mats[ 0 ] == NULL */
	for( SG_s_prog = SG_nstr = SG_nvar = 0, s = n_symop - 1; s >= 0; s-- ){
		symp = symops[ s ];
		if( !symp->s_active )
			continue;
		symp->s_addr = SG_n_mats;
		SG_n_mats++;
		SG_prog[ SG_s_prog ].i_op = getop( symp->s_op );
		for( pa = 0; pa < I_PARM_SIZE; pa++ )
			SG_prog[ SG_s_prog ].i_parms[ pa ] = -1;
		SG_prog[ SG_s_prog ].i_parms[ 0 ] = symp->s_addr;
		if( !strcmp( symp->s_op, "read" ) ){
			sp = symp->s_plist;
			SG_prog[ SG_s_prog ].i_parms[ 1 ] = sp->s_addr;
		}else if( !strcmp( symp->s_op, "merge" ) ){
			for( si = 0; si < S_MAX_INNER; si++ ){
				i_symp = symp->s_inner[ si ];
				if( i_symp )
					SG_prog[SG_s_prog].i_parms[1+si] =
						i_symp->s_addr;
				else
					SG_prog[SG_s_prog].i_parms[1+si] = 0;
			}
		}else if( !strcmp( symp->s_op, "multiply" ) ){
			for( si = 0; si < 2; si++ ){
				i_symp = symp->s_inner[ si ];
				SG_prog[SG_s_prog].i_parms[1+si] =
					i_symp->s_addr;
			}
		}else{
			i_symp = symp->s_inner[ 0 ];
			if( i_symp )
				SG_prog[SG_s_prog].i_parms[1] = i_symp->s_addr;
			else
				SG_prog[ SG_s_prog ].i_parms[ 1 ] = 0;
			for( pa = 2, sp = symp->s_plist; sp; sp = sp->s_next ){
				if( !strcmp( sp->s_name, "noid" ) ){
					sp->s_type = T_STR;
					sp->s_addr = SG_nstr;
					SG_prog[SG_s_prog].i_parms[pa]=SG_nstr;
					pa++;
					SG_nstr++;
				}else if( !strcmp( sp->s_name, "axestype" ) ){
					sp->s_type = T_STR;
					sp->s_addr = SG_nstr;
					SG_prog[SG_s_prog].i_parms[pa]=SG_nstr;
					pa++;
					SG_nstr++;
				}else if( !strcmp( sp->s_name, "center" ) ){
					sp->s_type = T_FLT;
					sp->s_addr = SG_nvar;
					SG_prog[SG_s_prog].i_parms[pa]=SG_nvar;
					pa++;
					SG_nvar += 3;
				}else if( !strncmp( sp->s_name, "axis", 4 ) ){
					sp->s_type = T_FLT;
					sp->s_addr = SG_nvar;
					SG_prog[SG_s_prog].i_parms[pa]=SG_nvar;
					pa++;
					SG_nvar += 3;
				}else if( !strncmp( sp->s_name, "angle", 5 ) ){
					sp->s_type = T_FLT;
					sp->s_addr = SG_nvar;
					SG_prog[SG_s_prog].i_parms[pa]=SG_nvar;
					pa++;
					SG_nvar++;
				}else if( !strcmp( sp->s_name, "dist" ) ){
					sp->s_type = T_FLT;
					sp->s_addr = SG_nvar;
					SG_prog[SG_s_prog].i_parms[pa]=SG_nvar;
					pa++;
					SG_nvar++;
				}else if( !strcmp( sp->s_name, "count" ) ){
					sp->s_type = T_INT;
					sp->s_addr = SG_nvar;
					SG_prog[SG_s_prog].i_parms[pa]=SG_nvar;
					pa++;
					SG_nvar++;
				}else if( strcmp( sp->s_name, "file" ) ){
					fprintf( stderr,
						"Unknown var type %s\n",
						sp->s_name );
				}
			}
		}
		SG_s_prog++;
	}

	/* reallocate mats: */

	SG_mats = ( MAT_T * )malloc( SG_n_mats * sizeof( MAT_T ) );
	for( s = 0; s < SG_n_mats; s++ ){
		SG_mats[ s ].nmats = 0;
		SG_mats[ s ].mats = NULL;
	}

	/* reallocate the state vars: */

	SG_act = ( int * )malloc( SG_nvar * sizeof( int ) );
	SG_val = ( double * )malloc( SG_nvar * sizeof( double ) );
	SG_lo = ( double * )malloc( SG_nvar * sizeof( double ) );
	SG_hi = ( double * )malloc( SG_nvar * sizeof( double ) );
	SG_inc = ( double * )malloc( SG_nvar * sizeof( double ) );
	SG_str = ( char ** )malloc( SG_nstr * sizeof( char * ) );

	if( set_vars( rfp ) )
		return( 1 );

/*
	dump_vars( stdout );
	dump_prog( stdout );
*/

	return( 0 );
}

static	int	set_vars( rfp )
FILE	*rfp;
{
	int	a, s, cnt;
	SYMOP_T	*symp;
	SYMPARM_T	*sp;
	double	x, y, z, ang;
	char	line[ 256 ];
	char	opname[ 32 ], pname[ 32 ];
	double	lv, hv, iv;
	int	addr, type;

	for( s = 0; s < SG_nvar; s++ ){
		SG_val[ s ] = 0.0;
		SG_lo [ s ] = 0.0;
		SG_hi [ s ] = 0.0;
		SG_inc[ s ] = 0.0;
	}
	for( s = 0; s < SG_nstr; s++ )
		SG_str[ s ] = NULL;

	for( s = n_symop - 1; s >= 0; s-- ){
		symp = symops[ s ];
		if( !symp->s_active )
			continue;
		if( ( sp = symp->s_plist ) == NULL )
			continue;
		for( ; sp; sp = sp->s_next ){
			if( !strcmp( sp->s_name, "noid" ) ){

				SG_str[sp->s_addr] = sp->s_value;

			}else if( !strcmp( sp->s_name, "axestype" ) ){

				SG_str[sp->s_addr] = sp->s_value;

			}else if( !strcmp( sp->s_name, "center" ) ||
				!strncmp( sp->s_name, "axis", 4 ) ){

				sscanf( sp->s_value, "%lf %lf %lf", &x, &y, &z );
				SG_val[sp->s_addr+0] = x;
				SG_val[sp->s_addr+1] = y;
				SG_val[sp->s_addr+2] = z;


			}else if( !strncmp( sp->s_name, "angle", 5 ) ||
				!strcmp( sp->s_name, "dist" ) ){

				sscanf( sp->s_value, "%lf", &ang );
				SG_val[sp->s_addr] = ang;

			}else if( !strcmp( sp->s_name, "count" ) ){

				sscanf( sp->s_value, "%d", &cnt );
				SG_val[sp->s_addr] = cnt;

			}
		}
	}

	for( SG_nact = 0; fgets( line, sizeof( line ), rfp ); ){
		if( *line == '#' )
			continue;
		sscanf( line, "%s %s %lf %lf %lf", opname, pname, &lv, &hv, &iv );
		if( findvar( opname, pname, &addr, &type ) ){
			SG_act[SG_nact] = addr;
			SG_lo [SG_act[SG_nact]] = lv;
			SG_hi [SG_act[SG_nact]] = hv;
			SG_inc[SG_act[SG_nact]] = iv;
			SG_nact++;
		}else{
			fprintf( stderr, "set_vars: op:parm %s:%s not found.\n",
				opname, pname );
			return( 1 );
		}
	}

	for( a = 0; a < SG_nact; a++ )
		SG_val[SG_act[a]] = SG_lo[SG_act[a]];

	return( 0 );
}

static	void	dump_vars( fp )
FILE	*fp;
{
	int	a, v, s;
	
	fprintf( fp, "%d active vars\n", SG_nact );
	for( a = 0; a < SG_nact; a++ ){
		fprintf( fp, "SG_act[%3d] = %3d:\n", a, SG_act[a] );
		fprintf( fp, "\tSG_val[%3d] = %lf\n",
			SG_act[a], SG_val[SG_act[a]] );
		fprintf( fp, "\tSG_lo [%3d] = %lf\n", SG_act[a],
			SG_lo [SG_act[a]] );
		fprintf( fp, "\tSG_hi [%3d] = %lf\n",
			SG_act[a], SG_hi [SG_act[a]] );
		fprintf( fp, "\tSG_inc[%3d] = %lf\n",
			SG_act[a], SG_inc[SG_act[a]] );
	}
	for( v = 0; v < SG_nvar; v++ )
		fprintf( fp, "\tSG_val[%3d] = %lf\n", v, SG_val[v] );
	for( s = 0; s < SG_nstr; s++ )
		fprintf( fp, "\tSG_str[%3d] = %s\n", s, SG_str[s] );
}

static	void	dump_prog( fp )
FILE	*fp;
{
	int	p, i;
	
	fprintf( fp, "%d inst\n", SG_s_prog );
	for( p = 0; p < SG_s_prog; p++ ){
		fprintf( fp, "SG_prog[%3d] = %3d:", p, SG_prog[p].i_op );
		for( i = 0; i < I_PARM_SIZE; i++ )
			fprintf( fp, " %4d", SG_prog[p].i_parms[i] );
		fprintf( fp, "\n" );
	}
}

static	int	findvar( opname, pname, addr, type )
char	opname[];
char	pname[];
int	*addr;
int	*type;
{
	int	s;
	SYMOP_T	*symp;
	SYMPARM_T	*sp;

	for( s = n_symop - 1; s >= 0; s-- ){
		symp = symops[ s ];
		if( !symp->s_active )
			continue;
		if( strcmp( opname, symp->s_name ) )
			continue;
		if( ( sp = symp->s_plist ) == NULL )
			continue;
		for( ; sp; sp = sp->s_next ){
			if( !strcmp( pname, sp->s_name ) ){
				*addr = sp->s_addr;
				*type = sp->s_type;
				return( 1 );
			}
		}
	}
	return( 0 );
}

int	tss_read( fp )
FILE	*fp;
{
	char	*lp, line[ 256 ], *wp, word[ 256 ];
	int	i, j;

	while( fgets( line, sizeof( line ), fp ) ){
		if( *line != '#' )
			break;
	}
	sscanf( line, "%d", &SG_s_prog );
	fgets( line, sizeof( line ), fp );
	sscanf( line, "%d", &SG_pc );
	for( i = 0; i < SG_s_prog; i++ ){
		fgets( line, sizeof( line ), fp );
		lp = line + strspn( line, " \t" );
		SG_prog[ i ].i_op = atoi( lp );
		for( j = 0; j < I_PARM_SIZE; j++ ){
			lp = strpbrk( lp, " \t" );
			lp += strspn( lp, " \t" );
			SG_prog[i].i_parms[j] = atoi( lp );
		}
	}

	fgets( line, sizeof( line ), fp );
	sscanf( line, "%d", &SG_n_mats );
	SG_mats = ( MAT_T * )malloc( SG_n_mats * sizeof( MAT_T ) );
	for( i = 0; i < SG_n_mats; i++ ){
		SG_mats[ i ].nmats = 0;
		SG_mats[ i ].mats = NULL;
	}

	fgets( line, sizeof( line ), fp );
	sscanf( line, "%d", &SG_ovfl );
	fgets( line, sizeof( line ), fp );
	sscanf( line, "%d", &SG_nact );
	SG_act = ( int * )malloc( SG_nact * sizeof( int ) );
	for( i = 0; i < SG_nact; i++ ){
		fgets( line, sizeof( line ), fp );
		sscanf( line, "%d", &SG_act[i] );
	}

	fgets( line, sizeof( line ), fp );
	sscanf( line, "%d", &SG_nvar );
	SG_val = ( double * )malloc( SG_nvar * sizeof( double ) );
	SG_lo  = ( double * )malloc( SG_nvar * sizeof( double ) );
	SG_hi  = ( double * )malloc( SG_nvar * sizeof( double ) );
	SG_inc = ( double * )malloc( SG_nvar * sizeof( double ) );
	for( i = 0; i < SG_nvar; i++ ){
		fgets( line, sizeof( line ), fp );
		sscanf( line, "%lf %lf %lf %lf",
			&SG_val[i], &SG_lo [i], &SG_hi [i], &SG_inc[i] );
	}

	fgets( line, sizeof( line ), fp );
	sscanf( line, "%d", &SG_nstr );
	SG_str = ( char ** )malloc( SG_nstr * sizeof( char * ) );
	for( i = 0; i < SG_nstr; i++ ){
		fgets( line, sizeof( line ), fp );
		sscanf( line, "%s", word );
		wp = ( char * )malloc( ( strlen( word )+1 ) * sizeof( char ) );
		strcpy( wp, word );
		SG_str[ i ] = wp;
	}

	return( 0 );
}

int	tss_write( fp )
FILE	*fp;
{
	int	i, j;
	int	*parms;
	
	fprintf( fp, "# sym program\n" );
	fprintf( fp, "%d inst\n", SG_s_prog );
	fprintf( fp, "%d pc\n", SG_pc );
	for( i = 0; i < SG_s_prog; i++ ){
		fprintf( fp, "%3d", SG_prog[ i ].i_op );
		parms = SG_prog[ i ].i_parms;
		for( j = 0; j < I_PARM_SIZE; j++ )
			fprintf( fp, " %4d", parms[ j ] );
		fprintf( fp, "\n" );
	}

	fprintf( fp, "%d matrices\n", SG_n_mats );
	fprintf( fp, "%d SG_ovfl\n", SG_ovfl );
	fprintf( fp, "%d SG_nact\n", SG_nact );
	for( i = 0; i < SG_nact; i++ ){
		fprintf( fp, "%d\n", SG_act[i] );
	}
	fprintf( fp, "%d SG_nvar\n", SG_nvar );
	for( i = 0; i < SG_nvar; i++ ){
		fprintf( fp, "%lg %lg %lg %lg\n",
			SG_val[i], SG_lo[i], SG_hi[i], SG_inc[i] );
	}
	fprintf( fp, "%d SG_nstr\n", SG_nstr );
	for( i = 0; i < SG_nstr; i++ )
		fprintf( fp, "%s\n", SG_str[i] );

	return( 0 );
}

int	tss_next( err, somats, nomats, omats, scur, ncur, cur )
int	*err;
int	*somats;
int	*nomats;
MATRIX_T omats[];
int	*scur;
int	*ncur;
double	cur[];
{
	int	a, carry;
	int	m;
	int	*parms;

	*err = 0;
	if( SG_ovfl )
		return( 0 );

	for( SG_pc = 0; SG_pc < SG_s_prog; SG_pc++ ){
		parms = SG_prog[ SG_pc ].i_parms;
		switch( SG_prog[ SG_pc ].i_op ){
		case OP_CYC_2 :
			SG_cyclic_2(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ] );
			break;
		case OP_CYC_3 :
			SG_cyclic_3(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ] );
			break;
		case OP_CYC_4 :
			SG_cyclic_4(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ] );
			break;
		case OP_CYC_5 :
			SG_cyclic_5(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ] );
			break;
		case OP_CYC_6 :
			SG_cyclic_6(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ] );
			break;
		case OP_CYC_N :
			SG_cyclic_N(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ],
				&SG_val [ parms[ 6 ] ],
				&SG_val [ parms[ 7 ] ] );
			break;
		case OP_CUBE :
			SG_cube(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ],
				&SG_val [ parms[ 6 ] ] );
			break;
		case OP_DIHED :
			SG_dihedral(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ],
				&SG_val [ parms[ 6 ] ],
				&SG_val [ parms[ 7 ] ] );
			break;
		case OP_HELIX :
			SG_helix(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ],
				&SG_val [ parms[ 6 ] ],
				&SG_val [ parms[ 7 ] ],
				&SG_val [ parms[ 8 ] ] );
			break;
		case OP_ICO :
			SG_ico(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ],
				&SG_val [ parms[ 6 ] ] );
			break;
		case OP_OCTA :
			SG_octa(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ],
				&SG_val [ parms[ 6 ] ] );
			break;
		case OP_ORIENT :
			SG_orient(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				&SG_val [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ],
				&SG_val [ parms[ 6 ] ],
				&SG_val [ parms[ 7 ] ],
				&SG_val [ parms[ 8 ] ],
				&SG_val [ parms[ 9 ] ] );
			break;
		case OP_ROTATE :
			SG_rotate(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				&SG_val [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ] );
			break;
		case OP_TETRA :
			SG_tetra(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				SG_str  [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ],
				&SG_val [ parms[ 6 ] ] );
			break;
		case OP_TRANS :
			SG_translate(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				SG_str  [ parms[ 2 ] ],
				&SG_val [ parms[ 3 ] ],
				&SG_val [ parms[ 4 ] ],
				&SG_val [ parms[ 5 ] ] );
			break;
		case OP_MERGE :
			SG_merge(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				&SG_mats[ parms[ 2 ] ],
				&SG_mats[ parms[ 3 ] ],
				&SG_mats[ parms[ 4 ] ] );
			break;
		case OP_MULT :
			SG_merge(
				&SG_mats[ parms[ 0 ] ],
				&SG_mats[ parms[ 1 ] ],
				&SG_mats[ parms[ 2 ] ] );
			break;
		}
	}

	*nomats = SG_mats[ SG_n_mats - 1 ].nmats;
	if( *nomats > *somats ){
		*err = 1;
		return( 0 );
	}

	for( m = 0; m < *nomats; m++ ){
		NAB_matcpy( omats[ m ], SG_mats[ SG_n_mats - 1 ].mats[ m ] );
	}

	*ncur = SG_nact;
	if( SG_nact > *scur ){
		*err = 2;
		return( 0 );
	}
	for( a = 0; a < SG_nact; a++ )
		cur[a] = SG_val[SG_act[a]];

	for( carry = 1, a = 0; carry && a < SG_nact; a++ ){
		if( SG_val[SG_act[a]] >= SG_hi[SG_act[a]] ){
			carry = 1;
			SG_val[SG_act[a]] = SG_lo[SG_act[a]];
		}else{
			carry = 0;
			SG_val[SG_act[a]] += SG_inc[SG_act[a]];
		}
	}
	SG_ovfl = carry;

	return( 1 );
}

int	tss_get( what, scur, ncur, cur )
int	*what;
int	*scur;
int	*ncur;
double	cur[];
{
	int a;

	if( SG_nact > *scur )
        	return( 2 );

	*ncur = SG_nact;
	switch( *what ){
	case W_CUR :
		for( a = 0; a < SG_nact; a++ )
			cur[a] = SG_val[SG_act[a]];
		break;
	case W_LO :
		for( a = 0; a < SG_nact; a++ )
			cur[a] = SG_lo[SG_act[a]];
		break;
	case W_HI :
		for( a = 0; a < SG_nact; a++ )
			cur[a] = SG_hi[SG_act[a]];
		break;
	case W_INC :
		for( a = 0; a < SG_nact; a++ )
			cur[a] = SG_inc[SG_act[a]];
		break;
	default :
		return( 1 );
	}
	return( 0 );
}

int	tss_set( what, ncur, cur )
int	*what;
int	*ncur;
double	cur[];
{
	int	a;

	if( SG_nact != *ncur ) 
		return( 2 );

	switch( *what ){
	case W_CUR :
		for( a = 0; a < SG_nact; a++ )
			SG_val[SG_act[a]] = cur[a];
		SG_ovfl = 0;
		break;
	case W_LO :
		for( a = 0; a < SG_nact; a++ )
			SG_lo[SG_act[a]] = cur[a];
		break;
	case W_HI:
		for( a = 0; a < SG_nact; a++ )
			SG_hi[SG_act[a]] = cur[a];
		break;
	case W_INC :
		for( a = 0; a < SG_nact; a++ )
			SG_inc[SG_act[a]] = cur[a];
		break;
	default :
		return( 1 );
	}
	return( 0 );
}

static	int	getop( s_op )
char	s_op[];
{

	if( !strcmp( s_op, "cyclic_2" ) )
		return( OP_CYC_2 );
	else if( !strcmp( s_op, "cyclic_3" ) )
		return( OP_CYC_3 );
	else if( !strcmp( s_op, "cyclic_4" ) )
		return( OP_CYC_4 );
	else if( !strcmp( s_op, "cyclic_5" ) )
		return( OP_CYC_5 );
	else if( !strcmp( s_op, "cyclic_6" ) )
		return( OP_CYC_6 );
	else if( !strcmp( s_op, "cyclic_N" ) )
		return( OP_CYC_N );
	else if( !strcmp( s_op, "cube" ) )
		return( OP_CUBE );
	else if( !strcmp( s_op, "dihedral" ) )
		return( OP_DIHED );
	else if( !strcmp( s_op, "helix" ) )
		return( OP_HELIX );
	else if( !strcmp( s_op, "ico" ) )
		return( OP_ICO );
	else if( !strcmp( s_op, "octa" ) )
		return( OP_OCTA );
	else if( !strcmp( s_op, "orient" ) )
		return( OP_ORIENT );
	else if( !strcmp( s_op, "rotate" ) )
		return( OP_ROTATE );
	else if( !strcmp( s_op, "tetra" ) )
		return( OP_TETRA );
	else if( !strcmp( s_op, "translate" ) )
		return( OP_TRANS );
	else if( !strcmp( s_op, "merge" ) )
		return( OP_MERGE );
	else if( !strcmp( s_op, "multiply" ) )
		return( OP_MULT );
	else if( !strcmp( s_op, "read" ) )
		return( OP_READ );
	else
		fprintf( stderr, "getop: Unknown operator: '%s'.\n",
			s_op );
	return( -1 );
}

static	int SG_cyclic_2( p_omat, p_imat, p_noid, p_axestype, p_center, p_axes )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes;
{
	double	p_angle;
	double	p_dist;
	double	p_count;

	p_angle = 180.0;
	p_dist = 0.0;
	p_count = 2;
	return( SG_helix( p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes,
		&p_angle, &p_dist, &p_count ) ); 
}

static	int SG_cyclic_3( p_omat, p_imat, p_noid, p_axestype, p_center, p_axes )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes;
{
	double	p_angle;
	double	p_dist;
	double	p_count;

	p_angle = 120.0;
	p_dist = 0.0;
	p_count = 3;
	return( SG_helix( p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes,
		&p_angle, &p_dist, &p_count ) ); 
}

static	int SG_cyclic_4( p_omat, p_imat, p_noid, p_axestype, p_center, p_axes )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes;
{
	double	p_angle;
	double	p_dist;
	double	p_count;

	p_angle = 90.0;
	p_dist = 0.0;
	p_count = 4;
	return( SG_helix( p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes,
		&p_angle, &p_dist, &p_count ) ); 
}

static	int SG_cyclic_5( p_omat, p_imat, p_noid, p_axestype, p_center, p_axes )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes;
{
	double	p_angle;
	double	p_dist;
	double	p_count;

	p_angle = 72.0;
	p_dist = 0.0;
	p_count = 5;
	return( SG_helix( p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes,
		&p_angle, &p_dist, &p_count ) ); 
}

static	int SG_cyclic_6( p_omat, p_imat, p_noid, p_axestype, p_center, p_axes )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes;
{
	double	p_angle;
	double	p_dist;
	double	p_count;

	p_angle = 60.0;
	p_dist = 0.0;
	p_count = 6;
	return( SG_helix( p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes,
		&p_angle, &p_dist, &p_count ) ); 
}

static	int SG_cyclic_N( p_omat, p_imat,
	p_noid, p_axestype, p_center, p_axes, p_angle, p_count )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes;
double	*p_angle;
double	*p_count;
{
	double	p_dist;

	p_dist = 0.0;
	return( SG_helix( p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes,
		p_angle, &p_dist, p_count ) ); 
}

static	int SG_cube( p_omat, p_imat, p_noid, p_axestype, p_center, p_axes1, p_axes2 )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes1;
double	*p_axes2;
{

	return( SG_polyhedron( "cube", p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes1, p_axes2, NULL ) ); 
}

static	int SG_dihedral( p_omat, p_imat,
	p_noid, p_axestype, p_center, p_axes1, p_axes2, p_count )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes1;
double	*p_axes2;
double	*p_count;
{

	return( SG_polyhedron( "dihedral", p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes1, p_axes2, p_count ) ); 
}

static	int SG_dodeca( p_omat, p_imat, p_noid, p_axestype, p_center, p_axes1, p_axes2 )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes1;
double	*p_axes2;
{

	return( SG_polyhedron( "dodeca", p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes1, p_axes2, NULL ) ); 
}

static	int SG_ico( p_omat, p_imat, p_noid, p_axestype, p_center, p_axes1, p_axes2 )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes1;
double	*p_axes2;
{
	return( SG_polyhedron( "ico", p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes1, p_axes2, NULL ) ); 
}

static	int SG_octa( p_omat, p_imat, p_noid, p_axestype, p_center, p_axes1, p_axes2 )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes1;
double	*p_axes2;
{
	return( SG_polyhedron( "octa", p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes1, p_axes2, NULL ) ); 
}

static	int SG_tetra( p_omat, p_imat, p_noid, p_axestype, p_center, p_axes1, p_axes2 )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes1;
double	*p_axes2;
{
	return( SG_polyhedron( "tetra", p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes1, p_axes2, NULL ) ); 
}

static	int SG_polyhedron( stype, p_omat, p_imat,
	p_noid, p_axestype, p_center, p_axes1, p_axes2, p_count )
char	stype[];
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes1;
double	*p_axes2;
double	*p_count;
{
	int	noid;
	POINT_T	pts[ 4 ];
	int	cnt;
	int	nmats;
	static	int	smats = 0;
	static	MATRIX_T	*mats = NULL;
	static	MAT_T	p_mat;

	noid = !strcmp( p_noid, "true" );

	pts[0][0] = p_center[0];
	pts[0][1] = p_center[1];
	pts[0][2] = p_center[2];

	pts[1][0] = p_axes1[0];
	pts[1][1] = p_axes1[1];
	pts[1][2] = p_axes1[2];

	pts[2][0] = p_axes2[0];
	pts[2][1] = p_axes2[1];
	pts[2][2] = p_axes2[2];

	if( !strcmp( p_axestype, "relative" ) ){
		NAB_ptadd( pts[1], pts[0], pts[1] );
		NAB_ptadd( pts[2], pts[0], pts[2] );
	}

	if( !strcmp( stype, "tetra" ) )
		cnt = 12;
	else if( !strcmp( stype, "cube" ) )
		cnt = 24;
	else if( !strcmp( stype, "octa" ) )
		cnt = 24;
	else if( !strcmp( stype, "dodeca" ) ){
		fprintf( stderr,
			"SG_polyhedron: dodeca not implemented.\n" );
		return( -1 );
	}else if( !strcmp( stype, "ico" ) )
		cnt = 60;
	else if( !strcmp( stype, "dihedral" ) )
		cnt = 2 * *p_count;

	nmats = cnt;
	if( nmats > smats ){
		if( smats > 0 )
			free( mats );
		mats = MAT_ALLOC( nmats );
		if( mats == NULL ){
			fprintf( stderr,
				"SG_polyhedron: allocation of mats failed.\n"
				);
			return( -1 );
		}else
			smats = nmats;
	}

	if( !strcmp( stype, "tetra" ) )
		MAT_tetra( pts, mats );
	else if( !strcmp( stype, "cube" ) )
		MAT_cube( pts, mats );
	else if( !strcmp( stype, "octa" ) )
		MAT_octa( pts, mats );
	else if( !strcmp( stype, "dodeca" ) ){
		fprintf( stderr,
			"SG_polyhedron: dodeca not implemented.\n" );
		return( -1 );
	}else if( !strcmp( stype, "ico" ) )
		MAT_ico( pts, mats );
	else if( !strcmp( stype, "dihedral" ) )
		MAT_dihedral( pts, cnt, mats );

	p_mat.nmats = nmats;
	p_mat.mats = mats;
	return( mk_omats( p_omat, p_imat, &p_mat, noid ) );
}

static	int SG_rotate( p_omat, p_imat, p_axestype, p_center, p_axes, p_angle )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_axestype;
double	*p_center;
double	*p_axes;
double	*p_angle;
{
	char	*p_noid;
	double	p_dist;
	double	p_count;

	p_noid = "true";
	p_dist = 0.0;
	p_count = 2;
	return( SG_helix( p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes,
		p_angle, &p_dist, &p_count ) );
}

static	int SG_translate( p_omat, p_imat, p_axestype, p_center, p_axes, p_dist )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_axestype;
double	*p_center;
double	*p_axes;
double	*p_dist;
{
	char	*p_noid;
	double	p_angle;
	double	p_count;

	p_noid = "true";
	p_angle = 0.0;
	p_count = 2;
	return( SG_helix( p_omat, p_imat,
		p_noid, p_axestype, p_center, p_axes,
		&p_angle, p_dist, &p_count ) );
}

static	int SG_helix( p_omat, p_imat,
	p_noid, p_axestype, p_center, p_axes, p_angle, p_dist, p_count )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_noid;
char	*p_axestype;
double	*p_center;
double	*p_axes;
double	*p_angle;
double	*p_dist;
double	*p_count;
{
	int	noid;
	POINT_T	pts[ 4 ];
	double	ang, dst;
	int	cnt;
	int	nmats;
	static	int	smats = 0;
	static	MATRIX_T	*mats = NULL;
	MAT_T	p_mat;

	noid = !strcmp( p_noid, "true" );

	pts[0][0] = p_center[0];
	pts[0][1] = p_center[1];
	pts[0][2] = p_center[2];

	pts[1][0] = p_axes[0];
	pts[1][1] = p_axes[1];
	pts[1][2] = p_axes[2];

	if( !strcmp( p_axestype, "relative" ) )
		NAB_ptadd( pts[1], pts[0], pts[1] );

	ang = *p_angle;
	dst = *p_dist;
	cnt = *p_count;

	nmats = cnt;
	if( nmats > smats ){
		if( smats > 0 )
			free( mats );
		mats = MAT_ALLOC( nmats );
		if( mats == NULL ){
			fprintf( stderr,
				"SG_helix: allocation of mats failed.\n" );
			return( -1 );
		}else
			smats = nmats;
	}

	MAT_helix( pts, &ang, &dst, &cnt, mats );

	p_mat.nmats = nmats;
	p_mat.mats = mats;
	return( mk_omats( p_omat, p_imat, &p_mat, noid ) );
}

static	int SG_orient( p_omat, p_imat, p_axestype, p_center,
	p_axes1, p_axes2, p_axes3, p_angle1, p_angle2, p_angle3 )
MAT_T	*p_omat;
MAT_T	*p_imat;
char	*p_axestype;
double	*p_center;
double	*p_axes1;
double	*p_axes2;
double	*p_axes3;
double	*p_angle1;
double	*p_angle2;
double	*p_angle3;
{
	double	center[3];
	double	pts[4][3];
	double	angs[3];
	MATRIX_T	mats[1];
	MAT_T	p_mat;

	pts[0][0] = p_center[0];
	pts[0][1] = p_center[1];
	pts[0][2] = p_center[2];

	pts[1][0] = p_axes1[0];
	pts[1][1] = p_axes1[1];
	pts[1][2] = p_axes1[2];

	pts[2][0] = p_axes2[0];
	pts[2][1] = p_axes2[1];
	pts[2][2] = p_axes2[2];

	pts[3][0] = p_axes3[0];
	pts[3][1] = p_axes3[1];
	pts[3][2] = p_axes3[2];

	if( !strcmp( p_axestype, "relative" ) ){
		NAB_ptadd( pts[1], pts[0], pts[1] );
		NAB_ptadd( pts[2], pts[0], pts[2] );
		NAB_ptadd( pts[3], pts[0], pts[3] );
	}

	angs[0] = *p_angle1;
	angs[1] = *p_angle2;
	angs[2] = *p_angle3;

	MAT_orient( pts, angs, mats );

	p_mat.nmats = 1;
	p_mat.mats = mats;
	return( mk_omats( p_omat, p_imat, &p_mat, 0 ) );
}

static	int SG_merge( p_omat, p_imat1, p_imat2, p_imat3, p_imat4 )
MAT_T	*p_omat;
MAT_T	*p_imat1;
MAT_T	*p_imat2;
MAT_T	*p_imat3;
MAT_T	*p_imat4;
{
	int	nomats, nimats1, nimats2, nimats3, nimats4;
	int	m, mo;
	MATRIX_T	*omats;

	nimats1 = p_imat1 ? p_imat1->nmats : 0;
	nimats2 = p_imat2 ? p_imat2->nmats : 0;
	nimats3 = p_imat3 ? p_imat3->nmats : 0;
	nimats4 = p_imat4 ? p_imat4->nmats : 0;

	nomats = nimats1 + nimats2 + nimats3 + nimats4;

	if( p_omat->nmats != nomats ){
		if( p_omat->nmats > 0 )
			free( p_omat->mats );
		omats = MAT_ALLOC( nomats );
		if( omats == NULL ){
			fprintf( stderr,
				"SG_merge: reallocation of omats failed.\n" );
			return( -1 );
		}
		p_omat->nmats = nomats;
		p_omat->mats = omats;
	}else
		omats = p_omat->mats;

	for( mo = 0, m = 0; m < nimats1; m++, mo++ )
		NAB_matcpy( omats[  mo ], p_imat1->mats[ m ] );
	for( m = 0; m < nimats2; m++, mo++ )
		NAB_matcpy( omats[  mo ], p_imat2->mats[ m ] );
	for( m = 0; m < nimats3; m++, mo++ )
		NAB_matcpy( omats[  mo ], p_imat3->mats[ m ] );
	for( m = 0; m < nimats4; m++, mo++ )
		NAB_matcpy( omats[  mo ], p_imat4->mats[ m ] );
		
	return( 0 );
}

static	int SG_multiply( p_omat, p_imat1, p_imat2 )
MAT_T	*p_omat;
MAT_T	*p_imat1;
MAT_T	*p_imat2;
{
	int	nomats, nimats1, nimats2;
	MATRIX_T	*imats1, *imats2, *omats;
	int	m, mo;
	
	nimats1 = p_imat1->nmats;
	nimats2 = p_imat2->nmats;

	nomats = nimats1 > nimats2 ? nimats1 : nimats2;

	if( p_omat->nmats != nomats ){
		if( p_omat->nmats > 0 )
			free( p_omat->mats );
		omats = MAT_ALLOC( nomats );
		if( omats == NULL ){
			fprintf( stderr,
				"SG_multiply: reallocation of omats failed.\n"
				);
			return( -1 );
		}
		p_omat->nmats = nomats;
		p_omat->mats = omats;
	}else
		omats = p_omat->mats;

	if( nimats1 > nimats2 ){
		for( mo = 0, m = 0; m < nimats2; m++, mo++ ){
			NAB_matcpy( omats[mo], 
				MAT_concat( p_imat1->mats[m],
				p_imat2->mats[m] ));
		}
		for( m = nimats2; m < nimats1; m++, mo++ ){
			NAB_matcpy( omats[mo], 
				MAT_concat( p_imat1->mats[m],
				p_imat2->mats[nimats2-1] ) );
		}
	}else{
		for( mo = 0, m = 0; m < nimats1; m++, mo++ ){
			NAB_matcpy( omats[mo], 
				MAT_concat( p_imat1->mats[m],
				p_imat2->mats[m] ) );
		}
		for( m = nimats1; m < nimats2; m++, mo++ ){
			NAB_matcpy( omats[mo], 
				MAT_concat( p_imat1->mats[nimats1-1],
				p_imat2->mats[m] ) );
		}
	}

	return( 0 );
}

static	int SG_read( p_omat, p_fname )
MAT_T	*p_omat;
double	*p_fname;
{

	return( 0 );
}

static	int	mk_omats( p_omat, p_imat, p_mat, noid )
MAT_T	*p_omat;
MAT_T	*p_imat;
MAT_T	*p_mat;
int	noid;
{
	int	nimats, nmats, nomats;
	int	m, ms, mi, mo;
	MATRIX_T	*imats, *mats, *omats;

	if( p_imat ){
		nimats = p_imat->nmats;
		imats = p_imat->mats;
	}else
		nimats = 0;

	nmats = p_mat->nmats;
	mats = p_mat->mats;

	if( nimats == 0 ){
		if( noid ){
			ms = 1;
			nomats = nmats - 1;
		}else{
			ms = 0;
			nomats = nmats;
		}
	}else{
		if( noid ){
			ms = 1;
			nomats = nimats * ( nmats - 1 );
		}else{
			ms = 0;
			nomats = nimats * nmats;
		}
	}

	if( p_omat->nmats != nomats ){
		if( p_omat->nmats > 0 )
			free( p_omat->mats );
		omats = MAT_ALLOC( nomats );
		if( omats == NULL ){
			fprintf( stderr,
				"mk_omats: reallocation of omats failed.\n" );
			return( -1 );
		}
		p_omat->nmats = nomats;
		p_omat->mats = omats;
	}else
		omats = p_omat->mats;

	if( nimats == 0 ){
		mo = 0;
		for( m = ms; m < nmats; m++, mo++ )
			NAB_matcpy( omats[ mo ], mats[ m ] );
	}else{
		mo = 0;
		for( m = ms; m < nmats; m++ ){
			for( mi = 0; mi < nimats; mi++, mo++ ){
				NAB_matcpy( omats[  mo ],
					MAT_concat( imats[ mi ], mats[ m ] ) );
			}
		}
	}

	return( 0 );
}
