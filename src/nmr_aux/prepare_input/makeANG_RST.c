#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int make_RST_restraint( int tnum, int rnum, char *rname, char *angname,
                        float lb, float ub );
int read_lib();
int residue_offset( char *a, char *b );

#define USAGE	"usage: %s -help\n       %s -pdb ambpdb_file [-con constraint] [-lib libfile] [-les lesfile] \n", argv[ 0 ], argv[ 0 ]
#define ERROR	(1)
#define SUCCESS	0

#define RAD(a)	((a)*.01745329251994)
#define DEG(a)	((a)/.01745329251994)

#define	DEFAULT_LIB "dat/tordef.lib"

#define REPORT	0
#define WORDSIZE	20

#define RNAME_SIZE	5
#define ANAME_SIZE	5
#define ANGNAME_SIZE	10
#define MAX_ATOMS	50000
#define MAX_TDEF	500
#define MAX_CONS	5000

#define PUCKER_AMP	38.5
#define R1_R4_RANGE	1.0
#define MIN_NU_RANGE	30.0
#define K_FORCE_LB	2.0
#define K_FORCE_UB	2.0

#define MAXFILENAME 128

int	torcnt=0;

typedef struct tordef {
	char    type[WORDSIZE+WORDSIZE];
	char    A[5][WORDSIZE];
	int     R[5], pseudo;
	float   J[3];
} TDEF;

TDEF	torlib[MAX_TDEF];
int	ndef;
int	pseudocnt=0;

typedef struct	atom	{
	int	rnum, anum;
	char	rname[RNAME_SIZE], aname[ANAME_SIZE];
}ATOM;

ATOM	atom_table[MAX_ATOMS];
int	nat;
char	line[100];

typedef struct restraint {
	int	rnum;
	ATOM	*ap[4];
	TDEF	*tdef;
	float	lb, ub;
}RST;

int   read_les = 0;
RST	rst[MAX_CONS];
FILE	*pdb=NULL, *cons=NULL, *lib=NULL, *les=NULL;

#define MAXCOPY 10
#define MAXATOMS 5000
int		invcnt[ MAXATOMS ];
int		invcop[ MAXATOMS ][ MAXCOPY ];

static char *help[] = {
"Usage: makeANG_RST -help\n",
"       makeANG_RST -pdb ambpdb_file [-con constraint] [-lib libfile]\n",
"                [-les lesfile ]\n",
"\n",
"This program takes as input a five column torsion angle constraint file\n",
"along with an AMBER pdb file of the molecule.\n",
"\n",
"It creates as output (to standard out) a list of constraints in RST format\n",
"that is readable by AMBER.\n",
"\n",
"--------------------------------------------------------------------------------\n",
"\n",
"The input torsion angle constraint file can be read from standard in or from\n",
"a file specified by the -con option on the command line.\n",
"\n",
"The input constraint file should look something like this:\n",
"1	GUA	PPA	111.5	144.0\n",
"2	CYT	EPSILN	20.9	100.0\n",
"2	CYT	PPA	115.9	134.2 \n",
"3	THY	ALPHA	20.4	35.6\n",
"4	ADE	GAMMA	54.7	78.8\n",
"5	GLY	PHI	30.5	60.3\n",
"6	ALA	CHI	20.0	50.0\n",
"....\n",
"\n",
"The first column is the residue number, the second is the residue name\n",
"(three letter code, or as defined in your own personal torsion library file).\n",
"Third is the angle name. Fourth is the lower bound.  Last is the upper bound.\n",
"\n",
"NOTE:  It is assumed that the lower bound and the upper bound define a\n",
"region of allowed conformation on the unit circle that is swept out in a\n",
"clockwise direction from lb->ub.  If the number in the lb column is > the\n",
"the number in the ub column, 360 will successively be subtracted from the lb\n",
"until lb < ub.  This preserves the clockwise definition of the allowed\n",
"conformation space, while also making the number that specifies the lower\n",
"bound less than the number that specifies the upper bound, as is required by\n",
"AMBER.  If this occurs, a warning message will be printed to stderr to\n",
"notify the user that the data has been modified.\n",
"\n",
"--------------------------------------------------------------------------------\n",
"\n",
"The angles that one can constrain in this manner are defined in the library\n",
"file that can be optionally specified on the command line with the -lib\n",
"flag, or the default library tordef.lib (written by Garry P. Gippert) will\n",
"be used.  If you wish to specify your own nomenclature, or add angles that\n",
"are not already defined in the default file, you should make a copy of this\n",
"file and modify it to suit your needs.\n",
"\n",
"The general format for an entry in the library is:\n",
"\n",
"LEU   PSI     N    CA   C    N+\n",
"\n",
"where the first column is the residue name, the second column is the angle\n",
"name that will appear in the input file when specifying this angle, and the\n",
"last four columns are the atom names that define the torsion angle.  When a\n",
"torsion angle contains atom(s) from a preceding or succeeding residue in the\n",
"structure, a '-' or '+' is appended to those atom names in the library,\n",
"thereby specifying that this is the case.  In the example above, the atoms\n",
"that define PSI for LEU residues are the N, CA, and C atoms of that same LEU\n",
"and the N atom of the residue after that LEU in the primary structure.  Note\n",
"that the order of atoms in the definition is important and should reflect\n",
"that the torsion angle rotates about the two central atoms as well as the\n",
"fact that the four atoms are bonded in the order that is specified in the\n",
"definition.\n",
"\n",
"If the first letter of the second field is 'J', this torsion is assumed\n",
"to be a J-coupling constraint.  In that case, three additional floats are\n",
"read at the end of the line, giving the A,B and C coefficients for the\n",
"Karplus relation for this torsion.  For example:\n",
"\n",
"ALA JHNA    H    N    CA   HA   9.5  -1.4  0.3\n",
"\n",
"will set up a J-coupling restraint for the HN-HA 3-bond coupling, assuming\n",
"a Karplus relation with A,B, C as 9.5, -1.4 and 0.3.  (These particular \n",
"values are from Bruschweiler and Case, JACS 116: 11199 (1994).)\n",
"\n",
"--------------------------------------------------------------------------------\n",
"\n",
"This program also supports pseudorotation phase angle constraints, and for\n",
"each of these will generate restraints for the 5 component angles which\n",
"correspond to the lb and ub values of the input psuedorotation constraint.\n",
"In the torsion library, a pseudorotation definition looks like:\n",
"\n",
"PSEUDO        CYT     PPA     NU0     NU1     NU2     NU3     NU4\n",
"CYT   NU0     C4'     O4'     C1'     C2'\n",
"CYT   NU1     O4'     C1'     C2'     C3'\n",
"CYT   NU2     C1'     C2'     C3'     C4'\n",
"CYT   NU3     C2'     C3'     C4'     O4'\n",
"CYT   NU4     C3'     C4'     O4'     C1'\n",
"\n",
"The first line describes that a PSEUDOrotation angle is to be defined for\n",
"CYT that is called PPA and is made up of the four angles NU0-NU4.  Then the\n",
"definition for NU0-NU4 should also apper in the file in the same format as\n",
"the example given above for LEU PSI.\n",
"\n",
"PPA stands for Pseudorotation Phase Angle and is the angle that should\n",
"appear in the input constraint file when using pseudorotation constraints.\n",
"The program then uses the definition of that PPA angle in the library file\n",
"to look for the 5 other angles (NU0-NU4 in this case) which it then\n",
"generates restraints for.  PPA for proline residues is included in the\n",
"standard library as well as for the DNA nucleotides.\n",
"\n",
"--------------------------------------------------------------------------------\n",
"\n",
"If the -les flag is set, the program will prepare torsion angle restraints\n",
"for multiple copies (LES) simulations.  In this case, the input pdb file\n",
"is one *without* LES copies, i.e. with just a single copy of the molecule.\n",
"The lesfile specified by this flag is created by the addles program, and\n",
"contains a mapping from original atom numbers into the copy numbers used\n",
"in the multiple-copies simulation.\n",
"\n",
"Torsion angle constraints defined here *cannot* span two different copy\n",
"sets, i.e., there cannot be some atoms of a particular torsion that are in\n",
"one multiple copy set, and other atoms from the same torsion that are in\n",
"other copy sets.  It *is* OK to have some atoms with single copies, and\n",
"others with multiple copies in the same torsion.  The program will create\n",
"as many duplicate torsions as there are copies.\n",
"\n",
"******************\n",
"makeANG_RST\n",
"jsmith@Scripps.EDU\n",
0
};

void print_RST()
{
	int	i, j, anum[4], bnum[4];
	int	ij, cntmax, icnt[4];
	float	r1, r2, r3, r4, rk2p;

	for( i=0; i<torcnt; i++ ){
		if( read_les ){           /*  loop over LES copies   */

			anum[0] = rst[i].ap[0]->anum; anum[1] = rst[i].ap[1]->anum;
			anum[2] = rst[i].ap[2]->anum; anum[3] = rst[i].ap[3]->anum;

			icnt[0] = invcnt[ anum[0] ]; cntmax = icnt[0];
			icnt[1] = invcnt[ anum[1] ]; if( icnt[1]> cntmax ) cntmax=icnt[1];
			icnt[2] = invcnt[ anum[2] ]; if( icnt[2]> cntmax ) cntmax=icnt[2];
			icnt[3] = invcnt[ anum[3] ]; if( icnt[3]> cntmax ) cntmax=icnt[3];

			r1 = rst[i].lb - R1_R4_RANGE; r4 = rst[i].ub + R1_R4_RANGE;
			r2 = rst[i].lb; r3 = rst[i].ub;

			for( ij=0; ij<cntmax; ij++ ){

				bnum[0] = icnt[0]==1 ? invcop[anum[0]][0] : invcop[anum[0]][ij];
				bnum[1] = icnt[1]==1 ? invcop[anum[1]][0] : invcop[anum[1]][ij];
				bnum[2] = icnt[2]==1 ? invcop[anum[2]][0] : invcop[anum[2]][ij];
				bnum[3] = icnt[3]==1 ? invcop[anum[3]][0] : invcop[anum[3]][ij];

				printf("# %d %s:  ", rst[i].rnum, rst[i].tdef->type);
				for( j=0; j<3; j++ ){
					printf("(%d %s %s)-", rst[i].ap[j]->rnum, rst[i].ap[j]->rname,
						rst[i].ap[j]->aname);
				}
				printf("(%d %s %s) %5.1f %5.1f\n",
					rst[i].ap[3]->rnum, rst[i].ap[3]->rname,
					rst[i].ap[3]->aname, rst[i].lb, rst[i].ub);
				printf(" &rst     iat = %5d, %5d, %5d, %5d,\n",
					bnum[0], bnum[1], bnum[2], bnum[3]);
				printf("	  r1 = %5.1f, r2 = %5.1f, r3 = %5.1f, r4 = %5.1f,\n", 
					r1, r2, r3, r4);
				if( rst[i].tdef->J[0] != 0.0 )
					printf("	  rjcoef = %.2f, %.2f, %.2f,\n", 
					rst[i].tdef->J[0], rst[i].tdef->J[1], rst[i].tdef->J[2] );
				rk2p = ij<cntmax-1 ? -K_FORCE_LB : K_FORCE_LB;
				printf("	  rk2 = %5.1f, ", rk2p);
				printf("rk3 = %5.1f,				/\n\n", K_FORCE_UB);
			}

		} else {                  /*  no LES copies to consider  */

			anum[0] = rst[i].ap[0]->anum; anum[1] = rst[i].ap[1]->anum;
			anum[2] = rst[i].ap[2]->anum; anum[3] = rst[i].ap[3]->anum;
			r1 = rst[i].lb - R1_R4_RANGE; r4 = rst[i].ub + R1_R4_RANGE;
			r2 = rst[i].lb; r3 = rst[i].ub;
			printf("# %d %s:  ", rst[i].rnum, rst[i].tdef->type);
			for( j=0; j<3; j++ ){
				printf("(%d %s %s)-", rst[i].ap[j]->rnum, rst[i].ap[j]->rname,
					rst[i].ap[j]->aname);
			}
			printf("(%d %s %s) %5.1f %5.1f\n",
				rst[i].ap[3]->rnum, rst[i].ap[3]->rname,
				rst[i].ap[3]->aname, rst[i].lb, rst[i].ub);
			printf(" &rst     iat = %5d, %5d, %5d, %5d,\n",
				anum[0], anum[1], anum[2], anum[3]);
			printf("	  r1 = %5.1f, r2 = %5.1f, r3 = %5.1f, r4 = %5.1f,\n", 
				r1, r2, r3, r4);
			if( rst[i].tdef->J[0] != 0.0 )
				printf("	  rjcoef = %.2f, %.2f, %.2f,\n", 
				rst[i].tdef->J[0], rst[i].tdef->J[1], rst[i].tdef->J[2] );
			if( i == 0 ){
				printf("	  rk2 = %5.1f, ", K_FORCE_LB);
				printf("rk3 = %5.1f,				/\n\n", K_FORCE_UB);
			}else{
				printf("	/\n\n");
			}
		}
	}
}

void read_cons()
{
	char	rname[RNAME_SIZE], angname[ANGNAME_SIZE];
	int	rnum, tors;
	float	lb, ub;

	while ( fgets(line, sizeof(line), cons) ){
		if ( strncmp( line, "#", 1 ) == 0 ) continue;
		if ( strlen( line ) < 5 ) continue; /* for blank-like lines */
		tors = 0;
		sscanf (line, "%d %s %s %f %f", &rnum, rname, angname, &lb, &ub);

		/* only use first three characters of rname, for compatibility with
			DIANA constraint format:   */
		rname[3] = '\0';

		if ((tors = make_RST_restraint( torcnt, rnum, rname, angname, lb, ub )))
			torcnt += tors;
	}
	fprintf(stderr, "\n%d constraints (%d torsion + %d pseudorotation) converted to RST format.\n\n",
		torcnt, torcnt-pseudocnt*5, pseudocnt);
}

void make_atom_table( FILE *pdb)
{
	for ( nat=0; (fgets(line, sizeof(line), pdb)) && (nat < MAX_ATOMS); nat++ ){
		if ( (strncmp( "ATOM", line, 4 ) == 0) || (strncmp( "HETATM", line, 6 ) == 0) ){
			sscanf(&line[ 6 ], "%d %s %s %d", &atom_table[ nat ].anum, atom_table[ nat ].aname,
							  atom_table[ nat ].rname, &atom_table[ nat ].rnum);
		}else nat--;
	}
	fclose( pdb );
}

void cmd_parse( int argc, char *argv[] )
{
	int	ac, i;
	char	pf[100], lf[100], cf[100], lesf[100];
	char	f_lib[MAXFILENAME];
	char	*ambh;

	if ( argc < 3 ){ 
		if ( (argc == 2) && (strcmp( argv[ 1 ], "-help") == 0) ){
			for ( i = 0; help[i] != 0; i++ ) fprintf( stderr, "%s", help[i] );
			exit( SUCCESS );
		}else{
			fprintf( stderr, USAGE );
			exit(ERROR);
		}
	}
	for ( ac=1; ac < argc; ac++){
		if ( strcmp( argv[ ac ], "-pdb" ) == 0 ){
			strcpy( pf, argv[ ac +1 ]);
			if ( ac == argc - 1 ){
				fprintf( stderr, USAGE );
				exit(ERROR);
			}else if ( (pdb = fopen( argv[ ac + 1 ], "r")) == NULL ) {
				fprintf( stderr, "Can't open pdb file %s\n", argv[ ac + 1 ] );
				exit (ERROR);
			}
		}else if ( (strcmp( argv[ ac ], "-con" )) == 0 ){
			strcpy( cf, argv[ ac + 1 ]);
			if ( ac == argc - 1 ){
				fprintf( stderr, USAGE );
				exit(ERROR);
			}else if ( (cons = fopen( argv[ ac + 1 ], "r")) == NULL ) {
				fprintf(stderr, "Can't open constraint file %s\n", argv[ ac + 1 ] );
				exit (ERROR);
			}
		}else if ( (strcmp( argv[ ac ], "-les" )) == 0 ){
			read_les = 1;
			strcpy( lesf, argv[ ac + 1 ]);
			if ( ac == argc - 1 ){
				fprintf( stderr, USAGE );
				exit(ERROR);
			}else if ( (les = fopen( argv[ ac + 1 ], "r")) == NULL ) {
				fprintf(stderr, "Can't open constraint file %s\n", argv[ ac + 1 ] );
				exit (ERROR);
			}
		}else if ( (strcmp( argv[ ac ], "-lib" )) == 0 ){
			strcpy( lf, argv[ ac +1 ]);
			if ( ac == argc - 1 ){
				fprintf( stderr, USAGE );
				exit(ERROR);
			}else if ( (lib = fopen( argv[ ac + 1 ], "r")) == NULL ) {
				fprintf(stderr, "Can't open torsion library file %s\n", argv[ ac + 1 ] );
				exit (ERROR);
			}
		}
	}
	if (cons == NULL){
		cons = stdin;
		fprintf(stderr, "Reading constraints from stdin.\n");
	}
	if (lib == NULL){
		if((ambh = getenv( "AMBERHOME" ))){
			strncpy( f_lib, ambh, MAXFILENAME-40 );
			strcat( f_lib, "/" );
		}else{
			strcpy( f_lib, "./" );
		}
		strcat(f_lib,DEFAULT_LIB);
		if ( (lib = fopen( f_lib, "r")) == 0 ){
			fprintf(stderr, "Can't open default torsion library file %s\n",
				f_lib );
			exit (ERROR);
		}else{
			ndef = read_lib();
			fprintf(stderr, "Read %d torsion definitions from default library file %s\n",
				ndef, f_lib );
		}
	}else{
		ndef = read_lib();
		fprintf(stderr, "Read %d torsion definitions from library file %s\n",
			ndef, lf );
	}
	if( pdb == NULL){
		fprintf( stderr, USAGE );
		exit( ERROR );
	}else{
		make_atom_table(pdb);
		if( nat > 0 ){
			printf( "# %d atoms read from pdb file %s.\n", nat, pf );
			fprintf(stderr, "Read %d atoms from %s for RST atom numbers:\n",
				nat, pf);
			fprintf(stderr, "Be sure that the atom numbers are the same ");
			fprintf(stderr, "as in the prmtop file that you will use.\n");
		}else{
			fprintf(stderr, "ERROR: Could not get any atoms from pdb file %s\n", pf);
			exit(ERROR);
		}
	}
	fprintf(stderr, "Reading constraints from %s\n", cf);
}



/*	NOTE:  These next two routines are borrowed from Garry Gippert's torvio program
	to ensure compatibility between the treatment of the torsion library in torvio
	and torconv.
*/
int read_lib()
{
	TDEF	*tdef;
	int	lineno, nd, i;
	float   J[3];
	char    rtyp[WORDSIZE], ttyp[WORDSIZE], word[WORDSIZE], line[100];
	char    A[5][WORDSIZE];
 
	lineno = 0;
	tdef = torlib;
	nd = 0; 
 
	while ( fgets( line, sizeof(line), lib ) ){
		lineno++;
		sscanf( line, "%s", word );
 		if ( strncmp( word, "#", 1) == 0 ) continue;
		if ( strcmp( word, "PSEUDO" ) == 0 ){
			if( 7 != sscanf( line, "%*s %s %s %s %s %s %s %s",
				rtyp, ttyp, A[0], A[1], A[2], A[3], A[4] ) ) fprintf( stderr, 
				"Trouble reading line %d from torsion library.\n", lineno);
			
			sprintf( tdef->type, "%s %s", rtyp, ttyp );
			for ( i = 0; i < 5; i++ )
			{
				strcpy( tdef->A[i], A[i] );
			}
			tdef->pseudo = 1;
		}else{
			if( 6 != sscanf( line, "%s %s %s %s %s %s",
				rtyp, ttyp, A[0], A[1], A[2], A[3] ) ) fprintf( stderr, 
				"Trouble reading line %d from torsion library.\n", lineno);
 
			sprintf( tdef->type, "%s %s", rtyp, ttyp );
			for ( i = 0; i < 4; i++ )
			{
				tdef->R[i] = residue_offset( A[i], tdef->A[i] );
				if ( REPORT )
				printf( "%s %s -> %s %d\n",
					tdef->type, A[i], tdef->A[i], tdef->R[i] );
			}
/*                           re-scan the line to pick up torsion constants: */
			if( ttyp[0] == 'J' ){
				if( 9 != sscanf( line, "%s %s %s %s %s %s %f %f %f",
					rtyp, ttyp, A[0], A[1], A[2], A[3], &J[0], &J[1], &J[2] ) )
					fprintf( stderr, 
					"Trouble reading line %d from torsion library.\n", lineno);
				tdef->J[0] = J[0];
				tdef->J[1] = J[1];
				tdef->J[2] = J[2];
			} else {
				tdef->J[0] = 0.0;
				tdef->J[1] = 0.0;
				tdef->J[2] = 0.0;
			}
		tdef->pseudo = 0;
		}
		nd++;
		tdef++;
	}
	fclose( lib );
	return(nd);
}

int residue_offset( char *a, char *b )
{
	int     l = strlen(a);
	char    s = a[l-1];
 
	strcpy( b, a );
	switch( s ) {
		case '-':
			b[l-1]='\0';
			return ( -1 );
		case '+':
			b[l-1]='\0';
			return ( 1 );
		default:
			return ( 0 );
	}
}

float *get_bounds( float *bvec, int j, float P_lower, float P_upper)
{
	float	nu;
	float	P, nu_min, nu_max;
	float	range, d;

	nu_min = ( RAD(PUCKER_AMP)*( cos(RAD(P_lower+(144*(j-2)))) ) );
	nu_max = ( RAD(PUCKER_AMP)*( cos(RAD(P_upper+(144*(j-2)))) ) );
	for( P = P_lower; P <= P_upper; P += 0.1 ){
		nu = ( RAD(PUCKER_AMP)*( cos(RAD(P+(144*(j-2)))) ) );
		if ( nu < nu_min ) nu_min = nu;
		if ( nu > nu_max ) nu_max = nu;
	}
	nu_min = DEG(nu_min);
	nu_max = DEG(nu_max);
	if ( (range = nu_max - nu_min) < MIN_NU_RANGE ){
		d = (MIN_NU_RANGE - range)/2;
		nu_min -= d;
		nu_max += d;
	}
	bvec[ 0 ] = nu_min;
	bvec[ 1 ] = nu_max;
	return( SUCCESS );
}

int make_RST_restraint( int tnum, int rnum, char *rname, char *angname,
                        float lb, float ub )
{
	int	i, acnt, dcnt=1, ares, pcnt=0;
	float	bvec[2];
	char	type[WORDSIZE];
	TDEF	*tdef;
	ATOM	*a;

	sprintf( type, "%s %s", rname, angname );
	tdef = torlib;

	while( (strcmp( type, tdef->type)) && (dcnt <= ndef) ){
		tdef++;
		dcnt++;
	}
	if( strcmp( type, tdef->type ) ){
		fprintf( stderr, "ERROR: %d %s: No definition for %s in torsion library.  Skipping...\n",
			rnum, type, type);
		return(0);
	}
	rst[tnum].tdef = tdef;
	rst[tnum].rnum = rnum;
	if ( lb > ub ){
		fprintf( stderr, "WARNING: %d %s %5.1f %5.1f:  lb > ub -->\n",
			rnum, tdef->type, lb, ub);
	}
	while ( lb > ub ){
		lb-=360;
		if ( lb < ub ){
			fprintf( stderr, "   lb reset to %5.1f for clockwise restraint definition.\n", lb);
		}
	}
	rst[tnum].lb = lb;
	rst[tnum].ub = ub;
	if( tdef->pseudo == 0 ){
		for( i = 0; i < 4; i++ ){
			ares = rnum + tdef->R[i];
			a = atom_table;
			acnt = 1;
			while( ( ares != a->rnum) && (acnt < nat) ){
				a++;
				acnt++;
			}
			if( ares != a->rnum ){
				fprintf( stderr, "ERROR: %d %s: Couldn't find residue %d in pdb file. Skipping...\n",
					rnum, type, ares );
				return(0);
			}
			while ( strcmp( tdef->A[i], a->aname ) ){ 
				if( ares != a->rnum ){
					fprintf( stderr, "%d %s: Couldn't find atom", rnum, type);
					fprintf( stderr, "%s in residue %d of pdb file.\n",
						tdef->A[i], ares );
					return(0);
				}
				a++;
			}
			rst[tnum].ap[i] = a;
		}
		return(1);
	}else{
		for( i = 0; i < 5; i++ ){
			get_bounds( bvec, i, lb, ub );
			pcnt += make_RST_restraint( tnum+i, rnum, rname, tdef->A[i], bvec[0], bvec[1] );
		}
		pseudocnt++;
		return( pcnt );
	}
}

int main( int argc, char *argv[] )
{
	int i, j, natles;

	fprintf( stderr, "\n******************************************************\n");
	fprintf( stderr, "makeANG_RST v 2.0, jsmith@scripps.edu\n"); 
	cmd_parse( argc, argv );
	fprintf( stderr, "******************************************************\n\n");
	if( read_les ){
		fscanf( les, "%d", &natles );
		if( natles != nat ){
			fprintf( stderr, "Bad number of atoms in les file: %d\n", natles );
			exit( 1 );
		}
		for( i=1; i<=nat; i++ ){
			fscanf( les, "%d", &invcnt[ i ] );
			for( j=0; j<invcnt[i]; j++ ){
				fscanf( les, "%d", &invcop[ i ][ j ] );
			}
		}
	}

	read_cons();
	print_RST();
	exit(0);
}

