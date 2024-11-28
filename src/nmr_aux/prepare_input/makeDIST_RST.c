/*
** makeDIST_RST converts a "simple" NMR constraint format to
**    the &rst namelist used by sander.
**
**    input files are very similar to the constraint formats used by
**      programs like DIANA and DISGEO;  output is designed primarily
**      for the SANDER module of AMBER.  Some options are used for other
**      purposes.
**
**    written by Garry P. Gippert, TSRI.
**    LES and AMBIGUOUS constraint modifications by D.A. Case.
*/

#include <ctype.h>
#include <stdlib.h>
#include "global.h"

static int lookupanum( int rnum, char *atyp );
static int lookuprtyp( int rnum, char *rtyp );
static int mapdefinition( char *atyp, char *rtyp );
static int read_data( FILE *fp, char *data, char *comment );
static int read_maplib( FILE *fp );
static int read_pdb_file( FILE *fp );
static int split( char *str, char wordptr[][WORDSIZE] );

#define DEFAULT_RST	("      rk2=20.0, rk3=20.0, ir6=1,")
#define DEFAULT_MAP ("map.DG-AMBER")

int		idum; 

struct molecule {
	char	name[WORDSIZE+WORDSIZE];
	int	anum;
	char	atyp[WORDSIZE];
	char	rtyp[WORDSIZE];
	int	rnum;
	float	xyz[MAXDIMEN];
	float	radii;
} mol[MAXATOMS];

int 	natoms;

#define MAXMAPLIB 13500

struct mapping {
	char	rtyp[WORDSIZE];            /* residue name */
	char	atyp[WORDSIZE];            /* atom or pseudo-atom name */
	char	atrn[MAXDEGEN][WORDSIZE];  /* expanded (real) atom names */
	int		rtrn[MAXDEGEN];            /* residue number of expanded atom;
		                                  if zero, use the number of the
		                                  resiude being considered   */      
	int	na;                            /* number of expanded atoms  */
	int maxamb;                        /* max. number of expanded atoms that
		                                  could participate; used in deciding
		                                  the upper bound; 1 <= maxamb <= na */
	float	pseudoc;                   /* pseudoatom distance correction */
} map[MAXMAPLIB];

int 	nmaps;

struct pseudoatom {
	int	rnum;
	char	rtyp[WORDSIZE];
	char	nrtyp[WORDSIZE];
	char	atyp[WORDSIZE];
	char	newatyp[MAXDEGEN][WORDSIZE];
	int	newanum[MAXDEGEN];
	int	na;
	char	comment[COMMENTSIZE];
	float   pseudoc;
} atma, atmb;

int		xpk,npk;
int		imix=1;
int		npeak=0;
int		maxres;  /* will hold number of residues in the system */
float	u_bound;
float	l_bound;
float	volume;
float	nu_bound;
float	nl_bound;

#define MAXCOPY 10
int		invcnt[ MAXATOMS ];
int		invcop[ MAXATOMS ][ MAXCOPY ];

static char *banner[] = {
	"# makeDIST_RST\n",
	0
};

static char *help[] = {
	"convert DIANA-like distance restraints to SANDER format\n",
	"\n",
	"input:\n",
	"   -upb <filename>: 7-col_NOE (necessary) distance bound file, OR\n",
	"   -ual <filename>: 8-col_UAL (necessary) upper/lower file, OR\n",
	"   -vol <filename>: 7-col_VOL (necessary) VOLUME file, one or more volumes\n",
	"   -pdb <filename>: PDBFILE using AMBER nomenclature and numbering\n",
	"   -map <filename>: MAP file  (default:map.DG-AMBER)\n",
	"\n",
	"output (all optional)\n",
	"   -dgm <filename>:  DGEOM95 restraint format\n",
	"   -rst <filename>:  SANDER restraint format\n",
	"   -rm6 <filename>:  DISGEO/NOEVIO format\n",
	"   -svf <filename>:  Sander Volume Format\n",
	"   \n",
	"other options:\n",
	"   -nocorr (do not correct upper bound for r**-6 averaging)\n",
	"   -altdis (use alternative form for the distance restraints)\n",
	"   -help   (gives you this explanation, overrides other parameters)\n",
	"   -report (gives you short runtime diagnostic output)\n",
	"   -imix   (for -vol input, says which mixing time; default is 1)\n",
	"\n",
	"errors come to stderr.\n",
	"\n",
	"Dg-nomenclature NOEs are translated to All-atom nomenclature\n",
	"according to mapping functions found in the maplib file.\n",
	0
};

#define MAXFILENAME 128

char	f_upb[MAXFILENAME];
char	f_les[MAXFILENAME];
char	f_ual[MAXFILENAME];
char	f_vol[MAXFILENAME];
char	f_pdb[MAXFILENAME];
char	f_map[MAXFILENAME];
char	f_dgm[MAXFILENAME];
char	f_rst[MAXFILENAME];
char	f_rm6[MAXFILENAME];
char	f_svf[MAXFILENAME];

int 	report=0;
int 	debug=0;
int 	read_upb=0;
int 	read_les=0;
int 	read_ual=0;
int 	read_vol=0;
int 	read_pdb=0;
int 	read_map=0;
int 	print_rst=0;
int 	print_dgm=0;
int 	print_rm6=0;
int 	print_svf=0;
int 	fatal_error=0;
int 	first_rst_output=1;
int		no_corr = 0;
int		altdis = 0;

FILE 	*fp_ual, *fp_les, *fp_upb, *fp_vol, *fp_pdb, *fp_map;
FILE 	*fp_dgm, *fp_rst, *fp_rm6, *fp_svf;


void write_dgm(FILE *fp, char *data, char *comment)
{
	int	i,j;
#ifdef DEBUG
	char 	*module="write_dgm";
	printf("%s: (%s) (%s)\n", module, data, comment);
#endif
	for (i=0;i<atma.na;i++) {
		for (j=0;j<atmb.na;j++) {
/*
**			if ( i==0 && j==0 ) {
**              fprintf(fp,"! %s%s%s%s",
**                data,
**               ( strlen( comment ) == 0 ? "" : " "),
**                comment,
**               ( strlen( comment ) == 0 ? "" : "\n"));
**			}
*/
			if (read_ual==0) {
				nl_bound = -1.0;
			} else {
				nl_bound = l_bound;
			}
			fprintf(fp," 1 %5d %-5s 1 %5d %-5s %8.2f %8.2f\n",
				atma.rnum, atma.newatyp[i],
				atmb.rnum, atmb.newatyp[j],
				nl_bound, nu_bound);
		}
	}
}

void write_svf( FILE *fp,  char *data,  char *comment)
{
#ifdef DEBUG
	char 	*module="write_svf";
#endif

	npeak++;
	fprintf(fp,
		" ihp(%d,%d) = %4d, jhp(%d,%d) = %4d, aexp(%d,%d) = %10.5f,\n",
		imix, npeak, atma.newanum[atma.na-1],
		imix, npeak, atmb.newanum[atmb.na-1],
		imix, npeak, volume);
}

void write_rm6(FILE *fp, char *data, char *comment)
{
	int 	i,j;
#ifdef DEBUG
	char 	*module="write_rm6";
	printf("%s: (%s) (%s)\n", module, data, comment);
#endif
	if ( (strncmp(data,"END OF SET",10)==0) ) {
		fprintf(fp,"%s",data);
		return;
	}
	for (i=0;i<atma.na;i++) {
		for (j=0;j<atmb.na;j++) {
			fprintf(fp,"%5d %-5s %-5s %5d %-5s %-5s %10.4f",
				atma.rnum, atma.nrtyp, atma.newatyp[i],
				atmb.rnum, atmb.nrtyp, atmb.newatyp[j],
				( (i==atma.na-1 && j==atmb.na-1) ? 
						nu_bound : 0.0 ));
			if ( i==0 && j==0 ) {
				fprintf(fp," # %s%s%s%s",
					data,
					( strlen( comment ) == 0 ? "" : " " ),
					comment,
					( strlen( comment ) == 0 ? "" : "\n" ));
			} else {
				fprintf(fp," #\n" );
			}
		}
	}
	return;
}

void write_rst( FILE *fp, char *data, char *comment)
{
	int 	i,j,k;
	int		ntot, itot, ktot[ MAXDEGEN ];
#ifdef DEBUG
	char 	*module="write_rst";
	printf("%s: (%s) (%s)\n",module, data, comment);
#endif

	fprintf(fp,"#\n");
	if ( strlen( comment ) == 0 )
		fprintf(fp,"# %s", data );
	else
		fprintf(fp,"# %s (%s)\n", data, comment);

/*	First check to see that this constraint never refers to itself,
	i.e., that all the atoms "a" are distanct from the atoms "b"    */

	for( i=0; i< atma.na; i++ ){
		for( j=0; j< atmb.na; j++ ){
			if( atma.newanum[i] == atmb.newanum[j] ){
				fprintf(fp,"#This appears to be self-referential: ignoring\n" );
				return;
			}
		}
	}

	fprintf(fp," &rst\n  ixpk= %d, nxpk= %d,", xpk,npk);
	if( read_les ) {

/*    here we translate the "original" atom numbers into LES copies   */
/*      (would need to set up invcnt and invcop somewhere)            */

		if( atma.na == 1 ){
			k = atma.newanum[0];
			atma.na = invcnt[ k ];
			for( i=0; i<atma.na; i++ ){
				atma.newanum[i] = invcop[ k ][ i ];
			}
		} else {
			ntot = 0;
			for( itot=0; itot<atma.na; itot++ ){
				ktot[ itot ] = atma.newanum[ itot ];
			}
			for( itot=0; itot<atma.na; itot++ ){
				k = ktot[itot];
				for( i=0; i<invcnt[k]; i++ ){
					atma.newanum[i+ntot] = invcop[ k ][ i ];
				}
				ntot += invcnt[k];
			}
			atma.na = ntot;
		}
		if( atmb.na == 1 ){
			k = atmb.newanum[0];
			atmb.na = invcnt[ k ];
			for( i=0; i<atmb.na; i++ ){
				atmb.newanum[i] = invcop[ k ][ i ];
			}
		} else {
			ntot = 0;
			for( itot=0; itot<atmb.na; itot++ ){
				ktot[ itot ] = atmb.newanum[ itot ];
			}
			for( itot=0; itot<atmb.na; itot++ ){
				k = ktot[itot];
				for( i=0; i<invcnt[k]; i++ ){
					atmb.newanum[i+ntot] = invcop[ k ][ i ];
				}
				ntot += invcnt[k];
			}
			atmb.na = ntot;
		}
	}

	/* adjust l_bound to be <= nu_bound, so distances are monotonic: */
    if( l_bound > nu_bound ) l_bound = nu_bound;

	fprintf(fp," iat=%4d,%4d, r1=%5.2f, r2=%5.2f, r3=%5.2f, r4=%5.2f,%s",
		(atma.na == 1 ? atma.newanum[0] : -1),
		(atmb.na == 1 ? atmb.newanum[0] : -1),
		(altdis ? 0.0 : l_bound-0.5), 
		l_bound, nu_bound, 
		(altdis ? nu_bound + 0.75 : nu_bound + 0.5),
		((atma.na>1 || atmb.na>1) ? "\n" : " "));
	if (first_rst_output) {
		fprintf(fp,"\n%s ialtd=%d,\n",DEFAULT_RST, altdis ? 1: 0 );
		first_rst_output=0;
	}
	if (atma.na>1) {
		fprintf(fp," igr1=");
		for (i=0;i<atma.na;i++){
			fprintf(fp,"%4d,", atma.newanum[i] );
			if( i%10 == 9 ) fprintf(fp, "\n      ");
		}
		fprintf(fp,"\n");
	}
	if (atmb.na>1) {
		fprintf(fp," igr2=");
		for (i=0;i<atmb.na;i++){
			fprintf(fp,"%4d,", atmb.newanum[i] );
			if( i%10 == 9 ) fprintf(fp, "\n      ");
		}
		fprintf(fp,"\n");
	}

	fprintf(fp," /\n");
}

void translate_atoms( char *data )
{
#ifdef	DEBUG
	char 	*module="translate_atoms";
#endif
	int 	mapa, mapb;
	int 	i, j, residue;


	if ( !lookuprtyp(atma.rnum,atma.nrtyp) ||
	     !lookuprtyp(atmb.rnum,atmb.nrtyp)) {
		fprintf( stderr, 
		    "ERROR residue number not in molecule: %s\n", data );
		exit( -1 );
	}

	if ((mapa=mapdefinition(atma.atyp,atma.rtyp))<0){
		fprintf( stderr, "ERROR no map function for %s %s :data= %s\n",
			atma.atyp, atma.rtyp, data );
		exit( -1 );
	}
	if ((mapb=mapdefinition(atmb.atyp,atmb.rtyp))<0){
		fprintf( stderr, "ERROR no map function for %s %s :data= %s\n",
			atmb.atyp, atmb.rtyp, data );
		exit( -1 );
	}

	strcpy(atma.rtyp,atma.nrtyp);
	atma.na=map[mapa].na;
	atma.pseudoc=map[mapa].pseudoc;
	for (i=0; i < atma.na; i++ ) {
		strcpy(atma.newatyp[i],map[mapa].atrn[i]);
		if (map[mapa].rtrn[i] == 0) residue = atma.rnum;
			else residue = map[mapa].rtrn[i];
		if ((j=lookupanum(residue,atma.newatyp[i]))<0) {
			fprintf( stderr, 
				"ERROR atom (%s) not in residue (%d) : %s\n",
				atma.newatyp[i], residue, data );
			exit( -1 );
		}
		atma.newanum[i]=j;
	}

	strcpy(atmb.rtyp,atmb.nrtyp);
	atmb.na=map[mapb].na;
	atmb.pseudoc=map[mapb].pseudoc;
	for (i=0; i < atmb.na; i++ ) {
		strcpy(atmb.newatyp[i],map[mapb].atrn[i]);
		if (map[mapb].rtrn[i] == 0) residue = atmb.rnum;
			else residue = map[mapb].rtrn[i];
		if ((j=lookupanum(residue,atmb.newatyp[i]))<0) {
			fprintf( stderr, 
				"ERROR atom (%s) not in residue (%d) : %s\n",
				atmb.newatyp[i], residue, data );
			exit( -1 );
		}
		atmb.newanum[i]=j;
	}
#ifdef	DEBUG
	printf("%s: %s",module, data);
	printf( " nu_bound=%f\n", nu_bound);
	printf("\t%4d %-5s %-5s (%d) (%-5s:", 
		atma.rnum, atma.rtyp, atma.atyp, atma.na, atma.nrtyp);
	for (i=0;i<MAXDEGEN;i++) 
		printf("%4d %-5s,",atma.newanum[i],atma.newatyp[i]); 
	printf(")\n");
	printf("\t%4d %-5s %-5s (%d) (%-5s:", 
		atmb.rnum, atmb.rtyp, atmb.atyp, atmb.na, atmb.nrtyp);
	for (i=0;i<MAXDEGEN;i++) 
		printf("%4d %-5s,",atmb.newanum[i],atmb.newatyp[i]); 
	printf(")\n");
#endif
}

void correct_upperbound(char *data)
{
	int 	mapa, mapb;
#ifdef	DEBUG
	char 	*module="correct_upperbound";
#endif
	float	arg;

	if ((mapa=mapdefinition(atma.atyp,atma.rtyp))<0){
		fprintf( stderr, "ERROR no map function for %s %s :data= %s\n",
			atma.atyp, atma.rtyp, data );
		exit( -1 );
	}
	if ((mapb=mapdefinition(atmb.atyp,atmb.rtyp))<0){
		fprintf( stderr, "ERROR no map function for %s %s :data= %s\n",
			atmb.atyp, atmb.rtyp, data );
		exit( -1 );
	}

	if( no_corr){
		nu_bound = u_bound;
	} else { 
		arg = atma.na * atmb.na;
		nu_bound = u_bound * exp( log( arg )/6. ) + atma.pseudoc + atmb.pseudoc;
	}

#ifdef DEBUG
	printf("%s: %s",module, data);
	printf( "%s  %s %8.3f\n", atma.rtyp, atmb.rtyp, u_bound );
	printf( " nu_bound=%8.3f %d %d\n", nu_bound, atma.na, atmb.na );
#endif
}


static int lookupanum( int rnum, char *atyp )
{
	int 	i;

	for (i=0; i<natoms; i++) {
		if (mol[i].rnum==rnum && ew(mol[i].atyp,atyp) )
			return (mol[i].anum);
	}
	return( -1 );
}

static int lookuprtyp( int rnum, char *rtyp )
{
	int 	i;

	if( rnum == 0 ){
		strcpy( rtyp, "AMB" );
		return( 1 );
	}
	for (i=0; i<natoms; i++) {
		if (mol[i].rnum==rnum) {
			strcpy(rtyp, mol[i].rtyp);
			return( 1 );
		}
	}
	return ( 0 );
}

#define FMT "%s: atyp=%s rtyp=%s map.rtyp=%s, rmatch=%d mat.atyp=%s amatch=%d\n"

static int mapdefinition( char *atyp, char *rtyp )
{
	char 	*module="mapdefinition";
	int 	i;

	for (i=0; i<nmaps; i++ ) {
		if (debug) {
			printf( FMT,
				module, atyp, rtyp, map[i].rtyp,
				strncmp(rtyp,map[i].rtyp,strlen(map[i].rtyp)),
				map[i].atyp,
				identical(atyp,map[i].atyp) );
		}
		if ( strncmp(rtyp,map[i].rtyp,strlen(map[i].rtyp))==0 &&
		     identical(atyp,map[i].atyp) )
			return( i );
		if ( strncmp("AMB",map[i].rtyp,strlen(map[i].rtyp))==0 &&
		     identical(atyp,map[i].atyp) )
			return( i );
	}
	return ( -1 );
}

void initialize_data( char *data, char *comment )
{
	char 	*module="initialize_data";
	int 	i;

	if (strlen(data)==0) {
		fprintf( stderr, "%s: ERROR NOE data length 0!: (%s) (%s)\n",
			module, data, comment );
		exit( -1 );
	}

	/* initialize atoms */
	atma.rnum= -1;				atmb.rnum= -1;
	strcpy(atma.rtyp,"");			strcpy(atmb.rtyp,"");
	strcpy(atma.nrtyp,"");			strcpy(atmb.nrtyp,"");
	strcpy(atma.atyp,"");			strcpy(atmb.atyp,"");
	for (i=0;i<MAXDEGEN;i++) {
		strcpy(atma.newatyp[i],"");	strcpy(atmb.newatyp[i],"");
		atma.newanum[i]= 0;		atmb.newanum[i]= 0;
	}
	atma.na = -1;				atmb.na = -1;
	strcpy(atma.comment,"");		strcpy(atmb.comment,"");

	u_bound = nu_bound = BIGDIST;
}

void parse_upb( char *data, char *comment )
{
	char 	*module="parse_upb";
	int 	c;

	c=sscanf( data, "%d %s %s %d %s %s %f",
		&atma.rnum, atma.rtyp, atma.atyp,
		&atmb.rnum, atmb.rtyp, atmb.atyp,
		&u_bound );

	if ( c!=7 ) {
		fprintf( stderr, 
		  "%s: ERROR NOE data not in 7(%d) column format: (%s) (%s)\n",
			module, c, data, comment );
		exit( -1 );
	}

	c = sscanf( comment, "%*s %d:%d", &xpk, &npk );
#ifdef DEBUG
	printf( "%s: %d %s %s %d %s %s %f %d\n",
		module,
		atma.rnum, atma.rtyp, atma.atyp,
		atmb.rnum, atmb.rtyp, atmb.atyp,
		u_bound, xpk);
#endif
}

void parse_ual( char *data, char *comment )
{
	char 	*module="parse_ual";
	int 	c;

	c=sscanf( data, "%d %s %s %d %s %s %f %f",
		&atma.rnum, atma.rtyp, atma.atyp,
		&atmb.rnum, atmb.rtyp, atmb.atyp,
		&l_bound, &u_bound );

	if ( c!=8 ) {
		fprintf( stderr, 
		  "%s: ERROR UAL data not in 8(%d) column format: (%s) (%s)\n",
			module, c, data, comment );
		exit( -1 );
	}

	c = sscanf( comment, "%*s %d:%d", &xpk, &npk );
#ifdef DEBUG
	printf( "%s: %d %s %s %d %s %s %f %f\n",
		module,
		atma.rnum, atma.rtyp, atma.atyp,
		atmb.rnum, atmb.rtyp, atmb.atyp,
		l_bound, u_bound);
#endif
}

void parse_vol( char *data, char *comment )
{
	char 	*module="parse_vol";
	int 	c;

	c=sscanf( data, "%d %s %s %d %s %s %f",
		&atma.rnum, atma.rtyp, atma.atyp,
		&atmb.rnum, atmb.rtyp, atmb.atyp,
		&volume );

	if ( c!=7 ) {
		fprintf( stderr, 
		  "%s: ERROR VOL data not in 7(%d) column format: (%s) (%s)\n",
			module, c, data, comment );
		exit( -1 );
	}
}

static int read_data( FILE *fp, char *data, char *comment )
{
#ifdef	DEBUG
	char 	*module="read_data";
#endif
	char 	line[LINESIZE], crud[LINESIZE]; 
	int 	slen;

	strcpy(data,"");
	strcpy(comment,"");
	while (strlen(data)==0){
		if (!getline) {
			if (report) printf("# no more data.\n");
			return(0);
		}
		if (lineis("#") || blankline || sscanf(line,"%s", crud )<=0){
			if ( print_rm6 && ! blankline )
				fprintf(fp_rm6, "%s", line );
			if (report) printf( "%s", line );
			continue;
		}
		sscanf( line, "%[^#]%[^\n]", data, comment );
		if (sscanf(data,"%s", crud )<=0){
			strcpy(data,"");
			if (report) printf( "%s", line );
			continue;
		}
	}
	return(1);
}

void leftcharwrap( char *str, int n )
{
	char 	t[WORDSIZE];
	int 	l = strlen(str);

	strncpy(t,str,n);	/* pulls the prefix off */
	t[n] = 0;		/* terminate, 0 or '\0' */
	strcpy(str,&str[n]);
	strcpy(&str[l-n],t);
}

static int read_pdb_file( FILE *fp )
{
	char 	*module="read_structure";
	int 	nat, atom_serial_number, residue_sequence_number, lineno = 0;
	char 	line[LINESIZE];
	char	atom[WORDSIZE], atom_name[WORDSIZE], residue_name[WORDSIZE];
	char 	l[WORDSIZE+WORDSIZE+WORDSIZE];
	float 	x, y, z;

#ifdef DEBUG
	printf( "%s: begin;\n", module );
#endif
	nat=0;
 	while ( getline ) {
		lineno++;
		line[21] = ' '; /* remove chain-id if it is present */
		if ( lineis("ATOM") || lineis("HETATM") ){
			sscanf( line, "%s %d %s %s %d %f %f %f",
				atom, &atom_serial_number, atom_name,
				residue_name, &residue_sequence_number,
				&x, &y, &z );
			mol[nat].anum=atom_serial_number;
			mol[nat].rnum=residue_sequence_number;
			sscanf( residue_name, "%s", mol[nat].rtyp );
			sscanf( atom_name, "%s", mol[nat].atyp );

			if (isdigit(*mol[nat].atyp))
				leftcharwrap(mol[nat].atyp,1);

			mol[nat].xyz[X]=x;
			mol[nat].xyz[Y]=y;
			mol[nat].xyz[Z]=z;
			sprintf(l, "%3d %-5s %-5s",
				mol[nat].rnum, mol[nat].rtyp, mol[nat].atyp );
			strcpy(mol[nat].name, l);

			if (debug) {
				printf( "%s: %-6s%5d %-4s %-3s  %4d ",
					module, atom, 
					mol[nat].anum, mol[nat].atyp,
					mol[nat].rtyp, mol[nat].rnum );
				printf( "   %8.3f%8.3f%8.3f\n",
					mol[nat].xyz[X],
					mol[nat].xyz[Y],
					mol[nat].xyz[Z]);
			}
			nat++;
			assert (nat < MAXATOMS);
		} else if ( lineis("HEADER") || lineis("COMPND") ||
		            lineis("SOURCE") || lineis("AUTHOR") ) {
			if (report) printf("%s",line);
		}
	}
	if (report) printf( "%s: %d atoms read;\n", module, nat );
	maxres = mol[nat].rnum;
	fclose( fp );
	return( nat );
}

static int read_maplib( FILE *fp )
{
	char 	*module="read_maplib";
	int 	i, j, nm, nw, lineno = 0;
	int		k, ij, npseudo, found;
	char	tmpatyp[WORDSIZE];
	char	tmprtyp[WORDSIZE];
	char 	line[LINESIZE];
	char 	rtyp[WORDSIZE];
	char 	word[MAXWORDS][WORDSIZE];

#ifdef DEBUG
	printf( "%s: begin;\n", module );
#endif
	nm=0;
 	while ( getline ) {
		lineno++;

		if ( lineis("RESIDUE") ){
			sscanf( line, "%*s %s", rtyp );
			continue;
		}

		if ( lineis("MAPPING") ){ 
			nw=split(line,word);
			assert( nw <= MAXWORDS );
			sscanf( rtyp,    "%s", map[nm].rtyp );
			sscanf( word[1], "%s", map[nm].atyp );
			for ( i = 0; i<nw-3; i++ ){
				sscanf(word[i+3], "%s", map[nm].atrn[i]);
				map[nm].rtrn[i] = 0;
			}
			map[nm].pseudoc=0.0;
			map[nm].na=nw-3;
#ifdef DEBUG
			printf ( "%s:MAP %-5s %-5s to ",
				module, map[nm].rtyp, map[nm].atyp );
			for ( i = 0; i<map[nm].na; i++ )
				printf( "%-5s", map[nm].atrn[i]);
			printf( "\n" );
#endif
			nm++; npseudo = nm;
			assert (nm < MAXMAPLIB);
		} else if ( lineis("PSEUDOC") ) {
			nw=split(line,word);
			sscanf( rtyp,    "%s", map[nm].rtyp );
			sscanf( word[1], "%s", map[nm].atyp );
			sscanf( word[3], "%s", map[nm].atrn[0] );
			map[nm].rtrn[0] = 0;
			sscanf( word[4], "%f", &map[nm].pseudoc );
			map[nm].na=1;
#ifdef DEBUG
			printf ( "%s:PSE %-5s %-5s to ",
				module, map[nm].rtyp, map[nm].atyp );
			for ( i = 0; i<map[nm].na; i++ )
				printf( "%-5s", map[nm].atrn[i]);
			printf( "%7.2f", map[nm].pseudoc );
			printf( "\n" );
#endif
			nm++; npseudo = nm;
			assert (nm < MAXMAPLIB);
		} else if ( lineis("AMBIG") ) {
			nw=split(line,word);
			assert( nw <= MAXWORDS );
			sscanf( "AMB",   "%s", map[nm].rtyp );
			sscanf( word[1], "%s", map[nm].atyp );
			for ( ij=0, i=0; i<(nw-3)/2; i++ ){
				sscanf(word[2*i+3], "%s", tmpatyp );
				sscanf(word[2*i+4], "%d", &idum );

/*              cycle over all "standard" mappings (which are defined
                by this point, and re-expand any keywords found          */

				if ( !lookuprtyp( idum, tmprtyp )) {
					fprintf( stderr, 
		    			"ERROR residue number not in molecule: %d\n", idum );
					exit( -1 );
				}
				for( found=0, j=0; j<npseudo; j++ ){
					if( strcmp( tmpatyp, map[j].atyp ) == 0 &&
					    strcmp( tmprtyp, map[j].rtyp ) == 0 ){
						for( k=0; k<map[j].na; k++ ){
							assert( ij < MAXDEGEN );
							strcpy( map[nm].atrn[ij], map[j].atrn[k] );
							map[nm].rtrn[ij] = idum;
							ij++; 
						}
						found = 1; break;
					}
				}
				if( !found ){
					assert( ij < MAXDEGEN );
					strcpy( map[nm].atrn[ij], tmpatyp );
					map[nm].rtrn[ij] = idum;
					ij++;
				}
			}
			map[nm].na = ij;
#ifdef DEBUG
			printf ( "%s:%-5s %-5s to ",
				module, map[nm].rtyp, map[nm].atyp );
			for ( i = 0; i<map[nm].na; i++ )
				printf( "(%-4s %3d); ", map[nm].atrn[i], map[nm].rtrn[i]);
			printf( "\n" );
#endif
			nm++;
			if (nm >= MAXMAPLIB){
				fprintf(stderr, 
					"ERROR Too many map definitions, max is %d\n",
					MAXMAPLIB );
				exit( -1 );
			}
		}
	}
	if (report) printf( "%s: %d mappings read;\n", module, nm );
	fclose( fp );
	return( nm );
}

static int split( char *str, char wordptr[][WORDSIZE] )
{
        int 	len, nwords;

        len = nwords = 0;
        while (*str != '\n' && *str != '\0') {
                if ( nwords > MAXWORDS) {
			fprintf(stderr, 
			"ERROR Too many words in split function, input\n::>%s<::\n",str);
			exit( -1 );
                        /*return(-1);*/
		}
                /* following takes care of leading white space */
                while (*str == ' ' || *str == '\t') str++;
                while (*str != ' ' && *str != '\n' && *str != '\0' 
				   && *str != '\t')
                        wordptr[nwords][len++] = *str++;
                wordptr[nwords++][len] = '\0';
                len = 0;
                /* following takes care of trailing white space, (maybe) */
                while (*str == ' ' || *str == '\t') str++;
        }
        return(nwords);
}

int main( int argc, char *argv[] )
{
	int 	i, j, natles;
	char 	comment[LINESIZE];
	char    data[LINESIZE];
	char	*ambh;

	if (argc<=1) {
		for ( i = 0; help[i] != 0; i++ ) fprintf( stderr, "%s", help[i] );
		exit( -1 );
	}

	if (argc>1 && strcmp(argv[1],"-help")==0 ) {
		for ( i = 0; banner[i] != 0; i++ ) fprintf( stderr, "%s", banner[i] );
		for ( i = 0; help[i] != 0; i++ ) fprintf( stderr, "%s", help[i] );
		exit( -1 );
	}

	for ( i = 0; banner[i] != 0; i++ ) fprintf( stderr, "%s", banner[i] );
	fprintf(stderr, "Currently configured for up to %d atoms\n", MAXATOMS);

	if((ambh = getenv( "AMBERHOME" ))){
		strncpy( f_map, ambh, MAXFILENAME-1 );
		strcat( f_map, "/dat/" );
	}else{
		strcpy( f_map, "./" );
	}
	strcat(f_map,DEFAULT_MAP);

	read_map=1;
	for ( i=1; i < argc; i++ ) {
		if( strcmp(argv[i],"-report")==0 ) {
			report=1;
		} else if( strcmp(argv[i],"-debug")==0 ) {
			debug=1;
		} else if( strcmp(argv[i],"-nocorr")==0 ) {
			no_corr=1;
		} else if( strcmp(argv[i],"-altdis")==0 ) {
			altdis=1;
		} else if( strcmp(argv[i],"-ual")==0 ) {
			i++;
			sscanf(argv[i],"%s",f_ual);
			if (report) printf( "ual file %s\n", f_ual );
			read_ual=1;
		} else if( strcmp(argv[i],"-les")==0 ) {
			i++;
			sscanf(argv[i],"%s",f_les);
			if (report) printf( "les file %s\n", f_les );
			read_les=1;
		} else if( strcmp(argv[i],"-upb")==0 ) {
			i++;
			sscanf(argv[i],"%s",f_upb);
			if (report) printf( "upb file %s\n", f_upb );
			read_upb=1;
		} else if( strcmp(argv[i],"-vol")==0 ) {
			i++;
			sscanf(argv[i],"%s",f_vol);
			if (report) printf( "vol file %s\n", f_vol );
			read_vol=1;
		} else if( strcmp(argv[i],"-pdb")==0 ) {
			i++;
			sscanf(argv[i],"%s",f_pdb);
			if (report) printf( "pdb file %s\n", f_pdb );
			read_pdb=1;
		} else if( strcmp(argv[i],"-map")==0 ) {
			i++;
			sscanf(argv[i],"%s",f_map);
			if (report) printf( "map file %s\n", f_map );
			read_map=1;
		} else if( strcmp(argv[i],"-rst")==0 ) {
			i++;
			sscanf(argv[i],"%s",f_rst);
			if (report) printf( "rst file %s\n", f_rst );
			print_rst=1;
		} else if( strcmp(argv[i],"-dgm")==0 ) {
			i++;
			sscanf(argv[i],"%s",f_dgm);
			if (report) printf( "dgm file %s\n", f_dgm );
			print_dgm=1;
		} else if( strcmp(argv[i],"-rm6")==0 ) {
			i++;
			sscanf(argv[i],"%s",f_rm6);
			if (report) printf( "rm6 file %s\n", f_rm6 );
			print_rm6=1;
		} else if( strcmp(argv[i],"-svf")==0 ) {
			i++;
			sscanf(argv[i],"%s",f_svf);
			if (report) printf( "svf file %s\n", f_svf );
			print_svf=1;
		} else if( strcmp(argv[i],"-imix")==0 ) {
			i++;
			sscanf(argv[i],"%d",&imix);
			if (report) printf( "imix= %d\n", imix );
		} else {
			fprintf( stderr, 
				"ERROR Argument (%s) not recognized.\n", 
				argv[i] );
			fatal_error++;
		}
	}

	/* check input file names for readability */
	if ( !read_upb && !read_vol && !read_ual ) {
		fprintf( stderr, "ERROR, NOE, UAL or VOL file not specified.\n" );
		fatal_error++;
	} else if ( read_upb && read_vol ) {
		fprintf( stderr, 
			"ERROR, cannot specify both NOE and VOL file.\n" );
		fatal_error++;
	} else if ( read_ual && (fp_ual=fopen(f_ual,"r"))==NULL) {
		fprintf( stderr, 
			"ERROR can't open ual file %s for reading.\n", f_ual );
		fatal_error++;
	} else if ( read_upb && (fp_upb=fopen(f_upb,"r"))==NULL) {
		fprintf( stderr, 
			"ERROR can't open noe file %s for reading.\n", f_upb );
		fatal_error++;
	} else if ( read_vol && (fp_vol=fopen(f_vol,"r"))==NULL) {
		fprintf( stderr, 
			"ERROR can't open vol file %s for reading.\n", f_vol );
		fatal_error++;
	}

	if (!read_pdb) {
		fprintf( stderr, "ERROR no PDB file specified.\n" );
		fatal_error++;
	} else if ( (fp_pdb=fopen(f_pdb,"r"))==NULL) {
		fprintf( stderr, 
			"ERROR can't open pdb file %s for reading.\n", f_pdb );
		fatal_error++;
	}

	if (read_les) {
		if ( (fp_les=fopen(f_les,"r"))==NULL) {
			fprintf( stderr, 
				"ERROR can't open les file %s for reading.\n", f_les );
			fatal_error++;
		}
	}

	if (!read_map) {
		fprintf( stderr, "ERROR no MAP file specified.\n" );
		fatal_error++;
	} else if ( (fp_map=fopen(f_map,"r"))==NULL) {
		fprintf( stderr, 
			"ERROR can't open map file %s for reading.\n", f_map );
		fatal_error++;
	}

	/* check output file names for writability */
	if (!print_rst) {
		fprintf( stderr, "Caution, output RST file not specified.\n" );
	} else if ( (fp_rst=fopen(f_rst,"w"))==NULL) {
		fprintf( stderr, 
			"ERROR can't open RST file %s for writing.\n", f_rst );
		fatal_error++;
	}

	if (!print_rm6) {
/*
		fprintf( stderr, 
			"Caution, output RM6NOE file not specified.\n" );
*/
	} else if ( (fp_rm6=fopen(f_rm6,"w"))==NULL) {
		fprintf( stderr, 
			"ERROR can't open RM6NOE file %s for writing.\n", 
			f_rm6 );
		fatal_error++;
	}

	if (!print_dgm) {
/*
		fprintf( stderr, 
			"Caution, output DGEOMIN file not specified.\n" );
*/
	} else if ( (fp_dgm=fopen(f_dgm,"w"))==NULL) {
		fprintf( stderr, 
			"ERROR can't open DGEOMIN file %s for writing.\n", 
			f_dgm );
		fatal_error++;
	}

	if ( read_vol && ! print_svf ) {
		fprintf( stderr, "ERROR VOL input requires SVF output.\n" );
		fatal_error++;
	} else if ( read_vol && (fp_svf=fopen(f_svf,"w"))==NULL ) {
		fprintf( stderr, 
			"ERROR can't open SVF file %s for writing.\n", f_svf );
		fatal_error++;
	}

	if (fatal_error) {
		exit( -1 );
	}

/*   Get here if command line all seems to be OK               */

	fprintf( stderr, "Using MAP file %s\n", f_map );

	natoms=read_pdb_file(fp_pdb);
	if (report) 
		printf( "Number of atoms read in PDB file %s = %d\n", 
						f_pdb, natoms );

	nmaps=read_maplib(fp_map);
	if (report) 
		printf("Number of mapping functions read in MAP file %s = %d\n",
						f_map, nmaps );

	if( read_les ){
		fscanf( fp_les, "%d", &natles );
		if( natles != natoms ){
			fprintf( stderr, "Bad number of atoms in les file: %d\n", natles );
			exit( 1 );
		}
		for( i=1; i<=natoms; i++ ){
			fscanf( fp_les, "%d", &invcnt[ i ] );
			for( j=0; j<invcnt[i]; j++ ){
				fscanf( fp_les, "%d", &invcop[ i ][ j ] );
			}
		}
	}

	if ( read_upb ) {
		while ( read_data(fp_upb,data,comment) ) {
			if ( !(strncmp(data,"END OF SET",10)==0) ) {
				initialize_data( data, comment );
				parse_upb(data,comment);
				translate_atoms(data);
				correct_upperbound(data);
				l_bound = 1.8;
				if (print_rst) write_rst(fp_rst,data,comment);
				if (print_rm6) write_rm6(fp_rm6,data,comment);
				if (print_dgm) write_dgm(fp_dgm,data,comment);
			} else {
				if (print_rm6) write_rm6(fp_rm6,data,comment);
			}
		}
	} else if ( read_ual ) {
		while ( read_data(fp_ual,data,comment) ) {
			if ( !(strncmp(data,"END OF SET",10)==0) ) {
				initialize_data( data, comment );
				parse_ual(data,comment);
				translate_atoms(data);
				correct_upperbound(data);
				if (print_rst) write_rst(fp_rst,data,comment);
				if (print_rm6) write_rm6(fp_rm6,data,comment);
				if (print_dgm) write_dgm(fp_dgm,data,comment);
			} else {
				if (print_rm6) write_rm6(fp_rm6,data,comment);
			}
		}
	} else if ( read_vol ) {
		fprintf( fp_svf, " &noeexp\n" );
		while ( read_data(fp_vol,data,comment) ) {
			initialize_data( data, comment );
			parse_vol( data, comment );
			translate_atoms(data);
			if (print_svf) write_svf(fp_svf,data,comment);
		}
		fprintf( fp_svf, " npeak(%d) = %d, -1\n /\n", imix, npeak );
		fprintf( fp_svf, " entire molecule is sub-molecule\n" );
		fprintf (fp_svf, "RES  1 %d\nEND\nEND\n", maxres );
		fprintf (fp_svf, " &noeexp   npeak=-1,   /\n" );
	}

	if (print_rst) fclose(fp_rst);
	if (print_rm6) fclose(fp_rm6);
	if (print_dgm) fclose(fp_dgm);
	if (print_svf) fclose(fp_svf);

	return(0);
}
/* END MAIN */
