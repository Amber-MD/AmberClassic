/*
** nef_to_RST converts a NMR-STAR (NEF) restraint list into the raw
**    files used by sander and pmemd.
**
**    written by Garry P. Gippert, TSRI.
**    LES and AMBIGUOUS constraint modifications by D.A. Case.
**    NEF ability added by DAC
*/

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "../../cifparse/cifparse.h"
#include "global.h"

#define DEFAULT_RST    ("      rk2=20.0, rk3=20.0, ir6=1,")
#define DEFAULT_MAP ("map.NEF-AMBER")

#define R1_R4_RANGE    1.0
#define K_FORCE_LB    2.0
#define K_FORCE_UB    2.0

int        idum; 

struct molecule {
    int    anum;
    char   atyp[WORDSIZE];
    char   rtyp[WORDSIZE];
    int    rnum;
} mol[MAXATOMS];

struct residues {
    int    rnum;
    char   cid;
    int    amber_rnum;
} residue_list[1000];

int     natoms;

#define MAXMAPLIB 13500

struct mapping {
    char    cid;                       /* chainId */
    char    rtyp[WORDSIZE];            /* residue name */
    char    atyp[WORDSIZE];            /* atom or pseudo-atom name */
    char    atrn[MAXDEGEN][WORDSIZE];  /* expanded (real) atom names */
    int        rtrn[MAXDEGEN];            /* residue number of expanded atom;
                                          if zero, use the number of the
                                          resiude being considered   */      
    int    na;                            /* number of expanded atoms  */
    float    pseudoc;                   /* pseudoatom distance correction */
} map[MAXMAPLIB];

int     nm;

struct pseudoatom {
    int    rnum;
    char    rtyp[WORDSIZE];
    char    nrtyp[WORDSIZE];
    char    atyp[WORDSIZE];
    char    newatyp[MAXDEGEN][WORDSIZE];
    int    newanum[MAXDEGEN];
    int    na;
    char    comment[COMMENTSIZE];
    float   pseudoc;
} atma, atmb;

int        xpk,npk;
int        imix=1;
int        npeak=0;
int        maxres;  /* will hold number of residues in the system */
int     npseudo; /* number of pseudoatom mappings, from MAP file */
float    u_bound;
float    l_bound;
float    volume;
float    nu_bound;
float    nl_bound;

#define MAXCOPY 10
int        invcnt[ MAXATOMS ];
int        invcop[ MAXATOMS ][ MAXCOPY ];

static char *banner[] = {
    "# nef_to_RST\n",
    0
};

static char *help[] = {
    "convert NEF restraints to Amber format\n",
    "\n",
    "input:\n",
    "   -nef <filename>: NEF file\n",
    "   -pdb <filename>: PDBFILE using AMBER nomenclature and numbering\n",
    "   -p4a <filename>: _renum file from pdb4amber\n",
    "   -map <filename>: MAP file  (default:map.NEF-AMBER)\n",
    "\n",
    "output\n",
    "   -rst <filename>:  SANDER DISANG format\n",
    "   -rdc <filename>:  SANDER DIP format\n",
    "   \n",
    "other options:\n",
    "   -nocorr (do not correct upper bound for r**-6 averaging)\n",
    "   -altdis (use alternative form for the distance restraints)\n",
    "   -help   (gives you this explanation, overrides other parameters)\n",
    "   -report (gives you short runtime diagnostic output)\n",
    "\n",
    "errors come to stderr.\n",
    "\n",
    0
};

#define MAXFILENAME 128

char    f_nef[MAXFILENAME];
char    f_pdb[MAXFILENAME];
char    f_map[MAXFILENAME];
char    f_rst[MAXFILENAME];
char    f_rdc[MAXFILENAME];
char    f_p4a[MAXFILENAME];   /* _renum file for pdb4amber  */

int     report=0;
int     debug=0;
int     read_nef=0;
int     read_map=0;
int     read_pdb=0;
int     read_p4a=0;
int     print_rst=0;
int     print_rdc=0;
int     fatal_error=0;
int     first_rst_output=1;
int        no_corr = 0;
int        altdis = 0;
int        nambig = 0;

FILE     *fp_nef, *fp_pdb, *fp_map, *fp_rst, *fp_rdc, *fp_p4a;

int lookupanum( int rnum, char *atyp )
{
    int     i;

    for (i=0; i<natoms; i++) {
        if (mol[i].rnum==rnum && ew(mol[i].atyp,atyp) )
            return (mol[i].anum);
    }
    return( -1 );
}

int lookuprtyp( int rnum, char *rtyp )
{
    int     i;

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

int residue_renumber( char cid, int rnum ){
    int found, ires, amber_rnum;

    found = 0;
    for (ires=1; ires<=maxres; ires++ ){
        if( cid == residue_list[ires].cid 
            && rnum == residue_list[ires].rnum ) { 
            found=1;
            amber_rnum = residue_list[ires].amber_rnum;
            break;
        }
    }
    assert( found );
    return( amber_rnum);
}

#define FMT "%s: atyp=%s rtyp=%s map.rtyp=%s, rmatch=%d mat.atyp=%s amatch=%d\n"

int mapdefinition( char *atyp, char *rtyp )
{
    char     *module="mapdefinition";
    int     i;

    for (i=0; i<nm; i++ ) {
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

int split( char *str, char wordptr[][WORDSIZE] )
{
        int     len, nwords;

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

void write_rst( FILE *fp, char *data, char *comment)
{
    int     i,j,k,dup;
#ifdef DEBUG
    char     *module="write_rst";
    printf("%s: (%s) (%s)\n",module, data, comment);
#endif

/*    First check to see that this constraint never refers to itself,
    i.e., that all the atoms "a" are distinct from the atoms "b"    */

    for( i=0; i< atma.na; i++ ){
        for( j=0; j< atmb.na; j++ ){
            if( atma.newanum[i] == atmb.newanum[j] ){
                fprintf(fp,"# This appears to be self-referential: ignoring\n" );
                fprintf(fp,"# details: i,j, atma[i], atmb[j]:  %d %d   %d %d\n\n",
                    i,j, atma.newanum[i],atmb.newanum[j] );
                return;
            }
        }
    }

    fprintf(fp," &rst\n  ixpk= %d, nxpk= %d,", xpk,npk);

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
            /* check for duplicates: */
            dup  = 0;
            for( k=0; k<i; k++ )
                if( atma.newanum[k] == atma.newanum[i] ) dup = 1;
            if( dup ) continue;
            fprintf(fp,"%4d,", atma.newanum[i] );
            if( i%10 == 9 ) fprintf(fp, "\n      ");
        }
        fprintf(fp,"\n");
    }
    if (atmb.na>1) {
        fprintf(fp," igr2=");
        for (i=0;i<atmb.na;i++){
            /* check for duplicates: */
            dup  = 0;
            for( k=0; k<i; k++ )
                if( atmb.newanum[k] == atmb.newanum[i] ) dup = 1;
            if( dup ) continue;
            fprintf(fp,"%4d,", atmb.newanum[i] );
            if( i%10 == 9 ) fprintf(fp, "\n      ");
        }
        fprintf(fp,"\n");
    }

    fprintf(fp," /\n\n");
}

void translate_atoms( char *data )
{
#ifdef    DEBUG
    char     *module="translate_atoms";
#endif
    int     mapa, mapb;
    int     i, j, residue;


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
#ifdef    DEBUG
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
    int     mapa, mapb;
#ifdef    DEBUG
    char     *module="correct_upperbound";
#endif
    float    arg;

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

void initialize_data( char *data, char *comment )
{
    char     *module="initialize_data";
    int     i;

    if (strlen(data)==0) {
        fprintf( stderr, "%s: ERROR NOE data length 0!: (%s) (%s)\n",
            module, data, comment );
        exit( -1 );
    }

    /* initialize atoms */
    atma.rnum= -1;                atmb.rnum= -1;
    strcpy(atma.rtyp,"");            strcpy(atmb.rtyp,"");
    strcpy(atma.nrtyp,"");            strcpy(atmb.nrtyp,"");
    strcpy(atma.atyp,"");            strcpy(atmb.atyp,"");
    for (i=0;i<MAXDEGEN;i++) {
        strcpy(atma.newatyp[i],"");    strcpy(atmb.newatyp[i],"");
        atma.newanum[i]= 0;        atmb.newanum[i]= 0;
    }
    atma.na = -1;                atmb.na = -1;
    strcpy(atma.comment,"");        strcpy(atmb.comment,"");

    u_bound = nu_bound = BIGDIST;
}

void parse_upb( char *data, char *comment )
{
    char     *module="parse_upb";
    int     c;

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

void parse_ambig( char *line )
{
    int nw, i, ij, j, k, idum, found;
    char     word[MAXWORDS][WORDSIZE];
    char    tmpatyp[WORDSIZE];
    char    tmprtyp[WORDSIZE];

    // fprintf( stderr, "in parse_ambig: nm=%d, %s\n", nm,line );

    nw=split(line,word);
    assert( nw <= MAXWORDS );
    sscanf( "AMB",   "%s", map[nm].rtyp );
    sscanf( word[1], "%s", map[nm].atyp );
    for ( ij=0, i=0; i<(nw-3)/2; i++ ){
        sscanf(word[2*i+3], "%s", tmpatyp );
        sscanf(word[2*i+4], "%d", &idum );

/*              cycle over all "standard" mappings (which are defined
        by this point), and re-expand any keywords found          */

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
    /* more info as comments in the output file: */
    fprintf ( fp_rst, "# %-5s %-5s to ",
        map[nm].rtyp, map[nm].atyp );
    for ( i = 0; i<map[nm].na; i++ )
        fprintf( fp_rst, "(%-4s %3d); ", map[nm].atrn[i], map[nm].rtrn[i]);
    fprintf( fp_rst, "\n" );
#endif

    nm++;
    if (nm >= MAXMAPLIB){
        fprintf(stderr, 
            "ERROR Too many map definitions, max is %d\n",
            MAXMAPLIB );
        exit( -1 );
    }
}

void leftcharwrap( char *str, int n )
{
    char     t[WORDSIZE];
    int     l = strlen(str);

    strncpy(t,str,n);    /* pulls the prefix off */
    t[n] = 0;        /* terminate, 0 or '\0' */
    strcpy(str,&str[n]);
    strcpy(&str[l-n],t);
}

int read_pdb_file( FILE *fp )
{
    char     *module="read_structure";
    int     nat, lineno = 0;
    char     line[LINESIZE];

#ifdef DEBUG
    printf( "%s: begin;\n", module );
#endif
    nat=0;
     while ( getline ) {
        lineno++;
        line[16] = ' '; /* remove alt conformer if it is present */
        line[26] = '\0'; /* this will end the residue number field */
        if ( lineis("ATOM") || lineis("HETATM") ){
            /* safe to use sscanf up to resdiue name, since blanks should
               always be present: */
            sscanf( &line[7], "%d %s %s",
                &mol[nat].anum, mol[nat].atyp, mol[nat].rtyp);
            mol[nat].rnum = atoi( &line[22] );
            if (isdigit(*mol[nat].atyp))
                leftcharwrap(mol[nat].atyp,1);

            if (debug) {
                printf( "%s: %5d %-4s %-3s  %4d ",
                    module, mol[nat].anum, mol[nat].atyp,
                    mol[nat].rtyp, mol[nat].rnum );
            }
            nat++;
            assert(nat < MAXATOMS);
        } else if ( lineis("HEADER") || lineis("COMPND") ||
                    lineis("SOURCE") || lineis("AUTHOR") ) {
            if (report) printf("%s",line);
        }
    }
    if (report) printf( "%s: %d atoms read;\n", module, nat );
    fclose( fp );
    return( nat );
}

void read_p4a_file( FILE *fp )
{
    char     *module="read_renum";
    int      lineno = 0;
    char     line[LINESIZE];

#ifdef DEBUG
    printf( "%s: begin;\n", module );
#endif
     while ( getline ) {
        lineno++;
        residue_list[lineno].cid = line[4];
        line[11] = '\0';
        residue_list[lineno].rnum = atoi( &line[8] );
        line[24] = '\0';
        residue_list[lineno].amber_rnum = atoi( &line[21] );
        //fprintf( stderr, "%c  %3d  %3d\n", residue_list[lineno].id,
        //      residue_list[lineno].rnum, residue_list[lineno].amber_rnum );
    }
    maxres = lineno; 
    if (report) printf( "%s: %d residue renumberings read;\n", module, lineno );
    fclose( fp );
}

int read_maplib( FILE *fp )
{
    char     *module="read_maplib";
    int     i, nw, lineno = 0;
    char     line[LINESIZE];
    char     rtyp[WORDSIZE];
    char     word[MAXWORDS][WORDSIZE];

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
            parse_ambig( line );
        }
    }
    if (report) printf( "%s: %d mappings read;\n", module, nm );
    fclose( fp );
    return( nm );
}

int main( int argc, char *argv[] )
{
    int     i, j, nmaps, anum[4];
    char     comment[LINESIZE], data[LINESIZE];

    char    *ambh;
    char    *ambig1=NULL;
    char    *ambig2=NULL;
    int iblock, stat, iCat, is_ambig, curr_id, next_id;
    int col1, col2, col3, col4, col5, col6, col7, col8, col9, 
        col10, col11, col12, col13, col14, col15, col16;
    int tot, dataset;
    int amber_rnum1, amber_rnum2, amber_rnum3, amber_rnum4;
    float r1, r2, r3, r4, dobs, gigj, dij;

    if (argc<=1) {
        for ( i = 0; help[i] != 0; i++ ) fprintf( stderr, "%s", help[i] );
        exit( -1 );
    }

    if (argc>1 && strcmp(argv[1],"-help")==0 ) {
        for ( i = 0; banner[i] != 0; i++ ) fprintf( stderr, "%s", banner[i] );
        for ( i = 0; help[i] != 0; i++ ) fprintf( stderr, "%s", help[i] );
        exit( -1 );
    }

    // for ( i = 0; banner[i] != 0; i++ ) fprintf( stderr, banner[i] );
    // fprintf(stderr, "Currently configured for up to %d atoms\n", MAXATOMS);

    if(( ambh = getenv( "AMBERHOME" ))){
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
        } else if( strcmp(argv[i],"-nef")==0 ) {
            i++;
            sscanf(argv[i],"%s",f_nef);
            if (report) printf( "nef file %s\n", f_nef );
            read_nef=1;
        } else if( strcmp(argv[i],"-pdb")==0 ) {
            i++;
            sscanf(argv[i],"%s",f_pdb);
            if (report) printf( "pdb file %s\n", f_pdb );
            read_pdb=1;
        } else if( strcmp(argv[i],"-p4a")==0 ) {
            i++;
            sscanf(argv[i],"%s",f_p4a);
            if (report) printf( "p4a file %s\n", f_p4a );
            read_p4a=1;
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
        } else if( strcmp(argv[i],"-rdc")==0 ) {
            i++;
            sscanf(argv[i],"%s",f_rdc);
            if (report) printf( "rdc file %s\n", f_rdc );
            print_rdc=1;
        } else {
            fprintf( stderr, 
                "ERROR Argument (%s) not recognized.\n", 
                argv[i] );
            fatal_error++;
        }
    }

    /* check input file names for readability */
    if ( !read_nef ) {
        fprintf( stderr, "ERROR: NEF file not specified.\n" );
        fatal_error++;
    } else if ( (fp_nef=fopen(f_nef,"r"))==NULL) {
        fprintf( stderr, 
            "ERROR can't open nef file %s for reading.\n", f_nef );
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

    if (!read_map) {
        fprintf( stderr, "ERROR no MAP file specified.\n" );
        fatal_error++;
    } else if ( (fp_map=fopen(f_map,"r"))==NULL) {
        fprintf( stderr, 
            "ERROR can't open map file %s for reading.\n", f_map );
        fatal_error++;
    }

    if (!read_p4a) {
        fprintf( stderr, "ERROR no P4A file specified.\n" );
        fatal_error++;
    } else if ( (fp_p4a=fopen(f_p4a,"r"))==NULL) {
        fprintf( stderr, 
            "ERROR can't open p4a file %s for reading.\n", f_p4a );
        fatal_error++;
    }

    /* check output file names for writability */
    if (!print_rst) {
        fprintf( stderr, "ERROR: output RST file not specified.\n" );
        fatal_error++;
    } else if ( (fp_rst=fopen(f_rst,"w"))==NULL) {
        fprintf( stderr, 
            "ERROR can't open RST file %s for writing.\n", f_rst );
        fatal_error++;
    }

    if (!print_rdc) {
        fprintf( stderr, "ERROR: output RDC file not specified.\n" );
        fatal_error++;
    } else if ( (fp_rdc=fopen(f_rdc,"w"))==NULL) {
        fprintf( stderr, 
            "ERROR can't open RDC file %s for writing.\n", f_rdc );
        fatal_error++;
    }

    if (fatal_error) exit( 1 );

    /*   Get here if command line all seems to be OK  */

    natoms=read_pdb_file(fp_pdb);
    if (report) 
        printf( "Number of atoms read in PDB file %s = %d\n", 
                        f_pdb, natoms );

    read_p4a_file(fp_p4a);

    nmaps=read_maplib(fp_map);
    if (report) 
        printf("Number of mapping functions read in MAP file %s = %d\n",
                        f_map, nmaps );

    /* Use the cifparse framework to extract information from an cif file: */
    ndb_cif_init();
    stat = ndb_cif_read_file(fp_nef);
    fclose( fp_nef );

#define CIFP cifFiles.datablocks[iblock].categories[iCat]

    /*  ----  First process the distance restraint lists:   ------ */

    for (iblock = 0; iblock < cifFiles.numDatablock; iblock++) {

        iCat = get_category_index(iblock, "nef_distance_restraint_list");
        if (iCat == -1) {
            continue;
        }

        col1 = get_column_index(iblock, iCat, "sf_category"); 
        assert( col1 >= 0 );
        assert( !strcmp( CIFP.rows[0].columns[col1],
          "nef_distance_restraint_list" ));
       
        iCat = get_category_index(iblock, "nef_distance_restraint");
        if (iCat == -1) {
            fprintf(stderr, "Target category %s is not in data block %d\n",
                    "nef_distance_restraint", iblock);
            exit(1);
        }

        fprintf( fp_rst, "#  converting NEF block %s\n# \n",
            cifFiles.datablocks[iblock].datablockName );

        col1 = get_column_index(iblock, iCat, "restraint_id"); /* id */
        assert( col1 >= 0 );
        col2 = get_column_index(iblock, iCat, "sequence_code_1");
        assert( col2 >= 0 );
        col3 = get_column_index(iblock, iCat, "residue_name_1");
        assert( col3 >= 0 );
        col4 = get_column_index(iblock, iCat, "atom_name_1");
        assert( col4 >= 0 );
        col6 = get_column_index(iblock, iCat, "sequence_code_2");
        assert( col6 >= 0 );
        col7 = get_column_index(iblock, iCat, "residue_name_2");
        assert( col7 >= 0 );
        col8 = get_column_index(iblock, iCat, "atom_name_2");
        assert( col8 >= 0 );
        col9 = get_column_index(iblock, iCat, "upper_limit");
        assert( col9 >= 0 );
        col10 = get_column_index(iblock, iCat, "chain_code_1");
        assert( col10 >= 0 );
        col11 = get_column_index(iblock, iCat, "chain_code_2");
        assert( col11 >= 0 );

#define R_ID atoi(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col1])
#define RESNO_1 atoi(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col2])
#define RESNO_2 atoi(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col6])
#define CID_1 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col10][0]
#define CID_2 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col11][0]
#define RESNAME_1 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col3]
#define RESNAME_2 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col7]
#define ATNAME_1 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col4]
#define ATNAME_2 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col8]
#define UPPERLIM atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col9])

        /*  --- need some better error processing here!  */

        is_ambig = 0;

        for (i = 0; i < CIFP.numRow; i++) {

            fprintf( fp_rst, "#  %d    %c %d %s %s     %c %d %s %s   %10.3f\n", 
                 R_ID, CID_1, RESNO_1, RESNAME_1, ATNAME_1,
                       CID_2, RESNO_2, RESNAME_2, ATNAME_2,  UPPERLIM );

            /* convert chainId/resno into Amber residue sequence number: */

            amber_rnum1 = residue_renumber( CID_1, RESNO_1 );
            amber_rnum2 = residue_renumber( CID_2, RESNO_2 );

            /* look ahead to see this is an ambiguous restraint: */

            curr_id = R_ID;
            if( i == CIFP.numRow - 1 ) next_id = -1;
            else next_id = atoi(CIFP.rows[i+1].columns[col1]);

            if( curr_id == next_id ){  /* need to store info for later */
                is_ambig++;
                if( is_ambig == 1 ){  /* this is a new restraint */
                    asprintf( &ambig1, "AMBIG %d:%d:1 = ", iblock, curr_id );
                    asprintf( &ambig2, "AMBIG %d:%d:2 = ", iblock, curr_id );
                }
                asprintf( &ambig1, "%s %s %d ", ambig1, 
                          ATNAME_1, amber_rnum1 );
                asprintf( &ambig2, "%s %s %d ", ambig2, 
                          ATNAME_2, amber_rnum2 );

            } else {

                if( is_ambig ){   /* we are at the end of an ambig restraint */
                    asprintf( &ambig1, "%s %s %d ", ambig1, 
                              ATNAME_1, amber_rnum1 );
                    asprintf( &ambig2, "%s %s %d ", ambig2, 
                              ATNAME_2, amber_rnum2 );
                    parse_ambig( ambig1 );
                    parse_ambig( ambig2 );
                    sprintf( data, "0 AMB %d:%d:1  0 AMB %d:%d:2  %8.3f",
                        iblock, curr_id, iblock, curr_id, UPPERLIM );

                } else {   /* this is an unamabigous restraint */

                    sprintf( data, "%d %s %s     %d %s %s   %10.3f\n", 
                          amber_rnum1, RESNAME_1, ATNAME_1,
                          amber_rnum2, RESNAME_2, ATNAME_2, UPPERLIM );
                }
                is_ambig = 0; 
                initialize_data( data, comment );
                parse_upb(data,comment);
                translate_atoms(data);
                correct_upperbound(data);
                l_bound = 1.8;
                write_rst(fp_rst,data,comment);
            }
        }
    }

    /*  ----  Now process the dihedral restraint lists:   ------ */

    for (iblock = 0; iblock < cifFiles.numDatablock; iblock++) {

        iCat = get_category_index(iblock, "nef_dihedral_restraint_list");
        if (iCat == -1) {
            continue;
        }

        col1 = get_column_index(iblock, iCat, "sf_category"); 
        assert( col1 >= 0 );
        assert( !strcmp( CIFP.rows[0].columns[col1],
          "nef_dihedral_restraint_list" ));
       
        iCat = get_category_index(iblock, "nef_dihedral_restraint");
        if (iCat == -1) {
            fprintf(stderr, "Target category %s is not in data block %d\n",
                    "nef_dihedral_restraint", iblock);
            exit(1);
        }

        fprintf( fp_rst, "#  converting NEF block %s\n# \n",
            cifFiles.datablocks[iblock].datablockName );

        col1 = get_column_index(iblock, iCat, "restraint_id"); /* id */
        assert( col1 >= 0 );
        col2 = get_column_index(iblock, iCat, "sequence_code_1");
        assert( col2 >= 0 );
        col3 = get_column_index(iblock, iCat, "atom_name_1");
        assert( col3 >= 0 );
        col4 = get_column_index(iblock, iCat, "sequence_code_2");
        assert( col4 >= 0 );
        col5 = get_column_index(iblock, iCat, "atom_name_2");
        assert( col5 >= 0 );
        col6 = get_column_index(iblock, iCat, "sequence_code_3");
        assert( col6 >= 0 );
        col7 = get_column_index(iblock, iCat, "atom_name_3");
        assert( col7 >= 0 );
        col8 = get_column_index(iblock, iCat, "sequence_code_4");
        assert( col8 >= 0 );
        col9 = get_column_index(iblock, iCat, "atom_name_4");
        assert( col9 >= 0 );
        col10 = get_column_index(iblock, iCat, "lower_limit");
        assert( col10 >= 0 );
        col11 = get_column_index(iblock, iCat, "upper_limit");
        assert( col11 >= 0 );
        col12 = get_column_index(iblock, iCat, "chain_code_1");
        assert( col12 >= 0 );
        col13 = get_column_index(iblock, iCat, "chain_code_2");
        assert( col13 >= 0 );
        col14 = get_column_index(iblock, iCat, "chain_code_3");
        assert( col14 >= 0 );
        col15 = get_column_index(iblock, iCat, "chain_code_4");
        assert( col15 >= 0 );
        col16 = get_column_index(iblock, iCat, "name");
        assert( col16 >= 0 );

#define DRESNO_1 atoi(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col2])
#define DRESNO_2 atoi(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col4])
#define DRESNO_3 atoi(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col6])
#define DRESNO_4 atoi(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col8])
#define DCID_1 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col12][0]
#define DCID_2 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col13][0]
#define DCID_3 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col14][0]
#define DCID_4 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col15][0]
#define DATNAME_1 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col3]
#define DATNAME_2 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col5]
#define DATNAME_3 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col7]
#define DATNAME_4 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col9]
#define DLOWERLIM atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col10])
#define DUPPERLIM atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col11])
#define DNAME cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col16]

        for (i = 0; i < CIFP.numRow; i++) {

            fprintf( fp_rst, "#  %d   %c %d %s | %c %d %s | %c %d %s | %c %d %s | %10.3f  %10.3f  %s\n", 
             R_ID, DCID_1, DRESNO_1, DATNAME_1,
                   DCID_2, DRESNO_2, DATNAME_2,
                   DCID_3, DRESNO_3, DATNAME_3,
                   DCID_4, DRESNO_4, DATNAME_4, DLOWERLIM, DUPPERLIM, DNAME );

            /* convert chainId/resno into Amber residue sequence number: */

            amber_rnum1 = residue_renumber( DCID_1, DRESNO_1 );
            amber_rnum2 = residue_renumber( DCID_2, DRESNO_2 );
            amber_rnum3 = residue_renumber( DCID_3, DRESNO_3 );
            amber_rnum4 = residue_renumber( DCID_4, DRESNO_4 );

            r2 = DLOWERLIM;
            r3 = DUPPERLIM;
            if( r2 > r3 ) r3 = r3 + 360.; /* try to ensure that r2<r3  */
            r1 = r2 - R1_R4_RANGE;  r4 = r3 + R1_R4_RANGE;

            for( j=0; j<4; j++ ) anum[j] = 0;
            for( j=0; j<natoms; j++ ){
                if( mol[j].rnum == amber_rnum1  && 
                    !strcmp( mol[j].atyp, DATNAME_1 ) )
                    anum[0] = j+1;
                if( mol[j].rnum == amber_rnum2  && 
                    !strcmp( mol[j].atyp, DATNAME_2 ) )
                    anum[1] = j+1;
                if( mol[j].rnum == amber_rnum3  && 
                    !strcmp( mol[j].atyp, DATNAME_3 ) )
                    anum[2] = j+1;
                if( mol[j].rnum == amber_rnum4  && 
                    !strcmp( mol[j].atyp, DATNAME_4 ) )
                    anum[3] = j+1;
            }
            assert( anum[0] > 0 );
            assert( anum[1] > 0 );
            assert( anum[2] > 0 );
            assert( anum[3] > 0 );


            fprintf( fp_rst, " &rst     iat = %5d, %5d, %5d, %5d,\n",
                anum[0], anum[1], anum[2], anum[3]);
            fprintf( fp_rst, "    r1 = %5.1f, r2 = %5.1f, r3 = %5.1f, r4 = %5.1f,\n", 
                r1, r2, r3, r4);
            if( i == 0 ){
                fprintf( fp_rst, "    rk2 = %5.1f, ", K_FORCE_LB);
                fprintf( fp_rst, "rk3 = %5.1f, /\n\n", K_FORCE_UB);
            }else{
                fprintf(fp_rst, " /\n\n");
            }
        }
    }


    /*  ----  Now process the rdc restraint lists:   ------ */

    dataset = 0; tot = 0;
    fprintf( fp_rdc, " &align\n" );

    for (iblock = 0; iblock < cifFiles.numDatablock; iblock++) {

        iCat = get_category_index(iblock, "nef_rdc_restraint_list");
        if (iCat == -1) {
            continue;
        }

        col1 = get_column_index(iblock, iCat, "sf_category"); 
        assert( col1 >= 0 );
        assert( !strcmp( CIFP.rows[0].columns[col1],
          "nef_rdc_restraint_list" ));
       
        iCat = get_category_index(iblock, "nef_rdc_restraint");
        if (iCat == -1) {
            fprintf(stderr, "Target category %s is not in data block %d\n",
                    "nef_rdc_restraint", iblock);
            exit(1);
        }

        fprintf( fp_rdc, "! \n!  converting NEF block %s\n! \n",
            cifFiles.datablocks[iblock].datablockName );
        dataset++;

        col1 = get_column_index(iblock, iCat, "restraint_id"); /* id */
        assert( col1 >= 0 );
        col2 = get_column_index(iblock, iCat, "sequence_code_1");
        assert( col2 >= 0 );
        col3 = get_column_index(iblock, iCat, "atom_name_1");
        assert( col3 >= 0 );
        col4 = get_column_index(iblock, iCat, "sequence_code_2");
        assert( col4 >= 0 );
        col5 = get_column_index(iblock, iCat, "atom_name_2");
        assert( col5 >= 0 );
        col6 = get_column_index(iblock, iCat, "target_value");
        assert( col6 >= 0 );
        col7 = get_column_index(iblock, iCat, "chain_code_1");
        assert( col7 >= 0 );
        col8 = get_column_index(iblock, iCat, "chain_code_2");
        assert( col8 >= 0 );

#define ACID_1 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col7][0]
#define ACID_2 cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col8][0]
#define ATARGET atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col6])

        for (i = 0; i < CIFP.numRow; i++) {

            tot++;
            fprintf( fp_rdc, "! \n!  %d    %c %d %s | %c %d %s  %10.3f\n", 
                 R_ID, ACID_1, DRESNO_1, DATNAME_1,
                       ACID_2, DRESNO_2, DATNAME_2, ATARGET );
            dobs = ATARGET;

            /* convert chainId/resno into Amber residue sequence number: */

            amber_rnum1 = residue_renumber( ACID_1, DRESNO_1 );
            amber_rnum2 = residue_renumber( ACID_2, DRESNO_2 );

            for( j=0; j<2; j++ ) anum[j] = 0;
            for( j=0; j<natoms; j++ ){
                if( mol[j].rnum == amber_rnum1  && 
                    !strcmp( mol[j].atyp, DATNAME_1 ) ){
                    anum[0] = j+1;
                    break;
                }
            }
            for( j=0; j<natoms; j++ ){
                if( mol[j].rnum == amber_rnum2  && 
                    !strcmp( mol[j].atyp, DATNAME_2 ) ){
                    anum[1] = j+1;
                    break;
                }
            }
            assert( anum[0] > 0 );
            assert( anum[1] > 0 );

            if ( !strcmp(DATNAME_1,"CA") ){  /* CA-HA */
                gigj = 1.4048;
                dij = 1.12;
                fprintf( fp_rst, "#  C angles for rdc restraint %d\n", tot );
                fprintf( fp_rst, " &rst   iat= %d, %d, %d, 0, r1=100., r2=109.5, r3=109.5, r4=120., rk2=100., rk3=100. /\n",
                    anum[1], anum[1]-1, anum[1]+1 ); /* HA-CA-CB */
                fprintf( fp_rst, " &rst   iat= %d, %d, %d, 0, r1=100., r2=109.5, r3=109.5, r4=120.,  /\n",
                    anum[1], anum[1]-1, anum[1]-3 ); /* HA-CA-N */
                /* fprintf( fp_rst, " &rst   iat= %d, %d, %d, 0, r1=100., r2=109.5, r3=109.5, r4=120.,  /\n",
                    anum[1]+1, anum[1], anum[1]-4 ); ??? */

            } else if (!strcmp(DATNAME_1,"N") ) {  /* N-H */
                gigj = -0.5663;
                dij = 1.04;
                fprintf( fp_rst, "#  N angles for rdc restraint %d\n", tot );
                fprintf( fp_rst, " &rst   iat= %d, %d, %d, 0, r1=110., r2=120., r3=120., r4=130., rk2=100., rk3=100. /\n",
                    anum[1], anum[1]-1, anum[1]-3 );  /* H-N-C' */
                fprintf( fp_rst, " &rst   iat= %d, %d, %d, 0, r1=110., r2=120., r3=120., r4=130.,  /\n",
                    anum[1], anum[1]-1, anum[1]+1 ); /* H-N-CA */

            } else {
                fprintf( stderr, "Only N or CA rdc restraints current supported\n" );
                exit(1);
            }

            fprintf( fp_rdc, "    id(%d)=%d,   jd(%d)=%d,  dobsl(%d)=%5.2f, dobsu(%d)=%5.2f,\n    gigj(%d)=%7.4f, dij(%d)=%4.2f, dataset(%d)=%d,\n",
                tot, anum[0], tot, anum[1], tot, dobs, tot, dobs, 
                tot, gigj, tot, dij, tot, dataset );
        }
    }
    fprintf( fp_rdc,  "   ndip=%d,  num_datasets=%d, dwt=%d*0.5, \n /\n\n",
        tot, dataset, tot );


    /* should free up memory cifparse has used; (but program is ending) */

    if( print_rst ) fclose(fp_rst);
    if( print_rdc ) fclose(fp_rdc);

    return(0);
}
