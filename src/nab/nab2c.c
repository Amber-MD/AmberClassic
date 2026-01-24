#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "nab.h"
#include "cgen.h"

extern	FILE	*yyin;
extern	int	yydebug;
extern int yyparse();

int	cg_dval = CGD_NONE;
int	cg_lineno = 1;
int	cg_emsg_lineno = 1;
char	cg_nfname[ 256 ] = "";
int	cg_noassert = 0;
int	cg_nodebug = 0;

static	int	setdebugval( char * );
static	int	init( char nfname[], int aopt );
static	int	split( char [], char *[], char * );
int	yyerror( void );

int main( int argc, char *argv[] )
{
	int	ac;
	int	sym;
	int	aopt = 0;
	char	nfname[ 256 ];

	*nfname = '\0';
	for( ac = 1; ac < argc; ac++ ){
		if( strcmp( argv[ ac ], "-nfname" ) == 0 ){
			ac++;
			strcpy( nfname, argv[ ac ] );
		}else if( strcmp( argv[ ac ], "-avs" ) == 0 )
			aopt = 1;
		else if( strncmp( argv[ ac ], "-noassert", 9 ) == 0 )
			cg_noassert = 1;
		else if( strncmp( argv[ ac ], "-cgdebug", 8 ) == 0 )
			cg_dval = setdebugval( argv[ ac ] );
		else if( strncmp( argv[ ac ], "-nodebug", 8 ) == 0 )
			cg_nodebug = 1;
	}

	if( cg_dval & CGD_PARSER )
		yydebug = 1;

	if( init( nfname, aopt ) )
		exit( 1 );

	while( (sym = yyparse()) ){
		fprintf( stderr, "%s:%d syntax error\n", 
			nfname, cg_lineno );
		CG_exit( 1 );
	}

	exit( 0 );
}

static	int	setdebugval( char *dval )
{
	char	*dp;
	int	f, nf, cgd;
	char	*dvals[ 10 ];

	if( (dp = strchr( dval, '=' )) ){
		dp++;
		nf = split( dp, dvals, "," );
		cgd = CGD_NONE;
		for( f = 0; f < nf; f++ ){
			if( !strcmp( dvals[ f ], "VAR" ) )
				cgd |= CGD_VARDECL;
			else if( !strcmp( dvals[ f ], "FIX" ) )
				cgd |= CGD_FIXEXPR;
			else if( !strcmp( dvals[ f ], "CHK" ) )
				cgd |= CGD_CHKEXPR;
			else if( !strcmp( dvals[ f ], "GEN" ) )
				cgd |= CGD_EXPRCODE;
			else if( !strcmp( dvals[ f ], "PAR" ) )
				cgd |= CGD_PARSER;
		}
		return( cgd );
	}else
		return( CGD_ALL );
}

static	int	init( char nfname[], int aopt )
{
	char	cfname[ 256 ];
	char	*np, *cp, *dp;

	if( !*nfname ){
		fprintf( stderr, "nab2c: Internal error: nfname required\n" );
		exit( 1 );
	}

	strcpy( cg_nfname, nfname );

	for( np = nfname, dp = NULL; *np; np++ ){
		if( *np == '.' )
			dp = np;
	}

	for( cp = cfname, np = nfname; np <= dp; np++ )
		*cp++ = *np;
	strcpy( cp, "c" );

	return( CG_init( cfname, aopt ) );
}

static	int	split( char str[], char *fields[], char *fsep )
{
	int	nf, flen, white;
	char	*sp, *fp, *efp, *nfp;

	if( !str )
		return( 0 );

	/* Use fsep of white space is special				*/ 
	/* although a \n is a white space character, a fsep of a	*/
	/* single \n is not considered white space, allowing multi line	*/
	/* strings to be broken into lines				*/

	if( !strcmp( fsep, "\n" ) )
		white = 0;
	else if( strspn( fsep, " \t\n" ) == strlen( fsep ) ){
		white = 1;
	}else
		white = 0;

	if( white ){
		for( nf = 0, sp = str; ; ){
			fp = sp + strspn( sp, fsep );
			if( !*fp )
				return( nf );
			if( (efp = strpbrk( fp, fsep )) ){
				if( !( flen = efp - fp ) )
					return( nf );
				nfp = (char *)malloc( (flen + 1) * sizeof( char ) );
				strncpy( nfp, fp, flen );
				nfp[ flen ] = '\0';
				fields[ nf ] = nfp;
				sp = efp;
				nf++;
			}else{
				flen = strlen( fp );
				nfp = (char *)malloc( (flen + 1) * sizeof( char ) );
				strcpy( nfp, fp );
				fields[ nf ] = nfp;
				nf++;
				return( nf );
			}
		}
	}else{
		for( nf = 0, sp = str; ; ){
			if( (fp = strchr( sp, *fsep )) ){
				flen = fp - sp;
				nfp = (char *)malloc( (flen + 1) * sizeof( char ) );
				strncpy( nfp, sp, flen );
				nfp[ flen ] = '\0';
				fields[ nf ] = nfp;
				nf++;
				sp = fp + 1;
			}else{
				flen = strlen( sp );
				nfp = (char *)malloc( (flen + 1) * sizeof( char ) );
				strcpy( nfp, sp );
				fields[ nf ] = nfp;
				nf++;
				return( nf );
			}
		}
	}
}

int	yyerror( void )
{

	return( 0 );
}
