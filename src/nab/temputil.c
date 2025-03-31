#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "nab.h"

int	*NAB_itemp( int );
REAL_T	*NAB_ftemp( REAL_T );
char	**NAB_ctemp( char * );
FILE	**NAB_fptemp( FILE * );

int	*NAB_itemp( int expr )
{
	int	*ip;

	if( ( ip = ( int * )malloc( sizeof( int ) ) ) == NULL ){
		fprintf( stderr, "NAB_itemp: can't allocate temporary\n" );
		exit( 1 );
	}
	*ip = expr;
	return( ip );
}

REAL_T	*NAB_ftemp( REAL_T expr )
{
	REAL_T	*fp;

	if( ( fp = ( REAL_T * )malloc( sizeof( REAL_T ) ) ) == NULL ){
		fprintf( stderr, "NAB_ftemp: can't allocate temporary\n" );
		exit( 1 );
	}
	*fp = expr;
	return( fp );
}

char	**NAB_ctemp( char *expr )
{
	char	**cp;
	int	len;

	if( ( cp = ( char ** )malloc( sizeof( char  * ) ) ) == NULL ){
		fprintf( stderr, "NAB_ftemp: can't allocate temporary 1\n" );
		exit( 1 );
	}
	len = strlen( expr ) + 1;
	if( ( *cp = ( char * )malloc( len * sizeof( char ) ) ) == NULL ){
		fprintf( stderr, "NAB_ftemp: can't allocate temporary 2\n" );
		exit( 1 );
	}
	strcpy( *cp, expr );
	return( cp );
}

FILE	**NAB_fptemp( FILE *expr )
{
	FILE	**fp;

	if( ( fp = ( FILE ** )malloc( sizeof( FILE * ) ) ) == NULL ){
		fprintf( stderr, "NAB_fptemp: can't allocate temporary\n" );
		exit( 1 );
	}
	*fp = expr;
	return( fp );
}
