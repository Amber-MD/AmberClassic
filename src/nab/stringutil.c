/*
 *	Set of routines used to execute nab string ops of =, +, substring
 *	and allocate a string for dest in sprinf() and *scanf()
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "nab.h"

char	*substr( char [], int, int );
char	*NAB_getline( FILE * );
int	split( char [], char *[], char * );
int	split_n( char [], int, char *[], char * );
int	NAB_index( char [], char [] );
char	*NAB_strcat( char [], char [] );
char	*NAB_strcpy( char **, char [] );
int	NAB_strcmp( char *, char * );
char	*NAB_readstring( char ** );
int	NAB_newstring( char **, int * );

char	NAB_rsbuf[ NAB_RSBUF_SIZE ];

char	*substr( char str[], int pos, int len )
{
	int	slen, alen;
	char	*ssp;

	slen = strlen( str );
	if( pos < 1 || pos > slen ){
		/* fprintf( stderr, "substr: pos (%d) not in string\n", pos ); */
		return( NULL );
	}else{
		alen = slen - pos + 1;
		if( len > alen ) len = alen;  /* just get as much as you can */
		if(( ssp = (char *)malloc( (len + 1) * sizeof(char))) == NULL){
			fprintf( stderr, "substr: can't allocate substring\n" );
			return( NULL );
		}
		strncpy( ssp, &str[ pos - 1 ], len );
		ssp[ len ] = '\0';
		return( ssp );
	}
}

char	*NAB_getline( FILE *fp )
{
	static	char	buf[ 1024 ];
	char	*bp;

	if( fgets( buf, sizeof( buf ), fp ) ){
		for( bp = buf; *bp; bp++ ){
			if( *bp == '\n' ){
				*bp = '\0';
				break;
			}
		}
		return( buf );
	}else
		return( NULL );
}

int	split( char str[], char *fields[], char *fsep )
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
			if( ( efp = strpbrk( fp, fsep ) ) ){
				if( !( flen = efp - fp ) )
					return( nf );
				nfp = (char *)malloc((flen + 1) * sizeof(char));
				strncpy( nfp, fp, flen );
				nfp[ flen ] = '\0';
				fields[ nf ] = nfp;
				sp = efp;
				nf++;
			}else{
				flen = strlen( fp );
				nfp = (char *)malloc((flen + 1) * sizeof(char));
				strcpy( nfp, fp );
				fields[ nf ] = nfp;
				nf++;
				return( nf );
			}
		}
	}else{
		for( nf = 0, sp = str; ; ){
			if( ( fp = strchr( sp, *fsep ) ) ){
				flen = fp - sp;
				nfp = (char *)malloc((flen + 1) * sizeof(char));
				strncpy( nfp, sp, flen );
				nfp[ flen ] = '\0';
				fields[ nf ] = nfp;
				nf++;
				sp = fp + 1;
			}else{
				flen = strlen( sp );
				nfp = (char *)malloc((flen + 1) * sizeof(char));
				strcpy( nfp, sp );
				fields[ nf ] = nfp;
				nf++;
				return( nf );
			}
		}
	}
}

int	split_n( char str[], int s_fields, char *fields[], char *fsep )
{
	int	nf, flen, white;
	char	*sp, *fp, *efp, *nfp;

	if( !str )
		return( 0 );

	/* Use fsep of white (tab,space,newline) space is special	*/ 
	/* although a \n is a white space character, an fsep of a	*/
	/* single \n is not considered white space, allowing multi line	*/
	/* strings to be broken into lines				*/
	/* Likewise an fsep of a \t is not white space allowing tab	*/
	/* separated files to contain NULL fields			*/

	if( !strcmp( fsep, "\n" ) )
		white = 0;
	else if( !strcmp( fsep, "\t" ) )
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
			if( ( efp = strpbrk( fp, fsep ) ) ){
				if( !( flen = efp - fp ) )
					return( nf );
				if( nf < s_fields ){
					nfp = (char *)
						malloc((flen+1)*sizeof(char));
					strncpy( nfp, fp, flen );
					nfp[ flen ] = '\0';
					fields[ nf ] = nfp;
				}
				sp = efp;
				nf++;
			}else{
				if( nf < s_fields ){
					flen = strlen( fp );
					nfp = (char *)
						malloc((flen+1)*sizeof(char));
					strcpy( nfp, fp );
					fields[ nf ] = nfp;
				}
				nf++;
				return( nf );
			}
		}
	}else{
		for( nf = 0, sp = str; ; ){
			if( ( fp = strchr( sp, *fsep ) ) ){
				if( nf < s_fields ){
					flen = fp - sp;
					nfp = (char *)
						malloc((flen+1)*sizeof(char));
					strncpy( nfp, sp, flen );
					nfp[ flen ] = '\0';
					fields[ nf ] = nfp;
				}
				nf++;
				sp = fp + 1;
			}else{
				if( nf < s_fields ){
					flen = strlen( sp );
					nfp = (char *)
						malloc((flen+1)*sizeof(char));
					strcpy( nfp, sp );
					fields[ nf ] = nfp;
				}
				nf++;
				return( nf );
			}
		}
	}
}

int	NAB_index( char s[], char t[] )
{
	char	*ip;

	if( !s || !t )
		return( 0 );
	if( !( ip = strchr( s, *t ) ) )
		return( 0 );
	else
		return( ip - s + 1 );
}

char	*NAB_strcat( char s[], char t[] )
{
	int	sl, tl;
	char	*np0, *np, *sp, *tp;

	if( !s && !t )
		return( NULL );

	sl = s ? strlen( s ) : 0;
	tl = t ? strlen( t ) : 0;
	if( !( np0 = ( char * )malloc( (sl + tl + 1) * sizeof(char) ) ) ){
		fprintf( stderr, "NAB_strcat: can't allocate new string\n" );
		exit( 1 );
	}

	if( s ){
		for( np = np0, sp = s; *sp; sp++ )
			*np++ = *sp;
	}else
		np = np0;
	if( t ){
		for( tp = t; ( *np++ = *tp++ ); )
			;
	}else
		*np0 = '\0';

	return( np0 );
}

char	*NAB_strcpy( char **s, char t[] )
{
	int	tl;

	if( !t ){
		if( !*s )
			free( *s );
		*s = NULL;
		return( NULL );
	}

	tl = strlen( t );

/* Memory leak fix via Andreas and Tom. */
/* *s is guaranteed to be either NULL or a previously returned pointer */
/* from malloc, et alia. */

	if( !(*s=(char*)realloc(*s,(tl+1)*sizeof(char)))){
		fprintf( stderr, "NAB_strcpy: can't allocate new string\n" );
		exit( 1 );	
	}
	strcpy(*s,t);
	return(*s);
}

int	NAB_strcmp( char *s, char *t )
{

	if( s && t )
		return( strcmp( s, t ) );
	else
		return( s - t );
}

char	*NAB_readstring( char **sp )
{
/* Caution:  hard coded maximum size of temporary strings !! */

#define MAXSTRINGLENGTH 1069
/* The first four-digit quasiall-even-digits non-quasi-repdigit emirp ! */

	if( ( *sp = ( char * )malloc( MAXSTRINGLENGTH * sizeof(char) ) )
			== NULL ){
		fprintf( stderr,
			"NAB_readstring: can't allocate new string\n" );
		exit( 1 );
	}
	**sp = '\0';

	return( *sp );
}

int	NAB_newstring( char **sp, int *size )
{

	return( ( *sp = ( char * )malloc( (*size) * sizeof(char) ) ) ? 0 : 1 );
}
