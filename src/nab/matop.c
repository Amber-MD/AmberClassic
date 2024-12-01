#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "nabcode.h"
#include "matop.h"

#define	BPN	14

#define	GSBUFSIZE	10000
static	char	gsbuf[ GSBUFSIZE ];
static	char	*gsb;


int	MAT_fprint( FILE *fp, int nmats, MATRIX_T mats[] )
{
	int	m, i;

	for( m = 0; m < nmats; m++ ){
		for( i = 0; i < 4; i++ ){
			fprintf( fp, "%13.6le %13.6le %13.6le %13.6le\n",
				mats[m][i][0], mats[m][i][1],
				mats[m][i][2], mats[m][i][3] );
		}
	}

	return( 0 );
}

int	MAT_sprint( char str[], int nmats, MATRIX_T mats[] )
{
	int	m, i;
	char	*sp;

	for( sp = str, m = 0; m < nmats; m++ ){
		for( i = 0; i < 4; i++ ){
			sprintf( sp, "%13.6le %13.6le %13.6le %13.6le",
				mats[m][i][0], mats[m][i][1],
				mats[m][i][2], mats[m][i][3] );
			sp += 4 * BPN - 1;
			*sp++ = ( i == 3 && m == nmats - 1 ) ? '\0' : '\n';
		}
	}

	return( 0 );
}

int	MAT_fscan( FILE *fp, int smats, MATRIX_T mats[] )
{
	char	*lp, line[ 1024 ];
	int	c, bl;
	int	r, nmats, nnum;

	gsb = gsbuf;
	for( nnum = 0, nmats = 0, r = 0; ( c = getc( fp ) ) != EOF; ){
		ungetc( c, fp );
		fgets( line, sizeof( line ), fp );
		if( *line == '#' ){
			strcpy( gsb, line );
			gsb += strlen( line );
			continue;
		}
		for( bl = 1, lp = line; *lp; lp++ ){
			if( !isspace( *lp ) ){
				bl = 0;
				break;
			}
		}
		if( bl )
			continue;
		if( nmats >= smats ){
			fprintf( stderr, "MAT_fscan: too many matrices.\n" );
			return( -1 );
		}
		sscanf( line, "%lf %lf %lf %lf",
			&mats[nmats][r][0], &mats[nmats][r][1],
			&mats[nmats][r][2], &mats[nmats][r][3] );
		nnum += 4;
		r++;
		if( r == 4 ){
			r = 0;
			nmats++;
		}
	}

	if( nnum % 16 ){
		fprintf( stderr, "MAT_fscan: incomplete matrix.\n" );
		return( -1 );
	}

	return( nmats );
}

int	MAT_sscan( char str[], int smats, MATRIX_T mats[] )
{
	char	*lp, *elp, *nlp;
	int	l, llen, bl;
	int	r, nmats, nnum;

	gsb = gsbuf;
	for( nnum = 0, nmats = 0, r = 0, lp = str; lp && *lp; lp = nlp ){
		if( (elp = strchr( lp, '\n' )) ){
			llen = elp - lp;
			nlp = elp + 1;
		}else{
			llen = strlen( lp );
			nlp = NULL;
		}
		if( *lp == '#' ){
			strncpy( gsb, lp, llen + 1 ); /* llen up to \n */
			gsb[ llen + 1 ] = '\0';
			gsb += llen + 1;
			continue;
		}
		for( bl = 1, l = 0; l < llen; l++ ){
			if( !isspace( lp[ l ] ) ){
				bl = 0;
				break;
			}
		}
		if( bl )
			continue;
		if( nmats >= smats ){
			fprintf( stderr, "MAT_sscan: too many matrices.\n" );
			return( -1 );
		}
		sscanf( lp, "%lf %lf %lf %lf",
			&mats[nmats][r][0], &mats[nmats][r][1], 
			&mats[nmats][r][2], &mats[nmats][r][3] );
		nnum += 4;
		r++;
		if( r == 4 ){
			r = 0;
			nmats++;
		}
	}

	if( nnum % 16 ){
		fprintf( stderr, "MAT_sscan: incomplete matrix.\n" );
		return( -1 );
	}

	return( nmats );
}

REF_MATRIX_T	MAT_concat( MATRIX_T ml, MATRIX_T mr )
{
	int	i, j, k;
	static	MATRIX_T m;
	
	for( i = 0; i < 4; i ++ ){
		for( j = 0; j < 4; j++ ){
			m[i][j] = 0.0;
			for( k = 0; k < 4; k++ )
				m[i][j] += ml[i][k] * mr[k][j];
		}
	}
	return( m );
}

int	MAT_count( char buf[] )
{
	char	*lp, *elp, *nlp;
	int	nnum;
	int	l, inspace, llen;

	for( nnum = 0, lp = buf; lp && *lp; lp = nlp ){
		if( (elp = strchr( lp, '\n' )) ){
			llen = elp - lp;
			nlp = elp + 1;
		}else{
			llen = strlen( lp );
			nlp = NULL;
		}
		if( *lp == '#' )
			continue;
		for( inspace = 1, l = 0; l < llen; l++ ){
			if( isspace( lp[ l ] ) ){
				inspace = 1;
			}else if( inspace ){
				inspace = 0;
				nnum++;
			}
		}
	}

	if( nnum % 16 ){
		fprintf( stderr, "MAT_count: last marix is incomplete.\n" );
		return( -1 );
	}else
		return( nnum / 16 );
}

char	*MAT_getsyminfo( void )
{
	char	*bp, *lp, *elp, *nlp;
	int	llen, size;
	char	*sbuf, *sbp;
		
	bp = gsbuf;
	for( size = 0, lp = bp; lp && *lp; lp = nlp ){
		if( (elp = strchr( lp, '\n' )) ){
			llen = elp - lp + 1;
			nlp = elp + 1;
		}else{
			llen = strlen( lp ) + 1;
			nlp = NULL;
		}
		if( *lp != '#' )
			continue;
		size += llen;
	}

	if( !( sbuf = ( char * )malloc( ( size + 1 ) * sizeof( char ) ) ) ){
		fprintf( stderr,
	"MAT_getsyminfo: can't allocate space for symmetry information.\n" );
		return( NULL );
	}

	*sbuf = '\0';
	for( sbp = sbuf, lp = bp; lp && *lp; lp = nlp ){
		if( (elp = strchr( lp, '\n' )) ){
			llen = elp - lp + 1;
			nlp = elp + 1;
		}else{
			llen = strlen( lp ) + 1;
			nlp = NULL;
		}
		if( *lp != '#' )
			continue;
		if( nlp ){
			strncpy( sbp, lp, llen );
			sbp[ llen ] = '\0';
			sbp += llen;
		}else{
			strcpy( sbp, lp );
			sbp[ llen ] = '\n';
			sbp[ llen + 1 ] = '\0';
		}
	}

	return( sbuf );
}

int	MAT_istrue( MATRIX_T m )
{
	int	i, j;

	for( i = 0; i < 4; i++ ){
		for( j = 0; j < 4; j++ )
			if( m[i][j] != 0.0 )
				return( 1 );
	}
	return( 0 );
}
