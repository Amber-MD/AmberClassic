#include <stdio.h>
#include <string.h>
#include "nab.h"

char	*getenv();

#define	D_NAME_SIZE	20
typedef	struct	distrec	{
	char	d_name[ D_NAME_SIZE ];
	int	d_count;
	REAL_T	d_mean;
	REAL_T	d_stddev;
	REAL_T	d_min;
	REAL_T	d_max;
} DISTREC;

#define	MAXDRECS	30000
static	DISTREC	drecs[ MAXDRECS ];
static	int	n_drecs;
static  char	savename[256];
static  char	filename[256];
static  char	nabhome[256];

int	getdist( char **name, char **lib, int *dcount,
	REAL_T *dmean, REAL_T *dstddev, REAL_T *dmin, REAL_T *dmax )

{
    char    *amberhome;
	FILE	*fp;
	char	line[ 100 ];
	DISTREC	*dp;
	int	i, j, k, cv;

    if( !( amberhome = (char *) getenv( "AMBERHOME" ) ) ){
         fprintf( stderr, "AMBERHOME not defined.\n" );
         return( 1 );
    }
	if ( ( strcmp ( *lib, savename ) != 0 ) ){
		n_drecs = 0;
		strcpy( savename, *lib );
		if( ( fp = fopen( *lib, "r" ) ) == NULL  ){
			strcpy( filename, *lib );
			sprintf( nabhome, "%s/dat/dgdb/%s", amberhome, *lib );
			if( ( fp = fopen( nabhome, "r" ) ) == NULL  ){
				fprintf( stderr, 
			"getdist: could not find database %s\n", *lib );
						return( 1 );
			}
		}
		dp = drecs;
		while( fgets( line, sizeof( line ), fp ) ){
			if( *line == '#' )
				continue;
#ifdef NAB_DOUBLE_PRECISION
			sscanf( line, "%s %d %lf %lf %lf %lf",
#else
			sscanf( line, "%s %d %f %f %f %f",
#endif
				dp->d_name, &dp->d_count,
				&dp->d_mean, &dp->d_stddev,
				&dp->d_min, &dp->d_max );
			dp++;
			n_drecs++;
		}
		fclose( fp );
	}

	for( i = 0, j = n_drecs - 1; i <= j; ){
		k = ( i + j ) / 2;
		dp = &drecs[ k ];
		if( ( cv = strcmp( dp->d_name, *name ) ) == 0 ){
			*dcount = dp->d_count;
			*dmean = dp->d_mean;
			*dstddev = dp->d_stddev;
			*dmin = dp->d_min;
			*dmax = dp->d_max;
			return( 0 );
		}else if( cv < 0 )
			i = k + 1;
		else
			j = k - 1;
	}
	*dcount = 0;
	*dmean = 0.0;
	*dstddev = 0.0;
	*dmin = 0.0;
	*dmax = 0.0;
	return( 1 );
}
