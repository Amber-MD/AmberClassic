#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "nabhome.h"
#   include <sys/wait.h>

typedef	unsigned long	NAB_SIZE_T;

static	void	n2c( int, char *, int, int, int, char *[], char [], char * );
static	void	cc( int, int, int, char [], int, char *[], char [] );

int main( int argc, char *argv[] )
{
	int	ac;
	int	copt = 0;
	int	cgdopt = 0;
	char	*cgdval = NULL;
	int	noassert = 0;
	int	nodebug = 0;
	int	oopt = 0;
	char	cppstring[ 1024 ];
	char	*nfmask;
	char	*dotp;
	char	ofname[ 256 ];
	
	cppstring[ 0 ] = '\0';
	if( argc == 1 ){
		fprintf( stderr,
"usage: %s [-c] [-Dstring] [-noassert] [-nodebug] [-o file] [-v] file(s)\n",
			argv[ 0 ] ); 
		exit( 1 );
	}

	if( ( nfmask = ( char * )malloc( (argc + 1)*sizeof(int) ) ) == NULL ){
		fprintf( stderr, "%s: can't allocate arg mask\n", argv[ 0 ] );
		exit( 1 );
	}
	
	/* get nab options:	*/
	strcpy( ofname, "a.out" );
	for( ac = 1; ac < argc; ac++ ){
		nfmask[ ac ] = 0;
		if( strcmp( argv[ ac ], "-c" ) == 0 ){
			copt = 1;
		}else if( strncmp( argv[ ac ], "-cgdebug", 8 ) == 0 ){
			cgdopt = 1;
			cgdval = argv[ ac ]; 
			argv[ ac ] = NULL;
		}else if( strncmp( argv[ ac ], "-D", 2 ) == 0 ){
			strcat( cppstring, argv[ ac ] );
			strcat( cppstring, " ");
		}else if( strncmp( argv[ ac ], "-I", 2 ) == 0 ){
			strcat( cppstring, argv[ ac ] );
			strcat( cppstring, " ");
		}else if( strncmp( argv[ ac ], "-noassert", 9 ) == 0 ){
			noassert = 1;
			argv[ ac ] = NULL;
		}else if( strncmp( argv[ ac ], "-nodebug", 8 ) == 0 ){
			nodebug = 1;
			argv[ ac ] = NULL;
		}else if( strcmp( argv[ ac ], "-o" ) == 0 ){
			oopt = 1;
			ac++;
			if( ac == argc ){
				fprintf( stderr, "%s: -o requires file name\n",
					argv[ 0 ] );
				exit( 1 );
			}else{
				strcpy( ofname, argv[ ac ] );
				nfmask[ ac ] = 0;
			}
		}else if( strcmp( argv[ ac ], "-v" ) == 0 ){
			cgdopt = 1;
			cgdval = "";
			argv[ ac ] = NULL;
		}else if( *argv[ ac ] != '-' ){
			if( (dotp = strrchr( argv[ ac ], '.' )) ){
				if( strcmp( dotp, ".nab" ) == 0 )
					nfmask[ ac ] = 1;
			}
		}
	}

	n2c( cgdopt, cgdval, noassert, nodebug,
		argc, argv, nfmask, cppstring );

	cc( copt, cgdopt, oopt, ofname, argc, argv, nfmask );

	exit( 0 );
}

static	void	n2c( int cgdopt, char *cgdval, int noassert, int nodebug,
	int argc, char *argv[], char nfmask[], char *cppstring )
/*
int	cgdopt;
char	*cgdval;
int	noassert;
int	nodebug;
int	argc;
char	*argv[];
char	nfmask[];
char	*cppstring;
*/
{
	int	ac;
	char	*n2c_ofname;
	char	*cpp_ofname;

	char	*cmd;
	int	    status;
    int     cpp_ofd, n2c_ofd;

    asprintf( &n2c_ofname, "/tmp/n2c_ofname_XXXXXX" );
    if( ( n2c_ofd = mkstemp( n2c_ofname ) ) < 0 ){
              perror( n2c_ofname );
              exit(1);
    }

	for( ac = 1; ac < argc; ac++ ){
		if( nfmask[ ac ] ){
			if( access( argv[ ac ], F_OK ) ){
				fprintf( stderr,
					"%s: %s: no such file (arg # %d)\n",
					argv[ 0 ], argv[ ac ], ac );
				unlink( n2c_ofname );
				exit( 1 );
			}

           /*  run cpp :  */

           asprintf( &cpp_ofname, "/tmp/cpp_ofname_XXXXXX" );
           if( (cpp_ofd = mkstemp( cpp_ofname ) ) < 0 ){
                 perror( cpp_ofname );
                 exit(1);
           }

            char *nabhome = NABHOME;

			asprintf( &cmd, "%s/bin/%s %s -I%s/include %s > %s",
				nabhome, "ucpp -l", cppstring, nabhome,
				argv[ ac ] ? argv[ ac ] : "", cpp_ofname );
			if( cgdopt ) fprintf( stderr, "cpp cmd: %s\n", cmd );
            status = system( cmd );
			if( status != 0 ){
				unlink( cpp_ofname );
				fprintf( stderr, "cpp failed!\n" ); exit(1);
			}
            lseek( cpp_ofd, 0L, L_SET );   /* needed?? */

            /*  next run nab2c:  */

#ifdef MPI
			asprintf( &cmd, "%s/bin/mpinab2c %s %s %s -nfname %s < %s > %s",
#else
			asprintf( &cmd, "%s/bin/nab2c %s %s %s -nfname %s < %s > %s",
#endif
				nabhome,
				cgdopt ? cgdval : "",
				noassert ? "-noassert" : "",
				nodebug ? "-nodebug" : "",
				argv[ ac ] ? argv[ ac ] : "",
                cpp_ofname, n2c_ofname );
			if( cgdopt ) fprintf( stderr, "nab2c cmd: %s\n", cmd );
            status = system( cmd );

			unlink( cpp_ofname );
			if( status != 0 ){
				unlink( n2c_ofname );
				fprintf( stderr, "nab2c failed!\n" ); exit(1);
			}
		}
	}
	unlink( n2c_ofname );
}

static	void	cc( int copt, int cgdopt, int oopt, char ofname[],
	int argc, char *argv[], char nfmask[] )
/*
int	copt;
int	cgdopt;
int	oopt;
char	ofname[];
int	argc;
char	*argv[];
char	nfmask[];
*/
{
	int	ac;
    char    *dotp;
	char	*cmd, word[ 1024 ];
	int     cmd_sz;
	int 	status;

    char *nabhome = NABHOME;
	cmd_sz = 1024;
	cmd = malloc(cmd_sz);
	sprintf( cmd, "%s -Wno-parentheses -I%s/include", CC, nabhome );
	for( ac = 1; ac < argc; ac++ ){
		word[0] = '\0';
		if( nfmask[ ac ] ){
			dotp = strrchr( argv[ ac ], '.' );
			strcpy( &dotp[ 1 ], "c" );
			sprintf( word, " %s", argv[ ac ] );
		}else if( argv[ ac ] ){
			sprintf( word, " %s", argv[ ac ] );
		}
		if (strlen(cmd) + strlen(word) + 1 > cmd_sz) {
		    cmd_sz += 1024;
		    cmd = realloc(cmd, cmd_sz);
		}
		if (strlen(word) > 0) {
		    strcat( cmd, word );
		}
	}
	if( !copt ){
		sprintf( word, " -L%s/lib -lnab -lcifparse", nabhome );
		if (strlen(cmd) + strlen(word) + strlen(FLIBS) + 6 > cmd_sz) {
		    cmd_sz += strlen(word) + strlen(FLIBS) + 6;
		    cmd = realloc(cmd, cmd_sz);
		}
		strcat( cmd, word );
		strcat( cmd, " " );
		strcat( cmd, FLIBS );
		strcat( cmd, " -lm" );
	}
	if( cgdopt ) fprintf( stderr, "cc cmd: %s\n", cmd );
	status = system( cmd ); 
	if (status != 0) {
		fprintf(stderr,"cc failed!\n"); exit (1);
	}
	free(cmd);
}
