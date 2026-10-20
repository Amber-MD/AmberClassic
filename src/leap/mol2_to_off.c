/*
 *	File:	utilTripos.c
 *
 ************************************************************************
 *                                                                      *
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *	Description:
 *		Read a TRIPOS file and write OFF files.
 */



#include	"basics.h"

#include	"classes.h"

#include	"library.h"

#include	"tripos.h"




/*
 *	main
 *
 *	Author:	Christian Schafmeister (1991)
 */
int main( argc, argv )
int             argc;
char*           argv[];
{
FILE*           fIn;
UNIT		uUnit;
STRING		sName, sLibrary;
LIBRARY		lLibrary;


    BasicsInitialize();

    if ( argc != 2 ) {
        fprintf( stderr, "Usage: %s TRIPOS_FILE\n", 
			argv[0] );
        exit(0);
    }

    fIn = FOPENCOMPLAIN( argv[1], "r" );
    if ( fIn == NULL ) {
        printf( "Could not open file: %s\n", argv[1] );
        exit(1);
    }


	/* Read a TRIPOS UNIT */

    while ( uUnit = uTriposReadUnit(fIn) ) {

		/* Rename the UNIT so that leap doesnt choke on it */

	strcpy( sName, "u" );
	strcat( sName, sContainerName(uUnit) );
	ContainerSetName( uUnit, sName );
	printf( "Read unit: %s\n", sContainerName(uUnit) );

	strcpy( sLibrary, sContainerName(uUnit) );
	strcat( sLibrary, ".lib" );

	lLibrary = lLibraryOpen( sLibrary, OPENREADWRITE );
	LibrarySave( lLibrary, sContainerName(uUnit), (OBJEKT) uUnit, NULL );
	LibraryClose( &lLibrary );
    }

    fclose( fIn );
    return 0;
}
