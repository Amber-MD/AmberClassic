#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabc.h"
FILE* nabout;

int main( int argc, char *argv[] )
{
	nabout = stdout; /*default*/

    MOLECULE_T *m;

    // get pdb file and number of atoms
    m = getpdb( argv[1], "" );
    int natm = countmolatoms( m, NULL );

    // for historical reasons, setxyz_from_mol() takes a POINT_T*,
    //    whereas sff routines like putxv()/getxv(), etc. take a double*

    POINT_T* xyzp = malloc( sizeof(POINT_T) * natm );
    setxyz_from_mol( &m, NULL, xyzp );

    char* title = "sample title";
    double start_time = 0.0;
    putxv( argv[2], title, natm, start_time, (double *)xyzp, (double *)xyzp );
}
