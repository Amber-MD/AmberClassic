#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabc.h"
FILE* nabout;

static MOLECULE_T *m, *m2;


int main( int argc, char *argv[] )
{
    nabout = stdout;

    m = getpdb( "gbrna.pdb", NULL );
    putcif( "gbrna.cif", "1DAC", m );

    m2 = getcif( "gbrna.cif", "1DAC" );
    putpdb( "gbrna2.pdb", m2, "-nobocc, -wwpdb" );
}
