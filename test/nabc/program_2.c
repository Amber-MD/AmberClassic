// Program 2 - Superimpose two DNA duplexes

#include "nabc.h"

FILE* nabout;

int main( int argc, char *argv[] )
{
	nabout = stdout; /*default*/

    MOLECULE_T *m, *mr;
    MATRIX_T  mat;
    double r;

    m = getpdb( "test.pdb", "" );
    mr = getpdb( "gcg10.pdb", "" );

    NAB_matcpy( mat, superimpose( m, "::C1'", mr, "::C1'" ) );
    transformmol( mat, m, NULL );
    putpdb( "test.sup.pdb", m, "" );
    // rmsd( m, "::C1'", mr, "::C1'", r );
    // printf( "rmsd = %8.3f\n", r );
    exit( 0 );
}
