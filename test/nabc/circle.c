//
//  nab Program to create duplex B-DNA circles.
//  Usage:
//      a.out nbp dlk %gc
//
//  NBP is the number of base pairs which must be a non zero multiple of 10;
//  DLK is the change in the linking number. It can be +, 0, or -. It is
//  the number of superhelical turns placed in the circle assuming that
//  DNA makes a complete helix every 10 basepairs.
//

#include "nabc.h"

int main( int argc, char* argv[] )
{

int     b, nbp, dlk;
double  gc, seed, rnum;
double       rad, rise, twist, ttw;
MOLECULE_T    *m, *m1;
char      *sbase, *abase;

if( argc < 3 || argc > 4){
    fprintf( stderr, "usage: %s nbp dlk [ %%gc ]\n", argv[ 0 ] );
    exit( 1 );
}

nbp = atoi( argv[ 1 ] );
if( !nbp || nbp % 10 ){
    fprintf( stderr,
    "%s: Num. of base pairs must be multiple of 10\\n", argv[ 1 ] );
    exit( 1 );
}

// Get the delta linking number
dlk = atoi( argv[ 2 ] );

// Get the %gc, use 60% if not given
if( argc < 4 )
    gc = 0.6;
else
    gc = atof( argv[ 3 ] ) / 100.0;

rise = 3.38;
twist = ( nbp / 10 + dlk ) * 360.0 / nbp;
rad = 0.5 * rise / SIN( 180.0 / nbp );

//printf( "nbp, dlk, %%gc, rise, twist, rad: %d %d %5.2f %5.2f %5.2f %5.2f\n",
//         nbp, dlk, gc, rise, twist, rad );

m = newmolecule();
addstrand( m, "A" );
addstrand( m, "B" );
ttw = 0.0;

for( b=1; b <= nbp; b++ ){

    //  Create 1 standard B-dna W/C base pair
    if( ( rnum = rand2() ) <= 0.5 * gc ){
        sbase = "c";
        abase = "g";
    }else if( rnum <= gc ){
        sbase = "g";
        abase = "c";
    }else if( rnum <= 1.0 - 0.5 * gc ){
        sbase = "t";
        abase = "a";
    }else{
        sbase = "a";
        abase = "t";
    }

    m1 = wc_helix(
        sbase, "dna", abase, "dna", 2.25, -4.96, 0.0, 0.0, "" );

    //  Twist the base pair;  1st base pair has 0 twist
    if( b > 1 )
        transformmol( newtransform(0.,0.,0.,0.,0.,ttw), m1, NULL );

    //  Rotate the base from the xy plane into the xz plane
    // transformmol( newtransform(0.0, 0.0, 0.0, 90.0, 0.0, 0.0), m1, NULL );

    //  Put the base axes on the circle at y = 0.
    transformmol( newtransform(rad, 0., 0., 0., 0., 0.), m1, NULL );

    //  rotate to proper place on circle;  1st bp has a rotation of 0
    if( b > 1 )
        transformmol( newtransform( 0., 0., 0., 0., -360.0*(b-1)/nbp, 0. ), 
       m1, NULL );

    //  Insert the new base pair into growing circle; connect
    //  2d & subsequent base pairs
    mergestr( m, "A", "last", m1, "sense", "first" );
    mergestr( m, "B", "first", m1, "anti", "last" );

    if( b > 1 ){
        connectres( m, "A", b - 1, "O3'", b, "P" );
        connectres( m, "B", 1, "O3'", 2, "P" );
    }

    ttw = ttw + twist;
    if( ttw >= 360.0 ) ttw = ttw - 360.0;
}

// Close the circle
connectres( m, "A", nbp, "O3'", 1, "P" );
connectres( m, "B", nbp, "O3'", 1, "P" );

putpdb( "circ.pdb", m, NULL );
putbnd( "circ.bnd", m );

}
