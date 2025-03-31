#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sff.h"
FILE* nabout;

int main( int argc, char* argv[] )
{

   nabout = stdout;

   PARMSTRUCT_T* prm = rdparm( argv[1] );
   int natm = prm->Natom;
   double* xyz = malloc( prm->Nat3 * (sizeof(double)) );
   double* v = malloc( prm->Nat3 * (sizeof(double)) );
   double start_time = 0.0;

   getxv( argv[2], natm, start_time, xyz, v );  // reads a restart file

   mm_options( "cut=999., diel=C, gb=1, dielc=1.0, rgbmax=999. ,kappa=.10395" );

   // nothing frozen or constrained for now:
   int* frozen = parseMaskString( "@ZZZ", prm, xyz, 2 );
   int* constrained = parseMaskString( "@ZZZ", prm, xyz, 2 );

   mme_init_sff( prm, frozen, constrained, NULL, NULL );

   int iter = -1;
   mme( xyz, v, &iter );

   int ier = nmode( xyz, prm->Nat3, mme2, 20, 1, 0.0, 0.0, 0);
   printf("nmode returns %d\n", ier );
   // mme2_timer();

}

