#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sff.h"
FILE* nabout;

int main( int argc, char *argv[] )
{

   PARMSTRUCT_T *prm;   //  struct to hold info from a prmtop file
   int natm, iter;
   double energy, grms;
   double start_time = 0.0;   // dummy, since this is minimization, not md

   nabout = stdout;    // change to redirect output (historical kludge)

//   read in the prmtop file and the coordinates:

   prm = rdparm( argv[1] );    // reads the prmtop file
   natm = prm->Natom;

   double *xyz = malloc( 3 * natm * (sizeof(double)) );
   double *grad = malloc( 3 * natm * (sizeof(double)) );
   getxv( argv[2], natm, start_time, xyz, grad );  // reads a restart file

//   test an atom expression:

   int* mask = parseMaskString( argv[3], prm, xyz, 2 );

   for( int i=0; i<natm; i++ ){
      printf( "%d   %s  %d  \n", mask[i], prm->AtomNames[i], prm->AtomRes[i] );
   }

}
