#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "sff.h"

#ifdef MPI
   int   mpierror( int );
   int   mpifinalize( void );
   int   mpiinit( int*, char**, int*, int* );
   int mytaskid;
   int numtasks;
#else
   int mytaskid = 0;
   int numtasks = 1;
#endif

FILE* nabout;

int main( int argc, char *argv[] )
{

//  Simple example of MD: argv[1] has the prmtop filename;
//                        argv[2] has a netcdf restart filename;
//     beware: no error checking is done -- this is just a sample driver

   nabout = stdout; //  default; otherwise set to a file pointer of your choice

#ifdef MPI
   if ( mpiinit(&argc, argv, &mytaskid, &numtasks) != 0 ) {
      printf("Error in mpiinit!");
      fflush(stdout);
      exit(1);
   }
#endif

//   read in the prmtop file and the coordinates:

   PARMSTRUCT_T* prm = rdparm( argv[1] );    // reads the prmtop file
   int natm = prm->Natom;
   int natm3 = 3 * natm;
   double* x = malloc( natm3 * (sizeof(double)) );
   double* f = malloc( natm3 * (sizeof(double)) );
   double* v = malloc( natm3 * (sizeof(double)) );

#if 1
   //   option: read a formatted restart file:
   int start_time;
   getxv( argv[2], natm, start_time, x, v );
#else
   //  (next lines could be bundled into a netcdfReadRestart() function)
   struct AmberNetcdf ain;
   int ier = netcdfLoad( &ain, argv[2] );
   // netcdfDebug( &ain );
   // netcdfInfo( &ain );
   netcdfGetFrame(  &ain, 0, x, NULL );
   netcdfGetVelocity( &ain, 0, v );
#endif

//   set up some md options:

   //  force-field related options:
   mm_options( "cut=20.0, rgbmax=20.0, ntpr=500, nsnb=10, gb=1, diel=C" );
   //  dynamics related options:
   mm_options( "gamma_ln=5.0, temp0=300.0, ntpr_md=100, rattle=1" );

   // nothing frozen or constrained for now:
   int frozen[natm], constrained[natm];
   memset( frozen, 0, natm * sizeof(int) );
   memset( constrained, 0, natm * sizeof(int) );

   mme_init_sff( prm, frozen, constrained, NULL, NULL );

//   run the md:

   md( natm3, 500, x, f, v, mme );
   // if( mytaskid == 0 )
   //    netcdfWriteRestart( "xfin.md2.x", natm, x, v, NULL, 1.0, 0.0 );

#ifdef MPI
   mpifinalize();
#endif

}
