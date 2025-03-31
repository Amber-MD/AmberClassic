#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "sff.h"
int mytaskid, numtasks;
FILE* nabout;

int main( int argc, char* argv[] )
{

   nabout = stdout;  // default; change if desired

   struct xmin_opt xo;  // no pointer to to allocate memory
   struct lmod_opt lo;

   int *lig_start,  *lig_end,  *lig_cent;
   double *xyz,  *grad,  *conflib,  *lmod_trajectory;
   double *tr_min,  *tr_max,  *rot_min,  *rot_max;

   lmod_opt_init(  &lo,  &xo );

   lo.niter = 3;
   lo.nconf = 10;
   lo.mc_option = 2;
   lo.nof_lmod_steps = 5;
   lo.random_seed = 99;
   lo.print_level = 2;

   xo.ls_maxatmov = 1.500000E-01;

   PARMSTRUCT_T* prm = rdparm( argv[1] );
   int natm =  prm->Natom;

   xyz =  malloc( 3 * natm * sizeof(double) );
   grad = malloc( 3 * natm * sizeof(double) );
   conflib = malloc( 3 * lo.nconf * natm * sizeof(double) );
   lmod_trajectory = malloc( 3 * natm * (lo.niter+1) * sizeof(double) );
   double start_time = 0.;
   getxv( argv[2], natm, start_time, xyz, grad );  // reads a restart file

   mm_options( "ntpr=99999, gb=0, cut=999.0, nsnb=9999, diel=R " );

   // nothing frozen or constrained for now:
   int frozen[natm], constrained[natm];
   memset( frozen, 0, natm * sizeof(int) );
   memset( constrained, 0, natm * sizeof(int) );

   mme_init_sff( prm, frozen, constrained, NULL, NULL );

   int iter = -1;   // historical flag to give more verbose output
   double energy = mme( xyz, grad, &iter );


   double glob_min_energy = lmod(  &natm, xyz, grad,  &energy, conflib, 
       lmod_trajectory, lig_start, lig_end, lig_cent, 
       tr_min, tr_max, rot_min, rot_max,  &xo,  &lo );

   if( mytaskid == 0 )
      printf( "\nGlob. min. E         = %12.3lf kcal/mol\n", glob_min_energy );

}
