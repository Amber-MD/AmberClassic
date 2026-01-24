#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "sff.h"
FILE* nabout;

//  Very primitive short rism minimizer, to illustrate the C-api; no argument-
//     checking, etc., etc.
//
//  Usage:  trismxmin  <parm-file>  <restart-file>
//       output to stdout


int main( int argc, char *argv[] )
{

   PARMSTRUCT_T *prm;   //  struct to hold info from a prmtop file
   XMIN_OPT_T xo;       //  options for the minimizer
   double *xyz,  *grad;
   double energy, grms;
   double start_time = 0.0;   // dummy, since this is minimization, not md

   nabout = stdout;    // change to redirect output (historical kludge)

//   options for the minimizer:

   xmin_opt_init( &xo );  // sets the default parameters;

   xo.maxiter = 5;                  // non-default minimization options:
   xo.grms_tol = 0.0005;
   xo.ls_maxatmov = 0.15;
   xo.print_level = 1;
   xo.method = 2;

//   read in the prmtop file and the coordinates:

   prm = rdparm( argv[1] );    // reads the prmtop file
   int natm = prm->Natom;
   xyz = malloc( 3 * natm * (sizeof(double)) );
   grad = malloc( 3 * natm * (sizeof(double)) );
   getxv( argv[2], natm, start_time, xyz, grad );  // reads a restart file

//   setup the force field parameters, and get an initial energy:

   mm_options( "ntpr=1, cut=99.0, diel=C " );
   mm_options( "ntpr=1, rism=1, closure=1, cut=9999.0" );
   mm_options( "xvvfile=../rism1d/spc-kh/spc.xvv.save" );
   mm_options( "verbose=0" );
   mm_options( "buffer=-1, ng=30,30,30, solvbox=15,15,15, solvcut=999" );
   mm_options( "tolerance=1e-11" );
   mm_options( "apply_rism_force=1, centering=2" );
   mm_options( "ntpr_rism=1, ntwrism=0" );
   mm_options( "uccoeff=-0.149818,-0.1136266,-0.00053163,0.0103954" );

   // nothing frozen or constrained for now:
   int* frozen = parseMaskString( "@ZZZ", prm, xyz, 2 );
   int* constrained = parseMaskString( "@ZZZ", prm, xyz, 2 );

   mme_init_sff( prm, frozen, constrained, NULL, NULL );
   int verbose = -1;   // historical flag to give more verbose output
   energy = mme( xyz, grad, &verbose );

//   run the minimization:

   char title[] = " non-rattle minimization";
   energy = xmin( mme,  &natm, xyz, grad,  &energy,  &grms,  &xo );
// putxv( argv[3], title, natm, start_time, xyz, xyz );
   energy = mme( xyz, grad, &verbose );

//  optional: get solvent distribution at the end:
//   mm_options( "guvfile=g.xmin, apply_rism_force=0, ntwrism=1" );
//   mme_init_sff( prm, frozen, constrained, NULL, NULL );
//   mme( xyz, grad, &verbose );

}

