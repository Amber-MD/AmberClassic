#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "sff.h"
#include "../src/rism-cuda/rism3d.h"
FILE* nabout;

//  Very primitive short rism-cuda minimizer, no argument-
//     checking, etc., etc.
//
//  Usage:  trism-cuda  <parm-file>  <restart-file>
//       output to stdout

REAL_T mme_rism( REAL_T * x, REAL_T * f, int *iter ){

   REAL_T energy;
   RISM3D * rsystem;

   energy = mme( x, f, iter );

   rsystem = new RISM3D;
   rsystem -> initialize( "trpcage_c", "trpcage_s", false );

   // here we need to copy the current coordinates (passed in as x) into
   //   the appropriate place inside rsystem:
   for (int iu = 0; iu < rsystem -> su -> num; ++iu) {
      int n = iu * 4;
      int nsff = iu * 3;
      rsystem -> su -> r[n  ] = x[nsff  ];
      rsystem -> su -> r[n+1] = x[nsff+1];
      rsystem -> su -> r[n+2] = x[nsff+2];
   }
   rsystem -> su -> setup_cuda();

   rsystem -> cal_potential();
   rsystem -> iterate(0);

   // get exchem, pmv, pressure:
   double pmv = rsystem -> cal_pmv();
   double pressure = rsystem -> cal_pressure();
   double * xmu = new double[rsystem -> sv -> natv * 2];
   rsystem -> cal_exchem(xmu);
   double ibeta = avogadoro * boltzmann * rsystem -> sv -> temper;
    
   double xmua = 0.0;
   for (int iv = 0; iv < rsystem -> sv -> natv; ++iv) { xmua += xmu[iv]; }   
   double erism = ibeta * xmua / kcal2J;
   printf( "emm = %10.5f  erism = %10.5f  total = %10.5f\n", 
      energy, erism, energy + erism );
   energy += erism;
   delete[] xmu;

   // get gradient:
   double * du;
   du = new double[rsystem -> su -> num * 3];
   rsystem -> cal_grad(du);
   double dv = rsystem -> ce -> dv / kcal2J;
   for (int iu = 0; iu < rsystem -> su -> num; ++iu) {
      int num = iu * 3;
      f[num] += du[num  ] * dv;
      f[num] += du[num+1] * dv;
      f[num] += du[num+2] * dv;
   }
   delete[] du;
   rsystem -> su -> free_cuda();
   return energy;
}

int main( int argc, char *argv[] )
{
   PARMSTRUCT_T *prm;   //  struct to hold info from a prmtop file
   XMIN_OPT_T xo;       //  options for the minimizer
   double *xyz,  *grad, *v;
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
   xyz = (double *) malloc( 3 * natm * (sizeof(double)) );
   grad = (double *) malloc( 3 * natm * (sizeof(double)) );
   v = (double *) malloc( 3 * natm * (sizeof(double)) );

   getxv( argv[2], natm, start_time, xyz, v );  // reads a restart file

//   setup the force field parameters, and get an initial energy:

   mm_options( (char *) "ntpr=1, cut=99.0," );
   // mm_options( "xvvfile=../rism1d/spc-kh/spc.xvv.save" );
   // mm_options( "verbose=0" );
   // mm_options( "buffer=-1, ng=30,30,30, solvbox=15,15,15, solvcut=999" );
   // mm_options( "tolerance=1e-7" );
   // mm_options( "apply_rism_force=1, centering=2" );
   // mm_options( "ntpr_rism=1, ntwrism=0" );
   // mm_options( "uccoeff=-0.149818,-0.1136266,-0.00053163,0.0103954" );
   mm_options( (char *) "ntpr_md=1" );

   // nothing frozen or constrained for now:
   int* frozen = parseMaskString( (char *) "@ZZZ", prm, xyz, 2 );
   int* constrained = parseMaskString( (char *) "@ZZZ", prm, xyz, 2 );

   mme_init_sff( prm, frozen, constrained, NULL, NULL );
   int dn=0;
   cudaSetDevice(dn);

// int verbose = -1;   // historical flag to give more verbose output
// energy = mme_rism( xyz, grad, &verbose );

//   run the simulation:

//    energy = xmin( mme_rism,  &natm, xyz, grad,  &energy,  &grms,  &xo );
// putxv( argv[3], title, natm, start_time, xyz, xyz );

   int natm3 = 3 * natm;
   md( natm3, 10, xyz, grad, v, mme_rism );


}

