//
//   Example 1: thioredoxin
//   Inputs:     2trx.parm7, 2trx.rst7
//   Process:    Minimization, Heating, Equilibration and Production
//   Output:     2trx.log, 2trx.traj
// 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sff.h"
FILE* nabout;

int main( int argc, char *argv[] )
{

   nabout = stdout;     // need to figure out how to not need this!
   PARMSTRUCT_T *prm;   //  struct to hold info from a prmtop file

   // load parameter-topology file
   prm = rdparm("2trx.parm7");

   // allocate storage for coord (xyz), force (f_xyz), vel (v), constraint (x0)
   double* xyz = malloc( sizeof(double) * 3 * prm->Natom );
   double* f_xyz = malloc( sizeof(double) * 3 * prm->Natom );
   double* v = malloc( sizeof(double) * 3 * prm->Natom );
   double* x0 = malloc( sizeof(double) * 3 * prm->Natom );

   // initialize values of coordinates
   getxv( "2trx.rst7", prm->Natom, 0.0, xyz, v );
   memcpy( x0, xyz, 3 * prm->Natom * sizeof(double));  // use initial coords as restaint target
   for( int i=0; i<10; i++ )
      fprintf( stderr, "%10.3f   %10.3f\n", xyz[i], x0[i] );

   // Step1: Minimization:
   mm_options("cut=15, hcp=1, hcp_h1=15, hcp_h2=999, hcp_h3=999, \
     diel=C, gb=5, gbsa=1, rgbmax=15, wcons=5.0, kappa=0.125, nsnb=1, ntpr=10");

   // nothing frozen, but all atoms constrained:
   int* frozen = parseMaskString( "@ZZZ", prm, xyz, 2 );
   int* constrained = parseMaskString( ":*", prm, xyz, 2 );

   mme_init_sff( prm, frozen, constrained, x0, NULL );

   int iter = -1;   // historical flag to give more verbose output
   mme( xyz, f_xyz, &iter );

   double dgrad=0.1;
   double fret = 0.0;
   double dfpred = 10.0;
   int maxiter = 200;
   int ier = conjgrad(xyz, &prm->Natom, &fret, mme, &dgrad, &dfpred, &maxiter); 
   printf("conjgrad return code is  %d\n", ier);


#if 0

/* Step2: Heating  */
ier = mm_options("cut=15, hcp=1, hcp_h1=15, hcp_h2=999, hcp_h3=999, \
   diel=C, gb=5, gbsa=1, rgbmax=15, wcons=1.0, kappa=0.125, tempi=0.0, temp0=300.0, \
   gamma_ln=50, dt=0.001, nsnb=1, ntpr=100, ntpr_md=100");
ier = setmol_from_xyz(m, NULL, xyz);
ier = setxyz_from_mol(m, "::", x0);
ier = mme_init(m, NULL, NULL, x0, NULL);
ier = md(3*m.natoms, 1000, xyz, f_xyz, v, mme);

/* Step3: Equilibration part 1 */
ier = mm_options("cut=15, hcp=1, hcp_h1=15, hcp_h2=999, hcp_h3=999, \
   diel=C, gb=5, gbsa=1, rgbmax=15, wcons=0.1, kappa=0.125, tempi=300.0, temp0=300.0, \
   gamma_ln=50, dt=0.001, nsnb=1, ntpr=100, ntpr_md=100");
ier = setmol_from_xyz(m, NULL, xyz);
ier = setxyz_from_mol(m, "::", x0);
ier = mme_init(m, NULL, NULL, x0, NULL);
ier = md(3*m.natoms, 1000, xyz, f_xyz, v, mme);

/* Step4: Equilibration part 2 */
ier = mm_options("cut=15, hcp=1, hcp_h1=15, hcp_h2=999, hcp_h3=999, \
   diel=C, gb=5, gbsa=1, rgbmax=15, wcons=0.01, kappa=0.125, tempi=300.0, temp0=300.0, \
   gamma_ln=50, dt=0.001, nsnb=1, ntpr=100, ntpr_md=100");
ier = setmol_from_xyz(m, NULL, xyz);
ier = setxyz_from_mol(m, "::", x0);
ier = mme_init(m, NULL, NULL, x0, NULL);
ier = md(3*m.natoms, 1000, xyz, f_xyz, v, mme);

/* Step5: Production */
ier = mm_options("cut=15, hcp=1, hcp_h1=15, hcp_h2=999, hcp_h3=999, \
   diel=C, gb=5, gbsa=1, rgbmax=15, wcons=0.0, kappa=0.125, tempi=300.0, temp0=300.0, \
   gamma_ln=50, dt=0.001, nsnb=1, nscm=100, ntpr=1000, ntpr_md=1000, ntwx=1000");
/* open trajectory file */
traj = fopen("2trx.traj", "w");
ier = mme_init(m,NULL,"::ZZZ",dummy,traj);
ier = md(3*m.natoms, 100000, xyz, f_xyz, v, mme);

#endif

}
