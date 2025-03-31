#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Grid.h"
#include "BSpline.h"
#include "Constants.h"
#include "CrdManip.h"
#include "pmeRecip.h"
#include "mdgxVector.h"
#include "Matrix.h"
#include "pmeRecip.h"
#include "Timings.h"
#include "fftw3.h"

#include "TrajectoryDS.h"
#include "TopologyDS.h"

//-----------------------------------------------------------------------------
// Factor2357: determine the smallest available integer that is greater than
//             the specified value x and is a multiple of 2, 3, 5, and 7.
//                                                                      
// Arguments:                                                           
//   x:     see function description above                              
//   min2:  the minimum number of factors of 2 that the result must have
//-----------------------------------------------------------------------------
int Factor2357(double x, int min2)
{
  int n2, n3, n5, n7, lim2, lim3, lim5, lim7;
  int fac2, fac3, fac5, fac7, gnum, mingnum;
  double lg;

  lg = log(x);
  lim2 = lg/log(2.0) + 1;
  lim3 = lg/log(3.0) + 1;
  lim5 = lg/log(5.0) + 1;
  lim7 = lg/log(7.0) + 1;
  mingnum = 210.0*x + 1.0e-6;
  for (n2 = min2; n2 < lim2; n2++) {
    fac2 = pow(2.0, n2)+1.0e-8;
    for (n3 = 0; n3 < lim3; n3++) {
      fac3 = pow(3.0, n3)+1.0e-8;
      for (n5 = 0; n5 < lim5; n5++) {
	fac5 = pow(5.0, n5)+1.0e-8;
	for (n7 = 0; n7 < lim7; n7++) {
	  fac7 = pow(7.0, n7)+1.0e-8;
	  gnum = fac2*fac3*fac5*fac7;
	  if (gnum >= x && gnum < mingnum) {
	    mingnum = gnum;
	  }
	}
      }
    }
  }

  return mingnum;
}

//-----------------------------------------------------------------------------
// SetMeshDims: set the grid dimensions according to the standard "as close to
//              1A as possible, without going over" rule.  Like PMEMD, there
//              is a requirement that at least one (in this program, the final)
//              dimension be a multiple of 2, and multiples of 2, 3, and 5 are
//              sought.                    
//                                                                      
// Arguments:                                                           
//   ng:      number of grid points in this dimension                   
//   gdim:    simulation cell dimensions                                
//-----------------------------------------------------------------------------
void SetMeshDims(int* ng, double* gdim)
{
  int i, min2;

  for (i = 0; i < 3; i++) {
    min2 = (i == 2) ? 1 : 0;
    ng[i] = (ng[i] > 0) ? Factor2357(ng[i], min2) : Factor2357(gdim[i], min2);
  }
}

//-----------------------------------------------------------------------------
// GammaSum: this function computes a quantity needed for the B mesh    
//           prefactors.  This method was contributed by Tom Darden as   
//           part of the SANDER code.                                    
//                                                                      
// Arguments:                                                           
//   m:    the distance from 0 assuming that the box wraps around (so, m
//         ranges from -32 to 31 on a grid of 64 points).               
//   ng:   the size of the grid in this dimension.                      
//   ordr: the interpolation order                                      
//-----------------------------------------------------------------------------
static double GammaSum(int m, int ng, int ordr)
{
  int i;
  double gsum, invng, x;

  gsum = 1.0;

  if (m != 0) {
    invng = 1.0/ng;
    x = PI*m*invng;
    for (i = 1; i <= 50; i++) {
      gsum = gsum + pow(x/(x + PI*i), ordr) + pow(x/(x - PI*i), ordr);
    }
  }

  return gsum;
}

//-----------------------------------------------------------------------------
// LoadPrefac: this routine loads prefactors for the B mesh, which can be used
//             whether the unit cell is orthorhombic or not.  In other words,
//             the B mesh is the same for a given size of grid and a given
//             interpolation order, regardless of the dimensions of the unit
//             cell.                             
//                                                                      
// Arguments:                                                           
//   ordr: the interpolation order (here, ordr is a single integer, as this
//         routine only load prefactors for B in one dimension)    
//   ng:   the size of the mesh in this dimension                       
//-----------------------------------------------------------------------------
double* LoadPrefac(int ordr, int ng)
{
  int i, j, ngc, ordr2;
  double twopi, invng, lambda, gsum, gsum2, sum1, sum2, arg;
  double bsparr[ordr];
  double* bspmod;
  double* bv;

  // Allocate memory
  bspmod = (double*)calloc(ng, sizeof(double));
  bv = (double*)calloc(ng, sizeof(double));

  // Load the coefficients for an on-point B-Spline
  for (i = 0; i < ordr; i++) {
    bsparr[i] = BSpln(i+1.0, ordr);
  }

  // Compute the moduli of the DFT
  twopi = 2.0*PI;
  invng = 1.0/ng;
  for (i = 0; i < ng; i++) {
    sum1 = 0.0;
    sum2 = 0.0;
    for (j = 0; j < ordr; j++) {
      arg = twopi*i*j*invng;
      sum1 = sum1 + bsparr[j]*cos(arg);
      sum2 = sum2 + bsparr[j]*sin(arg);
    }
    bspmod[i] = sum1*sum1 + sum2*sum2;
  }

  // Fix the case where Euler exp. spline interpolation fails
  for (i = 0; i < ng; i++) {
    if (bspmod[i] < 1.0e-7) {
      if (i > 0 && i < ng) {
	bspmod[i] = 0.5*(bspmod[i-1] + bspmod[i+1]);
      }
      else {

	// This should never happen, but if it does we set bspmod
	// very high so that the inverse will be very small.     
	bspmod[i] = 1.0e15;
	printf("LoadPrefac >> Warning!  DFT modulus at index %d out of %d "
	       "is near zero!\n", i, ng);
      }
    }
  }

  // Optimize the lambda coefficient
  ngc = ng/2;
  ordr2 = 2*ordr;
  for (i = 0; i < ng; i++) {
    j = (i > ngc) ? i-ng : i;
    gsum = GammaSum(j, ng, ordr);
    gsum2 = GammaSum(j, ng, ordr2);
    lambda = gsum/gsum2;
    bv[i] = lambda*lambda/bspmod[i];
  }

  // Free allocated memory
  free(bspmod);

  return bv;
}

//-----------------------------------------------------------------------------
// LoadMVec: load the m-vector that will be used in mesh computation.  The
//           m-vector is loaded without normalization by the unit cell
//           dimensions to work with the non-orthorhombic unit cell case.
//                                                                      
// Arguments:                                                           
//   ng:  the mesh size in this dimension                               
//-----------------------------------------------------------------------------
double* LoadMVec(int ng)
{
  int i;
  double* mv;

  mv = (double*)malloc(ng*sizeof(double));
  for (i = 0; i < ng; i++) {
    mv[i] = (i <= ng/2) ? i : i-ng;
  }

  return mv;
}

//-----------------------------------------------------------------------------
// LoadMVecShift: load the shifted m-vector that will be used in mesh   
//                computations.  The m-vector is loaded without normalization
//                by the unit cell dimensions to work with the non-orthorhombic
//                unit cell case.                  
//                                                                      
// Arguments:                                                           
//   ng:  the mesh size in this dimension                               
//-----------------------------------------------------------------------------
static double* LoadMVecShift(int ng)
{
  int i, hng, idx;
  double* mvs;

  mvs = (double*)malloc(ng*sizeof(double));
  hng = ng/2;
  if (2*hng < ng) {
    hng += 1;
  }
  for (i = 0; i < ng; i++) {
    idx = (ng - i) % ng;
    mvs[i] = (idx <= ng/2) ? idx : idx-ng;
  }

  return mvs;
}

//-----------------------------------------------------------------------------
// TestUnitCellOrtho: test whether the unit cell is orthorhombic.  This 
//                    routine returns 1 if the unit cell is orthorhombic
//                    and 0 otherwise.                                  
//                                                                      
// Arguments:                                                           
//   gdim:  the unit cell dimensions                                    
//-----------------------------------------------------------------------------
int TestUnitCellOrtho(double* gdim)
{
  int i, isorthog;

  isorthog = 1;
  for (i = 3; i < 6; i++) {
    if (fabs(gdim[i] - 0.5*PI) > 1.0e-8) {
      isorthog = 0;
      break;
    }
  }

  return isorthog;
}

//-----------------------------------------------------------------------------
// SelfEnergyCorrection: pre-compute the electrostatic self energy correction
//                       for all charges in the system.      
//                                                                      
// Arguments:                                                           
//   ewcoeff:  the Ewald coefficient                                    
//   tp:       the topology                                             
//-----------------------------------------------------------------------------
double SelfEnergyCorrection(double ewcoeff, prmtop *tp)
{
  int i;
  double ecorr;
  double *qtmp;

  ecorr = 0.0;
  qtmp = tp->Charges;
  for (i = 0; i < tp->natom; i++) {
    ecorr += qtmp[i]*qtmp[i];
  }
  ecorr *= -BIOQ*ewcoeff/sqrt(PI);

  return ecorr;
}

//-----------------------------------------------------------------------------
// CreateBCKit: create a kit for assembling the BC mesh (the Fourier transform
//              of the reciprocal space pair potential).  The naming
//              conventions used in this function follow those used in Essmann
//              et al., J. Chem. Phys 103:8577.  This kit is intended for use
//              with STANDARD Smooth Particle Mesh Ewald transformations.  The
//              transformations are done in-place, for one thing, whereas
//              Multi-Level Ewald requires that the initial charge density be
//              preserved as the potential is accumulated.
//
// Arguments:                                                           
//   rcinp:  reciprocal space calculation control data                  
//   Q:      three-dimensional charge array                             
//   crd:    particle coordinates                                       
//   tp:     systemp topology                                           
//-----------------------------------------------------------------------------
bckit CreateBCKit(reccon *rcinp, dbook *Q, coord *crd, prmtop *tp,
		  unsigned int plan)
{
  int i;
  int *ng, *ordr;
  double tpS, pivol;
  double *gdim;
  bckit PPk;

  // We use real-to-complex and complex-to-real in-place transforms
  if (plan != FFTW_UNALIGNED) {
    PPk.forwplan = fftw_plan_dft_r2c_3d(Q->pag, Q->row, Q->col, Q->data,
					Q->fdata, plan);
    PPk.backplan = fftw_plan_dft_c2r_3d(Q->pag, Q->row, Q->col, Q->fdata,
					Q->data, plan);
    PPk.plans = 1;
  }
  else {
    PPk.plans = 0;
  }

  // Unpack the input data structure
  ng = rcinp->ng;
  ordr = rcinp->ordr;
  gdim = crd->gdim;

  // Pre-compute the electrostatic self energy correction
  PPk.SelfEcorr = SelfEnergyCorrection(0.5/rcinp->S, tp);

  // Load B mesh prefactors
  PPk.Bx = LoadPrefac(ordr[0], ng[0]);
  PPk.By = LoadPrefac(ordr[1], ng[1]);
  PPk.Bz = LoadPrefac(ordr[2], ng[2]);

  // Load m values
  PPk.mx = LoadMVec(ng[0]);
  PPk.my = LoadMVec(ng[1]);
  PPk.mz = LoadMVec(ng[2]);
  PPk.mxs = LoadMVecShift(ng[0]);
  PPk.mys = LoadMVecShift(ng[1]);
  PPk.mzs = LoadMVecShift(ng[2]);

  // Allocate exponential tables
  PPk.Ex = (double*)malloc(ng[0]*sizeof(double));
  PPk.Ey = (double*)malloc(ng[1]*sizeof(double));
  PPk.Ez = (double*)malloc(ng[2]*sizeof(double));

  // Test unit cell orthogonality
  if (crd->isortho == 0) {
    return PPk;
  }

  // The rest of this code is optimizations for the orthorhombic
  // unit cell case.  If the cell were not orthorhombic, this   
  // function would have returned on the statement above.       

  // Load C mesh prefactors
  DVecMult(PPk.mx, ng[0], 1.0/gdim[0]);
  DVecMult(PPk.my, ng[1], 1.0/gdim[1]);
  DVecMult(PPk.mz, ng[2], 1.0/gdim[2]);
  pivol = BIOQ/(PI*gdim[0]*gdim[1]*gdim[2]);
  tpS = 2.0*PI*rcinp->S;
  tpS *= tpS;
  for (i = 0; i < ng[0]; i++) {
    PPk.Ex[i] = pivol*exp(-tpS*PPk.mx[i]*PPk.mx[i]);
  }
  for (i = 0; i < ng[1]; i++) {
    PPk.Ey[i] = exp(-tpS*PPk.my[i]*PPk.my[i]);
  }
  for (i = 0; i < ng[2]; i++) {
    PPk.Ez[i] = exp(-tpS*PPk.mz[i]*PPk.mz[i]);
  }

  return PPk;
}

//-----------------------------------------------------------------------------
// UpdateBCKit: update a kit for assembling the BC mesh (the Fourier transform
//              of the reciprocal space pair potential).  This routine builds
//              on CreateBCKit (above), and changes those values that must be
//              changed given a rescaling of the periodic box (for
//              constant-pressure simulations).       
//                                                                      
// Arguments:                                                           
//   Identical to CreateBCKit above, with the addition of the reciprocal
//   space pair potential kit PPk included as an input argument rather  
//   than an output.                                                    
//-----------------------------------------------------------------------------
void UpdateBCKit(reccon *rcinp, coord *crd, prmtop *tp, bckit *PPk)
{
  int i;
  int *ng;
  double tpS, pivol;
  double *gdim;

  // If the unit cell is non-orthorhombic,
  // nothing in fact needs to be done.    
  if (crd->isortho == 0) {
    return;
  }

  // Unpack the input data structure
  ng = rcinp->ng;
  gdim = crd->gdim;

  // If the unit cell is orthorhombic, update
  // the m values and exponential tables.    
  free(PPk->mx);
  free(PPk->my);
  free(PPk->mz);
  PPk->mx = LoadMVec(ng[0]);
  PPk->my = LoadMVec(ng[1]);
  PPk->mz = LoadMVec(ng[2]);
  DVecMult(PPk->mx, ng[0], 1.0/gdim[0]);
  DVecMult(PPk->my, ng[1], 1.0/gdim[1]);
  DVecMult(PPk->mz, ng[2], 1.0/gdim[2]);
  pivol = BIOQ/(PI*gdim[0]*gdim[1]*gdim[2]);
  tpS = 2.0*PI*rcinp->S;
  tpS *= tpS;
  for (i = 0; i < ng[0]; i++) {
    PPk->Ex[i] = pivol*exp(-tpS*PPk->mx[i]*PPk->mx[i]);
  }
  for (i = 0; i < ng[1]; i++) {
    PPk->Ey[i] = exp(-tpS*PPk->my[i]*PPk->my[i]);
  }
  for (i = 0; i < ng[2]; i++) {
    PPk->Ez[i] = exp(-tpS*PPk->mz[i]*PPk->mz[i]);
  }

  // Pre-compute the electrostatic self energy correction  
  // FIX ME!  This doesn't need to happen unless rcinp->S  
  // (the width of the Gaussian charge spread) has changed.
  PPk->SelfEcorr = SelfEnergyCorrection(0.5/rcinp->S, tp);
}

//-----------------------------------------------------------------------------
// DestroyBCKit: free all memory associated with the pair potential kit PPk
//-----------------------------------------------------------------------------
void DestroyBCKit(bckit *PPk)
{
  free(PPk->Bx);
  free(PPk->By);
  free(PPk->Bz);
  free(PPk->mx);
  free(PPk->my);
  free(PPk->mz);
  free(PPk->mxs);
  free(PPk->mys);
  free(PPk->mzs);
  free(PPk->Ex);
  free(PPk->Ey);
  free(PPk->Ez);
  if (PPk->plans == 1) {
    fftw_destroy_plan(PPk->forwplan);
    fftw_destroy_plan(PPk->backplan);
  }
}

//-----------------------------------------------------------------------------
// PrepPME: the analog of the PrepMLE function (see mleRecip.c), this work is
//          encapsulated here to keep the main function tidy.   
//                                                                      
// Arguments:                                                           
//   rcinp:      reciprocal space control information                   
//   crd:        coordinates of the system (for box size information)   
//-----------------------------------------------------------------------------
void PrepPME(cellgrid *CG, reccon *rcinp, coord *crd)
{
  int i;

  SetMeshDims(rcinp->ng, crd->gdim);
  if (CG->tid == 0) {
    rcinp->QL = (dbook*)malloc(CG->nthreads*sizeof(dbook));
    for (i = 1; i < CG->nthreads; i++) {
      rcinp->QL[i] = CreateDbook(rcinp->ng[0], rcinp->ng[1], rcinp->ng[2], 1);
    }
  }
  else {
    rcinp->QL = (dbook*)malloc(2*sizeof(dbook));
    rcinp->QL[1] = CreateDbook(rcinp->ng[0], rcinp->ng[1], rcinp->ng[2], 1);
  }
  rcinp->QL[0] = CreateDbook(rcinp->ng[0], rcinp->ng[1], rcinp->ng[2], 1);
}

//-----------------------------------------------------------------------------
// DestroyRecCon: free all memory associated with a reciprocal space    
//                control data structure.
//
// Arguments:
//   rcinp    : reciprocal space input control data
//   CG       : the cell grid
//-----------------------------------------------------------------------------
void DestroyRecCon(reccon *rcinp, cellgrid *CG)
{

  int i;
  
  if (rcinp->nlev == 1) {

    // The master process will have allocated     
    // buffers for grids from all other processes,
    // but other processes will have only one grid
    if (CG->tid == 0) {
      for (i = 1; i < CG->nthreads; i++) {
	DestroyDbook(&rcinp->QL[i]);
      }
    }
    else {
      DestroyDbook(&rcinp->QL[1]);
    }
    DestroyDbook(&rcinp->QL[0]);
  }
  if (rcinp->nlev > 1) {
    free(rcinp->SPrv);
    free(rcinp->SPcv);
    for (i = 0; i < rcinp->nlev; i++) {
      DestroyDbook(&rcinp->QL[i]);
      DestroyDbook(&rcinp->Urec[i]);
      fftw_destroy_plan(rcinp->forwplan[i]);
      fftw_destroy_plan(rcinp->backplan[i]);
    }
    DestroyDbook(&rcinp->QL[rcinp->nlev]);
    free(rcinp->Urec);
  }

  free(rcinp->QL);
  free(rcinp->ng);
}


//-----------------------------------------------------------------------------
// DestroyAdvancedRecCon: free memory associated with a reciprocal space
//                control data structure except ng which is only        
//                initialized once by InitBasicReccon                   
//
// Arguments:
//   rcinp    : reciprocal space input control data
//   CG       : the cell grid
//-----------------------------------------------------------------------------
void DestroyAdvancedRecCon(reccon *rcinp, cellgrid *CG)
{
  int i;
    if (CG->tid == 0) {
      for (i = 1; i < CG->nthreads; i++) {
		DestroyDbook(&rcinp->QL[i]);
      }
    }
    else { 
      DestroyDbook(&rcinp->QL[1]);
    }
    DestroyDbook(&rcinp->QL[0]);
  free(rcinp->QL);
} 


// BLOCK 1: Energy and virial are not required
#define NEEDENERGY 0
#define NEEDVIRIAL 0
#include "Convolution.c"
#undef NEEDVIRIAL
#undef NEEDENERGY
// END BLOCK 1

// BLOCK 2: Energy is required
#define NEEDENERGY 1

// BLOCK 2a: Virial is unecessary
#define NEEDVIRIAL 0
#include "Convolution.c"
#undef NEEDVIRIAL

// Block 2b: Virial is also needed
#define NEEDVIRIAL 1
#include "Convolution.c"
#undef NEEDVIRIAL

#undef NEEDENERGY
// END BLOCK 2
