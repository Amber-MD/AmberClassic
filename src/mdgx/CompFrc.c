#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CompFrc.h"
#include "Constants.h"
#include "Matrix.h"

#include "SpecialMath.h"

//-----------------------------------------------------------------------------
// SplineTab: this function tabulates a piecewise cubic spline for
//            interpolating the function u(r).                          
//                                                                      
// Arguments:                                                           
//   u:    values of the function to be interpolated, given at regular  
//         intervals (dr) of the argument (r)                           
//   du:   derivative values of the function to be interpolated         
//   nbin: the number of bins (size of the arrays u and du)             
//   dr:   the increment in the function argument r                     
//-----------------------------------------------------------------------------
CSpln* SplineTab(double* u, double* du, int nbin, double dr)
{
  int i;
  double r, r2, rpdr, rpdr2;
  double b[4];
  dmat A;
  CSpln* S;

  // Allocate memory
  S = (CSpln*)malloc(nbin*sizeof(CSpln));

  // Loop over all points and compute the spline
  A = CreateDmat(4, 4, 0);
  for (i = 0; i < nbin-1; i++) {

    // Set up the matrix equation
    r = i*dr;
    rpdr = r+dr;
    r2 = r*r;
    rpdr2 = rpdr*rpdr;
    b[0] = u[i];
    b[1] = u[i+1];
    b[2] = du[i];
    b[3] = du[i+1];
    A.map[0][0] = r2*r;
    A.map[0][1] = r2;
    A.map[0][2] = r;
    A.map[0][3] = 1.0;
    A.map[1][0] = rpdr2*rpdr;
    A.map[1][1] = rpdr2;
    A.map[1][2] = rpdr;
    A.map[1][3] = 1.0;
    A.map[2][0] = 3.0*r2;
    A.map[2][1] = 2.0*r;
    A.map[2][2] = 1.0;
    A.map[2][3] = 0.0;
    A.map[3][0] = 3.0*rpdr2;
    A.map[3][1] = 2.0*rpdr;
    A.map[3][2] = 1.0;
    A.map[3][3] = 0.0;

    // Solve the matrix equation
    AxbQRRxc(A, b, 0);
    BackSub(A, b);
    S[i].A = b[0];
    S[i].B = b[1];
    S[i].C = b[2];
    S[i].D = b[3];
  }

  // Free allocated memory
  DestroyDmat(&A);

  return S;
}

//-----------------------------------------------------------------------------
// SplineTab2: this function tabulates a piecewise cubic spline for     
//             interpolating the function u(r).  Unlike the SplineTab function
//             (above, however, this function merely seeks to minimize the
//             error in the approximation to u(r), not to produce a function
//             with a continuous derivative.         
//-----------------------------------------------------------------------------
CSpln* SplineTab2(double* u, double* umid, int nbin, double dr)
{
  int i;
  double r, r2, rpdr, rpdr2, rptdr, rptdr2, rpttdr, rpttdr2;
  double b[4];
  dmat A;
  CSpln* S;

  // Allocate memory
  S = (CSpln*)malloc(nbin*sizeof(CSpln));

  // Loop over all points and compute the spline
  A = CreateDmat(4, 4, 0);
  const double ONETHIRDdr = dr/3.0;
  const double TWOTHIRDSdr = 2.0*dr/3.0;
  for (i = 0; i < nbin-1; i++) {

    // Set up the matrix equation
    r = i*dr;
    rptdr = r + ONETHIRDdr;
    rpttdr = r + TWOTHIRDSdr;
    rpdr = r+dr;
    r2 = r*r;
    rptdr2 = rptdr*rptdr;
    rpttdr2 = rpttdr*rpttdr;
    rpdr2 = rpdr*rpdr;
    b[0] = u[i];
    b[1] = umid[2*i];
    b[2] = umid[2*i+1];
    b[3] = u[i+1];
    A.map[0][0] = r2*r;
    A.map[0][1] = r2;
    A.map[0][2] = r;
    A.map[0][3] = 1.0;
    A.map[1][0] = rptdr2*rptdr;
    A.map[1][1] = rptdr2;
    A.map[1][2] = rptdr;
    A.map[1][3] = 1.0;
    A.map[2][0] = rpttdr2*rpttdr;
    A.map[2][1] = rpttdr2;
    A.map[2][2] = rpttdr;
    A.map[2][3] = 1.0;
    A.map[3][0] = rpdr2*rpdr;
    A.map[3][1] = rpdr2;
    A.map[3][2] = rpdr;
    A.map[3][3] = 1.0;

    // Solve the matrix equation
    AxbQRRxc(A, b, 0);
    BackSub(A, b);
    S[i].A = b[0];
    S[i].B = b[1];
    S[i].C = b[2];
    S[i].D = b[3];
  }

  // Free allocated memory
  DestroyDmat(&A);

  return S;
}

//-----------------------------------------------------------------------------
// DirectSpaceR2: force tabulation for a Gaussian charge convolution scheme, as
//                used in the standard Smooth Particle Mesh Ewald method,
//                producing a cubic spline lookup table indexed by the squared
//                distance between two particles.
//
// Arguments:
//   range      : the range for the lookup table to span                
//   spc        : spacing of the lookup table                           
//   ewcoeff    : the Ewald coefficient, 0.5/alpha                      
//   contderiv  : flag to specify a continuous derivative (set to 1), or
//                not (set to zero)                                     
//   doerf      : flag to indicate that the erf(r)/r, and not erfc(r)/r,
//                should be splined.  Set to 1 to activate this feature;
//                default 0, do erfc(r)/r as needed by SPME.            
//-----------------------------------------------------------------------------
FrcTab DirectSpaceR2(double range, double spc, double ewcoeff, int contderiv,
		     int doerf)
{
  int i, j, ir;
  double r2, r, invr, invr2, ewr, rer, bexp, bfac, invdspc;
  double r2p, rp, invrp, ewrp, rerp, bexpp;
  double r2m, rm, invrm, ewrm, rerm, bexpm;
  double truefrc, truenrg, ntrpfrc, ntrpnrg;
  dmat u, du;
  FrcTab F;

  // Determine the number of bins and store other information
  F.rmax = range;
  F.dr = spc;
  F.ivdr = 1.0/spc;
  F.nbin = F.rmax*F.rmax*F.ivdr + 2;

  // Allocate scratch arrays
  u = CreateDmat(4, F.nbin, 0);
  if (contderiv == 1) {
    du = CreateDmat(2, F.nbin, 0);
  }
  else {
    du = CreateDmat(2, 2*F.nbin, 0);
  }

  // Pre-compute constants
  bfac = 2.0*ewcoeff/sqrt(PI);
  invdspc = 1.0/(0.002*spc);
  const double ONETHIRD = 1.0/3.0;
  const double TWOTHIRDS = 2.0/3.0;

  // Loop over all indices
  for (i = 0; i < F.nbin; i++) {

    // Basic information about the distance r
    r2 = (i*spc < 1.0e-8) ? 1.0e-8 : i*spc;
    r = sqrt(r2);
    invr = 1.0/r;

    // Electrostatic direct space force and potential
    invr2 = invr*invr;
    ewr = r*ewcoeff;
    if (doerf == 0) {
      rer = PerfcFast(ewr)*invr;
    }
    else {
      rer = (1.0-PerfcFast(ewr))*invr;
    }
    bexp = bfac*exp(-ewr*ewr);
    u.map[0][i] = -BIOQ*(rer + bexp)*invr2;
    u.map[1][i] = BIOQ*rer;

    // Make a spline with a continuous derivative
    if (contderiv == 1) {

      // r2 incremented by 0.001*spc
      r2p = r2+0.001*spc;
      rp = sqrt(r2p);
      invrp = 1.0/rp;
      ewrp = rp*ewcoeff;
      if (doerf == 0) {
	rerp = PerfcFast(ewrp)*invrp;
      }
      else {
	rerp = (1.0-PerfcFast(ewrp))*invrp;
      }
      bexpp = bfac*exp(-ewrp*ewrp);
      du.map[0][i] = -BIOQ*(rerp + bexpp)*invrp*invrp;
      du.map[1][i] = BIOQ*rerp;

      // r2 decremented by 0.001*spc
      if (i > 0) {
	r2m = r2-0.001*spc;
	rm = sqrt(r2m);
	invrm = 1.0/rm;
	ewrm = rm*ewcoeff;
	if (doerf == 0) {
	  rerm = PerfcFast(ewrm)*invrm;
	}
	else {
	  rerm = (1.0-PerfcFast(ewrm))*invrm;
	}
	bexpm = bfac*exp(-ewrm*ewrm);
	du.map[0][i] -= -BIOQ*(rerm + bexpm)*invrm*invrm;
	du.map[1][i] -= BIOQ*rerm;
      }
      else {
	for (j = 0; j < 2; j++) {
	  du.map[j][0] = 2.0*(du.map[j][0]-u.map[j][0]);
	}
      }
      for (j = 0; j < 2; j++) {
	du.map[j][i] *= invdspc;
      }
    }

    // Make the most accurate possible approximation of u(r)
    else {

      // The first target point is 1/3 of the way into the interval
      r2 = ((i+ONETHIRD)*spc < 1.0e-8) ? 1.0e-8 : (i+ONETHIRD)*spc;
      r = sqrt(r2);
      invr = 1.0/r;
      ewr = r*ewcoeff;
      if (doerf == 0) {
	rer = PerfcFast(ewr)*invr;
      }
      else {
	rer = (1.0-PerfcFast(ewr))*invr;
      }
      bexp = bfac*exp(-ewr*ewr);
      du.map[0][2*i] = -BIOQ*(rer + bexp)*invr*invr;
      du.map[1][2*i] = BIOQ*rer;

      // The second target point is 2/3 of the way into the interval
      r2 = ((i+TWOTHIRDS)*spc < 1.0e-8) ? 1.0e-8 : (i+TWOTHIRDS)*spc;
      r = sqrt(r2);
      invr = 1.0/r;
      ewr = r*ewcoeff;
      if (doerf == 0) {
	rer = PerfcFast(ewr)*invr;
      }
      else {
	rer = (1.0-PerfcFast(ewr))*invr;
      }
      bexp = bfac*exp(-ewr*ewr);
      du.map[0][2*i+1] = -BIOQ*(rer + bexp)*invr*invr;
      du.map[1][2*i+1] = BIOQ*rer;
    }
  }

  // Compute spline coefficients
  if (contderiv == 1) {
    F.dSD = SplineTab(u.map[0], du.map[0], F.nbin, F.dr);
    F.SD = SplineTab(u.map[1], du.map[1], F.nbin, F.dr);
  }
  else {
    F.dSD = SplineTab2(u.map[0], du.map[0], F.nbin, F.dr);
    F.SD = SplineTab2(u.map[1], du.map[1], F.nbin, F.dr);
  }

  // Estimate error
  r2 = 0.0;
  F.FitType = contderiv;
  F.fmaxerr = 0.0;
  F.fmaxrelerr = 0.0;
  F.fmaxerrloc = 0.0;
  F.fmaxrelerrloc = 0.0;
  F.umaxerr = 0.0;
  F.umaxrelerr = 0.0;
  F.umaxerrloc = 0.0;
  F.umaxrelerrloc = 0.0;
  for (i = 0; i < 10*(F.nbin - 2); i++) {
    r2 += 0.1*F.dr;
    if (r2 < 1.0) {
      continue;
    }
    r = sqrt(r2);
    ir = r2*F.ivdr;
    invr = 1.0/r;
    invr2 = invr*invr;
    ewr = r*ewcoeff;
    rer = PerfcFast(ewr)*invr;
    bexp = bfac*exp(-ewr*ewr);
    truefrc = -BIOQ*(rer + bexp)*invr2;
    truenrg = BIOQ*rer;
    ntrpfrc = (((F.dSD[ir].A*r2 + F.dSD[ir].B)*r2 + F.dSD[ir].C)*r2 +
	       F.dSD[ir].D);
    ntrpnrg = (((F.SD[ir].A*r2 + F.SD[ir].B)*r2 + F.SD[ir].C)*r2 +
               F.SD[ir].D);
    if (fabs(truefrc - ntrpfrc) > F.fmaxerr) {
      F.fmaxerr = fabs(truefrc - ntrpfrc);
      F.fmaxerrloc = r;
    }
    if (fabs(truefrc - ntrpfrc)/fabs(truefrc) > F.fmaxrelerr) {
      F.fmaxrelerr = fabs(truefrc - ntrpfrc)/fabs(truefrc);
      F.fmaxrelerrloc = r;
    }
    if (fabs(truenrg - ntrpnrg) > F.umaxerr) {
      F.umaxerr = fabs(truenrg - ntrpnrg);
      F.umaxerrloc = r;
    }
    if (fabs(truenrg - ntrpnrg)/fabs(truenrg) > F.umaxrelerr) {
      F.umaxrelerr = fabs(truenrg - ntrpnrg)/fabs(truenrg);
      F.umaxrelerrloc = r;
    }
  }

  // Free allocated memory
  DestroyDmat(&u);
  DestroyDmat(&du);

  return F;
}

//-----------------------------------------------------------------------------
// FreeFrcTab: free a force table.                                      
//                                                                      
// Arguments:                                                           
//   F:  the force table                                                
//-----------------------------------------------------------------------------
void FreeFrcTab(FrcTab *F)
{
  free(F->SD);
  free(F->dSD);
}
