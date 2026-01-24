#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "BSpline.h"

#include "pmeRecipDS.h"
#include "pmeDirectDS.h"
#include "CellManipDS.h"

//-----------------------------------------------------------------------------
// BSpln: this is the basic function for B-Spline computation.  It's not at all
//        optimal, but it is easy to code using the recursive definition of
//        B-Splines.                                      
//                                                                      
// Arguments:                                                           
//   x:    the value to compute on this spline                          
//   ordr: the order of the spline (the spline will have nonzero values 
//         for 0 < x < ordr)                                            
//-----------------------------------------------------------------------------
double BSpln(double x, int ordr)
{
  double om1, iom1, y;

  om1 = ordr-1.0;
  iom1 = 1.0/om1;
  if (ordr > 2) {
    y = (x*iom1)*BSpln(x,om1) + ((ordr-x)*iom1)*BSpln(x-1.0, ordr-1);
  }
  else if (ordr == 2) {
    if (x > 0.0 && x < 2.0) {
      y = 1.0 - fabs(x-1.0);
    }
    else {
      y = 0.0;
    }
  }
  else {
    printf("BSpln >> Error.  Order = %d.\n", ordr);
    exit(1);
  }

  return y;
}

//-----------------------------------------------------------------------------
// dBSpln: this is a basic function for B-Spline derivative computation.
//                                                                      
// Arguments:                                                           
//   x:    the value to compute on this spline                          
//   ordr: the order of the spline (the spline will have nonzero values 
//         for 0 <= x <= ordr)                                          
//-----------------------------------------------------------------------------
double dBSpln(double x, int ordr)
{
  double y;

  if (ordr > 2) {
    y = BSpln(x, ordr-1) - BSpln(x-1.0, ordr-1);
  }
  else {
    printf("dBSpln >> Error.  Order = %d.\n", ordr);
    exit(1);
  }

  return y;
}

//-----------------------------------------------------------------------------
// CreateBmap: create a data structure to hold all information necessary to
//             perform particle<-->mesh operations.  The structure is filled
//             by SplCoeff (see below).                          
//                                                                      
// Arguments:                                                           
//   rcinp:  the reciprocal space control data                          
//   natom:  the number of atoms in the system                          
//-----------------------------------------------------------------------------
bmap CreateBmap(reccon *rcinp, int natom)
{
  bmap pmmap;

  pmmap.natom = natom;
  pmmap.xcof = (bcof*)malloc(natom*rcinp->ordr[0]*sizeof(bcof));
  pmmap.ycof = (bcof*)malloc(natom*rcinp->ordr[1]*sizeof(bcof));
  pmmap.zcof = (bcof*)malloc(natom*rcinp->ordr[2]*sizeof(bcof));
  pmmap.atmid = (int*)malloc(natom*sizeof(int));

  return pmmap;
}

//-----------------------------------------------------------------------------
// DestroyBmap: destroy a bmap data structure.                          
//                                                                      
// Arguments:                                                           
//   A:  the B-spline map to destroy                                    
//-----------------------------------------------------------------------------
void DestroyBmap(bmap *A)
{
  free(A->atmid);
  free(A->xcof);
  free(A->ycof);
  free(A->zcof);
}

//-----------------------------------------------------------------------------
// FillQMap: fill the charge map for this charge in this dimension.     
//                                                                      
// Arguments:                                                           
//   SP:    the array of B-spline coefficients and charge mesh indices  
//          in this dimension (the dimension of interest)             
//   uc:    the result of ceil() acting on the scaled fractional        
//          coordinate of the particle in the dimension of interest   
//   ng:    the number of mesh points in the dimension of interest      
//   ordr:  the interpolation order in the dimension of interest        
//-----------------------------------------------------------------------------
static void FillQMap(bcof* SP, int uc, int ng, int ordr)
{
  int i;
  double invng;

  // Charge map
  invng = 1.0/ng;
  for (i = 0; i < ordr; i++) {
    SP[i].m = uc + i;
    SP[i].m = SP[i].m - ng*floor(SP[i].m*invng);
  }
}

//-----------------------------------------------------------------------------
// Extend4: extend third order B-spline coefficients to fourth order and
//          compute derivatives.                                        
//                                                                      
// Arguments:                                                           
//   SP:    the array of B-spline coefficients and charge mesh indices in this 
//          dimension (the dimension of interest)             
//   w:     one minus the non-integer part of the scaled fractional coordinate
//   uc:    the result of ceil() acting on the scaled fractional coordinate
//          of the particle in the dimension of interest   
//   ng:    the number of mesh points in the dimension of interest      
//-----------------------------------------------------------------------------
static void Extend4(bcof* SP, double w, int uc, int ng)
{
  // Fourth order spline derivatives
  SP[0].d = -SP[0].s;
  SP[1].d = SP[0].s - SP[1].s;
  SP[2].d = SP[1].s - SP[2].s;
  SP[3].d = SP[2].s;

  // Fourth order splines
  SP[3].s = w*SP[2].s/3.0;
  SP[2].s = ((w+1.0)*SP[1].s + (3.0-w)*SP[2].s)/3.0;
  SP[0].s = (1.0-w)*SP[0].s/3.0;
  SP[1].s = 1.0-SP[0].s-SP[2].s-SP[3].s;

  // Charge map
  FillQMap(SP, uc-4, ng, 4);
}

//-----------------------------------------------------------------------------
// Extend6: extend third order B-spline coefficients to sixth order and 
//          compute derivatives.                                        
//                                                                      
// Arguments:                                                           
//   See Extend4(), above.                                              
//-----------------------------------------------------------------------------
static void Extend6(bcof* SP, double w, int uc, int ng)
{
  // Fourth order splines
  SP[3].s = w*SP[2].s/3.0;
  SP[2].s = ((w+1.0)*SP[1].s + (3.0-w)*SP[2].s)/3.0;
  SP[0].s = (1.0-w)*SP[0].s/3.0;
  SP[1].s = 1.0-SP[0].s-SP[2].s-SP[3].s;

  // Fifth order splines
  SP[4].s = 0.25*w*SP[3].s;
  SP[3].s = 0.25*((w+1.0)*SP[2].s + (4.0-w)*SP[3].s);
  SP[2].s = 0.25*((w+2.0)*SP[1].s + (3.0-w)*SP[2].s);
  SP[1].s = 0.25*((w+3.0)*SP[0].s + (2.0-w)*SP[1].s);
  SP[0].s = 0.25*(1.0-w)*SP[0].s;

  // Sixth order spline derivatives
  SP[0].d = -SP[0].s;
  SP[1].d = SP[0].s - SP[1].s;
  SP[2].d = SP[1].s - SP[2].s;
  SP[3].d = SP[2].s - SP[3].s;
  SP[4].d = SP[3].s - SP[4].s;
  SP[5].d = SP[4].s;

  // Sixth order splines
  SP[5].s = 0.2*w*SP[4].s;
  SP[4].s = 0.2*((w+1.0)*SP[3].s + (5.0-w)*SP[4].s);
  SP[3].s = 0.2*((w+2.0)*SP[2].s + (4.0-w)*SP[3].s);
  SP[2].s = 0.2*((w+3.0)*SP[1].s + (3.0-w)*SP[2].s);
  SP[1].s = 0.2*((w+4.0)*SP[0].s + (2.0-w)*SP[1].s);
  SP[0].s = 0.2*(1.0-w)*SP[0].s;

  // Charge map
  FillQMap(SP, uc-6, ng, 6);
}

//-----------------------------------------------------------------------------
// ExtendN: extend (N-1)th order B-spline coefficients to Nth order.    
//                                                                      
// Arguments:                                                           
//   SP:    the array of B-spline coefficients and charge mesh indices in this
//          dimension (the dimension of interest)             
//   w:     one minus the non-integer part of the scaled fractional coordinate
//   ordr:  the interpolation order in the dimension of interest        
//-----------------------------------------------------------------------------
static void ExtendN(bcof* SP, double w, int ordr)
{
  int i, om1;
  double dv;

  om1 = ordr-1;
  dv = 1.0/om1;
  SP[om1].s = dv*w*SP[ordr-2].s;
  for (i = 1; i <= ordr-2; i++) {
    SP[om1-i].s = dv*((w+i)*SP[ordr-i-2].s + (ordr-i-w)*SP[om1-i].s);
  }
  SP[0].s = dv*(1.0-w)*SP[0].s;
}

//-----------------------------------------------------------------------------
// DiffN: obtain the derivatives of an N-th order B-Spline from (N-1)th order
//        B-Spline coefficients.                                  
//
// Arguments:                                                           
//   SP:    the array of B-spline coefficients and charge mesh indices in this 
//          dimension (the dimension of interest)             
//   ordr:  the interpolation order in the dimension of interest        
//-----------------------------------------------------------------------------
static void DiffN(bcof* SP, int ordr)
{
  int i;

  SP[0].d = -SP[0].s;
  for (i = 1; i <= ordr-2; i++) {
    SP[i].d = SP[i-1].s - SP[i].s;
  }
  SP[ordr-1].d = SP[ordr-2].s;
}

//-----------------------------------------------------------------------------
// CalcBSpln: this function provides the inner-workings of the B-spline 
//            computation loop used by SplCoeff and CellSplCoeff, making the
//            highly optimized and tedious code accessible to spline
//            computation in different data structures.          
//                                                                      
// Arguments:                                                           
//   u[x,y,z]:   the x, y, or z location of the point in mesh spacings  
//   ng[x,y,z]:  the x, y, or z mesh dimensions                         
//   iplan:      set to 1 if interpolation is (4,4,4), 2 if (6,4,4),    
//               3 if (4,4,6), 4 if (6,6,6), and 0 otherwise            
//   iordr:      interpolation order; used if iplan is 0                
//   tSP[x,y,z]: the b-spline coefficients to fill                      
//-----------------------------------------------------------------------------
static void CalcBSpln(double ux, double uy, double uz, int ngx, int ngy,
		      int ngz, int iplan, int* iordr, bcof* tSPx, bcof* tSPy,
		      bcof* tSPz)
{
  int i;
  double ucx, ucy, ucz, wx, wy, wz;

  // Convert to scaled fractional coordinates
  ucx = ceil(ux);
  ucy = ceil(uy);
  ucz = ceil(uz);

  // One pass to order 3
  wx = ux - ucx + 1.0;
  wy = uy - ucy + 1.0;
  wz = uz - ucz + 1.0;

  // Compute third order splines
  tSPx[2].s = 0.5*wx*wx;
  tSPy[2].s = 0.5*wy*wy;
  tSPz[2].s = 0.5*wz*wz;
  tSPx[0].s = 0.5*(1.0-wx)*(1.0-wx);
  tSPy[0].s = 0.5*(1.0-wy)*(1.0-wy);
  tSPz[0].s = 0.5*(1.0-wz)*(1.0-wz);
  tSPx[1].s = 1.0 - tSPx[0].s - tSPx[2].s;
  tSPy[1].s = 1.0 - tSPy[0].s - tSPy[2].s;
  tSPz[1].s = 1.0 - tSPz[0].s - tSPz[2].s;

  // Special case for isotropic 4th order interpolation
  if (iplan == 1) {
    Extend4(tSPx, wx, ucx, ngx);
    Extend4(tSPy, wy, ucy, ngy);
    Extend4(tSPz, wz, ucz, ngz);
  }

  // Special case for anisotropic 6-4-4 interpolation
  else if (iplan == 2) {
    Extend6(tSPx, wx, ucx, ngx);
    Extend4(tSPy, wy, ucy, ngy);
    Extend4(tSPz, wz, ucz, ngz);
  }

  // Special case for anisotropic 4-4-6 interpolation
  else if (iplan == 3) {
    Extend4(tSPx, wx, ucx, ngx);
    Extend4(tSPy, wy, ucy, ngy);
    Extend6(tSPz, wz, ucz, ngz);
  }

  // Special case for isotropic 6th order interpolation
  else if (iplan == 4) {
    Extend6(tSPx, wx, ucx, ngx);
    Extend6(tSPy, wy, ucy, ngy);
    Extend6(tSPz, wz, ucz, ngz);
  }

  // General case
  else {

    // Things in the X dimension
    if (iordr[0] >= 4) {
      for (i = 4; i < iordr[0]; i++) {
	ExtendN(tSPx, wx, i);
      }
      DiffN(tSPx, iordr[0]);
      ExtendN(tSPx, wx, iordr[0]);
    }
    else {
      tSPx[0].d = wx - 1.0;
      tSPx[1].d = 1.0 - wx - wx;
      tSPx[2].d = wx;
    }
    FillQMap(tSPx, ucx-iordr[0], ngx, iordr[0]);

    // Things in the Y dimension
    if (iordr[1] >= 4) {
      for (i = 4; i < iordr[1]; i++) {
	ExtendN(tSPy, wy, i);
      }
      DiffN(tSPy, iordr[1]);
      ExtendN(tSPy, wy, iordr[1]);
    }
    else {
      tSPy[0].d = wy - 1.0;
      tSPy[1].d = 1.0 - wy - wy;
      tSPy[2].d = wy;
    }
    FillQMap(tSPy, ucy-iordr[1], ngy, iordr[1]);

    // Things in the Z dimension
    if (iordr[2] >= 4) {
      for (i = 4; i < iordr[2]; i++) {
	ExtendN(tSPz, wz, i);
      }
      DiffN(tSPz, iordr[2]);
      ExtendN(tSPz, wz, iordr[2]);
    }
    else {
      tSPz[0].d = wz - 1.0;
      tSPz[1].d = 1.0 - wz - wz;
      tSPz[2].d = wz;
    }
    FillQMap(tSPz, ucz-iordr[2], ngz, iordr[2]);
  }
}

//-----------------------------------------------------------------------------
// SplCoeff: this is an advanced function for optimal computation of B-Spline
//           coefficients.                                     
//                                                                      
// Arguments:                                                           
//   crd:    the system coordinates                                     
//   pmmap:  the particle <--> mesh map                                 
//   rcinp:  the reciprocal space control data                          
//-----------------------------------------------------------------------------
void SplCoeff(coord *crd, bmap *pmmap, reccon *rcinp)
{
  int i, iplan, ngx, ngy, ngz;
  int *tordr;
  double ufacx, ufacy, ufacz;
  double *tloc, *dtmp;

  // Shortcuts
  ngx = rcinp->ng[0];
  ngy = rcinp->ng[1];
  ngz = rcinp->ng[2];
  if (crd->isortho == 1) {
    ufacx = ngx/crd->gdim[0];
    ufacy = ngy/crd->gdim[1];
    ufacz = ngz/crd->gdim[2];
  }
  else {
    dtmp = crd->U.data;
  }
  tloc = crd->loc;
  tordr = rcinp->ordr;
  iplan = 0;
  if (tordr[0] == 4 && tordr[1] == 4 && tordr[2] == 4) {
    iplan = 1;
  }
  else if (tordr[0] == 6 && tordr[1] == 4 && tordr[2] == 4) {
    iplan = 2;
  }
  else if (tordr[0] == 4 && tordr[1] == 4 && tordr[2] == 6) {
    iplan = 3;
  }
  else if (tordr[0] == 6 && tordr[1] == 6 && tordr[2] == 6) {
    iplan = 4;
  }

  // Loop over all atoms
  if (crd->isortho == 1) {
    for (i = 0; i < crd->natom; i++) {

      // Atom identification number
      pmmap->atmid[i] = crd->atmid[i];
      CalcBSpln(ufacx*tloc[3*i], ufacy*tloc[3*i+1], ufacz*tloc[3*i+2], ngx,
		ngy, ngz, iplan, tordr, &pmmap->xcof[tordr[0]*i],
		&pmmap->ycof[tordr[1]*i], &pmmap->zcof[tordr[2]*i]);
    }
  }
  else {
    for (i = 0; i < crd->natom; i++) {

      // Atom identification number
      pmmap->atmid[i] = crd->atmid[i];
      ufacx = dtmp[0]*tloc[3*i] + dtmp[1]*tloc[3*i+1] + dtmp[2]*tloc[3*i+2];
      ufacy = dtmp[3]*tloc[3*i] + dtmp[4]*tloc[3*i+1] + dtmp[5]*tloc[3*i+2];
      ufacz = dtmp[6]*tloc[3*i] + dtmp[7]*tloc[3*i+1] + dtmp[8]*tloc[3*i+2];
      CalcBSpln(ufacx*ngx, ufacy*ngy, ufacz*ngz, ngx, ngy, ngz, iplan, tordr,
                &pmmap->xcof[tordr[0]*i], &pmmap->ycof[tordr[1]*i],
		&pmmap->zcof[tordr[2]*i]);
    }
  }
}

//-----------------------------------------------------------------------------
// CellSplCoeff: this is an advanced function for optimal computation of
//               B-Spline coefficients based on a cell framework.       
//                                                                      
// Arguments:                                                           
//   C:      the cell of interest                                       
//   crd:    the system coordinates (for box dimensions)                
//   rcinp:  the reciprocal space control data                          
//-----------------------------------------------------------------------------
void CellSplCoeff(cell *C, coord *crd, reccon *rcinp)
{
  int i, iplan, ngx, ngy, ngz;
  int *tordr;
  double ufacx, ufacy, ufacz;
  double *dtmp;
  atomc *atmi;

  // Shortcuts
  ngx = rcinp->ng[0];
  ngy = rcinp->ng[1];
  ngz = rcinp->ng[2];
  if (crd->isortho == 1) {
    ufacx = ngx/crd->gdim[0];
    ufacy = ngy/crd->gdim[1];
    ufacz = ngz/crd->gdim[2];
  }
  else {
    dtmp = crd->U.data;
  }
  tordr = rcinp->ordr;
  iplan = 0;
  if (tordr[0] == 4 && tordr[1] == 4 && tordr[2] == 4) {
    iplan = 1;
  }
  else if (tordr[0] == 6 && tordr[1] == 4 && tordr[2] == 4) {
    iplan = 2;
  }
  else if (tordr[0] == 4 && tordr[1] == 4 && tordr[2] == 6) {
    iplan = 3;
  }
  else if (tordr[0] == 6 && tordr[1] == 6 && tordr[2] == 6) {
    iplan = 4;
  }

  // Loop over all atoms in the cell's primary sector
  if (crd->isortho == 1) {
    for (i = 0; i < C->nr[0]; i++) {
      atmi = &C->data[i];
      CalcBSpln(ufacx*atmi->loc[0], ufacy*atmi->loc[1], ufacz*atmi->loc[2],
		ngx, ngy, ngz, iplan, tordr, &C->xcof[tordr[0]*i],
		&C->ycof[tordr[1]*i], &C->zcof[tordr[2]*i]);
    }
  }
  else {
    for (i = 0; i < C->nr[0]; i++) {
      atmi = &C->data[i];
      ufacx = dtmp[0]*atmi->loc[0] + dtmp[1]*atmi->loc[1] +
	dtmp[2]*atmi->loc[2];
      ufacy = dtmp[3]*atmi->loc[0] + dtmp[4]*atmi->loc[1] +
	dtmp[5]*atmi->loc[2];
      ufacz = dtmp[6]*atmi->loc[0] + dtmp[7]*atmi->loc[1] +
	dtmp[8]*atmi->loc[2];
      CalcBSpln(ufacx*ngx, ufacy*ngy, ufacz*ngz, ngx, ngy, ngz, iplan, tordr,
                &C->xcof[tordr[0]*i], &C->ycof[tordr[1]*i],
                &C->zcof[tordr[2]*i]);
    }
  }
}
