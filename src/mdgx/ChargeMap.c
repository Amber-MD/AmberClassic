#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Grid.h"
#include "pmeRecip.h"
#include "mdgxVector.h"
#include "Constants.h"
#include "ChargeMap.h"
#include "Matrix.h"
#include "BSpline.h"

#include "BSplineDS.h"
#include "TopologyDS.h"
#include "CellManipDS.h"

//-----------------------------------------------------------------------------
// CellQ2Book: map a set of atoms from a cell's primary sector to a 
//             double-precision book for later use in various Particle Mesh
//             Ewald calculations.  Note that this routine will zero the book  
//             before mapping any charges.    
//                                                                      
// Arguments:                                                           
//   C:      the cell of interest                                       
//   rcinp:  the reciprocal space input parameters                      
//   Q:      the charge mesh                                            
//-----------------------------------------------------------------------------
void CellQ2Book(cell *C, reccon *rcinp, dbook *Q)
{
  int h, i, j, k, ordrX, ordrY ,ordrZ;
  double spx, spxy, atmq;
  double *dtmp;
  double **dtm2p;
  bcof *xcf, *ycf, *zcf;

  // Shortcuts
  ordrX = rcinp->ordr[0];
  ordrY = rcinp->ordr[1];
  ordrZ = rcinp->ordr[2];

  // Map all the atoms
  if (ordrZ == 4) {
    for (h = 0; h < C->nr[0]; h++) {
      xcf = &C->xcof[ordrX*h];
      ycf = &C->ycof[ordrY*h];
      zcf = &C->zcof[ordrZ*h];
      atmq = C->data[h].q;
      for (i = 0; i < ordrX; i++) {
        dtm2p = Q->map[xcf[i].m];
        spx = atmq*xcf[i].s;
        for (j = 0; j < ordrY; j++) {
	  dtmp = dtm2p[ycf[j].m];
	  spxy = spx*ycf[j].s;
	  dtmp[zcf[0].m] += spxy*zcf[0].s;
	  dtmp[zcf[1].m] += spxy*zcf[1].s;
	  dtmp[zcf[2].m] += spxy*zcf[2].s;
	  dtmp[zcf[3].m] += spxy*zcf[3].s;
	}
      }
    }
  }
  else if (ordrZ == 6) {
    for (h = 0; h < C->nr[0]; h++) {
      xcf = &C->xcof[ordrX*h];
      ycf = &C->ycof[ordrY*h];
      zcf = &C->zcof[ordrZ*h];
      atmq = C->data[h].q;
      for (i = 0; i < ordrX; i++) {
        dtm2p = Q->map[xcf[i].m];
        spx = atmq*xcf[i].s;
        for (j = 0; j < ordrY; j++) {
	  dtmp = dtm2p[ycf[j].m];
	  spxy = spx*ycf[j].s;
	  dtmp[zcf[0].m] += spxy*zcf[0].s;
	  dtmp[zcf[1].m] += spxy*zcf[1].s;
	  dtmp[zcf[2].m] += spxy*zcf[2].s;
	  dtmp[zcf[3].m] += spxy*zcf[3].s;
	  dtmp[zcf[4].m] += spxy*zcf[4].s;
	  dtmp[zcf[5].m] += spxy*zcf[5].s;
	}
      }
    }
  }
  else {
    for (h = 0; h < C->nr[0]; h++) {
      xcf = &C->xcof[ordrX*h];
      ycf = &C->ycof[ordrY*h];
      zcf = &C->zcof[ordrZ*h];
      atmq = C->data[h].q;
      for (i = 0; i < ordrX; i++) {
        dtm2p = Q->map[xcf[i].m];
        spx = atmq*xcf[i].s;
        for (j = 0; j < ordrY; j++) {
          dtmp = dtm2p[ycf[j].m];
          spxy = spx*ycf[j].s;
	  for (k = 0; k < ordrZ; k++) {
            dtmp[zcf[k].m] += spxy*zcf[k].s;
	  }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// CellIntrpFrc: interpolate forces from the mesh to atoms within cells.  This
//               routine will not zero forces; it continues to accumulate
//               forces on each atom.
//
// Arguments:
//   C:      the cell of interest                                       
//   crd:    contains simulation box dimensions                         
//   rcinp:  the reciprocal space input parameters                      
//   U:      the potential mesh                                         
//-----------------------------------------------------------------------------
void CellIntrpFrc(cell *C, coord *crd, reccon *rcinp, dbook *U)
{
  int h, i, j, k, natom, ordrX, ordrY, ordrZ, isortho;
  double SPx, SPy, SPz, SPdx, SPxy, SPdxy, SPxdy, fx, fy, fz, uk, atmq;
  double invgx, invgy, invgz, ngx, ngy, ngz, Zuk;
  double *utmp, *umd;
  double **utm2p;
  bcof *xcf, *ycf, *zcf;
  atomc *atmh;

  // Shortcuts
  natom = C->nr[0];
  ordrX = rcinp->ordr[0];
  ordrY = rcinp->ordr[1];
  ordrZ = rcinp->ordr[2];
  ngx = rcinp->ng[0];
  ngy = rcinp->ng[1];
  ngz = rcinp->ng[2];
  invgx = 1.0/crd->gdim[0];
  invgy = 1.0/crd->gdim[1];
  invgz = 1.0/crd->gdim[2];
  umd = crd->U.data;
  isortho = crd->isortho;

  // Loop over all atoms
  for (h = 0; h < natom; h++) {
    atmh = &C->data[h];
    xcf = &C->xcof[ordrX*h];
    ycf = &C->ycof[ordrY*h];
    zcf = &C->zcof[ordrZ*h];
    atmq = atmh->q;
    fx = 0.0;
    fy = 0.0;
    fz = 0.0;
    for (i = 0; i < ordrX; i++) {
      utm2p = U->map[xcf[i].m];
      SPx = atmq*xcf[i].s;
      SPdx = atmq*xcf[i].d;

      // Branches for 4, 6, and nth-order splines
      if (ordrZ == 4) {
	for (j = 0; j < ordrY; j++) {
	  utmp = utm2p[ycf[j].m];
	  SPy = ycf[j].s;
	  SPxy = SPx*SPy;
	  SPdxy = SPdx*SPy;
	  SPxdy = SPx*ycf[j].d;
	  Zuk = zcf[0].s*utmp[zcf[0].m] + zcf[1].s*utmp[zcf[1].m] +
	        zcf[2].s*utmp[zcf[2].m] + zcf[3].s*utmp[zcf[3].m];
	  fx -= SPdxy*Zuk;
	  fy -= SPxdy*Zuk;
	  fz -= SPxy*(zcf[0].d*utmp[zcf[0].m] + zcf[1].d*utmp[zcf[1].m] +
		      zcf[2].d*utmp[zcf[2].m] + zcf[3].d*utmp[zcf[3].m]);
	}
      }
      else if (ordrZ == 6) {
        for (j = 0; j < ordrY; j++) {
          utmp = utm2p[ycf[j].m];
          SPy = ycf[j].s;
          SPxy = SPx*SPy;
          SPdxy = SPdx*SPy;
          SPxdy = SPx*ycf[j].d;
          Zuk = zcf[0].s*utmp[zcf[0].m] + zcf[1].s*utmp[zcf[1].m] +
	        zcf[2].s*utmp[zcf[2].m] + zcf[3].s*utmp[zcf[3].m] +
	        zcf[4].s*utmp[zcf[4].m] + zcf[5].s*utmp[zcf[5].m];
          fx -= SPdxy*Zuk;
          fy -= SPxdy*Zuk;
          fz -= SPxy*(zcf[0].d*utmp[zcf[0].m] + zcf[1].d*utmp[zcf[1].m] +
                      zcf[2].d*utmp[zcf[2].m] + zcf[3].d*utmp[zcf[3].m] +
                      zcf[4].d*utmp[zcf[4].m] + zcf[5].d*utmp[zcf[5].m]);
        }
      }
      else {
	for (j = 0; j < ordrY; j++) {
	  utmp = utm2p[ycf[j].m];
	  SPy = ycf[j].s;
	  SPxy = SPx*SPy;
	  SPdxy = SPdx*SPy;
	  SPxdy = SPx*ycf[j].d;
	  for (k = 0; k < ordrZ; k++) {
	    uk = utmp[zcf[k].m];
	    SPz = zcf[k].s;
	    fx -= SPdxy*SPz*uk;
	    fy -= SPxdy*SPz*uk;
	    fz -= SPxy*zcf[k].d*uk;
	  }
	}
      }
    }

    // The total force must be normalized within the box.
    // Multiplying by the number of grid points and then 
    // the inverse grid dimension or the transformation  
    // into unit cell space (both in units of reciprocal 
    // Angstroms) produces forces in kcal/mol-A.          
    fx *= ngx;
    fy *= ngy;
    fz *= ngz;
    if (isortho == 1) {
      atmh->frc[0] += fx*invgx;
      atmh->frc[1] += fy*invgy;
      atmh->frc[2] += fz*invgz;
    }
    else {
      atmh->frc[0] += umd[0]*fx + umd[3]*fy + umd[6]*fz;
      atmh->frc[1] += umd[1]*fx + umd[4]*fy + umd[7]*fz;
      atmh->frc[2] += umd[2]*fx + umd[5]*fy + umd[8]*fz;
    }
  }
}
