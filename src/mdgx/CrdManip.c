#include <math.h>
#include <stdlib.h>
#include "Matrix.h"
#include "pmeRecip.h"
#include "CrdManip.h"
#include "mdgxVector.h"

#include "pmeDirectDS.h"
#include "TopologyDS.h"

//-----------------------------------------------------------------------------
// CreateCoord: create a coordinate structure.                          
//                                                                      
// Arguments:                                                           
//   natom:  the number of atoms to allocate for                        
//-----------------------------------------------------------------------------
coord CreateCoord(int natom)
{
  coord crd;

  crd.natom = natom;
  crd.atmid = (int*)calloc(natom, sizeof(int));
  crd.loc = (double*)calloc(3*natom, sizeof(double));
  crd.prvloc = (double*)calloc(3*natom, sizeof(double));
  crd.scrloc = (double*)calloc(3*natom, sizeof(double));
  crd.vel = (double*)calloc(3*natom, sizeof(double));
  crd.prvvel = (double*)calloc(3*natom, sizeof(double));
  crd.frc = (double*)calloc(3*natom, sizeof(double));
  crd.prvfrc = (double*)calloc(3*natom, sizeof(double));
  crd.scrfrc = (double*)calloc(3*natom, sizeof(double));
  crd.U = CreateDmat(3, 3, 0);
  crd.invU = CreateDmat(3, 3, 0);
  crd.fcnorm = CreateDmat(26, 3, 0);

  return crd;
}

//-----------------------------------------------------------------------------
// CopyCoord: copy coord struct crd into Xcrd.                          
//                                                                      
// Arguments:                                                           
//   crd:  the original coord struct                                    
//-----------------------------------------------------------------------------
coord CopyCoord(coord *crd)
{
  int i;
  coord Xcrd;

  Xcrd.natom  = crd->natom;
  Xcrd.isortho = crd->isortho;
  Xcrd.atmid  = CpyIVec(crd->atmid, crd->natom);
  Xcrd.loc    = CpyDVec(crd->loc, 3*crd->natom);
  Xcrd.prvloc = CpyDVec(crd->prvloc, 3*crd->natom);
  Xcrd.scrloc = CpyDVec(crd->scrloc, 3*crd->natom);
  Xcrd.vel    = CpyDVec(crd->vel, 3*crd->natom);
  Xcrd.prvvel = CpyDVec(crd->prvvel, 3*crd->natom);
  Xcrd.frc    = CpyDVec(crd->frc, 3*crd->natom);
  Xcrd.prvfrc = CpyDVec(crd->prvfrc, 3*crd->natom);
  Xcrd.scrfrc = CpyDVec(crd->scrfrc, 3*crd->natom);
  for (i = 0; i < 6; i++) {
    Xcrd.gdim[i] = crd->gdim[i];
  }
  for (i = 0; i < 3; i++) {
    Xcrd.hgdim[i] = crd->hgdim[i];
  }
  Xcrd.U = CreateDmat(3, 3, 0);
  Xcrd.invU = CreateDmat(3, 3, 0);
  CompXfrm(Xcrd.gdim, Xcrd.U, Xcrd.invU);
  CopyDmat(&Xcrd.fcnorm, &crd->fcnorm, 0);

  return Xcrd;
}

//-----------------------------------------------------------------------------
// DestroyCoord: free all memory associated with a coord struct.        
//                                                                      
// Arguments:                                                           
//   crd: the coord struct to destroy                                   
//-----------------------------------------------------------------------------
void DestroyCoord(coord *crd)
{
  free(crd->atmid);
  free(crd->loc);
  free(crd->prvloc);
  free(crd->scrloc);
  free(crd->vel);
  free(crd->prvvel);
  free(crd->frc);
  free(crd->prvfrc);
  free(crd->scrfrc);
  DestroyDmat(&crd->U);
  DestroyDmat(&crd->invU);
  DestroyDmat(&crd->fcnorm);
}

//-----------------------------------------------------------------------------
// TransCrd: translates a set of coordinates by step*tvec[].            
//
// Arguments:
//   crds:     the coordinates to translate, in the format X1, Y1, Z1, X2,...
//   natom:    the number of atoms
//   tvec:     vector along which to translate
//   step:     translate the coordinates by step*tvec
//-----------------------------------------------------------------------------
void TransCrd(double* crds, int natom, double* tvec, double step)
{
  int i;
  double nvec[3];

  for (i = 0; i < 3; i++) {
    nvec[i] = tvec[i]*step;
  }
  for (i = 0; i < natom; i++) {
    crds[3*i] += nvec[0];
    crds[3*i+1] += nvec[1];
    crds[3*i+2] += nvec[2];
  }
}

//-----------------------------------------------------------------------------
// RotateCrd: rotates a set of coordinates using matrix U.              
//
// Arguments:
//   crds:    the coordinates, in the format X1, Y1, Z1, X2, Y2, Z2, X3, ...
//   natom:   the number of atoms
//   U:       rotation matrix to apply
//-----------------------------------------------------------------------------
void RotateCrd(double* crds, int natom, dmat U)
{
  int i;
  double tmp_crds[3];
  double* dtmp;

  dtmp = &crds[0];
  for (i = 0; i < natom; i++) {
    tmp_crds[0] = U.data[0]*dtmp[0] + U.data[1]*dtmp[1] + U.data[2]*dtmp[2];
    tmp_crds[1] = U.data[3]*dtmp[0] + U.data[4]*dtmp[1] + U.data[5]*dtmp[2];
    tmp_crds[2] = U.data[6]*dtmp[0] + U.data[7]*dtmp[1] + U.data[8]*dtmp[2];
    dtmp[0] = tmp_crds[0];
    dtmp[1] = tmp_crds[1];
    dtmp[2] = tmp_crds[2];
    dtmp = &dtmp[3];
  }
}

//-----------------------------------------------------------------------------
// FindCoordCenter: find the center of a set of coordinates, assuming that the
//                  coordinates are all in the same image.
//
// Arguments:
//   C:      coordinates (double-precision array)                       
//   mass:   masses (double-precision array)                            
//   usem:   flag to activate use of masses                             
//   n:      the number of atoms                                        
//   cofm:   3-element array to store and return center of mass         
//-----------------------------------------------------------------------------
void FindCoordCenter(double* C, double* mass, int usem, int n, double* cofm)
{
  int i;
  double amass, tmass;

  cofm[0] = 0.0;
  cofm[1] = 0.0;
  cofm[2] = 0.0;
  if (usem == 1) {
    tmass = 0.0;
    for (i = 0; i < n; i++) {
      amass = mass[i];
      cofm[0] += amass*C[3*i];
      cofm[1] += amass*C[3*i+1];
      cofm[2] += amass*C[3*i+2];
      tmass += amass;
    }
    tmass = 1.0/tmass;
  }
  else {
    for (i = 0; i < n; i++) {
      cofm[0] += C[3*i];
      cofm[1] += C[3*i+1];
      cofm[2] += C[3*i+2];
    }
    tmass = 1.0/n;
  }
  cofm[0] *= tmass;
  cofm[1] *= tmass;
  cofm[2] *= tmass;
}

//-----------------------------------------------------------------------------
// BeardRotMat: Returns a rotation matrix r given specified rotations   
//              about the x, y, and z lab frame axes.  Here, we use a   
//              bias-free approach as derived in:                       
//                                                                      
//              Daniel A. Beard and Tamar Schlick. "Unbiased Rotational 
//              Moves for Rigid-Body Dynamics." Biophysical Journal     
//              85:2973-2976.  2003.                                    
//
// Arguments:
//   om_{x,y,z}:   rotation about the X, Y, or Z (lab frame) axes
//   r:            rotation matrix (computed within this function) that will
//                 accomplish the unbiased rotational move
//-----------------------------------------------------------------------------
void BeardRotMat(double om_x, double om_y, double om_z, dmat *r)
{
  double om, om2, om2_x, om2_y, om2_z, sin_om, cos_om, inv_om, inv_om2;

  // Convenient quantities
  om2_x = om_x*om_x;
  om2_y = om_y*om_y;
  om2_z = om_z*om_z;
  om2 = om2_x + om2_y + om2_z;
  om = sqrt(om2);
  inv_om = 1.0/om;
  inv_om2 = inv_om*inv_om;
  cos_om = cos(om);
  sin_om = sin(om);

  // Compute rotation matrix
  r->data[0] = ((om2_y + om2_z)*cos_om + om2_x) * inv_om2;
  r->data[1] = (om_x*om_y*inv_om2)*(1.0 - cos_om) + (om_z*inv_om)*sin_om;
  r->data[2] = (om_x*om_z*inv_om2)*(1.0 - cos_om) - (om_y*inv_om)*sin_om;
  r->data[3] = (om_x*om_y*inv_om2)*(1.0 - cos_om) - (om_z*inv_om)*sin_om;
  r->data[4] = ((om2_x + om2_z)*cos_om + om2_y) * inv_om2;
  r->data[5] = (om_y*om_z*inv_om2)*(1.0 - cos_om) + (om_x*inv_om)*sin_om;
  r->data[6] = (om_x*om_z*inv_om2)*(1.0 - cos_om) + (om_y*inv_om)*sin_om;
  r->data[7] = (om_y*om_z*inv_om2)*(1.0 - cos_om) - (om_x*inv_om)*sin_om;
  r->data[8] = ((om2_x + om2_y)*cos_om + om2_z) * inv_om2;
}

//-----------------------------------------------------------------------------
// CompXfrm: compute the transformation matrices U and invU that take a set of
//           coordinates into and out of box space.  The matrices U and invU
//           must be pre-allocated.
//
// Arguments:
//   cd:   the box dimensions
//   U:    transformation matrix INTO the box vector space
//   invU: transformation matrix OUT OF box vector space back into real space
//-----------------------------------------------------------------------------
void CompXfrm(double* cd, dmat U, dmat invU)
{
  double dx, dy;

  dx = (cos(cd[4])*cos(cd[5]) - cos(cd[3])) / (sin(cd[4]) * sin(cd[5]));
  dy = sqrt(1.0 - dx*dx);
  invU.data[0] = cd[0];
  invU.data[1] = cd[1]*cos(cd[5]);
  invU.data[2] = cd[2]*cos(cd[4]);
  invU.data[4] = cd[1]*sin(cd[5]);
  invU.data[5] = -cd[2]*sin(cd[4])*dx;
  invU.data[8] = cd[2]*sin(cd[4])*dy;
  ttInv(invU, U);
}

//-----------------------------------------------------------------------------
// OrthoReim: re-image a single set of three coordinates (such as a 
//            displacement) given an orthorhombic box.  No assumptions are made
//            are made about the initial displacements.
//
// Arguments:
//   dx:                                                                
//   dy:     Cartesian x, y, and z coordinates                          
//   dz:                                                                
//   [inv]U: transformation matrices [into] and out of unit cell space  
//-----------------------------------------------------------------------------
void OrthoReim(double *dx, double *dy, double *dz, dmat *U, dmat *invU)
{
  double ndx, ndy, ndz;
  double *dtmp;

  dtmp = U->data;
  ndx = dtmp[0]*(*dx);
  ndy = dtmp[4]*(*dy);
  ndz = dtmp[8]*(*dz);
  if (ndx >= 0.5) {
    ndx -= floor(ndx+0.5);
  }
  else if (ndx < -0.5) {
    ndx -= ceil(ndx-0.5);
  }
  if (ndy >= 0.5) {
    ndy -= floor(ndy+0.5);
  }
  else if (ndy < -0.5) {
    ndy -= ceil(ndy-0.5);
  }
  if (ndz >= 0.5) {
    ndz -= floor(ndz+0.5);
  }
  else if (ndz < -0.5) {
    ndz -= ceil(ndz-0.5);
  }
  dtmp = invU->data;
  *dx = dtmp[0]*ndx;
  *dy = dtmp[4]*ndy;
  *dz = dtmp[8]*ndz;
}

//-----------------------------------------------------------------------------
// NonOrthoReim: re-image a single set of three coordinates (such as a  
//               displacement) given an non-orthorhombic box and the
//               appropriate transformation matrices.  No assumptions are made
//               about the initial distances between coordinates.
//
// Arguments:
//   dx:                                                                
//   dy:     Cartesian x, y, and z coordinates                          
//   dz:                                                                
//   [inv]U: transformation matrices [into] and out of unit cell space  
//-----------------------------------------------------------------------------
void NonOrthoReim(double *dx, double *dy, double *dz, dmat *U, dmat *invU)
{
  double ndx, ndy, ndz;
  double *dtmp;

  dtmp = U->data;
  ndx = dtmp[0]*(*dx) + dtmp[1]*(*dy) + dtmp[2]*(*dz);
  ndy = dtmp[3]*(*dx) + dtmp[4]*(*dy) + dtmp[5]*(*dz);
  ndz = dtmp[6]*(*dx) + dtmp[7]*(*dy) + dtmp[8]*(*dz);
  if (ndx >= 0.5) {
    ndx -= floor(ndx+0.5);
  }
  else if (ndx < -0.5) {
    ndx -= ceil(ndx-0.5);
  }
  if (ndy >= 0.5) {
    ndy -= floor(ndy+0.5);
  }
  else if (ndy < -0.5) {
    ndy -= ceil(ndy-0.5);
  }
  if (ndz >= 0.5) {
    ndz -= floor(ndz+0.5);
  }
  else if (ndz < -0.5) {
    ndz -= ceil(ndz-0.5);
  }
  dtmp = invU->data;
  *dx = dtmp[0]*ndx + dtmp[1]*ndy + dtmp[2]*ndz;
  *dy = dtmp[3]*ndx + dtmp[4]*ndy + dtmp[5]*ndz;
  *dz = dtmp[6]*ndx + dtmp[7]*ndy + dtmp[8]*ndz;
}

//-----------------------------------------------------------------------------
// ImageBondedGroups: align atoms in each bonded group to be in the same
//                    periodic image.  This routine loops over all atoms in the
//                    group, starting with the first, and re-images other atoms
//                    to be in the minimum image. 
//
// Arguments:
//   crd:      the coordinates (this routine operates on a unified list of
//             coordinates, not a cell grid)                         
//   tp:       the topology (contains bonded lists of atoms)            
//-----------------------------------------------------------------------------
void ImageBondedGroups(coord *crd, prmtop *tp)
{
  int h, i, tgai3;
  double cenx, ceny, cenz, dx, dy, dz;
  lgrp *tg;

  // Take all atoms into box space
  RotateCrd(crd->loc, crd->natom, crd->U);

  // Loop over all groups
  for (h = 0; h < tp->ngrp; h++) {
    tg = &tp->lgrps[h];
    cenx = crd->loc[3*tg->atoms[0]];
    ceny = crd->loc[3*tg->atoms[0]+1];
    cenz = crd->loc[3*tg->atoms[0]+2];
    for (i = 1; i < tg->natom; i++) {
      tgai3 = 3*tg->atoms[i];
      dx = crd->loc[tgai3] - cenx;
      dy = crd->loc[tgai3+1] - ceny;
      dz = crd->loc[tgai3+2] - cenz;
      dx -= floor(dx+0.5);
      dy -= floor(dy+0.5);
      dz -= floor(dz+0.5);
      crd->loc[tgai3] = cenx + dx;
      crd->loc[tgai3+1] = ceny + dy;
      crd->loc[tgai3+2] = cenz + dz;
    }
  }

  // Take all atoms back to real space
  RotateCrd(crd->loc, crd->natom, crd->invU);
}

//-----------------------------------------------------------------------------
// QuatAlign: function for aligning sets of coordinates using a quaternion.
//            Assumes that the sets have been aligned in terms of geometric
//            center.  Returns as U the rotation which minimizes the RMSD
//            between the sets of coordinates.
//
// Arguments:
//   frameI:    the frame to align against frame II
//   frameII:   the reference frame
//   natom:     the number of atoms
//   mass:      optional vector of atomic / particle masses
//   m:         flag to have mass-weighted alignment used
//   U:         rotation matrix that will store the rotation needed to take
//              frameI into frameII
//-----------------------------------------------------------------------------
void QuatAlign(double* frameI, double* frameII, int natom, double* mass,
               int m, dmat *U)
{
  int i, itr, info;
  double maxeig, a, x, y, z, totmass, tmass;
  double aa, ab, ac, ba, bb, bc, ca, cb, cc, i1, i2, i3, g1, g2, g3;
  double diag[4];
  dmat R, V;

  // Compute the quaternion matrix
  R = CreateDmat(4, 4, 0);
  totmass = (m == 1) ? 1.0/DSum(mass, natom) : 1.0/natom;
  for (i = 0; i < natom; i++) {
    itr = 3*i;
    i1 = frameII[itr];
    i2 = frameII[itr+1];
    i3 = frameII[itr+2];
    g1 = frameI[itr];
    g2 = frameI[itr+1];
    g3 = frameI[itr+2];
    aa = i1*g1;
    ab = i1*g2;
    ac = i1*g3;
    ba = i2*g1;
    bb = i2*g2;
    bc = i2*g3;
    ca = i3*g1;
    cb = i3*g2;
    cc = i3*g3;
    if (m == 1) {
      tmass = mass[i];
      R.data[0] += tmass*(aa+bb+cc);
      R.data[1] += tmass*(cb-bc);
      R.data[2] += tmass*(ac-ca);
      R.data[3] += tmass*(ba-ab);
      R.data[5] += tmass*(aa-bb-cc);
      R.data[6] += tmass*(ab+ba);
      R.data[7] += tmass*(ca+ac);
      R.data[10] += tmass*(bb-cc-aa);
      R.data[11] += tmass*(bc+cb);
      R.data[15] += tmass*(cc-aa-bb);
    }
    else {
      R.data[0] += aa+bb+cc;
      R.data[1] += cb-bc;
      R.data[2] += ac-ca;
      R.data[3] += ba-ab;
      R.data[5] += aa-bb-cc;
      R.data[6] += ab+ba;
      R.data[7] += ca+ac;
      R.data[10] += bb-cc-aa;
      R.data[11] += bc+cb;
      R.data[15] += cc-aa-bb;
    }
  }
  R.data[4] = R.data[1];
  R.data[8] = R.data[2];
  R.data[12] = R.data[3];
  R.data[9] = R.data[6];
  R.data[13] = R.data[7];
  R.data[14] = R.data[11];
  for (i = 0; i < 16; i++) {
    R.data[i] *= totmass;
  }
  V = CreateDmat(4, 4, 0);
  info = jacobi(R.map, 4, diag, V.map, &i);
  eigsrt(diag, V.map, 4);
  a = V.data[0];
  x = V.data[4];
  y = V.data[8];
  z = V.data[12];
  
  // Construct the rotation matrix
  U->data[0] = a*a + x*x -y*y - z*z;
  U->data[1] = 2.0*(x*y + a*z);
  U->data[2] = 2.0*(z*x - a*y);
  U->data[3] = 2.0*(x*y - a*z);
  U->data[4] = a*a - x*x + y*y - z*z;
  U->data[5] = 2.0*(y*z + a*x);
  U->data[6] = 2.0*(z*x + a*y);
  U->data[7] = 2.0*(y*z - a*x);
  U->data[8] = a*a - x*x - y*y + z*z;

  DestroyDmat(&R);
}
