#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#define PI                     3.1415926535898

/***=======================================================================***/
/*** TransCrd: translates a set of coordinates by step*tvec[].             ***/
/***=======================================================================***/
void TransCrd(double* crds, int n_atoms, double* tvec, double step)
{
  int i;
  double nvec[3];

  for (i = 0; i < 3; i++) {
    nvec[i] = tvec[i]*step;
  }
  for (i = 0; i < n_atoms; i++) {
    crds[3*i] += nvec[0];
    crds[3*i+1] += nvec[1];
    crds[3*i+2] += nvec[2];
  }
}

/***=======================================================================***/
/*** RotateCrd: rotates a set of coordinates using matrix U.               ***/
/***=======================================================================***/
void RotateCrd(double* crds, int n_atoms, dmat U)
{
  int i;
  double tmp_crds[3];
  double* dtmp;

  dtmp = &crds[0];
  for (i = 0; i < n_atoms; i++) {
    tmp_crds[0] = U.data[0]*dtmp[0] + U.data[1]*dtmp[1] + U.data[2]*dtmp[2];
    tmp_crds[1] = U.data[3]*dtmp[0] + U.data[4]*dtmp[1] + U.data[5]*dtmp[2];
    tmp_crds[2] = U.data[6]*dtmp[0] + U.data[7]*dtmp[1] + U.data[8]*dtmp[2];
    dtmp[0] = tmp_crds[0];
    dtmp[1] = tmp_crds[1];
    dtmp[2] = tmp_crds[2];
    dtmp = &dtmp[3];
  }
}

/***=======================================================================***/
/*** BeardRotMat: Returns a rotation matrix r given specified rotations    ***/
/***              about the x, y, and z lab frame axes.  Here, we use a    ***/
/***              bias-free approach as derived in:                        ***/
/***                                                                       ***/
/***              Daniel A. Beard and Tamar Schlick. "Unbiased Rotational  ***/
/***              Moves for Rigid-Body Dynamics." Biophysical Journal      ***/
/***              85:2973-2976.  2003.                                     ***/
/***=======================================================================***/
void BeardRotMat(double om_x, double om_y, double om_z, dmat r)
{
  double om, om2, om2_x, om2_y, om2_z, sin_om, cos_om, inv_om, inv_om2;

  /*** Convenient quantities ***/
  om2_x = om_x*om_x;
  om2_y = om_y*om_y;
  om2_z = om_z*om_z;
  om2 = om2_x + om2_y + om2_z;
  om = sqrt(om2);
  inv_om = 1.0/om;
  inv_om2 = inv_om*inv_om;
  cos_om = cos(om);
  sin_om = sin(om);

  /*** Compute rotation matrix ***/
  r.data[0] = ((om2_y + om2_z)*cos_om + om2_x) * inv_om2;
  r.data[1] = (om_x*om_y*inv_om2)*(1.0 - cos_om) + (om_z*inv_om)*sin_om;
  r.data[2] = (om_x*om_z*inv_om2)*(1.0 - cos_om) - (om_y*inv_om)*sin_om;
  r.data[3] = (om_x*om_y*inv_om2)*(1.0 - cos_om) - (om_z*inv_om)*sin_om;
  r.data[4] = ((om2_x + om2_z)*cos_om + om2_y) * inv_om2;
  r.data[5] = (om_y*om_z*inv_om2)*(1.0 - cos_om) + (om_x*inv_om)*sin_om;
  r.data[6] = (om_x*om_z*inv_om2)*(1.0 - cos_om) + (om_y*inv_om)*sin_om;
  r.data[7] = (om_y*om_z*inv_om2)*(1.0 - cos_om) - (om_x*inv_om)*sin_om;
  r.data[8] = ((om2_x + om2_y)*cos_om + om2_z) * inv_om2;
}

/***=======================================================================***/
/*** ExtCrds: finds the extreme coordinates of a given vector, assuming    ***/
/***          it to be laid out as [x1 y1 z1, x2 y2 z2, ... xN yN zN]      ***/
/***=======================================================================***/
void ExtCrds(double* crds, int num_atoms, double* ext)
{
  int i, j;

  for (j = 0; j < 3; j++) {
    ext[j] = crds[j];
    ext[j+3] = crds[j];
  }
  for (i = 1; i < num_atoms; i++) {
    for (j = 0; j < 3; j++) {
      if (crds[3*i+j] < ext[j]) {
        ext[j] = crds[3*i+j];
      }
      else if (crds[3*i+j] > ext[j+3]) {
        ext[j+3] = crds[3*i+j];
      }
    }
  }
}

/***=======================================================================***/
/*** FindCrdCenter: finds the center of a set of coordinates.  Flag m == 1 ***/
/***                activates mass weighting.                              ***/
/***=======================================================================***/
void FindCrdCenter(double* crds, double* mass, int m, int n_atoms, double* cm)
{
  int i;
  double tmp_mass, sum_mass;

  sum_mass = 0;
  for (i = 0; i < 3; i++) {
    cm[i] = 0.0;
  }
  for (i = 0; i < n_atoms; i++) {
    tmp_mass = (m == 1) ? mass[i] : 1.0;
    cm[0] += tmp_mass*crds[3*i];
    cm[1] += tmp_mass*crds[3*i+1];
    cm[2] += tmp_mass*crds[3*i+2];
    sum_mass += tmp_mass;
  }
  sum_mass = 1.0/sum_mass;
  for (i = 0; i < 3; i++) {
    cm[i] *= sum_mass;
  }
}

/***=======================================================================***/
/*** CmpXfrm: compute the transformation matrices U and invU that take a   ***/
/***          set of coordinates into and out of box space.  The matrices  ***/
/***          U and invU must be pre-allocated.                            ***/
/***=======================================================================***/
void CmpXfrm(double* cd, dmat U, dmat invU)
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

/***=======================================================================***/
/*** ReImage: re-image n atoms into the primary cell.  This assumes that   ***/
/***          all coordinates have been translated into the crystal space. ***/
/***=======================================================================***/
void ReImage(double* crds, int n)
{
  int i, j;
  double *tmp_crd;

  for (i = 0; i < n; i++) {
    tmp_crd = &crds[3*i];
    for (j = 0; j < 3; j++) {
      if (tmp_crd[j] > 0.5) {
	tmp_crd[j] -= floor(tmp_crd[j]+0.5);
      }
      else if (tmp_crd[j] < -0.5) {
	tmp_crd[j] -= ceil(tmp_crd[j]-0.5);
      }
    }
  }
}

/***=======================================================================***/
/*** UndoOct: undo the ghastly things that AMBER does to octahedral boxes  ***/
/***          so that I can work with them more simply. This returns the   ***/
/***          necessary rotation matrix.                                   ***/
/***=======================================================================***/
dmat UndoOct()
{
  double cos1, sin1, cos2, sin2;
  dmat R;

  R = CreateDmat(3, 3);
  cos1 = cos(PI/4.0);
  sin1 = sin(PI/4.0);
  cos2 = sqrt(2.0)/sqrt(3.0);
  sin2 = 1.0/sqrt(3.0);
  R.data[0] = cos1*cos2;
  R.data[1] = -cos1*sin2;
  R.data[2] = sin1;
  R.data[3] = -sin1*cos2;
  R.data[4] = sin2*sin1;
  R.data[5] = cos1;
  R.data[6] = -sin2;
  R.data[7] = -cos2;
  R.data[8] = 0.0;

  return R;
}

/***=======================================================================***/
/*** SimpleReIm: re-image a coordinate assuming its within [ -1.0, 1.0].   ***/
/***=======================================================================***/
double SimpleReIm(double dx)
{
  if (dx > 0.5) {
    return dx - 1.0;
  }
  else if (dx < -0.5) {
    return dx + 1.0;
  }
  return dx;
}

/***=======================================================================***/
/*** RenderFrame: take a frame from .crd format into a double precision    ***/
/***              array.                                                   ***/
/***=======================================================================***/
int RenderFrame(double* crd, double* boxd, int natm, FILE *inp)
{
  int i, la, lc, tri_atm, ncon, lnpos;
  double npsw, runsum;
  double decnum[6], deckey[6];
  char tmpc;
  char line[128];

  /*** The decimal key ***/
  deckey[0] = 100000.0;
  deckey[1] = 10000.0;
  deckey[2] = 1000.0;
  deckey[3] = 100.0;
  deckey[4] = 10.0;
  deckey[5] = 1.0;

  /*** Obtain coordinates ***/
  tri_atm = natm*3;
  la = 0;
  npsw = 0.;
  for (i = 0; i < 128; i++) {
    line[i] = ' ';
  }
  while (la < tri_atm && fgets(line, 128, inp) != NULL) {

    lc = 0;
    lnpos = 0;
    while (lc < 10 && la < tri_atm) {

      /*** Obtain a number ***/
      for (i = 0; i < 6 ;i++) {
        decnum[i] = 0.0;
      }
      ncon = -1;
      while (line[lnpos] != '.') {
        tmpc = line[lnpos];
        if (tmpc > ' ') {
          if (ncon == -1) {
            if (tmpc == '-') {
              npsw = -1.0;
              ncon = 0;
            } else {
              npsw = 1.0;
              decnum[0] = tmpc - 48;
              ncon = 1;
            }
          }
          else {
            decnum[ncon++] = tmpc - 48;
          }
        }
        lnpos++;
      }

      /*** Compute the whole number ***/
      runsum = 0.0;
      for (i = 0; i < ncon; i++) {
        runsum += decnum[i]*deckey[6-ncon+i];
      }

      /*** Compute the decimal component ***/
      lnpos++;
      runsum += 0.1*(line[lnpos++] - 48);
      runsum += 0.01*(line[lnpos++] - 48);
      runsum += 0.001*(line[lnpos++] - 48);

      /*** Store the coordinates ***/
      crd[la] = runsum*npsw;
      lc++;
      la++;
    }
  }

  /*** Read the box information ***/
  if (la == tri_atm) {
    if (fgets(line, 128, inp) == NULL) {
      la = 0;
    }
    else {
      sscanf(line, "%lf%lf%lf", &boxd[0], &boxd[1], &boxd[2]);
      for (i = 3; i < 6; i++) {
        boxd[i] = 0.5*PI;
      }
    }
  }

  /*** If the frame is complete, then process it ***/
  return (la == tri_atm) ? 1 : 0;
}
