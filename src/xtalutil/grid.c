#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "grid.h"

/***=======================================================================***/
/*** CheckGridDims: check the requested dimensions of a grid to make sure  ***/
/***                that they look all right.                              ***/
/***=======================================================================***/
void CheckGridDims(int lx, int ly, int lz, double ox, double oy, double oz,
		   double mx, double my, double mz)
{
  /*** No negative values in the cell counts ***/
  if (lx <= 0) {
    printf("CheckGridDims >> Error.  lx = %d.\n", lx);
    exit(1);
  }
  if (ly <= 0) {
    printf("CheckGridDims >> Error.  ly = %d.\n", lx);
    exit(1);
  }
  if (lz <= 0) {
    printf("CheckGridDims >> Error.  lz = %d.\n", lx);
    exit(1);
  }
  if (mx - ox <= 0.0) {
    printf("CheckGridDims >> Error.  mx = %15.8e, ox = %15.8e.\n", mx, ox);
    exit(1);
  }
  if (my - oy <= 0.0) {
    printf("CheckGridDims >> Error.  my = %15.8e, oy = %15.8e.\n", my, oy);
    exit(1);
  }
  if (mz - oz <= 0.0) {
    printf("CheckGridDims >> Error.  mz = %15.8e, oz = %15.8e.\n", mz, oz);
    exit(1);
  }
}

/***=======================================================================***/
/*** CreateBoolGrid: create a boolean grid.  These grids are very compact  ***/
/***                 but are not as easy to access as grids of elemental   ***/
/***                 types.  Use WriteBoolGrid(i,j,k) to write to index    ***/
/***                 (i,j,k) of a boolean grid, and ReadBoolGrid(i,j,k) to ***/
/***                 read that element.                                    ***/
/***=======================================================================***/
bgrid CreateBoolGrid(int lx, int ly, int lz, double ox, double oy, double oz,
		     double mx, double my, double mz)
{
  size_t tcsize;
  bgrid tg;

  CheckGridDims(lx, ly, lz, ox, oy, oz, mx, my, mz);
  tg.lx = lx;
  tg.ly = ly;
  tg.lz = lz;
  tg.total_data = lx*ly*lz;
  tg.gx = (mx-ox)/lx;
  tg.gy = (my-oy)/ly;
  tg.gz = (mz-oz)/lz;
  tg.invgx = 1.0/tg.gx;
  tg.invgy = 1.0/tg.gy;
  tg.invgz = 1.0/tg.gz;
  tg.ox = ox;
  tg.oy = oy;
  tg.oz = oz;
  tg.mx = mx;
  tg.my = my;
  tg.mz = mz;

  /*** This grid is allocated by numbers of bits ***/
  tcsize = sizeof(unsigned char);
  tg.cdatasize = tg.lx*tg.ly*tg.lz/(tcsize*8)+1;
  tg.data = (unsigned char*)calloc(tg.cdatasize, tcsize);

  return tg;
}

/***=======================================================================***/
/*** WriteBoolGrid: write a value V to element i,j,k of boolean grid tg.   ***/
/***                A value of 0 is written iff V == 0, a value of 1 is    ***/
/***                written otherwise.                                     ***/
/***=======================================================================***/
void WriteBoolGrid(bgrid tg, int i, int j, int k, int V)
{
  int nbit;
  long long int idx, nbyte;

  idx = (i*tg.ly + j)*tg.lz + k;
  nbit = idx & 0x07;
  nbyte = idx >> 3;

  if (V == 0) {
    tg.data[nbyte] |= (0 << nbit);
  }
  else {
    tg.data[nbyte] |= (1 << nbit);
  }
}

/***=======================================================================***/
/*** ReadBoolGrid: read the value of bit (i,j,k) in boolean grid V, return ***/
/***               0 or 1.                                                 ***/
/***=======================================================================***/
int ReadBoolGrid(bgrid tg, int i, int j, int k)
{
  int nbit;
  long long int idx, nbyte;

  idx = (i*tg.ly + j)*tg.lz + k;
  nbit = idx & 0x07;
  nbyte = idx >> 3;
  if (tg.data[nbyte] & (1 << nbit)) {
    return 1;
  }
  return 0;
}

/***=======================================================================***/
/*** BiColorGrid: color a binary grid.                                     ***/
/***=======================================================================***/
int BiColorGrid(bgrid tg, double* crds, int natm, dmat U, dmat invU,
		double R, int updateUser)
{
  int h, h3, i, j, k, ibuff, jbuff, kbuff, ncol;
  int imin, jmin, kmin, imax, jmax, kmax;
  double ax, ay, az, dx, dy, dz, R2;
  double xyz0[3], ic1[3], ic2[3], ic3[3], thx[3], thy[3], thz[3];
  double px, py, pz;

  /*** Grid coloring counter ***/
  ncol = 0;

  /*** Determine limits ***/
  R2 = R*R;
  xyz0[0] = invU.data[0] + invU.data[1] + invU.data[2];
  xyz0[1] = invU.data[4] + invU.data[5];
  xyz0[2] = invU.data[8];
  for (i = 0; i < 3; i++) {
    ic1[i] = invU.map[i][0];
    ic2[i] = invU.map[i][1];
    ic3[i] = invU.map[i][2];
  }
  CrossP(ic2, ic3, thx);
  CrossP(ic1, ic3, thy);
  CrossP(ic1, ic2, thz);
  Normalize(thx, 3);
  Normalize(thy, 3);
  Normalize(thz, 3);
  ibuff = ceil(tg.lx/(fabs(DotP(thx, xyz0, 3))/R))+1.01;
  jbuff = ceil(tg.ly/(fabs(DotP(thy, xyz0, 3))/R))+1.01;
  kbuff = ceil(tg.lz/(fabs(DotP(thz, xyz0, 3))/R))+1.01;

  for (h = 0; h < natm; h++) {

    if (h % 64 == 0 && updateUser == 1) {
      fprintf(stderr, "\rBiColorGrid >> Coloring atom %6d", h);
      fflush(stderr);
    }

    h3 = 3*h;
    ax = crds[h3] - tg.ox;
    ay = crds[h3+1] - tg.oy;
    az = crds[h3+2] - tg.oz;
    i = ax / tg.gx;
    j = ay / tg.gy;
    k = az / tg.gz;
    imin = (i > ibuff) ? i - ibuff : 0;
    imax = (i + ibuff < tg.lx) ? i + ibuff : tg.lx;
    jmin = (j > jbuff) ? j - jbuff : 0;
    jmax = (j + jbuff < tg.ly) ? j + jbuff : tg.ly;
    kmin = (k > kbuff) ? k - kbuff : 0;
    kmax = (k + kbuff < tg.lz) ? k + kbuff : tg.lz;
    for (i = imin; i < imax; i++) {
      dx = ax - i*tg.gx;
      for (j = jmin; j < jmax; j++) {
	dy = ay - j*tg.gy;
	for (k = kmin; k < kmax; k++) {
	  if (ReadBoolGrid(tg, i, j, k) == 1) {
	    continue;
	  }
	  dz = az - k*tg.gz;
	  px = invU.data[0]*dx + invU.data[1]*dy + invU.data[2]*dz;
	  py = invU.data[3]*dx + invU.data[4]*dy + invU.data[5]*dz;
	  pz = invU.data[6]*dx + invU.data[7]*dy + invU.data[8]*dz;
	  if (px*px + py*py + pz*pz < R2) {

	    /*** Color this grid point ***/
	    WriteBoolGrid(tg, i, j, k, 1);
	    ncol++;
	  }
	}
      }
    }
  }

  return ncol;
}

/***=======================================================================***/
/*** CreateIgrid: create an integer grid.                                  ***/
/***=======================================================================***/
igrid CreateIgrid(long long int lx, long long int ly, long long int lz,
		  double gx, double gy, double gz, double ox, double oy,
		  double oz)
{
  int i, j;
  igrid tg;

  tg.lx = lx;
  tg.ly = ly;
  tg.lz = lz;
  tg.total_data = lx*ly*lz;
  tg.gx = gx;
  tg.gy = gy;
  tg.gz = gz;
  tg.invgx = 1.0/gx;
  tg.invgy = 1.0/gy;
  tg.invgz = 1.0/gz;
  tg.ox = ox;
  tg.oy = oy;
  tg.oz = oz;
  tg.data = (int*)calloc(lx*ly*lz, sizeof(int));
  tg.map = (int***)malloc(lx*sizeof(int**));
  for (i = 0; i < lx; i++) {
    tg.map[i] = (int**)malloc(lx*sizeof(int*));
    for (j = 0; j < ly; j++) {
      tg.map[i][j] = &tg.data[(i*ly + j)*lz];
    }
  }

  return tg;
}
