#ifndef CrdManipHeadings
#define CrdManipHeadings

#include "CrdManipDS.h"
#include "MatrixDS.h"
#include "TopologyDS.h"

coord CreateCoord(int natm);

coord CopyCoord(coord *crd);

void DestroyCoord(coord *crd);

void TransCrd(double* crds, int natom, double* tvec, double step);

void RotateCrd(double* crds, int natom, dmat U);

void FindCoordCenter(double* C, double* mass, int usem, int n, double* cofm);

void BeardRotMat(double om_x, double om_y, double om_z, dmat *r);

void CompXfrm(double* cd, dmat U, dmat invU);

void OrthoReim(double *dx, double *dy, double *dz, dmat *U, dmat *invU);

void NonOrthoReim(double *dx, double *dy, double *dz, dmat *U, dmat *invU);

void ImageBondedGroups(coord *crd, prmtop *tp);

void QuatAlign(double* frameI, double* frameII, int natom, double* mass,
               int m, dmat *U);

#endif
