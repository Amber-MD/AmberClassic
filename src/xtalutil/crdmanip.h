#include "matrix.h"

#ifndef CRDMANIP_FUNCS
#define CRDMANIP_FUNCS

void TransCrd(double* crds, int n_atoms, double* tvec, double step);

void RotateCrd(double* crds, int n_atoms, dmat U);

void BeardRotMat(double om_x, double om_y, double om_z, dmat r);

void ExtCrds(double* crds, int num_atoms, double* ext);

void FindCrdCenter(double* crds, double* mass, int m, int n_atoms, double* cm);

void CmpXfrm(double* cd, dmat U, dmat invU);

void ReImage(double* crds, int n);

dmat UndoOct();

double SimpleReIm(double dx);

int RenderFrame(double* crd, double* boxd, int natm, FILE *inp);

#endif
