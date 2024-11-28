#include "matrix.h"

#ifndef QUATALIGN_FUNCS
#define QUATALIGN_FUNCS

void TRED2(double** A, int n, double* d, double* e);

void TQLI(double* d, double* e, int n, double** z);

void QuatAlign(double* frameI, double* frameII, int num_atoms, double* mass,
	       int m, dmat *U);

double WeightRMSD(double* mola, double* molb, double* mass, int m, int natom);

#endif
