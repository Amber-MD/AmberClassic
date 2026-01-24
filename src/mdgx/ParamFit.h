#ifndef ParamFitHeadings
#define ParamFitHeadings

#include "ParamFitDS.h"
#include "TopologyDS.h"
#include "TrajectoryDS.h"

int str4cmp(char* A, char* B);

int strendswith(char* A, char* B);

int TypeCompare(char* T1a, char* T1b, char* T1c, char* T1d, char* T2a,
                char* T2b, char* T2c, char* T2d, int order, int strict);

double EvaluateCMAPEnergy(double val1, double val2, dmat *cmsurf,
                          dmat *cmdphi, dmat *cmdpsi, dmat *cmd2pp);

dmat CompExclMatrix(prmtop *tp, imat *partnb);

imat BuildMatchMatrix(mmsys *myconf, nmroper *myop);

void FitParams(prmset *mp, trajcon *tj);

#endif
