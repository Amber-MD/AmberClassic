#ifndef BarostatFunctions
#define BarostatFunctions

#include "TrajectoryDS.h"
#include "CompFrcDS.h"
#include "CrdManipDS.h"
#include "pmeDirectDS.h"
#include "pmeRecipDS.h"
#include "BSplineDS.h"
#include "CellManipDS.h"
#include "TimingsDS.h"

void PrepMCBarostat(trajcon *tj, coord *crd);

double CurrentSystemPressure(Energy *sysUV, coord *crd);

void GridScaleCoord(coord *crd, prmtop *tp, double* chi, dircon *dcinp,
		    reccon *rcinp, bckit *PPk, cellgrid *CG, trajcon *tj);

void BerendsenBarostat(coord *crd, prmtop *tp, trajcon *tj, Energy *sysUV,
                       dircon *dcinp, reccon *rcinp, bckit *PPk, cellgrid *CG);

void CheckPrimarySectors(cellgrid*CG);

void MonteCarloBarostat(coord *crd, prmtop *tp, trajcon *tj, Energy *sysUV,
			dircon *dcinp, FrcTab *Etab, FrcTab *EHtab,
			reccon *rcinp, bckit *PPk, cellgrid *CG,
			execon *etimers);

#endif
