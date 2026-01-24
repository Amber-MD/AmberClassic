#ifndef IPOLQ_FUNCS
#define IPOLQ_FUNCS

#include "CrdManipDS.h"
#include "CellManipDS.h"
#include "TopologyDS.h"
#include "pmeDirectDS.h"
#include "pmeRecipDS.h"
#include "TrajectoryDS.h"
#include "TimingsDS.h"
#include "AmberNetcdf.h"
#include "IPolQDS.h"

double* EvenSphere(int n, int verbose);

double AtomCode(double mass, int znum, double charge, char* tlet);

void PerformIPolQ(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
		  FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
		  trajcon *tj, Energy *sysUV, execon *etimers, cdftrj* Acdf,
		  ipqcon *ipqinp);

#endif
