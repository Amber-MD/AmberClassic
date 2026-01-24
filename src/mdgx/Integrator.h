#ifndef IntegratorHeadings
#define IntegratorHeadings

#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "CompFrcDS.h"
#include "CellManipDS.h"
#include "pmeRecipDS.h"
#include "pmeDirectDS.h"
#include "BSplineDS.h"
#include "TrajectoryDS.h"
#include "TimingsDS.h"
#include "AmberNetcdfDS.h"

void InitHistory(coord *crd);

void AtomForces(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
                FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
		Energy *sysUV, execon *etimers, trajcon *tj);

#ifdef MPI
void SumTotalEnergy(cellgrid *CG, Energy *sysUV);
#else
void SumTotalEnergy(Energy *sysUV);
#endif

void InitVelocities(coord* crd, cellgrid* CG, prmtop* tp, dircon *dcinp,
                    FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit* PPk,
		    trajcon *tj, Energy* sysUV, execon *etimers, cdftrj* Acdf,
		    int n);

void InitVelocitiesTI(coord* crd, cellgrid* CG, prmtop* tp, dircon *dcinp,
		      FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit* PPk,
		      trajcon *tj, Energy* sysUV, execon *etimers,
		      cdftrj* Acdf, int n);

void Dynamics(coord* crd, cellgrid* CG, prmtop* tp, dircon *dcinp,
	      FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit* PPk,
	      trajcon *tj, Energy* sysUV, execon *etimers, cdftrj* Acdf,
	      int n);

void DynamicsTI(coord* crd, cellgrid* CG, prmtop* tp, dircon *dcinp,
	        FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit* PPk,
	        trajcon *tj, Energy* sysUV, execon *etimers, cdftrj* Acdf,
		int n);

void Minimization(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
		  FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
		  trajcon *tj, Energy *sysUV, execon *etimers, cdftrj* Acdf);

void PrepForceReport(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
                     FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
		     trajcon *tj, execon *etimers);

#endif
