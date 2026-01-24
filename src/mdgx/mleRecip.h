#ifndef mleRecipHeadings
#define mleRecipHeadings

#include "mleRecipDS.h"
#include "GridDS.h"
#include "pmeRecipDS.h"
#include "CrdManipDS.h"
#include "TrajectoryDS.h"
#include "TopologyDS.h"
#include "TimingsDS.h"

void SpreadBook(dbook Qc, dbook Q, g2gmap SPr, g2gmap SPc);

void IntrpBook(dbook Q, dbook Qc, g2gmap SPr, g2gmap SPc);

void MeshLevelLimits(reccon *rcinp, int lcon, int *llim1, int *hlim1,
                     int *llim2, int *hlim2);

dbook CutBCMeshX(reccon *rcinp, dbook *Urec, int lcon, coord *crd,
		 prmtop *tp);

void PrepMLE(reccon *rcinp, coord *crd, prmtop *tp);

void mleConvQBC(reccon *rcinp, Energy *sysUV, bckit *PPk, execon *etimers);

#endif
