#ifndef mleOptHeadings
#define mleOptHeadings

#include "pmeRecipDS.h"
#include "MatrixDS.h"
#include "GridDS.h"
#include "mleOptDS.h"
#include "TopologyDS.h"

void PrepareLayerTargets(reccon *rcinp, dbook *Qact, dbook *Qcrs, int nlyr,
                         dmat* Utrg, dmat* Qtrg);

bssf LayerBasisFunctions(dmat *pLrec, int nlyr, int maxlyr, coord *crd);

void OptimizeLayer(reccon *rcinp, int lcon, dbook *Lrec, dbook *Qact,
                   dbook *Qcrs, int nlyr, coord *crd);

void OptimizeMesh(dbook *Uact, dbook *Lrec, reccon *rcinp, coord *crd,
                  prmtop *tp, int lcon);

#endif
