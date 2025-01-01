#ifndef ChargeFitFunctions
#define ChargeFitFunctions

#include "ChargeFitDS.h"
#include "pmeDirectDS.h"
#include "pmeRecipDS.h"
#include "TrajectoryDS.h"
#include "TopologyDS.h"
#include "TimingsDS.h"

fbook ReadEPotGrid(char* fname, prmtop *tp, coord *crd);

void FitCharges(fset *myfit, trajcon *tj, execon *etimers);

#endif
