#ifndef CommandHeadings
#define CommandHeadings

#include "pmeDirectDS.h"
#include "pmeRecipDS.h"
#include "TopologyDS.h"
#include "TrajectoryDS.h"
#include "CrdManipDS.h"
#include "ChargeFitDS.h"
#include "ParamFitDS.h"
#include "IPolQDS.h"
#include "ConfigSampDS.h"
#include "SinglePointEvalDS.h"
#include "LoopBuilderDS.h"
#include "PeptideDS.h"
#include "ParmEditDS.h"

void ReadCommFile(dircon *dcinp, reccon *rcinp, prmtop* tp, trajcon *tj,
                  fset *myfit, prmset *myparms, ipqcon *ipqinp,
		  configs *cfsinp, spdata *spelist, loopcon *lpbinp,
		  pepcon *ppctrl, parmed *pmedit, char* source);

void SetGridDims(reccon *rcinp, coord *crd);

#endif
