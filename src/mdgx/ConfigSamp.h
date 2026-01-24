#ifndef ConfigFunctions
#define ConfigFunctions

#include "ConfigSampDS.h"
#include "TrajectoryDS.h"
#include "CrdManipDS.h"

void InitConfigs(configs *cfsinp);

void VerticalToCoord(dmat *crdpop, int ncol, coord *tc);

void Config2Pdb(configs *cfsinp, coord *tc, trajcon *tj, prmtop *tp,
		Energy* sysUV, int cidx, int stlidx);

void SampleConfigurations(configs *cfsinp, trajcon *tj, prmtop *tp);

#endif
