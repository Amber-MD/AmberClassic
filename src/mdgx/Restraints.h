#ifndef RestraintFunctions
#define RestraintFunctions

#include "GridDS.h"
#include "RestraintsDS.h"
#include "TrajectoryDS.h"
#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "CellManipDS.h"

fbook ReadRestraintGrid(char* source, coord *crd);

void WriteRestraintGrid(char* fname, fbook *L, trajcon *tj);

void WriteRestraintGridDB(char* fname, dbook *L, coord *crd, trajcon *tj);

void CellPrivateGrids(fbook* LV, int nLV, cellgrid *CG);

fbook ScoreOnGrid(fbook *Lgrd, cling *cr, trajcon *tj, coord *crd);

void CreateBellyMask(trajcon *tj, prmtop *tp);

#endif
