#ifndef ChargeMapHeadings
#define ChargeMapHeadings

#include "GridDS.h"
#include "pmeRecipDS.h"
#include "BSplineDS.h"
#include "TopologyDS.h"
#include "CellManipDS.h"

void CellQ2Book(cell *C, reccon *rcinp, dbook *Q);

void CellIntrpFrc(cell *C, coord *crd, reccon *rcinp, dbook *U);

#endif
