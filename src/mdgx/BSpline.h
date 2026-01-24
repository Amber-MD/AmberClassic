#ifndef BSplineHeadings
#define BSplineHeadings

#include "pmeRecipDS.h"
#include "BSplineDS.h"
#include "CrdManipDS.h"
#include "CellManipDS.h"

double BSpln(double x, int ordr);

double dBSpln(double x, int ordr);

bmap CreateBmap(reccon *rcinp, int natom);

void DestroyBmap(bmap *A);

void SplCoeff(coord *crd, bmap *pmmap, reccon *rcinp);

void CellSplCoeff(cell *C, coord *crd, reccon *rcinp);

#endif
