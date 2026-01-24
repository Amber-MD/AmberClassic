#ifndef ParamReadFuncs
#define ParamReadFuncs

#include "ParamFitDS.h"
#include "TrajectoryDS.h"

int CrossRefAtomType(prmset *mp, char* aname);

void ReadParmFile(prmset *mp, trajcon *tj);

void ReadFrcmodFile(prmset *mp, trajcon *tj, int nfmod);

void BubbleChar4(char* a, char* b);

void RecastAtomTypes(prmset *mp);

void CleaveAtomTypes(prmset *mp);

void CreateNMROperation(nmroper *thisop);

void DestroyNMROperation(nmroper *thisop);

int CompareOperations(nmroper *opA, nmroper *opB, prmset *mp, int numerics);

void ReadNMROperationsFile(prmset *mp);

#endif
