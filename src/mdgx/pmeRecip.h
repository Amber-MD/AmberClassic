#ifndef pmeRecipHeadings
#define pmeRecipHeadings

#include "GridDS.h"
#include "pmeRecipDS.h"
#include "CrdManipDS.h"
#include "TrajectoryDS.h"
#include "TopologyDS.h"
#include "TimingsDS.h"
#include "CellManipDS.h"

int Factor2357(double x, int min2);

void SetMeshDims(int* ng, double* gdim);

double* LoadPrefac(int ordr, int ng);

double* LoadMVec(int ng);

int TestUnitCellOrtho(double* gdim);

double SelfEnergyCorrection(double ewcoeff, prmtop *tp);

bckit CreateBCKit(reccon *rcinp, dbook *Q, coord *crd, prmtop *tp,
		  unsigned int plan);

void UpdateBCKit(reccon *rcinp, coord *crd, prmtop *tp, bckit *PPk);

void DestroyBCKit(bckit *PPk);

void ConvQBC(reccon *rcinp, coord *crd, dbook *Q, bckit *PPk, execon *etimers);

void ConvQBCnrg(reccon *rcinp, coord *crd, dbook *Q, bckit *PPk, Energy *sysUV,
		execon *etimers);

void ConvQBCnrgvir(reccon *rcinp, coord *crd, dbook *Q, bckit *PPk,
		   Energy *sysUV, execon *etimers);

void PrepPME(cellgrid *CG, reccon *rcinp, coord *crd);

void DestroyRecCon(reccon *rcinp, cellgrid *CG);
 
void DestroyAdvancedRecCon(reccon *rcinp, cellgrid *CG);


#endif
