#ifndef MPIMapFunctions
#define MPIMapFunctions

#include "TrajectoryDS.h"
#include "CrdManipDS.h"
#include "pmeRecipDS.h"
#include "CellManipDS.h"
#include "MPIMapDS.h"

void SelectSystemsToTend(trajcon *tj);

void MeshSumSPME(reccon *rcinp, cellgrid *CG);

void MapProcessors(trajcon *tj);

void InitLoadBalance(cellgrid *CG, reccon *rcinp);

void MapProcessMeshFootprint(cellgrid *CG, reccon *rcinp, coord *crd);

#ifdef MPI
int InitMeshGatherSPME(reccon *rcinp, cellgrid *CG, MPI_Request* req,
		       int scatter);

void SendMeshPart(reccon *rcinp, cellgrid *CG, MPI_Request *req, int ncpu);

void RecvMeshPart(reccon *rcinp, cellgrid *CG, int ncpu);
#endif

void CellGridEquipComm(cellgrid *CG, trajcon *tj);

void MapProcessCellSharing(cellgrid *CG);

void DestroyAshr(ashr *A);

void AllocatePooledBuffers(cellgrid *CG);

void PlanCoordRedux(cellgrid *CG);

#ifdef MPI
void CoordinateReduction(cellgrid *CG, coord *crd, prmtop *tp);
#endif

#endif
