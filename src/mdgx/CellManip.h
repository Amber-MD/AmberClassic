#ifndef CellManipHeadings
#define CellManipHeadings

#include "CrdManipDS.h"
#include "pmeDirectDS.h"
#include "pmeRecipDS.h"
#include "TopologyDS.h"
#include "CellManipDS.h"
#include "CompFrcDS.h"
#include "TrajectoryDS.h"

cell CreateCell(int maxatom, int* ordr);

void HessianNorms(dmat *invU, double* cdepth);

void TakeCellGridDims(int* cdim, double* cdepth, coord *crd, dircon *dcinp);

int IsCentralAtom(double* corig, double* ccen, double* atmloc, double* Umat,
		  double* ng, int isortho);

cellgrid CreateCellGrid(coord *crd, dircon *dcinp, reccon *rcinp, prmtop *tp,
			trajcon *tj, int sysnum);

void DestroyCell(cell *C);

void DestroyCellGrid(cellgrid *CG);

void ComputeCellOrigins(cellgrid *CG, coord *crd);

void ComputeCellNormals(coord *crd);

void UpdateCellGPS(cellgrid *CG);

void UploadCellPosForce(cell *C, coord *crd);

void DownloadCellForces(cell *C, coord *crd);

void AtomsToCells(coord *crd, cellgrid *CG, prmtop *tp);

int SortAtomID(const void *atmA, const void *atmB);

void CopyCellContents(cellgrid *CG, int Oi, int Oj, int Ok, int Ni, int Nj,
		      int Nk, int mx, int my, int mz, int dreg, double Mcut,
		      coord *crd, prmtop *tp);

void ComputeDirectAll(cellgrid *CG, FrcTab *EFrc, dircon *dcinp, prmtop *tp);

void MixCellGrids(cellgrid *CGA, cellgrid *CGB, trajcon *tj);

void ZeroCellForces(cellgrid *CG);

#ifdef MPI
void ShareCoordinates(cellgrid *CG, trajcon *tj, double Mcut, coord *crd,
                      prmtop *tp, int ct);
#else
void ShareCoordinates(cellgrid *CG, double Mcut, coord *crd, prmtop *tp,
		      int ct);
#endif

void ShareTriple2(cellgrid *CG, double Mcut, coord *crd, prmtop *tp);

void CountMassive(cellgrid *CG, prmtop *tp);

void DirectTriple2(cell *C, coord *crd, FrcTab *Etab, dircon *dcinp,
		   prmtop *tp, Energy *sysUV);

void CellLRvdw(cell *C, coord *crd, prmtop *tp, Energy *sysUV);

void MapCellForcesToAtoms(cellgrid *CG, coord *crd);

void MapListForcesToCells(cellgrid *CG, coord *crd);

#ifdef MPI
void LinkCellGrid(cellgrid *CG, coord *crd, reccon *rcinp);

void MergeCellForces(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj);

void UpdateCells(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj);
#else
void LinkCellGrid(cellgrid *CG, reccon *rcinp);

void MergeCellForces(cellgrid *CG, coord *crd, prmtop *tp);

void UpdateCells(cellgrid *CG, coord *crd, prmtop *tp);
#endif

double CalcCellR2Cost(cell *C, cellgrid *CG);

#endif
