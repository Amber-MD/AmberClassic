#!/bin/bash

CPP=${1}

echo "#ifndef MDGX_H" > mdgx.h
echo "#define MDGX_H" >> mdgx.h

cat >> mdgx.h << EOF

#include <sys/time.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "fftw3.h"
EOF

cat Constants.h Macros.h >> mdgx.h

for DSFI in BSpline Matrix Grid MPIMap CellManip Topology ChargeFit \
            CompFrc CrdManip mleRecip ParamFit pmeDirect pmeRecip \
            Restraints Timings Trajectory ; do
  ${CPP} -C -DPREP_API ${DSFI}DS.h -o tmp
  cat tmp >> mdgx.h
  rm tmp
done

cat >> mdgx.h << EOG
/***=======================================================================***/
/*** uform: all the potential function structs wrapped up in one type.     ***/
/***=======================================================================***/
struct PotentialFunction {
  prmtop tp;        // Topology
  dircon dcinp;     // Direct-space controls
  FrcTab Etab;      // Direct-space standard (coarse) lookup table
  FrcTab EHtab;     // Direct-space high-resolution lookkup table
  reccon rcinp;     // Reciprocal space controls
  bckit PPk;        // Convolution support data
};
typedef struct PotentialFunction uform;

/***=======================================================================***/
/*** mdsys: all structs for a molecular dynamics trajectory in one type.   ***/
/***=======================================================================***/
struct MolecularDynamicsSystem {
  coord crd;        // Coordinates (phone book organization of the system)
  cellgrid CG;      // Cell grid (spatial reorganization of the system)
  Energy sysUV;     // System energy decomposition and virial
  execon etimers;   // Timing data (for this system only)
};
typedef struct MolecularDynamicsSystem mdsys;

void GetPrmTop(prmtop *tp, trajcon *tj, int adjbnd);

void FreeTopology(prmtop *tp);

coord ReadRst(prmtop *tp, char* source);

void DestroyCoord(coord *crd);

void LongRangeVDW(prmtop *tp, dircon *dcinp);

FrcTab DirectSpaceR2(double range, double spc, double ewcoeff, int contderiv);

void FreeFrcTab(FrcTab *F);

cellgrid CreateCellGrid(coord *crd, dircon *dcinp, reccon *rcinp, prmtop *tp,
                        trajcon *tj, int sysnum);

#ifdef MPI
void LinkCellGrid(cellgrid *CG, coord *crd, reccon *rcinp);
#else
void LinkCellGrid(cellgrid *CG, reccon *rcinp);
#endif

void DestroyCellGrid(cellgrid *CG);

void DestroyEnergyTracker(Energy *sysUV);

void PrepPME(cellgrid *CG, reccon *rcinp, coord *crd);

void DestroyRecCon(reccon *rcinp, cellgrid *CG);

bckit CreateBCKit(reccon *rcinp, dbook *Q, coord *crd, prmtop *tp, int DoPlan);

void DestroyBCKit(bckit *PPk);

void AtomForces(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
                FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
                Energy *sysUV, execon *etimers, trajcon *tj);

double KineticEnergy(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj);

double SystemTemperature(cellgrid *CG, coord *crd, prmtop *tp, Energy *sysUV,
                         trajcon *tj, int updateKE);

#ifdef MPI
void SumTotalEnergy(cellgrid *CG, Energy *sysUV);
#else
void SumTotalEnergy(Energy *sysUV);
#endif

void InitHistory(coord *crd);

void InitBasicDircon(dircon *dcinp, double NBcut);

void InitBasicReccon(reccon *rcinp, dircon *dcinp);

void PrimeBasicTopology(prmtop *tp);

void InitBasicTrajcon(trajcon *tj);

void InitExecon(execon *tm);

void InitializeEnergy(Energy *sysUV, trajcon *tj, prmtop *tp,
                      int allocBondUdc);

uform InitPotential(char* topsrc, double NBcut, trajcon *tj);

mdsys LoadCoordToGrid(char* crdname, uform *U, trajcon *tj);

void MMForceEnergy(uform *U, mdsys *MD, trajcon *tj);

void AtomsToCells(coord *crd, cellgrid *CG, prmtop *tp);

void ImageBondedGroups(coord *crd, prmtop *tp);

void DefineMPITypes(trajcon *tj);

void FreeMPITypes(trajcon *tj);

void TransCrd(double* crds, int natom, double* tvec, double step);

void RotateCrd(double* crds, int natom, dmat U);

void FindCoordCenter(double* C, double* mass, int usem, int n, double* cofm);

void QuatAlign(double* frameI, double* frameII, int natom, double* mass,
               int m, dmat *U);

imat CreateImat(int N, int M);

void DestroyImat(imat *A);

imat ReallocImat(imat *A, int M, int N);

dmat CreateDmat(int N, int M, int prepFFT);

void DestroyDmat(dmat *A);

dmat ReallocDmat(dmat *A, int M, int N);

void CopyDmat(dmat *Ac, dmat *A, int Acex);

void DestroyTrajCon(trajcon *tj);

void DestroyCmat(cmat *A);

#endif 
EOG

