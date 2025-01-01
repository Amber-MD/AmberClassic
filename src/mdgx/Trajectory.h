#ifndef TrajectoryHeadings
#define TrajectoryHeadings

#include "AmberNetcdfDS.h"

#include "CrdManipDS.h"
#include "CellManipDS.h"
#include "TopologyDS.h"
#include "TrajectoryDS.h"
#include "pmeDirectDS.h"
#include "pmeRecipDS.h"
#include "CompFrcDS.h"
#include "TimingsDS.h"

void UpdateStepNumber(trajcon *tj, int newstep);

coord ReadRst(prmtop *tp, char* source, double *starttime);

dmat ReadCrdTraj(prmtop *tp, char* source, int readbox);

dmat ReadCDFTraj(prmtop *tp, char* source, int readbox);

void SpliceFileName(trajcon *tj, char* base, char* suff, char* fname,
                    int dprc);

void InitializeEnergy(Energy *sysUV, trajcon *tj, prmtop *tp,
		      int allocBondUdc);

void DestroyEnergyTracker(Energy *sysUV);

void ExtendCoordinates(coord *tc, prmtop *tp);

int ValidCoord(char* fname);

coord InitCoords(prmtop *tp, trajcon *tj, int n);

void OpenDiagnosticsFile(trajcon *tj, prmtop *tp, dircon *dcinp, reccon *rcinp,
			 FrcTab *Etab, coord *crd, cellgrid *CG, int n);

void CloseDiagnosticsFile(trajcon *tj, prmtop *tp, reccon *rcinp,
			  Energy *sysUV, execon *etimers, cellgrid* CG,
			  int format);

void InitAverageEnergy(Energy *sysUV);

void AccAverageEnergy(Energy *sysUV);

void PrintStateInfo(Energy *sysUV, trajcon *tj, prmtop *tp, FILE *outp,
		    int format);

void WriteDiagnostics(trajcon *tj, prmtop *tp, dircon *dcinp, reccon *rcinp,
                      FrcTab *Etab, cellgrid *CG, coord *crd, Energy *sysUV,
		      execon *etimers, int n);

void WriteRst(cellgrid *CG, coord *tc, prmtop *tp, trajcon *tj, int n,
              int prbox);

void WriteCrd(cellgrid *CG, coord *crd, int cvf, trajcon *tj, prmtop *tp,
	      int n);

void WriteCDF(cellgrid *CG, coord *crd, int cvf, trajcon *tj, cdftrj *Acdf,
	      prmtop *tp, int n);

void ReprintInputFile(trajcon *tj, char* key1, char* key2, int printkey,
                      FILE *outp);

void DestroyTrajCon(trajcon *tj);

void SynchronizeCoordinates(cellgrid* CG, trajcon *tj);

dmat ReadListOfDoubles(char* fname);

#endif
