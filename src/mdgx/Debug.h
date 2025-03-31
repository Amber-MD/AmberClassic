#ifndef DebugFunctions
#define DebugFunctions

#include "CrdManipDS.h"
#include "CellManipDS.h"

int FindAtomInCell(cell *C, int atmid, int ireq);

void PrintCellContents(cellgrid *CG, char* outname, char* varname);

void FindAllInstances(cellgrid *CG, coord *crd, int aid);

void CheckCellContents(cellgrid *CG, coord *crd, prmtop *tp, int seekEP,
		       int chkbounds, int chkforces, int StopOnError,
		       char* announce);

#ifdef MPI
void PrintRecvInfo(int maxsize, char* tname, int recver, int sender, int tag,
		   MPI_Comm comm);

void PrintSendInfo(int maxsize, char* tname, int recver, int sender, int tag,
		   MPI_Comm comm, void* data, int verbose);

void DisplaySend(ashr *cshr, int isend, int offset, int msgsize, cellgrid *CG,
		 trajcon *tj, MPI_Datatype msgtype);
#endif

void Torque(prmtop *tp, coord *crd, int resid);

void PrintResidueExclusions(prmtop *tp, int rid);

void PrintResidueForces(cellgrid *CG, coord *crd, prmtop *tp, int rid);

void CellChecksum(cellgrid *CG, coord *crd, prmtop *tp, int nsctr, int doloc,
		  int dovel, int dofrc, int doploc, int dopvel, int dopfrc,
		  char* tagmsg, char* finmsg, double sleeptm);

cellgrid PrepCellGridCopy(cellgrid *CG);

void CopyCellGrid(cellgrid *CY, cellgrid *CG);

void Sector2Sector(cell *C, int sA, int sB, dircon *dcinp, prmtop *tp,
                   FrcTab *EFrc, Energy* sysUV, int dofrc, int donrg,
                   int dovir);

void ExhaustInteractions(coord *crd, prmtop *tp, FrcTab *EFrc, dircon *dcinp,
                         Energy *sysUV, int dofrc, int donrg, int dovir,
			 double neartol);

void ListBondedIntr(coord *crd, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
                    Energy *sysUV, int qform);

#endif
