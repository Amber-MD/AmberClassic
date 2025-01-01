#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mdgxVector.h"
#include "Topology.h"
#include "CellManip.h"
#include "CrdManip.h"
#include "Nonbonded.h"
#include "Bonded.h"
#include "Matrix.h"

//-----------------------------------------------------------------------------
// FindAtomInSector: find an atom within a cell sector based on its ID number
//                   in the master topology.                     
//
// Arguments:
//   C:      the cell of interest                                       
//   atmid:  the ID number of the atom of interest                      
//   ireq:   flag to indicate that the atom must be found               
//   sctr:   the sector of the cell to search                           
//-----------------------------------------------------------------------------
static int FindAtomInSector(cell *C, int atmid, int ireq, int sctr)
{
  int i, imin, imax, ihalf;
  atomc *catm;

  imin = 0;
  imax = C->nr[sctr]-1;
  ihalf = imin + (imax - imin)/2;
  catm = C->map[sctr];
  while (imin <= imax ) {

    // First, check to see if the atom
    // is just not in this sector     
    if (catm[imin].id > atmid || catm[imax].id < atmid) {
      break;
    }
    if (catm[ihalf].id > atmid) {
      imax = ihalf-1;
      ihalf = imin + (imax - imin)/2;
    }
    else if (catm[ihalf].id < atmid) {
      imin = ihalf+1;
      ihalf = imin + (imax - imin)/2;
    }
    if (catm[ihalf].id == atmid) {
      return sctr*C->maxatom + ihalf;
    }
    if (imax - imin <= 2) {
      for (i = imin; i <= imax; i++) {
        if (catm[i].id == atmid) {
          return sctr*C->maxatom + i;
	}
      }

      // If we're still here, the atom
      // was not found in this sector 
      break;
    }
  }

  // If we're still here the atom was not found at all.
  if (ireq == 1) {

    // Exit if the atom was needed
    printf("FindAtomInSector >> Error.  Atom %d could not be located.\n",
           atmid);
    printf("FindAtomInSector >> Sector contents:\n");
    printf("Sector %d:\n", sctr);
    for (i = 0; i < C->nr[sctr]; i++) {
      printf("%5d ", C->map[sctr][i].id);
    }
    printf("\n");
    exit(1);
  }

  return -1;
}

//-----------------------------------------------------------------------------
// FindAtomInCell: find an atom within a cell based on its ID number in the
//                 master topology.  All eight sectors, including the home cell
//                 and all of its imported regions, will be searched.  This
//                 routine returns an integer corresponding to the position of
//                 the atom in the data field of the cell struct.  If the
//                 requested atom cannot be found, this routine will cause the
//                 program to exit in error if the flag ireq is set to 1.
//
// Arguments:
//   C:      the cell of interest                                       
//   atmid:  the ID number of the atom of interest                      
//   ireq:   flag to indicate that the atom must be found               
//-----------------------------------------------------------------------------
int FindAtomInCell(cell *C, int atmid, int ireq)
{
  int i, j, jmin, jhalf, jmax;
  atomc *catm;

  // Search all eight sectors in this cell
  for (i = 0; i < 8; i++) {
    catm = C->map[i];
    jmin = 0;
    jmax = C->nr[i]-1;
    jhalf = jmin + (jmax - jmin)/2;
    while (jmin <= jmax ) {

      // First, check to see if the atom
      // is just not in this sector     
      if (catm[jmin].id > atmid || catm[jmax].id < atmid) {
        break;
      }
      if (catm[jhalf].id > atmid) {
        jmax = jhalf-1;
        jhalf = jmin + (jmax - jmin)/2;
      }
      else if (catm[jhalf].id < atmid) {
        jmin = jhalf+1;
        jhalf = jmin + (jmax - jmin)/2;
      }
      if (catm[jhalf].id == atmid) {
        return i*C->maxatom + jhalf;
      }
      if (jmax - jmin <= 2) {
        for (j = jmin; j <= jmax; j++) {
          if (catm[j].id == atmid) {
            return i*C->maxatom + j;
          }
        }

        // If we're still here, the atom
        // was not found in this sector 
        break;
      }
    }
  }

  // If we're still here the atom was not found at all.
  if (ireq == 1) {

    // Exit if the atom was needed
    printf("FindAtomInCell >> Error.  Atom %d could not be located.\n", atmid);
    printf("FindAtomInCell >> Cell contents:\n");
    for (i = 0; i < 8; i++) {
      printf("Sector %d:\n", i);
      jmin = 0;
      for (j = 0; j < C->nr[i]; j++) {
        printf("%5d ", C->map[i][j].id);
        jmin++;
        if (jmin == 12) {
          jmin = 0;
          printf("\n");
        }
      }
      printf("\n");
    }
    exit(1);
  }

  return -1;
}

//-----------------------------------------------------------------------------
// PrintCellContents: loop over all cells and print the contents of the primary
//                    sectors to a file.
//
// Arguments:                                                           
//   CG:       the cell grid                                            
//   outname:  the name of the output file (MatLab format)              
//   varname:  the name of the variable for accessing data in MatLab    
//-----------------------------------------------------------------------------
void PrintCellContents(cellgrid *CG, char* outname, char* varname)
{
  int i, j;
  cell *C;
  FILE *outp;

  outp = fopen(outname, "w");
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    fprintf(outp, "%% Cell %3d: %4d atoms\n", i, C->nr[0]);
    fprintf(outp, "%s(1:%d,:,%d) = [\n", varname, C->nr[0], i+1);
    for (j = 0; j < C->nr[0]; j++) {
      fprintf(outp, "%4d %4d %16.10lf %16.10lf %16.10lf %16.10lf\n",
	      C->data[j].id, C->data[j].lj, C->data[j].q, C->data[j].loc[0],
	      C->data[j].loc[1], C->data[j].loc[2]);
    }
    fprintf(outp, "];\n");
  }
  fclose(outp);
}

//-----------------------------------------------------------------------------
// FindAllInstances: find all instances of an atom in all cells.        
//                                                                      
// Arguments:                                                           
//   CG:      the cell grid                                             
//   crd:     the coordinates                                           
//   aid:     the atom number to search for                             
//-----------------------------------------------------------------------------
void FindAllInstances(cellgrid *CG, coord *crd, int aid)
{
  int ci, cj, ck, j, k;
  cell *C;
  atomc *atmt;

  // Print some basic information about the atom
  printf("Finding all instances of atom %d:\n", aid);
  printf("Most recent global location: [ %16.10lf %16.10lf %16.10lf ]\n",
	 crd->loc[3*aid], crd->loc[3*aid+1], crd->loc[3*aid+2]);

  // Loop over all cells and find instances in any sector
  for (ci = 0; ci < CG->ng[0]; ci++) {
    for (cj = 0; cj < CG->ng[1]; cj++) {
      for (ck = 0; ck < CG->ng[2]; ck++) {
	C = &CG->map[ci][cj][ck];
	for (j = 0; j < 8; j++) {
	  for (k = 0; k < C->nr[j]; k++) {
	    atmt = &C->map[j][k];
	    if (atmt->id == aid) {
	      printf("Atom %6d found in cell [%2d][%2d][%2d]:  ", aid, ci, cj,
		     ck);
	      printf("Sector %d -> Element %4d of %4d (max %4d)\n", j, k,
		     C->nr[j], C->maxatom);
	      printf("Cell Loc: [ %16.10lf %16.10lf %16.10lf ]\n",
		     atmt->loc[0], atmt->loc[1], atmt->loc[2]);
	      printf("Cell Lim: [ %16.10lf %16.10lf %16.10lf ] ->\n"
		     "          [ %16.10lf %16.10lf %16.10lf ]\n",
		     C->orig[0], C->orig[1], C->orig[2],
		     C->orig[0] + CG->celldim[0] + CG->celldim[3],
		     C->orig[1] + CG->celldim[1] + CG->celldim[3],
		     C->orig[2] + CG->celldim[2] + CG->celldim[3]);
	    }
	  }
	}
      }
    }
  }
}

//-----------------------------------------------------------------------------
// PrimeCCCErrMsg: initialize the error message for CheckCellContents() below.
//                 Critical information supplied includes the process, system,
//                 and atom identifications.
//
// Arguments:                                                           
//   CG:           the cell grid                                        
//   C:            the cell within the cell grid which has the error    
//   sctr:         the sector of the cell wherein the problem is found  
//   cidx:         the index of the troublesome atom within the cell (the local
//                 variable tpidx gives the "absolute" index within the actual
//                 topology)
//   tp:           the topology                                         
//   errname:      the name of the error                                
//-----------------------------------------------------------------------------
static int PrimeCCCErrMsg(cellgrid *CG, cell *C, int sctr, int cidx,
			  prmtop *tp, char* errname, char* ErrMsg)
{
  int tpidx, residx;
  char *errtmp;

  tpidx = (sctr >= 0) ? C->map[sctr][cidx].id : cidx;
  residx = LocateResID(tp, tpidx, 0, tp->nres);
  sprintf(ErrMsg, ">> %s on process %4d of system %3d:\n", errname, CG->tid,
	  CG->sysID);
  errtmp = &ErrMsg[strlen(ErrMsg)];
  sprintf(errtmp, ">> Atom %6d (%.4s %4d %.4s)", tpidx,
	  &tp->ResNames[4*residx], residx, &tp->AtomNames[4*tpidx]);
  if (sctr >= 0) {
    errtmp = &ErrMsg[strlen(ErrMsg)];
    sprintf(errtmp, " in cell [ %3d %3d %3d ][ %d ] -> %4d:", C->gbin[0],
	    C->gbin[1], C->gbin[2], sctr, cidx);
  }
  errtmp = &ErrMsg[strlen(ErrMsg)];
  sprintf(errtmp, "\n");

  return strlen(ErrMsg);
}

//-----------------------------------------------------------------------------
// CheckCellContents: this debugging function performs a complex assert,
//                    verifying that the contents of each cell sector really
//                    belong there.  This function will work in parallel
//                    implementations, and should serve as a useful debugging
//                    tool.
//
// Arguments:
//   CG:           the cell grid                                        
//   crd:          the coordinates                                      
//   tp:           the topology                                         
//   seekEP:       flag to activate seeking of extra points in certain  
//                 tests                                                
//   chkbounds:    flag to call for bounds checking on all cells--set to
//                 zero to turn this off if constraints have just been  
//                 applied or Verlet integration has just occured and   
//                 the UpdateCells() function of the CellManip library  
//                 has not yet been called                              
//   chkforces:    flag to call for force checking, which will report   
//                 any unusually large forces                           
//   StopOnError:  flag to have the program abort if errors are detected
//   announce:     message to print, giving an indication of where the  
//                 function was called                                  
//-----------------------------------------------------------------------------
void CheckCellContents(cellgrid *CG, coord *crd, prmtop *tp, int seekEP,
		       int chkbounds, int chkforces, int StopOnError,
		       char* announce)
{
  int i, j, k, m, cp3, problem, halt, msgidx;
  int* cellplace;
  double fmag;
  double climit[3], corig[3], bcrd[3], rlim[3];
  char* ErrMsg;
  char *errtmp;
  atomc *catm;
  cell *C, *CNew;

  // Print the announcement message
  if (announce[0] != '\0') {
    if (CG->tid == 0) {
      printf("CheckCellContents: %s\n", announce);
    }
  }

  // Allocate memory for the error message,
  // for duplicate atoms thus far located, 
  // and other needs                       
  ErrMsg = (char*)malloc(1024*sizeof(char));

  // Check to see that no atoms are outside of their cell boundaries
  halt = 0;
  if (chkbounds > 0) {
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];

      // Shortcut for cell limits
      DMatVecMult(&crd->U, C->midp, climit);
      DMatVecMult(&crd->U, C->orig, corig);
      for (j = 0; j < 3; j++) {
	climit[j] += 2.0*CG->celldim[4+j];
      }
      for (j = 0; j < chkbounds; j++) {
	catm = C->map[j];
	for (k = 0; k < C->nr[j]; k++) {

	  // Determine if the atom is outside the cell boundary
	  DMatVecMult(&crd->U, catm[k].loc, bcrd);
	  problem = 0;
	  for (m = 0; m < 3; m++) {
	    if (bcrd[m] < corig[m] || bcrd[m] >= climit[m] ||
		bcrd[m] != bcrd[m]) {
	      problem = 1;
	    }
	  }
	  if (problem == 1) {
	    msgidx = PrimeCCCErrMsg(CG, C, j, k, tp, "Bounds alert", ErrMsg);
	    errtmp = &ErrMsg[msgidx];
	    sprintf(errtmp, "  Real location: %10.5lf %10.5lf %10.5lf\n"
		    "  Unit cell loc: %10.5lf %10.5lf %10.5lf\n",
		    catm[k].loc[0], catm[k].loc[1], catm[k].loc[2], bcrd[0],
		    bcrd[1], bcrd[2]);
	    errtmp = &ErrMsg[strlen(ErrMsg)];
	    for (m = 0; m < 3; m++) {
	      rlim[m] = C->orig[m] + (C->midp[m]-C->orig[m]) *
		(1.0+2.0*CG->celldim[4+m]);
	    }
	    if (crd->isortho == 1) {
	      sprintf(errtmp, "  Real limits:   %10.5lf %10.5lf %10.5lf x\n"
		      "                 %10.5lf %10.5lf %10.5lf\n",
		      C->orig[0], C->orig[1], C->orig[2], rlim[0], rlim[1],
		      rlim[2]);
	    }
	    else {
	      sprintf(errtmp, "  Unit cell lim: %10.5lf %10.5lf %10.5lf x\n"
		      "                 %10.5lf %10.5lf %10.5lf\n",
		      corig[0], corig[1], corig[2], climit[0], climit[1],
		      climit[2]);
	    }
	    printf("%s\n", ErrMsg);
	    halt = 1;
	  }
	}
      }
    }
  }

  // Check to see that there are no unusually large forces
  if (chkforces == 1) {
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      for (j = 0; j < C->nr[0]; j++) {
	fmag = C->data[j].frc[0]*C->data[j].frc[0];
	fmag += C->data[j].frc[1]*C->data[j].frc[1];
	fmag += C->data[j].frc[2]*C->data[j].frc[2];
	fmag = sqrt(fmag);
	if (fmag > 10000.0 || fmag != fmag) {
	  msgidx = PrimeCCCErrMsg(CG, C, 0, j, tp, "Large force", ErrMsg);
	  errtmp = &ErrMsg[msgidx];
	  sprintf(errtmp, "   Force = [ %16.4lf %16.4lf %16.4lf ]\n",
		  C->data[j].frc[0], C->data[j].frc[1], C->data[j].frc[2]);
	  printf("%s\n", ErrMsg);
	  halt = 1;
	}
      }
    }
  }

  // Check to see that all atoms are accounted for
  cellplace = (int*)malloc(crd->natom*sizeof(int));
  SetIVec(cellplace, crd->natom, -1);
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    catm = C->data;
    for (m = 0; m < C->nr[0]; m++) {

      // First check if this atom has already been located
      if (catm[m].id < 0 || catm[m].id >= crd->natom) {
	printf("Error: atom [ %2d %2d %2d ]->[ 0, %4d ] ID = %d!\n",
	       i, j, k, m, catm[m].id);
	exit(1);
      }
      cp3 = catm[m].id;
      if (cellplace[cp3] >= 0) {
	CNew = &CG->data[cellplace[cp3]];
        msgidx = PrimeCCCErrMsg(CG, C, 0, k, tp, "Multiple cell locations",
				ErrMsg);
        errtmp = &ErrMsg[msgidx];
	sprintf(errtmp, "   Also located in [ %3d %3d %3d ]!\n", CNew->gbin[0],
		CNew->gbin[1], CNew->gbin[2]);
	printf("%s\n", ErrMsg);
	halt = 1;
      }

      // Mark this atom as located
      cellplace[cp3] = C->gbin[3];
    }
  }

  // Check to see that all cells contain properly ordered lists
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < 8; j++) {
      catm = C->map[j];
      for (k = 1; k < C->nr[j]; k++) {
	if (catm[k].id < catm[k-1].id) {
	  PrimeCCCErrMsg(CG, C, j, k, tp, "Improper order", ErrMsg);
	  printf("%s", ErrMsg);
	  PrimeCCCErrMsg(CG, C, j, k-1, tp, "   Out of order relative to",
			 ErrMsg);
	  printf("%s\n", ErrMsg);
	  halt = 1;
	}
      }
    }
  }

  // Check to see that all atoms are accounted for
#ifdef MPI
  int* buffplace;
  buffplace = (int*)malloc(crd->natom*sizeof(int));
  MPI_Reduce(cellplace, buffplace, crd->natom, MPI_INT, MPI_MAX, 0,
	     CG->dspcomm);
  ReflectIVec(cellplace, buffplace, crd->natom);
  free(buffplace);
#endif
  if (CG->tid == 0) {
    for (i = 0; i < tp->natom; i++) {
      if (cellplace[i] == -1 && (tp->Masses[i] > 1.0e-8 || seekEP == 1)) {
	PrimeCCCErrMsg(CG, C, -1, i, tp, "Nonassigned atom", ErrMsg);
	printf("%s", ErrMsg);
      }
    }
  }
  if (halt == 1 && StopOnError == 1) {
    if (CG->tid == 0 && announce[0] != '\0') {
      printf("Exiting due to previous errors.\n");
    }
    exit(1);
  }

#ifdef MPI
  // Sync processes
  MPI_Barrier(CG->dspcomm);
#endif

  // Free allocated memory
  free(cellplace);
  free(ErrMsg);
}

#ifdef MPI
//-----------------------------------------------------------------------------
// PrintRecvInfo: print the information relating to a recv posting.  This
//                function is intended to speed the debugging process by
//                encapsulating this sort of report into a form that can be
//                added conveniently anywhere in the main code.
//
// Arguments:                                                           
//   maxsize:     the maximum size of the incoming message              
//   tname:       string defining the data type (for example, "INT")    
//   recver:      the receiving process                                 
//   sender:      the process from which to expect the message          
//   tag:         the message tag                                       
//   comm:        the communicator over which the recv passes           
//-----------------------------------------------------------------------------
void PrintRecvInfo(int maxsize, char* tname, int recver, int sender, int tag,
		   MPI_Comm comm)
{
  printf("Posting recv %4d -> %4d (%10d): Max  %6d of %s over %10d\n", sender,
	 recver, tag, maxsize, tname, (int)comm);
}

//-----------------------------------------------------------------------------
// PrintSendInfo: print the information relating to a send posting.  This
//                function is intended to speed the debugging process by
//                encapsulating this sort of report into a form that can be
//                added conveniently anywhere in the main code.
//
// Arguments:                                                           
//   actsize:     the actual size of the incoming message               
//   tname:       string defining the data type (for example, "INT")    
//   recver:      the receiving process                                 
//   sender:      the process from which to expect the message          
//   tag:         the message tag                                       
//   comm:        the communicator over which the recv passes           
//   data:        buffer containing the message                         
//   verbose:     flag to activate printing of ALL message data         
//-----------------------------------------------------------------------------
void PrintSendInfo(int actsize, char* tname, int recver, int sender, int tag,
		   MPI_Comm comm, void* data, int verbose)
{
  int i, outsize, typID;
  char *ctmp;
  char* output;

  // Considerations for printing ATOMB type data
  outsize = 128;
  if (strcmp(tname, "ATOMB") == 0) {
    outsize += actsize*48;
    typID = 1;
  }
  else if (strcmp(tname, "ATOMBV") == 0) {
    outsize += actsize*80;
    typID = 2;
  }
  else if (strcmp(tname, "ATOMBX") == 0) {
    outsize += actsize*220;
    typID = 3;
  }
  else {

    // Verbose printing is not supported for this data type
    verbose = 0;
  }
  output = (char*)malloc(outsize*sizeof(char));
  sprintf(output, "Posting send %4d -> %4d (%10d): Size %6d of %s over %10d\n",
	  sender, recver, tag, actsize, tname, (int)comm);
  if (verbose == 0) {

    // Print message and free the memory before returning
    printf("%s", output);
    free(output);

    return;
  }

  // If we're still here, print the entire send
  if (typID == 1) {
    atomb *rdata;
    rdata = (atomb*)data;

    // It is assumed that, if the message is one of these
    // custom-defined atom-related types, that the first 
    // element is not occupied with actual data.         
    for (i = 1; i < actsize+1; i++) {
      ctmp = &output[strlen(output)];
      sprintf(ctmp, "  %6d %4d %10.4lf %10.4lf %10.4lf\n", rdata[i].id,
	      rdata[i].dreg, rdata[i].loc[0], rdata[i].loc[1],
	      rdata[i].loc[2]);
    }
  }
  if (typID == 2) {
    atombv *rdata;
    rdata = (atombv*)data;

    // It is assumed that, if the message is one of these
    // custom-defined atom-related types, that the first 
    // element is not occupied with actual data.         
    for (i = 1; i < actsize+1; i++) {
      ctmp = &output[strlen(output)];
      sprintf(ctmp, "  %6d %4d %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf"
	      "\n", rdata[i].id, rdata[i].dreg, rdata[i].loc[0],
	      rdata[i].loc[1], rdata[i].loc[2], rdata[i].vel[0],
	      rdata[i].vel[1], rdata[i].vel[2]);
    }
  }
  printf("%s", output);

  // Free allocated memory
  free(output);
}
#endif

//-----------------------------------------------------------------------------
// Torque: compute the torque for a residue.                            
//                                                                      
// Arguments:                                                           
//   tp:     the topology structure                                     
//   crd:    the system coordinates                                     
//   resid:  the residue number                                         
//-----------------------------------------------------------------------------
static void Torque(prmtop *tp, coord *crd, int resid)
{
  int i, j;
  double com[3], r[3], rxf[3], trq[3];

  // Compute center of mass
  com[0] = com[1] = com[2] = 0.0;
  for (i = tp->ResLims[resid]; i < tp->ResLims[resid+1]; i++) {
    for (j = 0; j < 3; j++) {
      com[j] += tp->Masses[i]*crd->loc[3*i+j];
    }
  }

  // Compute torque
  trq[2] = trq[1] = trq[0] = 0.0;
  for (i = tp->ResLims[resid]; i < tp->ResLims[resid+1]; i++) {
    for (j = 0; j < 3; j++) {
      r[j] = crd->loc[3*i+j] - com[j];
    }
    CrossP(r, &crd->frc[3*i], rxf);
    for (j = 0; j < 3; j++) {
      trq[j] += rxf[j];
    }
  }

  printf("Torque = [ %12.8lf %12.8lf %12.8lf ]\n", trq[0], trq[1], trq[2]);
}

//-----------------------------------------------------------------------------
// WriteElimPair: write an elimination pair.  This condenses the code in
//                PrintResidueEliminations below.                       
//                                                                      
// Arguments:                                                           
//   tp:        the topology                                            
//   A:         integer 1, 2, 3, or 4                                   
//   nA:        the number of 1:A eliminations                          
//   listA:     the list of 1:A eliminations                            
//-----------------------------------------------------------------------------
static void WriteElimPair(prmtop *tp, int A, int nA, nixpr* listA)
{
  int i, xres, yres;

  for (i = 0; i < nA; i++) {
    xres = LocateResID(tp, listA[i].atmX, 0, tp->nres);
    yres = LocateResID(tp, listA[i].atmY, 0, tp->nres);
    printf("  1:%d -> %6d / %6d ( %.4s %.4s / %.4s %.4s )\n", A, listA[i].atmX,
	   listA[i].atmY, &tp->AtomNames[4*listA[i].atmX],
	   &tp->ResNames[4*xres], &tp->AtomNames[4*listA[i].atmY],
	   &tp->ResNames[4*yres]);
  }
}

//-----------------------------------------------------------------------------
// PrintResidueExclusions: print the exclusion list, atom by atom, for the
//                         specified residue.
//
// Arguments:                                                           
//   tp:        the topology                                            
//   rid:       the residue ID                                          
//-----------------------------------------------------------------------------
void PrintResidueExclusions(prmtop *tp, int rid)
{
  int i, j, llim, hlim, jatm, jres;

  llim = tp->ResLims[rid];
  hlim = tp->ResLims[rid+1];
  for (i = llim; i < hlim; i++) {

    // Print exclusions
    printf("Exclusions for atom %6d ( %.4s %.4s ):\n", i, &tp->AtomNames[4*i],
	   &tp->ResNames[4*rid]);
    for (j = tp->ConExcl[i]; j < tp->ConExcl[i+1]; j++) {
      jatm = tp->ExclList[j];

      // Bail out if there are no exclusions
      if (jatm == -1) {
	continue;
      }
      jres = LocateResID(tp, jatm, 0, tp->nres);
      printf("  %6d %.4s %.4s\n", jatm, &tp->AtomNames[4*jatm],
	     &tp->ResNames[4*jres]);
    }

    // Print eliminations
    if (tp->EPInserted == 1) {
      printf("Eliminations for atom %6d ( %.4s %.4s ):\n", i,
	     &tp->AtomNames[4*i], &tp->ResNames[4*rid]);
      WriteElimPair(tp, 1, tp->ElimPair[i].n11, tp->ElimPair[i].list11);
      WriteElimPair(tp, 2, tp->ElimPair[i].n12, tp->ElimPair[i].list12);
      WriteElimPair(tp, 3, tp->ElimPair[i].n13, tp->ElimPair[i].list13);
      WriteElimPair(tp, 4, tp->ElimPair[i].n14, tp->ElimPair[i].list14);
      printf("\n");
    }
  }
}

//-----------------------------------------------------------------------------
// PrintResidueForces: print forces on all atoms of a particular residue to
//                     stdout.  This function searches the cell grid primary
//                     sectors for all atoms.
//
// Arguments:                                                           
//   CG:        the cell grid                                           
//   crd:       coordinates struct                                      
//   tp:        the topology                                            
//   rid:       the residue ID                                          
//-----------------------------------------------------------------------------
void PrintResidueForces(cellgrid *CG, coord *crd, prmtop *tp, int rid)
{
  int h, i, atmcid;
  double *ftmp;
  cell *C;

  printf("%% Residue %d:\n", rid);
  printf("%%  Force X       Force Y       Force Z   "
	 "  Velocity X    Velocity Y    Velocity Z   Cell Index "
	 "  Position X    Position Y    Position Z    Charge\n"
	 "------------- ------------- ------------- "
	 "------------- ------------- ------------- --- --- --- "
	 "------------- ------------- ------------- ----------\n");
  for (h = tp->ResLims[rid]; h < tp->ResLims[rid+1]; h++) {
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      atmcid = FindAtomInSector(C, h, 0, 0);
      if (atmcid >= 0) {
	ftmp = C->data[atmcid].frc;
	printf("%13.7lf %13.7lf %13.7lf %13.7lf %13.7lf %13.7lf", ftmp[0],
	       ftmp[1], ftmp[2], crd->vel[3*h], crd->vel[3*h+1],
	       crd->vel[3*h+2]);
	ftmp = C->data[atmcid].loc;
	printf(" %3d %3d %3d %13.7lf %13.7lf %13.7lf %10.7lf\n", C->gbin[0],
	       C->gbin[1], C->gbin[2], ftmp[0], ftmp[1], ftmp[2],
	       C->data[atmcid].q);
	break;
      }
    }
  }
  printf("\n");
}

//-----------------------------------------------------------------------------
// PrintAtomState: print location, velocity, forces, and previous values of all
//                 these numbers for a particular atom.  This function searches
//                 the cell grid primary sectors for all atoms.
//
// Arguments:                                                           
//   CG:        the cell grid                                           
//   crd:       coordinates struct                                      
//   tp:        the topology                                            
//   atmid:     the residue ID                                          
//-----------------------------------------------------------------------------
void PrintAtomState(cellgrid *CG, coord *crd, int atmid)
{
  int i, atmcid;
  double *ftmp;
  cell *C;

  // Update the atom GPS
  UpdateCellGPS(CG);

  printf("%% Atom %d:\n", atmid);
  printf("%%  Force X       Force Y       Force Z   "
         "  Velocity X    Velocity Y    Velocity Z   Cell Index "
         "  Position X    Position Y    Position Z    Charge\n"
         "------------- ------------- ------------- "
         "------------- ------------- ------------- --- --- --- "
         "------------- ------------- ------------- ----------\n");
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    atmcid = C->GPSptr[atmid];
    if (atmcid >= 0) {
      ftmp = C->data[atmcid].frc;
      printf("%13.7lf %13.7lf %13.7lf %13.7lf %13.7lf %13.7lf", ftmp[0],
	     ftmp[1], ftmp[2], crd->vel[3*atmid], crd->vel[3*atmid+1],
	     crd->vel[3*atmid+2]);
      ftmp = C->data[atmcid].loc;
      printf(" %3d %3d %3d %13.7lf %13.7lf %13.7lf\n", C->gbin[0], C->gbin[1],
	     C->gbin[2], ftmp[0], ftmp[1], ftmp[2]);
      ftmp = &crd->prvfrc[3*atmid];
      printf("%13.7lf %13.7lf %13.7lf %13.7lf %13.7lf %13.7lf", ftmp[0],
	     ftmp[1], ftmp[2], crd->prvvel[3*atmid], crd->prvvel[3*atmid+1],
	     crd->prvvel[3*atmid+2]);
      ftmp = &crd->prvloc[3*atmid];
      printf(" %3d %3d %3d %13.7lf %13.7lf %13.7lf\n", C->gbin[0], C->gbin[1],
	     C->gbin[2], ftmp[0], ftmp[1], ftmp[2]);
      break;
    }
  }
  printf("\n");
}

//-----------------------------------------------------------------------------
// CellChecksum: perform a checksum operation on all cells, over various values
//               depending on the context.
//
// Arguments:                                                           
//   CG:        the cell grid                                           
//   crd:       the coordinates struct                                  
//   tp:        the topology (needed for particle masses)               
//   nsctr:     the number of sectors for which to print results (if    
//              this is set to 1, only the primary sector results will  
//              be displayed)                                           
//   do[p]loc:  flag to activate (previous) atom location sums, set to  
//              any value 1, 2, or 3 to get x, y, or z component sums of
//              positions                                               
//   dofrc:     flag to activate force sums                             
//   do[p]vel:  flag to activate (previous) atom velocity sums          
//   tagmsg:    message tag for this checksum operation                 
//   finmsg:    final message to print, indicating completion           
//   sleeptm:   time to sleep between cell printouts, to ensure proper  
//              ordering in the output                                  
//-----------------------------------------------------------------------------
void CellChecksum(cellgrid *CG, coord *crd, prmtop *tp, int nsctr, int doloc,
		  int dovel, int dofrc, int doploc, int dopvel, int dopfrc,
		  char* tagmsg, char* finmsg, double sleeptm)
{
  int i, j, k, atmid, olen;
  double chksum[6];
  char syscall[64];
  char *ctmp;
  char* outline;
  cell *C;

  // Allocate memory for output message
  outline = (char*)malloc(32 + (doloc + dovel + dofrc + doploc +
				dopvel + dopfrc) * nsctr * 21 * sizeof(char));

#ifdef MPI
  // Sync up all processes
  MPI_Barrier(CG->dspcomm);
#endif
  if (CG->tid == 0) {
    printf("%s\n", tagmsg);
  }
  for (i = 0; i < CG->ncell; i++) {
#ifdef MPI
    MPI_Barrier(CG->dspcomm);
#endif
    if (CG->data[i].CGRank != CG->tid) {
      continue;
    }
    C = &CG->data[i];

    // Print the checksums for all sectors
    sprintf(outline, "%4d ", C->gbin[3]);
    olen = strlen(outline);
    ctmp = &outline[olen];
    for (j = 0; j < nsctr; j++) {

      // Accumulate checksums
      for (k = 0; k < 6; k++) {
	chksum[k] = 0.0;
      }
      for (k = 0; k < C->nr[j]; k++) {
	atmid = C->map[j][k].id;

	// Current coordinates
	if (doloc) {
	  chksum[0] += C->map[j][k].loc[doloc-1];
	}

	// Current velocities
	if (dovel && tp->Masses[atmid] > 1.0e-8) {
	  chksum[1] += crd->vel[3*atmid+dovel-1];
	}

	// Current forces
	if (dofrc) {
	  chksum[2] += C->map[j][k].frc[dofrc-1];
	}

	// Previous coordinates
	if (doploc) {
	  chksum[3] += crd->prvloc[3*atmid+doploc-1];
	}

	// Previous velocities
	if (dopvel && tp->Masses[atmid] > 1.0e-8) {
	  chksum[4] += crd->prvvel[3*atmid+dopvel-1];
	}

	// Previous forces
	if (dopfrc && tp->Masses[atmid] > 1.0e-8) {
	  chksum[5] += crd->prvfrc[3*atmid+dopfrc-1];
	}
      }

      // Print the checksums
      if (doloc) {
	sprintf(ctmp, "%20.12e ", chksum[0]);
	olen = strlen(outline);
	ctmp = &outline[olen];
      }
      if (dovel) {
	sprintf(ctmp, "%20.12e ", chksum[1]);
	olen = strlen(outline);
	ctmp = &outline[olen];
      }
      if (dofrc) {
	sprintf(ctmp, "%20.12e ", chksum[2]);
	olen = strlen(outline);
	ctmp = &outline[olen];
      }
      if (doploc) {
	sprintf(ctmp, "%20.12e ", chksum[3]);
	olen = strlen(outline);
	ctmp = &outline[olen];
      }
      if (dopvel) {
	sprintf(ctmp, "%20.12e ", chksum[4]);
	olen = strlen(outline);
	ctmp = &outline[olen];
      }
      if (dopfrc) {
	sprintf(ctmp, "%20.12e ", chksum[5]);
	olen = strlen(outline);
	ctmp = &outline[olen];
      }
    }
    sprintf(ctmp, "\n");
    if (sleeptm > 1.0e-4) {
      sprintf(syscall, "sleep %8.4lf\n", sleeptm);
      system(syscall);
    }
    printf("%s", outline);
  }
  if (CG->tid == 0) {
    printf("%s\n", finmsg);
  }
  if (sleeptm > 1.0e-4) {
    sprintf(syscall, "sleep %8.4lf\n", sleeptm);
    system(syscall);
  }

  // Free allocated memory
  free(outline);
}

#ifdef MPI
//-----------------------------------------------------------------------------
// DisplaySend: display a send for a message.  Works for MPI messages   
//              passed by the cell grid communication plans.            
//                                                                      
// Arguments:                                                           
//   cshr:     sharing packet
//   isend:    number of the send message
//   offset:   offset of the memory index
//   msgsize:  size of the message
//   CG:       the cell grid
//   tj:       trajectory control data
//   msgtype:  type of MPI message (single integer, double array, etc.)
//-----------------------------------------------------------------------------
void DisplaySend(ashr *cshr, int isend, int offset, int msgsize, cellgrid *CG,
		 trajcon *tj, MPI_Datatype msgtype)
{
  int j, k, xpcon;
  char *ctmp;
  char* outline;

  // Decide on the appropriate size for the output message
  outline = (char*)malloc((200 + msgsize*20 + cshr->send[isend].ncell*10) *
			  sizeof(char));
  sprintf(outline, "%% Posting send for msgID %10d, thread %d -> %d over comm "
	  "%10d\n%%   AtomIDs (cell, atoms):\n",
	  cshr->send[isend].BaseID + offset, CG->tid,
	  cshr->send[isend].partner, CG->dspcomm);
  xpcon = 0;
  for (j = 0; j < cshr->send[isend].ncell; j++) {
    ctmp = &outline[strlen(outline)];
    sprintf(ctmp, "%%   %4d -> ", cshr->send[isend].cellpt[j]);
    if (msgtype == tj->MPI_ATOMX) {
      for (k = 1; k < CG->Xexport[isend][xpcon].id; k++) {
	ctmp = &outline[strlen(outline)];
	sprintf(ctmp, "%6d ", CG->Xexport[isend][xpcon+k].id);
      }
      ctmp = &outline[strlen(outline)];
      sprintf(ctmp, "\n");
    }
    else if (msgtype == tj->MPI_ATOMB) {
      for (k = 1; k < CG->pexport[isend][xpcon].id; k++) {
	ctmp = &outline[strlen(outline)];
	sprintf(ctmp, "%6d ", CG->pexport[isend][xpcon+k].id);
      }
      ctmp = &outline[strlen(outline)];
      sprintf(ctmp, "\n");
    }
  }
  ctmp = &outline[strlen(outline)];
  sprintf(ctmp, "%% Message = [ ");
  if (msgtype == tj->MPI_ATOMX) {
    for (j = 0; j < msgsize; j++) {
      ctmp = &outline[strlen(outline)];
      sprintf(ctmp, "%6d ", CG->Xexport[isend][j].id);
    }
  }
  else if (msgtype == tj->MPI_ATOMB) {
    for (j = 0; j < msgsize; j++) {
      ctmp = &outline[strlen(outline)];
      sprintf(ctmp, "%6d ", CG->pexport[isend][j].id);
    }
  }
  ctmp = &outline[strlen(outline)];
  sprintf(ctmp, "];\n");

  // Print the send information
  printf("%s\n", outline);

  // Free allocated memory
  free(outline);
}
#endif

//-----------------------------------------------------------------------------
// PrepCellGridCopy: copy an existing cell grid, creating a blank template for
//                   cloning the cell grid at a later time.
//                                                                      
// Arguments:                                                           
//   CG:      the original cell grid                                    
//-----------------------------------------------------------------------------
cellgrid PrepCellGridCopy(cellgrid *CG)
{
  int i, j, k, m;
  cellgrid nCG;

  for (i = 0; i < 3; i++) {
    nCG.ng[i] = CG->ng[i];
    nCG.dbng[i] = CG->dbng[i];
  }
  for (i = 0; i < 7; i++) {
    nCG.celldim[i] = CG->celldim[i];
  }
  nCG.ncell = CG->ncell;
  nCG.data = (cell*)malloc(CG->ng[0]*CG->ng[1]*CG->ng[2]*sizeof(cell));
  nCG.sysID = CG->sysID;
  nCG.maxatom = CG->maxatom;
  nCG.map = (cell***)malloc(nCG.ng[0]*sizeof(cell**));
  for (i = 0; i < nCG.ng[0]; i++) {
    nCG.map[i] = (cell**)malloc(nCG.ng[1]*sizeof(cell*));
    for (j = 0; j < nCG.ng[1]; j++) {
      nCG.map[i][j] = &nCG.data[(i*nCG.ng[1] + j)*nCG.ng[2]];
      for (k = 0; k < nCG.ng[2]; k++) {
	nCG.map[i][j][k] = CreateCell(nCG.maxatom, CG->data[0].pmordr);
	for (m = 0; m < 3; m++) {
	  nCG.map[i][j][k].orig[m] = CG->map[i][j][k].orig[m];
	}
	for (m = 0; m < 4; m++) {
	  nCG.map[i][j][k].gbin[m] = CG->map[i][j][k].gbin[3];
	}
      }
    }
  }

  return nCG;
}

//-----------------------------------------------------------------------------
// CopyCellGrid: copy a cell grid into another, pre-allocated chunk of memory.
//               Useful for storing snapshots of cell grids at different times
//               in a run.
//
// Arguments:                                                           
//   CY:      the cell grid copy to compose                             
//   CG:      the original cell grid                                    
//-----------------------------------------------------------------------------
void CopyCellGrid(cellgrid *CY, cellgrid *CG)
{
  int i, j, k;
  cell *ccY, *ccG;

  for (i = 0; i < CG->ncell; i++) {
    ccG = &CG->data[i];
    ccY = &CY->data[i];
    ccY->CGRank = ccG->CGRank;
  }
  for (i = 0; i < CG->MyCellCount; i++) {
    ccG = &CG->data[CG->MyCellDomain[i]];
    ccY = &CY->data[CG->MyCellDomain[i]];
    for (j = 0; j < 8; j++) {
      ccY->nr[j] = ccG->nr[j];
      ccY->nr[j+8] = ccG->nr[j+8];
      for (k = 0; k < ccG->nr[j] + ccG->nr[j+8]; k++) {
	ccY->map[j][k] = ccG->map[j][k];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Sector2Sector: this routine is the simplest but least efficient form of the
//                cell-based direct space computation.  Forces, energies, and  
//                virials can be computed in this function just as in the
//                optimized routines found in Nonbonded and its associated
//                libraries, but this single function takes flags for all of
//                these options rather than being split into multiple functions
//                by pre-processor directives.  This function makes no special
//                cases for van-der Waals and electrostatic cutoffs; it merely
//                assumes them to be independent.                       
//
// Arguments:                                                           
//   C:        the cell of interest                                     
//   s[A,B]:   interactions from sector A to sector B will be computed  
//   dcinp:    direct space control information                         
//   tp:       the topology                                             
//   dofrc:    flag to activate force computation                       
//   donrg:    flag to activate energy computation                      
//   dovir:    flag to activate virial computation                      
//-----------------------------------------------------------------------------
void Sector2Sector(cell *C, int sA, int sB, dircon *dcinp, prmtop *tp,
		   FrcTab *EFrc, Energy* sysUV, int dofrc, int donrg,
		   int dovir)
{
  int i, j, jmin, aiid, ailjt, ajljt;
  double fx, fy, fz, fmag, aifx, aify, aifz, uelec, uvdwr6, uvdwr12, aiq;
  double *ljftmp, *ljutmp;
  double **ljfmap, **ljumap;
  atomc *listA, *listB;
  CSpln *Uspl, *Fspl;

  // Pre-compute constants and unpack structures
  const double uspc = EFrc->ivdr;
  const double Vcut2 = dcinp->Vcut*dcinp->Vcut;
  const double Ecut2 = dcinp->Ecut*dcinp->Ecut;
  listA = C->map[sA];
  listB = C->map[sB];
  if (dofrc == 1) {
    ljfmap = tp->LJftab.map;
    Fspl = EFrc->dSD;
  }
  if (donrg == 1) {
    ljumap = tp->LJutab.map;
    Uspl = EFrc->SD;
    uelec = 0.0;
    uvdwr6 = 0.0;
    uvdwr12 = 0.0;
  }
  for (i = 0; i < C->nr[sA]; i++) {
    const double atmx = listA[i].loc[0];
    const double atmy = listA[i].loc[1];
    const double atmz = listA[i].loc[2];
    ailjt = listA[i].lj;
    jmin = (sA == sB) ? i+1 : 0;
    aiq = listA[i].q;
    aiid = listA[i].id;
    if (dofrc == 1) {
      aifx = 0.0;
      aify = 0.0;
      aifz = 0.0;
      if (ailjt >= 0) {
	ljftmp = ljfmap[ailjt];
      }
    }
    if (donrg == 1 && ailjt >= 0) {
      ljutmp = ljumap[ailjt];
    }
    for (j = jmin; j < C->nr[sB]; j++) {
      const double dx = listB[j].loc[0] - atmx;
      const double dy = listB[j].loc[1] - atmy;
      const double dz = listB[j].loc[2] - atmz;
      const double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < Ecut2) {
        const int ir2 = r2*uspc;
	if (dofrc == 1) {
	  fmag = (((Fspl[ir2].A*r2 + Fspl[ir2].B)*r2 + Fspl[ir2].C)*r2 +
		  Fspl[ir2].D)*aiq*listB[j].q;
	}
	if (donrg == 1) {
	  uelec += (((Uspl[ir2].A*r2 + Uspl[ir2].B)*r2 + Uspl[ir2].C)*r2 +
		    Uspl[ir2].D)*aiq*listB[j].q;
	}
      }
      else {
	fmag = 0.0;
      }
      ajljt = listB[j].lj;
      if (ailjt >= 0 && ajljt >= 0 && r2 < Vcut2 &&
          (r2 > MINNB2 || TestBondedExclusion(aiid, listB[j].id, tp) == 0)) {
        const double invr2 = 1.0/r2;
        const double invr4 = invr2*invr2;
	if (dofrc == 1) {
	  fmag += invr4*invr4*(ljftmp[2*ajljt]*invr4*invr2 +
			       ljftmp[2*ajljt+1]);
	}
	if (donrg == 1) {
	  uvdwr12 += invr4*invr2*invr4*invr2*ljutmp[2*ajljt];
	  uvdwr6 += invr4*invr2*ljutmp[2*ajljt+1];
	}
      }
      if (dofrc == 1) {
	fx = fmag*dx;
	fy = fmag*dy;
	fz = fmag*dz;
	listB[j].frc[0] -= fx;
	listB[j].frc[1] -= fy;
	listB[j].frc[2] -= fz;
	aifx += fx;
	aify += fy;
	aifz += fz;
      }
      if (dovir == 1) {
	sysUV->Vir[0] += fx*dx;
	sysUV->Vir[1] += fx*dy;
	sysUV->Vir[2] += fx*dz;
	sysUV->Vir[3] += fy*dx;
	sysUV->Vir[4] += fy*dy;
	sysUV->Vir[5] += fy*dz;
	sysUV->Vir[6] += fz*dx;
	sysUV->Vir[7] += fz*dy;
	sysUV->Vir[8] += fz*dz;
      }
    }

    // Accumulate the force
    if (dofrc == 1) {
      listA[i].frc[0] += aifx;
      listA[i].frc[1] += aify;
      listA[i].frc[2] += aifz;
    }
  }

  // Accumulate energy
  if (donrg == 1) {
    sysUV->vdw12 += uvdwr12;
    sysUV->vdw6 += uvdwr6;
    sysUV->delec += uelec;
  }
}

//-----------------------------------------------------------------------------
// ExhaustInteractions: compute all interactions of a system, given a set of
//                      coordinates.  This will accomplish what any cell-based
//                      method would, but without any spatial decomposition.
//                      Close interactions will be checked.  It is extremely
//                      slow.
//
// Arguments:
//   crd:     system coordinates
//   tp:      system topology
//   EFrc:    electrostatic force between any two particles, to be scaled by
//            the particles' respective charges
//   dcinp:   direct space input control data
//   sysUV:   system energy and virial accumulators
//   dofrc:   flag to have this function compute (direct space) forces
//   donrg:   flag to have this function compute (direct space) energies
//   dovir:   flag to have this compute direct space virial contributions
//   neartol: tolerance for close interactions (report if this is violated)
//-----------------------------------------------------------------------------
void ExhaustInteractions(coord *crd, prmtop *tp, FrcTab *EFrc, dircon *dcinp,
			 Energy *sysUV, int dofrc, int donrg, int dovir,
			 double neartol)
{
  int i, j, ailjt, ajljt;
  double fx, fy, fz, aifx, aify, aifz, fmag, uelec, uvdwr6, uvdwr12, aiq;
  double dx, dy, dz;
  double *ljftmp, *ljutmp;
  double **ljumap, **ljfmap;
  CSpln *Uspl, *Fspl;

  // The force array of the coordinate structure
  // will store the forces computed in this     
  // function.  Zero these force accumulators.  
  if (dofrc == 1) {
    SetDVec(crd->frc, 3*crd->natom, 0.0);
  }

  // Square the close interaction tolerance
  neartol *= neartol;
  
  // Pre-compute constants and unpack structures
  const double uspc = EFrc->ivdr;
  const double Vcut2 = dcinp->Vcut*dcinp->Vcut;
  const double Ecut2 = dcinp->Ecut*dcinp->Ecut;
  if (dofrc == 1) {
    ljfmap = tp->LJftab.map;
    Fspl = EFrc->dSD;
  }
  if (donrg == 1) {
    ljumap = tp->LJutab.map;
    Uspl = EFrc->SD;
    uelec = 0.0;
    uvdwr6 = 0.0;
    uvdwr12 = 0.0;
  }
  for (i = 0; i < crd->natom-1; i++) {
    const double atmx = crd->loc[3*i];
    const double atmy = crd->loc[3*i+1];
    const double atmz = crd->loc[3*i+2];
    aiq = tp->Charges[i];
    ailjt = tp->LJIdx[i];
    if (dofrc == 1) {
      aifx = 0.0;
      aify = 0.0;
      aifz = 0.0;
      if (ailjt >= 0) {
        ljftmp = ljfmap[ailjt];
      }
    }
    if (donrg == 1 && ailjt >= 0) {
      ljutmp = ljumap[ailjt];
    }
    for (j = i+1; j < crd->natom; j++) {
      dx = crd->loc[3*j]-atmx;
      dy = crd->loc[3*j+1]-atmy;
      dz = crd->loc[3*j+2]-atmz;
      NonOrthoReim(&dx, &dy, &dz, &crd->U, &crd->invU);
      const double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < neartol && TestBondedExclusion(i, j, tp) == 0) {
	printf("ExhaustInteractions >> Near interaction between atoms "
	       "(%9.6lf)\nExhaustInteractions >>   %4.4s %5d %4.4s %6d and "
	       "%4.4s %5d %4.4s %6d\n", sqrt(r2),
	       &tp->ResNames[4*LocateResID(tp, i, 0, tp->nres)],
	       LocateResID(tp, i, 0, tp->nres), &tp->AtomNames[4*i], i,
	       &tp->ResNames[4*LocateResID(tp, j, 0, tp->nres)],
	       LocateResID(tp, j, 0, tp->nres), &tp->AtomNames[4*j], j);
      }
      if (r2 < Ecut2) {
        const int ir2 = r2*uspc;
	if (dofrc == 1) {
          fmag = (((Fspl[ir2].A*r2 + Fspl[ir2].B)*r2 + Fspl[ir2].C)*r2 +
                  Fspl[ir2].D)*aiq*tp->Charges[j];
        }
        if (donrg == 1) {
          uelec += (((Uspl[ir2].A*r2 + Uspl[ir2].B)*r2 + Uspl[ir2].C)*r2 +
                    Uspl[ir2].D)*aiq*tp->Charges[j];
        }
      }
      else {
        fmag = 0.0;
      }
      ajljt = tp->LJIdx[j];
      if (ailjt >= 0 && ajljt >= 0 && r2 < Vcut2 &&
          (r2 > MINNB2 || TestBondedExclusion(i, j, tp) == 0)) {
        const double invr2 = 1.0/r2;
	const double invr4 = invr2*invr2;
        if (dofrc == 1) {
          fmag += invr4*invr4*(ljftmp[2*ajljt]*invr4*invr2 +
                               ljftmp[2*ajljt+1]);
        }
        if (donrg == 1) {
          uvdwr12 += invr4*invr2*invr4*invr2*ljutmp[2*ajljt];
          uvdwr6 += invr4*invr2*ljutmp[2*ajljt+1];
        }
      }
      if (dofrc == 1) {
        fx = fmag*dx;
        fy = fmag*dy;
        fz = fmag*dz;
        crd->frc[3*j] -= fx;
        crd->frc[3*j+1] -= fy;
        crd->frc[3*j+2] -= fz;
        aifx += fx;
	aify += fy;
        aifz += fz;
      }
      if (dovir == 1) {
        sysUV->Vir[0] += fx*dx;
        sysUV->Vir[1] += fx*dy;
        sysUV->Vir[2] += fx*dz;
        sysUV->Vir[3] += fy*dx;
        sysUV->Vir[4] += fy*dy;
        sysUV->Vir[5] += fy*dz;
        sysUV->Vir[6] += fz*dx;
        sysUV->Vir[7] += fz*dy;
        sysUV->Vir[8] += fz*dz;
      }
    }

    // Accumulate the force
    if (dofrc == 1) {
      crd->frc[3*i] += aifx;
      crd->frc[3*i+1] += aify;
      crd->frc[3*i+2] += aifz;
    }
  }

  // Accumulate energy
  if (donrg == 1) {
    sysUV->vdw12 = uvdwr12;
    sysUV->vdw6 = uvdwr6;
    sysUV->delec = uelec;
  }
}

//-----------------------------------------------------------------------------
// ListBondedIntr: recurrecting an old function whereby bonded interactions
//                 are calculated using the list-based atom format.
//
// Arguments:
//   crd:     system coordinates
//   tp:      system topology
//   Cfrc:    electrostatic force between any two particles, to be scaled by
//            the particles' respective charges
//   Hfrc:    high-resolution short-ranged form of EFrc
//   sysUV:   system energy and virial accumulators
//   qform:   the format of the charges (point or Gaussian)
//-----------------------------------------------------------------------------
void ListBondedIntr(coord *crd, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
		    Energy *sysUV, int qform)
{
  int atmid, i, aid, bid, cid, did, eid;
  double *atmcrd, *atmfrc;
  bondlist *tblc;
  bondcomm *tbcom;
  angllist *talc;
  anglcomm *tacom;
  dihelist *thlc;
  dihecomm *thcom;
  cmaplist *tmlc;
  cmapcomm *tmcom;

  // Loop over all atoms
  for (atmid = 0; atmid < tp->natom; atmid++) {

    // Set pointers for the atom which may control bonded interactions
    atmcrd = &crd->loc[3*atmid];
    atmfrc = &crd->frc[3*atmid];
    
    // Bonds
    if (tp->BLC[atmid].nbond > 0) {
      tblc = &tp->BLC[atmid];
      for (i = 0; i < tblc->nbond; i++) {
        tbcom = &tblc->BC[i];
        bid = tbcom->b;
        if (sysUV->updateU == 2) {
          BondVir(atmcrd, &crd->loc[3*bid], atmfrc, &crd->frc[3*bid],
                  tbcom, tp, Cfrc, Hfrc, sysUV, qform);
        }
        else if (sysUV->updateU == 1) {
          BondFrcNrg(atmcrd, &crd->loc[3*bid], atmfrc,
                     &crd->frc[3*bid], tbcom, tp, Cfrc, Hfrc, sysUV, qform);
        }
        else if (sysUV->updateU == -1) {
          BondNrg(atmcrd, &crd->loc[3*bid], tbcom, tp, Cfrc, Hfrc,
                  sysUV, qform);
        }
        else {
          BondFrc(atmcrd, &crd->loc[3*bid], atmfrc, &crd->frc[3*bid],
                  tbcom, tp, Cfrc, Hfrc, qform);
        }
      }
    }

    // Angles
    if (tp->ALC[atmid].nangl > 0) {
      talc = &tp->ALC[atmid];
      for (i = 0; i < talc->nangl; i++) {
        tacom = &talc->AC[i];
        aid = tacom->a;
        cid = tacom->c;
        if (sysUV->updateU == 2) {
          AnglVir(&crd->loc[3*aid], atmcrd, &crd->loc[3*cid],
                  &crd->frc[3*aid], atmfrc, &crd->frc[3*cid],
                  tacom, tp, Cfrc, Hfrc, sysUV, qform);
        }
        else if (sysUV->updateU == 1) {
          AnglFrcNrg(&crd->loc[3*aid], atmcrd, &crd->loc[3*cid],
                     &crd->frc[3*aid], atmfrc, &crd->frc[3*cid], tacom, tp,
                     Cfrc, Hfrc, sysUV, qform);
        }
        else if (sysUV->updateU == -1) {
          AnglNrg(&crd->loc[3*aid], atmcrd, &crd->loc[3*cid],
                  tacom, tp, Cfrc, Hfrc, sysUV, qform);
        }
        else {
          AnglFrc(&crd->loc[3*aid], atmcrd, &crd->loc[3*cid],
                  &crd->frc[3*aid], atmfrc, &crd->frc[3*cid],
                  tacom, tp, Cfrc, Hfrc, qform);
        }
      }
    }

    // Dihedrals
    if (tp->HLC[atmid].ndihe > 0) {
      thlc = &tp->HLC[atmid];
      for (i = 0; i < thlc->ndihe; i++) {
        thcom = &thlc->HC[i];
        aid = thcom->a;
        cid = thcom->c;
        did = thcom->d;
        if (sysUV->updateU == 2) {
          DiheVir(&crd->loc[3*aid], atmcrd, &crd->loc[3*cid],
                  &crd->loc[3*did], &crd->frc[3*aid], atmfrc,
                  &crd->frc[3*cid], &crd->frc[3*did], thcom,
                  tp, Cfrc, Hfrc, sysUV, qform);
        }
        else if (sysUV->updateU == 1) {
          DiheFrcNrg(&crd->loc[3*aid], atmcrd, &crd->loc[3*cid],
                     &crd->loc[3*did], &crd->frc[3*aid], atmfrc,
                     &crd->frc[3*cid], &crd->frc[3*did], thcom,
                     tp, Cfrc, Hfrc, sysUV, qform);
        }
        else if (sysUV->updateU == -1) {
          DiheNrg(&crd->loc[3*aid], atmcrd, &crd->loc[3*cid],
                  &crd->loc[3*did], thcom, tp, Cfrc, Hfrc, sysUV, qform);
        }
        else if (sysUV->updateU == 0) {
          DiheFrc(&crd->loc[3*aid], atmcrd, &crd->loc[3*cid],
                  &crd->loc[3*did], &crd->frc[3*aid], atmfrc,
                  &crd->frc[3*cid], &crd->frc[3*did], thcom,
                  tp, Cfrc, Hfrc, qform);
        }
      }
    }

    // CMAPs
    if (tp->MLC[atmid].ncmap > 0) {
      tmlc = &tp->MLC[atmid];
      for (i = 0; i < tmlc->ncmap; i++) {
        tmcom = &tmlc->MC[i];
        aid = tmcom->a;
        bid = tmcom->b;
        did = tmcom->d;
        eid = tmcom->e;
        if (sysUV->updateU == 2) {
          CmapVir(&crd->loc[3*aid], &crd->loc[3*bid], atmcrd,
                  &crd->loc[3*did], &crd->loc[3*eid], &crd->frc[3*aid],
                  &crd->frc[3*bid], atmfrc, &crd->frc[3*did],
                  &crd->frc[3*eid], tmcom, tp, sysUV, qform);
        }
        else if (sysUV->updateU == 1) {
          CmapFrcNrg(&crd->loc[3*aid], &crd->loc[3*bid], atmcrd,
                     &crd->loc[3*did], &crd->loc[3*eid], &crd->frc[3*aid],
                     &crd->frc[3*bid], atmfrc, &crd->frc[3*did],
                     &crd->frc[3*eid], tmcom, tp, sysUV, qform);
        }
        else if (sysUV->updateU == -1) {
          CmapNrg(&crd->loc[3*aid], &crd->loc[3*bid], atmcrd,
                  &crd->loc[3*did], &crd->loc[3*eid], tmcom, tp, sysUV, qform);
        }
        else if (sysUV->updateU == 0) {
          CmapFrc(&crd->loc[3*aid], &crd->loc[3*bid], atmcrd,
                  &crd->loc[3*did], &crd->loc[3*eid], &crd->frc[3*aid],
                  &crd->frc[3*bid], atmfrc, &crd->frc[3*did],
                  &crd->frc[3*eid], tmcom, tp, qform);
        }
      }
    }
  }
}
