#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Random.h"
#include "mdgxVector.h"
#include "Thermostats.h"
#include "Barostats.h"
#include "CellManip.h"
#include "pmeRecip.h"
#include "CrdManip.h"
#include "Integrator.h"
#include "Macros.h"
#include "Trajectory.h"
#include "Constraints.h"
#include "Debug.h"

#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "pmeDirectDS.h"

//-----------------------------------------------------------------------------
// ScaleCell: this routine scales coordinates of atoms within a particular
//            cell by a factor of chi, which can be a three-element vector.
//            This routine scales all coordinates independently, so if there
//            are constraint groups involved their coordinates must be
//            prepared by PrepScaleCnst beforehand.
//
// Arguments:
//   C:       the cell to scale
//   crd:     the coordinates of the system
//   tp:      systemp topology
//   chi:     factor by which to scale the coordinates of atoms in the cell C
//-----------------------------------------------------------------------------
static void ScaleCell(cell *C, coord *crd, prmtop *tp, double* chi)
{
  int i;
  double *ltmp, *sltmp;

  for (i = 0; i < C->nr[0]; i++) {
    if (tp->MobileAtoms[C->data[i].id] == 0) {
      continue;
    }
    ltmp = C->data[i].loc;
    sltmp = &crd->loc[3*C->data[i].id];
    ltmp[0] *= chi[0];
    sltmp[0] *= chi[0];
    ltmp[1] *= chi[1];
    sltmp[1] *= chi[1];
    ltmp[2] *= chi[2];
    sltmp[2] *= chi[2];
  }
}

#ifdef MPI
//-----------------------------------------------------------------------------
// CountCGMsg: count the number of messages that this process must expect to
//             receive from others or send to others, based on the overlap of
//             cells in the grid.  In the first call, the function looks for
//             overlaps between cells of the original grid and the new one,
//             but in the second call it looks for the reverse to determine
//             receives and sends respectively.
//
// Arguments:                                                           
//   CG:    the original cell grid that has sent atoms to CGN           
//   CGN:   the newly created cell grid that has received atoms from CG,
//          and will now reorganize them                                
//   tj:    trajectory control information (also contains process rank) 
//   rgd:   the reference grid, the one whose atom counts are taken for 
//          buffer sizes                                                
//-----------------------------------------------------------------------------
static lgrp CountCGMsg(cellgrid *CG, cellgrid *CGN, prmtop *tp, int rgd)
{
  int i, m, minx, maxx, miny, maxy, minz, maxz, maxmsg, maxatom;
  int ci, cj, ck, cpuid, msgfound;
  lgrp lmsg;
  cell *Ccgn;

  // This routine makes use of the lgrp struct (see
  // the Topology library's data structures) but   
  // rather than a group of atoms it describes a   
  // group of messages, with the natom field giving
  // the number of messages and the atoms field    
  // keeping a list of descriptors for each message
  // [ process message partner ][ # of cells ]     
  const double ratiox = (double)(CG->ng[0])/(double)(CGN->ng[0]);
  const double ratioy = (double)(CG->ng[1])/(double)(CGN->ng[1]);
  const double ratioz = (double)(CG->ng[2])/(double)(CGN->ng[2]);
  const int CGtid = CG->tid;
  maxmsg = 8;
  lmsg.natom = 0;
  lmsg.atoms = (int*)malloc(2*maxmsg*sizeof(int));
  for (i = 0; i < CGN->MyCellCount; i++) {
    Ccgn = &CGN->data[CGN->MyCellDomain[i]];
    minx = Ccgn->gbin[0]*ratiox + 1.0e-8;
    maxx = (Ccgn->gbin[0]+1)*ratiox + 1.0e-8;
    miny = Ccgn->gbin[1]*ratioy + 1.0e-8;
    maxy = (Ccgn->gbin[1]+1)*ratioy + 1.0e-8;
    minz = Ccgn->gbin[2]*ratioz + 1.0e-8;
    maxz = (Ccgn->gbin[2]+1)*ratioz + 1.0e-8;
    if (maxx == CG->ng[0]) maxx = CG->ng[0]-1;
    if (maxy == CG->ng[1]) maxy = CG->ng[1]-1;
    if (maxz == CG->ng[2]) maxz = CG->ng[2]-1;
    if (fabs(ratiox-1.0) < 1.0e-8) maxx = minx;
    if (fabs(ratioy-1.0) < 1.0e-8) maxy = miny;
    if (fabs(ratioz-1.0) < 1.0e-8) maxz = minz;

    // Test the cells of CG that this cell of
    // CGN overlaps based on cell limits.  If
    // they are owned by different processes,
    // then there is a message to receive    
    // from another process.                 
    for (ci = minx; ci <= maxx; ci++) {
      for (cj = miny; cj <= maxy; cj++) {
	for (ck = minz; ck <= maxz; ck++) {
	  if (CG->map[ci][cj][ck].CGRank == CGtid) {
	    continue;
	  }
	  cpuid = CG->map[ci][cj][ck].CGRank;
	  msgfound = 0;
	  for (m = 0; m < lmsg.natom; m++) {
	    if (lmsg.atoms[2*m] == cpuid) {
	      lmsg.atoms[2*m+1] += 1;
	      msgfound = 1;
	    }
	  }
	  if (msgfound == 0) {
	    lmsg.atoms[2*lmsg.natom] = cpuid;
	    lmsg.atoms[2*lmsg.natom+1] = 1;
	    lmsg.natom += 1;
	    if (lmsg.natom == maxmsg) {
	      maxmsg += 8;
	      lmsg.atoms = (int*)realloc(lmsg.atoms, 2*maxmsg*sizeof(int));
	    }
	  }
	}
      }
    }
  }

  // Translate the cell counts into maximum atom buffer sizes
  maxatom = (rgd == 1) ? CG->maxatom : CGN->maxatom;
  for (i = 0; i < lmsg.natom; i++) {
    lmsg.atoms[2*i+1] = MIN(lmsg.atoms[2*i+1]*maxatom, tp->natom);
  }

  return lmsg;
}

//-----------------------------------------------------------------------------
// PostAllMsg: post recvs sends for communicating the domain of one process on
//             the current cell grid to other processes on a new cell grid.
//
// Arguments:                                                           
//   lrecv:   the list of receives                                      
//   lsend:   the list of sends                                         
//   tj:      trajectory control data (process rank and MPI data types) 
//   CGN:     the new cell grid (contains atoms from cells that this process
//            originally controlled, spread into primary sectons of its cells)
//   crd:     coordinates                                               
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
static void PostAllMsg(lgrp *lrecv, lgrp *lsend, trajcon *tj, cellgrid *CGN,
		       coord *crd, prmtop *tp)
{
  int i, j, k, a3k, nreq, atmid, natm;
  atombx** sbuff;
  atombx** rbuff;
  cell *C;
  atombx *atmX, *tmpX;
  MPI_Request* req;
  MPI_Status* stt;

  i = lrecv->natom+lsend->natom;
  req = (MPI_Request*)malloc(i*sizeof(MPI_Request));
  stt = (MPI_Status*)malloc(i*sizeof(MPI_Status));

  // Allocate memory and post receives
  rbuff = (atombx**)malloc((lrecv->natom+1)*sizeof(atombx*));
  for (i = 0; i < lrecv->natom; i++) {
    rbuff[i] = (atombx*)malloc(lrecv->atoms[2*i+1]*sizeof(atombx));
    MPI_Irecv(rbuff[i], lrecv->atoms[2*i+1], tj->MPI_ATOMX, lrecv->atoms[2*i],
	      6*(CGN->nthreads*CGN->tid + lrecv->atoms[2*i])*DSP_MSG_TYPES +
	      DSP_RESCALE, CGN->dspcomm, &req[i]);
  }
  nreq = lrecv->natom;

  // Allocate memory for sends
  sbuff = (atombx**)malloc((lsend->natom+1)*sizeof(atombx*));
  for (i = 0; i < lsend->natom; i++) {
    sbuff[i] = (atombx*)malloc(lsend->atoms[2*i+1]*sizeof(atombx));

    // As usual, the first "atom" is a blank
    // for storing the size of the send     
    sbuff[i][0].id = 1;
  }

  // Fill send buffers
  for (i = 0; i < CGN->ncell; i++) {
    if (CGN->data[i].CGRank != CGN->tid && CGN->data[i].nr[0] > 0) {
      C = &CGN->data[i];
      for (j = 0; j < lsend->natom; j++) {
	if (lsend->atoms[2*j] == C->CGRank) {
	  atmX = sbuff[j];
	  break;
	}
      }
      for (j = 0; j < C->nr[0]; j++) {
	tmpX = &atmX[atmX[0].id];
	atmid = C->data[j].id;
	tmpX->id = atmid;
	tmpX->dreg = C->gbin[3];
	for (k = 0; k < 3; k++) {
	  a3k = 3*atmid+k;
	  tmpX->loc[k] = C->data[j].loc[k];
	  tmpX->vel[k] = crd->vel[a3k];
	  tmpX->prvvel[k] = crd->prvvel[a3k];
	  tmpX->prvfrc[k] = crd->prvfrc[a3k];
	  tmpX->sysloc[k] = crd->loc[a3k];
	  tmpX->prvloc[k] = crd->prvloc[a3k];
	}
	atmX[0].id += 1;
      }
    }
  }
  for (i = 0; i < lsend->natom; i++) {
    MPI_Isend(sbuff[i], sbuff[i][0].id, tj->MPI_ATOMX, lsend->atoms[2*i],
	      6*(CGN->nthreads*lsend->atoms[2*i] + CGN->tid)*DSP_MSG_TYPES +
              DSP_RESCALE, CGN->dspcomm, &req[nreq]);
    nreq++;
  }

  // Wait for all messages to complete
  MPI_Waitall(nreq, req, stt);
  MPI_Barrier(CGN->dspcomm);

  // Unpack
  for (i = 0; i < lrecv->natom; i++) {
    atmX = rbuff[i];
    for (j = 1; j < atmX[0].id; j++) {
      tmpX = &atmX[j];
      C = &CGN->data[tmpX->dreg];
      natm = C->nr[0];
      atmid = tmpX->id;
      C->data[natm].id = atmid;
      C->data[natm].q = tp->Charges[atmid];
      C->data[natm].lj = tp->LJIdx[atmid];
      for (k = 0; k < 3; k++) {
	a3k = 3*atmid+k;
	C->data[natm].loc[k] = tmpX->loc[k];
	crd->vel[a3k] = tmpX->vel[k];
	crd->prvvel[a3k] = tmpX->prvvel[k];
	crd->prvfrc[a3k] = tmpX->prvfrc[k];
	crd->loc[a3k] = tmpX->sysloc[k];
	crd->prvloc[a3k] = tmpX->prvloc[k];
      }
      C->nr[0] += 1;
    }
  }

  // Free allocated memory
  for (i = 0; i < lrecv->natom; i++) {
    free(rbuff[i]);
  }
  free(rbuff);
  for (i = 0; i < lsend->natom; i++) {
    free(sbuff[i]);
  }
  free(sbuff);
  free(req);
  free(stt);
}
#endif

//-----------------------------------------------------------------------------
// RepartitionCellGrid: this function, which should only be called once in a
//                      blue moon, will change the way in which a cell grid is 
//                      partitioned.  It takes as a full cell grid and returns
//                      a new one containing the same data with a
//                      repartitioned scheme.          
//
// Arguments:                                                           
//   CG:      the original cell grid                                    
//   crd:     the coordinates in a unified list                         
//   dcinp:   direct-space control data (contains, the nonbonded cutoff 
//            which determines the minimum cell spacing)                
//   rcinp:   reciprocal space control data                             
//-----------------------------------------------------------------------------
static cellgrid RepartitionCellGrid(cellgrid *CG, coord *crd, dircon *dcinp,
				    reccon *rcinp, prmtop *tp, trajcon *tj)
{
  int i, j, cx, cy, cz;
  double invcelldim[3];
  double *utmp, *loctmp;
  cell *C, *Cn;
  cellgrid CGN;
#ifdef MPI
  lgrp lsend, lrecv;
  MPI_Comm tcomm;
#endif

  // The new cell grid will be created based on the 
  // current unit cell dimensions as found in crd.  
  // The old cell grid knows about the new unit cell
  // dimensions and all coordinates within it are   
  // correct, but the cells are either overgrown or 
  // too small for the direct-space cutoff.         
  CGN = CreateCellGrid(crd, dcinp, rcinp, tp, tj, CG->sysID);
#ifdef MPI
  tcomm = CGN.dspcomm;
  CGN.dspcomm = CG->dspcomm;
#endif
  for (i = 0; i < 3; i++) {
    invcelldim[i] = (crd->isortho == 1) ? 1.0/CGN.celldim[i] : CGN.ng[i];
  }

  // The two cell grids will work on the same processors,
  // but have different communicators.  This will require
  // communication between processes of the new grid's   
  // communicator only, as atoms in cells of the first   
  // grid owned by this process are transferred directly 
  // into cells of the new grid, irrespective of whether 
  // this process owns them.  It is then up to processes 
  // to talk over the second grid's communicator to move 
  // atoms around appropriately.                         
#ifdef MPI
  LinkCellGrid(&CGN, crd, rcinp);
#else
  LinkCellGrid(&CGN, rcinp);
#endif

  // Loop over all cells in the old grid and prepare
  // to ship them off to the new grid.  This routine
  // only considers atoms in the primary sectors of 
  // each cell, and ships them off appropriately.   
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    if (crd->isortho == 1) {
      for (j = 0; j < C->nr[0]; j++) {
	cx = C->data[j].loc[0]*invcelldim[0];
	cy = C->data[j].loc[1]*invcelldim[1];
	cz = C->data[j].loc[2]*invcelldim[2];
	Cn = &CGN.map[cx][cy][cz];
	Cn->data[Cn->nr[0]] = C->data[j];
	Cn->nr[0] += 1;
      }
    }
    else {
      for (j = 0; j < C->nr[0]; j++) {
	utmp = crd->U.data;
	loctmp = C->data[j].loc;

	// Here, invcelldim is not an inverse cell dimension
	// in the sense it is for the orthorhombic case, but
	// it does keep us from having to promote integers  
	// to doubles again and again.                      
	cx = (utmp[0]*loctmp[0] + utmp[1]*loctmp[1] + utmp[2]*loctmp[2]) *
	  invcelldim[0];
	cy = (utmp[3]*loctmp[0] + utmp[4]*loctmp[1] + utmp[5]*loctmp[2]) *
	  invcelldim[1];
	cz = (utmp[6]*loctmp[0] + utmp[7]*loctmp[1] + utmp[8]*loctmp[2]) *
	  invcelldim[2];
	Cn = &CGN.map[cx][cy][cz];
	Cn->data[Cn->nr[0]] = C->data[j];
	Cn->nr[0] += 1;
      }
    }
  }

  // At this point, each process will have mapped all   
  // cells which it owns in CG to its version of CGN.   
  // The various versions of CGN in each process must be
  // synchronized by process-to-process communication.  
#ifdef MPI
  if (CG->nthreads > 1) {
    lrecv = CountCGMsg(CG, &CGN, tp, 0);
    lsend = CountCGMsg(&CGN, CG, tp, 1);
    PostAllMsg(&lrecv, &lsend, tj, &CGN, crd, tp);
    free(lrecv.atoms);
    free(lsend.atoms);
  }
#endif

  // Free the original cell grid
#ifdef MPI
  CG->dspcomm = tcomm;
#endif
  DestroyCellGrid(CG);

  // The atoms of the new cell grid must now be re-ordered.
  for (i = 0; i < CGN.MyCellCount; i++) {
    Cn = &CGN.data[CGN.MyCellDomain[i]];
    qsort(Cn->data, Cn->nr[0], sizeof(atomc), SortAtomID);
  }

  return CGN;
}

//-----------------------------------------------------------------------------
// GridScaleCoord: this routine adjusts all system coordinates, taking into
//                 account rigidly constrained groups (rigidly constrained     
//                 groups are moved to rescale their centers of mass).
//                 Because constraints are again involved, there can be
//                 additional cell-to-cell communication.  In rare cases,
//                 rescaling may require that the cell grid dimensions change.
//                 If this happens, a new cell grid is allocated.  Coordinates
//                 from the old cell grid are always updated within the old
//                 cell grid, and then transferred to the new one if
//                 necessary.  WARNING: atoms in the primary sectors of cells
//                 emerge from this routine with their positions rescaled, but
//                 atoms in import regions do not.        
//
// Arguments:
//   crd:      the coordinates of the system
//   tp:       system topology
//   chi:      scaling factor to apply across the entire system, respecting
//             rigid molecules whose centers of mass get scaled instead
//   dcinp:    the direct space input control data
//   rcinp:    the reciprocal space input control data
//   PPk:      reciprocal space pair potential kit for grid-based interactions
//   CG:       the cell grid for decomposing direct space interactions
//   tj:       trajectory control data
//-----------------------------------------------------------------------------
void GridScaleCoord(coord *crd, prmtop *tp, double* chi, dircon *dcinp,
		    reccon *rcinp, bckit *PPk, cellgrid *CG,
		    trajcon *tj)
{
  int i, neednewgrid;
  int cdim[3];
  double cdepth[3];
  
  // If all coordinates stretch independently (i.e. there are 
  // no constraints such that every atom position is rescaled 
  // by the same factor) then the cell grid composition WILL  
  // NOT CHANGE during a volume rescaling operation.  However,
  // if there are constraints on the system then the various  
  // cells' atom contents may change during the rescaling.    
  // The approach here is to prepare the system so that the   
  // positions of all particles can be scaled independently.  
  // This is done by pre-adjusting the positions of all atoms 
  // in constraint groups (relative to the center of mass of  
  // the constraint group) such that they will end up in the  
  // right places when multiplied by the rescaling factor chi,
  // the performing the same cell-to-cell communication as is 
  // done for other position constraints ApplyGridCnst()      
  // function in the Constraints library.                     
  ApplyGridCnst(CG, crd, tp, tj, dcinp->Mcut, 2, chi);

  // Loop over all cells and scale positions of
  // particles in each cell's primary sector   
  for (i = 0; i < CG->MyCellCount; i++) {
    ScaleCell(&CG->data[CG->MyCellDomain[i]], crd, tp, chi);
  }

  // Update unit cell dimensions within the coord struct.
  for (i = 0; i < 3; i++) {
    crd->gdim[i] *= chi[i];
    crd->hgdim[i] = 0.5*crd->gdim[i];
  }
  CompXfrm(crd->gdim, crd->U, crd->invU);

  // Update the PME machinery to compensate for the new box size
  UpdateBCKit(rcinp, crd, tp, PPk);

  // Check the cell grid dimensions (note that the
  // TakeCellGridDims() function doesn't actually 
  // change anythign in CG--it just returns a lot 
  // of information that can be loaded into CG or 
  // another cell grid.                           
  TakeCellGridDims(cdim, cdepth, crd, dcinp);
  for (i = 0; i < 3; i++) {
    CG->celldim[i] = crd->gdim[i]/CG->ng[i];
    CG->celldim[4+i] = cdepth[i];
  }
  ComputeCellOrigins(CG, crd);

  // Now it is time to decide whether to keep the
  // old cell grid or make an entirely new one   
  neednewgrid = 0;
  for (i = 0; i < 3; i++) {
    if (cdim[i] != CG->ng[i]) {
      neednewgrid = 1;
    }
  }
  if (neednewgrid == 1) {

    // Allocate space for a new cell   
    // grid with the correct dimensions
    *CG = RepartitionCellGrid(CG, crd, dcinp, rcinp, tp, tj);
  }
}

//-----------------------------------------------------------------------------
// PrepMCBarostat: prepare volume rescaling factors for a Monte-Carlo   
//                 barostat.
//
// Arguments:                                                           
//   tj:     the trajectory control data                                
//   crd:    the system coordinates                                     
//-----------------------------------------------------------------------------
void PrepMCBarostat(trajcon *tj, coord *crd)
{
  int i;

  if ((tj->ntp == 1 && tj->MCBarostatFac[0] < 0.0) ||
      (tj->ntp == 2 && tj->MCBarostatFac[0] < 0.0 &&
       tj->MCBarostatFac[1] < 0.0 && tj->MCBarostatFac[2] < 0.0)) {
    printf("PrepMCBarostat >> Error.  Bad input to Monte Carlo barostat.\n"
	   "PrepMCBarostat >>   ntp = %d, mccomp = [ %7.4lf %7.4lf %7.4lf ]\n",
	   tj->ntp, tj->MCBarostatFac[0], tj->MCBarostatFac[1],
	   tj->MCBarostatFac[2]);
    exit(1);
  }

  // If one of the dimensions is unset, set it to zero
  for (i = 0; i < 3; i++) {
    if (tj->MCBarostatFac[i] < 0.0) {
      tj->MCBarostatFac[i] = 0.0;
    }
  }

  // Decide on the volume rescaling amount, an absolute
  // volume based on the system's initial volume       
  if (tj->ntp == 1) {
    tj->mcdVmax = ((1.0 + tj->MCBarostatFac[0]) * 
		   (1.0 + tj->MCBarostatFac[0]) *
		   (1.0 + tj->MCBarostatFac[0]) - 1.0) *
      (crd->invU.map[0][0] * crd->invU.map[1][1] * crd->invU.map[2][2]);
  }
  if (tj->ntp == 2) {
    tj->mcdVmax = ((1.0 + tj->MCBarostatFac[0]) * 
		   (1.0 + tj->MCBarostatFac[1]) *
		   (1.0 + tj->MCBarostatFac[2]) - 1.0) *
      (crd->invU.map[0][0] * crd->invU.map[1][1] * crd->invU.map[2][2]);
  }
}

//-----------------------------------------------------------------------------
// CurrentSystemPressure: compute the current pressure of the system.   
//
// Arguments:
//   sysUV:  the system energy and virial data structure
//   crd:    system coordinates
//-----------------------------------------------------------------------------
double CurrentSystemPressure(Energy *sysUV, coord *crd)
{
  double virial, volume;

  virial = sysUV->Vir[0] + sysUV->Vir[4] + sysUV->Vir[8];
  volume = crd->invU.data[0]*crd->invU.data[4]*crd->invU.data[8];

  return (2.0*sysUV->kine - virial) / (3.0 * volume);
}

//-----------------------------------------------------------------------------
// BerendsenBarostat: use a Berendsen piston to adjust the system volume
//                    in response to the instantaneous virial.          
//
// Arguments:
//   crd:     system coordinates
//   tp:      system topology
//   tj:      trajectory control data
//   sysUV:   system energy and virial information
//   dcinp:   direct space input control data
//   rcinp:   reciprocal space input control data
//   PPk:     reciprocal space pair potential toolkit
//   CG:      cell grid spanning the system
//-----------------------------------------------------------------------------
void BerendsenBarostat(coord *crd, prmtop *tp, trajcon *tj, Energy *sysUV,
                       dircon *dcinp, reccon *rcinp, bckit *PPk, cellgrid *CG)
{
  double chi, bptc, beta, pressure;
  double chivec[3];

  // Compute chi, the volume rescale
  // factor for isotropic rescaling 
  beta = tj->BerendsenPCoupl;
  bptc = tj->BerendsenPTime;
  sysUV->kine = KineticEnergy(CG, crd, tp, tj);
  pressure = CurrentSystemPressure(sysUV, crd);
  chi = 1.0 - beta * (tj->dt / bptc) * (tj->Ptarget - pressure*PCONVFAC);
  if (chi < 0.95) {
    printf("BerendsenBarostat >> Warning.  Large rescaling value clamped at "
	   "95%% of\nBerendsenBarostat >> instantaneous box size.\n");
    chi = 0.95;
  }
  else if (chi > 1.05) {

    // Terminate the run if this is truly a problem;
    // prevent mdgx from gobbling up system memory! 
    if (chi > 1.10) {
      printf("BerendsenBarostat >> Rescaling value %.4lf.  This is far too "
	     "large and\nBerendsenBarostat >> indicates a problem with the "
	     "system.\n", chi);

      // Run some quick diagnostics to try
      // and indicate what went wrong     
      int i, maxatm, maxatm2;
      double maxfrc, maxfrc2, ifrc, sfx, sfy, sfz;
      maxfrc = 0.0;
      maxfrc2 = 0.0;
      maxatm = -1;
      maxatm2 = -1;
      sfx = 0.0;
      sfy = 0.0;
      sfz = 0.0;
      for (i = 0; i < crd->natom; i++) {
	ifrc = sqrt(crd->frc[3*i]*crd->frc[3*i] +
		    crd->frc[3*i+1]*crd->frc[3*i+1] +
		    crd->frc[3*i+2]*crd->frc[3*i+2]);
	if (ifrc > maxfrc) {
	  maxfrc2 = maxfrc;
	  maxatm2 = maxatm;
	  maxfrc = ifrc;
	  maxatm = i;
	}
	if (ifrc > maxfrc2 && ifrc < maxfrc) {
	  maxfrc2 = ifrc;
	  maxatm2 = i;
	}
	sfx += crd->frc[3*i];
	sfy += crd->frc[3*i+1];
	sfz += crd->frc[3*i+2];
      }
      printf("BerendsenBarostat >> High absolute force: %.4lf on atom %d.\n",
	     maxfrc, maxatm); 
      printf("BerendsenBarostat >> High absolute force: %.4lf on atom %d.\n",
	     maxfrc2, maxatm2);
      printf("BerendsenBarostat >> Total force on system:\n"
	     "Berendsenbarostat >> [ %16.8lf %16.8lf %16.8lf ];\n",
	     sfx, sfy, sfz);
      printf("All forces:\n");
      for (i = 0; i < crd->natom; i++) {
	printf("%16.8lf%16.8lf%16.8lf\n", crd->frc[3*i], crd->frc[3*i+1],
	       crd->frc[3*i+2]);
      }
      exit(1);
    }

    // If the rescaling is not too severe, continue
    printf("BerendsenBarostat >> Warning.  Large rescaling value clamped at "
           "105%% of\nBerendsenBarostat >> instantaneous box size.\n");
    chi = 1.05;
  }

  // Rescale positions by chi ^ 1/3
  chi = pow(chi, 1.0 / 3.0);
  chivec[0] = chi;
  chivec[1] = chi;
  chivec[2] = chi;

  // Update positions
  GridScaleCoord(crd, tp, chivec, dcinp, rcinp, PPk, CG, tj);
}

//-----------------------------------------------------------------------------
// CalcPdV: compute (P dV) for a volume change of the unit cell based on the
//          initial coordinates and a set of rescaling factors.
//
// Arguments:                                                           
//   crd:    the system coordinates
//   chi:    the rescaling factors
//   tj:     the trajectory control data
//-----------------------------------------------------------------------------
static double CalcPdV(coord *crd, trajcon *tj, double *chi)
{
  int i;
  double Eorig, Enew, dE;
  double celldim[6];

  // Save the original unit cell parameters
  for (i = 0; i < 6; i++) {
    celldim[i] = crd->gdim[i];
  }

  // Compute the original pressure-volume energy
  Eorig = (crd->invU.map[0][0] * crd->invU.map[1][1] * crd->invU.map[2][2]) *
    tj->Ptarget;

  // Compute the new dimensions
  for (i = 0; i < 3; i++) {
    crd->gdim[i] *= chi[i];
  }
  CompXfrm(crd->gdim, crd->U, crd->invU);

  // Compute the new pressure-volume energy
  Enew = (crd->invU.map[0][0] * crd->invU.map[1][1] * crd->invU.map[2][2]) *
    tj->Ptarget;

  // Compute the difference in energies and restate it in kcal/mol
  dE = (Enew - Eorig) * 6.0221367e-2 / 4184.0;

  // Reinstate the original cell dimensions to leave them unchanged
  for (i = 0; i < 6; i++) {
    crd->gdim[i] = celldim[i];
  }
  CompXfrm(crd->gdim, crd->U, crd->invU);

  return dE;
}

//-----------------------------------------------------------------------------
// MonteCarloBarostat: this routine implements a Monte-Carlo barostat; while
//                     conceptually simple, requiring merely a coordinate
//                     rescaling followed by an energy recalculation and
//                     acceptance of the move based on the Metropolis
//                     criterion, the barostat may not be very efficient
//                     unless the new energy can be easily computed.  Energy
//                     recalculations are easy if there are no bond
//                     constraints and coordinate rescaling is isotropic.
//
// Arguments:                                                           
//   crd:    the system coordinates                                     
//   tp:     the system topology                                        
//   tj:     the trajectory control data                                
//-----------------------------------------------------------------------------
void MonteCarloBarostat(coord *crd, prmtop *tp, trajcon *tj, Energy *sysUV,
                        dircon *dcinp, FrcTab *Etab, FrcTab *EHtab,
                        reccon *rcinp, bckit *PPk, cellgrid *CG,
                        execon *etimers)
{
  int i;
  double cEtot, nEtot, cEtot2, nEtot2;
  double PVwork, beta, expfac;
  double chi[3], revchi[3], dVmove[4];
  Energy trialUV;

  // Decide on the move to make and compute acceptance probability
  if (tj->ntp == 1) {
    dVmove[0] = tj->mcdVmax*(ran2(&tj->rndcon) - 0.5) /
      (crd->invU.map[0][0]*crd->invU.map[1][1]*crd->invU.map[2][2]);
    dVmove[1] = 0.0;
    dVmove[2] = 0.0;
  }
  else {
    beta = 1.0 / (crd->invU.map[0][0]*crd->invU.map[1][1]*crd->invU.map[2][2]);
    dVmove[0] = tj->mcdVmax*(ran2(&tj->rndcon) - 0.5) * beta;
    dVmove[1] = tj->mcdVmax*(ran2(&tj->rndcon) - 0.5) * beta;
    dVmove[2] = tj->mcdVmax*(ran2(&tj->rndcon) - 0.5) * beta;
  }
  dVmove[3] = ran2(&tj->rndcon);
#ifdef MPI
  MPI_Bcast(&dVmove, 4, MPI_DOUBLE, 0, CG->dspcomm);
#endif
  if (tj->ntp == 1) {
    chi[0] = pow(1.0+dVmove[0], 1.0/3.0);
    chi[1] = chi[0];
    chi[2] = chi[0];
  }
  else {
    chi[0] = pow(1.0+dVmove[0], 1.0/3.0);
    chi[1] = pow(1.0+dVmove[1], 1.0/3.0);
    chi[2] = pow(1.0+dVmove[2], 1.0/3.0);
  }

  // Compute the pressure-volume work
  PVwork = CalcPdV(crd, tj, chi);

  // Compute current system energy
  if (sysUV->Esummed == 0) {
#ifdef MPI
    SumTotalEnergy(CG, sysUV);
#else
    SumTotalEnergy(sysUV);
#endif
  }
  if (tj->TI == 1 && sysUV[1].Esummed == 0) {
#ifdef MPI
    SumTotalEnergy(&CG[1], &sysUV[1]);
#else
    SumTotalEnergy(&sysUV[1]);
#endif
  }
  cEtot = sysUV->eptot;
  if (tj->TI == 1) {
    cEtot2 = sysUV[1].eptot;
  }

  // Buffer the current forces in the coord struct
  MapCellForcesToAtoms(CG, crd);
  if (tj->TI == 1) {
    MapCellForcesToAtoms(&CG[1], &crd[1]);
  }

  // Rescale the coordinates
  GridScaleCoord(crd, tp, chi, dcinp, rcinp, PPk, CG, tj);
  if (tj->TI == 1) {
    GridScaleCoord(&crd[1], &tp[1], chi, dcinp, rcinp, &PPk[1], &CG[1], tj);
  }

  // Recompute the energy, specifically setting the updateU  
  // field of the sysUV struct to 1 to signal that energy and
  // forces must be calculated for the new configuration.    
  InitializeEnergy(&trialUV, tj, tp, 1);
  trialUV.updateU = 1;
  AtomForces(crd, CG, tp, dcinp, Etab, EHtab, rcinp, PPk, &trialUV, etimers,
	     tj);
#ifdef MPI
  SumTotalEnergy(CG, &trialUV);
#else
  SumTotalEnergy(&trialUV);
#endif
  nEtot = trialUV.eptot;
  DestroyEnergyTracker(&trialUV);
  if (tj->TI == 1) {
    InitializeEnergy(&trialUV, tj, &tp[1], 1);
    trialUV.updateU = 1;
    AtomForces(&crd[1], &CG[1], &tp[1], dcinp, Etab, EHtab, rcinp, &PPk[1],
	       &trialUV, etimers, tj);
#ifdef MPI
    SumTotalEnergy(&CG[1], &trialUV);
#else
    SumTotalEnergy(&trialUV);
#endif
    nEtot2 = trialUV.eptot;
    MixCellGrids(&CG[0], &CG[1], tj);
    DestroyEnergyTracker(&trialUV);
  }

  // Now, compare the new energy (nEtot + PVwork) with the
  // old energy (cEtot) and decide whether to accept or   
  // reject the move based on the Metropolis criterion.   
  // The temperature used in the Metropolis criterion is  
  // that used in the system thermostat.                  
  beta = 1.0/(GASCNST*tj->Ttarget);
  if (tj->TI == 1) {
    expfac = (tj->mxA*nEtot + tj->mxB*nEtot2) -
      (tj->mxA*cEtot + tj->mxB*cEtot2) + PVwork -
      tp->nprtcl*log(chi[0]*chi[1]*chi[2])/beta;
  }
  else {
    expfac = nEtot - cEtot + PVwork -
      tp->nprtcl*log(chi[0]*chi[1]*chi[2])/beta;
  }

  // If exp(-beta*expfac) > 1.0, the acceptance is   
  // probability is within (0,1) so no MIN is needed.
  if (dVmove[3] > exp(-beta*expfac)) {

    // Reject the move; return coordinates to their original
    // positions and forces to their original values.       
    for (i = 0; i < 3; i++) {
      revchi[i] = 1.0/chi[i];
    }
    GridScaleCoord(crd, tp, revchi, dcinp, rcinp, PPk, CG, tj);
    MapListForcesToCells(CG, crd);
    if (tj->TI == 1) {
      GridScaleCoord(&crd[1], &tp[1], revchi, dcinp, rcinp, &PPk[1], &CG[1],
		     tj);
      MapListForcesToCells(&CG[1], &crd[1]);
    }
  }
}
