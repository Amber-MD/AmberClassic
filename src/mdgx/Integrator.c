#include <math.h>
#include <stdlib.h>
#include <time.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "Random.h"
#include "Nonbonded.h"
#include "Bonded.h"
#include "pmeRecip.h"
#include "mleRecip.h"
#include "CellManip.h"
#include "BSpline.h"
#include "mdgxVector.h"
#include "Grid.h"
#include "ChargeMap.h"
#include "Constraints.h"
#include "Integrator.h"
#include "Trajectory.h"
#include "Topology.h"
#include "Timings.h"
#include "VirtualSites.h"
#include "Thermostats.h"
#include "Barostats.h"
#include "CrdManip.h"
#include "Debug.h"
#include "Macros.h"
#include "MPIMap.h"

#include "CompFrcDS.h"
#include "pmeDirectDS.h"

//-----------------------------------------------------------------------------
// InitHistory: initialize arrays for coordinate and force history.  This
//              routine must be called after initial locations of atoms are
//              known and before position constraints are applied.
//                                                                      
// Arguments:                                                           
//   crd:     the coordinates struct (every processor will receive a complete
//            set of original coordinates as read from the input coordinates or
//            restart file, and initialize its history to the same values)
//-----------------------------------------------------------------------------
void InitHistory(coord *crd)
{
  int i;

  double *currl, *prvl, *prvf;

  currl = crd->loc;
  prvl = crd->prvloc;
  prvf = crd->prvfrc;
  for (i = 0; i < 3*crd->natom; i++) {
    prvl[i] = currl[i];
    prvf[i] = 0.0;
  }
}

//-----------------------------------------------------------------------------
// AtomForces: the complete force calculation routine, with branches for
//             special cases in which the energy or virial is desired.  
//                                                                      
// Arguments:                                                           
//   crd:      the coordinates (for velocity information)               
//   CG:       the cell grid (this is where coordinates and forces are really
//             being stored)                                     
//   tp:       the topology (for masses)                                
//   dcinp:    the direct space parameters                              
//   Etab:     the direct space electrostatic interaction spline table  
//   EHtab:    the high-resolution direct space electrostatic interaction
//             spline table                                 
//   rcinp:    the reciprocal space control parameters                  
//   PPk:      the reciprocal space pair potential mesh construction kit
//   sysUV:    state information data (potential energy and virial)     
//   etimers:  timing information                                       
//   tj:       trajectory control information (for MPI process maps)    
//-----------------------------------------------------------------------------
void AtomForces(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
		FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
		Energy *sysUV, execon *etimers, trajcon *tj)
{
  int i;
  cell *C;
#ifdef MPI
  int j, nreq;
  MPI_Request* req;
  MPI_Status* stt;
  
  // Allocate memory for MPI requests
  nreq = CG->nthreads;
  req = (MPI_Request*)malloc(nreq*sizeof(MPI_Request));
  stt = (MPI_Status*)malloc(nreq*sizeof(MPI_Status));
#endif

  // Initialize forces on all cells' primary  
  // sectors and share positions of particles;
  // place extra points                       
  ZeroCellForces(CG);
  etimers->nbInt += mdgxStopTimer(etimers);
#ifdef MPI
  ShareCoordinates(CG, tj, dcinp->Mcut, crd, tp, 0);
  ExtraPointLocations(tp, crd, CG, tj);
#else
  ShareCoordinates(CG, dcinp->Mcut, crd, tp, 0);
  ExtraPointLocations(tp, crd, CG);
#endif

  // Update the cells' atom GPS
  UpdateCellGPS(CG);
  
  // Much of this work counts as cell communication
  etimers->cellcomm += mdgxStopTimer(etimers);

  // Spline computations
  for (i = 0; i < CG->MyCellCount; i++) {
    CellSplCoeff(&CG->data[CG->MyCellDomain[i]], crd, rcinp);
  }
  etimers->nbBsp += mdgxStopTimer(etimers);

  // Particle -> mesh mapping
  i = (rcinp->QL[0].pfft == 1) ? 2*(rcinp->QL[0].col/2+1) : rcinp->QL[0].col;
  SetDVec(rcinp->QL[0].data, rcinp->QL[0].pag*rcinp->QL[0].row*i, 0.0);
  for (i = 0; i < CG->MyCellCount; i++) {
    CellQ2Book(&CG->data[CG->MyCellDomain[i]], rcinp, &rcinp->QL[0]);
  }
  etimers->nbPtM += mdgxStopTimer(etimers);

#ifdef MPI
  // Communicate the results of the splining
  // to processes that will perform FFTs    
  if (rcinp->nlev == 1) {
    nreq = InitMeshGatherSPME(rcinp, CG, req, 0);
    etimers->mpiMeshPack += mdgxStopTimer(etimers);
  }
  else {
    nreq = 0;
    if (CG->nthreads > 1) {
      printf("AtomForces >> Error!  Not yet ready for parallel MLE.\n");
      exit(1);
    }
  }
#endif
  
  // We must enter at least one iteration of
  // this loop, even if CG->MyCellCount is 0
  for (i = 0; i <= CG->MyCellCount; i++) {

    // Asynchronous communication for mesh merger 
    // should be finalized by halfway through this
    // processor's direct space calculations; the 
    // second half of the direct space calculation
    // will be used to mask asynchronous mesh     
    // redistribution.                            
    if (i == CG->MasterHalfLoad) {
      etimers->nbInt += mdgxStopTimer(etimers);
#ifdef MPI
      if (nreq > 0) {
        MPI_Waitall(nreq, req, stt);
      }
      MPI_Barrier(CG->dspcomm);
      etimers->mpiMeshPullWait += mdgxStopTimer(etimers);
#endif
      if (CG->tid == 0) {
	      if (rcinp->nlev == 1) {
#ifdef MPI
          // Unpack the mesh components from other processes
          for (j = 1; j < CG->nthreads; j++) {
            RecvMeshPart(rcinp, CG, j);
          }
          etimers->mpiMeshPack += mdgxStopTimer(etimers);
#endif
          // The usual Smooth Particle Mesh Ewald calculation
          if (sysUV->updateU == 2) {
            ConvQBCnrgvir(rcinp, crd, &rcinp->QL[0], PPk, sysUV, etimers);
          }
          else if (sysUV->updateU == 1 || sysUV->updateU == -1) {
            ConvQBCnrg(rcinp, crd, &rcinp->QL[0], PPk, sysUV, etimers);
          }
          else {
            ConvQBC(rcinp, crd, &rcinp->QL[0], PPk, etimers);
          }
        }
        else {

	  // Multi-Level Ewald calculation
          mleConvQBC(rcinp, sysUV, PPk, etimers);
        }
      }
#ifdef MPI
      InitMeshGatherSPME(rcinp, CG, req, 1);
      etimers->mpiMeshPack += mdgxStopTimer(etimers);
#endif
    }

    // Continue if we have already exceeded the cell count;
    // this is for the case where one process (i.e. the    
    // master process) actually controls zero cells.       
    if (i == CG->MyCellCount) {
      continue;
    }
    C = &CG->data[CG->MyCellDomain[i]];

    // Direct space force calculation ****/
    DirectTriple2(C, crd, Etab, dcinp, tp, sysUV);
    etimers->nbInt += mdgxStopTimer(etimers);

    // Bonded interactions
    C = &CG->data[CG->MyCellDomain[i]];
    CellBondedIntr(C, CG, crd, tp, Etab, EHtab, sysUV, rcinp->qform);
    etimers->bonds += mdgxStopTimer(etimers);
  }  
  
  // Apply long-ranged vdW correction
  if (dcinp->LRvdw == 1) {
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      CellLRvdw(C, crd, tp, sysUV);
    }
  }
  etimers->nbInt += mdgxStopTimer(etimers);

  // Interpolate electrostatics, mesh -> particle
  if (sysUV->updateU >= 0) {
#ifdef MPI
    if (nreq > 0) {
      MPI_Waitall(nreq, req, stt);
      etimers->mpiMeshPushWait += mdgxStopTimer(etimers);
      if (CG->tid != 0) {
	      RecvMeshPart(rcinp, CG, 0);
	      etimers->mpiMeshPack += mdgxStopTimer(etimers);
      }
    }
    MPI_Barrier(CG->dspcomm);
#endif
    for (i = 0; i < CG->MyCellCount; i++) {
      CellIntrpFrc(&CG->data[CG->MyCellDomain[i]], crd, rcinp, &rcinp->QL[0]);
    }
    etimers->nbMtP += mdgxStopTimer(etimers);
  }

  // Merge atom forces and map back to unified list
#ifdef MPI
  MergeCellForces(CG, crd, tp, tj);
#else
  MergeCellForces(CG, crd, tp);
#endif
  etimers->cellcomm += mdgxStopTimer(etimers);

#ifdef MPI
  // Free allocated memory
  free(req);
  free(stt);
#endif
}

//-----------------------------------------------------------------------------
// CellVerletC: velocity Verlet algorithm for updating positions.  This puts
//              the current coordinates of the atoms in the prvloc array of the
//              coord struct and the current set of forces in the prvfrc array
//              of the coord struct.  Works in a cell-to-cell framework,
//              releasing atoms into new cells as appropriate.
//
// Arguments:                                                           
//   CG:   the cell grid                                                
//   tp:   the topology                                                 
//   tj:   trajectory control information, including time step and flag 
//         for Langevin thermostat.  When defined with MPI, contains MPI
//         mapping information like process rank.                       
//   crd:  the coordinates (also contains velocity information)         
//-----------------------------------------------------------------------------
static void CellVerletC(cellgrid *CG, prmtop *tp, trajcon *tj, coord *crd)
{
  int i, j, gcon, g3con;
  int* moving;
  double pmovex, pmovey, pmovez, hmdt;
  double *ctmp, *ltmp, *vtmp, *ftmp, *mtmp, *pftmp;
  cell *C;

  // Determine if a belly mask is in effect
  const int usebelly = (tj->Leash.active == 1) ? tj->Leash.usebelly : 0;
  moving = tp->MobileAtoms;

  // Loop over all atoms in all cells, update positions,
  // and save the current set of coordinates and forces.
  mtmp = tp->InvMasses;
  if (tj->ntt == 3) {
    pftmp = crd->prvfrc;
  }
  const double dt = sqrt(418.4)*tj->dt;
  const double hdt = 0.5*dt;
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];

    // For a Langevin integrator, add perturbations to 
    // the forces here.  These perturbations will enter
    // into the new positions and get carried forward  
    // to the calculation of new velocities.           
    if (tj->ntt == 3) {
      long counter = tj->rndcon;
      const double c_explic = tj->lnth.c_explic;
      for (j = 0; j < C->nr[0]; j++) {
	gcon = C->data[j].id;
	hmdt = hdt*mtmp[gcon];
	ctmp = C->data[j].loc;
	ftmp = C->data[j].frc;
	g3con = 3*gcon;
	vtmp = &crd->vel[g3con];
	ltmp = &crd->loc[g3con];
	crd->prvloc[g3con] = ltmp[0];
	crd->prvloc[g3con+1] = ltmp[1];
	crd->prvloc[g3con+2] = ltmp[2];
	vtmp[0] = vtmp[0]*c_explic + hmdt*(ftmp[0] + pftmp[g3con]);
	vtmp[1] = vtmp[1]*c_explic + hmdt*(ftmp[1] + pftmp[g3con+1]);
	vtmp[2] = vtmp[2]*c_explic + hmdt*(ftmp[2] + pftmp[g3con+2]);
	pmovex = dt*vtmp[0];
	pmovey = dt*vtmp[1];
	pmovez = dt*vtmp[2];
	if (usebelly == 0 || moving[gcon] == 1) {
	  ctmp[0] += pmovex;
	  ltmp[0] += pmovex;
	  ctmp[1] += pmovey;
	  ltmp[1] += pmovey;
	  ctmp[2] += pmovez;
	  ltmp[2] += pmovez;
	}
      }
      tj->rndcon = counter;
    }
    else {
      for (j = 0; j < C->nr[0]; j++) {
	gcon = C->data[j].id;
	hmdt = hdt*mtmp[gcon];
	ctmp = C->data[j].loc;
	ftmp = C->data[j].frc;
	g3con = 3*gcon;
	vtmp = &crd->vel[g3con];
	ltmp = &crd->loc[g3con];
	crd->prvloc[g3con] = ltmp[0];
	crd->prvloc[g3con+1] = ltmp[1];
	crd->prvloc[g3con+2] = ltmp[2];
	vtmp[0] += hmdt*ftmp[0];
	vtmp[1] += hmdt*ftmp[1];
	vtmp[2] += hmdt*ftmp[2];
	pmovex = dt*vtmp[0];
	pmovey = dt*vtmp[1];
	pmovez = dt*vtmp[2];
	if (usebelly == 0 || moving[gcon] == 1) {
	  ctmp[0] += pmovex;
	  ltmp[0] += pmovex;
	  ctmp[1] += pmovey;
	  ltmp[1] += pmovey;
	  ctmp[2] += pmovez;
	  ltmp[2] += pmovez;
	}
      }
    }
  }

  // Now, atoms may have migrated out of cells.
  // Place them in new cells, as appropriate.  
#ifdef MPI
  UpdateCells(CG, crd, tp, tj);
#else
  UpdateCells(CG, crd, tp);
#endif
}

//-----------------------------------------------------------------------------
// CellVerletV: the velocity Verlet algorithm for updating velocities, in a
//              cell-based implementation.                         
//                                                                      
// Arguments:                                                           
//   CG:        the cell grid (coordinate information)                  
//   tp:        the topology (for masses)                               
//   tj:        trajectory control data (for timestep and Langevin dynamics
//              flag)                                          
//   crd:       the coordinates (for velocities)                        
//-----------------------------------------------------------------------------
static void CellVerletV(cellgrid *CG, prmtop *tp, trajcon *tj, coord *crd)
{
  int i, j, g3con;
  double dt;
  double *vtmp, *mtmp, *masstmp, *pftmp;
  cell *C;
  atomc *cdtmp;

  // Pointers to coordinate or topology structs
  dt = sqrt(418.4)*tj->dt;
  vtmp = crd->vel;
  mtmp = tp->InvMasses;

  // Simple model for Langevin dynamics, basically taken from     
  // Loncharich, Brooks and Pastor, Biopolymers 32:523-535 (1992),
  // Eq. 11.  (Note that the first term on the rhs of Eq. 11b     
  // should not be there.)                                        
  // For more info: see Pastor, Brooks & Szabo, Mol. Phys. 65:1409
  // 1988.  More detailed comments are in sff.c                   
  if (tj->ntt == 3) {
    pftmp = crd->prvfrc;
    masstmp = tp->Masses;
    long int counter = tj->rndcon;
    const double c_implic = tj->lnth.c_implic;
    const double sdfac = tj->lnth.sdfac;
    const double hdt = 0.5*dt;
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      for (j = 0; j < C->nr[0]; j++) {
	cdtmp = &C->data[j];
	const double rsd = sdfac * sqrt(masstmp[cdtmp->id]);
	const double hmdt = hdt*mtmp[cdtmp->id];	
	g3con = 3*cdtmp->id;
	pftmp[g3con] = rsd*GaussBoxMuller(&counter);
	pftmp[g3con+1] = rsd*GaussBoxMuller(&counter);
	pftmp[g3con+2] = rsd*GaussBoxMuller(&counter);
        vtmp[g3con] = (hmdt*(cdtmp->frc[0] + pftmp[g3con]) + 
		       vtmp[g3con]) * c_implic;
	vtmp[g3con+1] = (hmdt*(cdtmp->frc[1] + pftmp[g3con+1]) +
			 vtmp[g3con+1]) * c_implic;
	vtmp[g3con+2] = (hmdt*(cdtmp->frc[2] + pftmp[g3con+2]) +
			 vtmp[g3con+2]) * c_implic;
      }
    }
    tj->rndcon = counter;
  }
  else {
    const double hdt = 0.5*dt;
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      for (j = 0; j < C->nr[0]; j++) {
	cdtmp = &C->data[j];
	const double hmdt = hdt*mtmp[cdtmp->id];
	g3con = 3*cdtmp->id;
	vtmp[g3con] += hmdt*cdtmp->frc[0];
	vtmp[g3con+1] += hmdt*cdtmp->frc[1];
	vtmp[g3con+2] += hmdt*cdtmp->frc[2];
      }
    }
  }
}

#ifdef MPI
//-----------------------------------------------------------------------------
// MergeEnergy: when MPI is enabled, this function synchronizes energy  
//              accumulators across all cells.                          
//-----------------------------------------------------------------------------
static void MergeEnergy(cellgrid *CG, Energy *sysUV)
{
  int i, nelem;
  double* Eacc;
  double* Racc;

  // Copy energy components into a scratch array
  nelem = 7+sysUV->nUdc;
  Eacc = (double*)malloc(nelem*sizeof(double));
  Racc = (double*)calloc(nelem, sizeof(double));
  Eacc[0] = sysUV->bond;
  Eacc[1] = sysUV->angl;
  Eacc[2] = sysUV->dihe;
  Eacc[3] = sysUV->delec;
  Eacc[4] = sysUV->relec;
  Eacc[5] = sysUV->vdw12;
  Eacc[6] = sysUV->vdw6;
  for (i = 0; i < sysUV->nUdc; i++) {
    Eacc[i+7] = sysUV->BondUdc[i];
  }

  // Allreduce!  Yikes!  Necessary evil!
  MPI_Allreduce(Eacc, Racc, nelem, MPI_DOUBLE, MPI_SUM, CG->dspcomm);

  // Write back the reduced energies
  sysUV->bond = Racc[0];
  sysUV->angl = Racc[1];
  sysUV->dihe = Racc[2];
  sysUV->delec = Racc[3];
  sysUV->relec = Racc[4];
  sysUV->vdw12 = Racc[5];
  sysUV->vdw6 = Racc[6];
  for (i = 0; i < sysUV->nUdc; i++) {
    sysUV->BondUdc[i] = Racc[i+6];
  }

  // Free allocated memory
  free(Eacc);
  free(Racc);
}
#endif

//-----------------------------------------------------------------------------
// SumTotalEnergy: this routine computes the total potential energy as well as
//                 the total (kinetic plus potential) energy.  The kinetic
//                 energy is not updated inside this function, however, and so
//                 must be computed prior to calling this function.
//
// Arguments:                                                           
//   sysUV:  the energy / virial accumulator                            
//-----------------------------------------------------------------------------
#ifdef MPI
void SumTotalEnergy(cellgrid *CG, Energy *sysUV)
#else
void SumTotalEnergy(Energy *sysUV)
#endif
{
#ifdef MPI
  MergeEnergy(CG, sysUV);
#endif
  sysUV->eptot = sysUV->bond + sysUV->angl + sysUV->dihe + sysUV->cmap +
                 sysUV->delec + sysUV->relec + sysUV->vdw12 + sysUV->vdw6;
  sysUV->etot = sysUV->eptot + sysUV->kine;
  sysUV->elec = sysUV->delec + sysUV->relec;
  sysUV->Esummed = 1;
}

//-----------------------------------------------------------------------------
// WriteTrajFiles: write all the requested trajectory files if this step is one
//                 of interest.  The restart file is not written at this time,
//                 but rather earlier in the file, because the checkpointing
//                 information must be sufficient to continue the run as if it
//                 were continuous.  Restart files are therefore written before
//                 this routine is called, after the velocity update but before
//                 velocity corrections for any constraints.  Checkpoint data
//                 contains position, but not velocity corrections, and does
//                 not contain the latest thermostat or barostat corrections as
//                 are found in the trajectory data end the same step number.
//                 The TI variant will write twin coordinate files but only one
//                 unified output file.   
//
// Arguments:                                                           
//   tj:       trajectory control parameters                            
//   tp:       the topology (in the case of TI, a pair of topologies)   
//   dcinp:    direct space control data                                
//   rcinp:    reciprocal space control data                            
//   Etab:     the electrostatic direct space spline interpolation table
//   crd:      the coordinates                                          
//   CG:       the cell grid                                            
//   sysUV:    system energy and virial information                     
//   etimers:  timing data                                              
//   Acdf:     NetCDF control structs (coordinates, velocities, and / or
//             forces)                                                  
//   n:        data for the nth system is to be printed                 
//-----------------------------------------------------------------------------
static void WriteTrajFiles(trajcon *tj, prmtop *tp, dircon *dcinp,
                           reccon *rcinp, FrcTab *Etab, cellgrid *CG,
                           coord *crd, Energy *sysUV, execon *etimers,
                           cdftrj* Acdf, int n)
{
  // Output diagnostics
  if (tj->TI == 1 && n == 0) {
    if ((tj->ntpr != 0 && tj->currstep % tj->ntpr == 0) ||
	tj->currstep == tj->nstep ||
	(tj->nfistep > 0 && tj->currstep % tj->nfistep == 0)) {
      WriteDiagnostics(tj, tp, dcinp, rcinp, Etab, CG, crd, sysUV, etimers, n);
    }
  }
  else if (tj->TI != 1) {
    if ((tj->ntpr != 0 && tj->currstep % tj->ntpr == 0) ||
	tj->currstep == tj->nstep ||
        (tj->nfistep > 0 && tj->currstep % tj->nfistep == 0)) {
      WriteDiagnostics(tj, tp, dcinp, rcinp, Etab, CG, crd, sysUV, etimers, n);
    }
  }

  // Coordinates trajectory
  if (tj->ntwx != 0 && tj->currstep % tj->ntwx == 0) {
    if (tj->ioutfm == 0) {
      WriteCrd(CG, crd, 1, tj, tp, n);
    }
    else if (tj->ioutfm == 1) {
      WriteCDF(CG, crd, 1, tj, Acdf, tp, n);
    }
  }

  // Velocities trajectory
  if (tj->ntwv != 0 && tj->currstep % tj->ntwv == 0) {
    if (tj->ioutfm == 0) {
      WriteCrd(CG, crd, 2, tj, tp, n);
    }
    else if (tj->ioutfm == 1) {
      WriteCDF(CG, crd, 2, tj, Acdf, tp, n);
    }
  }

  // Forces trajectory
  if (tj->ntwf != 0 && tj->currstep % tj->ntwf == 0) {
    if (tj->ioutfm == 0) {
      WriteCrd(CG, crd, 3, tj, tp, n);
    }
    else if (tj->ioutfm == 1) {
      WriteCDF(CG, crd, 3, tj, Acdf, tp, n);
    }
  }
}

//-----------------------------------------------------------------------------
// ComputeDVDL: given two systems, their energies, and the value of the mixing
//              parameter, compute DV/DL.                        
//                                                                      
// Arguments:                                                           
//   sys[A,B]:   systems A and B energy state information               
//   tj:         trajectory control information                         
//-----------------------------------------------------------------------------
static double ComputeDVDL(Energy *sysA, Energy *sysB, trajcon *tj)
{
  return tj->dmxA*sysA->eptot + tj->dmxB*sysB->eptot;
}

#define NEED_TI 0
#include "IntegratorBranch.c"
#undef NEED_TI

#define NEED_TI 1
#include "IntegratorBranch.c"
#undef NEED_TI

//-----------------------------------------------------------------------------
// InitVelocities: initialize the velocities of all atoms in the system (or
//                 systems in the case of TI) to random values in a Maxwell
//                 distribution, if velocities are not yet assigned, or
//                 resumes the dynamics if velocities have been read from a
//                 restart file.                       
//
// Arguments:                                                           
//   crd:      the coordinates (already initialized by InitCoords)      
//   tp:       the topology (for masses)                                
//   T:        the temperature                                          
//   dcinp:    the direct space control parameters                      
//   Etab:     the direct space electrostatic interaction spline table  
//   EHtab:    the high-resolution direct space electrostatic           
//             interaction spline table                                 
//   rcinp:    the reciprocal space control parameters                  
//   Q:        the reciprocal space charge / potential mesh             
//   PPk:      the reciprocal space pair potential mesh construction kit
//   tj:       trajectory control information                           
//   sysUV:    energy and virial information                            
//   etimers:  timings data                                             
//   Acdf:     NetCDF control structs                                   
//   n:        the nth trajectory is being initialized                  
//-----------------------------------------------------------------------------
void InitVelocities(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
		    FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
		    trajcon *tj, Energy *sysUV, execon *etimers, cdftrj* Acdf,
		    int n)
{
  long long int currstep;

  // Group bonded atoms.  Initialize coordinates      
  // arrays not present in the checkpoint information.
  ImageBondedGroups(crd, tp);
  InitHistory(crd);
  if (tj->TI == 1) {
    ImageBondedGroups(&crd[1], &tp[1]);
    InitHistory(&crd[1]);
  }

  // The first order of business is to take coordinates
  // from the unified list into cell grids.            
  AtomsToCells(crd, CG, tp);
  if (tj->TI == 1) {
    AtomsToCells(&crd[1], &CG[1], &tp[1]);
  }

  // Random initial velocities
  if (tj->irest == 0 && tj->mode == 0) {
    AndersenThermostat(CG, crd, tp, tj, tj->Tinit);
  }
  etimers->Thermostat += mdgxStopTimer(etimers);
    
  // Constraints are applied to positions simply for the     
  // purpose of overcoming imprecision in the checkpoint data
  // if this is a restart from a simulation in progress.     
  ApplyGridCnst(CG, crd, tp, tj, dcinp->Mcut, 0, crd->prvfrc);
  if (tj->TI == 1) {
    ApplyGridCnst(&CG[1], &crd[1], &tp[1], tj, dcinp->Mcut, 0, crd[1].prvfrc);
  }
  etimers->Constraints += mdgxStopTimer(etimers);

  // Explicitly synchronize trajectory coordinates
  if (tj->TI == 1 && tj->nsynch > 0 && tj->currstep % tj->nsynch == 0) {
    SynchronizeCoordinates(CG, tj);
  }

  // Compute forces
  InitializeEnergy(sysUV, tj, tp, 0);
  AtomForces(crd, CG, tp, dcinp, Etab, EHtab, rcinp, PPk, sysUV, etimers, tj);
  if (tj->TI == 1) {
    InitializeEnergy(&sysUV[1], tj, &tp[1], 0);
    AtomForces(&crd[1], &CG[1], &tp[1], dcinp, Etab, EHtab, rcinp, &PPk[1],
	       &sysUV[1], etimers, tj);
    MixCellGrids(&CG[0], &CG[1], tj);
#ifdef MPI
    SumTotalEnergy(&CG[0], &sysUV[0]);
    SumTotalEnergy(&CG[1], &sysUV[1]);
#else
    SumTotalEnergy(&sysUV[0]);
    SumTotalEnergy(&sysUV[1]);
#endif
    sysUV[0].dVdL = ComputeDVDL(&sysUV[0], &sysUV[1], tj);
    sysUV[0].AVEdVdL += sysUV[0].dVdL;
    sysUV[0].RMSdVdL += sysUV[0].dVdL*sysUV[0].dVdL;
  }

  // Velocity Verlet, update velocities
  CellVerletV(CG, tp, tj, crd);
  if (tj->TI == 1) {
    CellVerletV(&CG[1], &tp[1], tj, &crd[1]);
  }
  etimers->Integ += mdgxStopTimer(etimers);

  // Velocity constraints
  currstep = tj->currstep;
  tj->currstep = 0;
  if (tj->ntp > 0 && tj->barostat == 1) {
    ReflectDVec(crd->prvvel, crd->vel, 3*tp->natom);
  }
  ApplyGridCnst(CG, crd, tp, tj, dcinp->Mcut, 1, crd->prvfrc);
  if (tj->ntp > 0 && tj->barostat == 1) {
    VVConstraintVirial(crd, tp, tj, sysUV);
  }
  if (tj->TI == 1) {
    if (tj->ntp > 0 && tj->barostat == 1) {
      ReflectDVec(crd[1].prvvel, crd[1].vel, 3*tp[1].natom);
    }
    ApplyGridCnst(&CG[1], &crd[1], &tp[1], tj, dcinp->Mcut, 1, crd->prvfrc);
    if (tj->ntp > 0 && tj->barostat == 1) {
      VVConstraintVirial(&crd[1], &tp[1], tj, &sysUV[1]);
    }
  }
  tj->currstep = currstep;
  etimers->Constraints += mdgxStopTimer(etimers);

  // Compute kinetic energy and temperature after the  
  // force computation in preparation for thermostating
  if (tj->ntt == 1 || (tj->ntt == 2 && tj->currstep % tj->vrand == 0)) {
    if (tj->TI == 1) {
      sysUV[0].T = SystemTemperatureTI(CG, crd, tp, sysUV, tj, 1);
      sysUV[1].kine = sysUV[0].kine;
      sysUV[1].T = sysUV[0].T;
    }
    else {
      sysUV->T = SystemTemperature(CG, crd, tp, sysUV, tj, 1);
    }
  }
  etimers->Thermostat += mdgxStopTimer(etimers);

  // Write trajectory files and outputs
  if (tj->irest == 0 || tj->currstep == 0) {
    WriteTrajFiles(tj, tp, dcinp, rcinp, Etab, CG, crd, sysUV, etimers,
		   Acdf, n);
    if (tj->TI == 1) {
      WriteTrajFiles(tj, tp, dcinp, rcinp, Etab, CG, crd, sysUV, etimers,
		     Acdf, 1);
    }
    etimers->Write += mdgxStopTimer(etimers);
  }

  // Berendsen thermostat
  if (tj->ntt == 1) {
    if (tj->TI == 1) {
      BerendsenThermostatTI(CG, crd, tj, sysUV);
    }
    else {
      BerendsenThermostat(CG, crd, tj, sysUV);
    }
  }

  // Andersen thermostat
  if (tj->ntt == 2 && tj->currstep % tj->vrand == 0) {
    AndersenThermostat(CG, crd, tp, tj, tj->Ttarget);
  }
  etimers->Thermostat += mdgxStopTimer(etimers);

  // Berendsen barostat
  if (tj->ntp > 0 && tj->barostat == 1) {
    BerendsenBarostat(crd, tp, tj, sysUV, dcinp, rcinp, PPk, CG);
  }

  // Monte-Carlo barostat
  if (tj->ntp > 0 && tj->barostat == 2 &&
      tj->currstep % tj->MCBarostatFreq == 0) {
    MonteCarloBarostat(crd, tp, tj, sysUV, dcinp, Etab, EHtab, rcinp, PPk, CG,
		       etimers);
  }
  etimers->Barostat += mdgxStopTimer(etimers);

  // Remove system momentum
  if (tj->RemoveMomentum > 0 && tj->currstep % tj->RemoveMomentum == 0) {
    if (tj->TI == 1) {
      RemoveMomentumTI(CG, crd, tp, tj);
    }
    else {
      RemoveMomentum(CG, crd, tp);
    }
  }
  etimers->Thermostat += mdgxStopTimer(etimers);
}

//-----------------------------------------------------------------------------
// Dynamics: routine to manage one step of dynamics, whether that be an 
//           ensemble of constant volume, temperature, or both.      
//
// Arguments:                                                           
//   crd:     the coordinates array                                     
//   CG:      the cell grid array                                       
//   tp:      the topology array                                        
//   dcinp:   the direct space control parameters                       
//   Etab:    the direct space electrostatic interaction spline table   
//   EHtab:   the high-resolution direct space electrostatic interaction
//            spline table                                              
//   rcinp:   the reciprocal space control parameters                   
//   Q:       the reciprocal space charge / potential mesh              
//   PPk:     the reciprocal space pair potential mesh construction kits
//   tj:      trajectory control parameters                             
//   sysUV:   energy and virial information                             
//   etimers: timing data                                               
//   Acdf:    NetCDF file/output pointers                               
//-----------------------------------------------------------------------------
void Dynamics(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
	      FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
	      trajcon *tj, Energy *sysUV, execon *etimers, cdftrj* Acdf,
	      int n)
{
  // Velocity Verlet, update positions
  CellVerletC(CG, tp, tj, crd);
  if (tj->TI == 1) {
    CellVerletC(&CG[1], &tp[1], tj, &crd[1]);
  }
  etimers->Integ += mdgxStopTimer(etimers);

  // Position constraints
  if (tj->ntp > 0 && tj->barostat == 1) {
    ReflectDVec(crd->prvvel, crd->vel, 3*tp->natom);
  }
  ApplyGridCnst(CG, crd, tp, tj, dcinp->Mcut, 0, crd->prvfrc);
  if (tj->ntp > 0 && tj->barostat == 1) {
    VVConstraintVirial(crd, tp, tj, sysUV);
  }
  if (tj->TI == 1) {
    if (tj->ntp > 0 && tj->barostat == 1) {
      ReflectDVec(crd[1].prvvel, crd[1].vel, 3*tp[1].natom);
    }
    ApplyGridCnst(&CG[1], &crd[1], &tp[1], tj, dcinp->Mcut, 0, crd->prvfrc);
    if (tj->ntp > 0 && tj->barostat == 1) {
      VVConstraintVirial(&crd[1], &tp[1], tj, &sysUV[1]);
    }
  }
  etimers->Constraints += mdgxStopTimer(etimers);

  // Explicitly synchronize trajectory coordinates
  if (tj->TI == 1 && tj->nsynch > 0 && tj->currstep % tj->nsynch == 0) {
    SynchronizeCoordinates(CG, tj);
  }

  // Compute forces
  InitializeEnergy(sysUV, tj, tp, 0);
  AtomForces(crd, CG, tp, dcinp, Etab, EHtab, rcinp, PPk, sysUV, etimers, tj);
  if (tj->TI == 1) {
    InitializeEnergy(&sysUV[1], tj, &tp[1], 0);
    AtomForces(&crd[1], &CG[1], &tp[1], dcinp, Etab, EHtab, rcinp, &PPk[1],
	       &sysUV[1], etimers, tj);
    MixCellGrids(&CG[0], &CG[1], tj);
    sysUV[0].dVdL = ComputeDVDL(&sysUV[0], &sysUV[1], tj);
    sysUV[0].AVEdVdL += sysUV[0].dVdL;
    sysUV[0].RMSdVdL += sysUV[0].dVdL*sysUV[0].dVdL;
  }

  // Velocity Verlet, update velocities
  CellVerletV(CG, tp, tj, crd);
  if (tj->TI == 1) {
    CellVerletV(&CG[1], &tp[1], tj, &crd[1]);
  }
  etimers->Integ += mdgxStopTimer(etimers);

  // Write restart file(s)
  if (tj->ntwr != 0 && tj->currstep % tj->ntwr == 0 &&
      tj->currstep >= tj->ntwr) {
    WriteRst(CG, crd, tp, tj, n, 1);
    if (tj->TI == 1) {
      WriteRst(&CG[1], &crd[1], &tp[1], tj, 1, 1);
    }
    etimers->Write += mdgxStopTimer(etimers);
  }

  // Velocity constraints
  if (tj->ntp > 0 && tj->barostat == 1) {
    ReflectDVec(crd->prvvel, crd->vel, 3*tp->natom);
  }
  ApplyGridCnst(CG, crd, tp, tj, dcinp->Mcut, 1, crd->prvfrc);
  if (tj->ntp > 0 && tj->barostat == 1) {
    VVConstraintVirial(crd, tp, tj, sysUV);
  }
  if (tj->TI == 1) {
    if (tj->ntp > 0 && tj->barostat == 1) {
      ReflectDVec(crd[1].prvvel, crd[1].vel, 3*tp[1].natom);
    }
    ApplyGridCnst(&CG[1], &crd[1], &tp[1], tj, dcinp->Mcut, 1, crd->prvfrc);
    if (tj->ntp > 0 && tj->barostat == 1) {
      VVConstraintVirial(&crd[1], &tp[1], tj, &sysUV[1]);
    }
  }
  etimers->Constraints += mdgxStopTimer(etimers);

  // Compute kinetic energy after the force computation  
  // in preparation for Andersen or Berendsen thermostats
  if (tj->ntt == 1 || tj->ntt == 2) {
    if (tj->TI == 1) {
      sysUV[0].T = SystemTemperatureTI(CG, crd, tp, sysUV, tj, 1);
      sysUV[1].T = sysUV[0].T;
    }
    else {
      sysUV->T = SystemTemperature(CG, crd, tp, sysUV, tj, 1);
    }
  }
  etimers->Thermostat += mdgxStopTimer(etimers);

  // Write trajectory files and outputs
  WriteTrajFiles(tj, tp, dcinp, rcinp, Etab, CG, crd, sysUV, etimers,
                 Acdf, n);
  if (tj->TI == 1) {
    WriteTrajFiles(tj, &tp[1], dcinp, rcinp, Etab, &CG[1], &crd[1], &sysUV[1],
		   etimers, &Acdf[3], 1);
  }
  etimers->Write += mdgxStopTimer(etimers);

  // Berendsen thermostat
  if (tj->ntt == 1) {
    if (tj->TI == 1) {
      BerendsenThermostatTI(CG, crd, tj, sysUV);
    }
    else {
      BerendsenThermostat(CG, crd, tj, sysUV);
    }
  }

  // Andersen thermostat
  if (tj->ntt == 2 && tj->currstep % tj->vrand == 0) {
    AndersenThermostat(CG, crd, tp, tj, tj->Ttarget);
  }
  etimers->Thermostat += mdgxStopTimer(etimers);

  // Berendsen barostat
  if (tj->ntp > 0 && tj->barostat == 1) {
    BerendsenBarostat(&crd[0], &tp[0], tj, &sysUV[0], dcinp, rcinp, &PPk[0],
		      &CG[0]);
  }

  // Monte-Carlo barostat
  if (tj->ntp > 0 && tj->barostat == 2 &&
      tj->currstep % tj->MCBarostatFreq == 0) {
    MonteCarloBarostat(crd, tp, tj, &sysUV[0], dcinp, Etab, EHtab,
		       rcinp, &PPk[0], &CG[0], etimers);
  }
  etimers->Barostat += mdgxStopTimer(etimers);

  // Remove system momentum
  if (tj->RemoveMomentum > 0 && tj->currstep % tj->RemoveMomentum == 0) {
    if (tj->TI == 1) {
      RemoveMomentumTI(CG, crd, tp, tj);
    }
    else {
      RemoveMomentum(CG, crd, tp);
    }
  }
  etimers->Thermostat += mdgxStopTimer(etimers);
}

#if 0
//-----------------------------------------------------------------------------
// XminAirlock: this function serves for transitions within the xminC() 
//              function of the NAB code (xminC.c in the directory      
//              ${AMBERHOME}/AmberTools/sff/ ), creating buffer arrays  
//              for particle positions and gradients and copying them   
//              into and out of the cell grid.                          
//                                                                      
// Arguments:                                                           
//-----------------------------------------------------------------------------
void XminAirlock(cellgrid *CG, coord *crd, prmtop *tp, dricon *dcinp,
		 FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
		 trajcon *tj, Energy *sysUV, execon *etimers, cdftrj* Acdf)
{
  int i, j, k, g3con;
  double* xpos;
  double* xgrad;
  cell *C;

  // Compute the energy and gradient
  InitializeEnergy(sysUV, tj, tp, 0);
  AtomForces(crd, CG, tp, dcinp, Etab, EHtab, rcinp, PPk, sysUV, etimers,
	     tj);

  // Copy cell grid information into contiguous arrays
  xpos = (double*)malloc(3*tp->natom*sizeof(double));
  xgrad = (double*)malloc(3*tp->natom*sizeof(double));
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < C->nr[0]; j++) {
      g3con = C->data[j].atmid * 3;
      for (k = 0; k < 3; k++) {
	xpos[g3con+k] = C->data[j].loc[k];
	xgrad[g3con+k] = -1.0*C->data[j].frc[k];
      }
    }
  }

  // Loop calling xminC, then performing tasks as appropriate
  XminDone = 0;
  XminErr = 0;
  while (XminDone == 0 && XminErr == 0) {

    // The arguments to xminC() are as follows:                
    //   xyz_min       :                                       
    //   minim_method  : the energy minimization method        
    //   maxiter       : the maximum number of iterations      
    //   grms_tol      : the gradient convergence tolerance    
    //   natm_ext      : the number of atoms in the system     
    //                   (multiplied by three)                 
    //   m_lbfgs       : depth of the LBFGS memory for LBFGS   
    //                   LBFGS minimization or preconditioning;
    //                   a value of 0 turns off preconditioning
    //   numdiff       : XMIN's finite difference Hv matrix-   
    //                   vector product method                 
    //   xyz_ext       : the coordinates                       
    //   enrg          : the energy of the system              
    //   grad_ext      : the gradient, negative of the forces  
    //   grms          : RMS deviation of the gradient relative
    //                   to the previous step                  
    //   iter:         : the number of this iteration          
    //   total_time    : the total time taken by minimization  
    //   print_level   : the level of verbosity                
    //   ls_method     :                                       
    //   ls_maxiter    :                                       
    //   ls_iter       :                                       
    //   ls_maxatomv   :                                       
    //   beta_armijo   : these are parameters for the Armijo   
    //                   optimization protocol                 

    MinE = xminC();
    if (xflg == DONE) {
      XminDone = 1;
      continue;
    }

    // Copy list coordinates back to cells
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      for (j = 0; j < C->nr[0]; j++) {
	g3con = C->data[j].atmid * 3;
	for (k = 0; k < 3; k++) {
	  C->data[j].loc[k] = xpos[g3con+k];
	}
      }
    }

    // Update cell populations
#ifdef MPI
    UpdateCells(CG, crd, tp, tj);
#else
    UpdateCells(CG, crd, tp);
#endif

    // Compute forces and energy, as needed
    InitializeEnergy(sysUV, tj, tp, 0);
    sysUV->updateU = (xflg == CALCENRG) ? -1 : (xflg == CALCGRAD) ? 0 : 1;
    AtomForces(crd, CG, tp, dcinp, Etab, EHtab, rcinp, PPk, sysUV, etimers,
	       tj);

    // Copy cell grid information into contiguous arrays
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      for (j = 0; j < C->nr[0]; j++) {
	g3con = C->data[j].atmid * 3;
	for (k = 0; k < 3; k++) {
	  xpos[g3con+k] = C->data[j].loc[k];
	  xgrad[g3con+k] = -1.0*C->data[j].frc[k];
	}
      }
    }
  }
}
#endif

//-----------------------------------------------------------------------------
// PrintTriLine: this function just helps me to print out numbers in triplets,
//               line by line.                                
//                                                                      
// Arguments:                                                           
//   val:   the values to print                                         
//   n:     the number of triplets                                      
//   outp:  the file to print values into                               
//-----------------------------------------------------------------------------
static void PrintTriLine(double* val, int n, FILE *outp)
{
  int i;

  for (i = 0; i < n; i++) {
    fprintf(outp, "%20.13e %20.13e %20.13e\n", val[3*i], val[3*i+1],
	    val[3*i+2]);
  }
}

//-----------------------------------------------------------------------------
// PrepForceReport: this function writes a complete report on forces acting on
//                  all atoms.  The breakdown of forces is similar to the
//                  breakdown of energies shown in the diagnostics / state
//                  information output file.        
//
// Arguments:                                                           
//   crd:     the coordinates (for velocity information)                
//   CG:      the cell grid                                             
//   tp:      the topology (for masses)                                 
//   dcinp:   the direct space parameters                               
//   Etab:    the direct space electrostatic interaction spline table   
//   EHtab:   the high-resolution direct space electrostatic interaction
//            spline table                                              
//   rcinp:   the reciprocal space control parameters                   
//   Q:       the reciprocal space charge / potential mesh              
//   PPk:     the reciprocal space pair potential mesh construction kit 
//   sysUV:   state information data (potential energy and virial)      
//-----------------------------------------------------------------------------
void PrepForceReport(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
                     FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
		     trajcon *tj, execon *etimers)
{
  int i, natm;
  double* bondVir;
  double* anglVir;
  double* diheVir;
  double* vdwVir;
  double* relecVir;
  double* delecVir;
  double *dtmp;
  double* bondfrc;
  double* anglfrc;
  double* dihefrc;
  double* relecfrc;
  double* delecfrc;
  double* vdwfrc;
  FILE *outp;
  prmtop tpX;
  Energy sysUV, sysUVt;
  time_t ct;

  // Begin printing the output file
  outp = fopen(tj->dumpname, "w");
  ct = time(&ct);
  fprintf(outp, "%% This is an mdgx force report.  It contains information on "
	  "all types of\n%% forces acting on the system's atoms, much like a "
	  "force dump file written by\n%% SANDER.\n\n%% Compilation date: %s"
	  "\n%% Input file text:\n\n", asctime(localtime(&ct)));
  for (i = 0; i < tj->inptext.row; i++) {
    fprintf(outp, "%% %s", tj->inptext.map[i]);
  }
  fprintf(outp, "\n");
  fclose(outp);

  // Forces, energies, and virials are all wanted
  InitializeEnergy(&sysUV, tj, tp, 1);
  InitializeEnergy(&sysUVt, tj, tp, 1);
  sysUV.updateU = 2;
  sysUVt.updateU = 2;

  // Make a copy of the topology
  tpX = CopyTopology(tp);
  natm = tp->natom;
  dtmp = crd->frc;

  // Remove nonbonded interactions from consideration
  SetDVec(tpX.Charges, natm, 0.0);
  SetDVec(tpX.LJftab.data, 2*tp->ntypes*tp->ntypes, 0.0);
  SetDVec(tpX.LJutab.data, 2*tp->ntypes*tp->ntypes, 0.0);

  // Compute bonded interactions
  bondfrc = (double*)calloc(3*natm, sizeof(double));
  crd->frc = bondfrc;
  for (i = 0; i < tpX.nBAH.ndihe; i++) {
    tpX.HParam[i].K = 0.0;
  }
  for (i = 0; i < tpX.nBAH.nangl; i++) {
    tpX.AParam[i].K = 0.0;
  }
  AtomsToCells(crd, CG, &tpX);
  ZeroCellForces(CG);
#ifdef MPI
  ShareCoordinates(CG, tj, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG, tj);
#else
  ShareCoordinates(CG, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG);
#endif
  for (i = 0; i < CG->ng[0]*CG->ng[1]*CG->ng[2]; i++) {
    CellBondedIntr(&CG->data[i], CG, crd, &tpX, Etab, EHtab, &sysUVt,
		   rcinp->qform);
  }
#ifdef MPI
  MergeCellForces(CG, crd, tp, tj);
#else
  MergeCellForces(CG, crd, tp);
#endif
  MapCellForcesToAtoms(CG, crd);
  sysUV.bond = sysUVt.bond;
  bondVir = CpyDVec(sysUVt.Vir, 9);
  SetDVec(sysUVt.Vir, 9, 0.0);

  // Compute angle interactions
  anglfrc = (double*)calloc(3*natm, sizeof(double));
  crd->frc = anglfrc;
  for (i = 0; i < tpX.nBAH.nbond; i++) {
    tpX.BParam[i].K = 0.0;
  }
  for (i = 0; i < tpX.nBAH.nangl; i++) {
    tpX.AParam[i].K = tp->AParam[i].K;
  }
  AtomsToCells(crd, CG, &tpX);
  ZeroCellForces(CG);
#ifdef MPI
  ShareCoordinates(CG, tj, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG, tj);
#else
  ShareCoordinates(CG, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG);
#endif
  for (i = 0; i < CG->ng[0]*CG->ng[1]*CG->ng[2]; i++) {
    CellBondedIntr(&CG->data[i], CG, crd, &tpX, Etab, EHtab, &sysUVt,
		   rcinp->qform);
  }
#ifdef MPI
  MergeCellForces(CG, crd, tp, tj);
#else
  MergeCellForces(CG, crd, tp);
#endif
  MapCellForcesToAtoms(CG, crd);
  sysUV.angl = sysUVt.angl;
  anglVir = CpyDVec(sysUVt.Vir, 9);
  SetDVec(sysUVt.Vir, 9, 0.0);

  // Compute dihedral interactions
  dihefrc = (double*)calloc(3*natm, sizeof(double));
  crd->frc = dihefrc;
  for (i = 0; i < tpX.nBAH.ndihe; i++) {
    tpX.HParam[i].K = tp->HParam[i].K;
  }
  for (i = 0; i < tpX.nBAH.nangl; i++) {
    tpX.AParam[i].K = 0.0;
  }
  AtomsToCells(crd, CG, &tpX);
  ZeroCellForces(CG);
#ifdef MPI
  ShareCoordinates(CG, tj, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG, tj);
#else
  ShareCoordinates(CG, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG);
#endif
  for (i = 0; i < CG->ng[0]*CG->ng[1]*CG->ng[2]; i++) {
    CellBondedIntr(&CG->data[i], CG, crd, &tpX, Etab, EHtab, &sysUVt,
		   rcinp->qform);
  }
#ifdef MPI
  MergeCellForces(CG, crd, tp, tj);
#else
  MergeCellForces(CG, crd, tp);
#endif
  MapCellForcesToAtoms(CG, crd);
  sysUV.dihe = sysUVt.dihe;
  diheVir = CpyDVec(sysUVt.Vir, 9);
  SetDVec(sysUVt.Vir, 9, 0.0);

  // Reciprocal space force calculation
  relecfrc = (double*)calloc(3*natm, sizeof(double));
  crd->frc = relecfrc;
  for (i = 0; i < natm; i++) {
    tpX.Charges[i] = tp->Charges[i];
  }
  AtomsToCells(crd, CG, &tpX);
  ZeroCellForces(CG);
#ifdef MPI
  ShareCoordinates(CG, tj, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG, tj);
#else
  ShareCoordinates(CG, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG);
#endif
  i = (rcinp->QL[0].pfft == 1) ? 2*(rcinp->QL[0].col/2+1) : rcinp->QL[0].col;
  SetDVec(rcinp->QL[0].data, rcinp->QL[0].pag*rcinp->QL[0].row*i, 0.0);
  for (i = 0; i < CG->ng[0]*CG->ng[1]*CG->ng[2]; i++) {
    CellSplCoeff(&CG->data[i], crd, rcinp);
    CellQ2Book(&CG->data[i], rcinp, &rcinp->QL[0]);
  }
  if (rcinp->nlev == 1) {
    ConvQBCnrgvir(rcinp, crd, &rcinp->QL[0], PPk, &sysUVt, etimers);
  }
  else {
    mleConvQBC(rcinp, &sysUVt, PPk, etimers);
  }
  for (i = 0; i < CG->ng[0]*CG->ng[1]*CG->ng[2]; i++) {
    CellIntrpFrc(&CG->data[i], crd, rcinp, &rcinp->QL[0]);
  }
#ifdef MPI
  MergeCellForces(CG, crd, tp, tj);
#else
  MergeCellForces(CG, crd, tp);
#endif
  MapCellForcesToAtoms(CG, crd);
  sysUV.relec = sysUVt.relec;
  relecVir = CpyDVec(sysUVt.Vir, 9);
  SetDVec(sysUVt.Vir, 9, 0.0);

  // Direct space electrostatic calculation
  delecfrc = (double*)calloc(3*natm, sizeof(double));
  crd->frc = delecfrc;
  for (i = 0; i < tpX.nBAH.ndihe; i++) {
    tpX.HParam[i].K = 0.0;
  }
  AtomsToCells(crd, CG, &tpX);
  ZeroCellForces(CG);
#ifdef MPI
  ShareCoordinates(CG, tj, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG, tj);
#else
  ShareCoordinates(CG, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG);
#endif
  for (i = 0; i < CG->ng[0]*CG->ng[1]*CG->ng[2]; i++) {
    DirectTriple2(&CG->data[i], crd, Etab, dcinp, &tpX, &sysUVt);
    CellBondedIntr(&CG->data[i], CG, crd, &tpX, Etab, EHtab, &sysUVt,
		   rcinp->qform);
  }
#ifdef MPI
  MergeCellForces(CG, crd, tp, tj);
#else
  MergeCellForces(CG, crd, tp);
#endif
  MapCellForcesToAtoms(CG, crd);
  sysUV.delec = sysUVt.delec;
  delecVir = CpyDVec(sysUVt.Vir, 9);
  SetDVec(sysUVt.Vir, 9, 0.0);

  // Direct space vdW force calculation
  vdwfrc = (double*)calloc(3*natm, sizeof(double));
  crd->frc = vdwfrc;
  SetDVec(tpX.Charges, natm, 0.0);
  for (i = 0; i < 2*tp->ntypes*tp->ntypes; i++) {
    tpX.LJftab.data[i] = tp->LJftab.data[i];
    tpX.LJutab.data[i] = tp->LJutab.data[i];
  }
  AtomsToCells(crd, CG, &tpX);
  ZeroCellForces(CG);
#ifdef MPI
  ShareCoordinates(CG, tj, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG, tj);
#else
  ShareCoordinates(CG, dcinp->Mcut, crd, &tpX, 0);
  ExtraPointLocations(tp, crd, CG);
#endif
  for (i = 0; i < CG->ng[0]*CG->ng[1]*CG->ng[2]; i++) {
    DirectTriple2(&CG->data[i], crd, Etab, dcinp, &tpX, &sysUVt);
    CellBondedIntr(&CG->data[i], CG, crd, &tpX, Etab, EHtab, &sysUVt,
		   rcinp->qform);
  }
#ifdef MPI
  MergeCellForces(CG, crd, tp, tj);
#else
  MergeCellForces(CG, crd, tp);
#endif
  MapCellForcesToAtoms(CG, crd);
  sysUV.vdw12 = sysUVt.vdw12;
  sysUV.vdw6 = sysUVt.vdw6;
  vdwVir = CpyDVec(sysUVt.Vir, 9);
  SetDVec(sysUVt.Vir, 9, 0.0);

  // Free the extra allocated topology, sum forces and virials
  FreeTopology(&tpX);
  DVec2VecAdd(sysUV.Vir, bondVir, 9);
  DVec2VecAdd(sysUV.Vir, anglVir, 9);
  DVec2VecAdd(sysUV.Vir, diheVir, 9);
  DVec2VecAdd(sysUV.Vir, delecVir, 9);
  DVec2VecAdd(sysUV.Vir, relecVir, 9);
  DVec2VecAdd(sysUV.Vir, vdwVir, 9);
  for (i = 0; i < 3*natm; i++) {
    dtmp[i] = bondfrc[i] + anglfrc[i] + dihefrc[i] + relecfrc[i] + 
      delecfrc[i] + vdwfrc[i];
  }
  crd->frc = dtmp;

  // Write all forces and energies
  outp = fopen(tj->dumpname, "a");
  fprintf(outp, "%s.Ubond = %20.13e; %% Bond interaction energy\n",
	  tj->DMPvar, sysUV.bond);
  fprintf(outp, "%s.Uangl = %20.13e; %% Angle interaction energy\n",
	  tj->DMPvar, sysUV.angl);
  fprintf(outp, "%s.Udihe = %20.13e; %% Dihedral interaction energy\n",
	  tj->DMPvar, sysUV.dihe);
  fprintf(outp, "%s.Udir  = %20.13e; %% Direct space sum electrostatic "
	  "energy\n", tj->DMPvar, sysUV.delec);
  fprintf(outp, "%s.Urec  = %20.13e; %% Reciprocal space sum "
	  "electrostatic energy\n", tj->DMPvar, sysUV.relec);
  fprintf(outp, "%s.Uvdw  = %20.13e; %% van-der Waals energy\n", tj->DMPvar,
	  sysUV.vdw12 + sysUV.vdw6);
  if (tj->DMPcrd == 1) {
    fprintf(outp, "\n%% Coordinates:\n%s.crd = [\n", tj->DMPvar);
    PrintTriLine(crd->loc, tp->natom, outp);
    fprintf(outp, "];\n");
    fprintf(outp, "\n%% Box dimensions (X, Y, Z, followed by a, b, g in "
	    "radians):\n%s.box = [\n", tj->DMPvar);
    PrintTriLine(crd->gdim, 2, outp);
    fprintf(outp, "];\n");
  }
  if (tj->DMPbond == 1) {
    fprintf(outp, "\n%% Virial contributions from bonds\n%s.bond_vir = [\n",
	    tj->DMPvar);
    PrintTriLine(bondVir, 3, outp);
    fprintf(outp, "];\n");
    fprintf(outp, "\n%% Forces due to bonds\n%s.bond_frc = [\n", tj->DMPvar);
    PrintTriLine(bondfrc, tp->natom, outp);
    fprintf(outp, "];\n");
  }
  if (tj->DMPangl == 1) {
    fprintf(outp, "\n%% Virial contributions from angles\n%s.angl_vir = [\n",
	    tj->DMPvar);
    PrintTriLine(anglVir, 3, outp);
    fprintf(outp, "];\n");
    fprintf(outp, "\n%% Forces due to angles\n%s.angl_frc = [\n", tj->DMPvar);
    PrintTriLine(anglfrc, tp->natom, outp);
    fprintf(outp, "];\n");
  }
  if (tj->DMPdihe == 1) {
    fprintf(outp, "\n%% Virial contributions from dihedrals\n%s.dihe_vir = "
	    "[\n", tj->DMPvar);
    PrintTriLine(diheVir, 3, outp);
    fprintf(outp, "];\n");
    fprintf(outp, "\n%% Forces due to dihedrals\n%s.dihe_frc = [\n",
	    tj->DMPvar);
    PrintTriLine(dihefrc, tp->natom, outp);
    fprintf(outp, "];\n");
  }
  if (tj->DMPdelec == 1) {
    fprintf(outp, "\n%% Virial contributions from direct space electrostatics"
	    "\n%s.delec_vir = [\n", tj->DMPvar);
    PrintTriLine(delecVir, 3, outp);
    fprintf(outp, "];\n");
    fprintf(outp, "\n%% Forces due to electrostatic direct space sum\n"
	    "%s.delec_frc = [\n", tj->DMPvar);
    PrintTriLine(delecfrc, tp->natom, outp);
    fprintf(outp, "];\n");
  }
  if (tj->DMPrelec == 1) {
    fprintf(outp, "\n%% Virial contributions from reciprocal space "
	    "electrostatics\n%s.relec_vir = [\n", tj->DMPvar);
    PrintTriLine(relecVir, 3, outp);
    fprintf(outp, "];\n");
    fprintf(outp, "\n%% Forces due to electrostatic reciprocal space sum\n"
	    "%s.relec_frc = [\n", tj->DMPvar);
    PrintTriLine(relecfrc, tp->natom, outp);
    fprintf(outp, "];\n");
  }
  if (tj->DMPvdw == 1) {
    fprintf(outp, "\n%% Virial contributions from van der Waals interactions"
	    "\n%s.vdw_vir = [\n", tj->DMPvar);
    PrintTriLine(vdwVir, 3, outp);
    fprintf(outp, "];\n");
    fprintf(outp, "\n%% Forces due to van-der Waals interactions\n"
	    "%s.vdw_frc = [\n", tj->DMPvar);
    PrintTriLine(vdwfrc, tp->natom, outp);
    fprintf(outp, "];\n");
  }
  if (tj->DMPall == 1) {
    fprintf(outp, "\n%% Total system virial\n%s.sum_vir = [\n", tj->DMPvar);
    PrintTriLine(sysUV.Vir, 3, outp);
    fprintf(outp, "];\n");
    fprintf(outp, "\n%% Total force on all atoms\n"
	    "%s.sum_frc = [\n", tj->DMPvar);
    PrintTriLine(crd->frc, tp->natom, outp);
    fprintf(outp, "];\n");
  }

  // Close up the output file
  fclose(outp);

  // Free allocated memory
  free(bondfrc);
  free(anglfrc);
  free(dihefrc);
  free(relecfrc);
  free(delecfrc);
  free(vdwfrc);
  free(bondVir);
  free(anglVir);
  free(diheVir);
  free(vdwVir);
  free(relecVir);
  free(delecVir);
  DestroyEnergyTracker(&sysUV);
  DestroyEnergyTracker(&sysUVt);
}
