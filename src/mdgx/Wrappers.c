#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "CellManip.h"
#include "CompFrc.h"
#include "CrdManip.h"
#include "Integrator.h"
#include "MPIMap.h"
#include "MPITypeCast.h"
#include "pmeDirect.h"
#include "pmeRecip.h"
#include "Thermostats.h"
#include "Topology.h"
#include "Trajectory.h"

//-----------------------------------------------------------------------------
// uform: all the potential function structs wrapped up in one type.    
//-----------------------------------------------------------------------------
struct PotentialFunction {
  prmtop tp;        // Topology
  dircon dcinp;     // Direct-space controls
  FrcTab Etab;      // Direct-space standard (coarse) lookup table
  FrcTab EHtab;     // Direct-space high-resolution lookkup table
  reccon rcinp;     // Reciprocal space controls
  bckit PPk;        // Convolution support data
};
typedef struct PotentialFunction uform;

//-----------------------------------------------------------------------------
// mdsys: all structs for a molecular dynamics trajectory in one type.  
//-----------------------------------------------------------------------------
struct MolecularDynamicsSystem {
  coord crd;       // Coordinates (phone book organization of the system)      
  cellgrid CG;     // Cell grid (spatial reorganization of the system)         
  Energy sysUV;    // System energy decomposition and virial                   
  execon etimers;  // Timing data (for this system only)
};
typedef struct MolecularDynamicsSystem mdsys;

//-----------------------------------------------------------------------------
// InitBasicDircon: initialize a basic direct-space control structure.  
//                  Allows the user to set the direct space cutoff.     
//                                                                      
// Arguments:                                                           
//   NBcut:    the non-bonded cutoff (for nuances involving different   
//             electrostatic and van-der Waals cutoffs, users can       
//             manually edit the struct after this routine)             
//-----------------------------------------------------------------------------
void InitBasicDircon(dircon *dcinp, double NBcut)
{
  dcinp->LRvdw = 1;
  dcinp->Ecut = NBcut;
  dcinp->Vcut = NBcut;
  dcinp->Mcut = NBcut;
  dcinp->MaxDens = 2.5;
  dcinp->invMcut = 1.0/dcinp->Mcut;
  dcinp->invEcut = 1.0/dcinp->Ecut;
  dcinp->Dtol = 1.0e-5;
  dcinp->ewcoeff = EwaldCoefficient(dcinp->Ecut, dcinp->Dtol);
  dcinp->sigma = 0.5/dcinp->ewcoeff;
  dcinp->lkpspc = 0.0625;
}

//-----------------------------------------------------------------------------
// InitBasicReccon: initialize a basic reciprocal-space control data    
//                  structure.  The default ewald parameters are enforced.
//                                                                      
// Arguments:                                                           
//   rcinp:    the reciprocal space control data struct                 
//   dcinp:    the direct space control data struct (for sigma)         
//-----------------------------------------------------------------------------
void InitBasicReccon(reccon *rcinp, dircon *dcinp)
{
  rcinp->nlev = 1;
  rcinp->ordr[0] = 4;
  rcinp->ordr[1] = 4;
  rcinp->ordr[2] = 4;
  rcinp->ng = (int*)malloc(3*sizeof(int));
  rcinp->ng[0] = -1;
  rcinp->ng[1] = -1;
  rcinp->ng[2] = -1;
  rcinp->S = dcinp->sigma;
}

//-----------------------------------------------------------------------------
// PrimeBasicTopology: this prepares a topology struct to be compatible with
//                     functions downstream in the flow of data.  This is not
//                     labeled an "Init" function because it does not precisely
//                     "initialize" a topology; there is no "default" topology
//                     with a fixed number of atoms or properties, but there
//                     are default values to some of the attributes.
//
// Arguments:                                                           
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
void PrimeBasicTopology(prmtop *tp)
{
  sprintf(tp->source, "prmtop");
  sprintf(tp->WaterName, "WAT ");
  tp->ljbuck = 0;
  tp->eprulesource[0] = '\0';
  tp->settle = 0;
  tp->rattle = 0;
  tp->lj14fac = 2.0;
  tp->elec14fac = 1.2;
}

//-----------------------------------------------------------------------------
// InitBasicTrajcon: initialize a trajectory control data structure.  This
//                   routine permits no user intervention at the moment.
//
// Arguments:                                                           
//   tj:     trajectory control data, just gets filled up with lots of  
//           default values to make other things able to run            
//-----------------------------------------------------------------------------
void InitBasicTrajcon(trajcon *tj)
{
#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &tj->tid);
  MPI_Comm_size(MPI_COMM_WORLD, &tj->nthreads);
  DefineMPITypes(&tj);
#else
  tj->tid = 0;
  tj->nthreads = 1;
#endif
  tj->mode = 0;
  tj->starttime = 0.0;
  tj->dt = 0.001;
  tj->Ptarget = 1.0;
  tj->rattletol = 1.0e-6;
  tj->MaxRattleIter = 100;
  tj->RemoveMomentum = 0;
  tj->ioutfm = 0;
  tj->ntt = 0;
  tj->ntp = 0;
  tj->barostat = 1;
  tj->vrand = 1000;
  tj->BerendsenTCoupl = 0.4;
  tj->BerendsenPTime = 1.0;
  tj->BerendsenPCoupl = 44.6;
  tj->MCBarostatFac[0] = 2.0e-3;
  tj->MCBarostatFac[1] = -1.0;
  tj->MCBarostatFac[2] = -1.0;
  tj->MCBarostatFreq = 100;
  tj->npth.TauT = 1.0;
  tj->npth.TauP = 1.0;
  tj->lnth.gamma_ln = 0.0;
  tj->nstep = 1;
  tj->nfistep = 0;
  tj->currfi = 0;
  tj->ntwr = 0;
  tj->ntwx = 0;
  tj->ntwv = 0;
  tj->ntwf = 0;
  tj->ntpr = 0;
  tj->irest = 0;
  tj->topchk = 1;
  tj->TI = 0;
  tj->mxorder = 1;
  tj->nsynch = 1000;
  tj->nsys = 1;
  tj->Leash.active = 0;

  // Initialize larger structs within the trajcon
  SelectSystemsToTend(tj);
}

//-----------------------------------------------------------------------------
// InitPotential: initialize real space components of the potential function,
//                including default direct space controls and force lookup
//                tables.                                  
//                                                                      
// Arguments:                                                           
//   topsrc:   the topology file source name (if NULL, uses "prmtop")   
//   NBcut:    the non-bonded cutoff                                    
//   tj:       trajectoty control data                                  
//-----------------------------------------------------------------------------
uform InitPotential(char* topsrc, double NBcut, trajcon *tj)
{
  double spc, ewc;
  uform thisU={0};

  PrimeBasicTopology(&thisU.tp);
  if (topsrc[0] != '\0') {
    strcpy(thisU.tp.source, topsrc);
  }
  
  GetPrmTop(&thisU.tp, tj, 1);
  
  InitBasicDircon(&thisU.dcinp, NBcut);
  InitBasicReccon(&thisU.rcinp, &thisU.dcinp);
  LongRangeVDW(&thisU.tp, &thisU.dcinp);
  spc = thisU.dcinp.lkpspc;
  ewc = thisU.dcinp.ewcoeff;
  thisU.Etab = DirectSpaceR2(NBcut, spc, ewc, 0, 0);
  thisU.EHtab = DirectSpaceR2(NBcut, spc, ewc, 0, 0);

  return thisU;
}

//-----------------------------------------------------------------------------
// LoadCoordToGrid: load coordinate information from a file (inpcrd or restart
//                  format), then create and load a cell grid for the atoms. 
//                  This routine will image all bonded groups, initialize
//                  prior coordinates and velocities, create the required
//                  reciprocal space support structures, and link the cell
//                  grid.                 
//                                                                      
// Arguments:                                                           
//   crdname:    the name of the coordinate file (input)                
//   U:          the potential function                                 
//   tj:         trajectory control, ensemble information (input)       
//-----------------------------------------------------------------------------
mdsys LoadCoordToGrid(char* crdname, uform *U, trajcon *tj)
{
  double tstamp;
  mdsys thisMD;

  // Read coordinates from disk
  tstamp = 0.0;
  thisMD.crd = (crdname[0] != '\0') ?
    ReadRst(&U->tp, crdname, &tstamp) : ReadRst(&U->tp, "inpcrd", &tstamp);
  ImageBondedGroups(&thisMD.crd, &U->tp);
  InitHistory(&thisMD.crd);

  // Create the cell grid and prepare reciprocal space support
  thisMD.CG = CreateCellGrid(&thisMD.crd, &U->dcinp, &U->rcinp, &U->tp, tj, 0);
  PrepPME(&thisMD.CG, &U->rcinp, &thisMD.crd);
  U->PPk = CreateBCKit(&U->rcinp, &U->rcinp.QL[0], &thisMD.crd, &U->tp,
		       FFTW_ESTIMATE);
#ifdef MPI
  LinkCellGrid(&thisMD.CG, &thisMD.crd, &U->rcinp);
#else
  LinkCellGrid(&thisMD.CG, &U->rcinp);
#endif

  // Load the cell grid
  AtomsToCells(&thisMD.crd, &thisMD.CG, &U->tp);

  return thisMD;
}

//-----------------------------------------------------------------------------
// MMForceEnergy: routine for computing the force and energy of a set of
//                coordinates given the topology and Hamiltonian.       
//                                                                      
// Arguments:                                                           
//   U:      the description of the Hamiltonian                         
//   MD:     the description of the system, including coordinates and   
//           energy                                                     
//   tj:     trajectory control data, essentially just passed down to   
//           other functions                                            
//-----------------------------------------------------------------------------
void MMForceEnergy(uform *U, mdsys *MD, trajcon *tj)
{
  InitializeEnergy(&MD->sysUV, tj, &U->tp, 1);
  MD->sysUV.updateU = 1;
  AtomForces(&MD->crd, &MD->CG, &U->tp, &U->dcinp, &U->Etab, &U->EHtab,
	     &U->rcinp, &U->PPk, &MD->sysUV, &MD->etimers, tj);
  MD->sysUV.kine = KineticEnergy(&MD->CG, &MD->crd, &U->tp, tj);
#ifdef MPI
  SumTotalEnergy(&MD->CG, &MD->sysUV);
#else
  SumTotalEnergy(&MD->sysUV);
#endif
}
