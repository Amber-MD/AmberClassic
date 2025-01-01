#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "Random.h"
#include "mdgxVector.h"
#include "Thermostats.h"
#include "CellManip.h"

#include "TrajectoryDS.h"
#include "TopologyDS.h"
#include "CrdManipDS.h"

#define NEED_TI 0
#include "ThermostatsBranch.c"
#undef NEED_TI

#define NEED_TI 1
#include "ThermostatsBranch.c"
#undef NEED_TI

//-----------------------------------------------------------------------------
// AndersenThermostat: a routine for re-assigning velocities from a     
//                     constant-temperature bath.  The TI variant of the
//                     function also assigns velocities to atoms in the second
//                     system either by copying from the original system (if
//                     there is a corresponding atom) or by computing new
//                     random numbers (if the atom is unique).
//                                                                      
// Arguments:                                                           
//   CG:      the cell grid                                             
//   crd:     the coordinates (cells reference coordinates for velocity 
//            information)                                              
//   tp:      the topology (for masses)                                 
//   tj:      trajectory control information                            
//   T:       the temperature target                                    
//-----------------------------------------------------------------------------
void AndersenThermostat(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj,
                        double T)
{
  int i, j, g3con;
  long counter;
  double efac;
  double *vtmp, *vntmp;
  cell *C;

  counter = tj->rndcon;
  const double ebeta = (T > 1.0e-8) ? sqrt(GASCNST*T) : 0.0;
  vtmp = crd->vel;
  if (tj->TI == 1) {
    vntmp = crd[1].vel;
  }

  // Special case: synchronize pseudo-random sequences
  if (tj->SyncRNG == 1) {
    for (i = 0; i < tp->natom; i++) {
      efac = ebeta*sqrt(tp->InvMasses[i]);
      vtmp[3*i] = efac*GaussBoxMuller(&counter);
      vtmp[3*i+1] = efac*GaussBoxMuller(&counter);
      vtmp[3*i+2] = efac*GaussBoxMuller(&counter);
    }
    if (tj->TI == 1 && tj->prc.uniB > 0) {
      for (i = 0; i < tp[1].natom; i++) {
        if (tj->prc.matchB[i] < 0) {
          efac = ebeta*sqrt(tp->InvMasses[i]);
          vntmp[3*i] = efac*GaussBoxMuller(&counter);
          vntmp[3*i+1] = efac*GaussBoxMuller(&counter);
          vntmp[3*i+2] = efac*GaussBoxMuller(&counter);
          j++;
        }
      }
    }
  }

  // Loop over all atoms in the first 
  // cell grid and reassign velocities
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    if (tj->SyncRNG == 0) {
      for (j = 0; j < C->nr[0]; j++) {
        efac = ebeta*sqrt(tp->InvMasses[C->data[j].id]);
        g3con = 3*C->data[j].id;
        vtmp[g3con] = efac*GaussBoxMuller(&counter);
        vtmp[g3con+1] = efac*GaussBoxMuller(&counter);
        vtmp[g3con+2] = efac*GaussBoxMuller(&counter);
      }
    }
    if (tj->TI == 0) {
      continue;
    }

    // If we're still here, TI is implemented
    int k, g3ncon;
    cell *Cn = &CG[1].data[CG->MyCellDomain[i]];
    prmcorr *prc;
    prc = &tj->prc;

    // Simplest case: no unique atoms in either topology
    // and perfect correspondence between atoms         
    if (prc->relate == 0) {
      for (j = 0; j < Cn->nr[0]; j++) {
        g3con = 3*C->data[j].id;
        vntmp[g3con] = vtmp[g3con];
        vntmp[g3con+1] = vtmp[g3con+1];
        vntmp[g3con+2] = vtmp[g3con+2];
      }
    }
    else {

      // Assign new velocities for unique atoms,   
      // unless it has already been done as part   
      // of the synchronized random number protocol
      if (tj->SyncRNG == 0) {
        for (j = 0; j < Cn->nr[0]; j++) {
          if (prc->matchB[Cn->data[j].id] < 0) {
            efac = ebeta*sqrt(tp[1].InvMasses[Cn->data[j].id]);
            g3con = 3*Cn->data[j].id;
            vntmp[g3con] = efac*GaussBoxMuller(&counter);
            vntmp[g3con+1] = efac*GaussBoxMuller(&counter);
            vntmp[g3con+2] = efac*GaussBoxMuller(&counter);
          }
        }
      }

      // Not-so-bad case: frame shift in atom correspondence
      if (prc->relate == 1) {
        k = 0;
        for (j = 0; j < C->nr[0]; j++) {
          if (prc->matchA[C->data[j].id] < 0) {
            continue;
          }
          while (k < Cn->nr[0] && prc->matchB[Cn->data[k].id] < 0) {
            k++;
          }
          g3con = 3*C->data[j].id;
          g3ncon = 3*Cn->data[k].id;
          vntmp[g3ncon] = vtmp[g3con];
          vntmp[g3ncon+1] = vtmp[g3con+1];
          vntmp[g3ncon+2] = vtmp[g3con+2];
          k++;
        }
      }

      // Really bad case: no order in the atom correspondence
      else {
        for (j = 0; j < C->nr[0]; j++) {
          if (prc->matchA[C->data[j].id] >= 0) {
            k = Cn->GPSptr[prc->matchA[C->data[j].id]];
            g3con = 3*C->data[j].id;
            g3ncon = 3*Cn->data[k].id;
            vntmp[g3ncon] = vtmp[g3con];
            vntmp[g3ncon+1] = vtmp[g3con+1];
            vntmp[g3ncon+2] = vtmp[g3con+2];
          }
        }
      }
    }
  }

  // On multiple processors this could result in serious      
  // errors.  The pseudo-random number sequence must either be
  // advanced by a set amount for each processor, or entirely 
  // desynchronized in order for this to work.                
  tj->rndcon = counter;
}

//-----------------------------------------------------------------------------
// CellNHThermostat: a compact function for a nasty series of expressions.
//                   This implements the repeatedly used thermostat in the
//                   Nose-Hoover thermo- (baro-)stat above.  Note that if there
//                   is no barsotat (barostat is 0), then this works on a 1/4
//                   timestep interval.  If there is, it works on a 1/8
//                   timestep interval.  
//                                                                      
// Arguments:                                                           
//   CG:       the cell grid                                            
//   crd:      the coordinates (cells reference coordinates for velocity
//             information)                                             
//   tp:       the topology                                             
//   tj:       trajectory control information                           
//   sysUV:    the system energy and virial                             
//   ETAt:     scaling factor relating to the barostat piston mass      
//   CHIt:     (modified and returned) exponential argument for velocity
//             rescaling, similar to chi in the Berendsen thermostat  
//   CHItp1e:  CHIt advanced 1/8 timestep                               
//   CHItp1q:  CHIt advanced 1/4 timestep                               
//   barostat: flag to indicate whether a Nose-Hoover barostat is active
//-----------------------------------------------------------------------------
void NHThermostat(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj,
		  Energy *sysUV, double ETAt, double *CHIt, double *CHItp1e,
		  double *CHItp1q, int barostat)
{
  int i, j, g3con;
  double dt, Ttarget, pmass, qmass, invqmass, sigma, vfac;
  double *vtmp;
  cell *C;

  // Unpack
  dt = (barostat == 1) ? 0.125*sqrt(418.4)*tj->dt : 0.25*sqrt(418.4)*tj->dt;
  Ttarget = (barostat == 1) ? GASCNST*tj->Ttarget : 0.0;
  pmass = (barostat == 1) ? ETAt*ETAt*tj->npth.pmass : 0.0;
  qmass = tj->npth.qmass;
  sigma = 2.0*tj->npth.sigma;
  invqmass = 1.0/qmass;

  // Compute the kinetic energy
  sysUV->kine = KineticEnergy(CG, crd, tp, tj);
  *CHItp1e = *CHIt + dt*(2.0*sysUV->kine + pmass - sigma - Ttarget)*invqmass;
  const int ncell = CG->ng[0]*CG->ng[1]*CG->ng[2];
  vtmp = crd->vel;
  vfac = exp(-2.0*dt*(*CHItp1e));
  for (i = 0; i < ncell; i++) {
    C = &CG->data[i];
    for (j = 0; j < C->nr[0]; j++) {
      g3con = 3*C->data[j].id;
      vtmp[g3con] *= vfac;
      vtmp[g3con+1] *= vfac;
      vtmp[g3con+2] *= vfac;
    }
  }
  sysUV->kine = KineticEnergy(CG, crd, tp, tj);
  *CHItp1q = *CHItp1e +
    dt*(2.0*sysUV->kine + pmass - sigma - Ttarget)*invqmass;
}

//-----------------------------------------------------------------------------
// PrepThermoBarostat: prepare constants for a thermostat and barostat, if
//                     needed.                                       
//                                                                      
// Arguments:                                                           
//   tp:      system topology                                           
//   tj:      trajectory control data (for lack of a better place to store
//            these constants)                                    
//-----------------------------------------------------------------------------
void PrepThermoBarostat(prmtop *tp, trajcon *tj)
{
  tj->npth.sigma = 0.5*tp->ndf*GASCNST*tj->Ttarget;
  tj->npth.qmass = 2.0*tj->npth.sigma*tj->npth.TauT*tj->npth.TauT;
  tj->npth.pmass = (tp->ndf + 3)*GASCNST*tj->Ttarget*pow(tj->npth.TauP, 2.0);
  tj->npth.chi = 0.0;
  tj->npth.eta = 0.0;
}

//-----------------------------------------------------------------------------
// PrepLangevinThermostat: prepare constants for a Langevin thermostat. 
//                                                                      
// Arguments:                                                           
//   tj:      trajectory control data (for lack of a better place to store
//            these constants)                                    
//-----------------------------------------------------------------------------
void PrepLangevinThermostat(trajcon *tj)
{
  double hdt, gammai;
  
  hdt = 0.5*sqrt(418.4)*tj->dt;
  gammai = tj->lnth.gamma_ln / sqrt(418.4);
  tj->lnth.c_implic = 1.0 / (1.0 + gammai*hdt);
  tj->lnth.c_explic = 1.0 - gammai*hdt;
  tj->lnth.sdfac = sqrt(gammai*GASCNST*tj->Ttarget/hdt);
}
