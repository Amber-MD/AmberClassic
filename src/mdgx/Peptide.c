#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "Topology.h"
#include "Trajectory.h"
#include "CrdManip.h"
#include "Manual.h"
#include "Matrix.h"
#include "Random.h"
#include "Timings.h"
#include "ParamFit.h"
#include "IPolQ.h"
#include "Peptide.h"
#include "Gpu.h"
#include "mdgxVector.h"
#include "MultiSimDS.h"
#ifdef CUDA
#  include "ArraySimulator.cuh"
#endif

//-----------------------------------------------------------------------------
// AssignFileByTopology: assign a name to an MD output or coordinates file
//                       based on the name of the system's topology.
//
// Arguments:
//   ppctrl:   peptide simulations control data (file names, directives for
//             certain enhanced sampling methods)
//   fname:    buffer to hold the file name
//   duptops:  flags to indicate that there are duplicate topologies for
//             multiple systems
//   sidx:     index of the system within ppctrl
//-----------------------------------------------------------------------------
static void AssignFileByTopology(pepcon *ppctrl, char* fname, int* duptops,
                                 int sidx)
{
  int i, slen;

  if (fname[0] == '\0') {
    i = 0;
    slen = strlen(ppctrl->tpinames.map[sidx]);
    while (i < slen && ppctrl->tpinames.map[sidx][i] != '.') {
      fname[i] = ppctrl->tpinames.map[sidx][i];
      i++;
    }
    if (duptops[sidx] > 0) {
      sprintf(&fname[i], "_S%d", duptops[sidx]);
    }
  }
}

//-----------------------------------------------------------------------------
// ExpandSystemIO: expand the IO and basic run conditions for each peptide
//                 simulation.  This is done after taking in all inputs from
//                 the &pptd as well as &cntrl namelists, as both may contain
//                 pertinent information.
//
// Arguments:
//   ppctrl:  peptide simulations control data (file names, directives for
//            certain enhanced sampling methods)
//   tj:      trajectory control data (contains &cntrl namelist information,
//            specifically the temperature if nothing is specified for any
//            given system)
//-----------------------------------------------------------------------------
static void ExpandSystemIO(pepcon *ppctrl, trajcon *tj)
{
  int i, j, k, m, slen, maxsys;
  int* toplen;
  int* duptops;

  // Fill in the &cntrl namelist temperature if none
  // has been provided for any of the systems.  Count
  // the total number of systems.
  maxsys = 0;
  for (i = 0; i < ppctrl->nsys; i++) {
    if (ppctrl->Tranges[i].x < 0.0) {
      ppctrl->Tranges[i].x = tj->Ttarget;
      ppctrl->Tranges[i].y = tj->Ttarget;
    }
    if (ppctrl->Sreplicas[i] > 1) {
      maxsys += ppctrl->Sreplicas[i];
    }
    else {
      maxsys += ppctrl->Treplicas[i] * ppctrl->Preplicas[i];
    }
  }

  // Check that there are not too many independent systems
  if (maxsys > 960 && tj->ioutfm == 1) {
    printf("ExpandSystemIO >> Too many systems (%d) have been requested.  The "
	   "maximum\nExpandSystemIO >> number of NetCDF trajectories that can "
	   "be written at once\nExpandSystemIO >> is 960.  Reduce the number "
	   "of systems (recommended), or\nExpandSystemIO >> switch to .crd "
	   "format (ioutfm = 0 in &cntrl namelist).\n", maxsys);
    exit(1);
  }

  // Fill in other details from the trajectory control data
  ppctrl->rattle = tj->rattle;

  // Detect systems with similar topology names
  toplen = (int*)malloc(ppctrl->nsys * sizeof(int));
  duptops = (int*)calloc(ppctrl->nsys, sizeof(int));
  for (i = 0; i < ppctrl->nsys; i++) {
    j = 0;
    slen = strlen(ppctrl->tpinames.map[i]);
    while (j < slen && ppctrl->tpinames.map[i][j] != '.') {
      j++;
    }
    toplen[i] = j;
  }
  for (i = 0; i < ppctrl->nsys; i++) {
    if (duptops[i] > 0) {
      continue;
    }
    k = 2;
    for (j = i + 1; j < ppctrl->nsys; j++) {
      if (toplen[i] == toplen[j] &&
          strncmp(ppctrl->tpinames.map[i],
                  ppctrl->tpinames.map[j], toplen[i]) == 0) {
        duptops[i] = 1;
        duptops[j] = k;
        k++;
      }
    }
  }

  // Fill in the output file names based on topologies,
  // if none have been supplied.
  for (i = 0; i < ppctrl->nsys; i++) {
    AssignFileByTopology(ppctrl, ppctrl->outbases.map[i], duptops, i);
    AssignFileByTopology(ppctrl, ppctrl->trjbases.map[i], duptops, i);
    AssignFileByTopology(ppctrl, ppctrl->rstrtbases.map[i], duptops, i);
  }

  // Expand the systems if replica exchange series have been specified
  ppctrl->tpinames   = ReallocCmat(&ppctrl->tpinames, maxsys, MAXNAME);
  ppctrl->tpfnames   = ReallocCmat(&ppctrl->tpfnames, maxsys, MAXNAME);
  ppctrl->crdnames   = ReallocCmat(&ppctrl->crdnames, maxsys, MAXNAME);
  ppctrl->outbases   = ReallocCmat(&ppctrl->outbases, maxsys, MAXNAME);
  ppctrl->trjbases   = ReallocCmat(&ppctrl->trjbases, maxsys, MAXNAME);
  ppctrl->rstrtbases = ReallocCmat(&ppctrl->rstrtbases, maxsys, MAXNAME);
  ppctrl->starttime  = (double*)malloc(maxsys * sizeof(double));
  SetDVec(ppctrl->starttime, maxsys, -1.0);
  ppctrl->Tmix = (double*)malloc(maxsys * sizeof(double));
  ppctrl->Pmix = (double*)malloc(maxsys * sizeof(double));
  ppctrl->T = (double*)malloc(maxsys * sizeof(double));
  m = maxsys - 1;
  for (i = ppctrl->nsys - 1; i >= 0; i--) {

    // Standard replicas, which will not exchange and run at the same
    // temperature, take precedence over any sort of REMD.
    if (ppctrl->Sreplicas[i] > 1) {
      for (j = ppctrl->Sreplicas[i] - 1; j >= 0; j--) {
        ppctrl->Tmix[m] = 0.0;
        ppctrl->Pmix[m] = 0.0;
        if (ppctrl->Tranges[i].x >= 0.0) {
          ppctrl->T[m] = ppctrl->Tranges[i].x;
        }
        else {
          ppctrl->T[m] = tj->Ttarget;
        }
        if (m != i) {
          strcpy(ppctrl->tpinames.map[m],   ppctrl->tpinames.map[i]);
          strcpy(ppctrl->tpfnames.map[m],   ppctrl->tpfnames.map[i]);
          strcpy(ppctrl->crdnames.map[m],   ppctrl->crdnames.map[i]);
          strcpy(ppctrl->outbases.map[m],   ppctrl->outbases.map[i]);
          strcpy(ppctrl->trjbases.map[m],   ppctrl->trjbases.map[i]);
          strcpy(ppctrl->rstrtbases.map[m], ppctrl->rstrtbases.map[i]);
        }
        sprintf(&ppctrl->outbases.map[m][strlen(ppctrl->outbases.map[m])],
                "_R%d", j + 1);
        sprintf(&ppctrl->trjbases.map[m][strlen(ppctrl->trjbases.map[m])],
                "_R%d", j + 1);
        sprintf(&ppctrl->rstrtbases.map[m][strlen(ppctrl->rstrtbases.map[m])],
                "_R%d", j + 1);
        m--;
      }
      continue;
    }

    // Make temperature and Hamiltonian replicas
    for (j = ppctrl->Treplicas[i] - 1; j >= 0; j--) {
      for (k = ppctrl->Preplicas[i] - 1; k >= 0; k--) {
        if (m != i) {
          strcpy(ppctrl->tpinames.map[m],   ppctrl->tpinames.map[i]);
          strcpy(ppctrl->tpfnames.map[m],   ppctrl->tpfnames.map[i]);
          strcpy(ppctrl->crdnames.map[m],   ppctrl->crdnames.map[i]);
          strcpy(ppctrl->outbases.map[m],   ppctrl->outbases.map[i]);
          strcpy(ppctrl->trjbases.map[m],   ppctrl->trjbases.map[i]);
          strcpy(ppctrl->rstrtbases.map[m], ppctrl->rstrtbases.map[i]);
        }
        if (ppctrl->Treplicas[i] > 1 && ppctrl->Preplicas[i] > 1) {
          sprintf(&ppctrl->outbases.map[m][strlen(ppctrl->outbases.map[m])],
                  "_T%d_P%d", j + 1, k + 1);
          sprintf(&ppctrl->trjbases.map[m][strlen(ppctrl->trjbases.map[m])],
                  "_T%d_P%d", j + 1, k + 1);
          slen = strlen(ppctrl->rstrtbases.map[m]);
          sprintf(&ppctrl->rstrtbases.map[m][slen], "_T%d_P%d", j + 1, k + 1);
        }
        else if (ppctrl->Treplicas[i] > 1) {
          sprintf(&ppctrl->outbases.map[m][strlen(ppctrl->outbases.map[m])],
                  "_T%d", j + 1);
          sprintf(&ppctrl->trjbases.map[m][strlen(ppctrl->trjbases.map[m])],
                  "_T%d", j + 1);
          slen = strlen(ppctrl->rstrtbases.map[m]);
          sprintf(&ppctrl->rstrtbases.map[m][slen], "_T%d", j + 1);
        }
        else if (ppctrl->Preplicas[i] > 1) {
          sprintf(&ppctrl->outbases.map[m][strlen(ppctrl->outbases.map[m])],
                  "_P%d", k + 1);
          sprintf(&ppctrl->trjbases.map[m][strlen(ppctrl->trjbases.map[m])],
                  "_P%d", k + 1);
          slen = strlen(ppctrl->rstrtbases.map[m]);
          sprintf(&ppctrl->rstrtbases.map[m][slen], "_P%d", k + 1);
        }
        if (ppctrl->Treplicas[i] == 1) {
          ppctrl->Tmix[m] = 0.0;
          if (ppctrl->Tranges[i].x >= 0.0) {
            ppctrl->T[m] = ppctrl->Tranges[i].x;
          }
          else {
            ppctrl->T[m] = tj->Ttarget;
          }
        }
        else {
          ppctrl->Tmix[m] = (double)(j) / (double)(ppctrl->Treplicas[i] - 1);
          ppctrl->T[m] = (ppctrl->Tmix[m] * ppctrl->Tranges[i].y) +
                         ((1.0 - ppctrl->Tmix[m]) * ppctrl->Tranges[i].x);
        }
        if (ppctrl->Preplicas[i] == 1) {
          ppctrl->Pmix[m] = 0.0;
        }
        else {
          ppctrl->Pmix[m] = (double)(k) / (double)(ppctrl->Preplicas[i] - 1);
        }
        m--;
      }
    }
  }
  ppctrl->nsys = maxsys;

  // Allocate CDF file handles
  if (tj->ioutfm == 1) {
    ppctrl->CDFHandles = (cdftrj*)malloc(ppctrl->nsys * sizeof(cdftrj));
  }

  // Make determinations based on &cntrl namelist input
  if (tj->ntwx == 0) {
    ppctrl->notraj = 1;
    if ((double)tj->nstep > 2.0e9) {
      printf("ExpandSystemIO >> Error.  No trajectory output was requested, "
	     "which is\nExpandSystemIO >> acceptable so long as the total "
	     "number of steps is less\nExpandSystemIO >> than 2^31 (~2 "
	     "billion).  The requested number of steps,\nExpandSystemIO >> "
	     "%lld, cannot be accommodated.\n", tj->nstep);
      exit(1);
    }
    tj->ntwx = tj->nstep;
  }
  else {
    ppctrl->notraj = 0;
  }
  if (ppctrl->rattle == 1 && ppctrl->nBondSteps > 1) {
    printf("ExpandSystemIO >> Warning.  Multiple time steps are incompatible "
	   "with RATTLE.\nExpandSystemIO >> Setting the number of bonded "
	   "interaction minor steps to 1.\n");
    ppctrl->nBondSteps = 1;
  }
  if (ppctrl->nBondSteps <= 0) {
    if (ppctrl->rattle == 1) {
      ppctrl->nBondSteps = 1;
    }
    else {
      ppctrl->nBondSteps = (int)((tj->dt + 9.999999e-4) / 0.001);
    }
    if (ppctrl->nBondSteps < 1) {
      ppctrl->nBondSteps = 1;
    }
  }
  if (ppctrl->rattle == 0 && tj->dt / (double)ppctrl->nBondSteps >= 0.001001) {
    printf("ExpandSystemIO >> Warning.  The number of bond minor steps (%d) "
           "will couple\nExpandSystemIO >> with the base time step (%.4lfps) "
           "to give a time step of\nExpandSystemIO >> %.4lfps for "
           "updating bonded interactions.  This is larger\nExpandSystemIO >> "
           "than is typically recommended.\n", ppctrl->nBondSteps, tj->dt,
           tj->dt / (double)ppctrl->nBondSteps);
  }
  else if (ppctrl->rattle == 1 && tj->dt > 0.004) {
    printf("ExpandSystemIO >> Warning.  Even with RATTLE (SHAKE equivalent "
	   "for mdgx's\nExpandSystemIO >> velocity Verlet integrator) "
	   "constraining bonds to Hydrogen,\nExpandSystemIO >> the maximum "
	   "recommended time step is 4fs.  With time steps\nExpandSystemIO >> "
	   "Longer than 2fs, hydrogen mass repartitioning is essential\n"
	   "ExpandSystemIO >> for energy conservation.\n");
  }
  if (tj->ntpr != 0) {
    if (tj->ntwx % tj->ntpr != 0) {
      tj->ntwx = ((tj->ntwx / tj->ntpr) + 1) * tj->ntpr;
      printf("ExpandSystemIO >> Warning.  The number diagnostic steps (mdout "
             "writing)\nExpandSystemIO >> must be a multiple of the number of "
             "trajectory writing\nExpandSystemIO >> steps.  Setting ntwx to "
             "%d.\n", tj->ntwx);
    }
  }
  if (tj->nstep % tj->ntwx != 0) {
    tj->nstep = ((tj->nstep / tj->ntwx) + 1) * tj->ntwx;
    printf("ExpandSystemIO >> Warning.  The total number of simulation steps "
           "must be a\nExpandSystemIO >> multiple of ntwx.  Setting nstlim / "
           "StepCount to %lld.\n", tj->nstep);
  }
  if (tj->ntwr != 0 && tj->ntwr % tj->ntwx != 0) {
    tj->ntwr = ((tj->ntwr / tj->ntwx) + 1) * tj->ntwx;
    printf("ExpandSystemIO >> Warning.  Restart files must be written at "
	   "multiples of the\nExpandSystemIO >> trajctory output frequency.  "
	   "Setting ntwr to %d.\n", tj->ntwr);
  }
  i = 0;
  while (tj->ntwr != 0 && tj->nstep % tj->ntwr != 0) {
    tj->ntwr += tj->ntwx;
    i = 1;
  }
  if (i == 1) {
    printf("ExpandSystemIO >> Warning.  The restart file writing frequency "
	   "must be a\nExpandSystemIO >> factor of the total number of "
	   "steps.  Setting ntwr to\nExpandSystemIO >> %d.\n", tj->ntwr);
  }
  ppctrl->nsegment = tj->nstep / tj->ntwx;
  ppctrl->nsgmdout = tj->ntwx  / tj->ntpr;

  // Free allocated memory
  free(toplen);
  free(duptops);
}

//-----------------------------------------------------------------------------
// InitIDArray: fill an ID array with blank information that will point to
//              atoms which exist and are not likely to generate severe
//              problems.
//
// Arguments:
//   idx:     output data array of indices to be initialized
//   nbexcl:  non-bonded exclusion matrix for the system
//   np:      the length of idx
//   order:   the order of indices to make (2 = bonds or pairs, 3 = angles,
//            4 = dihedrals)
//-----------------------------------------------------------------------------
static void InitIDArray(void* idx, dmat *nbexcl, int np, int order)
{
  int i, j, atmA, atmB, atmC, atmD;
  unsigned int *uidx;
  uint2 *u2idx;

  // Type-cast the output array
  if (order == 2 || order == 3) {
    uidx = (unsigned int*)idx;
  }
  else if (order == 4) {
    u2idx = (uint2*)idx;
  }

  // Loop over all items and find four unique atoms, preferably atoms that
  // all have valid, non-excluded non-bonded interactions with one another
  atmA = 0;
  for (i = 0; i < np; i++) {
    atmB = -1;
    j = atmA + 1;
    while (atmB == -1 && j < nbexcl->row) {
      if (nbexcl->map[atmA][j] > 1.0e-8) {
        atmB = j;
      }
      j++;
    }
    j = 0;
    while (atmB == -1 && j < atmA) {
      if (nbexcl->map[atmA][j] > 1.0e-8) {
        atmB = j;
      }
      j++;
    }
    j = 0;
    while (atmB == -1 && j < nbexcl->row) {
      if (j != atmA) {
        atmB = j;
      }
      j++;
    }
    if (order >= 3) {
      atmC = -1;
      j = atmB + 1;
      while (atmC == -1 && j < nbexcl->row) {
        if (nbexcl->map[atmA][j] > 1.0e-8 && j != atmA) {
          atmC = j;
        }
        j++;
      }
      j = 0;
      while (atmC == -1 && j < atmB) {
        if (nbexcl->map[atmA][j] > 1.0e-8 && j != atmA) {
          atmC = j;
        }
        j++;
      }
      j = 0;
      while (atmC == -1 && j < nbexcl->row) {
        if (j != atmA) {
          atmC = j;
        }
        j++;
      }
      if (order >= 4) {
        atmD = -1;
        j = atmC + 1;
        while (atmD == -1 && j < nbexcl->row) {
          if (nbexcl->map[atmA][j] > 1.0e-8 && j != atmA && j != atmB) {
            atmD = j;
          }
          j++;
        }
        j = 0;
        while (atmD == -1 && j < atmC) {
          if (nbexcl->map[atmA][j] > 1.0e-8 && j != atmA && j != atmB) {
            atmD = j;
          }
          j++;
        }
        j = 0;
        while (atmD == -1 && j < nbexcl->row) {
          if (j != atmA && j != atmB) {
            atmD = j;
          }
          j++;
        }
      }
    }
    if (order == 2) {
      uidx[i] = ((atmA << 16) | atmB);
    }
    else if (order == 3) {
      uidx[i] = ((atmA << 20) | (atmB << 10) | atmC);
    }
    else {
      u2idx[i].x = ((atmA << 20) | (atmB << 10) | atmC);
      u2idx[i].y = ((atmD << 16) | 0x1);
    }
    atmA++;
    if (atmA == nbexcl->row) {
      atmA = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// SetupRattleUnits: this complex procedure, similar in spirit to the bond
//                   work units in pmemd.cuda, gets its own routine.  It is
//                   called by InitGmsTopologies() below to set up data arrays
//                   for rapid execution of RATTLE constraints on the GPU.
//-----------------------------------------------------------------------------
static void SetupRattleUnits(pepcon *ppctrl, prmtop* tpbank, gpuMultiSim *gms)
{
  int h, i, j, k, warpcnst, nitem, atmIDs, cnstcount, ididx, atma, atmb;

  // Record the read start and read stop positions for each system.
  cnstcount = 0;
  gms->cnstReadLimits = CreateGpuInt2(ppctrl->nsys, 1);
  for (h = 0; h < ppctrl->nsys; h++) {

    // Log the read start position
    gms->cnstReadLimits.HostData[h].x = cnstcount;

    // Make warp instructions to fulfill all constraint groups
    warpcnst = 0;
    for (i = 0; i < tpbank[h].natom; i++) {
      if (tpbank[h].SHL[i].exe == 2) {
        nitem = tpbank[h].SHL[i].blist[0];
        if (nitem > GRID) {
          printf("SetupRattleUnits >> Error.  Atom %d controls %d "
                 "constraints.  This\nSetupRattleUnits >> exceeds the limit "
                 "permitted by the GPU code (%d).\n", i, nitem, GRID);
          exit(1);
        }
        if (warpcnst + nitem < GRID) {
          warpcnst += nitem;
        }
        else {
          cnstcount += GRID;
          warpcnst = nitem;
        }
      }
    }
    if (warpcnst > 0) {
      cnstcount += GRID;
    }

    // Begin building constraint instruction sets.  Each instruction set
    // is a series of up to 32 constraints, containing whole groups as
    // identified above.  If there are already 29 constraints in a set
    // and the next group contains 4 constraints, it has to go into the
    // next set.

    // Log the read stop position
    gms->cnstReadLimits.HostData[h].y = cnstcount;
  }

  // Build warp constraint instructions.  As with other things,
  // this will execute on the CPU in a way that emulates the GPU.
  gms->cnstInfo = CreateGpuUInt2(cnstcount, 1);
  cnstcount = 0;
  for (h = 0; h < ppctrl->nsys; h++) {
    warpcnst = 0;
    for (i = 0; i < tpbank[h].natom; i++) {
      if (tpbank[h].SHL[i].exe == 2) {
        nitem = tpbank[h].SHL[i].blist[0];
        if (warpcnst + nitem < GRID) {
          for (j = 0; j < nitem; j++) {
            ididx = (3 * tpbank[h].SHL[i].blist[0]) + 2;
            atma = tpbank[h].SHL[i].blist[3*j + 1];
            atmb = tpbank[h].SHL[i].blist[3*j + 2];
            atma = tpbank[h].SHL[i].blist[ididx + atma];
            atmb = tpbank[h].SHL[i].blist[ididx + atmb];
            atmIDs = ((atma << 16) | atmb);
            gms->cnstInfo.HostData[cnstcount].x = atmIDs;
            gms->cnstInfo.HostData[cnstcount].y =
              tpbank[h].SHL[i].blist[3*j + 3];
            cnstcount++;
          }
          warpcnst += nitem;
        }
        else {
          k = ((cnstcount + GRID_BITS_MASK) / GRID) * GRID;
          for (j = cnstcount; j < k; j++) {
            gms->cnstInfo.HostData[j].x = 0xffffffff;
            gms->cnstInfo.HostData[j].y = 0;
          }
          cnstcount = k;
          for (j = 0; j < nitem; j++) {
            ididx = (3 * tpbank[h].SHL[i].blist[0]) + 2;
            atma = tpbank[h].SHL[i].blist[3*j + 1];
            atmb = tpbank[h].SHL[i].blist[3*j + 2];
            atma = tpbank[h].SHL[i].blist[ididx + atma];
            atmb = tpbank[h].SHL[i].blist[ididx + atmb];
            atmIDs = ((atma << 16) | atmb);
            gms->cnstInfo.HostData[cnstcount].x = atmIDs;
            gms->cnstInfo.HostData[cnstcount].y =
              tpbank[h].SHL[i].blist[3*j + 3];
            cnstcount++;
          }
          warpcnst = nitem;
        }
      }
    }
    if (warpcnst > 0) {
      j = ((cnstcount + GRID_BITS_MASK) / GRID) * GRID;
      for (i = cnstcount; i < j; i++) {
        gms->cnstInfo.HostData[i].x = 0xffffffff;
        gms->cnstInfo.HostData[i].y = 0;
      }
      cnstcount = j;
    }
  }

  // Set device pointers
  gms->DVCcnstReadLimits = gms->cnstReadLimits.DevcData;
  gms->DVCcnstInfo       = gms->cnstInfo.DevcData;
}

//-----------------------------------------------------------------------------
// InitGmsConstants: the last of several routines for constructing the GPU
//                   Multi-Simulation struct.  This one sets constants within
//                   the struct so that the GPU understands the simulation
//                   conditions.
//
// Arguments:
//   ppctrl:  peptide simulations control data (file names, directives
//            for certain enhanced sampling methods)
//   tj:      trajectory control directives (many will be transferred to the
//            GPU Multi-Simulator)
//   gms:     the thing to initialize
//-----------------------------------------------------------------------------
static void InitGmsConstants(pepcon *ppctrl, trajcon *tj, gpuMultiSim *gms)
{
  int i;
  double gammai, hdt;

  // Output control information
  gms->ntwx         = tj->ntwx;
  gms->ntpr         = tj->ntpr;
  gms->nsgmdout     = ppctrl->nsgmdout;

  // The trajectory control data.  The Langevin thermostat is geared toward
  // the long timestep in the trajcon struct, not the shorter one in this
  // struct.  It will be applied only on the low-frequency steps.
  gms->nsys           = ppctrl->nsys;
  gms->dt             = tj->dt / (double)ppctrl->nBondSteps;
  gms->dtVF           = sqrt(418.4) * tj->dt / (double)ppctrl->nBondSteps;
  gms->invdtVF        = 1.0 / (sqrt(418.4) * tj->dt /
                               (double)ppctrl->nBondSteps);
  gms->igseed         = tj->igseed;
  gms->rndcon         = tj->rndcon;
  gms->slowFrcMult    = ppctrl->nBondSteps;
  gms->rattle         = ppctrl->rattle;
  gms->maxRattleIter  = tj->MaxRattleIter;
  gms->rattleTol      = tj->rattletol;
  gms->rtoldt         = tj->rattletol * (double)ppctrl->nBondSteps /
                        (sqrt(418.4) * tj->dt);
  gms->hdtm2invm      = FPSCALEfrc / (0.5 * gms->dt * sqrt(418.4));
  gms->hdtm2mass      = (0.5 * gms->dt * sqrt(418.4)) / FPSCALEfrc;

  // Set up the thermostat.  This will take into account
  // multiple time steps, if they are in effect.
  gms->Tstat.active   = (tj->ntt == 3);
  gms->Tstat.gamma_ln = tj->lnth.gamma_ln;
  gammai              = gms->Tstat.gamma_ln / sqrt(418.4);
  hdt                 = 0.5 * sqrt(418.4) * gms->dt;
  gms->Tstat.c_implic = 1.0 / (1.0 + gammai*hdt);
  gms->Tstat.c_explic = 1.0 - gammai*hdt;
  if (gms->Tstat.active == 1) {
    gms->Tstat.sdfac    = CreateGpuFloat(gms->nsys, 1);
    for (i = 0; i < gms->nsys; i++) {
      gms->Tstat.sdfac.HostData[i] = sqrt(gammai * GASCNST * ppctrl->T[i] /
                                          hdt);
    }

    // Set device pointer for the thermostat array
    gms->DVCsdfac     = gms->Tstat.sdfac.DevcData;
  }
  gms->Tstat.refresh  = ppctrl->prngRefresh;

  // Potential function control data, in particular GB parameters
  gms->igb          = ppctrl->igb;
  if (ppctrl->gboffset < -9.9e7) {
    gms->GBOffset = (float)0.09;
  }
  else {
    gms->GBOffset = ppctrl->gboffset;
  }
  if (gms->igb == 7) {
    gms->GBNeckScale = (float)0.361825;
  }
  else if (gms->igb == 8) {
    gms->GBNeckScale = (float)0.826836;
    if (ppctrl->gboffset < -9.9e7) {
      gms->GBOffset = (float)0.195141;
    }
  }
  gms->GBNeckCut    = GBNECKCUT;
  gms->dielectric   = ppctrl->dielectric;
  gms->kappa        = (float)0.0;
  gms->kscale       = KSCALE;
}

//-----------------------------------------------------------------------------
// InitGmsTopologies: one of several routines that constructs the GPU Multi-
//                    Simulation struct.  This one loads up topological
//                    information for all systems.
//
// Arguments:
//   ppctrl:  peptide simulations control data (file names, directives
//            for certain enhanced sampling methods)
//   tpbank:  list of all system topologies
//   gms:     carries instructions for computing forces on all systems, and
//            propagating dynamics (ported to the GPU if available)
//-----------------------------------------------------------------------------
static void InitGmsTopologies(pepcon *ppctrl, prmtop* tpbank, gpuMultiSim *gms)
{
  int h, i, j, k, m, ntile, katm, matm, nbcon, nattn, iatom, znumber;
  int atomcount, typecount, bondcount, anglcount, dihecount;
  unsigned int ubuff;
  uint2 u2buff;
  float2 f2buff;
  float4 f4buff;
  float val, val2;
  char tlet[HALFGRID];
  dmat* nbexcl;

  // Load Lennard-Jones interaction matrices
  // and all atoms' non-bonded descriptors
  nbcon = 0;
  atomcount = 0;
  gms->atomReadLimits = CreateGpuInt2(ppctrl->nsys, 1);
  for (i = 0; i < ppctrl->nsys; i++) {
    if (tpbank[i].ntypes > nbcon) {
      nbcon = tpbank[i].ntypes;
    }
    gms->atomReadLimits.HostData[i].x = atomcount;
    gms->atomReadLimits.HostData[i].y = atomcount + tpbank[i].natom;
    atomcount += ((tpbank[i].natom + GRID_BITS_MASK) / GRID) * GRID;
  }
  nbcon *= nbcon;
  gms->ljABoffset = nbcon;
  gms->atomLJID   = CreateGpuInt(atomcount, 1);
  gms->atomQ      = CreateGpuFloat(atomcount, 1);
  gms->atomMass   = CreateGpuFloat(atomcount, 1);
  gms->atomHDTM   = CreateGpuFloat(atomcount, 1);
  gms->ljFtab     = CreateGpuFloat2(nbcon * ppctrl->nsys, 1);
  gms->ljUtab     = CreateGpuFloat2(nbcon * ppctrl->nsys, 1);
  gms->typeCounts = CreateGpuInt(ppctrl->nsys, 1);
  atomcount = 0;
  for (i = 0; i < ppctrl->nsys; i++) {
    typecount = tpbank[i].ntypes;
    gms->typeCounts.HostData[i] = typecount;
    for (j = 0; j < typecount; j++) {
      for (k = 0; k < typecount; k++) {
        m = nbcon*i + j*typecount + k;
        gms->ljFtab.HostData[m].x = tpbank[i].LJftab.map[j][2*k];
        gms->ljFtab.HostData[m].y = tpbank[i].LJftab.map[j][2*k + 1];
        gms->ljUtab.HostData[m].x = tpbank[i].LJutab.map[j][2*k];
        gms->ljUtab.HostData[m].y = tpbank[i].LJutab.map[j][2*k + 1];
      }
    }
    for (j = 0; j < tpbank[i].natom; j++) {
      gms->atomQ.HostData[atomcount]      = sqrt(BIOQ) * tpbank[i].Charges[j];
      gms->atomMass.HostData[atomcount]   = tpbank[i].Masses[j];
      gms->atomHDTM.HostData[atomcount]   = 0.5 * gms->dt * sqrt(418.4) /
                                            (FPSCALEfrc * tpbank[i].Masses[j]);
      gms->atomLJID.HostData[atomcount]   = abs(tpbank[i].LJIdx[j]);
      atomcount++;
    }
    atomcount = ((atomcount + GRID_BITS_MASK) / GRID) * GRID;
  }

  // Log the number of system degrees of freedom in a convenient product
  gms->invNDF = CreateGpuDouble(ppctrl->nsys, 1);
  for (i = 0; i < ppctrl->nsys; i++) {
    gms->invNDF.HostData[i] = 1.0 / ((double)tpbank[i].ndf * GASCNST);
  }

  // Construct full non-bonded interaction matrices
  nbexcl = (dmat*)malloc(ppctrl->nsys * sizeof(dmat));
  gms->nbReadLimits = CreateGpuInt4(ppctrl->nsys, 1);
  for (i = 0; i < ppctrl->nsys; i++) {
    nbexcl[i] = CompExclMatrix(&tpbank[i], NULL);
  }

  // Allocate non-bonded pair masks
  j = 0;
  for (i = 0; i < ppctrl->nsys; i++) {
    ntile = (tpbank[i].natom + HALFGRID_MASK) / HALFGRID;
    gms->nbReadLimits.HostData[i].x = j;
    j += GRID * ntile * (ntile + 1) / 2;
    gms->nbReadLimits.HostData[i].y = j;
  }
#ifdef AMBER_PLATFORM_AMD
  gms->pairBitMasks = CreateGpuUInt2(j, 1);
#else
  gms->pairBitMasks = CreateGpuUInt(j, 1);
#endif

  // Quickly count the number of nonzero exclusions (attenuations)
  // so that a lit of them for all systems can be allocated
  nattn = 0;
  for (h = 0; h < ppctrl->nsys; h++) {
    gms->nbReadLimits.HostData[h].z = nattn;
    for (i = 1; i < tpbank[h].natom; i++) {
      for (j = 0; j < i; j++) {
        val = nbexcl[h].map[i][j];
        if (val > 1.0e-8 && val < 0.99999999) {
          nattn++;
        }
      }
    }
    gms->nbReadLimits.HostData[h].w = nattn;
    nattn = ((nattn + GRID_BITS_MASK) / GRID) * GRID;
  }
  gms->attnIDs = CreateGpuUInt(nattn, 1);
  gms->attnFactors = CreateGpuFloat2(nattn, 1);
  for (h = 0; h < ppctrl->nsys; h++) {
    i = gms->nbReadLimits.HostData[h].z;
    j = gms->nbReadLimits.HostData[h].w;
    InitIDArray(&gms->attnIDs.HostData[i], &nbexcl[h], j - i, 2);
  }

  // Loop over all 16 x 16 atom tiles, make a set of bit masks
  // describing the 16 exclusions of each of the 16 atoms, store
  // that set in sequence.
  nattn = 0;
  for (h = 0; h < ppctrl->nsys; h++) {
    ntile = (tpbank[h].natom + HALFGRID_MASK) / HALFGRID;

    // Nested loop over tiles.  This will set the ordering of the non-bonded
    // pair list on the GPU!
    nbcon = 0;
    for (i = 0; i < ntile; i++) {
      for (j = 0; j <= i; j++) {
        for (k = 0; k < HALFGRID; k++) {
#ifdef AMBER_PLATFORM_AMD
          u2buff.x = ((i << 24) | (j << 16));
          u2buff.y = 0;
#else
          ubuff = ((i << 24) | (j << 16));
#endif
          katm = i*HALFGRID + k;
          if (katm < tpbank[h].natom) {
            for (m = 0; m < HALFGRID; m++) {
              matm = j*HALFGRID + m;
              if (matm < katm) {
                val  = nbexcl[h].map[katm][matm];
                val2 = nbexcl[h].map[matm][katm];
#ifdef AMBER_PLATFORM_AMD
                u2buff.y |= ((val >= 0.99999999 && val2 >= 0.99999999) << m);
#else
                ubuff |= ((val >= 0.99999999 && val2 >= 0.99999999) << m);
#endif
                if ((val  > 1.0e-8 && val  < 0.99999999) ||
                    (val2 > 1.0e-8 && val2 < 0.99999999)) {

                  // Dihedrals will be calculated on the GPU much as they
                  // are in CPU standard MD: with attenuations subtracted
                  // from the full interactions post-facto.
                  gms->attnIDs.HostData[nattn] = ((katm << 16) | matm);
                  gms->attnFactors.HostData[nattn].x = val2;
                  gms->attnFactors.HostData[nattn].y = val;
                  nattn++;
                }
              }
            }
          }
          m = gms->nbReadLimits.HostData[h].x + nbcon;
#ifdef AMBER_PLATFORM_AMD
          gms->pairBitMasks.HostData[m] = u2buff;
#else
          gms->pairBitMasks.HostData[m] = ubuff;
#endif
          nbcon++;
        }

        // Second pass to make the m atom-perspective exclusion masks
        for (m = 0; m < HALFGRID; m++) {
        #ifdef AMBER_PLATFORM_AMD
          u2buff.x = ((i << 24) | (j << 16));
          u2buff.y = 0;
        #else
          ubuff = ((i << 24) | (j << 16));
        #endif
          matm = j*HALFGRID + m;
          if (matm < tpbank[h].natom) {
            for (k = 0; k < HALFGRID; k++) {
              katm = i*HALFGRID + k;
              if (katm > matm && katm < tpbank[h].natom) {
                val  = nbexcl[h].map[katm][matm];
                val2 = nbexcl[h].map[matm][katm];
#ifdef AMBER_PLATFORM_AMD
                u2buff.y |= ((val >= 0.99999999 && val2 >= 0.99999999) << k);
#else
                ubuff |= ((val >= 0.99999999 && val2 >= 0.99999999) << k);
#endif
              }
            }
          }
          k = gms->nbReadLimits.HostData[h].x + nbcon;
#ifdef AMBER_PLATFORM_AMD
          gms->pairBitMasks.HostData[k] = u2buff;
#else
          gms->pairBitMasks.HostData[k] = ubuff;
#endif
          nbcon++;
        }
      }
    }
    nattn = ((nattn + GRID_BITS_MASK) / GRID) * GRID;
  }

  // Allocate and fill GB neck arrays
  if (gms->igb != 6) {
    gms->NeckID    = CreateGpuInt(atomcount, 1);
    gms->rborn     = CreateGpuFloat(atomcount, 1);
    gms->reff      = CreateGpuFloat(atomcount, 1);
    gms->GBalpha   = CreateGpuFloat(atomcount, 1);
    gms->GBbeta    = CreateGpuFloat(atomcount, 1);
    gms->GBgamma   = CreateGpuFloat(atomcount, 1);
    gms->Fscreen   = CreateGpuFloat(atomcount, 1);
    for (h = 0; h < ppctrl->nsys; h++) {
      for (i = gms->atomReadLimits.HostData[h].x;
           i < gms->atomReadLimits.HostData[h].y; i++) {
        iatom = i - gms->atomReadLimits.HostData[h].x;
        gms->rborn.HostData[i] = tpbank[h].Radii[iatom];
        if (gms->igb != 7 && gms->igb != 8) {
          gms->Fscreen.HostData[i] = tpbank[h].Screen[iatom];
        }
        if (gms->igb == 2) {
          gms->GBalpha.HostData[i] = (float)0.8;
          gms->GBbeta.HostData[i]  = (float)0.0;
          gms->GBgamma.HostData[i] = (float)2.909125;
        }
        else if (gms->igb == 5) {
          gms->GBalpha.HostData[i] = (float)1.0;
          gms->GBbeta.HostData[i]  = (float)0.8;
          gms->GBgamma.HostData[i] = (float)4.85;
        }
        else if (gms->igb == 7) {
          gms->GBalpha.HostData[i] = (float)1.09511284;
          gms->GBbeta.HostData[i]  = (float)1.90792938;
          gms->GBgamma.HostData[i] = (float)2.50798245;
          znumber = AtomCode(gms->atomMass.HostData[i],
                             tpbank[h].ZNumber[iatom], gms->atomQ.HostData[i],
                             tlet) + 1.0e-8;
          switch (znumber) {
          case 1:  // Hydrogen
            gms->Fscreen.HostData[i] = (float)1.09085413633e0;
            break;
          case 6:  // Carbon
            gms->Fscreen.HostData[i] = (float)4.84353823306e-1;
            break;
          case 7:  // Nitrogen
            gms->Fscreen.HostData[i] = (float)7.00147318409e-1;
            break;
          case 8:  // Oxygen
            gms->Fscreen.HostData[i] = (float)1.06557401132e0;
            break;
          case 16: // Sulfur
            gms->Fscreen.HostData[i] = (float)6.02256336067e-1;
            break;
          default: // Not optimized
            gms->Fscreen.HostData[i] = (float)5.0e-1;
          }
        }
        else if (gms->igb == 8) {
          znumber = AtomCode(gms->atomMass.HostData[i],
                             tpbank[h].ZNumber[iatom], gms->atomQ.HostData[i],
                             tlet) + 1.0e-8;
          switch (znumber) {
          case 1: // Hydrogen
            gms->GBalpha.HostData[i] = GB8_ALPHA_H;
            gms->GBbeta.HostData[i]  = GB8_BETA_H;
            gms->GBgamma.HostData[i] = GB8_GAMMA_H;
            gms->Fscreen.HostData[i] = GB8_SCREEN_H;
            break;
          case 6: // Carbon
            gms->GBalpha.HostData[i] = GB8_ALPHA_C;
            gms->GBbeta.HostData[i]  = GB8_BETA_C;
            gms->GBgamma.HostData[i] = GB8_GAMMA_C;
            gms->Fscreen.HostData[i] = GB8_SCREEN_C;
            break;
          case 7: // Nitrogen
            gms->GBalpha.HostData[i] = GB8_ALPHA_N;
            gms->GBbeta.HostData[i]  = GB8_BETA_N;
            gms->GBgamma.HostData[i] = GB8_GAMMA_N;
            gms->Fscreen.HostData[i] = GB8_SCREEN_N;
            break;
          case 8: // Oxygen
            gms->GBalpha.HostData[i] = GB8_ALPHA_OS;
            gms->GBbeta.HostData[i]  = GB8_BETA_OS;
            gms->GBgamma.HostData[i] = GB8_GAMMA_OS;
            gms->Fscreen.HostData[i] = GB8_SCREEN_O;
            break;
          case 15: // Phosphorus
            gms->GBalpha.HostData[i] = GB8_ALPHA_P;
            gms->GBbeta.HostData[i]  = GB8_BETA_P;
            gms->GBgamma.HostData[i] = GB8_GAMMA_P;
            gms->Fscreen.HostData[i] = GB8_SCREEN_P;
            break;
          case 16: // Sulfur
            gms->GBalpha.HostData[i] = GB8_ALPHA_OS;
            gms->GBbeta.HostData[i]  = GB8_BETA_OS;
            gms->GBgamma.HostData[i] = GB8_GAMMA_OS;
            gms->Fscreen.HostData[i] = GB8_SCREEN_S;
            break;
          default: // Non-optimized defaults
            gms->GBalpha.HostData[i] = 1.0;
            gms->GBbeta.HostData[i]  = 0.8;
            gms->GBgamma.HostData[i] = 4.85;
            gms->Fscreen.HostData[i] = 5.0e-1;
          }
        }
        if (gms->igb == 7 || gms->igb == 8) {
          gms->NeckID.HostData[i] =
            (int)(((gms->rborn.HostData[i] - (float)1.0) * (float)20.0) +
                  (float)0.5);
          if (gms->NeckID.HostData[i] < 0 || gms->NeckID.HostData[i] > 20) {
            printf("InitGmsTopologies >> Error.  Atom %d in system %d is "
                   "outside the allowed\nInitGmsTopologies >> range of 1-2 "
                   "Angstroms for igb = %d.  Regenerate the\n"
                   "InitGmsTopologies >> topology with 'bondi' radii.", iatom,
                   h, gms->igb);
            exit(1);
          }
        }
      }
    }

    // Fill up the Neck GB parameters array
    gms->neckFactors = CreateGpuFloat2(24*24, 1);
    for (i = 0; i < 21; i++) {
      for (j = 0; j < 21; j++) {
        gms->neckFactors.HostData[21*i + j].x = neckMaxPos[i][j];
        gms->neckFactors.HostData[21*i + j].y = neckMaxVal[i][j];
      }
    }
  }

  // Allocate for bonded terms
  bondcount = 0;
  anglcount = 0;
  dihecount = 0;
  gms->bondReadLimits = CreateGpuInt2(ppctrl->nsys, 1);
  gms->anglReadLimits = CreateGpuInt2(ppctrl->nsys, 1);
  gms->diheReadLimits = CreateGpuInt2(ppctrl->nsys, 1);
  for (h = 0; h < ppctrl->nsys; h++) {
    gms->bondReadLimits.HostData[h].x = bondcount;
    gms->anglReadLimits.HostData[h].x = anglcount;
    gms->diheReadLimits.HostData[h].x = dihecount;
    bondcount += tpbank[h].withH.nbond + tpbank[h].woH.nbond;
    anglcount += tpbank[h].withH.nangl + tpbank[h].woH.nangl;
    for (i = 0; i < tpbank[h].natom; i++) {
      for (j = 0; j < tpbank[h].HLC[i].ndihe; j++) {
        dihecount += tpbank[h].HLC[i].HC[j].nt;
      }
    }
    gms->bondReadLimits.HostData[h].y = bondcount;
    gms->anglReadLimits.HostData[h].y = anglcount;
    gms->diheReadLimits.HostData[h].y = dihecount;

    // Pad the counts so that the last warp working on each
    // system can read blank entries
    bondcount = ((bondcount + GRID_BITS_MASK) / GRID) * GRID;
    anglcount = ((anglcount + GRID_BITS_MASK) / GRID) * GRID;
    dihecount = ((dihecount + GRID_BITS_MASK) / GRID) * GRID;
  }
  gms->bondIDs    = CreateGpuUInt(bondcount, 1);
  gms->anglIDs    = CreateGpuUInt(anglcount, 1);
  gms->diheIDs    = CreateGpuUInt2(dihecount, 1);
  gms->bondBasics = CreateGpuFloat2(bondcount, 1);
  gms->bondAugs   = CreateGpuFloat4(bondcount, 1);
  gms->anglInfo   = CreateGpuFloat2(anglcount, 1);
  gms->diheInfo   = CreateGpuFloat2(dihecount, 1);
  for (h = 0; h < ppctrl->nsys; h++) {
    i = gms->bondReadLimits.HostData[h].x;
    j = gms->bondReadLimits.HostData[h].y;
    j = ((j + GRID_BITS_MASK) / GRID) * GRID;
    InitIDArray(&gms->bondIDs.HostData[i], &nbexcl[h], j - i, 2);
    i = gms->anglReadLimits.HostData[h].x;
    j = gms->anglReadLimits.HostData[h].y;
    j = ((j + GRID_BITS_MASK) / GRID) * GRID;
    InitIDArray(&gms->anglIDs.HostData[i], &nbexcl[h], j - i, 3);
    i = gms->diheReadLimits.HostData[h].x;
    j = gms->diheReadLimits.HostData[h].y;
    j = ((j + GRID_BITS_MASK) / GRID) * GRID;
    InitIDArray(&gms->diheIDs.HostData[i], &nbexcl[h], j - i, 4);
  }

  // Fill up the bonded terms
  bondcount = 0;
  anglcount = 0;
  dihecount = 0;
  for (h = 0; h < ppctrl->nsys; h++) {
    for (i = 0; i < tpbank[h].natom; i++) {
      for (j = 0; j < tpbank[h].BLC[i].nbond; j++) {
        ubuff = ((i << 16) | tpbank[h].BLC[i].BC[j].b);
        gms->bondIDs.HostData[bondcount] = ubuff;
        k = tpbank[h].BLC[i].BC[j].t;
        f2buff.x = tpbank[h].BParam[k].K;
        f2buff.y = tpbank[h].BParam[k].l0;
        f4buff.x = tpbank[h].BParam[k].Kpull;
        f4buff.y = tpbank[h].BParam[k].lpull0;
        f4buff.z = tpbank[h].BParam[k].Kpress;
        f4buff.w = tpbank[h].BParam[k].lpress0;
        gms->bondBasics.HostData[bondcount] = f2buff;
        gms->bondAugs.HostData[bondcount] = f4buff;
        bondcount++;
      }
      for (j = 0; j < tpbank[h].ALC[i].nangl; j++) {
        ubuff = ((tpbank[h].ALC[i].AC[j].a << 20) | (i << 10) |
                 tpbank[h].ALC[i].AC[j].c);
        gms->anglIDs.HostData[anglcount] = ubuff;
        k = tpbank[h].ALC[i].AC[j].t;
        f2buff.x = tpbank[h].AParam[k].K;
        f2buff.y = tpbank[h].AParam[k].th0;
        gms->anglInfo.HostData[anglcount] = f2buff;
        anglcount++;
      }
      for (j = 0; j < tpbank[h].HLC[i].ndihe; j++) {
        u2buff.x = ((tpbank[h].HLC[i].HC[j].a << 20) | (i << 10) |
                    tpbank[h].HLC[i].HC[j].c);
        for (k = 0; k < tpbank[h].HLC[i].HC[j].nt; k++) {
          m = tpbank[h].HLC[i].HC[j].t[k];
          u2buff.y = ((tpbank[h].HLC[i].HC[j].d << 16) |
                      (int)(tpbank[h].HParam[m].N + 1.0e-8));
          gms->diheIDs.HostData[dihecount] = u2buff;
          f2buff.x = tpbank[h].HParam[m].K;
          f2buff.y = tpbank[h].HParam[m].Phi;
          gms->diheInfo.HostData[dihecount] = f2buff;
          dihecount++;
        }
      }
    }

    // Pad the counts to get to the starting point for the next system
    bondcount = ((bondcount + GRID_BITS_MASK) / GRID) * GRID;
    anglcount = ((anglcount + GRID_BITS_MASK) / GRID) * GRID;
    dihecount = ((dihecount + GRID_BITS_MASK) / GRID) * GRID;
  }

  // Set up RATTLE groups
  if (gms->rattle == 1) {
    SetupRattleUnits(ppctrl, tpbank, gms);
  }

  // Set shortcut pointers on the GPU Multi-Simulator struct
  gms->DVCatomReadLimits = gms->atomReadLimits.DevcData;
  gms->DVCbondReadLimits = gms->bondReadLimits.DevcData;
  gms->DVCbondIDs        = gms->bondIDs.DevcData;
  gms->DVCbondBasics     = gms->bondBasics.DevcData;
  gms->DVCbondAugs       = gms->bondAugs.DevcData;
  gms->DVCanglReadLimits = gms->anglReadLimits.DevcData;
  gms->DVCanglIDs        = gms->anglIDs.DevcData;
  gms->DVCanglInfo       = gms->anglInfo.DevcData;
  gms->DVCdiheReadLimits = gms->diheReadLimits.DevcData;
  gms->DVCdiheIDs        = gms->diheIDs.DevcData;
  gms->DVCdiheInfo       = gms->diheInfo.DevcData;
  gms->DVCpairBitMasks   = gms->pairBitMasks.DevcData;
  gms->DVCattnIDs        = gms->attnIDs.DevcData;
  gms->DVCattnFactors    = gms->attnFactors.DevcData;
  gms->DVCnbReadLimits   = gms->nbReadLimits.DevcData;
  gms->DVCatomQ          = gms->atomQ.DevcData;
  gms->DVCatomMass       = gms->atomMass.DevcData;
  gms->DVCatomHDTM       = gms->atomHDTM.DevcData;
  gms->DVCatomLJID       = gms->atomLJID.DevcData;
  gms->DVCtypeCounts     = gms->typeCounts.DevcData;
  gms->DVCljFtab         = gms->ljFtab.DevcData;
  gms->DVCljUtab         = gms->ljUtab.DevcData;
  if (gms->igb != 6) {
    gms->DVCNeckID         = gms->NeckID.DevcData;
    gms->DVCrborn          = gms->rborn.DevcData;
    gms->DVCreff           = gms->reff.DevcData;
    gms->DVCGBalpha        = gms->GBalpha.DevcData;
    gms->DVCGBbeta         = gms->GBbeta.DevcData;
    gms->DVCGBgamma        = gms->GBgamma.DevcData;
    gms->DVCFscreen        = gms->Fscreen.DevcData;
    gms->DVCneckFactors    = gms->neckFactors.DevcData;
  }

  // Free allocated memory
  for (i = 0; i < ppctrl->nsys; i++) {
    DestroyDmat(&nbexcl[i]);
  }
  free(nbexcl);
}

//-----------------------------------------------------------------------------
// InitGmsCoordinates: another of several routines for constructing the GPU
//                     Multi-Simulation struct.  This one loads coordinates
//                     from coord structs read from input files and puts them
//                     into arrays accessible to the GPU.  This will also
//                     allocate forces and energies.
//
// Arguments:
//   ppctrl:  peptide simulations control data (file names, directives
//            for certain enhanced sampling methods)
//   tpbank:  list of all system topologies
//   crdbank: array holding all starting coordinates read from restart files
//   gms:     carries instructions for computing forces on all systems, and
//            propagating dynamics (ported to the GPU if available)
//-----------------------------------------------------------------------------
static void InitGmsCoordinates(pepcon *ppctrl, prmtop* tpbank, coord* crdbank,
                               gpuMultiSim *gms)
{
  int i, j, k, atomcount, endidx;
  long int counter;

  // Get the total number of atoms to allocate
  atomcount = 0;
  for (i = 0; i < ppctrl->nsys; i++) {
    atomcount += ((tpbank[i].natom + GRID_BITS_MASK) / GRID) * GRID;
  }

  // Allocate arrays to hold coordinates and velocities
  gms->atomX = CreateGpuFloat(atomcount, 1);
  gms->atomY = CreateGpuFloat(atomcount, 1);
  gms->atomZ = CreateGpuFloat(atomcount, 1);
  gms->crdrX = CreateGpuFloat(atomcount, 1);
  gms->crdrY = CreateGpuFloat(atomcount, 1);
  gms->crdrZ = CreateGpuFloat(atomcount, 1);
  gms->velcX = CreateGpuFloat(atomcount, 1);
  gms->velcY = CreateGpuFloat(atomcount, 1);
  gms->velcZ = CreateGpuFloat(atomcount, 1);
  gms->velrX = CreateGpuFloat(atomcount, 1);
  gms->velrY = CreateGpuFloat(atomcount, 1);
  gms->velrZ = CreateGpuFloat(atomcount, 1);
  k = 0;
  for (i = 0; i < ppctrl->nsys; i++) {
    for (j = 0; j < tpbank[i].natom; j++) {
      gms->atomX.HostData[k] = crdbank[i].loc[3*j    ];
      gms->atomY.HostData[k] = crdbank[i].loc[3*j + 1];
      gms->atomZ.HostData[k] = crdbank[i].loc[3*j + 2];
      gms->velcX.HostData[k] = crdbank[i].vel[3*j    ];
      gms->velcY.HostData[k] = crdbank[i].vel[3*j + 1];
      gms->velcZ.HostData[k] = crdbank[i].vel[3*j + 2];
      k++;
    }
    k = ((k + GRID_BITS_MASK) / GRID) * GRID;
  }

  // Allocate arrays to hold forces
  gms->frcX  = CreateGpuInt(atomcount, 1);
  gms->frcY  = CreateGpuInt(atomcount, 1);
  gms->frcZ  = CreateGpuInt(atomcount, 1);
#ifdef CUDA
  counter = gms->rndcon;
#endif
  if (gms->Tstat.active == 1) {
    gms->langX          = CreateGpuInt(atomcount, 1);
    gms->langY          = CreateGpuInt(atomcount, 1);
    gms->langZ          = CreateGpuInt(atomcount, 1);
#ifdef CUDA
    gms->prngHeap       = CreateGpuFloat(ppctrl->totalPRNG, 1);
    gms->prngReadLimits = CreateGpuInt2(ppctrl->nsys, 1);
    j = 0;
    for (i = 0; i < ppctrl->nsys; i++) {
      gms->prngReadLimits.HostData[i].x = j;

      // Prepare the random number heap for the first
      // Velocity Verlet coordinates update
      for (k = 0; k < 3 * tpbank[i].natom; k++) {
	endidx = j + (3 * tpbank[i].natom * gms->Tstat.refresh) + k;
	gms->prngHeap.HostData[endidx] = GaussBoxMuller(&counter);
      }
      j += 3 * tpbank[i].natom * (gms->Tstat.refresh + 1);
      gms->prngReadLimits.HostData[i].y = j;
      j = ((j + GRID_BITS_MASK) / GRID) * GRID;
    }
#endif
  }
#ifdef CUDA
  gms->rndcon = counter;
#endif

  // Allocate arrays to hold energy terms
  gms->Uelec = CreateGpuDouble(ppctrl->nsys * ppctrl->nsgmdout, 1);
  gms->Uvdw  = CreateGpuDouble(ppctrl->nsys * ppctrl->nsgmdout, 1);
  gms->Usolv = CreateGpuDouble(ppctrl->nsys * ppctrl->nsgmdout, 1);
  gms->Ubond = CreateGpuDouble(ppctrl->nsys * ppctrl->nsgmdout, 1);
  gms->Uangl = CreateGpuDouble(ppctrl->nsys * ppctrl->nsgmdout, 1);
  gms->Udihe = CreateGpuDouble(ppctrl->nsys * ppctrl->nsgmdout, 1);
  gms->Ukine = CreateGpuDouble(ppctrl->nsys * ppctrl->nsgmdout, 1);
  gms->Temp  = CreateGpuDouble(ppctrl->nsys * ppctrl->nsgmdout, 1);

  // An array to hold progress trackers
  gms->systemPos = CreateGpuInt(2 * ppctrl->nsgmdout, 1);

  // Set shortcut pointers on the GPU Multi-Simulator struct
  gms->DVCatomX    = gms->atomX.DevcData;
  gms->DVCatomY    = gms->atomY.DevcData;
  gms->DVCatomZ    = gms->atomZ.DevcData;
  gms->DVCcrdrX    = gms->crdrX.DevcData;
  gms->DVCcrdrY    = gms->crdrY.DevcData;
  gms->DVCcrdrZ    = gms->crdrZ.DevcData;
  gms->DVCvelcX    = gms->velcX.DevcData;
  gms->DVCvelcY    = gms->velcY.DevcData;
  gms->DVCvelcZ    = gms->velcZ.DevcData;
  gms->DVCvelrX    = gms->velrX.DevcData;
  gms->DVCvelrY    = gms->velrY.DevcData;
  gms->DVCvelrZ    = gms->velrZ.DevcData;
  gms->DVCfrcX     = gms->frcX.DevcData;
  gms->DVCfrcY     = gms->frcY.DevcData;
  gms->DVCfrcZ     = gms->frcZ.DevcData;
  if (gms->Tstat.active == 1) {
    gms->DVClangX          = gms->langX.DevcData;
    gms->DVClangY          = gms->langY.DevcData;
    gms->DVClangZ          = gms->langZ.DevcData;
#ifdef CUDA
    gms->DVCprngHeap       = gms->prngHeap.DevcData;
    gms->DVCprngReadLimits = gms->prngReadLimits.DevcData;
#endif
  }
  gms->DVCuelec     = gms->Uelec.DevcData;
  gms->DVCuvdw      = gms->Uvdw.DevcData;
  gms->DVCusolv     = gms->Usolv.DevcData;
  gms->DVCubond     = gms->Ubond.DevcData;
  gms->DVCuangl     = gms->Uangl.DevcData;
  gms->DVCudihe     = gms->Udihe.DevcData;
  gms->DVCukine     = gms->Ukine.DevcData;
  gms->DVCtemp      = gms->Temp.DevcData;
  gms->DVCsystemPos = gms->systemPos.DevcData;
}

//-----------------------------------------------------------------------------
// UploadGmsAll: this function will upload the GPU Multi-Simulation struct,
//               including values of all pointers and contents of all data
//               arrays.  It should only need to be called once, after the
//               GPU Multi-Simulation struct has been assembled in host RAM.
//
// Arguments:
//   gms:     the GPU Multi-Simulation struct to upload
//-----------------------------------------------------------------------------
#ifdef CUDA
static void UploadGmsAll(gpuMultiSim *gms)
{
  // Takes care of putting all scalars, static arrays,
  // and pointers stored in gms up onto the GPU
  SetGmsImage(gms);

  // Upload the atom read limits
  UploadGpuInt2(&gms->atomReadLimits, NULL);

  // Upload coordinates
  UploadGpuFloat(&gms->atomX, NULL);
  UploadGpuFloat(&gms->atomY, NULL);
  UploadGpuFloat(&gms->atomZ, NULL);
  UploadGpuFloat(&gms->crdrX, NULL);
  UploadGpuFloat(&gms->crdrY, NULL);
  UploadGpuFloat(&gms->crdrZ, NULL);

  // Upload velocities
  UploadGpuFloat(&gms->velcX, NULL);
  UploadGpuFloat(&gms->velcY, NULL);
  UploadGpuFloat(&gms->velcZ, NULL);
  UploadGpuFloat(&gms->velrX, NULL);
  UploadGpuFloat(&gms->velrY, NULL);
  UploadGpuFloat(&gms->velrZ, NULL);

  // Upload forces
  UploadGpuInt(&gms->frcX, NULL);
  UploadGpuInt(&gms->frcY, NULL);
  UploadGpuInt(&gms->frcZ, NULL);
  if (gms->Tstat.active == 1) {
    UploadGpuInt(&gms->langX, NULL);
    UploadGpuInt(&gms->langY, NULL);
    UploadGpuInt(&gms->langZ, NULL);
    UploadGpuFloat(&gms->Tstat.sdfac, NULL);
#ifdef CUDA
    UploadGpuFloat(&gms->prngHeap, NULL);
    UploadGpuInt2(&gms->prngReadLimits, NULL);
#endif
  }

  // Upload energy terms
  UploadGpuDouble(&gms->Uelec, NULL);
  UploadGpuDouble(&gms->Uvdw, NULL);
  UploadGpuDouble(&gms->Usolv, NULL);
  UploadGpuDouble(&gms->Ubond, NULL);
  UploadGpuDouble(&gms->Uangl, NULL);
  UploadGpuDouble(&gms->Udihe, NULL);
  UploadGpuDouble(&gms->Ukine, NULL);
  UploadGpuDouble(&gms->Temp, NULL);

  // Upload tracking counters
  UploadGpuInt(&gms->systemPos, NULL);

  // Upload atom properties, exclusion masks, and partial attenuations
  UploadGpuFloat(&gms->atomQ, NULL);
  UploadGpuFloat(&gms->atomMass, NULL);
  UploadGpuFloat(&gms->atomHDTM, NULL);
  UploadGpuDouble(&gms->invNDF, NULL);
  UploadGpuInt(&gms->atomLJID, NULL);
#ifdef AMBER_PLATFORM_AMD
  UploadGpuUInt2(&gms->pairBitMasks, NULL);
#else
  UploadGpuUInt(&gms->pairBitMasks, NULL);
#endif
  UploadGpuInt4(&gms->nbReadLimits, NULL);
  UploadGpuUInt(&gms->attnIDs, NULL);
  UploadGpuFloat2(&gms->attnFactors, NULL);
  UploadGpuInt(&gms->typeCounts, NULL);
  UploadGpuFloat2(&gms->ljFtab, NULL);
  UploadGpuFloat2(&gms->ljUtab, NULL);

  // Upload GB attributes
  if (gms->igb != 6) {
    UploadGpuInt(&gms->NeckID, NULL);
    UploadGpuFloat(&gms->rborn, NULL);
    UploadGpuFloat(&gms->reff, NULL);
    UploadGpuFloat(&gms->GBalpha, NULL);
    UploadGpuFloat(&gms->GBbeta, NULL);
    UploadGpuFloat(&gms->GBgamma, NULL);
    UploadGpuFloat(&gms->Fscreen, NULL);
    UploadGpuFloat2(&gms->neckFactors, NULL);
  }

  // Upload bonded terms
  UploadGpuInt2(&gms->bondReadLimits, NULL);
  UploadGpuUInt(&gms->bondIDs, NULL);
  UploadGpuFloat2(&gms->bondBasics, NULL);
  UploadGpuFloat4(&gms->bondAugs, NULL);
  if (gms->rattle == 1) {
    UploadGpuInt2(&gms->cnstReadLimits, NULL);
    UploadGpuUInt2(&gms->cnstInfo, NULL);
  }
  UploadGpuInt2(&gms->anglReadLimits, NULL);
  UploadGpuUInt(&gms->anglIDs, NULL);
  UploadGpuFloat2(&gms->anglInfo, NULL);
  UploadGpuInt2(&gms->diheReadLimits, NULL);
  UploadGpuUInt2(&gms->diheIDs, NULL);
  UploadGpuFloat2(&gms->diheInfo, NULL);
}

//-----------------------------------------------------------------------------
// DownloadGmsAll: this function will download the GPU Multi-Simulation struct,
//               including values of all pointers and contents of all data
//               arrays.  It should only need to be called for debugging, after
//               the struct has been manipulated on the GPU.
//
// This is a debugging function.
//
// Arguments:
//   gms:     the GPU Multi-Simulation struct to upload
//-----------------------------------------------------------------------------
static void DownloadGmsAll(gpuMultiSim *gms)
{
  // Upload the atom read limits
  DownloadGpuInt2(&gms->atomReadLimits, NULL);

  // Download coordinates
  DownloadGpuFloat(&gms->atomX, NULL);
  DownloadGpuFloat(&gms->atomY, NULL);
  DownloadGpuFloat(&gms->atomZ, NULL);
  DownloadGpuFloat(&gms->crdrX, NULL);
  DownloadGpuFloat(&gms->crdrY, NULL);
  DownloadGpuFloat(&gms->crdrZ, NULL);

  // Download velocities
  DownloadGpuFloat(&gms->velcX, NULL);
  DownloadGpuFloat(&gms->velcY, NULL);
  DownloadGpuFloat(&gms->velcZ, NULL);
  DownloadGpuFloat(&gms->velrX, NULL);
  DownloadGpuFloat(&gms->velrY, NULL);
  DownloadGpuFloat(&gms->velrZ, NULL);

  // Download forces
  DownloadGpuInt(&gms->frcX, NULL);
  DownloadGpuInt(&gms->frcY, NULL);
  DownloadGpuInt(&gms->frcZ, NULL);
  if (gms->Tstat.active == 1) {
    DownloadGpuInt(&gms->langX, NULL);
    DownloadGpuInt(&gms->langY, NULL);
    DownloadGpuInt(&gms->langZ, NULL);
    DownloadGpuFloat(&gms->Tstat.sdfac, NULL);
#ifdef CUDA
    DownloadGpuFloat(&gms->prngHeap, NULL);
    DownloadGpuInt2(&gms->prngReadLimits, NULL);
#endif
  }

  // Download energy terms
  DownloadGpuDouble(&gms->Uelec, NULL);
  DownloadGpuDouble(&gms->Uvdw, NULL);
  DownloadGpuDouble(&gms->Usolv, NULL);
  DownloadGpuDouble(&gms->Ubond, NULL);
  DownloadGpuDouble(&gms->Uangl, NULL);
  DownloadGpuDouble(&gms->Udihe, NULL);
  DownloadGpuDouble(&gms->Ukine, NULL);
  DownloadGpuDouble(&gms->Temp, NULL);

  // Download tracking counters
  DownloadGpuInt(&gms->systemPos, NULL);

  // Download atom properties, exclusion masks, and partial attenuations
  DownloadGpuFloat(&gms->atomQ, NULL);
  DownloadGpuFloat(&gms->atomMass, NULL);
  DownloadGpuFloat(&gms->atomHDTM, NULL);
  DownloadGpuDouble(&gms->invNDF, NULL);
  DownloadGpuInt(&gms->atomLJID, NULL);
#ifdef AMBER_PLATFORM_AMD
  DownloadGpuUInt2(&gms->pairBitMasks, NULL);
#else
  DownloadGpuUInt(&gms->pairBitMasks, NULL);
#endif
  DownloadGpuInt4(&gms->nbReadLimits, NULL);
  DownloadGpuUInt(&gms->attnIDs, NULL);
  DownloadGpuFloat2(&gms->attnFactors, NULL);
  DownloadGpuInt(&gms->typeCounts, NULL);
  DownloadGpuFloat2(&gms->ljFtab, NULL);
  DownloadGpuFloat2(&gms->ljUtab, NULL);

  // Download GB attributes
  if (gms->igb != 6) {
    DownloadGpuInt(&gms->NeckID, NULL);
    DownloadGpuFloat(&gms->rborn, NULL);
    DownloadGpuFloat(&gms->reff, NULL);
    DownloadGpuFloat(&gms->GBalpha, NULL);
    DownloadGpuFloat(&gms->GBbeta, NULL);
    DownloadGpuFloat(&gms->GBgamma, NULL);
    DownloadGpuFloat(&gms->Fscreen, NULL);
    DownloadGpuFloat2(&gms->neckFactors, NULL);
  }

  // Download bonded terms
  DownloadGpuInt2(&gms->bondReadLimits, NULL);
  DownloadGpuUInt(&gms->bondIDs, NULL);
  DownloadGpuFloat2(&gms->bondBasics, NULL);
  DownloadGpuFloat4(&gms->bondAugs, NULL);
  if (gms->rattle == 1) {
    DownloadGpuInt2(&gms->cnstReadLimits, NULL);
    DownloadGpuUInt2(&gms->cnstInfo, NULL);
  }
  DownloadGpuInt2(&gms->anglReadLimits, NULL);
  DownloadGpuUInt(&gms->anglIDs, NULL);
  DownloadGpuFloat2(&gms->anglInfo, NULL);
  DownloadGpuInt2(&gms->diheReadLimits, NULL);
  DownloadGpuUInt2(&gms->diheIDs, NULL);
  DownloadGpuFloat2(&gms->diheInfo, NULL);
}
#endif

//-----------------------------------------------------------------------------
// CopyGms: copy a GPU Multi-Simulation struct, both on the host and on the
//          device.
//
// This is a debugging function.
//
// Arguments:
//   gms:      the GPU Multi-Simulation struct to copy
//   copySrc:  set to 0 to copy from host memory, 1 to copy from device memory
//-----------------------------------------------------------------------------
static gpuMultiSim CopyGms(gpuMultiSim *gms, int copySrc)
{
  int i;
  gpuMultiSim newgms;

  // Copy scalars (this will copy pointers, too, but those will be rewritten)
  newgms = *gms;

  // Copy atom reading limits
  newgms.atomReadLimits = CopyGpuInt2(&gms->atomReadLimits, copySrc);

  // Copy coordinates
  newgms.atomX = CopyGpuFloat(&gms->atomX, copySrc);
  newgms.atomY = CopyGpuFloat(&gms->atomY, copySrc);
  newgms.atomZ = CopyGpuFloat(&gms->atomZ, copySrc);
  newgms.crdrX = CopyGpuFloat(&gms->crdrX, copySrc);
  newgms.crdrY = CopyGpuFloat(&gms->crdrY, copySrc);
  newgms.crdrZ = CopyGpuFloat(&gms->crdrZ, copySrc);

  // Copy velocities
  newgms.velcX = CopyGpuFloat(&gms->velcX, copySrc);
  newgms.velcY = CopyGpuFloat(&gms->velcY, copySrc);
  newgms.velcZ = CopyGpuFloat(&gms->velcZ, copySrc);
  newgms.velrX = CopyGpuFloat(&gms->velrX, copySrc);
  newgms.velrY = CopyGpuFloat(&gms->velrY, copySrc);
  newgms.velrZ = CopyGpuFloat(&gms->velrZ, copySrc);

  // Copy Forces
  newgms.frcX  = CopyGpuInt(&gms->frcX, copySrc);
  newgms.frcY  = CopyGpuInt(&gms->frcY, copySrc);
  newgms.frcZ  = CopyGpuInt(&gms->frcZ, copySrc);
  if (gms->Tstat.active == 1) {
    newgms.langX          = CopyGpuInt(&gms->langX, copySrc);
    newgms.langY          = CopyGpuInt(&gms->langY, copySrc);
    newgms.langZ          = CopyGpuInt(&gms->langZ, copySrc);
    newgms.Tstat.sdfac    = CopyGpuFloat(&gms->Tstat.sdfac, copySrc);
#ifdef CUDA
    newgms.prngHeap       = CopyGpuFloat(&gms->prngHeap, copySrc);
    newgms.prngReadLimits = CopyGpuInt2(&gms->prngReadLimits, copySrc);
#endif
  }

  // Copy energy terms
  newgms.Uelec = CopyGpuDouble(&gms->Uelec, copySrc);
  newgms.Uvdw  = CopyGpuDouble(&gms->Uvdw, copySrc);
  newgms.Usolv = CopyGpuDouble(&gms->Usolv, copySrc);
  newgms.Ubond = CopyGpuDouble(&gms->Ubond, copySrc);
  newgms.Uangl = CopyGpuDouble(&gms->Uangl, copySrc);
  newgms.Udihe = CopyGpuDouble(&gms->Udihe, copySrc);
  newgms.Ukine = CopyGpuDouble(&gms->Ukine, copySrc);

  // Copy tracking counters
  newgms.systemPos = CopyGpuInt(&gms->systemPos, copySrc);

  // Copy atom properties, exclusion masks, and partial attenuations
  newgms.atomQ          = CopyGpuFloat(&gms->atomQ, copySrc);
  newgms.atomMass       = CopyGpuFloat(&gms->atomMass, copySrc);
  newgms.atomHDTM       = CopyGpuFloat(&gms->atomHDTM, copySrc);
  newgms.invNDF         = CopyGpuDouble(&gms->invNDF, copySrc);
  newgms.atomLJID       = CopyGpuInt(&gms->atomLJID, copySrc);
#ifdef AMBER_PLATFORM_AMD
  newgms.pairBitMasks   = CopyGpuUInt2(&gms->pairBitMasks, copySrc);
#else
  newgms.pairBitMasks   = CopyGpuUInt(&gms->pairBitMasks, copySrc);
#endif
  newgms.nbReadLimits   = CopyGpuInt4(&gms->nbReadLimits, copySrc);
  newgms.attnIDs        = CopyGpuUInt(&gms->attnIDs, copySrc);
  newgms.attnFactors    = CopyGpuFloat2(&gms->attnFactors, copySrc);
  newgms.typeCounts     = CopyGpuInt(&gms->typeCounts, copySrc);
  newgms.ljFtab         = CopyGpuFloat2(&gms->ljFtab, copySrc);
  newgms.ljUtab         = CopyGpuFloat2(&gms->ljUtab, copySrc);

  // Copy GB attributes
  if (gms->igb != 6) {
    newgms.NeckID      = CopyGpuInt(&gms->NeckID, copySrc);
    newgms.rborn       = CopyGpuFloat(&gms->rborn, copySrc);
    newgms.reff        = CopyGpuFloat(&gms->reff, copySrc);
    newgms.GBalpha     = CopyGpuFloat(&gms->GBalpha, copySrc);
    newgms.GBbeta      = CopyGpuFloat(&gms->GBbeta, copySrc);
    newgms.GBgamma     = CopyGpuFloat(&gms->GBgamma, copySrc);
    newgms.Fscreen     = CopyGpuFloat(&gms->Fscreen, copySrc);
    newgms.neckFactors = CopyGpuFloat2(&gms->neckFactors, copySrc);
  }

  // Copy bonded parameters
  newgms.bondReadLimits = CopyGpuInt2(&gms->bondReadLimits, copySrc);
  newgms.bondIDs        = CopyGpuUInt(&gms->bondIDs, copySrc);
  newgms.bondBasics     = CopyGpuFloat2(&gms->bondBasics, copySrc);
  newgms.bondAugs       = CopyGpuFloat4(&gms->bondAugs, copySrc);
  if (gms->rattle == 1) {
    newgms.cnstReadLimits = CopyGpuInt2(&gms->cnstReadLimits, copySrc);
    newgms.cnstInfo       = CopyGpuUInt2(&gms->cnstInfo, copySrc);
  }
  newgms.anglReadLimits = CopyGpuInt2(&gms->anglReadLimits, copySrc);
  newgms.anglIDs        = CopyGpuUInt(&gms->anglIDs, copySrc);
  newgms.anglInfo       = CopyGpuFloat2(&gms->anglInfo, copySrc);
  newgms.diheReadLimits = CopyGpuInt2(&gms->diheReadLimits, copySrc);
  newgms.diheIDs        = CopyGpuUInt2(&gms->diheIDs, copySrc);
  newgms.diheInfo       = CopyGpuFloat2(&gms->diheInfo, copySrc);

  // Copy pointers
  newgms.DVCatomReadLimits  = newgms.atomReadLimits.DevcData;
  newgms.DVCbondReadLimits  = newgms.bondReadLimits.DevcData;
  newgms.DVCbondIDs         = newgms.bondIDs.DevcData;
  newgms.DVCbondBasics      = newgms.bondBasics.DevcData;
  newgms.DVCbondAugs        = newgms.bondAugs.DevcData;
  if (gms->rattle == 1) {
    newgms.DVCcnstReadLimits  = newgms.cnstReadLimits.DevcData;
    newgms.DVCcnstInfo        = newgms.cnstInfo.DevcData;
  }
  newgms.DVCanglReadLimits  = newgms.anglReadLimits.DevcData;
  newgms.DVCanglIDs         = newgms.anglIDs.DevcData;
  newgms.DVCanglInfo        = newgms.anglInfo.DevcData;
  newgms.DVCdiheReadLimits  = newgms.diheReadLimits.DevcData;
  newgms.DVCdiheIDs         = newgms.diheIDs.DevcData;
  newgms.DVCdiheInfo        = newgms.diheInfo.DevcData;
  newgms.DVCpairBitMasks    = newgms.pairBitMasks.DevcData;
  newgms.DVCattnIDs         = newgms.attnIDs.DevcData;
  newgms.DVCattnFactors     = newgms.attnFactors.DevcData;
  newgms.DVCnbReadLimits    = newgms.nbReadLimits.DevcData;
  newgms.DVCatomQ           = newgms.atomQ.DevcData;
  newgms.DVCatomMass        = newgms.atomMass.DevcData;
  newgms.DVCatomHDTM        = newgms.atomHDTM.DevcData;
  newgms.DVCinvNDF          = newgms.invNDF.DevcData;
  newgms.DVCatomLJID        = newgms.atomLJID.DevcData;
  newgms.DVCtypeCounts      = newgms.typeCounts.DevcData;
  newgms.DVCljFtab          = newgms.ljFtab.DevcData;
  newgms.DVCljUtab          = newgms.ljUtab.DevcData;
  if (gms->igb !=6 ) {
    newgms.DVCNeckID          = newgms.NeckID.DevcData;
    newgms.DVCrborn           = newgms.rborn.DevcData;
    newgms.DVCreff            = newgms.reff.DevcData;
    newgms.DVCGBalpha         = newgms.GBalpha.DevcData;
    newgms.DVCGBbeta          = newgms.GBbeta.DevcData;
    newgms.DVCGBgamma         = newgms.GBgamma.DevcData;
    newgms.DVCFscreen         = newgms.Fscreen.DevcData;
    newgms.DVCneckFactors     = newgms.neckFactors.DevcData;
  }
  newgms.DVCatomX           = newgms.atomX.DevcData;
  newgms.DVCatomY           = newgms.atomY.DevcData;
  newgms.DVCatomZ           = newgms.atomZ.DevcData;
  newgms.DVCvelcX           = newgms.velcX.DevcData;
  newgms.DVCvelcY           = newgms.velcY.DevcData;
  newgms.DVCvelcZ           = newgms.velcZ.DevcData;
  newgms.DVCfrcX            = newgms.frcX.DevcData;
  newgms.DVCfrcY            = newgms.frcY.DevcData;
  newgms.DVCfrcZ            = newgms.frcZ.DevcData;
  if (gms->Tstat.active == 1) {
    newgms.DVClangX           = newgms.langX.DevcData;
    newgms.DVClangY           = newgms.langY.DevcData;
    newgms.DVClangZ           = newgms.langZ.DevcData;
#ifdef CUDA
    newgms.DVCprngHeap        = newgms.prngHeap.DevcData;
    newgms.DVCprngReadLimits  = newgms.prngReadLimits.DevcData;
#endif
  }
  newgms.DVCuelec           = newgms.Uelec.DevcData;
  newgms.DVCuvdw            = newgms.Uvdw.DevcData;
  newgms.DVCusolv           = newgms.Usolv.DevcData;
  newgms.DVCubond           = newgms.Ubond.DevcData;
  newgms.DVCuangl           = newgms.Uangl.DevcData;
  newgms.DVCudihe           = newgms.Udihe.DevcData;
  newgms.DVCukine           = newgms.Ukine.DevcData;
  newgms.DVCtemp            = newgms.Temp.DevcData;
  newgms.DVCsystemPos       = newgms.systemPos.DevcData;

  return newgms;
}

//-----------------------------------------------------------------------------
// DestroyGpuMultiSim: free all memory associated with a GPU Multi-Simulation
//                     struct.
//
// Arguments:
//   gms:     the multi-simulation struct to free
//-----------------------------------------------------------------------------
static void DestroyGpuMultiSim(gpuMultiSim *gms)
{
  DestroyGpuInt2(&gms->atomReadLimits);
  DestroyGpuInt2(&gms->bondReadLimits);
  DestroyGpuUInt(&gms->bondIDs);
  DestroyGpuFloat2(&gms->bondBasics);
  DestroyGpuFloat4(&gms->bondAugs);
  if (gms->rattle == 1) {
    DestroyGpuInt2(&gms->cnstReadLimits);
    DestroyGpuUInt2(&gms->cnstInfo);
  }
  DestroyGpuInt2(&gms->anglReadLimits);
  DestroyGpuUInt(&gms->anglIDs);
  DestroyGpuFloat2(&gms->anglInfo);
  DestroyGpuInt2(&gms->diheReadLimits);
  DestroyGpuUInt2(&gms->diheIDs);
  DestroyGpuFloat2(&gms->diheInfo);
  DestroyGpuUInt(&gms->pairBitMasks);
  DestroyGpuUInt(&gms->attnIDs);
  DestroyGpuFloat2(&gms->attnFactors);
  DestroyGpuInt4(&gms->nbReadLimits);
  DestroyGpuFloat(&gms->atomQ);
  DestroyGpuFloat(&gms->atomMass);
  DestroyGpuFloat(&gms->atomHDTM);
  DestroyGpuDouble(&gms->invNDF);
  DestroyGpuInt(&gms->atomLJID);
  DestroyGpuInt(&gms->typeCounts);
  DestroyGpuFloat2(&gms->ljFtab);
  DestroyGpuFloat2(&gms->ljUtab);
  if (gms->igb != 6) {
    DestroyGpuInt(&gms->NeckID);
    DestroyGpuFloat(&gms->rborn);
    DestroyGpuFloat(&gms->reff);
    DestroyGpuFloat(&gms->GBalpha);
    DestroyGpuFloat(&gms->GBbeta);
    DestroyGpuFloat(&gms->GBgamma);
    DestroyGpuFloat(&gms->Fscreen);
    DestroyGpuFloat2(&gms->neckFactors);
  }
  DestroyGpuDouble(&gms->Uelec);
  DestroyGpuDouble(&gms->Uvdw);
  DestroyGpuDouble(&gms->Usolv);
  DestroyGpuDouble(&gms->Ubond);
  DestroyGpuDouble(&gms->Uangl);
  DestroyGpuDouble(&gms->Udihe);
  DestroyGpuDouble(&gms->Ukine);
#if 0
  DestroyGpuDouble(&gms->Temp);
#endif
  DestroyGpuFloat(&gms->atomX);
  DestroyGpuFloat(&gms->atomY);
  DestroyGpuFloat(&gms->atomZ);
  DestroyGpuFloat(&gms->crdrX);
  DestroyGpuFloat(&gms->crdrY);
  DestroyGpuFloat(&gms->crdrZ);
  DestroyGpuFloat(&gms->velcX);
  DestroyGpuFloat(&gms->velcY);
  DestroyGpuFloat(&gms->velcZ);
  DestroyGpuFloat(&gms->velrX);
  DestroyGpuFloat(&gms->velrY);
  DestroyGpuFloat(&gms->velrZ);
  DestroyGpuInt(&gms->frcX);
  DestroyGpuInt(&gms->frcY);
  DestroyGpuInt(&gms->frcZ);
  if (gms->Tstat.active == 1) {
    DestroyGpuInt(&gms->langX);
    DestroyGpuInt(&gms->langY);
    DestroyGpuInt(&gms->langZ);
    DestroyGpuFloat(&gms->Tstat.sdfac);
#ifdef CUDA
    DestroyGpuFloat(&gms->prngHeap);
    DestroyGpuInt2(&gms->prngReadLimits);
#endif
  }
}

//-----------------------------------------------------------------------------
// AllocateSlmem: allocate data for an slmem struct to serve as temporary
//                scratch space during simulations, replacing the expensive
//                calls to allocate and free such memory during force
//                computations.
//-----------------------------------------------------------------------------
static slmem AllocateSlmem()
{
  slmem sta;

  sta.ixfrc     = (int*)malloc(1024 * sizeof(int));
  sta.iyfrc     = (int*)malloc(1024 * sizeof(int));
  sta.izfrc     = (int*)malloc(1024 * sizeof(int));
  sta.ixbuff    = (int*)malloc(1024 * sizeof(int));
  sta.iybuff    = (int*)malloc(1024 * sizeof(int));
  sta.izbuff    = (int*)malloc(1024 * sizeof(int));
  sta.active    = (int*)malloc(1024 * sizeof(int));
  sta.ljid      = (int*)malloc(1024 * sizeof(int));
  sta.neckid    = (int*)malloc(1024 * sizeof(int));
  sta.xcrd      = (float*)malloc(1024 * sizeof(float));
  sta.ycrd      = (float*)malloc(1024 * sizeof(float));
  sta.zcrd      = (float*)malloc(1024 * sizeof(float));
  sta.xprvcrd   = (float*)malloc(1024 * sizeof(float));
  sta.yprvcrd   = (float*)malloc(1024 * sizeof(float));
  sta.zprvcrd   = (float*)malloc(1024 * sizeof(float));
  sta.qval      = (float*)malloc(1024 * sizeof(float));
  sta.rborn     = (float*)malloc(1024 * sizeof(float));
  sta.reff      = (float*)malloc(1024 * sizeof(float));
  sta.sumdeijda = (float*)malloc(1024 * sizeof(float));
  sta.gbalpha   = (float*)malloc(1024 * sizeof(float));
  sta.gbbeta    = (float*)malloc(1024 * sizeof(float));
  sta.gbgamma   = (float*)malloc(1024 * sizeof(float));
  sta.fs        = (float*)malloc(1024 * sizeof(float));
  sta.psi       = (int*)malloc(1024 * sizeof(int));

  return sta;
}

//-----------------------------------------------------------------------------
// DestroySlmem: free data associated with a "shared local memory" mockup.
//
// Arguments:
//   sta:     the mockup to free
//-----------------------------------------------------------------------------
static void DestroySlmem(slmem *sta)
{
  free(sta->ixfrc);
  free(sta->iyfrc);
  free(sta->izfrc);
  free(sta->ixbuff);
  free(sta->iybuff);
  free(sta->izbuff);
  free(sta->active);
  free(sta->ljid);
  free(sta->xcrd);
  free(sta->ycrd);
  free(sta->zcrd);
  free(sta->xprvcrd);
  free(sta->yprvcrd);
  free(sta->zprvcrd);
  free(sta->qval);
  free(sta->neckid);
  free(sta->rborn);
  free(sta->reff);
  free(sta->sumdeijda);
  free(sta->gbalpha);
  free(sta->gbbeta);
  free(sta->gbgamma);
  free(sta->fs);
  free(sta->psi);
}

//-----------------------------------------------------------------------------
// CrossPf: function for finding the cross-product cr of vectors p and q in
//          fp32.  Note that vectors p and q are assumed to be
//          three-dimensional and only the first three components of these
//          vectors will be considered.
//-----------------------------------------------------------------------------
static void CrossPf(float* p, float* q, float* cr)
{
  cr[0] = p[1]*q[2] - p[2]*q[1];
  cr[1] = p[2]*q[0] - p[0]*q[2];
  cr[2] = p[0]*q[1] - p[1]*q[0];
}

//-----------------------------------------------------------------------------
// GmsFastStepEnergyAndForces: compute the energy and forces due to bond and
//                             angle terms, using data in host RAM and the CPU.
//                             This function pairs with the non-bonded and
//                             dihedral calculations to complete a multiple
//                             time-step dynamics engine.
//
// Arguments:
//   gms:       information for all systems, including coordinates and
//              topologies
//   sta:       System Temporary Allocations, to hold data on any one system's
//              coordinates, velocities, and forces.  This temporary scratch
//              space is pre-allocated once for the entirety of the simulation.
//   resetfrc:  flag to have forces reset (not used on the major step when
//              non-bonded and dihedral forces have been computed, but on the
//              minor step when bonds and angles are being updated alone)
//   stepidx:   step count within the segment (not total step count)
//-----------------------------------------------------------------------------
static void GmsFastStepEnergyAndForces(gpuMultiSim *gms, slmem sta,
                                       int resetfrc, int stepidx)
{
  int h, i, j, xtl, ytl, jatm, katm, matm, nrgidx;
  int natom, atmstart, ifx, ify, ifz, ljTabIdx, ntypes;
#ifdef AMBER_PLATFORM_AMD
  uint2 jmask;
#else
  unsigned int jmask;
#endif
  float atmx, atmy, atmz, fmag;
  float dx, dy, dz, dl, dlpu, dlpr, r2, r;
  float Leq, Keq, Lpull, Kpull, Lpress, Kpress;
  float bax, bay, baz, bcx, bcy, bcz;
  float theta, costheta, dtheta, mgba, mgbc, dA, mbabc, sqba, sqbc, invbabc;
  float adfx, adfy, adfz, cdfx, cdfy, cdfz;
  float elscale, ljscale, invr, invr2, invr4, invr6, invr8, qq;
  float2 ljterm;
  float ubond, uangl, uelec, uvdw;

  // Loop over all systems
  for (h = 0; h < gms->nsys; h++) {
    atmstart = gms->atomReadLimits.HostData[h].x;
    natom = gms->atomReadLimits.HostData[h].y - atmstart;
    ljTabIdx = h * gms->ljABoffset;
    ntypes = gms->typeCounts.HostData[h];

    // Import coordinates and properties into local arrays.
    // Initialize local force accumulators.
    for (i = 0; i < natom; i++) {
      sta.xcrd[i] = gms->atomX.HostData[i + atmstart];
      sta.ycrd[i] = gms->atomY.HostData[i + atmstart];
      sta.zcrd[i] = gms->atomZ.HostData[i + atmstart];
      sta.qval[i] = gms->atomQ.HostData[i + atmstart];
      sta.ljid[i] = gms->atomLJID.HostData[i + atmstart];
      if (resetfrc == 1) {
        sta.ixfrc[i] = 0;
        sta.iyfrc[i] = 0;
        sta.izfrc[i] = 0;
      }
      else {
        sta.ixfrc[i] = gms->frcX.HostData[i + atmstart];
        sta.iyfrc[i] = gms->frcY.HostData[i + atmstart];
        sta.izfrc[i] = gms->frcZ.HostData[i + atmstart];
      }
    }

    // Do non-bonded attenuations
    uelec = 0.0;
    uvdw = 0.0;
    for (i = gms->nbReadLimits.HostData[h].z;
         i < gms->nbReadLimits.HostData[h].w; i++) {
      jatm = (gms->attnIDs.HostData[i] >> 16);
      katm = (gms->attnIDs.HostData[i] & 0xffff);
      elscale = gms->attnFactors.HostData[i].x;
      ljscale = gms->attnFactors.HostData[i].y;

      // Compute the non-bonded interaction of these two atoms
      dx = sta.xcrd[katm] - sta.xcrd[jatm];
      dy = sta.ycrd[katm] - sta.ycrd[jatm];
      dz = sta.zcrd[katm] - sta.zcrd[jatm];
      r2 = dx*dx + dy*dy + dz*dz;
      invr2 = 1.0 / r2;
      invr4 = invr2 * invr2;
      invr6 = invr4 * invr2;
      invr8 = invr4 * invr4;
      r = sqrt(r2);
      invr = 1.0 / r;
      qq = elscale * sta.qval[katm] * sta.qval[jatm];
      ljterm = gms->ljFtab.HostData[ljTabIdx + sta.ljid[jatm]*ntypes +
                                    sta.ljid[katm]];
      fmag = -(qq * invr2 * invr) +
             (ljscale * invr8 * ((ljterm.x * invr6) + ljterm.y));
      fmag *= FPSCALEfrc;

      // Compute the energy
      ljterm = gms->ljUtab.HostData[ljTabIdx + sta.ljid[jatm]*ntypes +
                                    sta.ljid[katm]];
      uelec += qq * invr;
      uvdw  += ljscale * ((ljterm.x * invr6) + ljterm.y) * invr6;
      ifx = (int)(fmag * dx);
      ify = (int)(fmag * dy);
      ifz = (int)(fmag * dz);
      sta.ixfrc[jatm] += ifx;
      sta.iyfrc[jatm] += ify;
      sta.izfrc[jatm] += ifz;
      sta.ixfrc[katm] -= ifx;
      sta.iyfrc[katm] -= ify;
      sta.izfrc[katm] -= ifz;
    }

    // Loop over all bonded interactions
    ubond = 0.0;
    for (i = gms->bondReadLimits.HostData[h].x;
         i < gms->bondReadLimits.HostData[h].y; i++) {
      jatm = gms->bondIDs.HostData[i];
      katm = (jatm >> 16);
      jatm = (jatm & 0xffff);
      dx = sta.xcrd[katm] - sta.xcrd[jatm];
      dy = sta.ycrd[katm] - sta.ycrd[jatm];
      dz = sta.zcrd[katm] - sta.zcrd[jatm];
      r2 = dx*dx + dy*dy + dz*dz;
      r = sqrt(r2);
      Keq    = gms->bondBasics.HostData[i].x;
      Leq    = gms->bondBasics.HostData[i].y;
      Kpull  = gms->bondAugs.HostData[i].x;
      Lpull  = gms->bondAugs.HostData[i].y;
      Kpress = gms->bondAugs.HostData[i].z;
      Lpress = gms->bondAugs.HostData[i].w;
      dl = (Leq >= (float)0.0) ? Leq - r : r + Leq;
      dlpu = (r > Lpull)  ? Lpull  - r : (float)0.0;
      dlpr = (r < Lpress) ? Lpress - r : (float)0.0;
      ubond += Keq*dl*dl + Kpull*dlpu*dlpu + Kpress*dlpr*dlpr;
      fmag = -2.0 * FPSCALEfrc * (Keq*dl + Kpull*dlpu + Kpress*dlpr) / r;
      ifx = (int)(fmag * dx);
      ify = (int)(fmag * dy);
      ifz = (int)(fmag * dz);
      sta.ixfrc[jatm] += ifx;
      sta.iyfrc[jatm] += ify;
      sta.izfrc[jatm] += ifz;
      sta.ixfrc[katm] -= ifx;
      sta.iyfrc[katm] -= ify;
      sta.izfrc[katm] -= ifz;
    }

    // Loop over all angle interactions
    uangl = 0.0;
    for (i = gms->anglReadLimits.HostData[h].x;
         i < gms->anglReadLimits.HostData[h].y; i++) {
      jatm = gms->anglIDs.HostData[i];
      katm = ((jatm >> 10) & 0x3ff);
      matm = (jatm & 0x3ff);
      jatm = (jatm >> 20);
      bax = sta.xcrd[jatm] - sta.xcrd[katm];
      bay = sta.ycrd[jatm] - sta.ycrd[katm];
      baz = sta.zcrd[jatm] - sta.zcrd[katm];
      bcx = sta.xcrd[matm] - sta.xcrd[katm];
      bcy = sta.ycrd[matm] - sta.ycrd[katm];
      bcz = sta.zcrd[matm] - sta.zcrd[katm];
      mgba = bax*bax + bay*bay + baz*baz;
      mgbc = bcx*bcx + bcy*bcy + bcz*bcz;
      invbabc = 1.0/sqrt(mgba*mgbc);
      costheta = (bax*bcx + bay*bcy + baz*bcz) * invbabc;
      costheta = (costheta < (float)-1.0) ? (float)-1.0 :
                 (costheta > (float)1.0)  ? (float)1.0  : costheta;
      theta = acos(costheta);
      Keq = gms->anglInfo.HostData[i].x;
      Leq = gms->anglInfo.HostData[i].y;
      dtheta = theta - Leq;
      uangl += Keq * dtheta * dtheta;
      dA = (float)-2.0 * Keq * dtheta / sqrt((float)1.0 - costheta*costheta);
      sqba = dA / mgba;
      sqbc = dA / mgbc;
      mbabc = dA * invbabc;
      adfx = (bcx*mbabc - costheta*bax*sqba) * FPSCALEfrc;
      cdfx = (bax*mbabc - costheta*bcx*sqbc) * FPSCALEfrc;
      adfy = (bcy*mbabc - costheta*bay*sqba) * FPSCALEfrc;
      cdfy = (bay*mbabc - costheta*bcy*sqbc) * FPSCALEfrc;
      adfz = (bcz*mbabc - costheta*baz*sqba) * FPSCALEfrc;
      cdfz = (baz*mbabc - costheta*bcz*sqbc) * FPSCALEfrc;
      sta.ixfrc[jatm] -= adfx;
      sta.iyfrc[jatm] -= adfy;
      sta.izfrc[jatm] -= adfz;
      sta.ixfrc[katm] += adfx + cdfx;
      sta.iyfrc[katm] += adfy + cdfy;
      sta.izfrc[katm] += adfz + cdfz;
      sta.ixfrc[matm] -= cdfx;
      sta.iyfrc[matm] -= cdfy;
      sta.izfrc[matm] -= cdfz;
    }

    // Commit forces to "global"
    for (i = 0; i < natom; i++) {
      gms->frcX.HostData[i + atmstart] = sta.ixfrc[i];
      gms->frcY.HostData[i + atmstart] = sta.iyfrc[i];
      gms->frcZ.HostData[i + atmstart] = sta.izfrc[i];
    }

    // Commit energies to "global"
    if (stepidx % gms->ntpr == 0 && resetfrc == 0) {
      nrgidx = (h * gms->nsgmdout) + (stepidx / gms->ntpr);
      gms->Ubond.HostData[nrgidx]  = ubond;
      gms->Uangl.HostData[nrgidx]  = uangl;
      gms->Uelec.HostData[nrgidx] += uelec;
      gms->Uvdw.HostData[nrgidx]  += uvdw;
    }
  }
}

//-----------------------------------------------------------------------------
// GmsSlowStepEnergyAndForces: compute the non-bonded and dihedral energy of
//                             all systems, using data in host RAM and the CPU.
//                             The function GmsFastStepEnergyAndForces works
//                             to toss in the high-frequency degrees of freedom
//                             at much higher update rates.
//
// Arguments:
//   gms:      information for all systems, including coordinates and
//             topologies
//   sta:      System Temporary Allocations, to hold data on any one system's
//             coordinates, velocities, and forces.  This temporary scratch
//             space is pre-allocated once for the entirety of the simulation.
//   stepidx:  step count within the segment (not total step count)
//-----------------------------------------------------------------------------
static void GmsSlowStepEnergyAndForces(gpuMultiSim *gms, slmem sta,
                                       int stepidx)
{
  int h, i, j, k, m, xtl, ytl, jatm, katm, matm, patm, jljt;
  int natom, ntile, ntypes, atmstart, ljTabIdx, nrgidx;
  int ifx, ify, ifz;
#ifdef AMBER_PLATFORM_AMD
  uint2 jmask;
#else
  unsigned int jmask;
#endif
  float atmx, atmy, atmz, atmq, qq, fmag, elscale, ljscale, Keq, Leq;
  float dx, dy, dz, r2, r, invr, invr2, invr3, invr4, invr6, invr8;
  float invrad, invirad, invjrad, sumi, si, si2, sj, sj2, uij, tmpsd, dumbo;
  float sumda, thi, datmpj, datmpk, neck, fipsi;
  float mdist, mdist2, mdist6, atmrad, atmirad, atmjrad, rb2, atmq2h, atmqd2h;
  float expmkf, efac, fgbi, fgbk, dielfac, temp1, temp4, temp5, temp6;
  float theta, costheta, mgbc;
  float fa, fb1, fc1, fb2, fc2, fd, isinb2, isinc2, Npdc, sangle;
  float mgab, mgcd, invab, invbc, invcd, invabc, invbcd, cosb, cosc;
  float afmag, bfmag, cfmag, dfmag;
  float ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  float2 ljterm;
  double uelec, uvdw, udihe, usolv;

  // Loop over all systems
  for (h = 0; h < gms->nsys; h++) {
    atmstart = gms->atomReadLimits.HostData[h].x;
    natom = gms->atomReadLimits.HostData[h].y - atmstart;
    ljTabIdx = h * gms->ljABoffset;

    // Import coordinates and properties into local arrays.
    // Initialize local force accumulators.
    for (i = 0; i < natom; i++) {
      sta.xcrd[i]    = gms->atomX.HostData[i + atmstart];
      sta.ycrd[i]    = gms->atomY.HostData[i + atmstart];
      sta.zcrd[i]    = gms->atomZ.HostData[i + atmstart];
      sta.qval[i]    = gms->atomQ.HostData[i + atmstart];
      sta.ljid[i]    = gms->atomLJID.HostData[i + atmstart];
      sta.ixfrc[i]   = 0;
      sta.iyfrc[i]   = 0;
      sta.izfrc[i]   = 0;
    }
    if (gms->igb == 1 || gms->igb == 2 || gms->igb == 5 || gms->igb == 7 ||
        gms->igb == 8) {
      for (i = 0; i < natom; i++) {
        sta.rborn[i]   = gms->rborn.HostData[i + atmstart];
        sta.neckid[i]  = gms->NeckID.HostData[i + atmstart];
        sta.gbalpha[i] = gms->GBalpha.HostData[i + atmstart];
        sta.gbbeta[i]  = gms->GBbeta.HostData[i + atmstart];
        sta.gbgamma[i] = gms->GBgamma.HostData[i + atmstart];
        sta.fs[i]      = gms->Fscreen.HostData[i + atmstart];
      }
    }

    // Compute GB radii
    if (gms->igb != 6) {

      // Prime the psi accumulator array
      for (i = 0; i < natom; i++) {
        sta.psi[i] = 0;
      }

      // Loop over all pairs, reusing displacement and other information
      for (i = 1; i < natom; i++) {
        atmx = sta.xcrd[i];
        atmy = sta.ycrd[i];
        atmz = sta.zcrd[i];
        atmirad = sta.rborn[i] - gms->GBOffset;
        invirad = (float)1.0 / atmirad;
        for (j = 0; j < i; j++) {
          dx = atmx - sta.xcrd[j];
          dy = atmy - sta.ycrd[j];
          dz = atmz - sta.zcrd[j];
          atmjrad = sta.rborn[j] - gms->GBOffset;
          invjrad = (float)1.0 / atmjrad;
          r2   = dx*dx + dy*dy + dz*dz;
          invr = (float)1.0 / sqrt(r2);
          r    = r2 * invr;

          // First computation: i -> j
          sj = sta.fs[j] * atmjrad;
          sj2 = sj * sj;
          if (r > (float)4.0 * sj) {
            invr2  = invr * invr;
            tmpsd  = sj2 * invr2;
            dumbo  = TA + tmpsd*(TB + tmpsd*(TC + tmpsd*(TD + tmpsd*TDD)));
            sta.psi[i] -= (int)(sj * tmpsd * invr2 * dumbo * FPSCALErad);
          }
          else if (r > atmirad + sj) {
            sta.psi[i] -= (int)((float)0.5 * FPSCALErad *
                                ((sj / (r2 - sj2)) +
                                 ((float)0.5 * invr * log((r - sj) /
                                                          (r + sj)))));
          }
          else if (r > fabs(atmirad - sj)) {
            theta = (float)0.5 * invirad * invr * (r2 + atmirad*atmirad - sj2);
            uij = (float)1.0 / (r + sj);
            sta.psi[i] -= (int)((float)0.25 * FPSCALErad *
                                (invirad*((float)2.0 - theta) - uij +
                                 invr*log(atmirad * uij)));
          }
          else if (atmirad < sj) {
            sta.psi[i] -= (int)((float)0.5 * FPSCALErad *
                                ((sj / (r2 - sj2)) +
                                 ((float)2.0 * invirad) +
                                 ((float)0.5 * invr * log((sj - r) /
                                                          (sj + r)))));
          }

          // Second computation: j -> i
          si = sta.fs[i] * atmirad;
          si2 = si * si;
          if (r > (float)4.0 * si) {
            invr2  = invr * invr;
            tmpsd  = si2 * invr2;
            dumbo  = TA + tmpsd*(TB + tmpsd*(TC + tmpsd*(TD + tmpsd*TDD)));
            sta.psi[j] -= (int)(si * tmpsd * invr2 * dumbo * FPSCALErad);
          }
          else if (r > atmjrad + si) {
            sta.psi[j] -= (int)((float)0.5 * FPSCALErad *
                                ((si / (r2 - si2)) +
                                 ((float)0.5 * invr * log((r - si) /
                                                          (r + si)))));
          }
          else if (r > fabs(atmjrad - si)) {
            theta = (float)0.5 * invjrad * invr * (r2 + atmjrad*atmjrad - si2);
            uij = (float)1.0 / (r + si);
            sta.psi[j] -= (int)((float)0.25 * FPSCALErad *
                                (invjrad*((float)2.0 - theta) - uij +
                                 invr*log(atmjrad * uij)));
          }
          else if (atmjrad < si) {
            sta.psi[j] -= (int)((float)0.5 * FPSCALErad *
                                ((si / (r2 - si2)) +
                                 ((float)2.0 * invjrad) +
                                 ((float)0.5 * invr * log((si - r) /
                                                          (si + r)))));
          }

          // Neck GB contribution
          if ((gms->igb == 7 || gms->igb == 8) &&
              r < sta.rborn[i] + sta.rborn[j] + GBNECKCUT) {

            // First computation: i -> j
            mdist  = r - neckMaxPos[sta.neckid[i]][sta.neckid[j]];
            mdist2 = mdist * mdist;
            mdist6 = mdist2 * mdist2 * mdist2;
            neck   = neckMaxVal[sta.neckid[i]][sta.neckid[j]] /
                     ((float)1.0 + mdist2 + (float)0.3*mdist6);
            sta.psi[i] -= (int)(FPSCALErad * gms->GBNeckScale * neck);

            // Second computation: j -> i
            mdist  = r - neckMaxPos[sta.neckid[j]][sta.neckid[i]];
            mdist2 = mdist * mdist;
            mdist6 = mdist2 * mdist2 * mdist2;
            neck   = neckMaxVal[sta.neckid[j]][sta.neckid[i]] /
                     ((float)1.0 + mdist2 + (float)0.3*mdist6);
            sta.psi[j] -= (int)(FPSCALErad * gms->GBNeckScale * neck);
          }
        }
      }

      // Final pass to get effective Born radii
      for (i = 0; i < natom; i++) {
        fipsi = (float)(sta.psi[i]) / FPSCALErad;
        if (gms->igb == 1) {

          // Original (Hawkins-Craemer-Truhlar) effective radii
          invirad = (float)1.0 / (sta.rborn[i] - gms->GBOffset);
          sta.reff[i] = (float)1.0 / (invirad + fipsi);
          if (sta.reff[i] < (float)0.0) {
            sta.reff[i] = (float)30.0;
          }
        }
        else {

          // "GBAO" formulas
          atmirad = sta.rborn[i] - gms->GBOffset;
          invirad = (float)1.0 / atmirad;
          fipsi *= -atmirad;
          sta.reff[i] = (float)1.0 /
                        (invirad - tanh((sta.gbalpha[i] -
                                         (sta.gbbeta[i] * fipsi) +
                                         (sta.gbgamma[i] * fipsi * fipsi)) *
                                        fipsi) / sta.rborn[i]);
        }
      }
    }

    // Initialize for Generalized Born solvent forces
    usolv = 0.0;
    for (i = 0; i < natom; i++) {
      if (gms->igb != 6) {
        atmq = sta.qval[i];
        atmrad = sta.reff[i];
        expmkf = exp(-KSCALE * (gms->kappa) * atmrad) / gms->dielectric;
        dielfac = (float)1.0 - expmkf;
        atmq2h = (float)0.5 * atmq * atmq;
        atmqd2h = atmq2h * dielfac;
        usolv += -atmqd2h / atmrad;
        sta.sumdeijda[i] = atmqd2h - (KSCALE * gms->kappa * atmq2h *
                                      expmkf * atmrad);
      }
      else {
        sta.sumdeijda[i] = (float)0.0;
      }
    }

    // Loop over non-bonded tiles, reading from the pair list bit masks and
    // computing any non-bonded interactions that are not fully excluded.
    ntile = (natom + HALFGRID_MASK) / HALFGRID;
    ntypes = gms->typeCounts.HostData[h];
    xtl = 0;
    ytl = 0;
    uelec = 0.0;
    uvdw = 0.0;
    for (i = gms->nbReadLimits.HostData[h].x;
         i < gms->nbReadLimits.HostData[h].y; i += GRID) {
      for (j = 0; j < HALFGRID; j++) {
        jatm = xtl*HALFGRID + j;
        if (jatm >= natom) {
          continue;
        }

        // Shortcuts and pre-computations for this atom
        atmx = sta.xcrd[jatm];
        atmy = sta.ycrd[jatm];
        atmz = sta.zcrd[jatm];
        atmq = sta.qval[jatm];
        jljt = sta.ljid[jatm] * ntypes;
        jmask = gms->pairBitMasks.HostData[i + j];
        if (gms->igb != 6) {
          atmrad = sta.reff[jatm];
        }

        // Inner loop for this tile
        for (k = 0; k < HALFGRID; k++) {
          katm = ytl*HALFGRID + k;
          if (katm < jatm) {

            // Generalized Born computations (all pairs participate--
            // save the distance in case it is needed for Lennard-Jones
            // and Coulomb non-bonded interactions)
            dx = sta.xcrd[katm] - atmx;
            dy = sta.ycrd[katm] - atmy;
            dz = sta.zcrd[katm] - atmz;
            r2 = dx*dx + dy*dy + dz*dz;
            qq = sta.qval[katm] * atmq;
            if (gms->igb != 6) {
              rb2 = atmrad * sta.reff[katm];
              efac = exp(-r2 / ((float)4.0 * rb2));
              fgbi = (float)1.0 / sqrt(r2 + rb2 * efac);
              fgbk = -gms->kappa * KSCALE / fgbi;
              expmkf = exp(fgbk) / gms->dielectric;
              dielfac = (float)1.0 - expmkf;
              usolv += -qq * dielfac * fgbi;
              temp4 = fgbi * fgbi * fgbi;
              temp6 = qq * temp4 * (dielfac + fgbk * expmkf);
              fmag = temp6 * ((float)1.0 - (float)0.25 * efac);
              temp5 = (float)0.5 * efac * temp6 * (rb2 + (float)0.25*r2);
              sta.sumdeijda[jatm] += atmrad * temp5;
              sta.sumdeijda[katm] += sta.reff[katm] * temp5;
            }
            else {
              fmag = (float)0.0;
            }

            // Compute the non-bonded interaction of these two atoms
#ifdef AMBER_PLATFORM_AMD
              if (((jmask.y >> k) & 0x1) == 0x1) {
#else
              if (((jmask >> k) & 0x1) == 0x1) {
#endif
              invr = (float)1.0 / sqrt(r2);
              invr2 = invr  * invr;
              invr4 = invr2 * invr2;
              invr8 = invr4 * invr4;
              ljterm = gms->ljFtab.HostData[ljTabIdx + jljt + sta.ljid[katm]];
              fmag += -(qq * invr2 * invr) +
                      invr8*((ljterm.x * invr4 * invr2) + ljterm.y);

              // Compute the energy
              ljterm = gms->ljUtab.HostData[ljTabIdx + jljt + sta.ljid[katm]];
              uelec += qq * invr;
              uvdw  += ((ljterm.x * invr4 * invr2) + ljterm.y) * invr4 * invr2;
            }

            // Contribute forces
            fmag *= FPSCALEfrc;
            ifx = (int)(fmag * dx);
            ify = (int)(fmag * dy);
            ifz = (int)(fmag * dz);
            sta.ixfrc[jatm] += ifx;
            sta.iyfrc[jatm] += ify;
            sta.izfrc[jatm] += ifz;
            sta.ixfrc[katm] -= ifx;
            sta.iyfrc[katm] -= ify;
            sta.izfrc[katm] -= ifz;
          }
        }
      }
      if (ytl == xtl) {
        ytl = 0;
        xtl++;
      }
      else {
        ytl++;
      }
    }

    // Loop back over non-bonded tiles to fold in derivatives of
    // the effective Born radii to the forces on each atom
    if (gms->igb != 6) {

      // Work with the Born radii derivative array to save computations
      // in the inner loop
      for (i = 0; i < natom; i++) {
        atmirad = sta.rborn[i] - gms->GBOffset;
        fipsi = (float)(sta.psi[i]) / FPSCALErad;
        if (gms->igb != 1) {
          fipsi *= -atmirad;
        }
        thi = tanh((sta.gbalpha[i] - sta.gbbeta[i] * fipsi +
                    sta.gbgamma[i] * fipsi * fipsi) *
                   fipsi);
        sta.sumdeijda[i] *= (sta.gbalpha[i] -
                             ((float)2.0 * sta.gbbeta[i] * fipsi) +
                             ((float)3.0 * sta.gbgamma[i] * fipsi * fipsi)) *
                            ((float)1.0 - thi * thi) * atmirad / sta.rborn[i];
      }

      // Nested loop over tiles
      xtl = 0;
      ytl = 0;
      for (i = gms->nbReadLimits.HostData[h].x;
           i < gms->nbReadLimits.HostData[h].y; i += GRID) {
        for (j = 0; j < HALFGRID; j++) {
          jatm = xtl*HALFGRID + j;
          if (jatm >= natom) {
            continue;
          }

          // Shortcuts and pre-computations for this atom
          atmx = sta.xcrd[jatm];
          atmy = sta.ycrd[jatm];
          atmz = sta.zcrd[jatm];
          atmirad = sta.rborn[jatm] - gms->GBOffset;
          invirad = (float)1.0 / atmirad;
          sumda = sta.sumdeijda[jatm];

          // Inner loop for this tile
          for (k = 0; k < HALFGRID; k++) {
            katm = ytl*HALFGRID + k;
            if (katm < jatm) {
              dx = sta.xcrd[katm] - atmx;
              dy = sta.ycrd[katm] - atmy;
              dz = sta.zcrd[katm] - atmz;
              atmjrad = sta.rborn[katm] - gms->GBOffset;
              invjrad = (float)1.0 / atmjrad;
              r2 = dx*dx + dy*dy + dz*dz;
              invr = (float)1.0 / sqrt(r2);
              invr2 = invr * invr;
              r = r2 * invr;

              // First computation: i -> j
              sj = sta.fs[katm] * atmjrad;
              sj2 = sj * sj;
              if (r > (float)4.0 * sj) {
                tmpsd  = sj2 * invr2;
                dumbo  = TE + tmpsd*(TF + tmpsd*(TG + tmpsd*(TH + tmpsd*THH)));
                datmpj = tmpsd * sj * invr2 * invr2 * dumbo;
              }
              else if (r > atmirad + sj) {
                temp1  = (float)1.0 / (r2 - sj2);
                datmpj = temp1 * sj * ((float)-0.5 * invr2 + temp1) +
                         (float)0.25 * invr * invr2 * log((r - sj) /
                                                          (r + sj));
              }
              else if (r > fabs(atmirad - sj)) {
                temp1  = (float)1.0 / (r + sj);
                invr3  = invr2 * invr;
                datmpj = (float)-0.25 * (((float)-0.5 *
                                          (r2 - atmirad*atmirad + sj2) *
                                          invr3 * invirad * invirad) +
                                         (invr * temp1 * (temp1 - invr)) -
                                         (invr3 * log(atmirad * temp1)));
              }
              else if (atmirad < sj) {
                temp1  = (float)1.0 / (r2 - sj2);
                datmpj = (float)-0.5 * ((sj * invr2 * temp1) -
                                        ((float)2.0 * sj * temp1 * temp1) -
                                        ((float)0.5 * invr2 * invr *
                                         log((sj - r) / (sj + r))));
              }
              else {
                datmpj = (float)0.0;
              }

              // Second computation: j -> i
              si = sta.fs[jatm] * atmirad;
              si2 = si * si;
              if (r > (float)4.0 * si) {
                tmpsd  = si2 * invr2;
                dumbo  = TE + tmpsd*(TF + tmpsd*(TG + tmpsd*(TH + tmpsd*THH)));
                datmpk = tmpsd * si * invr2 * invr2 * dumbo;
              }
              else if (r > atmjrad + si) {
                temp1  = (float)1.0 / (r2 - si2);
                datmpk = temp1 * si * ((float)-0.5 * invr2 + temp1) +
                         (float)0.25 * invr * invr2 * log((r - si) /
                                                          (r + si));
              }
              else if (r > fabs(atmjrad - si)) {
                temp1  = (float)1.0 / (r + si);
                invr3  = invr2 * invr;
                datmpk = (float)-0.25 * (((float)-0.5 *
                                          (r2 - atmjrad*atmjrad + si2) *
                                          invr3 * invjrad * invjrad) +
                                         (invr * temp1 * (temp1 - invr)) -
                                         (invr3 * log(atmjrad * temp1)));
              }
              else if (atmjrad < si) {
                temp1  = (float)1.0 / (r2 - si2);
                datmpk = (float)-0.5 * ((si * invr2 * temp1) -
                                        ((float)2.0 * si * temp1 * temp1) -
                                        ((float)0.5 * invr2 * invr *
                                         log((si - r) / (si + r))));
              }
              else {
                datmpk = (float)0.0;
              }

              // Neck GB contributions
              if ((gms->igb == 7 || gms->igb == 8) &&
                  r < sta.rborn[jatm] + sta.rborn[katm] + GBNECKCUT) {

                // First computation: i -> j
                mdist = r - neckMaxPos[sta.neckid[jatm]][sta.neckid[katm]];
                mdist2 = mdist * mdist;
                mdist6 = mdist2 * mdist2 * mdist2;
                temp1 = (float)1.0 + mdist2 + (float)0.3*mdist6;
                temp1 = temp1 * temp1 * r;
                datmpj += (((float)2.0*mdist +
                            (float)1.8*mdist2*mdist2*mdist) *
                           neckMaxVal[sta.neckid[jatm]][sta.neckid[katm]] *
                           gms->GBNeckScale) / temp1;

                // Second computation: j -> i
                mdist = r - neckMaxPos[sta.neckid[katm]][sta.neckid[jatm]];
                mdist2 = mdist * mdist;
                mdist6 = mdist2 * mdist2 * mdist2;
                temp1 = (float)1.0 + mdist2 + (float)0.3*mdist6;
                temp1 = temp1 * temp1 * r;
                datmpk += (((float)2.0*mdist +
                            (float)1.8*mdist2*mdist2*mdist) *
                           neckMaxVal[sta.neckid[katm]][sta.neckid[jatm]] *
                           gms->GBNeckScale) / temp1;
              }

              // Contribute the derivatives to the force arrays
              fmag  = (datmpj * sumda) + (datmpk * sta.sumdeijda[katm]);
              fmag *= FPSCALEfrc;
              ifx = (int)(fmag * dx);
              ify = (int)(fmag * dy);
              ifz = (int)(fmag * dz);
              sta.ixfrc[jatm] -= ifx;
              sta.iyfrc[jatm] -= ify;
              sta.izfrc[jatm] -= ifz;
              sta.ixfrc[katm] += ifx;
              sta.iyfrc[katm] += ify;
              sta.izfrc[katm] += ifz;
            }
          }
        }
        if (ytl == xtl) {
          ytl = 0;
          xtl++;
        }
        else {
          ytl++;
        }
      }
    }

    // Loop over all dihedral interactions
    udihe = 0.0;
    for (i = gms->diheReadLimits.HostData[h].x;
         i < gms->diheReadLimits.HostData[h].y; i++) {
      jatm = gms->diheIDs.HostData[i].x;
      katm = ((jatm >> 10) & 0x3ff);
      matm = (jatm & 0x3ff);
      jatm = (jatm >> 20);
      patm = gms->diheIDs.HostData[i].y;
      Npdc = (float)(patm & 0xffff);
      patm = (patm >> 16);
      ab[0] = sta.xcrd[katm] - sta.xcrd[jatm];
      ab[1] = sta.ycrd[katm] - sta.ycrd[jatm];
      ab[2] = sta.zcrd[katm] - sta.zcrd[jatm];
      bc[0] = sta.xcrd[matm] - sta.xcrd[katm];
      bc[1] = sta.ycrd[matm] - sta.ycrd[katm];
      bc[2] = sta.zcrd[matm] - sta.zcrd[katm];
      cd[0] = sta.xcrd[patm] - sta.xcrd[matm];
      cd[1] = sta.ycrd[patm] - sta.ycrd[matm];
      cd[2] = sta.zcrd[patm] - sta.zcrd[matm];
      CrossPf(ab, bc, crabbc);
      CrossPf(bc, cd, crbccd);
      if (crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2] <
          (float)1.0e-4 ||
          crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2] <
          (float)1.0e-4) {
        continue;
      }
      costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] +
                 crabbc[2]*crbccd[2];
      costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] +
                        crabbc[2]*crabbc[2]) *
                       (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] +
                        crbccd[2]*crbccd[2]));
      CrossPf(crabbc, crbccd, scr);
      costheta = (costheta < (float)-1.0) ?
                 (float)-1.0 : (costheta > (float)1.0) ? (float)1.0 : costheta;
      if (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > (float)0.0) {
        theta = acosf(costheta);
      }
      else {
        theta = -acosf(costheta);
      }
      Keq = gms->diheInfo.HostData[i].x;
      Leq = gms->diheInfo.HostData[i].y;
      sangle = Npdc*theta - Leq;
      fmag = FPSCALEfrc * Keq * Npdc * sinf(sangle);
      udihe += Keq * ((float)1.0 + cos(sangle));
      mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      invab = (float)1.0 / mgab;
      mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      invbc = (float)1.0 / mgbc;
      mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
      invcd = (float)1.0 / mgcd;
      cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2])*invab*invbc;
      isinb2 = (cosb*cosb < (float)0.9999) ?
               (float)1.0/((float)1.0 - cosb*cosb) : (float)0.0;
      cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2])*invbc*invcd;
      isinc2 = (cosc*cosc < (float)0.9999) ?
               (float)1.0/((float)1.0 - cosc*cosc) : (float)0.0;
      isinb2 *= fmag;
      isinc2 *= fmag;
      invabc = invab*invbc;
      invbcd = invbc*invcd;
      for (j = 0; j < 3; j++) {
        crabbc[j] *= invabc;
        crbccd[j] *= invbcd;
      }

      // Transform the dihedral forces to cartesian coordinates
      fa = -invab * isinb2;
      fb1 = (mgbc - mgab*cosb) * invabc * isinb2;
      fb2 = cosc * invbc * isinc2;
      fc1 = (mgbc - mgcd*cosc) * invbcd * isinc2;
      fc2 = cosb * invbc * isinb2;
      fd = -invcd * isinc2;

      // Apply the dihedral forces
      sta.ixfrc[jatm] += crabbc[0] * fa;
      sta.ixfrc[katm] += fb1 * crabbc[0] - fb2 * crbccd[0];
      sta.ixfrc[matm] += -fc1 * crbccd[0] + fc2 * crabbc[0];
      sta.ixfrc[patm] += -fd * crbccd[0];
      sta.iyfrc[jatm] += crabbc[1] * fa;
      sta.iyfrc[katm] += fb1 * crabbc[1] - fb2 * crbccd[1];
      sta.iyfrc[matm] += -fc1 * crbccd[1] + fc2 * crabbc[1];
      sta.iyfrc[patm] += -fd * crbccd[1];
      sta.izfrc[jatm] += crabbc[2] * fa;
      sta.izfrc[katm] += fb1 * crabbc[2] - fb2 * crbccd[2];
      sta.izfrc[matm] += -fc1 * crbccd[2] + fc2 * crabbc[2];
      sta.izfrc[patm] += -fd * crbccd[2];
    }

    // Scale forces by the time step and commit them to "global"
    for (i = 0; i < natom; i++) {
      sta.ixfrc[i] *= gms->slowFrcMult;
      sta.iyfrc[i] *= gms->slowFrcMult;
      sta.izfrc[i] *= gms->slowFrcMult;
      gms->frcX.HostData[i + atmstart] = sta.ixfrc[i];
      gms->frcY.HostData[i + atmstart] = sta.iyfrc[i];
      gms->frcZ.HostData[i + atmstart] = sta.izfrc[i];
    }

    // Commit energies to "global"
    if (stepidx % gms->ntpr == 0) {
      nrgidx = (h * gms->nsgmdout) + (stepidx / gms->ntpr);
      gms->Uelec.HostData[nrgidx] = uelec;
      gms->Usolv.HostData[nrgidx] = usolv;
      gms->Uvdw.HostData[nrgidx]  = uvdw;
      gms->Udihe.HostData[nrgidx] = udihe;
    }
  }
}

//-----------------------------------------------------------------------------
// WriteMultiRestart: write checkpoint files for each system in the array.
//
// Arguments:
//   gms:       information for all systems, including coordinates and
//              topologies
//   tpbank:    bank of topologies for all systems
//   ppctrl:    peptide control data, contains the restart file names and
//              start times
//   tj:        trajectory control information, contains the file suffix
//   currstep:  the current step number
//-----------------------------------------------------------------------------
static void WriteMultiRestart(gpuMultiSim *gms, prmtop *tpbank, pepcon *ppctrl,
			      trajcon *tj, long long int currstep)
{
  int h, i, j, natom;
  coord tc;
  trajcon mockTJ;

  // Make a coordinates struct large enough to hold the largest system
  // (the number of atoms will be adjusted on the fly)
  j = 0;
  for (i = 0; i < ppctrl->nsys; i++) {
    natom = gms->atomReadLimits.HostData[i].y -
            gms->atomReadLimits.HostData[i].x;
    if (j < natom) {
      j = natom;
    }
  }
  tc = CreateCoord(j);

  // Create a mockup of the trajcon for managing the standard writer
  mockTJ.mode = 9;
  mockTJ.rstbase = CreateCmat(1, MAXNAME);
  mockTJ.rstsuff = CreateCmat(1, GRID);

  // Download coordinates, velocities, and forces from the GPU
#ifdef CUDA
  DownloadGpuFloat(&gms->crdrX, NULL);
  DownloadGpuFloat(&gms->crdrY, NULL);
  DownloadGpuFloat(&gms->crdrZ, NULL);
  DownloadGpuFloat(&gms->velrX, NULL);
  DownloadGpuFloat(&gms->velrY, NULL);
  DownloadGpuFloat(&gms->velrZ, NULL);
#endif

  // Loop over all systems
  for (h = 0; h < gms->nsys; h++) {

    // Put coordinates and velocities into the coordinate struct
    tc.natom = tpbank[h].natom;
    j = 0;
    for (i = gms->atomReadLimits.HostData[h].x;
	 i < gms->atomReadLimits.HostData[h].y; i++) {
      tc.loc[3*j    ] = gms->crdrX.HostData[i];
      tc.loc[3*j + 1] = gms->crdrY.HostData[i];
      tc.loc[3*j + 2] = gms->crdrZ.HostData[i];
      tc.vel[3*j    ] = gms->velrX.HostData[i];
      tc.vel[3*j + 1] = gms->velrY.HostData[i];
      tc.vel[3*j + 2] = gms->velrZ.HostData[i];
      j++;
    }

    // Velocities do not need to be re-wound as a snapshot was
    // taken at the correct moment in the code.
    mockTJ.currstep = 0;
    mockTJ.nfistep = 0;
    mockTJ.nstep = 1;
    mockTJ.currtime = ppctrl->starttime[h] + (tj->dt * (double)currstep);
    mockTJ.OverwriteOutput = tj->OverwriteOutput;
    strcpy(mockTJ.rstbase.map[0], ppctrl->rstrtbases.map[h]);
    strcpy(mockTJ.rstsuff.map[0], tj->rstsuff.map[0]);
    WriteRst(NULL, &tc, &tpbank[h], &mockTJ, 0, 0);
  }

  // Free allocated memory
  DestroyCmat(&mockTJ.rstbase);
  DestroyCmat(&mockTJ.rstsuff);
  DestroyCoord(&tc);
}

//-----------------------------------------------------------------------------
// SystemRattleC: positional corrections for RATTLE constraints in one of many
//                implicit solvent MD simulations.
//
// Arguments:
//   gms:       information for all systems, including coordinates and
//              topologies
//   sta:       System Temporary Allocations, to hold data on any one system's
//              coordinates, velocities, and forces.  This temporary scratch
//              space is pre-allocated once for the entirety of the simulation.
//   sysID:     the system for which sta has data (can't loop over all systems
//-----------------------------------------------------------------------------
static void SystemRattleC(gpuMultiSim *gms, slmem sta, int sysID)
{
  int i, j, niter, atma, atmb, violations, atmstart, atmend, natom;
  float delta, rma, rmb, term, xterm, yterm, zterm;
  double rx, ry, rz, rrefx, rrefy, rrefz, r2, l0, dot, xmove, ymove, zmove;

  // Treat coordinates as their initial positions (fp32) plus a fixed
  // precision pertubation.  The precision offered by a 32-bit integer
  // storing numbers ranging between -2 and +2 is enough to converge
  // RATTLE positions to very small tolerances.  Nine 32-bit arrays are
  // needed for this calculation: three floating-point arrays to hold
  // each atom's current X, Y, and Z coordinates, three more to hold
  // each atom's previous coordinates, and three integer arrays to hold
  // fixed-precision local adjustments.  At this stage, forces are not
  // needed (there's our three integer arrays) and, specifically on the
  // GPU, velocities can be stashed in registers to make room for the
  // previous coordinates.
  atmstart = gms->atomReadLimits.HostData[sysID].x;
  atmend   = gms->atomReadLimits.HostData[sysID].y;
  natom = atmend - atmstart;

  // Obtain positions from "global" memory and zero the perturbations
  j = 0;
  for (i = atmstart; i < atmend; i++) {
    sta.xcrd[j] = gms->atomX.HostData[i];
    sta.ycrd[j] = gms->atomY.HostData[i];
    sta.zcrd[j] = gms->atomZ.HostData[i];
    sta.ixfrc[j] = 0;
    sta.iyfrc[j] = 0;
    sta.izfrc[j] = 0;
    j++;
  }

  // Each warp will execute up to niter times, and all
  // RATTLE'd bonds in the warp will get processed for
  // the same number of iterations pending the slowest
  // one to converge.
  for (i = gms->cnstReadLimits.HostData[sysID].x;
       i < gms->cnstReadLimits.HostData[sysID].y; i += GRID) {
    niter = 0;
    violations = 1;
    while (niter < gms->maxRattleIter && violations > 0) {
      violations = 0;

      // Zero the elements of the buffer that matter
      SetIVec(sta.active, natom, 0);
      for (j = 0; j < GRID; j++) {
        atma = gms->cnstInfo.HostData[i + j].x;
        if (atma == 0xffffffff) {
          continue;
        }
        atmb = (atma & 0xffff);
        atma >>= 16;
        sta.ixbuff[atma] = 0;
        sta.iybuff[atma] = 0;
        sta.izbuff[atma] = 0;
        sta.ixbuff[atmb] = 0;
        sta.iybuff[atmb] = 0;
        sta.izbuff[atmb] = 0;
	sta.active[atma] = 1;
	sta.active[atmb] = 1;
      }
      for (j = 0; j < GRID; j++) {

        // Get the constraint parameters
        atma = gms->cnstInfo.HostData[i + j].x;
        if (atma == 0xffffffff) {
          continue;
        }
        atmb = (atma & 0xffff);
        atma >>= 16;
        l0 = (double)(gms->cnstInfo.HostData[i + j].y) / 1.0e8;
        l0 *= l0;
        rx =  (double)sta.xcrd[atmb] + (double)sta.ixfrc[atmb]*FPSCALEicn -
             ((double)sta.xcrd[atma] + (double)sta.ixfrc[atma]*FPSCALEicn);
        ry =  (double)sta.ycrd[atmb] + (double)sta.iyfrc[atmb]*FPSCALEicn -
             ((double)sta.ycrd[atma] + (double)sta.iyfrc[atma]*FPSCALEicn);
        rz =  (double)sta.zcrd[atmb] + (double)sta.izfrc[atmb]*FPSCALEicn -
             ((double)sta.zcrd[atma] + (double)sta.izfrc[atma]*FPSCALEicn);
        r2 = rx*rx + ry*ry + rz*rz;
        delta = (float)(l0 - r2);
        if (fabs(delta) > gms->rattleTol) {
          violations++;
          rrefx = (double)sta.xprvcrd[atmb] - (double)sta.xprvcrd[atma];
          rrefy = (double)sta.yprvcrd[atmb] - (double)sta.yprvcrd[atma];
          rrefz = (double)sta.zprvcrd[atmb] - (double)sta.zprvcrd[atma];
          dot = rx*rrefx + ry*rrefy + rz*rrefz;
          rma = gms->atomHDTM.HostData[atmstart + atma] * gms->hdtm2invm;
          rmb = gms->atomHDTM.HostData[atmstart + atmb] * gms->hdtm2invm;
          term = (float)1.2 * delta / ((float)2.0 * (float)dot * (rma + rmb));
          xterm = rrefx * term;
          yterm = rrefy * term;
          zterm = rrefz * term;
          sta.ixbuff[atma] -= (int)(xterm * rma * FPSCALEcn);
          sta.iybuff[atma] -= (int)(yterm * rma * FPSCALEcn);
          sta.izbuff[atma] -= (int)(zterm * rma * FPSCALEcn);
          sta.ixbuff[atmb] += (int)(xterm * rmb * FPSCALEcn);
          sta.iybuff[atmb] += (int)(yterm * rmb * FPSCALEcn);
          sta.izbuff[atmb] += (int)(zterm * rmb * FPSCALEcn);
        }
      }
      for (j = 0; j < natom; j++) {
	if (sta.active[j] == 1) {
	  sta.ixfrc[j] += sta.ixbuff[j];
	  sta.iyfrc[j] += sta.iybuff[j];
	  sta.izfrc[j] += sta.izbuff[j];
	}
      }
      niter++;
    }
    if (niter == gms->maxRattleIter && violations > 0) {
      printf("SystemRattleC >> Warning. %d iterations of RATTLE were unable "
             "to resolve\nSystemRattleC >> geometry in system %d.\n",
             gms->maxRattleIter, sysID);
    }
  }

  // Contribute the perturbations back to global
  for (i = 0; i < natom; i++) {
    xmove = (double)(sta.ixfrc[i]) * FPSCALEicn;
    ymove = (double)(sta.iyfrc[i]) * FPSCALEicn;
    zmove = (double)(sta.izfrc[i]) * FPSCALEicn;
    gms->atomX.HostData[atmstart + i] = sta.xcrd[i] + xmove;
    gms->atomY.HostData[atmstart + i] = sta.ycrd[i] + ymove;
    gms->atomZ.HostData[atmstart + i] = sta.zcrd[i] + zmove;
    gms->velcX.HostData[atmstart + i] += xmove * gms->invdtVF;
    gms->velcY.HostData[atmstart + i] += ymove * gms->invdtVF;
    gms->velcZ.HostData[atmstart + i] += zmove * gms->invdtVF;
  }
}

//-----------------------------------------------------------------------------
// MultiVelVerletC: the coordinate update for velocity Verlet integration,
//                  with and without Langevin thermostating.
//
// Arguments:
//   gms:       information for all systems, including coordinates and
//              topologies
//   sta:       System Temporary Allocations, to hold data on any one system's
//              coordinates, velocities, and forces.  This temporary scratch
//              space is pre-allocated once for the entirety of the simulation.
//-----------------------------------------------------------------------------
static void MultiVelVerletC(gpuMultiSim *gms, slmem sta)
{
  int h, i, atmstart, atmend;
  long int counter;
  float rsd, hmdt, dimass;
  const float c_explic = (gms->Tstat.active == 1) ? gms->Tstat.c_explic : 0.0;

  // Move atoms according to forces, with Langevin integration if requested
  if (gms->Tstat.active == 1) {
    counter = gms->rndcon;
  }
  const float dt = sqrt(418.4) * gms->dt;
  const float hdt = 0.5 * dt / FPSCALEfrc;
  for (h = 0; h < gms->nsys; h++) {
    atmstart = gms->atomReadLimits.HostData[h].x;
    atmend   = gms->atomReadLimits.HostData[h].y;

    // Store positions for geometric constraints
    if (gms->rattle == 1) {
      for (i = atmstart; i < atmend; i++) {
        sta.xprvcrd[i - atmstart] = gms->atomX.HostData[i];
        sta.yprvcrd[i - atmstart] = gms->atomY.HostData[i];
        sta.zprvcrd[i - atmstart] = gms->atomZ.HostData[i];
      }
    }

    // Atom movement
    if (gms->Tstat.active == 1) {
      for (i = atmstart; i < atmend; i++) {
        hmdt = gms->atomHDTM.HostData[i];
        gms->velcX.HostData[i] = (gms->velcX.HostData[i] * c_explic) +
                                 hmdt * (float)(gms->frcX.HostData[i] +
                                                gms->langX.HostData[i]);
        gms->velcY.HostData[i] = (gms->velcY.HostData[i] * c_explic) +
                                 hmdt * (float)(gms->frcY.HostData[i] +
                                                gms->langY.HostData[i]);
        gms->velcZ.HostData[i] = (gms->velcZ.HostData[i] * c_explic) +
                                 hmdt * (float)(gms->frcZ.HostData[i] +
                                                gms->langZ.HostData[i]);
        gms->atomX.HostData[i] += dt * gms->velcX.HostData[i];
        gms->atomY.HostData[i] += dt * gms->velcY.HostData[i];
        gms->atomZ.HostData[i] += dt * gms->velcZ.HostData[i];
      }
    }
    else {
      for (i = atmstart; i < atmend; i++) {
        hmdt = gms->atomHDTM.HostData[i];
        gms->velcX.HostData[i] += hmdt * (float)(gms->frcX.HostData[i]);
        gms->velcY.HostData[i] += hmdt * (float)(gms->frcY.HostData[i]);
        gms->velcZ.HostData[i] += hmdt * (float)(gms->frcZ.HostData[i]);
        gms->atomX.HostData[i] += dt * gms->velcX.HostData[i];
        gms->atomY.HostData[i] += dt * gms->velcY.HostData[i];
        gms->atomZ.HostData[i] += dt * gms->velcZ.HostData[i];
      }
    }

    // Apply geometric constraints
    if (gms->rattle == 1) {
      SystemRattleC(gms, sta, h);
    }
  }
  if (gms->Tstat.active == 1) {
    gms->rndcon = counter;
  }
}

//-----------------------------------------------------------------------------
// SystemRattleV: velocity corrections for RATTLE constraints in one of many
//                implicit solvent MD simulations.
//
// Arguments:
//   gms:       information for all systems, including coordinates and
//              topologies
//   sta:       System Temporary Allocations, to hold data on any one system's
//              coordinates, velocities, and forces.  This temporary scratch
//              space is pre-allocated once for the entirety of the simulation.
//-----------------------------------------------------------------------------
static void SystemRattleV(gpuMultiSim *gms, slmem sta, int sysID)
{
  int i, j, niter, atma, atmb, violations, atmstart, atmend, natom;
  float rma, rmb, term, xterm, yterm, zterm;
  double rx, ry, rz, vx, vy, vz, r2, l0, dot;

  atmstart = gms->atomReadLimits.HostData[sysID].x;
  atmend   = gms->atomReadLimits.HostData[sysID].y;
  natom    = atmend - atmstart;

  // Obtain positions from "global" memory.  Forces are
  // already stored in global memory, so no need to stash
  // them in the slmem struct.  The integer force arrays
  // will again serve to hold the perturbations.
  j = 0;
  for (i = atmstart; i < atmend; i++) {
    sta.xcrd[j] = gms->atomX.HostData[i];
    sta.ycrd[j] = gms->atomY.HostData[i];
    sta.zcrd[j] = gms->atomZ.HostData[i];
    sta.xprvcrd[j] = gms->velcX.HostData[i];
    sta.yprvcrd[j] = gms->velcY.HostData[i];
    sta.zprvcrd[j] = gms->velcZ.HostData[i];
    sta.ixfrc[j] = 0;
    sta.iyfrc[j] = 0;
    sta.izfrc[j] = 0;
    j++;
  }
  const float rtoldt = gms->rattleTol / (sqrt(418.4) * gms->dt);
  for (i = gms->cnstReadLimits.HostData[sysID].x;
       i < gms->cnstReadLimits.HostData[sysID].y; i += GRID) {
    niter = 0;
    violations = 1;
    while (niter < gms->maxRattleIter && violations > 0) {
      violations = 0;
      for (j = 0; j < GRID; j++) {
        atma = gms->cnstInfo.HostData[i + j].x;
        if (atma == 0xffffffff) {
          continue;
        }
        atmb = (atma & 0xffff);
        atma >>= 16;
        rx = (double)sta.xcrd[atmb] - (double)sta.xcrd[atma];
        ry = (double)sta.ycrd[atmb] - (double)sta.ycrd[atma];
        rz = (double)sta.zcrd[atmb] - (double)sta.zcrd[atma];
        vx =  (double)sta.xprvcrd[atmb] + (double)sta.ixfrc[atmb]*FPSCALEicnv -
             ((double)sta.xprvcrd[atma] + (double)sta.ixfrc[atma]*FPSCALEicnv);
        vy =  (double)sta.yprvcrd[atmb] + (double)sta.iyfrc[atmb]*FPSCALEicnv -
             ((double)sta.yprvcrd[atma] + (double)sta.iyfrc[atma]*FPSCALEicnv);
        vz =  (double)sta.zprvcrd[atmb] + (double)sta.izfrc[atmb]*FPSCALEicnv -
             ((double)sta.zprvcrd[atma] + (double)sta.izfrc[atma]*FPSCALEicnv);
        dot = rx*vx + ry*vy + rz*vz;
        rma = gms->atomHDTM.HostData[atmstart + atma] * gms->hdtm2invm;
        rmb = gms->atomHDTM.HostData[atmstart + atmb] * gms->hdtm2invm;
        l0  = (double)(gms->cnstInfo.HostData[i + j].y) / 1.0e8;
        term = (float)-1.2 * (float)dot / ((rma + rmb) * l0 * l0);
        if (fabs(term) > rtoldt) {
          violations++;
          xterm = rx * term;
          yterm = ry * term;
          zterm = rz * term;
          sta.ixfrc[atma] -= (int)(xterm * rma * FPSCALEcnv);
          sta.iyfrc[atma] -= (int)(yterm * rma * FPSCALEcnv);
          sta.izfrc[atma] -= (int)(zterm * rma * FPSCALEcnv);
          sta.ixfrc[atmb] += (int)(xterm * rmb * FPSCALEcnv);
          sta.iyfrc[atmb] += (int)(yterm * rmb * FPSCALEcnv);
          sta.izfrc[atmb] += (int)(zterm * rmb * FPSCALEcnv);
        }
      }
      niter++;
    }
    if (niter == gms->maxRattleIter && violations > 0) {
      printf("SystemRattleV >> Warning. %d iterations of RATTLE were unable "
             "to resolve\nSystemRattleV >> velocities in system %d.\n",
             gms->maxRattleIter, sysID);
    }
  }

  // Contribute the perturbations back to the velocity arrays
  for (i = 0; i < natom; i++) {
    gms->velcX.HostData[atmstart + i] = sta.xprvcrd[i] +
                                        (double)(sta.ixfrc[i]) * FPSCALEicnv;
    gms->velcY.HostData[atmstart + i] = sta.yprvcrd[i] +
                                        (double)(sta.iyfrc[i]) * FPSCALEicnv;
    gms->velcZ.HostData[atmstart + i] = sta.zprvcrd[i] +
                                        (double)(sta.izfrc[i]) * FPSCALEicnv;
  }
}

//-----------------------------------------------------------------------------
// MultiVelVerletV: the velocity update for velocity Verlet integration,
//                  with and without Langevin thermostating.
//
// Arguments:
//   gms:       information for all systems, including coordinates and
//              topologies
//   sta:       System Temporary Allocations, to hold data on any one system's
//              coordinates, velocities, and forces.  This temporary scratch
//              space is pre-allocated once for the entirety of the simulation.
//-----------------------------------------------------------------------------
static void MultiVelVerletV(gpuMultiSim *gms, slmem sta)
{
  int h, i, atmstart, atmend, lgX, lgY, lgZ;
  long int counter;
  float rsd, hmdt, dimass, hsdfac;
  const float c_implic = (gms->Tstat.active == 1) ? gms->Tstat.c_implic : 0.0;

  if (gms->Tstat.active == 1) {
    counter = gms->rndcon;
  }
  const float hdt = 0.5 * sqrt(418.4) * gms->dt / FPSCALEfrc;
  for (h = 0; h < gms->nsys; h++) {
    atmstart = gms->atomReadLimits.HostData[h].x;
    atmend   = gms->atomReadLimits.HostData[h].y;
    if (gms->Tstat.active == 1) {
      hsdfac = gms->Tstat.sdfac.HostData[h];
      for (i = atmstart; i < atmend; i++) {
        dimass = gms->atomMass.HostData[i];
        hmdt = gms->atomHDTM.HostData[i];
        rsd = hsdfac * sqrt(dimass);
        lgX = (int)(rsd * GaussBoxMuller(&counter) * FPSCALEfrc);
        lgY = (int)(rsd * GaussBoxMuller(&counter) * FPSCALEfrc);
        lgZ = (int)(rsd * GaussBoxMuller(&counter) * FPSCALEfrc);
        gms->langX.HostData[i] = lgX;
        gms->langY.HostData[i] = lgY;
        gms->langZ.HostData[i] = lgZ;
        gms->velcX.HostData[i] = (hmdt * (float)(gms->frcX.HostData[i] + lgX) +
                                  gms->velcX.HostData[i]) * c_implic;
        gms->velcY.HostData[i] = (hmdt * (float)(gms->frcY.HostData[i] + lgY) +
                                  gms->velcY.HostData[i]) * c_implic;
        gms->velcZ.HostData[i] = (hmdt * (float)(gms->frcZ.HostData[i] + lgZ) +
                                  gms->velcZ.HostData[i]) * c_implic;
      }
    }
    else {
      for (i = atmstart; i < atmend; i++) {
        hmdt = gms->atomHDTM.HostData[i];
        gms->velcX.HostData[i] += hmdt * (float)(gms->frcX.HostData[i]);
        gms->velcY.HostData[i] += hmdt * (float)(gms->frcY.HostData[i]);
        gms->velcZ.HostData[i] += hmdt * (float)(gms->frcZ.HostData[i]);
      }
    }

    // Velocity constraints
    if (gms->rattle == 1) {
      SystemRattleV(gms, sta, h);
    }
  }
  if (gms->Tstat.active == 1) {
    gms->rndcon = counter;
  }
}

//-----------------------------------------------------------------------------
// MultiKineticEnergy: compute the kinetic energy of multiple systems in
//                     flight.
//
// Arguments:
//   gms:       information for all systems, including coordinates and
//              topologies
//   cyclecon:  counter to track the position in the multiple time step cycle
//   stepidx:   step count within the segment (not total step count)
//-----------------------------------------------------------------------------
static void MultiKineticEnergy(gpuMultiSim *gms, int cyclecon, int stepidx)
{
  int h, i, atmstart, atmend, nrgidx;
  float result, vx, vy, vz;

  // Only do this on the first cycle of MTS
  if (cyclecon > 0) {
    return;
  }

  // Once again, nested loops over all systems, then all atoms
  for (h = 0; h < gms->nsys; h++) {
    atmstart = gms->atomReadLimits.HostData[h].x;
    atmend   = gms->atomReadLimits.HostData[h].y;
    result = (float)0.0;
    for (i = atmstart; i < atmend; i++) {
      vx = gms->velcX.HostData[i];
      vy = gms->velcY.HostData[i];
      vz = gms->velcZ.HostData[i];
      result += gms->atomMass.HostData[i] * (vx*vx + vy*vy + vz*vz);
    }
    if (gms->Tstat.active == 1) {
      result *= (float)1.0 + ((float)0.5 * gms->Tstat.gamma_ln *
                              gms->dt);
    }
    result *= (float)0.5;
    if (stepidx % gms->ntpr == 0) {
      nrgidx = (h * gms->nsgmdout) + (stepidx / gms->ntpr);
      gms->Ukine.HostData[nrgidx] = result;
      gms->Temp.HostData[nrgidx] = 2.0 * result * gms->invNDF.HostData[h];
    }
  }
}

//-----------------------------------------------------------------------------
// RecenterSystems: re-center all systems by translating their coordinates
//                  such that they become geometry-centered on the origin.
//                  This is helpful for keeping coordinates from becoming too
//                  large that precision is lost.  This routine uses fixed
//                  precision to maintain agreement between CPU and GPU code.
//
// Arguments:
//   gms:    information for all systems, including coordinates and topologies
//-----------------------------------------------------------------------------
static void RecenterSystemsDP(gpuMultiSim *gms)
{
  int h, i, atmstart, atmend;
  double np, geomcx, geomcy, geomcz;

  // Once again, nested loops over all systems, then all atoms
  for (h = 0; h < gms->nsys; h++) {
    atmstart = gms->atomReadLimits.HostData[h].x;
    atmend   = gms->atomReadLimits.HostData[h].y;
    geomcx = 0.0;
    geomcy = 0.0;
    geomcz = 0.0;
    for (i = atmstart; i < atmend; i++) {
      geomcx += gms->atomX.HostData[i];
      geomcy += gms->atomY.HostData[i];
      geomcz += gms->atomZ.HostData[i];
    }
    np = (double)(atmend - atmstart);
    geomcx /= np;
    geomcy /= np;
    geomcz /= np;
    for (i = atmstart; i < atmend; i++) {
      gms->atomX.HostData[i] -= geomcx;
      gms->atomY.HostData[i] -= geomcy;
      gms->atomZ.HostData[i] -= geomcz;
    }
  }
}

//-----------------------------------------------------------------------------
// RecenterSystems: re-center all systems by translating their coordinates
//                  such that they become geometry-centered on the origin.
//                  This is helpful for keeping coordinates from becoming too
//                  large that precision is lost.  This routine uses fixed
//                  precision to maintain agreement between CPU and GPU code.
//
// Arguments:
//   gms:    information for all systems, including coordinates and topologies
//-----------------------------------------------------------------------------
static void RecenterSystems(gpuMultiSim *gms)
{
  int h, i, atmstart, atmend, igeomcx, igeomcy, igeomcz;
  float np, geomcx, geomcy, geomcz;

  // Once again, nested loops over all systems, then all atoms
  for (h = 0; h < gms->nsys; h++) {
    atmstart = gms->atomReadLimits.HostData[h].x;
    atmend   = gms->atomReadLimits.HostData[h].y;
    igeomcx = 0;
    igeomcy = 0;
    igeomcz = 0;
    for (i = atmstart; i < atmend; i++) {
      igeomcx += (int)(gms->atomX.HostData[i] * RCSCALE);
      igeomcy += (int)(gms->atomY.HostData[i] * RCSCALE);
      igeomcz += (int)(gms->atomZ.HostData[i] * RCSCALE);
    }
    np = (float)(atmend - atmstart);
    geomcx = (float)igeomcx / (RCSCALE * np);
    geomcy = (float)igeomcy / (RCSCALE * np);
    geomcz = (float)igeomcz / (RCSCALE * np);
    for (i = atmstart; i < atmend; i++) {
      gms->atomX.HostData[i] -= geomcx;
      gms->atomY.HostData[i] -= geomcy;
      gms->atomZ.HostData[i] -= geomcz;
    }
  }
}

//-----------------------------------------------------------------------------
// MultiVelocityInit: initialize velocities for multiple systems according to
//                    a Maxwell distribution.
//
// Arguments:
//   gms:    information for all systems, including coordinates and topologies
//   tj:     trajectory control information (needed for irest)
//   ppctrl: peptide control data, used for storing system-specific notes on
//           initialization
//   sta:    scratch space for CPU-side computations on all systems.  If the
//           GPU is in use, then this struct will only be used for staging
//           velocities.
//   devspc: GPU specifications
//-----------------------------------------------------------------------------
static void MultiVelocityInit(gpuMultiSim *gms, trajcon *tj, pepcon *ppctrl,
			      slmem sta, gpuSpecs *devspc)
{
  int h, i, atmstart, atmend, hasvel, natom, rstore, tstore;
  long int counter;
  float efac;
#ifdef CUDA
  gpuMultiSim gmsCopy;
#endif

  // If the &cntrl namelist calls for fresh velocities, create
  // them.  If there are no initial velocities, provide some.
  counter = gms->rndcon;
  for (h = 0; h < gms->nsys; h++) {
    atmstart = gms->atomReadLimits.HostData[h].x;
    atmend   = gms->atomReadLimits.HostData[h].y;
    natom = atmend - atmstart;
    hasvel = 0;
    for (i = 0; i < natom; i++) {
      if (fabs(gms->velcX.HostData[i]) > (float)1.0e-4 ||
	  fabs(gms->velcY.HostData[i]) > (float)1.0e-4 ||
	  fabs(gms->velcZ.HostData[i]) > (float)1.0e-4) {
	hasvel = 1;
	break;
      }
    }
    if (hasvel == 0 || tj->irest == 0) {
      const float ebeta = (ppctrl->T[h] > 1.0e-8) ?
                           sqrt(GASCNST * ppctrl->T[h]) : 0.0;
      for (i = atmstart; i < atmend; i++) {
        efac = ebeta * sqrt(1.0 / gms->atomMass.HostData[i]);
        gms->velcX.HostData[i] = efac * GaussBoxMuller(&counter);
        gms->velcY.HostData[i] = efac * GaussBoxMuller(&counter);
        gms->velcZ.HostData[i] = efac * GaussBoxMuller(&counter);
      }
    }
  }
  gms->rndcon = counter;

  // Recenter the coordinates as they are first found
  // in the restart file using double precision.
  RecenterSystemsDP(gms);

  // Apply constraints to initial positions, simply
  // to overcome roundoff error in the restart file.
  if (gms->rattle == 1) {
    for (h = 0; h < gms->nsys; h++) {
      atmstart = gms->atomReadLimits.HostData[h].x;
      atmend   = gms->atomReadLimits.HostData[h].y;
      for (i = atmstart; i < atmend; i++) {
        sta.xprvcrd[i - atmstart] = gms->atomX.HostData[i];
        sta.yprvcrd[i - atmstart] = gms->atomY.HostData[i];
        sta.zprvcrd[i - atmstart] = gms->atomZ.HostData[i];
      }
      SystemRattleC(gms, sta, h);
    }
  }

  // Initialize the forces.  In order to co-opt the GPU dynamics routine
  // for computing forces of just the initial configuration, make a copy
  // of the simulation array and set all initial forces to zero.  Rewind
  // the positions half a time step (without constraining the resulting
  // particle positions) with the CPU routine so that the GPU kernel
  // can apply its own half step update to return to the correct initial
  // positions and then start computing forces and updating as usual.
#ifdef CUDA
  gmsCopy = CopyGms(gms, 0);
  gmsCopy.nsgmdout = 1;
  gmsCopy.ntpr = 1;
  gmsCopy.Tstat.active = 0;
  gmsCopy.dt = -(gms->dt);
  gmsCopy.rattle = 0;
  for (h = 0; h < gms->nsys; h++) {
    for (i = gms->atomReadLimits.HostData[h].x;
         i < gms->atomReadLimits.HostData[h].y; i++) {
      gmsCopy.frcX.HostData[i] = 0;
      gmsCopy.frcY.HostData[i] = 0;
      gmsCopy.frcZ.HostData[i] = 0;
    }
  }
  MultiVelVerletC(&gmsCopy, sta);
  gmsCopy.dt = gms->dt;
  gmsCopy.rattle = gms->rattle;
  UploadGmsAll(&gmsCopy);
  LaunchDynamics(&gmsCopy, ppctrl->blockDim, ppctrl->blocks, devspc);
  DownloadGpuFloat(&gmsCopy.atomX, gms->atomX.HostData);
  DownloadGpuFloat(&gmsCopy.atomY, gms->atomY.HostData);
  DownloadGpuFloat(&gmsCopy.atomZ, gms->atomZ.HostData);
  DownloadGpuFloat(&gmsCopy.velcX, gms->velcX.HostData);
  DownloadGpuFloat(&gmsCopy.velcY, gms->velcY.HostData);
  DownloadGpuFloat(&gmsCopy.velcZ, gms->velcZ.HostData);
  DownloadGpuInt(&gmsCopy.frcX, gms->frcX.HostData);
  DownloadGpuInt(&gmsCopy.frcY, gms->frcY.HostData);
  DownloadGpuInt(&gmsCopy.frcZ, gms->frcZ.HostData);
  DestroyGpuMultiSim(&gmsCopy);
#else
  // Rewind the coordinates for re-centering in line with the GPU
  gms->dt = -gms->dt;
  rstore = gms->rattle;
  gms->rattle = 0;
  tstore = gms->Tstat.active;
  gms->Tstat.active = 0;
  MultiVelVerletC(gms, sta);
  RecenterSystems(gms);
  gms->dt = -gms->dt;
  MultiVelVerletC(gms, sta);
  gms->rattle = rstore;
  gms->Tstat.active = tstore;

  // Compute the complete force.
  GmsSlowStepEnergyAndForces(gms, sta, 0);
  GmsFastStepEnergyAndForces(gms, sta, 0, 0);

  // Perform the first velocity half-kick on the CPU.  If
  // the GPU is in use, these results will get uploaded.
  MultiVelVerletV(gms, sta);

  // Perform subsequent intermediate steps on the CPU,
  // as this part is not computationally intensive.
  for (i = 1; i < gms->slowFrcMult; i++) {

    // Velocity Verlet, update positions
    MultiVelVerletC(gms, sta);

    // High-frequency force calculation
    GmsFastStepEnergyAndForces(gms, sta, 1, 0);

    // Velocity Verlet, update velocities
    MultiVelVerletV(gms, sta);
  }
#endif
}

//-----------------------------------------------------------------------------
// SnapshotPhaseSpace: take a snapshot of the positions and velocities in all
//                     systems for the purpose of writing a restart file.
//                     Immediately after a restart file is read, forces will be
//                     computed (both fast and slow components if MTS is
//                     active), velocities updated, and constraints applied.
//                     The snapshot therefore needs to happen right before
//                     computation of the slow forces on the system.
//
// Arguments:
//   gms:      information for all systems, including coordinates and
//             topologies
//-----------------------------------------------------------------------------
static void SnapshotPhaseSpace(gpuMultiSim *gms)
{
  int h, i;

  for (h = 0; h < gms->nsys; h++) {
    for (i = gms->atomReadLimits.HostData[h].x;
	 i < gms->atomReadLimits.HostData[h].y; i++) {
      gms->crdrX.HostData[i] = gms->atomX.HostData[i];
      gms->crdrY.HostData[i] = gms->atomY.HostData[i];
      gms->crdrZ.HostData[i] = gms->atomZ.HostData[i];
      gms->velrX.HostData[i] = gms->velcX.HostData[i];
      gms->velrY.HostData[i] = gms->velcY.HostData[i];
      gms->velrZ.HostData[i] = gms->velcZ.HostData[i];
    }
  }
}

//-----------------------------------------------------------------------------
// RunMultiDynamics: host-side code to run MD simulations of non-periodic
//                   systems.  This will use the many of the same data
//                   structures and concepts as the GPU device code to
//                   facilitate checking, prototyping, and learning.  This
//                   routine mimics Dynamics in the Integrator.c library for
//                   periodic systems in its Velocity Verlet implementation.
//
// Arguments:
//   gms:      information for all systems, including coordinates and
//             topologies
//   sta:      scratch space for force computations
//   stepidx:  step count within the segment (not total step count)
//-----------------------------------------------------------------------------
static void RunMultiDynamics(gpuMultiSim *gms, slmem sta, int stepidx)
{
  int i;

  // Re-center all systems. (Do not remove net momentum, as this likely
  // reflects Langevin thermostating and any residual momentum would be
  // replaced almost immediately.)
  if (stepidx % gms->ntpr == 0) {
    RecenterSystems(gms);
  }

  // Multiple time steps (all systems are done by each function inside)
  for (i = 0; i < gms->slowFrcMult; i++) {

    // Velocity Verlet, update positions
    MultiVelVerletC(gms, sta);

    // Take velocity snapshot for restart file writing
    if (i == 0 && stepidx == gms->ntwx - 1) {
      SnapshotPhaseSpace(gms);
    }

    // Compute long time-step forces if needed.  Always compute
    // short time-step forces, and compute intermediate forces
    // if appropriate.
    if (i == 0) {
      GmsSlowStepEnergyAndForces(gms, sta, stepidx);
      GmsFastStepEnergyAndForces(gms, sta, 0, stepidx);
    }
    else {
      GmsFastStepEnergyAndForces(gms, sta, 1, stepidx);
    }

    // Velocity Verlet, update velocities
    MultiVelVerletV(gms, sta);

    // Kinetic energy computation
    MultiKineticEnergy(gms, i, stepidx);
  }
}

//-----------------------------------------------------------------------------
// WriteMultiDiagnostics: write diagnostic output files for each system.  All
//                        files have already been initialized.
//
// Arguments:
//   ppctrl:  peptide simulations control data (file names, directives for
//            certain enhanced sampling methods)
//   tpbank:  the bank of topologies
//   tj:      trajectory control data (contains &cntrl namelist information)
//   gms:     information for all systems, including coordinates and
//            topologies
//   sysUV:   energy tracking for each system (holds running averages)
//   segidx:  segment index
//-----------------------------------------------------------------------------
static void WriteMultiDiagnostics(pepcon *ppctrl, prmtop* tpbank, trajcon *tj,
                                  gpuMultiSim *gms, Energy* sysUV, int segidx)
{
  int h, i, nrgidx;
  double uelec, usolv, uvdw, ubond, uangl, udihe, ukine;
  char fname[MAXNAME];
  FILE *outp;
  trajcon mockTJ;

  // Loop over all systems, then over all output iterations
  for (h = 0; h < ppctrl->nsys; h++) {
    SpliceFileName(tj, ppctrl->outbases.map[h], tj->outsuff, fname, 0);
    outp = fopen(fname, "a");
    for (i = 0; i < ppctrl->nsgmdout; i++) {

      // Prepare mock trajcon and Energy structs
      // to hold the information from each step
      mockTJ.currstep = (segidx * gms->ntwx) + (i * gms->ntpr);
      mockTJ.currtime = ppctrl->starttime[h] +
                        (mockTJ.currstep * gms->dt * (double)gms->slowFrcMult);
      nrgidx = (h * gms->nsgmdout) + i;
      sysUV[h].T = gms->Temp.HostData[nrgidx];
      uelec = gms->Uelec.HostData[nrgidx];
      usolv = gms->Usolv.HostData[nrgidx];
      uvdw  = gms->Uvdw.HostData[nrgidx];
      ubond = gms->Ubond.HostData[nrgidx];
      uangl = gms->Uangl.HostData[nrgidx];
      udihe = gms->Udihe.HostData[nrgidx];
      ukine = gms->Ukine.HostData[nrgidx];
      sysUV[h].kine  = ukine;
      sysUV[h].eptot = uelec + usolv + uvdw + ubond + uangl + udihe;
      sysUV[h].etot  = sysUV[h].eptot + ukine;
      sysUV[h].bond  = ubond;
      sysUV[h].angl  = uangl;
      sysUV[h].dihe  = udihe;
      sysUV[h].elec  = uelec;
      sysUV[h].vdw12 = uvdw;
      sysUV[h].relec = usolv;
      AccAverageEnergy(&sysUV[h]);

      // Direct the information to the mdout
      fprintf(outp, "\\\\\\\n");
      PrintStateInfo(&sysUV[h], &mockTJ, &tpbank[h], outp, 1);
    }
    fclose(outp);
  }
}

//-----------------------------------------------------------------------------
// WriteMultiTrajectory: write trajectories for each of the systems being
//                       simulated.
//
// Arguments:
//   ppctrl:  peptide simulations control data (file names, directives for
//            certain enhanced sampling methods)
//   tpbank:  bank of topologies for all systems
//   tj:      trajectory control data (contains &cntrl namelist information)
//   gms:     information for all systems, including coordinates and
//            topologies
//   segidx:  segment index
//   irstrt:  flag to signal that this is writing restart files
//-----------------------------------------------------------------------------
static void WriteMultiTrajectory(pepcon *ppctrl, prmtop* tpbank, trajcon *tj,
                                 gpuMultiSim *gms, int segidx, int irstrt)
{
  int i, j, k, natom;
  coord tc;
  trajcon mockTJ;

  // Make a coordinates struct large enough to hold the largest system
  // (the number of atoms will be adjusted on the fly)
  j = 0;
  for (i = 0; i < ppctrl->nsys; i++) {
    natom = gms->atomReadLimits.HostData[i].y -
            gms->atomReadLimits.HostData[i].x;
    if (j < natom) {
      j = natom;
    }
  }
  tc = CreateCoord(j);

  // Set trajectory details
  mockTJ.currstep = (segidx + 1) * tj->ntwx;
  mockTJ.nfistep = 0;
  mockTJ.nstep = tj->nstep;
  mockTJ.irest = tj->irest;
  mockTJ.ntpr = tj->ntpr;
  mockTJ.ntwx = tj->ntwx;
  mockTJ.trjbase = CreateCmat(1, MAXNAME);
  mockTJ.trjsuff = CreateCmat(1, GRID);
  mockTJ.OverwriteOutput = tj->OverwriteOutput;

  // Loop over all systems
  for (i = 0; i < ppctrl->nsys; i++) {

    // Set file names
    strcpy(mockTJ.trjbase.map[0], ppctrl->trjbases.map[i]);
    strcpy(mockTJ.trjsuff.map[0], tj->trjsuff.map[0]);

    // Build the coordinates from the gpuMultiSim data
    tc.natom = gms->atomReadLimits.HostData[i].y -
               gms->atomReadLimits.HostData[i].x;
    k = 0;
    for (j = gms->atomReadLimits.HostData[i].x;
         j < gms->atomReadLimits.HostData[i].y; j++) {
      tc.loc[3*k    ] = gms->atomX.HostData[j];
      tc.loc[3*k + 1] = gms->atomY.HostData[j];
      tc.loc[3*k + 2] = gms->atomZ.HostData[j];
      k++;
    }
    if (irstrt == 1) {
      k = 0;
      for (j = gms->atomReadLimits.HostData[i].x;
           j < gms->atomReadLimits.HostData[i].y; j++) {
        tc.vel[3*k    ] = gms->velcX.HostData[j];
        tc.vel[3*k + 1] = gms->velcY.HostData[j];
        tc.vel[3*k + 2] = gms->velcZ.HostData[j];
        k++;
      }
    }
    for (j = 0; j < 3; j++) {
      tc.gdim[j    ] = 0.0;
      tc.gdim[j + 3] = 90.0;
    }
    if (irstrt == 0) {
      if (tj->ioutfm == 0) {
        WriteCrd(NULL, &tc, 1, &mockTJ, &tpbank[i], 0);
      }
      else if (tj->ioutfm == 1) {
        WriteCDF(NULL, &tc, 1, &mockTJ, &ppctrl->CDFHandles[i], &tpbank[i], 0);
      }
    }
  }

  // Free allocated memory
  DestroyCmat(&mockTJ.trjbase);
  DestroyCmat(&mockTJ.trjsuff);
  DestroyCoord(&tc);
}

#ifdef CUDA
//-----------------------------------------------------------------------------
// PlanRandomUsage: plan the generation and usage of random numbers by the
//                  GPU.
//
// Arguments:
//   ppctrl:   peptide simulations control data (file names, directives for
//             certain enhanced sampling methods)
//   tpbank:   bank of topologies used by all systems
//   tj:       trajctory control information, contains ntpr for determining
//             the random number refresh rate
//-----------------------------------------------------------------------------
static void PlanRandomUsage(pepcon *ppctrl, prmtop *tpbank, trajcon *tj)
{
  int i, j, nprng, maxRefresh;

  // Count the total number of atoms in all systems
  nprng = 0;
  for (i = 0; i < ppctrl->nsys; i++) {
    nprng += 3 * tpbank[i].natom;
    nprng = ((nprng + GRID_BITS_MASK) / GRID) * GRID;
  }

  // The maximum number of randoms that will be allocated is dependent on the
  // amount of GPU RAM, or the maximum size of a single array.
  if (nprng >= 3 * 2097152) {
    printf("PlanRandomUsage >> There are too many distinct systems that would "
	   "require\nPlanRandomUsage >> unique random number sequences.  "
           "Separate this group of\nPlanRandomUsage >> simulations into "
	   "smaller batches.\n");
    exit(1);
  }

  // Random numbers will be refreshed in each system at a rate determined
  // by the size of the random array that can be allocated.
  ppctrl->prngRefresh = 25165824 / nprng;
  maxRefresh = ppctrl->prngRefresh;
  i = 1;
  for (j = 2; j <= 7; j++) {
    while (tj->ntpr % (i * j) == 0 && i * j < maxRefresh) {
      i *= j;
    }
  }
  ppctrl->prngRefresh = i;
  i = 1;
  for (j = 7; j >= 2; j--) {
    while (tj->ntpr % (i * j) == 0 && i * j < maxRefresh) {
      i *= j;
    }
  }
  if (i > ppctrl->prngRefresh) {
    ppctrl->prngRefresh = i;
  }
  ppctrl->prngRefresh *= ppctrl->nBondSteps;

  // Determine the PRNG heap allocation
  i = (ppctrl->prngRefresh + 1) * nprng;
  i = ((i + GRID_BITS_MASK) / GRID) * GRID;
  ppctrl->totalPRNG = i;
}

//-----------------------------------------------------------------------------
// PlanGpuUtilization: as the name implies, this will look at characteristics
//                     of the selected GPU and the simulations to be completed,
//                     then decide how best to run everything.
//
// Arguments:
//   ppctrl:   peptide simulations control data (file names, directives for
//             certain enhanced sampling methods)
//   tpbank:   bank of topologies used by all systems
//   tj:       trajectory control data, contains RATTLE setting
//   devspc:   decription of the GPU hardware
//-----------------------------------------------------------------------------
static void PlanGpuUtilization(pepcon *ppctrl, prmtop *tpbank, trajcon *tj,
                               gpuSpecs *devspc)
{
  int i, maxatom, msrounds;
  double xlOcc, lgOcc, mdOcc, smOcc;

  // Determine whether all systems fit within a certain size limit
  maxatom = 0;
  for (i = 0; i < ppctrl->nsys; i++) {
    if (tpbank[i].natom > maxatom) {
      maxatom = tpbank[i].natom;
    }
  }

  // Test that no system violates the maximum size limit
  if (maxatom > LG_ATOM_COUNT) {
    printf("mdgx >> One or more systems contains more than the limit of %d "
           "atoms.\n", LG_ATOM_COUNT);
    for (i = 0; i < ppctrl->nsys; i++) {
      if (tpbank[i].natom > LG_ATOM_COUNT) {
        printf("mdgx >>  %-50.50s :: %d atoms\n",
               ppctrl->tpinames.map[i], tpbank[i].natom);
      }
    }
    exit(1);
  }

  // If the user has specified a block size, try to accommodate it
  if (ppctrl->blockreq[0] != '\0') {
    if (strcmp(ppctrl->blockreq, "SMALL") == 0) {
      if (maxatom <= SM_ATOM_COUNT) {
#if defined(AMBER_PLATFORM_AMD)
    ppctrl->blockDim = 128;
    ppctrl->blocks = 8 * devspc->MPcount;
#else
	if (tj->rattle == 0) {
	  ppctrl->blockDim = (devspc->maxThrPerMP == 1024) ? 256 : 320;
	}
	else {
	  ppctrl->blockDim = (devspc->maxThrPerMP == 1024) ? 256 : 288;
	}
    ppctrl->blocks = 4 * devspc->MPcount;
#endif
	return;
      }
      else {
	printf("PlanGpuUtilization >> Warning.  Unable to accommodate request "
	       "for small blocks.\nPlanGpuUtilization >> Maximum atom size %d "
	       "> limit of %d.", maxatom, SM_ATOM_COUNT);
      }
    }
    else if (strcmp(ppctrl->blockreq, "MEDIUM") == 0) {
      if (maxatom <= MD_ATOM_COUNT) {
#if defined(AMBER_PLATFORM_AMD)
    ppctrl->blockDim = 256;
    ppctrl->blocks = 4 * devspc->MPcount;
#else
	if (tj->rattle == 0) {
	  ppctrl->blockDim = (devspc->maxThrPerMP == 1024) ? 512 : 640;
	}
	else {
	  ppctrl->blockDim = (devspc->maxThrPerMP == 1024) ? 512 : 576;
	}
    ppctrl->blocks = 2 * devspc->MPcount;
#endif
	return;
      }
      else {
	printf("PlanGpuUtilization >> Warning.  Unable to accommodate request "
	       "for medium-sized\nPlanGpuUtilization >> blocks.  Maximum atom "
	       "size %d > limit of %d.", maxatom, MD_ATOM_COUNT);
      }
    }
    else if (strcmp(ppctrl->blockreq, "LARGE") == 0) {
#if defined(AMBER_PLATFORM_AMD)
    ppctrl->blockDim = 512;
    ppctrl->blocks = 2 * devspc->MPcount;
#else
    ppctrl->blockDim = 1024;
    ppctrl->blocks = devspc->MPcount;
#endif
      return;
    }
  }

  // Seek the strategy that will best keep the GPU occupied.  Assume
  // Assume that a system will occupy smaller blocks slightly better
  // than larger ones, provided that the system does fit in either.
  msrounds = (ppctrl->nsys + devspc->MPcount - 1) / devspc->MPcount;
  lgOcc = (double)ppctrl->nsys / (double)(msrounds * devspc->MPcount);
  lgOcc *= 1.0 - (0.04 * LG_ATOM_COUNT / maxatom);
  if (maxatom < MD_ATOM_COUNT) {
    msrounds = (ppctrl->nsys + (2 * devspc->MPcount) - 1) /
               (2 * devspc->MPcount);
    mdOcc = (double)ppctrl->nsys / (double)(msrounds * 2 * devspc->MPcount);
  }
  else {
    mdOcc = 0.0;
  }
  mdOcc *= 1.0 - (0.04 * MD_ATOM_COUNT / maxatom);
  if (maxatom < SM_ATOM_COUNT) {
    msrounds = (ppctrl->nsys + (4 * devspc->MPcount) - 1) /
               (4 * devspc->MPcount);
    smOcc = (double)ppctrl->nsys / (double)(msrounds * 4 * devspc->MPcount);
  }
  else {
    smOcc = 0.0;
  }
  smOcc *= 1.0 - (0.04 * SM_ATOM_COUNT / maxatom);

  // Select the best utilization
  if (lgOcc > mdOcc && lgOcc > smOcc) {
#if defined(AMBER_PLATFORM_AMD)
    ppctrl->blockDim = 512;
    ppctrl->blocks = 2 * devspc->MPcount;
#else
    ppctrl->blockDim = 1024;
    ppctrl->blocks = devspc->MPcount;
#endif
  }
  else if (mdOcc > lgOcc && mdOcc > smOcc) {
#if defined(AMBER_PLATFORM_AMD)
    ppctrl->blockDim = 256;
    ppctrl->blocks = 4 * devspc->MPcount;
#else
    if (tj->rattle == 0) {
      ppctrl->blockDim = (devspc->maxThrPerMP == 1024) ? 512 : 640;
    }
    else {
      ppctrl->blockDim = (devspc->maxThrPerMP == 1024) ? 512 : 576;
    }
    ppctrl->blocks = 2 * devspc->MPcount;
#endif
  }
  else if (smOcc > lgOcc && smOcc > mdOcc) {
#if defined(AMBER_PLATFORM_AMD)
    ppctrl->blockDim = 128;
    ppctrl->blocks = 8 * devspc->MPcount;
#else
    if (tj->rattle == 0) {
      ppctrl->blockDim = (devspc->maxThrPerMP == 1024) ? 256 : 320;
    }
    else {
      ppctrl->blockDim = (devspc->maxThrPerMP == 1024) ? 256 : 288;
    }
    ppctrl->blocks = 4 * devspc->MPcount;
#endif
  }
}

//-----------------------------------------------------------------------------
// SortSystemSize: compare the sizes of two systems for quicksort to arrange
//                 the keys.
//
// Arguments:
//-----------------------------------------------------------------------------
static int SortSystemSize(const void *sysA, const void *sysB)
{
  int nA = ((int2*)sysA)[0].x;
  int nB = ((int2*)sysB)[0].x;

  if (nA < nB) {
    return 1;
  }
  else if (nA > nB) {
    return -1;
  }
  else {
    return 0;
  }
}

//-----------------------------------------------------------------------------
// OrganizeSystems: arrange systems in order of decreasing overall atom count,
//                  in a bid to put the most expensive simulations first on
//                  the GPU's docket.  Small simulations will backfill idle
//                  multiprocessors as the larger ones finish.
//
// Arguments:
//   ppctrl:   peptide simulations control data (contains information about
//             the simulations array layout on the GPU)
//   tpbank:   bank of topologies used by all systems
//   crdbank:  bank of coordinates for all systems
//-----------------------------------------------------------------------------
static void OrganizeSystems(pepcon *ppctrl, prmtop* tpbank, coord* crdbank)
{
  int i, j, nmol;
  double* THold;
  double* TmixHold;
  double* PmixHold;
  double* timeHold;
  prmtop* tpHold;
  coord* tcHold;
  cmat nameHold;

  // If there are more blocks than systems to run, don't change the order
  nmol = ppctrl->nsys;
  if (ppctrl->blocks >= nmol) {
    return;
  }

  // Make an array of all system sizes
  int2* sysatoms;
  sysatoms = (int2*)malloc(nmol * sizeof(int2));
  for (i = 0; i < nmol; i++) {
    sysatoms[i].x = tpbank[i].natom;
    sysatoms[i].y = i;
  }

  // Sort the keys
  qsort(sysatoms, nmol, sizeof(int2), SortSystemSize);

  // Re-arrange the list: topologies, coordinates, all file names,
  // starting times, plus temperature and Hamiltonian mixing values.
  tpHold = (prmtop*)malloc(nmol * sizeof(prmtop));
  tcHold = (coord*)malloc(nmol * sizeof(coord));
  nameHold = CreateCmat(6 * nmol, MAXNAME);
  THold = (double*)malloc(nmol * sizeof(double));
  TmixHold = (double*)malloc(nmol * sizeof(double));
  PmixHold = (double*)malloc(nmol * sizeof(double));
  timeHold = (double*)malloc(nmol * sizeof(double));
  for (i = 0; i < nmol; i++) {
    tpHold[i] = tpbank[i];
    tcHold[i] = crdbank[i];
    strcpy(nameHold.map[i         ], ppctrl->tpinames.map[i]);
    strcpy(nameHold.map[i +   nmol], ppctrl->tpfnames.map[i]);
    strcpy(nameHold.map[i + 2*nmol], ppctrl->crdnames.map[i]);
    strcpy(nameHold.map[i + 3*nmol], ppctrl->outbases.map[i]);
    strcpy(nameHold.map[i + 4*nmol], ppctrl->trjbases.map[i]);
    strcpy(nameHold.map[i + 5*nmol], ppctrl->rstrtbases.map[i]);
    THold[i] = ppctrl->T[i];
    TmixHold[i] = ppctrl->Tmix[i];
    PmixHold[i] = ppctrl->Pmix[i];
    timeHold[i] = ppctrl->starttime[i];
  }
  for (i = 0; i < nmol; i++) {
    tpbank[i] = tpHold[sysatoms[i].y];
    crdbank[i] = tcHold[sysatoms[i].y];
    strcpy(ppctrl->tpinames.map[i],   nameHold.map[sysatoms[i].y         ]);
    strcpy(ppctrl->tpfnames.map[i],   nameHold.map[sysatoms[i].y +   nmol]);
    strcpy(ppctrl->crdnames.map[i],   nameHold.map[sysatoms[i].y + 2*nmol]);
    strcpy(ppctrl->outbases.map[i],   nameHold.map[sysatoms[i].y + 3*nmol]);
    strcpy(ppctrl->trjbases.map[i],   nameHold.map[sysatoms[i].y + 4*nmol]);
    strcpy(ppctrl->rstrtbases.map[i], nameHold.map[sysatoms[i].y + 5*nmol]);
    ppctrl->T[i]         = THold[sysatoms[i].y];
    ppctrl->Tmix[i]      = TmixHold[sysatoms[i].y];
    ppctrl->Pmix[i]      = PmixHold[sysatoms[i].y];
    ppctrl->starttime[i] = timeHold[sysatoms[i].y];
  }

  // Free allocated memory
  free(tpHold);
  free(tcHold);
  free(THold);
  free(TmixHold);
  free(PmixHold);
  free(timeHold);
  DestroyCmat(&nameHold);
}

//-----------------------------------------------------------------------------
// PrintGpuDetails: adds an extra section to the mdout, just before results,
//                  to explain which GPU is in use and how it is being
//                  utilized.  This leaves the output file primed to print
//                  diagnostic state info.
//
// Arguments:
//   tj:       trajectory control data (contains &cntrl namelist information)
//   devspc:   device (GPU) specs
//   ppctrl:   peptide simulations control data (contains information about
//             the simulations array layout on the GPU)
//   sysid:    ID number of the system to rpint details for
//-----------------------------------------------------------------------------
static void PrintGpuDetails(trajcon *tj, gpuSpecs *devspc, pepcon *ppctrl,
                            int sysid)
{
  char fname[MAXNAME];
  FILE *outp;

  SpliceFileName(tj, ppctrl->outbases.map[sysid], tj->outsuff, fname, 1);
  outp = fopen(fname, "a");
  HorizontalRule(outp, 0);
  fprintf(outp, " (5.) GPU UTILIZATION\n\n");
  fprintf(outp, " - GPU Model: %s\n\n", devspc->name);
  fprintf(outp, " - GPU Architecture:           %d.%d\n", devspc->major,
          devspc->minor);
  fprintf(outp, " - Streaming Multiprocessors:  %d\n", devspc->MPcount);
  fprintf(outp, " - Max threads per SMP:        %d\n", devspc->maxThrPerMP);
  fprintf(outp, " - Simulation block size:      %d threads\n",
          ppctrl->blockDim);
  fprintf(outp, " - Block grid size:            %d\n", ppctrl->blocks);
  HorizontalRule(outp, 1);
  HorizontalRule(outp, 0);
  fprintf(outp, " (6.) RESULTS\n\n");
  fclose(outp);
}
#endif // CUDA-only functions

//-----------------------------------------------------------------------------
// PeptideSimulator: main function for running peptide simulations in implicit
//                   solvent.  If a valid &pptd namelist is present, mdgx will
//                   dive in here and then exit.
//
// Arguments:
//   ppctrl:   peptide simulations control data (file names, directives for
//             certain enhanced sampling methods)
//   tj:       trajectory control data (contains &cntrl namelist information)
//   etimers:  timings data
//-----------------------------------------------------------------------------
void PeptideSimulator(pepcon *ppctrl, trajcon *tj, char* epsrc,
                      execon *etimers)
{
  int i, j, k, nrgidx;
  long long int currstep;
  prmtop tpI, tpF;
  prmtop* tpbank;
  coord* crdbank;
  cellgrid mockCG;
  reccon mockRCINP;
  Energy* sysUV;
  slmem sta;
  gpuMultiSim gms;
  cmat chscr;
  gpuSpecs devspc;

  // Expand the systems so that each one has a single topology,
  // input coordinates file, trajectory, output, and restart
  // file to write.  Create hybrid topologies as needed based
  // on mixing factors.
  ExpandSystemIO(ppctrl, tj);

  // Get information about the GPU
#ifdef CUDA
  devspc = GetGpuSpecs(tj->Reckless);
#endif

  // Read in all system topologies and coordinates
  tpbank = (prmtop*)malloc(ppctrl->nsys * sizeof(prmtop));
  crdbank = (coord*)malloc(ppctrl->nsys * sizeof(coord));
  epsrc[0] = '\0';
  for (i = 0; i < ppctrl->nsys; i++) {

    // Get the topology, or a hybrid of two topologies
    if (ppctrl->Pmix[i] < 1.0e-8) {
      SetBasicTopology(&tpbank[i], ppctrl->tpinames.map[i], epsrc,
                       ppctrl->rattle, 0);
      GetPrmTop(&tpbank[i], tj, 1);
    }
    else {
      SetBasicTopology(&tpI, ppctrl->tpinames.map[i], epsrc,
                       ppctrl->rattle, 0);
      GetPrmTop(&tpI, tj, 1);
      SetBasicTopology(&tpF, ppctrl->tpfnames.map[i], epsrc,
                       ppctrl->rattle, 0);
      GetPrmTop(&tpF, tj, 1);
      tpbank[i] = InterpolateTopology(&tpI, &tpF, ppctrl->Pmix[i]);
      FreeTopology(&tpI);
      FreeTopology(&tpF);
    }

    // Get the coordinates
    crdbank[i] = ReadRst(&tpbank[i], ppctrl->crdnames.map[i],
                         &ppctrl->starttime[i]);
  }
#ifdef CUDA
  // Plan the attack
  PlanRandomUsage(ppctrl, tpbank, tj);
  PlanGpuUtilization(ppctrl, tpbank, tj, &devspc);

  // Re-arrange systems for best efficiency on the GPU
  OrganizeSystems(ppctrl, tpbank, crdbank);
#endif

  // Make non-bonded exclusion masks in CPU memory.  These will
  // start from the masks that other modules (parameter fitting,
  // configuration sampling) make, then compress the format
  // into something more streamlined for the GPU.
  InitGmsConstants(ppctrl, tj, &gms);
  InitGmsTopologies(ppctrl, tpbank, &gms);
  InitGmsCoordinates(ppctrl, tpbank, crdbank, &gms);
  sta = AllocateSlmem();
  MultiVelocityInit(&gms, tj, ppctrl, sta, &devspc);

  // If CUDA is compiled, stage this all on the GPU
#ifdef CUDA
  if (gms.Tstat.active == 1) {
    InitGpuPRNG(&gms, tj->igseed, ppctrl->blocks, ppctrl->blockDim);
  }
  UploadGmsAll(&gms);
#endif

  // Write the preambles to all mdout files
  tj->currstep = 0;
  chscr = CreateCmat(2, MAXNAME);
  strcpy(chscr.map[0], tj->ipcname.map[0]);
  strcpy(chscr.map[1], tj->outbase);
  for (i = 0; i < ppctrl->nsys; i++) {
    strcpy(tj->ipcname.map[0], ppctrl->crdnames.map[i]);
    strcpy(tj->outbase, ppctrl->outbases.map[i]);
    OpenDiagnosticsFile(tj, &tpbank[i], NULL, NULL, NULL, &crdbank[i], NULL,
                        0);
#ifdef CUDA
    PrintGpuDetails(tj, &devspc, ppctrl, i);
#endif
  }
  strcpy(tj->ipcname.map[0], chscr.map[0]);
  strcpy(tj->outbase, chscr.map[1]);

  // Initialize all trajectory files
  WriteMultiTrajectory(ppctrl, tpbank, tj, &gms, -1, 0);

  // Allocate to store energies and averages
  // that the output files can understand
  sysUV = (Energy*)malloc(ppctrl->nsys * sizeof(Energy));
  for (i = 0; i < ppctrl->nsys; i++) {
    InitializeEnergy(&sysUV[i], tj, &tpbank[i], 1);
    sysUV[i].P = 0.0;
    sysUV[i].V = 0.0;
    InitAverageEnergy(&sysUV[i]);
  }
  etimers->Write += mdgxStopTimer(etimers);

  // Dynamics loop
#ifdef CUDA
  for (i = 0; i < ppctrl->nsegment; i++) {

    // LaunchDynamics will cycle through "forces only" and
    // "forces with energies" kernels until the next coordinate
    // writing is at hand.
    LaunchDynamics(&gms, ppctrl->blockDim, ppctrl->blocks, &devspc);
    DownloadGpuFloat(&gms.atomX, NULL);
    DownloadGpuFloat(&gms.atomY, NULL);
    DownloadGpuFloat(&gms.atomZ, NULL);
    DownloadGpuDouble(&gms.Ubond, NULL);
    DownloadGpuDouble(&gms.Uangl, NULL);
    DownloadGpuDouble(&gms.Udihe, NULL);
    DownloadGpuDouble(&gms.Uelec, NULL);
    DownloadGpuDouble(&gms.Uvdw, NULL);
    DownloadGpuDouble(&gms.Usolv, NULL);
    DownloadGpuDouble(&gms.Ukine, NULL);

    // Compute temperature on the CPU side, based on reported kinetic energy
    for (j = 0; j < ppctrl->nsys; j++) {
      for (k = 0; k < gms.nsgmdout; k++) {
        nrgidx = (j * gms.nsgmdout) + k;
        gms.Temp.HostData[nrgidx] = 2.0 * gms.Ukine.HostData[nrgidx] *
                                    gms.invNDF.HostData[j];
      }
    }

    // Write mdout files
    WriteMultiDiagnostics(ppctrl, tpbank, tj, &gms, sysUV, i);

    // Write trajectory files
    if (ppctrl->notraj == 0) {
      WriteMultiTrajectory(ppctrl, tpbank, tj, &gms, i, 0);
    }

    // Write restart files with snapshotted data
    currstep = ((long long int)(i + 1) * (long long int)gms.ntpr *
                (long long int)gms.nsgmdout);
    if (currstep != 0 && tj->ntwr != 0) {
      if (currstep % tj->ntwr == 0) {
        WriteMultiRestart(&gms, tpbank, ppctrl, tj, currstep);
      }
    }
    etimers->Write += mdgxStopTimer(etimers);
  }
#else
  for (i = 0; i < ppctrl->nsegment; i++) {

    // Propagate dynamics
    for (j = 0; j < gms.ntwx; j++) {
      RunMultiDynamics(&gms, sta, j);
    }

    // Write mdout files
    WriteMultiDiagnostics(ppctrl, tpbank, tj, &gms, sysUV, i);

    // Write trajectory files
    WriteMultiTrajectory(ppctrl, tpbank, tj, &gms, i, 0);

    // Write restart files with snapshotted data
    currstep = ((long long int)(i + 1) * (long long int)gms.ntpr *
                (long long int)gms.nsgmdout);
    if (currstep != 0 && tj->ntwr != 0) {
      if (currstep % tj->ntwr == 0) {
        WriteMultiRestart(&gms, tpbank, ppctrl, tj, currstep);
      }
    }
    etimers->Write += mdgxStopTimer(etimers);
  }
#endif

  // Close all diagnostics files
  strcpy(chscr.map[1], tj->outbase);
  mockRCINP.nlev = 1;
#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mockCG.tid);
  mockCG.dspcomm = MPI_COMM_WORLD;
#else
  mockCG.tid = 0;
#endif
  for (i = 0; i < ppctrl->nsys; i++) {
    strcpy(tj->outbase, ppctrl->outbases.map[i]);
    CloseDiagnosticsFile(tj, &tpbank[i], &mockRCINP, &sysUV[i], etimers,
                         &mockCG, 1);
  }
  strcpy(tj->ipcname.map[0], chscr.map[0]);
  DestroyCmat(&chscr);
  etimers->Write += mdgxStopTimer(etimers);

  // Free allocated memory
  DestroyGpuMultiSim(&gms);
  for (i = 0; i < ppctrl->nsys; i++) {
    DestroyEnergyTracker(&sysUV[i]);
  }
  free(sysUV);
  DestroySlmem(&sta);

  // Exit
  exit(0);
}
