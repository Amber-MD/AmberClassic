#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <dirent.h>
#include <time.h>
#include "ConfigSamp.h"
#include "CrdManip.h"
#include "Topology.h"
#include "VirtualSites.h"
#include "Trajectory.h"
#include "Matrix.h"
#include "Manual.h"
#include "Parse.h"
#include "IPolQ.h"
#include "ParamFit.h"
#include "mdgxVector.h"
#include "Bonded.h"
#include "Macros.h"
#include "Random.h"

//-----------------------------------------------------------------------------
// InitConfigs: initialize a configuration sampling input data structure with
//              basic settings for file names and energy minimization
//              conditions.  This includes the qmset sub-structure QMsettings
//              (see ConfigSampDS.h for descriptions of both objects).
//
// Arguments:
//   cfsinp:    pointer to an uninitialized, unallocated configs struct
//-----------------------------------------------------------------------------
void InitConfigs(configs *cfsinp)
{
  cfsinp->nops = 0;
  cfsinp->ncombo = 0;
  cfsinp->ncyc = 100;
  cfsinp->maxcyc = 1000;
  cfsinp->nreopt = 0;
  cfsinp->nbelly = 0;
  cfsinp->atomlimit = 512;
  cfsinp->fconv = 1.0e-5;
  cfsinp->step0 = 0.01;
  cfsinp->stepconv = 1.0e-8;
  cfsinp->reshuffle = 0;
  cfsinp->Ereplace = 1.0e-4;
  cfsinp->proximity = 5.0;
  cfsinp->verbose = 0;
  cfsinp->freezeH = 0;
  cfsinp->showorigins = -1;
  cfsinp->rattle = 0;
  cfsinp->MaxRattleIter = 100;
  cfsinp->FreezePolarH = 0;
  cfsinp->PreventHBonds = 0;
  cfsinp->HBMargin = 2.5;
  cfsinp->SpinHydroxyls = 0;
  cfsinp->KSpinOH = 64.0;
  cfsinp->fbhwSpinOH = 60.0 * PI / 180.0;
  cfsinp->rattletol = 1.0e-6;
  cfsinp->ops = (cfigop*)malloc(32*sizeof(cfigop));
  cfsinp->combo = (coupler*)malloc(32*sizeof(coupler));
  cfsinp->outbase = CreateCmat(1, MAXNAME);
  cfsinp->outsuff = CreateCmat(1, MAXNAME);
  cfsinp->ostyle = CreateCmat(1, 32);
  cfsinp->belly = CreateCmat(1, MAXNAME);
  cfsinp->maxbstrn = 10.0;
  cfsinp->maxastrn = 10.0;
  cfsinp->strainlim = 200.0;
  cfsinp->rmsdtol = -1.0;
  cfsinp->simEtol = -1.0;
  cfsinp->QMsettings.ncpu = 1;
  cfsinp->QMsettings.spin = 1;
  cfsinp->QMsettings.MaxCore = 512;
  cfsinp->QMsettings.checkpoint = (char*)malloc(MAXNAME*sizeof(char));
  cfsinp->QMsettings.opt = 0;
  cfsinp->QMsettings.xferRst = 1;
  strcpy(cfsinp->outbase.map[0], "conf");
  strcpy(cfsinp->outsuff.map[0], "pdb");
  strcpy(cfsinp->ostyle.map[0], "PDB");
  sprintf(cfsinp->QMsettings.theory, "mp2");
  sprintf(cfsinp->QMsettings.basis, "cc-pVTZ");
  sprintf(cfsinp->shfdir, "down");
  sprintf(cfsinp->shuffletype, "jackknife");
  cfsinp->QMsettings.checkpoint[0] = '\0';
}

//-----------------------------------------------------------------------------
// DestroyConfigs: free the memory allocated for a configs data struct,
//                 including the QMsettings sub-structure.
//
// Arguments:
//   cfsinp:    pointer to an uninitialized, unallocated configs struct
//-----------------------------------------------------------------------------
void DestroyConfigs(configs *cfsinp)
{
  int i;

  for (i = 0; i < cfsinp->nops; i++) {
    DestroyCmat(&cfsinp->ops[i].atommasks);
    if (cfsinp->ops[i].order == 1) {
      DestroyDmat(&cfsinp->ops[i].refcrd);
    }
    free(cfsinp->ops[i].targets);
  }
  free(cfsinp->ops);
  free(cfsinp->combo);
  free(cfsinp->QMsettings.checkpoint);
  free(cfsinp->movable);
  DestroyCmat(&cfsinp->ostyle);
  DestroyCmat(&cfsinp->outbase);
  DestroyCmat(&cfsinp->outsuff);
  DestroyCmat(&cfsinp->origins);
  DestroyCmat(&cfsinp->belly);
}

//-----------------------------------------------------------------------------
// CreateExcltab: function for creating an exclusion table for a particular
//                atom.
//
// Arguments:
//   nexcl:    the number of exclusions to allocate
//-----------------------------------------------------------------------------
static excltab CreateExcltab(int nexcl)
{
  excltab extb;

  extb.nexcl = nexcl;
  extb.nblist = (int*)malloc(nexcl*sizeof(int));
  extb.ljval = (double*)malloc(nexcl*sizeof(double));
  extb.qqval = (double*)malloc(nexcl*sizeof(double));

  return extb;
}

//-----------------------------------------------------------------------------
// DestroyExcltab: function for freeing memory associated with an exclusions
//                 table.
//
// Arguments:
//   extb:     the exclusions table to free
//-----------------------------------------------------------------------------
static void DestroyExcltab(excltab *extb)
{
  free(extb->nblist);
  free(extb->ljval);
  free(extb->qqval);
}

//-----------------------------------------------------------------------------
// SortInteger: function for qsort to use when sorting integers into ascending
//              order.
//
// Arguments:
//   num[A,B]:   the integers to compare
//-----------------------------------------------------------------------------
int SortInteger(const void *numA, const void *numB)
{
  int intA = ((int*)numA)[0];
  int intB = ((int*)numB)[0];

  if (intA < intB) {
    return -1;
  }
  else if (intA > intB) {
    return 1;
  }
  else {
    return 0;
  }
}

//-----------------------------------------------------------------------------
// MakeExclList: a function for mapping all of the exclusions that each atom
//               of the topology will have.  This has to work reflexively for
//               both atoms on either side of a bond, angle, or dihedral.
//
// Arguments:
//   tp:      the system topology
//-----------------------------------------------------------------------------
static excltab* MakeExclList(prmtop *tp)
{
  int i, j, k, maxcount, found;
  int* xclcount;
  imat xclscratch;
  dmat qqvals;
  dmat ljvals;
  excltab* allex;

  // Allocate space to work in
  maxcount = 32;
  xclcount = (int*)calloc(tp->natom, sizeof(int));
  xclscratch = CreateImat(tp->natom, maxcount);
  qqvals = CreateDmat(tp->natom, maxcount, 0);
  ljvals = CreateDmat(tp->natom, maxcount, 0);
  allex = (excltab*)malloc(tp->natom*sizeof(excltab));

  // Count 1:2 and 1:3 exclusions
  for (i = 0; i < tp->natom; i++) {
    for (j = tp->ConExcl[i]; j < tp->ConExcl[i+1]; j++) {

      // Skip if there are no exclusions for this atom
      if (tp->ExclList[j] == -1) {
        continue;
      }

      // Catalog the 1:2 and 1:3 exclusions
      xclscratch.map[i][xclcount[i]] = tp->ExclList[j];
      qqvals.map[i][xclcount[i]] = 0.0;
      ljvals.map[i][xclcount[i]] = 0.0;
      xclscratch.map[tp->ExclList[j]][xclcount[tp->ExclList[j]]] = i;
      qqvals.map[tp->ExclList[j]][xclcount[tp->ExclList[j]]] = 0.0;
      ljvals.map[tp->ExclList[j]][xclcount[tp->ExclList[j]]] = 0.0;
      xclcount[i] += 1;
      xclcount[tp->ExclList[j]] += 1;
      if (xclcount[i] > maxcount || xclcount[tp->ExclList[j]] > maxcount) {
        xclscratch = ReallocImat(&xclscratch, tp->natom, maxcount);
        qqvals = ReallocDmat(&qqvals, tp->natom, maxcount);
        ljvals = ReallocDmat(&ljvals, tp->natom, maxcount);
      }
    }
  }

  // Count 1:4 exclusions
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->HLC[i].ndihe; j++) {
      if (tp->HLC[i].HC[j].eval14 == 1) {

        // Check to see whether this 1:4 is already present
        xclscratch.map[tp->HLC[i].HC[j].a][xclcount[tp->HLC[i].HC[j].a]] =
          tp->HLC[i].HC[j].d;
        qqvals.map[tp->HLC[i].HC[j].a][xclcount[tp->HLC[i].HC[j].a]] = 
          1.0 - tp->HLC[i].HC[j].scee;
        ljvals.map[tp->HLC[i].HC[j].a][xclcount[tp->HLC[i].HC[j].a]] = 
          1.0 - tp->HLC[i].HC[j].scnb;
        xclcount[tp->HLC[i].HC[j].a] += 1;
        xclscratch.map[tp->HLC[i].HC[j].d][xclcount[tp->HLC[i].HC[j].d]] =
          tp->HLC[i].HC[j].a;
        qqvals.map[tp->HLC[i].HC[j].d][xclcount[tp->HLC[i].HC[j].d]] =
          1.0 -tp->HLC[i].HC[j].scee;
        ljvals.map[tp->HLC[i].HC[j].d][xclcount[tp->HLC[i].HC[j].d]] =
          1.0 - tp->HLC[i].HC[j].scnb;
        xclcount[tp->HLC[i].HC[j].d] += 1;
        if (xclcount[tp->HLC[i].HC[j].a] >= maxcount || 
            xclcount[tp->HLC[i].HC[j].d] >= maxcount) {
          maxcount += 32;
          xclscratch = ReallocImat(&xclscratch, tp->natom, maxcount);
          qqvals = ReallocDmat(&qqvals, tp->natom, maxcount);
          ljvals = ReallocDmat(&ljvals, tp->natom, maxcount);
        }
      }
    }
  }

  // Allocate and fill each exclusion list
  for (i = 0; i < tp->natom; i++) {
    allex[i] = CreateExcltab(xclcount[i]);
    allex[i].idx = i;
    allex[i].llim = tp->natom;
    allex[i].hlim = 0;
    for (j = 0; j < xclcount[i]; j++) {
      allex[i].nblist[j] = xclscratch.map[i][j];
      allex[i].qqval[j] = qqvals.map[i][j];
      allex[i].ljval[j] = ljvals.map[i][j];
      if (xclscratch.map[i][j] < allex[i].llim) {
        allex[i].llim = xclscratch.map[i][j];
      }
      if (xclscratch.map[i][j] > allex[i].hlim) {
        allex[i].hlim = xclscratch.map[i][j];
      }
    }

    // Adjust the upper limit of the search to accomodate for loop behavior
    allex[i].hlim += 1;
  }

  // Free allocated memory
  DestroyImat(&xclscratch);
  DestroyDmat(&qqvals);
  DestroyDmat(&ljvals);
  free(xclcount);
  return allex;
}

//-----------------------------------------------------------------------------
// VectorRattle: fix the lengths of bonds between hydrogen atoms and
//               their parents to equilibrium values.
//
// Arguments:
//   tp:           topology for the system
//   cfsinp:       configuration sampling struct, containing copy counts and
//                 other minimization parameters
//   crdpop:       coordinates of the system in their new state (before bond
//                 length constraints have been applied)
//   crdsnapshot:  snapshot of coordinates prior to the move
//   nactive:      the numbe rof active conformations
//-----------------------------------------------------------------------------
static void VectorRattle(prmtop *tp, configs *cfsinp, dmat *crdpop,
                         dmat *crdsnapshot, int nactive)
{
  int i, j, k, ididx, nnode, niter, done, atma, atmb;
  int *myCC;
  int* atomids;
  double l0, delta, rx, ry, rz, rrefx, rrefy, rrefz;
  double r2, dot, term, xterm, yterm, zterm, rma, rmb;
  double rah[3];
  double* xref;
  double* invmass;
  double* xpos;
  
  // Allocate space to store reference coordinates
  xref    = (double*)malloc(3 * tp->RattleGrpMax * sizeof(double));
  xpos    = (double*)malloc(3 * tp->RattleGrpMax * sizeof(double));
  invmass = (double*)malloc(3 * tp->RattleGrpMax * sizeof(double));
  atomids = (int*)malloc(tp->RattleGrpMax * sizeof(int));

  // Shortcuts
  const double rtol = cfsinp->rattletol;
  const int maxiter = cfsinp->MaxRattleIter;
  
  // Loop over all atoms, find groups, iterate all instances of each group
  for (i = 0; i < tp->natom; i++) {

    // Does this atom control a RATTLE constraint group?
    if (tp->SHL[i].exe != 2) {
      continue;
    }
    
    // Get information about the constraint group
    myCC = tp->SHL[i].blist;
    ididx = 3*myCC[0] + 2;
    nnode = myCC[ididx - 1];
    for (j = 0; j < nnode; j++) {
      atomids[j] = myCC[ididx + j];
      invmass[j] = tp->InvMasses[atomids[j]];
    }
    
    // Loop over all copies (each will be done individually--there is not
    // much to be gained by trying to vectorize the iterative process of
    // converging the constraints).
    for (j = 0; j < nactive; j++) {  

      // Cache reference coordinates
      for (k = 0; k < nnode; k++) {
        xref[3*k    ] = crdsnapshot->map[3*atomids[k]    ][j];
        xref[3*k + 1] = crdsnapshot->map[3*atomids[k] + 1][j];
        xref[3*k + 2] = crdsnapshot->map[3*atomids[k] + 2][j];
        xpos[3*k    ] = crdpop->map[3*atomids[k]    ][j];
        xpos[3*k + 1] = crdpop->map[3*atomids[k] + 1][j];
        xpos[3*k + 2] = crdpop->map[3*atomids[k] + 2][j];
      }

      // Implement the constraints
      niter = 0;
      done = 0;
      while (niter < maxiter && done == 0) {
        done = 1;
        for (k = 0; k < myCC[0]; k++) {
          atma = myCC[3*k + 1];
          atmb = myCC[3*k + 2];
          rx = xpos[3*atmb    ] - xpos[3*atma    ];
          ry = xpos[3*atmb + 1] - xpos[3*atma + 1];
          rz = xpos[3*atmb + 2] - xpos[3*atma + 2];
          r2 = rx*rx + ry*ry + rz*rz;
          l0 = myCC[3*k + 3] / 1.0e8;
          l0 *= l0;
          delta = l0 - r2;
          if (fabs(delta) > rtol) {
            done = 0;
            rrefx = xref[3*atmb    ] - xref[3*atma    ];
            rrefy = xref[3*atmb + 1] - xref[3*atma + 1];
            rrefz = xref[3*atmb + 2] - xref[3*atma + 2];
            dot = rx*rrefx + ry*rrefy + rz*rrefz;
            term = 1.2 * delta / (2.0 * dot * (invmass[atma] + invmass[atmb]));
            xterm = rrefx * term;
            yterm = rrefy * term;
            zterm = rrefz * term;
            xpos[3*atma    ] -= xterm * invmass[atma];
            xpos[3*atma + 1] -= yterm * invmass[atma];
            xpos[3*atma + 2] -= zterm * invmass[atma];
            xpos[3*atmb    ] += xterm * invmass[atmb];
            xpos[3*atmb + 1] += yterm * invmass[atmb];
            xpos[3*atmb + 2] += zterm * invmass[atmb];
          }
        }
        niter++;
      }
      
      // Check for exceeding the maximum iteration count--this
      // produces a warning only, as strain may have induced this
      // problem and further minimization may resolve it.
      if (niter == maxiter && done == 0) {
        continue;
      }

      // Apply the constrained coordinates back to the active set
      for (k = 0; k < nnode; k++) {
        crdpop->map[3*atomids[k]    ][j] = xpos[3*k    ];
        crdpop->map[3*atomids[k] + 1][j] = xpos[3*k + 1];
        crdpop->map[3*atomids[k] + 2][j] = xpos[3*k + 2];
      }
    }
  }
  
  // Free allocated memory
  free(xref);
  free(xpos);
  free(invmass);
  free(atomids);
}

//-----------------------------------------------------------------------------
// VectorPlaceXpt: function for placing extra points onto all configurations
//                 using the frame information in the topology.
//
// Arguments:
//   tp:      the topology for the system
//   crdpop:  coordinates for all configurations
//   nactive:  the number of configurations actively being manipulated
//-----------------------------------------------------------------------------
static void VectorPlaceXpt(prmtop *tp, dmat *crdpop, int nactive)
{
  int i, j, k, epid, atma, atmb, atmc, atmd;
  double aloc[3], agloc[3], bloc[3], cloc[3], dloc[3], eploc[3], epgloc[3];
  expt *tmr;
  
  // Return immediately if there are no extra points
  if (tp->nxtrapt == 0) {
    return;
  }

  // Loop over each extra point, then over the coordinates.  Mock
  // arrays to hold the coordinates must be made to overcome the
  // rearrangement that vector processing of configurations entails.
  for (i = 0; i < tp->nxtrapt; i++) {

    // Information about the extra point
    tmr = &tp->xtrapts[i];
    epid = tmr->atomid;
    atma = tmr->fr1;
    atmb = tmr->fr2;
    atmc = (tmr->frstyle >= 2) ? tmr->fr3 : -1;
    atmd = (tmr->frstyle == 6 || tmr->frstyle == 9 || tmr->frstyle == 10) ?
           tmr->fr4 : -1;
    for (j = 0; j < nactive; j++) {

      // Copy atom locations into mock arrays
      for (k = 0; k < 3; k++) {        
        aloc[k] = crdpop->map[3*atma + k][j];
        bloc[k] = crdpop->map[3*atmb + k][j];
        agloc[k]  = 0.0;
        epgloc[k] = 0.0;
        if (atmc >= 0) {
          cloc[k] = crdpop->map[3*atmc + k][j];
        }
        if (atmd >= 0) {
          dloc[k] = crdpop->map[3*atmd + k][j];
        }
      }
      
      // Place the extra point
      XptLocator(aloc, agloc, bloc, cloc, dloc, eploc, epgloc, tmr);

      // Copy the extra point location back into the real coordinates
      for (k = 0; k < 3; k++) {
        crdpop->map[3*epid + k][j] = eploc[k];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// VectorCommuteXptForces: commute the forces from extra points onto their
//                         frame atoms for a set of configurations.
//
// Arguments:
//-----------------------------------------------------------------------------
static void VectorCommuteXptForces(prmtop *tp, dmat *crdpop, dmat *frcpop,
                                   int nactive)
{
  int i, j, k, atma, atmb, atmc, atmd, epid;
  double aloc[3], bloc[3], cloc[3], dloc[3], eploc[3];
  double afrc[3], bfrc[3], cfrc[3], dfrc[3], epfrc[3];
  expt *tmr;
  
  // Return immediately if there are no extra points
  if (tp->nxtrapt == 0) {
    return;
  }

  // Loop over each extra point, then over the coordinates.  Mock
  // arrays to hold the coordinates must be made to overcome the
  // rearrangement that vector processing of configurations entails.
  for (i = 0; i < tp->nxtrapt; i++) {

    // Information about the extra point
    tmr = &tp->xtrapts[i];
    epid = tmr->atomid;
    atma = tmr->fr1;
    atmb = tmr->fr2;
    atmc = (tmr->frstyle >= 2) ? tmr->fr3 : -1;
    atmd = (tmr->frstyle == 6 || tmr->frstyle == 9 || tmr->frstyle == 10) ?
      tmr->fr4 : -1;
    for (j = 0; j < nactive; j++) {
      
      // Copy atom locations and forces into mock arrays
      for (k = 0; k < 3; k++) {        
        aloc[k]  = crdpop->map[3*atma + k][j];
        bloc[k]  = crdpop->map[3*atmb + k][j];
        eploc[k] = crdpop->map[3*epid + k][j];
        afrc[k]  = 0.0;
        bfrc[k]  = 0.0;
        epfrc[k] = frcpop->map[3*epid + k][j];
        if (atmc >= 0) {
          cloc[k] = crdpop->map[3*atmc + k][j];
          cfrc[k] = 0.0;
        }
        if (atmd >= 0) {
          dloc[k] = crdpop->map[3*atmd + k][j];
          dfrc[k] = 0.0;
        }
      }
      
      // Commute the extra point's forces
      XptForceXfer(aloc, bloc, cloc, dloc, eploc, afrc, bfrc, cfrc, dfrc,
                   epfrc, tmr);

      // Copy the result back into the real forces.  Delete
      // forces on the extra point in the real array.
      for (k = 0; k < 3; k++) {
        frcpop->map[3*atma + k][j] += afrc[k];
        frcpop->map[3*atmb + k][j] += bfrc[k];
        if (atmc >= 0) {
          frcpop->map[3*atmc + k][j] += cfrc[k];
        }
        if (atmd >= 0) {
          frcpop->map[3*atmd + k][j] += dfrc[k];
        }
        frcpop->map[3*epid + k][j] = 0.0;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// InitCoordPopulation: initialize the "coordinate population," that is the
//                      set of initial structures from which all of the
//                      configurations will eventually be derived.  The total
//                      number of configurations may be specified by the user,
//                      but the places to find those initial coordinates are
//                      also user specified and certain inputs may cause mdgx
//                      to reassess the total number of configurations.  The
//                      user really only gets to say how many configurations to
//                      generate based on copies of a single template if the
//                      template is the only structure provided.  If the user
//                      supplies a list of files instead, or a trajectory with
//                      many frames, those will determine the ultimate number
//                      of configurations instead.
//
// Arguments:
//   cfsinp:  configuration sampling input data
//   tj:      trajectory control input data (the initial structure will be read
//            from here)
//   tp:      topology for the system (helps to determine whether a given file
//            is compatible as a coordinate set for an initial configuration)
//-----------------------------------------------------------------------------
dmat InitCoordPopulation(configs *cfsinp, trajcon *tj, prmtop *tp)
{
  int i, j, k, fitype, inptype;
  int ncoord, ndmat, nconf, hasbox, namelim;
  int* crdidx;
  int* dmatidx;
  int32_t flen;
  double lval, tstamp;
  double *dtmp, *dtm2p;
  FILE *finp;
  coord tc;
  coord* tcL;
  dmat crdpop, crdsnapshot;
  dmat *dmptr;
  dmat* cpopL;
  cmat crdnames, dmatnames, allnames;

  // First, a question: does the topology call for box information?
  hasbox = (tp->ifbox > 0) ? 1 : 0;

  // Next order of business: the user may have supplied a particular
  // coordinate file, in which case we will just read it in and make
  // as many copies of those coordinates as are needed, but the user
  // may alternatively have supplied a directory name, in which case
  // all files containing appropriate sets of coordinates in the
  // directory will be read in.  Or, the user may have supplied a
  // regular expression, in which case all files containing valid
  // set of coordinates matching the given template will be read in.
  inptype = FindFileType(tj->ipcname.map[0]);
  if (inptype == 0) {

    // Since we're just replicating one file, there is no need to wonder
    // where each new configuration came from.
    if (cfsinp->showorigins == -1) {
      cfsinp->showorigins = 0;
    }

    // First, check to see whether this is an input coordinates file
    // or perhaps a trajectory instead.
    fitype = ValidCoord(tj->ipcname.map[0]);
    if (fitype == 1 || fitype == 2 || fitype == 4) {

      // In this case, we will use the count given in the configs data
      // structure to know how many copies to make.  All configurations
      // will begin as copies of the one coordinate set.
      tstamp = 0.0;
      tc = ReadRst(tp, tj->ipcname.map[0], &tstamp);
      crdpop = CreateDmat(3*tp->natom, cfsinp->count, 0);
      for (i = 0; i < 3*tp->natom; i++) {
        lval = tc.loc[i];
        dtmp = crdpop.map[i];
        for (j = 0; j < cfsinp->count; j++) {
          dtmp[j] = lval;
        }
      }
      DestroyCoord(&tc);
    }
    else if (fitype == 3) {

      // In this case we will read the trajectory and take every frame,
      // ignoring whatever value the user might have given for the count.
      crdpop = ReadCrdTraj(tp, tj->ipcname.map[0], 0);
      cfsinp->count = crdpop.col;
    }
    else if (fitype == 5) {

      // Again, the value of count is reset by the
      // (NetCDF) trajectory frame count
      crdpop = ReadCDFTraj(tp, tj->ipcname.map[0], 0);
      cfsinp->count = crdpop.col;
    }

    // The origin of each configuration is the same
    cfsinp->origins = CreateCmat(cfsinp->count, MAXNAME);
    for (i = 0; i < cfsinp->count; i++) {
      strcpy(cfsinp->origins.map[i], tj->ipcname.map[0]);
    }
  }
  else if (inptype == 1 || inptype == 2) {

    // Since this will be a dump of all contents in the directory or
    // regular expression, we will want a clear description of where
    // each new configuration came from.
    if (cfsinp->showorigins == -1) {
      cfsinp->showorigins = 1;
    }

    // Open the directory and start parsing
    if (inptype == 1) {
      allnames = DirectoryFileSearch(tj->ipcname.map[0]);
    }

    // Parse the regular expression
    if (inptype == 2) {
      allnames = RegExpFileSearch(tj->ipcname.map[0]);
    }

    // Load up all coordinates
    tcL = (coord*)malloc(allnames.row*sizeof(coord));
    cpopL = (dmat*)malloc(allnames.row*sizeof(dmat));
    ncoord = 0;
    ndmat = 0;
    crdnames = CreateCmat(allnames.row, MAXNAME);
    dmatnames = CreateCmat(allnames.row, MAXNAME);
    for (i = 0; i < allnames.row; i++) {
      fitype = ValidCoord(allnames.map[i]);
      if (fitype == 1 || fitype == 2 || fitype == 4) {
        tstamp = 0.0;
        tcL[ncoord] = ReadRst(tp, allnames.map[i], &tstamp);
        strcpy(crdnames.map[ncoord], allnames.map[i]);
        ncoord++;
      }
      else if (fitype == 3 || fitype == 5) {
        if (fitype == 3) {
          cpopL[ndmat] = ReadCrdTraj(tp, allnames.map[i], 0);
        }
        else {
          cpopL[ndmat] = ReadCDFTraj(tp, allnames.map[i], 0);
        }
        strcpy(dmatnames.map[ndmat], allnames.map[i]);
        ndmat++;
      }
    }
    crdnames = ReallocCmat(&crdnames, ncoord, MAXNAME);
    dmatnames = ReallocCmat(&dmatnames, ndmat, MAXNAME);

    // Put some alphabetical ordering on all of the files, so that
    // a human can parse the results later.  Keep the individual
    // configurations at the front, but alphabetize those and the
    // trajectory files that follow.
    crdidx = AlphabetizeCmat(&crdnames);
    dmatidx = AlphabetizeCmat(&dmatnames);

    // Now collect all of the configurations that were found
    cfsinp->count = ncoord;
    for (i = 0; i < ndmat; i++) {
      cfsinp->count += cpopL[i].col;
    }
    cfsinp->origins = CreateCmat(cfsinp->count, MAXNAME);
    crdpop = CreateDmat(3*tp->natom, cfsinp->count, 0);
    for (i = 0; i < ncoord; i++) {
      for (j = 0; j < 3*tp->natom; j++) {
        crdpop.map[j][crdidx[i]] = tcL[i].loc[j];
      }
      strcpy(cfsinp->origins.map[crdidx[i]], crdnames.map[i]);
      DestroyCoord(&tcL[i]);
    }
    nconf = ncoord;
    for (i = 0; i < ndmat; i++) {
      for (j = 0; j < ndmat; j++) {
        if (dmatidx[j] == i) {
          dmptr = &cpopL[j];
        }
      }
      for (j = 0; j < dmptr->row; j++) {
        dtmp = &crdpop.map[j][nconf];
        dtm2p = dmptr->map[j];
        for (k = 0; k < dmptr->col; k++) {
          dtmp[k] = dtm2p[k];
        }
      }
      k = log10(dmptr->col);
      for (j = 0; j < dmptr->col; j++) {
        sprintf(cfsinp->origins.map[nconf+j], "%s, Frame %*d",
                dmatnames.map[i], k, j+1);
      }
      nconf += dmptr->col;
      DestroyDmat(dmptr);
    }

    // Free allocated memory
    free(crdidx);
    free(dmatidx);
    DestroyCmat(&crdnames);
    DestroyCmat(&dmatnames);
    DestroyCmat(&allnames);
  }
  
  // If requested, fix the lengths of bonds to hydrogen
  if (cfsinp->rattle == 1) {
    CopyDmat(&crdsnapshot, &crdpop, 0);
    VectorRattle(tp, cfsinp, &crdpop, &crdsnapshot, cfsinp->count);
    DestroyDmat(&crdsnapshot);
  }
  
  // Place extra points
  VectorPlaceXpt(tp, &crdpop, cfsinp->count);

  // Alert user as to progress
  if (cfsinp->verbose == 1) {
    printf("Configs >> Created %d initial configurations.\n", cfsinp->count);
  }

  // The matrix now contains a transposed array
  // of all configurations' coordinates.
  return crdpop;
}

// Include four branches of the VecMin.c library to
// get all of the necessary code, highly optimized.
#define DOENERGY 1
#include "VecMin.c"
#undef DOENERGY
#define DOFORCE 1
#include "VecMin.c"
#undef DOFORCE
#define DORESTRAINT 1
#include "VecMin.c"
#undef DORESTRAINT
#define DORESTE 1
#include "VecMin.c"
#undef DORESTE

//-----------------------------------------------------------------------------
// SetBaselineTargets: function to set all target values in a restraint to
//                     the minimum value that the restraint calls for, whether
//                     we are talking about a minimum value relative to the
//                     initial coordinates or an absolute range.
//
// Arguments:
//   cfsinp:  configuration sampling input data
//   myop:    pointer to the operation of interest
//   tp:      topology for the system (helps to determine whether a given file
//            is compatible as a coordinate set for an initial configuration)
//   crdpop:  matrix of coordinates, needed for seeing the coordinates of just
//            the first configuration actually
//-----------------------------------------------------------------------------
static void SetBaselineTargets(configs *cfsinp, cfigop *myop, prmtop *tp,
                               dmat *crdpop)
{
  int i;
  int *anchr;

  // Allocate either the array targets (for any of several one-dimensional
  // restraints) or refcrd (for three-dimensional, positional restraints).
  if (myop->order == 1) {
    myop->refcrd = CreateDmat(3, cfsinp->count, 0);
  }  
  myop->targets = (double*)malloc(cfsinp->count*sizeof(double));

  // Now set the baseline values, whether for an absolute range (i.e.
  // "GridSample" or "RandomSample" keywords) or for a relative perturbation
  // (the "GridPerturb" or "RandomPerturb" keywords).
  if (myop->absrange == 1) {
    if (myop->order == 1) {
      for (i = 0; i < cfsinp->count; i++) {
        myop->refcrd.map[0][i] += myop->minvx;
        myop->refcrd.map[1][i] += myop->minvy;
        myop->refcrd.map[2][i] += myop->minvz;
        myop->targets[i] = myop->distance;
      }
    }
    else {
      SetDVec(myop->targets, cfsinp->count, myop->minval);
    }
  }
  else {
    anchr = myop->anchors;
    if (myop->order == 1) {
      ReflectDVec(myop->refcrd.data, crdpop->map[3*anchr[0]], 3*cfsinp->count);
    }
    else if (myop->order == 2) {
      VectorBond(crdpop, 1, NULL, anchr[0], anchr[1], myop->targets,
                 0.0, 0.0, 0.0, NULL, 0.0, 0.0, cfsinp->count);
    }
    else if (myop->order == 3) {
      VectorAngle(crdpop, 1, NULL, anchr[0], anchr[1], anchr[2],
                  myop->targets, 0.0, NULL, cfsinp->count);
    }
    else {
      VectorTorsion(crdpop, anchr[0], anchr[1], anchr[2], anchr[3],
                    myop->targets, 1, NULL, 0, NULL, NULL, NULL,
                    cfsinp->count);
    }
    if (myop->order == 1) {

      // If this is a positional restraint on just one atom, it can be
      // restrained to sit a particular distance from the target position,
      // and the usual flat bottom half width can then be used to create
      // a permissible zone.
      for (i = 0; i < cfsinp->count; i++) {
       myop->refcrd.map[0][i] += myop->minvx;
       myop->refcrd.map[1][i] += myop->minvy;
       myop->refcrd.map[2][i] += myop->minvz;
       myop->targets[i] = myop->distance;
      }
    }
    else {
      for (i = 0; i < cfsinp->count; i++) {
        myop->targets[i] += myop->minval;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SetIncrementors: function called by InterpretOperations to set the values of
//                  four incrementors.  Either the first or the last three will
//                  be used to iteratively increase the values of restraint
//                  targets in the case of one- or three-dimensional
//                  restraints.
//
// Arguments:
//   myop:      the operation in question
//   finc:      the incrementor for one-dimensional restraints (returned)
//   fxinc:
//   fyinc:     incrementors for three-dimensional restraints (returned)
//   fzinc:
//   dimcount:  number of steps along a line, or total number of points to
//              sample out of a three-dimensional volume
//-----------------------------------------------------------------------------
static void SetIncrementors(cfigop *myop, double *finc, double *fxinc,
                           double *fyinc, double *fzinc, int dimcount)
{
  int dimc3;

  if (myop->order == 1) {
    if (myop->pegtype == 0) {
      dimc3 = pow(dimcount, 0.333);
      while(dimc3*dimc3*dimc3 < dimcount) {
        dimc3++;
      }
    }
    else {
      dimc3 = dimcount;
    }
    if (dimc3 == 0) {
      dimc3 = 1;
    }
    *fxinc = (myop->maxvx - myop->minvx)/dimc3;
    *fyinc = (myop->maxvy - myop->minvy)/dimc3;
    *fzinc = (myop->maxvz - myop->minvz)/dimc3;
    *finc = 0.0;
  }
  else {
    *finc = (myop->maxval - myop->minval)/dimcount;
  }
}

//-----------------------------------------------------------------------------
// SetSliderValue: set the value of a target along one or three dimensions.
//                 This encapsulates code that would otherwise be repeated in
//                 InterpretOperations below.
//
// Arguments:
//   myop:      the operation in question
//   finc:      the incrementor for one-dimensional restraints (returned)
//   fxinc:
//   fyinc:     incrementors for three-dimensional restraints (returned)
//   fzinc:
//   dimc3:     cube root of the number of steps along each dimension, invoked
//              if the order of the operation is 1 and its pegtype is 0 (spread
//              points over a whole grid volume)
//   step:      the number of the step along one axis of the grid
//   cidx:      index of the configuration for this particular target
//-----------------------------------------------------------------------------
static void SetSliderValue(cfigop *myop, double finc, double fxinc,
                          double fyinc, double fzinc, int dimc3, int step,
                           int cidx)
{
  int ix, iy, iz;

  if (myop->order == 1) {
    if (myop->pegtype == 0) {
      ix = MAX(0, step / (dimc3*dimc3));
      iy = MAX(0, (step - ix*dimc3*dimc3)/dimc3);
      iz = MAX(0, (step - (ix*dimc3 + iy)*dimc3));
      myop->refcrd.map[0][cidx] += (ix + 0.5)*fxinc;
      myop->refcrd.map[1][cidx] += (iy + 0.5)*fyinc;
      myop->refcrd.map[2][cidx] += (iz + 0.5)*fzinc;
    }
    else {
      ix = iy = iz = 0;   
      myop->refcrd.map[0][cidx] += (ix + 0.5)*fxinc;
      myop->refcrd.map[1][cidx] += (iy + 0.5)*fyinc;
      myop->refcrd.map[2][cidx] += (iz + 0.5)*fzinc;
    }
  }
  else {
    myop->targets[cidx] += (step + 0.5)*finc;
  }
}

//-----------------------------------------------------------------------------
// InterpretCombination: interpret a sinble combination.  This just tidies up
//                       the much more elaborate InterpretOperations function
//                       below by assigning integer indices to the labels that
//                       the combination entails.  Returns 1 if successful, 0
//                       if not.
//
// Arguments:
//   cfsinp:  configuration sampling input data
//   nc:      the number of the combination to interpret
//-----------------------------------------------------------------------------
static int InterpretCombination(configs *cfsinp, int nc)
{
  int i, j;

  for (i = 0; i < cfsinp->combo[nc].npcs; i++) {
    cfsinp->combo[nc].pcs[i] = -1;
    for (j = 0; j < cfsinp->nops; j++) {
      if (strcmp(cfsinp->combo[nc].labels[i], cfsinp->ops[j].label) == 0) {
        cfsinp->combo[nc].pcs[i] = j;
      }
    }
    if (cfsinp->combo[nc].pcs[i] < 0) {
      if (WordIsInteger(cfsinp->combo[nc].labels[i]) == 1) {
        cfsinp->combo[nc].pcs[i] = atoi(cfsinp->combo[nc].labels[i])-1;
      }
      if (cfsinp->combo[nc].pcs[i] >= cfsinp->nops ||
          cfsinp->combo[nc].pcs[i] < 0) {
        return 0;
      }
    }
  }

  return 1;
}

//-----------------------------------------------------------------------------
// IsPolarHydrogen: test if an atom qualifies as a polar hydrogen.  Its own
//                  atomic number must be 1, it must bond to only one other
//                  atom, and the atomic number of that atom must correspond
//                  to some electronegative atom.
//
// Arguments:
//   tp:        the topology
//   ai:        the atom to look at
//-----------------------------------------------------------------------------
static int IsPolarHydrogen(prmtop *tp, int ai)
{
  int i, j, nwithmass, nelec;

  // Make sure it's a hydrogen
  if (tp->ZNumber[ai] != 1) {
    return 0;
  }
  
  // If there are too many bonds to things with mass, this is not a polar H.
  // Make sure any massive atom it bonds to is electronegative.
  nwithmass = 0;
  nelec = 0;
  for (i = 0; i < tp->nb1234[ai].n12; i++) {
    j = tp->nb1234[ai].L12[i];
    if (tp->Masses[j] > 1.0e-8) {
      nwithmass++;
      if (tp->ZNumber[j] == 7 || tp->ZNumber[j] == 8 || tp->ZNumber[j] == 15 ||
          tp->ZNumber[j] == 16) {
        nelec++;
      }
    }
  }
  if (nwithmass != 1 || nelec == 0) {
    return 0;
  }
  
  // Looks like a polar hydrogen
  return 1;
}

//-----------------------------------------------------------------------------
// StopHydrogenBonding: create a series of restraints to limit hydrogen bonding
//                      among potential donors and acceptors.  Pairs which
//                      already have restraints defined will be ignored, and
//                      pairs which are already 1-2, 1-3, or 1-4 to one another
//                      will also be ignored.
//
// Arguments:
//   cfsinp:  configuration sampling input data
//   tp:      topology for the system (contains atom properties and bonding
//            patterns)
//-----------------------------------------------------------------------------
static void StopHydrogenBonding(configs *cfsinp, prmtop *tp)
{
  int i, j, k, iatm, jatm, nacc, ndon, nheavy, nbrs, nrst, resid, maxops;
  int bigatom;
  int* donors;
  int* acceptors;
  
  // Find donors and acceptors
  nacc = 0;
  ndon = 0;
  donors = (int*)malloc(tp->natom * sizeof(int));
  acceptors = (int*)malloc(tp->natom * sizeof(int));
  for (i = 0; i < tp->natom; i++) {
    if (IsPolarHydrogen(tp, i) == 1) {
      donors[ndon] = i;
      ndon++;
    }
    else if (tp->ZNumber[i] == 8 || tp->ZNumber[i] == 16) {
      acceptors[nacc] = i;
      nacc++;
    }
    else if ((tp->ZNumber[i] == 7 || tp->ZNumber[i] == 15) &&
             tp->nb1234[i].n12 < 4) {
      nheavy = 0;
      for (j = 0; j < tp->nb1234[i].n12; j++) {
        jatm = tp->nb1234[i].L12[j];
        nheavy += (tp->ZNumber[jatm] > 1);
      }
      if ((tp->ZNumber[i] == 7  && nheavy < 3) ||
          (tp->ZNumber[i] == 15 && nheavy < 4)){
        acceptors[nacc] = i;
        nacc++;
      }
    }
  }
  
  // Create restraints
  maxops = cfsinp->nops + 32;
  nrst = 1;
  cfsinp->ops = (cfigop*)realloc(cfsinp->ops, maxops * sizeof(cfigop));
  for (i = 0; i < nacc; i++) {
    iatm = acceptors[i];
    for (j = 0; j < ndon; j++) {
      jatm = donors[j];
      nbrs = 0;
      for (k = 0; k < tp->nb1234[iatm].n12; k++) {
        if (tp->nb1234[iatm].L12[k] == jatm) {
          nbrs = 1;
        }
      }
      for (k = 0; k < tp->nb1234[iatm].n13; k++) {
        if (tp->nb1234[iatm].L13[k] == jatm) {
          nbrs = 1;
        }        
      }
      for (k = 0; k < tp->nb1234[iatm].n14; k++) {
        if (tp->nb1234[iatm].L14[k] == jatm) {
          nbrs = 1;
        }        
      }
      for (k = 0; k < cfsinp->nops; k++) {
        if (cfsinp->ops[i].order == 2 &&
            ((cfsinp->ops[i].anchors[0] == iatm &&
              cfsinp->ops[i].anchors[1] == jatm) ||
             (cfsinp->ops[i].anchors[0] == jatm &&
              cfsinp->ops[i].anchors[1] == iatm))) {
          nbrs = 1;
        }
      }
      if (nbrs == 1) {
        continue;
      }

      // Create a new restraint operation using the flat-bottom,
      // half-width trick to penalize hydrogen bond formation.
      cfsinp->ops[cfsinp->nops].order      = 2;
      cfsinp->ops[cfsinp->nops].absrange   = 1;
      cfsinp->ops[cfsinp->nops].scatter    = 0;
      cfsinp->ops[cfsinp->nops].anchors[0] = iatm;
      cfsinp->ops[cfsinp->nops].anchors[1] = jatm;
      cfsinp->ops[cfsinp->nops].pegtype    = 0;
      cfsinp->ops[cfsinp->nops].Krst       = 64.0;
      cfsinp->ops[cfsinp->nops].quadwin    = 1.0;
      cfsinp->ops[cfsinp->nops].distance   = 0.0;
      bigatom = (tp->ZNumber[iatm] >= 15);
      cfsinp->ops[cfsinp->nops].minval     = 1000.0 + (0.25 * bigatom);
      cfsinp->ops[cfsinp->nops].maxval     = cfsinp->ops[cfsinp->nops].minval;
      cfsinp->ops[cfsinp->nops].fbhw       = cfsinp->ops[cfsinp->nops].minval -
                                             cfsinp->HBMargin;
      cfsinp->ops[cfsinp->nops].atommasks = CreateCmat(2, 32);
      sprintf(cfsinp->ops[cfsinp->nops].label, "AutoHB_%d\n", nrst);
      resid = LocateResID(tp, iatm, 0, tp->nres);
      sprintf(cfsinp->ops[cfsinp->nops].atommasks.map[0], ":%d@%-4.4s",
              resid, &tp->AtomNames[iatm]);
      resid = LocateResID(tp, jatm, 0, tp->nres);
      sprintf(cfsinp->ops[cfsinp->nops].atommasks.map[1], ":%d@%-4.4s",
              resid, &tp->AtomNames[jatm]);
      cfsinp->ops[cfsinp->nops].targets = (double*)malloc(cfsinp->count *
                                                          sizeof(double));
      SetDVec(cfsinp->ops[cfsinp->nops].targets, cfsinp->count,
              cfsinp->ops[cfsinp->nops].minval);
      nrst++;
      cfsinp->nops += 1;
      if (cfsinp->nops == maxops) {
        maxops += 32;
        cfsinp->ops = (cfigop*)realloc(cfsinp->ops, maxops * sizeof(cfigop));
      }
    }
  }
  cfsinp->ops = (cfigop*)realloc(cfsinp->ops, cfsinp->nops * sizeof(cfigop));
  
  // Free allocated memory
  free(donors);
  free(acceptors);
}

//-----------------------------------------------------------------------------
// SpinHydroxylGroups: induce random rotations of hydroxyl (and thiol) groups
//                     about the nearest dihedral axis.  Hydroxyl groups which
//                     are already affected by a dihedral restraint will be
//                     ignored, but distance restraints (such as those that
//                     prevent hydrogen bonding) will not interfere with the
//                     operations.
//
// Arguments:
//   cfsinp:  configuration sampling input data
//   tp:      topology for the system
//   tj:      trajectory control information, but needed for the random number
//            counter
//-----------------------------------------------------------------------------
static void SpinHydroxylGroups(configs *cfsinp, prmtop *tp, trajcon *tj)
{
  int i, j, k, m, jatm, katm, nhydro, resid, found, success;
  int* hydrogens;
  int* oxygens;
  int* heavy1;
  int* heavy2;

  nhydro = 0;
  hydrogens = (int*)malloc(tp->natom * sizeof(int));
  oxygens   = (int*)malloc(tp->natom * sizeof(int));
  heavy1    = (int*)malloc(tp->natom * sizeof(int));
  heavy2    = (int*)malloc(tp->natom * sizeof(int));
  for (i = 0; i < tp->natom; i++) {

    // Simplest check at the front
    if (tp->ZNumber[i] != 1 || tp->nb1234[i].n14 == 0) {
      continue;
    }
    
    // Skip if there are already torsion restraints acting on this
    found = 0;
    for (j = 0; j < cfsinp->nops; j++) {
      if (cfsinp->ops[j].order == 4) {
	if (cfsinp->ops[j].anchors[0] == i || cfsinp->ops[j].anchors[1] == i ||
	    cfsinp->ops[j].anchors[2] == i || cfsinp->ops[j].anchors[3] == i) {
	  found = 1;	  
	}
      }
    }
    if (found == 1) {
      continue;
    }

    // Look for a rotatable bond
    if (IsPolarHydrogen(tp, i) == 1 && tp->nb1234[i].n12 > 0) {
      if (tp->ZNumber[tp->nb1234[i].L12[0]] == 8 ||
          tp->ZNumber[tp->nb1234[i].L12[0]] == 16) {
	success = 0;
        for (j = 0; j < tp->nb1234[i].n13; j++) {
	  jatm = tp->nb1234[i].L13[j];
	  for (k = 0; k < tp->nb1234[jatm].n12; k++) {
	    katm = tp->nb1234[jatm].L12[k];
	    for (m = 0; m < tp->nb1234[i].n14; m++) {
	      if (tp->nb1234[i].L14[m] == katm) {
                hydrogens[nhydro] = i;
                oxygens[nhydro] = tp->nb1234[i].L12[0];
		if (success == 1) {
		  if (tp->ZNumber[heavy1[nhydro]] *
		      tp->ZNumber[heavy2[nhydro]] < tp->ZNumber[jatm] *
                                                    tp->ZNumber[katm]) {
                    heavy1[nhydro]  = jatm;
                    heavy2[nhydro]  = katm;
		  }
		}
		else {
                  heavy1[nhydro]  = jatm;
                  heavy2[nhydro]  = katm;
		}
		success = 1;
	      }
	    }
	  }
	}
	if (success == 1) {
	  nhydro++;
	}
      }
    }
  }

  // Set hydroxyl group restraints to spin them randomly into new rotameric
  // positions.  Do not try to force them into a non-equilibrium state.
  cfsinp->ops = (cfigop*)realloc(cfsinp->ops,
                                 (cfsinp->nops + nhydro) * sizeof(cfigop));
  for (i = 0; i < nhydro; i++) {

    // Create a new restraint operation using the flat-bottom,
    // half-width trick to penalize hydrogen bond formation.
    cfsinp->ops[cfsinp->nops].order      = 4;
    cfsinp->ops[cfsinp->nops].absrange   = 1;
    cfsinp->ops[cfsinp->nops].scatter    = 0;
    cfsinp->ops[cfsinp->nops].anchors[0] = hydrogens[i];
    cfsinp->ops[cfsinp->nops].anchors[1] = oxygens[i];
    cfsinp->ops[cfsinp->nops].anchors[2] = heavy1[i];
    cfsinp->ops[cfsinp->nops].anchors[3] = heavy2[i];
    cfsinp->ops[cfsinp->nops].pegtype    = 0;
    cfsinp->ops[cfsinp->nops].Krst       = cfsinp->KSpinOH;
    cfsinp->ops[cfsinp->nops].quadwin    = 1.0;
    cfsinp->ops[cfsinp->nops].fbhw       = cfsinp->fbhwSpinOH;
    cfsinp->ops[cfsinp->nops].distance   = 0.0;
    cfsinp->ops[cfsinp->nops].minval     = 2.0*ran2(&tj->rndcon) - 1.0;
    cfsinp->ops[cfsinp->nops].maxval     = cfsinp->ops[cfsinp->nops].minval;
    cfsinp->ops[cfsinp->nops].atommasks = CreateCmat(4, 32);
    sprintf(cfsinp->ops[cfsinp->nops].label, "SpinOH_%d\n", i);
    resid = LocateResID(tp, hydrogens[i], 0, tp->nres);
    sprintf(cfsinp->ops[cfsinp->nops].atommasks.map[0], ":%d@%-4.4s",
            resid, &tp->AtomNames[hydrogens[i]]);
    resid = LocateResID(tp, oxygens[i], 0, tp->nres);
    sprintf(cfsinp->ops[cfsinp->nops].atommasks.map[1], ":%d@%-4.4s",
            resid, &tp->AtomNames[oxygens[i]]);
    resid = LocateResID(tp, heavy1[i], 0, tp->nres);
    sprintf(cfsinp->ops[cfsinp->nops].atommasks.map[2], ":%d@%-4.4s",
            resid, &tp->AtomNames[heavy1[i]]);
    resid = LocateResID(tp, heavy2[i], 0, tp->nres);
    sprintf(cfsinp->ops[cfsinp->nops].atommasks.map[3], ":%d@%-4.4s",
            resid, &tp->AtomNames[heavy2[i]]);
    cfsinp->ops[cfsinp->nops].targets = (double*)malloc(cfsinp->count *
                                                        sizeof(double));
    for (j = 0; j < cfsinp->count; j++) {
      cfsinp->ops[cfsinp->nops].targets[j] =
        (2.0*ran2(&tj->rndcon) - 1.0) * PI;
    }
    cfsinp->nops += 1;
  }
  
  // Free allocated memory
  free(hydrogens);
  free(oxygens);
  free(heavy1);
  free(heavy2);
}

//-----------------------------------------------------------------------------
// InterpretOperations: with the topology and coordinates of a system, we can
//                      understand the meaning of ambmask strings pointing to
//                      atoms in the system.
//
// Arguments:
//   cfsinp:  configuration sampling input data
//   tp:      topology for the system (helps to determine whether a given file
//            is compatible as a coordinate set for an initial configuration)
//   crdpop:  matrix of coordinates, needed for seeing the coordinates of just
//            the first configuration actually
//   tj:      trajectory control data (needed for the random number counter)
//-----------------------------------------------------------------------------
static void InterpretOperations(configs *cfsinp, prmtop *tp, dmat *crdpop,
                                trajcon *tj)
{
  int i, j, k, m, jp, kp, mp, opj, opk, opm, msum;
  int dimcount, dimc3, cidx;
  int *anchr;
  int* amask;
  int* covered;
  double tbase, trange, trangex, trangey, trangez, dj, jinc, kinc, minc;
  double jxinc, jyinc, jzinc, kxinc, kyinc, kzinc, mxinc, myinc, mzinc;
  coord tc;

  // With the topology now understood, we can assign atoms to the restraints
  tc = CreateCoord(tp->natom);
  for (i = 0; i < 3*tp->natom; i++) {
    tc.loc[i] = crdpop->map[i][0];
  }
  for (i = 0; i < cfsinp->nops; i++) {
    for (j = 0; j < cfsinp->ops[i].order; j++) {
      amask = ParseAmbMask(cfsinp->ops[i].atommasks.map[j], tp, &tc);
      msum = ISum(amask, tp->natom);
      if (msum != 1) {
        printf("InterpretOperations >> Ambmask string %s matches %d atoms.\n",
               cfsinp->ops[i].atommasks.map[j], msum); 
        exit(1);
      }
      else {
        for (k = 0; k < tp->natom; k++) {
          if (amask[k] == 1) {
            cfsinp->ops[i].anchors[j] = k;
          }
        }
      }
      free(amask);
    }
  }

  // Combinatorial sampling is also an option: the restraints can be combined
  // to create orthogonal sampling in up to three dimensions at once.  If no
  // combinations are specified, then grid-based sampling proceeds along a
  // single line, with all restraint targets progressing from their minimum to
  // maximum ranges. (Random sampling will still cover two dimensions, but with
  // no guarantee of regularity!) If two or three operations are flagged as a
  // combo, then the grid sampling will cover a two- or three-dimensional
  // space with the same number of points.
  covered = (int*)calloc(cfsinp->nops, sizeof(int));
  for (i = 0; i < cfsinp->ncombo; i++) {

    // Interpret the combination.  Match every label it was given to the
    // label of a restraint operation.
    if (InterpretCombination(cfsinp, i) == 0) {
      continue;
    }

    // The number of samples along each dimension will be the nth
    // root of the total number of configurations, rounding up--it's
    // on the user to provide a good total number.  If this causes
    // there to be fewer configurations than complete gridded sampling
    // of the space would require, the space will not be fully sampled. 
    dimcount = pow((double)cfsinp->count, 1.0/cfsinp->combo[i].npcs);
    dimcount += 1.0 - 1.0e-8;

    // Mark all operations in this combination as covered
    for (j = 0; j < cfsinp->combo[i].npcs; j++) {
      covered[cfsinp->combo[i].pcs[j]] = 1;
    }

    // Any combo will have at least two, but maybe three, members
    opj = cfsinp->combo[i].pcs[0];
    SetIncrementors(&cfsinp->ops[opj], &jinc, &jxinc, &jyinc, &jzinc,
                   dimcount);
    SetBaselineTargets(cfsinp, &cfsinp->ops[opj], tp, crdpop);
    opk = cfsinp->combo[i].pcs[1];
    SetIncrementors(&cfsinp->ops[opk], &kinc, &kxinc, &kyinc, &kzinc,
                   dimcount);
    SetBaselineTargets(cfsinp, &cfsinp->ops[opk], tp, crdpop);
    dimc3 = pow(dimcount, 0.333);
    while(dimc3*dimc3*dimc3 < dimcount) {
      dimc3++;
    }
    cidx = 0;

    // Two-operation combo
    if (cfsinp->combo[i].npcs == 2) {
      for (j = 0; j < dimcount; j++) {
        for (k = 0; k < dimcount; k++) {

          // Bail out of the nested loop if we've
          // hit the allotted configuration count
          if (cidx >= cfsinp->count) {
            j = dimcount;
            k = dimcount;
            continue;
          }

          // Combinations of operations are exclusively doing grid sampling
          SetSliderValue(&cfsinp->ops[opj], jinc, jxinc, jyinc, jzinc, dimc3,
                         j, cidx);
          SetSliderValue(&cfsinp->ops[opk], kinc, kxinc, kyinc, kzinc, dimc3,
                         k, cidx);
          cidx++;
        }
      }
    }

    // Three-operation combo
    if (cfsinp->combo[i].npcs == 3) {
      opm = cfsinp->combo[i].pcs[2];
      SetIncrementors(&cfsinp->ops[opm], &minc, &mxinc, &myinc, &mzinc,
                      dimcount);
      SetBaselineTargets(cfsinp, &cfsinp->ops[opm], tp, crdpop);
      for (j = 0; j < dimcount; j++) {
        for (k = 0; k < dimcount; k++) {
          for (m = 0; m < dimcount; m++) {

            // Bail out of the nested loop if we've
            // hit the allotted configuration count
            if (cidx >= cfsinp->count) {
              j = dimcount;
              k = dimcount;
              m = dimcount;
              continue;
            }

            // Combinations of operations are exclusively doing grid sampling
            SetSliderValue(&cfsinp->ops[opj], jinc, jxinc, jyinc, jzinc, dimc3,
                           j, cidx);
            SetSliderValue(&cfsinp->ops[opk], kinc, kxinc, kyinc, kzinc, dimc3,
                           k, cidx);
            SetSliderValue(&cfsinp->ops[opm], minc, mxinc, myinc, mzinc, dimc3,
                           m, cidx);
            cidx++;
          }
        }
      }
    }
  }

  // Assign restraint targets for all remaining (singlet) operations.
  // Random sampling is permitted in singlet operations, but not in combos.
  for (i = 0; i < cfsinp->nops; i++) {
    if (covered[i] == 1) {
      continue;
    }
    if (cfsinp->ops[i].order == 1) {
      dimc3 = pow(cfsinp->count, 0.333);
      while (dimc3*dimc3*dimc3 < cfsinp->count) {
       dimc3++;
      }
    }
    else {
      dimc3 = 0;
    }
    SetBaselineTargets(cfsinp, &cfsinp->ops[i], tp, crdpop);

    // Random sampling
    if (cfsinp->ops[i].scatter == 1) {
      if (cfsinp->ops[i].order == 1) {
        trangex = cfsinp->ops[i].maxvx - cfsinp->ops[i].minvx;
        trangey = cfsinp->ops[i].maxvy - cfsinp->ops[i].minvy;
        trangez = cfsinp->ops[i].maxvz - cfsinp->ops[i].minvz;
      }
      else {
        trange = cfsinp->ops[i].maxval - cfsinp->ops[i].minval;
      }
      for (j = 0; j < cfsinp->count; j++) {
        if (cfsinp->ops[i].order == 1) {
          cfsinp->ops[i].refcrd.map[0][j] += trangex*ran2(&tj->rndcon);
          cfsinp->ops[i].refcrd.map[1][j] += trangey*ran2(&tj->rndcon);
          cfsinp->ops[i].refcrd.map[2][j] += trangez*ran2(&tj->rndcon);
        }
        else {
          cfsinp->ops[i].targets[j] += trange*ran2(&tj->rndcon);
        }
      }
    }

    // Grid sampling
    else {
      SetIncrementors(&cfsinp->ops[i], &jinc, &jxinc, &jyinc, &jzinc,
                     cfsinp->count);
      for (j = 0; j < cfsinp->count; j++) {
        SetSliderValue(&cfsinp->ops[i], jinc, jxinc, jyinc, jzinc,
                       dimc3, j, j);
      }
    }
  }

  // Assign general hydrogen-bond prevention restraints
  if (cfsinp->PreventHBonds == 1) {
    StopHydrogenBonding(cfsinp, tp);
  }

  // Assign general hydroxyl-spinning restraints
  if (cfsinp->SpinHydroxyls == 1) {
    SpinHydroxylGroups(cfsinp, tp, tj);
  }
  
  
  // Alert user as to progress
  if (cfsinp->verbose == 1) {
    printf("Configs >> Processed %d restraint operations.\n", cfsinp->nops);
  }

  // Free allocated memory
  DestroyCoord(&tc);
  free(covered);
}

//-----------------------------------------------------------------------------
// VectorBondContrib: execute the vectorized bonded interaction calculations
//                    for an array of configurations of a single system.  This
//                    will not process restraints--those are done in a separate
//                    function, VectorRestraints.
//
// Arguments:
//   tp:       system topology
//   mobile:   array of flags to state whether each atom is movable or not
//   crdpop:   matrix of coordinates, giving X1, Y1, Z1, X2, ..., Zn on
//             successive rows
//   frcpop:   matrix of forces, being accumulated before this function is
//             called and accepting contributions from this function
//   sysUV:    vector of energy accumulators for each conformation
//   compfrc:  flag to have forces computed if set to 1 (otherwise only the
//             energy will be computed)
//   nactive:  the number of configurations actively being manipulated
//-----------------------------------------------------------------------------
static void VectorBondContrib(prmtop *tp, int* mobile, dmat *crdpop,
                              dmat *frcpop, Energy* sysUV, int compfrc,
                              int nactive)
{
  int i, j, k, m;
  double phival;
  double Nperiod[6], Kamp[6];
  double* EqVal;
  double* theta;
  bondcomm *thisbc;
  anglcomm *thisac;
  dihecomm *thishc;
  cmapcomm *thismc;
  dmat Phase;

  // Allocate arrays to hold equilibrium lengths, angles, or phase angles
  // from the topology or restraints, as well as calculated angles from the
  // configurations.
  EqVal = (double*)malloc(nactive*sizeof(double));
  theta = (double*)malloc(nactive*sizeof(double));
  Phase = CreateDmat(nactive, 6, 0);

  // Loop over all bonds in the system
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->BLC[i].nbond; j++) {
      thisbc = &tp->BLC[i].BC[j];
      if (mobile[thisbc->a] + mobile[thisbc->b] == 0) {
        continue;
      }
      SetDVec(EqVal, nactive, tp->BParam[thisbc->t].l0);
      if (compfrc == 1) {
        VectorBondF(crdpop, frcpop, sysUV, thisbc->a, thisbc->b, theta,
                    tp->BParam[thisbc->t].K, tp->BParam[thisbc->t].Kpull,
                    tp->BParam[thisbc->t].Kpress, EqVal,
                    tp->BParam[thisbc->t].lpull0,
                    tp->BParam[thisbc->t].lpress0, nactive);
      }
      else {
        VectorBond(crdpop, 0, sysUV, thisbc->a, thisbc->b, theta,
                   tp->BParam[thisbc->t].K, tp->BParam[thisbc->t].Kpull,
                   tp->BParam[thisbc->t].Kpress, EqVal,
                   tp->BParam[thisbc->t].lpull0,
                   tp->BParam[thisbc->t].lpress0, nactive);
      }
    }
  }

  // Loop over all angles in the system
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->ALC[i].nangl; j++) {
      thisac = &tp->ALC[i].AC[j];
      if (mobile[thisac->a] + mobile[i] + mobile[thisac->c] == 0) {
        continue;
      }
      SetDVec(EqVal, nactive, tp->AParam[thisac->t].th0);
      if (compfrc == 1) {
        VectorAngleF(crdpop, frcpop, sysUV, thisac->a, i, thisac->c, theta,
                     tp->AParam[thisac->t].K, EqVal, nactive);
      }
      else {
        VectorAngle(crdpop, 0, sysUV, thisac->a, i, thisac->c, theta,
                    tp->AParam[thisac->t].K, EqVal, nactive);
      }
    }
  }

  // Loop over all dihedrals in the system
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->HLC[i].ndihe; j++) {
      thishc = &tp->HLC[i].HC[j];
      if (mobile[thishc->a] + mobile[i] +
          mobile[thishc->c] + mobile[thishc->d] == 0) {
        continue;
      }

      // Repack the dihedral into arrays that can be fed into the
      // vectorized calculation
      for (k = 0; k < thishc->nt; k++) {
        Kamp[k] = tp->HParam[thishc->t[k]].K;
        Nperiod[k] = tp->HParam[thishc->t[k]].N;
        phival = tp->HParam[thishc->t[k]].Phi;
        for (m = 0; m < nactive; m++) {
          Phase.map[m][k] = phival;
        }
      }
      if (compfrc == 1) {
        VectorTorsionF(crdpop, frcpop, thishc->a, i, thishc->c, thishc->d,
                       theta, sysUV, thishc->nt, Kamp, Nperiod, &Phase,
                       nactive);
      }
      else {
        VectorTorsion(crdpop, thishc->a, i, thishc->c, thishc->d, theta, 0,
                      sysUV, thishc->nt, Kamp, Nperiod, &Phase, nactive);
      }
    }
  }

  // Loop over all CMAPs in the system
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->MLC[i].ncmap; j++) {
      thismc = &tp->MLC[i].MC[j];
      if (mobile[thismc->a] + mobile[thismc->b] + mobile[i] +
          mobile[thismc->d] + mobile[thismc->e] == 0) {
        continue;
      }

      // There are no CMAP-like restraints, but a vectorized calculation
      // must still happen.
      if (compfrc == 1) {
        VectorCMAPF(crdpop, frcpop, thismc->a, thismc->b, i, thismc->d,
                    thismc->e, sysUV, thismc, tp, nactive);
      }
      else {
        VectorCMAP(crdpop, thismc->a, thismc->b, i, thismc->d,
                   thismc->e, sysUV, thismc, tp, nactive);
      }
    }
  }
  
  // Free allocated memory
  free(EqVal);
  free(theta);
  DestroyDmat(&Phase);
}

//-----------------------------------------------------------------------------
// VectorTrackRst: function for implementing track restraints.  A track
//                 restraint applies to one or two atoms, as specified, each of
//                 these real atoms in the system being tied to fictitious
//                 particles that travel strictly along a pre-defined track
//                 (the "track" particles).  Like the reference coordinates for
//                 point positional restraints, the track is pegged to the
//                 board.  When applying track restraints, the positions of the
//                 track particles are instantaneously equilibrated so as to
//                 minimize the energy of the (piecewise harmonic) restraints
//                 holding the track particles to their real atom partners and
//                 (if there are two, not just one fictitious particles), the
//                 restraint energy holding the track particles close to one
//                 another.  In some sense, a track restraint is a generalized
//                 point positional restraint.
//
// Arguments:
//-----------------------------------------------------------------------------
#if 0
static void VectorTrackRst(track paths, dmat *crdpop, dmat *frcpop,
                           dmat *sttspace, Energy* sysUV, Energy* trialUV,
                           int compfrc, int nactive)
{
  int i, j, ip1, ip2;
  double dp1, dp2, vp1, vp2, vp1fac, vp2fac, E0;
  double rp1[3], rp2[3], frp1[3], frp2[3];
  double *xtmp1, *xtmp2, *ytmp1, *ytmp2, *ztmp1, *ztmp2;
  dmat frctmp;

  // Set pointers in the dummy matrix struct
  frctmp.map = (double**)malloc(12*sizeof(double));
  frctmp.data = sttspace->map[12];
  for (i = 0; i < 12; i++) {
    frctmp.map[i] = sttspace->map[i+12];
  }
  
  // Get the atoms of the system that attach to the track particles.
  for (i = 0; i < nactive; i++) {
    vp1 = paths.p1[i];
    ip1 = vp1;
    dp1 = vp1 - ip1;
    if (paths[i].ntp == 2) {
      vp2 = paths.p2[i];
      ip2 = vp2;
      dp2 = vp2 - ip2;
    }
    for (j = 0; j < 3; j++) {
      sttspace->map[j][i] = crdpop->map[3*paths.atmA+j][i];
      frp1[j] = 0.0;
      if (paths[i].ntp == 2) {
        sttspace->map[j+9][i] = crdpop->map[3*paths.atmB+j][i];
        frp2[j] = 0.0;
      }
    }
  }

  // The initial locations of the track particles are given by their
  // progress stored in their respective paths.
  SetDVec(sttspace->map[3], 6*sttspace->col, 0.0);
  for (pwr = 0; pwr < 3; pwr++) {
    xtmp1 = paths.xcoef.map[3*ip1+pwr];
    ytmp1 = paths.ycoef.map[3*ip1+pwr];
    ztmp1 = paths.zcoef.map[3*ip1+pwr];
    vp1fac = 1.0;
    if (paths[i].ntp == 2) {
      xtmp2 = paths.xcoef.map[3*ip2+pwr];
      ytmp2 = paths.ycoef.map[3*ip2+pwr];
      ztmp2 = paths.zcoef.map[3*ip2+pwr];
      vp2fac = 1.0;
    }
    for (i = 0; i < 2-pwr; i++) {
      vp1fac *= vp1;
      if (paths[i].ntp == 2) {
        vp2fac *= vp2;
      }
    }
    for (i = 0; i < nactive; i++) {
      sttspace->map[3][i] += vp1fac*xtmp1[i];
      sttspace->map[4][i] += vp1fac*ytmp1[i];
      sttspace->map[5][i] += vp1fac*ztmp1[i];
    }
    if (paths[i].ntp == 2) {
      for (i = 0; i < nactive; i++) {
        sttspace->map[6][i] += vp2fac*xtmp2[i];
        sttspace->map[7][i] += vp2fac*ytmp2[i];
        sttspace->map[8][i] += vp2fac*ztmp2[i];
      }
    }
  }

  // Compute the track-to-real energy as a sum of one or two restraints.
  // Contribute the forces if necessary.
  VectorZeroEnergies(trialUV, nactive);
  VectorBondR(sttspace, &frctmp, trialUV, 0, 1, sttspace->map[24],
              cfsinp->ops[i].Krst, cfsinp->ops[i].targets,
              cfsinp->ops[i].fbhw, cfsinp->ops[i].quadwin, nactive);
  if (paths[i].ntp == 2) {
    VectorBondR(sttspace, &frctmp, trialUV, 2, 3, sttspace->map[24],
                cfsinp->ops[i].Krst, cfsinp->ops[i].targets,
                cfsinp->ops[i].fbhw, cfsinp->ops[i].quadwin, nactive);
  }

  // Compute the energy along the track 
}
#endif

//-----------------------------------------------------------------------------
// VectorApplyRestraints: this function applies restraints analogous to each
//                        potential form in VectorBondContrib above.  As such,
//                        it takes many of the same arguments, but it is kept
//                        separate from the other routine so that bonded
//                        energies can be computed separately.
//
// Arguments:
//   cfsinp:   configuration sampling input data
//   crdpop:   matrix of coordinates, giving X1, Y1, Z1, X2, ..., Zn on
//             successive rows
//   frcpop:   matrix of forces, being accumulated before this function is
//             called and accepting contributions from this function
//   sttspace: scratch space for coordinate and force storage, used in position
//             restraints
//   sysUV:    vector of energy accumulators for each conformation
//   compfrc:  flag to have forces computed if set to 1 (otherwise only the
//             energy will be computed)
//   nactive:  the number of configurations actively being manipulated
//-----------------------------------------------------------------------------
static void VectorApplyRestraints(configs *cfsinp, dmat *crdpop, dmat *frcpop,
                                  dmat *sttspace, Energy* sysUV, int compfrc,
                                  int nactive)
{
  int i, j, tmpex;
  int *anch;
  int *mobile;
  double* Kamp;
  double* theta;
  dmat Phase, frctmp, crdtmp;

  // Allocate memory for scratch computations
  theta = (double*)malloc(cfsinp->count*sizeof(double));
  Phase = CreateDmat(cfsinp->count, 1, 0);
  tmpex = 0;
  
  // Pointer to the movable atoms array
  mobile = cfsinp->movable;

  // Loop over all restraint operations and execute them.
  // (The equilibria of each restraint will have been
  // laid beforehand--the implementation is now done with
  // the constants for flat bottom half width, stiffness,
  // and maximum force / quadratic range.)
  for (i = 0; i < cfsinp->nops; i++) {
    anch = cfsinp->ops[i].anchors;
    if (cfsinp->ops[i].order == 1) {
      if (mobile[anch[0]] == 0) {
       continue;
      }

      // Set up some dummy variables
      if (tmpex == 0) {
       frctmp.map = (double**)malloc(6*sizeof(double*));
       frctmp.row = 6;
       frctmp.col = cfsinp->count;
       crdtmp.map = (double**)malloc(6*sizeof(double*));
       crdtmp.row = 6;
       crdtmp.col = cfsinp->count;
       tmpex = 1;
      }
      crdtmp.data = crdpop->map[3*anch[0]];
      frctmp.data = sttspace->data;
      for (j = 0; j < 3; j++) {
        crdtmp.map[j] = crdpop->map[3*anch[0]+j];
        crdtmp.map[j+3] = cfsinp->ops[i].refcrd.map[j];
        frctmp.map[j] = frcpop->map[3*anch[0]+j];
        frctmp.map[j+3] = sttspace->map[j];
      }
      if (compfrc == 1) {
       VectorBondR(&crdtmp, &frctmp, sysUV, 0, 1, theta,
                   cfsinp->ops[i].Krst, cfsinp->ops[i].targets,
                   cfsinp->ops[i].fbhw, cfsinp->ops[i].quadwin, nactive);
      }
      else {
       VectorBondRU(&crdtmp, sysUV, 0, 1, theta,
                    cfsinp->ops[i].Krst, cfsinp->ops[i].targets,
                    cfsinp->ops[i].fbhw, cfsinp->ops[i].quadwin, nactive);
      }
    }
    else if (cfsinp->ops[i].order == 2) {
      if (mobile[anch[0]] + mobile[anch[1]] == 0) {
        continue;
      }
      if (compfrc == 1) {
        VectorBondR(crdpop, frcpop, sysUV, anch[0], anch[1], theta,
                    cfsinp->ops[i].Krst, cfsinp->ops[i].targets,
                    cfsinp->ops[i].fbhw, cfsinp->ops[i].quadwin, nactive);
      }
      else {
        VectorBondRU(crdpop, sysUV, anch[0], anch[1], theta,
                     cfsinp->ops[i].Krst, cfsinp->ops[i].targets,
                     cfsinp->ops[i].fbhw, cfsinp->ops[i].quadwin, nactive);
      }
    }
    else if (cfsinp->ops[i].order == 3) {
      if (mobile[anch[0]] + mobile[anch[1]] + mobile[anch[2]] == 0) {
        continue;
      }
      if (compfrc == 1) {
        VectorAngleR(crdpop, frcpop, sysUV, anch[0], anch[1], anch[2],
                     theta, cfsinp->ops[i].Krst, cfsinp->ops[i].targets,
                     cfsinp->ops[i].fbhw, cfsinp->ops[i].quadwin, nactive);
      }
      else {
        VectorAngleRU(crdpop, sysUV, anch[0], anch[1], anch[2], theta,
                      cfsinp->ops[i].Krst, cfsinp->ops[i].targets,
                      cfsinp->ops[i].fbhw, cfsinp->ops[i].quadwin, nactive);
      }
    }
    else if (cfsinp->ops[i].order == 4) {
      if (mobile[anch[0]] + mobile[anch[1]] +
          mobile[anch[2]] + mobile[anch[3]] == 0) {
        continue;
      }
      ReflectDVec(Phase.data, cfsinp->ops[i].targets, cfsinp->count);
      if (compfrc == 1) {
        VectorTorsionR(crdpop, frcpop, anch[0], anch[1], anch[2], anch[3],
                       theta, sysUV, 1, &cfsinp->ops[i].Krst, &Phase,
                       cfsinp->ops[i].fbhw, cfsinp->ops[i].quadwin, nactive);
      }
      else {
        VectorTorsionRU(crdpop, anch[0], anch[1], anch[2], anch[3], theta,
                        sysUV, 1, &cfsinp->ops[i].Krst, &Phase,
                        cfsinp->ops[i].fbhw, cfsinp->ops[i].quadwin, nactive);
      }
    }
  }

  // Loop over all path restraints and apply them.
#if 0
  for (i = 0; i < cfsinp->npaths; i++) {

  }
#endif

  // Free allocated memory
  free(theta);
  DestroyDmat(&Phase);
  if (tmpex == 1) {
    free(frctmp.map);
    free(crdtmp.map);
  }
}

//-----------------------------------------------------------------------------
// CalcNBScaling: calculate the non-bonded interaction scaling (between zero,
//                nothing for an exclusion, and one, full).  This switches
//                between use of a fast table lookup for small molecules and a
//                much slower exclusion list check for larger ones.  Each
//                method draws upon a custom-built data structure.
//
// Arguments:
//   atmi:     the index of the ith atom
//   atmj:     the index of the jth atom
//   exclmeth: flag to decide which method determines the scaling factor, and
//             which scaling factor to return (0 = matrix method, 1 = list
//             method)
//   excl:     the matrix of exclusions (NULL if exclmeth == 1)
//   allex:    the list of exclusions (NULL if exclmeth == 0)
//   qqscl:    electrostatic scaling factor (returned)
//   ljscl:    Lennard-Jones scaling factor (returned)
//-----------------------------------------------------------------------------
static void CalcNBScaling(int atmi, int atmj, int exclmeth, dmat *excl,
                          excltab* allex, double *qqscl, double *ljscl)
{
  int i;

  // Elecrostatics
  if (exclmeth == 0) {
    *qqscl = excl->map[atmi][atmj];
    *ljscl = excl->map[atmj][atmi];
  }
  else {

    // Access the exclusions for the ith atom and see if j is within
    // the range of atoms that could be partly or wholly excluded.
    *qqscl = 1.0;
    *ljscl = 1.0;
    if (atmj >= allex[atmi].llim && atmj < allex[atmi].hlim) {
      for (i = 0; i < allex[atmi].nexcl; i++) {
        if (allex[atmi].nblist[i] == atmj) {
          *qqscl = allex[atmi].qqval[i];
          *ljscl = allex[atmi].ljval[i];
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// VectorNonbondContrib: a vectorized computation of non-bonded interactions.
//                       This assumes all-to-all interactions and makes use of
//                       a matrix to denote exclusions, borrowed from the
//                       parameter fitting routines.
//
// Arguments:
//   tp:       system topology
//   mobile:   array of flags to state whether each atom is movable or not
//   crdpop:   matrix of coordinates, giving X1, Y1, Z1, X2, ..., Zn on
//             successive rows
//   frcpop:   matrix of forces, being accumulated before this function is
//             called and accepting contributions from this function
//   excl:     exclusions and 1:4 attenuations matrix
//   sysUV:    vector of energy accumulators for each conformation
//   compfrc:  flag to have forces computed if set to 1 (otherwise only the
//             energy will be computed)
//   nactive:  the number of configurations actively being manipulated
//-----------------------------------------------------------------------------
static void VectorNonbondContrib(prmtop *tp, int* mobile, dmat *crdpop,
                                 dmat *frcpop, int exclmeth, dmat *excl,
                                 excltab* allex, Energy* sysUV, int compfrc,
                                 int nactive)
{
  int i, j, k, ailjt, ajljt;
  double qi, qij, ljAcoef, ljBcoef, dljAcoef, dljBcoef, qqscl, ljscl;
  double dx, dy, dz, invr, invr2, invr4, invr8, fmag;
  double *xcrdi, *xcrdj, *ycrdi, *ycrdj, *zcrdi, *zcrdj;
  double *xfrci, *xfrcj, *yfrci, *yfrcj, *zfrci, *zfrcj;

  // The familiar nested loop over all atoms, with an
  // additional inner loop for all configurations
  for (i = 0; i < tp->natom-1; i++) {
    qi = tp->Charges[i] * BIOQ;
    ailjt = tp->LJIdx[i];
    xcrdi = crdpop->map[3*i];
    ycrdi = crdpop->map[3*i+1];
    zcrdi = crdpop->map[3*i+2];
    xfrci = frcpop->map[3*i];
    yfrci = frcpop->map[3*i+1];
    zfrci = frcpop->map[3*i+2];
    for (j = i+1; j < tp->natom; j++) {

      // Check that at least one of the atoms is movable
      if (mobile[i] + mobile[j] == 0) {
        continue;
      }

      // Set pointers and prepare to work
      CalcNBScaling(i, j, exclmeth, excl, allex, &qqscl, &ljscl);
      qij = qqscl * qi * tp->Charges[j];
      ajljt = 2*tp->LJIdx[j];
      xcrdj = crdpop->map[3*j];
      ycrdj = crdpop->map[3*j+1];
      zcrdj = crdpop->map[3*j+2];
      xfrcj = frcpop->map[3*j];
      yfrcj = frcpop->map[3*j+1];
      zfrcj = frcpop->map[3*j+2];
      if (ailjt >= 0 && ajljt >= 0) {
        ljAcoef = ljscl * tp->LJutab.map[ailjt][ajljt];
        ljBcoef = ljscl * tp->LJutab.map[ailjt][ajljt+1];
        dljAcoef = ljscl * tp->LJftab.map[ailjt][ajljt];
        dljBcoef = ljscl * tp->LJftab.map[ailjt][ajljt+1];
      }
      else {
        ljAcoef = 0.0;
        ljBcoef = 0.0;
        dljAcoef = 0.0;
        dljBcoef = 0.0;
      }

      // Are forces needed, or just energies?
      if (compfrc == 1) {
        for (k = 0; k < nactive; k++) {
          dx = xcrdj[k] - xcrdi[k];
          dy = ycrdj[k] - ycrdi[k];
          dz = zcrdj[k] - zcrdi[k];
          invr2 = 1.0/(dx*dx + dy*dy + dz*dz);
          invr = sqrt(invr2);
          invr4 = invr2*invr2;
          invr8 = invr4*invr4;
          sysUV[k].elec += qij*invr;
          sysUV[k].vdw12 += ljAcoef*invr8*invr4;
          sysUV[k].vdw6 += ljBcoef*invr2*invr4;
          fmag = -qij*invr*invr2 + invr8*(dljAcoef*invr4*invr2 + dljBcoef);
          xfrci[k] += fmag*dx;
          yfrci[k] += fmag*dy;
          zfrci[k] += fmag*dz;
          xfrcj[k] -= fmag*dx;
          yfrcj[k] -= fmag*dy;
          zfrcj[k] -= fmag*dz;
        }
      }
      else {
        for (k = 0; k < nactive; k++) {
          dx = xcrdj[k] - xcrdi[k];
          dy = ycrdj[k] - ycrdi[k];
          dz = zcrdj[k] - zcrdi[k];
          invr2 = 1.0/(dx*dx + dy*dy + dz*dz);
          invr = sqrt(invr2);
          invr4 = invr2*invr2;
          sysUV[k].elec += qij*invr;
          sysUV[k].vdw12 += ljAcoef*invr4*invr4*invr4;
          sysUV[k].vdw6 += ljBcoef*invr2*invr4;
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// VectorZeroEnergies: a function that operates selectively on an array of
//                     energy trackers to initialize the relevant components
//                     of the molecular mechanics energy.  The Energy data
//                     structure (see TrajectoryDS.h) was designed for actual
//                     simulations of solitary, evolving configurations of a
//                     system, not energy minimizations of many configurations.
//                     It is, rather, co-opted for this purpose, so this
//                     function helps insulate the rest of the development
//                     from that fact.
//
// Arguments:
//   sysUV:        the array of energy trackers
//   nactive:      the number of configurations, or elements in sysUV
//-----------------------------------------------------------------------------
static void VectorZeroEnergies(Energy* sysUV, int nactive)
{
  int i;

  for (i = 0; i < nactive; i++) {
    sysUV[i].bond = 0.0;
    sysUV[i].angl = 0.0;
    sysUV[i].dihe = 0.0;
    sysUV[i].cmap = 0.0;
    sysUV[i].elec = 0.0;
    sysUV[i].relec = 0.0;
    sysUV[i].vdw6 = 0.0;
    sysUV[i].vdw12 = 0.0;
  }
}

//-----------------------------------------------------------------------------
// VectorSumEnergies: like VectorZeroEnergies, this insulates the rest of the
//                    development from the consequences of co-opting the
//                    Energy data structure for vectorized energy minimization.
//                    It merely sums the relevant components of the energy.
//
// Arguments:
//   sysUV:        the array of energy trackers
//   nactive:      the number of configurations, or elements in sysUV
//-----------------------------------------------------------------------------
static void VectorSumEnergies(Energy* sysUV, int nactive)
{
  int i;

  for (i = 0; i < nactive; i++) {
    sysUV[i].etot = sysUV[i].bond + sysUV[i].angl + sysUV[i].dihe +
                    sysUV[i].cmap + sysUV[i].elec + sysUV[i].relec +
                    sysUV[i].vdw12 + sysUV[i].vdw6;
  }
}

//-----------------------------------------------------------------------------
// CheckConstrainedBonds: check bonds that have been constrained by RATTLE, to
//                        verify that they have been held to their proper
//                        lengths.
//
// Arguments:
//   cfsinp:    the configuration sampling input data, contains the file
//              format, base name, and extension to write
//   tp:        the topology to work with
//   crdpop:    the coordinates of configurations found
//-----------------------------------------------------------------------------
static void CheckConstrainedBonds(configs *cfsinp, prmtop *tp, dmat *crdpop,
                                  ConvStat* cvs, FILE *outp)
{
  int i, j, k, ididx, nnode, npass, ncnst, atma, atmb;
  int *myCC;
  int* atomids;
  double rmsd, delta, r2, rx, ry, rz, l0;
  double* xpos;
  
  atomids = (int*)malloc(tp->RattleGrpMax * sizeof(int));
  xpos = (double*)malloc(3 * tp->RattleGrpMax * sizeof(double));

  rmsd = 0.0;
  for (i = 0; i < tp->natom; i++) {
    if (tp->SHL[i].exe != 2) {
      continue;
    }

    // Get information about the constraint group
    myCC = tp->SHL[i].blist;
    ididx = 3*myCC[0] + 2;
    nnode = myCC[ididx - 1];
    for (j = 0; j < nnode; j++) {
      atomids[j] = myCC[ididx + j];
    }
    for (j = 0; j < cfsinp->count; j++) {
      if (cvs[j].pass == 0) {
        continue;
      }

      // Cache reference coordinates
      for (k = 0; k < nnode; k++) {
        xpos[3*k    ] = crdpop->map[3*atomids[k]    ][j];
        xpos[3*k + 1] = crdpop->map[3*atomids[k] + 1][j];
        xpos[3*k + 2] = crdpop->map[3*atomids[k] + 2][j];
      }

      // Loop over each constraint and measure the distance precision
      for (k = 0; k < myCC[0]; k++) {
        atma = myCC[3*k + 1];
        atmb = myCC[3*k + 2];
        rx = xpos[3*atmb    ] - xpos[3*atma    ];
        ry = xpos[3*atmb + 1] - xpos[3*atma + 1];
        rz = xpos[3*atmb + 2] - xpos[3*atma + 2];
        r2 = rx*rx + ry*ry + rz*rz;
        l0 = myCC[3*k + 3] / 1.0e8;
        delta = l0 - sqrt(r2);
        rmsd += delta * delta;
      }
    }
  }
  free(atomids);
  free(xpos);
  // Normalize the deviations
  npass = 0;
  ncnst = 0;
  for (i = 0; i < cfsinp->count; i++) {
    npass += cvs[i].pass;
  }
  for (i = 0; i < tp->natom; i++) {
    if (tp->SHL[i].exe == 2) {
      ncnst += tp->SHL[i].blist[0];
    }
  }
  rmsd /= (double)ncnst * (double)npass;
  rmsd = sqrt(rmsd);
  fprintf(outp, " - Root mean-squared deviation of\n"
          "   constrained bond lengths:       %14.6e\n", rmsd);
}

//-----------------------------------------------------------------------------
// VectorLineMin: this will implement line minimization for every member of an
//                array of configurations using the initial coordinates,
//                forces, and step sizes (specific to each configuration).
//                The general strategy is to make one move, check the step size
//                to make sure it's reasonable, make a second move, check the
//                step size and direction of the move again, and finally make a
//                third move.  Energies computed for each configuration after
//                each move will determine what to do with the step size, and
//                the initial energy plus the energies after each of the three
//                moves will provide four points along the line for cubic
//                spline interpolation.  The minimum of the cubic spline
//                function (within the range defined by the four points) will
//                then define the ultimate move, and so long as those
//                configurations have lower energy they will replace the
//                original coordinates.  If a configuration's energy was not
//                reduced, its step size will be.  If a configuration's energy
//                was improved, its step size may be increased based on the
//                location of the minimum relative to the initial step size.
//
// Arguments:
//   cfsinp:     configuration sampling input data (holds restraint operations)
//   tp:         the system topology
//   exclmeth:   flag to indicate the method of detecting exclusions
//   excl:       matrix of exclusions
//   allex:      list of exclusions
//   crdpop:     coordinates of the population of configurations
//   frcpop:     forces acting on the population of configurations
//   sttspace:   scratch space for coordinate and force calculations when
//               positional restraints are in effect
//   step:       the initial step to take for each configuration (the energies
//               of trial configurations will be sampled at a number of points
//               based on crdpop plus multiples of step and frcpop)
//   nactive:    the number of active configurations (use the column count of
//               crdpop or frcpop to get the total number of configurations) 
//   rattleOn:   flag to have RATTLE constrained bonds
//-----------------------------------------------------------------------------
static void VectorLineMin(configs *cfsinp, prmtop *tp, int exclmeth,
                          dmat *excl, excltab* allex, dmat *crdpop,
                          dmat *frcpop, dmat *sttspace, dmat *crdtrial,
                          Energy* sysUV, Energy* trialUV, dmat *ecurve,
                          dmat *step, int nactive, int rattleOn)
{
  int i, j;
  double stepval, espline, minspline, maxstep, minstep;
  double *dtmp, *dtm2p, *dtm3p, *step1, *step2, *step3;
  dmat A;
  double b[4];
  
  // Load up the initial energies
  for (i = 0; i < nactive; i++) {
    ecurve->map[0][i] = sysUV[i].etot;
  }
  step1 = step->map[1];
  step2 = step->map[2];
  step3 = step->map[3];
  
  // Compute coordinates for each configuration after its first trial move,
  // then decide whether the step for that configuration was too big or too
  // small.  Make adjustments to the step size as needed.
  for (i = 0; i < 3*tp->natom; i++) {
    dtmp = crdtrial->map[i];
    dtm2p = crdpop->map[i];
    dtm3p = frcpop->map[i];
    for (j = 0; j < nactive; j++) {
      dtmp[j] = dtm2p[j] + step1[j]*dtm3p[j];
    }
  }

  // If requested, fix the lengths of bonds to hydrogen
  if (rattleOn == 1) {
    VectorRattle(tp, cfsinp, crdtrial, crdpop, nactive);
  }

  // Adjust extra points
  VectorPlaceXpt(tp, crdtrial, nactive);

  // These energy computations will not change the forces present in frcpop.
  // It's just passed in as a placeholder.
  VectorZeroEnergies(trialUV, nactive);
  VectorBondContrib(tp, cfsinp->movable, crdtrial, frcpop, trialUV, 0,
                    nactive);
  VectorApplyRestraints(cfsinp, crdtrial, frcpop, sttspace, trialUV, 0,
                        nactive);
  VectorNonbondContrib(tp, cfsinp->movable, crdtrial, frcpop, exclmeth, excl,
                       allex, trialUV, 0, nactive);
  VectorSumEnergies(trialUV, nactive);
  for (i = 0; i < nactive; i++) {
    ecurve->map[1][i] = trialUV[i].etot;

    // Decide whether the initial move was too large
    if (ecurve->map[1][i] > ecurve->map[0][i]) {
      step2[i] = 0.5*step1[i];
    }
    else {
      step2[i] = 2.0*step1[i];
    }
  }

  // Perform the next move and compute the energies
  for (i = 0; i < 3*tp->natom; i++) {
    dtmp = crdtrial->map[i];
    dtm2p = crdpop->map[i];
    dtm3p = frcpop->map[i];
    for (j = 0; j < nactive; j++) {
      dtmp[j] = dtm2p[j] + step2[j]*dtm3p[j];
    }
  }
  if (rattleOn == 1) {
    VectorRattle(tp, cfsinp, crdtrial, crdpop, nactive);
  }
  VectorPlaceXpt(tp, crdtrial, nactive);
  VectorZeroEnergies(trialUV, nactive);
  VectorBondContrib(tp, cfsinp->movable, crdtrial, frcpop, trialUV, 0,
                    nactive);
  VectorApplyRestraints(cfsinp, crdtrial, frcpop, sttspace, trialUV, 0,
                        nactive);
  VectorNonbondContrib(tp, cfsinp->movable, crdtrial, frcpop, exclmeth, excl,
                       allex, trialUV, 0, nactive);
  VectorSumEnergies(trialUV, nactive);
  for (i = 0; i < nactive; i++) {
    ecurve->map[2][i] = trialUV[i].etot;

    // Decide whether the second move was too large,
    // and what should be done for the third move.
    if (step2[i] < step1[i]) {

      // The first step was a failure, and the second step is too.
      if (ecurve->map[2][i] > ecurve->map[0][i]) {
        step3[i] = 0.25*step1[i];
      }

      // The second step improved the energy even though the first did not.
      else {
        step3[i] = 0.75*step1[i];
      }
    }
    else {

      // The first step was a success, but the second was not.
      if (ecurve->map[2][i] > ecurve->map[1][i]) {
        step3[i] = 1.5*step1[i];
      }

      // The second step was just as much a success as the first
      else {
        step3[i] = 3.0*step1[i];
      }
    }
  }

  // Perform the final move and compute the energies
  for (i = 0; i < 3*tp->natom; i++) {
    dtmp = crdtrial->map[i];
    dtm2p = crdpop->map[i];
    dtm3p = frcpop->map[i];
    for (j = 0; j < nactive; j++) {
      dtmp[j] = dtm2p[j] + step3[j]*dtm3p[j];
    }
  }
  if (rattleOn == 1) {
    VectorRattle(tp, cfsinp, crdtrial, crdpop, nactive);
  }
  VectorPlaceXpt(tp, crdtrial, nactive);
  VectorZeroEnergies(trialUV, nactive);
  VectorBondContrib(tp, cfsinp->movable, crdtrial, frcpop, trialUV, 0,
                    nactive);
  VectorApplyRestraints(cfsinp, crdtrial, frcpop, sttspace, trialUV, 0,
                        nactive);
  VectorNonbondContrib(tp, cfsinp->movable, crdtrial, frcpop, exclmeth, excl,
                       allex, trialUV, 0, nactive);
  VectorSumEnergies(trialUV, nactive);
  for (i = 0; i < nactive; i++) {
    ecurve->map[3][i] = trialUV[i].etot;
  }
  
  // Compute cubic splines for each configuration's energy profile, and
  // interpolate to find the minimum value within the available range.
  // The ecurve matrix gives the energies, and the step matrix gives the
  // position along the x-axis in each line minimization.
  A = CreateDmat(4, 4, 0);
  for (i = 0; i < nactive; i++) {

    // Construct and solve the matrix equation
    maxstep = 0.0;
    for (j = 0; j < 4; j++) {
      stepval = step->map[j][i];
      if (stepval > maxstep) {
        maxstep = stepval;
      }
      A.map[j][0] = stepval * stepval * stepval;
      A.map[j][1] = stepval * stepval;
      A.map[j][2] = stepval;
      A.map[j][3] = 1.0;
      b[j] = ecurve->map[j][i];
    }
    AxbQRRxc(A, b, 0);
    BackSub(A, b);

    // Use the solution to find the minimum value within the range
    minspline = b[3];
    minstep = 0.0;
    for (j = 0; j <= 200; j++) {
      stepval = 0.005*maxstep*j;
      espline = stepval*(stepval*(b[0]*stepval + b[1]) + b[2]) + b[3];
      if (espline < minspline) {
        minspline = espline;
        minstep = stepval;
      }
    }

    // Store the best step in the energy matrix (the ith column
    // has served its purpose for this line minimization)
    ecurve->map[0][i] = minstep;
  }
  DestroyDmat(&A);
  
  // Take the best steps for each configuration
  dtm3p = ecurve->map[0];
  for (i = 0; i < 3*tp->natom; i++) {
    dtmp = (rattleOn == 1) ? crdtrial->map[i] : crdpop->map[i];
    dtm2p = frcpop->map[i];
    for (j = 0; j < nactive; j++) {
      dtmp[j] += dtm2p[j]*dtm3p[j];
    }
  }
  if (rattleOn == 1) {
    VectorRattle(tp, cfsinp, crdtrial, crdpop, nactive);
    for (i = 0; i < 3*tp->natom; i++) {
      dtmp  = crdpop->map[i];
      dtm2p = crdtrial->map[i];
      for (j = 0; j < nactive; j++) {
        dtmp[j] = dtm2p[j];
      }
    }
  }
  VectorPlaceXpt(tp, crdpop, nactive);
  
  // Update the initial step sizes if the line minimization
  // was very successful, or if it couldn't do very much
  dtmp = step->map[1];
  dtm2p = step->map[2];
  dtm3p = step->map[3];
  for (i = 0; i < nactive; i++) {
    if (dtm3p[i] > 1.6*dtmp[i] && dtmp[i] < 0.1) {
      dtmp[i] *= 1.25;
    }
    else if (dtm3p[i] < dtmp[i]) {
      dtmp[i] *= 0.8;
    }
  }
}

//-----------------------------------------------------------------------------
// VectorConjugateGradient: compute the conjugate gradient moves from a matrix
//                          whose columns are arrays of forces on different
//                          configurations of a structure.
//
// Arguments:
//   tp:          the system topology
//   frcpop:      the population of force vectors
//   gg,dgg,gam:  scratch arrays pre-allocated and re-used for speed
//   nactive:     the number of active configurations
//-----------------------------------------------------------------------------
static void VectorConjugateGradient(prmtop *tp, dmat *frcpop, dmat *gdir,
                                    dmat *htemp, double* gg, double* dgg,
                                    double* gam, int nactive, int niter)
{
  int i, j;
  double *dtmp, *dtm2p, *dtm3p;

  // The implementation here may be a bit hard to follow.
  SetDVec(gg, nactive, 0.0);
  SetDVec(dgg, nactive, 0.0);
  for (i = 0; i < 3*tp->natom; i++) {
    dtmp = gdir->map[i];
    dtm2p = frcpop->map[i];
    for (j = 0; j < nactive; j++) {
      gg[j] += dtmp[j]*dtmp[j];
      dgg[j] += (dtm2p[j]-dtmp[j])*dtm2p[j];
    }
  }
  if (niter == 0) {
    SetDVec(gam, nactive, 0.0);
  }
  else {
    for (i = 0; i < nactive; i++) {
      gam[i] = dgg[i]/gg[i];
    }
  }
  for (i = 0; i < 3*tp->natom; i++) {
    dtmp = gdir->map[i];
    dtm2p = frcpop->map[i];
    dtm3p = htemp->map[i];
    for (j = 0; j < nactive; j++) {
      dtmp[j] = dtm2p[j];
      dtm3p[j] = dtmp[j] + gam[j]*dtm3p[j];
      dtm2p[j] = dtm3p[j];
    }
  }
}

//-----------------------------------------------------------------------------
// VectorDeleteFrozenAtomForces: even though force calculations are only done
//                               for terms that have at least one movable atom,
//                               forces will still accumulate on the frozen
//                               atoms to some degree.  This routine will erase
//                               those forces.
// 
// Arguments:
//   cfsinp:  configuration sampling input data
//   tp:      topology for the system
//   frcpop:  forces acting on the population of configurations
//-----------------------------------------------------------------------------
static void VectorDeleteFrozenAtomForces(configs *cfsinp, prmtop *tp,
                                         dmat *frcpop)
{
  int i;

  // Bail right out if there are no frozen atoms
  if (cfsinp->allmove == 1) {
    return;
  }

  // Loop over all atoms and delete forces on frozen ones
  for (i = 0; i < tp->natom; i++) {
    if (cfsinp->movable[i] == 0) {
      SetDVec(frcpop->map[3*i], cfsinp->count, 0.0);
      SetDVec(frcpop->map[3*i+1], cfsinp->count, 0.0);
      SetDVec(frcpop->map[3*i+2], cfsinp->count, 0.0);
    }
  }
}

//-----------------------------------------------------------------------------
// GenerateConfigurations: use conjugate-gradient energy minimization to form
//                         new configurations from the starting states.  The
//                         given restraints will differentiate or alter each
//                         structure.
//
// Arguments:
//   cfsinp:    configuration sampling input data
//   tp:        topology for the system (helps to determine whether a given
//              file is compatible as a coordinate set for an initial
//              configuration)
//   crd:       matrix of coordinates, giving X1, Y1, Z1, X2, ..., Zn on
//              successive rows
//   sysUV:     vector of energy accumulators for each conformation
//   cvs:       convergence statistics for the run
//   rattleOn:   flag to have RATTLE constrained bonds
//-----------------------------------------------------------------------------
void GenerateConfigurations(configs *cfsinp, prmtop *tp, dmat *crdpop,
                            Energy* sysUV, Energy* isysUV, ConvStat* cvs,
                            int rattleOn)
{
  int i, j, k, m, nconf, nactive, ndata, exclmeth, tmpnops;
  int *anch, *tmovable;
  int* indices;
  int* converged;
  double fconv2, fmag2, dval;
  double *dtmp, *dtm2p, *dtm3p;
  double* magfrc;
  double* gg;
  double* dgg;
  double* gam;
  double* totalE;
  dmat frcpop, crdtrial, crdxpose, excl, step;
  dmat ecurve, gdir, htemp, Phase, U, sttspace;
  cfigop *tmpops;
  bondcomm *thisbc;
  anglcomm *thisac;
  Energy tmpnrg;
  Energy* trialUV;
  Energy* frozenUV;
  ConvStat tmpcv;
  ConvStat* cvscopy;
  excltab* allex;

  // Compute the exclusions matrix
  if (tp->natom <= cfsinp->atomlimit) {
    excl = CompExclMatrix(tp, NULL);
    exclmeth = 0;
  }
  else {
    allex = MakeExclList(tp);
    exclmeth = 1;
  }

  // Allocate memory for forces and scratch computations
  nconf = cfsinp->count;
  frcpop = CreateDmat(3*tp->natom, nconf, 0);
  crdtrial = CreateDmat(3*tp->natom, nconf, 0);
  magfrc = (double*)malloc(nconf*sizeof(double));
  converged = (int*)malloc(nconf*sizeof(int));
  trialUV = (Energy*)malloc(nconf*sizeof(Energy));
  cvscopy = (ConvStat*)malloc(nconf*sizeof(ConvStat));
  sttspace = CreateDmat(25, nconf, 0);
  gg = (double*)malloc(nconf*sizeof(double));
  dgg = (double*)malloc(nconf*sizeof(double));
  gam = (double*)malloc(nconf*sizeof(double));
  gdir = CreateDmat(3*tp->natom, nconf, 0);
  htemp = CreateDmat(3*tp->natom, nconf, 0);
  step = CreateDmat(4, nconf, 0);
  SetDVec(step.map[1], nconf, cfsinp->step0);
  ecurve = CreateDmat(4, nconf, 0);
  ndata = 3*tp->natom*nconf;

  // Unpack minimization instructions from the configs input control data
  fconv2 = cfsinp->fconv * cfsinp->fconv;

  // Track the indices of configurations, even if they get re-arranged
  // when some become energy-minimized and others need more steps.
  indices = CountUp(nconf);
  nactive = nconf;

  // Initialize the convergence statistics
  for (i = 0; i < nconf; i++) {
    cvs[i].nstep = -1;
    cvs[i].pass = 1;
  }

  // Compute the energies of frozen atoms in each configuration
  if (cfsinp->allmove == 0) {

    // Compute the total energy of each configuration; hold on to the
    // actual set of moving atoms while the official pointer gets reset.
    tmovable = cfsinp->movable;
    cfsinp->movable = (int*)malloc(tp->natom*sizeof(int));
    SetIVec(cfsinp->movable, tp->natom, 1);
    frozenUV = (Energy*)malloc(nconf*sizeof(Energy));
    VectorZeroEnergies(frozenUV, nactive);
    SetDVec(frcpop.data, ndata, 0.0);
    VectorBondContrib(tp, cfsinp->movable, crdpop, &frcpop, frozenUV, 0,
                      nactive);
    VectorApplyRestraints(cfsinp, crdpop, &frcpop, &sttspace, frozenUV, 0,
                          nactive);
    VectorNonbondContrib(tp, cfsinp->movable, crdpop, &frcpop, exclmeth, &excl,
                         allex, frozenUV, 0, nactive);
    VectorCommuteXptForces(tp, crdpop, &frcpop, nactive);
    VectorSumEnergies(frozenUV, nactive);

    // Set the official pointer back to the correct list of moving atoms and
    // recompute the energies, now excluding interactions within the subset
    // of frozen atoms.  Subtract this energy from the previous result to 
    // get the energy due only to interactions among the frozen atoms.
    free(cfsinp->movable);
    cfsinp->movable = tmovable;
    VectorZeroEnergies(trialUV, nactive);
    SetDVec(frcpop.data, ndata, 0.0);
    VectorBondContrib(tp, cfsinp->movable, crdpop, &frcpop, trialUV, 0,
                      nactive);
    VectorApplyRestraints(cfsinp, crdpop, &frcpop, &sttspace, trialUV, 0,
                          nactive);
    VectorNonbondContrib(tp, cfsinp->movable, crdpop, &frcpop, exclmeth, &excl,
                         allex, trialUV, 0, nactive);
    VectorCommuteXptForces(tp, crdpop, &frcpop, nactive);
    VectorSumEnergies(trialUV, nactive);
    for (i = 0; i < nconf; i++) {
      frozenUV[i].bond -= trialUV[i].bond;
      frozenUV[i].angl -= trialUV[i].angl;
      frozenUV[i].dihe -= trialUV[i].dihe;
      frozenUV[i].cmap -= trialUV[i].cmap;
      frozenUV[i].elec -= trialUV[i].elec;
      frozenUV[i].vdw6 -= trialUV[i].vdw6;
      frozenUV[i].vdw12 -= trialUV[i].vdw12;
      frozenUV[i].relec -= trialUV[i].relec;
    }
  }

  // If requested, fix the lengths of bonds to hydrogen
  if (rattleOn == 1) {
    CopyDmat(&crdtrial, crdpop, 1);
    VectorRattle(tp, cfsinp, crdpop, &crdtrial, nactive);
  }

  // If there are no cycles, st some energies to avoid initialization failures
  if (cfsinp->maxcyc == 0) {
    VectorZeroEnergies(sysUV, nactive);
    SetDVec(frcpop.data, ndata, 0.0);
    VectorBondContrib(tp, cfsinp->movable, crdpop, &frcpop, sysUV, 1, nactive);
    VectorApplyRestraints(cfsinp, crdpop, &frcpop, &sttspace, sysUV, 1,
                          nactive);
    VectorNonbondContrib(tp, cfsinp->movable, crdpop, &frcpop, exclmeth, &excl,
                         allex, sysUV, 1, nactive);
    VectorCommuteXptForces(tp, crdpop, &frcpop, nactive);
    VectorDeleteFrozenAtomForces(cfsinp, tp, &frcpop);
    VectorSumEnergies(sysUV, nactive);
    for (j = 0; j < nconf; j++) {
      isysUV[j] = sysUV[j];
      cvs[j].lastE = sysUV[j].etot - sysUV[j].relec;
      cvs[j].lastER = sysUV[j].relec;
    }
  }
  
  // Generate configurations
  for (i = 0; i < cfsinp->maxcyc; i++) {
    
    // Alert the user as to progress
    if (cfsinp->verbose == 1 && i % 10 == 0) {
      fprintf(stderr, "\rConfigs >> Minimization step %5d of %5d (%4d active ",
              i, cfsinp->maxcyc, nactive);
      if (nactive > 1) {
        fprintf(stderr, "configurations)");
      }
      else {
        fprintf(stderr, "configuration) ");
      }
    }

    // Initialize and then compute the energies and forces.
    // It is legal, here, to zero forces on inactive configurations (using
    // ndata to span the entire matrix of 3*tp->natom x nconf numbers) as
    // these configurations will not be further manipulated.
    VectorZeroEnergies(sysUV, nactive);
    SetDVec(frcpop.data, ndata, 0.0);
    VectorBondContrib(tp, cfsinp->movable, crdpop, &frcpop, sysUV, 1, nactive);
    VectorApplyRestraints(cfsinp, crdpop, &frcpop, &sttspace, sysUV, 1,
                          nactive);
    VectorNonbondContrib(tp, cfsinp->movable, crdpop, &frcpop, exclmeth, &excl,
                         allex, sysUV, 1, nactive);
    VectorCommuteXptForces(tp, crdpop, &frcpop, nactive);
    VectorDeleteFrozenAtomForces(cfsinp, tp, &frcpop);
    VectorSumEnergies(sysUV, nactive);
    if (i == 0) {
      for (j = 0; j < nconf; j++) {
        isysUV[j] = sysUV[j];
        cvs[j].lastE = sysUV[j].etot - sysUV[j].relec;
        cvs[j].lastER = sysUV[j].relec;
      }
    }
    else {
      for (j = 0; j < nactive; j++) {
        cvs[j].lastdE = (sysUV[j].etot - sysUV[j].relec) - cvs[j].lastE;
        cvs[j].lastdER = sysUV[j].relec - cvs[j].lastER;
        cvs[j].lastE = sysUV[j].etot - sysUV[j].relec;
        cvs[j].lastER = sysUV[j].relec;
      }
    }

    // Check for convergence of forces in any configurations
    SetIVec(converged, nactive, 1);
    for (j = 0; j < nactive; j++) {
      cvs[j].maxF = 0.0;
    }
    for (j = 0; j < tp->natom; j++) {
      dtmp = frcpop.map[3*j];
      dtm2p = frcpop.map[3*j+1];
      dtm3p = frcpop.map[3*j+2];
      for (k = 0; k < nactive; k++) {
        fmag2 = dtmp[k]*dtmp[k] + dtm2p[k]*dtm2p[k] + dtm3p[k]*dtm3p[k];
        if (converged[k] == 1 && fmag2 > fconv2) {
          converged[k] = 0;
        }
        if (fmag2 > cvs[k].maxF) {
          cvs[k].maxF = fmag2;
          cvs[k].maxFatom = j;
        }
      }
    }
    for (j = 0; j < nactive; j++) {
      cvs[j].maxF = sqrt(cvs[j].maxF);
    }

    // Check for vanishing step sizes in any configurations
    dtmp = step.map[1];
    for (j = 0; j < nactive; j++) {
      if (dtmp[j] < cfsinp->stepconv) {
        converged[j] = 1;
      }
    }

    // Swap out converged configurations
    for (j = 0; j < nactive; j++) {
      if (converged[j] == 1) {

        // The jth configuration has reached convergence and should be taken
        // out of consideration.  Find the highest numbered non-converged
        // configuration so the two can be swapped.
        cvs[j].nstep = i;
        while (nactive > 0 && converged[nactive-1] == 1 && nactive > j) {
          nactive--;
        }
        if (j == nactive || nactive == 0) {

          // The last active configuration has reach convergence.
          break;
        }

        // If we are still here, swap the highest numbered non-converged
        // configuration with the lowest-numbered converged one.  The target
        // values of any restraints have to get swapped as well, as do the
        // entries in the step size and memory matrices for minimization.
        converged[j] = 0;
        converged[nactive-1] = 1;
        SWAP(indices[j], indices[nactive-1], k);
        SWAP(sysUV[j], sysUV[nactive-1], tmpnrg);
        if (cfsinp->allmove == 0) {
          SWAP(frozenUV[j], frozenUV[nactive-1], tmpnrg);
        }
        SWAP(cvs[j], cvs[nactive-1], tmpcv);
        SWAP(step.map[1][j], step.map[1][nactive-1], dval);
        SWAP(gg[j], gg[nactive-1], dval);
        SWAP(dgg[j], dgg[nactive-1], dval);
        SWAP(gam[j], gam[nactive-1], dval);
        for (k = 0; k < 3*tp->natom; k++) {
          SWAP(crdpop->map[k][j], crdpop->map[k][nactive-1], dval);
          SWAP(frcpop.map[k][j], frcpop.map[k][nactive-1], dval);
          SWAP(gdir.map[k][j], gdir.map[k][nactive-1], dval);
          SWAP(htemp.map[k][j], htemp.map[k][nactive-1], dval);
        }
        for (k = 0; k < cfsinp->nops; k++) {
          SWAP(cfsinp->ops[k].targets[j], cfsinp->ops[k].targets[nactive-1],
               dval);
        }
      }
    }
    if (nactive == 0) {
      break;
    }

    // Steepest descent or conjugate gradient optimization
    if (i >= cfsinp->ncyc) {
      VectorConjugateGradient(tp, &frcpop, &gdir, &htemp, gg, dgg, gam,
                              nactive, i-cfsinp->ncyc);
    }

    // Normalize the forces on each (remaining) active configuration
    SetDVec(magfrc, nactive, 0.0);
    for (j = 0; j < 3*tp->natom; j++) {
      dtmp = frcpop.map[j];
      for (k = 0; k < nactive; k++) {
        magfrc[k] += dtmp[k]*dtmp[k];
      }
    }
    for (j = 0; j < nactive; j++) {
      magfrc[j] = 1.0/sqrt(magfrc[j]);
    }
    for (j = 0; j < 3*tp->natom; j++) {
      dtmp = frcpop.map[j];
      for (k = 0; k < nactive; k++) {
        dtmp[k] *= magfrc[k];
      }
    }

    // Line minimization along the direction chosen
    VectorLineMin(cfsinp, tp, exclmeth, &excl, allex, crdpop, &frcpop,
                  &sttspace, &crdtrial, sysUV, trialUV, &ecurve, &step,
                  nactive, rattleOn);
    
    // If the minimization didn't fully converge, then we need
    // to recompute the final energies before bailing out.
    if (i == cfsinp->maxcyc - 1) {
      VectorZeroEnergies(sysUV, nactive);
      SetDVec(frcpop.data, ndata, 0.0);
      VectorBondContrib(tp, cfsinp->movable, crdpop, &frcpop, sysUV, 1,
                        nactive);
      VectorApplyRestraints(cfsinp, crdpop, &frcpop, &sttspace, sysUV, 1,
                            nactive);
      VectorNonbondContrib(tp, cfsinp->movable, crdpop, &frcpop, exclmeth,
                           &excl, allex, sysUV, 1, nactive);
      VectorCommuteXptForces(tp, crdpop, &frcpop, nactive);
      VectorDeleteFrozenAtomForces(cfsinp, tp, &frcpop);
      VectorSumEnergies(sysUV, nactive);
      for (j = 0; j < nactive; j++) {
        cvs[j].lastdE = (sysUV[j].etot - sysUV[j].relec) - cvs[j].lastE;
        cvs[j].lastdER = sysUV[j].relec - cvs[j].lastER;
        cvs[j].lastE = sysUV[j].etot - sysUV[j].relec;
        cvs[j].lastER = sysUV[j].relec;
      }
    }
  }

  // Final report on this minimization run
  if (cfsinp->verbose == 1 && nactive == 0) {
    fprintf(stderr, "\rConfigs >> Minimization step %5d of %5d (   0 "
            "active configurations)", i, cfsinp->maxcyc);
  }

  // Reorder the list of configurations (it got shuffled as configurations
  // were moved to the back after reaching convergence).  Start with
  // coordinates (crdpop), using crdtrial as scratch space.  Then do energies
  // (sysUV and frozenUV), each time using trialUV as a scratch space.  The
  // convergence statistics can also be reordered in one of these loops.
  // Finally, reorder the restraint targets so that the final coordinates of
  // each configuration actually reflect the restraints that built it.
  for (i = 0; i < 3*tp->natom; i++) {
    dtmp = crdtrial.map[i];
    dtm2p = crdpop->map[i];
    for (j = 0; j < nconf; j++) {
      dtmp[indices[j]] = dtm2p[j];
    }
  }
  ReflectDVec(crdpop->data, crdtrial.data, ndata);
  for (i = 0; i < nconf; i++) {
    trialUV[indices[i]] = sysUV[i];
    cvscopy[indices[i]] = cvs[i];
  }
  for (i = 0; i < nconf; i++) {
    sysUV[i] = trialUV[i];
    cvs[i] = cvscopy[i];
  }
  if (cfsinp->allmove == 0) {
    for (i = 0; i < nconf; i++) {
      trialUV[indices[i]] = frozenUV[i];
    }
    for (i = 0; i < nconf; i++) {
      frozenUV[i] = trialUV[i];
    }
  }
  for (i = 0; i < cfsinp->nops; i++) {
    for (j = 0; j < nconf; j++) {
      gg[indices[j]] = cfsinp->ops[i].targets[j];
    }
    ReflectDVec(cfsinp->ops[i].targets, gg, nconf);
    if (cfsinp->ops[i].order == 1) {
      for (j = 0; j < nconf; j++) {
        gg[indices[j]] = cfsinp->ops[i].refcrd.map[0][j];
      }
      ReflectDVec(cfsinp->ops[i].refcrd.map[0], gg, nconf);
      for (j = 0; j < nconf; j++) {
        gg[indices[j]] = cfsinp->ops[i].refcrd.map[1][j];
      }
      ReflectDVec(cfsinp->ops[i].refcrd.map[1], gg, nconf);
      for (j = 0; j < nconf; j++) {
        gg[indices[j]] = cfsinp->ops[i].refcrd.map[2][j];
      }
      ReflectDVec(cfsinp->ops[i].refcrd.map[2], gg, nconf);
    }
  }

  // Check the various types of strain on each configuration: if there is
  // too much strain on any one bond, any one angle, or too much overall
  // restraint energy, mark the configuration.

  // Check the bonds: re-use arrays allocated for
  // CG minimization to store equilibria
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->BLC[i].nbond; j++) {
      VectorZeroEnergies(trialUV, nconf);
      thisbc = &tp->BLC[i].BC[j];
      if (cfsinp->movable[thisbc->a] + cfsinp->movable[thisbc->b] == 0) {
        continue;
      }
      SetDVec(gg, nconf, tp->BParam[thisbc->t].l0);
      VectorBond(crdpop, 0, trialUV, thisbc->a, thisbc->b, dgg,
                 tp->BParam[thisbc->t].K, tp->BParam[thisbc->t].Kpull,
                 tp->BParam[thisbc->t].Kpress, gg,
                 tp->BParam[thisbc->t].lpull0, tp->BParam[thisbc->t].lpress0,
                 nconf);
      for (k = 0; k < nconf; k++) {
        if (trialUV[k].bond > cvs[k].maxbstrn) {
          cvs[k].maxbstrn = trialUV[k].bond;
          cvs[k].bstrnloc[0] = i;
          cvs[k].bstrnloc[1] = j;
          cvs[k].bstrnloc[2] = round(100000.0 * dgg[k]);
        }
      }
    }
  }

  // Check the angles: re-use arrays allocated for CG minimization
  // to store equilibria.  Skip angles with no moving atoms.
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->ALC[i].nangl; j++) {
      VectorZeroEnergies(trialUV, nconf);
      thisac = &tp->ALC[i].AC[j];
      if (cfsinp->movable[thisac->a] + cfsinp->movable[i] +
          cfsinp->movable[thisac->c] == 0) {
        continue;
      }
      SetDVec(gg, nconf, tp->AParam[thisac->t].th0);
      VectorAngle(crdpop, 0, trialUV, thisac->a, i, thisac->c, dgg,
                  tp->AParam[thisac->t].K, gg, nconf);
      for (k = 0; k < nconf; k++) {
        if (trialUV[k].angl > cvs[k].maxastrn) {
          cvs[k].maxastrn = trialUV[k].angl;
          cvs[k].astrnloc[0] = i;
          cvs[k].astrnloc[1] = j;
          cvs[k].astrnloc[2] = round(1000000.0 * dgg[k]);
        }
      }
    }
  }

  // Check the restraints: find the restraint in each configuration
  // that is furthest from equilibrium.  More re-use of dgg.  Again,
  // skip any restraints that apply to no mobile atoms.
  tmpnops = cfsinp->nops;
  tmpops = cfsinp->ops;
  cfsinp->nops = 1;
  for (i = 0; i < tmpnops; i++) {
    VectorZeroEnergies(trialUV, nconf);
    VectorApplyRestraints(cfsinp, crdpop, &frcpop, &sttspace, trialUV, 0,
                          nconf);
    for (j = 0; j < nconf; j++) {
      if (trialUV[j].relec > cvs[j].maxrst) {
        cvs[j].maxrst = trialUV[j].relec;
        cvs[j].maxrloc = i;
      }
    }
    if (i < tmpnops - 1) {
      cfsinp->ops = &cfsinp->ops[1];
    }
  }
  cfsinp->nops = tmpnops;
  cfsinp->ops = tmpops;
  
  // Record whether each configuration now fails
  for (i = 0; i < nconf; i++) {
    if (cvs[i].maxbstrn > cfsinp->maxbstrn) {
      cvs[i].pass = 0;
    }
    if (cvs[i].maxastrn > cfsinp->maxastrn) {
      cvs[i].pass = 0;
    }
    if (sysUV[i].relec > cfsinp->strainlim) {
      cvs[i].pass = 0;
    }
  }

  // With the tests now complete on the moveable atoms in the system,
  // the energies due to interactions among the frozen atoms can be
  // reincorporated.  In this manner, restraints applied strictly
  // between frozen atoms (which can do nothing to the structure),
  // cannot inadvertently cause energy minimizations which otherwise
  // worked to fail the sanity checks.
  if (cfsinp->allmove == 0) {
    for (i = 0; i < nconf; i++) {
      sysUV[i].bond += frozenUV[i].bond;
      sysUV[i].angl += frozenUV[i].angl;
      sysUV[i].dihe += frozenUV[i].dihe;
      sysUV[i].cmap += frozenUV[i].cmap;
      sysUV[i].elec += frozenUV[i].elec;
      sysUV[i].vdw6 += frozenUV[i].vdw6;
      sysUV[i].vdw12 += frozenUV[i].vdw12;
      sysUV[i].relec += frozenUV[i].relec;
    }
  }

  // Final check: filter conformations that are too similar in
  // positional RMSD if they are also too similar in energy.
  if (cfsinp->simEtol > 0.0) {
    if (cfsinp->rmsdtol > 0.0) {
      crdxpose = CreateDmat(nconf, 3*tp->natom, 0);
    }
    U = CreateDmat(3, 3, 0);
    totalE = (double*)malloc(nconf*sizeof(double));
    for (i = 0; i < nconf; i++) {
      totalE[i] = sysUV[i].bond + sysUV[i].angl + sysUV[i].dihe +
                  sysUV[i].cmap + sysUV[i].elec + sysUV[i].vdw6 +
                  sysUV[i].vdw12;
    }
    for (i = 0; i < 3*tp->natom; i++) {
      if (cfsinp->rmsdtol > 0.0) {
        for (j = 0; j < nconf; j++) {
          crdxpose.map[j][i] = crdpop->map[i][j];
        }
      }
    }
    for (i = 0; i < nconf-1; i++) {
      if (cvs[i].pass == 0) {
        continue;
      }
      for (j = i+1; j < nconf; j++) {
        if (cvs[j].pass == 0) {
          continue;
        }
        if (cfsinp->rmsdtol > 0.0) {
          QuatAlign(crdxpose.map[i], crdxpose.map[j],
                    tp->natom, tp->Masses, 1, &U);
          if (VecRMSD(crdxpose.map[i], crdxpose.map[j],
                      tp->natom) < cfsinp->rmsdtol &&
              fabs(totalE[i] - totalE[j]) < cfsinp->simEtol) {
            cvs[j].pass = 0;
          }
        }
        else if (fabs(totalE[i] - totalE[j]) < cfsinp->simEtol) {
          cvs[j].pass = 0;
        }
      }
    }
    free(totalE);
  }

  // Free allocated memory
  free(magfrc);
  free(converged);
  free(indices);
  free(gg);
  free(dgg);
  free(gam);
  free(trialUV);
  free(cvscopy);
  DestroyDmat(&frcpop);
  DestroyDmat(&crdtrial);
  DestroyDmat(&step);
  DestroyDmat(&ecurve);
  if (tp->natom <= cfsinp->atomlimit) {
    DestroyDmat(&excl);
  }
  else {
    for (i = 0; i < tp->natom; i++) {
      DestroyExcltab(&allex[i]);
    }
    free(allex);
  }
  DestroyDmat(&gdir);
  DestroyDmat(&htemp);
  DestroyDmat(&sttspace);
  if (cfsinp->allmove == 0) {
    free(frozenUV);
  }
}

//-----------------------------------------------------------------------------
// ShuffleConfigurations: shuffle the configurations with previously optimized
//                        results in order to re-attempt the same optimizations
//                        with different initial conditions.
//
// Arguments:
//   cfsinp:    the configuration sampling input data, contains the file
//              format, base name, and extension to write
//   tp:        the topology to work with
//   crdpop:    the coordinates of all configurations optimized with respect
//              to their restraint energies from a previous optimization
//   crdtrial:  scratch space for new configurations
//   tj:        trajectory control input data, used for the rndcon attribute
//              to drive the pseudo-random number generator
//-----------------------------------------------------------------------------
void ShuffleConfigurations(configs *cfsinp, prmtop *tp, dmat *crdpop,
                           dmat *crdtrial, trajcon *tj)
{
  int i, j, k, iswap, iselect, icount, nconf, nbelow;
  double dval, minE, threshE;
  double* restE;
  dmat Tsave, sttspace;
  Energy* trialUV;

  // Shortcuts
  nconf = cfsinp->count;

  // Random reshuffling, with replacement
  if (strcmp(cfsinp->shuffletype, "bootstrap") == 0) {
    for (i = 0; i < nconf; i++) {
      iswap = nconf*ran2(&tj->rndcon);
      for (j = 0; j < 3*tp->natom; j++) {
        crdtrial->map[j][i] = crdpop->map[j][iswap];
      }
    }
  }

  // Random reshuffling, without replacement
  if (strcmp(cfsinp->shuffletype, "jackknife") == 0) { 
    CopyDmat(crdtrial, crdpop, 1);
    for (i = 0; i < nconf; i++) {
      iswap = nconf*ran2(&tj->rndcon);
      for (j = 0; j < 3*tp->natom; j++) {
        SWAP(crdtrial->map[j][i], crdtrial->map[j][iswap], dval);
      }
    }
  }

  // Proximity reshuffling in the restrained degrees of freedom.
  // This is more complex than the other resampling plans.
  if (strcmp(cfsinp->shuffletype, "proximity") == 0) {

    // First, loop over all configurations and rescore EVERY configuration
    // with respect to the restraints of each of them.  Doing so requires
    // that we overwrite the restraint targets.  Save the originals first.
    Tsave = CreateDmat(cfsinp->nops, nconf, 0);
    sttspace = CreateDmat(25, nconf, 0);
    trialUV = (Energy*)malloc(nconf*sizeof(Energy));
    restE = (double*)malloc(nconf*sizeof(double));
    for (i = 0; i < cfsinp->nops; i++) {
      ReflectDVec(Tsave.map[i], cfsinp->ops[i].targets, nconf);
    }
    CopyDmat(crdtrial, crdpop, 1);
    for (i = 0; i < nconf; i++) {

      // Set restraints 
      for (j = 0; j < cfsinp->nops; j++) {
        SetDVec(cfsinp->ops[j].targets, nconf, Tsave.map[j][i]);
      }
      
      // Rescore all configurations
      VectorZeroEnergies(trialUV, nconf);
      VectorApplyRestraints(cfsinp, crdpop, NULL, &sttspace, trialUV, 0,
                            nconf);

      // Binary seach to find the threshold at which the right amount of
      // configurations have restraint energies below the desired value.
      for (j = 0; j < nconf; j++) {
        restE[j] = trialUV[j].relec;
      }
      nbelow = 0;
      minE = DExtreme(restE, nconf, 0);
      threshE = minE + cfsinp->proximity;
      for (j = 0; j < nconf; j++) {
        if (restE[j] < threshE) {
          nbelow++;
        }
      }

      // Select one of the admissible configurations to swap into the trial
      iselect = ran2(&tj->rndcon)*nbelow;
      icount = 0;
      for (j = 0; j < nconf; j++) {
        if (restE[j] < threshE) {
          if (icount == iselect) {
            for (k = 0; k < 3*tp->natom; k++) {
              crdtrial->map[k][i] = crdpop->map[k][j];
            }
            break;
          }
          icount++;
        }
      }
    }

    // Restore the restraints to their original settings
    for (i = 0; i < cfsinp->nops; i++) {
      ReflectDVec(cfsinp->ops[i].targets, Tsave.map[i], nconf);
    }

    // Free allocated memory
    DestroyDmat(&Tsave);
    DestroyDmat(&sttspace);
    free(restE);
  }
}

//-----------------------------------------------------------------------------
// MergeConfigurations: after reshuffling, better solutions to each restraint
//                      setup may have been found.  Merge those back into the
//                      set of configurations to report.
//
// Arguments:
//   cfsinp:    the configuration sampling input data, contains the file
//              format, base name, and extension to write
//   tp:        the topology to work with
//   crdpop:    the coordinates of all configurations optimized with respect
//              to their restraint energies from a previous optimization
//   crdtrial:  coordinates of new configurations
//   sysUV:     energies of the current set of configurations
//   trialUV:   energies of the newly solved trial configurations
//   cvs:       convergence statistics of the existing set of configurations
//   trialcvs:  convergence statistics of the trial configurations
//-----------------------------------------------------------------------------
void MergeConfigurations(configs *cfsinp, prmtop *tp, dmat *crdpop,
                         dmat *crdtrial, Energy* sysUV, Energy* trialUV,
                         Energy* rtlUV, Energy* rtlTrialUV, ConvStat* cvs,
                         ConvStat* trialcvs, ConvStat* rattlecvs,
                         ConvStat* rtlTrialcvs)
{
  int i, j, direc;
  double dE, eAorig, eAtrial;

  direc = (strcmp(cfsinp->shfdir, "up") == 0) ? 1 : 0;
  for (i = 0; i < cfsinp->count; i++) {
    if (cfsinp->rattle == 1) {

      // Take the non-restraint energy of the configuration
      eAorig = rtlUV[i].etot - rtlUV[i].relec;
      eAtrial = rtlTrialUV[i].etot - rtlTrialUV[i].relec;

      // If the original configuration does not pass sanity
      // checks and the new one does, take the new one.
      // If the new configuration is of lower (or higher)
      // non-restraint energy than the existing one, AND
      // the restraint energy is not going higher, take it
      // IF the new configuration passes sanity checks.
      if ((rattlecvs[i].pass == 0 && rtlTrialcvs[i].pass == 1) ||
          (((direc == 0 && eAtrial < eAorig - cfsinp->Ereplace) ||
           (direc == 1 && eAtrial > eAorig + cfsinp->Ereplace)) &&
           rtlTrialUV[i].relec < rtlUV[i].relec && rtlTrialcvs[i].pass == 1)) {
        for (j = 0; j < 3*tp->natom; j++) {
          crdpop->map[j][i] = crdtrial->map[j][i];
        }
        sysUV[i] = trialUV[i];
        rtlUV[i] = rtlTrialUV[i];
        cvs[i] = trialcvs[i];
        cvs[i].dEshuffle = cvs[i].dEshuffle +
                           (trialUV[i].etot - sysUV[i].etot);
        rattlecvs[i] = rtlTrialcvs[i];
        rattlecvs[i].dEshuffle = rattlecvs[i].dEshuffle +
                                 (rtlTrialUV[i].etot - rtlUV[i].etot);
        cfsinp->nreopt += 1;
      }
    }
    else {

      // Take the non-restraint energy of the configuration
      eAorig = sysUV[i].etot - sysUV[i].relec;
      eAtrial = trialUV[i].etot - trialUV[i].relec;

      // If the original configuration does not pass sanity
      // checks and the new one does, take the new one.
      // If the new configuration is of lower (or higher)
      // non-restraint energy than the existing one, AND
      // the restraint energy is not going higher, take it
      // IF the new configuration passes sanity checks.
      if ((cvs[i].pass == 0 && trialcvs[i].pass == 1) ||
          (((direc == 0 && eAtrial < eAorig - cfsinp->Ereplace) ||
           (direc == 1 && eAtrial > eAorig + cfsinp->Ereplace)) &&
           trialUV[i].relec < sysUV[i].relec && trialcvs[i].pass == 1)) {
        for (j = 0; j < 3*tp->natom; j++) {
          crdpop->map[j][i] = crdtrial->map[j][i];
        }
        dE = cvs[i].dEshuffle + (trialUV[i].etot - sysUV[i].etot);
        sysUV[i] = trialUV[i];
        cvs[i] = trialcvs[i];
        cvs[i].dEshuffle = dE;
        cfsinp->nreopt += 1;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// ArbitrateConfigurations: search through the convergence data, including
//                          high-energy bonds and angles, and determine how
//                          many passing configurations there are in the final
//                          set.
//
// Arguments:
//   cfsinp:    the configuration sampling input data, contains the file
//              format, base name, and extension to write
//   sysUV:     energies of the current set of configurations
//   cvs:       convergence statistics of the existing set of configurations
//-----------------------------------------------------------------------------
static void ArbitrateConfigurations(configs *cfsinp, Energy *sysUV,
                                    ConvStat* cvs)
{
  int i;

  // Check the bond strain, angle strain, and restraint energy
  // to see if each configuration will ultimately be printed.
  cfsinp->nbelim = 0;
  cfsinp->naelim = 0;
  cfsinp->nrelim = 0;
  cfsinp->npass = 0;
  for (i = 0; i < cfsinp->count; i++) {
    if (cvs[i].pass == 1) {
      cfsinp->npass += 1;
      continue;
    }
    if (cvs[i].maxbstrn > cfsinp->maxbstrn) {
      cfsinp->nbelim += 1;
    }
    if (cvs[i].maxastrn > cfsinp->maxastrn) {
      cfsinp->naelim += 1;
    }
    if (sysUV[i].relec > cfsinp->strainlim) {
      cfsinp->nrelim += 1;
    }
  }
}

//-----------------------------------------------------------------------------
// VerticalToCoord: transpose a column of the array containing the vectorized
//                  coordinates into a (pre-allocated) coord data structure.
//
// Arguments:
//   crdpop:    the coordinates of all configurations
//   ncol:      the column number of the configuration to print
//   tc:        coordinate data structure made to hold this configuration
//-----------------------------------------------------------------------------
void VerticalToCoord(dmat *crdpop, int ncol, coord *tc)
{
  int i;

  for (i = 0; i < crdpop->row; i++) {
    tc->loc[i] = crdpop->map[i][ncol];
  }
  for (i = 0; i < 3; i++) {
    tc->gdim[i] = 200.0;
    tc->gdim[i+3] = 0.5*PI;
  }
}

//-----------------------------------------------------------------------------
// GetOperationLimits: find the four limits of a flat-bottom harmonic well (the
//                     points at which the potential goes from linear to
//                     harmonic, harmonic to flat, flat to harmonic, and then
//                     harmonic to linar again) based on the nature of the
//                     restraint and the locally stored parameters.
//
// Arguments:
//   myop:       the operation to translate into sander/pmemd nmropt limits
//   cidx:       index of the configuration of interest--each restraint is
//               centered on a specific point for a given configuration
//   r{1,2,3,4}: values corresponding to sander/pmemd nmropt r1, r2, r3, and r4
//-----------------------------------------------------------------------------
static void GetOperationLimits(cfigop *myop, int cidx, double *r1, double *r2,
                               double *r3, double *r4)
{
  double cpt;

  cpt = myop->targets[cidx];
  *r1 = cpt - myop->fbhw - myop->quadwin;
  *r2 = cpt - myop->fbhw;
  *r3 = cpt + myop->fbhw;
  *r4 = cpt + myop->fbhw + myop->quadwin;
  if (myop->order >= 3) {
    *r1 *= 180.0/PI;
    *r2 *= 180.0/PI;
    *r3 *= 180.0/PI;
    *r4 *= 180.0/PI;
  }
}

//-----------------------------------------------------------------------------
// Config2Pdb: write a configuration to a PDB file.
//
// Arguments:
//   cfsinp:     configuration sampling input data
//   tc:         the coordinates of the configuration to write
//   tj:         trajectory control data (needed to tell whether it is all
//               right to overwrite pre-existing files)
//   tp:         system topology
//   sysUV:      array of energies for all configurations of this system (the
//               one of interest will be plucked out)
//   cidx:       index of the configuration to print
//   stlidx:     style of the configuration to print (index into ostyle)
//-----------------------------------------------------------------------------
void Config2Pdb(configs *cfsinp, coord *tc, trajcon *tj, prmtop *tp,
                Energy* sysUV, int cidx, int stlidx)
{
  int i, j, ires;
  double r1, r2, r3, r4;
  char* fname;
  FILE *foutp;

  // Concatenate the input file name
  fname = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(fname, "%s%d.%s", cfsinp->outbase.map[stlidx], cidx+1,
          cfsinp->outsuff.map[stlidx]);
  foutp = FOpenSafe(fname, tj->OverwriteOutput);
  fprintf(foutp, "HEADER  Printed by mdgx\n");
  fprintf(foutp, "REMARK%4d\n", 4);
  fprintf(foutp, "REMARK%4d Topology used to create this configuration:\n", 4);
  fprintf(foutp, "REMARK%4d   %s\n", 4, tp->source);
  fprintf(foutp, "REMARK%4d\n", 4);
  fprintf(foutp, "REMARK%4d Restraints used to create this configuration:"
          "\n", 4);
  for (i = 0; i < cfsinp->nops; i++) {
    fprintf(foutp, "REMARK%4d &rst\nREMARK%4d  iat=", 4, 4);
    for (j = 0; j < cfsinp->ops[i].order; j++) {
      fprintf(foutp, "%d,", cfsinp->ops[i].anchors[j]+1);
    }
    fprintf(foutp, "\n");
    GetOperationLimits(&cfsinp->ops[i], cidx, &r1, &r2, &r3, &r4);
    fprintf(foutp, "REMARK%4d  r1=%.4lf, r2=%.4lf, r3=%.4lf, r4=%.4lf,\n", 4,
            r1, r2, r3, r4);
    fprintf(foutp, "REMARK%4d  rk2=%.4lf, rk3=%.4lf,\nREMARK%4d &end\n", 4,
            cfsinp->ops[i].Krst, cfsinp->ops[i].Krst, 4);
  }
  fprintf(foutp, "REMARK%4d\n", 4);
  fprintf(foutp, "REMARK%4d Molecular mechanics energy decomposition "
          "(kcal/mol):\n", 4);
  fprintf(foutp, "REMARK%4d   Bond          %12.4lf\n", 4, sysUV[cidx].bond);
  fprintf(foutp, "REMARK%4d   Angle         %12.4lf\n", 4, sysUV[cidx].angl);
  fprintf(foutp, "REMARK%4d   Dihedral      %12.4lf (sum of proper and "
          "improper terms)\n", 4, sysUV[cidx].dihe);
  if (tp->ncmap > 0) {
    fprintf(foutp, "REMARK%4d   CMAP          %12.4lf\n", 4, sysUV[cidx].cmap);
  }
  fprintf(foutp, "REMARK%4d   van-der Waals %12.4lf (includes 1-4 "
          "interactions)\n", 4, sysUV[cidx].vdw6 + sysUV[cidx].vdw12);
  fprintf(foutp, "REMARK%4d   Electrostatic %12.4lf (includes 1-4 "
          "interactions)\n", 4, sysUV[cidx].elec);
  fprintf(foutp, "REMARK%4d   Restraints    %12.4lf\n", 4, sysUV[cidx].relec);
  fprintf(foutp, "REMARK%4d   Total energy  %12.4lf (includes restraint "
          "energy)\n", 4, sysUV[cidx].etot);
  fprintf(foutp, "REMARK%4d\n", 4);
  fprintf(foutp, "CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2lf P 1 1 1\n",
          tc->gdim[0], tc->gdim[1], tc->gdim[2], tc->gdim[3]*180.0/PI,
          tc->gdim[4]*180.0/PI, tc->gdim[5]*180.0/PI);
  for (i = 0; i < tp->natom; i++) {
    ires = LocateResID(tp, i, 0, tp->nres);
    fprintf(foutp, "ATOM %6d %.4s %.4sA%4d    %8.3lf%8.3lf%8.3lf\n",
            i+1, &tp->AtomNames[4*i], &tp->ResNames[4*ires], ires,
            tc->loc[3*i], tc->loc[3*i+1], tc->loc[3*i+2]);
  }
  fprintf(foutp, "END\n");
  fclose(foutp);

  // Free allocated memory
  free(fname);
}

//-----------------------------------------------------------------------------
// Config2QuantumInput: prints a Gaussian input file based on a configuration.
//                      This function is distinct from the printing capability
//                      in the IPolQ module because that is built with much
//                      more complexity, anticipating an external charge field
//                      and a different nomenclature.
//
// Arguments:
//   cfsinp:     configuration sampling input data
//   tc:         coordinates of the configuration to print (passed in because
//               this has already been created)
//   tp:         the topology to work with
//   tj:         trajectory control data structure (for its file overwriting
//               flag)
//   cidx:       index of the configuration to print
//   stlidx:     style of the configuration to print (index into ostyle)
//-----------------------------------------------------------------------------
static void Config2QuantumInput(configs *cfsinp, coord *tc, prmtop *tp,
                                trajcon *tj, int cidx, int stlidx)
{
  int i, j, itotq;
  double totq, znum;
  char tlet[8], optword[8];
  char* fname;
  FILE *foutp;

  // Construct the file name and open the file for writing
  fname = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(fname, "%s%d.%s", cfsinp->outbase.map[stlidx], cidx+1,
          cfsinp->outsuff.map[stlidx]);
  foutp = FOpenSafe(fname, tj->OverwriteOutput);

  // Write the header for the QM calculation
  if (strcmp(cfsinp->ostyle.map[stlidx], "GAUSSIAN") == 0) {
    fprintf(foutp, "--Link1--\n%%nproc=%d\n%%mem=%dMB\n%%chk=%s\n#"
            "# %s/%s maxdisk=16GB, pop=(minimal, MK), 5d,\n\n",
            cfsinp->QMsettings.ncpu, cfsinp->QMsettings.MaxCore,
            cfsinp->QMsettings.checkpoint, cfsinp->QMsettings.theory,
            cfsinp->QMsettings.basis);
  }
  else if (strcmp(cfsinp->ostyle.map[stlidx], "ORCA") == 0) {
    if (cfsinp->QMsettings.ncpu > 1) {
      fprintf(foutp, "! PAL%d\n", cfsinp->QMsettings.ncpu);
    }
    optword[0] = '\0';
    if (cfsinp->QMsettings.opt == 1) {
      sprintf(optword, "Opt");
    }
    fprintf(foutp, "! %s %s TightSCF %s NoKeepInts\n\n",
            cfsinp->QMsettings.theory, cfsinp->QMsettings.basis, optword);
    fprintf(foutp, "%%scf\n  MaxCore %d\nend\n", cfsinp->QMsettings.MaxCore);
    if (strcmp(cfsinp->QMsettings.theory, "mp2") == 0) {
      fprintf(foutp, "%%mp2\n  MaxCore %d\nend\n", cfsinp->QMsettings.MaxCore);
    }
  }

  // Print constraints
  if (cfsinp->QMsettings.opt == 1 && cfsinp->QMsettings.xferRst) {
    if (strcmp(cfsinp->ostyle.map[stlidx], "ORCA") == 0) {
      j = 0;
      for (i = 0; i < cfsinp->nops; i++) {
        j += (cfsinp->ops[i].order == 2 &&
              strncmp(cfsinp->ops[i].label, "AutoHB_", 7) != 0);
	j += (cfsinp->ops[i].order > 2);
      }
      if (j > 0) {
        fprintf(foutp, "%%geom\nConstraints\n");
        for (i = 0; i < cfsinp->nops; i++) {
          if (cfsinp->ops[i].order == 2 &&
              strncmp(cfsinp->ops[i].label, "AutoHB_", 7) != 0) {
            fprintf(foutp, "{ B %3d %3d C }\n", cfsinp->ops[i].anchors[0],
                    cfsinp->ops[i].anchors[1]);
          }
          else if (cfsinp->ops[i].order == 3) {
            fprintf(foutp, "{ A %3d %3d %3d C }\n", cfsinp->ops[i].anchors[0],
                    cfsinp->ops[i].anchors[1], cfsinp->ops[i].anchors[2]);
          }
          else if (cfsinp->ops[i].order == 4) {
            fprintf(foutp, "{ D %3d %3d %3d %3d C }\n",
                    cfsinp->ops[i].anchors[0], cfsinp->ops[i].anchors[1],
                    cfsinp->ops[i].anchors[2], cfsinp->ops[i].anchors[3]);
          }
	}
      }
      fprintf(foutp, "end\nend\n\n");
    }
  }

  // Print the molecular charge and spin multiplicity
  totq = DSum(tp->Charges, tp->natom);
  itotq = (totq < 0.0) ? totq - 0.5 : totq + 0.5;
  if (strcmp(cfsinp->ostyle.map[stlidx], "GAUSSIAN") == 0) {
    fprintf(foutp, "Generated by mdgx\n\n%d   %d\n", itotq,
            cfsinp->QMsettings.spin);
  }
  else if (strcmp(cfsinp->ostyle.map[stlidx], "ORCA") == 0) {
    fprintf(foutp, "*xyz  %d  %d\n", itotq, cfsinp->QMsettings.spin);
  }

  // Print the atomic elements and coordinates
  for (i = 0; i < tp->natom; i++) {
    if (tp->Masses[i] > 1.0e-8) {
      znum = AtomCode(tp->Masses[i], tp->ZNumber[i], tp->Charges[i], tlet);
      fprintf(foutp, "    %-4.4s  %12.7lf %12.7lf %12.7lf\n", tlet,
              tc->loc[3*i], tc->loc[3*i+1], tc->loc[3*i+2]);
    }
  }
  if (strcmp(cfsinp->ostyle.map[stlidx], "GAUSSIAN") == 0) {
    fprintf(foutp, "\n\n");
  }
  else if (strcmp(cfsinp->ostyle.map[stlidx], "ORCA") == 0) {
    fprintf(foutp, "*\n");
  }
  fclose(foutp);
  
  // Free allocated memory
  free(fname);
}

//-----------------------------------------------------------------------------
// WriteRestraintXferFile: write a file to transfer restraints used in the MM
//                         optimization to a QM optimization.  OpenFF format
//                         only.
//
// Arguments:
//   cfsinp:     configuration sampling input data
//   tc:         coordinates of the configuration to print (passed in because
//               this has already been created)
//   tp:         the topology to work with
//   tj:         trajectory control data structure (for its file overwriting
//               flag)
//   cidx:       index of the configuration to print
//   stlidx:     style of the configuration to print (index into ostyle)
//-----------------------------------------------------------------------------
static void WriteRestraintXferFile(configs *cfsinp, coord *tc, prmtop *tp,
                                   trajcon *tj, int cidx, int stlidx)
{
  int i;
  char* fname;
  FILE *foutp;

  // Write the restraints transfer file
  fname = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(fname, "%s%d.%s", cfsinp->outbase.map[stlidx], cidx+1,
          cfsinp->outsuff.map[stlidx]);
  foutp = FOpenSafe(fname, tj->OverwriteOutput);
  for (i = 0; i < cfsinp->nops; i++) {

    // Bond restraint potential
    if (cfsinp->ops[i].order == 2) {
      if (cfsinp->ops[i].fbhw > 1.0e-8) {

        // Split the force into two components
        fprintf(foutp, "fbnd_%dL = "
                "CustomBondForce(\"step(%9.5lf-r)*%9.5lf*(r-%9.5lf)^2\")\n",
                i, cfsinp->ops[i].targets[cidx] - cfsinp->ops[i].fbhw,
                cfsinp->ops[i].Krst,
                cfsinp->ops[i].targets[cidx] - cfsinp->ops[i].fbhw);
        fprintf(foutp, "fbnd_%dL.addBond(%d, %d)\n", i,
                cfsinp->ops[i].anchors[0], cfsinp->ops[i].anchors[1]);
        fprintf(foutp, "system.addForce(fbnd_%dL)\n", i);
        fprintf(foutp, "fbnd_%dH = "
                "CustomBondForce(\"step(r-%9.5lf)*%9.5lf*(r-%9.5lf)^2\")\n",
                i, cfsinp->ops[i].targets[cidx] + cfsinp->ops[i].fbhw,
                cfsinp->ops[i].Krst,
                cfsinp->ops[i].targets[cidx] + cfsinp->ops[i].fbhw);
        fprintf(foutp, "fbnd_%dH.addBond(%d, %d)\n", i, 
                cfsinp->ops[i].anchors[0], cfsinp->ops[i].anchors[1]);
        fprintf(foutp, "system.addForce(fbnd_%dH)\n", i);
      }
      else {

        // A harmonic force is easier to specify
        fprintf(foutp, "fbnd_%d = "
                "CustomBondForce(\"%9.5lf*(r-%9.5lf)^2\")\n", i,
                cfsinp->ops[i].Krst, cfsinp->ops[i].targets[cidx]);
        fprintf(foutp, "fbnd_%d.addBond(%d, %d)\n", i, 
                cfsinp->ops[i].anchors[0], cfsinp->ops[i].anchors[1]);
        fprintf(foutp, "system.addForce(fbnd_%d)\n", i);
      }
    }

    // Angle restraint potential
    else if (cfsinp->ops[i].order == 3) {
      if (cfsinp->ops[i].fbhw > 1.0e-8) {

        // Split the force into two components
        fprintf(foutp, "fang_%dL = "
                "CustomAngleForce(\"step(%9.5lf-theta)*%9.5lf*"
                "(theta-%9.5lf)^2\")\n", i,
                cfsinp->ops[i].targets[cidx] - cfsinp->ops[i].fbhw,
                cfsinp->ops[i].Krst,
                cfsinp->ops[i].targets[cidx] - cfsinp->ops[i].fbhw);
        fprintf(foutp, "fang_%dL.addAngle(%d, %d, %d)\n", i, 
                cfsinp->ops[i].anchors[0], cfsinp->ops[i].anchors[1],
                cfsinp->ops[i].anchors[2]);
        fprintf(foutp, "system.addForce(fang_%dL)\n", i);
        fprintf(foutp, "fang_%dH = "
                "CustomAngleForce(\"step(theta-%9.5lf)*%9.5lf*"
                "(theta-%9.5lf)^2\")\n", i,
                cfsinp->ops[i].targets[cidx] + cfsinp->ops[i].fbhw,
                cfsinp->ops[i].Krst,
                cfsinp->ops[i].targets[cidx] + cfsinp->ops[i].fbhw);
        fprintf(foutp, "fang_%dH.addAngle(%d, %d, %d)\n", i,
                cfsinp->ops[i].anchors[0], cfsinp->ops[i].anchors[1],
                cfsinp->ops[i].anchors[2]);
        fprintf(foutp, "system.addForce(fang_%dH)\n", i);
      }
      else {

        // A harmonic force is easier to specify
        fprintf(foutp, "fang_%d = "
                "CustomAngleForce(\"%9.5lf*(theta-%9.5lf)^2\")\n", i,
                cfsinp->ops[i].Krst, cfsinp->ops[i].targets[cidx]);
        fprintf(foutp, "fang_%d.addAngle(%d, %d, %d)\n", i,
                cfsinp->ops[i].anchors[0], cfsinp->ops[i].anchors[1],
                cfsinp->ops[i].anchors[2]);
        fprintf(foutp, "system.addForce(fang_%d)\n", i);
      }
    }

    // Torsion restraint potential
    else if (cfsinp->ops[i].order == 4) {
      if (cfsinp->ops[i].fbhw > 1.0e-8) {

        // Split the force into two components
        fprintf(foutp, "ftor_%dL = "
                "CustomTorsionForce(\"step(%9.5lf-theta)*%9.5lf*"
                "(theta-%9.5lf)^2\")\n", i,
                cfsinp->ops[i].targets[cidx] - cfsinp->ops[i].fbhw,
                cfsinp->ops[i].Krst,
                cfsinp->ops[i].targets[cidx] - cfsinp->ops[i].fbhw);
        fprintf(foutp, "ftor_%dL.addTorsion(%d, %d, %d, %d)\n", i,
                cfsinp->ops[i].anchors[0], cfsinp->ops[i].anchors[1],
                cfsinp->ops[i].anchors[2], cfsinp->ops[i].anchors[3]);
        fprintf(foutp, "system.addForce(ftor_%dL)\n", i);
        fprintf(foutp, "ftor_%dH = "
                "CustomTorsionForce(\"step(theta-%9.5lf)*%9.5lf*"
                "(theta-%9.5lf)^2\")\n", i,
                cfsinp->ops[i].targets[cidx] + cfsinp->ops[i].fbhw,
                cfsinp->ops[i].Krst,
                cfsinp->ops[i].targets[cidx] + cfsinp->ops[i].fbhw);
        fprintf(foutp, "ftor_%dH.addTorsion(%d, %d, %d, %d)\n", i,
                cfsinp->ops[i].anchors[0], cfsinp->ops[i].anchors[1],
                cfsinp->ops[i].anchors[2], cfsinp->ops[i].anchors[3]);
        fprintf(foutp, "system.addForce(ftor_%dH)\n", i);
      }
      else {

        // A harmonic force is easier to specify
        fprintf(foutp, "ftor_%d = "
                "CustomTorsionForce(\"%9.5lf*(theta-%9.5lf)^2\")\n", i,
                cfsinp->ops[i].Krst, cfsinp->ops[i].targets[cidx]);
        fprintf(foutp, "ftor_%d.addTorsion(%d, %d, %d, %d)\n", i,
                cfsinp->ops[i].anchors[0], cfsinp->ops[i].anchors[1],
                cfsinp->ops[i].anchors[2], cfsinp->ops[i].anchors[3]);
        fprintf(foutp, "system.addForce(ftor_%d)\n", i);
      }
    }
  }
  fclose(foutp);
  
  // Free allocated memory
  free(fname);
}

//-----------------------------------------------------------------------------
// CheckForAllBonds: check to see whether three atoms are all bonded to one.
//                   Return 1 if true, 0 if false.
//
// Arguments:
//   tp:         the topology to work with
//   ai:         the atom that might be the center
//   a(j,k,l):   the atoms that may be bonded
//-----------------------------------------------------------------------------
static int CheckForAllBonds(prmtop *tp, int ai, int aj, int ak, int al)
{
  int i, ntriad;

  ntriad = 0;
  for (i = 0; i < tp->nb1234[ai].n12; i++) {
    if (tp->nb1234[ai].L12[i] == aj || tp->nb1234[ai].L12[i] == ak ||
        tp->nb1234[ai].L12[i] == al) {
      ntriad++;
    }
  }

  return (ntriad == 3);
}

//-----------------------------------------------------------------------------
// PrintResidueAndAtom: print the name of a residue, its number, and the atom
//                      name of a specific topology index in a tight format.
//
// Arguments:
//   tp:         the topology to work with
//   ai:         the atom index
//   ostr:       the output string (pre-allocated to at least 32 bytes please)
//-----------------------------------------------------------------------------
static void PrintResidueAndAtom(prmtop *tp, int ai, char* ostr)
{
  int i, j, ri;

  ri = LocateResID(tp, ai, 0, tp->nres);
  j = 0;
  for (i = 0; i < 4; i++) {
    if (tp->ResNames[4*ri + i] != ' ') {
      ostr[j] = tp->ResNames[4*ri + i];
      j++;
    }
  }
  sprintf(&ostr[j], "(%d)-", ri);
  j = strlen(ostr);
  for (i = 0; i < 4; i++) {
    if (tp->AtomNames[4*ai + i] != ' ') {
      ostr[j] = tp->AtomNames[4*ai + i];
      j++;
    }
  }
  while (j < 13) {
    ostr[j] = ' ';
    j++;
  }
  ostr[13] = '\0';
}

//-----------------------------------------------------------------------------
// WriteGeomeTRICFreeze: write a file containing constraints for critical
//                       heavy-atom torsions that will effectively freeze a
//                       conformation in place, with local optimizations to
//                       bonds and bond angles being the only major
//                       rearrangements, when presented to the GeomeTRIC
//                       wrapper for OpenFF-approved quantum programs.
//
// Arguments:
//   cfsinp:     configuration sampling input data
//   tc:         coordinates of the configuration to print (passed in because
//               this has already been created)
//   tp:         the topology to work with
//   tj:         trajectory control data structure (for its file overwriting
//               flag)
//   cidx:       index of the configuration to print
//   stlidx:     style of the configuration to print (index into ostyle)
//-----------------------------------------------------------------------------
static void WriteGeomeTRICFreeze(configs *cfsinp, coord *tc, prmtop *tp,
                                 trajcon *tj, int cidx, int stlidx)
{
  int i, j, k, ip, iw, ik, totdihe;
  int* isImpCenter;
  int* notViable;
  char ostr[32];
  char* fname;
  FILE *foutp;

  // Make a list of torsions from the toplogy, subject to three rules:
  //   1.) Do not duplicate torsion angle constraints
  //   2.) If, in the torsion angle A-B-C-D, B and C are both improper centers,
  //       do not include the torsion
  //   3.) Do not count a torsion A-B-C-D if either A or D was a terminal atom
  //       (this is how most heavy atoms will be detected)
  isImpCenter = (int*)calloc(tp->natom, sizeof(int));
  totdihe = 0;
  for (i = 0; i < tp->natom; i++) {
    totdihe += tp->HLC[i].ndihe;
  }
  notViable = (int*)calloc(totdihe, sizeof(int));
  k = 0;
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->HLC[i].ndihe; j++) {
      
      // Check to see if all three atoms bond to this one
      if (CheckForAllBonds(tp, i, tp->HLC[i].HC[j].a, tp->HLC[i].HC[j].c,
                           tp->HLC[i].HC[j].d) == 1) {
        isImpCenter[i] = 1;
        notViable[k] = 1;
      }
      else if (CheckForAllBonds(tp, tp->HLC[i].HC[j].a, i, tp->HLC[i].HC[j].c,
                                tp->HLC[i].HC[j].d) == 1) {
        isImpCenter[tp->HLC[i].HC[j].a] = 1;
        notViable[k] = 1;
      }
      else if (CheckForAllBonds(tp, tp->HLC[i].HC[j].c, tp->HLC[i].HC[j].a, i,
                                tp->HLC[i].HC[j].d) == 1) {
        isImpCenter[tp->HLC[i].HC[j].c] = 1;
        notViable[k] = 1;
      }
      else if (CheckForAllBonds(tp, tp->HLC[i].HC[j].d, tp->HLC[i].HC[j].a,
                                tp->HLC[i].HC[j].c, i) == 1) {
        isImpCenter[tp->HLC[i].HC[j].d] = 1;
        notViable[k] = 1;
      }
      k++;
    }
  }

  // Accumulate important torsions
  k = 0;
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->HLC[i].ndihe; j++) {
      if (notViable[k] == 1) {
        k++;
        continue;
      }
      if (isImpCenter[tp->HLC[i].HC[j].c] == 1 && isImpCenter[i] == 1) {
        notViable[k] = 1;
      }
      else if (tp->nb1234[tp->HLC[i].HC[j].a].n12 == 1 ||
          tp->nb1234[tp->HLC[i].HC[j].d].n12 == 1) {
        if (cfsinp->FreezePolarH == 0 ||
            (IsPolarHydrogen(tp, tp->HLC[i].HC[j].a) == 0 &&
             IsPolarHydrogen(tp, tp->HLC[i].HC[j].d) == 0)) {
          notViable[k] = 1;
        }
      }

      // Final check: does this overlap with another torsion
      // that can freeze the same degree of freedom?
      ik = 0;
      for (ip = 0; ip < i; ip++) {
        for (iw = 0; iw < tp->HLC[ip].ndihe; iw++) {
          if (tp->HLC[ip].HC[iw].c == i && tp->HLC[i].HC[j].c == ip &&
              notViable[ik] == 0) {
            notViable[k] = 1;
          }
          else if (tp->HLC[ip].HC[iw].a == tp->HLC[i].HC[j].a &&
                   tp->HLC[ip].HC[iw].d == tp->HLC[i].HC[j].d &&
                   notViable[ik] == 0) {
            notViable[k] = 1;
          }
          else if (tp->HLC[ip].HC[iw].a == tp->HLC[i].HC[j].d &&
                   tp->HLC[ip].HC[iw].d == tp->HLC[i].HC[j].a &&
                   notViable[ik] == 0) {
            notViable[k] = 1;
          }
          ik++;
        }
      }
      for (iw = 0; iw < j; iw++) {
        if (tp->HLC[i].HC[iw].c == tp->HLC[i].HC[j].c && notViable[ik] == 0) {
          notViable[k] = 1;
        }
        ik++;
      }
      k++;
    }
  }
  
  // Write out the dihedrals that can be frozen
  if (ISum(notViable, totdihe) == totdihe) {
    return;
  }
  fname = (char*)malloc(MAXNAME * sizeof(char));
  sprintf(fname, "%s%d.%s", cfsinp->outbase.map[stlidx], cidx + 1,
          cfsinp->outsuff.map[stlidx]);
  foutp = FOpenSafe(fname, tj->OverwriteOutput);
  fprintf(foutp, "$freeze\n");
  k = 0;
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->HLC[i].ndihe; j++) {
      if (notViable[k] == 0) {
        fprintf(foutp, "%4d %4d %4d %4d  # ", tp->HLC[i].HC[j].a + 1, i + 1,
                tp->HLC[i].HC[j].c + 1, tp->HLC[i].HC[j].d + 1);
        PrintResidueAndAtom(tp, tp->HLC[i].HC[j].a, ostr);
        fprintf(foutp, "%s", ostr);
        PrintResidueAndAtom(tp, i, ostr);
        fprintf(foutp, "%s", ostr);
        PrintResidueAndAtom(tp, tp->HLC[i].HC[j].c, ostr);
        fprintf(foutp, "%s", ostr);
        PrintResidueAndAtom(tp, tp->HLC[i].HC[j].d, ostr);
        fprintf(foutp, "%s", ostr);
        fprintf(foutp, "\n");
      }
      k++;
    }
  }
  fclose(foutp);

  // Free allocated memory
  free(isImpCenter);
  free(notViable);
}

//-----------------------------------------------------------------------------
// WriteJSON: write a .json format file for OpenFF QCFractal to read as a 
//            quantum input deck.
//
// Arguments:
//   cfsinp:     configuration sampling input data
//   tc:         coordinates of the configuration to print (passed in because
//               this has already been created)
//   tp:         the topology to work with
//   tj:         trajectory control data structure (for its file overwriting
//               flag)
//   cidx:       index of the configuration to print
//   stlidx:     style of the configuration to print (index into ostyle)
//-----------------------------------------------------------------------------
#if 0
static void WriteJSON(configs *cfsinp, coord *tc, prmtop *tp, trajcon *tj,
                      int stlidx)
{
  int i;
  double dznum;
  char* fname, tlet[8];
  FILE *foutp;
  
  // Construct the file name and open the file for writing
  fname = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(fname, "%s%d.%s", cfsinp->outbase.map[stlidx], cidx+1,
          cfsinp->outsuff.map[stlidx]);
  foutp = FOpenSafe(fname, tj->OverwriteOutput);

  fprintf(foutp, "{\n  \"dataset_name\": \"%s\",\n", cfsinp->dsName);
  fprintf(foutp, "  \"dataset_tagline\": \"Data set created by mdgx\",\n");
  fprintf(foutp, "  \"dataset_type\": \"%s\",\n", cfsinp->dsType);
  fprintf(foutp, "  \"method\": \"%s\",\n", cfsinp->QMSettings.theory);
  fprintf(foutp, "  \"basis\": \"%s\",\n", cfsinp->QMSettings.basis);
  fprintf(foutp, "  \"program\": \"psi4\",\n");
  fprintf(foutp, "  \"maxiter\": 200,\n");
  fprintf(foutp, "  \"driver\": \"gradient\",\n");
  fprintf(foutp, "  \"scf_properties\": [\n    \"dipole\",\n"
	  "    \"quadrupole\",\n    \"wiberg_lowdin_indices\",\n"
	  "    \"mayer_indices\",\n  ],\n");
  fprintf(foutp, "  \"spec_name\": \"default\",\n");
  fprintf(foutp, "  \"spec_description\": \"OpenFF quantum chemistry input\","
	  "\n");
  fprintf(foutp, "  \"priority\": \"normal\",\n");
  fprintf(foutp, "  \"description\": \"mdgx conformational search outputs "
	  "converted to JSON format\",\n");
  fprintf(foutp, "  \"dataset_tags\": [\n    \"openff\"\n  ],\n");
  fprintf(foutp, "  \"compute_tag\": \"openff\",\n");
  fprintf(foutp, "  \"dataset\": {\n");
  fprintf(foutp, "    \"index\": \"");
  for (i = 0; i < tp->nres; i++) {
    for (j = 0; j < 4; j++) {
      if (tp->ResNames[4*i + j] != ' ') {
	fprintf(foutp, "%c", tp->ResNames[4*i + j]);
      }
    }
    if (i < tp->nres - 1) {
      fprintf(foutp, "_");
    }
  }
  fprintf(foutp, "\",\n");
  fprintf(foutp, "    \"initial_molecules\": [\n");
  for (i = 0; i < cfsinp->count; i++) {
    fprintf(foutp, "      {\n");
    fprintf(foutp, "        \"schema_name\": \"qcschema_molecule\",\n");
    fprintf(foutp, "        \"schema_version\": 2,\n");
    fprintf(foutp, "        \"validated\": true,\n");
    fprintf(foutp, "        \"symbols\": [\n");
    for (j = 0; j < tp->natom; j++) {
      if (tp->ZNumber[j] > 0) {
        dznum = AtomCode(0.0, tp->ZNumber[j], tp->Charges[j], tlet);
        fprintf(foutp, "          \"%s\",\n", tlet);
      }
    }
    fprintf(foutp, "        ],\n");
    fprintf(foutp, "        \"geometry\": [\n");
    for (j = 0; j < tp->natom; j++) {
      if (tp->ZNumber[j] > 0) {
        fprintf(foutp, "          %10.5lf,\n", tc->loc[3*j    ]);
        fprintf(foutp, "          %10.5lf,\n", tc->loc[3*j + 1]);
        fprintf(foutp, "          %10.5lf,\n", tc->loc[3*j + 2]);
      }
    }
    fprintf(foutp, "        ],\n");
    fprintf(foutp, "        \"name\": \"%s_%d\",\n",
            cfsinp->outbase.map[stlidx], cidx + 1);
    fprintf(foutp, "        \"molecular_charge\": %9.5lf,\n",
            DSum(tp->Charges, tp->natom));
    fprintf(foutp, "        \"molecular_multiplicity\": %d,\n",
	    
  }
  
  // Free allocated memory
  free(fname);
}
#endif
  
//-----------------------------------------------------------------------------
// PrintConfigurations: print each configuration found in a specified format.
//                      The configurations could be printed to a trajectory,
//                      a series of inpcrd files, PDB files, or a series of
//                      input files for various QM packages.
//
// Arguments:
//   cfsinp:     the configuration sampling input data, contains the file
//               format, base name, and extension to write
//   tp:         the topology to work with
//   crdpop:     the coordinates of configurations found
//   tj:         trajectory control data (for passing information down to
//               the standard trjaectory writing routines, so they can be
//               co-opted for this purpose)
//   sysUV:      energies of the configurations produced by restrained energy
//               minimization
//   cvs:        convergence statistics for the (initial, non-constrained) run
//-----------------------------------------------------------------------------
static void PrintConfigurations(configs *cfsinp, prmtop *tp, dmat *crdpop,
                                trajcon *tj, Energy* sysUV, ConvStat* cvs)
{
  int i, j, npass;
  char pdbname[MAXNAME];
  coord tc;
  cdftrj Acdf;

  // Allocate a temporary coordinates array to hold the results
  tc = CreateCoord(tp->natom);

  // Count the number of configurations that will pass--this is critical
  // for writing trajectories.  Scan through the outputs as well, to see
  // if trajectories need to be written, and prime them if they do.
  npass = 0;
  for (i = 0; i < cfsinp->count; i++) {
    if (cvs[i].maxbstrn <= cfsinp->maxbstrn &&
        cvs[i].maxastrn <= cfsinp->maxastrn &&
        sysUV[i].relec <= cfsinp->strainlim) {
      npass++;
    }
  }
  for (i = 0; i < cfsinp->ostyle.row; i++) {
    if (strcmp(cfsinp->ostyle.map[i], "CRD") == 0) {
      sprintf(tj->trjbase.map[0], "%s", cfsinp->outbase.map[i]);
      sprintf(tj->trjsuff.map[0], ".%s", cfsinp->outsuff.map[i]);
      tj->currstep = 0;
      WriteCrd(NULL, &tc, 1, tj, tp, 0);
    }
    else if (strcmp(cfsinp->ostyle.map[i], "CDF") == 0) {
      sprintf(tj->trjbase.map[0], "%s", cfsinp->outbase.map[i]);
      sprintf(tj->trjsuff.map[0], ".%s", cfsinp->outsuff.map[i]);
      tj->currstep = 0;
      WriteCDF(NULL, &tc, 1, tj, &Acdf, tp, 0);
    }
  }

  // Make sure that, when splicing file names, there will be no
  // attempt to add extra numbers based on uninitialized integers.
  // Also rig this so that every configuration will be written as
  // a "frame" of the trajectory to a single file.
  tj->nfistep = 0;
  tj->irest = 0;
  tj->ntwx = 1;
  tj->ntpr = -1;
  tj->nstep = npass+1;
  npass = 1;
  for (i = 0; i < cfsinp->count; i++) {

    // Check to ensure that this configuration passes muster
    if (cvs[i].pass == 0) {
      continue;
    }
    VerticalToCoord(crdpop, i, &tc);
    for (j = 0; j < cfsinp->ostyle.row; j++) {
      if (strcmp(cfsinp->ostyle.map[j], "INPCRD") == 0) {
        sprintf(tj->rstbase.map[0], "%s%d", cfsinp->outbase.map[j], i+1);
        sprintf(tj->rstsuff.map[0], ".%s", cfsinp->outsuff.map[j]);
        tj->currstep = npass;
        WriteRst(NULL, &tc, tp, tj, 0, -1);
      }
      else if (strcmp(cfsinp->ostyle.map[j], "PDB") == 0) {
        Config2Pdb(cfsinp, &tc, tj, tp, sysUV, i, j);
      }
      else if (strcmp(cfsinp->ostyle.map[j], "CRD") == 0) {
        sprintf(tj->trjbase.map[0], "%s", cfsinp->outbase.map[j]);
        sprintf(tj->trjsuff.map[0], ".%s", cfsinp->outsuff.map[j]);
        tj->currstep = npass;
        WriteCrd(NULL, &tc, 1, tj, tp, 0);
      }
      else if (strcmp(cfsinp->ostyle.map[j], "CDF") == 0) {
        sprintf(tj->trjbase.map[0], "%s", cfsinp->outbase.map[j]);
        sprintf(tj->trjsuff.map[0], ".%s", cfsinp->outsuff.map[j]);
        tj->currstep = npass;
        WriteCDF(NULL, &tc, 1, tj, &Acdf, tp, 0);
      }
      else if (strcmp(cfsinp->ostyle.map[j], "GAUSSIAN") == 0 ||
               strcmp(cfsinp->ostyle.map[j], "ORCA") == 0) {
        Config2QuantumInput(cfsinp, &tc, tp, tj, i, j);
      }
      else if (strcmp(cfsinp->ostyle.map[j], "GAMESS") == 0) {
  
      }
      else if (strcmp(cfsinp->ostyle.map[j], "MOLPRO") == 0) {
  
      }
      else if (strcmp(cfsinp->ostyle.map[j], "QRST") == 0) {
        WriteRestraintXferFile(cfsinp, &tc, tp, tj, i, j);
      }
      else if (strcmp(cfsinp->ostyle.map[j], "FREEZE") == 0) {
        WriteGeomeTRICFreeze(cfsinp, &tc, tp, tj, i, j);
      }
    }
    npass += 1;
  }

  // Free allocated memory
  DestroyCoord(&tc);
}

//-----------------------------------------------------------------------------
// ConfigEStats: a function for printing summary statistics on a vectorized
//               energy minimization run.
//
// Arguments:
//   cfsinp:    the configuration sampling input data, contains the file
//              format, base name, and extension to write
//   sysUV:     energies of the configurations produced by restrained energy
//              minimization
//   outp:      the output file being written
//-----------------------------------------------------------------------------
static void ConfigEStats(configs *cfsinp, Energy* sysUV, Energy* isysUV,
                         ConvStat* cvs, FILE *outp)
{
  int i;
  double* scratch;
  char* msg;

  // General explanation of this section
  HorizontalRule(outp, 0);
  PrintVADesc(0, "(1.)", 4, " ", 1, "Energies of the initial and final "
              "states.  Units are kcal/mol.  A number of statistics could be "
              "of interest here: this table tries to anticipate them.", 74, 0,
              outp);

  // Table header
  fprintf(outp, "\n %-24.24s %-12.12s %-12.12s %-12.12s %-12.12s\n",
          "Statistic", "Mean Value", " Std. Dev. ", "Min. Value",
          "Max. Value");
  fprintf(outp, " %-24.24s %-12.12s %-12.12s %-12.12s %-12.12s\n",
          " ", " (kcal/mol) ", " (kcal/mol) ", " (kcal/mol) ", " (kcal/mol) ");
  fprintf(outp, " ------------------------ ------------ ------------ "
          "------------ ------------\n");

  // Allocate scratch space to do computations
  scratch = (double*)malloc(cfsinp->count*sizeof(double));

  // Initial states
  for (i = 0; i < cfsinp->count; i++) {
    scratch[i] = isysUV[i].etot - isysUV[i].relec;
  }
  fprintf(outp, " %-24.24s %12.4lf %12.4lf %12.4lf %12.4lf\n",
          "Initial E(Model) (a)", DAverage(scratch, cfsinp->count),
          DStDev(scratch, cfsinp->count), DExtreme(scratch, cfsinp->count, 0),
          DExtreme(scratch, cfsinp->count, 1));
  for (i = 0; i < cfsinp->count; i++) {
    scratch[i] = isysUV[i].relec;
  }
  fprintf(outp, " %-24.24s %12.4lf %12.4lf %12.4lf %12.4lf\n",
          "Initial E(Restraint) (b)", DAverage(scratch, cfsinp->count),
          DStDev(scratch, cfsinp->count), DExtreme(scratch, cfsinp->count, 0),
          DExtreme(scratch, cfsinp->count, 1));

  // Final states
  for (i = 0; i < cfsinp->count; i++) {
    scratch[i] = sysUV[i].etot - sysUV[i].relec;
  }
  fprintf(outp, " %-24.24s %12.4lf %12.4lf %12.4lf %12.4lf\n",
          "Final E(Model) (c)", DAverage(scratch, cfsinp->count),
          DStDev(scratch, cfsinp->count), DExtreme(scratch, cfsinp->count, 0),
          DExtreme(scratch, cfsinp->count, 1));
  for (i = 0; i < cfsinp->count; i++) {
    scratch[i] = sysUV[i].relec;
  }
  fprintf(outp, " %-24.24s %12.4lf %12.4lf %12.4lf %12.4lf\n",
          "Final E(Restraint) (d)", DAverage(scratch, cfsinp->count),
          DStDev(scratch, cfsinp->count), DExtreme(scratch, cfsinp->count, 0),
          DExtreme(scratch, cfsinp->count, 1));

  // Deltas
  for (i = 0; i < cfsinp->count; i++) {
    scratch[i] = (sysUV[i].etot - sysUV[i].relec) -
                 (isysUV[i].etot - isysUV[i].relec);
  }
  fprintf(outp, " %-24.24s %12.4lf %12.4lf %12.4lf %12.4lf\n",
          "Delta E(Model) (e)", DAverage(scratch, cfsinp->count),
          DStDev(scratch, cfsinp->count), DExtreme(scratch, cfsinp->count, 0),
          DExtreme(scratch, cfsinp->count, 1));
  for (i = 0; i < cfsinp->count; i++) {
    scratch[i] = sysUV[i].relec - isysUV[i].relec;
  }
  fprintf(outp, " %-24.24s %12.4lf %12.4lf %12.4lf %12.4lf\n",
          "Delta E(Restraint) (f)", DAverage(scratch, cfsinp->count),
          DStDev(scratch, cfsinp->count), DExtreme(scratch, cfsinp->count, 0),
          DExtreme(scratch, cfsinp->count, 1));
  if (cfsinp->reshuffle > 0) {
    for (i = 0; i < cfsinp->count; i++) {
      scratch[i] = cvs[i].dEshuffle;
    }
    fprintf(outp, " %-24.24s %12.4lf %12.4lf %12.4lf %12.4lf\n",
            "Reoptimization (g)", DAverage(scratch, cfsinp->count),
            DStDev(scratch, cfsinp->count),
            DExtreme(scratch, cfsinp->count, 0),
            DExtreme(scratch, cfsinp->count, 1));
  }
  fprintf(outp, "\n");

  // Explain the table
  PrintVADesc(2, "(a)", 3, " ", 1, "The energies of the initial models, "
              "without the restraint energy (this is an indicator of how "
              "strained the configurations were before mdgx began to "
              "manipulate them).", 73, 0, outp);
  PrintVADesc(2, "(b)", 3, " ", 1, "When restraints were first applied to the "
              "initial models, these were the restraint energies.", 73, 0,
              outp);
  PrintVADesc(2, "(c)", 3, " ", 1, "The energies of the final models, again "
              "excluding restraint energy.  This is an indicator of how "
              "strained the molecular configurations became after "
              "manipulation.", 73, 0, outp);
  PrintVADesc(2, "(d)", 3, " ", 1, "Restraint energy needed to get the final "
              "configurations into place.", 73, 0, outp);
  PrintVADesc(2, "(e)", 3, " ", 1, "This is an indicator of the change in "
              "molecular mechanics energy.  E(model) will typically rise, "
              "indicating that structures became distorted in response to "
              "pressure from external restraints.", 73, 0, outp);
  PrintVADesc(2, "(f)", 3, " ", 1, "The change in restraint energy.  "
              "E(restraint) will typically fall as restraints force the "
              "structures into new configurations.", 73, 0, outp);
  if (cfsinp->reshuffle > 0) {
    msg = (char*)malloc(512*sizeof(char));
    if (strcmp(cfsinp->shuffletype, "jackknife") == 0 ||
        strcmp(cfsinp->shuffletype, "bootstrap") == 0) {
      sprintf(msg, "Reoptimization was performed by %s, randomly swapping "
              "solutions to different restraint sets and using them as the "
              "initial states in subsequent energy minimizations.  If this "
              "drove the energy %s relative to the previous solution, the "
              "result was retained.", cfsinp->shuffletype, cfsinp->shfdir);
    }
    if (strcmp(cfsinp->shuffletype, "proximity") == 0) {
      sprintf(msg, "Proximity reoptimization was performed for each state S "
              "by randomly selecting solutions that could be substituted for "
              "S without driving up the restraint energy by more than %.4lf " 
              "kcal/mol.  These substitutes were then put through a new round "
              "of energy minimization, and if the energy of the reoptimized "
              "structures went %s relative to the previous solution S, the "
              "result was kept.", cfsinp->proximity, cfsinp->shfdir);
    }
    PrintVADesc(2, "(g)", 3, " ", 1, msg, 73, 0, outp);
    free(msg);
  }
  HorizontalRule(outp, 1);

  // Free allocated memory
  free(scratch);
}

//-----------------------------------------------------------------------------
// ConvergenceStats: prints summary statistics on the convergence of energy
//                   minimization across all configurations
//
// Arguments:
//   cfsinp:     the configuration sampling input data, contains the file
//               format, base name, and extension to write
//   tp:         the topology to work with
//   crdpop:     the coordinates of finished configurations
//   cvs:        convergence statistics for the run
//   rattlecvs:  convergence statistics for the run
//   outp:       the output file being written
//-----------------------------------------------------------------------------
static void ConvergenceStats(configs *cfsinp, prmtop *tp, dmat *crdpop,
                          ConvStat* cvs, ConvStat* rattlecvs, FILE *outp)
{
  int i, j, k, nsuccess, nfail, avestep, stdstep, probatom, maxcount;
  int* floc;
  double* scratch;
  char* desc;

  // General explanation of this section
  HorizontalRule(outp, 0);
  desc = (char*)malloc(512*sizeof(char));
  sprintf(desc, "Convergence statistics for the run.  This block of "
          "information describes the energy minimization of ALL "
          "configurations, regardless of whether they were printed.  The next "
          "section presents counts of configurations passing the basic sanity "
         "checks.");
  if (cfsinp->reshuffle > 0) {
    sprintf(&desc[strlen(desc)], "  These convergence statistics pertain to "
            "the final round of optimization performed after %d rounds of "
            "reshuffling.", cfsinp->reshuffle);
  }
  if (cfsinp->rattle == 1) {
    sprintf(&desc[strlen(desc)], "  In addition to any restraints specified, "
           "constraints were applied to all bonds invovling hydrogen.  The "
           "RATTLE procedure does not mix with standard line minimization, "
           "even though harmonic restraints do.  Therefore, after each round "
           "of unCONSTRAINED minimization, additional cycles were performed "
           "with RATTLE interspersed.  The radius of convergence is very "
           "small (does not mix with standard line minimization), but the "
           "products are near local optimums and have bond lengths fixed at "
           "equilibrium values.");
  }
  PrintVADesc(0, "(2.)", 4, " ", 1, desc, 74, 0, outp);
  free(desc);

  // Allocate space for scratch computations
  scratch = (double*)malloc(cfsinp->count*sizeof(double));

  // Find the mean number of cycles required to converge
  // configurations that DID converge
  nsuccess = 0;
  for (i = 0; i < cfsinp->count; i++) {
    if (cvs[i].nstep >= 0) {
      scratch[nsuccess] = cvs[i].nstep;
      nsuccess++;
    }
  }
  avestep = DAverage(scratch, nsuccess);
  stdstep = DStDev(scratch, nsuccess);
  if (cfsinp->rattle == 1) {
    fprintf(outp, "                                   Standard Line   "
            "    RATTLE    \n");
    fprintf(outp, "                                   --------------  "
            "--------------");
  }
  fprintf(outp, "\n - Average steps to convergence:   %5d +/- %4d",
          avestep, stdstep);
  if (cfsinp->rattle == 1) {
    nsuccess = 0;
    for (i = 0; i < cfsinp->count; i++) {
      if (cvs[i].nstep >= 0) {
        scratch[nsuccess] = rattlecvs[i].nstep;
        nsuccess++;
      }
    }
    avestep = DAverage(scratch, nsuccess);
    stdstep = DStDev(scratch, nsuccess);
    fprintf(outp, "  %5d +/- %4d", avestep, stdstep);
  }
  fprintf(outp, "\n - Total converged configurations: %5d", nsuccess);
  fprintf(outp, "\n - Unconverged configurations:     %5d\n",
          cfsinp->count - nsuccess);
  if (cfsinp->reshuffle > 0) {
    fprintf(outp, " - Successful reoptimizations:     %5d\n", cfsinp->nreopt);
    j = 0;
    for (i = 0; i < cfsinp->count; i++) {
      if (cfsinp->rattle == 0 && fabs(cvs[i].dEshuffle) > 1.0e-8) {
        j++;
      }
      else if (cfsinp->rattle == 1 && fabs(rattlecvs[i].dEshuffle) > 1.0e-8) {
        j++;
      }
    }
    fprintf(outp, " - Configs with reoptimizations:   %5d\n", j);
  }
  fprintf(outp, "\n");

  // Describe the residual forces on non-converged configurations
  nfail = 0;
  for (i = 0; i < cfsinp->count; i++) {
    if (cvs[i].nstep < 0) {
      scratch[nfail] = cvs[i].maxF;
      nfail++;
    }
  }
  if (nfail > 0) {
    fprintf(outp, " - Average of the maximal residual force\n"
            "   on NON-converged configurations:      %14.6e "
            "kcal/mol-A\n", DAverage(scratch, nfail));
  }

  // Locate the maximum forces on any non-convergent
  // configurations and generalize the findings.
  if (nfail >= 5) {
    floc = (int*)calloc(tp->natom, sizeof(int));
    for (i = 0; i < cfsinp->count; i++) {
      if (cvs[i].nstep < 0) {
        floc[cvs[i].maxFatom] += 1;
      }
    }
    i = 0;
    probatom = 0;
    fprintf(outp, " - Average residual forces on frequently strained atoms:"
            "\n");
    while (i < 5 && probatom >= 0) {
      probatom = -1;
      maxcount = 0;
      for (j = 0; j < tp->natom; j++) {
        if (floc[j] > maxcount) {
          probatom = j;
          maxcount = floc[j];
        }
      }
      if (probatom >= 0) {
        floc[probatom] = 0;
        k = 0;
        for (j = 0; j < cfsinp->count; j++) {
          if (cvs[j].maxFatom == probatom) {
            scratch[k] = cvs[j].maxF;
            k++;
          }
        }
        fprintf(outp, "     Atom %4d, %4.4s (%2d configs):       %14.6e "
                "kcal/mol-A\n", probatom+1, &tp->AtomNames[4*probatom],
                maxcount, DAverage(scratch, maxcount));
      }
      i++;
    }
    free(floc);
    fprintf(outp, "\n");
  }

  // Print a metric on the bond length constraints
  if (cfsinp->rattle == 1) {
    CheckConstrainedBonds(cfsinp, tp, crdpop, cvs, outp);
  }
  HorizontalRule(outp, 1);

  // Free allocated memory
  free(scratch);
}

//-----------------------------------------------------------------------------
// SanityChecks: print a report on what might have gone wrong in some of the
//               configurations.  What restraints don't work for this molecule
//               and what gets too bent out of shape?
//
// Arguments:
//   cfsinp:    the configuration sampling input data, contains the file
//              format, base name, and extension to write
//   tp:        the topology to work with
//   cvs:       convergence statistics for the run
//   sysUV:     energy components for each configuration
//   outp:      the output file being written
//-----------------------------------------------------------------------------
static void SanityChecks(configs *cfsinp, prmtop *tp, ConvStat* cvs,
                         Energy* sysUV, FILE *outp)
{
  int i, j, atma, atmb, atmc;
  cfigop *myop;
  char* msg;

  // General explanation of this section
  HorizontalRule(outp, 0);
  msg = (char*)malloc(MAXLINE*sizeof(char));
  sprintf(msg, "Basic sanity checks.  If a configuration fails these tests, "
          "it was not printed out for further consideration but the reasons "
          "it failed will be summarized here.  The tests are that, first, no "
          "bond of a configuration may be strained more than %.2f kcal/mol.  "
          "Furthermore, no angle may be strained more than %.2f kcal/mol.  "
          "Finally, the maximum tolerated restraint energy (summed over all "
          "restraints) is %.2f kcal/mol. If any of these tests seem too "
          "stringent or not tight enough, the thresholds may be changed in "
          "the &configs namelist.", cfsinp->maxbstrn, cfsinp->maxastrn,
          cfsinp->strainlim);

  PrintVADesc(0, "(3.)", 4, " ", 1, msg, 74, 0, outp);
  fprintf(outp, "\n - Configurations passing all sanity checks:      %5d\n",
          cfsinp->npass);
  fprintf(outp, " - Configurations failing bond sanity check:      %5d\n",
          cfsinp->nbelim);
  fprintf(outp, " - Configurations failing angle sanity check:     %5d\n",
          cfsinp->naelim);
  fprintf(outp, " - Configurations failing restraint sanity check: %5d\n\n",
          cfsinp->nrelim);
  if (cfsinp->npass < cfsinp->count) {
    PrintVADesc(0, "", 0, " ", 1, "Descriptions of each problematic "
                "configuration follow, along with the names and numbers of "
                "atoms involved in each feature creating the problem.", 74, 0,
                outp);
  }

  // Free allocated memory
  free(msg);

  // Print the details of each failure.
  for (i = 0; i < cfsinp->count; i++) {
    if (cvs[i].pass == 1) {
      continue;
    }
    fprintf(outp, " - Configuration %d:\n", i+1);
    if (cvs[i].maxbstrn > cfsinp->maxbstrn) {
      atma = tp->BLC[cvs[i].bstrnloc[0]].BC[cvs[i].bstrnloc[1]].a;
      atmb = tp->BLC[cvs[i].bstrnloc[0]].BC[cvs[i].bstrnloc[1]].b;
      fprintf(outp, "     Strain is %9.4lf for bond  %4.4s-%4.4s      "
              "(%d, %d)\n", cvs[i].maxbstrn, &tp->AtomNames[4*atma],
              &tp->AtomNames[4*atmb], atma+1, atmb+1);
      j = tp->BLC[cvs[i].bstrnloc[0]].BC[cvs[i].bstrnloc[1]].t;
      fprintf(outp, "     Types %4.4s %4.4s, K = %9.5lf, Leq = %9.5lf, "
              "actual = %9.4lf\n", &tp->AtomTypes[4*atma],
              &tp->AtomTypes[4*atmb], tp->BParam[j].K, tp->BParam[j].l0,
              cvs[i].bstrnloc[2] / 100000.0);
    }
    if (cvs[i].maxastrn > cfsinp->maxastrn) {
      atma = tp->ALC[cvs[i].astrnloc[0]].AC[cvs[i].astrnloc[1]].a;
      atmb = cvs[i].astrnloc[0];
      atmc = tp->ALC[cvs[i].astrnloc[0]].AC[cvs[i].astrnloc[1]].c;
      fprintf(outp, "     Strain is %9.4lf for angle %4.4s-%4.4s-%4.4s "
              "(%d, %d, %d)\n", cvs[i].maxastrn, &tp->AtomNames[4*atma],
              &tp->AtomNames[4*atmb], &tp->AtomNames[4*atmc], atma+1, atmb+1,
              atmc+1);
      j = tp->ALC[cvs[i].astrnloc[0]].AC[cvs[i].astrnloc[1]].t;
      fprintf(outp, "     Types %4.4s %4.4s %4.4s, K = %8.4lf, "
              "Teq = %9.4lf, Actual = %9.4lf\n", &tp->AtomTypes[4*atma],
              &tp->AtomTypes[4*atmb], &tp->AtomTypes[4*atmc], tp->AParam[j].K,
              180.0 / PI * tp->AParam[j].th0,
              180.0 / PI * cvs[i].astrnloc[2] / 1000000.0);
    }
    if (sysUV[i].relec > cfsinp->strainlim) {
      fprintf(outp, "     Total restraint energy %9.4f\n", sysUV[i].relec);
    }
    if (cvs[i].maxrloc >= 0) {
      myop = &cfsinp->ops[cvs[i].maxrloc];
      fprintf(outp, "     Restraint binding atoms ");
      for (j = 0; j < myop->order; j++) {
        fprintf(outp, "%4.4s", &tp->AtomNames[4*myop->anchors[j]]);
        if (j < myop->order-1) {
          fprintf(outp, "-");
        }
      }
      fprintf(outp, "(");
      for (j = 0; j < myop->order; j++) {
        fprintf(outp, "%d", myop->anchors[j]+1);
        if (j < myop->order-1) {
          fprintf(outp, ", ");
        }
      }
      fprintf(outp, ") strained to %9.4f\n", cvs[i].maxrst);
    }
  }
  HorizontalRule(outp, 0);
}

//-----------------------------------------------------------------------------
// PrintOrigins: function to print an additional section to the report giving
//               the origin files of all configurations that were generated.
//               This can be important if reading from a directory, as the
//               directory parsing (dirent.h) does not read files in a way
//               that the average human would see order in.
//
// Arguments:
//   cfsinp:    the configuration sampling input data, contains the file
//              format, base name, and extension to write
//   outp:      the output file being written
//-----------------------------------------------------------------------------
static void PrintOrigins(configs *cfsinp, FILE *outp)
{
  int i;
  char* msg;

  // General explanation of this section
  fprintf(outp, "\n");
  HorizontalRule(outp, 0);
  msg = (char*)malloc(MAXLINE*sizeof(char));
  sprintf(msg, "Origins of each configuration.  The file from which the "
          "initial state was derived is listed next to the configuration "
          "number.");
  PrintVADesc(0, "(4.)", 4, " ", 1, msg, 74, 0, outp);
  fprintf(outp, "\n Config.   Original structure file\n");
  fprintf(outp, " -------   -----------------------\n");
  for (i = 0; i < cfsinp->count; i++) {
    fprintf(outp, "  %6d   %s\n", i+1, cfsinp->origins.map[i]);
  }
  HorizontalRule(outp, 0);
}

//-----------------------------------------------------------------------------
// PrintConfigsReport: print a report on what went on during configuration
//                     sampling.  
//
// Arguments:
//   cfsinp:     the configuration sampling input data, contains the file
//               format, base name, and extension to write
//   tp:         the topology to work with
//   crdpop:     the coordinates of configurations found
//   tj:         trajectory control data (for passing information down to
//               the standard trjaectory writing routines, so they can be
//               co-opted for this purpose)
//   sysUV:      energies of the configurations produced by restrained energy
//               minimization
//   cvs:        convergence statistics (in particular, information on whether
//               each configuration passes muster in the non-constrained run)
//   rattlecvs:  convergence statistics for the follow-on, constrained run
//-----------------------------------------------------------------------------
static void PrintConfigsReport(configs *cfsinp, prmtop *tp, dmat *crdpop,
                               trajcon *tj, Energy* sysUV, Energy* isysUV,
                               ConvStat* cvs, ConvStat* rattlecvs)
{
  int i;
  FILE *outp;
  time_t ct;

  // Input file header
  ct = time(&ct);
  outp = FOpenSafe(tj->outbase, tj->OverwriteOutput);
  PrintSplash(outp);
  fprintf(outp, "Run on %s", asctime(localtime(&ct)));
  HorizontalRule(outp, 0);
  fprintf(outp, "\nINPUT LINE TEXT:\n\n");
  PrintParagraph(tj->inpline, 79, NULL, outp);
  fprintf(outp, "\nINPUT FILE TEXT:\n\n");
  for (i = 0; i < tj->inptext.row; i++) {
    fprintf(outp, "%s", tj->inptext.map[i]);
  }
  HorizontalRule(outp, 1);

  // Print the energy results: what were the initial structure energies (sans
  // restraints), the initial restraint energies, the final structure energies,
  // the final restraint energies?
  ConfigEStats(cfsinp, sysUV, isysUV, cvs, outp);
  ConvergenceStats(cfsinp, tp, crdpop, cvs, rattlecvs, outp);
  SanityChecks(cfsinp, tp, cvs, sysUV, outp);
  if (cfsinp->showorigins == 1) {
    PrintOrigins(cfsinp, outp);
  }

  fclose(outp);
}

//-----------------------------------------------------------------------------
// AssignMobileAtoms: make a mask to assign atom mobility states based on a
//                    series of ambmask strings provided by the user.
//
// Arguments:
//   cfsinp:  configuration sampling input data
//   tp:      the system topology
//-----------------------------------------------------------------------------
static void AssignMobileAtoms(configs *cfsinp, prmtop *tp, dmat *crdpop)
{
  int i, j;
  int* amask;
  coord tc;

  // All atoms are mobile by default
  cfsinp->allmove = 0;
  if (cfsinp->nbelly == 0) {
    cfsinp->movable = (int*)malloc(tp->natom*sizeof(int));
    SetIVec(cfsinp->movable, tp->natom, 1);
    cfsinp->allmove = 1;
    return;
  }

  // If we're still here, assign mobile atoms based on ambmask strings
  // similar to how bellymasks work (only, we're taking as many as the
  // user likes to help reduce the problem of writing complex ambmasks)
  tc = CreateCoord(tp->natom);
  for (i = 0; i < 3*tp->natom; i++) {
    tc.loc[i] = crdpop->map[i][0];
  }
  cfsinp->movable = (int*)calloc(tp->natom, sizeof(int));
  for (i = 0; i < cfsinp->nbelly; i++) {
    amask = ParseAmbMask(cfsinp->belly.map[i], tp, &tc);
    IVec2VecAdd(cfsinp->movable, amask, tp->natom);
    free(amask);
  }
  cfsinp->allmove = (ISum(cfsinp->movable, tp->natom) == tp->natom) ? 1 : 0;
  DestroyCoord(&tc);
}

//-----------------------------------------------------------------------------
// SampleConfigurations: once the job control file has been read, the main
//                       function dives into this routine in order to manage
//                       the configuration sampling.
//
// Arguments:
//   cfsinp:  configuration sampling input data
//   tj:      trajectory control input data (the initial structure will be read
//            from here)
//   tp:      topology for the system (helps to determine whether a given file
//            is compatible as a coordinate set for an initial configuration)
//-----------------------------------------------------------------------------
void SampleConfigurations(configs *cfsinp, trajcon *tj, prmtop *tp)
{
  int i;
  dmat crdpop, crdtrial;
  Energy* sysUV;
  Energy* isysUV;
  Energy* rtlUV;
  Energy* irtlUV;
  Energy* trialUV;
  Energy* itrialUV;
  Energy* rtlTrialUV;
  Energy* irtlTrialUV;
  ConvStat* cvs;
  ConvStat* trialcvs;
  ConvStat* rattlecvs;
  ConvStat* rtlTrialcvs;
  
  // First order of business is to read the ONE topology
  // that the user has supplied.  This was not done as of
  // the time this function was called in main.
  tj->rattle = cfsinp->rattle;
  tp->rattle = cfsinp->rattle;
  GetPrmTop(tp, tj, 1);
  
  // Initialize the original coordinates for each configuration
  // and a corresponding array of energy trackers.
  crdpop = InitCoordPopulation(cfsinp, tj, tp);
  sysUV = (Energy*)malloc(cfsinp->count*sizeof(Energy));
  isysUV = (Energy*)malloc(cfsinp->count*sizeof(Energy));
  itrialUV = (Energy*)malloc(cfsinp->count*sizeof(Energy));
  cvs = (ConvStat*)calloc(cfsinp->count, sizeof(ConvStat));
  if (cfsinp->rattle == 1) {
    rtlUV = (Energy*)malloc(cfsinp->count*sizeof(Energy));
    rtlTrialUV = (Energy*)malloc(cfsinp->count*sizeof(Energy));
    irtlUV = (Energy*)malloc(cfsinp->count*sizeof(Energy));
    irtlTrialUV = (Energy*)malloc(cfsinp->count*sizeof(Energy));
    rattlecvs = (ConvStat*)calloc(cfsinp->count, sizeof(ConvStat));
  }
  if (cfsinp->reshuffle > 0) {
    crdtrial = CreateDmat(3*tp->natom, cfsinp->count, 0);
    trialUV = (Energy*)malloc(cfsinp->count*sizeof(Energy));
  }

  // Get the mask of movable atoms (default is to have all atoms move)
  AssignMobileAtoms(cfsinp, tp, &crdpop);

  // Assign anchor atoms for each operation;
  // decide on target values for each restraint.
  InterpretOperations(cfsinp, tp, &crdpop, tj);

  // Now the configurations can be generated.
  GenerateConfigurations(cfsinp, tp, &crdpop, sysUV, isysUV, cvs, 0);
  if (cfsinp->rattle == 1) {
    GenerateConfigurations(cfsinp, tp, &crdpop, rtlUV, irtlUV, rattlecvs, 1);
  }
  for (i = 0; i < cfsinp->reshuffle; i++) {
    trialcvs = (ConvStat*)calloc(cfsinp->count, sizeof(ConvStat));
    if (cfsinp->rattle == 1) {
      rtlTrialcvs = (ConvStat*)calloc(cfsinp->count, sizeof(ConvStat));
    }
    ShuffleConfigurations(cfsinp, tp, &crdpop, &crdtrial, tj);
    GenerateConfigurations(cfsinp, tp, &crdtrial, trialUV, itrialUV, trialcvs,
                           0);
    if (cfsinp->rattle == 1) {
      GenerateConfigurations(cfsinp, tp, &crdtrial, rtlTrialUV, irtlTrialUV,
                             rtlTrialcvs, 1);
    }
    MergeConfigurations(cfsinp, tp, &crdpop, &crdtrial, sysUV, trialUV, rtlUV,
                     rtlTrialUV, cvs, trialcvs, rattlecvs, rtlTrialcvs);
    free(trialcvs);
    if (cfsinp->rattle == 1) {
      free(rtlTrialcvs);
    }
  }
  if (cfsinp->rattle == 1) {
    ArbitrateConfigurations(cfsinp, sysUV, rattlecvs);
  }
  else {
    ArbitrateConfigurations(cfsinp, sysUV, cvs);
  }
  
  // Verbose output
  if (cfsinp->verbose == 1) {
    printf("\nConfigs >> Configuration optimization complete.\n");
  }

  // Print the configurations in the requested format.
  if (cfsinp->rattle == 0) {
    PrintConfigurations(cfsinp, tp, &crdpop, tj, sysUV, cvs);
  }
  else {
    PrintConfigurations(cfsinp, tp, &crdpop, tj, sysUV, rattlecvs);
  }

  // Print the report on this optimization.  
  PrintConfigsReport(cfsinp, tp, &crdpop, tj, sysUV, isysUV, cvs, rattlecvs);

  // Free what allocated memory we can
  free(sysUV);
  free(isysUV);
  free(itrialUV);
  if (cfsinp->rattle == 1) {
    free(rtlUV);
    free(irtlUV);
    free(rtlTrialUV);
    free(irtlTrialUV);
    free(rattlecvs);
  }
  DestroyDmat(&crdpop);
  DestroyConfigs(cfsinp);

  // We are done.
  exit(1);
}
