#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <dirent.h>
#include "Integrator.h"
#include "Trajectory.h"
#include "mdgxVector.h"
#include "Grid.h"
#include "Matrix.h"
#include "Manual.h"
#include "Parse.h"
#include "Random.h"
#include "ChargeMap.h"
#include "CrdManip.h"
#include "CellManip.h"
#include "pmeRecip.h"
#include "BSpline.h"
#include "Timings.h"
#include "Macros.h"
#include "Constants.h"
#include "BSpline.h"
#include "IPolQ.h"
#include "ChargeFit.h"
#include "Topology.h"

#ifdef MPI
#include "BroadcastCommand.h"
#include "MPIMap.h"
#endif

//-----------------------------------------------------------------------------
// EvenSphere: function for arranging n points at equal distances about 
//             the origin, on the unit spehere, using a 1/r repulsion  
//             potential.                                               
//                                                                      
// Arguments:                                                           
//   n:         the number of points to arrange                         
//   verbose:   flag to alert the user as to progress                   
//-----------------------------------------------------------------------------
double* EvenSphere(int n, int verbose)
{
  int i, j, npts, pts_remaining, pts_placed, improvement, nmove;
  double band_sep, radius, elevation, theta, pt_angle, ix, iy, iz, d, fmag;
  double f_totx, f_toty, f_totz, dx, dy, dz, fx, fy, fz, max_step, max_force;
  double new_potential, curr_potential;
  double* tmp;
  double* tm2p;
  double* fitmp;
  double* fjtmp;
  double* projf;
  double* evsph;
  double* new_coords;
  double* forces;

  // Allocate memory for the points
  evsph = (double*)malloc(3*n*sizeof(double));
  new_coords = (double*)malloc(3*n*sizeof(double));
  forces = (double*)malloc(3*n*sizeof(double));
  projf = (double*)malloc(3*sizeof(double));
  if (verbose == 1) {
    printf("\n");
  }

  // Place points in a semi-intelligent fashion
  evsph[0] = 0.0;
  evsph[1] = 0.0;
  evsph[2] = -1.0;
  band_sep = sqrt(4.0*PI/n);
  theta = PI-band_sep;
  radius = sin(theta);
  elevation = cos(theta);
  pts_remaining = n-1;
  pts_placed = 1;
  while(PI*radius/band_sep > 1.0 && theta > 0.0 && pts_remaining > 0) {

    // Place points around this band
    npts = MIN(2.0*PI*radius/band_sep + 1, pts_remaining);
    pt_angle = 2.0*PI/npts;
    for (i = pts_placed; i < pts_placed + npts; i++) {
      tmp = &evsph[3*i];
      tmp[2] = elevation;
      tmp[0] = radius*cos(i*pt_angle);
      tmp[1] = radius*sin(i*pt_angle);
    }
    pts_remaining -= npts;
    pts_placed += npts;
    theta -= band_sep;
    radius = sin(theta);
    elevation = cos(theta);
  }
  if (pts_remaining > 0) {
    if (theta <= 0.0) {
      theta += band_sep;
      theta *= 0.75;
    }
    if (PI*radius/band_sep <= 1.0) {
      theta += band_sep/2.0;
    }
    radius = sin(theta);
    elevation = cos(theta);
    pt_angle = 2.0*PI/pts_remaining;
    for (i = pts_placed; i < n; i++) {
      tmp = &evsph[3*i];
      tmp[2] = elevation;
      tmp[0] = radius*cos(i*pt_angle);
      tmp[1] = radius*sin(i*pt_angle);
    }
  }

  // Even out the distribution
  max_step = band_sep/10.0;
  curr_potential = -1.0;
  while(max_step > band_sep*0.0001 && n > 1) {
    for (i = 0; i < 3*n; i++) {
      forces[i] = 0.0;
    }

    // Compute forces
    for (i = 0; i < n-1; i++) {
      fitmp = &forces[3*i];
      tmp = &evsph[3*i];
      ix = tmp[0];
      iy = tmp[1];
      iz = tmp[2];
      f_totx = 0.0;
      f_toty = 0.0;
      f_totz = 0.0;
      for (j = i+1; j < n; j++) {
        tm2p = &evsph[3*j];
        dx = tm2p[0] - ix;
        dy = tm2p[1] - iy;
        dz = tm2p[2] - iz;
        d = dx*dx + dy*dy + dz*dz;
        if (d < 0.25) {
          fmag = 1.0/(d*sqrt(d));
          fx = fmag*dx;
          fy = fmag*dy;
          fz = fmag*dz;
          fjtmp = &forces[3*j];
          fjtmp[0] += fx;
          fjtmp[1] += fy;
          fjtmp[2] += fz;
          f_totx -= fx;
          f_toty -= fy;
          f_totz -= fz;
        }
      }
      fitmp[0] += f_totx;
      fitmp[1] += f_toty;
      fitmp[2] += f_totz;
    }

    // Isolate tangetial forces
    max_force = 0.0;
    for (i = 0; i < n; i++) {
      tmp = &evsph[3*i];
      fitmp = &forces[3*i];
      Project(fitmp, tmp, projf, 3);
      fitmp[0] -= projf[0];
      fitmp[1] -= projf[1];
      fitmp[2] -= projf[2];
      for (j = 0; j < 3; j++) {
        fitmp[j] -= projf[j];
      }
      if (fitmp[0]*fitmp[0] + fitmp[1]*fitmp[1] + fitmp[2]*fitmp[2] > 
          max_force) {
        max_force = fitmp[0]*fitmp[0] + fitmp[1]*fitmp[1] + fitmp[2]*fitmp[2];
      }
    }

    // Scale forces
    max_force = 1.0/sqrt(max_force);
    if (verbose == 1) {
      fprintf(stderr, "\rEvenSphere >> Max Force = %16.6lf  ||  Max Step = "
              "%16.6lf", 1.0/max_force, max_step);
      fflush(stderr);
    }
    for (i = 0; i < n; i++) {
      forces[3*i] *= max_force;
      forces[3*i+1] *= max_force;
      forces[3*i+2] *= max_force;
    }

    // Compute potential of current state
    if (curr_potential < 0.0) {
      curr_potential = 0.0;
      for (i = 0; i < n-1; i++) {
        tmp = &evsph[3*i];
        ix = tmp[0];
        iy = tmp[1];
        iz = tmp[2];
        for (j = i+1; j < n; j++) {
          tm2p = &evsph[3*j];
          dx = tm2p[0]-ix;
          dy = tm2p[1]-iy;
          dz = tm2p[2]-iz;
          d = dx*dx + dy*dy + dz*dz;
          if (d < 0.25) {
            curr_potential += 1.0/sqrt(d);
          }
        }
      }
    }

    // Improve the result
    improvement = 1;
    nmove = 0;
    while (improvement == 1) {

      // Compute a move along the gradient
      for (i = 0; i < n; i++) {
        tmp = &evsph[3*i];
        tm2p = &new_coords[3*i];
        fitmp = &forces[3*i];
        tm2p[0] = tmp[0] + fitmp[0]*max_step;
        tm2p[1] = tmp[1] + fitmp[1]*max_step;
        tm2p[2] = tmp[2] + fitmp[2]*max_step;
        d = 1.0/sqrt(tm2p[0]*tm2p[0] + tm2p[1]*tm2p[1] + tm2p[2]*tm2p[2]);
        tm2p[0] *= d;
        tm2p[1] *= d;
        tm2p[2] *= d;
      }

      // Compute the new potential
      new_potential = 0.0;
      for (i = 0; i < n-1; i++) {
        ix = new_coords[3*i];
        iy = new_coords[3*i+1];
        iz = new_coords[3*i+2];
        for (j = i+1; j < n; j++) {
          dx = new_coords[3*j] - ix;
          dy = new_coords[3*j+1] - iy;
          dz = new_coords[3*j+2] - iz;
          d = dx*dx + dy*dy + dz*dz;
          if (d < 0.25) {
            new_potential += 1.0/sqrt(d);
          }
        }
      }
      if (new_potential < curr_potential) {
        SWAP(evsph, new_coords, tmp);
        curr_potential = new_potential;
        nmove++;
      }
      else {
        improvement = 0;
      }
    }
    if (nmove == 0) {
      max_step *= 0.8;
    }
    else if (nmove < 5) {
      max_step *= 1.0 + 0.01*nmove;
    }
    else {
      max_step *= 1.05;
    }
  }
  if (verbose == 1) {
    printf("\n");
  }

  // Free allocated memory
  free(projf);
  free(forces);
  free(new_coords);

  return evsph;
}

//-----------------------------------------------------------------------------
// DefineBoundary: generates a set of equidistant points on a spherical surface
//                 and then uses that set of points to define a surface of
//                 points spread around a solute molecule at an arbitrary
//                 distance from at least one atom.        
//
// Arguments:
//   crd:     the coordinates of the molecule
//   nshell:  the number of shells to use (for N shells there will be N
//            discontinuous surfaces
//   r:       a vector of real numbers giving the radii of all shells
//   nbq:     the maximum number of points to place on each shell, per atom
//   tj:      trajectory control data (used here for the random number
//            generator handle)
//-----------------------------------------------------------------------------
static coord DefineBoundary(coord *crd, int nshell, double* r, int nbq,
                            trajcon *tj)
{
  int h, i, j, k, npt, isvalid, j3, j3p1, j3p2, nbqL;
  int* surflj;
  double dx, dy, dz, atmx, atmy, atmz, r2, r2min, rfac, ljt;
  double* sph;
  double* surf;
  dmat R;
  coord ts;

  R = CreateDmat(3, 3, 0);
  surf = (double*)malloc(3*nbq*nshell*crd->natom*sizeof(double));
  surflj = (int*)malloc(nbq*nshell*crd->natom*sizeof(int));
  npt = 0;
  rfac = sqrt(8.0*PI/nbq);
  for (h = 0; h < nshell; h++) {

    // A single surface point at the atom, or
    // multiple points on a sphere surface   
    if (r[h] < 0.0) {
      continue;
    }
    else if (r[h] < 1.0e-8) {
      sph = (double*)calloc(3, sizeof(double));
      nbqL = 1;
      r2min = -1.0;
      ljt = 0;
    }
    else {
      sph = EvenSphere(nbq, 0);
      nbqL = nbq;
      r2min = r[h]*r[h];
      ljt = 1;
    }
    for (i = 0; i < 3*nbqL; i++) {
      sph[i] *= r[h];
    }
    for (i = 0; i < crd->natom ; i++) {
      atmx = crd->loc[3*i];
      atmy = crd->loc[3*i+1];
      atmz = crd->loc[3*i+2];
      dx = rfac*(0.5-ran2(&tj->rndcon));
      dy = rfac*(0.5-ran2(&tj->rndcon));
      dz = rfac*(0.5-ran2(&tj->rndcon));
      BeardRotMat(dx, dy, dz, &R);
      for (j = 0; j < nbqL; j++) {
        j3 = 3*j;
        j3p1 = j3+1;
        j3p2 = j3+2;
        dx = R.data[0]*sph[j3] + R.data[1]*sph[j3p1] + R.data[2]*sph[j3p2];
        dy = R.data[3]*sph[j3] + R.data[4]*sph[j3p1] + R.data[5]*sph[j3p2];
        dz = R.data[6]*sph[j3] + R.data[7]*sph[j3p1] + R.data[8]*sph[j3p2];
        surf[3*npt] = dx + atmx;
        surf[3*npt+1] = dy + atmy;
        surf[3*npt+2] = dz + atmz;
        surflj[npt] = ljt;

        // Check to see if this point is valid
        isvalid = 1;
        for (k = 0; k < crd->natom; k++) {
          if (k == i) {
            continue;
          }
          dx = surf[3*npt] - crd->loc[3*k];
          dy = surf[3*npt+1] - crd->loc[3*k+1];
          dz = surf[3*npt+2] - crd->loc[3*k+2];
          r2 = dx*dx + dy*dy + dz*dz;
          if (r2 < r2min) {
            isvalid = 0;
            break;
          }
        }
        npt += isvalid;
      }
    }

    // Free allocated memory
    free(sph);
  }

  // Create a new coord struct with the surface points
  ts = CreateCoord(npt);
  ReflectDVec(ts.loc, surf, 3*npt);
  ReflectIVec(ts.atmid, surflj, npt);
  free(surf);
  free(surflj);
  CopyDmat(&ts.U, &crd->U, 1);
  CopyDmat(&ts.invU, &crd->invU, 1);
  ReflectDVec(ts.gdim, crd->gdim, 6);
  ts.isortho = crd->isortho;

  return ts;
}

//-----------------------------------------------------------------------------
// PrepIPolQ: allocate data structures and detect Orca executables in   
//            preparation for IPolQ calculations.                       
//                                                                      
// Arguments:                                                           
//   tj:       trajectory control data (for the step count, primarily)  
//   ipqinp:   IPolQ control data                                       
//   tp:       system topology                                          
//   crd:      system coordinates                                       
//   qswap:    the charge array to swap during IPolQ calculations       
//-----------------------------------------------------------------------------
static void PrepIPolQ(trajcon *tj, ipqcon *ipqinp, prmtop *tp, coord *crd,
                      double* qswap)
{
  int i, j;
  int* atommask;
  double qval, totq;
  coord solute;

  // Run through all charge swaps
  for (i = 0; i < ipqinp->nqmod; i++) {
    atommask = ParseAmbMask(ipqinp->QModMask.map[i], tp, crd);
    qval = ipqinp->QModVal[i];
    for (j = 0; j < tp->natom; j++) {
      if (atommask[j] == 1) {
        qswap[j] = qval;
      }
    }
  }

  // Check that a bellymask has been made
  if (tj->Leash.active == 0 || tj->Leash.usebelly == 0) {
    printf("PrepIPolQ >> A frozen solute molecule must be specified.\n");
    exit(1);
  }

  // Set the step count to zero and determine the minimum number
  // of steps needed to satisfy the sampling requirements.      
  if (ipqinp->neqstep % ipqinp->ntqs > 0) {
    ipqinp->neqstep += ipqinp->ntqs - ipqinp->neqstep % ipqinp->ntqs;
  }
  tj->nstep = ipqinp->neqstep + ipqinp->ntqs * ipqinp->nqframe;
  UpdateStepNumber(tj, 0);

  // Define the solute
  j = 0;
  totq = 0.0;
  for (i = 0; i < tp->natom; i++) {
    if (tp->MobileAtoms[i] == 0) {
      totq += tp->Charges[i];
      if (tp->Masses[i] > 1.0e-8) {
        j++;
      }
    }
  }
  ipqinp->stotq = lround(totq);
  solute = CreateCoord(j);
  j = 0;
  for (i = 0; i < tp->natom; i++) {
    if (tp->MobileAtoms[i] == 0 && tp->Masses[i] > 1.0e-8) {
      solute.atmid[j] = i;
      solute.loc[3*j] = crd->loc[3*i];
      solute.loc[3*j+1] = crd->loc[3*i+1];
      solute.loc[3*j+2] = crd->loc[3*i+2];
      j++;
    }
  }
  solute.natom = j;
  solute.isortho = crd->isortho;
  CopyDmat(&solute.U, &crd->U, 1);
  CopyDmat(&solute.invU, &crd->invU, 1);
  ReflectDVec(solute.gdim, crd->gdim, 6);
  ipqinp->Solute = solute;

  // Create the set of points that will describe the
  // electrostatic potential at and around the atom 
  // sites, then create the set of points at which  
  // to place charges to aid in the reproduction of 
  // that potential.  In parallel mode, compute the 
  // surfaces on the master and broadcast to slaves,
  // to keep the surfaces consistent as generating  
  // them does rely on random numbers.              
  if (tj->tid == 0) {
    ipqinp->Vsurf = DefineBoundary(&solute, ipqinp->nVshell, ipqinp->Vshell,
                                   ipqinp->nVphpt, tj);
    ipqinp->Qsurf = DefineBoundary(&solute, ipqinp->nQshell, ipqinp->Qshell,
                                   ipqinp->nQphpt, tj);
  }
#ifdef MPI
  int nVQ[2];
  if (tj->tid == 0) {
    nVQ[0] = ipqinp->Vsurf.natom;
    nVQ[1] = ipqinp->Qsurf.natom;
  }
  MPI_Bcast(nVQ, 2, MPI_INT, 0, tj->SysComm[0]);
  if (tj->tid != 0) {
    ipqinp->Vsurf = CreateCoord(nVQ[0]);
    ipqinp->Qsurf = CreateCoord(nVQ[1]);
  }
  BroadcastCoordinates(&ipqinp->Vsurf, tj, 0);
  BroadcastCoordinates(&ipqinp->Qsurf, tj, 0);
#endif

  // Prepare for block average accumulation of forces at points
  // of Vsurf.  nqframe is set to zero to become a counter.    
  ipqinp->Vfrc = (double*)calloc(3*ipqinp->Vsurf.natom, sizeof(double));
  ipqinp->SAfrc = CreateDmat(ipqinp->nqframe, 3*solute.natom, 0);
  ipqinp->nqframe = 0;

  // Prepare to accumulate the explicit charge density
  ipqinp->qnbrs = (int*)malloc(tp->natom*sizeof(int));
  ipqinp->nQcloud = 0;
  ipqinp->Qcloud = CreateDmat(tp->natom, 4, 0);
}

//-----------------------------------------------------------------------------
// Surf2Cells: this routine creates a fake "cell grid" to feed into the 
//             CellGridIntrp function so that test charges at the       
//             required surface points may sample electrostatic         
//             potentials and gradients.                                
//
// Arguments:
//   CG:     the cell grid of the original system (we will make a new cell grid
//           with as many cells in the same configuration as the grid for the
//           standard molecular system)
//   surf:   the surface of points to map into cells
//   rcinp:  reciprocal space control data
//-----------------------------------------------------------------------------
static cell* Surf2Cells(cellgrid *CG, coord *surf, reccon *rcinp)
{
  int i, j;
  int nc[3], ng[3];
  int* atomcount;
  int* gridbin;
  double bcrd[3];
  double* gridpos;
  cell* CGP;
  atomc *catm;

  // Count the surface points in each cell
  atomcount = (int*)calloc(CG->ncell, sizeof(int));
  gridpos = (double*)malloc(3*surf->natom*sizeof(double));
  gridbin = (int*)malloc(surf->natom*sizeof(int));
  ng[0] = CG->ng[0];
  ng[1] = CG->ng[1];
  ng[2] = CG->ng[2];
  for (i = 0; i < surf->natom; i++) {
    DMatVecMult(&surf->U, &surf->loc[3*i], bcrd);
    for (j = 0; j < 3; j++) {
      bcrd[j] -= floor(bcrd[j]);
      nc[j] = bcrd[j] * CG->ng[j];
    }
    DMatVecMult(&surf->invU, bcrd, &gridpos[3*i]);
    gridbin[i] = (nc[0]*ng[1] + nc[1])*ng[2] + nc[2];
    atomcount[gridbin[i]] += 1;
  }

  // Allocate cells
  CGP = (cell*)malloc(CG->ncell*sizeof(cell));
  for (i = 0; i < CG->ncell; i++) {
    CGP[i] = CreateCell(atomcount[i], rcinp->ordr);
    CGP[i].nr[0] = 0;
  }

  // Place atoms in cells
  for (i = 0; i < surf->natom; i++) {
    catm = &CGP[gridbin[i]].data[CGP[gridbin[i]].nr[0]];
    catm->id = i;
    catm->lj = surf->atmid[i];
    catm->q = 1.0;
    catm->loc[0] = gridpos[3*i];
    catm->loc[1] = gridpos[3*i+1];
    catm->loc[2] = gridpos[3*i+2];
    CGP[gridbin[i]].nr[0] += 1;
  }

  // Free allocated memory
  free(atomcount);
  free(gridpos);
  free(gridbin);

  return CGP;
}

//-----------------------------------------------------------------------------
// AllEleNeighbors: compute direct-space interactions between surface points
//                  and the actual atoms of the system, looping over all 27
//                  nearest neighbor cells.                 
//                                                                      
// Arguments:                                                           
//   CG:     the cell grid, containing locations of system atoms        
//   CGP:    array of cells which parallels the cell grid, contains locations
//           of surface points and buffers for accumulating forces on them
//   n:      the number of the cell (as indexed in this process's list  
//           of cells from CG) to start in                              
//   surf:   surface point coord struct (for transformation matrices)   
//   Etab:   tabulated coarse electrostatic potential                   
//   EHtab:  tabulated fine electrostatic potential                     
//-----------------------------------------------------------------------------
static void AllEleNeighbors(cellgrid *CG, cell* CGP, int n, coord *surf,
                            FrcTab *Etab, FrcTab *EHtab, dircon *dcinp,
                            ipqcon *ipqinp)
{
  int i, j, k, ai, aj, xCp, yCp, zCp, xbound, ybound, zbound, isortho;
  int ir2, reim, ilim, jlim, klim;
  int ivec[3], jvec[3], kvec[3];
  int* qnbrs;
  double cspc, hspc, ptx, pty, ptz, dx, dy, dz, r2min, r2, aiq, fmag;
  double shr2max;
  cell *C, *Cp, *CX;
  CSpln *CFspl, *HFspl;

  // Unpack
  C = &CG->data[n];
  Cp = &CGP[C->gbin[3]];
  isortho = surf->isortho;
  cspc = Etab->ivdr;
  hspc = EHtab->ivdr;
  CFspl = Etab->dSD;
  HFspl = EHtab->dSD;
  r2min = dcinp->Mcut*dcinp->Mcut;
  qnbrs = ipqinp->qnbrs;
  shr2max = ipqinp->Qshell[0] * ipqinp->Qshell[0];

  // First question: is this cell at a grid boundary?
  xCp = 0;
  yCp = 0;
  zCp = 0;
  if (C->gbin[0] == 0) xCp = -1;
  if (C->gbin[0] == CG->ng[0]-1) xCp = 1;
  if (C->gbin[1] == 0) yCp = -1;
  if (C->gbin[1] == CG->ng[1]-1) yCp = 1;
  if (C->gbin[2] == 0) zCp = -1;
  if (C->gbin[2] == CG->ng[2]-1) zCp = 1;

  // Set the cells to loop over
  ivec[0] = (C->gbin[0] > 0) ? C->gbin[0]-1 : CG->ng[0]-1;
  jvec[0] = (C->gbin[1] > 0) ? C->gbin[1]-1 : CG->ng[1]-1;
  kvec[0] = (C->gbin[2] > 0) ? C->gbin[2]-1 : CG->ng[2]-1;
  ivec[1] = C->gbin[0];
  jvec[1] = C->gbin[1];
  kvec[1] = C->gbin[2];
  if (CG->ng[0] > 2) {
    ivec[2] = (C->gbin[0] < CG->ng[0]-1) ? C->gbin[0]+1 : 0;
    ilim = 3;
  }
  else {
    ilim = 2;
  }
  if (CG->ng[1] > 2) {
    jvec[2] = (C->gbin[1] < CG->ng[1]-1) ? C->gbin[1]+1 : 0;
    jlim = 3;
  }
  else {
    jlim = 2;
  }
  if (CG->ng[2] > 2) {
    kvec[2] = (C->gbin[2] < CG->ng[2]-1) ? C->gbin[2]+1 : 0;
    klim = 3;
  }
  else {
    klim = 2;
  }
  for (i = 0; i < ilim; i++) {
    xbound = (ivec[i] == 0) ? -xCp : (ivec[i] == CG->ng[0]-1) ? xCp : 0;
    for (j = 0; j < jlim; j++) {
      ybound = (jvec[j] == 0) ? -yCp : (jvec[j] == CG->ng[1]-1) ? yCp : 0;
      for (k = 0; k < klim; k++) {
        zbound = (kvec[k] == 0) ? -zCp : (kvec[k] == CG->ng[2]-1) ? zCp : 0;

        // Set pointer to the other cell, but continue
        // if it is not owned by this process.        
        CX = &CG->map[ivec[i]][jvec[j]][kvec[k]];
        if (CX->CGRank != CG->tid) {
          continue;
        }
        reim = (xbound == -1 || ybound == -1 || zbound == -1) ? 1 : 0;
        for (ai = 0; ai < CX->nr[0]; ai++) {
          ptx = CX->data[ai].loc[0];
          pty = CX->data[ai].loc[1];
          ptz = CX->data[ai].loc[2];
          aiq = CX->data[ai].q;
          if (fabs(aiq) < 1.0e-8) {
            continue;
          }
          for (aj = 0; aj < Cp->nr[0]; aj++) {

            // Squared distance calculation
            dx = Cp->data[aj].loc[0] - ptx;
            dy = Cp->data[aj].loc[1] - pty;
            dz = Cp->data[aj].loc[2] - ptz;
            if (reim == 1) {
              if (isortho) {
                OrthoReim(&dx, &dy, &dz, &surf->U, &surf->invU);
              }
              else {
                NonOrthoReim(&dx, &dy, &dz, &surf->U, &surf->invU);
              }
            }
            r2 = dx*dx + dy*dy + dz*dz;

            // Force calculation and accumulation
            if (r2 < r2min && r2 >= MINNB2) {
              ir2 = r2*cspc;
              fmag = (((CFspl[ir2].A*r2 + CFspl[ir2].B)*r2 + 
                       CFspl[ir2].C)*r2 + CFspl[ir2].D)*aiq;
            }
            else if (r2 < MINNB2) {
              ir2 = r2*hspc;
              fmag = (((HFspl[ir2].A*r2 + HFspl[ir2].B)*r2 + 
                       HFspl[ir2].C)*r2 + HFspl[ir2].D)*aiq;
            }
            else {
              fmag = 0.0;
            }
            Cp->data[aj].frc[0] -= fmag*dx;
            Cp->data[aj].frc[1] -= fmag*dy;
            Cp->data[aj].frc[2] -= fmag*dz;

            // Atom ranging for charge accumulation
            if (Cp->data[aj].lj == 0 && r2 < shr2max) {
              qnbrs[CX->data[ai].id] = 1;
            }
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// UpdateChargePool: update the charge pool with charges and locations from the
//                   the latest snapshot.                          
//                                                                      
// Arguments:                                                           
//   ipqinp:     the IPolQ input parameter structure, where all the data is
//               also being stored                                   
//   tp:         system topology, for referencing atomic partial charges
//   CG:         the cell grid, where we get atomic coordinates         
//   qswap:      alternative charge array (used for SRFP calculations)  
//-----------------------------------------------------------------------------
static void UpdateChargePool(ipqcon *ipqinp, prmtop *tp, cellgrid *CG,
                             double* qswap)
{
  int i, j, nqc, maxqc;
  int* qnbtmp;
  double *dtmp;
  dmat *qcloud;
  cell *C;

  qnbtmp = ipqinp->qnbrs;
  qcloud = &ipqinp->Qcloud;
  nqc = ipqinp->nQcloud;
  maxqc = ipqinp->Qcloud.row;
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < C->nr[0]; j++) {
      if (qnbtmp[C->data[j].id] == 1) {
        dtmp = qcloud->map[nqc];
        dtmp[0] = C->data[j].loc[0];
        dtmp[1] = C->data[j].loc[1];
        dtmp[2] = C->data[j].loc[2];
        dtmp[3] = qswap[C->data[j].id];
        nqc++;
        if (nqc == maxqc) {
          maxqc += tp->natom;
          ipqinp->Qcloud = ReallocDmat(&ipqinp->Qcloud, maxqc, 4);
          qcloud = &ipqinp->Qcloud;
        }
      }
    }
  }
  ipqinp->nQcloud = nqc;
}

//-----------------------------------------------------------------------------
// AccSoluteForces: accumulate forces acting on solute atoms into a matrix.    
//                  The results can be analyzed later in terms of convergence.
//
// Arguments:
//   ipqinp:      the IPolQ input control data
//-----------------------------------------------------------------------------
static void AccSoluteForces(ipqcon *ipqinp)
{
  int i, j, nqframe;
  int *ljtmp;
  dmat *SAmat;

  // Stash the SRFP gradient on all solute atoms
  ljtmp = ipqinp->Vsurf.atmid;
  SAmat = &ipqinp->SAfrc;
  nqframe = ipqinp->nqframe;
  j = 0;
  for (i = 0; i < ipqinp->Vsurf.natom; i++) {
    if (ljtmp[i] == 0) {
      SAmat->map[nqframe][3*j] = ipqinp->Vsurf.frc[3*i];
      SAmat->map[nqframe][3*j+1] = ipqinp->Vsurf.frc[3*i+1];
      SAmat->map[nqframe][3*j+2] = ipqinp->Vsurf.frc[3*i+2];
      j++;
    }
  }

  // Update the sample counter
  ipqinp->nqframe += 1;
  if (ipqinp->nqframe == SAmat->row) {
    ipqinp->SAfrc = ReallocDmat(&ipqinp->SAfrc, ipqinp->nqframe+50,
                                ipqinp->SAfrc.col);
  }
}

//-----------------------------------------------------------------------------
// ComputeSRFP: compute the solvent reaction-field potential, defining the
//              "solvent" as any mobile atoms in the system and the "solute"
//              as any immobile atoms.                         
//                                                                      
// Arguments:                                                           
//   See PerformIPolQ function in this library                          
//-----------------------------------------------------------------------------
static void ComputeSRFP(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
                        FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
                        execon *etimers, ipqcon *ipqinp, double* qswap)
{
  int i, j, k;
  cell *C;
  cell* CGP;
#ifdef MPI
  int nreq;
  MPI_Request* req;
  MPI_Status* stt;

  // Allocate memory for MPI requests
  nreq = CG->nthreads;
  req = (MPI_Request*)malloc(nreq*sizeof(MPI_Request));
  stt = (MPI_Status*)malloc(nreq*sizeof(MPI_Status));
#endif

  // Swap charges with the SRFP generator array and
  // set partial charges of solute atoms to zero   
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < C->nr[0]; j++) {
      if (tp->MobileAtoms[C->data[j].id] == 0) {
        C->data[j].q = 0.0;
      }
      else {
        C->data[j].q = qswap[C->data[j].id];
      }
    }
  }

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
  nreq = InitMeshGatherSPME(rcinp, CG, req, 0);
  etimers->mpiMeshPack += mdgxStopTimer(etimers);
  if (nreq > 0) {
    MPI_Waitall(nreq, req, stt);
  }
  MPI_Barrier(CG->dspcomm);
  etimers->mpiMeshPullWait += mdgxStopTimer(etimers);
#endif
  if (CG->tid == 0) {
#ifdef MPI
    // Unpack the mesh components from other processes
    for (j = 1; j < CG->nthreads; j++) {
      RecvMeshPart(rcinp, CG, j);
    }
    etimers->mpiMeshPack += mdgxStopTimer(etimers);
#endif
    ConvQBC(rcinp, crd, &rcinp->QL[0], PPk, etimers);
  }
#ifdef MPI
  InitMeshGatherSPME(rcinp, CG, req, 1);
  etimers->mpiMeshPack += mdgxStopTimer(etimers);
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

  // At this point, we have an electrostatic potential grid
  // computed in the absence of solute partial charges.  We
  // can sample the grid at will to pull down electrostatic
  // potentials and forces on test charges.  We can compute
  // the direct space interactions using the existing cell 
  // decomposition and the primary sector contents.        
  CGP = Surf2Cells(CG, &ipqinp->Vsurf, rcinp);
  for (i = 0; i < CG->ncell; i++) {
    C = &CGP[i];
    for (j = 0; j < C->nr[0]; j++) {
      C->data[j].frc[0] = 0.0;
      C->data[j].frc[1] = 0.0;
      C->data[j].frc[2] = 0.0;
    }
  }
  SetIVec(ipqinp->qnbrs, tp->natom, 0);
  for (i = 0; i < CG->ncell; i++) {
    AllEleNeighbors(CG, CGP, i, &ipqinp->Vsurf, Etab, EHtab, dcinp, ipqinp);
  }
  UpdateChargePool(ipqinp, tp, CG, qswap);

  // It is also necessary to compute the effects of the 
  // reciprocal space electrostatics on the test charges
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CGP[CG->data[CG->MyCellDomain[i]].gbin[3]];
    CellSplCoeff(C, &ipqinp->Vsurf, rcinp);
    CellIntrpFrc(C, &ipqinp->Vsurf, rcinp, &rcinp->QL[0]);
  }

  // Upload forces to the list
  for (i = 0; i < CG->ncell; i++) {
    C = &CGP[i];
    for (j = 0; j < C->nr[0]; j++) {
      k = C->data[j].id;
      ipqinp->Vsurf.prvfrc[3*k] = C->data[j].frc[0];
      ipqinp->Vsurf.prvfrc[3*k+1] = C->data[j].frc[1];
      ipqinp->Vsurf.prvfrc[3*k+2] = C->data[j].frc[2];
    }
  }

#ifdef MPI
  // Reduce the forces to track statistics at each step.
  if (CG->tid == 0) {
    SetDVec(ipqinp->Vsurf.frc, 3*ipqinp->Vsurf.natom, 0.0);
  }
  MPI_Reduce(ipqinp->Vsurf.prvfrc, ipqinp->Vsurf.frc, 3*ipqinp->Vsurf.natom,
             MPI_DOUBLE, MPI_SUM, 0, CG->dspcomm);
#else
  double *dtmp;
  SWAP(ipqinp->Vsurf.prvfrc, ipqinp->Vsurf.frc, dtmp);
#endif

  // Accumulate the forces into block averages
  AccSoluteForces(ipqinp);
  DVec2VecAdd(ipqinp->Vfrc, ipqinp->Vsurf.frc, 3*ipqinp->Vsurf.natom);

  // Reset partial charges of solute atoms
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < C->nr[0]; j++) {
      C->data[j].q = tp->Charges[C->data[j].id];
    }
  }

  // Free allocated memory
  for (i = 0; i < CG->ncell; i++) {
    DestroyCell(&CGP[i]);
  }
  free(CGP);
#ifdef MPI
  free(req);
  free(stt);
#endif
}

#ifdef MPI
//-----------------------------------------------------------------------------
// PoolCharges: pool all charges from various processes involved in the MD.
//              The master process accumulates a complete charge cloud to
//              manipulate for quantum calculations.           
//                                                                      
// Arguments:                                                           
//   ipqinp:   IPolQ input data, containing tables of points and charges
//   CG:       cell grid, for communicator and thread tracking          
//-----------------------------------------------------------------------------
static void PoolCharges(ipqcon *ipqinp, cellgrid *CG)
{
  int i, j, llim, hlim, totalQc;
  int* nq;
  int* npq;
  double* Qc;

  // Allocate memory for MPI requests
  if (CG->nthreads == 1) {
    return;
  }

  // Reduce the charge counts
  nq = (int*)calloc(CG->nthreads, sizeof(int));
  npq = (int*)calloc(CG->nthreads, sizeof(int));
  nq[CG->tid] = ipqinp->nQcloud;
  MPI_Allreduce(nq, npq, CG->nthreads, MPI_INT, MPI_MAX, CG->dspcomm);

  // Pool the charge data
  totalQc = ISum(npq, CG->nthreads);
  ipqinp->Qcloud = ReallocDmat(&ipqinp->Qcloud, totalQc, 4);
  Qc = (double*)calloc(4*totalQc, sizeof(double));
  llim = 4*ISum(npq, CG->tid);
  hlim = llim + 4*ipqinp->nQcloud;
  j = 0;
  for (i = llim; i < hlim; i++) {
    Qc[i] = ipqinp->Qcloud.data[j];
    j++;
  }
  SetDVec(ipqinp->Qcloud.data, 4*totalQc, 0.0);
  MPI_Reduce(Qc, ipqinp->Qcloud.data, 4*totalQc, MPI_DOUBLE, MPI_SUM, 0,
             CG->dspcomm);
  if (CG->tid == 0) {
    ipqinp->nQcloud = totalQc;
  }

  // Free allocated memory
  free(nq);
  free(npq);
  free(Qc);
}
#endif

//-----------------------------------------------------------------------------
// SnapQVPoints: re-image all points in a list to a new origin.         
//                                                                      
// Arguments:                                                           
//   L:          the list of points                                     
//   n:          the number of atoms                                    
//   stride:     L[i*stride] accesses the ith point                     
//   c[x,y,z]:   location of the new origin                             
//   crd:        coordinates (for box information)                      
//-----------------------------------------------------------------------------
static void SnapQVPoints(double* L, int n, int stride, double cx, double cy,
                         double cz, coord *crd)
{
  int i;
  double *dtmp;

  for (i = 0; i < n; i++) {
    dtmp = &L[i*stride];
    dtmp[0] -= cx;
    dtmp[1] -= cy;
    dtmp[2] -= cz;
    if (crd->isortho == 1) {
      OrthoReim(&dtmp[0], &dtmp[1], &dtmp[2], &crd->U, &crd->invU);
    }
    else {
      NonOrthoReim(&dtmp[0], &dtmp[1], &dtmp[2], &crd->U, &crd->invU);
    }
    dtmp[0] += cx;
    dtmp[1] += cy;
    dtmp[2] += cz;
  }
}

//-----------------------------------------------------------------------------
// CenterQCloud: re-arrange the atoms of the charge cloud to place them all in
//               the same image.  This is accomplished by finding the center
//               of the solute coordinates and re-imaging all charge cloud
//               points to that.                           
//                                                                      
// Arguments:                                                           
//   ipqinp:   IPolQ input data, containing tables of points and charges
//   crd:      system coordinates (maintained by the dynamics routines) 
//-----------------------------------------------------------------------------
static void CenterQCloud(ipqcon *ipqinp, coord *crd)
{
  double scn[3];
  coord *solute;

  // Find the solute center of coordinates
  solute = &ipqinp->Solute;
  FindCoordCenter(solute->loc, solute->prvloc, 0, solute->natom, scn);

  // Re-image all solute atoms, charge cloud points, electrostatic
  // field samples, and all shell charge locations                
  SnapQVPoints(solute->loc, solute->natom, 3, scn[0], scn[1], scn[2], crd);
  SnapQVPoints(ipqinp->Qcloud.data, ipqinp->nQcloud, 4, scn[0], scn[1], scn[2],
               crd);
  SnapQVPoints(ipqinp->Vsurf.loc, ipqinp->Vsurf.natom, 3, scn[0], scn[1],
               scn[2], crd);
  SnapQVPoints(ipqinp->Qsurf.loc, ipqinp->Qsurf.natom, 3, scn[0], scn[1],
               scn[2], crd);
}

//-----------------------------------------------------------------------------
// AnalSRFPConv: analyze the convergence of the solvent reaction field  
//               potential (SRFP) at solute atom sites.                 
//                                                                      
// Arguments:                                                           
//   ipqinp:   IPolQ input data, containing tables of points and charges
//-----------------------------------------------------------------------------
static void AnalSRFPConv(ipqcon *ipqinp, FILE *outp)
{
  int h, i, j, blkcount;
  double *tVfrc;
  dmat blockave, tave;

  blkcount = ipqinp->nqframe / ipqinp->nblock;
  h = 0;
  blockave = CreateDmat(ipqinp->nblock, ipqinp->SAfrc.col, 0);
  for (i = 0; i < ipqinp->nblock; i++) {
    for (j = 0; j < blkcount; j++) {
      DVec2VecAdd(blockave.map[i], ipqinp->SAfrc.map[h], blockave.col);
      h++; 
    }
    DVecMult(blockave.map[i], blockave.col, 1.0/blkcount);
  }
  tave = CreateDmat(blockave.col, blockave.row, 0);
  for (i = 0; i < tave.row; i++) {
    for (j = 0; j < tave.col; j++) {
      tave.map[i][j] = blockave.map[j][i];
    }
  }
  fprintf(outp, "Convergence of the SRF at solute atom sites:\n");
  fprintf(outp, 
          "     Sim. X    Sim. Y    Sim. Z     Conv. X   Conv. Y   Conv. Z\n"
          "   --------- --------- ---------   --------- --------- ---------"
          "\n");
  tVfrc = ipqinp->Vfrc;
  for (i = 0; i < ipqinp->Solute.natom; i++) {
    fprintf(outp, "   %9.4lf %9.4lf %9.4lf   %9.4lf %9.4lf %9.4lf\n",
            tVfrc[3*i], tVfrc[3*i+1], tVfrc[3*i+2],
            DStDev(tave.map[3*i], tave.col), DStDev(tave.map[3*i+1], tave.col),
            DStDev(tave.map[3*i+2], tave.col));
  }

  // Free allocated memory
  DestroyDmat(&blockave);
  DestroyDmat(&tave);
}

//-----------------------------------------------------------------------------
// TrimQCloud: trim the charge cloud to eliminate charges that fall too 
//             close to the molecule.  This has to be done here because 
//             all of the solute atom positions must be known in order  
//             to compute it.                                           
//                                                                      
// Arguments:                                                           
//   ipqinp:   the IPolQ input control data
//-----------------------------------------------------------------------------
static void TrimQCloud(ipqcon *ipqinp)
{
  int i, j;
  int* expel;
  double atmx, atmy, atmz, dx, dy, dz, shr2min;
  double *dtmp;
  dmat *Qc;
  coord *S;

  // Bail out if there is no inner shell defined
  if (ipqinp->Qshell[3] <= 0.0) {
    return;
  }

  // Prune all points within the innermost shell
  Qc = &ipqinp->Qcloud;
  if (ipqinp->verbose == 1) {
    printf("\nIPolQ >> Trimming charge cloud from %d atoms to ", Qc->row);
  }
  dtmp = Qc->data;
  S = &ipqinp->Solute;
  expel = (int*)calloc(Qc->row, sizeof(int));
  shr2min = ipqinp->Qshell[3] * ipqinp->Qshell[3];
  for (i = 0; i < S->natom; i++) {
    atmx = S->loc[3*i];
    atmy = S->loc[3*i+1];
    atmz = S->loc[3*i+2];
    for (j = 0; j < Qc->row; j++) {
      dx = dtmp[4*j] - atmx;
      dy = dtmp[4*j+1] - atmy;
      dz = dtmp[4*j+2] - atmz;
      if (dx*dx + dy*dy + dz*dz < shr2min) {
        expel[j] = 1;
      }
    }
  }
  j = 0;
  for (i = 0; i < Qc->row; i++) {
    if (expel[i] == 0) {
      dtmp[4*j] = dtmp[4*i];
      dtmp[4*j+1] = dtmp[4*i+1];
      dtmp[4*j+2] = dtmp[4*i+2];
      dtmp[4*j+3] = dtmp[4*i+3];
      j++;
    }
  }
  free(expel);
  ipqinp->Qcloud = ReallocDmat(&ipqinp->Qcloud, j, 4);
  ipqinp->nQcloud = j;
  if (ipqinp->verbose == 1) {
    printf("%d.\n", j);
  }
}

//-----------------------------------------------------------------------------
// AtomCode: produce the element name based on its atomic number or mass in
//           the AMBER topology.
//                                                                      
// Arguments:                                                           
//   mass:     the mass of the atom                                     
//   znum:     the Z-number (returned, and inferred if given as zero)
//   charge:   atomic partial charge
//   tlet:     pre-allocated string to hold the result                  
//-----------------------------------------------------------------------------
double AtomCode(double mass, int znum, double charge, char* tlet)
{
  int i;

  // Go by atomic number if available
  if (znum == 1) {
    sprintf(tlet, "H");
  }
  if (znum == 2) {
    sprintf(tlet, "He");
  }
  if (znum == 3) {
    sprintf(tlet, "Li");
  }
  if (znum == 4) {
    sprintf(tlet, "Be");
  }
  if (znum == 5) {
    sprintf(tlet, "B");
  }
  if (znum == 6) {
    sprintf(tlet, "C");
  }
  if (znum == 7) {
    sprintf(tlet, "N");
  }
  if (znum == 8) {
    sprintf(tlet, "O");
  }
  if (znum == 9) {
    sprintf(tlet, "F");
  }
  if (znum == 10) {
    sprintf(tlet, "Ne");
  }
  if (znum == 11) {
    sprintf(tlet, "Na");
  }
  if (znum == 12) {
    sprintf(tlet, "Mg");
  }
  if (znum == 13) {
    sprintf(tlet, "Al");
  }
  if (znum == 14) {
    sprintf(tlet, "Si");
  }
  if (znum == 15) {
    sprintf(tlet, "P");
  }
  if (znum == 16) {
    sprintf(tlet, "S");
  }
  if (znum == 17) {
    sprintf(tlet, "Cl");
  }
  if (znum == 18) {
    sprintf(tlet, "Ar");
  }
  if (znum == 19) {
    sprintf(tlet, "K");
  }
  if (znum == 20) {
    sprintf(tlet, "Ca");
  }
  if (znum == 21) {
    sprintf(tlet, "Sc");
  }
  if (znum == 22) {
    sprintf(tlet, "Ti");
  }
  if (znum == 23) {
    sprintf(tlet, "V");
  }
  if (znum == 24) {
    sprintf(tlet, "Cr");
  }
  if (znum == 25) {
    sprintf(tlet, "Mn");
  }
  if (znum == 26) {
    sprintf(tlet, "Fe");
  }
  if (znum == 27) {
    sprintf(tlet, "Co");
  }
  if (znum == 28) {
    sprintf(tlet, "Ni");
  }
  if (znum == 29) {
    sprintf(tlet, "Cu");
  }
  if (znum == 30) {
    sprintf(tlet, "Zn");
  }
  if (znum == 31) {
    sprintf(tlet, "Ga");
  }
  if (znum == 32) {
    sprintf(tlet, "Ge");
  }
  if (znum == 33) {
    sprintf(tlet, "As");
  }
  if (znum == 34) {
    sprintf(tlet, "Se");
  }
  if (znum == 35) {
    sprintf(tlet, "Br");
  }
  if (znum == 36) {
    sprintf(tlet, "Kr");
  }
  if (znum == 47) {
    sprintf(tlet, "Ag");
  }
  if (znum == 53) {
    sprintf(tlet, "I");
  }
  if (znum == 54) {
    sprintf(tlet, "Xe");
  }
  if (znum == 55) {
    sprintf(tlet, "Cs");
  }
  if (znum > 0) {
    return (double)znum;
  }
  
  // If unsuccessful, go by mass and charge
  if (mass < 0.1) {
    sprintf(tlet, "EP");
    return 0.0;
  }
  if (fabs(mass - 1.008) < 0.1) {
    sprintf(tlet, "H");
    return 1.0;
  }
  if (fabs(mass - 4.003) < 0.1) {
    sprintf(tlet, "He");
    return 2.0;
  }
  if (fabs(mass - 6.940) < 0.1) {
    sprintf(tlet, "Li");
    return 3.0;
  }
  if (fabs(mass - 9.012) < 0.1) {
    sprintf(tlet, "Be");
    return 4.0;
  }
  if (fabs(mass - 10.810) < 0.1) {
    sprintf(tlet, "B");
    return 5.0;
  }
  if (fabs(mass - 12.011) < 0.1) {
    sprintf(tlet, "C");
    return 6.0;
  }
  if (fabs(mass - 14.007) < 0.1) {
    sprintf(tlet, "N");
    return 7.0;
  }
  if (fabs(mass - 15.999) < 0.1) {
    sprintf(tlet, "O");
    return 8.0;
  }
  if (fabs(mass - 18.998) < 0.1) {
    sprintf(tlet, "F");
    return 9.0;
  }
  if (fabs(mass - 20.180) < 0.1) {
    sprintf(tlet, "Ne");
    return 10.0;
  }
  if (fabs(mass - 22.990) < 0.1) {
    sprintf(tlet, "Na");
    return 11.0;
  }
  if (fabs(mass - 24.305) < 0.1) {
    sprintf(tlet, "Mg");
    return 12.0;
  }
  if (fabs(mass - 26.982) < 0.1) {
    sprintf(tlet, "Al");
    return 13.0;
  }
  if (fabs(mass - 28.085) < 0.1) {
    sprintf(tlet, "Si");
    return 14.0;
  }
  if (fabs(mass - 30.974) < 0.1) {
    sprintf(tlet, "P");
    return 15.0;
  }
  if (fabs(mass - 32.060) < 0.1) {
    sprintf(tlet, "S");
    return 16.0;
  }
  if (fabs(mass - 35.450) < 0.1) {
    sprintf(tlet, "Cl");
    return 17.0;
  }
  if (fabs(mass - 39.948) < 0.1 && fabs(charge) < 0.1) {
    sprintf(tlet, "Ar");
    return 18.0;
  }
  if (fabs(mass - 39.098) < 0.1) {
    sprintf(tlet, "K");
    return 19.0;
  }
  if (fabs(mass - 40.078) < 0.1 && fabs(charge - 2.0) < 1.0) {
    sprintf(tlet, "Ca");
    return 20.0;
  }
  if (fabs(mass - 44.956) < 0.1) {
    sprintf(tlet, "Sc");
    return 21.0;
  }
  if (fabs(mass - 47.867) < 0.1) {
    sprintf(tlet, "Ti");
    return 22.0;
  }
  if (fabs(mass - 50.942) < 0.1) {
    sprintf(tlet, "V");
    return 23.0;
  }
  if (fabs(mass - 51.996) < 0.1) {
    sprintf(tlet, "Cr");
    return 24.0;
  }
  if (fabs(mass - 54.938) < 0.1) {
    sprintf(tlet, "Mn");
    return 25.0;
  }
  if (fabs(mass - 55.845) < 0.1) {
    sprintf(tlet, "Fe");
    return 26.0;
  }
  if (fabs(mass - 58.933) < 0.1) {
    sprintf(tlet, "Co");
    return 27.0;
  }
  if (fabs(mass - 58.693) < 0.1) {
    sprintf(tlet, "Ni");
    return 28.0;
  }
  if (fabs(mass - 63.546) < 0.1) {
    sprintf(tlet, "Cu");
    return 29.0;
  }
  if (fabs(mass - 65.380) < 0.1) {
    sprintf(tlet, "Zn");
    return 30.0;
  }
  if (fabs(mass - 69.723) < 0.1) {
    sprintf(tlet, "Ga");
    return 31.0;
  }
  if (fabs(mass - 72.630) < 0.1) {
    sprintf(tlet, "Ge");
    return 32.0;
  }
  if (fabs(mass - 74.922) < 0.1) {
    sprintf(tlet, "As");
    return 33.0;
  }
  if (fabs(mass - 78.960) < 0.1) {
    sprintf(tlet, "Se");
    return 34.0;
  }
  if (fabs(mass - 79.904) < 0.1) {
    sprintf(tlet, "Br");
    return 35.0;
  }
  if (fabs(mass - 83.798) < 0.1) {
    sprintf(tlet, "Kr");
    return 36.0;
  }
  if (fabs(mass - 107.870) < 0.1) {
    sprintf(tlet, "Ag");
    return 47.0;
  }
  if (fabs(mass - 126.90) < 0.1) {
    sprintf(tlet, "I");
    return 53.0;
  }
  if (fabs(mass - 131.29) < 0.1) {
    sprintf(tlet, "Xe");
    return 54.0;
  }
  if (fabs(mass - 132.91) < 0.1) {
    sprintf(tlet, "Cs");
    return 55.0;
  }

  // Die because no element could be matched.
  printf("AtomCode >> Error.  Atom mass %8.4lf and charge %8.4lf corresponds "
         "to no\nAtomCode >> listed element.\n", mass, charge);
  exit(1);

  return 0.0;
}

//-----------------------------------------------------------------------------
// FitQShell: this function subtracts the contributions of explicit chargs to
//            the solvent reaction field potential at sites selected on and
//            near solute atoms, then computes values of shell charges.
//                                                                      
// Arguments:                                                           
//   ipqinp:   IPolQ input data, containing tables of points and charges
//   tj:       trajectory control data (for output file name)           
//   tp:       system topology (for atomic numbers)                     
//-----------------------------------------------------------------------------
static void FitQShell(ipqcon *ipqinp, trajcon *tj, prmtop *tp)
{
  int i, j;
  int *smtmp;
  double dx, dy, dz, qx, qy, qz, vx, vy, vz, dr2, rmin;
  double invdr, invdr2, fmag, invnf, qfac, Utot, Unuc, tcharge;
  double *Atmpx, *Atmpy, *Atmpz, *qltmp, *vltmp, *tVfrc, *vfptmp;
  double *sltmp;
  double* b;
  char tlet[8];
  char* fname;
  FILE *outp;
  coord *Vs, *Qs;
  dmat A;

  // Scale down explicit charges
  invnf = 1.0/ipqinp->nqframe;
  tcharge = 0.0;
  for (i = 0; i < ipqinp->nQcloud; i++) {
    ipqinp->Qcloud.map[i][3] *= invnf;
    tcharge += ipqinp->Qcloud.map[i][3];
  }
  tcharge += ipqinp->stotq;

  // Compute explicit charges' contributions.  Use the    
  // prvfrc array in the Vsurf struct to store the result.
  Vs = &ipqinp->Vsurf;
  tVfrc = ipqinp->Vfrc;
  for (i = 0; i < 3*Vs->natom; i++) {
    tVfrc[i] *= invnf;
  }
  vltmp = Vs->loc;
  vfptmp = Vs->prvfrc;
  SetDVec(vfptmp, 3*Vs->natom, 0.0);
  rmin = 9.0e2;
  for (i = 0; i < ipqinp->nQcloud; i++) {
    qx = ipqinp->Qcloud.map[i][0];
    qy = ipqinp->Qcloud.map[i][1];
    qz = ipqinp->Qcloud.map[i][2];
    qfac = BIOQ*ipqinp->Qcloud.map[i][3];
    for (j = 0; j < Vs->natom; j++) {
      dx = vltmp[3*j] - qx;
      dy = vltmp[3*j+1] - qy;
      dz = vltmp[3*j+2] - qz;
      dr2 = dx*dx + dy*dy + dz*dz;
      invdr2 = 1.0/dr2;
      invdr = sqrt(invdr2);
      if (dr2 < rmin) {
        rmin = dr2;
      }
      fmag = qfac*invdr2*invdr;
      vfptmp[3*j] += fmag*dx;
      vfptmp[3*j+1] += fmag*dy;
      vfptmp[3*j+2] += fmag*dz;
    }
  }
  rmin = sqrt(rmin);

  // Fill out the matrix for shell charges
  Qs = &ipqinp->Qsurf;
  qltmp = Qs->loc;
  A = CreateDmat(3*Vs->natom + Qs->natom + 1, Qs->natom, 0);
  b = (double*)malloc((3*Vs->natom + Qs->natom + 1)*sizeof(double));
  for (i = 0; i < Vs->natom; i++) {
    vx = vltmp[3*i];
    vy = vltmp[3*i+1];
    vz = vltmp[3*i+2];
    Atmpx = A.map[3*i];
    Atmpy = A.map[3*i+1];
    Atmpz = A.map[3*i+2];
    for (j = 0; j < Qs->natom; j++) {
      dx = qltmp[3*j] - vx;
      dy = qltmp[3*j+1] - vy;
      dz = qltmp[3*j+2] - vz;
      invdr2 = 1.0/(dx*dx + dy*dy + dz*dz);
      invdr = sqrt(invdr2);
      fmag = BIOQ*invdr2*invdr;
      Atmpx[j] = -fmag*dx;
      Atmpy[j] = -fmag*dy;
      Atmpz[j] = -fmag*dz;
    }
    b[3*i] = tVfrc[3*i] - vfptmp[3*i];
    b[3*i+1] = tVfrc[3*i+1] - vfptmp[3*i+1];
    b[3*i+2] = tVfrc[3*i+2] - vfptmp[3*i+2];
  }
  for (i = 0; i < Qs->natom; i++) {
    A.map[3*Vs->natom+i][i] = ipqinp->minqfac;
    b[3*Vs->natom+i] = 0.0;
    A.map[3*Vs->natom+Qs->natom][i] = 10000.0;
  }
  b[3*Vs->natom+Qs->natom] = -10000.0*tcharge;

  // Solve the matrix equation
  if (ipqinp->verbose == 1) {
    printf("IPolQ >> Fitting shell charges.\n");
  }
  AxbQRRxc(A, b, ipqinp->verbose);
  BackSub(A, b);
  if (ipqinp->verbose == 1) {
    printf("IPolQ >> Shell charge fitting complete.\n");
  }

  // Extend the charge cloud
  ipqinp->Qcloud = ReallocDmat(&ipqinp->Qcloud, ipqinp->nQcloud + Qs->natom,
                               4);
  j = ipqinp->nQcloud;
  for (i = 0; i < Qs->natom; i++) {
    ipqinp->Qcloud.map[j][0] = Qs->loc[3*i];
    ipqinp->Qcloud.map[j][1] = Qs->loc[3*i+1];
    ipqinp->Qcloud.map[j][2] = Qs->loc[3*i+2];
    ipqinp->Qcloud.map[j][3] = b[i];
    j++;
  }
  ipqinp->nQcloud += Qs->natom;

  // Analyze the result
  SetDVec(vfptmp, 3*Vs->natom, 0.0);
  for (i = 0; i < ipqinp->nQcloud; i++) {
    qx = ipqinp->Qcloud.map[i][0];
    qy = ipqinp->Qcloud.map[i][1];
    qz = ipqinp->Qcloud.map[i][2];
    qfac = BIOQ*ipqinp->Qcloud.map[i][3];
    for (j = 0; j < Vs->natom; j++) {
      dx = vltmp[3*j] - qx;
      dy = vltmp[3*j+1] - qy;
      dz = vltmp[3*j+2] - qz;
      invdr2 = 1.0/(dx*dx + dy*dy + dz*dz);
      invdr = sqrt(invdr2);
      fmag = qfac*invdr2*invdr;
      vfptmp[3*j] += fmag*dx;
      vfptmp[3*j+1] += fmag*dy;
      vfptmp[3*j+2] += fmag*dz;
    }
  }

  // Compute charge:charge interaction energy in
  // the cloud and between nuclei and the cloud.
  Utot = 0.0;
  Unuc = 0.0;
  qltmp = ipqinp->Qcloud.data;
  sltmp = ipqinp->Solute.loc;
  smtmp = ipqinp->Solute.atmid;
  for (i = 0; i < ipqinp->nQcloud-1; i++) {
    qx = qltmp[4*i];
    qy = qltmp[4*i+1];
    qz = qltmp[4*i+2];
    qfac = 332.0636*qltmp[4*i+3];
    for (j = i+1; j < ipqinp->nQcloud; j++) {
      dx = qltmp[4*j] - qx;
      dy = qltmp[4*j+1] - qy;
      dz = qltmp[4*j+2] - qz;
      Utot += qfac*qltmp[4*j+3]/sqrt(dx*dx + dy*dy + dz*dz);
    }
    for (j = 0; j < ipqinp->Solute.natom; j++) {
      dx = sltmp[3*j] - qx;
      dy = sltmp[3*j+1] - qy;
      dz = sltmp[3*j+2] - qz;
      Unuc += qfac*AtomCode(tp->Masses[smtmp[j]], tp->ZNumber[smtmp[j]],
                            tp->Charges[smtmp[j]], tlet) /
        sqrt(dx*dx + dy*dy + dz*dz);
    }
  }

  // Write to output
  fname = (char*)malloc(MAXNAME*sizeof(char));
  SpliceFileName(tj, tj->outbase, tj->outsuff, fname, 1);
  outp = fopen(fname, "a");
  HorizontalRule(outp, 0);
  fprintf(outp, " (7.) Solvent Reaction Field (Potential) (SRF(P))\n\n");
  fprintf(outp, " SRF evaluated at %d points, averaged over %d frames\n",
          Vs->natom, ipqinp->nqframe);
  fprintf(outp, " SRFP was approximated by: %7d scaled simulation charges\n"
          "                           %7d shell charges\n\n",
          ipqinp->nQcloud-Qs->natom, Qs->natom);
  fprintf(outp, " Closest explicit charge penetration to SRFP evaluation "
          "point: %8.4lf A\n\n", rmin);
  fprintf(outp, " Total charge-charge interaction energy in the cloud: "
          "%12.4lf kcal/mol\n", Utot);
  fprintf(outp, " Total nuclear interaction energy with the cloud:     "
          "%12.4lf kcal/mol\n", Unuc);
  fprintf(outp, " SRF approximated to %8.4lf kcal/(mol-A-e) RMS error "
          "(correlation %7.4lf)\n\n", VecRMSD(tVfrc, vfptmp, 3*Vs->natom),
          Pearson(tVfrc, vfptmp, 3*Vs->natom));
  fprintf(outp, " SRF at solute atom sites:\n");
  fprintf(outp, 
          "     Sim. X    Sim. Y    Sim. Z     Appr. X   Appr. Y   Appr. Z\n"
          "   --------- --------- ---------   --------- --------- ---------"
          "\n");
  for (i = 0; i < Vs->natom; i++) {
    if (Vs->atmid[i] == 0) {
      fprintf(outp, "   %9.4lf %9.4lf %9.4lf   %9.4lf %9.4lf %9.4lf\n",
              tVfrc[3*i], tVfrc[3*i+1], tVfrc[3*i+2], vfptmp[3*i],
              vfptmp[3*i+1], vfptmp[3*i+2]);
    }
  }
  fprintf(outp, "\n");
  AnalSRFPConv(ipqinp, outp);
  fclose(outp);

  // Free allocated memory
  DestroyDmat(&A);
  free(b);
  free(fname);
}

//-----------------------------------------------------------------------------
// XformChargeCloud: transform the charge cloud to follow the solute    
//                   coordinates.                                       
//                                                                      
// Arguments:                                                           
//   ipqinp:   IPolQ input data, containing tables of points and charges
//             as well as orca run parameters                           
//-----------------------------------------------------------------------------
static void XformChargeCloud(ipqcon *ipqinp, double* tI, double* tII, dmat U)
{
  int i, nQ;
  double* qcrd;
  dmat *Qc;

  nQ = ipqinp->nQcloud;
  Qc = &ipqinp->Qcloud;
  qcrd = (double*)malloc(3*nQ*sizeof(double));
  for (i = 0; i < nQ; i++) {
    qcrd[3*i] = Qc->map[i][0];
    qcrd[3*i+1] = Qc->map[i][1];
    qcrd[3*i+2] = Qc->map[i][2];
  }
  TransCrd(qcrd, nQ, tI, -1.0);
  RotateCrd(qcrd, nQ, U);
  TransCrd(qcrd, nQ, tII, 1.0);
  for (i = 0; i < nQ; i++) {
    Qc->map[i][0] = qcrd[3*i];
    Qc->map[i][1] = qcrd[3*i+1];
    Qc->map[i][2] = qcrd[3*i+2];
  }
  free(qcrd);
}

//-----------------------------------------------------------------------------
// WriteGaussianInput: write input for Gaussian calculations with IPolQ.
//                                                                      
// Arguments:                                                           
//   ipqinp:   IPolQ input data, containing tables of points and charges
//             as well as orca run parameters                           
//   tp:       the molecular topology (for atom masses and mobilities)  
//   env:      the environment (0 = vacuum, 1 = IPolQ solvent)          
//   tj:       trajectory control, for output overwriting permissions   
//-----------------------------------------------------------------------------
static void WriteGaussianInput(ipqcon *ipqinp, prmtop *tp, int env,
                               trajcon *tj)
{
  int i, j, nqmt;
  int *gdim;
  double origx, origy, origz;
  double Scofm[3], Sstdcofm[3];
  double *gspc;
  char inpname[MAXNAME], tlet[8], line[128], envcomm[64];
  FILE *outp, *vacinp;
  coord Sstd;
  coord *S;
  dmat U;
  dmat *Qc;

  // Open the input file for writing, with the proper extension
  if (env == 0) {
    sprintf(inpname, "%s.vacu", ipqinp->inpfile);
    sprintf(envcomm, ",");
  }
  else {
    sprintf(inpname, "%s.solv", ipqinp->inpfile);
    sprintf(envcomm, ", charge,");
  }
  outp = FOpenSafe(inpname, tj->OverwriteOutput);

  // Write the input file
  if (env == 0) {
    nqmt = tj->nthreads;
    if (ipqinp->nQMThreads > 0) {
      nqmt = ipqinp->nQMThreads;
    }
    fprintf(outp, "--Link1--\n%%nproc=%d\n%%mem=%dMB\n%%chk=%s.vacu.chk\n#"
            "%s/%s density=%s%s maxdisk=16GB, pop=(minimal, MK), 5d, "
            "punch=mo\n\n", nqmt, ipqinp->MaxCore, ipqinp->outfile,
            ipqinp->qmmeth, ipqinp->basis, ipqinp->qmmeth, envcomm);
  }
  else if (env == 1) {
    fprintf(outp, "--Link1--\n%%nproc=%d\n%%mem=%dMB\n%%chk=%s.solv.chk\n#"
            "%s/%s density=%s%s maxdisk=16GB, pop=(minimal, MK), 5d, "
            "punch=mo\n\n", tj->nthreads, ipqinp->MaxCore, ipqinp->outfile,
            ipqinp->qmmeth, ipqinp->basis, ipqinp->qmmeth, envcomm);
  }
  fprintf(outp, "Generated by mdgx\n\n%d   1\n", ipqinp->stotq);
  S = &ipqinp->Solute;
  if (env == 1) {

    // If we need a solvent charge environment, we must 
    // re-orient the charge field and solute coordinates
    // to match a standard orientation.                 
    Sstd = CreateCoord(S->natom);
    sprintf(inpname, "%s.vacu", ipqinp->outfile);
    if ((vacinp = fopen(inpname, "r")) == NULL) {
      printf("WriteGaussianInput >> Error.  Vacuum run output %s not found.\n",
             inpname);
      exit(1);
    }
    while (fgets(line, 128, vacinp) != NULL) {
      RemoveWhiteSpace(line, 0, ' ');
      if (line[0] == 'S' && strncmp(line, "Standard orientation:", 21) == 0) {
        for (i = 0; i < 4; i++) {
          fgets(line, 128, vacinp);
        }
        for (i = 0; i < Sstd.natom; i++) {
          fgets(line, 128, vacinp);
          sscanf(&line[31], "%lf%lf%lf", &Sstd.loc[3*i], &Sstd.loc[3*i+1],
                 &Sstd.loc[3*i+2]);
        }
      }
    }
    fclose(vacinp);
    FindCoordCenter(S->loc, S->prvloc, 0, S->natom, Scofm);
    FindCoordCenter(Sstd.loc, Sstd.prvloc, 0, Sstd.natom, Sstdcofm);
    TransCrd(S->loc, S->natom, Scofm, -1.0);
    TransCrd(Sstd.loc, Sstd.natom, Sstdcofm, -1.0);
    U = CreateDmat(3, 3, 0);
    QuatAlign(Sstd.loc, S->loc, S->natom, S->prvloc, 0, &U);
    RotateCrd(S->loc, S->natom, U);
    TransCrd(S->loc, S->natom, Sstdcofm, 1.0);

    // Manipulate the charge cloud appropriately
    XformChargeCloud(ipqinp, Scofm, Sstdcofm, U);
  }

  // Now print the solute and perhaps charge cloud coordinates
  for (i = 0; i < S->natom; i++) {

    // The coord struct's atmid array stores indices into the master topology
    j = S->atmid[i];
    AtomCode(tp->Masses[j], tp->ZNumber[j], tp->Charges[j], tlet);
    fprintf(outp, "    %s  %12.7lf %12.7lf %12.7lf\n", tlet,
            S->loc[3*i], S->loc[3*i+1], S->loc[3*i+2]);
  }
  fprintf(outp, "\n");
  if (env == 1) {
    Qc = &ipqinp->Qcloud;
    for (i = 0; i < ipqinp->nQcloud; i++) {
      fprintf(outp, "%11.6lf %11.6lf %11.6lf  %12.9lf\n", Qc->map[i][0],
              Qc->map[i][1], Qc->map[i][2], Qc->map[i][3]);
    }
  }
  fprintf(outp, "\n\n");
  fclose(outp);

  // Print electrostatic potential evaluation input
  if (env == 0) {
    if (ipqinp->CenterGrid == 1) {
      FindCoordCenter(S->loc, S->prvloc, 0, S->natom, Scofm);
    }
    else {
      Scofm[0] = 0.0;
      Scofm[1] = 0.0;
      Scofm[2] = 0.0;
    }
    gdim = ipqinp->gdim;
    gspc = ipqinp->gspc;
    origx = (gdim[0] % 2 == 0) ?
      Scofm[0] - (gdim[0]/2 - 0.5)*gspc[0] :
      Scofm[0] - (gdim[0]/2)*gspc[0];
    origy = (gdim[1] % 2 == 0) ?
      Scofm[1] - (gdim[1]/2 - 0.5)*gspc[1] :
      Scofm[1] - (gdim[1]/2)*gspc[1];
    origz = (gdim[2] % 2 == 0) ?
      Scofm[2] - (gdim[2]/2 - 0.5)*gspc[2] :
      Scofm[2] - (gdim[2]/2)*gspc[2];
    sprintf(inpname, "%s.eval", ipqinp->inpfile);
    outp = FOpenSafe(inpname, tj->OverwriteOutput);
    fprintf(outp, "-11, %9.4lf, %9.4lf, %9.4lf\n", origx, origy, origz);
    fprintf(outp, "%4d, %9.4lf, %9.4lf, %9.4lf\n", gdim[0], gspc[0], 0.0, 0.0);
    fprintf(outp, "%4d, %9.4lf, %9.4lf, %9.4lf\n", gdim[1], 0.0, gspc[1], 0.0);
    fprintf(outp, "%4d, %9.4lf, %9.4lf, %9.4lf\n", gdim[2], 0.0, 0.0, gspc[2]);
    fclose(outp);
  }

  // Free allocated memory
  if (env == 1) {
    DestroyCoord(&Sstd);
    DestroyDmat(&U);
  }
}

//-----------------------------------------------------------------------------
// WriteOrcaInput: write input for an orca calculation with IPolQ.      
//                                                                      
// Arguments:                                                           
//   Same as WriteGaussianInput above                                   
//-----------------------------------------------------------------------------
static void WriteOrcaInput(ipqcon *ipqinp, prmtop *tp, int env, trajcon *tj)
{
  int h, i, j, k, imin, imax;
  int *gdim;
  double dx, dy, dz, origx, origy, origz;
  double scofm[3];
  double *gspc;
  char inpname[MAXNAME], tlet[8], engrad[16];
  FILE *outp;
  coord *S;
  dmat *Qc;

  // Open the input file for writing, with the proper extension
  if (env == 0) {
    sprintf(inpname, "%s.vacu", ipqinp->inpfile);
  }
  else {
    sprintf(inpname, "%s.solv", ipqinp->inpfile);
  }
  outp = FOpenSafe(inpname, tj->OverwriteOutput);
  
  // Write the input file
  if (ipqinp->excitation == 1 || ipqinp->DoSinglePoint == 1) {
    sprintf(engrad, "EnGrad");
  }
  else {
    engrad[0] = '\0';
  }
  fprintf(outp, "! %s %s TightSCF %s KeepDens NoKeepInts\n\n",
          ipqinp->qmmeth, ipqinp->basis, engrad);
  if (ipqinp->excitation == 1) {
    fprintf(outp, "%%tddft\n  KeepInts false\n  NRoots %d\n  IRoot %d\n"
            " end\n\n", ipqinp->excitation+1, ipqinp->excitation);
  }
  if (tj->nthreads > 1 || ipqinp->nQMThreads > 0) {
    if (ipqinp->nQMThreads > 0) {
      fprintf(outp, "%%pal\n  nprocs %d\n end\n", ipqinp->nQMThreads);
    }
    else {
      fprintf(outp, "%%pal\n  nprocs %d\n end\n", tj->nthreads);
    }
  }
  if (strcmp(ipqinp->qmmeth, "MP2") == 0) {
    fprintf(outp, "%%scf\n  KeepInts false\n end\n");
    fprintf(outp, "%%mp2\n  Density relaxed\n  MaxCore %d\n end\n\n",
            ipqinp->MaxCore);
  }
  if (env == 1) {
    fprintf(outp, "%%pointcharges \"%s\"\n\n", ipqinp->ptqfile);
  }
  fprintf(outp, "*xyz  %d  1\n", ipqinp->stotq);
  S = &ipqinp->Solute;
  for (i = 0; i < S->natom; i++) {

    // The coord struct's atmid array stores indices into the master topology
    j = S->atmid[i];
    AtomCode(tp->Masses[j], tp->ZNumber[j], tp->Charges[j], tlet);
    fprintf(outp, "    %s  %12.7lf %12.7lf %12.7lf\n", tlet,
            S->loc[3*i], S->loc[3*i+1], S->loc[3*i+2]);
  }
  fprintf(outp, "*\n");
  fclose(outp);
  
  // Write input for the electrostatic potential evaluation
  if (env == 0) {
    if (ipqinp->DoSinglePoint == 0) {
      if (ipqinp->CenterGrid == 1) {
        FindCoordCenter(S->loc, S->prvloc, 0, S->natom, scofm);
      }
      else {
        scofm[0] = 0.0;
        scofm[1] = 0.0;
        scofm[2] = 0.0;
      }
      gdim = ipqinp->gdim;
      gspc = ipqinp->gspc;
      for (h = 0; h < tj->nthreads; h++) {

        // Determine the pages to write
        imin = h * gdim[0] / tj->nthreads;
        imax = (h == tj->nthreads-1) ? gdim[0] :
                                       (h+1) * gdim[0] / tj->nthreads;

        // Print points for each segment of the grid
        sprintf(inpname, "%s.eval.%d", ipqinp->inpfile, h);
        outp = FOpenSafe(inpname, tj->OverwriteOutput);
        fprintf(outp, "%d\n", (imax-imin) * gdim[1] * gdim[2]);
        origx = (gdim[0] % 2 == 0) ?
          scofm[0] - (gdim[0]/2 - 0.5)*gspc[0] :
          scofm[0] - (gdim[0]/2)*gspc[0];
        origy = (gdim[1] % 2 == 0) ?
          scofm[1] - (gdim[1]/2 - 0.5)*gspc[1] :
          scofm[1] - (gdim[1]/2)*gspc[1];
        origz = (gdim[2] % 2 == 0) ?
          scofm[2] - (gdim[2]/2 - 0.5)*gspc[2] :
          scofm[2] - (gdim[2]/2)*gspc[2];
        for (i = imin; i < imax; i++) {
          dx = (origx + i*gspc[0]) / B2ANG;
          for (j = 0; j < gdim[1]; j++) {
            dy = (origy + j*gspc[1]) / B2ANG;
            for (k = 0; k < gdim[2]; k++) {
              dz = (origz + k*gspc[2]) / B2ANG;
              fprintf(outp, "%10.6lf %10.6lf %10.6lf\n", dx, dy, dz);
            }
          }
        }
        fclose(outp);
      }
      ipqinp->gorig[0] = origx;
      ipqinp->gorig[1] = origy;
      ipqinp->gorig[2] = origz;
    }
    else {

      // Print points for the sites of the charge cloud.  Evaluate the
      // potential at these sites to get the exact energy contribution
      // due to these charges.
      for (h = 0; h < tj->nthreads; h++) {
        sprintf(inpname, "%s.eval.%d", ipqinp->inpfile, h);
        outp = FOpenSafe(inpname, tj->OverwriteOutput);
        fprintf(outp, "%d\n", ipqinp->nQcloud);
        for (i = 0; i < ipqinp->nQcloud; i++) {
          fprintf(outp, "%10.6lf %10.6lf %10.6lf\n",
                  ipqinp->Qcloud.map[i][0] / B2ANG,
                  ipqinp->Qcloud.map[i][1] / B2ANG,
                  ipqinp->Qcloud.map[i][2] / B2ANG);
        }
        fclose(outp);
      }
    }
    // Bail out if there is no solvent charge environment
    return;
  }
  
  outp = FOpenSafe(ipqinp->ptqfile, tj->OverwriteOutput);
  fprintf(outp, "%d\n", ipqinp->nQcloud);
  Qc = &ipqinp->Qcloud;
  for (i = 0; i < ipqinp->nQcloud; i++) {
    fprintf(outp, "%12.9lf  %11.6lf %11.6lf %11.6lf\n", Qc->map[i][3],
            Qc->map[i][0], Qc->map[i][1], Qc->map[i][2]);
  }
  fclose(outp);

}

//-----------------------------------------------------------------------------
// ReadOrcaVpot: when Orca performs the QM calculation, it will create a file
//               that is not in Gaussian cubegen format.  This function will
//               rectify that.                            
//                                                                      
// Arguments:                                                           
//   ipqinp:  IPolQ control parameters                                  
//   envstr:  the environment descriptor string ("vac" or "solv")       
//   tp:      system topology (for atom masses to get atomic numbers)   
//   tj:      trajectory control data (for file overwriting permissions)
//-----------------------------------------------------------------------------
static void ReadOrcaVpot(ipqcon *ipqinp, char* envstr, prmtop *tp, trajcon *tj)
{
  int h, i, j, k, npt, fullpt, idx, idy, idz;
  int *gdim, *itmp;
  double ptx, pty, ptz, dx, dy, dz, origx, origy, origz, val, nrg;
  double *gspc, *dtmp;
  char finame[MAXNAME], tlet[8];
  FILE *inp, *outp;
  dmat Ucloud;
  ibook Aocc;
  dbook A;
  coord *S;

  // If the calculation was a single point evaluation, then a grid is not
  // required and has not, infact, been computed.  Rather, all threads have
  // processed the same set of points corresponding to locations of the
  // charge density cloud.  Process these results into a human-readable
  // file that displays the charges, coordinates, and energies.
  if (ipqinp->DoSinglePoint == 1) {

    // This only happens for the condensed phase calculation
    if (strcmp(envstr, "solv") != 0) {
      return;
    }

    // Read the electrostatic potential evaluation
    sprintf(finame, "%s.%s.0", ipqinp->grdfile, envstr);
    if ((inp = fopen(finame, "r")) == NULL) {
      printf("ReadOrcaVpot >> Error.  QM grid file %s not found.\n", finame);
      exit(1);
    }
    fscanf(inp, "%d", &npt);
    Ucloud = CreateDmat(npt, 4, 0);
    for (i = 0; i < npt; i++) {
      fscanf(inp, "%lf%lf%lf%lf\n", &ptx, &pty, &ptz, &val);
      Ucloud.map[i][0] = ptx * B2ANG;
      Ucloud.map[i][1] = pty * B2ANG;
      Ucloud.map[i][2] = ptz * B2ANG;
      Ucloud.map[i][3] = val * H2KCAL;
    }
    fclose(inp);

    // Reprint the file in a custom format.  It still uses
    // the grid file name from user input.
    S = &ipqinp->Solute;
    sprintf(finame, "%s.uipq", ipqinp->outfile);
    outp = FOpenSafe(finame, tj->OverwriteOutput);
    fprintf(outp, "%% Generated by mdgx potential=%s\n", ipqinp->qmmeth);
    fprintf(outp, "%% Units of these data are Angstroms and kcal/mol-e\n\n");
    fprintf(outp, "%% Z-number (int, real) Crd X(A)    Crd Y       Crd Z\n");
    fprintf(outp, "coords = [\n");
    for (i = 0; i < S->natom; i++) {
      dx = AtomCode(tp->Masses[S->atmid[i]], tp->ZNumber[S->atmid[i]],
		    tp->Charges[S->atmid[i]], tlet);
      idx = dx + 0.01;
      fprintf(outp, "  %5d%12.6lf%12.6lf%12.6lf%12.6lf\n", idx, dx,
              S->loc[3*i], S->loc[3*i+1], S->loc[3*i+2]);
    }
    fprintf(outp, "];\n\n");
    fprintf(outp, "%% Total number of polarizing point charges: %d\n", npt);
    fprintf(outp, "%% Charge (e)       Crd X(A)   Crd Y      Crd Z        "
            "Uelec(QM)\n");
    fprintf(outp, "uelec = [\n");
    nrg = 0.0;
    for (i = 0; i < npt; i++) {
      fprintf(outp, " %13.9lf   %10.6lf %10.6lf %10.6lf %13.6lf\n",
              ipqinp->Qcloud.map[i][3], ipqinp->Qcloud.map[i][0],
              ipqinp->Qcloud.map[i][1], ipqinp->Qcloud.map[i][2],
              Ucloud.map[i][3]);
      nrg += ipqinp->Qcloud.map[i][3] * Ucloud.map[i][3];
    }
    fprintf(outp, "];\n\n");
    fclose(outp);

    // Write the total energy of charges in the
    // QM electrostatic field to mdout
    SpliceFileName(tj, tj->outbase, tj->outsuff, finame, 1);
    outp = fopen(finame, "a");
    fprintf(outp, "\n Potential energy of solvent charges in the QM field: "
            "%12.6lf kcal/mol\n", nrg);
    HorizontalRule(outp, 0);
    fclose(outp);

    return;
  }
  
  origx = ipqinp->gorig[0];
  origy = ipqinp->gorig[1];
  origz = ipqinp->gorig[2];
  gspc = ipqinp->gspc;
  gdim = ipqinp->gdim;
  npt = 0;
  for (h = 0; h < tj->nthreads; h++) {
    sprintf(finame, "%s.%s.%d", ipqinp->grdfile, envstr, h);
    if ((inp = fopen(finame, "r")) == NULL) {
      printf("ReadOrcaVpot >> Error.  QM grid file %s not found.\n", finame);
      exit(1);
    }
    fscanf(inp, "%d", &i);
    npt += i;
    fclose(inp);
  }
  fullpt = ipqinp->gdim[0] * ipqinp->gdim[1] * ipqinp->gdim[2];
  if (npt != fullpt) {
    printf("ReadOrcaVpot >> Expected %d points in the output file(s), found "
           "%d.\n", i, npt);
    exit(1);
  }
  A = CreateDbook(ipqinp->gdim[0], ipqinp->gdim[1], ipqinp->gdim[2], 0);
  Aocc = CreateIbook(ipqinp->gdim[0], ipqinp->gdim[1], ipqinp->gdim[2]);
  for (h = 0; h < tj->nthreads; h++) {
    sprintf(finame, "%s.%s.%d", ipqinp->grdfile, envstr, h);
    inp = fopen(finame, "r");
    fscanf(inp, "%d", &npt);
    for (i = 0; i < npt; i++) {
      fscanf(inp, "%lf%lf%lf%lf\n", &ptx, &pty, &ptz, &val);
      dx = (ptx * B2ANG - origx) / gspc[0] + 0.01;
      dy = (pty * B2ANG - origy) / gspc[1] + 0.01;
      dz = (ptz * B2ANG - origz) / gspc[2] + 0.01;
      idx = dx;
      idy = dy;
      idz = dz;
      A.map[idx][idy][idz] = val;
      if (Aocc.map[idx][idy][idz] == 1) {
        printf("ReadOrcaVpot >> Error.  Two values assigned to grid element "
               "%4d %4d %4d.\n", idx, idy, idz);
        exit(1);
      }
      Aocc.map[idx][idy][idz] += 1;
    }
    fclose(inp);
  }

  // Check to see that all grid points were covered
  for (i = 0; i < gdim[0]; i++) {
    for (j = 0; j < gdim[1]; j++) {
      itmp = Aocc.map[i][j];
      for (k = 0; k < gdim[2]; k++) {
        if (itmp[k] != 1) {
          printf("ReadOrcaVpot >> Error.  Grid point [ %4d %4d %4d ] was not "
                 "calculated.\n", i, j, k);
        }
      }
    }
  }

  // Reprint the file in cubegen format
  S = &ipqinp->Solute;
  sprintf(finame, "%s.%s", ipqinp->grdfile, envstr);
  outp = FOpenSafe(finame, tj->OverwriteOutput);
  fprintf(outp, " Generated by mdgx potential=%s\n", ipqinp->qmmeth);
  fprintf(outp, " Electrostatic potential from Total %s Density\n",
          ipqinp->qmmeth);
  fprintf(outp, "%5d%12.6lf%12.6lf%12.6lf\n", ipqinp->Solute.natom,
          origx / B2ANG, origy / B2ANG, origz / B2ANG);
  fprintf(outp, "%5d%12.6lf%12.6lf%12.6lf\n", gdim[0],
          gspc[0] / B2ANG, 0.0, 0.0);
  fprintf(outp, "%5d%12.6lf%12.6lf%12.6lf\n", gdim[1],
          0.0, gspc[1] / B2ANG, 0.0);
  fprintf(outp, "%5d%12.6lf%12.6lf%12.6lf\n", gdim[2],
          0.0, 0.0, gspc[2] / B2ANG);
  for (i = 0; i < S->natom; i++) {
    dx = AtomCode(tp->Masses[S->atmid[i]], tp->ZNumber[S->atmid[i]],
		  tp->Charges[S->atmid[i]], tlet);
    idx = dx + 0.01;
    fprintf(outp, "%5d%12.6lf%12.6lf%12.6lf%12.6lf\n", idx, dx,
            S->loc[3*i] / B2ANG, S->loc[3*i+1] / B2ANG, S->loc[3*i+2] / B2ANG);
  }
  h = 0;
  for (i = 0; i < A.pag; i++) {
    for (j = 0; j < A.row; j++) {
      dtmp = A.map[i][j];
      for (k = 0; k < A.col; k++) {
        fprintf(outp, " %12.5E", dtmp[k]);
        h++;
        if (h == 6) {
          fprintf(outp, "\n");
          h = 0;
        }
      }
      if (h != 0) {
        h = 0;
        fprintf(outp, "\n");
      }
    }
  }
  fclose(outp);

  // Remove files stored before cubegen format
  for (i = 0; i < tj->nthreads; i++) {
    sprintf(finame, "rm %s.%s.%d", ipqinp->grdfile, envstr, i);
    system(finame);
  }

  // Free allocated memory
  DestroyDbook(&A);
  DestroyIbook(&Aocc);
}

//-----------------------------------------------------------------------------
// ManageQMCalc: this function sets up a QM calculation, manages the system
//               call, and processes the results of the quantum calculation.
//                                                                      
// Arguments:                                                           
//   ipqinp:   IPolQ input data, containing tables of points and charges
//   tp:       the molecular topology (for atom masses and mobilities)  
//   tj:       trajectory control, for output overwriting permissions   
//-----------------------------------------------------------------------------
static void ManageQMCalc(ipqcon *ipqinp, prmtop *tp, trajcon *tj, int env)
{
  int i, slen, scrcreated;
  char envstr[8], rhoext[8];
  char* syscall;
  char* pwd;
  DIR* mydir;
  FILE *inp;
  
  // Determine the current directory so that
  // certain QM outputs may be saved        
  syscall = (char*)malloc(8192*sizeof(char));
  pwd = (char*)malloc(MAXNAME*sizeof(char));
  sprintf(syscall, "echo `pwd` >> %s", ipqinp->finfile);
  system(syscall);
  inp = fopen(ipqinp->finfile, "r");
  fscanf(inp, "%s", pwd);
  fclose(inp);
  sprintf(syscall, "rm -f %s", ipqinp->finfile);
  system(syscall);

  // Prepare quantum input files
  if (strcmp(ipqinp->qmprog, "orca") == 0) {
    WriteOrcaInput(ipqinp, tp, env, tj);
  }
  else if (strcmp(ipqinp->qmprog, "gaussian") == 0) {
    WriteGaussianInput(ipqinp, tp, env, tj);
  }
  if (env == 0) {
    sprintf(envstr, "vacu");
  }
  else {
    sprintf(envstr, "solv");
  }

  // Allocate for the system call and set out
  // lines to prepare for the QM calculation.
  syscall[0] = '\0';
  for (i = 0; i < ipqinp->prepcalls.row; i++) {
    slen = strlen(syscall);
    sprintf(&syscall[slen], "%s; ", ipqinp->prepcalls.map[i]);
  }

  // If a scratch directory has been specified, 
  // set out part of the system call to create it.
  if (ipqinp->scrdir[0] != '\0') {
    mydir = opendir(ipqinp->scrdir);
    if (mydir) {
      scrcreated = 0;
      closedir(mydir);
    }
    else {
      scrcreated = 1;
      slen = strlen(syscall);
      sprintf(&syscall[slen], "mkdir -p %s; ", ipqinp->scrdir);
    }
    if (strcmp(ipqinp->qmprog, "orca") == 0) {
      slen = strlen(syscall);
      sprintf(&syscall[slen], "cp %s.%s %s.eval* %s; ", ipqinp->inpfile,
              envstr, ipqinp->inpfile, ipqinp->scrdir);
      if (env == 1) {
        slen = strlen(syscall);
        sprintf(&syscall[slen], "cp %s %s; ", ipqinp->ptqfile, ipqinp->scrdir);
      }
    }
    else if (strcmp(ipqinp->qmprog, "gaussian") == 0) {
      slen = strlen(syscall);
      sprintf(&syscall[slen], "cp %s.%s %s.eval %s; ", ipqinp->inpfile, envstr,
              ipqinp->inpfile, ipqinp->scrdir);
    }
    slen = strlen(syscall);
    sprintf(&syscall[slen], "cd %s; ", ipqinp->scrdir);
  }
  
  // Part of the system call to run the QM calculation
  slen = strlen(syscall);
  if (strcmp(ipqinp->qmprog, "orca") == 0) {
    sprintf(&syscall[slen], "%s %s.%s > %s.%s; ", ipqinp->qmpath,
            ipqinp->inpfile, envstr, ipqinp->outfile, envstr);
    if (strcmp(ipqinp->qmmeth, "MP2") == 0) {
      sprintf(rhoext, "pmp2re");
    }
    else if (ipqinp->excitation > 0) {
      sprintf(rhoext, "cisp");
    }
    else {
      sprintf(rhoext, "scfp");
    }
    for (i = 0; i < tj->nthreads; i++) {
      slen = strlen(syscall);
      sprintf(&syscall[slen], "%s %s.%s.gbw %s.%s.%s %s.eval.%d %s.%s.%d "
              "> /dev/null & sleep 1; ", ipqinp->uvpath, ipqinp->inpfile,
              envstr, ipqinp->inpfile, envstr, rhoext, ipqinp->inpfile, i,
              ipqinp->grdfile, envstr, i);
    }
  }
  else if (strcmp(ipqinp->qmprog, "gaussian") == 0) {
    sprintf(&syscall[slen], "%s < %s.%s > %s.%s; ", ipqinp->qmpath,
            ipqinp->inpfile, envstr, ipqinp->outfile, envstr);
    slen = strlen(syscall);
    sprintf(&syscall[slen], "%s %s.%s.chk For-%s.%s.chk > /dev/null; %s 0 "
            "potential=%s For-%s.%s.chk %s.%s -1 h < %s.eval > /dev/null; ",
            ipqinp->fmpath, ipqinp->outfile, envstr, ipqinp->outfile, envstr,
            ipqinp->uvpath, ipqinp->qmmeth, ipqinp->outfile, envstr,
            ipqinp->grdfile, envstr, ipqinp->inpfile);
  }
  slen = strlen(syscall);
  sprintf(&syscall[slen], "wait; ");

  // Move files
  if (ipqinp->scrdir[0] != '\0') {
    slen = strlen(syscall);
    sprintf(&syscall[slen], "mv ");
    if (ipqinp->retqminp == 1) {
      slen = strlen(syscall);
      sprintf(&syscall[slen], "%s.%s ", ipqinp->inpfile, envstr);
    }
    if (ipqinp->retqmout == 1) {
      slen = strlen(syscall);
      sprintf(&syscall[slen], "%s.%s ", ipqinp->outfile, envstr);
    }
    if (ipqinp->retqmchk == 1) {
      slen = strlen(syscall);
      if (strcmp(ipqinp->qmprog, "gaussian") == 0) {
        sprintf(&syscall[slen], "%s.%s.chk ", ipqinp->outfile, envstr);
      }
      else if (strcmp(ipqinp->qmprog, "orca") == 0) {
        sprintf(&syscall[slen], "%s.%s.gbw %s.%s.%s ", ipqinp->inpfile,
                envstr, ipqinp->inpfile, envstr, rhoext);
      }
    }
    slen = strlen(syscall);
    sprintf(&syscall[slen], "%s.%s* %s; ", ipqinp->grdfile, envstr, pwd);
  }

  // Remove the scratch directory and point evaluation files
  if (ipqinp->scrdir[0] != '\0' && scrcreated == 1) {
    slen = strlen(syscall);
    sprintf(&syscall[slen], "cd %s; rm -rf %s; ", pwd, ipqinp->scrdir);
  }
  if (env == 1) {
    if (strcmp(ipqinp->qmprog, "orca") == 0 && ipqinp->retptfi == 0) {
      slen = strlen(syscall);
      sprintf(&syscall[slen], "rm %s; ", ipqinp->ptqfile);
    }
    slen = strlen(syscall);
    sprintf(&syscall[slen], "rm %s.eval*; ", ipqinp->inpfile);
  }
  
  // Additional commands in the system call
  // to clean up after the QM calculation     
  for (i = 0; i < ipqinp->postcalls.row; i++) {
    slen = strlen(syscall);
    sprintf(&syscall[slen], "%s; ", ipqinp->postcalls.map[i]);
  }
  system(syscall);

  // Convert file if necessary
  if (strcmp(ipqinp->qmprog, "orca") == 0) {
    ReadOrcaVpot(ipqinp, envstr, tp, tj);
  }

  // Free allocated memory
  free(syscall);
  free(pwd);
}

#ifdef MPI
//-----------------------------------------------------------------------------
// ObserveQMCalc: observe the quantum mechanical calculations managed by the
//                head node.  Check every few seconds for completion of the
//                vacuum quantum calculation, staggering each check to avoid
//                multiple processes all trying to look at the file at exactly
//                the same time.  When the vacuum calculation and electrostatic
//                potential evaluations look to be complete, stop idling.  This
//                prevents the mdgx slave processes from consuming CPU
//                resources while the master is driving other programs.
//
// Arguments:                                                           
//   ipqinp:  IPolQ control parameters (for file names)                 
//   CG:      cell grid, used to track thread IDs and total thread count
//-----------------------------------------------------------------------------
static void ObserveQMCalc(ipqcon *ipqinp, cellgrid *CG)
{
  int qmdone, lnmatch;
  char syscall[64], line[256], outname[MAXNAME];
  FILE *outp;

  // Staggered initial wait
  sprintf(syscall, "sleep %d", CG->tid+1);
  system(syscall);

  // Idle
  qmdone = 0;
  sprintf(syscall, "sleep %d", CG->nthreads + 5);
  sprintf(outname, "%s", ipqinp->finfile);
  while (qmdone == 0) {
    if ((outp = fopen(outname, "r")) == NULL) {
      qmdone = 0;
    }
    else {
      qmdone = 1;
      fclose(outp);
    }
    if (qmdone == 0) {
      system(syscall);
    }
  }
}
#endif

//-----------------------------------------------------------------------------
// ReadQMResult: read the results from a QM calculation and pick out the most
//               important pieces.  This is all in service of IPolQ single
//               point evaluations.
//
// Arguments:
//   fname:    the name of the QM output file
//   qmprog:   the name of the QM executable (not the path, the name)
//   tcrd:     pre-allocated coordinate struct to store forces
//   nrg:      on return, this will store the total QM potential energy
//-----------------------------------------------------------------------------
static void ReadQMResult(char* fname, char* qmprog, coord *tcrd, double *nrg)
{
  int i, j, k, natom, doread, reading, nw;
  cmat FC, lnw;

  FC = Ascii2Mem(fname,  256, 128, "Failed at ReadQMResult.\n");
  if (strcmp(qmprog, "orca") == 0) {
    for (i = 0; i < FC.row; i++) {
      nw = CountWords(FC.map[i]);

      // Check for total energy
      if (nw == 5) {
        lnw = ParseWords(FC.map[i]);
        if (strcmp(lnw.map[0], "FINAL") == 0 &&
            strcmp(lnw.map[1], "SINGLE") == 0 &&
            strcmp(lnw.map[2], "POINT") == 0 &&
            strcmp(lnw.map[3], "ENERGY") == 0 &&
	    WordIsNumber(lnw.map[4]) == 1) {
          *nrg = strtod(lnw.map[4], NULL) * H2KCAL;
        }
        DestroyCmat(&lnw);
      }
    
      // Check for forces
      if (nw == 2) {
        lnw = ParseWords(FC.map[i]);
        doread = (strcmp(lnw.map[0], "CARTESIAN") == 0 &&
                  strcmp(lnw.map[1], "GRADIENT") == 0) ? 1 : 0;	
        DestroyCmat(&lnw);
        if (doread == 1) {
          reading = 0;
          j = i + 1;
          while (j < FC.row && reading == 0) {
            nw = CountWords(FC.map[j]);	    
            if (nw == 6) {
              reading = 1;
            }
            if (reading == 0) {
              j++;
            }
          }
          k = 0;
          while (j < FC.row && reading == 1) {
            nw = CountWords(FC.map[j]);
            if (nw == 6) {
              lnw = ParseWords(FC.map[j]);
              if (WordIsNumber(lnw.map[3]) == 1 &&
                  WordIsNumber(lnw.map[4]) == 1 &&
                  WordIsNumber(lnw.map[5]) == 1) {
                if (k == tcrd->natom) {
                  printf("ReadQMResult >> Error.  Too many atomic forces have "
                         "been read from QM output\nReadQMResult >> %s.\n",
                         fname);
                  exit(1);
                }
                tcrd->frc[3*k    ] = strtod(lnw.map[3], NULL) * H2KCAL;
                tcrd->frc[3*k + 1] = strtod(lnw.map[4], NULL) * H2KCAL;
                tcrd->frc[3*k + 2] = strtod(lnw.map[5], NULL) * H2KCAL;
                k++;
              }
              else {
                reading = 0;
              }
              DestroyCmat(&lnw);
            }
            else {
              reading = 0;
            }
            j++;
          }
        }
      }
    }
  }
  else {
    printf("ReadQMResult >> Error.  Only ORCA outputs are currently supported "
	   "for this\nReadQMResult >> analysis.\n");
    exit(1);
  }
  DestroyCmat(&FC);
}

//-----------------------------------------------------------------------------
// PostProcessQMCalc: post-processing to determine information from the QM
//                    output.  The QM electrostatic potential is read and then
//                    analyzed with respect to the solvent charge density.
//
// Arguments:                                                           
//   ipqinp:  IPolQ control parameters
//   tp:      topology for the solvated system (ipqinp->Solute has indexing
//            into this structure)
//   tj:      trajectory control data (needed for file names)
//   env:     indicator of the environment we are in
//-----------------------------------------------------------------------------
static void PostProcessQMCalc(ipqcon *ipqinp, prmtop *tp, trajcon *tj, int env)
{
  int i, j;
  double nrg, dnrg;
  double scofm[3], tcofm[3];
  char* fname;
  char* gname;
  dmat U;
  dmat *Qc;
  fbook Ue;
  prmtop faketp;
  coord testcrd;
  coord *S;
  FILE *outp;

  // This is always getting done, but "post-processing" means something
  // different when an IPolQ single point calculation occurs.  In that
  // case, rather than read output from the electrostatic potential
  // evaluator (Gaussian cubegen, ORCA's orca_vpot), we read the output
  // from the QM files for both phases.  The energies, forces on all atoms,
  // and point charge locations are all stored in a custom, human-readable
  // format.
  S = &ipqinp->Solute;
  fname = (char*)malloc(MAXNAME*sizeof(char));
  gname = (char*)malloc(MAXNAME*sizeof(char));
  if (ipqinp->DoSinglePoint == 1) {
    testcrd = CreateCoord(S->natom);
    
    // Gas-phase part
    sprintf(gname, "%s.uipq", ipqinp->outfile);
    outp = fopen(gname, "a");
    sprintf(fname, "%s.vacu", ipqinp->outfile);
    ReadQMResult(fname, ipqinp->qmprog, &testcrd, &nrg);
    fprintf(outp, "fVacuum = [\n");
    for (i = 0; i < S->natom; i++) {
      j = LocateResID(tp, S->atmid[i], 0, tp->nres);
      fprintf(outp, "  %12.6lf  %12.6lf  %12.6lf  %% %4.4s %4.4s\n",
	      testcrd.frc[3*i], testcrd.frc[3*i + 1], testcrd.frc[3*i + 2],
	      &tp->AtomNames[4 * S->atmid[i]], &tp->ResNames[4*j]);
    }
    fprintf(outp, "];\n\n");
    fprintf(outp, "eVacuum = [ %14.4lf ];\n\n", nrg);

    // Condensed-phase part
    sprintf(fname, "%s.solv", ipqinp->outfile);
    ReadQMResult(fname, ipqinp->qmprog, &testcrd, &nrg);
    fprintf(outp, "fCondensed = [\n");
    for (i = 0; i < S->natom; i++) {
      j = LocateResID(tp, S->atmid[i], 0, tp->nres);
      fprintf(outp, "  %12.6lf  %12.6lf  %12.6lf  %% %4.4s %4.4s\n",
	      testcrd.frc[3*i], testcrd.frc[3*i + 1], testcrd.frc[3*i + 2],
	      &tp->AtomNames[4 * S->atmid[i]], &tp->ResNames[4*j]);
    }
    fprintf(outp, "];\n\n");
    fprintf(outp, "eCondensed = [ %14.4lf ];\n\n", nrg);
    fclose(outp);

    // Free allocated memory and return
    DestroyCoord(&testcrd);
    free(fname);
    free(gname);
    return;
  }
  
  // Create a fake topology and coordinate buffer for the grid reading
  sprintf(faketp.source, "NONE");
  faketp.natom = S->natom;
  faketp.EPInserted = 0;
  faketp.nxtrapt = 0;
  testcrd = CreateCoord(S->natom);

  // Read the grid
  if (env == 0) {
    sprintf(fname, "%s.vacu", ipqinp->grdfile);
  }
  else if (env == 1) {
    sprintf(fname, "%s.solv", ipqinp->grdfile);
  }
  Ue = ReadEPotGrid(fname, &faketp, &testcrd);
  if (env == 1) {
    FindCoordCenter(S->loc, S->prvloc, 0, S->natom, scofm);
    FindCoordCenter(testcrd.loc, S->prvloc, 0, S->natom, tcofm);
    TransCrd(S->loc, S->natom, scofm, -1.0);
    TransCrd(testcrd.loc, S->natom, tcofm, -1.0);
    U = CreateDmat(3, 3, 0);
    QuatAlign(testcrd.loc, S->loc, S->natom, S->prvloc, 0, &U);
    RotateCrd(S->loc, S->natom, U);
    TransCrd(S->loc, S->natom, tcofm, 1.0);
    TransCrd(testcrd.loc, S->natom, tcofm, 1.0);
    XformChargeCloud(ipqinp, scofm, tcofm, U);
    nrg = 0.0;
    Qc = &ipqinp->Qcloud;
    for (i = 0; i < ipqinp->nQcloud; i++) {
      TriInterp(Qc->map[i], &dnrg, &Ue, 1, 3);
      nrg += 332.0636 * dnrg * Qc->map[i][3];
    }
    SpliceFileName(tj, tj->outbase, tj->outsuff, fname, 1);
    outp = fopen(fname, "a");
    fprintf(outp, "\n Potential energy of solvent charges in the QM field: "
            "%12.6lf kcal/mol\n", nrg);
    HorizontalRule(outp, 0);
    fclose(outp);
  }
  
  // Free allocated memory
  free(fname);
  free(gname);
  DestroyDmat(&U);
  DestroyCoord(&testcrd);
  DestroyFbook(&Ue);
}

//-----------------------------------------------------------------------------
// PerformIPolQ: this function wraps all the necessary features for producing
//               quantum electrostatic potentials from a solute in water
//               inpcrd file.
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
//   ipqinp:  IPolQ control parameters                                  
//-----------------------------------------------------------------------------
void PerformIPolQ(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
                  FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
                  trajcon *tj, Energy *sysUV, execon *etimers, cdftrj* Acdf,
                  ipqcon *ipqinp)
{
  int ncq;
  double* qswap;
  char syscall[MAXNAME];
  FILE *outp;

  // Allocate space for charge array swapping on this process
  qswap = CpyDVec(tp->Charges, tp->natom);

  // Initial setup
  PrepIPolQ(tj, ipqinp, tp, crd, qswap);

  // Dynamics with periodic checks for convergence
  // of the solvent reaction field at solute sites
  if (ipqinp->verbose == 1 && CG->tid == 0) {
    printf("IPolQ >> Sampling solvent charge density.\n");
  }
  while (tj->currstep < tj->nstep) {
    UpdateStepNumber(tj, tj->currstep + 1);
    Dynamics(crd, CG, tp, dcinp, Etab, EHtab, rcinp, PPk, tj, sysUV, etimers,
             Acdf, 0);
    if (tj->currstep >= ipqinp->neqstep && tj->currstep % ipqinp->ntqs == 0) {
      ComputeSRFP(crd, CG, tp, dcinp, Etab, EHtab, rcinp, PPk, etimers,
                  ipqinp, qswap);
      if (ipqinp->verbose == 1) {
#ifdef MPI
        ncq = 0;
        MPI_Reduce(&ipqinp->nQcloud, &ncq, 1, MPI_INT, MPI_SUM, 0,
                   CG->dspcomm);
#else
        ncq = ipqinp->nQcloud;
#endif
        if (CG->tid == 0) {
          fprintf(stderr, "\rIPolQ >> Step %9lld | Frame %5d | %8d pts",
                  tj->currstep, ipqinp->nqframe, ncq);
          fflush(stderr);
        }
      }
    }
  }

  // The electrostatic potential is now accumulated.
  // The stored charge density can be augmented with
  // shell charges to produce a final charge array  
  // that leads into a quantum calculation.         
#ifdef MPI
  PoolCharges(ipqinp, CG);
  if (CG->tid == 0) {
#endif
    CenterQCloud(ipqinp, crd);
    TrimQCloud(ipqinp);
    FitQShell(ipqinp, tj, tp);
#ifdef MPI
  }
#endif

  // The charge array is ready; print input
  // files and run the quantum calculation.
  if (CG->tid == 0) {
    if (ipqinp->verbose == 1) {
      printf("IPolQ >> Running vacuum phase quantum calculation with %s.\n",
             ipqinp->qmprog);
    }
    ManageQMCalc(ipqinp, tp, tj, 0);
    if (ipqinp->verbose == 1) {
      printf("IPolQ >> Running solution phase quantum calculation with %s.\n",
             ipqinp->qmprog);
    }
    ManageQMCalc(ipqinp, tp, tj, 1);
    if (ipqinp->verbose == 1) {
      printf("IPolQ >> Quantum and electrostatic potential calculations "
             "completed.\n");
    }
    PostProcessQMCalc(ipqinp, tp, tj, 1);
    outp = fopen(ipqinp->finfile, "w");
    fprintf(outp, "Done.\n");
    fclose(outp);
  }
#ifdef MPI
  else {
    ObserveQMCalc(ipqinp, CG);
  }
  MPI_Barrier(CG->dspcomm);
#endif
  if (CG->tid == 0) {
    sprintf(syscall, "rm %s", ipqinp->finfile);
    system(syscall);
  }

  // Free allocated memory
  free(qswap);
}
