#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LoopBuilder.h"
#include "ConfigSamp.h"
#include "Matrix.h"
#include "Trajectory.h"
#include "CrdManip.h"
#include "Topology.h"
#include "mdgxVector.h"
#include "Parse.h"
#include "Random.h"

//-----------------------------------------------------------------------------
// InitLoopControl: function to initialize data for a loop prediction or
//                  completion run.
//
// Arguments:
//   lpbinp:     loop control data structure to fill
//-----------------------------------------------------------------------------
void InitLoopControl(loopcon *lpbinp)
{
  lpbinp->maxstretch = 2.0;
  lpbinp->lreconst = 8.0;
  lpbinp->belly = CreateCmat(1, MAXNAME);
  sprintf(lpbinp->belly.map[0], ":*");
}

//-----------------------------------------------------------------------------
// MoveNrg: make an incremental trial move of the points in the plane, then
//          compute the new energy.  This is specifically for two-dimensional
//          systems confined to the XY plane.
//
// Arguments:
//   Tpts:    locations of trial points
//   dpos:    the length of the move to make
//   frc:     forces on all particles / points
//-----------------------------------------------------------------------------
static double MoveNrg(dmat *Tpts, double dpos, dmat *frc)
{
  int i, j;
  double En, dx, dy, r2, r;

  for (i = 0; i < Tpts->row; i++) {
    Tpts->map[i][0] += dpos*frc->map[i][0];
    Tpts->map[i][1] += dpos*frc->map[i][1];
  }
  En = 0.0;
  for (i = 0; i < Tpts->row; i++) {
    for (j = 0; j < i; j++) {
      dx = Tpts->map[j][0] - Tpts->map[i][0];
      dy = Tpts->map[j][1] - Tpts->map[i][1];
      r2 = dx*dx + dy*dy;
      if (r2 < 1.0e-4) {
        r2 = 1.0e-4;
      }
      En += 0.1/sqrt(r2);
    }
    dx = Tpts->map[i][0];
    dy = Tpts->map[i][1];
    r2 = dx*dx + dy*dy;
    if (r2 > 1.0) {
      r = sqrt(r2);
      En += 100.0*(1.0-r)*(1.0-r);
    }
  }

  return En;
}

//-----------------------------------------------------------------------------
// EvenCircle: spread a group of points evenly around and within a circle.  The
//             points are laid out in the XY plane but can later be rotated in
//             three dimensional space.
//
// Arguments:
//   npts:     the number of points to set out
//-----------------------------------------------------------------------------
static dmat EvenCircle(int npts, long *rndcon)
{
  int i, j;
  double dx, dy, r2, r, invr, invr3, fx, fy, fmag;
  double dE, minE, stepval, E0, espline, minstep, maxstep;
  double dpos[4], b[4];
  dmat A, pts, frc, Tpts;

  pts = CreateDmat(npts, 3, 0);
  frc = CreateDmat(npts, 3, 0);
  Tpts = CreateDmat(npts, 3, 0);
  for (i = 0; i < npts; i++) {
    pts.map[i][0] = 2.0*(ran2(rndcon) - 0.5);
    pts.map[i][1] = 2.0*(ran2(rndcon) - 0.5);
    r2 = pts.map[i][0]*pts.map[i][0] + pts.map[i][1]*pts.map[i][1];
    if (r2 > 1.0) {
      invr = 1.0/sqrt(r2);
      pts.map[i][0] *= invr;
      pts.map[i][1] *= invr;
    }
  }
  dE = 1.0;
  dpos[0] = 0.0;
  dpos[1] = 0.001;
  A = CreateDmat(4, 4, 0);
  while (dpos[1] > 1.0e-6) {

    // Compute forces on all 'atoms' and the current energy
    E0 = 0.0;
    SetDVec(frc.data, 3*npts, 0.0);
    for (i = 0; i < npts; i++) {
      for (j = 0; j < i; j++) {
	dx = pts.map[j][0] - pts.map[i][0];
	dy = pts.map[j][1] - pts.map[i][1];
	r2 = dx*dx + dy*dy;
	if (r2 < 1.0e-4) {
	  r2 = 1.0e-4;
	}
	r = sqrt(r2);
	invr = 1.0/r;
	invr3 = invr*invr*invr;
	fx = 0.1*dx*invr3;
	fy = 0.1*dy*invr3;
	frc.map[i][0] -= fx;
	frc.map[i][1] -= fy;
	frc.map[j][0] += fx;
	frc.map[j][1] += fy;	
	E0 += 0.1*invr;
      }
      dx = pts.map[i][0];
      dy = pts.map[i][1];
      r2 = dx*dx + dy*dy;
      if (r2 > 1.0) {
	r = sqrt(r2);
	invr = 1.0/r;
	fx = 200.0*(1.0-r)*dx*invr;
	fy = 200.0*(1.0-r)*dy*invr;
	frc.map[i][0] += fx;
	frc.map[i][1] += fy;
        E0 += 100.0*(1.0-r)*(1.0-r);
      }
    }
    fmag = 0.0;
    for (i = 0; i < npts; i++) {
      fmag += frc.map[i][0]*frc.map[i][0] + frc.map[i][1]*frc.map[i][1];
    }
    fmag = 1.0/sqrt(fmag);
    for (i = 0; i < npts; i++) {
      frc.map[i][0] *= fmag;
      frc.map[i][1] *= fmag;
    }
    b[0] = E0;

    // Make moves and compute new energies
    for (i = 0; i < npts; i++) {
      Tpts.map[i][0] = pts.map[i][0];
      Tpts.map[i][1] = pts.map[i][1];
    }
    b[1] = MoveNrg(&Tpts, dpos[1], &frc);
    dpos[2] = (b[1] > b[0]) ? 0.5*dpos[1] : 2.0*dpos[1];
    b[2] = MoveNrg(&Tpts, dpos[2]-dpos[1], &frc);
    if (dpos[2] < dpos[1]) {
      dpos[3] = (b[2] > b[1]) ? 0.25*dpos[1] : 0.75*dpos[1];
    }
    else {
      dpos[3] = (b[2] > b[1]) ? 1.5*dpos[1] : 3.0*dpos[1];
    }
    b[3] = MoveNrg(&Tpts, dpos[3]-dpos[2], &frc);
    
    // Solve the cubic polynomial and find the minimum value
    maxstep = 0.0;
    for (i = 0; i < 4; i++) {
      A.map[i][0] = dpos[i]*dpos[i]*dpos[i];
      A.map[i][1] = dpos[i]*dpos[i];
      A.map[i][2] = dpos[i];
      A.map[i][3] = 1.0;
      if (dpos[i] > maxstep) {
	maxstep = dpos[i];
      }
    }
    AxbQRRxc(A, b, 0);
    BackSub(A, b);
    minE = b[3];
    minstep = 0.0;
    for (i = 0; i < 200; i++) {
      stepval = 0.005*maxstep*i;
      espline = stepval*(stepval*(b[0]*stepval + b[1]) + b[2]) + b[3];
      if (espline < minE) {
	minE = espline;
	minstep = stepval;
      }
    }
    dpos[1] = (minstep > dpos[1]*1.05) ? dpos[1]*1.05 : dpos[1]*0.95;
    for (i = 0; i < npts; i++) {
      pts.map[i][0] = pts.map[i][0] + minstep*frc.map[i][0];
      pts.map[i][1] = pts.map[i][1] + minstep*frc.map[i][1];
    }
    dE = E0 - minE;
  }

  // Free allocated memory
  DestroyDmat(&pts);
  DestroyDmat(&Tpts);
  DestroyDmat(&A);

  return frc;
}

//-----------------------------------------------------------------------------
// CalcPath: function for calculating the path in x, y, and z that connects a
//           series of points as a function of some fourth variable t.
//
// Arguments:
//   lpbinp:   the loop building control data (contains information about the
//             path size, etc.)
//   pidx:     index of the path to connect within lpbinp
//   dim:      index of the dimension to operate along (0 = x, 1 = y, 2 = z)
//   pathC:    scratch matrix for coefficients of the path splines in x, y, or
//             z (each direction is handled individually)
//   pathV:    coordinates of the path in x, y, or z (each direction is handled
//             separately in different function calls)
//-----------------------------------------------------------------------------
static void CalcPath(track *path, int dim, dmat *pathC, double *pathV)
{
  int i, nk;
  double di;
  double *ptmp;
  
  // The equations are: value at a knot, then derivative at that knot,
  // repeating for every knot.  The path makes a loop, requiring that
  // the value and derivative at the end points are identical.
  nk = path->nknot; 
  ptmp = path->crd;
  SetDVec(pathC->data, nk*nk, 0.0);
  for (i = 0; i < nk; i++) {
    di = i;
    
    // Value of path function at t = i
    pathC->map[i][3*i] = di*di;
    pathC->map[i][3*i+1] = di;
    pathC->map[i][3*i+2] = 1.0;
    pathV[i] = ptmp[3*i+dim];

    // Value of path function at t = i+1
    pathC->map[i+1][3*i] = (di+1.0)*(di+1.0);
    pathC->map[i+1][3*i+1] = di+1.0;
    pathC->map[i+1][3*i+2] = 1.0;
    pathV[i+1] = (i == nk-1) ? ptmp[dim] : ptmp[3*(i+1)+dim];

    // Value of path derivative at t = i
    pathC->map[i+2][3*i] = 2.0*di;
    pathC->map[i+2][3*i+1] = 1.0;
    pathC->map[i+2][3*i+2] = 0.0;
    if (i == nk-1) {
      pathC->map[i+2][2] = -1.0;
    }
    else {
      pathC->map[i+2][3*(i+1)] = -2.0*di;
      pathC->map[i+2][3*(i+1)+1] = -1.0;
    }
    pathV[i+2] = 0.0;
  }

  // Solve the matrix equation
  AxbQRRxc(*pathC, pathV, 0);
  BackSub(*pathC, pathV);

  // Store the results
  if (dim == 0) {
    ReflectDVec(path->xcoef, pathV, 3*nk);
  }
  else if (dim == 1) {
    ReflectDVec(path->ycoef, pathV, 3*nk);
  }
  else if (dim == 2) {
    ReflectDVec(path->zcoef, pathV, 3*nk);
  }
}

//-----------------------------------------------------------------------------
// TraceConnection: if there is a long bond in a structure, there could be many
//                  different results based on the way in which it eventually
//                  connects.  Different results could be achieved by pulling
//                  the bond together over different paths.  This routine will
//                  outline a series of paths to try.
//
// Arguments:
//   
//-----------------------------------------------------------------------------
static void TraceConnections(loopcon *lpbinp, long *rndcon, imat *LB,
			     dmat *crdpop)
{
  int i, j, k, nlb, nnp;
  double disp;
  double dr[3], dv[3], zvec[3], endpt[3];
  double* pathV;
  dmat EC, U, ECrt, pathC;
  
  // Make the stencil of points
  EC = EvenCircle(lpbinp->npath, rndcon);
  ECrt = CreateDmat(3, 3, 0);

  // Prepare scratch space for path spline computations
  pathC = CreateDmat(2*lpbinp->pathlength-2, 2*lpbinp->pathlength-2, 0);
  pathV = (double*)malloc((2*lpbinp->pathlength-2)*sizeof(double));
  
  // Line up layers of the stencil along the axis
  // between the ends of the long bond
  nlb = 0;
  for (i = 0; i < LB->row; i++) {
    for (j = 0; j < 3; j++) {
      endpt[j] = crdpop->map[3*LB->map[i][0]+j][0];
      dr[j] = crdpop->map[3*LB->map[i][1]+j][0] - endpt[j];
    }
    if (sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]) < lpbinp->lreconst) {
      nlb++;
    }
  }
  lpbinp->paths = (track*)malloc(lpbinp->npath*nlb*sizeof(track));
  lpbinp->nlongbond = nlb;
  nlb = 0;
  for (i = 0; i < LB->row; i++) {
    for (j = 0; j < 3; j++) {
      endpt[j] = crdpop->map[3*LB->map[i][0]+j][0];
      dr[j] = crdpop->map[3*LB->map[i][1]+j][0] - endpt[j];
    }
    if (sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]) < lpbinp->lreconst) {
      continue;
    }

    // This is a long bond and paths can now be traced that have a chance to
    // connect the two ends.  Make a set of paths that will then be used as
    // tracks to bring the two ends of the bond together.
    zvec[0] = 0.0;
    zvec[1] = 0.0;
    zvec[2] = 1.0;
    U = RotateA2B(zvec, dr);
    SetDVec(ECrt.data, 9, 0.0);
    for (j = 0; j < lpbinp->npath; j++) {
      ECrt.map[j][0] = U.data[0]*EC.map[j][0] + U.data[1]*EC.map[j][1] +
	               U.data[2]*EC.map[j][2] + endpt[0];
      ECrt.map[j][1] = U.data[0]*EC.map[j][0] + U.data[1]*EC.map[j][1] +
	               U.data[2]*EC.map[j][2] + endpt[1];
      ECrt.map[j][2] = U.data[0]*EC.map[j][0] + U.data[1]*EC.map[j][1] +
	               U.data[2]*EC.map[j][2] + endpt[2];
    }
    disp = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]) / lpbinp->pathlength;
    nnp = lpbinp->npath*nlb;

    // Lay out the points for each path.
    for (j = 0; j < lpbinp->npath; j++) {
      lpbinp->paths[j+nnp].atmA = LB->map[i][0];
      lpbinp->paths[j+nnp].atmB = LB->map[i][1];
      for (k = 0; k < 3; k++) {
        lpbinp->paths[j+nnp].crd[k] = endpt[k];
        lpbinp->paths[j+nnp].crd[3*lpbinp->npath-3+k] = endpt[k] + dr[k];
      }
    }
    for (j = 1; j < lpbinp->pathlength-1; j++) {
      dv[0] = j*disp*dr[0];
      dv[1] = j*disp*dr[1];
      dv[2] = j*disp*dr[2];
      for (k = 0; k < lpbinp->npath; k++) {
	lpbinp->paths[k+nnp].crd[j] = ECrt.map[k][0] + dv[0];
	lpbinp->paths[k+nnp].crd[j+1] = ECrt.map[k][1] + dv[1];
	lpbinp->paths[k+nnp].crd[j+2] = ECrt.map[k][2] + dv[2];
      }
    }

    // Find the piecewise splines defining each path
    for (j = 0; j < lpbinp->npath; j++) {
      for (k = 0; k < 3; k++) {
        CalcPath(&lpbinp->paths[j+nnp], k, &pathC, pathV);
      }
    }
    nlb++;
  }

  // Free allocated memory
  DestroyDmat(&EC);
  DestroyDmat(&ECrt);
  DestroyDmat(&pathC);
  free(pathV);
}

//-----------------------------------------------------------------------------
// FindLongBonds: find all long bonds within a coordinate set, defined as those
//                which are stretched beyond a certain limit.
//
// Arguments:
//   lpbinp:  loop building controls
//   tc:      the coordinate set to look at
//   tp:      the topology to check
//-----------------------------------------------------------------------------
static imat FindLongBonds(loopcon *lpbinp, coord *tc, prmtop *tp)
{
  int i, j, nlb, maxlb;
  double dx, dy, dz;
  bondcomm *thisbc;
  imat LB;

  maxlb = 32;
  LB = CreateImat(maxlb, 2);
  nlb = 0;
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->BLC[i].nbond; j++) {
      thisbc = &tp->BLC[i].BC[j];
      if (lpbinp->movable[thisbc->a] + lpbinp->movable[thisbc->b] == 0) {
	continue;
      }
      dx = tc->loc[3*thisbc->b] - tc->loc[3*thisbc->a];
      dy = tc->loc[3*thisbc->b+1] - tc->loc[3*thisbc->a+1];
      dz = tc->loc[3*thisbc->b+2] - tc->loc[3*thisbc->a+2];
      if (dx*dx + dy*dy + dz*dz >
	  tp->BParam[thisbc->t].l0 + lpbinp->maxstretch) {
	LB.map[nlb][0] = thisbc->a;
	LB.map[nlb][1] = thisbc->b;
	nlb++;
	if (nlb >= maxlb) {
	  maxlb += 32;
	  LB = ReallocImat(&LB, maxlb, 2);
	}
      }
    }
  }
  LB = ReallocImat(&LB, nlb, 2);

  return LB;
}

//-----------------------------------------------------------------------------
// DetectMobileAtoms: this function is much like its counterpart in the
//                    configuration sampling module, but it works off of a
//                    different struct.
//
// Arguments:
//   lpbinp:     loop building input data
//   tp:         the system topology
//   tc:         system coordinates (not a transposed array, a coord struct)
//-----------------------------------------------------------------------------
static void DetectMobileAtoms(loopcon *lpbinp, prmtop *tp, coord *tc)
{
  int i;
  int* amask;
  
  lpbinp->movable = (int*)calloc(tp->natom, sizeof(int));
  for (i = 0; i < lpbinp->nbelly; i++) {
    amask = ParseAmbMask(lpbinp->belly.map[i], tp, tc);
    IVec2VecAdd(lpbinp->movable, amask, tp->natom);
    free(amask);
  }
  lpbinp->allmove = (ISum(lpbinp->movable, tp->natom) == tp->natom) ? 1 : 0;
}

//-----------------------------------------------------------------------------
// LoopBuilder: a function for creating loops, or more precisely adjusting the
//              conformations of pieces of a larger structure to minimize the
//              overall energy.  This routine invokes functions from the
//              configuration sampling module in a script-like fashion, with
//              its own set of inputs.
//-----------------------------------------------------------------------------
void LoopBuilder(loopcon *lpbinp, trajcon *tj, prmtop *tp)
{
  int i, nactive, nconf;
  int* mobile;
  double tstamp;
  imat LB;
  dmat crdpop, frcpop;
  coord tc;

  // As it was with configuration sampling, the topology will not have been
  // read in main by the time the program dives into this routine.  Read the
  // topology now.
  GetPrmTop(tp, tj, 1);

  // Initialize the coordinate population.  This will be done using a single
  // set of coordinates.
  tstamp = 0.0;
  tc = ReadRst(tp, tj->ipcname.map[0], &tstamp);
  crdpop = CreateDmat(3*tp->natom, lpbinp->count, 0);
  for (i = 0; i < 3*tp->natom; i++) {
    SetDVec(crdpop.map[i], lpbinp->count, tc.loc[i]);
  }
  frcpop = CreateDmat(3*tp->natom, lpbinp->count, 0);
  
  // Search for long bonds and find mobile atoms before
  // destroying the original coord struct
  DetectMobileAtoms(lpbinp, tp, &tc);
  LB = FindLongBonds(lpbinp, &tc, tp);
  DestroyCoord(&tc);
  TraceConnections(lpbinp, &tj->rndcon, &LB, &crdpop);
  
  // The topology will have already been read in, as well as the
  // original set of coordinates.  Allocate for the specified number of
  // trial conformations (the true number could be much higher, but
  // only a certain number will be actively worked on at once).
  nconf = lpbinp->count;

  // Compute the energy of the original structure once in order
  // to understand where there might be significant strain.  Construct
  // bond splints to heal any bonds that are terribly stretched.

  // Once the inital structure is reasonably minimized, make all of
  // the other conformations copies of that.  Proceed to minimize
  // batches of conformations toward particular targets.  A target
  // consists of one or more desired locations for particular atoms
  // as well as backbone dihedral conformations.  This continues
  // until all targets have been hit from the initial conformation
  // or a previously computed solution from a nearby target. 
  // Promising solutions, whether low in energy or simply adequate
  // and far from previously encountered solutions, will be saved.

  // Make additional attempts to reach each target from previous
  // successes.  Keep the lowest energy solution to any one target.

  // Report the best results.
}
