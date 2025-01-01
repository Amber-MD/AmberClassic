#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mdgxVector.h"
#include "CrdManip.h"
#include "Constraints.h"
#include "CellManip.h"
#include "Topology.h"
#include "Debug.h"

#include "CrdManipDS.h"
#include "TrajectoryDS.h"
#include "pmeDirectDS.h"
#include "pmeRecipDS.h"

//-----------------------------------------------------------------------------
// Vec2PrimeC: function needed by settlec.  Essentially computes a
//             matrix-vector multiplication, where n1, n2, and n0 are the rows
//             of the matrix.                                  
//
// Arguments:                                                           
//   V:      the vector to transform                                    
//   n[120]: the rows of the transformation matrix                      
//-----------------------------------------------------------------------------
static void Vec2PrimeC(double* V, double* n1, double* n2, double* n0)
{
  double Vp[3];

  Vp[0] = Dot3(n1, V);
  Vp[1] = Dot3(n2, V);
  Vp[2] = Dot3(n0, V);
  V[0] = Vp[0];
  V[1] = Vp[1];
  V[2] = Vp[2];
}

//-----------------------------------------------------------------------------
// settlec: the first stage of SETTLE (analytic bond length constraints for
//          rigid, three-point water molecules).                    
//                                                                      
// This function is adapted from code developed by the Theoretical and  
// Computational Biophysics Group in the Beckman Institute for Advanced 
// Science and Technology at the University of Illinois at              
// Urbana-Champaign.                                                    
//                                                                      
// Arguments:                                                           
//   ref:     the initial coordinates of all three atoms after the most recent
//            application of constraints, in the order oxygen, hydrogen 1,
//            hydrogen 2                                    
//   m[OH]:   the masses of the oxygen (O) and hydrogen (H) atoms       
//   pos:     on input, the perturbed coordinates (as after a position update,
//            prompted by forces acting on the atoms); on output, the
//            constrained coordinates                       
//   r[abc]:  the canonical positions of water atoms (see the diagram in
//            Constraints.h)                                            
//-----------------------------------------------------------------------------
static void settlec(const double* ref, double* pos, settleparm *FWtab)
{
  int i;
  double sinphi, cosphi, sinpsi, cospsi, rbphi, tmp, tmp1, tmp2;
  double alpha, beta, gamma, a2b2, sintheta, costheta, mO, mH, ra, rb, rc;
  double b0[3], c0[3], d0[3], a1[3], b1[3], c1[3], a2[3], b2[3], c2[3];
  double a3[3], b3[3], c3[3], m0[3], m1[3], m2[3], n0[3], n1[3], n2[3];

  // Unpack SETTLE parameters
  mO = FWtab->mO;
  mH = FWtab->mH;
  ra = FWtab->ra;
  rb = FWtab->rb;
  rc = FWtab->rc;

  for (i = 0; i < 3; i++) {

    // Vectors in the plane of the original positions
    b0[i] = ref[3+i]-ref[i];
    c0[i] = ref[6+i]-ref[i];

    // Compute and store the new center of mass
    d0[i] = (pos[i]*mO + pos[3+i]*mH + pos[6+i]*mH)/(mO+mH+mH);

    // Translate the system's center of mass to the origin
    a1[i] = pos[i] - d0[i];
    b1[i] = pos[3+i] - d0[i];
    c1[i] = pos[6+i] - d0[i];
  }

  // Vectors describing transformation from the original coordinate
  // system to the 'primed' coordinate system as in the diagram.   
  CrossP(b0, c0, n0);
  CrossP(a1, n0, n1);
  CrossP(n0, n1, n2);
  UnitVector3(n0);
  UnitVector3(n1);
  UnitVector3(n2);
  Vec2PrimeC(b0, n1, n2, n0);
  Vec2PrimeC(c0, n1, n2, n0);
  Vec2PrimeC(a1, n1, n2, n0);
  Vec2PrimeC(b1, n1, n2, n0);
  Vec2PrimeC(c1, n1, n2, n0);

  // Compute positions of canonical water
  sinphi = a1[2]/ra;
  tmp = 1.0-sinphi*sinphi;
  cosphi = sqrt(tmp);
  sinpsi = (b1[2] - c1[2])/(2.0*rc*cosphi);
  tmp = 1.0-sinpsi*sinpsi;
  cospsi = sqrt(tmp);

  rbphi = -rb*cosphi;
  tmp1 = rc*sinpsi*sinphi;
  tmp2 = rc*sinpsi*cosphi;
 
  a2[0] = 0.0;
  a2[1] = ra*cosphi;
  a2[2] = ra*sinphi;
  b2[0] = -rc*cospsi;
  b2[1] = rbphi - tmp1;
  b2[2] = -rb*sinphi + tmp2;
  c2[0] = rc*cosphi;
  c2[1] = rbphi+tmp1;
  c2[2] = -rb*sinphi - tmp2;

  // There are no a0 terms because a0 was subtracted off
  // when b0 and c0 were first defined.                 
  alpha = b2[0]*(b0[0] - c0[0]) + b0[1]*b2[1] + c0[1]*c2[1];
  beta  = b2[0]*(c0[1] - b0[1]) + b0[0]*b2[1] + c0[0]*c2[1];
  gamma  = b0[0]*b1[1] - b1[0]*b0[1] + c0[0]*c1[1] - c1[0]*c0[1];
 
  a2b2 = alpha*alpha + beta*beta;
  sintheta = (alpha*gamma - beta*sqrt(a2b2 - gamma*gamma))/a2b2;
  costheta = sqrt(1.0 - sintheta*sintheta);

  a3[0] = -a2[1]*sintheta;
  a3[1] = a2[1]*costheta;
  a3[2] = a1[2];
  b3[0] = b2[0]*costheta - b2[1]*sintheta;
  b3[1] = b2[0]*sintheta + b2[1]*costheta;
  b3[2] = b1[2];
  c3[0] = -b2[0]*costheta - c2[1]*sintheta;
  c3[1] = -b2[0]*sintheta + c2[1]*costheta;
  c3[2] = c1[2];

  // Undo the transformation; generate new normal vectors
  // from the transpose.                                 
  m1[0] = n1[0];
  m1[1] = n2[0];
  m1[2] = n0[0];
  m2[0] = n1[1];
  m2[1] = n2[1];
  m2[2] = n0[1];
  m0[0] = n1[2];
  m0[1] = n2[2];
  m0[2] = n0[2];
  Vec2PrimeC(a3, m1, m2, m0);
  Vec2PrimeC(b3, m1, m2, m0);
  Vec2PrimeC(c3, m1, m2, m0);
 
  // Translate the system's center of mass back to its proper position
  for (i = 0; i < 3; i++) {
    pos[i] = a3[i] + d0[i];
    pos[3+i] = b3[i] + d0[i];
    pos[6+i] = c3[i] + d0[i];
  }
}

//-----------------------------------------------------------------------------
// settlev: velocity constraints portion of the SETTLE algorithm.       
//                                                                      
// This function is adapted from code developed by the Theoretical and  
// Computational Biophysics Group in the Beckman Institute for Advanced 
// Science and Technology at the University of Illinois at              
// Urbana-Champaign.                                                    
//                                                                      
// Arguments:                                                           
//   pos:    the constrained positions of the three atoms (as output by 
//           settlec, above)                                            
//   m[OH]:  masses of oxygen (O) and hydrogen (H) atoms                
//   vel:    on input, the unconstrained velocities, given in the same order
//           as the atoms (O, H1, H2); on output, the constrained velocities
//-----------------------------------------------------------------------------
static void settlev(const double* pos, double* vel, settleparm *FWtab)
{
  int i;
  double cosA, cosB, cosC, vab, vbc, vca, mOH, d, tab, tbc, tca, mO, mH;
  double AB[3], BC[3], CA[3], vtAB[3], vtBC[3], vtCA[3], ga[3], gb[3], gc[3];

  // Unpack SETTLE parameters
  mO = FWtab->mO;
  mH = FWtab->mH;

  for (i = 0; i < 3; i++) {
    AB[i] = pos[3+i]-pos[i];
    BC[i] = pos[6+i]-pos[3+i];
    CA[i] = pos[i]-pos[6+i];
  }

  UnitVector3(AB);
  UnitVector3(BC);
  UnitVector3(CA);

  cosA = -Dot3(AB, CA);
  cosB = -Dot3(BC, AB);
  cosC = -Dot3(CA, BC);

  for (i = 0; i < 3; i++) {
    vtAB[i] = (vel[3+i]-vel[i]);
    vtBC[i] = (vel[6+i]-vel[3+i]);
    vtCA[i] = (vel[i]-vel[6+i]);
  }
  vab = Dot3(AB, vtAB);
  vbc = Dot3(BC, vtBC);
  vca = Dot3(CA, vtCA);

  mOH = mO+mH;
  d = mH/(mOH*mOH + mH*cosA*(mO*cosB*cosC - mH*cosA) -
	  0.5*mO*mOH*(cosB*cosB + cosC*cosC));
  tab = (vab*(2.0*mOH - mO*cosC*cosC) + 
	 vbc*(mH*cosC*cosA - mOH*cosB) +
	 vca*(mO*cosB*cosC - 2.0*mH*cosA))*mO*d;
  tbc = (vbc*(mOH*mOH - mH*mH*cosA*cosA) +
	 vca*mO*(mH*cosA*cosB - mOH*cosC) +
	 vab*mO*(mH*cosC*cosA - mOH*cosB))*d;
  tca = (vca*(2*mOH - mO*cosB*cosB) +
	 vab*(mO*cosB*cosC - 2*mH*cosA) +
	 vbc*(mH*cosA*cosB - mOH*cosC))*mO*d;

  mO = 0.5/mO;
  mH = 0.5/mH;
  for (i = 0; i < 3; i++) {
    ga[i] = tab*AB[i] - tca*CA[i];
    gb[i] = tbc*BC[i] - tab*AB[i];
    gc[i] = tca*CA[i] - tbc*BC[i];
    vel[i] += mO*ga[i];
    vel[3+i] += mH*gb[i];
    vel[6+i] += mH*gc[i];
  }
}

#define CNSTTYPE 0
#include "ConstrainPos.c"
#undef CNSTTYPE

#define CNSTTYPE 1
#include "ConstrainPos.c"
#undef CNSTTYPE

//-----------------------------------------------------------------------------
// LoadVelCnstResult: this routine makes the decision whether to enter the
//                    results of the application of velocity constraints
//                    directly or prepare them for export to another cell which
//                    actually controls the atoms.  There are a lot of
//                    arguments to this little function, but it is likely to be
//                    inlined by a compiler and makes for much simpler coding
//                    as many types of velocity constraints can make use of it.
//
// Arguments:
//   C:        the cell which applied the constraint                    
//   expsctr:  8-elemenet array, zero if the sector does NOT need to be 
//             exported and 1 if it does                                
//   atmreg:   the region of the cell in which the atom resides         
//   atmid:    the ID number of the atom                                
//   atmvel:   the constrained velocity computed for this atom          
//   sysvel:   pointer to the location of the atom's velocities in the  
//             associated coord struct (atmvel gets put into sysvel if  
//             expsctr[atmreg] is 0)                                    
//-----------------------------------------------------------------------------
static void LoadVelCnstResult(cell *C, int* expsctr, int atmreg, int atmid,
			      double *atmvel, double *sysvel)
{
  int nexp;

  // This atom will need to be exported
  if (expsctr[atmreg]) {
    nexp = C->nexp;
    C->Vexport[nexp].id = atmid;
    C->Vexport[nexp].dreg = atmreg;
    C->Vexport[nexp].vel[0] = atmvel[0];
    C->Vexport[nexp].vel[1] = atmvel[1];
    C->Vexport[nexp].vel[2] = atmvel[2];
    C->nexp += 1;
  }

  // This atom is controlled by the cell
  else {
    sysvel[0] = atmvel[0];
    sysvel[1] = atmvel[1];
    sysvel[2] = atmvel[2];
  }
}

//-----------------------------------------------------------------------------
// RattleGroupVel: velocity constraints stage of RATTLE algorithm for rigid
//                 bonds in the context of Velocity-Verlet integration.
//
// This function is adapted from code in the tinker software package and
// the AMBER sff program.                                               
//
// Arguments:                                                           
//   CG:       the cell grid                                            
//   crd:      the coordinates (this is where velocity information      
//             residues)                                                
//   tj:       trajectory control information                           
//   tp:       the topology                                             
//   cloc:     location of the cell of interest                         
//   myCC:     vector containing codified information about the atoms   
//             and connectivity in this RATTLE group                    
//   cellidx:  pre-allocated array to hold the indices of all atoms in  
//             this constraint group, indices into the cell's data array
//   xref:     pre-allocated array holding velocities of atoms          
//   xpos:     pre-allocated array holding coordinates of atoms         
//   xposc:    pre-allocated array holding copy of xpos to compare after
//             constraints are applied                                  
//-----------------------------------------------------------------------------
static void RattleGroupVel(coord *crd, trajcon *tj, prmtop *tp, cell *C,
			   int* myCC, int* cellidx, double* xpos, double* xvel,
			   int* expsctr)
{
  int i, j, i3, ididx, atma, atmb, atma3, atmb3, nnode;
  int done, niter;
  double l0, rx, ry, rz, vx, vy, vz, dot;
  double rma, rmb, term, xterm, yterm, zterm;
  atomc *Cdidx;

  // First, this function must find all of
  // the atoms in the constraint group.   
  ididx = 3*myCC[0]+2;
  nnode = myCC[ididx-1];
  for (i = 0; i < nnode; i++) {
    atma = myCC[ididx+i];
    cellidx[i] = C->GPSptr[atma];
    i3 = 3*i;
    Cdidx = &C->data[cellidx[i]];
    for (j = 0; j < 3; j++) {
      xpos[i3+j] = Cdidx->loc[j];
      xvel[i3+j] = crd->vel[3*atma+j];
    }
  }

  // Apply the RATTLE constraints to velocities
  const double rtoldt = tj->rattletol / (sqrt(418.4)*tj->dt);
  niter = 0;
  done = 0;
  const int maxiter = tj->MaxRattleIter;
  while (done == 0 && niter < maxiter) {
    done = 1;
    for (i = 0; i < myCC[0]; i++) {
      atma = myCC[3*i+1];
      atmb = myCC[3*i+2];
      atma3 = 3*atma;
      atmb3 = 3*atmb;
      rx = xpos[atmb3] - xpos[atma3];
      ry = xpos[atmb3+1] - xpos[atma3+1];
      rz = xpos[atmb3+2] - xpos[atma3+2];
      vx = xvel[atmb3] - xvel[atma3];
      vy = xvel[atmb3+1] - xvel[atma3+1];
      vz = xvel[atmb3+2] - xvel[atma3+2];
      dot = rx * vx + ry * vy + rz * vz;
      rma = tp->InvMasses[myCC[ididx + atma]];
      rmb = tp->InvMasses[myCC[ididx + atmb]];
      l0 = myCC[3*i+3]/1.0e8;
      term = -1.2 * dot / ((rma + rmb) * l0 * l0);
      if (fabs(term) > rtoldt) {
	done = 0;
	xterm = rx * term;
	yterm = ry * term;
	zterm = rz * term;
	xvel[atma3] -= xterm*rma;
	xvel[atma3+1] -= yterm*rma;
	xvel[atma3+2] -= zterm*rma;
	xvel[atmb3] += xterm*rmb;
	xvel[atmb3+1] += yterm*rmb;
	xvel[atmb3+2] += zterm*rmb;
      }
    }
    niter++;
  }

  // Test to make sure that the maximum  
  // number of iterations was not reached
  if (niter == maxiter && done == 0) {
    printf("RattleGroupVel >> Error.  Maximum RATTLE iterations exceeded.\n\n"
	   "Group = [\n");
    for (i = 0; i < myCC[0]; i++) {

      // The use of the variables atma and atmb changes here.
      // The program is about to exit in error, so there is  
      // no risk of confusing the variables later on.        
      atma = myCC[ididx + myCC[3*i+1]];
      atmb = myCC[ididx + myCC[3*i+2]];
      printf("%4d %.4s %.4s -> %4d %.4s %.4s\n", atma,
	     &tp->ResNames[4*LocateResID(tp, atma, 0, tp->nres)],
	     &tp->AtomNames[4*atma], atmb,
	     &tp->ResNames[4*LocateResID(tp, atmb, 0, tp->nres)],
	     &tp->AtomNames[4*atmb]);
    }
    printf("];\n");
    exit(1);
  }

  // Prepare to export the velocity adjustments
  for (i = 0; i < nnode; i++) {
    atma = myCC[ididx+i];
    LoadVelCnstResult(C, expsctr, cellidx[i]/C->maxatom, atma, &xvel[3*i],
		      &crd->vel[3*atma]);
  }
}

//-----------------------------------------------------------------------------
// CellVelocityCnst: this is a wrapper for SETTLE, RATTLE, and other constraint
//                   functions in a framework of cells.  The function utilizes
//                   velocity information that should already be present for
//                   relevant atoms when the function is called.
//
// Arguments:
//   CG:    the cell grid (this is where position information and assignments
//          as to which cell will perform the constraint comes from)
//   crd:   the coordinates (this is where velocity information resides)
//   tp:    the topology                                                
//   cloc:  location of the cell of interest                            
//-----------------------------------------------------------------------------
static void CellVelocityCnst(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj,
			     int cgid)
{
  int h, i, j, k, ni, nj, nk, atmid, atmid3, nh1, nh1id, nh1reg, nh1id3;
  int nh2, nh2id, nh2reg, nh2id3, ncpu;
  int expsctr[8];
  int* cellidx;
  double* xpos;
  double* xvel;
  double *ccen, *dbng;
  double cloc[3];
  atomc *atmhi, *atmnh1, *atmnh2;
  cnstcomm *tSHL;
  cell *C;

  // Pointers
  C = &CG->data[cgid];
  cloc[0] = C->gbin[0];
  cloc[1] = C->gbin[1];
  cloc[2] = C->gbin[2];
  tSHL = tp->SHL;
  ccen = &CG->celldim[4];
  dbng = CG->dbng;

  // Pre-allocate arrays for storing this
  // group's coordinates and velocities  
  xpos = (double*)malloc(tp->RattleGrpMax*sizeof(double));
  xvel = (double*)malloc(tp->RattleGrpMax*sizeof(double));
  cellidx = (int*)malloc(tp->RattleGrpMax/3*sizeof(int));

  // Determine which sectors will require export;
  h = 0;
  ncpu = C->CGRank;
  for (k = cloc[2]; k < cloc[2]+2; k++) {
    nk = (k == CG->ng[2]) ? 0 : k;
    for (j = cloc[1]; j < cloc[1]+2; j++) {
      nj = (j == CG->ng[1]) ? 0 : j;
      for (i = cloc[0]; i < cloc[0]+2; i++) {
	ni = (i == CG->ng[0]) ? 0 : i;
	expsctr[h] = (CG->map[ni][nj][nk].CGRank != ncpu);
	h++;
      }
    }
  }

  // Loop over all eight cell sectors and 
  // find atoms in the central region that
  // own constraint groups.               
  C->nexp = 0;
  for (h = 0; h < 8; h++) {
    for (i = 0; i < C->nr[h]; i++) {
      atmhi = &C->map[h][i];
      if (tSHL[atmhi->id].exe == 1 &&
	  IsCentralAtom(C->orig, ccen, atmhi->loc, crd->U.data,
			dbng, crd->isortho) >= 0) {
	atmid = atmhi->id;
        nh1id = tSHL[atmid].blist[0];
        nh2id = tSHL[atmid].blist[1];
	nh1 = C->GPSptr[nh1id];
	nh2 = C->GPSptr[nh2id];
	nh1reg = nh1 / C->maxatom;
	nh2reg = nh2 / C->maxatom;
	atmnh1 = &C->data[nh1];
	atmnh2 = &C->data[nh2];
	atmid3 = 3*atmid;
	nh1id3 = 3*nh1id;
	nh2id3 = 3*nh2id;
        for (j = 0; j < 3; j++) {
          xpos[j] = atmhi->loc[j];
          xpos[3+j] = atmnh1->loc[j];
          xpos[6+j] = atmnh2->loc[j];
          xvel[j] = crd->vel[atmid3+j];
          xvel[3+j] = crd->vel[nh1id3+j];
          xvel[6+j] = crd->vel[nh2id3+j];
        }

	// Because all of the positions come from atoms    
	// within the cell struct, no re-imaging is needed.
        settlev(xpos, xvel, &tp->FWtab);

	// Add this atom to the export table, with velocities in 
	// the place of the position field.  The Vexport table   
	// is also employed to keep track of which sector of this
	// cell each atom with a velocity update started in.     
	LoadVelCnstResult(C, expsctr, h, atmid, xvel, &crd->vel[atmid3]);
	LoadVelCnstResult(C, expsctr, nh1reg, nh1id, &xvel[3],
			  &crd->vel[nh1id3]);
	LoadVelCnstResult(C, expsctr, nh2reg, nh2id, &xvel[6],
			  &crd->vel[nh2id3]);
      }
      else if (tSHL[atmhi->id].exe == 2 &&
	       IsCentralAtom(C->orig, ccen, atmhi->loc,
			     crd->U.data, dbng, crd->isortho) >= 0) {
	RattleGroupVel(crd, tj, tp, C, tSHL[atmhi->id].blist, cellidx, xpos,
		       xvel, expsctr);
      }
    }
  }

  // Free allocated memory
  free(xpos);
  free(xvel);
  free(cellidx);
}

//-----------------------------------------------------------------------------
// BroadcastVelocities: this function broadcasts velocity updates due to
//                      constraints.  After a position update and subsequent
//                      atom imports, cells know velocities for all atoms in
//                      all their sectors.  Most of the velocity information
//                      is in fact unnecessary, but in some cases it will be
//                      updated by constraints and then only certain cells
//                      will know the true velocities of certain atoms.  This
//                      routine will share that knowledge appropriately.
//
// Arguments:
//   CG:      the cell grid                                             
//   cellid:  the location of the cell in the cell grid                 
//   imove:   the direction of the move (0 for +X, 1 for +Y, 2 for +Z)  
//   nmsg:    the number of the process-to-process send                 
//-----------------------------------------------------------------------------
static void BroadcastVelocities(cellgrid *CG, int cellid, int imove, int nmsg)
{
  int i, oreg, aexp, iC1, iC2, iC3, iC4, iD1, iD2, iD3, iD4;
  int dloc[3];
  cell *C, *D;
  atomb *texport;

  C = &CG->data[cellid];
  dloc[0] = C->gbin[0];
  dloc[1] = C->gbin[1];
  dloc[2] = C->gbin[2];
  dloc[imove] = (C->gbin[imove] == CG->ng[imove]-1) ? 0 : C->gbin[imove]+1;
  D = &CG->map[dloc[0]][dloc[1]][dloc[2]];

  // Determine the sector-to-sector transfers
  if (imove == 0) {
    iC1 = 1;    iC2 =  3;    iC3 =  5;    iC4 =  7;
    iD1 = 0;    iD2 =  2;    iD3 =  4;    iD4 =  6;
  }
  else if (imove == 1) {
    iC1 = 2;    iC2 = -3;    iC3 =  6;    iC4 = -7;
    iD1 = 0;    iD2 =  1;    iD3 =  4;    iD4 =  5;
  }
  else {
    iC1 = 4;    iC2 = -5;    iC3 = -6;    iC4 = -7;
    iD1 = 0;    iD2 =  1;    iD3 =  2;    iD4 =  3;
  }

  // Scan all atoms whose velocities were changed 
  // in C and load the appropriate ones for export
  aexp = 1;
  texport = C->pexport;
  C->pexport = (C->CGRank == D->CGRank) ? D->import :
    &CG->pexport[nmsg][CG->nexp[nmsg]];
  for (i = 0; i < C->nexp; i++) {
    oreg = C->Vexport[i].dreg;
    if (oreg == iC1 || oreg == iC2 || oreg == iC3 || oreg == iC4) {
      C->pexport[aexp].id = C->Vexport[i].id;
      C->pexport[aexp].loc[0] = C->Vexport[i].vel[0];
      C->pexport[aexp].loc[1] = C->Vexport[i].vel[1];
      C->pexport[aexp].loc[2] = C->Vexport[i].vel[2];
      C->pexport[aexp].dreg = (oreg == iC1) ? iD1 : (oreg == iC2) ? iD2 :
	(oreg == iC3) ? iD3 : iD4;
      aexp++;
    }
  }
  C->pexport[0].id = aexp;

  // Replace the export pointer in C
  C->pexport = texport;
  if (C->CGRank != D->CGRank) {
    CG->nexp[nmsg] += aexp;
  }
}

//-----------------------------------------------------------------------------
// UnpackVelocities: unpacks the velocities transferred during velocity
//                   constraint updates.
//
// Arguments: 
//   D:       destination cell into which velocities have been delivered
//   crd:     the system coordiantes, at least to the extent they are known by
//            this process
//-----------------------------------------------------------------------------
static void UnpackVelocities(cell *D, coord *crd)
{
  int i, atmid, aexp;
  double *vtmp;

  // Update new velocities in D
  vtmp = crd->vel;
  for (i = 1; i < D->import[0].id; i++) {

    // If this atom is in the cell's primary sector, update  
    // the velocity officially.  Otherwise, it's just passing
    // through, so contribute it to the Vexport buffer.      
    if (D->import[i].dreg == 0) {
      atmid = D->import[i].id;
      vtmp[3*atmid] = D->import[i].loc[0];
      vtmp[3*atmid+1] = D->import[i].loc[1];
      vtmp[3*atmid+2] = D->import[i].loc[2];
    }
    else {

      // Here, aexp refers not to the number of atomb structs
      // that D is about to export but instead to the number 
      // of atombv structs that D has accumulated and WILL   
      // filter for export on the next BroadcastVelocities().
      // D->nexp was initialized by CellPositionCnst() above,
      // and holds the number of atoms in constraint groups  
      // whose velocities have been changed.  D->nexp is NOT 
      // padded by 1 as it is not exported.                  
      aexp = D->nexp;
      D->Vexport[aexp].id = D->import[i].id;
      D->Vexport[aexp].dreg = D->import[i].dreg;
      D->Vexport[aexp].vel[0] = D->import[i].loc[0];
      D->Vexport[aexp].vel[1] = D->import[i].loc[1];
      D->Vexport[aexp].vel[2] = D->import[i].loc[2];
      D->nexp = aexp + 1;
    }
  }
}

//-----------------------------------------------------------------------------
// C2CProcessCnstUpdate: after constraints have been applied, there may be
//                       position and velocity updates in each cell's Vexport
//                       array that will need to be transferred selectively to 
//                       neighboring cells by the export array.  The Vexport
//                       array stores data on all modified coordinates which
//                       the cell does not own, because the Vexport array is
//                       needed as a buffer for specifc exports to each
//                       receiving cell; velocity constraint broadcasting makes
//                       a similar use of the next larger data type.  This
//                       function handles communications from each cell to its
//                       seven neighbors in need of updates.  In the position
//                       constraint form, ???? is "Cnst" and Positions are
//                       transferred as relative moves only, obviating the need
//                       for re-imaging.  In the rescaling form, ???? is "Rscl"
//                       and absolute positions are communicated.
//
// Arguments:
//   CG:           the cell grid                                        
//   cloc:         the location of the cell in the cell grid            
//   imove:        the dimension in which the move is occuring          
//   nmsg:         the number of the process-to-process send message    
//-----------------------------------------------------------------------------
void C2CProcessCnstUpdate(cellgrid *CG, int cellid, int imove, int nmsg)
{
  int i, oreg, aexp;
  int iC1, iC2, iC3, iC4, iD1, iD2, iD3, iD4;
  int dloc[3];
  cell *C, *D;
  atomb *cRdata;
  atombv *cVdata;

  C = &CG->data[cellid];
  dloc[0] = C->gbin[0];
  dloc[1] = C->gbin[1];
  dloc[2] = C->gbin[2];
  dloc[imove] = (C->gbin[imove] == CG->ng[imove]-1) ? 0 : C->gbin[imove]+1;
  D = &CG->map[dloc[0]][dloc[1]][dloc[2]];

  // Determine the sector-to-sector transfers
  if (imove == 0) {
    iC1 = 1;    iC2 = -3;    iC3 = -5;    iC4 = -7;
    iD1 = 0;    iD2 = -2;    iD3 = -4;    iD4 = -6;
  }
  else if (imove == 1) {
    iC1 = 2;    iC2 =  3;    iC3 = -6;    iC4 = -7;
    iD1 = 0;    iD2 =  1;    iD3 = -4;    iD4 = -5;
  }
  else {
    iC1 = 4;    iC2 =  5;    iC3 =  6;    iC4 =  7;
    iD1 = 0;    iD2 =  1;    iD3 =  2;    iD4 =  3;
  }

  // Filter the export data; transfer directly to the D cell  
  // import region if both cells are owned by the same process
  cVdata = C->Vexport;
  cRdata = (C->CGRank == D->CGRank) ? D->import :
    &CG->pexport[nmsg][CG->nexp[nmsg]];
  aexp = 1;
  for (i = 0; i < C->nexp; i++) {
    oreg = cVdata[i].dreg;
    if (oreg == iC1 || oreg == iC2 || oreg == iC3 || oreg == iC4) {
      cRdata[aexp].id = cVdata[i].id;
      cRdata[aexp].loc[0] = cVdata[i].loc[0];
      cRdata[aexp].loc[1] = cVdata[i].loc[1];
      cRdata[aexp].loc[2] = cVdata[i].loc[2];
      cRdata[aexp].dreg = (oreg == iC1) ? iD1 : (oreg == iC2) ? iD2 :
        (oreg == iC3) ? iD3 : iD4;
      aexp++;
    }

  }
  cRdata[0].id = aexp;

  // If this is a process-to-process message,
  // update the pooled export count          
  if (C->CGRank != D->CGRank) {
    CG->nexp[nmsg] += aexp;
  }
}

//-----------------------------------------------------------------------------
// ApplyGridCnst: this is a wrapper for Cell(Velocity/Position)Cnst, which in
//                turn wrap the SETTLE and RATTLE constraint functions.  The
//                style of action is determined by the app flag as described
//                below.
//
// Arguments:
//   CG:    the cell grid (this is where positions come from)           
//   crd:   the coordinates (this is where velocities come from)        
//   tp:    the topology                                                
//   tj:    trajectory control information (contains the time step, MPI 
//          maps, and MPI data type definitions)                        
//   Mcut:  the maximum direct space cutoff                             
//   app:   the application (0 = position constraints, 1 = velocity     
//          constraints, 2 = position rescaling)                        
//   chi:   vector of rescaling factors (calls to this function in the  
//          InitVelocities() and Dynamics() routines pass the rarely    
//          used prvfrc field of a coord struct to this function, as the
//          chi argument is not used by those function calls at all)    
//-----------------------------------------------------------------------------
void ApplyGridCnst(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj,
		   double Mcut, int app, double* chi)
{
  int i, j, imove, istr, ilim, iinc;
  double invdt;
  ashr *cshr;

#ifdef MPI
  MPI_Request* req;
  MPI_Status* stt;
  req = (MPI_Request*)malloc((CG->nsend + CG->nrecv)*sizeof(MPI_Request));
  stt = (MPI_Status*)malloc((CG->nsend + CG->nrecv)*sizeof(MPI_Status));
#endif

  // Bail right out if there are no constraints
  // and this is not a rescaling               
  if (app < 2 && tp->settle == 0 && tp->rattle == 0) {
    return;
  }

  // In the case of position constraint application, the  
  // positions of atoms in primary sectors have just been 
  // updated and must be shared amongst neighboring cells.
  // In the case of velocity constraint application,      
  // current velocities are known only for atoms in the   
  // primary sector and therefore must be shared.  In the 
  // case of coordinate rescaling, current positions of   
  // atoms in ALL sectors are already known, so no sharing
  // is needed.                                           
  if (app < 2) {
#ifdef MPI
    ShareCoordinates(CG, tj, Mcut, crd, tp, 2-app);
#else
    ShareCoordinates(CG, Mcut, crd, tp, 2-app);
#endif
  }

  // Update the atom GPS
  UpdateCellGPS(CG);
  for (i = 0; i < CG->MyCellCount; i++) {

    // Apply constraints
    if (app == 0) {
      CellPositionCnst(CG, CG->MyCellDomain[i], crd, tp, tj);
    }
    else if (app == 1) {
      CellVelocityCnst(CG, crd, tp, tj, CG->MyCellDomain[i]);
    }
    else if (app == 2) {
      CellPositionRscl(CG, CG->MyCellDomain[i], crd, tp, chi);
    }
  }

  // Positions or velocities of atoms in all cells have been
  // updated by constraints, but the updates are only valid 
  // in the cells which control the atoms.  Export tables   
  // were accumulated for each cell, so make the appropriate
  // data transfers to apply all constraints in the primary 
  // sectors of the atoms' home cells.                      
  if (app == 0 || app == 2) {
    istr = 2;
    iinc = -1;
    ilim = -1;
    invdt = 1.0/(sqrt(418.4)*tj->dt);
  }
  else {
    istr = 0;
    iinc = 1;
    ilim = 3;
  }
  for (imove = istr; imove != ilim; imove += iinc) {
    cshr = &CG->DirCommPlan.frcmg[imove];

#ifdef MPI
    // Post receives
    for (i = 0; i < cshr->nrecv; i++) {
      MPI_Irecv(CG->import[i], CG->maximp[i], tj->MPI_ATOMB,
		cshr->recv[i].partner,
		cshr->recv[i].BaseID + DSP_CONSTRAIN + app,
		CG->dspcomm, &req[i]);
    }
    int nreq = cshr->nrecv;
#endif

    // Post sends
    for (i = 0; i < cshr->nsend; i++) {
      CG->nexp[i] = 0;
      if (app == 0 || app == 2) {
	for (j = 0; j < cshr->send[i].ncell; j++) {
	  C2CProcessCnstUpdate(CG, cshr->send[i].cellpt[j], imove, i);
	}
      }
      else if (app == 1) {
	for (j = 0; j < cshr->send[i].ncell; j++) {
	  BroadcastVelocities(CG, cshr->send[i].cellpt[j], imove, i);
	}
      }
#ifdef MPI
      if (cshr->send[i].partner != CG->tid) {
        MPI_Isend(CG->pexport[i], CG->nexp[i], tj->MPI_ATOMB,
                  cshr->send[i].partner,
		  cshr->send[i].BaseID + DSP_CONSTRAIN + app,
                  CG->dspcomm, &req[nreq]);
        nreq++;
      }
#endif
    }

#ifdef MPI
    // Wait on sends
    if (nreq > 0) {
      MPI_Waitall(nreq, req, stt);
    }
    MPI_Barrier(CG->dspcomm);

    // Unpack
    for (i = 0; i < cshr->nrecv; i++) {
      int nimp = 0;
      atomb *timport;
      for (j = 0; j < cshr->recv[i].ncell; j++) {
        cell *D;
        D = &CG->data[cshr->recv[i].cellpt[j]];
        timport = D->import;
        D->import = &CG->import[i][nimp];
	if (app == 0) {
	  UnpackProcessCnstUpdate(D, crd, invdt);
	}
	else if (app == 1) {
	  UnpackVelocities(D, crd);
	}
	else if (app == 2) {
	  UnpackProcessRsclUpdate(D, crd);
	}
        nimp += D->import[0].id;
        D->import = timport;
      }
    }
#endif
    for (i = 0; i < cshr->selfrecv.ncell; i++) {
      if (app == 0) {
	UnpackProcessCnstUpdate(&CG->data[cshr->selfrecv.cellpt[i]], crd,
				invdt);
      }
      else if (app == 1) {
	UnpackVelocities(&CG->data[cshr->selfrecv.cellpt[i]], crd);
      }
      else if (app == 2) {
        UnpackProcessRsclUpdate(&CG->data[cshr->selfrecv.cellpt[i]], crd);
      }
    }
  }

  // If the application is for positions, atom migration may occur
  if (app == 0 || app == 2) {

    // At this point, position and velocity updates due to the 
    // constraints have been applied across cells.  However,   
    // atoms may have migrated out of their primary cells.  It 
    // is therefore necessary to have another round of updates,
    // just as it was at the end of CellVerletC.               
#ifdef MPI
    UpdateCells(CG, crd, tp, tj);
#else
    UpdateCells(CG, crd, tp);
#endif
  }

#ifdef MPI
  // Free allocated memory
  free(req);
  free(stt);
#endif
}

//-----------------------------------------------------------------------------
// VVConstraintVirial: compute the virial contributions from the SETTLE bond
//                     constraints in either the first or second
//                     Velocity-Verlet integration step.
//                
// Arguments:                                                           
//   crd:     the coordinates (for velocity information)                
//   tp:      the topology (for masses)                                 
//   tj:      trajectory control information                            
//   sysUV:   energy and virial information                             
//-----------------------------------------------------------------------------
void VVConstraintVirial(coord *crd, prmtop *tp, trajcon *tj, Energy *sysUV)
{
  int i, j, k, atjk3;
  int atnum[24];
  double invdt, acc, jmass;
  double* oldcrd;
  double* newcrd;
  double* oldvel;
  double* newvel;

  // Shortcuts
  invdt = (tj->dt > 1.0e-12) ? 1.0 / (sqrt(418.4)*tj->dt) : 0.0;
  newcrd = crd->loc;
  oldcrd = crd->prvloc;
  newvel = crd->vel;
  oldvel = crd->prvvel;

  // If this is the first step, there are special considerations
  if (tj->currstep == 0) {

    // Loop over all residues
    for (i = 0; i < tp->natom; i++) {

      // Does this atom control a SETTLE constraint group?
      if (tp->SHL[i].exe == 1) {
	atnum[0] = i;
	atnum[1] = tp->SHL[i].blist[0];
	atnum[2] = tp->SHL[i].blist[1];

	// Loop over the oxygen and both hydrogen atoms
	for (j = 0; j < 3; j++) {
	  jmass = tp->Masses[atnum[j]];
	  for (k = 0; k < 3; k++) {
	    atjk3 = 3*atnum[j]+k;
            acc = (newcrd[atjk3] - oldcrd[atjk3])*invdt - oldvel[atjk3];
            sysUV->Vir[4*k] -= acc*oldcrd[atjk3]*jmass;
          }
        }
      }

      // Does this atom control other constraint groups?
    }
  }
  else {

    // Loop over all residues
    for (i = 0; i < tp->natom; i++) {

      // Does this atom control a SETTLE constraint group?
      if (tp->SHL[i].exe == 1) {
	atnum[0] = i;
	atnum[1] = tp->SHL[i].blist[0];
	atnum[2] = tp->SHL[i].blist[1];

	// Loop over the oxygen and both hydrogen atoms
	for (j = 0; j < 3; j++) {
	  jmass = tp->Masses[atnum[j]];
	  for (k = 0; k < 3; k++) {
	    atjk3 = 3*atnum[j]+k;
            acc = newvel[atjk3] - oldvel[atjk3];
            sysUV->Vir[4*k] -= acc*newcrd[atjk3]*jmass;
          }
        }
      }

      // Does this atom control other constraint groups?
    }
  }
}
