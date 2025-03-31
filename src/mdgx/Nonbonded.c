#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Constants.h"
#include "Macros.h"
#include "SpecialMath.h"
#include "mdgxVector.h"
#include "Nonbonded.h"

#include "CellManipDS.h"
#include "pmeDirectDS.h"
#include "CompFrcDS.h"
#include "TopologyDS.h"
#include "TrajectoryDS.h"

//-----------------------------------------------------------------------------
// TestBondedExclusion: this function tests whether two atoms share a bonded
//                      exclusion by checking the topology.  This is a fairly
//                      time-consuming process, but it only happens if the
//                      distance between two atoms is extremely small.
//
// Arguments:                                                           
//   a1:  the first atom                                                
//   a2:  the second atom                                               
//   tp:  the topology                                                  
//-----------------------------------------------------------------------------
int TestBondedExclusion(int a1, int a2, prmtop *tp)
{
  int i, llim, hlim, isexcl;
  int* exatm;

  // Make sure that a1 < a2 for reference in the exclusion list
  if (a1 > a2) {
    SWAP(a1, a2, i);
  }

  // Determine limits and loop
  llim = tp->ConExcl[a1];
  hlim = tp->ConExcl[a1+1];
  exatm = &tp->ExclList[llim];
  isexcl = 0;
  for (i = 0; i < hlim-llim; i++) {
    if (exatm[i] == a2) {
      isexcl = 1;
      break;
    }
  }
  
  return isexcl;
}

//-----------------------------------------------------------------------------
// CalcAtomR2: calculate r2 between one atom and a list of others, then store
//             the results for later use.                         
//
// Arguments:
//   lim2:      the number of atoms that we are dealing with (the arrays ordr
//              and atmJ are passed in as pointers from the start positions of
//              the atom sequence in the event that this is a loop over
//              interactions within a cell's primary region)
//   ordr:      the order in which to access atoms from atmJ            
//   atmJ:      the atom candidates for interactions                    
//   p[x,y,z]:  the xyz coordinates of the atom being interacted with   
//   IDbuff:    buffer for positive result atom IDs into the cell list  
//   r2buff:    buffer for positive result displacements                
//   cut2:      the squared cutoff                                      
//-----------------------------------------------------------------------------
static int CalcAtomR2(int lim2, int *ordr, atomc *atmJ, const double px,
		      const double py, const double pz, int* IDbuff,
		      rngbuff* r2buff, double cut2)
{
  int j, npair;
  int nj[8];
  double dx[8], dy[8], dz[8], r2[8];

  const int l8lim2 = 8*(lim2/8);
  npair = 0;
  for (j = 0; j < l8lim2; j += 8) {

    // Compute eight r2 calculations
    nj[0] = ordr[j];
    dx[0] = atmJ[nj[0]].loc[0] - px;
    dy[0] = atmJ[nj[0]].loc[1] - py;
    dz[0] = atmJ[nj[0]].loc[2] - pz;
    r2[0] = dx[0]*dx[0] + dy[0]*dy[0] + dz[0]*dz[0];
    nj[1] = ordr[j+1];
    dx[1] = atmJ[nj[1]].loc[0] - px;
    dy[1] = atmJ[nj[1]].loc[1] - py;
    dz[1] = atmJ[nj[1]].loc[2] - pz;
    r2[1] = dx[1]*dx[1] + dy[1]*dy[1] + dz[1]*dz[1];
    nj[2] = ordr[j+2];
    dx[2] = atmJ[nj[2]].loc[0] - px;
    dy[2] = atmJ[nj[2]].loc[1] - py;
    dz[2] = atmJ[nj[2]].loc[2] - pz;
    r2[2] = dx[2]*dx[2] + dy[2]*dy[2] + dz[2]*dz[2];
    nj[3] = ordr[j+3];
    dx[3] = atmJ[nj[3]].loc[0] - px;
    dy[3] = atmJ[nj[3]].loc[1] - py;
    dz[3] = atmJ[nj[3]].loc[2] - pz;
    r2[3] = dx[3]*dx[3] + dy[3]*dy[3] + dz[3]*dz[3];
    nj[4] = ordr[j+4];
    dx[4] = atmJ[nj[4]].loc[0] - px;
    dy[4] = atmJ[nj[4]].loc[1] - py;
    dz[4] = atmJ[nj[4]].loc[2] - pz;
    r2[4] = dx[4]*dx[4] + dy[4]*dy[4] + dz[4]*dz[4];
    nj[5] = ordr[j+5];
    dx[5] = atmJ[nj[5]].loc[0] - px;
    dy[5] = atmJ[nj[5]].loc[1] - py;
    dz[5] = atmJ[nj[5]].loc[2] - pz;
    r2[5] = dx[5]*dx[5] + dy[5]*dy[5] + dz[5]*dz[5];
    nj[6] = ordr[j+6];
    dx[6] = atmJ[nj[6]].loc[0] - px;
    dy[6] = atmJ[nj[6]].loc[1] - py;
    dz[6] = atmJ[nj[6]].loc[2] - pz;
    r2[6] = dx[6]*dx[6] + dy[6]*dy[6] + dz[6]*dz[6];
    nj[7] = ordr[j+7];
    dx[7] = atmJ[nj[7]].loc[0] - px;
    dy[7] = atmJ[nj[7]].loc[1] - py;
    dz[7] = atmJ[nj[7]].loc[2] - pz;
    r2[7] = dx[7]*dx[7] + dy[7]*dy[7] + dz[7]*dz[7];

    // Store positive results
    if (r2[0] < cut2) {
      IDbuff[npair] = nj[0];
      r2buff[npair].dx = dx[0];
      r2buff[npair].dy = dy[0];
      r2buff[npair].dz = dz[0];
      r2buff[npair].r2 = r2[0];
      npair++;
    }
    if (r2[1] < cut2) {
      IDbuff[npair] = nj[1];
      r2buff[npair].dx = dx[1];
      r2buff[npair].dy = dy[1];
      r2buff[npair].dz = dz[1];
      r2buff[npair].r2 = r2[1];
      npair++;
    }
    if (r2[2] < cut2) {
      IDbuff[npair] = nj[2];
      r2buff[npair].dx = dx[2];
      r2buff[npair].dy = dy[2];
      r2buff[npair].dz = dz[2];
      r2buff[npair].r2 = r2[2];
      npair++;
    }
    if (r2[3] < cut2) {
      IDbuff[npair] = nj[3];
      r2buff[npair].dx = dx[3];
      r2buff[npair].dy = dy[3];
      r2buff[npair].dz = dz[3];
      r2buff[npair].r2 = r2[3];
      npair++;
    }
    if (r2[4] < cut2) {
      IDbuff[npair] = nj[4];
      r2buff[npair].dx = dx[4];
      r2buff[npair].dy = dy[4];
      r2buff[npair].dz = dz[4];
      r2buff[npair].r2 = r2[4];
      npair++;
    }
    if (r2[5] < cut2) {
      IDbuff[npair] = nj[5];
      r2buff[npair].dx = dx[5];
      r2buff[npair].dy = dy[5];
      r2buff[npair].dz = dz[5];
      r2buff[npair].r2 = r2[5];
      npair++;
    }
    if (r2[6] < cut2) {
      IDbuff[npair] = nj[6];
      r2buff[npair].dx = dx[6];
      r2buff[npair].dy = dy[6];
      r2buff[npair].dz = dz[6];
      r2buff[npair].r2 = r2[6];
      npair++;
    }
    if (r2[7] < cut2) {
      IDbuff[npair] = nj[7];
      r2buff[npair].dx = dx[7];
      r2buff[npair].dy = dy[7];
      r2buff[npair].dz = dz[7];
      r2buff[npair].r2 = r2[7];
      npair++;
    }
  }
  if (j <= lim2-4) {

    // Compute remaining r2 calculations
    nj[0] = ordr[j];
    dx[0] = atmJ[nj[0]].loc[0] - px;
    dy[0] = atmJ[nj[0]].loc[1] - py;
    dz[0] = atmJ[nj[0]].loc[2] - pz;
    r2[0] = dx[0]*dx[0] + dy[0]*dy[0] + dz[0]*dz[0];
    nj[1] = ordr[j+1];
    dx[1] = atmJ[nj[1]].loc[0] - px;
    dy[1] = atmJ[nj[1]].loc[1] - py;
    dz[1] = atmJ[nj[1]].loc[2] - pz;
    r2[1] = dx[1]*dx[1] + dy[1]*dy[1] + dz[1]*dz[1];
    nj[2] = ordr[j+2];
    dx[2] = atmJ[nj[2]].loc[0] - px;
    dy[2] = atmJ[nj[2]].loc[1] - py;
    dz[2] = atmJ[nj[2]].loc[2] - pz;
    r2[2] = dx[2]*dx[2] + dy[2]*dy[2] + dz[2]*dz[2];
    nj[3] = ordr[j+3];
    dx[3] = atmJ[nj[3]].loc[0] - px;
    dy[3] = atmJ[nj[3]].loc[1] - py;
    dz[3] = atmJ[nj[3]].loc[2] - pz;
    r2[3] = dx[3]*dx[3] + dy[3]*dy[3] + dz[3]*dz[3];

    // Store positive results
    if (r2[0] < cut2) {
      IDbuff[npair] = nj[0];
      r2buff[npair].dx = dx[0];
      r2buff[npair].dy = dy[0];
      r2buff[npair].dz = dz[0];
      r2buff[npair].r2 = r2[0];
      npair++;
    }
    if (r2[1] < cut2) {
      IDbuff[npair] = nj[1];
      r2buff[npair].dx = dx[1];
      r2buff[npair].dy = dy[1];
      r2buff[npair].dz = dz[1];
      r2buff[npair].r2 = r2[1];
      npair++;
    }
    if (r2[2] < cut2) {
      IDbuff[npair] = nj[2];
      r2buff[npair].dx = dx[2];
      r2buff[npair].dy = dy[2];
      r2buff[npair].dz = dz[2];
      r2buff[npair].r2 = r2[2];
      npair++;
    }
    if (r2[3] < cut2) {
      IDbuff[npair] = nj[3];
      r2buff[npair].dx = dx[3];
      r2buff[npair].dy = dy[3];
      r2buff[npair].dz = dz[3];
      r2buff[npair].r2 = r2[3];
      npair++;
    }
    j += 4;
  }
  if (j <= lim2-2) {

    // Compute remaining r2 calculations
    nj[0] = ordr[j];
    dx[0] = atmJ[nj[0]].loc[0] - px;
    dy[0] = atmJ[nj[0]].loc[1] - py;
    dz[0] = atmJ[nj[0]].loc[2] - pz;
    r2[0] = dx[0]*dx[0] + dy[0]*dy[0] + dz[0]*dz[0];
    nj[1] = ordr[j+1];
    dx[1] = atmJ[nj[1]].loc[0] - px;
    dy[1] = atmJ[nj[1]].loc[1] - py;
    dz[1] = atmJ[nj[1]].loc[2] - pz;
    r2[1] = dx[1]*dx[1] + dy[1]*dy[1] + dz[1]*dz[1];

    // Store positive results
    if (r2[0] < cut2) {
      IDbuff[npair] = nj[0];
      r2buff[npair].dx = dx[0];
      r2buff[npair].dy = dy[0];
      r2buff[npair].dz = dz[0];
      r2buff[npair].r2 = r2[0];
      npair++;
    }
    if (r2[1] < cut2) {
      IDbuff[npair] = nj[1];
      r2buff[npair].dx = dx[1];
      r2buff[npair].dy = dy[1];
      r2buff[npair].dz = dz[1];
      r2buff[npair].r2 = r2[1];
      npair++;
    }
    j += 2;
  }
  if (j < lim2) {

    // Compute remaining r2 calculation
    nj[0] = ordr[j];
    dx[0] = atmJ[nj[0]].loc[0] - px;
    dy[0] = atmJ[nj[0]].loc[1] - py;
    dz[0] = atmJ[nj[0]].loc[2] - pz;
    r2[0] = dx[0]*dx[0] + dy[0]*dy[0] + dz[0]*dz[0];
    
    // Store positive result
    if (r2[0] < cut2) {
      IDbuff[npair] = nj[0];
      r2buff[npair].dx = dx[0];
      r2buff[npair].dy = dy[0];
      r2buff[npair].dz = dz[0];
      r2buff[npair].r2 = r2[0];
      npair++;
    }
  }

  return npair;
}

//-----------------------------------------------------------------------------
// CullDualRange: further cull the list in the event that there are two cutoffs
//                (the squared electrostatic cutoff Ecut2 is required to be
//                less than or equal to the Lennard-Jones cutoff Vcut2).  The
//                results are stored in qIDbuff and qr2buff; ljIDbuff and
//                ljr2buff are unchanged.         
//                                                                      
// Arguments:                                                           
//   nljpair:   the number of Lennard-Jones successes reported from     
//              CalcAtomR2 above                                        
//   Ecut2:     the squared electrostatic cutoff                        
//   ljIDbuff:  buffer of Lennard-Jones success IDs into the cell list  
//   ljr2buff:  buffer of Lennard-Jones success displacements           
//   qIDbuff:   buffer of electrostatic success IDs into the cell list  
//   qr2buff:   buffer of electrostatic success displacements           
//-----------------------------------------------------------------------------
static int CullDualRange(int nljpair, double Ecut2, int* ljIDbuff,
                         rngbuff* ljr2buff, int* qIDbuff, rngbuff* qr2buff)
{
  int j, nqpair;

  nqpair = 0;
  for (j = 0; j < nljpair; j++) {
    if (ljr2buff[j].r2 < Ecut2) {
      qIDbuff[nqpair] = ljIDbuff[j];
      qr2buff[nqpair] = ljr2buff[j];
      nqpair++;
    }
  }

  return nqpair;
}

//-----------------------------------------------------------------------------
// CopyAtomSuccess: directly copy the Lennard-Jones success list in the event
//                  that the electrostatic cutoff is the same as the
//                  Lennard-Jones cutoff.  The results are stored in qIDbuff
//                  and qr2buff; ljIDbuff and ljr2buff are unchanged.
//
// Arguments:                                                           
//   nljpair:   the number of Lennard-Jones successes reported from     
//              CalcAtomR2 above                                        
//   ljIDbuff:  buffer of Lennard-Jones success IDs into the cell list  
//   ljr2buff:  buffer of Lennard-Jones success displacements           
//   qIDbuff:   buffer of electrostatic success IDs into the cell list  
//   qr2buff:   buffer of electrostatic success displacements           
//-----------------------------------------------------------------------------
static int CopyAtomSuccess(int nljpair, int* ljIDbuff, rngbuff* ljr2buff,
			   int* qIDbuff, rngbuff* qr2buff)
{
  int j;

  for (j = 0; j < nljpair; j++) {
    qIDbuff[j] = ljIDbuff[j];
    qr2buff[j] = ljr2buff[j];
  }

  return nljpair;
}

//-----------------------------------------------------------------------------
// CullCloseLJ: further cull the Lennard-Jones success list by pairs that are
//              in close contact and part of bonded exclusions.
//
// Arguments:
//   nljpair:   the number of Lennard-Jones successes reported from     
//              CalcAtomR2 above                                        
//   atmJ:      atoms in the cell available for interactions            
//   ljIDbuff:  buffer of Lennard-Jones success IDs into the cell list  
//   ljr2buff:  buffer of Lennard-Jones success displacements           
//   aiid:      the ID number of the first atom in the pair (index into 
//              the topology, not the cell atom list)                   
//   tp:        the topology                                            
//-----------------------------------------------------------------------------
static int CullCloseLJ(int nljpair, atomc* atmJ, int* ljIDbuff,
		       rngbuff* ljr2buff, int aiid, prmtop *tp)
{
  int j, nljpairWC;

  nljpairWC = 0;
  for (j = 0; j < nljpair; j++) {
    if (ljr2buff[j].r2 >= MINNB2 ||
        TestBondedExclusion(aiid, atmJ[ljIDbuff[j]].id, tp) == 0) {
      ljIDbuff[nljpairWC] = ljIDbuff[j];
      ljr2buff[nljpairWC] = ljr2buff[j];
      nljpairWC++;
    }
  }

  return nljpairWC;
}

// BLOCK 1: forces are required
#define NEEDFORCE 1

// BLOCK 1: energy (and virial) are not needed
#define NEEDENERGY 0
#define NEEDVIRIAL 0
#define CHKRANGE 1
#define LISTSAME 1
#include "NBOnlyLJ.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBOnlyLJ.c"
#undef LISTSAME
#undef CHKRANGE

#define CHKRANGE 0
#define LISTSAME 0
#include "NBOnlyLJ.c"
#undef LISTSAME
#undef CHKRANGE
#undef NEEDVIRIAL
#undef NEEDENERGY

#define NEEDENERGY 1
#define NEEDVIRIAL 0
#define CHKRANGE 1
#define LISTSAME 1
#include "NBOnlyLJ.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBOnlyLJ.c"
#undef LISTSAME
#undef CHKRANGE

#define CHKRANGE 0
#define LISTSAME 0
#include "NBOnlyLJ.c"
#undef LISTSAME
#undef CHKRANGE
#undef NEEDVIRIAL

#define NEEDVIRIAL 1
#define CHKRANGE 1
#define LISTSAME 1
#include "NBOnlyLJ.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBOnlyLJ.c"
#undef LISTSAME
#undef CHKRANGE

#define CHKRANGE 0
#define LISTSAME 0
#include "NBOnlyLJ.c"
#undef LISTSAME
#undef CHKRANGE
#undef NEEDVIRIAL
#undef NEEDENERGY
#undef NEEDFORCE

// BLOCK 2: the forces are not required
#define NEEDFORCE 0
#define NEEDENERGY 1 
#define NEEDVIRIAL 0
#define CHKRANGE 1
#define LISTSAME 1
#include "NBOnlyLJ.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBOnlyLJ.c"
#undef LISTSAME
#undef CHKRANGE

#define CHKRANGE 0
#define LISTSAME 0
#include "NBOnlyLJ.c"
#undef LISTSAME
#undef CHKRANGE
#undef NEEDVIRIAL
#undef NEEDENERGY
#undef NEEDFORCE

// BLOCK 1: the forces are required
#define NEEDFORCE 1
#define NEEDENERGY 0
#define NEEDVIRIAL 0
#define LISTSAME 1
#include "NBOnlyQQ.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBOnlyQQ.c"
#undef LISTSAME
#undef NEEDVIRIAL
#undef NEEDENERGY

#define NEEDENERGY 1
#define NEEDVIRIAL 0
#define LISTSAME 1
#include "NBOnlyQQ.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBOnlyQQ.c"
#undef LISTSAME
#undef NEEDVIRIAL

#define NEEDVIRIAL 1
#define LISTSAME 1
#include "NBOnlyQQ.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBOnlyQQ.c"
#undef LISTSAME
#undef NEEDVIRIAL
#undef NEEDENERGY
#undef NEEDFORCE

// BLOCK 2: the forces are not required
#define NEEDFORCE 0
#define NEEDENERGY 1
#define NEEDVIRIAL 0
#define LISTSAME 1
#include "NBOnlyQQ.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBOnlyQQ.c"
#undef LISTSAME
#undef NEEDENERGY
#undef NEEDVIRIAL
#undef NEEDFORCE

// BLOCKS 1 and 2: the forces are required
#define NEEDFORCE 1

// BLOCK 1: energy (and virial) are not needed
#define NEEDENERGY 0
#define NEEDVIRIAL 0
#define VGTE 1
#define CHKRANGE 1
#define LISTSAME 1
#include "NBLoop.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 1
#define CHKRANGE 0
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 0
#define CHKRANGE 1
#define LISTSAME 1
#include "NBLoop.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 0
#define CHKRANGE 0
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE
#undef NEEDVIRIAL
#undef NEEDENERGY
// END BLOCK 1

// BLOCK 2a: energy (but not virial) is needed,
//           in addition to forces             
#define NEEDENERGY 1
#define NEEDVIRIAL 0
#define VGTE 1
#define CHKRANGE 1
#define LISTSAME 1
#include "NBLoop.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 1
#define CHKRANGE 0
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 0
#define CHKRANGE 1
#define LISTSAME 1
#include "NBLoop.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 0
#define CHKRANGE 0
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE
#undef NEEDVIRIAL
// END BLOCK 2a

// BLOCK 2b: energy, virial, and forces are all needed
#define NEEDVIRIAL 1
#define VGTE 1
#define CHKRANGE 1
#define LISTSAME 1
#include "NBLoop.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 1
#define CHKRANGE 0
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 0
#define CHKRANGE 1
#define LISTSAME 1
#include "NBLoop.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 0
#define CHKRANGE 0
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE
#undef NEEDVIRIAL
#undef NEEDENERGY
// END BLOCK 2b

#undef NEEDFORCE
// END BLOCKS 1 and 2

// BLOCK 3: energy alone, but not forces or virials, are needed
#define NEEDFORCE 0
#define NEEDENERGY 1
#define NEEDVIRIAL 0
#define VGTE 1
#define CHKRANGE 1
#define LISTSAME 1
#include "NBLoop.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 1
#define CHKRANGE 0
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 0
#define CHKRANGE 1
#define LISTSAME 1
#include "NBLoop.c"
#undef LISTSAME
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE

#define VGTE 0
#define CHKRANGE 0
#define LISTSAME 0
#include "NBLoop.c"
#undef LISTSAME
#undef VGTE
#undef CHKRANGE
#undef NEEDVIRIAL
#undef NEEDENERGY
#undef NEEDFORCE
// END BLOCK 3
