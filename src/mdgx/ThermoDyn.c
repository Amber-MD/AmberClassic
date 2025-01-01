#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Macros.h"
#include "mdgxVector.h"
#include "VirtualSites.h"
#include "Topology.h"
#include "CellManip.h"
#include "Integrator.h"

#include "CrdManipDS.h"
#include "pmeRecipDS.h"
#include "pmeDirectDS.h"
#include "TimingsDS.h"

//-----------------------------------------------------------------------------
// CompareAtoms: compare two atoms in different topologies to determine whether
//               they are congruent.  Congruence is defined by the atoms having
//               identical positions and names.        
//                                                                      
// Arguments:                                                           
//   crdA     coordinates for system A                                  
//   crdB     coordinates for system B                                  
//   tpA:     topology for system A                                     
//   tpB:     topology for system B                                     
//   na:      atom from system A                                        
//   nb:      atom from system B                                        
//-----------------------------------------------------------------------------
static int CompareAtoms(coord *crdA, coord *crdB, prmtop *tpA, prmtop *tpB,
			int na, int nb)
{
  int resa, resb;
  double ax, ay, az, bx, by, bz;

  // First, compare the masses
  if (fabs(tpA->Masses[na] - tpB->Masses[nb]) > 1.0e-8) {
    return 0;
  }

  // Next, compare the locations
  ax = crdA->loc[3*na];
  ay = crdA->loc[3*na+1];
  az = crdA->loc[3*na+2];
  bx = crdB->loc[3*nb];
  by = crdB->loc[3*nb+1];
  bz = crdB->loc[3*nb+2];
  if (fabs(ax-bx) < 5.0e-8 && fabs(ay-by) < 5.0e-8 &&
      fabs(az-bz) < 5.0e-8) {

    // These atoms are in the same locations,
    // so do they have similar names?        
    if (strncmp(&tpA->AtomNames[4*na], &tpB->AtomNames[4*nb], 4) == 0) {

      // Finally, do these atoms have similar residue names?
      resa = LocateResID(tpA, na, 0, tpA->nres);
      resb = LocateResID(tpB, nb, 0, tpB->nres);
      if (strncmp(&tpA->ResNames[4*resa], &tpB->ResNames[4*resb], 4) == 0) {

	// Return true; these atoms correspond (although
	// that is not to say they precisely match)     
	return 1;
      }
    }
  }

  // If the two atoms did not match on       
  // at least one of the levels, return false
  return 0;
}

//-----------------------------------------------------------------------------
// MarkUniqueBondedTerm: this function will loop over bonds controlled by an
//                       atom and determine whether an equivalent bond is
//                       present in each of two topologies.  If the bond is
//                       unique to the first topology, it will be catalogged
//                       in the unique bonded terms list.
//                                                                      
// Arguments:                                                           
//   tp1:     topology for the original system                          
//   tp2:     topology for the second system (this function determines  
//            whether a bond is unique to tp1, not tp2--for the reverse 
//            situation, flip tp1 and tp2 in the arguments and call this    
//            function again)                                      
//   n1:      atom from the original system                             
//   n2:      corresponding atom from the second system                 
//   match1:  the topology correspondence map for atoms in system 1 matching
//            atoms in system 2                                
//   uterms:  the list of unique bonded terms                           
//   ordr:    order of the bonded term (2 = bonds, 3 = angles, 4 = dihedrals)
//-----------------------------------------------------------------------------
static void MarkUniqueBondedTerm(prmtop *tp1, prmtop *tp2, int n1, int n2,
				 int* match1, ublist *uterms, int ordr)
{
  int i, j, k, ilim, found, nt, eval14;
  int ratmA, ratmB, ratmC, ratmD;
  double keq, l0, th0, scee, scnb;
  bondlist *bond1, *bond2;
  angllist *angl1, *angl2;
  dihelist *dihe1, *dihe2;
  dihedef *ddef1, *ddef2;

  if (ordr == 2) {
    bond1 = &tp1->BLC[n1];
    bond2 = &tp2->BLC[n2];
    ilim = bond1->nbond;
  }
  else if (ordr == 3) {
    angl1 = &tp1->ALC[n1];
    angl2 = &tp2->ALC[n2];
    ilim = angl1->nangl;
  }
  else if (ordr == 4) {
    dihe1 = &tp1->HLC[n1];
    dihe2 = &tp2->HLC[n2];
    ilim = dihe1->ndihe;
  }
  for (i = 0; i < ilim; i++) {

    // First, make pointers to the bond / angle / dihedral list
    // in topology A.  Then, loop over all bonds controlled by 
    // this atom's counterpart in topology B.  If there is a   
    // match, then this bond has an equivalent in each topology
    // and does not need further consideration.                
    found = 0;
    if (ordr == 2) {
      ratmB = bond1->BC[i].b;
      l0 = tp1->BParam[bond1->BC[i].t].l0;
      keq = tp1->BParam[bond1->BC[i].t].K;
      for (j = 0; j < bond2->nbond; j++) {
	if (match1[ratmB] == bond2->BC[j].b &&
	    fabs(l0 - tp2->BParam[bond2->BC[j].t].l0) < 1.0e-8 &&
	    fabs(keq - tp2->BParam[bond2->BC[j].t].K) < 1.0e-8) {
	  found = 1;
	}
      }
    }
    else if (ordr == 3) {
      ratmA = angl1->AC[i].a;
      ratmC = angl1->AC[i].c;
      th0 = tp1->AParam[angl1->AC[i].t].th0;
      keq = tp1->AParam[angl1->AC[i].t].K;
      for (j = 0; j < angl2->nangl; j++) {
        if (match1[ratmA] == angl2->AC[j].a &&
	    match1[ratmC] == angl2->AC[j].c &&
            fabs(th0 - tp2->AParam[angl2->AC[j].t].th0) < 1.0e-8 &&
	    fabs(keq - tp2->AParam[angl2->AC[j].t].K) < 1.0e-8) {
          found = 1;
        }
      }
    }
    else if (ordr == 4) {
      ratmA = dihe1->HC[i].a;
      ratmC = dihe1->HC[i].c;
      ratmD = dihe1->HC[i].d;
      nt = dihe1->HC[i].nt;
      scnb = dihe1->HC[i].scnb;
      scee = dihe1->HC[i].scee;
      eval14 = dihe1->HC[i].eval14;
      for (j = 0; j < dihe2->ndihe; j++) {
        if (match1[ratmA] == dihe2->HC[j].a &&
            match1[ratmC] == dihe2->HC[j].c &&
            match1[ratmD] == dihe2->HC[j].d &&
	    nt == dihe2->HC[j].nt && eval14 == dihe2->HC[j].eval14 &&
	    fabs(scnb - dihe2->HC[j].scnb) < 1.0e-8 &&
	    fabs(scee - dihe2->HC[j].scee) < 1.0e-8) {
          found = 1;

	  // So far, so good, but there is
	  // a fourier series to compare. 
	  for (k = 0; k < nt; k++) {
	    ddef1 = &tp1->HParam[dihe1->HC[j].t[k]];
	    ddef2 = &tp2->HParam[dihe2->HC[j].t[k]];
	    if (fabs(ddef1->K - ddef2->K) > 1.0e-8 ||
		fabs(ddef1->N - ddef2->N) > 1.0e-8 ||
		fabs(ddef1->Phi - ddef2->Phi) > 1.0e-8) {
	      found = 0;
	    }
	  }
        }
      }
    }
    if (found == 0) {
      j = uterms->nitem;
      uterms->items = (int*)realloc(uterms->items, (j+1)*sizeof(int));
      uterms->items[j] = i;
      uterms->nitem = j+1;
    }
  }
}

//-----------------------------------------------------------------------------
// CrossRefAtoms: look more closely at two atoms which are known to be  
//                comparable: do they have similar bonding patterns, or 
//                identical nonbonded charateristics?  If so, note this 
//                information in the topology correlation struct.       
//                                                                      
// Arguments:                                                           
//   tpA:     topology for system A                                     
//   tpB:     topology for system B                                     
//   na:      atom from system A (matches atom nb from system B)        
//   nb:      atom from system B (matches atom na from system A)        
//   prc:     the topology correlation map                              
//-----------------------------------------------------------------------------
static void CrossRefAtoms(prmtop *tpA, prmtop *tpB, int na, int nb,
			  prmcorr *prc)
{
  // First, check nonbonded characteristics; if these match, then
  // it is likely that forces can be computed once for the first 
  // system and only slightly modified in the second system.     
  if (fabs(tpA->Charges[na] - tpB->Charges[nb]) < 1.0e-8) {
    prc->corrA[na] += 1;
    prc->corrB[nb] += 1;
    prc->dQ[nb] = 0.0;
  }
  else {

    // Store the difference in charge that would be necessary
    // in order to bring these topologies into agreement     
    prc->dQ[nb] = tpB->Charges[nb] - tpA->Charges[na];
  }

  // Next, check van-der Waals characteristics; there is an
  // implicit assumption that the same type names imply the
  // same type:type interactions in each topology.         
  if (strncmp(&tpA->AtomTypes[4*na], &tpB->AtomTypes[4*nb], 4) == 0) {
    prc->corrA[na] += 2;
    prc->corrB[nb] += 2;
  }

  // Now, check bonding patterns.  If a bond (or angle, or dihedral)
  // exists in topology A but does not exist in topology B, then its
  // effects will need to be removed from the forces and energies   
  // computed for topology A before they are added to the forces and
  // energies computed for topology B.  Bonded terms that change in 
  // stiffness or equilibrium length between the two topologies will
  // likewise need to be subtracted from A and added to B.          
  MarkUniqueBondedTerm(tpA, tpB, na, nb, prc->matchA, &prc->SubBondA[na], 2);
  MarkUniqueBondedTerm(tpB, tpA, nb, na, prc->matchB, &prc->AddBondB[nb], 2);
  MarkUniqueBondedTerm(tpA, tpB, na, nb, prc->matchA, &prc->SubAnglA[na], 3);
  MarkUniqueBondedTerm(tpB, tpA, nb, na, prc->matchB, &prc->AddAnglB[nb], 3);
  MarkUniqueBondedTerm(tpA, tpB, na, nb, prc->matchA, &prc->SubDiheA[na], 4);
  MarkUniqueBondedTerm(tpB, tpA, nb, na, prc->matchB, &prc->AddDiheB[nb], 4);
  if (prc->SubBondA[na].nitem == 0 && prc->AddBondB[nb].nitem == 0 &&
      prc->SubAnglA[na].nitem == 0 && prc->AddAnglB[nb].nitem == 0 &&
      prc->SubDiheA[na].nitem == 0 && prc->AddDiheB[nb].nitem == 0) {
    prc->corrA[na] += 4;
    prc->corrB[nb] += 4;
  }
}

//-----------------------------------------------------------------------------
// CompTpCorr: compute the correspondence between topologies, atom by atom.
//             Changes in charge, changes in Lennard-Jones atom type, changes
//             in bond structure, correspondences between atoms.
//
// Arguments:
//   tpA    : topology of system A
//   tpB    : topology of system B
//   crdA   : coordinates of system A
//   crdB   : coordinates of system B
//-----------------------------------------------------------------------------
prmcorr CompTpCorr(prmtop *tpA, prmtop *tpB, coord *crdA, coord *crdB)
{
  int i, j, lmnB, matchfound;
  prmcorr prc;

  // Store pointers to each topology that
  // this correspondence map will relate 
  prc.tpA = tpA;
  prc.tpB = tpB;

  // Initialize all atoms in each topology to have
  // no correspondence in the other topology      
  prc.matchA = (int*)malloc(tpA->natom*sizeof(int));
  prc.matchB = (int*)malloc(tpB->natom*sizeof(int));
  SetIVec(prc.matchA, tpA->natom, -1);
  SetIVec(prc.matchB, tpB->natom, -1);
  prc.corrA = (int*)calloc(tpA->natom, sizeof(int));
  prc.corrB = (int*)calloc(tpB->natom, sizeof(int));
  prc.dQ = (double*)calloc(tpB->natom, sizeof(double));
  prc.SubBondA = (ublist*)calloc(tpA->natom, sizeof(ublist));
  prc.AddBondB = (ublist*)calloc(tpB->natom, sizeof(ublist));
  prc.SubAnglA = (ublist*)calloc(tpA->natom, sizeof(ublist));
  prc.AddAnglB = (ublist*)calloc(tpB->natom, sizeof(ublist));
  prc.SubDiheA = (ublist*)calloc(tpA->natom, sizeof(ublist));
  prc.AddDiheB = (ublist*)calloc(tpB->natom, sizeof(ublist));

  // Loop over all atoms, find counterparts
  // with the same coordinates and names   
  lmnB = 0;
  for (i = 0; i < tpA->natom; i++) {

    // Start the search for a corresponding atom
    // from the last matching atom found in the 
    // second topology / coordinate set         
    matchfound = 0;
    for (j = lmnB; j < tpB->natom; j++) {
      if (CompareAtoms(crdA, crdB, tpA, tpB, i, j) == 1) {
	prc.matchA[i] = j;
	prc.matchB[j] = i;
	matchfound = 1;
	lmnB = j;
	break;
      }
    }

    // Continue on to the next atom if
    // a match has already been found 
    if (matchfound == 1) {
      continue;
    }

    // Wrap around and search the second topology /
    // coordinates set from the beginning.         
    for (j = 0; j < lmnB; j++) {
      if (CompareAtoms(crdA, crdB, tpA, tpB, i, j) == 1) {
	prc.matchA[i] = j;
	prc.matchB[j] = i;
	break;
      }
    }
  }

  // Thus far, we have created a list of atom correspondences
  // for each topology.  If these topologies contain the same
  // atoms, it may make things much easier.  Even if they do 
  // not, it is necessary to understand how the topologies   
  // relate to each other.                                   
  for (i = 0; i < tpA->natom; i++) {
    if (prc.matchA[i] >= 0) {
      CrossRefAtoms(tpA, tpB, i, prc.matchA[i], &prc);
    }
  }

  // Count the number of unique atoms in systems A and
  // B.  Set the "relate" field based on the totality 
  // of correspondence between systems A and B.       
  prc.uniA = 0;
  for (i = 0; i < tpA->natom; i++) {
    if (prc.matchA[i] < 0) {
      prc.uniA += 1; 
    }
  }
  prc.comAB = tpA->natom - prc.uniA;
  prc.uniB = 0;
  for (i = 0; i < tpB->natom; i++) {
    if (prc.matchB[i] < 0) {
      prc.uniB += 1; 
    }
  }
  prc.relate =  (prc.uniA > 0 || prc.uniB > 0) ? 1 : 0;
  if (prc.relate > 0) {
    matchfound = -1;
    for (i = 0; i < tpA->natom; i++) {
      if (prc.matchA[i] >= 0) {
	if (prc.matchA[i] < matchfound) {
	  prc.relate = 2;
	  break;
	}
	matchfound = prc.matchA[i];
      }
    }
  }

  // Determine the number of mutual degrees of freedom
  prc.nxdf = tpA->ndf;
  prc.nxprtcl = tpA->nprtcl;
  for (i = 0; i < tpB->natom; i++) {
    if (prc.matchB[i] < 0 && tpB->Masses[i] > 1.0e-8) {
      prc.nxdf += 3;
      prc.nxprtcl += 1;
      if (tpB->SHL[i].exe == 1) {
	prc.nxdf -= 3;
	prc.nxprtcl -= 2;
      }
    }
  }

  return prc;
}

//-----------------------------------------------------------------------------
// ForceReeval: re-evaluate forces, energies, and virials with a second 
//              topology for the purpose of computing a hybrid trajectory in
//              which the coordinates are propagated by a weighted average of
//              the forces due to each of the two topologies.  This propagation
//              of a single trajectory by a mixture of two Hamiltonians is a
//              major part of certain thermodynamic integration protocols.
//                                                                      
// Arguments:                                                           
//   crd      : system coordinates
//   CG       : the cell grid for this system, decomposing all direct space
//              interactions
//   tp       : system topology
//   dcinp    : dierct space input control data
//   Etab     : (low-density) electrostatic direct space energy lookup table
//   EHtab    : high-density electrostatic direct space energy lookup table
//   rcinp    : reciprocal space input control data
//   PPk      : reciprocal space mesh-to-mesh pair potential tool kit
//   sysUV    : system energy and virial accumulators
//   etimers  : timers and time accumulators for each part of the calculation
//   prc      : topology correlation structure
//   tj       : trajectory input control data
//-----------------------------------------------------------------------------
void ForceReeval(coord* crd, cellgrid* CG, prmtop* tp, dircon* dcinp,
		 FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit* PPk,
		 Energy* sysUV, execon *etimers, prmcorr *prc, trajcon *tj)
{

  // Topologies A and B contain many of the same
  // atoms, but may contain different parameters
  if (prc->vdw == 0 && prc->elec == 0) {

    // The non-bonded forces from topologies A and B are
    // consistent, so forces can be copied over directly

    // 1:4 interactions may have changed,
    // though, so account for that       
    if (prc->vdw14 == 1 || prc->elec14 == 1) {

    }
  }
  else {

    // Non-bonded interactions have changed; the    
    // easiest thing is just to recompute everything
    AtomForces(&crd[1], &CG[1], &tp[1], dcinp, Etab, EHtab, rcinp, &PPk[1],
	       &sysUV[1], etimers, tj);
  }

  // Add the interactions of appearing atoms in topology B 
  // with one another (by standard non-bonded potentials)  
  // and with atoms already existing in topology A (by     
  // modified non-bonded potentials)  The general strategy 
  // is to replace r in 1/(r^n) expressions with a lambda- 
  // dependent function 1/((b*L + r)^n) where b is a simple
  // function (or even a constant) and L is lambda or      
  // 1-lambda if atoms are disappearing or appearing,      
  // respectively.  The modified potentials do not go to   
  // singularities even as r approaches zero.              
  if (prc->appear == 1) {

  }
  if (prc->disappear == 1) {

  }
}
