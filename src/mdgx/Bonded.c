#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mdgxVector.h"
#include "Bonded.h"
#include "CrdManip.h"
#include "CellManip.h"
#include "Macros.h"
#include "Debug.h"

#include "TopologyDS.h"
#include "TrajectoryDS.h"

//-----------------------------------------------------------------------------
// ElimInteraction: atoms X and Y may interact, but they should not contribute
//                  to the total energy and forces because of exclusions. 
//                  This routine will test whether the interaction of atoms X
//                  and Y has contributed to energies, forces, or virials,
//                  then determine what needs to be done to remove that
//                  contribution.
//
// Arguments:                                                           
//   C:      the cell of interest
//   elist:  the list of eliminations to perform                        
//   nelim:  the length of elist                                        
//   tp:     the topology                                               
//   iordr:  the order of the bonded exclusion, interactions are 1:iordr
//   Cfrc:   the "coarse" electrostatic force/energy spline table       
//   Hfrc:   the "fine" electrostatic force/energy spline table         
//   sysUV:  the system energy and virial (state information)           
//-----------------------------------------------------------------------------
static void ElimInteraction(cell *C, nixpr* elist, int nelim, prmtop *tp,
			    int iordr, FrcTab *Cfrc, FrcTab *Hfrc,
			    Energy *sysUV, int qform)
{
  int i, xid, yid;
  double vdwscl, elecscl;
  double *xloc, *yloc;

  if (iordr == 4) {
    vdwscl = tp->lj14fac;
    elecscl = tp->elec14fac;
  }
  else {
    vdwscl = 1.0;
    elecscl = 1.0;
  }
  for (i = 0; i < nelim; i++) {
    xid = C->GPSptr[elist[i].atmX];
    yid = C->GPSptr[elist[i].atmY];
    xloc = C->data[xid].loc;
    yloc = C->data[yid].loc;
    if (sysUV->updateU == 2) {
      AttenuatePairVir(tp, Cfrc, Hfrc, elist[i].atmX, elist[i].atmY,
		       yloc[0]-xloc[0], yloc[1]-xloc[1], yloc[2]-xloc[2], 0.0,
		       C->data[xid].frc, C->data[yid].frc, sysUV,
		       elecscl, vdwscl, qform);
    }
    else if (sysUV->updateU == 1) {
      AttenuatePairFrcNrg(tp, Cfrc, Hfrc, elist[i].atmX, elist[i].atmY,
		          yloc[0]-xloc[0], yloc[1]-xloc[1], yloc[2]-xloc[2],
			  0.0, C->data[xid].frc, C->data[yid].frc, sysUV,
			  elecscl, vdwscl, qform);
    }
    else if (sysUV->updateU == 0) {
      AttenuatePairFrc(tp, Cfrc, Hfrc, elist[i].atmX, elist[i].atmY,
		       yloc[0]-xloc[0], yloc[1]-xloc[1], yloc[2]-xloc[2],
		       0.0, C->data[xid].frc, C->data[yid].frc, elecscl,
		       vdwscl, qform);
    }
    else if (sysUV->updateU == -1) {
      AttenuatePairNrg(tp, Cfrc, Hfrc, elist[i].atmX, elist[i].atmY,
		       yloc[0]-xloc[0], yloc[1]-xloc[1], yloc[2]-xloc[2],
		       sysUV, elecscl, vdwscl, qform);
    }
  }
}

//-----------------------------------------------------------------------------
// CellBondedIntr: compute bonded interactions for a cell decomposition of atom
//                 coordinates.  This routine can be called only after the
//                 "ShareTriple2" function in the CellManip library has shared
//                 atoms between cells.  This routine searches for atoms within
//                 a region defined as the same size as the home cell C, but
//                 translated by one half the nonbonded cutoff in all
//                 directions.         
//
// Arguments:                                                           
//   C:      the cell of interest                                       
//   CG:     the cell grid (used for dimensions)                        
//   tp:     the topology                                               
//   crd:    the coordinates (used for transformation matrices)         
//   Cfrc:   the "coarse" electrostatic force/energy spline table       
//   Hfrc:   the "fine" electrostatic force/energy spline table         
//   sysUV:  the system energy and virial (state information)           
//   qform:  format of the charges (points, or spherical Gaussians)     
//-----------------------------------------------------------------------------
void CellBondedIntr(cell *C, cellgrid *CG, coord *crd, prmtop *tp,
		    FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV, int qform)
{
  int i, j, k, atmid, aid, bid, cid, did, eid;
  double *atmcrd, *atmfrc;
  bondlist *tblc;
  bondcomm *tbcom;
  angllist *talc;
  anglcomm *tacom;
  dihelist *thlc;
  dihecomm *thcom;
  cmaplist *tmlc;
  cmapcomm *tmcom;
  auxelim *eptmp;

#ifdef MPI
  // Set counters for load-balancing
  C->nbond = 0;
  C->nangl = 0;
  C->ndihe = 0;
  C->ncmap = 0;
  C->nelim = 0;
#endif

  // Loop over all regions of this cell
  for (i = 0; i < 8; i++) {

    // Loop over all atoms within the cell region
    for (j = 0; j < C->nr[i]; j++) {
      if (IsCentralAtom(C->orig, &CG->celldim[4], C->map[i][j].loc,
			crd->U.data, CG->dbng, crd->isortho) == -1) {
	continue;
      }

      // If we are still here, then this atom may control
      // bonded interactions that should be computed.    
      atmid = C->map[i][j].id;
      atmcrd = C->map[i][j].loc;
      atmfrc = C->map[i][j].frc;

      // Bonds
      if (tp->BLC[atmid].nbond > 0) {
	tblc = &tp->BLC[atmid];
#ifdef MPI
        C->nbond += tblc->nbond;
#endif
	for (k = 0; k < tblc->nbond; k++) {
	  tbcom = &tblc->BC[k];
	  bid = C->GPSptr[tbcom->b];
	  if (sysUV->updateU == 2) {
	    BondVir(atmcrd, C->data[bid].loc, atmfrc, C->data[bid].frc,
		    tbcom, tp, Cfrc, Hfrc, sysUV, qform);
	  }
	  else if (sysUV->updateU == 1) {
	    BondFrcNrg(atmcrd, C->data[bid].loc, atmfrc,
		       C->data[bid].frc, tbcom, tp, Cfrc, Hfrc, sysUV, qform);
	  }
	  else if (sysUV->updateU == -1) {
	    BondNrg(atmcrd, C->data[bid].loc, tbcom, tp, Cfrc, Hfrc,
		    sysUV, qform);
	  }
	  else {
	    BondFrc(atmcrd, C->data[bid].loc, atmfrc, C->data[bid].frc,
		    tbcom, tp, Cfrc, Hfrc, qform);
	  }
	}
      }

      // Angles
      if (tp->ALC[atmid].nangl > 0) {
	talc = &tp->ALC[atmid];
#ifdef MPI
        C->nangl += talc->nangl;
#endif
	for (k = 0; k < talc->nangl; k++) {
	  tacom = &talc->AC[k];
	  aid = C->GPSptr[tacom->a];
	  cid = C->GPSptr[tacom->c];
	  if (sysUV->updateU == 2) {
	    AnglVir(C->data[aid].loc, atmcrd, C->data[cid].loc,
		    C->data[aid].frc, atmfrc, C->data[cid].frc,
		    tacom, tp, Cfrc, Hfrc, sysUV, qform);
	  }
	  else if (sysUV->updateU == 1) {
	    AnglFrcNrg(C->data[aid].loc, atmcrd, C->data[cid].loc,
		       C->data[aid].frc, atmfrc, C->data[cid].frc, tacom, tp,
		       Cfrc, Hfrc, sysUV, qform);
	  }
	  else if (sysUV->updateU == -1) {
	    AnglNrg(C->data[aid].loc, atmcrd, C->data[cid].loc,
		    tacom, tp, Cfrc, Hfrc, sysUV, qform);
	  }
	  else {
	    AnglFrc(C->data[aid].loc, atmcrd, C->data[cid].loc,
		    C->data[aid].frc, atmfrc, C->data[cid].frc,
		    tacom, tp, Cfrc, Hfrc, qform);
	  }
	}
      }

      // Dihedrals
      if (tp->HLC[atmid].ndihe > 0) {
	thlc = &tp->HLC[atmid];
#ifdef MPI
        C->ndihe += thlc->ndihe;
#endif
	for (k = 0; k < thlc->ndihe; k++) {
	  thcom = &thlc->HC[k];
	  aid = C->GPSptr[thcom->a];
	  cid = C->GPSptr[thcom->c];
	  did = C->GPSptr[thcom->d];
	  if (sysUV->updateU == 2) {
	    DiheVir(C->data[aid].loc, atmcrd, C->data[cid].loc,
		    C->data[did].loc, C->data[aid].frc, atmfrc,
		    C->data[cid].frc, C->data[did].frc, thcom,
		    tp, Cfrc, Hfrc, sysUV, qform);
	  }
	  else if (sysUV->updateU == 1) {
	    DiheFrcNrg(C->data[aid].loc, atmcrd, C->data[cid].loc,
		       C->data[did].loc, C->data[aid].frc, atmfrc,
		       C->data[cid].frc, C->data[did].frc, thcom,
		       tp, Cfrc, Hfrc, sysUV, qform);
	  }
	  else if (sysUV->updateU == -1) {
	    DiheNrg(C->data[aid].loc, atmcrd, C->data[cid].loc,
		    C->data[did].loc, thcom, tp, Cfrc, Hfrc, sysUV, qform);
	  }
	  else if (sysUV->updateU == 0) {
	    DiheFrc(C->data[aid].loc, atmcrd, C->data[cid].loc,
		    C->data[did].loc, C->data[aid].frc, atmfrc,
		    C->data[cid].frc, C->data[did].frc, thcom,
		    tp, Cfrc, Hfrc, qform);
	  }
	}
      }
      
      // CMAPs
      if (tp->MLC[atmid].ncmap > 0) {
	tmlc = &tp->MLC[atmid];
#ifdef MPI
        C->ncmap += tmlc->ncmap;
#endif
	for (k = 0; k < tmlc->ncmap; k++) {
	  tmcom = &tmlc->MC[k];
	  aid = C->GPSptr[tmcom->a];
	  bid = C->GPSptr[tmcom->b];
	  did = C->GPSptr[tmcom->d];
	  eid = C->GPSptr[tmcom->e];
	  if (sysUV->updateU == 2) {
	    CmapVir(C->data[aid].loc, C->data[bid].loc, atmcrd,
		    C->data[did].loc, C->data[eid].loc, C->data[aid].frc,
		    C->data[bid].frc, atmfrc, C->data[did].frc,
		    C->data[eid].frc, tmcom, tp, sysUV, qform);
	  }
	  else if (sysUV->updateU == 1) {
	    CmapFrcNrg(C->data[aid].loc, C->data[bid].loc, atmcrd,
		       C->data[did].loc, C->data[eid].loc, C->data[aid].frc,
		       C->data[bid].frc, atmfrc, C->data[did].frc,
		       C->data[eid].frc, tmcom, tp, sysUV, qform);
	  }
	  else if (sysUV->updateU == -1) {
	    CmapNrg(C->data[aid].loc, C->data[bid].loc, atmcrd,
		    C->data[did].loc, C->data[eid].loc, tmcom, tp,
		    sysUV, qform);
	  }
	  else if (sysUV->updateU == 0) {
	    CmapFrc(C->data[aid].loc, C->data[bid].loc, atmcrd,
		    C->data[did].loc, C->data[eid].loc, C->data[aid].frc,
		    C->data[bid].frc, atmfrc, C->data[did].frc,
		    C->data[eid].frc, tmcom, tp, qform);
	  }
	}
      }

      // If there have been no extra points added
      // to the topology, we can just continue
      if (tp->nxtrapt == 0) {
	continue;
      }
      eptmp = &tp->ElimPair[atmid];
      ElimInteraction(C, eptmp->list11, eptmp->n11, tp, 1, Cfrc, Hfrc, sysUV,
		      qform);
      ElimInteraction(C, eptmp->list12, eptmp->n12, tp, 2, Cfrc, Hfrc, sysUV,
		      qform);
      ElimInteraction(C, eptmp->list13, eptmp->n13, tp, 3, Cfrc, Hfrc, sysUV,
		      qform);
      ElimInteraction(C, eptmp->list14, eptmp->n14, tp, 4, Cfrc, Hfrc, sysUV,
		      qform);
#ifdef MPI
      C->nelim += eptmp->n11 + eptmp->n12 + eptmp->n13 + eptmp->n14;
#endif
    }
  }
}

// BLOCKS 1 and 2: the forces are required
#define NEEDFORCE 1

// BLOCK 1: energy (and virial) are not needed
#define NEEDENERGY 0
#define NEEDVIRIAL 0
#include "BondContrib.c"
#undef NEEDVIRIAL
#undef NEEDENERGY

// BLOCK 2a: energy (but not virial) is needed,
//           in addition to forces             
#define NEEDENERGY 1
#define NEEDVIRIAL 0
#include "BondContrib.c"
#undef NEEDVIRIAL

// BLOCK 2b: energy, virial, and forces are all needed
#define NEEDVIRIAL 1
#include "BondContrib.c"
#undef NEEDVIRIAL
#undef NEEDENERGY

#undef NEEDFORCE
// END BLOCKS 1 and 2

// BLOCK 3: energy alone, but not forces or virials, are needed
#define NEEDFORCE 0
#define NEEDENERGY 1
#define NEEDVIRIAL 0
#include "BondContrib.c"
#undef NEEDVIRIAL
#undef NEEDENERGY
#undef NEEDFORCE
// END BLOCK 3
