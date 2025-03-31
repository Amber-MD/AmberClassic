#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mdgxVector.h"
#include "Matrix.h"
#include "ParamFit.h"
#include "Topology.h"
#include "CrdManip.h"
#include "Parse.h"
#include "ChargeFit.h"
#include "VirtualSites.h"
#include "Manual.h"
#include "Nonbonded.h"
#include "Trajectory.h"
#include "Macros.h"
#include "ParamOut.h"
#include "ParamRead.h"

#include "CrdManipDS.h"

//-----------------------------------------------------------------------------
// str4cmp: a faster 4-character string comparison function, because the
//          TypeCompare function is REALLY slow when order == 4.  Returns 0 if
//          the two strings are identical, 1 if not.       
//
// Arguments:
//   A, B:  the strings to compare, assumed to be no more and no less than
//          four characters each
//-----------------------------------------------------------------------------
int str4cmp(char* A, char* B)
{
  if (A[0] == B[0] && A[1] == B[1] && A[2] == B[2] && A[3] == B[3]) {
    return 0;
  }
  else {
    return 1;
  }
}

//-----------------------------------------------------------------------------
// strendswith: Checks if one string ends with another.  Returns 1 if true, 0
//              otherwise (opposite of how strcmp, strncmp work).
//
// Arguments:
//   A, B:  the strings to compare
//-----------------------------------------------------------------------------
int strendswith(char* A, char* B)
{
  int A_len, B_len;

  if (A == NULL || B == NULL) {
    return 0;
  }
  A_len = strlen(A);
  B_len = strlen(B);

  if (B_len > A_len) {
    return 0;
  }
  else {
    return (0 == strncmp(&A[A_len-B_len], B, B_len));
  }
}

//-----------------------------------------------------------------------------
// TypeCompare: compare the types of two pairs of atoms, returning 2 if the
//              types can be connected directly, 1 if they can be connected
//              by a wildcard, and 0 if they cannot connect.  
//                                                                      
// Arguments:                                                           
//   T[1,2][a-d]:   the first and second sets' various types a-d (note that the
//                  second set should come from the parameter library file,
//                  i.e. parm99.dat)            
//   order:         order of the comparison (2 = bond, 3 = angle, 4 = torsion)
//   strict:        if set to 1, this will require that the atom type names
//                  strictly match in either forward or reverse order; no
//                  wildcards!                                
//-----------------------------------------------------------------------------
int TypeCompare(char* T1a, char* T1b, char* T1c, char* T1d, char* T2a,
                char* T2b, char* T2c, char* T2d, int order, int strict)
{
  int aworks, bworks, cworks, dworks;
  int awild, bwild, cwild, dwild, nwild, ataken, btaken, ctaken, dtaken;

  if (order == 2) {
    if ((str4cmp(T1a, T2a) == 0 && str4cmp(T1b, T2b) == 0) ||
        (str4cmp(T1a, T2b) == 0 && str4cmp(T1b, T2a) == 0)) {
      return 2;
    }
    if (strict == 1) {
      return 0;
    }
    if ((str4cmp(T1a, T2a) == 0 || str4cmp(T2a, "X   ") == 0) && 
        (str4cmp(T1b, T2b) == 0 || str4cmp(T2b, "X   ") == 0)) {
      return 1;
    }
    if ((str4cmp(T1a, T2b) == 0 || str4cmp(T2b, "X   ") == 0) && 
        (str4cmp(T1b, T2a) == 0 || str4cmp(T2a, "X   ") == 0)) {
      return 1;
    }
    if (str4cmp(T2a, "X   ") == 0 && str4cmp(T2b, "X   ") == 0) {
      return 1;
    }
  }
  else if (order == 3) {
    if ((str4cmp(T1a, T2a) == 0 && str4cmp(T1b, T2b) == 0 &&
         str4cmp(T1c, T2c) == 0) ||
        (str4cmp(T1a, T2c) == 0 && str4cmp(T1b, T2b) == 0 &&
         str4cmp(T1c, T2a) == 0)) {
      return 2;
    }
    if (strict == 1) {
      return 0;
    }

    // Check for parallel correspondence, with wildcards
    aworks = (str4cmp(T2a, "X   ") == 0 || str4cmp(T1a, T2a) == 0) ? 1 : 0;
    bworks = (str4cmp(T2b, "X   ") == 0 || str4cmp(T1b, T2b) == 0) ? 1 : 0;
    cworks = (str4cmp(T2c, "X   ") == 0 || str4cmp(T1c, T2c) == 0) ? 1 : 0;
    if (aworks == 1 && bworks == 1 && cworks == 1) {
      return 1;
    }

    // Check for anti-parallel correspondence, with wildcards
    cworks = (str4cmp(T2a, "X   ") == 0 || str4cmp(T1c, T2a) == 0) ? 1 : 0;
    bworks = (str4cmp(T2b, "X   ") == 0 || str4cmp(T1b, T2b) == 0) ? 1 : 0;
    aworks = (str4cmp(T2c, "X   ") == 0 || str4cmp(T1a, T2c) == 0) ? 1 : 0;
    if (aworks == 1 && bworks == 1 && cworks == 1) {
      return 1;
    }
  }
  else if (order == 4 || order == 5) {

    // Check for explicit correspondence in the types
    if ((str4cmp(T1a, T2a) == 0 && str4cmp(T1b, T2b) == 0 &&
         str4cmp(T1c, T2c) == 0 && str4cmp(T1d, T2d) == 0) ||
        (str4cmp(T1a, T2d) == 0 && str4cmp(T1b, T2c) == 0 &&
         str4cmp(T1c, T2b) == 0 && str4cmp(T1d, T2a) == 0)) {
      return 2;
    }
    if (strict == 1) {
      return 0;
    }

    // Check for parallel correspondence, with wildcards
    aworks = (str4cmp(T2a, "X   ") == 0 || str4cmp(T1a, T2a) == 0) ? 1 : 0;
    bworks = (str4cmp(T2b, "X   ") == 0 || str4cmp(T1b, T2b) == 0) ? 1 : 0;
    cworks = (str4cmp(T2c, "X   ") == 0 || str4cmp(T1c, T2c) == 0) ? 1 : 0;
    dworks = (str4cmp(T2d, "X   ") == 0 || str4cmp(T1d, T2d) == 0) ? 1 : 0;
    if (aworks == 1 && bworks == 1 && cworks == 1 && dworks == 1) {
      return 1;
    }

    // Check for anti-parallel correspondence, with wildcards
    dworks = (str4cmp(T2a, "X   ") == 0 || str4cmp(T1d, T2a) == 0) ? 1 : 0;
    cworks = (str4cmp(T2b, "X   ") == 0 || str4cmp(T1c, T2b) == 0) ? 1 : 0;
    bworks = (str4cmp(T2c, "X   ") == 0 || str4cmp(T1b, T2c) == 0) ? 1 : 0;
    aworks = (str4cmp(T2d, "X   ") == 0 || str4cmp(T1a, T2d) == 0) ? 1 : 0;
    if (aworks == 1 && bworks == 1 && cworks == 1 && dworks == 1) {
      return 1;
    }

    // Last chance to match an improper dihedral with wildcards
    if (order == 5) {

      // The third atom type of the second set must be  
      // the center, but it could still be a wildcard.  
      // The goal is now to match the non-wildcard atoms
      // of the second set (the template type sequence) 
      // uniquely to atom types in the first set.       
      awild = (str4cmp(T2a, "X   ") == 0) ? 1 : 0;
      bwild = (str4cmp(T2b, "X   ") == 0) ? 1 : 0;
      cwild = (str4cmp(T2c, "X   ") == 0) ? 1 : 0;
      dwild = (str4cmp(T2d, "X   ") == 0) ? 1 : 0;
      nwild = awild + bwild + dwild;
      if (nwild == 0) {
        return 0;
      }

      // First assume that the third atom    
      // type of the first set is the center.
      if (cwild || str4cmp(T1c, T2c) == 0) {
        aworks = 0;
        bworks = 0;
        dworks = 0;
        ataken = 0;
        btaken = 0;
        dtaken = 0;
        if (awild == 0) {
          if (str4cmp(T1a, T2a) == 0) {
            aworks = 1;
            ataken = 1;
          }
          else if (str4cmp(T1b, T2a) == 0) {
            aworks = 1;
            btaken = 1;
          }
          else if (str4cmp(T1d, T2a) == 0) {
            aworks = 1;
            dtaken = 1;
          }
        }
        else {
          aworks = 1;
        }
        if (bwild == 0) {
          if (ataken == 0 && str4cmp(T1a, T2b) == 0) {
            bworks = 1;
            ataken = 1;
          }
          else if (btaken == 0 && str4cmp(T1b, T2b) == 0) {
            bworks = 1;
            btaken = 1;
          }
          else if (dtaken == 0 && str4cmp(T1d, T2b) == 0) {
            bworks = 1;
            dtaken = 1;
          }
        }
        else {
          bworks = 1;
        }
        if (dwild == 0) {
          if (ataken == 0 && str4cmp(T1a, T2d) == 0) {
            dworks = 1;
            ataken = 1;
          }
          else if (btaken == 0 && str4cmp(T1b, T2d) == 0) {
            dworks = 1;
            btaken = 1;
          }
          else if (dtaken == 0 && str4cmp(T1d, T2d) == 0) {
            dworks = 1;
            dtaken = 1;
          }
        }
        else {
          dworks = 1;
        }

        // It was implict that C works, by matching
        // the third atom of the first sequence.   
        if (aworks == 1 && bworks == 1 && dworks == 1) {
          return 1;
        }
      }

      // If that didn't work, assume that the second
      // type of the first set is the center.       
      if (cwild || str4cmp(T1b, T2c) == 0) {
        aworks = 0;
        bworks = 0;
        dworks = 0;
        ataken = 0;
        ctaken = 0;
        dtaken = 0;
        if (awild == 0) {
          if (str4cmp(T1a, T2a) == 0) {
            aworks = 1;
            ataken = 1;
          }
          else if (str4cmp(T1c, T2a) == 0) {
            aworks = 1;
            ctaken = 1;
          }
          else if (str4cmp(T1d, T2a) == 0) {
            aworks = 1;
            dtaken = 1;
          }
        }
        else {
          aworks = 1;
        }
        if (bwild == 0) {
          if (ataken == 0 && str4cmp(T1a, T2b) == 0) {
            bworks = 1;
            ataken = 1;
          }
          else if (ctaken == 0 && str4cmp(T1c, T2b) == 0) {
            bworks = 1;
            ctaken = 1;
          }
          else if (dtaken == 0 && str4cmp(T1d, T2b) == 0) {
            bworks = 1;
            dtaken = 1;
          }
        }
        else {
          bworks = 1;
        }
        if (dwild == 0) {
          if (ataken == 0 && str4cmp(T1a, T2d) == 0) {
            dworks = 1;
            ataken = 1;
          }
          else if (ctaken == 0 && str4cmp(T1c, T2d) == 0) {
            dworks = 1;
            ctaken = 1;
          }
          else if (dtaken == 0 && str4cmp(T1d, T2d) == 0) {
            dworks = 1;
            dtaken = 1;
          }
        }
        else {
          dworks = 1;
        }

        // It was implict that C works, by matching
        // the second atom of the first sequence.  
        if (aworks == 1 && bworks == 1 && dworks == 1) {
          return 1;
        }        
      }

      // If neither of those cases worked, this match is lost.
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// CmapResidueOverlap: function to determine whether two CMAPs apply to some
//                     or all of the same residues, and return a coded value
//                     indicating the result.
//
// Arguments:
//   
//-----------------------------------------------------------------------------
static int CmapResidueOverlap(xcmapdef *cmapA, xcmapdef *cmapB)
{
  int i, j, coverage;

  coverage = 0;
  for (i = 0; i < cmapA->nres; i++) {
    for (j = 0; j < cmapB->nres; j++) {
      if (str4cmp(cmapA->resnames.map[i], cmapB->resnames.map[j]) == 0) {
        coverage++;
      }
    }
  }
  
  return (coverage == cmapA->nres || coverage == cmapB->nres);
}

//-----------------------------------------------------------------------------
// DestroyXCmapDef: free data associated with an extended CMAP definition.
//                  It's not obvious where in the lbirary to place this
//                  function, but this will need to be done in at least two
//                  parts of the code, so an encapsulating function is
//                  appropriate.
//
// Arguments:
//   xcm:     the struct to destroy
//-----------------------------------------------------------------------------
static void DestroyXCmapDef(xcmapdef *xcm)
{
  free(xcm->title);
  DestroyDmat(&xcm->esrf);
  DestroyDmat(&xcm->dphi);
  DestroyDmat(&xcm->dpsi);
  DestroyDmat(&xcm->d2pp);
  DestroyDmat(&xcm->esrfFit);
  DestroyDmat(&xcm->dphiFit);
  DestroyDmat(&xcm->dpsiFit);
  DestroyDmat(&xcm->d2ppFit);
  DestroyDmat(&xcm->stencil);
  if (xcm->rstw > 0.0) {
    DestroyDmat(&xcm->target);
  }
  DestroyCmat(&xcm->resnames);
}

//-----------------------------------------------------------------------------
// PruneDuplicates: prune any duplicate definitions that can be found in a
//                  list.                                             
//                                                                      
// Arguments:                                                           
//   defs:         the list of definitions                              
//   count:        the initial number of definitions                    
//   abscount:     the absolute number of definitions in a longer list (for
//                 dihedrals, when impropers are a separate section of the
//                 same list)                                    
//   order:        the order of the bonded terms (2 = bonds, 3 = angles,
//                 4 = dihedrals, 5 = impropers, 6 = NMR operations, 7 = CMAPs)
//-----------------------------------------------------------------------------
static int PruneDuplicates(void* defs, int count, int abscount, int order)
{
  int i, j, k, m, nc, match, nappl;
  int* duplicate;
  xbonddef *bonds;
  xangldef *angls;
  torterm *torsions;
  xcmapdef *cmaps;
  
  // Return immediately if there are one or zero items
  if (count <= 1) {
    return count;
  }

  // Recast the list
  if (order == 2) {
    bonds = (xbonddef*)defs;
  }
  else if (order == 3) {
    angls = (xangldef*)defs;
  }
  else if (order == 4 || order == 5) {
    torsions = (torterm*)defs;
  }
  else if (order == 7) {
    cmaps = (xcmapdef*)defs;
  }
  
  // Prepare to mark duplicates
  duplicate = (int*)calloc(count, sizeof(int));
  for (i = 0; i < count-1; i++) {
    if (duplicate[i] > 0) {
      continue;
    }
    for (j = i+1; j < count; j++) {
      if (duplicate[j] > 0) {
        continue;
      }
      if (order == 2) {
        duplicate[j] = TypeCompare(bonds[i].atype, bonds[i].btype, "    ",
                                   "    ", bonds[j].atype, bonds[j].btype,
                                   "    ", "    ", order, 1);
      }
      else if (order == 3) {
        duplicate[j] = TypeCompare(angls[i].atype, angls[i].btype,
                                   angls[i].ctype, "    ", angls[j].atype,
                                   angls[j].btype, angls[j].ctype, "    ",
                                   order, 1);
      }
      else if (order == 4 || order == 5) {
        duplicate[j] = TypeCompare(torsions[i].atype, torsions[i].btype,
                                   torsions[i].ctype, torsions[i].dtype,
                                   torsions[j].atype, torsions[j].btype,
                                   torsions[j].ctype, torsions[j].dtype,
                                   order, 1);
        if (duplicate[j] > 0) {
          if (torsions[i].impr != torsions[j].impr ||
              fabs(torsions[i].phase - torsions[j].phase) > 1.0e-4 ||
              fabs(torsions[i].pn - torsions[j].pn) > 1.0e-4) {
            duplicate[j] = 0;
          }
        }
      }
      else if (order == 7) {
        match = 0;
        for (k = 0; k < 5; k++) {
          if (cmaps[i].relResID[k] == cmaps[j].relResID[k]) {
            match++;
          }
        }
        if (match < 5) {
          continue;
        }
        if (CmapResidueOverlap(&cmaps[i], &cmaps[j]) == 0) {
          continue;
        }

        // If we are still here, then these two CMAPs have the same residue
        // applicability, or at least one of them subsumes the other.
        // Check the atom correspondence, and if that's a match, then the
        // two CMAPs are duplicates.
        if (strncmp(cmaps[i].atoms, cmaps[j].atoms, 20) == 0) {
          duplicate[j] = 1;
        }
      }
    }
  }

  // Preserve only the unique parameter specifications
  for (i = 0; i < count; i++) {
    if (duplicate[i] > 0) {
      if (order == 2) {
        free(bonds[i].comment);
      }
      else if (order == 3) {
        free(angls[i].comment);
      }
      else if (order == 4 || order == 5) {
        free(torsions[i].comment);
      }
      else if (order == 7) {
        DestroyXCmapDef(&cmaps[i]);
      }
    }
  }
  nc = 0;
  for (i = 0; i < count; i++) {
    if (duplicate[i] == 0) {
      if (order == 2) {
        bonds[nc] = bonds[i];
      }
      else if (order == 3) {
        angls[nc] = angls[i];
      }
      else if (order == 4 || order == 5) {
        torsions[nc] = torsions[i];
      }
      else if (order == 7) {
        cmaps[nc] = cmaps[i];
      }
      nc++;
    }
  }

  // If this is about proper dihedrals, shift things up
  if (order == 4) {
    j = nc;
    for (i = count; i < abscount; i++) {
      torsions[j] = torsions[i];
      j++;
    }
  }

  // Free allocated memory
  free(duplicate);

  return nc;
}

//-----------------------------------------------------------------------------
// ReadSystemConf: read a system conformation and store it in memory, placing
//                 extra points if needed.                      
//                                                                      
// Arguments:                                                           
//   mp:      the master parameter set                                  
//   nc:      the number of the conformation being read                 
//   mt:      the maximum number of topologies the master parameter set can
//            currently hold                                        
//   tj:      trajectory control information                            
//-----------------------------------------------------------------------------
static void ReadSystemConf(int *nc, prmset *mp, int *mt, trajcon *tj)
{
  int i, j, needprmtop, inptype, fitype, atmcount, frnum;
  double tstamp;
  char *ctmp, *ctm2p;
  prmtop *tp;
  mmsys *myconf;
  dmat traj, evalues;
  cmat allnames;

  // Decide if a new topology needs to be read
  myconf = &mp->conf[*nc];
  needprmtop = 1;
  for (i = 0; i < mp->nunisys; i++) {
    if (strcmp(mp->tpencyc[i].source, myconf->tpsrc) == 0) {
      myconf->tp = &mp->tpencyc[i];
      myconf->GroupNum = i;
      needprmtop = 0;
    }
  }

  // Read a new topology
  if (needprmtop == 1) {

    // Augment the topology with universal settings
    tp = &mp->tpencyc[mp->nunisys];
    strcpy(tp->source, myconf->tpsrc);
    strcpy(tp->WaterName, mp->WaterName);
    tp->ljbuck = mp->ljbuck;
    if (mp->ep[0] != '\0') {
      strcpy(tp->eprulesource, mp->ep);
    }
    else {
      tp->eprulesource[0] = '\0';
    }
  
    // Settle and Rattle are disabled in     
    // parameter fitting.  Nothing is moving. 
    tp->settle = 0;
    tp->rattle = 0;
    tp->lj14fac = mp->lj14fac;
    tp->elec14fac = mp->elec14fac;

    // Read the topology from disk
    GetPrmTop(tp, tj, 1);

    // Placeholders for certain topology fields;   
    // they are not used in this module but do need
    // to be allocated so they can later be freed. 
    tp->lVDWc = (double*)calloc(tp->ntypes, sizeof(double));
    tp->rattlemask = (char*)calloc(MAXNAME, sizeof(char));
    tp->norattlemask = (char*)calloc(MAXNAME, sizeof(char));

    // Set this conformation's topology pointer
    myconf->tp = &mp->tpencyc[mp->nunisys];
    myconf->GroupNum = mp->nunisys;

    // Increment the number of unique systems
    mp->nunisys += 1;

    // Allocate more space if needed
    if (mp->nunisys == *mt) {
      *mt += 32;
      mp->tpencyc = (prmtop*)realloc(mp->tpencyc, (*mt)*sizeof(prmtop));
      for (i = 0; i <= *nc; i++) {
        mp->conf[i].tp = &mp->tpencyc[mp->conf[i].GroupNum];
      }
    }
  }

  // Read the coordinates
  inptype = FindFileType(myconf->crdsrc);
  if (inptype == 0) {
    fitype = ValidCoord(myconf->crdsrc);
    if (fitype == 1 || fitype == 2 || fitype == 4) {
      tstamp = 0.0;
      myconf->crd = ReadRst(myconf->tp, myconf->crdsrc, &tstamp);

      // The energy MUST be given as a number, not as a file name here.
      if (WordIsNumber(myconf->esrc) == 1) {
        myconf->etrg = atof(myconf->esrc);
      }
      else {
        printf("ReadSystemConf >> Error.  %s\nReadSystemConf >> is not a "
               "valid energy for a single structure.\n", myconf->esrc);
        exit(1);
      }
    }
    else if (fitype == 3 || fitype == 5) {

      // Parse the trajectory
      if (fitype == 3) {
        if (mp->verbose) {
          printf("ReadSystemConf >> Warning.  It is not recommended to fit "
                 "parameters based on\nReadSystemConf >> coordinates with "
                 "only three decimal places of precision.\nReadSystemConf >> "
                 "Consider replacing %s\nReadSystemConf >> for system %s.\n",
                 myconf->crdsrc, myconf->tpsrc);
        }
        traj = ReadCrdTraj(myconf->tp, myconf->crdsrc, 0);
      }
      else {
        traj = ReadCDFTraj(myconf->tp, myconf->crdsrc, 0);
      }

      // The array of conformations must now be extended, as this
      // trajectory will likely contain more than just one conformation.
      mp->conf = (mmsys*)realloc(mp->conf,
                                 (mp->nconf+traj.col-1)*sizeof(mmsys));
      myconf = &mp->conf[*nc];
      mp->nconf += traj.col - 1;
      for (i = mp->nconf-1; i > (*nc)+traj.col-1; i--) {
        mp->conf[i] = mp->conf[i-traj.col+1];
      }

      // Save pointers to the original coordinate and topology paths
      ctmp = myconf->crdsrc;
      ctm2p = myconf->tpsrc;

      // Read in the energies for all conformations of this trajectory.
      evalues = ReadListOfDoubles(myconf->esrc);

      // Add new conformations for this trajectory.
      atmcount = myconf->tp->natom;
      for (i = *nc; i < (*nc)+traj.col; i++) {
        mp->conf[i].crdsrc = (char*)malloc(MAXNAME*sizeof(char));
        mp->conf[i].tpsrc = (char*)malloc(MAXNAME*sizeof(char));
        mp->conf[i].esrc = (char*)malloc(MAXNAME*sizeof(char));
        mp->conf[i].tp = myconf->tp;
        mp->conf[i].GroupNum = myconf->GroupNum;

        // Each frame comes from one trajectory, but should be kept distinct.
        // The topology of each conformation from this trajectory is the same.
        sprintf(mp->conf[i].crdsrc, "%s, Frame %d", ctmp, i-(*nc)+1);
        strcpy(mp->conf[i].tpsrc, ctm2p);
        mp->conf[i].etrg = evalues.map[0][i-(*nc)];
        mp->conf[i].crd = CreateCoord(atmcount);
        frnum = i - (*nc);
        for (j = 0; j < 3*atmcount; j++) {
          mp->conf[i].crd.loc[j] = traj.map[j][frnum];
        }
      }
      *nc += traj.col - 1;

      // Free memory before it is lost
      free(ctmp);
      free(ctm2p);
      DestroyDmat(&evalues);
    }
  }
  else if (inptype == 1 || inptype == 2) {

    // Parse the regular expression or directory
    if (inptype == 1) {
      allnames = DirectoryFileSearch(myconf->crdsrc);
    }
    else {
      allnames = RegExpFileSearch(myconf->crdsrc);
    }

    // Files in these locations must be of the energy + coordinates type.
    // No reliable list of energies can be made to match a regular
    // expression search or directory listing, so every file must bundle
    // coordinates and energies together.  The directories could hold
    // other files, however, so be silent and skip files that do not
    // conform to the requirements.  The "energy source" field is taken
    // to be a secondary topology containing the vacuum charge set, in
    // this case.
    for (i = 0; i < allnames.row; i++) {
      fitype = ValidCoord(allnames.map[i]);
      if (fitype != 6) {
        continue;
      }
      
    }
  }
}

//-----------------------------------------------------------------------------
// FreeConformation: free all memory associated with a conformation.    
//                                                                      
// Arguments:                                                           
//   conf:   the conformation to free                                   
//-----------------------------------------------------------------------------
static void FreeConformation(mmsys *conf)
{
  free(conf->crdsrc);
  free(conf->tpsrc);
  free(conf->esrc);
  DestroyCoord(&conf->crd);
  DestroyDmat(&conf->excl);
  DestroyDmat(&conf->nbnrg);
}

//-----------------------------------------------------------------------------
// FindInstanceRanges: loop through the conformations and find the first and
//                     last instances of each system / topology.    
//                                                                      
// Arguments:                                                           
//   mp:      the fitting data (contains a list of all systems)         
//-----------------------------------------------------------------------------
static void FindInstanceRanges(prmset *mp)
{
  int i;

  mp->FirstConf = (int*)malloc(mp->nunisys*sizeof(int));
  mp->LastConf = (int*)malloc(mp->nunisys*sizeof(int));
  SetIVec(mp->FirstConf, mp->nunisys, -1);
  for (i = 0; i < mp->nconf; i++) {
    if (mp->FirstConf[mp->conf[i].GroupNum] == -1) {
      mp->FirstConf[mp->conf[i].GroupNum] = i;
    }
    if (mp->FirstConf[mp->conf[i].GroupNum] >= 0) {
      mp->LastConf[mp->conf[i].GroupNum] = i;
    }
  }
  for (i = 0; i < mp->nunisys; i++) {
    if (mp->FirstConf[i] == -1) {
      printf("mdgx >> Error.  No coordinates for topology %s\n"
             "mdgx >> were specified.\n", mp->tpencyc[i].source);
      exit(1);
    }
  }
}

//-----------------------------------------------------------------------------
// ConfUnitConversion: This function converts the energies to internal units
//                     of kcal/mol, if they are not already in that format.
//
// Arguments:                                                           
//   mp:      the fitting data (contains a list of all systems)         
//-----------------------------------------------------------------------------
static void ConfUnitConversion(prmset *mp)
{
  int i, j;
  char* msg;

  // Determine the energy units
  j = strlen(mp->NrgUnits);
  for (i = 0; i < j; i++) {
    mp->NrgUnits[i] = ToUpper(mp->NrgUnits[i]);
  }
  if (strcmp(mp->NrgUnits, "KCAL") == 0 ||
      strcmp(mp->NrgUnits, "KILOCALORIES") == 0) {
    return;
  }
  if (strcmp(mp->NrgUnits, "HARTREE") == 0 ||
      strcmp(mp->NrgUnits, "ATOMIC") == 0) {
    for (i = 0; i < mp->nconf; i++) {
      mp->conf[i].etrg *= 627.509469;
    }
  }
  else if (strcmp(mp->NrgUnits, "KJ") == 0 ||
           strcmp(mp->NrgUnits, "KILOJOULES") == 0) {
    for (i = 0; i < mp->nconf; i++) {
      mp->conf[i].etrg /= 4.184;
    }
  }
  else if (strcmp(mp->NrgUnits, "J") == 0 ||
           strcmp(mp->NrgUnits, "JOULES") == 0) {
    for (i = 0; i < mp->nconf; i++) {
      mp->conf[i].etrg /= 4184.0;
    }
  }
  else {
    msg = (char*)malloc(512*sizeof(char));
    sprintf(msg, "Error.  Invalid energy units of fitting data %s.",
            mp->NrgUnits);
    printf("\n");
    PrintParagraph(msg, 79, "mdgx >>", stdout);
    free(msg);
  }
}

//-----------------------------------------------------------------------------
// CompEnorm: compute enorm for all systems, a normalized target energy 
//            adjusted to put the average of all target energies for each
//            system at ther average molecular mechanics energy according to
//            the original model.                          
//                                                                      
// Arguments:                                                           
//   mp:     the fitting data set                                       
//-----------------------------------------------------------------------------
static void CompEnorm(prmset *mp)
{
  int i, j, k, ngrp, maxsize;
  int* grpsize;
  double gspread, gmean, gmin, gmax, totwt;
  double* eave;
  double* eadj;
  imat allid;
  dmat allval;

  // Count the populations of each system
  mp->GroupCount = (int*)calloc(mp->nunisys, sizeof(int));
  for (i = 0; i < mp->nconf; i++) {
    mp->GroupCount[mp->conf[i].GroupNum] += 1;
  }

  eave = (double*)calloc(mp->nunisys, sizeof(double));
  eadj = (double*)calloc(mp->nunisys, sizeof(double));
  grpsize = (int*)calloc(mp->nunisys, sizeof(int));
  for (i = 0; i < mp->nconf; i++) {
    ngrp = mp->conf[i].GroupNum;
    eave[ngrp] += mp->conf[i].etrg;
    eadj[ngrp] += mp->conf[i].eorig;
    grpsize[ngrp] += 1;
  }
  maxsize = 0;
  for (i = 0; i < mp->nunisys; i++) {
    eave[i] /= grpsize[i];
    eadj[i] /= grpsize[i];
    maxsize = (grpsize[i] > maxsize) ? grpsize[i] : maxsize;
    grpsize[i] = 0;
  }

  // Compute enorm for all groups
  for (i = 0; i < mp->nconf; i++) {
    ngrp = mp->conf[i].GroupNum;
    mp->conf[i].enorm = mp->conf[i].etrg - eave[ngrp] + eadj[ngrp];
  }

  // Check for outliers and warn if they exist
  allval = CreateDmat(mp->nunisys, maxsize, 0);
  allid = CreateImat(mp->nunisys, maxsize);
  for (i = 0; i < mp->nconf; i++) {
    ngrp = mp->conf[i].GroupNum;
    allval.map[ngrp][grpsize[ngrp]] = mp->conf[i].enorm;
    allid.map[ngrp][grpsize[ngrp]] = i;
    grpsize[ngrp] += 1;
  }
  for (i = 0; i < mp->nunisys; i++) {
    gspread = DStDev(allval.map[i], grpsize[i]);
    gmean = DAverage(allval.map[i], grpsize[i]);
    for (j = 0; j < grpsize[i]; j++) {
      if (fabs(allval.map[i][j] - gmean) / gspread > mp->esigtol) {
        if (mp->verbose == 1) {
          printf("mdgx >> Warning.  Conformation %6d (system %s)\n"
                 "mdgx >> energy is %7.4lf sigma from the mean.\n"
                 "mdgx >> Sigma = %9.4lf, target energy = %9.4lf (mean "
                 "%9.4lf)\n\n", allid.map[i][j],
                 mp->conf[allid.map[i][j]].tp->source,
                 fabs(allval.map[i][j]-gmean) / gspread, gspread,
                 allval.map[i][j], gmean);
        }

        // Destroy this offending conformation
        if (mp->RemoveOutliers == 1) {
          FreeConformation(&mp->conf[allid.map[i][j]]);
          for (k = allid.map[i][j]; k < mp->nconf-1; k++) {
            mp->conf[k] = mp->conf[k+1];
          }
          for (k = 0; k < mp->nunisys*maxsize; k++) {
            if (allid.data[k] > allid.map[i][j]) {
              allid.data[k] -= 1;
            }
          }
          for (k = j; k < grpsize[i]-1; k++) {
            allid.map[i][k] = allid.map[i][k+1];
            allval.map[i][k] = allval.map[i][k+1];
          }
          grpsize[i] -= 1;
          mp->nconf -= 1;
        }
      }
    }
  }

  // Assign weights based on the target energies
  for (i = 0; i < mp->nunisys; i++) {
    if (grpsize[i] == 1) {
      mp->conf[allid.map[i][0]].wt = 1.0;
      continue;
    }
    gmax = DExtreme(allval.map[i], grpsize[i], 1);
    gmin = DExtreme(allval.map[i], grpsize[i], -1);
    gspread = 1.0/(gmax-gmin);
    totwt = 0.0;
    for (j = 0; j < grpsize[i]; j++) {
      mp->conf[allid.map[i][j]].wt = 
        fabs(gmax - allval.map[i][j])*gspread + mp->wtfloor;
      totwt += mp->conf[allid.map[i][j]].wt;
    }
    totwt = grpsize[i]/totwt;
    for (j = 0; j < grpsize[i]; j++) {
      mp->conf[allid.map[i][j]].wt *= totwt;
    }
  }

  // Free allocated memory
  free(eave);
  free(eadj);
  free(grpsize);
  DestroyDmat(&allval);
  DestroyImat(&allid);
}

//-----------------------------------------------------------------------------
// AllocateTermKeyRepositories: allocate the repositories for all term keys.
//                              Each conformation will put its keys in one
//                              line of each of these matrices.
//
// Arguments:
//   mp:       the parameter set containing all conformations and unique
//             topologies
//-----------------------------------------------------------------------------
static void AllocateTermKeyRepositories(prmset *mp)
{
  mp->BondMapIDBucket = (bidx**)malloc(mp->nunisys * sizeof(bidx*));
  mp->AnglMapIDBucket = (aidx**)malloc(mp->nunisys * sizeof(aidx*));
  mp->DiheMapIDBucket = (hidx**)malloc(mp->nunisys * sizeof(hidx*));
  mp->NOpsMapIDBucket = (ridx**)malloc(mp->nunisys * sizeof(ridx*));
  mp->SplnMapIDBucket = (midx**)malloc(mp->nunisys * sizeof(midx*));
  mp->BondMapValBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->BondMapUxBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->BondMapUyBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->BondMapUzBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->BondMapUwBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->BondMapCntrbBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->AnglMapValBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->AnglMapUxBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->AnglMapUyBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->AnglMapCntrbBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->DiheMapValBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->DiheMapUBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->DiheMapCntrbBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->NOpsMapValBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->NOpsMapKernelBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->NOpsMapCntrbBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->SplnMapValBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
  mp->SplnMapCntrbBucket = (dmat*)malloc(mp->nunisys * sizeof(dmat));
}

//-----------------------------------------------------------------------------
// AllocateTermKey: allocate memory for a conformation's term key.
//
// Arguments:
//   mp:       the parameter set (needed for the case of CMAPs, when the map
//             resolution may not be known
//   confID    the conformation ID number (index into mp->conf)
//   order:    the order of the terms, 2 for bonds, 3 for angles, 4 for
//             dihedral terms, 6 for NMR operations, and 7 for CMAPs
//   nterm:    the number of terms found with the specified order
//   nconf:    the number of conformations that share the same topology and
//             thus will have related term keys
//   doalloc:  do the memory allocation for the requisite matrix.  Each row of
//             this matrix will serve as the term key (of appropriate order)
//             for one conformation out of all the conformations that share a
//             particular topology.  If this flag is set to zero, the matrix
//             is assumed to have already been allocated.  The matrices for
//             all orders of terms and all topologies are stored in arrays of
//             the parameter set.
//-----------------------------------------------------------------------------
static void AllocateTermKey(prmset *mp, int confID, int order, int nterm,
                            int nconf, int confPos)
{
  int i, topID, ksize, mapdim, tmpterm;
  mmsys *conf;

  // The conformation shows where these allocations and data will go
  conf = &mp->conf[confID];
  topID = conf->GroupNum;
  
  // Return if there are no terms to allocate, after making the appropriate
  // marks on the matrix struct.
  tmpterm = nterm;
  nterm += (nterm == 0);
  if (confPos == 0) {
    if (order == 2) {
      mp->BondMapIDBucket[topID] = (bidx*)malloc(nterm * nconf * sizeof(bidx));
      mp->BondMapValBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->BondMapUxBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->BondMapUyBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->BondMapUzBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->BondMapUwBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->BondMapCntrbBucket[topID] = CreateDmat(nconf, nterm, 0);
    }
    else if (order == 3) {
      mp->AnglMapIDBucket[topID] = (aidx*)malloc(nterm * nconf * sizeof(aidx));
      mp->AnglMapValBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->AnglMapUxBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->AnglMapUyBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->AnglMapCntrbBucket[topID] = CreateDmat(nconf, nterm, 0);
    }
    else if (order == 4) {
      mp->DiheMapIDBucket[topID] = (hidx*)malloc(nterm * nconf * sizeof(hidx));
      mp->DiheMapValBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->DiheMapUBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->DiheMapCntrbBucket[topID] = CreateDmat(nconf, nterm, 0);
    }
    else if (order == 6) {
      mp->NOpsMapIDBucket[topID] = (ridx*)malloc(nterm * nconf * sizeof(ridx));
      mp->NOpsMapValBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->NOpsMapKernelBucket[topID] = CreateDmat(nconf, nterm, 0);
      mp->NOpsMapCntrbBucket[topID] = CreateDmat(nconf, nterm, 0);
    }
    else if (order == 7) {
      mp->SplnMapIDBucket[topID] = (midx*)malloc(nterm * nconf * sizeof(midx));
      mp->SplnMapValBucket[topID] = CreateDmat(nconf, 2 * nterm, 0);
      mp->SplnMapCntrbBucket[topID] = CreateDmat(nconf, nterm, 0);
    }
  }
  nterm = tmpterm;
  
  // Switch to deal with each kind of term
  if (order == 2) {
    conf->bmap.id = &mp->BondMapIDBucket[topID][nterm * confPos];
    conf->bmap.val = mp->BondMapValBucket[topID].map[confPos];
    conf->bmap.UkernelX = mp->BondMapUxBucket[topID].map[confPos];
    conf->bmap.UkernelY = mp->BondMapUyBucket[topID].map[confPos];
    conf->bmap.UkernelZ = mp->BondMapUzBucket[topID].map[confPos];
    conf->bmap.UkernelW = mp->BondMapUwBucket[topID].map[confPos];
    conf->bmap.Ucontrib = mp->BondMapCntrbBucket[topID].map[confPos];
    conf->bmap.nbond = nterm;
  }
  else if (order == 3) {
    conf->amap.id = &mp->AnglMapIDBucket[topID][nterm * confPos];
    conf->amap.val = mp->AnglMapValBucket[topID].map[confPos];
    conf->amap.UkernelX = mp->AnglMapUxBucket[topID].map[confPos];
    conf->amap.UkernelY = mp->AnglMapUyBucket[topID].map[confPos];
    conf->amap.Ucontrib = mp->AnglMapCntrbBucket[topID].map[confPos];
    conf->amap.nangl = nterm;
  }
  else if (order == 4) {
    conf->hmap.id = &mp->DiheMapIDBucket[topID][nterm * confPos];
    conf->hmap.val = mp->DiheMapValBucket[topID].map[confPos];
    conf->hmap.Ukernel = mp->DiheMapUBucket[topID].map[confPos];
    conf->hmap.Ucontrib = mp->DiheMapCntrbBucket[topID].map[confPos];
    conf->hmap.ntterm = nterm;
  }
  else if (order == 6) {
    conf->rmap.id = &mp->NOpsMapIDBucket[topID][nterm * confPos];
    conf->rmap.val = mp->NOpsMapValBucket[topID].map[confPos];
    conf->rmap.Ukernel = mp->NOpsMapKernelBucket[topID].map[confPos];
    conf->rmap.Ucontrib = mp->NOpsMapCntrbBucket[topID].map[confPos];
    conf->rmap.nops = nterm;
  }
  else if (order == 7) {
    conf->mmap.id = &mp->SplnMapIDBucket[topID][nterm * confPos];
    conf->mmap.val = mp->SplnMapValBucket[topID].map[confPos];
    conf->mmap.Ucontrib = mp->SplnMapCntrbBucket[topID].map[confPos];
    conf->mmap.ncmap = nterm;
  }
}

//-----------------------------------------------------------------------------
// ReflectTermKey: copy all values of one conformation's term key map to
//                 another conformation.                                
//                                                                      
// Arguments:                                                           
//   mp:     the set of systems and parameter guidelines                
//   confi:  the conformation with a complete key                       
//   confj:  the conformation with an allocated but empty key           
//   order:  the order of the terms, 2 for bonds, 3 for angles, and 4   
//           for dihedral terms                                         
//-----------------------------------------------------------------------------
static void ReflectTermKey(prmset *mp, mmsys *confi, mmsys *confj, int order)
{
  int i, nterm;

  if (order == 2) {
    nterm = confi->bmap.nbond;
    for (i = 0; i < nterm; i++) {
      confj->bmap.id[i] = confi->bmap.id[i];
      mp->bonds[confi->bmap.id[i].key].ninst += 1;
    }
  }
  else if (order == 3) {
    nterm = confi->amap.nangl;
    for (i = 0; i < nterm; i++) {
      confj->amap.id[i] = confi->amap.id[i];
      mp->angls[confi->amap.id[i].key].ninst += 1;
    }
  }
  else if (order == 4) {
    nterm = confi->hmap.ntterm;
    for (i = 0; i < nterm; i++) {
      confj->hmap.id[i] = confi->hmap.id[i];
      mp->torsions[confi->hmap.id[i].key].ninst += 1;
    }
  }
  else if (order == 6) {
    nterm = confi->rmap.nops;
    for (i = 0; i < nterm; i++) {
      confj->rmap.id[i] = confi->rmap.id[i];
      mp->nmrops[confi->rmap.id[i].key].ninst += 1;
    }
  }
  else if (order == 7) {
    nterm = confi->mmap.ncmap;
    for (i = 0; i < nterm; i++) {
      confj->mmap.id[i] = confi->mmap.id[i];
      mp->cmaps[confi->mmap.id[i].key].ninst += 1;
    }
  }
}

//-----------------------------------------------------------------------------
// FindCentralAtom: function for finding the central atom of an improper
//                  dihedral taken from a topology file.  The problem is that
//                  improper dihedrals in topologies can store the central
//                  atom in either the second or third position.
//                                                                      
// Arguments:                                                           
//   tp:     the topology                                               
//   atidx:  the index of the atom in the topology that controls the    
//           improper dihedral of interest (this may or may not actually
//           be the central atom--this integer is passed in just to     
//           locate the dihedral).                                      
//   impidx: the index of the improper dihedral in the list of them     
//           controlled by the atom                                     
//-----------------------------------------------------------------------------
static int FindCentralAtom(prmtop *tp, int atidx, int impidx)
{
  int i, j, k;
  int aidx[4], bondsto[4][4];

  // Unpack atoms in the improper dihedral
  aidx[0] = tp->HLC[atidx].HC[impidx].a;
  aidx[1] = atidx;
  aidx[2] = tp->HLC[atidx].HC[impidx].c;
  aidx[3] = tp->HLC[atidx].HC[impidx].d;

  // Tally up which of these four atoms bond to others
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      bondsto[i][j] = 0;
    }
  }
  for (i = 0; i < 4; i++) {
    for (j = 0; j < tp->BLC[aidx[i]].nbond; j++) {
      for (k = 0 ; k < 4; k++) {
        if (tp->BLC[aidx[i]].BC[j].b == aidx[k]) {
          bondsto[i][k] = 1;
          bondsto[k][i] = 1;
        }
      }
    }
  }
  for (i = 0; i < 4; i++) {
    if (bondsto[i][0] + bondsto[i][1] + bondsto[i][2] + bondsto[i][3] == 3) {
      return i;
    }
  }

  // Error message
  printf("FindCentralAtom >> Warning.  No atom in dihedral %4.4s %4.4s %4.4s "
         "%4.4s\nFindCentralAtom >> could be identified as the center.\n",
         &tp->AtomNames[4*tp->HLC[atidx].HC[impidx].a],
         &tp->AtomNames[4*atidx],
         &tp->AtomNames[4*tp->HLC[atidx].HC[impidx].c],
         &tp->AtomNames[4*tp->HLC[atidx].HC[impidx].d]);

  return 2;
}

//-----------------------------------------------------------------------------
// PruneMatchMatrix: remove identical rows from the matrix of matches.  This is
//                   done for counting the number of actual NMR operations.
//
// Arguments:
//   M:       the match matrix
//-----------------------------------------------------------------------------
static void PruneMatchMatrix(imat *M)
{
  int i, j, k, match;
  int* preserve;

  preserve = (int*)malloc(M->row*sizeof(int));
  SetIVec(preserve, M->row, 1);
  for (i = 0; i < M->row; i++) {
    if (preserve[i] == 0) {
      continue;
    }
    if (ISum(M->map[i], M->col) == 0) {
      preserve[i] = 0;
      continue;
    }
    for (j = i+1; j < M->row; j++) {
      if (preserve[j] == 0) {
        continue;
      }
      match = 0;
      if (M->col == 2) {
        if ((M->map[i][0] == M->map[j][0] && M->map[i][1] == M->map[j][1]) ||
            (M->map[j][0] == M->map[i][1] && M->map[j][1] == M->map[i][0])) {
          match = 1;
        }
      }
      else if (M->col == 3 && M->map[i][1] == M->map[j][1]) {
        if ((M->map[i][0] == M->map[j][0] && M->map[i][2] == M->map[j][2]) ||
            (M->map[j][0] == M->map[i][2] && M->map[j][2] == M->map[i][0])) {
          match = 1;
        }
      }
      else if (M->col == 4) {
        if ((M->map[i][0] == M->map[j][0] && M->map[i][1] == M->map[j][1] &&
             M->map[i][2] == M->map[j][2] && M->map[i][3] == M->map[j][3]) ||
            (M->map[i][0] == M->map[j][3] && M->map[i][1] == M->map[j][2] &&
             M->map[i][2] == M->map[j][1] && M->map[i][3] == M->map[j][0])) {
          match = 1;
        }
      }
      if (match == 1) {
        preserve[j] = 0;
      }
    }
  }
  j = 0;
  for (i = 0; i < M->row; i++) {
    if (preserve[i] == 1) {
      for (k = 0; k < M->col; k++) {
        M->map[j][k] = M->map[i][k];
      }
      j++;
    }
  }
  if (j == 0) {
    *M = ReallocImat(M, 1, M->col);
    M->row = 0;
  }
  else {
    *M = ReallocImat(M, j, M->col);
  }
}

//-----------------------------------------------------------------------------
// BuildMatchMatrix: build the matrix of matches for a particular system and
//                   a particulat NMR restraint.
//
// Arguments:
//   myconf: conformation within the fitting data
//   myop:   the operation of interest
//-----------------------------------------------------------------------------
imat BuildMatchMatrix(mmsys *myconf, nmroper *myop)
{
  int h, i, j, k, nmatch;
  int* amask;
  int* bmask;
  int* cmask;
  int* dmask;
  imat matches;
  prmtop *tp;

  // Map each atom mask and sibling
  tp = myconf->tp;
  matches = CreateImat(16, myop->order);
  nmatch = 0;
  for (h = 0; h < myop->nsibling; h++) {
    amask = ParseAmbMask(myop->amask.map[h], tp, &myconf->crd);
    bmask = ParseAmbMask(myop->bmask.map[h], tp, &myconf->crd);
    if (myop->order >= 3) {
      cmask = ParseAmbMask(myop->cmask.map[h], tp, &myconf->crd);
    }
    if (myop->order == 4) {
      dmask = ParseAmbMask(myop->dmask.map[h], tp, &myconf->crd);
    }

    // Loop over all bonded parameters and try to find matches among the
    // active atoms in each of the masks.
    for (i = 0; i < tp->natom; i++) {
      if (myop->order == 2) {
        for (j = 0; j < tp->BLC[i].nbond; j++) {
          if ((amask[tp->BLC[i].BC[j].a] == 1 &&
               bmask[tp->BLC[i].BC[j].b] == 1) ||
              (amask[tp->BLC[i].BC[j].b] == 1 &&
               bmask[tp->BLC[i].BC[j].a] == 1)) {
            matches.map[nmatch][0] = tp->BLC[i].BC[j].a;
            matches.map[nmatch][1] = tp->BLC[i].BC[j].b;
            nmatch++;
            if (nmatch == matches.row) {
              matches = ReallocImat(&matches, nmatch*2, myop->order);
            }
          }
        }
        for (j = 0; j < tp->ALC[i].nangl; j++) {
          if ((amask[tp->ALC[i].AC[j].a] == 1 &&
               bmask[tp->ALC[i].AC[j].c] == 1) ||
              (amask[tp->ALC[i].AC[j].c] == 1 &&
               bmask[tp->ALC[i].AC[j].a] == 1)) {
            matches.map[nmatch][0] = tp->ALC[i].AC[j].a;
            matches.map[nmatch][1] = tp->ALC[i].AC[j].c;
            nmatch++;
            if (nmatch == matches.row) {
              matches = ReallocImat(&matches, nmatch*2, myop->order);
            }
          }
        }
        for (j = 0; j < tp->HLC[i].ndihe; j++) {
          if (tp->HLC[i].HC[j].impr == 1 &&
              ((amask[tp->HLC[i].HC[j].a] == 1 &&
                bmask[tp->HLC[i].HC[j].d] == 1) ||
               (amask[tp->HLC[i].HC[j].d] == 1 &&
                bmask[tp->HLC[i].HC[j].a] == 1))) {
            matches.map[nmatch][0] = tp->HLC[i].HC[j].a;
            matches.map[nmatch][1] = tp->HLC[i].HC[j].d;
            nmatch++;
            if (nmatch == matches.row) {
              matches = ReallocImat(&matches, nmatch*2, myop->order);
            }
          }
        }
      }
      else if (myop->order == 3) {
        for (j = 0; j < tp->ALC[i].nangl; j++) {
          if (bmask[i] == 1 && (amask[tp->ALC[i].AC[j].a] == 1 &&
                                cmask[tp->ALC[i].AC[j].c] == 1) ||
                               (amask[tp->ALC[i].AC[j].c] == 1 &&
                                cmask[tp->ALC[i].AC[j].a] == 1)) {
            matches.map[nmatch][0] = tp->ALC[i].AC[j].a;
            matches.map[nmatch][1] = i;
            matches.map[nmatch][2] = tp->ALC[i].AC[j].c;
            nmatch++;
            if (nmatch == matches.row) {
              matches = ReallocImat(&matches, nmatch*2, myop->order);
            }
          }
        }
      }
      else if (myop->order == 4) {
        for (j = 0; j < tp->HLC[i].ndihe; j++) {
          if (tp->HLC[i].HC[j].impr == 1) {
            continue;
          }
          if ((amask[tp->HLC[i].HC[j].a] == 1 && bmask[i] == 1 &&
               cmask[tp->HLC[i].HC[j].c] == 1 &&
               dmask[tp->HLC[i].HC[j].d] == 1) ||
              (amask[tp->HLC[i].HC[j].d] == 1 &&
               bmask[tp->HLC[i].HC[j].c] == 1 && cmask[i] == 1 &&
               dmask[tp->HLC[i].HC[j].a] == 1)) {
            matches.map[nmatch][0] = tp->HLC[i].HC[j].a;
            matches.map[nmatch][1] = i;
            matches.map[nmatch][2] = tp->HLC[i].HC[j].c;
            matches.map[nmatch][3] = tp->HLC[i].HC[j].d;
            nmatch++;
            if (nmatch == matches.row) {
              matches = ReallocImat(&matches, nmatch*2, myop->order);
            }
          }
        }
      }
    }
  }
  PruneMatchMatrix(&matches);

  // Organize the match matrix to put the lowest index first in each row
  j = matches.col-1;
  for (i = 0; i < matches.row; i++) {
    if (matches.map[i][0] > matches.map[i][j]) {
      SWAP(matches.map[i][0], matches.map[i][j], k);

      // Nothing more needs to be done for order = 3, but for
      // order = 4 there is another pair of numbers to swap.
      if (myop->order == 4) {
        SWAP(matches.map[i][1], matches.map[i][j-1], k);
      }
    }
  }

  // Free allocated memory
  free(amask);
  free(bmask);
  if (myop->order >= 3) {
    free(cmask);
  }
  if (myop->order == 4) {
    free(dmask);
  }

  return matches;
}

//-----------------------------------------------------------------------------
// CountNMROperations: count the number of optimizable NMR operations in a
//                     system for the purposes of parameter fitting. 
//
// Arguments:
//   mp:       the set of systems and parameter guidelines                
//   confidx:  the index of the conformation
//-----------------------------------------------------------------------------
static int CountNMROperations(prmset *mp, int confidx)
{
  int i, nterm, nmatch;
  imat matches;

  // Unpack this system
  nterm = 0;

  // Loop over all NMR operations
  for (i = 0; i < mp->nops; i++) {
    matches = BuildMatchMatrix(&mp->conf[confidx], &mp->nmrops[i]);
    nterm += matches.row;
    DestroyImat(&matches);
  }

  return nterm;
}

//-----------------------------------------------------------------------------
// MapNMROperations: map the various NMR operations that influence the energy
//                   of a conformation.
//
// Arguments:
//   mp:       the set of systems and parameter guidelines
//   confidx:  the index of the conformation currently under examination
//-----------------------------------------------------------------------------
static void MapNMROperations(prmset *mp, int confidx)
{
  int i, j, nterm, nmatch;
  imat matches;
  prmtop *tp;
  mmsys *cf;

  cf = &mp->conf[confidx];
  tp = cf->tp;
  nterm = 0;
  for (i = 0; i < mp->nops; i++) {
    matches = BuildMatchMatrix(&mp->conf[confidx], &mp->nmrops[i]);
    for (j = 0; j < matches.row; j++) {
      cf->rmap.id[nterm].order = mp->nmrops[i].order;
      cf->rmap.id[nterm].key = i;
      cf->rmap.id[nterm].a = matches.map[j][0];
      cf->rmap.id[nterm].b = matches.map[j][1];
      if (mp->nmrops[i].order >= 3) {
        cf->rmap.id[nterm].c = matches.map[j][2];
      }
      if (mp->nmrops[i].order == 4) {
        cf->rmap.id[nterm].d = matches.map[j][3];
      }
      nterm++;
    }
    DestroyImat(&matches);
  }
}

//-----------------------------------------------------------------------------
// MapCmapTerms: map the CMAPs to the original, unique surfaces read in from
//               parameter (i.e. parm##.dat and frcmod) files.
//
// Arguments:
//   mp:     the set of systems and parameter guidelines
//   confidx:  the index of the conformation currently under examination
//-----------------------------------------------------------------------------
static void MapCmapTerms(prmset *mp, int confidx)
{
  int i, j, k, nterm, centralResID, tpResScore, found;
  int tpResID[5];
  char atmnames[24];
  prmtop *tp;
  mmsys *cf;

  cf = &mp->conf[confidx];
  tp = cf->tp;
  nterm = 0;
  for (i = 0; i < tp->ncmap; i++) {

    // The atoms are easy enough, and since the atom types are
    // irrelevant to the CMAP there are no such type comparisons.
    cf->mmap.id[i].a   = tp->Mterms[i].a;
    cf->mmap.id[i].b   = tp->Mterms[i].b;
    cf->mmap.id[i].c   = tp->Mterms[i].c;
    cf->mmap.id[i].d   = tp->Mterms[i].d;
    cf->mmap.id[i].e   = tp->Mterms[i].e;

    // However, finding the exact CMAP surface to which this
    // corresponds will present its own minor challenges.  Get
    // the residue number profile, then the name of the central
    // residue, then the names of all atoms.
    tpResID[0] = LocateResID(tp, tp->Mterms[i].a, 0, tp->nres);
    tpResID[1] = LocateResID(tp, tp->Mterms[i].b, 0, tp->nres);
    tpResID[2] = LocateResID(tp, tp->Mterms[i].c, 0, tp->nres);
    tpResID[3] = LocateResID(tp, tp->Mterms[i].d, 0, tp->nres);
    tpResID[4] = LocateResID(tp, tp->Mterms[i].e, 0, tp->nres);
    centralResID = tpResID[2];
    for (j = 0; j < 5; j++) {
      tpResID[j] -= centralResID;
    }

    // Loop over all CMAPs found in parameter files.  Compare the
    // residue profiles, central residue names, and atom names to
    // determine a match.
    found = 0;
    for (j = 0; j < mp->ncmaps; j++) {

      // If the residue profiles are not a 100% match, bail out.
      tpResScore = 0;
      for (k = 0; k < 5; k++) {
        tpResScore += (tpResID[k] ==  mp->cmaps[j].relResID[k]);
      }
      if (tpResScore < 5) {
        continue;
      }

      // If the central residue name is not found in the list for this
      // CMAP, bail out.
      tpResScore = 0;
      for (k = 0; k < mp->cmaps[j].resnames.row; k++) {
        if (str4cmp(&tp->ResNames[4 * centralResID],
                    mp->cmaps[j].resnames.map[k]) == 0) {
          tpResScore = 1;
        }
      }
      if (tpResScore == 0) {
        continue;
      }

      // If the atom names match, then this is the CMAP we want.
      tpResScore = str4cmp(&tp->AtomNames[4 * tp->Mterms[i].a],
                           &mp->cmaps[j].atoms[0]) +
                   str4cmp(&tp->AtomNames[4 * tp->Mterms[i].b],
                           &mp->cmaps[j].atoms[4]) +
                   str4cmp(&tp->AtomNames[4 * tp->Mterms[i].c],
                           &mp->cmaps[j].atoms[8]) +
                   str4cmp(&tp->AtomNames[4 * tp->Mterms[i].d],
                           &mp->cmaps[j].atoms[12]) +
                   str4cmp(&tp->AtomNames[4 * tp->Mterms[i].e],
                           &mp->cmaps[j].atoms[16]);
      if (tpResScore != 0) {
        continue;
      }

      // Mark this CMAP's origin and indicate that it has been found
      mp->conf[confidx].mmap.id[nterm].key = j;
      mp->cmaps[j].ninst += 1;
      found = 1;
    }
    if (found == 0) {
      printf("MakeTermKey >> Error.  A CMAP spanning atoms %4.4s %4.4s %4.4s "
             "%4.4s %4.4s\nMakeTermKey >> with central residue %4.4s is not "
             "in the parameter files.\n", &tp->AtomNames[4 * tp->Mterms[i].a],
             &tp->AtomNames[4 * tp->Mterms[i].b],
             &tp->AtomNames[4 * tp->Mterms[i].c],
             &tp->AtomNames[4 * tp->Mterms[i].d],
             &tp->AtomNames[4 * tp->Mterms[i].e],
             &tp->ResNames[4 * centralResID]);
    }
    mp->conf[confidx].mmap.id[nterm].a = tp->Mterms[i].a;
    mp->conf[confidx].mmap.id[nterm].b = tp->Mterms[i].b;
    mp->conf[confidx].mmap.id[nterm].c = tp->Mterms[i].c;
    mp->conf[confidx].mmap.id[nterm].d = tp->Mterms[i].d;
    mp->conf[confidx].mmap.id[nterm].e = tp->Mterms[i].e;
    nterm++;
  }
}

//-----------------------------------------------------------------------------
// MakeTermKey: make a table of unique bonded interactions across all systems,
//              and in so doing create a table of bonded terms for each system
//              with indexing into the master key.      
//
// Arguments:                                                           
//   mp:     the set of systems and parameter guidelines
//   order:  the order of the terms, 2 for bonds, 3 for angles, 4 for
//           dihedral terms, 7 for CMAPs, and 6 for NMR operations
//-----------------------------------------------------------------------------
static void MakeTermKey(prmset *mp, int order)
{
  int i, j, k, m, n, ii, jmin, klim, nmpterm, nterm, found, norder, isimpr;
  int atmA, atmB, atmC, atmD, resA, resB, resC, resD;
  int compval, maxcomp;
  int* havekey;
  double diheK, dihePhi, diheN;
  char typeA[8], typeB[8], typeC[8], typeD[8], imprflag[32];
  char otypeA[8], otypeB[8], otypeC[8], otypeD[8];
  prmtop *tp;

  // Count the number of bonded parameters in the fitting set's master key
  if (order == 2) {
    nmpterm = mp->nbond;
    for (i = 0; i < nmpterm; i++) {
      mp->bonds[i].ninst = 0;
    }
  }
  else if (order == 3) {
    nmpterm = mp->nangl;
    for (i = 0; i < nmpterm; i++) {
      mp->angls[i].ninst = 0;
    }
  }
  else if (order == 4) {
    nmpterm = mp->ntor;
    for (i = 0; i < nmpterm; i++) {
      mp->torsions[i].ninst = 0;
    }
  }
  else if (order == 6) {
    nmpterm = mp->nops;
    for (i = 0; i < nmpterm; i++) {
      mp->nmrops[i].ninst = 0;
    }
  }
  else if (order == 7) {
    nmpterm = mp->ncmaps;
    for (i = 0; i < nmpterm; i++) {
      mp->cmaps[i].ninst = 0;
    }
  }

  // Allocate an array to track which conformations have keys
  havekey = (int*)calloc(mp->nconf, sizeof(int));

  // Allocate space for keys in all conformations
  for (i = 0; i < mp->nconf; i++) {

    // Skip if this key is already allocated
    if (havekey[i] == 1) {
      continue;
    }

    // Unpack this system
    tp = mp->conf[i].tp;
    nterm = 0;
    if (order == 2) {
      for (j = 0; j < tp->natom; j++) {
        nterm += tp->BLC[j].nbond;
      }
    }
    else if (order == 3) {
      for (j = 0; j < tp->natom; j++) {
        nterm += tp->ALC[j].nangl;
      }
    }
    else if (order == 4) {
      for (j = 0; j < tp->natom; j++) {
        for (k = 0; k < tp->HLC[j].ndihe; k++) {
          nterm += tp->HLC[j].HC[k].nt;
        }
      }
    }
    else if (order == 6) {
      nterm = CountNMROperations(mp, i);
    }
    else if (order == 7) {
      nterm += tp->ncmap;
    }

    // Search for other conformations of the same topology
    found = mp->conf[i].GroupNum;
    k = 1;    
    for (j = i+1; j < mp->nconf; j++) {
      k += (havekey[j] == 0 && mp->conf[j].GroupNum == found);
    }
    AllocateTermKey(mp, i, order, nterm, k, 0);
    havekey[i] = 1;

    // Mark other conformations and redirect their term keys
    // to lines of the matrix 
    found = mp->conf[i].GroupNum;
    k = 1;
    for (j = i+1; j < mp->nconf; j++) {
      if (havekey[j] == 0 && mp->conf[j].GroupNum == found) {
        AllocateTermKey(mp, j, order, nterm, k, k);
        k++;
        havekey[j] = 1;
      }
    }
  }

  // Fill out the keys
  for (i = 0; i < mp->nconf; i++) {

    // Skip if this conformation already has a key
    if (havekey[i] == 2) {
      continue;
    }

    // Unpack this system; skip the loop if NMR operations or CMAPs
    // are the subject (these work on names of atoms or residues,
    // not atom types)
    tp = mp->conf[i].tp;
    nterm = 0;
    jmin = (order == 6 || order == 7) ? tp->natom : 0;
    for (j = jmin; j < tp->natom; j++) {
      klim = (order == 2) ? tp->BLC[j].nbond :
        (order == 3) ? tp->ALC[j].nangl : tp->HLC[j].ndihe;
      for (k = 0; k < klim; k++) {
        if (order == 2) {
          atmA = tp->BLC[j].BC[k].a;
          atmB = tp->BLC[j].BC[k].b;
        }
        else if (order == 3) {
          atmA = tp->ALC[j].AC[k].a;
          atmB = j;
          atmC = tp->ALC[j].AC[k].c;
        }
        else if (order == 4) {
          atmA = tp->HLC[j].HC[k].a;
          atmB = j;
          atmC = tp->HLC[j].HC[k].c;
          atmD = tp->HLC[j].HC[k].d;
        }
        strncpy(typeA, &tp->AtomTypes[4*atmA], 4);
        strncpy(typeB, &tp->AtomTypes[4*atmB], 4);
        if (order > 2) {
          strncpy(typeC, &tp->AtomTypes[4*atmC], 4);
        }
        if (order == 4) {
          strncpy(typeD, &tp->AtomTypes[4*atmD], 4);
        }

        // Reference this bonded term against the master key
        found = 0;
        maxcomp = 0;
        if (order == 2) {
          for (m = 0; m < nmpterm; m++) {
            compval = TypeCompare(typeA, typeB, "    ", "    ",
                                  mp->bonds[m].atype, mp->bonds[m].btype,
                                  "    ", "    ", 2, 0);
            if (compval > maxcomp) {

              // Undo the previous assignment if one exists
              if (maxcomp == 1) {
                mp->bonds[mp->conf[i].bmap.id[nterm].key].ninst -= 1;
              }
              maxcomp = compval;

              // Assign or reassign the bond to its place in the key
              mp->conf[i].bmap.id[nterm].key = m;
              mp->bonds[m].ninst += 1;
              found = 1;
            }
          }
          if (found == 0) {

            // This bond is not present in the parameter set
            printf("MakeTermKey >> Error.  A bond between %.4s and %.4s is "
                   "not found in the\nMakeTermKey >> parameter files.\n",
                   typeA, typeB);
            exit(1);
          }
          mp->conf[i].bmap.id[nterm].a = atmA;
          mp->conf[i].bmap.id[nterm].b = atmB;
          nterm++;
        }
        else if (order == 3) {
          for (m = 0; m < nmpterm; m++) {
            compval = TypeCompare(typeA, typeB, typeC, "    ",
                                   mp->angls[m].atype, mp->angls[m].btype,
                                  mp->angls[m].ctype, "    ", 3, 0);
            if (compval > maxcomp) {

              // Undo the previous assignment if one exists
              if (maxcomp == 1) {
                mp->angls[mp->conf[i].amap.id[nterm].key].ninst -= 1;
              }
              maxcomp = compval;

              // Assign or reassign the angle to its place in the key
              mp->conf[i].amap.id[nterm].key = m;
              mp->angls[m].ninst += 1;
              found = 1;
            }
          }
          if (found == 0) {

            // This angle is not present in the parameter set
            printf("MakeTermKey >> Error.  An angle between %.4s, %.4s, and "
                   "%.4s is not found in the\nMakeTermKey >> parameter files."
                   "\n", typeA, typeB, typeC);
            exit(1);
          }
          mp->conf[i].amap.id[nterm].a = atmA;
          mp->conf[i].amap.id[nterm].b = atmB;
          mp->conf[i].amap.id[nterm].c = atmC;
          nterm++;
        }
        else if (order == 4) {
          isimpr = tp->HLC[j].HC[k].impr;
          norder = order + isimpr;
          for (ii = 0; ii < 4; ii++) {
            otypeA[ii] = typeA[ii];
            otypeB[ii] = typeB[ii];
            otypeC[ii] = typeC[ii];
            otypeD[ii] = typeD[ii];
          }

          // Find the central atom of the four, and order the    
          // atoms apprpriately, if this is an improper dihedral.
          if (isimpr == 1) {
            m = FindCentralAtom(tp, j, k);

            // Order the atom types A, B, and D if this is to
            // be compared to an improper.                   
            if (m == 0) {
              for (ii = 0; ii < 4; ii++) {
                otypeC[ii] = typeA[ii];
                otypeA[ii] = typeC[ii];
              }
            }
            else if (m == 1) {
              for (ii = 0; ii < 4; ii++) {
                otypeB[ii] = typeC[ii];
                otypeC[ii] = typeB[ii];
              }
            }
            else if (m == 3) {
              for (ii = 0; ii < 4; ii++) {
                otypeC[ii] = typeD[ii];
                otypeD[ii] = typeC[ii];
              }
            }
            BubbleChar4(otypeA, otypeB);
            BubbleChar4(otypeB, otypeD);
            BubbleChar4(otypeA, otypeB);
          }

          // Loop over all terms
          for (m = 0; m < tp->HLC[j].HC[k].nt; m++) {

            // For dihedrals, there is another layer of indexing.
            // We can identify the stiffness K, phase angle Phi, 
            // and periodicity N of the mth torsion term of this 
            // kth dihedral fourier series controlled by atom j  
            // of system i.                                       
            diheK = tp->HParam[tp->HLC[j].HC[k].t[m]].K;
            dihePhi = tp->HParam[tp->HLC[j].HC[k].t[m]].Phi;
            diheN = tp->HParam[tp->HLC[j].HC[k].t[m]].N;

            // Now scan through the list of unique
            // torsion terms identified thus far. 
            found = 0;
            maxcomp = 0;
            for (n = 0; n < nmpterm; n++) {
              if (mp->torsions[n].impr != isimpr) {
                continue;
              }

              // Finally: make the comparison of atom types.
              compval = TypeCompare(otypeA, otypeB, otypeC, otypeD,
                                    mp->torsions[n].atype,
                                    mp->torsions[n].btype,
                                    mp->torsions[n].ctype,
                                    mp->torsions[n].dtype, norder, 0);
              if (compval > maxcomp) {

                // This conditional is nested for clarity
                if (fabs(mp->torsions[n].phase - dihePhi) < 1.0e-4 &&
                    fabs(mp->torsions[n].pn - diheN) < 1.0e-4) {

                  // Undo the previous assignment if one exists
                  if (maxcomp == 1) {
                    mp->torsions[mp->conf[i].hmap.id[nterm].key].ninst -= 1;
                  }
                  maxcomp = compval;

                  // Assign the torsion to its place in the key
                  mp->conf[i].hmap.id[nterm].key = n;
                  mp->torsions[n].ninst += 1;
                  found = 1;
                }
              }
            }
            if (found == 0) {

              // This torsion is not present in the parameter set
              if (isimpr == 1) {
                sprintf(imprflag, "n improper");
              }
              else {
                imprflag[0] = '\0';
              }
              printf("MakeTermKey >> Error.  A%s torsion term for %.4s %.4s "
                     "%.4s %.4s is not\nMakeTermKey >> found in the parameter "
                     "files.  Amplitude %9.4lf, periodicity\nMakeTermKey >> "
                     "%9.4lf, and phase angle %9.4lf were required.  Data "
                     "point\nMakeTermKey >> %s is not completely represented "
                     "by the parameter files.\n", imprflag, typeA, typeB,
                     typeC, typeD, diheK, diheN, dihePhi, mp->conf[i].crdsrc);
              resA = LocateResID(tp, atmA, 0, tp->nres);
              resB = LocateResID(tp, atmB, 0, tp->nres);
              resC = LocateResID(tp, atmC, 0, tp->nres);
              resD = LocateResID(tp, atmD, 0, tp->nres);
              printf("MakeTermKey >> Atoms in this dihedral:\n"
                     "MakeTermKey >>    %.4s %.4s\n"
                     "MakeTermKey >>    %.4s %.4s\n"
                     "MakeTermKey >>    %.4s %.4s\n"
                     "MakeTermKey >>    %.4s %.4s\n", &tp->ResNames[4*resA],
                     &tp->AtomNames[4*atmA], &tp->ResNames[4*resB],
                     &tp->AtomNames[4*atmB], &tp->ResNames[4*resC],
                     &tp->AtomNames[4*atmC], &tp->ResNames[4*resD],
                     &tp->AtomNames[4*atmD]);
              exit(1);
            }

            // The number of torsion terms in this
            // system must be incremented here    
            mp->conf[i].hmap.id[nterm].a = atmA;
            mp->conf[i].hmap.id[nterm].b = atmB;
            mp->conf[i].hmap.id[nterm].c = atmC;
            mp->conf[i].hmap.id[nterm].d = atmD;
            nterm++;
          }
        }
      }
    }

    // Keys for NMR operations and CMAPs are not part of the big loop above.
    // They are handled in separate functions.
    if (order == 6) {
      MapNMROperations(mp, i);
    }
    else if (order == 7) {
      MapCmapTerms(mp, i);
    }

    // Increment the flag to state that this system has not only been
    // assigned to its term key repositories, but that they have also
    // been filled out.
    havekey[i] = 2;

    // Search for other conformations
    found = mp->conf[i].GroupNum;
    for (j = i+1; j < mp->nconf; j++) {
      if (havekey[j] == 1 && mp->conf[j].GroupNum == found) {
        ReflectTermKey(mp, &mp->conf[i], &mp->conf[j], order);
        havekey[j] = 2;
      }
    }
  }
  
  // Free allocated memory
  free(havekey);
}

//-----------------------------------------------------------------------------
// SetRestraint: encapsulates decisions about how to set the restraint  
//               equations for a given parameter system.
//
// Arguments:                                                           
//   mp:         the set of systems and parameter guidelines            
//   order:      the order of the term (bond=2, angle=3, torsion=4)     
//   prmidx:     the index of the parameter in the bond, angle, or      
//               torsion arrays of mp                                   
//-----------------------------------------------------------------------------
static void SetRestraint(prmset *mp, int order, int prmidx)
{
  int i, compval;

  // Initialize things with default parameters
  if (order == 2) {
    mp->bonds[prmidx].rstwK = mp->grstB;
    mp->bonds[prmidx].rstwl0 = mp->grstBcpl * mp->grstB;
    mp->bonds[prmidx].targK = mp->bonds[prmidx].K;
    mp->bonds[prmidx].targl0 = mp->bonds[prmidx].l0;
  }
  else if (order == 3) {
    mp->angls[prmidx].rstwK = mp->grstA;
    mp->angls[prmidx].rstwTh0 = mp->grstAcpl * mp->grstA;
    mp->angls[prmidx].targK = mp->angls[prmidx].K;
    mp->angls[prmidx].targTh0 = mp->angls[prmidx].th0;
  }
  else if (order == 4) {
    mp->torsions[prmidx].rstw = mp->grstH;
    mp->torsions[prmidx].target = 0.0;
  }
  else if (order == 6) {
    mp->nmrops[prmidx].rstw = mp->grstR;
    mp->nmrops[prmidx].target = 1.0;
  }
  else if (order == 7) {
    mp->cmaps[prmidx].rstw = mp->grstM;
    mp->cmaps[prmidx].target = CreateDmat(mp->cmaps[prmidx].esrfFit.row,
                                          mp->cmaps[prmidx].esrfFit.row, 0);
  }
  
  // Check through lists of specific parameters with 
  // user-specified restraint weights to find a match
  // to the one in question                          
  if (order == 2) {
    for (i = 0; i < mp->nuserrstB; i++) {
      if (TypeCompare(mp->bonds[prmidx].atype, mp->bonds[prmidx].btype,
                      "    ", "    ", mp->userrstB[i].atype,
                      mp->userrstB[i].btype, "    ", "    ", 2, 0) == 2) {
        if (mp->userrstB[i].has_wK == 1) {
          mp->bonds[prmidx].rstwK = mp->userrstB[i].rstwK;
        }
        if (mp->userrstB[i].has_wl0 == 1) {
          mp->bonds[prmidx].rstwl0 = mp->userrstB[i].rstwl0;
        }
        if (mp->userrstB[i].has_tK == 1) {
          mp->bonds[prmidx].targK = mp->userrstB[i].targK;
        }
        if (mp->userrstB[i].has_tl0 == 1) {
          mp->bonds[prmidx].targl0 = mp->userrstB[i].targl0;
        }
      }
    }
  }
  else if (order == 3) {
    for (i = 0; i < mp->nuserrstA; i++) {
      if (TypeCompare(mp->angls[prmidx].atype, mp->angls[prmidx].btype,
                      mp->angls[prmidx].ctype, "    ", mp->userrstA[i].atype,
                      mp->userrstA[i].btype, mp->userrstA[i].ctype, "    ", 3,
                      0) == 2) {
        if (mp->userrstA[i].has_wK == 1) {
          mp->angls[prmidx].rstwK = mp->userrstA[i].rstwK;
        }
        if (mp->userrstA[i].has_wTh0 == 1) {
          mp->angls[prmidx].rstwTh0 = mp->userrstA[i].rstwTh0;
        }
        if (mp->userrstA[i].has_tK == 1) {
          mp->angls[prmidx].targK = mp->userrstA[i].targK;
        }
        if (mp->userrstA[i].has_tTh0 == 1) {
          mp->angls[prmidx].targTh0 = mp->userrstA[i].targTh0;
        }
      }
    }
  }
  else if (order == 4) {
    for (i = 0; i < mp->nuserrstH; i++) {
      if (TypeCompare(mp->torsions[prmidx].atype, mp->torsions[prmidx].btype,
                      mp->torsions[prmidx].ctype, mp->torsions[prmidx].dtype,
                      mp->userrstH[i].atype, mp->userrstH[i].btype,
                      mp->userrstH[i].ctype, mp->userrstH[i].dtype,
                      4, 0) == 2) {
        if (fabs(mp->userrstH[i].pn) < 1.0e-4 ||
            fabs(fabs(mp->userrstH[i].pn) -
                 fabs(mp->torsions[prmidx].pn)) < 1.0e-4) {
          mp->torsions[prmidx].rstw = mp->userrstH[i].rstw;
          mp->torsions[prmidx].target = mp->userrstH[i].target;
        }
      }
    }
  }
  else if (order == 6) {
    for (i = 0; i < mp->nuserrstR; i++) {
      compval = CompareOperations(&mp->userrstR[i],
                                  &mp->nmrops[prmidx], mp, 0);
      if (compval == 0 || compval == 2) {
        mp->nmrops[prmidx].rstw = mp->userrstR[i].rstw;
        mp->nmrops[prmidx].target = mp->userrstR[i].target;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// FindAdjustableTerms: this routine parses the lists of bond, angle, torsion,
//                      and scee / scnb terms to obtain a list of adjustable
//                      parameters.  Adjustable parameters are suggested by the
//                      user in the form of atom types or names of atoms making
//                      up each bonded term or subsystem containing a scaling
//                      factor.  By matching these types or names to bonds in
//                      the list of structures, the program identifies the
//                      contributors to a linear least-squares fit and, if
//                      necessary, iteratively refines that fit by adjustment
//                      of nonlinear contributions such as bond length or
//                      dihedral phase angle.
//
// Arguments:                                                           
//   mp:     the set of systems and parameter guidelines                
//-----------------------------------------------------------------------------
static void FindAdjustableTerms(prmset *mp)
{
  int i, j, found, ncol;
  char* msg;

  // The fitting matrix column counter
  ncol = 0;

  // All bonds, angles, and torsions default
  // to no place in the fitting matrix      
  mp->nbvar = 0;
  mp->navar = 0;
  mp->nhvar = 0;
  mp->nmvar = 0;
  for (i = 0; i < mp->nbond; i++) {
    mp->bonds[i].fitcolX = -1;
    mp->bonds[i].fitcolY = -1;
    mp->bonds[i].fitcolZ = -1;
    mp->bonds[i].fitcolW = -1;
  }
  for (i = 0; i < mp->nangl; i++) {
    mp->angls[i].fitcolX = -1;
    mp->angls[i].fitcolY = -1;
  }
  for (i = 0; i < mp->ntor; i++) {
    mp->torsions[i].fitcol = -1;
  }
  for (i = 0; i < mp->nops; i++) {
    mp->nmrops[i].fitcol = -1;
  }
  for (i = 0; i < mp->ncmaps; i++) {
    mp->cmaps[i].fitcoli = -1;
    mp->cmaps[i].fitcolf = -1;
  }

  // Allocate memory for printing warning messages
  msg = (char*)malloc(MAXLINE*sizeof(char));

  // If all bonds are selected, enumerate them
  if (mp->FitAllBonds == 1) {
    for (i = 0; i < mp->nbond; i++) {
      if (mp->bonds[i].ninst > 0) {
        SetRestraint(mp, 2, i);
        if (mp->FitBondEq == 1) {
          mp->bonds[i].lbasisX = mp->bonds[i].l0 - mp->lpost;
          mp->bonds[i].lbasisY = mp->bonds[i].l0 + mp->lpost;
          mp->bonds[i].fitcolX = ncol;
          ncol++;
          mp->bonds[i].fitcolY = ncol;
          ncol++;
          mp->nbvar += 2;
        }
        else {
          mp->bonds[i].fitcolX = ncol;
          ncol++;
          mp->nbvar += 1;
        }
        if (mp->FitAllAugs == 1 ||
            (mp->FitCurrentAugs == 1 && mp->bonds[i].isAug == 1)) {
          mp->bonds[i].fitcolZ = ncol;
          ncol++;
          mp->bonds[i].fitcolW = ncol;
          ncol++;
          mp->nbvar += 2;
        }
      }
    }
  }

  // Check adjustable bonds
  else {
    for (i = 0; i < mp->nbadj; i++) {
      found = 0;
      for (j = 0; j < mp->nbond; j++) {
        if (TypeCompare(mp->badj[i].atype, mp->badj[i].btype, "    ",
                        "    ", mp->bonds[j].atype, mp->bonds[j].btype, "    ",
                        "    ", 2, 1) == 2 && mp->bonds[j].ninst > 0) {
          found = 1;
          if (mp->bonds[j].fitcolX >= 0) {
            continue;
          }
          SetRestraint(mp, 2, j);
          if (mp->FitBondEq == 1) {
            mp->bonds[j].lbasisX = mp->bonds[j].l0 - mp->lpost;
            mp->bonds[j].lbasisY = mp->bonds[j].l0 + mp->lpost;
            mp->bonds[j].fitcolX = ncol;
            ncol++;
            mp->bonds[j].fitcolY = ncol;
            ncol++;
            mp->nbvar += 2;
          }
          else {
            mp->bonds[j].fitcolX = ncol;
            ncol++;
            mp->nbvar += 1;
          }
          if (mp->FitAllAugs == 1 ||
              (mp->FitCurrentAugs == 1 && mp->bonds[i].isAug == 1)) {
            mp->bonds[j].fitcolZ = ncol;
            ncol++;
            mp->bonds[j].fitcolW = ncol;
            ncol++;
            mp->nbvar += 2;
          }
        }
      }
      if (found == 0) {
        sprintf(msg, "Bond type %.4s %.4s marked for optimization but not "
                "found in any systems.", mp->badj[i].atype, mp->badj[i].btype);
        PrintParagraph(msg, 79, "mdgx >>", stdout);
      }
    }
  }

  // If all angles are selected, enumerate them
  if (mp->FitAllAngles == 1) {
    for (i = 0; i < mp->nangl; i++) {
      if (mp->angls[i].ninst > 0) {
        SetRestraint(mp, 3, i);
        if (mp->FitAnglEq == 1) {
          mp->angls[i].tbasisX = mp->angls[i].th0 - mp->thpost;
          mp->angls[i].tbasisY = mp->angls[i].th0 + mp->thpost;
          mp->angls[i].fitcolX = ncol;
          ncol++;
          mp->angls[i].fitcolY = ncol;
          ncol++;
          mp->navar += 2;
        }
        else {
          mp->angls[i].fitcolX = ncol;
          ncol++;
          mp->navar += 1;
        }
      }
    }
  }

  // Check adjustable angles
  else {
    for (i = 0; i < mp->naadj; i++) {
      found = 0;
      for (j = 0; j < mp->nangl; j++) {
        if (TypeCompare(mp->aadj[i].atype, mp->aadj[i].btype,
                        mp->aadj[i].ctype, "    ", mp->angls[j].atype,
                        mp->angls[j].btype, mp->angls[j].ctype, "    ",
                        3, 1) == 2 && mp->angls[j].ninst > 0) {
          found = 1;
          if (mp->angls[j].fitcolX >= 0) {
            continue;
          }
          SetRestraint(mp, 3, j);
          if (mp->FitAnglEq == 1) {
            mp->angls[j].tbasisX = mp->angls[j].th0 - mp->thpost;
            mp->angls[j].tbasisY = mp->angls[j].th0 + mp->thpost;
            mp->angls[j].fitcolX = ncol;
            ncol++;
            mp->angls[j].fitcolY = ncol;
            ncol++;
            mp->navar += 2;
          }
          else {
            mp->angls[j].fitcolX = ncol;
            ncol++;
            mp->navar += 1;
          }
        }
      }
      if (found == 0) {
        sprintf(msg, "Angle type %.4s %.4s %.4s marked for optimization but "
               "not found in any systems.", mp->aadj[i].atype,
                mp->aadj[i].btype, mp->aadj[i].ctype);
        PrintParagraph(msg, 79, "mdgx >>", stdout);
      }
    }
  }

  // If all torsions are selected, enumerate them
  if (mp->FitAllTorsions == 1) {
    for (i = 0; i < mp->ntor; i++) {
      if (mp->torsions[i].ninst > 0 && mp->torsions[i].impr == 0) {
        SetRestraint(mp, 4, i);
        mp->torsions[i].fitcol = ncol;
        mp->nhvar += 1;
        ncol++;
      }
    }
  }

  // Check adjustable torsions
  else {
    for (i = 0; i < mp->nhadj; i++) {
      found = 0;
      for (j = 0; j < mp->ntor; j++) {
        if ((TypeCompare(mp->hadj[i].atype, mp->hadj[i].btype,
                         mp->hadj[i].ctype, mp->hadj[i].dtype,
                         mp->torsions[j].atype, mp->torsions[j].btype,
                         mp->torsions[j].ctype, mp->torsions[j].dtype,
                         4, 1) == 2 ||
             TypeCompare(mp->hadj[i].atype, mp->hadj[i].btype,
                         mp->hadj[i].ctype, mp->hadj[i].dtype,
                         mp->torsions[j].atype, mp->torsions[j].btype,
                         mp->torsions[j].ctype, mp->torsions[j].dtype,
                         5, 1) == 2) && mp->torsions[j].ninst > 0) {
          found = 1;
          if (mp->torsions[j].fitcol >= 0) {
            continue;
          }
          SetRestraint(mp, 4, j);
          mp->torsions[j].fitcol = ncol;
          mp->nhvar += 1;
          ncol++;
        }
      }
      if (found == 0) {
        sprintf(msg, "Torsion type %.4s %.4s %.4s %.4s marked for "
                "optimization but not found in any systems.\n",
                mp->hadj[i].atype, mp->hadj[i].btype, mp->hadj[i].ctype,
                mp->hadj[i].dtype);
        PrintParagraph(msg, 79, "mdgx >>", stdout);
      }
    }
  }

  // If all CMAP surfaces are optimizable, enumerate them
  if (mp->FitAllCMAPs == 1) {
    for (i = 0; i < mp->ncmaps; i++) {
      if (mp->cmaps[i].ninst > 0) {
        SetRestraint(mp, 7, i);
        mp->cmaps[i].fitcoli = ncol;
        mp->nmvar += 1;
        ncol += mp->cmaps[i].esrfFit.row * mp->cmaps[i].esrfFit.row;
        mp->cmaps[i].fitcolf = ncol;
      }
    }
  }

  // Check adjustable CMAP surfaces
  else {
    for (i = 0; i < mp->nmadj; i++) {
      found = 0;
      for (j = 0; j < mp->ncmaps; j++) {
        if (CmapResidueOverlap(&mp->madj[i], &mp->cmaps[j]) == 1) {
          found = 1;
          if (mp->cmaps[j].fitcoli >= 0) {
            continue;
          }
          SetRestraint(mp, 7, j);
        }
      }
    }
  }
  
  // If all NMR operations are optimizable, enumerate them
  if (mp->FitAllNMROps == 1) {
    for (i = 0; i < mp->nops; i++) {
      if (mp->nmrops[i].ninst > 0) {
        SetRestraint(mp, 6, i);
        mp->nmrops[i].fitcol = ncol;
        mp->nrvar += 1;
        ncol++;
      }
    }
  }

  // Check adjustable NMR operations
  else {
    for (i = 0; i < mp->nradj; i++) {

      // If there is a label specified, try to match that
      found = 0;
      for (j = 0; j < mp->nops; j++) {
        if (mp->radj[i].label[0] != '\0' &&
            strcmp(mp->radj[i].label, mp->nmrops[i].label) == 0) {
          found = 1;
        }
      }

      // If there are atom masks specified, try to match those
      // by comparison across all systems
      for (j = 0; j < mp->nops; j++) {

      }
      if (found == 0) {
        sprintf(msg, "NMR operation %s marked for optimization but not "
                "found.", mp->radj[i].label);
        PrintParagraph(msg, 79, "mdgx >>", stdout);
      }
    }
  }

  // Store the number of adjustable parameters
  mp->nparm = ncol;

  // Check whether 1:4 scaling factors are in play
  mp->nparm += mp->fitscnb + mp->fitscee;

  // Free allocated memory
  free(msg);
}

//-----------------------------------------------------------------------------
// CheckGeometryRestraints: check that the sums of various angles are valid
//                          restraints, containing atom types and angles
//                          actually present in the training set.
//                                                                      
// Arguments:                                                           
//   mp:     the set of systems and parameter guidelines                
//-----------------------------------------------------------------------------
static void CheckGeometryRestraints(prmset *mp)
{
  int i, j, k, m, found, nadj, engage;
  int* keep;
  cmat *gtype;
  geomrst tmpsum;

  // Bail out if no such restraints
  if (mp->nanglsum == 0) {
    return;
  }

  // Examine each restraint.  Verify that its angles  
  // are present and that at least two are adjustable.     
  keep = (int*)malloc(mp->nanglsum*sizeof(int));
  for (i = 0; i < mp->nanglsum; i++) {
    engage = 1;
    nadj = 0;
    gtype = &mp->anglsum[i].atmtype;
    for (j = 0; j < mp->anglsum[i].nvar; j++) {
      found = 0;
      for (k = 0; k < mp->nangl; k++) {
        if (mp->angls[k].ninst > 0 &&
            TypeCompare(gtype->map[3*j], gtype->map[3*j+1], gtype->map[3*j+2],
                        "    ", mp->angls[k].atype, mp->angls[k].btype,
                        mp->angls[k].ctype, "    ", 3, 1) == 2) {
          found = 1;
          if (mp->angls[k].fitcolY >= 0) {
            mp->anglsum[i].fitcol[2*j] = mp->angls[k].fitcolX;
            mp->anglsum[i].fitcol[2*j+1] = mp->angls[k].fitcolY;
            mp->anglsum[i].basis[2*j] = mp->angls[k].tbasisX;
            mp->anglsum[i].basis[2*j+1] = mp->angls[k].tbasisY;
            nadj++;
          }
          else {
            mp->anglsum[i].target -= mp->angls[k].th0;
          }
        }
      }
      if (found == 0) {
        engage = 0;
      }
    }
    if (nadj < 2) {
      engage = 0;
    }

    // Warn the user if this contraint is not applicable
    if (engage == 0 && mp->verbose == 1) {
      printf("mdgx >> Warning.  Angle sum restraint cannot be engaged due to "
             "some angles not\nmdgx >> being present, or not enough angles in "
             "the group being adjustable:\n");
      for (j = 0; j < mp->anglsum[i].nvar; j++) {
        printf("mdgx >>   %4.4s %4.4s %4.4s  ",
               mp->anglsum[i].atmtype.map[3*j],
               mp->anglsum[i].atmtype.map[3*j+1],
               mp->anglsum[i].atmtype.map[3*j+2]);
        if (mp->anglsum[i].fitcol[2*j] > 0) {
          printf("-> Adjustable\n");
        }
        else {
          printf("\n");
        }
      }
    }

    // Tidy up this constraint so that it
    // contains only adjustable angles   
    if (engage == 1) {
      for (j = 0; j < mp->anglsum[i].nvar; j++) {
        if (mp->anglsum[i].fitcol[2*j] == -1) {
          for (k = j; k < mp->anglsum[i].nvar-1; k++) {
            mp->anglsum[i].fitcol[2*k] = mp->anglsum[i].fitcol[2*(k+1)];
            mp->anglsum[i].fitcol[2*k+1] = mp->anglsum[i].fitcol[2*(k+1)+1];
            mp->anglsum[i].basis[2*k] = mp->anglsum[i].basis[2*(k+1)];
            mp->anglsum[i].basis[2*k+1] = mp->anglsum[i].basis[2*(k+1)+1];
            strncpy(mp->anglsum[i].atmtype.map[k],
                    mp->anglsum[i].atmtype.map[k+1], 4);
          }
          j--;
          mp->anglsum[i].nvar--;
        }
      }
    }
    keep[i] = engage;
  }
  j = 0;
  for (i = 0; i < mp->nanglsum; i++) {
    if (keep[i] == 1) {
      if (j != i) {
        SWAP(mp->anglsum[j], mp->anglsum[i], tmpsum);
      }
      j++;
    }
  }
  for (i = j; i < mp->nanglsum; i++) {
    DestroyCmat(&mp->anglsum[i].atmtype);
    free(mp->anglsum[i].fitcol);
    free(mp->anglsum[i].rstrow);
    free(mp->anglsum[i].basis);
  }
  mp->nanglsum = ISum(keep, mp->nanglsum);

  // Free allocated memory
  free(keep);
}

//-----------------------------------------------------------------------------
// ComputeDistance: compute the distance between points A and B.
//
// Arguments:
//   crd:           the coordinate set to work from
//   [a,b]idx:      indices of the atoms A, B, and C in the coordinate set
//-----------------------------------------------------------------------------
static double ComputeDistance(coord *crd, int aidx, int bidx)
{
  int i;
  double r, dx, dy, dz;
  double *aptr, *bptr;

  aptr = &crd->loc[3*aidx];
  bptr = &crd->loc[3*bidx];
  dx = bptr[0] - aptr[0];
  dy = bptr[1] - aptr[1];
  dz = bptr[2] - aptr[2];
  r = sqrt(dx*dx + dy*dy + dz*dz);

  return r;
}

//-----------------------------------------------------------------------------
// ComputeAngle: compute the value of an angle A-B-C.  This is encapsulated
//               here, but not made more widely available, because only in
//               these fitting routines are the individual displacements no
//               longer necessary.
//
// Arguments:
//   crd:           the coordinate set to work from
//   [a,b,c]idx:    indices of the atoms A, B, and C in the coordinate set
//-----------------------------------------------------------------------------
static double ComputeAngle(coord *crd, int aidx, int bidx, int cidx)
{
  int i;
  double mgba, mgbc, invbabc, costheta, theta;
  double ba[3], bc[3];
  double *aptr, *bptr, *cptr;

  aptr = &crd->loc[3*aidx];
  bptr = &crd->loc[3*bidx];
  cptr = &crd->loc[3*cidx];
  for (i = 0; i < 3; i++) {
    ba[i] = aptr[i] - bptr[i];
    bc[i] = cptr[i] - bptr[i];
  }
  mgba = ba[0]*ba[0] + ba[1]*ba[1] + ba[2]*ba[2];
  mgbc = bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2];
  invbabc = 1.0/sqrt(mgba*mgbc);
  costheta = (ba[0]*bc[0] + ba[1]*bc[1] + ba[2]*bc[2]) * invbabc;
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  theta = acos(costheta);

  return theta;
}

//-----------------------------------------------------------------------------
// ComputeDiehdral: compute the value of a dihedral angle A-B-C-D.  This is
//                  encapsulated here for the same reason as ComputeAngle
//                  above.
//
// Arguments:
//   crd:           the coordinate set to work from
//   [a,b,c,d]idx:  indices of the atoms A, B, C, and D in the coordinate set
//-----------------------------------------------------------------------------
static double ComputeDihedral(coord *crd, int aidx, int bidx, int cidx,
                              int didx)
{
  int i;
  double costheta, theta;
  double ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  double *aptr, *bptr, *cptr, *dptr;

  aptr = &crd->loc[3*aidx];
  bptr = &crd->loc[3*bidx];
  cptr = &crd->loc[3*cidx];
  dptr = &crd->loc[3*didx];
  for (i = 0; i < 3; i++) {
    ab[i] = bptr[i] - aptr[i];
    bc[i] = cptr[i] - bptr[i];
    cd[i] = dptr[i] - cptr[i];
  }
  CrossP(ab, bc, crabbc);
  CrossP(bc, cd, crbccd);
  costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  costheta /=
    sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
         (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  CrossP(crabbc, crbccd, scr);
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  if (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) {
    theta = acos(costheta);
  }
  else {
    theta = -acos(costheta);
  }

  return theta;
}

//-----------------------------------------------------------------------------
// ComputeCMAPPosition: compute the two values of a CMAP, that is dihedrals
//                      A-B-C-D and B-C-D-E.  This draws on ComputeDihedral
//                      above.
//
// Arguments:
//   crd:             the coordinate set to work from
//   [a,b,c,d,e]idx:  indices of the atoms A through E in the coordinate set
//   val[1,2]:        first and second angle values of the CMAP (returned)
//-----------------------------------------------------------------------------
static void ComputeCMAPPosition(coord *crd, int aidx, int bidx, int cidx,
                                int didx, int eidx, double *val1, double *val2)
{
  *val1 = ComputeDihedral(crd, aidx, bidx, cidx, didx) + PI;
  *val2 = ComputeDihedral(crd, bidx, cidx, didx, eidx) + PI;
  if (*val1 < 0.0 || *val1 >= 2.0 * PI) {
    *val1 -= floor(*val1 / (2.0 * PI)) * 2.0 * PI;
  }
  if (*val2 < 0.0 || *val2 >= 2.0 * PI) {
    *val2 -= floor(*val2 / (2.0 * PI)) * 2.0 * PI;
  }
}

//-----------------------------------------------------------------------------
// EvaluateCMAPEnergy: evaluate the energy of a CMAP given a surface to work
//                     from.
//
// Arguments:
//   val[1,2]:  first and second values of the CMAP position (values in
//              radians--they will be converted to degrees when needed)
//   cmsurf:    the CMAP energy surface
//   cmdphi:    the CMAP derivative with respect to its first angle dimension 
//   cmdpsi:    the CMAP derivative with respect to its second angle dimension 
//   cmd2pp:    the CMAP mixed partial derivative
//-----------------------------------------------------------------------------
double EvaluateCMAPEnergy(double val1, double val2, dmat *cmsurf,
                          dmat *cmdphi, dmat *cmdpsi, dmat *cmd2pp)
{
  int binx, biny, binxp, binyp, nbin;
  double bval1, bval2, dx, dy;
  double nrg00, nrg01, nrg10, nrg11, dphi00, dphi01, dphi10, dphi11;
  double dpsi00, dpsi01, dpsi10, dpsi11, d2pp00, d2pp01, d2pp10, d2pp11;
  double a00, a10, a20, a30, a01, a11, a21, a31, a02, a12, a22, a32;
  double a03, a13, a23, a33, Ecmap;
  
  // Get the bin and displacements within the bin
  nbin = cmsurf->row;
  bval1 = (double)nbin * (val1 / (2.0 * PI));
  bval2 = (double)nbin * (val2 / (2.0 * PI));
  binx = (int)bval1;
  biny = (int)bval2;
  dx = bval1 - (double)binx;
  dy = bval2 - (double)biny;
  binx += nbin * (binx < 0);
  biny += nbin * (biny < 0);
  binx -= nbin * (binx >= nbin);
  biny -= nbin * (biny >= nbin);
  binxp = binx + 1;
  binyp = biny + 1;
  binxp -= (binxp >= nbin) * nbin;
  binyp -= (binyp >= nbin) * nbin;
 
  // Read off grid energies and derivatives
  nrg00  = cmsurf->map[binx ][biny ];
  nrg01  = cmsurf->map[binx ][binyp];
  nrg10  = cmsurf->map[binxp][biny ];
  nrg11  = cmsurf->map[binxp][binyp];
  dphi00 = cmdphi->map[binx ][biny ];
  dphi01 = cmdphi->map[binx ][binyp];
  dphi10 = cmdphi->map[binxp][biny ];
  dphi11 = cmdphi->map[binxp][binyp];
  dpsi00 = cmdpsi->map[binx ][biny ];
  dpsi01 = cmdpsi->map[binx ][binyp];
  dpsi10 = cmdpsi->map[binxp][biny ];
  dpsi11 = cmdpsi->map[binxp][binyp];
  d2pp00 = cmd2pp->map[binx ][biny ];
  d2pp01 = cmd2pp->map[binx ][binyp];
  d2pp10 = cmd2pp->map[binxp][biny ];
  d2pp11 = cmd2pp->map[binxp][binyp];
  
  // Evaluate the energy at the point
  a00        =                    nrg00;
  a10        =                   dphi00;
  a20        = (-3.0 * nrg00 ) + ( 3.0 * nrg10) -
               ( 2.0 * dphi00) -                   dphi10;
  a30        = ( 2.0 * nrg00 ) - ( 2.0 * nrg10) +
                       dphi00  +                   dphi10;
  a01        =                   dpsi00;
  a11        =                   d2pp00;
  a21        = (-3.0 * dpsi00) + ( 3.0 * dpsi10) -
               ( 2.0 * d2pp00) -                   d2pp10;
  a31        = ( 2.0 * dpsi00) - ( 2.0 * dpsi10) +
                       d2pp00  +                   d2pp10;
  a02        = (-3.0 * nrg00 ) + ( 3.0 * nrg01 ) -
               ( 2.0 * dpsi00) -                   dpsi01;
  a12        = (-3.0 * dphi00) + ( 3.0 * dphi01) -
               ( 2.0 * d2pp00) -                   d2pp01;
  a22        = ( 9.0 * (nrg00  - nrg10  - nrg01  + nrg11 ) ) +
               ( 6.0 * dphi00) + ( 3.0 * dphi10) -
               ( 6.0 * dphi01) - ( 3.0 * dphi11) +
               ( 6.0 * dpsi00) - ( 6.0 * dpsi10) +
               ( 3.0 * dpsi01) - ( 3.0 * dpsi11) +
               ( 4.0 * d2pp00) + ( 2.0 * d2pp10) +
               ( 2.0 * d2pp01) +                   d2pp11;
  a32        = (-6.0 * (nrg00  - nrg10  - nrg01  + nrg11 ) ) +
               (-3.0 * (dphi00 + dphi10 - dphi01 - dphi11) ) +
               (-4.0 * dpsi00) + ( 4.0 * dpsi10) -
                       d2pp01  -                   d2pp11  +
               (-2.0 * (d2pp00 + d2pp10 + dpsi01 - dpsi11) );
  a03        = ( 2.0 * nrg00 ) - ( 2.0 * nrg01 ) +
                       dpsi00  +                   dpsi01;
  a13        = ( 2.0 * dphi00) - ( 2.0 * dphi01) +
                       d2pp00  +                   d2pp01;
  a23        = (-6.0 * (nrg00  - nrg10  - nrg01  + nrg11 ) ) +
               (-2.0 * (d2pp00 + dphi10 + d2pp01 - dphi11) ) +
               (-3.0 * (dpsi00 - dpsi10 + dpsi01 - dpsi11) ) +
               (-4.0 * dphi00) - d2pp10 +
               ( 4.0 * dphi01) - d2pp11;
  a33        = ( 4.0 * (nrg00  - nrg10  - nrg01  + nrg11 ) ) +
               ( 2.0 * (dphi00 + dphi10 - dphi01 - dphi11) ) +
               ( 2.0 * (dpsi00 - dpsi10 + dpsi01 - dpsi11) ) +
                        d2pp00 +                   d2pp10  +
                        d2pp01 +                   d2pp11;

  // Evaluate the final energy
  Ecmap = ((a33*dy + a32)*dy + a31)*dy + a30;
  Ecmap = ((a23*dy + a22)*dy + a21)*dy + a20 + dx*Ecmap;
  Ecmap = ((a13*dy + a12)*dy + a11)*dy + a10 + dx*Ecmap;
  Ecmap = ((a03*dy + a02)*dy + a01)*dy + a00 + dx*Ecmap;

  return Ecmap;
}

//-----------------------------------------------------------------------------
// StencilScalingFactor: encapsulates the tedious computations of stencil
//                       contributions to the influence of CMAP grid points on
//                       the final energy.  This takes advantage of a concise
//                       statement of the bicubic spline interpolation formula
//                       that uses only 4x4 matrices:
//
//                                    [ a00 a01 a02 a03 ]   [  1   ]
//       Ecmap = [ 1 dx dx^2 dx^3 ] * [ a10 a11 a12 a13 ] * [  dy  ]
//                                    [ a20 a21 a22 a23 ]   [ dy^2 ]
//                                    [ a30 a31 a32 a33 ]   [ dy^3 ]
//
// where the values a00, a01, etc. are familiar from the code.  The matrix of
// values a00, a01, ..., a33 may be obtained from a compact 4x4 matrix
// multiplication:
//
//       [ a00 a01 a02 a03 ]         [  nrg00  nrg01 dpsi00 dpsi01 ]
//       [ a10 a11 a12 a13 ] = C0' * [  nrg10  nrg11 dpsi10 dpsi11 ] * C0
//       [ a20 a21 a22 a23 ]         [ dphi00 dphi01 d2pp00 dp2201 ]
//       [ a30 a31 a32 a33 ]         [ dphi10 dphi11 d2pp10 dp2211 ]
//
// where the matrix C0 is [ 1 0 0 0 ; 0 0 1 0 ; -3 3 -2 -1 ; 2 -2 1 1 ].
// This form can drastically reduce the number of computations required to
// accumulate the influence matrix.
//
// Arguments:
//   dx, dy:      fractional coordinates within the grid bin for the arbitrary
//                position on the CMAP
//   isx, isy:    the matrix element of the [ nrg00, ..., d2pp11 ] matrix to
//                fold into the result
//   S, V:        scratch matrices
//   C0, C0t:     right and left matrices in the equation above
//-----------------------------------------------------------------------------
static double StencilScalingFactor(double dx, double dy, int isx, int isy,
                                   dmat *S, dmat *V, dmat *C0, dmat *C0t)
{
  int i, j;
  double factor, dxv, dyv;

  // Set the component of the S (ngr00, ..., d2pp11) matrix
  SetDVec(S->data, 16, 0.0);
  S->map[isx][isy] = 1.0;

  // Get the matrix [ a00, ..., a33 ] for this energy / derivative
  DMatMult(*S, *C0, *V);
  DMatMult(*C0t, *V, *S);

  // Accumulate the answer
  factor = 0.0;
  dxv = 1.0;
  for (i = 0; i < 4; i++) {
    dyv = 1.0;
    for (j = 0; j < 4; j++) {
      factor += S->map[i][j] * dxv * dyv;
      dyv *= dy;
    }
    dxv *= dx;
  }

  return factor;
}

//-----------------------------------------------------------------------------
// ApplyStencil: apply the stencil to obtain the energy weight matrix of a
//               CMAP's first partial derivative or a mixed partial derivative
//               for the grid point of interest.
//
// Arguments:
//   s:      the stencil to apply
//   M:      the matrix to fill out (this will be zero'ed prior to applying the
//           stencil)
//   xpt:    the column to start on
//   ypt:    the row to start on
//   dir:    direction (0 for x/phi, 1 for y/psi)
//   order:  0 for the undifferentiated function, 1 for first partial
//           derivative, 2 for mixed partial derivative
//-----------------------------------------------------------------------------
void ApplyStencil(dmat *s, dmat *M, int xpt, int ypt, int dir, int order)
{
  int i, j, ip, jp, npts;

  npts = s->col;
  SetDVec(M->data, M->row * M->col, 0.0);
  if (order == 0) {
    M->map[xpt][ypt] = 1.0;
  }
  else if (dir == 0 && order == 1) {

    // Partial derivative in phi
    for (i = 0; i < npts; i++) {
      ip = i + xpt;
      ip -= (ip >= npts) * npts;
      M->map[ip][ypt] = s->data[i];
    }
  }
  else if (dir == 1 && order == 1) {

    // Partial derivative in psi
    for (i = 0; i < npts; i++) {
      ip = i + ypt;
      ip -= (ip >= npts) * npts;
      M->map[xpt][ip] = s->data[i];
    }
  }
  else {

    // Anything else gets the mixed partial derivative
    for (i = 0; i < npts; i++) {
      ip = i + xpt;
      ip -= (ip >= npts) * npts;
      for (j = 0 ; j < npts; j++) {
        jp = j + ypt;
        jp -= (jp >= npts) * npts;
        M->map[ip][jp] = s->data[i] * s->data[j];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// GetCMAPFootprint: given a pair of angle coordinates (phi and psi), produce
//                   the contributions (scaling factors on the values of the
//                   CMAP surface) for all CMAP grid points.  This function
//                   computes how the values of the CMAP derivatives at four
//                   grid points surrounding the coordinates in question
//                   depend on contributions from the rest of the grid, then
//                   looks at how the CMAP values and derivatives at those
//                   four grid points contribute to the final result.  A
//                   central relationship is expressed by the code above, or
//                   in long form:
//
// Ecmap =
//         a03*dy*dy*dy +          a02*dy*dy +          a01*dy +          a00 +
//      dx*a13*dy*dy*dy +       dx*a12*dy*dy +       dx*a11*dy +       dx*a10 +
//   dx*dx*a23*dy*dy*dy +    dx*dx*a22*dy*dy +    dx*dx*a21*dy +    dx*dx*a20 +
//dx*dx*dx*a33*dy*dy*dy + dx*dx*dx*a32*dy*dy + dx*dx*dx*a31*dy + dx*dx*dx*a30
//
//                   While this is expensive to compute, it saves a great deal
//                   of memory and only needs to be done when composing the
//                   fitting matrix or evaluating the solution.
//
// Arguments:
//   val1:     the first angle coordinate
//   val2:     the second angle coordinate
//   stencil:  dependencies of the first derivative of the CMAP surface on
//             neighboring grid points.  This applies in either direction,
//             and since the map is square the outer product of the stencil
//             with itself gives the weights for the mixed partial derivative
//             at some grid point.  Also contains the CMAP dimension.
//   tprof:    scratch space for generating the contributions matrix
//   result:   the proportionate contributions of all grid points for this
//             phi, psi coordinate (returned)
//-----------------------------------------------------------------------------
static void GetCMAPFootPrint(double val1, double val2, dmat *stencil,
                             dmat *tprof, dmat *result)
{
  int i, j, nbin, binx, biny, binxp, binyp;
  double bval1, bval2, dx, dy;
  dmat C0, C0t, S, V;
  
  // Compute the bin and displacement  
  nbin = stencil->col;
  bval1 = (double)nbin * (val1 / (2.0 * PI));
  bval2 = (double)nbin * (val2 / (2.0 * PI));
  binx = (int)bval1;
  biny = (int)bval2;
  dx = bval1 - (double)binx;
  dy = bval2 - (double)biny;
  binx += nbin * (binx < 0);
  biny += nbin * (biny < 0);
  binx -= nbin * (binx >= nbin);
  biny -= nbin * (biny >= nbin);
  binxp = binx + 1;
  binyp = biny + 1;
  binxp -= (binxp >= nbin) * nbin;
  binyp -= (binyp >= nbin) * nbin;

  // Clear the result
  SetDVec(result->data, result->row * result->col, 0.0);

  // Set the matrices C0 and C0t (see the StencilScalingFactor function
  // above).
  S   = CreateDmat(4, 4, 0);
  V   = CreateDmat(4, 4, 0);
  C0  = CreateDmat(4, 4, 0);
  C0t = CreateDmat(4, 4, 0);
  C0.map[0][0] =  1.0;
  C0.map[0][2] = -3.0;
  C0.map[0][3] =  2.0;
  C0.map[1][2] =  3.0;
  C0.map[1][3] = -2.0;
  C0.map[2][1] =  1.0;
  C0.map[2][2] = -2.0;
  C0.map[2][3] =  1.0;
  C0.map[3][2] = -1.0;
  C0.map[3][3] =  1.0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      C0t.map[i][j] = C0.map[j][i];
    }
  }

  // Add the raw energy contributions from the surrounding points
  result->map[binx ][biny ] = StencilScalingFactor(dx, dy, 0, 0,
                                                   &S, &V, &C0, &C0t);
  result->map[binxp][biny ] = StencilScalingFactor(dx, dy, 1, 0,
                                                   &S, &V, &C0, &C0t);
  result->map[binx ][binyp] = StencilScalingFactor(dx, dy, 0, 1,
                                                   &S, &V, &C0, &C0t);
  result->map[binxp][binyp] = StencilScalingFactor(dx, dy, 1, 1,
                                                   &S, &V, &C0, &C0t);
  
  // Compute the dphi profiles and add multiples of them to the result.
  ApplyStencil(stencil, tprof,  binx,  biny, 0, 1);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 2, 0,
                                                   &S, &V, &C0, &C0t));
  ApplyStencil(stencil, tprof,  binx, binyp, 0, 1);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 2, 1,
                                                   &S, &V, &C0, &C0t));
  ApplyStencil(stencil, tprof, binxp,  biny, 0, 1);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 3, 0,
                                                   &S, &V, &C0, &C0t));
  ApplyStencil(stencil, tprof, binxp, binyp, 0, 1);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 3, 1,
                                                   &S, &V, &C0, &C0t));

  // Compute the dpsi profiles and add multiples of them to the result.
  ApplyStencil(stencil, tprof,  binx,  biny, 1, 1);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 0, 2,
                                                   &S, &V, &C0, &C0t));
  ApplyStencil(stencil, tprof,  binx, binyp, 1, 1);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 0, 3,
                                                   &S, &V, &C0, &C0t));
  ApplyStencil(stencil, tprof, binxp,  biny, 1, 1);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 1, 2,
                                                   &S, &V, &C0, &C0t));
  ApplyStencil(stencil, tprof, binxp, binyp, 1, 1);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 1, 3,
                                                   &S, &V, &C0, &C0t));
  
  // Compute the d2pp profiles and add multiples of them to the result.
  ApplyStencil(stencil, tprof,  binx,  biny, 0, 2);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 2, 2,
                                                   &S, &V, &C0, &C0t));
  ApplyStencil(stencil, tprof,  binx, binyp, 0, 2);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 2, 3,
                                                   &S, &V, &C0, &C0t));
  ApplyStencil(stencil, tprof, binxp,  biny, 0, 2);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 3, 2,
                                                   &S, &V, &C0, &C0t));
  ApplyStencil(stencil, tprof, binxp, binyp, 0, 2);
  AddDmat(result, tprof, 1.0, StencilScalingFactor(dx, dy, 3, 3,
                                                   &S, &V, &C0, &C0t));
}

//-----------------------------------------------------------------------------
// EvalNMROperation: function to encapsulate the evaluation of NMR operations,
//                   so that this meticulous code need not be recoded many
//                   times.
//
// Arguments:
//   myop:     the operation of interest
//   r:        the coordinate along the operation
//-----------------------------------------------------------------------------
static double EvalNMROperation(nmroper *myop, double r)
{
  int i;
  double dl, val;

  for (i = 0; i < myop->nlayer; i++) {
    if (r < myop->r1[i]) {

      // Linear where it will smoothly transition into rk2 harmonic
      dl = (myop->r2[i] - myop->r1[i]);
      val = (myop->rk2[i])*dl*dl + (2.0*myop->rk2[i])*dl*(myop->r1[i] - r);
    }
    else if (r < myop->r2[i]) {

      // Parabolic, heading towards zero
      dl = (myop->r2[i] - r);
      val = (myop->rk2[i])*dl*dl;
    }
    else if (r < myop->r3[i]) {
      val = 0.0;
    }
    else if (r < myop->r4[i]) {

      // Parabolic, heading towards some asymptote
      dl = (myop->r4[i] - r);
      val = (myop->rk3[i])*dl*dl;
    }
    else {

      // Linear with a slope where rk3 harmonic left off
      dl = (myop->r4[i] - myop->r3[i]);
      val = (myop->rk3[i])*dl*dl + (2.0*myop->rk3[i])*dl*(r - myop->r4[i]);
    }
  }

  return val;
}

//-----------------------------------------------------------------------------
// BondTermParse: compute the contributions due to a bonded term, parsed into
//                the kernel and coefficient-scaled values.
//
// Arguments:                                                           
//   conf:   the molecular mechanics system, complete with topology,    
//           coordinates, and energy tables                             
//   idx:    the index of the bonded term in the conformation
//   order:  the order of the bonded term (2 = bond, 3 = angle...)      
//   mp:     the parameter set, containing master tables of terms       
//   val:    returned value of the bonded term (i.e. bond length, dihedral
//           angle)
//   krnlX:  the "X" kernel of the bond contribution (returned) (for    
//           torsions, this is the only contribution--for special cases 
//           of bonds and angles there may also be a "Y" contribution)  
//   krnlY:  the "Y" kernel of the bond contribution (returned)         
//   krnlZ:  the "Z" kernel of the bond contribution, describing the pulling
//           augmentation for bonds (returned)
//   krnlW:  the "W" kernel of the bond contribution, describing the
//           compression augmentation for bonds (returned)
//   cntrb:  the contribution, krnl*stiffness (returned)                
//-----------------------------------------------------------------------------
static int BondTermParse(mmsys *conf, int idx, int order, prmset *mp,
                         double *val, double *krnlX, double *krnlY,
                         double *krnlW, double *krnlZ, double *cntrb)
{
  int natom;
  double r, dl, sangle, theta;
  double dtheta, dlX, dlY, dlZ, dlW, dthetaX, dthetaY;
  xbonddef *bndtmp;
  xangldef *angtmp;
  torterm *tortmp;
  nmroper *myop;

  // Bond
  if (order == 2) {
    r = ComputeDistance(&conf->crd, conf->bmap.id[idx].a,
                        conf->bmap.id[idx].b);
    bndtmp = &mp->bonds[conf->bmap.id[idx].key];
    dl = r - bndtmp->l0;
    dlZ = (r > bndtmp->lpull0)  ? r - bndtmp->lpull0 : 0.0;
    dlW = (r < bndtmp->lpress0) ? r - bndtmp->lpress0 : 0.0;
    *val = dl;
    if (mp->FitBondEq == 1) {
      dlX = r - bndtmp->lbasisX;
      dlY = r - bndtmp->lbasisY;
      *krnlX = dlX * dlX;
      *krnlY = dlY * dlY;
      *krnlZ = dlZ * dlZ;
      *krnlW = dlW * dlW;
    }
    else {
      *krnlX = dl * dl;
      *krnlZ = dlZ * dlZ;
      *krnlW = dlW * dlW;
    }
    *cntrb = (bndtmp->K * dl * dl) + (bndtmp->Kpull * dlZ * dlZ) +
             (bndtmp->Kpress * dlW * dlW);
  }

  // Angle
  else if (order == 3) {
    theta = ComputeAngle(&conf->crd, conf->amap.id[idx].a,
                         conf->amap.id[idx].b, conf->amap.id[idx].c);
    angtmp = &mp->angls[conf->amap.id[idx].key];
    dtheta = theta - angtmp->th0;
    *val = dtheta;
    if (mp->FitAnglEq == 1) {
      dthetaX = theta - angtmp->tbasisX;
      dthetaY = theta - angtmp->tbasisY;
      *krnlX = dthetaX*dthetaX;
      *krnlY = dthetaY*dthetaY;
    }
    else {
      *krnlX = dtheta*dtheta;
    }
    *cntrb = angtmp->K*dtheta*dtheta;
  }

  // Torsion
  else if (order == 4) {
    theta = ComputeDihedral(&conf->crd, conf->hmap.id[idx].a,
                            conf->hmap.id[idx].b, conf->hmap.id[idx].c,
                            conf->hmap.id[idx].d);
    tortmp = &mp->torsions[conf->hmap.id[idx].key];
    sangle = tortmp->pn*theta - tortmp->phase;
    *val = theta;
    *krnlX = 1.0 + cos(sangle);
    *cntrb = tortmp->K * (1.0 + cos(sangle));
  }

  // CMAP
  else if (order == 7) {
    ComputeCMAPPosition(&conf->crd, conf->mmap.id[idx].a,
                        conf->mmap.id[idx].b, conf->mmap.id[idx].c,
                        conf->mmap.id[idx].d, conf->mmap.id[idx].e,
                        &val[0], &val[1]);
    *cntrb = EvaluateCMAPEnergy(val[0], val[1],
                                &mp->cmaps[conf->mmap.id[idx].key].esrf,
                                &mp->cmaps[conf->mmap.id[idx].key].dphi,
                                &mp->cmaps[conf->mmap.id[idx].key].dpsi,
                                &mp->cmaps[conf->mmap.id[idx].key].d2pp);
  }
  
  // NMR operation
  else if (order == 6) {

    // Compute the coordinate--whether it is distance, angle, or dihedral
    myop = &mp->nmrops[conf->rmap.id[idx].key];
    natom = myop->order;
    if (natom == 2) {
      r = ComputeDistance(&conf->crd, conf->rmap.id[idx].a,
                          conf->rmap.id[idx].b);
    }
    else if (natom == 3) {
      r = ComputeAngle(&conf->crd, conf->rmap.id[idx].a,
                       conf->rmap.id[idx].b, conf->rmap.id[idx].c);
    }
    else if (natom == 4) {
      r = ComputeDihedral(&conf->crd, conf->rmap.id[idx].a,
                          conf->rmap.id[idx].b, conf->rmap.id[idx].c,
                          conf->rmap.id[idx].d);
    }

    // Every NMR restraint now works the same, until you take a derivative,
    // which we don't (yet).
    *val = r;
    *krnlX = EvalNMROperation(myop, r);
    *cntrb = *krnlX;
  }

  // Alert if the contribution exceeds a specified tolerance
  return (*cntrb > mp->mmtol) ? 1 : 0;
}

//-----------------------------------------------------------------------------
// AllBondedTerms: compute the molecular mechanics contributions due to bonded
//                 terms.
//
// Arguments:                                                           
//   conf:   the molecular mechanics system, complete with topology,    
//           coordinates, and energy tables                             
//   mp:     the parameter set, containing master tables of bonded terms
//   strain: indicator that there may be a bad interaction in this      
//           configuration, grounds to explore alternative geometries   
//-----------------------------------------------------------------------------
static double AllBondedTerms(mmsys *conf, prmset *mp, int *strain)
{
  int i, scon;
  double mmnrg, tmpval;

  // Unpack the system; only the coordinates are needed, as all
  // topology information is now stored in the prmset struct   
  mmnrg = 0.0;

  // Flag for bad interactions, which may     
  // indicate a misplaced atom in the molecule
  scon = 0;

  // Bonds
  for (i = 0; i < conf->bmap.nbond; i++) {
    scon += BondTermParse(conf, i, 2, mp, &conf->bmap.val[i],
                          &conf->bmap.UkernelX[i], &conf->bmap.UkernelY[i],
                          &conf->bmap.UkernelZ[i], &conf->bmap.UkernelW[i],
                          &conf->bmap.Ucontrib[i]);
    mmnrg += conf->bmap.Ucontrib[i];
  }

  // Angles
  for (i = 0; i < conf->amap.nangl; i++) {
    scon += BondTermParse(conf, i, 3, mp, &conf->amap.val[i],
                          &conf->amap.UkernelX[i], &conf->amap.UkernelY[i],
                          &tmpval, &tmpval, &conf->amap.Ucontrib[i]);
    mmnrg += conf->amap.Ucontrib[i];
  }

  // Dihedral fourier terms
  for (i = 0; i < conf->hmap.ntterm; i++) {
    scon += BondTermParse(conf, i, 4, mp, &conf->hmap.val[i],
                          &conf->hmap.Ukernel[i], &tmpval, &tmpval, &tmpval,
                          &conf->hmap.Ucontrib[i]);
    mmnrg += conf->hmap.Ucontrib[i];
  }

  // CMAP terms
  for (i = 0; i < conf->mmap.ncmap; i++) {
    scon += BondTermParse(conf, i, 7, mp, &conf->mmap.val[2*i], &tmpval,
                          &tmpval, &tmpval, &tmpval, &conf->mmap.Ucontrib[i]);
  }

  // NMR operations
  for (i = 0; i < conf->rmap.nops; i++) {
    scon += BondTermParse(conf, i, 6, mp, &conf->rmap.val[i],
                          &conf->rmap.Ukernel[i], &tmpval, &tmpval, &tmpval,
                          &conf->rmap.Ucontrib[i]);
    mmnrg += conf->rmap.Ucontrib[i];
  }

  // Are there any strained interactions?
  *strain = scon;

  return mmnrg;
}

//-----------------------------------------------------------------------------
// PlaceAtoms: this function references a list of atoms to rearrange a molecule
//             A partial molecule is built, unless the list contains aliases
//             for all atoms.                          
//                                                                      
// Arguments:                                                           
//   currvec:   the current state vector of the system... currvec[k] contains
//              the identity of atmorder[k] in the real system 
//   atmorder:  the atom ordering vector... it is not simply 0, 1, 2, ... n
//              because the molecule must be filled out one atom at a time,
//              always with a means of checking whether the latest addition
//              has introduced some impossible situation
//   asmts:     contains a more intelligible list of which atom is in what
//              slot; asmts[k] is the number of the atom from the original
//              coordinate which is now being used as atom k according to the
//              topology.  This array is initialized to -1 at the beginning of
//              this routine, to carry additional information about which
//              atoms have NOT been placed.     
//   conf:      the molecular system, contains topology and coordinates 
//   ref:       the reference coordinates                               
//-----------------------------------------------------------------------------
static void PlaceAtoms(int* currvec, int* atmorder, int* asmts, int nassigned,
                       mmsys *conf, coord *ref)
{
  int i, j;

  SetIVec(asmts, conf->tp->natom, -1);
  for (i = 0; i < nassigned; i++) {
    for (j = 0; j < 3; j++) {
      conf->crd.loc[3*atmorder[i]+j] = ref->loc[3*currvec[i]+j];
    }
    asmts[atmorder[i]] = currvec[i];
  }
}

//-----------------------------------------------------------------------------
// IncrementState: increment the state of the guess held in currvec.  Returns 1
//                 if the highest possible state has been reached.
//
// Arguments:                                                           
//   currvec:   the current state vector of the system... currvec[k]    
//              contains the identity of atmorder[k] in the real system 
//-----------------------------------------------------------------------------
static int IncrementState(int* currvec, int *nassigned, int natom)
{
  int i, pivot, inci;

  // Increment the ith counter of the state vector
  // and then check backwards to verify that no   
  // duplicate assignments have been made.  Lower-
  // numbered counters take precedence if any     
  // duplicates exist.                            
  pivot = *nassigned-1;
  inci = 1;
  while (inci == 1 && pivot >= 0) {
    inci = 0;
    currvec[pivot] += 1;

    // If we have hit the roof, erase the most recently
    // placed atoms and back-track until counters can  
    // be incremented again.                           
    if (currvec[pivot] == natom) {
      currvec[pivot] = -1;
      inci = 1;
      pivot--;
      continue;
    }

    // If we have successfully incremented, check to be
    // sure that there are no duplicate atoms.         
    for (i = 0; i < pivot; i++) {
      if (currvec[i] == currvec[pivot]) {
        inci = 1;
      }
    }
  }

  // The pivot may have shifted backwards
  *nassigned = pivot+1;

  // Incrementing was successful and the new   
  // state is a candidate for energy evaluation      
  if (inci == 0) {
    return 0;
  }

  // The pivot had to be pushed all the way back to 
  // -1, indicating that all states have been tested
  else {
    return 1;
  }
}

//-----------------------------------------------------------------------------
// WriteInpcrd: write an inpcrd file in AMBER format.  This routine is here
//              because it is not currently needed anywhere else in the
//              program.
//
// Arguments:                                                           
//   crd:      the coordinates                                          
//   bdim:     the box dimensions                                       
//   finame:   the name of the inpcrd file to write                     
//   tj:       the trajectory control data (for overwrite flag)         
//-----------------------------------------------------------------------------
static void WriteInpcrd(coord *crd, double* bdim, char* finame, trajcon *tj)
{
  int h, i;
  FILE *outp;

  outp = FOpenSafe(finame, tj->OverwriteOutput);
  fprintf(outp, "Coordinates generated by mdgx\n%6d\n", crd->natom);
  h = 0;
  for (i = 0; i < crd->natom; i++) {
    fprintf(outp, "%12.7lf%12.7lf%12.7lf", crd->loc[3*i], crd->loc[3*i+1],
            crd->loc[3*i+2]);
    h++;
    if (h == 2) {
      fprintf(outp, "\n");
      h = 0;
    }
  }
  if (h != 0) {
    fprintf(outp, "\n");
  }
  for (i = 0; i < 6; i++) {
    fprintf(outp, "%12.7lf", bdim[i]);
  }
  fprintf(outp, "\n");
  fclose(outp);
}

//-----------------------------------------------------------------------------
// RearrangeMol: rearrange the atoms of a molecule to minimize molecular
//               mechanics bonded term contributions.                   
//                                                                      
// Arguments:                                                           
//   conf:   the molecular mechanics system, complete with topology,    
//           coordinates, and energy tables                             
//   mp:     the parameter set, containing master tables of bonded terms
//   tj:     trajectory control data (for overwrite flag)               
//-----------------------------------------------------------------------------
static void RearrangeMol(mmsys *conf, prmset *mp, trajcon *tj)
{
  int i, j, atma, atmb, natom, nfound, nadd, nassigned, finished, strain;
  int iter, relaxed;
  int* atmorder;
  int* atmfound;
  int* currvec;
  int* bestvec;
  int* assignments;
  double dw, dx, dy, dz, currnrg;
  double bdim[6];
  double* bestnrg;
  dmat bndmat, dstmat, stiffmat;
  coord locbuff;
  coord *crd;
  xbonddef *bndtmp;

  // Make a matrix of bonded terms
  natom = conf->tp->natom;
  bndmat = CreateDmat(natom, natom, 0);
  stiffmat = CreateDmat(natom, natom, 0);
  for (i = 0; i < conf->bmap.nbond; i++) {
    atma = conf->bmap.id[i].a;
    atmb = conf->bmap.id[i].b;
    bndtmp = &mp->bonds[conf->bmap.id[i].key];
    bndmat.map[atma][atmb] = bndtmp->l0;
    bndmat.map[atmb][atma] = bndtmp->l0;
    stiffmat.map[atma][atmb] = bndtmp->K;
    stiffmat.map[atmb][atma] = bndtmp->K;
  }

  // The best energy after adding each atom
  bestnrg = (double*)malloc(natom*sizeof(double));
  SetDVec(bestnrg, natom, 1.0e30);

  // Compute the inter-atomic distances in the coordinate set
  crd = &conf->crd;
  dstmat = CreateDmat(natom, natom, 0);
  for (i = 0; i < natom; i++) {
    for (j = 0; j < natom; j++) {
      dx = crd->loc[3*j] - crd->loc[3*i];
      dy = crd->loc[3*j+1] - crd->loc[3*i+1];
      dz = crd->loc[3*j+2] - crd->loc[3*i+2];
      dstmat.map[i][j] = sqrt(dx*dx + dy*dy + dz*dz);
    }
  }

  // Plot a course through the molecule: in      
  // what order shall we add atoms to the system?
  atmorder = (int*)malloc(natom*sizeof(int));
  atmfound = (int*)calloc(natom, sizeof(int));
  atmorder[0] = 0;
  atmfound[0] = 1;
  nfound = 1;
  while (nfound < natom) {
    nadd = 0;
    for (i = 0; i < nfound; i++) {
      for (j = 0; j < conf->bmap.nbond; j++) {
        if (conf->bmap.id[j].a == atmorder[i] &&
          atmfound[conf->bmap.id[j].b] == 0) {
          atmfound[conf->bmap.id[j].b] = 1;
          atmorder[nfound+nadd] = conf->bmap.id[j].b;
          nadd++;
        }
        if (conf->bmap.id[j].b == atmorder[i] &&
          atmfound[conf->bmap.id[j].a] == 0) {
          atmfound[conf->bmap.id[j].a] = 1;
          atmorder[nfound+nadd] = conf->bmap.id[j].a;
          nadd++;
        }
      }
    }
    nfound += nadd;
  }

  // Copy the coordinates for safe keeping
  locbuff = CopyCoord(&conf->crd);

  // Find putative identities for each atom based on bonding.
  // The vectors currvec and bestvec record the current and  
  // best configurations found thus far.  These vectors do   
  // NOT hold the identity of kth atom in position k, rather 
  // they hold the identity of the atmorder[k]th atom at     
  // position k.  The molecule is filled out according to    
  // the array atmorder.                                     
  currvec = (int*)malloc(natom*sizeof(int));
  bestvec = (int*)malloc(natom*sizeof(int));
  assignments = (int*)calloc(natom, sizeof(int));
  SetIVec(assignments, natom, -1);
  SetIVec(currvec, natom, -1);
  nassigned = 2;
  currvec[0] = 0;
  currvec[1] = 1;
  finished = 0;
  iter = 0;
  relaxed = 0;
  while (finished == 0) {

    // Test this configuration to the extent 
    // that atom positions have been assigned
    PlaceAtoms(currvec, atmorder, assignments, nassigned, conf, &locbuff);
    currnrg = 0.0;
    strain = 0;
    for (i = 0; i < conf->bmap.nbond; i++) {
      if (assignments[conf->bmap.id[i].a] >= 0 &&
          assignments[conf->bmap.id[i].b] >= 0) {
        strain += BondTermParse(conf, i, 2, mp, &dz, &dx, &dw, &dw, &dw, &dy);
        currnrg += dy;
      }
    }
    for (i = 0; i < conf->amap.nangl; i++) {
      if (assignments[conf->amap.id[i].a] >= 0 &&
          assignments[conf->amap.id[i].b] >= 0 &&
          assignments[conf->amap.id[i].c] >= 0) {
        strain += BondTermParse(conf, i, 3, mp, &dz, &dx, &dw, &dw, &dw, &dy);
        currnrg += dy;
      }
    }
    for (i = 0; i < conf->hmap.ntterm; i++) {
      if (assignments[conf->hmap.id[i].a] >= 0 &&
          assignments[conf->hmap.id[i].b] >= 0 &&
          assignments[conf->hmap.id[i].c] >= 0 &&
          assignments[conf->hmap.id[i].d] >= 0) {
        strain += BondTermParse(conf, i, 4, mp, &dz, &dx, &dw, &dw, &dw, &dy);
        currnrg += dy;
      }
    }

    // Now, the energy has been summed; is it
    // acceptable, and if so is it the best?
    if (strain == 0) {
      if (currnrg < bestnrg[nassigned-1]) {
        bestnrg[nassigned-1] = currnrg;
        if (nassigned == natom) {
          relaxed = 1;
          ReflectIVec(bestvec, currvec, natom);
        }
      }
      if (nassigned < natom) {
        currvec[nassigned] = -1;
        nassigned++;
      }
    }
    finished = IncrementState(currvec, &nassigned, natom);
    iter++;
  }
  if (relaxed == 0) {
    printf("mdgx >>   Optimization unsuccessful after %d iterations.\n", iter);
    conf->PassedEtol = 0;
  }
  else {
    printf("mdgx >>   Optimization complete after %d iterations.\n", iter);
    printf("mdgx >>   Best energy %12.6lf\n", bestnrg[natom-1]);
    conf->PassedEtol = 2;
    for (i = 0; i < natom; i++) {
      for (j = 0; j < 3; j++) {
        conf->crd.loc[3*atmorder[i]+j] = locbuff.loc[3*bestvec[i]+j];
      }
    }
    AllBondedTerms(conf, mp, &strain);
    for (i = 0; i < 3; i++) {
      bdim[i] = 100.0;
      bdim[i+3] = 90.0;
    }
    if (tj->OverwriteOutput == 1) {
      WriteInpcrd(&conf->crd, bdim, conf->crdsrc, tj);
    }
  }

  // Free allocated memory
  free(atmfound);
  free(bestnrg);
  free(currvec);
  free(bestvec);
  DestroyCoord(&locbuff);
}

//-----------------------------------------------------------------------------
// CompBondedMME: compute molecular mechanics energy due to bonded
//                interactions.                                         
//                                                                      
// Arguments:                                                           
//   conf:   the molecular mechanics system, complete with topology,    
//           coordinates, and energy tables                             
//   mp:     the parameter set, containing master tables of bonded terms
//   tj:     trajectory control data (for passing overwrite flag to     
//           RearrangeMol if needed)                                    
//-----------------------------------------------------------------------------
static void CompBondedMME(mmsys *conf, prmset *mp, trajcon *tj)
{
  int strain;
  double mmnrg;

  // Check for high molecular mechanics energy
  mmnrg = AllBondedTerms(conf, mp, &strain);
  if (strain > 0) {
    printf("\nmdgx >> Rearrangement required for conformation:\nmdgx >> "
           "  Topology %s / Coordinates %s.\n"
           "mdgx >>   %d strained interactions, total bonded energy "
           "%12.6lf\n", conf->tp->source, conf->crdsrc, strain, mmnrg);
    RearrangeMol(conf, mp, tj);
  }
  else {
    conf->PassedEtol = 1;
  }
}

//-----------------------------------------------------------------------------
// CompExclMatrix: compute the exclusion and scaling matrix for a small system,
//                 with electrostatic attenuations in the upper triangle (row
//                 index i is greater than column index j) and Lennard-Jones
//                 attenuations in the lower triangle.
//
// Arguments:
//   tp:      the system topology
//   partnb:  the (pre-allocated) non-bonded partial interaction matrix, which
//            reads 1 if there is a 1:4 interaction between two atoms and zero
//            otherwise.  If the pointer is fed in as NULL this will not be
//            computed.
//-----------------------------------------------------------------------------
dmat CompExclMatrix(prmtop *tp, imat *partnb)
{
  int i, j, k, m, p, ni, nj, parni, parnj, rootAtom;
  double fscee, fscnb;
  dmat excl;

  // Assume all interactions are complete to start with
  excl = CreateDmat(tp->natom, tp->natom, 0);
  SetDVec(excl.data, tp->natom*tp->natom, 1.0);

  // Completely eliminate 1:2 and 1:3 interactions,
  // and zero out 1:4 interactions as well.
  for (i = 0; i < tp->natom-1; i++) {
    for (j = i+1; j < tp->natom; j++) {
      if (TestBondedExclusion(i, j, tp) == 1) {
        excl.map[i][j] = 0.0;
        excl.map[j][i] = 0.0;
      }
    }
  }

  // Mark 1:4 attenuations
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->HLC[i].ndihe; j++) {
      if (tp->HLC[i].HC[j].eval14 == 1) {
        ni = tp->HLC[i].HC[j].a;
        nj = tp->HLC[i].HC[j].d;
        if (ni < nj) {
          excl.map[ni][nj] = 1.0 - tp->HLC[i].HC[j].scee;
          excl.map[nj][ni] = 1.0 - tp->HLC[i].HC[j].scnb;
        }
        else {
          excl.map[nj][ni] = 1.0 - tp->HLC[i].HC[j].scee;
          excl.map[ni][nj] = 1.0 - tp->HLC[i].HC[j].scnb;
        }
        if (partnb != NULL) {
          partnb->map[ni][nj] = 1;
          partnb->map[nj][ni] = 1;
        }
      }
    }
  }

  // Mark any other 1:1, 1:2, 1:3, and 1:4 eliminations
  // (due to virtual sites, mainly)
  if (tp->EPInserted == 1 || tp->ExclMarked == 1) {
    for (i = 0; i < tp->natom; i++) {
      for (j = 0; j < tp->ElimPair[i].n11; j++) {
        ni = tp->ElimPair[i].list11[j].atmX;
        nj = tp->ElimPair[i].list11[j].atmY;
        excl.map[ni][nj] = 0.0;
        excl.map[nj][ni] = 0.0;
      }
      for (j = 0; j < tp->ElimPair[i].n12; j++) {
        ni = tp->ElimPair[i].list12[j].atmX;
        nj = tp->ElimPair[i].list12[j].atmY;
        excl.map[ni][nj] = 0.0;
        excl.map[nj][ni] = 0.0;
      }
      for (j = 0; j < tp->ElimPair[i].n13; j++) {
        ni = tp->ElimPair[i].list13[j].atmX;
        nj = tp->ElimPair[i].list13[j].atmY;
        excl.map[ni][nj] = 0.0;
        excl.map[nj][ni] = 0.0;
      }
      for (j = 0; j < tp->ElimPair[i].n14; j++) {
        ni = tp->ElimPair[i].list14[j].atmX;
        nj = tp->ElimPair[i].list14[j].atmY;

        // Loop over atoms ni and nj, determine whether they are
        // massive atoms themselves or which massive atoms are
        // their parents, then determine whether there is a
        // particular 1:4 scaling factor at work between those
        // two massive atoms and apply that scaling factor to
        // the exclusions matrix.
        parni = (tp->Masses[ni] > 1.0e-8) ? ni : tp->nb1234[ni].L11[0];
        parnj = (tp->Masses[nj] > 1.0e-8) ? nj : tp->nb1234[nj].L11[0];
        fscee = tp->elec14fac;
        fscnb = tp->lj14fac;
        for (k = 0; k < tp->nb1234[parni].n14; k++) {
          if (tp->nb1234[parni].L14[k] == parnj) {

            // This is a 1:4 interaction.  Assign it to the standard
            // 1:4 scaling factors, then look for dihedrals that could
            // indicate specific factors.
            for (m = 0; m < tp->nb1234[parni].n12; m++) {
              rootAtom = tp->nb1234[parni].L12[m];
              for (p = 0; p < tp->HLC[rootAtom].ndihe; p++) {
                if ((tp->HLC[rootAtom].HC[p].a == parni &&
                     tp->HLC[rootAtom].HC[p].d == parnj) ||
                    (tp->HLC[rootAtom].HC[p].a == parnj &&
                     tp->HLC[rootAtom].HC[p].d == parni)) {
                  fscee = tp->HLC[rootAtom].HC[p].scee;
                  fscnb = tp->HLC[rootAtom].HC[p].scnb;
                }
              }
            }
            for (m = 0; m < tp->nb1234[parnj].n12; m++) {
              rootAtom = tp->nb1234[parnj].L12[m];
              for (p = 0; p < tp->HLC[rootAtom].ndihe; p++) {
                if ((tp->HLC[rootAtom].HC[p].a == parni &&
                     tp->HLC[rootAtom].HC[p].d == parnj) ||
                    (tp->HLC[rootAtom].HC[p].a == parnj &&
                     tp->HLC[rootAtom].HC[p].d == parni)) {
                  fscee = tp->HLC[rootAtom].HC[p].scee;
                  fscnb = tp->HLC[rootAtom].HC[p].scnb;
                }
              }
            }
          }
        }
        if (ni < nj) {
          excl.map[ni][nj] = 1.0 - fscee;
          excl.map[nj][ni] = 1.0 - fscnb;
        }
        else {
          excl.map[nj][ni] = 1.0 - fscee;
          excl.map[ni][nj] = 1.0 - fscnb;
        }
      }
    }
  }

  return excl;
}

//-----------------------------------------------------------------------------
// CompNBMatrix: compute the matrix of non-bonded interactions for this 
//               molecule.  After this computation, the matrix will be parsed
//               for exclusions and 1:4 interactions.            
//                                                                      
// Arguments:                                                           
//   mp:     the parameter set, containing fitting directives           
//   conf:   the molecular mechanics system, complete with topology,
//           coordinates, and energy tables                             
//-----------------------------------------------------------------------------
static void CompNBMatrix(prmset *mp, mmsys *conf)
{
  int i, j, ni, nj, natom;
  double atmx, atmy, atmz, atmq, dx, dy, dz, invr, invr2, invr6, invr12;
  double Aval, Bval;
  double* ljAtmp;
  imat partnb;
  dmat excl, nrg;
  prmtop *tp;
  coord *crd;

  // Unpack the system
  tp = conf->tp;
  crd = &conf->crd;
  natom = tp->natom;

  // Allocate tables
  nrg = CreateDmat(natom, natom, 0);
  partnb = CreateImat(natom, natom);

  // Allocate the matrix.  Electrostatic interactions
  // go above the diagonal, van-der Waals below it.  
  for (i = 0; i < tp->natom-1; i++) {
    atmx = crd->loc[3*i];
    atmy = crd->loc[3*i+1];
    atmz = crd->loc[3*i+2];
    atmq = BIOQ*tp->Charges[i];
    if (tp->LJIdx[i] >= 0) {
      ljAtmp = tp->LJutab.map[tp->LJIdx[i]];
    }
    for (j = i+1; j < tp->natom; j++) {
      dx = crd->loc[3*j] - atmx;
      dy = crd->loc[3*j+1] - atmy;
      dz = crd->loc[3*j+2] - atmz;
      invr2 = 1.0/(dx*dx + dy*dy + dz*dz);
      invr = sqrt(invr2);
      invr6 = invr2*invr2*invr2;
      invr12 = invr6*invr6;

      // Compute electrostatic interaction
      nrg.map[i][j] = atmq*tp->Charges[j]*invr;

      // Compute Lennard-Jones interaction
      nrg.map[j][i] = 0.0;
      if (tp->LJIdx[i] >= 0 && tp->LJIdx[j] >= 0) {
        Aval = ljAtmp[2*tp->LJIdx[j]];
        Bval = ljAtmp[2*tp->LJIdx[j]+1];
        nrg.map[j][i] = Aval*invr12 + Bval*invr6;
      }
    }
  }

  // Compute the exclusion matrix for this system,
  // and fill the partial 1:4 interactions table.
  excl = CompExclMatrix(tp, &partnb);

  // Compute the non-bonded electrostatic and Lennard-Jones kernels
  conf->EEkernel = 0.0;
  conf->LJkernel = 0.0;
  conf->EEnonfit = 0.0;
  conf->LJnonfit = 0.0;
  for (i = 0; i < tp->natom-1; i++) {
    for (j = i+1; j < tp->natom; j++) {
      if (partnb.map[i][j] == 1) {
        if (mp->fitscee == 1) {
          conf->EEkernel += nrg.map[i][j];
        }
        else {
          conf->EEnonfit += nrg.map[i][j] * excl.map[i][j];
        }
        if (mp->fitscnb == 1) {
          conf->LJkernel += nrg.map[j][i];
        }
        else {
          conf->LJnonfit += nrg.map[j][i] * excl.map[j][i];
        }
      }
      else {
        conf->EEnonfit += nrg.map[i][j] * excl.map[i][j];
        conf->LJnonfit += nrg.map[j][i] * excl.map[j][i];
      }
    }
  }

  // Commit tables to the system
  conf->nbnrg = nrg;
  conf->excl = excl;

  // Free allocated memory
  DestroyImat(&partnb);
}

//-----------------------------------------------------------------------------
// SumEnergyMM: sum the energy according to a molecular mechanics model.
//                                                                      
// Arguments:                                                           
//   conf:     the system of interest                                   
//-----------------------------------------------------------------------------
static double SumEnergyMM(mmsys *conf)
{
  int i, j, natom;
  double eQQ, eLJ, etot;

  eQQ = 0.0;
  eLJ = 0.0;
  natom = conf->tp->natom;
  for (i = 0; i < natom-1; i++) {
    for (j = i+1; j < natom; j++) {
      eQQ += conf->nbnrg.map[i][j] * conf->excl.map[i][j];
      eLJ += conf->nbnrg.map[j][i] * conf->excl.map[j][i];
    }
  }
  etot = DSum(conf->bmap.Ucontrib, conf->bmap.nbond);
  etot += DSum(conf->amap.Ucontrib, conf->amap.nangl);
  etot += DSum(conf->hmap.Ucontrib, conf->hmap.ntterm);
  etot += eQQ + eLJ;

  return etot;
}

//-----------------------------------------------------------------------------
// GetInstanceCounts: function to produce a vector of integers recording the
//                    number of times that a variable for any column in the
//                    fitting matrix is expected to have a non-zero entry in
//                    that matrix.
//
// Arguments:
//   mp:       fitting control data (contains atom types and examples of each
//             atom in various topologies)
//-----------------------------------------------------------------------------
static int* GetInstanceCounts(prmset *mp)
{
  int i, j, k, fcon, colidx, colfin;
  int* result;

  result = (int*)calloc(mp->nparm + mp->nunisys, sizeof(int));
  for (i = 0; i < mp->nunisys; i++) {
    fcon = mp->FirstConf[i];
    for (j = 0; j < mp->conf[fcon].bmap.nbond; j++) {
      colidx = mp->bonds[mp->conf[fcon].bmap.id[j].key].fitcolX;
      if (colidx >= 0) {
	result[colidx] += mp->GroupCount[i];
      }
      colidx = mp->bonds[mp->conf[fcon].bmap.id[j].key].fitcolY;
      if (colidx >= 0) {
	result[colidx] += mp->GroupCount[i];
      }
    }
    for (j = 0; j < mp->conf[fcon].amap.nangl; j++) {
      colidx = mp->angls[mp->conf[fcon].amap.id[j].key].fitcolX;
      if (colidx >= 0) {
	result[colidx] += mp->GroupCount[i];
      }
      colidx = mp->angls[mp->conf[fcon].amap.id[j].key].fitcolY;
      if (colidx >= 0) {
	result[colidx] += mp->GroupCount[i];
      }      
    }
    for (j = 0; j < mp->conf[fcon].hmap.ntterm; j++) {
      colidx = mp->torsions[mp->conf[fcon].hmap.id[j].key].fitcol;
      if (colidx >= 0) {
	result[colidx] += mp->GroupCount[i];
      }
    }
    for (j = 0; j < mp->conf[fcon].rmap.nops; j++) {
      colidx = mp->nmrops[mp->conf[fcon].rmap.id[j].key].fitcol;
      if (colidx >= 0) {
	result[colidx] += mp->GroupCount[i];
      }
    }
    for (j = 0; j < mp->conf[fcon].mmap.ncmap; j++) {
      colidx = mp->cmaps[mp->conf[fcon].mmap.id[j].key].fitcoli;
      if (colidx >= 0) {
	colfin = mp->cmaps[mp->conf[fcon].mmap.id[j].key].fitcolf;
	for (k = colidx; k < colfin; k++) {
          result[k] += mp->GroupCount[i];
	}
      }
    }
  }

  return result;
}

//-----------------------------------------------------------------------------
// CountInstances: count the instances of a variable occurring in the fitting
//                 matrix.                                      
//                                                                      
// Arguments:                                                           
//   A:        the fitting matrix, in its original form before the QR   
//             decomposition                                            
//   rstart:   the starting row of the fitting matrix, when we switch from data
//             to restraints                                  
//   ridx:     the (column) index of the variable of interest           
//-----------------------------------------------------------------------------
static int CountInstances(dmat *A, int rstart, int ridx)
{
  int i, ninst;

  ninst = 0;
  for (i = 0; i < rstart; i++) {
    if (fabs(A->map[i][ridx]) > 1.0e-12) {
      ninst += 1;
    }
  }
  if (ninst == 0) {
    ninst = 1;
  }

  return ninst;
}

//-----------------------------------------------------------------------------
// FitMatrixRestraints: add constraints to the fitting matrix of bonds, angles,
//                      and dihedral potentials.                
//                                                                      
// Arguments:
//   mp:         fitting control data (contains atom types and examples of each
//               atom in various topologies)
//   A:          the fitting matrix, in its original form before the QR
//               decomposition
//   b:          the target vector for the equation Ax = b
//   rstart:     the starting row of the fitting matrix, when we switch from
//               data to restraints
//   InstTally:  the number of times each variable is expected to appear in the
//               fitting matrix
//-----------------------------------------------------------------------------
static void FitMatrixRestraints(prmset *mp, dmat *A, double* b, int rstart,
				int* InstTally)
{
  int h, i, j, ninst, vcon, ngrd, colA, colB;
  double Acnfac, xdiff, ydiff;

  // The matrix row at which to start adding restraints
  vcon = rstart;

  // Harmonic restraints to bond stiffnesses and equilibria
  for (h = 0; h < mp->nbond; h++) {
    if (mp->bonds[h].fitcolX < 0 || mp->bonds[h].rstwK < 1.0e-12) {
      continue;
    }
    ninst = (mp->RestScaleMeth == 0) ? InstTally[mp->bonds[h].fitcolX] :
            CountInstances(A, rstart, mp->bonds[h].fitcolX);
    Acnfac = ninst*mp->bonds[h].rstwK;
    A->map[vcon][mp->bonds[h].fitcolX] = Acnfac;
    b[vcon] = Acnfac*mp->bonds[h].targK;
    mp->bonds[h].restKeq = vcon;
    if (mp->FitBondEq == 1) {

      // The sum of stiffnesses from each basis function 
      // produces the final stiffness.  Restrain the sum.
      A->map[vcon][mp->bonds[h].fitcolY] = Acnfac;

      // Second constraint on the ratio of the two stiffness
      // constants to move the solution towards the desired 
      // equilibrium value.                                 
      vcon++;
      xdiff = mp->bonds[h].targl0 - mp->bonds[h].lbasisX;
      ydiff = mp->bonds[h].targl0 - mp->bonds[h].lbasisY;
      Acnfac = ninst*mp->bonds[h].rstwl0;
      A->map[vcon][mp->bonds[h].fitcolX] = Acnfac*xdiff;
      A->map[vcon][mp->bonds[h].fitcolY] = Acnfac*ydiff;
      b[vcon] = 0.0;
      mp->bonds[h].restLeq = vcon;
    }
    vcon++;
  }

  // Restraints to angle stiffnesses and equilibria
  for (h = 0; h < mp->nangl; h++) {
    if (mp->angls[h].fitcolX < 0 || mp->angls[h].rstwK < 1.0e-12) {
      continue;
    }
    ninst = (mp->RestScaleMeth == 0) ? InstTally[mp->angls[h].fitcolX] :
            CountInstances(A, rstart, mp->angls[h].fitcolX);
    Acnfac = ninst*mp->angls[h].rstwK;
    A->map[vcon][mp->angls[h].fitcolX] = Acnfac;
    b[vcon] = Acnfac*mp->angls[h].targK;
    mp->angls[h].restKeq = vcon;
    if (mp->FitAnglEq == 1) {

      // The sum of stiffnesses from each basis function 
      // produces the final stiffness.  Restrain the sum.
      A->map[vcon][mp->angls[h].fitcolY] = Acnfac;

      // Second constraint on the ratio of the two stiffness
      // constants to move the solution towards the desired 
      // equilibrium value.                                 
      vcon++;
      xdiff = mp->angls[h].targTh0 - mp->angls[h].tbasisX;
      ydiff = mp->angls[h].targTh0 - mp->angls[h].tbasisY;
      Acnfac = ninst*mp->angls[h].rstwTh0;
      A->map[vcon][mp->angls[h].fitcolX] = Acnfac*xdiff;
      A->map[vcon][mp->angls[h].fitcolY] = Acnfac*ydiff;
      b[vcon] = 0.0;
      mp->angls[h].restTh0 = vcon;
    }
    vcon++;
  }

  // Restraints to torsion amplitudes
  for (h = 0; h < mp->ntor; h++) {
    if (mp->torsions[h].fitcol < 0 || mp->torsions[h].rstw < 1.0e-12) {
      continue;
    }
    ninst = (mp->RestScaleMeth == 0) ? InstTally[mp->torsions[h].fitcol] :
            CountInstances(A, rstart, mp->torsions[h].fitcol);
    Acnfac = ninst*mp->torsions[h].rstw;
    if (mp->torsions[h].impr == 0) {
      A->map[vcon][mp->torsions[h].fitcol] = Acnfac;
      b[vcon] = Acnfac * mp->torsions[h].target;
    }
    else {
      A->map[vcon][mp->torsions[h].fitcol] = Acnfac;
      b[vcon] = Acnfac * mp->torsions[h].K;
    }
    mp->torsions[h].restKval = vcon;
    vcon++;
  }

  // Restraints to 1-4 scaling factors
  if (mp->fitscee == 1 && mp->grst14 > 1.0e-12) {
    h = mp->nbvar + mp->navar + mp->nhvar;
    ninst = CountInstances(A, rstart, h);
    A->map[vcon][h] = ninst*mp->grst14;
    b[vcon] = ninst*mp->grst14*mp->elec14fac;
    vcon++;
  }
  if (mp->fitscnb == 1 && mp->grst14 > 1.0e-12) {
    h = mp->nbvar + mp->navar + mp->nhvar + mp->fitscee;
    ninst = CountInstances(A, rstart, h);
    A->map[vcon][h] = ninst*mp->grst14;
    b[vcon] = ninst*mp->grst14*mp->lj14fac;
    vcon++;
  }

  // Restraints to CMAP surfaces
  for (h = 0; h < mp->ncmaps; h++) {
    ngrd = mp->cmaps[h].esrfFit.row;
    if (mp->cmaps[h].fitcoli < 0 || mp->cmaps[h].rstw < 1.0e-12) {
      continue;
    }
    colA = mp->cmaps[h].fitcoli;
    ninst = (mp->RestScaleMeth == 0) ? InstTally[colA] :
            CountInstances(A, rstart, colA);
    A->map[vcon][colA] = ninst * mp->cmaps[h].rstw;
    b[vcon] = 0.0;
    vcon++;
    for (i = 0; i < ngrd; i++) {
      for (j = 0; j < ngrd; j++) {
        ninst = (mp->RestScaleMeth == 0) ? InstTally[colA] :
                CountInstances(A, rstart, colA);
	colB = (j < ngrd-1) ? colA + 1 : colA - ngrd + 1;
	A->map[vcon][colA] =  ninst * mp->cmaps[h].rstw;
	A->map[vcon][colB] = -ninst * mp->cmaps[h].rstw;
	b[vcon] = 0.0;
	vcon++;
	colB = (i < ngrd-1) ? colA + ngrd : colA - ngrd*(ngrd - 1);
	A->map[vcon][colA] =  ninst * mp->cmaps[h].rstw;
	A->map[vcon][colB] = -ninst * mp->cmaps[h].rstw;
	b[vcon] = 0.0;
	vcon++;
	colA++;
      }
    }
  }
  
  // Restraints to NMR optimization potentials
  for (h = 0; h < mp->nops; h++) {
    if (mp->nmrops[h].fitcol < 0 || mp->nmrops[h].rstw < 1.0e-12) {
      continue;
    }
    ninst = InstTally[mp->nmrops[h].fitcol];
    Acnfac = ninst*mp->nmrops[h].rstw;
    A->map[vcon][mp->nmrops[h].fitcol] = Acnfac;
    b[vcon] = Acnfac * mp->nmrops[h].target;
    mp->nmrops[h].restKval = vcon;
    vcon++;
  }
}

//-----------------------------------------------------------------------------
// UnrelatedDihedrals: given two columns of the fitting matrix, return 1 if
//                     they describe the actions of two unrelated dihedrals. 
//                     Return 0 otherwise.                  
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains history and number of fitting
//             points)                                          
//   n[i,j]:   the numbers of the fitting matrix columns                
//-----------------------------------------------------------------------------
static int UnrelatedDihedrals(prmset *mp, int ni, int nj)
{
  int i, idihe, jdihe;

  if (ni == nj) {
    return 0;
  }
  idihe = -1;
  jdihe = -1;
  for (i = 0; i < mp->nhvar; i++) {
    if (mp->torsions[i].fitcol == ni) {
      idihe = i;
    }
    if (mp->torsions[i].fitcol == nj) {
      jdihe = i;
    }
  }
  if (idihe < 0 || jdihe < 0) {
    return 0;
  }
  if (TypeCompare(mp->torsions[idihe].btype, mp->torsions[idihe].ctype, "    ",
                  "    ", mp->torsions[jdihe].btype, mp->torsions[jdihe].ctype,
                  "    ", "    ", 2, 1) == 2) {
    return 0;
  }

  return 1;
}

//-----------------------------------------------------------------------------
// FillFittingMatrix: fill up the fitting matrix given the set of adjustable
//                    parameters, known evergy targets, and molecular mechanics
//                    properties of each system.  Also loads energy targets
//                    into a pre-allocated array.
//
// Arguments:
//   mp:         fitting control data (contains history and number of fitting
//               points)
//   A:          the pre-allocated fitting matrix
//   b:          the pre-allocated array for computing Ax = B
//   InstTally:  the number of times each variable is expected to appear in the
//               fitting matrix
//-----------------------------------------------------------------------------
static void FillFittingMatrix(prmset *mp, dmat *A, double* b, int* InstTally)
{
  int h, i, j, k, cidx;
  double colsum, ctest, sq, dtj;
  double *dtmp, *dtm2p;
  double* vec_a;
  double* vec_a2;
  xbonddef *bondkey;
  xangldef *anglkey;
  torterm *torkey;
  xcmapdef *cmapkey;
  nmroper *NMRkey;
  dmat vec_ab, cmfprnt, cmtprof;

  // If there are CMAP terms, allocate for the biggest of them
  if (mp->ncmaps > 0) {
    h = 0;
    for (i = 0; i < mp->ncmaps; i++) {
      if (mp->cmaps[i].esrfFit.row > h) {
        h = mp->cmaps[i].esrfFit.row;
      }
    }
    cmtprof = CreateDmat(h, h, 0);
    cmfprnt = CreateDmat(h, h, 0);
  }
  
  // Loop over all systems
  for (h = 0; h < mp->nconf; h++) {

    // Start accumulating the MM energy components
    // which are not subject to fitting           
    mp->conf[h].nonfitmm = 0.0;

    // Bonds
    for (i = 0; i < mp->conf[h].bmap.nbond; i++) {
      bondkey = &mp->bonds[mp->conf[h].bmap.id[i].key];
      cidx = bondkey->fitcolX;
      if (cidx < 0) {
        mp->conf[h].nonfitmm += mp->conf[h].bmap.Ucontrib[i];
      }
      else {
        A->map[h][cidx] += mp->conf[h].bmap.UkernelX[i];
        if (mp->FitBondEq == 1) {
          cidx = bondkey->fitcolY;
          A->map[h][cidx] += mp->conf[h].bmap.UkernelY[i];
        }
        cidx = bondkey->fitcolZ;
        if (cidx >= 0) {
          A->map[h][cidx] += mp->conf[h].bmap.UkernelZ[i];
        }
        cidx = bondkey->fitcolW;
        if (cidx >= 0) {
          A->map[h][cidx] += mp->conf[h].bmap.UkernelW[i];
        }
      }
    }

    // Angles
    for (i = 0; i < mp->conf[h].amap.nangl; i++) {
      anglkey = &mp->angls[mp->conf[h].amap.id[i].key];
      cidx = anglkey->fitcolX;
      if (cidx < 0) {
        mp->conf[h].nonfitmm += mp->conf[h].amap.Ucontrib[i];
      }
      else {
        A->map[h][cidx] += mp->conf[h].amap.UkernelX[i];
        if (mp->FitAnglEq == 1) {
          cidx = anglkey->fitcolY;
          A->map[h][cidx] += mp->conf[h].amap.UkernelY[i];
        }
      }
    }

    // Torsions
    for (i = 0; i < mp->conf[h].hmap.ntterm; i++) {
      torkey = &mp->torsions[mp->conf[h].hmap.id[i].key];
      cidx = torkey->fitcol;
      if (cidx < 0) {
        mp->conf[h].nonfitmm += mp->conf[h].hmap.Ucontrib[i];
      }
      else {
        A->map[h][cidx] += mp->conf[h].hmap.Ukernel[i];
      }
    }

    // CMAP terms
    for (i = 0; i < mp->conf[h].mmap.ncmap; i++) {
      cmapkey = &mp->cmaps[mp->conf[h].mmap.id[i].key];
      cidx = cmapkey->fitcoli;
      if (cidx < 0) {
	mp->conf[h].nonfitmm += mp->conf[h].mmap.Ucontrib[i];
      }
      else {
        GetCMAPFootPrint(mp->conf[h].mmap.val[2*i],
                         mp->conf[h].mmap.val[2*i + 1], &cmapkey->stencil,
                         &cmtprof, &cmfprnt);
        for (j = 0; j < cmapkey->stencil.col; j++) {
          for (k = 0; k < cmapkey->stencil.col; k++) {
	    A->map[h][cidx] += cmfprnt.map[j][k];
            cidx++;
	  }
        }
      }
    }
    
    // NMR restraints
    for (i = 0; i < mp->conf[h].rmap.nops; i++) {
      NMRkey = &mp->nmrops[mp->conf[h].rmap.id[i].key];
      cidx = NMRkey->fitcol;
      if (cidx < 0) {
        mp->conf[h].nonfitmm += mp->conf[h].rmap.Ucontrib[i];
      }
      else {
        A->map[h][cidx] += mp->conf[h].rmap.Ukernel[i];
      }
    }

    // Non-bonded energy
    mp->conf[h].nonfitmm += mp->conf[h].EEnonfit;
    mp->conf[h].nonfitmm += mp->conf[h].LJnonfit;
    if (mp->fitscee == 1) {
      cidx = mp->nbvar + mp->navar + mp->nhvar;
      A->map[h][cidx] = mp->conf[h].EEkernel;
    }
    if (mp->fitscnb == 1) {
      cidx = mp->nbvar + mp->navar + mp->nhvar + 1;
      A->map[h][cidx] = mp->conf[h].LJkernel;
    }

    // Constant for energy offset between MM and QM models
    A->map[h][mp->nparm + mp->conf[h].GroupNum] = 1.0;

    // Energy target
    if (mp->zeroNonfit == 1) {
      mp->conf[h].nonfitmm = 0.0;
    }
    b[h] = mp->conf[h].enorm - mp->conf[h].nonfitmm;

    // Report to the user
    if (mp->verbose == 1 && (h % 1000) == 0) {
      fprintf(stderr, "\rmdgx >> Evaluating conformation %7d of %7d", h,
	      mp->nconf);
      fflush(stderr);
    }
  }

  // Add restraint conditions
  if (mp->verbose == 1) {
    printf("\nmdgx >> Assigning restraints.\n");
  }
  FitMatrixRestraints(mp, A, b, mp->nconf, InstTally);

  // Check the matrix A for zero and highly correlated columns
  mp->nzerocol = 0;
  mp->zerocol = (int*)malloc(A->col*sizeof(int));
  mp->ncorrcol = 0;
  mp->corrcol = (int*)malloc((A->col+1)*A->col*sizeof(int));
  mp->corrval = (double*)malloc(((A->col+1)*A->col/2)*sizeof(double));

  // Pearson correlations computed en masse to save time
  vec_a = (double*)calloc(mp->nparm, sizeof(double));
  vec_a2 = (double*)calloc(mp->nparm, sizeof(double));
  vec_ab = CreateDmat(mp->nparm, mp->nparm, 0);
  for (i = 0; i < mp->nconf; i++) {
    dtmp = A->map[i];
    for (j = 0; j < mp->nparm; j++) {
      dtj = dtmp[j];
      vec_a[j] += dtj;
      vec_a2[j] += dtj*dtj;
      dtm2p = vec_ab.map[j];
      for (k = j+1; k < mp->nparm; k++) {
        dtm2p[k] += dtj*dtmp[k];
      }
    }
  }

  // Report stupidly high correlations
  for (i = 0; i < mp->nparm; i++) {
    if (vec_a2[i] < 1.0e-12) {
      mp->zerocol[mp->nzerocol] = i;
      mp->nzerocol += 1;
    }
    for (j = i+1; j < mp->nparm; j++) {
      sq = (mp->nconf*vec_a2[i] - vec_a[i]*vec_a[i]) *
        (mp->nconf*vec_a2[j] - vec_a[j]*vec_a[j]);
      if (sq < 1.0e-12) {
        ctest = 0.0;
      }
      else {
        ctest = (mp->nconf*vec_ab.map[i][j] - vec_a[i]*vec_a[j]) / sqrt(sq);
      }
      if (1.0 - fabs(ctest) < 1.0e-4 ||
          (1.0 - fabs(ctest) < 1.0e-1 && UnrelatedDihedrals(mp, i, j) == 1)) {
        mp->corrcol[2*mp->ncorrcol] = i;
        mp->corrcol[2*mp->ncorrcol+1] = j;
        mp->corrval[mp->ncorrcol] = ctest;
        mp->ncorrcol += 1;
      }
    }
  }

  // Free allocated memory
  if (mp->ncmaps > 0) {
    DestroyDmat(&cmtprof);
    DestroyDmat(&cmfprnt);
  }
  DestroyDmat(&vec_ab);
  free(vec_a);
  free(vec_a2);
}

//-----------------------------------------------------------------------------
// WeightFittingMatrix: this function will add weights to the rows of the
//                      fitting matrix for the purpose of raising or lowering
//                      the importance of each conformation.   
//                                                                      
// Arguments:                                                           
//   mp:       fitting control data (contains history and number of fitting
//             points)                                          
//   A:        the pre-allocated fitting matrix                         
//   b:        the pre-allocated array for computing Ax = B             
//-----------------------------------------------------------------------------
static void WeightFittingMatrix(prmset *mp, dmat *A, double *b)
{
  int i, j;
  double sclfac;
  double *dtmp;

  // Bail out if the minimum weight is so  
  // high as to make reweighting pointless.
  if (mp->wtfloor > 100.0) {
    return;
  }

  // Scale each row and the corresponding target value
  for (i = 0; i < mp->nconf; i++) {
    dtmp = A->map[i];
    sclfac = mp->conf[i].wt;
    for (j = 0; j < A->col; j++) {
      dtmp[j] *= sclfac;
    }
    b[i] *= sclfac;
  }
}

//-----------------------------------------------------------------------------
// CountRestraints: function for totaling up the number of restraints on bonds,
//                  angles, and dihedrals.  General restraints on all these
//                  terms may be placed, as well as specific restraints on
//                  individual fitted variables. 
//
// Arguments:                                                           
//   mp:  fitting control data                                          
//-----------------------------------------------------------------------------
static int CountRestraints(prmset *mp)
{
  int i, rcount, nrestraint;

  // Initialize the counter
  nrestraint = 0;

  // Bonds
  for (i = 0; i < mp->nbond; i++) {
    if (mp->bonds[i].fitcolX >= 0 && mp->bonds[i].rstwK > 1.0e-12) {
      nrestraint++;
    }
    if (mp->bonds[i].fitcolY >= 0 && mp->bonds[i].rstwl0 > 1.0e-12) {
      nrestraint++;
    }
  }

  // Angles
  for (i = 0; i < mp->nangl; i++) {
    if (mp->angls[i].fitcolX >= 0 && mp->angls[i].rstwK > 1.0e-12) {
      nrestraint++;
    }
    if (mp->angls[i].fitcolY >= 0 && mp->angls[i].rstwTh0 > 1.0e-12) {
      nrestraint++;
    }
  }

  // Dihedrals
  for (i = 0; i < mp->ntor; i++) {
    if (mp->torsions[i].fitcol >= 0 && mp->torsions[i].rstw > 1.0e-12) {
      nrestraint++;
    }
  }

  // CMAPS
  for (i = 0; i < mp->ncmaps; i++) {
    if (mp->cmaps[i].fitcoli >= 0 && mp->cmaps[i].rstw > 1.0e-12) {
      rcount = mp->cmaps[i].esrfFit.row;
      nrestraint += (rcount * rcount * 2) + 1;
    }
  }

  // 1:4 fitting terms
  if (mp->grst14 > 1.0e-12) nrestraint += mp->fitscnb + mp->fitscee;

  // NMR operations
  for (i = 0; i < mp->nops; i++) {
    if (mp->nmrops[i].fitcol >= 0 && mp->nmrops[i].rstw > 1.0e-12) {
      nrestraint++;
    }
  }

  return nrestraint;
}

//-----------------------------------------------------------------------------
// DestroyPrmset: detroy a prmset structure and all data in it; this may not
//                free data associated with conformations which have been
//                rejected from the fitting set, but those are freed as they
//                are rejected.                           
//                                                                      
// Arguments:                                                           
//   mp:  fitting control data                                          
//-----------------------------------------------------------------------------
static void DestroyPrmset(prmset *mp)
{
  int i;

  free(mp->zerocol);
  free(mp->corrcol);
  free(mp->corrval);
  free(mp->GroupCount);
  free(mp->FirstConf);
  free(mp->LastConf);
  for (i = 0; i < mp->nconf; i++) {
    FreeConformation(&mp->conf[i]);
  }
  for (i = 0; i < mp->natom; i++) {
    free(mp->atoms[i].comment);
  }
  for (i = 0; i < mp->nbond; i++) {
    free(mp->bonds[i].comment);
  }
  for (i = 0; i < mp->nangl; i++) {
    free(mp->angls[i].comment);
  }
  for (i = 0; i < mp->ntor; i++) {
    free(mp->torsions[i].comment);
  }
  for (i = 0; i < mp->nbadj; i++) {
    free(mp->badj[i].comment);
  }
  for (i = 0; i < mp->naadj; i++) {
    free(mp->aadj[i].comment);
  }
  for (i = 0; i < mp->nhadj; i++) {
    free(mp->hadj[i].comment);
  }
  for (i = 0; i < mp->nunisys; i++) {
    FreeTopology(&mp->tpencyc[i]);
  }
  for (i = 0; i < mp->nspectralH; i++) {
    DestroyCmat(&mp->spectralH[i].atmtypes);
  }
  if (mp->ncmaps > 0) {
    for (i = 0; i < mp->ncmaps; i++) {
      DestroyXCmapDef(&mp->cmaps[i]);
    }
    free(mp->cmaps);
  }
  free(mp->recast);
  free(mp->cleave);
  free(mp->conf);
  free(mp->badj);
  free(mp->aadj);
  free(mp->hadj);
  free(mp->atoms);
  free(mp->bonds);
  free(mp->angls);
  free(mp->torsions);
  free(mp->tpencyc);
  free(mp->ititl);
  free(mp->icomm);
  free(mp->NMROpsFile);
  free(mp->NMRParmFile);
  free(mp->spectralH);
  for (i = 0; i < mp->nunisys; i++) {
    free(mp->BondMapIDBucket[i]);
    free(mp->AnglMapIDBucket[i]);
    free(mp->DiheMapIDBucket[i]);
    free(mp->NOpsMapIDBucket[i]);
    free(mp->SplnMapIDBucket[i]);
    DestroyDmat(&mp->BondMapValBucket[i]);
    DestroyDmat(&mp->BondMapUxBucket[i]);
    DestroyDmat(&mp->BondMapUyBucket[i]);
    DestroyDmat(&mp->BondMapUzBucket[i]);
    DestroyDmat(&mp->BondMapUwBucket[i]);
    DestroyDmat(&mp->BondMapCntrbBucket[i]);
    DestroyDmat(&mp->AnglMapValBucket[i]);
    DestroyDmat(&mp->AnglMapUxBucket[i]);
    DestroyDmat(&mp->AnglMapUyBucket[i]);
    DestroyDmat(&mp->AnglMapCntrbBucket[i]);
    DestroyDmat(&mp->DiheMapValBucket[i]);
    DestroyDmat(&mp->DiheMapUBucket[i]);
    DestroyDmat(&mp->DiheMapCntrbBucket[i]);
    DestroyDmat(&mp->NOpsMapValBucket[i]);
    DestroyDmat(&mp->NOpsMapKernelBucket[i]);
    DestroyDmat(&mp->NOpsMapCntrbBucket[i]);
    DestroyDmat(&mp->SplnMapValBucket[i]);
    DestroyDmat(&mp->SplnMapCntrbBucket[i]);
  }
  free(mp->BondMapIDBucket);
  free(mp->AnglMapIDBucket);
  free(mp->DiheMapIDBucket);
  free(mp->NOpsMapIDBucket);
  free(mp->SplnMapIDBucket);
  free(mp->BondMapValBucket);
  free(mp->BondMapUxBucket);
  free(mp->BondMapUyBucket);
  free(mp->BondMapUzBucket);
  free(mp->BondMapUwBucket);
  free(mp->BondMapCntrbBucket);
  free(mp->AnglMapValBucket);
  free(mp->AnglMapUxBucket);
  free(mp->AnglMapUyBucket);
  free(mp->AnglMapCntrbBucket);
  free(mp->DiheMapValBucket);
  free(mp->DiheMapUBucket);
  free(mp->DiheMapCntrbBucket);
  free(mp->NOpsMapValBucket);
  free(mp->NOpsMapKernelBucket);
  free(mp->NOpsMapCntrbBucket);
  free(mp->SplnMapValBucket);
  free(mp->SplnMapCntrbBucket);
}

//-----------------------------------------------------------------------------
// GeometryConstraints: this is a function to wrap a lot of iterative solution
//                      of Ax = b, but subject to more restraint forces.  It
//                      is a lot like spectral resampling, but different enough
//                      to warrant a new function.  Here we will push and pull
//                      on selected groups of angles until we add or subtract
//                      just enough to get the optimal sum of all of them.
//
// Arguments:                                                           
//   mp:     fitting control data, contains topology and coordinate file
//           names, as well as keys for correlating parameters across   
//           systems.  Also contains the frcmod output file name.       
//   A :     original coefficient matrix for bonded parameter fitting   
//   b :     original target vector for bonded parameter fitting        
//   x :     global optimum solution vector to Ax = b                   
//-----------------------------------------------------------------------------
static void GeometryConstraints(prmset *mp, dmat *A, double* b, double* x)
{
  int h, i, j, pivot, nvar, fcX, fcY;
  double bssX, bssY, distrib, range;
  double* basecombo;
  double* currcombo;
  double* bwork;
  double* btest;
  dmat Awork;

  // Alloate scratch memory
  Awork = CreateDmat(A->row, A->col, 0);
  bwork = (double*)malloc(A->row*sizeof(double));
  btest = (double*)malloc(A->row*sizeof(double));

  // Loop over all geometry constraints
  for (h = 0; h < mp->nanglsum; h++) {

    // The angles associated with this constraint have already
    // been found, as have the rows of the restraints.  Decide
    // on a combination of the angle equilibria and set the   
    // restraints accordingly.                                
    nvar = mp->anglsum[h].nvar;
    range = mp->anglsum[h].srange;
    basecombo = (double*)malloc(nvar*sizeof(double));
    for (i = 0; i < nvar; i++) {
      fcX = mp->anglsum[h].fitcol[2*i];
      fcY = mp->anglsum[h].fitcol[2*i+1];
      bssX = mp->anglsum[h].basis[2*i];
      bssY = mp->anglsum[h].basis[2*i+1];
      basecombo[i] = (x[fcX]*bssX + x[fcY]*bssY) / (x[fcX] + x[fcY]);

      // Grab the restraint row for this angle,
      // now that it is known.                 
      for (j = 0; j < mp->nangl; j++) {
        if (mp->angls[j].fitcolX == mp->anglsum[h].fitcol[2*i]) {
          mp->anglsum[h].rstrow[i] = mp->angls[j].restTh0;
        }
      }
    }

    // Satisfy the angle sum constraint.
    distrib = mp->anglsum[h].target - DSum(basecombo, nvar);
    distrib /= nvar;
    for (i = 0; i < mp->anglsum[h].nvar; i++) {
      basecombo[i] += distrib;
    }

    // Roll through all variable angles, adjust them by up to
    // five degrees from their base values (the base values  
    // are the equilibria found in the optimal solution plus 
    // whatever adjustment satisfies the constraint).        
    // Re-evaluate the entire solution for all combinations  
    // in which all angles are within five degrees of their  
    // base values.                                          
    pivot = 0;
    currcombo = (double*)malloc(nvar*sizeof(double));
    for (i = 0; i < nvar; i++) {
      currcombo[i] = basecombo[i] - range;
    }
    if (mp->verbose == 1) {
      printf("mdgx >> Evaluating geometry constraint %d.\n", h);
    }
    while(pivot < nvar-1) {

      // Set the angle variable in the pivot 
      // and increment the pivot if necessary
      currcombo[pivot] += 0.25*PI/180.0;
      while (pivot < nvar-1 &&
             currcombo[pivot] > basecombo[pivot] + range) {
        currcombo[pivot] = basecombo[pivot] - range;
        pivot++;
        currcombo[pivot] += 0.25*PI/180.0;
      }

      // The final angle variable satisfies the constraint
      currcombo[nvar-1] = mp->anglsum[h].target -
        DSum(currcombo, nvar-1);

      // Check to see whether all angles
      // are within the correct range.  
      if (currcombo[nvar-1] < basecombo[nvar-1] - range ||
          currcombo[nvar-1] > basecombo[nvar-1] + range) {
        continue;
      }

      // If we are still here, change the restraints
      // and evaluate a new solution.               
      printf("We are still here with: ");
      for (i = 0; i < nvar; i++) {
        printf("%7.2lf ", currcombo[i]*180.0/PI);
      }
      printf("\n");
      CopyDmat(&Awork, A, 1);
      ReflectDVec(bwork, b, A->row);
      for (i = 0; i < nvar; i++) {
        fcX = mp->anglsum[h].fitcol[2*i];
        fcY = mp->anglsum[h].fitcol[2*i+1];
        bssX = mp->anglsum[h].basis[2*i];
        bssY = mp->anglsum[h].basis[2*i+1];
        j = mp->anglsum[h].rstrow[i];
        Awork.map[j][fcX] = 0.1*A->row*(currcombo[i] - bssX);
        Awork.map[j][fcY] = 0.1*A->row*(currcombo[i] - bssY);
      }
      AxbQRRxc(Awork, bwork, mp->verbose);
      BackSub(Awork, bwork);
      printf("  -> Solution: ");
      for (i = 0; i < nvar; i++) {
        fcX = mp->anglsum[h].fitcol[2*i];
        fcY = mp->anglsum[h].fitcol[2*i+1];
        bssX = mp->anglsum[h].basis[2*i];
        bssY = mp->anglsum[h].basis[2*i+1];
        printf("%7.2lf %7.2lf   ", (bwork[fcX]*bssX + bwork[fcY]*bssY) /
               (bwork[fcX] + bwork[fcY]) * 180.0/PI, bwork[fcX] + bwork[fcY]);
      }

      DMatVecMult(A, bwork, btest);
      printf("               :: %9.4lf\n", VecRMSD(btest, b, mp->nconf));

    }

    // Free allocated memory
    free(basecombo);
    free(currcombo);
  }

  // Free allocated memory
  DestroyDmat(&Awork);
  free(bwork);
  free(btest);
}

//-----------------------------------------------------------------------------
// InitSpectrumControl: this function will initialize a spectrum control
//                      structure.  See the SpectrumControl structure   
//                      declaration in ParamFitDS.h for descriptions of the
//                      eponymous input variables.                  
//-----------------------------------------------------------------------------
static void InitSpectrumControl(speccon *sc, int level, int order, int colX,
                                int colY, double llimX, double hlimX,
                                double llimY, double hlimY, double gspcX,
                                double gspcY)
{
  sc->level = level;
  sc->order = order;
  sc->colXmain = colX;
  sc->colYmain = colY;
  sc->llimX = llimX;
  sc->hlimX = hlimX;
  sc->llimY = llimY;
  sc->hlimY = hlimY;
  sc->gspcX = gspcX;
  sc->gspcY = gspcY;
}

//-----------------------------------------------------------------------------
// MarkSpectrumVariables: this function takes a spectrum request and finds
//                        appropriate dihedrals to meet it.       
//                                                                      
// Arguments:                                                           
//   mp:     fitting control data, contains topology and coordinate file names,
//           as well as keys for correlating parameters across systems.  Also
//           contains the frcmod output file name.       
//   sr:     spectrum request                                           
//   marks:  vector of SpectrumControl structs to enumerate dihedrals for the 
//           spectrum (allocated when it comes in, initialized in this
//           function)                                          
//-----------------------------------------------------------------------------
static void MarkSpectrumVariables(prmset *mp, specreq *sr, speccon* marks)
{
  int i, j, match;

  if (sr->order == 4) {
    for (i = 0; i < mp->ntor; i++) {

      // Continue if this is not an adjustable dihedral
      if (mp->torsions[i].fitcol < 0) {
        continue;
      }

      // Check to see that all types enumerated in the
      // spectrum request are found in the dihedral.  
      match = 0;
      if (sr->atmtypes.row == 1) {
        if (str4cmp(sr->atmtypes.map[0], mp->torsions[i].atype) == 0 ||
            str4cmp(sr->atmtypes.map[0], mp->torsions[i].btype) == 0 ||
            str4cmp(sr->atmtypes.map[0], mp->torsions[i].ctype) == 0 ||
            str4cmp(sr->atmtypes.map[0], mp->torsions[i].dtype) == 0) {
          match = 1;
        }
      }
      if (sr->atmtypes.row == 2) {
        if (TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
                        mp->torsions[i].ctype, mp->torsions[i].dtype,
                        sr->atmtypes.map[0], sr->atmtypes.map[1], "X   ",
                        "X   ", sr->order, 0) > 0 ||
            TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
                        mp->torsions[i].ctype, mp->torsions[i].dtype, "X   ",
                        sr->atmtypes.map[0], sr->atmtypes.map[1], "X   ",
                        sr->order, 0) > 0 ||
            TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
                        mp->torsions[i].ctype, mp->torsions[i].dtype,
                        "X   ", "X   ", sr->atmtypes.map[0],
                        sr->atmtypes.map[1], sr->order, 0) > 0 ||
            TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
                        mp->torsions[i].ctype, mp->torsions[i].dtype,
                        sr->atmtypes.map[1], sr->atmtypes.map[0], "X   ",
                        "X   ", sr->order, 0) > 0 ||
            TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
                        mp->torsions[i].ctype, mp->torsions[i].dtype,
                        "X   ", sr->atmtypes.map[1], sr->atmtypes.map[0],
                        "X   ", sr->order, 0) > 0 ||
            TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
                        mp->torsions[i].ctype, mp->torsions[i].dtype, "X   ",
                        "X   ", sr->atmtypes.map[1], sr->atmtypes.map[0],
                        sr->order, 0) > 0) {
          match = 1;
        }
      }
      else if (sr->atmtypes.row == 3) {
        if (TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
                        mp->torsions[i].ctype, mp->torsions[i].dtype,
                        sr->atmtypes.map[0], sr->atmtypes.map[1],
                        sr->atmtypes.map[2], "X   ", sr->order, 0) > 0 ||
            TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
                        mp->torsions[i].ctype, mp->torsions[i].dtype,
                        "X   ", sr->atmtypes.map[0], sr->atmtypes.map[1],
                        sr->atmtypes.map[2], sr->order, 0) > 0 ||
            TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
                        mp->torsions[i].ctype, mp->torsions[i].dtype,
                        sr->atmtypes.map[2], sr->atmtypes.map[1],
                        sr->atmtypes.map[0], "X   ", sr->order, 0) > 0 ||
            TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
                        mp->torsions[i].ctype, mp->torsions[i].dtype,
                        "X   ", sr->atmtypes.map[2], sr->atmtypes.map[1],
                        sr->atmtypes.map[0], sr->order, 0) > 0) {
          match = 1;
        }
      }
      else if (sr->atmtypes.row == 4) {
        if (TypeCompare(mp->torsions[i].atype, mp->torsions[i].btype,
                        mp->torsions[i].ctype, mp->torsions[i].dtype,
                        sr->atmtypes.map[0], sr->atmtypes.map[1],
                        sr->atmtypes.map[2], sr->atmtypes.map[3],
                        sr->order, 1) > 0) {
          match = 1;
        }
      }
      if (match == 0) {
        continue;
      }

      // This torsion meets the criteria for spectral reoptimization.
      j = mp->torsions[i].fitcol;
      if (marks[j].level < sr->level) {
        InitSpectrumControl(&marks[j], sr->level, sr->order, j, -1, sr->minval,
                            sr->maxval, 0.0, 0.0, sr->gspc, 0.1);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// FindRestraintEquations: a function for filling out the rest of the
//                         SpectrumControl data struct by finding which rows of
//                         each matrix appear to serve as restraint equations
//                         for particular columns.
//
// Arguments:                                                           
//   A    : the main or reduced fitting matrix                          
//   sc   : the spectrum control structure                              
//   red  : flag to indicate whether this is the main (set to 0) or     
//          reduced matrix (set to 1)                                   
//-----------------------------------------------------------------------------
static void FindRestraintEquations(dmat *A, speccon *sc, int red)
{
  int i, j, rcandidate, colX, colY, rstX, rstY;
  double *dtmp;

  // Bail out if there's nothing to find
  if (sc->order == 0) {
    if (red == 0) {
      sc->rstXmain = -1;
      sc->rstYmain = -1;
    }
    else {
      sc->rstXred = -1;
      sc->rstYred = -1;
    }
    return;
  }

  // Initialize the rows of restraint equations to -1,
  // and initialize the columns of the variables based
  // on which matrix we are dealing with.             
  rstX = -1;
  rstY = -1;
  if (red == 0) {
    colX = sc->colXmain;
    colY = sc->colYmain;
  }
  else {
    colX = sc->colXred;
    colY = sc->colYred;
  }

  // Find restraint equations in the matrix
  for (i = A->row-1; i >= 0; i--) {
    dtmp = A->map[i];
    rcandidate = 1;
    for (j = 0; j < A->col; j++) {
      if (fabs(dtmp[j]) > 1.0e-12 && j != colX && j != colY) {
        rcandidate = 0;
        break;
      }
    }
    if (rcandidate == 1) {
      if (sc->order == 4 && fabs(dtmp[colX]) > 1.0e-12) {        
        rstX = i;
      }
      if ((sc->order == 3 || sc->order == 2) &&
          fabs(dtmp[colX]) > 1.0e-12) {
        if (sc->rstYmain < 0) {
          rstY = i;
        }
        else {
          rstX = i;
        }
      }
    }
  }

  // Store the restraint equation row indices
  if (red == 0) {
    sc->rstXmain = rstX;
    sc->rstYmain = rstY;
  }
  else {
    sc->rstXred = rstX;
    sc->rstYred = rstY;
  }
}

//-----------------------------------------------------------------------------
// ComputeSysRMSE: a function for computing the root mean squared error for
//                 each unique system, given a large set of data in which 
//                 conformations of separate systems may be interspersed.
//
// Arguments:                                                           
//   sysrmse  : vector of doubles, as many indices as there are unique systems.
//              This vector will be initialized within this function.
//   sysidx   : map of system indices, a vector of integers             
//   target   : the target energies of each conformation in the data set
//   guess    : the energy guesses for each conformation                
//   mp       : fitting control data, contains counts of conformations (in
//              total and on a per-system basis) as well as the number of
//              unique systems                                
//-----------------------------------------------------------------------------
static void ComputeSysRMSE(double* sysrmse, int* sysidx, double* target,
                           double* guess, prmset *mp)
{
  int i;
  double dv;

  SetDVec(sysrmse, mp->nunisys, 0.0);
  for (i = 0; i < mp->nconf; i++) {
    dv = guess[i] - target[i];
    sysrmse[sysidx[i]] += dv*dv;
  }
  for (i = 0; i < mp->nunisys; i++) {
    sysrmse[i] = sqrt(sysrmse[i] / mp->GroupCount[i]);
  }
}

//-----------------------------------------------------------------------------
// SpectrumSolution: a function for solving a reduced form of the matrix
//                   equation many times, with varying restraints on selected
//                   parameters to test the quality of alternate solutions.
//                   Alternate solutions which come within a tolerance of the
//                   normally restrained optimum will be sorted, and
//                   sufficiently divergent alternatives will be reported in
//                   separate files.  This function finds the alternative
//                   solutions for evaluation later.
//
// Arguments:                                                           
//   A:     a copy of the original, unweighted fitting matrix.  An extra copy
//          of A will be made so that at least one pristine copy of A can be
//          preserved for scoring and posterior evaluations.   
//   b:     a copy of the original target vector.  Like the fitting matrix
//          itself, another copy will be produced so that the first copy can
//          be saved for posterior evaluations.          
//   x:     the original solution vector, essential for subtracting the 
//          contributions of fixed, solved parameters when generating   
//          alternative solutions.                                      
//   mp:    fitting control data, contains topology and coordinate file names,
//          as well as keys for correlating parameters across systems.  Also
//          contains the frcmod output file name.        
//   nsoln: the number of alternate solutions produced (an array of them will
//          be returned, but the count goes back in the form of a pointer to
//          integer)                                         
//-----------------------------------------------------------------------------
static sresamp SpectrumSolution(dmat *A, double* b, double* x, prmset *mp)
{
  int i, j, ired, jred, nredcol, nredrow, nredconf, redo;
  int maxsoln, isoln, firstpass, passorder;
  int* retaincol;
  int* retainrow;
  int* ninst;
  int* sysidx;
  double rx, ry, rweight, pinwt;
  double *dtmp;
  double* xc;
  double* btest;
  double* bredtest;
  double* bred;
  double* bredc;
  double* sysrmse;
  double* sysrmsetest;
  dmat Ared, Aredc;
  speccon* spv;
  sresamp ms;

  // If spectrum resampling is not requested, just fill in
  // the necessary components of the resampling structure 
  // and bail right out.                                  
  ms.A = A;
  ms.x = x;
  ms.b = b;
  if (mp->spectrum == 0) {
    ms.naltsoln = 0;
    ms.specvar = (speccon*)calloc(1, sizeof(speccon));
    return ms;
  }

  // Reduce the fitting matrix to its bare essentials, the
  // columns containing variables we want to manipulate.  
  // Initializing spv with calloc sets all of its         
  // fields to zero.                                      
  spv = (speccon*)calloc(A->col, sizeof(speccon));
  retaincol = (int*)malloc(A->col*sizeof(int));
  for (i = 0; i < mp->nspectralH; i++) {
    MarkSpectrumVariables(mp, &mp->spectralH[i], spv);
  }
  for (i = mp->nparm; i < A->col; i++) {
    InitSpectrumControl(&spv[i], 1, 0, i, -1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    spv[i].colXmain = i;
    spv[i].colYmain = -1;
  }
  nredcol = 0;
  for (i = 0; i < A->col; i++) {
    retaincol[i] = -1;
    if (spv[i].level > 0) {
      retaincol[i] = nredcol;
      nredcol++;
    }
  }

  // The mapping of the main to the reduced matrix stored in
  // retaincol now needs to be imprinted on the spectrum    
  // control variables.                                     
  for (i = 0; i < A->col; i++) {
    if (spv[i].level > 0) {
      if (spv[i].colXmain >= 0) {
        spv[i].colXred = retaincol[spv[i].colXmain];
      }
      if (spv[i].colYmain >= 0) {
        spv[i].colYred = retaincol[spv[i].colYmain];
      }
    }
  }

  // Create the reduced matrix and solution vector.
  // Keep a tally of how many actual conformations 
  // are part of the reduced matrix.  Make a map of
  // the systems to which each row corresponds.    
  retainrow = (int*)calloc(A->row, sizeof(int));
  sysidx = (int*)malloc(A->row*sizeof(int));
  nredrow = 0;
  nredconf = 0;
  for (i = 0; i < A->row; i++) {
    retainrow[i] = -1;
    dtmp = A->map[i];
    for (j = 0; j < A->col; j++) {
      if (spv[j].level > 0 && fabs(dtmp[j]) > 1.0e-12) {
        retainrow[i] = nredrow;
        if (i < mp->nconf) {
          nredconf++;
        }
        nredrow++;
        break;
      }
    }
    for (j = mp->nparm; j < A->col; j++) {
      if (dtmp[j] > 1.0e-1) {
        sysidx[i] = j - mp->nparm;
      }
    }
  }
  Ared = CreateDmat(nredrow, nredcol, 0);
  bred = (double*)calloc(nredrow, sizeof(double));
  ired = 0;
  for (i = 0; i < A->row; i++) {
    if (retainrow[i] == -1) {
      continue;
    }
    dtmp = A->map[i];
    jred = 0;
    bred[ired] = b[i];
    for (j = 0; j < A->col; j++) {
      if (spv[j].level == 0) {
        bred[ired] -= dtmp[j]*x[j];
      }
      else {
        Ared.map[ired][jred] = dtmp[j];
        jred++;
      }
    }
    ired++;
  }

  // Rather than add complexity to the parameter structures
  // and try to track the restraint equations, we will just
  // seek out patterns in the reduced fitting matrix which 
  // look like restraint equations for each parameter.     
  for (i = 0; i < A->col; i++) {
    FindRestraintEquations(A, &spv[i], 0);
    FindRestraintEquations(&Ared, &spv[i], 1);
  }

  // Allocate for reduced solutions, assess the global
  // minimum solution, and get an approximate count of
  // instances so that new restraints can be designed.
  btest = (double*)malloc(A->row*sizeof(double));
  ninst = (int*)malloc(Ared.col*sizeof(int));
  for (i = 0; i < Ared.col; i++) {
    ninst[i] = CountInstances(&Ared, nredconf, i);
  }
  DMatVecMult(A, x, btest);
  sysrmse = (double*)malloc(mp->nunisys*sizeof(double));
  sysrmsetest = (double*)malloc(mp->nunisys*sizeof(double));
  ComputeSysRMSE(sysrmse, sysidx, b, btest, mp);

  // Allocate space for storing the solutions we are about
  // to create.  This space may need to be expanded later.
  // The first "alternate" solution is the global optimum.
  maxsoln = 32;
  ms.xalt = (altsol*)malloc(maxsoln*sizeof(altsol));
  ms.xalt[0].ndim = A->col;
  ms.xalt[0].x = CpyDVec(x, A->col);
  ms.xalt[0].b = (double*)malloc(A->row*sizeof(double));
  DMatVecMult(A, x, ms.xalt[0].b);
  ms.xalt[0].rstpos = CreateImat(1, 2);
  ms.xalt[0].rstval = CreateDmat(1, 2, 0);
  ms.xalt[0].rstpos.map[0][0] = -1;
  isoln = 1;

  // Loop over the spectrum variables, selecting one or two at 
  // a time.  Reduce the matrix problem using the original     
  // solution vector to fill in settings for variables outside 
  // of the immediate reoptimization.  Energy offset columns   
  // corresponding to each system will be retained from the    
  // original matrix.  Ghost rows, those containing all zeros  
  // after removing unnecessary columns, will be eliminated.   
  // Within an inner loop or loops, constrain each spectrum    
  // variable to a specific value, solve the reduced problem,  
  // and store the unique solution.                            
  passorder = -1;
  for (i = 0; i < A->col; i++) {
    if (spv[i].level != 2) {
      continue;
    }

    // We are now about to resample a parameter.  In order to do
    // so most efficiently, we need to decide on a restraint    
    // weight.  If this is a new class of parameters (bonds,    
    // angles, dihedrals), then we need to determine the right  
    // strength of restraint, and maybe adjust it in flight.    
    // When encountering a new class of parameter, start at     
    // 10.0x the typical value, adjust upwards by 20% increments
    // if this isn't strong enough, then adjust up or down by 5%
    // increments.                                              
    if (spv[i].order != passorder) {
      firstpass = 1;
      passorder = spv[i].order;
      pinwt = 50.0;
    }
    if (pinwt > 100.0) {
      pinwt = 100.0;
    }

    // Start with a sensible guess for how tightly we must restrain
    // the target parameter.  This restraint multiplier will evolve
    // over the range of trial values for this target parameter.   
    for (rx = spv[i].llimX; rx <= spv[i].hlimX; rx += spv[i].gspcX) {
      for (ry  = spv[i].llimY; ry <= spv[i].hlimY; ry += spv[i].gspcY) {
        redo = 1;
        while (redo == 1) {

          // Make the reduced matrix and target, with new restraints
          CopyDmat(&Aredc, &Ared, 0);
          bredc = CpyDVec(bred, Ared.row);

            // FIX ME!!!  Need to implement for bond and angle sampling,
          // and that will take some more bookkeeping.
          if (spv[i].order == 4 && spv[i].rstXred >= 0) {
            rweight = ninst[spv[i].colXred] * mp->grstH * pinwt;
            Aredc.map[spv[i].rstXred][spv[i].colXred] = rweight;
            bredc[spv[i].rstXred] = rx * rweight;
          }
          AxbQRRxc(Aredc, bredc, 0);
          BackSub(Aredc, bredc);

          // Check if the parameter has stayed within its
          // predefined range.  If not, ratchet up the   
          // restraint and try again.                    
          if (spv[i].order == 4 && spv[i].rstXred >= 0) {
            if (fabs(bredc[spv[i].colXred] - rx) < spv[i].gspcX) {
              redo = 0;
              pinwt *= 0.95;
            }
            else {
              if (firstpass == 1) {
                pinwt /= 0.80;
              }
              else {
                pinwt /= 0.95;
              }
              DestroyDmat(&Aredc);
              free(bredc);
            }
          }
        }

        // Success--there is now a new solution to the reduced  
        // system of equations which satisfies the conditions   
        // placed on the parameter of interest.  Map the reduced
        // solution back to the complete system of equations    
        // and compute its approximate score across all systems 
        // in the training set.                                 
        xc = CpyDVec(x, A->col);
        for (j = 0; j < A->col; j++) {
          if (spv[j].level > 0) {
            if (spv[j].order == 2 || spv[j].order == 3) {
              if (spv[j].colXmain >= 0 && spv[j].colXred >= 0) {
                xc[spv[j].colXmain] = bredc[spv[j].colXred];
              }
              if (spv[j].colYmain >= 0 && spv[j].colYred >= 0) {
                xc[spv[j].colYmain] = bredc[spv[j].colYred];
              }
            }
            else if (spv[j].order == 4 || spv[j].order == 0) {
              if (spv[j].colXmain >= 0 && spv[j].colXred >= 0) {
                      xc[spv[j].colXmain] = bredc[spv[j].colXred];
              }
            }
          }
        }
        bredtest = (double*)malloc(A->row*sizeof(double));
        DMatVecMult(A, xc, bredtest);
        ComputeSysRMSE(sysrmsetest, sysidx, b, bredtest, mp);
        for (j = 0; j < mp->nunisys; j++) {
          sysrmsetest[j] -= sysrmse[j] * mp->spvtol;
        }

        // If the sub-optimal solution is sufficiently accurate
        // across all systems in the training set, store it for
        // filtering later, along with a compact statement of  
        // the particular restraints that created it.          
        if (DExtreme(sysrmsetest, mp->nunisys, 1) < 0.0) {
          ms.xalt[isoln].ndim = A->col;
          ms.xalt[isoln].x = xc;
          ms.xalt[isoln].b = bredtest;

          // FIX ME!!!  Need to add support for bonds and angles,
          // mirroring what is needed above.
          if (spv[i].order == 4) {
            ms.xalt[isoln].rstpos = CreateImat(1, 2);
            ms.xalt[isoln].rstval = CreateDmat(1, 2, 0);
            if (spv[i].rstXred >= 0) {
              ms.xalt[isoln].rstpos.map[0][0] = spv[i].rstXmain;
              ms.xalt[isoln].rstpos.map[0][1] = spv[i].colXmain;
              ms.xalt[isoln].rstval.map[0][0] = rweight;
              ms.xalt[isoln].rstval.map[0][1] = rx * rweight;
            }
            else {
              ms.xalt[isoln].rstpos.map[0][0] = -1;
            }
          }
          isoln += 1;
          if (isoln == maxsoln) {
            maxsoln += 32;
            ms.xalt = (altsol*)realloc(ms.xalt, maxsoln*sizeof(altsol));
          }
        }
        else {
          free(xc);
          free(bredtest);
        }
        DestroyDmat(&Aredc);
        free(bredc);
      }
    }

    // Whatever we were working on, the first pass of it is complete.
    firstpass = 0;
  }

  // Free allocated memory
  DestroyDmat(&Ared);
  free(btest);
  free(bred);
  free(retainrow);
  free(retaincol);
  free(ninst);
  free(sysrmse);
  free(sysidx);
  free(sysrmsetest);

  // Return alternate solutions as part of a resampling structure
  ms.naltsoln = isoln;
  ms.specvar = spv;

  return ms;
}

//-----------------------------------------------------------------------------
// SpectrumEvaluation: evaluate alternative solutions to the (reduced) least
//                     squares problem in terms of their relationships to one
//                     another and to the optimum solution.  Those with
//                     sufficient differences are recomputed according to the
//                     full least squares problem.  The performance of each
//                     such solution is evaluated in detail.
//
// Arguments:                                                           
//   mp:  fitting control data, contains topology and coordinate file   
//        names, as well as keys for correlating parameters across      
//        systems.  Also contains the frcmod output file name.          
//-----------------------------------------------------------------------------
static void SpectrumEvaluation(prmset *mp, sresamp *ms)
{
  int i, j, k, qualcount, bestloc, skip, ninst;
  int* qualified;
  double rmsd, bestrmsd;
  double* bcpy;
  double* btest;
  imat *rstp;
  dmat ssd, Acpy;
  dmat *rstv;
  altsol tmpalt;

  // Allocate space to store the indices of
  // the most mutually distant solutions   
  qualified = (int*)calloc(mp->nspectralx, sizeof(int));

  // Step through all solutions and compute mutual parameter RMSDs.
  // This will get filtered into one of the matrices stored in ms. 
  ssd = CreateDmat(ms->naltsoln, ms->naltsoln, 0);
  for (i = 0; i < ms->naltsoln; i++) {
    for (j = i+1; j < ms->naltsoln; j++) {
      ssd.map[i][j] = VecRMSD(ms->xalt[i].b, ms->xalt[j].b, mp->nconf);
      ssd.map[j][i] = ssd.map[i][j];
    }
  }

  // The first qualifying solution is 0, the global optimum.
  // Find a number of additional solutions based on maximum 
  // mutual distance.  Find a total of nqual solutions.     
  for (qualcount = 1; qualcount < mp->nspectralx; qualcount++) {

    // Loop over all solutions, skipping those already in
    // the list, and assess the mutual distance to any   
    // that are in the list.                             
    bestrmsd = 0.0;
    for (j = 0; j < ms->naltsoln; j++) {
      skip = 0;
      for (k = 0; k < qualcount; k++) {
        if (qualified[k] == j) {
          skip = 1;
        }
      }
      if (skip) {
        continue;
      }
      rmsd = 0.0;
      for (k = 0; k < qualcount; k++) {
        rmsd += ssd.map[j][k] * ssd.map[j][k];
      }
      rmsd = sqrt(rmsd/qualcount);
      if (rmsd > bestrmsd) {
        bestrmsd = rmsd;
        bestloc = j;
      }
    }
    qualified[qualcount] = bestloc;
  }

  // With some of the most distant solutions selected, 
  // solve them again, this time using the full matrix.
  // Store the way in which the initial solutions      
  // relate to one another for future comparison.      
  ms->ssd1 = CreateDmat(mp->nspectralx, mp->nspectralx, 0);
  for (i = 0; i < mp->nspectralx; i++) {
    for (j = i+1; j < mp->nspectralx; j++) {
      ms->ssd1.map[i][j] = ssd.map[qualified[i]][qualified[j]];
      ms->ssd1.map[j][i] = ms->ssd1.map[i][j];
    }
  }
  for (i = 0; i < mp->nspectralx; i++) {
    CopyDmat(&Acpy, ms->A, 0);
    bcpy = CpyDVec(ms->b, ms->A->row);
    rstp = &ms->xalt[qualified[i]].rstpos;
    rstv = &ms->xalt[qualified[i]].rstval;
    for (j = 0; j < rstp->row; j++) {
      if (rstp->map[j][0] >= 0) {
        for (k = 1; k < rstp->col; k++) {
          Acpy.map[rstp->map[j][0]][rstp->map[j][k]] = rstv->map[j][k-1];
        }
        bcpy[rstp->map[j][0]] = rstv->map[j][rstp->col-1];
      }
    }
    AxbQRRxc(Acpy, bcpy, 0);
    BackSub(Acpy, bcpy);
    btest = (double*)malloc(ms->A->row*sizeof(double));
    DMatVecMult(ms->A, bcpy, btest);
    ReflectDVec(ms->xalt[qualified[i]].x, bcpy, ms->A->col);
    ReflectDVec(ms->xalt[qualified[i]].b, btest, ms->A->col);
    DestroyDmat(&Acpy);
    free(bcpy);
    free(btest);
  }

  // With the new and maximally distant solutions
  // known, the rest can now be discarded.       
  for (i = 0; i < mp->nspectralx; i++) {
    SWAP(ms->xalt[i], ms->xalt[qualified[i]], tmpalt);
  }
  for (i = mp->nspectralx; i < ms->naltsoln; i++) {
    free(ms->xalt[i].b);
    free(ms->xalt[i].x);
    DestroyImat(&ms->xalt[i].rstpos);
    DestroyDmat(&ms->xalt[i].rstval);
  }
  ms->xalt = (altsol*)realloc(ms->xalt, mp->nspectralx*sizeof(altsol));
  ms->ssd2 = CreateDmat(mp->nspectralx, mp->nspectralx, 0);
  for (i = 0; i < mp->nspectralx; i++) {
    for (j = i+1; j < mp->nspectralx; j++) {
      ms->ssd2.map[i][j] = VecRMSD(ms->xalt[i].b, ms->xalt[j].b, mp->nconf);
      ms->ssd2.map[j][i] = ms->ssd2.map[i][j];
    }
  }

  // The restraint values include factors of the number of    
  // instances in which the variable being restrained appears.
  for (i = 0; i < mp->nspectralx; i++) {
    rstp = &ms->xalt[i].rstpos;
    rstv = &ms->xalt[i].rstval;
    for (j = 0; j < rstp->row; j++) {
      if (rstp->map[j][0] >= 0) {
        ninst = CountInstances(ms->A, mp->nconf, rstp->map[j][1]);
        for (k = 1; k <= rstp->col; k++) {
          rstv->map[j][k-1] /= ninst;
        }
      }
    }
  }

  // Free allocated memory
  DestroyDmat(&ssd);
  free(qualified);
}

//-----------------------------------------------------------------------------
// DetectIPolQData: run through the conformations and check to see whether any
//                  of them invoke IPolQ energy and coordinate bundles.  If
//                  any of them do, then all of them must.  Set a flag to
//                  indicate that IPolQ data fitting is in effect, which will
//                  activate additional calculations when it comes time to
//                  compute the potential energy surfaces.  Exit if only part
//                  of the data is energy and coordinate bundles.
//
// Arguments:
//   mp:  fitting control data, contains topology and coordinate file names, as
//        well as keys for correlating parameters across systems.
//-----------------------------------------------------------------------------
static int DetectIPolQData(prmset *mp)
{
  int i, inptype, fitype, hasIpq, hasNonIpq;

  hasIpq = 0;
  hasNonIpq = 0;
  for (i = 0; i < mp->nconf; i++) {
    inptype = FindFileType(mp->conf[i].crdsrc);
    if (inptype == 0) {
      fitype = ValidCoord(mp->conf[i].crdsrc);
      if (fitype == 6) {
        hasIpq = 1;
      }
      else {
        hasNonIpq = 1;
      }
    }
    else if (inptype == 1 || inptype == 2) {

      // The only way to get directories or regular expressions
      // into the fitting data is to have them be IPolQ bundles.
      hasIpq = 1;
    }
  }
  if (hasIpq == 1 && hasNonIpq == 1) {
    printf("DetectIPolQData >> Error.  A mixture of files containing IPolQ "
           "energy and\nDetectIPolQData >> coordinate bundles has been "
           "detected, and this is not\nDetectIPolQData >> valid.\n");
    exit(1);
  }

  // If the IPolQ data fitting framework is in effect, say so
  return hasIpq;
}

//-----------------------------------------------------------------------------
// FitParams: the main function for optimizing parameters to fit a given
//            potential energy surface.                                 
//                                                                      
// Arguments:                                                           
//   mp:  fitting control data, contains topology and coordinate file names, as
//        well as keys for correlating parameters across systems.  Also
//        contains the frcmod output file name.          
//   tj:  the trajectory control data, contains the statistics output file name
//-----------------------------------------------------------------------------
void FitParams(prmset *mp, trajcon *tj)
{
  int i, j, naltsoln;
  int* InstTally;
  double* b;
  double* bcpy;
  dmat A, Acpy;
  sresamp ms;

  // Sort the command input information
  mp->nbadj = PruneDuplicates(mp->badj, mp->nbadj, mp->nbadj, 2);
  mp->naadj = PruneDuplicates(mp->aadj, mp->naadj, mp->naadj, 3);
  mp->nhadj = PruneDuplicates(mp->hadj, mp->nhadj, mp->nhadj, 4);
  mp->nmadj = PruneDuplicates(mp->madj, mp->nmadj, mp->nmadj, 7);
  
  // Read the conformations and topologies
  if (mp->verbose == 1) {
    printf("\nmdgx >> Unique adjustable parameters enumerated.\n");
  }
  mp->nunisys = 0;
  mp->tpencyc = (prmtop*)malloc(32*sizeof(prmtop));
  j = 32;
  mp->ipolq = DetectIPolQData(mp);
  for (i = 0; i < mp->nconf; i++) {
    ReadSystemConf(&i, mp, &j, tj);
    if (mp->verbose == 1 && (mp->nconf < 100 || i % 100 == 0)) {
      fprintf(stderr, "\rmdgx >> Read conformation %6d of %6d.", i, mp->nconf);
      fflush(stderr);
    }
  }
  ConfUnitConversion(mp);
  FindInstanceRanges(mp);
  if (mp->verbose == 1) {
    printf("\rmdgx >> %6d conformations read.             \n", mp->nconf);
  }
  
  // Read the underlying parameters
  ReadParmFile(mp, tj);
  for (i = 0;i < tj->fmodfile.row; i++) {
    ReadFrcmodFile(mp, tj, i);
  }
  ReadNMROperationsFile(mp);
  if (mp->verbose == 1) {
    printf("mdgx >> Parameter keys read.\n");
  }

  // Insert new atom types by cloning
  // or rewriting existing terms     
  RecastAtomTypes(mp);
  CleaveAtomTypes(mp);

  // Cloning of bond, angle, or dihedral terms may have resulted in
  // duplication.  Remove the duplicates.  CMAPs would not be affected.
  mp->nbond = PruneDuplicates(mp->bonds, mp->nbond, mp->nbond, 2);
  mp->nangl = PruneDuplicates(mp->angls, mp->nangl, mp->nangl, 3);
  mp->ndihe = PruneDuplicates(mp->torsions, mp->ndihe, mp->ntor, 4);
  mp->nimpr = PruneDuplicates(&mp->torsions[mp->ndihe], mp->nimpr, mp->nimpr,
                              5);
  mp->ntor = mp->ndihe + mp->nimpr;
  mp->ncmaps = PruneDuplicates(mp->cmaps, mp->ncmaps, mp->ncmaps, 7);

  // Correlate bond, angle, torsion, and CMAP terms to each topology
  AllocateTermKeyRepositories(mp);  
  MakeTermKey(mp, 2);
  MakeTermKey(mp, 3);
  MakeTermKey(mp, 4);
  MakeTermKey(mp, 7);

  // If NMR restraints are specified as part of the energy surface,
  // they will come in already pruned of most duplicates.  The only
  // question is whether any of them are applicable to the systems
  // at hand, and whether some of the remaining restraints actually
  // apply strictly to the same atoms every time.
  MakeTermKey(mp, 6);

  // Identify terms that are adjustable by
  // linear least-squares fitting, and tie
  // them to fitting matrix columns       
  FindAdjustableTerms(mp);
  CheckGeometryRestraints(mp);
  if (mp->verbose == 1) {
    printf("mdgx >> %d adjustable terms mapped.\n", mp->nparm);
  }

  // Compute molecular mechanics energies and   
  // contributions from various adjustable terms
  for (i = 0; i < mp->nconf; i++) {

    // Bonded interactions
    CompBondedMME(&mp->conf[i], mp, tj);

    // Nonbonded interactions
    CompNBMatrix(mp, &mp->conf[i]);

    // Sum the original molecular mechanics energy
    mp->conf[i].eorig = (mp->zeroNonfit==0) ? SumEnergyMM(&mp->conf[i]) : 0.0;
    if (mp->verbose == 1 && (mp->nconf < 100 || i % 100 == 0)) {
      fprintf(stderr, "\rmdgx >> Computed energy for conformation %6d "
              "of %6d.", i, mp->nconf);
      fflush(stderr);
    }
  }
  if (mp->verbose == 1) {
    printf("\rmdgx >> %6d conformations' energies computed.             "
           "   \n", mp->nconf);
  }

  // Prune bad conformations
  for (i = 0; i < mp->nconf; i++) {
    if (mp->conf[i].PassedEtol == 0) {
      FreeConformation(&mp->conf[i]);
      for (j = i+1; j < mp->nconf; j++) {
        mp->conf[j-1] = mp->conf[j];
      }
      mp->nconf -= 1;
      i--;
    }
  }
  if (mp->verbose == 1) {
    printf("mdgx >> High-energy conformations pruned.\n");
  }

  // Adjust the target energies to meet the averages
  // of the original MM energies for the surviving  
  // conformations of each system                   
  CompEnorm(mp);
  if (mp->verbose == 1) {
    printf("mdgx >> Energies normalized, conformations weighted.\n"
           "mdgx >> Accumulating and checking fitting matrix.\n");
  }
  
  // Create the fitting matrix and solve
  InstTally = GetInstanceCounts(mp);
  mp->nrestraint = CountRestraints(mp);
  A = CreateDmat(mp->nconf + mp->nrestraint, mp->nparm + mp->nunisys, 0);
  b = (double*)malloc((mp->nconf + mp->nrestraint)*sizeof(double));
  FillFittingMatrix(mp, &A, b, InstTally);
  CopyDmat(&Acpy, &A, 0);
  bcpy = CpyDVec(b, A.row);
  WeightFittingMatrix(mp, &A, b);
  WeightFittingMatrix(mp, &Acpy, bcpy);
  
  // Solve the fitting matrix for the global, restrained minimum
  AxbQRRxc(A, b, mp->verbose);
  BackSub(A, b);
  if (mp->verbose == 1) {
    printf("mdgx >> Matrix decomposition complete.\n");
  }

  // Satisfy angle sum contraints
  GeometryConstraints(mp, &Acpy, bcpy, b);

  // Reduce the fitting matrix and make a spectrum of solutions
  ms = SpectrumSolution(&Acpy, bcpy, b, mp);
  if (mp->spectrum == 1) {
    SpectrumEvaluation(mp, &ms);
  }

  // Report the best parameters and the fit
  ParamReport(mp, &ms, tj);
  if (mp->verbose == 1) {
    printf("mdgx >> Parameters written to file %s.\n", tj->dumpname);
    printf("mdgx >> Analysis written to file %s.\n", tj->outbase);
  }

  // Free allocated memory
  free(b);
  free(bcpy);
  free(InstTally);
  DestroyDmat(&A);
  DestroyDmat(&Acpy);
  DestroyPrmset(mp);

  // Exit!  We are done.
  exit(0);
}
