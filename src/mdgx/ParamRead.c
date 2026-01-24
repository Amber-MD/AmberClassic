#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Matrix.h"
#include "mdgxVector.h"
#include "Topology.h"
#include "Parse.h"
#include "Macros.h"
#include "ParamFit.h"

//-----------------------------------------------------------------------------
// CrossRefAtomType: function to find the index of a given atom type among all
//                   of those stored in a parameter set.      
//                                                                      
// Arguments:                                                           
//   mp:             the parameter set                                  
//   aname:          the name of the atom type                          
//-----------------------------------------------------------------------------
int CrossRefAtomType(prmset *mp, char* aname)
{
  int i;
  char typeA[8];

  Type4Char(typeA, aname);
  for (i = 0; i < mp->natom; i++) {
    if (strncmp(mp->atoms[i].atype, typeA, 4) == 0) {
      return i;
    }
  }

  return -1;
}

//-----------------------------------------------------------------------------
// CrossRefBondType: function to find the index of a given atom type among all
//                   of those stored in a parameter set.      
//                                                                      
// Arguments:                                                           
//   mp:             the parameter set                                  
//   [a,b]name:      the names of the atom types                         
//-----------------------------------------------------------------------------
static int CrossRefBondType(prmset *mp, char* aname, char* bname)
{
  int i;
  char typeA[8], typeB[8];

  Type4Char(typeA, aname);
  Type4Char(typeB, bname);
  for (i = 0; i < mp->nbond; i++) {
    if ((strncmp(mp->bonds[i].atype, typeA, 4) == 0 &&
         strncmp(mp->bonds[i].btype, typeB, 4) == 0) ||
        (strncmp(mp->bonds[i].atype, typeB, 4) == 0 &&
         strncmp(mp->bonds[i].btype, typeA, 4) == 0)) {
      return i;
    }
  }

  return -1;
}

//-----------------------------------------------------------------------------
// CrossRefHBondType: function to find the index of a given hydrogen bond
//                    potential among all of those store in a parameter set.
//
// Arguments:
//   mp:             the parameter set
//   [a,b]name:      the names of the atom types
//-----------------------------------------------------------------------------
static int CrossRefHBondType(prmset *mp, char* aname, char* bname)
{
  int i;
  char typeA[8], typeB[8];

  Type4Char(typeA, aname);
  Type4Char(typeB, bname);
  for (i = 0; i < mp->nhb1012; i++) {
    if ((strncmp(mp->hb1012[i].atype, typeA, 4) == 0 &&
         strncmp(mp->hb1012[i].btype, typeB, 4) == 0) ||
        (strncmp(mp->hb1012[i].atype, typeB, 4) == 0 &&
         strncmp(mp->hb1012[i].btype, typeA, 4) == 0)) {
      return i;
    }
  }

  return -1;
}

//-----------------------------------------------------------------------------
// CrossRefAnglType: function to find the index of a given atom type among all
//                   of those stored in a parameter set.      
//                                                                      
// Arguments:                                                           
//   mp:             the parameter set                                  
//   [a,b,c]name:    the names of the atom types                          
//-----------------------------------------------------------------------------
static int CrossRefAnglType(prmset *mp, char* aname, char* bname, char* cname)
{
  int i;
  char typeA[8], typeB[8], typeC[8];

  Type4Char(typeA, aname);
  Type4Char(typeB, bname);
  Type4Char(typeC, cname);
  for (i = 0; i < mp->nangl; i++) {
    if ((strncmp(mp->angls[i].atype, typeA, 4) == 0 &&
         strncmp(mp->angls[i].btype, typeB, 4) == 0 && 
         strncmp(mp->angls[i].ctype, typeC, 4) == 0) ||
        (strncmp(mp->angls[i].atype, typeC, 4) == 0 &&
         strncmp(mp->angls[i].btype, typeB, 4) == 0 &&
         strncmp(mp->angls[i].ctype, typeA, 4) == 0)) {
      return i;
    }
  }

  return -1;
}

//-----------------------------------------------------------------------------
// CrossRefDiheType: function to find the index of a given atom type among all
//                   of those stored in a parameter set.      
//
// Arguments:                                                           
//   mp:             the parameter set                                  
//   [a,b,c,d]name:  the names of the atom types                          
//-----------------------------------------------------------------------------
static int CrossRefDiheType(prmset *mp, cmat lwords, int order)
{
  int i;
  double K, phase, pn;
  char typeA[8], typeB[8], typeC[8], typeD[8];

  Type4Char(typeA, lwords.map[0]);
  Type4Char(typeB, lwords.map[1]);
  Type4Char(typeC, lwords.map[2]);
  Type4Char(typeD, lwords.map[3]);
  for (i = 0; i < mp->ntor; i++) {

    // Skip entries that don't match the type of
    // dihedral we seek (proper or improper).   
    if ((order == 4 && mp->torsions[i].impr == 1) ||
        (order == 5 && mp->torsions[i].impr == 0)) {
      continue;
    }
    if ((strncmp(mp->torsions[i].atype, typeA, 4) == 0 && 
         strncmp(mp->torsions[i].btype, typeB, 4) == 0 &&
         strncmp(mp->torsions[i].ctype, typeC, 4) == 0 &&
         strncmp(mp->torsions[i].dtype, typeD, 4) == 0) ||
        (strncmp(mp->torsions[i].atype, typeD, 4) == 0 &&
         strncmp(mp->torsions[i].btype, typeC, 4) == 0 &&
         strncmp(mp->torsions[i].ctype, typeB, 4) == 0 &&
         strncmp(mp->torsions[i].dtype, typeA, 4) == 0)) {

      // The types match, but what about other    
      // aspects of this fourier term declaration?
      if (order == 4) {
        pn = fabs(atof(lwords.map[7]));
      }
      else if (order == 5) {
        pn = fabs(atof(lwords.map[6]));
      }
      if (fabs(mp->torsions[i].pn - pn) < 1.0e-4) {
        return i;
      }
    }
  }

  return -1;
}

//-----------------------------------------------------------------------------
// ParmFileComment: record a comment from a parameter file.             
//                                                                      
// Arguments:                                                           
//   lwords:    the list of words which may form a comment              
//   wstart:    the starting word                                       
//-----------------------------------------------------------------------------
static char* ParmFileComment(cmat *lwords, int wstart)
{
  int i, j;
  char* commtext;

  commtext = (char*)malloc(MAXLINE*sizeof(char));
  commtext[0] = '\0';
  if (wstart >= lwords->row) {
    return commtext;
  }
  j = 0;
  for (i = wstart; i < lwords->row; i++) {
    sprintf(&commtext[j], "%s ", lwords->map[i]);
    j = strlen(commtext);
  }
  commtext[j-1] = '\0';

  return commtext;
}

//-----------------------------------------------------------------------------
// RunOfWords: check to see that there is a run of atom types or numbers in
//             selected positions of a word list.                    
//                                                                      
// Arguments:                                                           
//   lwords:    the list of words                                       
//   istart:    the start position of the atom types run                
//   iend:      the final position of the atom types run                
//   style:     set to 0 for atom types, 1 for numbers                  
//-----------------------------------------------------------------------------
static int RunOfWords(cmat *lwords, int istart, int iend, int style)
{
  int i;

  for (i = istart; i <= iend; i++) {
    if (i >= 0 && i < lwords->row) {
      if (style == 0 && WordIsAtomType(lwords->map[i]) == 0) {
        return 0;
      }
      if (style == 1 && WordIsNumber(lwords->map[i]) == 0) {
        return 0;
      }
    }
  }

  return 1;
}

//-----------------------------------------------------------------------------
// FormatTest: function to test the format of a line for matches to each type
//             of parameter file data.                             
//                                                                      
// Arguments:                                                           
//   lwords:    the line, converted to words                            
//   order:     integer denoting the type of input needed               
//-----------------------------------------------------------------------------
static int FormatTest(cmat *lwords, int order)
{
  // Atoms
  if (order == 1 || order == 8) {
    if (lwords->row < 2) {
      return 0;
    }
    if (WordIsAtomType(lwords->map[0]) == 1 &&
        WordIsNumber(lwords->map[1]) == 1) {
      if (order == 8 && WordIsNumber(lwords->map[2]) == 0) {
        return 0;
      }
      return 1;
    }
  }

  // Bonds, hydrogen bonds
  if (order == 2 || order == 9) {
    if (lwords->row < 4) {
      return 0;
    }
    if (RunOfWords(lwords, 0, 1, 0) == 1 && RunOfWords(lwords, 2, 3, 1) == 1) {
      return 1;
    }
  }

  // Angles
  if (order == 3) {
    if (lwords->row < 5) {
      return 0;
    }
    if (RunOfWords(lwords, 0, 2, 0) == 1 && RunOfWords(lwords, 3, 4, 1) == 1) {
      return 1;
    }
  }

  // Dihedrals
  if (order == 4) {
    if (lwords->row < 8) {
      return 0;
    }
    if (RunOfWords(lwords, 0, 3, 0) == 1 && RunOfWords(lwords, 4, 7, 1) == 1) {
      return 1;
    }
  }

  // Improper Dihedrals
  if (order == 5) {
    if (lwords->row < 7) {
      return 0;
    }
    if (RunOfWords(lwords, 0, 3, 0) == 1 && RunOfWords(lwords, 4, 6, 1) == 1) {
      return 1;
    }
  }

  // Line of atom types
  if (order == 6 || order == 7) {
    if (lwords->row == 0) {
      return 0;
    }
    if (RunOfWords(lwords, 0, lwords->row-1, 0) == 1) {
      return 1;
    }
  }
  
  // Nope, the line didn't conform
  return 0;
}

//-----------------------------------------------------------------------------
// ScanDeclarations: scan through the parameter or frcmod file for lines that
//                   conform to the required format.         
//                                                                      
// Arguments:                                                           
//   mp:      the parameter set                                         
//   cfi:     the parameter input file, converted to a character matrix 
//   order:   the order of the comparison                               
//   nentry:  the number of entries already logged, to which this scan will add
//   strline: the line on which to start the scan                       
//-----------------------------------------------------------------------------
static int ScanDeclarations(prmset *mp, cmat *cif, int order, int nentry,
                            int strline, char* tercode)
{
  int i, j, k, m, llen, ndash, represent, entryID, newadd;
  int segter, formatmatch, matchfound, maxentry;
  char line[MAXLINE], typetmp[8];
  cmat lwords, workingtypes;

  // Allocate memory
  maxentry = nentry + 32;
  if (nentry == 0) {
    if (order == 1) {
      mp->atoms = (xatomdef*)malloc(maxentry*sizeof(xatomdef));
    }
    else if (order == 2) {
      mp->bonds = (xbonddef*)malloc(maxentry*sizeof(xbonddef));
    }
    else if (order == 3) {
      mp->angls = (xangldef*)malloc(maxentry*sizeof(xangldef));
    }
    else if (order == 4) {
      mp->torsions = (torterm*)malloc(maxentry*sizeof(torterm));
    }
    else if (order == 6) {
      mp->Hydrophilics = CreateCmat(1, 5);
      mp->Hydrophilics.row = 0;
    }
    else if (order == 7) {
      mp->eqgroups = (eagrp*)malloc(maxentry*sizeof(eagrp));
    }
    else if (order == 9) {
      mp->hb1012 = (xhb1012def*)malloc(maxentry*sizeof(xhb1012def));
    }
  }
  else {
    if (order == 1) {
      mp->atoms = (xatomdef*)realloc(mp->atoms, maxentry*sizeof(xatomdef));
    }
    else if (order == 2) {
      mp->bonds = (xbonddef*)realloc(mp->bonds, maxentry*sizeof(xbonddef));
    }
    else if (order == 3) {
      mp->angls = (xangldef*)realloc(mp->angls, maxentry*sizeof(xangldef));
    }
    else if (order == 4 || order == 5) {
      mp->torsions = (torterm*)realloc(mp->torsions, maxentry*sizeof(torterm));
    }
    else if (order == 9) {
      mp->hb1012 = (xhb1012def*)realloc(mp->hb1012,
                                        maxentry*sizeof(xhb1012def));
    }
  }
  if (order == 4) {
    workingtypes = CreateCmat(4, 8);
  }

  // Loop through the available portions of the input file
  matchfound = 0;
  segter = 0;
  for (i = strline; i < cif->row; i++) {

    // Copy the line so it can be modified
    strcpy(line, cif->map[i]);
    RemoveWhiteSpace(line, 0, ' ');
    llen = strlen(line);

    // Remove dashes
    if (order >= 2 && order <= 5) {
      ndash = 0;
      for (j = 0; j < llen; j++) {
        if (line[j] == '-') {
          ndash++;
          line[j] = ' ';
        }

        // Break if the maximum number of
        // dashes have been encountered  
        if (ndash == order-1) {
          break;
        }
      }
    }

    // Scan the line into words
    lwords = ParseWords(line);

    // Determine whether the format matches
    formatmatch = FormatTest(&lwords, order);
    if (formatmatch == 0 && order == 4) {

      // This would fail the format test, but there is one more chance. 
      // If the words on this line are numbers such that four atom types
      // could be inserted to make it work out, fix the line and roll   
      // with that.                                                     
      if (RunOfWords(&lwords, 0, 3, 1) == 1 &&
          workingtypes.map[0][0] != '\0') {
        j = lwords.row;
        lwords = ReallocCmat(&lwords, j+4, lwords.col);
        for (k = lwords.row-1; k >= 4; k--) {
          strcpy(lwords.map[k], lwords.map[k-4]);
        }
        for (k = 0; k < 4; k++) {
          strncpy(lwords.map[k], workingtypes.map[k], 2);
          lwords.map[k][2] = '\0';
        }
      }
      formatmatch = FormatTest(&lwords, order);
    }

    // If the format does match, proceed accordingly
    if (formatmatch == 1) {
      matchfound = 1;             
      if (order == 1 || order == 8) {
        entryID = CrossRefAtomType(mp, lwords.map[0]);
      }
      else if (order == 2) {
        entryID = CrossRefBondType(mp, lwords.map[0], lwords.map[1]);
      }
      else if (order == 3) {
        entryID = CrossRefAnglType(mp, lwords.map[0], lwords.map[1],
                                   lwords.map[2]);
      }
      else if (order == 4 || order == 5) {
        entryID = CrossRefDiheType(mp, lwords, order);
      }
      else if (order == 9) {
        entryID = CrossRefHBondType(mp, lwords.map[0], lwords.map[1]);
      }
      else {
        entryID = -1;
      }
      if (entryID == -1) {
        newadd = 1;
        entryID = nentry;
      }
      else {
        newadd = 0;
      }
      if (order == 1) {
        Type4Char(mp->atoms[entryID].atype, lwords.map[0]);
        mp->atoms[entryID].mass = atof(lwords.map[1]);
        if (lwords.row > 2) {
          if (WordIsNumber(lwords.map[2]) == 1) {
            mp->atoms[entryID].apol = atof(lwords.map[2]);
            mp->atoms[entryID].comment = ParmFileComment(&lwords, 3);
          }
          else {
            mp->atoms[entryID].apol = 0.0;
            mp->atoms[entryID].comment = ParmFileComment(&lwords, 2);
          }
        }
        else {
          mp->atoms[entryID].apol = 0.0;
          mp->atoms[entryID].comment = ParmFileComment(&lwords, 2);
        }
        mp->atoms[entryID].ljsig = 0.0;
        mp->atoms[entryID].ljeps = 0.0;
      }
      else if (order == 2) {
        Type4Char(mp->bonds[entryID].atype, lwords.map[0]);
        Type4Char(mp->bonds[entryID].btype, lwords.map[1]);
        mp->bonds[entryID].K = atof(lwords.map[2]);
        mp->bonds[entryID].l0 = atof(lwords.map[3]);
        if (lwords.row >= 8 && WordIsNumber(lwords.map[4]) == 1 &&
            WordIsNumber(lwords.map[5]) == 1 &&
            WordIsNumber(lwords.map[6]) == 1 &&
            WordIsNumber(lwords.map[7]) == 1) {
          mp->bonds[entryID].Kpull   = atof(lwords.map[4]);
          mp->bonds[entryID].lpull0  = atof(lwords.map[5]);
          mp->bonds[entryID].Kpress  = atof(lwords.map[6]);
          mp->bonds[entryID].lpress0 = atof(lwords.map[7]);
          mp->bonds[entryID].isAug   = 1;
          mp->bonds[entryID].comment = ParmFileComment(&lwords, 8);
        }
        else {
          mp->bonds[entryID].Kpull   = 0.0;
          mp->bonds[entryID].lpull0  = 100.0;
          mp->bonds[entryID].Kpress  = 0.0;
          mp->bonds[entryID].lpress0 = 100.0;
          mp->bonds[entryID].isAug   = 0;
          mp->bonds[entryID].comment = ParmFileComment(&lwords, 4);
        }
      }
      else if (order == 3) {
        Type4Char(mp->angls[entryID].atype, lwords.map[0]);
        Type4Char(mp->angls[entryID].btype, lwords.map[1]);
        Type4Char(mp->angls[entryID].ctype, lwords.map[2]);
        mp->angls[entryID].K = atof(lwords.map[3]);
        mp->angls[entryID].th0 = atof(lwords.map[4])*PI/180.0;
        mp->angls[entryID].comment = ParmFileComment(&lwords, 5);
      }
      else if (order == 4 || order == 5) {
        Type4Char(mp->torsions[entryID].atype, lwords.map[0]);
        Type4Char(mp->torsions[entryID].btype, lwords.map[1]);
        Type4Char(mp->torsions[entryID].ctype, lwords.map[2]);
        Type4Char(mp->torsions[entryID].dtype, lwords.map[3]);
        if (order == 4) {
          for (j = 0; j < 4; j++) {
            strncpy(workingtypes.map[j], lwords.map[j], 2);
            workingtypes.map[j][2] = '\0';
          }
          mp->torsions[entryID].K = atof(lwords.map[5]) / atof(lwords.map[4]);
          mp->torsions[entryID].phase = atof(lwords.map[6])*PI/180.0;
          mp->torsions[entryID].pn = atof(lwords.map[7]);
          mp->torsions[entryID].impr = 0;
          mp->torsions[entryID].comment = ParmFileComment(&lwords, 8);
        }
        else if (order == 5) {
          mp->torsions[entryID].K = atof(lwords.map[4]);
          mp->torsions[entryID].phase = atof(lwords.map[5])*PI/180.0;
          mp->torsions[entryID].pn = atof(lwords.map[6]);
          mp->torsions[entryID].impr = 1;
          mp->torsions[entryID].comment = ParmFileComment(&lwords, 7);
        }
        if (mp->torsions[entryID].pn < 0.0) {
          mp->torsions[entryID].singlet = -1;
          mp->torsions[entryID].pn = fabs(mp->torsions[entryID].pn);
        }
        else {
          mp->torsions[entryID].singlet = 1;
        }

      }
      else if (order == 6) {
        k = mp->Hydrophilics.row;
        mp->Hydrophilics = ReallocCmat(&mp->Hydrophilics, k + lwords.row, 5);
        for (j = 0; j < lwords.row; j++) {
          Type4Char(mp->Hydrophilics.map[k+j], lwords.map[j]);
        }
      }
      else if (order == 7) {
        mp->eqgroups = (eagrp*)realloc(mp->eqgroups,
                                       (entryID+1)*sizeof(eagrp));
        mp->eqgroups[entryID].natom = lwords.row;
        mp->eqgroups[entryID].types =
          (char*)malloc((4*lwords.row+1)*sizeof(char));
        for (j = 0; j < lwords.row; j++) {
          Type4Char(&mp->eqgroups[entryID].types[4*j], lwords.map[j]);
        }
      }
      else if (order == 8) {

        // In this case, the objective is not to add
        // to a growing list, but to contribute more
        // information to what already exists.
        mp->atoms[entryID].ljsig = atof(lwords.map[1]);
        mp->atoms[entryID].ljeps = atof(lwords.map[2]);
        Type4Char(typetmp, lwords.map[0]);
        for (j = 0; j < mp->neqgroups; j++) {
          represent = 0;
          for (k = 0; k < mp->eqgroups[j].natom; k++) {
            if (strncmp(&mp->eqgroups[j].types[4*k], typetmp, 4) == 0) {
              represent = 1;
            }
          }
          if (represent == 1) {
            for (k = 0; k < mp->eqgroups[j].natom; k++) {
              m = CrossRefAtomType(mp, &mp->eqgroups[j].types[4*k]);
              if (m == -1) {
                printf("ScanDeclarations >> Error.  Cross-reference to atom "
                       "type %.4s failed.\n", &mp->eqgroups[j].types[4*k]);
                exit(1);
              }
              mp->atoms[m].ljsig = atof(lwords.map[1]);
              mp->atoms[m].ljeps = atof(lwords.map[2]);
            }
          }
        }
      }
      else if (order == 9) {
        Type4Char(mp->hb1012[entryID].atype, lwords.map[0]);
        Type4Char(mp->hb1012[entryID].btype, lwords.map[1]);
        mp->hb1012[entryID].Aterm = atof(lwords.map[2]);
        mp->hb1012[entryID].Bterm = atof(lwords.map[3]);
        mp->hb1012[entryID].comment = ParmFileComment(&lwords, 4);
      }

      // Increment the number of entries
      if (newadd == 1) {
        nentry += 1;
        if (order == 1) {
          mp->natom = nentry;
        }
        else if (order == 2) {
          mp->nbond = nentry;
        }
        else if (order == 3) {
          mp->nangl = nentry;
        }
        else if (order == 4 || order == 5) {
          mp->ntor = nentry;
        }
        else if (order == 7) {
          mp->neqgroups = nentry;
        }
        else if (order == 9) {
          mp->nhb1012 = nentry;
        }
      }
      if (nentry >= maxentry) {
        maxentry += 32;
        if (order == 1) {
          mp->atoms = (xatomdef*)realloc(mp->atoms, maxentry*sizeof(xatomdef));
        }
        else if (order == 2) {
          mp->bonds = (xbonddef*)realloc(mp->bonds, maxentry*sizeof(xbonddef));
        }
        else if (order == 3) {
          mp->angls = (xangldef*)realloc(mp->angls, maxentry*sizeof(xangldef));
        }
        else if (order == 4 || order == 5) {
          mp->torsions = (torterm*)realloc(mp->torsions,
                                           maxentry*sizeof(torterm));
        }
      }
    }

    // Consider the possibility that this is a segment terminator,
    // or that this is the end of all parameter declarations.     
    if ((strcmp(tercode, "BLANKLINE") == 0 && line[0] == '\n') ||
        strcmp(tercode, "ONELINE") == 0) {
      segter = 1;
      i++;
    }
    if (strncmp(line, "END", 3) == 0 || strncmp(line, "end", 3) == 0) {
      segter = 1;
      i = cif->row;
    }

    // Free allocated memory
    DestroyCmat(&lwords);

    // Terminate the search
    if (segter == 1) {
      break;
    }
  }

  // Free allocated memory
  if (order == 4) {
    DestroyCmat(&workingtypes);
  }

  // What about the next line of the file?
  if (matchfound == 0) {
    return strline;
  }
  else {
    return i;
  }
}

//-----------------------------------------------------------------------------
// ReadFrcmodCmapData: encapsulates input for CMAP parameters in a frcmod
//                     file, which have a different flavor than other frcmod
//                     components.
//
// Arguments:
//   mp:       the growing parameter set
//   cfi:      the parameter input file, converted to a character matrix 
//   strline:  the line on which to start reading CMAP data
//-----------------------------------------------------------------------------
static void ReadFrcmodCmapData(prmset *mp, cmat *cfi, int strline)
{
  int i, j, k, countseen, mapcon, rc, cc, prmcon, cmres, ndataln;
  int initCmapCount;
  double* b;
  double* tvec;
  double* spln;
  imat termlimits;
  dmat A;
  cmat lwords, dtwords;
  
  // Search for CMAP_COUNT flags within the remainder of the frcmod file.
  // These will give an indication of the range over which each CMAP
  // term's data can be found.
  initCmapCount = mp->ncmaps;
  if (mp->ncmaps == 0) {
    mp->cmaps = (xcmapdef*)malloc(sizeof(xcmapdef));
  }
  termlimits = CreateImat(2, 64);
  termlimits.map[0][0] = strline;
  for (i = strline; i < cfi->row; i++) {

    // Scan the line into words
    lwords = ParseWords(cfi->map[i]);
    if (lwords.row == 3) {
      if (strcmp(lwords.map[0], "%FLAG") == 0 &&
	  strcmp(lwords.map[1], "CMAP_COUNT") == 0) {
	termlimits.map[1][mp->ncmaps] = i;
	mp->ncmaps += 1;
        countseen = strtol(lwords.map[2], NULL, 10);
	if (countseen != mp->ncmaps - initCmapCount) {
	  if (initCmapCount == 0) {
            printf("ReadFrcmodCmapData >> Warning.  CMAP count may be out of "
                   "order.\nReadFrcmodCmapData >> Read map %3d but expected "
		   "%3d.\n", countseen, mp->ncmaps);
	  }
	  else {
            printf("ReadFrcmodCmapData >> Warning.  CMAP count may be out of "
                   "order.  Read map\nReadFrcmodCmapData >> %3d but expected "
                   "%3d (after reading %3d maps previously).\n", countseen,
		   mp->ncmaps, initCmapCount);
	  }
	}
	if (mp->ncmaps == termlimits.col) {
	  termlimits = ReallocImat(&termlimits, 2, termlimits.col + 32);
	}
	termlimits.map[0][mp->ncmaps] = i;
      }
    }
    DestroyCmat(&lwords);
  }
  if (mp->verbose == 1) {
    if (initCmapCount == 0) {
      printf("mdgx >> Detected %d CMAP surfaces.\n", mp->ncmaps);
    }
    else {
      printf("mdgx >> Detected %d additional CMAP surfaces.\n",
	     mp->ncmaps - initCmapCount);
    }
  }
  
  // Pre-allocate space for the CMAPs, then look for their components
  mp->cmaps = (xcmapdef*)realloc(mp->cmaps, mp->ncmaps * sizeof(xcmapdef));
  for (i = initCmapCount; i < mp->ncmaps; i++) {
    mp->cmaps[i].title = (char*)malloc(MAXLINE * sizeof(char));
  }

  // Detect the resolution of each map (exactly one CMAP_RESOLUTION
  // keyword must be present for every map).  Fill in default values.
  mapcon = initCmapCount;
  for (i = strline; i < cfi->row; i++) {

    // Scan the line into words
    lwords = ParseWords(cfi->map[i]);
    if (lwords.row == 3) {
      if (strcmp(lwords.map[0], "%FLAG") == 0 &&
          strcmp(lwords.map[1], "CMAP_RESOLUTION") == 0) {
        cmres = strtol(lwords.map[2], NULL, 10);
        mp->cmaps[mapcon].esrf = CreateDmat(cmres, cmres, 0);
        mp->cmaps[mapcon].dphi = CreateDmat(cmres, cmres, 0);
        mp->cmaps[mapcon].dpsi = CreateDmat(cmres, cmres, 0);
        mp->cmaps[mapcon].d2pp = CreateDmat(cmres, cmres, 0);
	if (mp->cmapFitResln > 0) {
	  j = mp->cmapFitResln;
          mp->cmaps[mapcon].esrfFit = CreateDmat(j, j, 0);
          mp->cmaps[mapcon].dphiFit = CreateDmat(j, j, 0);
          mp->cmaps[mapcon].dpsiFit = CreateDmat(j, j, 0);
          mp->cmaps[mapcon].d2ppFit = CreateDmat(j, j, 0);
	}
	else {
          mp->cmaps[mapcon].esrfFit = CreateDmat(cmres, cmres, 0);
          mp->cmaps[mapcon].dphiFit = CreateDmat(cmres, cmres, 0);
          mp->cmaps[mapcon].dpsiFit = CreateDmat(cmres, cmres, 0);
          mp->cmaps[mapcon].d2ppFit = CreateDmat(cmres, cmres, 0);
	}
        sprintf(mp->cmaps[mapcon].atoms, "C   N   CA  C   N   ");
        mp->cmaps[mapcon].relResID[0] = -1;
        mp->cmaps[mapcon].relResID[1] =  0;
        mp->cmaps[mapcon].relResID[2] =  0;
        mp->cmaps[mapcon].relResID[3] =  0;
        mp->cmaps[mapcon].relResID[4] =  1;
	mp->cmaps[mapcon].rstw = -1.0;
	mapcon++;
      }
    }
  }

  // Read the parameters and other information
  mapcon = initCmapCount - 1;
  for (i = strline; i < cfi->row; i++) {

    // Scan the line into words
    lwords = ParseWords(cfi->map[i]);
    if (lwords.row == 3) {
      if (strcmp(lwords.map[0], "%FLAG") == 0 &&
	  strcmp(lwords.map[1], "CMAP_COUNT") == 0) {

	// Increment the counter to the map number that will now be filled.
	// (The CMAP_COUNT keyword must appear first in the frcmod file.)
	mapcon++;
      }
      else if (strcmp(lwords.map[0], "%FLAG") == 0) {
	if (strcmp(lwords.map[1], "CMAP_RESLIST") == 0) {
	  mp->cmaps[mapcon].nres = strtol(lwords.map[2], NULL, 10);
	  dtwords = ParseWords(cfi->map[i+1]);
	  if (dtwords.row != mp->cmaps[mapcon].nres) {
	    printf("ReadFrcmodCmapData >> Warning.  Expected %d residues for "
		   "CMAP\nReadFrcmodCmapData >> %s.\nReadFrcmodCmapData >> "
		   "Found %d candidate names on the next line.\n",
		   mp->cmaps[mapcon].nres, mp->cmaps[mapcon].title,
		   dtwords.row);
	    if (dtwords.row < mp->cmaps[mapcon].nres) {
	      printf("ReadFrcmodCmapData >> This map will apply to only %d "
		     "residues.\n", dtwords.row);
	      mp->cmaps[mapcon].nres = dtwords.row;
	    }
	    else {
	      printf("ReadFrcmodCmapData >> This map will apply to the first "
		     "%d named residues.\n", mp->cmaps[mapcon].nres);
	    }
	  }
	  mp->cmaps[mapcon].resnames = CreateCmat(mp->cmaps[mapcon].nres, 5);
	  for (j = 0; j < mp->cmaps[mapcon].nres; j++) {
	    TouchupAtomName(dtwords.map[j]);
	    for (k = 0; k < 4; k++) {
	      mp->cmaps[mapcon].resnames.map[j][k] = dtwords.map[j][k];
	    }
	  }
	}
      }
    }
    else if (lwords.row == 2) {
      if (strcmp(lwords.map[0], "%FLAG") == 0) {
	if (strcmp(lwords.map[1], "CMAP_TITLE") == 0 && i < cfi->row - 1) {
	  j = 0;
	  while (cfi->map[i+1][j] != '\n') {
	    mp->cmaps[mapcon].title[j] = cfi->map[i+1][j];
	    j++;
	  }
	  mp->cmaps[mapcon].title[j] = '\0';
	}
	else if (strcmp(lwords.map[1], "CMAP_PARAMETER") == 0) {
	  cmres = mp->cmaps[mapcon].esrf.row;
	  ndataln = (cmres*cmres + 7) / 8;
	  rc = 0;
	  cc = 0;
	  for (j = i+1; j < i + 1 + ndataln; j++) {
            dtwords = ParseWords(cfi->map[j]);
	    for (k = 0; k < dtwords.row; k++) {
	      if (WordIsNumber(dtwords.map[k]) == 1 && rc < cmres) {
		mp->cmaps[mapcon].esrf.map[rc][cc] = strtod(dtwords.map[k],
                                                            NULL);
		cc++;
		if (cc == mp->cmaps[mapcon].esrf.col) {
		  cc = 0;
		  rc++;
		}
	      }
	      else if (rc < cmres) {
		printf("ReadFrcmodCmapData >> Error.  Expected CMAP value "
                       "(%2d,%2d) but found\nReadFrcmodCmapData >>\"%s\".\n",
		       rc, cc, dtwords.map[k]);
	      }
	      else {
		printf("ReadFrcmodCmapData >> Error.  Expected end of CMAP "
                       "values, but found\nReadFrcmodCmapData >>\"%s\".\n",
		       dtwords.map[k]);
	      }
	    }
	    DestroyCmat(&dtwords);
	  }
	}
      }
    }
    else if (lwords.row == 7) {
      if (strcmp(lwords.map[0], "%FLAG") == 0) {
	if (strcmp(lwords.map[1], "CMAP_ATMLIST") == 0) {
	  for (j = 2; j < 7; j++) {
	    if (strlen(lwords.map[j]) > 4) {
	      printf("ReadFrcmodCmapData >> Error.  Atom name %s has too "
                     "many characters.\n", lwords.map[j]);
	      exit(1);
	    }
	    for (k = 0; k < strlen(lwords.map[j]); k++) {
	      mp->cmaps[mapcon].atoms[4*(j-2) + k] = lwords.map[j][k];
	    }
	    TouchupAtomName(&mp->cmaps[mapcon].atoms[4*(j-2)]);
	  }
	}
	if (strcmp(lwords.map[1], "CMAP_RESIDX") == 0) {
	  for (j = 2; j < 7; j++) {
            mp->cmaps[mapcon].relResID[j-2] =
              strtol(lwords.map[j], NULL, 10);
	  }
	}
      }
    }
    DestroyCmat(&lwords);
  }
  if (mp->verbose == 1) {
    printf("mdgx >> CMAP parameters read.\n");
  }
  
  // Compute derivatives (use the linear method only).  Check that
  // all CMAPs are of consistent size.
  for (i = initCmapCount; i < mp->ncmaps; i++) {
    if (mp->verbose == 1) {
      fprintf(stderr, "\rmdgx >> Computing derivatives for CMAP %4d", i);
      fflush(stderr);
    }
    ComputeLinearCmapDerivatives(&mp->cmaps[i].esrf, &mp->cmaps[i].dphi, 
				 &mp->cmaps[i].dpsi, &mp->cmaps[i].d2pp);
  }
  if (mp->verbose == 1) {
    printf("\n");
  }
  
  // Compute the stencils.  These apply to the fitted surface, which may be
  // of different dimension than the surface read from the frcmod file.
  // The purpose of allowing a different size on the fitted surface is to
  // coarsen the map.  It's overkill to make a 24-point surface with
  // bicubic spline interpolation.  When the derivatives are computed
  // properly (in periodic boundary conditions), an eight-point map should
  // cover the same level of detail as a cosine series with a periodicity
  // as high as 4, and cut the number of fitted surface points by a factor
  // of 9.
  if (mp->verbose == 1) {
    printf("mdgx >> Computing CMAP stencils.\n");
  }
  for (i = initCmapCount; i < mp->ncmaps; i++) {
    mp->cmaps[i].stencil = Compute1DStencil(mp->cmaps[i].esrfFit.row,
                                            mp->cmaps[i].esrf.map[0]);
  }
}

//-----------------------------------------------------------------------------
// ReadParmFile: read a force field parameter file, to help identify wildcard
//               parameters and thus pare down the number of variables to fit.
//
// Arguments:                                                           
//   mp:      the fitting data (contains a list of all systems)         
//   tj:      trajectory control data (contains the frcmod file name)   
//-----------------------------------------------------------------------------
void ReadParmFile(prmset *mp, trajcon *tj)
{
  int nextline, tmpi;
  cmat cfi;

  // Test file existence
  if (tj->parmfile[0] == '\0') {
    printf("ReadParmFile >> Error.  Force field parameter file not specified."
           "\n");
    exit(1);
  }
  cfi = Ascii2Mem(tj->parmfile, 256, 2,
                  "Force field parameter file not found.");

  // Read the atom declarations
  mp->natom = 0;
  nextline = ScanDeclarations(mp, &cfi, 1, mp->natom, 1, "BLANKLINE");

  // Record the hydrophilic atoms
  nextline = ScanDeclarations(mp, &cfi, 6, 0, nextline, "ONELINE");

  // Read bond, angle, dihedral, and improper declarations
  mp->nbond = 0;
  nextline = ScanDeclarations(mp, &cfi, 2, mp->nbond, nextline, "BLANKLINE");
  mp->nangl = 0;
  nextline = ScanDeclarations(mp, &cfi, 3, mp->nangl, nextline, "BLANKLINE");
  mp->ntor = 0;
  nextline = ScanDeclarations(mp, &cfi, 4, mp->ntor, nextline, "BLANKLINE");
  nextline = ScanDeclarations(mp, &cfi, 5, mp->ntor, nextline, "BLANKLINE");
  mp->nhb1012 = 0;
  nextline = ScanDeclarations(mp, &cfi, 9, mp->nhb1012, nextline, "BLANKLINE");

  // Read atom equivalencies
  mp->neqgroups = 0;
  nextline = ScanDeclarations(mp, &cfi, 7, 0, nextline, "BLANKLINE");

  // Read non-bonded parameters and    
  // cross-reference with equivalencies
  nextline = ScanDeclarations(mp, &cfi, 8, 0, nextline, "BLANKLINE");

  // Look for CMAPs--probably won't find any, but
  // the array must be initialzed here.
  mp->ncmaps = 0;
  if (DetectInCmat(&cfi, "CMAP", &nextline, &tmpi, nextline, 0) == 1) {
    ReadFrcmodCmapData(mp, &cfi, nextline);
  }
  
  // Free allocated memory
  DestroyCmat(&cfi);
}

//-----------------------------------------------------------------------------
// ReadFrcmodFile: read a force field modification file, to add to the  
//                 parameters already detected by ReadParmFile() above. 
//
// Arguments:                                                           
//   mp:      the fitting data (contains a list of all systems)         
//   tj:      trajectory control data (contains the frcmod file name)   
//   nfmod:   the number of the frcmod file to read (there could be multiple)
//-----------------------------------------------------------------------------
void ReadFrcmodFile(prmset *mp, trajcon *tj, int nfmod)
{
  int i, nextline;
  cmat cfi;

  // Test file existence
  if (tj->fmodfile.map[nfmod][0] == '\0') {
    return;
  }
  cfi = Ascii2Mem(tj->fmodfile.map[nfmod], 256, 2,
                  "Force field modification file not found.");

  // Read bond, angle, dihedral, and improper declarations
  if (DetectInCmat(&cfi, "MASS", &nextline, &i, 0, 0) == 1) {
    nextline = ScanDeclarations(mp, &cfi, 1, mp->natom, nextline, "BLANKLINE");
  }
  if (DetectInCmat(&cfi, "BOND", &nextline, &i, nextline, 0) == 1) {
    nextline = ScanDeclarations(mp, &cfi, 2, mp->nbond, nextline, "BLANKLINE");
  }
  if (DetectInCmat(&cfi, "ANGL", &nextline, &i, nextline, 0) == 1) {
    nextline = ScanDeclarations(mp, &cfi, 3, mp->nangl, nextline, "BLANKLINE");
  }
  if (DetectInCmat(&cfi, "DIHE", &nextline, &i, nextline, 0) == 1) {
    nextline = ScanDeclarations(mp, &cfi, 4, mp->ntor, nextline, "BLANKLINE");
  }
  if (DetectInCmat(&cfi, "IMPR", &nextline, &i, nextline, 0) == 1) {
    nextline = ScanDeclarations(mp, &cfi, 5, mp->ntor, nextline, "BLANKLINE");
  }

  // Modifications to hydrogen bonding parameters
  if (DetectInCmat(&cfi, "HBON", &nextline, &i, nextline, 0) == 1) {
    nextline = ScanDeclarations(mp, &cfi, 9, mp->nhb1012, nextline,
                                "BLANKLINE");
  }

  // Read modified nobonded parameters
  if (DetectInCmat(&cfi, "NONB", &nextline, &i, nextline, 0) == 1) {
    nextline = ScanDeclarations(mp, &cfi, 8, 0, nextline, "BLANKLINE");
  }
  
  // Read CMAP parameters
  if (DetectInCmat(&cfi, "CMAP", &nextline, &i, nextline, 0) == 1) {
    ReadFrcmodCmapData(mp, &cfi, nextline);
  }
}

//-----------------------------------------------------------------------------
// RecastInComment: this function tries to find instances of the former atoms
//                  in comments, and replace them.                
//
// Arguments:                                                           
//   comment:    the comment to parse                                   
//   recast:     the type replacement information                       
//-----------------------------------------------------------------------------
static void RecastInComment(char* comment, typeswitch recast)
{
  int i, slen;

  slen = strlen(comment);
  for (i = 0; i < slen-1; i++) {
    if (comment[i] == recast.orig[0] && comment[i+1] == recast.orig[1]) {

      // This is a possible hit
      if ((i == 0 || comment[i-1] == ' ' || comment[i-1] == '-' ||
           comment[i-1] == '(') &&
          (i == slen-2 || comment[i+2] == ' ' || comment[i+2] == '-' ||
           comment[i+2] == ')')) {

        // This is a hit; replace the atom type name
        comment[i] = recast.pnew[0];
        comment[i+1] = recast.pnew[1];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// BubbleChar4: bubble sort operation for two arrays of 4 characters.   
//                                                                      
// Arguments:                                                           
//   [a,b]:    the two strings to compare and possibly swap             
//-----------------------------------------------------------------------------
void BubbleChar4(char* a, char* b)
{
  int i;
  char c;

  if ((str4cmp(b, "X   ") == 0 && str4cmp(a, "X   ") != 0) ||
      (str4cmp(a, "X   ") != 0 && strncmp(a, b, 4) > 0)) {
    for (i = 0; i < 4; i++) {
      SWAP(a[i], b[i], c);
    }
  }
}

//-----------------------------------------------------------------------------
// AlphabetizeImproper: this function will re-order the A, B and D atom types
//                      in an improper dihedral to ensure that they appear in
//                      alphabetical order.              
//                                                                      
// Arguments:                                                           
//   ti:      the torsion definition                                    
//-----------------------------------------------------------------------------
static void AlphabetizeImproper(torterm *ti)
{
  // First, check to see that this really is an improper
  if (ti->impr == 0) {
    return;
  }
  BubbleChar4(ti->atype, ti->btype);
  BubbleChar4(ti->btype, ti->dtype);
  BubbleChar4(ti->atype, ti->btype);
}

//-----------------------------------------------------------------------------
// ReplaceAtomType: rewrite an atom type name if it matches.            
//                                                                      
// Arguments:                                                           
//   ctype:   the current atom type                                     
//   orig:    the target atom type to change                            
//   new:     the name of the new atom type                             
//-----------------------------------------------------------------------------
static void ReplaceAtomType(char* ctype, char* orig, char* new)
{
  if (str4cmp(ctype, orig) == 0) {
    strncpy(ctype, new, 4);
  }
}

//-----------------------------------------------------------------------------
// RecordAtomTypeChange: record the changes made to the atom types of specific
//                       atoms in each topology.               
//                                                                      
// Arguments:                                                           
//   mp:     the master parameter set                                   
//   tp:     the particular topology that is being changed              
//   atnum:  the number of the atom in the topology                     
//   origt:  the original atom type                                     
//   newt:   the new atom type                                          
//   style:  the style of change (recast, branch)                       
//   src:    ambmask string dictating the change                        
//-----------------------------------------------------------------------------
static void RecordAtomTypeChange(prmset *mp, prmtop *tp, int atnum,
                                 char* origt, char* newt, char* style,
                                 char* src)
{
  int resid;

  resid = LocateResID(tp, atnum, 0, tp->nres);
  sprintf(mp->ChangeLog.map[mp->nchng],
          " %-16.16s  %.4s %2d %.4s  %.2s -> %.2s  %.6s  %-27.27s\n",
          tp->source, &tp->AtomNames[4*atnum], resid, &tp->ResNames[4*resid],
          origt, newt, style, src);
  mp->nchng += 1;
  if (mp->nchng == mp->ChangeLog.row) {
    mp->ChangeLog = ReallocCmat(&mp->ChangeLog, mp->nchng+32, 80);
  }
}

//-----------------------------------------------------------------------------
// RecastAtomTypes: rename all instances of a particular atom type in the
//                  parameter set.  This function performs a hard renaming of
//                  the atom type wherever it can find it, and then checks
//                  only the alphabetical ordering of improper dihedrals.  No
//                  parameters are cloned.      
//
// Arguments:                                                           
//   mp:      the fitting data (contains a list of all parameters)      
//-----------------------------------------------------------------------------
void RecastAtomTypes(prmset *mp)
{
  int i, j, k;
  char *origt, *newt;
  prmtop *tp;

  // Initialize the change log
  mp->nchng = 0;
  mp->ChangeLog = CreateCmat(32, 80);

  // Change the names of the atom type in all circumstances
  for (i = 0; i < mp->nrecast; i++) {

    // Warn the user what we are about to do
    printf("mdgx >> Recasting atom type %.2s as %.2s.\n", mp->recast[i].orig,
           mp->recast[i].pnew);
    printf("mdgx >>\nmdgx >> In order to use the resulting parameter / frcmod "
           "file, change the\nmdgx >> declaration of type %.2s in your leaprc "
           "file and replace all instances\nmdgx >> of this atom type in the "
           "associated library files.\nmdgx >>\n", mp->recast[i].orig);

    // Loop over all atoms, bonds, angles, and dihedrals
    origt = mp->recast[i].orig;
    newt = mp->recast[i].pnew;
    for (j = 0; j < mp->natom; j++) {
      ReplaceAtomType(mp->atoms[j].atype, origt, newt);
      RecastInComment(mp->atoms[j].comment, mp->recast[i]);
    }
    for (j = 0; j < mp->nbond; j++) {
      ReplaceAtomType(mp->bonds[j].atype, origt, newt);
      ReplaceAtomType(mp->bonds[j].btype, origt, newt);
      RecastInComment(mp->bonds[j].comment, mp->recast[i]);
    }
    for (j = 0; j < mp->nangl; j++) {
      ReplaceAtomType(mp->angls[j].atype, origt, newt);
      ReplaceAtomType(mp->angls[j].btype, origt, newt);
      ReplaceAtomType(mp->angls[j].ctype, origt, newt);
      RecastInComment(mp->angls[j].comment, mp->recast[i]);
    }
    for (j = 0; j < mp->ntor; j++) {
      ReplaceAtomType(mp->torsions[j].atype, origt, newt);
      ReplaceAtomType(mp->torsions[j].btype, origt, newt);
      ReplaceAtomType(mp->torsions[j].ctype, origt, newt);
      ReplaceAtomType(mp->torsions[j].dtype, origt, newt);
      RecastInComment(mp->torsions[j].comment, mp->recast[i]);
    }

    // Repair any impropers if we messed up the alphabetical order
    for (j = 0; j < mp->ntor; j++) {
      if (mp->torsions[j].impr == 1) {
        AlphabetizeImproper(&mp->torsions[j]);
      }
    }

    // Loop over all topologies and change atom type names
    for (j = 0; j < mp->nunisys; j++) {
      tp = &mp->tpencyc[i];
      for (k = 0; k < tp->natom; k++) {
        if (strncmp(&tp->AtomTypes[4*k], origt, 4) == 0) {
          strncpy(&tp->AtomTypes[4*k], newt, 4);
          RecordAtomTypeChange(mp, tp, k, origt, newt, "Recast", " ");
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// RepackWildcards: repack wildcard terms at the front of a set of dihedrals.
//
// Arguments:                                                           
//   tset:      the list of torsion terms                               
//   tbuff:     buffer for sorting torsion terms                        
//   nt:        the length of tset                                      
//-----------------------------------------------------------------------------
static void RepackWildcards(torterm *tset, torterm* tbuff, int nt)
{
  int i, nwild, nplain;

  nwild = 0;
  nplain = 0;
  for (i = 0; i < nt; i++) {
    if (str4cmp(tset[i].atype, "X   ") != 0 &&
        str4cmp(tset[i].btype, "X   ") != 0 &&
        str4cmp(tset[i].ctype, "X   ") != 0 &&
        str4cmp(tset[i].dtype, "X   ") != 0) {
      tbuff[nplain] = tset[i];
      nplain++;
    }
    else {
      tset[nwild] = tset[i];
      nwild++;
    }
  }
  for (i = 0; i < nplain; i++) {
    tset[nwild+i] = tbuff[i];
  }
}

//-----------------------------------------------------------------------------
// PermuteBondBranches: run through all permutations of the branching atom
//                      type in a bond.
//
// Arguments:                                                           
//   mp:        the parameter set                                       
//   bondid:    the ID of the bond we're cloning                        
//   branchid:  the ID of the branch command governing the cloning      
//-----------------------------------------------------------------------------
static void PermuteBondBranches(prmset *mp, int bondid, int branchid)
{
  int i, j;
  int cclim[2];
  xbonddef *tbond;
  typebranch *tbranch;

  tbond = &mp->bonds[bondid];
  tbranch = &mp->cleave[branchid];
  cclim[0] = 1 - str4cmp(tbond->atype, tbranch->orig);
  cclim[1] = 1 - str4cmp(tbond->btype, tbranch->orig);
  for (i = 0; i <= cclim[0]; i++) {
    for (j = 0; j <= cclim[1]; j++) {
      if (i + j == 0) {
        continue;
      }

      // Make a new copy of this bond
      mp->bonds[mp->nbond] = *tbond;
      mp->bonds[mp->nbond].dup = 1;
      if (i == 1) {
        ReplaceAtomType(mp->bonds[mp->nbond].atype, tbranch->orig,
                        tbranch->pnew);
      }
      if (j == 1) {
        ReplaceAtomType(mp->bonds[mp->nbond].btype, tbranch->orig,
                        tbranch->pnew);
      }
      mp->bonds[mp->nbond].comment = (char*)malloc(MAXLINE*sizeof(char));
      sprintf(mp->bonds[mp->nbond].comment, "Branched from %.2s-%.2s",
              tbond->atype, tbond->btype);
      mp->nbond += 1;
    }
  }
}

//-----------------------------------------------------------------------------
// PermuteAngleBranches: run through all permutations of the branching atom
//                       type in an angle.                         
//
// Arguments:                                                           
//   mp:        the parameter set                                       
//   anglid:    the ID of the bond we're cloning                        
//   branchid:  the ID of the branch command governing the cloning      
//-----------------------------------------------------------------------------
static void PermuteAngleBranches(prmset *mp, int anglid, int branchid)
{
  int i, j, k;
  int cclim[3];
  xangldef *tangl;
  typebranch *tbranch;

  tangl = &mp->angls[anglid];
  tbranch = &mp->cleave[branchid];
  cclim[0] = 1 - str4cmp(tangl->atype, tbranch->orig);
  cclim[1] = 1 - str4cmp(tangl->btype, tbranch->orig);
  cclim[2] = 1 - str4cmp(tangl->ctype, tbranch->orig);
  for (i = 0; i <= cclim[0]; i++) {
    for (j = 0; j <= cclim[1]; j++) {
      for (k = 0; k <= cclim[2]; k++) {
        if (i + j + k == 0) {
          continue;
        }

        // Make a new copy of this bond
        mp->angls[mp->nangl] = *tangl;
        mp->angls[mp->nangl].dup = 1;
        if (i == 1) {
          ReplaceAtomType(mp->angls[mp->nangl].atype, tbranch->orig,
                          tbranch->pnew);
        }
        if (j == 1) {
          ReplaceAtomType(mp->angls[mp->nangl].btype, tbranch->orig,
                          tbranch->pnew);
        }
        if (k == 1) {
          ReplaceAtomType(mp->angls[mp->nangl].ctype, tbranch->orig,
                          tbranch->pnew);
        }
        mp->angls[mp->nangl].comment = (char*)malloc(MAXLINE*sizeof(char));
        sprintf(mp->angls[mp->nangl].comment, "Branched from %.2s-%.2s-%.2s",
                tangl->atype, tangl->btype, tangl->ctype);
        mp->nangl += 1;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// PermuteTorsionBranches: run through all permutations of the branching atom
//                         type in an angle.                       
//
// Arguments:                                                           
//   mp:        the parameter set                                       
//   torid:     the ID of the bond we're cloning                        
//   branchid:  the ID of the branch command governing the cloning      
//-----------------------------------------------------------------------------
static void PermuteTorsionBranches(prmset *mp, int torid, int branchid)
{
  int i, j, k, m;
  int cclim[4];
  torterm *ttor;
  typebranch *tbranch;

  ttor = &mp->torsions[torid];
  tbranch = &mp->cleave[branchid];
  cclim[0] = 1 - str4cmp(ttor->atype, tbranch->orig);
  cclim[1] = 1 - str4cmp(ttor->btype, tbranch->orig);
  cclim[2] = 1 - str4cmp(ttor->ctype, tbranch->orig);
  cclim[3] = 1 - str4cmp(ttor->dtype, tbranch->orig);
  for (i = 0; i <= cclim[0]; i++) {
    for (j = 0; j <= cclim[1]; j++) {
      for (k = 0; k <= cclim[2]; k++) {
        for (m = 0; m <= cclim[3]; m++) {
          if (i + j + k + m == 0) {
            continue;
          }

          // Make a new copy of this bond
          mp->torsions[mp->ntor] = *ttor;
          mp->torsions[mp->ntor].dup = 1;
          if (i == 1) {
            ReplaceAtomType(mp->torsions[mp->ntor].atype, tbranch->orig,
                            tbranch->pnew);
          }
          if (j == 1) {
            ReplaceAtomType(mp->torsions[mp->ntor].btype, tbranch->orig,
                            tbranch->pnew);
          }
          if (k == 1) {
            ReplaceAtomType(mp->torsions[mp->ntor].ctype, tbranch->orig,
                            tbranch->pnew);
          }
          if (m == 1) {
            ReplaceAtomType(mp->torsions[mp->ntor].dtype, tbranch->orig,
                            tbranch->pnew);
          }
          mp->torsions[mp->ntor].comment = (char*)malloc(MAXLINE*sizeof(char));
          sprintf(mp->torsions[mp->ntor].comment,
                  "Branched from %.2s-%.2s-%.2s-%.2s", ttor->atype,
                  ttor->btype, ttor->ctype, ttor->dtype);
          mp->ntor += 1;
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// CleaveAtomTypes: rename particular instances of an atom type in the  
//                  parameter set, clone parameters, and prepare to retain
//                  cloned parameters (with independently fitted values) in
//                  the output.                              
//                                                                      
// Arguments:                                                           
//   mp:      the fitting data (contains a list of all parameters)      
//-----------------------------------------------------------------------------
void CleaveAtomTypes(prmset *mp)
{
  int i, j, k, origid, ncopy, nimpr, ndihe;
  int* atmmask;
  char *origt, *newt;
  torterm* torbuff;
  prmtop *tp;

  // Set duplication flags
  for (i = 0; i < mp->natom; i++) {
    mp->atoms[i].dup = 0;
  }
  for (i = 0; i < mp->nbond; i++) {
    mp->bonds[i].dup = 0;
  }
  for (i = 0; i < mp->nangl; i++) {
    mp->angls[i].dup = 0;
  }
  for (i = 0; i < mp->ntor; i++) {
    mp->torsions[i].dup = 0;
  }
  for (i = 0; i < mp->nhb1012; i++) {
    mp->hb1012[i].dup = 0;
  }

  // Loop over all atom cleaves
  for (i = 0; i < mp->ncleave; i++) {
    origt = mp->cleave[i].orig;
    newt = mp->cleave[i].pnew;

    // Make sure this type does not exist  
    // already, then create a new atom type
    origid = -1;
    for (j = 0; j < mp->natom; j++) {
      if (str4cmp(mp->atoms[j].atype, mp->cleave[i].pnew) == 0) {
        printf("mdgx >> Error.  Branching requested in type %.2s to create "
               "type %.2s, but\nmdgx >> the new type already exists.\n",
               mp->cleave[i].orig, mp->cleave[i].pnew);
        exit(1);
      }
      if (str4cmp(mp->atoms[j].atype, mp->cleave[i].orig) == 0) {
        origid = j;
      }
    }
    if (origid == -1) {
      printf("mdgx >> Error.  Branching requested in type %.2s to create "
             "type %.2s, but\nmdgx >> the original type does not exist.\n",
             mp->cleave[i].orig, mp->cleave[i].pnew);
    }
    mp->atoms = (xatomdef*)realloc(mp->atoms, (mp->natom+1)*sizeof(xatomdef));
    mp->atoms[mp->natom] = mp->atoms[origid];
    mp->atoms[mp->natom].dup = 1;
    strncpy(mp->atoms[mp->natom].atype, mp->cleave[i].pnew, 4);
    mp->atoms[mp->natom].comment = (char*)malloc(MAXLINE*sizeof(char));
    sprintf(mp->atoms[mp->natom].comment, "Branched from type %.2s in "
            "'%s'", mp->cleave[i].orig, mp->cleave[i].instances);
    mp->natom += 1;

    // Duplicate bond terms
    ncopy = 0;
    for (j = 0; j < mp->nbond; j++) {
      if (str4cmp(mp->bonds[j].atype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->bonds[j].btype, mp->cleave[i].orig) == 0) {
        ncopy += 3;
      }
    }
    mp->bonds = (xbonddef*)realloc(mp->bonds, (mp->nbond+ncopy) *
                                   sizeof(xbonddef));
    origid = mp->nbond;
    for (j = 0; j < origid; j++) {
      if (str4cmp(mp->bonds[j].atype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->bonds[j].btype, mp->cleave[i].orig) == 0) {
        PermuteBondBranches(mp, j, i);
      }
    }

    // Duplicate angle terms
    ncopy = 0;
    for (j = 0; j < mp->nangl; j++) {
      if (str4cmp(mp->angls[j].atype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->angls[j].btype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->angls[j].ctype, mp->cleave[i].orig) == 0) {
        ncopy += 7;
      }
    }
    mp->angls = (xangldef*)realloc(mp->angls, (mp->nangl+ncopy) *
                                   sizeof(xangldef));
    origid = mp->nangl;
    for (j = 0; j < origid; j++) {
      if (str4cmp(mp->angls[j].atype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->angls[j].btype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->angls[j].ctype, mp->cleave[i].orig) == 0) {
        PermuteAngleBranches(mp, j, i);
      }
    }

    // Duplicate dihedral terms
    ncopy = 0;
    for (j = 0; j < mp->ntor; j++) {
      if (str4cmp(mp->torsions[j].atype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->torsions[j].btype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->torsions[j].ctype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->torsions[j].dtype, mp->cleave[i].orig) == 0) {
        ncopy += 15;
      }
    }
    mp->torsions = (torterm*)realloc(mp->torsions, (mp->ntor+ncopy) *
                                     sizeof(torterm));
    origid = mp->ntor;
    for (j = 0; j < origid; j++) {
      if (str4cmp(mp->torsions[j].atype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->torsions[j].btype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->torsions[j].ctype, mp->cleave[i].orig) == 0 ||
          str4cmp(mp->torsions[j].dtype, mp->cleave[i].orig) == 0) {
        PermuteTorsionBranches(mp, j, i);
      }
    }
  }

  // Repair any impropers if we messed up the alphabetical order
  for (i = 0; i < mp->ntor; i++) {
    if (mp->torsions[i].impr == 1) {
      AlphabetizeImproper(&mp->torsions[i]);
    }
  }

  // Restore the order of dihedrals, the impropers
  // in the torsions array, then place wildcards  
  // ahead of specific terms.                     
  torbuff = (torterm*)malloc(mp->ntor*sizeof(torterm));
  ndihe = 0;
  nimpr = 0;
  for (i = 0; i < mp->ntor; i++) {
    if (mp->torsions[i].impr == 0) {
      mp->torsions[ndihe] = mp->torsions[i];
      ndihe++;
    }
    else {
      torbuff[nimpr] = mp->torsions[i];
      nimpr++;
    }
  }
  for (i = 0; i < nimpr; i++) {
    mp->torsions[ndihe+i] = torbuff[i];
  }
  mp->ndihe = ndihe;
  mp->nimpr = nimpr;
  mp->ntor = mp->ndihe + mp->nimpr;
  RepackWildcards(mp->torsions, torbuff, ndihe);
  RepackWildcards(&mp->torsions[ndihe], torbuff, nimpr);

  // Parse contexts of each atom type change,   
  // loop over all topologies, apply the changes
  for (i = 0; i < mp->nunisys; i++) {
    tp = &mp->tpencyc[i];
    for (j = 0; j < mp->ncleave; j++) {
      origt = mp->cleave[j].orig;
      newt = mp->cleave[j].pnew;

      // Loop over all conformations of this system.  
      // The mask will be built when the first example
      // of the system is encountered.                
      for (k = 0; k < mp->nconf; k++) {
        if (mp->conf[k].GroupNum == i) {
          atmmask = ParseAmbMask(mp->cleave[j].instances, tp,
                                 &mp->conf[k].crd);
          break;
        }
      }
      for (k = 0; k < tp->natom; k++) {
        if (atmmask[k] == 1) {
          if (str4cmp(&tp->AtomTypes[4*k], origt) != 0) {
            printf("mdgx >> Error.  The atoms specified by ambmask %s\n"
                   "mdgx >> contain atoms of type %.2s when they should all "
                   "be of type %.2s.\nmdgx >> Fix the atom type branching and "
                   "rerun mdgx.\n", mp->cleave[j].instances,
                   &tp->AtomTypes[4*k], origt);
            exit(1);
          }
          strncpy(&tp->AtomTypes[4*k], newt, 4);
          RecordAtomTypeChange(mp, tp, k, origt, newt, "Branch", 
                               mp->cleave[j].instances);
        }
      }
      free(atmmask);
    }
  }

  // Free allocated memory
  free(torbuff);
}

//-----------------------------------------------------------------------------
// CreateNMROperation: initialize settings and allocate memory for an NMR-like
//                     operation as used by parameter fitting.  This is not the
//                     same as a NMR-like operation used by configuration
//                     sampling.
//
// Arguments:
//   thisop:    pointer to the operation to allocate
//-----------------------------------------------------------------------------
void CreateNMROperation(nmroper *thisop)
{
  int i;

  for (i = 0; i < 6; i++) {
    thisop->r1[i] = 0.0;
    thisop->r2[i] = 0.0;
    thisop->r3[i] = 0.0;
    thisop->r4[i] = 0.0;
    thisop->rk2[i] = 0.0;
    thisop->rk3[i] = 0.0;
  }
  thisop->target = 0.0;
  thisop->amask = CreateCmat(1, MAXNAME);
  thisop->bmask = CreateCmat(1, MAXNAME);
  thisop->cmask = CreateCmat(1, MAXNAME);
  thisop->dmask = CreateCmat(1, MAXNAME);
  thisop->style = (char*)malloc(MAXNAME*sizeof(char));
  thisop->label = (char*)malloc(MAXNAME*sizeof(char));
  thisop->label[0] = '\0';
  sprintf(thisop->style, "halfcup");
  thisop->stylecode = 0;
  thisop->nsibling = 1;
  thisop->usetypes = 0;
}

//-----------------------------------------------------------------------------
// DestroyNMROperation: free memory associated with a parameter fitting
//                      NMR restraint potential.
//
// Arguments:
//   thisop:    pointer to the operation to allocate
//-----------------------------------------------------------------------------
void DestroyNMROperation(nmroper *thisop)
{
  DestroyCmat(&thisop->amask);
  DestroyCmat(&thisop->bmask);
  DestroyCmat(&thisop->cmask);
  DestroyCmat(&thisop->dmask);
  free(thisop->style);
  free(thisop->label);
}

//-----------------------------------------------------------------------------
// CompareOperations: a function for taking the match matrice of two NMR
//                    operations across all systems and thereby finding
//                    whether they are truly identical or not.  This function
//                    works much like strcmp() in that it returns 0 if there
//                    are no differences, 1 if the two are different and
//                    distinct, and 2 if one of the operations subsumes the
//                    other.
//
// Arguments:
//   op[A,B]:    the two operations to compare
//   mp:         fitting data (contains conformations and topologies)
//   numerics:   flag to have numbers (not just labels or atom masks) checked
//-----------------------------------------------------------------------------
int CompareOperations(nmroper *opA, nmroper *opB, prmset *mp, int numerics)
{
  int i, j, k, nrow, nmatchA, nmatchB, order, found, mval;
  imat matchA, matchB;

  // If the labels are identical, this means that the two operations should
  // be taken as the same irrespective of their other properties.  Return an
  // unambiguous value to signify this case.
  if (strcmp(opA->label, opB->label) == 0) {
    return 4;
  }

  // Proceed to test the simples things first, then expand to more and 
  // more complex comparisons until every aspect of the operations has
  // been compared.  Return 0 if the operations are identical, 1 if not.
  if (numerics == 1) {

    // Basic numerical checks
    if (opA->usetypes != opB->usetypes ||
        opA->order != opB->order || opA->stylecode != opB->stylecode) {
      return 1;
    }

    // Now check the knots of the basis functions--could they look the same?
    for (i = 0; i < 6; i++) {
      if (fabs(opA->r1[i] - opB->r1[i]) > 1.0e-8 ||
          fabs(opA->r2[i] - opB->r2[i]) > 1.0e-8 ||
          fabs(opA->r3[i] - opB->r3[i]) > 1.0e-8 ||
          fabs(opA->r4[i] - opB->r4[i]) > 1.0e-8) {
        return 1;
      }
    }
  }
  order = opA->order;

  // The atom mask strings offer a chance to identify a match
  // before having to do even harder things.
  if ((order == 2 &&
       ((CompareCmat(&opA->amask, &opB->amask) == 0 &&
         CompareCmat(&opA->bmask, &opB->bmask) == 0) ||
        (CompareCmat(&opA->amask, &opB->bmask) == 0 &&
         CompareCmat(&opA->bmask, &opB->amask) == 0))) ||
      (order == 3 &&
       ((CompareCmat(&opA->amask, &opB->amask) == 0 &&
         CompareCmat(&opA->bmask, &opB->bmask) == 0 &&
         CompareCmat(&opA->cmask, &opB->cmask) == 0) ||
        (CompareCmat(&opA->amask, &opB->cmask) == 0 &&
         CompareCmat(&opA->bmask, &opB->bmask) == 0 &&
         CompareCmat(&opA->cmask, &opB->amask) == 0))) ||
      (order == 4 &&
       ((CompareCmat(&opA->amask, &opB->amask) == 0 &&
         CompareCmat(&opA->bmask, &opB->bmask) == 0 &&
         CompareCmat(&opA->cmask, &opB->cmask) == 0 &&
         CompareCmat(&opA->dmask, &opB->dmask) == 0) ||
        (CompareCmat(&opA->amask, &opB->dmask) == 0 &&
         CompareCmat(&opA->bmask, &opB->cmask) == 0 &&
         CompareCmat(&opA->cmask, &opB->bmask) == 0 &&
         CompareCmat(&opA->dmask, &opB->amask) == 0)))) {
    return 0;
  }

  // If it's still unclear, then we must loop over all systems and
  // look at where the NMR operations will apply in order to decide
  // whether they are actually the same or different.
  for (i = 0; i < mp->nunisys; i++) {
    matchA = BuildMatchMatrix(&mp->conf[mp->FirstConf[i]], opA);
    matchB = BuildMatchMatrix(&mp->conf[mp->FirstConf[i]], opB);
    mval = DistinguishImat(&matchA, &matchB);
    DestroyImat(&matchA);
    DestroyImat(&matchB);
    if (mval > 0) {
      return mval;
    }
  }

  // If this point has been reached, complete correspondence was found
  return 0;
}

//-----------------------------------------------------------------------------
// ReadNMROperationsFile: read the NMR operations from a separate parameter
//                        file, with a unique format.  The operations specified
//                        in this file can and will be written in a format that
//                        the standard Amber packages can read, but there are
//                        additional details in this file format that make it
//                        convenient for users to specify operations that will
//                        do more complex things.
//
// Arguments:
//   mp:      the fitting data (contains a list of all parameters and the
//            name of the additional input file)
//-----------------------------------------------------------------------------
void ReadNMROperationsFile(prmset *mp)
{
  int i, j, opcount, na, nb, nc, nd;
  int* unique;
  double rbegin, rend, rkval;
  char* direction;
  imat compval;
  cmat cfi, Lwords, atomnames, glossary;
  nmlgroup nmlg;
  nmroper *myop;

  // Bail out if there is nothing to do
  if (mp->NMROpsFile[0] == '\0') {
    mp->nops = 0;
    return;
  }

  // Read the file verbatim into memory,
  // then extract the proper namelists
  nmlg = ReadNamelistsFromFile(mp->NMROpsFile, 1, "nmropt");

  // Allocate memory for NMR restraint operations
  // in the context of data fitting
  mp->nmrops = (nmroper*)malloc(nmlg.count*sizeof(nmroper));
  direction = (char*)malloc(MAXNAME*sizeof(char));
  glossary = CreateCmat(32, 64);
  AddToGlossary(&glossary, 10, "Atom1", "atm1", "Atom2", "atm2", "Atom3",
                "atm3", "Atom4", "atm4", "Atoms", "iatm");
  AddToGlossary(&glossary, 12, "Operation", "style", "PenaltyGrows", "growing",
                "r1", "r2", "r3", "r4",        "StartParabolic", "begin",
                "StopParabolic", "end");
  AddToGlossary(&glossary, 2, "Stiffness", "K");
  opcount = 0;
  for (i = 0; i < nmlg.count; i++) {

    // Defaults for this restraint
    atomnames = CreateCmat(4, MAXNAME);
    direction[0] = '\0';
    rbegin = 0.0;
    rend = 0.0;
    myop = &mp->nmrops[opcount];

    // Allocate needed memory
    CreateNMROperation(myop);

    // Read input data
    SeekReal(nmlg.nml[i], &myop->r1[0], "r1", "r1");
    SeekReal(nmlg.nml[i], &myop->r2[0], "r2", "r2");
    SeekReal(nmlg.nml[i], &myop->r3[0], "r3", "r3");
    SeekReal(nmlg.nml[i], &myop->r4[0], "r4", "r4");
    SeekReal(nmlg.nml[i], &myop->rk2[0], "rk2", "rk2");
    SeekReal(nmlg.nml[i], &myop->rk3[0], "rk3", "rk3");
    SeekReal(nmlg.nml[i], &rkval, "Stiffness", "K");
    SeekReal(nmlg.nml[i], &rbegin, "StartParabolic", "begin");
    SeekReal(nmlg.nml[i], &rend, "StopParabolic", "end");
    SeekMultiString(nmlg.nml[i], &myop->amask, "Atom1", "atm1", &glossary, 1);
    SeekMultiString(nmlg.nml[i], &myop->bmask, "Atom2", "atm2", &glossary, 1);
    SeekMultiString(nmlg.nml[i], &myop->cmask, "Atom3", "atm3", &glossary, 1);
    SeekMultiString(nmlg.nml[i], &myop->dmask, "Atom4", "atm4", &glossary, 1);
    SeekMultiString(nmlg.nml[i], &atomnames, "Atoms", "iatm", &glossary, 1);
    SeekString(nmlg.nml[i], myop->style, "Operation", "style");
    SeekString(nmlg.nml[i], direction, "PenaltyGrows", "growing");
    SeekString(nmlg.nml[i], myop->label, "Label", "label");

    // Check if there are multiple copies
    na = 0;
    for (j = 0; j < myop->amask.row; j++) {
      if (myop->amask.map[j][0] != '\0') {
        na++;
      }
    }
    nb = 0;
    for (j = 0; j < myop->bmask.row; j++) {
      if (myop->bmask.map[j][0] != '\0') {
        nb++;
      }
    }
    nc = 0;
    for (j = 0; j < myop->cmask.row; j++) {
      if (myop->cmask.map[j][0] != '\0') {
        nc++;
      }
    }
    nd = 0;
    for (j = 0; j < myop->dmask.row; j++) {
      if (myop->dmask.map[j][0] != '\0') {
        nd++;
      }
    }
    j = (na > nb) ? na : nb;
    j = (j > nc) ? j : nc;
    j = (j > nd) ? j : nd;
    if (j > 0) {
      myop->nsibling = j;
    }
    myop->amask = ReallocCmat(&myop->amask, myop->nsibling, MAXNAME);
    myop->bmask = ReallocCmat(&myop->bmask, myop->nsibling, MAXNAME);
    myop->cmask = ReallocCmat(&myop->cmask, myop->nsibling, MAXNAME);
    myop->dmask = ReallocCmat(&myop->dmask, myop->nsibling, MAXNAME);

    // Fill in atom names if they were specified in a particular way
    if (myop->amask.map[0][0] == '\0' && atomnames.row >= 1 &&
        atomnames.map[0][0] != '\0') {
      strcpy(myop->amask.map[0], atomnames.map[0]);
    }
    if (myop->bmask.map[0][0] == '\0' && atomnames.row >= 2) {
      strcpy(myop->bmask.map[0], atomnames.map[1]);
    }
    if (myop->cmask.map[0][0] == '\0' && atomnames.row >= 3) {
      strcpy(myop->cmask.map[0], atomnames.map[2]);
    }
    if (myop->dmask.map[0][0] == '\0' && atomnames.row >= 4) {
      strcpy(myop->dmask.map[0], atomnames.map[3]);
    }
    if (myop->amask.map[0][0] == '\0' || myop->bmask.map[0][0] == '\0') {
      DestroyNMROperation(myop);
      continue;
    }
    myop->order = (myop->cmask.map[0][0] == '\0') ? 2 :
                  (myop->dmask.map[0][0] == '\0') ? 3 : 4;

    // Did the user specify a standard NMR restraint?
    if ((fabs(myop->rk2[i]) > 1.0e-8 || fabs(myop->rk3[i]) > 1.0e-8) &&
             (myop->r4[i] >= myop->r3[i] && myop->r3[i] >= myop->r2[i] &&
              myop->r2[i] >= myop->r1[i]) &&
             (myop->r4[i] - myop->r3[i] >= 1.0e-4 ||
              myop->r2[i] - myop->r1[i] > 1.0e-4)) {
      sprintf(myop->style, "standard");
      myop->stylecode = 0;
      myop->nlayer = 1;
      opcount++;
    }

    // Did the user specify a half-parabolic restraint
    // using the vernacular keywords?
    else if (CaselessStrcmp(myop->style, "halfcup", -1) == 0 ||
             CaselessStrcmp(myop->style, "halfparabola", -1) == 0) {
      if (direction[0] != '\0') {
        if (CaselessStrcmp(direction, "right", -1) == 0 &&
            fabs(rend-rbegin) > 1.0e-8) {
          myop->r3[0] = rbegin;
          myop->r4[0] = rend;
          myop->rk2[0] = 0.0;
          if (fabs(myop->rk3[0]) < 1.0e-8) {
            myop->rk3[0] = 1.0;
          }
        }
        if (CaselessStrcmp(direction, "left", -1) == 0 &&
            fabs(rend-rbegin) > 1.0e-8) {
          myop->r2[0] = rbegin;
          myop->r1[0] = rend;
          if (fabs(myop->rk2[0]) < 1.0e-8) {
            myop->rk2[0] = 1.0;
          }
          myop->rk3[0] = 0.0;
        }
      }
      else if (fabs(rbegin) > 1.0e-8 && fabs(rend) > 1.0e-8 &&
               fabs(rend - rbegin) > 1.0e-4) {
        if (rbegin > rend + 1.0e-8) {
          myop->r2[0] = rbegin;
          myop->r1[0] = rend;
          myop->rk2[0] = 0.0;
          if (fabs(myop->rk3[0]) < 1.0e-8) {
            myop->rk3[0] = 1.0;
          }
        }
        else if (rend > rbegin + 1.0e-8) {
          myop->r3[0] = rbegin;
          myop->r4[0] = rend;
          myop->rk3[0] = 0.0;
          if (fabs(myop->rk2[0]) < 1.0e-8) {
            myop->rk2[0] = 1.0;
          }
        }
      }
      else if (myop->r4[0] >= myop->r3[0] && myop->r3[0] >= myop->r2[0] &&
               myop->r2[0] >= myop->r1[0]) {
        if (myop->r4[0] - myop->r3[0] >= 1.0e-4 &&
            myop->r2[0] - myop->r1[0] < 1.0e-4 &&
            fabs(myop->rk3[0]) < 1.0e-8) {
          myop->rk3[0] = 1.0;
        }
        if (myop->r4[0] - myop->r3[0] < 1.0e-4 &&
            myop->r2[0] - myop->r1[0] >= 1.0e-4 &&
            fabs(myop->rk2[0]) < 1.0e-8) {
          myop->rk2[0] = 1.0;
        }
      }

      // Check the input: at least r1 and r2, or r3 and r4, must be zero while
      // the others are not.  Increment the operation count if successful.
      if ((fabs(myop->r1[0]) < 1.0e-8 && fabs(myop->r2[0]) < 1.0e-8 &&
           fabs(myop->r4[0]) > 1.0e-8 && myop->r3[0] < myop->r4[0]) ||
          (fabs(myop->r1[0]) > 1.0e-8 && myop->r2[0] > myop->r1[0] &&
           fabs(myop->r3[0]) < 1.0e-8 && fabs(myop->r4[0]) < 1.0e-8)) {
        sprintf(myop->style, "halfcup");
        myop->stylecode = 1;
        myop->nlayer = 1;
        opcount++;
      }
      else {
        DestroyNMROperation(myop);
      }
    }

    // Did the user specify an S-function restraint (really a combination of
    // NMR restraints) using the vernacular keywords?
    else if (CaselessStrcmp(myop->style, "sfunc", -1) == 0 ||
             CaselessStrcmp(myop->style, "riser", -1) == 0) {

      // Check the validity of the input
      if (myop->r4[0] < myop->r3[0] && myop->r3[0] <= myop->r2[0] &&
          myop->r2[0] < myop->r1[0]) {
        SWAP(myop->r4[0], myop->r1[0], rend);
        SWAP(myop->r3[0], myop->r2[0], rend);
      }
      if (myop->r1[0] < myop->r2[0] && myop->r2[0] <= myop->r3[0] &&
          myop->r3[0] < myop->r4[0]) {

        // This is a bit tricky: the user will specify r1, r2, r3, and r4 in
        // increasing order along the coordinate of interest (if the user
        // gives those quantities in decreasing order, that's fine too, they
        // will just be rearranged).  However, for an S-function, there are
        // really two half-parabola NMR restraints, both increasing as the
        // coordinate goes right.  But, the user is supposed to be able to
        // order up an S-function with just one NMR namelist.  The first
        // half parabola increases in magnitude (it does not decrease) from
        // r1 to r2.  The second half parabola increases in magnitude from
        // r3 to r4 BUT has opposite sign to the first. 
        myop->rk2[0] = 1.0;
        myop->rk3[0] = -(myop->r2[0] - myop->r1[0]) /
                        (myop->r4[0] - myop->r3[0]);
        sprintf(myop->style, "sfunc");
        myop->stylecode = 2;
        myop->nlayer = 2;
        opcount++;
      }
      else {
        DestroyNMROperation(myop);
      }
    }

    // The user didn't make any sense.
    else {
      DestroyNMROperation(myop);
    }

    // Free allocated memory
    DestroyCmat(&atomnames);
  }
  mp->nops = opcount;

  // Add labels if any operations are missing them
  for (i = 0; i < mp->nops; i++) {
    if (mp->nmrops[i].label[0] == '\0') {
      sprintf(mp->nmrops[i].label, "%d", i+1);
    }
  }

  // Prune duplicate operations
  unique = (int*)malloc(mp->nops*sizeof(int));
  SetIVec(unique, mp->nops, 1);
  compval = CreateImat(mp->nops, mp->nops);
  for (i = 0; i < mp->nops-1; i++) {
    for (j = i+1; j < mp->nops; j++) {
      compval.map[i][j] = CompareOperations(&mp->nmrops[i],
                                            &mp->nmrops[j], mp, 1);
      compval.map[j][i] = compval.map[i][j];
    }
  }

  // Now, several passes over the comparison matrix.  First we must
  // eliminate operations that are superceded by later ones.  This
  // is done by seeking out comparisons with a value of 3 (unique
  // applications of a later operation, when j > i), or comparisons
  // with a value of 0 (the operations are identical), or comparisons
  // with a value of 4 (the two operations have the same label).
  for (i = 0; i < mp->nops; i++) {
    for (j = i+1; j < mp->nops; j++) {
      if (compval.map[i][j] == 3 || compval.map[i][j] == 0 ||
          compval.map[i][j] == 4) {
        unique[i] = 0;
        break;
      }
    }
  }

  // Now, forward eliminations can proceed.  If the comparison got
  // a value of 2 (the first operation i, again i < j, had unique
  // applications), others down the line that it comprises may be
  // removed.
  for (i = 0; i < mp->nops; i++) {
    if (unique[i] == 0) {
      continue;
    }
    for (j = i+1; j < mp->nops; j++) {
      if (compval.map[i][j] == 2) {
        unique[j] = 0;
      }
    }
  }

  // Free memory associated with unused operations
  for (i = 0; i < mp->nops; i++) {
    if (unique[i] == 0) {
      DestroyNMROperation(&mp->nmrops[i]);
    }
  }

  // Repack the operations array
  j = 0;
  for (i = 0; i < mp->nops; i++) {
    if (unique[i] == 1) {
      mp->nmrops[j] = mp->nmrops[i];
      j++;
    }
  }
  mp->nops = j;

  // Make sure everything has a label, and convert
  // units from degrees to radians
  for (i = 0; i < mp->nops; i++) {
    if (mp->nmrops[i].order == 3 || mp->nmrops[i].order == 4) {
      for (j = 0; j < 6; j++) {
        mp->nmrops[i].r1[j] *= PI/180.0;
        mp->nmrops[i].r2[j] *= PI/180.0;
        mp->nmrops[i].r3[j] *= PI/180.0;
        mp->nmrops[i].r4[j] *= PI/180.0;
      }
    }
  }

  // Set the number of NMR operations in the parameter set
  mp->nmrops = (nmroper*)realloc(mp->nmrops, mp->nops*sizeof(nmroper));

  // Alert the user as to progress
  if (mp->verbose == 1) {
    printf("mdgx >> Read %d piecewise basis functions.\n", mp->nops);
  }

  // Free allocated memory
  DestroyImat(&compval);
  free(direction);
  free(unique);
}
