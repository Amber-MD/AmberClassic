#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Parse.h"
#include "mdgxVector.h"
#include "Matrix.h"
#include "Constraints.h"
#include "VirtualSites.h"
#include "Topology.h"
#include "CellManip.h"
#include "CrdManip.h"
#include "Macros.h"
#include "Debug.h"
#include "Manual.h"

//-----------------------------------------------------------------------------
// ExpandLJTables: expands the list of Lennard-Jones parameters to include a
//                 newly added type, then recompile any derivative tables.
//
// Arguments:
//   tp:      the system topology
//   sig:     sigma parameter to add to the LJ table
//   eps:     epsilon parameter to add to the LJ table
//-----------------------------------------------------------------------------
void ExpandLJTables(prmtop *tp, double sig, double eps)
{
  int h, i, j;
  imat S;
  double lja, ljb, tsig, teps;
  double* Vsig;
  double* Veps;

  // Get the current vectors of pure Lennard-Jones parameters
  Vsig = (double*)malloc((tp->ntypes+1)*sizeof(double));
  Veps = (double*)malloc((tp->ntypes+1)*sizeof(double));
  for (i = 0; i < tp->ntypes; i++) {
    lja = tp->LJA[tp->NBParmIdx[(tp->ntypes+1)*i]];
    ljb = tp->LJB[tp->NBParmIdx[(tp->ntypes+1)*i]];
    Vsig[i] = (lja > 1.0e-04 && ljb > 1.0e-04) ? pow(lja/ljb, 1.0/6.0) : 0.0;
    Veps[i] = (Vsig[i] > 1.0e-04) ? 0.25*ljb/pow(Vsig[i], 6.0) : 0.0;
  }
  Vsig[tp->ntypes] = sig;
  Veps[tp->ntypes] = eps;

  // The nonbonded parameter index can
  // be inscribed into a square matrix
  S = CreateImat(tp->ntypes+1, tp->ntypes+1);
  h = 0;
  for (i = 0; i < tp->ntypes; i++) {
    for (j = 0; j < tp->ntypes; j++) {
      S.map[i][j] = tp->NBParmIdx[h];
      h++;
    }
  }
  h = tp->ntypes*(tp->ntypes+1)/2;
  for (i = 0; i < tp->ntypes+1; i++) {
    S.map[tp->ntypes][i] = h;
    S.map[i][tp->ntypes] = h;
    h++;
  }
  tp->NBParmIdx = (int*)realloc(tp->NBParmIdx,
				(tp->ntypes+1)*(tp->ntypes+1)*sizeof(int));
  h = 0;
  for (i = 0; i < tp->ntypes+1; i++) {
    for (j = 0; j < tp->ntypes+1; j++) {
      tp->NBParmIdx[h] = S.map[i][j];
      h++;
    }
  }

  // The Lennard-Jones A and B coefficient
  // arrays can now be re-compiled        
  tp->LJA = (double*)realloc(tp->LJA,
			     (tp->ntypes+1)*(tp->ntypes+2)/2*sizeof(double));
  tp->LJB = (double*)realloc(tp->LJB,
			     (tp->ntypes+1)*(tp->ntypes+2)/2*sizeof(double));
  for (i = 0; i < tp->ntypes+1; i++) {
    for (j = 0; j < tp->ntypes+1; j++) {
      h = S.map[i][j];
      if (h >= 0) {
	tsig = 0.5*(Vsig[i]+Vsig[j]);
	teps = 4.0*sqrt(Veps[i]*Veps[j]);
	tp->LJA[h] = teps*pow(tsig, 12.0);
	tp->LJB[h] = teps*pow(tsig,  6.0);
      }
      else {
	h = -h;
	tp->LJA[h] = 0.0;
        tp->LJB[h] = 0.0;
      }
    }
  }

  // Increment the number of Lennard-Jones types
  tp->ntypes += 1;

  // Finally, recompile the force and energy tables
  DestroyDmat(&tp->LJftab);
  DestroyDmat(&tp->LJutab);
  OrderLJParameters(tp);

  // Free allocated memory
  free(Vsig);
  free(Veps);
  DestroyImat(&S);
}

//-----------------------------------------------------------------------------
// IgnoredEPKeyword: function to let the user know that a given keyword is not
//                   needed by a certain rule type and will be ignored.
//
// Arguments:
//   itest:      flag indicating that a keyword has been specified
//   keyword:    the keyword
//   alias:      keyword alias
//   frtype:     the frame type
//-----------------------------------------------------------------------------
static void IgnoredEPKeyword(int itest, char* keyword, char* alias,
                             char* frtype)
{
  char* msg;
  
  // If the offending keyword was never specified, there is nothing to do
  if (itest == 0) {
    return;
  }

  // Construct the message
  msg = (char*)malloc(MAXLINE * sizeof(char));
  sprintf(msg, "Warning.  A value for %s / %s was given, but frame type %s "
	  "takes no such keyword.  Value will be ignored.", keyword, alias,
	  frtype);
  PrintParagraph(msg, 79, "ReadEPRuleFile", stdout);
  
  // Free allocated memory
  free(msg);
}

//-----------------------------------------------------------------------------
// MissingEPKeyword: function to let the user know that a given keyword needed
//                   by a certain rule type is missing.
//
// Arguments:
//   itest:      flag indicating that a keyword has been specified
//   atest:      flag indicating that an alternative was supplied  
//   keyword:    the keyword
//   alias:      keyword alias
//   frtype:     the frame type
//-----------------------------------------------------------------------------
static void MissingEPKeyword(int itest, int atest, char* keyword, char* alias,
                             char* frtype)
{
  char* msg;
  
  // If the keyword was specified, there is nothing to do.
  // Check that any alternative keyword was also not supplied.
  if (itest == 1 || atest == 1) {
    return;
  }

  // Construct the message
  msg = (char*)malloc(MAXLINE * sizeof(char));
  sprintf(msg, "Error.  A value for %s / %s is required for frame type %s.",
	  keyword, alias, frtype);
  PrintParagraph(msg, 79, "ReadEPRuleFile", stdout);
  
  // Free allocated memory
  free(msg);

  // Haha, we are going to exit anyway.
  exit(1);
}

//-----------------------------------------------------------------------------
// MistakenEPKeyword: function to let the user know that a given keyword needed
//                    by a certain rule type is missing, but that a proxy has
//                    been specified that will suffice.
//
// Arguments:
//   itest:       flag indicating that a keyword has been specified
//   atest:       flag indicating that an alternative was supplied  
//   keyword:     the keyword
//   alias:       keyword alias
//   altkeyword:  the alternative keyword
//   altalias:    the alternative keyword alias
//   frtype:      the frame type
//-----------------------------------------------------------------------------
static void MistakenEPKeyword(int itest, int atest, char* keyword, char* alias,
                              char* altkeyword, char* altalias, char* frtype)
{
  char* msg;
  
  // If the offending keyword was never specified, there is nothing to do
  if (atest == 0) {
    return;
  }

  // Construct the message
  msg = (char*)malloc(MAXLINE * sizeof(char));

  // Case where the true keyword is unspecified
  if (itest == 0) {
    sprintf(msg, "Warning.  A value for %s / %s is required for frame type "
	    "%s, but a value for %s / %s was found instead.  The alternative "
	    "keyword will suffice to create this frame type, but may have "
	    "been specified with different units than was intended.", keyword,
	    alias, frtype, altkeyword, altalias);
  }

  // Case where both the keyword and its alternative are specified
  if (itest == 1) {
    sprintf(msg, "Warning.  A value for %s / %s is required for frame type "
	    "%s, but a value for %s / %s was also found.  The value given in "
	    "the latter of these two directives will be taken.\n", keyword,
            alias, frtype, altkeyword, altalias);
  }

  // Print the result
  PrintParagraph(msg, 79, "ReadEPRuleFile", stdout);
  
  // Free allocated memory
  free(msg);
}

//-----------------------------------------------------------------------------
// ReadEPRuleFile: read a file containing user-specified extra point rules.
//                 Rules specified in this file will override built-in rules
//                 within the program, such as four-point water rules.
//
// Arguments:
//   tp:      the system topology
//-----------------------------------------------------------------------------
void ReadEPRuleFile(prmtop *tp)
{
  int i, j, neprule, maxrule, collect, found;
  int v12given, v13given, v23given, v1egiven, thetagiven, epgiven, atomgiven;
  int v12x13given, parentgiven, frame2given, frame3given, frame4given;
  int qgiven, ljgiven;
  double realtmp, sig, eps, aval, bval;
  char line[MAXLINE], tmpepn[64], tmpres[64], tmpstyle[64];
  FILE *inp;
  cmat lwords, tmpfr;
  eprule *tmr;

  // Initialize variables relating to extra points
  tp->EPInserted = 0;
  tp->norigatom = tp->natom;

  // If there's no file specified, we're done.
  if (tp->eprulesource[0] == '\0') {
    tp->neprule = 0;
    return;
  }
  maxrule = 10;
  tp->eprules = (eprule*)malloc(maxrule*sizeof(eprule));
  if ((inp = fopen(tp->eprulesource, "r")) == NULL) {
    printf("ReadEPRuleFile >> Error.  Extra points rule file %s not found.\n",
           tp->eprulesource);
    exit(1);
  }
  neprule = 0;
  tmpfr = CreateCmat(4, MAXNAME);
  while (AdvanceToSegment(inp, "rule", 0) != 0) {
    collect = 1;

    // Prepare to record dimension specifications
    epgiven     = 0;
    atomgiven   = 0;
    v12given    = 0;
    v13given    = 0;
    v23given    = 0;
    v1egiven    = 0;
    thetagiven  = 0;
    v12x13given = 0;
    parentgiven = 0;
    frame2given = 0;
    frame3given = 0;
    frame4given = 0;
    qgiven      = 0;
    ljgiven     = 0;
    
    // Initialize this rule
    tmr = &tp->eprules[neprule];
    tmpfr.map[0][0] = '\0';
    tmpfr.map[1][0] = '\0';
    tmpfr.map[2][0] = '\0';
    tmpfr.map[3][0] = '\0';
    tmpepn[0] = '\0';
    tmpres[0] = '\0';
    tmpstyle[0] = '\0';
    tmr->frstyle = -1;
    tmr->excl2 = 0;
    tmr->excl3 = 0;
    tmr->excl4 = 0;
    tmr->d1 = 0.0;
    tmr->d2 = 0.0;
    tmr->d3 = 0.0;
    tmr->Charge = 1.1e30;
    tmr->sig = -1.0;
    tmr->eps = -1.0;
    tmr->Join = 0;
    tmr->Rotat = 0;
    tmr->modElec = 0;
    tmr->modLJ = 0;
    for (i = 0; i < 5; i++) {
      tmpepn[i] = '\0';
      tmpres[i] = '\0';
    }

    // The final initialization detail is the Lennard-Jones type.
    // By default, an extra point should have no Lennard-Jones
    // properties, but it is not enough to say -1 because if the
    // topology is written out this -1 will get converted to the
    // type +1 with whatever properties atom type index 1 entails
    // (and that's C indexing--Fortran indexing in the actual
    // file would increment the type to +2).  So search for a
    // type in the topology with zero Lennard-Jones sigma and
    // epsilon.  If none can be found, create a new index with
    // those properties and extend the tables.
    found = 0;
    for (i = 0; i < tp->ntypes; i++) {
      sig = 0.0;
      eps = 0.0;
      if (tp->NBParmIdx[i*tp->ntypes + i] >= 0) {
        aval = tp->LJA[tp->NBParmIdx[i*tp->ntypes + i]];
        bval = tp->LJB[tp->NBParmIdx[i*tp->ntypes + i]];
        if (aval > 1.0e-8) {
          sig = pow(aval/bval, 1.0/6.0);
          eps = bval/pow(sig, 6.0);
        }
      }
      if (sig < 1.0e-8 && eps < 1.0e-8) {
	tmr->LJIdx = i;
	found = 1;
	break;
      }
    }
    if (found == 0) {
      tmr->LJIdx = -(tp->ntypes);
      ExpandLJTables(tp, 0.0, 0.0);
    }
    while (collect == 1) {

      // Information about each rule is collected in much the same way  
      // as information for any of the namelists in the command section.
      collect = ReadNamelistLine(line, &lwords, "ReadEPRuleFile", inp);
      if (collect == 0) {
        continue;
      }
      SeekString(lwords, tmpepn, "ExtraPoint", "epname");
      SeekString(lwords, tmpepn, "AtomName", "atom");
      SeekInt(lwords, &tmr->frstyle, "FrameStyle", "style");
      SeekInt(lwords, &tmr->excl2, "Exclude2", "excl2");
      SeekInt(lwords, &tmr->excl3, "Exclude3", "excl3");
      SeekInt(lwords, &tmr->excl4, "Exclude4", "excl4");
      SeekString(lwords, tmpfr.map[0], "ParentAtom", "parent");
      SeekString(lwords, tmpfr.map[0], "FrameAtom1", "frame1");
      SeekString(lwords, tmpfr.map[1], "FrameAtom2", "frame2");
      SeekString(lwords, tmpfr.map[2], "FrameAtom3", "frame3");
      SeekString(lwords, tmpfr.map[3], "FrameAtom4", "frame4");
      SeekReal(lwords, &tmr->d1, "Vector12", "v12");
      SeekReal(lwords, &tmr->d1, "Vector1E", "v1e");
      SeekReal(lwords, &tmr->d2, "Vector13", "v13");
      SeekReal(lwords, &tmr->d2, "Theta", "theta");
      SeekReal(lwords, &tmr->d3, "Vector23", "v23");
      SeekReal(lwords, &tmr->d3, "Vector12x13", "v12x13");
      SeekReal(lwords, &tmr->Charge, "Charge", "q");
      SeekReal(lwords, &tmr->sig, "Sigma", "sig");
      SeekReal(lwords, &tmr->eps, "Epsilon", "eps");
      SeekString(lwords, tmpres, "ResidueName", "residue");
      SeekString(lwords, tmpstyle, "FrameType", "eptype");

      // Record which dimensions have been specified
      if (WordInGlossary(&lwords, "Vector12") == 1 ||
	  WordInGlossary(&lwords, "v12") == 1) {
	v12given = 1;
      }
      else if (WordInGlossary(&lwords, "Vector13") == 1 ||
               WordInGlossary(&lwords, "v13") == 1) {
	v13given = 1;
      }
      else if (WordInGlossary(&lwords, "Vector23") == 1 ||
               WordInGlossary(&lwords, "v23") == 1) {
        v23given = 1;
      }
      else if (WordInGlossary(&lwords, "Vector1E") == 1 ||
               WordInGlossary(&lwords, "v1e") == 1) {
       v1egiven = 1;
      }
      else if (WordInGlossary(&lwords, "Theta") == 1 ||
               WordInGlossary(&lwords, "theta") == 1) {
        thetagiven = 1;
      }
      else if (WordInGlossary(&lwords, "Vector12x13") == 1 ||
               WordInGlossary(&lwords, "v12x13") == 1) {
        v12x13given = 1;
      }
      else if (WordInGlossary(&lwords, "ParentAtom") == 1 ||
               WordInGlossary(&lwords, "parent") == 1) {
        parentgiven = 1;
      }
      else if (WordInGlossary(&lwords, "FrameAtom1") == 1 ||
               WordInGlossary(&lwords, "frame1") == 1) {
        parentgiven = 1;
      }
      else if (WordInGlossary(&lwords, "FrameAtom2") == 1 ||
               WordInGlossary(&lwords, "frame2") == 1) {
        frame2given = 1;
      }
      else if (WordInGlossary(&lwords, "FrameAtom3") == 1 ||
               WordInGlossary(&lwords, "frame3") == 1) {
        frame3given = 1;
      }
      else if (WordInGlossary(&lwords, "FrameAtom4") == 1 ||
               WordInGlossary(&lwords, "frame4") == 1) {
        frame4given = 1;
      }
      else if (WordInGlossary(&lwords, "ExtraPoint") == 1 ||
               WordInGlossary(&lwords, "epname") == 1) {
        epgiven = 1;
      }
      else if (WordInGlossary(&lwords, "AtomName") == 1 ||
               WordInGlossary(&lwords, "atom") == 1) {
        atomgiven = 1;
      }
      else if (WordInGlossary(&lwords, "Charge") == 1 ||
               WordInGlossary(&lwords, "q") == 1) {
        qgiven = 1;
      }
      else if (WordInGlossary(&lwords, "Sigma") == 1 ||
               WordInGlossary(&lwords, "sig") == 1) {
        ljgiven |= 0x1;
      }
      else if (WordInGlossary(&lwords, "Epsilon") == 1 ||
               WordInGlossary(&lwords, "eps") == 1) {
        ljgiven |= 0x10;
      }

      // Free allocated memory
      DestroyCmat(&lwords);
    }

    // Check that the residue and extra point names are properly specified.
    if (tmpepn[0] == '\0') {
      printf("ReadEPRuleFile >> Error.  Extra point name unspecified.  Use "
             "the ExtraPoint /\nReadEPRuleFile >> epname or AtomName / atom "
	     "keywords, as appropriate.\n");
      exit(1);
    }
    if (tmpres[0] == '\0') {
      printf("ReadEPRuleFile >> Error.  Residue name unspecified.  Use the "
             "ResidueName /\nReadEPRuleFile >> residue keyword.\n");
      exit(1);
    }

    // Check that this rule actually sets non-bonded properties.
    if (qgiven == 0 && ljgiven == 0) {
      printf("ReadEPRuleFile >> Error.  Rule specifies no non-bonded "
	     "properties.  Use the\nReadEPRuleFile >> Charge / q, Sigma / "
	     "sig, or Epsilon / eps keywords to\nReadEPRuleFile >> specify "
             "charge, sigma, or epsilon for this new particle.\n");
      exit(1);
    }
    if (ljgiven == 0x1) {
      printf("ReadEPRuleFile >> Error.  Lennard-Jones Sigma was specified, "
	     "but no Epsilon\nReadEPRuleFile >> was given.  Use the Epsilon / "
	     "eps keyword to fix this.\n");
      exit(1);
    }
    else if (ljgiven == 0x10) {
      printf("ReadEPRuleFile >> Error.  Lennard-Jones Epsilon was specified, "
	     "but no Sigma\nReadEPRuleFile >> was given.  Use the Sigma / sig "
	     "keyword to fix this.\n");
      exit(1);
    }	
    
    // Set the frame style based on keywords
    if (tmpstyle[0] != '\0') {
      if (strcmp(tmpstyle, "Flex-2") == 0 ||
          strcmp(tmpstyle, "FlexDis2") == 0) {
	tmr->frstyle = 1;
      }
      else if (strcmp(tmpstyle, "FD-2") == 0 ||
               strcmp(tmpstyle, "FixedDis2") == 0) {
	tmr->frstyle = 11;
      }
      else if (strcmp(tmpstyle, "Flex-3") == 0 ||
               strcmp(tmpstyle, "FlexDis3") == 0) {
	tmr->frstyle = 2;
      }
      else if (strcmp(tmpstyle, "FD-3") == 0 ||
               strcmp(tmpstyle, "FixedDis3") == 0) {
	tmr->frstyle = 3;
      }
      else if (strcmp(tmpstyle, "FAD-3") == 0 ||
               strcmp(tmpstyle, "FixAnglDis") == 0) {
	tmr->frstyle = 4;
      }
      else if (strcmp(tmpstyle, "Out-3") == 0 ||
               strcmp(tmpstyle, "OutOfPlane") == 0) {
	tmr->frstyle = 5;
      }
      else if (strcmp(tmpstyle, "FD-4") == 0 ||
               strcmp(tmpstyle, "FourPoint") == 0) {
	tmr->frstyle = 6;
      }
      else if (strcmp(tmpstyle, "None") == 0 ||
	       strcmp(tmpstyle, "NONE") == 0 ||
	       strcmp(tmpstyle, "ChangeQ") == 0) {
	tmr->frstyle = 0;
      }
      else {
	printf("ReadEPRuleFile >> Acceptable values for the frame include:\n\n"
	       "  NONE         There is no frame--rather, this rule describes "
	       "a real atom\n"
	       "  ChangeQ      whose partial charge is to change, probably to "
	       "accommodate\n"
	       "               extra points added elsewhere.\n\n"
	       "  Flex-2       Frame consists of atoms parent and frame2, "
	       "extra point along\n"
	       "  FlexDis2     the line beween them, at a displacement of v12 "
	       "times the\n"
	       "               vector parent->frame2, starting at parent.  A "
	       "value of 0.5\n"
	       "               will place an extra point at the bond "
	       "midpoint.  Negative\n"
	       "               values of v12 can model sigma holes halogen "
	       "parent atoms.\n\n"
	       "  FD-2         Frame consists of atoms parent and frame2, "
	       "extra point along\n"
	       "  FixedDis2    the line beween them, at an exact distance v1e "
	       "from parent.\n"
	       "               Negative values of v1e can again model sigma "
	       "holes for parent\n"
	       "               halogen atoms, and fixed dipoles can be "
	       "modeled with any\n"
	       "               value of v1e (use a non-neglible value, "
	       "perhaps >= 0.25, for\n"
	       "               best results on a normal PME grid).\n\n"
	       "  Flex-3       Frame consists of atoms parent, frame2, and "
	       "frame3, extra\n"
	       "  FlexDis3     point in the plane of the frame, at a "
	       "displacement of v12\n"
	       "               times the distance between parent and frame2 "
	       "plus v13 times\n"
	       "               the distance between parent and frame3.  Small "
	       "values of v12\n"
	       "               and v13 can create water molecules in the "
	       "TIP-4P family where\n"
	       "               oxygen is the parent atom, and will not "
	       "stretch if the frame\n"
	       "               itself is rigid.\n\n"
	       "  FD-3         Frame consists of atoms parent, frame2, and "
	       "frame3, extra point\n"
	       "  FixedDis3    in the plane of the frame, at an exact "
	       "distance v1e from the\n"
	       "               parent atom such that the extra point lies "
	       "along a line defined\n"
	       "               by the parent atom and a point Y.  The point Y "
	       "lies one the\n"
	       "               line between frame2 and frame3, at a "
	       "proportion v23 of the way\n"
	       "               from frame2 to frame3.  This is a planar "
	       "analog of frame FD-2.\n\n"
	       "  FAD-3        Frame consists of atoms parent, frame2, and "
	       "frame3, extra point\n"
	       "  FixAnglDis   in the plane of the frame, at an exact "
	       "distance v1e from the\n"
	       "               parent atom such that the extra point, parent, "
	       "and frame2 atoms\n"
	       "               form an angle theta (specified in radians).  "
	       "This is perhaps\n"
	       "               the most rigid planar frame, as both v1e and "
	       "theta are\n"
	       "               invariant to flexion of the frame itself.\n\n"
	       "  Out-3        Frame consists of atoms parent, frame2, and "
	       "frame3, extra point\n"
	       "  OutOfPlane   out of the plane of the frame.  Proportions of "
	       "displacements\n"
	       "               between the parent and frame2 or frame3 atoms, "
	       "v12 and v13,\n"
	       "               define a point Y in the plane.  Starting at Y, "
	       "the extra point\n"
	       "               is located on a plane normal at a proportion "
	       "v12x13 of the\n"
	       "               cross product of the two displacements (not "
	       "taking into account\n"
	       "               the proportions v12 and v13).\n\n"
	       "  FD-4         Frame consists of atoms parent, frame2, "
	       "frame3, and frame4,\n"
	       "  FourPoint    extra point placed atop a pyramidal structure "
	       "determined by the\n"
	       "               cross products of displacements between "
	       "frame2, frame3, and\n"
	       "               frame4.\n\n");
	       
	exit(1);
      }
    }
    
    // Frame style 0 is a special case--it means that this rule
    // merely modifies the nonbonded interactions of this atom.
    if (tmr->frstyle == 0) {
      tmr->nfratm = 0;

      // Check the charge.  The (impossible) default value of    
      // 1.1e30 means that electrostatics should not be modified.
      if (tmr->Charge < 1.0e30) {
        tmr->modElec = 1;
      }

      // Check the Lennard-Jones sigma and epsilon.  
      // Then, declare these parameters to be a new  
      // Lennard-Jones type, and expand the topology 
      // parameter lists in order to accommodate the 
      // new type.                                   
      if (tmr->eps >= 1.0e-8 && tmr->sig >= 1.0e-8) {
        tmr->LJIdx = tp->ntypes;
        tmr->modLJ = 1;
        ExpandLJTables(tp, tmr->sig, tmr->eps);
      }
    }

    // Check for duplicate specifiers
    if (epgiven == 1 && atomgiven == 1) {
      printf("ReadEPRuleFile >> Warning.  Both the ExtraPoint / epname and "
	     "the AtomName /\nReadEPRuleFile >> atom keywords were specified."
	     "These keywords point to the\nReadEPRuleFile >> same directive "
	     "and are intended for use in different\nReadEPRuleFile >> "
	     "contexts.\n");
    }
    
    // Check frame dimension specifications
    if (tmr->frstyle == 0) {
      IgnoredEPKeyword(v12given, "Vector12", "v12", "NONE");
      IgnoredEPKeyword(v13given, "Vector13", "v13", "NONE");
      IgnoredEPKeyword(v23given, "Vector23", "v23", "NONE");
      IgnoredEPKeyword(v1egiven, "Vector1E", "v1e", "NONE");
      IgnoredEPKeyword(thetagiven, "Theta", "theta", "NONE");
      IgnoredEPKeyword(v12x13given, "Vector13x23", "v13x23", "NONE");
      IgnoredEPKeyword(parentgiven, "ParentAtom", "parent", "NONE");
      IgnoredEPKeyword(frame2given, "FrameAtom2", "frame2", "NONE");
      IgnoredEPKeyword(frame3given, "FrameAtom3", "frame3", "NONE");
      IgnoredEPKeyword(frame4given, "FrameAtom4", "frame4", "NONE");
      if (epgiven == 1 && atomgiven == 0) {
	printf("ReadEPRuleFile >> This rule attempts to change the non-bonded "
	       "properties of\nReadEPRuleFile >> an atom, but invokes the "
               "ExtraPoint / epname keyword.  This\nReadEPRuleFile >> will "
	       "work, but the AtomName / atom keyword should be used\n"
	       "ReadEPRuleFile >> instead.\n");
      }      
    }
    else if (tmr->frstyle == 1) {
      MissingEPKeyword(v12given, v1egiven, "Vector12", "v12", "Flex-2");
      MissingEPKeyword(parentgiven, 0, "ParentAtom", "parent", "Flex-2");
      MissingEPKeyword(frame2given, 0, "FrameAtom2", "frame2", "Flex-2");
      MistakenEPKeyword(v12given, v1egiven, "Vector12", "v12", "Vector1E",
			"v1e", "Flex-2");
      IgnoredEPKeyword(v13given, "Vector13", "v13", "Flex-2");
      IgnoredEPKeyword(v23given, "Vector23", "v23", "Flex-2");
      IgnoredEPKeyword(v1egiven, "Vector1E", "v1e", "Flex-2");
      IgnoredEPKeyword(thetagiven, "Theta", "theta", "Flex-2");
      IgnoredEPKeyword(v12x13given, "Vector13x23", "v13x23", "Flex-2");
      IgnoredEPKeyword(frame3given, "FrameAtom3", "frame3", "Flex-2");
      IgnoredEPKeyword(frame4given, "FrameAtom4", "frame4", "Flex-2");
    }
    else if (tmr->frstyle == 11) {
      MissingEPKeyword(v1egiven, v12given, "Vector1E", "v1e", "FD-2");
      MissingEPKeyword(parentgiven, 0, "ParentAtom", "parent", "FD-2");
      MissingEPKeyword(frame2given, 0, "FrameAtom2", "frame2", "FD-2");
      MistakenEPKeyword(v1egiven, v12given, "Vector1E", "v1e", "Vector12",
			"v12", "FD-2");
      IgnoredEPKeyword(v13given, "Vector13", "v13", "FD-2");
      IgnoredEPKeyword(v23given, "Vector23", "v23", "FD-2");
      IgnoredEPKeyword(thetagiven, "Theta", "theta", "FD-2");
      IgnoredEPKeyword(v12x13given, "Vector13x23", "v13x23", "FD-2");
      IgnoredEPKeyword(frame3given, "FrameAtom3", "frame3", "FD-2");
      IgnoredEPKeyword(frame4given, "FrameAtom4", "frame4", "FD-2");
    }
    else if (tmr->frstyle == 2) {
      MissingEPKeyword(v12given, v1egiven, "Vector12", "v12", "Flex-3");
      MissingEPKeyword(v13given, thetagiven, "Vector13", "v13", "Flex-3");
      MissingEPKeyword(parentgiven, 0, "ParentAtom", "parent", "Flex-3");
      MissingEPKeyword(frame2given, 0, "FrameAtom2", "frame2", "Flex-3");
      MissingEPKeyword(frame3given, 0, "FrameAtom3", "frame3", "Flex-3");
      MistakenEPKeyword(v12given, v1egiven, "Vector12", "v12", "Vector1E",
			"v1e", "Flex-3");
      MistakenEPKeyword(v13given, thetagiven, "Vector13", "v13", "Theta",
			"theta", "Flex-3");
      IgnoredEPKeyword(v23given, "Vector23", "v23", "Flex-3");
      IgnoredEPKeyword(v1egiven, "Vector1E", "v1e", "Flex-3");
      IgnoredEPKeyword(v12x13given, "Vector13x23", "v13x23", "Flex-3");
      IgnoredEPKeyword(frame4given, "FrameAtom4", "frame4", "Flex-3");
    }
    else if (tmr->frstyle == 3) {
      MissingEPKeyword(v1egiven, v12given, "Vector1E", "v1e", "FD-3");
      MissingEPKeyword(v23given, v12x13given, "Vector23", "v23", "FD-3");
      MissingEPKeyword(parentgiven, 0, "ParentAtom", "parent", "FD-3");
      MissingEPKeyword(frame2given, 0, "FrameAtom2", "frame2", "FD-3");
      MissingEPKeyword(frame3given, 0, "FrameAtom3", "frame3", "FD-3");
      MistakenEPKeyword(v1egiven, v12given, "Vector1E", "v1e", "Vector12",
			"v12", "FD-3");
      MistakenEPKeyword(v23given, v12x13given, "Vector23", "v23",
			"Vector12x13", "v12x13", "FD-3");
      IgnoredEPKeyword(v13given, "Vector13", "v13", "FD-3");
      IgnoredEPKeyword(thetagiven, "Theta", "theta", "FD-3");
      IgnoredEPKeyword(frame4given, "FrameAtom4", "frame4", "FD-3");
    }
    else if (tmr->frstyle == 4) {
      MissingEPKeyword(v1egiven, v12given, "Vector1E", "v1e", "FAD-3");
      MissingEPKeyword(thetagiven, v13given, "Theta", "theta", "FAD-3");
      MissingEPKeyword(parentgiven, 0, "ParentAtom", "parent", "FAD-3");
      MissingEPKeyword(frame2given, 0, "FrameAtom2", "frame2", "FAD-3");
      MissingEPKeyword(frame3given, 0, "FrameAtom3", "frame3", "FAD-3");
      MistakenEPKeyword(v1egiven, v12given, "Vector1E", "v1e", "Vector12",
                        "v12", "FAD-3");
      MistakenEPKeyword(thetagiven, v13given, "Theta", "theta", "Vector13",
                        "v13", "FAD-3");
      IgnoredEPKeyword(v12x13given, "Vector12x13", "v12x13", "FAD-3");
      IgnoredEPKeyword(v23given, "Vector23", "v23", "FAD-3");
      IgnoredEPKeyword(frame4given, "FrameAtom4", "frame4", "FAD-3");
    }
    else if (tmr->frstyle == 5) {
      MissingEPKeyword(v12given, v1egiven, "Vector12", "v12", "Out-3");
      MissingEPKeyword(v13given, thetagiven, "Vector13", "v13", "Out-3");
      MissingEPKeyword(v12x13given, v23given, "Vector12x13", "v12x13",
                       "Out-3");
      MissingEPKeyword(parentgiven, 0, "ParentAtom", "parent", "Out-3");
      MissingEPKeyword(frame2given, 0, "FrameAtom2", "frame2", "Out-3");
      MissingEPKeyword(frame3given, 0, "FrameAtom3", "frame3", "Out-3");
      MistakenEPKeyword(v12given, v1egiven, "Vector12", "v12", "Vector1E",
                        "v1e", "Out-3");
      MistakenEPKeyword(v13given, thetagiven, "Vector13", "v13", "Theta",
                        "theta", "Out-3");
      MistakenEPKeyword(v12x13given, v23given, "Vector12x13", "v12x13",
			"Vector23", "v23", "Out-3");
      IgnoredEPKeyword(frame4given, "FrameAtom4", "frame4", "Out-3");
    }
    else if (tmr->frstyle == 6) {
      MissingEPKeyword(v1egiven, v12given, "Vector1E", "v1e", "FD-4");
      MissingEPKeyword(v13given, thetagiven, "Vector13", "v13", "FD-4");
      MissingEPKeyword(v23given, v12x13given, "Vector23", "v23", "FD-4");
      MissingEPKeyword(parentgiven, 0, "ParentAtom", "parent", "FD-4");
      MissingEPKeyword(frame2given, 0, "FrameAtom2", "frame2", "FD-4");
      MissingEPKeyword(frame3given, 0, "FrameAtom3", "frame3", "FD-4");
      MissingEPKeyword(frame4given, 0, "FrameAtom4", "frame4", "FD-4");
      MistakenEPKeyword(v1egiven, v12given, "Vector1E", "v1e", "Vector12",
                        "v12", "FD-4");
      MistakenEPKeyword(v13given, thetagiven, "Vector13", "v13", "Theta",
                        "theta", "FD-4");
      MistakenEPKeyword(v23given, v12x13given, "Vector23", "v23",
			"Vector12x13", "v12x13", "FD-4");

      // The ordering of parameters in frame style 6 needs adjustment
      // to fit mdgx's internal conventions.  From input, we've taken
      // d1 as the v1e parameter, d2 as the v13, and d3 as the v23.
      // What we need is for d1 to be what the user entered for v13,
      // d2 what the user entered for v23, and d3 what the user
      // entered for v1e.  Once these values are reordered, they can
      // be written to a topology correctly.  Likewise, reading these
      // values from the topology as opposed to the &rule namelist
      // will have them in the correct order.
      realtmp = tmr->d1;
      tmr->d1 = tmr->d2;
      tmr->d2 = tmr->d3;
      tmr->d3 = realtmp;
    }
    
    // Fix the impossible charge value if the user did not include one
    if (tmr->Charge > 1.0e30) {
      tmr->Charge = 0.0;
    }

    // Check the style of this rule.  Type 1 is a two-atom     
    // frame, all others save for type 6 are three-atom frames.
    if (tmr->frstyle < 0) {
      printf("ReadEPRuleFile >> Error.  No frame type was specified.  Use "
	     "keywords\nReadEPRuleFile >> FrameType / eptype to describe the "
	     "frame, or enter the\nReadEPRuleFile >> numerical identifier "
	     "with the FrameStyle / style\nReadEPRuleFile >> keywords (the "
	     "former input is recommended).\n");
      exit(1);
    }
    if (tmr->frstyle > 11) {
      printf("ReadEPRuleFile >> Error.  Frame style %d is not allowed.\n",
             tmr->frstyle);
      exit(1);
    }
    if (tmr->frstyle == 1 || tmr->frstyle == 11) {
      tmr->nfratm = 2;
    }
    else if (tmr->frstyle < 6) {
      tmr->nfratm = 3;
    }
    else if (tmr->frstyle == 6 || tmr->frstyle == 9 || tmr->frstyle == 10) {
      tmr->nfratm = 4;
    }
    else {
      printf("ReadEPRuleFile >> Error.  Frame style %d is invalid.\n",
             tmr->frstyle);
      exit(1);
    }

    // Check that the atom and residue names are properly specified.
    if (strlen(tmpfr.map[0]) > 4 || strlen(tmpfr.map[1]) > 4 ||
        strlen(tmpfr.map[2]) > 4 || strlen(tmpfr.map[3]) > 4 ||
        strlen(tmpepn) > 4) {
      sprintf(line, "%s %s %s %s %s", tmpfr.map[0], tmpfr.map[1], tmpfr.map[2],
              tmpfr.map[3], tmpepn);
      printf("ReadEPRuleFile >> Error.  Atom names [ %s ] invalid (4 "
             "characters max\nReadEPRuleFile >> per name).\n", line);
      exit(1);
    }
    if (strlen(tmpres) > 4) {
      printf("ReadEPRuleFile >> Error.  Residue name [ %s ] invalid (4 "
             "characters max).\n", tmpres);
      exit(1);
    }
    for (i = 0; i < 4; i++) {
      if (tmpfr.map[i][0] != '\0') {
        for (j = strlen(tmpfr.map[i]); j < 4; j++) {
          tmpfr.map[i][j] = ' ';
        }
      }
    }
    for (i = strlen(tmpres); i < 4; i++) {
      tmpres[i] = ' ';
    }
    for (i = strlen(tmpepn); i < 4; i++) {
      tmpepn[i] = ' ';
    }
    for (j = 0; j < 4; j++) {
      tmr->resname[j] = tmpres[j];
      tmr->epname[j] = tmpepn[j];
      if (tmr->frstyle > 0) {
        tmr->fr1[j] = tmpfr.map[0][j];
        tmr->fr2[j] = tmpfr.map[1][j];
      }
      if (tmr->frstyle > 1) {
        tmr->fr3[j] = tmpfr.map[2][j];
      }
      if (tmr->frstyle == 6 || tmr->frstyle == 9 || tmr->frstyle == 10) {
        tmr->fr4[j] = tmpfr.map[3][j];
      }
    }
    if ((tmr->frstyle == 1 && (strncmp(tmr->fr1, tmr->fr2, 4) == 0 ||
                               strncmp(tmr->epname, tmr->fr1, 4) == 0 ||
                               strncmp(tmr->epname, tmr->fr2, 4) == 0)) ||
        (tmr->frstyle > 1 && (strncmp(tmr->fr1, tmr->fr3, 4) == 0 ||
                              strncmp(tmr->fr2, tmr->fr3, 4) == 0 ||
                              strncmp(tmr->epname, tmr->fr3, 4) == 0)) ||
        ((tmr->frstyle == 6 || tmr->frstyle == 9 || tmr->frstyle == 10) &&
         (strncmp(tmr->fr2, tmr->fr4, 4) == 0 ||
          strncmp(tmr->fr2, tmr->fr4, 4) == 0 ||
          strncmp(tmr->epname, tmr->fr4, 4) == 0))) {
      sprintf(line, "%s %s %s %s", tmpfr.map[0], tmpfr.map[1], tmpfr.map[2],
              tmpfr.map[3]);
      printf("ReadEPRuleFile >> Error.  One or more frame atoms is not "
             "properly specified.\nReadEPRuleFile >> Frame atoms are [ %s ]."
             "\n", line);
      exit(1);
    }

    // Increment the number of known rules
    neprule++;
    if (neprule == maxrule) {
      maxrule += 10;
      tp->eprules = (eprule*)realloc(tp->eprules, maxrule*sizeof(eprule));
    }
  }

  // Adjust the length of the extra point rule list to 
  // accommodate only as many rules as have been found.
  tp->eprules = (eprule*)realloc(tp->eprules, neprule*sizeof(eprule));
  tp->neprule = neprule;

  // Free allocated memory
  DestroyCmat(&tmpfr);
}

//-----------------------------------------------------------------------------
// CountBonds: count all bonds extending to and from atom ai according to
//             topology tp.                                          
//
// Arguments:
//   tp:     the system topology
//   ai:     the atom for which to count bonds
//-----------------------------------------------------------------------------
static int CountBonds(prmtop *tp, int ai)
{
  int i, j, nbnd;

  nbnd = 0;
  for (i = 0; i < tp->natom; i++) {
    if (i == ai) {
      nbnd += tp->BLC[i].nbond;
    }
    else {
      for (j = 0; j < tp->BLC[i].nbond; j++) {
        if (tp->BLC[i].BC[j].b == ai) {
          nbnd++;
        }
      }
    }
  }

  return nbnd;
}

//-----------------------------------------------------------------------------
// EPRuleDefined: this routine searches through a list of rules in the  
//                topology tp to see if any rules matching the residue name can
//                be found.  Although this function will return the number of a
//                matching rule if one is found, there may be more matching
//                rules, so the list is re-checked in the ImplementEPRule
//                routine.                       
//
// Arguments:
//   tp:      system topology
//   resid:   index of the residue to check
//   atmid:   index of the atom within the residue to check
//-----------------------------------------------------------------------------
static int EPRuleDefined(prmtop *tp, int resid, int atmid)
{
  int i, epid;

  epid = -1;
  for (i = 0; i < tp->neprule; i++) {
    if (strncmp(tp->eprules[i].resname, &tp->ResNames[4*resid], 4) == 0 &&
        strncmp(tp->eprules[i].epname, &tp->AtomNames[4*atmid], 4) == 0) {
      epid = i;
      break;
    }
  }

  return epid;
}

//-----------------------------------------------------------------------------
// ReallocateXtraPtsArray: reallocate a topology's extra point array in the
//                         event that more space is or will be needed.
//
// Arguments:
//   tp:      the topology at hand
//   safety:  the amount of space that must be available (or else array length
//            will be doubled)
//-----------------------------------------------------------------------------
static void ReallocateXtraPtsArray(prmtop *tp, int safety)
{
  int i;

  if (tp->nxtrapt == 0 || tp->nxtrapt > tp->maxtrapt - safety) {
    tp->maxtrapt = (tp->nxtrapt < 32) ? tp->maxtrapt + 32 : tp->maxtrapt * 2;
    tp->xtrapts = (expt*)realloc(tp->xtrapts, tp->maxtrapt*sizeof(expt));
    for (i = tp->nxtrapt; i < tp->maxtrapt; i++) {
      tp->xtrapts[i].atomid = -1;
      tp->xtrapts[i].frstyle = -1;
      tp->xtrapts[i].nfratm = 0;
      tp->xtrapts[i].fr1 = -1;
      tp->xtrapts[i].fr2 = -1;
      tp->xtrapts[i].fr3 = -1;
      tp->xtrapts[i].fr4 = -1;
      tp->xtrapts[i].atm11 = 0;
      tp->xtrapts[i].d1 = 0.0;
      tp->xtrapts[i].d2 = 0.0;
      tp->xtrapts[i].d3 = 0.0;
    }
  }  
}

//-----------------------------------------------------------------------------
// ImplementEPRule: implement a user-defined extra point rule.          
//                                                                      
// Arguments:                                                           
//   tp:     the topology structure                                     
//   resid:  the residue number                                         
//   epid:   the extra point rule number                                
//-----------------------------------------------------------------------------
static void ImplementEPRule(prmtop *tp, int resid, int epid)
{
  int i, llim, hlim;
  eprule *ert;

  llim = tp->ResLims[resid];
  hlim = tp->ResLims[resid+1];

  // Transfer information from the rule to the particle
  ert = &tp->eprules[epid];
  tp->xtrapts[tp->nxtrapt].frstyle = ert->frstyle;
  tp->xtrapts[tp->nxtrapt].nfratm = ert->nfratm;
  tp->xtrapts[tp->nxtrapt].d1 = ert->d1;
  tp->xtrapts[tp->nxtrapt].d2 = ert->d2;
  tp->xtrapts[tp->nxtrapt].d3 = ert->d3;

  // Find the extra point and frame atoms
  for (i = llim; i < hlim; i++) {
    if (strncmp(ert->epname, &tp->AtomNames[4*i], 4) == 0) {
      tp->xtrapts[tp->nxtrapt].atomid = i;
    }
    if (strncmp(ert->fr1, &tp->AtomNames[4*i], 4) == 0) {
      tp->xtrapts[tp->nxtrapt].fr1 = i;
    }
    if (strncmp(ert->fr2, &tp->AtomNames[4*i], 4) == 0) {
      tp->xtrapts[tp->nxtrapt].fr2 = i;
    }
    if (ert->nfratm > 2) {
      if (strncmp(ert->fr3, &tp->AtomNames[4*i], 4) == 0) {
        tp->xtrapts[tp->nxtrapt].fr3 = i;
      }
    }
    if (ert->nfratm > 3) {
      if (strncmp(ert->fr4, &tp->AtomNames[4*i], 4) == 0) {
        tp->xtrapts[tp->nxtrapt].fr4 = i;
      }
    }
  }

  // Increment the number of extra points
  tp->nxtrapt += 1;
  ReallocateXtraPtsArray(tp, 1);
}

//-----------------------------------------------------------------------------
// AmberFrame2mdgx: convert the frame information for a custom virtual site
//                  in an Amber prmtop into mdgx's internal conventions.  Most
//                  of the custom virtual sites in the rest of Amber were in
//                  mdgx long beforehand, but a couple of new frame types
//                  were added and the numbering scheme used in pmemd cannot
//                  match that used in mdgx.  This routine encapsulates the
//                  conversion.
//
// Arguments:
//   tp:     the topology structure
//   frmID:  the topology frame to convert
//-----------------------------------------------------------------------------
static void AmberFrame2mdgx(prmtop *tp, int frmID)
{
  int i, amberStyle;
  expt *tpxtra;
  
  amberStyle = tp->customVSi[6*frmID + 5];
  tpxtra = &tp->xtrapts[tp->nxtrapt];
  if (amberStyle == 4) {

    // Amber style 4 is mdgx style 1, in-line flexible two atom frame
    tpxtra->frstyle = 1;
    tpxtra->nfratm  = 2;
    tpxtra->atomid  = tp->customVSi[6*frmID] - 1;
    tpxtra->fr1     = tp->customVSi[6*frmID + 1] - 1;
    tpxtra->fr2     = tp->customVSi[6*frmID + 2] - 1;
    tpxtra->d1      = tp->customVSd[3*frmID];
  }
  else if (amberStyle == 5) {

    // Amber style 5 is mdgx style 2, in-plane flexible three atom frame
    tpxtra->frstyle = 2;
    tpxtra->nfratm  = 3;
    tpxtra->atomid  = tp->customVSi[6*frmID] - 1;
    tpxtra->fr1     = tp->customVSi[6*frmID + 1] - 1;
    tpxtra->fr2     = tp->customVSi[6*frmID + 2] - 1;
    tpxtra->fr3     = tp->customVSi[6*frmID + 3] - 1;
    tpxtra->d1      = tp->customVSd[3*frmID];
    tpxtra->d2      = tp->customVSd[3*frmID + 1];
  }
  else if (amberStyle == 6) {

    // Amber style 6 is mdgx style 3, in-plane fixed-distance three atom frame
    tpxtra->frstyle = 3;
    tpxtra->nfratm  = 3;
    tpxtra->atomid  = tp->customVSi[6*frmID] - 1;
    tpxtra->fr1     = tp->customVSi[6*frmID + 1] - 1;
    tpxtra->fr2     = tp->customVSi[6*frmID + 2] - 1;
    tpxtra->fr3     = tp->customVSi[6*frmID + 3] - 1;
    tpxtra->d1      = tp->customVSd[3*frmID];
    tpxtra->d3      = tp->customVSd[3*frmID + 1];
  }
  else if (amberStyle == 7) {

    // Amber style 6 is mdgx style 4, in-plane fixed angle and distance
    // three atom frame
    tpxtra->frstyle = 4;
    tpxtra->nfratm  = 3;
    tpxtra->atomid  = tp->customVSi[6*frmID] - 1;
    tpxtra->fr1     = tp->customVSi[6*frmID + 1] - 1;
    tpxtra->fr2     = tp->customVSi[6*frmID + 2] - 1;
    tpxtra->fr3     = tp->customVSi[6*frmID + 3] - 1;
    tpxtra->d1      = tp->customVSd[3*frmID];
    tpxtra->d2      = tp->customVSd[3*frmID + 1];
  }
  else if (amberStyle == 8) {

    // Amber style 8 is mdgx style 5, out-of-plane flexible three atom frame
    tpxtra->frstyle = 5;
    tpxtra->nfratm  = 3;
    tpxtra->atomid  = tp->customVSi[6*frmID] - 1;
    tpxtra->fr1     = tp->customVSi[6*frmID + 1] - 1;
    tpxtra->fr2     = tp->customVSi[6*frmID + 2] - 1;
    tpxtra->fr3     = tp->customVSi[6*frmID + 3] - 1;
    tpxtra->d1      = tp->customVSd[3*frmID];
    tpxtra->d2      = tp->customVSd[3*frmID + 1];
    tpxtra->d3      = tp->customVSd[3*frmID + 2];
  }
  else if (amberStyle == 9) {

    // Amber style 9 is mdgx style 11, in-line fixed distance two atom frame
    tpxtra->frstyle = 11;
    tpxtra->nfratm  = 2;
    tpxtra->atomid  = tp->customVSi[6*frmID] - 1;
    tpxtra->fr1     = tp->customVSi[6*frmID + 1] - 1;
    tpxtra->fr2     = tp->customVSi[6*frmID + 2] - 1;
    tpxtra->d1      = tp->customVSd[3*frmID];
  }
  else if (amberStyle == 10) {

    // Amber style 10 is mdgx style 6, the four atom frame
    tpxtra->frstyle = 6;
    tpxtra->nfratm  = 4;
    tpxtra->atomid  = tp->customVSi[6*frmID] - 1;
    tpxtra->fr1     = tp->customVSi[6*frmID + 1] - 1;
    tpxtra->fr2     = tp->customVSi[6*frmID + 2] - 1;
    tpxtra->fr3     = tp->customVSi[6*frmID + 3] - 1;
    tpxtra->fr4     = tp->customVSi[6*frmID + 4] - 1;
    tpxtra->d1      = tp->customVSd[3*frmID];
    tpxtra->d2      = tp->customVSd[3*frmID + 1];
    tpxtra->d3      = tp->customVSd[3*frmID + 2];
  }
  
  // Increment the number of extra points
  tp->nxtrapt += 1;
  ReallocateXtraPtsArray(tp, 1);
}

//-----------------------------------------------------------------------------
// FindBondLength: find the length of a bond connecting two atoms A and B,
//                 first searching bonds controlled by atom A and then bonds
//                 controlled by atom B.                     
//                                                                      
// Arguments:                                                           
//   tp:     the topology structure                                     
//   A:      the first atom                                             
//   B:      the first atom                                             
//-----------------------------------------------------------------------------
static double FindBondLength(prmtop *tp, int A, int B)
{
  int i;

  // Check bonds controlled by atom A for recipient atom B
  for (i = 0; i < tp->BLC[A].nbond; i++) {
    if (tp->BLC[A].BC[i].b == B) {
      return tp->BParam[tp->BLC[A].BC[i].t].l0;
    }
  }

  // Check bonds controlled by atom B for recipient atom A
  for (i = 0; i < tp->BLC[B].nbond; i++) {
    if (tp->BLC[B].BC[i].b == A) {
      return tp->BParam[tp->BLC[B].BC[i].t].l0;
    }
  }

  // If we're still here, there was a problem, no bond was found.
  printf("FindBondLength >> Error.  No bond could be located between "
         "the extra\nFindBondLength >> point %.4s (atom %d) and the parent "
         "atom %.4s (atom %d).\n", &tp->AtomNames[4*A], A + 1,
	 &tp->AtomNames[4*B], B + 1);
  exit(1);
}

//-----------------------------------------------------------------------------
// Implement4PWaterRule: implement a rule for the family of four-point rigid
//                       waters.                                  
//                                                                      
// Arguments:                                                           
//   tp:     the topology structure                                     
//   resid:  the residue number                                         
//-----------------------------------------------------------------------------
static void Implement4PWaterRule(prmtop *tp, int resid)
{
  int i, nxt, llim, hlim;
  double maxmass, repo, rho, rhh, rise;

  // Before continuing in this function, check to make
  // sure that we won't overrun the extra points array
  ReallocateXtraPtsArray(tp, 1);
  
  // The frame style is "2"
  nxt = tp->nxtrapt;
  tp->xtrapts[nxt].frstyle = 2;
  tp->xtrapts[nxt].nfratm = 3;
  tp->xtrapts[nxt].fr2 = -1;
  tp->xtrapts[nxt].fr3 = -1;

  // Find the extra point and frame atoms by mass.
  // The most massive atom must be the oxygen, and
  // the extra point must have negligible mass of 
  // less than 1.0e-8.                            
  llim = tp->ResLims[resid];
  hlim = tp->ResLims[resid+1];
  maxmass = 0.0;
  for (i = llim; i < hlim; i++) {
    if (tp->Masses[i] < 1.0e-8) {
      tp->xtrapts[nxt].atomid = i;
    }
    if (tp->Masses[i] > maxmass) {
      maxmass = tp->Masses[i];
      tp->xtrapts[nxt].fr1 = i;
    }
  }
  for (i = llim; i < hlim; i++) {
    if (i != tp->xtrapts[nxt].fr1 && i != tp->xtrapts[nxt].atomid) {
      if (tp->xtrapts[nxt].fr2 == -1) {
        tp->xtrapts[nxt].fr2 = i;
      }
      else {
        tp->xtrapts[nxt].fr3 = i;
      }
    }
  }

  // Find the distance form the extra point to the oxygen
  repo = FindBondLength(tp, tp->xtrapts[nxt].atomid, tp->xtrapts[nxt].fr1);
  rho = FindBondLength(tp, tp->xtrapts[nxt].fr1, tp->xtrapts[nxt].fr2);
  rhh = FindBondLength(tp, tp->xtrapts[nxt].fr2, tp->xtrapts[nxt].fr3);
  rise = 2.0 * sqrt(rho*rho - 0.25*rhh*rhh);
  tp->xtrapts[nxt].d1 = repo / rise;
  tp->xtrapts[nxt].d2 = repo / rise;

  // The extra point is 1:1 connected only to frame atom 1
  tp->xtrapts[nxt].atm11 = 1;

  // Increment the extra point counter
  tp->nxtrapt += 1;
  ReallocateXtraPtsArray(tp, 1);
}

//-----------------------------------------------------------------------------
// Implement5PWaterRule: implement a rule for the family of five-point rigid
//                       waters.                                  
//                                                                      
// Arguments:                                                           
//   tp:     the topology structure                                     
//   resid:  the residue number                                         
//-----------------------------------------------------------------------------
static void Implement5PWaterRule(prmtop *tp, int resid)
{
  int i, nepfound, nxt, llim, hlim;
  double maxmass, roep, rhh, roh, ahoh;
  double bunnyrun, bunnyrise, hydrorun, cpmag, dfac;

  // Before continuing in this function, check to make
  // sure that we won't overrun the extra points array
  ReallocateXtraPtsArray(tp, 2);

  // The frame style is "5"
  nxt = tp->nxtrapt;
  tp->xtrapts[nxt].frstyle = 5;
  tp->xtrapts[nxt].nfratm = 3;
  tp->xtrapts[nxt].fr2 = -1;
  tp->xtrapts[nxt].fr3 = -1;
  tp->xtrapts[nxt+1].frstyle = 5;
  tp->xtrapts[nxt+1].nfratm = 3;
  tp->xtrapts[nxt+1].fr2 = -1;
  tp->xtrapts[nxt+1].fr3 = -1;

  // Find the extra points and frame atoms by mass.
  // The most massive atom must be the oxygen, and 
  // the extra point must have negligible mass of  
  // less than 1.0e-8.                             
  llim = tp->ResLims[resid];
  hlim = tp->ResLims[resid+1];
  maxmass = 0.0;
  nepfound = 0;
  for (i = llim; i < hlim; i++) {
    if (tp->Masses[i] < 1.0e-8) {
      tp->xtrapts[nxt+nepfound].atomid = i;
      nepfound++;
    }
    if (tp->Masses[i] > maxmass) {
      maxmass = tp->Masses[i];
      tp->xtrapts[nxt].fr1 = i;
      tp->xtrapts[nxt+1].fr1 = i;
    }
  }
  for (i = llim; i < hlim; i++) {
    if (i != tp->xtrapts[nxt].fr1 &&
        i != tp->xtrapts[nxt].atomid &&
        i != tp->xtrapts[nxt+1].atomid) {
      if (tp->xtrapts[nxt].fr2 == -1) {
        tp->xtrapts[nxt].fr2 = i;
        tp->xtrapts[nxt+1].fr2 = i;
      }
      else {
        tp->xtrapts[nxt].fr3 = i;
        tp->xtrapts[nxt+1].fr3 = i;
      }
    }
  }

  // Find the distances from the extra points to the oxygen
  roep = FindBondLength(tp, tp->xtrapts[nxt].atomid, tp->xtrapts[nxt].fr1);
  ahoh = FindBondLength(tp, tp->xtrapts[nxt+1].atomid, tp->xtrapts[nxt+1].fr1);
  if (fabs(roep - ahoh) > 1.0e-8) {
    printf("Implement5PWaterRule >> Error.  Extra point distances to "
           "oxygen atom differ \nImplement5PWaterRule >> by %13.10lf\n",
           fabs(roep - ahoh));
    exit(1);
  }
  roh = FindBondLength(tp, tp->xtrapts[nxt].fr1, tp->xtrapts[nxt].fr2);
  rhh = FindBondLength(tp, tp->xtrapts[nxt].fr2, tp->xtrapts[nxt].fr3);

  // The angle between the extra point, oxygen, and
  // hydrogen is hard-wired to be 109.47122 degrees   
  bunnyrun  = roep*cos(0.5*TETANG);
  bunnyrise = roep*sin(0.5*TETANG);
  hydrorun = sqrt(roh*roh - 0.25*rhh*rhh);
  ahoh = PI - 2.0*acos((0.5*rhh) / roh);
  dfac = -0.5*(bunnyrun/hydrorun);
  tp->xtrapts[nxt].d1 = dfac;
  tp->xtrapts[nxt].d2 = dfac;
  tp->xtrapts[nxt+1].d1 = dfac;
  tp->xtrapts[nxt+1].d2 = dfac;

  // Now the cross-product of the O-H1 and O-H2 bonds
  // must be computed in order to determine d3       
  cpmag = roh*roh*sin(ahoh);
  tp->xtrapts[nxt].d3 = bunnyrise/cpmag;
  tp->xtrapts[nxt+1].d3 = -bunnyrise/cpmag;

  // The extra points are 1:1 connected only to frame atom 1
  tp->xtrapts[nxt].atm11 = 1;
  tp->xtrapts[nxt+1].atm11 = 1;

  // Increment the extra point counter.  No further  
  // expansion of the extra points array is needed as
  // that was done at the beginning of the function. 
  tp->nxtrapt += 2;
  ReallocateXtraPtsArray(tp, 2);
}

//-----------------------------------------------------------------------------
// ImplementType1EP: implement a type I extra point.  This means planar
//                   geometry if there is only one extra point bonded to the
//                   parent atom, tetrahedral if there are two (there's no
//                   other number of points that we know how to do).
//
// Arguments:
//   tp:     the topology structure
//   epID:   the index of the extra point in the topology
//   pos:    the position of the extra point in the list for the same parent
//           atom
//   nep:    the number of extra points bound to the same parent atom
//-----------------------------------------------------------------------------
static void ImplementType1EP(prmtop *tp, int epID, int pos, int nep)
{
  int i, j, parAtom, nxt;
  double dL;
  
  // Before continuing in this function, check to make
  // sure that we won't overrun the extra points array
  ReallocateXtraPtsArray(tp, 1);
  parAtom = tp->nb1234[epID].L11[0];
  nxt = tp->nxtrapt;
  tp->xtrapts[nxt].atomid = epID;
  tp->xtrapts[nxt].frstyle = (nep == 1) ? 7 : 8;
  tp->xtrapts[nxt].nfratm = 3;
  tp->xtrapts[nxt].fr1 = parAtom;
  j = 0;
  for (i = 0; i < tp->nb1234[tp->xtrapts[nxt].fr1].n12; i++) {
    if (tp->nb1234[tp->xtrapts[nxt].fr1].L12[i] != epID) {
      if (j == 0) {
        tp->xtrapts[nxt].fr2 = tp->nb1234[tp->xtrapts[nxt].fr1].L12[i];
      }
      else {
        tp->xtrapts[nxt].fr3 = tp->nb1234[tp->xtrapts[nxt].fr1].L12[i];
      }
      j++;
    }
  }
  tp->xtrapts[nxt].atm11 = 1;

  // Handle the idiosyncracies of one or two EP
  if (nep == 1) {
    tp->xtrapts[nxt].d1 = FindBondLength(tp, epID, tp->xtrapts[nxt].fr1);
  }
  else {
    dL = FindBondLength(tp, epID, tp->xtrapts[nxt].fr1);
    tp->xtrapts[nxt].d1 = dL * cos(0.5*TETANG);
    tp->xtrapts[nxt].d2 = dL * sin(0.5*TETANG);
    if (pos == 0) {
      tp->xtrapts[nxt].d2 *= -1.0;
    }
  }

  // Increment the extra point count
  tp->nxtrapt += 1;
  ReallocateXtraPtsArray(tp, 1);
}

//-----------------------------------------------------------------------------
// ImplementType2EP: implement a type II extra point.  This means planar
//                   geometry with the points in the plane of the 
//                   parent atom.  The geometry is linear if there's only one
//                   atom (i.e. a sigma hole on a halogen), or trigonal planar
//                   if there are two extra points.
//
// Arguments:
//   tp:     the topology structure
//   epID:   the index of the extra point
//   pos:    the position of the extra point in the list for the same parent
//           atom
//   nep:    the number of extra points bound to the same parent atom
//-----------------------------------------------------------------------------
static void ImplementType2EP(prmtop *tp, int epID, int pos, int nep)
{
  int i, j, parAtom, secAtom, nxt;
  double dL;
  
  // Before continuing in this function, check to make
  // sure that we won't overrun the extra points array
  ReallocateXtraPtsArray(tp, 1);
  parAtom = tp->nb1234[epID].L11[0];
  nxt = tp->nxtrapt;
  tp->xtrapts[nxt].atomid = epID;
  tp->xtrapts[nxt].frstyle = (nep == 1) ? 9 : 10;
  tp->xtrapts[nxt].nfratm = 4;
  tp->xtrapts[nxt].fr1 = parAtom;
  for (i = 0; i < tp->nb1234[tp->xtrapts[nxt].fr1].n12; i++) {
    if (tp->nb1234[tp->xtrapts[nxt].fr1].L12[i] != epID) {
      tp->xtrapts[nxt].fr2 = tp->nb1234[tp->xtrapts[nxt].fr1].L12[i];
      break;
    }
  }
  j = 0;
  for (i = 0; i < tp->nb1234[tp->xtrapts[nxt].fr2].n12; i++) {
    if (tp->nb1234[tp->xtrapts[nxt].fr2].L12[i] != tp->xtrapts[nxt].fr1 &&
        tp->Masses[tp->nb1234[tp->xtrapts[nxt].fr2].L12[i]] > 1.0e-8) {
      if (j == 0) {
        tp->xtrapts[nxt].fr3 = tp->nb1234[tp->xtrapts[nxt].fr2].L12[i];
      }
      else {
        tp->xtrapts[nxt].fr4 = tp->nb1234[tp->xtrapts[nxt].fr2].L12[i];
      }
      j++;
    }
  }
  tp->xtrapts[nxt].atm11 = 1;
  if (nep == 1) {
    tp->xtrapts[nxt].d1 = FindBondLength(tp, epID, tp->xtrapts[nxt].fr1);
  }
  else {
    dL = FindBondLength(tp, epID, tp->xtrapts[nxt].fr1);
    tp->xtrapts[nxt].d1 = 0.5*dL;
    tp->xtrapts[nxt].d2 = 0.5*sqrt(3.0)*dL;
    if (pos == 1) {
      tp->xtrapts[nxt].d2 *= -1.0;
    }
  }
  
  // Increment the extra point count
  tp->nxtrapt += 1;
  ReallocateXtraPtsArray(tp, 1);
}

//-----------------------------------------------------------------------------
// ShivForInt: insert an integer into an array at a particular point, shifting
//             all other elements of the array one step down.  The original
//             array is reallocated (lengthened) and then returned.
//                                                                      
// Arguments:                                                           
//   V:     the original array of integers (also returned)              
//   n:     the original size of the array V                            
//   p:     the position at which to add the new integer (all elements of the
//          original array p, p+1, p+2, ..., n-1 will be shifted to positions
//          p+1, p+2, p+3, ..., n)                         
//   a:     the value of the integer to add                             
//-----------------------------------------------------------------------------
static int* ShivForInt(int* V, int n, int p, int a)
{
  int i;

  // Reallocate the array
  V = (int*)realloc(V, (n+1)*sizeof(int));

  // Shift existing elements
  for (i = n; i > p; i--) {
    V[i] = V[i-1];
  }
  V[p] = a;

  return V;
}

//-----------------------------------------------------------------------------
// ShivForDouble: insert a double-precision real number into an array at a
//                particular point, shifting all other elements of the array
//                one step down.  The original array is reallocated
//                (lengthened) and then returned.           
//                                                                      
// Arguments:                                                           
//   V:     the original array of doubles (also returned)               
//   n:     the original size of the array V                            
//   p:     the position at which to add the new double (all elements of the
//          original array p, p+1, p+2, ..., n-1 will be shifted to positions
//          p+1, p+2, p+3, ..., n)                            
//   a:     the value of the double to add                              
//-----------------------------------------------------------------------------
static double* ShivForDouble(double* V, int n, int p, double a)
{
  int i;

  // Reallocate the array
  V = (double*)realloc(V, (n+1)*sizeof(double));

  // Shift existing elements
  for (i = n; i > p; i--) {
    V[i] = V[i-1];
  }
  V[p] = a;

  return V;
}

//-----------------------------------------------------------------------------
// ShivForChar4: insert a string of four characters into an array at a  
//               particular point, shifting all other elements of the array one
//               step down.  The original array is reallocated (lengthened) and
//               then returned.                        
//                                                                      
// Arguments:                                                           
//   V:     the original array of characters (also returned)            
//   n:     n x 4 char is the original size of the array V              
//   p:     the position at which to add the new 4-character word       
//   a:     the 4-character word to add                                 
//-----------------------------------------------------------------------------
static char* ShivForChar4(char* V, int n, int p, char* a)
{
  int i;

  // Reallocate the array
  V = (char*)realloc(V, 4*(n+1)*sizeof(char));

  // Shift existing elements
  for (i = n; i > p; i--) {
    V[4*i] = V[4*(i-1)];
    V[4*i+1] = V[4*(i-1)+1];
    V[4*i+2] = V[4*(i-1)+2];
    V[4*i+3] = V[4*(i-1)+3];
  }
  for (i = 0; i < 4; i++) {
    V[4*p+i] = a[i];
  }

  return V;
}

//-----------------------------------------------------------------------------
// ThresholdIncrementor: increment integer n if it is greater than or   
//                       equal to the threshold thr; return the result. 
//-----------------------------------------------------------------------------
static int ThresholdIncrementor(int n, int thr)
{
  if (abs(n) >= thr) {
    if (n >= 0) {
      n += 1;
    }
    else {
      n -= 1;
    }
  }

  return n;
}

//-----------------------------------------------------------------------------
// UpdateBondIndices: update the atom indices for lists of bonds, angles, and
//                    dihedrals.
//
// Arguments:
//   B:      the list of bonds
//   A:      the list of angles
//   H:      the list of dihedrals
//   M:      the list of CMAPs
//   nb:     the size of the list of bonds
//   na:     the size of the list of angles
//   nh:     the size of the list of dihedrals
//   nm:     the size of the list of CMAPs
//   atmpos: the position of the extra point (atom) being inserted      
//-----------------------------------------------------------------------------
static void UpdateBondIndices(bond* B, int nb, angle* A, int na, dihedral* H,
                              int nh, cmapterm* M, int nm, int atmpos)
{
  int i;

  // Bonds
  for (i = 0; i < nb; i++) {
    B[i].a = ThresholdIncrementor(B[i].a, atmpos);
    B[i].b = ThresholdIncrementor(B[i].b, atmpos);
  }

  // Angles
  for (i = 0; i < na; i++) {
    A[i].a = ThresholdIncrementor(A[i].a, atmpos);
    A[i].b = ThresholdIncrementor(A[i].b, atmpos);
    A[i].c = ThresholdIncrementor(A[i].c, atmpos);
  }

  // Dihedrals
  for (i = 0; i < nh; i++) {
    H[i].a = ThresholdIncrementor(H[i].a, atmpos);
    H[i].b = ThresholdIncrementor(H[i].b, atmpos);
    H[i].c = ThresholdIncrementor(H[i].c, atmpos);
    H[i].d = ThresholdIncrementor(H[i].d, atmpos);
  }

  // CMAPs
  for (i = 0; i < nm; i++) {
    M[i].a = ThresholdIncrementor(M[i].a, atmpos);
    M[i].b = ThresholdIncrementor(M[i].b, atmpos);
    M[i].c = ThresholdIncrementor(M[i].c, atmpos);
    M[i].d = ThresholdIncrementor(M[i].d, atmpos);
    M[i].e = ThresholdIncrementor(M[i].e, atmpos);
  }
}

//-----------------------------------------------------------------------------
// UpdateBondArrays: updates the atom indices for lists of bonds, angles
//                   and dihedrals in the manner of UpdateBondIndices   
//                   above, but operating on lists used during dynamics.
//                                                                      
// Arguments:                                                           
//   tp:     the topology structure                                     
//   atmpos: the position of the extra point (atom) being inserted      
//-----------------------------------------------------------------------------
static void UpdateBondArrays(prmtop *tp, int atmpos)
{
  int i, j;

  // Lengthen the bond, angle, and dihedral arrays. 
  // The added extra points control no bonded terms.
  tp->BLC = (bondlist*)realloc(tp->BLC, tp->natom * sizeof(bondlist));
  tp->ALC = (angllist*)realloc(tp->ALC, tp->natom * sizeof(angllist));
  tp->HLC = (dihelist*)realloc(tp->HLC, tp->natom * sizeof(dihelist));
  tp->MLC = (cmaplist*)realloc(tp->MLC, tp->natom * sizeof(cmaplist));
  for (i = tp->natom-1; i > atmpos; i--) {
    tp->BLC[i] = tp->BLC[i-1];
    tp->ALC[i] = tp->ALC[i-1];
    tp->HLC[i] = tp->HLC[i-1];
    tp->MLC[i] = tp->MLC[i-1];
  }
  tp->BLC[atmpos].nbond = 0;
  tp->BLC[atmpos].BC = (bondcomm*)malloc(sizeof(bondcomm));
  tp->ALC[atmpos].nangl = 0;
  tp->ALC[atmpos].AC = (anglcomm*)malloc(sizeof(anglcomm));
  tp->HLC[atmpos].ndihe = 0;
  tp->HLC[atmpos].HC = (dihecomm*)malloc(sizeof(dihecomm));
  tp->MLC[atmpos].ncmap = 0;
  tp->MLC[atmpos].MC = (cmapcomm*)malloc(sizeof(cmapcomm));

  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->BLC[i].nbond; j++) {
      tp->BLC[i].BC[j].a = ThresholdIncrementor(tp->BLC[i].BC[j].a, atmpos);
      tp->BLC[i].BC[j].b = ThresholdIncrementor(tp->BLC[i].BC[j].b, atmpos);
    }
    for (j = 0; j < tp->ALC[i].nangl; j++) {
      tp->ALC[i].AC[j].a = ThresholdIncrementor(tp->ALC[i].AC[j].a, atmpos);
      tp->ALC[i].AC[j].c = ThresholdIncrementor(tp->ALC[i].AC[j].c, atmpos);
    }
    for (j = 0; j < tp->HLC[i].ndihe; j++) {
      tp->HLC[i].HC[j].a = ThresholdIncrementor(tp->HLC[i].HC[j].a, atmpos);
      tp->HLC[i].HC[j].c = ThresholdIncrementor(tp->HLC[i].HC[j].c, atmpos);
      tp->HLC[i].HC[j].d = ThresholdIncrementor(tp->HLC[i].HC[j].d, atmpos);
    }
    for (j = 0; j < tp->MLC[i].ncmap; j++) {
      tp->MLC[i].MC[j].a = ThresholdIncrementor(tp->MLC[i].MC[j].a, atmpos);
      tp->MLC[i].MC[j].b = ThresholdIncrementor(tp->MLC[i].MC[j].b, atmpos);
      tp->MLC[i].MC[j].d = ThresholdIncrementor(tp->MLC[i].MC[j].d, atmpos);
      tp->MLC[i].MC[j].e = ThresholdIncrementor(tp->MLC[i].MC[j].e, atmpos);
    }
  }
}

//-----------------------------------------------------------------------------
// UpdateExtraPointIndices: update the array of extra point atom IDs and frame
//                          atom IDs.                             
//                                                                      
// Arguments:                                                           
//   tp:         the topology structure                                 
//   atmpos:     the atom number of the extra point                     
//-----------------------------------------------------------------------------
static void UpdateExtraPointIndices(prmtop *tp, int atmpos)
{
  int i;

  // The loop runs over all but the most recently added extra point,
  // as the most recently added extra point's atom position is what 
  // causes atoms downstream to need their ID numbers incremented   
  // Any frame atoms of the most recently added extra point, by     
  // virtue of being within the same residue, will not be downstream
  // of the newly added extra point because the extra point is      
  // placed as the last atom in the residue.                        
  for (i = 0; i < tp->nxtrapt-1; i++) {
    tp->xtrapts[i].atomid = ThresholdIncrementor(tp->xtrapts[i].atomid,
                                                 atmpos);
    tp->xtrapts[i].fr1 = ThresholdIncrementor(tp->xtrapts[i].fr1, atmpos);
    tp->xtrapts[i].fr2 = ThresholdIncrementor(tp->xtrapts[i].fr2, atmpos);
    if (tp->xtrapts[i].frstyle > 1) {
      tp->xtrapts[i].fr3 = ThresholdIncrementor(tp->xtrapts[i].fr3, atmpos);
    }
    if (tp->xtrapts[i].frstyle == 6 || tp->xtrapts[i].frstyle == 9 ||
	tp->xtrapts[i].frstyle == 10) {
      tp->xtrapts[i].fr4 = ThresholdIncrementor(tp->xtrapts[i].fr4, atmpos);
    }
  }
}

//-----------------------------------------------------------------------------
// ExtendMap1234: extend a map1234 struct to catalog an additional 1:1, 1:2,
//                1:3, or 1:4 interaction.                         
//                                                                      
// Arguments:                                                           
//   imap:  the map1234 structure                                       
//   naatm: the atom which has a neighbor, nbatm
//   c:     indicates that a 1:c interaction is being added to the list 
//   nbatm: the atom that is a neighbor in a 1:c sense
//   epchk: flag to have exclusions checked for whether they have extra points
//          (turn this off to log exclusions involving extra point--that's a
//          separate pass over all atoms once exclusions involving atoms with
//          mass have been done)
//-----------------------------------------------------------------------------
static void ExtendMap1234(prmtop *tp, int naatm, int c, int nbatm, int epchk)
{
  int i;
  map1234 *imap;

  imap = &tp->nb1234[naatm];
  
  // Intercept extra points, things that have no mass.  They are, in fact, 1:1
  // to their parent atoms, 1:2 to atoms their parents are bonded to, etc.
  // Additional exclusions will be added to the list to accommodate these
  // particles.
  if (epchk == 1 &&
      (tp->Masses[naatm] < 1.0e-8 || tp->Masses[nbatm] < 1.0e-8)) {
    return;
  }
  
  // Check to see that this atom is not already a 1:1
  // neighbor, then extend the 1:1 list if requested
  for (i = 0; i < imap->n11; i++) {
    if (imap->L11[i] == nbatm) {
      return;
    }
  }
  if (c == 1) {
    imap->L11[imap->n11] = nbatm;
    imap->n11 += 1;
    if (imap->n11 == imap->mx11) {
      imap->mx11 += 5;
      imap->L11 = (int*)realloc(imap->L11, imap->mx11*sizeof(int));
    }
  }

  // Check to see that this atom is not already a 1:2
  // neighbor, then extend the 1:2 list if requested
  for (i = 0; i < imap->n12; i++) {
    if (imap->L12[i] == nbatm) {
      return;
    }
  }
  if (c == 2) {
    imap->L12[imap->n12] = nbatm;
    imap->n12 += 1;
    if (imap->n12 == imap->mx12) {
      imap->mx12 += 5;
      imap->L12 = (int*)realloc(imap->L12, imap->mx12*sizeof(int));
    }
  }

  // Check to see that this atom is not already a 1:3
  // neighbor, then extend the 1:3 list if requested
  for (i = 0; i < imap->n13; i++) {
    if (imap->L13[i] == nbatm) {
      return;
    }
  }
  if (c == 3) {
    imap->L13[imap->n13] = nbatm;
    imap->n13 += 1;
    if (imap->n13 == imap->mx13) {
      imap->mx13 += 5;
      imap->L13 = (int*)realloc(imap->L13, imap->mx13*sizeof(int));
    }
  }

  // Check to see that this atom is not already a 1:4
  // neighbor, then extend the 1:4 list if requested
  for (i = 0; i < imap->n14; i++) {
    if (imap->L14[i] == nbatm) {
      return;
    }
  }
  if (c == 4) {
    imap->L14[imap->n14] = nbatm;
    imap->n14 += 1;
    if (imap->n14 == imap->mx14) {
      imap->mx14 += 5;
      imap->L14 = (int*)realloc(imap->L14, imap->mx14*sizeof(int));
    }
  }
}

//-----------------------------------------------------------------------------
// AllocateMap1234: allocate a map1234 struct.                          
//-----------------------------------------------------------------------------
static void AllocateMap1234(map1234 *tmapi)
{
  tmapi->n11 = 0;
  tmapi->n12 = 0;
  tmapi->n13 = 0;
  tmapi->n14 = 0;
  tmapi->mx11 = 5;
  tmapi->mx12 = 10;
  tmapi->mx13 = 20;
  tmapi->mx14 = 40;
  tmapi->L11 = (int*)malloc(5*sizeof(int));
  tmapi->L12 = (int*)malloc(10*sizeof(int));
  tmapi->L13 = (int*)malloc(20*sizeof(int));
  tmapi->L14 = (int*)malloc(40*sizeof(int));
}

//-----------------------------------------------------------------------------
// Map13Neighbors: record 1:3 neighbors of atma in topology tp.         
//                                                                      
// Arguments:                                                           
//   atma:    the atom of interest
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
static void Map13Neighbors(int atma, prmtop *tp)
{
  int i, j, atmb, atmc;

  for (i = 0; i < tp->nb1234[atma].n12; i++) {
    atmb = tp->nb1234[atma].L12[i];
    for (j = 0; j < tp->nb1234[atmb].n12; j++) {
      atmc = tp->nb1234[atmb].L12[j];
      if (atmc != atma) {
        ExtendMap1234(tp, atma, 3, atmc, 1);
        ExtendMap1234(tp, atmc, 3, atma, 1);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Map14Neighbors: record 1:4 neighbors of atma in topology tp.         
//                                                                      
// Arguments:                                                           
//   atma:    the atom of interest
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
static void Map14Neighbors(int atma, prmtop *tp)
{
  int i, j, atmc, atmd;

  for (i = 0; i < tp->nb1234[atma].n13; i++) {
    atmc = tp->nb1234[atma].L13[i];
    for (j = 0; j < tp->nb1234[atmc].n12; j++) {
      atmd = tp->nb1234[atmc].L12[j];
      if (atmd != atma) {
        ExtendMap1234(tp, atma, 4, atmd, 1);
        ExtendMap1234(tp, atmd, 4, atma, 1);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Share1234Neighbors: record atoms A and B as being 1:1 neighbors, and 
//                     let atom B inherit all of atom A's 1:2, 1:3, and 
//                     1:4 neighbors, and add atom B to the appropriate 
//                     neighbor lists of all of atom A's neighbors.     
//                                                                      
// Arguments:                                                           
//   tp:      the topology                                              
//   A, B:    the atoms whose neighbors must now be shared, as one is an
//            extra point bound to the other
//-----------------------------------------------------------------------------
static void Share1234Neighbors(prmtop *tp, int A, int B)
{
  int i, An11, An12, An13, An14;
  int *AL11, *AL12, *AL13, *AL14;

  // Make atoms A and B 1:1 neighbors
  ExtendMap1234(tp, A, 1, B, 0);
  ExtendMap1234(tp, B, 1, A, 0);

  // Share 1:1 neighbors
  An11 = tp->nb1234[A].n11;
  AL11 = tp->nb1234[A].L11;
  for (i = 0; i < An11; i++) {
    if (AL11[i] != B) {
      ExtendMap1234(tp, B, 1, AL11[i], 0);
      ExtendMap1234(tp, AL11[i], 1, B, 0);
    }
  }

  // Share 1:2 neighbors
  An12 = tp->nb1234[A].n12;
  AL12 = tp->nb1234[A].L12;
  for (i = 0; i < An12; i++) {
    if (AL12[i] != B) {
      ExtendMap1234(tp, B, 2, AL12[i], 0);
      ExtendMap1234(tp, AL12[i], 2, B, 0);
    }
  }

  // Share 1:3 neighbors
  An13 =  tp->nb1234[A].n13;
  AL13 = tp->nb1234[A].L13;
  for (i = 0; i < An13; i++) {
    if (AL13[i] != B) {
      ExtendMap1234(tp, B, 3, AL13[i], 0);
      ExtendMap1234(tp, AL13[i], 3, B, 0);
    }
  }

  // Share 1:4 neighbors
  An14 =  tp->nb1234[A].n14;
  AL14 = tp->nb1234[A].L14;
  for (i = 0; i < An14; i++) {
    if (AL14[i] != B) {
      ExtendMap1234(tp, B, 4, AL14[i], 0);
      ExtendMap1234(tp, AL14[i], 4, B, 0);
    }
  }
}

//-----------------------------------------------------------------------------
// Connect1234: make a map of the topology, indicating for every atom   
//              what its 1:1, 1:2, 1:3, and 1:4 neighbors are.  The map 
//              is symmetric and reflexive, unlike the lists of bonds,  
//              angles, or dihedrals controlled by each atom.           
//                                                                      
// Arguments:                                                           
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
void Connect1234(prmtop *tp)
{
  int i, j, k;

  // Allocate
  tp->nb1234 = (map1234*)calloc(tp->natom, sizeof(map1234));
  for (i = 0; i < tp->natom; i++) {
    AllocateMap1234(&tp->nb1234[i]);
  }

  // Arrays have been created for bonds controlled by each atom, but
  // fall back on the original bond lists to maintain agreement with
  // Amber atom ordering, particularly when assigning EP frames.
  for (i = 0; i < tp->woH.nbond; i++) {
    ExtendMap1234(tp, tp->BNoH[i].a, 2, tp->BNoH[i].b, 1);
    ExtendMap1234(tp, tp->BNoH[i].b, 2, tp->BNoH[i].a, 1);
  }
  for (i = 0; i < tp->withH.nbond; i++) {
    ExtendMap1234(tp, tp->BIncH[i].a, 2, tp->BIncH[i].b, 1);
    ExtendMap1234(tp, tp->BIncH[i].b, 2, tp->BIncH[i].a, 1);
  }
  
  // Loop over bonds for 1:3 neighbors
  for (i = 0; i < tp->natom; i++) {
    Map13Neighbors(i, tp);
  }

  // Loop over bonds for 1:4 neighbors
  for (i = 0; i < tp->natom; i++) {
    Map14Neighbors(i, tp);
  }

  // Loop over all extra points in the topology and extend the exclusions list
  for (i = 0; i < tp->natom; i++) {
    if (fabs(tp->Masses[i]) < 1.0e-8) {
      for (j = 0; j < tp->natom; j++) {
        for (k = 0; k < tp->BLC[j].nbond; k++) {
          if (tp->BLC[j].BC[k].b == i) {
            Share1234Neighbors(tp, tp->BLC[j].BC[k].a, i);
          }
          if (tp->BLC[j].BC[k].a == i) {
            Share1234Neighbors(tp, tp->BLC[j].BC[k].b, i);
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// AddNewElimPair: function for adding a new 1:1, 1:2, 1:3, or 1:4 pair (whose
//                 interaction must be eliminated) to the growing list.
//                                                                      
// Arguments:                                                           
//   tp:         the topology structure                                 
//   atmH:       atom which will be in charge of the new elimination    
//   atm[XY]:    atoms with a pairwise interaction to eliminate (atmX is
//               going to be an extra point, by the way this routine    
//               gets called)                                           
//   ncon:       values 1 to 4 to specify 1:1, 1:2, 1:3, or 1:4         
//-----------------------------------------------------------------------------
void AddNewElimPair(prmtop *tp, int atmH, int atmX, int atmY, int ncon)
{
  int nXY;

  if (ncon == 1) {
    nXY = tp->ElimPair[atmH].n11;
    tp->ElimPair[atmH].list11 = (nixpr*)realloc(tp->ElimPair[atmH].list11,
                                               (nXY+1)*sizeof(nixpr));
    tp->ElimPair[atmH].list11[nXY].atmX = atmX;
    tp->ElimPair[atmH].list11[nXY].atmY = atmY;
    tp->ElimPair[atmH].n11 = nXY+1;
  }
  else if (ncon == 2) {
    nXY = tp->ElimPair[atmH].n12;
    tp->ElimPair[atmH].list12 = (nixpr*)realloc(tp->ElimPair[atmH].list12,
                                               (nXY+1)*sizeof(nixpr));
    tp->ElimPair[atmH].list12[nXY].atmX = atmX;
    tp->ElimPair[atmH].list12[nXY].atmY = atmY;
    tp->ElimPair[atmH].n12 = nXY+1;
  }
  else if (ncon == 3) {
    nXY = tp->ElimPair[atmH].n13;
    tp->ElimPair[atmH].list13 = (nixpr*)realloc(tp->ElimPair[atmH].list13,
                                               (nXY+1)*sizeof(nixpr));
    tp->ElimPair[atmH].list13[nXY].atmX = atmX;
    tp->ElimPair[atmH].list13[nXY].atmY = atmY;
    tp->ElimPair[atmH].n13 = nXY+1;
  }
  else if (ncon == 4) {
    nXY = tp->ElimPair[atmH].n14;
    tp->ElimPair[atmH].list14 = (nixpr*)realloc(tp->ElimPair[atmH].list14,
                                               (nXY+1)*sizeof(nixpr));
    tp->ElimPair[atmH].list14[nXY].atmX = atmX;
    tp->ElimPair[atmH].list14[nXY].atmY = atmY;
    tp->ElimPair[atmH].n14 = nXY+1;
  }
}

//-----------------------------------------------------------------------------
// IncrementEliminationList: this function increments the numbers in a list of
//                           pair eliminations, mindful of the fact that the
//                           integers in the list may be set to negative values
//                           for purposes of indicating whether the
//                           interactions involve electrostatics, van-der
//                           Waals, or both.    
//                                                                      
// Arguments:                                                           
//   prtmp:     the list of pair eliminations                          
//   n:         the number of elimiations in the list                  
//   atmpos:     the atom number of the newly added extra point         
//-----------------------------------------------------------------------------
static void IncrementEliminationList(nixpr *prtmp, int n, int atmpos)
{
  int i;

  for (i = 0; i < n; i++) {
    prtmp[i].atmX = ThresholdIncrementor(prtmp[i].atmX, atmpos);
    prtmp[i].atmY = ThresholdIncrementor(prtmp[i].atmY, atmpos);
  }
}

//-----------------------------------------------------------------------------
// Update1234Eliminations: function for listing new 1:2, 1:3, and 1:4   
//                         eliminations that may need to be applied when an
//                         extra point is added to a topology.  The exclusions
//                         are there, but electrostatic and van-der Waals
//                         interactions may still need to nullified post-hoc.
//                         These eliminations are accumulated immediately once
//                         new extra points are placed in order to avoid having
//                         to search for them all later and preventing double-
//                         counting of these eliminations.              
//                                                                      
// Arguments:                                                           
//   tp:         the topology structure                                 
//   atmpos:     the atom number of the extra point                     
//   ifr1:       the ID number of the extra point's parent atom         
//-----------------------------------------------------------------------------
static void Update1234Eliminations(prmtop *tp, int atmpos, int ifr1)
{
  int i, j, k, vdWTrue, vdWPossible, elecTrue, elecPossible, mpfound;
  int atmB, atmC, atmD;
  auxelim *tpi;

  // Update the eliminations list as all other lists have been
  tp->ElimPair = (auxelim*)realloc(tp->ElimPair, tp->natom*sizeof(auxelim));
  for (i = tp->natom-1; i > atmpos; i--) {
    tp->ElimPair[i] = tp->ElimPair[i-1];
  }
  tp->ElimPair[atmpos].n11 = 0;
  tp->ElimPair[atmpos].list11 = (nixpr*)malloc(sizeof(nixpr));
  tp->ElimPair[atmpos].n12 = 0;
  tp->ElimPair[atmpos].list12 = (nixpr*)malloc(sizeof(nixpr));
  tp->ElimPair[atmpos].n13 = 0;
  tp->ElimPair[atmpos].list13 = (nixpr*)malloc(sizeof(nixpr));
  tp->ElimPair[atmpos].n14 = 0;
  tp->ElimPair[atmpos].list14 = (nixpr*)malloc(sizeof(nixpr));
  for (i = 0; i < tp->natom; i++) {
    if (i == atmpos) {
      continue;
    }
    tpi = &tp->ElimPair[i];
    IncrementEliminationList(tpi->list11, tpi->n11, atmpos);
    IncrementEliminationList(tpi->list12, tpi->n12, atmpos);
    IncrementEliminationList(tpi->list13, tpi->n13, atmpos);
    IncrementEliminationList(tpi->list14, tpi->n14, atmpos);
  }

  // Now we can look at the eliminations for this extra point
  vdWPossible = (tp->LJIdx[atmpos] > 0);
  elecPossible = (fabs(tp->Charges[atmpos]) > 1.0e-8);

  // 1:1 "parent / sibling" interactions
  for (i = 0; i < tp->nb1234[atmpos].n11; i++) {
    atmB = tp->nb1234[atmpos].L11[i];
    elecTrue = (elecPossible && fabs(tp->Charges[atmB]) > 1.0e-8);
    vdWTrue = (vdWPossible && tp->LJIdx[atmB] >= 0);

    // The parent of the extra point shall control the eliminations
    if (elecTrue || vdWTrue) {
      AddNewElimPair(tp, ifr1, atmpos, atmB, 1);
    }
  }

  // 1:2 "bonded" interactions
  for (i = 0; i < tp->nb1234[atmpos].n12; i++) {
    atmB = tp->nb1234[atmpos].L12[i];
    elecTrue = (elecPossible && fabs(tp->Charges[atmB]) > 1.0e-8);
    vdWTrue = (vdWPossible && tp->LJIdx[atmB] >= 0);

    // The parent of the extra point shall control the eliminations
    if (elecTrue || vdWTrue) {
      AddNewElimPair(tp, ifr1, atmpos, atmB, 2);
    }
  }

  // 1:3 "angle neighbor" interactions
  for (i = 0; i < tp->nb1234[atmpos].n13; i++) {
    atmC = tp->nb1234[atmpos].L13[i];
    elecTrue = (elecPossible && fabs(tp->Charges[atmC]) > 1.0e-8);
    vdWTrue = (vdWPossible && tp->LJIdx[atmC] >= 0);
    if (vdWTrue || elecTrue) {

      // Find an atom which is 1:2 to atmpos and 1:2 to atmC,
      // a midpoint to which the interaction can be assigned.
      mpfound = 0;
      for (j = 0; j < tp->nb1234[atmpos].n12; j++) {
        atmB = tp->nb1234[atmpos].L12[j];
        for (k = 0; k < tp->nb1234[atmB].n12; k++) {
          if (tp->nb1234[atmB].L12[k] == atmC) {
            mpfound = 1;
            break;
          }
        }
        if (mpfound == 1) {
          AddNewElimPair(tp, atmB, atmpos, atmC, 3);
          break;
        }
      }
      if (mpfound == 0) {
        printf("Update1234Eliminations >> Error.  No intermediate atom found "
               "for 1:4 pair\nUpdate1234Eliminations >> %d and %d.\n", atmpos,
               atmC);
        exit(1);
      }
    }
  }

  // 1:4 nonbonded interactions
  for (i = 0; i < tp->nb1234[atmpos].n14; i++) {
    atmD = tp->nb1234[atmpos].L14[i];
    elecTrue = (elecPossible && fabs(tp->Charges[atmD]) > 1.0e-8);
    vdWTrue = (vdWPossible && tp->LJIdx[atmD] >= 0);
    if (vdWTrue || elecTrue) {

      // Find an atom which is 1:2 to atmpos and 1:3 to atmD,
      // a midpoint to which the interaction can be assigned.
      mpfound = 0;
      for (j = 0; j < tp->nb1234[atmpos].n12; j++) {
        atmB = tp->nb1234[atmpos].L12[j];
        for (k = 0; k < tp->nb1234[atmB].n13; k++) {
          if (tp->nb1234[atmB].L13[k] == atmD) {
            mpfound = 1;
            break;
          }
        }
        if (mpfound == 1) {
          AddNewElimPair(tp, atmB, atmpos, atmD, 4);
          break;
        }
      }
      if (mpfound == 0) {
        printf("ListNew14Interactions >> Error.  No intermediate atom found "
               "for 1:4 pair\nListNew14Interactions >> %d and %d.\n", atmpos,
               atmD);
        exit(1);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// UpdateConnectivity: update the list of map1234 structs assigned to each atom
//                     detailing 1:1, 1:2, 1:3, and 1:4 interactions.
//                                                                      
// Arguments:                                                           
//   tp:         the topology structure                                 
//   atmpos:     the atom number of the extra point                     
//   nxt:        the number of the extra point in the extra points list 
//   ert:        the extra point rule governing this extra point        
//-----------------------------------------------------------------------------
static void UpdateConnectivity(prmtop *tp, int atmpos, int nxt, eprule *ert)
{
  int i, j;

  // First, extend the connectivity array
  tp->nb1234 = (map1234*)realloc(tp->nb1234, tp->natom*sizeof(map1234));
  for (i = tp->natom-1; i > atmpos; i--) {
    tp->nb1234[i] = tp->nb1234[i-1];
  }

  // Fill in the connectivity for the
  // newcomer, the new extra point   
  AllocateMap1234(&tp->nb1234[atmpos]);

  // Loop over all atom neighbors and increment the numbering
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->nb1234[i].n12; j++) {
      if (tp->nb1234[i].L12[j] >= atmpos) {
        tp->nb1234[i].L12[j] += 1;
      }
    }
    for (j = 0; j < tp->nb1234[i].n13; j++) {
      if (tp->nb1234[i].L13[j] >= atmpos) {
        tp->nb1234[i].L13[j] += 1;
      }
    }
    for (j = 0; j < tp->nb1234[i].n14; j++) {
      if (tp->nb1234[i].L14[j] >= atmpos) {
        tp->nb1234[i].L14[j] += 1;
      }
    }
  }

  // The extra point being added is, by
  // definition, 1:1 to frame atom 1.  
  Share1234Neighbors(tp, tp->xtrapts[nxt].fr1, atmpos);

  // If frame atoms 2, 3, or 4 are also 1:1 to the 
  // extra point, add their neighbor lists as well.
  if (ert->excl2 == 1) {
    Share1234Neighbors(tp, tp->xtrapts[nxt].fr2, atmpos);
  }
  if (ert->excl3 == 1) {
    Share1234Neighbors(tp, tp->xtrapts[nxt].fr3, atmpos);
  }
  if (ert->excl4 == 1) {
    Share1234Neighbors(tp, tp->xtrapts[nxt].fr4, atmpos);
  }
}

//-----------------------------------------------------------------------------
// UpdateConstraintCommands: inserts a new, blank constraint command for this
//                           extra point.  Extra points cannot, by definition,
//                           control any constraint groups. 
//                                                                      
// Arguments:                                                           
//   tp:      the topology structure                                     
//   atmpos:  the position of the extra point now being added
//-----------------------------------------------------------------------------
static void UpdateConstraintCommands(prmtop *tp, int atmpos)
{
  int i, j;
  int *itmp;

  // Update all constraint group indices
  for (i = 0; i < tp->natom-1; i++) {
    if (tp->SHL[i].exe == 1) {
      itmp = tp->SHL[i].blist;
      for (j = 0; j < 2; j++) {
        itmp[j] = ThresholdIncrementor(itmp[j], atmpos);
      }
    }
    else if (tp->SHL[i].exe == 2) {
      itmp = &tp->SHL[i].blist[3*tp->SHL[i].blist[0]+2];
      for (j = 0; j < itmp[-1]; j++) {
        itmp[j] = ThresholdIncrementor(itmp[j], atmpos);
      }
    }
  }

  // Reallocate
  tp->SHL = (cnstcomm*)realloc(tp->SHL, tp->natom*sizeof(cnstcomm));
  for (i = tp->natom-1; i > atmpos; i--) {
    tp->SHL[i] = tp->SHL[i-1];
  }

  // Blank constraint command
  tp->SHL[atmpos].exe = 0;
  tp->SHL[atmpos].blist = (int*)calloc(1, sizeof(int)); 
}

//-----------------------------------------------------------------------------
// UpdateAtomGroups: increments the numbering of all atom groups in the 
//                   topology, and adds the newly inserted extra point to the
//                   appropriate group.  This is for purposes of reconnecting
//                   groups right before writing a restart file.
//
// Arguments:                                                           
//   tp:      the topology structure                                     
//   atmpos:  the position of the extra point now being added
//   fr1:     frame atom 1 of the newly added extra point
//-----------------------------------------------------------------------------
static void UpdateAtomGroups(prmtop *tp, int atmpos, int fr1)
{
  int i, j, grp4xpt;
  int *itmp;
  lgrp *tpg;

  // Increment group numbering
  grp4xpt = -1;
  for (i = 0; i < tp->ngrp; i++) {
    itmp = tp->lgrps[i].atoms;
    for (j = 0; j < tp->lgrps[i].natom; j++) {
      itmp[j] = ThresholdIncrementor(itmp[j], atmpos);

      // If the atom in this group is the extra point's  
      // frame atom 1, it will not have been incremeneted
      // but it should be recorded so that the new extra 
      // point can be added to this group.               
      if (itmp[j] == fr1) {
        grp4xpt = i;
      }
    }
  }

  // Add the extra point to the appropriate group
  if (grp4xpt < 0) {
    return;
  }
  tpg = &tp->lgrps[grp4xpt];
  j = -1;
  for (i = 0; i < tpg->natom; i++) {
    if (tpg->atoms[i] > atmpos) {
      j = i;
      tpg->atoms = ShivForInt(tpg->atoms, tpg->natom, i, atmpos);
      break;
    }
  }
  if (j == -1) {
    tpg->atoms = ShivForInt(tpg->atoms, tpg->natom, tpg->natom, atmpos);
  }
  tpg->natom += 1;
}

//-----------------------------------------------------------------------------
// AllocateElimList: allocate a pair eliminations list with blanks everywhere.
//                   Encapsulated in its own function to allow
//                   FulfillExclusions() in the Topology library to call it.
//                                                                      
// Arguments:                                                           
//   tp:         the topology structure                                 
//-----------------------------------------------------------------------------
void AllocateElimList(prmtop *tp)
{
  int i;

  tp->ElimPair = (auxelim*)malloc(tp->natom*sizeof(auxelim));
  for (i = 0; i < tp->natom; i++) {
    tp->ElimPair[i].n11 = 0;
    tp->ElimPair[i].list11 = (nixpr*)malloc(sizeof(nixpr));
    tp->ElimPair[i].n12 = 0;
    tp->ElimPair[i].list12 = (nixpr*)malloc(sizeof(nixpr));
    tp->ElimPair[i].n13 = 0;
    tp->ElimPair[i].list13 = (nixpr*)malloc(sizeof(nixpr));
    tp->ElimPair[i].n14 = 0;
    tp->ElimPair[i].list14 = (nixpr*)malloc(sizeof(nixpr));
  }
}

//-----------------------------------------------------------------------------
// PrepForNewEP: prepare a prmtop struct for new extra points.  This requires
//               allocation of certain auxiliary arrays, one of which is
//               retained throughout the dynamics, and setting the EPAdded
//               field of the prmtop to 1 to indicate that the prmtop has been
//               made ready for new extra points.   
//
// Arguments:                                                           
//   tp:         the topology structure                                 
//-----------------------------------------------------------------------------
static void PrepForNewEP(prmtop *tp)
{
  tp->EPInserted = 1;
  AllocateElimList(tp);
  
  // Record the original atom numbering
  tp->OldAtomNum = CountUp(tp->natom);
  tp->norigatom = tp->natom;
}

//-----------------------------------------------------------------------------
// InsertExtraPoint: inserts an extra point into a topology struct--no simple
//                   task.  The extra point will be inserted as prescribed by
//                   an extra point rule, and all appropriate exclusions will
//                   be added.  An insertion requires significant migration of
//                   information in the topology struct, as each new extra
//                   point will be inserted as the last atom of an existing
//                   residue which then requires all atom ID numbers in
//                   subsequent residues to be incremented.
//
// Arguments:
//   tp:           the topology data structure
//   resid:        index of the residue receiving the extra point
//   ruleid:       the index of the rule that the extra point will follow
//   molAtomLims:  atom index limits of each molecule 
//-----------------------------------------------------------------------------
static void InsertExtraPoint(prmtop *tp, int resid, int ruleid,
                             int* molAtomLims)
{
  int i, atmpos, ifr1, ifr2, ifr3, ifr4, nxt, nmol, stride, dosearch;
  eprule *ert;

  // Indicate that extra points have been inserted.
  // This happens exactly once, when the first new 
  // extra point is inserted.                      
  if (tp->EPInserted == 0) {
    PrepForNewEP(tp);
  }

  // Set a pointer to the extra point rule
  ert = &tp->eprules[ruleid];

  // Locate the frame atoms.  These can be located now and
  // not updated, as the frame atoms are all required to  
  // be within the residue and the extra point is added at
  // the end of the residue, ensuring that the frame atom 
  // numbers will be unchanged.                           
  ifr1 = FindAtom(tp, tp->ResLims[resid], tp->ResLims[resid+1], ert->fr1);
  ifr2 = FindAtom(tp, tp->ResLims[resid], tp->ResLims[resid+1], ert->fr2);
  if (ert->frstyle > 1 && ert->frstyle != 11) {
    ifr3 = FindAtom(tp, tp->ResLims[resid], tp->ResLims[resid+1], ert->fr3);
  }
  if (ert->frstyle == 6 || ert->frstyle == 9 || ert->frstyle == 10) {
    ifr4 = FindAtom(tp, tp->ResLims[resid], tp->ResLims[resid+1], ert->fr4);
  }

  // The new extra point counts as a new atom
  // tacked on at the end of the residue.    
  atmpos = tp->ResLims[resid+1];

  // The new extra point can now be marked
  // for processing during dynamics.      
  nxt = tp->nxtrapt;
  tp->xtrapts[nxt].frstyle = ert->frstyle;
  tp->xtrapts[nxt].fr1 = ifr1;
  tp->xtrapts[nxt].fr2 = ifr2;
  tp->xtrapts[nxt].d1 = ert->d1;
  tp->xtrapts[nxt].atomid = atmpos;
  if (ert->frstyle == 1 || ert->frstyle == 11) {

    // If the frame style is 1 or 11, it's just two
    // atoms and the one distance parameter. 
    tp->xtrapts[nxt].nfratm = 2;
  }
  else {

    // There is guaranteed to be a second distance or
    // angle parameter, and at least a third atom.   
    tp->xtrapts[nxt].d2 = ert->d2;
    tp->xtrapts[nxt].fr3 = ifr3;
    if (ert->frstyle < 6) {

      // There are three atoms in the frame.
      tp->xtrapts[nxt].nfratm = 3;

      // A third distance parameter is needed.
      if (ert->frstyle == 3 || ert->frstyle == 5) {
        tp->xtrapts[nxt].d3 = ert->d3;
      }
    }
    else {

      // A fourth atom is required, as well
      // as a fourth distance parameter.   
      tp->xtrapts[nxt].fr4 = ifr4;
      tp->xtrapts[nxt].nfratm = 4;
      tp->xtrapts[nxt].d3 = ert->d3;
    }
  }
  tp->nxtrapt += 1;
  ReallocateXtraPtsArray(tp, 1);
  
  // Update the number of atoms in the topology.
  tp->natom += 1;

  // Update the residue limits
  for (i = resid+1; i <= tp->nres; i++) {
    tp->ResLims[i] += 1;
  }
  
  // Update the molecular atom index limits
  if (tp->ifbox > 0) {
    stride = tp->nspm / 4;
    if (stride < 1) {
      stride = 1;
    }
    nmol = tp->nspm / 2;
    dosearch = 1;
    while (dosearch == 1) {
      if (atmpos <= molAtomLims[nmol]) {
        nmol -= stride;
      }
      else if (atmpos > molAtomLims[nmol+1]) {
        nmol += stride;
      }
      else {
        dosearch = 0;
      }
      stride /= 2;
      if (stride < 1) {
        stride = 1;
      }
    }
    for (i = nmol + 1; i <= tp->nspm; i++) {
      molAtomLims[i] += 1;
    }
    tp->Nsp[nmol] += 1;
  }
  
  // Reallocate integer arrays as needed
  tp->LJIdx = ShivForInt(tp->LJIdx, tp->natom - 1, atmpos, ert->LJIdx);
  tp->Join = ShivForInt(tp->Join, tp->natom - 1, atmpos, ert->Join);
  tp->Rotat = ShivForInt(tp->Rotat, tp->natom - 1, atmpos, ert->Rotat);
  tp->ZNumber = ShivForInt(tp->ZNumber, tp->natom - 1, atmpos, 0);  
  
  // Reallocate double precision arrays as needed
  tp->Charges = ShivForDouble(tp->Charges, tp->natom - 1, atmpos, ert->Charge);
  tp->OrigCharges = ShivForDouble(tp->OrigCharges, tp->natom - 1, atmpos,
                                  ert->Charge * 18.2223);
  tp->Masses = ShivForDouble(tp->Masses, tp->natom - 1, atmpos, 0.0);
  tp->InvMasses = ShivForDouble(tp->InvMasses, tp->natom - 1, atmpos, 0.0);
  tp->Radii = ShivForDouble(tp->Radii, tp->natom - 1, atmpos, 0.0);
  tp->Screen = ShivForDouble(tp->Screen, tp->natom - 1, atmpos, 0.0);

  // Reallocate character arrays as needed
  tp->AtomNames = ShivForChar4(tp->AtomNames, tp->natom - 1, atmpos,
                               ert->epname);
  tp->AtomTypes = ShivForChar4(tp->AtomTypes, tp->natom - 1, atmpos, "EP  ");
  tp->TreeSymbols = ShivForChar4(tp->TreeSymbols, tp->natom - 1, atmpos,
                                 "BLA ");

  // Increment the atom numbering
  tp->OldAtomNum = ShivForInt(tp->OldAtomNum, tp->natom - 1, atmpos, -1);

  // Increment bond identification numbers
  UpdateBondIndices(tp->BIncH, tp->withH.nbond, tp->AIncH, tp->withH.nangl,
                    tp->HIncH, tp->withH.ndihe, tp->Mterms, tp->ncmap, atmpos);
  UpdateBondIndices(tp->BNoH, tp->woH.nbond, tp->ANoH, tp->woH.nangl,
                    tp->HNoH, tp->woH.ndihe, NULL, 0, atmpos);
  UpdateBondArrays(tp, atmpos);

  // Increment extra point atom ID numbers and frames
  UpdateExtraPointIndices(tp, atmpos);

  // Increment the constraint records
  UpdateConstraintCommands(tp, atmpos);

  // Increment atom groups
  UpdateAtomGroups(tp, atmpos, ifr1);

  // Increment connectivity records
  UpdateConnectivity(tp, atmpos, nxt, ert);

  // Accumulate lists of 1:1, 1:2, 1:3, and 1:4 eliminations
  Update1234Eliminations(tp, atmpos, ifr1);
}

//-----------------------------------------------------------------------------
// ModifyAtomProperties: alter the properties of an atom in the original
//                       topology.  This routine does not add extra     
//                       points.                                        
//
// Arguments:
//   tp:      the topology
//   resid:   index of the residue receiving the extra point
//   ruleid:  the index of the rule that the extra point will follow
//-----------------------------------------------------------------------------
static void ModifyAtomProperties(prmtop *tp, int resid, int ruleid)
{
  int i, atmid;
  eprule *ert;

  // Locate the atom within the residue
  ert = &tp->eprules[ruleid];
  atmid = -1;
  for (i = tp->ResLims[resid]; i < tp->ResLims[resid+1]; i++) {
    if (strncmp(&tp->AtomNames[4*i], ert->epname, 4) == 0) {
      atmid = i;
      break;
    }
  }
  if (atmid == -1) {
    printf("ModifyAtomProperties >> Error.  Atom %.4s not found in residue "
           "%d.\n", ert->epname, resid);
  }
  if (ert->modElec == 1) {
    tp->Charges[atmid] = ert->Charge;
    tp->OrigCharges[atmid] = ert->Charge * 18.2223;
  }
  if (ert->modLJ == 1) {
    tp->LJIdx[atmid] = ert->LJIdx;
  }
}

//-----------------------------------------------------------------------------
// DetectType12EP: detect "type I" or "type II" extra points according to
//                 Amber conventions.  If there is one point, it will lie in
//                 the plane of the three frame atoms, along the line defined
//                 by the parent atom and the bisector of the other two frame
//                 atoms.
//
// Arguments:
//   tp:     the topology structure                                     
//   epID:   the ID number of the extra point
//-----------------------------------------------------------------------------
static int DetectType12EP(prmtop *tp, int epID)
{
  int i, nWithMass, parAtom;

  // Find out the parent atom.
  nWithMass = 0;
  for (i = 0; i < tp->nb1234[epID].n11; i++) {
    nWithMass += (tp->Masses[tp->nb1234[epID].L11[i]] > 1.0e-8);
  }
  if (nWithMass > 1) {
    printf("DetectType12EP >> Error.  Extra point %d bonds to %d massive "
           "atoms.\n", epID + 1, nWithMass);
    exit(1);
  }
  parAtom = tp->nb1234[epID].L11[0];

  // If the parent atom has bonds to two other atoms with mass, that will
  // classify this extra point as type I.  If there is only one other atom,
  // that classifies this extra point as type II.
  if (tp->nb1234[parAtom].n12 == 2) {
    return 1;
  }
  else if (tp->nb1234[parAtom].n12 == 1) {
    return 2;
  }
  else {
    printf("DetectType12EP >> Error.  Extra point %d is bonded to an atom "
           "with\nDetectType12EP >> %d surrogates.\n", epID + 1,
           tp->nb1234[parAtom].n12);
    exit(1);
    
    return 0;
  }
}

//-----------------------------------------------------------------------------
// DetermineEPFrames: determine the frame atoms, and frame type, for all
//                    massless atoms in the system.  The basic rules are
//                    that a frame must be found for every extra point, 
//                    and that the frame atoms must lie within the same 
//                    residue.  While extra points are initialized on a 
//                    residue-by-residue basis, the frames are stored on
//                    an atom-by-atom basis so that calculations do not 
//                    need to reference whole residues.                 
//                                                                      
// Arguments:                                                           
//   tp:     the topology structure                                     
//-----------------------------------------------------------------------------
void DetermineEPFrames(prmtop *tp)
{
  int i, j, k, epid, nep, nresatm, epfound, thereAreTwoEP, thisIsLast;
  int* eplist;
  int* rulefound;
  int* molAtomLims;

  // Initialize data pertaining to extra points
  tp->nxtrapt = 0;
  tp->maxtrapt = 32;
  tp->xtrapts = (expt*)malloc(tp->maxtrapt * sizeof(expt));
  ReallocateXtraPtsArray(tp, 1);
  
  // Loop over all residues and search for    
  // existing extra points or applicable rules
  for (i = 0; i < tp->nres; i++) {

    // Count any extra points that this residue may already have
    nep = 0;
    nresatm = tp->ResLims[i+1] - tp->ResLims[i];
    eplist = (int*)calloc(nresatm, sizeof(int));
    rulefound = (int*)calloc(nresatm, sizeof(int));
    for (j = tp->ResLims[i]; j < tp->ResLims[i+1]; j++) {
      if (tp->Masses[j] < 1.0e-8) {
        eplist[nep] = j;
        nep++;
	ReallocateXtraPtsArray(tp, 1);
      }
    }

    // To begin, we check to see if there is a user-defined rule
    // for dealing with this extra point.  If there is no such  
    // rule, then the program's built-in rules take effect.  In 
    // this manner, mdgx allows users to define their own extra 
    // point rules that can override the default behavior.      
    for (j = 0; j < nep; j++) {
      epid = EPRuleDefined(tp, i, eplist[j]);
      if (epid > -1) {
        ImplementEPRule(tp, i, epid);
        rulefound[j] = 1;
      }
    }

    // Check the custom virtual site specifications that came
    // with the topology.  If one describes the extra point,
    // translate it into mdgx's internal conventions.
    if (tp->nCustomVS > 0) {
      for (j = 0; j < nep; j++) {
        for (k = 0; k < tp->nCustomVS; k++) {
          if (tp->customVSi[6*k] - 1 == eplist[j]) {
	    AmberFrame2mdgx(tp, k);
	    rulefound[j] = 1;
	  }
	}
      }
    }
    
    // If there is only a bond between the extra point and the
    // oxygen atom, then we will assume that the extra point  
    // lies along the bisector between the O-H1 and O-H2 bonds
    // and that this is a TIP4P-like water molecule.          
    else if (nresatm == 4 && nep == 1 && rulefound[0] == 0 &&
        CountBonds(tp, eplist[0]) == 1 &&
        strncmp(&tp->ResNames[4*i], tp->WaterName, 4) == 0) {
      Implement4PWaterRule(tp, i);
      rulefound[0] = 1;
    }

    // If there are two extra points, five atoms in the residue,
    // the residue is what the topology identifies as water, and
    // this extra point has only one bond, then we will assume  
    // that this extra point is part of a TIP5P-like water      
    // molecule.                                                
    else if (nresatm == 5 && nep == 2 && rulefound[0] == 0 &&
             rulefound[1] == 0 && CountBonds(tp, eplist[0]) == 1 &&
             CountBonds(tp, eplist[1]) == 1 &&
             strncmp(&tp->ResNames[4*i], tp->WaterName, 4) == 0) {
      Implement5PWaterRule(tp, i);
      rulefound[0] = 1;
      rulefound[1] = 1;
    }

    // Otherwise, the extra point must be ascertained by the number
    // and types of atoms that its parent atom binds.
    else {
      for (j = 0; j < nep ; j++) {
        if (rulefound[j] == 1) {
          continue;
        }

        // Before implementing the rule, it is important to see if
        // there are other extra points bonded to the same atom.
        // There can be at most two under the current Amber rules.
        thereAreTwoEP = 0;
        thisIsLast = 0;
        for (k = 0; k < nep; k++) {
	  if (tp->nb1234[eplist[j]].n11 == 0) {
	    continue;
	  }
          if (tp->nb1234[eplist[j]].L11[0] == tp->nb1234[eplist[k]].L11[0]) {
            thereAreTwoEP++;
            if (k < j) {
              thisIsLast++;
            }
          }
        }
	if (thereAreTwoEP == 0) {

	  // This is a sign that the extra point didn't bond to anything
	  // else, and is therefore probably a rogue particle in an
	  // incomplete topology.
	  printf("DetermineEPFrames >> Warning.  No atoms were found bonding "
		 "to atom %4d (%4.4s).\nDetermineEPFrames >> This is often a "
		 "sign of an incomplete topology, which\nDetermineEPFrames >> "
		 "has extra points in need of custom frame specifications.\n",
		 eplist[j], &tp->AtomNames[4*eplist[j]]);
	  continue;
	}
        thereAreTwoEP = (thereAreTwoEP == 2) ? 1 : 0;
        if (DetectType12EP(tp, eplist[j]) == 1) {
          ImplementType1EP(tp, eplist[j], thereAreTwoEP * thisIsLast,
                           thereAreTwoEP + 1);
        }
        else {
          ImplementType2EP(tp, eplist[j], thereAreTwoEP * thisIsLast,
                           thereAreTwoEP + 1);
        }
        rulefound[j] = 1;
      }
    }
    
    // If all extra points have not been accounted
    // for, it is unclear what to do              
    if (ISum(rulefound, nep) < nep) {
      printf("DetermineEPFrames >> Error.  No rules describe %d extra "
             "points in a residue:\n", nep - ISum(rulefound, nep));
      for (j = tp->ResLims[i]; j < tp->ResLims[i+1]; j++) {
        printf("DetermineEPFrames >>     %6d %.4s %5d %.4s\n", j,
               &tp->AtomNames[4*j], i, &tp->ResNames[4*i]);
      }
      exit(1);
    }

    // Free allocated memory
    free(eplist);
    free(rulefound);
  }
  
  // In addition to extra points that may already be included in the topology
  // file, any extra point rules provided by the user must be implemented.
  // The next step is to go through each rule applicable to residue i, and add
  // extra points as required by the rules.  In this manner, the mdgx program
  // has the capability of expanding a topology after reading it from a file.
  // Start by making an array of all molecule limits (similar to tp->ResLims)
  // to help track which molecule atom counts to increment as new extra points
  // get added.
  if (tp->ifbox > 0) {
    molAtomLims = (int*)calloc(tp->nspm + 1, sizeof(int));
    j = 0;
    for (i = 0; i < tp->nspm; i++) {
      molAtomLims[i] = j;
      j += tp->Nsp[i];
    }
    molAtomLims[tp->nspm] = tp->natom;
  }
  for (i = 0; i < tp->nres; i++) {
    for (j = 0; j < tp->neprule; j++) {

      // If this extra point rule does not apply to the residue, skip
      if (strncmp(tp->eprules[j].resname, &tp->ResNames[4*i], 4) != 0) {
        continue;
      }

      // Check to see if this extra point is already present
      epfound = 0;
      for (k = tp->ResLims[i]; k < tp->ResLims[i+1]; k++) {
        if (strncmp(tp->eprules[j].epname, &tp->AtomNames[4*k], 4) == 0) {
          epfound = 1;
          break;
        }
      }
      if (epfound == 1 && tp->eprules[j].frstyle > 0) {
        continue;
      }

      // Modify atom parameters
      else if (epfound == 1 && tp->eprules[j].frstyle == 0) {
        ModifyAtomProperties(tp, i, j);
      }

      // Add the extra point if it is not already present
      else {
        InsertExtraPoint(tp, i, j, molAtomLims);
      }
    }
  }
  if (tp->ifbox > 0) {
    free(molAtomLims);
  }

  // Trim the extra points array to fit the actual size needed.
  // The number of extra points has been accumulated during the
  // loops over all residues, but the total allocated memory is
  // almost certainly more than is needed.                     
  if (tp->nxtrapt < tp->maxtrapt) {
    if (tp->nxtrapt > 0) {
      tp->xtrapts = (expt*)realloc(tp->xtrapts, tp->nxtrapt*sizeof(expt));
    }
    else {
      tp->xtrapts = (expt*)realloc(tp->xtrapts, sizeof(expt));
    }
  }

  // Also store the number of extra points in the slot that the rest of
  // Amber will recognize
  tp->numextra = tp->nxtrapt;
}

//-----------------------------------------------------------------------------
// XptLocator: encapsulates the function for placing an extra point given the
//             frame atoms and frame type.                    
//                                                                      
// Arguments:                                                           
//   [a,b,c,d]loc:  the locations of frame atoms                        
//   agloc:         the global instance of the first frame atom         
//   ep[g]loc:      the locations of the extra point and its global     
//                  instance                                            
//   tmr:           the frame data                                      
//-----------------------------------------------------------------------------
void XptLocator(double *aloc, double *agloc, double *bloc, double *cloc,
                double *dloc, double *eploc, double *epgloc, expt *tmr)
{
  int i;
  double dvx, dvy, dvz, magdvec, invab, invab2, magrP2, abbcOabab, delfr1;
  double rabfac, rPfac, magrab, magrac;
  double rab[3], rac[3], rad[3], rbc[3], rbd[3], abac[3], rPerp[3];
  double rja[3], rjb[3], rm[3];

  // Compute displacements
  for (i = 0; i < 3; i++) {
    rab[i] = bloc[i] - aloc[i];
  }
  if (tmr->frstyle > 1 && tmr->frstyle < 9) {
    for (i = 0; i < 3; i++) {
      rac[i] = cloc[i] - aloc[i];
      rbc[i] = cloc[i] - bloc[i];
    }
  }
  if (tmr->frstyle == 9 || tmr->frstyle == 10) {
    for (i = 0; i < 3; i++) {
      rbc[i] = 0.5*(cloc[i] + bloc[i]);
      rbd[i] = 0.5*(dloc[i] + bloc[i]);
      rab[i] = rbc[i] - aloc[i];
      rac[i] = rbd[i] - aloc[i];
    }
  }
  if (tmr->frstyle == 6) {
    for (i = 0; i < 3; i++) {
      rad[i] = dloc[i] - aloc[i];
      rja[i] = (tmr->d1 * rac[i]) - rab[i];
      rjb[i] = (tmr->d2 * rad[i]) - rab[i];
    }
  }
  
  // Place extra point.  Try the GROMACS styles first.
  if (tmr->frstyle == 1) {
    for (i = 0; i < 3; i++) {
      delfr1 = rab[i]*tmr->d1;
      eploc[i] = aloc[i] + delfr1;
      epgloc[i] = agloc[i] + delfr1;
    }
  }
  else if (tmr->frstyle == 11) {
    invab = tmr->d1 / sqrt(rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2]);
    for (i = 0; i < 3; i++) {
      eploc[i] = aloc[i] + (rab[i] * invab);
      epgloc[i] = agloc[i] + (rab[i] * invab);
    }
  }
  else if (tmr->frstyle == 2) {
    for (i = 0; i < 3; i++) {
      delfr1 = rab[i]*tmr->d1 + rac[i]*tmr->d2;
      eploc[i] = aloc[i] + delfr1;
      epgloc[i] = agloc[i] + delfr1;
    }
  }
  else if (tmr->frstyle == 3) {
    dvx = rab[0] + tmr->d3*rbc[0];
    dvy = rab[1] + tmr->d3*rbc[1];
    dvz = rab[2] + tmr->d3*rbc[2];
    magdvec = 1.0/sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
    for (i = 0; i < 3; i++) {
      delfr1 = tmr->d1*(rab[i] + tmr->d3*rbc[i])*magdvec;
      eploc[i] = aloc[i] + delfr1;
      epgloc[i] = agloc[i] + delfr1;
    }
  }
  else if (tmr->frstyle == 4) {
    invab2 = 1.0/Dot3(rab, rab);
    abbcOabab = Dot3(rab, rbc) * invab2;
    rPerp[0] = rbc[0] - abbcOabab*rab[0];
    rPerp[1] = rbc[1] - abbcOabab*rab[1];
    rPerp[2] = rbc[2] - abbcOabab*rab[2];
    magrP2 = Dot3(rPerp, rPerp);
    rabfac = tmr->d1*cos(tmr->d2)*sqrt(invab2);
    rPfac = tmr->d1*sin(tmr->d2)/sqrt(magrP2);
    for (i = 0; i < 3; i++) {
      delfr1 = rabfac*rab[i] + rPfac*rPerp[i];
      eploc[i] = aloc[i] + delfr1;
      epgloc[i] = agloc[i] + delfr1;
    }
  }
  else if (tmr->frstyle == 5) {
    CrossP(rab, rac, abac);
    for (i = 0; i < 3; i++) {
      delfr1 = rab[i]*tmr->d1 + rac[i]*tmr->d2 + abac[i]*tmr->d3;
      eploc[i] = aloc[i] + delfr1;
      epgloc[i] = agloc[i] + delfr1;
    }
  }
  else if (tmr->frstyle == 6) {
    CrossP(rja, rjb, rm);
    magdvec = tmr->d3 / sqrt(rm[0]*rm[0] + rm[1]*rm[1] + rm[2]*rm[2]);
    for (i = 0; i < 3; i++) {
      eploc[i] = aloc[i] + (magdvec * rm[i]);
      epgloc[i] = agloc[i] + (magdvec * rm[i]);
    }
  }

  // Amber EP frame styles cover chemical space, but are less general.
  else if (tmr->frstyle >= 7 && tmr->frstyle <= 10) {

    // This basic operation is common to all Amber frame types, once
    // the "B" and "C" particles have been adapted for frames 9 and 10.
    magrab = 1.0 / sqrt(rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2]);
    magrac = 1.0 / sqrt(rac[0]*rac[0] + rac[1]*rac[1] + rac[2]*rac[2]);
    for (i = 0; i < 3; i++) {
      rab[i] *= magrab;
      rac[i] *= magrac;
      abac[i] = -0.5*(rab[i] + rac[i]);
    }
    delfr1 = tmr->d1 /
             sqrt(abac[0]*abac[0] + abac[1]*abac[1] + abac[2]*abac[2]);

    // A single, planar type I or type II extra point can now be placed.
    for (i = 0; i < 3; i++) {
      eploc[i] = aloc[i] + delfr1*abac[i];      
      epgloc[i] = agloc[i] + delfr1*abac[i];
    }

    // For tetrahedral type I extra points, move a distance d1 along
    // the local bisector in the plane of the three frame atoms, then
    // a distance d2 along the plane normal.  The vectors rab and rac
    // have already been normalized, so their cross product is also a
    // unit vector.
    if (tmr->frstyle == 8) {
      CrossP(rab, rac, abac);
      magdvec = sqrt(abac[0]*abac[0] + abac[1]*abac[1] + abac[2]*abac[2]);
      abac[0] /= magdvec;
      abac[1] /= magdvec;
      abac[2] /= magdvec;
      for (i = 0; i < 3; i++) {
        eploc[i] += tmr->d2 * abac[i];
        epgloc[i] += tmr->d2 * abac[i];
      }
    }

    // For twin type II extra points, move perpendicular to the
    // previous axis, that is parallel to the direction of B->C.
    else if (tmr->frstyle == 10) {
      for (i = 0; i < 3; i++) {
        abac[i] = rac[i] - rab[i];
      }
      delfr1 = tmr->d2 /
        sqrt(abac[0]*abac[0] + abac[1]*abac[1] + abac[2]*abac[2]);
      for (i = 0; i < 3; i++) {
        eploc[i] += delfr1 * abac[i];
        epgloc[i] += delfr1 * abac[i];        
      }
    }
  }

  // Failure to place the extra point.  Exit with error.
  else {
    printf("XptLocator >> Error.  Not prepared to handle frame type %d.\n",
           tmr->frstyle);
    exit(1);
  }
}

//-----------------------------------------------------------------------------
// MolPlaceXpt: place extra points around a molecule given its topology.
//
// Arguments:
//   crd:      coordinates of the molecule (extra point sites are
//             included but all set to zero)
//   tp:       topology decribing the molecule
//-----------------------------------------------------------------------------
void MolPlaceXpt(coord *crd, prmtop *tp)
{
  int i;
  double *loctmp, *aloc, *bloc, *cloc, *dloc, *eploc;
  expt *tmr;

  // Loop over all extra points and place them
  for (i = 0; i < tp->nxtrapt; i++) {
    tmr = &tp->xtrapts[i];
    loctmp = crd->loc;
    aloc = &loctmp[3*tmr->fr1];
    bloc = &loctmp[3*tmr->fr2];
    if (tmr->frstyle > 1) {
      cloc = &loctmp[3*tmr->fr3];
    }
    if (tmr->frstyle == 6 || tmr->frstyle == 9 || tmr->frstyle == 10) {
      dloc = &loctmp[3*tmr->fr4];
    }
    eploc = &loctmp[3*tmr->atomid];
    XptLocator(aloc, aloc, bloc, cloc, dloc, eploc, eploc, tmr);
  }
}


//-----------------------------------------------------------------------------
// CellPlaceXpt: places an extra point based on information about its frame.
//               This function operates on atoms found within cell structs.
//                                                                      
// Arguments:                                                           
//   tp:     the topology structure                                     
//   crd:    the system coordinates                                     
//   C:      the cell of interest                                       
//   CG:     the cell grid of interest                                  
//-----------------------------------------------------------------------------
static void CellPlaceXpt(prmtop *tp, coord *crd, cell *C, cellgrid *CG)
{
  int h, i, j, k, bid, cid, did, belem, celem, delem, csector;
  int breg, creg, dreg;
  double dvx, dvy, dvz, magdvec, invab2, magrP2, abbcOabab, delfr1;
  double rabfac, rPfac;
  double rab[3], eploc[3];
  double *aloc, *bloc, *cloc, *dloc, *ccen, *agloc, *epgloc, *dbng;
  expt *tmr;
  atomc *atmhi, *atmjk;
  lgrp *hiFR1;

  // Prepare this cell to receive extra points
  SetIVec(&C->nr[8], 8, 0);

  // Pointers
  ccen = &CG->celldim[4];
  dbng = CG->dbng;

  // Loop over all cell regions 
  for (h = 0; h < 8; h++) {
    for (i = 0; i < C->nr[h]; i++) {

      // Test this atom to see if it is    
      // frame atom 1 of any extra point(s)
      atmhi = &C->map[h][i];
      hiFR1 = &tp->FR1Idx[atmhi->id];
      if (hiFR1->natom == 0 ||
          IsCentralAtom(C->orig, ccen, atmhi->loc,
                        crd->U.data, dbng, crd->isortho) < 0) {
        continue;
      }

      // This atom is frame atom 1 to one or 
      // more extra points; place them one by
      // one and then add them to the proper 
      // regions of this cell.               
      for (j = 0; j < hiFR1->natom; j++) {

        // Pointer to extra point structure
        tmr = &tp->xtrapts[hiFR1->atoms[j]];

        // Pointers to atom locations
        aloc = atmhi->loc;
        bid = C->GPSptr[tmr->fr2];
        breg = bid/C->maxatom;
        belem = bid - breg*C->maxatom;
        bloc = C->map[breg][belem].loc;
        if (tmr->frstyle != 1) {
          cid = C->GPSptr[tmr->fr3];
          creg = cid/C->maxatom;
          celem = cid - creg*C->maxatom;
          cloc = C->map[creg][celem].loc;
        }
        if (tmr->frstyle == 6 || tmr->frstyle == 9 || tmr->frstyle == 10) {
          did = C->GPSptr[tmr->fr4];
          dreg = did/C->maxatom;
          delem = did - dreg*C->maxatom;
          dloc = C->map[dreg][delem].loc;
        }        
        agloc = &crd->loc[3*atmhi->id];
        epgloc = &crd->loc[3*tmr->atomid];
        XptLocator(aloc, agloc, bloc, cloc, dloc, eploc, epgloc, tmr);
        
        // The location of the extra point is now known.   
        // Determine what sector of the cell it lies in.   
        // Add the extra point as a new site in that secor.
        // The extra 8 indices of cell C's nr array are    
        // used to store the number of extra points added  
        // to each sector in this manner.                  
        for (k = 0; k < 3; k++) {
          rab[k] = eploc[k] - C->orig[k];
        }
        RotateCrd(rab, 1, crd->U);
        for (k = 0; k < 3; k++) {
          rab[k] = floor(rab[k]*dbng[k]);
        }
        csector = rab[0] + 2.0*rab[1] + 4.0*rab[2] + 1.0e-5;
        atmjk = &C->map[csector][C->nr[csector]+C->nr[8+csector]];
        atmjk->id = tmr->atomid;
        atmjk->q = tp->Charges[tmr->atomid];
        atmjk->lj = tp->LJIdx[tmr->atomid];
        for (k = 0; k < 3; k++) {
          atmjk->loc[k] = eploc[k];
          atmjk->frc[k] = 0.0;
        }
        C->nr[8+csector] += 1;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// LoadEPForExport: this function loads extra points (or, in principle, any
//                  atoms) from one of a cell's eight sectors into its export
//                  buffer.                                  
//                                                                      
// Arguments:                                                           
//   Aorig:        the original atom found in one of the cell's sectors 
//   Aexp:         pointer an element of the cell's export array        
//   dreg:         the destination region of the cell                   
//   reim[x,y,z]:  reimaging coordinate displacements                   
//-----------------------------------------------------------------------------
static void LoadEP4Export(atomc *Aorig, atomb *Aexp, int dreg, double* reimv)
{
  Aexp->id = Aorig->id;
  Aexp->dreg = dreg;
  Aexp->loc[0] = Aorig->loc[0] + reimv[0];
  Aexp->loc[1] = Aorig->loc[1] + reimv[1];
  Aexp->loc[2] = Aorig->loc[2] + reimv[2];
}

//-----------------------------------------------------------------------------
// PackageCellEP: this function wraps the four loops for sharing extra points
//                of one cell's four sectors with another cell.  
//                                                                      
// Arguments:                                                           
//   C:   the cell                                                      
//   i[C,D][1,2,3,4]: the sectors of cells O and N sharing extra points 
//   reimv:           the reimaging vector                              
//-----------------------------------------------------------------------------
static void PackageCellEP(cell* C, int iC1, int iC2, int iC3, int iC4,
                          int iD1, int iD2, int iD3, int iD4, double* reimv)
{
  int i, nexp;

  nexp = 1;
  for (i = C->nr[iC1]; i < C->nr[iC1] + C->nr[iC1+8]; i++) {
    LoadEP4Export(&C->map[iC1][i], &C->pexport[nexp], iD1, reimv);
    nexp++;
  }
  for (i = C->nr[iC2]; i < C->nr[iC2] + C->nr[iC2+8]; i++) {
    LoadEP4Export(&C->map[iC2][i], &C->pexport[nexp], iD2, reimv);
    nexp++;
  }
  for (i = C->nr[iC3]; i < C->nr[iC3] + C->nr[iC3+8]; i++) {
    LoadEP4Export(&C->map[iC3][i], &C->pexport[nexp], iD3, reimv);
    nexp++;
  }
  for (i = C->nr[iC4]; i < C->nr[iC4] + C->nr[iC4+8]; i++) {
    LoadEP4Export(&C->map[iC4][i], &C->pexport[nexp], iD4, reimv);
    nexp++;
  }
  C->pexport[0].id = nexp;
  C->nexp = nexp;
}

//-----------------------------------------------------------------------------
// UnpackCellEP: this function will unpack the imported extra points in a cell.
//                                                                      
// Arguments:                                                           
//   C:    the cell                                                     
//   tp:   the topology                                                 
//   fcn:  the face unit normal vectors                                 
//-----------------------------------------------------------------------------
static void UnpackCellEP(cell *C, prmtop *tp, dmat *fcn)
{
  int i, dreg, treg;
  int nadd[8];
  double fx[3], fy[3], fz[3], atmloc[3];
  atomb *tmpab;
  atomc *genatm;

  // Determine the number of atoms imported into each region
  SetIVec(nadd, 8, 0);
  for (i = 1; i < C->import[0].id; i++) {
    nadd[C->import[i].dreg] += 1;
  }
  for (i = 0; i < 3; i++) {
    fx[i] = fcn->map[0][i];
    fy[i] = fcn->map[1][i];
    fz[i] = fcn->map[2][i];
  }
  for (i = 1; i < C->import[0].id; i++) {
    tmpab = &C->import[i];
    atmloc[0] = tmpab->loc[0] - C->midp[0];
    atmloc[1] = tmpab->loc[1] - C->midp[1];
    atmloc[2] = tmpab->loc[2] - C->midp[2];
    dreg = tmpab->dreg;

    // Determine displacements from the planes of interest
    treg = (atmloc[0]*fx[0] + atmloc[1]*fx[1] + atmloc[2]*fx[2] >= 0.0) +
      (atmloc[0]*fy[0] + atmloc[1]*fy[1] + atmloc[2]*fy[2] >= 0.0)*2 +
      (atmloc[0]*fz[0] + atmloc[1]*fz[1] + atmloc[2]*fz[2] >= 0.0)*4;
    if (treg != dreg) {
      continue;
    }

    // If we're still here, this atom is a keeper
    genatm = &C->map[dreg][C->nr[dreg]+C->nr[dreg+8]];
    genatm->id = tmpab->id;
    genatm->loc[0] = tmpab->loc[0];
    genatm->loc[1] = tmpab->loc[1];
    genatm->loc[2] = tmpab->loc[2];
    genatm->frc[0] = 0.0;
    genatm->frc[1] = 0.0;
    genatm->frc[2] = 0.0;
    genatm->q = tp->Charges[tmpab->id];
    genatm->lj = tp->LJIdx[tmpab->id];
    C->nr[dreg+8] += 1;
  }
}

//-----------------------------------------------------------------------------
// UnpackCellGridEP: this routine handles unpacking of extra points from both
//                   cell local import buffers and cell grid pooled import
//                   buffers.                                    
//                                                                      
// Arguments:                                                           
//   cshr:      the communication plan                                  
//   CG:        the cell grid                                           
//   tp:        the topology                                            
//   fcn:       unit normals to cell faces                              
//-----------------------------------------------------------------------------
static void UnpackCellGridEP(ashr *cshr, cellgrid *CG, prmtop *tp, dmat *fcn)
{
  int i, j, nimp;
  atomb *timport;
  cell *N;

  // Unpack received, pooled messages
  for (i = 0; i < cshr->nrecv; i++) {
    nimp = 0;
    for (j = 0; j < cshr->recv[i].ncell; j++) {
      N = &CG->data[cshr->recv[i].cellpt[j]];
      timport = N->import;
      N->import = &CG->import[i][nimp];
      UnpackCellEP(N, tp, fcn);
      nimp += N->import[0].id;
      N->import = timport;
    }
  }

  // Unpack the "self receive"
  for (i = 0; i < cshr->selfrecv.ncell; i++) {
    UnpackCellEP(&CG->data[cshr->selfrecv.cellpt[i]], tp, fcn);
  }

#ifdef MPI
  // MPI implementations with multiple threads require a barrier 
  // here to avoid over-writing the contents of the import tables
  // before they can be transferred into the cells' atom tables  
  MPI_Barrier(CG->dspcomm);
#endif
}

//-----------------------------------------------------------------------------
// BroadcastSectors: broadcast extra points between individual sectors of
//                   cells.                                          
//                                                                      
// Arguments:                                                           
//   O,N:             the cells sharing extra points                    
//   i[O,N][1,2,3,4]: the sectors of cells O and N sharing extra points 
//   msign:           the sign of the move (+1 for moving forward, as in
//                    sharing atoms, -1 for reverse, as in merging forces)
//   reimv:           the reimaging vector for going O->N               
//   nmsg:            the number of the send message                    
//-----------------------------------------------------------------------------
static void BroadcastSectors(cell *O, cell *N, int iO1, int iO2, int iO3,
                             int iO4, int iN1, int iN2, int iN3, int iN4,
                             int msign, double* reimv, cellgrid *CG, int nmsg)
{
  atomb *texport;

  // Load the O cell's exports into the pooled export table
  if (O->CGRank != N->CGRank) {
    texport = O->pexport;
    O->pexport = &CG->pexport[nmsg][CG->nexp[nmsg]];
  }

  // Transfer directly to the N cell's import table
  else if (msign == 1) {
    texport = O->pexport;
    O->pexport = N->import;
  }

  // Add extra points from the requested
  // sectors to cell O's export table   
  PackageCellEP(O, iO1, iO2, iO3, iO4, iN1, iN2, iN3, iN4, reimv);

  // Return the O cell export pointer to its original location
  if (O->CGRank != N->CGRank) {
    CG->nexp[nmsg] += O->pexport[0].id;
    O->pexport = texport;
  }
  else if (msign == 1) {
    O->pexport = texport;
  }
}

//-----------------------------------------------------------------------------
// BroadcastExtraPoints: after extra points have been placed from within a cell
//                       framework, they must be broadcast to all cells.  This
//                       is done by having each cell make a list, for each of
//                       its regions, detailing the extra points for which it
//                       knows the locations.  The lists begin when a cell
//                       places extra points for atoms within its central
//                       region, the patch of territory that is at least 0.5
//                       cell lengths from any cell face.  The extra points can
//                       end up in any sector of the cell, however.  Because of
//                       this, a cell's list of extra points spans all eight
//                       sectors.                             
//                                                                      
// Arguments:                                                           
//   cellid:    the ID number of the cell in the grid                   
//   imove:     the axis of the move (0 = X, 1 = Y, 2 = Z)              
//   msign:     the direction of the move along the axis (1 moves as    
//              coordinates would be shared amongst cells, -1 moves as forces
//              would be merged among cells)                     
//   CG:        the cell grid                                           
//   crd:       the coordinates (for grid dimensions and orthonormality)
//   tp:        the topology                                            
//-----------------------------------------------------------------------------
static void BroadcastExtraPoints(int cellid, int imove, int msign, int nmsg,
                                 cellgrid *CG, coord *crd)
{
  int i;
  int nL[3];
  double reimv[3];
  cell *O, *N;

  // Origin and destination cells, reimaging considerations
  O = &CG->data[cellid];
  for (i = 0; i < 3; i++) {
    reimv[i] = 0.0;
    nL[i] = O->gbin[i];
  }
  nL[imove] = (msign == 1) ? (O->gbin[imove] - 1) : (O->gbin[imove] + 1);
  if (nL[imove] == -1) {
    nL[imove] = CG->ng[imove] - 1;
    if (crd->isortho == 1) {
      reimv[imove] = crd->gdim[imove];
    }
    else {
      reimv[imove] = 1.0;
      RotateCrd(reimv, 1, crd->invU);
    }
  }
  else if (nL[imove] == CG->ng[imove]) {
    nL[imove] = 0;
    if (crd->isortho == 1) {
      reimv[imove] = -crd->gdim[imove];
    }
    else {
      reimv[imove] = -1.0;
      RotateCrd(reimv, 1, crd->invU);
    }
  }
  N = &CG->map[nL[0]][nL[1]][nL[2]];

  // X-direction broadcast
  if (imove == 0) {
    if (msign == 1) {
      BroadcastSectors(O, N, 0, 2, 4, 6, 1, 3, 5, 7, msign, reimv, CG, nmsg);
    }
    else {
      BroadcastSectors(O, N, 1, 3, 5, 7, 0, 2, 4, 6, msign, reimv, CG, nmsg);
    }
  }

  // Y-direction broadcast
  else if (imove == 1) {
    if (msign == 1) {
      BroadcastSectors(O, N, 0, 1, 4, 5, 2, 3, 6, 7, msign, reimv, CG, nmsg);
    }
    else {
      BroadcastSectors(O, N, 2, 3, 6, 7, 0, 1, 4, 5, msign, reimv, CG, nmsg);
    }
  }

  // Z-direction broadcast
  else if (imove == 2) {
    if (msign == 1) {
      BroadcastSectors(O, N, 0, 1, 2, 3, 4, 5, 6, 7, msign, reimv, CG, nmsg);
    }
    else {
      BroadcastSectors(O, N, 4, 5, 6, 7, 0, 1, 2, 3, msign, reimv, CG, nmsg);
    }
  }
}

//-----------------------------------------------------------------------------
// CellMergeSort: merge the two parts of the atomc data array in the primary
//                sector of a cell.  This function implements a quicksort to
//                sort the elements of the second part of the array (the extra
//                points, added in eight blocks and therefore not in much
//                discernable order) and then uses a merge sort to combine the
//                two parts of the array as the first part is already sorted.
//                                                                      
// Arguments:                                                           
//   C:       the cell                                                  
//-----------------------------------------------------------------------------
static void CellMergeSort(cell *C)
{
  int i, atmcon, epcon;
  atomc *tatm;

  // Bail out if there are no extra points to merge
  if (C->nr[8] == 0) {
    return;
  }

  // Sort the second part of the array
  qsort(&C->data[C->nr[0]], C->nr[8], sizeof(atomc), SortAtomID);

  // Set counters.  The strategy is to start counters at the     
  // beginning of the list of atoms (atmcom) and at the beginning
  // of the list of extra points (epcon).  If the extra point at 
  // epcon has a lower atom ID number that the atom at atmcon,   
  // then the extra point replaces the atom and the atom gets put
  // in a holding array.                                         
  atmcon = C->nr[0]-1;
  const int natm = C->nr[0] + C->nr[8];
  const int cnr0 = C->nr[0];
  epcon = natm-1;
  tatm = C->data;
  for (i = natm-1; i >= 0; i--) {
    if (atmcon < 0 || (epcon >= cnr0 && tatm[epcon].id > tatm[atmcon].id)) {
      C->atmscr[i] = tatm[epcon];
      epcon--;
    }
    else {
      C->atmscr[i] = tatm[atmcon];
      atmcon--;
    }
  }
  for (i = 0; i < natm; i++) {
    C->data[i] = C->atmscr[i];
  }

  // The extra points are now merged into the primary sector
  C->nr[0] += C->nr[8];
}

//-----------------------------------------------------------------------------
// RemoveXpt: this function will loop over all cells and remove any extra
//            points in each sector.  This is done to avoid the duplication
//            and buildup of extra points, which would otherwise occur in the
//            course of normal dynamics as extra points are not subject to the
//            usual atom migration rules. 
//                                                                      
// Arguments:                                                           
//   CG:      the cell grid                                             
//   tp:      the system topology
//-----------------------------------------------------------------------------
static void RemoveXpt(cellgrid *CG, prmtop *tp)
{
  int i, j, k, nmassive;
  double *mass;
  cell *C;
  atomc *catm;

  mass = tp->Masses;
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < 8; j++) {
      catm = C->map[j];
      nmassive = 0;
      for (k = 0; k < C->nr[j]; k++) {
        if (mass[catm[k].id] > 1.0e-8) {
          if (nmassive < k) {
            catm[nmassive] = catm[k];
          }
          nmassive++;
        }
      }
      C->nr[j] = nmassive;
    }
  }
}

//-----------------------------------------------------------------------------
// ExtraPointLocations: this function wraps the complex series of calls to
//                      CellPlaceXpt and BroadcastExtraPoints in order to
//                      manage the placement of extra points in the context of
//                      a cell grid.                     
//                                                                      
// Arguments:                                                           
//   tp:         the topology                                           
//   crd:        the coordinates                                        
//   CG:         the cell grid                                          
//   tj:         in MPI implementation, defines MPI custom types        
//-----------------------------------------------------------------------------
#ifdef MPI
void ExtraPointLocations(prmtop *tp, coord *crd, cellgrid *CG, trajcon *tj)
#else
void ExtraPointLocations(prmtop *tp, coord *crd, cellgrid *CG)
#endif
{
  int i, j, k, imove;
  cell *O, *N;
  ashr *cshr, *cshr2;

#ifdef MPI
  MPI_Request* req;
  MPI_Request* req2;
  MPI_Status* stt;
  MPI_Status* stt2;
  req = (MPI_Request*)malloc((CG->nsend + CG->nrecv)*sizeof(MPI_Request));
  stt = (MPI_Status*)malloc((CG->nsend + CG->nrecv)*sizeof(MPI_Status));
  req2 = (MPI_Request*)malloc((CG->nsend + CG->nrecv)*sizeof(MPI_Request));
  stt2 = (MPI_Status*)malloc((CG->nsend + CG->nrecv)*sizeof(MPI_Status));
#endif

  // Bail right out if there are no extra points
  if (tp->nxtrapt == 0) {
    return;
  }

  // Make sure that there are no extra points 
  // currently in the cell grid, and also that
  // all cells which this process does not own
  // have zero atom count.                    
  RemoveXpt(CG, tp);

  // Update the atom GPS in preparation
  // for lots of referencing.          
  UpdateCellGPS(CG);

  // Set the extra point counters in each cell to zero.  
  // Then, loop over all cells and place each extra point
  // as indicated by the location of its frame atom 1.   
  for (i = 0; i < CG->MyCellCount; i++) {
    CellPlaceXpt(tp, crd, &CG->data[CG->MyCellDomain[i]], CG);
  }

  // Each cell might now have some extra points in its 
  // various sectors.  While they are still in easily  
  // accessible lists of known length at the back of   
  // each sector's atom list, broadcast them to sectors
  // of other cells which need to know their locations.
  // First, broadcast extra points in sectors 1, 3, 5, 
  // and 7 of cell (i,j,k) to sectors 0, 2, 4, and 6 of
  // cell (i+1,j,k), and vice-versa.  Then, broadcast  
  // extra points in sectors 2, 3, 6, and 7 of cell    
  // (i,j,k) to sectors 0, 1, 4, and 5 of cell         
  // (i,j+1,k).  Finally, broadcast extra points in    
  // sectors 5, 6, 7, and 8 of cell (i,j,k) to sectors 
  // 0, 1, 2, and 3 of cell (i,j,k+1).  The process is 
  // more complex than other transfers, as the pair of 
  // cells must swap their coordinates and there is no 
  // simple position metric by which to judge which    
  // extra points were originally placed by particular 
  // cells.  Cells therefore transfer to one another's 
  // import buffers, then load up their own export     
  // buffers for the second pulse, unpack their import 
  // buffers, and communicate the second pulse.        
  for (imove = 0; imove < 3; imove++) {
    cshr = &CG->DirCommPlan.mvshr[imove];

#ifdef MPI
    // Post receives for the first pulse
    for (i = 0; i < cshr->nrecv; i++) {
      MPI_Irecv(CG->import[i], CG->maximp[i], tj->MPI_ATOMB,
                cshr->recv[i].partner, cshr->recv[i].BaseID + DSP_EXP_LOC1,
                CG->dspcomm, &req[i]);
    }
    int nreq = cshr->nrecv;
#endif

    // Post sends for the first pulse
    for (i = 0; i < cshr->nsend; i++) {
      CG->nexp[i] = 0;
      for (j = 0; j < cshr->send[i].ncell; j++) {
        BroadcastExtraPoints(cshr->send[i].cellpt[j], imove, 1, i, CG, crd);
      }
#ifdef MPI
      if (cshr->send[i].partner != CG->tid) {
        MPI_Isend(CG->pexport[i], CG->nexp[i], tj->MPI_ATOMB,
                  cshr->send[i].partner, cshr->send[i].BaseID + DSP_EXP_LOC1,
                  CG->dspcomm, &req[nreq]);
        nreq++;
      }
#endif
    }

    // Set the cell:cell communication plan
    // in preparation for the second pulse 
    cshr2 = &CG->DirCommPlan.frcmg[imove];

#ifdef MPI
    // Wait for messages to complete in the first pulse
    MPI_Waitall(nreq, req, stt);

    // Post receives for the second pulse
    for (i = 0; i < cshr2->nrecv; i++) {
      MPI_Irecv(CG->import[i], CG->maximp[i], tj->MPI_ATOMB,
                cshr2->recv[i].partner, cshr2->recv[i].BaseID + DSP_EXP_LOC2,
                CG->dspcomm, &req2[i]);
    }
    int nreq2 = cshr2->nrecv;
#endif

    // Prepare sends for the second pulse
    for (i = 0; i < cshr2->nsend; i++) {
      CG->nexp[i] = 0;
      for (j = 0; j < cshr2->send[i].ncell; j++) {
        BroadcastExtraPoints(cshr2->send[i].cellpt[j], imove, -1, i, CG, crd);
      }
    }

    // Unpack from the first pulse
    UnpackCellGridEP(cshr, CG, tp, &crd->fcnorm);

    // Post sends from the second pulse
    for (i = 0; i < cshr2->nsend; i++) {
#ifdef MPI
      if (cshr2->send[i].partner != CG->tid) {
        MPI_Isend(CG->pexport[i], CG->nexp[i], tj->MPI_ATOMB,
                  cshr2->send[i].partner, cshr2->send[i].BaseID + DSP_EXP_LOC2,
                  CG->dspcomm, &req2[nreq2]);
        nreq2++;
      }
      else {
#endif
        for (j = 0; j < cshr2->send[i].ncell; j++) {
          O = &CG->data[cshr2->send[i].cellpt[j]];
          N = &CG->data[cshr2->selfrecv.cellpt[j]];
          for (k = 0; k < O->pexport[0].id; k++) {
            N->import[k] = O->pexport[k];
          }
        }
#ifdef MPI
      }
#endif
    }
#ifdef MPI
    // Wait for messages to complete in second pulse
    MPI_Waitall(nreq2, req2, stt2);
#endif

    // Unpack from the second pulse
    UnpackCellGridEP(cshr2, CG, tp, &crd->fcnorm);
  }

  // Now, all cells have the right extra points in them,
  // but the increasing order of atom IDs within all    
  // sectors of all cells must be maintained.  There are
  // essentially two lists in each sector of each cell  
  // with atoms arranged in order of increasing atom ID 
  // number: atoms with mass, and atoms without mass    
  // (extra points).                                    
  for (i = 0; i < CG->MyCellCount; i++) {
    N = &CG->data[CG->MyCellDomain[i]];
    CellMergeSort(N);
    for (j = 1; j < 8; j++) {
      qsort(N->map[j], N->nr[j] + N->nr[j+8], sizeof(atomc), SortAtomID);
      N->nr[j] += N->nr[j+8];
    }
  }

#ifdef MPI
  // Free allocated memory
  free(req);
  free(req2);
  free(stt);
  free(stt2);
#endif
}

//-----------------------------------------------------------------------------
// XptForceXfer: transfer the force on an extra point to as many as four frame
//               atoms.
//
// Arguments:
//   [a,b,c,d,ep]loc:    locations of frame atoms and the extra point
//   [a,b,c,d,ep]frc:    forces on the frame atoms and the extra point
//-----------------------------------------------------------------------------
void XptForceXfer(double *aloc, double *bloc, double *cloc, double *dloc,
		  double *eploc, double *afrc, double *bfrc, double *cfrc,
		  double *dfrc, double *epfrc, expt *tmr)
{
  int i;
  double dvx, dvy, dvz, g01, g02, g12, gamma, pbfac, magdvec, dotdu, dotdv;
  double magab, magac, magPr, magab2, magPr2, invab, invPr, invab2, invPr2;
  double F1, F2, F3, f1fac, f2fac, nf1fac, nf2fac, abbcOabab, magup, magvp;
  double dotuv, sfac, uvdis, vudis, dphidu, dphidv, dphidw, fproj, magdvec2;
  double cfx, cfy, cfz;
  double rab[3], rac[3], rad[3], rbc[3], fb[3], fc[3], fd[3];
  double raEP[3], rPerp[3];
  double pbold[3], rja[3], rjb[3], rjab[3], rm[3], rt[3], dr[3], up[3], vp[3];
  
  // Compute displacements for higher-complexity frames
  if (tmr->frstyle > 2) {
    if (tmr->frstyle == 9 || tmr->frstyle == 10) {
      for (i = 0; i < 3; i++) {
        rab[i] = 0.5*(bloc[i] + cloc[i]) - aloc[i];
        rac[i] = 0.5*(bloc[i] + dloc[i]) - aloc[i];
      }
    }
    else {
      for (i = 0; i < 3; i++) {
        rab[i] = bloc[i] - aloc[i];
        rac[i] = cloc[i] - aloc[i];
        rbc[i] = cloc[i] - bloc[i];
      }
    }
  }

  // Transfer forces
  if (tmr->frstyle == 1) {
    for (i = 0; i < 3; i++) {
      afrc[i] += (1.0-tmr->d1) * epfrc[i];
      bfrc[i] += tmr->d1 * epfrc[i];
      epfrc[i] = 0.0;
    }
  }
  else if (tmr->frstyle == 11) {
    invab = 1.0 / sqrt(rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2]);
    fproj = 0.0;
    for (i = 0; i < 3; i++) {
      raEP[i] = eploc[i] - aloc[i];
      fproj = fproj + raEP[i]*epfrc[i];
    }
    fproj /= raEP[0]*raEP[0] + raEP[1]*raEP[1] + raEP[2]*raEP[2];
    for (i = 0; i < 3; i++) {
      fb[i] = invab * tmr->d1 * (epfrc[i] - fproj*raEP[i]);
      afrc[i] += epfrc[i] - fb[i];
      bfrc[i] += fb[i];
      epfrc[i] = 0.0;
    }
  }
  else if (tmr->frstyle == 2) {
    for (i = 0; i < 3; i++) {
      afrc[i] += (1.0-tmr->d1-tmr->d2) * epfrc[i];
      bfrc[i] += tmr->d1 * epfrc[i];
      cfrc[i] += tmr->d2 * epfrc[i];
      epfrc[i] = 0.0;
    }
  }
  else if (tmr->frstyle == 3) {
    dvx = rab[0] + tmr->d3*rbc[0];
    dvy = rab[1] + tmr->d3*rbc[1];
    dvz = rab[2] + tmr->d3*rbc[2];
    gamma = tmr->d1/sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
    raEP[0] = eploc[0] - aloc[0];
    raEP[1] = eploc[1] - aloc[1];
    raEP[2] = eploc[2] - aloc[2];
    pbfac = Dot3(raEP, epfrc) / Dot3(raEP, raEP);
    pbold[0] = epfrc[0] - pbfac*raEP[0];
    pbold[1] = epfrc[1] - pbfac*raEP[1];
    pbold[2] = epfrc[2] - pbfac*raEP[2];
    for (i = 0; i < 3; i++) {
      afrc[i] += epfrc[i]      - gamma * pbold[i];
      bfrc[i] += (1.0-tmr->d3) * gamma * pbold[i];
      cfrc[i] +=      tmr->d3  * gamma * pbold[i];
      epfrc[i] = 0.0;
    }
  }
  else if (tmr->frstyle == 4) {
    magdvec = Dot3(rab, rbc) / Dot3(rab, rab);
    rPerp[0] = rbc[0] - magdvec*rab[0];
    rPerp[1] = rbc[1] - magdvec*rab[1];
    rPerp[2] = rbc[2] - magdvec*rab[2];
    magab2 = Dot3(rab, rab);
    magab = sqrt(magab2);
    invab = 1.0/magab;
    invab2 = invab*invab;
    magPr2 = Dot3(rPerp, rPerp);
    magPr = sqrt(magPr2);
    invPr = 1.0/magPr;
    invPr2 = invPr*invPr;
    f1fac = Dot3(rab, epfrc)*invab2;
    f2fac = Dot3(rPerp, epfrc)*invPr2;
    nf1fac = tmr->d1*cos(tmr->d2)*invab;
    nf2fac = tmr->d1*sin(tmr->d2)*invPr;
    abbcOabab = Dot3(rab, rbc)*invab2;
    for (i = 0; i < 3; i++) {
      F1 = epfrc[i] - f1fac*rab[i];
      F2 = F1 - f2fac*rPerp[i];
      F3 = f1fac*rPerp[i];
      afrc[i] += epfrc[i]  - nf1fac*F1 + nf2fac*(abbcOabab*F2 + F3);
      bfrc[i] += nf1fac*F1 - nf2fac*(F2 + abbcOabab*F2 + F3);
      cfrc[i] += nf2fac*F2;
      epfrc[i] = 0.0;
    }
  }
  else if (tmr->frstyle == 5) {
    g01 = tmr->d3*rac[2];
    g02 = tmr->d3*rac[1];
    g12 = tmr->d3*rac[0];
    fb[0] = tmr->d1*epfrc[0] -     g01*epfrc[1] +     g02*epfrc[2];
    fb[1] =     g01*epfrc[0] + tmr->d1*epfrc[1] -     g12*epfrc[2];
    fb[2] =    -g02*epfrc[0] +     g12*epfrc[1] + tmr->d1*epfrc[2];
    g01 = tmr->d3*rab[2];
    g02 = tmr->d3*rab[1];
    g12 = tmr->d3*rab[0];
    fc[0] = tmr->d2*epfrc[0] +     g01*epfrc[1] -     g02*epfrc[2];
    fc[1] =    -g01*epfrc[0] + tmr->d2*epfrc[1] +     g12*epfrc[2];
    fc[2] =     g02*epfrc[0] -     g12*epfrc[1] + tmr->d2*epfrc[2];
    for (i = 0; i < 3; i++) {
      afrc[i] += epfrc[i] - fb[i] - fc[i];
      bfrc[i] += fb[i];
      cfrc[i] += fc[i];
      epfrc[i] = 0.0;
    }
  }
  else if (tmr->frstyle == 6) {
    for (i = 0; i < 3; i++) {
      rad[i]  = dloc[i] - aloc[i];
      rja[i]  = (tmr->d1 * rac[i]) - rab[i];
      rjb[i]  = (tmr->d2 * rad[i]) - rab[i];
      rjab[i] = rjb[i] - rja[i];
    }
    CrossP(rja, rjb, rm);
    magdvec = 1.0 / sqrt(rm[0]*rm[0] + rm[1]*rm[1] + rm[2]*rm[2]);
    magdvec2 = magdvec * magdvec;
    cfx = tmr->d3 * magdvec * epfrc[0];
    cfy = tmr->d3 * magdvec * epfrc[1];
    cfz = tmr->d3 * magdvec * epfrc[2];
    CrossP(rm, rjab, rt);
    rt[0] = rt[0] * magdvec2;
    rt[1] = rt[1] * magdvec2;
    rt[2] = rt[2] * magdvec2;
    fb[0] = -rm[0]*rt[0]*cfx + (rjab[2] - rm[1]*rt[0])*cfy -
            (rjab[1] + rm[2]*rt[0])*cfz;
    fb[1] = -(rjab[2] + rm[0]*rt[1])*cfx - rm[1]*rt[1]*cfy +
            (rjab[0] - rm[2]*rt[1])*cfz;
    fb[2] = (rjab[1] - rm[0]*rt[2])*cfx - (rjab[0] + rm[1]*rt[2])*cfy -
            rm[2]*rt[2]*cfz;
    rt[0] = (rjb[1]*rm[2] - rjb[2]*rm[1]) * magdvec2 * tmr->d1;
    rt[1] = (rjb[2]*rm[0] - rjb[0]*rm[2]) * magdvec2 * tmr->d1;
    rt[2] = (rjb[0]*rm[1] - rjb[1]*rm[0]) * magdvec2 * tmr->d1;
    fc[0] = -rm[0]*rt[0]*cfx - (tmr->d1*rjb[2] + rm[1]*rt[0])*cfy +
            (tmr->d1*rjb[1] - rm[2]*rt[0])*cfz;
    fc[1] = (tmr->d1*rjb[2] - rm[0]*rt[1])*cfx - rm[1]*rt[1]*cfy -
            (tmr->d1*rjb[0] + rm[2]*rt[1])*cfz;
    fc[2] = -(tmr->d1*rjb[1] + rm[0]*rt[2])*cfx +
            (tmr->d1*rjb[0] - rm[1]*rt[2])*cfy - rm[2]*rt[2]*cfz;
    rt[0] = (rm[1]*rja[2] - rm[2]*rja[1]) * magdvec2 * tmr->d2;
    rt[1] = (rm[2]*rja[0] - rm[0]*rja[2]) * magdvec2 * tmr->d2;
    rt[2] = (rm[0]*rja[1] - rm[1]*rja[0]) * magdvec2 * tmr->d2;
    fd[0] = -rm[0]*rt[0]*cfx + (tmr->d2*rja[2] - rm[1]*rt[0])*cfy -
            (tmr->d2*rja[1] + rm[2]*rt[0])*cfz;
    fd[1] = -(tmr->d2*rja[2] + rm[0]*rt[1])*cfx - rm[1]*rt[1]*cfy +
            (tmr->d2*rja[0] - rm[2]*rt[1])*cfz;
    fd[2] = (tmr->d2*rja[1] - rm[0]*rt[2])*cfx +
            -(tmr->d2*rja[0] + rm[1]*rt[2])*cfy - rm[2]*rt[2]*cfz;
    for (i = 0; i < 3; i++) {
      afrc[i] += epfrc[i] - fb[i] - fc[i] - fd[i];
      bfrc[i] += fb[i];
      cfrc[i] += fc[i];
      dfrc[i] += fd[i];
      epfrc[i] = 0.0;
    }
  }
  else if (tmr->frstyle >= 7 && tmr->frstyle <= 10) {
    double rx = eploc[0] - aloc[0];
    double ry = eploc[1] - aloc[1];
    double rz = eploc[2] - aloc[2];
    double torqueX = ry*epfrc[2] - rz*epfrc[1];
    double torqueY = rz*epfrc[0] - rx*epfrc[2];
    double torqueZ = rx*epfrc[1] - ry*epfrc[0];
    magab = 1.0 / sqrt(rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2]);
    magac = 1.0 / sqrt(rac[0]*rac[0] + rac[1]*rac[1] + rac[2]*rac[2]);
    for (i = 0; i < 3; i++) {
      rab[i] *= magab;
      rac[i] *= magac;
    }
    CrossP(rab, rac, rPerp);
    magdvec = 1.0 / sqrt(rPerp[0]*rPerp[0] + rPerp[1]*rPerp[1] +
                         rPerp[2]*rPerp[2]);
    rPerp[0] *= magdvec;
    rPerp[1] *= magdvec;
    rPerp[2] *= magdvec;
    dr[0] = rac[0] - rab[0];
    dr[1] = rac[1] - rab[1];
    dr[2] = rac[2] - rab[2];
    dotdu = rab[0]*dr[0] + rab[1]*dr[1] + rab[2]*dr[2];
    dotdv = rac[0]*dr[0] + rac[1]*dr[1] + rac[2]*dr[2];
    for (i = 0; i < 3; i++) {
      up[i] = dr[i] - dotdu*rab[i];
      vp[i] = dr[i] - dotdv*rac[i];
    }
    magup = 1.0 / sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
    magvp = 1.0 / sqrt(vp[0]*vp[0] + vp[1]*vp[1] + vp[2]*vp[2]);
    for (i = 0; i < 3; i++) {
      up[i] *= magup;
      vp[i] *= magvp;
    }
    dotuv = rab[0]*rac[0] + rab[1]*rac[1] + rab[2]*rac[2];
    sfac = 1.0 / sqrt(1.0 - dotuv*dotuv);
    uvdis = magab * sfac;
    vudis = magac * sfac;
    dphidu = -(torqueX*rab[0]   + torqueY*rab[1]   + torqueZ*rab[2]);
    dphidv = -(torqueX*rac[0]   + torqueY*rac[1]   + torqueZ*rac[2]);
    dphidw = -(torqueX*rPerp[0] + torqueY*rPerp[1] + torqueZ*rPerp[2]);
    dphidv *= uvdis;
    dphidu *= vudis;
    magab *= 0.5 * dphidw;
    magac *= 0.5 * dphidw;
    for (i = 0; i < 3; i++) {
      fb[i] = -rPerp[i]*dphidv + up[i]*magab;
      fc[i] =  rPerp[i]*dphidu + vp[i]*magac;
      if (tmr->frstyle < 9) {
        bfrc[i] -= fb[i];
        cfrc[i] -= fc[i];
      }
      else {
        bfrc[i] -= 0.5 * (fb[i] + fc[i]);
        cfrc[i] -= 0.5 * (fb[i]);
        dfrc[i] -= 0.5 * (fc[i]);
      }
      afrc[i] += fb[i] + fc[i] + epfrc[i];
      epfrc[i] = 0.0;
    }
  }
  else {
    printf("XptForceXfer >> Error.  Not prepared to handle frame "
           "type %d.\n", tmr->frstyle);
    exit(1);
  }
}

//-----------------------------------------------------------------------------
// CellXferEPForces: transfer forces from extra points within a cell to their
//                   frame atoms.  Only forces on extra points in the cell's
//                   central region will be transferred, but the transfered
//                   forces must then be propagated over neighboring cells. 
//                   This function is a counterpart of CellPlaceXpt.
//                                                                      
// Arguments:                                                           
//   C:      the cell of interest                                       
//   tp:     the topology structure                                     
//   crd:    the system coordinates                                     
//-----------------------------------------------------------------------------
void CellXferEPForces(cell *C, prmtop *tp, coord *crd, cellgrid *CG)
{
  int h, i, j, k, epid, bid, cid, did;
  int epreg, breg, creg, dreg, epelem, belem, celem, delem;
  double *aloc, *bloc, *cloc, *dloc, *eploc;
  double *afrc, *bfrc, *cfrc, *dfrc, *epfrc, *ccen, *dbng;
  expt *tmr;
  atomc *atmhi;
  lgrp *hiFR1;

  // Pointers
  ccen = &CG->celldim[4];
  dbng = CG->dbng;

  // Loop over all sectors
  for (h = 0; h < 8; h++) {
    for (i = 0; i < C->nr[h]; i++) {
      atmhi = &C->map[h][i];
      hiFR1 = &tp->FR1Idx[atmhi->id];
      if (hiFR1->natom == 0 ||
          IsCentralAtom(C->orig, ccen, atmhi->loc,
                        crd->U.data, dbng, crd->isortho) < 0) {
        continue;
      }

      // This atom is frame atom 1 to one or 
      // more extra points; place them one by
      // one and then add them to the proper 
      // regions of this cell.               
      for (j = 0; j < hiFR1->natom; j++) {

        // Pointer to extra point structure
        tmr = &tp->xtrapts[hiFR1->atoms[j]];

        // Pointers to atom locations
        aloc = atmhi->loc;
        afrc = atmhi->frc;
        bid = C->GPSptr[tmr->fr2];
        breg = bid/C->maxatom;
        belem = bid - breg*C->maxatom;
        bloc = C->map[breg][belem].loc;
        bfrc = C->map[breg][belem].frc;
        epid = C->GPSptr[tmr->atomid];
        epreg = epid/C->maxatom;
        epelem = epid - epreg*C->maxatom;
        eploc = C->map[epreg][epelem].loc;
        epfrc = C->map[epreg][epelem].frc;
        if (tmr->frstyle != 1) {
          cid = C->GPSptr[tmr->fr3];
          creg = cid/C->maxatom;
          celem = cid - creg*C->maxatom;
          cloc = C->map[creg][celem].loc;
          cfrc = C->map[creg][celem].frc;
        }
        if (tmr->frstyle == 6 || tmr->frstyle == 9 || tmr->frstyle == 10) {
          did = C->GPSptr[tmr->fr4];
          dreg = did/C->maxatom;
          delem = did - dreg*C->maxatom;
          dloc = C->map[dreg][delem].loc;
          dfrc = C->map[dreg][delem].frc;
        }

        // Compute displacements for higher-complexity frames
        XptForceXfer(aloc, bloc, cloc, dloc, eploc, afrc, bfrc, cfrc, dfrc,
                     epfrc, tmr);
      }
    }
  }
}
