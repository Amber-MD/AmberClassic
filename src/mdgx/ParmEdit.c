#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ParmEdit.h"
#include "Topology.h"
#include "Trajectory.h"
#include "VirtualSites.h"
#include "Parse.h"
#include "Manual.h"
#include "mdgxVector.h"

//-----------------------------------------------------------------------------
// ReportExtraPoints: append a report of extra points, added and inherited, in
//                    the recently printed topology.
//
// Arguments:
//   tp:       the (now edited) topology
//   outp:     output file where report shall be written
//-----------------------------------------------------------------------------
static void ReportExtraPoints(prmtop *tp, FILE *outp)
{
  int i, j, found, nunique, nbpidx, ngiven;
  int epOrigPop[12];
  int epInsPop[12];
  int* epDescribed;
  int* examples;
  double tsig, teps;
  char exword[16], spcword[16];
  expt thispt;

  if (tp->EPInserted == 1) {
    fprintf(outp, "\n Extra points were added to this topology.\n");
  }
  else {
    fprintf(outp, "\n Extra points were inherited from the original topology."
	    "\n");
  }
  SetIVec(epOrigPop, 12, 0);
  SetIVec(epInsPop, 12, 0);
  for (i = 0; i < tp->nxtrapt; i++) {
    if (tp->EPInserted) {
      if (tp->OldAtomNum[tp->xtrapts[i].atomid] == -1) {
        epInsPop[tp->xtrapts[i].frstyle] += 1;
      }
      else {
        epOrigPop[tp->xtrapts[i].frstyle] += 1;
      }
    }
    else {
      epOrigPop[tp->xtrapts[i].frstyle] += 1;
    }
  }
  fprintf(outp, "\n     Frame Type    Type Alias    Inherited     Added\n");
  fprintf(outp, "     ----------   ------------   ---------   ---------\n");
  for (i = 1; i < 12; i++) {
    if (epOrigPop[i] > 0 || epInsPop[i] > 0) {
      if (i == 1) {
	fprintf(outp, "       Flex-2       FlexDis2  ");
      }
      else if (i == 11) {
	fprintf(outp, "        FD-2        FixedDis2 ");
      }
      else if (i == 2) {
	fprintf(outp, "       Flex-3       FlexDis3  ");
      }
      else if (i == 3) {
	fprintf(outp, "        FD-3        FixedDis3 ");
      }
      else if (i == 4) {
	fprintf(outp, "       FAD-3       FixAnglDis ");
      }
      else if (i == 5) {
	fprintf(outp, "       Out-3       OutOfPlane ");
      }
      else if (i == 6) {
	fprintf(outp, "        FD-4        FourPoint ");
      }
      else if (i == 7) {
	fprintf(outp, "     Inf-Bis-1    InferredBis1");
      }
      else if (i == 8) {
	fprintf(outp, "     Inf-Bis-2    InferredBis2");
      }
      else if (i == 9) {
	fprintf(outp, "     Inf-InP-1    InferredInP1");
      }
      else if (i == 10) {
	fprintf(outp, "     Inf-InP-2    InferredInP2");
      }
      fprintf(outp, "    %7d     %7d\n", epOrigPop[i], epInsPop[i]);
    }
  }
  fprintf(outp, "\n Details of added extra points follow.\n"
	  "   Example format:\n"
	  "     [Extra Point] <Residue> <Res. No.> -- <Atom Name> <Atom No.>\n"
	  "     [Parent Atom] <Residue> <Res. No.> -- <Atom Name> <Atom No.>\n"
	  "\n");
  epDescribed = (int*)calloc(tp->nxtrapt, sizeof(int));
  examples = (int*)calloc(tp->nxtrapt, sizeof(int));
  nunique = 0;
  while (ISum(epDescribed, tp->nxtrapt) < tp->nxtrapt) {
    i = 0;
    while (epDescribed[i] == 1) {
      i++;
    }
    thispt = tp->xtrapts[i];
    SetIVec(examples, tp->nxtrapt, 0);
    examples[i] = 1;
    epDescribed[i] = 1;
    for (j = i+1; j < tp->nxtrapt; j++) {
      if (tp->xtrapts[j].frstyle == thispt.frstyle &&
	  fabs(tp->Charges[tp->xtrapts[j].atomid] -
               tp->Charges[thispt.atomid]) < 1.0e-6 &&
	  fabs(tp->xtrapts[j].d1 - thispt.d1) < 1.0e-6 &&
	  fabs(tp->xtrapts[j].d2 - thispt.d2) < 1.0e-6 &&
	  fabs(tp->xtrapts[j].d3 - thispt.d3) < 1.0e-6 &&
          tp->LJIdx[tp->xtrapts[j].atomid] == tp->LJIdx[thispt.atomid]) {
	examples[j] = 1;
        epDescribed[j] = 1;
      }
    }
    fprintf(outp, " Virtual Site ID %4d (%4d instances):\n", nunique,
	    ISum(examples, tp->nxtrapt));
    if (thispt.frstyle == 1) {
      fprintf(outp, "   Frame type:   %s / %s\n", "Flex-2",
	      "FlexDis2");
    }
    else if (thispt.frstyle == 11) {
      fprintf(outp, "   Frame type:   %s / %s\n", "FD-2",
	      "FixedDis2");
    }
    else if (thispt.frstyle == 2) {
      fprintf(outp, "   Frame type:   %s / %s\n", "Flex-3",
	      "FlexDis3");
    }
    else if (thispt.frstyle == 3) {
      fprintf(outp, "   Frame type:   %s / %s\n", "FD-3",
	      "FixedDis3");
    }
    else if (thispt.frstyle == 4) {
      fprintf(outp, "   Frame type:   %s / %s\n", "FAD-3",
	      "FixAnglDis");
    }
    else if (thispt.frstyle == 5) {
      fprintf(outp, "   Frame type:   %s / %s\n", "Out-3",
	      "OutOfPlane");
    }
    else if (thispt.frstyle == 6) {
      fprintf(outp, "   Frame type:   %s / %s\n", "FD-4",
	      "FourPoint");
    }
    else if (thispt.frstyle == 7) {
      fprintf(outp, "   Frame type:   %s / %s\n", "Inf-Bis-1",
	      "InferredBis1");
    }
    else if (thispt.frstyle == 8) {
      fprintf(outp, "   Frame type:   %s / %s\n", "Inf-Bis-2",
	      "InferredBis2");
    }
    else if (thispt.frstyle == 9) {
      fprintf(outp, "   Frame type:   %s / %s\n", "Inf-InP-1",
	      "InferredInP1");
    }
    else if (thispt.frstyle == 10) {
      fprintf(outp, "   Frame type:   %s / %s\n", "Inf-InP-2",
	      "InferredInP2");
    }
    if (tp->LJIdx[thispt.atomid] >= 0) {
      nbpidx = tp->NBParmIdx[tp->LJIdx[thispt.atomid] * (tp->ntypes + 1)];
      if (tp->LJA[nbpidx] > 1.0e-8 && tp->LJB[nbpidx] > 1.0e-8) {
        tsig = pow(tp->LJA[nbpidx] / tp->LJB[nbpidx], 1.0 / 6.0);
        teps = 0.25 * tp->LJA[nbpidx] / pow(tsig, 12.0);
      }
      else {
        tsig = 0.0;
        teps = 0.0;
      }
    }
    else {
      tsig = 0.0;
      teps = 0.0;
    }
    fprintf(outp, "   Charge:     %9.5lf\n   LJ-Sigma:   %9.5lf\n"
	    "   LJ-Epsilon: %9.5lf\n", tp->Charges[thispt.atomid], tsig, teps);
    if (ISum(examples, tp->nxtrapt) == 1) {
      strcpy(exword, "Example: ");
    }
    else {
      strcpy(exword, "Examples:");
    }
    strcpy(spcword, "         ");
    ngiven = 0;
    for (i = 0; i < tp->nxtrapt; i++) {
      if (examples[i] == 1) {

	// Check to see that no similar example, by residue name, has
	// been given before.
	found = 0;
	for (j = 0; j < i; j++) {
	  if (examples[j] == 2) {
	    found = 1;
	  }
	}
	if (found) {
	  continue;
	}

	// Print residue, name, and ID for the extra point plus its parent atom
	if (ngiven == 0) {
	  fprintf(outp, "   %s   ", exword);
	}
	else {
	  fprintf(outp, "   %s   ", spcword);
	}
        j = tp->xtrapts[i].atomid;
        fprintf(outp, "[ %4.4s %5d -- %4.4s %6d ] -> ",
                &tp->ResNames[4*LocateResID(tp, j, 0, tp->nres)],
                LocateResID(tp, j, 0, tp->nres), &tp->AtomNames[4*j], j);
        j = tp->xtrapts[i].fr1;
        fprintf(outp, "[ %4.4s %5d -- %4.4s %6d ]\n",
                &tp->ResNames[4*LocateResID(tp, j, 0, tp->nres)],
                LocateResID(tp, j, 0, tp->nres), &tp->AtomNames[4*j], j);
	ngiven++;
	examples[i] = 2;
      }
    }
  }

  // Free allocated memory
  free(epDescribed);
  free(examples);
}

//-----------------------------------------------------------------------------
// CheckTopologyCharges: check each residue of a topology for non-integral
//                       charges.
//
// Arguments:
//   tp:       the topology to check
//   outp:     report file (mdout)
//-----------------------------------------------------------------------------
static void CheckTopologyCharges(prmtop *tp,  FILE *outp, double tol)
{
  int i, j, epOffRes, resid, parid, nprob;
  int* epResAddress;
  int* resIsEP;
  double rchg;

  // Check for extra points that are their own residues
  epOffRes = 0;
  epResAddress = (int*)malloc(tp->nxtrapt * sizeof(int));
  resIsEP = (int*)calloc(tp->nres, sizeof(int));
  for (i = 0; i < tp->nxtrapt; i++) {
    resid = LocateResID(tp, tp->xtrapts[i].atomid, 0, tp->nres);
    if (tp->ResLims[resid+1] - tp->ResLims[resid] == 1) {

      // Mark this residue as a solitary EP and find its parent residue
      resIsEP[resid] = 1;
      parid = LocateResID(tp, tp->xtrapts[i].fr1, 0, tp->nres);
      epResAddress[i] = parid;
      if (parid != resid) {
        epOffRes = 1;
      }
    }
  }
  
  // Check each residue, or the unified extra point + residue of parent atom
  nprob = 0;
  fprintf(outp, " - Checking net charges per residue:\n");
  for (i = 0; i < tp->nres; i++) {
    if (resIsEP[i] == 1) {
      continue;
    }
    rchg = DSum(&tp->Charges[tp->ResLims[i]],
                tp->ResLims[i+1] - tp->ResLims[i]);
    if (epOffRes == 1) {
      for (j = 0; j < tp->nxtrapt; j++) {
	if (epResAddress[j] == i) {
	  rchg += tp->Charges[tp->xtrapts[j].atomid];
	}
      }
    }
    if (fabs(rchg - round(rchg)) > tol) {
      fprintf(outp, " - Residue %4.4s %5d has net charge %9.5lf\n",
	      &tp->ResNames[4*i], i, rchg);
      nprob++;
    }
  }
  if (nprob == 0) {
    fprintf(outp, " - All residues have integral net charges.\n");
  }
  if (epOffRes == 1) {
    fprintf(outp, " - Extra points were found to constitute their own "
	    "residues, and to integrate\n   into other residues.  Net charges "
	    "were compiled accounting for this fact.\n");
  }

  // Finally, check the total charge on the system
  rchg = DSum(tp->Charges, tp->natom);
  if (fabs(rchg - round(rchg)) > 5.0 * tol) {
    fprintf(outp, " - The net charge on the system is %9.5lf\n", rchg);
  }
  
  // Free allocated memory
  free(epResAddress);
  free(resIsEP);
}

//-----------------------------------------------------------------------------
// EditPrmtop: edit a topology file and reprint the result.  Primarily, this is
//             for adding extra points (virtual sites) with a series of &rule
//             namelists in an eprulefile, but there could be other uses.
//
// Arguments:
//   tp:       the (now edited) topology
//   pmedit:   data structure with details as to what shall be printed
//-----------------------------------------------------------------------------
void EditPrmtop(prmtop *tp, trajcon *tj, parmed *pmedit)
{
  int EPins, iwarning;
  char suffbuffer[32];
  char* desc;
  char* basebuffer;
  FILE *outp;
  time_t ct;
  coord tc;
  
  // Create the new coordinates
  tc = InitCoords(tp, tj, 0);
  MolPlaceXpt(&tc, tp);
  
  // Print the output file as part of the main routine, since there are only
  // a handful of things to do and calls to other functions in this library
  // can handle other options.
  ct = time(&ct);
  outp = FOpenSafe(tj->outbase, tj->OverwriteOutput);
  PrintSplash(outp);
  fprintf(outp, "Run on %s", asctime(localtime(&ct)));
  
  // Make a new coordinates file, temporarily telling the writer that
  // there were no inserted extra points (that would otherwise be omitted
  // from output).  Let any user-specified restart file name override
  // the default edited coordinate file name.
  EPins = tp->EPInserted;
  tp->EPInserted = 0;
  basebuffer = (char*)malloc(MAXNAME * sizeof(char));
  strcpy(basebuffer, tj->rstbase.map[0]);
  strcpy(suffbuffer, tj->rstsuff.map[0]);
  if (strcmp(pmedit->crdfile, "inpcrd.edit") != 0) {
    strcpy(tj->rstbase.map[0], pmedit->crdfile);
    tj->rstsuff.map[0][0] = '\0';
  }
  WriteRst(NULL, &tc, tp, tj, 0, -1);
  tp->EPInserted = EPins;
  strcpy(tj->rstbase.map[0], basebuffer);
  strcpy(tj->rstsuff.map[0], suffbuffer);
  
  // Reprint the input file
  ReprintInputFile(tj, "NONE", "NONE", 0, outp);
  
  // Print the new topology
  HorizontalRule(outp, 1);
  desc = (char*)malloc(MAXLINE * sizeof(char));
  PrintVADesc(0, "(1.)", 4, " ", 1, "Printing System Topology.\n", 77, 0,
	      outp);
  if (pmedit->topchk == 1) {
    fprintf(outp, " - Checking chirality of standard amino acids.\n");
    iwarning = ProteinChiralityCheck(tp, &tc, outp);
    fprintf(outp, " - Checking for omitted disulfides.\n");
    iwarning += FindDisulfides(tp, &tc, outp);
    fprintf(outp, " - Checking for (nonstandard) Lennard-Jones combining "
            "rules.\n");
    iwarning += CheckLJRules(tp, outp);
    if (iwarning == 0) {
      fprintf(outp, " - No problems with the topology detected.\n\n");
    }
    CheckTopologyCharges(tp,  outp, pmedit->qtol);
  }
  else {
    fprintf(outp, " - Check skipped.\n");
  }
  if (tp->nxtrapt > 0) {
    ReportExtraPoints(tp, outp);
  }
  HorizontalRule(outp, 1);
  fclose(outp);
  PutPrmTop(tp, pmedit->prmfile, pmedit->prmtitle);

  // Free allocated memory
  free(desc);
  free(basebuffer);
  
  // Exit, we are done.
  exit(0);
}
