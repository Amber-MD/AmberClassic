#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "Matrix.h"
#include "Parse.h"
#include "CrdManip.h"
#include "Topology.h"
#include "Trajectory.h"
#include "SinglePointEval.h"
#include "mdgxVector.h"
#include "Manual.h"
#include "Constants.h"

//-----------------------------------------------------------------------------
// DetectQMOutput: detect the output of a quantum mechanics program.  Return
//                 an integer code to indicate the result.  0 = not QM output,
//                 1 = ORCA output, 2 = Gaussian output, 3 = MolPro output,
//                 4 = GAMESS output, 5 = simple list of numbers.
//
// Arguments:
//   fname:    the name of the file to test
//-----------------------------------------------------------------------------
static int DetectQMOutput(char* fname)
{
  int i, j, ftype;
  cmat cfi;

  // Read the file into memory
  cfi = Ascii2Mem(fname, 256, 8, "File not found");

  // Tests for each QM output type
  ftype = 0;
  if (DetectInCmat(&cfi, "* O   R   C   A *", &i, &j, 0, 0) == 1) {
    ftype = 1;
  }
  if (DetectInCmat(&cfi, "This software contains proprietary and confidential",
		   &i, &j, 0, 0) == 1 &&
      DetectInCmat(&cfi, "including trade secrets, belonging to Gaussian, Inc",
		   &i, &j, i, 0) == 1) {
    ftype = 2;
  }

  // Free allocated memory
  DestroyCmat(&cfi);

  // If there hasn't yet been a file type detected,
  // try looking at it as just a plain list of numbers. 
  if (ftype == 0 && ValidRealData(fname) == 1) {
    ftype = 5;
  }

  // Return the code that was matched
  return ftype;
}

//-----------------------------------------------------------------------------
// ReadOrcaSP: read an ORCA output file describing a single point energy
//             calculation.
//
// Arguments:
//   fname:    the name of the ORCA output file
//   molcrd:   coordinates to use (if NULL, the coordinates in the file will
//             be used instead)
//   tp:       system topology
//   info:     data structure to store the results, including the coordinates
//             and energy
//-----------------------------------------------------------------------------
static void ReadOrcaSP(char* fname, double* molcrd, prmtop *tp, qmresult *info)
{
  int i, j, k, triali, trialj;
  cmat cfi, Lwords;
  
  // Read in the file
  cfi = Ascii2Mem(fname, 256, 8, "File could not be read.");

  // Default settings for the QM calculation
  info->pass = 1;
  info->nbfunc = -1;
  sprintf(info->prog, "orca");
  
  // Checks to ensure that the file is complete
  info->pass = (DetectInCmat(&cfi, "****ORCA TERMINATED NORMALLY****",
			     &i, &j, 0, 0) == 1) ? info->pass : 0;

  // Get information from the file that may be of interest
  info->pass = (DetectInCmat(&cfi, "# of primitive gaussian functions",
			     &i, &j, 0, 0) == 1) ? info->pass : 0;
  if (info->pass) {
    if (DetectInCmat(&cfi, "...", &i, &j, 0, 0) == 1) {
      sscanf(&cfi.map[i][j], "%d", &info->nbfunc);
    }
  }
  if (DetectInCmat(&cfi, "ORCA SCF", &i, &j, 0, 0) == 1) {

    // The method will currently be logged as "SCF", but this may change
    // if a post Hartree-Fock method is called later 
    sprintf(info->method, "scf");
    if (DetectInCmat(&cfi, "Nuclear Repulsion      ENuc", &i, &j, i, 0) == 1) {
      DoubleFromLine(cfi.map[i], 4, &info->nucrepE);
    }
    if (DetectInCmat(&cfi, "Energy Change          TolE", &i, &j, i, 0) == 1) {
      DoubleFromLine(cfi.map[i], 4, &info->scftol);
    }
  }
  if (DetectInCmat(&cfi, "ORCA MP2 CALCULATION", &i, &j, 0, 0) == 1) {
    sprintf(info->method, "mp2");
    info->pass = (DetectInCmat(&cfi, "MP2 CORRELATION ENERGY   :",
			       &i, &j, i, 0) == 1) ? info->pass : 0;
    if (info->pass == 1) {
      sscanf(&cfi.map[i][j], "%lf", &info->mp2corrE);
    }
  }

  // This is the big deal--the single point energy
  info->pass = (DetectInCmat(&cfi, "FINAL SINGLE POINT ENERGY",
			     &i, &j, 0, 0) == 1) ? info->pass : 0;
  if (info->pass == 1) {
    while (DetectInCmat(&cfi, "FINAL SINGLE POINT ENERGY",
                        &triali, &trialj, i+1, 0) == 1) {
      i = triali;
      j = trialj;
    }
    sscanf(&cfi.map[i][j], "%lf", &info->spE);
  }

  // Now, do we want the coordinates?
  if (molcrd == NULL) {
    info->pass = (DetectInCmat(&cfi, "CARTESIAN COORDINATES (ANGSTROEM)",
			       &i, &j, 0, 0) == 1) ? info->pass : 0;
    info->crd = (double*)malloc(3*tp->natom*sizeof(double));
    if (info->pass == 1) {
      while (DetectInCmat(&cfi, "CARTESIAN COORDINATES (ANGSTROEM)",
                          &triali, &trialj, i+1, 0) == 1) {
        i = triali;
        j = trialj;
      }
      for (j = 0; j < tp->natom; j++) {
	Lwords = ParseWords(cfi.map[i+2+j]);
	if (Lwords.row != 4) {
	  info->pass = 0;
	}
	for (k = 0; k < 3; k++) {
	  if (WordIsNumber(Lwords.map[k+1])) {
	    info->crd[3*j+k] = atof(Lwords.map[k+1]);
	  }
	  else {
	    info->pass = 0;
	  }
	}
	DestroyCmat(&Lwords);
      }

      // Test the next line: it should be blank
      if (CountWords(cfi.map[i+2+tp->natom]) != 0) {
	info->pass = 0;
      }
    }
  }
  else {
    info->crd = CpyDVec(molcrd, 3*tp->natom);
  }

  // Free allocated memory
  DestroyCmat(&cfi);
}

//-----------------------------------------------------------------------------
// ReadGaussianSP: read a Gaussian output file describing a single point
//                 energy calculation.
//
// Arguments:
//   fname:    the name of the Gaussian output file
//   molcrd:   coordinates to use (if NULL, the coordinates in the file will
//             be used instead)
//   tp:       system topology
//   info:     data structure to store the results, including the coordinates
//             and energy
//-----------------------------------------------------------------------------
static void ReadGaussianSP(char* fname, double* molcrd, prmtop *tp,
			   qmresult *info)
{
  int i, j, k;
  char* lbuff;
  cmat cfi, Lwords;

  // Read in the file
  cfi = Ascii2Mem(fname, 256, 8, "File could not be read.");
  lbuff = (char*)malloc(MAXLINE*sizeof(char));

  // Default settings for the QM calculation
  info->pass = 1;
  info->nbfunc = -1;
  sprintf(info->prog, "gaussian");

  // Checks to ensure that the file is complete
  info->pass = (DetectInCmat(&cfi, " Normal termination of Gaussian", &i, &j,
                             cfi.row-1, 1) == 1) ? info->pass : 0;

  // Get information from the file that may be of interest
  info->pass = (DetectInCmat(&cfi, "symmetry adapted cartesian basis func",
                             &i, &j, 0, 0) == 1) ? info->pass : 0;
  if (info->pass) {
    sscanf(&cfi.map[i][10], "%d", &info->nbfunc);
  }
  if (DetectInCmat(&cfi, " SCF Done:  E(", &i, &j, 0, 0) == 1) {

    // The method will currently be logged as "SCF", but this may change
    // if a post Hartree-Fock method is called later
    sprintf(info->method, "scf");
    if (DetectInCmat(&cfi, "EE=", &i, &j, i, 0) == 1) {
      strcpy(lbuff, cfi.map[i]);
      EqualSpace(lbuff);
      DoubleFromLine(lbuff, 8, &info->nucrepE);
    }
    if (DetectInCmat(&cfi, "Delta-E=", &i, &j, i, 1) == 1) {
      strcpy(lbuff, cfi.map[i]);
      EqualSpace(lbuff);
      DoubleFromLine(cfi.map[i], 5, &info->scftol);
    }
  }
  if (DetectInCmat(&cfi, "EUMP2 =", &i, &j, 0, 0) == 1) {
    sprintf(info->method, "mp2");
    strcpy(lbuff, cfi.map[i]);
    EqualSpace(lbuff);
    info->pass = (DoubleFromLine(lbuff, 2,
                                 &info->mp2corrE) == 0) ? info->pass : 0;
    info->pass = (DoubleFromLine(lbuff, 5,
                                 &info->spE) == 0) ? info->pass : 0;
  }

  // Now, do we want the coordinates?
  if (molcrd == NULL) {
    sprintf(lbuff, "Number     Number       Type             X           Y");
    info->pass = (DetectInCmat(&cfi, lbuff,
                               &i, &j, 0, 0) == 1) ? info->pass : 0;
    info->crd = (double*)malloc(3*tp->natom*sizeof(double));
    if (info->pass == 1) {
      for (j = 0; j < tp->natom; j++) {
        for (k = 0; k < 3; k++) {
	  if (DoubleFromLine(cfi.map[i+2+j], k+3, &info->crd[3*j+k]) == 1) {
            info->pass = 0;
	  }
        }
      }

      // Test the next line: it should be a bunch of dashes
      if (CountWords(cfi.map[i+2+tp->natom]) != 1) {
	info->pass = 0;
      }
    }
  }
  else {
    info->crd = CpyDVec(molcrd, 3*tp->natom);
  }

  // Free allocated memory
  free(lbuff);
}

//-----------------------------------------------------------------------------
// FilterSPData: this function will clean up a lot of things that could have
//               gone wrong with the aggregate data set.  Checks on individual
//               input files were performed as they were read in, but now that
//               all of the data is in memory we can check to see that it all
//               makes sense together.
//
// Arguments:
//   spelist:     list of single point energies to parse, along with the
//                structures they describe
//   tp:          the system topology
//-----------------------------------------------------------------------------
static void FilterSPData(spdata *spelist, prmtop *tp)
{
  int i, j, k, np, nclst, maxloc;
  int* itemsrc;
  int* ptsrc;
  double clstcenter;
  double* energies;

  // Count the total number of data points.  Make a
  // list of all of them and where they came from.
  np = 0;
  for (i = 0; i < spelist->nitems; i++) {
    np += spelist->items[i].npass;
    for (j = 0; j < spelist->items[i].npts; j++) {
      spelist->items[i].points[j].ccremoved = 0;
    }
  }
  if (np == 0) {
    return;
  }
  energies = (double*)malloc(np*sizeof(double));
  itemsrc = (int*)malloc(np*sizeof(int));
  ptsrc = (int*)malloc(np*sizeof(int));
  k = 0;
  for (i = 0; i < spelist->nitems; i++) {
    for (j = 0; j < spelist->items[i].npts; j++) {
      if (spelist->items[i].points[j].pass == 0) {
	continue;
      }

      // The energy from the single point calculations will
      // oftgen be in Hartrees.
      if (strcmp(spelist->items[i].points[j].prog, "orca") == 0 ||
	  strcmp(spelist->items[i].points[j].prog, "gaussian") == 0 ||
	  strcmp(spelist->items[i].points[j].prog, "molpro") == 0 ||
	  strcmp(spelist->items[i].points[j].prog, "gamess") == 0) {
	energies[k] = spelist->items[i].points[j].spE * 627.509;
      }
      itemsrc[k] = i;
      ptsrc[k] = j;
      k++;
    }
  }

  // First check: does the data have a reasonable distribution of energies?
  if (spelist->enfCluster == 1) {

    // Starting from each data point, determine how much of the
    // data set can be fit into a cluster of the specified width
    // using that data point as the center.
    spelist->maxclst = 0;
    maxloc = -1;
    for (i = 0; i < np; i++) {

      // First method: clusters centerd on actual data points
      nclst = 0;
      for (j = 0; j < np; j++) {
	if (fabs(energies[j] - energies[i]) <= spelist->clstWidth) {
	  nclst++;
	}
      }
      if (nclst > spelist->maxclst) {
	spelist->maxclst = nclst;
        maxloc = i;
	clstcenter = energies[i];
      }

      // Next method: clusters with an upper bound of an actual data point
      nclst = 0;
      for (j = 0; j < np; j++) {
        if (j == i || (energies[i] - energies[j] >= 0.0 &&
                       energies[i] - energies[j] <= 2.0*spelist->clstWidth)) {
          nclst++;
        }
      }
      if (nclst > spelist->maxclst) {
        spelist->maxclst = nclst;
        maxloc = i;
	clstcenter = energies[i] - spelist->clstWidth;
      }

      // Final method: clusters with a lower bound of an actual data point
      nclst = 0;
      for (j = 0; j < np; j++) {
        if (j == i || (energies[j] - energies[i] >= 0.0 &&
                       energies[j] - energies[i] <= 2.0*spelist->clstWidth)) {
          nclst++;
        }
      }
      if (nclst > spelist->maxclst) {
        spelist->maxclst = nclst;
        maxloc = i;
	clstcenter = energies[i] + spelist->clstWidth;
      }
    }

    // Points not within the best cluster no longer pass
    spelist->clstCull = 0;
    for (i = 0; i < np; i++) {
      if (fabs(energies[i] - clstcenter) > spelist->clstWidth) {
	spelist->items[itemsrc[i]].points[ptsrc[i]].pass = 0;
	spelist->items[itemsrc[i]].points[ptsrc[i]].ccremoved = 1;
	spelist->clstCull += 1;
      }
    }
  }

  // Free allocated memory
  free(itemsrc);
  free(ptsrc);
  free(energies);
}

//-----------------------------------------------------------------------------
// PrintSPStructures: this function works very much like PrintConfigurations in
//                    the ConfigSamp.c library, but it is coded separately due
//                    to the different layout of the structures to be printed
//                    and the specific nature of the output.
//
// Arguments:
//   spelist:     list of single point energies to parse, along with the
//                structures they describe
//   tj:          trajectory control input data (this contains the output file
//                name and a flag to tell whether it's OK to overwrite output)
//   tp:          the system topology
//-----------------------------------------------------------------------------
static void PrintSPStructures(spdata *spelist, trajcon *tj, prmtop *tp)
{
  int i, j, nfrm;
  char* trajname;
  char* msg;
  FILE *nrgout;
  coord tc;
  cdftrj Acdf;

  // Count the number of structures to print
  nfrm = 0;
  for (i = 0; i < spelist->nitems; i++) {
    nfrm += spelist->items[i].npass;
  }
  
  // Make sure that, when splicing file names, there will be no
  // attempt to add extra numbers based on uninitialized integers.
  // Also rig this so that every configuration will be written as
  // a "frame" of the trajectory to a single file.
  tj->nfistep = 0;
  tj->irest = 0;
  tj->ntwx = 1;
  tj->ntpr = -1;
  tj->currstep = 0;
  tj->nstep = nfrm;

  // Loop over all single point evaluation tasks, and within
  // each task loop over the data points that passed.
  nrgout = FOpenSafe(tj->dumpname, tj->OverwriteOutput);
  trajname = (char*)malloc(MAXNAME*sizeof(char));
  SpliceFileName(tj, tj->trjbase.map[0], tj->trjsuff.map[0], trajname, 0); 
  msg = (char*)malloc(MAXLINE*sizeof(char));
  sprintf(msg, "The following energies apply to conformations described "
	  "by the topology %s.  Structures are given in trajectory "
	  "%s, count: %d.", tp->source, trajname, nfrm);
  PrintParagraph(msg, 79, "%", nrgout);
  free(msg);
  free(trajname);

  // Print the trajectory in NetCDF format to preserve precision
  tc = CreateCoord(tp->natom);
  for (i = 0; i < 3; i++) {
    tc.gdim[i] = 128.0;
    tc.gdim[i+3] = 0.5*PI;
  }
  WriteCDF(NULL, &tc, 1, tj, &Acdf, tp, 0);
  tj->currstep = 1;
  for (i = 0; i < spelist->nitems; i++) {
    if (spelist->items[i].npass == 0) {
      continue;
    }
    fprintf(nrgout, "%% Item %d:\n", i);
    if (spelist->items[i].takeqmcrd == 0) {
      fprintf(nrgout, "%%   Coords -> %s\n", spelist->items[i].crdsrc);
    }
    fprintf(nrgout, "%%   Energy -> %s\n", spelist->items[i].esrc);    
    for (j = 0; j < spelist->items[i].npts; j++) {
      if (spelist->items[i].points[j].pass != 1) {
	continue;
      }
      ReflectDVec(tc.loc, spelist->items[i].points[j].crd, 3*tp->natom);
      WriteCDF(NULL, &tc, 1, tj, &Acdf, tp, 0);
      fprintf(nrgout, "%16.8lf\n", spelist->items[i].points[j].spE);
      tj->currstep += 1;
    }
  }
  fclose(nrgout);

  // Free allocated memory
  DestroyCoord(&tc);
}

//-----------------------------------------------------------------------------
// PrintSPReport: print a report on the single point evaluations that have been
//                done, highlighting failures and any rearrangements that had
//                to be done.
//
// Arguments:
//   spelist:     list of single point energies to parse, along with the
//                structures they describe
//   tj:          trajectory control input data (this contains the output file
//                name and a flag to tell whether it's OK to overwrite output)
//   tp:          the system topology
//-----------------------------------------------------------------------------
static void PrintSPReport(spdata *spelist, trajcon *tj, prmtop *tp)
{
  int i, j, totalread, totalpass;
  char* msg;
  char* trajname;
  FILE *outp;
  time_t ct;

  // Input file header
  ct = time(&ct);
  outp = FOpenSafe(tj->outbase, tj->OverwriteOutput);
  PrintSplash(outp);
  fprintf(outp, "Run on %s", asctime(localtime(&ct)));
  ReprintInputFile(tj, "data", "DataPoint", spelist->PrintItems, outp);

  // Print a summary of the successes and failures
  totalread = 0;
  totalpass = 0;
  for (i = 0; i < spelist->nitems; i++) {
    totalread += spelist->items[i].npts;
    for (j = 0; j < spelist->items[i].npts; j++) {
      if (spelist->items[i].points[j].pass == 1) {
	totalpass++;
      }
    }
  }
  HorizontalRule(outp, 0);
  msg = (char*)malloc(MAXLINE*sizeof(char));
  trajname = (char*)malloc(MAXNAME*sizeof(char));
  SpliceFileName(tj, tj->trjbase.map[0], tj->trjsuff.map[0], trajname, 0);
  sprintf(msg, "Summary of successes and failures.  A total of %d structures "
          "and energies were read in, of which %d passed inspection and were "
          "then written to the output trajectory %s and the energy data file "
          "%s.\n", totalread, totalpass, trajname, tj->dumpname);
  PrintVADesc(0, "(1.)", 4, " ", 1, msg, 74, 0, outp);
  HorizontalRule(outp, 1);

  // Print a summary of actions taken to filter the results
  if (spelist->enfCluster == 1) {
    HorizontalRule(outp, 0);
    sprintf(msg, "Structures and energies were checked to ensure that the "
	    "collection of data points stays within a sensible range.  A "
	    "clustering of structures whose energies all fell within %.2lf "
	    "kcal/mol were taken, while structures that could not be included "
	    "in this group were left out.  IF THE RESULTS SEEM STRANGE, "
            "please check the units of the single point energies and note "
            "that mdgx is trying to make appropriate adjustments but may fail "
            "if the source of the data is unfamiliar.", spelist->clstWidth);
    PrintVADesc(0, "(2.)", 4, " ", 1, msg, 74, 0, outp);
    if (spelist->clstCull > 0) {
      fprintf(outp, "\n Data points that were culled due to this filter:\n");
      for (i = 0; i < spelist->nitems; i++) {
        for (j = 0; j < spelist->items[i].npts; j++) {
	  if (spelist->items[i].points[j].ccremoved == 1) {
	    fprintf(outp, "  - %s (energy %12.4lf)\n",
		    spelist->items[i].points[j].source,
		    spelist->items[i].points[j].spE);
	  }
	}
      }
    }
    else {
      fprintf(outp, " All structures's energies could be grouped within the "
	      "specified range.\n");
    }
    HorizontalRule(outp, 1);
  }

  // Close the output file
  fclose(outp);

  // Free allcoated memory
  free(msg);
  free(trajname);
}

//-----------------------------------------------------------------------------
// ReadSinglePointEnergies: this central function for the small &speval module
//                          looks in a series of files or directories
//                          containing quantum outputs.  The quantum energies
//                          in these directories will be paired against the
//                          expressed coordinates corresponding to them as a
//                          trajectory and list of energies for later use as
//                          inputs to the parameter fitting routines.
//
// Arguments:
//   spelist:     list of single point energies to parse, along with the
//                structures they describe
//   tj:          trajectory control input data (this contains the output file
//                name and a flag to tell whether it's OK to overwrite output)
//   tp:          the system topology
//-----------------------------------------------------------------------------
void ReadSinglePointEnergies(spdata *spelist, trajcon *tj, prmtop *tp)
{
  int i, j, k, nconf, problem, usecoord, usedmat, ndata, fitype;
  int crdtype, nrgtype;
  double tstamp;
  cmat allnames;
  coord tc;
  dmat ttraj, nrgvals;

  i = 0;
  problem = 0;
  ndata = 0;
  while (i < spelist->nitems) {

    // Report issues
    if (problem == 1) {
      printf("ReadSinglePointEnergies >> If specifying coordinates "
	     "separate from energies,\nReadSinglePointEnergies >> there "
	     "must be precisely one energy and one set\n"
	     "ReadSinglePointEnergies >> of coordinates specified, or a "
	     "series of\nReadSinglePointEnergies >> conformations in a "
	     "single trajectory file given\nReadSinglePointEnergies >> "
	     "alongside a file containing a list with the\n"
	     "ReadSinglePointEnergies >> same number of energies.  The "
	     "\"data point\" pairing\nReadSinglePointEnergies >> %s and %s "
	     "violates this rule.\n", spelist->items[i].crdsrc,
	     spelist->items[i].esrc);
      problem = 0;
      i++;
    }

    // Determine the file types
    if (spelist->items[i].takeqmcrd == 0) {
      crdtype = FindFileType(spelist->items[i].crdsrc);
    }
    nrgtype = FindFileType(spelist->items[i].esrc);

    // Case 1: take one set of coordinates and the
    //         energy directly from a QM output file.
    if (spelist->items[i].takeqmcrd == 1 && nrgtype == 0) {
      spelist->items[i].npts = 1;
      spelist->items[i].points = (qmresult*)malloc(sizeof(qmresult));
      nrgtype = DetectQMOutput(spelist->items[i].esrc);
      if (nrgtype == 1) {
	ReadOrcaSP(spelist->items[i].esrc, NULL, tp,
		   &spelist->items[i].points[0]);
      }
      spelist->items[i].points[0].source = (char*)malloc(MAXNAME*sizeof(char));
      sprintf(spelist->items[i].points[0].source, "%s", spelist->items[i].esrc);
    }

    // Case 2: take sets of coordinates and corresponding
    //         energies from a set of QM output files.
    else if (spelist->items[i].takeqmcrd == 1 &&
	     (nrgtype == 1 || nrgtype == 2)) {
      if (nrgtype == 1) {
        allnames = DirectoryFileSearch(spelist->items[i].esrc);
      }
      else {
        allnames = RegExpFileSearch(spelist->items[i].esrc);
      }
      spelist->items[i].npts = allnames.row;
      spelist->items[i].points = (qmresult*)malloc(allnames.row *
						   sizeof(qmresult));
      for (j = 0; j < allnames.row; j++) {
	nrgtype = DetectQMOutput(allnames.map[j]);
        if (nrgtype == 1) {
          ReadOrcaSP(allnames.map[j], NULL, tp, &spelist->items[i].points[j]);
	}
	else if (nrgtype == 2) {
	  ReadGaussianSP(allnames.map[j], NULL, tp,
			 &spelist->items[i].points[j]);
	}
	else {

	  // Blank or dud file	  
	  spelist->items[i].points[j].pass = 0;
	  spelist->items[i].points[j].crd = (double*)malloc(2*sizeof(double));
	}
        spelist->items[i].points[j].source =
          (char*)malloc(MAXNAME*sizeof(char));
        sprintf(spelist->items[i].points[j].source, "%s", allnames.map[j]);
      }
      DestroyCmat(&allnames);
    }
    else if (spelist->items[i].takeqmcrd == 0 &&
	     crdtype == 0 && nrgtype == 0) {
      crdtype = ValidCoord(spelist->items[i].crdsrc);
      nrgtype = DetectQMOutput(spelist->items[i].esrc);

      // Case 3: take coordinates from a trajectory, inpcrd,
      //         or restart file in ascii or NetCDF format).
      //         Read energies from an ascii list.
      if (crdtype > 0 && nrgtype == 5) {
	nrgvals = ReadRealData(spelist->items[i].esrc);
	if (crdtype == 1 || crdtype == 2 || crdtype == 4) {
	  tstamp = 0.0;
	  tc = ReadRst(tp, spelist->items[i].crdsrc, &tstamp);
	  spelist->items[i].npts = 1;
	  spelist->items[i].points = (qmresult*)malloc(sizeof(qmresult));
	  if (nrgvals.col != 1) {
	    problem = 1;
	    spelist->items[i].npts = 0;
	    DestroyCoord(&tc);
	    continue;
	  }
          sprintf(spelist->items[i].points[0].prog, "amber");
          sprintf(spelist->items[i].points[0].method, "unknown");
	  spelist->items[i].points[0].spE = nrgvals.data[0];
	  spelist->items[i].points[0].crd = CpyDVec(tc.loc, 3*tp->natom);
	  DestroyCoord(&tc);
	  spelist->items[i].points[0].pass = 1;
          spelist->items[i].points[0].source =
            (char*)malloc(MAXNAME*sizeof(char));
          sprintf(spelist->items[i].points[0].source, "%s", spelist->items[i].esrc);
	}
	else if (crdtype == 3 || crdtype == 5) {
	  if (crdtype == 3) {
	    ttraj = ReadCrdTraj(tp, spelist->items[i].crdsrc, 0);
	  }
	  else {
	    ttraj = ReadCDFTraj(tp, spelist->items[i].crdsrc, 0);
	  }
          spelist->items[i].npts = ttraj.row;
          spelist->items[i].points = (qmresult*)malloc(ttraj.row *
						       sizeof(qmresult));
          if (nrgvals.col != ttraj.col) {
            problem = 1;
            spelist->items[i].npts = 0;
	    DestroyDmat(&ttraj);
	    continue;
          }
	  for (j = 0; j < ttraj.col; j++) {
            sprintf(spelist->items[i].points[j].prog, "amber");
            sprintf(spelist->items[i].points[j].method, "unknown");
	    spelist->items[i].points[j].spE = nrgvals.data[j];
	    spelist->items[i].points[j].crd = (double*)malloc(3*tp->natom *
							      sizeof(double));
	    for (k = 0; k < 3*tp->natom; k++) {
	      spelist->items[i].points[j].crd[k] = ttraj.map[k][j];
	    }
	    spelist->items[i].points[j].pass = 1;
            spelist->items[i].points[j].source =
              (char*)malloc(MAXNAME*sizeof(char));
            sprintf(spelist->items[i].points[j].source, "%s / %s, frame %d",
                    spelist->items[i].esrc, spelist->items[i].crdsrc, j);
	  }
	  DestroyDmat(&ttraj);
	}
      }

      // Case 4: take coordinates from an inpcrd or restart
      //         file (ascii or NetCDF format).  Read the
      //         energy from a QM output file.
      if ((crdtype == 1 || crdtype == 2 || crdtype == 4) &&
	  (nrgtype >= 1 && nrgtype <= 4)) {
	tstamp = 0.0;
	tc = ReadRst(tp, spelist->items[i].crdsrc, &tstamp);
	spelist->items[i].npts = 1;
	spelist->items[i].points = (qmresult*)malloc(sizeof(qmresult));
	if (nrgtype == 1) {
	  ReadOrcaSP(spelist->items[i].crdsrc, tc.loc, tp,
		     &spelist->items[i].points[0]);
	}
        spelist->items[i].points[0].source =
          (char*)malloc(MAXNAME*sizeof(char));
        sprintf(spelist->items[i].points[0].source, "%s / %s",
                spelist->items[i].esrc, spelist->items[i].crdsrc);
	DestroyCoord(&tc);
      }
    }

    // Increment the loop control variable
    i++;
  }

  // Count the number of frames that pass muster
  for (i = 0; i < spelist->nitems; i++) {
    spelist->items[i].npass = 0;
    for (j = 0; j < spelist->items[i].npts; j++) {
      if (spelist->items[i].points[j].pass == 1) {
	spelist->items[i].npass += 1;
      }
    }
  }

  // Filter the results: in a lot of cases there could be outliers, and
  // it is tedious to have to remove them by hand.
  FilterSPData(spelist, tp);

  // Print a trajectory with all of the coordinates and a file with all of
  // the energies for conformations with energies that passed.
  PrintSPStructures(spelist, tj, tp);

  // Print the output report and exit
  PrintSPReport(spelist, tj, tp);
  exit(1);
}
