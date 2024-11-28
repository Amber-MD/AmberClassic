#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "pdbRead.h"
#include "crdmanip.h"
#include "matrix.h"
#define PI                     3.1415926535898

int main(int argc, char *argv[])
{
  int i;
  char tag[64], outname[256];
  pdb p, t;
  symT S;

  if (argc == 1) {
    printf("\nUnitCell >> A program for recreating a crystallographic unit "
	   "cell from a PDB\n""UnitCell >> structure.\n\n"
	   "Options:\n"
           "  -p   : the structure to reassemble (PDB format)\n"
           "  -o   : the output structure (PDB format)\n");
    exit(0);
  }

  /*** Input ***/
  p.source[0] = '\0';
  outname[0] = '\0';
  for (i = 0; i < argc-2; i += 2) {
    strcpy(tag, *++argv);
    if (strcmp(tag, "-p") == 0) {
      strcpy(p.source, *++argv);
    }
    else if (strcmp(tag, "-o") == 0) {
      strcpy(outname, *++argv);
    }
    else {
      printf("UnitCell >> Error.  Tag %s not recognized.\n", tag);
      exit(1);
    }
  }

  /*** Check input ***/
  if (p.source[0] == '\0') {
    printf("UnitCell >> Error.  Original PDB file not specified!\n");
    exit(1);
  }
  if (outname[0] == '\0') {
    printf("UnitCell >> Error.  Output PDB file not specified!\n");
    exit(1);
  }

  /*** Get the original PDB ***/
  GetPDB(&p, 0, 0);

  /*** Scan PDB file for symmetry information ***/
  S = LoadSymmetryData(p.source);

  /*** Create the new PDB ***/
  t = ExpandPDBSym(&p, &S, 1);
  strcpy(t.source, outname);

  /*** Print PDB ***/
  ModPdbRA(&t);
  PutPDB(&t, &S, t.source, "STANDARD", "HEADER  ", 0);

  return 0;
}
