#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "topRead.h"

int main(int argc, char *argv[])
{
  int i, j, checkBub, checkChir;
  double gspc, bubdi;
  char rmapname[MAXNAME], tag[MAXNAME], tmptag[MAXNAME];
  coord tc;
  prmtop tp;

  /*** Options ***/
  if (argc == 1) {
    printf("\nCheckMD >> A program for automated checking of an MD simulation."
	   "\n\nOptions:\n"
	   "  -p        : the topology file (new prmtop format)\n"
	   "  -c        : the coordinate file (default format .crd)\n"
	   "  -fc       : the coordinate file format (\"CRD\" or \"RST\")\n"
	   "  -V        : check for vaccuum bubbles (not done by default)\n"
           "  -G        : grid spacing to use when checking for vacuum "
	   "bubbles\n"
           "  -B        : bubble radius to test when checking for vacuum "
	   "bubbles (default 2.0A)\n"
	   "  -Nochiral : skip chirality check\n"
           "  -Rmap     : residue map file name (default 'resID.txt')\n\n");
    exit(1);
  }

  tp.source[0] = '\0';
  tc.source[0] = '\0';
  tc.rst = 0;
  sprintf(rmapname, "resID.txt");
  checkBub = 0;
  checkChir = 1;
  gspc = 0.2;
  bubdi = 2.0;

  for (i = 0; i < argc-1; i += 2) {
    strcpy(tag, *++argv);
    if (strcmp(tag, "-p") == 0) {
      strcpy(tp.source, *++argv);
    }
    else if (strcmp(tag, "-c") == 0) {
      strcpy(tc.source, *++argv);
    }
    else if (strcmp(tag, "-fc") == 0) {
      strcpy(tmptag, *++argv);
      for (j = 0; j < strlen(tmptag); j++) {
	tmptag[j] = toupper(tmptag[j]);
      }
      if (strcmp(tmptag, "CRD") == 0) {
	tc.rst = 0;
      }
      else if (strcmp(tmptag, "RST") == 0) {
	tc.rst = 1;
      }
      else {
	printf("CheckMD >> Error.  Coordinate file type %s unkonwn.\n\n",
	       tmptag);
	exit(1);
      }
    }
    else if (strcmp(tag, "-V") == 0) {
      checkBub = 1;
      i--;
    }
    else if (strcmp(tag, "-G") == 0) {
      gspc = atof(*++argv);
    }
    else if (strcmp(tag, "-B") == 0) {
      bubdi = atof(*++argv);
    }
    else if (strcmp(tag, "-NoChiral") == 0) {
      checkChir = 0;
    }
    else if (strcmp(tag, "-Rmap") == 0) {
      strcpy(rmapname, *++argv);
    }
    else {
      printf("CheckMD >> Error.  Identifier %s unknown.\n\n", tag);
      exit(1); 
    }
  }

  /*** Check input ***/
  if (tc.source[0] == '\0') {
    printf("CheckMD >> Error.  Coordinate file not specified!\n");
  }
  if (tp.source[0] == '\0') {
    printf("CheckMD >> Error.  Topology file not specified!\n");
  }

  GetPrmTop(&tp, 1, 1);
  GetRst(&tc, &tp);
  MapResidues(&tp, &tc, rmapname);
  FindDisulfides(&tp, &tc);
  if (checkBub == 1) {
    VacuumBubble(&tp, &tc, gspc, bubdi);
  }
  if (checkChir == 1) {
    ProteinChiralityCheck(&tp, &tc);
  }

  return 0;
}
