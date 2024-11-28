#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "pdbRead.h"
#include "matrix.h"

/***=======================================================================***/
/*** GetNames: reads a file containing a list of file names (or, in        ***/
/***           general, words, one per line).  Any line beginning with '%' ***/
/***           is assumed to be a comment.  All other lines are assumed to ***/
/***           contain one character string.  If more than one is present, ***/
/***           only the first string is read.  If the flag "cex" is set to ***/
/***           1, each file name will be tested for existence.             ***/
/***=======================================================================***/
cmat GetNames(char* ptsfile, int cex)
{
  int i;
  char testname[MAXNAME], line[MAXLINE];
  cmat names;
  FILE *pts, *test;

  names.row = 0;
  if ((pts = fopen(ptsfile, "r")) == NULL) {
    printf("GetNames >> Error.  Calibration Points File %s not found!\n",
	   ptsfile);
  }
  while(fgets(line, MAXLINE, pts) != NULL) {
    remove_whitespace(line, MAXLINE);
    if (line[0] != '%') {
      names.row += 1;
      if (cex == 1) {
        sscanf(line, "%s", testname);
        if ((test = fopen(testname, "r")) != NULL) {
	  fclose(test);
        }
        else {
	  printf("GetNames >> Error.  File %s not found!\n", testname);
	  exit(1);
	}
      }
    }
  }
  rewind(pts);

  /*** Read calibration point names ***/
  i = 0;
  names = CreateCmat(names.row, 512);
  while(fgets(line, MAXLINE, pts) != NULL) {
    remove_whitespace(line, MAXLINE);
    if (line[0] != '%') {
      sscanf(line, "%s", names.map[i]);
      i++;
    }
  }
  fclose(pts);

  printf("GetNames >> Found %d names in %s\n\n", names.row, ptsfile);
  return names;
}

/***=======================================================================***/
/*** GetDoubles: reads a list of double-precision numbers, one per line.   ***/
/***=======================================================================***/
double* GetDoubles(char* fname, int *nval)
{
  int i, maxi;
  double* d;
  char line[256];
  FILE *inp;

  if ((inp = fopen(fname, "r")) == NULL) {
    printf("GetDoubles >> Error.  File %s not found.\n", fname);
  }
  i = 0;
  maxi = 4096;
  d = (double*)malloc(maxi*sizeof(double));
  while(fgets(line, 256, inp) != NULL) {
    remove_whitespace(line, 256);
    if (line[0] != '%') {
      sscanf(line, "%lf", &d[i]);
      i++;
    }
    if (i == maxi-1) {
      maxi += 4096;
      d = (double*)realloc(d, maxi*sizeof(double));
    }
  }
  fclose(inp);
  d = (double*)realloc(d, i*sizeof(double));
  *nval = i;

  return d;
}
