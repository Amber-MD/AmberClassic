#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stringDefs.h"

int main(int argc, char* argv[])
{
  int i, natm, nline;
  double boxd[6];
  char iname[MAXNAME], oname[MAXNAME], line[MAXLINE], tag[MAXNAME];
  char** mfi;
  FILE *inp, *outp;

  /*** Options ***/
  if (argc == 1) {
    printf("\nChBox >> A program for changing the box dimensions of an AMBER "
	   "restart file.\n\nOptions:\n"
	   "  -c  : the original coordinate file (.RST format, coordinates "
	   "only)\n"
	   "  -o  : the output coordinate file (.RST format)\n"
           "  -al : box alpha angle\n"
           "  -bt : box beta angle\n"
           "  -gm : box gamma angle\n"
           "  -X  : X dimension of the box\n"
           "  -Y  : Y dimension of the box\n"
           "  -Z  : Z dimension of the box\n\n");
    exit(1);
  }

  iname[0] = '\0';
  oname[0] = '\0';
  boxd[0] = -1.0;
  boxd[1] = -1.0;
  boxd[2] = -1.0;
  boxd[3] = 90.0;
  boxd[4] = 90.0;
  boxd[5] = 90.0;
  for (i = 0; i < argc-2; i += 2) {
    strcpy(tag, *++argv);
    if (strcmp(tag, "-c") == 0) {
      strcpy(iname, *++argv);
    }
    else if (strcmp(tag, "-o") == 0) {
      strcpy(oname, *++argv);
    }
    else if (strcmp(tag, "-X") == 0) {
      boxd[0] = atof(*++argv);
    }
    else if (strcmp(tag, "-Y") == 0) {
      boxd[1] = atof(*++argv);
    }
    else if (strcmp(tag, "-Z") == 0) {
      boxd[2] = atof(*++argv);
    }
    else if (strcmp(tag, "-al") == 0) {
      boxd[3] = atof(*++argv);
    }
    else if (strcmp(tag, "-bt") == 0) {
      boxd[4] = atof(*++argv);
    }
    else if (strcmp(tag, "-gm") == 0) {
      boxd[5] = atof(*++argv);
    }
    else {
      printf("ChBox >> Error.  Identifier %s unkonwn.\n\n", tag);
      exit(1);
    }
  }

  /*** Check ***/
  if (iname[0] == '\0') {
    printf("ChBox >> Error.  Input file not specified.\n");
  }
  if (oname[0] == '\0') {
    printf("ChBox >> Error.  Output file not specified.\n");
  }

  /*** Start reading ***/
  if ((inp = fopen(iname, "r")) == NULL) {
    printf("ChBox >> Error.  File %s not found!\n", iname);
    exit(1);
  }
  fgets(line, 128, inp);
  fscanf(inp, "%d", &natm);
  if (natm < 0 || natm > 200000000) {
    printf("ChBox.c >> Error.  Silly number of atoms %d.\n", natm);
    exit(1);
  }
  nline = 3*(natm+1)/6;
  mfi = (char**)malloc(nline*sizeof(char*));
  fgets(line, MAXLINE, inp);
  for (i = 0; i < nline; i++) {
    fgets(line, MAXLINE, inp);
    mfi[i] = (char*)malloc(MAXLINE*sizeof(char));
    strcpy(mfi[i], line);
  }
  fclose(inp);

  /*** Start writing ***/
  outp = fopen(oname, "w");
  fprintf(outp, "    \n%6d\n", natm);
  for (i = 0; i < nline; i++) {
    fprintf(outp, "%s", mfi[i]);
  }
  for (i = 0; i < 6; i++) {
    fprintf(outp, "%12.7lf", boxd[i]);
  }
  fprintf(outp, "\n");
  fclose(outp);

  return 0;
}
