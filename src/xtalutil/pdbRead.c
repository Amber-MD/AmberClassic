#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "pdbRead.h"
#define PI                     3.1415926535898
#include "crdmanip.h"

/***=======================================================================***/
/*** GetPDB: function for taking in the information within a standard-     ***/
/***         format PDB file.  The information is returned within a pdb    ***/
/***         type structure.                                               ***/
/***=======================================================================***/
void GetPDB(pdb *thispdb, int verbosity, int pqr_format)
{
  int i, j, t, curr_res, rnumrescue, anumrescue;
  int aninum, anisoufound;
  char line[MAXLINE];
  FILE *pdbsource;

  /*** Scan the PDB file for atoms ***/
  if ((pdbsource = fopen(thispdb->source, "r")) == NULL) {
    printf("GetPDB >> Error.  PDB File %s not found!\n", thispdb->source);
    exit(1);
  }
  thispdb->n_atoms = 0;
  thispdb->n_ter = 0;
  while(fgets(line, MAXLINE, pdbsource) != NULL) {
    if ((line[0] == 'A' && strncmp(line, "ATOM", 4) == 0) ||
        (line[0] == 'H' && strncmp(line, "HETATM", 6) == 0)) {
      thispdb->n_atoms += 1;
    }
    else if (strncmp(line, "TER", 3) == 0) {
      thispdb->n_ter += 1;
    }
  }
  if (verbosity == 1) {
    printf("GetPDB >> Found %d atoms in %s.\n", thispdb->n_atoms, 
      thispdb->source);
  }
  rewind(pdbsource);

  /*** Allocate memory for the molecule ***/
  thispdb->res_nums = (int*)malloc(thispdb->n_atoms*sizeof(int));
  thispdb->atom_nums = (int*)malloc(thispdb->n_atoms*sizeof(int));
  thispdb->chain = (char*)malloc(thispdb->n_atoms*sizeof(char));
  thispdb->icode = (char*)malloc(thispdb->n_atoms*sizeof(char));
  thispdb->alt = (char*)malloc(thispdb->n_atoms*sizeof(char));
  thispdb->element = (char*)malloc(2*thispdb->n_atoms*sizeof(char));
  thispdb->atom_names = (char*)malloc(4*thispdb->n_atoms*sizeof(char));
  thispdb->res_names = (char*)malloc(4*thispdb->n_atoms*sizeof(char));
  thispdb->crds = (double*)calloc(3*thispdb->n_atoms, sizeof(double));
  thispdb->termini = (int*)malloc(thispdb->n_ter*sizeof(int));
  thispdb->anisou = (int*)calloc(6*thispdb->n_atoms, sizeof(int));
  thispdb->pqr = pqr_format;
  thispdb->aug = 0;
  thispdb->AnisouRec = 0;
  thispdb->charges = (double*)malloc(thispdb->n_atoms*sizeof(double));
  thispdb->radii = (double*)malloc(thispdb->n_atoms*sizeof(double));

  /*** Read in the PDB information ***/
  i = 0;
  t = 0;
  rnumrescue = 0;
  anumrescue = 0;
  while(fgets(line, MAXLINE, pdbsource) != NULL) {
    if ((line[0] == 'A' && strncmp(line, "ATOM", 4) == 0) ||
        (line[0] == 'H' && strncmp(line, "HETATM", 6) == 0)) {
      if (line[6] != '*') {
        sscanf(&line[6], "%d", &thispdb->atom_nums[i]);
        anumrescue = thispdb->atom_nums[i];
      }
      else {
        thispdb->atom_nums[i] = anumrescue;
        anumrescue++;
      }
      for (j = 0; j < 4; j++) {
        thispdb->atom_names[i*4+j] = line[j+12];
        thispdb->res_names[i*4+j] = line[j+17];
      }
      if (line[22] != '*') {
        if (thispdb->pqr == 0) {
          sscanf(&line[16], "%1c",&thispdb->alt[i]);
          sscanf(&line[21], "%1c",&thispdb->chain[i]);
          sscanf(&line[22], "%4d",&thispdb->res_nums[i]);
          sscanf(&line[26], "%1c",&thispdb->icode[i]);
          sscanf(&line[30], "%8lf",&thispdb->crds[i*3]);
          sscanf(&line[38], "%8lf",&thispdb->crds[i*3+1]);
          sscanf(&line[46], "%8lf",&thispdb->crds[i*3+2]);
          sscanf(&line[54], "%6lf",&thispdb->charges[i]);
          sscanf(&line[60], "%6lf",&thispdb->radii[i]);
          sscanf(&line[76], "%1c",&thispdb->element[i*2]);
          sscanf(&line[77], "%1c",&thispdb->element[i*2+1]);
        } else {
          sscanf(&line[21], "%c%d%lf%lf%8lf%8lf%8lf", &thispdb->chain[i],
                 &thispdb->res_nums[i], &thispdb->crds[i*3],
                 &thispdb->crds[i*3+1], &thispdb->crds[i*3+2],
                 &thispdb->charges[i], &thispdb->radii[i]);
       }
       rnumrescue = thispdb->res_nums[i]+1;
      } else {
       thispdb->alt[i] = line[16];
       thispdb->chain[i] = line[21];
        if (thispdb->pqr == 0) {
          sscanf(&line[30], "%8lf",&thispdb->crds[i*3]);
          sscanf(&line[38], "%8lf",&thispdb->crds[i*3+1]);
          sscanf(&line[46], "%8lf",&thispdb->crds[i*3+2]);
        } 
        else
          sscanf(&line[30], "%lf%lf%lf%lf%lf", &thispdb->crds[i*3],
                 &thispdb->crds[i*3+1], &thispdb->crds[i*3+2],
                &thispdb->charges[i], &thispdb->radii[i]);
        thispdb->res_nums[i] = rnumrescue;
      }

      /*** Don't filter the character inputs: assume input PDB is correct ***/
      // touchup(&thispdb->atom_names[i*4], 4);
      // touchup(&thispdb->res_names[i*4], 4);
      // alphanumeric(&thispdb->chain[i]);

      i++;
    }
    else if (line[0] == 'A' && strncmp(line, "ANISOU", 6) == 0) {
      thispdb->AnisouRec = 1;
      if (line[6] != '*') {
        sscanf(&line[6], "%d", &aninum);
        anisoufound = 0;
        for (j = 0; j < i; j++) {
          if (aninum == thispdb->atom_nums[j]) {
            sscanf(&line[28], "%d%d%d%d%d%d", &thispdb->anisou[6*j],
                   &thispdb->anisou[6*j+1], &thispdb->anisou[6*j+2],
                   &thispdb->anisou[6*j+3], &thispdb->anisou[6*j+4],
                   &thispdb->anisou[6*j+5]);
            anisoufound = 1;
          }
        }
        if (anisoufound == 0) {
          printf("GetPDB >> Error.  Anisotropy records for atom %d cannot be "
                 "parsed.\nGetPDB >> Exiting.\n", aninum);
          exit(1);
        }
      }
      else {
        printf("GetPDB >> Error.  Anisotropy records cannot contain **** "
               "numerical entries.\nGetPDB >> Exiting.\n");
        exit(1);
      }
    }
    else if (strncmp(line, "TER", 3) == 0) {
      thispdb->termini[t] = i;
      t++;
      rnumrescue++;
    }
  }
  fclose(pdbsource);

  /*** Determine the residue limits and normalize residue indices ***/
  thispdb->n_res = 1;
  curr_res = thispdb->res_nums[0];
  j = 0;
  for (i = 0; i < thispdb->n_atoms; i++) {
    if (thispdb->res_nums[i] != curr_res) {
      thispdb->n_res += 1;
      curr_res = thispdb->res_nums[i];
      j++;
    }
    thispdb->res_nums[i] = j;
  }
  thispdb->res_lims = (int*)malloc((thispdb->n_res+1)*sizeof(int));
  curr_res = thispdb->res_nums[0];
  thispdb->res_lims[0] = 0;
  j = 1;
  for (i = 0; i < thispdb->n_atoms; i++) {
    if (thispdb->res_nums[i] != curr_res) {
      curr_res = thispdb->res_nums[i];
      thispdb->res_lims[j] = i;
      j++;
    }
  }
  thispdb->res_lims[thispdb->n_res] = thispdb->n_atoms;

  /*** Make sure TER cards fall in between residues ***/
  for (i = 0; i < thispdb->n_ter; i++) {
    t = 0;
    for (j = i; j < thispdb->n_res+1; j++) {
      if (thispdb->res_lims[j] == thispdb->termini[i]) {
        t = 1;
        break;
      }
    }
    if (t == 0 && verbosity == 1) {
      printf("GetPDB >> Warning.  TER card after atom %d does not match "
             "residue limits!\n", thispdb->termini[i]);
    }
  }

  /*** No TER card after the last residue ***/
  if (thispdb->n_ter > 0 &&
      thispdb->termini[thispdb->n_ter-1] == thispdb->n_atoms) {
    thispdb->n_ter -= 1;
  }

  if (verbosity == 1) {
    printf("GetPDB >> Bracketed %d residues.\n\n", thispdb->n_res);
  }
}

/***=======================================================================***/
/*** PutPDB: function for writing PDB files after Unitcell and AddToBox;   ***/
/***         improved version of the original and a companion to GetPDB.   ***/
/***=======================================================================***/
void PutPDB(pdb *thispdb, symT *S, char* filename, char* pdb_form,
            char* custom_header, int verbosity)
{
  int i, j;
  int* itmp;
  double* dtmp;
  char* ctmp;
  char* ctm2p;
  char line[MAXLINE];
  FILE *pdbtxt;
  
  S->boxd[3]*=180/PI;
  S->boxd[4]*=180/PI;
  S->boxd[5]*=180/PI;

  pdbtxt = fopen(filename, "w");
  fprintf(pdbtxt, "HEADER  Generated by AmberTools\n%s\n", custom_header);
  fprintf(pdbtxt, "CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n", S->boxd[0], S->boxd[1],S->boxd[2], S->boxd[3], S->boxd[4], S->boxd[5]);
  for (i = 0; i < thispdb->n_atoms; i++) {
    ctmp = &thispdb->atom_names[i*4];
    ctm2p = &thispdb->res_names[i*4];
    dtmp = &thispdb->crds[i*3];
    itmp = &thispdb->anisou[i*6];

    /*** The basic information ***/
    if (thispdb->pqr == 0) {
#if 0
      printf( "\n%.4s\n", ctmp );
      printf( "%c\n", thispdb->alt[i] );
      printf( "%.4s\n", ctm2p );
      printf( "%6.2lf\n", thispdb->charges[i] );
      printf( "%6.2lf\n", thispdb->radii[i] );
      printf( "%c\n", thispdb->element[2*i] );
      printf( "%c\n", thispdb->element[2*i+1] );
      exit(1);
#endif
      sprintf(line, 
   "ATOM %6d %.4s%c%.4s%c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %c%c\n",
          (thispdb->atom_nums[i]+1) % 100000, ctmp, 
          thispdb->alt[i], ctm2p, thispdb->chain[i],
          (thispdb->res_nums[i]+1) % 10000, 
          dtmp[0], dtmp[1], dtmp[2],
          thispdb->charges[i], thispdb->radii[i], 
          thispdb->element[2*i], thispdb->element[2*i+1]);
    }
    else {
      sprintf(line, "ATOM %6d %.4s %.4s %4d    %8.3lf%8.3lf%8.3lf"
       "%6.2lf %6.2lf\n", (thispdb->atom_nums[i]+1) % 100000, ctmp, ctm2p,
       (thispdb->res_nums[i]+1) % 10000, dtmp[0], dtmp[1], dtmp[2],
       thispdb->charges[i], thispdb->radii[i]);
    }
    // realign_name(&line[12]);
    fprintf(pdbtxt, "%s", line);

    if (thispdb->AnisouRec == 1) {
   sprintf(line, "ANISOU%5d %.4s%c%.4s%c%4d  %7d%7d%7d%7d%7d%7d      %c%c\n",
        (thispdb->atom_nums[i]+1) % 100000, ctmp, 
        thispdb->alt[i], ctm2p, thispdb->chain[i], thispdb->res_nums[i]+1,
        itmp[0], itmp[1], itmp[2], itmp[3], itmp[4], itmp[5],
        thispdb->element[2*i], thispdb->element[2*i+1]);
        // realign_name(&line[12]);
        fprintf(pdbtxt, "%s", line);
    }

    /*** Print TER card ***/
    for (j = 0; j < thispdb->n_ter; j++) {
      if (thispdb->termini[j] == i+1) {
        fprintf(pdbtxt, "TER\n");
        break;
      }
      if (thispdb->termini[j] > i+1) {
        break;
      }
    }
  }
  fprintf(pdbtxt, "END\n");
  fclose(pdbtxt);
  if (verbosity == 1) {
    printf("PutPDB >> Molecule printed to file %s in %s format\n\n",
           filename, pdb_form);
  }
}

/***=======================================================================***/
/*** PutPDBPro: function for writing PDB files after running PropPDB, an   ***/
/***         improved version of the original and a companion to GetPDB.   ***/
/***=======================================================================***/
void PutPDBPro(pdb *thispdb, symT *S, int* irep, char* filename, char* pdb_form,
            char* custom_header, int verbosity)
{
  int i, j;
  int* itmp;
  double* dtmp;
  char* ctmp;
  char* ctm2p;
  char line[MAXLINE];
  FILE *pdbtxt;
  
  S->boxd[3]*=180/PI;
  S->boxd[4]*=180/PI;
  S->boxd[5]*=180/PI;

  pdbtxt = fopen(filename, "w");
  fprintf(pdbtxt, "HEADER  Generated by AmberTools\n%s\n", custom_header);
  fprintf(pdbtxt, "CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n", S->boxd[0]*irep[0], S->boxd[1]*irep[1],S->boxd[2]*irep[2], S->boxd[3], S->boxd[4], S->boxd[5]);
  for (i = 0; i < thispdb->n_atoms; i++) {
    ctmp = &thispdb->atom_names[i*4];
    ctm2p = &thispdb->res_names[i*4];
    dtmp = &thispdb->crds[i*3];
    itmp = &thispdb->anisou[i*6];

    /*** The basic information ***/
    if (thispdb->pqr == 0) {
      sprintf(line, 
   "ATOM %6d %.4s%c%.4s%c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %c%c\n",
          (thispdb->atom_nums[i]+1) % 100000, ctmp, 
          thispdb->alt[i], ctm2p, thispdb->chain[i],
          (thispdb->res_nums[i]+1) % 10000, dtmp[0], dtmp[1], dtmp[2],
          thispdb->charges[i], thispdb->radii[i], 
          thispdb->element[2*i], thispdb->element[2*i+1]);
    }
    else {
      sprintf(line, "ATOM %6d %.4s %.4s %4d    %8.3lf%8.3lf%8.3lf"
        "%6.2lf %6.2lf\n", (thispdb->atom_nums[i]+1) % 100000, ctmp, ctm2p,
        (thispdb->res_nums[i]+1) & 10000, dtmp[0], dtmp[1], dtmp[2],
        thispdb->charges[i], thispdb->radii[i]);
    }
    // realign_name(&line[12]);
    fprintf(pdbtxt, "%s", line);
    if (thispdb->AnisouRec == 1) {
   sprintf(line, "ANISOU%5d %.4s%c%.4s%c%4d  %7d%7d%7d%7d%7d%7d      %c%c\n",
        (thispdb->atom_nums[i]+1) % 100000, ctmp, 
        thispdb->alt[i], ctm2p, thispdb->chain[i], thispdb->res_nums[i]+1,
        itmp[0], itmp[1], itmp[2], itmp[3], itmp[4], itmp[5],
        thispdb->element[2*i], thispdb->element[2*i+1]);
        // realign_name(&line[12]);
        fprintf(pdbtxt, "%s", line);
    }

    /*** Print TER card ***/
    for (j = 0; j < thispdb->n_ter; j++) {
      if (thispdb->termini[j] == i+1) {
        fprintf(pdbtxt, "TER\n");
        break;
      }
      if (thispdb->termini[j] > i+1) {
        break;
      }
    }
  }
  fprintf(pdbtxt, "END\n");
  fclose(pdbtxt);
  if (verbosity == 1) {
    printf("PutPDB >> Molecule printed to file %s in %s format\n\n",
           filename, pdb_form);
  }
}

/***=======================================================================***/
/*** FreePDB: function for freeing a pdb structure.  Easy as pie!          ***/
/***=======================================================================***/
void FreePDB(pdb *thispdb)
{
  free(thispdb->res_nums);
  free(thispdb->atom_nums);
  free(thispdb->res_lims);
  free(thispdb->chain);
  free(thispdb->atom_names);
  free(thispdb->res_names);
  free(thispdb->crds);
  free(thispdb->termini);
  free(thispdb->anisou);
  if (thispdb->pqr == 1) {
    free(thispdb->charges);
    free(thispdb->radii);
  }
  if (thispdb->aug == 1) {
    free(thispdb->z);
    free(thispdb->mass);
  }
  thispdb->pqr = 0;
  thispdb->aug = 0;
  thispdb->AnisouRec = 0;
}

/***=======================================================================***/
/*** SSCAN_ADV: scans a string from a larger string, then advances to the  ***/
/***            last character scanned in the larger line.  Ignores any    ***/
/***            whitespace.                                                ***/
/***=======================================================================***/
char* sscan_adv(char *ctmp, char* cword)
{
  while (ctmp[0] == ' ') {
    ctmp = &ctmp[1];
  }
  sscanf(ctmp, "%s", cword);
  ctmp = &ctmp[strlen(cword)];

  return ctmp;
}

/***=======================================================================***/
/*** TOUCHUP: combines the functions append_whitespace, remove_whitespace, ***/
/***          and alphanumeric into one convenient package.                ***/
/***=======================================================================***/
void touchup(char* a, int asize)
{
  int i;

  append_whitespace(a, asize);
  remove_whitespace(a, asize);
  for (i = 0; i < asize; i++) {
    if (a[i] == '+' || a[i] == '-') {
      continue;
    }
    alphanumeric(&a[i]);
  }
}

/***=======================================================================***/
/*** REMOVE_WHITESPACE: function for removing white space at the beginning ***/
/***                    of a string.                                       ***/
/***=======================================================================***/
void remove_whitespace(char* a, int asize)
{
  int i, j, num_rem;

  num_rem = 0;
  while (a[num_rem] == ' ' && num_rem < asize) {
    num_rem++;
  }
  j = 0;
  for (i = num_rem; i < asize; i++) {
    a[j] = a[i];
    j++;
  }
  for (i = asize-num_rem; i < asize; i++) {
    a[i] = ' ';
  }
}

/***=======================================================================***/
/*** APPEND_WHITESPACE: function for adding whitespace to the end of a     ***/
/***                    string, staring at the (first) instance of '\0'    ***/
/***=======================================================================***/
void append_whitespace(char* a, int asize)
{
  int i, whiteout;
  
  whiteout = 0;
  for (i = 0; i < asize; i++) {
    if (a[i] == '\0') {
      whiteout = 1;
    }
    if (whiteout == 1) {
      a[i] = ' ';
    }
  }
}

/***=======================================================================***/
/*** ALPHANUMERIC: make sure a character is alphanumeric; if not, make it  ***/
/***               whitespace.                                             ***/
/***=======================================================================***/
void alphanumeric(char *a)
{
  if (isalphanum(*a) == 0) {
    *a = ' ';
  }
}

/***=======================================================================***/
/*** ISALPHANUM: return 0 if a character is not 'alphanumeric'             ***/
/***=======================================================================***/
int isalphanum(char a)
{
  if ((a < 'a' || a > 'z') && (a < 'A' || a > 'Z') && (a < '0' || a > '9') &&
      (a != '*') && (a != '\'')) {
    return 0;
  }
  else {
    return 1;
  }
}

/***=======================================================================***/
/*** REALIGN_NAME: function for realigning an atom name of four characters ***/
/***               to look like the original PDB format (one space goes in ***/
/***               front, unless the first character is a letter or the    ***/
/***               name has four characters).                              ***/
/***=======================================================================***/
void realign_name(char* atom_name)
{
  int i;

  if (atom_name[0] != ' ' && atom_name[3] == ' ' &&
      (atom_name[0] < '0' || atom_name[0] > '9')) {
    for (i = 3; i > 0; i--) {
      atom_name[i] = atom_name[i-1];
    }
    atom_name[0] = ' ';
  }
}

/***=======================================================================***/
/*** CopyAtoms: copies atoms from one PDB to another.                      ***/
/***=======================================================================***/
void CopyAtoms(pdb *add, int astart, int aend, pdb *base, int bstart, int chn)
{
  int i, j, ib;

  for (i = astart; i < aend; i++) {
    ib = (i-astart) + bstart;
    for (j = 0; j < 4; j++) {
      base->atom_names[4*ib+j] = add->atom_names[4*i+j];
      base->res_names[4*ib+j] = add->res_names[4*i+j];
    }
    for (j = 0; j < 2; j++) {
      base->element[2*ib+j] = add->element[2*i+j];
    }
    for (j = 0; j < 3; j++) {
      base->crds[3*ib+j] = add->crds[3*i+j];
    }
    base->charges[ib] = add->charges[i];
    base->radii[ib] = add->radii[i];
    base->atom_nums[ib] = add->atom_nums[i];
    base->res_nums[ib] = add->res_nums[i];
    base->chain[ib] = (chn > 0) ? chn : add->chain[i];
    base->alt[ib] = add->alt[i];
  }
}

/***=======================================================================***/
/*** DeleteResidue: delete a residue from a pdb structure, moving termini  ***/
/***                and residue limits as appropriate.                     ***/
/***=======================================================================***/
void DeleteResidue(pdb *mol, int r)
{
  int i, j, i3, i4, rlen, r3len, r4len;

  /*** Advance all atoms ***/
  rlen = mol->res_lims[r+1] - mol->res_lims[r];
  r4len = 4*rlen;
  r3len = 3*rlen;
  for (i = mol->res_lims[r+1]; i < mol->n_atoms; i++) {
    i3 = 3*i;
    i4 = 4*i;
    for (j = 0; j < 4; j++) {
      mol->atom_names[i4-r4len+j] = mol->atom_names[i4+j];
      mol->res_names[i4-r4len+j] = mol->res_names[i4+j];
    }
    for (j = 0; j < 3; j++) {
      mol->crds[i3-r3len+j] = mol->crds[i3+j];
    }
    mol->atom_nums[i-rlen] = mol->atom_nums[i];
    mol->res_nums[i-rlen] = mol->res_nums[i];
    mol->chain[i-rlen] = mol->chain[i];
  }

  /*** Update termini ***/
  for (i = 0; i < mol->n_ter; i++) {
    if (mol->termini[i] > mol->res_lims[r] &&
        mol->termini[i] <= mol->res_lims[r+1]) {
      for (j = i+1; j < mol->n_ter; j++) {
        mol->termini[j-1] = mol->termini[j];
      }
      mol->n_ter -= 1;
    }
    if (mol->termini[i] >= mol->res_lims[r+1]) {
      mol->termini[i] -= rlen;
    }
  }

  /*** Update residue limits ***/
  for (i = r+1; i < mol->n_res-1; i++) { 
    mol->res_lims[i] = mol->res_lims[i+1] - rlen;
  }

  /*** Update the numbers of atoms and residues ***/
  mol->n_atoms -= rlen;
  mol->n_res -= 1;
}

/***=======================================================================***/
/*** ModPdbRA: modify the numbers of residues and atoms so that they do    ***/
/***           not overload the fields.                                    ***/
/***=======================================================================***/
void ModPdbRA(pdb *w)
{
  int i, j;

  for (i = 0; i < w->n_res; i++) {
    for (j = w->res_lims[i]; j < w->res_lims[i+1]; j++) {
      w->res_nums[j] = i;
    }
  }
  for (i = 0; i < w->n_atoms; i++) {
    w->atom_nums[i] = i;
  }
  for (i = 0; i < w->n_atoms; i++) {
    w->atom_nums[i] = w->atom_nums[i] % 1000000;
    w->res_nums[i] = w->res_nums[i] % 10000;
  }
}

/***=======================================================================***/
/*** AugmentPDB: augment a PDB structure with additional information.      ***/
/***=======================================================================***/
void AugmentPDB(pdb *thispdb)
{
  int i, j, k;

  thispdb->aug = 1;
  thispdb->z = (int*)calloc(thispdb->n_atoms, sizeof(int));
  thispdb->mass = (double*)calloc(thispdb->n_atoms, sizeof(double));

  /*** Atomic numbers ***/
  for (i = 0; i < thispdb->n_atoms; i++) {
    j = 4*i;
    k = thispdb->atom_names[j];
    while(k > 47 && k < 58) {
      j++;
      k = thispdb->atom_names[j];
    }
    if (k == 72) {
      thispdb->z[i] = 1;
      thispdb->mass[i] = 1.0078;
    }
    if (k == 67) {
      thispdb->z[i] = 6;
      thispdb->mass[i] = 12.0000;
    }
    if (k == 79) {
      thispdb->z[i] = 8;
      thispdb->mass[i] = 15.9996;
    }
    if (k == 78) {
      thispdb->z[i] = 7;
      thispdb->mass[i] = 14.0067;
    }
    if (k == 83) {
      thispdb->z[i] = 16;
      thispdb->mass[i] = 32.0650;
    }
  }
}

/***=======================================================================***/
/*** LoadSymmetryData: load symmetry information from a PDB file.          ***/
/***=======================================================================***/
symT LoadSymmetryData(char* source)
{
  int i, nm;
  double m1, m2, m3, t;
  char line[128], cword[128];
  FILE *inp;
  symT S;
  int is_CRYST1;

  is_CRYST1 = 0;
  if ((inp = fopen(source, "r")) == NULL) {
    printf("GetPDB >> Error.  PDB File %s not found!\n", source);
    exit(1);
  }
  S.nsym = 0;
  while(fgets(line, 128, inp) != NULL) {
    if (strncmp(line, "REMARK", 6) == 0) {
      sscanf(line, "%s%s%s", cword, cword, cword);
      if (strncmp(cword, "SMTRY1", 6) == 0) {
        S.nsym += 1;
      }
    }
  }
  S.tmat = (dmat*)malloc(S.nsym*sizeof(dmat));
  for (i = 0; i < S.nsym; i++) {
    S.tmat[i] = CreateDmat(3, 3);
  }
  S.tvec = (double*)malloc(3*S.nsym*sizeof(double));
  rewind(inp);
  while(fgets(line, 128, inp) != NULL) {
    if (strncmp(line, "REMARK", 6) == 0) {
      sscanf(line, "%s%s%s", cword, cword, cword);
      if (strncmp(cword, "SMTRY", 5) == 0) {
        sscanf(line, "%s%s%s%d%lf%lf%lf%lf\n", cword, cword, cword, &nm, &m1,
               &m2, &m3, &t);
        cword[6] = '\0';
        i = atoi(&cword[5])-1;
        nm -= 1;
        S.tmat[nm].data[3*i] = m1;
        S.tmat[nm].data[3*i+1] = m2;
        S.tmat[nm].data[3*i+2] = m3;
        S.tvec[3*nm+i] = t;
      }
    }
    if (strncmp(line, "CRYST1", 6) == 0) {
      sscanf(&line[6], "%lf%lf%lf%lf%lf%lf\n", &S.boxd[0], &S.boxd[1],
                        &S.boxd[2], &S.boxd[3], &S.boxd[4], &S.boxd[5]);
      for (i = 3; i < 6; i++) {
        S.boxd[i] *= PI/180.0;
      }
      is_CRYST1 = 1;
    }
  }
  if (is_CRYST1 == 0) {
    printf("pdbRead >> Error. PDB file must contain CRYST1 card!\n");
    exit(1);
  }
  fclose(inp);
  
  return S;
}

/***=======================================================================***/
/*** ExpandPDBSym: expand a PDB file based on the symmetry operations of   ***/
/***               the crystal space group.                                ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   p:       the pdb struct to expand                                   ***/
/***   S:       the symmetry information about the crystal space group     ***/
/***   stride:  flag to activate atom and residue number offsets in the    ***/
/***            new pdb                                                    ***/
/***=======================================================================***/
pdb ExpandPDBSym(pdb *p, symT *S, int stride)
{
  int h, i, m, n, rskip, tm;
  double x, y, z;
  double *dtmp;
  pdb t;

  /*** Now, make the new PDB ***/
  rskip = p->n_res;
  t.n_atoms = S->nsym*p->n_atoms;
  t.n_res = S->nsym*p->n_res;
  t.n_ter = S->nsym*(p->n_ter+1);
  t.aug = 0;
  t.pqr = 0;
  t.AnisouRec = p->AnisouRec;
  t.res_nums = (int*)malloc(t.n_atoms*sizeof(int));
  t.atom_nums = (int*)malloc(t.n_atoms*sizeof(int));
  t.res_lims = (int*)malloc((t.n_res+1)*sizeof(int));
  t.termini = (int*)malloc(t.n_ter*sizeof(int));
  t.alt = (char*)malloc(t.n_atoms*sizeof(char));
  t.element = (char*)malloc(2*t.n_atoms*sizeof(char));
  t.chain = (char*)malloc(t.n_atoms*sizeof(char));
  t.atom_names = (char*)malloc(4*t.n_atoms*sizeof(char));
  t.res_names = (char*)malloc(4*t.n_atoms*sizeof(char));
  t.crds = (double*)malloc(3*t.n_atoms*sizeof(double));
  t.charges = (double*)malloc(t.n_atoms*sizeof(double));
  t.radii = (double*)malloc(t.n_atoms*sizeof(double));
  t.anisou = (int*)malloc(6*t.n_atoms*sizeof(int));
  

  /*** Construct symmetry-related copies ***/
  h = 0;
  for (i = 0; i < S->nsym; i++) {
    for (m = 0; m < p->n_atoms; m++) {
      tm = m + h*p->n_atoms;
      if (stride == 0) {
        t.res_nums[tm] = p->res_nums[m];
        t.atom_nums[tm] = p->atom_nums[m];
      }
      else {
        t.res_nums[tm] = p->res_nums[m] + h*rskip;
        t.atom_nums[tm] = p->atom_nums[m] + h*p->n_atoms;
      }
      t.alt[tm] = p->alt[m];
      t.chain[tm] = p->chain[m];
      t.element[2*tm] = p->element[2*m];
      t.element[2*tm+1] = p->element[2*m+1];
      for (n = 0; n < 4; n++) {
        t.res_names[4*tm+n] = p->res_names[4*m+n];
        t.atom_names[4*tm+n] = p->atom_names[4*m+n];
      }
      dtmp = S->tmat[i].data;
      x = p->crds[3*m];
      y = p->crds[3*m+1];
      z = p->crds[3*m+2];
      t.crds[3*tm] = dtmp[0]*x + dtmp[1]*y + dtmp[2]*z + S->tvec[3*i];
      t.crds[3*tm+1] = dtmp[3]*x + dtmp[4]*y + dtmp[5]*z + S->tvec[3*i+1];
      t.crds[3*tm+2] = dtmp[6]*x + dtmp[7]*y + dtmp[8]*z + S->tvec[3*i+2];
      for (n = 0; n < 6; n++) {
        t.anisou[6*tm+n] = p->anisou[6*m+n];
      }
      t.charges[tm] = p->charges[m];
      t.radii[tm] = p->radii[m];
    }
    for (m = 0; m < p->n_res; m++) {
      tm = m + h*p->n_res;
      t.res_lims[tm] = p->res_lims[m] + h*p->n_atoms;
    }
    for (m = 0; m < p->n_ter; m++) {
      tm = m + h*(p->n_ter+1);
      t.termini[tm] = p->termini[m] + h*p->n_atoms;
    }
    t.termini[m + h*(p->n_ter+1)] = (h+1)*p->n_atoms;
    h++;
  }
  t.res_lims[t.n_res] = t.n_atoms;

  return t;
}

/***=======================================================================***/
/*** TilePDB: expand a PDB file by tiling the unit cell in space.          ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   p:       the pdb struct to expand                                   ***/
/***   boxd:    6-element vector of the box dimensions                     ***/
/***   irep:    the number of tiles along the three box vectors            ***/
/***   stride:  flag to activate atom and residue number offsets in the    ***/
/***            new pdb                                                    ***/
/***=======================================================================***/
pdb TilePDB(pdb *p, symT *S, int* irep, int stride)
{
  int h, i, j, k, m, n, nrep, rskip, tm;
  dmat U, invU;
  pdb t;

  /*** get boxd from symT S, Mei ***/
  double* boxd;
  boxd = &S->boxd[0];

  /*** Now, make the new PDB ***/
  nrep = irep[0]*irep[1]*irep[2];
  rskip = p->n_res;
  t.n_atoms = nrep*p->n_atoms;
  t.n_res = nrep*p->n_res;
  t.n_ter = nrep*(p->n_ter+1);
  t.aug = 0;
  t.pqr = 0;
  t.AnisouRec = p->AnisouRec;
  t.res_nums = (int*)malloc(t.n_atoms*sizeof(int));
  t.atom_nums = (int*)malloc(t.n_atoms*sizeof(int));
  t.res_lims = (int*)malloc((t.n_res+1)*sizeof(int));
  t.termini = (int*)malloc(t.n_ter*sizeof(int));
  t.alt = (char*)malloc(t.n_atoms*sizeof(char));
  t.element = (char*)malloc(2*t.n_atoms*sizeof(char));
  t.chain = (char*)malloc(t.n_atoms*sizeof(char));
  t.atom_names = (char*)malloc(4*t.n_atoms*sizeof(char));
  t.res_names = (char*)malloc(4*t.n_atoms*sizeof(char));
  t.crds = (double*)malloc(3*t.n_atoms*sizeof(double));
  t.charges = (double*)malloc(t.n_atoms*sizeof(double));
  t.radii = (double*)malloc(t.n_atoms*sizeof(double));
  t.anisou = (int*)malloc(6*t.n_atoms*sizeof(int));

  /*** Go into box space ***/
  U = CreateDmat(3, 3);
  invU = CreateDmat(3, 3);
  CmpXfrm(boxd, U, invU);
  RotateCrd(p->crds, p->n_atoms, U);

  /*** Propagate ***/
  for (m = 0; m < p->n_atoms; m++) {
    t.res_nums[m] = p->res_nums[m];
    t.atom_nums[m] = p->atom_nums[m];
    t.chain[m] = p->chain[m];
    t.alt[m] = p->alt[m];
    t.element[2*m] = p->element[2*m];
    t.element[2*m+1] = p->element[2*m+1];
    for (n = 0; n < 4; n++) {
      t.res_names[4*m+n] = p->res_names[4*m+n];
      t.atom_names[4*m+n] = p->atom_names[4*m+n];
    }
    for (n = 0; n < 3; n++) {
      t.crds[3*m+n] = p->crds[3*m+n];
    }
    t.charges[m] = p->charges[m];
    t.radii[m] = p->radii[m];
  }
  for (m = 0; m < p->n_res; m++) {
    t.res_lims[m] = p->res_lims[m];
  }
  for (m = 0; m < p->n_ter; m++) {
    t.termini[m] = p->termini[m];
  }
  t.termini[m] = p->n_atoms;

  h = 0;
  for (i = 0; i < irep[0]; i++) {
    for (j = 0; j < irep[1]; j++) {
      for (k = 0; k < irep[2]; k++) {
        for (m = 0; m < p->n_atoms; m++) {
          tm = m + h*p->n_atoms;
          t.res_nums[tm] = p->res_nums[m] + h*rskip;
          t.atom_nums[tm] = p->atom_nums[m] + h*p->n_atoms;
          t.chain[tm] = p->chain[m];
          t.alt[tm] = p->alt[m];
          t.element[2*tm] = p->element[2*m];
          t.element[2*tm+1] = p->element[2*m+1];
          for (n = 0; n < 4; n++) {
            t.res_names[4*tm+n] = p->res_names[4*m+n];
            t.atom_names[4*tm+n] = p->atom_names[4*m+n];
          }
          t.crds[3*tm] = p->crds[3*m] + i;
          t.crds[3*tm+1] = p->crds[3*m+1] + j;
          t.crds[3*tm+2] = p->crds[3*m+2] + k;
          for (n = 0; n < 6; n++) {
            t.anisou[6*tm+n] = p->anisou[6*m+n];
          }
          t.charges[tm] = p->charges[m];
          t.radii[tm] = p->radii[m];
        }
        for (m = 0; m < p->n_res; m++) {
          tm = m + h*p->n_res;
          t.res_lims[tm] = p->res_lims[m] + h*p->n_atoms;
        }
        for (m = 0; m < p->n_ter; m++) {
          tm = m + h*(p->n_ter+1);
          t.termini[tm] = p->termini[m] + h*p->n_atoms;
        }
        t.termini[m + h*(p->n_ter+1)] = (h+1)*p->n_atoms;
        h++;
      }
    }
  }
  t.res_lims[t.n_res] = t.n_atoms;

  /*** Back to real space ***/
  RotateCrd(t.crds, t.n_atoms, invU);

  return t;
}

