#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "Topology.h"
#include "mdgxVector.h"
#include "Matrix.h"
#include "Grid.h"
#include "CrdManip.h"
#include "Constants.h"
#include "Macros.h"
#include "CompFrc.h"
#include "Manual.h"
#include "Parse.h"
#include "VirtualSites.h"
#include "Buckingham.h"
#include "BroadcastCommand.h"
#include "Debug.h"
#include "Restraints.h"
#include "ParamFit.h"
 
#include "pmeDirectDS.h"
#include "TrajectoryDS.h"

//-----------------------------------------------------------------------------
// RoundTotalCharge: force the total system charge to be an integer value if
//                   it is close to an integer multiple of the proton charge
//                   to begin with.                       
//                                                                      
// Arguments:                                                           
//   tp:     the system topology                                        
//-----------------------------------------------------------------------------
static void RoundTotalCharge(prmtop *tp)
{
  int i, nwithq;
  double tq;

  // Compute the sum of all charges
  tq = DSum(tp->Charges, tp->natom);
  tp->initq = tq;
  if (tq - floor(tq) > 0.1 && ceil(tq) - tq > 0.1 &&
      fabs(tq - floor(tq)) > 1.0e-8 && fabs(ceil(tq) - tq) > 1.0e-8) {
    return;
  }

  // If we're still here, the charge is close to an
  // integer value and so should be neutralized.   
  // The excess charge is added only to atoms whose
  // charges are not close to zero, to ensure that 
  // atoms with zero charge can continue to be a   
  // speed advantage for the dynamics routines.    
  tq = (tq - floor(tq) < 0.1) ? floor(tq) - tq : ceil(tq) - tq;
  nwithq = 0;
  for (i = 0; i < tp->natom; i++) {
    if (fabs(tp->Charges[i]) >= 1.0e-8) {
      nwithq++;
    }
  }
  tq /= nwithq;
  for (i = 0; i < tp->natom; i++) {
    if (fabs(tp->Charges[i]) >= 1.0e-8) {
      tp->Charges[i] += tq;
    }
  }
}

//-----------------------------------------------------------------------------
// RecursiveBondSeek: recursive routine for seeking out, and labeling, all
//                    atoms connected together by bonds.  This routine contains
//                    an internal check to ensure that atoms do not belong to
//                    multiple groups.           
//                                                                      
// Arguments:                                                           
//   oatm:   the origin atom                                            
//   labl:   the number with which to label all atoms connected to oatm 
//   glist:  the groups list, storing labels for all atoms              
//   tp:     the system topology                                        
//-----------------------------------------------------------------------------
static void RecursiveBondSeek(int oatm, int labl, int* glist, prmtop *tp)
{
  int i, xatm;

  // If the atom is already counted in this group, return
  if (glist[oatm] == labl) {
    return;
  }

  // If the atom is already in another group, we're in trouble
  if (glist[oatm] >= 0) {
    printf("RecursiveBondSeek >> Error.  Atom %d is in more than one "
           "group!\n", oatm);
    exit(1);
  }
  glist[oatm] = labl;
  for (i = 0; i < tp->nb1234[oatm].n11; i++) {
    xatm = tp->nb1234[oatm].L11[i];
    glist[xatm] = labl;
    RecursiveBondSeek(xatm, labl, glist, tp);
  }
  for (i = 0; i < tp->nb1234[oatm].n12; i++) {
    xatm = tp->nb1234[oatm].L12[i];
    glist[xatm] = labl;
    RecursiveBondSeek(xatm, labl, glist, tp);
  }
}

//-----------------------------------------------------------------------------
// DefineAtomChains: loop over all atoms, recursively listing all atoms 
//                   which are linked by bonds.                         
//
// Arguments:
//   tp      : topology to operate on
//-----------------------------------------------------------------------------
static void DefineAtomChains(prmtop *tp)
{
  int i, j, k, ngrp;
  int* glist;

  // Set a list of atoms that have been linked
  glist = (int*)malloc(tp->natom*sizeof(int));
  SetIVec(glist, tp->natom, -1);
  ngrp = 0;
  for (i = 0; i < tp->natom; i++) {

    // If this atom is already in a group, continue
    if (glist[i] >= 0) {
      continue;
    }
    
    // Recursively find all atoms linked to this one via bonds
    RecursiveBondSeek(i, ngrp, glist, tp);
    ngrp++;
  }

  // Write the topology key for linked groups of atoms
  tp->ngrp = ngrp;
  tp->lgrps = (lgrp*)malloc(tp->ngrp*sizeof(lgrp));
  for (i = 0; i < tp->ngrp; i++) {
    tp->lgrps[i].natom = 0;
    for (j = 0; j < tp->natom; j++) {
      if (glist[j] == i) {
        tp->lgrps[i].natom += 1;
      }
    }
    tp->lgrps[i].atoms = (int*)malloc(tp->lgrps[i].natom*sizeof(int));
    k = 0;
    for (j = 0; j < tp->natom; j++) {
      if (glist[j] == i) {
        tp->lgrps[i].atoms[k] = j;
        k++;
      }
    }
  }

  // Free allocated memory
  free(glist);
}

//-----------------------------------------------------------------------------
// Ascii2Mem: scan an ascii text file verbatim into memory; this routine is
//            particularly useful for scanning topology files for rapid
//            parsing.  The output of this function is a cmat character matrix.
//
// Arguments:                                                           
//   fname:     the flag to search for                                  
//   maxw:      the maximum length of a line in the file                
//   stride:    the number of blocks to allocate pre-emptively, each time the
//              last block is filled
//   errmsg:    an error message to help interpret failures in this function
//-----------------------------------------------------------------------------
cmat Ascii2Mem(char* fname, int maxw, int stride, char* errmsg)
{
  int i, j, k, llen, lnum, nblock, maxblock, nline;
  char *ctmp, *cjtmp;
  char* line;
  FILE *inp;
  cmat C;
  cmat* Cblock;

  // Open the file
  if ((inp = fopen(fname, "r")) == NULL) {
    printf("Ascii2Mem >> Error.  File %s not found!\nAscii2Mem >> %s\n", fname,
           errmsg);
    exit(1);
  }

  // Allocate the matrix array Ct with 2 blocks initially.
  // A block is 128 lines of length maxw.                 
  line = (char*)malloc(maxw*sizeof(char));
  maxblock = stride;
  Cblock = (cmat*)malloc(maxblock*sizeof(cmat));
  for (i = 0; i < stride; i++) {
    Cblock[i] = CreateCmat(128, maxw);
  }
  nblock = 0;
  lnum = 0;
  while (fgets(line, maxw, inp) != NULL) {
    llen = strlen(line);
    for (i = 0; i < llen; i++) {
      Cblock[nblock].map[lnum][i] = line[i];
    }
    lnum++;
    if (lnum == 128) {
      nblock++;
      lnum = 0;
      if (nblock == maxblock) {
        maxblock += stride;
        Cblock = (cmat*)realloc(Cblock, maxblock*sizeof(cmat));
        for (i = nblock; i < maxblock; i++) {
          Cblock[i] = CreateCmat(128, maxw);
        }
      }
    }
  }

  // At this point, the entire file should be in memory.
  // Now, translate it into a single cmat struct.       
  nline = nblock*128 + lnum;
  C = CreateCmat(nline, maxw);
  for (i = 0; i < nblock; i++) {
    for (j = 0; j < 128; j++) {
      ctmp = C.map[i*128+j];
      cjtmp = Cblock[i].map[j];
      for (k = 0; k < maxw; k++) {
        ctmp[k] = cjtmp[k];
      }
    }
  }
  for (j = 0; j < lnum; j++) {
    ctmp = C.map[nblock*128+j];
    cjtmp = Cblock[nblock].map[j];
    for (k = 0; k < maxw; k++) {
      ctmp[k] = cjtmp[k];
    }
  }

  // Free allocated memory
  free(line);
  for (i = 0; i < maxblock; i++) {
    DestroyCmat(&Cblock[i]);
  }
  free(Cblock);
  fclose(inp);

  return C;
}

//-----------------------------------------------------------------------------
// ScanToFlag: function that scans a cmat struct representing an ascii file
//             read into memory, seeking a keyword.  After finding the keyword,
//             ScanToFlag verifies that the next line is a comment and then
//             returns the number of the second line after the one containing
//             the keyword.                    
//                                                                      
// Arguments:                                                           
//   Ctop:       the topology, in ascii format                          
//   keyword:    the keyword to search for                              
//   Lstart:     the suggested starting line                            
//   reqinfo:    flag to indicate that the keyword is vital; the program will
//               abort if the keyword is not found                 
//-----------------------------------------------------------------------------
static int ScanToFlag(cmat *Ctop, char* keyword, int Lstart, int reqinfo)
{
  int i, klen;

  // Get the length of the keyword
  klen = strlen(keyword);

  // Scan from the starting line to the end
  i = (Lstart >= 0) ? Lstart : 0;
  while (i < Ctop->row) {
    if (Ctop->map[i][0] == '%' &&
        strncmp(&Ctop->map[i][6], keyword, klen) == 0) {
      if (Ctop->map[i+1][0] != '%') {
        printf("ScanToFlag >> Error.  No comment line after keyword %s.\n",
               keyword);
        exit(1);
      }
      return i+2;
    }
    i++;
  }

  // Scan from the beginning
  for (i = 0; i < Lstart; i++) {
    if (Ctop->map[i][0] == '%' &&
        strncmp(&Ctop->map[i][6], keyword, klen) == 0) {
      if (Ctop->map[i+1][0] != '%') {
        printf("ScanToFlag >> Error.  No comment line after keyword %s.\n",
               keyword);
        exit(1);
      }
      return i+2;
    }
  }

  // Return -1 or exit if the information was vital
  if (reqinfo == 1) {
    printf("ScanToFlag >> Error.  Could not locate keyword %s.\n", keyword);
    exit(1);
  }
  return -1;
}

//-----------------------------------------------------------------------------
// ReadXiL: read a set of numbers from a memory-commited topology file 
//          assuming the format (X=6,10)i8.  This function returns an array of
//          integers.                                               
//                                                                      
// Arguments:                                                           
//   Ctop:       the topology, in ascii format                          
//   lnum:       going into the function, this is the line on which the 
//               integers start; coming out of the function this is the last 
//               line containing integers, or possibly the next line after
//               that if there were 8 integers on the final line from which
//               integers were read--either way, it is a good point to start
//               looking for the next keyword     
//   N:          the number of integers to read in total
//   X:          the number of integers to read per line
//   L:          the length of character string alotted to each integer
//-----------------------------------------------------------------------------
static int* ReadXiL(cmat *Ctop, int *lnum, int N, int X, int L)
{
  int h, i, j, lcurr, endSeen;
  int* iread;
  char nbuf[32];

  // Allocate the array to be returned
  iread = (int*)malloc(N * sizeof(int));

  // Prepare the buffer array
  for (i = 0; i < 32; i++) {
    nbuf[i] = '\0';
  }
  h = 0;
  lcurr = *lnum;
  for (i = 0; i < N; i++) {
    endSeen = 0;
    for (j = 0; j < L; j++) {
      if (Ctop->map[lcurr][h*L+j] == '\0' || Ctop->map[lcurr][h*L+j] == '\n') {
	endSeen = 1;
	break;
      }
      else {
        nbuf[j] = Ctop->map[lcurr][h*L+j];
      }
    }
    iread[i] = (endSeen == 0) ? atoi(nbuf) : 0;
    h++;
    if (h == X) {
      h = 0;
      lcurr += 1;
      if (lcurr >= Ctop->row && i < N-1) {
        printf("ReadXiL >> Error.  Seeking %d integers, encountered end of "
               "file.\n", N);
        exit(1);
      }
    }
  }
  *lnum = lcurr;

  return iread;
}

//-----------------------------------------------------------------------------
// Read10i8Test: read a set of integers from a memory-commited topology file
//               assuming the format 10i8.  The length of the list is not
//               known when this function is called; rather, it will read
//               so long as there are integers in this format, until a new
//               section or the end of the file is encountered.
//
// Arguments:
//   Ctop:       the topology, in ascii format
//   lnum:       going into the function, this is the line on which the
//               integers start; coming out of the function this is the last
//               line containing integers, or possibly the next line after
//               that if there were 8 integers on the final line from which
//               integers were read--either way, it is a good point to start
//               looking for the next keyword
//   N:          the number of integers that have been read (returned)
//-----------------------------------------------------------------------------
static int* Read10i8Test(cmat *Ctop, int *lnum, int *N)
{
  int i, j, nints, maxints, lnterm, currln, currpos, valid;
  int* iread;
  char nbuf[16];

  // Allocate the array to be returned
  maxints = 32;
  iread = (int*)malloc(maxints * sizeof(int));

  // Find the last line
  lnterm = *lnum;
  while (lnterm < Ctop->row && strncmp(Ctop->map[lnterm], "%FLAG", 5) != 0) {
    lnterm++;
  }

  // Start reading numbers
  currln = *lnum;
  currpos = 0;
  valid = 1;
  nints = 0;
  while (valid == 1) {
    for (i = 0; i < 8; i++) {
      nbuf[i] = Ctop->map[currln][currpos+i];
    }
    nbuf[8] = '\0';

    // Test whether the word just read is an integer 
    if (WordIsInteger(nbuf) == 0) {
      valid = 0;
    }
    else {
      iread[nints] = atoi(nbuf);
      nints++;
      if (nints == maxints) {
	maxints *= 2;
	iread = (int*)realloc(iread, maxints * sizeof(int));
      }
    }
    currpos += 8;
    if (currpos == 80) {
      currpos = 0;
      currln++;
      if (currln == lnterm) {
	valid = 0;
      }
    }
  }
  if (currln < lnterm - 1) {
    printf("Read10i8Test >> Warning.  Unexpected end of integer content at "
	   "line %d.\nRead10i8Test >> Expected integers to continue until "
	   "line %d.\n", currln, lnterm);
  }

  // Return results
  *N = nints;
  *lnum = currln;
  iread = (int*)realloc(iread, nints * sizeof(int));
  return iread;
}

//-----------------------------------------------------------------------------
// Read5e16: read a set of numbers from a memory-commited topology file 
//           assuming the format 5e16.  This function returns an array of
//           doubles.                                                
//                                                                      
// Arguments:                                                           
//   Ctop:       the topology, in ascii format                          
//   lnum:       going into the function, this is the line on which the 
//               integers start; coming out of the function this is the last 
//               line containing integers, or possibly the next line after
//               that if there were 8 integers on the final line from which
//               integers were read--either way, it is a good point to start
//               looking for the next keyword     
//   N:          the number of doubles to read                          
//-----------------------------------------------------------------------------
static double* Read5e16(cmat *Ctop, int *lnum, int N)
{
  int h, i, j, lcurr;
  double* dread;
  char nbuf[17];

  // Allocate the array to be returned
  dread = (double*)malloc(N*sizeof(double));

  // Prepare the buffer array
  nbuf[16] = '\0';
  h = 0;
  lcurr = *lnum;
  for (i = 0; i < N; i++) {
    for (j = 0; j < 16; j++) {
      nbuf[j] = Ctop->map[lcurr][h*16+j];
    }
    dread[i] = atof(nbuf);
    h++;
    if (h == 5) {
      h = 0;
      lcurr += 1;
      if (lcurr >= Ctop->row && i < N-1) {
        printf("Read5e16 >> Error.  Seeking %d reals, encountered end of "
               "file.\n", N);
        exit(1);
      }
    }
  }
  *lnum = lcurr;

  return dread;
}

//-----------------------------------------------------------------------------
// Read20a4: read a set of characters in the format 20a4 and place them in one
//           long continuous array.  Starting from a character matrix this is
//           actually quite simple.                      
//                                                                      
// Arguments:                                                           
//   Ctop:       the topology, in ascii format                          
//   lnum:       going into the function, this is the line on which the 
//               character sets start; coming out of the function this is the
//               last line containing character sets, or possibly the next line
//               after that if there were 8 integers on the final line from
//               which integers were read--either way, it is a good point to
//               start looking for the next keyword     
//   N:          the number of character sets to read                   
//-----------------------------------------------------------------------------
static char* Read20a4(cmat *Ctop, int *lnum, int N)
{
  int h, i, j, lcurr;
  char* cread;

  // Allocate array for return
  cread = (char*)malloc(4*N*sizeof(char));
  h = 0;
  lcurr = *lnum;
  for (i = 0; i < N; i++) {
    for (j = 0; j < 4; j++) {
      cread[4*i+j] = Ctop->map[lcurr][h*4+j];
    }
    h++;
    if (h == 20) {
      h = 0;
      lcurr += 1;
      if (lcurr >= Ctop->row && i < N-1) {
        printf("Read20a4 >> Error.  Seeking %d character sets, encountered "
               "end of file.\n", N);
        exit(1);
      }
    }
  }
  *lnum = lcurr;

  return cread;
}

//-----------------------------------------------------------------------------
// PrmTopPreamble: this dissects the preamble to a topology by first reading
//                 an array of 32 integers after the keyword "POINTERS" and
//                 then picking it apart.                
//                                                                      
// Arguments:                                                           
//   Ctop:       the topology, in ascii format                          
//   tp:         the topology struct                                    
//-----------------------------------------------------------------------------
static void PrmTopPreamble(cmat *Ctop, prmtop *tp, int *npline)
{
  int* pointers;

  // Check to see that this is a LEaP prmtop
  // format file, not a chamber prmtop.     
  *npline = ScanToFlag(Ctop, "CTITLE", 0, 0);
  if (*npline >= 0) {
    printf("PrmTopPreamble >> Error.  Unable to read chamber prmtop format"
           ".\n");
    exit(1);
  }

  *npline = ScanToFlag(Ctop, "POINTERS", 0, 1);
  pointers = ReadXiL(Ctop, npline, 32, 10, 8);
  tp->natom = pointers[0];
  tp->ntypes = pointers[1];
  tp->withH.nbond = pointers[2];
  tp->woH.nbond = pointers[3];
  tp->withH.nangl = pointers[4];
  tp->woH.nangl = pointers[5];
  tp->withH.ndihe = pointers[6];
  tp->woH.ndihe = pointers[7];
  tp->nhparm = pointers[8];
  tp->nparm = pointers[9];
  tp->tnexcl = pointers[10];
  tp->nres = pointers[11];
  tp->woHC.nbond = pointers[12]; 
  tp->woHC.nangl = pointers[13];
  tp->woHC.ndihe = pointers[14];
  tp->nBAH.nbond = pointers[15];
  tp->nBAH.nangl = pointers[16];
  tp->nBAH.ndihe = pointers[17];
  tp->natyp = pointers[18];
  tp->nphb = pointers[19];
  tp->ifpert = pointers[20];
  tp->pert.nbond = pointers[21];
  tp->pert.nangl = pointers[22];
  tp->pert.ndihe = pointers[23];
  tp->wpert.nbond = pointers[24];
  tp->wpert.nangl = pointers[25];
  tp->wpert.ndihe = pointers[26];
  tp->ifbox = pointers[27];
  tp->nmxrs = pointers[28];
  tp->ifcap = pointers[29];
  tp->numextra = pointers[30];
  tp->ncopy = pointers[31];

  // Set some defaults for this prmtop struct here
  tp->EPInserted = 0;

  // Free allcoated memory
  free(pointers);
}

//-----------------------------------------------------------------------------
// EnumerateBonds: this function loops through all bond lists, creating 
//                 atom-specific lists of bonds based on the Ath atom in the
//                 A-B sequence.  The Ath atom "owns" the bond; cells which
//                 control the Ath atom (cells which have the Ath atom in their
//                 primary sector) are responsible for figuring out where other
//                 atoms in the bond are and then computing forces.
//                                                                      
//                 This function should be called once, after the system
//                 topology has been read, as is done at the end of GetPrmTop()
//                 below.                                   
//                                                                      
// Arguments:                                                           
//   tp:  the topology                                                  
//-----------------------------------------------------------------------------
static void EnumerateBonds(prmtop *tp)
{
  int i, j, k;

  // Allocate the array for all atoms
  tp->BLC = (bondlist*)calloc(tp->natom, sizeof(bondlist));

  // One pass to figure out how many bonds each atom controls.
  for (i = 0; i < tp->withH.nbond; i++) {
    j = tp->BIncH[i].a;
    if (j < 0) {
      printf("EnumerateBonds >> Error.  Bond anchor index value is %d.\n"
             "EnumerateBonds >> No anchor index value should be less than "
             "zero.\n", j);
      exit(1);
    }
    tp->BLC[j].nbond += 1;
  }
  for (i = 0; i < tp->woH.nbond; i++) {
    j = tp->BNoH[i].a;
    if (j < 0) {
      printf("EnumerateBonds >> Error.  Bond anchor index value is %d.\n"
             "EnumerateBonds >> No anchor index value should be less than "
             "zero.\n", j);
      exit(1);
    }
    tp->BLC[j].nbond += 1;
  }

  // Allocate the bond arrays for each individual atom
  for (i = 0; i < tp->natom; i++) {
    tp->BLC[i].BC = (bondcomm*)malloc(tp->BLC[i].nbond*sizeof(bondcomm));
    tp->BLC[i].nbond = 0;
  }

  // A second pass to index the bonds by anchor atoms
  for (i = 0; i < tp->withH.nbond; i++) {
    j = tp->BIncH[i].a;
    k = tp->BLC[j].nbond;
    tp->BLC[j].BC[k].a = j;
    tp->BLC[j].BC[k].b = tp->BIncH[i].b;
    tp->BLC[j].BC[k].t = tp->BIncH[i].idx;
    tp->BLC[j].BC[k].H = 1;
    tp->BLC[j].nbond = k+1;
  }
  for (i = 0; i < tp->woH.nbond; i++) {
    j = tp->BNoH[i].a;
    k = tp->BLC[j].nbond;
    tp->BLC[j].BC[k].a = j;
    tp->BLC[j].BC[k].b = tp->BNoH[i].b;
    tp->BLC[j].BC[k].t = tp->BNoH[i].idx;
    tp->BLC[j].BC[k].H = 0;
    tp->BLC[j].nbond = k+1;
  }
}

//-----------------------------------------------------------------------------
// EnumerateAngles: this function loops through all angle lists, making 
//                  atom-specific lists of angles based on the Bth atom in the
//                  A-B-C sequence.  The Bth atom "owns" the angle.
//                                                                      
//                  This function should be called once, after the system
//                  topology has been read, as is done at the end of
//                  GetPrmTop() below.                               
//                                                                      
// Arguments:                                                           
//   tp:  the topology                                                  
//-----------------------------------------------------------------------------
static void EnumerateAngles(prmtop *tp)
{
  int i, j, k, atma, atmc;

  // Allocate the array for all atoms
  tp->ALC = (angllist*)calloc(tp->natom, sizeof(angllist));

  // One pass to figure out how many angles each atom controls.
  for (i = 0; i < tp->withH.nangl; i++) {
    j = tp->AIncH[i].b;
    if (j < 0) {
      printf("EnumerateAngles >> Error.  Angle anchor index value is %d.\n"
             "EnumerateAngles >> No anchor index value should be less than "
             "zero.\n", j);
      exit(1);
    }
    tp->ALC[j].nangl += 1;
  }
  for (i = 0; i < tp->woH.nangl; i++) {
    j = tp->ANoH[i].b;
    if (j < 0) {
      printf("EnumerateAngles >> Error.  Angle anchor index value is %d.\n"
             "EnumerateAngles >> No anchor index value should be less than "
             "zero.\n", j);
      exit(1);
    }
    tp->ALC[j].nangl += 1;
  }

  // Allocate the angle arrays for each individual atom
  for (i = 0; i < tp->natom; i++) {
    tp->ALC[i].AC = (anglcomm*)malloc(tp->ALC[i].nangl*sizeof(anglcomm));
    tp->ALC[i].nangl = 0;
  }

  // A second pass to index the angles by anchor atoms
  for (i = 0; i < tp->withH.nangl; i++) {
    j = tp->AIncH[i].b;
    k = tp->ALC[j].nangl;
    tp->ALC[j].AC[k].a = tp->AIncH[i].a;
    tp->ALC[j].AC[k].c = tp->AIncH[i].c;
    tp->ALC[j].AC[k].t = tp->AIncH[i].idx;
    tp->ALC[j].nangl = k+1;
  }
  for (i = 0; i < tp->woH.nangl; i++) {
    j = tp->ANoH[i].b;
    k = tp->ALC[j].nangl;
    tp->ALC[j].AC[k].a = tp->ANoH[i].a;
    tp->ALC[j].AC[k].c = tp->ANoH[i].c;
    tp->ALC[j].AC[k].t = tp->ANoH[i].idx;
    tp->ALC[j].nangl = k+1;
  }

  // A final pass to determine whether to exclude 1:3 interactions
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < tp->ALC[i].nangl; j++) {

      // By default, 1:3 interactions are excluded
      tp->ALC[i].AC[j].excl = 1;
      atma = tp->ALC[i].AC[j].a;
      atmc = tp->ALC[i].AC[j].c;
      for (k = 0; k < tp->BLC[atma].nbond; k++) {
        if (tp->BLC[atma].BC[k].b == atmc) {
          tp->ALC[i].AC[j].excl = 0;
          break;
        }
      }
      if (tp->ALC[i].AC[j].excl == 1) {
        for (k = 0; k < tp->BLC[atmc].nbond; k++) {
          if (tp->BLC[atmc].BC[k].b == atma) {
            tp->ALC[i].AC[j].excl = 0;
            break;
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// EnumerateDihedrals: this function loops through all dihedral lists, creating
//                     atom-specific lists of dihedrals based on the Bth atom
//                     in the A-B-C-D sequence.  The Bth atom "owns" the
//                     dihedral; cells which control the Bth atom for bonded
//                     interactions are responsible for figuring out where the
//                     other atoms in the dihedral are and then computing
//                     forces.          
//                                                                      
//                     This function should be called once, after the system
//                     topology has been read, as is done at the end of
//                     GetPrmTop() below.                        
//                                                                      
// Arguments:                                                           
//   tp:  the topology                                                  
//-----------------------------------------------------------------------------
static void EnumerateDihedrals(prmtop *tp)
{
  int i, j, k, maxdihe, dihefound;
  dihelist scrd;

  // Allocate the array for all atoms
  tp->HLC = (dihelist*)calloc(tp->natom, sizeof(dihelist));

  // One pass to figure out how many dihedrals each atom controls.
  // The distinction between dihedrals with and without hydrogen  
  // is a historical artifact, so we lump them all together here. 
  for (i = 0; i < tp->withH.ndihe; i++) {
    j = tp->HIncH[i].b;
    if (j < 0) {
      printf("EnumerateDihedrals >> Error.  Dihedral anchor index value is "
             "%d.\nEnumerateDihedrals >> No anchor index value should be less "
             "than zero.\n", j);
      exit(1);
    }
    tp->HLC[j].ndihe += 1;
  }
  for (i = 0; i < tp->woH.ndihe; i++) {
    j = tp->HNoH[i].b;
    if (j < 0) {
      printf("EnumerateDihedrals >> Error.  Dihedral anchor index value is "
             "%d.\nEnumerateDihedrals >> No anchor index value should be less "
             "than zero.\n", j);
      exit(1);
    }
    tp->HLC[j].ndihe += 1;
  }

  // Allocate the dihedral arrays for each individual atom
  for (i = 0; i < tp->natom; i++) {
    tp->HLC[i].HC = (dihecomm*)malloc(tp->HLC[i].ndihe*sizeof(dihecomm));
    tp->HLC[i].ndihe = 0;
  }

  // Another pass to assign dihedrals to individual atoms
  for (i = 0; i < tp->withH.ndihe; i++) {
    j = tp->HIncH[i].b;
    k = tp->HLC[j].ndihe;
    tp->HLC[j].HC[k].impr = (tp->HIncH[i].d < 0) ? 1 : 0;
    tp->HLC[j].HC[k].a = abs(tp->HIncH[i].a);
    tp->HLC[j].HC[k].eval14 = (tp->HIncH[i].c < 0) ? 0 : 1;
    tp->HLC[j].HC[k].scee = tp->scee[tp->HIncH[i].idx];
    tp->HLC[j].HC[k].scnb = tp->scnb[tp->HIncH[i].idx];
    tp->HLC[j].HC[k].c = abs(tp->HIncH[i].c);
    tp->HLC[j].HC[k].d = abs(tp->HIncH[i].d);
    tp->HLC[j].HC[k].nt = tp->HIncH[i].idx;
    tp->HLC[j].ndihe = k+1;
  }
  for (i = 0; i < tp->woH.ndihe; i++) {
    j = tp->HNoH[i].b;
    k = tp->HLC[j].ndihe;
    tp->HLC[j].HC[k].impr = (tp->HNoH[i].d < 0) ? 1 : 0;
    tp->HLC[j].HC[k].a = abs(tp->HNoH[i].a);
    tp->HLC[j].HC[k].eval14 = (tp->HNoH[i].c < 0) ? 0 : 1;
    tp->HLC[j].HC[k].scee = tp->scee[tp->HNoH[i].idx];
    tp->HLC[j].HC[k].scnb = tp->scnb[tp->HNoH[i].idx];
    tp->HLC[j].HC[k].c = abs(tp->HNoH[i].c);
    tp->HLC[j].HC[k].d = abs(tp->HNoH[i].d);
    tp->HLC[j].HC[k].nt = tp->HNoH[i].idx;
    tp->HLC[j].ndihe = k+1;
  }

  // Now a pass to condense dihedrals that share the same atoms
  maxdihe = 0;
  for (i = 0; i < tp->natom; i++) {
    if (tp->HLC[i].ndihe > maxdihe) {
      maxdihe = tp->HLC[i].ndihe;
    }
  }
  scrd.HC = (dihecomm*)malloc(maxdihe*sizeof(dihecomm));
  for (i = 0; i < maxdihe; i++) {
    scrd.HC[i].t = (int*)malloc(maxdihe*sizeof(int));
  }
  for (i = 0; i < tp->natom; i++) {

    // Initialize the buffer dihedral command structure
    scrd.ndihe = 0;
    for (j = 0; j < maxdihe; j++) {
      scrd.HC[j].nt = 0;
    }

    // Loop over all dihedrals controlled by this atom and
    // contribute them to the buffer structure.           
    for (j = 0; j < tp->HLC[i].ndihe; j++) {
      dihefound = 0;
      for (k = 0; k < scrd.ndihe; k++) {
        if (scrd.HC[k].a == tp->HLC[i].HC[j].a &&
            scrd.HC[k].c == tp->HLC[i].HC[j].c &&
            scrd.HC[k].d == tp->HLC[i].HC[j].d &&
            scrd.HC[k].impr == tp->HLC[i].HC[j].impr) {
          dihefound = 1;
          scrd.HC[k].t[scrd.HC[k].nt] = tp->HLC[i].HC[j].nt;
          scrd.HC[k].nt += 1;
          if (tp->HLC[i].HC[j].eval14 == 1) {
            scrd.HC[k].eval14 = 1;
          }
          break;
        }
      }

      // Make a new dihedral for this atom
      if (dihefound == 0) {
        scrd.HC[scrd.ndihe].a = tp->HLC[i].HC[j].a;
        scrd.HC[scrd.ndihe].c = tp->HLC[i].HC[j].c;
        scrd.HC[scrd.ndihe].d = tp->HLC[i].HC[j].d;
        scrd.HC[scrd.ndihe].impr = tp->HLC[i].HC[j].impr;
        scrd.HC[scrd.ndihe].eval14 = tp->HLC[i].HC[j].eval14;
        scrd.HC[scrd.ndihe].scee = tp->HLC[i].HC[j].scee;
        scrd.HC[scrd.ndihe].scnb = tp->HLC[i].HC[j].scnb;
        scrd.HC[scrd.ndihe].t[0] = tp->HLC[i].HC[j].nt;
        scrd.HC[scrd.ndihe].nt = 1;
        scrd.ndihe += 1;
      }
    }

    // Copy the buffer dihedral command structure into the topology
    tp->HLC[i].ndihe = scrd.ndihe;
    free(tp->HLC[i].HC);
    tp->HLC[i].HC = (dihecomm*)malloc(scrd.ndihe*sizeof(dihecomm));
    for (j = 0; j < scrd.ndihe; j++) {
      tp->HLC[i].HC[j].a = scrd.HC[j].a;
      tp->HLC[i].HC[j].c = scrd.HC[j].c;
      tp->HLC[i].HC[j].d = scrd.HC[j].d;
      tp->HLC[i].HC[j].impr = scrd.HC[j].impr;
      tp->HLC[i].HC[j].eval14 = scrd.HC[j].eval14;
      tp->HLC[i].HC[j].scee = scrd.HC[j].scee;
      tp->HLC[i].HC[j].scnb = scrd.HC[j].scnb;
      tp->HLC[i].HC[j].nt = scrd.HC[j].nt;
      tp->HLC[i].HC[j].t = CpyIVec(scrd.HC[j].t, scrd.HC[j].nt);
    }
  }

  // Free allocated memory
  for (i = 0; i < maxdihe; i++) {
    free(scrd.HC[i].t);
  }
  free(scrd.HC);
}

//-----------------------------------------------------------------------------
// EnumerateCmaps: this function loops through the list of CMAP terms to assign
//                 each CMAP to the C atom in its sequence.  The Amber CMAP
//                 terms are really a special case of the more general CHARMM
//                 CMAPs that can involve any two dihedrals without any atoms
//                 shared between them.  In Amber CMAPs, the C atom is the
//                 third atom of the first (phi) dihedral and the second atom
//                 of the second (psi) dihedral.
//
// Arguments:
//   tp:  the topology
//-----------------------------------------------------------------------------
static void EnumerateCmaps(prmtop *tp)
{
  int i, j, atmidx, cidx;
  int* cmcount;
  
  // Allocate the array for all atoms
  tp->MLC = (cmaplist*)calloc(tp->natom, sizeof(cmaplist));
  
  // One pass to find out how many terms each atom controls
  cmcount = (int*)calloc(tp->natom, sizeof(int));
  for (i = 0; i < tp->ncmap; i++) {
    cmcount[tp->Mterms[i].c] += 1;
  }

  // Allocate the CMAP commands for all atoms, then
  // assign pointers to subsections of the array.
  if (tp->ncmap == 0) {
    tp->MLC[0].MC = (cmapcomm*)malloc(sizeof(cmapcomm));
  }
  else {
    tp->MLC[0].MC = (cmapcomm*)malloc(tp->ncmap * sizeof(cmapcomm));
  }
  j = 0;
  for (i = 0; i < tp->natom; i++) {
    tp->MLC[i].MC = &tp->MLC[0].MC[j];
    tp->MLC[i].ncmap = cmcount[i];
    j += cmcount[i];
    cmcount[i] = 0;
  }

  // Assign CMAPs to atoms
  for (i = 0; i < tp->ncmap; i++) {
    atmidx =  tp->Mterms[i].c;
    cidx = cmcount[atmidx];
    tp->MLC[atmidx].MC[cidx].a   = tp->Mterms[i].a;
    tp->MLC[atmidx].MC[cidx].b   = tp->Mterms[i].b;
    tp->MLC[atmidx].MC[cidx].d   = tp->Mterms[i].d;
    tp->MLC[atmidx].MC[cidx].e   = tp->Mterms[i].e;
    tp->MLC[atmidx].MC[cidx].idx = tp->Mterms[i].idx;
    tp->MLC[i].MC[cidx].eval14i  = 0;
    tp->MLC[i].MC[cidx].eval14ii = 0;
    tp->MLC[i].MC[cidx].scnbi    = 0.0;
    tp->MLC[i].MC[cidx].scnbii   = 0.0;
    tp->MLC[i].MC[cidx].sceei    = 0.0;
    tp->MLC[i].MC[cidx].sceeii   = 0.0;
    cmcount[tp->Mterms[i].c] = cidx + 1;
  }

  // Free allocated memory
  free(cmcount);
}

//-----------------------------------------------------------------------------
// OrderLJParameters: this routine re-organizes van-der Waals parameters into
//                    things that are more intuitive to read.      
//                                                                      
// Arguments:                                                           
//   tp:  the topology                                                  
//-----------------------------------------------------------------------------
void OrderLJParameters(prmtop *tp)
{
  int i, j, haslj;
  double aval, bval, sig, eps, gam, cval, dval;
  double *ftmp, *utmp;

  // Build the Lennard-Jones tables
  if (tp->ljbuck == 0) {
    tp->LJftab = CreateDmat(tp->ntypes, 2*tp->ntypes, 0);
    tp->LJutab = CreateDmat(tp->ntypes, 2*tp->ntypes, 0);
    for (i = 0; i < tp->ntypes; i++) {
      ftmp = tp->LJftab.map[i];
      utmp = tp->LJutab.map[i];
      for (j = 0; j < tp->ntypes; j++) {

        // There is no support for the 10-12 hydrogen bonding potential.
        // Atoms whose van-der Waals types are less than 0 are assumed  
        // to have no van-der Waals properties whatsoever.              
        if (tp->NBParmIdx[i*tp->ntypes+j] >= 0) {
          aval = tp->LJA[tp->NBParmIdx[i*tp->ntypes+j]];
          bval = tp->LJB[tp->NBParmIdx[i*tp->ntypes+j]];
          if (aval > 1.0e-8) {
            sig = pow(aval/bval, 1.0/6.0);
            eps = bval/pow(sig, 6.0);
            ftmp[2*j] = -12.0*eps*pow(sig, 12.0);
            ftmp[2*j+1] = 6.0*eps*pow(sig, 6.0);
            utmp[2*j] = eps*pow(sig, 12.0);
            utmp[2*j+1] = -eps*pow(sig, 6.0);
          }
          else {
            ftmp[2*j] = 0.0;
            ftmp[2*j+1] = 0.0;
            utmp[2*j] = 0.0;
            utmp[2*j+1] = 0.0;
          }
        }
        else {
          ftmp[2*j] = 0.0;
          ftmp[2*j+1] = 0.0;
          utmp[2*j] = 0.0;
          utmp[2*j+1] = 0.0;
        }
      }
    }
  }
  else if (tp->ljbuck == 1) {
    tp->LJftab = CreateDmat(tp->ntypes, 4*tp->ntypes, 0);
    tp->LJutab = CreateDmat(tp->ntypes, 4*tp->ntypes, 0);
    for (i = 0; i < tp->ntypes; i++) {
      ftmp = tp->LJftab.map[i];
      utmp = tp->LJutab.map[i];
      for (j = 0; j < tp->ntypes; j++) {

        // This modified Buckingham exp-6 potential uses a (1/r)^12 term
        // to continue the positive divergence of the potential as r    
        // goes to zero, such that there isn't even an inflection point.
        // This (1/r)^12 terms requires its own scaling factor, which is
        // computed here based on the exp terms (aval and bval) and the 
        // (1/r)^6 term (cval).  The potential is, succinctly:          
        //                                                              
        //      U(r) = aval*exp(-bval*r) - cval/(r^6) + dval/(r^12)     
        //                                                              
        // The terms aval, bval, and cval are computed based on the eps,
        // gam, and sig parameters which are actually read in from the  
        // topology file.  This is in contrast to the Lennard-Jones     
        // potential, where it is aval and bval which are read from the 
        // topology file and eps and sig which are back-calculated.     
        if (tp->NBParmIdx[i*tp->ntypes+j] >= 0) {
          eps = tp->LJA[tp->NBParmIdx[i*tp->ntypes+j]];
          gam = tp->LJB[tp->NBParmIdx[i*tp->ntypes+j]];
          sig = tp->LJC[tp->NBParmIdx[i*tp->ntypes+j]];
          aval = 6.0*eps*exp(eps)/(gam-6.0);
          bval = gam/sig;
          cval = eps*pow(sig, 6.0)/(1.0-6.0/gam);
          dval = DetermineDCoef(eps, gam, sig);
          utmp[4*j] = aval;
          utmp[4*j+1] = bval;
          utmp[4*j+2] = cval;
          utmp[4*j+3] = dval;
          ftmp[4*j] = -aval*bval;
          ftmp[4*j+1] = bval;
          ftmp[4*j+2] = -6.0*cval;
          ftmp[4*j+3] = -12.0*dval;
        }
        else {
          ftmp[4*j] = 0.0;
          ftmp[4*j+1] = 0.0;
          ftmp[4*j+2] = 0.0;
          ftmp[4*j+3] = 0.0;
          utmp[4*j] = 0.0;
          utmp[4*j+1] = 0.0;
          utmp[4*j+2] = 0.0;
          utmp[4*j+3] = 0.0;
        }
      }
    }
  }
  else {
    printf("OrderParameters >> Error.  Unrecognized van-der Waals model type "
           "%d.\n", tp->ljbuck);
    exit(1);
  }

  // Identify atoms that have no Lennard-Jones properties
  for (i = 0; i < tp->ntypes; i++) {
    haslj = 0;
    utmp = tp->LJutab.map[i];
    for (j = 0; j < tp->ntypes; j++) {
      if (fabs(utmp[2*j]) > 1.0e-8 || fabs(utmp[2*j+1]) > 1.0e-8) {
        haslj = 1;
        break;
      }
    }
    if (haslj == 0) {
      for (j = 0; j < tp->natom; j++) {
        if (tp->LJIdx[j] == i) {
          tp->LJIdx[j] = -i;
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// OrderBondParameters: this routine re-organizes bond, angle, and dihedral
//                      parameters into things that are more intuitive to read.
//                                                                      
// Arguments:                                                           
//   tp:  the topology                                                  
//-----------------------------------------------------------------------------
static void OrderBondParameters(prmtop *tp)
{
  int i;

  // Organize bond parameters
  tp->BParam = (bonddef*)malloc(tp->nBAH.nbond*sizeof(bonddef));
  for (i = 0; i < tp->nBAH.nbond; i++) {
    tp->BParam[i].K       = tp->BondK[i];
    tp->BParam[i].l0      = tp->BondEq[i];
    tp->BParam[i].Kpull   = tp->BondPullK[i];
    tp->BParam[i].lpull0  = tp->BondPullEq[i];
    tp->BParam[i].Kpress  = tp->BondPressK[i];
    tp->BParam[i].lpress0 = tp->BondPressEq[i];
  }

  // Organize angle parameters
  tp->AParam = (angldef*)malloc(tp->nBAH.nangl*sizeof(angldef));
  for (i = 0; i < tp->nBAH.nangl; i++) {
    tp->AParam[i].K = tp->AnglK[i];
    tp->AParam[i].th0 = tp->AnglEq[i];
  }

  // Organize dihedral parameters
  tp->HParam = (dihedef*)malloc(tp->nBAH.ndihe*sizeof(dihedef));
  for (i = 0; i < tp->nBAH.ndihe; i++) {
    tp->HParam[i].K = tp->DiheK[i];
    tp->HParam[i].N = tp->DiheN[i];
    tp->HParam[i].Phi = tp->DihePhi[i];
  }
}

//-----------------------------------------------------------------------------
// FindThreePointWaters: locate all three-point water molecules within the
//                       topology, with the user-specified water name as the
//                       reference.                         
//                                                                      
// Arguments:                                                           
//   tp:  the topology                                                  
//-----------------------------------------------------------------------------
static void FindThreePointWaters(prmtop *tp)
{
  int i, j, nox, nh1, nh2, nwat;

  // Allocate space for the array of fast water indices
  nwat = 0;
  tp->SHL = (cnstcomm*)calloc(tp->natom, sizeof(cnstcomm));

  // If SETTLE is not requested, return
  if (tp->settle == 0) {
    tp->nwat = nwat;
    return;
  }

  // Locate all fast waters
  for (i = 0; i < tp->nres; i++) {
    if (strncmp(&tp->ResNames[4*i], tp->WaterName, 4) == 0) {
      nox = -1;
      nh1 = -1;
      nh2 = -1;
      for (j = tp->ResLims[i]; j < tp->ResLims[i+1]; j++) {
        if (tp->Masses[j] > 15.9 && tp->Masses[j] < 16.1) {            
          nox = j;
        }
        else if (tp->Masses[j] < 1.1 && tp->Masses[j] > 0.9 && nh1 < 0) {
          nh1 = j;
        }
        else if (tp->Masses[j] < 1.1 && tp->Masses[j] > 0.9 && nh2 < 0) {
          nh2 = j;
        }
      }
      if (nox < 0) {
        printf("FindThreePointWaters >> Error.  Oxygen not found for residue "
               "%d.\n", i);
        exit(1);
      }
      if (nh1 < 0 || nh2 < 0) {
        printf("FindThreePointWaters >> Error.  Hydrogen not found for "
               "residue %d.\n", i);
        exit(1);
      }

      // The oxygen atom is flagged with a +1 SHL execution value,
      // to signify that it controls a three-point water style    
      // constraint group.                                        
      tp->SHL[nox].exe = 1;
      tp->SHL[nox].blist = (int*)malloc(2*sizeof(int));
      tp->SHL[nox].blist[0] = nh1;
      tp->SHL[nox].blist[1] = nh2;

      // The hydrogens are flagged with -1 SHL execution values,
      // to signify that they are part of a constraint group    
      // but not in control of it.                              
      tp->SHL[nh1].exe = -1;
      tp->SHL[nh2].exe = -1;
      nwat++;
    }
  }
  tp->nwat = nwat;
}

//-----------------------------------------------------------------------------
// CompSettleGeom: compute the geometric parameters for SETTLE on the   
//                 user-specified water molecules.                      
//                                                                      
// Arguments:                                                           
//   tp:  the topology                                                  
//-----------------------------------------------------------------------------
static void CompSettleGeom(prmtop *tp)
{
  int i, j, k, resfound, nox, atma, atmb, nh1, nh2;
  double dOH, dHH, massO, massH, t1;

  // Find an example of a rigid water molecule
  resfound = 0;
  i = 0;
  while (resfound == 0 && i < tp->nres) {
    for (j = tp->ResLims[i]; j < tp->ResLims[i+1]; j++) {
      if (tp->SHL[j].exe == 1) {
        resfound = 1;
        nox = j;
        nh1 = tp->SHL[j].blist[0];
        nh2 = tp->SHL[j].blist[1];
      }
    }
    if (resfound == 0) {
      i++;
      continue;
    }

    // Get bond length for O-H
    dOH = -1.0;
    dHH = -1.0;
    for (j = tp->ResLims[i]; j < tp->ResLims[i+1]; j++) {
      for (k = 0; k < tp->BLC[j].nbond; k++) {
        atma = tp->BLC[j].BC[k].a;
        atmb = tp->BLC[j].BC[k].b;
        if (tp->Masses[atma] > 14.9 &&
            tp->Masses[atmb] < 3.1 && tp->Masses[atmb] > 0.5) {
          dOH = tp->BParam[tp->BLC[j].BC[k].t].l0;
        }
        else if (tp->Masses[atma] < 3.1 && tp->Masses[atma] > 0.5 &&
                 tp->Masses[atmb] > 14.9) {
          dOH = tp->BParam[tp->BLC[j].BC[k].t].l0;
        }
        else if (tp->Masses[atma] < 3.1 && tp->Masses[atma] > 0.5 &&
                 tp->Masses[atmb] < 3.1 && tp->Masses[atmb] > 0.5) {
          dHH = tp->BParam[tp->BLC[j].BC[k].t].l0;
        }
      }
    }

    // Check that bond lengths are known
    if (dOH < 0.0) {
      printf("CompSettleGeom >> Error.  Unable to determine O-H bond "
             "length.\n");
      exit(1);
    }
    if (dHH < 0.0) {
      printf("CompSettleGeom >> Error.  Unable to determine H-H bond "
             "length.\n");
      exit(1);
    }

    massO = tp->Masses[nox];
    massH = tp->Masses[nh1];
    if (fabs(tp->Masses[nh2] - tp->Masses[nh1]) > 1.0e-9) {
      printf("CompSettleGeom >> Error.  Masses of hydrogen atoms do not "
             "agree.\n");
      exit(1);
    }
    t1 = 0.5 * massO / massH;
    tp->FWtab.rc = 0.5* dHH;
    tp->FWtab.ra = sqrt(dOH*dOH - tp->FWtab.rc*tp->FWtab.rc) / (t1 + 1.0);
    tp->FWtab.rb = t1 * tp->FWtab.ra;
    tp->FWtab.mO = massO;
    tp->FWtab.mH = massH;
  }
  if (resfound == 0) {
    printf("CompSettleGeom >> Error.  Unable to locate a rigid water "
           "molecule.\n");
    exit(1);
  }
}

//-----------------------------------------------------------------------------
// LabelSettleGroups: this function steps through the topology once to  
//                    enumerate all SETTLE-constrained groups of atoms.  Doing
//                    so prevents the code from having to loop back over all
//                    atoms every time a constrained pair must be tested.
//
// Arguments:
//   tp:      the system topology
//-----------------------------------------------------------------------------
static int* LabelSettleGroups(prmtop *tp)
{
  int i, ngrp;
  int* settlegrp;

  settlegrp = (int*)malloc(tp->natom*sizeof(int));
  SetIVec(settlegrp, tp->natom, -1);
  ngrp = 0;
  for (i = 0; i < tp->natom; i++) {
    if (tp->SHL[i].exe == 1) {
      settlegrp[i] = ngrp;
      settlegrp[tp->SHL[i].blist[0]] = ngrp;
      settlegrp[tp->SHL[i].blist[1]] = ngrp;
      ngrp++;
    }
  }

  return settlegrp;
}

//-----------------------------------------------------------------------------
// IsInSettleGroup: this function determines whether a pair of atoms which will
//                  be subjected to a RATTLE constraint is already part of a
//                  SETTLE group and therefore should not be RATTLE'd.
//                                                                      
// Arguments:                                                           
//   ratlpr:     the pair of atoms subject to RATTLE                      
//   settlegrp:  indices for a group of atoms subject to SETTLE
//-----------------------------------------------------------------------------
static int IsInSettleGroup(nixpr ratlpr, int* settlegrp)
{
  if (settlegrp[ratlpr.atmX] >= 0 &&
      settlegrp[ratlpr.atmX] == settlegrp[ratlpr.atmY]) {
    return 1;
  }

  return 0;
}

//-----------------------------------------------------------------------------
// AddRattleMask: this function interprets any rattlemasks from the input file
//                and flags any bonds between two atoms found in the SAME mask
//                for constraint by RATTLE.            
//                                                                      
// Arguments:                                                           
//   tp:         the topology                                           
//   ratlpr:     the pair of atoms subject to RATTLE                    
//   settlegrp:  groupings for atoms whose positions are constrained by the
//               SETTLE algorithm                                   
//-----------------------------------------------------------------------------
static int AddRattleMask(prmtop *tp, nixpr* ratlbond, int* settlegrp)
{
  int i, j, nratl;
  int* atmmask;
  coord crd;

  // Bail right out if there is no rattlemask
  if (tp->rattlemask[0] == '\0') {
    return 0;
  }

  // Initialize the number of RATTLE-constrained bonds.
  // The atom IDs of RATTLE-constrained bonds populate 
  // a growing list of pairs that is allocated before  
  // this function is called.  This function is where  
  // the first entries to the list might be placed.    
  nratl = 0;

  // Parse the ambmask string
  crd = CreateCoord(tp->natom);
  atmmask = ParseAmbMask(tp->rattlemask, tp, &crd);
  for (i = 0; i < tp->natom; i++) {
    if (atmmask[i] == 0) {
      continue;
    }

    // If we're still here, atmmask[i] is 1.  Loop over  
    // all bonds controlled by this atom to see if any   
    // of the bonding partners are also part of the mask.
    for (j = 0; j < tp->BLC[i].nbond; j++) {
      if (atmmask[tp->BLC[i].BC[j].b] == 0) {
        continue;
      }

      // This bond has been called for constraint
      ratlbond[nratl].atmX = i;
      ratlbond[nratl].atmY = tp->BLC[i].BC[j].b;
      if (IsInSettleGroup(ratlbond[nratl], settlegrp) == 0) {
        nratl++;
      }
    }
  }

  // Free allocated memory
  free(atmmask);
  DestroyCoord(&crd);

  return nratl;
}

//-----------------------------------------------------------------------------
// FindBondsToHydrogen: find all bonds to hydrogen atoms.  In mdgx, these bonds
//                      are defined as those connecting an atom of mass 2.0
//                      amu or less to any other atom. 
//                                                                      
// Arguments:                                                           
//   tp:         the topology                                           
//   ratlbond:   the list of bonds to be RATTLE'd                       
//   nratl:      the number of bonds to be RATTLE'd (also returned)     
//   settlegrp:  groupings for atoms whose positions are constrained by the
//               SETTLE algorithm                                   
//-----------------------------------------------------------------------------
static int FindBondsToHydrogen(prmtop *tp, nixpr* ratlbond, int nratl,
                               int* settlegrp)
{
  int i, j, natm;

  // Loop over all atoms
  for (i = 0; i < tp->natom; i++) {

    // Count the number of RATTLE'd bonds to hydrogen
    // this atom controls.  If the atom itself is a  
    // light / hydrogen atom, then all bonds count.
    if (tp->ZNumber[i] == 1 ||
        (tp->Masses[i] < 2.0 && tp->Masses[i] > 1.0e-8)) {
      for (j = 0; j < tp->BLC[i].nbond; j++) {
        ratlbond[nratl].atmX = i;
        ratlbond[nratl].atmY = tp->BLC[i].BC[j].b;
        if (IsInSettleGroup(ratlbond[nratl], settlegrp) == 0) {
          nratl++;
        }
      }
    }
    else {
      for (j = 0; j < tp->BLC[i].nbond; j++) {
        natm = tp->BLC[i].BC[j].b;
        if (tp->ZNumber[natm] == 1 ||
            (tp->Masses[natm] < 2.0 && tp->Masses[natm] > 1.0e-8)) {
          ratlbond[nratl].atmX = i;
          ratlbond[nratl].atmY = tp->BLC[i].BC[j].b;
          if (IsInSettleGroup(ratlbond[nratl], settlegrp) == 0) {
            nratl++;
          }
        }
      }
    }
  }

  return nratl;
}

//-----------------------------------------------------------------------------
// SortRattle1: function called by quicksort for comparing the first atoms in
//              RATTLE'd bonds into increasing order of their ID numbers.
//
// Arguments:
//   bnd{A,B}  : the bonded atom pair
//-----------------------------------------------------------------------------
static int SortRattle1(const void *bndA, const void *bndB)
{
  int idA = ((nixpr*)bndA)[0].atmX;
  int idB = ((nixpr*)bndB)[0].atmX;

  if (idA < idB) {
    return -1;
  }
  else if (idA > idB) {
    return 1;
  }
  else {
    return 0;
  }
}

//-----------------------------------------------------------------------------
// SortRattle2: function called by quicksort for comparing the second   
//              atoms in RATTLE'd bonds into increasing order of their  
//              ID numbers.                                             
//
// Arguments:
//   bnd{A,B}  : the bonded atom pair
//-----------------------------------------------------------------------------
static int SortRattle2(const void *bndA, const void *bndB)
{
  int idA = ((nixpr*)bndA)[0].atmY;
  int idB = ((nixpr*)bndB)[0].atmY;

  if (idA < idB) {
    return -1;
  }
  else if (idA > idB) {
    return 1;
  }
  else {
    return 0;
  }
}

//-----------------------------------------------------------------------------
// PruneRattleMask: this function reads the "norattlemask" string and prunes
//                  the constraints involving any mask atoms.    
//                                                                      
// Arguments:                                                           
//   tp:         the topology                                           
//   ratlbond:   the list of bonds to be RATTLE'd                       
//   nratl:      the number of bonds to be RATTLE'd (also returned)     
//-----------------------------------------------------------------------------
static int PruneRattleMask(prmtop *tp, nixpr* ratlbond, int nratl)
{
  int i, j, llim, hlim, cratl;
  int* elimbnd;
  int* nrmask;
  coord crd;

  // First sort the list of RATTLE'd bonds
  // in order to prune duplicates.        
  elimbnd = (int*)calloc(nratl, sizeof(int));
  qsort(ratlbond, nratl, sizeof(nixpr), SortRattle1);
  llim = 0;
  while (llim < nratl-1) {
    hlim = llim+1;
    while (hlim < nratl && ratlbond[hlim].atmX == ratlbond[llim].atmX) {
      hlim++;
    }
    for (i = llim; i < hlim-1; i++) {
      for (j = i+1; j < hlim; j++) {
        if (ratlbond[j].atmX == ratlbond[i].atmX &&
            ratlbond[j].atmY == ratlbond[i].atmY) {
          elimbnd[j] = 1;
        }
      }
    }
    llim = hlim;
  }

  // Read the no-RATTLE mask
  if (tp->norattlemask[0] != '\0') {
    crd = CreateCoord(tp->natom);
    nrmask = ParseAmbMask(tp->rattlemask, tp, &crd);

    // Loop over all RATTLE'd bonds
    for (i = 0; i < nratl; i++) {
      if (nrmask[ratlbond[i].atmX] == 1 || nrmask[ratlbond[i].atmY] == 1) {
        elimbnd[i] = 1;
      }
    }

    // Free allocated memory
    DestroyCoord(&crd);
    free(nrmask);
  }

  // Condense the bond list
  cratl = 0;
  for (i = 0; i < nratl; i++) {
    if (elimbnd[i] == 0) {
      if (cratl < i) {
        ratlbond[cratl] = ratlbond[i];
      }
      cratl++;
    }
  }

  // Free allocated memory
  free(elimbnd);

  return cratl;
}

//-----------------------------------------------------------------------------
// SeekBondMatch: find a bond matching the key within a specified array.
//                                                                      
// Arguments:                                                           
//   A:       the constraint graph data struct                          
//   key:     the key to find/match                                     
//   lnum:    the number of the list to find the match within           
//   xory:    flag to indicate whether to match in the X or Y component 
//            (1 for X, 2 for Y)                                        
//   ordr:    the order of the match (1 for only one atom or 2 for both atoms)
//-----------------------------------------------------------------------------
static int SeekBondMatch(cnstgrp *A, nixpr key, int lnum, int xory, int ordr)
{
  int j, llim, hlim, hlf, pt, ptN;
  nixpr *svec;

  // Set the vector to look in based on the lnum argument; this
  // feature basically doubles the utility of this function.   
  if (lnum == 1) {
    svec = A->ratlbond;
  }
  else {
    svec = A->ratlbond2;
  }

  // Binary search to find a matching element of the array in
  // the slot specified by the xory argument; this ability to
  // chose between X and Y arguments also doubles utility.   
  llim = 0;
  hlim = A->nratl;
  pt = -1;
  while (llim < hlim-2 && pt < 0) {
    hlf = llim + (hlim-llim)/2;
    if (xory == 1) {
      if (key.atmX < svec[hlf].atmX) {
        hlim = hlf;
      }
      else if (key.atmX == svec[hlf].atmX) {
        pt = hlf;
      }
      else {
        llim = hlf;
      }
    }
    else {
      if (key.atmY < svec[hlf].atmY) {
        hlim = hlf;
      }
      else if (key.atmY == svec[hlf].atmY) {
        pt = hlf;
      }
      else {
        llim = hlf;
      }
    }
  }
  if (pt < 0) {
    for (j = llim; j < hlim; j++) {
      if (xory == 1 && key.atmX == svec[j].atmX) {
        pt = j;
        break;
      }
      else if (xory == 2 && key.atmY == svec[j].atmY) {
        pt = j;
        break;
      }
    }
  }

  // We're almost done if only the one atom needs to be matched.
  // The ordr argument allows us to match just the one slot, or 
  // to continue and match both slots.  If only one slot is to  
  // be matched, then we need to find the lowest index element  
  // of the array that matches, but if both slots must match    
  // then there is only one element of the array that suffices. 
  if (ordr == 1) {

    // Find the minimum number of bond that matches
    j = pt;
    if (xory == 1) {
      while (j >= 0 && svec[j].atmX == key.atmX) {
        pt = j;
        j--;
      }
    }
    if (xory == 2) {
      while (j >= 0 && svec[j].atmY == key.atmY) {
        pt = j;
        j--;
      }
    }
    return pt;
  }

  // This bond element matches in one slot, but the second slot 
  // must also match.  Count backwards and then forwards until a
  // complete match is found.                                   
  if (xory == 2) {
    if (svec[pt].atmX != key.atmX) {
      ptN = pt-1;
      while (ptN >= 0 && key.atmY == svec[ptN].atmY &&
             svec[ptN].atmX != key.atmX) {
        ptN--;
      }
      if (ptN >= 0 && key.atmY == svec[ptN].atmY &&
          svec[ptN].atmX == key.atmX) {
        pt = ptN;
      }
      else {
        ptN = pt+1;
        while (ptN < A->nratl && key.atmY == svec[ptN].atmY &&
               svec[ptN].atmX != key.atmX) {
          ptN++;
        }
        if (ptN < A->nratl && key.atmY == svec[ptN].atmY &&
            svec[ptN].atmX == key.atmX) {
          pt = ptN;
        }
        else {
          printf("SeekBondMatch >> Error.  Searched high and low, no "
                 "match for [ %4d %4d ].\n", key.atmX, key.atmY);
          exit(1);
        }
      }
    }
  }
  else {
    if (svec[pt].atmY != key.atmY) {
      ptN = pt-1;
      while (ptN >= 0 && key.atmX == svec[ptN].atmX &&
             svec[ptN].atmY != key.atmY) {
        ptN--;
      }
      if (ptN >= 0 && key.atmX == svec[ptN].atmX &&
          svec[ptN].atmY == key.atmY) {
        pt = ptN;
      }
      else {
        ptN = pt+1;
        while (ptN < A->nratl && key.atmX == svec[ptN].atmX &&
               svec[ptN].atmY != key.atmY) {
          ptN++;
        }
        if (ptN < A->nratl && key.atmX == svec[ptN].atmX &&
            svec[ptN].atmY == key.atmY) {
          pt = ptN;
        }
        else {
          printf("SeekBondMatch >> Error.  Searched high and low, no "
                 "match for [ %4d %4d ].\n", key.atmX, key.atmY);
          exit(1);
        }
      }
    }
  }

  return pt;
}

//-----------------------------------------------------------------------------
// MapBondCorrespondence: computes the correspondence between ordered lists of
//                        bonds in a ConstraintGraph struct.   
//                                                                      
// Arguments:                                                           
//   A:       the constraint graph data struct                          
//-----------------------------------------------------------------------------
static void MapBondCorrespondence(cnstgrp *A)
{
  int i, pt;

  // Loop over all bonds in the first group
  // and find their matches in the second  
  for (i = 0; i < A->nratl; i++) {
    pt = SeekBondMatch(A, A->ratlbond[i], 2, 2, 2);

    // The ith element of the first array matches the matchpt
    // element of the second array.  Record that result.     
    A->map12[i] = pt;
    A->map21[pt] = i;
  }
}

//-----------------------------------------------------------------------------
// FindRattleLinks: find RATTLE-constrained all bonds (pairs of atoms) which
//                  share one atom in common with the key.        
//                                                                      
// Arguments:                                                           
//   A:       the constraint graph data struct                          
//   key:     the key to link against                                   
//-----------------------------------------------------------------------------
static void FindRattleLinks(cnstgrp *A, nixpr key)
{
  int i, ptX, ptY;

  // Find all bonds which match the key in
  // either their first or second slots   
  ptX = SeekBondMatch(A, key, 1, 1, 1);
  ptY = SeekBondMatch(A, key, 2, 2, 1);
  i = ptX;
  while (i >= 0 && i < A->nratl && A->ratlbond[i].atmX == key.atmX) {
    if (A->rsrv[i] == 1) {
      i++;
      continue;
    }
    A->ratlgrp[A->grpsize] = A->ratlbond[i];
    A->rsrv[i] = 1;
    A->grpsize += 1;
  }
  i = ptY;
  while (i >= 0 && i < A->nratl && A->ratlbond2[i].atmY == key.atmY) {
    if (A->rsrv[A->map21[i]] == 1) {
      i++;
      continue;
    }
    A->ratlgrp[A->grpsize] = A->ratlbond2[i];
    A->rsrv[A->map21[i]] = 1;
    A->grpsize += 1;
  }
}

//-----------------------------------------------------------------------------
// GraphRattleGroup: assembles a graph of atoms participating in RATTLE 
//                   constraints.  The node of the graph with the most edges
//                   indicates which atom shall control the RATTLE constraint
//                   group.                                  
//                                                                      
// Arguments:                                                           
//   A:        the constraint graph data struct                         
//   nstatrt:  the initial edge for constructing the graph              
//-----------------------------------------------------------------------------
static void GraphRattleGroup(cnstgrp *A, int nstart)
{
  int i;
  nixpr key;
  // Initialize the group of RATTLE'd bonds
  // and the open-ended bonds in the group 
  SetIVec(A->openend, A->nratl, 0);
  SetIVec(A->nodeorder, A->nratl, 0);

  // Iteratively find all of the links for every  
  // outstanding addition to the growing list of  
  // RATTLE'd bonds in this group.  The first     
  // member of the RATTLE group is the RATTLE'd   
  // bond at position nstart of the ratlbond      
  // array (not the ratlbond2 array), and it is   
  // initially labeled as "open-ended" pending a  
  // search by the FindRattleLinks function above.
  A->openend[0] = 1;
  A->ratlgrp[0] = A->ratlbond[nstart];
  A->rsrv[nstart] = 1;
  A->grpsize = 1;
  while (ISum(A->openend, A->grpsize) > 0) {
    for (i = 0; i < A->grpsize; i++) {
      key = A->ratlgrp[i];
      FindRattleLinks(A, key);
      key.atmX = A->ratlgrp[i].atmY;
      key.atmY = A->ratlgrp[i].atmX;
      FindRattleLinks(A, key);
      A->openend[i] = 0;
    }
  }
}

//-----------------------------------------------------------------------------
// FindIDInList: this function exhaustively searches the growing list of unique
//               atom ID numbers in the RATTLE'd constraint group for the ID
//               number key and returns true (1) or false (0) if the key is or
//               is not found.                         
//                                                                      
// Arguments:                                                           
//   uniqatm:   the current list of unique atoms                        
//   nuniq:     the current length of n                                 
//   key:       the key to seach for                                    
//-----------------------------------------------------------------------------
static int FindIDInList(int* uniqatm, int nuniq, int key)
{
  int i;

  for (i = 0; i < nuniq; i++) {
    if (uniqatm[i] == key) {
      return i;
    }
  }

  return -1;
}

//-----------------------------------------------------------------------------
// AssignRattleGroup: this function will assign a group of RATTLE'd bonds to a
//                    particular atom in the topology based on the maximum
//                    number of bonds connecting to any single atom in the
//                    RATTLE-constrained group.      
//                                                                      
// Arguments:                                                           
//   tp:      the topology                                              
//   A:       the constraint group (contains a list of all bonds subject
//            to RATTLE as well as sub-lists specific to the group)     
//-----------------------------------------------------------------------------
static void AssignRattleGroup(prmtop *tp, cnstgrp *A)
{
  int i, j, maxcon, curratm, currcon, hub, nuniq;
  int* uniqatm;
  double l0;
  bondlist *blc;

  // Order the RATTLE group by the X atoms
  // and count the maximum connectivity.  
  uniqatm = (int*)malloc((A->grpsize+1)*sizeof(int));
  qsort(A->ratlgrp, A->grpsize, sizeof(nixpr), SortRattle1);
  curratm = A->ratlgrp[0].atmX;
  currcon = 0;
  maxcon = 0;
  hub = curratm;
  uniqatm[0] = curratm;
  nuniq = 1;
  for (i = 0; i < A->grpsize; i++) {
    if (A->ratlgrp[i].atmX == curratm) {
      currcon++;
    }
    if (currcon > maxcon ||
        (currcon == maxcon && tp->Masses[curratm] > tp->Masses[hub])) {
      maxcon = currcon;
      hub = curratm;
    }
    if (A->ratlgrp[i].atmX != curratm) {
      currcon = 0;
      curratm = A->ratlgrp[i].atmX;
      uniqatm[nuniq] = curratm;
      nuniq++;
    }
  }

  // Order the RATTLE group by the Y atoms
  // and count the maximum connectivity.  
  qsort(A->ratlgrp, A->grpsize, sizeof(nixpr), SortRattle2);
  curratm = A->ratlgrp[0].atmY;
  currcon = 0;
  if (FindIDInList(uniqatm, nuniq, curratm) < 0) {
    uniqatm[nuniq] = curratm;
    nuniq++;
  }
  for (i = 0; i < A->grpsize; i++) {
    if (A->ratlgrp[i].atmY == curratm) {
      currcon++;
    }
    if (currcon > maxcon ||
        (currcon == maxcon && tp->Masses[curratm] > tp->Masses[hub])) {
      maxcon = currcon;
      hub = curratm;
    }
    if (A->ratlgrp[i].atmY != curratm) {
      currcon = 0;
      curratm = A->ratlgrp[i].atmY;
      if (FindIDInList(uniqatm, nuniq, curratm) < 0) {
        uniqatm[nuniq] = curratm;
        nuniq++;
      }
    }
  }
  
  // The atom of maximum connectivity is now found.  Make  
  // the constraint group.  The constraints are codified in
  // an integer array, the first element of which is the   
  // number of constraints in the constraint group.        
  // Thereafter, the elements of the array read atom A and 
  // atom B of the bond followed by an integer denoting the
  // constrained bond length times 100,000,000. (Note that,
  // assuming a 32-bit integer, this limits the length of  
  // bond constraints to about 20A; also, atoms A and B are
  // codified according to their order of appearance in the
  // list of atoms involved in this constraint group, not  
  // their official atom IDs as found in the topology.)    
  // Finally, after the list of bond atoms and lengths, the
  // number of unique atoms in this constraint group is    
  // given followed by a list of the topology ID numbers of
  // all unique atoms in this constraint group.              
  tp->SHL[hub].exe = 2;
  tp->SHL[hub].blist = (int*)malloc((3*A->grpsize+2+nuniq)*sizeof(int));
  tp->SHL[hub].blist[0] = A->grpsize;
  for (i = 0; i < A->grpsize; i++) {
    tp->SHL[hub].blist[3*i+1] = FindIDInList(uniqatm, nuniq,
                                             A->ratlgrp[i].atmX);
    tp->SHL[hub].blist[3*i+2] = FindIDInList(uniqatm, nuniq,
                                             A->ratlgrp[i].atmY);
    l0 = -1.0;
    blc = &tp->BLC[A->ratlgrp[i].atmX];
    for (j = 0; j < blc->nbond; j++) {
      if (blc->BC[j].b == A->ratlgrp[i].atmY) {
        l0 = tp->BParam[blc->BC[j].t].l0;
        break;
      }
    }
    if (l0 < 0.0) {
      blc = &tp->BLC[A->ratlgrp[i].atmY];
      for (j = 0; j < blc->nbond; j++) {
        if (blc->BC[j].b == A->ratlgrp[i].atmX) {
          l0 = tp->BParam[blc->BC[j].t].l0;
          break;
        }
      }
    }
    tp->SHL[hub].blist[3*i+3] = l0*1.0e8 + 0.001;
    if (A->ratlgrp[i].atmX != hub) {
      tp->SHL[A->ratlgrp[i].atmX].exe = -1;
    }
    if (A->ratlgrp[i].atmY != hub) {
      tp->SHL[A->ratlgrp[i].atmY].exe = -1;
    }
  }
  tp->SHL[hub].blist[3*A->grpsize+1] = nuniq;
  for (i = 0; i < nuniq; i++) {
    tp->SHL[hub].blist[3*A->grpsize+2+i] = uniqatm[i];
  }

  // Free allocated memory
  free(uniqatm);
}

//-----------------------------------------------------------------------------
// CodifyRattle: there is now a list of bonds that must be subjected to RATTLE,
//               but it is disjoint and must be properly organized so that, in
//               such cases as multiple RATTLE'd bonds sharing the same atom,
//               the iterative constraint applications can be interwoven during
//               dynamics.        
//                                                                      
// Arguments:                                                           
//   tp:         the topology                                           
//   ratlbond:   the list of bonds to be RATTLE'd                       
//   nratl:      the number of bonds to be RATTLE'd                     
//-----------------------------------------------------------------------------
static void CodifyRattle(prmtop *tp, cnstgrp *A)
{
  int i;

  // The list of bonds must be ordered to avoid an O(N^2) search.
  // The list of bonds in therefore ordered according to the atom
  // ID numbers of the first atom in each bond.  A second copy of
  // the bond list is made and ordered according to the atom ID  
  // of the second atom in the bond.                             
  for (i = 0; i < A->nratl; i++) {
    A->ratlbond2[i] = A->ratlbond[i];
  }
  qsort(A->ratlbond, A->nratl, sizeof(nixpr), SortRattle1);
  qsort(A->ratlbond2, A->nratl, sizeof(nixpr), SortRattle2);

  // Map the correspondence between bonds
  // in the first and second lists       
  MapBondCorrespondence(A);

  // Step along the first list; identify common atoms with either 
  // the first or second atom in the RATTLE'd bond.  If there are 
  // no other constrained bonds that share either atom, assign the
  // RATTLE constraint to the heavier of the two atoms, or the    
  // first if they are of the same mass.  If there are other      
  // RATTLE'd bonds that share the same atoms, assign them all to 
  // the common atom.  Then there's the case of chains or rings.  
  // This becomes yucky, but in all cases the solution is to draw 
  // a graph that connects all of the bonds and count the number  
  // of connections to each node.                                 
  SetIVec(A->rsrv, A->nratl, 0);
  for (i = 0; i < A->nratl; i++) {
    if (A->rsrv[i] == 1) {
      continue;
    }
    GraphRattleGroup(A, i);

    // A group of RATTLE-constrained bonds has now been found.
    // Determine the atoms which the bonds attach, and then   
    // find the atom with the highest attachment count that   
    // will ultimately control the bond in the topology.      
    AssignRattleGroup(tp, A);

    // Record the maximum size of a RATTLE-constrained
    // group for reference during dynamics            
    tp->RattleGrpMax = MAX(tp->RattleGrpMax, 3*(A->grpsize+1));
  }
}

//-----------------------------------------------------------------------------
// AllocateConstraintGraph: allocate memory for a constraint graph strcuture.
//                                                                      
// Arguments:                                                           
//   n:       the maximum number of constraints that might be in any one
//            graph, and thus the size of various allocated arrays      
//-----------------------------------------------------------------------------
static cnstgrp AllocateConstraintGraph(int n)
{
  cnstgrp A;

  A.openend = (int*)malloc((n+1)*sizeof(int));
  A.nodeorder = (int*)malloc(n*sizeof(int));
  A.rsrv = (int*)malloc(n*sizeof(int));
  A.map12 = (int*)malloc(n*sizeof(int));
  A.map21 = (int*)malloc(n*sizeof(int));
  A.ratlbond = (nixpr*)malloc(n*sizeof(nixpr));
  A.ratlbond2 = (nixpr*)malloc(n*sizeof(nixpr));
  A.ratlgrp = (nixpr*)malloc(n*sizeof(nixpr));

  return A;
}

//-----------------------------------------------------------------------------
// DeallocateConstraintGraph: allocate memory for a constraint graph structure.
//                                                                      
// Arguments:                                                           
//   n:       the maximum number of constraints that might be in any one
//            graph, and thus the size of various allocated arrays      
//-----------------------------------------------------------------------------
static void DeallocateConstraintGraph(cnstgrp A)
{
  free(A.openend);
  free(A.nodeorder);
  free(A.rsrv);
  free(A.map12);
  free(A.map21);
  free(A.ratlbond);
  free(A.ratlbond2);
  free(A.ratlgrp);
}

//-----------------------------------------------------------------------------
// CheckOffExclusion: check off an exclusion that is satisfied by some bond,
//                    angle, or dihedral term.                    
//                                                                      
// Arguments:                                                           
//   tp:       the topology                                             
//   aatm:     the first atom of the bond / angle / dihedral            
//   batm:     the second atom of the bond / angle / diehdral           
//   exclscr:  scratch array of exclusions                              
//-----------------------------------------------------------------------------
static void CheckOffExclusion(prmtop *tp, int aatm, int batm, int* exclscr)
{
  int i;

  for (i = tp->ConExcl[aatm]; i < tp->ConExcl[aatm+1]; i++) {
    if (tp->ExclList[i] == batm) {
      exclscr[i] = 1;
    }
  }
  for (i = tp->ConExcl[batm]; i < tp->ConExcl[batm+1]; i++) {
    if (tp->ExclList[i] == aatm) {
      exclscr[i] = 1;
    }
  }
}

//-----------------------------------------------------------------------------
// FulfillExclusions: make sure that all exclusions are covered by 1:2, 1:3,
//                    or 1:4 bonded interactions.  If not, allocate the pair
//                    elimination list and fill it out with any exclusions
//                    which would otherwise be missed.  A flag in the topology
//                    (ExclMarked) signifies that certain exclusions must be
//                    managed in this manner.
//                                                                      
// Arguments:                                                           
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
static void FulfillExclusions(prmtop *tp)
{
  int i, j, k, elj, elimtype;
  int* exclscr;
  bondlist *BLi;
  angllist *ALi;
  dihelist *HLi;
  auxelim *tAEi;

  // Loop over all atoms, allocating a list of  
  // exclusions and checking off whether they   
  // are handled by bonds, angles, or dihedrals.
  exclscr = (int*)calloc(tp->ConExcl[tp->natom], sizeof(int));
  
  // Blanks in the exclusion list do not count as exclusions, 
  // they merely mean that a particular atom has no exclusions
  for (i = 0; i < tp->ConExcl[tp->natom]; i++) {
    if (tp->ExclList[i] == -1) {
      exclscr[i] = 1;
    }
  }

  // Bonds, angles, and dihedrals
  for (i = 0; i < tp->natom; i++) {
    BLi = &tp->BLC[i];
    for (j = 0; j < BLi->nbond; j++) {
      CheckOffExclusion(tp, BLi->BC[j].a, BLi->BC[j].b, exclscr);
    }
    ALi = &tp->ALC[i];
    for (j = 0; j < ALi->nangl; j++) {
      CheckOffExclusion(tp, ALi->AC[j].a, ALi->AC[j].c, exclscr);
    }
    HLi = &tp->HLC[i];
    for (j = 0; j < HLi->ndihe; j++) {
      if (HLi->HC[j].impr == 0) {
        CheckOffExclusion(tp, HLi->HC[j].a, HLi->HC[j].d, exclscr);
      }
    }
  }

  // Exclusions that already exist for added extra points
  if (tp->EPInserted == 1) {
    for (i = 0; i < tp->natom; i++) {
      tAEi = &tp->ElimPair[i];
      for (j = 0; j < tAEi->n11; j++) {
        CheckOffExclusion(tp, tAEi->list11[j].atmX, tAEi->list11[j].atmY,
                          exclscr);
      }
      for (j = 0; j < tAEi->n12; j++) {
        CheckOffExclusion(tp, tAEi->list12[j].atmX, tAEi->list12[j].atmY,
                          exclscr);
      }
      for (j = 0; j < tAEi->n13; j++) {
        CheckOffExclusion(tp, tAEi->list13[j].atmX, tAEi->list13[j].atmY,
                          exclscr);
      }
      for (j = 0; j < tAEi->n14; j++) {
        CheckOffExclusion(tp, tAEi->list14[j].atmX, tAEi->list14[j].atmY,
                          exclscr);
      }
    }
  }

  // Scan for any exclusions unaccounted for
  tp->ExclMarked = 0;
  i = 0;
  while (i < tp->natom && tp->ExclMarked == 0) {
    for (j = tp->ConExcl[i]; j < tp->ConExcl[i+1]; j++) {
      if (exclscr[j] == 0) {
        tp->ExclMarked = 1;
        break;
      }
    }
    i++;
  }

  // If we must perform exclusions, allocate 
  // the elimination list and enumerate them.
  // The first atom of each exclusion counts 
  // as the controlling atom.                
  if (tp->ExclMarked == 1 && tp->EPInserted == 0) {
    AllocateElimList(tp);
  }
  for (i = 0; i < tp->natom; i++) {
    for (j = tp->ConExcl[i]; j < tp->ConExcl[i+1]; j++) {
      if (exclscr[j] == 0) {

	// The proximity is important.  1:1, 1:2, and 1:3 are
	// all strict eliminations, but 1:4 are attenuations.
	elimtype = 0;
	elj = tp->ExclList[j];
	for (k = 0; k < tp->nb1234[i].n11; k++) {
	  if (tp->nb1234[i].L11[k] == elj) {
	    elimtype = 1;
	  }
	}
	for (k = 0; k < tp->nb1234[i].n12; k++) {
	  if (tp->nb1234[i].L12[k] == elj) {
	    elimtype = 2;
	  }
	}
	for (k = 0; k < tp->nb1234[i].n13; k++) {
	  if (tp->nb1234[i].L13[k] == elj) {
	    elimtype = 3;
	  }
	}
	for (k = 0; k < tp->nb1234[i].n14; k++) {
	  if (tp->nb1234[i].L14[k] == elj) {
	    elimtype = 4;
	  }
	}
	if (elimtype == 0) {
	  printf("FulfillExclusions >> Error.  Exclusion %d - %d cannot be "
		 "traced.\n", i + 1, elj + 1);
	  printf("FulfillExclusions >> Atoms: %4.4s %4.4s in residues "
		 "%4.4s %4.4s\n", &tp->AtomNames[4*i], &tp->AtomNames[4*elj],
		 &tp->ResNames[4*LocateResID(tp, i, 0, tp->nres)],
		 &tp->ResNames[4*LocateResID(tp, elj, 0, tp->nres)]);
	  exit(1);
	}
	AddNewElimPair(tp, i, i, elj, elimtype);
      }
    }
  }

  // Free allocated memory
  free(exclscr);
}

//-----------------------------------------------------------------------------
// FindRattleBonds: find all bonds that require RATTLE constraints.     
//                                                                      
// Arguments:                                                           
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
static void FindRattleBonds(prmtop *tp)
{
  int i, maxcnst;
  int* settlegrp;
  cnstgrp mygrf;
  
  // Bail out if RATTLE is not wanted
  if (tp->rattle == 0) {
    return;
  }

  // Make labels of all SETTLE groups to avoid double-constraining
  settlegrp = LabelSettleGroups(tp);

  // Allocate memory to store all possible bond constraints
  maxcnst = 0;
  for (i = 0; i < tp->natom; i++) {
    maxcnst += tp->BLC[i].nbond;
  }
  mygrf = AllocateConstraintGraph(maxcnst);

  // First, any special RATTLE bonds get added
  mygrf.nratl = AddRattleMask(tp, mygrf.ratlbond, settlegrp);

  // Next, seek and RATTLE all hydrogens
  mygrf.nratl = FindBondsToHydrogen(tp, mygrf.ratlbond, mygrf.nratl,
                                    settlegrp);

  // Finally, prune constraints in the norattlemask string
  mygrf.nratl = PruneRattleMask(tp, mygrf.ratlbond, mygrf.nratl);

  // Codify the RATTLE'd bonds into the prmtop constraints
  CodifyRattle(tp, &mygrf);

  // Free allocated memory
  free(settlegrp);
  DeallocateConstraintGraph(mygrf);
}

//-----------------------------------------------------------------------------
// RebuildExclusionList: rebuild an exclusion list based on connectivity
//                       information.                                   
//                                                                      
// Arguments:                                                           
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
static void RebuildExclusionList(prmtop *tp)
{
  int i, j, k, totalexcl, atmexcl;
  int nswap, ibuff;
  
  // Reallocate the exclusion indexing arrays
  tp->NExcl = (int*)realloc(tp->NExcl, tp->natom*sizeof(int));
  tp->ConExcl = (int*)realloc(tp->ConExcl, (tp->natom+1)*sizeof(int));
  
  // Count the number of exclusions
  totalexcl = 0;
  tp->ConExcl[0] = 0;
  for (i = 0; i < tp->natom; i++) {
    atmexcl = 0;
    for (j = 0; j < tp->nb1234[i].n11; j++) {
      if (tp->nb1234[i].L11[j] > i) {
        atmexcl++;
      }
    }
    for (j = 0; j < tp->nb1234[i].n12; j++) {
      if (tp->nb1234[i].L12[j] > i) {
        atmexcl++;
      }
    }
    for (j = 0; j < tp->nb1234[i].n13; j++) {
      if (tp->nb1234[i].L13[j] > i) {
        atmexcl++;
      }
    }
    for (j = 0; j < tp->nb1234[i].n14; j++) {
      if (tp->nb1234[i].L14[j] > i) {
        atmexcl++;
      }
    }
    tp->NExcl[i] = atmexcl;
    if (atmexcl == 0) {
      atmexcl = 1;
    }
    tp->ConExcl[i+1] = tp->ConExcl[i] + atmexcl;
    totalexcl = totalexcl+atmexcl;
  }

  // Update the total number of exclusions in case it is needed for printing
  tp->tnexcl = totalexcl;

  // Allocate the exclusions array
  tp->ExclList = (int*)realloc(tp->ExclList, totalexcl*sizeof(int));
  totalexcl = 0;
  for (i = 0; i < tp->natom; i++) {
    atmexcl = totalexcl;
    for (j = 0; j < tp->nb1234[i].n11; j++) {
      if (tp->nb1234[i].L11[j] > i) {
        tp->ExclList[totalexcl] = tp->nb1234[i].L11[j];
        totalexcl++;
      }
    }
    for (j = 0; j < tp->nb1234[i].n12; j++) {
      if (tp->nb1234[i].L12[j] > i) {
        tp->ExclList[totalexcl] = tp->nb1234[i].L12[j];
        totalexcl++;
      }
    }
    for (j = 0; j < tp->nb1234[i].n13; j++) {
      if (tp->nb1234[i].L13[j] > i) {
        tp->ExclList[totalexcl] = tp->nb1234[i].L13[j];
        totalexcl++;
      }
    }
    for (j = 0; j < tp->nb1234[i].n14; j++) {
      if (tp->nb1234[i].L14[j] > i) {
        tp->ExclList[totalexcl] = tp->nb1234[i].L14[j];
        totalexcl++;
      }
    }

    // If this atom has no exclusions it gets a -1 entry
    if (atmexcl == totalexcl) {
      tp->ExclList[totalexcl] = -1;
      totalexcl++;
    }
    
    // Sort the exclusions for this atom (go bubble sort!)
    for (j = tp->ConExcl[i]; j < tp->ConExcl[i+1]-2; j++) {
      nswap = 0;
      for (k = j; k < tp->ConExcl[i+1]-1; k++) {
	if (tp->ExclList[k] > tp->ExclList[k+1]) {
	  SWAP(tp->ExclList[k], tp->ExclList[k+1], ibuff);
	  nswap++;
	}
      }
      if (nswap == 0) {
        break;
      }
      else {
	j--;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// FinalizeNewEP: finalize a prmtop struct after adding new extra points by
//                tidying up a few things.                           
//                                                                      
// Arguments:                                                           
//   tp:         the topology structure                                 
//-----------------------------------------------------------------------------
static void FinalizeNewEP(prmtop* tp)
{
  int i, j;

  // Record the numbers of extra points' first frame atoms in an 
  // atom-by-atom array.                                         
  tp->FR1Idx = (lgrp*)calloc(tp->natom, sizeof(lgrp));
  for (i = 0; i < tp->nxtrapt; i++) {
    tp->FR1Idx[tp->xtrapts[i].fr1].natom += 1;
  }
  for (i = 0; i < tp->natom; i++) {
    j = (tp->FR1Idx[i].natom == 0) ? 1 : tp->FR1Idx[i].natom;
    tp->FR1Idx[i].atoms = (int*)malloc(j*sizeof(int));
    tp->FR1Idx[i].natom = 0;
  }
  for (i = 0; i < tp->nxtrapt; i++) {
    j = tp->xtrapts[i].fr1;
    tp->FR1Idx[j].atoms[tp->FR1Idx[j].natom] = i;
    tp->FR1Idx[j].natom += 1;
  }

  // Bail out if no extra points are present
  if (tp->nxtrapt == 0) {
    return;
  }

  // Use the connectivity array to rebuild the exclusion list
  RebuildExclusionList(tp);
}

//-----------------------------------------------------------------------------
// PadArrayByHalf: fill an array of twice the length of the input vector with
//                 [ second half ] [ full array ] [ first half ] to complete a
//                 double image of the input vector.  The input vector is a
//                 row or column of a double-precision matrix, as indicated by
//                 a flag.
//
// Arguments:
//   D:        the double-precision matrix
//   i:        the index of the row or column to read.
//   rc:       flag to indicate whether a row (0) or column (1) is to be read
//   V:        pre-allocated vector to fill
//-----------------------------------------------------------------------------
static void PadArrayByHalf(dmat *D, int i, int rc, double* V)
{
  int j, nval, hnval;

  nval = (rc == 0) ? D->row : D->col;
  hnval = nval / 2;
  for (j = 0; j < nval; j++) {
    if (rc == 0) {
      V[j + hnval] = D->map[i][j];
      if (j < hnval) {
	V[hnval + nval + j] = D->map[i][j];
      }
      else {
	V[j - hnval       ] = D->map[i][j];
      }
    }
    else {
      V[j + hnval] = D->map[j][i];
      if (j < hnval) {
	V[hnval + nval + j] = D->map[j][i];
      }
      else {
	V[j - hnval       ] = D->map[j][i];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// GenerateCubicSpline: the first stage of cubic spline production.  This
//                      method ensures continuous derivatives at the boundaries
//                      between all segments, even though I suspect there is a
//                      better way to fit the splines with smaller derivatives.
//
// Arguments:
//   npts:      the number of points in the spline (twice the number of
//              points in a CMAP along any one dimension)
//   stepsize:  the step size for the CMAP, in degrees
//   tvec:      vector of energy values from one strip of a CMAP
//   tvec2:     first pass at filling the cubic spline
//-----------------------------------------------------------------------------
static void GenerateCubicSpline(int npts, double stepsize, double* tvec,
				double* tvec2)
{
  int i;
  double p;
  double* tmp;

  tmp = (double*)calloc(npts, sizeof(double));
  SetDVec(tvec2, npts, 0.0);
  for (i = 1; i < npts-1; i++) {
    p = 0.5 * tvec2[i-1] + 2.0;
    tvec2[i] = -0.5 / p;
    tmp[i] = (tvec[i+1] - 2.0*tvec[i] + tvec[i-1]) / stepsize;
    tmp[i] = ((3.0 * tmp[i] / stepsize) - 0.5*tmp[i-1]) / p;
  }
  tvec2[npts-1] = 0.0;
  for (i = npts-2; i >= 0; i--) {
    tvec2[i] = (tvec2[i] * tvec2[i+1]) + tmp[i];
  }

  // Free allocated memory
  free(tmp);
}

//-----------------------------------------------------------------------------
// EvaluateCubicSpline: the second stage of cubic spline production.
//
// Arguments:
//   npts:      the number of points in the spline (only the number of points
//              in a CMAP along any one dimension, not twice)
//   stepsize:  the CMAP discretization, in degrees
//   tvec:      the energy vector
//   tvec2:     first stage spline production vector
//   spln:      the final product cubic spline, again pre-allocated
//-----------------------------------------------------------------------------
static void EvaluateCubicSpline(int npts, double stepsize, double* tvec,
                                double* tvec2, double* spln)
{
  int i, ip;
  
  SetDVec(spln, npts, 0.0);
  for (i = 0; i < npts; i++) {
    ip = i + npts/2;
    spln[i] = ((tvec[ip+1] - tvec[ip]) / stepsize) -
      (stepsize * ((1.0 / 3.0)*tvec2[ip] + (1.0 / 6.0)*tvec2[ip+1]));
  }
}

//-----------------------------------------------------------------------------
// LinearSplineCoeff: solve the cubic spline coefficients using a linear
//                    system of equations.  This gives an exact solution to
//                    the problem of continuous derivatives, whereas the pmemd
//                    method is only approximate.
//
// Arguments:
//   npts:      the number of points in the spline (only the number of points
//              in a CMAP along any one dimension, not twice)
//   tvec:      the energy (or first partial derivative) vector
//   dspln:     the final product derivatives, again pre-allocated
//   A:         pre-allocated matrix to set up and solve the system
//   b:         pre-allocated array to set up and solve the system, Ax = b
//-----------------------------------------------------------------------------
static void LinearSplineCoeff(int npts, double* tvec, double* dspln, dmat *A,
			      double* b)
{
  int i, ip;
  double stepsize, phi, phip, phix;

  // Initialize the matrix and vectors
  SetDVec(A->data, A->row * A->col, 0.0);
  SetDVec(b, A->row, 0.0);
  
  // Set up the system of equations
  stepsize = 360.0 / (double)npts;  
  for (i = 0; i < npts; i++) {
    phi = stepsize * (double)i;
    phip = phi + stepsize;
    phix = phip;
    ip = i + 1;
    if (ip == npts) {
      ip = 0;
      phix = 0.0;
    }

    // Energy values for the piecewise spline must meet
    A->map[4*i    ][4*i    ] = phi * phi * phi;
    A->map[4*i    ][4*i + 1] = phi * phi;
    A->map[4*i    ][4*i + 2] = phi;
    A->map[4*i    ][4*i + 3] = 1.0;
    b[4*i    ]               = tvec[i];
    A->map[4*i + 1][4*i    ] = phip * phip * phip;
    A->map[4*i + 1][4*i + 1] = phip * phip;
    A->map[4*i + 1][4*i + 2] = phip;
    A->map[4*i + 1][4*i + 3] = 1.0;
    b[4*i + 1]               = tvec[ip];

    // First derivatives of each spline segment must
    // be equal at the boundary between two segments
    A->map[4*i + 2][4*i     ] =  3.0 * phip * phip;
    A->map[4*i + 2][4*i  + 1] =  2.0 * phip;
    A->map[4*i + 2][4*i  + 2] =  1.0;
    A->map[4*i + 2][4*ip    ] = -3.0 * phix * phix;
    A->map[4*i + 2][4*ip + 1] = -2.0 * phix;
    A->map[4*i + 2][4*ip + 2] = -1.0;
    b[4*i + 2]                = 0.0;

    // Second derivatives of each spline segment must
    // be equal at the boundary between two segments,
    // despite the fact that second derivatives don't
    // carry over into the bicubic spline table.
    A->map[4*i + 3][4*i     ] =  6.0 * phip;
    A->map[4*i + 3][4*i  + 1] =  2.0;
    A->map[4*i + 3][4*ip    ] = -6.0 * phix;
    A->map[4*i + 3][4*ip + 1] = -2.0;
    b[4*i + 3]                = 0.0;
  }

  // Solve the system
  AxbQRRxc(*A, b, 0);
  BackSub(*A, b);
  
  // Translate the first derivatives into the returned array
  for (i = 0; i < npts; i++) {
    phi = stepsize * (double)i;
    dspln[i] = ((3.0*b[4*i]*phi) + (2.0*b[4*i + 1]))*phi + b[4*i + 2];
  }
}

//-----------------------------------------------------------------------------
// ComputeLinearCmapDerivatives: function for computing derivatives of a CMAP
//                               by a linear system of equations.  This gives
//                               results that are nearly identical to what's
//                               in pmemd and CHARMM, but are in fact more
//                               correct (there are no boundary artifacts in
//                               this method).
//
// Arguments:
//   esrf:    the energy surface (input)
//   dphi:    the first dimension partial derivative (returned)
//   dpsi:    the second dimension partial derivative (returned)
//   d2pp:    the mixed partial derivative
//   
//-----------------------------------------------------------------------------
void ComputeLinearCmapDerivatives(dmat *esrf, dmat *dphi, dmat *dpsi,
				  dmat *d2pp)
{
  int i, j, npts;
  double stepsize;
  double* b;
  double* tvec;
  double* spln;
  dmat A;

  // Allocate scratch space
  npts = esrf->row;
  stepsize = 360.0 / (double)npts;
  A = CreateDmat(4 * npts, 4 * npts, 0);
  b = (double*)malloc(4 * npts * sizeof(double));
  tvec = (double*)malloc(npts * sizeof(double));
  spln = (double*)malloc(npts * sizeof(double));
  
  // Compute dPhi one column at a time
  for (i = 0; i < npts; i++) {
    for (j = 0; j < npts; j++) {
      tvec[j] = esrf->map[j][i];
    }
    LinearSplineCoeff(npts, tvec, spln, &A, b);
    for (j = 0; j < npts; j++) {
      dphi->map[j][i] = spln[j] * stepsize;
    }
  }

  // Compute dPsi on row at a time
  for (i = 0; i < npts; i++) {
    for (j = 0; j < npts; j++) {
      tvec[j] = esrf->map[i][j];
    }
    LinearSplineCoeff(npts, tvec, spln, &A, b);
    for (j = 0; j < npts; j++) {
      dpsi->map[i][j] = spln[j] * stepsize;
    }
  }

  // Compute dPhi/dPsi
  for (i = 0; i < npts; i++) {
    for (j = 0; j < npts; j++) {
      tvec[j] = dphi->map[i][j];
    }
    LinearSplineCoeff(npts, tvec, spln, &A, b);
    for (j = 0; j < npts; j++) {
      d2pp->map[i][j] = spln[j] * stepsize;
    }
  }

  // Free allocated memory
  DestroyDmat(&A);
  free(b);
  free(tvec);
  free(spln);
}

//-----------------------------------------------------------------------------
// Compute1DStencil: compute the one-dimensional stencil for expressing the
//                   first partial derivative of a piecewise cubic spline as a
//                   linear combination of the other points on the spline.
//
// Arguments:
//   npts:    the number of points in the spline.  That's all that matters.
//-----------------------------------------------------------------------------
dmat Compute1DStencil(int npts, double* tvec)
{
  int i, j, k, ip, pivot;
  double stepsize, tmpdb, phi, phip, phix, factor, vrmsd, maxval;
  double* b;
  double* bcpy;
  double* blin;
  dmat A, Acpy, T, dT, meandT;
  
  // Initialize a mock matrix and a tracking matrix.  The tracking matrix
  // records what multiple of the value of each grid point in the spline
  // contributes to the evolving solution vector as the matrix undergoes
  // Gaussian elimination.
  A = CreateDmat(4 * npts, 4 * npts, 0);
  T = CreateDmat(4 * npts, npts, 0);
  b = (double*)calloc(4 * npts, sizeof(double));
  stepsize = 360.0 / (double)npts;
  for (i = 0; i < npts; i++) {
    phi = 0.0;
    phip = phi + stepsize;
    ip = i + 1;
    if (ip == npts) {
      ip = 0;
    }

    // Energy values for the piecewise spline must meet
    A.map[4*i    ][4*i    ] = phi * phi * phi;
    A.map[4*i    ][4*i + 1] = phi * phi;
    A.map[4*i    ][4*i + 2] = phi;
    A.map[4*i    ][4*i + 3] = 1.0;
    T.map[4*i    ][i      ] = 1.0;
    A.map[4*i + 1][4*i    ] = phip * phip * phip;
    A.map[4*i + 1][4*i + 1] = phip * phip;
    A.map[4*i + 1][4*i + 2] = phip;
    A.map[4*i + 1][4*i + 3] = 1.0;
    T.map[4*i + 1][ip     ] = 1.0;
    
    // First derivatives of each spline segment must
    // be equal at the boundary between two segments
    A.map[4*i + 2][4*i     ] =  3.0 * phip * phip;
    A.map[4*i + 2][4*i  + 1] =  2.0 * phip;
    A.map[4*i + 2][4*i  + 2] =  1.0;
    A.map[4*i + 2][4*ip    ] = -3.0 * phi * phi;
    A.map[4*i + 2][4*ip + 1] = -2.0 * phi;
    A.map[4*i + 2][4*ip + 2] = -1.0;

    // Second derivatives of each spline segment must
    // be equal at the boundary between two segments,
    // despite the fact that second derivatives don't
    // carry over into the bicubic spline table.
    A.map[4*i + 3][4*i     ] =  6.0 * phip;
    A.map[4*i + 3][4*i  + 1] =  2.0;
    A.map[4*i + 3][4*ip    ] = -6.0 * phi;
    A.map[4*i + 3][4*ip + 1] = -2.0;

    // Set up the solution vector just as a check on the final result
    b[4*i    ]               = tvec[i];
    b[4*i + 1]               = tvec[ip];
  }

  // Create cop[ies of the matrix and vector to solve by QR decomposition
  CopyDmat(&Acpy, &A, 0);
  bcpy = CpyDVec(b, 4 * npts);

  // Proceed with Gaussian elimination, make a lower triangular matrix
  for (i = 0; i < A.row; i++) {
    
    // Find the row with the largest (absolute) non-zero value in this column
    maxval = fabs(A.map[i][i]);
    pivot = i;
    for (j = i+1; j < A.row; j++) {
      if (fabs(A.map[j][i]) > maxval) {
	maxval = fabs(A.map[j][i]);
	pivot = j;
      }
    }
    
    // Make this the ith row
    if (pivot != i) {
      for (k = 0; k < A.col; k++) {
	SWAP(A.map[pivot][k], A.map[i][k], tmpdb);
      }
      for (k = 0; k < T.col; k++) {
	SWAP(T.map[pivot][k], T.map[i][k], tmpdb);
      }
      SWAP(b[pivot], b[i], tmpdb);
    }

    // Kill everything below the ith row
    for (j = i+1; j < A.row; j++) {
      if (fabs(A.map[j][i]) > 0.0) {
        factor = A.map[j][i] / A.map[i][i];
	for (k = i; k < A.col; k++) {
	  A.map[j][k] -= factor * A.map[i][k];
	}
	for (k = 0; k < T.col; k++) {
	  T.map[j][k] -= factor * T.map[i][k];
	}
	b[j] -= factor * b[i];
      }
    }
  }

  // Make a diagonal matrix
  for (i = A.row-1; i >= 0; i--) {

    // Normalize each row
    factor = 1.0 / A.map[i][i];
    for (j = 0; j < T.col; j++) {
      T.map[i][j] *= factor;
    }
    A.map[i][i] = 1.0;
    b[i] *= factor;

    // Use it to kill everything above
    for (j = i-1; j >= 0; j--) {
      factor = A.map[j][i];
      A.map[j][i] -= factor * A.map[i][i];
      for (k = 0; k < T.col; k++) {
	T.map[j][k] -= factor * T.map[i][k];
      }
      b[j] -= factor * b[i];
    }
  }

  // Perform QR decomposition to compare results (check the numerical
  // stability of our Gaussian elimination)
  AxbQRRxc(Acpy, bcpy, 0);
  BackSub(Acpy, bcpy);
  vrmsd = VecRMSD(b, bcpy, 4 * npts);
  if (vrmsd > 1.0e-6) {
    printf("Compute1DStencil >> Warning.  Error in Partial-Pivot Gaussian "
           "Elimination vs\nCompute1DStencil >> QR decomposition is %12.8lf.  "
	   "This could indicate an invalid\nCompute1DStencil >> formula "
           "estimate of CMAP derivative dependencies.\n", vrmsd);
  }
  blin = (double*)calloc(4 * npts, sizeof(double));
  for (i = 0; i < 4*npts; i++) {
    for (j = 0; j < npts; j++) {
      blin[i] += T.map[i][j] * tvec[j];
    }
  }
  vrmsd = VecRMSD(b, blin, 4 * npts);
  if (vrmsd > 1.0e-6) {
    printf("Compute1DStencil >> Warning.  Error in forumlaic derivative "
           "derivation vs\nCompute1DStencil >> linear solution is %12.8lf.  "
	   "This could indicate an invalid\nCompute1DStencil >> formula "
           "estimate of CMAP derivative dependencies.\n", vrmsd);
  }
      
  // Take the formulas for the first derivative at every grid point and
  // average them (they should be the same, but shifted).  The average
  // formula for the first derivative, which should be extremely close
  // to any one of the values in the formula matrix, is the stencil.
  // Scale the stencil by the step size to get the derivative into the
  // proper units of kcal/mol - degree.
  dT = CreateDmat(npts, npts, 0);
  for (i = 0; i < npts; i++) {
    for (j = 0; j < npts; j++) {
      k = i + j;
      k -= (k >= npts) * npts;
      dT.map[i][j] = T.map[4*i + 2][k];
    }
  }
  meandT = CreateDmat(1, npts, 0);
  for (j = 0; j < dT.col; j++) {
    for (i = 0; i < dT.row; i++) {
      meandT.data[j] += dT.map[i][j];
    }
    meandT.data[j] *= stepsize / dT.row;
  }
  
  // Free allocated memory
  DestroyDmat(&A);
  DestroyDmat(&Acpy);
  DestroyDmat(&T);
  DestroyDmat(&dT);
  free(b);
  free(bcpy);
  free(blin);

  // Return the stencil
  return meandT;
}

//-----------------------------------------------------------------------------
// ComputeCmapDerivatives: compute the derivatives for a particular CMAP.
//                         This unpacks the energy surface matrix to get all
//                         the information it needs.  This function assumes
//                         that CMAPs span 360 degrees in both directions.
//
// Arguments:
//   tcmap:   the CMAP struct to operate on
//-----------------------------------------------------------------------------
static void ComputeCmapDerivatives(cmapdef *tcmap)
{
  int i, j, k, kp, m, mp, npts, hnpts;
  double stepsize, dE, dPhiEst, dPsiEst, dPhiPsiEst;
  double* tvec;
  double* tvec2;
  double* spln;

  // Convenient numbers
  npts = tcmap->esrf.row;
  if (npts % 2 != 0) {
    printf("ComputeCmapDerivatives >> Error.  The dimension of any CMAP must "
           "be a\nComputeCmapDerivatives >> multiple of 2.\n");
    exit(1);
  }
  npts = tcmap->esrf.row;
  hnpts = npts / 2;
  stepsize = 360.0 / tcmap->esrf.row;

  // Allocate scratch arrays
  tvec  = (double*)malloc(2 * npts * sizeof(double));
  tvec2 = (double*)malloc(2 * npts * sizeof(double));
  spln  = (double*)malloc(npts * sizeof(double));

  // Compute dPhi
  for (i = 0; i < npts; i++) {
    PadArrayByHalf(&tcmap->esrf, i, 1, tvec);
    GenerateCubicSpline(2 * npts, stepsize, tvec, tvec2);
    EvaluateCubicSpline(npts, stepsize, tvec, tvec2, spln);
    for (j = 0; j < npts; j++) {
      tcmap->dphi.map[j][i] = spln[j] * stepsize;
    }
  }

  // Compute dPsi
  for (i = 0; i < npts; i++) {
    PadArrayByHalf(&tcmap->esrf, i, 0, tvec);
    GenerateCubicSpline(2 * npts, stepsize, tvec, tvec2);
    EvaluateCubicSpline(npts, stepsize, tvec, tvec2, spln);
    for (j = 0; j < npts; j++) {
      tcmap->dpsi.map[i][j] = spln[j] * stepsize;
    }
  }

  // Compute dPhi/dPsi
  for (i = 0; i < npts; i++) {
    PadArrayByHalf(&tcmap->dphi, i, 0, tvec);
    GenerateCubicSpline(2 * npts, stepsize, tvec, tvec2);
    EvaluateCubicSpline(npts, stepsize, tvec, tvec2, spln);
    for (j = 0; j < npts; j++) {
      tcmap->d2pp.map[i][j] = spln[j] * stepsize;
    }
  }

  // Compute dPhi by solving a linear system of equations
  ComputeLinearCmapDerivatives(&tcmap->esrf, &tcmap->Ldphi, &tcmap->Ldpsi,
			       &tcmap->Ld2pp);
  
  // Check for differences in the exact and approximate solutions.
  // The approximate solution is used by pmemd.
  tcmap->muePhi = 0.0;
  tcmap->muePsi = 0.0;
  tcmap->mue2PP = 0.0;
  for (i = 0; i < npts; i++) {
    for (j = 0; j < npts; j++) {
      tcmap->muePhi += fabs(tcmap->Ldphi.map[i][j] - tcmap->dphi.map[i][j]);
      tcmap->muePsi += fabs(tcmap->Ldpsi.map[i][j] - tcmap->dpsi.map[i][j]);
      tcmap->mue2PP += fabs(tcmap->Ld2pp.map[i][j] - tcmap->d2pp.map[i][j]);
    }
  }
  tcmap->muePhi /= (double)(npts * npts);
  tcmap->muePsi /= (double)(npts * npts);
  tcmap->mue2PP /= (double)(npts * npts);
  
  // Make the influence stencil appropriate for this CMAP.  The influence
  // stencil expresses dPhi or dPsi as a linear combination of the grid
  // points in the column or row.  It can be read as a giant formula.
  // The outer product of the stencil with itself gives the influence
  // stencil for computing the dPhi/dPsi mixed partial derivative.
  for (i = 0; i < npts; i++) {
    tvec[i] = tcmap->esrf.map[i][0];
  }
  tcmap->stencil = Compute1DStencil(npts, tvec);

  // Check dPhi, dPsi, and dPhi_dPsi estimates against the exact values
  tcmap->mueFPhi = 0.0;
  tcmap->mueFPsi = 0.0;
  tcmap->mueF2PP = 0.0;
  for (i = 0; i < npts; i++) {
    for (j = 0; j < npts; j++) {
      dPhiEst = 0.0;
      for (k = 0; k < npts; k++) {
	kp = i + k;
	kp -= (kp >= npts) * npts;
	dPhiEst += tcmap->esrf.map[kp][j] * tcmap->stencil.data[k];
      }
      tcmap->mueFPhi += fabs(tcmap->Ldphi.map[i][j] - dPhiEst);
    }
  }
  for (i = 0; i < npts; i++) {
    for (j = 0; j < npts; j++) {
      dPsiEst = 0.0;
      for (k = 0; k < npts; k++) {
	kp = j + k;
	kp -= (kp >= npts) * npts;
	dPsiEst += tcmap->esrf.map[i][kp] * tcmap->stencil.data[k];
      }
      tcmap->mueFPsi += fabs(tcmap->Ldpsi.map[i][j] - dPsiEst);
    }
  }
  for (i = 0; i < npts; i++) {
    for (j = 0; j < npts; j++) {
      dPhiPsiEst = 0.0;
      for (k = 0; k < npts; k++) {
	kp = i + k;
	kp -= (kp >= npts) * npts;
	for (m = 0; m < npts; m++) {
	  mp = j + m;
	  mp -= (mp >= npts) * npts;
          dPhiPsiEst += tcmap->esrf.map[kp][mp] * tcmap->stencil.data[k] *
                        tcmap->stencil.data[m];
	}
      }
      tcmap->mueF2PP += fabs(tcmap->Ld2pp.map[i][j] - dPhiPsiEst);
    }
  }
  tcmap->mueFPhi /= (double)(tcmap->esrf.row * tcmap->esrf.col);
  tcmap->mueFPsi /= (double)(tcmap->esrf.row * tcmap->esrf.col);
  tcmap->mueF2PP /= (double)(tcmap->esrf.row * tcmap->esrf.col);
  
  // Free allocated memory
  free(tvec);
  free(tvec2);
  free(spln);
}

//-----------------------------------------------------------------------------
// ReadCmapInfo: read CMAP information from an Amber topology file.  Prepare
//               the eneregy and derivative matrices for bicubic spline
//               interpolation.
//
// Arguments:
//   Ctop:       character matrix encoding the prmtop file
//   npline:     the line to start on
//   tp:         pointer to the developing topology
//-----------------------------------------------------------------------------
static void ReadCmapInfo(cmat *Ctop, int npline, prmtop *tp)
{
  int i, j, k, m, nmaps, dataoffset;
  int *itmp;
  char* code;
  
  itmp = ReadXiL(Ctop, &npline, 2, 2, 8);
  tp->ncmap     = itmp[0];
  tp->ncmapsurf = itmp[1];
  free(itmp);
  npline = ScanToFlag(Ctop, "CMAP_RESOLUTION", npline, 1);
  itmp = ReadXiL(Ctop, &npline, tp->ncmapsurf, 20, 4);
  tp->MParam = (cmapdef*)malloc(tp->ncmapsurf * sizeof(cmapdef));
  code = (char*)malloc(MAXNAME*sizeof(char));
  for (i = 0; i < tp->ncmapsurf; i++) {
    sprintf(code, "CMAP_PARAMETER_%02d", i + 1);
    npline = ScanToFlag(Ctop, code, npline, 1);
    tp->MParam[i].note    = (char*)malloc(MAXLINE * sizeof(char));
    tp->MParam[i].esrf    = CreateDmat(itmp[i], itmp[i], 0);
    tp->MParam[i].dphi    = CreateDmat(itmp[i], itmp[i], 0);
    tp->MParam[i].dpsi    = CreateDmat(itmp[i], itmp[i], 0);
    tp->MParam[i].d2pp    = CreateDmat(itmp[i], itmp[i], 0);
    tp->MParam[i].Ldphi   = CreateDmat(itmp[i], itmp[i], 0);
    tp->MParam[i].Ldpsi   = CreateDmat(itmp[i], itmp[i], 0);
    tp->MParam[i].Ld2pp   = CreateDmat(itmp[i], itmp[i], 0);
    tp->MParam[i].stencil = CreateDmat(2, itmp[i], 0);
    if (strncmp(Ctop->map[npline - 1], "%COMMENT", 8) == 0) {
      j = 8;
      while (Ctop->map[npline - 1][j] == ' ') {
        j++;
      }
      k = 0;
      while (Ctop->map[npline - 1][j] != '\n') {
        tp->MParam[i].note[k] = Ctop->map[npline - 1][j];
        k++;
        j++;
      }
      tp->MParam[i].note[k] = '\0';
      dataoffset = 0;
    }
    else {
      tp->MParam[i].note[0] = '\0';
      dataoffset = -1;
    }
    if (strncmp(Ctop->map[npline + dataoffset], "%FORMAT(8F9.5)", 14) == 0) {
      npline += dataoffset + 1;
      m = 0;
      for (j = 0; j < tp->MParam[i].esrf.row; j++) {
        for (k = 0; k < tp->MParam[i].esrf.col; k++) {	  
          tp->MParam[i].esrf.map[j][k] =
            RealXpYf(&Ctop->map[npline][9*m], 9, 5);
          m++;
          if (m == 8) {
            m = 0;
            npline++;
          }
        }
      }
    }
    else {
      printf("GetPrmTop >> Error.  Topology file %s\nGetPrmTop >> does not "
             "conform to CMAP format conventions.  Up to one line of \n"
             "GetPrmTop >> comments are permitted, but data must then follow "
             "in 8F9.5 format.\n", tp->source);
      printf("GetPrmTop >>\n");
      printf("GetPrmTop >> Line contents:\nGetPrmTop >> %s\n",
             Ctop->map[npline + dataoffset]);
      exit(1);
    }
  }
  free(code);
  free(itmp);
  npline = ScanToFlag(Ctop, "CMAP_INDEX", npline, 0);
  tp->Mterms = (cmapterm*)ReadXiL(Ctop, &npline, 6 * tp->ncmap, 6, 8);

  // Adjustment of the CMAP array is not the same
  // as adjustment of other bonded term arrays
  for (i = 0; i < tp->ncmap; i++) {
    tp->Mterms[i].a   -= 1;
    tp->Mterms[i].b   -= 1;
    tp->Mterms[i].c   -= 1;
    tp->Mterms[i].d   -= 1;
    tp->Mterms[i].e   -= 1;
    tp->Mterms[i].idx -= 1;
  }

  // Compute the CMAP derivatives
  for (i = 0; i < tp->ncmapsurf; i++) {
    ComputeCmapDerivatives(&tp->MParam[i]);
  }
}

//-----------------------------------------------------------------------------
// SetBasicTopology: set important parameters in a topology to common values.
//                   This routine is intended for specific purposes, and will
//                   allocate some arrays of its own.  This could generate
//                   modest memory leaks if the function is used for other
//                   purposes.
//
// Arguments:
//   tp:       blank topology struct (no allocations anticipated)
//   dsource:  name of the topology file that will be read
//   epsource: name of the extra points rule file to use
//   dorattle: flag to have this topology initialized with RATTLE on
//   dosettle: flag to have this topology initialized with RATTLE on
//-----------------------------------------------------------------------------
void SetBasicTopology(prmtop *tp, char* dsource, char* epsource, int dorattle,
                      int dosettle)
{
  tp->lj14fac = 0.5;
  tp->elec14fac = 1.0 - 1.0/1.2;
  tp->rattle = dorattle;
  tp->settle = dosettle;
  tp->ljbuck = 0;
  tp->nCustomVS = 0;
  sprintf(tp->WaterName, "WAT");
  tp->norattlemask = (char*)calloc(MAXNAME, sizeof(char));
  tp->rattlemask = (char*)calloc(MAXNAME, sizeof(char));
  tp->norattlemask[0] = '\0';
  tp->rattlemask[0] = '\0';
  strcpy(tp->source, dsource);
  strcpy(tp->eprulesource, epsource);
  tp->lVDWc = (double*)calloc(32, sizeof(double));  
}

//-----------------------------------------------------------------------------
// AdjustBondArray: the bonded arrays are tricky to adjust, because the atom
//                  numbers are given as pointer indices and the bond constants
//                  are in need of the typical Fortran-C array indexing
//                  adjustment.  We assume array A to have N sets of P numbers
//                  each, of which the first P-1 are atom numbers.
//-----------------------------------------------------------------------------
static void AdjustBondArray(int* A, int N, int P)
{
  int h, i, j;

  h = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < P-1; j++) {
      A[h] = A[h]/3;
      h++;
    }
    A[h] -= 1;
    h++;
  }
}

//-----------------------------------------------------------------------------
// ResetBondArray: undo the operation of AdjustBondArray.
//
// Arguments:
//   A:       integer array representation of the bond array
//   N:       the number of terms in the bond array (the length of A is N*P)
//   P:       the size of each term
//-----------------------------------------------------------------------------
static void ResetBondArray(int* A, int N, int P)
{
  int h, i, j;

  h = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < P-1; j++) {
      A[h] = A[h] * 3;
      h++;
    }
    A[h] += 1;
    h++;
  }
}

//-----------------------------------------------------------------------------
// GetPrmTop: load an AMBER 7 prmtop topology file and return it as a prmtop
//            data structure.                                    
//                                                                      
// Arguments:
//   tp:         the topology structure                                 
//   tj:         trajectory control information (contains the thread ID for
//               MPI mode, as well as belly mask or RATTLE information)
//   adjbnd:     flag to activate number index adjustments, to make the 
//               Fortran-designed AMBER prmtop readable by a C program  
//-----------------------------------------------------------------------------
void GetPrmTop(prmtop *tp, trajcon *tj, int adjbnd)
{
  int i, j, npline, varscee, varscnb, nvsfrm;
  cmat Ctop;
  
  // Read the topology into memory as an ascii character matrix
#ifdef MPI
  int tlen[2];
  if (tj->tid == 0) {
    Ctop = Ascii2Mem(tp->source, 96, 8, "Missing Topology file.");
    tlen[0] = Ctop.row;
    tlen[1] = Ctop.col;
  }
  MPI_Bcast(tlen, 2, MPI_INT, 0, MPI_COMM_WORLD);
  if (tj->tid != 0) {
    Ctop = CreateCmat(tlen[0], tlen[1]);
  }
  MPI_Valgrind_bcast(Ctop.data, tlen[0]*tlen[1], MPI_CHAR, 0, MPI_COMM_WORLD,
                     IDEBUG);
#else
  Ctop = Ascii2Mem(tp->source, 96, 8, "Missing Topology file.");
#endif

  // Version stamp
  if (Ctop.row == 0) {
    printf("GetPrmTop >> Error.  Topology file %s contains no data.\n",
           tp->source);
    exit(1);
  }
  strcpy(tp->vstamp, Ctop.map[0]);
  
  // First, the pointer information
  PrmTopPreamble(&Ctop, tp, &npline);

  // Now, read each atom name, charge, mass, type, type name, and class
  npline = ScanToFlag(&Ctop, "ATOM_NAME", npline, 1);
  tp->AtomNames = Read20a4(&Ctop, &npline, tp->natom);
  npline = ScanToFlag(&Ctop, "CHARGE", npline, 1);
  tp->Charges = Read5e16(&Ctop, &npline, tp->natom);
  tp->OrigCharges = CpyDVec(tp->Charges, tp->natom);
  DVecMult(tp->Charges, tp->natom, 1.0/18.2223);
  npline = ScanToFlag(&Ctop, "ATOMIC_NUMBER", npline, 0);
  if (npline > 0) {
    tp->ZNumber = ReadXiL(&Ctop, &npline, tp->natom, 10, 8);
  }
  else {
    tp->ZNumber = (int*)calloc(tp->natom, sizeof(int));
  }
  npline = ScanToFlag(&Ctop, "MASS", npline, 1);
  tp->Masses = Read5e16(&Ctop, &npline, tp->natom);
  npline = ScanToFlag(&Ctop, "ATOM_TYPE_INDEX", npline, 1);
  tp->LJIdx = ReadXiL(&Ctop, &npline, tp->natom, 10, 8);
  if (adjbnd == 1) {
    IVecAdd(tp->LJIdx, tp->natom, -1);
  }
  npline = ScanToFlag(&Ctop, "AMBER_ATOM_TYPE", npline, 1);
  tp->AtomTypes = Read20a4(&Ctop, &npline, tp->natom);
  npline = ScanToFlag(&Ctop, "TREE_CHAIN_CLASSIFICATION", npline, 1);
  tp->TreeSymbols = Read20a4(&Ctop, &npline, tp->natom);

  // Read the number of excluded atoms and make lists
  npline = ScanToFlag(&Ctop, "NUMBER_EXCLUDED_ATOMS", npline, 1);
  tp->NExcl = ReadXiL(&Ctop, &npline, tp->natom, 10, 8);
  tp->ConExcl = (int*)calloc(tp->natom+1, sizeof(int));
  j = 0;
  for (i = 0; i < tp->natom; i++) {
    tp->ConExcl[i] = j;
    j = j + tp->NExcl[i];
  }
  tp->ConExcl[tp->natom] = j;
  npline = ScanToFlag(&Ctop, "EXCLUDED_ATOMS_LIST", npline, 1);
  tp->ExclList = ReadXiL(&Ctop, &npline, tp->ConExcl[tp->natom], 10, 8);
  if (adjbnd == 1) {
    IVecAdd(tp->ExclList, tp->ConExcl[tp->natom], -1);
  }
  
  // Now get the force constants for all sorts of interactions
  npline = ScanToFlag(&Ctop, "NONBONDED_PARM_INDEX", npline, 1);
  tp->NBParmIdx = ReadXiL(&Ctop, &npline, tp->ntypes*tp->ntypes, 10, 8);
  if (adjbnd == 1) {
    IVecAdd(tp->NBParmIdx, tp->ntypes*tp->ntypes, -1);
  }
  npline = ScanToFlag(&Ctop, "BOND_FORCE_CONSTANT", npline, 1);
  tp->BondK = Read5e16(&Ctop, &npline, tp->nBAH.nbond);
  npline = ScanToFlag(&Ctop, "BOND_EQUIL_VALUE", npline, 1);
  tp->BondEq = Read5e16(&Ctop, &npline, tp->nBAH.nbond);
  npline = ScanToFlag(&Ctop, "BOND_STIFFNESS_PULL_ADJ", npline, 0);
  if (npline > 0) {
    tp->BondPullK = Read5e16(&Ctop, &npline, tp->nBAH.nbond);
  }
  else {
    tp->BondPullK = (double*)calloc(tp->nBAH.nbond, sizeof(double));
  }
  npline = ScanToFlag(&Ctop, "BOND_EQUIL_PULL_ADJ", npline, 0);
  if (npline > 0) {
    tp->BondPullEq = Read5e16(&Ctop, &npline, tp->nBAH.nbond);
  }
  else {
    tp->BondPullEq = (double*)calloc(tp->nBAH.nbond, sizeof(double));
  }
  npline = ScanToFlag(&Ctop, "BOND_STIFFNESS_PRESS_ADJ", npline, 0);
  if (npline > 0) {
    tp->BondPressK = Read5e16(&Ctop, &npline, tp->nBAH.nbond);
  }
  else {
    tp->BondPressK = (double*)calloc(tp->nBAH.nbond, sizeof(double));
  }
  npline = ScanToFlag(&Ctop, "BOND_EQUIL_PRESS_ADJ", npline, 0);
  if (npline > 0) {
    tp->BondPressEq = Read5e16(&Ctop, &npline, tp->nBAH.nbond);
  }
  else {
    tp->BondPressEq = (double*)calloc(tp->nBAH.nbond, sizeof(double));
  }
  npline = ScanToFlag(&Ctop, "ANGLE_FORCE_CONSTANT", npline, 1);
  tp->AnglK = Read5e16(&Ctop, &npline, tp->nBAH.nangl);
  npline = ScanToFlag(&Ctop, "ANGLE_EQUIL_VALUE", npline, 1);
  tp->AnglEq = Read5e16(&Ctop, &npline, tp->nBAH.nangl);
  npline = ScanToFlag(&Ctop, "DIHEDRAL_FORCE_CONSTANT", npline, 1);
  tp->DiheK = Read5e16(&Ctop, &npline, tp->nBAH.ndihe);
  npline = ScanToFlag(&Ctop, "DIHEDRAL_PERIODICITY", npline, 1);
  tp->DiheN = Read5e16(&Ctop, &npline, tp->nBAH.ndihe);
  npline = ScanToFlag(&Ctop, "DIHEDRAL_PHASE", npline, 1);
  tp->DihePhi = Read5e16(&Ctop, &npline, tp->nBAH.ndihe);

  // These factors for variable 1-4 scaling may or may not
  // be present in the topology.  If they are not present,
  // the program will not abort.                          
  npline = ScanToFlag(&Ctop, "SCEE_SCALE_FACTOR", npline, 0);
  if (npline > 0) {
    tp->scee = Read5e16(&Ctop, &npline, tp->nBAH.ndihe);
    varscee = 1;
  }
  else {
    varscee = 0;
  }
  npline = ScanToFlag(&Ctop, "SCNB_SCALE_FACTOR", npline, 0);
  if (npline > 0) {
    tp->scnb = Read5e16(&Ctop, &npline, tp->nBAH.ndihe);
    varscnb = 1;
  }
  else {
    varscnb = 0;
  }

  // The SOLTY field is currently unused, but reserved.  Handle it anyway.
  npline = ScanToFlag(&Ctop, "SOLTY", npline, 0);
  if (npline > 0) {
    tp->solty = Read5e16(&Ctop, &npline, tp->natyp);
  }
  else {
    tp->solty = (double*)calloc(tp->natyp, sizeof(double));
  }
  
  // More force constants and nonbonded parameters
  npline = ScanToFlag(&Ctop, "LENNARD_JONES_ACOEF", npline, 1);
  tp->LJA = Read5e16(&Ctop, &npline, tp->ntypes*(tp->ntypes+1)/2);
  npline = ScanToFlag(&Ctop, "LENNARD_JONES_BCOEF", npline, 1);
  tp->LJB = Read5e16(&Ctop, &npline, tp->ntypes*(tp->ntypes+1)/2);
  if (tp->ljbuck == 1) {
    npline = ScanToFlag(&Ctop, "LENNARD_JONES_CCOEF", npline, 1);
    tp->LJC = Read5e16(&Ctop, &npline, tp->ntypes*(tp->ntypes+1)/2);
  }
  npline = ScanToFlag(&Ctop, "HBOND_ACOEF", npline, 1);
  tp->SolA = Read5e16(&Ctop, &npline, tp->nphb*(tp->nphb+1)/2);
  npline = ScanToFlag(&Ctop, "HBOND_BCOEF", npline, 1);
  tp->SolB = Read5e16(&Ctop, &npline, tp->nphb*(tp->nphb+1)/2);
  npline = ScanToFlag(&Ctop, "HBCUT", npline, 1);
  tp->HBCut = Read5e16(&Ctop, &npline, tp->nphb*(tp->nphb+1)/2);

  // The GB radii set description can be multiple words
  npline = ScanToFlag(&Ctop, "RADIUS_SET", npline, 0);
  if (npline >= 0) {
    i = 0;
    j = 0;
    while (i < Ctop.col && Ctop.map[npline][i] != '\n' &&
	   Ctop.map[npline][i] != '\0') {
      if (Ctop.map[npline][i] != ' ') {
	j = i+1;
      }
      i++;
    }
    for (i = 0; i < j; i++) {
      tp->RadSet[i] = Ctop.map[npline][i];
    }
    tp->RadSet[i] = '\0';
  }
  else {
    tp->RadSet[0] = '\0';
  }
  npline = ScanToFlag(&Ctop, "RADII", npline, 0);
  if (npline >= 0) {
    tp->Radii = Read5e16(&Ctop, &npline, tp->natom);
  }
  else {
    tp->Radii = (double*)calloc(tp->natom, sizeof(double));
  }
  npline = ScanToFlag(&Ctop, "SCREEN", npline, 0);
  if (npline >= 0) {
    tp->Screen = Read5e16(&Ctop, &npline, tp->natom);
  }
  else {
    tp->Screen = (double*)calloc(tp->natom, sizeof(double));
  }
 
  // Get residue limits and labels
  npline = ScanToFlag(&Ctop, "RESIDUE_LABEL", npline, 1);
  tp->ResNames = Read20a4(&Ctop, &npline, tp->nres);
  npline = ScanToFlag(&Ctop, "RESIDUE_POINTER", npline, 1);
  tp->ResLims = ReadXiL(&Ctop, &npline, tp->nres, 10, 8);
  if (adjbnd == 1) {
    IVecAdd(tp->ResLims, tp->nres, -1);
  }
  tp->ResLims = (int*)realloc(tp->ResLims, (tp->nres+1)*sizeof(int));
  tp->ResLims[tp->nres] = tp->natom;

  // Bond, angle, and dihedral lists
  npline = ScanToFlag(&Ctop, "BONDS_INC_HYDROGEN", npline, 1);
  tp->BIncH = (bond*)ReadXiL(&Ctop, &npline, 3*tp->withH.nbond, 10, 8);
  if (adjbnd == 1) {
    AdjustBondArray((int*)tp->BIncH, tp->withH.nbond, 3);
  }
  npline = ScanToFlag(&Ctop, "BONDS_WITHOUT_HYDROGEN", npline, 1);
  tp->BNoH = (bond*)ReadXiL(&Ctop, &npline, 3*tp->woH.nbond, 10, 8);
  if (adjbnd == 1) {
    AdjustBondArray((int*)tp->BNoH, tp->woH.nbond, 3);
  }
  npline = ScanToFlag(&Ctop, "ANGLES_INC_HYDROGEN", npline, 1);
  tp->AIncH = (angle*)ReadXiL(&Ctop, &npline, 4*tp->withH.nangl, 10, 8);
  if (adjbnd == 1) {
    AdjustBondArray((int*)tp->AIncH, tp->withH.nangl, 4);
  }
  npline = ScanToFlag(&Ctop, "ANGLES_WITHOUT_HYDROGEN", npline, 1);
  tp->ANoH = (angle*)ReadXiL(&Ctop, &npline, 4*tp->woH.nangl, 10, 8);
  if (adjbnd == 1) {
    AdjustBondArray((int*)tp->ANoH, tp->woH.nangl, 4);
  }
  npline = ScanToFlag(&Ctop, "DIHEDRALS_INC_HYDROGEN", npline, 1);
  tp->HIncH = (dihedral*)ReadXiL(&Ctop, &npline, 5*tp->withH.ndihe, 10, 8);
  if (adjbnd == 1) {
    AdjustBondArray((int*)tp->HIncH, tp->withH.ndihe, 5);
  }
  npline = ScanToFlag(&Ctop, "DIHEDRALS_WITHOUT_HYDROGEN", npline, 1);
  tp->HNoH = (dihedral*)ReadXiL(&Ctop, &npline, 5*tp->woH.ndihe, 10, 8);
  if (adjbnd == 1) {
    AdjustBondArray((int*)tp->HNoH, tp->woH.ndihe, 5);
  }

  // CMAP information
  tp->ncmap = 0;
  tp->ncmapsurf = 0;
  npline = ScanToFlag(&Ctop, "CMAP_COUNT", npline, 0);
  if (npline > 0) {
    ReadCmapInfo(&Ctop, npline, tp);
  }
  
  // Miscellaneous
  npline = ScanToFlag(&Ctop, "JOIN_ARRAY", npline, 1);
  tp->Join = ReadXiL(&Ctop, &npline, tp->natom, 10, 8);
  npline = ScanToFlag(&Ctop, "IROTAT", npline, 1);
  tp->Rotat = ReadXiL(&Ctop, &npline, tp->natom, 10, 8);

  // For periodic boxes
  if (tp->ifbox > 0) {
    npline = ScanToFlag(&Ctop, "SOLVENT_POINTERS", npline, 1);
    sscanf(Ctop.map[npline], "%d%d%d", &tp->iptres, &tp->nspm, &tp->nspsol);
    npline = ScanToFlag(&Ctop, "ATOMS_PER_MOLECULE", npline, 1);
    tp->Nsp = ReadXiL(&Ctop, &npline, tp->nspm, 10, 8);
    npline = ScanToFlag(&Ctop, "BOX_DIMENSIONS", npline, 1);
    sscanf(Ctop.map[npline], "%lf%lf%lf%lf", &tp->gdim[4], &tp->gdim[0],
           &tp->gdim[1], &tp->gdim[2]);
  }

  // For customized virtual sites
  if (tp->numextra > 0) {
    npline = ScanToFlag(&Ctop, "VIRTUAL_SITE_FRAMES", npline, 0);
    if (npline >= 0) {
      tp->customVSi = Read10i8Test(&Ctop, &npline, &nvsfrm);

      // Check the sanity of custom frame specifications
      if (nvsfrm % 6 != 0) {
	printf("GetPrmTop >> Warning.  A virtual site frame definition "
               "appears to be incomplete.\n");
      }
      tp->nCustomVS = nvsfrm / 6;
      npline = ScanToFlag(&Ctop, "VIRTUAL_SITE_FRAME_DETAILS", npline, 0);
      tp->customVSd = Read5e16(&Ctop, &npline, tp->nCustomVS * 3);      
    }
    else {
      tp->nCustomVS = 0;
    }
  }
  else {
    tp->nCustomVS = 0;
  }
  
  // For polarizable topologies (just read IPOL for now)
  npline = ScanToFlag(&Ctop, "IPOL", npline, 0);
  if (npline >= 0) {
    tp->ipol = atoi(Ctop.map[npline]);
  }
  else {
    tp->ipol = 0;
  }
    
  // Delete the memory copy of the input file
  DestroyCmat(&Ctop);

  // Allocate arrays for the 1-4 scaling factors if they do not  
  // already exist, and perform the same internal conversions as 
  // were done for the global 1-4 scaling factors in the         
  // GetCntrlNamelist function of the Command library if they do.
  if (varscee == 0) {
    tp->scee = (double*)malloc(tp->nBAH.ndihe*sizeof(double));
    SetDVec(tp->scee, tp->nBAH.ndihe, tp->elec14fac);
  }
  else {
    for (i = 0; i < tp->nBAH.ndihe; i++) {
      if (fabs(tp->scee[i]) > 1.0e-8) {
        tp->scee[i] = 1.0 - 1.0/tp->scee[i];
      }
      else {
        tp->scee[i] = 0.0;
      }
    }
  }
  if (varscnb == 0) {
    tp->scnb = (double*)malloc(tp->nBAH.ndihe*sizeof(double));
    SetDVec(tp->scnb, tp->nBAH.ndihe, tp->lj14fac);
  }
  else {
    for (i = 0; i < tp->nBAH.ndihe; i++) {
      if (fabs(tp->scnb[i]) > 1.0e-8) {
        tp->scnb[i] = 1.0 - 1.0/tp->scnb[i];
      }
      else {
        tp->scnb[i] = 0.0;
      }
    }
  }

  // Pre-computations for masses or other constants on a per-atom basis
  tp->InvMasses = (double*)malloc(tp->natom*sizeof(double));
  tp->TotalMass = 0.0;
  for (i = 0; i < tp->natom; i++) {
    tp->InvMasses[i] = (tp->Masses[i] < 1.0e-8) ? 0.0 : 1.0/tp->Masses[i];
    tp->TotalMass += tp->Masses[i];
  }

  // Currently no support for IFCAP > 0, IFPERT > 0, and/or IPOL == 1

  // Re-organize some of the convoluted topology information
  OrderLJParameters(tp);
  OrderBondParameters(tp);
  EnumerateBonds(tp);
  EnumerateAngles(tp);
  EnumerateDihedrals(tp);
  EnumerateCmaps(tp);

  // Define groups of linked atoms
  Connect1234(tp);
  DefineAtomChains(tp);

  // Prepare for SETTLE; the RattleGrpMax field is set to a minimum
  // of 9 here because SETTLE constraints make use of some of the  
  // same arrays later allocated in the CellPositionCnst and       
  // CnstPositionRscl functions (see ConstrainPos.c)               
  tp->RattleGrpMax = 9;
  FindThreePointWaters(tp);
  if (tp->nwat > 0) {
    CompSettleGeom(tp);
  }

  // Prepare for RATTLE
  FindRattleBonds(tp);

  // Process information regarding extra points
#ifdef MPI
  if (tj->tid == 0) {
    ReadEPRuleFile(tp);
  }
  BroadcastEPInfo(tp, tj);
#else
  ReadEPRuleFile(tp);
#endif
  DetermineEPFrames(tp);
  FinalizeNewEP(tp);

  // If there are any exclusions as yet
  // unaccounted for, fulfill them now.
  FulfillExclusions(tp);

  // Round the total amount of charge on the system
  RoundTotalCharge(tp);

  // Set a bellymask
  CreateBellyMask(tj, tp);

  // Determine the number of system degrees of freedom and particles
  tp->ncnst = 0;
  tp->nprtcl = tp->natom - tp->nxtrapt;
  if (tp->rattle == 1 || tp->settle == 1) {
    for (i = 0; i < tp->natom; i++) {
      if (tp->SHL[i].exe == 1) {
        tp->ncnst += 3;
        tp->nprtcl -= 2;
      }
      if (tp->SHL[i].exe == 2) {
        tp->ncnst += tp->SHL[i].blist[0];
        tp->nprtcl -= tp->SHL[i].blist[3*tp->SHL[i].blist[0]+1] - 1;
      }
    }
  }
  j = 0;
  for (i = 0; i < tp->natom; i++) {
    if (tp->MobileAtoms[i] == 1) {
      j++;
    }
  }
  tp->ndf = 3*(j - 1) - tp->ncnst;
}

//-----------------------------------------------------------------------------
// CopyTopology: copy a topology structure.                             
//                                                                      
// Arguments:                                                           
//   tp:   the original topology                                        
//-----------------------------------------------------------------------------
prmtop CopyTopology(prmtop *tp)
{
  int i, j, natm, nnix, grpsize, nuniq;
  prmtop tpC;

  // General copying
  tpC = *tp;

  // Copy integer arrays
  natm = tp->natom; 
  tpC.LJIdx = CpyIVec(tp->LJIdx, natm);
  tpC.NExcl = CpyIVec(tp->NExcl, natm);
  tpC.ConExcl = CpyIVec(tp->ConExcl, natm+1);
  tpC.ExclList = CpyIVec(tp->ExclList, tp->ConExcl[natm]);
  tpC.NBParmIdx = CpyIVec(tp->NBParmIdx, tp->ntypes*tp->ntypes);
  tpC.ResLims = CpyIVec(tp->ResLims, tp->nres + 1);
  tpC.Join = CpyIVec(tp->Join, natm);
  tpC.Rotat = CpyIVec(tp->Rotat, natm);
  if (tp->ifbox > 0) {
    tpC.Nsp = CpyIVec(tp->Nsp, tp->nspm);
  }
  
  // Copy double-precision real arrays
  tpC.Charges = CpyDVec(tp->Charges, natm);
  tpC.OrigCharges = CpyDVec(tp->OrigCharges, natm);
  tpC.ZNumber = CpyIVec(tp->ZNumber, natm);
  tpC.Masses = CpyDVec(tp->Masses, natm);
  tpC.InvMasses = CpyDVec(tp->InvMasses, natm);
  tpC.BondK = CpyDVec(tp->BondK, tp->nBAH.nbond);
  tpC.BondEq = CpyDVec(tp->BondEq, tp->nBAH.nbond);
  tpC.BondPullK = CpyDVec(tp->BondPullK, tp->nBAH.nbond);
  tpC.BondPullEq = CpyDVec(tp->BondPullEq, tp->nBAH.nbond);
  tpC.BondPressK = CpyDVec(tp->BondPressK, tp->nBAH.nbond);
  tpC.BondPressEq = CpyDVec(tp->BondPressEq, tp->nBAH.nbond);
  tpC.AnglK = CpyDVec(tp->AnglK, tp->nBAH.nangl);
  tpC.AnglEq = CpyDVec(tp->AnglEq, tp->nBAH.nangl);
  tpC.DiheK = CpyDVec(tp->DiheK, tp->nBAH.ndihe);
  tpC.DiheN = CpyDVec(tp->DiheN, tp->nBAH.ndihe);
  tpC.DihePhi = CpyDVec(tp->DihePhi, tp->nBAH.ndihe);
  tpC.LJA = CpyDVec(tp->LJA, tp->ntypes*(tp->ntypes+1)/2);
  tpC.LJB = CpyDVec(tp->LJB, tp->ntypes*(tp->ntypes+1)/2);
  tpC.SolA = CpyDVec(tp->SolA, tp->nphb*(tp->nphb+1)/2);
  tpC.SolB = CpyDVec(tp->SolB, tp->nphb*(tp->nphb+1)/2);
  tpC.HBCut = CpyDVec(tp->HBCut, tp->nphb*(tp->nphb+1)/2);
  tpC.Radii = CpyDVec(tp->Radii, natm);
  tpC.Screen = CpyDVec(tp->Screen, natm);
  tpC.scee = CpyDVec(tp->scee, tp->nBAH.ndihe);
  tpC.scnb = CpyDVec(tp->scnb, tp->nBAH.ndihe);
  tpC.solty = CpyDVec(tp->solty, tp->natyp);
  tpC.lVDWc = CpyDVec(tp->lVDWc, tp->ntypes);
  tpC.MobileAtoms = CpyIVec(tp->MobileAtoms, tp->natom);

  // Copy character arrays
  tpC.AtomNames = (char*)malloc(4 * natm * sizeof(char));
  tpC.AtomTypes = (char*)malloc(4 * natm * sizeof(char));
  tpC.TreeSymbols = (char*)malloc(4 * natm * sizeof(char));
  tpC.ResNames = (char*)malloc(4 * tp->nres * sizeof(char));
  for (i = 0; i < 4*natm; i++) {
    tpC.AtomNames[i] = tp->AtomNames[i];
    tpC.AtomTypes[i] = tp->AtomTypes[i];
    tpC.TreeSymbols[i] = tp->TreeSymbols[i];
  }
  for (i = 0; i < 4*tp->nres; i++) {
    tpC.ResNames[i] = tp->ResNames[i];
  }
  tpC.norattlemask = (char*)calloc(MAXNAME, sizeof(char));
  tpC.rattlemask = (char*)calloc(MAXNAME, sizeof(char));
  strcpy(tpC.norattlemask, tp->norattlemask);
  strcpy(tpC.rattlemask, tp->rattlemask);

  // Copy arrays of structures
  tpC.BIncH = (bond*)malloc(tp->withH.nbond * sizeof(bond));
  for (i = 0; i < tp->withH.nbond; i++) {
    tpC.BIncH[i] = tp->BIncH[i];
  }
  tpC.BNoH = (bond*)malloc(tp->woH.nbond * sizeof(bond));
  for (i = 0; i < tp->woH.nbond; i++) {
    tpC.BNoH[i] = tp->BNoH[i];
  }
  tpC.AIncH = (angle*)malloc(tp->withH.nangl * sizeof(angle));
  for (i = 0; i < tp->withH.nangl; i++) {
    tpC.AIncH[i] = tp->AIncH[i];
  }
  tpC.ANoH = (angle*)malloc(tp->woH.nangl * sizeof(angle));
  for (i = 0; i < tp->woH.nangl; i++) {
    tpC.ANoH[i] = tp->ANoH[i];
  }
  tpC.HIncH = (dihedral*)malloc(tp->withH.ndihe * sizeof(dihedral));
  for (i = 0; i < tp->withH.ndihe; i++) {
    tpC.HIncH[i] = tp->HIncH[i];
  }
  tpC.HNoH = (dihedral*)malloc(tp->woH.ndihe * sizeof(dihedral));
  for (i = 0; i < tp->woH.ndihe; i++) {
    tpC.HNoH[i] = tp->HNoH[i];
  }
  tpC.BLC = (bondlist*)malloc(natm * sizeof(bondlist));
  for (i = 0; i < tp->natom; i++) {
    tpC.BLC[i] = tp->BLC[i];
    tpC.BLC[i].BC = (bondcomm*)malloc(tpC.BLC[i].nbond * sizeof(bondcomm));
    for (j = 0; j < tpC.BLC[i].nbond; j++) {
      tpC.BLC[i].BC[j] = tp->BLC[i].BC[j];
    }
  }
  tpC.ALC = (angllist*)malloc(natm * sizeof(angllist));
  for (i = 0; i < tp->natom; i++) {
    tpC.ALC[i] = tp->ALC[i];
    tpC.ALC[i].AC = (anglcomm*)malloc(tpC.ALC[i].nangl * sizeof(anglcomm));
    for (j = 0; j < tpC.ALC[i].nangl; j++) {
      tpC.ALC[i].AC[j] = tp->ALC[i].AC[j];
    }
  }
  tpC.HLC = (dihelist*)malloc(natm * sizeof(dihelist));
  for (i = 0; i < tp->natom; i++) {
    tpC.HLC[i] = tp->HLC[i];
    tpC.HLC[i].HC = (dihecomm*)malloc(tpC.HLC[i].ndihe * sizeof(dihecomm));
    for (j = 0; j < tpC.HLC[i].ndihe; j++) {
      tpC.HLC[i].HC[j] = tp->HLC[i].HC[j];
      tpC.HLC[i].HC[j].t = CpyIVec(tp->HLC[i].HC[j].t, tp->HLC[i].HC[j].nt);
    }
  }
  tpC.MLC = (cmaplist*)malloc(natm * sizeof(cmaplist));
  for (i = 0; i < tp->natom; i++) {
    tpC.MLC[i] = tp->MLC[i];
    tpC.MLC[i].MC = (cmapcomm*)malloc(tpC.MLC[i].ncmap * sizeof(cmapcomm));
    for (j = 0; j < tpC.MLC[i].ncmap; j++) {
      tpC.MLC[i].MC[j] = tp->MLC[i].MC[j];      
    }
  }
  tpC.BParam = (bonddef*)malloc(tp->nBAH.nbond * sizeof(bonddef));
  for (i = 0; i < tp->nBAH.nbond; i++) {
    tpC.BParam[i] = tp->BParam[i];
  }
  tpC.AParam = (angldef*)malloc(tp->nBAH.nangl * sizeof(angldef));
  for (i = 0; i < tp->nBAH.nangl; i++) {
    tpC.AParam[i] = tp->AParam[i];
  }
  tpC.HParam = (dihedef*)malloc(tp->nBAH.ndihe * sizeof(dihedef));
  for (i = 0; i < tp->nBAH.ndihe; i++) {
    tpC.HParam[i] = tp->HParam[i];
  }
  tpC.Mterms = (cmapterm*)malloc(tp->ncmap * sizeof(cmapterm));
  for (i = 0; i < tp->ncmap; i++) {
    tpC.Mterms[i] = tp->Mterms[i];
  }
  tpC.MParam = (cmapdef*)malloc(tp->ncmapsurf * sizeof(cmapdef));
  for (i = 0; i < tp->ncmapsurf; i++) {

    // The copying happens field by field, because I'm not sure what
    // happens if CopyDmat is applied and asked to create new data
    // when its two dmat pointers point to the same thing.
    tpC.MParam[i].note = (char*)malloc(MAXLINE * sizeof(char));
    strcpy(tpC.MParam[i].note, tp->MParam[i].note);
    CopyDmat(&tpC.MParam[i].esrf, &tp->MParam[i].esrf, 0);
    CopyDmat(&tpC.MParam[i].dphi, &tp->MParam[i].dphi, 0);
    CopyDmat(&tpC.MParam[i].dpsi, &tp->MParam[i].dpsi, 0);
    CopyDmat(&tpC.MParam[i].d2pp, &tp->MParam[i].d2pp, 0);
    CopyDmat(&tpC.MParam[i].Ldphi, &tp->MParam[i].Ldphi, 0);
    CopyDmat(&tpC.MParam[i].Ldpsi, &tp->MParam[i].Ldpsi, 0);
    CopyDmat(&tpC.MParam[i].Ld2pp, &tp->MParam[i].Ld2pp, 0);
    CopyDmat(&tpC.MParam[i].stencil, &tp->MParam[i].stencil, 0);
    tpC.MParam[i].muePhi = tp->MParam[i].muePhi;
    tpC.MParam[i].muePsi = tp->MParam[i].muePsi;
    tpC.MParam[i].mue2PP = tp->MParam[i].mue2PP;
    tpC.MParam[i].mueFPhi = tp->MParam[i].mueFPhi;
    tpC.MParam[i].mueFPsi = tp->MParam[i].mueFPsi;
    tpC.MParam[i].mueF2PP = tp->MParam[i].mueF2PP;
  }
  tpC.SHL = (cnstcomm*)malloc(natm * sizeof(cnstcomm));
  tpC.FR1Idx = (lgrp*)malloc(natm * sizeof(lgrp));
  for (i = 0; i < natm; i++) {

    // Copy constraint groups
    tpC.SHL[i] = tp->SHL[i];
    if (tpC.SHL[i].exe == 1) {
      tpC.SHL[i].blist = CpyIVec(tp->SHL[i].blist, 2);
    }
    if (tpC.SHL[i].exe == 2) {
      grpsize = tp->SHL[i].blist[0];
      nuniq = tp->SHL[i].blist[3*grpsize + 1];
      tpC.SHL[i].blist = CpyIVec(tp->SHL[i].blist, 3*grpsize + 2 + nuniq);
    }

    // Copy frame atom 1 data
    tpC.FR1Idx[i].natom = tp->FR1Idx[i].natom;
    tpC.FR1Idx[i].atoms = CpyIVec(tp->FR1Idx[i].atoms, tp->FR1Idx[i].natom);
  }

  // Copy matrices
  tpC.LJftab = CreateDmat(tp->ntypes, 2*tp->ntypes, 0);
  tpC.LJutab = CreateDmat(tp->ntypes, 2*tp->ntypes, 0);
  for (i = 0; i < 2*tp->ntypes*tp->ntypes; i++) {
    tpC.LJftab.data[i] = tp->LJftab.data[i];
    tpC.LJutab.data[i] = tp->LJutab.data[i];
  }

  // Copy the connectivity array
  tpC.nb1234 = (map1234*)malloc(natm * sizeof(map1234));
  for (i = 0; i < natm; i++) {
    tpC.nb1234[i] = tp->nb1234[i];
    tpC.nb1234[i].L11 = CpyIVec(tp->nb1234[i].L11, tp->nb1234[i].n11);
    tpC.nb1234[i].L12 = CpyIVec(tp->nb1234[i].L12, tp->nb1234[i].n12);
    tpC.nb1234[i].L13 = CpyIVec(tp->nb1234[i].L13, tp->nb1234[i].n13);
    tpC.nb1234[i].L14 = CpyIVec(tp->nb1234[i].L14, tp->nb1234[i].n14);
  }

  // Copy group lists
  tpC.ngrp = tp->ngrp;
  tpC.lgrps = (lgrp*)malloc(tpC.ngrp * sizeof(lgrp));
  for (i = 0; i < tpC.ngrp; i++) {
    tpC.lgrps[i].natom = tp->lgrps[i].natom;
    tpC.lgrps[i].atoms = CpyIVec(tp->lgrps[i].atoms, tpC.lgrps[i].natom);
  }

  // Copy data related to new extra-points
  if (tpC.nxtrapt > 0) {
    tpC.xtrapts = (expt*)malloc(tpC.nxtrapt * sizeof(expt));
  }
  else {
    tpC.xtrapts = (expt*)malloc(sizeof(expt));
  }
  for (i = 0; i < tpC.nxtrapt; i++) {
    tpC.xtrapts[i] = tp->xtrapts[i];
  }
  if (tp->EPInserted == 1 || tp->ExclMarked == 1) {
    tpC.OldAtomNum = CpyIVec(tp->OldAtomNum, natm);
    tpC.ElimPair = (auxelim*)malloc(tpC.natom * sizeof(auxelim));
    for (i = 0; i < tpC.natom; i++) {
      tpC.ElimPair[i] = tp->ElimPair[i];
      nnix = MAX(tpC.ElimPair[i].n11, 1);
      tpC.ElimPair[i].list11 = (nixpr*)malloc(nnix * sizeof(nixpr));
      for (j = 0; j < nnix; j++) {
        tpC.ElimPair[i].list11[j] = tp->ElimPair[i].list11[j];
      }
      nnix = MAX(tpC.ElimPair[i].n12, 1);
      tpC.ElimPair[i].list12 = (nixpr*)malloc(nnix * sizeof(nixpr));
      for (j = 0; j < nnix; j++) {
        tpC.ElimPair[i].list12[j] = tp->ElimPair[i].list12[j];
      }
      nnix = MAX(tpC.ElimPair[i].n13, 1);
      tpC.ElimPair[i].list13 = (nixpr*)malloc(nnix * sizeof(nixpr));
      for (j = 0; j < nnix; j++) {
        tpC.ElimPair[i].list13[j] = tp->ElimPair[i].list13[j];
      }
      nnix = MAX(tpC.ElimPair[i].n14, 1);
      tpC.ElimPair[i].list14 = (nixpr*)malloc(nnix * sizeof(nixpr));
      for (j = 0; j < nnix; j++) {
        tpC.ElimPair[i].list14[j] = tp->ElimPair[i].list14[j];
      }
    }
  }
  if (tp->nCustomVS > 0) {
    tpC.customVSi = CpyIVec(tp->customVSi, tp->nCustomVS * 6);
    tpC.customVSd = CpyDVec(tp->customVSd, tp->nCustomVS * 3);
  }

  return tpC;
}

//-----------------------------------------------------------------------------
// InterpolateTopology: interpolate a new topology by linear combination of
//                      two other topologies.                        
//                                                                      
// Arguments:                                                           
//   tp[A,B]:    the topology endpoints for interpolation               
//   lambda:     the interpolation coordinate                           
//-----------------------------------------------------------------------------
prmtop InterpolateTopology(prmtop *tpA, prmtop *tpB, double lambda)
{
  int i, iprob;
  double amix, bmix;
  prmtop tpN;

  // Mixing factors for each topology
  amix = 1.0 - lambda;
  bmix = lambda;
  
  // Check for congruence between the topologies
  iprob = 0;
  if (tpA->natom != tpB->natom) {
    printf("InterpolateTopology >> Atom count mismatch (%d vs %d)\n",
           tpA->natom, tpB->natom);
    iprob = 1;
  }
  if (tpA->nres != tpB->nres) {
    printf("InterpolateTopology >> Residue count mismatch (%d vs %d)\n",
           tpA->nres, tpB->nres);
    iprob = 1;
  }
  if (tpA->ntypes != tpB->ntypes) {
    printf("InterpolateTopology >> Atom type count mismatch (%d vs %d)\n",
           tpA->ntypes, tpB->ntypes);
    iprob = 1;
  }
  for (i = 0; i < (1 - iprob)*tpA->natom; i++) {
    if (tpA->LJIdx[i] != tpB->LJIdx[i]) {
      printf("InterpolateTopology >> Mismatch in Lennard-Jones indexing (atom "
             "%d)\n", i);
      iprob = 1;
    }
    if (str4cmp(&tpA->AtomNames[4*i], &tpB->AtomNames[4*i]) != 0) {
      printf("InterpolateTopology >> Atom name mismatch (atom %d)\n", i + 1);
      iprob = 1;
    }
  }
  for (i = 0; i < (1 - iprob)*tpA->nres; i++) {
    if (str4cmp(&tpA->ResNames[4*i], &tpB->ResNames[4*i]) != 0) {
      printf("InterpolateTopology >> Residue name mismatch (residue %d)\n",
             i + 1);
      iprob = 1;
    }
  }
  if (tpA->nBAH.nbond != tpB->nBAH.nbond) {
    printf("InterpolateTopology >> Unique bond count mismatch (%d vs %d)\n",
           tpA->nBAH.nbond, tpB->nBAH.nbond);
    iprob = 1;
  }
  if (tpA->nBAH.nangl != tpB->nBAH.nangl) {
    printf("InterpolateTopology >> Unique angle count mismatch (%d vs %d)\n",
           tpA->nBAH.nangl, tpB->nBAH.nangl);
    iprob = 1;
  }
  if (tpA->nBAH.ndihe != tpB->nBAH.ndihe) {
    printf("InterpolateTopology >> Unique dihedral count mismatch "
           "(%d vs %d)\n", tpA->nBAH.ndihe, tpB->nBAH.ndihe);
    iprob = 1;
  }
  for (i = 0; i < tpA->nBAH.ndihe; i++) {
    if ((int)(tpA->HParam[i].N + 1.0e-8) != (int)(tpB->HParam[i].N + 1.0e-8)) {
      printf("InterpolateTopology >> Dihedral periodicity mismatch "
             "(%d vs %d)\n", (int)(tpA->HParam[i].N + 1.0e-8),
             (int)(tpB->HParam[i].N + 1.0e-8));
      iprob = 1;
    }
  }
  if (tpA->ncmap != tpB->ncmap) {
    printf("InterpolateTopology >> Unique CMAP count mismatch (%d vs %d)\n",
	   tpA->ncmap, tpB->ncmap);
    iprob = 1;
  }
  for (i = 0; i < tpA->ncmapsurf; i++) {
    if (tpA->MParam[i].esrf.row != tpB->MParam[i].esrf.row) {
      printf("InterpolateTopology >> CMAP resolution mismatch (%d vs %d)\n",
	     tpA->MParam[i].esrf.row, tpB->MParam[i].esrf.row);
      iprob = 1;
    }
  }
  if (tpA->ipol != tpB->ipol) {
    iprob = 1;
  }
  if (iprob == 1) {
    printf("InterpolateTopology >>   Topology 1: %s\n", tpA->source);
    printf("InterpolateTopology >>   Topology 2: %s\n", tpB->source);
    exit(1);
  }

  // Start making a new topology
  tpN = CopyTopology(tpA);
  for (i = 0; i < tpA->natom; i++) {
    tpN.Charges[i] = (amix * tpA->Charges[i]) +
                     (bmix * tpB->Charges[i]);
    tpN.OrigCharges[i] = (amix * tpA->OrigCharges[i]) +
                         (bmix * tpB->OrigCharges[i]);
  }
  for (i = 0; i < tpA->ntypes * (tpA->ntypes + 1) / 2; i++) {
    tpN.LJA[i] = (amix * tpA->LJA[i]) + (bmix * tpB->LJA[i]);
    tpN.LJB[i] = (amix * tpA->LJB[i]) + (bmix * tpB->LJB[i]);
  }
  CopyDmat(&tpN.LJftab, &tpA->LJftab, 0);
  CopyDmat(&tpN.LJutab, &tpA->LJutab, 0);
  AddDmat(&tpN.LJftab, &tpB->LJftab, 1.0 - lambda, lambda);
  AddDmat(&tpN.LJutab, &tpB->LJutab, 1.0 - lambda, lambda);
  for (i = 0; i < tpA->nBAH.nbond; i++) {
    tpN.BParam[i].K =  (amix * tpA->BParam[i].K)  + (bmix * tpB->BParam[i].K);
    tpN.BParam[i].l0 = (amix * tpA->BParam[i].l0) + (bmix * tpB->BParam[i].l0);
    tpN.BParam[i].Kpull = (amix * tpA->BParam[i].Kpull) +
                          (bmix * tpB->BParam[i].Kpull);
    tpN.BParam[i].lpull0 = (amix * tpA->BParam[i].lpull0) +
                           (bmix * tpB->BParam[i].lpull0);
    tpN.BParam[i].Kpress = (amix * tpA->BParam[i].Kpress) +
                          (bmix * tpB->BParam[i].Kpress);
    tpN.BParam[i].lpress0 = (amix * tpA->BParam[i].lpress0) +
                           (bmix * tpB->BParam[i].lpress0);
  }
  for (i = 0; i < tpA->nBAH.nangl; i++) {
    tpN.AParam[i].K =   (amix * tpA->AParam[i].K)   +
                        (bmix * tpB->AParam[i].K);
    tpN.AParam[i].th0 = (amix * tpA->AParam[i].th0) +
                        (bmix * tpB->AParam[i].th0);
  }
  for (i = 0; i < tpA->nBAH.ndihe; i++) {
    tpN.HParam[i].K =   (amix * tpA->HParam[i].K)   +
                        (bmix * tpB->HParam[i].K);
    tpN.HParam[i].Phi = (amix * tpA->HParam[i].Phi) +
                        (bmix * tpB->HParam[i].Phi);
    tpN.scee[i] = (amix * tpA->scee[i]) + (bmix * tpB->scee[i]);
    tpN.scnb[i] = (amix * tpA->scnb[i]) + (bmix * tpB->scnb[i]);
  }
  for (i = 0; i < tpA->ncmapsurf; i++) {
    tpN.MParam[i].note = (char*)malloc(MAXLINE * sizeof(char));
    sprintf(tpN.MParam[i].note, "%s (interpolation %9.5f %9.5lf)",
	    tpA->MParam[i].note, amix, bmix);
    CopyDmat(&tpN.MParam[i].esrf, &tpA->MParam[i].esrf, 0);
    CopyDmat(&tpN.MParam[i].dphi, &tpA->MParam[i].dphi, 0);
    CopyDmat(&tpN.MParam[i].dpsi, &tpA->MParam[i].dpsi, 0);
    CopyDmat(&tpN.MParam[i].d2pp, &tpA->MParam[i].d2pp, 0);
    CopyDmat(&tpN.MParam[i].Ldphi, &tpA->MParam[i].Ldphi, 0);
    CopyDmat(&tpN.MParam[i].Ldpsi, &tpA->MParam[i].Ldpsi, 0);
    CopyDmat(&tpN.MParam[i].Ld2pp, &tpA->MParam[i].Ld2pp, 0);
    CopyDmat(&tpN.MParam[i].stencil, &tpA->MParam[i].stencil, 0);
    AddDmat(&tpN.MParam[i].esrf, &tpB->MParam[i].esrf, amix, bmix);
    AddDmat(&tpN.MParam[i].dphi, &tpB->MParam[i].dphi, amix, bmix);
    AddDmat(&tpN.MParam[i].dpsi, &tpB->MParam[i].dpsi, amix, bmix);
    AddDmat(&tpN.MParam[i].d2pp, &tpB->MParam[i].d2pp, amix, bmix);
    AddDmat(&tpN.MParam[i].Ldphi, &tpB->MParam[i].Ldphi, amix, bmix);
    AddDmat(&tpN.MParam[i].Ldpsi, &tpB->MParam[i].Ldpsi, amix, bmix);
    AddDmat(&tpN.MParam[i].Ld2pp, &tpB->MParam[i].Ld2pp, amix, bmix);
    AddDmat(&tpN.MParam[i].stencil, &tpB->MParam[i].stencil, amix, bmix);
    tpN.MParam[i].muePhi = (amix * tpA->MParam[i].muePhi) +
                           (bmix * tpB->MParam[i].muePhi);
    tpN.MParam[i].muePsi = (amix * tpA->MParam[i].muePsi) +
                           (bmix * tpB->MParam[i].muePsi);
    tpN.MParam[i].mue2PP = (amix * tpA->MParam[i].mue2PP) +
                           (bmix * tpB->MParam[i].mue2PP);
    tpN.MParam[i].mueFPhi = (amix * tpA->MParam[i].mueFPhi) +
                            (bmix * tpB->MParam[i].mueFPhi);
    tpN.MParam[i].mueFPsi = (amix * tpA->MParam[i].mueFPsi) +
                            (bmix * tpB->MParam[i].mueFPsi);
    tpN.MParam[i].mueF2PP = (amix * tpA->MParam[i].mueF2PP) +
                            (bmix * tpB->MParam[i].mueF2PP);
  }  
  for (i = 0; i < tpA->natyp; i++) {
    tpN.solty[i] = (amix * tpA->solty[i]) + (bmix * tpB->solty[i]);
  }
  for (i = 0; i < 3 * tpA->nCustomVS; i++) {
    tpN.customVSd[i] = (amix * tpA->customVSd[i]) + (bmix * tpB->customVSd[i]);
  }
  for (i = 0; i < tpA->ncmapsurf; i++) {

    // Interpolate all of the CMAP matrices--don't bother recomputing
    // errors or discrepancies in the derivations.  The stencil is
    // set for a given size of CMAP, and does not need interpolation.
    AddDmat(&tpN.MParam[i].esrf, &tpB->MParam[i].esrf, amix, bmix);
    AddDmat(&tpN.MParam[i].dphi, &tpB->MParam[i].dphi, amix, bmix);
    AddDmat(&tpN.MParam[i].dpsi, &tpB->MParam[i].dpsi, amix, bmix);
    AddDmat(&tpN.MParam[i].d2pp, &tpB->MParam[i].d2pp, amix, bmix);
    AddDmat(&tpN.MParam[i].Ldphi, &tpB->MParam[i].Ldphi, amix, bmix);
    AddDmat(&tpN.MParam[i].Ldpsi, &tpB->MParam[i].Ldpsi, amix, bmix);
    AddDmat(&tpN.MParam[i].Ld2pp, &tpB->MParam[i].Ld2pp, amix, bmix);
  }

  return tpN;
}

//-----------------------------------------------------------------------------
// FreeTopology: free a topology structure.                             
//                                                                      
// Arguments:                                                           
//   tp:  the topology                                                  
//-----------------------------------------------------------------------------
void FreeTopology(prmtop *tp)
{
  int i, j;

  // Free arrays of elemental types
  free(tp->LJIdx);
  free(tp->NExcl);
  free(tp->ConExcl);
  free(tp->ExclList);
  free(tp->NBParmIdx);
  free(tp->ResLims);
  free(tp->Join);
  free(tp->Rotat);
  if (tp->ifbox > 0) {
    free(tp->Nsp);
  }
  free(tp->Charges);
  free(tp->OrigCharges);
  free(tp->ZNumber);
  free(tp->Masses);
  free(tp->InvMasses);
  free(tp->BondK);
  free(tp->BondEq);
  free(tp->BondPullK);
  free(tp->BondPullEq);
  free(tp->BondPressK);
  free(tp->BondPressEq);
  free(tp->AnglK);
  free(tp->AnglEq);
  free(tp->DiheK);
  free(tp->DiheN);
  free(tp->DihePhi);
  free(tp->LJA);
  free(tp->LJB);
  free(tp->SolA);
  free(tp->SolB);
  free(tp->HBCut);
  free(tp->Radii);
  free(tp->Screen);
  free(tp->AtomNames);
  free(tp->ResNames);
  free(tp->AtomTypes);
  free(tp->TreeSymbols);
  free(tp->scee);
  free(tp->scnb);
  free(tp->solty);
  free(tp->lVDWc);
  free(tp->norattlemask);
  free(tp->rattlemask);
  free(tp->MobileAtoms);

  // Free arrays of structs
  free(tp->BIncH);
  free(tp->BNoH);
  free(tp->AIncH);
  free(tp->ANoH);
  free(tp->HIncH);
  free(tp->HNoH);
  for (i = 0; i < tp->natom; i++) {
    free(tp->BLC[i].BC);
    free(tp->ALC[i].AC);
    for (j = 0; j < tp->HLC[i].ndihe; j++) {
      free(tp->HLC[i].HC[j].t);
    }
    free(tp->HLC[i].HC);
    free(tp->FR1Idx[i].atoms);
    if (tp->SHL[i].exe == 1) {
      free(tp->SHL[i].blist);
    }
  }
  free(tp->MLC[0].MC);
  free(tp->BLC);
  free(tp->ALC);
  free(tp->HLC);
  free(tp->MLC);
  free(tp->BParam);
  free(tp->AParam);
  free(tp->HParam);
  for (i = 0; i < tp->ncmapsurf; i++) {
    free(tp->MParam[i].note);
    DestroyDmat(&tp->MParam[i].esrf);
    DestroyDmat(&tp->MParam[i].dphi);
    DestroyDmat(&tp->MParam[i].dpsi);
    DestroyDmat(&tp->MParam[i].d2pp);
    DestroyDmat(&tp->MParam[i].Ldphi);
    DestroyDmat(&tp->MParam[i].Ldpsi);
    DestroyDmat(&tp->MParam[i].Ld2pp);
    DestroyDmat(&tp->MParam[i].stencil);
  }
  if (tp->ncmapsurf > 0) {
    free(tp->MParam);
  }
  if (tp->ncmap > 0) {
    free(tp->Mterms);
  }
  free(tp->SHL);
  for (i = 0; i < tp->ngrp; i++) {
    free(tp->lgrps[i].atoms);
  }
  free(tp->lgrps);
  free(tp->FR1Idx);
  free(tp->xtrapts);
  if (tp->neprule > 0) {
    free(tp->eprules);
  }

  // Free some pre-computed parameter tables
  DestroyDmat(&tp->LJftab);
  DestroyDmat(&tp->LJutab);

  // Destroy the connectivity array
  for (i = 0; i < tp->natom; i++) {
    free(tp->nb1234[i].L11);
    free(tp->nb1234[i].L12);
    free(tp->nb1234[i].L13);
    free(tp->nb1234[i].L14);
  }
  free(tp->nb1234);

  // Free structures for added extra points
  if (tp->ExclMarked == 1 || tp->EPInserted == 1) {
    for (i = 0; i < tp->natom; i++) {
      free(tp->ElimPair[i].list11);
      free(tp->ElimPair[i].list12);
      free(tp->ElimPair[i].list13);
      free(tp->ElimPair[i].list14);
    }
    free(tp->ElimPair);
  }
  if (tp->EPInserted == 1) {
    free(tp->OldAtomNum);
  }
  if (tp->nCustomVS > 0) {
    free(tp->customVSi);
    free(tp->customVSd);
  }
}

//-----------------------------------------------------------------------------
// RebuildCustomVS: rebuild the custom virtual site arrays ahead of topology
//                  printing.  This will ensure that all sites, both native to
//                  the topology as it was read and added since by some &rule
//                  namelist, are included in one standard format.  The older
//                  Amber extra point styles, those defined by inference from
//                  the bonding structure around the parent atom, are not
//                  counted as custom extra points / virtual sites.
//
// Arguments:
//   tp:     the system topology
//-----------------------------------------------------------------------------
static void RebuildCustomVS(prmtop *tp)
{
  int i, ncvs;

  // Free current data if custom extra points were already present
  if (tp->nCustomVS > 0) {
    free(tp->customVSi);
    free(tp->customVSd);
  }
  
  // Count the number of custom extra point frames and reallocate
  ncvs = 0;
  for (i = 0; i < tp->nxtrapt; i++) {
    ncvs += (tp->xtrapts[i].frstyle > 0 &&
	     (tp->xtrapts[i].frstyle < 7 || tp->xtrapts[i].frstyle > 10));
  }
  tp->nCustomVS = ncvs;
  tp->customVSi = (int*)malloc(6 * ncvs * sizeof(int));
  SetIVec(tp->customVSi, 6 * ncvs, -1);
  tp->customVSd = (double*)calloc(3 * ncvs, sizeof(double));

  // Repopulate the arrays
  ncvs = 0;
  for (i = 0; i < tp->nxtrapt; i++) {
    if (tp->xtrapts[i].frstyle == 0 ||\
	(tp->xtrapts[i].frstyle >= 7 && tp->xtrapts[i].frstyle <= 10)) {
      continue;
    }
    tp->customVSi[6*i]        = tp->xtrapts[i].atomid + 1;
    tp->customVSi[6*ncvs + 1] = tp->xtrapts[i].fr1 + 1;
    tp->customVSi[6*ncvs + 2] = tp->xtrapts[i].fr2 + 1;
    tp->customVSd[3*ncvs]     = tp->xtrapts[i].d1;
    if (tp->xtrapts[i].frstyle == 1) {
      tp->customVSi[6*ncvs + 5] = 4;
    }
    else if (tp->xtrapts[i].frstyle == 2) {
      tp->customVSi[6*ncvs + 3] = tp->xtrapts[i].fr3 + 1;
      tp->customVSi[6*ncvs + 5] = 5;
      tp->customVSd[3*ncvs + 1] = tp->xtrapts[i].d2;      
    }
    else if (tp->xtrapts[i].frstyle == 3) {
      tp->customVSi[6*ncvs + 3] = tp->xtrapts[i].fr3 + 1;
      tp->customVSi[6*ncvs + 5] = 6;
      tp->customVSd[3*ncvs + 1] = tp->xtrapts[i].d3;      
    }
    else if (tp->xtrapts[i].frstyle == 4) {
      tp->customVSi[6*ncvs + 3] = tp->xtrapts[i].fr3 + 1;
      tp->customVSi[6*ncvs + 5] = 7;
      tp->customVSd[3*ncvs + 1] = tp->xtrapts[i].d2;      
    }
    else if (tp->xtrapts[i].frstyle == 5) {
      tp->customVSi[6*ncvs + 3] = tp->xtrapts[i].fr3 + 1;
      tp->customVSi[6*ncvs + 5] = 8;
      tp->customVSd[3*ncvs + 1] = tp->xtrapts[i].d2;
      tp->customVSd[3*ncvs + 2] = tp->xtrapts[i].d3;
    }
    else if (tp->xtrapts[i].frstyle == 6) {
      tp->customVSi[6*ncvs + 3] = tp->xtrapts[i].fr3 + 1;
      tp->customVSi[6*ncvs + 4] = tp->xtrapts[i].fr4 + 1;
      tp->customVSi[6*ncvs + 5] = 10;
      tp->customVSd[3*ncvs + 1] = tp->xtrapts[i].d2;
      tp->customVSd[3*ncvs + 2] = tp->xtrapts[i].d3;
    }
    else if (tp->xtrapts[i].frstyle == 11) {
      tp->customVSi[6*ncvs + 5] = 9;
    }
    ncvs++;
  }
}

//-----------------------------------------------------------------------------
// FindAtom: finds an atom in a list of names, between limits il and ih.
//
// Arguments:
//   tp:      the system topology
//   i{l,h}:  limts of the atom names list (a character array) to search
//   aname:   the atom name to search for
//-----------------------------------------------------------------------------
int FindAtom(prmtop *tp, int il, int ih, char* aname)
{
  int i;

  for (i = il; i < ih; i++) {
    if (strncmp(aname, &tp->AtomNames[4*i], 4) == 0) {
      return i;
    }
  }

  // Error message
  printf("FindAtom >> Error.  Could not locate atom %.s\n", aname);
  exit(1);

  return -1;
}

//-----------------------------------------------------------------------------
// Format5e16: print a chunk of topology file with this numerical       
//             format.                                                  
//
// Arguments:
//   cname:     the subsection header (i.e. "CHARGES")
//   values:    array of real numbers to print in 5e16 format
//   N:         the length of the values array
//   outp:      output file pointer
//-----------------------------------------------------------------------------
static void Format5e16(char* cname, double* values, int N, FILE *outp)
{
  int h, i, j, slen;
  char nstr[MAXNAME];

  fprintf(outp, "%%FLAG %-74.74s\n%-80.80s\n", cname, "%FORMAT(5E16.8)");
  h = 0;
  for (i = 0; i < N; i++) {
    sprintf(nstr, "%16.8e", values[i]);
    slen = strlen(nstr);
    for (j = 0; j < slen; j++) {
      if (nstr[j] == 'e') {
        nstr[j] = 'E';
      }
    }
    fprintf(outp, "%s", nstr);
    h++;
    if (h == 5) {
      h = 0;
      fprintf(outp, "\n");
    }
  }
  if (h != 0 || N == 0) {
    fprintf(outp, "\n");
  }
}

//-----------------------------------------------------------------------------
// FormatXiY: print a chunk of topology file with this numerical format.
//
// Arguments:
//   cname:     the subsection header (i.e. "CHARGES")
//   values:    array of integers to print in 10i8 format
//   N:         the length of the values array
//   outp:      output file pointer
//   x:         the number of integers per line
//   y:         the maximum length of each integer
//-----------------------------------------------------------------------------
static void FormatXiY(char* cname, int* values, int N, FILE *outp, int x,
                      int y)
{
  int h, i;
  char fmtstring[32];

  sprintf(fmtstring, "%%FORMAT(%dI%d)", x, y);
  fprintf(outp, "%%FLAG %-74.74s\n%-80.80s\n", cname, fmtstring);
  h = 0;
  for (i = 0; i < N; i++) {
    fprintf(outp, "%*d", y, values[i]);
    h++;
    if (h == x) {
      h = 0;
      fprintf(outp, "\n");
    }
  }
  if (h != 0 || N == 0) {
    fprintf(outp, "\n");
  }
}

//-----------------------------------------------------------------------------
// FormatXfYpZ: print a chunk of topology file with this numerical format.
//
// Arguments:
//   cname:     the subsection header (i.e. "CHARGES")
//   values:    array of real numbers to print in X * %Y.Zlf format
//   N:         the length of the values array
//   comment:   optional comment to print just before the format line (if set
//              to NULL no comment will be printed)
//   outp:      output file pointer
//   x:         the number of real numbers per line
//   y:         the overall format length of each real number
//   z:         the number of digits after the decimal
//-----------------------------------------------------------------------------
static void FormatXfYpZ(char* cname, double* values, int N, char* comment,
                        FILE *outp, int x, int y, int z)
{
  int h, i;
  char fmtstring[32];

  sprintf(fmtstring, "%%FORMAT(%dF%d.%d)", x, y, z);
  if (comment == NULL) {
    fprintf(outp, "%%FLAG %-74.74s\n%-80.80s\n", cname, fmtstring);
  }
  else {
    fprintf(outp, "%%FLAG %-74.74s\n%%COMMENT  %-70.70s\n%-80.80s\n", cname,
            comment, fmtstring);
  }
  h = 0;
  for (i = 0; i < N; i++) {
    fprintf(outp, "%*.*lf", y, z, values[i]);
    h++;
    if (h == x) {
      h = 0;
      fprintf(outp, "\n");
    }
  }
  if (h != 0 || N == 0) {
    fprintf(outp, "\n");
  }
}

//-----------------------------------------------------------------------------
// Format20a4: print a chunk of a topology file with this character format.
//
// Arguments:
//   cname:     the subsection header (i.e. "CHARGES")
//   values:    array of characters to print in 20a4 format
//   N:         the length of the values array, divided by 4 (there are N
//              groups of four cahracters each)
//   outp:      output file pointer
//-----------------------------------------------------------------------------
static void Format20a4(char* cname, char* values, int N, FILE *outp)
{
  int h, i;

  fprintf(outp, "%%FLAG %-74.74s\n%-80.80s\n", cname, "%FORMAT(20a4)");
  h = 0;
  for (i = 0; i < N; i++) {
    fprintf(outp, "%.4s", &values[4*i]);
    h++;
    if (h == 20) {
      h = 0;
      fprintf(outp, "\n");
    }
  }
  if (h != 0 || N == 0) {
    fprintf(outp, "\n");
  }
}

//-----------------------------------------------------------------------------
// PutPrmTop: put a topology into a file.                               
//
// Arguments:
//   tp:      system topology
//   fname:   name of the file to write
//   title:   the title of this prmtop to be written on the second line (right
//            after the version stamp)
//-----------------------------------------------------------------------------
void PutPrmTop(prmtop *tp, char* fname, char* title)
{
  int i;
  int* itmp;
  double baugsum;
  double* dtmp;
  double* solty;
  char numstr[32];
  FILE *outp = fopen(fname, "w");
  
  fprintf(outp, "%s%-80.80s\n%-80.80s\n%-80.80s\n", tp->vstamp, "%FLAG TITLE",
	  "%FORMAT(20a4)", title);
  fprintf(outp, "%-80.80s\n%-80.80s\n", "%FLAG POINTERS", "%FORMAT(10I8)");
  fprintf(outp, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n"
	  "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n"
	  "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n%8d", tp->natom,
	  tp->ntypes, tp->withH.nbond, tp->woH.nbond, tp->withH.nangl,
	  tp->woH.nangl, tp->withH.ndihe, tp->woH.ndihe, tp->nhparm,
	  tp->nparm, tp->tnexcl, tp->nres, tp->woHC.nbond, tp->woHC.nangl,
	  tp->woHC.ndihe, tp->nBAH.nbond, tp->nBAH.nangl, tp->nBAH.ndihe,
	  tp->natyp, tp->nphb, tp->ifpert, tp->pert.nbond, tp->pert.nangl,
	  tp->pert.ndihe, tp->wpert.nbond, tp->wpert.nangl, tp->wpert.ndihe,
	  tp->ifbox, tp->nmxrs, tp->ifcap, tp->numextra);
  if (tp->ncopy > 0) {
    fprintf(outp, "%8d\n", tp->ncopy);
  }
  else {
    fprintf(outp, "\n");
  }
  
  // Print atom names, charges, masses, and atomic numbers  
  Format20a4("ATOM_NAME", tp->AtomNames, tp->natom, outp);
  Format5e16("CHARGE", tp->OrigCharges, tp->natom, outp);
  FormatXiY("ATOMIC_NUMBER", tp->ZNumber, tp->natom, outp, 10, 8);
  Format5e16("MASS", tp->Masses, tp->natom, outp);

  // Print the atom type index, reverting to Fortran numbering
  itmp = CpyIVec(tp->LJIdx, tp->natom);
  for (i = 0; i < tp->natom; i++) {
    if (itmp[i] < 0) {
      itmp[i] = -itmp[i];
    }
  }
  IVecAdd(itmp, tp->natom, 1);
  FormatXiY("ATOM_TYPE_INDEX", itmp, tp->natom, outp, 10, 8);
  free(itmp);

  // Print the number of excluded atoms
  itmp = CpyIVec(tp->NExcl, tp->natom);
  for (i = 0; i < tp->natom; i++) {
    if (itmp[i] == 0) {
      itmp[i] = 1;
    }
  }
  FormatXiY("NUMBER_EXCLUDED_ATOMS", itmp, tp->natom, outp, 10, 8);
  free(itmp);

  // Print the non-bonded parameter index, again adjusting back to Fortran
  itmp = CpyIVec(tp->NBParmIdx, tp->ntypes*tp->ntypes);
  IVecAdd(itmp, tp->ntypes*tp->ntypes, 1);
  FormatXiY("NONBONDED_PARM_INDEX", itmp, tp->ntypes*tp->ntypes, outp, 10, 8);
  free(itmp);

  // Print residue labels and pointers, adjusting the latter for Fortran
  Format20a4("RESIDUE_LABEL", tp->ResNames, tp->nres, outp);
  itmp = CpyIVec(tp->ResLims, tp->nres);
  IVecAdd(itmp, tp->nres, 1);
  FormatXiY("RESIDUE_POINTER", itmp, tp->nres, outp, 10, 8);
  free(itmp);

  // Print valence parameter information, testing some fields for relevance
  Format5e16("BOND_FORCE_CONSTANT", tp->BondK, tp->nBAH.nbond, outp);
  Format5e16("BOND_EQUIL_VALUE", tp->BondEq, tp->nBAH.nbond, outp);
  if (DAbsSum(tp->BondPullK, tp->nBAH.nbond) > 1.0e-8) {
    Format5e16("BOND_STIFFNESS_PULL_ADJ", tp->BondPullK, tp->nBAH.nbond, outp);
    Format5e16("BOND_EQUIL_PULL_ADJ", tp->BondPullEq, tp->nBAH.nbond, outp);
  }
  if (DAbsSum(tp->BondPressK, tp->nBAH.nbond) > 1.0e-8) {
    Format5e16("BOND_STIFFNESS_PRESS_ADJ", tp->BondPressK, tp->nBAH.nbond,
               outp);
    Format5e16("BOND_EQUIL_PRESS_ADJ", tp->BondPressEq, tp->nBAH.nbond, outp);
  }
  Format5e16("ANGLE_FORCE_CONSTANT", tp->AnglK, tp->nBAH.nangl, outp);
  Format5e16("ANGLE_EQUIL_VALUE", tp->AnglEq, tp->nBAH.nangl, outp);
  Format5e16("DIHEDRAL_FORCE_CONSTANT", tp->DiheK, tp->nBAH.ndihe, outp);
  Format5e16("DIHEDRAL_PERIODICITY", tp->DiheN, tp->nBAH.ndihe, outp);
  Format5e16("DIHEDRAL_PHASE", tp->DihePhi, tp->nBAH.ndihe, outp);

  // The 1:4 exclusion values must be adjuste back towards their original
  // prmtop values, rather than the mdgx pre-computations done in the array.
  dtmp = (double*)malloc(tp->nBAH.ndihe * sizeof(double));
  for (i = 0; i < tp->nBAH.ndihe; i++) {
    dtmp[i] = (fabs(tp->scee[i]) > 1.0e-8) ? 1.0 / (1.0 - tp->scee[i]) : 0.0;
  }
  Format5e16("SCEE_SCALE_FACTOR", dtmp, tp->nBAH.ndihe, outp);
  for (i = 0; i < tp->nBAH.ndihe; i++) {
    dtmp[i] = (fabs(tp->scnb[i]) > 1.0e-8) ? 1.0 / (1.0 - tp->scnb[i]) : 0.0;
  }
  Format5e16("SCNB_SCALE_FACTOR", dtmp, tp->nBAH.ndihe, outp);
  free(dtmp);
  
  // The SOLTY field is reserved for future use
  Format5e16("SOLTY", tp->solty, tp->natyp, outp);

  // Print non-bonded parameters
  Format5e16("LENNARD_JONES_ACOEF", tp->LJA, tp->ntypes*(tp->ntypes+1)/2,
             outp);
  Format5e16("LENNARD_JONES_BCOEF", tp->LJB, tp->ntypes*(tp->ntypes+1)/2,
	     outp);

  // Print valence term indexing, making adjustments as needed to return
  // to the original, Fortran-centric, Amber prmtop format
  itmp = CpyIVec((int*)tp->BIncH, 3*tp->withH.nbond);
  ResetBondArray(itmp, tp->withH.nbond, 3);
  FormatXiY("BONDS_INC_HYDROGEN", itmp, 3*tp->withH.nbond, outp, 10, 8);
  free(itmp);
  itmp = CpyIVec((int*)tp->BNoH, 3*tp->woH.nbond);
  ResetBondArray(itmp, tp->woH.nbond, 3);
  FormatXiY("BONDS_WITHOUT_HYDROGEN", itmp, 3*tp->woH.nbond, outp, 10, 8);
  free(itmp);
  itmp = CpyIVec((int*)tp->AIncH, 4*tp->withH.nangl);
  ResetBondArray(itmp, tp->withH.nangl, 4);
  FormatXiY("ANGLES_INC_HYDROGEN", itmp, 4*tp->withH.nangl, outp, 10, 8);
  free(itmp);
  itmp = CpyIVec((int*)tp->ANoH, 4*tp->woH.nangl);
  ResetBondArray(itmp, tp->woH.nangl, 4);
  FormatXiY("ANGLES_WITHOUT_HYDROGEN", itmp, 4*tp->woH.nangl, outp, 10, 8);
  free(itmp);
  itmp = CpyIVec((int*)tp->HIncH, 5*tp->withH.ndihe);
  ResetBondArray(itmp, tp->withH.ndihe, 5);
  FormatXiY("DIHEDRALS_INC_HYDROGEN", itmp, 5*tp->withH.ndihe, outp, 10, 8);
  free(itmp);
  itmp = CpyIVec((int*)tp->HNoH, 5*tp->woH.ndihe);
  ResetBondArray(itmp, tp->woH.ndihe, 5);
  FormatXiY("DIHEDRALS_WITHOUT_HYDROGEN", itmp, 5*tp->woH.ndihe, outp, 10, 8);
  free(itmp);

  // Print the excluded atoms list, again adjusting back to Fortran numbering
  itmp = CpyIVec(tp->ExclList, tp->ConExcl[tp->natom]);
  IVecAdd(itmp, tp->ConExcl[tp->natom], 1);
  FormatXiY("EXCLUDED_ATOMS_LIST", itmp, tp->ConExcl[tp->natom], outp, 10, 8);
  free(itmp);
  Format5e16("HBOND_ACOEF", tp->SolA, tp->nphb*(tp->nphb+1)/2, outp);
  Format5e16("HBOND_BCOEF", tp->SolB, tp->nphb*(tp->nphb+1)/2, outp);
  Format5e16("HBCUT", tp->HBCut, tp->nphb*(tp->nphb+1)/2, outp);
  Format20a4("AMBER_ATOM_TYPE", tp->AtomTypes, tp->natom, outp);
  Format20a4("TREE_CHAIN_CLASSIFICATION", tp->TreeSymbols, tp->natom, outp);
  FormatXiY("JOIN_ARRAY", tp->Join, tp->natom, outp, 10, 8);
  FormatXiY("IROTAT", tp->Rotat, tp->natom, outp, 10, 8);
  if (tp->ifbox > 0) {
    fprintf(outp, "%%FLAG %-74.74s\n%-80.80s\n", "SOLVENT_POINTERS",
	    "%FORMAT(3I8)");
    fprintf(outp, "%8d%8d%8d\n", tp->iptres, tp->nspm, tp->nspsol);
    FormatXiY("ATOMS_PER_MOLECULE", tp->Nsp, tp->nspm, outp, 10, 8);
    solty = (double*)malloc(4 * sizeof(double));
    solty[0] = tp->gdim[4];
    solty[1] = tp->gdim[0];
    solty[2] = tp->gdim[1];
    solty[3] = tp->gdim[2];
    Format5e16("BOX_DIMENSIONS", solty, 4, outp);
    free(solty);
  }
  fprintf(outp, "%-80.80s\n%-80.80s\n%-80.80s\n", "%FLAG RADIUS_SET",
	  "%FORMAT(1a80)", tp->RadSet);
  Format5e16("RADII", tp->Radii, tp->natom, outp);
  Format5e16("SCREEN", tp->Screen, tp->natom, outp);

  // Regenerate the virtual site frame arrays, if needed
  if (tp->numextra > 0) {
    RebuildCustomVS(tp);
    FormatXiY("VIRTUAL_SITE_FRAMES", tp->customVSi, 6 * tp->nCustomVS, outp,
              10, 8);
    Format5e16("VIRTUAL_SITE_FRAME_DETAILS", tp->customVSd, 3 * tp->nCustomVS,
               outp);
  }
  fprintf(outp, "%-80.80s\n%-80.80s\n%8d\n", "%FLAG IPOL", "%FORMAT(1I8)",
	  tp->ipol);

  // CMAP information goes here
  if (tp->ncmap > 0) {
    sprintf(numstr, "%8d%8d", tp->ncmap, tp->ncmapsurf); 
    fprintf(outp, "%-80.80s\n%-80.80s\n%-80.80s\n", "%FLAG CMAP_COUNT",
	    "%FORMAT(2I8)", numstr);
    itmp = (int*)malloc(tp->ncmapsurf * sizeof(int));
    for (i = 0; i < tp->ncmapsurf; i++) {
      itmp[i] = tp->MParam[i].esrf.row;
    }
    FormatXiY("CMAP_RESOLUTION", itmp, tp->ncmapsurf, outp, 20, 4);
    free(itmp);
    for (i = 0; i < tp->ncmapsurf; i++) {
      sprintf(numstr, "CMAP_PARAMETER_%02d", i + 1);
      FormatXfYpZ(numstr, tp->MParam[i].esrf.data,
                  tp->MParam[i].esrf.row * tp->MParam[i].esrf.col,
                  tp->MParam[i].note, outp, 8, 9, 5);
    }
    itmp = (int*)malloc(6 * tp->ncmap * sizeof(int));
    for (i = 0; i < tp->ncmap; i++) {
      itmp[6*i    ] = tp->Mterms[i].a + 1;
      itmp[6*i + 1] = tp->Mterms[i].b + 1;
      itmp[6*i + 2] = tp->Mterms[i].c + 1;
      itmp[6*i + 3] = tp->Mterms[i].d + 1;
      itmp[6*i + 4] = tp->Mterms[i].e + 1;
      itmp[6*i + 5] = tp->Mterms[i].idx + 1;
    }
    FormatXiY("CMAP_INDEX", itmp, 6 * tp->ncmap, outp, 6, 8);
  }

  // Close the file
  fclose(outp);
}

//-----------------------------------------------------------------------------
// Chirality: get the chirality of a center, given the number CA of the 
//            atom at the center and four atoms around it.  We expect A 
//            to be the first atom in the list, B the second, C the     
//            third, and D the fourth atom in order of greatest to least
//            chiral dominance.                                         
//
// Arguments:
//   tc:      coordinates to analyze
//   CA:      index of the possibly chiral atom
//   A,B,C,D: atoms to which the atom CA is bonded
//-----------------------------------------------------------------------------
static int Chirality(coord *tc, int CA, int A, int B, int C, int D)
{
  int i;
  double cax[3], ax[3], bx[3], cx[3], dx[3], ra[3], rb[3], rc[3], rd[3];
  double prja[3], prjb[3], prjc[3], acrb[3], bcrc[3];

  // Get coordinates
  for (i = 0; i < 3; i++) {
    cax[i] = tc->loc[3*CA+i];
    ax[i] = tc->loc[3*A+i];
    bx[i] = tc->loc[3*B+i];
    cx[i] = tc->loc[3*C+i];
    dx[i] = tc->loc[3*D+i];
    ra[i] = ax[i] - cax[i];
    rb[i] = bx[i] - cax[i];
    rc[i] = cx[i] - cax[i];
    rd[i] = dx[i] - cax[i];
  }

  // Remove projection
  Project(ra, rd, prja, 3);
  Project(rb, rd, prjb, 3);
  Project(rc, rd, prjc, 3);
  for (i = 0; i < 3; i++) {
    ra[i] -= prja[i];
    rb[i] -= prjb[i];
    rc[i] -= prjc[i];
  }

  // Find chirality via cross products
  CrossP(ra, rb, acrb);
  CrossP(rb, rc, bcrc);

  // L chirality is when ACRB and BCRC point opposite RD
  if (DotP(acrb, rd, 3) < 0.0 && DotP(bcrc, rd, 3) < 0.0) {
    return 1;
  }

  // Otherwise, this has D chirality
  return -1;
}

//-----------------------------------------------------------------------------
// ProteinChiralityCheck: checks the chiral centers of amino acids.  Any
//                        residues with D chirality are reported.       
//
// Arguments:
//   tp:      system topology
//   tc:      system coordinates
//   outp:    output file being written (mdout)
//-----------------------------------------------------------------------------
int ProteinChiralityCheck(prmtop *tp, coord *tc, FILE *outp)
{
  int i, nca, nn, nc, ncb, nha, iwarning;
  char* ctmp;

  iwarning = 0;
  for (i = 0; i < tp->nres; i++) {
    ctmp = &tp->ResNames[4*i];
    if (!(strncmp(ctmp, "ALA", 3) == 0 || strncmp(ctmp, "ARG", 3) == 0 ||
          strncmp(ctmp, "ASN", 3) == 0 || strncmp(ctmp, "ASP", 3) == 0 ||
          strncmp(ctmp, "ASH", 3) == 0 || strncmp(ctmp, "GLH", 3) == 0 ||
          strncmp(ctmp, "CYS", 3) == 0 || strncmp(ctmp, "CYX", 3) == 0 ||
          strncmp(ctmp, "GLN", 3) == 0 || strncmp(ctmp, "GLU", 3) == 0 ||
          strncmp(ctmp, "HID", 3) == 0 || strncmp(ctmp, "HIE", 3) == 0 ||
          strncmp(ctmp, "HIP", 3) == 0 || strncmp(ctmp, "HIS", 3) == 0 ||
          strncmp(ctmp, "ILE", 3) == 0 || strncmp(ctmp, "LEU", 3) == 0 ||
          strncmp(ctmp, "LYS", 3) == 0 || strncmp(ctmp, "MET", 3) == 0 ||
          strncmp(ctmp, "PHE", 3) == 0 || strncmp(ctmp, "PRO", 3) == 0 ||
          strncmp(ctmp, "SER", 3) == 0 || strncmp(ctmp, "THR", 3) == 0 ||
          strncmp(ctmp, "TRP", 3) == 0 || strncmp(ctmp, "TYR", 3) == 0 ||
          strncmp(ctmp, "VAL", 3) == 0)) {
      continue;
    }
    
    // If we're still here, we've got some amino acid
    nca = FindAtom(tp, tp->ResLims[i], tp->ResLims[i+1], "CA  ");
    nn = FindAtom(tp, tp->ResLims[i], tp->ResLims[i+1], "N   ");
    nc = FindAtom(tp, tp->ResLims[i], tp->ResLims[i+1], "C   ");
    ncb = FindAtom(tp, tp->ResLims[i], tp->ResLims[i+1], "CB  ");
    nha = FindAtom(tp, tp->ResLims[i], tp->ResLims[i+1], "HA  ");

    // Now, get the chirality of these four atoms
    if (Chirality(tc, nca, nn, nc, ncb, nha) == -1) {
      fprintf(outp, " Warning.  Residue %.4s %d has D chirality at its CA "
              "atom.\n", &tp->ResNames[4*i], i);
      iwarning++;
    }
  }

  return iwarning;
}

//-----------------------------------------------------------------------------
// FindDisulfides: find sulfur atoms on CYS residues that are very close
//                 together and prompt the user if they are not bonded. 
//
// Arguments:
//   tp:      system topology
//   tc:      system coordinates
//   outp:    output file being written (mdout)
//-----------------------------------------------------------------------------
int FindDisulfides(prmtop *tp, coord *tc, FILE *outp)
{
  int i, j, k, iresno, jresno, bonded, ndss, iwarning;
  double dx[3];
  char wrnmsg[256];

  ndss = 0;
  iwarning = 0;
  for (i = 0; i < tp->natom; i++) {
    if ((tp->Masses[i] <= 32.0 || tp->Masses[i] >= 33.0) &&
        tp->AtomNames[4*i] != 'S') {
      continue;
    }

    // This looks like sulfur, but is it part of a cysteine?
    iresno = LocateResID(tp, i, 0, tp->nres);
    if (strncmp(&tp->ResNames[4*iresno], "CY", 2) != 0 &&
        strncmp(&tp->ResNames[4*iresno+1], "CY", 2) != 0) {
      continue;
    }
    for (j = i+1; j < tp->natom; j++) {
      if ((tp->Masses[j] <= 32.0 || tp->Masses[j] >= 33.0) &&
          tp->AtomNames[4*j] != 'S') {
        continue;
      }
      jresno = LocateResID(tp, j, 0, tp->nres);
      if (strncmp(&tp->ResNames[4*jresno], "CY", 2) != 0 &&
          strncmp(&tp->ResNames[4*jresno+1], "CY", 2) != 0) {
        continue;
      }

      // So, we've got two atoms that look interesting
      dx[0] = tc->loc[3*i] - tc->loc[3*j];
      dx[1] = tc->loc[3*i+1] - tc->loc[3*j+1];
      dx[2] = tc->loc[3*i+2] - tc->loc[3*j+2];
      if (dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] < 16.0) {

        // They're close together... are they bonded?
        bonded = 0;
        for (k = 0; k < tp->woH.nbond; k++) {
          if ((tp->BNoH[k].a == i && tp->BNoH[k].b == j) ||
              (tp->BNoH[k].b == i && tp->BNoH[k].a == j)) {
            bonded = 1;
            ndss++;
            break;
          }
        }
        if (bonded == 0) {
          sprintf(wrnmsg, " Warning.  Atoms %.4s %4d %.4s and %.4s %4d %.4s "
                  "may be part of an unfulfilled disulfide bridge.",
                  &tp->ResNames[4*iresno], iresno+1, &tp->AtomNames[4*i],
                  &tp->ResNames[4*jresno], jresno+1, &tp->AtomNames[4*j]);
          PrintParagraph(wrnmsg, 79, NULL, outp);
          iwarning++;
        }
      }
    }
  }
  if (ndss > 0) {
    fprintf(outp, " There are %d disulfide bonds in the system.\n", ndss);
  }

  return iwarning;
}

//-----------------------------------------------------------------------------
// PrintMatchingTypeNames: print the type names of atoms matching a particular
//                         index.                            
//                                                                      
// Arguments:                                                           
//   tp:        the topology                                            
//   idx:       the type index to match against                         
//   outp:      the output file currently being written                 
//-----------------------------------------------------------------------------
static void PrintMatchingTypeNames(prmtop *tp, int idx, FILE *outp)
{
  int i, j, nmatch, maxmatch, atmex;
  cmat C;

  C = CreateCmat(32, 4);
  nmatch = 0;
  maxmatch = 32;
  for (i = 0; i < tp->natom; i++) {
    if (tp->LJIdx[i] != idx) {
      continue;
    }
    atmex = 0;
    for (j = 0; j < nmatch; j++) {
      if (str4cmp(C.map[j], &tp->AtomTypes[4*i]) == 0) {
        atmex = 1;
        break;
      }
    }
    if (atmex == 0) {
      strncpy(C.map[nmatch], &tp->AtomTypes[4*i], 4);
      nmatch++;
      if (nmatch == maxmatch-1) {
        maxmatch += 32;
        C = ReallocCmat(&C, maxmatch, 4);
      }
    }
  }
  fprintf(outp, "%s\n", C.data);

  // Free allocated memory
  DestroyCmat(&C);
}

//-----------------------------------------------------------------------------
// CheckLJRules: check Lennard-Jones combining rules, report any that do not
//               conform to standard Lorentz-Berthelot or geometric combining
//               rules.                                       
//                                                                      
// Arguments:                                                           
//   tp:        the topology                                            
//   outp:      the output file (currently writing the opening segment) 
//-----------------------------------------------------------------------------
int CheckLJRules(prmtop *tp, FILE *outp)
{
  int i, j, iwarn, lbrule, gmrule;
  int* idcon;
  double A, B, sig, eps, esig, eeps;
  double* sigV;
  double* epsV;

  // No problems initially
  iwarn = 0;

  // Scratch arrays for type parameters
  sigV = (double*)malloc(tp->ntypes*sizeof(double));
  epsV = (double*)malloc(tp->ntypes*sizeof(double));

  // Count the occurences of each type
  idcon = (int*)calloc(tp->ntypes, sizeof(int));
  for (i = 0; i < tp->natom; i++) {
    if (tp->LJIdx[i] >= 0) {
      idcon[tp->LJIdx[i]] += 1;
    }
  }

  // Scan through each atom type
  for (i = 0; i < tp->ntypes; i++) {

    // Determine sigma and epsilon for this type with itself
    if (tp->NBParmIdx[i*tp->ntypes+i] >= 0) {
      A = tp->LJA[tp->NBParmIdx[i*tp->ntypes+i]];
      B = tp->LJB[tp->NBParmIdx[i*tp->ntypes+i]];
      sigV[i] = (B > 1.0e-04) ? pow(A/B, 1.0/6.0) : 0.0;
      epsV[i] = (sigV[i] > 1.0e-04) ?
        0.25*tp->LJB[tp->NBParmIdx[i*tp->ntypes+i]]/pow(sigV[i], 6.0) : 0.0;
    }
  }

  // Now check through each pair
  lbrule = 0;
  gmrule = 0;
  for (i = 0; i < tp->ntypes-1; i++) {
    if (idcon[i] <= 0) {
      continue;
    }
    for (j = i+1; j < tp->ntypes; j++) {
      if (idcon[j] <= 0) {
        continue;
      }
      if (tp->NBParmIdx[i*tp->ntypes+j] >= 0) {
        A = tp->LJA[tp->NBParmIdx[i*tp->ntypes+j]];
        B = tp->LJB[tp->NBParmIdx[i*tp->ntypes+j]];
        sig = (B > 1.0e-04) ? pow(A/B, 1.0/6.0) : 0.0;
        eps = (sig > 1.0e-04) ?
          0.25*tp->LJB[tp->NBParmIdx[i*tp->ntypes+j]]/pow(sig, 6.0) : 0.0;
      }
      else {
        continue;
      }
      if (fabs(sig - 0.5*(sigV[i]+sigV[j])) < 1.0e-6 &&
          fabs(eps - sqrt(epsV[i]*epsV[j])) < 1.0e-6) {

        // This is a Lorentz-Berthelot pair
        lbrule++;
      }
      else if (fabs(sig - sqrt(sigV[i]*sigV[j])) < 1.0e-6 &&
               fabs(eps - sqrt(epsV[i]*epsV[j])) < 1.0e-6) {

        // This is a geometric pair
        gmrule++;
      }
      else {

        // This pair breaks known combining rules
        fprintf(outp, " - Nonstandard LJ pair detected.\n");
        fprintf(outp, "     Atom Types (i): ");
        PrintMatchingTypeNames(tp, i, outp);
        fprintf(outp, "     Atom Types (j): ");
        PrintMatchingTypeNames(tp, j, outp);
        fprintf(outp, "     Effective [ Sigma, Epsilon ]:  %9.6lf, %9.6lf\n",
                sig, eps);
        if (lbrule > 0) {
          esig = 0.5*(sigV[i]+sigV[j]);
          eeps = sqrt(epsV[i]*epsV[j]);
        }
        else if (gmrule > 0) {
          esig = sqrt(sigV[i]*sigV[j]);
          eeps = sqrt(epsV[i]*epsV[j]);
        }
        else {
          esig = -1.0;
          eeps = -1.0;
        }
        if (esig >= 0.0 && eeps >= 0.0) {
          fprintf(outp, "     Expected  [ Sigma, Epsilon ]:  %9.6lf, %9.6lf\n",
                  esig, eeps);
        }
        else {
          fprintf(outp, "     No expectation value can yet be computed.\n");
        }
        iwarn = 1;
      }
    }
  }

  // Check that geometric and LB combining rules are not competing
  if (gmrule > 0 && lbrule > 0) {
    fprintf(outp, " - Both geometric and Lorentz-Berthelot combining rules "
            "seem to be in use.\n");
    iwarn = 1;
  }
  else if (gmrule > 0) {
    fprintf(outp, " - Geometric combining rules are in effect.\n");
  }
  else if (lbrule > 0) {
    fprintf(outp, " - Lorentz-Berthelot combining rules are in effect.\n");
  }

  // Free allocated memory
  free(sigV);
  free(epsV);
  free(idcon);

  return iwarn;
}

//-----------------------------------------------------------------------------
// LocateResID: locate the number of a residue containing atom ai within
//              topology tp.  This function uses a binary search.       
//                                                                      
// Arguments:                                                           
//   tp:        the topology                                            
//   ai:        the atom number                                         
//   [l,h]lim:  the lower and upper limits of the residues to search    
//-----------------------------------------------------------------------------
int LocateResID(prmtop *tp, int ai, int llim, int hlim)
{
  int cres;

  // First, test the middle residue
  while (hlim > llim) {
    cres = (hlim+llim)/2;
    if (ai >= tp->ResLims[cres] && ai < tp->ResLims[cres+1]) {
      return cres;
    }
    else if (ai < tp->ResLims[cres]) {
      hlim = (cres > llim) ? cres-1 : cres;
    }
    else if (ai >= tp->ResLims[cres+1]) {
      llim = cres+1;
    }
  }
  if (ai < tp->ResLims[llim] || ai >= tp->ResLims[hlim+1]) {
    printf("LocateResID >> Error.  Atom %d not located.\n", ai);
  }
  if (ai >= tp->ResLims[llim] && ai < tp->ResLims[llim+1]) {
    return llim;
  }
  return hlim;
}

//-----------------------------------------------------------------------------
// LongRangeVDW: integrate, over a distance from the non-bonded cutoff to
//               infinity, the long-range van-der Waals correction appropriate
//               to each atom type in this topology.        
//                                                                      
// Arguments:                                                           
//   tp:    the topology                                                
//   dcinp: the direct space control data                               
//-----------------------------------------------------------------------------
void LongRangeVDW(prmtop *tp, dircon *dcinp)
{
  int i, j, k, Nij;
  double lVDW6, lVDW12, eps, sig, virial, sfac, sig6, sig12, invr3, invr9;
  double lja, ljb;

  // Allocate the array of correction factors
  tp->lVDWc = (double*)malloc(3*tp->ntypes*sizeof(double));
  invr3 = 1.0/pow(dcinp->Vcut, 3.0);
  invr9 = invr3*invr3*invr3;
  for (i = 0; i < tp->ntypes; i++) {

    // Compute the interaction correction, sans unit cell volume
    lVDW6 = 0.0;
    lVDW12 = 0.0;
    virial = 0.0;
    for (j = 0; j < tp->ntypes; j++) {
      if (tp->NBParmIdx[i*tp->ntypes+j] >= 0) {
        lja = tp->LJA[tp->NBParmIdx[i*tp->ntypes+j]];
        ljb = tp->LJB[tp->NBParmIdx[i*tp->ntypes+j]];
        sig = (ljb > 1.0e-04) ? pow(lja/ljb, 1.0/6.0) : 0.0;
        eps = (sig > 1.0e-04) ? 
          0.25*tp->LJB[tp->NBParmIdx[i*tp->ntypes+j]]/pow(sig, 6.0) : 0.0;
      }
      else {
        continue;
      }

      // If we're still here, there are significant interactions
      Nij = 0;
      for (k = 0; k < tp->natom; k++) {
        if (tp->LJIdx[k] == j) {
          Nij++;
        }
      }
      sfac = 8.0*PI*eps*Nij;
      sig6 = pow(sig, 6.0);
      sig12 = sig6*sig6;
      lVDW6 -= sfac*sig6*invr3/3.0;
      lVDW12 += sfac*sig12*invr9/9.0;
      virial += sfac*(2.0*sig6*invr3 - (4.0/3.0)*sig12*invr9);
    }
    tp->lVDWc[3*i] = lVDW6;
    tp->lVDWc[3*i+1] = lVDW12;
    tp->lVDWc[3*i+2] = virial;
  }
}

