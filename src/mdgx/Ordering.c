//-----------------------------------------------------------------------------
// Divider(style): this highly optimized routine will go over any region within
//                 a cell and compute how its atoms participate in up to four
//                 types of interactions, depending on the variant of the
//                 function.  Counters for the atom lists accumulated by this
//                 function are initialized before calling it.  Variants with
//                 2, 4, and 7 interactions are compiled.
//
// Arguments:                                                           
//   C:             the cell in which we're operating                   
//   sctrS:         the cell sector in which we're operating            
//   fcn[A...G]:    the interface unit normal vectors                   
//   cploc:         the cell center point location                      
//   invcut:        metrics for scaling the atom-to-plane distances and 
//   invEEcut:        thus obtaining a bin number                       
//   metric:                                                            
//   cordr:         the accumulating atom order array                   
//   cnsr[A...G]:   the lengths of lists stored in cordr; note that     
//                  cnsr[A...G] point to elements of the array cnsr in  
//                  the Ranger function.  Elements of the larger cnsr   
//                  array are graduated by the maximum number of atoms  
//                  per cell so that it indexes directly into cordr     
//   offset:        the number to add to the cell partition index;      
//                  because cell partitions points to regular segments  
//                  of the cell data aray, the kth atom of partition P  
//                  is the (k+maxatom)th atom of partition (P-1) and the
//                  (maxatom-k)th atom of partition P+1.  Adjusting the 
//                  atom numbering in this manner allows multiple       
//                  partitions to be treated as one.                    
//-----------------------------------------------------------------------------
#if VGTE == 1
  #if MULTIPLE == 1
static void Divider1VgtE(cell *C, int sctrS, double *fcnA, double *cploc,
			 double invcut, double invEEcut, int* cordr,
			 int *cnsrA, int offset)
  #elif MULTIPLE == 2
static void Divider2VgtE(cell *C, int sctrS, double *fcnA, double *fcnB,
			 double *cploc, double invcut, double invEEcut,
			 int* cordr, int *cnsrA, int *cnsrB, int offset)
  #elif MULTIPLE == 4
static void Divider4VgtE(cell *C, int sctrS, double *fcnA, double *fcnB,
			 double *fcnC, double *fcnD, double *cploc,
			 double invcut, double invEEcut, int* cordr,
			 int *cnsrA, int *cnsrB, int *cnsrC, int *cnsrD,
			 int offset)
  #else
static void Divider7VgtE(cell *C, int sctrS, double *fcnA, double *fcnB,
			 double *fcnC, double *fcnD, double *fcnE,
                         double *fcnF, double *fcnG, double *cploc,
			 double invcut, double invEEcut, int* cordr,
			 int *cnsrA, int *cnsrB, int *cnsrC, int *cnsrD,
			 int *cnsrE, int *cnsrF, int *cnsrG, int offset)
  #endif
#else
  #if MULTIPLE == 1
static void Divider1VeqE(cell *C, int sctrS, double *fcnA, double *cploc,
			 double metric, int* cordr, int *cnsrA, int offset)
  #elif MULTIPLE == 2
static void Divider2VeqE(cell *C, int sctrS, double *fcnA, double *fcnB,
			 double *cploc, double metric, int* cordr, int *cnsrA,
			 int *cnsrB, int offset)
  #elif MULTIPLE == 4
static void Divider4VeqE(cell *C, int sctrS, double *fcnA, double *fcnB,
			 double *fcnC, double *fcnD, double *cploc,
			 double metric, int* cordr, int *cnsrA, int *cnsrB,
			 int *cnsrC, int *cnsrD, int offset)
  #else
static void Divider7VeqE(cell *C, int sctrS, double *fcnA, double *fcnB,
			 double *fcnC, double *fcnD, double *fcnE,
                         double *fcnF, double *fcnG, double *cploc,
			 double metric, int* cordr, int *cnsrA, int *cnsrB,
			 int *cnsrC, int *cnsrD, int *cnsrE, int *cnsrF,
			 int *cnsrG, int offset)
  #endif
#endif
{
  int i, ioff, idRA;
  double dx, dy, dz, dRA;
#if VGTE == 1
  double metric;
#endif
#if MULTIPLE >= 2
  int idRB;
  double dRB;
#endif
#if MULTIPLE >= 4
  int idRC, idRD;
  double dRC, dRD;
#endif
#if MULTIPLE == 7
  int idRE, idRF, idRG;
  double dRE, dRF, dRG;
#endif
  atomc *cmapS;

  cmapS = C->map[sctrS];
  ioff = offset;
  for (i = 0; i < C->nr[sctrS]; i++) {

    // Switch metrics if necessary
#if VGTE == 1
    metric = (cmapS[i].lj < 0) ? invEEcut : invcut;
#endif

    // Atom location relative to cell midpoint
    dx = cploc[0] - cmapS[i].loc[0];
    dy = cploc[1] - cmapS[i].loc[1];
    dz = cploc[2] - cmapS[i].loc[2];

    // Displacements from planes of interest
    dRA = (fcnA[0]*dx + fcnA[1]*dy + fcnA[2]*dz) * metric;
#if MULTIPLE >= 2
    dRB = (fcnB[0]*dx + fcnB[1]*dy + fcnB[2]*dz) * metric;
#endif
#if MULTIPLE >= 4
    dRC = (fcnC[0]*dx + fcnC[1]*dy + fcnC[2]*dz) * metric;
    dRD = (fcnD[0]*dx + fcnD[1]*dy + fcnD[2]*dz) * metric;
#endif
#if MULTIPLE == 7
    dRE = (fcnE[0]*dx + fcnE[1]*dy + fcnE[2]*dz) * metric;
    dRF = (fcnF[0]*dx + fcnF[1]*dy + fcnF[2]*dz) * metric;
    dRG = (fcnG[0]*dx + fcnG[1]*dy + fcnG[2]*dz) * metric;
#endif
    if (dRA < 4.0) {
      idRA = dRA;
      cordr[cnsrA[idRA]] = ioff;
      cnsrA[idRA] += 1;
    }
#if MULTIPLE >= 2
    if (dRB < 4.0) {
      idRB = dRB;
      cordr[cnsrB[idRB]] = ioff;
      cnsrB[idRB] += 1;
    }
#endif
#if MULTIPLE >= 4
    if (dRC < 4.0) {
      idRC = dRC;
      cordr[cnsrC[idRC]] = ioff;
      cnsrC[idRC] += 1;
    }
    if (dRD < 4.0) {
      idRD = dRD;
      cordr[cnsrD[idRD]] = ioff;
      cnsrD[idRD] += 1;
    }
#endif
#if MULTIPLE >= 7
    if (dRE < 4.0) {
      idRE = dRE;
      cordr[cnsrE[idRE]] = ioff;
      cnsrE[idRE] += 1;
    }
    if (dRF < 4.0) {
      idRF = dRF;
      cordr[cnsrF[idRF]] = ioff;
      cnsrF[idRF] += 1;
    }
    if (dRG < 4.0) {
      idRG = dRG;
      cordr[cnsrG[idRG]] = ioff;
      cnsrG[idRG] += 1;
    }
#endif
    ioff++;
  }
}
