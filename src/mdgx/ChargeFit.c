#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mdgxVector.h"
#include "Matrix.h"
#include "Grid.h"
#include "CellManip.h"
#include "VirtualSites.h"
#include "Macros.h"
#include "CrdManip.h"
#include "Parse.h"
#include "Manual.h"
#include "ChargeFit.h"
#include "Trajectory.h"
#include "Topology.h"
#include "ParamFit.h"
#include "Random.h"
#include "Timings.h"

//-----------------------------------------------------------------------------
// AssignGridTopologies: when fitting electrostatic potential grids, it may be
//                       that some grids pertain to different models than the
//                       original topology.  In that case, the topologies must
//                       be read into an array and accessed with each grid in
//                       the fit.  
//
// Arguments:                                                           
//   myfit:   the fitting control data (contains matrices of names for  
//            topologies and extra points rules files, number of grids) 
//   tj:      trajectory control data                                   
//-----------------------------------------------------------------------------
static void AssignGridTopologies(fset *myfit, trajcon *tj)
{
  int i, j, k, atomsmatch, nchk;
  int* tmpidx;
  prmtop *tpi, *tpj;

  // Read all topologies
  tmpidx = (int*)malloc(myfit->ngrd * sizeof(int));
  SetIVec(tmpidx, myfit->ngrd, -1);
  myfit->TPbank = (prmtop*)malloc(myfit->ngrd * sizeof(prmtop));
  for (i = 0; i < myfit->ngrd; i++) {
    
    // Data that gets added to the topology from input file
    SetBasicTopology(&myfit->TPbank[i], myfit->tpname.map[i], myfit->eprule,
		     0, 0);
    GetPrmTop(&myfit->TPbank[i], tj, 1);
    tmpidx[i] = i;
  }

  // Correlate topologies.  Even if they have different names, two
  // topologies that have the same residues, atom names, and atom
  // types should be considered equivalent.
  for (i = 0; i < myfit->ngrd; i++) {
    if (tmpidx[i] != i) {
      continue;
    }
    tpi = &myfit->TPbank[i];
    for (j = i+1; j < myfit->ngrd; j++) {
      if (tmpidx[j] != j) {
        continue;
      }
      tpj = &myfit->TPbank[j];
      if (tpi->natom != tpj->natom || tpi->nres != tpj->nres) {
        continue;
      }
      atomsmatch = 1;
      nchk = 4*tpi->nres;
      k = 0;
      while (k < nchk && atomsmatch == 1) {
        if (tpi->ResNames[k] != tpj->ResNames[k]) {
          atomsmatch = 0;
        }
        k++;
      }
      nchk = 4*tpi->natom;
      k = 0;
      while (k < nchk && atomsmatch == 1) {
        if (tpi->AtomNames[k] != tpj->AtomNames[k]) {
          atomsmatch = 0;
        }
        k++;
      }
      k = 0;
      nchk = tpi->natom;
      while (k < nchk && atomsmatch == 1) {
        if (tpi->LJIdx[k] != tpj->LJIdx[k]) {
          atomsmatch = 0;
        }
        k++;
      }
      if (atomsmatch == 1) {
        tmpidx[j] = i;
      }
    }
  }

  // Re-read topologies, taking only the unique ones.
  for (i = 0; i < myfit->ngrd; i++) {
    FreeTopology(&myfit->TPbank[i]);
  }
  myfit->tpcount = 0;
  myfit->tpidx = (int*)malloc(myfit->ngrd * sizeof(int));
  for (i = 0; i < myfit->ngrd; i++) {
    if (tmpidx[i] == i) {
      
      // Data that gets added to the topology from input file
      k = myfit->tpcount;
      SetBasicTopology(&myfit->TPbank[k], myfit->tpname.map[i], myfit->eprule,
		       0, 0);
      GetPrmTop(&myfit->TPbank[k], tj, 1);
      for (j = i; j < myfit->ngrd; j++) {
        if (tmpidx[j] == i) {
          myfit->tpidx[j] = k;
        }
      }
      myfit->tpcount = k + 1;
    }
  }
  free(tmpidx);
  if (myfit->verbose == 1) {
    printf("mdgx >> Identified %d unique topologies out of %d conformations.\n",
           myfit->tpcount, myfit->ngrd);
  }
  
  // If there are multiple topologies in this fit, take the    
  // total charge constraints from the topologies themselves.  
  // Charges will be refit, but their sum should stay the same.
  myfit->totalq = (double*)malloc(myfit->tpcount*sizeof(double));
  for (i = 0; i < myfit->tpcount; i++) {
    myfit->totalq[i] = DSum(myfit->TPbank[i].Charges, myfit->TPbank[i].natom);
  }
}

//-----------------------------------------------------------------------------
// ColumnsForAtoms: there are a certain number of topologies involved in the
//                  fit, implying a number of unique atoms.  This routine will
//                  determine which atoms are unique, with consideration to 
//                  user-specified charge equalization constraints, and return
//                  a matrix of integers that assigns each atom to a column of
//                  the fitting matrix.
//
// Arguments:                                                           
//   myfit:      the fitting control data, contains a bank of topologies
//-----------------------------------------------------------------------------
static imat ColumnsForAtoms(fset *myfit)
{
  int i, j, k, maxatm, icol, natm, offset, maxcol;
  int* atmmask;
  int* atmresv;
  int* colscr;
  imat AC;
  coord crd;

  // Allocate the correspondence matrix
  icol = 0;
  maxatm = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    if (myfit->TPbank[i].natom > maxatm) {
      maxatm = myfit->TPbank[i].natom;
    }
    icol += myfit->TPbank[i].natom;
  }
  AC = CreateImat(myfit->tpcount, maxatm);
  maxcol = icol;

  // Determine the charge equalization restraints, which may
  // span multiple systems, and count the number of columns 
  atmresv = (int*)calloc(maxcol, sizeof(int));
  for (i = 0; i < myfit->nqeq; i++) {
    natm = 0;
    for (j = 0; j < myfit->tpcount; j++) {
      crd = CreateCoord(myfit->TPbank[j].natom);
      atmmask = ParseAmbMask(myfit->qeq[i].maskstr, &myfit->TPbank[j], &crd);
      natm += ISum(atmmask, myfit->TPbank[j].natom);
      free(atmmask);
      DestroyCoord(&crd);
    }
    myfit->qeq[i].atoms = (int*)malloc(natm*sizeof(int));
    myfit->qeq[i].natm = natm;
    natm = 0;
    offset = 0;
    for (j = 0; j < myfit->tpcount; j++) {
      crd = CreateCoord(myfit->TPbank[j].natom);
      atmmask = ParseAmbMask(myfit->qeq[i].maskstr, &myfit->TPbank[j], &crd);
      for (k = 0; k < myfit->TPbank[j].natom; k++) {
        if (atmmask[k] == 1) {
          myfit->qeq[i].atoms[natm] = k + offset;
          if (atmresv[k+offset] == 1) {
            printf("ColumnsForAtoms >> Error.  Overlapping charge "
                   "equalization constraints\nColumnsForAtoms >> detected.  "
                   "Terminating on restraint %s.\n", myfit->qeq[i].maskstr);
            exit(1);
          }
          atmresv[k+offset] = 1;
          natm++;
        }
      }
      offset += myfit->TPbank[j].natom;
      free(atmmask);
      DestroyCoord(&crd);
    }

    // The total number of columns is decreased by the
    // number of atoms in this restraint minus one    
    if (natm >= 2) {
      icol -= natm-1;
    }
  }
  free(atmresv);

  // Fill out the correspondence array
  natm = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    for (j = 0; j < myfit->TPbank[i].natom; j++) {
      AC.map[i][j] = natm;
      natm++;
    }
  }
  colscr = CountUp(maxcol);
  for (i = 0; i < myfit->nqeq; i++) {
    for (j = 1; j < myfit->qeq[i].natm; j++) {
      colscr[myfit->qeq[i].atoms[j]] = myfit->qeq[i].atoms[0];
    }
  }

  // Adjust the correspondence array to be sequential
  k = 0;
  atmmask = (int*)malloc(maxcol*sizeof(int));
  SetIVec(atmmask, maxcol, 0);
  for (i = 0; i < maxcol; i++) {
    if (atmmask[i] == 1) {
      continue;
    }
    offset = colscr[i];
    for (j = i; j < maxcol; j++) {
      if (atmmask[j] == 0 && colscr[j] == offset) {
        colscr[j] = k;
        atmmask[j] = 1;
      }
    }
    k++;
  }
  free(atmmask);

  // Atom correspondence back into the matrix
  k = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    for (j = 0; j < myfit->TPbank[i].natom; j++) {
      AC.map[i][j] = colscr[k];
      k++;
    }
  }

  // Free allocated memory
  free(colscr);

  myfit->q2fit = icol;
  return AC;
}

//-----------------------------------------------------------------------------
// CheckDegeneracy: procedure for checking through all fitting matrix columns
//                  to see whether any two columns pertain to atoms with
//                  overlapping names and types.
//
// Arguments:
//   myfit:      the fitting control data, contains a bank of topologies
//   Atm2Col:    a map of which atoms' partial charges will be solved in each
//               column of the fitting matrix
//-----------------------------------------------------------------------------
static void CheckDegeneracy(fset *myfit, imat *Atm2Col)
{
  int i, j, k, m, ires, jres, irna, jrna, nitype, njtype, flagged;
  double ilja, iljb, jlja, jljb;
  char *ianm, *janm, *irnm, *jrnm;
  imat dgs;
  prmtop *tpi, *tpj;

  // Map for any flagged atoms
  j = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    if (myfit->TPbank[i].natom > j) {
      j = myfit->TPbank[i].natom;
    }
  }
  dgs = CreateImat(myfit->tpcount, j);
  for (i = 0; i < myfit->tpcount; i++) {
    tpi = &myfit->TPbank[i];
    nitype = tpi->ntypes;
    for (k = 0; k < tpi->natom; k++) {
      if (dgs.map[i][k] == 1) {
	continue;
      }
      
      // Identify this atom's residue, name, and other features
      ires = LocateResID(tpi, k, 0, tpi->nres);
      irna = tpi->ResLims[ires + 1] - tpi->ResLims[ires];
      ianm = &tpi->AtomNames[4*k];
      irnm = &tpi->ResNames[4*ires];
      if (tpi->LJIdx[k] >= 0) {
        ilja = tpi->LJA[tpi->NBParmIdx[(nitype + 1)*tpi->LJIdx[k]]];
        iljb = tpi->LJB[tpi->NBParmIdx[(nitype + 1)*tpi->LJIdx[k]]];
      }
      else {
        ilja = 0.0;
        iljb = 0.0;
      }
      flagged = 0;
      
      // Loop over this topology and all those following it
      for (j = i; j < myfit->tpcount; j++) {
        tpj = &myfit->TPbank[j];
        njtype = tpj->ntypes;
        for (m = 0; m < myfit->TPbank[j].natom; m++) {
          if (dgs.map[j][m] == 1) {
            continue;
          }

          // Identify this atom's residue, name, and other features
          jres = LocateResID(tpj, m, 0, tpj->nres);
          jrna = tpj->ResLims[jres + 1] - tpj->ResLims[jres];
          janm = &tpj->AtomNames[4*m];
          jrnm = &tpj->ResNames[4*jres];
          if (tpj->LJIdx[m] >= 0) {
            jlja = tpj->LJA[tpj->NBParmIdx[(njtype + 1)*tpj->LJIdx[m]]];
            jljb = tpj->LJB[tpj->NBParmIdx[(njtype + 1)*tpj->LJIdx[m]]];
          }
          else {
            jlja = 0.0;
            jljb = 0.0;
          }
          
          // Check for correspondence
          if (ires != jres && irna == jrna &&
              str4cmp(ianm, janm) == 0 && str4cmp(irnm, jrnm) == 0 &&
	      fabs(ilja - jlja) <= 1.0e-4 * fabs(ilja + 1.0e-4) &&
              fabs(iljb - jljb) <= 1.0e-4 * fabs(iljb + 1.0e-4)) {
            if (flagged == 0) {
              printf("mdgx >> Warning.  Atom %4.4s %2d %4.4s of topology %s\n"
                     "mdgx >>   and atom     %4.4s %2d %4.4s of topology %s\n"
                     "mdgx >>   are fitted independently, but may need to be "
                     "linked.\n", irnm, ires + 1, ianm, tpi->source, jrnm,
		     jres + 1, janm, tpj->source);
	      flagged++;
	      dgs.map[i][k] = 1;
	      dgs.map[j][m] = 1;
            }
	    else {
	      flagged++;
	      dgs.map[i][k] = 1;
	      dgs.map[j][m] = 1;
	    }
          }
        }
      }
    }
  }

  // Free allocated memory
  DestroyImat(&dgs);
}

//-----------------------------------------------------------------------------
// ReadEPotGrid: read an electrostatic potential grid as output by the Gaussian
//               cubegen utility.  Only formatted cubegen outputs are read.
//
// Arguments:                                                           
//   fname:     the name of the cubegen output file                     
//   tp:        the name of the topology used to interpret the molecular
//              coordiantes at the head of the cubegen output file      
//   crd:       a coord struct to hold the coordinates at the head of   
//              the cubegen output file                                 
//-----------------------------------------------------------------------------
fbook ReadEPotGrid(char* fname, prmtop *tp, coord *crd)
{
  int h, i, j, k, iinit, jinit, kinit, headerfound;
  int slen, nnum, irexp, formgss;
  double rsign, rsum;
  double orig[3];
  double* nucchg;
  double* eptab;
  float *ftmp;
  char line[MAXLINE];
  char *ctmp;
  FILE *inp;
  dmat rL;
  cmat lwords;
  fbook Ue;

  // Lookup table for exponents with strictly formatted input
  eptab = (double*)malloc(199*sizeof(double));
  eptab[99] = 1.0;
  for (i = 1; i <= 99; i++) {
    eptab[99+i] = 10.0*eptab[99+i-1];
    eptab[99-i] = 0.1*eptab[99-i+1];
  }

  // Open the cubegen file
  if ((inp = fopen(fname, "r")) == NULL) {
    printf("ReadEPotGrid >> Error.  File %s not found.\n", fname);
    exit(1);
  }

  // Process the cubegen output
  headerfound = 1;
  while (headerfound == 1) {
    fgets(line, 128, inp);
    j = strlen(line);
    headerfound = 0;
    for (i = 0; i < strlen(line); i++) {
      if ((line[i] < '0' || line[i] > '9') && line[i] != '.' &&
          line[i] != '-' && line[i] != ' ' && line[i] != '\n') {
        headerfound = 1;
        break;
      }
    }
  }
  sscanf(line, "%d%lf%lf%lf", &crd->natom, &orig[0], &orig[1], &orig[2]);
  if (crd->natom != tp->natom - tp->nxtrapt) {
    printf("ReadEPotGrid >> Error.  %d real atoms present in %s,"
           "\nReadEPotGrid >> %d in the associated topology %s.\n", crd->natom,
           fname, tp->natom, tp->source);
    exit(1);
  }
  rL = CreateDmat(3, 3, 0);
  fgets(line, MAXLINE, inp);
  sscanf(line, "%d%lf%lf%lf", &i, &rL.data[0], &rL.data[1], &rL.data[2]);
  fgets(line, MAXLINE, inp);
  sscanf(line, "%d%lf%lf%lf", &j, &rL.data[3], &rL.data[4], &rL.data[5]);
  fgets(line, MAXLINE, inp);
  sscanf(line, "%d%lf%lf%lf", &k, &rL.data[6], &rL.data[7], &rL.data[8]);
  Ue = CreateFbook(i, j, k, &rL, orig);
  DestroyDmat(&rL);
  nucchg = (double*)malloc(crd->natom*sizeof(double));
  for (i = 0; i < crd->natom; i++) {
    fgets(line, MAXLINE, inp);
    sscanf(line, "%d%lf%lf%lf%lf", &j, &nucchg[i], &crd->loc[3*i],
           &crd->loc[3*i+1], &crd->loc[3*i+2]);
  }
  if (tp->nxtrapt > 0) {

    // Treat all extra points as "added" particles and assume that the QM
    // electrostatic grid only contains information about real atoms.
    j = 0;
    k = tp->norigatom;
    for (i = 0; i < tp->natom; i++) {
      j += (tp->Masses[i] >= 1.0e-8);
    }
    tp->norigatom = j;
    if (tp->EPInserted == 0) {
      tp->OldAtomNum = (int*)malloc(tp->natom * sizeof(int));
    }
    j = 0;
    for (i = 0; i < tp->natom; i++) {
      if (tp->Masses[i] < 1.0e-8) {
        tp->OldAtomNum[i] = -1;
      }
      else {
        tp->OldAtomNum[i] = j;
        j++;
      }
    }
    ExtendCoordinates(crd, tp);
    
    // Set the topology back to its original state
    if (tp->EPInserted == 0) {
      free(tp->OldAtomNum);
    }
    tp->norigatom = k;
  }
  if (crd->natom != tp->natom) {
    printf("ReadEPotGrid >> Error.  Atom count on grid %d, %d in topology "
           "%s.\n", crd->natom, tp->natom, tp->source);
    exit(1);
  }

  // Scan the grid data into memory, with routines to
  // boost performance when reading a standard format
  i = 0;
  j = 0;
  k = 0;
  formgss = 1;
  ftmp = Ue.map[i][j];
  while (formgss == 1 && fgets(line, MAXLINE, inp) != NULL && i < Ue.pag) {
    nnum = 0;
    slen = strlen(line);
    for (h = 0; h < slen; h++) {
      if (line[h] == '.') {
        nnum++;
      }
    }
    for (h = 0; h < nnum; h++) {
      if (line[13*h] != ' ' || line[13*h+3] != '.' ||
          !(line[13*h+9] == 'E' || line[13*h+9] == 'e')) {
        formgss = 0;
      }
    }
    if (formgss == 0) {
      printf("Format broken!\n");
      lwords = ParseWords(line);
      for (h = 0; h < lwords.row; h++) {
        ftmp[k] = atof(lwords.map[h]);
        k++;
        if (k == Ue.col) {
          k = 0;
          j++;
          if (j == Ue.row) {
            j = 0;
            i++;
          }
          if (i < Ue.pag) {
            ftmp = Ue.map[i][j];
          }
        }
      }
      DestroyCmat(&lwords);
      continue;
    }

    // The format appears to be upheld; read the numbers
    for (h = 0; h < nnum; h++) {
      ctmp = &line[13*h];
      rsign = (ctmp[1] == '-') ? -1.0 : 1.0;
      rsum = (ctmp[2] - '0') + 0.1*(ctmp[4] - '0') + 0.01*(ctmp[5] - '0') +
        0.001*(ctmp[6] - '0') + 0.0001*(ctmp[7] - '0') +
        0.00001*(ctmp[8] - '0');
      irexp = 10*(ctmp[11] - '0') + ctmp[12] - '0';
      if (ctmp[10] == '-') {
        irexp = -irexp;
      }
      ftmp[k] = rsign*rsum*eptab[irexp+99];
      k++;
      if (k == Ue.col) {
        k = 0;
        j++;
        if (j == Ue.row) {
          j = 0;
          i++;
        }
        if (i < Ue.pag) {
          ftmp = Ue.map[i][j];
        }
      }
    }
  }

  // Bail out of strict formatting and pick up where we left off
  if (formgss == 0) {
    iinit = i;
    jinit = j;
    kinit = k;
    for (i = iinit; i < Ue.pag; i++) {
      for (j = jinit; j < Ue.row; j++) {
        ftmp = Ue.map[i][j];
        for (k = kinit; k < Ue.col; k++) {
          fscanf(inp, "%f", &ftmp[k]);
          iinit = 0;
          jinit = 0;
          kinit = 0;
        }
      }
    }
  }
  fclose(inp);

  // Free allocated memory
  free(nucchg);
  free(eptab);

  const int nelem = Ue.pag*Ue.row*Ue.col;
  for (j = 0; j < nelem; j++) {
    Ue.data[j] *= H2KCAL;
  }
  for (j = 0; j < 9; j++) {
    Ue.L.data[j] *= B2ANG;
  }
  for (j = 0; j < 3; j++) {
    Ue.orig[j] *= B2ANG;
  }
  for (j = 0; j < 3*crd->natom; j++) {
    crd->loc[j] *= B2ANG;
  }
  MolPlaceXpt(crd, tp);
  
  return Ue;
}

//-----------------------------------------------------------------------------
// PrepUPot: prepares an electrostatic potential grid read from a file for RESP
//           fitting.  Units of the grid potential, scale, and molecular
//           coordinates are fitted.                          
//
// Arguments:                                                           
//   UPot:    the electrostatic potential grid                          
//   crd:     molecular coordinates associated with UPot                
//   tp:      the topology                                              
//   Uflag:   grid of indices describing the accessibility and availability of
//            points in UPot                            
//   Uminr:   grid recording the minimum range of grid points to an atom
//   myfit:   RESP fitting control data                                 
//-----------------------------------------------------------------------------
static void PrepUPot(fbook *UPot, coord *crd, prmtop *tp, cbook *Uflag,
                     dbook *Uminr, fset *myfit)
{
  int h, i, j, k, ib, jb, kb, m, ljt;
  int mini, minj, mink, maxi, maxj, maxk, buffi, buffj, buffk;
  double r2, invr2, invr4, sig, eps, Afac, Bfac, dx, dy, dz, di, dj, dk;
  double atmloc[3], ptloc[3], orig[3], cdepth[3];
  double *dtmp, *dtm2p, *Ldata;
  char *ctmp;
  dmat *Utbl;
  dbook Ulj;

  Ulj = CreateDbook(UPot->pag, UPot->row, UPot->col, 0);

  // Compute the Lennard-Jones energy of each grid point
  Utbl = &tp->LJutab;
  Ldata = UPot->L.data;
  orig[0] = UPot->orig[0];
  orig[1] = UPot->orig[1];
  orig[2] = UPot->orig[2];
  SetDVec(Uminr->data, Uminr->pag * Uminr->row * Uminr->col, 1.0e12);
  for (h = 0; h < tp->natom; h++) {

    // Do not include extra points in this analysis.  They could
    // alter the way data is selected, leading to biased results.
    if (tp->Masses[h] < 1.0e-8) {
      continue;
    }
    for (m = 0; m < 3; m++) {
      atmloc[m] = crd->loc[3*h+m];
    }
    ljt = tp->LJIdx[h];
    if (ljt < 0) {
      sig = 0.0;
      eps = 0.0;
    }
    else {
      sig = pow(-1.0*Utbl->map[ljt][2*ljt]/Utbl->map[ljt][2*ljt+1], 1.0/6.0);
      eps = -0.25 * Utbl->map[ljt][2*ljt+1] / pow(sig, 6.0);
      sig = 0.5*(sig + myfit->psig);
      eps = sqrt(eps*myfit->peps);
      Afac = 4.0*eps*pow(sig, 12.0);
      Bfac = 4.0*eps*pow(sig, 6.0);
    }
    for (i = 0; i < Ulj.pag; i++) {
      di = i;
      for (j = 0; j < Ulj.row; j++) {
        dtmp = Ulj.map[i][j];
        dtm2p = Uminr->map[i][j];
        dj = j;
        ptloc[0] = Ldata[0]*di + Ldata[1]*dj + orig[0] - atmloc[0];
        ptloc[1] = Ldata[3]*di + Ldata[4]*dj + orig[1] - atmloc[1];
        ptloc[2] = Ldata[6]*di + Ldata[7]*dj + orig[2] - atmloc[2];
        for (k = 0; k < Ulj.col; k++) {
          r2 = ptloc[0]*ptloc[0] + ptloc[1]*ptloc[1] + ptloc[2]*ptloc[2];
          ptloc[0] += Ldata[2];
          ptloc[1] += Ldata[5];
          ptloc[2] += Ldata[8];
          invr2 = 1.0/(r2);
          invr4 = invr2*invr2;
          dtmp[k] += Afac*invr4*invr4*invr4 - Bfac*invr4*invr2;
          if (r2 < dtm2p[k]) {
            dtm2p[k] = r2;
          }
        }
      }
    }
  }

  // Compute the grid buffer region for accessibility
  HessianNorms(&UPot->L, cdepth);
  buffi = myfit->prbarm/cdepth[0] + 1;
  buffj = myfit->prbarm/cdepth[1] + 1;
  buffk = myfit->prbarm/cdepth[2] + 1;
  for (i = 0; i < Ulj.pag; i++) {
    for (j = 0; j < Ulj.row; j++) {
      dtmp = Ulj.map[i][j];
      ctmp = Uflag->map[i][j];
      for (k = 0; k < Ulj.col; k++) {
        
        // If this point is accessible to the probe,
        // declare it so and mark all surrounding   
        // points accessible if they are within the 
        // probe arm's reach.                       
        if (dtmp[k] < myfit->stericlim) {
          ctmp[k] = 1;
          continue;
        }

        // If we're still here, this point is not 
        // accessible to probe but may nonetheless
        // be accessible to the probe arm.        
        mini = MAX(i-buffi, 0);
        minj = MAX(j-buffj, 0);
        mink = MAX(k-buffk, 0);
        maxi = MIN(i+buffi, Ulj.pag-1);
        maxj = MIN(j+buffj, Ulj.row-1);
        maxk = MIN(k+buffk, Ulj.col-1);
        for (ib = mini; ib <= maxi; ib++) {
          di = ib - i;
          for (jb = minj; jb <= maxj; jb++) {
            dj = jb - j;
            dtm2p = Ulj.map[ib][jb];
            dk = mink-k-1;
            ptloc[0] = di*Ldata[0] + dj*Ldata[1] + dk*Ldata[2];
            ptloc[1] = di*Ldata[3] + dj*Ldata[4] + dk*Ldata[5];
            ptloc[2] = di*Ldata[6] + dj*Ldata[7] + dk*Ldata[8];
            for (kb = mink; kb <= maxk; kb++) {
              ptloc[0] += Ldata[2];
              ptloc[1] += Ldata[5];
              ptloc[2] += Ldata[8];

              // In order for this point to be worth investigating,
              // it has to be possible for a probe to be situated  
              // on this point with an acceptable steric energy and
              // then reach its arm out to the point that we've    
              // thus far declared inaccessible.                   
              if (dtm2p[kb] > myfit->stericlim) {
                continue;
              }
              if (ptloc[0]*ptloc[0] + ptloc[1]*ptloc[1] + ptloc[2]*ptloc[2] <
                  myfit->prbarm) {

                // Mark the point accessible and 
                // bail out of these nested loops
                ctmp[k] = 1;
                ib = maxi+1;
                jb = maxj+1;
                kb = maxk+1;
              }
            }
          }
        }
      }
    }
  }

  // Free allocated memory
  DestroyDbook(&Ulj);
}

//-----------------------------------------------------------------------------
// MaskAllAtoms: restraints may span multiple systems.  This function examines
//               the restraint mask and applies it to each system individually
//               to determine how many unique atoms it involves.
//
// Arguments:                                                           
//   myfit:       the fitting control data                              
//   qmask:       the mask of all charges (length the number of columns 
//                in the fitting matrix)                                
//   maskstr:     the mask string (ambmask format)                      
//   Atm2Col:     matrix describing the mapping of atoms in each system 
//                to columns of the fitting constraint matrix           
//   stack:       flag to record counts of how many atoms share each    
//                unique fitted charge state; set to zero for boolean   
//                indications as to whether each fitted charge is       
//                present in the mask, 1 to record counts               
//   spec:        do the mask only for a specific topology; set to -1 to
//                make the mask for all topologies                      
//-----------------------------------------------------------------------------
static int MaskAllAtoms(fset *myfit, int* qmask, char* maskstr,
                        imat *Atm2Col, int stack, int spec)
{
  int i, j, natm;
  int *itmp;
  int* atmmask;
  prmtop *tp;
  coord crd;

  SetIVec(qmask, myfit->q2fit, 0);
  for (i = 0; i < myfit->tpcount; i++) {
    if (spec >= 0 && i != spec) {
      continue;
    }
    tp = &myfit->TPbank[i];
    crd = CreateCoord(tp->natom);
    atmmask = ParseAmbMask(maskstr, tp, &crd);
    itmp = Atm2Col->map[i];
    for (j = 0; j < tp->natom; j++) {
      if (atmmask[j] == 1) {
        qmask[itmp[j]] = (stack == 1) ? qmask[itmp[j]] + 1 : 1;
      }
    }
    free(atmmask);
    DestroyCoord(&crd);
  }
  natm = ISum(qmask, myfit->q2fit);

  return natm;
}

//-----------------------------------------------------------------------------
// MakeCnstMatrix: make a constraint matrix for charge fitting based on a
//                 series of ambmask strings and other data specified by the
//                 user.  This constraint matrix has n+2 columns for fitting n
//                 charges; column n+1 contains the target, to go on the right
//                 hand side of the linear least squares problem (the "b"
//                 vector in Ax=b); column n+2 contains the constraint type as
//                 follows:  
//                 0 = total charge constraint (always present),        
//                 1 = minimize q constraint, 2 = sum q constraint,     
//                 3 = tether q constraint, 4 = solvated q constraint   
//                 5 = vacuum q constraint.                             
//                                                                      
// Arguments:                                                           
//   myfit:    the fitting control data                                 
//   Atm2Col:  matrix describing the mapping of atoms in each system to 
//             columns of the fitting constraint matrix                 
//   fitmat:   the augmented fitting matrix [ A b ]                     
//-----------------------------------------------------------------------------
static dmat MakeCnstMatrix(fset *myfit, imat *Atm2Col, dmat *fitmat)
{
  int i, j, k, m, nrst, rstcon, natm, tpkey, maxatom;
  int* allqmask;
  int* fixqmask;
  int* tmpqmask;
  int* nsamp;
  double rstfac, totq, ninst;
  double *dtmp;
  dmat cnstmat;
  
  // The obligatory total charge restraints, one per system
  cnstmat = CreateDmat(myfit->tpcount, myfit->q2fit+2, 0);
  maxatom = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    maxatom += myfit->TPbank[i].natom;
    for (j = 0; j < myfit->TPbank[i].natom; j++) {
      cnstmat.map[i][Atm2Col->map[i][j]] += 1.0e7;
    }
    cnstmat.map[i][myfit->q2fit] = (1.0e7)*myfit->totalq[i];
  }
  nrst = myfit->tpcount;
  rstcon = myfit->tpcount;

  // Count the number of times each atom is sampled in the
  // fitting matrix, the number of equations involved, in 
  // order to properly scale the restraint weights.       
  nsamp = (int*)calloc(myfit->q2fit, sizeof(int));
  for (i = 0; i < fitmat->row; i++) {
    dtmp = fitmat->map[i];
    for (j = 0; j < myfit->q2fit; j++) {
      if (fabs(dtmp[j]) > 1.0e-8) {
        nsamp[j] += 1;
      }
    }
  }

  // Store the sampling for later analysis (and perhaps use)
  myfit->nsamp = nsamp;

  // Charge minimization restraints
  allqmask = (int*)malloc(maxatom*sizeof(int));
  for (i = 0; i < myfit->nqmin; i++) {

    // This restraint may span multiple systems.  We must   
    // determine how many unique atoms are being restrained. 
    natm = MaskAllAtoms(myfit, allqmask, myfit->qmin[i].maskstr,
                        Atm2Col, 0, -1);
    if (natm == 0) {
      continue;
    }

    // Add new restraints for each specified charge
    nrst += natm;
    cnstmat = ReallocDmat(&cnstmat, nrst, myfit->q2fit+2);
    for (j = 0; j < myfit->q2fit; j++) {
      if (allqmask[j] == 1) {
        dtmp = cnstmat.map[rstcon];
        dtmp[j] = myfit->qminwt*nsamp[j];
        dtmp[myfit->q2fit+1] = 1.01;
        rstcon++;
      }
    }
  }

  // Charge group sum restraints
  tmpqmask = (int*)malloc(maxatom*sizeof(int));
  for (i = 0; i < myfit->nqsum; i++) {
    
    // Check to see that this charge group
    // does exist in some system.                           
    natm = 0;
    for (j = 0; j < myfit->tpcount; j++) {
      natm = MaskAllAtoms(myfit, allqmask, myfit->qsum[i].maskstr,
                          Atm2Col, 1, j);
      if (natm > 0) {
        tpkey = j;
        break;
      }
    }
    if (natm == 0) {
      continue;
    }

    // It is now known that this charge group exists in
    // one or more of the systems, in whole or in part.
    // Now we must determine whether the charge group  
    // exists entirely or not at all in every system.  
    // If this condition is not met, it would result in
    // very strange and contradictory restraints, so   
    // stop the program if a bad case is detected.     
    for (j = 0; j < myfit->tpcount; j++) {
      natm = MaskAllAtoms(myfit, tmpqmask, myfit->qsum[i].maskstr,
                          Atm2Col, 1, j);
      if (ISum(tmpqmask, natm) == 0) {
        continue;
      }
      for (k = 0; k < myfit->q2fit; k++) {
        if (tmpqmask[k] != allqmask[k]) {
          printf("MakeCnstMatrix >> Error.  Mismatch in masks generated for "
                 "topologies\nMakeCnstMatrix >> %s and %s.\n",
                 myfit->TPbank[tpkey].source, myfit->TPbank[j].source);
          printf("MakeCnstMatrix >> This has occurred in charge sum restraint "
                 "mask processing\nMakeCnstMatrix >> and probably indicates a "
                 "missing charge equivalence between\nMakeCnstMatrix >> atoms "
                 "of the same residue existing in multiple systems.\n"
                 "MakeCnstMatrix >> For example, if you have ACE-ALA-NME and "
                 "ACE-ARG-NME in your\nMakeCnstMatrix >> data set and "
                 "specified that residues ACE and NME have\nMakeCnstMatrix >> "
                 "zero net charge, you need to specify that every instance "
                 "of\nMakeCnstMatrix >> each atom in ACE or NME be taken as "
                 "equivalent across all\nMakeCnstMatrix >> systems.  "
                 "Specifying sumq ':ACE' 0.0 and sumq ':NME' 0.0\n"
                 "MakeCnstMatrix >> would need to be accompanied by equalq "
                 "directives for each\nMakeCnstMatrix >> atom in ACE and NME."
                 "\n");
          exit(1);
        }
      }
    }

    // Add the new restraint for this group
    nrst += 1;
    cnstmat = ReallocDmat(&cnstmat, nrst, myfit->q2fit+2);
    dtmp = cnstmat.map[rstcon];
    for (j = 0; j < myfit->q2fit; j++) {
      if (allqmask[j] > 0) {
        dtmp[j] = allqmask[j] * 1.0e7;
      }
    }
    dtmp[myfit->q2fit] = 1.0e7 * myfit->qsum[i].target;
    dtmp[myfit->q2fit+1] = 2.01;
    rstcon++;
  }

  // Charge tethering
  if (myfit->model == 0 && myfit->tether == 1) {
    cnstmat = ReallocDmat(&cnstmat, nrst+myfit->q2fit, myfit->q2fit+2);
    for (i = 0; i < myfit->q2fit; i++) {
      rstfac = myfit->qtthwt*nsamp[i];
      cnstmat.map[rstcon][i] = rstfac;
      ninst = 0.0;
      totq = 0.0;
      for (j = 0; j < myfit->tpcount; j++) {
        for (k = 0; k < myfit->TPbank[j].natom; k++) {
          if (Atm2Col->map[j][k] == i) {
            totq += myfit->TPbank[j].Charges[k];
            ninst += 1.0;
          }
        }
      }
      cnstmat.map[rstcon][myfit->q2fit] = rstfac * totq / ninst;
      cnstmat.map[rstcon][myfit->q2fit+1] = 3.01;
      rstcon++;
    }
    nrst += myfit->q2fit;
  }

  // Charge fixing
  fixqmask = (int*)malloc(maxatom*sizeof(int));
  for (i = 0; i < myfit->nqfix; i++) {

    // This restraint may span multiple systems.  We must   
    // determine how many unique atoms are being restrained.
    natm = MaskAllAtoms(myfit, fixqmask, myfit->qfixsolv[i].maskstr,
                        Atm2Col, 0, -1);
    if (natm == 0) {
      continue;
    }

    // Add new restraints for each specified charge
    nrst += 2;
    cnstmat = ReallocDmat(&cnstmat, nrst, myfit->q2fit+2);
    for (j = 0; j < myfit->q2fit; j++) {
      if (fixqmask[j] == 1) {
        cnstmat.map[rstcon][j] = fixqmask[j] * 1.0e7;
        cnstmat.map[rstcon][myfit->q2fit] = 1.0e7 * myfit->qfixsolv[i].target *
                                            natm;
        cnstmat.map[rstcon][myfit->q2fit+1] = 4.01;
        cnstmat.map[rstcon+1][j] = fixqmask[j] * 1.0e7;
        cnstmat.map[rstcon+1][myfit->q2fit] = 1.0e7 *
                                              myfit->qfixvacu[i].target * natm;
        cnstmat.map[rstcon+1][myfit->q2fit+1] = 5.01;
      }
    }
    rstcon+=2;
  }
  if (myfit->verbose == 1) {
    printf("mdgx >> Constraint matrix composed.\n");
  }

  // Free allocated memory
  free(allqmask);
  free(fixqmask);
  free(tmpqmask);

  return cnstmat;
}

//-----------------------------------------------------------------------------
// SortPtDistance: function called by quicksort for comparing minimum distances
//                 between a grid point and the molecule.     
//                                                                      
// Arguments:                                                           
//   pt[A,B]:    the atomc structs                                      
//-----------------------------------------------------------------------------
static int SortPtDistance(const void *ptA, const void *ptB)
{
  double disA = ((fitpt*)ptA)[0].minr;
  double disB = ((fitpt*)ptB)[0].minr;

  if (disA < disB) {
    return -1;
  }
  else if (disA > disB) {
    return 1;
  }
  else {
    return 0;
  }
}

//-----------------------------------------------------------------------------
// ContributeFitPt: contribute this fitting point to the matrix of candidates.
//
// Arguments:                                                           
//   myfit:    the fitting set                                          
//   nsys:     the system number                                        
//   A:        the fitting matrix                                       
//   nln:      the line number                                          
//   grdpt:    pointer to grid point catalog                            
//   crd:      coordinates                                              
//   tp:       topology                                                 
//   fitline:  the next line of the growing fitting matrix to be written
//   fitrec:   record of the fitting points, kept for user interest
//-----------------------------------------------------------------------------
static void ContributeFitPt(fset *myfit, int nsys, dmat *A, int nln,
                            fitpt *grdpt, coord *crd, fbook *UPotA,
                            fbook *UPotB, imat *Atm2Col, FILE* fitrec)
{
  int i, ncol;
  int *colmap;
  double dx, dy, dz, r, wt, minr;
  double dijk[3];
  double *fitline;

  // Unpack structures and perform simple catalogging
  wt = myfit->wt[nsys];
  colmap = Atm2Col->map[myfit->tpidx[nsys]];
  ncol = myfit->q2fit;
  fitline = A->map[nln];
  myfit->FPtOrigins[nln] = myfit->tpidx[nsys];

  // Compute the sum of kq/r interactions
  dijk[0] = grdpt->ix;
  dijk[1] = grdpt->iy;
  dijk[2] = grdpt->iz;
  RotateCrd(dijk, 1, UPotA->L);
  dijk[0] += UPotA->orig[0];
  dijk[1] += UPotA->orig[1];
  dijk[2] += UPotA->orig[2];
  minr = 1.0e30;
  for (i = 0; i < crd->natom; i++) {
    dx = crd->loc[3*i] - dijk[0];
    dy = crd->loc[3*i+1] - dijk[1];
    dz = crd->loc[3*i+2] - dijk[2];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    fitline[colmap[i]] += wt*BIOQ/r;
    if (fitrec != NULL && r < minr) {
      minr = r;
    }
  }
  fitline[ncol] = wt*0.5*(UPotA->map[grdpt->ix][grdpt->iy][grdpt->iz] +
                          UPotB->map[grdpt->ix][grdpt->iy][grdpt->iz]);
  fitline[ncol+1] = wt*UPotA->map[grdpt->ix][grdpt->iy][grdpt->iz];
  if (fitrec != NULL) {
    fprintf(fitrec, "  %9.4lf %9.4lf %9.4lf   %9.4lf %9.4lf  %9.4lf\n",
            dijk[0], dijk[1], dijk[2], minr, fitline[ncol], fitline[ncol+1]);
  }
}

//-----------------------------------------------------------------------------
// FlagFitPt: flag all points within a short cutoff around a fitting point that
//            was just committed to the matrix of candidate fitting points.
//            These flags, which prevent points from being used in the fit,
//            make it so that extremely closely positioned points do not both
//            contribute to the fitting matrix.
//
// Arguments:
//   grdpts:   array of candidate fitting points                        
//   ptid:     ID number of the fitting point within the grdpts catalog 
//   npt:      the total number of points catalogged in grdpts          
//   UPot:     the electrostatic potential grid                         
//   UIdx:     indexes for whether the grid points are flagged or not   
//   catIdx:   maps the points in UPot/UIdx to their entries in grdpts  
//-----------------------------------------------------------------------------
static void FlagFitPt(fset *myfit, fitpt* grdpts, int ptid, fbook *UPot,
                      cbook *Uflag, ibook *UIdx, int buffi, int buffj,
                      int buffk)
{
  int i, j, k, imin, imax, jmin, jmax, kmin, kmax, homei, homej, homek;
  double r2, dx, dy, dz;
  double dijk[3];
  double *Ldata;
  int *itmp;
  char *ctmp;

  homei = grdpts[ptid].ix;
  homej = grdpts[ptid].iy;
  homek = grdpts[ptid].iz;
  imin = MAX(0, homei - buffi);
  imax = MIN(UPot->pag, homei + buffi);
  jmin = MAX(0, homej - buffj);
  jmax = MIN(UPot->pag, homej + buffj);
  kmin = MAX(0, homek - buffk);
  kmax = MIN(UPot->pag, homek + buffk);
  Ldata = UPot->L.data;
  for (i = imin; i < imax; i++) {
    for (j = jmin; j < jmax; j++) {
      ctmp = Uflag->map[i][j];
      itmp = UIdx->map[i][j];
      for (k = kmin; k < kmax; k++) {
        if (ctmp[k] == 0) {
          continue;
        }
        dx = i - homei;
        dy = j - homej;
        dz = k - homek;
        dijk[0] = Ldata[0]*dx + Ldata[1]*dy + Ldata[2]*dz;
        dijk[1] = Ldata[3]*dx + Ldata[4]*dy + Ldata[5]*dz;
        dijk[2] = Ldata[6]*dx + Ldata[7]*dy + Ldata[8]*dz;
        r2 = dijk[0]*dijk[0] + dijk[1]*dijk[1] + dijk[2]*dijk[2];
        if (r2 < myfit->flimit) {
          ctmp[k] = 0;
          grdpts[itmp[k]].flagged = 1;
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// HistogramFitPt: place a fitting point in a histogram based on its distance
//                 from the molecule.                          
//                                                                      
// Arguments:                                                           
//   myfit:    the fitting control structure (contains the pre-allocated
//             histogram)                                               
//   gpt:      the fitting point                                        
//-----------------------------------------------------------------------------
static void HistogramFitPt(fset *myfit, fitpt *gpt)
{
  int ir;

  ir = gpt->minr/myfit->fhistbin;
  myfit->fitpthist[ir] += 1;
}

//-----------------------------------------------------------------------------
// WriteConformation: write a conformation of the molecule, for purposes of
//                    visualizing the charge distribution, in PDB format.
//
// Arguments:                                                           
//   h:       the number of the conformation                            
//   crd:     the coordinates of the conformation                       
//   tp:      the topology                                              
//   mfit:    the fitting set                                           
//   tj:      trajectory control data (directive to overwrite outputs)  
//-----------------------------------------------------------------------------
static void WriteConformation(int h, double* crd, prmtop *tp, fset *myfit,
                              trajcon *tj)
{
  int i, ires;
  char outname[MAXNAME];
  FILE *outp;

  // Only do this for the first instance of each system
  for (i = 0; i < h; i++) {
    if (myfit->tpidx[i] == myfit->tpidx[h]) {
      return;
    }
  }

  // Print the conformation iin PDB format
  sprintf(outname, "%s.%s", tp->source, myfit->confext);
  outp = FOpenSafe(outname, tj->OverwriteOutput);
  for (i = 0; i < tp->natom; i++) {
    ires = LocateResID(tp, i, 0, tp->nres);
    fprintf(outp, "ATOM %6d %.4s %.4s%c%4d    %8.3lf%8.3lf%8.3lf\n",
            i, &tp->AtomNames[4*i], &tp->ResNames[4*ires], 'A', 1,
            crd[3*i], crd[3*i+1], crd[3*i+2]);
  }
  fprintf(outp, "END\n");
  fclose(outp);
}

//-----------------------------------------------------------------------------
// SelectFitPoint: introduce a criterion to limit the number of points far from
//                 the molecular surface which enter the fit.  The criterion is
//                 that, after a certain cutoff Rc, the probability of
//                 accepting a point at a distance r from the molecular surface
//                 drops off such that the number of points accepted at r would
//                 be equal to the number at Rc if the molecule were perfectly
//                 spherical.
//
// Arguments:                                                           
//   gpt:      the grid point                                           
//   myfit:    the fitting control structure                            
//   tj:       trajectory control data (contains random number counter) 
//-----------------------------------------------------------------------------
static int SelectFitPoint(fitpt *gpt, fset *myfit, trajcon *tj)
{
  double r;

  // Accept if the point is within Rc of the molecule
  if (gpt->minr < myfit->Rc) {
    return 1;
  }

  // Accept a roughly constant number of points for
  // any given distance from the molecule beyond Rc
  r = myfit->Rc / gpt->minr;
  if (gpt->minr < myfit->Rmax && ran2(&tj->rndcon) < r*r) {
    return 1;
  }

  // Reject the point
  return 0;
}

//-----------------------------------------------------------------------------
// MakeFittingMatrix: compute the fitting matrix based on points taken from all
//                    grids, sieved by various user-specified cutoffs.
//
// Arguments:
//   myfit:    the fitting control structure                            
//   Atm2Col:  correspondence between atoms and matrix columns          
//   allcrd:   matrix to accumulate coordinates of all molecular        
//             conformations                                            
//   tj:       trajectory control information (contains random number   
//             generator counter)                                       
//   etimers:  the execution timings data                               
//-----------------------------------------------------------------------------
static dmat MakeFittingMatrix(fset *myfit, imat *Atm2Col, dmat *allcrd,
                              trajcon *tj, execon *etimers)
{
  int h, i, j, k, m, buffi, buffj, buffk, gcon, ngpt, namelen;
  int congruent, setfitpt, totalfitpt;
  long long int totalmem;
  double dx, dy, dz, r, minr;
  double cdepth[3], dijk[3];
  double *dtmp;
  char varname[MAXNAME];
  char *ctmp;
  dmat fitmat;
  fitpt* grdpts;
  ibook UIdx;
  dbook Uminr;
  fbook UPot, auxUPot;
  cbook Uflag;
  coord crd, auxcrd;
  prmtop *tp;
  FILE* fitrec;
  
  // If a file was requested to store the locations and intensities
  // of each fitting point, start that here.
  fitrec = NULL;
  if (myfit->ptrecfile[0] != '\0') {
    fitrec = FOpenSafe(myfit->ptrecfile, tj->OverwriteOutput);
    fprintf(fitrec, "%% Fitting points, ordered by system and conformation\n\n");
  }
  
  // Check to ensure that there is sufficient
  // room in memory for the fitting matrix   
  totalmem = myfit->ngrd;
  totalmem *= myfit->nfitpt;
  totalmem *= (myfit->q2fit+1) * 8;
  if (myfit->model == 0) {
    totalmem *= 2;
  }
  else {
    totalmem *= 5;
  }
  if (totalmem > myfit->MaxMem) {
    printf("MakeFittingMatrix >> Fitting matrix would require %lld bytes of "
           "memory.\nMakeFittingMatrix >> This exceeds the allowed %lld byte "
           "limit.  Increase\nMakeFittingMatrix >> the maximum allowed memory "
           "with the maxmem keyword in the\nMakeFittingMatrix >> &fit "
           "namelist or reduce the number of fitting points.\n", totalmem,
           myfit->MaxMem);
    exit(1);
  }
  if (myfit->verbose == 1) {
    printf("mdgx >> Fitting matrices will occupy up to %lld bytes of memory.\n",
           totalmem);
    namelen = 0;
    for (i = 0; i < myfit->ngrd; i++) {
      namelen = MAX(namelen, strlen(myfit->gname.map[i]));
    }
  }

  // Allocate space for the fitting matrix
  totalfitpt = 0;
  fitmat = CreateDmat(myfit->ngrd*myfit->nfitpt, myfit->q2fit+2, 0);
  myfit->FPtOrigins = (int*)malloc(myfit->ngrd*myfit->nfitpt*sizeof(int));

  // Allocate memory for the fitting histogram
  myfit->fitpthist = (int*)calloc(10.0/myfit->fhistbin+1, sizeof(int));
  etimers->Setup += mdgxStopTimer(etimers);

  // Loop over all grids
  for (h = 0; h < myfit->ngrd; h++) {

    // Titillate the user
    if (myfit->verbose == 1) {
      fprintf(stderr, "\rmdgx >> Composing fit for %s", myfit->gname.map[h]);
      for (i = 0; i < namelen - strlen(myfit->gname.map[h]); i++) {
        fprintf(stderr, " ");
      }
      fflush(stderr);
    }
    if (fitrec != NULL) {
      sprintf(varname, "%s", myfit->TPbank[myfit->tpidx[h]].source);
      for (i = 0; i < strlen(varname); i++) {
        if (varname[i] == '.') {
          varname[i] = '\0';
        }
      }
      fprintf(fitrec, "%s_%d = [\n", varname, h);
    }

    // Allocate coordinates for this grid's system
    tp = &myfit->TPbank[myfit->tpidx[h]];
    crd = CreateCoord(tp->natom);
    auxcrd = CreateCoord(tp->natom);
    etimers->Setup += mdgxStopTimer(etimers);

    // Order the list of points as a function of distance from  
    // the solute.  The list will then be searched in increasing
    // order of distance for candidate fitting points.
    UPot = ReadEPotGrid(myfit->gname.map[h], tp, &crd);
    etimers->Integ += mdgxStopTimer(etimers);
    Uflag = CreateCbook(UPot.pag, UPot.row, UPot.col);
    UIdx = CreateIbook(UPot.pag, UPot.row, UPot.col);
    Uminr = CreateDbook(UPot.pag, UPot.row, UPot.col, 0);
    etimers->Setup += mdgxStopTimer(etimers);
    PrepUPot(&UPot, &crd, tp, &Uflag, &Uminr, myfit);
    etimers->cellcomm += mdgxStopTimer(etimers);
    if (myfit->model == 1) {
      auxUPot = ReadEPotGrid(myfit->auxgname.map[h], tp, &auxcrd);
      etimers->Integ += mdgxStopTimer(etimers);
      congruent = 1;
      for (i = 0; i < 9; i++) {
        if (fabs(auxUPot.L.data[i] - UPot.L.data[i]) > 1.0e-8) {
          congruent = 0;
        }
      }
      for (i = 0; i < 3; i++) {
        if (fabs(auxUPot.orig[i] - UPot.orig[i]) > 1.0e-8) {
          congruent = 0;
        }
      }
      if (congruent == 0) {
        printf("MakeFittingMatrix >> Error.  Non-correspondence in grid "
               "dimensions for\nMakeFittingMatrix >> %s and %s.\n",
               myfit->auxgname.map[h], myfit->gname.map[h]);
        exit(1);
      }
      for (i = 0; i < 3*crd.natom; i++) {
        if (fabs(auxcrd.loc[i] - crd.loc[i]) > 1.0e-5) {
          congruent = 0;
        }
      }
      if (congruent == 0) {
        printf("MakeFittingMatrix >> Error.  Non-correspondence in atom "
               "positions for\nMakeFittingMatrix >> %s and %s.\n",
               myfit->auxgname.map[h], myfit->gname.map[h]);
        exit(1);
      }

      // If we're still here, this auxiliary grid is
      // clean and compatible with the original grid
    }
    etimers->Setup += mdgxStopTimer(etimers);

    // Order the fitting data
    grdpts = (fitpt*)malloc(UPot.pag*UPot.row*UPot.col*sizeof(fitpt));
    gcon = 0;
    for (i = 0; i < UPot.pag; i++) {
      for (j = 0; j < UPot.row; j++) {
        ctmp = Uflag.map[i][j];
        dtmp = Uminr.map[i][j];
        for (k = 0; k < UPot.col; k++) {

          // Skip this point if it is inaccessible
          if (ctmp[k] == 0) {
            continue;
          }

          // Catalog this grid point
          grdpts[gcon].ix = i;
          grdpts[gcon].iy = j;
          grdpts[gcon].iz = k;
          grdpts[gcon].flagged = 0;
          grdpts[gcon].minr = sqrt(dtmp[k]);
          gcon++;
        }
      }
    }
    ngpt = gcon;

    // Sort fitting points based on distance to solute
    qsort(grdpts, ngpt, sizeof(fitpt), SortPtDistance);
    etimers->nbFFT += mdgxStopTimer(etimers);

    // Label points in the index map based on the new catalog order
    for (i = 0; i < ngpt; i++) {
      UIdx.map[grdpts[i].ix][grdpts[i].iy][grdpts[i].iz] = i;
    }

    // Determine the buffer region for testing exclusions
    HessianNorms(&UPot.L, cdepth);
    buffi = myfit->prbarm/cdepth[0] + 1;
    buffj = myfit->prbarm/cdepth[1] + 1;
    buffk = myfit->prbarm/cdepth[2] + 1;

    // Loop over each grid point, selecting points for fitting 
    // when they are accessible.  Selected points will mask out
    // others nearby acccording to a minimum distance (flimit) 
    // specified in the fitting struct.                        
    setfitpt = 0;
    for (i = 0; i < ngpt; i++) {

      // Bail out if the quota for this set is reached
      if (setfitpt == myfit->nfitpt) {
        break;
      }

      // Continue if this point is already flagged
      if (grdpts[i].flagged == 1) {
        continue;
      }

      // Accept the point with some probability based
      // on its distance from the molecular surface  
      if (SelectFitPoint(&grdpts[i], myfit, tj) == 1) {
        if (myfit->model == 0) {
          ContributeFitPt(myfit, h, &fitmat, totalfitpt, &grdpts[i], &crd,
                          &UPot, &UPot, Atm2Col, fitrec);
        }
        else {
          ContributeFitPt(myfit, h, &fitmat, totalfitpt, &grdpts[i], &crd,
                          &UPot, &auxUPot, Atm2Col, fitrec);
        }
        FlagFitPt(myfit, grdpts, i, &UPot, &Uflag, &UIdx, buffi, buffj, buffk);
        HistogramFitPt(myfit, &grdpts[i]);
        totalfitpt++;
        setfitpt++;
      }
      else {

        // Even if the point is rejected, the region around
        // it must be flagged in order to achieve the point
        // spread that is desired.                         
        FlagFitPt(myfit, grdpts, i, &UPot, &Uflag, &UIdx, buffi, buffj, buffk);
      }
    }
    etimers->nbBsp += mdgxStopTimer(etimers);
    if (fitrec != NULL) {
      fprintf(fitrec, "];\n\n");
    }

    // Commit coordinates to buffer for printing if necessary
    for (i = 0; i < 3*tp->natom; i++) {
      allcrd->map[h][i] = crd.loc[i];
    }

    // Print coordinates to conformation file
    WriteConformation(h, crd.loc, tp, myfit, tj);

    // Free allocated memory
    DestroyFbook(&UPot);
    if (myfit->model == 1) {
      DestroyFbook(&auxUPot);
    }
    DestroyIbook(&UIdx);
    DestroyDbook(&Uminr);
    DestroyCbook(&Uflag);
    free(grdpts);
    DestroyCoord(&crd);
    DestroyCoord(&auxcrd);
  }

  // Update the size of the fitting matrix
  if (totalfitpt < fitmat.row) {
    fitmat = ReallocDmat(&fitmat, totalfitpt, fitmat.col);
  }
  
  // Titillate the user
  if (myfit->verbose == 1) {
    printf("\nmdgx >> Fitting matrix composed.\n");
    totalmem = fitmat.row;
    totalmem *= fitmat.col * 8;
    if (myfit->model == 0) {
      totalmem *= 2;
    }
    else {
      totalmem *= 5;
    }
    printf("mdgx >> Refined guess for matrix memory requirement: %lld bytes\n"
           "mdgx >> This does not include constraint equations and other "
           "incidental\nmdgx >> costs, but can serve as a better estimate "
           "for cases that may require\nmdgx >> large memory allocations.\n",
           totalmem);
  }
  if (fitrec != NULL) {
    fclose(fitrec);
  }
  
  return fitmat;
}

//-----------------------------------------------------------------------------
// CalcGroupQ: calculate the total charge of a group of atoms in a system; it
//             could be the entire system.
//
// Arguments:                                                           
//   q:        the array of charge variables                            
//   mask:     the mask of atoms defining this group                    
//   Atm2Col:  the key by which atoms map to the charges in q           
//   natom:    the number of atoms in the system (length of mask)       
//-----------------------------------------------------------------------------
static double CalcGroupQ(double* q, int* mask, int* Atm2Col, int natom)
{
  int i;
  double totq;

  totq = 0.0;
  for (i = 0; i < natom; i++) {
    totq += q[Atm2Col[i]] * (double)(mask[i]);
  }

  return totq;
}

//-----------------------------------------------------------------------------
// EvalGroupQSums: evaluate all charge group sums.
//
// Arguments:
//   q:          the array of charge variables                            
//   myfit:      fitting control data
//   Atm2Col:    matrix of atom indices into the fitting matrix columns   
//   qsummasks:  pre-evaluated masks for all charge group sum constraints
//-----------------------------------------------------------------------------
static double EvalGroupQSums(double* q, fset *myfit, imat *Atm2Col,
                             imat *qsummasks, int* qsumexample)
{
  int i, j, ngrp;
  int* atmmask;
  double dq, rmsdq;
  prmtop *tp;
  coord crd;

  ngrp = 0;
  rmsdq = 0.0;

  // Evaluate dq for systems
  for (i = 0; i < myfit->tpcount; i++) {
    tp = &myfit->TPbank[i];

    dq = myfit->totalq[i];
    for (j = 0; j < tp->natom; j++) {
      dq -= q[Atm2Col->map[i][j]];
    }
    rmsdq += dq*dq;
    ngrp++;
  }
  
  // Evaluate dq for any relevant charge group sum restraints
  for (i = 0; i < myfit->nqsum; i++) {
    j = qsumexample[i];

    // Bail out if there was no instance of this charge sum restraint
    if (j < 0) {
      continue;
    }
    tp = &myfit->TPbank[j];
    if (ISum(qsummasks[j].map[i], tp->natom) > 0) {
      dq = CalcGroupQ(q, qsummasks[j].map[i], Atm2Col->map[j], tp->natom) -
           myfit->qsum[i].target;
      rmsdq += dq*dq;
      ngrp++;
    }
  }

  return sqrt(rmsdq / ngrp);
}

//-----------------------------------------------------------------------------
// SnapFittedQ: snap the fitted set of charges to the nearest 1.0e-5 proton
//              units, ensure that all charges that were intended to be
//              equalized are perfectly equalized, ensure that the sum of all
//              charges adds to the requested value, and ensure that all
//              charges that were fixed are fixed.      
//                                                                      
// Arguments:                                                           
//   q:         the fitted charges                                      
//   myfit:     fitting control data                                    
//   Atm2Col:   matrix of atom indices into the fitting matrix columns   
//   qfix:      charge fix restraints                                    
//   npert:     the number of charge perturbations to attempt at one time
//-----------------------------------------------------------------------------
static int SnapFittedQ(double* q, fset *myfit, imat *Atm2Col, nail* qfix,
                       int npert)
{
  int i, j, k, niter, npass;
  int* atmmask;
  int* fixqmask;
  int* qsumexample;
  double chksum, Sbest, Scurr, qiorig, qjorig;
  double* qbest;
  imat* qsummasks;
  prmtop *tp;
  coord crd;

  // Titillate the user
  if (myfit->verbose == 1) {
    printf("mdgx >> Snapping fitted charges, rank %d.\n", npert);
  }
  
  // Round all charges to the nearest hundred thousandth of a unit
  for (i = 0; i < myfit->q2fit; i++) {
    if (q[i] >= 0.0) {
      j = q[i]*1.0e5;
      q[i] = j*1.0e-5;
    }
    else {
      j = -q[i]*1.0e5;
      q[i] = -j*1.0e-5;
    }
  }

  // Make the atom masks for all systems and all qsum restraints just once
  qsummasks = (imat*)malloc(myfit->tpcount * sizeof(imat));
  qsumexample = (int*)malloc(myfit->nqsum * sizeof(int));
  SetIVec(qsumexample, myfit->nqsum, -1);
  for (i = 0; i < myfit->tpcount; i++) {
    tp = &myfit->TPbank[i];
    qsummasks[i] = CreateImat(myfit->nqsum, tp->natom);
    crd = CreateCoord(tp->natom);
    for (j = 0; j < myfit->nqsum; j++) {
      atmmask = ParseAmbMask(myfit->qsum[j].maskstr, tp, &crd);
      ReflectIVec(qsummasks[i].map[j], atmmask, tp->natom);
      if (ISum(atmmask, tp->natom) > 0) {
        qsumexample[j] = i;
      }
      free(atmmask);
    }
    DestroyCoord(&crd);
  }  
  
  // Freeze fixed charges
  fixqmask = (int*)calloc(myfit->q2fit, sizeof(int));
  for (i = 0; i < myfit->nqfix; i++) {
    for (j = 0; j < myfit->tpcount; j++) {
      tp = &myfit->TPbank[j];
      atmmask = ParseAmbMask(qfix[i].maskstr, tp, &crd);
      for (k = 0; k < tp->natom; k++) {
        if (atmmask[k] == 1) {
          q[Atm2Col->map[j][k]] = qfix[i].target;
          fixqmask[Atm2Col->map[j][k]] = 1;
        }
      }
      free(atmmask);
    }
  }

  // Count the number of charge sums that already pass
  npass = 0;
  atmmask = (int*)calloc(myfit->nqsum, sizeof(int));
  for (i = 0; i < myfit->tpcount; i++) {
    tp = &myfit->TPbank[i];
    for (j = 0; j < myfit->nqsum; j++) {
      if (ISum(qsummasks[i].map[j], tp->natom) > 0) {
        chksum = CalcGroupQ(q, qsummasks[i].map[j], Atm2Col->map[i],
                            tp->natom) - myfit->qsum[j].target;
        if (fabs(chksum) < 1.0e-10) {
          atmmask[j] = 1;
        }
      }
    }
  }
  if (myfit->verbose == 1) {
    printf("mdgx >>   %4d / %4d sums pass.  ", ISum(atmmask, myfit->nqsum),
           myfit->nqsum);
  }
  free(atmmask);
  
  // Work on the charge variables, one at a time.  Use steepest
  // descent optimization in set increments of 1.0e-5 proton   
  // charges to minimize (hopefully eliminate) any deviations  
  // from the ideal charge sum values.                         
  qbest = CpyDVec(q, myfit->q2fit);
  Sbest = EvalGroupQSums(q, myfit, Atm2Col, qsummasks, qsumexample);
  if (myfit->verbose == 1) {
    printf("Average error: %9.6lfe\n", Sbest);
  }
  niter = 0;
  if (npert == 1) {
    while (niter < myfit->maxsnap && Sbest > 1.0e-10) {
      ReflectDVec(qbest, q, myfit->q2fit);
      for (i = 0; i < myfit->q2fit; i++) {

        // Skip fixed charges
        if (fixqmask[i] == 1) {
          continue;
        }

        // Initial loop perturbing one charge at a time
        qiorig = q[i];
        q[i] = qiorig - 1.0e-5;
        Scurr = EvalGroupQSums(q, myfit, Atm2Col, qsummasks, qsumexample);
        if (Scurr < Sbest - 1.0e-10) {
          Sbest = Scurr;
          ReflectDVec(qbest, q, myfit->q2fit);
          continue;
        }
        q[i] = qiorig + 1.0e-5;
        Scurr = EvalGroupQSums(q, myfit, Atm2Col, qsummasks, qsumexample);
        if (Scurr < Sbest - 1.0e-10) {
          Sbest = Scurr;
          ReflectDVec(qbest, q, myfit->q2fit);
          continue;
        }
        q[i] = qiorig;
      }
      ReflectDVec(q, qbest, myfit->q2fit);
      niter++;
    }
  }
  else {
    while (niter < myfit->maxsnap && Sbest > 1.0e-10) {
      ReflectDVec(qbest, q, myfit->q2fit);
      for (i = 0; i < myfit->q2fit; i++) {

        // Skip fixed charges
        if (fixqmask[i] == 1) {
          continue;
        }
        for (j = i+1; j < myfit->q2fit; j++) {
      
          // Skip fixed charges
          if (fixqmask[j] == 1) {
            continue;
          }

          // A pair of perturbable charges; test all combinations
          qiorig = q[i];
          qjorig = q[j];
          q[i] = qiorig - 1.0e-5;
          q[j] = qjorig - 1.0e-5;
          Scurr = EvalGroupQSums(q, myfit, Atm2Col, qsummasks, qsumexample);
          if (Scurr < Sbest - 1.0e-10) {
            Sbest = Scurr;
            ReflectDVec(qbest, q, myfit->q2fit);
            continue;
          }
          q[i] = qiorig - 1.0e-5;
          q[j] = qjorig + 1.0e-5;
          Scurr = EvalGroupQSums(q, myfit, Atm2Col, qsummasks, qsumexample);
          if (Scurr < Sbest - 1.0e-10) {
            Sbest = Scurr;
            ReflectDVec(qbest, q, myfit->q2fit);
            continue;
          }
          q[i] = qiorig + 1.0e-5;
          q[j] = qjorig - 1.0e-5;
          Scurr = EvalGroupQSums(q, myfit, Atm2Col, qsummasks, qsumexample);
          if (Scurr < Sbest - 1.0e-10) {
            Sbest = Scurr;
            ReflectDVec(qbest, q, myfit->q2fit);
            continue;
          }
          q[i] = qiorig + 1.0e-5;
          q[j] = qjorig + 1.0e-5;
          Scurr = EvalGroupQSums(q, myfit, Atm2Col, qsummasks, qsumexample);
          if (Scurr < Sbest - 1.0e-10) {
            Sbest = Scurr;
            ReflectDVec(qbest, q, myfit->q2fit);
            continue;
          }
          q[i] = qiorig;
          q[j] = qjorig;
        }
      }
      ReflectDVec(q, qbest, myfit->q2fit);
      niter++;
    }
  }

  // Store the best charge set
  ReflectDVec(q, qbest, myfit->q2fit);

  // Free allocated memory
  free(qbest);
  free(fixqmask);
  for (i = 0; i < myfit->tpcount; i++) {
    DestroyImat(&qsummasks[i]);
  }
  free(qsummasks);
  free(qsumexample);
  
  if (niter == myfit->maxsnap && Sbest > 1.0e-10) {
    return -1;
  }
  else {
    return niter;
  }
}

//-----------------------------------------------------------------------------
// SystemQTest: test the accuracy of the electrostatics for each system.
//                                                                      
// Arguments:                                                           
//   btar:       the vector of target results                           
//   bpred:      predicted electrostatic potential values               
//   npt:        the number of elements in btar, bpred                  
//   myfit:      the fitting command data                               
//   row:        row of the results matrices to write                   
//-----------------------------------------------------------------------------
static void SystemQTest(double* btar, double* bpred, int npt, fset *myfit,
                        int row)
{
  int i, j, nspt;
  double* ttar;
  double* tpred;

  // Allocate scratch arrays
  ttar = (double*)malloc(npt*sizeof(double));
  tpred = (double*)malloc(npt*sizeof(double));

  // Loop over all systems
  for (i = 0; i < myfit->tpcount; i++) {
    nspt = 0;
    for (j = 0; j < npt; j++) {
      if (myfit->FPtOrigins[j] == i) {
        ttar[nspt] = btar[j];
        tpred[nspt] = bpred[j];
        nspt++;
      }
    }
    myfit->QScorr.map[row][i] = Pearson(ttar, tpred, nspt);
    myfit->QSrmsd.map[row][i] = VecRMSD(ttar, tpred, nspt);
  }

  // Free allocated memory
  free(ttar);
  free(tpred);
}

//-----------------------------------------------------------------------------
// SolvRESP: solve a linear least squares problem and report statistics on the
//           results.  The fitting matrix is stored in the first n-1 columns of
//           the fitting and constraint matrices; the target vector is stored
//           in the final column, for simplified data passing.
//
// Arguments:                                                           
//   fitmat:     the fitting matrix and target vector derived from      
//               electrostatic potential data                           
//   cnstmat:    the constraints matrix and target vector               
//   myfit:      the fitting command data                               
//   etimers:    the execution timings data                             
//-----------------------------------------------------------------------------
static dmat SolveRESP(dmat *fitmat, dmat *cnstmat, fset *myfit, imat *Atm2Col,
                      execon *etimers)
{
  int i, j, ip, tval, nqrest, nqrest2;
  int* replaceGeneralCnst;
  double* bfit;
  double* btst;
  double* bpred;
  dmat Afit, Atst, qres;

  // Titillate the user
  if (myfit->verbose == 1) {
    printf("mdgx >> Solving linear least-squares problem for %d independent "
           "charges.\n", myfit->q2fit);
  }

  // Fit new charges
  if (myfit->model == 0) {
    i = 1;
    j = 1;
  }
  else {
    i = 4;
    j = 3;
  }
  qres = CreateDmat(i, myfit->q2fit, 0);
  myfit->QScorr = CreateDmat(j, myfit->tpcount, 0);
  myfit->QSrmsd = CreateDmat(j, myfit->tpcount, 0);

  // Count the number of restraints that apply to this fit:  
  // total q, minimize q, sum q, tether q, and fix solvated q
  nqrest = 0;
  for (i = 0; i < cnstmat->row; i++) {
    tval = cnstmat->map[i][myfit->q2fit+1];
    if (tval == 0 || tval == 1 || tval == 2 || tval == 3 || tval == 4) {
      nqrest++;
    }
  }
  Afit = CreateDmat(fitmat->row+nqrest, myfit->q2fit, 0);
  bfit = (double*)malloc((fitmat->row+nqrest)*sizeof(double));
  for (i = 0; i < fitmat->row; i++) {
    for (j = 0; j < myfit->q2fit; j++) {
      Afit.map[i][j] = fitmat->map[i][j];
    }
    bfit[i] = fitmat->map[i][myfit->q2fit];
  }
  ip = fitmat->row;
  for (i = 0; i < cnstmat->row; i++) {
    tval = cnstmat->map[i][myfit->q2fit+1];
    if (tval == 0 || tval == 1 || tval == 2 || tval == 3 || tval == 4) {
      for (j = 0; j < myfit->q2fit; j++) {
        Afit.map[ip][j] = cnstmat->map[i][j];
      }
      bfit[ip] = cnstmat->map[i][myfit->q2fit];
      ip++;
    }
  }
  AxbQRRxc(Afit, bfit, myfit->verbose);
  BackSub(Afit, bfit);
  DestroyDmat(&Afit);
  etimers->nbPtM += mdgxStopTimer(etimers);

  // Test the result
  Atst = CreateDmat(fitmat->row, myfit->q2fit, 0);
  btst = (double*)malloc(fitmat->row*sizeof(double));
  bpred = (double*)malloc(fitmat->row*sizeof(double));
  for (i = 0; i < fitmat->row; i++) {
    for (j = 0; j < myfit->q2fit; j++) {
      Atst.map[i][j] = fitmat->map[i][j];
    }
    btst[i] = fitmat->map[i][myfit->q2fit];
  }
  DMatVecMult(&Atst, bfit, bpred);
  SystemQTest(btst, bpred, fitmat->row, myfit, 0);

  // Store this set of charges along with its statistics
  ReflectDVec(qres.map[0], bfit, myfit->q2fit);
  myfit->SnapCount[0] = SnapFittedQ(qres.map[0], myfit, Atm2Col,
                                    myfit->qfixsolv, 1);
  if (myfit->SnapCount[0] < 0) {
    myfit->SnapCount[0] = SnapFittedQ(qres.map[0], myfit, Atm2Col,
                                      myfit->qfixsolv, 2);
  }
  myfit->Qcorr[0] = Pearson(btst, bpred, Atst.row);
  myfit->Qrmsd[0] = VecRMSD(btst, bpred, Atst.row);

  // Free allocated memory
  DestroyDmat(&Atst);
  free(bfit);
  free(btst);
  free(bpred);
  etimers->nbCnv += mdgxStopTimer(etimers);

  // If this is IPolQ, try a different approach
  if (myfit->model == 1) {

    // Count the number of restraints that apply to solvated and vacuum
    // fits:                                                           
    //
    // vacuum q: total q, minimize q, sumq, tether q, and fix vacuum q
    // delta q: total q, sumq, and fix solvated q
    //
    // Also add additional, modest minimization restraints on each delta q
    // and amplify them if the corresponding vacuum q is held to a special
    // minimum value by some particular amount or if the vacuum q is
    // tethered.
    nqrest = 0;
    nqrest2 = 0;
    replaceGeneralCnst = (int*)calloc(myfit->q2fit, sizeof(int));
    for (i = 0; i < cnstmat->row; i++) {
      tval = cnstmat->map[i][myfit->q2fit+1];
      if (tval == 0 || tval == 1 || tval == 2 || tval == 3 || tval == 5) {
        nqrest++;
      }
      if (tval == 0 || tval == 2 || tval == 4) {
        nqrest2++;
      }
      if (tval == 1 || tval == 3) {
        for (j = 0; j < myfit->q2fit; j++) {
          if (fabs(cnstmat->map[i][j]) > 1.0e-4) {
            replaceGeneralCnst[j] = 1;
          }
        }
      }
    }
    nqrest2 += myfit->q2fit;
    Afit = CreateDmat(2*fitmat->row + nqrest + nqrest2, 2*myfit->q2fit, 0);
    bfit = (double*)malloc((2*fitmat->row + nqrest + nqrest2) *
                           sizeof(double));
    for (i = 0; i < fitmat->row; i++) {
      for (j = 0; j < myfit->q2fit; j++) {
        Afit.map[i][j] = fitmat->map[i][j];
      }
      bfit[i] = fitmat->map[i][myfit->q2fit+1];
    }
    ip = fitmat->row;
    for (i = 0; i < fitmat->row; i++) {
      for (j = 0; j < myfit->q2fit; j++) {
        Afit.map[ip][j] = fitmat->map[i][j];
        Afit.map[ip][j + myfit->q2fit] = fitmat->map[i][j];
      }
      bfit[ip] = fitmat->map[i][myfit->q2fit];
      ip++;
    }
    for (i = 0; i < cnstmat->row; i++) {
      tval = cnstmat->map[i][myfit->q2fit+1];
      if (tval == 0 || tval == 1 || tval == 2 || tval == 3 || tval == 5) {
        for (j = 0; j < myfit->q2fit; j++) {
          Afit.map[ip][j] = cnstmat->map[i][j];
        }
        bfit[ip] = cnstmat->map[i][myfit->q2fit];
        ip++;
      }
    }
    for (i = 0; i < cnstmat->row; i++) {
      tval = cnstmat->map[i][myfit->q2fit+1];
      if (tval == 0 || tval == 2 || tval == 4) {
        for (j = 0; j < myfit->q2fit; j++) {
          Afit.map[ip][j] = cnstmat->map[i][j];
          Afit.map[ip][j + myfit->q2fit] = cnstmat->map[i][j];
        }
        bfit[ip] = cnstmat->map[i][myfit->q2fit];
        ip++;
      }
    }

    // Restraints to keep charge deltas small, or squash them entirely
    for (i = 0; i < myfit->q2fit; i++) {
      if (replaceGeneralCnst[i] == 1) {
        Afit.map[ip][i+myfit->q2fit] = 1.0e7 * myfit->nsamp[i];
      }
      else {
        Afit.map[ip][i+myfit->q2fit] = myfit->qtthwt * myfit->nsamp[i];
      }
      bfit[ip] = 0.0;
      ip++;
    }
    AxbQRRxc(Afit, bfit, myfit->verbose);
    BackSub(Afit, bfit);
    DestroyDmat(&Afit);
    etimers->nbMtP += mdgxStopTimer(etimers);

    // Free allocated memory
    free(replaceGeneralCnst);

    // Store the vacuum, polarized, and delta charge sets
    ReflectDVec(qres.map[1], bfit, myfit->q2fit);
    myfit->SnapCount[1] = SnapFittedQ(qres.map[1], myfit, Atm2Col,
                                      myfit->qfixvacu, 1);
    if (myfit->SnapCount[1] < 0) {
      myfit->SnapCount[1] = SnapFittedQ(qres.map[1], myfit, Atm2Col,
                                        myfit->qfixvacu, 2);
    }
    ReflectDVec(qres.map[2], bfit, myfit->q2fit);
    DVec2VecAdd(qres.map[2], &bfit[myfit->q2fit], myfit->q2fit);
    myfit->SnapCount[2] = SnapFittedQ(qres.map[2], myfit, Atm2Col,
                                      myfit->qfixsolv, 1);
    if (myfit->SnapCount[2] < 0) {
      myfit->SnapCount[2] = SnapFittedQ(qres.map[2], myfit, Atm2Col,
                                        myfit->qfixsolv, 2);
    }
    ReflectDVec(qres.map[3], &bfit[myfit->q2fit], myfit->q2fit);

    // Test the result
    Atst = CreateDmat(fitmat->row, myfit->q2fit, 0);
    btst = (double*)malloc(fitmat->row*sizeof(double));
    bpred = (double*)malloc(fitmat->row*sizeof(double));
    for (i = 0; i < fitmat->row; i++) {
      for (j = 0; j < myfit->q2fit; j++) {
        Atst.map[i][j] = fitmat->map[i][j];
      }
      btst[i] = fitmat->map[i][myfit->q2fit+1];
    }
    DMatVecMult(&Atst, bfit, bpred);
    myfit->Qcorr[1] = Pearson(btst, bpred, Atst.row);
    myfit->Qrmsd[1] = VecRMSD(btst, bpred, Atst.row);
    SystemQTest(btst, bpred, fitmat->row, myfit, 1);
    DestroyDmat(&Atst);
    Atst = CreateDmat(fitmat->row, 2*myfit->q2fit, 0);
    for (i = 0; i < fitmat->row; i++) {
      for (j = 0; j < myfit->q2fit; j++) {
        Atst.map[i][j] = fitmat->map[i][j];
        Atst.map[i][j+myfit->q2fit] = fitmat->map[i][j];
      }
      btst[i] = fitmat->map[i][myfit->q2fit];
    }
    DMatVecMult(&Atst, bfit, bpred);
    myfit->Qcorr[2] = Pearson(btst, bpred, Atst.row);
    myfit->Qrmsd[2] = VecRMSD(btst, bpred, Atst.row);
    SystemQTest(btst, bpred, fitmat->row, myfit, 2);
    DestroyDmat(&Atst);
    etimers->nbCnv += mdgxStopTimer(etimers);
  }

  // Titillate the user
  if (myfit->verbose == 1) {
    printf("mdgx >> Fit complete.\n");
  }

  return qres;
}

//-----------------------------------------------------------------------------
// CalculateDipoles: this function takes the best fitted charge set and puts it
//                   on each of the fitting conformations to compute dipole
//                   moments.  Generates additional text in the output file.
//
// Arguments:                                                           
//   myfit:    the fitting data and input variables                     
//   qres:     the fitted charges                                       
//   allcrd:   matrix of all coordinates for all systems                
//   Atm2Col:  matrix of atom indices into the fitting matrix columns   
//   outp:     the output file being written                            
//-----------------------------------------------------------------------------
static void CalculateDipoles(fset *myfit, dmat *qres, dmat *allcrd,
                             imat *Atm2Col, FILE *outp)
{
  int h, i, j, k, i3, gnamelen, tpnamelen, nrep;
  int *qidx;
  double DPdev;
  double* sqscr;
  dmat dp;
  dmat sqdp;
  prmtop *tp;

  // Allocate memory and compute dipoles
  dp = CreateDmat(3, 3*myfit->ngrd, 0);
  sqdp = CreateDmat(3, 3*myfit->ngrd, 0);
  sqscr = (double*)malloc(3*myfit->ngrd*sizeof(double));
  for (i = 0; i < myfit->ngrd; i++) {

    // The dipole is not relevant if the molecule is charged
    if (fabs(myfit->totalq[myfit->tpidx[i]]) > 1.0e-8) {
      continue;
    }
    tp = &myfit->TPbank[myfit->tpidx[i]];
    qidx = Atm2Col->map[myfit->tpidx[i]];
    for (j = 0; j < tp->natom; j++) {
      for (k = 0; k < 3; k++) {
        dp.map[0][3*i+k] += qres->map[0][qidx[j]]*allcrd->map[i][3*j+k];
        if (myfit->model == 1) {
          dp.map[1][3*i+k] += qres->map[1][qidx[j]]*allcrd->map[i][3*j+k];
          dp.map[2][3*i+k] += qres->map[2][qidx[j]]*allcrd->map[i][3*j+k];
        }
      }
    }
  }

  // Print outputs
  HorizontalRule(outp, 0);
  fprintf(outp, "(4.) Dipole moments on all conformations, (units of Debye)"
          "\n\n");
  tpnamelen = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    j = strlen(myfit->TPbank[i].source);
    if (tpnamelen < j) {
      tpnamelen = j;
    }
  }
  for (i = 0; i < myfit->ngrd; i++) {
    i3 = 3*i;
    for (j = 0; j < 3; j++) {
      sqdp.map[j][i] =
      sqrt(dp.map[j][i3]*dp.map[j][i3] + dp.map[j][i3+1]*dp.map[j][i3+1] +
           dp.map[j][i3+2]*dp.map[j][i3+2]) * EA2DEBYE;
    }
  }
  gnamelen = 0;
  for (i = 0; i < myfit->ngrd; i++) {
    j = strlen(myfit->gname.map[i]);
    if (gnamelen < j) {
      gnamelen = j;
    }
  }
  if (gnamelen < 9) {
    gnamelen = 9;
  }
  if (myfit->DispAllDP == 1) {
    if (myfit->model == 0) {
      fprintf(outp, "  %-*s  %-*s IPolQ,orig    Vacuum  IPolQ,pert\n",
              gnamelen, " ", tpnamelen, " ");
    }
    fprintf(outp, "  %-*s  %-*s", gnamelen, "Grid File", tpnamelen, "System");
    fprintf(outp, "  ");
    for (i = 0; i < gnamelen; i++) {
      fprintf(outp, "-");
    }
    fprintf(outp, "  ");
    for (i = 0; i < tpnamelen; i++) {
      fprintf(outp, "-");
    }
    fprintf(outp, "\n");
    for (i = 0; i < myfit->ngrd; i++) {
      fprintf(outp, "  %-*s  %-*s", gnamelen, myfit->gname.map[i],
              tpnamelen, myfit->TPbank[myfit->tpidx[i]].source);
      for (j = 0; j < 3; j++) {
        fprintf(outp, "  %9.5lf", sqdp.map[j][i]);
      }
      fprintf(outp, "\n");
    }
    fprintf(outp, "\n");
  }
  if (myfit->model == 0) {
    fprintf(outp, "  %-*s  %-*s       Result\n", gnamelen, " ",
            tpnamelen, " ");
  }
  else if (myfit->model == 1) {
    fprintf(outp, "  %-*s     IPolQ,orig        Vacuum        IPolQ,pert\n",
            tpnamelen, " ");
  }
  fprintf(outp, " System");
  for (i = 0; i < tpnamelen-6; i++) {
    fprintf(outp, " ");
  }
  nrep = (myfit->model == 0) ? 1 : 3;
  for (i = 0; i < nrep; i++) {
    fprintf(outp, "    Mean   Stdev.");
  }
  fprintf(outp, "\n ");
  for (i = 0; i < tpnamelen; i++) {
    fprintf(outp, "-");
  }
  for (i = 0; i < nrep; i++) {
    fprintf(outp, "   ------  ------");
  }
  fprintf(outp, "\n");
  for (h = 0; h < myfit->tpcount; h++) {
    fprintf(outp, " %-*s", tpnamelen, myfit->TPbank[h].source);
    for (i = 0; i < nrep; i++) {
      k = 0;
      for (j = 0; j < myfit->ngrd; j++) {
        if (myfit->tpidx[j] == h) {
          sqscr[k] = sqdp.map[i][j];
          k++;
        }
      }
      DPdev = (k >= 2) ? DStDev(sqscr, k) : 0.0;
      fprintf(outp, "   %6.2lf  %6.2lf", DAverage(sqscr, k), DPdev);
    }
    fprintf(outp, "\n");
  }
  HorizontalRule(outp, 1);

  // Free allocated memory
  DestroyDmat(&dp);
  DestroyDmat(&sqdp);
  free(sqscr);
}

//-----------------------------------------------------------------------------
// QFitTimingsData: print the timings data from this charge fitting run.
//
// Arguments:                                                           
//   etimers:     the execution timings struct                          
//   outp:        the mdout file currently being written                
//-----------------------------------------------------------------------------
static void QFitTimingsData(execon *etimers, FILE *outp)
{
  double tt;

  tt = etimers->Setup +     // Setup time (memory allocation, constraint calc)
       etimers->cellcomm +  // Grid coloring time for fitting data
       etimers->nbPtM +     // Small REsP problem time
       etimers->Integ +     // Grid reading time
       etimers->nbFFT +     // Point sorting time
       etimers->nbBsp;      // Point marking time
  if (etimers->nbMtP > 1.0e-8) {
    tt += etimers->nbMtP;   // Large IPolQ problem time
  }
  tt = 100.0/tt;

  HorizontalRule(outp, 0);
  fprintf(outp, "(6.) Timings for the fitting problem.\n\n");
  fprintf(outp, " Segment                 Time / Percentage\n"
          " ----------------------  ---------  ------\n");
  fprintf(outp, " Setup                   %9.2lf  %6.2lf\n", etimers->Setup,
          etimers->Setup*tt);
  fprintf(outp, " Grid Reading            %9.2lf  %6.2lf\n", etimers->Integ,
          etimers->Integ*tt);
  fprintf(outp, " Grid Processing         %9.2lf  %6.2lf\n", etimers->cellcomm,
          etimers->cellcomm*tt);
  fprintf(outp, " Point Sorting           %9.2lf  %6.2lf\n", etimers->nbFFT,
          etimers->nbFFT*tt);
  fprintf(outp, " Point Selection         %9.2lf  %6.2lf\n", etimers->nbBsp,
          etimers->nbBsp*tt);
  fprintf(outp, " Least Squares, REsP     %9.2lf  %6.2lf\n", etimers->nbPtM,
          etimers->nbPtM*tt);
  if (etimers->nbMtP > 1.0e-8) {
    fprintf(outp, " Least Squares, IPolQ    %9.2lf  %6.2lf\n", etimers->nbMtP,
            etimers->nbMtP*tt);
  }
  fprintf(outp, " Output Printing         %9.2lf  %6.2lf\n", etimers->Write,
          etimers->Write*tt);
  fprintf(outp, " Total CPU Time          %9.2lf  %6.2lf\n", 100.0/tt, 100.0);
  HorizontalRule(outp, 1);
}

//-----------------------------------------------------------------------------
// PrintEPRules: print a file of extra point rules that will modify each
//               topology and add any needed extra points.              
//                                                                      
// Arguments:                                                           
//   myfit:      fitting control data                                   
//   tpnum:      the number of the topology for which to print rules (the
//               fitted charges for each topology are now stored in the
//               topology structs themselves)                       
//   tj:         trajectory control information                         
//-----------------------------------------------------------------------------
static void PrintEPRules(fset *myfit, int tpnum, trajcon *tj)
{
  int i, j, natm, nres;
  eprule *tmr;
  char outname[MAXNAME];
  FILE *outp;
  prmtop *tp;

  tp = &myfit->TPbank[tpnum];
  sprintf(outname, "%s.%s", tp->source, myfit->epext);
  outp = FOpenSafe(outname, tj->OverwriteOutput);
  fprintf(outp, "%% Extra points rules for charge model fitted as described "
          "in\n%% %s.\n%%\n\n", tj->outbase);
  fprintf(outp, "%% The following rules will modify the charges of\n%% atoms "
          "given in the topology %s.\n", tp->source);
  for (i = 0; i < tp->natom; i++) {
    nres = LocateResID(tp, i, 0, tp->nres);
    if (tp->EPInserted == 0 || tp->OldAtomNum[i] >= 0) {
      fprintf(outp, "&rule\n  ResidueName  %.4s\n  ExtraPoint   %.4s\n  "
              "FrameStyle   0\n  Charge       %9.5lf\n&end\n\n",
              &tp->ResNames[4*nres], &tp->AtomNames[4*i],
              tp->Charges[i]);
    }
  }
  if (tp->EPInserted == 0) {
    return;
  }
  fprintf(outp, "%% The following rules will add charged extra points\n%% "
          "to the original topology.\n");
  for (i = 0; i < tp->neprule; i++) {
    tmr = &tp->eprules[i];
    for (j = 0; j < tp->natom; j++) {
      if (strncmp(&tp->AtomNames[4*j], tmr->epname, 4) == 0) {
        natm = j;
      }
    }
    fprintf(outp, "&rule\n  ResidueName  %.4s\n  ExtraPoint   %.4s\n  "
             "FrameStyle   %d\n", tmr->resname, tmr->epname, tmr->frstyle);
    fprintf(outp, "  FrameAtom1   %.4s\n  FrameAtom2   %.4s\n", tmr->fr1,
            tmr->fr2);
    if (tmr->frstyle > 1) {
      fprintf(outp, "  FrameAtom3   %.4s\n", tmr->fr3);
    }
    if (tmr->frstyle == 6 || tmr->frstyle == 9 || tmr->frstyle == 10) {
      fprintf(outp, "  FrameAtom4   %.4s\n", tmr->fr4);
    }
    fprintf(outp, "  Vector12     %16.10lf\n", tmr->d1);
    if (tmr->frstyle == 2 || tmr->frstyle == 5 || tmr->frstyle == 6) {
      fprintf(outp, "  Vector13     %16.10lf\n", tmr->d2);
    }
    else if (tmr->frstyle == 4) {
      fprintf(outp, "  Theta        %16.10lf\n", tmr->d2);
    }
    if (tmr->frstyle == 3) {
      fprintf(outp, "  Vector23     %16.10lf\n", tmr->d3);
    }
    else if (tmr->frstyle == 5) {
      fprintf(outp, "  Vector12x13  %16.10lf\n", tmr->d3);
    }
    fprintf(outp, "  Charge       %9.5lf\n", tp->Charges[natm]);
    if (tmr->sig >= 0.0) {
      fprintf(outp, "  Sigma        %16.10lf\n", tmr->sig);
    }
    if (tmr->eps >= 0.0) {
      fprintf(outp, "  Epsilon      %16.10lf\n", tmr->eps);
    }
    fprintf(outp, "&end\n\n");
  }
  fclose(outp);
}

//-----------------------------------------------------------------------------
// PrintRespHistogram: print a histogram of the RESP fitting points, indicating
//                     their minimum distance to the solute. 
//
// Arguments:                                                           
//   myfit:      the fitting command data, used here to determine the   
//               number of jackknife fits and their style               
//-----------------------------------------------------------------------------
static void PrintRespHistogram(fset *myfit, trajcon *tj)
{
  int i;
  FILE *outp;

  outp = FOpenSafe(myfit->histfile, tj->OverwriteOutput);
  fprintf(outp, "%% Molecule:fit point distance histogram for charge model "
          "fitted\n%% as described in %s.\n%%\n%% Counts are normalized by "
          "the total number of fitting points.\n\n%% Bin Center     Count\n"
          "%% ----------   ---------\n", tj->outbase);
  for (i = 0; i < 10.0/myfit->fhistbin; i++) {
    fprintf(outp, "  %10.6lf   %9.6lf\n", (i+0.5)*myfit->fhistbin,
            (double)myfit->fitpthist[i]/(myfit->ngrd*myfit->nfitpt));
  }
  fclose(outp);
}

//-----------------------------------------------------------------------------
// UnitResidueMatch: ensure that the residue name matches the unit name.  Some
//                   degree of leniency is granted if the user specifies.
//
// Arguments:
//  
//-----------------------------------------------------------------------------
static int UnitResidueMatch(fset *myfit, char* resname, char* unitname)
{
  int i, lres, lunit;
  char unitbuff[8], resbuff[8];

  Type4Char(unitbuff, unitname);
  for (i = 0; i < 4; i++) {
    resbuff[i] = resname[i];
  }
  lres = 4;
  lunit = 4;
  for (i = 3; i >= 0; i--) {
    if (resbuff[i] == ' ') {
      lres = i;
    }
    if (unitbuff[i] == ' ') {
      lunit = i;
    }
  }
  if (strcmp(myfit->qlibmatch, "PERMISSIVE") == 0) {
    for (i = 0; i <= lunit - lres; i++) {
      if (strncmp(&unitname[i], resname, lres) == 0) {
        return 1;
      }
    }
  }
  else {
    if (lunit == lres && strncmp(unitname, resname, lunit) == 0) {
      return 1;
    }
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
// OutputAmbLib: function to print a new Amber library file containing the
//               fitted charges.
//
// Arguments:
//   myfit:      the fitting control data
//   tj:         trajectory control information (needed for output file name)
//   outname:    name of the library output to write
//   tmpname:    name of the template library file to base this on
//   qres:       matrix filled with putative solutions to the RESP fit  
//   qtype:      type of charges to log (0 = REsP / original IPolQ, 1 = vacuum
//               phase IPolQ, 2 = new IPolQ, 3 = hyperpolarized phase IPolQ)
//   Atm2Col:    correspondence of atoms in each system and columns in the
//               fitting matrix                                     
//   outp:       output report file
//-----------------------------------------------------------------------------
static void OutputAmbLib(fset *myfit, trajcon *tj, char* outname,
                         char* tmpname, dmat *qres, int qtype, imat *Atm2Col,
                         FILE *outp)
{
  int i, j, k, m, n, p, ajn, anchored, resid, ierr, atommatched, resmatched;
  double* unitcharges;
  char line[MAXLINE], aname[8], atype[8];
  FILE *foutp;
  cmat tmplt, lwords, ewords, unitatoms, unittypes;
  prmtop *tp;  

  // Read the original, template library file
  tmplt = Ascii2Mem(tmpname, 384, 16, "Amber library template not found.");
  fprintf(outp, " Notes printing %s, based on\n template library file %s:\n",
          outname, tmpname);

  // Step through the template line by line, and edit
  for (i = 0; i < tmplt.row; i++) {
    if (strncmp(tmplt.map[i], "!entry.", 7) == 0) {
      lwords = ParseWords(tmplt.map[i]);
      for (j = 0; j < strlen(lwords.map[0]); j++) {
        if (lwords.map[0][j] == '.') {
          lwords.map[0][j] = ' ';
        }
      }
      ewords = ParseWords(lwords.map[0]);
      DestroyCmat(&lwords);
      if (strcmp(ewords.map[2], "unit") == 0 &&
          strcmp(ewords.map[3], "atoms") == 0) {

        // Seek out all atoms in this unit and match them to a residue
        // within one of the topologies
        j = i+1;
        k = 0;
        while (strncmp(tmplt.map[j], "!entry.", 7) != 0) {
          k++;
          j++;
        }
        unitatoms = CreateCmat(k, 5);
        unittypes = CreateCmat(k, 5);
        unitcharges = (double*)malloc(k * sizeof(double));
        for (j = i + 1; j <= i + k; j++) {
          lwords = ParseWords(tmplt.map[j]);
          if (WordIsAtomType(lwords.map[1]) == 0) {
            printf("OutputAmbLib >> Error.  %s (atom %s in unit %s) is not "
                   "recognized as an\nOutputAmbLib >> Amber atom type.\n",
                   lwords.map[1], lwords.map[0], ewords.map[1]);
            ierr = 1;
            return;
          }
          Type4Char(unitatoms.map[j - (i+1)], lwords.map[0]);
          Type4Char(unittypes.map[j - (i+1)], lwords.map[1]);
          DestroyCmat(&lwords);
        }
        
        // Check that there are no duplicates within the atom names (this
        // should be unnecessary, but it bulletproofs what comes next)
        for (j = 0; j < unitatoms.row; j++) {
          for (k = j+1; k < unitatoms.row; k++) {
            if (str4cmp(unitatoms.map[j], unitatoms.map[k]) == 0) {
              printf("OutputAmbLib >> Error.  Duplicate atom names found in "
                     "template library\nOutputAmbLib >> file (%4.4s and "
                     "%4.4s)\n", unitatoms.map[j], unitatoms.map[k]);
              ierr = 1;
              return;
            }
          }
        }
        
        // With the atom names (and types) known, compare them to the
        // names and types present in each topology and try to find a match.
        for (j = 0; j < myfit->tpcount; j++) {
          tp = &myfit->TPbank[j];

          // Loop through all atoms of the topology and see if one
          // matches any atom of the set we are looking for.  Once
          // found, find the residue of this matching atom in the
          // topology, and then search that residue to verify that
          // all atoms in the set are contained in the residue, and
          // that the residue contains exactly as many atoms as the
          // set we are trying to match.
          resmatched = -1;
          for (k = 0; k < tp->natom; k++) {
            anchored = 0;
            for (m = 0; m < unitatoms.row; m++) {
              if (str4cmp(unitatoms.map[m], &tp->AtomNames[4*k]) == 0) {
                anchored = 1;
              }
            }
            if (anchored == 0) {
              continue;
            }
            resid = LocateResID(tp, k, 0, tp->nres);
            if (tp->ResLims[resid+1] - tp->ResLims[resid] != unitatoms.row) {
              continue;
            }

            // Match the residue and unit names, invoking a measure of
            // permissiveness if the user allows it.
            if (UnitResidueMatch(myfit, &tp->ResNames[4*resid],
                                 ewords.map[1]) == 0) {
              continue;
            }
            resmatched = resid;
            for (m = tp->ResLims[resid]; m < tp->ResLims[resid+1]; m++) {
              atommatched = 0;
              for (n = 0; n < unitatoms.row; n++) {
                if (str4cmp(unitatoms.map[n], &tp->AtomNames[4*m]) == 0) {
                  atommatched = 1;
                }
              }
              if (atommatched == 0) {
                resmatched = -1;
              }
            }
            if (resmatched >= 0) {

              // Check that the atom types are consistent, fill in the
              // correct atom types if not, and change the charges.
              for (m = 0; m < unitatoms.row; m++) {
                for (n = tp->ResLims[resmatched];
                     n < tp->ResLims[resmatched + 1]; n++) {
                  if (str4cmp(unitatoms.map[m], &tp->AtomNames[4*n]) == 0) {
                    if (str4cmp(unittypes.map[m], &tp->AtomTypes[4*n]) != 0) {
                      fprintf(outp, " - Atom type update in %s / %s: %4.4s "
                              "-> %4.4s\n", ewords.map[1], unitatoms.map[m],
                              unittypes.map[m], &tp->AtomTypes[4*n]);
                    }
                    for (p = 0; p < 4; p++) {
                      unittypes.map[m][p] = tp->AtomTypes[4*n + p];
                    }
                    if (qtype == 0 || qtype == 1 || qtype == 2) {
                      unitcharges[m] = qres->map[qtype][Atm2Col->map[j][n]];
                    }
                    else if (qtype == 3) {
                      ajn = Atm2Col->map[j][n];
                      unitcharges[m] = qres->map[1][ajn] +
                                       (2.0 * qres->map[3][ajn]);
                    }                      
                  }
                }
              }
              fprintf(outp, " - Printing unit %s with %9.5lf net charge.\n",
                      ewords.map[1], DSum(unitcharges, unitatoms.row));
              
              // Write the results back into the scanned character matrix for
              // the Amber library file.
              for (m = i + 1; m <= i + unitatoms.row; m++) {
                lwords = ParseWords(tmplt.map[m]);
                aname[0] = '"';
                for (n = 0; n < 4; n++) {
                  if (unitatoms.map[m - (i+1)][n] != ' ') {
                    aname[n + 1] = unitatoms.map[m - (i+1)][n];
                  }
                  else {
                    break;
                  }
                }
                aname[n + 1] = '"';
                aname[n+2] = '\0';
                atype[0] = '"';
                for (n = 0; n < 4; n++) {
                  if (unittypes.map[m - (i+1)][n] != ' ') {
                    atype[n + 1] = unittypes.map[m - (i+1)][n];
                  }
                  else {
                    break;
                  }
                }
                atype[n + 1] = '"';
                atype[n+2] = '\0';                
                sprintf(tmplt.map[m], "    %-6.6s  %-5.5s %3d %3d %8d %3d %3d "
                        "%8.5lf\n", aname, atype, atoi(lwords.map[2]),
                        atoi(lwords.map[3]), atoi(lwords.map[4]),
                        atoi(lwords.map[5]), atoi(lwords.map[6]),
                        unitcharges[m - (i+1)]);
                DestroyCmat(&lwords);
              }

              // Break if the residue has been matched
              break;
            }
          }

          // Break if the residue has already been matched
          if (resmatched >= 0) {
            break;
          }
        }
        DestroyCmat(&unitatoms);
        DestroyCmat(&unittypes);
        free(unitcharges);
      }
      DestroyCmat(&ewords);
    }
  }

  // Spacing in the output
  fprintf(outp, "\n");

  // Write the output based on the modified data in memory
  if ((foutp = fopen(outname, "w")) == NULL) {
    printf("OutputAmbLib >> Error.  Could not open %s for writing.\n",
           outname);
    exit(1);
  }
  for (i = 0; i < tmplt.row; i++) {
    fprintf(foutp, "%s", tmplt.map[i]);
  }
  fclose(foutp);
}

//-----------------------------------------------------------------------------
// OutputRESP: produce formatted output to present the results of the RESP
//             calculation.  If an extra points file was specified, then a new
//             extra points file with the appropriate names and charges
//             assigned to each point will be written.  A set of scaled charges
//             will be written in 5 x %16.8e format for inclusion in the
//             topology file.               
//
// Arguments:
//   myfit:      the fitting command data, used here to determine the number of
//               jackknife fits and their style               
//   tj:         trajectory control data, used here for directives on output
//               overwriting and input file reprinting           
//   qres:       matrix filled with putative solutions to the RESP fit  
//   allcrd:     coordinates of all fitting conformations               
//   Atm2Col:    correspondence of atoms in each system and columns in the
//               fitting matrix                                     
//   etimers:    the execution timings data                             
//-----------------------------------------------------------------------------
static void OutputRESP(fset *myfit, trajcon *tj, dmat *qres, dmat *allcrd,
                       imat *Atm2Col, execon *etimers)
{
  int i, j, k, nchg;
  double qval;
  FILE *outp;
  time_t ct;
  prmtop *tp;

  // Open the output file
  ct = time(&ct);
  outp = FOpenSafe(tj->outbase, tj->OverwriteOutput);
  fprintf(outp, "Run on %s", asctime(localtime(&ct)));

  // Reprint the input file
  HorizontalRule(outp, 0);
  fprintf(outp, "\nINPUT LINE TEXT:\n\n");
  PrintParagraph(tj->inpline, 79, NULL, outp);
  fprintf(outp, "\nINPUT FILE TEXT:\n\n");
  for (i = 0; i < tj->inptext.row; i++) {
    fprintf(outp, "%s", tj->inptext.map[i]);
  }
  HorizontalRule(outp, 1);

  // Print the best resulting charge set(s)
  HorizontalRule(outp, 0);
  fprintf(outp, "(1.) Overall accuracy of the derived charges.\n\n");
  if (myfit->model == 0) {
    fprintf(outp, " Correlation of fitted to original electrostatic "
            "potential:       %8.5lf\n", myfit->Qcorr[0]);
    fprintf(outp, " Error fitted versus original electrostatic potential "
            "(kcal/mol): %8.5lf\n", myfit->Qrmsd[0]);
    fprintf(outp, " Snap iterations required:                           "
            "                 %4d\n\n", myfit->SnapCount[0]);
    PrintVADesc(0, " ", 1, " ", 1, "The electrostatic potential is easier to "
                "fit for some systems and targets than others.  The following "
                "table gives values for individual systems.\n", 77, 0,outp);
  }
  else {
    PrintVADesc(0, " ", 1, " ", 1, "The electrostatic potential was fitted in "
                "two ways.  First, a traditional REsP fit was performed to "
                "obtain a set of IPolQ charges by fitting them directly to "
                "the average of the vacuum and condensed-phase potentials.  "
                "Next, a two-component REsP was performed to fit a set of "
                "charges to the vacuum phase potential only, restrained under "
                "the same conditions as the first set of IPolQ charges, and "
                "to simultaneously fit a second set of IPolQ charges composed "
                "of the vacuum charge set plus a minimal perturbation.  This "
                "second approach to fitting IPolQ charges produces a result "
                "which can be related to charges appropriate for potential "
                "energy surfaces of systems in vacuum, which suggests a means "
                "of obtaining an appropriate set of torsion potential Fourier "
                "terms for IPolQ charges.\n", 77, 0, outp);
    fprintf(outp, "                    Accuracy to Aggregate Target\n"
            "   Charge Set      Correlation   RMS Error   Snap\n"
            " -------------     -----------  -----------  ----\n"
            " IPolQ (orig)        %9.4lf  %11.4lf  %4d\n"
            " Vacuum              %9.4lf  %11.4lf  %4d\n"
            " IPolQ (pert)        %9.4lf  %11.4lf  %4d\n\n",
            myfit->Qcorr[0], myfit->Qrmsd[0], myfit->SnapCount[0],
            myfit->Qcorr[1], myfit->Qrmsd[1], myfit->SnapCount[1],
            myfit->Qcorr[2], myfit->Qrmsd[2], myfit->SnapCount[2]);
  }
  PrintVADesc(0, " ", 1, " ", 1, "The electrostatic potential is easier to "
              "fit for some systems and targets than others.  The following "
              "table gives values for individual systems.\n", 77, 0,outp);
  if (myfit->model == 1) {
    fprintf(outp, "                                   IPolQ (orig) "
            "     Vacuum      IPolQ (pert)\n");
  }
  fprintf(outp, " System                            Corr   RMSE  ");
  if (myfit->model == 0) {
    fprintf(outp, "\n");
  }
  else {
    fprintf(outp, "  Corr   RMSE    Corr   RMSE\n");
  }
  fprintf(outp, " ------------------------------   ------ ------  ");
  if (myfit->model == 0) {
    fprintf(outp, "\n");
  }
  else {
    fprintf(outp, "------ ------  ------ ------\n");
  }
  for (i = 0; i < myfit->tpcount; i++) {
    fprintf(outp, " %-30.30s   %6.3lf %6.2lf", myfit->TPbank[i].source,
            myfit->QScorr.map[0][i], myfit->QSrmsd.map[0][i]);
    if (myfit->model == 0) {
      fprintf(outp, "\n");
    }
    else {
      fprintf(outp, "  %6.3lf %6.2lf  %6.3lf %6.2lf\n",
              myfit->QScorr.map[1][i], myfit->QSrmsd.map[1][i],
              myfit->QScorr.map[2][i], myfit->QSrmsd.map[2][i]);
    }
  }

  HorizontalRule(outp, 1);
  HorizontalRule(outp, 0);
  fprintf(outp, "(2.) Charges on all atoms, proton units\n");
  for (i = 0; i < myfit->tpcount; i++) {
    if (myfit->model == 0) {
      fprintf(outp, "\n System %d: %s\n Atom  Residue    Charge  \n"
              " ----  -------  ----------\n", i+1,
              myfit->TPbank[i].source);
    }
    else {
      fprintf(outp, "\n System %d: %s\n"
              " Atom  Residue  IPolQ,orig   Vacuum Q    Delta Q    "
              "IPolQ,pert\n"
              " ----  -------  ----------  ----------  ----------  "
              "----------\n",
              i + 1, myfit->TPbank[i].source);
    }
    tp = &myfit->TPbank[i];
    for (j = 0; j < tp->natom; j++) {
      nchg = Atm2Col->map[i][j];
      k = LocateResID(tp, j, 0, tp->nres);
      if (myfit->model == 0) {
        fprintf(outp, " %.4s    %4.4s   %10.5lf\n", &tp->AtomNames[4*j],
                &tp->ResNames[4*k], qres->map[0][nchg]);
      }
      else {
        fprintf(outp, " %.4s    %4.4s   %10.5lf  %10.5lf  %10.5lf  %10.5lf\n",
                &tp->AtomNames[4*j], &tp->ResNames[4*k], qres->map[0][nchg],
                qres->map[1][nchg], qres->map[3][nchg], qres->map[2][nchg]);
      }
    }
  }
  HorizontalRule(outp, 1);

  // Print the best resulting charge set in prmtop format
  HorizontalRule(outp, 0);
  fprintf(outp, "(3.) Charges on all atoms, prmtop format (proton units "
          "scaled by 18.2223)\n");
  for (i = 0; i < myfit->tpcount; i++) {
    tp = &myfit->TPbank[i];
    fprintf(outp, "\n System %d: %s\n", i+1, myfit->TPbank[i].source);
    k = 0;
    for (j = 0; j < tp->natom; j++) {
      qval = (myfit->model == 0) ? qres->map[0][Atm2Col->map[i][j]] :
        qres->map[2][Atm2Col->map[i][j]];
      fprintf(outp, "%16.8e", qval*sqrt(BIOQ));
      k++;
      if (k == 5) {
        k = 0;
        fprintf(outp, "\n");
      }
    }
    if (k > 0) {
      fprintf(outp, "\n");
    }
  }
  HorizontalRule(outp, 1);

  // Compute and print dipole moments
  CalculateDipoles(myfit, qres, allcrd, Atm2Col, outp);

  // Print the library file
  if (myfit->namblib > 0) {
    HorizontalRule(outp, 0);
    if (myfit->namblib == 1) {
      fprintf(outp, "(5.) Amber-format Library File\n\n");
    }
    else {
      fprintf(outp, "(5.) Amber-format Library Files\n\n");
    }
    for (i = 0; i < myfit->namblib; i++) {
      OutputAmbLib(myfit, tj, myfit->qlibout.map[i], myfit->qlibtmp.map[i],
                   qres, myfit->qlibstyle.map[0][i], Atm2Col, outp);
    }
    HorizontalRule(outp, 1);
  }
  
  // Print timings data
  etimers->Write = mdgxStopTimer(etimers);
  QFitTimingsData(etimers, outp);
  fclose(outp);

  // Print an extra points file
  if (myfit->epext[0] != '\0') {
    for (i = 0; i < myfit->tpcount; i++) {
      PrintEPRules(myfit, i, tj);
    }
  }

  // Print the histogram of fitting points
  if (myfit->histfile[0] != '\0') {
    PrintRespHistogram(myfit, tj);
  }
}

//-----------------------------------------------------------------------------
// DestroyFitData: cleanup for all data related to charge fitting.      
//                                                                      
// Arguments:                                                           
//   myfit:  the fitting control data                                   
//-----------------------------------------------------------------------------
static void DestroyFitData(fset *myfit)
{
  int i;

  // Simple frees
  free(myfit->fitpthist);
  free(myfit->tpidx);
  free(myfit->totalq);
  free(myfit->wt);
  free(myfit->nsamp);
  free(myfit->FPtOrigins);
  
  // Constraint structs
  for (i = 0; i < myfit->nqeq; i++) {
    free(myfit->qeq[i].atoms);
    free(myfit->qeq[i].maskstr);
  }
  free(myfit->qeq);
  for (i = 0; i < myfit->nqmin; i++) {
    free(myfit->qmin[i].maskstr);
  }
  free(myfit->qmin);
  for (i = 0; i < myfit->nqsum; i++) {
    free(myfit->qsum[i].maskstr);
  }
  free(myfit->qsum);

  // File names
  DestroyCmat(&myfit->gname);
  DestroyCmat(&myfit->auxgname);
  DestroyCmat(&myfit->tpname);

  // Topologies
  for (i = 0; i < myfit->tpcount; i++) {
    FreeTopology(&myfit->TPbank[i]);
  }
  free(myfit->TPbank);
  free(myfit->eprule);

  // Libraries
  if (myfit->namblib > 0) {
    DestroyCmat(&myfit->qlibout);
    DestroyCmat(&myfit->qlibtmp);
    DestroyImat(&myfit->qlibstyle);
  }
}

//-----------------------------------------------------------------------------
// FitCharges: the main function for parameter optimization.            
//                                                                      
// Arguments:                                                           
//   myfit:    the fitting control data                                 
//   tj:       trajectory control information                           
//   etimers:  the execution timings                                    
//-----------------------------------------------------------------------------
void FitCharges(fset *myfit, trajcon *tj, execon *etimers)
{
  int i, maxatm;
  imat Atm2Col;
  dmat cnstmat, fitmat, qres;
  dmat allcrd;

  // Match up topology names with grids
  AssignGridTopologies(myfit, tj);

  // Match up atoms in all topologies with matrix columns
  Atm2Col = ColumnsForAtoms(myfit);
  CheckDegeneracy(myfit, &Atm2Col);
  
  // Allocate space for coordinate buffer
  maxatm = 0;
  for (i = 0; i < myfit->tpcount; i++) {
    if (myfit->TPbank[i].natom > maxatm) {
      maxatm = myfit->TPbank[i].natom;
    }
  }
  allcrd = CreateDmat(myfit->ngrd, 3*maxatm, 0);
  etimers->Setup += mdgxStopTimer(etimers);

  // Create the matrix of all fitting points
  fitmat = MakeFittingMatrix(myfit, &Atm2Col, &allcrd, tj, etimers);

  // Create the matrix of constraints
  cnstmat = MakeCnstMatrix(myfit, &Atm2Col, &fitmat);
  etimers->Setup += mdgxStopTimer(etimers);

  // Solve the linear least squares problem
  qres = SolveRESP(&fitmat, &cnstmat, myfit, &Atm2Col, etimers);
  etimers->Setup += mdgxStopTimer(etimers);

  // Re-evaluate each system in light of the fitted charges
  //  SystemReeval(myfit, tj, &qres, &Atm2Col, etimers);

  // Output results
  OutputRESP(myfit, tj, &qres, &allcrd, &Atm2Col, etimers);
  
  // Free allocated memory
  DestroyDmat(&fitmat);
  DestroyDmat(&cnstmat);
  DestroyDmat(&allcrd);
  DestroyImat(&Atm2Col);
  DestroyFitData(myfit);

  // Exit!  We are done.
  exit(1);
}
