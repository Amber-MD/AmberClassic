#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef MPI
#include <mpi.h>
#endif

#include "Restraints.h"
#include "mdgxVector.h"
#include "Matrix.h"
#include "Grid.h"
#include "CrdManip.h"
#include "Parse.h"
#include "Topology.h"
#include "BSpline.h"

#include "CellManipDS.h"
#include "TrajectoryDS.h"

//-----------------------------------------------------------------------------
// ReadRestraintGrid: read a grid-based restraint file in XPLOR format. 
//                                                                      
// Arguments:                                                           
//   fname:   the name of the grid file                                 
//   crd:     the coordinates of the system (used for a sanity check, to ensure
//            that that the dimensions of the grid are in rough agreement with
//            the system itself)                         
//-----------------------------------------------------------------------------
fbook ReadRestraintGrid(char* source, coord *crd)
{
  int i, j, k;
  int glim[6];
  double gbox[6];
  char line[128];
  FILE *inp;
  dmat L;
  fbook rgrd;

  // Check for the existence of the grid
  if ((inp = fopen(source, "r")) == NULL) {
    printf("ReadRestraintGrid >> Error.  Grid file %s not found.\n", source);
    exit(1);
  }

  // Read the preamble; require that a single file contain the entire
  // grid and expect alpha, beta, and gamma angles in degrees.       
  fgets(line, 128, inp);
  fgets(line, 128, inp);
  sscanf(line, "%d", &j);
  for (i = 0; i < j; i++) {
    fgets(line, 128, inp);
  }
  fgets(line, 128, inp);
  sscanf(line, "%d%d%d%d%d%d%d%d%d", &rgrd.pag, &glim[0], &glim[1], &rgrd.row,
	 &glim[2], &glim[3], &rgrd.col, &glim[4], &glim[5]);
  if (glim[0] != 1 || glim[1] != rgrd.pag || glim[2] != 1 ||
      glim[3] != rgrd.row || glim[4] != 1 || glim[5] != rgrd.col) {
    printf("ReadRestraintGrid >> Error.  Grid file %s indices invalid.\n"
	   "ReadRestraintGrid >> [ %4d %4d ] x [ %4d %4d ] x [ %4d %4d ]\n",
	   source, glim[0], glim[1], glim[2], glim[3], glim[4], glim[5]);
    exit(1);
  }
  fgets(line, 128, inp);
  for (i = 0; i < 6; i++) {
    sscanf(&line[12*i], "%lf", &gbox[i]);
    if (i >= 3) {
      gbox[i] *= PI/180.0;
    }
  }
  for (i = 0; i < 6; i++) {
    if (fabs(crd->gdim[i]-gbox[i])/crd->gdim[i] > 1.0e-4) {
      printf("ReadRestraintGrid >> Error.  Grid dimensions do not match "
	     "coordinates.\nReadRestraintGrid >> %10.4lf (grid %d) versus "
	     "%10.4lf (coordinate %d).\n", gbox[i], i, crd->gdim[i], i);
      exit(1);
    }

    // Overwrite gbox indices once they've been checked
    gbox[i] = 0.0;
  }

  // Allocate grid
  L = CreateDmat(3, 3, 0);
  for (i = 0; i < 3; i++) {
    L.map[0][i] = crd->U.map[0][i] * rgrd.pag;
    L.map[1][i] = crd->U.map[1][i] * rgrd.row;
    L.map[2][i] = crd->U.map[2][i] * rgrd.col;
  }
  rgrd = CreateFbook(rgrd.pag, rgrd.row, rgrd.col, &L, gbox);

  // Read each section, with k varying slowest (ZXY in XPLOR format)
  fgets(line, 128, inp);
  if (strncmp(line, "ZYX", 3) != 0) {
    printf("ReadRestraintGrid >> Error.  Expected ZXY flag, got:\n%s", line);
    exit(1);
  }
  for (k = 0; k < glim[5]; k++) {
    fscanf(inp, "%d", &i);
    if (i != k+1) {
      printf("ReadRestraintGrid >> Error.  Expected section %d, got:\n%s",
	     k+1, line);
      exit(1);
    }
    for (i = 0; i < glim[3]; i++) {
      for (j = 0; j < glim[1]; j++) {
	fscanf(inp, "%E", &rgrd.map[j][i][k]);
      }
    }
  }

  // Read the final landmark
  fscanf(inp, "%d", &i);
  if (i != -9999) {
    printf("ReadRestraintGrid >> Error.  -9999 landmark expected, found:\n%s",
	   line);
    exit(1);
  }
  fclose(inp);

  // Free allocated memory
  DestroyDmat(&L);

  return rgrd;
}

//-----------------------------------------------------------------------------
// WriteRestraintGrid: write a grid-based restraint file (XPLOR format).
//                                                                      
// Arguments:                                                           
//   fname:   the name of the grid file                                 
//   L:       the name of the grid file to write                        
//   crd:     coordinates (for box information)                         
//   tj:      trajectory control (for output overwriting directives)    
//-----------------------------------------------------------------------------
void WriteRestraintGrid(char* fname, fbook *L, trajcon *tj)
{
  int h, i, j, k;
  double x, y, z, a, b, g, invz;
  double zres[3], zax[3];
  FILE *outp;

  // Determine grid dimensions.
  x = L->invU.map[0][0];
  if (fabs(L->invU.map[0][1]) < 1.0e-8) {
    g = 0.5*PI;
  }
  else {
    g = atan(L->invU.map[1][1] / L->invU.map[0][1]);
  }
  y = L->invU.map[0][1] / cos(g);
  zax[0] = 0.0;
  zax[1] = 0.0;
  zax[2] = 1.0;
  DMatVecMult(&L->invU, zax, zres);
  z = sqrt(zres[0]*zres[0] + zres[1]*zres[1] + zres[2]*zres[2]);
  b = acos(L->invU.map[0][2] / z);
  invz = 1.0/z;

  // It is a very poorly numerically conditioned
  // problem to pull the angle alpha (a) out of 
  // the inverse transformation matrix.  Use the
  // law of cosines instead.                    
  zres[0] *= invz;
  zres[1] *= invz;
  zres[2] *= invz;
  zres[0] -= 1.0;
  a = acos(0.5*(2.0 - zres[0]*zres[0] - zres[1]*zres[1] - zres[2]*zres[2]));

  outp = FOpenSafe(fname, tj->OverwriteOutput);
  fprintf(outp, "\n");
  fprintf(outp, "%8d\n", 4);
  fprintf(outp, "Written by mdgx\nThis density map is designed to guide "
	  "simulated trajectories into ensembles\nwhich better resemble the "
	  "crystallographic electron density.\nUse with electron density "
	  "rules in %s.\n", tj->Leash.GridDefsFile);
  fprintf(outp, "%8d%8d%8d%8d%8d%8d%8d%8d%8d\n", L->pag, 1, L->pag, L->row, 1,
	  L->row, L->col, 1, L->col);
  fprintf(outp, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\nZYX\n",
	  x, y, z, a*180.0/PI, b*180.0/PI, g*180.0/PI);
  for (k = 0; k < L->col; k++) {
    fprintf(outp, "%8d\n", k+1);
    h = 0;
    for (j = 0; j < L->row; j++) {
      for (i = 0; i < L->pag; i++) {
        fprintf(outp, "%12.5E", L->map[i][j][k]);
	h++;
	if (h == 6) {
	  h = 0;
	  fprintf(outp, "\n");
	}
      }
    }
    if (h != 0) {
      fprintf(outp, "\n");
    }
  }
  fprintf(outp, "-9999\n%12.4lf%12.4lf\n", 1.0, 0.5);
  fclose(outp);
}

//-----------------------------------------------------------------------------
// WriteRestraintGrid: write a grid-based restraint file (XPLOR format).
//                                                                      
// Arguments:                                                           
//   fname:   the name of the grid file                                 
//   L:       the name of the grid file to write                        
//   crd:     coordinates (for box information)                         
//   tj:      trajectory control (for output overwriting directives)    
//-----------------------------------------------------------------------------
void WriteRestraintGridDB(char* fname, dbook *L, coord *crd, trajcon *tj)
{
  int h, i, j, k;
  FILE *outp;

  outp = FOpenSafe(fname, tj->OverwriteOutput);
  fprintf(outp, "\n");
  fprintf(outp, "%8d\n", 4);
  fprintf(outp, "Written by mdgx\nThis density map is designed to guide "
	  "simulated trajectories into ensembles\nwhich better resemble the "
	  "crystallographic electron density.\nUse with electron density "
	  "rules in %s.\n", tj->Leash.GridDefsFile);
  fprintf(outp, "%8d%8d%8d%8d%8d%8d%8d%8d%8d\n", L->pag, 1, L->pag, L->row, 1,
	  L->row, L->col, 1, L->col);
  fprintf(outp, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\nZYX\n", crd->gdim[0],
	  crd->gdim[1], crd->gdim[2], crd->gdim[3]*180.0/PI,
	  crd->gdim[4]*180.0/PI, crd->gdim[5]*180.0/PI);
  for (k = 0; k < L->col; k++) {
    fprintf(outp, "%8d\n", k+1);
    h = 0;
    for (j = 0; j < L->row; j++) {
      for (i = 0; i < L->pag; i++) {
        fprintf(outp, "%12.5E", L->map[i][j][k]);
	h++;
	if (h == 6) {
	  h = 0;
	  fprintf(outp, "\n");
	}
      }
    }
    if (h != 0) {
      fprintf(outp, "\n");
    }
  }
  fprintf(outp, "-9999\n%12.4lf%12.4lf\n", 1.0, 0.5);
  fclose(outp);
}

//-----------------------------------------------------------------------------
// ReadGridDefinitions: parse a file containing the definitions of how atoms
//                      interact with the grid.                   
//                                                                      
// Arguments:                                                           
//   source:    file of definitions by which atoms cling to the grid    
//   tp:        the system topology; rules for atoms clinging to the    
//              grid will be written here just like EP frame rules      
//-----------------------------------------------------------------------------
void ReadGridDefinitions(char* source, prmtop *tp)
{
  int collect;
  double radius;
  double gss[5], wt[5];
  char* atommask;
  char* line;
  cmat lwords;
  FILE *inp;

  // Open the input file
  if ((inp = fopen(source, "r")) == NULL) {
    printf("ReadGridDefinitions >> Error.  Definitions file %s not found.\n",
           source);
    exit(1);
  }

  // Prepare to take in cling rules
  tp->nclingrule = 0;

  // Repeatedly scan for a &sticky namelist
  line = (char*)malloc(MAXLINE*sizeof(char));
  atommask = (char*)malloc(MAXLINE*sizeof(char));
  while (AdvanceToSegment(inp, "cling", 0) != 0) {
    collect = 1;
    while (collect == 1) {
      collect = ReadNamelistLine(line, &lwords, "ReadGridDefinitions", inp);
      if (collect == 0) {
        continue;
      }
      SeekString(lwords, atommask, "AtomMask", "atoms");
      SeekReal(lwords, &gss[0], "Gaussian1", "gss1");
      SeekReal(lwords, &gss[1], "Gaussian2", "gss2");
      SeekReal(lwords, &gss[2], "Gaussian3", "gss3");
      SeekReal(lwords, &gss[3], "Gaussian4", "gss4");
      SeekReal(lwords, &wt[0], "Weight1", "wt1");
      SeekReal(lwords, &wt[1], "Weight2", "wt2");
      SeekReal(lwords, &wt[2], "Weight3", "wt3");
      SeekReal(lwords, &wt[3], "Weight4", "wt4");
      SeekReal(lwords, &wt[4], "WeightF", "flat");
      SeekReal(lwords, &radius, "Extent", "rad");
    }

    // Encode this cling rule
    tp->nclingrule += 1;

    DestroyCmat(&lwords);
  }

  // Free allocated memory
  free(line);
  free(atommask);
}

//-----------------------------------------------------------------------------
// ImposeESymmetry: with the potential (or restraint potential) due to the
//                  asymmetric unit plotted, spread the potential throughout
//                  the unit cell according to symmetry operations.
//
// Arguments:                                                           
//   kgrd:     the original grid, which will be consumed after the new  
//             grid is fully computed and ready for return              
//   xfrm:     the symmetry operations                                  
//-----------------------------------------------------------------------------
fbook ImposeESymmetry(fbook *kgrd, symop *xfrm)
{
  int h, i, j, k;
  double prx, pry, prz, ptwt, ptwtk;
  double sr[3];
  double *spctmp, *Utmp;
  fbook sgrd;

  // Create a new grid to accumulate the symmetry-related potential
  sgrd = CreateFbook(kgrd->pag, kgrd->row, kgrd->col, &kgrd->L, kgrd->orig);

  // Loop over all grid points and project them onto the new grid
  ptwt = 1.0/xfrm->nt;
  for (i = 0; i < kgrd->pag; i++) {
    for (j = 0; j < kgrd->row; j++) {
      for (k = 0; k < kgrd->col; k++) {

	// Compute the original grid point location
	spctmp = kgrd->invL.data;
	prx = i*spctmp[0] + j*spctmp[1] + k*spctmp[2] + kgrd->orig[0];
	pry = i*spctmp[3] + j*spctmp[4] + k*spctmp[5] + kgrd->orig[1];
	prz = i*spctmp[6] + j*spctmp[7] + k*spctmp[8] + kgrd->orig[2];

	// Compute the symmetry-related points
	ptwtk = ptwt*kgrd->map[i][j][k];
	for (h = 0; h < xfrm->nt; h++) {
	  Utmp = xfrm->rmat[h].data;
	  sr[0] = prx*Utmp[0] + pry*Utmp[1] + prz*Utmp[2] + xfrm->tvec[3*h];
	  sr[1] = prx*Utmp[3] + pry*Utmp[4] + prz*Utmp[5] + xfrm->tvec[3*h+1];
	  sr[2] = prx*Utmp[6] + pry*Utmp[7] + prz*Utmp[8] + xfrm->tvec[3*h+2];

	  // For each symmetry-related point,
	  // compute the interpolation domain
	  TriInterp(sr, &ptwtk, &sgrd, 0, 1);
	}
      }
    }
  }

  // Destroy the original grid
  DestroyFbook(kgrd);

  return sgrd;
}

//-----------------------------------------------------------------------------
// PrivateGridOrigin: determine the origin and extent of a cell's private grid
//                    integer indices into the main grid.  
//                                                                      
// Arguments:                                                           
//   kgrd:    the main potential grid spanning the entire simulation box
//   C:       the cell of interest                                      
//   CG:      the cell grid                                             
//   idir:    the axis along which to find the origin                   
//   llim:    the lower grid limit (returned)                           
//   hlim:    the upper grid limit (returned)                           
//-----------------------------------------------------------------------------
static void PrivateGridOrigin(fbook *kgrd, cell *C, cellgrid *CG, int idir,
			      int *llim, int *hlim)
{
  double axfac, axfacp1;

  // Find the voxel containing the cell origin, then
  // find the voxel containing the cell's far limit.
  // Finally add 4 for fourth-order interpolation.  
  axfac = (double)C->gbin[idir]/CG->dbng[idir];
  axfacp1 = ((double)C->gbin[idir]+1.0)/CG->dbng[idir];
  if (idir == 0) {
    *llim = kgrd->pag * axfac;
    *hlim = kgrd->pag * axfacp1 + 4;
  }
  else if (idir == 1) {
    *llim = kgrd->row * axfac;
    *hlim = kgrd->row * axfacp1 + 4;
  }
  else {
    *llim = kgrd->col * axfac;
    *hlim = kgrd->col * axfacp1 + 4;
  }
  *hlim -= *llim;
}

//-----------------------------------------------------------------------------
// FillPrivateGrid: step over indices of the large grid spanning the entire
//                  simulation cell and take the piece that is important to
//                  this cell.  The piece may be formed by wrapping voxels in
//                  the large grid.                  
//                                                                      
// Arguments:                                                           
//   pgrd:     the private grid to fill; the origin of this grid is set with
//             values that, when converted to integers, produce the starting
//             indices for reading Lgrd                        
//   Lgrd:     the larger grid from which to extract a piece            
//   n[i,j,k]: the dimensions of the private grid piece                 
//   lstr:     the point on the larger grid on which to start copying   
//-----------------------------------------------------------------------------
void FillPrivateGrid(fbook *pgrd, fbook *Lgrd, int ni, int nj, int nk,
		     int* lstr)
{
  int i, j, k;
  int idxi;
  int idxj;
  int* idxk;
  float *ptmp, *ltmp;
  
  // Array to store the inner loop indices
  idxk = (int*)malloc(nk*sizeof(int));
  i = lstr[2];
  for (k = 0; k < nk; k++) {
    if (i < 0) {
      i += Lgrd->col;
    }
    if (i >= Lgrd->col) {
      i -= Lgrd->col;
    }
    idxk[k] = i;
    i++;
  }

  // Loop over all pgrd indices
  idxi = lstr[0];
  for (i = 0; i < ni; i++) {
    if (idxi < 0) {
      idxi += Lgrd->pag;
    }
    else if (idxi >= Lgrd->pag) {
      idxi -= Lgrd->pag;
    }
    idxj = lstr[1];
    for (j = 0; j < nj; j++) {
      if (idxj < 0) {
	idxj += Lgrd->row;
      }
      else if (idxj >= Lgrd->row) {
	idxj -= Lgrd->row;
      }
      ptmp = pgrd->map[i][j];
      ltmp = Lgrd->map[idxi][idxj];
      for (k = 0; k < nk; k++) {
	ptmp[k] = ltmp[idxk[k]];
      }
      idxj++;
    }
    idxi++;
  }

  // Free allocated memory
  free(idxk);
}

//-----------------------------------------------------------------------------
// CellPrivateGrids: take a big grid and give each of the simulation cells
//                   their own private piece of it.  In parallel mode, the
//                   master process will read a restraint grid and then send
//                   pieces of it as needed to each slave process through this
//                   function.                     
//                                                                      
// Arguments:                                                           
//   LV:     the array of large grids from which all private grids are derived
//   nLV:    the number of large grids to distribute over cells         
//   CG:     the cell grid to which all private grids are assigned      
//-----------------------------------------------------------------------------
void CellPrivateGrids(fbook* LV, int nLV, cellgrid *CG)
{
  int h, i, j;
  int iRorig[3], iExtent[3];
  double dRorig[3], nRorig[3];
  cell *C;
  fbook epbook;
#ifdef MPI
  MPI_Status stt;
#endif

  // Loop over all cells and allocate data
  for (i = 0; i < CG->ncell; i++) {
    C = &CG->data[i];
    if (C->CGRank == CG->tid) {
      C->nFscr = nLV;
      C->Fscr = (fbook*)malloc(nLV * sizeof(fbook));
    }
  }

  // Loop over all grids and distribute data
  for (h = 0; h < nLV; h++) {

    // Loop over all cells; create a subgrid    
    // for each cell and export it if necessary.
    for (i = 0; i < CG->ncell; i++) {
      C = &CG->data[i];

      // If this process owns the cell,
      // allocate memory for the grid. 
      for (j = 0; j < 3; j++) {
	PrivateGridOrigin(&LV[h], C, CG, j, &iRorig[j], &iExtent[j]);
	dRorig[j] = iRorig[j];
      }
      DMatVecMult(&LV[h].invL, dRorig, nRorig);
      if (C->CGRank == CG->tid) {
	C->Fscr[h] = CreateFbook(iExtent[0], iExtent[1], iExtent[2], &LV[h].L,
				 nRorig);
      }
      else if (CG->tid == 0) {
	epbook = CreateFbook(iExtent[0], iExtent[1], iExtent[2], &LV[h].L,
			     nRorig);
      }

      // If this is process owns the cell, fill the
      // private grid; otherwise post a receive.   
      if (CG->tid == 0) {
	if (C->CGRank == 0) {
	  FillPrivateGrid(&C->Fscr[h], &LV[h], iExtent[0], iExtent[1],
			  iExtent[2], iRorig);
	}
#ifdef MPI
	else {
	  FillPrivateGrid(&epbook, &LV[h], iExtent[0], iExtent[1], iExtent[2],
			  iRorig);
	  MPI_Send(epbook.data, iExtent[0]*iExtent[1]*iExtent[2], MPI_FLOAT,
		   C->CGRank, i, CG->dspcomm);
	}
#endif
      }
#ifdef MPI
      else {
	MPI_Recv(C->Fscr[h].data, iExtent[0]*iExtent[1]*iExtent[2], MPI_FLOAT,
		 0, i, CG->dspcomm, &stt);
      }
#endif
    }
  }
}

//-----------------------------------------------------------------------------
// ScoreOnGrid: pre-compute the score that an atom of a particular type 
//              (defined by a restraint grid cling rule) would obtain for being
//              at a particular grid point.  Copy the basic grid and compute
//              such scores for all points, storing the results in a new grid.
//              Return the new grid, a lookup table for that grid definition
//              that is suitable for B-Spline interpolation.
//                                                                      
// Arguments:
//   Lgrd    : the three-dimensional histogram
//   cr      : a restraint
//   tj      : trajectory input control data
//   crd     : system coordinates
//-----------------------------------------------------------------------------
fbook ScoreOnGrid(fbook *Lgrd, cling *cr, trajcon *tj, coord *crd)
{
  int h, i, j, k, ixgrfac, ii, jj, kk;
  double gamp, gsig, ampfac, expfac, r2, xgrfac, invxgrfac;
  double spln4ij, invcd, nnormfac, gtot, ngreal, ngimag;
  double gvec[3], gloc[3], spln4[4];
  double *dtmp, *d2tmp;
  float *ftmp;
  fftw_complex *ntmp, *gtmp;
  dbook G, Ngrd, Pgrd;
  fbook scrL;
  fftw_plan gforw, nforw, nback;

  // First, promote the grid to a higher resolution if requested
  ixgrfac = tj->Leash.XpandGrid;
  xgrfac = ixgrfac;
  invxgrfac = 1.0/xgrfac;
  Ngrd = PromoteFbook(Lgrd, tj->Leash.XpandGrid, 1);

  // Loop over each Gaussian and prepare to make a convolution
  G = CreateDbook(Ngrd.pag, Ngrd.row, Ngrd.col, 1);
  for (h = 0; h < cr->ngss; h++) {
    gamp = cr->gss[h].amp;
    gsig = cr->gss[h].sig;
    ampfac = 2.0*PI*gsig*gsig;
    ampfac *= ampfac*ampfac;
    ampfac = gamp/sqrt(ampfac);
    expfac = 0.5/(gsig*gsig);
    for (i = 0; i < G.pag; i++) {
      gvec[0] = (i <= G.pag/2) ? i*invxgrfac : (i-G.pag)*invxgrfac;
      for (j = 0; j < G.row; j++) {
	gvec[1] = (j <= G.row/2) ? j*invxgrfac : (j-G.row)*invxgrfac;
	dtmp = G.map[i][j];
	for (k = 0; k < G.col; k++) {
	  gvec[2] = (k <= G.col/2) ? k*invxgrfac : (k-G.col)*invxgrfac;
	  DMatVecMult(&Lgrd->invL, gvec, gloc);
	  r2 = gloc[0]*gloc[0] + gloc[1]*gloc[1] + gloc[2]*gloc[2];
	  dtmp[k] += ampfac*exp(-r2*expfac);
	}
      }
    }
  }

  // Normalize the Gaussian grid
  gtot = 0.0;
  for (i = 0; i < G.pag; i++) {
    for (j = 0; j < G.row; j++) {
      dtmp = G.map[i][j];
      for (k = 0; k < G.col; k++) {
	gtot += dtmp[k];
      }
    }
  }
  gtot = 1.0/gtot;
  for (i = 0; i < G.pag; i++) {
    for (j = 0; j < G.row; j++) {
      dtmp = G.map[i][j];
      for (k = 0; k < G.col; k++) {
        dtmp[k] *= gtot;
      }
    }
  }

  // Perform convolution of the Gaussian density
  // with the density-based restraint potential.
  gforw = fftw_plan_dft_r2c_3d(G.pag, G.row, G.col, G.data, G.fdata,
			       FFTW_ESTIMATE);
  nforw = fftw_plan_dft_r2c_3d(Ngrd.pag, Ngrd.row, Ngrd.col, Ngrd.data,
			       Ngrd.fdata, FFTW_ESTIMATE);
  nback = fftw_plan_dft_c2r_3d(Ngrd.pag, Ngrd.row, Ngrd.col, Ngrd.fdata,
			       Ngrd.data, FFTW_ESTIMATE);
  fftw_execute(gforw);
  fftw_execute(nforw);
  j = Ngrd.pag * Ngrd.row * (Ngrd.col/2+1);
  gtmp = G.fdata;
  ntmp = Ngrd.fdata;
  for (i = 0; i < j; i++) {
    ngreal = ntmp[i][0]*gtmp[i][0] - ntmp[i][1]*gtmp[i][1];
    ngimag = ntmp[i][0]*gtmp[i][1] + ntmp[i][1]*gtmp[i][0];
    ntmp[i][0] = ngreal;
    ntmp[i][1] = ngimag;
  }
  fftw_execute(nback);
  fftw_destroy_plan(gforw);
  fftw_destroy_plan(nforw);
  fftw_destroy_plan(nback);

  // Free allocated memory
  DestroyDbook(&G);

  // Write results to a new grid, the same size as the
  // original Lgrd.  No interpolation is done at this 
  // step, as the previous grid expansion was done    
  // simply to create a finer step size for numerical 
  // integration of the Gaussian electron density.    
  Pgrd = PromoteFbook(Lgrd, 1, 1);
  nnormfac = 1.0 / (Ngrd.pag * Ngrd.row * Ngrd.col);
  ii = 0;
  for (i = 0; i < Ngrd.pag; i+=ixgrfac) {
    jj = 0;
    for (j = 0; j < Ngrd.row; j+=ixgrfac) {
      kk = 0;
      dtmp = Pgrd.map[ii][jj];
      d2tmp = Ngrd.map[i][j];
      for (k = 0; k < Ngrd.col; k+=ixgrfac) {
	dtmp[kk] = d2tmp[k] * nnormfac;
	kk++;
      }
      jj++;
    }
    ii++;
  }

  // Free allocate memory
  DestroyDbook(&Ngrd);

  // De-convolute the new lookup table
  // (Pgrd) with a 4th-order B-spline.
  G = CreateDbook(Lgrd->pag, Lgrd->row, Lgrd->col, 1);
  spln4[0] = BSpln(0.0, 4);
  spln4[1] = BSpln(1.0, 4);
  spln4[2] = BSpln(2.0, 4);
  spln4[3] = BSpln(3.0, 4);
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      spln4ij = spln4[i]*spln4[j];
      for (k = 0; k < 4; k++) {
	G.map[i][j][k] = spln4ij*spln4[k];
      }
    }
  }
  gforw = fftw_plan_dft_r2c_3d(G.pag, G.row, G.col, G.data, G.fdata,
                               FFTW_ESTIMATE);
  nforw = fftw_plan_dft_r2c_3d(Pgrd.pag, Pgrd.row, Pgrd.col, Pgrd.data,
                               Pgrd.fdata, FFTW_ESTIMATE);
  nback = fftw_plan_dft_c2r_3d(Pgrd.pag, Pgrd.row, Pgrd.col, Pgrd.fdata,
                               Pgrd.data, FFTW_ESTIMATE);
  fftw_execute(gforw);
  fftw_execute(nforw);
  j = Pgrd.pag * Pgrd.row * (Pgrd.col/2+1);
  gtmp = G.fdata;
  ntmp = Pgrd.fdata;
  for (i = 0; i < j; i++) {
    invcd = 1.0 / (gtmp[i][0]*gtmp[i][0] + gtmp[i][1]*gtmp[i][1]);
    ngreal = (ntmp[i][0]*gtmp[i][0] + ntmp[i][1]*gtmp[i][1]) * invcd;
    ngimag = (ntmp[i][1]*gtmp[i][0] - ntmp[i][0]*gtmp[i][1]) * invcd;
    ntmp[i][0] = ngreal;
    ntmp[i][1] = ngimag;
  }
  fftw_execute(nback);
  fftw_destroy_plan(gforw);
  fftw_destroy_plan(nforw);
  fftw_destroy_plan(nback);

  // Free allocated memory
  DestroyDbook(&G);

  // Copy the deconvoluted potential into the fbook-format lookup
  scrL = CreateFbook(Lgrd->pag, Lgrd->row, Lgrd->col, &Lgrd->L, Lgrd->orig);
  nnormfac = 1.0 / (Pgrd.pag * Pgrd.row * Pgrd.col);
  for (i = 0; i < scrL.pag; i++) {
    for (j = 0; j < scrL.row; j++) {
      ftmp = scrL.map[i][j];
      dtmp = Pgrd.map[i][j];
      for (k = 0; k < scrL.col; k++) {
	ftmp[k] = dtmp[k]*nnormfac;
      }
    }
  }

  // Free allocated memory
  DestroyDbook(&Pgrd);

  return scrL;
}

//-----------------------------------------------------------------------------
// CellIntrpRstr: cell-based interpolation of grid restraint potentials.  This
//                function replicates a small amount of code from
//                CellIntrpFrc() in the ChargeMap library to explicitly specify
//                fourth-order interpolation (higher orders will not buy any
//                more accuracy in this case) and to operate on fbook structs
//                kept by each cell individually.      
//                                                                      
// Arguments:
//   C      : the cell on which to operate
//   tp     : the system topology
//-----------------------------------------------------------------------------
void CellIntrpRstr(cell *C, prmtop *tp)
{
  int h, i, j, natom;
  double fx, fy, fz, atmq, SPx, SPy, SPxy, SPdx, SPdxy, SPxdy, Zuk;
  float *utmp;
  float **utm2p;
  bcof *xcf, *ycf, *zcf;

  // Shortcuts
  natom = C->nr[0];

  // FIX ME!!! This code is totally not done.
  // Gotta figure out the units and info that
  // will be available in the cell struct,   
  // then come back to this.        FIX ME!!!
  for (h = 0; h < natom; h++) {
    xcf = &C->xcof[4*h];
    ycf = &C->ycof[4*h];
    zcf = &C->zcof[4*h];
    atmq = C->data[h].q;
    fx = 0.0;
    fy = 0.0;
    fz = 0.0;
    for (i = 0; i < 4; i++) {

      // FIX ME!!!  How will I get the grids packed into the cell?
      //      utm2p = U->map[xcf[i].m];
      SPx = atmq*xcf[i].s;
      SPdx = atmq*xcf[i].d;
      for (j = 0; j < 4; j++) {
	utmp = utm2p[ycf[j].m];
	SPy = ycf[j].s;
	SPxy = SPx*SPy;
	SPdxy = SPdx*SPy;
	SPxdy = SPx*ycf[j].d;
	Zuk = zcf[0].s*utmp[zcf[0].m] + zcf[1].s*utmp[zcf[1].m] +
	  zcf[2].s*utmp[zcf[2].m] + zcf[3].s*utmp[zcf[3].m];
	fx -= SPdxy*Zuk;
	fy -= SPxdy*Zuk;
	fz -= SPxy*(zcf[0].d*utmp[zcf[0].m] + zcf[1].d*utmp[zcf[1].m] +
		    zcf[2].d*utmp[zcf[2].m] + zcf[3].d*utmp[zcf[3].m]);

      }
    }
  }
}

//-----------------------------------------------------------------------------
// CreateBellyMask: create a belly mask for atoms that will be permitted to be
//                  mobile in the simulation.  This routine will also check to
//                  make sure that no constraints are evaluated if they involve
//                  the frozen atoms.         
//
// Arguments:                                                           
//   tj:    trajectory control data (contains restraint information)    
//   tp:    topology for this system                                    
//-----------------------------------------------------------------------------
void CreateBellyMask(trajcon *tj, prmtop *tp)
{
  int i, j, jlim, resi, resj, jatm, ismobile, grpsize;
  coord crd;
  char moveword[16];

  // If there is no mask, just set all    
  // massive atoms as mobile and bail out.
  if (tj->Leash.active == 0 || tj->Leash.usebelly == 0) {
    tp->MobileAtoms = (int*)malloc(tp->natom*sizeof(int));
    for (i = 0; i < tp->natom; i++) {
      tp->MobileAtoms[i] = (tp->Masses[i] > 1.0e-8) ? 1 : 2;
    }
    return;
  }

  // Create a fake set of coordinates
  crd = CreateCoord(tp->natom);

  // Evaluate the belly mask
  if (tj->Leash.BellyMask[0] != '\0') {
    tp->MobileAtoms = ParseAmbMask(tj->Leash.BellyMask, tp, &crd);
  }
  else {
    tp->MobileAtoms = ParseAmbMask(tj->Leash.FrozenMask, tp, &crd);
    for (i = 0; i < tp->natom; i++) {
      tp->MobileAtoms[i] = 1 - tp->MobileAtoms[i];
    }
  }

  // Mark any extra points with a 2 (if they would
  // otherwise be part of the mobile group).      
  for (i = 0; i < tp->natom; i++) {
    if (tp->Masses[i] < 1.0e-8) {
      tp->MobileAtoms[i] = (tp->MobileAtoms[i] == 1) ? 2 : 0;
    }
  }

  // Analyze each constraint group
  if (tp->rattle == 1 || tp->settle == 1) {
    for (i = 0; i < tp->natom; i++) {
      if (tp->SHL[i].exe == 1) {
	jlim = 2;
      }
      else if (tp->SHL[i].exe == 2) {
	grpsize = tp->SHL[i].blist[0];
        jlim = tp->SHL[i].blist[3*grpsize + 1];
      }
      else {
	continue;
      }

      // If we're still here, this atom controls a constraint group.
      // Make sure the entire group is frozen or mobile to match this  
      // controlling atom.  Note RATTLE-guided groups (as opposed to
      // SETTLE).
      ismobile = tp->MobileAtoms[i];
      if (ismobile == 0) {
	tp->SHL[i].exe = 0;
      }
      for (j = 0; j < jlim; j++) {
	if (tp->SHL[i].exe == 1) {
	  jatm = tp->SHL[i].blist[j];
	}
	else {
	  jatm = tp->SHL[i].blist[3*grpsize + 2 + j];
	}
	if (tp->MobileAtoms[jatm] != ismobile && tp->MobileAtoms[jatm] != 2) {
	  if (ismobile == 0) {
	    sprintf(moveword, "immobilized");
	  }
	  else {
	    sprintf(moveword, "mobile");
	  }
	  resi = LocateResID(tp, i, 0, tp->nres);
	  resj = LocateResID(tp, jatm, 0, tp->nres);
	  printf("CreateBellyMask >> Error.  Atom %d (residue %.4s %d %.4s)\n"
		 "CreateBellyMask >> controls a constraint group and is "
		 "%s, but\nCreateBellyMask >> atom %d (residue %.4s "
		 "%d %.4s), part of the\nCreateBellyMask >> constraint group, "
		 "is not %s.\n", i, &tp->ResNames[4*resi], resi,
		 &tp->AtomNames[4*i], moveword, jatm, &tp->ResNames[4*resj],
		 resj, &tp->AtomNames[4*jatm], moveword);
	  exit(1);
	}
      }
    }
  }

  // Free allocated memory
  DestroyCoord(&crd);
}
