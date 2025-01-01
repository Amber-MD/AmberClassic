#include <math.h>
#include <stdlib.h>
#include "Grid.h"
#include "BSpline.h"
#include "mdgxVector.h"
#include "pmeRecip.h"
#include "mleRecip.h"
#include "Timings.h"
#include "Matrix.h"
#include "fftw3.h"

//-----------------------------------------------------------------------------
// GrdCofR: create a set of B-Spline coefficients for interpolation between
//          two grids.                                          
//                                                                      
// Arguments:                                                           
//   rcinp:   reciprocal space control information                      
//   lcon:    the reciprocal space mesh level that this map describes   
//   ndim:    each map describes interpolation in one dimension (Y or Z for
//            slab MLE); ndim = 0 implies interpolation in X, ndim = 1 implies
//            interpolation in Y, and ndim = 2 implies interpolation in Z
//-----------------------------------------------------------------------------
static g2gmap GrdCofR(reccon *rcinp, int lcon, int ndim)
{
  int i, j, ggordr;
  g2gmap G;
  bmap mml;
  coord crl;
  reccon rcnew;

  // Allocate tables
  G.ng = rcinp->ng[ndim];
  ggordr = rcinp->ggordr;
  G.ggordr = ggordr;
  G.s = CreateDmat(G.ng, ggordr, 0);
  G.m = CreateImat(G.ng, ggordr);

  // Create a temporary spline map, as if this were
  // a 3D particle <--> mesh interpolation problem.
  // The system is said to be orthorhombic, because
  // grid-to-grid interpolation occurs between to  
  // sets of points that are already in the same   
  // vector space.                                 
  rcnew.ng = (int*)malloc(3*sizeof(int));
  for (i = 0; i < 3; i++) {
    rcnew.ng[i] = rcinp->ng[3*lcon+i];
    rcnew.ordr[i] = ggordr;
    crl.gdim[i] = rcinp->ng[i];
  }
  crl.natom = G.ng;
  crl.atmid = (int*)calloc(G.ng, sizeof(int));
  crl.loc = (double*)calloc(3*G.ng, sizeof(double));
  for (i = 0; i < rcinp->ng[ndim]; i++) {
    crl.loc[3*i+ndim] = i;
  }
  crl.isortho = 1;
  mml = CreateBmap(&rcnew, G.ng);
  SplCoeff(&crl, &mml, &rcnew);

  // Copy the particle <--> mesh maps into the
  // grid <--> grid map structure             
  for (i = 0; i < G.ng; i++) {
    for (j = 0; j < ggordr; j++) {
      G.s.map[i][j] = (ndim == 0) ? mml.xcof[i*ggordr+j].s :
	(ndim == 1) ? mml.ycof[i*ggordr+j].s : mml.zcof[i*ggordr+j].s;
      G.m.map[i][j] = (ndim == 0) ? mml.xcof[i*ggordr+j].m :
        (ndim == 1) ? mml.ycof[i*ggordr+j].m : mml.zcof[i*ggordr+j].m;
    }
  }

  // Free allocated memory
  DestroyBmap(&mml);

  return G;
}

//-----------------------------------------------------------------------------
// SpreadBook: a function for interpolating a coarse book from a fine book,
//             for later use in Multi-Level Ewald calculations.   
//-----------------------------------------------------------------------------
void SpreadBook(dbook Qc, dbook Q, g2gmap SPr, g2gmap SPc)
{
  int i, j, k, m, ggordr;
  int* tmap;
  double qpt;
  double *dtmp, *dctmp, *tcof;
  double* cgbuff;
  double **dtm2p, **dctm2p;

  // Check input
  ggordr = SPr.ggordr;

  // Allocate buffer memory
  cgbuff = (double*)malloc(Qc.col*sizeof(double));

  // Wipe clean the book that is to be interpolated
  if (Qc.pfft == 1) {
    SetDVec(Qc.data, 2*Qc.pag*Qc.row*(Qc.col/2+1), 0.0);
  }
  else {
    SetDVec(Qc.data, Qc.pag*Qc.row*Qc.col, 0.0);
  }

  // Loop over all pages
  for (i = 0; i < Q.pag; i++) {
    dtm2p = Q.map[i];
    dctm2p = Qc.map[i];

    // Loop over all rows
    for (j = 0; j < Q.row; j++) {
      dtmp = dtm2p[j];

      // First, spread along the columns
      SetDVec(cgbuff, Qc.col, 0.0);
      for (k = 0; k < Q.col; k++) {
        qpt = dtmp[k];
        tmap = SPc.m.map[k];
        tcof = SPc.s.map[k];
        for (m = 0; m < ggordr; m++) {
          cgbuff[tmap[m]] += qpt*tcof[m];
        }
      }

      // Now, spread along the rows
      tmap = SPr.m.map[j];
      tcof = SPr.s.map[j];
      for (m = 0; m < ggordr; m++) {
        dctmp = dctm2p[tmap[m]];
        qpt = tcof[m];
        for (k = 0; k < Qc.col; k++) {
          dctmp[k] += cgbuff[k]*qpt;
        }
      }
    }
  }

  // Free allocated memory
  free(cgbuff);
}

//-----------------------------------------------------------------------------
// IntrpBook: interpolate a fine book from a coarse book, for use in
//            Multi-Level Ewald calculations.  This function will NOT zero the
//            book Q, so that a fine book can be accumulated from many coarse
//            books.                                   
//
// Arguments:
//   Q:    the fine book to be interpolated / accumulated               
//   Qc:   the coarse book                                              
//   SPr:  interpolation coefficients for rows                          
//   SPc:  interpolation coefficients for columns                       
//-----------------------------------------------------------------------------
void IntrpBook(dbook Q, dbook Qc, g2gmap SPr, g2gmap SPc)
{
  int i, j, k, m, ggordr;
  int* tmap;
  double qpt;
  double *dtmp, *dctmp, *tcof;
  double* cgbuff;
  double **dtm2p, **dctm2p;

  // Check input
  ggordr = SPr.ggordr;

  // Allocate buffer memory
  cgbuff = (double*)malloc(Qc.col*sizeof(double));

  // Loop over all pages
  for (i = 0; i < Q.pag; i++) {
    dtm2p = Q.map[i];
    dctm2p = Qc.map[i];

    // Loop over all rows
    for (j = 0; j < Q.row; j++) {
      dtmp = dtm2p[j];

      // Sum the rows
      tmap = SPr.m.map[j];
      tcof = SPr.s.map[j];
      SetDVec(cgbuff, Qc.col, 0.0);
      for (m = 0; m < ggordr; m++) {
        dctmp = dctm2p[tmap[m]];
        qpt = tcof[m];
        for (k = 0; k < Qc.col; k++) {
          cgbuff[k] += dctmp[k]*qpt;
        }
      }

      // Now, interpolate the fine grid points ****/
      for (k = 0; k < Q.col; k++) {
        tmap = SPc.m.map[k];
        tcof = SPc.s.map[k];
        qpt = 0.0;
        for (m = 0; m < ggordr; m++) {
          qpt += cgbuff[tmap[m]]*tcof[m];
        }
        dtmp[k] += qpt;
      }
    }
  }

  // Free allocated memory
  free(cgbuff);
}

//-----------------------------------------------------------------------------
// MakeRealBCMesh: pre-compute the BC mesh in real space for performing MLE
//                 calculations.                                    
//
// Arguments:
//   rcinp:  reciprcoal space control data                              
//   crd:    coordinates (for simulation box information)               
//   PPK:    the reciprocal space pair potential kit                    
//-----------------------------------------------------------------------------
static dbook MakeRealBCMesh(reccon *rcinp, coord *crd, prmtop *tp)
{
  dbook Urec;
  bckit PPk;
  execon etimers;

  // Create a charge mesh with a single +1 charge at (0,0,0)
  Urec = CreateDbook(rcinp->ng[0], rcinp->ng[1], rcinp->ng[2], 1);

  // Create a pair potential kit
  PPk = CreateBCKit(rcinp, &Urec, crd, tp, FFTW_MEASURE);
  Urec.map[0][0][0] = 1.0;

  // Convolute this mesh with the reciprocal space pair potential
  // to obtain the reciprocal space pair potential in real space 
  InitExecon(&etimers);

  // DEBUG IDEA: copy this function to a static function in this  
  // library, then find the lines which cause the slew of valgrind
  // errors after the function is called.                         
  ConvQBC(rcinp, crd, &Urec, &PPk, &etimers);

  // Free allocated memory
  DestroyBCKit(&PPk);

  return Urec;
}

//-----------------------------------------------------------------------------
// MeshLevelLimits: determine which layers of the reciprocal space pair 
//                  potential this mesh is responsible for.             
//                                                                      
// Arguments:                                                           
//   rcinp:  reciprcoal space control data                              
//   lcon:   the level at which to cut the mesh (0 is the first level of
//           coarsening, typically by a factor of 2; 1 is the next      
//           level, typically by a factor of 4; and so on)              
//-----------------------------------------------------------------------------
void MeshLevelLimits(reccon *rcinp, int lcon, int *llim1, int *hlim1,
		     int *llim2, int *hlim2)
{
  int i;
  int* ng;

  ng = rcinp->ng;
  *llim1 = (lcon == 0) ? 0 : 1;
  *hlim1 = 1;
  for (i = 0; i < lcon; i++) {
    *llim1 += rcinp->PadYZ[i];
    *hlim1 += rcinp->PadYZ[i];
  }
  *hlim1 += rcinp->PadYZ[lcon];
  if (lcon == rcinp->nlev - 1) {
    *hlim1 = rcinp->ng[0] - *llim1 + 1;
    *llim2 = 0;
    *hlim2 = 0;
  }
  else {
    *llim2 = rcinp->ng[0];
    for (i = 0; i <= lcon; i++) {
      *llim2 -= rcinp->PadYZ[i];
    }
    *hlim2 = (lcon == 0) ? ng[0] : ng[0] - *llim1 + 1;
  }
}

//-----------------------------------------------------------------------------
// CutBCMeshYZ: cut the reciprocal space pair potential along the YZ plane,
//              tangential to the source.                        
//                                                                      
// Arguments:                                                           
//   rcinp:  reciprcoal space control data                              
//   Urec:   the original reciprocal space pair potential               
//   lcon:   the level at which to cut the mesh (0 is the first level of
//           coarsening, typically by a factor of 2; 1 is the next      
//           level, typically by a factor of 4; and so on)              
//-----------------------------------------------------------------------------
static dbook CutBCMeshYZ(reccon *rcinp, dbook *Urec, int lcon)
{
  int i, j, k, hj, hk, llim1, llim2, hlim1, hlim2;
  int *ng;
  double bfac, dj, dk;
  double *rectmp, *cpytmp;
  dbook Ucpy, Lrec;
  fftw_plan forwp;

  // Extract the relevant part of the reciprocal space pair potential
  ng = rcinp->ng;
  Ucpy = CreateDbook(ng[0], ng[1], ng[2], 0); 

  // Determine the layers of padding that have
  // already been assigned to other meshes    
  MeshLevelLimits(rcinp, lcon, &llim1, &hlim1, &llim2, &hlim2);

  // Copy the layers of padding
  for (i = llim1; i < hlim1; i++) {
    for (j = 0; j < ng[1]; j++) {
      rectmp = Urec->map[i][j];
      cpytmp = Ucpy.map[i][j];
      for (k = 0; k < ng[2]; k++) {
	cpytmp[k] = rectmp[k];
      }
    }
  }
  for (i = llim2; i < hlim2; i++) {
    for (j = 0; j < ng[1]; j++) {
      rectmp = Urec->map[i][j];
      cpytmp = Ucpy.map[i][j];
      for (k = 0; k < ng[2]; k++) {
	cpytmp[k] = rectmp[k];
      }
    }
  }

  // Ucpy now contains the part of the potential that is of interest. 
  // Ucpy may need to be coarsened in order to perform MLE.           
  Lrec = CreateDbook(ng[3*lcon], ng[3*lcon+1], ng[3*lcon+2], 1);
  forwp = fftw_plan_dft_r2c_3d(Lrec.pag, Lrec.row, Lrec.col,
			       Lrec.data, Lrec.fdata, FFTW_ESTIMATE);

  // cfac must be an integer, or interpolation
  // will fail to produce an accurate result. 
  if (ng[1] % ng[3*lcon+1] != 0 || ng[2] % ng[3*lcon+2] != 0) {
    printf("CutBCMeshXY >> Error.  Higher level mesh dimensions must be "
	   "factors\nCutBCMeshXY >> of the original mesh dimensions.\n");
    exit(1);
  }
  for (i = 0; i < Ucpy.pag; i++) {
    hj = 0;
    const double cjfac = ((double)ng[1])/((double)ng[3*lcon+1]);
    for (dj = 1.0e-6; dj < Ucpy.row; dj += cjfac) {
      j = dj;
      rectmp = Ucpy.map[i][j];
      cpytmp = Lrec.map[i][hj];
      hk = 0;
      const double ckfac = ((double)ng[2])/((double)ng[3*lcon+2]);
      for (dk = 1.0e-6; dk < Ucpy.col; dk += ckfac) {
	k = dk;
	cpytmp[hk] = rectmp[k];
	hk++;
      }
      hj++;
    }
  }
  fftw_execute(forwp);
  fftw_destroy_plan(forwp);

  // Normalization of the FFT by the number of elements
  bfac = 1.0/(Lrec.pag*Lrec.row*Lrec.col);
  const int nelem = 2 * Lrec.pag * Lrec.row * (Lrec.col/2+1);
  rectmp = Lrec.data;
  for (i = 0; i < nelem; i++) {
    rectmp[i] *= bfac;
  }

  // If Urec is coarsened, fold in the B-spline reciprocal
  // space factors.  Otherwise, return.                   
  if (lcon == 0) {
    return Lrec;
  }

  double* By;
  double* Bz;
  By = LoadPrefac(rcinp->ggordr, Lrec.row);
  Bz = LoadPrefac(rcinp->ggordr, Lrec.col);
  for (i = 0; i < Lrec.pag; i++) {
    for (j = 0; j < Lrec.row; j++) {
      rectmp = Lrec.map[i][j];
      for (k = 0; k < Lrec.col/2 + 1; k++) {
	bfac = By[j]*Bz[k];
	rectmp[2*k] *= bfac;
	rectmp[2*k+1] *= bfac;
      }
    }
  }
  free(By);
  free(Bz);

  return Lrec;
}

//-----------------------------------------------------------------------------
// PrepMLE: prepare for MLE by allocating grid-to-grid interpolation    
//          coefficients, making the reciprocal space pair potential, and
//          cutting the pair potential into "short" and "long" ranged
//          components.                                          
//                                                                      
// Arguments:                                                           
//   rcinp:  reciprocal space control data                              
//   crd:    system coordinates (needed for box dimensions)             
//   tp:     system topology                                            
//-----------------------------------------------------------------------------
void PrepMLE(reccon *rcinp, coord *crd, prmtop *tp)
{
  int i, cg0, cg1, cg2;
  int *ng;

  // Set / adjust the finest-resolution mesh dimensions,
  // then compute dimensions for coraser meshes         
  rcinp->ng = (int*)realloc(rcinp->ng, 3*rcinp->nlev*sizeof(int));
  ng = rcinp->ng;
  for (i = 0; i < rcinp->nlev; i++) {
    ng[3*i] = ng[0];
    ng[3*i+1] = ng[1]/rcinp->cfac[i] + 1.0e-8;
    ng[3*i+2] = ng[2]/rcinp->cfac[i] + 1.0e-8;
    SetMeshDims(ng, crd->gdim);
  }

  // Make maps of coefficients for every level, other than the
  // the finest level for which no coefficients are needed, of
  // the reciprocal space calculation                         
  rcinp->SPrv = (g2gmap*)malloc(rcinp->nlev*sizeof(g2gmap));
  rcinp->SPcv = (g2gmap*)malloc(rcinp->nlev*sizeof(g2gmap));
  for (i = 1; i < rcinp->nlev; i++) {
    rcinp->SPrv[i] = GrdCofR(rcinp, i, 1);
    rcinp->SPcv[i] = GrdCofR(rcinp, i, 2);
  }

  // Create charge maps and prepare for convolutions on those
  rcinp->forwplan = (fftw_plan*)malloc(rcinp->nlev*sizeof(fftw_plan));
  rcinp->backplan = (fftw_plan*)malloc(rcinp->nlev*sizeof(fftw_plan));
  rcinp->QL = (dbook*)malloc((rcinp->nlev+1)*sizeof(dbook));
  rcinp->QL[rcinp->nlev] = CreateDbook(ng[0], ng[1], ng[2], 0);
  for (i = 0; i < rcinp->nlev; i++) {
    cg0 = ng[3*i];
    cg1 = ng[3*i+1];
    cg2 = ng[3*i+2];
    rcinp->QL[i] = CreateDbook(cg0, cg1, cg2, 1);
    rcinp->forwplan[i] =
      fftw_plan_dft_r2c_3d(cg0, cg1, cg2, rcinp->QL[i].data,
			   rcinp->QL[i].fdata, FFTW_MEASURE);
    rcinp->backplan[i] =
      fftw_plan_dft_c2r_3d(cg0, cg1, cg2, rcinp->QL[i].fdata,
			   rcinp->QL[i].data, FFTW_MEASURE);
  }

  // Create the reciprocal space pair potential.  This potential  
  // applies to the original simulation box volume.  If necessary,
  // the potential can be rescaled to accommodate a new simulation
  // box volume.                                                  
  rcinp->Urec = (dbook*)malloc((rcinp->nlev+1)*sizeof(dbook));
  rcinp->Urec[rcinp->nlev] = MakeRealBCMesh(rcinp, crd, tp);
  for (i = 0; i < rcinp->nlev; i++) {
    rcinp->Urec[i] = CutBCMeshYZ(rcinp, &rcinp->Urec[rcinp->nlev], i);
  }
}

//-----------------------------------------------------------------------------
// mleConvQBC: this function convolutes the charge mesh Q with the reciprocal
//             space pair potential stored in rcinp to arrive at the reciprocal
//             space electrostatic potential.  Because the electrostatic
//             potential energy is obtained in real space, it does not require
//             additional commands inserted into the middle of the convolution
//             loop.  Therefore, this routine does not have a separate version
//             to produce the energy.
//
// Arguments:
//   rcinp:      reciprocal space input control data
//   sysUV:      energy and virial accumulators
//   PPk:        reciprocal space mesh-to-mesh pair potential toolkit
//   etimers:    timing data for the dynamics run
//-----------------------------------------------------------------------------
void mleConvQBC(reccon *rcinp, Energy *sysUV, bckit *PPk, execon *etimers)
{
  int i, j, k, npc;
  double qureal, quimag;
  double *dtmp, *dtm2p;
  fftw_complex *qtmp, *utmp;

  // If energy must be computed, retain the charge mesh
  if (sysUV->updateU == 1) {
    for (i = 0; i < rcinp->QL[0].pag; i++) {
      for (j = 0; j < rcinp->QL[0].row; j++) {
	dtmp = rcinp->QL[0].map[i][j];
	dtm2p = rcinp->QL[rcinp->nlev].map[i][j];
	for (k = 0; k < rcinp->QL[0].col; k++) {
	  dtm2p[k] = dtmp[k];
	}
      }
    }
  }
  etimers->nbCnv += mdgxStopTimer(etimers);

  // Spread charge to higher level meshes
  for (i = 1; i < rcinp->nlev; i++) {
    SpreadBook(rcinp->QL[i], rcinp->QL[0], rcinp->SPrv[i], rcinp->SPcv[i]);
  }
  etimers->nbMtM += mdgxStopTimer(etimers);

  // Convolute each charge mesh with the appropriate potential
  for (i = 0; i < rcinp->nlev; i++) {
    fftw_execute(rcinp->forwplan[i]);
    etimers->nbFFT += mdgxStopTimer(etimers);
    npc = rcinp->QL[i].pag * rcinp->QL[i].row * (rcinp->QL[i].col/2+1);
    qtmp = rcinp->QL[i].fdata;
    utmp = rcinp->Urec[i].fdata;
    for (j = 0; j < npc; j++) {
      qureal = qtmp[j][0]*utmp[j][0] - qtmp[j][1]*utmp[j][1];
      quimag = qtmp[j][0]*utmp[j][1] + qtmp[j][1]*utmp[j][0];
      qtmp[j][0] = qureal;
      qtmp[j][1] = quimag;
    }
    etimers->nbCnv += mdgxStopTimer(etimers);
    fftw_execute(rcinp->backplan[i]);
  }
  etimers->nbFFT += mdgxStopTimer(etimers);

  // Interpolate the fine electrostatic potential
  for (i = 1; i < rcinp->nlev; i++) {
    IntrpBook(rcinp->QL[0], rcinp->QL[i], rcinp->SPrv[i], rcinp->SPcv[i]);
  }
  etimers->nbMtM += mdgxStopTimer(etimers);

  // Compute the energy
  if (sysUV->updateU == 1) {
    qureal = 0.0;
    for (i = 0; i < rcinp->QL[0].pag; i++) {
      for (j = 0; j < rcinp->QL[0].row; j++) {
	dtmp = rcinp->QL[rcinp->nlev].map[i][j];
	dtm2p = rcinp->QL[0].map[i][j];
	for (k = 0; k < rcinp->QL[0].col; k++) {
	  qureal += dtmp[k]*dtm2p[k];
	}
      }
    }
    sysUV->relec = 0.5*qureal + PPk->SelfEcorr;
  }
  etimers->nbCnv += mdgxStopTimer(etimers);
}
