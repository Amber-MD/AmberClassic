#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "CellManip.h"
#include "mdgxVector.h"
#include "Matrix.h"
#include "Grid.h"
#include "Nonbonded.h"
#include "CrdManip.h"
#include "VirtualSites.h"
#include "Macros.h"
#include "Debug.h"
#include "Timings.h"
#include "MPIMap.h"
#include "Topology.h"

#include "pmeDirectDS.h"
#include "pmeRecipDS.h"
#include "CellManipDS.h"
#include "CompFrcDS.h"
#include "TrajectoryDS.h"
#include "BSplineDS.h"

//-----------------------------------------------------------------------------
// CreateCell: this routine creates a cell, the basic unit of input for the
//             tower and plate or triple2 methods.                  
//                                                                      
// Arguments:                                                           
//   maxatom:    the maximum number of atoms that this cell is to handle
//   ordr:       the order of particle <-> mesh interpolation (needed to
//               allocate B-spline coefficient arrays, and stored in the cell
//               for later use)                                    
//-----------------------------------------------------------------------------
cell CreateCell(int maxatom, int* ordr)
{
  int i;
  cell C;

  // The amount of memory that must be allocated is roughly
  // 3kb per atom, plus some incidentals related to the    
  // cell itself and maps within the cell.                 
  C.maxatom = maxatom;
  C.pmordr[0] = ordr[0];
  C.pmordr[1] = ordr[1];
  C.pmordr[2] = ordr[2];
  C.nr = (int*)calloc(16, sizeof(int));
  C.nsr = (int*)calloc(104, sizeof(int));
  C.xcof = (bcof*)malloc(maxatom*MAX(ordr[0],4)*sizeof(bcof));
  C.ycof = (bcof*)malloc(maxatom*MAX(ordr[1],4)*sizeof(bcof));
  C.zcof = (bcof*)malloc(maxatom*MAX(ordr[2],4)*sizeof(bcof));
  C.atmscr = (atomc*)malloc(maxatom*sizeof(atomc));
  C.ordr = CreateImat(104, maxatom);
  C.supordr = CreateImat(2, maxatom);
  C.qIDbuff = (int*)malloc(maxatom*sizeof(int));
  C.ljIDbuff = (int*)malloc(maxatom*sizeof(int));
  C.qr2buff = (rngbuff*)malloc(maxatom*sizeof(rngbuff));
  C.ljr2buff = (rngbuff*)malloc(maxatom*sizeof(rngbuff));
  C.data = (atomc*)malloc(8*maxatom*sizeof(atomc));
  C.map = (atomc**)malloc(8*sizeof(atomc*));
  for (i = 0; i < 8; i++) {
    C.map[i] = &C.data[i*maxatom];
  }
  C.import = (atomb*)malloc(4*maxatom*sizeof(atomb));
  C.pexport = (atomb*)malloc(4*maxatom*sizeof(atomb));
  C.Vimport = (atombv*)malloc(8*maxatom*sizeof(atombv));
  C.Vexport = (atombv*)malloc(8*maxatom*sizeof(atombv));
  C.Ximport = (atombx*)malloc(8*maxatom*sizeof(atombx));
  C.Xexport = (atombx*)malloc(8*maxatom*sizeof(atombx));

  // Grid-based restraints are not turned on by default
  C.nFscr = 0;

  return C;
}

//-----------------------------------------------------------------------------
// HessianNorms: this routine computes the distances between box faces given
//               the inverse box transformation matrix.           
//                                                                      
// Arguments:                                                           
//   invU:    the inverse box transformation matrix                     
//   cdepth:  upon return, stores the distances between box faces       
//-----------------------------------------------------------------------------
void HessianNorms(dmat *invU, double* cdepth)
{
  int i;
  double xyz0[3], thx[3], thy[3], thz[3], ic1[3], ic2[3], ic3[3];

  xyz0[0] = invU->data[0] + invU->data[1] + invU->data[2];
  xyz0[1] = invU->data[4] + invU->data[5];
  xyz0[2] = invU->data[8];
  for (i = 0; i < 3; i++) {
    ic1[i] = invU->map[i][0];
    ic2[i] = invU->map[i][1];
    ic3[i] = invU->map[i][2];
  }
  CrossP(ic2, ic3, thx);
  CrossP(ic1, ic3, thy);
  CrossP(ic1, ic2, thz);
  Normalize(thx, 3);
  Normalize(thy, 3);
  Normalize(thz, 3);
  cdepth[0] = fabs(DotP(thx, xyz0, 3));
  cdepth[1] = fabs(DotP(thy, xyz0, 3));
  cdepth[2] = fabs(DotP(thz, xyz0, 3));
}

//-----------------------------------------------------------------------------
// TakeCellGridDims: compute the cell grid dimensions based on direct space
//                   cutoffs and the box dimensions.              
//                                                                      
// Arguments:                                                           
//   cdim:   the dimensions of the cell grid (a 3-element vector)       
//   cdepth: returns the origin of the cell midpoint, expressed in terms of
//           cell length                                             
//   crd:    the coordinates                                            
//   dcinp:  the direct space parameters                                
//-----------------------------------------------------------------------------
void TakeCellGridDims(int* cdim, double* cdepth, coord *crd, dircon *dcinp)
{
  int i;
  double invMcut;
  dmat *invU;

  // Determine cell dimensions and the number of cells.  This uses the
  // Hessian normal form to compute the thickness of each cell (the   
  // distance between two cell faces), which works generally for      
  // non-orthorhombic as well as orthorhombic simulation boxes.       
  invU = &crd->invU;
  HessianNorms(invU, cdepth);
  invMcut = 1.0/dcinp->Mcut;
  for (i = 0; i < 3; i++) {
    cdim[i] = floor(cdepth[i]*invMcut);
    cdepth[i] /= cdim[i];
    cdepth[i] = 0.5/(cdepth[i]*invMcut);
  }
}

//-----------------------------------------------------------------------------
// IsCentralAtom: determine if this atom is in the center of the cell, at least
//                0.5 box lengths from any face of the extended cell.  Returns
//                0-7 (the sector number in which the atom is located) if true,
//                or -1 if false.                
//                                                                      
// Arguments:                                                           
//   corig:   the origin of the cell of interest                        
//   ccen:    the location of the origin of the cell central region     
//   atmloc:  the location of the atom within cell C                    
//   U:       the transformation matrix for converting to fractional
//            coordinates                                               
//   ng:      the number of grid cells in each dimension within the simulation
//-----------------------------------------------------------------------------
int IsCentralAtom(double* corig, double* ccen, double* atmloc, double* Umat,
		  double* ng, int isortho)
{
  int i;
  double dorig[3];

  if (isortho) {
    dorig[0] = (atmloc[0] - corig[0])*ng[0]*Umat[0];
    if (dorig[0] < ccen[0] || dorig[0] >= ccen[0]+1.0) {
      return -1;
    }
    dorig[1] = (atmloc[1] - corig[1])*ng[1]*Umat[4];
    if (dorig[1] < ccen[1] || dorig[1] >= ccen[1]+1.0) {
      return -1;
    }
    dorig[2] = (atmloc[2] - corig[2])*ng[2]*Umat[8];
    if (dorig[2] < ccen[2] || dorig[2] >= ccen[2]+1.0) {
      return -1;
    }
  }
  else {
    dorig[0] = atmloc[0] - corig[0];
    dorig[1] = atmloc[1] - corig[1];
    dorig[2] = atmloc[2] - corig[2];
    dorig[0] = dorig[0]*Umat[0] + dorig[1]*Umat[1] + dorig[2]*Umat[2];
    dorig[1] = dorig[0]*Umat[3] + dorig[1]*Umat[4] + dorig[2]*Umat[5];
    dorig[2] = dorig[0]*Umat[6] + dorig[1]*Umat[7] + dorig[2]*Umat[8];
    dorig[0] *= ng[0];
    if (dorig[0] < ccen[0] || dorig[0] >= ccen[0]+1.0) {
      return -1;
    }
    dorig[1] *= ng[1];
    if (dorig[1] < ccen[1] || dorig[1] >= ccen[1]+1.0) {
      return -1;
    }
    dorig[2] *= ng[2];
    if (dorig[2] < ccen[2] || dorig[2] >= ccen[2]+1.0) {
      return -1;
    }
  }
  dorig[0] = floor(dorig[0]);
  dorig[1] = floor(dorig[1]);
  dorig[2] = floor(dorig[2]);
  i = dorig[0] + 2.0*dorig[1] + 4.0*dorig[2] + 1.0e-5;

  return i;
}

//-----------------------------------------------------------------------------
// GetMaxDensityConfig: get the maximum density configuration of this system
//                      by determining massive atoms with the highest site to
//                      mass ratio and considering how many might fit in a
//                      single cell.                
//                                                                      
// Arguments:                                                           
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
static int GetMaxDensityConfig(prmtop *tp, cellgrid *CG, dmat *invU,
			       double dtarget)
{
  int i, j, nsp, maxatom;
  double rmass, maxsp, cvol;

  // If there are fewer than 100 atoms, we can just bail right out
  if (tp->natom <= 100) {
    return tp->natom;
  }

  // Determine the maximum number of residue sites per gram 
  maxsp = 0.0;
  for (i = 0; i < tp->nres; i++) {
    nsp = tp->ResLims[i+1] - tp->ResLims[i];
    rmass = 0.0;
    for (j = tp->ResLims[i]; j < tp->ResLims[i+1]; j++) {
      rmass += tp->Masses[j];
    }
    maxsp = MAX(maxsp, nsp/rmass);
  }

  // Determine the volume of a cell in mL
  cvol = invU->data[0]*invU->data[4]*invU->data[8] / 
    ((CG->ng[0]*CG->ng[1]*CG->ng[2])*1.0e24/AVOGADRO);

  // This is the maximum number of sites that will
  // fit into a cell at the maximum stated density
  maxatom = dtarget*cvol*maxsp;
  maxatom = MIN(maxatom, tp->natom+1);

  return maxatom;
}

//-----------------------------------------------------------------------------
// CreateCellGrid: this routine creates a cell grid, for storing atoms that
//                 will be involved in direct space interactions.  
//                                                                      
// Arguments:                                                           
//   crd:     the coordinates                                           
//   dcinp:   the direct space parameters                               
//   rcinp:   the reciprocal space parameters                           
//   CG:      the cell grid                                             
//   tj:      defined with MPI, the trajectory control information for access
//            to MPI mapping information                         
//   sysnum:  the number of the system this cell grid applies to        
//-----------------------------------------------------------------------------
cellgrid CreateCellGrid(coord *crd, dircon *dcinp, reccon *rcinp, prmtop *tp,
			trajcon *tj, int sysnum)
{
  int i, j, k, ix, iy, iz, nbx, nby, nbz;
  int cdim[3];
  double cdepth[3];
  cellgrid CG;

  // Compute the cell grid dimensions and allocate the cell grid
  TakeCellGridDims(cdim, cdepth, crd, dcinp);
  for (i = 0; i < 3; i++) {
    CG.ng[i] = cdim[i];
    CG.dbng[i] = CG.ng[i];
    CG.celldim[i] = crd->gdim[i]/CG.ng[i];
    CG.celldim[i+4] = cdepth[i];
  }
  CG.celldim[3] = dcinp->Mcut;
  CG.ncell = CG.ng[0]*CG.ng[1]*CG.ng[2];
  CG.data = (cell*)malloc(CG.ng[0]*CG.ng[1]*CG.ng[2]*sizeof(cell));
  CG.sysID = sysnum;
  
  // Determine a suitable number of atoms that each cell
  // can be expected to hold.  Note that this number    
  // CANNOT be increased later in the simulation.       
  CG.maxatom = GetMaxDensityConfig(tp, &CG, &crd->invU, dcinp->MaxDens);

  // Allocate GPS buffers to rapidly cross-reference atoms
  // in the master topology and in the cell data arrays.  
  // Link each cell to one of the GPS buffers; the buffers
  // and cells that point to them are arranged such that  
  // no one GPS buffer (each row of the matrix AtmGPS) can
  // have more than one copy of any atom.    
  nbx = (CG.ng[0] % 2 == 0) ? 2 : 3;
  nby = (CG.ng[1] % 2 == 0) ? 2 : 3;
  nbz = (CG.ng[2] % 2 == 0) ? 2 : 3;
  CG.AtmGPS = CreateImat(nbx*nby*nbz, tp->natom);

  // Allocate each cell and set its GPS pointer
  CG.map = (cell***)malloc(CG.ng[0]*sizeof(cell**));
  for (i = 0; i < CG.ng[0]; i++) {
    ix = (nbx == 3 && i == CG.ng[0]-1) ? 2 : i % 2;
    CG.map[i] = (cell**)malloc(CG.ng[1]*sizeof(cell*));
    for (j = 0; j < CG.ng[1]; j++) {
      iy = (nby == 3 && j == CG.ng[1]-1) ? 2 : j % 2;
      CG.map[i][j] = &CG.data[(i*CG.ng[1] + j)*CG.ng[2]];
      for (k = 0; k < CG.ng[2]; k++) {
	iz = (nbz == 3 && k == CG.ng[2]-1) ? 2 : k % 2;
	CG.map[i][j][k] = CreateCell(CG.maxatom, rcinp->ordr);
	CG.map[i][j][k].GPSptr = CG.AtmGPS.map[(ix*nby + iy)*nbz + iz];
      }
    }
  }

  // Compute cell origins
  ComputeCellOrigins(&CG, crd);

  // Compute the normal vectors to each cell
  // face for movement in x, y, and z       
  ComputeCellNormals(crd);

  // Processes that will tend this cell grid.
  // The communicator is created here, and   
  // messaging plans are created next in     
  // the LinkCellGrid() function below.      
  CellGridEquipComm(&CG, tj);

  // Check to see that the cell grid is of acceptable dimensions
  for (i = 0; i < 3; i++) {
    if (cdim[i] < 2) {
      printf("CreateCellGrid >> Error.  Cell grid must be at least two "
	     "cells in all\nCreateCellGrid >> dimensions such that the "
	     "simulation cell is at least\nCreateCellGrid >> twice the "
	     "nonbonded cutoff in all dimensions.  This likely\n"
	     "CreateCellGrid >> occurred because the simulation system is "
	     "too small.\n");
      exit(1);
    }
  }

  return CG;
}

//-----------------------------------------------------------------------------
// LinkCellGrid: after the cell grid communicator has been created above in
//               the CreateCellGrid() function, some additional work may need
//               to be done in order to determine properties of the reciprocal
//               space mesh.                             
//                                                                      
// Arguments:                                                           
//   CG:       the cell grid to link                                    
//   crd:      system coordinates                                       
//   rcinp:    reciprocal space control information                     
//-----------------------------------------------------------------------------
#ifdef MPI
void LinkCellGrid(cellgrid *CG, coord *crd, reccon *rcinp)
#else
void LinkCellGrid(cellgrid *CG, reccon *rcinp)
#endif
{
  InitLoadBalance(CG, rcinp);
#ifdef MPI
  MapProcessMeshFootprint(CG, rcinp, crd);
#endif
  MapProcessCellSharing(CG);
  AllocatePooledBuffers(CG);
  PlanCoordRedux(CG);
}

//-----------------------------------------------------------------------------
// DestroyCell: this routine destroys a cell.                           
//                                                                      
// Arguments:                                                           
//   C:  the cell                                                       
//-----------------------------------------------------------------------------
void DestroyCell(cell *C)
{
  int i;

  free(C->nr);
  free(C->nsr);
  free(C->xcof);
  free(C->ycof);
  free(C->zcof);
  free(C->map);
  free(C->data);
  free(C->import);
  free(C->pexport);
  free(C->Vimport);
  free(C->Vexport);
  free(C->Ximport);
  free(C->Xexport);
  free(C->atmscr);
  free(C->qr2buff);
  free(C->qIDbuff);
  free(C->ljr2buff);
  free(C->ljIDbuff);
  DestroyImat(&C->ordr);
  DestroyImat(&C->supordr);
  if (C->nFscr > 0) {
    for (i = 0; i < C->nFscr; i++) {
      DestroyFbook(&C->Fscr[i]);
    }
    free(C->Fscr);
  }
}

//-----------------------------------------------------------------------------
// DestroyCellGrid: this routine destroys a cell grid.                  
//                                                                      
// Arguments:                                                           
//   CG:  the cell grid                                                 
//-----------------------------------------------------------------------------
void DestroyCellGrid(cellgrid *CG)
{
  int i;

  for (i = 0; i < CG->ng[0]*CG->ng[1]*CG->ng[2]; i++) {
    DestroyCell(&CG->data[i]);
  }
  for (i = 0; i < CG->ng[0]; i++) {
    free(CG->map[i]);
  }
  free(CG->map);
  free(CG->data);

  // Free pooled import / export buffers
  for (i = 0; i < CG->nsend; i++) {
    free(CG->pexport[i]);
    free(CG->Vexport[i]);
    free(CG->Xexport[i]);
  }
  for (i = 0; i < CG->nrecv; i++) {
    free(CG->import[i]);
    free(CG->Vimport[i]);
    free(CG->Ximport[i]);
  }
  free(CG->maximp);
  free(CG->maxexp);
  free(CG->pexport);
  free(CG->import);
  free(CG->Vexport);
  free(CG->Vimport);
  free(CG->Xexport);
  free(CG->Ximport);
  free(CG->nexp);

#ifdef MPI
  // Free coordinate reduction buffers
  int nthr, rank, planstr, planlim;
  MPI_Comm_rank(CG->dspcomm, &rank);
  MPI_Comm_size(CG->dspcomm, &nthr);
  if (nthr > 1) {
    if (rank == 0) {
      for (i = 0; i < nthr-1; i++) {
	free(CG->CrdPool[i]);
      }
    }
    else {
      free(CG->CrdPool[0]);
    }
    free(CG->CrdPool);
    free(CG->CrdPoolSize);
  }
  planstr = (CG->tid == 0) ? 1 : 0;
  planlim = (CG->tid == 0) ? CG->nthreads : 1;
  for (i = planstr; i < planlim; i++) {
    free(CG->MeshCommPlan[i].pagl);
    free(CG->MeshCommPlan[i].coll);
    free(CG->MeshCommPlan[i].pcl);
  }
  free(CG->MeshCommPlan);
  MPI_Comm_free(&CG->dspcomm);
#endif

  // Free direct space communication plans
  for (i = 0; i < 3; i++) {
    DestroyAshr(&CG->DirCommPlan.mvshr[i]);
    DestroyAshr(&CG->DirCommPlan.frcmg[i]);
  }

  // Free cell domain
  free(CG->MyCellDomain);

  // Free atom GPS
  DestroyImat(&CG->AtmGPS);
}

//-----------------------------------------------------------------------------
// ComputeCellOrigins: compute the origin of each cell, to act as a reference
//                     point during cell-cell transfers and other operations.
//                                                                      
// Arguments:                                                           
//   CG:     the cell grid                                              
//   crd:    the coordinates                                            
//-----------------------------------------------------------------------------
void ComputeCellOrigins(cellgrid *CG, coord *crd)
{
  int h, i, j, k, nx, ny, nz;
  double mmx, mmy, mmz, mmxp1, mmyp1, mmzp1, invnx, invny, invnz;
  double *utmp;
  cell** ctm2p;
  cell* ctmp;

  utmp = crd->invU.data;
  nx = CG->ng[0];
  ny = CG->ng[1];
  nz = CG->ng[2];

  // The orthorhombic case
  if (crd->isortho == 1) {
    h = 0;
    for (i = 0; i < nx; i++) {
      ctm2p = CG->map[i];
      mmx = i*CG->celldim[0];
      mmxp1 = mmx+CG->celldim[0];
      for (j = 0; j < ny; j++) {
	ctmp = ctm2p[j];
	mmy = j*CG->celldim[1];
	mmyp1 = mmy+CG->celldim[1];
	for (k = 0; k < nz; k++) {
	  mmz = k*CG->celldim[2];
	  ctmp[k].orig[0] = mmx;
	  ctmp[k].orig[1] = mmy;
	  ctmp[k].orig[2] = mmz;
	  ctmp[k].midp[0] = mmxp1;
	  ctmp[k].midp[1] = mmyp1;
	  ctmp[k].midp[2] = mmz+CG->celldim[2];
	  ctmp[k].gbin[0] = i;
	  ctmp[k].gbin[1] = j;
	  ctmp[k].gbin[2] = k;
	  ctmp[k].gbin[3] = h;
	  h++;
	}
      }
    }
    return;
  }

  // The non-orthorhombic case
  invnx = 1.0/nx;
  invny = 1.0/ny;
  invnz = 1.0/nz;
  h = 0;
  for (i = 0; i < nx; i++) {
    ctm2p = CG->map[i];
    mmx = i*invnx;
    mmxp1 = mmx+invnx;
    for (j = 0; j < ny; j++) {
      ctmp = ctm2p[j];
      mmy = j*invny;
      mmyp1 = mmy+invny;
      for (k = 0; k < nz; k++) {
	mmz = k*invnz;
	mmzp1 = mmz+invnz;
	ctmp[k].orig[0] = utmp[0]*mmx + utmp[1]*mmy + utmp[2]*mmz;
	ctmp[k].orig[1] = utmp[3]*mmx + utmp[4]*mmy + utmp[5]*mmz;
	ctmp[k].orig[2] = utmp[6]*mmx + utmp[7]*mmy + utmp[8]*mmz;
        ctmp[k].midp[0] = utmp[0]*mmxp1 + utmp[1]*mmyp1 + utmp[2]*mmzp1;
        ctmp[k].midp[1] = utmp[3]*mmxp1 + utmp[4]*mmyp1 + utmp[5]*mmzp1;
        ctmp[k].midp[2] = utmp[6]*mmxp1 + utmp[7]*mmyp1 + utmp[8]*mmzp1;
	ctmp[k].gbin[0] = i;
	ctmp[k].gbin[1] = j;
	ctmp[k].gbin[2] = k;
	ctmp[k].gbin[3] = h;
	h++;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SetPoint: set a point P to row n of the three-column matrix A, or to the
//           average of rows n and m (if m is non-negative).        
//
// Arguments:
//   P:        the point
//   A:        coordinate matrix
//   n:        row of the matrix A
//   m:        second, optional row of the matrix A
//-----------------------------------------------------------------------------
static void SetPoint(double* P, dmat *A, int n, int m)
{
  P[0] = A->map[n][0];
  P[1] = A->map[n][1];
  P[2] = A->map[n][2];
  if (m >= 0) {
    P[0] = 0.5*(P[0] + A->map[m][0]);
    P[1] = 0.5*(P[1] + A->map[m][1]);
    P[2] = 0.5*(P[2] + A->map[m][2]);
  }
}

//-----------------------------------------------------------------------------
// ComputeCellNormals: compute vectors normal to planes dividing cells along
//                     the axis of each cell:cell interaction.    
//                                                                      
// Arguments:                                                           
//   crd:    the coordinates                                            
//-----------------------------------------------------------------------------
void ComputeCellNormals(coord *crd)
{
  int i;
  double P1[3], P2[3], P3[3];
  double *dtmp, *dtm2p;
  dmat nvec;

  // Compute six unit vectors, representing
  // each of the box vectors.              
  nvec = CreateDmat(6, 3, 0);
  for (i = 0; i < 3; i++) {
    nvec.map[i][0] = crd->invU.map[0][i];     // Unit cell (+)a-vector
    nvec.map[i][1] = crd->invU.map[1][i];     // Unit cell (+)b-vector
    nvec.map[i][2] = crd->invU.map[2][i];     // Unit cell (+)g-vector
    nvec.map[3+i][0] = -1.0*nvec.map[i][0];   // Unit cell (-)a-vector
    nvec.map[3+i][1] = -1.0*nvec.map[i][1];   // Unit cell (-)b-vector
    nvec.map[3+i][2] = -1.0*nvec.map[i][2];   // Unit cell (-)g-vector
    UnitVector3(nvec.map[i]);
    UnitVector3(nvec.map[3+i]);
  }

  // For all plane computations from here on out, P1 is the origin
  P1[0] = 0.0;
  P1[1] = 0.0;
  P1[2] = 0.0;

  // Normals to each cell face are stored in the ROWS of fcnorm. 
  // If box vectors are labelled a, b, and c, then the rows of   
  // fcnorm are normal to the bc, ac, and ab faces, respectively.
  UnitNormal2Plane(crd->fcnorm.map[0], P1, nvec.map[2], nvec.map[1]);
  UnitNormal2Plane(crd->fcnorm.map[1], P1, nvec.map[0], nvec.map[2]);
  UnitNormal2Plane(crd->fcnorm.map[2], P1, nvec.map[1], nvec.map[0]);

  // Normals to each cell edgewise divider are stored in the rows
  // of ednorm, which is not a square matrix.  The first through 
  // sixth rows pertain to interactions between sectors 0->3,    
  // 0->5, 0->6, 1->2, 1->4, and 2->4, respectively.             
  SetPoint(P2, &nvec, 5, -1);
  SetPoint(P3, &nvec, 1, 3);
  UnitNormal2Plane(crd->fcnorm.map[3], P1, P3, P2);
  SetPoint(P2, &nvec, 4, -1);
  SetPoint(P3, &nvec, 2, 3);
  UnitNormal2Plane(crd->fcnorm.map[4], P1, P2, P3);
  SetPoint(P2, &nvec, 3, -1);
  SetPoint(P3, &nvec, 4, 2);
  UnitNormal2Plane(crd->fcnorm.map[5], P1, P3, P2);
  SetPoint(P2, &nvec, 5, -1);
  SetPoint(P3, &nvec, 0, 1);
  UnitNormal2Plane(crd->fcnorm.map[6], P1, P2, P3);
  SetPoint(P2, &nvec, 4, -1);
  SetPoint(P3, &nvec, 2, 0);
  UnitNormal2Plane(crd->fcnorm.map[7], P1, P3, P2);
  SetPoint(P2, &nvec, 3, -1);
  SetPoint(P3, &nvec, 1, 2);
  UnitNormal2Plane(crd->fcnorm.map[8], P1, P2, P3);

  // Normals to each cell vertex divider are stored in the rows 
  // of vtnorm, which is a square matrix but is probably not    
  // convenient for any rotations.  The first through third rows
  // pertain to interactions between sectors 0->7, 1->6, 2->5,  
  // and 4->3, respectively.                                    
  SetPoint(P2, &nvec, 4, 0);
  SetPoint(P3, &nvec, 2, 3);
  UnitNormal2Plane(crd->fcnorm.map[9], P1, P2, P3);
  SetPoint(P2, &nvec, 4, 3);
  SetPoint(P3, &nvec, 2, 0);
  UnitNormal2Plane(crd->fcnorm.map[10], P1, P3, P2);
  SetPoint(P2, &nvec, 4, 3);
  SetPoint(P3, &nvec, 5, 4);
  UnitNormal2Plane(crd->fcnorm.map[11], P1, P3, P2);
  SetPoint(P2, &nvec, 4, 0);
  SetPoint(P3, &nvec, 3, 5);
  UnitNormal2Plane(crd->fcnorm.map[12], P1, P3, P2);

  // Negative unit normals
  dtmp = crd->fcnorm.data;
  dtm2p = &crd->fcnorm.data[39];
  for (i = 0; i < 39; i++) {
    dtm2p[i] = -1.0*dtmp[i];
  }

  // Free allocated memory
  DestroyDmat(&nvec);
}

//-----------------------------------------------------------------------------
// UpdateCellGPS: this function will update the coord struct's AtmGPS array
//                with the sectors and indices where each of the cell's atoms
//                can be found.                            
//                                                                      
// Arguments:                                                           
//   CG:        the cell grid, containing its own GPS guide             
//-----------------------------------------------------------------------------
void UpdateCellGPS(cellgrid *CG)
{
  int i, j, k, maxatom, jmax;
  cell *C;
  atomc *catm;

  // Loop over all cells and set the number of instances
  // of each atom encountered in each one to zero.      
  maxatom = CG->maxatom;
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < 8; j++) {
      catm = C->map[j];
      jmax = maxatom * j;
      for (k = 0; k < C->nr[j]; k++) {
        C->GPSptr[catm[k].id] = jmax + k;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SortAtomID: function called by quicksort for comparing the ID numbers of two
//             atomc structs.                                    
//                                                                      
// Arguments:                                                           
//   atm[A,B]:   the atomc structs                                      
//-----------------------------------------------------------------------------
int SortAtomID(const void *atmA, const void *atmB)
{
  int idA = ((atomc*)atmA)[0].id;
  int idB = ((atomc*)atmB)[0].id;

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
// AtomsToCells: preliminary import of atoms into cells.                
//                                                                      
// Arguments:                                                           
//   crd:    the coordinates                                            
//   CG:     the cell grid                                              
//   tp:     the topology                                               
//-----------------------------------------------------------------------------
void AtomsToCells(coord *crd, cellgrid *CG, prmtop *tp)
{
  int h, i, j;
  int nc[3];
  double celldim[3], gdim[3], invgdim[3], invcelldim[3], tmppos[3], bcrd[3];
  cell *C;

  // Pre-compute cell dimensions if this is a rectangular system
  if (crd->isortho == 1) {
    for (i = 0; i < 3; i++) {
      gdim[i] = crd->gdim[i];
      celldim[i] = CG->celldim[i];
      invgdim[i] = 1.0/gdim[i];
      invcelldim[i] = 1.0/celldim[i];
    }
  }

  // Set the number of atoms in all cells to zero
  for (i = 0; i < CG->ng[0]*CG->ng[1]*CG->ng[2]; i++) {
    for (j = 0; j < 8; j++) {
      CG->data[i].nr[j] = 0;
    }
  }

  // Loop over all atoms and assign each to a cell
  for (h = 0; h < crd->natom; h++) {

    // Skip virtual sites / extra points
    if (tp->Masses[h] < 1.0e-8) {
      continue;
    }
    if (crd->isortho == 0) {

      // In a triclinic system, we must convert coordinates
      // into box space before assigning them to cells     
      DMatVecMult(&crd->U, &crd->loc[3*h], bcrd);
      for (i = 0; i < 3; i++) {
	bcrd[i] -= floor(bcrd[i]);
	nc[i] = bcrd[i] * CG->ng[i];
      }
      DMatVecMult(&crd->invU, bcrd, tmppos);
    }
    else {

      // Is the system is rectangular, then we can
      // use the pre-computed cell dimensions     
      for (i = 0; i < 3; i++) {
	tmppos[i] = crd->loc[3*h+i];
      }
      for (i = 0; i < 3; i++) {
	tmppos[i] = tmppos[i]*invgdim[i];
	tmppos[i] = gdim[i]*(tmppos[i] - floor(tmppos[i]));
	nc[i] = tmppos[i]*invcelldim[i];
	if (nc[i] >= CG->ng[i] &&
	    tmppos[i]*invcelldim[i] - 1.0e-8 < CG->ng[i]) {
	  nc[i] = tmppos[i]*invcelldim[i] - 1.0e-8;
	}
      }
    }
    C = &CG->map[nc[0]][nc[1]][nc[2]];
    j = C->nr[0];
    C->map[0][j].id = h;
    C->map[0][j].lj = tp->LJIdx[h];
    C->map[0][j].q = tp->Charges[h];
    for (i = 0; i < 3; i++) {
      C->map[0][j].loc[i] = tmppos[i];
      C->map[0][j].frc[i] = 0.0;
    }
    j++;
    C->nr[0] = j;
    if (j == C->maxatom && tp->natom / C->maxatom > 8) {
      printf("AtomsToCellsOrtho >> Warning.  Atom count in cell [ %2d %2d "
	     "%2d ] has met\nAtomsToCellsOrtho >> the maximum of %d.  It is "
	     "likely that the run will need to be restarted\n"
	     "AtomsToCellsOrtho >> with a higher value of rho.\n",
	     nc[0], nc[1], nc[2], C->maxatom);
    }
  }
}

//-----------------------------------------------------------------------------
// MixCells: mix the forces in two cells.  This function is extracted from  
//           MixCellGrids because it may be convenient to just bail out of it.
//
// Arguments:
//   cA:     the first cell
//   cB:     the scond cell
//   prc:    stores correlations among topology structures
//   mxA:    mixing factor for the first cell cA
//   mxB:    mixing factor for the second cell cB
//-----------------------------------------------------------------------------
static void MixCells(cell *cA, cell *cB, prmcorr *prc, double mxA, double mxB)
{
  int i, j, acon, bcon, aseek;
  double df;
  atomc *atmA, *atmB;

  // Simplest case: all atoms in one topology 
  // correspond directly to atoms in the other
  if (prc->relate == 0) {
    for (i = 0; i < cA->nr[0]; i++) {
      for (j = 0; j < 3; j++) {
	df = mxA*cA->data[i].frc[j] + mxB*cB->data[i].frc[j];
	cA->data[i].frc[j] = df;
	cB->data[i].frc[j] = df;
      }
    }
  }

  // Not-so-bad case: atoms in one topology correspond
  // to atoms in the other, but there may be skips    
  else if (prc->relate == 1) {
    acon = 0;
    bcon = 0;
    while (acon < cA->nr[0] && bcon < cB->nr[0]) {

      // Set pointers to atoms in each cell
      atmA = &cA->data[acon];
      atmB = &cB->data[bcon];

      // If the atom in cell A does not correspond
      // to anything in topology B and therefore  
      // does not correspond to anything in cell B
      if (prc->matchA[atmA->id] == -1) {
	acon++;
	continue;
      }

      // Same as above, but for atom B
      if (prc->matchB[atmB->id] == -1) {
	bcon++;
	continue;
      }

      // Atoms A and B had better be partners
      for (j = 0; j < 3; j++) {
	df = mxA*atmA->frc[j] + mxB*atmB->frc[j];
	atmA->frc[j] = df;
	atmB->frc[j] = df;
      }
      acon++;
      bcon++;
    }
  }

  // Really bad case: atoms in one topology may or may not  
  // correspond to atoms in the other, and there is no order
  else if (prc->relate == 2) {
    for (i = 0; i < cA->nr[0]; i++) {
      aseek = prc->matchA[cA->data[i].id];
      if (aseek == -1) {
	continue;
      }
      atmA = &cA->data[i];
      atmB = &cB->data[cB->GPSptr[aseek]];
      for (j = 0; j < 3; j++) {
	df = mxA*atmA->frc[j] + mxB*atmB->frc[j];
	atmA->frc[j] = df;
	atmB->frc[j] = df;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// MixCellGrids: in cases where two cell grids have corresponding atoms, in
//               order to synchronize the trajectories forces on the 
//               corresponding atoms must be mixed.                     
//                                                                      
// CG[A,B]:   the cell grids                                            
// tj:        contains mixing constant information                      
//-----------------------------------------------------------------------------
void MixCellGrids(cellgrid *CGA, cellgrid *CGB, trajcon *tj)
{
  int i;
  double mxA, mxB;
  prmcorr *prc;

  prc = &tj->prc;

  // mxB is the mixing factor for the second   
  // set for forces due to the second topology.
  // mxA is the mixing factor for the first.   
  mxA = 1.0;
  for (i = 0; i < tj->mxorder; i++) {
    mxA *= (1.0 - tj->lambda);
  }
  mxB = 1.0 - mxA;

  // Loop over primary sectors of all cells 
  // in CG and mix forces with atoms in CGN.
  // While the two cells grids may not have 
  // all the same atoms, they will be of    
  // similar sizes and most atoms will have 
  // counterparts in each grid.             
  for (i = 0; i < CGA->MyCellCount; i++) {
    int j = CGA->MyCellDomain[i];
    MixCells(&CGA->data[j], &CGB->data[j], prc, mxA, mxB);
  }
}

//-----------------------------------------------------------------------------
// ZeroCellForces: zero forces on atoms in the primary sectors of all cells.
//                                                                      
// Arguments:                                                           
//   CG:     the cell grid                                              
//-----------------------------------------------------------------------------
void ZeroCellForces(cellgrid *CG)
{
  int i, j;
  cell *C;

  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < C->nr[0]; j++) {
      C->data[j].frc[0] = 0.0;
      C->data[j].frc[1] = 0.0;
      C->data[j].frc[2] = 0.0;
    }
  }
}

//-----------------------------------------------------------------------------
// PackageCell: this function packages atoms in a cell sector for export to
//              another cell.  This routine is called by ShareCellContents,
//              and as such it can only package one type of information about
//              atoms such as location, velocity or force.  Works for
//              orthorhombic unit cells.  
//                                                                      
// Arguments:                                                           
//   C:         current cell, owned by this process                     
//   crd:       the coordinates, as understood by this process          
//   imove:     the direction of the move from C -> D                   
//   needreim:  flag to indicate whether reimaging is needed            
//   reimv:     the re-imaging vector                                   
//   Mcut:      the maximum direct space cutoff                         
//   i[C,D]:    regions within cells C and D that are being transferred 
//   ct:        flag to indicate what is being transferred              
//   tp:        the topology                                            
//-----------------------------------------------------------------------------
static void PackageCell(cell *C, coord *crd, int imove, int needreim,
			double* reimv, double Mcut, int iC, int iD,
			int ct, prmtop *tp)
{
  int i, nexp;
  double clim;
  double *vtmp, *fcnorm;
  double* crange;
  atomc *catm;

  // Counter for the current number of atoms to export
  nexp = C->nexp;
  catm = C->map[iC];

  // The relevant limits
  crange = (double*)malloc(C->nr[iC]*sizeof(double));
  fcnorm = crd->fcnorm.map[imove];
  if (crd->isortho == 0) {
    for (i = 0; i < C->nr[iC]; i++) {
      crange[i] = fcnorm[0]*(catm[i].loc[0] - C->orig[0]) +
        fcnorm[1]*(catm[i].loc[1] - C->orig[1]) +
        fcnorm[2]*(catm[i].loc[2] - C->orig[2]);
    }
    clim = Mcut;
  }
  else {
    for (i = 0; i < C->nr[iC]; i++) {
      crange[i] = catm[i].loc[imove];
    }
    clim = C->orig[imove] + Mcut;
  }

  // Transfer locations
  if (ct == 0) {
    if (needreim) {
      for (i = 0; i < C->nr[iC]; i++) {
	if (crange[i] < clim && tp->Masses[catm[i].id] > 1.0e-8) {
	  C->pexport[nexp].id = catm[i].id;
	  C->pexport[nexp].dreg = iD;
	  C->pexport[nexp].loc[0] = catm[i].loc[0] + reimv[0];
	  C->pexport[nexp].loc[1] = catm[i].loc[1] + reimv[1];
	  C->pexport[nexp].loc[2] = catm[i].loc[2] + reimv[2];
	  nexp++;
	}
      }
    }
    else {
      for (i = 0; i < C->nr[iC]; i++) {
	if (crange[i] < clim && tp->Masses[catm[i].id] > 1.0e-8) {
	  C->pexport[nexp].id = catm[i].id;
	  C->pexport[nexp].dreg = iD;
	  C->pexport[nexp].loc[0] = catm[i].loc[0];
	  C->pexport[nexp].loc[1] = catm[i].loc[1];
	  C->pexport[nexp].loc[2] = catm[i].loc[2];
	  nexp++;
	}
      }
    }
  }

  // Transfer velocities
  else if (ct == 1) {
    for (i = 0; i < C->nr[iC]; i++) {
      if (crange[i] < clim) {
	C->pexport[nexp].id = catm[i].id;
	C->pexport[nexp].dreg = iD;
	vtmp = &crd->vel[3*catm[i].id];
	C->pexport[nexp].loc[0] = vtmp[0];
	C->pexport[nexp].loc[1] = vtmp[1];
	C->pexport[nexp].loc[2] = vtmp[2];
	nexp++;
      }
    }
  }

  // Transfer current and previous positions
  else if (ct == 2) {
    int atmid3;
    if (needreim) {
      for (i = 0; i < C->nr[iC]; i++) {
	if (crange[i] < clim && tp->Masses[catm[i].id] > 1.0e-8) {
	  atmid3 = 3*catm[i].id;
	  C->Vexport[nexp].id = catm[i].id;
	  C->Vexport[nexp].dreg = iD;
	  C->Vexport[nexp].loc[0] = catm[i].loc[0] + reimv[0];
	  C->Vexport[nexp].loc[1] = catm[i].loc[1] + reimv[1];
	  C->Vexport[nexp].loc[2] = catm[i].loc[2] + reimv[2];
	  C->Vexport[nexp].vel[0] = crd->prvloc[atmid3];
	  C->Vexport[nexp].vel[1] = crd->prvloc[atmid3+1];
	  C->Vexport[nexp].vel[2] = crd->prvloc[atmid3+2];
	  nexp++;	  
	}
      }
    }
    else {
      for (i = 0; i < C->nr[iC]; i++) {
	if (crange[i] < clim && tp->Masses[catm[i].id] > 1.0e-8) {
	  atmid3 = 3*catm[i].id;
	  C->Vexport[nexp].id = catm[i].id;
	  C->Vexport[nexp].dreg = iD;
	  C->Vexport[nexp].loc[0] = catm[i].loc[0];
	  C->Vexport[nexp].loc[1] = catm[i].loc[1];
	  C->Vexport[nexp].loc[2] = catm[i].loc[2];
	  C->Vexport[nexp].vel[0] = crd->prvloc[atmid3];
	  C->Vexport[nexp].vel[1] = crd->prvloc[atmid3+1];
	  C->Vexport[nexp].vel[2] = crd->prvloc[atmid3+2];
	  nexp++;
	}
      }
    }
  }

  // Update the export count
  C->nexp = nexp;

  // Free allocated memory
  free(crange);
}

//-----------------------------------------------------------------------------
// UnpackCellContents: this function, much simpler than the related     
//                     PackageCell(Ortho,Nonortho) functions above, will unpack
//                     imported atoms in the destination cell.  Because this
//                     function is only called with MPI compilation, it assumes
//                     that messages passed between cells contain the number of
//                     atoms passed as the atom ID field of the first atomb
//                     struct, and then the list of real atoms thereafter.
//
// Arguments:                                                           
//   D:      the destination cell                                       
//   ct:     the type of information to transfer (0 = positions,        
//           1 = velocties)                                             
//   tp:     the topology                                               
//   crd:    the coordinates                                            
//-----------------------------------------------------------------------------
static void UnpackCellContents(cell *D, int ct, prmtop *tp, coord *crd)
{
  int i, dreg;
  double *vtmp;
  atomc *datm;
  atomb *dimp;
  atombv *dVimp;

  // Get the size of the import from the first "atom"
  D->nimp = (ct < 2) ? D->import[0].id : D->Vimport[0].id;

  // Loop over all imported atoms, with behavior 
  // dependent on the type of data being imported
  if (ct == 0) {
    for (i = 1; i < D->nimp; i++) {
      dimp = &D->import[i];
      dreg = dimp->dreg;
      datm = &D->map[dreg][D->nr[dreg]];
      datm->id = dimp->id;
      datm->q = tp->Charges[dimp->id];
      datm->lj = tp->LJIdx[dimp->id];
      datm->loc[0] = dimp->loc[0];
      datm->loc[1] = dimp->loc[1];
      datm->loc[2] = dimp->loc[2];
      datm->frc[0] = 0.0;
      datm->frc[1] = 0.0;
      datm->frc[2] = 0.0;
      D->nr[dreg] += 1;
    }
  }
  else if (ct == 1) {
    for (i = 1; i < D->nimp; i++) {
      dimp = &D->import[i];

      // Note that this merely sets the velocity of the atom 
      // in the associated coord struct on the receiving cell
      vtmp = &crd->vel[3*dimp->id];
      vtmp[0] = dimp->loc[0];
      vtmp[1] = dimp->loc[1];
      vtmp[2] = dimp->loc[2];
    }
  }
  else if (ct == 2) {
    for (i = 1; i < D->nimp; i++) {
      dVimp = &D->Vimport[i];
      dreg = dVimp->dreg;
      datm = &D->map[dreg][D->nr[dreg]];
      datm->id = dVimp->id;
      datm->q = tp->Charges[dVimp->id];
      datm->lj = tp->LJIdx[dVimp->id];
      vtmp = &crd->prvloc[3*dVimp->id];
      datm->loc[0] = dVimp->loc[0];
      datm->loc[1] = dVimp->loc[1];
      datm->loc[2] = dVimp->loc[2];
      vtmp[0] = dVimp->vel[0];
      vtmp[1] = dVimp->vel[1];
      vtmp[2] = dVimp->vel[2];
      datm->frc[0] = 0.0;
      datm->frc[1] = 0.0;
      datm->frc[2] = 0.0;
      D->nr[dreg] += 1;
    }
  }
}

//-----------------------------------------------------------------------------
// HandoffCell: this function performs what PackageCell and UnpackCellContents
//              would in the case that the same process owns both the sending
//              and receiving cells.  This function is only used in the case of
//              position information transfer, as velocity information will
//              already be consistent between all cells sharing the same
//              process.
//
// Arguments:                                                           
//   C,D:       cells C (current) and D (destination)                   
//   i[C,D]:    regions within cells C and D that are being transferred 
//   imove:     the direction of the move from C -> D                   
//   needreim:  flag to indicate whether reimaging is needed            
//   reimv:     the re-imaging vector                                   
//   Mcut:      the maximum direct space cutoff                         
//   tp:        the topology                                            
//   isortho:   flag to tell whether the simulation box is orthorhombic 
//   fcnorm:     vector normal to the face between the two cells        
//-----------------------------------------------------------------------------
static void HandoffCell(cell *C, cell *D, int iC, int iD, int imove,
			int needreim, double* reimv, double Mcut,
			prmtop *tp, int isortho, double* fcnorm)
{
  int i, nexp;
  double clim;
  double* crange;
  atomc *catm, *datm;

  // Pointers to sending and receiving cell atom arrays
  catm = C->map[iC];
  datm = D->map[iD];

  // The relevant limits
  crange = (double*)malloc(C->nr[iC]*sizeof(double));
  if (isortho == 0) {
    for (i = 0; i < C->nr[iC]; i++) {
      crange[i] = fcnorm[0]*(catm[i].loc[0] - C->orig[0]) +
	fcnorm[1]*(catm[i].loc[1] - C->orig[1]) +
	fcnorm[2]*(catm[i].loc[2] - C->orig[2]);
    }
    clim = Mcut;
  }
  else {
    for (i = 0; i < C->nr[iC]; i++) {
      crange[i] = catm[i].loc[imove];
    }
    clim = C->orig[imove] + Mcut;
  }

  // Count of transferred atoms
  nexp = 0;
  if (needreim) {
    if (tp->nxtrapt > 0) {
      for (i = 0; i < C->nr[iC]; i++) {
	if (crange[i] < clim && tp->Masses[catm[i].id] > 1.0e-8) {
	  datm[nexp] = catm[i];
	  datm[nexp].loc[0] += reimv[0];
	  datm[nexp].loc[1] += reimv[1];
	  datm[nexp].loc[2] += reimv[2];
	  datm[nexp].frc[0] = 0.0;
	  datm[nexp].frc[1] = 0.0;
	  datm[nexp].frc[2] = 0.0;
	  nexp++;
	}
      }
    }
    else {
      for (i = 0; i < C->nr[iC]; i++) {
        if (crange[i] < clim) {
          datm[nexp] = catm[i];
	  datm[nexp].loc[0] += reimv[0];
	  datm[nexp].loc[1] += reimv[1];
	  datm[nexp].loc[2] += reimv[2];
          datm[nexp].frc[0] = 0.0;
          datm[nexp].frc[1] = 0.0;
          datm[nexp].frc[2] = 0.0;
          nexp++;
        }
      }
    }
  }
  else {
    if (tp->nxtrapt > 0) {
      for (i = 0; i < C->nr[iC]; i++) {
	if (crange[i] < clim && tp->Masses[catm[i].id] > 1.0e-8) {
          datm[nexp] = catm[i];
	  datm[nexp].frc[0] = 0.0;
	  datm[nexp].frc[1] = 0.0;
	  datm[nexp].frc[2] = 0.0;
	  nexp++;
	}
      }
    }
    else {
      for (i = 0; i < C->nr[iC]; i++) {
        if (crange[i] < clim) {
          datm[nexp] = catm[i];
          datm[nexp].frc[0] = 0.0;
          datm[nexp].frc[1] = 0.0;
          datm[nexp].frc[2] = 0.0;
          nexp++;
        }
      }
    }
  }

  // The population of the receiving cell sector
  // is the number of transferred atoms         
  D->nr[iD] = nexp;

  // Free allocated memory
  free(crange);
}

//-----------------------------------------------------------------------------
// ShareCellContents: this function shares the contents of up to four sectors  
//                    in a cell with the relevant neighboring cell.  This
//                    function makes what would otherwise be a set of seven
//                    messages for importing atom coordinates into just three.
//                                                                      
// Arguments:                                                           
//   CG:           the cell grid                                        
//   crd:          the coordinates                                      
//   cloc:         the location of the sending cell                     
//   imove:        the dimension in which transfer occurs (0 = -X, 1 = -Y,
//                 2 = -Z)                                      
//   Mcut:         the maximum direct space cutoff (electrostatic and van-der  
//                 Waals forces)                                
//   tp:           the topology                                         
//   iC[1,2,3,4]:  the sectors of the sending cell to package up        
//   iD[1,2,3,4]:  the sectors of the receiving cell that well get information
//                 from sectors iC[1,2,3,4], respectively   
//   ct:           the type of information to transfer (0 = positions,  
//                 1 = velocities)                                      
//-----------------------------------------------------------------------------
static void ShareCellContents(cellgrid *CG, coord *crd, cell *C, int imove,
			      double Mcut, prmtop *tp, int iC1, int iC2,
			      int iC3, int iC4, int iD1, int iD2, int iD3,
			      int iD4, int ct, int nmsg)
{
  int i, needreim;
  int dloc[3];
  double reimv[3];
  double *nvec;
  cell *D;
  atomb *texport;
  atombv *tVexport;

  // Determine the destination cell and reimaging considerations
  for (i = 0; i < 3; i++) {
    dloc[i] = C->gbin[i];
    reimv[i] = 0.0;
  }
  if (C->gbin[imove] == 0) {
    dloc[imove] = CG->ng[imove]-1;
    reimv[imove] = 1.0;
    needreim = 1;
    RotateCrd(reimv, 1, crd->invU);
  }
  else {
    dloc[imove] = C->gbin[imove]-1;
    needreim = 0;
  }
  D = &CG->map[dloc[0]][dloc[1]][dloc[2]];
  nvec = crd->fcnorm.map[imove];
  if (C->CGRank == D->CGRank) {

    // If ct is 1 (velocities, not positions, are   
    // being transferred), skip this entirely.      
    // Previous positions are stored in arrays at   
    // the coordinate struct level; like velocities,
    // they do not need to be transferred.          
    if (ct == 1) {
      return;
    }
    if (iC1 >= 0) HandoffCell(C, D, iC1, iD1, imove, needreim, reimv,
			      Mcut, tp, crd->isortho, nvec);
    if (iC2 >= 0) HandoffCell(C, D, iC2, iD2, imove, needreim, reimv,
			      Mcut, tp, crd->isortho, nvec);
    if (iC3 >= 0) HandoffCell(C, D, iC3, iD3, imove, needreim, reimv,
			      Mcut, tp, crd->isortho, nvec);
    if (iC4 >= 0) HandoffCell(C, D, iC4, iD4, imove, needreim, reimv,
			      Mcut, tp, crd->isortho, nvec);

    // This cell:cell communication is complete
    return;
  }

  // If we're still here, this entails communication between     
  // processes.  Loop over four sectors and package the contents.
  // Exports are shunted to a cellgrid-wide array for grouped    
  // transfer.                                                   
  C->nexp = 1;
  if (ct < 2) {
    texport = C->pexport;
    C->pexport = &CG->pexport[nmsg][CG->nexp[nmsg]];
  }
  else {
    tVexport = C->Vexport;
    C->Vexport = &CG->Vexport[nmsg][CG->nexp[nmsg]];
  }
  if (iC1 >= 0) PackageCell(C, crd, imove, needreim, reimv, Mcut, iC1,
			    iD1, ct, tp);
  if (iC2 >= 0) PackageCell(C, crd, imove, needreim, reimv, Mcut, iC2,
			    iD2, ct, tp);
  if (iC3 >= 0) PackageCell(C, crd, imove, needreim, reimv, Mcut, iC3,
			    iD3, ct, tp);
  if (iC4 >= 0) PackageCell(C, crd, imove, needreim, reimv, Mcut, iC4,
			    iD4, ct, tp);

  // Store the number of exported atoms in the first
  // "atom" to avoid a call to MPI_Get_count later  
  if (ct < 2) {
    C->pexport[0].id = C->nexp;
    CG->nexp[nmsg] += C->nexp;
    C->pexport = texport;
  }
  else {
    C->Vexport[0].id = C->nexp;
    CG->nexp[nmsg] += C->nexp;
    C->Vexport = tVexport;
  }
}

//-----------------------------------------------------------------------------
// ShareCoordinates: have all cells share coordinates (which could mean 
//                   positions, but also velocities or even forces) for atoms
//                   which each cell owns with other cells which need the
//                   information for atoms in their import regions.
//
// Arguments:                                                           
//   CG:   the cell grid                                                
//   tj:   trjaectory controls (contains custom MPI data types)         
//   Mcut: the maximum direct space cutoff (electrostatic and van-der Waals  
//         forces)                                                
//   crd:  the coordinates                                              
//   tp:   the topology                                                 
//   ct:   the type of information to transfer (0 = positions, 1 = velocities,
//         2 = positions and previous positions)            
//-----------------------------------------------------------------------------
#ifdef MPI
void ShareCoordinates(cellgrid *CG, trajcon *tj, double Mcut, coord *crd,
		      prmtop *tp, int ct)
#else
void ShareCoordinates(cellgrid *CG, double Mcut, coord *crd, prmtop *tp,
		      int ct)
#endif
{
  int imove, h, i;
  cell *C;
  int j, nimp;
  atomb *timport;
  atombv *tVimport;
  ashr *cshr;
#ifdef MPI
  int ntag, nreq;
  MPI_Request* req;
  MPI_Status* stt;
#endif

  // The number of cells should be computed in all cases
#ifdef MPI
  req = (MPI_Request*)malloc((CG->nrecv+CG->nsend)*sizeof(MPI_Request));
  stt = (MPI_Status*)malloc((CG->nrecv+CG->nsend)*sizeof(MPI_Status));
#endif

  // If we are sending positions, all import sectors must
  // be cleared first; set the atom counts in all sectors
  // (except the primary sector) of all cells to zero.   
  if (ct == 0 || ct == 2) {
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      for (h = 1; h < 8; h++) {
	      C->nr[h] = 0;
      }
    }
  }

  // Loop over all cells with pulses in -X, -Y, and -Z
  for (imove = 0; imove < 3; imove++) {

    // Pointer to the communication plan
    cshr = &CG->DirCommPlan.mvshr[imove];

#ifdef MPI
    // Post receives
    for (i = 0; i < cshr->nrecv; i++) {
      ntag = cshr->recv[i].BaseID + DSP_CRD_SHARE + ct;
      if (ct < 2) {
	      MPI_Irecv(CG->import[i], CG->maximp[i], tj->MPI_ATOMB,
		              cshr->recv[i].partner, ntag, CG->dspcomm, &req[i]);
      }
      else {
        MPI_Irecv(CG->Vimport[i], CG->maximp[i], tj->MPI_ATOMV,
		              cshr->recv[i].partner, ntag, CG->dspcomm, &req[i]);
      }
    }

    // Post sends
    nreq = cshr->nrecv;
#endif
    for (i = 0; i < cshr->nsend; i++) {
      CG->nexp[i] = 0;
      for (j = 0; j < cshr->send[i].ncell; j++) {
	      C = &CG->data[cshr->send[i].cellpt[j]];
	      if (imove == 0) {
          ShareCellContents(CG, crd, C, imove, Mcut, tp, 0, -1, -1, -1,
			                      1, -1, -1, -1, ct, i);
        }
	      else if (imove == 1) {
	        ShareCellContents(CG, crd, C, imove, Mcut, tp, 0, 1, -1, -1,
			                      2, 3, -1, -1, ct, i);
        }
	      else if (imove == 2) {
	        ShareCellContents(CG, crd, C, imove, Mcut, tp, 0, 1, 2, 3,
			                      4, 5, 6, 7, ct, i);
	      }
      }
#ifdef MPI
      if (cshr->send[i].partner != CG->tid) {
	ntag = cshr->send[i].BaseID + DSP_CRD_SHARE + ct;
	if (ct < 2) {
	  MPI_Isend(CG->pexport[i], CG->nexp[i], tj->MPI_ATOMB,
	            cshr->send[i].partner, ntag, CG->dspcomm, &req[nreq]);
      	}
	else {
	  MPI_Isend(CG->Vexport[i], CG->nexp[i], tj->MPI_ATOMV,
	            cshr->send[i].partner, ntag, CG->dspcomm, &req[nreq]);
        }
        nreq++;
      }
#endif
    }
#ifdef MPI
    // Wait until all receives from this pulse are finished
    if (nreq > 0) {
      MPI_Waitall(nreq, req, stt);
    }
    MPI_Barrier(CG->dspcomm);
#endif

    // Unpack the results, if any
    for (i = 0; i < cshr->nrecv; i++) {
      nimp = 0;
      for (j = 0; j < cshr->recv[i].ncell; j++) {
	      C = &CG->data[cshr->recv[i].cellpt[j]];
	      if (ct < 2) {
	        timport = C->import;
	        C->import = &CG->import[i][nimp];
	        nimp += C->import[0].id;
	      }
	      else {
	        tVimport = C->Vimport;
	        C->Vimport = &CG->Vimport[i][nimp];
	        nimp += C->Vimport[0].id;
	      }
	      UnpackCellContents(C, ct, tp, crd);
        if (ct < 2) {
          C->import = timport;
        }
        else {
          C->Vimport = tVimport;
        }
      }
    }
  }

#ifdef MPI
  // Free allocate memory
  free(req);
  free(stt);
#endif
}

//-----------------------------------------------------------------------------
// CompactOrder: this function concatenates the orderings of atoms in one of
//               four sectors with those in other sectors.       
//                                                                      
// Arguments:                                                           
//   cordr: the padded list, containing gaps so that various sectors can hold
//          a potentially large number of atoms                    
//   astr:  the starting point of the first sector in cordr             
//   bstr:  the starting point of the second sector in cordr            
//   n:     the number of atoms in this segment                         
//-----------------------------------------------------------------------------
static void CompactOrder(int* cordr, int astr, int bstr, int n)
{
  int i;

  for (i = 0; i < n; i++) {
    cordr[astr+i] = cordr[bstr+i];
  }
}

#define VGTE 1
#define MULTIPLE 1
#include "Ordering.c"
#undef MULTIPLE
#define MULTIPLE 2
#include "Ordering.c"
#undef MULTIPLE
#define MULTIPLE 4
#include "Ordering.c"
#undef MULTIPLE
#define MULTIPLE 7
#include "Ordering.c"
#undef MULTIPLE
#undef VGTE

#define VGTE 0
#define MULTIPLE 1
#include "Ordering.c"
#undef MULTIPLE
#define MULTIPLE 2
#include "Ordering.c"
#undef MULTIPLE
#define MULTIPLE 4
#include "Ordering.c"
#undef MULTIPLE
#define MULTIPLE 7
#include "Ordering.c"
#undef MULTIPLE
#undef VGTE

//-----------------------------------------------------------------------------
// AtomRanger: an efficient method for determining, precisely once per atom in
//             the simulation, all of the bins an atom may occupy for
//             interactions between cell sectors.            
//                                                                      
// Arguments:                                                           
//   C:        the cell whose regions are being ordered                 
//   crd:      coordinate struct (for interface normals)                
//   dcinp:    the direct space control information (for cutoffs)       
//-----------------------------------------------------------------------------
static void AtomRanger(cell *C, coord *crd, dircon *dcinp)
{
  int i, j, i4, i4p1, i4p2, i8, maxatom;
  int *cnsr, *itmp, *itm2p, *cordr;
  double invcut, invEEcut;
  double **fcn;

  // Prepare counters for the numbers of atoms in each sub-region
  fcn = crd->fcnorm.map;
  maxatom = C->maxatom;
  cordr = C->ordr.data;
  cnsr = C->nsr;
  for (i = 0; i < 104; i++) {
    cnsr[i] = i*maxatom;
  }

  // Cutoff markers: invcut is the inverse direct space cutoff,
  // multiplied by 4 so that this factor can be used to turn   
  // displacements directly into bin numbers.  Similarly,      
  // invEEcut pertains to atoms with no Lennard-Jones          
  // properties when dcinp->Vcut > dcinp->Ecut.                
  invcut = 4.0*dcinp->invMcut;
  if (dcinp->Vcut > dcinp->Ecut) {
    invEEcut = 4.0*dcinp->invEcut;
    Divider7VgtE(C, 0, fcn[0], fcn[1], fcn[2], fcn[3], fcn[4], fcn[5], fcn[9],
		 C->midp, invcut, invEEcut, cordr, cnsr, &cnsr[8], &cnsr[16],
		 &cnsr[24], &cnsr[32], &cnsr[40], &cnsr[72], 0);
    Divider4VgtE(C, 1, fcn[13], fcn[6], fcn[7], fcn[10], C->midp, invcut,
		 invEEcut, cordr, &cnsr[4], &cnsr[48], &cnsr[56], &cnsr[80],
		 0);
    Divider4VgtE(C, 2, fcn[14], fcn[19], fcn[8], fcn[11], C->midp, invcut,
		 invEEcut, cordr, &cnsr[12], &cnsr[52], &cnsr[64],
		 &cnsr[88], 0);
    Divider4VgtE(C, 4, fcn[15], fcn[20], fcn[21], fcn[12], C->midp, invcut,
		 invEEcut, cordr, &cnsr[20], &cnsr[60], &cnsr[68],
		 &cnsr[96], 0);
    Divider2VgtE(C, 3, fcn[16], fcn[25], C->midp, invcut, invEEcut, cordr,
		 &cnsr[28], &cnsr[100], 0);
    Divider2VgtE(C, 5, fcn[17], fcn[24], C->midp, invcut, invEEcut, cordr,
		 &cnsr[36], &cnsr[92], 0);
    Divider2VgtE(C, 6, fcn[18], fcn[23], C->midp, invcut, invEEcut, cordr,
		 &cnsr[44], &cnsr[84], 0);
    Divider1VgtE(C, 7, fcn[22], C->midp, invcut, invEEcut, cordr, &cnsr[76],
		 0);
  }
  else {
    Divider7VeqE(C, 0, fcn[0], fcn[1], fcn[2], fcn[3], fcn[4], fcn[5], fcn[9],
		 C->midp, invcut, cordr, cnsr, &cnsr[8], &cnsr[16], &cnsr[24],
		 &cnsr[32], &cnsr[40], &cnsr[72], 0);
    Divider4VeqE(C, 1, fcn[13], fcn[6], fcn[7], fcn[10], C->midp, invcut,
		 cordr, &cnsr[4], &cnsr[48], &cnsr[56], &cnsr[80], 0);
    Divider4VeqE(C, 2, fcn[14], fcn[19], fcn[8], fcn[11], C->midp, invcut,
		 cordr, &cnsr[12], &cnsr[52], &cnsr[64], &cnsr[88], 0);
    Divider4VeqE(C, 4, fcn[15], fcn[20], fcn[21], fcn[12], C->midp, invcut,
		 cordr, &cnsr[20], &cnsr[60], &cnsr[68], &cnsr[96], 0);
    Divider2VeqE(C, 3, fcn[16], fcn[25], C->midp, invcut, cordr, &cnsr[28],
		 &cnsr[100], 0);
    Divider2VeqE(C, 5, fcn[17], fcn[24], C->midp, invcut, cordr, &cnsr[36],
		 &cnsr[92], 0);
    Divider2VeqE(C, 6, fcn[18], fcn[23], C->midp, invcut, cordr, &cnsr[44],
		 &cnsr[84], 0);
    Divider1VeqE(C, 7, fcn[22], C->midp, invcut, cordr, &cnsr[76], 0);
  }

  // Compact the atom orderings
  for (i = 1; i < 104; i++) {
    cnsr[i] -= i*maxatom;
  }
  for (i = 0; i < 26; i++) {
    i4 = 4*i;
    i4p1 = i4+1;
    i4p2 = i4+2;
    cnsr[i4p1] += cnsr[i4];
    cnsr[i4p2] += cnsr[i4p1];
    cnsr[i4+3] += cnsr[i4p2];
    CompactOrder(C->ordr.map[i4], cnsr[i4], maxatom, cnsr[i4p1]-cnsr[i4]);
    CompactOrder(C->ordr.map[i4], cnsr[i4p1], 2*maxatom,
                 cnsr[i4p2]-cnsr[i4p1]);
    CompactOrder(C->ordr.map[i4], cnsr[i4p2], 3*maxatom,
                 cnsr[i4+3]-cnsr[i4p2]);
  }

  // Copy the numberings of sector B's 
  // contents for fast retrieval later.
  for (i = 0; i < 13; i++) {
    i8 = 8*i;
    itmp = C->ordr.map[i8+4];
    itm2p = C->ordr.map[i8+5];
    for (j = 0; j < cnsr[i8+6]; j++) {
      itm2p[j] = itmp[j];
    }
    itm2p = C->ordr.map[i8+6];
    for (j = 0; j < cnsr[i8+5]; j++) {
      itm2p[j] = itmp[j];
    }
    itm2p = C->ordr.map[i8+7];
    for (j = 0; j < cnsr[i8+4]; j++) {
      itm2p[j] = itmp[j];
    }
  }
}

//-----------------------------------------------------------------------------
// ManageIntr: after partitioning atoms in each of two sectors into four
//             different sub-regions, we can exploit the known distances
//             between each sub-region to reduce the number of squared  
//             distances that must be evaluated.                        
//                                                                      
// Arguments:                                                           
//   C:        the cell for which interactions are computed             
//   s[A,B]:   the A and B sectors, interactions are computed A -> B    
//   o[A,B]:   indices into C->ordr array based on choices of sA and sB 
//   EFrc:     spline table of electrostatic forces (this and subsequent
//             parameters are passed down to various CompIntr functions)
//   dcinp:    direct space control parameters                          
//   tp:       the topology
//   sysUV:    energy and virial quantities
//   isortho:  flag to indicate whether the system is                   
//-----------------------------------------------------------------------------
static void ManageIntr(cell *C, int sA, int sB, int oA, int oB, FrcTab *EFrc,
		       dircon *dcinp, prmtop *tp, Energy *sysUV)
{
  int i, cstr;
  int *cordrA, *cordrB, *cnsr;

  // Set pointers
  cordrA = C->ordr.map[oA];
  cordrB = C->ordr.map[oB];
  cnsr = &C->nsr[oA];

  // A1:B(1,2,3,4), A2:B(1,2,3), A3:B(1,2), A4:B1
  for (i = 0; i < 4; i++) {
    cstr = (i == 0) ? 0 : cnsr[i-1];

    // Compute energies, virials, and forces
    if (sysUV->updateU == 2) {
      if (dcinp->Vcut > dcinp->Ecut) {
	if (i*0.25*dcinp->Vcut >= MINNB) {
	  CompIntrVgtEnrgvir(C, C->map[sA], C->map[sB], &cordrA[cstr],
			     &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i],
			     EFrc, dcinp, tp, sysUV);
	}
	else {
	  CompIntrVgtECRnrgvir(C, C->map[sA], C->map[sB], &cordrA[cstr],
			       &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i],
			       EFrc, dcinp, tp, sysUV);
	}
      }
      else {
	if (i*0.25*dcinp->Vcut >= MINNB) {
	  CompIntrVeqEnrgvir(C, C->map[sA], C->map[sB], &cordrA[cstr],
			     &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i],
			     EFrc, dcinp, tp, sysUV);
	}
	else {
	  CompIntrVeqECRnrgvir(C, C->map[sA], C->map[sB], &cordrA[cstr],
			       &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i],
			       EFrc, dcinp, tp, sysUV);
	}
      }
      continue;
    }

    // Compute energies and forces
    else if (sysUV->updateU == 1) {
      if (dcinp->Vcut > dcinp->Ecut) {
	if (i*0.25*dcinp->Vcut >= MINNB) {
	  CompIntrVgtEnrg(C, C->map[sA], C->map[sB], &cordrA[cstr],
			  &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i], EFrc,
			  dcinp, tp, sysUV);
	}
	else {
	  CompIntrVgtECRnrg(C, C->map[sA], C->map[sB], &cordrA[cstr],
			    &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i],
			    EFrc, dcinp, tp, sysUV);
	}
      }
      else {
	if (i*0.25*dcinp->Vcut >= MINNB) {
	  CompIntrVeqEnrg(C, C->map[sA], C->map[sB], &cordrA[cstr],
			  &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i], EFrc,
			  dcinp, tp, sysUV);
	}
	else {
	  CompIntrVeqECRnrg(C, C->map[sA], C->map[sB], &cordrA[cstr],
			    &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i],
			    EFrc, dcinp, tp, sysUV);
	}
      }
      continue;
    }

    // Compute energies alone
    else if (sysUV->updateU == -1) {
      if (dcinp->Vcut > dcinp->Ecut) {
	if (i*0.25*dcinp->Vcut >= MINNB) {
	  CompIntrVgtEnrgx(C, C->map[sA], C->map[sB], &cordrA[cstr],
			   &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i],
			   EFrc, dcinp, tp, sysUV);
	}
	else {
	  CompIntrVgtECRnrgx(C, C->map[sA], C->map[sB], &cordrA[cstr],
			     &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i],
			     EFrc, dcinp, tp, sysUV);
	}
      }
      else {
	if (i*0.25*dcinp->Vcut >= MINNB) {
	  CompIntrVeqEnrgx(C, C->map[sA], C->map[sB], &cordrA[cstr],
			   &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i],
			   EFrc, dcinp, tp, sysUV);
	}
	else {
	  CompIntrVeqECRnrgx(C, C->map[sA], C->map[sB], &cordrA[cstr],
			     &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i],
			     EFrc, dcinp, tp, sysUV);
	}
      }
      continue;
    }

    // If we're still here, only forces are necessary
    else {
      if (dcinp->Vcut > dcinp->Ecut) {
        if (i*0.25*dcinp->Vcut >= MINNB) {
	  CompIntrVgtE(C, C->map[sA], C->map[sB], &cordrA[cstr],
		       &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i], EFrc,
		       dcinp, tp);
	}
	else {
	  CompIntrVgtECR(C, C->map[sA], C->map[sB], &cordrA[cstr],
			 &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i], EFrc,
			 dcinp, tp);
	}
      }
      else {
	if (i*0.25*dcinp->Vcut >= MINNB) {
	  CompIntrVeqE(C, C->map[sA], C->map[sB], &cordrA[cstr],
		       &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i], EFrc,
		       dcinp, tp);
	}
	else {
	  CompIntrVeqECR(C, C->map[sA], C->map[sB], &cordrA[cstr],
			 &cordrB[i*C->maxatom], cnsr[i]-cstr, cnsr[7-i], EFrc,
			 dcinp, tp);
	}
      }
    }
  }
}

//-----------------------------------------------------------------------------
// DirectTriple2: perform a direct space nonbonded calculation for all atoms
//                in this cell using the "triple2," a.k.a. the "1/8th Shell"
//                method.                                 
//                                                                      
// Arguments:                                                           
//   C:       the cell                                                  
//   crd:     the coordinate array (for unit normal vectors)            
//   celldim: the dimensions of the cell (taken from the overarching cellgrid
//            structure)                                       
//   EFrc:    the electrostatic force spline table (also contains spline
//            coefficients for the potential, if desired)               
//   dcinp:   the direct space input parameters                         
//   tp:      the topology (passed in for rare cases in which we must refer to
//            the nonbonded exclusion list to determine whether an interaction
//            should actually be computed)               
//-----------------------------------------------------------------------------
void DirectTriple2(cell *C, coord *crd, FrcTab *EFrc, dircon *dcinp,
		   prmtop *tp, Energy *sysUV)
{
  // Compute tables for cell:cell interactions
  AtomRanger(C, crd, dcinp);
  if (dcinp->Vcut > dcinp->Ecut) {

    // First, do interactions within the primary sector/region
    if (sysUV->updateU == 2) {
      CompSameVgtECRnrgvir(C, C->data, C->nr[0], EFrc, dcinp, tp, sysUV);
    }
    else if (sysUV->updateU == 1) {
      CompSameVgtECRnrg(C, C->data, C->nr[0], EFrc, dcinp, tp, sysUV);
    }
    else if (sysUV->updateU == -1) {
      CompSameVgtECRnrgx(C, C->data, C->nr[0], EFrc, dcinp, tp, sysUV);
    }
    else {
      CompSameVgtECR(C, C->data, C->nr[0], EFrc, dcinp, tp);
    }

    // Next, do interactions across faces (3)
    ManageIntr(C, 0, 1,  0,   4, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 0, 2,  8,  12, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 0, 4, 16,  20, EFrc, dcinp, tp, sysUV);

    // Next, do interactions across edges (6)
    ManageIntr(C, 0, 3, 24,  28, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 0, 5, 32,  36, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 0, 6, 40,  44, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 1, 2, 48,  52, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 1, 4, 56,  60, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 2, 4, 64,  68, EFrc, dcinp, tp, sysUV);

    // Finally, do interactions across the center point (4)
    ManageIntr(C, 0, 7, 72,  76, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 1, 6, 80,  84, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 2, 5, 88,  92, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 4, 3, 96, 100, EFrc, dcinp, tp, sysUV);
  }
  else {

    // First, do interactions within the primary sector/region
    if (sysUV->updateU == 2) {
      CompSameVeqECRnrgvir(C, C->data, C->nr[0], EFrc, dcinp, tp, sysUV);
    }
    else if (sysUV->updateU == 1) {
      CompSameVeqECRnrg(C, C->data, C->nr[0], EFrc, dcinp, tp, sysUV);
    }
    else if (sysUV->updateU == -1) {
      CompSameVeqECRnrgx(C, C->data, C->nr[0], EFrc, dcinp, tp, sysUV);
    }
    else {
      CompSameVeqECR(C, C->data, C->nr[0], EFrc, dcinp, tp);
    }

    // Next, do interactions across faces (3)
    ManageIntr(C, 0, 1,  0,   4, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 0, 2,  8,  12, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 0, 4, 16,  20, EFrc, dcinp, tp, sysUV);

    // Next, do interactions across edges (6)
    ManageIntr(C, 0, 3, 24,  28, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 0, 5, 32,  36, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 0, 6, 40,  44, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 1, 2, 48,  52, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 1, 4, 56,  60, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 2, 4, 64,  68, EFrc, dcinp, tp, sysUV);

    // Finally, do interactions across the center point (4)
    ManageIntr(C, 0, 7, 72,  76, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 1, 6, 80,  84, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 2, 5, 88,  92, EFrc, dcinp, tp, sysUV);
    ManageIntr(C, 4, 3, 96, 100, EFrc, dcinp, tp, sysUV);
  }
}

//-----------------------------------------------------------------------------
// CellLRvdw: compute long-ranged van-der Waals interactions for all atoms
//            within a cell's primary sector.                     
//                                                                      
// Arguments:                                                           
//   C:       the cell                                                  
//-----------------------------------------------------------------------------
void CellLRvdw(cell *C, coord *crd, prmtop *tp, Energy *sysUV)
{
  int i, ityp;
  double invV;

  if (sysUV->updateU == 0) {
    return;
  }

  invV = 1.0/(crd->invU.map[0][0]*crd->invU.map[1][1]*crd->invU.map[2][2]);
  for (i = 0; i < C->nr[0]; i++) {
    ityp = C->data[i].lj;
    if (ityp >= 0) {
      sysUV->vdw6 += tp->lVDWc[3*ityp]*invV;
      sysUV->vdw12 += tp->lVDWc[3*ityp+1]*invV;
      if (sysUV->updateU == 2) {
	sysUV->Vir[0] += tp->lVDWc[3*ityp+2]*invV/3.0;
	sysUV->Vir[4] += tp->lVDWc[3*ityp+2]*invV/3.0;
	sysUV->Vir[8] += tp->lVDWc[3*ityp+2]*invV/3.0;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// LoadForceForExport: this funciton loads forces for an atoms from one of a
//                     cell's eight sectors for export to a sector of another
//                     cell.                                 
//                                                                      
// Arguments:                                                           
//   Aorig:        the original atom found in one of the cell's sectors 
//   Aexp:         pointer an element of the cell's export array        
//   dreg:         the destination region of the cell                   
//-----------------------------------------------------------------------------
static void LoadForceForExport(atomc *Aorig, atomb *Aexp, int dreg)
{
  Aexp->id = Aorig->id;
  Aexp->dreg = dreg;
  Aexp->loc[0] = Aorig->frc[0];
  Aexp->loc[1] = Aorig->frc[1];
  Aexp->loc[2] = Aorig->frc[2];
}

//-----------------------------------------------------------------------------
// ForcePkgLoop: this function wraps the decision to loop over all atoms
//               within a region of a cell, or just extra points.       
//                                                                      
// Arguments:                                                           
//   C:        the cell                                                 
//   iC:       the sector of the cell to be packaged, negative if only extra
//             points are to be packaged                          
//   iD:       the destination region of some other cell where the atom forces
//             will be going                                     
//   nexp:     the number of atoms packaged thus far (also returned)    
//   tp:       the topology                                             
//-----------------------------------------------------------------------------
static int ForcePkgLoop(cell *C, int iC, int iD, int nexp, prmtop *tp)
{
  int i;

  // Normal case: package atoms with mass only
  if (iC >= 0) {
    for (i = 0; i < C->nr[iC]; i++) {
      if (tp->Masses[C->map[iC][i].id] >= 1.0e-8) {
	LoadForceForExport(&C->map[iC][i], &C->pexport[nexp], iD);
	nexp++;
      }
    }
  }

  // Special cases: package extra points only
  else {
    iC *= -1;

    // Extra points only for sector 0 (note that all 
    // points would normally be packaged in sector 0)
    if (iC > 7) {
      iC = 0;
    }

    for (i = 0; i < C->nr[iC]; i++) {
      if (tp->Masses[C->map[iC][i].id] < 1.0e-8) {
	LoadForceForExport(&C->map[iC][i], &C->pexport[nexp], iD);
	nexp++;
      }
    }
  }

  return nexp;
}

//-----------------------------------------------------------------------------
// PackageCellForces: this function wraps four loops for sharing forces on
//                    atoms between sectors of two adjacent cells.  This
//                    routine intercepts iCN >= 8, which implies that iCN ->
//                    iDN transfer (where N = 1, 2, 3, or 4) is to be skipped
//                    entirely.                        
//                                                                      
// Arguments:                                                           
//   C:           the cell in which atoms are being packaged            
//   iC[1,2,3,4]: the sectors of cell C to package                      
//   iD[1,2,3,4]: the sectors of cell D receiving forces from sectors   
//                iC[1,2,3,4] of cell C, respectively                   
//   tp:          the topology                                          
//-----------------------------------------------------------------------------
static void PackageCellForces(cell *C, int iC1, int iC2, int iC3, int iC4,
			      int iD1, int iD2, int iD3, int iD4, prmtop *tp)
{
  if (iC1 < 8) {
    C->nexp = ForcePkgLoop(C, iC1, iD1, 1, tp);
  }
  if (iC2 < 8) {
    C->nexp = ForcePkgLoop(C, iC2, iD2, C->nexp, tp);
  }
  if (iC3 < 8) {
    C->nexp = ForcePkgLoop(C, iC3, iD3, C->nexp, tp);
  }
  if (iC4 < 8) {
    C->nexp = ForcePkgLoop(C, iC4, iD4, C->nexp, tp);
  }
}

//-----------------------------------------------------------------------------
// ForceUnpackLoop: this function wraps the loop over a single region of a cell
//                  for unpacking imported forces, and increments the import
//                  counter appropriately.                   
//                                                                      
// Arguments:                                                           
//   D:           the cell in which atoms are being unpacked            
//   iD:          the region of the cell which is being unpacked        
//-----------------------------------------------------------------------------
static int ForceUnpackLoop(cell *D, int iD, int impcon)
{
  int i;
  atomb *Dimp;
  atomc *Datms;

  // There may be no imported atoms; return
  // immediately if that is the case        
  if (D->nimp == 1) {
    return impcon;
  }

  // Set pointer to imported atom struct; it will be
  // updated less frequently than the atom counter  
  Dimp = &D->import[impcon];

  // Bail out if the destination region of the current
  // atom is not the region we're looping over        
  if (Dimp->dreg != iD) {
    return impcon;
  }

  // Loop over all atoms in the region, compare them
  // to the current imported atom, and increment the
  // imported atom counter as appropriate           
  Datms = D->map[iD];
  for (i = 0; i < D->nr[iD]; i++) {

    // It is possible that the import buffer will contain 
    // atoms that are not present in the region.  However,
    // both the import region and the cell's atom array   
    // are arranged in increasing order of atom ID.       
    while (Datms[i].id > Dimp->id) {
      impcon++;
      Dimp = &D->import[impcon];
      if (impcon == D->nimp || Dimp->dreg != iD) {
	return impcon;
      }
    }
    if (Datms[i].id == Dimp->id) {
      Datms[i].frc[0] += Dimp->loc[0];
      Datms[i].frc[1] += Dimp->loc[1];
      Datms[i].frc[2] += Dimp->loc[2];
      impcon++;
      Dimp = &D->import[impcon];
      if (impcon == D->nimp || Dimp->dreg != iD) {
	return impcon;
      }
    }
  }

  // All atoms in this cell region have been searched.
  // If there are any more imports for this region,   
  // they cannot be placed and so must be skipped.    
  while (impcon < D->nimp && Dimp->dreg == iD) {
    impcon++;
    Dimp = &D->import[impcon];
  }

  return impcon;
}

//-----------------------------------------------------------------------------
// UnpackCellForces: this function wraps four loops for sharing forces on atoms
//                   between sectors of two adjacent cells; it performs the
//                   reciprocal task of unpacking forces.  
//                                                                      
// Arguments:                                                           
//   D:           the cell in which atoms are being unpacked            
//   iD[1,2,3,4]: the sectors of cell D receiving forces from entries in
//                the cell's input buffer                               
//-----------------------------------------------------------------------------
static void UnpackCellForces(cell *D, int iD1, int iD2, int iD3, int iD4)
{
  int impcon;

  D->nimp = D->import[0].id;
  impcon = ForceUnpackLoop(D, iD1, 1);
  if (impcon == D->nimp) {
    return;
  }
  impcon = ForceUnpackLoop(D, iD2, impcon);
  if (impcon == D->nimp) {
    return;
  }
  impcon = ForceUnpackLoop(D, iD3, impcon);
  if (impcon == D->nimp) {
    return;
  }
  impcon = ForceUnpackLoop(D, iD4, impcon);
}

//-----------------------------------------------------------------------------
// ForceUnpackSwitch: this function wraps the UnpackCellForces function.
//                                                                      
// Arguments:                                                           
//   N:        the cell receiving information                           
//   imove:    the X, Y, or Z direction of the move                     
//   spread:   0 for +imove (focusing), 1 for -imove (spreading EP)     
//-----------------------------------------------------------------------------
static void ForceUnpackSwitch(cell *N, int imove, int spread)
{
  if (imove == 0) {
    if (spread == 0) {
      UnpackCellForces(N, 0, 2, 4, 6);
    }
    else {
      UnpackCellForces(N, 1, 3, 5, 7);
    }
  }
  else if (imove == 1) {
    if (spread == 0) {
      UnpackCellForces(N, 0, 1, 4, 5);
    }
    else {
      UnpackCellForces(N, 2, 3, 6, 7);
    }
  }
  else if (imove == 2) {
    if (spread == 0) {
      UnpackCellForces(N, 0, 1, 2, 3);
    }
    else {
      UnpackCellForces(N, 4, 5, 6, 7);
    }
  }
}

//-----------------------------------------------------------------------------
// GatherForces: forces on extra points must first be broadcast to all cell
//               sectors so that they can be redistributed to atoms with mass
//               before the final force merger.               
//                                                                      
// Arguments:                                                           
//   cellid:  number of the cell of interest, implies its location in   
//            the cell grid                                             
//   imove:   the direction of the move (0 = X, 1 = Y, 2 = Z)           
//   CG:      the cell grid                                             
//   tp:      the topology                                              
//   tj:      for MPI compilation, contains MPI process map information 
//   spread:  0 if this routine is to collect forces into the primary   
//            sector, 1 if this routine is to propagate forces into all 
//            sectors                                                   
//   EPonly:  1 if only extra points are to be considered, 0 if extra   
//            are to be excluded from the merger                        
//-----------------------------------------------------------------------------
static void GatherForces(int cellid, int imove, cellgrid *CG, prmtop *tp,
			 int spread, int EPonly, int nmsg)
{
  int i;
  int oL[3], nL[3];
  cell *O, *N;
  atomb *texport;

  // The original cell
  O = &CG->data[cellid];

  // Destination cell and reimaging considerations
  for (i = 0; i < 3; i++) {
    oL[i] = O->gbin[i];
    nL[i] = oL[i];
  }
  if (spread == 0) {
    if (oL[imove] == CG->ng[imove] - 1) {
      nL[imove] = 0;
    }
    else {
      nL[imove] = oL[imove] + 1;
    }
  }
  else {
    if (oL[imove] == 0) {
      nL[imove] = CG->ng[imove] - 1;
    }
    else {
      nL[imove] = oL[imove] - 1;
    }
  }
  N = &CG->map[nL[0]][nL[1]][nL[2]];

  // If the destination cell is on a different process,
  // swap the O cell's export pointer for a segment of 
  // the cell grid's pooled export buffer.  Otherwise, 
  // set the O cell's export pointer to the N cell's   
  // import buffer and perform the transfer implicitly.
  texport = O->pexport;
  if (N->CGRank != O->CGRank) {
    O->pexport = &CG->pexport[nmsg][CG->nexp[nmsg]];
  }
  else {
    O->pexport = N->import;
  }

  // X-direction message
  if (imove == 0) {
    if (spread == 0) {
      if (EPonly == 1) {
	PackageCellForces(O, -1, -3, -5, -7, 0, 2, 4, 6, tp);
      }
      else {
	PackageCellForces(O, 1, 3, 5, 7, 0, 2, 4, 6, tp);
      }
    }
    else {
      PackageCellForces(O, -8, 8, 8, 8, 1, 3, 5, 7, tp);
    }
  }

  // Y-direction message
  else if (imove == 1) {
    if (spread == 0) {
      if (EPonly == 1) {
	PackageCellForces(O, -2, 8, -6, 8, 0, 1, 4, 5, tp);
      }
      else {
	PackageCellForces(O, 2, 8, 6, 8, 0, 1, 4, 5, tp);
      }
    }
    else {
      PackageCellForces(O, -8, -1, 8, 8, 2, 3, 6, 7, tp);
    }
  }

  // Z-direction message
  else if (imove == 2) {
    if (spread == 0) {
      if (EPonly == 1) {
	PackageCellForces(O, -4, 8, 8, 8, 0, 1, 2, 3, tp);
      }
      else {
	PackageCellForces(O, 4, 8, 8, 8, 0, 1, 2, 3, tp);
      }
    }
    else {
      PackageCellForces(O, -8, -1, -2, -3, 4, 5, 6, 7, tp);
    }
  }

  // Record the number of atom forces being transferred
  O->pexport[0].id = O->nexp;

  // If the origin and destination cells are on different
  // processes, increment the cell grid's pooled buffer  
  // count and return the O cell's export pointer to its 
  // local export buffer.  Otherwise, the transfer has   
  // already occurred and unpacking can proceed at once. 
  O->pexport = texport;
  if (N->CGRank != O->CGRank) {
    CG->nexp[nmsg] += O->nexp;
  }
  else {
    ForceUnpackSwitch(N, imove, spread);
  }
}

//-----------------------------------------------------------------------------
// GatherForcePulse: encapsulates a lot of code that is used in three different
//                   variants throughout the force gathering routine.                                           
//                                                                      
// Arguments:                                                           
//   CG:      the cell grid                                             
//   tp:      the topology                                              
//   tj:      trajectory control information (holds the process rank and
//            defined MPI types)                                        
//   spread:  0 if this routine is to collect forces into the primary   
//            sector, 1 if this routine is to propagate forces into all 
//            sectors                                                   
//   EPonly:  1 if only extra points are to be considered, 0 if extra   
//            points are to be excluded from the merger                 
//-----------------------------------------------------------------------------
#ifdef MPI
static void GatherForcePulse(cellgrid *CG, prmtop *tp, trajcon *tj, int spread,
			     int EPonly, MPI_Request* req, MPI_Status* stt)
#else
static void GatherForcePulse(cellgrid *CG, prmtop *tp, int spread, int EPonly)
#endif
{
  int i, j, imove;
  ashr *cshr;
  cell *N;

  for (imove = 0; imove < 3; imove++) {
    cshr = (spread == 1) ? &CG->DirCommPlan.mvshr[imove] :
      &CG->DirCommPlan.frcmg[imove];

    // Post receives
#ifdef MPI
    for (i = 0; i < cshr->nrecv; i++) {
      MPI_Irecv(CG->import[i], CG->maximp[i], tj->MPI_ATOMB,
		cshr->recv[i].partner,
		cshr->recv[i].BaseID + DSP_FRC_MERGE + 2*spread + EPonly,
		CG->dspcomm, &req[i]);
    }
    int nreq = cshr->nrecv;
#endif

    // Post sends
    for (i = 0; i < cshr->nsend; i++) {
      CG->nexp[i] = 0;
      for (j = 0; j < cshr->send[i].ncell; j++) {
	GatherForces(cshr->send[i].cellpt[j], imove, CG, tp, spread, EPonly,
		     i);
      }
#ifdef MPI
      if (cshr->send[i].partner != CG->tid) {
	MPI_Isend(CG->pexport[i], CG->nexp[i], tj->MPI_ATOMB,
		  cshr->send[i].partner,
		  cshr->send[i].BaseID + DSP_FRC_MERGE + 2*spread + EPonly,
		  CG->dspcomm, &req[nreq]);
	nreq++;
      }
#endif
    }
#ifdef MPI
    if (nreq > 0) {
      MPI_Waitall(nreq, req, stt);
    }
    MPI_Barrier(CG->dspcomm);
#endif
    for (i = 0; i < cshr->nrecv; i++) {
      int nimp = 0;
      atomb *timport;
      for (j = 0; j < cshr->recv[i].ncell; j++) {
	N = &CG->data[cshr->recv[i].cellpt[j]];
	timport = N->import;
	N->import = &CG->import[i][nimp];
	ForceUnpackSwitch(N, imove, spread);
	nimp += N->import[0].id;
	N->import = timport;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// MergeCellForces: merge the results of a force calculation.  After this is
//                  done, the processors should be ready to integrate the
//                  equations of motion for particles in each cell's primary
//                  sector.                         
//                                                                      
// Arguments:                                                           
//   CG:      the cell grid                                             
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
#ifdef MPI
void MergeCellForces(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj)
#else
void MergeCellForces(cellgrid *CG, coord *crd, prmtop *tp)
#endif
{
  int i, j, k;
  cell *C;
  atomc *catm;

#ifdef MPI
  MPI_Request* req;
  MPI_Status* stt;
  req = (MPI_Request*)malloc(2*CG->nthreads*sizeof(MPI_Request));
  stt = (MPI_Status*)malloc(2*CG->nthreads*sizeof(MPI_Status));
#endif

  // Broadcast extra point forces
  if (tp->nxtrapt > 0) {

    // Pulses in X, Y, and Z to accumulate forces on 
    // extra points at the primary sector of one cell
#ifdef MPI
    GatherForcePulse(CG, tp, tj, 0, 1, req, stt);
#else
    GatherForcePulse(CG, tp, 0, 1);
#endif

    // Zero extra point forces on all but the primary sectors
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      for (j = 1; j < 8; j++) {
	catm = C->map[j];
	for (k = 0; k < C->nr[j]; k++) {
	  if (tp->Masses[catm[k].id] < 1.0e-8) {
	    catm[k].frc[0] = 0.0;
	    catm[k].frc[1] = 0.0;
	    catm[k].frc[2] = 0.0;
	  }
	}
      }
    }

    // Pulses in X, Y, and Z to propagate forces on
    // extra points across all sectors of all cells
#ifdef MPI
    GatherForcePulse(CG, tp, tj, 1, 1, req, stt);
#else
    GatherForcePulse(CG, tp, 1, 1);
#endif

    // Update the atom GPS in preparation
    // for lots of referencing.          
    UpdateCellGPS(CG);

    // Loop over all cells and transfer forces
    // from extra points to frame atoms       
    for (i = 0; i < CG->MyCellCount; i++) {
      CellXferEPForces(&CG->data[CG->MyCellDomain[i]], tp, crd, CG);
    }
  }

  // Roundup of all forces on atoms with mass
#ifdef MPI
  GatherForcePulse(CG, tp, tj, 0, 0, req, stt);
#else
  GatherForcePulse(CG, tp, 0, 0);
#endif

#ifdef MPI
  // Free allocated memory
  free(req);
  free(stt);
#endif
}

//-----------------------------------------------------------------------------
// MapCellForcesToAtoms: this function, which should not be called during     
//                       typical MD steps, takes forces from the atom lists in 
//                       cells and maps them back to the ordered list of atoms.
//
// Arguments:                                                           
//   CG:   the cell grid                                                
//   crd:  the coordinates structure (the ordered list of atoms)        
//-----------------------------------------------------------------------------
void MapCellForcesToAtoms(cellgrid *CG, coord *crd)
{
  int i, j, k, aid;
  cell *C;

  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < C->nr[0]; j++) {
      aid = 3*C->data[j].id;
      for (k = 0; k < 3; k++) {
	crd->frc[aid+k] = C->data[j].frc[k];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// MapListForcesToCells: this function, which should not be called during     
//                       typical MD steps, takes forces from the ordered lists 
//                       of atoms and maps them back to the cell grid.  It is
//                       used primarily when the coord struct is buffering
//                       forces for the cell grid.
//
// Arguments:                                                           
//   CG:   the cell grid                                                
//   crd:  the coordinates structure (the ordered list of atoms)        
//-----------------------------------------------------------------------------
void MapListForcesToCells(cellgrid *CG, coord *crd)
{
  int i, j, k, aid;
  cell *C;

  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < C->nr[0]; j++) {
      aid = 3*C->data[j].id;
      for (k = 0; k < 3; k++) {
	C->data[j].frc[k] = crd->frc[aid+k];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SingleAtomOrder: in a cell's atom map, there may be a single atom at the end
//                  of the primary sector atom list that needs to be placed to
//                  maintain an ascending order of atom IDs.  This function
//                  makes things right again.       
//
// Arguments:
//   D:       the cell                                                  
//-----------------------------------------------------------------------------
static void SingleAtomOrder(cell *D)
{
  int i, lastid, insertpt;
  atomc abuf;
  atomc *datm;

  const int natm = D->nr[0];
  datm = D->data;  
  lastid = datm[natm-1].id;

  // Bail out if the list is already in order
  if (natm == 1 || lastid > datm[natm-2].id) {
    return;
  }

  // Find the insertion point
  insertpt = natm-1;
  for (i = 0; i < natm-1; i++) {
    if (lastid < datm[i].id) {
      insertpt = i;
      break;
    }
  }
  
  // Shift the list
  abuf = datm[natm-1];
  for (i = natm-1; i > insertpt; i--) {
    datm[i] = datm[i-1];
  }
  datm[insertpt] = abuf;
}

//-----------------------------------------------------------------------------
// ExpungeAtomBlanks: after atoms migrate out of this cell, there will be holes
//                    left in the primary sector atom list.  Find and eliminate
//                    the holes.                     
//
// Arguments:
//   C:      the cell to operate on
//-----------------------------------------------------------------------------
static void ExpungeAtomBlanks(cell *C)
{
  int h, i;
  atomc *catm;

  catm = C->data;
  h = 0;

  for (i = 0; i < C->nr[0]; i++) {
    if (catm[i].id >= 0) {
      if (h < i) {
	catm[h] = catm[i];
      }
      h++;
    }
  }
  C->nr[0] = h;
}

//-----------------------------------------------------------------------------
// C2CDirectTransfer: transfer atoms that have migrated out of cell C into cell
//                    D.  This assumes that the same process controls cells C
//                    and D.                           
//                                                                      
// Arguments:                                                           
//   C,D:       the cells of interest                                   
//   mflag:     direction of the move                                   
//   vv:        dimension of the move (X = 0, Y = 1, Z = 2)             
//   clim:      the relevant boundary of the cell's primary region in   
//              the dimension of the move                               
//-----------------------------------------------------------------------------
static void C2CDirectTransfer(cell *C, cell *D, int mflag, int vv, double clim,
			      double *reimv, prmtop *tp)
{
  int i;
  atomc *catm, *datm;

  // Loop over all atoms in cell C and
  // transfer those leaving to cell D 
  catm = C->data;
  for (i = 0; i < C->nr[0]; i++) {
    if (tp->Masses[catm[i].id] < 1.0e-8) {
      continue;
    }
    if ((mflag == -1 && catm[i].loc[vv] < clim) ||
	(mflag == 1 && catm[i].loc[vv] >= clim)) {
      datm = &D->data[D->nr[0]];
      *datm = catm[i];
      datm->loc[0] += reimv[0];
      datm->loc[1] += reimv[1];
      datm->loc[2] += reimv[2];
      D->nr[0] += 1;
      SingleAtomOrder(D);

      // Mark this atom as having left
      catm[i].id = -1;
    }
  }

  // Clean up in cell C
  ExpungeAtomBlanks(C);
}

//-----------------------------------------------------------------------------
// C2CProcessTransfer: transfer atoms that have migrated out of cell C into
//                     cell D.  Different processes control cells C and D.
//
// Arguments:                                                           
//   C,D:       the cells of interest                                   
//   mflag:     direction of the move                                   
//   vv:        dimension of the move (X = 0, Y = 1, Z = 2)             
//   clim:      the relevant boundary of the cell's primary region in   
//              the dimension of the move                               
//   reimv:     the reimaging vector                                    
//   crd:       coordinates array (for absolute system locations of all 
//              sites, velocities, previous positions and velocities,   
//              and previous forces on atoms)                           
//   tp:        the topology                                            
//   tj:        trajectory control information for MPI process maps     
//   req:       pointer to single MPI request for asynch communication  
//   nreq:      number of MPI requests thus far on this process         
//-----------------------------------------------------------------------------
static void C2CProcessTransfer(cell *C, int mflag, int vv, double clim,
			       double *reimv, coord *crd)
{
  int i, j, ja3, nexp, atmid;
  atomc *catm;
  atombx *cXexp;

  // Loop over all atoms in cell C
  // and load the export buffer   
  nexp = 1;
  catm = C->data;
  cXexp = C->Xexport;
  for (i = 0; i < C->nr[0]; i++) {
    if ((mflag == -1 && catm[i].loc[vv] < clim) ||
        (mflag == 1 && catm[i].loc[vv] >= clim)) {

      // Load export buffer
      atmid = catm[i].id;
      cXexp[nexp].id = atmid;
      cXexp[nexp].dreg = 0;
      for (j = 0; j < 3; j++) {
	cXexp[nexp].loc[j] = catm[i].loc[j] + reimv[j];
	ja3 = 3*atmid+j;
	cXexp[nexp].vel[j] = crd->vel[ja3];
	cXexp[nexp].sysloc[j] = crd->loc[ja3];
	cXexp[nexp].prvloc[j] = crd->prvloc[ja3];
	cXexp[nexp].prvfrc[j] = crd->prvfrc[ja3];
	cXexp[nexp].prvvel[j] = crd->prvvel[ja3];
      }
      nexp++;

      // Mark this atom as having left
      catm[i].id = -1;
    }
  }

  // Clean up in cell C
  ExpungeAtomBlanks(C);

  // Transfer the buffer
  C->nexp = nexp;
  cXexp[0].id = nexp;
}

#ifdef MPI
//-----------------------------------------------------------------------------
// UnpackAtomAllInfo: unpacks information from an imported atombx array.  This
//                    finishes the process of atom migration from one cell to
//                    another.                              
//
// Arguments:
//   D:       the destination cell with imports                         
//   crd:     the coordinates (receives a great deal of the information)
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
static void UnpackAtomAllInfo(cell *D, coord *crd, prmtop *tp)
{
  int i, j, ja3, atmid;
  atomc *datm, *datmp;
  atombx *dXimp;

  // Loop over the imports to cell D and load it up
  datm = D->data;
  for (i = 1; i < D->Ximport[0].id; i++) {
    datmp = &datm[D->nr[0]];
    dXimp = &D->Ximport[i];
    atmid = dXimp->id;
    datmp->id = atmid;
    datmp->lj = tp->LJIdx[atmid];
    datmp->q = tp->Charges[atmid];
    for (j = 0; j < 3; j++) {
      datmp->loc[j] = dXimp->loc[j];
      ja3 = 3*atmid+j;
      crd->vel[ja3] = dXimp->vel[j];
      crd->loc[ja3] = dXimp->sysloc[j];
      crd->prvloc[ja3] = dXimp->prvloc[j];
      crd->prvfrc[ja3] = dXimp->prvfrc[j];
      crd->prvvel[ja3] = dXimp->prvvel[j];
    }
    D->nr[0] += 1;
    SingleAtomOrder(D);
  }
}
#endif

//-----------------------------------------------------------------------------
// ReleaseAtoms: when atoms move out of a cell, they are released to other
//               cells.  This could require that each cell communicate with
//               twenty-six nearest neighbors, but by releasing atoms in X, Y,
//               and Z sequentially there are only six messages per cell.
//                                                                      
// Arguments:                                                           
//   CG:      the cell grid                                             
//   crd:     the coordinates                                           
//   tp:      the topology                                              
//   cellid:  the number of the cell within the cell grid               
//   vv:      the direction along which to release (+X = 1, +Y = 2,     
//              +Z = 3, -X = -1, -Y = -2, -Z = -3)                      
//-----------------------------------------------------------------------------
static void ReleaseAtoms(cellgrid *CG, coord *crd, prmtop *tp, int cellid,
			 int vv, int nmsg)
{
  int i, mflag;
  int cijk[3], dijk[3];
  double clim;
  double reimv[3];
  double *celldim;
  cell *C, *D;
  atombx *tXexport;

  // Determine the cell limits
  celldim = CG->celldim;
  cijk[0] = CG->data[cellid].gbin[0];
  cijk[1] = CG->data[cellid].gbin[1];
  cijk[2] = CG->data[cellid].gbin[2];
  C = &CG->map[cijk[0]][cijk[1]][cijk[2]];
  if (vv < 0) {
    vv *= -1;
    vv -= 1;
    mflag = -1;
    if (crd->isortho == 0) {
      clim = cijk[vv]/CG->dbng[vv];
    }
    else {
      clim = C->orig[vv];
    }
  }
  else {
    vv -= 1;
    mflag = 1;
    if (crd->isortho == 0) {
      clim = (cijk[vv]+1.0)/CG->dbng[vv];
    }
    else {
      clim = C->orig[vv] + celldim[vv];
    }
  }

  // Determine the destination cell
  dijk[0] = cijk[0];
  dijk[1] = cijk[1];
  dijk[2] = cijk[2];
  dijk[vv] += mflag;
  if (dijk[vv] == CG->ng[vv]) {
    dijk[vv] = 0;
  }
  else if (dijk[vv] == -1) {
    dijk[vv] = CG->ng[vv] - 1;
  }
  D = &CG->map[dijk[0]][dijk[1]][dijk[2]];

  // Determine re-imaging considerations
  for (i = 0; i < 3; i++) {
    reimv[i] = 0.0;
    if (cijk[i] == 0 && dijk[i] == CG->ng[i]-1 && mflag == -1) {
      reimv[i] = 1.0;
    }
    else if (cijk[i] == CG->ng[i]-1 && dijk[i] == 0 && mflag == 1) {
      reimv[i] = -1.0;
    }
  }
  if (crd->isortho == 1) {
    RotateCrd(reimv, 1, crd->invU);
  }

  // Loop over all atoms in the cell's primary sector:
  // which ones have moved outside the cell's limits? 
  if (C->CGRank == D->CGRank) {
    C2CDirectTransfer(C, D, mflag, vv, clim, reimv, tp);
  }
  else {
    tXexport = C->Xexport;
    C->Xexport = &CG->Xexport[nmsg][CG->nexp[nmsg]];
    C2CProcessTransfer(C, mflag, vv, clim, reimv, crd);
    CG->nexp[nmsg] += C->Xexport[0].id;
    C->Xexport = tXexport;
  }
}

//-----------------------------------------------------------------------------
// UpdateCells: this function wraps the ReleaseAtoms and MigrateAtoms   
//              functions (above).                                      
//                                                                      
// Arguments:                                                           
//   CG:      the cell grid                                             
//   crd:     the coordinates (for velocities, system locations, etc.)  
//   tp:      the topology                                              
//-----------------------------------------------------------------------------
#ifdef MPI
void UpdateCells(cellgrid *CG, coord *crd, prmtop *tp, trajcon*tj)
#else
void UpdateCells(cellgrid *CG, coord *crd, prmtop *tp)
#endif
{
  int h, i, j, imove, cidx;
  double dx, dy, dz;
  double *Utmp;
  ashr* cshr;
  atomc *catm;
  cell *C;

#ifdef MPI
  int maxmsg;
  int nreq[2];
  MPI_Request* req;
  MPI_Request *tmpreq;
  MPI_Status* stt;
  MPI_Status *tmpstt;
  maxmsg = CG->nrecv + CG->nsend;
  req = (MPI_Request*)malloc(2*maxmsg*sizeof(MPI_Request));
  stt = (MPI_Status*)malloc(2*maxmsg*sizeof(MPI_Status));
#endif

  // If this is not an orthorhombic cell, do one   
  // pass over all cells to transform coordinates  
  // of atoms in the primary sectors into box space
  if (crd->isortho == 0) {
    Utmp = crd->U.data;
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      catm = C->data;
      for (j = 0; j < C->nr[0]; j++) {
	dx = catm[j].loc[0]*Utmp[0] + catm[j].loc[1]*Utmp[1] +
	  catm[j].loc[2]*Utmp[2];
	dy = catm[j].loc[0]*Utmp[3] + catm[j].loc[1]*Utmp[4] +
	  catm[j].loc[2]*Utmp[5];
	dz = catm[j].loc[0]*Utmp[6] + catm[j].loc[1]*Utmp[7] +
	  catm[j].loc[2]*Utmp[8];
	catm[j].loc[0] = dx;
	catm[j].loc[1] = dy;
	catm[j].loc[2] = dz;
      }
    }
  }

  // Loop over all cells and move atoms
  cshr = (ashr*)malloc(2*sizeof(ashr));
  for (imove = 1; imove < 4; imove++) {
    cshr[0] = CG->DirCommPlan.frcmg[imove-1];
    cshr[1] = CG->DirCommPlan.mvshr[imove-1];
    cidx = 0;
    for (h = 1; h > -2; h -= 2) {
#ifdef MPI
      // Post receives
      tmpreq = &req[cidx*maxmsg];
      tmpstt = &stt[cidx*maxmsg];
      for (i = 0; i < cshr[cidx].nrecv; i++) {
	MPI_Irecv(CG->Ximport[i], CG->maximp[i], tj->MPI_ATOMX,
		  cshr[cidx].recv[i].partner,
		  cshr[cidx].recv[i].BaseID + DSP_UPDATE + cidx,
		  CG->dspcomm, &tmpreq[i]);
      }
      nreq[cidx] = cshr[cidx].nrecv;
#endif
      // Post sends
      for (i = 0; i < cshr[cidx].nsend; i++) {
	CG->nexp[i] = 0;
	for (j = 0; j < cshr[cidx].send[i].ncell; j++) {
	  ReleaseAtoms(CG, crd, tp, cshr[cidx].send[i].cellpt[j], h*imove, i);
	}
#ifdef MPI
	if (cshr[cidx].send[i].partner != CG->tid) {
	  MPI_Isend(CG->Xexport[i], CG->nexp[i], tj->MPI_ATOMX,
		    cshr[cidx].send[i].partner,
		    cshr[cidx].send[i].BaseID + DSP_UPDATE + cidx,
		    CG->dspcomm, &tmpreq[nreq[cidx]]);
	  nreq[cidx] += 1;
	}
#endif
      }
#ifdef MPI
      MPI_Waitall(nreq[cidx], tmpreq, tmpstt);
      for (i = 0; i < cshr[cidx].nrecv; i++) {
	int nimp = 0;
	atombx *tXimport;
	cell *D;
	for (j = 0; j < cshr[cidx].recv[i].ncell; j++) {
	  D = &CG->data[cshr[cidx].recv[i].cellpt[j]];
	  tXimport = D->Ximport;
	  D->Ximport = &CG->Ximport[i][nimp];
	  nimp += D->Ximport[0].id;
	  UnpackAtomAllInfo(D, crd, tp);
	  D->Ximport = tXimport;
	}
      }
      MPI_Barrier(CG->dspcomm);
#endif
      cidx++;
    }
  }

  // Free allocated memory
#ifdef MPI
  free(req);
  free(stt);
#endif
  free(cshr);

  // If this is not an orthorhombic cell, do a pass  
  // over all cells to transform coordinates of atoms
  // in the primary sectors back into real space     
  if (crd->isortho == 0) {
    Utmp = crd->invU.data;
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      catm = C->data;
      for (j = 0; j < C->nr[0]; j++) {
        dx = catm[j].loc[0]*Utmp[0] + catm[j].loc[1]*Utmp[1] +
          catm[j].loc[2]*Utmp[2];
        dy = catm[j].loc[0]*Utmp[3] + catm[j].loc[1]*Utmp[4] +
          catm[j].loc[2]*Utmp[5];
        dz = catm[j].loc[0]*Utmp[6] + catm[j].loc[1]*Utmp[7] +
          catm[j].loc[2]*Utmp[8];
        catm[j].loc[0] = dx;
        catm[j].loc[1] = dy;
        catm[j].loc[2] = dz;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// R2Est: estimate the number of R2 interactions that will have to be computed
//        for a sector-to-sector interaction based on the sizes of subregions
//        being compared as well as the density and composition of either
//        sector.                                 
//                                                                      
// Arguments:
//   nr:      the number of atoms (in the first eight elements) and extra
//            points (in the second eight elements) in the cell's various
//            sectors.
//   nrEE:    the number of atoms in each sector with only electrostatic
//            properties
//   nrLJ:    the number of atoms in each sector with only Lennard-Jones
//            properties
//-----------------------------------------------------------------------------
static double R2Est(int* nr, int* nrEE, int* nrLJ, double* cdens, int sA,
		    int sB, double compfac, double regsize)
{
  double r2c;

  // The prefactor for subregion comparisons, attentuated 
  // by the failure rate of interactions due to mismatches
  // between atom properties.
  if (nr[sA] * nr[sB] > 1.0e-6) {
    r2c = compfac *
      (double)(nr[sA]*nr[sB] - nrLJ[sA]*nrEE[sB] - nrEE[sA]*nrLJ[sB]) /
      (double)(nr[sA]*nr[sB]);
  }
  else {
    r2c = 0.0;
  }

  // The interaction rate is now known.  All that must
  // be done is fold in the number of atoms in each of
  // the two regions that will actually participate.  
  r2c *= cdens[sA] * regsize * cdens[sB] * regsize;

  return r2c;
}

//-----------------------------------------------------------------------------
// CalcCellR2Cost: estimate the number of pairwise distances that must be
//                 computed for a particular cell.
//
// Arguments:
//   C:     the cell of interest                                        
//   CG:    the cell grid                                               
//-----------------------------------------------------------------------------
double CalcCellR2Cost(cell *C, cellgrid *CG)
{
  int i, j;
  int nrEE[8], nrLJ[8];
  double nR2;
  double cdens[8], thickness[3];
  atomc *catm;

  // Count the numbers of atoms in each region with  
  // electrostatic, Lennard-Jones, or both properties
  for (i = 0; i < 8; i++) {
    catm = C->map[i];
    nrLJ[i] = 0;
    nrEE[i] = 0;
    for (j = 0; j < C->nr[i]; j++) {
      if (fabs(catm[j].q) > 1.0e-8 && catm[j].lj < 0) {
        nrEE[i] += 1;
      }
      else if (fabs(catm[j].q) < 1.0e-8 && catm[j].lj >= 0) {
        nrLJ[i] += 1;
      }
    }
  }

  // Estimate the density of each sector
  thickness[0] = 0.5/CG->celldim[4];
  thickness[1] = 0.5/CG->celldim[5];
  thickness[2] = 0.5/CG->celldim[6];
  cdens[0] = thickness[0] * thickness[1] * thickness[2];
  cdens[1] = thickness[1] * thickness[2];
  cdens[2] = thickness[0] * thickness[2];
  cdens[3] = thickness[2];
  cdens[4] = thickness[0] * thickness[1];
  cdens[5] = thickness[1];
  cdens[6] = thickness[0];
  cdens[7] = 1.0;
  for (i = 0; i < 8; i++) {
    cdens[i] = C->nr[i] / cdens[i];
  }
  
  // Estimate the number of distances to test
  nR2 = C->nr[0] * (C->nr[0] - 1) / 2 - nrEE[0] * nrLJ[0];
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 0, 1, 0.625,
               thickness[1]*thickness[2]);
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 0, 2, 0.625,
               thickness[0]*thickness[2]);
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 0, 4, 0.625,
               thickness[0]*thickness[1]);
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 0, 3, 70.0/256.0, thickness[2]);
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 0, 5, 70.0/256.0, thickness[1]);
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 0, 6, 70.0/256.0, thickness[0]);
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 1, 4, 70.0/256.0, thickness[1]);
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 2, 4, 70.0/256.0, thickness[0]);
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 1, 2, 70.0/256.0, thickness[2]);
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 0, 7, 442.0/4096.0, 0.5*sqrt(3.0));
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 1, 6, 442.0/4096.0, 0.5*sqrt(3.0));
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 2, 5, 442.0/4096.0, 0.5*sqrt(3.0));
  nR2 += R2Est(C->nr, nrEE, nrLJ, cdens, 3, 4, 442.0/4096.0, 0.5*sqrt(3.0));
  
  return nR2;
}
