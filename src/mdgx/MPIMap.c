#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "Constants.h"
#include "mdgxVector.h"
#include "Matrix.h"
#include "MPIMap.h"
#include "Grid.h"
#include "Macros.h"
#include "BSpline.h"
#include "CrdManip.h"
#include "VirtualSites.h"

#include "TrajectoryDS.h"
#include "TopologyDS.h"
#include "pmeRecipDS.h"
#include "CellManipDS.h"

//-----------------------------------------------------------------------------
// SelectSystemsToTend: select the systems that this process will tend.  This
//                      routine takes the processor map identifying the various
//                      clusters and makes sense of out what systems to place
//                      on which processors.           
//
// Arguments:                                                           
//   tj:       the trajectory control information                       
//-----------------------------------------------------------------------------
void SelectSystemsToTend(trajcon *tj)
{
  int i, j, iclus, isys, icpu, sysCPU;
  lgrp* GroupMap;

  sysCPU = (tj->TI == 1) ? tj->nthreads : tj->nthreads / tj->nsys;
  tj->SystemCPUs = CreateImat(tj->nsys, MAX(sysCPU, 1));

  // Special cases: just one system, or TI engaged
  if (tj->nsys == 1 || tj->TI == 1) {
    for (i = 0; i < tj->nsys; i++) {
      for (j = 0; j < tj->nthreads; j++) {
	tj->SystemCPUs.map[i][j] = j;
      }
    }
  }

  // Case I: there are as many or more threads than processes
  else if (sysCPU >= 1) {

    // Detect if there are going to be spare threads
    if (sysCPU*tj->nsys < tj->nthreads) {
      printf("SelectSystemsToTend >> Warning.  Dividing %d systems among %d "
	     "threads will\nSelectSystemToTend >> leave %d threads idle.\n",
	     tj->nsys, tj->nthreads, tj->nthreads - sysCPU*tj->nsys);
    }

    // Make a map to remember which processors are occupied
    GroupMap = (lgrp*)malloc(tj->nCPUcluster*sizeof(lgrp));
    for (i = 0; i < tj->nCPUcluster; i++) {
      GroupMap[i].natom = tj->CPUcluster[i].natom;
      GroupMap[i].atoms = (int*)calloc(GroupMap[i].natom, sizeof(int));
    }

    // Start placing systems on clusters
    isys = 0;
    iclus = 0;
    while (isys < tj->nsys && iclus < tj->nCPUcluster) {

      // Loop over all processes in this cluster, assigning blocks
      // of sysCPU processors to each system so long as a block of
      // that size is available.                                  
      icpu = 0;
      for (i = 0; i < tj->CPUcluster[iclus].natom / sysCPU; i++) {
	for (j = 0; j < sysCPU; j++) {
	  tj->SystemCPUs.map[isys][j] = tj->CPUcluster[iclus].atoms[icpu];
	  GroupMap[iclus].atoms[icpu] = 1;
	  icpu++;
	}
	isys++;
	if (isys == tj->nsys) {
	  break;
	}
      }
      iclus++;
    }

    // Finish placing systems
    iclus = 0;
    i = 0;
    while (isys < tj->nsys) {
      icpu = 0;
      while (i < sysCPU && icpu < tj->CPUcluster[iclus].natom) {
	if (GroupMap[iclus].atoms[icpu] == 0) {
	  tj->SystemCPUs.map[isys][i] = tj->CPUcluster[iclus].atoms[icpu];
	  i++;
	}
	icpu++;
      }
      if (i == sysCPU) {
	isys++;
	i = 0;
      }
    }

    // Free allocated memory
    for (i = 0; i < tj->nCPUcluster; i++) {
      free(GroupMap[i].atoms);
    }
    free(GroupMap);
  }

  // Case II: allocate systems to each process and loop as needed
  else {
    iclus = 0;
    icpu = 0;
    for (isys = 0; isys < tj->nsys; isys++) {
      tj->SystemCPUs.map[isys][0] = tj->CPUcluster[iclus].atoms[icpu];
      icpu++;
      if (icpu == tj->CPUcluster[iclus].natom) {
	iclus++;
	icpu = 0;
      }
      if (iclus == tj->nCPUcluster) {
	iclus = 0;
      }
    }
  }

  // Enumerate the systems that this process must worry about
  tj->MySystemCount = 0;
  for (i = 0; i < tj->nsys; i++) {
    for (j = 0; j < tj->SystemCPUs.col; j++) {
      if (tj->SystemCPUs.map[i][j] == tj->tid) {
	tj->MySystemCount += 1;
	break;
      }
    }
  }
  tj->MySystemDomain = (int*)malloc(tj->MySystemCount*sizeof(int));
  isys = 0;
  for (i = 0; i < tj->nsys; i++) {
    for (j = 0; j < tj->SystemCPUs.col; j++) {
      if (tj->SystemCPUs.map[i][j] == tj->tid) {
	tj->MySystemDomain[isys] = i;
	isys++;
	break;
      }
    }
  }

#ifdef MPI
  // Create communicators for each system
  MPI_Group cpugrp, gworld;

  MPI_Comm_group(MPI_COMM_WORLD, &gworld);
  tj->SysComm = (MPI_Comm*)malloc(tj->nsys*sizeof(MPI_Comm));
  for (i = 0; i < tj->nsys; i++) {
    MPI_Group_incl(gworld, tj->SystemCPUs.col, tj->SystemCPUs.map[i], &cpugrp);
    MPI_Comm_create(MPI_COMM_WORLD, cpugrp, &tj->SysComm[i]);
  }
  MPI_Group_free(&gworld);
#endif
}

#ifdef MPI
//-----------------------------------------------------------------------------
// InitMeshGather: begin the process of gathering grid information onto the
//                 master process for the SPME reciprocal space convolution.
//                 All processes other than the master will post asynchronous
//                 sends, the master will post asynchronous receives for all
//                 sub-processes, and direct-space calculations will proceed
//                 until all messages have been received.
//
// Arguments:                                                           
//   rcinp:     reciprocal space information (where the mesh grids are) 
//   CG:        the cell grid (for MPI process maps)                    
//   req:       array of MPI requests                                   
//   scatter:   0 for gathering, 1 for distributing                     
//-----------------------------------------------------------------------------
int InitMeshGatherSPME(reccon *rcinp, cellgrid *CG, MPI_Request* req,
		       int scatter)
{
  int i, maxelem;

  // Bail out if there is only one thread
  if (CG->nthreads == 1) {
    return 0;
  }

  // Master process posts receives or sends
  maxelem = rcinp->ng[0]*rcinp->ng[1]*rcinp->ng[2];
  if (CG->tid == 0) {
    for (i = 1; i < CG->nthreads; i++) {
      if (scatter == 0) {
	MPI_Irecv(rcinp->QL[i].data, maxelem, MPI_DOUBLE, i, i, CG->dspcomm,
		  &req[i-1]);
      }
      else {
	SendMeshPart(rcinp, CG, &req[i-1], i);
      }
    }
    return CG->nthreads-1;
  }

  // All other processes post sends or receives
  else {
    if (scatter == 0) {
      SendMeshPart(rcinp, CG, &req[0], 0);
    }
    else {
      MPI_Irecv(rcinp->QL[1].data, maxelem, MPI_DOUBLE, 0, CG->tid,
		CG->dspcomm, &req[0]);
    }
    return 1;
  }
}

//-----------------------------------------------------------------------------
// MapProcessMeshFootprint: determine what territory on the reciprocal space
//                          mesh the direct space decomposition cells owned by
//                          a particular process project onto.
//
// Arguments:                                                           
//   CG:       the cell grid (identifies cells owned by this process)   
//   rcinp:    reciprocal space control information (dictates the mesh  
//             dimensions)                                              
//   crd:      the atom coordinates                                     
//-----------------------------------------------------------------------------
void MapProcessMeshFootprint(cellgrid *CG, reccon *rcinp, coord *crd)
{
  int h, i, j, k, ci, cj, ck, mini, minj, mink, maxi, maxj, maxk;
  int pagefull, colfull, colsum, onpc;
  int* coltmp;
  int** pagtmp;
  coord ptcrd;
  bmap pmmap;
  cell *C;
  ibook mygrid;

  // If this is the master process, post receives
  if (CG->tid == 0) {
    CG->MeshCommPlan = (gsplc*)malloc(CG->nthreads*sizeof(gsplc));
    int* ibuff;
    MPI_Request* req;
    MPI_Status* stt;
    ibuff = (int*)malloc(3*CG->nthreads*sizeof(int));
    req = (MPI_Request*)malloc(CG->nthreads*sizeof(MPI_Request));
    stt = (MPI_Status*)malloc(CG->nthreads*sizeof(MPI_Status));
    for (i = 1; i < CG->nthreads; i++) {
      MPI_Irecv(&ibuff[3*i], 3, MPI_INT, i, i, CG->dspcomm, &req[i-1]);
    }
    MPI_Waitall(CG->nthreads-1, req, stt);
    for (i = 1; i < CG->nthreads; i++) {
      CG->MeshCommPlan[i].npag = ibuff[3*i];
      CG->MeshCommPlan[i].ncol = ibuff[3*i+1];
      CG->MeshCommPlan[i].npc = ibuff[3*i+2];
      CG->MeshCommPlan[i].pagl = (int*)malloc((ibuff[3*i]+1)*sizeof(int));
      CG->MeshCommPlan[i].coll = (int*)malloc((2*ibuff[3*i+1]+1)*sizeof(int));
      CG->MeshCommPlan[i].pcl = (int*)malloc((4*ibuff[3*i+2]+1)*sizeof(int));
    }
    for (i = 1; i < CG->nthreads; i++) {
      MPI_Irecv(CG->MeshCommPlan[i].pagl, ibuff[3*i], MPI_INT, i, i,
		CG->dspcomm, &req[i-1]);
    }
    MPI_Waitall(CG->nthreads-1, req, stt);
    for (i = 1; i < CG->nthreads; i++) {
      MPI_Irecv(CG->MeshCommPlan[i].coll, 2*ibuff[3*i+1], MPI_INT, i, i,
		CG->dspcomm, &req[i-1]);
    }
    MPI_Waitall(CG->nthreads-1, req, stt);
    for (i = 1; i < CG->nthreads; i++) {
      MPI_Irecv(CG->MeshCommPlan[i].pcl, 4*ibuff[3*i+2], MPI_INT, i, i,
		CG->dspcomm, &req[i-1]);
    }
    MPI_Waitall(CG->nthreads-1, req, stt);

    // Free allocated memory
    free(req);
    free(stt);
    free(ibuff);

    return;
  }

  // Allocate space for fake coordinates and P->M maps
  ptcrd = CreateCoord(2);
  pmmap = CreateBmap(rcinp, 2);
  ptcrd.isortho = crd->isortho;
  mygrid = CreateIbook(rcinp->ng[0], rcinp->ng[1], rcinp->ng[2]);
  for (i = 0; i < 6; i++) {
    ptcrd.gdim[i] = crd->gdim[i];
  }
  CompXfrm(ptcrd.gdim, ptcrd.U, ptcrd.invU);

  // Loop over all cells in the grid, work
  // on the one owned by this process     
  for (h = 0; h < CG->MyCellCount; h++) {
    C = &CG->data[CG->MyCellDomain[h]];

    // Coordinates for "atoms" at the extrema of this cell
    for (i = 0; i < 3; i++) {
      ptcrd.loc[i] = C->orig[i];
      ptcrd.loc[i+3] = C->midp[i] - 1.0e-8;
    }
    SplCoeff(&ptcrd, &pmmap, rcinp);
    
    // Color the fake grid
    mini = pmmap.xcof[0].m;
    maxi = pmmap.xcof[2*rcinp->ordr[0]-1].m;
    minj = pmmap.ycof[0].m;
    maxj = pmmap.ycof[2*rcinp->ordr[1]-1].m;
    mink = pmmap.zcof[0].m;
    maxk = pmmap.zcof[2*rcinp->ordr[2]-1].m;
    if (mini > maxi) {
      mini -= rcinp->ng[0];
    }
    if (minj > maxj) {
      minj -= rcinp->ng[1];
    }
    if (mink > maxk) {
      mink -= rcinp->ng[2];
    }
    for (i = mini; i <= maxi; i++) {
      ci = (i < 0) ? i + rcinp->ng[0] : i;
      for (j = minj; j <= maxj; j++) {
	cj = (j < 0) ? j + rcinp->ng[1] : j;
	for (k = mink; k <= maxk; k++) {
	  ck = (k < 0) ? k + rcinp->ng[2] : k;
	  mygrid.map[ci][cj][ck] = 1;
	}
      }
    }
  }

  // Make the grid communication map
  gsplc gcomm;
  for (h = 0; h < 2; h++) {

    // First pass determines the number of each data
    // type, second pass identifies the data itself 
    if (h == 1) {
      gcomm.pagl = (int*)malloc(gcomm.npag*sizeof(int));
      gcomm.coll = (int*)malloc(2*gcomm.ncol*sizeof(int));
      gcomm.pcl = (int*)malloc(4*gcomm.npc*sizeof(int));
    }
    gcomm.npag = 0;
    gcomm.ncol = 0;
    gcomm.npc = 0;
    for (i = 0; i < mygrid.pag; i++) {

      // Determine whether this process controls a complete page
      pagefull = 1;
      pagtmp = mygrid.map[i];
      for (j = 0; j < mygrid.row; j++) {
	coltmp = pagtmp[j];
	for (k = 0; k < mygrid.col; k++) {
	  if (coltmp[k] == 0) {
	    pagefull = 0;
	    k = mygrid.col;
	    j = mygrid.row;
	  }
	}
      }
      if (pagefull == 1) {
	if (h == 1) {
	  gcomm.pagl[gcomm.npag] = i;
	}
	gcomm.npag += 1;
	continue;
      }

      // Determine whether this process controls full columns  
      // or just pieces, and count how many pieces there may be
      for (j = 0; j < mygrid.row; j++) {
	colfull = 1;
	colsum = 0;
	coltmp = pagtmp[j];
	for (k = 0; k < mygrid.col; k++) {
	  if (coltmp[k] == 0) {
	    colfull = 0;
	  }
	  else {
	    colsum += 1;
	  }
	}
	if (colfull == 1) {
	  if (h == 1) {
	    gcomm.coll[2*gcomm.ncol] = i;
	    gcomm.coll[2*gcomm.ncol+1] = j;
	  }
	  gcomm.ncol += 1;
	}
	else if (colsum > 0) {
	  onpc = 0;
	  for (k = 0; k < mygrid.col; k++) {
	    if (coltmp[k] != onpc) {
	      if (onpc == 0) {
		if (h == 1) {
		  gcomm.pcl[4*gcomm.npc] = i;
		  gcomm.pcl[4*gcomm.npc+1] = j;
		  gcomm.pcl[4*gcomm.npc+2] = k;
		}
	      }
	      else {
		if (h == 1) {
		  gcomm.pcl[4*gcomm.npc+3] = k;
		}
		gcomm.npc += 1;
	      }
	      onpc = coltmp[k];
	    }
	  }
	  if (onpc == 1) {
	    if (h == 1) {
	      gcomm.pcl[4*gcomm.npc+3] = mygrid.col;
	    }
	    gcomm.npc += 1;
	  }
	}
      }
    }
  }

  // Free allocated memory
  DestroyCoord(&ptcrd);
  DestroyBmap(&pmmap);
  DestroyIbook(&mygrid);

  // Send this plan to the master process
  int icount[3];
  MPI_Request myreq;
  MPI_Status mystt;
  icount[0] = gcomm.npag;
  icount[1] = gcomm.ncol;
  icount[2] = gcomm.npc;
  MPI_Isend(icount, 3, MPI_INT, 0, CG->tid, CG->dspcomm, &myreq);
  MPI_Wait(&myreq, &mystt);
  MPI_Isend(gcomm.pagl, gcomm.npag, MPI_INT, 0, CG->tid, CG->dspcomm, &myreq);
  MPI_Wait(&myreq, &mystt);
  MPI_Isend(gcomm.coll, 2*gcomm.ncol, MPI_INT, 0, CG->tid, CG->dspcomm,
	    &myreq);
  MPI_Wait(&myreq, &mystt);
  MPI_Isend(gcomm.pcl, 4*gcomm.npc, MPI_INT, 0, CG->tid, CG->dspcomm, &myreq);
  MPI_Wait(&myreq, &mystt);

  // Remember this plan for this process
  CG->MeshCommPlan = (gsplc*)malloc(sizeof(gsplc));
  CG->MeshCommPlan[0] = gcomm;
}

//-----------------------------------------------------------------------------
// SendMeshPart: package a part of the reciprocal space mesh for export to the
//               master process using information compiled during the creation
//               of the cell grid to determine the minimum amount of mesh data.
//
// Arguments:                                                           
//   rcinp:  the reciprocal space control information                   
//   CG:     the cell grid (contains mesh communication instructions)   
//   req:    array of MPI requests                                      
//   ncpu:   the process to which the mesh is being sent--0 if the      
//           meshes are being focused onto the master, > 0 if the master
//           is redistributing them to process ID "scatter"             
//-----------------------------------------------------------------------------
void SendMeshPart(reccon *rcinp, cellgrid *CG, MPI_Request *req, int ncpu)
{
  int i, j, k, npt;
  int ng[3];
  double *xdata, *coltmp;
  double **pagtmp;
  double ***qmesh;
  gsplc *gcomm;

  // Pointers and shortcuts
  gcomm = &CG->MeshCommPlan[ncpu];
  const int npag = gcomm->npag;
  const int ncol = gcomm->ncol;
  const int npc = gcomm->npc;
  ng[0] = rcinp->ng[0];
  ng[1] = rcinp->ng[1];
  ng[2] = rcinp->ng[2];
  xdata = (ncpu == 0) ? rcinp->QL[1].data : rcinp->QL[ncpu].data;
  qmesh = rcinp->QL[0].map;
  npt = 0;

  // Package the full pages
  for (i = 0; i < npag; i++) {
    pagtmp = qmesh[gcomm->pagl[i]];
    for (j = 0; j < ng[1]; j++) {
      coltmp = pagtmp[j];
      for (k = 0; k < ng[2]; k++) {
	xdata[npt] = coltmp[k];
	npt++;
      }
    }
  }

  // Package the full columns
  for (i = 0; i < ncol; i++) {
    coltmp = qmesh[gcomm->coll[2*i]][gcomm->coll[2*i+1]];
    for (j = 0; j < ng[2]; j++) {
      xdata[npt] = coltmp[j];
      npt++;
    }
  }

  // Package the column pieces
  for (i = 0; i < npc; i++) {
    coltmp = qmesh[gcomm->pcl[4*i]][gcomm->pcl[4*i+1]];
    for (j = gcomm->pcl[4*i+2]; j < gcomm->pcl[4*i+3]; j++) {
      xdata[npt] = coltmp[j];
      npt++;
    }
  }

  // Send the message via MPI
  if (ncpu == 0) {
    MPI_Isend(xdata, npt, MPI_DOUBLE, 0, CG->tid, CG->dspcomm, req);
  }
  else {
    MPI_Isend(xdata, npt, MPI_DOUBLE, ncpu, ncpu, CG->dspcomm, req);
  }
}

//-----------------------------------------------------------------------------
// RecvMeshPart: unpack parts of the reciprocal space mesh on the master
//               process using information compiled during the creation of the
//               cell grid to determine the minimum amount of mesh data.
//
// Arguments:                                                           
//   rcinp:  the reciprocal space control information                   
//   CG:     the cell grid (contains mesh communication instructions)   
//   ncpu:   the process from which this mesh data was received--0 if   
//           the master sent this (electrostatic potential) data, > 0 if
//           one of the slaves sent this (charge mesh) data             
//-----------------------------------------------------------------------------
void RecvMeshPart(reccon *rcinp, cellgrid *CG, int ncpu)
{
  int i, j, k, npt;
  int ng[3];
  double *xdata, *coltmp;
  double **pagtmp;
  double ***qmesh;
  gsplc *gcomm;

  // Pointers and shortcuts
  gcomm = &CG->MeshCommPlan[ncpu];
  const int npag = gcomm->npag;
  const int ncol = gcomm->ncol;
  const int npc = gcomm->npc;
  ng[0] = rcinp->ng[0];
  ng[1] = rcinp->ng[1];
  ng[2] = rcinp->ng[2];
  xdata = (ncpu == 0) ? rcinp->QL[1].data : rcinp->QL[ncpu].data;
  qmesh = rcinp->QL[0].map;
  npt = 0;

  // If electrostatic potential data is coming from
  // the master, unpack it directly into the mesh  
  if (ncpu == 0) {

    // Unpack the full pages
    for (i = 0; i < npag; i++) {
      pagtmp = qmesh[gcomm->pagl[i]];
      for (j = 0; j < ng[1]; j++) {
	coltmp = pagtmp[j];
	for (k = 0; k < ng[2]; k++) {
	  coltmp[k] = xdata[npt];
	  npt++;
	}
      }
    }

    // Unpack the full columns
    for (i = 0; i < ncol; i++) {
      coltmp = qmesh[gcomm->coll[2*i]][gcomm->coll[2*i+1]];
      for (j = 0; j < ng[2]; j++) {
	coltmp[j] = xdata[npt];
	npt++;
      }
    }

    // Unpack the column pieces
    for (i = 0; i < npc; i++) {
      coltmp = qmesh[gcomm->pcl[4*i]][gcomm->pcl[4*i+1]];
      for (j = gcomm->pcl[4*i+2]; j < gcomm->pcl[4*i+3]; j++) {
	coltmp[j] = xdata[npt];
	npt++;
      }
    }
  }
  else {

    // Unpack the full pages
    for (i = 0; i < gcomm->npag; i++) {
      pagtmp = qmesh[gcomm->pagl[i]];
      for (j = 0; j < ng[1]; j++) {
	coltmp = pagtmp[j];
	for (k = 0; k < ng[2]; k++) {
	  coltmp[k] += xdata[npt];
	  npt++;
	}
      }
    }

    // Unpack the full columns
    for (i = 0; i < ncol; i++) {
      coltmp = qmesh[gcomm->coll[2*i]][gcomm->coll[2*i+1]];
      for (j = 0; j < ng[2]; j++) {
	coltmp[j] += xdata[npt];
	npt++;
      }
    }

    // Unpack the column pieces
    for (i = 0; i < npc; i++) {
      coltmp = qmesh[gcomm->pcl[4*i]][gcomm->pcl[4*i+1]];
      for (j = gcomm->pcl[4*i+2]; j < gcomm->pcl[4*i+3]; j++) {
	coltmp[j] += xdata[npt];
	npt++;
      }
    }
  }
}
#endif

//-----------------------------------------------------------------------------
// MapProcessors: this function examines the host names of processors in the
//                list to determine how many separate nodes are involved in
//                the calculation.  Processes with the same host name are
//                assumed to be tightly connected.        
//                                                                      
// Arguments:                                                           
//   tj:      trajectory control information                            
//-----------------------------------------------------------------------------
void MapProcessors(trajcon *tj)
{
#ifdef MPI
  int i, j, k, ngrp;
  int* grpcount;
  int* grplist;
  char* thishost;
  cmat pnames;

  // Gather the names of all processes
  thishost = (char*)malloc(MAXTITL*sizeof(char));
  gethostname(thishost, MAXTITL);
  pnames = CreateCmat(tj->nthreads, MAXTITL);
  MPI_Gather(thishost, MAXTITL, MPI_CHAR, pnames.data, MAXTITL, MPI_CHAR,
	     0, MPI_COMM_WORLD);

  // Determine which processors are tightly connected
  grplist = (int*)malloc(tj->nthreads*sizeof(int));
  grpcount = (int*)malloc(tj->nthreads*sizeof(int));
  if (tj->tid == 0) {
    ngrp = 0;
    SetIVec(grplist, tj->nthreads, -1);
    i = 0;
    ngrp = 0;
    while (i < tj->nthreads) {
      if (grplist[i] == -1) {

	// Find all members of this group
	grplist[i] = ngrp;
	grpcount[ngrp] = 1;
	for (j = 1; j < tj->nthreads; j++) {
	  if (strcmp(pnames.map[j], pnames.map[i]) == 0) {
	    grplist[j] = ngrp;
	    grpcount[ngrp] += 1;
	  }
	}
	ngrp++;
      }
      i++;
    }
  }

  // Broadcast the group information
  MPI_Bcast(&ngrp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(grplist, tj->nthreads, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(grpcount, ngrp, MPI_INT, 0, MPI_COMM_WORLD);

  // Store the group information in the trajectory struct
  tj->nCPUcluster = ngrp;
  tj->CPUcluster = (lgrp*)malloc(ngrp*sizeof(lgrp));
  for (i = 0; i < ngrp; i++) {
    tj->CPUcluster[i].natom = grpcount[i];
    tj->CPUcluster[i].atoms = (int*)malloc(grpcount[i]*sizeof(int));
    k = 0;
    for (j = 0; j < tj->nthreads; j++) {
      if (grplist[j] == i) {
	tj->CPUcluster[i].atoms[k] = j;
	k++;
      }
    }
  }

  // Free allocated memory
  free(thishost);
  free(grpcount);
  free(grplist);
  DestroyCmat(&pnames);
#else

  // Dummy allocation of CPU clusters
  tj->nCPUcluster = 1;
  tj->CPUcluster = (lgrp*)malloc(sizeof(lgrp));
  tj->CPUcluster[0].natom = 1;
  tj->CPUcluster[0].atoms = (int*)malloc(sizeof(int));
  tj->CPUcluster[0].atoms[0] = 0;
#endif
}

//-----------------------------------------------------------------------------
// BalanceSmallSPME: initialize the load distribution for a small SPME  
//                   calculation.  The load balancing in this case is fairly
//                   straightforward: divide direct space and particle <-->
//                   mesh work amongs all processors and save one for computing
//                   FFTs.  Assign modest direct space and particle <--> mesh
//                   work to the processor tasked with crunching the FFT, and
//                   back it off dynamically as needed to optimize performance.
//
// Arguments:                                                           
//   CG:     the (newly formed) cell grid                               
//   tj:     trajectory control information                             
//-----------------------------------------------------------------------------
static void BalanceSmallSPME(cellgrid *CG)
{
  int i, j, k, npmin, nleft, ccpu, cproc, jdir, kdir;
  int* cellcount;

  // The master process gets the reciprocal space FFT work;    
  // estimate that in the typical case the FFT will take 1/10th
  // as long as the direct space nonbonded calculation.        
  const int ncell = CG->ncell;
  CG->MyCellDomain = (int*)malloc(ncell*sizeof(int));

  // First, determine how many cells each
  // process will have, at a minimum     
  nleft = ncell % CG->nthreads;
  npmin = (ncell - nleft) / CG->nthreads;
  cellcount = (int*)malloc(CG->nthreads*sizeof(int));
  SetIVec(cellcount, CG->nthreads, npmin);
  for (i = CG->nthreads-1; i >= CG->nthreads-nleft; i--) {
    cellcount[i] += 1;
  }

  // Assume that the reciprocal space convolution
  // will take 1/10th of the direct space work    
  nleft = ncell/10;
  nleft = MIN(nleft, cellcount[0]);
  if (CG->nthreads > 1) {
    cellcount[0] -= nleft;
    while (nleft > 0) {
      for (i = 1; i < CG->nthreads; i++) {
	if (nleft > 0) {
	  cellcount[i] += 1;
	  nleft--;
	}
      }
    }
  }

  // Snake through the cell grid, assigning
  // cells to CPUs in an orderly fashion.  
  ccpu = CG->nthreads-1;
  cproc = cellcount[ccpu];
  jdir = 1;
  kdir = 1;
  i = 0;
  j = 0;
  k = 0;
  while (ccpu >= 0 && cproc > 0) {

    // Assign this cell to the current CPU
    CG->map[i][j][k].CGRank = ccpu;
    cproc--;
    if (cproc == 0) {

      // Move on to a new processor
      ccpu--;
      if (ccpu >= 0) {
	cproc = cellcount[ccpu];
      }
    }

    // Increment the counters
    if ((kdir == 1 && k < CG->ng[2]-1) || (kdir == -1 && k > 0)) {
      k += kdir;
    }
    else {
      kdir *= -1;
      if ((jdir == 1 && j < CG->ng[1]-1) || (jdir == -1 && j > 0)) {
	j += jdir;
      }
      else {
	jdir *= -1;
	i++;
      }
    }
  }

  // List out all cells associated with each process,
  // and determine the master process's load.        
  j = 0;
  k = 0;
  for (i = 0; i < ncell; i++) {
    if (CG->data[i].CGRank == CG->tid) {
      CG->MyCellDomain[j] = i;
      j++;
    }
    if (CG->data[i].CGRank == 0) {
      k++;
    }
  }
  CG->MyCellCount = j;
  CG->MasterHalfLoad = k / 2;

  // Free allocated memory
  free(cellcount);
}

//-----------------------------------------------------------------------------
// CellGridEquipComm: the cellgrid struct is assigned its communicator, which
//                    as defined at the end of this library's function
//                    SelectSystemsToTend() is somewhat akin to MPI_COMM_WORLD
//                    for the system that the cellgrid struct describes.
//
// Arguments:                                                           
//   CG:     the (newly formed) cell grid                               
//   tj:     trajectory control information                             
//-----------------------------------------------------------------------------
void CellGridEquipComm(cellgrid *CG, trajcon *tj)
{
  // Store the number of threads working on this cell grid
  CG->nthreads = tj->SystemCPUs.col;

#ifdef MPI
  // Transfer the communicator
  MPI_Comm_dup(tj->SysComm[CG->sysID], &CG->dspcomm);
  MPI_Comm_rank(CG->dspcomm, &CG->tid);
#else
  CG->tid = 0;
#endif
}

//-----------------------------------------------------------------------------
// InitLoadBalance: initialize the load balancing for a new cell grid.  The
//                  plan is not necessarily to have either direct or reciprocal
//                  space processors "own" the others, but to set up an initial
//                  plan for distributing the workload in a balanced manner.
//                  Plans created by this routine will be modified later.
//
// Arguments:                                                           
//   CG:     the (newly formed) cell grid                               
//   rcinp:  the reciprocal space control information                   
//   tj:     trajectory control information                             
//-----------------------------------------------------------------------------
void InitLoadBalance(cellgrid *CG, reccon *rcinp)
{
  int i;

  // Easy if there's only one process
  if (CG->nthreads == 1) {
    CG->tid = 0;
    CG->nthreads = 1;
    for (i = 0; i < CG->ncell; i++) {
      CG->data[i].CGRank = 0;
    }
    CG->MyCellCount = CG->ncell;
    CG->MasterHalfLoad = CG->ncell / 2;
    CG->MyCellDomain = CountUp(CG->ncell);
    return;
  }

  // For a standard SPME calculation
  if (rcinp->nlev == 1) {
    if (CG->nthreads <= 8) {
      BalanceSmallSPME(CG);
    }
    else {
      printf("InitLoadBalance >> Error.  Not yet prepared for parallelism "
	     "beyond 8 CPUs.\n");
      exit(1);
    }
  }

  // FIX ME!!!  MLE deserves better!
  if (rcinp->nlev > 1) {
    if (CG->nthreads <= 8) {
      BalanceSmallSPME(CG);
    }
  }
}

//-----------------------------------------------------------------------------
// ConstructCellPulse: determine the communications between (and within)
//                     processes that must occur in order for cells to  
//                     communicate with their neighbors as part of an X, Y,
//                     or Z pulse.                                   
//                                                                      
// Arguments:                                                           
//   CG:     the (newly formed) cell grid                               
//   imove:  the direction of the pulse (0 = X, 1 = Y, 2 = Z)           
//   rev:    0 for a forward pulse (coordinate sharing), 1 for a reverse
//           pulse (force merger)                                       
//-----------------------------------------------------------------------------
static ashr ConstructCellPulse(cellgrid *CG, int imove, int rev)
{
  int i, j, sendcon, recvcon, scell, rcell;
  int dloc[3];
  int* procpush;
  int* procpull;
  int* procself;
  int* cellpush;
  int* otgsend;
  int* increcv;
  ashr cshr;
  cell *C, *D;

  // Scratch space
  procpush = (int*)malloc(CG->ncell*sizeof(int));
  procpull = (int*)malloc(CG->ncell*sizeof(int));
  procself = (int*)malloc(CG->ncell*sizeof(int));
  cellpush = (int*)malloc(CG->ncell*sizeof(int));
  increcv = (int*)malloc(CG->nthreads*sizeof(int));
  otgsend = (int*)malloc(CG->nthreads*sizeof(int));

  // Label the cells as they push
  // or pull from one another    
  for (i = 0; i < CG->ncell; i++) {

    // Determine the origin (C) and destination (D) cells
    C = &CG->data[i];
    for (j = 0; j < 3; j++) {
      dloc[j] = C->gbin[j];
    }
    if (rev == 0) {
      dloc[imove] = (C->gbin[imove] == 0) ? CG->ng[imove]-1 : C->gbin[imove]-1;
    }
    else {
      dloc[imove] = (C->gbin[imove] == CG->ng[imove]-1) ? 0 : C->gbin[imove]+1;
    }
    D = &CG->map[dloc[0]][dloc[1]][dloc[2]];

    // Make notes of the pushing and pulling
    procself[i] = C->CGRank;
    procpush[i] = D->CGRank;
    cellpush[i] = D->gbin[3];
    procpull[D->gbin[3]] = C->CGRank;
  }

  // Count all of the sends and receives
  for (i = 0; i < CG->nthreads; i++) {
    otgsend[i] = 0;
    increcv[i] = 0;
  }
  for (i = 0; i < CG->ncell; i++) {
    if (procself[i] == CG->tid) {
      otgsend[procpush[i]] += 1;
      increcv[procpull[i]] += 1;
    }
  }
  sendcon = 0;
  recvcon = 0;
  scell = 0;
  rcell = 0;
  for (i = 0; i < CG->nthreads; i++) {
    if (otgsend[i] > 0) {
      sendcon++;
      scell += otgsend[i];
    }
    if (increcv[i] > 0) {
      if (i != CG->tid) {
	recvcon++;
      }
      rcell += increcv[i];
    }
  }

  // Allocate memory for sends and receives
  cshr.nsend = sendcon;
  cshr.nrecv = recvcon;
  cshr.send = (pcgrp*)malloc(MAX(sendcon, 1)*sizeof(pcgrp));
  cshr.recv = (pcgrp*)malloc(MAX(recvcon, 1)*sizeof(pcgrp));
  cshr.slist = (int*)malloc(MAX(scell, 1)*sizeof(int));
  cshr.rlist = (int*)malloc(MAX(rcell, 1)*sizeof(int));

  // Detail all sends and receives between processes
  sendcon = 0;
  recvcon = 0;
  scell = 0;
  rcell = 0;
  for (i = 0; i < CG->nthreads; i++) {
    if (otgsend[i] != 0 && i != CG->tid) {
      cshr.send[sendcon].ncell = otgsend[i];
      cshr.send[sendcon].partner = i;
      cshr.send[sendcon].BaseID =
	(6*(CG->tid*CG->nthreads + cshr.send[sendcon].partner) + 
	 (3*rev + imove))*DSP_MSG_TYPES;
      cshr.send[sendcon].cellpt = &cshr.slist[scell];
      for (j = 0; j < CG->ncell; j++) {
	if (procself[j] == CG->tid && procpush[j] == i) {
	  cshr.slist[scell] = j;
	  scell++;
	}
      }
      sendcon++;
    }
    if (increcv[i] != 0 && i != CG->tid) {
      cshr.recv[recvcon].ncell = increcv[i];
      cshr.recv[recvcon].partner = i;
      cshr.recv[recvcon].BaseID =
	(6*(i*CG->nthreads + CG->tid) + (3*rev + imove))*DSP_MSG_TYPES;
      cshr.recv[recvcon].cellpt = &cshr.rlist[rcell];
      for (j = 0; j < CG->ncell; j++) {
	if (procself[j] == CG->tid && procpull[j] == i) {
	  cshr.rlist[rcell] = j;
	  rcell++;
	}
      }
      recvcon++;
    }
  }

  // Detail the "self send" of each process
  // (which never actually goes over MPI)  
  if (otgsend[CG->tid] != 0) {
    cshr.send[sendcon].ncell = otgsend[CG->tid];
    cshr.send[sendcon].partner = CG->tid;
    cshr.send[sendcon].cellpt = &cshr.slist[scell];
    cshr.selfrecv.ncell = otgsend[CG->tid];
    cshr.selfrecv.partner = CG->tid;
    cshr.selfrecv.cellpt = &cshr.rlist[rcell];
    for (j = 0; j < CG->ncell; j++) {
      if (procself[j] == CG->tid) {
	if (procpush[j] == CG->tid) {
	  cshr.slist[scell] = j;
	  cshr.rlist[rcell] = cellpush[j];
	  scell++;
	  rcell++;
	}
      }
    }
  }
  else {
    cshr.selfrecv.ncell = 0;
  }

  // Free allocated memory
  free(procpush);
  free(procpull);
  free(procself);
  free(cellpush);
  free(increcv);
  free(otgsend);

  return cshr;
}

//-----------------------------------------------------------------------------
// MapProcessCellSharing: analyze the direct space communications needed for
//                        coordinate and force communication, and prepare plans
//                        for performing these operations during dynamics.
//
// Arguments:                                                           
//   CG:     the (newly formed) cell grid                               
//   tj:     trajectory control information                             
//-----------------------------------------------------------------------------
void MapProcessCellSharing(cellgrid *CG)
{
  int imove;

  // List communications for ShareCoordinates()
  for (imove = 0; imove < 3; imove++) {
    CG->DirCommPlan.mvshr[imove] = ConstructCellPulse(CG, imove, 0);
  }

  // List communications for MergeCellForces()
  for (imove = 0; imove < 3; imove++) {
    CG->DirCommPlan.frcmg[imove] = ConstructCellPulse(CG, imove, 1);
  }
}

//-----------------------------------------------------------------------------
// DestroyAshr: destroy a cell sharing plan.                            
//-----------------------------------------------------------------------------
void DestroyAshr(ashr *A)
{
  free(A->slist);
  free(A->rlist);
  free(A->send);
  free(A->recv);
}

//-----------------------------------------------------------------------------
// AllocatePooledBuffers: allocate space for buffers that will send and receive
//                        pooled messages between processes for a cell grid.
//
// Arguments:                                                           
//   CG:      the (recently created) cell grid                          
//-----------------------------------------------------------------------------
void AllocatePooledBuffers(cellgrid *CG)
{
  int i, j;
  int nsend, nrecv;
  ashr *cshr;

  // Find the largest number of sends and receives
  nsend = 0;
  nrecv = 0;
  for (i = 0; i < 3; i++) {
    cshr = &CG->DirCommPlan.mvshr[i];
    nsend = MAX(nsend, cshr->nsend);
    nrecv = MAX(nrecv, cshr->nrecv);
    cshr = &CG->DirCommPlan.frcmg[i];
    nsend = MAX(nsend, cshr->nsend);
    nrecv = MAX(nrecv, cshr->nrecv);
  }
  CG->nsend = nsend;
  CG->nrecv = nrecv;

  // Find the biggest import and export
  CG->maxexp = (int*)calloc(nsend+1, sizeof(int));
  CG->maximp = (int*)calloc(nrecv+1, sizeof(int));
  for (i = 0; i < 3; i++) {
    cshr = &CG->DirCommPlan.mvshr[i];
    for (j = 0; j < cshr->nsend; j++) {
      CG->maxexp[j] = MAX(CG->maxexp[j], cshr->send[j].ncell);
    }
    for (j = 0; j < cshr->nrecv; j++) {
      CG->maximp[j] = MAX(CG->maximp[j], cshr->recv[j].ncell);
    }
    cshr = &CG->DirCommPlan.frcmg[i];
    for (j = 0; j < cshr->nsend; j++) {
      CG->maxexp[j] = MAX(CG->maxexp[j], cshr->send[j].ncell);
    }
    for (j = 0; j < cshr->nrecv; j++) {
      CG->maximp[j] = MAX(CG->maximp[j], cshr->recv[j].ncell);
    }
  }
  for (i = 0; i < nsend; i++) {
    CG->maxexp[i] *= 4*CG->maxatom;
  }
  for (i = 0; i < nrecv; i++) {
    CG->maximp[i] *= 4*CG->maxatom;
  }

  // Allocate for the biggest import and export
  CG->import = (atomb**)malloc((nrecv+1)*sizeof(atomb*));
  CG->pexport = (atomb**)malloc((nsend+1)*sizeof(atomb*));
  CG->Vimport = (atombv**)malloc((nrecv+1)*sizeof(atombv*));
  CG->Vexport = (atombv**)malloc((nsend+1)*sizeof(atombv*));
  CG->Ximport = (atombx**)malloc((nrecv+1)*sizeof(atombx*));
  CG->Xexport = (atombx**)malloc((nsend+1)*sizeof(atombx*));
  CG->nexp = (int*)malloc(nsend*sizeof(int));
  for (i = 0; i < nsend; i++) {
    CG->pexport[i] = (atomb*)malloc(CG->maxexp[i]*sizeof(atomb));
    CG->Vexport[i] = (atombv*)malloc(CG->maxexp[i]*sizeof(atombv));
    CG->Xexport[i] = (atombx*)malloc(CG->maxexp[i]*sizeof(atombx));
  }
  for (i = 0; i < nrecv; i++) {
    CG->import[i] = (atomb*)malloc(CG->maximp[i]*sizeof(atomb));
    CG->Vimport[i] = (atombv*)malloc(CG->maximp[i]*sizeof(atombv));
    CG->Ximport[i] = (atombx*)malloc(CG->maximp[i]*sizeof(atombx));
  }
}

//-----------------------------------------------------------------------------
// PlanCoordRedux: plan the communication needed to inform the master process
//                 with information on the coordinates of all atoms.  This is
//                 performed with different buffers than cell-to-cell
//                 communication due to the uncertainty of whether either the
//                 import or export buffers would be large enough to
//                 accommodate each process's message.  
//
// Arguments:                                                           
//   CG:    the cell grid (for atom assignment information)             
//-----------------------------------------------------------------------------
void PlanCoordRedux(cellgrid *CG)
{
  int i, j, npcell;

  // Bail right out if there is only one thread
  if (CG->nthreads == 1) {
    return;
  }

  // Allocate buffers
  if (CG->tid == 0) {
    CG->CrdPoolSize = (int*)malloc((CG->nthreads - 1)*sizeof(int));
    CG->CrdPool = (atombx**)malloc((CG->nthreads - 1)*sizeof(atombx*));
    for (i = 0; i < CG->nthreads - 1; i++) {

      // Count the number of cells associated with this process
      npcell = 0;
      for (j = 0; j < CG->ncell; j++) {
	if (CG->data[j].CGRank == i+1) {
	  npcell++;
	}
      }
      CG->CrdPoolSize[i] = npcell*CG->maxatom;
      CG->CrdPool[i] = (atombx*)malloc(npcell*CG->maxatom*sizeof(atombx));
    }
  }
  else {
    CG->CrdPoolSize = (int*)malloc(sizeof(int));
    CG->CrdPoolSize[0] = CG->MyCellCount*CG->maxatom;
    CG->CrdPool = (atombx**)malloc(sizeof(atombx*));
    CG->CrdPool[0] = (atombx*)malloc(CG->CrdPoolSize[0]*sizeof(atombx));
  }
}

#ifdef MPI
//-----------------------------------------------------------------------------
// CoordinateReduction: reduce the coordinates to the master process of this
//                      cell grid.                                 
//
// Arguments:                                                           
//   CG:    the cell grid (for atom assignment information)             
//   crd:   coordinates (where velocities and other information is      
//          ultimately kept)                                            
//-----------------------------------------------------------------------------
void CoordinateReduction(cellgrid *CG, coord *crd, prmtop *tp)
{
  int i, j, k, atmid;
  int* AtomsOwned;
  int* cellcat;
  int* mycellcat;
  double* scratch;
  cell *C;
  atomc *tatm;

  // Return immediately if there is only one thread
  if (CG->nthreads == 1) {
    return;
  }

  // Map forces, velocities, and positions to the coordinate arrays,
  // but zero out coordinates for anything that this process doesn't
  // own, in preparation for MPI_Reduce.                            
  AtomsOwned = (int*)calloc(crd->natom, sizeof(int));
  cellcat = (int*)calloc(crd->natom, sizeof(int));
  mycellcat = (int*)calloc(crd->natom, sizeof(int));
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < C->nr[0]; j++) {
      atmid = C->data[j].id;
      AtomsOwned[atmid] = 1;
      mycellcat[atmid] = CG->MyCellDomain[i];
      for (k = 0; k < 3; k++) {
	crd->frc[3*atmid+k] = C->data[j].frc[k];
      }
    }
  }
  for (i = 0; i < crd->natom; i++) {
    if (AtomsOwned[i] == 0) {
      for (j = 3*i; j < 3*i+3; j++) {
	crd->vel[j] = 0.0;
	crd->frc[j] = 0.0;
	crd->loc[j] = 0.0;
      }
    }
  }

  // MPI_Reduce to get atom information to the master process
  MPI_Reduce(mycellcat, cellcat, crd->natom, MPI_INT, MPI_SUM, 0, CG->dspcomm);
  scratch = (double*)calloc(3*crd->natom, sizeof(double));
  MPI_Reduce(crd->loc, scratch, 3*crd->natom, MPI_DOUBLE, MPI_SUM, 0,
	     CG->dspcomm);
  if (CG->tid == 0) {
    ReflectDVec(crd->loc, scratch, 3*crd->natom);
  }
  SetDVec(scratch, 3*crd->natom, 0.0);
  MPI_Reduce(crd->vel, scratch, 3*crd->natom, MPI_DOUBLE, MPI_SUM, 0,
	     CG->dspcomm);
  if (CG->tid == 0) {
    ReflectDVec(crd->vel, scratch, 3*crd->natom);
  }
  SetDVec(scratch, 3*crd->natom, 0.0);
  MPI_Reduce(crd->frc, scratch, 3*crd->natom, MPI_DOUBLE, MPI_SUM, 0,
	     CG->dspcomm);
  if (CG->tid == 0) {
    ReflectDVec(crd->frc, scratch, 3*crd->natom);
    for (i = 0; i < CG->ncell; i++) {
      C = &CG->data[i];
      if (C->CGRank != 0) {
	C->nr[0] = 0;
      }
    }
    for (i = 0; i < crd->natom; i++) {
      if (AtomsOwned[i] == 1) {
	continue;
      }
      C = &CG->data[cellcat[i]];
      tatm = &C->data[C->nr[0]];
      C->nr[0] += 1;
      tatm->id = i;
      for (j = 0; j < 3; j++) {
	tatm->frc[j] = crd->frc[3*i+j];
      }
    }

    // If there are extra points, their positions are NOT
    // guaranteed to be current.  However, we do have all
    // of the locations of all atoms in a convenient XYZ 
    // list, so the extra points may be replaced now.    
    double *aloc, *bloc, *cloc, *dloc, *eploc;
    expt *tmr;
    for (i = 0; i < tp->nxtrapt; i++) {
      tmr = &tp->xtrapts[i];
      aloc = &crd->loc[3*tmr->fr1];
      bloc = &crd->loc[3*tmr->fr2];
      if (tmr->frstyle > 1) {
	cloc = &crd->loc[3*tmr->fr3];
      }
      if (tmr->frstyle == 6 || tmr->frstyle == 9 || tmr->frstyle == 10) {
	dloc = &crd->loc[3*tmr->fr4];
      }
      eploc = &crd->loc[3*tmr->atomid];
      XptLocator(aloc, aloc, bloc, cloc, dloc, eploc, eploc, tmr);
    }
  }

  // Free allocated memory
  free(scratch);
  free(cellcat);
  free(mycellcat);
  free(AtomsOwned);
}

//-----------------------------------------------------------------------------
// P2PCellRebalance: this routine manages the work of transferring cells
//                   between two processes to even out the workload.  It
//                   examines wait times experienced by each process in 
//                   the direct-space pool and decides which processes  
//                   will give and receive work units (cells).          
//-----------------------------------------------------------------------------


#endif
