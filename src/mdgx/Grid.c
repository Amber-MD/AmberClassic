#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Matrix.h"
#include "Grid.h"
#include "CrdManip.h"
#include "Constants.h"

//-----------------------------------------------------------------------------
// CreateIbook: create an integer book of dimensions M x N x P.         
//-----------------------------------------------------------------------------
ibook CreateIbook(int M, int N, int P)
{
  int i, j;
  ibook A;

  // Allocate memory
  A.data = (int*)calloc(M*N*P, sizeof(int));

  // Make the map
  A.map = (int***)malloc(M*sizeof(int**));
  for (i = 0; i < M; i++) {
    A.map[i] = (int**)malloc(N*sizeof(int*));
    for (j = 0; j < N; j++) {
      A.map[i][j] = &A.data[i*N*P + j*P];
    }
  }

  // Store the dimensionality
  A.pag = M;
  A.row = N;
  A.col = P;

  return A;
}

//-----------------------------------------------------------------------------
// CreateDbook: create a double-precision real book with M x N x P indices. 
//              There is added functionality for preparing the array for 1, 2,
//              and 3-dimensional real-to-complex (R2C) Discrete Fourier
//              Transforms (DFTs).                     
//                                                                      
// Arguments:                                                           
//   M:       the number of pages (the most slowly varying index)       
//   N:       the number of rows (the second most slowly varying index) 
//   P:       the number of columns (the most rapidly varying index)    
//   prepFFT: flag to request preparation for R2C transforms--setting prepFFT=1
//            pads along the 3rd dimension.  This will prepare the array for a
//            single 3D R2C DFT, as well as a series of 2D or 1D R2C DFTs.
//-----------------------------------------------------------------------------
dbook CreateDbook(int M, int N, int P, int prepFFT)
{
  int i, j, Pp;
  dbook A;

  // Prepare this book for 3D FFTs
  A.pfft = prepFFT;
  Pp = (prepFFT == 1) ? 2*(P/2+1) : P;

  // Allocate memory
  A.data = (double*)calloc(M*N*Pp, sizeof(double));

  // This is pointer-array abuse, but I know of no other way to do it
  if (prepFFT > 0) {
    A.fdata = (fftw_complex*)A.data;
  }

  // Make the map
  A.map = (double***)malloc(M*sizeof(double**));
  for (i = 0; i < M; i++) {
    A.map[i] = (double**)malloc(N*sizeof(double*));
    for (j = 0; j < N; j++) {
      A.map[i][j] = &A.data[(i*N+j)*Pp];
    }
  }

  // Set up the complex maps
  if (prepFFT == 1) {
    A.fmap = (fftw_complex***)malloc(M*sizeof(fftw_complex**));
    for (i = 0; i < M; i++) {
      A.fmap[i] = (fftw_complex**)malloc(N*sizeof(fftw_complex*));
      for (j = 0; j < N; j++) {
	A.fmap[i][j] = &A.fdata[(i*N+j)*Pp/2];
      }
    }
  }

  // Store the dimensionality
  A.pag = M;
  A.row = N;
  A.col = P;

  return A;
}

//-----------------------------------------------------------------------------
// CreateFbook: create a single-precision real book with M x N x P indices.
//              This is the format for grids read from files, or grids that
//              represent potential fields for which there is not a ton of
//              memory.
//                                                                      
// Arguments:                                                           
//   M:       the number of pages (the most slowly varying index)       
//   N:       the number of rows (the second most slowly varying index) 
//   P:       the number of columns (the most rapidly varying index)    
//   voxel:   matrix describing the side lengths of a grid bin / voxel  
//   orig:    the grid origin (3-element vector)                        
//-----------------------------------------------------------------------------
fbook CreateFbook(int M, int N, int P, dmat *voxel, double* orig)
{
  int i, j;
  fbook A;

  // Allocate memory
  A.data = (float*)calloc(M*N*P, sizeof(float));

  // Make the map
  A.map = (float***)malloc(M*sizeof(float**));
  for (i = 0; i < M; i++) {
    A.map[i] = (float**)malloc(N*sizeof(float*));
    for (j = 0; j < N; j++) {
      A.map[i][j] = &A.data[(i*N+j)*P];
    }
  }

  // Store the dimensionality
  A.pag = M;
  A.row = N;
  A.col = P;

  // Allocate the spacing
  A.L = CreateDmat(3, 3, 0);
  A.invL = CreateDmat(3, 3, 0);
  A.U = CreateDmat(3, 3, 0);
  A.invU = CreateDmat(3, 3, 0);

  // Compute the transformation matrices
  CopyDmat(&A.L, voxel, 1);
  ttInv(A.L, A.invL);
  for (i = 0; i < 3; i++) {
    A.U.map[0][i] = A.L.map[0][i] / M;
    A.U.map[1][i] = A.L.map[1][i] / N;
    A.U.map[2][i] = A.L.map[2][i] / P;
    A.orig[i] = orig[i];
  }
  ttInv(A.U, A.invU);

  // Test for orthogonality of the grid vectors
  if (fabs(A.invU.data[1]) < 1.0e-8 && fabs(A.invU.data[2]) < 1.0e-8 &&
      fabs(A.invU.data[5]) < 1.0e-8) {
    A.isortho = 1;
  }
  else {
    A.isortho = 0;
  }

  return A;
}

//-----------------------------------------------------------------------------
// PromoteFbook: promote an fbook grid to higher resolution and double  
//               precision.                                             
//                                                                      
// Arguments:                                                           
//   A:         the original grid                                       
//   pfac:      promotion factor, a value of 1 simply promotes precision
//   prepfft:   prepare for FFT                                         
//-----------------------------------------------------------------------------
dbook PromoteFbook(fbook *A, int pfac, int prepfft)
{
  int i, j, k;
  double invpfac;
  double crd[3], gvec[3];
  dbook D;

  D = CreateDbook(pfac*A->pag, pfac*A->row, pfac*A->col, prepfft);  
  invpfac = 1.0/pfac;
  for (i = 0; i < D.pag; i++) {
    gvec[0] = i * invpfac;
    for (j = 0; j < D.row; j++) {
      gvec[1] = j * invpfac;
      for (k = 0; k < D.col; k++) {
	gvec[2] = k * invpfac;
	DMatVecMult(&A->invL, gvec, crd); 
	TriInterp(crd, &D.map[i][j][k], A, 1, 3);
      }
    }
  }

  return D;
}

//-----------------------------------------------------------------------------
// CreateCbook: create a character (or unsigned 8-bit integer) book of  
//              dimensions M x N x P.                                   
//-----------------------------------------------------------------------------
cbook CreateCbook(int M, int N, int P)
{
  int i, j;
  cbook A;

  // Allocate memory
  A.data = (char*)calloc(M*N*P, sizeof(char));

  // Make the map
  A.map = (char***)malloc(M*sizeof(char**));
  for (i = 0; i < M; i++) {
    A.map[i] = (char**)malloc(N*sizeof(char*));
    for (j = 0; j < N; j++) {
      A.map[i][j] = &A.data[i*N*P + j*P];
    }
  }

  // Store the dimensionality
  A.pag = M;
  A.row = N;
  A.col = P;

  return A;
}

//-----------------------------------------------------------------------------
// DestroyIbook: destroy an integer book.                               
//-----------------------------------------------------------------------------
void DestroyIbook(ibook *A)
{
  int i;

  // Free the map
  for (i = 0; i < A->pag; i++) {
    free(A->map[i]);
  }
  free(A->map);

  // Free the data
  free(A->data);
}

//-----------------------------------------------------------------------------
// DestroyFbook: destroy a single-precision real book.                  
//-----------------------------------------------------------------------------
void DestroyFbook(fbook *A)
{
  int i;

  // Free the map
  for (i = 0; i < A->pag; i++) {
    free(A->map[i]);
  }
  free(A->map);

  // Free the grid spacing matrix
  DestroyDmat(&A->L);
  DestroyDmat(&A->invL);
  DestroyDmat(&A->U);
  DestroyDmat(&A->invU);

  // Free the data
  free(A->data);
}

//-----------------------------------------------------------------------------
// DestroyDbook: destroy a double-precision real book.                  
//-----------------------------------------------------------------------------
void DestroyDbook(dbook *A)
{
  int i;

  // Free the map
  for (i = 0; i < A->pag; i++) {
    free(A->map[i]);
  }
  free(A->map);
  if (A->pfft == 1) {
    for (i = 0; i < A->pag; i++) {
      free(A->fmap[i]);
    }
    free(A->fmap);
  }

  // Free the data
  free(A->data);
}

//-----------------------------------------------------------------------------
// DestroyCbook: destroy a character book.                              
//-----------------------------------------------------------------------------
void DestroyCbook(cbook *A)
{
  int i;

  // Free the map
  for (i = 0; i < A->pag; i++) {
    free(A->map[i]);
  }
  free(A->map);

  // Free the data
  free(A->data);
}

//-----------------------------------------------------------------------------
// CompressQ: convert a 64-bit double-precision real book to a 32-bit integer
//            book, with a list of exceptions needing special treatment.  The
//            integers range from roughly -2.1e9 to +2.1e9, so charges are
//            stored as integers in the range -2.1 to +2.1, with nine decimal
//            places of precision.  The vast majority grid points will fall in
//            this category.  A separate list enumerates the few exceptional
//            cases.       
//                                                                     
// Arguments:                                                           
//   cQ: the compressed charge book                                     
//   Q:  the uncompressed charge book                                   
//                                                                      
// Note that the compressed charge book is assumed to be pre-allocated, 
// as this routine is expected to be used many times to convert the same
// charge book into a compressed version.  See also DecompressQ in this 
// library.                                                             
//-----------------------------------------------------------------------------
void CompressQ(qbook *cQ, dbook *Q)
{
  int i, j, M, N, P;
  int *Qcomp;
  double *Qorig;

  // Set the exception counter to zero
  cQ->nexcp = 0;

  // Set pointers and abbreviate constants
  M = cQ->qdata.pag;
  N = cQ->qdata.row;
  P = cQ->qdata.col;
  Qcomp = cQ->qdata.data;
  Qorig = Q->data;

  // Checkbook dimensions
  if (M != Q->pag || N != Q->row || P != Q->col) {
    printf("CompressQ >> Error.  Book dimensions are inconsistent.\n"
	   "CompressQ >> Original book   = [ %4d pages of %4d x %4d ]\n"
	   "CompressQ >> Compressed book = [ %4d pages of %4d x %4d ]\n",
	   M, N, P, Q->pag, Q->row, Q->col);
    exit(1);
  }

  // Loop over all data
  for (i = 0; i < M*N*P; i++) {
    if (Qorig[i] >= -2.14748364 && Qorig[i] <= 2.14748364) {
      Qcomp[i] = Qorig[i]*1.0e9;
    }
    else {
      Qcomp[i] = 0.0;
      j = cQ->nexcp;
      cQ->eidx[j] = i;
      cQ->eval[j] = Qorig[i];
      cQ->nexcp = ++j;
      if (j > cQ->maxexcp) {
	cQ->maxexcp += 1024;
	cQ->eidx = (int*)realloc(cQ->eidx, cQ->maxexcp*sizeof(int));
	cQ->eval = (double*)realloc(cQ->eval, cQ->maxexcp*sizeof(double));
      }
    }
  }
}

//-----------------------------------------------------------------------------
// DecompressQ: convert a 32-bit integer book back into a 64-bit
//              double-precision real book, abiding a list of exceptions
//              showing grid points that contain lots of charge.        
//                                                                      
// Arguments:                                                           
//   Q:  the uncompressed charge book                                   
//   cQ: the compressed charge book                                     
//                                                                      
// As in CompressQ, the uncompressed charge book must be pre-allocated; this
// routine is designed to be used many times with the same data structures.
//-----------------------------------------------------------------------------
void DecompressQ(dbook *Q, qbook *cQ)
{
  int i, M, N, P;
  int *Qcomp;
  double *Qnew;

  // Set pointers and abbreviate constants
  M = cQ->qdata.row;
  N = cQ->qdata.col;
  P = cQ->qdata.pag;
  Qcomp = cQ->qdata.data;
  Qnew = Q->data;

  // Loop over all data
  for (i = 0; i < M*N*P; i++) {
    Qnew[i] = Qcomp[i];
  }

  // Loop over exceptions
  for (i = 0; i < cQ->nexcp; i++) {
    Qnew[cQ->eidx[i]] = cQ->eval[i];
  }
}

//-----------------------------------------------------------------------------
// ExtractDpage: extract a page from a book and copy it into a dmat struct.
//                                                                      
// Arguments:                                                           
//   B:   the book of double-precision real numbers                     
//   P:   the matrix of double-precision real numbers                   
//   n:   the page number to extract                                    
//   ll:  flag to tell whether the page is pre-allocated (1) or must be 
//          allocated (0)                                               
//-----------------------------------------------------------------------------
void ExtractDpage(dbook *B, dmat *P, int n, int ll)
{
  int i, j;
  double *btmp, *ptmp;

  // Allocate the matrix (page) if necessary.  Note that the page will
  // be not be prepared for FFTs.
  if (ll == 0) {
    *P = CreateDmat(B->row, B->col, 0);
  }

  // Check
  if (B->row > P->row) {
    printf("ExtractDpage >> Error.  %d rows in the book page, %d rows "
	   "available in the matrix.\n", B->row, P->row);
    exit(1);
  }
  if (B->col > P->col) {
    printf("ExtractDpage >> Error.  %d columns in the book page, %d columns "
	   "available in the matrix.\n", B->col, P->col);
    exit(1);
  }

  // Copy values
  for (i = 0; i < B->row; i++) {
    btmp = B->map[n][i];
    ptmp = P->map[i];
    for (j = 0; j < B->col; j++) {
      ptmp[j] = btmp[j];
    }
  }
}

//-----------------------------------------------------------------------------
// ScribeDpage: copy a double-precision matrix into a dbook struct.     
//                                                                      
// Arguments:                                                           
//   P:   the matrix of double-precision real numbers                   
//   B:   the book of double-precision real numbers                     
//   n:   the page number to which the matrix values should be written  
//-----------------------------------------------------------------------------
void ScribeDpage(dmat *P, dbook *B, int n)
{
  int i, j;
  double *btmp, *ptmp;

  // Check
  if (P->row > B->row) {
    printf("ScribeDpage >> Error.  %d rows in the book page, %d rows "
           "available in the matrix.\n", P->row, B->row);
    exit(1);
  }
  if (P->col > B->col) {
    printf("ScribeDpage >> Error.  %d columns in the book page, %d columns "
           "available in the matrix.\n", P->col, B->col);
    exit(1);
  }

  // Copy values
  for (i = 0; i < P->row; i++) {
    ptmp = P->map[i];
    btmp = B->map[n][i];
    for (j = 0; j < P->col; j++) {
      btmp[j] = ptmp[j];
    }
  }
}

//-----------------------------------------------------------------------------
// Dmat2Dbook: promote a double-precision M x N real matrix into a 1 x M x N
//             double precision real book.  FFT preparations are not covered.
//                                                                      
// Arguments:                                                           
//   A:  the matrix to promote                                          
//-----------------------------------------------------------------------------
dbook Dmat2Dbook(dmat *A)
{
  dbook Ab;

  Ab.pag = 1;
  Ab.row = A->row;
  Ab.col = A->col;
  Ab.data = A->data;
  Ab.pfft = 0;
  Ab.map = (double***)malloc(sizeof(double**));
  Ab.map[0] = A->map;
  
  return Ab;
}

//-----------------------------------------------------------------------------
// SumDbook: compute the sum of a double-precision real book, accounting for
//           possible 3D-FFT padding.                               
//                                                                      
// Arguments:                                                           
//   A:  the grid to sum
//-----------------------------------------------------------------------------
double SumDbook(dbook *A)
{
  int i, j, k;
  double qs;
  double *dtmp;

  qs = 0.0;
  const int klim = A->col;
  for (i = 0; i < A->pag; i++) {
    for (j = 0; j < A->row; j++) {
      dtmp = A->map[i][j];
      for (k = 0; k < klim; k++) {
	qs += dtmp[k];
      }
    }
  }

  return qs;
}

//-----------------------------------------------------------------------------
// TriCubic2Point: tricubic interpolation from a grid to a point.       
//                                                                      
// Arguments:                                                           
//   crd:        an array of coordinates for exactly one point
//   A:          a single-precision real grid
//   {i,j,k}vec: arrays of indices for X, Y, and Z coordinates of the point
//   idir:       the type of interpolation.  Set to 1 to return the value of
//               the grid at the location of the point, otherwise this function
//               will return the derivative of the field at the point in the
//               crd array, overwriting the coordinates.
//-----------------------------------------------------------------------------
static double Tricubic2Point(double* crd, fbook *A, int* ivec, int* jvec,
			     int* kvec, int idir)
{
  int i, j;
  float *atmp;
  double x, x2, x3, y, y2, y3, z, z2, z3, val;

  // These arrays are declared statically
  // to avoid having to reallocate them. 
  static double s[4], tcof[3], t[16], ucof[3], u[4], pcof[3];

  x = crd[0];
  y = crd[1];
  z = crd[2];
  x2 = x*x;
  x3 = x2*x;
  y2 = y*y;
  y3 = y2*y;
  z2 = z*z;
  z3 = z2*z;

  // Compute the interpolants in the k direction
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      atmp = A->map[ivec[i]][jvec[j]];
      s[0] = atmp[kvec[0]];
      s[1] = atmp[kvec[1]];
      s[2] = atmp[kvec[2]];
      s[3] = atmp[kvec[3]];
      tcof[0] = ONE_SIXTH*(s[3] - s[0]) + 0.5*(s[1] - s[2]);
      tcof[1] = 0.5*(s[0]+s[2]) - s[1];
      tcof[2] = -ONE_THIRD*s[0] - 0.5*s[1] + s[2] - ONE_SIXTH*s[3];
      t[i*4 + j] = tcof[0]*z3 + tcof[1]*z2 + tcof[2]*z + s[1];
    }
  }

  // Compute the interpolants in the j direction
  for (i = 0; i < 4; i++) {
    s[0] = t[4*i];
    s[1] = t[4*i+1];
    s[2] = t[4*i+2];
    s[3] = t[4*i+3];
    ucof[0] = ONE_SIXTH*(s[3] - s[0]) + 0.5*(s[1] - s[2]);
    ucof[1] = 0.5*(s[0]+s[2]) - s[1];
    ucof[2] = -ONE_THIRD*s[0] - 0.5*s[1] + s[2] - ONE_SIXTH*s[3];
    u[i] = ucof[0]*y3 + ucof[1]*y2 + ucof[2]*y + s[1];
  }

  // Compute the interpolant in the i direction
  pcof[0] = ONE_SIXTH*(u[3] - u[0]) + 0.5*(u[1] - u[2]);
  pcof[1] = 0.5*(u[0]+u[2]) - u[1];
  pcof[2] = -ONE_THIRD*u[0] - 0.5*u[1] + u[2] - ONE_SIXTH*u[3];
  val = pcof[0]*x3 + pcof[1]*x2 + pcof[2]*x + u[1];
  if (idir == 1) {
    return val;
  }

  // If we're still here, we need to compute the forces
  crd[0] = 3.0*pcof[0]*x2 + 2.0*pcof[1]*x + pcof[2];

  // Compute interpolants in the i direction
  for (i = 0; i < 4; i++) {
    s[0] = t[i];
    s[1] = t[i+4];
    s[2] = t[i+8];
    s[3] = t[i+12];
    ucof[0] = ONE_SIXTH*(s[3] - s[0]) + 0.5*(s[1] - s[2]);
    ucof[1] = 0.5*(s[0]+s[2]) - s[1];
    ucof[2] = -ONE_THIRD*s[0] - 0.5*s[1] + s[2] - ONE_SIXTH*s[3];
    u[i] = ucof[0]*x3 + ucof[1]*x2 + ucof[2]*x + s[1];
  }

  // Compute the interpolant in the j direction
  pcof[0] = ONE_SIXTH*(u[3] - u[0]) + 0.5*(u[1] - u[2]);
  pcof[1] = 0.5*(u[0]+u[2]) - u[1];
  pcof[2] = -ONE_THIRD*u[0] - 0.5*u[1] + u[2] - ONE_SIXTH*u[3];
  if (fabs(val - (pcof[0]*y3 + pcof[1]*y2 + pcof[2]*y + u[1])) > 1.0e-8) {
    printf("Oops!  Val = %12.8lf and %12.8lf\n", val,
	   pcof[0]*y3 + pcof[1]*y2 + pcof[2]*y + u[1]);
    exit(1);
  }
  crd[1] = 3.0*pcof[0]*y2 + 2.0*pcof[1]*y + pcof[2];

  // Compute 16 new interpolants in the i direction
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      s[0] = A->map[ivec[0]][jvec[i]][kvec[j]];
      s[1] = A->map[ivec[1]][jvec[i]][kvec[j]];
      s[2] = A->map[ivec[2]][jvec[i]][kvec[j]];
      s[3] = A->map[ivec[3]][jvec[i]][kvec[j]];
      tcof[0] = ONE_SIXTH*(s[3] - s[0]) + 0.5*(s[1] - s[2]);
      tcof[1] = 0.5*(s[0]+s[2]) - s[1];
      tcof[2] = -ONE_THIRD*s[0] - 0.5*s[1] + s[2] - ONE_SIXTH*s[3];
      t[i*4 + j] = tcof[0]*x3 + tcof[1]*x2 + tcof[2]*x + s[1];
    }
  }

  // Compute interpolants in the j direction
  for (i = 0; i < 4; i++) {
    s[0] = t[i];
    s[1] = t[i+4];
    s[2] = t[i+8];
    s[3] = t[i+12];
    ucof[0] = ONE_SIXTH*(s[3] - s[0]) + 0.5*(s[1] - s[2]);
    ucof[1] = 0.5*(s[0]+s[2]) - s[1];
    ucof[2] = -ONE_THIRD*s[0] - 0.5*s[1] + s[2] - ONE_SIXTH*s[3];
    u[i] = ucof[0]*y3 + ucof[1]*y2 + ucof[2]*y + s[1];
  }

  // Compute the interpolant in the k direction
  pcof[0] = ONE_SIXTH*(u[3] - u[0]) + 0.5*(u[1] - u[2]);
  pcof[1] = 0.5*(u[0]+u[2]) - u[1];
  pcof[2] = -ONE_THIRD*u[0] - 0.5*u[1] + u[2] - ONE_SIXTH*u[3];
  if (fabs(val - (pcof[0]*z3 + pcof[1]*z2 + pcof[2]*z + u[1])) > 1.0e-8) {
    printf("Oops!  Val = %12.8lf and %12.8lf\n", val,
           pcof[0]*z3 + pcof[1]*z2 + pcof[2]*z + u[1]);
    exit(1);
  }
  crd[2] = 3.0*pcof[0]*z2 + 2.0*pcof[1]*z + pcof[2];


  return val;
}

//-----------------------------------------------------------------------------
// TriInterp: tri(linear,cubic) interpolation between points and a grid.
//                                                                      
// Arguments:                                                           
//   crd:    the coordinates of the point (if idir is set to -1, the function
//           derivative at the point will be returned in crd)  
//   val:    the value of the function at the point                     
//   A:      the three-dimensional grid, complete with origin and spacing
//           information                                        
//   idir:   the direction of interpolation (0 = particle -> grid,      
//           1 = grid -> particle, -1 = grid->particle with derivatives)
//   order:  the order of interpolation (1 = linear, 3 = cubic)         
//-----------------------------------------------------------------------------
void TriInterp(double* crd, double *val, fbook *A, int idir, int order)
{
  int i, j, k, ip1, jp1, kp1;
  double dx, dy, dz, ldx, ldy, ldz;
  double *ltmp, *utmp, *iutmp;

  // Determine the grid coordinates
  dx = (crd[0] - A->orig[0]);
  dy = (crd[1] - A->orig[1]);
  dz = (crd[2] - A->orig[2]);
  ltmp = A->L.data;
  utmp = A->U.data;
  iutmp = A->invU.data;
  if (A->isortho == 1) {
    dx *= utmp[0];
    dy *= utmp[4];
    dz *= utmp[8];
    dx -= floor(dx);
    dy -= floor(dy);
    dz -= floor(dz);
    dx *= iutmp[0];
    dy *= iutmp[4];
    dz *= iutmp[8];
    dx *= ltmp[0];
    dy *= ltmp[4];
    dz *= ltmp[8];
  }
  else {
    ldx = utmp[0]*dx + utmp[1]*dy + utmp[2]*dz;
    ldy = utmp[3]*dx + utmp[4]*dy + utmp[5]*dz;
    ldz = utmp[6]*dx + utmp[7]*dy + utmp[8]*dz;
    ldx -= floor(ldx);
    ldy -= floor(ldy);
    ldz -= floor(ldz);
    dx = iutmp[0]*ldx + iutmp[1]*ldy + iutmp[2]*ldz;
    dy = iutmp[3]*ldx + iutmp[4]*ldy + iutmp[5]*ldz;
    dz = iutmp[6]*ldx + iutmp[7]*ldy + iutmp[8]*ldz;
    ldx = ltmp[0]*dx + ltmp[1]*dy + ltmp[2]*dz;
    ldy = ltmp[3]*dx + ltmp[4]*dy + ltmp[5]*dz;
    ldz = ltmp[6]*dx + ltmp[7]*dy + ltmp[8]*dz;
    dx = ldx;
    dy = ldy;
    dz = ldz;
  }
  i = dx;
  ip1 = (i == A->pag-1) ? 0 : i+1;
  dx -= i;
  j = dy;
  jp1 = (j == A->row-1) ? 0 : j+1;
  dy -= j;
  k = dz;
  kp1 = (k == A->col-1) ? 0 : k+1;
  dz -= k;
  if (idir == 0) {

    // Contribute to grid points
    if (order == 1) {

      // This is where the point value is folded into the sum.
      ldx = (*val)*(1.0-dx);
      ldy = 1.0-dy;
      ldz = 1.0-dz;
      dx *= (*val);
      A->map[i][j][k] += ldx*ldy*ldz;
      A->map[ip1][j][k] += dx*ldy*ldz;
      A->map[i][jp1][k] += ldx*dy*ldz;
      A->map[ip1][jp1][k] += dx*dy*ldz;
      A->map[i][j][kp1] += ldx*ldy*dz;
      A->map[ip1][j][kp1] += dx*ldy*dz;
      A->map[i][jp1][kp1] += ldx*dy*dz;
      A->map[ip1][jp1][kp1] += dx*dy*dz;
    }
    else {
      printf("TriInterp >> Error.  Linear interpolation of points to grids "
	     "only.\n");
      exit(1); 
    }
  }
  else {
    if (order == 1) {
      ldx = 1.0-dx;
      ldy = 1.0-dy;
      ldz = 1.0-dz;
      *val =  ((A->map[i][j][k] * ldx + A->map[ip1][j][k] * dx)*ldy +
	       (A->map[i][jp1][k] * ldx + A->map[ip1][jp1][k] * dx)*dy)*ldz + 
	((A->map[i][j][kp1] * ldx + A->map[ip1][j][kp1] * dx)*ldy +
	 (A->map[i][jp1][kp1] * ldx + A->map[ip1][jp1][kp1] * dx)*dy)*dz;
    }
    else {

      // Load all coefficients from grid; compute
      // first sixteen cubic interpolants.       
      int ivec[4], jvec[4], kvec[4];
      double loc[3];
      ivec[0] = (i == 0) ? A->pag-1 : i-1;
      jvec[0] = (j == 0) ? A->row-1 : j-1;
      kvec[0] = (k == 0) ? A->col-1 : k-1;
      ivec[1] = i;
      jvec[1] = j;
      kvec[1] = k;
      ivec[2] = ip1;
      jvec[2] = jp1;
      kvec[2] = kp1;
      ivec[3] = (ip1 < A->pag-1) ? ip1+1 : 0;
      jvec[3] = (jp1 < A->row-1) ? jp1+1 : 0;
      kvec[3] = (kp1 < A->col-1) ? kp1+1 : 0;
      loc[0] = dx;
      loc[1] = dy;
      loc[2] = dz;

      // If requested, the derivative will be   
      // computed and returned in the crd array.
      *val = Tricubic2Point(loc, A, ivec, jvec, kvec, idir);
      if (idir == -1) {
	crd[0] = loc[0];
	crd[1] = loc[1];
	crd[2] = loc[2];
      }
    }
  }
}
