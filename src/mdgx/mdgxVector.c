#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mdgxVector.h"

//-----------------------------------------------------------------------------
// CpyIVec: copy an integer vector V of length n.                       
//-----------------------------------------------------------------------------
int* CpyIVec(int* V, int n)
{
  int i;
  int* C;

  C = (int*)malloc(n*sizeof(int));
  for (i = 0; i < n; i++) {
    C[i] = V[i];
  }

  return C;
}

//-----------------------------------------------------------------------------
// CpyDVec: copies a double-precision real vector V of length n into C; C is
//          allocated and returned.                                
//-----------------------------------------------------------------------------
double* CpyDVec(double* V, int n)
{
  int i;
  double* C;

  C = (double*)malloc(n*sizeof(double));
  for (i = 0; i < n; i++) {
    C[i] = V[i];
  }

  return C;
}

//-----------------------------------------------------------------------------
// ReflectIVec: copies an integer vector V into C, assuming that C has  
//              already been allocated.  Nothing is returned.           
//-----------------------------------------------------------------------------
void ReflectIVec(int* C, int* V, int n)
{
  int i;

  for (i = 0; i < n; i++) {
    C[i] = V[i];
  }
}

//-----------------------------------------------------------------------------
// ReflectDVec: copies a double-precision real vector V into C, assuming
//              that C has already been allocated.  Nothing is returned.
//-----------------------------------------------------------------------------
void ReflectDVec(double* C, double* V, int n)
{
  int i;

  for (i = 0; i < n; i++) {
    C[i] = V[i];
  }
}

//-----------------------------------------------------------------------------
// SetIVec: set all N elements of an integer vector V to s.             
//-----------------------------------------------------------------------------
void SetIVec(int* V, int N, int s)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] = s;
  }
}

//-----------------------------------------------------------------------------
// SetDVec: set all N elements of a double-precision real vector V to s.
//-----------------------------------------------------------------------------
void SetDVec(double* V, int N, double s)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] = s;
  }
}

//-----------------------------------------------------------------------------
// CountUp: produces an integer vector of N entries 0, 1, ..., N-1      
//-----------------------------------------------------------------------------
int* CountUp(int N)
{
  int i;
  int * C;

  C = (int*)malloc(N*sizeof(int));
  for (i = 0; i < N; i++) {
    C[i] = i;
  }

  return C;
}

//-----------------------------------------------------------------------------
// ISum: computes the sum of an integer vector V of length n.           
//-----------------------------------------------------------------------------
int ISum(int* V, int n)
{
  int i, s;

  s = 0;
  for (i = 0; i < n; i++) {
    s += V[i];
  }

  return s;
}

//-----------------------------------------------------------------------------
// DetectInIVec: detect an integer im in vector V of length n.  Returns the
//               index if the integer is found, -1 if the integer cannot be
//               found.                                       
//-----------------------------------------------------------------------------
int DetectInIVec(int im, int* V, int n)
{
  int i;

  for (i = 0; i < n; i++) {
    if (V[i] == im) {
      return i;
    }
  }

  return -1;
}

//-----------------------------------------------------------------------------
// DSum: computes the sum of a double-precision real vector V of length n.
//-----------------------------------------------------------------------------
double DSum(double* V, int n)
{
  int i;
  double s;

  s = 0.0;
  for (i = 0; i < n; i++) {
    s += V[i];
  }

  return s;
}

//-----------------------------------------------------------------------------
// DAbsSum: computes the sum of absolute values of a double-precision real
//          vector V of length n.
//-----------------------------------------------------------------------------
double DAbsSum(double* V, int n)
{
  int i;
  double s;

  s = 0.0;
  for (i = 0; i < n; i++) {
    s += fabs(V[i]);
  }

  return s;
}

//-----------------------------------------------------------------------------
// Normalize: normalize a double-precision real vector V of length N.   
//-----------------------------------------------------------------------------
void Normalize(double* V, int N)
{
  int i;
  double m;

  m = 1.0/DMag(V, N);
  for (i = 0; i < N; i++) {
    V[i] *= m;
  }
}

//-----------------------------------------------------------------------------
// DAverage: computes the mean value of a double-precision real vector. 
//-----------------------------------------------------------------------------
double DAverage(double* V, int n)
{
  return DSum(V, n)/n;
}

//-----------------------------------------------------------------------------
// DStDev: computes the standard deviation of a double-precision real vector.
//-----------------------------------------------------------------------------
double DStDev(double* V, int n)
{
  int i;
  double m, m2;

  m = DAverage(V, n);
  m2 = 0.0;
  for (i = 0; i < n; i++) {
    m2 += (V[i]-m)*(V[i]-m);
  }
  if (m2 < 1.0e-8) {
    return 0.0;
  }

  return sqrt(m2/(n-1));
}

//-----------------------------------------------------------------------------
// pythag: compute the pythagorean hypotenuse length associated with a and b,
//         in a manner that is not prone to catastrophic error.
//-----------------------------------------------------------------------------
double pythag(double a, double b)
{
  double absa, absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) {
    return absa*sqrt(1.0+(absb/absa)*(absb/absa));
  }
  else {
    return absb*sqrt(1.0+(absa/absb)*(absa/absb));
  }
}

//-----------------------------------------------------------------------------
// IExtreme: integer extremum of a vector V of length N.  Setting iflag to 1
//           gets the maximum, anything else the minimum.          
//-----------------------------------------------------------------------------
int IExtreme(int* V, int N, int iflag)
{
  int i, etr;

  etr = V[0];
  if (iflag == 1) {
    for (i = 1; i < N; i++) {
      if (V[i] > etr) {
        etr = V[i];
      }
    }
  }
  else {
    for (i = 1; i < N; i++) {
      if (V[i] < etr) {
        etr = V[i];
      }
    }
  }

  return etr;
}

//-----------------------------------------------------------------------------
// DExtreme: real extremum of a vector V of length N.  Setting iflag to 1 gets
//           the maximum, anything else the minimum.          
//-----------------------------------------------------------------------------
double DExtreme(double* V, int N, int iflag)
{
  int i;
  double etr;

  etr = V[0];
  if (iflag == 1) {
    for (i = 1; i < N; i++) {
      if (V[i] > etr) {
        etr = V[i];
      }
    }
  }
  else {
    for (i = 1; i < N; i++) {
      if (V[i] < etr) {
        etr = V[i];
      }
    }
  }

  return etr;
}

//-----------------------------------------------------------------------------
// DotP: Dot product of two double-precision real vectors V1 and V2 of  
//       length N.                                                      
//-----------------------------------------------------------------------------
double DotP(double* V1, double* V2, int N)
{
  int i;
  double dp;

  dp = 0.0;
  for (i = 0; i < N; i++) {
    dp += V1[i]*V2[i];
  }

  return dp;
}

//-----------------------------------------------------------------------------
// Dot3: Dot product of two double-precision real vectors V1 and V2 of length
//       3.  Useful for coordinate or displacement dot products. 
//-----------------------------------------------------------------------------
double Dot3(double* V1, double* V2)
{
  return V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2];
}

//-----------------------------------------------------------------------------
// DMag: takes the magnitude of a double-precision real vector.         
//-----------------------------------------------------------------------------
double DMag(double* V, int n)
{
  int i;
  double m;

  m = 0.0;
  for (i = 0; i < n; i++) {
    m += V[i]*V[i];
  }

  return sqrt(m);
}

//-----------------------------------------------------------------------------
// UnitVector3: normalize a vector V of length 3.                       
//-----------------------------------------------------------------------------
void UnitVector3(double* V)
{
  double mg;

  mg = 1.0/sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
  V[0] *= mg;
  V[1] *= mg;
  V[2] *= mg;
}

//-----------------------------------------------------------------------------
// IVecAdd: add S to an integer vector V of length N.                   
//-----------------------------------------------------------------------------
void IVecAdd(int* V, int N, int S)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] += S;
  }
}

//-----------------------------------------------------------------------------
// DVecAdd: add S to a double-precision real vector V of length N.      
//-----------------------------------------------------------------------------
void DVecAdd(double* V, int N, double S)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] += S;
  }
}

//-----------------------------------------------------------------------------
// IVec2VecAdd: add vector W to vector V, assuming that V is of length N
//              and W is of length N or larger.                         
//-----------------------------------------------------------------------------
void IVec2VecAdd(int* V, int* W, int N)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] += W[i];
  }
}

//-----------------------------------------------------------------------------
// DVec2VecAdd: add vector W to vector V, assuming that V is of length N and W
//              is of length N or larger.                         
//-----------------------------------------------------------------------------
void DVec2VecAdd(double* V, double* W, int N)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] += W[i];
  }
}

//-----------------------------------------------------------------------------
// IVecMult: multiply an integer vector V of length N by S.             
//-----------------------------------------------------------------------------
void IVecMult(int* V, int N, int S)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] *= S;
  }
}

//-----------------------------------------------------------------------------
// DVecMult: multiply a double-precision real vector V of length N by S.
//-----------------------------------------------------------------------------
void DVecMult(double* V, int N, double S)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] *= S;
  }
}

//-----------------------------------------------------------------------------
// TrimAcos: trim the acos argument to be within (-1, 1).               
//-----------------------------------------------------------------------------
double TrimAcos(double acosarg)
{
  if (acosarg <= -1.0+1.0e-12 || acosarg >= 1.0-1.0e-12) {
    if (acosarg < -1.00001 || acosarg > 1.00001) {
      printf("TRIM_ACOS >> Error: ACOS argument is %12.4lf\nTRIM_ACOS >> "
             "This is erroneously large!\n", acosarg);
      exit(1);
    }
    if (acosarg >= 1.0-1.0e-12) {
      return 1.0-1.0e-12;
    }
    else {
      return -1.0+1.0e-12;
    }
  }

  return acosarg;
}

//-----------------------------------------------------------------------------
// CrossP: function for finding the cross-product cr of vectors p and q.  Note
//         that vectors p and q are assumed to be three-dimensional and only
//         the first three components of these vectors will be considered.
//-----------------------------------------------------------------------------
void CrossP(double* p, double* q, double* cr)
{
  cr[0] = p[1]*q[2] - p[2]*q[1];
  cr[1] = p[2]*q[0] - p[0]*q[2];
  cr[2] = p[0]*q[1] - p[1]*q[0];
}

//-----------------------------------------------------------------------------
// Project: function for projecting a vector r of length dim onto q.    
//          The result is returned as p.                                
//-----------------------------------------------------------------------------
void Project(double* r, double* q, double* p, int dim)
{
  int i;
  double dpval, mgq; 

  mgq = DMag(q, dim);
  dpval = DotP(r, q, dim)/(mgq*mgq);
  for (i = 0; i < dim; i++) {
    p[i] = q[i]*dpval;
  }
}

//-----------------------------------------------------------------------------
// Dihedral: function for determining a torsional angle given four atom 
//           locations in the order A B C D:                            
//                                                 D                    
//                                                /                     
//                                        --- B--C -->                  
//                                           /                          
//                                          A                           
//                                                                      
//           Above, looking down the BC axis, A and D are taken to have a
//           dihedral angle of 180 degrees.  Rotating the D atom into the 
//           screen generates takes the angle from 180 -> -180 in a negative
//           direction (170, 160, etc.).  Rotating D out of the screen
//           generates a positive trend in the dihedral angle
//           (-180, -170, -160, -150, etc.).                            
//                                                                      
// Arguments:                                                           
// [a,b,c,d]locs: cartesian coordinates of the atoms                    
//                                                                      
// Internal Variables:                                                  
// ba,cb,cd: vectors joining the named points, i.e. b --> a             
// crabcb,crcbcd,scr: vectors derived from ba,cb,cd                     
// theta: value of the dihedral angle (returned)                        
//-----------------------------------------------------------------------------
double Dihedral(double* alocs, double* blocs, double* clocs, double* dlocs)
{
  int i;
  double theta, acosarg;
  double ab[3], bc[3], cd[3], scr[3], crabbc[3], crbccd[3];

  for (i = 0; i < 3; i++) {
    ab[i] = blocs[i] - alocs[i];
    bc[i] = clocs[i] - blocs[i];
    cd[i] = dlocs[i] - clocs[i];
  }
  CrossP(ab, bc, crabbc);
  CrossP(bc, cd, crbccd);
  acosarg = DotP(crabbc, crbccd, 3)/(DMag(crabbc, 3)*DMag(crbccd, 3));
  CrossP(crabbc, crbccd, scr);
  if (DotP(scr, bc, 3) > 0.0) {
    theta = acos(TrimAcos(acosarg));
  }
  else {
    theta = -acos(TrimAcos(acosarg));
  }

  return theta;
}

//-----------------------------------------------------------------------------
// Pearson: function for determining (Pearson) correlation coefficient between
//          two vectors vec_a and vec_b of rank num_t.          
//                                                                      
// Internal Variables:                                                  
// i:                  loop control variable
// sum_[a,b,ab,aa,bb]: the sum of the vectors vec_a and vec_b, or the product
//                     thereof
// pcorr:              correlation coefficient between vectors vec_a and vec_b
//-----------------------------------------------------------------------------
double Pearson(double* vec_a, double* vec_b, int num_t)
{
  int i;
  double sum_a, sum_b, sum_ab, sum_aa, sum_bb, pcorr, va, vb, sq;

  sum_a = 0.0;
  sum_b = 0.0;
  sum_ab = 0.0;
  sum_aa = 0.0;
  sum_bb = 0.0;
  for (i = 0; i < num_t; i++) {
    va = vec_a[i];
    vb = vec_b[i];
    sum_a += va;
    sum_b += vb;
    sum_ab += va*vb;
    sum_aa += va*va;
    sum_bb += vb*vb;
  }

  sq = (num_t*sum_aa - sum_a*sum_a)*(num_t*sum_bb - sum_b*sum_b);
  if (sq < 1.0e-12) {
    return 0.0;
  }
  pcorr = (num_t*sum_ab - sum_a*sum_b) / sqrt(sq);

  return pcorr;
}

//-----------------------------------------------------------------------------
// SumTrace3: takes the sum of the trace of a 3x3 matrix A, represented by a
//            vector of 9 elements which give the matrix elements in the order
//            (1,1), (1,2), (1,3), (2,1,) ...              
//-----------------------------------------------------------------------------
double SumTrace3(double* A)
{
  return A[0] + A[4] + A[8];
}

//-----------------------------------------------------------------------------
// PascalTriangle: compute Pascal's triangle of a specified order.  Returns an
//                 array with the required coefficients.     
//
// Arguments:                                                           
//   n:      the order of the triangle to compute                       
//-----------------------------------------------------------------------------
double* PascalTriangle(int n)
{
  int i, j, pvalj, pvaljm1;
  double* pval;

  pval = (double*)calloc(n, sizeof(double));
  pval[0] = 1.0;
  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++) {
      pvalj = pval[j];
      pval[j] = (j == 0) ? pval[0] : pval[j] + pvaljm1;
      pvaljm1 = pvalj;
    }
  }

  return pval;
}

//-----------------------------------------------------------------------------
// VecRMSD: function for computing the RMSD between two vectors of size n,
//          returning the result as a double.                        
//-----------------------------------------------------------------------------
double VecRMSD(double* va, double* vb, int n)
{
  int i;
  double dx, rmsdev;

  rmsdev = 0.0;
  for (i = 0; i < n; i++) {
    dx = vb[i] - va[i];
    rmsdev += dx*dx;
  }
  rmsdev = rmsdev / n;

  return sqrt(rmsdev);
}

//-----------------------------------------------------------------------------
// UnitNormal2Plane: compute a unit normal vector unv to a plane given points
//                   A, B, and C in the plane.                   
//-----------------------------------------------------------------------------
void UnitNormal2Plane(double* unv, double* A, double* B, double* C)
{
  double BA[3], BC[3];

  // Compute the vectors BA and BC
  BA[0] = A[0] - B[0];
  BA[1] = A[1] - B[1];
  BA[2] = A[2] - B[2];
  BC[0] = C[0] - B[0];
  BC[1] = C[1] - B[1];
  BC[2] = C[2] - B[2];

  // Take the cross product and normalize
  CrossP(BA, BC, unv);
  UnitVector3(unv);
}
