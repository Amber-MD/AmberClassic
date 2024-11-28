#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"

/***=======================================================================***/
/*** CpyIVec: copy an integer vector V of length n.                        ***/
/***=======================================================================***/
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

/***=======================================================================***/
/*** CountUp: produces an integer vector of length N with entries 1,2,..,N ***/
/***=======================================================================***/
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

/***=======================================================================***/
/*** ISum: computes the sum of an integer vector V of length n.            ***/
/***=======================================================================***/
int ISum(int* V, int n)
{
  int i, s;

  s = 0;
  for (i = 0; i < n; i++) {
    s += V[i];
  }

  return s;
}

/***=======================================================================***/
/*** DSum: computes the sum of a double-precision real vector V of length  ***/
/***       n.                                                              ***/
/***=======================================================================***/
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

/***=======================================================================***/
/*** Normalize: normalize a double-precision real vector V of length N.    ***/
/***=======================================================================***/
void Normalize(double* V, int N)
{
  int i;
  double m;

  m = 1.0/DMag(V, N);
  for (i = 0; i < N; i++) {
    V[i] *= m;
  }
}

/***=======================================================================***/
/*** DAverage: computes the mean value of a double-precision real vector.  ***/
/***=======================================================================***/
double DAverage(double* V, int n)
{
  return DSum(V, n)/n;
}

/***=======================================================================***/
/*** DStDev: computes the standard deviation of a double-precision real    ***/
/***         vector.                                                       ***/
/***=======================================================================***/
double DStDev(double* V, int n)
{
  int i;
  double m, m2;

  m = DAverage(V, n);
  m2 = 0.0;
  for (i = 0; i < n; i++) {
    m2 += (V[i]-m)*(V[i]-m);
  }

  return sqrt(m2/(n-1));
}

/***=======================================================================***/
/*** CpyDVec: copies a double-precision real vector V into C.              ***/
/***=======================================================================***/
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

/***=======================================================================***/
/*** IExtreme: integer extremum of a vector V of length N.  Setting iflag  ***/
/***           to 1 gets the maximum, anything else the minimum.           ***/
/***=======================================================================***/
int IExtreme(int* V, int N, int iflag)
{
  int i, etr;

  etr = V[0];
  for (i = 1; i < N; i++) {
    if (iflag == 1 && V[i] > etr) {
      etr = V[i];
    }
    else if (V[i] < etr) {
      etr = V[i];
    }
  }

  return etr;
}

/***=======================================================================***/
/*** DExtreme: real extremum of a vector V of length N.  Setting iflag     ***/
/***           to 1 gets the maximum, anything else the minimum.           ***/
/***=======================================================================***/
double DExtreme(double* V, int N, int iflag)
{
  int i, etr;

  etr = V[0];
  for (i = 1; i < N; i++) {
    if (iflag == 1 && V[i] > etr) {
      etr = V[i];
    }
    else if (V[i] < etr) {
      etr = V[i];
    }
  }

  return etr;
}

/***=======================================================================***/
/*** SetVec: set all N elements of an integer vector V to s.               ***/
/***=======================================================================***/
void SetIVec(int* V, int N, int s)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] = s;
  }
}

/***=======================================================================***/
/*** DotP: Dot product of two double-precision real vectors V1 and V2 of   ***/
/***       length N.                                                       ***/
/***=======================================================================***/
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

/***=======================================================================***/
/*** DMag: takes the magnitude of a double-precision real vector.          ***/
/***=======================================================================***/
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

/***=======================================================================***/
/*** AddToIVec: add S to an integer vector V of length N.                  ***/
/***=======================================================================***/
void AddToIVec(int* V, int N, int S)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] += S;
  }
}

/***=======================================================================***/
/*** TrimAcos: trim the acos argument to be within (-1, 1).                ***/
/***=======================================================================***/
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

/***=======================================================================***/
/*** CrossP: function for finding the cross-product cr of vectors p and q. ***/
/***         Note that vectors p and q are assumed to be three-dimensional ***/
/***         and only the first three components of these vectors will be  ***/
/***         considered.                                                   ***/
/***=======================================================================***/
void CrossP(double p[3], double q[3], double cr[3])
{
  cr[0] = p[1]*q[2]-p[2]*q[1];
  cr[1] = p[2]*q[0]-p[0]*q[2];
  cr[2] = p[0]*q[1]-p[1]*q[0];
}

/***=======================================================================***/
/*** Project: function for projecting a vector r of length dim onto q.     ***/
/***          The result is returned as p.                                 ***/
/***=======================================================================***/
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

/***=======================================================================***/
/*** Dihedral: function for determining a torsional angle given four atom  ***/
/***           locations in the order A B C D:                             ***/
/***                                                 D                     ***/
/***                                                /                      ***/
/***                                        --- B--C -->                   ***/
/***                                           /                           ***/
/***                                          A                            ***/
/***                                                                       ***/
/***           Above, looking down the BC axis, A and D are taken to have  ***/
/***           a dihedral angle of 180 degrees.  Rotating the D atom into  ***/
/***           the screen generates takes the angle from 180 -> -180 in a  ***/
/***           negative direction (170, 160, etc.).  Rotating D out of     ***/
/***           the screen generates a positive trend in the dihedral angle ***/
/***           (-180, -170, -160, -150, etc.).                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/*** [a,b,c,d]locs: cartesian coordinates of the atoms                     ***/
/***                                                                       ***/
/*** Internal Variables:                                                   ***/
/*** ba,cb,cd: vectors joining the named points, i.e. b --> a              ***/
/*** crabcb,crcbcd,scr: vectors derived from ba,cb,cd                      ***/
/*** theta: value of the dihedral angle (returned)                         ***/
/***=======================================================================***/
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

/***=======================================================================***/
/*** Pearson: function for determining (Pearson) correlation coefficient   ***/
/***          between two vectors vec_a and vec_b of rank num_t.           ***/
/***                                                                       ***/
/*** Internal Variables:                                                   ***/
/*** i: loop control variable                                              ***/
/*** sum_[a,b,ab,aa,bb]: the sum of the vectors vec_a, vec_b, or a product ***/
/***                     thereof                                           ***/
/*** pcorr: the correlation coefficient between vectors vec_a and vec_b    ***/
/***=======================================================================***/
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
