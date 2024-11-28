/*
 * newton.c: Newton-Raphson minimizer.
 *
 * Provided by Russell A. Brown (russ.brown@sun.com)
 *
 * Line search added by David A. Case (case@scripps.edu)
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "sff.h"
#include "memutil.h"
#include "timer.h"

#if defined(MPI) || defined(SCALAPACK)
#include "mpi.h"
#endif

/* Here are offsets into the ScaLAPACK descriptor array. */

#define DTYPE_ (0)
#define CTXT_ (1)
#define M_ (2)
#define N_ (3)
#define MB_ (4)
#define NB_ (5)
#define RSRC_ (6)
#define CSRC_ (7)
#define LLD_ (8)
#define DLEN_ (9)

REAL_T seconds(void);

INT_T get_blocksize(void);

INT_T get_mytaskid(void);

INT_T get_numtasks(void);

INT_T get_nr_debug(void);

void dgemm_(char *, char *, INT_T *, INT_T *, INT_T *, REAL_T *,
            REAL_T *, INT_T *, REAL_T *, INT_T *, REAL_T *, REAL_T *,
            INT_T *);

void dposv_(char *, INT_T *, INT_T *, REAL_T *,
            INT_T *, REAL_T *, INT_T *, INT_T *);

INT_T dsyev_(char *, char *, INT_T *, REAL_T *, INT_T *,
             REAL_T *, REAL_T *, INT_T *, INT_T *);

REAL_T dnrm2_(INT_T *, REAL_T *, INT_T *);

REAL_T ddot_(INT_T *, REAL_T *, INT_T *, REAL_T *, INT_T *);

void dcopy_(INT_T *, REAL_T *, INT_T *, REAL_T *, INT_T *);

void daxpy_(INT_T *, REAL_T *, REAL_T *, INT_T *, REAL_T *, INT_T *);

void sl_init_(INT_T *, INT_T *, INT_T *);

void blacs_gridinfo_(INT_T *, INT_T *, INT_T *, INT_T *, INT_T *);

void blacs_barrier_(INT_T *, char *);

void blacs_gridexit_(INT_T *);

INT_T numroc_(INT_T *, INT_T *, INT_T *, INT_T *, INT_T *);

void descinit_(INT_T *, INT_T *, INT_T *, INT_T *, INT_T *,
               INT_T *, INT_T *, INT_T *, INT_T *, INT_T *);

void pdgemr2d_(INT_T *, INT_T *,
               REAL_T *, INT_T *, INT_T *, INT_T *,
               REAL_T *, INT_T *, INT_T *, INT_T *, INT_T *);

void pdgemm_(char *, char *, INT_T *, INT_T *, INT_T *, REAL_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *, REAL_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *);

void pdposv_(char *, INT_T *, INT_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *, INT_T *);

void pdsyev_(char *, char *, INT_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *, REAL_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *,
             REAL_T *, INT_T *, INT_T *);

INT_T MPI_Finalize(void);

void exit(int);

INT_T conjgrad(REAL_T *, INT_T *, REAL_T *,
               REAL_T(*func) (REAL_T *, REAL_T *, INT_T *),
               REAL_T *, REAL_T *, INT_T *);

REAL_T mme(REAL_T *, REAL_T *, INT_T *);

INT_T myroc(int, int, int, int);

size_t adr1d(int, int, int);

size_t adr2d(int, int, int, int, int, int, int);

REAL_T *ptr1d(REAL_T *, int *, int);

REAL_T *ptr2d(REAL_T *, int *, int, int);


/***********************************************************************
                            LEVEL()
************************************************************************/

/*
 * Prepare a matrix to adjust the Hessian matrix to put the overall
 * translations and rotations at a relatively high frequency.
 *
 * See Nguyen and Case, J. Phys. Chem. 89, 4020 (1985).
 *
 * Calling parameters are as follows:
 *
 * h      the Hessian matrix to be modified
 * x[n]   the atomic (x,y,z) coordinates
 * m      the atomic masses
 * d      unit vector matrix
 * n      the number of coordinates
 * natom  the number of atoms
 */

static
void level(REAL_T * x, REAL_T * m, REAL_T * d, INT_T n, INT_T natom)
{
   REAL_T xx, yy, zz, mass, xsum, ysum, zsum, msum;
   REAL_T vnor1, vnor2, vnor3, vnor4;
   INT_T i, iat;

   /* Normalize the translational and rotational degrees of freedom. */

   xsum = ysum = zsum = msum = 0.0;

   for (i = 1; i <= natom; i++) {

      xx = x[3 * (i - 1) + 1];
      yy = x[3 * (i - 1) + 2];
      zz = x[3 * (i - 1) + 3];
      mass = m[i];

      xsum += xx * xx * mass;
      ysum += yy * yy * mass;
      zsum += zz * zz * mass;
      msum += mass;
   }

   vnor1 = sqrt(xsum + ysum);
   vnor2 = sqrt(ysum + zsum);
   vnor3 = sqrt(xsum + zsum);
   vnor4 = sqrt(msum);

   /*
    * Set up the translational and rotational unit vectors.
    * See Eckart, Phys. Rev. 47, 552 (1935).
    */

   for (i = 1; i <= 6 * n; i++) {
      d[i] = 0.0;
   }

   iat = 1;
   for (i = 1; i <= n; i += 3) {

      xx = x[3 * (iat - 1) + 1];
      yy = x[3 * (iat - 1) + 2];
      zz = x[3 * (iat - 1) + 3];
      mass = sqrt(m[iat]);

      /*
       * For reference, here is the code for d as a matrix:
       *
       *      d[i][1] = mass / vnor4;
       *      d[i+1][2] = mass / vnor4;
       *      d[i+2][3] = mass / vnor4;
       *
       *      d[i][4] = -mass * yy / vnor1;
       *      d[i+1][4] = mass * xx / vnor1;
       *
       *      d[i+1][5] = -mass * zz / vnor2;
       *      d[i+2][5] = mass * yy / vnor2;
       *
       *      d[i][6] = mass * zz / vnor3;
       *      d[i+2][6] = -mass * xx / vnor3;
       */

      d[i] = mass / vnor4;
      d[n + i + 1] = mass / vnor4;
      d[2 * n + i + 2] = mass / vnor4;

      d[3 * n + i] = -mass * yy / vnor1;
      d[3 * n + i + 1] = mass * xx / vnor1;

      d[4 * n + i + 1] = -mass * zz / vnor2;
      d[4 * n + i + 2] = mass * yy / vnor2;

      d[5 * n + i] = mass * zz / vnor3;
      d[5 * n + i + 2] = -mass * xx / vnor3;

      iat++;
   }
}

/***********************************************************************
                            LINEMIN()
************************************************************************/

/*
 *  Find the value of alpha such that the energy at position x + alpha*d
 *  is minimized.  On input, one needs the current position, x[], the
 *  search direction d[], an initial guess for alpha, a pointer to the
 *  function to be used, a maximum number of evaluations to allow,
 *  and relative (tol[0]) and absolute (tol[1]) tolerances.
 *
 *  On output, x[] will be updated to the best position, alpha will
 *  have its optimal value, and evlmax will give the number of searches
 *  actually carried out.  The other parameters will be unchanged.
 *
 *  Based on a routine by H.T. Lau, in "A Numerical Library in C for
 *  Scientists and Engineers", (Boca Raton: CRC Press, 1995).
 */

void linemin(INT_T * n, REAL_T x[], REAL_T d[], REAL_T * alpha,
             REAL_T(*func) (REAL_T *, REAL_T *, int *),
             INT_T * evlmax, REAL_T tol[])
{
   int ncopy, one, evl, notinint, minustwo, mytaskid;
   REAL_T *x0, *g;
   REAL_T alpha0, y, nd, reltol, abstol, f0, df0, oldf, olddf, w, z, f1,
       df1, eps, aid, f, df;

   mytaskid = get_mytaskid();
   one = 1;
   minustwo = -2;
   ncopy = *n;
   x0 = vector(0, ncopy);
   g = vector(0, ncopy);
   nd = dnrm2_(n, d, &one);
   reltol = tol[0];
   abstol = tol[1];
   evl = 0;
   alpha0 = 0.0;
   f0 = (*func) (x, g, &minustwo);
   if (mytaskid == 0) {
      printf("For alpha = %10.5f energy = %20.10f\n", alpha0, f0);
      fflush(stdout);
   }
   df0 = ddot_(n, d, &one, g, &one);
   oldf = f0;
   olddf = df0;
   y = *alpha;
   dcopy_(n, x, &one, x0, &one);
   daxpy_(n, alpha, d, &one, x0, &one);
   f1 = (*func) (x0, g, &minustwo);
   if (mytaskid == 0) {
      printf("For alpha = %10.5f energy = %20.10f\n", y, f1);
      fflush(stdout);
   }
   df1 = ddot_(n, d, &one, g, &one);
   notinint = 1;
   dcopy_(n, x, &one, x0, &one);
   eps = (dnrm2_(n, x, &one) * reltol + abstol) / nd;

   while (1) {
      if (notinint)
         notinint = (df1 < 0.0);
      aid = *alpha;
      if (df1 >= 0.0) {
         /* cubic interpolation */
         z = 3.0 * (oldf - f1) / (*alpha) + olddf + df1;
         w = sqrt(z * z - olddf * df1);
         *alpha *= 1.0 - (df1 + w - z) / (df1 - olddf + w + w);
         if (*alpha < eps)
            *alpha = eps;
         else if (aid - (*alpha) < eps)
            *alpha = aid - eps;
      } else {
         if (notinint) {
            alpha0 = *alpha = y;
            olddf = df1;
            oldf = f1;
         } else {
            *alpha *= 0.5;
         }
      }
      y = (*alpha) + alpha0;
      dcopy_(n, x0, &one, x, &one);
      daxpy_(n, &y, d, &one, x, &one);
      eps = (dnrm2_(n, x, &one) * reltol + abstol) / nd;
      f = (*func) (x, g, &minustwo);
      if (mytaskid == 0) {
         printf("For alpha = %10.5f energy = %20.10f\n", y, f);
         fflush(stdout);
      }
      df = ddot_(n, d, &one, g, &one);
      evl++;
      if (evl > *evlmax)
         break;
      if (notinint || df > 0.0) {
         df1 = df;
         f1 = f;
      } else {
         alpha0 = y;
         *alpha = aid - (*alpha);
         olddf = df;
         oldf = f;
      }
      if (*alpha < 2.0 * eps)
         break;
   }
   *alpha = y;
   *evlmax = evl;
   free_vector(g, 0, ncopy);
   free_vector(x0, 0, ncopy);
   return;
}

/***********************************************************************
                            NEWTON()
************************************************************************/

/*
 * Newton-Raphson minimizer
 *
 * Calling parameters are as follows:
 *
 * x[n]    contains the variables to be updated
 * n       is the number of variables
 * f       will contain the value of the final value of the function f
 * func    pointer to function that computes the energy, gradient and hessian
 * rmsgrad quit when the rms of the gradient is less than rmsgrad
 * nradd   floating point value to add to the diagonal of the hessian
 * maxiter quit when iter exceeds maxiter
 *
 *
 *
 *     return codes:
 *             >0    converged, final iteration number
 *             -1    exceeded the maximum number of iterations
 *             -2    Hessian was not positive definite; quit, since
 *                   bad things are likely to happen.  The last structure
 *                   examined will be in the x array.
 *             -3    the eigenvalue calculation did not converge, so quit
 *
 */

INT_T newton(REAL_T x[], INT_T * n, REAL_T * f,
             REAL_T(*func1) (REAL_T *, REAL_T *, int *),
             REAL_T(*func2) (REAL_T *, REAL_T *, REAL_T *, REAL_T *,
                             REAL_T *, INT_T *, INT_T *, INT_T *,
                             INT_T *, INT_T *, INT_T *, INT_T *,
                             INT_T *, INT_T *, char *),
             REAL_T * rmsgrad, REAL_T * nradd, INT_T * maxiter)
{
   REAL_T *d = NULL, *g = NULL, *m = NULL, *h = NULL, *grad = NULL;
   REAL_T sumg, dgrad;
   REAL_T tnewton1, t1, t2;
   int niter, i, j, natom, ret_val, mcopy, numtasks, mytaskid;
   int evlmax, maxit, lwork, nr_debug, gridim;
   int one = 1, info = 0, six = 6;
   int context_PxQ = -1, context_1x1 = -1, context_Nx1 = -1;
   int descH_PxQ[DLEN_], descG_PxQ[DLEN_], descG_1x1[DLEN_];
   char uplo, jobz, transa, transb;
   REAL_T scale = 3500.0, dblone = 1.0;
   REAL_T alphal, toll[2], sqrn, rmsg, fret, dfpred;
   size_t ncopy, in;
   char *name;

#ifdef SCALAPACK
   int zero = 0;
   int myrow, mycol, nprow, npcol;
   int bs3, ierror, nb, np, sizesytrd;
   int lldH_PxQ, lldG_PxQ, lldG_1x1, lldD_PxQ, lldD_1x1;
   size_t locpH_PxQ, locpG_PxQ, locpG_1x1, locpD_PxQ, locpD_1x1;
   size_t locqH_PxQ, locqG_PxQ, locqG_1x1, locqD_PxQ, locqD_1x1;
   size_t sizeH_PxQ, sizeG_PxQ, sizeG_1x1, sizeD_PxQ, sizeD_1x1;
   int descD_PxQ[DLEN_], descD_1x1[DLEN_];
   REAL_T *ptr, *reductarr = NULL, *dee = NULL, *work = NULL;
#endif

/*
 *   If DSYEV is defined, then lapack routines are used to get the low
 *   eigenvalues of the Hessian, so that it can be updated to ensure that
 *   is it positive definite.  This could be a part of a more robust NR
 *   minimization procedure, but it still needs work.
 *
 */

#undef DSYEV
#ifdef DSYEV
   REAL_T *h2, *roots;
#endif

   t1 = seconds();
   tnewton1 = t1;
   name = (char *) malloc((*n) * sizeof(char));
   assert( name );

   /* Get mytaskid and, if SCALAPACK is defined, numtasks. */

   mytaskid = get_mytaskid();

#ifdef SCALAPACK
   numtasks = get_numtasks();
#endif

   /* Get the print control variable. */

   nr_debug = get_nr_debug();


   /* All arrays, even the x array, will be indexed from 1 to *n. */

   --x;

   /*
    * Allocate some dynamic vectors and matrices.  Full copies of
    * the m and d arrays are allocated for each task.
    */

   mcopy = *n;
   ncopy = *n;
   m = vector(1, ncopy);        /* atomic masses */

#ifndef SCALAPACK

   /* If SCALAPACK is not defined, allocate full copies of g, d, grad and h. */

   g = vector(1, ncopy);        /* gradient vector */
   grad = vector(1, ncopy);     /* copy of gradient vector */
   d = vector(1, 6 * ncopy);    /* scratch for level and dsyev_ */
   h = vector(0, ncopy * ncopy);        /* hessian (linearized) */

#else

   /*
    * If SCALAPACK is defined, allocate distributed copies of g, d, grad
    * and h, as well as a vector to be used by the MPI_Allreduce function.
    *
    * Create a general context.  Although context_PxQ does comprise all
    * of the processes, it appears that the general context must be
    * distinct from the context(s) of the matrices that participate
    * in the redistribution via pdgemr2d.
    *
    * The topologic layout of the general context does not appear
    * to be important, so this general context is row cyclic.
    *
    * It appears to be important that the most general context be created
    * first, followed by successively less general contexts.  For this
    * code the correct order of creation is Nx1 then PxQ then 1x1.  Each
    * context is subsumed by the earlier created contexts.  I don't know
    * what to do about overlapping contexts where one does not sumsume
    * the other.  Failure to the most general context first leads to
    * synchronization errors.
    */

   sl_init_(&context_Nx1, &numtasks, &one);

   /* Calculate the dimensions of the largest possible square process grid. */

   gridim = (int) (sqrt((REAL_T) numtasks) + 0.5);
   if (gridim * gridim > numtasks) {
      gridim -= 1;
   }
   if (gridim <= 0) {
      gridim = 1;
   }

   /*
    * Initialize the process grid for block cyclic distribution of matrices.
    * Note that a returned context of -1 indicates that the task is not
    * active on the process grid.
    */

   sl_init_(&context_PxQ, &gridim, &gridim);

   /* Initialize the process grid for a single (1x1) process. */

   sl_init_(&context_1x1, &one, &one);

   /*
    * Get the number of rows and columns on the block cyclic (PxQ) process
    * grid, as well as this task's row and column on the grid.
    */

   blacs_gridinfo_(&context_PxQ, &nprow, &npcol, &myrow, &mycol);

   /*
    * Get the blocksize for a square block.  Because in egb2 the
    * loop index i selects three rows of the Hessian, multiply
    * the blocksize by three.
    */

   bs3 = 3 * get_blocksize();

   /*
    * If this task is on the process grid, set up the array descriptors.
    * If this task isn't on the process grid, set descD_PxQ[CTXT_],
    * descG_PxQ[CTXT_] and descH_PxQ[CTXT_] to -1.  These values will
    * be used by pdgemr2d to determine activity on the grid.
    */

   if (context_PxQ >= 0) {

      /*
       * Allocate then distribute the Hessian matrix h, the shift
       * matrix dee and the gradient vector g on the block cyclic
       * process grid.  The gradient vector g does not need to be
       * allocated for columns other than column 0, but the descinit_
       * function must be called for all columns.
       *
       * The numroc_ function is used to calculate the number of matrix
       * elements that are distributed across a PxQ processor grid.
       */

      locpD_PxQ = numroc_(&mcopy, &bs3, &myrow, &zero, &nprow);
      locqD_PxQ = numroc_(&six, &bs3, &mycol, &zero, &npcol);
      sizeD_PxQ = locpD_PxQ * locqD_PxQ;
      lldD_PxQ = locpD_PxQ;
      descinit_(descD_PxQ, &mcopy, &six, &bs3, &bs3,
                &zero, &zero, &context_PxQ, &lldD_PxQ, &info);
      dee = vector(1, sizeD_PxQ);

      locpG_PxQ = numroc_(&mcopy, &bs3, &myrow, &zero, &nprow);
      locqG_PxQ = numroc_(&one, &bs3, &mycol, &zero, &npcol);
      sizeG_PxQ = locpG_PxQ * locqG_PxQ;
      lldG_PxQ = locpG_PxQ;
      descinit_(descG_PxQ, &mcopy, &one, &bs3, &bs3,
                &zero, &zero, &context_PxQ, &lldG_PxQ, &info);
      if (mycol == 0) {
         g = vector(1, sizeG_PxQ);
      }

      locpH_PxQ = numroc_(&mcopy, &bs3, &myrow, &zero, &nprow);
      locqH_PxQ = numroc_(&mcopy, &bs3, &mycol, &zero, &npcol);
      sizeH_PxQ = locpH_PxQ * locqH_PxQ;
      lldH_PxQ = locpH_PxQ;
      descinit_(descH_PxQ, &mcopy, &mcopy, &bs3, &bs3,
                &zero, &zero, &context_PxQ, &lldH_PxQ, &info);
      h = vector(0, sizeH_PxQ);

   } else {
      descD_PxQ[CTXT_] = -1;
      descG_PxQ[CTXT_] = -1;
      descH_PxQ[CTXT_] = -1;
   }

   /*
    * Get the number of rows and columns on the single process grid,
    * as well as this task's row and column on the grid.
    */

   blacs_gridinfo_(&context_1x1, &nprow, &npcol, &myrow, &mycol);

   /*
    * If this task is on the process grid, set up the array descriptors.
    * If this task isn't on the process grid, set descD_1x1[CTXT_]
    * descG_1x1[CTXT_] and descH_1x1[CTXT_] to -1.  These values will
    * be used by pdgemr2d to determine activity on the grid.
    */

   if (context_1x1 >= 0) {

      /*
       * Allocate then distribute the matrix d, the gradient vector grad
       * and the hessian matrix hess on the single process grid.  The
       * descinit_ function is called for these arrays for only the task
       * that is active on the grid.
       *
       * Also, for the other tasks that are not active, the gradient vector
       * grad is allocated because the copy of grad from the single process
       * grid will be broadcast to the other copies.  Thesizes of the gradient
       * vector grad for these other tasks is determined by the ncopy variable.
       *
       * The numroc_ function is used to calculate the number of matrix
       * elements that are distributed across a 1x1 processor grid.
       */

      locpD_1x1 = numroc_(&mcopy, &bs3, &myrow, &zero, &nprow);
      locqD_1x1 = numroc_(&six, &bs3, &mycol, &zero, &npcol);
      sizeD_1x1 = locpD_1x1 * locqD_1x1;
      lldD_1x1 = locpD_1x1;
      descinit_(descD_1x1, &mcopy, &six, &bs3, &bs3,
                &zero, &zero, &context_1x1, &lldD_1x1, &info);
      d = vector(1, sizeD_1x1);

      locpG_1x1 = numroc_(&mcopy, &bs3, &myrow, &zero, &nprow);
      locqG_1x1 = numroc_(&one, &bs3, &mycol, &zero, &npcol);
      sizeG_1x1 = locpG_1x1 * locqG_1x1;
      lldG_1x1 = locpG_1x1;
      descinit_(descG_1x1, &mcopy, &one, &bs3, &bs3,
                &zero, &zero, &context_1x1, &lldG_1x1, &info);
      grad = vector(1, sizeG_1x1);

   } else {
      grad = vector(1, ncopy);
      descD_1x1[CTXT_] = -1;
      descG_1x1[CTXT_] = -1;
   }

   /* Allocate the temporary array that is used in reductions. */

   reductarr = vector(1, ncopy);

#endif                          /* ifdef SCALAPACK */


   /* Initialize rms cutoff and iteration counter. */

   sqrn = sqrt((float) ncopy);
   dgrad = (*rmsgrad) * (*rmsgrad) * ncopy;
   niter = 0;

   /* Here is the minimization loop. */

   do {

      /* Increment the iteration counter. */

      ++niter;

      /*
       * For non-ScaLAPACK execution, set some variables to values that
       * will select the proper sections of code below.
       */

#ifndef SCALAPACK

      gridim = 1;
      context_PxQ = 0;
      context_1x1 = 1;

#endif

      t2 = seconds();
      *tnewtonOther += t2 - t1;
      t1 = t2;
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

      /*
       *   Compute the function value in "f" and its
       *   gradient (with respect to x) in "grad".
       */

      *f = (*func1) (&x[1], &grad[1], &niter);

      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

      t2 = seconds();
      *tnewtonMME += t2 - t1;
      t1 = t2;
      if (nr_debug) {
         if (mytaskid == 0) {
            printf("\nmme time = %10.2f\n\n", t2 - t1);
            fflush(stdout);
         }
      }

      /*
       * Calculate the rms value (squared) of the gradient vector,
       * and quit if that value is below the rms cutoff.
       */

      sumg = 0.0;
      for (i = 1; i <= ncopy; i++) {
         sumg += grad[i] * grad[i];
      }

      if (sumg < dgrad) {
         break;
      }

      t2 = seconds();
      *tnewtonOther += t2 - t1;
      t1 = t2;

      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

      /*
       * Compute the function value in "f",
       * its gradient (with respect to x) in "g", and
       * its hessian (with respect to x) in "h".
       * The atomic masses are returned in "m".
       * The number of atoms is returned in "natom".
       *
       * The grad, descG_PxQ, descG_1x1, descH_PxQ,
       * calling parameters supply ScaLAPACK information,
       * or are dummy arguments for non-ScaLAPACK execution.
       *
       * The &x[1], &g[1] and &m[1] calling parameters
       * map from 1..n indexing in this newton function
       * to 0..n-1 indexing in *func (which defaults to mme2).
       * This technique is not used for h because it must be
       * indexed from 1 to n.
       */

      *f = (*func2) (&x[1], &g[1], h, &m[1],
                     &grad[1], descG_PxQ, descG_1x1, descH_PxQ,
                     &context_PxQ, &context_1x1, &context_Nx1,
                     &gridim, &natom, &niter, name);

      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

      t2 = seconds();
      *tnewtonMME2 += t2 - t1;
      if (nr_debug) {
         if (mytaskid == 0) {
            printf("\nmme2 time = %10.2f\n\n", t2 - t1);
            fflush(stdout);
         }
      }
      t1 = t2;

      /*
       * "Level-shift" the eigenvalues that correspond to overall
       * translation and rotation of the molecule to high frequencies.
       * The Hessian matrix h is not altered by the level function;
       * rather, a matrix multiply is performed after return from level.
       *
       * Before d was converted from a matrix to a vector, the syntax was:
       *
       *      for (i = 0; i < ncopy; i++)
       *        for (j = 0; j < ncopy; j++)
       *          for (ic = 0; ic < 6; ic++)
       *            h[i][j] += 3500.0 * d[i][ic] * d[j][ic];
       *
       * For ScaLAPACK, calculate the level shift matrix only for context_1x1
       * because any task that is not on the 1x1 grid does not have valid data.
       *
       * For non-ScaLAPACK execution, the context_1x1 variable has been set to 1
       * above so that the level function will be called.
       */

      if (context_1x1 >= 0) {
         level(x, m, d, ncopy, natom);
         t2 = seconds();
         *tnewtonLevel += t2 - t1;
         if (nr_debug) {
            if (mytaskid == 0) {
               printf("level time = %10.2f\n\n", t2 - t1);
               fflush(stdout);
            }
         }
         t1 = t2;
      }

      /*
       * The ScaLAPACK pdgemm_ function appears to quit unexpectedly
       * for large matrices on a 1x1 process grid, so bypass ScaLAPACK
       * and use the LAPACK dgemm_ function instead.  Also use dposv_
       * and dsyev_ instead of pdposv_ and pdsyev_ under these conditions.
       *
       * The correct test is (gridim == 1), not (nprow == 1 && npcol == 1)
       * since processes that aren't on the 1x1 grid have nprow == npcol == -1,
       * which would direct control to pdgemr2d_ (below) that would hang because
       * it would not be called from all processes, specifically not from
       * the process that is on the 1x1 grid and has nprow == npcol == 1.
       *
       * Non-ScaLAPACK execution will select this section of code because
       * gridim and context_PxQ were set to 1 and 0, respectively, above.
       */

      if (gridim == 1) {
         if (context_PxQ >= 0) {

            /* Calculate the level shift matrix. */

            level(x, m, d, ncopy, natom);

            t2 = seconds();
            *tnewtonLevel += t2 - t1;
            if (nr_debug) {
               if (mytaskid == 0) {
                  printf("level time = %10.2f\n\n", t2 - t1);
                  fflush(stdout);
               }
            }
            t1 = t2;

            /* Multiply the Hessian by the level shift matrix d. */

            transa = 'N';
            transb = 'T';
            dgemm_(&transa, &transb, &mcopy, &mcopy, &six, &scale,
                   &d[1], &mcopy, &d[1], &mcopy, &dblone, h, &mcopy);

            if (nr_debug) {
               if (mytaskid == 0) {
                  printf("dgemm time = %10.2f\n\n", t2 - t1);
                  fflush(stdout);
               }
            }
#ifndef DSYEV
            /* Modify the diagonal elements of the Hessian. */
            /*  Do this here using the input value, rather than trying
               to use dsyev() to figure out the proper value   */

            if (mytaskid == 0) {
               printf(" adding %10.5f to diagonal of the hessian\n",
                      *nradd);
               fflush(stdout);
            }
            for (i = 0; i < ncopy; i++) {
               h[i + ncopy * i] += *nradd;
            }
#endif

            t2 = seconds();
            *tnewtonOther += t2 - t1;
            t1 = t2;

            /* Solve the linear system via Cholesky factorization. */

            uplo = 'L';
            dposv_(&uplo, n, &one, h, n, &g[1], n, &info);

            t2 = seconds();
            *tnewtonCholesky += t2 - t1;
            if (nr_debug) {
               if (mytaskid == 0) {
                  printf("dposv time = %10.2f\n\n", t2 - t1);
                  fflush(stdout);
               }
            }
            t1 = t2;

            if (info < 0) {
               if (mytaskid == 0) {
                  printf("dposv returns %d; giving up...\n\n", info);
                  fflush(stdout);
               }
               ret_val = -2;
               goto cleanup;
            } else if (info > 0) {

#ifdef DSYEV

               /*
                *  Add a constant to the diagonal of the Hessian, ostensibly to
                *  make it positive definite..
                */

               if (mytaskid == 0) {
                  printf("dposv returns %d; attempting dsyev...\n\n",
                         info);
                  fflush(stdout);
               }

               /* Try modifying the Hessian with eigenvalues from dsyev_. */

               h2 = vector(0, ncopy * ncopy);   /* copy of hessian */

#pragma omp parallel for

               for (i = 0; i < ncopy; i++) {
                  in = ncopy * i;
                  for (j = 0; j < ncopy; j++) {
                     h2[in + j] = h[in + j];
                  }
               }

               roots = vector(0, ncopy);        /* for eigenvalues */

               /* Call dsyev_ using d as the work array. */

               uplo = 'L';
               jobz = 'N';
               lwork = 6 * ncopy;
               dsyev_(&jobz, &uplo, n, h2, n, roots, &d[1], &lwork, &info);

               t2 = seconds();
               *tnewtonDSYEV += t2 - t1;
               if (nr_debug) {
                  if (mytaskid == 0) {
                     printf("dsyev time = %10.2f\n\n", t2 - t1);
                     fflush(stdout);
                  }
               }
               t1 = t2;

               if (info) {
                  if (mytaskid == 0) {
                     printf("dspev returns %d; giving up...\n\n", info);
                     fflush(stdout);
                  }
                  ret_val = -3;
                  goto cleanup;
               }

               *nradd = 0.01 - roots[0];
               *nradd = *nradd < 0.0 ? 0.0 : *nradd;
               if (mytaskid == 0) {
                  printf("finding nradd: %10.5f %10.5f %10.5f\n",
                         roots[0], roots[1], *nradd);
                  fflush(stdout);
               }
               free_vector(roots, 0, ncopy);
               free_vector(h2, 0, ncopy * ncopy);

               t2 = seconds();
               if (nr_debug) {
                  if (mytaskid == 0) {
                     printf("dsyev time = %10.2f\n\n", t2 - t1);
                     fflush(stdout);
                  }
               }
               *tnewtonDSYEV += t2 - t1;
               t1 = t2;

               /* Modify the diagonal elements of the Hessian. */

               for (i = 0; i < ncopy; i++) {
                  h[i + ncopy * i] += *nradd;
               }

               t2 = seconds();
               *tnewtonOther += t2 - t1;
               t1 = t2;

               /* Another attempt at LAPACK solution to the NR equations. */

               uplo = 'L';
               dposv_(&uplo, n, &one, h, n, &g[1], n, &info);

               t2 = seconds();
               if (nr_debug) {
                  if (mytaskid == 0) {
                     printf("dposv time = %10.2f\n\n", t2 - t1);
                     fflush(stdout);
                  }
               }

               if (info && mytaskid == 0) {
                  printf("dposv returns %d; giving up...\n\n", info);
                  fflush(stdout);
                  ret_val = -2;
                  goto cleanup;
               }
#endif                          /* ifdef DSYEV  */

            }

            /*
             * Only one process is on the grid, so copy g into grad
             * instead of executing pdgemr2d_ below.
             */

            for (i = 1; i <= ncopy; i++) {
               grad[i] = g[i];
            }
         }
      } else {

#ifdef SCALAPACK

         /*
          * Here is the ScaLAPACK code that will execute if gridim != 1,
          * i.e., for more than one process on the process grid.
          *
          * Redistribute the level shift matrix d into the matrix dee.
          * The function pdgemr2d_ will hang unless called from all
          * tasks.
          */

         pdgemr2d_(&mcopy, &six,
                   d, &one, &one, descD_1x1,
                   dee, &one, &one, descD_PxQ, &context_Nx1);


         t2 = seconds();
         if (nr_debug) {
            if (mytaskid == 0) {
               printf("pdgemr2d time = %10.2f\n\n", t2 - t1);
               fflush(stdout);
            }
         }

         /* Multiply the Hessian by the level shift matrix dee. */

         if (context_PxQ >= 0) {
            transa = 'N';
            transb = 'T';
            pdgemm_(&transa, &transb, &mcopy, &mcopy, &six, &scale,
                    &dee[1], &one, &one, descD_PxQ,
                    &dee[1], &one, &one, descD_PxQ, &dblone,
                    h, &one, &one, descH_PxQ);
         }

         t2 = seconds();
	 *tnewtonOther += t2 - t1;
         if (nr_debug) {
            if (mytaskid == 0) {
               printf("pdgemm time = %10.2f\n\n", t2 - t1);
               fflush(stdout);
            }
         }
         t1 = t2;

         /* Solve the linear system via Cholesky factorization. */

         if (context_PxQ >= 0) {
            uplo = 'L';
            pdposv_(&uplo, n, &one,
                    h, &one, &one, descH_PxQ,
                    &g[1], &one, &one, descG_PxQ, &info);

            t2 = seconds();
            if (nr_debug) {
               if (mytaskid == 0) {
                  printf("pdposv time = %10.2f\n\n", t2 - t1);
                  fflush(stdout);
               }
            }
            *tnewtonCholesky += t2 - t1;
            t1 = t2;
         }

         if (context_PxQ >= 0 && info < 0) {
            if (mytaskid == 0) {
               printf("pdposv returns %d; giving up...\n\n", info);
               fflush(stdout);
            }
            ret_val = -2;
            goto cleanup;
         } else if (context_PxQ >= 0 && info > 0) {

#ifdef DSYEV

            /*
             *  Add a constant to the diagonal of the Hessian, ostensibly to
             *  make it positive definite..
             */

            if (mytaskid == 0) {
               printf("pdposv returns %d; attempting pdsyev...\n\n", info);
               fflush(stdout);
            }

            if (context_PxQ >= 0) {

               /* Calculate minimum size for the work array (from pdsyev_ source). */

               blacs_gridinfo_(&context_PxQ, &nprow, &npcol, &myrow,
                               &mycol);
               nb = descH_PxQ[NB_];
               np = numroc_(&mcopy, &nb, &myrow, &zero, &nprow);
               sizesytrd =
                   (3 * nb > nb * (np + 1)) ? 3 * nb : nb * (np + 1);
               lwork = 5 * ncopy + sizesytrd + 1;

               /* Allocate some temporary arrays. */

               h2 = vector(0, sizeH_PxQ);       /* copy of hessian */
               work = vector(0, lwork); /* work array (change size?) */
               roots = vector(0, ncopy);        /* for eigenvalues */

               for (i = 0; i < sizeH_PxQ; i++) {
                  h2[i] = h[i];
               }
            }

            /* Call pdsyev_ using d as the work array. */

            if (context_PxQ >= 0) {
               uplo = 'L';
               jobz = 'N';
               pdsyev_(&jobz, &uplo, n,
                       h2, &one, &one, descH_PxQ, roots,
                       h2, &one, &one, descH_PxQ, work, &lwork, &info);

               t2 = seconds();
               if (nr_debug) {
                  if (mytaskid == 0) {
                     printf("pdsyev time = %10.2f\n\n", t2 - t1);
                     fflush(stdout);
                  }
               }
               *tnewtonDSYEV += t2 - t1;
               t1 = t2;

               if (info) {
                  if (mytaskid == 0) {
                     printf("pdsyev returns %d; giving up...\n\n", info);
                     fflush(stdout);
                  }
                  ret_val = -3;
                  goto cleanup;
               }

               /* The roots vector is valid only within context_PxQ. */

               *nradd = 0.01 - roots[0];
               *nradd = *nradd < 0.0 ? 0.0 : *nradd;
               if (mytaskid == 0) {
                  printf("finding nradd: %10.5f %10.5f %10.5f\n",
                         roots[0], roots[1], *nradd);
                  fflush(stdout);
               }

               free_vector(h2, 0, ncopy * ncopy);
               free_vector(work, 0, lwork);
               free_vector(roots, 0, ncopy);
            }

            /*
             * Get the processor row and column in context_PxQ, as well as
             * the local leading dimension of the Hessian, and  modify the
             * diagonal elements of the Hessian.  If the task is not active
             * in context_PxQ, the myroc function returns 0.
             */

            if (context_PxQ >= 0) {
               for (i = 0; i < ncopy; i++) {
                  ptr = ptr2d(h, descH_PxQ, i, i);
                  if (ptr != NULL)
                     *ptr += *nradd;
               }
            }

            t2 = seconds();
            *tnewtonOther += t2 - t1;
            t1 = t2;

            /* Another attempt at ScaLAPACK solution to the NR equations. */

            if (context_PxQ >= 0) {
               uplo = 'L';
               pdposv_(&uplo, n, &one,
                       h, &one, &one, descH_PxQ,
                       &g[1], &one, &one, descG_PxQ, &info);
            }

            t2 = seconds();
            *tnewtonCholesky += t2 - t1;
            if (nr_debug) {
               if (mytaskid == 0) {
                  printf("pdposv time = %10.2f\n\n", t2 - t1);
                  fflush(stdout);
               }
            }
            t1 = t2;

            if (context_PxQ >= 0 && info) {
               if (mytaskid == 0) {
                  printf("pdposv returns %d; giving up...\n\n", info);
                  fflush(stdout);
               }
               ret_val = -2;
               goto cleanup;
            }
#endif                          /* ifdef DSYEV */

         }

         /*
          * More than one process is active on the grid, so gather the
          * distributed g vectors to the non-distributed grad vector.
          */

         pdgemr2d_(&mcopy, &one,
                   &g[1], &one, &one, descG_PxQ,
                   &grad[1], &one, &one, descG_1x1, &context_Nx1);

         t2 = seconds();
	 *tnewtonOther += t2 - t1;
         if (nr_debug) {
            if (mytaskid == 0) {
               printf("pdgemr2d time = %10.2f\n\n", t2 - t1);
               fflush(stdout);
            }
         }
         t1 = t2;

#endif                          /* ifdef SCALAPACK */

      }

      /*
       * Each MPI task requires its own copy of the grad vector, which might
       * be accomplished by broadcasting grad from the process that is active
       * on the 1x1 grid to all other tasks, including those tasks that are
       * not active on the PxQ grid.  However, the MPI_Bcast function requires
       * that the broadcasting task be specified.  Because I am unsure that in
       * all implementations of ScaLAPACK the process that is active on the
       * 1x1 grid is task 0, I will use the MPI_Allreduce function instead of
       * MPI_Bcast.  MPI_Allreduce does not require specification of the
       * broadcasting task.
       */

#ifdef SCALAPACK

      if (numtasks > 1) {

         /*
          * Zero the grad vector for all processes except the process that
          * is active on the 1x1 grid then call MPI_Allreduce, but only if
          * more than one process is executing.
          */

         if (context_1x1 < 0) {
            for (i = 1; i <= ncopy; i++) {
               grad[i] = 0.0;
            }
         }

         t2 = seconds();
	 *tnewtonOther += t2 - t1;
         t1 = t2;

         /*
          * Because the grad vector is nonzero for one process only,
          * MPI_Allreduce accomplishes a broadcast.
          */

         ierror = MPI_Allreduce(&grad[1], &reductarr[1], mcopy,
                                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         if (ierror != MPI_SUCCESS) {
            printf("Error in newton grad reduction, error = %d  "
                   "mytaskid = %d\n", ierror, mytaskid);
         }

         t2 = seconds();
	 *tnewtonOther += t2 - t1;
         if (nr_debug) {
            if (mytaskid == 0) {
               printf("MPI_Allreduce time = %10.2f\n\n", t2 - t1);
               fflush(stdout);
            }
         }
         t1 = t2;

         for (i = 0; i <= ncopy; i++) {
            grad[i] = reductarr[i];
         }
      }
#endif

      /* update the coordinates:  */

      alphal = dnrm2_(n, &grad[1], &one) / sqrn;
      if (mytaskid == 0) {
         printf("rms of search direction: %10.7f\n", alphal);
         fflush(stdout);
      }
      if (alphal < 0.0001) {

         /* If we are pretty close to a stationary point, rely on NR:  */

         for (i = 1; i <= ncopy; i++)
            x[i] -= grad[i];

      } else {

         /*
          * Since the second derivatives are much more expensive than first
          * derivatives, here we try to get the most out of the NR procedure.
          * First, do a line search along the NR search direction:
          */

         for (i = 1; i <= ncopy; i++)
            grad[i] = -grad[i];
         alphal = 1.0;
         evlmax = 5;
         toll[0] = 0.00000001;
         toll[1] = 0.00000001;
         linemin(n, &x[1], &grad[1], &alphal, func1, &evlmax, toll);

#if 0
         /*
          * Next, do a few steps of CG minimization for all variables, 
          * which should do a little bit of fixing along directions
          * perpendicular to the NR search direction.
          */

         rmsg = 0.00001;
         dfpred = 0.01;         /* universal value? get from linemin?  */
         maxit = 25;
         conjgrad(&x[1], n, &fret, func1, &rmsg, &dfpred, &maxit);
         if (mytaskid == 0) {
            printf("After CG minimization, energy = %20.10f\n", fret);
            fflush(stdout);
         }
#endif
      }
   } while (niter < *maxiter);

   /* Prepare return value. */

   if (niter >= *maxiter)
      ret_val = -1;
   else
      ret_val = niter;

   /*
    * Free the dynamic vectors and matrices, including those of egb2
    * which is called to free static vectors when *func is called with
    * an iteration count equal to -3.
    */

 cleanup:
   niter = -3;
   (*func2) (&x[1], &g[1], h, &m[1],
             &grad[1], descG_PxQ, descG_1x1, descH_PxQ,
             &context_PxQ, &context_1x1, &context_Nx1,
             &gridim, &natom, &niter, name);

   if (m != NULL)
      free_vector(m, 1, ncopy);

#ifndef SCALAPACK

   if (d != NULL)
      free_vector(d, 1, 6 * ncopy);
   if (g != NULL)
      free_vector(g, 1, ncopy);
   if (grad != NULL)
      free_vector(grad, 1, ncopy);
   if (h != NULL)
      free_vector(h, 0, ncopy * ncopy);

#else

   if (reductarr != NULL)
      free_vector(reductarr, 1, ncopy);
   if (context_PxQ >= 0) {
      if (dee != NULL)
         free_vector(dee, 1, sizeD_PxQ);
      if (g != NULL)
         free_vector(g, 1, sizeG_PxQ);
      if (h != NULL)
         free_vector(h, 0, sizeH_PxQ);
   }
   if (context_1x1 >= 0) {
      if (d != NULL)
         free_vector(d, 1, sizeD_1x1);
      if (grad != NULL)
         free_vector(grad, 1, sizeG_1x1);
   } else {
      if (grad != NULL)
         free_vector(grad, 1, ncopy);
   }

   /* Exit the process grid only for tasks that are on the process grid. */

   if (context_PxQ >= 0) {
      blacs_gridexit_(&context_PxQ);
   }
   if (context_1x1 >= 0) {
      blacs_gridexit_(&context_1x1);
   }

   if (context_Nx1 >= 0) {
      blacs_gridexit_(&context_Nx1);
   }
#endif                          /* ifndef SCALAPACK */

   t2 = seconds();
   *tnewtonOther += t2 - t1;
   *tnewton += t2 - tnewton1;
   return ret_val;
}
