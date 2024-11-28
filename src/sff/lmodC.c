/***************************************************************************/
/*      LMOD: Written by Istvan Kolossvary     time stamp: 06/03/2016      */
/***************************************************************************/

/*  Note: there are two "public" routines here:
 *  lmodC():  a fairly low-level interface (currently used by sander), which
 *            requires its own driver that understands the reverse
 *            communication needed; this name may get mangled for Fortran.
 *  lmod():   a higher-level, NAB-like interface, which takes the system
 *            to be simulated as arguments, and handles the reverse
 *            communication internally.
 */

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>           /* For gettimeofday */

#if defined(MPI) || defined(SCALAPACK)
#include <mpi.h>
#endif

#ifdef SQM
#  define lmodC lmodc_
#else
#  include "sff.h"
#endif

#define DBL_EPSILON 2.2204460492503131e-16

#define SQR(a) ((a)*(a))
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))


#define DONE            0
#define CALCGRAD_NEWNBL 1
#define CALCBOTH_NEWNBL 2
#define CALCGRAD_OLDNBL 3
#define CALCBOTH_OLDNBL 4
#define MINIMIZE        5
#define RELAX           6

#define PARAMS_ERROR    -1
#define ILLEGAL_STATUS  -2
#define MALLOC_ERROR    -3
#define FILE_ERROR      -4
#define DIVZERO_ERROR   -5
#define ARPACK_ERROR    -6
#define CONFLIB_ERROR   -7


#define ZERO            0.0
#define ONE             1.0
#define BIG             1e10
#define TINY            1e-20
#define YES             1
#define NO              0
#define TRUE            1
#define FALSE           0
#define PI              3.141592653589793
#define RAD2DEG         (180./PI)
#define DEG2RAD         (PI/180.)


#define METROPOLIS      1
#define TOTAL_QUENCHING 2
#define QUICK_QUENCHING 3
#define MAX_MC_TRY      10000
#define MAX_ZIGZAG_ITER 20
#define FULLY_MINIMIZED 10.0    /* 0.5 */
#define CONFLIB_UPDATE  10

#define HESSVEC_FORWARD
#undef HESSVEC_CENTRAL

/* ARPACK: */
#define ARPK_ITER       1000000
#define ARPK_TOL        0.0001  /* 0.0 = machine precision */
#define ARPK_SM         1
#define ARPK_SA         2
#define ARPK_LM         3
#define ARPK_LA         4
#define ARPK_BE         5


static int DEBUG_ARPACK = NO;
static int PRINT_LMOD = NO;
static int DEBUG_LMOD = NO;


static struct archive {
   double energy;
   double rad;
   int mult;
   double *xyz;
   struct archive *next;
} *conflib_archive = NULL;

#ifdef SQM
static FILE *nabout;
#endif

/*
    Defined elsewhere: (../arpack/dsarpack.f, parm.c)
*/
extern int dsarpack_();
extern int get_mytaskid();      /* for MPI */
extern int setseed();           /* for use of NAB's own    */
extern int rseed();             /* random number genarator */
extern double rand2();          /* see rand2.c             */


/*
    Private to liblmod.c:
static void absort();
static void archive_structure();
static void arpack();
static void calc_centroid();
static double calc_rad();
static void calc_rot_matrix();
static void destroy_archive();
static void diag();
static void hessdump();
static void load_archive();
static void my_free();
static void *my_malloc();
static int read_archive();
static int read_and_sort_archive();
static double rmsfit();
static void restart_lmod();
static void rot_ligand();
static void separate_close_pairs();
static void trans_ligand();
static int update_archive();
static void write_archive();
*/

static void absort(int n, double *crr, int *brr, int *error_flag);
static void archive_structure(int ndim, double *xyz, double energy, double rad, int *error_flag);
static void arpack(int *ndim, int *nof_requested_modes, int *nof_computed_modes,
                   int *arpk_dim, double *eigvals, double *eigvecs, int *spectrum,
                   int *want_eigvecs, double *xyz, double *grad, int *return_flag, int *label);
static double calc_rad(int natm, double *xyz);
static void calc_centroid(double *xyz, int start, int end, double *cent_x,
                          double *cent_y, double *cent_z);
static void calc_rot_matrix(double alpha, double rotax_x, double rotax_y,
                            double rotax_z, double rotmat[3][3]);
static void destroy_archive();
static void diag(int n, double aa[3][3], double *d, double x[3][3], int *error_flag);
/* static void hessdump(int *ndim, double *xyz, double *grad, int * *return_flag, int *label); */
static void load_archive(int *nconf, int ndim, double *xyz_array, int *error_flag);
static void my_free(void *poi);
static void *my_malloc(void *(*malloc_method) (size_t), const char *s,
                       size_t nmemb, size_t size, int *error_flag);
/* static int read_archive(char *fname, int ndim, int *error_flag); */
static int read_and_sort_archive(char *fname, int ndim, int *error_flag);
static void restart_lmod(int ndim, int nmax, double *xyz);
static double rmsfit(int ndim, double *ref, double *xyz, char *which, int *error_flag);
static void separate_close_pairs(int do_all_pairs, int do_ligs_only, int ndim,
                     double *xyz, double *grad, int nligs, int *lig_start, int *lig_end);
static void trans_ligand(double *xyz, int start, int end, double dx, double dy, double dz);
static int update_archive(int ndim_ext, int ndim_local, double *energy_window, double *separation_rms, int *error_flag);
static void write_archive(char *fname, int ndim, int *error_flag);


/*****
        Dynamic memory allocation:
*****/
static void *my_malloc(void *(*malloc_method) (size_t), const char *s,
                       size_t nmemb, size_t size, int *error_flag)
/*
    malloc_method is a pointer to a particular malloc function.
*/
{
   void *poi;
   *error_flag = FALSE;
   if ((poi = (void *) (*malloc_method) (nmemb * size)) == NULL) {
      perror(s);
      fflush(stderr);
      *error_flag = MALLOC_ERROR;
      return NULL;
   }
   memset(poi, 0, nmemb * size);        /* clear memory */
 /*printf("\n Allocated %10d bytes at %p for %s\n",(int)(nmemb*size),poi,s+10);
   fflush(stdout); */
   return poi;
}


static void my_free(void *poi)
{
   if (poi != NULL)
      free(poi);
 /*printf("\n Deallocated                   %p\n",poi); */
   fflush(stdout); 
}


/*****
        Linked list archiving:
*****/
static void
archive_structure(int ndim, double *xyz, double energy, double rad,
                  int *error_flag)
{
   int last;
   struct archive *curr, *prev;
   if (conflib_archive == NULL) {       /* open list */
      conflib_archive = (struct archive *)
          my_malloc(malloc,
                    "\nERROR in archive_structure/my_malloc(struct archive *conflib_archive)",
                    1, sizeof(struct archive), error_flag);
      if (*error_flag) {
         conflib_archive = NULL;
         return;
      }
      conflib_archive->energy = energy;
      conflib_archive->rad = rad;
      conflib_archive->mult = 1;
      conflib_archive->xyz = (double *)
          my_malloc(malloc,
                    "\nERROR in archive_structure/my_malloc(double *conflib_archive->xyz)",
                    ndim, sizeof(double), error_flag);
      if (*error_flag) {
         my_free(conflib_archive);
         conflib_archive = NULL;
         return;
      }
      memcpy(conflib_archive->xyz, xyz, ndim * sizeof(double));
      conflib_archive->next = NULL;
   } else {                     /* sort existing list in ascending energetic order */
      for (curr = conflib_archive, prev = NULL, last = YES; curr != NULL;
           prev = curr, curr = curr->next) {
         if (energy <= curr->energy) {  /* found slot right before which new structure is to be inserted */
            if (prev == NULL) { /* new global minimum goes to the top of the list */
               conflib_archive = (struct archive *)
                   my_malloc(malloc,
                             "\nERROR in archive_structure/my_malloc(struct archive *conflib_archive)",
                             1, sizeof(struct archive), error_flag);
               if (*error_flag) {
                  conflib_archive = curr;
                  destroy_archive();
                  return;
               }
               conflib_archive->energy = energy;
               conflib_archive->rad = rad;
               conflib_archive->mult = 1;
               conflib_archive->xyz = (double *)
                   my_malloc(malloc,
                             "\nERROR in archive_structure/my_malloc(double *conflib_archive->xyz)",
                             ndim, sizeof(double), error_flag);
               if (*error_flag) {
                  my_free(conflib_archive);
                  conflib_archive = curr;
                  destroy_archive();
                  return;
               }
               memcpy(conflib_archive->xyz, xyz, ndim * sizeof(double));
               conflib_archive->next = curr;
            } else {
               prev->next = (struct archive *)
                   my_malloc(malloc,
                             "\nERROR in archive_structure/my_malloc(struct archive *prev->next)",
                             1, sizeof(struct archive), error_flag);
               if (*error_flag) {
                  prev->next = curr;
                  destroy_archive();
                  return;
               }
               prev->next->energy = energy;
               prev->next->rad = rad;
               prev->next->mult = 1;
               prev->next->xyz = (double *)
                   my_malloc(malloc,
                             "\nERROR in archive_structure/my_malloc(double *prev->next->xyz)",
                             ndim, sizeof(double), error_flag);
               if (*error_flag) {
                  my_free(prev->next);
                  prev->next = curr;
                  destroy_archive();
                  return;
               }
               memcpy(prev->next->xyz, xyz, ndim * sizeof(double));
               prev->next->next = curr;
            }
            last = NO;
            break;
         }
      }
      if (last) {               /* add currently highest energy structure to the bottom of the list */
         prev->next = (struct archive *)
             my_malloc(malloc,
                       "\nERROR in archive_structure/my_malloc(struct archive *prev->next)",
                       1, sizeof(struct archive), error_flag);
         if (*error_flag) {
            prev->next = NULL;
            destroy_archive();
            return;
         }
         prev->next->energy = energy;
         prev->next->rad = rad;
         prev->next->mult = 1;
         prev->next->xyz = (double *)
             my_malloc(malloc,
                       "\nERROR in archive_structure/my_malloc(double *prev->next->xyz)",
                       ndim, sizeof(double), error_flag);
         if (*error_flag) {
            my_free(prev->next);
            prev->next = NULL;
            destroy_archive();
            return;
         }
         memcpy(prev->next->xyz, xyz, ndim * sizeof(double));
         prev->next->next = NULL;
      }
   }
}


static int update_archive(int ndim_ext, int ndim_local, double *energy_window, double *separation_rms, int *error_flag)
/* Archive 'conflib_archive' must be sorted by increasing energy! */
{
   int ndat, first_here;
   struct archive *curr, *curr2, *next_curr, *hold, *last;
   double *tmp, rms;
   tmp = (double *)
       my_malloc(malloc,
                 "\nERROR in update_archive/my_malloc(double *tmp)", ndim_ext,
                 sizeof(double), error_flag);
   if (*error_flag)
      return(0);

   for (curr = conflib_archive; curr != NULL; last = curr, curr = curr->next) { /* walk through linked list */
      if (fabs(curr->energy - conflib_archive->energy) > *energy_window) {      /* cut off bottom of list   */
         last->next = NULL;     /* close top part of list */
         do {
            hold = curr->next;
            my_free(curr->xyz);
            my_free(curr);
            curr = hold;
         } while (curr != NULL);
         break;
      }
   }

   for (curr = conflib_archive; curr->next != NULL; curr = next_curr) {

      for (first_here = TRUE, curr2 = curr; curr2->next != NULL;) {
      /*if(get_mytaskid()==0)
        printf("curr=%p  curr->xyz=%p  curr->energy=%f  curr2->next=%p  curr2->next->xyz=%p  curr2->next->energy=%f\n",
        curr,curr->xyz,curr->energy,curr2->next,curr2->next->xyz,curr2->next->energy);*/
         memcpy(tmp, curr2->next->xyz, ndim_ext * sizeof(double));
         rms = rmsfit(ndim_ext, curr->xyz, tmp, "RMS", error_flag);
         rms *= sqrt((double)ndim_ext / (double)ndim_local);  /* frozen atoms' contribution is zero,       */
         /*if(get_mytaskid()==0)printf("rms=%f\n",rms);*/     /* denominator in rms shouldn't include them */
         if (*error_flag) {
            my_free(tmp);
            destroy_archive();
            return(0);
         }
         if (rms <= *separation_rms) {              /* Confs identical? */
            curr->mult += curr2->next->mult;        /* Keep only curr   */
            hold = curr2->next->next;
            my_free(curr2->next->xyz);
            my_free(curr2->next);
            curr2->next = hold;
            if (curr2->next == NULL) {/* deleted last entry in conf list*/
               break;
            }
         }
         else {
            if (first_here) {
               next_curr = curr2->next;
               first_here = FALSE;
            }
            curr2 = curr2->next;
         }
      }
      if (first_here == TRUE) { /* no more unique conf found */
         curr->next = NULL;
         break;
      }
   }
   for (ndat=0, curr = conflib_archive; curr != NULL; curr = curr->next) {ndat++;}

   my_free(tmp);
   return ndat;
}


static void write_archive(char *fname, int ndim, int *error_flag)
{
   struct archive *walk;
   FILE *file;
   if ((file = fopen(fname, "wb")) == NULL) {
      fprintf(stderr, "\nERROR in write_archive(): Cannot open %s.\n",
              fname);
      fflush(stderr);
      *error_flag = FILE_ERROR;
      return;
   }
   for (walk = conflib_archive; walk != NULL; walk = walk->next) {
      fwrite(&walk->energy, sizeof(walk->energy), 1, file);
      fwrite(&walk->rad, sizeof(walk->rad), 1, file);
      fwrite(&walk->mult, sizeof(walk->mult), 1, file);
      fwrite(walk->xyz, sizeof(double), ndim, file);
   }
   fclose(file);
}


static void destroy_archive()
{
   struct archive *walk, *hold;
   for (walk = conflib_archive; walk != NULL; walk = hold) {
      hold = walk->next;
      my_free(walk->xyz);
      my_free(walk);
   }
   conflib_archive = NULL;      /* label archive as empty */
}

#if 0
static int read_archive(char *fname, int ndim, int *error_flag)
{
   int ndat;
   struct archive *curr, *prev;
   FILE *file;
   if (conflib_archive != NULL)
      destroy_archive();
   conflib_archive = (struct archive *)
       my_malloc(malloc,
                 "\nERROR in read_archive/my_malloc(struct archive *conflib_archive)",
                 1, sizeof(struct archive), error_flag);
   if (*error_flag) {
      conflib_archive = NULL;
      return(0);
   }
   if ((file = fopen(fname, "rb")) == NULL) {
      fprintf(stderr, "\nERROR in read_archive(): Cannot open %s.\n",
              fname);
      fflush(stderr);
      *error_flag = FILE_ERROR;
      return(0);
   }
   for (ndat = 0, prev = NULL, curr = conflib_archive;;
        prev = curr, curr = curr->next, ndat++) {
      if (fread(&curr->energy, sizeof(curr->energy), 1, file) < 1)
         break;
      if (fread(&curr->rad, sizeof(curr->rad), 1, file) < 1)
         break;
      if (fread(&curr->mult, sizeof(curr->mult), 1, file) < 1)
         break;
      curr->xyz = (double *)
          my_malloc(malloc,
                    "\nERROR in read_archive/my_malloc(double *curr->xyz)",
                    ndim, sizeof(double), error_flag);
      if (*error_flag) {
         if (prev != NULL) {
            my_free(curr);      /* free extra slot */
            prev->next = NULL;  /* close archive   */
            destroy_archive();
         } else {
            my_free(conflib_archive);
            conflib_archive = NULL;
         }
         return(0);
      }
      if (fread(curr->xyz, sizeof(double), ndim, file) < ndim)
         break;
      curr->next = (struct archive *)
          my_malloc(malloc,
                    "\nERROR in read_archive/my_malloc(struct archive *curr->next)",
                    1, sizeof(struct archive), error_flag);
      if (*error_flag) {
         curr->next = NULL;     /* close archive   */
         destroy_archive();
         return(0);
      }
   }
   if (ferror(file)) {
      fprintf(stderr, "\nERROR in read_archive(): Input error in %s.\n",
              fname);
      fflush(stderr);
      *error_flag = FILE_ERROR;
      return(0);
   }
   my_free(curr);               /* free extra slot */
   prev->next = NULL;           /* close archive   */
   fclose(file);
   return ndat;
}
#endif

static int read_and_sort_archive(char *fname, int ndim, int *error_flag) /* ndim = ndim_ext = 3 * natm_ext */
{
   FILE *file;
   int ndat, mult;
   double energy, rad, *xyz;

   if ((file = fopen(fname, "rb")) == NULL) {
      fprintf(stderr, "\nERROR in read_and_sort_archive(): Cannot open %s.\n",
              fname);
      fflush(stderr);
      *error_flag = FILE_ERROR;
      return(0);
   }
   xyz = (double *)
       my_malloc(malloc,
                 "\nERROR in read_and_sort_archive/my_malloc(double *xyz)",
                 ndim, sizeof(double), error_flag);
   if (*error_flag) {
      return(0);
   }

   for (ndat = 0;; ndat++) {

      if (fread(&energy, sizeof(energy), 1, file) < 1)
         break;
      if (fread(&rad, sizeof(rad), 1, file) < 1)
         break;
      if (fread(&mult, sizeof(mult), 1, file) < 1)
         break;
      if (fread(xyz, sizeof(double), ndim, file) < ndim)
         break;

      /* add conf to conflib */
      archive_structure(ndim, xyz, energy, rad, error_flag); /* sorting happens here */
      if (*error_flag) {
         return(0);
      }
   }
   if (ferror(file)) {
      fprintf(stderr, "\nERROR in read_archive(): Input error in %s.\n",
              fname);
      fflush(stderr);
      *error_flag = FILE_ERROR;
      return(0);
   }

   fclose(file);
   my_free(xyz);
   return ndat;
}


static void
load_archive(int *nconf, int ndim, double *xyz_array, int *error_flag)
{
   int n;
   struct archive *walk;
   if (*nconf <= 0)
      return;
   if (conflib_archive == NULL) {
      fprintf(stderr,
              "\nERROR in load_archive(): conflib database is empty.\n");
      fflush(stderr);
      *error_flag = CONFLIB_ERROR;
      return;
   }
   for (n = 0, walk = conflib_archive; walk != NULL; walk = walk->next) {
      memcpy(xyz_array + n * ndim, walk->xyz, ndim * sizeof(double));
      if (++n == *nconf)
         break;
   }
   if (n < *nconf)
      *nconf = n;
}


static void restart_lmod(int ndim, int nmax, double *xyz)
{
   int n, nrand;
   struct archive *walk;
   if (conflib_archive == NULL) {       /* do nothing */
      return;
   }
   nrand = 1 + nmax * rand2();
   for (n = 0, walk = conflib_archive; walk != NULL; walk = walk->next) {
      memcpy(xyz, walk->xyz, ndim * sizeof(double));
      if (++n == nrand)
         break;
   }
}


/*****
        Miscellaneous:
*****/
static double calc_rad(int natm, double *xyz)
/* Calculate radius of giration. */
{
   int i;
   double avg = ZERO, var = ZERO, dist;
   double center[3];
   center[0] = center[1] = center[2] = ZERO;
   for (i = 0; i < natm; i++) {
      center[0] += xyz[i * 3];
      center[1] += xyz[i * 3 + 1];
      center[2] += xyz[i * 3 + 2];
   }
   center[0] /= natm;
   center[1] /= natm;
   center[2] /= natm;
   for (i = 0; i < natm; i++) {
      dist = SQR(xyz[i * 3] - center[0])
          + SQR(xyz[i * 3 + 1] - center[1])
          + SQR(xyz[i * 3 + 2] - center[2]);
      var += dist;
      dist = sqrt(dist);
      avg += dist;
   }
   avg /= natm;
   var = var / natm - SQR(avg);
   return sqrt(MAX(ZERO, var));
}


static void
diag(int n, double aa[3][3], double *d, double x[3][3], int *error_flag)
{
   int i, j, k, l, np, n2, j1;
   double *e, *a, ha, ep, tl, h, g, f, s, h1, ze, b, el, r, p, dx, c, ei;
   e = (double *) my_malloc(malloc,
                            "\nERROR in diag/my_malloc(double *e)", n,
                            sizeof(double), error_flag);
   if (*error_flag)
      return;
   a = (double *) my_malloc(malloc,
                            "\nERROR in diag/my_malloc(double *a)", n * n,
                            sizeof(double), error_flag);
   if (*error_flag) {
      my_free(e);
      return;
   }
   for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         a[i * n + j] = aa[i][j];
   tl = 1.e-29;
   ha = .5;
   np = n - 1;
   n2 = np;
   ep = 1.e-08;
   ze = ZERO;
   for (i = 0; i <= np; i++)
      for (j = i + 1; j <= np; j++)
         a[i * n + j] = ZERO;
   if (np > 1)
      for (i = np; i >= 2; i--) {
         l = i - 2;
         h = 0;
         g = a[i * n + i - 1];
         for (k = 0; k <= l; k++) {
            f = a[i * n + k];
            h += f * f;
         }
         if (h > 0) {
            s = h + g * g;
            if (s >= tl) {
               l++;
               f = g;
               g = sqrt(s);
               g = (g + s / g) * ha;
               if (f > 0.)
                  g = -g;
               h = s - f * g;
               a[i * n + i - 1] = f - g;
               f = 0.;
               h1 = 1. / h;
               for (j = 0; j <= l; j++) {
                  a[j * n + i] = a[i * n + j] * h1;
                  s = 0.;
                  for (k = 0; k <= j; k++)
                     s += a[j * n + k] * a[i * n + k];
                  j1 = j + 1;
                  if (j1 <= l)
                     for (k = j1; k <= l; k++)
                        s += a[k * n + j] * a[i * n + k];
                  e[j] = s * h1;
                  f += s * a[j * n + i];
               }
               h1 *= f * ha;
               for (j = 0; j <= l; j++) {
                  f = a[i * n + j];
                  s = e[j];
                  s -= h1 * f;
                  e[j] = s;
                  for (k = 0; k <= j; k++)
                     a[j * n + k] -= f * e[k] + a[i * n + k] * s;
               }
            }
         }
         d[i] = h;
         e[i - 1] = g;
      }
   e[0] = a[n];
   d[0] = a[0];
   a[0] = 1.;
   for (i = 1; i <= np; i++) {
      l = i - 1;
      if (d[i] > 0.) {
         for (j = 0; j <= l; j++) {
            s = 0.0;
            for (k = 0; k <= l; k++)
               s += a[k * n + j] * a[i * n + k];
            for (k = 0; k <= l; k++)
               a[k * n + j] -= s * a[k * n + i];
         }
      }
      d[i] = a[i * (n + 1)];
      a[i * (n + 1)] = 1.;
      for (j = 0; j <= l; j++)
         a[i * n + j] = a[j * n + i] = ze;
   }
   b = 0.0;
   f = 0.0;
   e[np] = 0.0;
   for (l = 0; l <= np; l++) {
      el = e[l];
      r = fabs(el);
      h = ep * (fabs(d[l]) + r);
      if (h > b)
         b = h;
      if (r > b) {
         for (i = l + 1; i <= np; i++)
            if (fabs(e[i]) <= b) {
               j = i;
               i = np + 1;
            }
         do {
            p = (d[l + 1] - d[l]) * ha / el;
            dx = p * p + 1.;
            r = sqrt(dx);
            r = (r + dx / r) * ha;
            if (p < 0.0)
               p -= r;
            else
               p += r;
            h = d[l] - el / p;
            for (i = l; i <= np; i++)
               d[i] -= h;
            f += h;
            p = d[j];
            c = 1.;
            s = 0.0;
            j1 = j - 1;
            for (i = j1; i >= l; i--) {
               ei = e[i];
               g = c * ei;
               h = c * p;
               if (fabs(p) < fabs(ei)) {
                  c = p / ei;
                  dx = c * c + 1.;
                  r = sqrt(dx);
                  r = (r + dx / r) * ha;
                  e[i + 1] = s * ei * r;
                  s = 1. / r;
                  c /= r;
               } else {
                  c = ei / p;
                  dx = c * c + 1.;
                  r = sqrt(dx);
                  r = (r + dx / r) * ha;
                  e[i + 1] = s * p * r;
                  s = c / r;
                  c = 1. / r;
               }
               p = c * d[i] - s * g;
               d[i + 1] = h + s * (c * g + s * d[i]);
               for (k = 0; k <= np; k++) {
                  h = a[k * n + i + 1];
                  g = a[k * n + i];
                  a[k * n + i + 1] = g * s + h * c;
                  a[k * n + i] = g * c - h * s;
               }
            }
            el = s * p;
            e[l] = el;
            d[l] = c * p;
         } while (fabs(el) > b);
      }
      d[l] += f;
   }
   for (i = 0; i <= n2; i++) {
      k = i;
      p = d[i];
      for (j = i + 1; j <= np; j++)
         if (d[j] > p) {
            k = j;
            p = d[j];
         }
      if (k != i) {
         d[k] = d[i];
         d[i] = p;
         for (j = 0; j <= np; j++) {
            p = a[j * n + i];
            a[j * n + i] = a[j * n + k];
            a[j * n + k] = p;
         }
      }
   }
   for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         x[i][j] = a[i * n + j];
   my_free(e);
   my_free(a);
}


static double
rmsfit(int ndim, double *ref, double *xyz, char *which, int *error_flag)
/*
    Based on   W. Kabsch; Acta Cryst. (1976). A32, 922-923,
               W. Kabsch; Acta Cryst. (1978). A34, 827-828.
*/
{
   int i, j, n, natm;
   double rmatrx[3][3], rtrmat[3][3], d[3], z[3][3], b[3][3], cntr_ref[3],
       cntr_xyz[3], umat[3][3];
   double fac, xsq, ysq, zsq, dsq, dmax, savex, savey, savez, rms;
   double *save_ref, *save_xyz;  /* only compute RMS in place, do not apply transformation */
   natm = ndim / 3;
   save_ref = (double *)
       my_malloc(malloc,
                 "\nERROR in rmsfit/my_malloc(double *save_ref)", ndim,
                 sizeof(double), error_flag);
   save_xyz = (double *)
       my_malloc(malloc,
                 "\nERROR in rmsfit/my_malloc(double *save_xyz)", ndim,
                 sizeof(double), error_flag);
   if (*error_flag)
      return(0.0);
   memcpy(save_ref, ref, ndim * sizeof(double));
   memcpy(save_xyz, xyz, ndim * sizeof(double));
  for (i = 0; i < 3; i++)
      cntr_ref[i] = cntr_xyz[i] = ZERO;
   for (n = 0; n < natm; n++) {
      cntr_ref[0] += ref[n * 3];
      cntr_ref[1] += ref[n * 3 + 1];
      cntr_ref[2] += ref[n * 3 + 2];
      cntr_xyz[0] += xyz[n * 3];
      cntr_xyz[1] += xyz[n * 3 + 1];
      cntr_xyz[2] += xyz[n * 3 + 2];
   }
   cntr_ref[0] /= natm;
   cntr_ref[1] /= natm;
   cntr_ref[2] /= natm;
   cntr_xyz[0] /= natm;
   cntr_xyz[1] /= natm;
   cntr_xyz[2] /= natm;
   for (n = 0; n < natm; n++) {
      ref[n * 3] -= cntr_ref[0];
      ref[n * 3 + 1] -= cntr_ref[1];
      ref[n * 3 + 2] -= cntr_ref[2];
      xyz[n * 3] -= cntr_xyz[0];
      xyz[n * 3 + 1] -= cntr_xyz[1];
      xyz[n * 3 + 2] -= cntr_xyz[2];
   }
   for (j = 0; j < 3; j++)
      for (i = 0; i < 3; i++)
         rmatrx[i][j] = ZERO;
   for (n = 0; n < natm; n++) {
      rmatrx[0][0] += 0.1 * ref[n * 3] * xyz[n * 3];
      rmatrx[1][0] += 0.1 * ref[n * 3 + 1] * xyz[n * 3];
      rmatrx[2][0] += 0.1 * ref[n * 3 + 2] * xyz[n * 3];
      rmatrx[0][1] += 0.1 * ref[n * 3] * xyz[n * 3 + 1];
      rmatrx[1][1] += 0.1 * ref[n * 3 + 1] * xyz[n * 3 + 1];
      rmatrx[2][1] += 0.1 * ref[n * 3 + 2] * xyz[n * 3 + 1];
      rmatrx[0][2] += 0.1 * ref[n * 3] * xyz[n * 3 + 2];
      rmatrx[1][2] += 0.1 * ref[n * 3 + 1] * xyz[n * 3 + 2];
      rmatrx[2][2] += 0.1 * ref[n * 3 + 2] * xyz[n * 3 + 2];
   }
   for (j = 0; j < 3; j++)
      for (i = 0; i < 3; i++)
         rtrmat[j][i] = rmatrx[0][i] * rmatrx[0][j]
             + rmatrx[1][i] * rmatrx[1][j]
             + rmatrx[2][i] * rmatrx[2][j];
   diag(3, rtrmat, d, z, error_flag);
   if (*error_flag)
      return(0.0);
   z[0][2] = z[1][0] * z[2][1] - z[2][0] * z[1][1];
   z[1][2] = z[2][0] * z[0][1] - z[0][0] * z[2][1];
   z[2][2] = z[0][0] * z[1][1] - z[1][0] * z[0][1];
   for (j = 0; j < 2; j++) {
      for (i = 0; i < 3; i++)
         b[i][j] = rmatrx[i][0] * z[0][j]
             + rmatrx[i][1] * z[1][j]
             + rmatrx[i][2] * z[2][j];
      fac =
          sqrt(b[0][j] * b[0][j] + b[1][j] * b[1][j] + b[2][j] * b[2][j]);
      for (i = 0; i < 3; i++)
         b[i][j] /= MAX(fac, TINY);
   }
   b[0][2] = b[1][0] * b[2][1] - b[2][0] * b[1][1];
   b[1][2] = b[2][0] * b[0][1] - b[0][0] * b[2][1];
   b[2][2] = b[0][0] * b[1][1] - b[1][0] * b[0][1];
   for (j = 0; j < 3; j++)
      for (i = 0; i < 3; i++)
         umat[i][j] = b[i][0] * z[j][0]
             + b[i][1] * z[j][1]
             + b[i][2] * z[j][2];
   savex = cntr_xyz[0];
   savey = cntr_xyz[1];
   savez = cntr_xyz[2];
   for (n = 0, dsq = dmax = ZERO; n < natm; n++) {
      savex = xyz[n * 3];
      savey = xyz[n * 3 + 1];
      savez = xyz[n * 3 + 2];
      xyz[n * 3] = umat[0][0] * savex
          + umat[0][1] * savey + umat[0][2] * savez;
      xsq = SQR(xyz[n * 3] - ref[n * 3]);
      dsq += xsq;
      xyz[n * 3 + 1] = umat[1][0] * savex
          + umat[1][1] * savey + umat[1][2] * savez;
      ysq = SQR(xyz[n * 3 + 1] - ref[n * 3 + 1]);
      dsq += ysq;
      xyz[n * 3 + 2] = umat[2][0] * savex
          + umat[2][1] * savey + umat[2][2] * savez;
      zsq = SQR(xyz[n * 3 + 2] - ref[n * 3 + 2]);
      dsq += zsq;
      if ((xsq + ysq + zsq) > dmax)
         dmax = xsq + ysq + zsq;
   }
   /* move ref and xyz structures back to their original position: */
   memcpy(ref, save_ref, ndim * sizeof(double));
   memcpy(xyz, save_xyz, ndim * sizeof(double));
   my_free(save_ref);
   my_free(save_xyz);
   if (!strcmp(which, "RMS")) {
      rms = sqrt(dsq / natm);
      return rms;
   } else {
      dmax = sqrt(dmax);
      return dmax;
   }
}


/*****
        Trans/Rot ligand:
*****/
static void
trans_ligand(double *xyz, int start, int end, double dx, double dy,
             double dz)
{
   int i, x, y, z;
   for (i = start; i <= end; i++) {
      x = 3 * i;
      y = x + 1;
      z = y + 1;
      xyz[x] += dx;
      xyz[y] += dy;
      xyz[z] += dz;
   }
}


static void
calc_rot_matrix(double alpha, double rotax_x, double rotax_y,
                double rotax_z, double rotmat[3][3])
{
   double ca, sa;
   ca = cos(alpha);
   sa = sin(alpha);
   rotmat[0][0] = rotax_x * rotax_x + ca * (ONE - rotax_x * rotax_x);
   rotmat[0][1] = (ONE - ca) * rotax_x * rotax_y - sa * rotax_z;
   rotmat[0][2] = (ONE - ca) * rotax_x * rotax_z + sa * rotax_y;
   rotmat[1][0] = (ONE - ca) * rotax_y * rotax_x + sa * rotax_z;
   rotmat[1][1] = rotax_y * rotax_y + ca * (ONE - rotax_y * rotax_y);
   rotmat[1][2] = (ONE - ca) * rotax_y * rotax_z - sa * rotax_x;
   rotmat[2][0] = (ONE - ca) * rotax_z * rotax_x - sa * rotax_y;
   rotmat[2][1] = (ONE - ca) * rotax_z * rotax_y + sa * rotax_x;
   rotmat[2][2] = rotax_z * rotax_z + ca * (ONE - rotax_z * rotax_z);
}


static void
rot_ligand(double *xyz, int start, int end, double cent_x, double cent_y,
           double cent_z, double rotmat[3][3])
{
   int i, x, y, z;
   double temp_x, temp_y, temp_z;
   for (i = start; i <= end; i++) {
      x = 3 * i;
      y = x + 1;
      z = y + 1;
      temp_x = xyz[x] - cent_x;
      temp_y = xyz[y] - cent_y;
      temp_z = xyz[z] - cent_z;
      xyz[x] = cent_x + rotmat[0][0] * temp_x
          + rotmat[0][1] * temp_y + rotmat[0][2] * temp_z;
      xyz[y] = cent_y + rotmat[1][0] * temp_x
          + rotmat[1][1] * temp_y + rotmat[1][2] * temp_z;
      xyz[z] = cent_z + rotmat[2][0] * temp_x
          + rotmat[2][1] * temp_y + rotmat[2][2] * temp_z;
      if (fabs(xyz[x]) < 0.0000007)
         xyz[x] = ZERO;
      if (fabs(xyz[y]) < 0.0000007)
         xyz[y] = ZERO;
      if (fabs(xyz[z]) < 0.0000007)
         xyz[z] = ZERO;
   }
}



static void
calc_centroid(double *xyz, int start, int end, double *cent_x,
              double *cent_y, double *cent_z)
{
   int i, x, y, z, n;
   *cent_x = *cent_y = *cent_z = ZERO;
   for (i = start; i <= end; i++) {
      x = 3 * i;
      y = x + 1;
      z = y + 1;
      *cent_x += xyz[x];
      *cent_y += xyz[y];
      *cent_z += xyz[z];
   }
   n = (end - start) + 1;
   *cent_x /= n;
   *cent_y /= n;
   *cent_z /= n;
}


/*****
        ARPACK related functions:
*****/
static void absort(int n, double *crr, int *brr, int *error_flag)
{
   int i, j, b;
   double *arr, a;
   arr = (double *)
       my_malloc(malloc,
                 "\nERROR in absort/my_malloc(double *arr)", n,
                 sizeof(double), error_flag);
   if (*error_flag)
      return;
   for (i = 0; i < n; i++)
      arr[i] = crr[i];
   for (j = 1; j < n; j++) {
      a = arr[j];
      b = brr[j];
      i = j - 1;
      while (i >= 0 && fabs(arr[i]) > fabs(a)) {
         arr[i + 1] = arr[i];
         brr[i + 1] = brr[i];
         i--;
      }
      arr[i + 1] = a;
      brr[i + 1] = b;
   }
   my_free(arr);
}


static void
arpack(int *ndim, int *nof_requested_modes, int *nof_computed_modes,
       int *arpk_dim, double *eigvals, double *eigvecs, int *spectrum,
       int *want_eigvecs, double *xyz, double *grad, int *return_flag,
       int *label)
{
/* ARPACK: */
   static int arpk_debug_flag, arpk_iter, arpk_error_flag;
   static double arpk_tol;
   static double *v = NULL, *workl = NULL, *workd = NULL, *d =
       NULL, *resid = NULL, *ax = NULL;
   static int *select = NULL;
/* Local: */
   static int allocated, status_flag, error_flag;
   switch (*label) {
   case 0:
      arpk_debug_flag = DEBUG_ARPACK;
      arpk_iter = ARPK_ITER - 1;
      arpk_tol = ARPK_TOL;
      *arpk_dim = MIN(*ndim, *arpk_dim);
      allocated = NO;
      status_flag = 0;
      error_flag = FALSE;
      if (!allocated) {
         v = (double *)
             my_malloc(malloc,
                       "\nERROR in arpack/my_malloc(double *v)",
                       (*ndim) * (*arpk_dim), sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         workl = (double *)
             my_malloc(malloc,
                       "\nERROR in arpack/my_malloc(double *workl)",
                       (*arpk_dim) * ((*arpk_dim) + 8), sizeof(double),
                       &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         workd = (double *)
             my_malloc(malloc,
                       "\nERROR in arpack/my_malloc(double *workd)",
                       3 * (*ndim), sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         d = (double *)
             my_malloc(malloc,
                       "\nERROR in arpack/my_malloc(double *d)",
                       (*arpk_dim) * 2, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         resid = (double *)
             my_malloc(malloc,
                       "\nERROR in arpack/my_malloc(double *resid)",
                       (*ndim), sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         ax = (double *)
             my_malloc(malloc,
                       "\nERROR in arpack/my_malloc(double *ax)",
                       (*ndim), sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         select = (int *)
             my_malloc(malloc,
                       "\nERROR in arpack/my_malloc(int *select)",
                       (*arpk_dim), sizeof(int), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         allocated = YES;
      }
      goto L00;
   case 1:
      goto L01;
   default:
      fprintf(stderr, "\nERROR in arpack(): Illegal status.\n");
      fflush(stderr);
      if (allocated)
         allocated = NO;
      *label = ILLEGAL_STATUS;
      goto error_cleanup;
   }
 L00:
   for (status_flag = 0;;) {
    L01:
      dsarpack_(ndim, nof_requested_modes, nof_computed_modes, arpk_dim,
                &arpk_iter, &arpk_tol, eigvals, eigvecs, spectrum,
                want_eigvecs, &arpk_error_flag, &arpk_debug_flag, v,
                workl, workd, d, resid, ax, select, xyz, grad,
                return_flag, &status_flag);
      if (status_flag > 0) {
         *label = 1;
         return;                /* dsarpack continue */
      } else if (status_flag < 0) {
         if (allocated)
            allocated = NO;
         *label = status_flag;
         goto error_cleanup;    /* dsarpack error    */
      } else
         break;                 /* dsarpack done     */
   }
   if (allocated) {             /* deallocate local arrays */
      my_free(v);
      my_free(workl);
      my_free(workd);
      my_free(d);
      my_free(resid);
      my_free(ax);
      my_free(select);
      allocated = NO;
   }
   *label = 0;                  /* arpack() done */
   return;

 error_cleanup:
   my_free(v);
   my_free(workl);
   my_free(workd);
   my_free(d);
   my_free(resid);
   my_free(ax);
   my_free(select);
}


#ifdef HESSVEC_CENTRAL
int
hessvec_(int *ndim, double *vec_in, double *vec_out, double *xyz,
         double *grad, int *return_flag, int *label)
/*
    H(xyz)*v = ( grad(xyz+tiny_step*v) - grad(xyz-tiny_step*v) ) / (2*tiny_step)
*/
{
   static int i, n;
   static double *xyz_save = NULL, *grad_save = NULL, *grad_orig =
       NULL, sqrt_epsmach, tiny_step, dot, xyz_norm, vec_in_norm, max,
       vec_in_max;
   static int allocated, error_flag;
   switch (*label) {
   case 0:
      allocated = NO;
      error_flag = FALSE;
      n = (*ndim);
      if (!allocated) {
         xyz_save = (double *)
             my_malloc(malloc,
                       "\nERROR in hessvec/my_malloc(double *xyz_save)",
                       n, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         grad_save = (double *)
             my_malloc(malloc,
                       "\nERROR in hessvec/my_malloc(double *grad_save)",
                       n, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         grad_orig = (double *)
             my_malloc(malloc,
                       "\nERROR in hessvec/my_malloc(double *grad_orig)",
                       n, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         allocated = YES;
      }
      sqrt_epsmach = sqrt(DBL_EPSILON);
      memcpy(grad_orig, grad, n * sizeof(double));
      goto L00;
   case 1:
      goto L01;
   case 2:
      goto L02;
   default:
      fprintf(stderr, "\nERROR in hessvec(): Illegal status.\n");
      fflush(stderr);
      if (allocated)
         allocated = NO;
      *label = ILLEGAL_STATUS;
      goto error_cleanup;
   }
 L00:
   memcpy(xyz_save, xyz, n * sizeof(double));
   for (i = 0, dot = ZERO; i < n; i++)
      dot += SQR(xyz_save[i]);
   xyz_norm = sqrt(dot);
   for (i = 0, dot = ZERO; i < n; i++)
      dot += SQR(vec_in[i]);
   vec_in_norm = sqrt(dot);
   for (i = 0, vec_in_max = ZERO; i < n; i++)
      if ((max = fabs(vec_in[i])) >= vec_in_max)
         vec_in_max = max;
   /* Derreumaux, Zhang, Schlick, Brooks J. Comput. Chem. 15, 532-552 (1994), p. 541: */
   tiny_step =
       MIN((2. * sqrt_epsmach * (ONE + xyz_norm) / vec_in_norm),
           (sqrt_epsmach / vec_in_max));
   for (i = 0; i < n; i++)
      xyz[i] -= tiny_step * vec_in[i];
   *return_flag = CALCGRAD_OLDNBL;
   *label = 1;
   return;
 L01:
   memcpy(grad_save, grad, n * sizeof(double));
   memcpy(xyz, xyz_save, n * sizeof(double));
   for (i = 0; i < n; i++)
      xyz[i] += tiny_step * vec_in[i];
   *return_flag = CALCGRAD_OLDNBL;
   *label = 2;
   return;
 L02:
   for (i = 0; i < n; i++)
      vec_out[i] = (grad[i] - grad_save[i]) / (2. * tiny_step); /* load vec_out[] */
   memcpy(xyz, xyz_save, n * sizeof(double));   /* restore  xyz[] */
   memcpy(grad, grad_orig, n * sizeof(double)); /* restore grad[] */
   if (allocated) {
      my_free(xyz_save);
      my_free(grad_save);
      my_free(grad_orig);
      allocated = NO;
   }
   *label = 0;                  /* hessvec_() done */
   return;

 error_cleanup:
   my_free(xyz_save);
   my_free(grad_save);
   my_free(grad_orig);
}


#else
int
hessvec_(int *ndim, double *vec_in, double *vec_out, double *xyz,
         double *grad, int *return_flag, int *label)
/*
    H(xyz)*v = ( grad(xyz+tiny_step*v) - grad(xyz) ) / tiny_step
    
    !!! Make sure that on entry grad[] is up-to-date wrt xyz[], because it is not calculated here. !!!
*/
{
   static int i, n;
   static double *xyz_save = NULL, *grad_save =
       NULL, sqrt_epsmach, tiny_step, dot, xyz_norm, vec_in_norm;
   static int allocated, error_flag;
   switch (*label) {
   case 0:
      allocated = NO;
      error_flag = FALSE;
      n = (*ndim);
      if (!allocated) {
         xyz_save = (double *)
             my_malloc(malloc,
                       "\nERROR in hessvec/my_malloc(double *xyz_save)",
                       n, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         grad_save = (double *)
             my_malloc(malloc,
                       "\nERROR in hessvec/my_malloc(double *grad_save)",
                       n, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         allocated = YES;
      }
      sqrt_epsmach = sqrt(DBL_EPSILON);
      goto L00;
   case 1:
      goto L01;
   default:
      fprintf(stderr, "\nERROR in hessvec(): Illegal status.\n");
      fflush(stderr);
      if (allocated)
         allocated = NO;
      *label = ILLEGAL_STATUS;
      goto error_cleanup;
   }
 L00:
   memcpy(grad_save, grad, n * sizeof(double));
   memcpy(xyz_save, xyz, n * sizeof(double));
   for (i = 0, dot = ZERO; i < n; i++)
      dot += SQR(xyz_save[i]);
   xyz_norm = sqrt(dot);
   for (i = 0, dot = ZERO; i < n; i++)
      dot += SQR(vec_in[i]);
   vec_in_norm = sqrt(dot);
   /* Derreumaux, Zhang, Schlick, Brooks J. Comput. Chem. 15, 532-552 (1994), p. 541: */
   tiny_step = 2. * sqrt_epsmach * (ONE + xyz_norm) / vec_in_norm;
   for (i = 0; i < n; i++)
      xyz[i] += tiny_step * vec_in[i];
   *return_flag = CALCGRAD_OLDNBL;
   *label = 1;
   return(0);
 L01:
   for (i = 0; i < n; i++)
      vec_out[i] = (grad[i] - grad_save[i]) / tiny_step;        /* load vec_out[] */
   memcpy(xyz, xyz_save, n * sizeof(double));   /* restore  xyz[] */
   memcpy(grad, grad_save, n * sizeof(double)); /* restore grad[] */
   if (allocated) {
      my_free(xyz_save);
      my_free(grad_save);
      allocated = NO;
   }
   *label = 0;                  /* hessvec_() done */
   return(0);

 error_cleanup:
   my_free(xyz_save);
   my_free(grad_save);
   return(0);
}
#endif

#if 0
static void
hessdump(int *ndim, double *xyz, double *grad, int *return_flag,
         int *label)
{
   static int i, j, n;
   static double *unit_vector = NULL, *hessian_row = NULL;
   static FILE *file;
   static int allocated, status_flag, error_flag;
   switch (*label) {
   case 0:
      allocated = NO;
      status_flag = 0;
      error_flag = FALSE;
      n = (*ndim);
      if (!allocated) {
         unit_vector = (double *)
             my_malloc(malloc,
                       "\nERROR in hessdump/my_malloc(double *unit_vector)",
                       n, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         hessian_row = (double *)
             my_malloc(malloc,
                       "\nERROR in hessdump/my_malloc(double *hessian_row)",
                       n, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         allocated = YES;
      }
      if ((file = fopen("hessian.mat", "w")) == NULL) {
         fprintf(stderr, "\nERROR in hessdump(): Cannot open %s.\n",
                 "hessian.mat");
         fflush(stderr);
         if (allocated)
            allocated = NO;
         *label = FILE_ERROR;
         goto error_cleanup;
      }
      goto L00;
   case 1:
      goto L01;
   default:
      fprintf(stderr, "\nERROR in hessdump(): Illegal status.\n");
      fflush(stderr);
      if (allocated)
         allocated = NO;
      *label = ILLEGAL_STATUS;
      goto error_cleanup;
   }
 L00:
   for (i = 0; i < n; i++, fprintf(file, "\n")) {
      for (j = 0; j < n; j++)
         if (j == i)
            unit_vector[j] = ONE;
         else
            unit_vector[j] = ZERO;
      for (status_flag = 0;;) {
       L01:
         hessvec_(ndim, unit_vector, hessian_row,
                  xyz, grad, return_flag, &status_flag);
         if (status_flag > 0) {
            *label = 1;
            return;             /* hessvec continue */
         } else if (status_flag < 0) {
            if (allocated)
               allocated = NO;
            *label = status_flag;
            goto error_cleanup; /* hessvec error    */
         } else
            break;              /* hessvec done     */
      }                         /* end hessvec() */
      for (j = 0; j < n; j++)
         fprintf(file, "  %13.5e", hessian_row[j]);
   }
   fclose(file);
   if (allocated) {
      my_free(unit_vector);
      my_free(hessian_row);
      allocated = NO;
   }
   *label = 0;                  /* hessdump() done */
   return;

 error_cleanup:
   my_free(unit_vector);
   my_free(hessian_row);
}
#endif

static void
separate_close_pairs(int do_all_pairs, int do_ligs_only, int ndim,
                     double *xyz, double *grad, int nligs, int *lig_start,
                     int *lig_end)
/* Atoms closer than 1 A will be moved at least 1.2 A away
   along their connecting axis. (after Peter Shenkin) */
{
   int i, ii, j, k, natm = ndim / 3, problem, max_separ_cycle =
       5, i_is_frozen, j_is_frozen;
   double dx, dy, dz, dist, dist2, frac_separ, min_separ_dist2 = 1.0;
   if (do_all_pairs) {
      for (k = 0, problem = FALSE; k < max_separ_cycle;
           k++, problem = FALSE) {
         for (i = 0; i < natm; i++) { /* grad[] need not be up-to-date, it is
                                         only used for catching frozen atoms */
            if (grad[i * 3] == ZERO && grad[i * 3 + 1] == ZERO
                && grad[i * 3 + 2] == ZERO)
               i_is_frozen = YES;
            else
               i_is_frozen = NO;
            for (j = i + 1; j < natm; j++) {
               if (grad[j * 3] == ZERO && grad[j * 3 + 1] == ZERO
                   && grad[j * 3 + 2] == ZERO)
                  j_is_frozen = YES;
               else
                  j_is_frozen = NO;
               if (i_is_frozen && j_is_frozen)
                  continue;
               dx = xyz[i * 3    ] - xyz[j * 3];
               dy = xyz[i * 3 + 1] - xyz[j * 3 + 1];
               dz = xyz[i * 3 + 2] - xyz[j * 3 + 2];
               dist2 = SQR(dx) + SQR(dy) + SQR(dz);
               if (dist2 < min_separ_dist2) {
                  problem = TRUE;
                  if (dist2 == ZERO)
                     if (i_is_frozen || j_is_frozen)
                        dx = dy = dz = sqrt(2. / 5.);
                     else
                        dx = dy = dz = sqrt(1. / 5.);
                  else {
                     dist = sqrt(dist2);
                     frac_separ = (1.2 - dist) / (2.0 * dist);
                     if (i_is_frozen || j_is_frozen)
                        frac_separ *= 2.0;
                     dx *= frac_separ;
                     dy *= frac_separ;
                     dz *= frac_separ;
                  }
                  if (!i_is_frozen) {
                     xyz[i * 3    ] += dx;
                     xyz[i * 3 + 1] += dy;
                     xyz[i * 3 + 2] += dz;
                  }
                  if (!j_is_frozen) {
                     xyz[j * 3    ] -= dx;
                     xyz[j * 3 + 1] -= dy;
                     xyz[j * 3 + 2] -= dz;
                  }
               }
            }
         }
         if (!problem)
            break;
      }
   }
   if (do_ligs_only) {
      for (k = 0, problem = FALSE; k < max_separ_cycle;
           k++, problem = FALSE) {
         for (ii = 0; ii < nligs; ii++)
            for (i = lig_start[ii] - 1; i < lig_end[ii]; i++) {
               if (grad[i * 3] == ZERO && grad[i * 3 + 1] == ZERO
                   && grad[i * 3 + 2] == ZERO)
                  i_is_frozen = YES;
               else
                  i_is_frozen = NO;
               for (j = 0; j < lig_start[ii] - 1; j++) {
                  if (grad[j * 3] == ZERO && grad[j * 3 + 1] == ZERO
                      && grad[j * 3 + 2] == ZERO)
                     j_is_frozen = YES;
                  else
                     j_is_frozen = NO;
                  if (i_is_frozen && j_is_frozen)
                     continue;
                  dx = xyz[i * 3    ] - xyz[j * 3];
                  dy = xyz[i * 3 + 1] - xyz[j * 3 + 1];
                  dz = xyz[i * 3 + 2] - xyz[j * 3 + 2];
                  dist2 = SQR(dx) + SQR(dy) + SQR(dz);
                  if (dist2 < min_separ_dist2) {
                     problem = TRUE;
                     if (dist2 == ZERO)
                        if (i_is_frozen || j_is_frozen)
                           dx = dy = dz = sqrt(2. / 5.);
                        else
                           dx = dy = dz = sqrt(1. / 5.);
                     else {
                        dist = sqrt(dist2);
                        frac_separ = (1.2 - dist) / (2.0 * dist);
                        if (i_is_frozen || j_is_frozen)
                           frac_separ *= 2.0;
                        dx *= frac_separ;
                        dy *= frac_separ;
                        dz *= frac_separ;
                     }
                     if (!i_is_frozen) {
                        xyz[i * 3    ] += dx;
                        xyz[i * 3 + 1] += dy;
                        xyz[i * 3 + 2] += dz;
                     }
                     if (!j_is_frozen) {
                        xyz[j * 3    ] -= dx;
                        xyz[j * 3 + 1] -= dy;
                        xyz[j * 3 + 2] -= dz;
                     }
                  }
               }
               for (j = lig_end[ii]; j < natm; j++) {
                  if (grad[j * 3] == ZERO && grad[j * 3 + 1] == ZERO
                      && grad[j * 3 + 2] == ZERO)
                     j_is_frozen = YES;
                  else
                     j_is_frozen = NO;
                  if (i_is_frozen && j_is_frozen)
                     continue;
                  dx = xyz[i * 3    ] - xyz[j * 3];
                  dy = xyz[i * 3 + 1] - xyz[j * 3 + 1];
                  dz = xyz[i * 3 + 2] - xyz[j * 3 + 2];
                  dist2 = SQR(dx) + SQR(dy) + SQR(dz);
                  if (dist2 < min_separ_dist2) {
                     problem = TRUE;
                     if (dist2 == ZERO)
                        if (i_is_frozen || j_is_frozen)
                           dx = dy = dz = sqrt(2. / 5.);
                        else
                           dx = dy = dz = sqrt(1. / 5.);
                     else {
                        dist = sqrt(dist2);
                        frac_separ = (1.2 - dist) / (2.0 * dist);
                        if (i_is_frozen || j_is_frozen)
                           frac_separ *= 2.0;
                        dx *= frac_separ;
                        dy *= frac_separ;
                        dz *= frac_separ;
                     }
                     if (!i_is_frozen) {
                        xyz[i * 3    ] += dx;
                        xyz[i * 3 + 1] += dy;
                        xyz[i * 3 + 2] += dz;
                     }
                     if (!j_is_frozen) {
                        xyz[j * 3    ] -= dx;
                        xyz[j * 3 + 1] -= dy;
                        xyz[j * 3 + 2] -= dz;
                     }
                  }
               }
            }
         if (!problem)
            break;
      }
   }
}


/*****
        Reverse communication LMOD function:
*****/
double
lmodC(int *nlmodit, int *nmod, int *kmod, int *rotran, int *natm_ext,
      double *xyz_ext, double *enrg, double *grad_ext, int *nconf,
      double *enrg_win, double *separ_rms, double *conflib, double *trajectory,
      int *arpk_recalc, int *arpk_dim, int *lmod_restart, int *topten,
      int *mc_option, double *rtemp, double *lmod_step_min,
      double *lmod_step_max, int *nof_lmod_step, int *nligs,
      int *lig_start, int *lig_end, int *try_now, int *ntry, double *trmin,
      double *trmax, int *lig_rot_cent, double *angmin, double *angmax,
      int *seed, int *print_level, double *lmod_time, double *extn_time,
      int *return_flag, int *label)
{
/* ARPACK: */
   static int nof_requested_modes, nof_computed_modes;
   static int spectrum, want_eigvecs;
   static double *eigvals = NULL, *eigvecs = NULL;
/* Local: */
   static int seed3;            /* see rand2.c */
   static int allocated, status_flag, error_flag;
   static int lmod_iter, restart, zigzag_iter, max_zigzag_iter,
       i, j, k, kk, l, n, cnt;
   static int barrier_crossing_test_on, do_all, do_ligs;
   static int *index = NULL;
   static int ikk, lig_start_kk, lig_end_kk, lig_rot_cent_kk;
   static double ref_energy, energy, energy_old, min_energy,
       glob_min_energy, rad, sum, max_atmov, scale, lmod_step, rms,
       rms_old, grad_rms;
   static double xtrans, ytrans, ztrans, trnorm, trscale, xcent, ycent,
       zcent, xrot, yrot, zrot, rotnorm, rotang, rotmat[3][3];
   static double *xyz_cent = NULL, *lmod_vec = NULL, *pboltz = NULL,
       *xyz_hold = NULL, *lmod_move = NULL, *trot_move = NULL;
   static struct archive *walk;
   static int sizeof_archive;
   static clock_t lmod_time_stamp, extn_time_stamp;
   static int natm_local, ndim_local, ndim_ext, *atm_indx = NULL;
   static double *xyz_local = NULL, *grad_local = NULL;

#ifdef SQM
   nabout = stdout;
#endif

   switch (*label) {
   case 0:
      *lmod_time = ZERO;
      *extn_time = ZERO;
      lmod_time_stamp = clock();
/*
                    Parameter check list:
*/
      if (*nlmodit < 0) {
         fprintf(stderr,
         "\nERROR in lmod(): Requested number of LMOD iterations negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*nmod < 1) {
         fprintf(stderr,
         "\nERROR in lmod(): Wrong number of vibrational modes.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*kmod < 1) {
         fprintf(stderr,
         "\nERROR in lmod(): Wrong number of variable modes.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*kmod > *nmod) {
         fprintf(stderr,
         "\nERROR in lmod(): Number of variable modes > total");
         fprintf(stderr, " number of modes.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*rotran != 6 && *rotran != 3 && *rotran != 1 && *rotran != 0) {
         fprintf(stderr,
         "\nERROR in lmod(): Wrong number of trans/rot modes.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*natm_ext < 0) {
         fprintf(stderr, "\nERROR in lmod(): Number of atoms negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*natm_ext < 2) {
         fprintf(stderr, "\nERROR in lmod(): Too few atoms.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*nconf <= 0) {
         fprintf(stderr,
         "\nERROR in lmod(): Requested number of confs zero or negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*enrg_win < ZERO) {
         fprintf(stderr, "\nERROR in lmod(): Energy window negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*separ_rms <= ZERO) {
         fprintf(stderr, 
             "\nERROR in lmod(): Conf separation RMS zero or negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*arpk_recalc <= 0) {
         fprintf(stderr,
         "\nERROR in lmod(): Eigvec update frequency zero or negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*arpk_dim < 0) {
         fprintf(stderr,
         "\nERROR in lmod(): Requested number of Arnoldi basis");
         fprintf(stderr, " vectors negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*arpk_dim > 0 && *arpk_dim < (*nmod + 1)) {
         fprintf(stderr,
         "\nERROR in lmod(): Number of Arnoldi basis vectors");
         fprintf(stderr, " should be at least one\n");
         fprintf(stderr,
         "                 more than the requested number");
         fprintf(stderr, " of vibrational modes.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*lmod_restart <= 0) {
         fprintf(stderr,
         "\nERROR in lmod(): LMOD restart frequency zero or negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*topten <= 0) {
         fprintf(stderr,
         "\nERROR in lmod(): Pool of structures for LMOD restart empty.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*mc_option < 1 || *mc_option > 3) {
         fprintf(stderr,
         "\nERROR in lmod(): Unknown Monte Carlo option.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*rtemp <= ZERO) {
         fprintf(stderr, "\nERROR in lmod(): RT should be positive.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*lmod_step_min <= ZERO) {
         fprintf(stderr,
         "\nERROR in lmod(): Min. LMOD step size should be positive.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*lmod_step_max <= ZERO) {
         fprintf(stderr,
         "\nERROR in lmod(): Max. LMOD step size should be positive.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*lmod_step_min > *lmod_step_max) {
         fprintf(stderr,
          "\nERROR in lmod(): LMOD step min/max reversed.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*nof_lmod_step < 0) {
         fprintf(stderr,
         "\nERROR in lmod(): Requested number of LMOD steps negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*nligs < 0) {
         fprintf(stderr,
         "\nERROR in lmod(): User-defined number of ligands negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      for (i = 0; i < (*nligs); i++) {
         if (lig_start[i] < 1 || lig_start[i] > (*natm_ext)
             || lig_end[i] < 1 || lig_end[i] > (*natm_ext)
             || lig_start[i] > lig_end[i]) {
            fprintf(stderr,
            "\nERROR in lmod(): Ligand %d's start/end out of range.\n",
            i + 1);
            fflush(stderr);
            *label = PARAMS_ERROR;
            return(0.0);
         }
         if (trmin[i] > trmax[i]) {
            fprintf(stderr,
            "\nERROR in lmod(): Ligand %d's trans min/max reversed.\n",
            i + 1);
            fflush(stderr);
            *label = PARAMS_ERROR;
            return(0.0);
         }
         if (angmin[i] > angmax[i]) {
            fprintf(stderr,
            "\nERROR in lmod(): Ligand %d's rot min/max reversed.\n",
            i + 1);
            fflush(stderr);
            *label = PARAMS_ERROR;
            return(0.0);
         }
         if (lig_rot_cent[i] != 0
             && (lig_rot_cent[i] < 0 || lig_rot_cent[i] < lig_start[i]
                 || lig_rot_cent[i] > lig_end[i])) {
            fprintf(stderr,
            "\nERROR in lmod(): Ligand %d's rot center out of range.\n",
            i + 1);
            fflush(stderr);
            *label = PARAMS_ERROR;
            return(0.0);
         }
      }
      if (*try_now <= 0) {
         fprintf(stderr,
         "\nERROR in lmod(): LMOD iter vs. ligand trans/rot ");
         fprintf(stderr, "frequency must be positive.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*ntry < 0) {
         fprintf(stderr,
         "\nERROR in lmod(): Requested number of ligand trans/rot");
         fprintf(stderr, " moves negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*seed < 0) {
         fprintf(stderr,
         "\nERROR in lmod(): User-defined random seed negative.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
      if (*print_level == 0) {
         DEBUG_ARPACK = NO;
         PRINT_LMOD = NO;
         DEBUG_LMOD = NO;
      } else if (*print_level == 1) {
         DEBUG_ARPACK = NO;
         PRINT_LMOD = YES;
         DEBUG_LMOD = NO;
      } else if (*print_level == 2) {
         DEBUG_ARPACK = NO;
         PRINT_LMOD = YES;
         DEBUG_LMOD = YES;
      } else if (*print_level == 3) {
         DEBUG_ARPACK = YES;
         PRINT_LMOD = YES;
         DEBUG_LMOD = YES;
      } else if (*print_level == 4) {
         DEBUG_ARPACK = YES;
         PRINT_LMOD = NO;
         DEBUG_LMOD = NO;
      } else {
         fprintf(stderr, "\nERROR in lmod(): Print level out of range.\n");
         fflush(stderr);
         *label = PARAMS_ERROR;
         return(0.0);
      }
/*
        Check for frozen atoms and create a local
        environment with only moving atoms:
*/
      natm_local = (*natm_ext);
      for (i=0; i<(*natm_ext); i++) {
         if (grad_ext[i*3]==0 && grad_ext[i*3+1]==0 && grad_ext[i*3+2]==0)
            natm_local--;
      }
      ndim_local = 3 * natm_local;
      ndim_ext   = 3 *(*natm_ext);
/*
	Parameters OK, start calculation:
*/
      allocated = NO;
      status_flag = 0;
      error_flag = FALSE;
      restart = NO;
      if (*arpk_dim == 0)
         *arpk_dim = ndim_local;
      else {
         *arpk_dim += (*rotran);
         *arpk_dim = MIN(*arpk_dim, ndim_local);
      }
      spectrum = ARPK_SM;
      want_eigvecs = YES;
      nof_requested_modes = MIN((*nmod), (ndim_local - 1));
      nof_requested_modes = MIN(nof_requested_modes + (*rotran),
                                ndim_local - 1);  /* plus ro-trans modes */
      if ((*nof_lmod_step) == 0) {
         max_zigzag_iter = MAX_ZIGZAG_ITER;
         barrier_crossing_test_on = YES;
      } else {
         max_zigzag_iter = (*nof_lmod_step);
         barrier_crossing_test_on = NO;
      }
      if (!allocated) {
         eigvals = (double *)
             my_malloc(malloc,
                       "\nERROR in lmod/my_malloc(double *eigvals)",
                       nof_requested_modes, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         eigvecs = (double *)
             my_malloc(malloc,
                       "\nERROR in lmod/my_malloc(double *eigvecs)",
                       ndim_local * nof_requested_modes, sizeof(double),
                       &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         index = (int *)
             my_malloc(malloc,
                       "\nERROR in lmod/my_malloc(int *index)", ndim_local,
                       sizeof(int), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         xyz_cent = (double *)
             my_malloc(malloc,
                       "\nERROR in lmod/my_malloc(double *xyz_cent)",
                       ndim_local, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         lmod_vec = (double *)
             my_malloc(malloc,
                       "\nERROR in lmod/my_malloc(double *lmod_vec)",
                       ndim_local, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         pboltz = (double *)
             my_malloc(malloc,
                       "\nERROR in lmod/my_malloc(double *pboltz)",
                       2 * nof_requested_modes, sizeof(double),
                       &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         xyz_hold = (double *)
             my_malloc(malloc,
                       "\nERROR in lmod/my_malloc(double *xyz_hold)",
                       ndim_local, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         lmod_move = (double *)
             my_malloc(malloc,
                       "\nERROR in lmod/my_malloc(double *lmod_move)",
                       2 * nof_requested_modes * ndim_local, sizeof(double),
                       &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         trot_move = (double *)
             my_malloc(malloc,
                       "\nERROR in lmod/my_malloc(double *trot_move)",
                       ndim_local, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         xyz_local = (double *)
             my_malloc(malloc,
                     "\nERROR in lmod/my_malloc(double *xyz_local)",
                     ndim_local, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         grad_local = (double *)
             my_malloc(malloc,
                     "\nERROR in lmod/my_malloc(double *grad_local)",
                     ndim_local, sizeof(double), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         atm_indx = (int *)
             my_malloc(malloc,
                     "\nERROR in lmod/my_malloc(int *atm_indx)",
                     natm_local, sizeof(int), &error_flag);
         if (error_flag) {
            *label = error_flag;
            goto error_cleanup;
         }
         allocated = YES;
      }
      for (i=j=0; i<(*natm_ext); i++) {  /* generate local -> ext mapping */
         if (grad_ext[i*3]!=0 || grad_ext[i*3+1]!=0 || grad_ext[i*3+2]!=0)
            atm_indx[j++] = i;
      }
      for (i=0; i<natm_local; i++) {  /* load xyz_local[] */
         j = atm_indx[i];
         xyz_local[3*i  ] = xyz_ext[3*j  ];
         xyz_local[3*i+1] = xyz_ext[3*j+1];
         xyz_local[3*i+2] = xyz_ext[3*j+2];
      }
      /* Make sure that all MPI processes use the same seed */
      if (get_mytaskid() == 0) {
         if (*seed > 0)
            seed3 = -(*seed);   /* see rand2.c */
         else
            seed3 = rseed();    /* see rand2.c */
      }
#if defined(MPI) || defined(SCALAPACK)
      MPI_Bcast(&seed3, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
      setseed(&seed3);          /* see rand2.c */
      goto L00;
   case 1:
      lmod_time_stamp = clock();
      *extn_time += (double) (lmod_time_stamp - extn_time_stamp);
      goto L01;
   case 2:
      lmod_time_stamp = clock();
      *extn_time += (double) (lmod_time_stamp - extn_time_stamp);
      goto L02;
   case 3:
      lmod_time_stamp = clock();
      *extn_time += (double) (lmod_time_stamp - extn_time_stamp);
      goto L03;
   case 4:
      lmod_time_stamp = clock();
      *extn_time += (double) (lmod_time_stamp - extn_time_stamp);
      goto L04;
   case 5:
      lmod_time_stamp = clock();
      *extn_time += (double) (lmod_time_stamp - extn_time_stamp);
      goto L05;
   case 6:
      lmod_time_stamp = clock();
      *extn_time += (double) (lmod_time_stamp - extn_time_stamp);
      goto L06;
   case 7:
      lmod_time_stamp = clock();
      *extn_time += (double) (lmod_time_stamp - extn_time_stamp);
      goto L07;
   default:
      fprintf(stderr, "\nERROR in lmod(): Illegal status.\n");
      fflush(stderr);
      if (allocated)
         allocated = NO;
      *label = ILLEGAL_STATUS;
      goto error_cleanup;
   }
 L00:
/*    for( status_flag=0 ;; ){  */
 L01:
#if 0
        /* Update grad_local[], hessdump needs new grad on input:*/
        for (i=0; i<natm_local; i++) {
           j = atm_indx[i];
           grad_local[3*i  ] = grad_ext[3*j  ];
           grad_local[3*i+1] = grad_ext[3*j+1];
           grad_local[3*i+2] = grad_ext[3*j+2];
        }
        hessdump( &ndim_local, xyz_local, grad_local,
                  return_flag, &status_flag );
        if( status_flag > 0 ){
            /* compute grad at new xyz: */
            for (i=0; i<natm_local; i++) {  /* update xyz_ext[] */
               j = atm_indx[i];
               xyz_ext[ 3*j  ] = xyz_local[ 3*i  ];
               xyz_ext[ 3*j+1] = xyz_local[ 3*i+1];
               xyz_ext[ 3*j+2] = xyz_local[ 3*i+2];
            }
            extn_time_stamp = clock();
            *lmod_time += (double)(extn_time_stamp - lmod_time_stamp);
            *label = 1;
            return;                     /* hessdump continue */
        }
        else if( status_flag < 0 ){
            if( allocated )
                allocated = NO;
            *label = status_flag;
            goto error_cleanup;         /* hessdump error    */
        }
        else
            break;                      /* hessdump done     */
    } /* end hessdump() */
#endif

/*
    Initialize conflib if "conflib.dat" found:
*/
   if( fopen( "conflib.dat", "rb" ) != NULL )    {

      if (get_mytaskid() == 0) {
         fprintf( nabout, "\n*****  W A R N I N G  *****\n\n" );
         fprintf( nabout, "It appears that results from a previous LMOD search were found\n" );
         fprintf( nabout, "in conflib.dat in the launch directory.  LMOD will assume that\n" );
         fprintf( nabout, "those molecules in conflib.dat are low-energy conformations of\n" );
         fprintf( nabout, "the current input molecule of the present LMOD run and will be\n" );
         fprintf( nabout, "used to seed the new search with multiple conformations.\n\n" );
         fprintf( nabout, "If conflib.dat is incompatible with the current search, *KILL*\n" );
         fprintf( nabout, "this LMOD job, remove conflib.dat, and restart LMOD.\n\n" );
      }

      read_and_sort_archive( "conflib.dat", ndim_ext, &error_flag );
      if (error_flag) {
         if (allocated)
            allocated = NO;
         *label = error_flag;
         goto error_cleanup;
      }
   }

/*
    Minimize initial structure:
*/
   for (i=0; i<natm_local; i++) {  /* update xyz_ext[] */
      j = atm_indx[i];
      xyz_ext[ 3*j  ] = xyz_local[ 3*i  ];
      xyz_ext[ 3*j+1] = xyz_local[ 3*i+1];
      xyz_ext[ 3*j+2] = xyz_local[ 3*i+2];
   }
   *return_flag = MINIMIZE;
   extn_time_stamp = clock();
   *lmod_time += (double) (extn_time_stamp - lmod_time_stamp);
   *label = 2;
   return(0.0);
 L02:
   for (i=0; i<natm_local; i++) {  /* update xyz_local[] w/ minimized coords */
      j = atm_indx[i];
      xyz_local[3*i  ] = xyz_ext[3*j  ];
      xyz_local[3*i+1] = xyz_ext[3*j+1];
      xyz_local[3*i+2] = xyz_ext[3*j+2];
   }
   energy = *enrg;
   rad = calc_rad(natm_local, xyz_local);
   /* archive structure in conflib: */
   archive_structure(ndim_ext, xyz_ext, energy, rad, &error_flag);
   if (error_flag) {
      if (allocated)
         allocated = NO;
      *label = error_flag;
      goto error_cleanup;
   }
   /* start recording LMOD trajectory: */
   memcpy(trajectory, xyz_ext, ndim_ext * sizeof(double));
   if (get_mytaskid() == 0) {
      if (PRINT_LMOD) {
         fprintf( nabout, "__________________________________________________");
         if (DEBUG_LMOD)
            fprintf( nabout, "_________________________");
         fprintf(nabout, "\n                  Low-Mode Simulation\n");
         if (!DEBUG_LMOD)
            fprintf( nabout, "-----------------------------------------------\n");
      }
   }
/*
    M A I N  I T E R A T I O N  L O O P:
*/
   for (lmod_iter = 0; lmod_iter < (*nlmodit); lmod_iter++) {
/*
        Calc energy and gradient:
*/
      for (i=0; i<natm_local; i++) {  /* update xyz_ext[] */
         j = atm_indx[i];
         xyz_ext[ 3*j  ] = xyz_local[ 3*i  ];
         xyz_ext[ 3*j+1] = xyz_local[ 3*i+1];
         xyz_ext[ 3*j+2] = xyz_local[ 3*i+2];
      }
      *return_flag = CALCBOTH_NEWNBL;   /* grad[] is required for forward
                                           hessvec formula as reference! */
      extn_time_stamp = clock();
      *lmod_time += (double) (extn_time_stamp - lmod_time_stamp);
      *label = 3;
      return(0.0);
    L03:
      for (i=0; i<natm_local; i++) {  /* update grad_local[] */
         j = atm_indx[i];
         grad_local[3*i  ] = grad_ext[3*j  ];
         grad_local[3*i+1] = grad_ext[3*j+1];
         grad_local[3*i+2] = grad_ext[3*j+2];
      }
      ref_energy = *enrg;
      for (i = 0, grad_rms = ZERO; i < ndim_local; i++)
         grad_rms += SQR(grad_local[i]);
      grad_rms = sqrt(grad_rms / ndim_local);
/*
        Calculate low-mode eigenvectors:
*/
      if ((lmod_iter % (*arpk_recalc)) && !restart)
         goto SKIP_ARPK;
      restart = NO;
      for (status_flag = 0;;) {
       L04:
         for (i=0; i<natm_local; i++) {  /* update grad_local[] w/ new grad
                                            ARPACK needs grad for computing
                                            Hess * some_vector in hessvec() */
            j = atm_indx[i];
            grad_local[3*i  ] = grad_ext[3*j  ];
            grad_local[3*i+1] = grad_ext[3*j+1];
            grad_local[3*i+2] = grad_ext[3*j+2];
         }
         arpack(&ndim_local, &nof_requested_modes, &nof_computed_modes,
                arpk_dim, eigvals, eigvecs, &spectrum, &want_eigvecs,
                xyz_local, grad_local, return_flag, &status_flag);
         if (status_flag > 0) {
            for (i=0; i<natm_local; i++) {  /* update xyz_ext for
                                               computing new grad */
               j = atm_indx[i];
               xyz_ext[ 3*j  ] = xyz_local[ 3*i  ];
               xyz_ext[ 3*j+1] = xyz_local[ 3*i+1];
               xyz_ext[ 3*j+2] = xyz_local[ 3*i+2];
            }
            extn_time_stamp = clock();
            *lmod_time += (double) (extn_time_stamp - lmod_time_stamp);
            *label = 4;
            return(0.0);             /* arpack continue */
         } else if (status_flag < 0) {
            if (allocated)
               allocated = NO;
            *label = status_flag;
            goto error_cleanup; /* arpack error    */
         } else {
#if 0
             if (get_mytaskid() == 0 ) {
               fprintf( nabout, "\nEigenvalues\n\n");
               for(i=0; i<nof_computed_modes; i++)
                   fprintf(nabout, "%12.5g ",eigvals[i]);
               fprintf(nabout, "\n\nEigenvectors\n\n");
               for(i=0; i<ndim; i++){
                   for(j=0; j<nof_computed_modes; j++)
                       fprintf(nabout, "%12.5g ",eigvecs[j*ndim+i]);
                   fprintf(nabout, "\n");
               }
            }
#endif
            break;              /* arpack done     */
         }
      }                         /* end arpack() */
      for (i = 0; i < nof_computed_modes; i++)
         index[i] = i;
      /* sort eigvals in ascending absolute value: */
      if (error_flag) {
      absort(nof_computed_modes, eigvals, index, &error_flag);
         if (allocated)
            allocated = NO;
         *label = error_flag;
         goto error_cleanup;
      }
    SKIP_ARPK:
      /* save coords of central structure: */
      memcpy(xyz_cent, xyz_local, ndim_local * sizeof(double));
      if (get_mytaskid() == 0) {
         if (DEBUG_LMOD) {
            fprintf (nabout, "--------------------------------------");
            fprintf (nabout, "-------------------------------------\n");
         }
         if (PRINT_LMOD) {
            fprintf(nabout, "%6d   %cE =%13.3f (%6.3f)  Rg = %8.3f\n",
                   lmod_iter + 1, ((*nligs)
                                   && (!(lmod_iter % (*arpk_recalc)))
                                   && lmod_iter) ? '*' : ' ', ref_energy,
                   grad_rms, calc_rad(natm_local, xyz_local));
            fflush(nabout);
         }
      }
/*
        Cycle through vibrational modes:
*/
      for (k = (*rotran), kk = 0, cnt = 0, min_energy = BIG;
           k < nof_computed_modes; k++) {
         /* randomly select (*kmod) vib. modes to be varied: */
         if (cnt == *kmod)
            break;
         if (*kmod < (nof_computed_modes - *rotran)
             && (*kmod - cnt) < (nof_computed_modes - k)
             && rand2() >=
             ((double) (*kmod) / (double) (nof_computed_modes - *rotran)))
            continue;
         cnt++;
         for (i = 0; i < ndim_local; i++)
            /* pull out eigvec: */
            lmod_vec[i] = eigvecs[index[k] * ndim_local + i];
         for (l = 1; l <= 2; l++) {
/*
   Climb energy barrier utilizing the LMOD ZIG-ZAG algorithm:

   A single LMOD move inherently involves excessive bond
   stretching and bond angle bending in Cartesian space. Therefore the primarily
   torsional trajectory drawn by the low modes of vibration on the PES is
   severely contaminated by this naive, linear approximation and, therefore, the
   actual Cartesian LMOD trajectory often misses its target by climbing walls
   rather than crossing over into neighboring valleys at not too high altitudes.
   The ZIG-ZAG algorithm consists of a series of alternating short LMOD moves
   along the low-mode eigenvector (ZIG) followed by a few steps of minimization
   (ZAG), which is expected to relax excessive stretches and bends more than
   reversing the torsional move.  Therefore, it is expected that such a ZIG-ZAG
   trajectory will eventually be dominated by concerted torsional movements and
   will carry the molecule over the energy barrier in a way that is not too
   different from finding a saddle point and crossing over into the next valley
   like passing through a mountain pass.
*/
            /* calculate maximum distance any single
               atom would move along lmod_vec[]: */
            for (i = 0, max_atmov = ZERO; i < natm_local; i++) {
               for (j = 0, sum = ZERO; j < 3; j++)
                  sum += SQR(lmod_vec[i * 3 + j]);
               if (sum > max_atmov)
                  max_atmov = sum;
            }
            max_atmov = sqrt(max_atmov);
            /* pick random distance between user-defined min/max values: */
            lmod_step =
                *lmod_step_min + (*lmod_step_max -
                                  *lmod_step_min) * rand2();
            /* fastest moving atom moves exactly lmod_step distance unit: */
            scale = lmod_step / max_atmov;
            for (zigzag_iter = 0, energy_old = ref_energy, rms_old =
                 ZERO; zigzag_iter < max_zigzag_iter; zigzag_iter++) {
               for (i = 0; i < ndim_local; i++)
                  xyz_local[i] += scale * lmod_vec[i];      /* "ZIG" move */
               /* close pair separation operates on xyz_ext[]: */
               for (i=0; i<natm_local; i++) {  /* load xyz_ext[] */
                  j = atm_indx[i];
                  xyz_ext[ 3*j  ] = xyz_local[ 3*i  ];
                  xyz_ext[ 3*j+1] = xyz_local[ 3*i+1];
                  xyz_ext[ 3*j+2] = xyz_local[ 3*i+2];
               }
               do_all = TRUE;
               do_ligs = FALSE;
               separate_close_pairs(do_all, do_ligs, ndim_ext, xyz_ext,
                                    grad_ext, *nligs, lig_start, lig_end);
               for (i=0; i<natm_local; i++) {  /* load xyz_local[] */
                  j = atm_indx[i];
                  xyz_local[3*i  ] = xyz_ext[3*j  ];
                  xyz_local[3*i+1] = xyz_ext[3*j+1];
                  xyz_local[3*i+2] = xyz_ext[3*j+2];
               }
/*
                    Relax structure externally:
*/
               for (i=0; i<natm_local; i++) {  /* update xyz_ext[]
                                                  it is up-to-date, but it
                                                  is better to be explicit */
                  j = atm_indx[i];
                  xyz_ext[ 3*j  ] = xyz_local[ 3*i  ];
                  xyz_ext[ 3*j+1] = xyz_local[ 3*i+1];
                  xyz_ext[ 3*j+2] = xyz_local[ 3*i+2];
               }
               *return_flag = RELAX;                        /* "ZAG" move */
               extn_time_stamp = clock();
               *lmod_time += (double) (extn_time_stamp - lmod_time_stamp);
               *label = 5;
               return(0.0);
             L05:
               for (i=0; i<natm_local; i++) {  /* update xyz_local[]
                                                  w/ relaxed coords */
                  j = atm_indx[i];
                  xyz_local[3*i  ] = xyz_ext[3*j  ];
                  xyz_local[3*i+1] = xyz_ext[3*j+1];
                  xyz_local[3*i+2] = xyz_ext[3*j+2];
               }
               if (barrier_crossing_test_on) {
                  memcpy(xyz_hold, xyz_local, ndim_local * sizeof(double));
                  rms = rmsfit(ndim_local, xyz_cent, xyz_hold,
                               "RMS", &error_flag);
                  if (error_flag) {
                     if (allocated)
                        allocated = NO;
                     *label = error_flag;
                     goto error_cleanup;
                  }
/*
   IF the following Boolean holds:

   1.  The current endpoint of the zigzag trajectory is lower than ref_energy.
           - OR -
   2.  The endpoint is at least lower than it was in the previous
       zigzag iteration step.
           - AND -
       The molecule has also moved farther away from xyz_cent[]
       in terms of all-atom superposition RMS.

   THEN

       The LMOD ZIG-ZAG trajectory has crossed an energy barrier.  :-)
*/
                  if (*enrg < ref_energy
                      || (*enrg < energy_old && rms > rms_old))
                     break;     /* break out from zigzag */
                  else {
                     energy_old = *enrg;
                     rms_old = rms;
                  }
               }
            }                   /* end zigzag_iter */
/*
                Minimize structure after LMOD move:
*/
            for (i=0; i<natm_local; i++) {  /* update xyz_ext[]
                                               it is up-to-date, but it
                                               is better to be explicit */
               j = atm_indx[i];
               xyz_ext[ 3*j  ] = xyz_local[ 3*i  ];
               xyz_ext[ 3*j+1] = xyz_local[ 3*i+1];
               xyz_ext[ 3*j+2] = xyz_local[ 3*i+2];
            }
            *return_flag = MINIMIZE;
            extn_time_stamp = clock();
            *lmod_time += (double) (extn_time_stamp - lmod_time_stamp);
            *label = 6;
            return(0.0);
          L06:
            for (i=0; i<natm_local; i++) {  /* update xyz/grad_local[] */
               j = atm_indx[i];
               xyz_local[ 3*i  ] = xyz_ext[ 3*j  ];
               xyz_local[ 3*i+1] = xyz_ext[ 3*j+1];
               xyz_local[ 3*i+2] = xyz_ext[ 3*j+2];
               grad_local[3*i  ] = grad_ext[3*j  ];
               grad_local[3*i+1] = grad_ext[3*j+1];
               grad_local[3*i+2] = grad_ext[3*j+2];
            }
            energy = *enrg;
            /* assign Boltzmann factor: */
            for (i = 0, grad_rms = ZERO; i < ndim_local; i++)
               grad_rms += SQR(grad_local[i]);
            grad_rms = sqrt(grad_rms / ndim_local);
            if (grad_rms <= FULLY_MINIMIZED) {  /* structure fully minimized */
               if (energy <= ref_energy)
                  pboltz[kk++] = ONE;
               else
                  pboltz[kk++] = exp(-(energy - ref_energy) / (*rtemp));
            } else {
               pboltz[kk++] = ZERO;
            }
            rad = calc_rad(natm_local, xyz_local);
            /* archive structure in conflib */
            archive_structure(ndim_ext, xyz_ext, energy, rad, &error_flag);
            if (error_flag) {
               if (allocated)
                  allocated = NO;
               *label = error_flag;
               goto error_cleanup;
            }
            for (i = 0; i < ndim_local; i++)  /* store total LMOD ZIGZAG move */
               lmod_move[(kk - 1) * ndim_local + i] =
                                                 xyz_local[i] - xyz_cent[i];
            if (get_mytaskid() == 0) {
               if (DEBUG_LMOD) {
                  memcpy(xyz_hold, xyz_local, ndim_local * sizeof(double));
                  rms =
                      rmsfit(ndim_local, xyz_cent, xyz_hold,
                             "RMS", &error_flag);
                  if (error_flag) {
                     if (allocated)
                        allocated = NO;
                     *label = error_flag;
                     goto error_cleanup;
                  }
                  fprintf( nabout,
           "%3d  /%2d %cE =%13.3f (%6.3f)  Rg = %8.3f  rmsd=%7.3f  p=%7.4f\n",
           k - (*rotran) + 1, MIN(zigzag_iter + 1,
                                  max_zigzag_iter),
           (grad_rms <= FULLY_MINIMIZED) ? ' ' : '!', energy,
           grad_rms, calc_rad(natm_local, xyz_local), rms,
           pboltz[kk - 1]);
                  fflush(nabout);
               }
            }
/*
   Total quenching:
   Next point on LMOD trajectory will be
    the lowest energy neighbor of xyz_cent[].
   Quick quenching:
   Next point on LMOD trajectory will be
    the first minimum lower than xyz_cent[].
*/
            if (energy <= min_energy) {
               min_energy = energy;
               n = kk - 1;
               if ((*mc_option == QUICK_QUENCHING)
                   && (min_energy < ref_energy))
                  goto SKIP_METROPOLIS;
            }
            /* restore central structure: */
            memcpy(xyz_local, xyz_cent, ndim_local * sizeof(double));
            for (i = 0; i < ndim_local; i++)
               lmod_vec[i] = -lmod_vec[i];  /* reverse direction of eigvec  */
         }
      }                         /* end k vibrational modes */
/*
        Select next point on the LMOD trajectory using Metropolis criterion:
*/
      if (*mc_option == METROPOLIS) {
         cnt = 0;
         do {
            n = kk * rand2();
         } while (rand2() > pboltz[n] && ++cnt < MAX_MC_TRY);
      }
    SKIP_METROPOLIS:
      for (i = 0; i < ndim_local; i++)
         /* cross over barrier: */
         xyz_local[i] = xyz_cent[i] + lmod_move[n * ndim_local + i];
      for (i=0; i<natm_local; i++) {  /* update xyz_ext[] */
         j = atm_indx[i];
         xyz_ext[ 3*j  ] = xyz_local[ 3*i  ];
         xyz_ext[ 3*j+1] = xyz_local[ 3*i+1];
         xyz_ext[ 3*j+2] = xyz_local[ 3*i+2];
      }
      memcpy(trajectory + (lmod_iter + 1) * ndim_ext, xyz_ext,
             ndim_ext * sizeof(double));  /* add structure to LMOD trajectory */
/*
        Update archive and restart simulation periodically:
*/
      if ((lmod_iter % (*lmod_restart) == (*lmod_restart) - 1)
          || (lmod_iter % CONFLIB_UPDATE == CONFLIB_UPDATE - 1)) {
         /* update/cleanup conflib archive: */
         sizeof_archive = update_archive(ndim_ext, ndim_local, enrg_win, separ_rms, &error_flag);
         if (error_flag) {
            if (allocated)
               allocated = NO;
            *label = error_flag;
            goto error_cleanup;
         }
         if (get_mytaskid() == 0) {
            /* save archive to disk: */
            write_archive("conflib.dat", ndim_ext, &error_flag);
            if (error_flag) {
               if (allocated)
                  allocated = NO;
               *label = error_flag;
               goto error_cleanup;
            }
         }
         if (get_mytaskid() == 0) {
            if (PRINT_LMOD) {
               fprintf( nabout,
                   "__________________________________________________\n");
               fprintf( nabout, " Top Ten:\n");
               for (cnt = 0, walk = conflib_archive; walk != NULL;
                    walk = walk->next) {
                  fprintf(nabout, "   %5d  E =%13.3f /%2d  Rg = %8.3f\n", ++cnt,
                         walk->energy, walk->mult, walk->rad);
                  if (cnt == 10)
                     break;
               }
               fprintf( nabout,
                   "--------------------------------------------------\n");
               fprintf( nabout, " Total Nof confs = %5d\n", sizeof_archive );
               fprintf( nabout,
                   "________________________\n");
               fflush(nabout);
            }
         }
         if (lmod_iter % (*lmod_restart) == (*lmod_restart) - 1) {
            /* restart simulation with a low-energy structure: */
            restart_lmod(ndim_ext, MIN(*topten, sizeof_archive), xyz_ext);
            restart = YES;
            for (i=0; i<natm_local; i++) {  /* update xyz_local[]
                                               w/ restart coords */
               j = atm_indx[i];
               xyz_local[3*i  ] = xyz_ext[3*j  ];
               xyz_local[3*i+1] = xyz_ext[3*j+1];
               xyz_local[3*i+2] = xyz_ext[3*j+2];
            }
         }
      }
/*
        At every (*try_now)-th iteration apply explicit trans/rot to ligand(s):
*/
      if ((*nligs) && (lmod_iter < (*nlmodit) - 1)
          && !((lmod_iter + 1) % (*try_now))) {
         /* save xyz[] before applying trans/rot: */
         memcpy(xyz_hold, xyz_local, ndim_local * sizeof(double));
         for (k = 0, min_energy = BIG; k < (*ntry); k++) {
            for (kk = 0; kk < (*nligs); kk++) {
               do {
                  xtrans = 2 * rand2() - 1;
                  ytrans = 2 * rand2() - 1;
                  ztrans = 2 * rand2() - 1;
               } while ((trnorm =
                         sqrt(SQR(xtrans) + SQR(ytrans) +
                              SQR(ztrans))) < TINY);
               trscale = trmin[kk] + (trmax[kk] - trmin[kk]) * rand2();
               xtrans *= trscale / trnorm;
               ytrans *= trscale / trnorm;
               ztrans *= trscale / trnorm;
               for (ikk = 0; ikk < natm_local; ikk++ ) {
                 if ( atm_indx[ikk] == (lig_start[kk] -1) ) {  /* lig_start[] contains external PDB atom numbers */
                   lig_start_kk = ikk;
                   break;
                 }
               }
               for (ikk = 0; ikk < natm_local; ikk++ ) {
                 if ( atm_indx[ikk] == (lig_end[kk] -1) ) {  /* lig_end[] contains external PDB atom numbers */
                   lig_end_kk = ikk;
                   break;
                 }
               }
               trans_ligand(xyz_local, lig_start_kk, lig_end_kk, xtrans,
                            ytrans, ztrans);
               do {
                  xrot = 2 * rand2() - 1;
                  yrot = 2 * rand2() - 1;
                  zrot = 2 * rand2() - 1;
               } while ((rotnorm =
                         sqrt(SQR(xrot) + SQR(yrot) + SQR(zrot))) < TINY);
               xrot /= rotnorm;
               yrot /= rotnorm;
               zrot /= rotnorm;
               rotang = angmin[kk] + (angmax[kk] - angmin[kk]) * rand2();
               calc_rot_matrix(DEG2RAD * rotang, xrot, yrot, zrot, rotmat);
               if (lig_rot_cent[kk]) {
                 for (ikk = 0; ikk < natm_local; ikk++ ) {
                   if ( atm_indx[ikk] == (lig_rot_cent[kk] -1)   ) {  /* lig_rot_cent[] contains external PDB atom numbers */
                     lig_rot_cent_kk = ikk;
                     break;
                   }
                 }
                 xcent = xyz_local[ lig_rot_cent_kk * 3    ];
                 ycent = xyz_local[ lig_rot_cent_kk * 3 + 1];
                 zcent = xyz_local[ lig_rot_cent_kk * 3 + 2];
               } else
                 calc_centroid(xyz_local, lig_start_kk, lig_end_kk,
                                &xcent, &ycent, &zcent);
               rot_ligand(xyz_local, lig_start_kk, lig_end_kk, xcent,
                          ycent, zcent, rotmat);
            }
            /* close pair separation operates on xyz_ext[]: */
            for (i=0; i<natm_local; i++) {  /* update xyz_ext[] */
               j = atm_indx[i];
               xyz_ext[ 3*j  ] = xyz_local[ 3*i  ];
               xyz_ext[ 3*j+1] = xyz_local[ 3*i+1];
               xyz_ext[ 3*j+2] = xyz_local[ 3*i+2];
            }
            do_all = FALSE;
            do_ligs = TRUE;
            separate_close_pairs(do_all, do_ligs, ndim_ext, xyz_ext,
                                 grad_ext, *nligs, lig_start, lig_end);
            for (i=0; i<natm_local; i++) { /* update xyz_local[] */
               j = atm_indx[i];
               xyz_local[3*i  ] = xyz_ext[3*j  ];
               xyz_local[3*i+1] = xyz_ext[3*j+1];
               xyz_local[3*i+2] = xyz_ext[3*j+2];
            }
/*
            Minimize structure after trans/rot move:
*/
            for (i=0; i<natm_local; i++) {  /* update xyz_ext[]
                                               it is up-to-date, but it
                                               is better to be explicit */
               j = atm_indx[i];
               xyz_ext[ 3*j  ] = xyz_local[ 3*i  ];
               xyz_ext[ 3*j+1] = xyz_local[ 3*i+1];
               xyz_ext[ 3*j+2] = xyz_local[ 3*i+2];
            }
            *return_flag = MINIMIZE;
            extn_time_stamp = clock();
            *lmod_time += (double) (extn_time_stamp - lmod_time_stamp);
            *label = 7;
            return(0.0);
          L07:
            for (i=0; i<natm_local; i++) {  /* update xyz_local[]
                                               strictly not necessary,
                                               see memcpy below */
               j = atm_indx[i];
               xyz_local[3*i  ] = xyz_ext[3*j  ];
               xyz_local[3*i+1] = xyz_ext[3*j+1];
               xyz_local[3*i+2] = xyz_ext[3*j+2];
            }
/*
                Save lowest energy trans/rot-ed ligand(s)' configuration:
*/
            if ((*enrg) <= min_energy) {
               min_energy = (*enrg);
               for (i = 0; i < ndim_local; i++)
                  trot_move[i] = xyz_local[i] - xyz_hold[i];
            }
            /* restore xyz[] before another trans/rot: */
            memcpy(xyz_local, xyz_hold, ndim_local * sizeof(double));
         }                      /* end trans/rot loop */
         for (i = 0; i < ndim_local; i++)
            /* apply best trans/rot move: */
            xyz_local[i] = xyz_hold[i] + trot_move[i];
      }                         /* end nligs  iter */
   }                            /* end lmod() iter */
/*
    Do final update/cleanup on conflib archive:
*/
   sizeof_archive = update_archive(ndim_ext, ndim_local, enrg_win, separ_rms, &error_flag);
   if (error_flag) {
      if (allocated)
         allocated = NO;
      *label = error_flag;
      goto error_cleanup;
   }
   if (get_mytaskid() == 0) {
      write_archive("conflib.dat", ndim_ext, &error_flag);
      if (error_flag) {
         if (allocated)
            allocated = NO;
         *label = error_flag;
         goto error_cleanup;
      }
   }
   if (get_mytaskid() == 0) {
      if (PRINT_LMOD || DEBUG_LMOD) {
         fprintf(nabout, "__________________________________________________\n");
         if (!DEBUG_LMOD)
            fprintf(nabout, " Top Ten:\n");
         else
            fprintf(nabout, " Full list:\n");
         for (cnt = 0, walk = conflib_archive; walk != NULL;
              walk = walk->next) {
            fprintf(nabout, "   %5d  E =%13.3f /%2d  Rg = %8.3f\n", ++cnt,
                   walk->energy, walk->mult, walk->rad);
            if (!DEBUG_LMOD && cnt == 10)
               break;
            if (DEBUG_LMOD && cnt == *nconf)
               break;
         }
         fprintf( nabout,
             "--------------------------------------------------\n");
         fprintf( nabout, " Total Nof confs = %5d\n", sizeof_archive );
         fprintf( nabout,
                   "________________________\n\n");
         fflush(nabout);
      }
   }
/*
    Load archive into conflib[] for use in parent program:
*/
   load_archive(nconf, ndim_ext, conflib, &error_flag);
   if (error_flag) {
      if (allocated)
         allocated = NO;
      *label = error_flag;
      goto error_cleanup;
   }
   glob_min_energy = conflib_archive->energy;
   /* update outgoing xyz_ext[] w/ glob min coords: */
   memcpy(xyz_ext, conflib_archive->xyz, ndim_ext * sizeof(double));
/*
    LMOD done, return to parent program:
*/
   if (allocated) {
      my_free(eigvals);
      my_free(eigvecs);
      my_free(index);
      my_free(xyz_cent);
      my_free(lmod_vec);
      my_free(pboltz);
      my_free(xyz_hold);
      my_free(lmod_move);
      my_free(trot_move);
      my_free(xyz_local);
      my_free(grad_local);
      my_free(atm_indx);
      allocated = NO;
   }
   if (conflib_archive != NULL)
      destroy_archive();
   *lmod_time += (double) (clock() - lmod_time_stamp);
   *lmod_time /= CLOCKS_PER_SEC;
   *extn_time /= CLOCKS_PER_SEC;
   *return_flag = DONE;
   *label = 0;                  /* lmod() done */
   return glob_min_energy;

 error_cleanup:
   my_free(eigvals);
   my_free(eigvecs);
   my_free(index);
   my_free(xyz_cent);
   my_free(lmod_vec);
   my_free(pboltz);
   my_free(xyz_hold);
   my_free(lmod_move);
   my_free(trot_move);
   my_free(xyz_local);
   my_free(grad_local);
   my_free(atm_indx);
   if (conflib_archive != NULL)
      destroy_archive();
   return(0.0);
}

#ifndef SQM

/*
    The following used to be in ../nab/lmod.nab, but here is translated to
    C to further separate the nab language from the sff routines.

    For now, the function to be used continues to be mme(); in the future,
    it could be passed as a pointer, as in xminC.c.  But the need to run
    lmod optimization on a generic function is less clear than for
    minimization itself.
*/

//       libLMOD reverse communication molecular simulation package.
//                      Written by Istvan Kolossvary.

INT_T   lmod_opt_init( struct lmod_opt *lo, struct xmin_opt *xo )
{
        lo->minim_grms          = 0.1;
        lo->niter               = 10;
        lo->nmod                = 5;
        lo->kmod                = 3;
        lo->nrotran_dof         = 6;
        lo->nconf               = 10;
        lo->energy_window       = 50.0;
        lo->conf_separation_rms = 2.0;
        lo->eig_recalc          = 5;
        lo->ndim_arnoldi        = 0;
        lo->lmod_restart        = 10;
        lo->n_best_struct       = 10;
        lo->mc_option           = 1;
        lo->rtemp               = 1.5;
        lo->lmod_step_size_min  = 2.0;
        lo->lmod_step_size_max  = 5.0;
        lo->nof_lmod_steps      = 0;
        lo->lmod_relax_grms     = 1.0;
        lo->nlig                = 0;
        lo->apply_rigdock       = 2;
        lo->nof_poses_to_try    = 10;
        lo->random_seed         = 314159;
        lo->print_level         = 0;

        xo->mol_struct_opt      = 1;
        xo->maxiter             = 20000;
/*      xo->grms_tol            = The default is handled within lmod() below
                                 using lo.minim_grms and lo.lmod_relax_grms. */
        xo->method              = 2;
        xo->numdiff             = 1;
        xo->m_lbfgs             = 3;
        xo->ls_method           = 2;
        xo->ls_maxiter          = 20;
        xo->ls_maxatmov         = 0.5;
        xo->beta_armijo         = 0.5;
        xo->c_armijo            = 0.4;
        xo->mu_armijo           = 1.0;
        xo->ftol_wolfe          = 0.0001;
        xo->gtol_wolfe          = 0.9;
        xo->print_level         = 0;
        return(0);
}

REAL_T lmod(INT_T *natm, REAL_T *xyz, REAL_T *grad, REAL_T *energy,
           REAL_T *conflib, REAL_T *lmod_trajectory,
           INT_T *lig_start, INT_T *lig_end, INT_T *lig_cent, 
           REAL_T *tr_min, REAL_T *tr_max, REAL_T *rot_min, REAL_T *rot_max,
           struct xmin_opt *xo, struct lmod_opt *lo )

{

/*
 *  This package can be used with any program that can calculate the energy
 *  and the gradient for a particular [x, y, z] atomic configuration. There is
 *  no size limit, but the xyz[], grad[], conflib[], lmod_trajectory[],
 *  lig_start[], lig_end[], lig_cent[], tr_min[], tr_max[], rot_min[], and
 *  rot_max[] arrays must be allocated by the calling program.
 *
 *  Input params:  Number of atoms, xyz[] and grad[] arrays, conflib[] and
 *  lmod_trajectory[] arrays, number of LMOD simulation steps, number of
 *  low-frequency vibrational modes used by LMOD, number of conformations to
 *  be stored within a user-specified energy window, the number of ligands to
 *  be docked with their identification lig_start[] and lig_end[] as well as
 *  rotran parameters tr_min[], tr_max[], rot_min[] and rot_max[] to specify
 *  how the ligand(s) move around during the simulation (see details below).
 *  
 *  Output params:  CPU time spent in the lmod() routine, and possibly an
 *  error message (see below).
 *  
 *  On exit, the low-energy conflib library will be loaded into the conflib[]
 *  array and the LMOD simulation path (trajectory) is loaded into
 *  lmod_trajectory[]. The conflib libray is also stored and periodically
 *  updated in a binary file called "conflib.dat". Each block in conflib.dat
 *  represents a single conformation sorted by increasing energy (glob. min.
 *  first). A single block consists of a header followed by the x, y, z
 *  coordinates of the atoms. The header holds three numbers: float (8 bytes)
 *  energy of the conformation, float (8 bytes) radius of gyration and int (4
 *  bytes) multiplicity, i.e., the number of times that particular
 *  conformation was found during the LMOD simulation.  The return value of
 *  lmod() is the global minimum energy and the glob. min. structure is
 *  loaded into xyz[].
 *
 *
 *  Params:
 *  
 *  natm            Number of atoms included in the calculation.
 *
 *  xyz             Allocated array of (x, y, z) atomic coordinates (dimension
 *                  = 3*natm).
 *
 *  grad            Allocated array of the gradient (dimension = 3*natm).
 *
 *  minim_grms      RMS gradient convergence criterion for minimization.
 *                  Should be <= 0.1.
 *
 *  energy          Energy value.
 *
 *  niter           Number of LMOD iterations.
 *
 *  nmod            Total number of low-frequency vibrational modes utilized
 *                  by LMOD.
 *
 *  kmod            Number of modes out of nmod to be explored in each LMOD
 *                  iteration.
 *
 *  nrotran_dof     Number of external/trivial/rotranslational degrees of
 *                  freedom.  (Allowed values depending on the number of
 *                  frozen/tethered atoms in the system: zero atoms dof=6, 1
 *                  atom dof=3, 2 atoms dof=1, >=3 atoms dof=0.)
 *
 *  nconf           Number of conformations or docking modes/poses to be
 *                  stored in conflib[].
 *
 *  energy_window   Max energy gap above global minimum for stored structures.
 *                  Energy unit defined in calling program.
 *
 *  conf_separation_rms  Min pairwise all-atom RMSD separation for
 *                       stored structures.
 *
 *  conflib         Allocated array of (x, y, z) atomic coordinates of conflib
 *                  structures (dimension = nconf).
 *
 *  lmod_trajectory Allocated array of (x, y, z) atomic coordinates of
 *                  consecutive structures visited along the LMOD
 *                  path/trajectory (dimension = niter+1).
 *
 *  eig_recalc      Frequency by which lmod() recalculates the low-mode
 *                  eigenvectors. eig_recalc=1 means that eigenvectors are
 *                  recalculated at every LMOD iteration, eig_recalc=5 means
 *                  they are only recalculated at every five iterations, etc.
 *
 *  ndim_arnoldi    This is the dimension of the Arnoldi factorization.
 *                  Basically, the ARPACK package used for the eigenvector
 *                  calculations solves multiple "small" eigenproblems instead
 *                  of a single "large" problem, which is the diagonalization
 *                  of the 3xnatm by 3xnatm Hessian matrix. This parameter is
 *                  the user- specified dimension of the "small" problem.
 *                  Allowed range is nmod+1 <= ndim_arnoldi <= 3xnatm.
 *                  ndim_arnoldi=0 translates to ndim_arnoldi=3xnatm which
 *                  means that the "small" problem and the "large" problem are
 *                  identical. This is the preferred/fastest calculation for
 *                  samll to medium size systems, because ARPACK is guaranteed
 *                  to converge in a single iteration. The ARPACK calculation
 *                  scales with 3xnatm x ndim_arnoldi^2 and, therefore, for
 *                  larger molecules there is an optimal ndim_arnoldi <<
 *                  3xnatm that converges much faster in multiple iterations
 *                  (possibly thousands or tens of thousands of iterations).
 *                  For proteins, ndim_arnoldi=1000 is generally a good value.
 *
 *  lmod_restart    Frequency by which conflib[] is updated and LMOD
 *                  simulation is restarted with a randomly chosen structure
 *                  among the n_best_struct lowest energy structures found so
 *                  far. A good value for lmod_restart is 10, which means that
 *                  conflib[] is updated (and written to the binary file
 *                  conflib.dat) after every 10th LMOD iteration and the
 *                  simulation is restarted with one of the lowest-energy
 *                  structures.  If lmod_restart >= niter (number of LMOD
 *                  iterations), the simulation will never restart.
 *
 *  n_best_struct   Number of the lowest-energy structures found so far at a
 *                  particular LMOD restart point.  The structure to be used
 *                  for the restart will be chosen randomly from this pool.
 *                  n_best_struct = 10 is generally a good choice.
 *                  n_best_struct = 1 allows the user to explore the
 *                  neighborhood of the then current global minimum.
 *
 *  mc_option       Monte Carlo option. Allowed values: '1' Metropolis Monte
 *                  Carlo, '2' "total quenching" (the LMOD trajectory always
 *                  proceeds towards the lowest lying neighbor of a particular
 *                  energy well found after exhaustive search along all of the
 *                  low modes), and '3' "quick quenching" (the LMOD trajectory
 *                  proceeds towards the first neighbor found, which is lower
 *                  in energy than the current point on the path, without
 *                  exploring the remaining modes).
 *
 *  rtemp           The value of RT at a particular, user-defined temperature,
 *                  given in the energy unit used in the calling program.
 *                  rtemp is utilized in the Metropolis criterion.
 *
 *  The basic tenet of LMOD is climbing energy barriers with ease. In libLMOD
 *  this is done by utilizing the LMOD ZIG-ZAG algorithm:
 *
 *  A single LMOD move inherently involves excessive bond stretching and bond
 *  angle bending in Cartesian space. Therefore the primarily torsional
 *  trajectory drawn by the low modes of vibration on the PES is severely
 *  contaminated by this naive, linear approximation and, therefore, the
 *  actual Cartesian LMOD trajectory often misses its target by climbing walls
 *  rather than crossing over into neighboring valleys at not too high
 *  altitudes.  The ZIG-ZAG algorithm consists of a series of alternating
 *  short LMOD moves along the low-mode eigenvector (ZIG) followed by a few
 *  steps of minimization (ZAG), which is expected to relax excessive
 *  stretches and bends more than reversing the torsional move.  Therefore, it
 *  is expected that such a ZIG-ZAG trajectory will eventually be dominated by
 *  concerted torsional movements and will carry the molecule over the energy
 *  barrier in a way that is not too different from finding a saddle point and
 *  crossing over into the next valley like passing through a mountain pass.
 *                  
 *  lmod_step_size_min  Minimum lengths of a single LMOD ZIG move specified in
 *                      the distance unit used in the calling program. A
 *                      generally good value is 2 Angs.
 *
 *  lmod_step_size_max  Maximum lengths of a single LMOD ZIG move specified 
 *                      in the distance unit used in the calling program. A
 *                      generally good value is 5 Angs.  The actual length of
 *                      the ZIG move will be chosen randomly between the min
 *                      and max values.
 *
 *  nof_lmod_steps  Total number of ZIG-ZAG moves. nof_lmod_steps=0 means that
 *                  the number of ZIG-ZAG moves is not pre-defined, instead
 *                  LMOD will attempt to cross the barrier in as many ZIG-ZAG
 *                  moves as it is necessary. The criterion of crossing a
 *                  barrier is as follows:
 *
 *                  IF the following Boolean holds:
 *
 *                  1.  The current endpoint of the zigzag trajectory is lower
 *                      than the bottom of the current energy well from where
 *                      the zigzag trajectory starts.
 *
 *                              - OR -
 *
 *                  2.  The endpoint is at least lower than it was in the
 *                       previous zigzag iteration step.
 *                              - AND -
 *                       The molecule has also moved farther away from the
 *                       starting point in terms of all-atom superposition
 *                       RMS.
 *
 *                  THEN
 *
 *                  The LMOD ZIG-ZAG trajectory has crossed an energy 
 *                  barrier :-)
 *
 *                  nof_lmod_steps > 0 means that multiple barriers may be
 *                  crossed and LMOD can carry the molecule to a large
 *                  distance on the PES without severely distorting the
 *                  geometry.
 *
 *  lmod_relax_grms This is the endpoint criterion of the ZIG
 *                   relaxation/minimization in terms of gradient RMS.
 *                   lmod_relax_grms = 1 is generally a good value for
 *                   crossing the closest barrier. For larger moves a less
 *                   stringent criterion can be adequate.
 *
 *  nlig            Number of ligands considered for flexible docking.
 *
 *  lig_start/      It is assumed that ligand atoms form a consecutive list
 *  lig_end         (no gaps). lig_start[i] is the lowest atom number in
 *                  ligand 'i' and lig_end[i] is the highest atom number in
 *                  ligand 'i'.  lig_start[] and lig_end[] must be allocated
 *                  in the calling program (dimension = nlig).
 *
 *  lig_cent        lig_cent[i] defines the center of rotation of ligand 'i'.
 *                  lig_start[i] <= lig_cent[i] <= lig_end[i] specifies a
 *                  particular atom as the center of rotation whereas
 *                  lig_cent[i] = 0 means that the center of rotation will be
 *                  the center of gravity (geometric centroid) of the ligand.
 *                  lig_cent[] should also be allocated by the calling program
 *                  to nlig dimensions.
 *
 *  apply_rigdock   Frequency by which lmod() applies "nof_poses_to_try"
 *                  rigid-body explicit trans/rot to the ligand(s).
 *                  apply_rigdock=1 means that such trans/rot takes place at
 *                  every LMOD iteration, apply_rigdock=5 means trans/rot is
 *                  applied at every five iterations, etc.
 *
 *  nof_poses_to_try Number of explicit rotations/translations applied to the
 *                   ligand(s) after each apply_rigdock-th LMOD iteration.
 *
 *  tr_min/         Range of random translation of ligand 'i' between a
 *  tr_max          minimum distance tr_min[i] and a maximum distance
 *                  tr_max[i]. The distance must be specified in the distance
 *                  units used in the calling program. tr_min[] and tr_max[]
 *                  must be allocated in the calling program (dimension = nlig).
 *
 *  rot_min/        Range of random rotation of ligand 'i' between a minimum
 *  rot_max         angle rot_min[i] and a maximum angle rot_max[i] about the
 *                  origin specified by lig_cent[i]. The angle is given in
 *                  degrees.  rot_min[] and rot_max[] also must be allocated
 *                  in the calling program (dimension = nlig).
 *
 *  random_seed     Random seed used to reproduce previous runs. random_seed =
 *                  0 uses a hardware seed.
 *
 *  print_level     Verbosity: 0= none, 1= some details, 2= more details, 3=
 *                  everything incl. ARPACK.
 *
 *  lmod_time       Total time spent in the lmod() routine in CPU sec.
 *
 *  aux_time        Total time spent in the calling program for calculating
 *                  energies/gradients and do energy minimizations when
 *                  requested by lmod().
 *
 *  XMIN parameters controllable through lmod:
 *  xo.maxiter    Max number of minimization steps.
 *  xo.ls_maxiter Max number of line search steps per minimization step.
 *  xo.method     Minimization method: 1= PRCG, 2= LBFGS, 
 *                                     3= LBFGS-preconditioned TNCG.
 *  xo.ls_method  Line search method:  1= modified Armijo.
 *                                     2= Wolfe (J. J. More', D. J. Thuente).
 *  xo.numdiff    Method used in finite difference Hv matrix-vector products:
 *                  1= forward difference, 2= central difference.
 *  xo.m_lbfgs    Depth of LBFGS memory for LBFGS minimization or TNCG 
 *                  preconditioning. Suggested value 3.
 *                  m_lbfgs=0 with TNCG minimization turns off preconditioning.
 *  xo.ls_maxatmov Max atomic coord movement allowed in line search, range > 0.
 *  xo.beta_armijo Armijo beta param, range (0, 1).
 *  xo.c_armijo    Armijo c param,    range (0, 0.5 ).
 *  xo.mu_armijo   Armijo mu param,   range (0, 2).
 *  xo.ftol_wolfe  Wolfe ftol param,  range (0, 0.5 ).
 *  xo.gtol_wolfe  Wolfe gtol param,  range (xo.ftol_wolfe, 1 ).
 *  xo.print_level  Verbosity: 0= none, 1= minim details, 2= minim and line
 *                             search details plus CG details in TNCG.
 *  error_flag      Error flag. lmod() will print a descriptive error message.
 *
 */

	INT_T return_flag, status_flag;
	REAL_T glob_min_energy, xmin_grms;
    INT_T NEG_TWO = -2, NEG_FOUR = -4;
    REAL_T mme( REAL_T*, REAL_T*, INT_T*);

//  Call LMOD:
	for (status_flag = lo->error_flag = 0;;) {
		glob_min_energy = lmodC(&(lo->niter), &(lo->nmod), &(lo->kmod),
               &(lo->nrotran_dof), natm, xyz, energy, grad,
			   &(lo->nconf), &(lo->energy_window), &(lo->conf_separation_rms),
			   conflib,lmod_trajectory, &(lo->eig_recalc), &(lo->ndim_arnoldi),
			   &(lo->lmod_restart), &(lo->n_best_struct), &(lo->mc_option),
			   &(lo->rtemp), &(lo->lmod_step_size_min),
			   &(lo->lmod_step_size_max), &(lo->nof_lmod_steps), &(lo->nlig),
			   lig_start, lig_end, &(lo->apply_rigdock),
			   &(lo->nof_poses_to_try), tr_min, tr_max, lig_cent,
			   rot_min, rot_max, &(lo->random_seed), &(lo->print_level),
			   &(lo->lmod_time), &(lo->aux_time), &return_flag,
			   &status_flag);

//      Finished LMOD simulation:
		if (return_flag == 0) {
			if (status_flag >= 0)
				lo->error_flag = 0;
			else
				lo->error_flag = status_flag;
		}

//      Force NB list update by passing '-4' to mme():
		else if (return_flag == 1 || return_flag == 2) {
			if (status_flag >= 0)
				*energy = mme(xyz, grad, &NEG_FOUR);
			else
				lo->error_flag = status_flag;
		}

//      Prevent NB list update by passing '-2' to mme():
		else if (return_flag == 3 || return_flag == 4) {
			if (status_flag >= 0)
				*energy = mme(xyz, grad, &NEG_TWO);
			else
				lo->error_flag = status_flag;
		}

//      Minimize current xyz[] structure externally:
		else if (return_flag == 5) {
			if (status_flag >= 0) {
                xo->grms_tol = lo->minim_grms;
				*energy = xmin(&mme, natm, xyz, grad, energy,
                                              &xmin_grms, xo );
			} else
				lo->error_flag = status_flag;
		}

//      Relax current xyz[] structure externally:
		else if (return_flag == 6) {
			if (status_flag >= 0) {
                xo->grms_tol = lo->lmod_relax_grms;
				*energy = xmin(&mme, natm, xyz, grad, energy,
                                              &xmin_grms, xo );
			} else
				lo->error_flag = status_flag;
		} else {
			fprintf(stderr, "\n LMOD ERROR: return_flag corrupted.\n");
			lo->error_flag = -100;
			return (0.0);
		}
		if (lo->error_flag || !return_flag)
			break;
	}
	if (lo->error_flag) {
		fprintf(stderr, "\n LMOD ERROR: %d\n", status_flag);
		return (0.0);
	}
	return (glob_min_energy);
}
#endif
