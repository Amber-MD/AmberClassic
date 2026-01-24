/*
 *
 *	Collection of routines to:
 *		constuct a distance matrix from the bounds matrix;
 *      create the metric matrix
 *      diagonalize this and constuct the coordinates
 *
 */
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define NDEBUG

#include "nab.h"
#include "errormsg.h"
#include "memutil.h"
#include "traceback.h"
#include "chirvol.h"

#define	B_NAME_SIZE	32
typedef struct bdbent {
    char b_name[B_NAME_SIZE];
    REAL_T b_mean;
    REAL_T b_stddev;
    REAL_T b_scale;
    REAL_T b_min;
    REAL_T b_max;
} BDBENT;

#define	IBSIZE		40
#define	JBSIZE		10

#define	CLOSEST		0.5
#define	FURTHEST	10000.0

#define	EPSILON		1e-8

#define	NEARZERO(a)	( fabs(a) < 10.0 * EPSILON ? 0.0 : (a) )

    /* molecule stuff:      */
static int natoms;

    /* varaibles set in options string:         */

static int eamax = 10;          /* maximum number of embed attempts */
static REAL_T randpair = 0.0;   /* use randomized pairwise metrization  */
static int ddm = 0;             /* dump distance matrix once created */
static int rdm = 0;             /* read distance matrix */
static int dmm = 0;             /* dump metric matrix once created */
static int rmm = 0;             /* read metric matrix */
static int mxp = 0;             /* use aexp for metrization */
static int majorize_opt = 0;    /* use Guttman transform majorization */
static char ddmn[256];          /* dump distance matrix filename  */
static char rdmn[256];          /* read distance matrix filename  */
static char dmmn[256];          /* dump metric matrix filename  */
static char rmmn[256];          /* read metric matrix filename */
static char mxpn[256];          /* initial metrization atom expression */
static int gdist = 0;           /* if true, use Gaussian distances  */
static REAL_T k4d = 1.0;        /* force constant for squeezing out 4th d */
static int ntpr = 10;           /* print frequency in db_viol  */
static REAL_T kchi = 1.0;       /* chirality weight     */
static int niter = 20;          /* number of iterations in majorization */
static REAL_T lbpen = 3.5;      /* penalty in db_viol for lower bounds viol. */
static int sqviol = 0;          /* use straight parabolas for db_viol   */
static int pembed = 0;          /* use embeding scheme described in de Groot,
                                   et al., Proteins 29:240-251 (1997)   */
static int shuffle = 1;         /* randomize starting coords for pembed scheme
                                   inside a box 0..rbox, 0..rbox, 0..rbox  */
static REAL_T rbox = 20.;       /* box dimension for random embed    */
static int riter = 1000;        /* max. interations for random pembed  */
static REAL_T slearn = 1.0;     /* starting learning parameter for pembed  */

    /* other DG variables:              */

static REAL_T **bmat;           /* bounds matrix    */
static REAL_T *dmat;            /* distance matrix, in linear form  */
static REAL_T **vec;            /* biggest three eigenvectors */
static REAL_T **mmat;           /* metric matrix    */
static REAL_T *com2;            /* squared com distances */

static REAL_T *e_vals;
static REAL_T *e_work;          /* scratch diag array used by diagonalization */
static int nchi;                /* number of chi ctrs   */
static CHIRAL_T *chi;           /* chirality tetrads    */
static REAL_T mostneg;          /* most negative diagonal element of metric matrix */

static int seed = -1;
static REAL_T pencut = -1.0;    /* if >= zero, print energy contributions */

REAL_T rand2();
REAL_T gauss2(REAL_T *, REAL_T *);
REAL_T vdw_radius( BOUNDS_T *, int );

    /* the user interface */
int dg_options(BOUNDS_T *, char *);
int embed(BOUNDS_T *, REAL_T *);
int metrize(BOUNDS_T * bp, int nmetrize);
static int findPaths(BOUNDS_T *, int);
static int makeBoundList(BOUNDS_T *);
int geodesics(BOUNDS_T *);
int tsmooth(BOUNDS_T *, REAL_T);
REAL_T db_viol(REAL_T *, REAL_T *, int *);
REAL_T db_viol3(REAL_T *, REAL_T *, int *);

    /* routines used for distance geom:     */

static void mk_dmat(BOUNDS_T *, REAL_T *, REAL_T(*)());
static int mk_mmat(BOUNDS_T *, REAL_T *, REAL_T **);
static int rd_mmat(char[], BOUNDS_T *, REAL_T **);
static void trifix(int, int, int, REAL_T **);

    /* create distances */

static REAL_T u_rdist(REAL_T, REAL_T);
static REAL_T g_rdist(REAL_T, REAL_T);

static void dumpmat(FILE *, BOUNDS_T *, REAL_T **);
void mkdgname(BOUNDS_T *, int, char[]);

#include "lex.dgo.c"


int embed(BOUNDS_T * bp, REAL_T * nxyz)
{
    int rval = 0;
    int i, j, k, iter, miter, nea;
    int iters[3];
    REAL_T d1, d2;
    REAL_T f, r, lb, ub;
    REAL_T(*rdist) ();
    REAL_T *grad;
    REAL_T dx, dy, dz, fiter, friter, learn;

    rdist = u_rdist;
    if (gdist)
        rdist = g_rdist;

    if (pembed) {

/*  ---- Use radomized embed scheme of de Groot et al.  -----------------   */

        if ((grad = vector(0, 4 * natoms - 1)) == NULL) {
            rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "grad");
            exit(1);
        }

        if (shuffle) {
            /*  first, set up random coordinates  */

            for (i = 0; i < natoms; i++) {
                nxyz[4 * i] = rbox * rand2();
                nxyz[4 * i + 1] = rbox * rand2();
                nxyz[4 * i + 2] = rbox * rand2();
                nxyz[4 * i + 3] = 0.002 * rand2() - 0.001;
            }
        }

        /*  Now, randomly choose pairs of atoms, and adjust their positions:  */

        learn = slearn;
        friter = riter;
        for (iter = 1; iter <= riter; iter++) {
            for (miter = 0; miter < 1000; miter++) {
                /* count interations in K  */
                i = natoms * rand2();
                j = natoms * rand2();
                if (i == j)
                    continue;
                if (i > j) {
                    k = i;
                    i = j;
                    j = k;
                }               /* (now i<j)  */
                lb = bmat[j][i];
                ub = bmat[i][j];

                dx = nxyz[4 * j] - nxyz[4 * i];
                dy = nxyz[4 * j + 1] - nxyz[4 * i + 1];
                dz = nxyz[4 * j + 2] - nxyz[4 * i + 2];
                r = sqrt(dx * dx + dy * dy + dz * dz);
                if (r < lb || r > ub) {

                    /* get the new distance: the original de Groot et al idea
                       places this randomly between the upper and lower bound;
                       Agrafiotis (J. Comput. Chem. 24:1215, 2003), in a 
                       slighly different context, suggested damping the move
                       by a learning parameter, which decreases from 1 to 0.  */

                    d1 = rdist(lb, ub);
                    d2 = learn * (d1 - r) / (2. * r);

                    nxyz[4 * j] += d2 * dx;
                    nxyz[4 * j + 1] += d2 * dy;
                    nxyz[4 * j + 2] += d2 * dz;
                    nxyz[4 * i] -= d2 * dx;
                    nxyz[4 * i + 1] -= d2 * dy;
                    nxyz[4 * i + 2] -= d2 * dz;
                }
            }

            /*  "learning" parameter: goes from 1 to 0  */
            fiter = iter;
            learn = slearn * (1. - 0.95 * fiter / friter);
            if (iter == 1 || iter % ntpr == 0)
                f = db_viol(nxyz, grad, &iter);
        }
        free_vector(grad, 0, natoms - 1);

    } else {

/*  ---- Diagonalize the metric matrix to get coords.   -----------------   */

        /*   allocate memory   */

        if ((dmat = vector(0, natoms * (natoms - 1) / 2)) == NULL) {
            rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "dmat");
            exit(1);
        }
        if ((mmat = matrix(0, natoms - 1, 0, natoms - 1)) == NULL) {
            rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "mmat");
            exit(1);
        }
        if ((vec = matrix(0, natoms - 1, 0, 2)) == NULL) {
            rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "mmat");
            exit(1);
        }
        if ((com2 = vector(0, natoms - 1)) == NULL) {
            rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "com2");
            exit(1);
        }
        if ((e_vals = vector(0, natoms - 1)) == NULL) {
            rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "e_vals");
            exit(1);
        }
        if ((e_work = vector(0, natoms - 1)) == NULL) {
            rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "e_work");
            exit(1);
        }

        for (nea = 0; nea < eamax; nea++) {

            printf("embed attempt %d\n", nea + 1);

            if (rmm)
                rd_mmat(rmmn, bp, mmat);
            else {
                mk_dmat(bp, dmat, rdist);
                mk_mmat(bp, dmat, mmat);
            }

            /*  embed the metric matrix;  */

            /*  Use power method with deflation    */

            miter = 400;
            for (j = 0; j < 3; j++) {
                for (i = 0; i < natoms; i++) {
                    vec[i][j] = 1.0;
                }
            }

            for (k = 0; k < 3; k++) {
                e_vals[k] = d1 = 0.0;
                for (i = 0; i < natoms; i++) {
                    e_work[i] = 0.0;
                    for (j = 0; j < natoms; j++)
                        e_work[i] += mmat[i][j] * vec[j][k];
                    d1 += e_work[i] * e_work[i];
                }

                if (d1 < natoms * 1.e-5) {
                    for (i = 0; i < natoms; i++)
                        e_work[i] = rand2();
                }

                for (iter = 0; iter < miter; iter++) {
                    d1 = d2 = 0.0;
                    for (i = 0; i < natoms; i++) {
                        vec[i][k] = 0.0;
                        for (j = 0; j < natoms; j++)
                            vec[i][k] += mmat[i][j] * e_work[j];
                        d1 += vec[i][k] * vec[i][k];
                        d2 += vec[i][k] * e_work[i];
                    }
                    r = fabs((e_vals[k] - d2) / d2);
                    e_vals[k] = d2;
                    assert(d1 >= 0.);
                    d1 = sqrt(d1);
                    for (i = 0; i < natoms; i++) {
                        vec[i][k] /= d1;
                        e_work[i] = vec[i][k];
                    }
                    if (r < 1.e-5)
                        break;
                }
                iters[k] = iter;

                for (i = 0; i < natoms; i++) {
                    for (j = 0; j < natoms; j++) {
                        mmat[j][i] -= e_vals[k] * vec[i][k] * vec[j][k];
                    }
                }
            }

            rval = 0;
            for (k = 0; k < 3; k++) {
                printf("eigenvector %d: %8.3f  iter: %3d\n",
                       k + 1, e_vals[k], iters[k]);
                if (e_vals[k] < 0.0) {
                    rval = 1;
                    break;
                }
                e_vals[k] += 2. * mostneg;
                e_vals[k] = NEARZERO(e_vals[k]);
                assert(e_vals[k] >= 0.);
                f = sqrt(e_vals[k]);
                for (i = 0; i < natoms; i++)
                    nxyz[4 * i + k] = f * vec[i][k];
            }
            printf("\n");

            /* set 4th dimension coord to random value between
               -0.001 and 0.001: */

            for (i = 0; i < natoms; i++)
                nxyz[4 * i + 3] = 0.002 * rand2() - 0.001;

            if (majorize_opt) {
                printf("Guttman transform majorization of coordinates\n");
                printf("  -- this option not implemented in this version!\n");
                exit(1);
/*			majorize( nxyz );    */
            }

            if (!rval)
                break;
        }

        if (rval) {
            printf("Unable to successfully embed\n");
            exit(1);
        }

        /*  free memory   */

        free_vector(e_work, 0, natoms - 1);
        free_vector(e_vals, 0, natoms - 1);
        free_vector(com2, 0, natoms - 1);
        free_matrix(mmat, 0, natoms - 1, 0, natoms - 1);
        free_vector(dmat, 0, natoms * (natoms - 1) / 2);
        free_matrix(vec, 0, natoms - 1, 0, 2);
    }

    return (0);
}

int metrize(BOUNDS_T * bp, int nmetrize)
/*
BOUNDS_T	*bp;
int		nmetrize;
*/
/*   Like the metrize option in mk_dmat, but this updates the bounds
**   matrix (rather than a copy of it), and can be called directly from
**   nab programs to tighten a bounds matrix.
*/
{
    int i, j, k, m, natoms;
    REAL_T lb, ub;
    REAL_T **bmat;

    bmat = bp->b_bmat;
    natoms = bp->b_natoms;

    printf("random metrization for %d pairs\n", nmetrize);

    for (i = 0; i < nmetrize; i++) {
        k = natoms * rand2();
        m = natoms * rand2();
        if (k > m) {
            j = k;
            k = m;
            m = j;              /* (now k <= m ) */
        }
        lb = bmat[m][k];
        ub = bmat[k][m];
        if (ub - lb < 2.0 * EPSILON)
            continue;
        bmat[k][m] = bmat[m][k] = u_rdist(lb, ub);
        trifix(k, m, natoms, bmat);
    }
    return (0);
}

REAL_T db_viol(REAL_T * pos, REAL_T * grad, int *iter)
/*   distance-bounds violations:
**
**    if sqviol is set:
**
**         if d > ub:  e = (d-ub)**2
**         if d < lb:  e = lbpen*(d-lb)**2
**
**    else:
**
**         if d > ub:  e = (d**2/ub**2 -1)
**         if d < lb:  e = lbpen* (2lb**2/(d**2+lb**2) -1)**2
**
**    [where d = distance, ub=upper bound, lb=lower bound]
**    [lbpen = 3.5 default factor comes from dgeom95 ]
**
**    also: add 0.5* k4d * w**2  
*/
{
    int i, j, ca;
    REAL_T etotf;
    REAL_T etot, etotb, etotk, etot4, f, eterm;
    REAL_T dx, dy, dz, dw, d2, gsq, wmax, lb2, ub2, temp, inv, err;
    REAL_T gx, gy, gz, gw, lbpen2, dis;
    CHIRAL_T *cp;
    REAL_T dvol[12], vol, kchired;

    lbpen2 = 2. * lbpen;
    etotb = 0.;
    etotk = 0.;
    etot4 = 0.;
    for (i = 0; i < 4 * natoms; i++)
        grad[i] = 0.0;

    for (i = 0; i < natoms; i++) {
        err = pos[4 * i + 3];
        grad[4 * i + 3] = k4d * err;
        etot4 += 0.5 * k4d * err * err;
    }

    for (i = 0; i < natoms - 1; i++) {
        gx = gy = gz = gw = 0.0;
        for (j = i + 1; j < natoms; j++) {

            /*  skip uninteresting bounds:  */
            if (bmat[j][i] == CLOSEST && bmat[i][j] == FURTHEST)
                continue;

            dx = pos[4 * i + 0] - pos[4 * j + 0];
            dy = pos[4 * i + 1] - pos[4 * j + 1];
            dz = pos[4 * i + 2] - pos[4 * j + 2];
            dw = pos[4 * i + 3] - pos[4 * j + 3];
            d2 = dx * dx + dy * dy + dz * dz + dw * dw;
            ub2 = bmat[i][j] * bmat[i][j];
            lb2 = bmat[j][i] * bmat[j][i];

            if (d2 > ub2) {
                if (sqviol) {
                    dis = sqrt(d2);
                    err = dis - bmat[i][j];
                    eterm = err * err;
                    etotb += eterm;
                    f = 2. * err / dis;
                } else {
                    inv = 1. / ub2;
                    eterm = d2 * inv - 1.;
                    etotb += eterm;
                    f = 2. * inv;
                }
            } else if (d2 < lb2) {
                if (sqviol) {
                    dis = sqrt(d2);
                    err = dis - bmat[j][i];
                    eterm = lbpen * err * err;
                    etotb += eterm;
                    f = lbpen2 * err / dis;
                } else {
                    inv = 1. / (d2 + lb2);
                    temp = lb2 * inv;
                    err = 2. * temp - 1.;
                    f = -8. * lbpen * err * temp * inv;
                    eterm = lbpen * err * err;
                    etotb += eterm;
/*  DAC changes
					inv = 1./d2;
					eterm = lbpen*(lb2*inv - 1.);
					etotb += eterm;
					f = -lbpen2*lb2*inv*inv;
*/
                }
            } else {
                continue;
            }
            if (pencut >= 0.0 && eterm >= pencut) {
                printf("\tdistance : %5d %5d %9.2f %9.2f %9.2f %9.2f\n",
                       i + 1, j + 1, eterm, dis, bmat[j][i], bmat[i][j]);
            }
            dx = f * dx;
            dy = f * dy;
            dz = f * dz;
            dw = f * dw;
            gx += dx;
            gy += dy;
            gz += dz;
            gw += dw;
            grad[4 * j] -= dx;
            grad[4 * j + 1] -= dy;
            grad[4 * j + 2] -= dz;
            grad[4 * j + 3] -= dw;
        }
        grad[4 * i] += gx;
        grad[4 * i + 1] += gy;
        grad[4 * i + 2] += gz;
        grad[4 * i + 3] += gw;
    }

    /*   chirality violations:    */

    for (i = 0; i < nchi; i++) {
        cp = &chi[i];
        chirvol(4, cp->c_anum[0], cp->c_anum[1], cp->c_anum[2],
                cp->c_anum[3], pos, dvol, &vol);
        kchired = kchi / (fabs(cp->c_dist) > 1.0 ? fabs(cp->c_dist) : 1.0);
        eterm = kchired * (vol - cp->c_dist) * (vol - cp->c_dist);
        etotk += eterm;
        if (pencut >= 0.0 && eterm >= pencut) {
            printf("\tchirality: %5d %9.2f %5d %5d %5d %5d %9.2f %9.2f\n",
                   i + 1, eterm, cp->c_anum[0] + 1, cp->c_anum[1] + 1,
                   cp->c_anum[2] + 1, cp->c_anum[3] + 1, vol, cp->c_dist);
        }
        for (j = 0; j < 4; j++) {
            ca = cp->c_anum[j];
            grad[4 * ca + 0] +=
                kchired * dvol[3 * j + 0] * (vol - cp->c_dist) * 2.;
            grad[4 * ca + 1] +=
                kchired * dvol[3 * j + 1] * (vol - cp->c_dist) * 2.;
            grad[4 * ca + 2] +=
                kchired * dvol[3 * j + 2] * (vol - cp->c_dist) * 2.;
        }
    }

    etot = etotb + etotk + etot4;
    if ((ntpr != 0) && (*iter % ntpr == 0 || *iter == 1)) {
        for (gsq = 0.0, i = 0; i < 4 * natoms; i++)
            gsq += grad[i] * grad[i];
        gsq = sqrt(gsq / (4 * natoms));
        for (wmax = 0.0, i = 0; i < natoms; i++)
            wmax =
                fabs(pos[4 * i + 3]) > wmax ? fabs(pos[4 * i + 3]) : wmax;
        if (*iter == 1)
            printf( "               tot       bounds       chi       "
                    "e4d       gsq       wmax\n");
        printf("db:%7d %-10.4g %-10.4g %-10.4g %-10.4g %-10.4g %-10.4g\n",
               *iter, etot, etotb, etotk, etot4, gsq, wmax);
    }
    etotf = etot;
    return (etotf);
}

REAL_T db_viol3(REAL_T * pos, REAL_T * grad, int *iter)
/*   distance-bounds violations:
**
**    if d > ub:  e = (d**2/ub**2 -1)
**
**    if d < lb:  e = lbpen* (2lb**2/(d**2+lb**2) -1)**2
**
**    [where d = distance, ub=upper bound, lb=lower bound]
**    [lbpen = 3.5 default factor comes from dgeom95 ]
**
*/
{
    int i, j, ca;
    REAL_T etotf;
    REAL_T etot, etotb, etotk, f, eterm;
    REAL_T dx, dy, dz, d2, gsq, lb2, ub2, temp, inv, err;
    REAL_T gx, gy, gz, lbpen2, dis;
    CHIRAL_T *cp;
    REAL_T dvol[12], vol, kchired;

    lbpen2 = 2. * lbpen;
    etotb = 0.;
    etotk = 0.;
    for (i = 0; i < 3 * natoms; i++)
        grad[i] = 0.0;

    for (i = 0; i < natoms - 1; i++) {
        gx = gy = gz = 0.0;
        for (j = i + 1; j < natoms; j++) {

            /*  skip uninteresting bounds:  */
            if (bmat[j][i] == CLOSEST && bmat[i][j] == FURTHEST)
                continue;

            dx = pos[3 * i + 0] - pos[3 * j + 0];
            dy = pos[3 * i + 1] - pos[3 * j + 1];
            dz = pos[3 * i + 2] - pos[3 * j + 2];
            d2 = dx * dx + dy * dy + dz * dz;
            ub2 = bmat[i][j] * bmat[i][j];
            lb2 = bmat[j][i] * bmat[j][i];

            if (d2 > ub2) {
                if (sqviol) {
                    dis = sqrt(d2);
                    err = dis - bmat[i][j];
                    eterm = err * err;
                    etotb += eterm;
                    f = 2. * err / dis;
                } else {
                    inv = 1. / ub2;
                    eterm = d2 * inv - 1.;
                    etotb += eterm;
                    f = 2. * inv;
                }
            } else if (d2 < lb2) {
                if (sqviol) {
                    dis = sqrt(d2);
                    err = dis - bmat[j][i];
                    eterm = lbpen * err * err;
                    etotb += eterm;
                    f = lbpen2 * err / dis;
                } else {
                    inv = 1. / (d2 + lb2);
                    temp = lb2 * inv;
                    err = 2. * temp - 1.;
                    f = -8. * lbpen * err * temp * inv;
                    eterm = lbpen * err * err;
                    etotb += eterm;
/*  DAC changes
					inv = 1./d2;
					eterm = lbpen*(lb2*inv - 1.);
					etotb += eterm;
					f = -lbpen2*lb2*inv*inv;
*/
                }
            } else {
                continue;
            }
            if (pencut >= 0.0 && eterm >= pencut) {
                printf("\tdistance : %5d %5d %9.2f\n",
                       i + 1, j + 1, eterm);
            }
            dx = f * dx;
            dy = f * dy;
            dz = f * dz;
            gx += dx;
            gy += dy;
            gz += dz;
            grad[3 * j] -= dx;
            grad[3 * j + 1] -= dy;
            grad[3 * j + 2] -= dz;
        }
        grad[3 * i] += gx;
        grad[3 * i + 1] += gy;
        grad[3 * i + 2] += gz;
    }

    /*   chirality violations:    */

    for (etotk = 0.0, i = 0; i < nchi; i++) {
        cp = &chi[i];
        chirvol(3, cp->c_anum[0], cp->c_anum[1], cp->c_anum[2],
                cp->c_anum[3], pos, dvol, &vol);
        kchired = kchi / (fabs(cp->c_dist) > 1.0 ? fabs(cp->c_dist) : 1.0);
        eterm = kchired * (vol - cp->c_dist) * (vol - cp->c_dist);
        etotk += eterm;
        if (pencut >= 0.0 && eterm >= pencut) {
            printf("\tchirality: %5d       %9.2f\n", i + 1, eterm);
        }
        for (j = 0; j < 3; j++) {
            ca = cp->c_anum[j];
            grad[3 * ca + 0] +=
                kchired * dvol[3 * j + 0] * (vol - cp->c_dist) * 2.;
            grad[3 * ca + 1] +=
                kchired * dvol[3 * j + 1] * (vol - cp->c_dist) * 2.;
            grad[3 * ca + 2] +=
                kchired * dvol[3 * j + 2] * (vol - cp->c_dist) * 2.;
        }
    }

    etot = etotb + etotk;
    if ((ntpr != 0) && (*iter % ntpr == 0 || *iter == 1)) {
        for (gsq = 0.0, i = 0; i < 3 * natoms; i++)
            gsq += grad[i] * grad[i];
        gsq = sqrt(gsq / (3 * natoms));
        if (*iter == 1)
            printf("              tot     bounds      chi      gsq\n");
        printf("db:%7d %-10.4g %-10.4g %-10.4g %-10.4g\n",
               *iter, etot, etotb, etotk, gsq);
    }
    etotf = etot;
    return (etotf);
}

static void mk_dmat(BOUNDS_T * bp, REAL_T * dmat, REAL_T(*rdist) ())
{
    int i, j, k, m, npair, nmetrize, ij;
    REAL_T lb, ub;
    REAL_T **bmat, **bmatc;
    REAL_T gap;

    bmat = bp->b_bmat;

    npair = natoms * (natoms - 1) / 2;
    nmetrize = 0.01 * randpair * npair;
    if (randpair > 0.0) {
        /* use partial randomized pairwise distance metrization */

        if ((bmatc = matrix(0, natoms - 1, 0, natoms - 1)) == NULL) {
            rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "bmatc");
            exit(1);
        }

        for (i = 0; i < natoms; i++) {
            for (j = 0; j < natoms; j++) {
                bmatc[i][j] = bmat[i][j];
            }
        }
        printf("random metrization for %d pairs\n", nmetrize);

        for (i = 0; i < nmetrize; i++) {
            k = natoms * rand2();
            m = natoms * rand2();
            if (k > m) {
                j = k;
                k = m;
                m = j;          /* (now k <= m ) */
            }
            lb = bmatc[m][k];
            ub = bmatc[k][m];
            if (ub - lb < 2.0 * EPSILON)
                continue;
            bmatc[k][m] = bmatc[m][k] = (*rdist) (lb, ub);
            trifix(k, m, natoms, bmatc);
            if (i % 5000 == 0 && i != 0)
                fprintf(stderr, "completed %5d calls to trifix\n", i);
        }

        gap = 0.0;
        ij = 0;
        for (i = 0; i < natoms - 1; i++) {
            for (j = i + 1; j < natoms; j++) {
                lb = bmatc[j][i];
                ub = bmatc[i][j];
                gap += ub - lb;
                if (ub - lb < 2.0 * EPSILON) {
                    dmat[ij] = 0.5 * (ub + lb);
                    ij++;
                } else {
                    dmat[ij] = (*rdist) (lb, ub);
                    ij++;
                }
            }
        }
        free_matrix(bmatc, 0, natoms - 1, 0, natoms - 1);

    } else {
        /* just select distances, don't metrize:  */

        gap = 0.0;
        ij = 0;
        for (i = 0; i < natoms - 1; i++) {
            for (j = i + 1; j < natoms; j++) {
                lb = bmat[j][i];
                ub = bmat[i][j];
                gap += ub - lb;
                if (ub - lb < 2.0 * EPSILON) {
                    dmat[ij++] = 0.5 * (ub + lb);
                } else {
                    dmat[ij++] = (*rdist) (lb, ub);
                }
            }
        }
    }
    gap /= (npair - nmetrize);
    printf("average bounds gap = %8.3f\n", gap);

}

static int rd_mmat(char op[], BOUNDS_T * bp, REAL_T ** mmat)
{
    char *mmfname;
    FILE *mmfp;
    int i, j;
    char line[256];
    int ri, rj;
    char ainame[20], riname[20];
    char ajname[20], rjname[20];
    char name[20], l_name[20];
    REAL_T d;

    if (!(mmfname = strchr(op, '='))) {
        fprintf(stderr, "rd_mmat: syntax: %s\n", op);
        exit(1);
    }
    mmfname++;
    if (!*mmfname) {
        fprintf(stderr, "rd_mmat: syntax: %s\n", op);
        exit(1);
    }
    if ((mmfp = fopen(mmfname, "r")) == NULL) {
        fprintf(stderr, "rd_mmat: can't open metric mat file %s\n",
                mmfname);
        exit(1);
    }
    for (i = 0; i < bp->b_natoms; i++)
        mmat[i][i] = 0.0;
    *l_name = '\0';
    for (i = 0; fgets(line, sizeof(line), mmfp);) {
#ifdef NAB_DOUBLE_PRECISION
        sscanf(line, "%d %s %s %d %s %s %lf",
#else
        sscanf(line, "%d %s %s %d %s %s %f",
#endif
               &ri, riname, ainame, &rj, rjname, ajname, &d);
        sprintf(name, "%d %s %s", ri, riname, ainame);
        if (strcmp(l_name, name)) {
            i++;
            j = 1;
        }
        mmat[i][j] = d;
        mmat[j][i] = d;
        strcpy(l_name, name);
        j++;
    }
    fclose(mmfp);
    return (0);
}

static int mk_mmat(BOUNDS_T * bp, REAL_T * dmat, REAL_T ** mmat)
{
    int i, j, n, ij;
    int nneg;
    REAL_T d2, rog2, mij;
    FILE *mfp;

    n = natoms;

    mostneg = 0.;
    for (i = 0; i < n; i++)
        com2[i] = 0.0;

    rog2 = 0.0;
    ij = 0;
    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            d2 = dmat[ij] * dmat[ij];
            ij++;
            rog2 += d2;
            com2[i] += d2;
            com2[j] += d2;
        }
    }
    rog2 /= n;

    for (nneg = 0, i = 0; i < n; i++) {
        com2[i] -= rog2;
        com2[i] /= n;
        if (com2[i] < 0.0) {
            nneg++;
            mostneg = com2[i] < mostneg ? com2[i] : mostneg;
        }
    }

    if (nneg > 0) {
        fprintf(stderr, "mk_mmat: %d neg COM dist: %8.3f\n", nneg,
                mostneg);
    }

    for (i = 0; i < n; i++)
        mmat[i][i] = com2[i] - 2. * mostneg;
    ij = 0;
    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            mij = 0.5 * (com2[i] + com2[j] - dmat[ij] * dmat[ij]);
            ij++;
            if (mij * mij / (com2[i] * com2[j]) > 1.) {
                assert(com2[i] * com2[j] >= 0.);
                mij = (mij > 0. ? 0.95 : -0.95) * sqrt(com2[i] * com2[j]);
            }
            mmat[j][i] = mmat[i][j] = mij;
        }
    }

    if (dmm) {
        if ((mfp = fopen(dmmn, "w"))) {
            dumpmat(mfp, bp, mmat);
            fclose(mfp);
        } else {
            fprintf(stderr,
                    "mk_mmat: can't write metric matrix file %s\n", dmmn);
            exit(1);
        }
    }

    return (0);
}

static void trifix(int p, int q, int n, REAL_T ** b)
/*
int		p,q,n;
REAL_T	**b;
*/
/*
**   update a bounds matrix b (of dimension n x n) to reflect changes
**   that are implied by a recent change in the bounds connecting atoms
**   p and q.
**
**   Algorithm here was developed by Jay Ponder; see J. Mol. Biol.
**     264: 585-602 (1996).
*/
{

    int i, k, ip, iq, np, nq;
    int *pt, *qt, *pc, *qc;
    REAL_T eps, d_ip_min, d_ip_max, d_iq_min, d_iq_max;
    REAL_T *pmin, *pmax, *qmin, *qmax;


    if ((pt = ivector(0, n - 1)) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "pt");
        exit(1);
    }

    if ((qt = ivector(0, n - 1)) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "qt");
        exit(1);
    }

    if ((pc = ivector(0, n - 1)) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "pc");
        exit(1);
    }

    if ((qc = ivector(0, n - 1)) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "qc");
        exit(1);
    }

    if ((pmin = vector(0, n - 1)) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "pmin");
        exit(1);
    }

    if ((pmax = vector(0, n - 1)) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "pmax");
        exit(1);
    }

    if ((qmin = vector(0, n - 1)) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "qmin");
        exit(1);
    }

    if ((qmax = vector(0, n - 1)) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "qmax");
        exit(1);
    }

    eps = 1.0e-4;
    np = nq = 0;
    for (i = 0; i < n; i++)
        pc[i] = qc[i] = 1;

    for (i = 0; i < p; i++) {
        pmin[i] = b[p][i];
        pmax[i] = b[i][p];
    }
    for (i = p; i < n; i++) {
        pmin[i] = b[i][p];
        pmax[i] = b[p][i];
    }
    for (i = 0; i < q; i++) {
        qmin[i] = b[q][i];
        qmax[i] = b[i][q];
    }
    for (i = q; i < n; i++) {
        qmin[i] = b[i][q];
        qmax[i] = b[q][i];
    }

/* check upper bounds:  */

    for (i = 0; i < n; i++) {
        d_ip_max = qmax[p] + qmax[i];
        if (pmax[i] > d_ip_max + eps) {
            pt[np++] = i;
            pmax[i] = d_ip_max;
            pc[i] = 0;
        }
        d_iq_max = pmax[q] + pmax[i];
        if (qmax[i] > d_iq_max + eps) {
            qt[nq++] = i;
            qmax[i] = d_iq_max;
            qc[i] = 0;
        }
    }

    for (ip = 0; ip < np; ip++) {
        i = pt[ip];
        d_ip_max = pmax[i];
        for (iq = 0; iq < nq; iq++) {
            k = qt[iq];
            if (i < k) {
                if (b[i][k] > d_ip_max + pmax[k])
                    b[i][k] = d_ip_max + pmax[k];
            } else {
                if (b[k][i] > d_ip_max + pmax[k])
                    b[k][i] = d_ip_max + pmax[k];
            }
        }
    }

/* check lower bounds:  */

    for (i = 0; i < n; i++) {
        d_ip_min = qmin[p] - qmax[i] > qmin[i] - qmax[p] ?
            qmin[p] - qmax[i] : qmin[i] - qmax[p];
        if (pmin[i] < d_ip_min - eps) {
            if (pc[i])
                pt[np++] = i;
            pmin[i] = d_ip_min;
        }

        d_iq_min = pmin[q] - pmax[i] > pmin[i] - pmax[q] ?
            pmin[q] - pmax[i] : pmin[i] - pmax[q];
        if (qmin[i] < d_iq_min - eps) {
            if (qc[i])
                qt[nq++] = i;
            qmin[i] = d_iq_min;
        }
    }

    for (ip = 0; ip < np; ip++) {
        i = pt[ip];
        d_ip_min = pmin[i];
        d_ip_max = pmax[i];
        for (iq = 0; iq < nq; iq++) {
            k = qt[iq];
            if (i < k) {
                if (d_ip_min - pmax[k] > b[k][i])
                    b[k][i] = d_ip_min - pmax[k];
                if (pmin[k] - d_ip_max > b[k][i])
                    b[k][i] = pmin[k] - d_ip_max;
            } else {
                if (d_ip_min - pmax[k] > b[i][k])
                    b[i][k] = d_ip_min - pmax[k];
                if (pmin[k] - d_ip_max > b[i][k])
                    b[i][k] = pmin[k] - d_ip_max;
            }
        }
    }

/*  some debug output: 

for( i=0; i<n; i++ ){
	if( pmin[i] - eps > pmax[i] )
		fprintf( stderr, "triout: %3d %3d %8.3f %8.3f %8.3f %8.3f\n", 
			i,p,pmin[i],pmax[i], b[i][p], b[p][i] );
	if( qmin[i] - eps > qmax[i] )
		fprintf( stderr, "triout: %3d %3d %8.3f %8.3f %8.3f %8.3f\n", 
			i,q,qmin[i],qmax[i], b[i][q], b[q][i] );
}

*/

    for (i = 0; i < p; i++) {
        b[p][i] = pmin[i];
        b[i][p] = pmax[i];
    }
    for (i = p; i < n; i++) {
        b[i][p] = pmin[i];
        b[p][i] = pmax[i];
    }
    for (i = 0; i < q; i++) {
        b[q][i] = qmin[i];
        b[i][q] = qmax[i];
    }
    for (i = q; i < n; i++) {
        b[i][q] = qmin[i];
        b[q][i] = qmax[i];
    }

    free_ivector(pt, 0, n - 1);
    free_ivector(qt, 0, n - 1);
    free_ivector(pc, 0, n - 1);
    free_ivector(qc, 0, n - 1);
    free_vector(pmin, 0, n - 1);
    free_vector(pmax, 0, n - 1);
    free_vector(qmin, 0, n - 1);
    free_vector(qmax, 0, n - 1);

    return;

}

static REAL_T u_rdist(REAL_T lb, REAL_T ub)
{
    REAL_T d;

    d = (ub - lb) * rand2() + lb;
    return (d);
}

static REAL_T g_rdist(REAL_T lb, REAL_T ub)
{
    REAL_T d, mu, sigma;

    mu = lb + 0.65 * (ub - lb);
    sigma = (ub - lb) / 4.0;
    d = gauss2(&mu, &sigma);
    if (d < lb)
        d = lb;
    if (d > ub)
        d = ub;
    return (d);
}

static void dumpmat(FILE * fp, BOUNDS_T * bp, REAL_T ** mat)
/*
FILE	*fp;
BOUNDS_T	*bp;
REAL_T	**mat;
*/
{
    int a1, a2;
    char a1nm[20], a2nm[20];

    for (a1 = 0; a1 < bp->b_natoms - 1; a1++) {
        mkdgname(bp, a1, a1nm);
        for (a2 = a1 + 1; a2 < bp->b_natoms; a2++) {
            mkdgname(bp, a2, a2nm);
            fprintf(fp, "%-12s %-12s %9.3f\n", a1nm, a2nm, mat[a1][a2]);
        }
    }
}

static int makeBoundList(BOUNDS_T * bp)
/*
  make a list of bounds from each atom to every higher atom number in which
  there exists a lower bound != van der waals or an upper bound < infinity.
  The bound list is stored in a natoms*natoms matrix boundList, which is
  attached to the bounds struct.
*/
{
    int i = 0, j = 0, boundCounter;
    int **boundList, natoms;
    int big, small;
    REAL_T vdw, delta = 1e-9;

    natoms = bp->b_natoms;
    if ((boundList = (int **) imatrix(0, natoms, 0, natoms)) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "boundList");
        exit(1);
    }
    bp->b_boundList = boundList;

    for (i = 0; i < natoms; i++) {
        boundCounter = 0;
        for (j = 0; j < natoms; j++) {
            bp->b_boundList[i][j] = -1;
            if (j == i)
                continue;
            vdw = vdw_radius(bp, i) + vdw_radius(bp, j);
            if (CLOSEST > vdw)
                vdw = CLOSEST;
            if (j > i) {
                big = j;
                small = i;
            } else {
                big = i;
                small = j;
            }
            if (((bp->b_bmat[big][small] < vdw - delta) ||
                 (bp->b_bmat[big][small] > vdw + delta)) &&
                (bp->b_bmat[big][small] > delta)) {
                bp->b_boundList[i][boundCounter++] = j;
            } else if (bp->b_bmat[small][big] < FURTHEST - delta) { /* upper bound */
                bp->b_boundList[i][boundCounter++] = j;
            }
        }
    }

    return (0);
}

static int findPaths(BOUNDS_T * bp, int root)

/*
  Routine to find triangle smoothed upper and lower bounds of each
  atom to a specified root atom using a sparse variant of the
  Bellman-Ford shortest path algorithm. routine mimicks pseudocode
  provided in Crippen & Havel's "Distance Geometry and Molecular
  Conformation", Algorithm 6.9, p. 301, with several differences: 
  (1) where the reference algorithm uses one while loop for the atom
  queue, we use two - one for the lower bound path lengths and one
  for the upper bound path lengths.  (2) where the reference algorithm
  checks to see if PathLength[ j ] - LowerBound[ j, -k ] > PathLength[ k ], 
  we use the absolute value to reflect the possibility that the
  PathLength[j] may be less than the LowerBound[j,-k].
*/
{
    int *queue, *inQueue, i, j, k, l, head, tail, natoms, enter;
    REAL_T *lowerBoundPathLength, *upperBoundPathLength, upperBound,
        lowerBound;
    REAL_T big, small, smaller, vdw, delta = 1e-9;

    natoms = bp->b_natoms;
    if ((queue = (int *) malloc(natoms * sizeof(int))) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "queue");
        exit(1);
    }
    if ((inQueue = (int *) malloc(natoms * sizeof(int))) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "inQueue");
        exit(1);
    }
    if ((lowerBoundPathLength =
         (REAL_T *) malloc(natoms * sizeof(REAL_T))) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "lowerBoundPathLength");
        exit(1);
    }
    if ((upperBoundPathLength =
         (REAL_T *) malloc(natoms * sizeof(REAL_T))) == NULL) {
        rt_errormsg_s(TRUE, E_NOMEM_FOR_S, "upperBoundPathLength");
        exit(1);
    }
    for (i = 0; i < natoms; i++) {  /* initialize loop */
        queue[i] = -1;
        inQueue[i] = FALSE;
        upperBoundPathLength[i] = FURTHEST;
        lowerBoundPathLength[i] = 0.0;
    }
    head = root;
    tail = root;
    upperBoundPathLength[root] = 0.0;
    while (head != -1) {
        j = head;
        inQueue[j] = FALSE;
        head = queue[head];
        for (l = 0; ((bp->b_boundList[j][l] != -1) && (l < natoms - 1));
             l++) {
            k = bp->b_boundList[j][l];
            if (k == root)
                continue;
            else if (k < j) {
                lowerBound = bp->b_bmat[j][k];
                upperBound = bp->b_bmat[k][j];
            } else {
                lowerBound = bp->b_bmat[k][j];
                upperBound = bp->b_bmat[j][k];
            }
            big = upperBoundPathLength[j] + upperBound;
            small = lowerBound - upperBoundPathLength[j];
            smaller = lowerBoundPathLength[j] - upperBound;
            if (small < smaller)
                small = smaller;
            enter = FALSE;
            if (upperBoundPathLength[k] > big + delta) {
                upperBoundPathLength[k] = big;
                if (!inQueue[k])
                    enter = TRUE;
            }
            if (lowerBoundPathLength[k] < small - delta) {
                if (root > k)
                    vdw = bp->b_bmat[root][k];
                else
                    vdw = bp->b_bmat[k][root];
                if (small < vdw)
                    small = vdw;
                lowerBoundPathLength[k] = small;
                if (!inQueue[k])
                    enter = TRUE;
            }
            if (enter) {
                inQueue[k] = TRUE;
                if (head == -1)
                    head = k;
                else
                    queue[tail] = k;
                queue[k] = -1;
                tail = k;
            }
        }
    }
    bp->b_lowerBoundPathLength = lowerBoundPathLength;
    bp->b_upperBoundPathLength = upperBoundPathLength;
    return (0);
}

int geodesics(BOUNDS_T * bp)
/*
  smooths the bounds matrix via the triangle inequality using a sparse
  matrix version of a shortest path algorithm.  See Crippen & Havel,
  "Distance Geometry and Molecular Conformation, Algorithm 6.9, p. 301.
*/
{
    int i, m;
    REAL_T lb, ub, lowerBoundPathLength_to_m, upperBoundPathLength_to_m,
        natoms;
    char inm[40], mnm[40];

    makeBoundList(bp);
    natoms = bp->b_natoms;
    for (i = 0; i < natoms; i++) {
        findPaths(bp, i);
        for (m = 0; m < i; m++) {
            lowerBoundPathLength_to_m = bp->b_lowerBoundPathLength[m];
            upperBoundPathLength_to_m = bp->b_upperBoundPathLength[m];
            lb = bp->b_bmat[i][m];
            ub = bp->b_bmat[m][i];
            if (lowerBoundPathLength_to_m > lb) {
                bp->b_bmat[i][m] = lowerBoundPathLength_to_m;
            }
            if (upperBoundPathLength_to_m < ub) {
                bp->b_bmat[m][i] = upperBoundPathLength_to_m;
            }
            if (ub < lb) {
                mkdgname(bp, i, inm);
                mkdgname(bp, m, mnm);
                fprintf(stderr,
                        "BOUNDS VIOLATION ( inverting lb, ub ): %s - %s ( %f - %f )\n",
                        inm, mnm, lb, ub);
                lb = bp->b_bmat[m][i];  /* bounds violations: */
                ub = bp->b_bmat[i][m];  /* ub must have shortest path */
                bp->b_bmat[i][m] = lb;  /* switch ub and lb */
                bp->b_bmat[m][i] = ub;
            }
        }
        for (m = i + 1; m < natoms; m++) {
            lowerBoundPathLength_to_m = bp->b_lowerBoundPathLength[m];
            upperBoundPathLength_to_m = bp->b_upperBoundPathLength[m];
            lb = bp->b_bmat[m][i];
            ub = bp->b_bmat[i][m];
            if (lowerBoundPathLength_to_m > lb) {
                bp->b_bmat[m][i] = lowerBoundPathLength_to_m;
            }
            if (upperBoundPathLength_to_m < ub) {
                bp->b_bmat[i][m] = upperBoundPathLength_to_m;
            }
            if (ub < lb) {
                mkdgname(bp, i, inm);
                mkdgname(bp, m, mnm);
                fprintf(stderr,
                        "BOUNDS VIOLATION ( inverting lb, ub ): %s - %s ( %f - %f )\n",
                        inm, mnm, lb, ub);
                lb = bp->b_bmat[i][m];  /* bounds violations: */
                ub = bp->b_bmat[m][i];  /* ub must have shortest path */
                bp->b_bmat[m][i] = lb;  /* switch ub and lb */
                bp->b_bmat[i][m] = ub;
            }
        }
        free(bp->b_lowerBoundPathLength);
        free(bp->b_upperBoundPathLength);
    }
    return (0);
}

int tsmooth(BOUNDS_T * bp, REAL_T delta)
{
    int i, j, k, ikmin, ikmax, jstart, nviol;
    REAL_T temp, q, l_ik, u_ik, l_ij, u_ij, l_jk, u_jk;
    REAL_T **bmat;

    bmat = bp->b_bmat;
    q = 0.;
    nviol = 0;

    for (k = 0; k < bp->b_natoms; k++) {

        if (k % 500 == 0 && k != 0)
            fprintf(stderr, "ready to triangle smooth atom %5d\n", k);

        for (i = 0; i < bp->b_natoms - 1; i++) {

            if (i == k)
                continue;
            ikmin = i < k ? i : k;
            ikmax = i > k ? i : k;
            l_ik = bmat[ikmax][ikmin];
            u_ik = bmat[ikmin][ikmax];

            /*  if a lower bound is greater than an upper bound, 
               we must already have reported this; don't continue
               to tsmooth from what are unreliable bounds
             */
            if (l_ik - u_ik > delta)
                continue;

            for (j = i + 1; j < k; j++) {

                l_ij = bmat[j][i];
                u_ij = bmat[i][j];
                l_jk = bmat[k][j];
                u_jk = bmat[j][k];

                /*  if a lower bound is greater than an upper bound, 
                   we must already have reported this; don't continue
                   to tsmooth from what are unreliable bounds
                 */
                if (l_ij - u_ij > delta || l_jk - u_jk > delta)
                    continue;

                if (u_ij > u_ik + u_jk) {
                    bmat[i][j] = u_ik + u_jk;
#					ifdef TS_DEBUG
                    fprintf(stderr, "udate: %4d %4d %8.3f\n", i + 1, j + 1,
                            u_ik + u_jk);
                    fprintf(stderr,
                            "            %4d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                            k + 1, l_ij, u_ij, l_ik, u_ik, l_jk, u_jk);
#					endif
                }
                temp =
                    l_ik - u_jk > l_jk - u_ik ? l_ik - u_jk : l_jk - u_ik;
                if (l_ij < temp) {
                    bmat[j][i] = temp;
#					ifdef TS_DEBUG
                    fprintf(stderr, "ldate: %4d %4d %8.3f\n", i + 1, j + 1,
                            temp);
                    fprintf(stderr,
                            "            %4d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                            k + 1, l_ij, u_ij, l_ik, u_ik, l_jk, u_jk);
#					endif
                }

                q = bmat[j][i] - bmat[i][j];
                if (q > delta) {
                    nviol++;
                    printf( "T1err: %4d %4d %4d %7.2f %7.2f %7.2f %7.2f "
                            "%7.2f %7.2f %7.2f %7.2f\n",
                         i + 1, j + 1, k + 1, bmat[j][i], bmat[i][j], l_ij,
                         u_ij, l_ik, u_ik, l_jk, u_jk);
                }

            }
            jstart = i > k ? i + 1 : k + 1;
            for (j = jstart; j < bp->b_natoms; j++) {

                l_ij = bmat[j][i];
                u_ij = bmat[i][j];
                l_jk = bmat[j][k];
                u_jk = bmat[k][j];

                /*  if a lower bound is greater than an upper bound, 
                   we must already have reported this; don't continue
                   to tsmooth from what are unreliable bounds
                 */
                if (l_ij - u_ij > delta || l_jk - u_jk > delta)
                    continue;

                if (u_ij > u_ik + u_jk) {
                    bmat[i][j] = u_ik + u_jk;
#					ifdef TS_DEBUG
                    fprintf(stderr, "udate: %4d %4d %8.3f\n", i + 1, j + 1,
                            u_ik + u_jk);
                    fprintf(stderr,
                            "            %4d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                            k + 1, l_ij, u_ij, l_ik, u_ik, l_jk, u_jk);
#					endif
                }
                temp =
                    l_ik - u_jk > l_jk - u_ik ? l_ik - u_jk : l_jk - u_ik;
                if (l_ij < temp) {
                    bmat[j][i] = temp;
#					ifdef TS_DEBUG
                    fprintf(stderr, "ldate: %4d %4d %8.3f\n", i + 1, j + 1,
                            temp);
                    fprintf(stderr,
                            "            %4d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                            k + 1, l_ij, u_ij, l_ik, u_ik, l_jk, u_jk);
#					endif
                }

                q = bmat[j][i] - bmat[i][j];
                if (q > delta) {
                    nviol++;
                    printf( "T2err: %4d %4d %4d (%7.2f %7.2f)(%7.2f %7.2f)"
                            "(%7.2f %7.2f)(%7.2f %7.2f)\n",
                         i + 1, j + 1, k + 1, bmat[j][i], bmat[i][j], l_ij,
                         u_ij, l_ik, u_ik, l_jk, u_jk);
                }
            }
        }

        if (nviol) {            /* need to sweep through bounds: swap ub
                                   and lb if bounds are out of order   */

            for (i = 0; i < bp->b_natoms - 1; i++) {
                for (j = i + 1; j < bp->b_natoms; j++) {
                    if (bmat[j][i] > bmat[i][j]) {
                        temp = bmat[j][i];
                        bmat[j][i] = bmat[i][j];
                        bmat[i][j] = temp;
                    }
                }
            }
        }

    }
    return (0);
}

