/* 
   Based on a Fortran routine by S. Sridharan.  See S. Sridharan, A.
   Nicholls and K.A. Sharp, "A rapid method for calculating derivatives of
   solvent accessible surface areas of molecules", J. Computat. Chem.
   16: 1038-1044 (1995).

   Adapted for NAB by DAC, 3/98 

   Input: atomic coordinates, radii (vdw+water probe), the number of atoms,
   and a surface free energy (surface tension value in kcal/A**2)

   Output: forces in kcal/mole/A  into array "fat"
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sff.h"

#define MXG 31
#define MXACUM 250000
#define MXPR 250000
#define NVER 260
#define NEDGE 520

int	sasad( REAL_T *, REAL_T *,  REAL_T *, int, REAL_T );
static	int	cube( REAL_T, int, REAL_T *, REAL_T *, REAL_T *, REAL_T *, REAL_T *, REAL_T *,
	int *, int *, int * );

int	sasad( REAL_T *crd, REAL_T *r0,  REAL_T *fat, int natm, REAL_T gamma )
/*
REAL_T *crd, *r0, *fat;
int natm;
REAL_T gamma;
*/
{

    REAL_T tiny = (REAL_T)1e-6;

    int i__2, i__3;
    REAL_T r__1, r__2;

    int *cbal;
    int edgv[2*NEDGE]	/* was [2][NEDGE] */;
    int ilvl, nlvl;
    REAL_T rvmg;
    int nprx, i, j, k;
    REAL_T hfact, d2, x1, x2, x3;
    int ie, ii, ne;
    REAL_T *r02;
    int jj, ip;
    REAL_T rm[9];
    int iv;
    REAL_T tm;
    int nv;
    REAL_T xm, ym;
    int st[NEDGE];
    REAL_T cf1, cf2, cf3;
    int ic1, ie1, ie2, ic2, ic3;
    REAL_T r0a;
    int ia1, ia2;
    REAL_T ds2, dx1, dx2, dx3;
    int iv1, iv2, ix1, ix2, ix3;
    REAL_T sm1, sm2, *rs2, xo1, xo2, xo3, rv1, rv2, dy1, dy2, dy3;
    int edg[NEDGE];
    REAL_T rad, dmg, frc, ctf, rdn, rij, rsa, pre, csp, tta, vmg, 
	    cst;
    int *oti;
    REAL_T *ver;
    int *pls, nvi, npr, nvo;
    REAL_T snt;
    int nst, nxv, *cbn1, *cbn2;
    REAL_T frc1, frc2, frc3, ctf2, tij1, tij2, tij3, cbai;

    /* Parameter adjustments */
    fat -= 4;
    --r0;
    crd -= 4;

	/* allocate space */
	if( ( cbn1 = ( int * )malloc( MXG*MXG*MXG*sizeof(int) ) ) == NULL ){
		fprintf( stderr, "Unable to allocate space for cbn1: %d\n", MXG );
		exit( 1 );
	}
	if( ( cbn2 = ( int * )malloc( MXG*MXG*MXG*sizeof(int) ) ) == NULL ){
		fprintf( stderr, "Unable to allocate space for cbn2: %d\n", MXG );
		exit( 1 );
	}
	if( ( pls = ( int * )malloc( 2*MXPR*sizeof(int) ) ) == NULL ){
		fprintf( stderr, "Unable to allocate space for pls: %d\n", MXPR );
		exit( 1 );
	}
	if( ( ver = ( REAL_T * )malloc( 3*NVER*sizeof(REAL_T) ) ) == NULL ){
		fprintf( stderr, "Unable to allocate space for ver: %d\n", NVER );
		exit( 1 );
	}
	if( ( oti = ( int * )malloc( NVER*sizeof(int) ) ) == NULL ){
		fprintf( stderr, "Unable to allocate space for oti: %d\n", NVER );
		exit( 1 );
	}
	if( ( cbal = ( int * )malloc( MXACUM*sizeof(int) ) ) == NULL ){
		fprintf( stderr, "Unable to allocate space for cbal: %d\n", MXACUM );
		exit( 1 );
	}
	if( ( r02 = ( REAL_T * )malloc( natm*sizeof(REAL_T) ) ) == NULL ){
		fprintf( stderr, "Unable to allocate space for r02: %d\n", natm );
		exit( 1 );
	}
	if( ( rs2 = ( REAL_T * )malloc( natm*sizeof(REAL_T) ) ) == NULL ){
		fprintf( stderr, "Unable to allocate space for rs2: %d\n", natm );
		exit( 1 );
	}

    /* Function Body */

/* nvi=initial no. of points ; nlvl = no. of levels in the hierarchy */
/* final no. of points on unit circle=nvi*2**nlvl */
/*thus nvi=8,10,12 and nlvl=2 will lead to 32,40 and 48 points on the circle*/
/* the above point densities seem to be the optimum numbers for nlvl=2; if */
/*higher point densities are desired it is better to increase nlvl than nvi*/

    nlvl = 4;
    nvi = 8;

    cube(2.0, natm, &crd[4], &r0[1], &cbai, &xo1, &xo2, &xo3, cbn1, cbn2, cbal);

    tta = (REAL_T)360. / (REAL_T) nvi;
    for (i = 1; i <= nvi; ++i) {
		rdn = (i - 1) * tta;
		ver[i * 3 - 3] = cos(rdn);
		ver[i * 3 - 2] = sin(rdn);
		ver[i * 3 - 1] = 0.;
		j = i + 1;
		if (i == nvi) {
	    	j = 1;
		}
		edgv[(i << 1) - 2] = i;
		edgv[(i << 1) - 1] = j;
    }
    nv = nvi;
    ne = nvi;
    ie2 = 0;

    for (ilvl = 1; ilvl <= nlvl; ++ilvl) {
		ie1 = ie2 + 1;
		ie2 = ne;
		for (ie = ie1; ie <= ie2; ++ie) {
	    	iv1 = edgv[(ie << 1) - 2];
	    	iv2 = edgv[(ie << 1) - 1];
	    	xm = ver[iv1 * 3 - 3] + ver[iv2 * 3 - 3];
	    	ym = ver[iv1 * 3 - 2] + ver[iv2 * 3 - 2];
	    	vmg = sqrt(xm * xm + ym * ym);
	    	++nv;
	    	ver[nv * 3 - 3] = xm / vmg;
	    	ver[nv * 3 - 2] = ym / vmg;
	    	ver[nv * 3 - 1] = (REAL_T)0.;
	    	++ne;
	    	edg[ie - 1] = ne;
	    	edgv[(ne << 1) - 2] = iv1;
	    	edgv[(ne << 1) - 1] = nv;
	    	++ne;
	    	edgv[(ne << 1) - 2] = nv;
	    	edgv[(ne << 1) - 1] = iv2;
		}
    }
    ne = ie2;
    for (ie = ie1; ie <= ie2; ++ie) edg[ie - 1] = -1;

    if (nv > NVER || ne > 520) {
		fprintf( stderr, "nv or ne exceeds max value: %d %d\n", nv, ne );
		exit(1);
    }
    hfact = gamma * 3.14159265 / nv;

    for (i = 1; i <= natm; ++i) {
		r0a = r0[i];
		rsa = r0[i] * .99999;
		r02[i - 1] = r0a * r0a;
		rs2[i - 1] = rsa * rsa;
    }
    npr = 0;
    for (i = 1; i <= natm; ++i) {
		rad = r0[i];
		x1 = crd[i * 3 + 1];
		x2 = crd[i * 3 + 2];
		x3 = crd[i * 3 + 3];
		ix1 = (x1 - xo1) * cbai;
		ix2 = (x2 - xo2) * cbai;
		ix3 = (x3 - xo3) * cbai;
		i__2 = cbn2[ix1 + (ix2 + ix3 * MXG) * MXG];
		for (jj = cbn1[ix1 + (ix2 + ix3 * MXG) * MXG]; jj <= i__2; ++jj) {
	    	j = cbal[jj - 1];
	    	if (j > i) {
				ctf = rad + r0[j];
				ctf2 = ctf * ctf;
				dx1 = crd[j * 3 + 1] - x1;
				dx2 = crd[j * 3 + 2] - x2;
				dx3 = crd[j * 3 + 3] - x3;
				d2 = dx1 * dx1 + dx2 * dx2 + dx3 * dx3;
				if (d2 < ctf2) {
		    		++npr;
		    		pls[(npr << 1) - 2] = i;
		    		pls[(npr << 1) - 1] = j;
				}
	    	}
		}
    }

    if (npr > MXPR) {
		fprintf( stderr, " # pairs exceeds max value: %d\n", npr );
		exit(1);
    }

    cube(1.0, natm, &crd[4], &r0[1], &cbai, &xo1, &xo2, &xo3, cbn1, cbn2, cbal);

    nprx = 0;
    nxv = 0;

    for (ip = 1; ip <= npr; ++ip) {
		i = pls[(ip << 1) - 2];
		j = pls[(ip << 1) - 1];
		dx1 = crd[j * 3 + 1] - crd[i * 3 + 1];
		dx2 = crd[j * 3 + 2] - crd[i * 3 + 2];
		dx3 = crd[j * 3 + 3] - crd[i * 3 + 3];
		d2 = dx1 * dx1 + dx2 * dx2 + dx3 * dx3;
		dmg = sqrt(d2);
		pre = (r02[i - 1] - r02[j - 1]) / d2 + (REAL_T)1.;
		tij1 = crd[i * 3 + 1] + pre * (REAL_T).5 * dx1;
		tij2 = crd[i * 3 + 2] + pre * (REAL_T).5 * dx2;
		tij3 = crd[i * 3 + 3] + pre * (REAL_T).5 * dx3;

		r__1 = r0[i] + r0[j];
		r__2 = r0[i] - r0[j];
		rij = sqrt(r__1 * r__1 - d2) * (REAL_T).5 * sqrt(d2 - r__2 * r__2) / dmg;
		rvmg = sqrt(dx1 * dx1 + dx2 * dx2);
		rv1 = -dx2 / (rvmg + tiny);
		rv2 = dx1 / (rvmg + tiny);
		cst = dx3 / dmg;
		snt = sqrt((REAL_T)1. - cst * cst);
		csp = (REAL_T)1. - cst;
		tm = csp * rv1;
		sm1 = snt * rv1;
		sm2 = snt * rv2;
		rm[0] = tm * rv1 + cst;
		rm[3] = tm * rv2;
		rm[6] = sm2;
		rm[1] = tm * rv2;
		rm[4] = csp * rv2 * rv2 + cst;
		rm[7] = -sm1;
		rm[2] = -sm2;
		rm[5] = sm1;
		rm[8] = cst;
		nvo = 0;
		for (iv = 1; iv <= nvi; ++iv) {
	    	cf1 = rm[0] * ver[iv * 3 - 3] + rm[3] * ver[iv * 3 - 2] + rm[6] * 
		    	ver[iv * 3 - 1];
	    	cf2 = rm[1] * ver[iv * 3 - 3] + rm[4] * ver[iv * 3 - 2] + rm[7] * 
		    	ver[iv * 3 - 1];
	    	cf3 = rm[2] * ver[iv * 3 - 3] + rm[5] * ver[iv * 3 - 2] + rm[8] * 
		    	ver[iv * 3 - 1];
	    	cf1 = tij1 + rij * cf1;
	    	cf2 = tij2 + rij * cf2;
	    	cf3 = tij3 + rij * cf3;
	    	ic1 = (cf1 - xo1) * cbai;
	    	ic2 = (cf2 - xo2) * cbai;
	    	ic3 = (cf3 - xo3) * cbai;
	    	i__3 = cbn2[ic1 + (ic2 + ic3 * MXG) * MXG];
	    	for (ii = cbn1[ic1 + (ic2 + ic3 * MXG) * MXG]; ii <= i__3; ++ii) {
				k = cbal[ii - 1];
				dy1 = crd[k * 3 + 1] - cf1;
				dy2 = crd[k * 3 + 2] - cf2;
				dy3 = crd[k * 3 + 3] - cf3;
				ds2 = dy1 * dy1 + dy2 * dy2 + dy3 * dy3;
				if (ds2 < rs2[k - 1]) {
		    		oti[iv - 1] = k;
		    		goto L10;
				}
	    	}
	    	++nvo;
	    	oti[iv - 1] = 0;
L10: 		;
		}
		nst = 0;
		for (ie = 1; ie <= nvi; ++ie) {
	    	ia1 = oti[edgv[(ie << 1) - 2] - 1];
	    	ia2 = oti[edgv[(ie << 1) - 1] - 1];
	    	if (ia1 > 0 && ia1 == ia2) continue;
	    	++nst;
	    	st[nst - 1] = ie;
		}
		while (nst > 0) {

	    	ie = st[nst - 1];
	    	--nst;
	    	ia1 = oti[edgv[(ie << 1) - 2] - 1];
	    	ia2 = oti[edgv[(ie << 1) - 1] - 1];
	    	iv = ie + nvi;
	    	cf1 = rm[0] * ver[iv * 3 - 3] + rm[3] * ver[iv * 3 - 2] + rm[6] * 
		    	ver[iv * 3 - 1];
	    	cf2 = rm[1] * ver[iv * 3 - 3] + rm[4] * ver[iv * 3 - 2] + rm[7] * 
		    	ver[iv * 3 - 1];
	    	cf3 = rm[2] * ver[iv * 3 - 3] + rm[5] * ver[iv * 3 - 2] + rm[8] * 
		    	ver[iv * 3 - 1];
	    	cf1 = tij1 + rij * cf1;
	    	cf2 = tij2 + rij * cf2;
	    	cf3 = tij3 + rij * cf3;
	    	if (ia1 != 0) {
	    		dy1 = crd[ia1 * 3 + 1] - cf1;
	    		dy2 = crd[ia1 * 3 + 2] - cf2;
	    		dy3 = crd[ia1 * 3 + 3] - cf3;
	    		ds2 = dy1 * dy1 + dy2 * dy2 + dy3 * dy3;
	    		if (ds2 < rs2[ia1 - 1]) {
					oti[iv - 1] = ia1;
					if (edg[ie - 1] > 0) {
		    			++nst;
		    			st[nst - 1] = edg[ie - 1] + 1;
					}
					continue;
	    		}
			}
	    	if (ia2 != 0) {
	    		dy1 = crd[ia2 * 3 + 1] - cf1;
	    		dy2 = crd[ia2 * 3 + 2] - cf2;
	    		dy3 = crd[ia2 * 3 + 3] - cf3;
		 		ds2 = dy1 * dy1 + dy2 * dy2 + dy3 * dy3;
	    		if (ds2 < rs2[ia2 - 1]) {
					oti[iv - 1] = ia2;
					if (edg[ie - 1] > 0) {
		    			++nst;
		    			st[nst - 1] = edg[ie - 1];
					}
					continue;
	    		}
			}
	    	ic1 = (cf1 - xo1) * cbai;
	    	ic2 = (cf2 - xo2) * cbai;
	    	ic3 = (cf3 - xo3) * cbai;
	    	i__2 = cbn2[ic1 + (ic2 + ic3 * MXG) * MXG];
	    	for (ii = cbn1[ic1 + (ic2 + ic3 * MXG) * MXG]; ii <= i__2; ++ii) {
				k = cbal[ii - 1];
				dy1 = crd[k * 3 + 1] - cf1;
				dy2 = crd[k * 3 + 2] - cf2;
				dy3 = crd[k * 3 + 3] - cf3;
				ds2 = dy1 * dy1 + dy2 * dy2 + dy3 * dy3;
				if (ds2 < rs2[k - 1]) {
		    		oti[iv - 1] = k;
		    		if (edg[ie - 1] > 0) {
						++nst;
						st[nst - 1] = edg[ie - 1];
						++nst;
						st[nst - 1] = edg[ie - 1] + 1;
		    		}
		    		goto L60;
				}
	    	}
	    	++nvo;
	    	oti[iv - 1] = 0;
	    	if (edg[ie - 1] > 0) {
				++nst;
				st[nst - 1] = edg[ie - 1];
				++nst;
				st[nst - 1] = edg[ie - 1] + 1;
	    	}
L60:		;
		}

		if (nvo > 0) {
	    	nxv += nvo;
	    	++nprx;
	    	frc = nvo * hfact / d2 * (r0[i] * (d2 + r02[j - 1] - r02[i - 1]) 
		    	+ r0[j] * (d2 + r02[i - 1] - r02[j - 1]));
	    	frc1 = frc * dx1 / dmg;
	    	frc2 = frc * dx2 / dmg;
	    	frc3 = frc * dx3 / dmg;
	    	fat[i * 3 + 1] += frc1;
	    	fat[i * 3 + 2] += frc2;
	    	fat[i * 3 + 3] += frc3;
	    	fat[j * 3 + 1] -= frc1;
	    	fat[j * 3 + 2] -= frc2;
	    	fat[j * 3 + 3] -= frc3;
		}
    }

	free( cbn1 );
	free( cbn2 );
	free( pls );
	free( ver );
	free( oti );
	free( cbal );
	free( r02 );
	free( rs2 );

    return 0;
} /* sasad_ */

#define MXCB 40

static	int	cube( REAL_T fac, int natm,
	REAL_T *crd, REAL_T *r0, REAL_T *cbai, REAL_T *xo, REAL_T *yo, REAL_T *zo,
	int *cbn1, int *cbn2, int *cbal )
/*
REAL_T fac;
int natm;
REAL_T *crd, *r0, *cbai, *xo, *yo, *zo;
int *cbn1, *cbn2, *cbal;
*/
{
    int i__4;
    REAL_T cbln;
    int iatm, *cbal2;
	int i, j, k, m, n;
    REAL_T x, y, z;
    int nmxcb, ic, il, im, in, ix, iy, iz, jx, jy, jz;
    REAL_T xp, yp, zp, bl1, bl2, bl3;
    int *cbn;
    REAL_T off, fac2;

    /* Parameter adjustments */
    --cbal;
    --r0;
    crd -= 4;

	/*  allocations  */
	if( ( cbal2 = ( int * )malloc( MXCB*MXG*MXG*MXG*sizeof(int) ) ) == NULL ){
		fprintf( stderr, "Unable to allocate space for cbal2: %d\n", MXCB*MXG*MXG*MXG );
		exit( 1 );
	}
	if( ( cbn = ( int * )malloc( MXG*MXG*MXG*sizeof(int) ) ) == NULL ){
		fprintf( stderr, "Unable to allocate space for cbn: %d\n", MXG*MXG*MXG );
		exit( 1 );
	}

    cbln = (REAL_T)0.;
    for (iatm = 1; iatm <= natm; ++iatm) {
		if (r0[iatm] > cbln) cbln = r0[iatm];
    }
    cbln = fac * cbln;
    fac2 = (REAL_T)2.;
    if (fac == (REAL_T)2.) fac2 = (REAL_T)1.;
    *xo = *yo = *zo = 1000.;
    xp = yp = zp = -1000.;
    off = 0.1;

    for (iatm = 1; iatm <= natm; ++iatm) {
		if (crd[iatm * 3 + 1] < *xo) *xo = crd[iatm * 3 + 1];
		if (crd[iatm * 3 + 1] >  xp)  xp = crd[iatm * 3 + 1];
		if (crd[iatm * 3 + 2] < *yo) *yo = crd[iatm * 3 + 2];
		if (crd[iatm * 3 + 2] >  yp)  yp = crd[iatm * 3 + 2];
		if (crd[iatm * 3 + 3] < *zo) *zo = crd[iatm * 3 + 3];
		if (crd[iatm * 3 + 3] >  zp)  zp = crd[iatm * 3 + 3];
    }

    *xo = *xo - fac2 * cbln - off;
    *yo = *yo - fac2 * cbln - off;
    *zo = *zo - fac2 * cbln - off;
    xp = xp + fac2 * cbln + off;
    yp = yp + fac2 * cbln + off;
    zp = zp + fac2 * cbln + off;
    bl1 = xp - *xo;
    bl2 = yp - *yo;
    bl3 = zp - *zo;
    il = bl1 / cbln;
    im = bl2 / cbln;
    in = bl3 / cbln;
    *cbai = (REAL_T)1. / cbln;

    for (i = 0; i <= il; ++i) {
		for (j = 0; j <= im; ++j) {
	    	for (k = 0; k <= in; ++k) {
				cbn[i + (j + k * MXG) * MXG] = 0;
				cbn1[i + (j + k * MXG) * MXG] = 1;
				cbn2[i + (j + k * MXG) * MXG] = 0;
	    	}
		}
    }

    for (i = 1; i <= natm; ++i) {
		x = (crd[i * 3 + 1] - *xo) / cbln;
		ix = (int) x;
		y = (crd[i * 3 + 2] - *yo) / cbln;
		iy = (int) y;
		z = (crd[i * 3 + 3] - *zo) / cbln;
		iz = (int) z;
		++cbn[ix + (iy + iz * MXG) * MXG];
		cbal2[cbn[ix + (iy + iz * MXG) * MXG] + (ix + (iy + iz * MXG) * MXG)
* MXCB 
			- 1] = i;
    }
    n = 0;
    nmxcb = 0;
    for (jx = 1; jx <= il-1; ++jx) {
	for (jy = 1; jy <= im-1; ++jy) {
	    for (jz = 1; jz <= in-1; ++jz) {
		if (cbn[jx + (jy + jz * MXG) * MXG] > nmxcb) {
		    nmxcb = cbn[jx + (jy + jz * MXG) * MXG];
		}
		if (cbn[jx + (jy + jz * MXG) * MXG] > MXCB) {
		    fprintf( stderr, "cubing neighbor limit mxcb exceeded: %d\n", MXCB);
			exit(1);
		}
		m = 0;
		ix = jx;
		iy = jy;
		iz = jz;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* -1,0,0 */
		--ix;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 1,0,0 */
		ix += 2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 0,-1,0 */
		--ix;
		--iy;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 0,1,0 */
		iy += 2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 0,0,-1 */
		--iy;
		--iz;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 0,0,1 */
		iz += 2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* nn=2 */
/* 1,0,1 */
		++ix;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* -1,0,1 */
		ix += -2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 0,1,1 */
		++ix;
		++iy;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 0,-1,1 */
		iy += -2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* -1,-1,0 */
		--iz;
		--ix;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 1,-1,0 */
		ix += 2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 1,1,0 */
		iy += 2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* -1,1,0 */
		ix += -2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* -1,0,-1 */
		--iz;
		--iy;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 1,0,-1 */
		ix += 2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 0,1,-1 */
		--ix;
		++iy;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 0,-1,-1 */
		iy += -2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* nn=3 */
/* -1,-1,-1 */
		--ix;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 1,-1,-1 */
		ix += 2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 1,1,-1 */
		iy += 2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* -1,1,-1 */
		ix += -2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* -1,1,1 */
		iz += 2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 1,1,1 */
		ix += 2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* 1,-1,1 */
		iy += -2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
/* -1,-1,1 */
		ix += -2;
		i__4 = cbn[ix + (iy + iz * MXG) * MXG];
		for (ic = 1; ic <= i__4; ++ic) {
		    ++m;
		    cbal[n + m] = cbal2[ic + (ix + (iy + iz * MXG) * MXG) * MXCB 
			    - 1];
		}
		if (m > 0) {
		    cbn1[jx + (jy + jz * MXG) * MXG] = n + 1;
		    cbn2[jx + (jy + jz * MXG) * MXG] = n + m;
		    n += m;
		}
	    }
	}
    }
    fprintf( stderr, "max cubing neighbors: %d\n", nmxcb );
    if (n > MXACUM) {
		fprintf( stderr, "max number of cubing neighbors exceeeded: %d\n ",
			MXACUM );
		exit(1);
    }

	free( cbal2 );
	free( cbn );
    return 0;
} /* cube_ */
