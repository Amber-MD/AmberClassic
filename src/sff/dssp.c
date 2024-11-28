/* dsspdac.f -- translated by f2c (version 20030306).
*/

/* Common Block Declarations */

struct {
    REAL_T v0, wv, ron0, won, q;
} dsspvars;


REAL_T fdssp(REAL_T * x, REAL_T * f)
{
    REAL_T edssp;
    /* System generated locals */
    int i1, i2;
    REAL_T d1, d2, d3;

    /* Builtin functions */
    REAL_T sqrt(REAL_T), exp(REAL_T);

    /* Local variables */
    static int i, j, k;
    static REAL_T u, v, fr, xc[3], fv, xh[3], xn[3], xo[3], rch, xch[3],
        rcn, xcn[3], roh, xoh[3], ron, xon[3], rch2, rch3, rcn2, rcn3,
        roh2, roh3, ron2, ron3;
    static int ipepc[500];
    static REAL_T dfrdr;
    static int npepc;
    static REAL_T dvdxc[3], dfvdv, dvdxh[3], drdxo[3], drdxn[3], dvdxn[3],
        dvdxo[3], constf;
    extern int udssp(REAL_T *, REAL_T *, REAL_T *, REAL_T *, REAL_T *,
                     REAL_T *, REAL_T *, REAL_T *);

    /* Parameter adjustments */
    --f;
    --x;

    /* Function Body */
    edssp = 0.0;
    dsspvars.v0 = 1.;
    dsspvars.wv = 1.;
    dsspvars.ron0 = 1.;
    dsspvars.won = 1.;
    dsspvars.q = 80.;
/*        loop over all pairs of peptide groups */
    i1 = npepc - 1;
    for (i = 1; i <= i1; ++i) {
/*       get the coordinates needed for the function: carbonyl of peptide "i" */
        for (k = 1; k <= 3; ++k) {
            xc[k - 1] = x[ipepc[i - 1] * 3 - 3 + k];
            xo[k - 1] = x[(ipepc[i - 1] + 1) * 3 - 3 + k];
        }
        i2 = npepc;
        for (j = i + 1; j <= i2; ++j) {
/*              NH group of peptide "j" */
            for (k = 1; k <= 3; ++k) {
                xn[k - 1] = x[(ipepc[j - 1] + 2) * 3 - 3 + k];
                xh[k - 1] = x[(ipepc[j - 1] + 3) * 3 - 3 + k];
                xcn[k - 1] = xc[k - 1] - xn[k - 1];
                xch[k - 1] = xc[k - 1] - xh[k - 1];
                xon[k - 1] = xo[k - 1] - xn[k - 1];
                xoh[k - 1] = xo[k - 1] - xh[k - 1];
            }
            d1 = xcn[0];
            d2 = xcn[1];
            d3 = xcn[2];
            rcn2 = d1 * d1 + d2 * d2 + d3 * d3;
            d1 = xch[0];
            d2 = xch[1];
            d3 = xch[2];
            rch2 = d1 * d1 + d2 * d2 + d3 * d3;
            d1 = xon[0];
            d2 = xon[1];
            d3 = xon[2];
            ron2 = d1 * d1 + d2 * d2 + d3 * d3;
            d1 = xoh[0];
            d2 = xoh[1];
            d3 = xoh[2];
            roh2 = d1 * d1 + d2 * d2 + d3 * d3;
            rcn = sqrt(rcn2);
            rch = sqrt(rch2);
            ron = sqrt(ron2);
            roh = sqrt(roh2);
            d1 = rcn;
            rcn3 = d1 * (d1 * d1);
            d1 = rch;
            rch3 = d1 * (d1 * d1);
            d1 = ron;
            ron3 = d1 * (d1 * d1);
            d1 = roh;
            roh3 = d1 * (d1 * d1);
            udssp(&ron, &rch, &roh, &rcn, &u, &v, &fv, &fr);
            edssp += u;

/*   update the forces:  compute DSSP force on all atoms in the peptide pair */
/*   (i) C-O...H-N (j) */
/*   since U = V * F(V) * F(ron), dU/dx = ( dV/dx ) * F(V) * F(ron) + */
/*             V * { [ dF(V)/dV ] * [ dV/dx ] } * F(ron) + */
/*             V * F(V) * [ dF(ron)/dx ] */
/*             add V * [ dF(V)/dx ] * F(r) term */
/*             evaluate dV/dx [ = ( dV/dr ) * ( dr/dx ) ] */

            for (k = 1; k <= 3; ++k) {
                dvdxc[k - 1] =
                    dsspvars.q * (-xch[k - 1] / rch3 + xcn[k - 1] / rcn3);
                dvdxo[k - 1] =
                    dsspvars.q * (-xon[k - 1] / ron3 + xoh[k - 1] / roh3);
                dvdxn[k - 1] =
                    dsspvars.q * (xon[k - 1] / ron3 - xcn[k - 1] / rcn3);
                dvdxh[k - 1] =
                    dsspvars.q * (xch[k - 1] / rch3 - xoh[k - 1] / roh3);
            }

/*      evaluate dF(V)/dV */

            dfvdv = exp((v - dsspvars.v0) / dsspvars.wv) / dsspvars.wv;
            dfvdv = -dfvdv * fv * fv;
            constf = v * dfvdv * fr;
            for (k = 1; k <= 3; ++k) {
                f[ipepc[i - 1] * 3 - 3 + k] -= constf * dvdxc[k - 1];
                f[(ipepc[i - 1] + 1) * 3 - 3 + k] -= constf * dvdxo[k - 1];
                f[(ipepc[j - 1] + 2) * 3 - 3 + k] -= constf * dvdxn[k - 1];
                f[(ipepc[j - 1] + 3) * 3 - 3 + k] -= constf * dvdxh[k - 1];
            }
/*    add V * F(V) * [ dF(r)/dx ] term:  since F(r) depends only on ron, */
/*    only xo and xn derivatives of ron are needed */
/*    evaluate dF(r)/dr and dr/dx */

            dfrdr =
                exp((ron - dsspvars.ron0) / dsspvars.won) / dsspvars.won;
            dfrdr = -dfrdr * fr * fr;
            for (k = 1; k <= 3; ++k) {
                drdxo[k - 1] = xon[k - 1] / ron;
                drdxn[k - 1] = -xon[k - 1] / ron;
            }
            constf = v * fv * dfrdr;
            for (k = 1; k <= 3; ++k) {
                f[(ipepc[i - 1] + 1) * 3 - 3 + k] -= constf * drdxo[k - 1];
                f[(ipepc[j - 1] + 2) * 3 - 3 + k] -= constf * drdxn[k - 1];
            }

/* add (dV/dx) * F(V) * F(r) term */

            constf = fv * fr;
            for (k = 1; k <= 3; ++k) {
                f[ipepc[i - 1] * 3 - 3 + k] -= constf * dvdxc[k - 1];
                f[(ipepc[i - 1] + 1) * 3 - 3 + k] -= constf * dvdxo[k - 1];
                f[(ipepc[j - 1] + 2) * 3 - 3 + k] -= constf * dvdxn[k - 1];
                f[(ipepc[j - 1] + 3) * 3 - 3 + k] -= constf * dvdxh[k - 1];
            }
        }
    }
    return edssp;
}

int udssp(REAL_T * ron, REAL_T * rch, REAL_T * roh,
          REAL_T * rcn, REAL_T * u, REAL_T * v, REAL_T * fv, REAL_T * fr)
{
    extern int fermi(REAL_T *, REAL_T *, REAL_T *, REAL_T *),
        vdssp(REAL_T *, REAL_T *, REAL_T *, REAL_T *, REAL_T *);

/* DSSP potential; separately evaluates cutoff (Fermi) functions in V and ron */
/*     F(V) and F(ron) */

    vdssp(ron, rch, roh, rcn, v);
    fermi(v, &dsspvars.v0, &dsspvars.wv, fv);
    fermi(ron, &dsspvars.ron0, &dsspvars.won, fr);
    *u = *v * *fv * *fr;
    return 0;

}

int vdssp(REAL_T * ron, REAL_T * rch, REAL_T * roh,
          REAL_T * rcn, REAL_T * v)
{

/*     DSSP potential (aside from cutoff functions in V and ron) */
/*     V = 332 * Qc * Qn * ( 1/ron + 1/rch - 1/roh - 1/rcn ) */

    *v = 1. / *ron + 1. / *rch - 1. / *roh - 1. / *rcn;
    *v = dsspvars.q * *v;
    return 0;
}

int fermi(REAL_T * x, REAL_T * x0, REAL_T * w, REAL_T * f)
{
    /* Builtin functions */
    REAL_T exp(REAL_T);

/*     Fermi function: { exp[ (x-x0)/w ] + 1 }^-1 */

    if (*w > 0.) {
        *f = exp((*x - *x0) / *w) + 1.;
        *f = 1. / *f;
    } else if (*w <= 0.) {
        if (*x <= *x0) {
            *f = 1.;
        } else if (*x > *x0) {
            *f = 0.;
        }
    }
    return 0;
}
