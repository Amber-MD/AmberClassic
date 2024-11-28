
/* rattle and rattle2 come more */
/* or less  (rather more) from */
/* the tinker package */
/* ($tinkerhome/source/rattle.f) */

#define FAC 1.2
#define EPS 0.0000001

int rattle(REAL_T, REAL_T *, REAL_T *, REAL_T *, REAL_T *);
int rattle2(REAL_T, REAL_T *, REAL_T *, REAL_T *);

int rattle(REAL_T dt, REAL_T * x, REAL_T * xold, REAL_T * v, REAL_T * minv)
/*
REAL_T dt;
REAL_T *x, *xold, *v, *minv;
*/
{
   int maxiter;
   int niter, i;
   int done;
   float eps;
   int nratH, nrat;
   int at1, at2;
   int atyp;
   REAL_T xr, yr, zr, xo, yo, zo, xterm, yterm, zterm;
   REAL_T dist2, dot;
   REAL_T krat, term, delta, rma, rmb;

   eps = EPS;
   maxiter = 1000;
   niter = 0;
   done = 0;
   nratH = prm->Nbonh;
   nrat = prm->Mbona;
   if( irattle == 2 ) nrat=0;

   while ((done == 0) && (niter < maxiter)) {
      niter++;
      done = 1;

      for (i = 0; i < nratH + nrat; i++) {
         if (i < nratH) {
            at1 = prm->BondHAt1[i];
            at2 = prm->BondHAt2[i];
            atyp = prm->BondHNum[i] - 1;
            krat = prm->Req[atyp];
         } else {
            at1 = prm->BondAt1[i - nratH];
            at2 = prm->BondAt2[i - nratH];
            atyp = prm->BondNum[i - nratH] - 1;
            krat = prm->Req[atyp];
         }

         xr = x[at2] - x[at1];
         yr = x[at2 + 1] - x[at1 + 1];
         zr = x[at2 + 2] - x[at1 + 2];
         dist2 = xr * xr + yr * yr + zr * zr;

         delta = krat * krat - dist2;

         if (fabs(delta) > eps) {
            done = 0;
            xo = xold[at2] - xold[at1];
            yo = xold[at2 + 1] - xold[at1 + 1];
            zo = xold[at2 + 2] - xold[at1 + 2];

            dot = xr * xo + yr * yo + zr * zo;
            rma = minv[at1];
            rmb = minv[at2];

#define FAC 1.2
            term = FAC * delta / (2.0 * (rma + rmb) * dot);

            xterm = xo * term;
            yterm = yo * term;
            zterm = zo * term;

            x[at1] -= xterm * rma;
            x[at1 + 1] -= yterm * rma;
            x[at1 + 2] -= zterm * rma;

            x[at2] += xterm * rmb;
            x[at2 + 1] += yterm * rmb;
            x[at2 + 2] += zterm * rmb;

            rma = rma / dt;
            rmb = rmb / dt;
            v[at1] -= xterm * rma;
            v[at1 + 1] -= yterm * rma;
            v[at1 + 2] -= zterm * rma;

            v[at2] += xterm * rmb;
            v[at2 + 1] += yterm * rmb;
            v[at2 + 2] += zterm * rmb;
         }
      }
   }

   if (niter >= maxiter - 1){
      fprintf(nabout, "ERROR in RATTLE\n");
      exit(1);
   }

   return (niter);
}

int rattle2(REAL_T dt, REAL_T * x, REAL_T * v, REAL_T * minv)
/*
REAL_T dt;
REAL_T *x, *v, *minv;
*/
{
   int maxiter, niter, done, i;
   REAL_T eps, xr, yr, zr, xv, yv, zv, xterm, yterm, zterm;
   REAL_T dot, rma, rmb, krat, term;
   int nratH, nrat;
   int at1, at2, atyp;

   eps = EPS / dt;
   maxiter = 1000;
   niter = 0;
   nratH = prm->Nbonh;
   nrat = prm->Mbona;
   if( irattle == 2 ) nrat=0;
   done = 0;

   while ((done == 0) && (niter < maxiter)) {
      niter++;
      done = 1;

      for (i = 0; i < nratH + nrat; i++) {
         if (i < nratH) {
            at1 = prm->BondHAt1[i];
            at2 = prm->BondHAt2[i];
            atyp = prm->BondHNum[i] - 1;
            krat = prm->Req[atyp];
         } else {
            at1 = prm->BondAt1[i - nratH];
            at2 = prm->BondAt2[i - nratH];
            atyp = prm->BondNum[i - nratH] - 1;
            krat = prm->Req[atyp];
         }

         xr = x[at2] - x[at1];
         yr = x[at2 + 1] - x[at1 + 1];
         zr = x[at2 + 2] - x[at1 + 2];

         xv = v[at2] - v[at1];
         yv = v[at2 + 1] - v[at1 + 1];
         zv = v[at2 + 2] - v[at1 + 2];

         dot = xr * xv + yr * yv + zr * zv;
         rma = minv[at1];
         rmb = minv[at2];

         term = -dot * FAC / ((rma + rmb) * krat * krat);

         if (fabs(term) > eps) {
            done = 0;
            xterm = xr * term;
            yterm = yr * term;
            zterm = zr * term;

            v[at1] -= xterm * rma;
            v[at1 + 1] -= yterm * rma;
            v[at1 + 2] -= zterm * rma;

            v[at2] += xterm * rmb;
            v[at2 + 1] += yterm * rmb;
            v[at2 + 2] += zterm * rmb;
         }

      }
   }
   if (niter >= maxiter - 1){
      fprintf(nabout, "ERROR in RATTLE2\n");
      exit(1);
   }
   return (niter);
}
