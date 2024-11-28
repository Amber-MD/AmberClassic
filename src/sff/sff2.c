/*
 * sff2.c:  simple force field: implement routines to read a "prmtop"
 * file from AMBER, and calculate the energy.  Implements only
 * some of AMBER's functionality: bonds, angles, dihedrals, and
 * nonbonded interactions with a distance-dependent dielectric.
 * Does not (yet) include support for period systems.
 *
 * Main interface is through routines "mme2" and "mme_init".
 *
 * Second derivatives, in particular the egb2 function, were
 * implemented by Russ Brown (russ.brown@sun.com).
 *
 * Changes to difang, eangl2 and ephi2 by David A. Case (case@scripps.edu).
 *
 * For the OpenMP implementation, only the egb2 function is parallelized
 * so the econs2, ebond2, eangl2, ephi2 and nbond2 functions execute in
 * single-threaded mode.
 *
 * For the ScaLAPACK implementation, only the egb2 function is parallelized
 * so the econs2, ebond2, eangl2, ephi2 and nbond2 functions all compute
 * each gradient vector and Hessian matrix element.  However, because the
 * gradient and Hessian are distributed, the ptr1d, ptr1f, ptr2d and ptr2f
 * functions are used to write the subvector and submatrix elements by only
 * the process that owns the subvector or submatrix.  The SUPERD_T type is
 * used to facilitate argument transmission to the ptr2f function.
 */

#ifdef SCALAPACK

typedef struct superd {
   INT_T dtype, ctxt, m, n, mb, nb, rsrc, csrc, lld, row, col;
   size_t size;
   INT_T *divmb, *divnb, *modmb, *modnb, *divrw, *divcl, *modrw, *modcl;
} SUPERD_T;

#endif

int numroc_(int *, int *, int *, int *, int *);

void descinit_(int *, int *, int *, int *, int *,
               int *, int *, int *, int *, int *);

void pdgemr2d_(INT_T *, INT_T *,
               REAL_T *, INT_T *, INT_T *, INT_T *,
               REAL_T *, INT_T *, INT_T *, INT_T *, INT_T *);

void pdgemm_(char *, char *, INT_T *, INT_T *, INT_T *, REAL_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *, REAL_T *,
             REAL_T *, INT_T *, INT_T *, INT_T *);

void dgemm_(char *, char *, INT_T *, INT_T *, INT_T *, REAL_T *, REAL_T *,
            INT_T *, REAL_T *, INT_T *, REAL_T *, REAL_T *, INT_T *);

#ifdef SCALAPACK

/***********************************************************************
                            SUPERDESC()
************************************************************************/

/* Create and initialize a superdescriptor. */

SUPERD_T *superdesc(INT_T * desc)
{
   SUPERD_T *sd;
   int nprow, npcol, myrow, mycol, i, zero = 0;
   size_t locp, locq;

   /* Malloc a superdescriptor and copy values from the descriptor. */

   if ((sd = (SUPERD_T *) malloc(sizeof(SUPERD_T))) == NULL) {
      printf("Error allocating superdescriptor in superdesc!\n");
      fflush(stdout);
      mpierror(-1);
   }

   sd->dtype = desc[DTYPE_];
   sd->ctxt = desc[CTXT_];
   sd->m = desc[M_];
   sd->n = desc[N_];
   sd->mb = desc[MB_];
   sd->nb = desc[NB_];
   sd->rsrc = desc[RSRC_];
   sd->csrc = desc[CSRC_];
   sd->lld = desc[LLD_];

   if (gb2_debug) {
      printf("superdesc: context = %d\n", sd->ctxt);
      fflush(stdout);
   }

   /* Initialize the division and modulus array pointers to NULL. */

   sd->divmb = sd->modmb = sd->divrw = sd->modrw = NULL;
   sd->divnb = sd->modnb = sd->divcl = sd->modcl = NULL;

   /* Return NULL if the context is invalid. */

   if (sd->lld < 0)
      return (NULL);

   /* 
    * Calculate the maximum size of the submatrix array.
    * Store the row and column for this process.
    */

   blacs_gridinfo_(&(sd->ctxt), &nprow, &npcol, &myrow, &mycol);
   locp = numroc_(&(sd->m), &(sd->mb), &myrow, &zero, &nprow);
   locq = numroc_(&(sd->n), &(sd->nb), &mycol, &zero, &npcol);

   sd->size = locp * locq;
   sd->row = myrow;
   sd->col = mycol;

   /* Return NULL if the size or the matrix dimensions are invalid. */

   if (sd->size < 0 || sd->m < 0 || sd->n < 0)
      return (NULL);

   /* Create and fill the division and modulus arrays. */

   sd->divmb = ivector(0, sd->m);
   sd->modmb = ivector(0, sd->m);
   sd->divrw = ivector(0, sd->m);
   sd->modrw = ivector(0, sd->m);

   sd->divnb = ivector(0, sd->n);
   sd->modnb = ivector(0, sd->n);
   sd->divcl = ivector(0, sd->n);
   sd->modcl = ivector(0, sd->n);

   for (i = 0; i < sd->m; i++) {
      (sd->divmb)[i] = i / (sd->mb);
      (sd->modmb)[i] = i % (sd->mb);
      (sd->divrw)[i] = i / nprow;
      (sd->modrw)[i] = i % nprow;
   }

   for (i = 0; i < sd->n; i++) {
      (sd->divnb)[i] = i / (sd->nb);
      (sd->modnb)[i] = i % (sd->nb);
      (sd->divcl)[i] = i / npcol;
      (sd->modcl)[i] = i % npcol;
   }

   /* Return the superdescriptor. */

   return (sd);

}

/***********************************************************************
                            FREE_SUPERDESC()
************************************************************************/

/* Free a superdescriptor including all of its dynamic arrays. */

void free_superdesc(SUPERD_T * sd)
{
   if (sd == NULL)
      return;

   if (sd->divmb != NULL)
      free_ivector(sd->divmb, 0, sd->m);
   if (sd->modmb != NULL)
      free_ivector(sd->modmb, 0, sd->m);
   if (sd->divrw != NULL)
      free_ivector(sd->divrw, 0, sd->m);
   if (sd->modrw != NULL)
      free_ivector(sd->modrw, 0, sd->m);

   if (sd->divnb != NULL)
      free_ivector(sd->divnb, 0, sd->n);
   if (sd->modnb != NULL)
      free_ivector(sd->modnb, 0, sd->n);
   if (sd->divcl != NULL)
      free_ivector(sd->divcl, 0, sd->n);
   if (sd->modcl != NULL)
      free_ivector(sd->modcl, 0, sd->n);

   free(sd);
}

/***********************************************************************
                            ADR1D()
************************************************************************/

/* Perform the 1D address mapping from the global to the local address. */

static
size_t adr1d(int i, int mb, int nprow)
{
   return ((i / (nprow * mb)) * mb + i % mb);
}

/***********************************************************************
                            ADR2D()
************************************************************************/

/* Perform the 2D address mapping from the global to the local address. */

static
size_t adr2d(int i, int j, int mb, int nb, int nprow, int npcol, int lld)
{
   size_t lbi, lbj, bai, baj;

   /* Calculate the local block address (lbi, lbj) */

   lbi = i / (nprow * mb);
   lbj = j / (npcol * nb);

   /* Calculate the offset (bai, baj) within the local block. */

   bai = i % mb;
   baj = j % nb;

   /* Calculate the local array address for column-major addressing. */

   return (lld * (nb * lbj + baj) + mb * lbi + bai);
}

/***********************************************************************
                            PTR1D()
************************************************************************/

/* Return a pointer to a submatrix element or NULL if undefined. */

REAL_T *ptr1d(REAL_T * a, int *desc, int i)
{
   int nprow, npcol, myrow, mycol, locp, locq, zero = 0;
   size_t adr;

   /* Return NULL if the matrix isn't active in the context. */

   if (desc[CTXT_] < 0)
      return (NULL);

   /* Return NULL if the global index is out of range. */

   if (i > desc[M_]) {
      printf("Global vector address overflow, i=%d  imax=%d  context=%d\n",
             i, desc[M_], desc[CTXT_]);
      fflush(stdout);
      return (NULL);
   }

   /* Get the process grid information. */

   blacs_gridinfo_(&desc[CTXT_], &nprow, &npcol, &myrow, &mycol);

   /*
    * Return NULL if the submatrix doesn't exist in this process.
    * Perform the column test first because the row test is likely
    * to be true for egb2 because a row test is performed at the
    * level of the loop index.
    */

   if (mycol != 0 || !myroc(i, desc[MB_], nprow, myrow))
      return (NULL);

   /* Calculate the submatrix address and return NULL if out of range. */

   locp = numroc_(&desc[M_], &desc[MB_], &myrow, &zero, &nprow);
   locq = numroc_(&desc[N_], &desc[NB_], &mycol, &zero, &npcol);
   adr = adr1d(i, desc[MB_], nprow);
   if (adr > (size_t) (locp * locq)) {
      printf
          ("Local subvector address overflow, adr=%lld  adrmax=%lld  context=%d\n",
           adr, (size_t) (locp * locq), desc[CTXT_]);
      fflush(stdout);
      return (NULL);
   }

   /* If the vector array isn't allocated return NULL. */

   if (a == NULL) {
      printf("Local subvector is not allocated, context=%d\n",
             desc[CTXT_]);
      fflush(stdout);
      return (NULL);
   }

   /* Otherwise, return the address of the matrix element. */

   return (&a[adr]);
}

/***********************************************************************
                            PTR1F()
************************************************************************/

/* Return a pointer to a submatrix element or NULL if undefined. */

REAL_T *ptr1f(REAL_T * a, SUPERD_T * sd, int i)
{
   int lbi, bai, adr;

   /* Return NULL if the superdescriptor is NULL. */

   if (sd == NULL)
      return (NULL);

   /* Return NULL if the matrix isn't active in the context. */

   if (sd->ctxt < 0)
      return (NULL);

   /* Return NULL if the global index is out of range. */

   if (i < 0 || i > sd->m) {
      printf
          ("Global vector address overflow, i=%d  imax=%lld  context=%d\n",
           i, sd->m, sd->ctxt);
      fflush(stdout);
      return (NULL);
   }

   /*
    * Return NULL if the submatrix doesn't exist in this process.
    * Perform the column test first because the row test is likely
    * to be true for egb2 because a row test is performed at the
    * level of the loop index.
    */

   lbi = (sd->divmb)[i];
   if (sd->col != 0 || (sd->modrw)[lbi] != sd->row)
      return (NULL);

   /* Calculate the submatrix address and return NULL if out of range. */

   lbi = (sd->divrw)[lbi];
   bai = (sd->modmb)[i];
   adr = (sd->mb) * lbi + bai;
   if (adr < 0 || adr > sd->size) {
      printf
          ("Local subvector address overflow, adr=%lld  adrmax=%lld  context=%d\n",
           adr, sd->size, sd->ctxt);
      fflush(stdout);
      return (NULL);
   }

   /* If the vector array isn't allocated return NULL. */

   if (a == NULL) {
      printf("Local subvector is not allocated, context=%d\n", sd->ctxt);
      fflush(stdout);
      return (NULL);
   }

   /* Otherwise, return the address of the matrix element. */

   return (&a[adr]);
}

/***********************************************************************
                            PTR2D()
************************************************************************/

/* Return a pointer to a submatrix element or NULL if undefined. */

REAL_T *ptr2d(REAL_T * a, int *desc, int i, int j)
{
   int nprow, npcol, myrow, mycol, locp, locq, zero = 0;
   size_t adr;

   /* Return NULL if the matrix isn't active in the context. */

   if (desc[CTXT_] < 0)
      return (NULL);

   /* Return NULL if the global indices are out of range. */

   if (i > desc[M_]) {
      printf
          ("Global matrix address overflow, i=%d  imax=%lld  context=%d\n",
           i, desc[M_], desc[CTXT_]);
      fflush(stdout);
      return (NULL);
   }
   if (j > desc[N_]) {
      printf
          ("Global matrix address overflow, j=%d  jmax=%lld  context=%d\n",
           j, desc[N_], desc[CTXT_]);
      fflush(stdout);
      return (NULL);
   }

   /* Get the process grid information. */

   blacs_gridinfo_(&desc[CTXT_], &nprow, &npcol, &myrow, &mycol);

   /* Return NULL if the submatrix doesn't exist in this process. */

   if (!myroc(i + desc[RSRC_], desc[MB_], nprow, myrow) ||
       !myroc(j + desc[CSRC_], desc[NB_], npcol, mycol))
      return (NULL);

   /* Calculate the submatrix address and return NULL if out of range. */

   locp = numroc_(&desc[M_], &desc[MB_], &myrow, &zero, &nprow);
   locq = numroc_(&desc[N_], &desc[NB_], &mycol, &zero, &npcol);
   adr = adr2d(i, j, desc[MB_], desc[NB_], nprow, npcol, locp);
   if (adr > (size_t) (locp * locq)) {
      printf
          ("Local submatrix address overflow, adr=%lld  adrmax=%lld  context=%d\n",
           adr, (size_t) (locp * locq), desc[CTXT_]);
      fflush(stdout);
      return (NULL);
   }

   /* If the matrix array isn't allocated return NULL. */

   if (a == NULL) {
      printf("Local submatrix is not allocated, context=%d\n",
             desc[CTXT_]);
      fflush(stdout);
      return (NULL);
   }

   /* Otherwise, return the address of the matrix element. */

   return (&a[adr]);
}

/***********************************************************************
                            PTR2F()
************************************************************************/

/*
 * Return a pointer to a submatrix element or NULL if undefined.
 * Always print an error message for 'hard' errors, such as address
 * overflows, but only print an error message for normal NULL,
 * such as 'matrix isn't active in context' if gb2_debug != 0.
 */

REAL_T *ptr2f(REAL_T * a, SUPERD_T * sd, int i, int j, int *info)
{
   size_t lbi, lbj, bai, baj, adr;

   /* Return NULL if the superdescriptor is NULL. */

   if (sd == NULL) {
      if (gb2_debug) {
         printf("ptr2f: superdescriptor is NULL\n");
         fflush(stdout);
      }
      *info = -1;
      return (NULL);
   }

   /* Return NULL if the matrix isn't active in the context. */

   if (sd->ctxt < 0) {
      if (gb2_debug) {
         printf("ptr2f: matrix isn't active in context %d\n", sd->ctxt);
         fflush(stdout);
      }
      *info = -2;
      return (NULL);
   }

   /* Return NULL if the global indices are out of range. */

   if (i < 0) {
      printf
          ("ptr2f: matrix row address underflow, i=%d  imax=%d  context=%d\n",
           i, sd->m, sd->ctxt);
      fflush(stdout);
      *info = -3;
      return (NULL);
   }
   if (i > sd->m) {
      printf
          ("ptr2f: matrix row address overflow, i=%d  imax=%d  context=%d\n",
           i, sd->m, sd->ctxt);
      fflush(stdout);
      *info = -4;
      return (NULL);
   }
   if (j < 0) {
      printf
          ("ptr2f: matrix column address underflow, j=%d  jmax=%d  context=%d\n",
           j, sd->n, sd->ctxt);
      fflush(stdout);
      *info = -5;
      return (NULL);
   }
   if (j > sd->n) {
      printf
          ("ptr2f: matrix column address overflow, j=%d  jmax=%d  context=%d\n",
           j, sd->n, sd->ctxt);
      fflush(stdout);
      *info = -6;
      return (NULL);
   }

   /*
    * Return NULL if the submatrix doesn't exist in this process.
    * Perform the column test first because the row test is likely
    * to be true for egb2 because a row test is performed at the
    * level of the loop index.
    *
    * The test is (csrc+j/nb)%npcol != mycol || (rsrc+i/mb)%nprow != myrow
    */

   lbi = (sd->divmb)[i];
   lbj = (sd->divnb)[j];
   if ((sd->modcl)[sd->csrc + lbj] != sd->col ||
       (sd->modrw)[sd->rsrc + lbi] != sd->row) {
      if (gb2_debug) {
         printf
             ("ptr2f: submatrix doesn't exist in process, row=%d  col=%d\n",
              sd->row, sd->col);
         fflush(stdout);
      }
      *info = -7;
      return (NULL);
   }

   /* Calculate the local block address (lbi=i/mb/nprow, lbj=j/nb/npcol) */

   lbi = (sd->divrw)[lbi];
   lbj = (sd->divcl)[lbj];

   /* Calculate the offset (bai=i%mb, baj=j%nb) within the local block. */

   bai = (sd->modmb)[i];
   baj = (sd->modnb)[j];

   /* Calculate the submatrix address and return NULL if out of range. */

   adr = (sd->lld) * ((sd->nb) * lbj + baj) + (sd->mb) * lbi + bai;
   if (adr < 0 || adr > sd->size) {
      printf
          ("ptr2f: submatrix address overflow, adr=%lld  adrmax=%lld  context=%d\n",
           adr, sd->size, sd->ctxt);
      fflush(stdout);
      *info = -8;
      return (NULL);
   }

   /* If the matrix array isn't allocated return NULL. */

   if (a == NULL) {
      printf("ptr2f: submatrix is not allocated, context=%d\n", sd->ctxt);
      fflush(stdout);
      *info = -9;
      return (NULL);
   }

   /* Otherwise, return the address of the matrix element. */

   *info = 0;
   return (&a[adr]);
}

#endif

/***********************************************************************
                            ECONS2()
************************************************************************/

/* Calculate the constrained energy and derivatives. */

static
REAL_T econs2(REAL_T * x, REAL_T * f, REAL_T * g,
              int context_PxQ, int *descF_PxQ, int *descG_PxQ)
{
   int i;
   size_t n;
   REAL_T e_cons, rx, ry, rz;

#ifdef SCALAPACK
   REAL_T *ptr;
#endif

   n = 3 * prm->Natom;
   e_cons = 0.0;
   for (i = 0; i < prm->Natom; i++) {
      if (constrained[i]) {
         rx = x[3 * i] - x0[3 * i];
         ry = x[3 * i + 1] - x0[3 * i + 1];
         rz = x[3 * i + 2] - x0[3 * i + 2];
         e_cons += wcons * (rx * rx + ry * ry + rz * rz);

#ifdef SCALAPACK

         ptr = ptr1d(f, descF_PxQ, 3 * i + 0);
         if (ptr != NULL)
            *ptr += 2. * wcons * rx;
         ptr = ptr1d(f, descF_PxQ, 3 * i + 1);
         if (ptr != NULL)
            *ptr += 2. * wcons * ry;
         ptr = ptr1d(f, descF_PxQ, 3 * i + 2);
         if (ptr != NULL)
            *ptr += 2. * wcons * rz;

         ptr = ptr2d(g, descG_PxQ, 3 * i + 0, 3 * i + 0);
         if (ptr != NULL)
            *ptr += 2. * wcons;
         ptr = ptr2d(g, descG_PxQ, 3 * i + 1, 3 * i + 1);
         if (ptr != NULL)
            *ptr += 2. * wcons;
         ptr = ptr2d(g, descG_PxQ, 3 * i + 1, 3 * i + 1);
         if (ptr != NULL)
            *ptr += 2. * wcons;
#else
         f[3 * i] += 2. * wcons * rx;
         f[3 * i + 1] += 2. * wcons * ry;
         f[3 * i + 2] += 2. * wcons * rz;

         g[3 * i + 0 + n * (3 * i + 0)] += 2. * wcons;
         g[3 * i + 1 + n * (3 * i + 1)] += 2. * wcons;
         g[3 * i + 2 + n * (3 * i + 2)] += 2. * wcons;
#endif
      }
   }
   return (e_cons);
}

/***********************************************************************
                            EBOND2()
************************************************************************/

/* Calculate bond stretching energy and derivatives. */

static
REAL_T ebond2(int nbond, int *a1, int *a2, int *atype,
              REAL_T * Rk, REAL_T * Req, REAL_T * x, REAL_T * f,
              REAL_T * g, int context_PxQ, int *descF_PxQ, int *descG_PxQ)
{
   int i, k, l, at1, at2, atyp;
   size_t n;
   REAL_T e_bond, r, r2, rinv, r2inv, rx, ry, rz, dfx, dfy, dfz, db, df,
       dg, e;
   REAL_T di[3], dj[3], d2ii[3][3], d2jj[3][3], d2ij[3][3], ijterm;

#ifdef SCALAPACK
   REAL_T *ptr;
#endif

   n = 3 * prm->Natom;
   e_bond = 0.0;
   for (i = 0; i < nbond; i++) {
      at1 = a1[i];
      at2 = a2[i];
      atyp = atype[i] - 1;

      /* Calculate the interatomic distance and its (x,y,z) components. */

      rx = x[at1] - x[at2];
      ry = x[at1 + 1] - x[at2 + 1];
      rz = x[at1 + 2] - x[at2 + 2];
      r2 = rx * rx + ry * ry + rz * rz;
      r2inv = 1.0 / r2;
      r = sqrt(r2);
      rinv = r * r2inv;

      /*
       * Calculate the first and second derivatives of the interatomic
       * distance Dij with respect to the cartesian coordinates of atoms
       * i and j.  The results are placed into five arrays:
       *
       *   di[] for the first derivatives with respect to atom i
       *   dj[] for the first derivatives with respect to atom j
       *   d2ii[] for the second derivatives with respect to atom i
       *   d2jj[] for the second derivatives with respect to atom j
       *   d2ij[] for the second derivatives with respect to atoms i and j
       *
       * Some useful symmetry obtains.  The d2ii, d2jj and d2ij arrays
       * are symmetric in that their lower and upper triangles are
       * the transposes of one another.  Also, d2jj equals d2ii, and
       * d2ij is the negative of d2ii.
       *
       * As was exploited by David Case and his colleagues, the vector dj
       * is the negative of the vector di.
       */

      di[0] = rx;
      di[1] = ry;
      di[2] = rz;

      dj[0] = -rx;
      dj[1] = -ry;
      dj[2] = -rz;

      /* Load the upper triangle of d2ii. */

      d2ii[0][0] = 1.0 - rx * rx * r2inv;
      d2ii[0][1] = -rx * ry * r2inv;
      d2ii[0][2] = -rx * rz * r2inv;

      d2ii[1][1] = 1.0 - ry * ry * r2inv;
      d2ii[1][2] = -ry * rz * r2inv;

      d2ii[2][2] = 1.0 - rz * rz * r2inv;

      /* Finish loading the rest of all of the matrices. */

      for (k = 0; k < 3; k++) {

         /* Load the upper triangles of d2jj and d2ij. */

         for (l = k; l < 3; l++) {
            d2jj[k][l] = d2ii[k][l];
            d2ij[k][l] = -d2ii[k][l];
         }

         /* Load the symmetric elements of d2ii, d2jj and d2ij. */

         for (l = k + 1; l < 3; l++) {
            d2ii[l][k] = d2ii[k][l];
            d2jj[l][k] = d2jj[k][l];
            d2ij[l][k] = d2ij[k][l];
         }
      }

      /*
       * Calculate the bond stretching energy as:
       *
       *      e = Rk*(r - Req)^2
       */

      db = r - Req[atyp];
      df = Rk[atyp] * db;
      e = df * db;
      e_bond += e;

      /*
       * For the calculation of the 1st derivatives of the energy, the
       * 1st derivatives of the bond distance are multiplied by df:
       *
       *      de/dx = df*(dr/dx)
       *
       * where ds/dx represents the 1st derivative of the bond distance, and: 
       *
       *      df = de/dr = 2*Rk*(r - Req)
       *
       * because the chain rule gives:
       *
       *      de/dx = (de/dr)*(dr/dx) = df*(dr/dx)
       *
       * Divide df by r because this factor is not present in dr/dx.
       */

      df *= 2.0 * rinv;

      dfx = df * rx;
      dfy = df * ry;
      dfz = df * rz;

#ifdef SCALAPACK

      ptr = ptr1d(f, descF_PxQ, at1 + 0);
      if (ptr != NULL)
         *ptr += dfx;
      ptr = ptr1d(f, descF_PxQ, at1 + 1);
      if (ptr != NULL)
         *ptr += dfy;
      ptr = ptr1d(f, descF_PxQ, at1 + 2);
      if (ptr != NULL)
         *ptr += dfz;

      ptr = ptr1d(f, descF_PxQ, at2 + 0);
      if (ptr != NULL)
         *ptr -= dfx;
      ptr = ptr1d(f, descF_PxQ, at2 + 1);
      if (ptr != NULL)
         *ptr -= dfy;
      ptr = ptr1d(f, descF_PxQ, at2 + 2);
      if (ptr != NULL)
         *ptr -= dfz;

#else
      f[at1] += dfx;
      f[at1 + 1] += dfy;
      f[at1 + 2] += dfz;

      f[at2] -= dfx;
      f[at2 + 1] -= dfy;
      f[at2 + 2] -= dfz;
#endif

      /*
       * For the calculation of the 2nd derivatives of the energy, the 1st and 2nd
       * derivatives of the bond distance are multiplied by dg and df as follows:
       *
       *      d2e/dx2dx1 = df*(d2r/dx2dx1) + dg*(dr/dx2)*(dr/dx1)
       *
       * where d2s/dx2dx1 represents a 2nd derivative of the bond distance,
       * dr/dx1 and dr/dx2 represent 1st derivatives of the bond distance, and:
       *
       *      dg = (d/dr)df = (d/dr)(de/dr) = 2*Rk
       *
       * because the chain rule gives:
       *
       *      d2e/dx2dx1 = (d/dx2)(de/dx1)
       *
       *                 = (d/dx2)[(de/dr)*(dr/dx1)]
       *
       *                 = (d/dx2)[df*(dr/dx1)]
       *
       *                 = [(d/dx2)df]*(dr/dx1) + df*[(d/dx2)(dr/dx1)]
       *
       *                 = [(d/dx2)(de/dr)]*(dr/dx1) + df*(d2r/dx2dx1)
       *
       *                 = {[(d/dr)(de/dr)]*(dr/dx2)}*(dr/dx1) + df*(d2r/dx2dx1)
       *
       *                 = dg*(dr/dx2)*(dr/dx1) + df*(d2r/dx2dx1)
       *
       * Divide dg by r^2 because this factor is not present in (dr/dx1)*(dr/dx2).
       */

      dg = 2.0 * Rk[atyp] * r2inv;

      for (k = 0; k < 3; k++) {
         for (l = 0; l < 3; l++) {
            ijterm = df * d2ij[k][l] + dg * di[k] * dj[l];

#ifdef SCALAPACK

            ptr = ptr2d(g, descG_PxQ, at1 + l, at1 + k);
            if (ptr != NULL)
               *ptr += df * d2ii[k][l] + dg * di[k] * di[l];
            ptr = ptr2d(g, descG_PxQ, at2 + l, at2 + k);
            if (ptr != NULL)
               *ptr += df * d2jj[k][l] + dg * dj[k] * dj[l];
            ptr = ptr2d(g, descG_PxQ, at2 + l, at1 + k);
            if (ptr != NULL)
               *ptr += ijterm;
            ptr = ptr2d(g, descG_PxQ, at1 + k, at2 + l);
            if (ptr != NULL)
               *ptr += ijterm;
#else
            g[at1 + k + n * (at1 + l)] +=
                df * d2ii[k][l] + dg * di[k] * di[l];
            g[at2 + k + n * (at2 + l)] +=
                df * d2jj[k][l] + dg * dj[k] * dj[l];
            g[at1 + k + n * (at2 + l)] += ijterm;
            g[at2 + l + n * (at1 + k)] += ijterm;
#endif
         }
      }
   }

   if (e_debug)
      EXPR("%9.3f", e_bond);
   return (e_bond);
}

/***********************************************************************
                            DIFANG()
************************************************************************/

/*
 * This function sets up second derivatives for angle-like potentials,
 * and is used by the eangl2 and ephi2 functions.
 *
 * Calling parameters are as follows:
 *
 * cst - cosine of theta
 * s - bond vectors ij and kj
 * dc - first der. of costheta w/respect to the cartesian differences in s
 * ddc - output of second derivative of costheta w/respect to the cartesian
 *       differences
 */

static
void difang(REAL_T cst, REAL_T * s, REAL_T * dc, REAL_T ** ddc)
{
   REAL_T bij, boi2, boij3, boji3, bkj, boj2, br, bt;
   int i, i1, j;

   /*     ----- first set up distances needed -----   */

   bij = s[1] * s[1] + s[2] * s[2] + s[3] * s[3];
   if ( bij < 0.01 ) bij = 0.01;  /* should cutoff be smaller? */
   boi2 = 1.0 / bij;
   bij = sqrt(bij);
   bkj = s[4] * s[4] + s[5] * s[5] + s[6] * s[6];
   if ( bkj < 0.01 ) bkj = 0.01;  /* should cutoff be smaller? */
   boj2 = 1.0 / bkj;
   bkj = sqrt(bkj);
   br = boi2;

   /*
    *     ----- set up ddc first as second derivative of costheta w/respect
    *           to cartesian differences -----
    */

   for (i = 1; i <= 6; i++) {
      if (i == 4)
         br = boj2;
      ddc[i][i] =
          -br * (2.0 * s[i] * dc[i] + cst * (1.0 - s[i] * s[i] * br));
   }

   for (i = 1; i <= 3; i++) {
      ddc[i][i + 3] = -boi2 * s[i] * dc[i + 3] - boj2 * s[i + 3] * dc[i]
          - boi2 * boj2 * (s[i] * s[i + 3] * cst - bij * bkj);
   }

   for (i = 1; i <= 2; i++) {
      i1 = i + 1;
      for (j = i1; j <= 3; j++) {
         ddc[i][j] =
             -boi2 * (s[i] * dc[j] + s[j] * dc[i] -
                      boi2 * cst * s[i] * s[j]);
      }
   }

   for (i = 4; i <= 5; i++) {
      i1 = i + 1;
      for (j = i1; j <= 6; j++) {
         ddc[i][j] =
             -boj2 * (s[i] * dc[j] + s[j] * dc[i] -
                      boj2 * cst * s[i] * s[j]);
      }
   }

   boij3 = boj2 / (bij * bkj);
   boji3 = boi2 / (bij * bkj);
   bt = boi2 * boj2 * cst;
   for (i = 1; i <= 3; i++) {
      for (j = 4; j <= 6; j++) {
         if (j != i + 3)
            ddc[i][j] =
                -boij3 * s[i + 3] * s[j] - boji3 * s[i] * s[j - 3] +
                bt * s[i] * s[j];
      }
   }

   /*    ----- symmetrize the ddc array -----   */

   for (i = 1; i <= 6; i++) {
      for (j = i; j <= 6; j++) {
         ddc[j][i] = ddc[i][j];
      }
   }
   return;
}

/***********************************************************************
                            EANGL2()
************************************************************************/

/* Calculate the bond bending energy and derivatives.*/

static
REAL_T eangl2(int nang, int *a1, int *a2, int *a3, int *atype,
              REAL_T * Tk, REAL_T * Teq, REAL_T * x, REAL_T * f,
              REAL_T * g, int context_PxQ, int *descF_PxQ, int *descG_PxQ)
{
   int i, j, ia, atyp, at1, at2, at3, ist, iof, inew, jst, jof, jnew;
   size_t n;
   REAL_T x1, y1, z1, x2, y2, z2, x3, y3, z3;
   REAL_T c, st, st2c, dtheta, theta, df, ddf, e, e_theta;
   REAL_T rij, rkj, rij2, rkj2, rrik;

   REAL_T s[7], dc[7], dr[10];
   REAL_T **ddc, **ddr;

#ifdef SCALAPACK
   REAL_T *ptr;
#endif

   ddc = matrix(1, 6, 1, 6);
   ddr = matrix(1, 9, 1, 9);

   n = 3 * prm->Natom;
   e_theta = 0.0;
   for (ia = 0; ia < nang; ia++) {
      at1 = a1[ia];
      at2 = a2[ia];
      at3 = a3[ia];
      atyp = atype[ia] - 1;

      /* Set up the variables x1..x3, y1..y3 and z1..z3 for use later. */

      x1 = x[at1];
      y1 = x[at1 + 1];
      z1 = x[at1 + 2];

      x2 = x[at2];
      y2 = x[at2 + 1];
      z2 = x[at2 + 2];

      x3 = x[at3];
      y3 = x[at3 + 1];
      z3 = x[at3 + 2];

      s[1] = x1 - x2;
      s[2] = y1 - y2;
      s[3] = z1 - z2;
      s[4] = x3 - x2;
      s[5] = y3 - y2;
      s[6] = z3 - z2;
      rij2 = s[1] * s[1] + s[2] * s[2] + s[3] * s[3];
      rkj2 = s[4] * s[4] + s[5] * s[5] + s[6] * s[6];
      rij = sqrt(rij2);
      rkj = sqrt(rkj2);
      rrik = rij * rkj;
      c = (s[1] * s[4] + s[2] * s[5] + s[3] * s[6]) / rrik;

      c = c > 1.0 ? 1.0 : c;
      c = c < -1.0 ? -1.0 : c;
      theta = acos(c);
      dtheta = theta - Teq[atyp];

      /* df and ddf are derivatives of E with respect to c == costheta: */
      df = dtheta * Tk[atyp];
      ddf = 2. * Tk[atyp];
      e = df * dtheta;
      e_theta += e;

      /*
       * Calculate the sine then limit small values of the sine to the range
       * -10**-3..+10**-3 so as to avoid division by zero.
       */

      st = sin(theta);
      if (st > 0 && st < 1.e-3)
         st = 1.e-3;
      else if (st < 0 && st > -1.e-3)
         st = -1.e-3;
      df *= 2.0;                /* check this!!  */

      /* dc = derivative of c with respect to cartesian differences: */

      for (i = 1; i <= 3; i++) {
         dc[i] = (s[i + 3] / rkj - c * s[i] / rij) / rij;
         dc[i + 3] = (s[i] / rij - c * s[i + 3] / rkj) / rkj;
      }

      /* get ddc = second derivative of c with respect to 
         cartesian differences:  */
      difang(c, s, dc, ddc);

      /* change ddc to second. der. of theta with respect to 
         cartesian differences:  */
      st2c = c / (st * st);
      for (i = 1; i <= 6; i++) {
         for (j = i; j <= 6; j++) {
            ddc[i][j] = -(ddc[i][j] + dc[i] * dc[j] * st2c) / st;
            ddc[j][i] = ddc[i][j];
         }
      }

      /* change dc to derivative of theta w/ respect to cartesian
         differences:  */
      for (i = 1; i <= 6; i++) {
         dc[i] = -dc[i] / st;
      }

      /* dr will hold -derivates of theta w/ respect to cartesians: */
      for (i = 1; i <= 3; i++) {
         dr[i] = -dc[i];
         dr[i + 6] = -dc[i + 3];
         dr[i + 3] = dc[i] + dc[i + 3];
      }

      /* update the forces:  */

#ifdef SCALAPACK

      ptr = ptr1d(f, descF_PxQ, at1 + 0);
      if (ptr != NULL)
         *ptr -= df * dr[1];
      ptr = ptr1d(f, descF_PxQ, at1 + 1);
      if (ptr != NULL)
         *ptr -= df * dr[2];
      ptr = ptr1d(f, descF_PxQ, at1 + 2);
      if (ptr != NULL)
         *ptr -= df * dr[3];

      ptr = ptr1d(f, descF_PxQ, at2 + 0);
      if (ptr != NULL)
         *ptr -= df * dr[4];
      ptr = ptr1d(f, descF_PxQ, at2 + 1);
      if (ptr != NULL)
         *ptr -= df * dr[5];
      ptr = ptr1d(f, descF_PxQ, at2 + 2);
      if (ptr != NULL)
         *ptr -= df * dr[6];

      ptr = ptr1d(f, descF_PxQ, at3 + 0);
      if (ptr != NULL)
         *ptr -= df * dr[7];
      ptr = ptr1d(f, descF_PxQ, at3 + 1);
      if (ptr != NULL)
         *ptr -= df * dr[8];
      ptr = ptr1d(f, descF_PxQ, at3 + 2);
      if (ptr != NULL)
         *ptr -= df * dr[9];

#else

      f[at1] -= df * dr[1];
      f[at1 + 1] -= df * dr[2];
      f[at1 + 2] -= df * dr[3];

      f[at2] -= df * dr[4];
      f[at2 + 1] -= df * dr[5];
      f[at2 + 2] -= df * dr[6];

      f[at3] -= df * dr[7];
      f[at3 + 1] -= df * dr[8];
      f[at3 + 2] -= df * dr[9];

#endif

      /* ddr will hold second derivatives of theta with respect to carts: */
      for (i = 1; i <= 3; i++) {
         for (j = i; j <= 3; j++) {
            ddr[i][j] = ddc[i][j];
            ddr[i + 6][j + 6] = ddc[i + 3][j + 3];
         }
         for (j = 1; j <= 3; j++) {
            ddr[i][j + 6] = ddc[i][j + 3];
            ddr[i][j + 3] = -ddc[i][j] - ddc[i][j + 3];
            ddr[i + 3][j + 6] = -ddc[i + 3][j] - ddc[i + 3][j + 3];
         }
      }
      for (i = 1; i <= 3; i++) {
         for (j = i; j <= 3; j++) {
            ddr[i + 3][j + 3] = -ddr[i][j + 3] - ddr[i + 3][j + 6];
         }
      }
      /*
         for( i=1; i<=9; i++ ){
         for( j=i; j<=9; j++ ){
         ddr[j][i] = ddr[i][j];
         }
         }
       */

      for (i = 1; i <= 9; i++) {
         ist = at3;
         if (i <= 6)
            ist = at2;
         if (i <= 3)
            ist = at1;
         iof = i % 3;
         if (iof == 0)
            iof = 3;
         inew = ist + iof - 1;

         for (j = i; j <= 9; j++) {
            jst = at3;
            if (j <= 6)
               jst = at2;
            if (j <= 3)
               jst = at1;
            jof = j % 3;
            if (jof == 0)
               jof = 3;
            jnew = jst + jof - 1;

#ifdef SCALAPACK

            ptr = ptr2d(g, descG_PxQ, jnew, inew);
            if (ptr != NULL)
               *ptr += ddf * dr[i] * dr[j] + df * ddr[i][j];
            if (inew != jnew) {
               ptr = ptr2d(g, descG_PxQ, inew, jnew);
               if (ptr != NULL)
                  *ptr += ddf * dr[i] * dr[j] + df * ddr[i][j];
            }
#else
            g[inew + n * jnew] += ddf * dr[i] * dr[j] + df * ddr[i][j];
            if (inew != jnew)
               g[jnew + n * inew] += ddf * dr[i] * dr[j] + df * ddr[i][j];
#endif
         }
      }
   }

   free_matrix(ddc, 1, 6, 1, 6);
   free_matrix(ddr, 1, 6, 1, 6);

   if (e_debug)
      EXPR("%9.3f", e_theta);
   return (e_theta);
}

/***********************************************************************
                            EPHI2()
************************************************************************/

/* Calculate the dihedral torsion energy and derivatives. */

static
REAL_T ephi2(int nphi, int *a1, int *a2,
             int *a3, int *a4, int *atype,
             REAL_T * Pk, REAL_T * Pn, REAL_T * Phase,
             REAL_T * x, REAL_T * f, REAL_T * g,
             int context_PxQ, int *descF_PxQ, int *descG_PxQ)
{
   int i, j, k, l, at1, at2, at3, at4, atyp, iper, iphi;
   size_t n;
   int inew, ist, iof, jnew, jst, jof;
   REAL_T e, df, ddf, e_tors;
   REAL_T xij, yij, zij, xkj, ykj, zkj, xkl, ykl, zkl;
   REAL_T dx, dy, dz, gx, gy, gz, bi, bk, ct, z1, z2;
   REAL_T ux, uy, uz, vx, vy, vz, dx1, dy1, dz1, delta;
   REAL_T phi, yy;
   REAL_T ct2, arg, z11, z12, z22;
   REAL_T t[7], dc[7], dr[13];
   REAL_T dum, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12;
   REAL_T **dtx, **ddr, **ddc;

#ifdef SCALAPACK
   REAL_T *ptr;
#endif

   ddc = matrix(1, 6, 1, 6);
   dtx = matrix(1, 6, 1, 12);
   ddr = matrix(1, 12, 1, 12);

   n = 3 * prm->Natom;
   e_tors = 0.0;
   for (iphi = 0; iphi < nphi; iphi++) {

      /*
       * Get the atom numbers.  Atom 3 may be negative to specify that
       * the end group interactions are to be ignored.  Atom 4 may be 
       * negative to specify that the dihedral angle is an "improper
       * torsion".  However, neither of these conditions is actually
       * treated in this function, so just ignore negative values.
       */

      at1 = a1[iphi];
      at2 = a2[iphi];
      at3 = abs(a3[iphi]);
      at4 = abs(a4[iphi]);
      atyp = atype[iphi] - 1;

      /*     ----- calculation of the dihedral -----     */

      xij = x[at1] - x[at2];
      yij = x[at1 + 1] - x[at2 + 1];
      zij = x[at1 + 2] - x[at2 + 2];

      xkj = x[at3] - x[at2];
      ykj = x[at3 + 1] - x[at2 + 1];
      zkj = x[at3 + 2] - x[at2 + 2];

      xkl = x[at3] - x[at4];
      ykl = x[at3 + 1] - x[at4 + 1];
      zkl = x[at3 + 2] - x[at4 + 2];

      dx = yij * zkj - zij * ykj;
      dy = zij * xkj - xij * zkj;
      dz = xij * ykj - yij * xkj;

      gx = zkj * ykl - ykj * zkl;
      gy = xkj * zkl - zkj * xkl;
      gz = ykj * xkl - xkj * ykl;

      bi = dx * dx + dy * dy + dz * dz;
      bk = gx * gx + gy * gy + gz * gz;
      ct = dx * gx + dy * gy + dz * gz;

      /*    ----- approximate if linear dihedral 
            ----- assumes force constant is eventually zero ----- */

      if ( bi <= 0.01 ) bi = 0.01;
      if ( bk <= 0.01 ) bk = 0.01;

      bi = sqrt(bi);
      bk = sqrt(bk);
      z1 = 1. / bi;
      z2 = 1. / bk;
      ct = ct * z1 * z2;
      ct = ct > 1.0 ? 1.0 : ct;
      ct = ct < -1.0 ? -1.0 : ct;

      /*     ----- ct value above is actually -cosphi; here we change
       *           its sign -----
       */
      ct = -ct;

      /*     ----- calculate the energy and derivatives -----


       *     ----- get df = first der. of potential w/respect to cosphi; and
       *           ddf = second der. of potntial w/respect to cosphi -----

       *           the torsional potential is assumed to have the form:
       *            e = pk(ic) * (1.0+phase*cos(pn(ic)*phi)
       *            where phase = 1.0 or -1.0, and pn = 1,2,3,4, or 6

       *     ----- energy terms for dihedrals are expressed in terms of
       *           cosphi in order to eliminate problems for planar angles ----
       */

    multi_term:
         if ((fabs(Phase[atyp] - 3.142) > 0.01) && (fabs(Phase[atyp]) > 0.01)){

            /*   here we have special code for phases that are not zero or pi */

            phi = acos(ct);

            /*
               now calculate sin(phi) because cos(phi) is symmetric, so
               we can decide between +-phi.
             */

            ux = -yij * zkj + zij * ykj;
            uy = -zij * xkj + xij * zkj;
            uz = -xij * ykj + yij * xkj;

            vx = -ykj * zkl + zkj * ykl;
            vy = -zkj * xkl + xkj * zkl;
            vz = -xkj * ykl + ykj * xkl;

            dx1 = uy * vz - uz * vy;
            dy1 = uz * vx - ux * vz;
            dz1 = ux * vy - uy * vx;

            dx1 = dx1*xkj + dy1*ykj + dz1*zkj;
            if (dx1 < 0.0)
               phi = -phi;

            delta = Pn[atyp]*phi - Phase[atyp];
            e = Pk[atyp]*(1.0 + cos(delta));
            e_tors += e;
            yy = sin(phi);

            if (fabs(yy) > 0.001) {
               df = Pk[atyp]*Pn[atyp]*sin(delta)/yy;
               ddf = -Pk[atyp]*Pn[atyp]*( Pn[atyp]*cos(delta) -
                       ct*sin(delta)/yy )/(yy*yy);
            } else {
               /* if sin(phi) happens to be at zero or pi, adjust it slightly */
               if (phi > -1. && phi < 1.) {  /* set sin(phi) = 0.001 */
                  df = Pk[atyp]*Pn[atyp]*sin(delta) * 1000.;
                  ddf = -Pk[atyp]*Pn[atyp]*1000000. * ( Pn[atyp]*cos(delta) -
                          ct*1000.*sin(delta) );
               } else {  /* set sin(phi) = -0.001  */
                  df = -Pk[atyp]*Pn[atyp]*sin(delta) * 1000.;
                  ddf = -Pk[atyp]*Pn[atyp]*1000000. * ( Pn[atyp]*cos(delta) +
                          ct*1000.*sin(delta) );
               }
            }

         } else {   /* here is "usual" case, where the phase is 0 or pi */

      iper = (int) fabs(Pn[atyp]) + 0.0001;
      assert(iper != (fabs(Pn[atyp]) - 0.0001));
      assert(iper >= 1 && iper <= 6);
      switch (iper) {

      case 1:
         e = ct;
         df = 1.0;
         ddf = 0.0;
         break;

      case 2:
         e = 2.0*ct*ct - 1.0;
         df = 4.0*ct;
         ddf = 4.0;
         break;

      case 3:
         ct2 = ct*ct;
         e = ct*(4.0*ct2 - 3.0);
         df = 12.0*ct2 - 3.0;
         ddf = 24.0*ct;
         break;

      case 4:
         ct2 = ct*ct;
         e = 1.0 + ct2*8.0*(ct2 - 1.0);
         df = 32.0*ct2*ct - 16.*ct;
         ddf = 96.*ct2 -16.;
         break;

      case 5:
         ct2 = ct*ct;
         e = 16.*ct2*ct2*ct - 20.*ct2*ct + 5.*ct;
         df = 80.*ct2*ct2 - 60.*ct2 + 5.;
         ddf = 320.*ct2*ct - 120.*ct;
         break;

      case 6:
         ct2 = ct * ct;
         e = ct2 * (ct2 * (ct2 * 32.0 - 48.0) + 18.0) - 1.0;
         df = ct * (ct2 * (ct2 * 192.0 - 192.0) + 36.0);
         ddf = ct2 * (ct2 * 960.0 - 576.0) + 36.0;
         break;

      default:
         fprintf(stderr, "bad periodicity: %d\n", iper);
         exit(1);
      }

      if (fabs(Phase[atyp] - 3.142) < 0.01)
         arg = -Pk[atyp];
      else
         arg = Pk[atyp];

      e = Pk[atyp] + arg * e;
      df = df * arg;
      ddf = ddf * arg;
      e_tors += e;
      }

      t[1] = dx;
      t[2] = dy;
      t[3] = dz;

      t[4] = -gx;
      t[5] = -gy;
      t[6] = -gz;

      /*
       *     ----- now, set up array dc = first der. of cosphi w/respect
       *           to the cartesian differences t -----
       */

      z11 = z1 * z1;
      z12 = z1 * z2;
      z22 = z2 * z2;
      for (i = 1; i <= 3; i++) {
         dc[i] = t[i + 3] * z12 - ct * t[i] * z11;
         dc[i + 3] = t[i] * z12 - ct * t[i + 3] * z22;
      }

      /*
       *     ----- subroutine difang will now create array ddc which is second
       *           derivative of cosphi with respect to the t s -----
       */

      difang(ct, t, dc, ddc);

      /* ---now set up array s, given on page 118 of cff book-- */

      s1 = xij;
      s2 = yij;
      s3 = zij;
      s4 = xkj;
      s5 = ykj;
      s6 = zkj;
      s7 = -xkj;
      s8 = -ykj;
      s9 = -zkj;
      s10 = -xkl;
      s11 = -ykl;
      s12 = -zkl;

      /*
       *     ----- set up dtx[i][j] = derivative of t[i] w/respect to x[j]
       *           see p. 120 of cff book -----
       */

      for (i = 1; i <= 6; i++) {
         for (j = 1; j <= 12; j++) {
            dtx[i][j] = 0.0;
         }
      }

      dtx[1][2] = s6;
      dtx[1][3] = -s5;
      dtx[1][5] = s3 - s6;
      dtx[1][6] = s5 - s2;
      dtx[1][8] = -s3;
      dtx[1][9] = s2;
      dtx[2][1] = -s6;
      dtx[2][3] = s4;
      dtx[2][4] = s6 - s3;
      dtx[2][6] = s1 - s4;
      dtx[2][7] = s3;
      dtx[2][9] = -s1;
      dtx[3][1] = s5;
      dtx[3][2] = -s4;
      dtx[3][4] = s2 - s5;
      dtx[3][5] = s4 - s1;
      dtx[3][7] = -s2;
      dtx[3][8] = s1;
      dtx[4][5] = s12;
      dtx[4][6] = -s11;
      dtx[4][8] = s9 - s12;
      dtx[4][9] = s11 - s8;
      dtx[4][11] = -s9;
      dtx[4][12] = s8;
      dtx[5][4] = -s12;
      dtx[5][6] = s10;
      dtx[5][7] = s12 - s9;
      dtx[5][9] = s7 - s10;
      dtx[5][10] = s9;
      dtx[5][12] = -s7;
      dtx[6][4] = s11;
      dtx[6][5] = -s10;
      dtx[6][7] = s8 - s11;
      dtx[6][8] = s10 - s7;
      dtx[6][10] = -s8;
      dtx[6][11] = s7;

      /*
       *     ----- set up dr array, containing -first derivative of cosphi with
       *           respect to cartesians -----
       */

      for (i = 1; i <= 12; i++) {
         dum = 0.0;
         for (j = 1; j <= 6; j++) {
            dum += dc[j] * dtx[j][i];
         }
         dr[i] = -dum;
      }


      /*     ----- update the force array ----- */

#ifdef SCALAPACK

      ptr = ptr1d(f, descF_PxQ, at1 + 0);
      if (ptr != NULL)
         *ptr -= df * dr[1];
      ptr = ptr1d(f, descF_PxQ, at1 + 1);
      if (ptr != NULL)
         *ptr -= df * dr[2];
      ptr = ptr1d(f, descF_PxQ, at1 + 2);
      if (ptr != NULL)
         *ptr -= df * dr[3];

      ptr = ptr1d(f, descF_PxQ, at2 + 0);
      if (ptr != NULL)
         *ptr -= df * dr[4];
      ptr = ptr1d(f, descF_PxQ, at2 + 1);
      if (ptr != NULL)
         *ptr -= df * dr[5];
      ptr = ptr1d(f, descF_PxQ, at2 + 2);
      if (ptr != NULL)
         *ptr -= df * dr[6];

      ptr = ptr1d(f, descF_PxQ, at3 + 0);
      if (ptr != NULL)
         *ptr -= df * dr[7];
      ptr = ptr1d(f, descF_PxQ, at3 + 1);
      if (ptr != NULL)
         *ptr -= df * dr[8];
      ptr = ptr1d(f, descF_PxQ, at3 + 2);
      if (ptr != NULL)
         *ptr -= df * dr[9];

      ptr = ptr1d(f, descF_PxQ, at4 + 0);
      if (ptr != NULL)
         *ptr -= df * dr[10];
      ptr = ptr1d(f, descF_PxQ, at4 + 1);
      if (ptr != NULL)
         *ptr -= df * dr[11];
      ptr = ptr1d(f, descF_PxQ, at4 + 2);
      if (ptr != NULL)
         *ptr -= df * dr[12];

#else

      f[at1] -= df * dr[1];
      f[at1 + 1] -= df * dr[2];
      f[at1 + 2] -= df * dr[3];

      f[at2] -= df * dr[4];
      f[at2 + 1] -= df * dr[5];
      f[at2 + 2] -= df * dr[6];

      f[at3] -= df * dr[7];
      f[at3 + 1] -= df * dr[8];
      f[at3 + 2] -= df * dr[9];

      f[at4] -= df * dr[10];
      f[at4 + 1] -= df * dr[11];
      f[at4 + 2] -= df * dr[12];

#endif

      /*
       *     ----- now set up the ddr array = second der. of cosphi w/respect
       *           to cartesians; first we take the first term of last formula
       *           on p. 113 of cff book -----
       */

      for (i = 1; i <= 12; i++) {
         for (j = i; j <= 12; j++) {
            ddr[i][j] = 0.0;
            for (k = 1; k <= 6; k++) {
               for (l = 1; l <= 6; l++) {
                  ddr[i][j] =
                      ddr[i][j] + ddc[k][l] * dtx[k][i] * dtx[l][j];
               }
            }
         }
      }

      /*     ----- now do the second term of this equation -----  */

      ddr[2][9] += dc[1];
      ddr[3][5] += dc[1];
      ddr[6][8] += dc[1];
      ddr[2][6] -= dc[1];
      ddr[3][8] -= dc[1];
      ddr[5][9] -= dc[1];
      ddr[1][6] += dc[2];
      ddr[3][7] += dc[2];
      ddr[4][9] += dc[2];
      ddr[1][9] -= dc[2];
      ddr[3][4] -= dc[2];
      ddr[6][7] -= dc[2];
      ddr[1][8] += dc[3];
      ddr[2][4] += dc[3];
      ddr[5][7] += dc[3];
      ddr[1][5] -= dc[3];
      ddr[2][7] -= dc[3];
      ddr[4][8] -= dc[3];
      ddr[5][12] += dc[4];
      ddr[6][8] += dc[4];
      ddr[9][11] += dc[4];
      ddr[5][9] -= dc[4];
      ddr[6][11] -= dc[4];
      ddr[8][12] -= dc[4];
      ddr[4][9] += dc[5];
      ddr[6][10] += dc[5];
      ddr[7][12] += dc[5];
      ddr[4][12] -= dc[5];
      ddr[6][7] -= dc[5];
      ddr[9][10] -= dc[5];
      ddr[4][11] += dc[6];
      ddr[5][7] += dc[6];
      ddr[8][10] += dc[6];
      ddr[4][8] -= dc[6];
      ddr[5][10] -= dc[6];
      ddr[7][11] -= dc[6];

      /*     ----- Now form the second derivative matrix -----  */

      for (i = 1; i <= 12; i++) {
         ist = at4;
         if (i <= 9)
            ist = at3;
         if (i <= 6)
            ist = at2;
         if (i <= 3)
            ist = at1;
         iof = i % 3;
         if (iof == 0)
            iof = 3;
         inew = ist + iof - 1;
         for (j = i; j <= 12; j++) {
            jst = at4;
            if (j <= 9)
               jst = at3;
            if (j <= 6)
               jst = at2;
            if (j <= 3)
               jst = at1;
            jof = j % 3;
            if (jof == 0)
               jof = 3;
            jnew = jst + jof - 1;

#ifdef SCALAPACK

            ptr = ptr2d(g, descG_PxQ, jnew, inew);
            if (ptr != NULL)
               *ptr += ddf * dr[i] * dr[j] + df * ddr[i][j];
            if (inew != jnew) {
               ptr = ptr2d(g, descG_PxQ, inew, jnew);
               if (ptr != NULL)
                  *ptr += ddf * dr[i] * dr[j] + df * ddr[i][j];
            }
#else
            g[inew + n * jnew] += ddf * dr[i] * dr[j] + df * ddr[i][j];
            if (inew != jnew)
               g[jnew + n * inew] += ddf * dr[i] * dr[j] + df * ddr[i][j];
#endif
         }
      }


#ifdef PRINT_EPHI
      printf("%4d %4d %4d %4d %4d %9.4f\n", i + 1, at1 / 3, at2 / 3,
             at3 / 3, at4 / 3, e);
#endif

      /* A negative value of Pn means that more terms are needed in the expansion. */

      if (Pn[atyp] < 0) {
         atyp++;
         goto multi_term;
      }
   }
   if (e_debug)
      EXPR("%9.3f", e_tors);

   free_matrix(dtx, 1, 6, 1, 12);
   free_matrix(ddr, 1, 12, 1, 12);
   free_matrix(ddc, 1, 6, 1, 6);

   return (e_tors);
}

/***********************************************************************
                            NBOND2()
************************************************************************/

/* 
 * Calculate the non-bonded energy, 1st and 2nd derivatives.
 * This function is complicated by the fact that it must
 * process two forms of pair lists: the 1-4 pair list and
 * the non-bonded pair list.  The non-bonded pair list
 * must be modified by the excluded atom list whereas the
 * 1-4 pair list is used unmodified.  Also, the non-bonded
 * pair list comprises lower and upper triangles whereas
 * the 1-4 pair list comprises an upper triangle only.
 *
 * Calling parameters are as follows:
 *
 * lpears - the number of pairs on the lower triangle pair list
 * upears - the number of pairs on the upper trianble pair list
 * pearlist - either the 1-4 pair list or the non-bonded pair list
 * N14 - set to 0 for the non-bonded pair list, 1 for the 1-4 pair list
 * x - the atomic coordinate array
 * f - the gradient vector
 * g - the Hessian matrix mapped onto a linear vector
 * enb - Van der Waals energy return value, passed by reference
 * eel - Coulombic energy return value, passed by reference
 * scnblist - list of 1 / scale factors for vdw energy (1-4 interactions)
 * sceelist - list of 1 / scale factors for Coulombic energy (1-4 interactions)
 * context_PxQ - input: the distributed vector and matrix context for ScaLAPACK
 * descF_PxQ - input: the ScaLAPACK descriptor for vector f in context_PxQ
 * descG_PxQ - input: the ScaLAPACK descriptor for matrix g in context_PxQ
 *
 * Note, scnblist and sceelist are NULL when N14 == 0, so assume scale factors
 * are 1 in that case.
 */

static
int nbond2(int *lpears, int *upears, int **pearlist, int N14,
           REAL_T * x, REAL_T * f, REAL_T * g, REAL_T * enb, REAL_T * eel,
           REAL_T **scnblist, REAL_T **sceelist,
           int context_PxQ, int *descF_PxQ, int *descG_PxQ)
{
   int i, j, k, l, jn, ic, npr, lpair, iaci;
   size_t n;
   int *iexw;
   REAL_T r, r2, r2inv, df2, dg2, r6, r10, f1, f2, df, dg;
   REAL_T dis, kij, d0, diff, rinv, rs, rssq, eps1, epsi, cgijr, pow, cgi, cgj;
   REAL_T dumx, dumy, dumz, enbfaci, eelfaci;
   int nhbpair, ibig, isml;
   REAL_T xi, yi, zi, xij, yij, zij;
   REAL_T di[3], dj[3], d2ii[3][3], d2jj[3][3], d2ij[3][3], ijterm;

#ifdef SCALAPACK
   REAL_T *ptr;
#endif

#if 0
   REAL_T hbener, nb14;
   nhbpair = 0;
   hbener = 0.0;
   nb14 = 0.0;
#endif

#define SIG 0.3
#define DIW 78.0
#define C1 38.5

   /*
    * For the various methods of calculating the nonbonded energy
    * both the energy and its 1st and 2nd derivatives are calculated.
    *
    * For the calculation of the 1st derivatives of the energy, the
    * 1st derivatives the nonbonded distance r are multiplied by df:
    *
    *      de/dx = df*(dr/dx)
    *
    * where dr/dx represents the 1st derivative of the distance, and 
    * 
    *      df = de/dr
    *
    * is required by the chain rule:
    *
    *      de/dx = (de/dr)*(dr/dx) = df*(dr/dx)
    *
    * For the calculation of the 2nd derivatives of the energy, the 1st and
    * 2nd derivatives of the distance are multiplied by dg and df as follows:
    *
    *      d2e/dxjdxi = df*(d2r/dxjdxi) + dg*(dr/dxj)*(dr/dxi)
    *
    * where d2r/dxjdxi represents a 2nd derivative of the distance,
    * dr/dxi and dr/dxj represent 1st derivatives of the distance, and:
    *
    *      dg = (d/dr)df
    *
    * is required by the chain rule:
    *
    *      d2e/dxjdxi = (d/dxj)(de/dxi)
    *
    *                 = (d/dxj)[(de/dr)*(dr/dxi)]
    *
    *                 = (d/dxj)[df*(dr/dxi)]
    *
    *                 = [(d/dxj)df]*(dr/dxi) + df*[(d/dxj)(dr/dxi)]
    *
    *                 = [(d/dxj)(de/dr)]*(dr/dxi) + df*(d2r/dxjdxi)
    *
    *                 = {[(d/dr)(de/dr)]*(dr/dxj)}*(dr/dxi) + df*(d2r/dxjdxi)
    *
    *                 = dg*(dr/dxj)*(dr/dxi) + df*(d2r/dxjdxi)
    *
    */

   n = 3 * prm->Natom;
   *enb = 0.;
   *eel = 0.;
   // Default scale factors are 1.0 and will be re-set for 1-4 interactions from
   // sceelist and scnblist
   enbfaci = 1.0;
   eelfaci = 1.0;

   /*
    * Allocate and initialize the iexw array used for skipping excluded
    * atoms.  Note that because of the manner in which iexw is used, it
    * is necessary to initialize it before only the first iteration of
    * the following loop.
    */

   iexw = ivector(-1, prm->Natom);
   for (i = -1; i < prm->Natom; i++) {
      iexw[i] = -1;
   }

   /* Loop over all atoms i except for the final atom i. */

   for (i = 0; i < prm->Natom - 1; i++) {

      /* Check whether there are any atoms j on the pair list of atom i. */

      npr = upears[i];
      if (npr <= 0)
         continue;

      iaci = prm->Ntypes * (prm->Iac[i] - 1);
      dumx = 0.0;
      dumy = 0.0;
      dumz = 0.0;
      xi = x[3 * i];
      yi = x[3 * i + 1];
      zi = x[3 * i + 2];
      cgi = prm->Charges[i];

      /*
       * Expand the excluded list into the iexw array by storing i
       * at array address j.
       */

      for (j = 0; j < prm->Iblo[i]; j++) {
         iexw[IexclAt[i][j] - 1] = i;
      }

      /*
       * If the 'N14' calling parameter is clear, use the beginning
       * address of the upper triangle pair list, which happens
       * to be the number of atoms on the lower triangle pair list.
       * If the 'N14' calling parameter is set, the beginning
       * address is zero because no lower triangle pair list is
       * used for the N14 interactions.
       */

      if (N14 == 0) {
         lpair = lpears[i];
      } else {
         lpair = 0;
      }

      /* Select atom j from the pair list.  Non-graceful error handling. */

      for (jn = 0; jn < npr; jn++) {

         if (pearlist[i] == NULL) {
            printf("NULL pair list entry in nbond2, taskid = %d\n",
                   mytaskid);
            fflush(stdout);
         }
         j = pearlist[i][lpair + jn];

         /* If the 'N14' calling parameter is set, grab eelfaci and enbfaci from
          * the sceelist and scnblist, respectively.  The proper factors are set
          * up in exactly the same way as the 1-4 pairlist so indices can be
          * re-used (and they are already inverted)
          */

         if (N14) {
            enbfaci = scnblist[i][lpair + jn];
            eelfaci = sceelist[i][lpair + jn];
         }

         cgj = prm->Charges[j] * eelfaci;

         /*
          * If the 'N14' calling parameter is clear, check whether
          * this i,j pair is exempted by the excluded atom list.
          */

         if (N14 != 0 || iexw[j] != i) {
            xij = xi - x[3 * j];
            yij = yi - x[3 * j + 1];
            zij = zi - x[3 * j + 2];
            r2 = xij * xij + yij * yij + zij * zij;
            r2inv = 1.0 / r2;
            r = sqrt(r2);
            rinv = r * r2inv;

            /*
             * Calculate the first and second derivatives of the interatomic
             * distance Dij with respect to the cartesian coordinates of atoms
             * i and j.  The results are placed into five arrays:
             *
             *   di[] for the first derivatives with respect to atom i
             *   dj[] for the first derivatives with respect to atom j
             *   d2ii[] for the second derivatives with respect to atom i
             *   d2jj[] for the second derivatives with respect to atom j
             *   d2ij[] for the second derivatives with respect to atoms i and j
             *
             * Some useful symmetry obtains.  The d2ii, d2jj and d2ij arrays
             * are symmetric in that their lower and upper triangles are
             * the transposes of one another.  Also, d2jj equals d2ii, and
             * d2ij is the negative of d2ii.
             *
             * As was exploited by David Case and his colleagues, the vector dj
             * is the negative of the vector di.
             */

            di[0] = xij;
            di[1] = yij;
            di[2] = zij;

            dj[0] = -xij;
            dj[1] = -yij;
            dj[2] = -zij;

            /* Load the upper triangle of d2ii. */

            d2ii[0][0] = 1.0 - xij * xij * r2inv;
            d2ii[0][1] = -xij * yij * r2inv;
            d2ii[0][2] = -xij * zij * r2inv;

            d2ii[1][1] = 1.0 - yij * yij * r2inv;
            d2ii[1][2] = -yij * zij * r2inv;

            d2ii[2][2] = 1.0 - zij * zij * r2inv;

            /* Finish loading the rest of all of the matrices. */

            for (k = 0; k < 3; k++) {

               /* Load the upper triangles of d2jj and d2ij. */

               for (l = k; l < 3; l++) {
                  d2jj[k][l] = d2ii[k][l];
                  d2ij[k][l] = -d2ii[k][l];
               }

               /* Load the symmetric elements of d2ii, d2jj and d2ij. */

               for (l = k + 1; l < 3; l++) {
                  d2ii[l][k] = d2ii[k][l];
                  d2jj[l][k] = d2jj[k][l];
                  d2ij[l][k] = d2ij[k][l];
               }
            }

            /* Calculate the energy and derivatives according to dield. */

            if (dield == -3) {

               /* special code Ramstein & Lavery dielectric, 94 force field */

               rs = SIG * r;
               rssq = rs * rs;
               pow = exp(-rs);
               eps1 = rssq + rs + rs + 2.0;
               epsi = 1.0 / (DIW - C1 * pow * eps1);
               cgijr = cgi * cgj * rinv * epsi;
               *eel += cgijr;
               df2 = -cgijr * (1.0 + C1 * pow * rs * rssq * epsi);
               dg2 =
                   cgijr * (epsi *
                            (C1 * pow * SIG *
                             (rssq *
                              (2.0 * (epsi * C1 * pow * SIG * rssq + rinv)
                               + SIG) - 2.0 * SIG * rs)) + 2.0 * r2inv);
               ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
               if (ic >= 0) {
                  r6 = r2inv * r2inv * r2inv;
                  f2 = prm->Cn2[ic] * r6;
                  f1 = prm->Cn1[ic] * r6 * r6;
                  *enb += (f1 - f2) * enbfaci;
                  df = (df2 + (6.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
                  dg = dg2 + (156.0 * f1 - 42.0 * f2 * enbfaci) * r2inv;
               } else {
                  df = df2 * rinv;
                  dg = dg2;
               }

            } else if (dield == -4) {

               /* distance-dependent dielectric code, 94 ff */
               /* epsilon = r  */

               rs = cgi * cgj * r2inv;
               df2 = -2.0 * rs;
               dg2 = 6.0 * rs * r2inv;
               *eel += rs;
               ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
               if (ic >= 0) {
                  r6 = r2inv * r2inv * r2inv;
                  f2 = prm->Cn2[ic] * r6;
                  f1 = prm->Cn1[ic] * r6 * r6;
                  *enb += (f1 - f2) * enbfaci;
                  df = (df2 + (6.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
                  dg = dg2 + (156.0 * f1 - 42.0 * f2 * enbfaci) * r2inv;
               } else {
                  df = df2 * rinv;
                  dg = dg2;
               }

            } else if (dield == -5) {

               /* non-bonded term from yammp  */

               dis = r;
               ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
               d0 = prm->Cn2[ic];
               if (dis < d0) {
                  kij = prm->Cn1[ic];
                  diff = dis - d0;
                  *enb += kij * diff * diff;
                  df = 2.0 * kij * diff;
                  dg = 2.0 * kij;
               } else {
                  df = 0.0;
                  dg = 0.0;
               }
            } else {

               /*
                * Code for various dielectric models.
                * The df2 variable should hold r(dV/dr).
                */

               if (dield == 0) {

                  /* epsilon = r  */

                  rs = cgi * cgj * r2inv;
                  df2 = -2.0 * rs;
                  dg2 = 6.0 * rs * r2inv;
                  *eel += rs;

               } else if (dield == 1) {

                  /* epsilon = 1  */

                  rs = cgi * cgj * rinv;
                  df2 = -rs;
                  dg2 = 2.0 * rs * r2inv;
                  *eel += rs;

               } else if (dield == -2) {

                  /* Ramstein & Lavery dielectric, PNAS 85, 7231 (1988). */

                  rs = SIG * r;
                  rssq = rs * rs;
                  pow = exp(-rs);
                  eps1 = rssq + rs + rs + 2.0;
                  epsi = 1.0 / (DIW - C1 * pow * eps1);
                  cgijr = cgi * cgj * rinv * epsi;
                  *eel += cgijr;
                  df2 = -cgijr * (1.0 + C1 * pow * rs * rssq * epsi);
                  dg2 =
                      cgijr * (epsi *
                               (C1 * pow * SIG *
                                (rssq *
                                 (2.0 *
                                  (epsi * C1 * pow * SIG * rssq + rinv)
                                  + SIG) - 2.0 * SIG * rs)) + 2.0 * r2inv);
               }

               /* Calculate either Van der Waals or hydrogen bonded term. */

               ic = prm->Cno[iaci + prm->Iac[j] - 1];
               if (ic > 0 || enbfaci != 1.0) {
                  if (ic > 0) {
                     ic--;
                  } else {
                     ibig = prm->Iac[i] > prm->Iac[j] ?
                         prm->Iac[i] : prm->Iac[j];
                     isml = prm->Iac[i] > prm->Iac[j] ?
                         prm->Iac[j] : prm->Iac[i];
                     ic = ibig * (ibig - 1) / 2 + isml - 1;
                  }
                  r6 = r2inv * r2inv * r2inv;
                  f2 = prm->Cn2[ic] * r6;
                  f1 = prm->Cn1[ic] * r6 * r6;
                  *enb += (f1 - f2) * enbfaci;
                  df = (df2 + (6.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
                  dg = dg2 + (156.0 * f1 - 42.0 * f2) * enbfaci * r2inv;
#if 0
                  if (enbfaci != 1.0)
                     nb14 += (f1 - f2) * enbfaci;
#endif
               } else {
                  ic = -ic - 1;
                  r10 = r2inv * r2inv * r2inv * r2inv * r2inv;
                  f2 = prm->HB10[ic] * r10;
                  f1 = prm->HB12[ic] * r10 * r2inv;
                  *enb += (f1 - f2) * enbfaci;
                  df = (df2 + (10.0 * f2 - 12.0 * f1) * enbfaci) * rinv;
                  dg = dg2 + (156.0 * f1 - 110.0 * f2 * enbfaci) * r2inv;
                  nhbpair++;
#if 0
                  hbener += (f1 - f2) * enbfaci;
#endif
               }
            }

            /*
             * The df term includes Dij in the denominator, and the dg term
             * includes Dij^2 in the denominator; hence, terms such as d2ii
             * and di to not include Dij.
             */

            df *= rinv;
            dg *= r2inv;

            /* Update the gradient and the Jacobian. */

            dumx += df * di[0];
            dumy += df * di[1];
            dumz += df * di[2];

#ifdef SCALAPACK

            ptr = ptr1d(f, descF_PxQ, 3 * j + 0);
            if (ptr != NULL)
               *ptr += df * dj[0];
            ptr = ptr1d(f, descF_PxQ, 3 * j + 1);
            if (ptr != NULL)
               *ptr += df * dj[1];
            ptr = ptr1d(f, descF_PxQ, 3 * j + 2);
            if (ptr != NULL)
               *ptr += df * dj[2];
#else
            f[3 * j] += df * dj[0];
            f[3 * j + 1] += df * dj[1];
            f[3 * j + 2] += df * dj[2];
#endif

            /*
             * Calculate the second derivatives of Dij for the Van der Waals
             * and Coulombic interactions.
             */

            for (k = 0; k < 3; k++) {
               for (l = 0; l < 3; l++) {
                  ijterm = df * d2ij[k][l] + dg * di[k] * dj[l];

#ifdef SCALAPACK

                  ptr = ptr2d(g, descG_PxQ, 3 * i + l, 3 * i + k);
                  if (ptr != NULL)
                     *ptr += df * d2ii[k][l] + dg * di[k] * di[l];
                  ptr = ptr2d(g, descG_PxQ, 3 * j + l, 3 * j + k);
                  if (ptr != NULL)
                     *ptr += df * d2jj[k][l] + dg * dj[k] * dj[l];
                  ptr = ptr2d(g, descG_PxQ, 3 * j + l, 3 * i + k);
                  if (ptr != NULL)
                     *ptr += ijterm;
                  ptr = ptr2d(g, descG_PxQ, 3 * i + k, 3 * j + l);
                  if (ptr != NULL)
                     *ptr += ijterm;
#else
                  g[3 * i + k + n * (3 * i + l)] +=
                      df * d2ii[k][l] + dg * di[k] * di[l];
                  g[3 * j + k + n * (3 * j + l)] +=
                      df * d2jj[k][l] + dg * dj[k] * dj[l];
                  g[3 * i + k + n * (3 * j + l)] += ijterm;
                  g[3 * j + l + n * (3 * i + k)] += ijterm;
#endif
               }
            }
         }
      }

      /* For atom i, the gradient is updated in the i-loop only. */

#ifdef SCALAPACK

      ptr = ptr1d(f, descF_PxQ, 3 * i + 0);
      if (ptr != NULL)
         *ptr += dumx;

      ptr = ptr1d(f, descF_PxQ, 3 * i + 1);
      if (ptr != NULL)
         *ptr += dumy;

      ptr = ptr1d(f, descF_PxQ, 3 * i + 2);
      if (ptr != NULL)
         *ptr += dumz;

#else

      f[3 * i] += dumx;
      f[3 * i + 1] += dumy;
      f[3 * i + 2] += dumz;

#endif

   }

   /* Deallocate the iexw array. */

   free_ivector(iexw, -1, prm->Natom);

   return (0);
}

/***********************************************************************
                            EGB2()
************************************************************************/

/*
 * Calculate the generalized Born energy, 1st and 2nd derivatives.
 *
 * Calling parameters are as follows:
 *
 * lpears - the number of pairs on the lower triangle pair list
 * upears - the number of pairs on the upper trianble pair list
 * pearlist - the pair list, contiguous for the upper and lower triangles
 * lpearsnp - number of pairs on the non-polar lower triangle pair list
 * upearsnp - number of pairs on the non-polar upper trianble pair list
 * pearlistnp - non-polar pair list, contiguous for upper & lower triangles
 * x - input: the atomic (x,y,z) coordinates
 * f - updated: the gradient vector
 * g - updated: the Hessian matrix
 * grad - updated: the gradient vector in context_1x1
 * fs - input: overlap parameters
 * rborn - input: atomic radii
 * q - input: atomic charges
 * kappa - input: inverse of the Debye-Huckel length
 * diel_ext - input: solvent dielectric constant
 * enb - updated: Lennard-Jones energy
 * eelt - updated: gas-phase electrostatic energy
 * esurf - updated: nonpolar surface area solvation free energy
 * enp - updated: nonpolar van der Waals solvation free energy
 * context_Nx1 - input: the general context for ScaLAPACK
 * context_1x1 - input: the non-distributed vector context for ScaLAPACK
 * context_PxQ - input: the distributed vector and matrix context for ScaLAPACK
 * desc_3Nx1 - input: ScaLAPACK descriptor for the gradient vector in context_PxQ
 * desc_1x1 - input: ScaLAPACK descriptor for the gradient vector in context_1x1
 * desc_3Nx3N - input: ScaLAPACK descriptor for the Hessian matrix in context_PxQ
 * gridim - input: the ScaLAPACK process grid dimension (=1 for single task)
 * freevectors - input: if !=0 free the static vectors and return
 */

static
REAL_T egb2(INT_T * lpears, INT_T * upears, INT_T ** pearlist,
            INT_T * lpearsnp, INT_T * upearsnp, INT_T ** pearlistnp,
            REAL_T * x, REAL_T * f, REAL_T * g, REAL_T * grad,
            REAL_T * fs, REAL_T * rborn, REAL_T * q, REAL_T * kappa,
            REAL_T * diel_ext, REAL_T * enb, REAL_T * eelt, REAL_T * esurf,
            REAL_T * enp,
            INT_T context_Nx1, INT_T context_1x1, INT_T context_PxQ,
            INT_T * desc_3Nx1, INT_T * desc_1x1, INT_T * desc_3Nx3N,
            INT_T gridim, INT_T freevectors)

  /*FGB taylor coefficients follow */
  /* from A to L :                 */
  /* 1/3 , 2/5 , 3/7 , 4/9 , 5/11  */
  /* 4/3 , 12/5 , 24/7 , 40/9 , 60/11 */
  /* 20/3, 84/5, 216/7, 440/9, 780/11 */
#define TA 0.33333333333333333333
#define TB 0.4
#define TC 0.42857142857142857143
#define TD 0.44444444444444444444
#define TDD 0.45454545454545454545
#define TE 1.33333333333333333333
#define TF 2.4
#define TG 3.42857142857142857143
#define TH 4.44444444444444444444
#define THH 5.45454545454545454545
#define TI 6.66666666666666666667
#define TJ 16.8
#define TK 30.8571428571428571429
#define TL 48.8888888888888888889
#define TLL 70.9090909090909090909
{

#ifdef SCALAPACK
   int zero = 0, one = 1, info = 0;
   int myrow, mycol, nprow, npcol;
   int bs, bs3, ierror, lld;
   static int desc_NxN[DLEN_], desc_Nx3N[DLEN_];
   static SUPERD_T *super_3Nx3N = NULL, *super_NxN = NULL, *super_Nx3N =
       NULL;
   static REAL_T *reductarr = NULL, *gg = NULL;
   static REAL_T *sg = NULL, *sh = NULL, *sr = NULL, *sn = NULL;
   REAL_T *ptr;
   size_t ii, locp;
#endif

   int i, j, k, l, m, n, m3;
   int npairs, ic, iaci, threadnum, numthreads;
   int *iexw;
   size_t in, i3, in3, jn, j3, jn3, nn, n3, nn3;
   REAL_T elec, evdw, epol, dielfac, qi, qj, qiqj, fgbi, fgbk, rb2, expmkf;
   REAL_T xi, yi, zi, xj, yj, zj, xij, yij, zij;
   REAL_T dumbo, tmpsd;
   REAL_T df, dg, dx, dy, dz;
   REAL_T dedx, dedy, dedz;
   REAL_T dij1i, dij2i, dij3i, temp1, temp2;
   REAL_T tgb21, t1, t2;
   REAL_T qi2h, qid2h, datmp, da2tmp;
   REAL_T theta, ri1i;

   REAL_T dij, sumi, ijterm;
   REAL_T eel, f6, f12, rinv, r2inv, r4inv, r6inv;
   REAL_T r2, ri, ri2, riinv, rj, rjinv, rirjinv, rj2, sj, sj2, thi;
   REAL_T uij, efac, temp4, temp5, temp6, temp7, sumda;
   REAL_T d_dij, e_dij, f_dij, g_dij, h_dij, i_dij;
   REAL_T j_dij, l_dij, m_dij, n_dij;
   REAL_T a_dij, b_dij, c_dij, p_dij, r_dij;
   REAL_T temp8, temp9, temp10, temp11, temp12, temp13;
   REAL_T temp14, temp15, temp16, temp17, temp18;
   REAL_T temp20, temp21, temp22, temp23;

   REAL_T rgbmax1i, rgbmax2i, rgbmaxpsmax2;
   REAL_T dxj, dyj, dzj;
   REAL_T di[3], dj[3], d2ii[3][3], d2jj[3][3], d2ij[3][3];

   char transa, transb;
   REAL_T dblone = 1.0;

   static REAL_T *reff = NULL, *dreff = NULL, *psi = NULL;
   static REAL_T *sumdeijda = NULL, *sumdeijdn = NULL;
   static REAL_T *sumdeijdg = NULL, *sumdeijdh = NULL;

   REAL_T evdwnp, vdwdenom, vdwterm;

   t1 = seconds();
   tgb21 = t1;
   /* If PRINT_EGB2_TIMES is defined print some calculation times. */

#undef PRINT_EGB2_TIMES


   /*
    * If ROW_CYCLIC is defined allocate the Hessian and some temporary
    * matrices in a row cyclic manner via the "first touch" method.
    */

#define ROW_CYCLIC

   if (gb2_debug) {
      printf("enter egb2, n = %d\n", prm->Natom);
#ifdef SCALAPACK
      printf("egb2: context_Nx1=%d  context_1x1=%d  context_PxQ=%d\n",
             context_Nx1, context_1x1, context_PxQ);
      printf
          ("egb2: descF_PxQ[CTXT_]=%d  descF_1x1[CTXT_]=%d  descG_PxQ[CTXT_]=%d\n",
           desc_3Nx1[CTXT_], desc_1x1[CTXT_], desc_3Nx3N[CTXT_]);
#endif
      fflush(stdout);
   }


   /* Make sure that we have the second derivatives available */
   if( gb>1 ){
      fprintf( stderr, "Cannot do second derivatives with gb>1\n" );
      exit(1);
   }

   /* Compute these because they are used a fair amount. */

   n = prm->Natom;
   nn = n * n;
   m3 = 3 * n;
   n3 = m3;
   nn3 = n * n3;

   /*
    * If freevectors != 0, deallocate the static arrays that have been
    * previously allocated and return.  Note that when ScaLAPACK is
    * defined, the size of the arrays is determined by the 'size' field
    * of the associated superdescriptor, and so in principle the array
    * should be free before its associated superdescriptor.  In practice,
    * however, the free_vector function doesn't use the size argument,
    * so any dummy value could be supplied.
    */

   if (freevectors != 0) {

#ifndef SCALAPACK

      free_vector(sumdeijdn, 0, nn);
      sumdeijdn = NULL;

      free_vector(sumdeijdg, 0, nn3);
      sumdeijdg = NULL;

      free_vector(sumdeijdh, 0, nn3);
      sumdeijdh = NULL;

      free_vector(dreff, 0, nn3);
      dreff = NULL;

#else

      if (context_PxQ >= 0) {

         free_superdesc(super_3Nx3N);
         super_3Nx3N = NULL;

         free_vector(sumdeijdn, 0, super_NxN->size);
         sumdeijdn = NULL;

         free_superdesc(super_NxN);
         super_NxN = NULL;

         free_vector(sumdeijdg, 0, super_Nx3N->size);
         sumdeijdg = NULL;

         free_vector(sumdeijdh, 0, super_Nx3N->size);
         sumdeijdh = NULL;

         free_vector(dreff, 0, super_Nx3N->size);
         dreff = NULL;

         free_superdesc(super_Nx3N);
         super_Nx3N = NULL;

         free_vector(reductarr, 0, 3 * n3);
         reductarr = NULL;
      }

      free_vector(gg, 0, 3 * n3);
      gg = NULL;

      free_vector(sg, 0, n3);
      sg = NULL;

      free_vector(sh, 0, n3);
      sh = NULL;

      free_vector(sr, 0, n3);
      sr = NULL;

      free_vector(sn, 0, n);
      sn = NULL;

#endif

      free_vector(reff, 0, n);
      reff = NULL;

      free_vector(sumdeijda, 0, n);
      sumdeijda = NULL;

      free_vector(psi, 0, n);
      psi = NULL;

      if (gb2_debug) {
         printf("freevectors\n");
         fflush(stdout);
      }

      return (0.0);
   }

   /* Allocate some static arrays if they have not been allocated already. */

   if (reff == NULL)
      reff = vector(0, n);
   if (sumdeijda == NULL)
      sumdeijda = vector(0, n);
   if ((psi == NULL) && (gb == 2 || gb == 5))
      psi = vector(0, n);

#ifndef SCALAPACK

   /*
    * If SCALAPACK is not defined, allocate sumdeijdn, sumdeijdg, sumdeijdh
    * and dreff as NAB vectors.
    */

   if (sumdeijdn == NULL)
      sumdeijdn = vector(0, nn);
   if (sumdeijdg == NULL)
      sumdeijdg = vector(0, nn3);
   if (sumdeijdh == NULL)
      sumdeijdh = vector(0, nn3);
   if (dreff == NULL)
      dreff = vector(0, nn3);

#else

   /*
    * If SCALAPACK is defined, initialize matrix descriptors for the
    * sumdeindn, sumdeijdg, sumdeijdh and dreff matrices, then create
    * superdescriptors that contain the matrix size (for each process),
    * as well as precomputed division and modulus arrays to accelerate
    * the mapping from a global matrix (row, column) address to a local
    * submatrix array address.  Once the superdescriptors are created,
    * allocate the submatrix arrays.
    *
    * Get the block size for the square block and multiply by three.
    */

   bs = blocksize;
   bs3 = 3 * bs;

   /*
    * Get the number of rows and columns on the block cyclic (PxQ) process
    * grid, as well as this task's row and column on the grid.
    */

   blacs_gridinfo_(&context_PxQ, &nprow, &npcol, &myrow, &mycol);

   /*
    * Allocate the reduction array in all processes because MPI_Allreduce
    * will use MPI_COMM_WORLD as its communicator.  It is unclear that
    * it is possible to use a ScaLAPACK context as an MPI communicator,
    * so the default MPI_COMM_WORLD is used instead.
    */

   if (reductarr == NULL)
      reductarr = vector(0, 3 * n3);

   /*
    * If this task is on the process grid, set up the array descriptors.
    * If this task isn't on the process grid, set desc_NxN[CTXT_] and
    * desc_Nx3N[CTXT_] to -1 in case any of the associated matrices are
    * called by the pdgemr2d_ function, and set the superdescriptors
    * to NULL.
    */

   /*
    * Distribute the sumdeijdg, sumdeijdh and dreff matrices on the
    * block cyclic process grid.  Because the number of columns is
    * three times the number of rows for these matrices, use a row
    * block size of blocksize and a column block size of 3*blocksize
    * so that a test of the row index i and the column index j (not 3*j)
    * may be used to identify the grid process to which a global matrix
    * address maps.
    *
    * Distribute the sumdeijdn matrix on the block cyclic process
    * grid.  Use row and column block sizes of blocksize because
    * the number of rows equals the number of columns for this matrix.
    *
    * The numroc_ function is used to calculate the number of matrix
    * elements that are distributed across a PxQ processor grid.
    *
    * Once the descriptor has been initialized by the descinit function,
    * create a superdescriptor from the descriptor.
    *
    * Allocate the gg, sg, sh, sr and sn vectors for the diagonal
    * elements of the g, sumdeijdg, sumdeijdh, sumdeijdr and
    * sumdeijdn matrices.
    */

   if (context_PxQ >= 0) {

      if (super_3Nx3N == NULL) {
         if ((super_3Nx3N = superdesc(desc_3Nx3N)) == NULL) {
            printf("egb2: superdescriptor 3Nx3N == NULL, mytaskid = %d\n",
                   mytaskid);
            fflush(stdout);
         }
      }

      if (super_NxN == NULL) {
         locp = numroc_(&n, &bs, &myrow, &zero, &nprow);
         lld = locp;
         descinit_(desc_NxN, &n, &n, &bs, &bs, &zero, &zero,
                   &context_PxQ, &lld, &info);
         if (info) {
            printf
                ("egb2: error in desc_NxN initialization, mytaskid = %d\n",
                 mytaskid);
            fflush(stdout);
         }
         if ((super_NxN = superdesc(desc_NxN)) == NULL) {
            printf("egb2: superdescriptor NxN == NULL, mytaskid = %d\n",
                   mytaskid);
            fflush(stdout);
         }
      }
      if (sumdeijdn == NULL && super_NxN->size > 0) {
         sumdeijdn = vector(0, super_NxN->size);
      }

      if (super_Nx3N == NULL) {
         locp = numroc_(&n, &bs, &myrow, &zero, &nprow);
         lld = locp;
         descinit_(desc_Nx3N, &n, &m3, &bs, &bs3, &zero, &zero,
                   &context_PxQ, &lld, &info);
         if (info) {
            printf
                ("egb2: error in desc_Nx3N initialization, mytaskid = %d\n",
                 mytaskid);
            fflush(stdout);
         }
         if ((super_Nx3N = superdesc(desc_Nx3N)) == NULL) {
            printf("egb2: superdescriptor Nx3N == NULL, mytaskid = %d\n",
                   mytaskid);
            fflush(stdout);
         }
      }
      if (sumdeijdg == NULL && super_Nx3N->size > 0) {
         sumdeijdg = vector(0, super_Nx3N->size);
      }
      if (sumdeijdh == NULL && super_Nx3N->size > 0) {
         sumdeijdh = vector(0, super_Nx3N->size);
      }
      if (dreff == NULL && super_Nx3N->size > 0) {
         dreff = vector(0, super_Nx3N->size);
      }
   } else {
      desc_NxN[CTXT_] = desc_Nx3N[CTXT_] = -1;
      super_3Nx3N = super_NxN = super_Nx3N = NULL;
   }

   gg = vector(0, 3 * n3);
   sr = vector(0, n3);
   sg = vector(0, n3);
   sh = vector(0, n3);
   sn = vector(0, n);

   if (gb2_debug) {
      printf("vector allocation, mytaskid = %d\n", mytaskid);
      fflush(stdout);
   }
#endif

   /*
    * Smooth "cut-off" in calculating GB effective radii.
    * Implementd by Andreas Svrcek-Seiler and Alexey Onufriev.
    * The integration over solute is performed up to rgbmax and includes
    * parts of spheres; that is an atom is not just "in" or "out", as
    * with standard non-bonded cut.  As a result, calclated effective
    * radii are less than rgbmax. This saves time, and there is no
    * discontinuity in dReff/drij.
    *
    * Only the case rgbmax > 5*max(sij) = 5*fsmax ~ 9A is handled; this is
    * enforced in mdread().  Smaller values would not make much physical
    * sense anyway.
    *
    * Note: rgbmax must be less than or equal to cut so that the pairlist
    * generated from cut may be applied to calculation of the effective
    * radius and its derivatives.
    */

   if (rgbmax > cut) {
      printf("Error in egb2: rgbmax = %f is greater than cutoff = %f\n",
             rgbmax, cut);
      exit(1);
   }
   rgbmax1i = 1.0 / rgbmax;
   rgbmax2i = rgbmax1i * rgbmax1i;
   rgbmaxpsmax2 = (rgbmax + prm->Fsmax) * (rgbmax + prm->Fsmax);

#if 0
   if (gb2_debug) {
      printf("Effective Born radii:\n");
   }
#endif

   /* 
    * Get the "effective" Born radii via the approximate pairwise method
    * Use Eqs 9-11 of Hawkins, Cramer, Truhlar, J. Phys. Chem. 100:19824
    * (1996).
    *
    * For SCALAPACK or MPI, initialize all elements of the reff array
    * to 0 because although each task will calculate only a subset of the
    * elements, a reduction is used to combine the results from all tasks.
    * Initialize this array even for tasks that are not on the process grid
    * because MPI_Allreduce will use the default MPI_COMM_WORLD.
    *
    * For OpenMP, the combination of all threads loads all elements
    * of the reff array so there is no need to initialize it to 0.
    *
    * For single-threaded execution all elements of the reff array are
    * loaded so there is no need to initialize it to 0.
    */

#if defined(SCALAPACK) || defined(MPI)

   for (i = 0; i < n; i++) {
      reff[i] = 0.0;
   }

   if (gb2_debug) {
      printf("reff initialized, mytaskid = %d\n", mytaskid);
      fflush(stdout);
   }
#endif

#pragma omp parallel \
  private (i, xi, yi, zi, ri, ri1i, sumi, j, k, xij, yij, zij, \
           r2, dij1i, dij, sj, sj2, uij, dij2i, tmpsd, dumbo, theta, \
           threadnum, numthreads )
   {
      /*
       * Get the thread number and the number of threads for multi-threaded
       * execution under OpenMP.  For ScaLAPACK, use the process column and
       * the number of process columns, respectively.  For single-threaded
       * execution and MPI there is no need to set these variables.
       */

#ifdef OPENMP

      threadnum = omp_get_thread_num();
      numthreads = omp_get_num_threads();

#elif defined(SCALAPACK)

      threadnum = mycol;
      numthreads = npcol;

#endif

      /*
       * Loop over all atoms i.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       * 
       * OpenMP utilizes all threads.  ScaLAPACK utilizes all processes
       * on the process grid.  MPI is not parallelized.
       *
       * The reff array is written in the following loops.  It is necessary
       * to synchronize the OpenMP threads that execute these loops after
       * loop execution so that a race condition does not exist for reading
       * the reff array before it is written.  Even if all subsequent loops
       * use loop index to thread mapping that is identical to that of the
       * following loop, elements of the reff array are indexed by other
       * than such a loop indices, so synchronization is necessary.
       *
       * OpenMP synchronization is accomplished by the implied barrier
       * at the end of this parallel region.  ScaLAPACK synchronization is
       * accomplished by MPI_Allreduce.
       */

      for (i = 0; i < n; i++) {

#if defined(OPENMP) || defined(SCALAPACK)

         if (!myroc(i, blocksize, numthreads, threadnum))
            continue;
#endif

         xi = x[3 * i];
         yi = x[3 * i + 1];
         zi = x[3 * i + 2];
         ri = rborn[i] - gboffset;
         ri1i = 1. / ri;
         sumi = 0.0;

         /* Select atom j from the pair list.  Non-graceful error handling. */

         for (k = 0; k < lpears[i] + upears[i]; k++) {

            if (pearlist[i] == NULL) {
               printf
                   ("NULL pair list entry in egb2 loop 1, mytaskid = %d\n",
                    mytaskid);
               fflush(stdout);
            }
            j = pearlist[i][k];

#ifdef SCALAPACK

            if (gb2_debug) {
               if (!myroc(j, blocksize, nprow, myrow)) {
                  printf("pair list error!\n");
                  fflush(stdout);
               }
            }
#endif

            xij = xi - x[3 * j];
            yij = yi - x[3 * j + 1];
            zij = zi - x[3 * j + 2];
            r2 = xij * xij + yij * yij + zij * zij;
            if (r2 > rgbmaxpsmax2)
               continue;
            dij1i = 1.0 / sqrt(r2);
            dij = r2 * dij1i;
            sj = fs[j] * (rborn[j] - gboffset);
            sj2 = sj * sj;

            /*
             * The following equations are from the Appendix of Schaefer and Froemmel,
             * JMB 216:1045-1066, 1990.  Taylor series expansion for d>>s is by Andreas
             * Svrcek-Seiler.  Smooth rgbmax idea is from Andreas Svrcek-Seiler and
             * Alexey Onufriev.
             */

            if (dij > rgbmax + sj)
               continue;

            if (dij > rgbmax - sj) {

               uij = 1. / (dij - sj);

               sumi -= 0.125 * dij1i * (1.0 + 2.0 * dij * uij +
                                        rgbmax2i * (r2 -
                                                    4.0 * rgbmax * dij -
                                                    sj2) -
                                        2.0 * log(rgbmax * uij));

            } else if (dij > 4.0 * sj) {

               dij2i = dij1i * dij1i;
               tmpsd = sj2 * dij2i;
               dumbo =
                   TA + tmpsd * (TB +
                                 tmpsd * (TC +
                                          tmpsd * (TD + tmpsd * TDD)));

               sumi -= sj2 * sj * dij2i * dij2i * dumbo;

            } else if (dij > ri + sj) {

               sumi -= 0.5 * (sj / (r2 - sj2) +
                              0.5 * dij1i * log((dij - sj) / (dij + sj)));

            } else if (dij > fabs(ri - sj)) {

               theta = 0.5 * ri1i * dij1i * (r2 + ri * ri - sj2);
               uij = 1. / (dij + sj);

               sumi -= 0.25 * (ri1i * (2. - theta) - uij +
                               dij1i * log(ri * uij));

            } else if (ri < sj) {

               sumi -= 0.5 * (sj / (r2 - sj2) + 2. * ri1i +
                              0.5 * dij1i * log((sj - dij) / (sj + dij)));

            }

         }

#define GBALPHA 0.8
#define GBBETA 0.0
#define GBGAMMA 2.909125

         /*
          * For ScaLAPACK store sums instead of reciprocal sums
          * because more terms will be contributed to each sum
          * from other process columns via the subsequent reduction.
          */

#ifdef SCALAPACK

         if (gb == 1) {

            /* "standard" (HCT) effective radii:  */

            reff[i] = sumi;

         } else {

            /* "gbao" formulas:  */

            psi[i] = -ri * sumi;

            reff[i] =
                -tanh((GBALPHA +
                       GBGAMMA * psi[i] * psi[i]) * psi[i]) / rborn[i];
         }
#else

         if (gb == 1) {

            /* "standard" (HCT) effective radii:  */

            reff[i] = 1.0 / (ri1i + sumi);

         } else {

            /* "gbao" formulas:  */

            psi[i] = -ri * sumi;

            reff[i] = ri1i - tanh((GBALPHA + GBGAMMA * psi[i] * psi[i])
                                  * psi[i]) / rborn[i];
         }
#endif

      }
   }

   if (gb2_debug) {
      printf("reff calculated, mytaskid = %d\n", mytaskid);
      fflush(stdout);
   }

   /*
    * The ScaLAPACK synchronization is accomplished via reduction
    * of the reff array.
    *
    * No reduction is accomplished for MPI because egb2 is not
    * parallelized for MPI.
    */

#ifdef SCALAPACK

   ierror = MPI_Allreduce(reff, reductarr, n,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (ierror != MPI_SUCCESS) {
      printf
          ("Error in egb reff reduction, error = %d  myrow = %d  mycol = %d\n",
           ierror, myrow, mycol);
   }

   if (gb2_debug) {
      printf("MPI_Allreduce reff, mytaskid = %d\n", mytaskid);
      fflush(stdout);
   }

   for (i = 0; i < n; i++) {
      reff[i] = reductarr[i];
      if (reff[i] == 0.0) {
         printf("reff[%d]=0!, myrow = %d  mycol = %d\n", i, myrow, mycol);
         fflush(stdout);
      }
   }

   /* For ScaLAPACK reciprocate the reduced sums in reff. */

   for (i = 0; i < n; i++) {
      ri = rborn[i] - gboffset;
      ri1i = 1. / ri;
      reff[i] = 1.0 / (ri1i + reff[i]);
   }

#endif

   *esurf = 0.0;

   /*
    * Compute the generalized Born energy for all pair interactions, as
    * well as the Van der Waals and Coulombic energy for all pairs that
    * are not on the excluded list.  Compute the 1st and 2nd derivatives
    * of the of the Van der Waals and Coulombic energy for all non-excluded
    * pairs.
    *
    * In addition, compute the 1st derivative of the Born energy for all
    * pairs of atoms (selected by the indices i and j).  David Case and
    * his colleagues have shown that the first derivative for an atom pair
    * may be calculated by first building the following equations in
    * terms of Ri (the effective radius of atom i), Rj (the effective
    * radius of atom j) and Dij (the distance between atoms i and j):
    *
    *      rb2 = Ri * Rj
    *
    *      efac = exp(Dij^2 / (4 * rb2))
    *
    *      fgbi = 1 / sqrt(Dij^2 + rb2 * efac)
    *
    *      fgbk = -kappa * KSCALE / fgbi
    *
    *      expmkf = exp(fgbk) / DIELECTRIC_CONSTANT
    *
    *      dielfac = 1 - expmkf
    *
    *      temp4 = fgbi * fgbi * fgbi
    *
    *      temp6 = Qi * Qj * temp4 * (dielfac + fgbk * expmkf)
    *
    *      temp5 = efac * temp6 * (rb2 + 0.25 * r2) / 2
    *
    * From the above equations it is possible to express the first
    * derivative of the Born energy with respect to any cartesian coordinate
    * (designated by Xl) by summing contributions from the three following
    * equations:
    *
    * 1.   d(Dij)/d(Xl) * A(Dij) where
    *
    *      A(Dij) = (1 - efac / 4) * temp6
    *
    * 2.   d(Ri)/d(Xl) * B(Dij) where
    *
    *      B(Dij) = Ri * temp5
    *
    * 3.   d(Rj)/d(Xl) * C(Dij) where
    *
    *      C(Dij) = Rj * temp5
    *
    * Equation 1 is the simplest of the three to evaluate.  One calculates
    * all of the derivatives of Dij with respect to cartesian coordinates,
    * multiplies these derivatives by A(Dij), and sums the resulting products
    * to the gradient vector.  There are six derivatives involving atoms
    * i and j:
    *
    *      d(Dij)/d(Xi), d(Dij)/d(Yi), d(Dij)/d(Zi)
    *
    *      d(Dij)/d(Xj), d(Dij)/d(Yj), d(Dij)/d(Zj)
    *
    * The derivatives of Dij with respect to the cartesian coordinates
    * contain a factor of Dij in the denominator.  This factor is absorbed
    * by equation 1 (and all of the other equations that incorporate the
    * derivatives) rather than carried with the individual derivatives.
    *
    * Equations 2 and 3 are more complicated because the derivatives of the
    * effective radii Ri and Rj with respect to the cartesian coordinates are
    * sometimes series, but sometimes single terms.  For example:
    *
    *      d(Ri)/d(Xi) = d[f(Dij)]/d(Dij) * d(Dij)/d(Xi)
    *
    *                  + d[f(Dik)]/d(Dik) * d(Dik)/d(Xi) + ...
    *
    *                  + d[f(Din)]/d(Din) * d(Din)/d(Xi)
    *
    *      d(Ri)/d(Xj) = d[f(Dij)]/d(Dij) * d(Dij)/d(Xj)
    *
    *      d(Rj)/d(Xj) = d[f(Dij)]/d(Dij) * d(Dij)/d(Xj)
    *
    *                  + d[f(Djk)]/d(Djk) * d(Djk)/d(Xj) + ...
    *
    *                  + d[f(Djn)]/d(Djn) * d(Djn)/d(Xj)
    *
    *      d(Rj)/d(Xi) = d[f(Dij)]/d(Dij) * d(Dij)/d(Xi)
    *
    * where n is the number of atoms and f() designates a function of the
    * interatomic distance.
    *
    * Moreover, when all of the atom pairs are considered, one can refactor
    * the computation to obtain terms such as:
    *
    *      d(Ri)/d(Xi) * [B(Dij) + B(Dik) + ... + B(Din)]
    *
    *      d(Rj)/d(Xj) * [C(Dij) + C(Djk) + ... + C(Djn)]
    *
    * The square brackets [] indicate that B(Dij) + B(Dik) + ... + B(Din) are
    * summed prior to multiplication by the derivatives of the effective radius.
    *
    * David Case and his colleagues have discovered an algorithm that permits
    * evaluation of the first derivatives while requiring that the first
    * derivatives of the effective radii be evaluated once only.  They
    * use an array called "sumdeijda" into which the above series are summed.
    * For example, the series B(Dij) + B(Dik) + ... + B(Din) is summed into
    * sumdeijda[i], and the series C(Dij) + C(Djk) + ... + C(Din) is summed
    * into sumdeijda[j].  One can think of the contents of sumdeijda[i]
    * as the sum of everything that gets multiplied by d(Ri)/d(Xl), where
    * d(Ri)/d(Xl) represents all of the derivatives of Ri with respect to the
    * cartesian coordinates.  Similarly, one can think of the contents of
    * sumdeijda[j] as the sum of everything that gets multiplied by d(Rj)/d(Xl),
    * where d(Rj)/d(Xl) represents all of the derivatives of Rj with respect
    * to the cartesian coordinates.
    *
    * Once the sumdeijda array is loaded, the d(Ri)/d(Xl) and d(Rj)/d(Xl) are
    * calculated and multiplied by the respective elements of sumdeijda.  For
    * example, for d(Ri)/d(Xl) the individual terms of the d(Ri)/d(Xi) series
    * are calculated.  Each term has an x-, y- and z-component.  The components
    * of the jth term are:
    *
    *      d[f(Dij)]/d(Dij) * d(Dij)/dXi
    *
    *      d[f(Dij)]/d(Dij) * d(Dij)/dYi
    *
    *      d[f(Dij)]/d(Dij) * d(Dij)/dZi 
    *
    *      d[f(Dij)]/d(Dij) * d(Dij)/dXj
    *
    *      d[f(Dij)]/d(Dij) * d(Dij)/dYj
    *
    *      d[f(Dij)]/d(Dij) * d(Dij)/dZj
    *
    * Note that the jth term contains derivatives with respect to the cartesian
    * coordinates of atoms i and j.  The derivatives with respect to the cartesian
    * coordinates of atom j are stored in da[j].  The derivatives with respect to
    * the cartesian coordinates of atom i are summed into da[i] with the result
    * that da[i] contains derivatives from f(Dij), f(Dik), ..., f(Din).
    *
    * Each element of the da array is then multiplied by sumdeijda[i] to form:
    *
    * [ d(Rij)/d(Xl) + d(Rik)/d(Xl) + ... ] * [f(Dij) + f(Dik) + ... ]

    * where the series d(Rij)/d(Xl) + d(Rik)/d(Xl) + ... contains n terms each
    * of which has an (x,y,z) coordinate.  Therefore this result is a vector of
    * length 3*n where n is the number of atoms.  This vector represents
    * the derivatives of Ri with respect to all of the cartesian coordinates.
    * It is summed into the gradient vector.
    *
    * The process is then repeated for Rj, Rk, ..., Rn.
    *
    * There is one further refinement that relates to the first derivatives
    * of the "diagonal" terms, i.e., those energy terms that are functions
    * of the effective radius Ri but not of the interatomic distance Dij.
    * One begins by building up the following equations:
    *
    *      expmkf = exp(-kappa * KSCALE * Ri) / DIELECTRIC_CONSTANT
    *
    *      dielfac = 1 - expmkf
    *
    *      qi2h = Qi^2 / 2
    *
    *      qid2h = di2h * dielfac
    *
    * Then the first derivative with respect to Ri is the equation:
    *
    * 4.   d(Ri)/d(Xl) * R(Dij) where
    *
    *      R(Dij) = qid2h - kappa * KSCALE * qi2h * expmkf * Ri
    *
    *
    * Equation 4, like equation 2, is multiplied by d(Ri)/d(Xl) and
    * therefore one can add these two equations and store the result
    * in the ith element of the sumdeijda array.
    *
    * Note that equations 2, 3 and 4 include a factor of either -Ri*Ri
    * or -Rj*Rj which arises from the derivative of Ri or Rj, respectively.
    * The negation is introduced when  sumda is multiplied by an element
    * from the dreff array.
    *
    * The above algorithm may be extended to the calculation of the second
    * derivatives.  One proceeds by first building the following equations
    * in terms of Ri, Rj and Dij:
    *
    *      RiInv = 1 / Ri
    *
    *      RjInv = 1 / Rj
    *
    *      RiRjInv = RiInv * RjInv
    *
    *      temp4 = fgbi * fgbi * fgbi
    *
    *      temp5 = Qi * Qj * temp4 * (dielfac + fgbk * expmkf) / 2
    *
    *      temp6 = fgbi * fgbi * temp5
    *
    *      temp7 = Qi * Qj * kappa^2 * KSCALE^2 * temp4 * expmkf
    *
    *      temp8 = 1 - efac / 4
    *
    *      temp9 = temp7 / 4 - 3 * temp6 / 2
    *
    *      temp10 = Dij^2 / 4
    *
    *      temp11 = Ri * Rj + temp10
    *
    *      temp12 = Ri^2 * Rj^2
    *
    *      temp13 = Ri * temp11
    *
    *      temp14 = Rj * temp11
    *
    *      temp15 = temp10 * RjInv
    *
    *      temp16 = temp10 * RiInv
    *
    *      temp17 = temp9 * temp13 * efac
    *
    *      temp18 = temp9 * temp14 * efac
    *
    *      temp20 = temp8 * temp9
    *
    *      temp21 = temp20 * temp13
    *
    *      temp22 = temp20 * temp14
    *
    *      temp23 = 2 * efac;
    *
    *
    * From the above equations it is possible to express the second
    * derivative of the Born energy with respect to two cartesian coordinates
    * (designated by Xl and Xm) by summing contributions from the following
    * equations:
    *
    * 5.   d(Ri)/d(Xm) * [ d(Dij)/d(Xl) * D(Dij) ] where
    *
    *      D(Dij) = temp23 * (temp21 - temp15 * temp5 / 4)
    *
    *
    * 6.   d(Rj)/d(Xm) * [ d(Dij)/d(Xl) * E(Dij) ] where
    *
    *      E(Dij) = temp23 * (temp22 - temp16 * temp5 / 4)
    *
    *
    * 7.   d(Ri)/d(Xl) * [ d(Dij)/d(Xm) * J(Dij) ] where
    *
    *      J(Dij) = temp23 * ((Ri - RjInv * temp11) * temp5 / 4 + temp21)
    *
    *
    * 8.   d(Rj)/d(Xl) * [ d(Dij)/d(Xm) * N(Dij) ] where
    *
    *      N(Dij) = temp23 * ((Rj - RiInv * temp11) * temp5 / 4 + temp22)
    *
    *
    * 9.   d(Ri)/d(Xm) * [ d(Ri)/d(Xl) * H(Dij) ] where
    *
    *      H(Dij) = efac * (temp10 * (temp13 * RjInv - Ri^2) * temp5 + temp13 * temp17)
    *
    *
    * 10.  d(Rj)/d(Xm) * [ d(Rj)/d(Xl) * M(Dij) ] where
    *
    *      M(Dij) = efac * (temp10 * (temp14 * RiInv - Rj^2) * temp5 + temp14 * temp18)
    *
    *
    * 11.  d(Rj)/d(Xm) * [ d(Ri)/d(Xl) * I(Dij) ] where
    *
    *      I(Dij) = efac * ((temp12 + temp16 * temp13) * temp5 + temp13 * temp18)
    *
    *
    * 12.  d(Ri)/d(Xm) * [ d(Rj)/d(Xl) * L(Dij) ] where
    *
    *      L(Dij) = efac * ((temp12 + temp15 * temp14) * temp5 + temp14 * temp17)
    *
    *
    * 15.  d2(Dij)/d(Xm)d(Xl) * F(Dij) where
    *
    *      F(Dij) = 2 * temp8 * temp5
    *
    *
    * 16.  d(Dij)/d(Xm) * [ d(Dij)/d(Xl) * G(Dij) ] where
    *
    *      G(Dij) = 4 * temp8 * temp20
    *
    *               + (2 * temp8 + efac * temp10 * RiRjInv) * temp5 / Dij^2
    *
    *
    * Using the above equations 5-16 it is possible to calculate the
    * second derivatives of the Born energy in an optimized manner.
    * Note, however, that computing the second derivatives of the effective
    * radius with respect to the cartesian coordinates involves some
    * complexity.  This fact can be appreciated by beginning with the
    * following expression for the effective radius Ri of atom i:
    *
    *
    * Ri = 1 / [1/Rborn(i) + f(Dij) + f(Dik) + ... + f(Din)]
    *
    *
    * where Rborn(i) represents the Born radius of atom i and f() respesents
    * a function of the interatomic distance.  Differentiating with respect to
    * a cartesian coordinate Xl gives:
    *
    *
    * d(Ri)/d(Xl) = -Ri^2 * { d[f(Dij)]/d(Dij) * d(Dij)/d(Xl) +
    * 
    *                         d[f(Dik)]/d(Dik) * d(Dik)/d(Xl) + ... +
    *
    *                         d[f(Din)]/d(Din) * d(Din)/d(Xl) }
    *
    *
    * The above expression is the general case.  Many of the terms are
    * zero, for example, d(Dij)/d(Xl) is zero when l equals neither i nor j.
    *
    * Differentiating again with respect to a cartesian coordinate Xm gives:
    *
    *
    * d2(Ri)/d(Xm)d(Xl) =
    *
    * - Ri^2 * { d[f(Dij)]/d(Dij) * d2(Dij)/d(Xm)d(Xl) +
    *
    *            d[f(Dik)]/d(Dik) * d2(Dik)/d(Xm)d(Xl) + ... +
    *
    *            d[f(Din)]/d(Din) * d2(Din)/d(Xm)d(Xl) }
    *
    *
    * - Ri^2 * { d2[f(Dij)]/d(Dij)2 * d(Dij)/d(Xm) * d(Dij)/d(Xl) +
    *
    *            d2[f(Dik)]/d(Dik)2 * d(Dik)/d(Xm) * d(Dik)/d(Xl) + ... +
    *
    *            d2[f(Din)]/d(Din)2 * d(Din)/d(Xm) * d(Din)/d(Xl) }
    *
    *
    * + 2Ri^3 * { d[f(Dij)]/d(Dij) * d(Dij)/d(Xm) * d[f(Dij)]/d(Dij) * d(Dij)/d(Xl) +
    *
    *             d[f(Dij)]/d(Dij) * d(Dij)/d(Xm) * d[f(Dik)]/d(Dik) * d(Dik)/d(Xl) + ... +
    *
    *             d[f(Dij)]/d(Dij) * d(Dij)/d(Xm) * d[f(Din)]/d(Din) * d(Din)/d(Xl) +
    *
    * 
    *             d[f(Dik)]/d(Dik) * d(Dik)/d(Xm) * d[f(Dij)]/d(Dij) * d(Dij)/d(Xl) +
    *
    *             d[f(Dik)]/d(Dik) * d(Dik)/d(Xm) * d[f(Dik)]/d(Dik) * d(Dik)/d(Xl) + ... +
    *
    *             d[f(Dik)]/d(Dik) * d(Dik)/d(Xm) * d[f(Din)]/d(Din) * d(Din)/d(Xl) + ... +
    *
    *
    *             d[f(Din)]/d(Din) * d(Din)/d(Xm) * d[f(Dij)]/d(Dij) * d(Dij)/d(Xl) +
    *
    *             d[f(Din)]/d(Din) * d(Din)/d(Xm) * d[f(Dik)]/d(Dik) * d(Dik)/d(Xl) + ... +
    *
    *             d[f(Din)]/d(Din) * d(Din)/d(Xm) * d[f(Din)]/d(Din) * d(Din)/d(Xl) }
    *
    *
    * Many of the terms in the above expression are zero, or example, when neither
    * l nor m equals i or j.  However, many non-zero terms remain, and in particular
    * there exists a large number of the terms that are multiplied by 2Ri^3.  But
    * fortunately these particular terms may be computed by simply multiplying
    * the derivative of the Ri with respect to all cartesian coordinates by itself.
    *
    * As can be appreciated from the above analysis, differentiation gives rise to
    * factors of -Ri^2, Ri^3, -Rj^2 and Rj^3.  A factor of Ri^2 or Rj^2 is included
    * equations 5-16, as well as equations 17-18 (below).
    *
    * Square brackets [ ] have been added to many of equations 5-18
    * to indicate the order of evaluation.  The derivative inside of the
    * brackets will be called the inner derivative.  The derivative outside
    * of the brackets will be called the outer derivative.  The inner
    * derivative will be evaluated prior to the outer derivative.  The
    * inner derivative that is multiplied by the accompanying functions
    * of Dij, e.g., H(Dij) of equation 7, produces a gradient vector.
    *
    * Equations 15 and 16 are analogous to equation 1 and are simple to
    * evaluate because they do not involve derivatives of the effective
    * radii Ri and Rj and so do not comprise series.  Instead, they comprise
    * the derivatives of Dij with respect to the cartesian coordinates
    * of atoms i and j only.
    *
    * Equations 13 and 14 are analogous to equations 2 and 3 except that
    * they involve second derivatives of the effective radii Ri and Rj instead
    * of first derivatives.  However, the second derivatives of the effective
    * radius Ri include terms such as d2(Ri)/d(Xl)d(Xm) which affect many Hessian
    * elements when neither l nor m equals i.
    *
    * Equations 9, 10, 11 and 12 imply two applications of the approach used
    * to evaluate equations 2 and 3.  However, the result of the multiplication
    * by the inner derivative (which is enclosed in [ ]) is not a single term
    * but rather is a gradient vector.  Therefore, multiplication by the outer
    * derivative involves multiplying each derivative of Ri or Rj by all of the
    * terms in the gradient vector.  Moreover, this process requires a
    * two-dimensional array, i.e., an array of gradient vectors which is named
    * either "sumdeijdg" or "sumdeijdh", and is analogous to the sumdeijda array
    * that is used to evaluate equations 2, 3 and 4.  The order of evaluation of
    * the derivatives appears to be unimportant for equations 9 and 10, i.e., the
    * inner derivatives may be swappped  with the outer derivatives without
    * affecting the results.
    *
    * However, there is a subtlety involving equations 11 and 12 which requires
    * that if the inner derivatives are swapped with the outer derivatives in
    * equation 11, they must be similarly swapped in equation 12.  This subtlety
    * arises because the the sumdeijdg array is used to collect all inner
    * derivatives with respect to Xl that will subsequently be multiplied by
    * outer derivatives with respect to Xm  An analogous sumdeijdh array is used
    * to collect all inner derivatives with respect to Xm that will subsequently
    * be multiplied by outer derivatives with respect to Xl.  The inner derivatives
    * of equations 11 and 12 are stored in the sumdeijdg array, which implies
    * that for these equations the inner derivatives are taken with respect to Xl.
    * Swapping the inner and outer derivatives for both of these equations preserves
    * the symmetry of the Hessian matrix, but swapping the inner and outer derivatives
    * for only one of these equations breaks the symmetry and leads to erroneous
    * results.  For example, if one were to swap the inner and outer derivatives
    * in equation 11, but not in equation 12, the following pair of equations would
    * obtain:
    *
    * 11. d(Ri)/d(Xm) * [ d(Rj)/d(Xl) * I(Dij) ]
    *
    * 12.  d(Ri)/d(Xm) * [ d(Rj)/d(Xl) * L(Dij) ]
    *
    * This pair of equations is not symmetric and is therefore incorrect.
    *
    * Equations 5, 6, 7 and 8 are evaluated by first multiplying by
    * an inner derivative of Dij and then by multiplying by an outer
    * derivative of Ri or Rj.  The inner derivatives of Dij are
    * collected into either the sumdeijdg array or sumdeijdh array.
    * The derivatives of a particular Dij produce a gradient vector
    * that is stored in one row of either sumdeijdg or sumdeijdh.
    * Subsequent multiplication by the outer derivative of Ri or Rj
    * involves multiplying each derivative of Ri or Rj by all of the
    * terms in the gradient vector.  This process involves
    * a subtlety.  In equations 5 and 6 the inner derivative is with
    * respect to Xl whereas in equations 7 and 8 the inner derivative
    * is with respect to Xm.  For this reason, the inner derivatives
    * of equations 5 and 6 are collected into the sumdeijdg array,
    * whereas the inner derivatives of equations 7 and 8 are collected
    * into the sumdeijdh array, following the convention that was
    * discussed above for equations 11 and 12.  The outer derivatives
    * taken with respect to Xm are computed using the sumdeijdg array,
    * whereas the outer derivatives taken with respect to Xl are computed
    * using the sumdeijdh array.  These two sets of derivatives are
    * written into transpose-related elements of the Hessian matrix in
    * order to preserve the symmetry of the Hessian.
    *
    * On another topic, it is possible to generate equations for the second
    * derivatives of the non-pairwise ("diagonal") terms.  One begins by
    * building up the following equations:
    *
    *      expmkf = exp(-kappa * KSCALE * Ri) / DIELECTRIC_CONSTANT
    *
    *      dielfac = 1 - expmkf
    *
    *      temp18 = Ri * (kappa * KSCALE * expmkf * Ri * (2 + kappa * KSCALE * Ri)
    *
    *                     - 2 * dielfac)
    *
    *
    * Using the above equations it is possible to calculate the second
    * derivatives of the non-pairwise terms via the two following equations:
    *
    * 17.  d2(Ri)/d(Xm)d(Xl) * R(Dij) where R(Dij) is defined in equation 4.
    *
    *
    * 18.  d(Ri)/d(Xm) * [ d(Ri)/d(Xl) * P(Dij) ] where
    *
    *      P(Dij) = Qi^2 * temp18 / 2
    *
    * Equations 17 and 18 are handled in a manner analogous to equation 4.
    *
    *
    * Calculate the Born, van der Waals and Coulombic energies.  Complete
    * equations 15 and 16, and set up the remaining equations.  For OpenMP
    * in this collection of loops sum data relative only to atom i or to
    * pairs ii and ij so that OpenMP can be used to parallelize the outer loop.
    * For ScaLAPACK sum data relative to atoms i and j, as well as to pairs
    * ii and jj.  When in single-threaded mode sum data relative to atoms
    * i and j, as well as to pairs ii, ij, ji and jj.
    *
    * Note: the iexw vector is allocated and deallocated
    * within the following parallel region.  It may be useful to declare
    * a global array outside of the parallel region as is done for reff
    * sumdeijda and f.
    *
    * Initialize the energy by all OpenMP threads or ScaLAPACK processes.
    */

   epol = elec = evdw = evdwnp = 0.0;

#ifndef __PGI
#pragma omp parallel reduction (+: epol, elec, evdw, evdwnp) \
  private (i, i3, in, in3, xi, yi, zi, qi, qj, ri, ri2, riinv, expmkf, dielfac, \
           qi2h, qid2h, temp17, temp18, p_dij, r_dij, iaci, j, j3, jn, jn3, \
           xij, yij, zij, r2, r2inv, qiqj, rj, rj2, rjinv, rirjinv, rb2, efac, \
           fgbi, fgbk, temp4, temp6, a_dij, temp5, b_dij, c_dij, di, dj, d2ii, d2ij, \
           k, l, df, dg, rinv, r4inv, eel, ic, r6inv, f6, f12, dedx, dedy, dedz, \
           temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, \
           temp16, temp20, temp21, temp22, temp23, d_dij, e_dij, f_dij, \
           g_dij, h_dij, i_dij, j_dij, l_dij, m_dij, n_dij, \
           npairs, xj, yj, zj, thi, threadnum, numthreads, \
           d2jj, ri1i, temp1, dij3i, dx, dy, dz, m, iexw, \
           datmp, tmpsd, dumbo, temp2, da2tmp, sj, sj2, dij1i, dij2i, dij, \
           dxj, dyj, dzj, sumda, ijterm, vdwdenom, vdwterm)
#endif
   {

      /*
       * Get the thread number and the number of threads for multi-threaded
       * execution under OpenMP.  For ScaLAPACK, use the process column and
       * the number of process columns, respectively.  These variables are
       * not used for MPI or single-threaded execution.
       */

#ifdef OPENMP
      threadnum = omp_get_thread_num();
      numthreads = omp_get_num_threads();
#elif defined(SCALAPACK)
      threadnum = mycol;
      numthreads = npcol;
#endif

      /*
       * For OpenMP, initialize the sumdeijda, sumdeijdg, sumdeijdh,
       * sumdeijdn and dreff arrays inside of the parallel region.
       *
       * The "first touch" memory allocation strategy will locate
       * the elements of these arrays that are initialized by
       * a particular CPU in memory that is local to that CPU.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       *
       * It is not necessary to synchronize OpenMP threads following
       * this loop because the next loops use identical loop index to
       * thread mapping.  Hence the sumdeijda array is guaranteed to
       * be initialized above before it is updated below.
       *
       * For ScaLAPACK, each process maintains a full copy of the sumdeijda
       * and grad arrays, so initialize the entirety of these arrays for each
       * task.  Note that because each process maintains a full copy of these
       * arrays they may be updated in i,j loops only and not in j,i loops.
       *
       * For ScaLAPACK initialize as well the gg, sg, sh, sr and sn vectors.
       *
       * Because at present MPI_Allreduce uses the default MPI_COMM_WORLD,
       * initialization must be performed for all tasks, not just the tasks
       * on the process grid.
       *
       * Each process on the process grid maintains a submatrix of the
       * distributed sumdeijdg, sumdeijdh, sumdeijdn and dreff matrices,
       * so initialize the submatrix that is maintained by each process.
       */

#ifdef SCALAPACK

      for (i = 0; i < n; i++) {
         sumdeijda[i] = sn[i] = 0.0;
      }

      for (i = 0; i < n3; i++) {
         grad[i] = sg[i] = sh[i] = sr[i] = 0.0;
      }

      for (i = 0; i < 3 * n3; i++) {
         gg[i] = 0.0;
      }

      if (context_PxQ >= 0) {
         for (i = 0; i < super_NxN->size; i++) {
            sumdeijdn[i] = 0.0;
         }
         for (i = 0; i < super_Nx3N->size; i++) {
            sumdeijdg[i] = 0.0;
         }
         for (i = 0; i < super_Nx3N->size; i++) {
            sumdeijdh[i] = 0.0;
         }
         for (i = 0; i < super_Nx3N->size; i++) {
            dreff[i] = 0.0;
         }
      }

      if (gb2_debug) {
         printf("arrays initialized\n");
         fflush(stdout);
      }
#else

#ifndef ROW_CYCLIC
#pragma omp for
#endif

      for (i = 0; i < n; i++) {

#if defined(ROW_CYCLIC) && defined(OPENMP)

         if (!myroc(i, blocksize, numthreads, threadnum))
            continue;
#endif

         in = n * i;
         in3 = 3 * in;

         sumdeijda[i] = 0.0;

         for (j = 0; j < n; j++) {
            sumdeijdn[in + j] = 0.0;
         }

         for (j = 0; j < n3; j++) {
            sumdeijdg[in3 + j] = sumdeijdh[in3 + j] = dreff[in3 + j] = 0.0;
         }
      }

#endif

      /*
       * Allocate and initialize the iexw array used for skipping excluded
       * atoms.  Note that because of the manner in which iexw is used, it
       * is necessary to initialize it before only the first iteration of
       * the following loop.
       */

      iexw = ivector(-1, n);
      for (i = -1; i < n; i++) {
         iexw[i] = -1;
      }

      /*
       * Loop over all atoms i.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       *
       * For OpenMP the threads are arranged in a row cyclic manner so each
       * iteration of the following loops is assigned to a different thread.
       *
       * For ScaLAPACK the processes are arranged in a block cyclic manner
       * so each iteration of the following loops is assigned to all process
       * columns of a given process row.  The mapping of iterations onto
       * process rows is accomplished such that the submatrices maintained
       * by the processes of the row are indexed by the global matrix row
       * index that is specified by the loop index.  Before a matrix element
       * is updated, a check is performed to select the correct process
       * column as well.
       *
       * The Hessian matrix g has dimensionality 3Nx3N.  The gradient
       * vector f has dimensionality 3Nx1.  The sumdeijdn array has
       * dimensionality NxN.  The sumdeijdg, sumdeijdh and dreff arrays
       * have dimensionality Nx3N.  To facilitate selection of the correct
       * process row and column for all of these matrices, the block
       * dimensions have been chosen to match the dimensionality, that is,
       * a dimensionality of 3N requires a block dimension of 3*blocksize,
       * whereas a dimensionality of N requires a block dimension of
       * blocksize.  This approach permits the use of the bs variable
       * instead of the bs3 variable in the calls to the myroc function
       * below.
       *
       * It is not necessary to synchronize OpenMP threads following
       * this loop because the next loop uses identical loop index to
       * thread mapping.  Hence there will be no race condition associated
       * with updates of the sumdeijda, sumdeijdn, sumdeijdg, sumdeijdh,
       * gradient and Hessian arrays.  In particlar, this (i,j) loop
       * updates i elements of the sumdeijda and f gradient vectors, and
       * (i,i) and (i,j) elements of the sumdeijdg, sumdeijdh and Hessian
       * matrices, whereas the next (j,i) loop updates j elements of the
       * sumdeijda and f gradient vectors,and (j,j) and (j,i) elements
       * of the sumdeijdg, sumdeijdh and Hessian matrices, so each loop
       * updates the same array elements.
       *
       * For ScaLAPACK this loop updates i and j elements of the sumdeijda
       * and (grad) gradient vectors, so there is no need to update these
       * vectors in the next (j,i) loop.  For ScaLAPACK this loop updates
       * matrices in the same manner as OpenMP, so the following (j,i) loop
       * is necessary.
       *
       * For MPI and single-threaded execution this loop updates all elements
       * of all vectors and matrices, so the following (j,i) loop is not
       * necessary.
       */

      for (i = 0; i < n; i++) {

#if defined(OPENMP) || defined(SCALAPACK)

         if (!myroc(i, blocksize, numthreads, threadnum))
            continue;

#endif

         in = n * i;

         qi = q[i];
         ri = reff[i];

         /*
          * If atom i is not frozen, compute the "diagonal" energy that
          * is a function of only the effective radius Ri but not of the
          * interatomic distance Dij.  Compute also the contribution of
          * the diagonal energy term to the sum by which the derivatives
          * of Ri will be multiplied.
          */

         if (!frozen[i]) {
            expmkf = exp(-KSCALE * (*kappa) * ri) / (*diel_ext);
            dielfac = 1.0 - expmkf;
            qi2h = 0.5 * qi * qi;
            qid2h = qi2h * dielfac;

            vdwterm = 0.0;
            vdwdenom = 1.0;

            /*
             * For ScaLAPACK, update scalars such as epol by only one
             * process row (such as row 0) of a given process column
             * in this i-loop, but update scalars by all process rows
             * in the i,j and j,i nested loops.
             */

#ifdef SCALAPACK

            if (myrow == 0) {
               epol += -qid2h / ri;
               evdwnp += vdwterm;
            }
#else
            epol += -qid2h / ri;
            evdwnp += vdwterm;
#endif

            vdwterm *= -3.0 * vdwdenom * ri * ri;
            r_dij =
                qid2h - KSCALE * (*kappa) * qi2h * expmkf * ri + vdwterm;

            vdwterm *= -4.0 * vdwdenom * ri * ri;
            temp18 =
                ((*kappa) * KSCALE * expmkf * ri *
                 (2.0 + (*kappa) * KSCALE * ri)
                 - 2.0 * dielfac) * ri;
            p_dij = qi2h * temp18 + vdwterm;

            /*
             * Each of the ScaLAPACK process rows performs the following
             * update to sumdeijda[i] because it occurs in the i-loop not
             * the j-loop.  Therefore, disallow all but one update.
             */

#ifdef SCALAPACK

            if (myrow == 0) {
               sumdeijda[i] += r_dij;
            }
#else
            sumdeijda[i] += r_dij;
#endif

            /*
             * For ScaLAPACK, the diagonal elements of matrices
             * such as sumdeijdn are not updated; rather, an
             * associated vector such as sn is updated instead.
             * This approach is necessary because the diagonal
             * elements (i,i) and (j,j) receive contributions
             * from all (i,j) atom pairs; however, because the
             * pair list is distributed across process columns,
             * a process that calculates the contribution from
             * atom pair (i,j) will be able to access element
             * (i,j) of the matrix but that process may not be
             * able to access element (i,i) or (j,j).  Hence
             * a vector such as sn exists in each process for
             * the purpose of updating the diagonal elements of
             * these matrices.
             *
             * Also, as discussed above for sumdeijda[i], disallow
             * all but one update of sn[i] in the i-loop.
             */

#ifdef SCALAPACK

            if (myrow == 0) {
               sn[i] += p_dij;
            }
#else
            sumdeijdn[in + i] += p_dij;
#endif
         }

         /*
          * Skip the pair calculations if there are no atoms j on the
          * pair list of atom i.
          */

         npairs = upears[i];
         if (npairs <= 0)
            continue;

         i3 = 3 * i;
         in3 = 3 * in;

         xi = x[i3];
         yi = x[i3 + 1];
         zi = x[i3 + 2];

         iaci = prm->Ntypes * (prm->Iac[i] - 1);

         ri2 = ri * ri;
         riinv = 1.0 / ri;

         /* Initialize the first derivative accumulators. */

         dx = dy = dz = 0.0;

         /*
          * Expand the excluded atom list into the iexw array by storing i
          * at array address j.
          */

         for (j = 0; j < prm->Iblo[i]; j++) {
            iexw[IexclAt[i][j] - 1] = i;
         }

         /* Select atom j from the pair list.  Non-graceful error handling. */

         for (m = lpears[i]; m < lpears[i] + npairs; m++) {

            if (pearlist[i] == NULL) {
               printf("NULL pair list entry in egb2 loop 3, taskid = %d\n",
                      mytaskid);
               fflush(stdout);
            }
            j = pearlist[i][m];

            j3 = 3 * j;
            jn = n * j;
            jn3 = 3 * jn;

            /* Continue computing the non-diagonal energy term. */

            xij = xi - x[j3];
            yij = yi - x[j3 + 1];
            zij = zi - x[j3 + 2];
            r2 = xij * xij + yij * yij + zij * zij;

            /*
             * Because index j is retrieved from the pairlist array it is
             * not constrained to a particular range of values; therefore,
             * the threads that have loaded the reff array must be
             * synchronized prior to the use of reff below.
             */

            r2inv = 1.0 / r2;
            qiqj = qi * q[j];
            rj = reff[j];
            rj2 = rj * rj;
            rjinv = 1.0 / rj;
            rirjinv = riinv * rjinv;
            rb2 = ri * rj;
            efac = exp(-r2 / (4.0 * rb2));
            fgbi = 1.0 / sqrt(r2 + rb2 * efac);
            fgbk = -(*kappa) * KSCALE / fgbi;

            /*
             * Calculate the "non-diagonal" energy term, i.e., the term that is a
             * function of the effective radius Ri and the interatomic distance Dij.
             */

            expmkf = exp(fgbk) / (*diel_ext);
            dielfac = 1.0 - expmkf;

            epol += -qiqj * dielfac * fgbi;

            temp4 = fgbi * fgbi * fgbi;
            temp6 = qiqj * temp4 * (dielfac + fgbk * expmkf);

            /* The following term is required by equation 1. */

            a_dij = temp6 * (1.0 - 0.25 * efac);

            /*
             * Prepare to compute the 1st derivatives of the effective radii.
             * The following code performs the sums indicated by equations 2 and 3.
             * Compute the contribution of the non-diagonal energy term to the sum
             * by which the derivatives of Ri and Rj will be multiplied. For OpenMP
             * sum the contribution for the derivative of Ri only; otherwise, sum
             * the contribution for the derivatives of Ri and Rj.
             */

            temp5 = 0.5 * efac * temp6 * (rb2 + 0.25 * r2);
            b_dij = ri * temp5;
            c_dij = rj * temp5;

            sumdeijda[i] += b_dij;

#if !defined(OPENMP) && !defined(SCALAPACK)

            sumdeijda[j] += c_dij;
#endif


            /*
             * Calculate the first and second derivatives of the interatomic
             * distance Dij with respect to the cartesian coordinates of atoms
             * i and j.  The results are placed into five arrays:
             *
             *   di[] for the first derivatives with respect to atom i
             *   dj[] for the first derivatives with respect to atom j
             *   d2ii[] for the second derivatives with respect to atom i
             *   d2jj[] for the second derivatives with respect to atom j
             *   d2ij[] for the second derivatives with respect to atoms i and j
             *
             * Some useful symmetry obtains.  The d2ii and d2ij arrays
             * are symmetric in that their lower and upper triangles are
             * the transposes of one another.  Also, d2jj equals d2ii, and
             * d2ij is the negative of d2ii.
             *
             * Note: the additional factor of the interatomic distance Dij that
             * ought to be included in the denominator of the following equations
             * is provided by the absence of Dij in equations 15 and 16, as well as
             * by one less power of Dij in the denominator of the equations for the
             * derivatives of the Van der Waals and Coulombic energy terms.
             */

            di[0] = xij;
            di[1] = yij;
            di[2] = zij;

            dj[0] = -xij;
            dj[1] = -yij;
            dj[2] = -zij;

            /* Load the upper triangle of d2ii. */

            d2ii[0][0] = 1.0 - xij * xij * r2inv;
            d2ii[0][1] = -xij * yij * r2inv;
            d2ii[0][2] = -xij * zij * r2inv;

            d2ii[1][1] = 1.0 - yij * yij * r2inv;
            d2ii[1][2] = -yij * zij * r2inv;

            d2ii[2][2] = 1.0 - zij * zij * r2inv;

            /* Finish loading the rest of all of the matrices. */

            for (k = 0; k < 3; k++) {

               /* Load the upper triangles of d2jj and d2ij. */

               for (l = k; l < 3; l++) {

#if !defined(OPENMP) && !defined(SCALAPACK)

                  d2jj[k][l] = d2ii[k][l];
#endif
                  d2ij[k][l] = -d2ii[k][l];
               }

               /* Load the symmetric elements of d2ii, d2jj and d2ij. */

               for (l = k + 1; l < 3; l++) {
                  d2ii[l][k] = d2ii[k][l];

#if !defined(OPENMP) && !defined(SCALAPACK)

                  d2jj[l][k] = d2jj[k][l];
#endif
                  d2ij[l][k] = d2ij[k][l];
               }
            }

            /*
             * Compute the Van der Waals and Coulomb energy here because
             * if egb2 is called nbond2 is not called.  Remember that the
             * derivatives of Dij with respect to x, y and z supply one more
             * power of 1/Dij but that power of 1/Dij is included in the
             * equations below, instead of with the derivatives of Dij.
             *
             * Compute the Van der Waals and Coulombic energies for only
             * those pairs which are not on the excluded atom list.  Any
             * pair on the excluded atom list will have atom i stored at
             * address j of the iexw array.  It is not necessary to reset
             * the elements of the iexw array to -1 between successive
             * iterations in i because an i,j pair is uniquely identified
             * by atom i stored at array address j.  Thus for example, the
             * i+1,j pair would be stored at the same address as the i,j
             * pair but after the i,j pair were used.
             *
             * Initialize df here instead of within the following if statement
             * so that df is initialized unconditionally.
             */

            df = 0.0;

            if (iexw[j] != i) {

               /* Initialize dg. */

               dg = 0.0;

               rinv = 1. / sqrt(r2);
               r4inv = r2inv * r2inv;

               /*  gas-phase Coulomb energy:  */

               eel = qiqj * rinv;

               elec += eel;

               df -= eel;
               dg += 2. * eel;

               /* Lennard-Jones energy:   */

               ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
               if (ic >= 0) {
                  r6inv = r2inv * r2inv * r2inv;
                  f6 = prm->Cn2[ic] * r6inv;
                  f12 = prm->Cn1[ic] * r6inv * r6inv;

                  evdw += f12 - f6;

                  df += 6. * f6 - 12. * f12;
                  dg += 156. * f12 - 42. * f6;
               }

               /*
                * The df term includes Dij in the denominator, and the dg term
                * includes Dij^2 in the denominator; hence, terms such as d2ii
                * and di to not include Dij.
                */

               df *= r2inv;
               dg *= r4inv;

               /*
                * Calculate the second derivatives of Dij for the Van der Waals
                * and Coulombic interactions.  This calculation uses df and dg
                * for these interactions only.
                *
                * When in single-threaded mode sum to the Hessian the ii, ij, ji
                * and jj terms.  For OpenMP sum to the ii and ij terms only.
                * For ScaLAPACK sum to the ji term only and sum to the gg vector
                * instead of to the ii term.
                */

               for (k = 0; k < 3; k++) {
                  for (l = 0; l < 3; l++) {
                     ijterm = df * d2ij[k][l] + dg * di[k] * dj[l];

#ifdef SCALAPACK

                     gg[3 * (i3 + k) + l] +=
                         df * d2ii[k][l] + dg * di[k] * di[l];

                     ptr = ptr2f(g, super_3Nx3N, j3 + l, i3 + k, &info);
                     if (ptr != NULL)
                        *ptr += ijterm;
#else
                     g[n3 * (i3 + k) + i3 + l] +=
                         df * d2ii[k][l] + dg * di[k] * di[l];
                     g[n3 * (i3 + k) + j3 + l] += ijterm;
#endif

#if !defined(OPENMP) && !defined(SCALAPACK)

                     g[n3 * (j3 + l) + i3 + k] += ijterm;
                     g[n3 * (j3 + k) + j3 + l] +=
                         df * d2jj[k][l] + dg * dj[k] * dj[l];
#endif
                  }
               }
            }

            /*
             * Add to df the contribution to the gradient vector from
             * equation 1 so that it will be added to the gradient vector
             * along with the contributions from the Van der Waals and
             * Coulombic terms.
             */

            df += a_dij;

            /*
             * Complete the first derivatives of Dij for the Van der Waals
             * and Coulombic interactions, as well as for the component of
             * the first derivatives that is contributed by equation 1.
             *
             * Sum into the first derivative accumulators (dx, dy, dz)
             * the derivatives that are computed relative to the Cartesian
             * coordiantes of atom i.
             *
             * For ScaLAPACK sum into the gradient vector grad the derivative
             * that is computed relative to the Cartesian coordinates of atom
             * j.  For all other cases except OpenMP sum this derivative into
             * the gradient vector f.
             */

            dedx = df * xij;
            dedy = df * yij;
            dedz = df * zij;

            dx += dedx;
            dy += dedy;
            dz += dedz;

#ifdef SCALAPACK

            grad[j3] -= dedx;
            grad[j3 + 1] -= dedy;
            grad[j3 + 2] -= dedz;

#elif !defined(OPENMP)

            f[j3] -= dedx;
            f[j3 + 1] -= dedy;
            f[j3 + 2] -= dedz;
#endif

            /* The following code sets up the sums for equations 5-16. */

            temp5 = 0.5 * temp6;
            temp6 = fgbi * fgbi * temp5;
            temp7 =
                qiqj * (*kappa) * (*kappa) * KSCALE * KSCALE * temp4 *
                expmkf;
            temp8 = 1 - 0.25 * efac;
            temp9 = 0.25 * temp7 - 1.5 * temp6;
            temp10 = 0.25 * r2;
            temp11 = rb2 + temp10;
            temp12 = ri2 * rj2;
            temp13 = ri * temp11;
            temp14 = rj * temp11;
            temp15 = temp10 * rjinv;
            temp16 = temp10 * riinv;
            temp17 = temp9 * temp13 * efac;
            temp18 = temp9 * temp14 * efac;
            temp20 = temp8 * temp9;
            temp21 = temp20 * temp13;
            temp22 = temp20 * temp14;
            temp23 = 2.0 * efac;

            d_dij = temp23 * (temp21 - 0.25 * temp15 * temp5);

            e_dij = temp23 * (temp22 - 0.25 * temp16 * temp5);

            f_dij = 2.0 * temp8 * temp5;

            g_dij = 4.0 * temp8 * temp20
                + (2.0 * temp8 + efac * temp10 * rirjinv) * r2inv * temp5;

            h_dij =
                efac * (temp10 * (temp13 * rjinv - ri2) * temp5 +
                        temp13 * temp17);

            i_dij =
                efac * ((temp12 + temp16 * temp13) * temp5 +
                        temp13 * temp18);

            j_dij =
                temp23 * ((ri - rjinv * temp11) * 0.25 * temp5 + temp21);

            l_dij =
                efac * ((temp12 + temp15 * temp14) * temp5 +
                        temp14 * temp17);

            m_dij =
                efac * (temp10 * (temp14 * riinv - rj2) * temp5 +
                        temp14 * temp18);

            n_dij =
                temp23 * ((rj - riinv * temp11) * 0.25 * temp5 + temp22);


            /*
             * Evaluate equations 15 and 16 and sum the results to the Hessian
             * matrix.  These equations involve the derivatives d(Dij)/d(Xl),
             * d(Dij)/d(Xm) and d2(Dij)/d(Xl)d(Xm), which derivatives exist
             * only for values of k==i, k==j, l==i and l==j.
             *
             * When in single-threaded mode sum to the Hessian the ii, ij, ji
             * and jj terms.  For OpenMP sum to the ii and ij terms only.
             * For ScaLAPACK sum to the ji term only and sum to the gg vector
             * instead of to the ii term.
             *
             * The df term includes Dij in the denominator, and the dg term
             * includes Dij^2 in the denominator; hence, terms such as d2ii
             * and di to not include Dij.
             */

            df = f_dij;
            dg = g_dij;

            for (k = 0; k < 3; k++) {
               for (l = 0; l < 3; l++) {
                  ijterm = df * d2ij[k][l] + dg * di[k] * dj[l];

#ifdef SCALAPACK

                  gg[3 * (i3 + k) + l] +=
                      df * d2ii[k][l] + dg * di[k] * di[l];

                  ptr = ptr2f(g, super_3Nx3N, j3 + l, i3 + k, &info);
                  if (ptr != NULL)
                     *ptr += ijterm;
#else
                  g[n3 * (i3 + k) + i3 + l] +=
                      df * d2ii[k][l] + dg * di[k] * di[l];
                  g[n3 * (i3 + k) + j3 + l] += ijterm;
#endif

#if !defined(OPENMP) && !defined(SCALAPACK)

                  g[n3 * (j3 + l) + i3 + k] += ijterm;
                  g[n3 * (j3 + k) + j3 + l] +=
                      df * d2jj[k][l] + dg * dj[k] * dj[l];
#endif
               }
            }

            /*
             * Setup the inner differentiations indicated by equations 9 and 12.
             * The column index selects the first derivative to be calculated,
             * i.e., d(Ri) or d(Rj).  The row index is i because we are
             * setting up the outer multiplication by d(Ri).
             */

#ifdef SCALAPACK

            sn[i] += h_dij;

            ptr = ptr2f(sumdeijdn, super_NxN, j, i, &info);
            if (ptr != NULL)
               *ptr += i_dij;
#else
            sumdeijdn[in + i] += h_dij;
            sumdeijdn[in + j] += l_dij;
#endif

            /*
             * Setup the inner differentiations indicated by equations 10 and 11.
             * The column index selects the first derivative to be calculated,
             * i.e., d(Ri) or d(Rj).  The row index is j because we are
             * setting up the outer multiplication by d(Rj).  When in single-
             * threaded mode sum both terms.  When in multi-threaded mode sum
             * neither term.
             */

#if !defined(OPENMP) && !defined(SCALAPACK)

            sumdeijdn[jn + i] += i_dij;
            sumdeijdn[jn + j] += m_dij;
#endif

            /*
             * Perform the inner differentiation indicated by equations 5, 6, 7 and 8.
             * The row or first index is selected by i because we are setting up the
             * outer differentiation by the derivatives of Ri.  The column or second
             * index is selected by i or j because the Dij term affects the derivatives
             * of atoms i and j.  Both the sumdeijdg and sumdeijdh arrays are used
             * because equations 5 and 6 perform inner differentiation with respect to
             * Xl, whereas equations 7 and 8 perform inner differentiation with respect
             * to Xm.  Distinguishing these two cases is accomplished by transposing
             * the Hessian addresses, as will be done in the final two calls to the
             * dgemm() function below.
             *
             * When in single-threaded mode sum into the sumdeigdg and sumdeijdh arrays
             * all terms.  When in multi-threaded mode sum only those terms that comprise
             * d(Ri) times d(Dij).
             *
             * For ScaLAPACK, update matrices such as sumdeijdg
             * by the correct process column.  The correct process
             * row has already been selected by the loop index i.
             * Note that because sumdeijdg has dimensionality Nx3N
             * both blocksize (bs) and 3*blocksize (bs3) are used in
             * calculation of the local submatrix address.
             */

            for (k = 0; k < 3; k++) {

#ifdef SCALAPACK

               sg[i3 + k] += d_dij * di[k];
               sh[i3 + k] += j_dij * di[k];

               ptr = ptr2f(sumdeijdg, super_Nx3N, j, i3 + k, &info);
               if (ptr != NULL)
                  *ptr += e_dij * di[k];

               ptr = ptr2f(sumdeijdh, super_Nx3N, j, i3 + k, &info);
               if (ptr != NULL)
                  *ptr += n_dij * di[k];

#else
               sumdeijdg[in3 + i3 + k] += d_dij * di[k];
               sumdeijdg[in3 + j3 + k] += d_dij * dj[k];
               sumdeijdh[in3 + i3 + k] += j_dij * di[k];
               sumdeijdh[in3 + j3 + k] += j_dij * dj[k];
#endif

#if !defined(OPENMP) && !defined(SCALAPACK)

               sumdeijdg[jn3 + i3 + k] += e_dij * di[k];
               sumdeijdg[jn3 + j3 + k] += e_dij * dj[k];
               sumdeijdh[jn3 + i3 + k] += n_dij * di[k];
               sumdeijdh[jn3 + j3 + k] += n_dij * dj[k];
#endif
            }
         }

         /* Sum the first derivative accumulators into the gradient vector. */

#ifdef SCALAPACK

         grad[i3] += dx;
         grad[i3 + 1] += dy;
         grad[i3 + 2] += dz;

#else

         f[i3] += dx;
         f[i3 + 1] += dy;
         f[i3 + 2] += dz;

#endif

      }

      if (gb2_debug) {
         printf("end loop 1\n");
         fflush(stdout);
      }

      /*
       * Do not execute the following nest of loops for single-threads or MPI.
       * For OpenMP update the j elements of the gradient vector f and of the
       * sumdeijda array.  For OpenMP and ScaLAPACK update the j,i and j,j
       * elements of the sumdeijdg, sumdeijdh and Hessian matrices.
       */

#if defined(OPENMP) || defined(SCALAPACK)

      /*
       * Initialize the iexw array used for skipping excluded atoms.
       * Note that because of the manner in which iexw is used, it
       * is necessary to initialize it before only the first iteration
       * of the following loop.
       */

      for (i = -1; i < n; i++) {
         iexw[i] = -1;
      }

      /*
       * Loop over all atoms j.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       */

      for (j = 0; j < n; j++) {

         if (!myroc(j, blocksize, numthreads, threadnum))
            continue;

         /*
          * Skip the pair calculations if there are no atoms i on the
          * pair list of atom j.
          */

         npairs = lpears[j];
         if (npairs <= 0)
            continue;

         j3 = 3 * j;
         jn = n * j;
         jn3 = n * j3;

         xj = x[j3];
         yj = x[j3 + 1];
         zj = x[j3 + 2];

         qj = q[j];
         rj = reff[j];
         rj2 = rj * rj;
         rjinv = 1.0 / rj;

         /* Initialize the first derivative accumulators. */

         dx = dy = dz = 0.0;

         /*
          * Expand the excluded list into the iexw array by storing j
          * at array address i.
          */

         for (i = 0; i < Jblo[j]; i++) {
            iexw[JexclAt[j][i] - 1] = j;
         }

         /* Select atom i from the pair list.  Non-graceful error handling. */

         for (m = 0; m < npairs; m++) {

            if (pearlist[j] == NULL) {
               printf("NULL pair list entry in egb2 loop 4, taskid = %d\n",
                      mytaskid);
               fflush(stdout);
            }
            i = pearlist[j][m];

            in = n * i;
            i3 = 3 * i;
            in3 = n * i3;

            xij = x[i3] - xj;
            yij = x[i3 + 1] - yj;
            zij = x[i3 + 2] - zj;
            r2 = xij * xij + yij * yij + zij * zij;
            r2inv = 1.0 / r2;

            iaci = prm->Ntypes * (prm->Iac[i] - 1);

            /*
             * Because index i is retrieved from the pairlist array it is
             * not constrained to a particular range of values; therefore,
             * the threads that have loaded the reff array must be
             * synchronized prior to the use of reff below.
             */

            qiqj = q[i] * qj;
            ri = reff[i];
            ri2 = ri * ri;
            riinv = 1.0 / ri;
            rirjinv = riinv * rjinv;
            rb2 = ri * rj;
            efac = exp(-r2 / (4.0 * rb2));
            fgbi = 1.0 / sqrt(r2 + rb2 * efac);
            fgbk = -(*kappa) * KSCALE / fgbi;

            /*
             * Calculate the "non-diagonal" energy term, i.e., the term that is a
             * function of the effective radius Ri and the interatomic distance Dij.
             */

            expmkf = exp(fgbk) / (*diel_ext);
            dielfac = 1.0 - expmkf;

            temp4 = fgbi * fgbi * fgbi;
            temp6 = qiqj * temp4 * (dielfac + fgbk * expmkf);

            /* The following term is required by equation 1. */

            a_dij = temp6 * (1.0 - 0.25 * efac);

            /*
             * Prepare to compute the 1st derivatives of the effective radii.
             * The following code performs the sums indicated by equations 2 and 3.
             * Compute the contribution of the non-diagonal energy term to the sum
             * by which the derivative of Rj will be multiplied.
             */

            temp5 = 0.5 * efac * temp6 * (rb2 + 0.25 * r2);
            c_dij = rj * temp5;

            /* Update the j element of the sumdeijda array. */

            sumdeijda[j] += c_dij;

            /*
             * Calculate the first and second derivatives of the interatomic
             * distance Dij with respect to the cartesian coordinates of atoms
             * i and j.  The results are placed into two arrays:
             *
             *   di[] for the first derivatives with respect to atom i
             *   dj[] for the first derivatives with respect to atom j
             *   d2jj[] for the second derivatives with respect to atom j
             *   d2ij[] for the second derivatives with respect to atoms i and j
             *
             * Some useful symmetry obtains.  The d2ii and d2ij arrays
             * are symmetric in that their lower and upper triangles are
             * the transposes of one another.  Also, d2jj equals d2ii, and
             * d2ij is the negative of d2ii.
             *
             * Note: the additional factor of the interatomic distance Dij that
             * ought to be included in the denominator of the following equations
             * is provided by the absence of Dij in equations 15 and 16, as well as
             * by one less power of Dij in the denominator of the equations for the
             * derivatives of the Van der Waals and Coulombic energy terms.
             */

            di[0] = xij;
            di[1] = yij;
            di[2] = zij;

            dj[0] = -xij;
            dj[1] = -yij;
            dj[2] = -zij;

            /* Load the upper triangle of d2jj. */

            d2jj[0][0] = 1.0 - xij * xij * r2inv;
            d2jj[0][1] = -xij * yij * r2inv;
            d2jj[0][2] = -xij * zij * r2inv;

            d2jj[1][1] = 1.0 - yij * yij * r2inv;
            d2jj[1][2] = -yij * zij * r2inv;

            d2jj[2][2] = 1.0 - zij * zij * r2inv;

            /* Finish loading the rest of all of the matrices. */

            for (k = 0; k < 3; k++) {

               /* Load the upper triangle of d2ij. */

               for (l = k; l < 3; l++) {
                  d2ij[k][l] = -d2jj[k][l];
               }

               /* Load the symmetric elements of d2jj and d2ij. */

               for (l = k + 1; l < 3; l++) {
                  d2jj[l][k] = d2jj[k][l];
                  d2ij[l][k] = d2ij[k][l];
               }
            }

            /*
             * Compute the Van der Waals and Coulomb energy here because
             * if egb2 is called nbond2 is not called.  Remember that the
             * derivatives of Dij with respect to x, y and z supply one more
             * power of 1/Dij but that power of 1/Dij is included in the
             * equations below, instead of with the derivatives of Dij.
             *
             * Compute the Van der Waals and Coulombic energies for only
             * those pairs which are not on the excluded atom list.  Any
             * pair on the excluded atom list will have atom i stored at
             * address j of the iexw array.  It is not necessary to reset
             * the elements of the iexw array to -1 between successive
             * iterations in i because an i,j pair is uniquely identified
             * by atom i stored at array address j.  Thus for example, the
             * i,j+1 pair would be stored at the same address as the i,j
             * pair but after the i,j pair were used.
             *
             * Initialize df here instead of within the following if statement
             * so that df is initialized unconditionally.
             */

            df = 0.0;

            if (iexw[i] != j) {

               /* Initialize dg. */

               dg = 0.0;

               rinv = 1. / sqrt(r2);
               r4inv = r2inv * r2inv;

               /*  gas-phase Coulomb energy:  */

               eel = qiqj * rinv;

               df -= eel;
               dg += 2. * eel;

               /* Lennard-Jones energy:   */

               ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
               if (ic >= 0) {
                  r6inv = r2inv * r2inv * r2inv;
                  f6 = prm->Cn2[ic] * r6inv;
                  f12 = prm->Cn1[ic] * r6inv * r6inv;

                  df += 6. * f6 - 12. * f12;
                  dg += 156. * f12 - 42. * f6;
               }

               /*
                * The df term includes Dij in the denominator, and the dg term
                * includes Dij^2 in the denominator; hence, terms such as d2ii
                * and di to not include Dij.
                */

               df *= r2inv;
               dg *= r4inv;

               /*
                * Calculate the second derivatives of Dij for the Van der Waals
                * and Coulombic interactions.  This calculation uses df and dg
                * for these interactions only.
                *
                * For OpenMP sum the ji and jj terms into the Hessian.
                * For ScaLAPACK sum the ij term into the Hessian but
                * sum the jj term into the gg vector.
                */

               for (k = 0; k < 3; k++) {
                  for (l = 0; l < 3; l++) {

#ifdef SCALAPACK

                     ptr = ptr2f(g, super_3Nx3N, i3 + k, j3 + l, &info);
                     if (ptr != NULL)
                        *ptr += df * d2ij[k][l] + dg * di[k] * dj[l];

                     gg[3 * (j3 + k) + l] +=
                         df * d2jj[k][l] + dg * dj[k] * dj[l];

#else

                     g[n3 * (j3 + l) + i3 + k] +=
                         df * d2ij[k][l] + dg * di[k] * dj[l];
                     g[n3 * (j3 + k) + j3 + l] +=
                         df * d2jj[k][l] + dg * dj[k] * dj[l];
#endif
                  }
               }
            }

            /*
             * Add to df the contribution to the gradient vector from
             * equation 1 so that it will be added to the gradient vector
             * along with the contributions from the Van der Waals and
             * Coulombic terms.
             */

            df += a_dij;

            /*
             * Complete the first derivatives of Dij for the Van der Waals
             * and Coulombic interactions, as well as for the component of
             * the first derivatives that is contributed by equation 1.
             *
             * Sum the derivatives into the derivative accumulators.
             */

            dedx = df * xij;
            dedy = df * yij;
            dedz = df * zij;

            dx += dedx;
            dy += dedy;
            dz += dedz;

            /* The following code sets up the sums for equations 5-16. */

            temp5 = 0.5 * temp6;
            temp6 = fgbi * fgbi * temp5;
            temp7 =
                qiqj * (*kappa) * (*kappa) * KSCALE * KSCALE * temp4 *
                expmkf;
            temp8 = 1 - 0.25 * efac;
            temp9 = 0.25 * temp7 - 1.5 * temp6;
            temp10 = 0.25 * r2;
            temp11 = rb2 + temp10;
            temp12 = ri2 * rj2;
            temp13 = ri * temp11;
            temp14 = rj * temp11;
            temp15 = temp10 * rjinv;
            temp16 = temp10 * riinv;
            temp17 = temp9 * temp13 * efac;
            temp18 = temp9 * temp14 * efac;
            temp20 = temp8 * temp9;
            temp21 = temp20 * temp13;
            temp22 = temp20 * temp14;
            temp23 = 2.0 * efac;

            d_dij = temp23 * (temp21 - 0.25 * temp15 * temp5);

            e_dij = temp23 * (temp22 - 0.25 * temp16 * temp5);

            f_dij = 2.0 * temp8 * temp5;

            g_dij = 4.0 * temp8 * temp20
                + (2.0 * temp8 + efac * temp10 * rirjinv) * r2inv * temp5;

            h_dij =
                efac * (temp10 * (temp13 * rjinv - ri2) * temp5 +
                        temp13 * temp17);

            i_dij =
                efac * ((temp12 + temp16 * temp13) * temp5 +
                        temp13 * temp18);

            j_dij =
                temp23 * ((ri - rjinv * temp11) * 0.25 * temp5 + temp21);

            l_dij =
                efac * ((temp12 + temp15 * temp14) * temp5 +
                        temp14 * temp17);

            m_dij =
                efac * (temp10 * (temp14 * riinv - rj2) * temp5 +
                        temp14 * temp18);

            n_dij =
                temp23 * ((rj - riinv * temp11) * 0.25 * temp5 + temp22);


            /*
             * Evaluate equations 15 and 16 and sum the results to the Hessian
             * matrix.  These equations involve the derivatives d(Dij)/d(Xl),
             * d(Dij)/d(Xm) and d2(Dij)/d(Xl)d(Xm), which derivatives exist
             * only for values of k==i, k==j, l==i and l==j.
             *
             * For OpenMP sum the ji and jj terms into the Hessian.
             * For ScaLAPACK sum the ij term into the Hessian but
             * sum the jj term into the gg vector.
             *
             * The df term includes Dij in the denominator, and the dg term
             * includes Dij^2 in the denominator; hence, terms such as d2ii
             * and di to not include Dij.
             */

            df = f_dij;
            dg = g_dij;

            for (k = 0; k < 3; k++) {
               for (l = 0; l < 3; l++) {

#ifdef SCALAPACK

                  ptr = ptr2f(g, super_3Nx3N, i3 + k, j3 + l, &info);
                  if (ptr != NULL)
                     *ptr += df * d2ij[k][l] + dg * di[k] * dj[l];

                  gg[3 * (j3 + k) + l] +=
                      df * d2jj[k][l] + dg * dj[k] * dj[l];

#else

                  g[n3 * (j3 + l) + i3 + k] +=
                      df * d2ij[k][l] + dg * di[k] * dj[l];
                  g[n3 * (j3 + k) + j3 + l] +=
                      df * d2jj[k][l] + dg * dj[k] * dj[l];
#endif
               }
            }

            /*
             * Setup the inner differentiations indicated by equations 10 and 11.
             * The column index selects the first derivative to be calculated,
             * i.e., d(Rj).  The row index is j because we are setting up the
             * outer multiplication by d(Rj).
             */

#ifdef SCALAPACK

            ptr = ptr2f(sumdeijdn, super_NxN, i, j, &info);
            if (ptr != NULL)
               *ptr += l_dij;

            sn[j] += m_dij;
#else
            sumdeijdn[jn + i] += i_dij;
            sumdeijdn[jn + j] += m_dij;
#endif

            /*
             * Perform the inner differentiation indicated by equations 5, 6, 7 and 8.
             * The row or first index is selected by j because we are setting up the
             * outer differentiation by the derivatives of Rj.  The column or second
             * index is selected by i or j because the Dij term affects the derivatives
             * of atoms i and j.  Both the sumdeijdg and sumdeijdh arrays are used
             * because equations 5 and 6 perform inner differentiation with respect to
             * Xl, whereas equations 7 and 8 perform inner differentiation with respect
             * to Xm.  Distinguishing these two cases is accomplished by transposing
             * the Hessian addresses, as will be done in the final two calls to the
             * dgemm() function below.
             *
             * Sum those terms that comprise d(Rj) times d(Dij).
             */

            for (k = 0; k < 3; k++) {

#ifdef SCALAPACK

               sg[j3 + k] += e_dij * dj[k];
               sh[j3 + k] += n_dij * dj[k];

               ptr = ptr2f(sumdeijdg, super_Nx3N, i, j3 + k, &info);
               if (ptr != NULL)
                  *ptr += d_dij * dj[k];

               ptr = ptr2f(sumdeijdh, super_Nx3N, i, j3 + k, &info);
               if (ptr != NULL)
                  *ptr += j_dij * dj[k];
#else
               sumdeijdg[jn3 + i3 + k] += e_dij * di[k];
               sumdeijdg[jn3 + j3 + k] += e_dij * dj[k];
               sumdeijdh[jn3 + i3 + k] += n_dij * di[k];
               sumdeijdh[jn3 + j3 + k] += n_dij * dj[k];
#endif
            }
         }

         /*
          * For OpenMP but not for ScaLAPACK, sum the first derivative
          * accumulators into the gradient vector f.
          */

#ifdef OPENMP

         f[j3] -= dx;
         f[j3 + 1] -= dy;
         f[j3 + 2] -= dz;

#endif

      }

      if (gb2_debug) {
         printf("end loop 2\n");
         fflush(stdout);
      }
#endif                          /* if defined(OPENMP) || defined(SCALAPACK) */


      /*
       * There is no pragma directive associated with the above 'for' loop
       * so explicit thread synchronization is performed here.
       */

#pragma omp barrier

      /*
       * The SCALAPACK synchronization is accomplised via reduction
       * of the sumdeijda array.  Each element of this array should
       * be non-zero for one process only.
       *
       * No reduction is needed for MPI because no parallelization is performed.
       *
       * Actually, for both OpenMP and SCALAPACK it is possible to defer
       * synchronization until after the following (i,j) nest of loops
       * because the mapping of loop indices to threads or tasks is
       * identical to that of the two previous nests of loops, and because
       * only the i elements of the sumdeijda array are read in the following
       * nest of loops.  However, synchronization is possible at this point,
       * and may serve to better balance the workload in the next two nests
       * of loops.
       */

#ifdef SCALAPACK

      ierror = MPI_Allreduce(sumdeijda, reductarr, n,
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (ierror != MPI_SUCCESS) {
         printf
             ("Error in egb sumdeijda reduction, error = %d  mytaskid = %d\n",
              ierror, mytaskid);
      }

      if (gb2_debug) {
         printf("MPI_Allreduce sumdeijda\n");
         fflush(stdout);
      }

      for (i = 0; i < n; i++) {
         sumdeijda[i] = reductarr[i];
      }

#endif

      /* Free the iexw array within this potentially parallel region of code. */

      free_ivector(iexw, -1, n);

      /*
       * Compute the first and second derivatives of the effective radius.
       * Update the gradient vector using the first derivatives.  Complete
       * equations 2, 3 and 4.  Begin equations 13, 14 and 17 which update
       * the Hessian matrix.
       *
       * For OpenMP this loop updates the i elements of the gradient
       * vector f as well as the (i,i) and (i,j) elements of the dreff and
       * Hessian matrices.
       *
       * For ScaLAPACK this loop updates the i and j elements of the gradient
       * vector grad as well as the (i,i) and (i,j) elements of the dreff and
       * Hessian matrices.
       *
       * For MPI and single-threaded execution this loop updates all elements
       * of all vectors and matrices, so the following (j,i) loop is not
       * necessary.
       *
       * Loop over all atoms i.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       *
       * It is not necessary to synchronize OpenMP threads following
       * this loop because the next loop uses identical loop index to
       * thread mapping.  Hence there will be no race condition associated
       * with updates of the dreff, gradient and Hessian arrays.
       */

      for (i = 0; i < n; i++) {

#if defined(OPENMP) || defined(SCALAPACK)

         if (!myroc(i, blocksize, numthreads, threadnum))
            continue;

#endif

         /*
          * Don't calculate derivatives of the effective radius of atom i
          * if atom i is frozen or if there are no pair atoms j associated
          * with atom i.
          */

         npairs = lpears[i] + upears[i];
         if (frozen[i] || (npairs <= 0))
            continue;

         /* Initialize the first derivative accumulators. */

         dx = dy = dz = 0.0;

         i3 = 3 * i;
         in = n * i;
         in3 = 3 * in;

         /*
          * The following captures all terms that are to be multiplied by the
          * first derivative of the effective radius (equations 2, 3 and 4).
          */

         sumda = sumdeijda[i];

         if (gb > 1) {

            ri = rborn[i] - gboffset;
            thi =
                tanh((gbalpha[i] - gbbeta[i] * psi[i] +
                      gbgamma[i] * psi[i] * psi[i]) * psi[i]);
            sumda *=
                (gbalpha[i] - 2.0 * gbbeta[i] * psi[i] +
                 3.0 * gbgamma[i] * psi[i] * psi[i])
                * (1.0 - thi * thi) * ri / rborn[i];
         }

         /*
          * The following captures all terms that are to be multiplied by the
          * second derivative of the effective radius (equations 13, 14 and 17).
          */

         sumda = sumdeijda[i];

         /* Prepare to compute the derivatives of the effective radius Ri. */

         xi = x[i3];
         yi = x[i3 + 1];
         zi = x[i3 + 2];
         ri = rborn[i] - gboffset;
         ri1i = 1. / ri;

         /* Select atom j from the pair list.  Non-graceful error handling. */

         for (m = 0; m < npairs; m++) {

            if (pearlist[i] == NULL) {
               printf("NULL pair list entry in egb2 loop 5, taskid = %d\n",
                      mytaskid);
               fflush(stdout);
            }
            j = pearlist[i][m];

            j3 = 3 * j;
            jn3 = n * j3;

            xij = xi - x[j3];
            yij = yi - x[j3 + 1];
            zij = zi - x[j3 + 2];
            r2 = xij * xij + yij * yij + zij * zij;

            /* Ignore the ij atom pair if their separation exceeds the GB cutoff. */

            if (r2 > rgbmaxpsmax2)
               continue;

            dij1i = 1.0 / sqrt(r2);
            dij2i = dij1i * dij1i;
            dij = r2 * dij1i;
            sj = fs[j] * (rborn[j] - gboffset);
            sj2 = sj * sj;

            /*
             * The following are the numerator of the first derivatives of the
             * equations from the Appendix of Schaefer and Froemmel and of the
             * Taylor series expansion for d>>s by Andreas Svrcek-Seiler.  Smooth
             * rgbmax idea is from Andreas Svrcek-Seiler and Alexey Onufriev.  The
             * derivatives are with respect to Dij.  The complete derivative is
             * formed by multiplying the numerator by -Ri*Ri.  The factor of Ri*Ri
             * has been moved to equations 5-18.  The negation is deferred until later.
             * When the chain rule is used to form the first derivatives of the
             * effective radius with respect to the cartesian coordinates, an
             * additional factor of Dij appears in the denominator.  That factor
             * is included in the following expressions.
             */

            if (dij > rgbmax + sj)
               continue;

            if (dij > rgbmax - sj) {

               temp1 = 1. / (dij - sj);
               dij3i = dij1i * dij2i;

               datmp = 0.125 * dij3i * ((r2 + sj2) *
                                        (temp1 * temp1 - rgbmax2i) -
                                        2.0 * log(rgbmax * temp1));

            } else if (dij > 4.0 * sj) {

               tmpsd = sj2 * dij2i;
               dumbo =
                   TE + tmpsd * (TF +
                                 tmpsd * (TG +
                                          tmpsd * (TH + tmpsd * THH)));
               datmp = tmpsd * sj * dij2i * dij2i * dumbo;

            } else if (dij > ri + sj) {

               temp1 = 1. / (r2 - sj2);
               datmp = temp1 * sj * (-0.5 * dij2i + temp1)
                   + 0.25 * dij1i * dij2i * log((dij - sj) / (dij + sj));

            } else if (dij > fabs(ri - sj)) {

               temp1 = 1. / (dij + sj);
               dij3i = dij2i * dij1i;
               datmp =
                   -0.25 * (-0.5 * (r2 - ri * ri + sj2) * dij3i * ri1i *
                            ri1i + dij1i * temp1 * (temp1 - dij1i)
                            - dij3i * log(ri * temp1));

            } else if (ri < sj) {

               temp1 = 1. / (r2 - sj2);
               datmp =
                   -0.5 * (sj * dij2i * temp1 - 2. * sj * temp1 * temp1 -
                           0.5 * dij2i * dij1i * log((sj - dij) /
                                                     (sj + dij)));

            } else
               continue;

            /*
             * Negate the first derivative of Ri with respect to the
             * interatomic distance in order to account for the
             * minus sign of -Ri*Ri.  The factor of Ri*Ri is omitted
             * because it was included in equations 5-18.  These
             * derivatives contain the appropriate factor of Dij
             * in the denominator because it is supplied by datmp.
             */

            df = -datmp;

            /*
             * The following are the numerator of the second derivatives of the
             * equations from the Appendix of Schaefer and Froemmel and of the
             * Taylor series expansion for d>>s by Andreas Svrcek-Seiler.  Smooth
             * rgbmax idea is from Andreas Svrcek-Seiler and Alexey Onufriev.  The
             * derivatives are with respect to Dij.
             */

            if (dij > rgbmax + sj)
               continue;

            if (dij > rgbmax - sj) {

               temp1 = 1. / (dij - sj);
               temp2 = temp1 * temp1;
               dij3i = dij1i * dij2i;

               da2tmp =
                   -0.25 * (dij3i * (sj2 * rgbmax2i - 1.0) +
                            dij1i * temp2 +
                            2.0 * (log(rgbmax * temp1) * dij3i +
                                   temp1 * (dij2i - temp2)));

            } else if (dij > 4.0 * sj) {

               tmpsd = sj2 * dij2i;
               dumbo =
                   TI + tmpsd * (TJ +
                                 tmpsd * (TK +
                                          tmpsd * (TL + tmpsd * TLL)));
               da2tmp = tmpsd * sj * dij2i * dij2i * dumbo;


            } else if (dij > ri + sj) {

               temp1 = 2. * dij * sj;
               temp2 = (dij - sj) * (dij + sj);
               dij3i = dij2i * dij1i;

               da2tmp = (0.5 * log((dij - sj) / (dij + sj)) +
                         dij * sj * ((r2 + sj2) * (r2 - sj2) +
                                     temp1 * temp1) / (temp2 * temp2 *
                                                       temp2)) * dij3i;

            } else if (dij > fabs(ri - sj)) {

               temp1 = dij + sj;
               temp2 = dij * temp1;
               ri2 = ri * ri;

               da2tmp =
                   ((2. * ri2 * log(ri / temp1) +
                     sj2) * temp1 * temp1 * temp1 +
                    ri2 * sj * (2. * dij * dij - dij * sj -
                                sj2)) / (4. * ri2 * temp2 * temp2 * temp2);

            } else if (ri < sj) {

               temp1 = sj2 - r2;
               temp2 = 2. * sj * dij;
               dij3i = dij2i * dij1i;

               da2tmp = (0.5 * log((sj - dij) / (sj + dij)) + sj * dij *
                         (temp1 * (sj2 + r2) - temp2 * temp2) /
                         (temp1 * temp1 * temp1)) * dij3i;

            } else
               continue;

            /* 
             * Calculate the second derivative of the effective radius
             * with respect to the interatomic distance.  A factor of
             * the Ri*Ri has been omitted from the second derivative
             * because it was included in equations 5-18.
             *
             * When the chain rule is used to form the second derivatives of the
             * effective radius with respect to the cartesian coordinates, an
             * additional factor of Dij*Dij appears in the denominator.  That
             * factor is included in the following expression.
             */

            dg = da2tmp * dij2i;

            /*
             * Calculate the first and second derivatives of the interatomic
             * distance Dij with respect to the cartesian coordinates of atoms
             * i and j.  The results are placed into five arrays:
             *
             *   di[] for the first derivatives with respect to atom i
             *   dj[] for the first derivatives with respect to atom j
             *   d2ii[] for the second derivatives with respect to atom i
             *   d2jj[] for the second derivatives with respect to atom j
             *   d2ij[] for the second derivatives with respect to atoms i and j
             *
             * Some useful symmetry obtains.  The d2ii and d2ij arrays
             * are symmetric in that their lower and upper triangles are
             * the transposes of one another.  Also, d2jj equals d2ii, and
             * d2ij is the negative of d2ii.
             *
             *
             * Note: the additional factor of the interatomic distance Dij that
             * ought to be included in the denominator of the following equations
             * is included instead in the first and second derivatives of the
             * effective radius with respect to Dij (df and dg) that were calculated
             * above.
             */

            di[0] = xij;
            di[1] = yij;
            di[2] = zij;

            dj[0] = -xij;
            dj[1] = -yij;
            dj[2] = -zij;

            /* Load the upper triangle of d2ii. */

            r2inv = 1.0 / r2;

            d2ii[0][0] = 1.0 - xij * xij * r2inv;
            d2ii[0][1] = -xij * yij * r2inv;
            d2ii[0][2] = -xij * zij * r2inv;

            d2ii[1][1] = 1.0 - yij * yij * r2inv;
            d2ii[1][2] = -yij * zij * r2inv;

            d2ii[2][2] = 1.0 - zij * zij * r2inv;

            /* Finish loading the rest of all of the matrices. */

            for (k = 0; k < 3; k++) {

               /* Load the upper triangles of d2jj and d2ij. */

               for (l = k; l < 3; l++) {

#if !defined(OPENMP) && !defined(SCALAPACK)

                  d2jj[k][l] = d2ii[k][l];
#endif
                  d2ij[k][l] = -d2ii[k][l];
               }

               /* Load the symmetric elements of d2ii, d2jj and d2ij. */

               for (l = k + 1; l < 3; l++) {
                  d2ii[l][k] = d2ii[k][l];

#if !defined(OPENMP) && !defined(SCALAPACK)

                  d2jj[l][k] = d2jj[k][l];
#endif
                  d2ij[l][k] = d2ij[k][l];
               }
            }

            /*
             * Store the 1st derivatives of Ri with respect to (xj,yj,zj)
             * if SCALAPACK is not defined.  If SCALAPACK is defined, store
             * these derivatives in the next nest of loops.
             */

            dxj = df * dj[0];
            dyj = df * dj[1];
            dzj = df * dj[2];

#ifndef SCALAPACK

            dreff[in3 + j3] = dxj;
            dreff[in3 + j3 + 1] = dyj;
            dreff[in3 + j3 + 2] = dzj;
#endif

            /*
             * If OPENMP is not defined, update the gradient vector f
             * with d(Ri)/d(Xj).  If SCALAPACK is defined update the
             * gradient vector grad instead of the gradient vector f.
             */

#ifdef SCALAPACK

            grad[j3] += sumda * dxj;
            grad[j3 + 1] += sumda * dyj;
            grad[j3 + 2] += sumda * dzj;

#elif !defined(OPENMP)

            f[j3] += sumda * dxj;
            f[j3 + 1] += sumda * dyj;
            f[j3 + 2] += sumda * dzj;
#endif

            /*
             * Sum the 1st derivatives of Ri with respect to (xi,yi,zi).
             * Exploit the fact that d(Dij)/d(Xi) = -d(Dij)/d(Xj).
             */

            dx -= dxj;
            dy -= dyj;
            dz -= dzj;

            /*
             * Form the second derivatives of the effective radius with
             * respect to the cartesian coordinates using the chain rule,
             * multiply the result by sumda and sum to the Hessian matrix.
             * This step implements equations 13, 14 and 17, except that
             * terms such as d2(Ri)/d(Dij)d(Dik) are not handled here
             * but instead are handled in the first call to dgemm() below.
             *
             * When in single-threaded mode update the Hessian using the ii, ij
             * ji and jjterms.  For OpenMP update the Hessian using the ii and ij
             * terms only.  For ScaLAPACK update the Hessian using the ji term
             * only and update the i-element of the gg vector.
             */

            for (k = 0; k < 3; k++) {
               for (l = 0; l < 3; l++) {
                  ijterm = sumda * (df * d2ij[k][l] + dg * di[k] * dj[l]);

#ifdef SCALAPACK

                  gg[3 * (i3 + k) + l] +=
                      sumda * (df * d2ii[k][l] + dg * di[k] * di[l]);

                  ptr = ptr2f(g, super_3Nx3N, j3 + l, i3 + k, &info);
                  if (ptr != NULL)
                     *ptr += ijterm;
#else
                  g[n3 * (i3 + k) + i3 + l] +=
                      sumda * (df * d2ii[k][l] + dg * di[k] * di[l]);
                  g[n3 * (i3 + k) + j3 + l] += ijterm;
#endif

#if !defined(OPENMP) && !defined(SCALAPACK)

                  g[n3 * (j3 + l) + i3 + k] += ijterm;
                  g[n3 * (j3 + k) + j3 + l] +=
                      sumda * (df * d2jj[k][l] + dg * dj[k] * dj[l]);
#endif
               }
            }
         }

         /*
          * Store the sum of first derivatives of Ri with respect to (xi,yi,zi)
          * in the diagonal elements of the dreff array.  These derivatives
          * include the appropriate factor of Dij but not a factor of Ri*Ri.
          * They include the minus sign of -Ri*Ri that was supplied in df.
          */

#ifdef SCALAPACK

         sr[i3] = dx;
         sr[i3 + 1] = dy;
         sr[i3 + 2] = dz;

#else

         dreff[in3 + i3] = dx;
         dreff[in3 + i3 + 1] = dy;
         dreff[in3 + i3 + 2] = dz;

#endif

         /* Update the gradient vector with d(Ri)/d(Xi). */

#ifdef SCALAPACK

         grad[i3] += sumda * dx;
         grad[i3 + 1] += sumda * dy;
         grad[i3 + 2] += sumda * dz;

#else

         f[i3] += sumda * dx;
         f[i3 + 1] += sumda * dy;
         f[i3 + 2] += sumda * dz;

#endif

      }

      if (gb2_debug) {
         printf("end loop 3\n");
         fflush(stdout);
      }

      /*
       * This nest of loops updates the (j,i) and (j,j) elements of the
       * Hessian matrix for ScaLAPACK or OpenMP.  For OpenMP update
       * the j elements of the gradient vector f.  For ScaLAPACK do
       * not update the gradient vector grad.  This nest of loop is
       * not necessary for single-threaded execution or for MPI.
       * 
       * The control flow in this nest of loops must be identical to 
       * that of the preceeding nest of loops.  Perform the minimum
       * calculation necessary to reconstitute the conditions of
       * the preceeding nest of loops.
       *
       * Because no computation follows the following 'for' loop within
       * this parallel region, the implied barrier at the end of the
       * parallel region enforces the necessary thread synchronization
       * for OpenMP threads.
       */

#if defined(OPENMP) || defined(SCALAPACK)

      /*
       * Loop over all atoms j.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       */

      for (j = 0; j < n; j++) {

         if (!myroc(j, blocksize, numthreads, threadnum))
            continue;

         /*
          * Don't calculate derivatives with respect to atom j if there
          * are no atoms i associated with atom j.
          */

         npairs = lpears[j] + upears[j];
         if (npairs <= 0)
            continue;

         j3 = 3 * j;
         jn3 = n * j3;

         xj = x[j3];
         yj = x[j3 + 1];
         zj = x[j3 + 2];

         sj = fs[j] * (rborn[j] - gboffset);
         sj2 = sj * sj;

         /* Initialize the first derivative accumulators. */

         dx = dy = dz = 0.0;

         /* Select atom i from the pair list.  Non-graceful error handling. */

         for (m = 0; m < npairs; m++) {

            if (pearlist[j] == NULL) {
               printf("NULL pair list entry in egb2 loop 6, taskid = %d\n",
                      mytaskid);
               fflush(stdout);
            }
            i = pearlist[j][m];

            /*
             * Don't calculate derivatives of the effective radius of atom i
             * if atom i is frozen.
             */

            if (frozen[i])
               continue;

            i3 = 3 * i;
            in = (n * i);

            xij = x[i3] - xj;
            yij = x[i3 + 1] - yj;
            zij = x[i3 + 2] - zj;
            r2 = xij * xij + yij * yij + zij * zij;

            /* Ignore the ij atom pair if their separation exceeds the GB cutoff. */

            if (r2 > rgbmaxpsmax2)
               continue;

            dij1i = 1.0 / sqrt(r2);
            dij2i = dij1i * dij1i;
            dij = r2 * dij1i;

            ri = rborn[i] - gboffset;
            ri1i = 1. / ri;

            /*
             * The following are the numerator of the first derivatives of the
             * equations from the Appendix of Schaefer and Froemmel and of the
             * Taylor series expansion for d>>s by Andreas Svrcek-Seiler.  Smooth
             * rgbmax idea is from Andreas Svrcek-Seiler and Alexey Onufriev.  The
             * derivatives are with respect to Dij.  The complete derivative is
             * formed by multiplying the numerator by -Ri*Ri.  The factor of Ri*Ri
             * has been moved to equations 5-18.  The negation is deferred until later.
             * When the chain rule is used to form the first derivatives of the
             * effective radius with respect to the cartesian coordinates, an
             * additional factor of Dij appears in the denominator.  That factor
             * is included in the following expressions.
             */

            if (dij > rgbmax + sj)
               continue;

            if (dij > rgbmax - sj) {

               temp1 = 1. / (dij - sj);
               dij3i = dij1i * dij2i;

               datmp = 0.125 * dij3i * ((r2 + sj2) *
                                        (temp1 * temp1 - rgbmax2i) -
                                        2.0 * log(rgbmax * temp1));

            } else if (dij > 4.0 * sj) {

               tmpsd = sj2 * dij2i;
               dumbo =
                   TE + tmpsd * (TF +
                                 tmpsd * (TG +
                                          tmpsd * (TH + tmpsd * THH)));
               datmp = tmpsd * sj * dij2i * dij2i * dumbo;

            } else if (dij > ri + sj) {

               temp1 = 1. / (r2 - sj2);
               datmp = temp1 * sj * (-0.5 * dij2i + temp1)
                   + 0.25 * dij1i * dij2i * log((dij - sj) / (dij + sj));

            } else if (dij > fabs(ri - sj)) {

               temp1 = 1. / (dij + sj);
               dij3i = dij2i * dij1i;
               datmp =
                   -0.25 * (-0.5 * (r2 - ri * ri + sj2) * dij3i * ri1i *
                            ri1i + dij1i * temp1 * (temp1 - dij1i)
                            - dij3i * log(ri * temp1));

            } else if (ri < sj) {

               temp1 = 1. / (r2 - sj2);
               datmp =
                   -0.5 * (sj * dij2i * temp1 - 2. * sj * temp1 * temp1 -
                           0.5 * dij2i * dij1i * log((sj - dij) /
                                                     (sj + dij)));

            } else
               continue;

            /*
             * Negate the first derivative of Ri with respect to the
             * interatomic distance in order to account for the
             * minus sign of -Ri*Ri.  The factor of Ri*Ri is omitted
             * because it was included in equations 5-18.  These
             * derivatives contain the appropriate factor of Dij
             * in the denominator because it is supplied by datmp.
             */

            df = -datmp;

            /*
             * The following are the numerator of the second derivatives of the
             * equations from the Appendix of Schaefer and Froemmel and of the
             * Taylor series expansion for d>>s by Andreas Svrcek-Seiler.  Smooth
             * rgbmax idea is from Andreas Svrcek-Seiler and Alexey Onufriev.  The
             * derivatives are with respect to Dij.
             */

            if (dij > rgbmax + sj)
               continue;

            if (dij > rgbmax - sj) {

               temp1 = 1. / (dij - sj);
               temp2 = temp1 * temp1;
               dij3i = dij1i * dij2i;

               da2tmp =
                   -0.25 * (dij3i * (sj2 * rgbmax2i - 1.0) +
                            dij1i * temp2 +
                            2.0 * (log(rgbmax * temp1) * dij3i +
                                   temp1 * (dij2i - temp2)));

            } else if (dij > 4.0 * sj) {

               tmpsd = sj2 * dij2i;
               dumbo =
                   TI + tmpsd * (TJ +
                                 tmpsd * (TK +
                                          tmpsd * (TL + tmpsd * TLL)));
               da2tmp = tmpsd * sj * dij2i * dij2i * dumbo;


            } else if (dij > ri + sj) {

               temp1 = 2. * dij * sj;
               temp2 = (dij - sj) * (dij + sj);
               dij3i = dij2i * dij1i;

               da2tmp = (0.5 * log((dij - sj) / (dij + sj)) +
                         dij * sj * ((r2 + sj2) * (r2 - sj2) +
                                     temp1 * temp1) / (temp2 * temp2 *
                                                       temp2)) * dij3i;

            } else if (dij > fabs(ri - sj)) {

               temp1 = dij + sj;
               temp2 = dij * temp1;
               ri2 = ri * ri;

               da2tmp =
                   ((2. * ri2 * log(ri / temp1) +
                     sj2) * temp1 * temp1 * temp1 +
                    ri2 * sj * (2. * dij * dij - dij * sj -
                                sj2)) / (4. * ri2 * temp2 * temp2 * temp2);

            } else if (ri < sj) {

               temp1 = sj2 - r2;
               temp2 = 2. * sj * dij;
               dij3i = dij2i * dij1i;

               da2tmp = (0.5 * log((sj - dij) / (sj + dij)) + sj * dij *
                         (temp1 * (sj2 + r2) - temp2 * temp2) /
                         (temp1 * temp1 * temp1)) * dij3i;

            } else
               continue;

            /* 
             * Calculate the second derivative of the effective radius
             * with respect to the interatomic distance.  A factor of
             * the Ri*Ri has been omitted from the second derivative
             * because it was included in equations 5-18.
             *
             * When the chain rule is used to form the second derivatives of the
             * effective radius with respect to the cartesian coordinates, an
             * additional factor of Dij*Dij appears in the denominator.  That
             * factor is included in the following expression.
             */

            dg = da2tmp * dij2i;

            /*
             * Calculate the first and second derivatives of the interatomic
             * distance Dij with respect to the cartesian coordinates of atoms
             * i and j.  The results are placed into four arrays:
             *
             *   di[] for the first derivatives with respect to atom i
             *   dj[] for the first derivatives with respect to atom j
             *   d2jj[] for the second derivatives with respect to atom j
             *   d2ij[] for the second derivatives with respect to atoms i and j
             *
             * Note: the additional factor of the interatomic distance Dij that
             * ought to be included in the denominator of the following equations
             * is provided by the absence of Dij in equations 15 and 16, as well as
             * by one less power of Dij in the denominator of the equations for the
             * derivatives of the Van der Waals and Coulombic energy terms.
             */

            di[0] = xij;
            di[1] = yij;
            di[2] = zij;

            dj[0] = -xij;
            dj[1] = -yij;
            dj[2] = -zij;

            /* Load the upper triangle of d2jj. */

            r2inv = 1.0 / r2;

            d2jj[0][0] = 1.0 - xij * xij * r2inv;
            d2jj[0][1] = -xij * yij * r2inv;
            d2jj[0][2] = -xij * zij * r2inv;

            d2jj[1][1] = 1.0 - yij * yij * r2inv;
            d2jj[1][2] = -yij * zij * r2inv;

            d2jj[2][2] = 1.0 - zij * zij * r2inv;

            /* Finish loading the rest of all of the matrices. */

            for (k = 0; k < 3; k++) {

               /* Load the upper triangle of d2ij. */

               for (l = k; l < 3; l++) {
                  d2ij[k][l] = -d2jj[k][l];
               }

               /* Load the symmetric elements of d2ij and d2jj. */

               for (l = k + 1; l < 3; l++) {
                  d2ij[l][k] = d2ij[k][l];
                  d2jj[l][k] = d2jj[k][l];
               }
            }

            /*
             * The following captures all terms that are to be multiplied by the
             * second derivative of the effective radius (equations 13, 14 and 17).
             *
             * Because index i is retrieved from the pairlist array it is
             * not constrained to a particular range of values; therefore,
             * the threads that have loaded the sumdeijda array have been
             * synchronized above prior to the use of sumdeijda below.
             */

            sumda = sumdeijda[i];

            /*
             * Form the second derivatives of the effective radius with
             * respect to the cartesian coordinates using the chain rule,
             * multiply the result by sumda and sum to the Hessian matrix.
             * This step implements equations 13, 14 and 17, except that
             * terms such as d2(Ri)/d(Dij)d(Dik) are not handled here
             * but instead are handled in the first call to dgemm() below.
             *
             * For OpenMP update the Hessian using the ji and jj terms only.
             * For ScaLAPACK update the Hessian using the ij term only, and
             * update the j element of the ff vector.
             */

            for (k = 0; k < 3; k++) {
               for (l = 0; l < 3; l++) {

#ifdef SCALAPACK

                  ptr = ptr2f(g, super_3Nx3N, i3 + k, j3 + l, &info);
                  if (ptr != NULL)
                     *ptr +=
                         sumda * (df * d2ij[k][l] + dg * di[k] * dj[l]);

                  gg[3 * (j3 + k) + l] +=
                      sumda * (df * d2jj[k][l] + dg * dj[k] * dj[l]);

#else

                  g[n3 * (j3 + l) + i3 + k] +=
                      sumda * (df * d2ij[k][l] + dg * di[k] * dj[l]);
                  g[n3 * (j3 + k) + j3 + l] +=
                      sumda * (df * d2jj[k][l] + dg * dj[k] * dj[l]);
#endif
               }
            }

            /*
             * The following captures all terms that are to be multiplied by the
             * first derivative of the effective radius (equations 2, 3 and 4).
             */

            if (gb > 1) {

               ri = rborn[i] - gboffset;
               thi =
                   tanh((gbalpha[i] - gbbeta[i] * psi[i] +
                         gbgamma[i] * psi[i] * psi[i]) * psi[i]);
               sumda *=
                   (gbalpha[i] - 2.0 * gbbeta[i] * psi[i] +
                    3.0 * gbgamma[i] * psi[i] * psi[i])
                   * (1.0 - thi * thi) * ri / rborn[i];
            }

            /* Store the 1st derivatives of Ri with respect to (xj,yj,zj). */

            dxj = df * dj[0];
            dyj = df * dj[1];
            dzj = df * dj[2];

#ifdef SCALAPACK

            ptr = ptr2f(dreff, super_Nx3N, i, j3 + 0, &info);
            if (ptr != NULL)
               *ptr = dxj;

            ptr = ptr2f(dreff, super_Nx3N, i, j3 + 1, &info);
            if (ptr != NULL)
               *ptr = dyj;

            ptr = ptr2f(dreff, super_Nx3N, i, j3 + 2, &info);
            if (ptr != NULL)
               *ptr = dzj;
#endif

            /* Update the first derivative accumulators with d(Ri)/d(Xj). */

            dx += sumda * dxj;
            dy += sumda * dyj;
            dz += sumda * dzj;
         }

         /*
          * Update the gradient vector f from the first derivative accumulators,
          * but only for OpenMP, not for ScaLAPACK.
          */

#ifdef OPENMP

         f[j3] += dx;
         f[j3 + 1] += dy;
         f[j3 + 2] += dz;

#endif

      }

      if (gb2_debug) {
         printf("end loop 4\n");
         fflush(stdout);
      }
#endif                          /* if defined(OPENMP) || defined(SCALAPACK) */

   }

#ifdef SCALAPACK

   /*
    * Add the gradient vector grad to the gradient vector f.
    * This step is straightforward for a 1x1 process grid
    * but requires multiple steps for a larger process grid.
    */

   if (gridim == 1) {
      if (context_PxQ >= 0) {
         for (i = 0; i < m3; i++) {
            f[i] += grad[i];
         }

         if (gb2_debug) {
            printf("f += grad\n");
            fflush(stdout);
         }
      }
   } else {

      /*
       * Perform a reduction over the gradient vector grad.
       * Leave the result in reductarr which exists in all processes.
       */

      ierror = MPI_Allreduce(grad, reductarr, m3,
                             MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (ierror != MPI_SUCCESS) {
         printf("Error in egb grad reduction, error = %d  mytaskid = %d\n",
                ierror, mytaskid);
      }

      if (gb2_debug) {
         printf("MPI_Allreduce grad\n");
         fflush(stdout);
      }

      /* Gather f from context_PxQ into grad in context_1x1. */

      pdgemr2d_(&m3, &one,
                f, &one, &one, desc_3Nx1,
                grad, &one, &one, desc_1x1, &context_Nx1);

      if (gb2_debug) {
         printf("pgemr2d f\n");
         fflush(stdout);
      }

      /* Add reductarr to grad. */

      if (context_1x1 >= 0) {
         for (i = 0; i < m3; i++) {
            grad[i] += reductarr[i];
         }

         if (gb2_debug) {
            printf("f += reductarr\n");
            fflush(stdout);
         }
      }

      /* Distribute grad from context_1x1 into f in context_PxQ. */

      pdgemr2d_(&m3, &one,
                grad, &one, &one, desc_1x1,
                f, &one, &one, desc_3Nx1, &context_Nx1);

      if (gb2_debug) {
         printf("pgemr2d grad\n");
         fflush(stdout);
      }
   }

   /*
    * The diagonal elements of the Hessian, sumdeijdg, sumdeijgh,
    * sumdeijdn and dreff matrices are stored in the gg, sg, sh,
    * sn and sr vectors.  Perform a reduction across each vector
    * and write the result into corresponding matrix.
    *
    * Here is the code for the Hessian.
    */

   ierror = MPI_Allreduce(gg, reductarr, 3 * m3,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (ierror != MPI_SUCCESS) {
      printf("Error in egb gg reduction, error = %d  mytaskid = %d\n",
             ierror, mytaskid);
   }

   if (gb2_debug) {
      printf("MPI_Allreduce gg\n");
      fflush(stdout);
   }

   for (i = 0; i < n; i++) {
      i3 = 3 * i;
      for (k = 0; k < 3; k++) {
         for (l = 0; l < 3; l++) {
            ptr = ptr2f(g, super_3Nx3N, i3 + k, i3 + l, &info);
            if (ptr != NULL)
               *ptr += reductarr[3 * (i3 + k) + l];
         }
      }
   }

   /* Here is the code for sumdeijdg. */

   ierror = MPI_Allreduce(sg, reductarr, m3,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (ierror != MPI_SUCCESS) {
      printf("Error in egb sg reduction, error = %d  mytaskid = %d\n",
             ierror, mytaskid);
   }

   if (gb2_debug) {
      printf("MPI_Allreduce sg\n");
      fflush(stdout);
   }

   for (i = 0; i < n; i++) {
      i3 = 3 * i;
      for (k = 0; k < 3; k++) {
         ptr = ptr2f(sumdeijdg, super_Nx3N, i, i3 + k, &info);
         if (ptr != NULL)
            *ptr += reductarr[i3 + k];
      }
   }

   /* Here is the code for sumdeijdh. */

   ierror = MPI_Allreduce(sh, reductarr, m3,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (ierror != MPI_SUCCESS) {
      printf("Error in egb sh reduction, error = %d  mytaskid = %d\n",
             ierror, mytaskid);
   }

   if (gb2_debug) {
      printf("MPI_Allreduce sh\n");
      fflush(stdout);
   }

   for (i = 0; i < n; i++) {
      i3 = 3 * i;
      for (k = 0; k < 3; k++) {
         ptr = ptr2f(sumdeijdh, super_Nx3N, i, i3 + k, &info);
         if (ptr != NULL)
            *ptr += reductarr[i3 + k];
      }
   }

   /* Here is the code for dreff. */

   ierror = MPI_Allreduce(sr, reductarr, m3,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (ierror != MPI_SUCCESS) {
      printf("Error in egb sr reduction, error = %d  mytaskid = %d\n",
             ierror, mytaskid);
   }

   if (gb2_debug) {
      printf("MPI_Allreduce sr\n");
      fflush(stdout);
   }

   for (i = 0; i < n; i++) {
      i3 = 3 * i;
      for (k = 0; k < 3; k++) {
         ptr = ptr2f(dreff, super_Nx3N, i, i3 + k, &info);
         if (ptr != NULL)
            *ptr += reductarr[i3 + k];
      }
   }

   /* Here is the code for sumdeijdn. */

   ierror = MPI_Allreduce(sn, reductarr, n,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (ierror != MPI_SUCCESS) {
      printf("Error in egb sn reduction, error = %d  mytaskid = %d\n",
             ierror, mytaskid);
   }

   if (gb2_debug) {
      printf("MPI_Allreduce sn\n");
      fflush(stdout);
   }

   for (i = 0; i < n; i++) {
      ptr = ptr2f(sumdeijdn, super_NxN, i, i, &info);
      if (ptr != NULL)
         *ptr += reductarr[i];
   }

   /*
    * Add to the diagonal elements of sumdeijdn the d2(Ri)/d(Dij)d(Dik)
    * terms weighted by Ri and sumdeijda[i].
    */

   for (i = 0; i < n; i++) {
      ptr = ptr2f(sumdeijdn, super_NxN, i, i, &info);
      if (ptr != NULL) {
         *ptr += 2.0 * sumdeijda[i] * reff[i];
      }
   }

   t2 = seconds();
   *tgb2Other += t2 - t1;
   if (gb2_debug) {
      printf("\nmatrix preparation time = %10.2f  n = %d  m3 = %d\n\n",
             t2 - t1, n, m3);
      fflush(stdout);
   }
   t1 = t2;

   /*
    * Perform the inner multiplications of equations 9, 10, 11, 12
    * and 18 in a manner analogous to the multiplications of
    * of equations 2, 3 and 4.  Also, multiply each element of each
    * row of the dreff array by the d2(Ri)/d(Dij)d(Dik) terms.
    *
    * Next perform the outer multiplication for equations 5, 6, 9, 10,
    * 11 and 12.
    *
    * Finally, perform the outer multiplication for equations 7 and 8.
    *
    * Some ScaLAPACK functions appear to quit unexpectedly
    * for large matrices on a 1x1 process grid, e.g. the pdgemm_
    * function in newton.c, so bypass the ScaLAPACK pdgemm_ and
    * use the LAPACK dgemm_ function instead.
    *
    * The correct test is (gridim == 1), not (nprow == 1 && npcol == 1)
    * since processes that aren't on the 1x1 grid have nprow == npcol == -1,
    * which would direct control to pdgemr2d_ (below) that would hang because
    * it would not be called from all processes, specifically not from
    * the process that is on the 1x1 grid and has nprow == npcol == 1.
    */

   if (gridim == 1) {
      if (context_PxQ >= 0) {

         transa = 'T';
         transb = 'N';
         dgemm_(&transa, &transb, &n, &m3, &n, &dblone,
                sumdeijdn, &n, dreff, &n, &dblone, sumdeijdg, &n);

         t2 = seconds();
         if (gb2_debug) {
            printf("dgemm1 time = %10.2f\n\n", t2 - t1);
            fflush(stdout);
         }
         *tgb2dgemm1 += t2 - t1;
         t1 = t2;

         transa = 'T';
         transb = 'N';
         dgemm_(&transa, &transb, &m3, &m3, &n, &dblone,
                dreff, &n, sumdeijdg, &n, &dblone, g, &m3);

         t2 = seconds();
         if (gb2_debug) {
            printf("dgemm2 time = %10.2f\n\n", t2 - t1);
            fflush(stdout);
         }
         *tgb2dgemm2 += t2 - t1;
         t1 = t2;

         dgemm_(&transa, &transb, &m3, &m3, &n, &dblone,
                sumdeijdh, &n, dreff, &n, &dblone, g, &m3);

         t2 = seconds();
         if (gb2_debug) {
            printf("dgemm3 time = %10.2f\n\n", t2 - t1);
            fflush(stdout);
         }
         *tgb2dgemm3 += t2 - t1;
      }
   } else {
      if (context_PxQ >= 0) {
         transa = 'N';
         transb = 'N';
         pdgemm_(&transa, &transb, &n, &m3, &n, &dblone,
                 sumdeijdn, &one, &one, desc_NxN,
                 dreff, &one, &one, desc_Nx3N, &dblone,
                 sumdeijdg, &one, &one, desc_Nx3N);

         t2 = seconds();
         if (gb2_debug) {
            printf("pdgemm1 time = %10.2f\n\n", t2 - t1);
            fflush(stdout);
         }
         *tgb2dgemm1 += t2 - t1;
         t1 = t2;

         transa = 'T';
         transb = 'N';
         pdgemm_(&transa, &transb, &m3, &m3, &n, &dblone,
                 dreff, &one, &one, desc_Nx3N,
                 sumdeijdg, &one, &one, desc_Nx3N, &dblone,
                 g, &one, &one, desc_3Nx3N);

         t2 = seconds();
         if (gb2_debug) {
            printf("pdgemm2 time = %10.2f\n\n", t2 - t1);
            fflush(stdout);
         }
         *tgb2dgemm2 += t2 - t1;
         t1 = t2;

         transa = 'T';
         transb = 'N';
         pdgemm_(&transa, &transb, &m3, &m3, &n, &dblone,
                 sumdeijdh, &one, &one, desc_Nx3N,
                 dreff, &one, &one, desc_Nx3N, &dblone,
                 g, &one, &one, desc_3Nx3N);

         t2 = seconds();
         if (gb2_debug) {
            printf("pdgemm3 time = %10.2f\n\n", t2 - t1);
            fflush(stdout);
         }
         *tgb2dgemm3 += t2 - t1;
      }
   }

#else

   /*
    * Here is the non-ScaLAPACK section of the code.  Note that the
    * sumdeijdg, sumdeijdh and dreff matrices are built in a transposed
    * manner relative to the manner in which they are built for
    * ScaLAPACK execution because the ScaLAPACK submatrices are built
    * in column major order.
    *
    * Add to the diagonal elements of sumdeijdn the d2(Ri)/d(Dij)d(Dik)
    * terms weighted by Ri and sumdeijda[i].
    */

   for (i = 0; i < n; i++) {
      in = (n * i);
      sumdeijdn[in + i] += 2.0 * sumdeijda[i] * reff[i];
   }

   t2 = seconds();
   if (gb2_debug) {
      printf("\nmatrix preparation time = %10.2f  n = %d  m3 = %d\n\n",
             t2 - t1, n, m3);
      fflush(stdout);
   }
   *tgb2Other += t2 - t1;
   t1 = t2;

   /*
    * Perform the inner multiplications of equations 9, 10, 11, 12
    * and 18 in a manner analogous to the multiplications of
    * of equations 2, 3 and 4.  Also, multiply each element of each
    * row of the dreff array by the d2(Ri)/d(Dij)d(Dik) terms.
    *
    * Next perform the outer multiplication for equations 5, 6, 9, 10,
    * 11 and 12.
    *
    * Finally, perform the outer multiplication for equations 7 and 8.
    */

   transa = 'N';
   transb = 'N';
   dgemm_(&transa, &transb, &m3, &n, &n, &dblone,
          dreff, &m3, sumdeijdn, &n, &dblone, sumdeijdg, &m3);

   t2 = seconds();
   if (gb2_debug) {
      printf("dgemm1 time = %10.2f\n\n", t2 - t1);
      fflush(stdout);
   }
   *tgb2dgemm1 += t2 - t1;
   t1 = t2;

   transa = 'N';
   transb = 'T';



   dgemm_(&transa, &transb, &m3, &m3, &n, &dblone,
          sumdeijdg, &m3, dreff, &m3, &dblone, g, &m3);


   t2 = seconds();
   if (gb2_debug) {
      printf("dgemm2 time = %10.2f\n\n", t2 - t1);
      fflush(stdout);
   }
   *tgb2dgemm2 += t2 - t1;
   t1 = t2;

   transa = 'N';
   transb = 'T';
   dgemm_(&transa, &transb, &m3, &m3, &n, &dblone,
          dreff, &m3, sumdeijdh, &m3, &dblone, g, &m3);

   t2 = seconds();
   if (gb2_debug) {
      printf("dgemm3 time = %10.2f\n\n", t2 - t1);
      fflush(stdout);
   }
   *tgb2dgemm3 += t2 - t1;

#endif

   if (gb2_debug) {
      printf("exit egb2\n");
      fflush(stdout);
   }

   /*
    * Return elec, evdw and evdwnp through the parameters eelt, enb and enp.
    * These variables are computed in parallel.
    */

   *eelt = elec;
   *enb = evdw;
   *enp = evdwnp;
   t2 = seconds();
   *tgb2Other += t2 - t1;
   *tgb2 += t2 - tgb21;
   return (epol);
}

/***********************************************************************
                            MME2()
************************************************************************/

/*
 * Compute the energy and first and second derivatives.
 *
 * Calling parameters are as follows:
 *
 * x - input: the atomic (x,y,z) coordinates
 * f - updated: the gradient vector
 * g - updated: the Hessian matrix
 * m - returned: the atomic masses
 * grad - updated: the gradient vector in context_1x1
 * descF_PxQ - input: the ScaLAPACK descriptor for vector f in context_PxQ
 * descF_1x1 - input: the ScaLAPACK descriptor for vector grad in context_1x1
 * descG_PxQ - input: the ScaLAPACK descriptor for matrix g in context_PxQ
 * context_PxQ - input: the distributed vector and matrix context for ScaLAPACK
 * context_1x1 - input: the non-distributed vector context for ScaLAPACK
 * context_Nx1 - input: the general context for ScaLAPACK
 * gridim - input: the ScaLAPACK process grid dimension (=1 for single task)
 * natom - returned: the number of atoms
 * iter - input: the iteration counter, which if negative selects the following:
 *        -1 print some energy values
 *        -3 call egb to deallocate static arrays, then deallocate grad
 *        -(any other value) normal execution
 */

REAL_T mme2(REAL_T * x, REAL_T * f, REAL_T * g, REAL_T * m, REAL_T * grad,
            INT_T * descF_PxQ, INT_T * descF_1x1, INT_T * descG_PxQ,
            INT_T * context_PxQ, INT_T * context_1x1, INT_T * context_Nx1,
            INT_T * gridim, INT_T * natom, INT_T * iter, char *name)
{
   REAL_T ebh, eba, eth, eta, eph, epa, enb, eel, enb14, eel14;
   REAL_T ecn, esurf, evdwnp, frms;
   REAL_T ene[20];
   REAL_T tmme21, t1, t2;
   int i, j, k, n, i3, threadnum, numthreads, maxthreads;
   size_t in9, ii, jj, kk;
   char atsymb;

#ifdef SCALAPACK
   int myrow, mycol, nprow, npcol, mb, nb;
   size_t locpF_PxQ, locqF_PxQ, locpG_PxQ, locqG_PxQ;
   int zero = 0, one = 1;
   REAL_T *ptr, reductarr[20];
   char ctrl;
#endif

   /* Enforce 3D (and disallow 4D) for second derivatives. */

   t1 = seconds();
   tmme21 = t1;
   dim = 3;
   n = 3 * prm->Natom;

   /* Error exit if this is a periodic system: second derivatives are not
      supported there:  */
   if( prm->IfBox ){
      fprintf( stderr, "Second derivatives are not supported for periodic systems; exiting.\n" );
      exit(1);
   }

   for (i = 0; i < prm->Natom; i++) {
      atsymb = prm->AtomNames[i][0];
      /* if (atsymb == 'H') {printf("yo dude %d %d",i,prm->Natom);} */
      /* first letter atom name */
      name[i] = atsymb;
   }

   /*
    * If the iteration count equals -3, call egb2 to deallocate the
    * static arrays, then return; otherwise, simply return.
    */

   if (*iter == -3) {
      egb2(lpairs, upairs, pairlist, lpairs, upairs, pairlist,
           x, f, g, grad, prm->Fs,
           prm->Rborn, prm->Charges, &kappa, &epsext, &enb, &eel, &esurf,
           &evdwnp, *context_Nx1, *context_1x1, *context_PxQ,
           descF_PxQ, descF_1x1, descG_PxQ, *gridim, 1);
      return (0.0);
   }

   /*
    * If the iteration count is zero or one, print the energy header
    * for task 0 only.
    */

   if ((*iter == 0 || *iter == 1) && mytaskid == 0) {
      fprintf(nabout, "      iter    Total       bad      vdW     elect"
              "   nonpolar   genBorn      frms\n");
      fflush(nabout);
   }

   /*
    * Write the checkpoint file every nchk2 iterations if the chknm
    * variable is non-NULL.
    */

   if (chknm != NULL && (*iter > 0 && *iter % nchk2 == 0)) {
      checkpoint(chknm, prm->Natom, x, *iter);
   }

   /*
    * Build the non-bonded pair list if it hasn't already been built;
    * rebuild it every nsnb iterations.  This pair list will be used for
    * first derivative calculations (for line minimization within
    * Newton-Raphson), for non-14 nonbonded interactions (i.e.,
    * when nbond2 is called instead of egb2) and for Born electrostatic
    * and nonpolar second derivative calculations in both single-threaded
    * and multi-threaded execution under OpenMP.  This pair list is global
    * and fully populated for OpenMP, and local and partially populated
    * for MPI.
    *
    * If Generalized Born surface area calculations are selected,
    * build the non-polar pair list if it hasn't already been built;
    * rebuild it every nsnp iterations.  This pair list will be used for
    * non-polar first derivative calculations (for line minimization
    * within Newton-Raphson) and for non-polar egb interactions. This
    * pair list is global and fully populated for OpenMP, and local and 
    * fully populated for MPI.  It is fully populated for MPI so that
    * experiments can be performed with loop scheduling and block sizes.
    * See the comment in the mme34 function of the eff.c file.
    *
    * If MPI is defined, build a second pair list that is to be used
    * by egb2 and nbond2.  This pair list will be fully populated for
    * each task because neither the Born second derivatives nor the
    * non-14 nonbonded second derivative calculations are parallelized.
    *
    * If SCALAPACK is defined, build a second pair list that is to be
    * used by egb2 for the calculation of second derivatives for the
    * Born electrostatic term and the Born van der Waals nonpolar term.
    * This pair list will conform to the block cyclic process grid.
    *
    * If SCALAPACK is defined and GB surface area is selected, build a
    * pair list that is to be used by egb2 for the calculation of the
    * second derivatives of the Born surface area nonpolar term.  This pair
    * list will be fully populated for each MPI process because calculation
    * of the Born surface area nonpolar term is not parallelized.
    * This pair list is always provided to egb2 but it will contain
    * all NULL entries and will not be used unless Generalized Born
    * surface area calculations are selected.
    */

   if (nb_pairs < 0 || (*iter > 0 && *iter % nsnb == 0)) {
      t2 = seconds();
      *tmme2Other += t2 - t1;
      t1 = t2;
      nb_pairs = nblist(lpairs, upairs, pairlist, x, *context_PxQ, 1, cut,
                        prm->Natom, dim, frozen);
      t2 = seconds();
      *tmme2Pair += t2 - t1;
      t1 = t2;
   }

#ifdef SCALAPACK

   if (gb) {
      if (nb_pairs2 < 0 || (*iter > 0 && *iter % nsnb == 0)) {
	t2 = seconds();
	*tmme2Other += t2 - t1;
	t1 = t2;
	nb_pairs2 =
	  nblist(lpairs2, upairs2, pairlist2, x, *context_PxQ, -1, cut,
		 prm->Natom, dim, frozen);
	t2 = seconds();
	*tmme2Pair += t2 - t1;
	t1 = t2;
      }

   } else {
      if (nb_pairs2 < 0 || (*iter > 0 && *iter % nsnb == 0)) {
	t2 = seconds();
	*tmme2Other += t2 - t1;
	t1 = t2;
	nb_pairs2 =
	  nblist(lpairs2, upairs2, pairlist2, x, *context_PxQ, 0, cut,
		 prm->Natom, dim, frozen);
	t2 = seconds();
	*tmme2Pair += t2 - t1;
	t1 = t2;
      }
   }

#elif defined(MPI)

   if (nb_pairs2 < 0 || (*iter > 0 && *iter % nsnb == 0)) {
     t2 = seconds();
     *tmme2Other += t2 - t1;
     t1 = t2;
     nb_pairs2 =
       nblist(lpairs2, upairs2, pairlist2, x, *context_PxQ, 0, cut,
	      prm->Natom, dim, frozen);
     t2 = seconds();
     *tmme2Pair += t2 - t1;
     t1 = t2;
   }

#endif

#ifdef SCALAPACK

   /*
    * Get the number of rows and columns on the row cyclic (PxQ) process grid,
    * as well as this task's row and column on the grid.
    */

   blacs_gridinfo_(context_PxQ, &nprow, &npcol, &myrow, &mycol);

   /*
    * Only processes that are active on the PxQ grid initialize
    * the gradient vector and Hessian matrix.
    */

   if (*context_PxQ >= 0) {

      /*
       * The numroc_ function is used to calculate the number of vector or
       * matrix elements that are distributed across a PxQ processor grid
       * for the gradient vector or the Hessian matrix.  Only processes
       * in grid column zero initialize the gradient vector.
       */

      mb = descF_PxQ[MB_];
      nb = descF_PxQ[NB_];

      locpF_PxQ = numroc_(&n, &mb, &myrow, &zero, &nprow);
      locqF_PxQ = numroc_(&one, &nb, &mycol, &zero, &npcol);

      if (mycol == 0) {
         for (i = 0; i < locpF_PxQ * locqF_PxQ; i++) {
            f[i] = 0.0;
         }
      }

      /* Get the block sizes for the Hessian matrix. */

      mb = descG_PxQ[MB_];
      nb = descG_PxQ[NB_];

      locpG_PxQ = numroc_(&n, &mb, &myrow, &zero, &nprow);
      locqG_PxQ = numroc_(&n, &nb, &mycol, &zero, &npcol);

      for (ii = 0; ii < locpG_PxQ * locqG_PxQ; ii++) {
         g[ii] = 0.0;
      }
   }
#else

#pragma omp parallel private (i, j, i3, ii, jj, in9, threadnum, numthreads)
   {

      /*
       * Get the thread number, number of threads and maximum number of threads
       * for multi-threaded execution.  The maxthreads variables is used only
       * for allocation of the gradient array.
       */

#ifdef OPENMP
      threadnum = omp_get_thread_num();
      numthreads = omp_get_num_threads();
      maxthreads = omp_get_max_threads();
#else
      maxthreads = 1;
#endif

      /*
       * Initialize the gradient and Hessian arrays inside of the
       * OpenMP parallel region.
       *
       * The "first touch" memory allocation strategy will locate
       * the elements of these arrays that are initialized by
       * a particular CPU in memory that is local to that CPU.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       *
       * Because no computation follows the following 'for' loop within
       * this parallel region, the implied barrier at the end of the
       * parallel region enforces the necessary thread synchronization.
       */

#ifndef ROW_CYCLIC
#pragma omp for
#endif

      for (i = 0; i < prm->Natom; i++) {

#if defined(ROW_CYCLIC) && defined(OPENMP)
         if (!myroc(i, blocksize, numthreads, threadnum))
            continue;
#endif

         i3 = 3 * i;
         ii = (size_t) i;
         in9 = (9 * prm->Natom * ii);

         f[i3] = f[i3 + 1] = f[i3 + 2] = 0.0;

         for (j = 0; j < 9 * prm->Natom; j++) {
            jj = (size_t) j;
            g[in9 + jj] = 0.0;
         }
      }
   }

#endif

   t2 = seconds();
   *tmme2Other += t2 - t1;
   t1 = t2;

   ebh = ebond2(prm->Nbonh, prm->BondHAt1, prm->BondHAt2,
                prm->BondHNum, prm->Rk, prm->Req, x, f, g,
                *context_PxQ, descF_PxQ, descG_PxQ);
   eba = ebond2(prm->Mbona, prm->BondAt1, prm->BondAt2,
                prm->BondNum, prm->Rk, prm->Req, x, f, g,
                *context_PxQ, descF_PxQ, descG_PxQ);
   ene[3] = ebh + eba;
   t2 = seconds();
   *tmme2Bond += t2 - t1;
   t1 = t2;

   eth = eangl2(prm->Ntheth, prm->AngleHAt1, prm->AngleHAt2,
                prm->AngleHAt3, prm->AngleHNum,
                prm->Tk, prm->Teq, x, f, g,
                *context_PxQ, descF_PxQ, descG_PxQ);
   eta = eangl2(prm->Ntheta, prm->AngleAt1, prm->AngleAt2,
                prm->AngleAt3, prm->AngleNum,
                prm->Tk, prm->Teq, x, f, g,
                *context_PxQ, descF_PxQ, descG_PxQ);
   ene[4] = eth + eta;
   t2 = seconds();
   *tmme2Angl += t2 - t1;
   t1 = t2;

   eph = ephi2(prm->Nphih, prm->DihHAt1, prm->DihHAt2,
               prm->DihHAt3, prm->DihHAt4, prm->DihHNum,
               prm->Pk, prm->Pn, prm->Phase, x, f, g,
               *context_PxQ, descF_PxQ, descG_PxQ);
   epa = ephi2(prm->Mphia, prm->DihAt1, prm->DihAt2,
               prm->DihAt3, prm->DihAt4, prm->DihNum,
               prm->Pk, prm->Pn, prm->Phase, x, f, g,
               *context_PxQ, descF_PxQ, descG_PxQ);
   ene[5] = eph + epa;
   t2 = seconds();
   *tmme2Phi += t2 - t1;
   t1 = t2;

   ene[6] = 0.0;                /*  hbond term not in Amber-94 force field */

   nbond2(lpairs, prm->N14pairs, N14pearlist, 1, x, f, g, &enb14, &eel14,
          N14scnbpearlist, N14sceepearlist, *context_PxQ, descF_PxQ, descG_PxQ);
   t2 = seconds();
   *tmme2Nonb += t2 - t1;
   t1 = t2;

   if (e_debug) {
      EXPR("%9.3f", enb14);
      EXPR("%9.3f", eel14);
   }
   ene[7] = enb14;
   ene[8] = eel14;

   if (nconstrained) {
      ecn = econs2(x, f, g, *context_PxQ, descF_PxQ, descG_PxQ);
      t2 = seconds();
      *tmme2Cons += t2 - t1;
      t1 = t2;
   } else
      ecn = 0.0;
   ene[9] = ecn;

   /* Calculate the generalized Born energy and derivatives. */

   if (gb) {

#if defined(SCALAPACK) || defined(MPI)

      ene[10] =
          egb2(lpairs2, upairs2, pairlist2, lpairs2np, upairs2np,
               pairlist2np, x, f, g, grad, prm->Fs, prm->Rborn,
               prm->Charges, &kappa, &epsext, &enb, &eel, &esurf, &evdwnp,
               *context_Nx1, *context_1x1, *context_PxQ, descF_PxQ,
               descF_1x1, descG_PxQ, *gridim, 0);
#else

      ene[10] =
          egb2(lpairs, upairs, pairlist, lpairsnp, upairsnp, pairlistnp, x,
               f, g, grad, prm->Fs, prm->Rborn, prm->Charges, &kappa,
               &epsext, &enb, &eel, &esurf, &evdwnp, *context_Nx1,
               *context_1x1, *context_PxQ, descF_PxQ, descF_1x1, descG_PxQ,
               *gridim, 0);
#endif

      t2 = seconds();
      *tmme2Born += t2 - t1;
      t1 = t2;
      ene[1] = enb;
      ene[2] = eel;
      ene[11] = esurf;
      ene[12] = evdwnp;
   } else {

#if defined(SCALAPACK) || defined(MPI)
      nbond2(lpairs2, upairs2, pairlist2, 0, x, f, g, &enb, &eel, NULL, NULL,
             *context_PxQ, descF_PxQ, descG_PxQ);
#else
      nbond2(lpairs, upairs, pairlist, 0, x, f, g, &enb, &eel, NULL, NULL,
             *context_PxQ, descF_PxQ, descG_PxQ);
#endif

      t2 = seconds();
      *tmme2Nonb += t2 - t1;
      t1 = t2;
      ene[1] = enb;
      ene[2] = eel;
      ene[10] = 0.0;
      ene[11] = 0.0;
      ene[12] = 0.0;
   }

   /*
    * Zero out the frozen forces.  Use the local leading dimension of
    * the distributed Hessian matrix for the ScaLAPACK implementation.
    */

#pragma omp parallel private (k, j, kk, jj, threadnum, numthreads)
   {
      /*
       * Get the thread number and the number of threads for multi-threaded
       * execution.  Use default values of 0 and 1, respectively, for
       * single-threaded execution.
       */

#ifdef SCALAPACK
      threadnum = myrow;
      numthreads = nprow;
#elif defined(OPENMP)
      threadnum = omp_get_thread_num();
      numthreads = omp_get_num_threads();
#else
      threadnum = 0;
      numthreads = 1;
#endif

      /*
       * Use blocksize to assign the loop index i to OpenMP threads
       * or ScaLAPACK process rows.  This approach is identical to
       * the approach used in the egb2 function.
       */

#pragma omp for

      for (k = 0; k < prm->Natom; k++) {

#if defined(OPENMP) || defined(SCALAPACK)
         if (!myroc(k, blocksize, numthreads, threadnum))
            continue;
#endif

         kk = (size_t) k;

         if (frozen[k]) {

#ifndef SCALAPACK
            f[3 * k] = 0.0;
            f[3 * k + 1] = 0.0;
            f[3 * k + 2] = 0.0;
#else
            /*
             * Only processes in grid column 0 access the gradient vector.
             * The correct grid row was selected by the above continue statement,
             * which selects three rows at a time because the block size of
             * the gradient vector is three times blocksize.
             */

            if (mycol == 0) {
               ptr = ptr1d(f, descF_PxQ, 3 * k + 0);
               if (ptr != NULL)
                  *ptr = 0.0;
               ptr = ptr1d(f, descF_PxQ, 3 * k + 1);
               if (ptr != NULL)
                  *ptr = 0.0;
               ptr = ptr1d(f, descF_PxQ, 3 * k + 2);
               if (ptr != NULL)
                  *ptr = 0.0;
            }
#endif
         }
      }

      /* Zero the [k,j] elements of the Hessian matrix. */

      for (k = 0; k < prm->Natom; k++) {

         if (!myroc(k, blocksize, numthreads, threadnum))
            continue;
         kk = (size_t) k;

         if (frozen[k]) {
            for (j = 0; j < prm->Natom; j++) {
               jj = (size_t) j;

#ifndef SCALAPACK
               g[3 * jj + 0 + n * (3 * kk + 0)] = 0.0;
               g[3 * jj + 1 + n * (3 * kk + 1)] = 0.0;
               g[3 * jj + 2 + n * (3 * kk + 2)] = 0.0;
#else
               /*
                * The correct grid row was selected by the above continue statement.
                * Determine the correct processor grid column for Hessian matrix
                * access using the ptr2d function.
                */

               ptr = ptr2d(g, descG_PxQ, 3 * kk + 0, 3 * jj + 0);
               if (ptr != NULL)
                  *ptr = 0.0;
               ptr = ptr2d(g, descG_PxQ, 3 * kk + 1, 3 * jj + 1);
               if (ptr != NULL)
                  *ptr = 0.0;
               ptr = ptr2d(g, descG_PxQ, 3 * kk + 2, 3 * jj + 2);
               if (ptr != NULL)
                  *ptr = 0.0;
#endif
            }
         }
      }

#ifdef SCALAPACK

      /* This step is unnecessary for a square process grid. */

      threadnum = mycol;
      numthreads = npcol;

#endif

      /* Zero the [j,k] elements of the Hessian matrix. */

      for (k = 0; k < prm->Natom; k++) {

         if (!myroc(k, blocksize, numthreads, threadnum))
            continue;
         kk = (size_t) k;

         if (frozen[k]) {

            for (j = 0; j < prm->Natom; j++) {
               jj = (size_t) j;

#ifndef SCALAPACK
               g[3 * kk + 0 + n * (3 * jj + 0)] = 0.0;
               g[3 * kk + 1 + n * (3 * jj + 1)] = 0.0;
               g[3 * kk + 2 + n * (3 * jj + 2)] = 0.0;
#else
               /*
                * The correct grid column was selected by the above continue statement.
                * Determine the correct processor grid row for Hessian matrix
                * access using the ptr2d function.
                */

               ptr = ptr2d(g, descG_PxQ, 3 * jj + 0, 3 * kk + 0);
               if (ptr != NULL)
                  *ptr = 0.0;
               ptr = ptr2d(g, descG_PxQ, 3 * jj + 1, 3 * kk + 1);
               if (ptr != NULL)
                  *ptr = 0.0;
               ptr = ptr2d(g, descG_PxQ, 3 * jj + 2, 3 * kk + 2);
               if (ptr != NULL)
                  *ptr = 0.0;
#endif
            }
         }
      }
   }

   /*
    * Calculate the RMS error of the gradient vector.
    *
    * For non-ScaLAPACK execution, set some variables to values that
    * will select the proper sections of code below.
    */

#ifndef SCALAPACK

   *gridim = 1;
   *context_PxQ = 0;

#endif

   /*
    * Some ScaLAPACK functions appear to quit unexpectedly
    * for large matrices on a 1x1 process grid, e.g. the pdgemm_
    * function in newton.c, so bypass the ScaLAPACK pdgemr2d_
    * function and calculate the RMS error using f instead of grad.
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

   if (*gridim == 1) {
      if (*context_PxQ >= 0) {
         frms = 0.0;
         for (i = 0; i < n; i++)
            frms += f[i] * f[i];
         frms = sqrt(frms / ((REAL_T) n));
      }
   } else {

#ifdef SCALAPACK

      /*
       * For ScaLAPACK, begin by gathering the distributed gradient vector f
       * into the non-distributed gradient vector grad.  The function pdgemr2d_
       * will hang unless called from all processes.
       */

      pdgemr2d_(&n, &one,
                f, &one, &one, descF_PxQ,
                grad, &one, &one, descF_1x1, context_Nx1);

      /*
       * Calculate the RMS error for the gradient vector grad in context_1x1
       * where grad exists.
       */

      if (*context_1x1 >= 0) {
         frms = 0.0;
         for (i = 0; i < n; i++) {
            frms += grad[i] * grad[i];
         }
         frms = sqrt(frms / ((REAL_T) n));
      }
#endif

   }

   /* Calculate the total energy. */

   ene[0] = 0.0;
   for (k = 1; k <= 12; k++) {
      ene[0] += ene[k];
   }

   /*
    * If SCALAPACK is defined perform a reduction of the ene array,
    * but only the vdW, Coulombic, Born and nonpolar vdW energies have
    * been computed in parallel by egb2, not by nbond2 which does
    * not compute in parallel, and similarly not the surface area
    * term of the non-polar Born energy that was not computed in
    * parallel.
    */

#ifdef SCALAPACK
   if (gb) {
      MPI_Allreduce(ene, reductarr, 13, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      ene[1] = reductarr[1];
      ene[2] = reductarr[2];
      ene[10] = reductarr[10];
      ene[12] = reductarr[12];
   }
#endif

   /*
    * Print the energies and rms gradient but only for positive values
    * of the iteration counter, and if SCALAPACK is defined, only for
    * context_1x1 where the gradient vector grad is valid.
    */

#ifdef SCALAPACK
   if (*context_1x1 >= 0)
#elif defined(MPI)
   if (mytaskid == 0)
#endif
   {
      if (*iter > -1 && (*iter == 1 || *iter % ntpr == 0)) {
         printf("ff:%6d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2e\n",
                *iter, ene[0], ene[3] + ene[4] + ene[5],
                ene[1] + ene[7], ene[2] + ene[8],
                ene[9] + ene[11] + ene[12], ene[10], frms);
         fflush(stdout);
      }
   }

   /*
    * Load the atomic masses into the array m because the calling
    * function (newton) needs these masses to "level shift" to high
    * frequencies the eigenvalues associated with overall translation
    * and rotation of the molecule.
    *
    * Load the number of atoms as well.
    */

   for (i = 0; i < prm->Natom; i++) {
      m[i] = prm->Masses[i];
   }

   *natom = prm->Natom;

   /* A value of -1 for the iteration counter is reserved for printing. */

   if (*iter == -1) {
      printf("Inside mme2:\n");
      printf("     bond:  %15.9f\n", ene[3]);
      printf("    angle:  %15.9f\n", ene[4]);
      printf(" dihedral:  %15.9f\n", ene[5]);
      printf("    enb14:  %15.9f\n", ene[7]);
      printf("    eel14:  %15.9f\n", ene[8]);
      printf("      enb:  %15.9f\n", ene[1]);
      printf("      eel:  %15.9f\n", ene[2]);
      printf("      egb:  %15.9f\n", ene[10]);
      printf("    econs:  %15.9f\n", ene[9]);
      printf("    esurf:  %15.9f\n", ene[11]);
      printf("    Total:  %15.9f\n", ene[0]);
   }

   t2 = seconds();
   *tmme2Other += t2 - t1;

   /*
    * The following code is useful in debugging the gradient and Hessian.
    *
    * Compare the 1st derivatives from egb to the 1st derivatives from egb2.
    */

#undef GRADIENT_DEBUG
#ifdef GRADIENT_DEBUG

   REAL_T *f1 = vector(0, n);
   int iter_dum = 1;
   mme34(x, f1, &iter_dum);

   frms = 0.0;

   for (i = 0; i < n; i++) {
      frms += (f[i] - f1[i]) * (f[i] - f1[i]);
   }
   printf("rms delta df = %9.2e\n", sqrt(frms) / ((REAL_T) n));

   free_vector(f1, 0, n);

#endif


   /* Construct the 2nd derivatives by finite difference from the 1st derivatives. */

#undef HESSIAN_DEBUG
#ifdef HESSIAN_DEBUG

#define DELX (0.00001)
#define MULX (50000.0)

   ntpr = 2;
   int iter_dumm = 3;
   REAL_T eplus, eminus;
   REAL_T *e = vector(0, n);
   REAL_T *h = vector(0, n * n);
   REAL_T *d = vector(0, n * n);
   REAL_T *fplus = vector(0, n);
   REAL_T *fminus = vector(0, n);
   REAL_T maxdev = 0.0;
   REAL_T mindev = 1.0e10;
   REAL_T dev;
   for (i = 0; i < n; i++) {
      x[i] += DELX;
      eplus = mme34(x, fplus, &iter_dumm);
      x[i] -= 2.0 * DELX;
      eminus = mme34(x, fminus, &iter_dumm);
      x[i] += DELX;
      e[i] = MULX * (eplus - eminus);
      for (j = 0; j < n; j++) {
         h[n * i + j] = MULX * (fplus[j] - fminus[j]);
      }
   }

   /*
    * Calculate the rms difference between the first derivatives calculated
    * analytically and via finite differences.
    */

   maxdev = 0.0;
   mindev = 1.0e10;
   k = 0;
   frms = 0.0;
   for (i = 0; i < n; i++) {
      k++;
      dev = (e[i] - f[i]) * (e[i] - f[i]);
      frms += dev;
      if (dev > maxdev) {
         maxdev = dev;
      }
      if (dev < mindev) {
         mindev = dev;
      }
   }
   printf("rms delta df = %9.2e   maxdev = %9.2e   mindev = %9.2e\n",
          sqrt(frms) / ((REAL_T) k), maxdev, mindev);

   /*
    * Calculate the rms difference between the second derivatives calculated
    * analytically and via finite differences.
    */

   maxdev = 0.0;
   mindev = 1.0e10;
   k = 0;
   frms = 0.0;
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         k++;
         dev =
             (h[n * i + j] - g[n * i + j]) * (h[n * i + j] - g[n * i + j]);
         frms += dev;
         dev = sqrt(dev);
         d[n * i + j] = dev;
         if (dev > maxdev) {
            maxdev = dev;
         }
         if (dev < mindev) {
            mindev = dev;
         }
      }
   }
   printf("rms delta dg = %9.2e   maxdev = %9.2e   mindev = %9.2e\n",
          sqrt(frms) / ((REAL_T) k), maxdev, mindev);

   /*
    * Calculate the rms difference between the second derivatives
    * for the upper and lower triangles of the Hessian.
    */

   maxdev = 0.0;
   mindev = 1.0e10;
   k = 0;
   frms = 0.0;
   for (i = 0; i < n; i++) {
      for (j = i + 1; j < n; j++) {
         k++;
         dev =
             (h[n * i + j] - g[n * i + j]) * (h[n * i + j] - g[n * i + j]);
         frms += dev;
         dev = sqrt(dev);
         d[n * i + j] = dev;
         if (dev > maxdev) {
            maxdev = dev;
         }
         if (dev < mindev) {
            mindev = dev;
         }
      }
   }
   printf("rms delta dg = %9.2e   maxdev = %9.2e   mindev = %9.2e\n",
          sqrt(frms) / ((REAL_T) k), maxdev, mindev);

   maxdev = 0.0;
   mindev = 1.0e10;
   k = 0;
   frms = 0.0;
   for (j = 0; j < n; j++) {
      for (i = 0; i < j; i++) {
         k++;
         dev =
             (h[n * i + j] - g[n * i + j]) * (h[n * i + j] - g[n * i + j]);
         frms += dev;
         dev = sqrt(dev);
         d[n * i + j] = dev;
         if (dev > maxdev) {
            maxdev = dev;
         }
         if (dev < mindev) {
            mindev = dev;
         }
      }
   }
   printf("rms delta dg = %9.2e   maxdev = %9.2e   mindev = %9.2e\n",
          sqrt(frms) / ((REAL_T) k), maxdev, mindev);

   /* Use the finite-difference Hessian not the analytic Hessian. */

#undef USE_FINITE_DIFFERENCE
#ifdef USE_FINITE_DIFFERENCE
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         g[n * i + j] = h[n * i + j];
      }
   }
#endif

#undef PRINT_DEVIATION
#ifdef PRINT_DEVIATION
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         printf("deviation[%d][%d] = %9.2e\n", i, j, d[i][j]);
      }
      printf("\n");
   }
   printf("\n");
#endif

   ntpr = 1;
   free_vector(fminus, 0, 3 * prm->Natom);
   free_vector(fplus, 0, 3 * prm->Natom);
   free_vector(h, 0, n * n);
   free_vector(d, 0, n * n);
   free_vector(e, 0, n);
#endif

   /* Check that the Hessian is symmetric. */

#undef TEST_SYMMETRY
#ifdef TEST_SYMMETRY
   k = 0;
   frms = 0.0;
   for (i = 0; i < n; i++) {
      for (j = i + 1; j < n; j++) {
         k++;
         frms +=
             (g[n * i + j] - g[n * j + i]) * (g[n * i + j] - g[n * j + i]);
      }
   }
   printf("rms symmetry = %9.2e\n", sqrt(frms / ((REAL_T) k)));
#endif

   /* Return the total energy as the value of mme2. */

   t2 = seconds();
   *tmme2Other += t2 - t1;
   *tmme2 += t2 - tmme21;
   return (ene[0]);
}

