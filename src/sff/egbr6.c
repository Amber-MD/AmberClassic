
/***********************************************************************
                            EGBr6()
************************************************************************/

/*
 * Calculate the generalized Born energy and first derivatives.
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
 * fs - input: overlap parameters
 * rborn - input: atomic radii
 * q - input: atomic charges
 * kappa - input: inverse of the Debye-Huckel length
 * diel_ext - input: solvent dielectric constant
 * enb - updated: Lennard-Jones energy
 * eelt - updated: gas-phase electrostatic energy
 * esurf - updated: nonpolar surface area solvation free energy
 * enp - updated: nonpolar van der Waals solvation free energy
 * freevectors - if !=0 free the static vectors and return
 */

#define GOFFSET 0.0
static
REAL_T egbr6(INT_T * lpears, INT_T * upears, INT_T ** pearlist,
             INT_T * lpearsnp, INT_T * upearsnp, INT_T ** pearlistnp,
             REAL_T * x, REAL_T * f, REAL_T * fs, REAL_T * rborn,
             REAL_T * q, REAL_T * kappa, REAL_T * diel_ext, REAL_T * enb,
             REAL_T * eelt, REAL_T * esurf, REAL_T * enp,
             INT_T freevectors)
{
#if defined(MPI) || defined(SCALAPACK)
   int ierror;
   static REAL_T *reductarr = NULL;
#endif

   static REAL_T *reff = NULL, *sumdeijda = NULL, *psi = NULL;

   int i, i34, j, j34, k, threadnum, numthreads, maxthreads, foff, soff;
   int npairs, ic, iaci;
   int *iexw;
   size_t n;
   REAL_T epol, dielfac, qi, qiqj, fgbi, fgbk, rb2, expmkf;
   REAL_T elec, evdw, sumda, daix, daiy, daiz;
   REAL_T xi, yi, zi, xij, yij, zij;
   REAL_T dedx, dedy, dedz, de;
   REAL_T dij1i, dij3i, temp1;
   REAL_T qi2h, qid2h, datmp;
   REAL_T ri1i, dij2i;
   REAL_T theta;

   REAL_T dij, sumi;
   REAL_T eel, f6, f12, rinv, r2inv, r6inv;
   REAL_T r2, ri, rj, sj, sj2;
   REAL_T uij, efac, temp4, temp5, temp6, reffi2;

   /*
    * Determine the size of the sumdeijda array.  If OPENMP is defined,
    * a copy of this array must be allocated for each thread; otherwise
    * only one copy is allocated.
    */

#ifdef OPENMP
   maxthreads = omp_get_max_threads();
#else
   maxthreads = 1;
#endif

   n = (size_t) prm->Natom;

   *enp = 0.0;

   /*
    * If freevectors != 0, deallocate the static arrays that have been
    * previously allocated and return.
    */

   if (freevectors != 0) {
      if (reff != NULL)
         free_vector(reff, 0, n);
      reff = NULL;
      if (sumdeijda != NULL)
         free_vector(sumdeijda, 0, maxthreads * n);
      sumdeijda = NULL;
      if (psi != NULL)
         free_vector(psi, 0, n);
      psi = NULL;
#if defined(MPI) || defined(SCALAPACK)
      if (reductarr != NULL)
         free_vector(reductarr, 0, n);
      reductarr = NULL;
#endif
      return (0.0);
   }


   /* Allocate some static arrays if they have not been allocated already. */

   if (reff == NULL)
      reff = vector(0, n);
   if (sumdeijda == NULL)
      sumdeijda = vector(0, maxthreads * n);
   if ((psi == NULL) && (gb == 2 || gb == 5))
      psi = vector(0, n);
#if defined(MPI) || defined(SCALAPACK)
   if (reductarr == NULL)
      reductarr = vector(0, n);
#endif

   if (gb_debug)
      fprintf(nabout, "Effective Born radii:\n");


   /*
    * For MPI or ScaLAPACK, initialize all elements of the reff array.
    * Although each task will calculate only a subset of the elements,
    * a reduction is used to combine the results from all tasks.
    * If a gather were used instead of a reduction, no initialization
    * would be necessary.
    */

#if defined(MPI) || defined(SCALAPACK)
   for (i = 0; i < prm->Natom; i++) {
      reff[i] = 0.0;
   }
#endif

#pragma omp parallel \
  private (i, xi, yi, zi, ri, ri1i, sumi, j, k, xij, yij, zij, \
           r2, dij1i, dij, sj, sj2, uij, dij2i, theta, temp1, \
           threadnum, numthreads)
   {

      /*
       * Get the thread number and the number of threads for multi-threaded
       * execution under OpenMP.  For all other cases, including ScaLAPACK,
       * MPI and single-threaded execution, use the values that have been
       * stored in mytaskid and numtasks, respectively.
       */

#if defined(OPENMP)
      threadnum = omp_get_thread_num();
      numthreads = omp_get_num_threads();
#else
      threadnum = mytaskid;
      numthreads = numtasks;
#endif

      /*
       * Loop over all atoms i.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       *
       * The reff array is written in the following loops.  It is necessary to
       * synchronize the OpenMP threads or MPI tasks that execute these loops
       * following loop execution so that a race condition does not exist for
       * reading the reff array before it is written.  Even if all subsequent
       * loops use loop index to thread or task mapping that is identical to
       * that of the following loop, elements of the reff array are indexed by
       * other loop indices, so synchronization is necessary.
       *
       * OpenMP synchronization is accomplished by the implied barrier
       * at the end of this parallel region.  MPI synchronization is
       * accomplished by MPI_Allreduce.
       */

      for (i = 0; i < prm->Natom; i++) {

#if defined(OPENMP) || defined(MPI) || defined(SCALAPACK)
         if (!myroc(i, blocksize, numthreads, threadnum))
            continue;
#endif

         xi = x[3 * i];
         yi = x[3 * i + 1];
         zi = x[3 * i + 2];

         ri = rborn[i] - GOFFSET;
         ri1i = 1. / ri;
         sumi = 0.0;

         /* Select atom j from the pair list.  Non-graceful error handling. */

         for (k = 0; k < lpears[i] + upears[i]; k++) {

            if (pearlist[i] == NULL) {
               fprintf(nabout,
                       "NULL pair list entry in egb loop 1, taskid = %d\n",
                       mytaskid);
               fflush(nabout);
            }
            j = pearlist[i][k];

            xij = xi - x[3 * j];
            yij = yi - x[3 * j + 1];
            zij = zi - x[3 * j + 2];
            r2 = xij * xij + yij * yij + zij * zij;

            dij1i = 1.0 / sqrt(r2);
            dij = r2 * dij1i;
            sj = fs[j] * (rborn[j] - GOFFSET);
            sj2 = sj * sj;

            /* Contributions to the inverse cube effective radii; see
               Eqs. 10-12 of Tjong & Zhou, JPCB 111:3055, 2007        */

            if (dij > ri + sj) {        /* Case I */

               temp1 = sj / (r2 - sj2);
               sumi -= temp1 * temp1 * temp1;

            } else if (dij > fabs(ri - sj)) {   /* Case II */
               temp1 = 1. / (dij + sj);
               sumi -=
                   0.0625 * dij1i * (8. * dij *
                                     (1. / (ri * ri * ri) -
                                      temp1 * temp1 * temp1)
                                     - 3. * (r2 - sj2) * (pow(ri, -4.) -
                                                          pow(temp1, 4.))
                                     - 6. * (1. / (ri * ri) -
                                             temp1 * temp1));

            } else if (ri < sj) {       /* Case III */
               temp1 = sj / (r2 - sj2);
               sumi -= dij1i / r2 + temp1 * temp1 * temp1;

            }

         }

         temp1 = ri1i * ri1i * ri1i + sumi;
         if (temp1 <= 0.0) {
            reff[i] = 30.0;
         } else {
            reff[i] = pow(temp1, (-1. / 3.));
         }

         if (gb_debug)
            fprintf(nabout, "%d\t%15.7f\t%15.7f\n", i + 1, rborn[i],
                    reff[i]);
      }
   }

   /* The MPI synchronization is accomplished via reduction of the reff array. */

#if defined(MPI) || defined(SCALAPACK)

   ierror = MPI_Allreduce(reff, reductarr, prm->Natom,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (ierror != MPI_SUCCESS) {
      fprintf(nabout,
              "Error in egb reff reduction, error = %d  mytaskid = %d\n",
              ierror, mytaskid);
   }
   for (i = 0; i < prm->Natom; i++) {
      reff[i] = reductarr[i];
   }

#endif

   *esurf = 0.0;                /* no LCPO in this version  */

   /* Compute the GB, Coulomb and Lennard-Jones energies and derivatives. */

   epol = elec = evdw = 0.0;

#pragma omp parallel reduction (+: epol, elec, evdw) \
  private (i, i34, ri, qi, expmkf, dielfac, qi2h, qid2h, iaci, \
           xi, yi, zi, k, j, j34, xij, yij, zij, r2, qiqj, \
           rj, rb2, efac, fgbi, fgbk, temp4, temp6, eel, de, temp5, \
           rinv, r2inv, ic, r6inv, f6, f12, dedx, dedy, dedz, \
           iexw, threadnum, numthreads, foff, soff, \
           sumda, ri1i, dij1i, datmp, daix, daiy, daiz, \
           dij2i, dij, sj, sj2, temp1, dij3i, npairs)
   {

      /*
       * Get the thread number and the number of threads for multi-threaded
       * execution under OpenMP.  For all other cases, including ScaLAPACK,
       * MPI and single-threaded execution, use the values that have been
       * stored in mytaskid and numtasks, respectively.
       */

#if defined(OPENMP)
      threadnum = omp_get_thread_num();
      numthreads = omp_get_num_threads();
#else
      threadnum = mytaskid;
      numthreads = numtasks;
#endif

      /*
       * Compute offsets into the gradient and sumdeijda arrays for this
       * thread, but only if OPENMP is defined.
       */

#ifdef OPENMP
      soff = prm->Natom * threadnum;
      foff = 3 * soff;
#else
      soff = 0;
      foff = 0;
#endif

      /*
       * Initialize the sumdeijda array inside of the parallel region.
       *
       * For OpenMP, the "first touch" memory allocation strategy will
       * locate the copy of the sumdeijda array that is initialized by
       * a particular CPU in memory that is local to that CPU.
       *
       * It is not necessary to synchronize OpenMP threads following
       * this loop because each thread initializes the particular copy
       * of the sumdeijda array that it subsequently updates.  A similar
       * argument applies for MPI tasks.
       */

      for (i = 0; i < prm->Natom; i++) {
         sumdeijda[soff + i] = 0.0;
      }

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

      /*
       * Loop over all atoms i.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       *
       * Synchronization of OpenMP threads will occur following this loop
       * because the parallel region ends after this loop.  Following
       * synchronization, a reduction of the sumdeijda array will be
       * performed.
       *
       * Synchronization of MPI tasks will occur via the MPI_Reduce function.
       */

      for (i = 0; i < prm->Natom; i++) {

#if defined(OPENMP) || defined(MPI) || defined(SCALAPACK)
         if (!myroc(i, blocksize, numthreads, threadnum))
            continue;
#endif

         ri = reff[i];
         qi = q[i];

         /*
          * If atom i is not frozen, compute the "diagonal" energy that
          * is a function of only the effective radius Ri but not of the
          * interatomic distance Dij.  Compute also the contribution of
          * the diagonal energy term to the sum by which the derivative
          * of Ri will be multiplied.
          */

         if (!frozen[i]) {
            expmkf = exp(-KSCALE * (*kappa) * ri) / (*diel_ext);
            dielfac = 1.0 - expmkf;
            qi2h = 0.5 * qi * qi;
            qid2h = qi2h * dielfac;
            epol += -qid2h / ri;

            sumdeijda[soff + i] +=
                qid2h - KSCALE * (*kappa) * qi2h * expmkf * ri;
         }

         /*
          * Skip the pair calculations if there are no atoms j on the
          * pair list of atom i.
          */

         npairs = upears[i];
         if (!npairs)
            continue;

         i34 = 3 * i;

         xi = x[i34];
         yi = x[i34 + 1];
         zi = x[i34 + 2];

         iaci = prm->Ntypes * (prm->Iac[i] - 1);

         /*
          * Expand the excluded atom list into the iexw array by storing i
          * at array address j.
          */

         for (j = 0; j < prm->Iblo[i]; j++) {
            iexw[IexclAt[i][j] - 1] = i;
         }

         /* Initialize the derivative accumulators. */

         daix = daiy = daiz = 0.0;

         /* Select atoms j from the pair list.  Non-graceful error handling. */

         for (k = lpears[i]; k < lpears[i] + npairs; k++) {

            if (pearlist[i] == NULL) {
               fprintf(nabout,
                       "NULL pair list entry in egb loop 3, taskid = %d\n",
                       mytaskid);
               fflush(nabout);
            }
            j = pearlist[i][k];

            j34 = 3 * j;

            /* Continue computing the non-diagonal energy term. */

            xij = xi - x[j34];
            yij = yi - x[j34 + 1];
            zij = zi - x[j34 + 2];
            r2 = xij * xij + yij * yij + zij * zij;

            /*
             * Because index j is retrieved from the pairlist array it is
             * not constrained to a particular range of values; therefore,
             * the threads that have loaded the reff array must be
             * synchronized prior to the use of reff below.
             */

            qiqj = qi * q[j];
            rj = reff[j];
            rb2 = ri * rj;
#ifdef GRYCUK_FGB
            fgbi = 1.0 / sqrt(r2 + rb2);
#else
            efac = exp(-r2 / (4.0 * rb2));
            fgbi = 1.0 / sqrt(r2 + rb2 * efac);
#endif
            fgbk = -(*kappa) * KSCALE / fgbi;

            expmkf = exp(fgbk) / (*diel_ext);
            dielfac = 1.0 - expmkf;

            epol += -qiqj * dielfac * fgbi;

            temp4 = fgbi * fgbi * fgbi;
            temp6 = qiqj * temp4 * (dielfac + fgbk * expmkf);

            /* de holds (1/r) dE/dr, where r is the bond distance  */
#ifdef GRYCUK_FGB
            de = temp6;
#else
            de = temp6 * (1.0 - 0.25 * efac);
#endif

            /* temp5 holds -(1/ri) dE/dai,  where ai = 1/ri; this is
               also equal to -(1/rj) dE/daj                         */
#ifdef GRYCUK_FGB
            temp5 = 0.5 * temp6 * rb2;
#else
            temp5 = 0.5 * efac * temp6 * (rb2 + 0.25 * r2);
#endif

            /*
             * sumdeijda[] will hold minus the derivative of the GB
             * energy with respect to the inverse effective radii.
             */

            sumdeijda[soff + i] += ri * temp5;
            sumdeijda[soff + j] += rj * temp5;

            /*
             * Compute the Van der Waals and Coulombic energies for only
             * those pairs that are not on the excluded atom list.  Any
             * pair on the excluded atom list will have atom i stored at
             * address j of the iexw array.  It is not necessary to reset
             * the elements of the iexw array to -1 between successive
             * iterations in i because an i,j pair is uniquely identified
             * by atom i stored at array address j.  Thus for example, the
             * i+1,j pair would be stored at the same address as the i,j
             * pair but after the i,j pair were used.
             */

            if (iexw[j] != i) {

               rinv = 1. / sqrt(r2);
               r2inv = rinv * rinv;

               /*  gas-phase Coulomb energy:  */

               eel = qiqj * rinv;
               elec += eel;
               de -= eel * r2inv;

               /* Lennard-Jones energy:   */

               ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
               if (ic >= 0) {
                  r6inv = r2inv * r2inv * r2inv;
                  f6 = prm->Cn2[ic] * r6inv;
                  f12 = prm->Cn1[ic] * r6inv * r6inv;
                  evdw += f12 - f6;
                  de -= (12. * f12 - 6. * f6) * r2inv;
               }
            }

            /*
             * Sum to the gradient vector the derivatives of Dij that are
             * computed relative to the Cartesian coords of atoms i and j.
             */

            dedx = de * xij;
            dedy = de * yij;
            dedz = de * zij;

            daix += dedx;
            daiy += dedy;
            daiz += dedz;

            f[foff + j34] -= dedx;
            f[foff + j34 + 1] -= dedy;
            f[foff + j34 + 2] -= dedz;

         }

         /* Update the i elements of the gradient. */

         f[foff + i34] += daix;
         f[foff + i34 + 1] += daiy;
         f[foff + i34 + 2] += daiz;

      }

      /* Free the iexw array within this potentially parallel region of code. */

      free_ivector(iexw, -1, prm->Natom);
   }
#if defined(MPI) || defined(SCALAPACK)

   /* Perform a reduction of sumdeijda */

   ierror = MPI_Allreduce(sumdeijda, reductarr, prm->Natom,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (ierror != MPI_SUCCESS) {
      fprintf(nabout,
              "Error in egb sumdeijda reduction, error = %d  mytaskid = %d\n",
              ierror, mytaskid);
   }
   for (i = 0; i < prm->Natom; i++) {
      sumdeijda[i] = reductarr[i];
   }
#endif
#ifdef OPENMP

   /*
    * Perform a reduction of sumdeijda.
    * In the following, the (j,i) loop nesting is more efficient
    * than (i,j) loop nesting.
    */

#pragma omp parallel for private (i) schedule(static, 8)
   for (j = 0; j < prm->Natom; j++) {
      for (i = 1; i < maxthreads; i++) {
         sumdeijda[j] += sumdeijda[prm->Natom * i + j];
      }
   }
#endif

#pragma omp parallel \
  private (i, i34, ri, qi, expmkf, dielfac, qi2h, qid2h, iaci, \
           xi, yi, zi, k, j, j34, xij, yij, zij, r2, qiqj, \
           rj, rb2, efac, fgbi, fgbk, temp4, temp6, eel, de, temp5, \
           rinv, r2inv, ic, r6inv, f6, f12, dedx, dedy, dedz, \
           iexw, threadnum, numthreads, foff, \
           sumda, ri1i, dij1i, datmp, daix, daiy, daiz, \
           dij2i, dij, sj, sj2, temp1, dij3i, reffi2, npairs)
   {
      /*
       * Get the thread number and the number of threads for multi-threaded
       * execution under OpenMP.  For all other cases, including ScaLAPACK,
       * MPI and single-threaded execution, use the values that have been
       * stored in mytaskid and numtasks, respectively.
       */

#if defined(OPENMP)
      threadnum = omp_get_thread_num();
      numthreads = omp_get_num_threads();
#else
      threadnum = mytaskid;
      numthreads = numtasks;
#endif

      /*
       * Compute an offset into the gradient array for this thread,
       * but only if OPENMP is defined.  There is no need to
       * compute an offset into the sumdeijda array because all
       * copies of this array have been reduced into copy zero.
       */

#ifdef OPENMP
      foff = prm->Natom * 3 * threadnum;
#else
      foff = 0;
#endif

      /*
       * Compute the derivatives of the effective radius Ri of atom i
       * with respect to the cartesian coordinates of each atom j.  Sum
       * all of these derivatives into the gradient vector.
       *
       * Loop over all atoms i.
       *
       * Synchronization of OpenMP threads will occur following this loop
       * because the parallel region ends after this loop.  A reduction
       * of the gradient array will occur in the mme34 function, either
       * for OpenMP or MPI.  This reduction will synchronize the MPI
       * tasks, so an explicit barrier is not necessary at the end of
       * this loop.
       *
       * Explicitly assign threads to loop indices for the following loop,
       * in a manner equivalent to (static, N) scheduling with OpenMP, and
       * identical to the manner in which threads are assigned in nblist.
       */

      for (i = 0; i < prm->Natom; i++) {

#if defined(OPENMP) || defined(MPI) || defined(SCALAPACK)
         if (!myroc(i, blocksize, numthreads, threadnum))
            continue;
#endif

         /*
          * Don't calculate derivatives of the effective radius of atom i
          * if atom i is frozen or if there are no pair atoms j associated
          * with atom i.
          */

         npairs = lpears[i] + upears[i];
         if (frozen[i] || !npairs)
            continue;

         i34 = 3 * i;

         xi = x[i34];
         yi = x[i34 + 1];
         zi = x[i34 + 2];

         ri = rborn[i] - GOFFSET;
         ri1i = 1. / ri;
         reffi2 = reff[i] * reff[i];

         sumda = sumdeijda[i];

         /* Initialize the derivative accumulators. */

         daix = daiy = daiz = 0.0;

         /* Select atom j from the pair list.  Non-graceful error handling. */

         for (k = 0; k < npairs; k++) {

            if (pearlist[i] == NULL) {
               fprintf(nabout,
                       "NULL pair list entry in egb loop 4, taskid = %d\n",
                       mytaskid);
               fflush(nabout);
            }
            j = pearlist[i][k];

            j34 = 3 * j;

            xij = xi - x[j34];
            yij = yi - x[j34 + 1];
            zij = zi - x[j34 + 2];
            r2 = xij * xij + yij * yij + zij * zij;

            dij1i = 1.0 / sqrt(r2);
            dij2i = dij1i * dij1i;
            dij = r2 * dij1i;
            sj = fs[j] * (rborn[j] - GOFFSET);
            sj2 = sj * sj;

            /*  datmp will hold -(1/d) D[ ai, d ], where ai is the
               inverse of the effective radius, and d is the distance
               between atoms i and j.                                 */

            if (dij > ri + sj) {        /* Case I */

               datmp = reffi2 * 2. * sj * sj2 * pow((r2 - sj2), -4.);

            } else if (dij > fabs(ri - sj)) {   /* Case II  */

               temp1 = 1. / (dij + sj);
               temp1 = temp1 * temp1 * temp1 * temp1;
               dij3i = dij1i * dij1i * dij1i;
               datmp = 0.0625 * reffi2 * dij3i *
                   ((r2 + sj2 - 2. * ri * ri) *
                    ri1i * ri1i * ri1i * ri1i +
                    (r2 + 4. * dij * sj + sj2) * temp1);

            } else if (ri < sj) {       /* Case III */

               datmp = reffi2 * 2. * sj * sj2 * pow((r2 - sj2), -4.);

            } else {
               datmp = 0.;
            }

            /* Sum the derivatives into daix, daiy and daiz. */

            daix += xij * datmp;
            daiy += yij * datmp;
            daiz += zij * datmp;

            datmp *= sumda;
            f[foff + j34] += xij * datmp;
            f[foff + j34 + 1] += yij * datmp;
            f[foff + j34 + 2] += zij * datmp;

         }

         f[foff + i34] -= sumda * daix;
         f[foff + i34 + 1] -= sumda * daiy;
         f[foff + i34 + 2] -= sumda * daiz;

      }
   }

   /* Return elec and evdw through the reference parameters eelt and enb. */

   *eelt = elec;
   *enb = evdw;

   return (epol);
}
