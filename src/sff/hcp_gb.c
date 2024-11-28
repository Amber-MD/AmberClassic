/**********************************************************************
* implementation of HCPO-GB and HCP-GB                                *
* Izadi, Aandakrishnan and Onufriev, JCTC 2016 (in prep)              *
* Aandakrishnan and Onufriev, JCTC, 2011                              *
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#define KSCALE (0.73)

/*FGB taylor coefficients follow */
/* from A to H :                 */
/* 1/3 , 2/5 , 3/7 , 4/9 , 5/11  */
/* 4/3 , 12/5 , 24/7 , 40/9 , 60/11 */

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

#define PACK_FACT 1.0    /* 1.1604 = 1/cube_root(0.64) */
#define LAMBDA 1.0       /* 1.4 makes it worse */

#define HCP_RAD 1        /* 0=vol based, 1=ao, 2=as */

/***********************************************************************
 *                          egb_hcp_rbondi()
 *
 * Calculate radius for higher level hcp components
 *
 * Calling Parameters:
 * prm    - parameter/topology structure    (input)
 * r_hcpx - (x=1/2) radius of residue/chain (update)
 *
 ************************************************************************/

int egb_hcp_rbondi(PARMSTRUCT_T *prm, REAL_T * r_hcp1, REAL_T * r_hcp2, REAL_T * r_hcp3) 
{

        int c, s, r, a, q;                 /* index for strand, residue, atom, charge */
        int a_from, a_to;                  /* for atoms within a residue */
        int r_from, r_to;		   /* for residues within a strand */
        int s_from, s_to;		   /* for strands within a complex */
        float temp_r;

        for (c = 0; c < prm->Ncomplex; c++) /* for each strand */
        {
                /* initialize complex radius */
                for (q = 0; q < hcp; q++)
                {
                        r_hcp3[c * hcp + q] = 0.0;
                }
                s_from = prm->Ipcomplex[c] - 1;
                if (c < prm->Ncomplex - 1) { s_to = prm->Ipcomplex[c + 1] - 1; }
                else { s_to = prm->Nstrand; }
                for (s = s_from; s < s_to; s++) /* for each strand */
                {
                        /* initialize strand radius */
                        for (q = 0; q < hcp; q++)
                        {
                                r_hcp2[s * hcp + q] = 0.0;
                        }
                        r_from = prm->Ipstrand[s] - 1;
                        if (s < prm->Nstrand - 1) { r_to = prm->Ipstrand[s + 1] - 1; }
                        else { r_to = prm->Nres; }
                        for (r = r_from; r < r_to; r++)         /* for each residue */
                        {
                                /* initialize residue radius */
                                for (q = 0; q < hcp; q++)
                                {
                                        r_hcp1[r * hcp + q] = 0.0;
                                }

                                a_from = prm->Ipres[r];
                                if (r + 1 < prm->Nres) { a_to = prm->Ipres[r + 1]; }
                                else { a_to = prm->Natom; }
                                for (a = a_from; a < a_to; a++)    /* for each atom */
                                {
                                        q = q2idx(prm->Charges[a]);
                                        temp_r = (prm->Rborn[a] - gboffset) * prm->Fs[a]; 
                                        r_hcp1[r * hcp + q] += (temp_r * temp_r * temp_r);
                                }  /* end for each atoms */

                                for (q = 0; q < hcp; q++)
                                {
                                        /* roll up to next level */
                                        r_hcp2[s * hcp + q] += r_hcp1[r * hcp + q];

                                        /* compute residue radius */
                                        r_hcp1[r * hcp + q] = cbrt(r_hcp1[r * hcp + q]) * PACK_FACT;
                                }

                        }  /* end for each residue */

                        for (q = 0; q < hcp; q++)
                        {
                                /* roll up to next level */
                                r_hcp3[c * hcp + q] += r_hcp2[s * hcp + q];
                                /* compute strand radius */
                                r_hcp2[s * hcp + q] = cbrt(r_hcp2[s * hcp + q]) * PACK_FACT;
                        }

                }  /* end for each strand */

                for (q = 0; q < hcp; q++)
                {
                        /* compute complex radius */
                        r_hcp3[c * hcp + q] = cbrt(r_hcp3[c * hcp + q]) * PACK_FACT;
                }
        }
        return(0);

}  /* end egb_hcp_rbondi() */



/***********************************************************************
 *                        egb_hcp_approx_reff()
 *
 * Calculate effective Born radii for hcp components
 *
 * Calling parameters are as follows:
 *
 * j             - HCP component j                   (input)
 * r2            - dist^2 between i and j            (input)
 * sj            - adjusted component radius for j   (input)
 * sumi          - contribution to Born radius       (output)
 *
 ************************************************************************/

static REAL_T ebg_hcp_approx_reff(int j, REAL_T r2, REAL_T sj)
{
        REAL_T dij1i, dij2i, sj2, tmpsd, dumbo, sumi;

       // why sqrt?
       // why j, inputed to the routine, is not used anywhere?

        dij2i = 1.0 /r2;
       // dij1i = 1.0 /sqrt(r2);
       // dij2i = dij1i * dij1i;

        sj2 = sj * sj;
    

        tmpsd = sj2 * dij2i;
        dumbo = TA + tmpsd * (TB + tmpsd * (TC +   tmpsd * (TD + tmpsd * TDD)));

        sumi = sj * tmpsd * dij2i * dumbo;

        return(sumi);
}


/***********************************************************************
 *                        egb_hcp_atom_reff()
 *
 * Calculate effective Born radii for atoms
 *
 * Calling parameters are as follows:
 *
 * x             - the atomic (x,y,z) coordinates         (input)
 * x_hcpx        - (x=1/2) geom center residue/chain      (input)
 * rborn         - atom bondi radii                       (input)
 * fs            - atom scaling (overlap) factor          (input)
 * r_hcpx        - (x=1/2) residue/chain radius           (input)
 * reff_hcp0     - atom effective born radii              (update) 
 * psi           - interim computed values for reuse      (update)
 *
 ************************************************************************/

int egb_hcp_atom_reff(REAL_T * x, REAL_T * x_hcp1, REAL_T * x_hcp2, REAL_T * x_hcp3,
                REAL_T * q_hcp1, REAL_T * q_hcp2, REAL_T * q_hcp3,
                REAL_T * rborn, REAL_T * fs, REAL_T * r_hcp1, REAL_T * r_hcp2, REAL_T * r_hcp3,
                REAL_T * reff_hcp0, REAL_T * psi)
{

        int i, q;
        int threadnum, numthreads, maxthreads, numcopies;
        REAL_T xi, yi, zi;
        REAL_T dij1i;
        REAL_T uij, theta, ri1i, dij2i;

        REAL_T dij, sumi, t1, t2;
        REAL_T r2, ri, sj, sj2;
        REAL_T dumbo, tmpsd;

        int c, s, r, a;      /* complex, strand, residue, atom index */
        int a_from, a_to;    /* for atoms within a residue */
        int r_from, r_to;    /* for residues within a strand */
        int s_from, s_to;    /* for strands within a complex */
        REAL_T dist2, dist2_hcp1, dist2_hcp2, dist2_hcp3;

        /* square of threshold distances */
        dist2_hcp1 = hcp_h1 * hcp_h1;
        dist2_hcp2 = hcp_h2 * hcp_h2;
        dist2_hcp3 = hcp_h3 * hcp_h3;


#if defined(MPI) || defined(SCALAPACK)
        int ierror;
        REAL_T *reductarr = NULL;
        MPI_Status status;
#endif

        /*
         * Determine the size of the iexw array.  If OPENMP is
         * defined, a copy of this array must be allocated for
         * each thread; otherwise, only one copy is allocated.
         */

        maxthreads = 1;
        numcopies = 1;

        /* 
         * Get the "effective" Born radii via the approximate pairwise method.
         * Use Eqs 9-11 of Hawkins, Cramer, Truhlar, J. Phys. Chem. 100:19824
         * (1996).
         *
         * For MPI or ScaLAPACK, initialize all elements of the reff array.
         * Although each task will calculate only a subset of the elements,
         * a reduction is used to combine the results from all tasks.
         * If a gather were used instead of a reduction, no initialization
         * would be necessary.
         */


        /*
         * Get the thread number and the number of threads for multi-threaded
         * execution under OpenMP.  For all other cases, including ScaLAPACK,
         * MPI and single-threaded execution, use the values that have been
         * stored in mytaskid and numtasks, respectively.
         */

        threadnum = mytaskid;
        numthreads = numtasks;

        /*
         * Loop over all atoms i.
         *
         * For MPI or ScaLAPACK, explicitly assign tasks to loop indices
         * for the following loop in a manner equivalent to (static, N)
         * scheduling for OpenMP.  For OpenMP use (dynamic, N) scheduling.
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
         * at the end of this 'pragma omp for'.  MPI synchronization is
         * accomplished by MPI_Allreduce.
         */

       // if (gb_debug) 
        //        printf("Effective Born Radii:\n");

        for (i = 0; i < prm->Natom; i++) 
        {

                reff_hcp0[i] = 0.0;

#if defined(MPI)
                if (!myroc(i, blocksize, numthreads, threadnum))
                        continue;
#endif

                xi = x[dim * i];
                yi = x[dim * i + 1];
                zi = x[dim * i + 2];

                ri = rborn[i] - gboffset;
                ri1i = 1. / ri;
                sumi = 0.0;

                for (c = 0; c < prm->Ncomplex; c++)   /* for each complex in the molecule */
                {
                        dist2 = calc_dist2(xi, yi, zi,
                                        x_hcp3[3*c+0], x_hcp3[3*c+1], x_hcp3[3*c+2]);
                        if (dist2 > dist2_hcp3)
                        {
                                for (q = 0; q < hcp; q++)
                                {
                                        dist2 = calc_dist2(xi, yi, zi, q_hcp3[4*(c*hcp+q)], q_hcp3[4*(c*hcp+q)+1], q_hcp3[4*(c*hcp+q)+2]);
                                        sumi -= ebg_hcp_approx_reff(c, dist2, r_hcp3[c * hcp + q]);
                                }
                        }
                        else
                        {
                                s_from = prm->Ipcomplex[c] - 1;
                                if (c + 1 < prm->Ncomplex) { s_to = prm->Ipcomplex[c + 1] - 1; }
                                else { s_to = prm->Nstrand; }

                                for (s = s_from; s < s_to; s++)   /* for each strand in the complex */
                                {
                                        dist2 = calc_dist2(xi, yi, zi,
                                                        x_hcp2[3*s+0], x_hcp2[3*s+1], x_hcp2[3*s+2]);
                                        if (dist2 > dist2_hcp2)
                                        {
                                                for (q = 0; q < hcp; q++)
                                                {
                                                        dist2 = calc_dist2(xi, yi, zi, q_hcp2[4*(s*hcp+q)], q_hcp2[4*(s*hcp+q)+1], q_hcp2[4*(s*hcp+q)+2]);
                                                        sumi -= ebg_hcp_approx_reff(s, dist2, r_hcp2[s * hcp + q]);
                                                }
                                        }
                                        else
                                        {
                                                r_from = prm->Ipstrand[s] - 1;
                                                if (s + 1 < prm->Nstrand) { r_to = prm->Ipstrand[s + 1] - 1; }
                                                else { r_to = prm->Nres; }

                                                for (r = r_from; r < r_to; r++)         /* for each residue in the strand */
                                                {
                                                        dist2 = calc_dist2(xi, yi, zi,
                                                                        x_hcp1[3*r+0], x_hcp1[3*r+1], x_hcp1[3*r+2]);
                                                        if (dist2 > dist2_hcp1)
                                                        {
                                                                for (q = 0; q < hcp; q++)
                                                                {
                                                                        dist2 = calc_dist2(xi, yi, zi, q_hcp1[4*(r*hcp+q)], q_hcp1[4*(r*hcp+q)+1], q_hcp1[4*(r*hcp+q)+2]);
                                                                        sumi -= ebg_hcp_approx_reff(r, dist2, r_hcp1[r * hcp + q]);
                                                                }
                                                        }
                                                        else
                                                        {
                                                                a_from = prm->Ipres[r];
                                                                if (r + 1 < prm->Nres) { a_to = prm->Ipres[r + 1]; }
                                                                else { a_to = prm->Natom; }

                                                                for (a = a_from; a < a_to; a++)     /* for each atom in the residue */
                                                                {
                                                                        if (a != i)   /* not same atom */
                                                                        {
                                                                                r2 = calc_dist2(xi, yi, zi, x[dim*a], x[dim*a+1], x[dim*a+2]);

                                                                                dij1i = 1.0 / sqrt(r2);
                                                                                dij = r2 * dij1i;
                                                                                sj = fs[a] * (rborn[a] - gboffset);
                                                                                sj2 = sj * sj;

                                                                                /*
                                                                                 * ---following are from the Appendix of Schaefer and Froemmel,
                                                                                 * JMB 216:1045-1066, 1990;  Taylor series expansion for d>>s
                                                                                 * is by Andreas Svrcek-Seiler;
                                                                                 */

                                                                                if (dij > 4.0 * sj) {
                                                                                        dij2i = dij1i * dij1i;
                                                                                        tmpsd = sj2 * dij2i;
                                                                                        dumbo = TA + tmpsd * (TB + tmpsd * (TC +
                                                                                                                tmpsd * (TD + tmpsd * TDD)));
                                                                                        sumi -= sj * tmpsd * dij2i * dumbo;

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

                                                                        }   /* end if not same atom */
                                                                }   /* end for each atom a */
                                                        }   /* end else in residue threshold dist */
                                                }   /* end for each residue r */
                                        }   /* end else in strand threshod dist */
                                }   /* end for each strand s */
                        }   /* end else in complex threshold dist */
                }   /* end for each complex */


                if (gb == 1) {

                        /* "standard" (HCT) effective radii:  */
                        reff_hcp0[i] = 1.0 / (ri1i + sumi);
                        if (reff_hcp0[i] < 0.0)
                                reff_hcp0[i] = 30.0;

                } else {

                        /* "gbao" formulas:  */
                        psi[i] = -ri * sumi;
                        reff_hcp0[i] = 1.0 / (ri1i - tanh((gbalpha[i] - gbbeta[i] * psi[i] +
                                                        gbgamma[i] * psi[i] * psi[i]) * psi[i]) / rborn[i]);
                }
                   //  reff_hcp0[i] = 0;
        }

        /* The MPI synchronization is accomplished via reduction of the reff array. */
#if defined(MPI) || defined(SCALAPACK)
        t1 = seconds();

        reductarr = vector(0, prm->Natom);

        ierror = MPI_Allreduce(reff_hcp0, reductarr, prm->Natom,
                        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (ierror != MPI_SUCCESS) {
                fprintf(nabout,
                                "Error in egb reff reduction, error = %d  mytaskid = %d\n",
                                ierror, mytaskid);
        }
        for (i = 0; i < prm->Natom; i++) {
                reff_hcp0[i] = reductarr[i];
        }

        free_vector(reductarr, 0, prm->Natom);

        /* Update the reduction time. */
        t2 = seconds();
        *treduceegb += t2 - t1;
        t1 = t2;

#endif

        return(0);

}  /* end function egb_hcp_atom_reff() */


/***********************************************************************
 *                        calc_sumi_hcp()
 *
 * Calculate contribution of each component to effective Born radii
 *
 * Calling parameters are as follows:
 *
 * ri            - component radius of i             (input)
 * r2            - dist^2 between i and j            (input)
 * sj            - component radius for j            (input)
 * sumi          - contribution to Born radius       (output)
 *
 ************************************************************************/

REAL_T calc_sumi_hcp(REAL_T ri, REAL_T r2, REAL_T sj)
{

        REAL_T ri1i, dij1i, dij, dij2i, sj2, tmpsd, dumbo, theta, uij;
        REAL_T sumi = 0;

        //   ri *= PACK_FACT;
        //   sj *= PACK_FACT;

        ri1i = 1. / ri;
        dij1i = 1.0 / sqrt(r2);
        dij = r2 * dij1i;

        sj2 = sj * sj;

        if (dij > 4.0 * sj) {

                dij2i = dij1i * dij1i;
                tmpsd = sj2 * dij2i;
                dumbo = TA + tmpsd * (TB + tmpsd * (TC + tmpsd * (TD + tmpsd * TDD)));
                sumi = sj * tmpsd * dij2i * dumbo;

        } else if (dij > ri + sj) {

                sumi = 0.5 * (sj / (r2 - sj2) +
                                0.5 * dij1i * log((dij - sj) / (dij + sj)));

        } else if (dij > fabs(ri - sj)) {

                theta = 0.5 * ri1i * dij1i * (r2 + ri * ri - sj2);
                uij = 1. / (dij + sj);
                sumi = 0.25 * (ri1i * (2. - theta) - uij + dij1i * log(ri * uij));

        } else if (ri < sj) {

                sumi = 0.5 * (sj / (r2 - sj2) + 2. * ri1i +
                                0.5 * dij1i * log((sj - dij) / (sj + dij)));

        }

        return(sumi * LAMBDA);

}

/***********************************************************************
 *                        egb_hcp_reff_new()
 *
 * Calculate effective Born radii for hcp components
 *
 * Calling parameters are as follows:
 *
 * x_hcpx        - HCP component coordinates         (input)
 * r_hcpx        - HCP component radius              (input)
 * reff_hcpx     - HCP component Born radius         (output)
 *
 ************************************************************************/


void egb_hcp_reff_new(REAL_T * x_hcp1, REAL_T * q_hcp1, REAL_T * r_hcp1, REAL_T * reff_hcp1,
                REAL_T * x_hcp2, REAL_T * q_hcp2, REAL_T * r_hcp2, REAL_T * reff_hcp2,
                REAL_T * x_hcp3, REAL_T * q_hcp3, REAL_T * r_hcp3, REAL_T * reff_hcp3,
                REAL_T hcp_h1, REAL_T hcp_h2, REAL_T hcp_h3, PARMSTRUCT_T * prm)
{
        int i, qi, qj, c, s, r, s_from, s_to, r_from, r_to;
        REAL_T xi, yi, zi;
        REAL_T ri, ri1i, sumi, r2, psi;
        REAL_T dist2_hcp2, dist2_hcp3;

        /* square of threshold distances */
        dist2_hcp2 = hcp_h2*hcp_h2;   
        dist2_hcp3 = hcp_h3*hcp_h3;

        /* calculate residue Born radii */
        for (i = 0; i < prm->Nres; i++)
        {

                for (qi = 0; qi < hcp; qi++)
                {
                        xi = q_hcp1[4*(i * hcp + qi)];
                        yi = q_hcp1[4*(i * hcp + qi) + 1];
                        zi = q_hcp1[4*(i * hcp + qi) + 2];

                        ri = r_hcp1[i * hcp + qi];
                        ri1i = 1. / ri;
                        sumi = 0.0;


                        for (c = 0; c < prm->Ncomplex; c++)
                        {
                                r2 = calc_dist2(xi, yi, zi, 
                                                x_hcp3[3*c], x_hcp3[3*c+1], x_hcp3[3*c+2]);

                                if (r2 > dist2_hcp3)
                                {
                                        for (qj = 0; qj < hcp; qj++)
                                        {
                                                r2 = calc_dist2(xi, yi, zi, 
                                                                q_hcp3[4*(c*hcp+qj)], q_hcp3[4*(c*hcp+qj)+1], q_hcp3[4*(c*hcp+qj)+2]);
                                                sumi += calc_sumi_hcp(ri, r2, r_hcp3[c*hcp+qj]);
                                        }
                                }
                                else
                                { 
                                        s_from = prm->Ipcomplex[c] - 1;
                                        if (c < prm->Ncomplex - 1) { s_to = prm->Ipcomplex[c + 1] - 1; }
                                        else { s_to = prm->Nstrand; }
                                        for (s = s_from; s < s_to; s++) /* for each strand */
                                        {
                                                r2 = calc_dist2(xi, yi, zi, 
                                                                x_hcp2[3*s], x_hcp2[3*s+1], x_hcp2[3*s+2]);

                                                if (r2 > dist2_hcp2)
                                                {
                                                        for (qj = 0; qj < hcp; qj++)
                                                        {
                                                                r2 = calc_dist2(xi, yi, zi, 
                                                                                q_hcp2[4*(s*hcp+qj)], q_hcp2[4*(s*hcp+qj)+1], q_hcp2[4*(s*hcp+qj)+2]);
                                                                sumi += calc_sumi_hcp(ri, r2, r_hcp2[s*hcp+qj]);
                                                        }
                                                }
                                                else
                                                { 
                                                        r_from = prm->Ipstrand[s] - 1;
                                                        if (s < prm->Nstrand - 1) { r_to = prm->Ipstrand[s + 1] - 1; }
                                                        else { r_to = prm->Nres; }
                                                        for (r = r_from; r < r_to; r++)         /* for each residue */
                                                        {

                                                                for (qj = 0; qj < hcp; qj++)
                                                                {
                                                                        if ((i != r) || (qi != qj))
                                                                        {
                                                                                r2 = calc_dist2(xi, yi, zi, 
                                                                                                q_hcp1[4*(r*hcp+qj)], q_hcp1[4*(r*hcp+qj)+1], q_hcp1[4*(r*hcp+qj)+2]);
                                                                                sumi += calc_sumi_hcp(ri, r2, r_hcp1[r*hcp+qj]);
                                                                        }
                                                                }   /* end if not same residue */
                                                        }   /* end for each residue r */
                                                }   /* end if inside strand threshold dist */
                                        }   /* end for each strand s */
                                }   /* end if inside complex threshold distance */
                        }   /* end for each complex c */

                        if (gb == 1) {

                                /* "standard" (HCT) effective radii:  */
                                reff_hcp1[i*hcp+qi] = 1.0 / (ri1i + sumi);
                                if (reff_hcp1[i*hcp+qi] < 0.0)
                                        reff_hcp1[i*hcp+qi] = 30.0;

                        } else {

                                /* "gbao" formulas:  */
                                psi = -ri * sumi;
                                reff_hcp1[i*hcp+qi] = 1.0 / (ri1i - tanh((gbalpha[i] - gbbeta[i] * psi +
                                                                gbgamma[i] * psi * psi) * psi) / ri);
                        }

                       // if (gb_debug)
                        //        printf("Residue %d (%d)\t%15.7f\t%15.7f\n", i+1, qi + 1, r_hcp1[i*hcp+qi], reff_hcp1[i*hcp+qi]); 
                }   /* end for each hcp charge qi */
        }   /* end for each residue i */

        /* calculate strand Born radii */
        for (i = 0; i < prm->Nstrand; i++)
        {

                for (qi = 0; qi < hcp; qi++)
                {
                        xi = q_hcp2[4*(i * hcp + qi)];
                        yi = q_hcp2[4*(i * hcp + qi) + 1];
                        zi = q_hcp2[4*(i * hcp + qi) + 2];

                        ri = r_hcp2[i * hcp + qi];
                        ri1i = 1. / ri;
                        sumi = 0.0;

                        for (c = 0; c < prm->Ncomplex; c++)
                        {
                                r2 = calc_dist2(xi, yi, zi, 
                                                x_hcp3[3*c], x_hcp3[3*c+1], x_hcp3[3*c+2]);

                                if (r2 > dist2_hcp3)
                                {
                                        for (qj = 0; qj < hcp; qj++)
                                        {
                                                r2 = calc_dist2(xi, yi, zi, 
                                                                q_hcp3[4*(c*hcp+qj)], q_hcp3[4*(c*hcp+qj)+1], q_hcp3[4*(c*hcp+qj)+2]);
                                                sumi += calc_sumi_hcp(ri, r2, r_hcp3[c*hcp+qj]);
                                        }
                                }
                                else
                                { 
                                        s_from = prm->Ipcomplex[c] - 1;
                                        if (c < prm->Ncomplex - 1) { s_to = prm->Ipcomplex[c + 1] - 1; }
                                        else { s_to = prm->Nstrand; }
                                        for (s = s_from; s < s_to; s++) /* for each strand */
                                        {
                                                for (qj = 0; qj < hcp; qj++)
                                                {
                                                        if ((i != s) && (qi != qj))   /* not same strand */
                                                        {
                                                                r2 = calc_dist2(xi, yi, zi, 
                                                                                q_hcp2[4*(s*hcp+qj)], q_hcp2[4*(s*hcp+qj)+1], q_hcp2[4*(s*hcp+qj)+2]);
                                                                sumi += calc_sumi_hcp(ri, r2, r_hcp2[s*hcp+qj]);
                                                        }
                                                }
                                        }   /* end for each strand s */
                                }   /* end if inside complex threshold distance */
                        }   /* end for each complex c */

                        if (gb == 1) {

                                /* "standard" (HCT) effective radii:  */
                                reff_hcp2[i*hcp+qi] = 1.0 / (ri1i + sumi);
                                if (reff_hcp2[i*hcp+qi] < 0.0)
                                        reff_hcp2[i*hcp+qi] = 30.0;

                        } else {

                                /* "gbao" formulas:  */
                                psi = -ri * sumi;
                                reff_hcp2[i*hcp+qi] = 1.0 / (ri1i - tanh((gbalpha[i] - gbbeta[i] * psi +
                                                                gbgamma[i] * psi * psi) * psi) / ri);
                        }

                       // if (gb_debug)
                         //       printf("Strand %d (%d)\t%15.7f\t%15.7f\n", i + 1, qi + 1, r_hcp2[i*hcp+qi], reff_hcp2[i*hcp+qi]);
                }   /* end for each hcp charge qi */
        }   /* end for each strand i */

        /* calculate complex Born radii */
        for (i = 0; i < prm->Ncomplex; i++)
        {

                for (qi = 0; qi < hcp; qi++)
                {
                        xi = q_hcp3[4*(i * hcp + qi)];
                        yi = q_hcp3[4*(i * hcp + qi) + 1];
                        zi = q_hcp3[4*(i * hcp + qi) + 2];

                        ri = r_hcp3[i * hcp + qi];
                        ri1i = 1. / ri;
                        sumi = 0.0;

                        for (c = 0; c < prm->Ncomplex; c++)
                        {
                                for (qj = 0; qj < hcp; qj++)
                                {
                                        if ((i != c) && (qi != qj))  /* not same complex */
                                        {
                                                r2 = calc_dist2(xi, yi, zi, 
                                                                q_hcp3[4*(c*hcp+qj)], q_hcp3[4*(c*hcp+qj)+1], q_hcp3[4*(c*hcp+qj)+2]);
                                                sumi += calc_sumi_hcp(ri, r2, r_hcp3[c*hcp+qj]);
                                        }
                                }
                        }   /* end for each complex c */

                        if (gb == 1) {

                                /* "standard" (HCT) effective radii:  */
                                reff_hcp3[i*hcp+qi] = 1.0 / (ri1i + sumi);
                                if (reff_hcp3[i*hcp+qi] < 0.0)
                                        reff_hcp3[i*hcp+qi] = 30.0;

                        } else {

                                /* "gbao" formulas:  */
                                psi = -ri * sumi;
                                reff_hcp3[i*hcp+qi] = 1.0 / (ri1i - tanh((gbalpha[i] - gbbeta[i] * psi +
                                                                gbgamma[i] * psi * psi) * psi) / ri);
                        }

                       // if (gb_debug)
                         //       printf("Complex %d (%d)\t%15.7f\t%15.7f\n", i + 1, qi + 1, r_hcp3[i*hcp+qi], reff_hcp3[i*hcp+qi]);
                }   /* end for each hcp charge qi */
        }   /* end for each complex i */

}




/***********************************************************************
 *                         egb_hcp_reff()
 *
 * Calculate effective Born radii for residues
 *
 * Calling parameters are as follows:
 * prm       - parameter/topology structure           (input)
 * reff_hcp0 - atom effective born radii              (input)
 * reff_hcpx - (x=1/2) residue/chain born radii       (update)
 *
 ************************************************************************/
int egb_hcp_reff(PARMSTRUCT_T * prm,
                REAL_T * reff_hcp0, REAL_T * reff_hcp1, REAL_T * reff_hcp2, REAL_T * reff_hcp3)
{

        int c, s, r, a, m;      /* complex, strand, residue, atom index */
        int a_from, a_to;       /* for atoms within a residue */
        int r_from, r_to;	   /* for residues within a strand */
        REAL_T q, hcp1_q[3], hcp2_q[3], hcp3_q[3];


        for (c = 0; c < prm->Ncomplex; c++)
        {
                for (m = 0; m < hcp; m++)
                {
                        reff_hcp3[c * hcp + m] = 0;
                        hcp3_q[m] = 0;
                }

                for (s = 0; s < prm->Nstrand; s++)
                {
                        for (m = 0; m < hcp; m++)
                        {
                                reff_hcp2[s * hcp + m] = 0;
                                hcp2_q[m] = 0;
                        }

                        r_from = prm->Ipstrand[s] - 1;
                        if (s < prm->Nstrand - 1) { r_to = prm->Ipstrand[s + 1] - 1; }
                        else { r_to = prm->Nres; }
                        for (r = r_from; r < r_to; r++)
                        {
                                for (m = 0; m < hcp; m++)
                                {
                                        reff_hcp1[r * hcp + m] = 0;
                                        hcp1_q[m] = 0;
                                }

                                a_from = prm->Ipres[r];
                                if (r + 1 < prm->Nres) {a_to = prm->Ipres[r + 1]; }
                                else                   {a_to = prm->Natom; }
                                for (a = a_from; a < a_to; a++)
                                {
                                        q = prm->Charges[a];
                                        m = q2idx(q);

                                       // if (HCP_RAD == 2)    /* archontis-simonson */
                                       // {
                                       //         reff_hcp1[r * hcp + m] += (q * q / reff_hcp0[a]);
                                       //         hcp1_q[m] += (q * q);
                                      //  }
                                       // else			/* anandakrishnan-onufriev */
                                      //  {
                                                reff_hcp1[r * hcp + m] += (q / sqrt(reff_hcp0[a]));
                                                hcp1_q[m] += q;
                                      //  }
                                }

                                for (m = 0; m < hcp; m++)
                                {
                                        /* rollup totals to next higher level */
                                        reff_hcp2[s * hcp + m] += reff_hcp1[r * hcp + m];
                                        hcp2_q[m] += hcp1_q[m];
                                        /* compute residue level effective born radii */
                                      //  if (HCP_RAD == 2)
                                     //   {
                                      //          reff_hcp1[r * hcp + m] = hcp1_q[m] / reff_hcp1[r * hcp +  m];
                                      //  }
                                      //  else
                                      //  {
                                                if (fabs(reff_hcp1[r * hcp +m]) < 1)
                                                {
                                                        reff_hcp1[r * hcp + m] = 2.3;   /* for partial charge sum << 1 */
                                                }
                                                else
                                                {
                                                        reff_hcp1[r * hcp + m] = (hcp1_q[m] * hcp1_q[m]) / 
                                                                (reff_hcp1[r * hcp +  m] * reff_hcp1[r * hcp + m]);
                                                }
                                       // }
                                      //  if (gb_debug)
                                        //        printf("Residue %d (%d)\t%15.7f\n", r+1, m+1, reff_hcp1[r * hcp + m]);
                                }
                        }

                        for (m = 0; m < hcp; m++)
                        {
                                /* rollup totals to next higher level */
                                reff_hcp3[c * hcp + m] += reff_hcp2[s * hcp + m];
                                hcp3_q[m] += hcp2_q[m];
                                /* compute strand level effective born radii */
                                if (HCP_RAD == 2)
                                        reff_hcp2[s * hcp + m] = hcp2_q[m] / reff_hcp2[s * hcp +  m];
                                else
                                        reff_hcp2[s * hcp + m] = (hcp2_q[m] * hcp2_q[m]) / 
                                                (reff_hcp2[s * hcp + m] * reff_hcp2[s * hcp + m]);
                                //if (gb_debug)
                                  //      printf("Strand %d (%d)\t%15.7f\n", s+1, m+1, reff_hcp2[s * hcp + m]);
                        }
                }

                for (m = 0; m < hcp; m++)
                {
                        /* compute complex level effective born radii */
                      //  if (HCP_RAD == 2)
                      //          reff_hcp3[c * hcp + m] = hcp3_q[m] / reff_hcp3[c * hcp +  m];
                     //   else
                                reff_hcp3[c * hcp + m] = (hcp3_q[m] * hcp3_q[m]) / 
                                        (reff_hcp3[c * hcp + m] * reff_hcp3[c * hcp + m]);
                       // if (gb_debug)
                         //       printf("complex %d (%d)\t%15.7f\n", c+1, m+1, reff_hcp3[c * hcp + m]);
                }
        }	
        return(0);
}

/***********************************************************************
 *                         egb_hcpo_reff_avg()
 * Izadi, Anandakrishnan and Onufriev, JCTC 2016 (in prep)
 * 
 * Calculate effective Born radii for residues, modified for new hcp (HCPO)
 * Calling parameters are as follows:
 * prm       - parameter/topology structure           (input)
 * reff_hcp0 - atom effective born radii              (input)
 * reff_hcpx - (x=1/2) residue/chain born radii       (update)
 * q_hcpi     - (i=1,2,3) approximate charges/coordinates for residue/strand/complex (input)  
 * x_hcpi     - (i=1,2,3) coordinates of center of geometry for residue/strand/complex (input)
 ************************************************************************/
int egb_hcpo_reff_avg(PARMSTRUCT_T * prm,
                REAL_T * reff_hcp0, REAL_T * reff_hcp1, REAL_T * reff_hcp2, REAL_T * reff_hcp3, 
                 REAL_T * q_hcp1, REAL_T * q_hcp2, REAL_T * q_hcp3,
                 REAL_T * x_hcp1, REAL_T * x_hcp2, REAL_T * x_hcp3, REAL_T * x)
{

        int c, s, r, a, m, n;      /* complex, strand, residue, atom index */
        int a_from, a_to;       /* for atoms within a residue */
        int r_from, r_to;	   /* for residues within a strand */
        REAL_T q, hcp1_q[3], hcp2_q[3], hcp3_q[3];
        REAL_T dummy=0, rdis2l1, ql1, rdis2l2, ql2, rdis2l3, ql3; 
        int ja; 


        for (c = 0; c < prm->Ncomplex; c++)
        {
                for (m = 0; m < hcp; m++)
                {
                        reff_hcp3[c * hcp + m] = 0;
                        hcp3_q[m] = 0;
                }

                for (s = 0; s < prm->Nstrand; s++)
                {
                        for (m = 0; m < hcp; m++)
                        {
                                reff_hcp2[s * hcp + m] = 0;
                                hcp2_q[m] = 0;
                        }

                        r_from = prm->Ipstrand[s] - 1;
                        if (s < prm->Nstrand - 1) { r_to = prm->Ipstrand[s + 1] - 1; }
                        else { r_to = prm->Nres; }
                        for (r = r_from; r < r_to; r++)
                        {
                                for (m = 0; m < hcp; m++)
                                {
                                        reff_hcp1[r * hcp + m] = 0;
                                        hcp1_q[m] = 0;
                                }

                                a_from = prm->Ipres[r];
                                if (r + 1 < prm->Nres) {a_to = prm->Ipres[r + 1]; }
                                else                   {a_to = prm->Natom; }
                                for (a = a_from; a < a_to; a++)
                                {
                                        q = (prm->Charges[a]);
                                         m = q2idx(q);  
                                 // the avergaing in performed over positive and negative charges separately, 
                                 // there is a sum over total charges and apparantly 
                                 // summing up negative and positive charges makes the total 
                                 // charge close to zero (for neutral cases) and create issues, 
                                 // so +/- charges are averaged separtely.  
                                        ja = 3 * a;     
                                        rdis2l1=(x[ja+0]-q_hcp1[4*(r*hcp+m)+0])*(x[ja+0]-q_hcp1[4*(r*hcp+m)+0])+
                                                (x[ja+1]-q_hcp1[4*(r*hcp+m)+1])*(x[ja+1]-q_hcp1[4*(r*hcp+m)+1])+
                                                (x[ja+2]-q_hcp1[4*(r*hcp+m)+2])*(x[ja+2]-q_hcp1[4*(r*hcp+m)+2]); 
                                        ql1 =fabs(q) / sqrt(rdis2l1); // scale charges for new averaging method
                                        reff_hcp1[r * hcp + m] += (ql1 / sqrt(reff_hcp0[a]));
                                        hcp1_q[m] += ql1;        
                                } // end loob over a
                                   
                                for (m = 0; m < hcp; m++)
                                {
                                      if (fabs(reff_hcp1[r * hcp +m]) < 1)
                                      {
                                        reff_hcp1[r * hcp + m] = 2.3;   /* for partial charge sum << 1 */
                                      }
                                      else
                                      {
                                        reff_hcp1[r * hcp + m] = (hcp1_q[m] * hcp1_q[m]) / 
                                                      (reff_hcp1[r * hcp +  m] * reff_hcp1[r * hcp + m]);
                                      }
                                 }
                                // since averaging is done for positive and negative charges 
                                // independently, we should make sure reff is assigned to the right charges
                                if (q_hcp1[4*(r*hcp+0)+3] > q_hcp1[4*(r*hcp+1)+3]) 
                                {
                                   dummy=reff_hcp1[r * hcp + 0];
                                   reff_hcp1[r * hcp + 0] = reff_hcp1[r * hcp + 1];
                                   reff_hcp1[r * hcp + 1] = dummy;
                                } 
			                    for (m = 0; m < hcp; m++)
			                 	{	
                                   if (q_hcp1[4*(r*hcp+m)+3]>0){
                                     n = 1;
                                   } else {
                                     n = 0;
                                   }  
                                   rdis2l2=(q_hcp1[4*(r*hcp+m)+0]-q_hcp2[4*(s*hcp+n)+0])*(q_hcp1[4*(r*hcp+m)+0]-q_hcp2[4*(s*hcp+n)+0])+
                                           (q_hcp1[4*(r*hcp+m)+1]-q_hcp2[4*(s*hcp+n)+1])*(q_hcp1[4*(r*hcp+m)+1]-q_hcp2[4*(s*hcp+n)+1])+
                                           (q_hcp1[4*(r*hcp+m)+2]-q_hcp2[4*(s*hcp+n)+2])*(q_hcp1[4*(r*hcp+m)+2]-q_hcp2[4*(s*hcp+n)+2]); 
                                   ql2 = q_hcp1[4*(r*hcp+m)+3] / sqrt(rdis2l2);
			                       reff_hcp2[s * hcp + n] += (ql2 / sqrt(reff_hcp1[r * hcp + m]));
			                       hcp2_q[n] += ql2;         
			                    }	   
                        } // end loop over r
                      
                        for (m = 0; m < hcp; m++)
                        {
                            if (fabs(reff_hcp2[s * hcp +m]) < 1)
                            {
                               reff_hcp2[s * hcp + m] = 2.3;   /* for partial charge sum << 1 */
                            }
                            else
                            {
                               reff_hcp2[s * hcp + m] = (hcp2_q[m] * hcp2_q[m]) / 
                                         (reff_hcp2[s * hcp + m] * reff_hcp2[s * hcp + m]);
                            }
                        }
                     /*   if (q_hcp2[4*(s*hcp+0)+3] > q_hcp2[4*(s*hcp+1)+3]) 
                        {
                          dummy=reff_hcp2[s * hcp + 0];
                          reff_hcp2[s * hcp + 0] = reff_hcp2[s * hcp + 1];
                          reff_hcp2[s * hcp + 1] = dummy;
                        } */
			            for (m = 0; m < hcp; m++)
			            {	
                          if (q_hcp2[4*(s*hcp+m)+3]>0){
                             n = 1;
                          } else {
                             n = 0;
                          }
                           rdis2l3=(q_hcp2[4*(s*hcp+m)+0]-q_hcp3[4*(c*hcp+n)+0])*(q_hcp2[4*(s*hcp+m)+0]-q_hcp3[4*(c*hcp+n)+0])+
                                   (q_hcp2[4*(s*hcp+m)+1]-q_hcp3[4*(c*hcp+n)+1])*(q_hcp2[4*(s*hcp+m)+1]-q_hcp3[4*(c*hcp+n)+1])+
                                   (q_hcp2[4*(s*hcp+m)+2]-q_hcp3[4*(c*hcp+n)+2])*(q_hcp2[4*(s*hcp+m)+2]-q_hcp3[4*(c*hcp+n)+2]); 
                           ql3 = q_hcp2[4*(s*hcp+m)+3] / sqrt(rdis2l3);
			               reff_hcp3[c * hcp + n] += (ql3 / sqrt(reff_hcp2[s * hcp + m]));
			               hcp3_q[n] += ql3;         
                        }
                          
                }  // ending the loop over strands (s)

                for (m = 0; m < hcp; m++)
                {
                  if (fabs(reff_hcp3[c * hcp +m]) < 1)
                  {  
                       reff_hcp3[c * hcp + m] = 2.3;   /* for partial charge sum << 1 */
                  }
                  else
                  {  
                       reff_hcp3[c * hcp + m] = (hcp3_q[m] * hcp3_q[m]) / 
                                     (reff_hcp3[c * hcp + m] * reff_hcp3[c * hcp + m]);
                  }
                }
              //  if (q_hcp3[4*(c*hcp+0)+3] > q_hcp3[4*(c*hcp+1)+3]) 
              //  {
              //     dummy=reff_hcp3[c * hcp + 0];
               //    reff_hcp3[c * hcp + 0] = reff_hcp3[c * hcp + 1];
               //    reff_hcp3[c * hcp + 1] = dummy;
               // }
               //printf("newavg2newavg2"); 
               
        }	
        return(0);
}

/***********************************************************************
 *                         HCP_GBSA()
 *
 * Calculate effective Born radii for residues
 *
 * Calling parameters are as follows:
 * npairs_hcpx - (x=0/1/2) number of pairs in each pairlist (input)
 * pairs_hcpx  - (x=0/1/2) list of atom/residue/chains      (input)
 * Various parameters from egb_nbond_hcp()		    (input)

 ************************************************************************/ 

void hcp_gbsa(REAL_T *x, REAL_T *f, PARMSTRUCT_T *prm, REAL_T *esurf)
{
        /*
         * Calculate either the volume or the surface area component of the AGB
         * nonpolar hydration free energy.  Use Eqs 2-15 of Gallicchio and Levy,
         * J. Comput. Chem. 25: 479 (2004).
         *
         */

        /* LCPO stuff follows */
        int count, count2, icount;
        REAL_T si, sumAij, sumAjk, sumAijAjk, sumdAijddijdxi;
        REAL_T sumdAijddijdyi, sumdAijddijdzi, sumdAijddijdxiAjk;
        REAL_T sumdAijddijdyiAjk, sumdAijddijdziAjk, rij, tmpaij, Aij, dAijddij;
        REAL_T dAijddijdxj, dAijddijdyj, dAijddijdzj;
        REAL_T sumdAjkddjkdxj, sumdAjkddjkdyj, sumdAjkddjkdzj, p3p4Aij;
        REAL_T xk, yk, zk, rjk2, djk1i, rjk, vdw2dif, tmpajk, Ajk, sumAjk2, dAjkddjk;
        REAL_T dAjkddjkdxj, dAjkddjkdyj, dAjkddjkdzj, lastxj, lastyj, lastzj;
        REAL_T dAidxj, dAidyj, dAidzj, Ai, dAidxi, dAidyi, dAidzi;
        REAL_T totsasa;

        int i, j, k;
        REAL_T xi, yi, zi, xij, yij, zij, xj, yj, zj;
        REAL_T dij1i;

        REAL_T dij;
        REAL_T r2;

        int c, s, r, a;      /* complex, strand, residue, atom index */
        int a_from, a_to;    /* for atoms within a residue */
        int r_from, r_to;    /* for residues within a strand */
        int s_from, s_to;    /* for strands within a complex */
        REAL_T dist2, dist2_hcp1, dist2_hcp2, dist2_hcp3;

        /* square of threshold distances */
        dist2_hcp1 = hcp_h1 * hcp_h1;
        dist2_hcp2 = hcp_h2 * hcp_h2;
        dist2_hcp3 = hcp_h3 * hcp_h3;

        *esurf = 0.0;

        /*
         * Main LCPO stuff follows.
         */

        /* Loop over all atoms i. */

        count = 0;
        for (i = 0; i < prm->Natom; i++) {

#if defined(MPI) 
                if (!myroc(i, blocksize, numtasks, mytaskid))
                        continue;
#endif

                xi = x[3 * i];
                yi = x[3 * i + 1];
                zi = x[3 * i + 2];

                for (c=0; c < prm->Ncomplex; c++)   /* for each complex in the molecule */
                {
                        dist2 = calc_dist2(xi, yi, zi,
                                        x_hcp3[3*c+0], x_hcp3[3*c+1], x_hcp3[3*c+2]);
                        if (dist2 < dist2_hcp3)
                        {
                                s_from = prm->Ipcomplex[c] - 1;
                                if (c + 1 < prm->Ncomplex) { s_to = prm->Ipcomplex[c + 1] - 1; }
                                else { s_to = prm->Nstrand; }

                                for (s=s_from; s < s_to; s++)   /* for each strand in the complex */
                                {
                                        dist2 = calc_dist2(xi, yi, zi,
                                                        x_hcp2[3*s+0], x_hcp2[3*s+1], x_hcp2[3*s+2]);
                                        if (dist2 < dist2_hcp2)
                                        {
                                                r_from = prm->Ipstrand[s] - 1;
                                                if (s + 1 < prm->Nstrand) { r_to = prm->Ipstrand[s + 1] - 1; }
                                                else { r_to = prm->Nres; }

                                                for (r = r_from; r < r_to; r++)         /* for each residue in the strand */
                                                {
                                                        dist2 = calc_dist2(xi, yi, zi,
                                                                        x_hcp1[3*r+0], x_hcp1[3*r+1], x_hcp1[3*r+2]);
                                                        if (dist2 < dist2_hcp1)
                                                        {
                                                                a_from = prm->Ipres[r];
                                                                if (r + 1 < prm->Nres) { a_to = prm->Ipres[r + 1]; }
                                                                else { a_to = prm->Natom; }

                                                                for (a = a_from; a < a_to; a++)     /* for each atom in the residue */
                                                                {
                                                                        if (a != i)   /* not same atom */
                                                                        {
                                                                                xij = xi - x[3 * a];
                                                                                yij = yi - x[3 * a + 1];
                                                                                zij = zi - x[3 * a + 2];
                                                                                r2 = xij * xij + yij * yij + zij * zij;

                                                                                dij1i = 1.0 / sqrt(r2);
                                                                                dij = r2 * dij1i;

                                                                                if ((P0[i] + P0[a]) > dij) {
                                                                                        if (P0[i] > 2.5 && P0[a] > 2.5) {

                                                                                                ineighbor[count] = a + 1;
                                                                                                /* this +1 is VERY important */
                                                                                                count += 1;
                                                                                        }
                                                                                }
                                                                        }   /* end if not same atom */
                                                                }   /* end for each atom a */
                                                        }   /* end if in threshold dist h1 */
                                                }   /* end for each residue r */
                                        }   /* end if in threshold dist h2 */
                                }   /* end for each strand s */
                        }   /* end if in threshold dist h3 */
                }   /* end for each complex c */

                ineighbor[count] = 0;
                count = count + 1;
        }

        totsasa = 0.0;
        count = 0;
        for (i = 0; i < prm->Natom; i++) {
#if defined(MPI) 
                if (!myroc(i, blocksize, numtasks, mytaskid))
                        continue;
#endif
                if (ineighbor[count] == 0) {
                        count = count + 1;
                } else {

                        si = PI*4 * P0[i] * P0[i];
                        sumAij = 0.0;
                        sumAjk = 0.0;
                        sumAjk2 = 0.0;
                        sumAijAjk = 0.0;
                        sumdAijddijdxi = 0.0;
                        sumdAijddijdyi = 0.0;
                        sumdAijddijdzi = 0.0;
                        sumdAijddijdxiAjk = 0.0;
                        sumdAijddijdyiAjk = 0.0;
                        sumdAijddijdziAjk = 0.0;

                        icount = count;

L70:j = ineighbor[count] - 1;
    /* mind the -1 , fortran again */
    xi = x[3 * i];
    yi = x[3 * i + 1];
    zi = x[3 * i + 2];

    xj = x[3 * j];
    yj = x[3 * j + 1];
    zj = x[3 * j + 2];

    r2 = (xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) + (zi -
                    zj) *
            (zi - zj);
    dij1i = 1. / sqrt(r2);
    rij = r2 * dij1i;
    tmpaij =
            P0[i] - rij * 0.5 - (P0[i] * P0[i] -
                            P0[j] * P0[j]) * 0.5 * dij1i;
    Aij = PI*2 * P0[i] * tmpaij;

    dAijddij =
            PI * P0[i] * (dij1i * dij1i *
                            (P0[i] * P0[i] - P0[j] * P0[j]) - 1.0);

    dAijddijdxj = dAijddij * (xj - xi) * dij1i;
    dAijddijdyj = dAijddij * (yj - yi) * dij1i;
    dAijddijdzj = dAijddij * (zj - zi) * dij1i;

    sumAij = sumAij + Aij;

    count2 = icount;
    sumAjk2 = 0.0;
    sumdAjkddjkdxj = 0.0;
    sumdAjkddjkdyj = 0.0;
    sumdAjkddjkdzj = 0.0;

    p3p4Aij = -surften * (P3[i] + P4[i] * Aij);

L80:k = ineighbor[count2] - 1;
    /*same as above, -1 ! */
    if (j == k)
            goto L85;

    xk = x[3 * k];
    yk = x[3 * k + 1];
    zk = x[3 * k + 2];

    rjk2 =
            (xj - xk) * (xj - xk) + (yj - yk) * (yj - yk) + (zj -
                            zk) *
            (zj - zk);

    djk1i = 1.0 / sqrt(rjk2);
    rjk = rjk2 * djk1i;

    if (P0[j] + P0[k] > rjk) {
            vdw2dif = P0[j] * P0[j] - P0[k] * P0[k];
            tmpajk = 2.0 * P0[j] - rjk - vdw2dif * djk1i;

            Ajk = PI * P0[j] * tmpajk;

            sumAjk = sumAjk + Ajk;
            sumAjk2 = sumAjk2 + Ajk;

            dAjkddjk =
                    PI * P0[j] * djk1i * (djk1i * djk1i * vdw2dif - 1.0);

            dAjkddjkdxj = dAjkddjk * (xj - xk);
            dAjkddjkdyj = dAjkddjk * (yj - yk);
            dAjkddjkdzj = dAjkddjk * (zj - zk);

            f[3 * k] = f[3 * k] + dAjkddjkdxj * p3p4Aij;
            f[3 * k + 1] = f[3 * k + 1] + dAjkddjkdyj * p3p4Aij;
            f[3 * k + 2] = f[3 * k + 2] + dAjkddjkdzj * p3p4Aij;

            sumdAjkddjkdxj = sumdAjkddjkdxj + dAjkddjkdxj;
            sumdAjkddjkdyj = sumdAjkddjkdyj + dAjkddjkdyj;
            sumdAjkddjkdzj = sumdAjkddjkdzj + dAjkddjkdzj;

    }

L85:count2 = count2 + 1;
    if (ineighbor[count2] != 0) {
            goto L80;
    } else {
            count2 = icount;
    }


    sumAijAjk = sumAijAjk + Aij * sumAjk2;

    sumdAijddijdxi = sumdAijddijdxi - dAijddijdxj;
    sumdAijddijdyi = sumdAijddijdyi - dAijddijdyj;
    sumdAijddijdzi = sumdAijddijdzi - dAijddijdzj;
    sumdAijddijdxiAjk = sumdAijddijdxiAjk - dAijddijdxj * sumAjk2;
    sumdAijddijdyiAjk = sumdAijddijdyiAjk - dAijddijdyj * sumAjk2;
    sumdAijddijdziAjk = sumdAijddijdziAjk - dAijddijdzj * sumAjk2;

    lastxj = dAijddijdxj * sumAjk2 + Aij * sumdAjkddjkdxj;
    lastyj = dAijddijdyj * sumAjk2 + Aij * sumdAjkddjkdyj;
    lastzj = dAijddijdzj * sumAjk2 + Aij * sumdAjkddjkdzj;

    dAidxj = surften * (P2[i] * dAijddijdxj +
                    P3[i] * sumdAjkddjkdxj + P4[i] * lastxj);
    dAidyj =
            surften * (P2[i] * dAijddijdyj +
                            P3[i] * sumdAjkddjkdyj + P4[i] * lastyj);
    dAidzj =
            surften * (P2[i] * dAijddijdzj +
                            P3[i] * sumdAjkddjkdzj + P4[i] * lastzj);

    f[3 * j] = f[3 * j] + dAidxj;
    f[3 * j + 1] = f[3 * j + 1] + dAidyj;
    f[3 * j + 2] = f[3 * j + 2] + dAidzj;

    count = count + 1;
    if (ineighbor[count] != 0) {
            goto L70;
    } else {
            count = count + 1;
    }

    Ai = P1[i] * si + P2[i] * sumAij + P3[i] * sumAjk +
            P4[i] * sumAijAjk;

    dAidxi =
            surften * (P2[i] * sumdAijddijdxi +
                            P4[i] * sumdAijddijdxiAjk);
    dAidyi =
            surften * (P2[i] * sumdAijddijdyi +
                            P4[i] * sumdAijddijdyiAjk);
    dAidzi =
            surften * (P2[i] * sumdAijddijdzi +
                            P4[i] * sumdAijddijdziAjk);

    f[3 * i] = f[3 * i] + dAidxi;
    f[3 * i + 1] = f[3 * i + 1] + dAidyi;
    f[3 * i + 2] = f[3 * i + 2] + dAidzi;

    totsasa = totsasa + Ai;



                }
        }

        *esurf = totsasa * surften;

}


/***********************************************************************
 *                         EGB_CALC_ENERGY_APPROX()
 *
 * Calculate effective Born radii for residues
 *
 * Calling parameters are as follows:
 * npairs_hcpx - (x=0/1/2) number of pairs in each pairlist (input)
 * pairs_hcpx  - (x=0/1/2) list of atom/residue/chains      (input)
 * reff_hcpx                                                (input)  
 * q_hcpx      - (x=1/2) approximate charges                (input)
 * Various parameters from egb_nbond_hcp()		    (input)

 ************************************************************************/ 

void egb_calc_energy_approx(int j, REAL_T * q_hcp, REAL_T * reff_hcp, int hcp, int i,
                REAL_T qi, REAL_T ri, REAL_T *kappa, REAL_T *diel_ext, 
                REAL_T *elec, REAL_T *evdw, REAL_T *epol, REAL_T *sumdeijda, 
                REAL_T xi, REAL_T yi, REAL_T zi, 
                REAL_T *daix, REAL_T *daiy, REAL_T *daiz)
{

        int  j34, k;
        REAL_T xij, yij, zij, r2, rinv, r2inv, qiqj, rj, rb2, efac, fgbi, fgbk, expmkf, dielfac, eel, temp4, temp5, temp6;
        REAL_T de, dedx, dedy, dedz;


        for(k = 0; k < hcp; k++)
        {
                j34 = 4 * (j * hcp + k);

                if (fabs(q_hcp[j34 + 3]) > MIN_Q)
                { 
                        /* Continue computing the non-diagonal energy term. */

                        xij = xi - q_hcp[j34];
                        yij = yi - q_hcp[j34 + 1];
                        zij = zi - q_hcp[j34 + 2];
                        r2 = xij * xij + yij * yij + zij * zij;

                        qiqj = qi * q_hcp[j34 + 3];
                        rj = reff_hcp[j * hcp + k];
                        rb2 = ri * rj;
                        efac = exp(-r2 / (4.0 * rb2));
                        fgbi = 1.0 / sqrt(r2 + rb2 * efac);
                        fgbk = -(*kappa) * KSCALE / fgbi;

                        expmkf = exp(fgbk) / (*diel_ext);
                        dielfac = 1.0 - expmkf;

                        *epol += -qiqj * dielfac * fgbi;

                        temp4 = fgbi * fgbi * fgbi;
                        temp6 = qiqj * temp4 * (dielfac + fgbk * expmkf);

                        de = temp6 * (1.0 - 0.25 * efac);

                        temp5 = 0.5 * efac * temp6 * (rb2 + 0.25 * r2);

                        /*
                         * Compute the contribution of the non-diagonal energy term to the
                         * sum by which the derivatives of Ri and Rj will be multiplied.
                         */

                        sumdeijda[i] += (ri * temp5);

                        rinv = 1. / sqrt(r2);
                        r2inv = rinv * rinv;

                        /*  gas-phase Coulomb energy:  */
                        eel = qiqj * rinv;

                        *elec += eel;
                        de -= eel * r2inv;

                        /*
                         * Sum to the gradient vector the derivatives of Dij that are
                         * computed relative to the cartesian coordinates of atoms i and j.
                         */

                        dedx = de * xij;
                        dedy = de * yij;
                        dedz = de * zij;

                        *daix += dedx;
                        *daiy += dedy;
                        *daiz += dedz;

                } /* end if q > MIN_Q */

        }  /* end for each hcp charge */

}

/***********************************************************************
 *                         EGB_CALC_ENERGY_ATOM()
 *
 * Calculate effective Born radii for residues
 *
 * Calling parameters are as follows:
 * npairs_hcpx - (x=0/1/2) number of pairs in each pairlist (input)
 * pairs_hcpx  - (x=0/1/2) list of atom/residue/chains      (input)
 * reff_hcpx                                                (input)  
 * x           - co-ordinate array                          (input)
 * Various parameters from egb_nbond_hcp()		    (input)
 *
 ************************************************************************/ 

void egb_calc_energy_atom(int j, REAL_T * x, REAL_T * reff_hcp, int *iexw, int hcp, int i, 
                int iaci, REAL_T qi, REAL_T ri, REAL_T *f, REAL_T *kappa, REAL_T *diel_ext, 
                REAL_T *elec, REAL_T *epol, REAL_T *evdw, REAL_T *sumdeijda, 
                REAL_T xi, REAL_T yi, REAL_T zi, REAL_T wi, 
                REAL_T *daix, REAL_T *daiy, REAL_T *daiz, REAL_T *daiw)
{

        int  j34, ic;
        REAL_T xij, yij, zij, r2, rinv, r2inv, r6inv, f6, f12, qiqj, rj, rb2, efac, fgbi, fgbk, expmkf, dielfac, eel;
        REAL_T temp4, temp5, temp6, de, dedx, dedy, dedz;


        j34 = dim * j;

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

        qiqj = qi * prm->Charges[j];
        rj = reff_hcp[j];
        rb2 = ri * rj;
        efac = exp(-r2 / (4.0 * rb2));
        fgbi = 1.0 / sqrt(r2 + rb2 * efac);
        fgbk = -(*kappa) * KSCALE / fgbi;

        expmkf = exp(fgbk) / (*diel_ext);
        dielfac = 1.0 - expmkf;
        *epol += -qiqj * dielfac * fgbi;

        temp4 = fgbi * fgbi * fgbi;
        temp6 = qiqj * temp4 * (dielfac + fgbk * expmkf);
        de = temp6 * (1.0 - 0.25 * efac);

        temp5 = 0.5 * efac * temp6 * (rb2 + 0.25 * r2);

        /*
         * Compute the contribution of the non-diagonal energy term to the
         * sum by which the derivatives of Ri and Rj will be multiplied.
         */

        sumdeijda[i] += ri * temp5;

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
         *
         * The de term contains one more factor of Dij in the denominator
         * so that terms such as dedx do not need to include 1/Dij. 
         */

        if (iexw[j] != i) {

                rinv = 1. / sqrt(r2);
                r2inv = rinv * rinv;
                /*  gas-phase Coulomb energy:  */

                eel = qiqj * rinv;
                *elec += eel;
                de -= eel * r2inv;

                /* Lennard-Jones energy:   */

                ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
                if (ic >= 0) {
                        r6inv = r2inv * r2inv * r2inv;
                        f6 = prm->Cn2[ic] * r6inv;
                        f12 = prm->Cn1[ic] * r6inv * r6inv;
                        /* elect only */
                        *evdw += f12 - f6;
                        de -= (12. * f12 - 6. * f6) * r2inv;
                        /* */
                }
        }

        /*
         * Sum to the gradient vector the derivatives of Dij that are
         * computed relative to the cartesian coordinates of atoms i and j.
         */


        dedx = de * xij;
        dedy = de * yij;
        dedz = de * zij;

        *daix += dedx;
        *daiy += dedy;
        *daiz += dedz;

}

/***********************************************************************
 *                         EGB_CALC_DERIV_APPROX()
 *
 * Calculate effective Born radii for residues
 *
 * Calling parameters are as follows:
 * npairs_hcpx - (x=0/1/2) number of pairs in each pairlist (input)
 * pairs_hcpx  - (x=0/1/2) list of atom/residue/chains      (input)
 * r_hcpx                                                   (input)  
 * x_hcpx      - (x=1/2) approximate charges                (input)
 * Various parameters from egb_nbond_hcp()		    (input)

 ************************************************************************/ 

void egb_calc_deriv_approx(int j, REAL_T *q_hcp, REAL_T *r_hcp, int i,
                REAL_T *f,   
                REAL_T xi, REAL_T yi, REAL_T zi, REAL_T wi, 
                REAL_T *daix, REAL_T *daiy, REAL_T *daiz, REAL_T *daiw, 
                REAL_T sumda, REAL_T ri, REAL_T ri1i)
{

        int j34, q;
        REAL_T xij,yij,zij, r2, dij, dij1i, dij2i, sj, sj2, tmpsd, dumbo, datmp;
        REAL_T dedx, dedy, dedz, temp1, dij3i;


        for (q = 0; q < hcp; q++)
        {
                j34 = 4 * (j * hcp + q);

                xij = xi - q_hcp[j34];
                yij = yi - q_hcp[j34 + 1];
                zij = zi - q_hcp[j34 + 2];
                r2 = xij * xij + yij * yij + zij * zij;

                dij1i = 1.0 / sqrt(r2);
                dij2i = dij1i * dij1i;

                dij = r2 * dij1i;

                sj = r_hcp[j * hcp + q];
                sj2 = sj * sj;


                if (dij > 4.0 * sj) {

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
                                -0.25 * (-0.5 * (r2 - ri * ri + sj2) * dij3i
                                                * ri1i * ri1i + dij1i * temp1 * (temp1 - dij1i)
                                                - dij3i * log(ri * temp1));

                } else if (ri < sj) {

                        temp1 = 1. / (r2 - sj2);
                        datmp =
                                -0.5 * (sj * dij2i * temp1 - 2. * sj * temp1 * temp1 -
                                                0.5 * dij2i * dij1i * log((sj - dij) /
                                                        (sj + dij)));

                } else {
                        datmp = 0.;
                }


                datmp *= sumda * 2;   /* to account for fji due to component */
                dedx = xij * datmp;
                dedy = yij * datmp;
                dedz = zij * datmp;

                /* Sum the derivatives into daix, daiy, daiz and daiw. */
                *daix += dedx; 
                *daiy += dedy; 
                *daiz += dedz; 

        }   /* end for each hcp charge q */


}

/***********************************************************************
 *                         EGB_CALC_DERIV_ATOM()
 *
 * Calculate effective Born radii for residues
 *
 * Calling parameters are as follows:
 * npairs_hcpx - (x=0/1/2) number of pairs in each pairlist (input)
 * pairs_hcpx  - (x=0/1/2) list of atom/residue/chains      (input)
 * f                                                        (input)  
 * x           - (x=1/2) approximate charges                (input)
 * Various parameters from egb_nbond_hcp()		    (input)

 ************************************************************************/ 

void egb_calc_deriv_atom(int j, REAL_T *x, REAL_T *f, REAL_T ri, REAL_T ri1i, REAL_T sumda,
                REAL_T xi, REAL_T yi, REAL_T zi, REAL_T wi, 
                REAL_T *daix, REAL_T *daiy, REAL_T *daiz, REAL_T *daiw)
{

        int j34;
        REAL_T xij,yij,zij, r2, dij, dij1i, dij2i, dij3i, sj, sj2, tmpsd, temp1, dumbo, datmp;
        REAL_T dedx, dedy, dedz;


        j34 = 3 * j;

        xij = xi - x[j34];
        yij = yi - x[j34 + 1];
        zij = zi - x[j34 + 2];
        r2 = xij * xij + yij * yij + zij * zij;

        dij1i = 1.0 / sqrt(r2);
        dij2i = dij1i * dij1i;
        dij = r2 * dij1i;

        sj = prm->Fs[j] * (prm->Rborn[j] - gboffset);
        sj2 = sj * sj;

        /*
         * The following are the numerator of the first derivatives of the
         * effective radius Ri with respect to the interatomic distance Dij.
         * They are derived from the equations from the Appendix of Schaefer
         * and Froemmel as well as from the Taylor series expansion for d>>s
         * by Andreas Svrcek-Seiler.  The smooth rgbmax idea is from Andreas
         * Svrcek-Seiler and Alexey Onufriev.  The complete derivative is
         * formed by multiplying the numerator by -Ri*Ri.  The factor of Ri*Ri
         * has been moved to the terms that are multiplied by the derivative.
         * The negation is deferred until later.  When the chain rule is used
         * to form the first derivatives of the effective radius with respect
         * to the cartesian coordinates, an additional factor of Dij appears
         * in the denominator.  That factor is included in the following
         * expressions.
         */

        if (dij > 4.0 * sj) {

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
                        -0.25 * (-0.5 * (r2 - ri * ri + sj2) * dij3i
                                        * ri1i * ri1i + dij1i * temp1 * (temp1 - dij1i)
                                        - dij3i * log(ri * temp1));

        } else if (ri < sj) {

                temp1 = 1. / (r2 - sj2);
                datmp =
                        -0.5 * (sj * dij2i * temp1 - 2. * sj * temp1 * temp1 -
                                        0.5 * dij2i * dij1i * log((sj - dij) /
                                                (sj + dij)));

        } else {
                datmp = 0.;
        }

        datmp *= sumda;
        dedx = xij * datmp;
        dedy = yij * datmp;
        dedz = zij * datmp;

        /* Sum the derivatives into daix, daiy, daiz and daiw. */

        *daix += dedx;
        *daiy += dedy;
        *daiz += dedz;

        /*
         * Sum the derivatives relative to atom j (weighted by -sumdeijda[i])
         * into the gradient vector.  For example, f[j34 + 2] contains the
         * derivatives of Ri with respect to the z-coordinate of atom j.
         */

        f[j34] += dedx;
        f[j34 + 1] += dedy;
        f[j34 + 2] += dedz;


}

/***********************************************************************
 *                          egb_hcp_nbond()
 *
 * Calculate the generalized Born energy and first derivatives.
 *
 * Calling parameters are as follows:
 *
 * npairs_hcpx   - (x=0/1/2) number of pairs in pairlist  (input)
 * pairlist_hcpx - (x=0/1/2) pairlist for atoms/res/chain (input)
 * Iblo_hcp      - number of excluded atoms               (input)
 * IexclAt       - excluded atoms list                    (input)
 * x             - atom cordinates                        (input)
 * x_hcpx        - (x=1/2) residue/chain geom centers     (input)
 * q_hcpx        - (x=1/2) residue/chain hcp charges      (input)
 * r_hcpx        - (x=1/2) residue/chain radii            (input)
 * prm           - parameter/topology structure           (input)
 * reff_hcpx     - (x=0/1/2) effective born radii         (input)
 * psi           - intermediate computational values      (input)
 * f             - gradient vector                        (update)
 * kappa         - GB parameter                           (input)
 * diel_ext      - external dielectric                    (input)
 * enb           - VDW energy                             (update)
 * eelt          - Coulomb energy                         (update)
 * esurf         - GBSA energy                            (update)
 * enp           -
 *
 * Return value: GB energy
 *
 ************************************************************************/
/*
 * Compute the non-polar contributions, depending on the value of gbsa.
 */
REAL_T egb_hcp_nbond(
                INT_T * Iblo_hcp, INT_T ** IexclAt_hcp, REAL_T * x, REAL_T * x_hcp1, REAL_T * x_hcp2, REAL_T * x_hcp3, 
                REAL_T * q_hcp1, REAL_T * q_hcp2, REAL_T * q_hcp3, REAL_T * r_hcp1, REAL_T * r_hcp2, REAL_T * r_hcp3,
                PARMSTRUCT_T * prm, REAL_T * reff_hcp0, REAL_T * reff_hcp1, REAL_T * reff_hcp2, REAL_T * reff_hcp3,
                REAL_T * psi, REAL_T * f, REAL_T * kappa, REAL_T * diel_ext, 
                REAL_T * enb, REAL_T * eelt, REAL_T * esurf, REAL_T * enp)
{

        REAL_T *sumdeijda = NULL;
        int *iexw = NULL;

        int i, i34, j;
        int iaci;
        REAL_T epol, dielfac, qi, expmkf;
        REAL_T elec, evdw, sumda, daix, daiy, daiz, daiw;
        REAL_T xi, yi, zi, wi; 
        REAL_T qi2h, qid2h; 
        REAL_T ri1i;

        REAL_T ri, thi;

        REAL_T evdwnp;

        /* hcp structure index */
        int a, r, s, c, a_from, a_to, r_from, r_to, s_from, s_to;
        REAL_T dist2, dist2_hcp1, dist2_hcp2, dist2_hcp3;

        /* square of threshold distances */
        dist2_hcp1 = hcp_h1 * hcp_h1;
        dist2_hcp2 = hcp_h2 * hcp_h2;
        dist2_hcp3 = hcp_h3 * hcp_h3;


        if (sumdeijda == NULL) {
                sumdeijda = vector(0, prm->Natom);
        }


        /* Compute the GB, Coulomb and Lennard-Jones energies and derivatives. */

        epol = elec = evdw = evdwnp = 0.0;

        /*
         * Initialize the sumdeijda array inside of the parallel region.
         *
         * For MPI and ScaLAPACK, each process has its own copy of the
         * array which must be initialized in its entirety because a
         * call to MPI_Allreduce will be used to reduce the array.
         * The MPI reduction will synchronize the processes.
         */

        //      for (i = 0; i < prm->Natom; i++) {
        //	sumdeijda[i] = 0.0;
        //      }

        /* iexw stores the list of excluded atoms of each atom */
        iexw = ivector(-1, prm->Natom+1);
        for (i = -1; i < prm->Natom+1; i++) {
                iexw[i] = -1;
        }


        /*Step 1 - Calculate energy */

        for (i = 0; i < prm->Natom; i++) {

#if defined(MPI) || defined(SCALAPACK)
                if (!myroc(i, blocksize, numtasks, mytaskid))
                        continue;
#endif

                sumdeijda[i] = 0.0;

                ri = reff_hcp0[i];
                qi = prm->Charges[i];

                /*
                 * If atom i is not frozen, compute the "diagonal" energy that
                 * is a function of only the effective radius Ri but not of the
                 * interatomic distance Dij.  Compute also the contribution of
                 * the diagonal energy term to the sum by which the derivative
                 * of Ri will be multiplied.  If requested, calculate the van
                 * der Waals component of the non-polar solvation free energy
                 * and its derivatives.
                 */

                if (!frozen[i]) {
                        expmkf = exp(-KSCALE * (*kappa) * ri) / (*diel_ext);
                        dielfac = 1.0 - expmkf;
                        qi2h = 0.5 * qi * qi;
                        qid2h = qi2h * dielfac;
                        epol += -qid2h * 2.0 / ri; /* HCP: multiply by 2.0 */

                        sumdeijda[i] +=
                             qid2h - KSCALE * (*kappa) * qi2h * expmkf * ri;
                }

                i34 = dim * i;

                xi = x[i34];
                yi = x[i34 + 1];
                zi = x[i34 + 2];

                iaci = prm->Ntypes * (prm->Iac[i] - 1);

                /*
                 * Expand the excluded atom list into the iexw array by storing i
                 * at array address j.
                 */

                for (j = 0; j < Iblo_hcp[i]; j++) {
                        iexw[IexclAt_hcp[i][j] - 1] = i;
                }

                /* Initialize the derivative accumulators. */

                daix = daiy = daiz = daiw = 0.0;

                for (c=0; c < prm->Ncomplex; c++)   /* for each complex in the molecule */
                {
                        dist2 = calc_dist2(xi, yi, zi,
                                        x_hcp3[3*c+0], x_hcp3[3*c+1], x_hcp3[3*c+2]);
                        if (dist2 > dist2_hcp3)
                        {
                                egb_calc_energy_approx(c, q_hcp3, reff_hcp3, hcp, i, qi, ri, kappa, diel_ext,
                                                &elec, &evdw, &epol, sumdeijda, xi, yi, zi,
                                                &daix, &daiy, &daiz);
                        }
                        else
                        {
                                s_from = prm->Ipcomplex[c] - 1;
                                if (c + 1 < prm->Ncomplex) { s_to = prm->Ipcomplex[c + 1] - 1; }
                                else { s_to = prm->Nstrand; }

                                for (s = s_from; s < s_to; s++)   /* for each strand in the complex */
                                {
                                        dist2 = calc_dist2(xi, yi, zi,
                                                        x_hcp2[3*s+0], x_hcp2[3*s+1], x_hcp2[3*s+2]);
                                        if (dist2 > dist2_hcp2)
                                        {
                                                egb_calc_energy_approx(s, q_hcp2, reff_hcp2, hcp, i, qi, ri, kappa, diel_ext,
                                                                &elec, &evdw, &epol, sumdeijda, xi, yi, zi,
                                                                &daix, &daiy, &daiz);
                                        }
                                        else
                                        {
                                                r_from = prm->Ipstrand[s] - 1;
                                                if (s + 1 < prm->Nstrand) { r_to = prm->Ipstrand[s + 1] - 1; }
                                                else { r_to = prm->Nres; }

                                                for (r = r_from; r < r_to; r++)         /* for each residue in the strand */
                                                {

                                                        dist2 = calc_dist2(xi, yi, zi,
                                                                        x_hcp1[3*r+0], x_hcp1[3*r+1], x_hcp1[3*r+2]);
                                                        if (dist2 > dist2_hcp1)
                                                        {
                                                                egb_calc_energy_approx(r, q_hcp1, reff_hcp1, hcp, i, qi, ri, kappa, diel_ext,
                                                                                &elec, &evdw, &epol, sumdeijda, xi, yi, zi,
                                                                                &daix, &daiy, &daiz);
                                                        }
                                                        else
                                                        {
                                                                a_from = prm->Ipres[r];
                                                                if (r + 1 < prm->Nres) { a_to = prm->Ipres[r + 1]; }
                                                                else { a_to = prm->Natom; }

                                                                for (a = a_from; a < a_to; a++)     /* for each atom in the atom */
                                                                {
                                                                        if (a != i)   /* not same atom */
                                                                        {
                                                                                egb_calc_energy_atom(a, x, reff_hcp0, iexw, hcp, i, iaci, qi, ri, f, kappa, diel_ext,
                                                                                                &elec, &epol, &evdw, sumdeijda, xi, yi, zi, wi, &daix, &daiy, &daiz, &daiw);
                                                                        }   /* end if not same atom */
                                                                }   /* end for each atom a */
                                                        }   /* end else in residue threshold dist */
                                                }   /* end for each residue r */
                                        }   /* end else in strand threshod dist */
                                }   /* end for each strand s */
                        }   /* end else in complex threshold dist */
                }   /* end for each complex c */


                /* Update the i elements of the gradient array */
                f[i34] += daix;
                f[i34 + 1] += daiy;
                f[i34 + 2] += daiz;

        }  /* end for each atom i */


        /*
         * Compute the derivatives of the effective radius Ri of atom i
         * with respect to the cartesian coordinates of each atom j.  Sum
         * all of these derivatives into the gradient vector.
         *
         * Loop over all atoms i.
         *
         * Synchronization of OpenMP threads will occur following this
         * loop nest because of the '#pragma omp for'.
         * 
         * A reduction of the gradient array will occur in the mme34 function,
         * either for OpenMP or MPI.  This reduction will synchronize the MPI
         * tasks, so an explicit barrier is not necessary following this loop.
         *
         * For MPI or ScaLAPACK, explicitly assign tasks to loop indices
         * for the following loop in a manner equivalent to (static, N)
         * scheduling for OpenMP.  For OpenMP use (dynamic, N) scheduling.
         */


        /* Step 2 - calculate derivative of solvation energy */

        for (i = 0; i < prm->Natom; i++) {

#if defined(MPI) 
                if (!myroc(i, blocksize, numtasks, mytaskid))
                        continue;
#endif

                /*
                 * Don't calculate derivatives of the effective radius of atom i
                 * if atom i is frozen or if there are no pair atoms j associated
                 * with atom i.
                 */

                if ( frozen[i] )
                        continue;

                i34 = dim * i;

                xi = x[i34];
                yi = x[i34 + 1];
                zi = x[i34 + 2];

                ri = prm->Rborn[i] - gboffset;
                ri1i = 1. / ri;

                sumda = sumdeijda[i];

                if (gb > 1) {

                        ri = prm->Rborn[i] - gboffset;
                        thi =
                                tanh((gbalpha[i] - gbbeta[i] * psi[i] +
                                                        gbgamma[i] * psi[i] * psi[i]) * psi[i]);
                        sumda *=
                                (gbalpha[i] - 2.0 * gbbeta[i] * psi[i] +
                                 3.0 * gbgamma[i] * psi[i] * psi[i])
                                * (1.0 - thi * thi) * ri / prm->Rborn[i];
                }
                /* Initialize the derivative accumulators. */

                daix = daiy = daiz = daiw = 0.0;

                for (c=0; c < prm->Ncomplex; c++)   /* for each complex in the molecule */
                {
                        dist2 = calc_dist2(xi, yi, zi,
                                        x_hcp3[3*c+0], x_hcp3[3*c+1], x_hcp3[3*c+2]);
                        if (dist2 > dist2_hcp3)
                        {
                                egb_calc_deriv_approx(c, q_hcp3, r_hcp3, i, f, 
                                                xi, yi, zi, wi, &daix, &daiy, &daiz, &daiw, sumda, ri, ri1i);
                        }
                        else
                        {
                                s_from = prm->Ipcomplex[c] - 1;
                                if (c + 1 < prm->Ncomplex) { s_to = prm->Ipcomplex[c + 1] - 1; }
                                else { s_to = prm->Nstrand; }

                                for (s = s_from; s < s_to; s++)   /* for each strand in the complex */
                                {
                                        dist2 = calc_dist2(xi, yi, zi,
                                                        x_hcp2[3*s+0], x_hcp2[3*s+1], x_hcp2[3*s+2]);
                                        if (dist2 > dist2_hcp2)
                                        {
                                                egb_calc_deriv_approx(s, q_hcp2, r_hcp2, i, f, 
                                                                xi, yi, zi, wi, &daix, &daiy, &daiz, &daiw, sumda, ri, ri1i);
                                        }
                                        else
                                        {
                                                r_from = prm->Ipstrand[s] - 1;
                                                if (s + 1 < prm->Nstrand) { r_to = prm->Ipstrand[s + 1] - 1; }
                                                else { r_to = prm->Nres; }

                                                for (r = r_from; r < r_to; r++)         /* for each residue in the strand */
                                                {
                                                        dist2 = calc_dist2(xi, yi, zi,
                                                                        x_hcp1[3*r+0], x_hcp1[3*r+1], x_hcp1[3*r+2]);
                                                        if (dist2 > dist2_hcp1)
                                                        {
                                                                egb_calc_deriv_approx(r, q_hcp1, r_hcp1, i, f, 
                                                                                xi, yi, zi, wi, &daix, &daiy, &daiz, &daiw, sumda, ri, ri1i);
                                                        }
                                                        else
                                                        {
                                                                a_from = prm->Ipres[r];
                                                                if (r + 1 < prm->Nres) { a_to = prm->Ipres[r + 1]; }
                                                                else { a_to = prm->Natom; }

                                                                for (a = a_from; a < a_to; a++)     /* for each atom in the residue */
                                                                {
                                                                        if (a != i)   /* not same atom */
                                                                        {
                                                                                egb_calc_deriv_atom(a, x, f, ri, ri1i, sumda,
                                                                                                xi, yi, zi, wi, &daix, &daiy, &daiz, &daiw);
                                                                        }   /* end if not same atom */
                                                                }   /* end for each atom a */
                                                        }   /* end else in residue threshold dist */
                                                }   /* end for each residue r */
                                        }   /* end else in strand threshod dist */
                                }   /* end for each strand s */
                        }   /* end else in complex threshold dist */
                }   /* end for each complex c */


                /*
                 * Update the gradient vector with the sums of derivatives of the
                 * effective radius Ri with respect to the cartesian coordinates.
                 * For example, f[i34 + 1] contains the sum of derivatives of Ri
                 * with respect to the y-coordinate of each atom.  Multiply by
                 * -sumdeijda[i] here (instead of merely using datmp multiplied by
                 * -sumdeijda) in order to distribute the product across the sum of
                 * derivatives in an attempt to obtain greater numeric stability.
                 */

                f[i34] -= daix;
                f[i34 + 1] -= daiy;
                f[i34 + 2] -= daiz;


        }   /* end for each atom i */

        /*
         * Return elec, evdw and evdwnp through the parameters eelt, enb and enp.
         * These variables are computed in parallel.
         */

        *eelt = elec * 0.5;
        *enb = evdw * 0.5;
        *enp = evdwnp * 0.5;

        return (epol * 0.5);
}


/***********************************************************************
 *                           egb_hcp() 
 *
 * Calculate the non-bonded energy and first derivatives using GB with HCP.
 * The non-bonded pair list must be modified by the excluded atom list 
 * 
 * Calling parameters are as follows:
 *
 * npairs_hcpx   - (x=0/1/2) number of pairs in pairlist  (input)
 * pairlist_hcpx - (x=0/1/2) pairlist for atoms/res/chain (input)
 * Iblo_hcp      - number of excluded atoms               (input)
 * IexclAt       - excluded atoms list                    (input)
 * x             - atom cordinates                        (input)
 * x_hcpx        - (x=1/2) residue/chain geom centers     (input)
 * q_hcpx        - (x=1/2) residue/chain hcp charges      (input)
 * r_hcpx        - (x=1/2) residue/chain radii            (input)
 * prm           - parameter/topology structure           (input)
 * f             - gradient vector                        (update)
 * kappa         - GB parameter                           (input)
 * epsext        - external dielectric                    (input)
 * enb           - VDW energy                             (update)
 * eel           - Coulomb energy                         (update)
 * esurf         - GBSA energy                            (update)
 * evdwnp        -
 *
 * Return value: GB energy
 *
 *************************************************************************/

REAL_T egb_hcp(int * Iblo_hcp, int ** IexclAt_hcp, 
                REAL_T * x, REAL_T * x_hcp1, REAL_T * x_hcp2, REAL_T * x_hcp3,
                REAL_T * q_hcp1, REAL_T * q_hcp2, REAL_T * q_hcp3,
                REAL_T * r_hcp1, REAL_T * r_hcp2, REAL_T * r_hcp3, 
                REAL_T * reff_hcp0, REAL_T * reff_hcp1, REAL_T * reff_hcp2, REAL_T * reff_hcp3, 
                PARMSTRUCT_T * prm, REAL_T * f, REAL_T * kappa, REAL_T * epsext, 
                REAL_T * enb, REAL_T * eel, REAL_T * esurf, REAL_T * evdwnp)
{
        REAL_T *psi = NULL;
        REAL_T e_gb;

        /* calculate approximate charges and geometric centers */
        if (hcpo==0) {
          calc_approx_q(x, x_hcp1, x_hcp2, x_hcp3, q_hcp1, q_hcp2, q_hcp3, hcp, prm);
        } else {
            calc_approx_q_hcpo(x, x_hcp1, x_hcp2, x_hcp3, q_hcp1, q_hcp2, q_hcp3, hcp, prm);
        }


        /* calculate effective born radii for atoms */
        reff_hcp0 = vector(0, prm->Natom); 
        psi = vector(0, prm->Natom); 

        egb_hcp_atom_reff(x, x_hcp1, x_hcp2, x_hcp3,
                        q_hcp1, q_hcp2, q_hcp3,
                        prm->Rborn, prm->Fs, r_hcp1, r_hcp2, r_hcp3, reff_hcp0, psi);

        /* calcluate effective born radii for higher level groupings */
        reff_hcp1 = vector(0, prm->Nres * hcp);
        reff_hcp2 = vector(0, prm->Nstrand * hcp);
        reff_hcp3 = vector(0, prm->Ncomplex * hcp);
       // if (HCP_RAD > 0)
       // {
       //    printf("hcpo= %d  , hcp=   %d, dhcp=   %f\n", hcpo, hcp, dhcp);
            if (hcpo==0) 
            {
                egb_hcp_reff(prm, reff_hcp0, reff_hcp1, reff_hcp2, reff_hcp3);  /*averaging based on AO methodss*/
            }
            else
            {
                 egb_hcpo_reff_avg(prm, reff_hcp0, reff_hcp1, reff_hcp2, reff_hcp3, q_hcp1, q_hcp2, q_hcp3,
                                x_hcp1, x_hcp2, x_hcp3, x); /* new averaging for hcpo, Izadi et al JCTC 2016 */
            } 
        //}
        //else
       // {
        //        egb_hcp_reff_new(x_hcp1, q_hcp1, r_hcp1, reff_hcp1, 
       //                         x_hcp2, q_hcp2, r_hcp2, reff_hcp2, 
        //                        x_hcp3, q_hcp2, r_hcp3, reff_hcp3, 
        //                        hcp_h1, hcp_h2, hcp_h3, prm);
       // }

        /* calculate nonbonded energy and force using GB */
        e_gb = egb_hcp_nbond(
           Iblo_hcp, IexclAt_hcp, x, x_hcp1, x_hcp2, x_hcp3,
           q_hcp1, q_hcp2, q_hcp3, r_hcp1, r_hcp2, r_hcp3,
           prm, reff_hcp0, reff_hcp1, reff_hcp2, reff_hcp3,
           psi, f, kappa, epsext, enb, eel, esurf, evdwnp);

        /* calculate surface energy and force using GBSA */
        if (gbsa)
        {
                hcp_gbsa(x, f, prm, esurf);
        }

        free_vector(reff_hcp0, 0, prm->Natom);
        free_vector(reff_hcp1, 0, prm->Nres * hcp);
        free_vector(reff_hcp2, 0, prm->Nstrand * hcp);
        free_vector(reff_hcp3, 0, prm->Ncomplex * hcp);
        free_vector(psi, 0, prm->Natom);

        return(e_gb);
}


