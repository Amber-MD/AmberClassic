#define PIx4 (12.5663706143591724639918)
#define PIx2 ( 6.2831853071795862319959)
#define PIx1 ( 3.1415926535897931159979)
   /*
    * Calculate either the volume or the surface area component of the AGB
    * nonpolar hydration free energy.  Use Eqs 2-15 of Gallicchio and Levy,
    * J. Comput. Chem. 25: 479 (2004).
    *
    * The gbsa option is interpreted as follows:
    * 0 - no nonpolar calculation
    * 1 - original LCPO approximation
    *
    * Calculate the cavity (volume or surface area) component of the nonpolar
    * GB energy.  Note that for OpenMP execution, the totalvolume and surfacearea
    * variables are reduction variables.  For MPI, the reduction occurs in the
    * mme34 function via reduction of the ene array.
    */

   if (gbsa == 1) {

      /*
       * Main LCPO stuff follows.
       */

      /* Loop over all atoms i. */

      count = 0;
      for (i = 0; i < prm->Natom; i++) {

         xi = x[3 * i];
         yi = x[3 * i + 1];
         zi = x[3 * i + 2];
         ri = rborn[i] - gboffset;
         ri1i = 1. / ri;
         sumi = 0.;

         /* Select atom j from the pair list. */

         for (k = 0; k < lpears[i] + upears[i]; k++) {

            j = pearlist[i][k];

            xij = xi - x[3 * j];
            yij = yi - x[3 * j + 1];
            zij = zi - x[3 * j + 2];
            r2 = xij * xij + yij * yij + zij * zij;
            if (r2 > rgbmaxpsmax2)
               continue;
            dij1i = 1.0 / sqrt(r2);
            dij = r2 * dij1i;

            if ((P0[i] + P0[j]) > dij) {
               if (P0[i] > 2.5 && P0[j] > 2.5) {

                  ineighbor[count] = j + 1;
                  /* this +1 is VERY important */
                  count += 1;
               }
            }
         }

         ineighbor[count] = 0;
         count = count + 1;
      }

      totsasa = 0.0;
      count = 0;
      for (i = 0; i < prm->Natom; i++) {
         if (ineighbor[count] == 0) {
            count = count + 1;
         } else {

            si = PIx4 * P0[i] * P0[i];
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
            Aij = PIx2 * P0[i] * tmpaij;

            dAijddij =
                PIx1 * P0[i] * (dij1i * dij1i *
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

               Ajk = PIx1 * P0[j] * tmpajk;

               sumAjk = sumAjk + Ajk;
               sumAjk2 = sumAjk2 + Ajk;

               dAjkddjk =
                   PIx1 * P0[j] * djk1i * (djk1i * djk1i * vdw2dif - 1.0);

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
            /*  printf("totsasa up to %d %f\n",i,totsasa); */
         }
      }
      /*      printf("SASA %f , ESURF %f \n",totsasa,totsasa*surften); */
      *esurf = totsasa * surften;

   } else {
      *esurf = 0.0;
   }
