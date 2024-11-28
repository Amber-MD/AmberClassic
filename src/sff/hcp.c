/**********************************************************************
 * implementation of Hierarchical Charge Partitioning (HCP):          *
 * Izadi, Aandakrishnan and Onufriev, JCTC 2016 (in prep)             *
 * Aandakrishnan and Onufriev, Journal of Computational Chemistry     *
 * July 2009                                                          *
 **********************************************************************/

#define SIG 0.3
#define DIW 78.0
#define C1 38.5

#define MIN_Q   0.00001  /* assume 0 if less than this */


/**********************************************************************
 * 
 *    Calculate eigenvalues and eigenvectors for OPCA
 *
 * ********************************************************************/
void eigen_decomposition(double mQ[3][3], double eigvec[3][3], double eigval[3]) {
  int i,imax,imin,imean;
  double A[3],B[3],C[3],eigvals[3];
  double k0,k1,k2,p,q,s,alpha,pi,arc,arc1,AA,BB,CC,eigvecmag;

   pi=3.14159265;  
  // Analytical equations to find eigen values and eigen vetors for a 3X3 matrix (i.e. the quadrupole)
  k0=mQ[0][0]*mQ[1][1]*mQ[2][2]+2.0*mQ[0][1]*mQ[1][2]*mQ[0][2]- 
     mQ[0][0]*mQ[1][2]*mQ[1][2]-mQ[1][1]*mQ[0][2]*mQ[0][2]-
     mQ[2][2]*mQ[0][1]*mQ[0][1];
  k1=-mQ[0][0]*mQ[1][1]-mQ[1][1]*mQ[2][2]-mQ[0][0]*mQ[2][2] 
     +mQ[0][1]*mQ[0][1]+mQ[1][2]*mQ[1][2]+mQ[0][2]*mQ[0][2];
  k2=mQ[0][0]+mQ[1][1]+mQ[2][2];

  p=-k1-k2*k2/3.0;
  q=-k0-2.0*k2*k2*k2/27.0-k1*k2/3.0;

  if ((-p)<0){printf("ERROR: eigen_decomposition, sqrt(negative value)");}
  s=-2.0*sqrt(-p/3.0);
  arc=-4.0*q/(s*s*s);

  alpha=acos(arc);
  eigvals[0]= s*cos(alpha/3.0)       +k2/3.0;
  eigvals[1]=-s*cos((alpha+pi)/3.0)  +k2/3.0;
  eigvals[2]=-s*cos((alpha-pi)/3.0)  +k2/3.0;

  eigval[0]=eigvals[0];
  eigval[1]=eigvals[1];
  eigval[2]=eigvals[2];

//if (fabs(fabs(eigvals[2])- fabs(eigvals[0])) > 0.5*fabs(eigvals[0])) {  
//  if (fabs(eigvals[2]) > fabs(eigvals[0])) {
//    eigval[0]=eigvals[2];
//    eigval[2]=eigvals[0];
//   }
//}

  for ( i = 0; i < 3; i++) { 

    A[i]=mQ[0][0]-eigval[i];
    B[i]=mQ[1][1]-eigval[i];
    C[i]=mQ[2][2]-eigval[i];
    AA=(mQ[0][2]*mQ[0][1]-A[i]*mQ[1][2]);
    BB=(mQ[0][1]*mQ[1][2]-B[i]*mQ[0][2]);
    CC=(mQ[0][2]*mQ[1][2]-C[i]*mQ[0][1]);
    eigvec[0][i]=BB*CC;
    eigvec[1][i]=CC*AA;
    eigvec[2][i]=AA*BB;
    eigvecmag=sqrt(eigvec[0][i]*eigvec[0][i]+eigvec[1][i]*eigvec[1][i]
             +eigvec[2][i]*eigvec[2][i]);
    if (eigvecmag != 0) {
      eigvec[0][i]=eigvec[0][i]/eigvecmag;
      eigvec[1][i]=eigvec[1][i]/eigvecmag;
      eigvec[2][i]=eigvec[2][i]/eigvecmag;
    }
  }

}



/***********************************************************************
 *                          calc_dist2()
 *
 * Calculate distance square
 *
 * Calling Parameters:
 * xi, yi, zi - point i
 * xj, yj, zj - point j
 * 
 * return: dist2
 ************************************************************************/

static REAL_T calc_dist2(REAL_T xi, REAL_T yi, REAL_T zi, REAL_T xj, REAL_T yj, REAL_T zj)
{
	REAL_T dist2, xij, yij, zij;

	xij = xi - xj;
	yij = yi - yj;
	zij = zi - zj;
	dist2 = xij * xij + yij * yij + zij * zij;

	return(dist2);
}

/***********************************************************************
 *                          q2idx()
 *
 * Determine index for grouping charges for hcp approximation
 *
 * Calling Parameters: q (charge)
 * 
 * return: idx
 ************************************************************************/

int q2idx(REAL_T q)
{
	int idx = 0;

	switch (hcp)
	{
		/* cutoff - should never get here */
		case 0:
			break;
			/* 1 charge (idx=0) */
		case 1:
			break;
			/* 2 charges (pos:idx=1 & neg:idx=0) */
		case 2:
			if (q > 0.0) { idx = 1; }
			break;
			/* 3 charges (idx=0/1/2) */
		case 3:
			if (q > -0.333)
			{
				idx = 1;
				if (q > 0.333) { idx = 2; }
			}
			break;
	} 

	return(idx);
}




/************************************************************************
*  
* Izadi, Anandakrishnan and Onufriev, JCTC 2016 (in prep)
* 
*   FUNCTION: OPCA_K2 -- CALCULATES OPCA K=2 FOR CHARGED & NUETRAL CASES
*INPUTS:   MULTIPOLE MOMENTS: MONOPOLE, DIPOLE, QUADRUPOLE
*OUTPUT:   updated with OPCA point charges
*
*************************************************************************/
 void OPCA_K2( int res, float qtotal, float p[3], float Q[3][3], float rmax, float rcog[3],
                float* qbar1, float rbar1[3], float* qbar2, float rbar2[3], float R1, float R2, float coef2) {

   int     a, i, j;                    /* loop counters */
   float   rcoq[3],rcoq2,prcoq;        /* to transfer the quadrupole */
   double A[3][3], V[3][3], d[3];      /* FOR the eigen value calculations*/
   float cent_of_dip[3], cdotq;
   float dis = 0.0;
   float R1s, R2s, dis2, denom, dis_mod;
   float rbar, coef=1.5;  
   float cent_of_dip_mod[3];
   float sumqpp=0.0, mpsq, dk[3], sumqp;
   float dr, q;
   float rbar1abs=0.0, rbar2abs=0.0;


   *qbar1=0.0;
   *qbar2=0.0;
   rbar1[0] = 0.0;
   rbar1[1] = 0.0;
   rbar1[2] = 0.0;
   rbar2[0] = 0.0;
   rbar2[1] = 0.0;
   rbar2[2] = 0.0;
 
   if (fabs(qtotal) > 0.6*18.2223) {   // This part calculates the 2-charge OPCA for charged distributions


      /* calculate the center of charge */
      rcoq[0]=p[0]/qtotal;
      rcoq[1]=p[1]/qtotal;
      rcoq[2]=p[2]/qtotal;

      //  Shifting the center of expansion to the center of charge for quadrupole
      rcoq2= rcoq[0] * rcoq[0] +rcoq[1] * rcoq[1] +rcoq[2] * rcoq[2];
      prcoq= p[0] * rcoq[0]+ p[1] * rcoq[1]+ p[2] * rcoq[2];
      for (i=0; i<3; i++) {
        for (j=i; j<3; j++) {
                Q[i][j] = 0.5 * ( 2.0 * Q[i][j] - (3.0 * p[i] * rcoq[j] + 3.0 * rcoq[i] *  p[j] )
                        +   qtotal * (3.0 * rcoq[i] * rcoq[j]));
        }
      }
      Q[1][0] =  Q[0][1];
      Q[2][0] =  Q[0][2];
      Q[2][1] =  Q[1][2];

      for (i=0; i<3; i++) {
                Q[i][i] = Q[i][i] +  prcoq - 0.5 * qtotal *rcoq2;
      }
      for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
                A[i][j]=Q[i][j];
        }
      }

      // Find principal directions and values
      eigen_decomposition( A, V, d);  /* eigen value calculations */


      rbar = coef * rmax;
      if (((rbar*rbar*qtotal)+d[0]) == 0.0) {rbar=0.0;}  /* to avoid divide by zero, 
                                                          switch to 1-q OPCA if denominator is zero. */

      (*qbar1)=qtotal * d[0] / ((rbar * rbar * qtotal) + d[0]);
      (*qbar2)=qtotal-(*qbar1);
      rbar1[0] = rbar * V[0][0];
      rbar1[1] = rbar * V[1][0];
      rbar1[2] = rbar * V[2][0];

      if ((*qbar2) != 0.0) {
         rbar2[0]=-(*qbar1)/ (*qbar2)*rbar1[0];
         rbar2[1]=-(*qbar1)/ (*qbar2)*rbar1[1];
         rbar2[2]=-(*qbar1)/ (*qbar2)*rbar1[2];
      }
      else {
         rbar2[0]=0.0;
         rbar2[1]=0.0;
         rbar2[2]=0.0;
      }

     /* Selecting the dipole that its center is closer to the center of geometry*/
      cent_of_dip[0] = (rbar1[0]+rbar2[0])/2.0;
      cent_of_dip[1] = (rbar1[1]+rbar2[1])/2.0;
      cent_of_dip[2] = (rbar1[2]+rbar2[2])/2.0;

      cdotq = cent_of_dip[0] * (rcoq[0]-rcog[0])
            + cent_of_dip[1] * (rcoq[1]-rcog[1])
            + cent_of_dip[2] * (rcoq[2]-rcog[2]);

      for (i=0; i<3; i++) {
        if (cdotq > 0) {
          rbar1[i]=-rbar1[i];
          rbar2[i]=-rbar2[i];
        }
        rbar1[i]=rbar1[i]+rcoq[i]-rcog[i];
        rbar2[i]=rbar2[i]+rcoq[i]-rcog[i];
      }

   } /* end for (qtotal > 0.6) */
   else {   /* calculating 2-charge OPCA for neutral distributions */

     // calculating dk (the center of dipole) for OPCA neutral 
      for ( i=0; i<3; i++) {
        for (j=0; j<3; j++) {
          sumqpp = sumqpp + Q[i][j]*p[i]*p[j];
        }
      }
      mpsq=p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
      memset(dk, 0, sizeof(dk));
      for (i=0; i<3; i++) {
        sumqp=0.0;
        for (j=0; j<3; j++) {
            sumqp=sumqp+Q[i][j]*p[j];
         }
         dk[i]=dk[i]+2.0/(3.0*mpsq)*(sumqp-(sumqpp/(4.0*mpsq)*p[i]));
      }

      dr = coef2 * rmax;
      q=sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2)) / dr;
      *qbar1= q;
      rbar1[0]= +p[0]/(2.0*q) + dk[0] - rcog[0];
      rbar1[1]= +p[1]/(2.0*q) + dk[1] - rcog[1];
      rbar1[2]= +p[2]/(2.0*q) + dk[2] - rcog[2];
      *qbar2= -q;
      rbar2[0]= -p[0]/(2.0*q) + dk[0] - rcog[0];
      rbar2[1]= -p[1]/(2.0*q) + dk[1] - rcog[1];
      rbar2[2]= -p[2]/(2.0*q) + dk[2] - rcog[2];

   }   /* end for else for (q<=0.6) */

   // constraining the charges within a certain threshold (<R1<R2)
   rbar1abs=sqrt(rbar1[0]*rbar1[0]+rbar1[1]*rbar1[1]+rbar1[2]*rbar1[2]);
   rbar2abs=sqrt(rbar2[0]*rbar2[0]+rbar2[1]*rbar2[1]+rbar2[2]*rbar2[2]);  

   R1s=R1*R1;
   R2s=R2*R2;
   denom=(R1s-R2s)*(R1s-R2s)*(R1s-R2s);
   dis=rbar1abs;
   dis2=dis*dis;
   dis_mod =0.0;

   if (dis < R1) {
     dis_mod = dis;
   } else if (dis > R2) {
     dis_mod = R2;
   } else {
     dis_mod =(R1s-dis2)*(R1s-dis2)*(R1s+2.0*dis*dis-3.0*R2s)/denom*(R2-dis)+dis;
   }
   for (i=0; i<3; i++) {
      rbar1[i] = dis_mod/dis*rbar1[i]+rcog[i];
   }
   dis=rbar2abs;
   dis2=dis*dis;
   dis_mod =0.0;
   if (dis < R1) {
      dis_mod = dis;
   } else if (dis > R2) {
      dis_mod = R2;
   } else {
      dis_mod =(R1s-dis2)*(R1s-dis2)*(R1s+2.0*dis*dis-3.0*R2s)/denom*(R2-dis)+dis;
   }
   for (i=0; i<3; i++) {
      rbar2[i] = dis_mod/dis*rbar2[i]+rcog[i];
   } 
    
}  /* end for OPCA_k2 function*/


/***********************************************************************
 *                          calc_approx_q()
 *
 * Calculate approximate charges and geometric center of components
 *
 * Calling Parameters:
 * x      - atom coordinates                 (input)
 * x_hcp1 - geometric centers of residues    (update)
 * x_hcp2 - geometric centers of strands     (update)
 * q_hcp1 - approximate charges for residues (update)
 * q_hcp2 - approximate charges for strands  (update)
 * prm    - parmeter/topology structure      (input)
 ************************************************************************/
void calc_approx_q(REAL_T * x, REAL_T * x_hcp1, REAL_T * x_hcp2, REAL_T * x_hcp3,
		REAL_T * q_hcp1, REAL_T * q_hcp2, REAL_T * q_hcp3, INT_T hcp, PARMSTRUCT_T * prm)
{

	int c, s, r, a, m, j3, j4, j3a, j4a;           		/* index for strand, residue, atom, charge */
	int a_from, a_to;                  			/* for atoms within a residue */
	int r_from, r_to;		      			/* for residues within a strand */
	int s_from, s_to;		      			/* for strands within a complex */
	int natoms_hcp1, natoms_hcp2, natoms_hcp3;      	/* total number of atoms */
	REAL_T q;                          			/* charge */

	for (c = 0; c < prm->Ncomplex; c++) /* for each complex */
	{
		/* initialize complexes */
		natoms_hcp3 = 0;
		j3 = 3 * c;
		x_hcp3[j3+0] = x_hcp3[j3+1] = x_hcp3[j3+2] = 0.0;

		for (m = 0; m < hcp; m++) 
		{ 
			j4 = 4 * (c * hcp + m);
			q_hcp3[j4] = q_hcp3[j4+1] = q_hcp3[j4+2] = q_hcp3[j4+3] = 0.0;
		}
		s_from = prm->Ipcomplex[c] - 1;
		if (c < prm->Ncomplex - 1) { s_to = prm->Ipcomplex[c + 1] - 1; }
		else { s_to = prm->Nstrand; }
		for (s = s_from; s < s_to; s++) /* for each strand */
		{
			/* initialize strands */
			natoms_hcp2 = 0;
			j3 = 3 * s;
			x_hcp2[j3+0] = x_hcp2[j3+1] = x_hcp2[j3+2] = 0.0;

			for (m = 0; m < hcp; m++) 
			{ 
				j4 = 4 * (s * hcp + m);
				q_hcp2[j4+0] = q_hcp2[j4+1] = 
					q_hcp2[j4+2] = q_hcp2[j4+3] = 0.0;
			}

			r_from = prm->Ipstrand[s] - 1;
			if (s < prm->Nstrand - 1) { r_to = prm->Ipstrand[s + 1] - 1; }
			else { r_to = prm->Nres; }
			for (r = r_from; r < r_to; r++)         /* for each residue */
			{
				/* initialize residues */
				natoms_hcp1 = 0;
				j3 = 3 * r;
				x_hcp1[j3+0] = x_hcp1[j3+1] = x_hcp1[j3+2] = 0.0;

				for (m = 0; m < hcp; m++) 
				{ 
					j4 = 4 * (r * hcp + m);
					q_hcp1[j4+0] = q_hcp1[j4+1] = 
						q_hcp1[j4+2] = q_hcp1[j4+3] = 0.0;
                             //   printf("charges   %f    %f\n", q_hcp1[4 * (0 * hcp + 0)+0], q_hcp1[4 * (0 * hcp + 1)+0]); 
				}

				a_from = prm->Ipres[r];
				if (r + 1 < prm->Nres) { a_to = prm->Ipres[r + 1]; }
				else { a_to = prm->Natom; }
				for (a = a_from; a < a_to; a++)    /* for each atoms */
				{
					q = prm->Charges[a];
					m = q2idx(q);
					++natoms_hcp1;

					/* add to geometric center */
					j3 = 3 * r;
					j3a = 3 * a;
					x_hcp1[j3+0] += x[j3a+0];
					x_hcp1[j3+1] += x[j3a+1];
					x_hcp1[j3+2] += x[j3a+2];
					/* add to center of charge */
					j4 = 4 * (r * hcp + m);
					q_hcp1[j4+0] += x[j3a+0] * q;
					q_hcp1[j4+1] += x[j3a+1] * q;
					q_hcp1[j4+2] += x[j3a+2] * q;
					q_hcp1[j4+3] += q;

				}  /* end for each atoms */

				/* roll up to next level */
				natoms_hcp2 += natoms_hcp1;
				j3 = 3 * s;
				j3a = 3 * r;
				x_hcp2[j3+0] += x_hcp1[j3a+0];
				x_hcp2[j3+1] += x_hcp1[j3a+1];
				x_hcp2[j3+2] += x_hcp1[j3a+2];
				/* compute geometric center*/
				if (natoms_hcp1 > 0)
				{
					x_hcp1[j3a+0] /= natoms_hcp1;
					x_hcp1[j3a+1] /= natoms_hcp1;
					x_hcp1[j3a+2] /= natoms_hcp1;
				}
				for (m = 0; m < hcp; m++) 
				{ 
					/* roll up to next level */
					j4 = 4 * (s * hcp + m);
					j4a = 4 * (r * hcp + m);
					q_hcp2[j4+0] += q_hcp1[j4a+0];
					q_hcp2[j4+1] += q_hcp1[j4a+1];
					q_hcp2[j4+2] += q_hcp1[j4a+2];
					q_hcp2[j4+3] += q_hcp1[j4a+3];
					/* compute approx charge loc */
					q = q_hcp1[j4a+3];
					if (fabs(q) > MIN_Q)
					{  
						q_hcp1[j4a+0] /= q;
						q_hcp1[j4a+1] /= q;
						q_hcp1[j4a+2] /= q;
					}
					else
					{
						q_hcp1[j4a+0] = x_hcp1[j3a+0];
						q_hcp1[j4a+1] = x_hcp1[j3a+1];
						q_hcp1[j4a+2] = x_hcp1[j3a+2];
					}

				}
                                  
			}  /* end for each residue */

			/* roll up to next level */
			natoms_hcp3 += natoms_hcp2;
			j3 = 3 * c;
			j3a = 3 * s;
			x_hcp3[j3+0] += x_hcp2[j3a+0];
			x_hcp3[j3+1] += x_hcp2[j3a+1];
			x_hcp3[j3+2] += x_hcp2[j3a+2];
			/* compute geometric center */
			if (natoms_hcp2 > 0)
			{
				x_hcp2[j3a+0] /= natoms_hcp2;
				x_hcp2[j3a+1] /= natoms_hcp2;
				x_hcp2[j3a+2] /= natoms_hcp2;
			}
			for (m = 0; m < hcp; m++) 
			{ 
				/* roll up to next level */
				j4 = 4 * (c * hcp + m);
				j4a = 4 * (s * hcp + m);
				q_hcp3[j4+0] += q_hcp2[j4a+0];
				q_hcp3[j4+1] += q_hcp2[j4a+1];
				q_hcp3[j4+2] += q_hcp2[j4a+2];
				q_hcp3[j4+3] += q_hcp2[j4a+3];
				/* compute approx charge loc */
				q = q_hcp2[j4a+3];
				if (fabs(q) > MIN_Q)
				{  
					q_hcp2[j4a+0] /= q;
					q_hcp2[j4a+1] /= q;
					q_hcp2[j4a+2] /= q;
				}
				else
				{
					q_hcp2[j4a+0] = x_hcp2[j3a+0];
					q_hcp2[j4a+1] = x_hcp2[j3a+1];
					q_hcp2[j4a+2] = x_hcp2[j3a+2];
				}
			}
		}  /* end for each strand */
		/* compute geometric center and approx charges */
		if (natoms_hcp3 > 0)
		{
			j3a = 3 * c;
			x_hcp3[j3a+0] /= natoms_hcp3;
			x_hcp3[j3a+1] /= natoms_hcp3;
			x_hcp3[j3a+2] /= natoms_hcp3;
		}
		for (m = 0; m < hcp; m++) 
		{ 
			j4a = 4 * (c * hcp + m);
			q = q_hcp3[j4a+3];
			if (fabs(q) > MIN_Q)
			{  
				q_hcp3[j4a+0] /= q;
				q_hcp3[j4a+1] /= q;
				q_hcp3[j4a+2] /= q;
			}
			else
			{
				q_hcp3[j4a+0] = x_hcp3[j3a+0];
				q_hcp3[j4a+1] = x_hcp3[j3a+1];
				q_hcp3[j4a+2] = x_hcp3[j3a+2];
			}
		}
	} /* end for each complex */   

}   /*end calc_approx_q() */


/***********************************************************************
 *                          calc_approx_q_hcpo()
 *
 * Calculate geometric center of components and Optimal Point Charges
 *
 * Calling Parameters:
 * x      - atom coordinates                 (input)
 * x_hcp1 - geometric centers of residues    (update)
 * x_hcp2 - geometric centers of strands     (update)
 * q_hcp1 - approximate charges for residues (update)
 * q_hcp2 - approximate charges for strands  (update)
 * prm    - parmeter/topology structure      (input)
 ************************************************************************/
void calc_approx_q_hcpo(REAL_T * x, REAL_T * x_hcp1, REAL_T * x_hcp2, REAL_T * x_hcp3,
		REAL_T * q_hcp1, REAL_T * q_hcp2, REAL_T * q_hcp3, INT_T hcp, PARMSTRUCT_T * prm)
{

	int c, s, r, a, m, j3, j4, j3a, j4a, ja, i, j;          /* index for strand, residue, atom, charge */
	int a_from, a_to;                  			/* for atoms within a residue */
	int r_from, r_to;		      			/* for residues within a strand */
	int s_from, s_to;		      			/* for strands within a complex */
	int natoms_hcp1, natoms_hcp2, natoms_hcp3;      	/* total number of atoms */
	REAL_T q,qh;                          			/* charge */
        float  qt_hcp1,qt_hcp2,qt_hcp3;                         /* monopole */  
        float  p_hcp1[3], p_hcp2[3], p_hcp3[3];                 /* dipole */       
        float  Q_hcp1[3][3], Q_hcp2[3][3], Q_hcp3[3][3];        /* quadrupole */   
        float  atom_r2,rmax,rmax2,rmax2i, rcog[3];              /* size of the distribution,center of geometry*/     
        float qbar1,qbar2,rbar1[3],rbar2[3];                    /* optimal point charges, postions and values */
        float R1, R2;              /* thershold distances to constraint charges*/
        float R1h1=1.8, R2h1=2.0;  /* for hcp1: R1=R1h1*rmax, R2=R2h1*rmax*/
        float R1h2=0.8, R2h2=0.9;  /* for hcp2: R1=R1h2*rmax, R2=R2h2*rmax*/
        float R1h3=0.8, R2h3=0.9;  /* for hcp3: R1=R1h3*rmax, R2=R2h3*rmax*/
       // float dhcp = 0.25;         /* determines distance between opca charges in charged cases in hcp */ 
    // c complex
    // s strand
    // r residue 
    // a atom
    // j3a > pointer to the original charge distribution x,y,z
    // j3 > pointer to the hcp approximation (center of geometry):
    //           complex >  (x_hcp3,y_hcp3,z_hcp3)
    //           strand  >  (x_hcp2,y_hcp2,z_hcp2)
         //           residue >  (x_hcp1,y_hcp1,z_hcp1)
         // hcp > number of point charges for OPCA 
        // j4 pointer to the hcp approximations (center of charge)
        //             complex > (q_hcp3, 1 to 4, and m=0 to hcp) 
        //             strand  > (q_hcp2, 1 to 4, ") 
        //             residue > (q_hcp1, 1 to 4, ") 

        if (gb>0) {
          R1h1 = 0.8; // times Rmax gives the first threshold for the smoothing function for hcp1
          R2h1 = 0.9; // times Rmax gives the second threshold for the smoorhting function for hcp1 
        } 

	for (c = 0; c < prm->Ncomplex; c++) /* for each complex */
	{
		/* initialize complexes */
		natoms_hcp3 = 0;
		j3 = 3 * c;
		x_hcp3[j3+0] = x_hcp3[j3+1] = x_hcp3[j3+2] = 0.0;
        	qt_hcp3 = 0.0;                                       // monopole
        	p_hcp3[0]    = p_hcp3[1]    = p_hcp3[2]    = 0.0;    // dipole
        	Q_hcp3[0][0] = Q_hcp3[0][1] = Q_hcp3[0][2] = 0.0;    // quadrupole
        	Q_hcp3[1][0] = Q_hcp3[1][1] = Q_hcp3[1][2] = 0.0; 
        	Q_hcp3[2][0] = Q_hcp3[2][1] = Q_hcp3[2][2] = 0.0; 
		for (m = 0; m < hcp; m++) 
		{ 
			j4 = 4 * (c * hcp + m);
			q_hcp3[j4] = q_hcp3[j4+1] = q_hcp3[j4+2] = q_hcp3[j4+3] = 0.0;
		}
		s_from = prm->Ipcomplex[c] - 1;
		if (c < prm->Ncomplex - 1) { s_to = prm->Ipcomplex[c + 1] - 1; }
		else { s_to = prm->Nstrand; }
		for (s = s_from; s < s_to; s++) /* for each strand */
		{
			/* initialize strands */
			natoms_hcp2 = 0;
			j3 = 3 * s;
			x_hcp2[j3+0] = x_hcp2[j3+1] = x_hcp2[j3+2] = 0.0;
        	        qt_hcp2 = 0.0;                                    // monopole
                        p_hcp2[0]    = p_hcp2[1]    = p_hcp2[2] = 0.0;    // dipole
                        Q_hcp2[0][0] = Q_hcp2[0][1] = Q_hcp2[0][2] = 0.0; // quadrupole
                        Q_hcp2[1][0] = Q_hcp2[1][1] = Q_hcp2[1][2] = 0.0; 
                        Q_hcp2[2][0] = Q_hcp2[2][1] = Q_hcp2[2][2] = 0.0; 
			for (m = 0; m < hcp; m++) 
			{ 
				j4 = 4 * (s * hcp + m);
				q_hcp2[j4+0] = q_hcp2[j4+1] = 
				q_hcp2[j4+2] = q_hcp2[j4+3] = 0.0;
			}

			r_from = prm->Ipstrand[s] - 1;
			if (s < prm->Nstrand - 1) { r_to = prm->Ipstrand[s + 1] - 1; }
			else { r_to = prm->Nres; }
			for (r = r_from; r < r_to; r++)         /* for each residue */
			{
				/* initialize residues */
				natoms_hcp1 = 0;
				j3 = 3 * r;
				x_hcp1[j3+0] = x_hcp1[j3+1] = x_hcp1[j3+2] = 0.0;
                		qt_hcp1 = 0.0;                                    // monopole
                		p_hcp1[0]    = p_hcp1[1]    = p_hcp1[2] = 0.0;    // dipole
                		Q_hcp1[0][0] = Q_hcp1[0][1] = Q_hcp1[0][2] = 0.0; // quadrupole
                		Q_hcp1[1][0] = Q_hcp1[1][1] = Q_hcp1[1][2] = 0.0; 
                		Q_hcp1[2][0] = Q_hcp1[2][1] = Q_hcp1[2][2] = 0.0; 

				for (m = 0; m < hcp; m++) 
				{ 
					j4 = 4 * (r * hcp + m);
					q_hcp1[j4+0] = q_hcp1[j4+1] = 
					q_hcp1[j4+2] = q_hcp1[j4+3] = 0.0;
				}

				a_from = prm->Ipres[r];
				if (r + 1 < prm->Nres) { a_to = prm->Ipres[r + 1]; }
				else { a_to = prm->Natom; }
				for (a = a_from; a < a_to; a++)    /* for each atoms */
				{
					q = prm->Charges[a];
					m = q2idx(q);
					++natoms_hcp1;

					/* add to geometric center */
					j3 = 3 * r;
					j3a = 3 * a;
					x_hcp1[j3+0] += x[j3a+0];
					x_hcp1[j3+1] += x[j3a+1];
					x_hcp1[j3+2] += x[j3a+2];
					/* add to center of charge */
					j4 = 4 * (r * hcp + m);
                    			/* add to monopole*/
                   			qt_hcp1 += q;             
                    			/* add to dipole*/
                  			p_hcp1[0] += x[j3a+0] * q; 
                    			p_hcp1[1] += x[j3a+1] * q; 
                   			p_hcp1[2] += x[j3a+2] * q; 
                    			/* add to quadrupole */
                   			atom_r2=x[j3a+0]*x[j3a+0]+x[j3a+1]*x[j3a+1]
                                        +x[j3a+2]*x[j3a+2];
                    			qh = 0.5*q; 
                    			Q_hcp1[0][0] += qh*(3.0*x[j3a+0]*x[j3a+0]-atom_r2); 
                    			Q_hcp1[1][1] += qh*(3.0*x[j3a+1]*x[j3a+1]-atom_r2); 
                    			Q_hcp1[2][2] += qh*(3.0*x[j3a+2]*x[j3a+2]-atom_r2); 
                    			Q_hcp1[0][1] += qh*(3.0*x[j3a+0]*x[j3a+1]); 
                    			Q_hcp1[0][2] += qh*(3.0*x[j3a+0]*x[j3a+2]); 
                    			Q_hcp1[1][2] += qh*(3.0*x[j3a+1]*x[j3a+2]); 
               
				}  /* end for each atoms */
                		Q_hcp1[1][0] =  Q_hcp1[0][1]; 
                		Q_hcp1[2][0] =  Q_hcp1[0][2]; 
                		Q_hcp1[2][1] =  Q_hcp1[1][2]; 

                 		/* acumulate the multipoles to next level */
                 		qt_hcp2 += qt_hcp1; 
                		for (i=0; i<3; i++) {
                        	   p_hcp2[i] += p_hcp1[i];
                 		}
                		for (i=0; i<3; i++) {
                  	  		for (j=i; j<3; j++) {
                                   		Q_hcp2[i][j] += Q_hcp1[i][j];
                    	   		}
                 		}
                		Q_hcp2[1][0] =  Q_hcp2[0][1];
                 		Q_hcp2[2][0] =  Q_hcp2[0][2];
                 		Q_hcp2[2][1] =  Q_hcp2[1][2];
                                        

				/* roll up to next level */
				natoms_hcp2 += natoms_hcp1;
				j3 = 3 * s;
				j3a = 3 * r;
				x_hcp2[j3+0] += x_hcp1[j3a+0];
				x_hcp2[j3+1] += x_hcp1[j3a+1];
				x_hcp2[j3+2] += x_hcp1[j3a+2];
				/* compute geometric center*/
				if (natoms_hcp1 > 0)
				{
					x_hcp1[j3a+0] /= natoms_hcp1;
					x_hcp1[j3a+1] /= natoms_hcp1;
					x_hcp1[j3a+2] /= natoms_hcp1;
				}
                		rcog[0] = x_hcp1[j3a+0];
                		rcog[1] = x_hcp1[j3a+1];
                		rcog[2] = x_hcp1[j3a+2];
                		rmax2 = 0.0;
                		rmax2i = 0.0;
                		for (a = a_from; a < a_to; a++) {
                          		ja = 3 * a;
                         		 rmax2 = (x[ja+0]-rcog[0])*(x[ja+0]-rcog[0])+
                                  	(x[ja+1]-rcog[1])*(x[ja+1]-rcog[1])+
                                  	(x[ja+2]-rcog[2])*(x[ja+2]-rcog[2]);
                          		if (rmax2 > rmax2i) { rmax2i = rmax2 ;}
                		}
                		rmax = sqrt(rmax2i);
                		R1 = R1h1*rmax; 
               			R2 = R2h1*rmax;
                		OPCA_K2(r, qt_hcp1, p_hcp1, Q_hcp1, rmax, rcog, 
                               &qbar1, rbar1, &qbar2, rbar2, R1, R2, dhcp);
                		j4a = 4 * (r * hcp + 0);
                		q_hcp1[j4a+3] = qbar1;
               			q_hcp1[j4a+0] = rbar1[0];
                		q_hcp1[j4a+1] = rbar1[1];
                		q_hcp1[j4a+2] = rbar1[2];
               			j4a = 4 * (r * hcp + 1);
                		q_hcp1[j4a+3] = qbar2; 
                		q_hcp1[j4a+0] = rbar2[0];
                		q_hcp1[j4a+1] = rbar2[1];
                		q_hcp1[j4a+2] = rbar2[2];
			}  /* end for each residue */
                                
			/* roll up to next level */
			natoms_hcp3 += natoms_hcp2;
			j3 = 3 * c;
			j3a = 3 * s;
			x_hcp3[j3+0] += x_hcp2[j3a+0];
			x_hcp3[j3+1] += x_hcp2[j3a+1];
			x_hcp3[j3+2] += x_hcp2[j3a+2];
			/* compute geometric center */
			if (natoms_hcp2 > 0)
			{
				x_hcp2[j3a+0] /= natoms_hcp2;
				x_hcp2[j3a+1] /= natoms_hcp2;
				x_hcp2[j3a+2] /= natoms_hcp2;
			}
           		qt_hcp3 += qt_hcp2; 
            		for (i=0; i<3; i++) {
                		 p_hcp3[i] += p_hcp2[i];
            		}
            		for (i=0; i<3; i++) {
                		for (j=i; j<3; j++) {
                      			Q_hcp3[i][j] += Q_hcp2[i][j];
                		}    
            		}
            		Q_hcp3[1][0] =  Q_hcp3[0][1];
            		Q_hcp3[2][0] =  Q_hcp3[0][2];
            		Q_hcp3[2][1] =  Q_hcp3[1][2];

            		// In this part the hcp=2 charges are replaced by OPCA2, level 2
            		rcog[0] = x_hcp2[j3a+0];
            		rcog[1] = x_hcp2[j3a+1];
            		rcog[2] = x_hcp2[j3a+2];
            		rmax2 = 0.0;
            		rmax2i = 0.0;
            		for (r = r_from; r < r_to; r++) {
                    		j3a = 3 * r; 
                    		rmax2 = (x_hcp1[j3a+0]-rcog[0])*(x_hcp1[j3a+0]-rcog[0])+
                            		(x_hcp1[j3a+1]-rcog[1])*(x_hcp1[j3a+1]-rcog[1])+
                            		(x_hcp1[j3a+2]-rcog[2])*(x_hcp1[j3a+2]-rcog[2]);
                   		if (rmax2 > rmax2i) { rmax2i = rmax2 ;}
             		}
             		rmax = sqrt(rmax2i);
             		R1 = R1h2*rmax; 
             		R2 = R2h2*rmax; 
             		OPCA_K2(s, qt_hcp2, p_hcp2, Q_hcp2, rmax, rcog, 
                        &qbar1, rbar1, &qbar2, rbar2, R1, R2, dhcp);
             		j4a = 4 * (s * hcp + 0);
            		q_hcp2[j4a+3] = qbar1;
	                q_hcp2[j4a+0] = rbar1[0];
 	                q_hcp2[j4a+1] = rbar1[1];
            	        q_hcp2[j4a+2] = rbar1[2];
           	        j4a = 4 * (s * hcp + 1);
             		q_hcp2[j4a+3] = qbar2; 
             		q_hcp2[j4a+0] = rbar2[0];
             		q_hcp2[j4a+1] = rbar2[1];
             		q_hcp2[j4a+2] = rbar2[2];
             		j4a = 4 * (s * hcp + 0);
             		j4a = 4 * (s * hcp + 1);
                                
		}  /* end for each strand */
		/* compute geometric center and approx charges */
		if (natoms_hcp3 > 0)
		{
			j3a = 3 * c;
			x_hcp3[j3a+0] /= natoms_hcp3;
			x_hcp3[j3a+1] /= natoms_hcp3;
			x_hcp3[j3a+2] /= natoms_hcp3;
		}

       		 // In this part the hcp=2 charges are replaced by OPCA2, level 3
        	rcog[0] = x_hcp3[j3a+0];
        	rcog[1] = x_hcp3[j3a+1];
        	rcog[2] = x_hcp3[j3a+2];
        	rmax2 = 0.0;
        	rmax2i = 0.0;
 
        	for (s = s_from; s < s_to; s++) {
                	j3a = 3 * s; 
               		rmax2 = (x_hcp2[j3a+0]-rcog[0])*(x_hcp2[j3a+0]-rcog[0])+
                        (x_hcp2[j3a+1]-rcog[1])*(x_hcp2[j3a+1]-rcog[1])+
                        (x_hcp2[j3a+2]-rcog[2])*(x_hcp2[j3a+2]-rcog[2]);
                	if (rmax2 > rmax2i) { rmax2i = rmax2 ;}
         	}
         	rmax = sqrt(rmax2i);
        	R1 = R1h3*rmax; 
        	R2 = R2h3*rmax; 
         	OPCA_K2(c, qt_hcp3, p_hcp3, Q_hcp3, rmax, rcog, 
                 &qbar1, rbar1, &qbar2, rbar2, R1, R2, dhcp);
 		 j4a = 4 * (c * hcp + 0);
         	q_hcp3[j4a+3] = qbar1;
         	q_hcp3[j4a+0] = rbar1[0];
         	q_hcp3[j4a+1] = rbar1[1];
         	q_hcp3[j4a+2] = rbar1[2];
		j4a = 4 * (c * hcp + 1);
         	q_hcp3[j4a+3] = qbar2; 
         	q_hcp3[j4a+0] = rbar2[0];
         	q_hcp3[j4a+1] = rbar2[1];
         	q_hcp3[j4a+2] = rbar2[2];
	} /* end for each complex */   

}   /*end calc_approx_q_hcpo() */



/***********************************************************************
 *                          CALC_ENERGY_APPROX() 
 *
 * Calculate force/energy for various levels of approximation
 *
 * Calling parameters:
 * npairs_hcpx - (x=0/1/2) number of pairs in each pairlist (input)
 * pairs_hcpx  - (x=0/1/2) list of atom/residue/chains      (input)
 * x           - atom coordinates                           (input)
 * q_hcpx      - (x=1/2) approximate charges                (input)
 * Various parameters from nbond_hcp()			    (input)

 ************************************************************************/

void calc_energy_approx(int j, REAL_T *q_hcp, REAL_T *x,int hcp, int i, REAL_T cgi, REAL_T *df, REAL_T *elec, REAL_T *evdw, REAL_T xi, REAL_T yi, REAL_T zi, REAL_T wi, 
		REAL_T *dumx, REAL_T *dumy, REAL_T *dumz, REAL_T *dumw)
{

	int k;
	REAL_T cgj, xij, yij, zij, r2, r2inv, r, rinv, rs, rssq, pow, eps1, epsi, cgijr, df2, dedx, dedy, dedz;

	for (k = 0; k < hcp; k++)
	{

		cgj = q_hcp[4*(j*hcp+k)+3];
		if (fabs(cgj) > MIN_Q)
		{ 
			xij = xi - q_hcp[4*(j*hcp+k)+0];
			yij = yi - q_hcp[4*(j*hcp+k)+1];
			zij = zi - q_hcp[4*(j*hcp+k)+2];
			r2 = xij * xij + yij * yij + zij * zij;

			r2inv = 1.0 / r2;
			r = sqrt(r2);
			rinv = r * r2inv;

			/* Calculate the energy and derivatives according to dield. */

			if (dield == -3) {

				/* special code Ramstein & Lavery dielectric, 94 force field */
				rs = SIG * r;
				rssq = rs * rs;
				pow = exp(-rs);
				eps1 = rssq + rs + rs + 2.0;
				epsi = 1.0 / (DIW - C1 * pow * eps1);
				cgijr = cgi * cgj * rinv * epsi;
				*elec += cgijr;
				*df = df2 * rinv;

			} else if (dield == -4) {

				/* distance-dependent dielectric code, 94 ff */
				/* epsilon = r  */
				rs = (cgi * cgj * r2inv);
				df2 = -2.0 * rs;
				*elec += rs;
				*df = df2 * rinv;


			} else if (dield == -5) {

				/* non-bonded term from yammp  */
				*df = 0.0;

			} else {

				/*
				 * Code for various dielectric models.
				 * The df2 variable should hold r(dV/dr).
				 */

				if (dield == 0) {

					/* epsilon = r  */

					rs = cgi * cgj * r2inv;
					df2 = -2.0 * rs;
					*elec += rs;

				} else if (dield == 1) {

					/* epsilon = 1  */
					rs = (cgi * cgj * rinv);
					df2 = -rs;
					*elec += rs;

				} else if (dield == -2) {

					/* Ramstein & Lavery dielectric, PNAS 85, 7231 (1988). */

					rs = SIG * r;
					rssq = rs * rs;
					pow = exp(-rs);
					eps1 = rssq + rs + rs + 2.0;
					epsi = 1.0 / (DIW - C1 * pow * eps1);
					cgijr = cgi * cgj * rinv * epsi;
					*elec += cgijr;
					df2 = -cgijr * (1.0 + C1 * pow * rs * rssq * epsi);
				}

				/* no vdw force calculated for residue/chain appr */
				*df = df2 * rinv;
			}


			/*
			 * The de term contains one more factor of Dij in the denominator
			 * so that terms such as dedx do not need to include 1/Dij. 
			 *
			 * Update the gradient for atom j.
			 */
			*df *= rinv;

			dedx = *df * xij;
			dedy = *df * yij;
			dedz = *df * zij;

			*dumx += dedx;
			*dumy += dedy;
			*dumz += dedz;

		}  /* end-if q > MIN_Q */

	}  /* end of for each approximate charge */
}

/***********************************************************************
 *                          CALC_ENERGY_ATOM() 
 *
 * Calculate force/energy for various levels of approximation
 *
 * Calling parameters:
 * npairs_hcpx - (x=0/1/2) number of pairs in each pairlist (input)
 * pairs_hcpx  - (x=0/1/2) list of atom/residue/chains      (input)
 * x           - atom coordinates                           (input)
 * iexw        - excluded pairlist                          (input)
 * Various parameters from nbond_hcp()			    (input)

 ************************************************************************/
void calc_energy_atom(int j, int *iexw, PARMSTRUCT_T *prm, REAL_T *x, int i, int iaci,
		REAL_T cgi, REAL_T *df, REAL_T *elec, REAL_T *evdw, 
		REAL_T xi, REAL_T yi, REAL_T zi, REAL_T wi, REAL_T *dumx, REAL_T *dumy, REAL_T *dumz, REAL_T *dumw)
{

	int ic, ibig, isml;
	REAL_T xij, yij, zij, r2, r2inv, r, rinv, rs, rssq, pow, eps1, epsi, cgijr, df2, r6, r10, f1, f2, dis, d0, kij, diff, dedx, dedy, dedz;

	if (iexw[j] != i)
	{
		xij = xi - x[dim * j + 0];
		yij = yi - x[dim * j + 1];
		zij = zi - x[dim * j + 2];

		/*             FROM HERE TO END OF IF MOVE TO A SUBROUTINE! */

		r2 = xij * xij + yij * yij + zij * zij;

		r2inv = 1.0 / r2;
		r = sqrt(r2);
		rinv = r * r2inv;

		/* Calculate the energy and derivatives according to dield. */

		if (dield == -3) {

			/* special code Ramstein & Lavery dielectric, 94 force field */

			rs = SIG * r;
			rssq = rs * rs;
			pow = exp(-rs);
			eps1 = rssq + rs + rs + 2.0;
			epsi = 1.0 / (DIW - C1 * pow * eps1);
			cgijr = cgi * prm->Charges[j] * rinv * epsi;
			*elec += cgijr;
			df2 = -cgijr * (1.0 + C1 * pow * rs * rssq * epsi);
			ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
			if (ic >= 0) {
				r6 = r2inv * r2inv * r2inv;
				f2 = prm->Cn2[ic] * r6;
				f1 = prm->Cn1[ic] * r6 * r6;
				*evdw += (f1 - f2);
				*df = (df2 + (6.0 * f2 - 12.0 * f1)) * rinv;
			} else {
				*df = df2 * rinv;
			}

		} else if (dield == -4) {

			/* distance-dependent dielectric code, 94 ff */
			/* epsilon = r  */

			rs = cgi * prm->Charges[j] * r2inv;
			df2 = -2.0 * rs;
			*elec += rs;
			ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
			if (ic >= 0) {
				r6 = r2inv * r2inv * r2inv;
				f2 = prm->Cn2[ic] * r6;
				f1 = prm->Cn1[ic] * r6 * r6;
				*evdw += (f1 - f2);
				*df = (df2 + (6.0 * f2 - 12.0 * f1)) * rinv;
			} else {
				*df = df2 * rinv;
			}

		} else if (dield == -5) {

			/* non-bonded term from yammp  */       

			dis = r;
			ic = prm->Cno[iaci + prm->Iac[j] - 1] - 1;
			d0 = prm->Cn2[ic];
			if (dis < d0) {
				kij = prm->Cn1[ic];
				diff = dis - d0;
				*evdw += kij * diff * diff;
				*df = 2.0 * kij * diff;
			} else {
				*df = 0.0;
			}
		} else {

			/*
			 * Code for various dielectric models.
			 * The df2 variable should hold r(dV/dr).
			 */

			if (dield == 0) {

				/* epsilon = r  */

				rs = cgi * prm->Charges[j] * r2inv;
				df2 = -2.0 * rs;
				*elec += rs;


			} else if (dield == 1) {

				/* epsilon = 1  */
				rs = (cgi * prm->Charges[j] * rinv);
				df2 = -rs;
				*elec += rs;

			} else if (dield == -2) {

				/* Ramstein & Lavery dielectric, PNAS 85, 7231 (1988). */

				rs = SIG * r;
				rssq = rs * rs; 
				pow = exp(-rs);
				eps1 = rssq + rs + rs + 2.0;
				epsi = 1.0 / (DIW - C1 * pow * eps1);
				cgijr = cgi * prm->Charges[j] * rinv * epsi;
				*elec += cgijr;
				df2 = -cgijr * (1.0 + C1 * pow * rs * rssq * epsi);
			}

			/* Calculate either Van der Waals or hydrogen bonded term. */
			ic =prm->Cno[iaci + prm->Iac[j] - 1];
			if (ic > 0) 
			{

				if (ic > 0) 
				{
					ic--;
				} 
				else 
				{
					ibig = prm->Iac[i] > prm->Iac[j] ?
						prm->Iac[i] : prm->Iac[j];
					isml = prm->Iac[i] > prm->Iac[j] ?
						prm->Iac[j] : prm->Iac[i];
					ic = ibig * (ibig - 1) / 2 + isml - 1;
				}
				r6 = r2inv * r2inv * r2inv;
				f2 = prm->Cn2[ic] * r6;
				f1 = prm->Cn1[ic] * r6 * r6;
				*evdw += (f1 - f2);
				*df = (df2 + (6.0 * f2 - 12.0 * f1)) * rinv;
			} 
			else {
				ic = -ic - 1;
				r10 = r2inv * r2inv * r2inv * r2inv * r2inv;
				f2 = prm->HB10[ic] * r10;
				f1 = prm->HB12[ic] * r10 * r2inv;
				*evdw += (f1 - f2);
				*df = (df2 + (10.0 * f2 - 12.0 * f1)) * rinv;
#if 0
				hbener += (f1 - f2);
#endif
			}
		}

		/*
		 * The de term contains one more factor of Dij in the denominator
		 * so that terms such as dedx do not need to include 1/Dij. 
		 *
		 * Update the gradient for atom j.
		 */

		*df *= rinv;

		dedx = *df * xij;
		dedy = *df * yij;
		dedz = *df * zij;

		*dumx += dedx;
		*dumy += dedy;
		*dumz += dedz;

	}  /* end of if not on exlcuded atom list */
}

/***********************************************************************
 *                          NBOND_HCP() 
 *
 * Calculate non-bonded energy and first derivative
 *
 * Calling parameters:
 * npairs_hcpx - (x=0/1/2) number of pairs in each pairlist (input)
 * pairs_hcpx  - (x=0/1/2) list of atom/residue/chains      (input)
 * Iblo_hcp    - number of excluded atoms                   (input)
 * IexclAt_hcp - list of excluded atoms                     (input)
 * x           - atom coordinates                           (input)
 * q_hcpx      - (x=1/2) approximate charges                (input)
 * f           - first derivative (gradient) vector         (update)
 * enb         - Van der Waals energy                       (update)
 * eel         - Coulombic energy                           (update)
 * 
 ************************************************************************/

int nbond_hcp(int * Iblo_hcp, int ** IexclAt_hcp, REAL_T * x, REAL_T *x_hcp1, REAL_T *x_hcp2, REAL_T *x_hcp3, REAL_T * q_hcp1, REAL_T * q_hcp2, REAL_T * q_hcp3, REAL_T * f, REAL_T * enb, REAL_T * eel)
{
	int i, j, iaci, foff, threadnum, numthreads;	
	int *iexw;
	REAL_T dumx, dumy, dumz, dumw, cgi;
	REAL_T df, evdw, elec;
	REAL_T xi, yi, zi, wi;
	REAL_T dist2;
        REAL_T dist2_hcp1, dist2_hcp2, dist2_hcp3;
	int c, s, res, a;  	          				/* index for complex, strand, residue and atom */
	int a_from, a_to;		          			/* for atoms within a residue */
	int r_from, r_to;			      			/* for residues withnin a strand */
	int s_from, s_to;			      			/* for strands within a complex */			
	dist2_hcp1 = hcp_h1*hcp_h1;
	dist2_hcp2 = hcp_h2*hcp_h2;   
	dist2_hcp3 = hcp_h3*hcp_h3;

	evdw = 0.0;
	elec = 0.0;

	/*
	 * thread number and the number of threads for multi-threaded
	 * execution under OpenMP.  For all other cases, including ScaLAPACK,
	 * MPI and single-threaded execution, use the values that have been
	 * stored in mytaskid and numtasks, respectively.
	 */

	threadnum = mytaskid;
	numthreads = numtasks;
	foff = 0;

	/* calculate approximate charges and geometric centers */
    	if (hcpo==0) {
       	   calc_approx_q(x, x_hcp1, x_hcp2, x_hcp3, q_hcp1, q_hcp2, q_hcp3, hcp, prm);
	} else {
	   calc_approx_q_hcpo(x, x_hcp1, x_hcp2, x_hcp3, q_hcp1, q_hcp2, q_hcp3, hcp, prm);
    	}

	/* iexw stores the list of excluded atoms of each atom */
	iexw = ivector(-1, prm->Natom+1);
	for (i = -1; i < prm->Natom+1; i++) {
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
	 * Synchronization of MPI tasks will occur via the MPI_Allreduce
	 * function that is called from within mme34.
	 */

	for (i = 0; i < prm->Natom ; i++) {

#if defined(MPI)
		if (!myroc(i, blocksize, numthreads, threadnum))
			continue;
#endif

		iaci = prm->Ntypes * (prm->Iac[i] - 1);
		dumx = dumy = dumz = 0.0;
		xi = x[dim * i + 0];
		yi = x[dim * i + 1];
		zi = x[dim * i + 2];

		cgi = prm->Charges[i];

		/*
		 * Expand the excluded list into the iexw array 
		 * by storing i at array address j
		 */

		for (j = 0; j < Iblo_hcp[i]; j++) {
			iexw[IexclAt_hcp[i][j] - 1] = i;
		} 


		for(c=0; c<prm->Ncomplex; c++)
		{
			dist2 = calc_dist2(xi, yi, zi, x_hcp3[3*c+0], x_hcp3[3*c+1], x_hcp3[3*c+2]);
			if(dist2 > dist2_hcp3)
			{
				calc_energy_approx(c, q_hcp3, x, hcp, i, cgi, &df, &elec, &evdw, xi, yi, zi, wi, &dumx, &dumy, &dumz, &dumw);
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
						calc_energy_approx(s, q_hcp2, x, hcp, i, cgi, &df, &elec, &evdw, xi, yi, zi, wi, &dumx, &dumy, &dumz, &dumw);
					}
					else
					{
						r_from = prm->Ipstrand[s] - 1;
						if (s + 1 < prm->Nstrand) { r_to = prm->Ipstrand[s + 1] - 1; }
						else { r_to = prm->Nres; }

						for (res = r_from; res < r_to; res++)         /* for each residue in the strand */
						{
							dist2 = calc_dist2(xi, yi, zi,
									x_hcp1[3*res+0], x_hcp1[3*res+1], x_hcp1[3*res+2]);
							if (dist2 > dist2_hcp1)
							{
								calc_energy_approx(res, q_hcp1, x, hcp, i, cgi, &df, &elec, &evdw, xi, yi, zi, wi, &dumx, &dumy, &dumz, &dumw);
							}
							else
							{
								a_from = prm->Ipres[res];
								if (res + 1 < prm->Nres) { a_to = prm->Ipres[res + 1]; }
								else { a_to = prm->Natom; }

								for (a = a_from; a < a_to; a++)     /* for each atom in the atom */
								{
									if (a != i)   /* not same atom */
									{
										calc_energy_atom(a, iexw, prm, x, i, iaci, cgi, &df, &elec, &evdw, xi, yi, zi, wi, &dumx, &dumy, &dumz, &dumw);
									}   /* end if not same atom */
								}   /* end for each atom a */
							}   /* end else in residue threshold dist */
						}   /* end for each residue r */
					}   /* end else in strand threshod dist */
				}   /* end for each strand s */
			}   /* end else in complex threshold dist */
		}   /* end for each complex c */

		f[foff + dim * i + 0] += dumx;
		f[foff + dim * i + 1] += dumy;
		f[foff + dim * i + 2] += dumz;

	}  /* end-for each atom */

	/* Deallocate the iexw array within this potentially parallel region. */
	free_ivector(iexw, -1, prm->Natom+1);

	/* Return evdw and elec through by-reference calling parameters. */

	*enb = evdw/2.0; 
	*eel = elec/2.0; 
    
	return (0);
}

