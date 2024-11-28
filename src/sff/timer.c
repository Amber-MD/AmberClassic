/***********************************************************************
                            INIT_TIMERS()
************************************************************************/
void init_timers(void)
{
  int i;
  for(i=0; i<NTIMER; i++)
    timearr[i] = 0;

  tinit = &(timearr[0]);
  tmd = &(timearr[1]);
  tmdOther = &(timearr[2]);
  tmdIO = &(timearr[3]);
  tnmode = &(timearr[4]);
  tnmodeOther = &(timearr[5]);
  tnmodeFact = &(timearr[6]);
  tnmodeEigenvec = &(timearr[7]);
  tnmodeInvdiag = &(timearr[8]);
  tnmodeAmat = &(timearr[9]);
  tnmodeSort = &(timearr[10]);
  tnmodeNorm = &(timearr[11]);
  tnmodeLan = &(timearr[12]);
  tnmodeHessian = &(timearr[13]);
  tnmodeDSYEV = &(timearr[14]);
  tnmodeDSYEVD = &(timearr[15]);
  tnmodeDiagA = &(timearr[16]);
  tconjgrad = &(timearr[17]);
  tconjgradMME = &(timearr[18]);
  tconjgradOther = &(timearr[19]);
  tnewton = &(timearr[20]);
  tnewtonLevel = &(timearr[21]);
  tnewtonCholesky = &(timearr[22]);
  tnewtonMME = &(timearr[23]);
  tnewtonMME2 = &(timearr[24]);
  tnewtonOther = &(timearr[25]);
  tnewtonDSYEV = &(timearr[26]);
  tmme = &(timearr[27]);
  tmmeCons = &(timearr[28]);
  tmmeNonb = &(timearr[29]);
  tmmePair = &(timearr[30]);
  tmmeBond = &(timearr[31]);
  tmmeAngl = &(timearr[32]);
  tmmePhi = &(timearr[33]);
  tmmeBorn = &(timearr[34]);
  tmmeOther = &(timearr[35]);
  tmmePB = &(timearr[36]);
  tmmeRism = &(timearr[37]);
  tmme2 = &(timearr[38]);
  tmme2Cons = &(timearr[39]);
  tmme2Nonb = &(timearr[40]);
  tmme2Pair = &(timearr[41]);
  tmme2Bond = &(timearr[42]);
  tmme2Angl = &(timearr[43]);
  tmme2Phi = &(timearr[44]);
  tmme2Born = &(timearr[45]);
  tmme2Other = &(timearr[46]);
  tgb2 = &(timearr[47]);
  tgb2dgemm1 = &(timearr[48]);
  tgb2dgemm2 = &(timearr[49]);
  tgb2dgemm3 = &(timearr[50]);
  tgb2Other = &(timearr[51]);
  treduceegb = &(timearr[52]);
  treducemme = &(timearr[53]);
}

/***********************************************************************
                            MME_TIMER()
************************************************************************/

/* Print a timing summary but only for task zero. */

INT_T mme_timer(void)
{
   /* Use the maximum time from all MPI tasks or SCALAPACK processes. */

  int labelLen = 20;
#if defined(MPI) || defined(SCALAPACK)
  REAL_T reductarr[NTIMER];

  MPI_Allreduce(timearr, reductarr, NTIMER, MPI_DOUBLE, MPI_MAX,
		MPI_COMM_WORLD);
  
   *tinit = reductarr[0];
   *tmd = reductarr[1];
   *tmdOther = reductarr[2];
   *tmdIO = reductarr[3];
   *tnmode = reductarr[4];
   *tnmodeOther = reductarr[5];
   *tnmodeFact = reductarr[6];
   *tnmodeEigenvec = reductarr[7];
   *tnmodeInvdiag = reductarr[8];
   *tnmodeAmat = reductarr[9];
   *tnmodeSort = reductarr[10];
   *tnmodeNorm = reductarr[11];
   *tnmodeLan = reductarr[12];
   *tnmodeHessian = reductarr[13];
   *tnmodeDSYEV = reductarr[14];
   *tnmodeDSYEVD = reductarr[15];
   *tnmodeDiagA = reductarr[16];
   *tconjgrad = reductarr[17];
   *tconjgradMME = reductarr[18];
   *tconjgradOther = reductarr[19];
   *tnewton = reductarr[20];
   *tnewtonLevel = reductarr[21];
   *tnewtonCholesky = reductarr[22];
   *tnewtonMME = reductarr[23];
   *tnewtonMME2 = reductarr[24];
   *tnewtonOther = reductarr[25];
   *tnewtonDSYEV = reductarr[26];
   *tmme = reductarr[27];
   *tmmeCons = reductarr[28];
   *tmmeNonb = reductarr[29];
   *tmmePair = reductarr[30];
   *tmmeBond = reductarr[31];
   *tmmeAngl = reductarr[32];
   *tmmePhi = reductarr[33];
   *tmmeBorn = reductarr[34];
   *tmmeOther = reductarr[35];
   *tmmePB = reductarr[36];
   *tmmeRism = reductarr[37];
   *tmme2 = reductarr[38];
   *tmme2Cons = reductarr[39];
   *tmme2Nonb = reductarr[40];
   *tmme2Pair = reductarr[41];
   *tmme2Bond = reductarr[42];
   *tmme2Angl = reductarr[43];
   *tmme2Phi = reductarr[44];
   *tmme2Born = reductarr[45];
   *tmme2Other = reductarr[46];
   *tgb2 = reductarr[47];
   *tgb2dgemm1 = reductarr[48];
   *tgb2dgemm2 = reductarr[49];
   *tgb2dgemm3 = reductarr[50];
   *tgb2Other = reductarr[51];
   *treduceegb = reductarr[52];
   *treducemme = reductarr[53];
#endif

   if (mytaskid == 0) {
      fprintf(nabout, "\n|Timing summary:\n");
      /*Major sections with timings include md, normal mode, conjgrad, newton, lmod and xmin*/
      /* Top level */
      fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Initialize", *tinit);
      fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Molec. Dyn.", *tmd);
      fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Normal Mode", *tnmode);
      fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Conj. Grad.", *tconjgrad);
      fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Newton", *tnewton);
      /* fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Lmod", lmod_opt.lmod_time); */
      /* fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Xmin", xmin_opt.xmin_time); */
      fprintf(nabout, "|-------------------------\n");
      /* fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Total", *tinit + *tmd + *tnmode  */
      /* 	      + *tconjgrad + *tnewton + lmod_opt.lmod_time + xmin_opt.xmin_time); */
      fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Total", *tinit + *tmd + *tnmode 
	      + *tconjgrad + *tnewton);

      if(*tmd > 0){
	fprintf(nabout,"\n|Molecular Dynamics Timing Summary:\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "1st Derivative", *tmd-*tmdIO-*tmdOther);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "I/O", *tmdIO);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Other", *tmdOther);
	fprintf(nabout, "|-------------------------\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Total", *tmd);
      }
      
      if(*tnmode > 0){
	fprintf(nabout,"\n|Normal Mode Timing Summary:\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Build Hessian", 
		*tnmodeHessian);
	if(*tnmodeDSYEV > 0)
	  fprintf(nabout, "|    %*s  %10.3f\n", labelLen, 
		  "diagonalize w/ dsyev", *tnmodeDSYEV);
	if(*tnmodeDSYEVD > 0)
	  fprintf(nabout, "|    %*s  %10.3f\n", labelLen, 
		  "diagonalize w/ dsyevd", *tnmodeDSYEVD);
	if(*tnmodeFact > 0){
	  fprintf(nabout, "|    %*s  %10.3f\n", labelLen, 
		  "Cholesky factorization", *tnmodeFact);
	  fprintf(nabout, "|    %*s  %10.3f\n", labelLen, 
		  "diagonalize inverted Hessian", *tnmodeInvdiag);
	  fprintf(nabout, "|    %*s  %10.3f\n", labelLen, 
		  "Hessian eigenvectors", *tnmodeEigenvec);
	}
	if(*tnmodeAmat > 0){
	  fprintf(nabout, "|    %*s  %10.3f\n", labelLen, 
		  "Setup 'a' matrix", *tnmodeAmat);
	  fprintf(nabout, "|    %*s  %10.3f\n", labelLen, 
		  "Diagonalize 'a' matrix", *tnmodeDiagA);
	  fprintf(nabout, "|    %*s  %10.3f\n", labelLen, 
		  "Sort eigenvalues", *tnmodeSort);
	  fprintf(nabout, "|    %*s  %10.3f\n", labelLen, 
		  "Normalize 'L' matrix", *tnmodeNorm);
	  fprintf(nabout, "|    %*s  %10.3f\n", labelLen, 
		  "Langevin modes", *tnmodeLan);
	}
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Other", *tnmodeOther);
	fprintf(nabout, "|-------------------------\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Total", *tnmode);
	
      }

      if(*tconjgrad > 0){
	fprintf(nabout, "\n|Conjugate gradient timing summary:\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "1st Derivative", *tconjgradMME);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Other", *tconjgradOther);
	fprintf(nabout, "|-------------------------\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Total", *tconjgrad);
	/* fprintf(nabout, "|    %*s  %10.3e\n", labelLen, "Relative Error", */
	/* 	(*tconjgrad-*tconjgradMME-*tconjgradOther)/ *tconjgrad); */
      }
      
      if(*tnewton > 0){
	fprintf(nabout, "\n|Newton-Rhapson timing summary:\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "1st Derivative", *tnewtonMME);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "1st & 2nd Derivative", *tnewtonMME2);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Level", *tnewtonLevel);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "dposv", *tnewtonCholesky);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "dsyev", *tnewtonDSYEV);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Other", *tnewtonOther);
	fprintf(nabout, "|-------------------------\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Total", *tnewton);
	/* fprintf(nabout, "|    %*s  %10.3e\n", labelLen, "Relative Error",  */
	/* 	(*tnewton-*tnewtonMME-*tnewtonMME2-*tnewtonLevel-*tnewtonCholesky-*tnewtonDSYEV-*tnewtonOther) / *tnewton); */
      }

      /* energy and gradients */
      /* mme */
      if(*tmme > 0){
	fprintf(nabout, "\n|1st derivative timing summary:\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "constraints", *tmmeCons);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "bonds", *tmmeBond);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "angles", *tmmeAngl);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "torsions", *tmmePhi);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "pairlist", *tmmePair);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "nonbond", *tmmeNonb);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "gen. Born", *tmmeBorn);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Poisson Boltzmann", *tmmePB);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "3D-RISM", *tmmeRism);
#ifdef MPI
	/* PBSA is serial only so wait time on idle nodes ends up in 
	   tmmeOther.  As a simple work around we take the difference*/ 
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Other", *tmmeOther-*tmmePB);
#else
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Other", *tmmeOther);
#endif
	fprintf(nabout, "|-------------------------\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Total", *tmme);
	/* fprintf(nabout, "|    %*s  %10.3e\n", labelLen, "Relative Error",  */
	/* 	(*tmme-*tmmeCons-*tmmeBond-*tmmeAngl-*tmmePhi-*tmmePair-*tmmeNonb-*tmmeBorn-*tmmeRism-*tmmeOther) / *tmme); */
      }
      
      /* mme2 */
      if(*tmme2 > 0){
	fprintf(nabout, "\n|1st & 2nd derivative timing summary:\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "constraints", *tmme2Cons);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "bonds", *tmme2Bond);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "angles", *tmme2Angl);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "torsions", *tmme2Phi);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "pairlist", *tmme2Pair);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "nonbond", *tmme2Nonb);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "gen. Born", *tmme2Born);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Other", *tmme2Other);
	fprintf(nabout, "|-------------------------\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Total", *tmme2);
	/* fprintf(nabout, "|    %*s  %10.3e\n", labelLen, "Relative Error",  */
	/* 	(*tmme2-*tmme2Cons-*tmme2Bond-*tmme2Angl-*tmme2Phi-*tmme2Pair-*tmme2Nonb-*tmme2Born-*tmme2Other) / *tmme2); */
      }

      /* EGB2 */
      if(*tgb2 > 0){
	fprintf(nabout, "\n|Generalized Born 1st & 2nd derivative timing summary:\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "dgemm1", *tgb2dgemm1);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "dgemm2", *tgb2dgemm2);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "dgemm3", *tgb2dgemm3);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Other", *tgb2Other);
	fprintf(nabout, "|-------------------------\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "Total", *tgb2);
	/* fprintf(nabout, "|    %*s  %10.3e\n", labelLen, "Relative Error", (*tgb2-*tgb2dgemm1-*tgb2dgemm2-*tgb2dgemm3-*tgb2Other) / *tgb2); */
      }

      /* MPI */
      if(*treduceegb + *treducemme > 0){
	fprintf(nabout, "\n|MPI timing summary:\n");
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "egb reduct.", *treduceegb);
	fprintf(nabout, "|    %*s  %10.3f\n", labelLen, "mme reduct.", *treducemme);
      }
      fflush(nabout);
   }
   /* 3D-RISM */
#ifdef RISMSFF   
   if(*tmmeRism > 0 && rismData.rism){
     /* use the built in timer report */
     rism_printtimer_(); 
   }
#endif /*RISMSFF*/
   return (0);
}
