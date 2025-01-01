#if NEEDFORCE == 1
  #if NEEDENERGY == 0
    #if LISTSAME == 1
      #define ONLYQQINTR OnlyQQSame
    #else
      #define ONLYQQINTR OnlyQQIntr
    #endif
  #else
    #if NEEDVIRIAL == 0
      #if LISTSAME == 1
        #define ONLYQQINTR OnlyQQSamenrg
      #else
        #define ONLYQQINTR OnlyQQIntrnrg
      #endif
    #else
      #if LISTSAME == 1
        #define ONLYQQINTR OnlyQQSamenrgvir
      #else
        #define ONLYQQINTR OnlyQQIntrnrgvir
      #endif
    #endif
  #endif
#else
  #if LISTSAME == 1
    #define ONLYQQINTR OnlyQQSamenrgx
  #else
    #define ONLYQQINTR OnlyQQIntrnrgx
  #endif
#endif

#if NEEDFORCE == 1
#if NEEDENERGY == 0
#if LISTSAME == 1
//-----------------------------------------------------------------------------
// ONLYQQINTR: this specialized routine will compute only electrostatic
//             interactions between particles.  See the BOTHINTR function
//             description in NBLoop.c for definitions of its arguments.
//-----------------------------------------------------------------------------
static void ONLYQQINTR(atomc *atm1, int *ordr1,  int str1,  int str2,
		       int lim1,  int lim2, int* qIDbuff, rngbuff* qr2buff,
		       dircon *dcinp, FrcTab *EFrc)
#else
static void ONLYQQINTR(atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		       int str1,  int str2,  int lim1, int lim2, int* qIDbuff,
		       rngbuff* qr2buff, dircon *dcinp, FrcTab *EFrc)
#endif
#else
#if LISTSAME == 1
static void ONLYQQINTR(atomc *atm1, int *ordr1,  int str1,  int str2,
		       int lim1,  int lim2, int* qIDbuff, rngbuff* qr2buff,
		       dircon *dcinp, FrcTab *EFrc, Energy *sysUV)
#else
static void ONLYQQINTR(atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		       int str1,  int str2,  int lim1, int lim2, int* qIDbuff,
		       rngbuff* qr2buff, dircon *dcinp, FrcTab *EFrc,
		       Energy *sysUV)
#endif
#endif
#else
#if LISTSAME == 1
static void ONLYQQINTR(atomc *atm1, int *ordr1,  int str1,  int str2,
		       int lim1,  int lim2, int* qIDbuff, rngbuff* qr2buff,
		       dircon *dcinp, FrcTab *EFrc, Energy *sysUV)
#else
static void ONLYQQINTR(atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		       int str1,  int str2,  int lim1, int lim2, int* qIDbuff,
		       rngbuff* qr2buff, dircon *dcinp, FrcTab *EFrc,
		       Energy *sysUV)
#endif
#endif
{
  int i, j, nqpair;

  /*** Pre-compute a bunch of constants and unpack structures ***/
  const double uspc = EFrc->ivdr;
  const double Ecut2 = dcinp->Ecut*dcinp->Ecut;
#if NEEDFORCE == 1
  double aifx, aify, aifz;
  double fxA, fyA, fzA, fxB, fyB, fzB, fxC, fyC, fzC, fmagA, fmagB, fmagC;
  double fxD, fyD, fzD, fmagD;
  CSpln* Fspl;
  Fspl = EFrc->dSD;
#endif
#if NEEDENERGY == 1
  double uelec;
  CSpln* Uspl;
  Uspl = EFrc->SD;
  uelec = 0.0;
#endif

#if LISTSAME == 0
  /*** The inner loop should be the longer of the two ***/
  if (lim1-str1 > lim2-str2) {
    atomc *atmT;
    int strT, limT;
    int* ordrT;
    SWAP(atm1, atm2, atmT);
    SWAP(str1, str2, strT);
    SWAP(lim1, lim2, limT);
    SWAP(ordr1, ordr2, ordrT);
  }
#endif

  /*** The outer direct space loop ***/
  for (i = str1; i < lim1; i++) {
    const int ni = ordr1[i];
    const double atmx = atm1[ni].loc[0];
    const double atmy = atm1[ni].loc[1];
    const double atmz = atm1[ni].loc[2];

    /*** Properties of this atom ***/
    const double aiq = atm1[ni].q;
#if NEEDFORCE == 1
    aifx = 0.0;
    aify = 0.0;
    aifz = 0.0;
#endif

    /*** Inner loop for caching r2 values and successes ***/
#if LISTSAME == 1
    const int jmin = MAX(i+1, str2);
    nqpair = CalcAtomR2(lim2 - jmin, &ordr1[jmin], atm1, atmx, atmy, atmz,
			qIDbuff, qr2buff, Ecut2);
#else
    nqpair = CalcAtomR2(lim2 - str2, &ordr2[str2], atm2, atmx, atmy, atmz,
			qIDbuff, qr2buff, Ecut2);
#endif

    /*** Inner direct space loop for electrostatics ***/
    const int npairl4 = 4*(nqpair/4);
    for (j = 0; j < npairl4; j += 4) {
      const int jp1 = j+1;
      const int jp2 = j+2;
      const int jp3 = j+3;
      const int njA = qIDbuff[j];
      const int njB = qIDbuff[jp1];
      const int njC = qIDbuff[jp2];
      const int njD = qIDbuff[jp3];
      const double r2A = qr2buff[j].r2;
      const double r2B = qr2buff[jp1].r2;
      const double r2C = qr2buff[jp2].r2;
      const double r2D = qr2buff[jp3].r2;
      const int irA = r2A*uspc;
      const int irB = r2B*uspc;
      const int irC = r2C*uspc;
      const int irD = r2D*uspc;
#if NEEDFORCE == 1
#if LISTSAME == 1
      fmagA = (((Fspl[irA].A*r2A + Fspl[irA].B)*r2A + Fspl[irA].C)*r2A +
	       Fspl[irA].D)*aiq*atm1[njA].q;
      fmagB = (((Fspl[irB].A*r2B + Fspl[irB].B)*r2B + Fspl[irB].C)*r2B +
      	       Fspl[irB].D)*aiq*atm1[njB].q;
      fmagC = (((Fspl[irC].A*r2C + Fspl[irC].B)*r2C + Fspl[irC].C)*r2C +
      	       Fspl[irC].D)*aiq*atm1[njC].q;
      fmagD = (((Fspl[irD].A*r2D + Fspl[irD].B)*r2D + Fspl[irD].C)*r2D +
      	       Fspl[irD].D)*aiq*atm1[njD].q;
#else
      fmagA = (((Fspl[irA].A*r2A + Fspl[irA].B)*r2A + Fspl[irA].C)*r2A +
	       Fspl[irA].D)*aiq*atm2[njA].q;
      fmagB = (((Fspl[irB].A*r2B + Fspl[irB].B)*r2B + Fspl[irB].C)*r2B +
      	       Fspl[irB].D)*aiq*atm2[njB].q;
      fmagC = (((Fspl[irC].A*r2C + Fspl[irC].B)*r2C + Fspl[irC].C)*r2C +
      	       Fspl[irC].D)*aiq*atm2[njC].q;
      fmagD = (((Fspl[irD].A*r2D + Fspl[irD].B)*r2D + Fspl[irD].C)*r2D +
      	       Fspl[irD].D)*aiq*atm2[njD].q;
#endif
      fxA = fmagA*qr2buff[j].dx;
      fyA = fmagA*qr2buff[j].dy;
      fzA = fmagA*qr2buff[j].dz;
      fxB = fmagB*qr2buff[jp1].dx;
      fyB = fmagB*qr2buff[jp1].dy;
      fzB = fmagB*qr2buff[jp1].dz;
      fxC = fmagC*qr2buff[jp2].dx;
      fyC = fmagC*qr2buff[jp2].dy;
      fzC = fmagC*qr2buff[jp2].dz;
      fxD = fmagD*qr2buff[jp3].dx;
      fyD = fmagD*qr2buff[jp3].dy;
      fzD = fmagD*qr2buff[jp3].dz;
      aifx += fxA + fxB + fxC + fxD;
      aify += fyA + fyB + fyC + fyD;
      aifz += fzA + fzB + fzC + fzD;
#if LISTSAME == 1
      atm1[njA].frc[0] -= fxA;
      atm1[njA].frc[1] -= fyA;
      atm1[njA].frc[2] -= fzA;
      atm1[njB].frc[0] -= fxB;
      atm1[njB].frc[1] -= fyB;
      atm1[njB].frc[2] -= fzB;
      atm1[njC].frc[0] -= fxC;
      atm1[njC].frc[1] -= fyC;
      atm1[njC].frc[2] -= fzC;
      atm1[njD].frc[0] -= fxD;
      atm1[njD].frc[1] -= fyD;
      atm1[njD].frc[2] -= fzD;
#else
      atm2[njA].frc[0] -= fxA;
      atm2[njA].frc[1] -= fyA;
      atm2[njA].frc[2] -= fzA;
      atm2[njB].frc[0] -= fxB;
      atm2[njB].frc[1] -= fyB;
      atm2[njB].frc[2] -= fzB;
      atm2[njC].frc[0] -= fxC;
      atm2[njC].frc[1] -= fyC;
      atm2[njC].frc[2] -= fzC;
      atm2[njD].frc[0] -= fxD;
      atm2[njD].frc[1] -= fyD;
      atm2[njD].frc[2] -= fzD;
#endif
#endif
#if NEEDENERGY == 1
#if LISTSAME == 1
      uelec += (((Uspl[irA].A*r2A + Uspl[irA].B)*r2A + Uspl[irA].C)*r2A +
		Uspl[irA].D)*aiq*atm1[njA].q;
      uelec += (((Uspl[irB].A*r2B + Uspl[irB].B)*r2B + Uspl[irB].C)*r2B +
      		Uspl[irB].D)*aiq*atm1[njB].q;
      uelec += (((Uspl[irC].A*r2C + Uspl[irC].B)*r2C + Uspl[irC].C)*r2C +
      		Uspl[irC].D)*aiq*atm1[njC].q;
      uelec += (((Uspl[irD].A*r2D + Uspl[irD].B)*r2D + Uspl[irD].C)*r2D +
      		Uspl[irD].D)*aiq*atm1[njD].q;
#else
      uelec += (((Uspl[irA].A*r2A + Uspl[irA].B)*r2A + Uspl[irA].C)*r2A +
                Uspl[irA].D)*aiq*atm2[njA].q;
      uelec += (((Uspl[irB].A*r2B + Uspl[irB].B)*r2B + Uspl[irB].C)*r2B +
                Uspl[irB].D)*aiq*atm2[njB].q;
      uelec += (((Uspl[irC].A*r2C + Uspl[irC].B)*r2C + Uspl[irC].C)*r2C +
                Uspl[irC].D)*aiq*atm2[njC].q;
      uelec += (((Uspl[irD].A*r2D + Uspl[irD].B)*r2D + Uspl[irD].C)*r2D +
                Uspl[irD].D)*aiq*atm2[njD].q;
#endif
#if NEEDVIRIAL == 1
      sysUV->Vir[0] += fxA*qr2buff[j].dx + fxB*qr2buff[jp1].dx + 
	fxC*qr2buff[jp2].dx + fxD*qr2buff[jp3].dx;
      sysUV->Vir[1] += fxA*qr2buff[j].dy + fxB*qr2buff[jp1].dy +
	fxC*qr2buff[jp2].dy + fxD*qr2buff[jp3].dy;
      sysUV->Vir[2] += fxA*qr2buff[j].dz + fxB*qr2buff[jp1].dz +
	fxC*qr2buff[jp2].dz + fxD*qr2buff[jp3].dz;
      sysUV->Vir[3] += fyA*qr2buff[j].dx + fyB*qr2buff[jp1].dx +
	fyC*qr2buff[jp2].dx + fyD*qr2buff[jp3].dx;
      sysUV->Vir[4] += fyA*qr2buff[j].dy + fyB*qr2buff[jp1].dy +
	fyC*qr2buff[jp2].dy + fyD*qr2buff[jp3].dy;
      sysUV->Vir[5] += fyA*qr2buff[j].dz + fyB*qr2buff[jp1].dz +
	fyC*qr2buff[jp2].dz + fyD*qr2buff[jp3].dz;
      sysUV->Vir[6] += fzA*qr2buff[j].dx + fzB*qr2buff[jp1].dx +
	fzC*qr2buff[jp2].dx + fzD*qr2buff[jp3].dx;
      sysUV->Vir[7] += fzA*qr2buff[j].dy + fzB*qr2buff[jp1].dy +
	fzC*qr2buff[jp2].dy + fzD*qr2buff[jp3].dy;
      sysUV->Vir[8] += fzA*qr2buff[j].dz + fzB*qr2buff[jp1].dz +
	fzC*qr2buff[jp2].dz + fzD*qr2buff[jp3].dz;
#endif
#endif
    }
    for (j = npairl4; j < nqpair; j++) {
      const int njA = qIDbuff[j];
      const double r2A = qr2buff[j].r2;
      const int irA = r2A*uspc;
#if NEEDFORCE == 1
#if LISTSAME == 1
      fmagA = (((Fspl[irA].A*r2A + Fspl[irA].B)*r2A + Fspl[irA].C)*r2A +
	       Fspl[irA].D)*aiq*atm1[njA].q;
#else
      fmagA = (((Fspl[irA].A*r2A + Fspl[irA].B)*r2A + Fspl[irA].C)*r2A +
	       Fspl[irA].D)*aiq*atm2[njA].q;
#endif
      fxA = fmagA*qr2buff[j].dx;
      fyA = fmagA*qr2buff[j].dy;
      fzA = fmagA*qr2buff[j].dz;
      aifx += fxA;
      aify += fyA;
      aifz += fzA;
#if LISTSAME == 1
      atm1[njA].frc[0] -= fxA;
      atm1[njA].frc[1] -= fyA;
      atm1[njA].frc[2] -= fzA;
#else
      atm2[njA].frc[0] -= fxA;
      atm2[njA].frc[1] -= fyA;
      atm2[njA].frc[2] -= fzA;
#endif
#endif
#if NEEDENERGY == 1
#if LISTSAME == 1
      uelec += (((Uspl[irA].A*r2A + Uspl[irA].B)*r2A + Uspl[irA].C)*r2A +
		Uspl[irA].D)*aiq*atm1[njA].q;
#else
      uelec += (((Uspl[irA].A*r2A + Uspl[irA].B)*r2A + Uspl[irA].C)*r2A +
                Uspl[irA].D)*aiq*atm2[njA].q;
#endif
#if NEEDVIRIAL == 1
      sysUV->Vir[0] += fxA*qr2buff[j].dx;
      sysUV->Vir[1] += fxA*qr2buff[j].dy;
      sysUV->Vir[2] += fxA*qr2buff[j].dz;
      sysUV->Vir[3] += fyA*qr2buff[j].dx;
      sysUV->Vir[4] += fyA*qr2buff[j].dy;
      sysUV->Vir[5] += fyA*qr2buff[j].dz;
      sysUV->Vir[6] += fzA*qr2buff[j].dx;
      sysUV->Vir[7] += fzA*qr2buff[j].dy;
      sysUV->Vir[8] += fzA*qr2buff[j].dz;
#endif
#endif
    }
#if NEEDFORCE == 1

    /*** Accumulate the force ***/
    atm1[ni].frc[0] += aifx;
    atm1[ni].frc[1] += aify;
    atm1[ni].frc[2] += aifz;
#endif
  }
#if NEEDENERGY == 1
  sysUV->delec += uelec;
#endif
}

#undef ONLYQQINTR
