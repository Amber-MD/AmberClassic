#if NEEDFORCE == 1
  #if NEEDENERGY == 0
    #if CHKRANGE == 1
      #if LISTSAME == 1
        #define ONLYLJINTR OnlyLJSameCR
      #else
        #define ONLYLJINTR OnlyLJIntrCR
      #endif
    #else
      #define ONLYLJINTR OnlyLJIntr
    #endif
  #else
    #if NEEDVIRIAL == 0
      #if CHKRANGE == 1
        #if LISTSAME == 1
          #define ONLYLJINTR OnlyLJSameCRnrg
        #else
          #define ONLYLJINTR OnlyLJIntrCRnrg
        #endif
      #else
        #define ONLYLJINTR OnlyLJIntrnrg
      #endif
    #else
      #if CHKRANGE == 1
        #if LISTSAME == 1
          #define ONLYLJINTR OnlyLJSameCRnrgvir
        #else
          #define ONLYLJINTR OnlyLJIntrCRnrgvir
        #endif
      #else
        #define ONLYLJINTR OnlyLJIntrnrgvir
      #endif
    #endif
  #endif
#else
  #if CHKRANGE == 1
    #if LISTSAME == 1
      #define ONLYLJINTR OnlyLJSameCRnrgx
    #else
      #define ONLYLJINTR OnlyLJIntrCRnrgx
    #endif
  #else
    #define ONLYLJINTR OnlyLJIntrnrgx
  #endif
#endif

#if NEEDENERGY == 0
#if LISTSAME == 1
//-----------------------------------------------------------------------------
// ONLYLJINTR: this specialized routine will compute only Lennard-Jones
//             interactions between particles.  See the BOTHINTR function
//             description in NBLoop.c for definitions of its arguments.
//-----------------------------------------------------------------------------
static void ONLYLJINTR(atomc *atm1, int *ordr1, const int str1, const int str2,
		       const int lim1, const int lim2, int* ljIDbuff,
		       rngbuff* ljr2buff, dircon *dcinp, prmtop *tp)
#else
static void ONLYLJINTR(atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		       const int str1, const int str2, const int lim1,
		       const int lim2, int* ljIDbuff, rngbuff* ljr2buff,
		       dircon *dcinp, prmtop *tp)
#endif
#else
#if LISTSAME == 1
static void ONLYLJINTR(atomc *atm1, int *ordr1, const int str1, const int str2,
		       const int lim1, const int lim2, int* ljIDbuff,
		       rngbuff* ljr2buff, dircon *dcinp, prmtop *tp,
		       Energy *sysUV)
#else
static void ONLYLJINTR(atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		       const int str1, const int str2, const int lim1,
		       const int lim2, int* ljIDbuff, rngbuff* ljr2buff,
		       dircon *dcinp, prmtop *tp, Energy *sysUV)
#endif
#endif
{
  int i, j, npair;
#if NEEDFORCE == 1
  double aifx, aify, aifz, fx, fy, fz, fmag;
  double *ljftmp;
  double **ljfmap;
#endif
#if NEEDENERGY == 1
  double uvdwr6, uvdwr12;
  double *ljutmp;
  double **ljumap;
#endif

  /*** Pre-compute a bunch of constants and unpack structures ***/
  const double Vcut2 = dcinp->Vcut*dcinp->Vcut;
#if NEEDFORCE == 1
  ljfmap = tp->LJftab.map;
#endif
#if NEEDENERGY == 1
  ljumap = tp->LJutab.map;
  uvdwr6 = 0.0;
  uvdwr12 = 0.0;
#endif

  /*** The outer direct space loop ***/
  for (i = str1; i < lim1; i++) {
    const int ni = ordr1[i];
    const double atmx = atm1[ni].loc[0];
    const double atmy = atm1[ni].loc[1];
    const double atmz = atm1[ni].loc[2];

    /*** Properties of this atom ***/
    const int ailjt = atm1[ni].lj;
#if CHKRANGE == 1
    const int aiid = atm1[ni].id;
#endif
#if NEEDFORCE == 1
    aifx = 0.0;
    aify = 0.0;
    aifz = 0.0;
    ljftmp = ljfmap[ailjt];
#endif
#if NEEDENERGY == 1
    ljutmp = ljumap[ailjt];
#endif

    /*** Inner loop for caching r2 values and successes ***/
#if LISTSAME == 1
    const int jmin = MAX(i+1, str2);
    npair = CalcAtomR2(lim2 - jmin, &ordr1[jmin], atm1, atmx, atmy, atmz,
		       ljIDbuff, ljr2buff, Vcut2);
#else
    npair = CalcAtomR2(lim2 - str2, &ordr2[str2], atm2, atmx, atmy, atmz,
		       ljIDbuff, ljr2buff, Vcut2);
#endif
#if CHKRANGE == 1
#if LISTSAME == 1
    npair = CullCloseLJ(npair, atm1, ljIDbuff, ljr2buff, aiid, tp);
#else
    npair = CullCloseLJ(npair, atm2, ljIDbuff, ljr2buff, aiid, tp);
#endif
#endif

    /*** Inner direct space loop for LJ interactions ***/
    for (j = 0; j < npair; j++) {
      ljr2buff[j].r2 = 1.0/ljr2buff[j].r2;
    }
    for (j = 0; j < npair; j++) {
      const int nj = ljIDbuff[j];
#if LISTSAME == 1
      const int ajljt = 2*atm1[nj].lj;
#else
      const int ajljt = 2*atm2[nj].lj;
#endif
      const double invr2 = ljr2buff[j].r2;
      const double invr4 = invr2 * invr2;
#if NEEDFORCE == 1

      /*** Accumulate the force ***/
      fmag = invr4*invr4*(ljftmp[ajljt]*invr4*invr2 + ljftmp[ajljt+1]);
      fx = fmag*ljr2buff[j].dx;
      fy = fmag*ljr2buff[j].dy;
      fz = fmag*ljr2buff[j].dz;
      aifx += fx;
      aify += fy;
      aifz += fz;
#if LISTSAME == 1
      atm1[nj].frc[0] -= fx;
      atm1[nj].frc[1] -= fy;
      atm1[nj].frc[2] -= fz;
#else
      atm2[nj].frc[0] -= fx;
      atm2[nj].frc[1] -= fy;
      atm2[nj].frc[2] -= fz;
#endif
#endif
#if NEEDENERGY == 1
      uvdwr12 += invr4*invr2*invr4*invr2*ljutmp[ajljt];
      uvdwr6 += invr4*invr2*ljutmp[ajljt+1];
#if NEEDVIRIAL == 1
      sysUV->Vir[0] += fx*ljr2buff[j].dx;
      sysUV->Vir[1] += fx*ljr2buff[j].dy;
      sysUV->Vir[2] += fx*ljr2buff[j].dz;
      sysUV->Vir[3] += fy*ljr2buff[j].dx;
      sysUV->Vir[4] += fy*ljr2buff[j].dy;
      sysUV->Vir[5] += fy*ljr2buff[j].dz;
      sysUV->Vir[6] += fz*ljr2buff[j].dx;
      sysUV->Vir[7] += fz*ljr2buff[j].dy;
      sysUV->Vir[8] += fz*ljr2buff[j].dz;
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
  sysUV->vdw12 += uvdwr12;
  sysUV->vdw6 += uvdwr6;
#endif
}

#undef ONLYLJINTR
