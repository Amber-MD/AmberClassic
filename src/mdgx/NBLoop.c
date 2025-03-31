#if NEEDFORCE == 1
  #if NEEDENERGY == 0
    #if VGTE == 1 && CHKRANGE == 1
      #if LISTSAME == 1
        #define COMPINTRVSSE CompSameVgtECR
        #define BOTHINTR     DoAllSameVgtECR
      #else
        #define COMPINTRVSSE CompIntrVgtECR
        #define BOTHINTR     DoAllIntrVgtECR
      #endif
    #elif VGTE == 1 && CHKRANGE == 0
      #define COMPINTRVSSE CompIntrVgtE
      #define BOTHINTR     DoAllIntrVgtE
    #elif VGTE == 0 && CHKRANGE == 1
      #if LISTSAME == 1
        #define COMPINTRVSSE CompSameVeqECR
        #define BOTHINTR     DoAllSameVeqECR
      #else
        #define COMPINTRVSSE CompIntrVeqECR
        #define BOTHINTR     DoAllIntrVeqECR
      #endif
    #else
      #define COMPINTRVSSE CompIntrVeqE
      #define BOTHINTR     DoAllIntrVeqE
    #endif
  #else
    #if NEEDVIRIAL == 0
      #if VGTE == 1 && CHKRANGE == 1
        #if LISTSAME == 1
          #define COMPINTRVSSE CompSameVgtECRnrg
          #define BOTHINTR     DoAllSameVgtECRnrg
        #else
          #define COMPINTRVSSE CompIntrVgtECRnrg
          #define BOTHINTR     DoAllIntrVgtECRnrg
        #endif
      #elif VGTE == 1 && CHKRANGE == 0
        #define COMPINTRVSSE CompIntrVgtEnrg
        #define BOTHINTR     DoAllIntrVgtEnrg
      #elif VGTE == 0 && CHKRANGE == 1
        #if LISTSAME == 1
          #define COMPINTRVSSE CompSameVeqECRnrg
          #define BOTHINTR     DoAllSameVeqECRnrg
        #else
          #define COMPINTRVSSE CompIntrVeqECRnrg
          #define BOTHINTR     DoAllIntrVeqECRnrg
        #endif
      #else
        #define COMPINTRVSSE CompIntrVeqEnrg
        #define BOTHINTR     DoAllIntrVeqEnrg
      #endif
    #else
      #if VGTE == 1 && CHKRANGE == 1
        #if LISTSAME == 1
          #define COMPINTRVSSE CompSameVgtECRnrgvir
          #define BOTHINTR     DoAllSameVgtECRnrgvir
        #else
          #define COMPINTRVSSE CompIntrVgtECRnrgvir
          #define BOTHINTR     DoAllIntrVgtECRnrgvir
        #endif
      #elif VGTE == 1 && CHKRANGE == 0
        #define COMPINTRVSSE CompIntrVgtEnrgvir
        #define BOTHINTR     DoAllIntrVgtEnrgvir
      #elif VGTE == 0 && CHKRANGE == 1
        #if LISTSAME == 1
          #define COMPINTRVSSE CompSameVeqECRnrgvir
          #define BOTHINTR     DoAllSameVeqECRnrgvir
        #else
          #define COMPINTRVSSE CompIntrVeqECRnrgvir
          #define BOTHINTR     DoAllIntrVeqECRnrgvir
        #endif
      #else
        #define COMPINTRVSSE CompIntrVeqEnrgvir
        #define BOTHINTR     DoAllIntrVeqEnrgvir
      #endif
    #endif
  #endif
#else
  #if VGTE == 1 && CHKRANGE == 1
    #if LISTSAME == 1
      #define COMPINTRVSSE CompSameVgtECRnrgx
      #define BOTHINTR     DoAllSameVgtECRnrgx
    #else
      #define COMPINTRVSSE CompIntrVgtECRnrgx
      #define BOTHINTR     DoAllIntrVgtECRnrgx
    #endif
  #elif VGTE == 1 && CHKRANGE == 0
    #define COMPINTRVSSE CompIntrVgtEnrgx
    #define BOTHINTR     DoAllIntrVgtEnrgx
  #elif VGTE == 0 && CHKRANGE == 1
    #if LISTSAME == 1
      #define COMPINTRVSSE CompSameVeqECRnrgx
      #define BOTHINTR     DoAllSameVeqECRnrgx
    #else
      #define COMPINTRVSSE CompIntrVeqECRnrgx
      #define BOTHINTR     DoAllIntrVeqECRnrgx
    #endif
  #else
    #define COMPINTRVSSE CompIntrVeqEnrgx
    #define BOTHINTR     DoAllIntrVeqEnrgx
  #endif
#endif

//-----------------------------------------------------------------------------
// BOTHINTR: the many guises of this function compute both electrostatic and
//           Lennard-Jones interactions; it assumes that all atoms it gets have
//           nonzero properties of both types.
//
// Arguments:
//   atm1:     the first list of atoms in a cell
//   atm2:     the second list of atoms in a cell--if only interactions
//             WITHIN a cell are to be computed, only the first list is given
//   ordr1:    the order of access into the atm1 array.  This stores all of the
//             atoms in atm1 numbered 0 to (number of atoms in the cell), but
//             rearranged according to the regions of the cell as they are
//             needed.
//   ordr2:    the order of access into the atm2 array--if only interactions
//             WITHIN a cell are to be computed, only ordr1 is given
//   str{1,2}: the starting points in the ordering arrays to get to the
//             regions of interest
//   lim{1,2}: the end points of the ordering arrays at which to stop (note
//             that str1 and str2 (as well as lim1 and lim2) will both be given
//             for a cell's self-interactions, as there may still be two
//             regions within the cell to compute interactions for)
//   qIDbuff:  buffer to store the indices of atoms which will make
//             electrostatic interactions
//   ljIDbuff: buffer to store the indices of atoms which will make
//             Lennard-Jones interactions
//   qr2buff:  buffer to store the distances at which successful electrostatic
//             (direct space) interactions will occur
//   ljr2buff: buffer to store the distances at which successful Lennard-Jones
//             interactions will occur
//   dcinp:    the direct space input control data
//   rcinp:    the reciprocal space input control data
//   tp:       the system topology
//   EFrc:     the electrostatic pre-tabulated force lookup table
//-----------------------------------------------------------------------------
#if NEEDENERGY == 0
#if LISTSAME == 1
static void BOTHINTR(atomc *atm1, int *ordr1, int str1, int str2,
		     int lim1, int lim2, int* qIDbuff, int* ljIDbuff,
		     rngbuff* qr2buff, rngbuff* ljr2buff, dircon *dcinp,
		     prmtop *tp, FrcTab *EFrc)
#else
static void BOTHINTR(atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		     int str1, int str2, int lim1, int lim2, int* qIDbuff,
		     int* ljIDbuff, rngbuff* qr2buff, rngbuff* ljr2buff,
		     dircon *dcinp, prmtop *tp, FrcTab *EFrc)
#endif
#else
#if NEEDFORCE == 1
#if LISTSAME == 1
static void BOTHINTR(atomc *atm1, int *ordr1, int str1, int str2, int lim1,
		     int lim2, int* qIDbuff, int* ljIDbuff, rngbuff* qr2buff,
		     rngbuff* ljr2buff, dircon *dcinp, prmtop *tp,
		     FrcTab *EFrc, Energy *sysU)
#else
static void BOTHINTR(atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		     int str1, int str2, int lim1, int lim2, int* qIDbuff,
		     int* ljIDbuff, rngbuff* qr2buff, rngbuff* ljr2buff,
		     dircon *dcinp, prmtop *tp, FrcTab *EFrc, Energy *sysU)
#endif
#else
#if LISTSAME == 1
static void BOTHINTR(atomc *atm1, int *ordr1, int str1, int str2, int lim1,
		     int lim2, int* qIDbuff, int* ljIDbuff, rngbuff* qr2buff,
		     rngbuff* ljr2buff, dircon *dcinp, prmtop *tp,
		     FrcTab *EFrc, Energy *sysU)
#else
static void BOTHINTR(atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		     int str1, int str2, int lim1, int lim2, int* qIDbuff,
		     int* ljIDbuff, rngbuff* qr2buff, rngbuff* ljr2buff,
		     dircon *dcinp, prmtop *tp, FrcTab *EFrc, Energy *sysU)
#endif
#endif
#endif
{
  int i, j, nqpair, nljpair;
#if NEEDFORCE == 1
  double aifx, aify, aifz;
  double fxA, fyA, fzA, fmagA, fxB, fyB, fzB, fmagB;
  double fxC, fyC, fzC, fmagC, fxD, fyD, fzD, fmagD;
  double *ljftmp;
  double **ljfmap;
  CSpln* Fspl;
#endif
#if NEEDENERGY == 1
  double uvdwr6, uvdwr12, uelec;
  double *ljutmp;
  double **ljumap;
  CSpln* Uspl;
#endif

  // Pre-compute a bunch of constants and unpack structures
  const double uspc = EFrc->ivdr;
  const double Vcut2 = dcinp->Vcut*dcinp->Vcut;
#if VGTE == 1
  const double Ecut2 = dcinp->Ecut*dcinp->Ecut;
#endif
#if NEEDFORCE == 1
  ljfmap = tp->LJftab.map;
  Fspl = EFrc->dSD;
#endif
#if NEEDENERGY == 1
  ljumap = tp->LJutab.map;
  Uspl = EFrc->SD;
  uelec = 0.0;
  uvdwr6 = 0.0;
  uvdwr12 = 0.0;
#endif

#if LISTSAME == 0
  // The inner loop should be the longer of the two
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

  // The outer direct space loop
  for (i = str1; i < lim1; i++) {
    const int ni = ordr1[i];
    const double atmx = atm1[ni].loc[0];
    const double atmy = atm1[ni].loc[1];
    const double atmz = atm1[ni].loc[2];

    // Properties of this atom
    const double aiq = atm1[ni].q;
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

    // Cache r2 values and successes
#if LISTSAME == 1
    const int jmin = MAX(i+1, str2);
    nljpair = CalcAtomR2(lim2-jmin, &ordr1[jmin], atm1, atmx, atmy, atmz,
			 ljIDbuff, ljr2buff, Vcut2);
#else
    nljpair = CalcAtomR2(lim2-str2, &ordr2[str2], atm2, atmx, atmy, atmz,
			 ljIDbuff, ljr2buff, Vcut2);
#endif
#if VGTE == 1
    nqpair = CullDualRange(nljpair, Ecut2, ljIDbuff, ljr2buff, qIDbuff,
			   qr2buff);
    const int nsljpair = nljpair;
#else
    nqpair = CopyAtomSuccess(nljpair, ljIDbuff, ljr2buff, qIDbuff, qr2buff);
#endif
#if CHKRANGE == 1
#if LISTSAME == 1
    nljpair = CullCloseLJ(nljpair, atm1, ljIDbuff, ljr2buff, aiid, tp);
#else
    nljpair = CullCloseLJ(nljpair, atm2, ljIDbuff, ljr2buff, aiid, tp);
#endif
#endif

    // Very likely that the electrostatics are a 
    // complete subset of the Lennard-Jones pairs
#if VGTE == 0
    if (nljpair == nqpair) {
#else
    if (nljpair == nqpair && nsljpair == nljpair) {
#endif

      // Inner direct space loop for combined
      // LJ and electrostatic interactions   
      const int npairl4 = 4*(nljpair/4);
      for (j = 0; j < npairl4; j += 4) {
	const int jp1 = j+1;
	const int jp2 = j+2;
	const int jp3 = j+3;
	const int njA = ljIDbuff[j];
	const int njB = ljIDbuff[jp1];
	const int njC = ljIDbuff[jp2];
	const int njD = ljIDbuff[jp3];
        const double r2A = ljr2buff[j].r2;
        const double r2B = ljr2buff[jp1].r2;
        const double r2C = ljr2buff[jp2].r2;
        const double r2D = ljr2buff[jp3].r2;
        const int irA = r2A*uspc;
        const int irB = r2B*uspc;
        const int irC = r2C*uspc;
        const int irD = r2D*uspc;
#if LISTSAME == 1
	const int ajljtA = atm1[njA].lj;
	const int ajljtB = atm1[njB].lj;
	const int ajljtC = atm1[njC].lj;
	const int ajljtD = atm1[njD].lj;
#else
	const int ajljtA = atm2[njA].lj;
	const int ajljtB = atm2[njB].lj;
	const int ajljtC = atm2[njC].lj;
	const int ajljtD = atm2[njD].lj;
#endif
	const double invr2A = 1.0/r2A;
	const double invr4A = invr2A*invr2A;
	const double invr2B = 1.0/r2B;
	const double invr4B = invr2B*invr2B;
	const double invr2C = 1.0/r2C;
	const double invr4C = invr2C*invr2C;
	const double invr2D = 1.0/r2D;
	const double invr4D = invr2D*invr2D;
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
	fmagA += invr4A*invr4A*(ljftmp[2*ajljtA]*invr4A*invr2A +
				ljftmp[2*ajljtA+1]);
	fmagB += invr4B*invr4B*(ljftmp[2*ajljtB]*invr4B*invr2B +
				ljftmp[2*ajljtB+1]);
	fmagC += invr4C*invr4C*(ljftmp[2*ajljtC]*invr4C*invr2C +
				ljftmp[2*ajljtC+1]);
	fmagD += invr4D*invr4D*(ljftmp[2*ajljtD]*invr4D*invr2D +
				ljftmp[2*ajljtD+1]);
	fxA = fmagA*ljr2buff[j].dx;
	fyA = fmagA*ljr2buff[j].dy;
	fzA = fmagA*ljr2buff[j].dz;
	fxB = fmagB*ljr2buff[jp1].dx;
	fyB = fmagB*ljr2buff[jp1].dy;
	fzB = fmagB*ljr2buff[jp1].dz;
	fxC = fmagC*ljr2buff[jp2].dx;
	fyC = fmagC*ljr2buff[jp2].dy;
	fzC = fmagC*ljr2buff[jp2].dz;
	fxD = fmagD*ljr2buff[jp3].dx;
	fyD = fmagD*ljr2buff[jp3].dy;
	fzD = fmagD*ljr2buff[jp3].dz;
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
	uvdwr12 += invr4A*invr2A*invr4A*invr2A*ljutmp[2*ajljtA];
	uvdwr6 += invr4A*invr2A*ljutmp[2*ajljtA+1];
	uvdwr12 += invr4B*invr2B*invr4B*invr2B*ljutmp[2*ajljtB];
	uvdwr6 += invr4B*invr2B*ljutmp[2*ajljtB+1];
	uvdwr12 += invr4C*invr2C*invr4C*invr2C*ljutmp[2*ajljtC];
	uvdwr6 += invr4C*invr2C*ljutmp[2*ajljtC+1];
	uvdwr12 += invr4D*invr2D*invr4D*invr2D*ljutmp[2*ajljtD];
	uvdwr6 += invr4D*invr2D*ljutmp[2*ajljtD+1];
#if NEEDVIRIAL == 1
	sysU->Vir[0] += fxA*ljr2buff[j].dx + fxB*ljr2buff[jp1].dx +
	  fxC*ljr2buff[jp2].dx + fxD*ljr2buff[jp3].dx;
	sysU->Vir[1] += fxA*ljr2buff[j].dy + fxB*ljr2buff[jp1].dy +
	  fxC*ljr2buff[jp2].dy + fxD*ljr2buff[jp3].dy;
	sysU->Vir[2] += fxA*ljr2buff[j].dz + fxB*ljr2buff[jp1].dz +
	  fxC*ljr2buff[jp2].dz + fxD*ljr2buff[jp3].dz;
	sysU->Vir[3] += fyA*ljr2buff[j].dx + fyB*ljr2buff[jp1].dx +
	  fyC*ljr2buff[jp2].dx + fyD*ljr2buff[jp3].dx;
	sysU->Vir[4] += fyA*ljr2buff[j].dy + fyB*ljr2buff[jp1].dy +
	  fyC*ljr2buff[jp2].dy + fyD*ljr2buff[jp3].dy;
	sysU->Vir[5] += fyA*ljr2buff[j].dz + fyB*ljr2buff[jp1].dz +
	  fyC*ljr2buff[jp2].dz + fyD*ljr2buff[jp3].dz;
	sysU->Vir[6] += fzA*ljr2buff[j].dx + fzB*ljr2buff[jp1].dx +
	  fzC*ljr2buff[jp2].dx + fzD*ljr2buff[jp3].dx;
	sysU->Vir[7] += fzA*ljr2buff[j].dy + fzB*ljr2buff[jp1].dy +
	  fzC*ljr2buff[jp2].dy + fzD*ljr2buff[jp3].dy;
	sysU->Vir[8] += fzA*ljr2buff[j].dz + fzB*ljr2buff[jp1].dz +
	  fzC*ljr2buff[jp2].dz + fzD*ljr2buff[jp3].dz;
#endif
#endif
      }
      if (j <= nljpair - 2) {
	const int jp1 = j+1;
	const int njA = ljIDbuff[j];
	const int njB = ljIDbuff[jp1];
        const double r2A = ljr2buff[j].r2;
        const double r2B = ljr2buff[jp1].r2;
        const int irA = r2A*uspc;
        const int irB = r2B*uspc;
#if LISTSAME == 1
	const int ajljtA = atm1[njA].lj;
	const int ajljtB = atm1[njB].lj;
#else
	const int ajljtA = atm2[njA].lj;
	const int ajljtB = atm2[njB].lj;
#endif
	const double invr2A = 1.0/r2A;
	const double invr4A = invr2A*invr2A;
	const double invr2B = 1.0/r2B;
	const double invr4B = invr2B*invr2B;
#if NEEDFORCE == 1
#if LISTSAME == 1
        fmagA = (((Fspl[irA].A*r2A + Fspl[irA].B)*r2A + Fspl[irA].C)*r2A +
                 Fspl[irA].D)*aiq*atm1[njA].q;
        fmagB = (((Fspl[irB].A*r2B + Fspl[irB].B)*r2B + Fspl[irB].C)*r2B +
                 Fspl[irB].D)*aiq*atm1[njB].q;
#else
        fmagA = (((Fspl[irA].A*r2A + Fspl[irA].B)*r2A + Fspl[irA].C)*r2A +
                 Fspl[irA].D)*aiq*atm2[njA].q;
        fmagB = (((Fspl[irB].A*r2B + Fspl[irB].B)*r2B + Fspl[irB].C)*r2B +
                 Fspl[irB].D)*aiq*atm2[njB].q;
#endif
	fmagA += invr4A*invr4A*(ljftmp[2*ajljtA]*invr4A*invr2A +
				ljftmp[2*ajljtA+1]);
	fmagB += invr4B*invr4B*(ljftmp[2*ajljtB]*invr4B*invr2B +
				ljftmp[2*ajljtB+1]);
	fxA = fmagA*ljr2buff[j].dx;
	fyA = fmagA*ljr2buff[j].dy;
	fzA = fmagA*ljr2buff[j].dz;
	fxB = fmagB*ljr2buff[jp1].dx;
	fyB = fmagB*ljr2buff[jp1].dy;
	fzB = fmagB*ljr2buff[jp1].dz;
	aifx += fxA + fxB;
	aify += fyA + fyB;
	aifz += fzA + fzB;
#if LISTSAME == 1
	atm1[njA].frc[0] -= fxA;
	atm1[njA].frc[1] -= fyA;
	atm1[njA].frc[2] -= fzA;
	atm1[njB].frc[0] -= fxB;
	atm1[njB].frc[1] -= fyB;
	atm1[njB].frc[2] -= fzB;
#else
	atm2[njA].frc[0] -= fxA;
	atm2[njA].frc[1] -= fyA;
	atm2[njA].frc[2] -= fzA;
	atm2[njB].frc[0] -= fxB;
	atm2[njB].frc[1] -= fyB;
	atm2[njB].frc[2] -= fzB;
#endif
#endif
#if NEEDENERGY == 1
#if LISTSAME == 1
        uelec += (((Uspl[irA].A*r2A + Uspl[irA].B)*r2A + Uspl[irA].C)*r2A +
                  Uspl[irA].D)*aiq*atm1[njA].q;
        uelec += (((Uspl[irB].A*r2B + Uspl[irB].B)*r2B + Uspl[irB].C)*r2B +
                  Uspl[irB].D)*aiq*atm1[njB].q;
#else
        uelec += (((Uspl[irA].A*r2A + Uspl[irA].B)*r2A + Uspl[irA].C)*r2A +
                  Uspl[irA].D)*aiq*atm2[njA].q;
        uelec += (((Uspl[irB].A*r2B + Uspl[irB].B)*r2B + Uspl[irB].C)*r2B +
                  Uspl[irB].D)*aiq*atm2[njB].q;
#endif
	uvdwr12 += invr4A*invr2A*invr4A*invr2A*ljutmp[2*ajljtA];
	uvdwr6 += invr4A*invr2A*ljutmp[2*ajljtA+1];
	uvdwr12 += invr4B*invr2B*invr4B*invr2B*ljutmp[2*ajljtB];
	uvdwr6 += invr4B*invr2B*ljutmp[2*ajljtB+1];
#if NEEDVIRIAL == 1
	sysU->Vir[0] += fxA*ljr2buff[j].dx + fxB*ljr2buff[jp1].dx;
	sysU->Vir[1] += fxA*ljr2buff[j].dy + fxB*ljr2buff[jp1].dy;
	sysU->Vir[2] += fxA*ljr2buff[j].dz + fxB*ljr2buff[jp1].dz;
	sysU->Vir[3] += fyA*ljr2buff[j].dx + fyB*ljr2buff[jp1].dx;
	sysU->Vir[4] += fyA*ljr2buff[j].dy + fyB*ljr2buff[jp1].dy;
	sysU->Vir[5] += fyA*ljr2buff[j].dz + fyB*ljr2buff[jp1].dz;
	sysU->Vir[6] += fzA*ljr2buff[j].dx + fzB*ljr2buff[jp1].dx;
	sysU->Vir[7] += fzA*ljr2buff[j].dy + fzB*ljr2buff[jp1].dy;
	sysU->Vir[8] += fzA*ljr2buff[j].dz + fzB*ljr2buff[jp1].dz;
#endif
#endif
	j+= 2;
      }
      if (j != nljpair) {
        const int njA = ljIDbuff[j];
        const double r2A = ljr2buff[j].r2;
        const int irA = r2A*uspc;
#if LISTSAME == 1
        const int ajljtA = atm1[njA].lj;
#else
        const int ajljtA = atm2[njA].lj;
#endif
        const double invr2A = 1.0/r2A;
        const double invr4A = invr2A*invr2A;
#if NEEDFORCE == 1
#if LISTSAME == 1
        fmagA = (((Fspl[irA].A*r2A + Fspl[irA].B)*r2A + Fspl[irA].C)*r2A +
                 Fspl[irA].D)*aiq*atm1[njA].q;
#else
        fmagA = (((Fspl[irA].A*r2A + Fspl[irA].B)*r2A + Fspl[irA].C)*r2A +
                 Fspl[irA].D)*aiq*atm2[njA].q;
#endif
        fmagA += invr4A*invr4A*(ljftmp[2*ajljtA]*invr4A*invr2A +
                                ljftmp[2*ajljtA+1]);
        fxA = fmagA*ljr2buff[j].dx;
        fyA = fmagA*ljr2buff[j].dy;
        fzA = fmagA*ljr2buff[j].dz;
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
        uvdwr12 += invr4A*invr2A*invr4A*invr2A*ljutmp[2*ajljtA];
        uvdwr6 += invr4A*invr2A*ljutmp[2*ajljtA+1];
#if NEEDVIRIAL == 1
        sysU->Vir[0] += fxA*ljr2buff[j].dx;
        sysU->Vir[1] += fxA*ljr2buff[j].dy;
        sysU->Vir[2] += fxA*ljr2buff[j].dz;
        sysU->Vir[3] += fyA*ljr2buff[j].dx;
        sysU->Vir[4] += fyA*ljr2buff[j].dy;
        sysU->Vir[5] += fyA*ljr2buff[j].dz;
        sysU->Vir[6] += fzA*ljr2buff[j].dx;
        sysU->Vir[7] += fzA*ljr2buff[j].dy;
        sysU->Vir[8] += fzA*ljr2buff[j].dz;
#endif
#endif
      }
    }
    else {

      // Inner direct space loop for LJ interactions
      const int npairl4 = (nljpair/4)*4;
      for (j = 0; j < npairl4; j += 4) {
	const int jp1 = j+1;
	const int jp2 = j+2;
	const int jp3 = j+3;
	const int njA = ljIDbuff[j];
	const int njB = ljIDbuff[jp1];
	const int njC = ljIDbuff[jp2];
	const int njD = ljIDbuff[jp3];
#if LISTSAME == 1
	const int ajljtA = atm1[njA].lj;
	const int ajljtB = atm1[njB].lj;
	const int ajljtC = atm1[njC].lj;
	const int ajljtD = atm1[njD].lj;
#else
	const int ajljtA = atm2[njA].lj;
	const int ajljtB = atm2[njB].lj;
	const int ajljtC = atm2[njC].lj;
	const int ajljtD = atm2[njD].lj;
#endif
	const double invr2A = 1.0/ljr2buff[j].r2;
	const double invr4A = invr2A*invr2A;
	const double invr2B = 1.0/ljr2buff[jp1].r2;
	const double invr4B = invr2B*invr2B;
	const double invr2C = 1.0/ljr2buff[jp2].r2;
	const double invr4C = invr2C*invr2C;
	const double invr2D = 1.0/ljr2buff[jp3].r2;
	const double invr4D = invr2D*invr2D;
#if NEEDFORCE == 1
	fmagA = invr4A*invr4A*(ljftmp[2*ajljtA]*invr4A*invr2A +
			       ljftmp[2*ajljtA+1]);
	fmagB = invr4B*invr4B*(ljftmp[2*ajljtB]*invr4B*invr2B +
			       ljftmp[2*ajljtB+1]);
	fmagC = invr4C*invr4C*(ljftmp[2*ajljtC]*invr4C*invr2C +
			       ljftmp[2*ajljtC+1]);
	fmagD = invr4D*invr4D*(ljftmp[2*ajljtD]*invr4D*invr2D +
			       ljftmp[2*ajljtD+1]);
	fxA = fmagA*ljr2buff[j].dx;
	fyA = fmagA*ljr2buff[j].dy;
	fzA = fmagA*ljr2buff[j].dz;
	fxB = fmagB*ljr2buff[jp1].dx;
	fyB = fmagB*ljr2buff[jp1].dy;
	fzB = fmagB*ljr2buff[jp1].dz;
	fxC = fmagC*ljr2buff[jp2].dx;
	fyC = fmagC*ljr2buff[jp2].dy;
	fzC = fmagC*ljr2buff[jp2].dz;
	fxD = fmagD*ljr2buff[jp3].dx;
	fyD = fmagD*ljr2buff[jp3].dy;
	fzD = fmagD*ljr2buff[jp3].dz;
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
	uvdwr12 += invr4A*invr2A*invr4A*invr2A*ljutmp[2*ajljtA];
	uvdwr6 += invr4A*invr2A*ljutmp[2*ajljtA+1];
	uvdwr12 += invr4B*invr2B*invr4B*invr2B*ljutmp[2*ajljtB];
	uvdwr6 += invr4B*invr2B*ljutmp[2*ajljtB+1];
	uvdwr12 += invr4C*invr2C*invr4C*invr2C*ljutmp[2*ajljtC];
	uvdwr6 += invr4C*invr2C*ljutmp[2*ajljtC+1];
	uvdwr12 += invr4D*invr2D*invr4D*invr2D*ljutmp[2*ajljtD];
	uvdwr6 += invr4D*invr2D*ljutmp[2*ajljtD+1];
#if NEEDVIRIAL == 1
	sysU->Vir[0] += fxA*ljr2buff[j].dx + fxB*ljr2buff[jp1].dx +
	  fxC*ljr2buff[jp2].dx + fxD*ljr2buff[jp3].dx;
	sysU->Vir[1] += fxA*ljr2buff[j].dy + fxB*ljr2buff[jp1].dy +
	  fxC*ljr2buff[jp2].dy + fxD*ljr2buff[jp3].dy;
	sysU->Vir[2] += fxA*ljr2buff[j].dz + fxB*ljr2buff[jp1].dz +
	  fxC*ljr2buff[jp2].dz + fxD*ljr2buff[jp3].dz;
	sysU->Vir[3] += fyA*ljr2buff[j].dx + fyB*ljr2buff[jp1].dx +
	  fyC*ljr2buff[jp2].dx + fyD*ljr2buff[jp3].dx;
	sysU->Vir[4] += fyA*ljr2buff[j].dy + fyB*ljr2buff[jp1].dy +
	  fyC*ljr2buff[jp2].dy + fyD*ljr2buff[jp3].dy;
	sysU->Vir[5] += fyA*ljr2buff[j].dz + fyB*ljr2buff[jp1].dz +
	  fyC*ljr2buff[jp2].dz + fyD*ljr2buff[jp3].dz;
	sysU->Vir[6] += fzA*ljr2buff[j].dx + fzB*ljr2buff[jp1].dx +
	  fzC*ljr2buff[jp2].dx + fzD*ljr2buff[jp3].dx;
	sysU->Vir[7] += fzA*ljr2buff[j].dy + fzB*ljr2buff[jp1].dy +
	  fzC*ljr2buff[jp2].dy + fzD*ljr2buff[jp3].dy;
	sysU->Vir[8] += fzA*ljr2buff[j].dz + fzB*ljr2buff[jp1].dz +
	  fzC*ljr2buff[jp2].dz + fzD*ljr2buff[jp3].dz;
#endif
#endif
      }
      if (j <= nljpair-2) {
	const int jp1 = j+1;
	const int njA = ljIDbuff[j];
	const int njB = ljIDbuff[jp1];
#if LISTSAME == 1
	const int ajljtA = atm1[njA].lj;
	const int ajljtB = atm1[njB].lj;
#else
	const int ajljtA = atm2[njA].lj;
	const int ajljtB = atm2[njB].lj;
#endif
	const double invr2A = 1.0/ljr2buff[j].r2;
	const double invr4A = invr2A*invr2A;
	const double invr2B = 1.0/ljr2buff[jp1].r2;
	const double invr4B = invr2B*invr2B;
#if NEEDFORCE == 1
	fmagA = invr4A*invr4A*(ljftmp[2*ajljtA]*invr4A*invr2A +
			       ljftmp[2*ajljtA+1]);
	fmagB = invr4B*invr4B*(ljftmp[2*ajljtB]*invr4B*invr2B +
			       ljftmp[2*ajljtB+1]);
	fxA = fmagA*ljr2buff[j].dx;
	fyA = fmagA*ljr2buff[j].dy;
	fzA = fmagA*ljr2buff[j].dz;
	fxB = fmagB*ljr2buff[jp1].dx;
	fyB = fmagB*ljr2buff[jp1].dy;
	fzB = fmagB*ljr2buff[jp1].dz;
	aifx += fxA + fxB;
	aify += fyA + fyB;
	aifz += fzA + fzB;
#if LISTSAME == 1
	atm1[njA].frc[0] -= fxA;
	atm1[njA].frc[1] -= fyA;
	atm1[njA].frc[2] -= fzA;
	atm1[njB].frc[0] -= fxB;
	atm1[njB].frc[1] -= fyB;
	atm1[njB].frc[2] -= fzB;
#else
	atm2[njA].frc[0] -= fxA;
	atm2[njA].frc[1] -= fyA;
	atm2[njA].frc[2] -= fzA;
	atm2[njB].frc[0] -= fxB;
	atm2[njB].frc[1] -= fyB;
	atm2[njB].frc[2] -= fzB;
#endif
#endif
#if NEEDENERGY == 1
	uvdwr12 += invr4A*invr2A*invr4A*invr2A*ljutmp[2*ajljtA];
	uvdwr6 += invr4A*invr2A*ljutmp[2*ajljtA+1];
	uvdwr12 += invr4B*invr2B*invr4B*invr2B*ljutmp[2*ajljtB];
	uvdwr6 += invr4B*invr2B*ljutmp[2*ajljtB+1];
#if NEEDVIRIAL == 1
	sysU->Vir[0] += fxA*ljr2buff[j].dx + fxB*ljr2buff[jp1].dx;
	sysU->Vir[1] += fxA*ljr2buff[j].dy + fxB*ljr2buff[jp1].dy;
	sysU->Vir[2] += fxA*ljr2buff[j].dz + fxB*ljr2buff[jp1].dz;
	sysU->Vir[3] += fyA*ljr2buff[j].dx + fyB*ljr2buff[jp1].dx;
	sysU->Vir[4] += fyA*ljr2buff[j].dy + fyB*ljr2buff[jp1].dy;
	sysU->Vir[5] += fyA*ljr2buff[j].dz + fyB*ljr2buff[jp1].dz;
	sysU->Vir[6] += fzA*ljr2buff[j].dx + fzB*ljr2buff[jp1].dx;
	sysU->Vir[7] += fzA*ljr2buff[j].dy + fzB*ljr2buff[jp1].dy;
	sysU->Vir[8] += fzA*ljr2buff[j].dz + fzB*ljr2buff[jp1].dz;
#endif
#endif
	j += 2;
      }
      if (j != nljpair) {
	const int njA = ljIDbuff[j];
#if LISTSAME == 1
	const int ajljtA = atm1[njA].lj;
#else
	const int ajljtA = atm2[njA].lj;
#endif
	const double invr2A = 1.0/ljr2buff[j].r2;
	const double invr4A = invr2A*invr2A;
#if NEEDFORCE == 1
	fmagA = invr4A*invr4A*(ljftmp[2*ajljtA]*invr4A*invr2A +
			       ljftmp[2*ajljtA+1]);
	fxA = fmagA*ljr2buff[j].dx;
	fyA = fmagA*ljr2buff[j].dy;
	fzA = fmagA*ljr2buff[j].dz;
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
	uvdwr12 += invr4A*invr2A*invr4A*invr2A*ljutmp[2*ajljtA];
	uvdwr6 += invr4A*invr2A*ljutmp[2*ajljtA+1];
#if NEEDVIRIAL == 1
	sysU->Vir[0] += fxA*ljr2buff[j].dx;
	sysU->Vir[1] += fxA*ljr2buff[j].dy;
	sysU->Vir[2] += fxA*ljr2buff[j].dz;
	sysU->Vir[3] += fyA*ljr2buff[j].dx;
	sysU->Vir[4] += fyA*ljr2buff[j].dy;
	sysU->Vir[5] += fyA*ljr2buff[j].dz;
	sysU->Vir[6] += fzA*ljr2buff[j].dx;
	sysU->Vir[7] += fzA*ljr2buff[j].dy;
	sysU->Vir[8] += fzA*ljr2buff[j].dz;
#endif
#endif
      }

      // Inner direct space loop for electrostatics
      const int nqpairl4 = 4*(nqpair/4);
      for (j = 0; j < nqpairl4; j += 4) {
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
        sysU->Vir[0] += fxA*ljr2buff[j].dx + fxB*ljr2buff[jp1].dx +
          fxC*ljr2buff[jp2].dx + fxD*ljr2buff[jp3].dx;
        sysU->Vir[1] += fxA*ljr2buff[j].dy + fxB*ljr2buff[jp1].dy +
          fxC*ljr2buff[jp2].dy + fxD*ljr2buff[jp3].dy;
        sysU->Vir[2] += fxA*ljr2buff[j].dz + fxB*ljr2buff[jp1].dz +
          fxC*ljr2buff[jp2].dz + fxD*ljr2buff[jp3].dz;
        sysU->Vir[3] += fyA*ljr2buff[j].dx + fyB*ljr2buff[jp1].dx +
          fyC*ljr2buff[jp2].dx + fyD*ljr2buff[jp3].dx;
        sysU->Vir[4] += fyA*ljr2buff[j].dy + fyB*ljr2buff[jp1].dy +
          fyC*ljr2buff[jp2].dy + fyD*ljr2buff[jp3].dy;
        sysU->Vir[5] += fyA*ljr2buff[j].dz + fyB*ljr2buff[jp1].dz +
          fyC*ljr2buff[jp2].dz + fyD*ljr2buff[jp3].dz;
        sysU->Vir[6] += fzA*ljr2buff[j].dx + fzB*ljr2buff[jp1].dx +
          fzC*ljr2buff[jp2].dx + fzD*ljr2buff[jp3].dx;
        sysU->Vir[7] += fzA*ljr2buff[j].dy + fzB*ljr2buff[jp1].dy +
          fzC*ljr2buff[jp2].dy + fzD*ljr2buff[jp3].dy;
        sysU->Vir[8] += fzA*ljr2buff[j].dz + fzB*ljr2buff[jp1].dz +
          fzC*ljr2buff[jp2].dz + fzD*ljr2buff[jp3].dz;
#endif
#endif
      }
      if (j <= nqpair-2) {
	const int jp1 = j+1;
	const int njA = qIDbuff[j];
	const int njB = qIDbuff[jp1];
	const double r2A = qr2buff[j].r2;
	const double r2B = qr2buff[jp1].r2;
	const int irA = r2A*uspc;
	const int irB = r2B*uspc;
#if NEEDFORCE == 1
#if LISTSAME == 1
	fmagA = (((Fspl[irA].A*r2A + Fspl[irA].B)*r2A + Fspl[irA].C)*r2A +
		 Fspl[irA].D)*aiq*atm1[njA].q;
	fmagB = (((Fspl[irB].A*r2B + Fspl[irB].B)*r2B + Fspl[irB].C)*r2B +
		 Fspl[irB].D)*aiq*atm1[njB].q;
#else
	fmagA = (((Fspl[irA].A*r2A + Fspl[irA].B)*r2A + Fspl[irA].C)*r2A +
		 Fspl[irA].D)*aiq*atm2[njA].q;
	fmagB = (((Fspl[irB].A*r2B + Fspl[irB].B)*r2B + Fspl[irB].C)*r2B +
		 Fspl[irB].D)*aiq*atm2[njB].q;
#endif
	fxA = fmagA*qr2buff[j].dx;
	fyA = fmagA*qr2buff[j].dy;
	fzA = fmagA*qr2buff[j].dz;
	fxB = fmagB*qr2buff[jp1].dx;
	fyB = fmagB*qr2buff[jp1].dy;
	fzB = fmagB*qr2buff[jp1].dz;
	aifx += fxA + fxB;
	aify += fyA + fyB;
	aifz += fzA + fzB;
#if LISTSAME == 1
	atm1[njA].frc[0] -= fxA;
	atm1[njA].frc[1] -= fyA;
	atm1[njA].frc[2] -= fzA;
	atm1[njB].frc[0] -= fxB;
	atm1[njB].frc[1] -= fyB;
	atm1[njB].frc[2] -= fzB;
#else
	atm2[njA].frc[0] -= fxA;
	atm2[njA].frc[1] -= fyA;
	atm2[njA].frc[2] -= fzA;
	atm2[njB].frc[0] -= fxB;
	atm2[njB].frc[1] -= fyB;
	atm2[njB].frc[2] -= fzB;
#endif
#endif
#if NEEDENERGY == 1
#if LISTSAME == 1
	uelec += (((Uspl[irA].A*r2A + Uspl[irA].B)*r2A + Uspl[irA].C)*r2A +
		  Uspl[irA].D)*aiq*atm1[njA].q;
	uelec += (((Uspl[irB].A*r2B + Uspl[irB].B)*r2B + Uspl[irB].C)*r2B +
		  Uspl[irB].D)*aiq*atm1[njB].q;
#else
	uelec += (((Uspl[irA].A*r2A + Uspl[irA].B)*r2A + Uspl[irA].C)*r2A +
		  Uspl[irA].D)*aiq*atm2[njA].q;
	uelec += (((Uspl[irB].A*r2B + Uspl[irB].B)*r2B + Uspl[irB].C)*r2B +
		  Uspl[irB].D)*aiq*atm2[njB].q;
#endif
#if NEEDVIRIAL == 1
	sysU->Vir[0] += fxA*qr2buff[j].dx + fxB*qr2buff[jp1].dx;
	sysU->Vir[1] += fxA*qr2buff[j].dy + fxB*qr2buff[jp1].dy;
	sysU->Vir[2] += fxA*qr2buff[j].dz + fxB*qr2buff[jp1].dz;
	sysU->Vir[3] += fyA*qr2buff[j].dx + fyB*qr2buff[jp1].dx;
	sysU->Vir[4] += fyA*qr2buff[j].dy + fyB*qr2buff[jp1].dy;
	sysU->Vir[5] += fyA*qr2buff[j].dz + fyB*qr2buff[jp1].dz;
	sysU->Vir[6] += fzA*qr2buff[j].dx + fzB*qr2buff[jp1].dx;
	sysU->Vir[7] += fzA*qr2buff[j].dy + fzB*qr2buff[jp1].dy;
	sysU->Vir[8] += fzA*qr2buff[j].dz + fzB*qr2buff[jp1].dz;
#endif
#endif
	j += 2;
      }
      if (j != nqpair) {
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
	sysU->Vir[0] += fxA*qr2buff[j].dx;
	sysU->Vir[1] += fxA*qr2buff[j].dy;
	sysU->Vir[2] += fxA*qr2buff[j].dz;
	sysU->Vir[3] += fyA*qr2buff[j].dx;
	sysU->Vir[4] += fyA*qr2buff[j].dy;
	sysU->Vir[5] += fyA*qr2buff[j].dz;
	sysU->Vir[6] += fzA*qr2buff[j].dx;
	sysU->Vir[7] += fzA*qr2buff[j].dy;
	sysU->Vir[8] += fzA*qr2buff[j].dz;
#endif
#endif
      }
    }
#if NEEDFORCE == 1

    // Accumulate the force
    atm1[ni].frc[0] += aifx;
    atm1[ni].frc[1] += aify;
    atm1[ni].frc[2] += aifz;
#endif
  }
#if NEEDENERGY == 1
  sysU->vdw12 += uvdwr12;
  sysU->vdw6 += uvdwr6;
  sysU->delec += uelec;
#endif
}

//-----------------------------------------------------------------------------
// CompIntrVSSE: compute all interactions for two lists of particles.  In the
//               event that the lists contain the same atoms, distinct versions
//               of the function will adjust the inner loop limits to avoid
//               double-counting.  In the "VgtE" version, the van-der Waals
//               cutoff is greater than the electrostatic direct-space cutoff.
//               If the two cutoffs are in fact identical, a more efficient
//               "VeqE" version is used.  The "CR" version implements range
//               checking van-der Waals interactions that might take place at
//               very close range, but this check is not needed for many cases.
//               The "nrg" extension of the function is used when a potential
//               energy computation is needed, the "nrgvir" extension is used
//               when virial calculations are required, and the "nrgx"
//               extension is used when energy, and no other property, is
//               needed.                      
//                                                                      
//               This function is really just a wrapper for three highly
//               optimized classes of functions, computing interactions
//               between particles with only Lennard-Jones properties, only
//               electrostatic properties, or both properties.     
//                                                                      
// Arguments:                                                           
//   atm1/2:  atom list 1 or 2                                          
//   ordr1/2: the order in which to access atoms from list 1            
//   natm1/2: the number of atoms in list 1                             
//   EFrc:    the electrostatic potential spline table (the force is    
//            obtained as the derivative of the spline)                 
//   dcinp:   the direct space input parameters                         
//-----------------------------------------------------------------------------
#if NEEDENERGY == 0
#if CHKRANGE == 1
#if LISTSAME == 1
void COMPINTRVSSE(cell *C, atomc *atm1, const int natm1, FrcTab *EFrc,
		  dircon *dcinp, prmtop *tp)
#else
void COMPINTRVSSE(cell *C, atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		  const int natm1, const int natm2, FrcTab *EFrc,
		  dircon *dcinp, prmtop *tp)
#endif
#else
void COMPINTRVSSE(cell *C, atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
		  const int natm1, const int natm2, FrcTab *EFrc,
		  dircon *dcinp, prmtop *tp)
#endif
#else
#if CHKRANGE == 1
#if LISTSAME == 1
void COMPINTRVSSE(cell *C, atomc *atm1, const int natm1, FrcTab *EFrc,
		  dircon *dcinp, prmtop *tp, Energy *sysUV)
#else
void COMPINTRVSSE(cell *C, atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
                  const int natm1, const int natm2, FrcTab *EFrc,
		  dircon *dcinp, prmtop *tp, Energy *sysUV)
#endif
#else
void COMPINTRVSSE(cell *C, atomc *atm1, atomc *atm2, int *ordr1, int *ordr2,
                  const int natm1, const int natm2, FrcTab *EFrc,
		  dircon *dcinp, prmtop *tp, Energy *sysUV)
#endif
#endif
{
  int i;
  int *sordr1;

  // The first order of business is to sort the lists
  // of atoms into atoms with only electrostatic     
  // properties, both electrostatic and Lennard-Jones
  // properties, and only Lennard-Jones properties.  
  int nqq1 = 0;
  int nal1 = 0;
  int nlj1 = 0;
  sordr1 = C->supordr.map[0];
  for (i = 0; i < natm1; i++) {
#if LISTSAME == 1
    if (atm1[i].lj < 0) {
      sordr1[nqq1] = i;
#else
    if (atm1[ordr1[i]].lj < 0) {
      sordr1[nqq1] = ordr1[i];
#endif
      nqq1++;
    }
  }
  for (i = 0; i < natm1; i++) {
#if LISTSAME == 1
    if (atm1[i].lj >= 0) {
      if (fabs(atm1[i].q) > 1.0e-8) {
	sordr1[nqq1+nal1] = i;
#else
    if (atm1[ordr1[i]].lj >= 0) {
      if (fabs(atm1[ordr1[i]].q) > 1.0e-8) {
	sordr1[nqq1+nal1] = ordr1[i];
#endif
	nal1++;
      }
      else {
	nlj1++;
#if LISTSAME == 1
	sordr1[natm1-nlj1] = i;
#else
	sordr1[natm1-nlj1] = ordr1[i];
#endif
      }
    }
  }
#if LISTSAME == 0
  int nqq2 = 0;
  int nal2 = 0;
  int nlj2 = 0;
  int *sordr2;
  sordr2 = C->supordr.map[1];
  for (i = 0; i < natm2; i++) {
    if (atm2[ordr2[i]].lj < 0) {
      sordr2[nqq2] = ordr2[i];
      nqq2++;
    }
  }
  for (i = 0; i < natm2; i++) {
    if (atm2[ordr2[i]].lj >= 0) {
      if (fabs(atm2[ordr2[i]].q) > 1.0e-8) {
	sordr2[nqq2+nal2] = ordr2[i];
	nal2++;
      }
      else {
	nlj2++;
	sordr2[natm2-nlj2] = ordr2[i];
      }
    }
  }
#endif

  // Now we peform a series of loops over the various 
  // classes of atoms in this cell.  First up: those  
  // with only electrostatic interactions to consider.
#if LISTSAME == 1
  if (nqq1 > 0) {
#else
  if (nqq1+nal1 > 0 && nqq2+nal2 > 0) {
#endif
#if NEEDFORCE == 1
  #if NEEDENERGY == 0
    #if LISTSAME == 1
    OnlyQQSame(atm1, sordr1, 0, 0, nqq1, nqq1+nal1, C->qIDbuff, C->qr2buff,
	       dcinp, EFrc);
    #else
    OnlyQQIntr(atm1, atm2, sordr1, sordr2, 0, 0, nqq1, nqq2+nal2, C->qIDbuff,
	       C->qr2buff, dcinp, EFrc);
    OnlyQQIntr(atm1, atm2, sordr1, sordr2, nqq1, 0, nqq1+nal1, nqq2,
	       C->qIDbuff, C->qr2buff, dcinp, EFrc);
    #endif
  #else
    #if NEEDVIRIAL == 0
      #if LISTSAME == 1
    OnlyQQSamenrg(atm1, sordr1, 0, 0, nqq1, nqq1+nal1, C->qIDbuff, C->qr2buff,
		  dcinp, EFrc, sysUV);
      #else
    OnlyQQIntrnrg(atm1, atm2, sordr1, sordr2, 0, 0, nqq1, nqq2+nal2,
		  C->qIDbuff, C->qr2buff, dcinp, EFrc, sysUV);
    OnlyQQIntrnrg(atm1, atm2, sordr1, sordr2, nqq1, 0, nqq1+nal1, nqq2,
		  C->qIDbuff, C->qr2buff, dcinp, EFrc, sysUV);
      #endif
    #else
      #if LISTSAME == 1
    OnlyQQSamenrgvir(atm1, sordr1, 0, 0, nqq1, nqq1+nal1, C->qIDbuff,
		     C->qr2buff, dcinp, EFrc, sysUV);
      #else
    OnlyQQIntrnrgvir(atm1, atm2, sordr1, sordr2, 0, 0, nqq1, nqq2+nal2,
		     C->qIDbuff, C->qr2buff, dcinp, EFrc, sysUV);
    OnlyQQIntrnrgvir(atm1, atm2, sordr1, sordr2, nqq1, 0, nqq1+nal1, nqq2,
		     C->qIDbuff, C->qr2buff, dcinp, EFrc, sysUV);
      #endif
    #endif
  #endif
#else
  #if LISTSAME == 1
    OnlyQQSamenrgx(atm1, sordr1, 0, 0, nqq1, nqq1+nal1, C->qIDbuff,
		   C->qr2buff, dcinp, EFrc, sysUV);
  #else
    OnlyQQIntrnrgx(atm1, atm2, sordr1, sordr2, 0, 0, nqq1, nqq2+nal2,
		   C->qIDbuff, C->qr2buff, dcinp, EFrc, sysUV);
    OnlyQQIntrnrgx(atm1, atm2, sordr1, sordr2, nqq1, 0, nqq1+nal1, nqq2,
		   C->qIDbuff, C->qr2buff, dcinp, EFrc, sysUV);
  #endif
#endif
  }

  // Now we loop over atoms with both electrostatic
  // and Lennard-Jones properties.                 
#if LISTSAME == 1
  if (nal1 > 0) {
#else
  if (nal1 > 0 && nal2 > 0) {
#endif
#if NEEDFORCE == 1
  #if NEEDENERGY == 0
    #if VGTE == 1 && CHKRANGE == 1
      #if LISTSAME == 1
    DoAllSameVgtECR(atm1, sordr1, nqq1, nqq1, nqq1+nal1, nqq1+nal1, C->qIDbuff,
		    C->ljIDbuff, C->qr2buff, C->ljr2buff, dcinp, tp, EFrc);
      #else
    DoAllIntrVgtECR(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
		    nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
		    C->ljr2buff, dcinp, tp, EFrc);
      #endif
    #elif VGTE == 1 && CHKRANGE == 0
    DoAllIntrVgtE(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1, nqq2+nal2,
		  C->qIDbuff, C->ljIDbuff, C->qr2buff, C->ljr2buff, dcinp,
		  tp, EFrc);
    #elif VGTE == 0 && CHKRANGE == 1
      #if LISTSAME == 1
    DoAllSameVeqECR(atm1, sordr1, nqq1, nqq1, nqq1+nal1, nqq1+nal1, C->qIDbuff,
                    C->ljIDbuff, C->qr2buff, C->ljr2buff, dcinp, tp, EFrc);
      #else
    DoAllIntrVeqECR(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
		    nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
		    C->ljr2buff, dcinp, tp, EFrc);
      #endif
    #else
    DoAllIntrVeqE(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1, nqq2+nal2,
		  C->qIDbuff, C->ljIDbuff, C->qr2buff, C->ljr2buff, dcinp,
		  tp, EFrc);
    #endif
  #else
    #if NEEDVIRIAL == 0
      #if VGTE == 1 && CHKRANGE == 1
        #if LISTSAME == 1
    DoAllSameVgtECRnrg(atm1, sordr1, nqq1, nqq1, nqq1+nal1, nqq1+nal1,
		       C->qIDbuff, C->ljIDbuff, C->qr2buff, C->ljr2buff, dcinp,
		       tp, EFrc, sysUV);
        #else
    DoAllIntrVgtECRnrg(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
		       nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
		       C->ljr2buff, dcinp, tp, EFrc, sysUV);
        #endif
      #elif VGTE == 1 && CHKRANGE == 0
    DoAllIntrVgtEnrg(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
		     nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
		     C->ljr2buff, dcinp, tp, EFrc, sysUV);
      #elif VGTE == 0 && CHKRANGE == 1
        #if LISTSAME == 1
    DoAllSameVeqECRnrg(atm1, sordr1, nqq1, nqq1, nqq1+nal1, nqq1+nal1,
                       C->qIDbuff, C->ljIDbuff, C->qr2buff, C->ljr2buff, dcinp,
                       tp, EFrc, sysUV);
        #else
    DoAllIntrVeqECRnrg(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
		       nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
		       C->ljr2buff, dcinp, tp, EFrc, sysUV);
        #endif
      #else
    DoAllIntrVeqEnrg(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
		     nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
		     C->ljr2buff, dcinp, tp, EFrc, sysUV);
      #endif
    #else
      #if VGTE == 1 && CHKRANGE == 1
        #if LISTSAME == 1
    DoAllSameVgtECRnrgvir(atm1, sordr1, nqq1, nqq1, nqq1+nal1, nqq1+nal1,
			  C->qIDbuff, C->ljIDbuff, C->qr2buff, C->ljr2buff,
			  dcinp, tp, EFrc, sysUV);
        #else
    DoAllIntrVgtECRnrgvir(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
			  nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
			  C->ljr2buff, dcinp, tp, EFrc, sysUV);
        #endif
      #elif VGTE == 1 && CHKRANGE == 0
    DoAllIntrVgtEnrgvir(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
			nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
			C->ljr2buff, dcinp, tp, EFrc, sysUV);
      #elif VGTE == 0 && CHKRANGE == 1
        #if LISTSAME == 1
    DoAllSameVeqECRnrgvir(atm1, sordr1, nqq1, nqq1, nqq1+nal1, nqq1+nal1,
                          C->qIDbuff, C->ljIDbuff, C->qr2buff, C->ljr2buff,
                          dcinp, tp, EFrc, sysUV);
        #else
    DoAllIntrVeqECRnrgvir(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
			  nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
			  C->ljr2buff, dcinp, tp, EFrc, sysUV);
        #endif
      #else
    DoAllIntrVeqEnrgvir(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
                        nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
                        C->ljr2buff, dcinp, tp, EFrc, sysUV);
      #endif
    #endif
  #endif
#else
  #if VGTE == 1 && CHKRANGE == 1
    #if LISTSAME == 1
    DoAllSameVgtECRnrgx(atm1, sordr1, nqq1, nqq1, nqq1+nal1, nqq1+nal1,
			C->qIDbuff, C->ljIDbuff, C->qr2buff, C->ljr2buff,
			dcinp, tp, EFrc, sysUV);
    #else
    DoAllIntrVgtECRnrgx(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
                        nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
                        C->ljr2buff, dcinp, tp, EFrc, sysUV);
    #endif
  #elif VGTE == 1 && CHKRANGE == 0
    DoAllIntrVgtEnrgx(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
		      nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
		      C->ljr2buff, dcinp, tp, EFrc, sysUV);
  #elif VGTE == 0 && CHKRANGE == 1
    #if LISTSAME == 1
    DoAllSameVeqECRnrgx(atm1, sordr1, nqq1, nqq1, nqq1+nal1, nqq1+nal1,
                        C->qIDbuff, C->ljIDbuff, C->qr2buff, C->ljr2buff,
                        dcinp, tp, EFrc, sysUV);
    #else
    DoAllIntrVeqECRnrgx(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
			nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
			C->ljr2buff, dcinp, tp, EFrc, sysUV);
    #endif
  #else
    DoAllIntrVeqEnrgx(atm1, atm2, sordr1, sordr2, nqq1, nqq2, nqq1+nal1,
		      nqq2+nal2, C->qIDbuff, C->ljIDbuff, C->qr2buff,
		      C->ljr2buff, dcinp, tp, EFrc, sysUV);
  #endif
#endif
  }

  // Finally, we loop over interactions that
  // involve only Lennard-Jones interactions
#if LISTSAME == 1
  if (nlj1 > 0) {
#else
  if (nlj1+nal1 > 0 && nlj2+nal2 > 0) {
#endif
#if NEEDFORCE == 1
  #if NEEDENERGY == 0
    #if CHKRANGE == 1
      #if LISTSAME == 1
    OnlyLJSameCR(atm1, sordr1, nqq1, nqq1+nal1, natm1, natm1, C->ljIDbuff,
		 C->ljr2buff, dcinp, tp);
      #else
    OnlyLJIntrCR(atm1, atm2, sordr1, sordr2, nqq1, nqq2+nal2, natm1, natm2,
		 C->ljIDbuff, C->ljr2buff, dcinp, tp);
    OnlyLJIntrCR(atm1, atm2, sordr1, sordr2, nqq1+nal1, nqq2, natm1, nqq2+nal2,
                 C->ljIDbuff, C->ljr2buff, dcinp, tp);
      #endif
    #else
      #if LISTSAME == 1
    OnlyLJSame(atm1, sordr1, nqq1, nqq1+nal1, natm1, natm1, C->ljIDbuff,
	       C->ljr2buff, dcinp, tp);
      #else
    OnlyLJIntr(atm1, atm2, sordr1, sordr2, nqq1, nqq2+nal2, natm1, natm2,
	       C->ljIDbuff, C->ljr2buff, dcinp, tp);
    OnlyLJIntr(atm1, atm2, sordr1, sordr2, nqq1+nal1, nqq2, natm1, nqq2+nal2,
	       C->ljIDbuff, C->ljr2buff, dcinp, tp);
      #endif
    #endif
  #else
    #if NEEDVIRIAL == 0
      #if CHKRANGE == 1
        #if LISTSAME == 1
    OnlyLJSameCRnrg(atm1, sordr1, nqq1, nqq1+nal1, natm1, natm1, C->ljIDbuff,
		    C->ljr2buff, dcinp, tp, sysUV);
        #else
    OnlyLJIntrCRnrg(atm1, atm2, sordr1, sordr2, nqq1, nqq2+nal2, natm1, natm2,
		    C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
    OnlyLJIntrCRnrg(atm1, atm2, sordr1, sordr2, nqq1+nal1, nqq2, natm1,
		    nqq2+nal2, C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
        #endif
      #else
        #if LISTSAME == 1
    OnlyLJSamenrg(atm1, sordr1, nqq1, nqq1+nal1, natm1, natm1, C->ljIDbuff,
		  C->ljr2buff, dcinp, tp, sysUV);
        #else
    OnlyLJIntrnrg(atm1, atm2, sordr1, sordr2, nqq1, nqq2+nal2, natm1, natm2,
		  C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
    OnlyLJIntrnrg(atm1, atm2, sordr1, sordr2, nqq1+nal1, nqq2, natm1,
		  nqq2+nal2, C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
        #endif
      #endif
    #else
      #if CHKRANGE == 1
        #if LISTSAME == 1
    OnlyLJSameCRnrgvir(atm1, sordr1, nqq1, nqq1+nal1, natm1, natm1,
		       C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
        #else
    OnlyLJIntrCRnrgvir(atm1, atm2, sordr1, sordr2, nqq1, nqq2+nal2, natm1,
		       natm2, C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
    OnlyLJIntrCRnrgvir(atm1, atm2, sordr1, sordr2, nqq1+nal1, nqq2, natm1,
		       nqq2+nal2, C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
        #endif
      #else
        #if LISTSAME == 1
    OnlyLJSamenrgvir(atm1, sordr1, nqq1, nqq1+nal1, natm1, natm1,
		     C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
        #else
    OnlyLJIntrnrgvir(atm1, atm2, sordr1, sordr2, nqq1, nqq2+nal2, natm1,
		     natm2, C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
    OnlyLJIntrnrgvir(atm1, atm2, sordr1, sordr2, nqq1+nal1, nqq2, natm1,
		     nqq2+nal2, C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
        #endif
      #endif
    #endif
  #endif
#else
  #if CHKRANGE == 1
    #if LISTSAME == 1
    OnlyLJSameCRnrgx(atm1, sordr1, nqq1, nqq1+nal1, natm1, natm1,
		     C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
    #else
    OnlyLJIntrCRnrgx(atm1, atm2, sordr1, sordr2, nqq1, nqq2+nal2, natm1,
                     natm2, C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
    OnlyLJIntrCRnrgx(atm1, atm2, sordr1, sordr2, nqq1+nal1, nqq2, natm1,
                     nqq2+nal2, C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
                    
    #endif
  #else
    #if LISTSAME == 1
    OnlyLJSamenrgx(atm1, sordr1, nqq1, nqq1+nal1, natm1, natm1,
		   C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
    #else
    OnlyLJIntrnrgx(atm1, atm2, sordr1, sordr2, nqq1, nqq2+nal2, natm1,
		   natm2, C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
    OnlyLJIntrnrgx(atm1, atm2, sordr1, sordr2, nqq1+nal1, nqq2, natm1,
		   nqq2+nal2, C->ljIDbuff, C->ljr2buff, dcinp, tp, sysUV);
    #endif
  #endif
#endif
  }
}

#undef COMPINTRVSSE
#undef BOTHINTR
