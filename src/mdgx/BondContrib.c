//-----------------------------------------------------------------------------
// Like the nonbonded loop, there are different versions of the bond, angle,
// and dihedral routines.  The different variants can calculate forces alone,
// energies alone, forces and energies, or forces, energies, and virial
// contributions.                                  
//-----------------------------------------------------------------------------
#if NEEDFORCE == 1
  #if NEEDENERGY == 1
    #if NEEDVIRIAL == 1 
      #define ATTNPFRC AttenuatePairVir
      #define BONDCALC BondVir
      #define ANGLCALC AnglVir
      #define DIHECALC DiheVir
      #define CMAPCALC CmapVir
    #else
      #define ATTNPFRC AttenuatePairFrcNrg
      #define BONDCALC BondFrcNrg
      #define ANGLCALC AnglFrcNrg
      #define DIHECALC DiheFrcNrg
      #define CMAPCALC CmapFrcNrg
    #endif
  #else
    #define ATTNPFRC AttenuatePairFrc
    #define BONDCALC BondFrc
    #define ANGLCALC AnglFrc
    #define DIHECALC DiheFrc
    #define CMAPCALC CmapFrc
  #endif
#else
  #define ATTNPFRC AttenuatePairNrg
  #define BONDCALC BondNrg
  #define ANGLCALC AnglNrg
  #define DIHECALC DiheNrg
  #define CMAPCALC CmapNrg
#endif

//-----------------------------------------------------------------------------
// AttenuatePairForce: routine to modify the pair force between atoms in order
//                     to account for exclusion of basic van-der Waals or 
//                     electrostatic interactions.  After any necessary  
//                     corrections are added, the forces, energies, and
//                     virials are accumulated.           
//
// Arguments:                                                           
//   tp        : topology                                               
//   Cfrc      : coarse spline table for electrostatic interactions     
//   Hfrc      : fine-grained spline for electrostatic interactions     
//   atm[A,B]  : atoms A and B that are interacting with some degree of 
//               attenuation / exclusion                                
//   d[x,y,z]  : Cartesian displacement between atoms A and B           
//   fmag      : magnitude of the force (already computed before calling
//               this function)                                         
//   [a,b]frc  : pointers to forces on atoms A and B                    
//   sysUV     : system energy and virial parameters                    
//   elec14fac : electrostatic scaling factor, which could also be zero 
//               for an explicit exclusion                              
//   lj14fac   : Lennard-Jones 1-4 scaling factor                       
//-----------------------------------------------------------------------------
#if NEEDFORCE == 1
#if NEEDENERGY == 1
void ATTNPFRC(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA, int atmB,
	      double dx, double dy, double dz, double fmag, double* afrc,
	      double* bfrc, Energy *sysUV, double elec14fac, double lj14fac,
	      int qform)
#else
void ATTNPFRC(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA, int atmB,
	      double dx, double dy, double dz, double fmag, double* afrc,
	      double* bfrc, double elec14fac, double lj14fac, int qform)
#endif
#else
void ATTNPFRC(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA, int atmB,
	      double dx, double dy, double dz, Energy *sysUV, double elec14fac,
	      double lj14fac, int qform)
#endif
{
  int irc, irh;
  double r, r2, invr, invr2, qq;
#if NEEDFORCE == 1
  CSpln *CdSD, *HdSD;
#endif
#if NEEDENERGY == 1
  CSpln *CSD, *HSD;
#endif

  // (Re)compute r and related quantities
  r2 = dx*dx + dy*dy + dz*dz;
  r = sqrt(r2);  
  invr = 1.0/r;
  invr2 = invr*invr;

  // Correct the electrostatics.  In the case of standard point    
  // charges, all electrostatic interactions are counted in the    
  // nonbonded routine (direct and reciprocal space), so exclusions
  // or 1:4 screened interactions must be addressed.  However,     
  // electrostatic direct-space interactions are computed with a   
  // coarse lookup table, which is inaccurate at short range.  If  
  // an electrostatic interaction was calculated at short range, we
  // subtract off the contributions from the coarse lookup table   
  // and replace them with calculations from a more accurate fine  
  // lookup table before subtracting the analytic qq/r interaction.
  // For Gaussian charges, there is another twist: no direct space 
  // interactions were calculated, but the spline tables have been 
  // populated with erf(r)/r, not (1-erf(r))/r.                    
  qq = tp->Charges[atmA] * tp->Charges[atmB];
  if (qform == 0) {
    if (r2 < MINNB2) {
      irc = r2*Cfrc->ivdr;
#if NEEDFORCE == 1
      CdSD = Cfrc->dSD;
      fmag -= qq*(((CdSD[irc].A*r2 + CdSD[irc].B)*r2 + CdSD[irc].C)*r2 +
		  CdSD[irc].D);
#endif
#if NEEDENERGY == 1
      CSD = Cfrc->SD;
      sysUV->delec -= qq*(((CSD[irc].A*r2 + CSD[irc].B)*r2 + CSD[irc].C)*r2 +
			  CSD[irc].D);
#endif
      irh = r2*Hfrc->ivdr;
#if NEEDFORCE == 1
      HdSD = Hfrc->dSD;
      fmag += qq*(((HdSD[irh].A*r2 + HdSD[irh].B)*r2 + HdSD[irh].C)*r2 +
		  HdSD[irh].D);
#endif
#if NEEDENERGY == 1
      HSD = Hfrc->SD;
      sysUV->delec += qq*(((HSD[irh].A*r2 + HSD[irh].B)*r2 + HSD[irh].C)*r2 +
			  HSD[irh].D);
#endif
    }
#if NEEDFORCE == 1
    fmag += elec14fac*BIOQ*qq*invr*invr2;
#endif
#if NEEDENERGY == 1
    sysUV->delec -= elec14fac*BIOQ*qq*invr;
#endif
  }
  else {
    irc = r2*Cfrc->ivdr;
#if NEEDFORCE == 1
    CdSD = Cfrc->dSD;
    fmag += qq*(((CdSD[irc].A*r2 + CdSD[irc].B)*r2 + CdSD[irc].C)*r2 +
                CdSD[irc].D);
#endif
#if NEEDENERGY == 1
    CSD = Cfrc->SD;
    sysUV->delec += qq*(((CSD[irc].A*r2 + CSD[irc].B)*r2 + CSD[irc].C)*r2 +
                        CSD[irc].D);
#endif
  }

  // Correct van-der Waals eliminations that have not already been addressed,
  // or van-der Waals attenuations that have been eliminated
  int ljA, ljB;
  double invr4, invr6;
  if (r2 <= MINNB2) {
    lj14fac = lj14fac - 1.0;
  }
  ljA = tp->LJIdx[atmA];
  ljB = tp->LJIdx[atmB];
  if (ljA >= 0 && ljB >= 0) {
    invr4 = invr2*invr2;
    invr6 = invr4*invr2;
#if NEEDFORCE == 1
    fmag -= lj14fac*invr4*invr4*(tp->LJftab.map[ljA][2*ljB]*invr6 +
                                 tp->LJftab.map[ljA][2*ljB+1]);
#endif
#if NEEDENERGY == 1
    sysUV->vdw12 -= lj14fac*invr6*invr6*tp->LJutab.map[ljA][2*ljB];
    sysUV->vdw6 -= lj14fac*invr6*tp->LJutab.map[ljA][2*ljB+1];
#endif
  }

#if NEEDFORCE == 1
  // Accumulate the force
  afrc[0] += fmag*dx;
  afrc[1] += fmag*dy;
  afrc[2] += fmag*dz;
  bfrc[0] -= fmag*dx;
  bfrc[1] -= fmag*dy;
  bfrc[2] -= fmag*dz;

#if NEEDVIRIAL == 1
  // Accumulate the stress tensor
  sysUV->Vir[0] += dx*fmag*dx;
  sysUV->Vir[1] += dx*fmag*dy;
  sysUV->Vir[2] += dx*fmag*dz;
  sysUV->Vir[3] += dy*fmag*dx;
  sysUV->Vir[4] += dy*fmag*dy;
  sysUV->Vir[5] += dy*fmag*dz;
  sysUV->Vir[6] += dz*fmag*dx;
  sysUV->Vir[7] += dz*fmag*dy;
  sysUV->Vir[8] += dz*fmag*dz;
#endif
#endif
}

//-----------------------------------------------------------------------------
// BondFrc: computes the force on two atoms due to a bonded interaction.
//                                                                      
// Arguments:                                                           
//   [ab]ptr:    the coordinates for atoms A or B                       
//   [ab]frc:    the forces for atoms A or B                            
//   bcom:       the bond command                                       
//   tp:         the topology                                           
//   Cfrc:       the "coarse" electrostatic force/energy spline table   
//   Hfrc:       the "fine" electrostatic force/energy spline table     
//   sysUV:      the system energy and virial (state information)       
//-----------------------------------------------------------------------------
#if NEEDFORCE == 1
#if NEEDENERGY == 1
void BONDCALC(double *aptr, double *bptr, double *afrc, double *bfrc,
              bondcomm *bcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
              Energy *sysUV, int qform)
#else
void BONDCALC(double *aptr, double *bptr, double *afrc, double *bfrc,
              bondcomm *bcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
              int qform)
#endif
#else
void BONDCALC(double *aptr, double *bptr, bondcomm *bcom, prmtop *tp,
              FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV, int qform)
#endif
{
  double r, r2, dx, dy, dz, dl, dlpu, dlpr;
#if NEEDFORCE == 1
  double fmag;
#endif

  // Compute displacement
  dx = bptr[0] - aptr[0];
  dy = bptr[1] - aptr[1];
  dz = bptr[2] - aptr[2];

  // Accumulate the bond force and energy
  r2 = dx*dx + dy*dy + dz*dz;
  r = sqrt(r2);
  if (tp->BParam[bcom->t].l0 >= 0.0) {
    dl = tp->BParam[bcom->t].l0 - r;
  }
  else {
    dl = r + tp->BParam[bcom->t].l0;
  }
  if (r > tp->BParam[bcom->t].lpull0) {
    dlpu = tp->BParam[bcom->t].lpull0 - r;
  }
  else {
    dlpu = 0.0;
  }
  if (r < tp->BParam[bcom->t].lpress0) {
    dlpr = tp->BParam[bcom->t].lpress0 - r;
  }
  else {
    dlpr = 0.0;
  }
#if NEEDFORCE == 1  
  fmag = -2.0*((tp->BParam[bcom->t].K      * dl) +
	       (tp->BParam[bcom->t].Kpull  * dlpu) +
	       (tp->BParam[bcom->t].Kpress * dlpr))/r;
#endif
#if NEEDENERGY == 1
  sysUV->bond += (tp->BParam[bcom->t].K      * dl   * dl) +
                 (tp->BParam[bcom->t].Kpull  * dlpu * dlpu) +
                 (tp->BParam[bcom->t].Kpress * dlpr * dlpr);
  sysUV->BondUdc[3*bcom->t] += r*r;

  // WARNING: Not sure what this BondUdc is really about.  Is it correct to
  //          simply add the contributions of the extended terms like this?
  sysUV->BondUdc[3*bcom->t+1] += r*tp->BParam[bcom->t].l0 +
                                 r*tp->BParam[bcom->t].lpull0 +
                                 r*tp->BParam[bcom->t].lpress0;
  sysUV->BondUdc[3*bcom->t+2] += pow(tp->BParam[bcom->t].l0, 2.0) +
                                 pow(tp->BParam[bcom->t].lpull0, 2.0) +
                                 pow(tp->BParam[bcom->t].lpress0, 2.0);
#endif

#if NEEDFORCE == 1
#if NEEDVIRIAL == 0
#if NEEDENERGY == 1
  AttenuatePairFrcNrg(tp, Cfrc, Hfrc, bcom->a, bcom->b, dx, dy, dz, fmag,
                      afrc, bfrc, sysUV, 1.0, 1.0, qform);
#else
  AttenuatePairFrc(tp, Cfrc, Hfrc, bcom->a, bcom->b, dx, dy, dz, fmag,
                   afrc, bfrc, 1.0, 1.0, qform);
#endif
#else
  AttenuatePairVir(tp, Cfrc, Hfrc, bcom->a, bcom->b, dx, dy, dz, fmag,
                   afrc, bfrc, sysUV, 1.0, 1.0, qform);
#endif
#else
  AttenuatePairNrg(tp, Cfrc, Hfrc, bcom->a, bcom->b, dx, dy, dz, sysUV, 1.0,
                   1.0, qform);
#endif
}

//-----------------------------------------------------------------------------
// AnglFrc: computes the force on two atoms due to an angle interaction.
//                                                                      
// Arguments:                                                           
//   [abc]ptr:   the coordinates for atoms A, B, or C                   
//   [abc]frc:   the forces for atoms A, B, or C                        
//   acom:       the angle command                                      
//   tp:         the topology                                           
//   Cfrc:       the "coarse" electrostatic force/energy spline table   
//   Hfrc:       the "fine" electrostatic force/energy spline table     
//   sysUV:      the system energy and virial (state information)       
//-----------------------------------------------------------------------------
#if NEEDFORCE == 1
#if NEEDENERGY == 1
void ANGLCALC(double *aptr, double *bptr, double *cptr, double *afrc,
              double *bfrc, double *cfrc, anglcomm *acom, prmtop *tp,
              FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV, int qform)
#else
void ANGLCALC(double *aptr, double *bptr, double *cptr, double *afrc,
              double *bfrc, double *cfrc, anglcomm *acom, prmtop *tp,
              FrcTab *Cfrc, FrcTab *Hfrc, int qform)
#endif
#else
void ANGLCALC(double *aptr, double *bptr, double *cptr, anglcomm *acom,
              prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV, int qform)
#endif
{
  int i;
  double costheta, theta, mgba, mgbc;
#if NEEDFORCE == 1
  double dA, mbabc, sqba, sqbc, adf, cdf;
#endif
  double invbabc, dtheta;
  double ac[3], ba[3], bc[3];

  // Compute displacements
  for (i = 0; i < 3; i++) {
    ba[i] = aptr[i] - bptr[i];
    bc[i] = cptr[i] - bptr[i];
    ac[i] = cptr[i] - aptr[i];
  }

  // On to the angle force computation
  mgba = ba[0]*ba[0] + ba[1]*ba[1] + ba[2]*ba[2];
  mgbc = bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2];
  invbabc = 1.0/sqrt(mgba*mgbc);
  costheta = (ba[0]*bc[0] + ba[1]*bc[1] + ba[2]*bc[2]) * invbabc;
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  theta = acos(costheta);
  dtheta = theta - tp->AParam[acom->t].th0;
#if NEEDFORCE == 1
  dA = -2.0*tp->AParam[acom->t].K*dtheta / sqrt(1.0 - costheta*costheta);
  sqba = dA/mgba;
  sqbc = dA/mgbc;
  mbabc = dA * invbabc;

  // Accumulate the angle forces and stress tensor
  for (i = 0; i < 3; i++) {
    adf = bc[i]*mbabc - costheta*ba[i]*sqba;
    cdf = ba[i]*mbabc - costheta*bc[i]*sqbc;
    afrc[i] -= adf;
    bfrc[i] += adf + cdf;
    cfrc[i] -= cdf;
#if NEEDVIRIAL == 1
    sysUV->Vir[0+i] += ba[0]*adf + bc[0]*cdf;
    sysUV->Vir[3+i] += ba[1]*adf + bc[1]*cdf;
    sysUV->Vir[6+i] += ba[2]*adf + bc[2]*cdf;
#endif
  }
#endif

#if NEEDENERGY == 1
  // Accumulate the angle energy
  sysUV->angl += tp->AParam[acom->t].K*dtheta*dtheta;
#endif

  // Bail out if this angle's A:C interactions do not need to be
  // excluded (i.e. they were already excluded by some bond)    
  if (acom->excl == 0) {
    return; 
  }

#if NEEDFORCE == 1
#if NEEDVIRIAL == 0
#if NEEDENERGY == 1
  AttenuatePairFrcNrg(tp, Cfrc, Hfrc, acom->a, acom->c, ac[0], ac[1], ac[2],
                      0.0, afrc, cfrc, sysUV, 1.0, 1.0, qform);
#else
  AttenuatePairFrc(tp, Cfrc, Hfrc, acom->a, acom->c, ac[0], ac[1], ac[2], 0.0,
                   afrc, cfrc, 1.0, 1.0, qform);
#endif
#else
  AttenuatePairVir(tp, Cfrc, Hfrc, acom->a, acom->c, ac[0], ac[1], ac[2], 0.0,
                   afrc, cfrc, sysUV, 1.0, 1.0, qform);
#endif
#else
  AttenuatePairNrg(tp, Cfrc, Hfrc, acom->a, acom->c, ac[0], ac[1], ac[2],
                   sysUV, 1.0, 1.0, qform);
#endif
}

//-----------------------------------------------------------------------------
// DiheFrc: computes the force on four atoms due to a dihedral motion.  This
//          also applies to improper dihedrals.  Standard and improper
//          dihedrals are treated using the following notation:
//                                                                      
//             Standard Dihedral          Improper Dihedral             
//                                                                      
//                       D                           D                  
//                      /                           /                   
//                  B--C                        B--C                    
//                 /                                \                   
//                A                                  A                  
//                                                                      
// Arguments:                                                           
//   [abcd]ptr:   the coordinates for atoms A, B, C, or D               
//   [abcd]frc:   the forces for atoms A, B, C, or D                    
//   hcom:        the dihedral command, giving details on all Fourier series   
//                terms and a flag to tell whether this is an improper or a   
//                standard dihedral                       
//   hdef:        the array of dihedral definitions                     
//-----------------------------------------------------------------------------
#if NEEDFORCE == 1
#if NEEDENERGY == 1
void DIHECALC(double *aptr, double *bptr, double *cptr, double *dptr,
              double *afrc, double *bfrc, double *cfrc, double *dfrc,
              dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
              Energy *sysUV, int qform)
#else
void DIHECALC(double *aptr, double *bptr, double *cptr, double *dptr,
              double *afrc, double *bfrc, double *cfrc, double *dfrc,
              dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
              int qform)
#endif
#else
void DIHECALC(double *aptr, double *bptr, double *cptr, double *dptr,
              dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
              Energy *sysUV, int qform)
#endif
{
  int i, idx;
  double theta, costheta;
#if NEEDFORCE == 1
  double fr, fa, fb1, fc1, fb2, fc2, fd, isinb2, isinc2;
  double mgab, mgbc, mgcd, invab, invbc, invcd, invabc, invbcd, cosb, cosc;
#endif
  double sangle;
  double ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  dihedef *hdef;

  // Compute displacements
  for (i = 0; i < 3; i++) {
    ab[i] = bptr[i] - aptr[i];
    bc[i] = cptr[i] - bptr[i];
    cd[i] = dptr[i] - cptr[i];
  }
  CrossP(ab, bc, crabbc);
  CrossP(bc, cd, crbccd);
  costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  costheta /=
    sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
         (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  CrossP(crabbc, crbccd, scr);
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  if (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) {
    theta = acos(costheta);
  }
  else {
    theta = -acos(costheta);
  }

  // Compute the magnitude of the force and accumulate the energy
#if NEEDFORCE == 1
  fr = 0.0;
#endif
  hdef = tp->HParam;
  for (i = 0; i < hcom->nt; i++) {
    idx = hcom->t[i];
    sangle = hdef[idx].N*theta - hdef[idx].Phi;
#if NEEDFORCE == 1
    fr += hdef[idx].K * hdef[idx].N * sin(sangle);
#endif
#if NEEDENERGY == 1
    sysUV->dihe += hdef[idx].K * (1.0 + cos(sangle));
#endif
  }

#if NEEDFORCE == 1
  // Other pre-computations
  mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
  invab = 1.0/mgab;
  mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
  invbc = 1.0/mgbc;
  mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
  invcd = 1.0/mgcd;
  cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2])*invab*invbc;
  isinb2 = (cosb*cosb < 0.9999) ? 1.0/(1.0 - cosb*cosb) : 0.0;
  cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2])*invbc*invcd;
  isinc2 = (cosc*cosc < 0.9999) ? 1.0/(1.0 - cosc*cosc) : 0.0;
  isinb2 *= fr;
  isinc2 *= fr;
  invabc = invab*invbc;
  invbcd = invbc*invcd;
  for (i = 0; i < 3; i++) {
    crabbc[i] *= invabc;
    crbccd[i] *= invbcd;
  }

  // Transform the dihedral forces to cartesian coordinates
  fa = -invab * isinb2;
  fb1 = (mgbc - mgab*cosb) * invabc * isinb2;
  fb2 = cosc * invbc * isinc2;
  fc1 = (mgbc - mgcd*cosc) * invbcd * isinc2;
  fc2 = cosb * invbc * isinb2;
  fd = -invcd * isinc2;
  for (i = 0; i < 3; i++) {
    afrc[i] += crabbc[i] * fa;
    bfrc[i] += fb1 * crabbc[i] - fb2 * crbccd[i];
    cfrc[i] += -fc1 * crbccd[i] + fc2 * crabbc[i];
    dfrc[i] += -fd * crbccd[i];
  }
#endif

  // Evaluate 1-4 interactions
  if (hcom->eval14 == 1) {
#if NEEDFORCE == 1
#if NEEDVIRIAL == 0
#if NEEDENERGY == 1
  AttenuatePairFrcNrg(tp, Cfrc, Hfrc, hcom->a, hcom->d, dptr[0] - aptr[0],
                      dptr[1] - aptr[1], dptr[2] - aptr[2], 0.0, afrc, dfrc,
                      sysUV, hcom->scee, hcom->scnb, qform);
#else
  AttenuatePairFrc(tp, Cfrc, Hfrc, hcom->a, hcom->d, dptr[0] - aptr[0],
                   dptr[1] - aptr[1], dptr[2] - aptr[2], 0.0, afrc, dfrc,
                   hcom->scee, hcom->scnb, qform);
#endif
#else
  AttenuatePairVir(tp, Cfrc, Hfrc, hcom->a, hcom->d, dptr[0] - aptr[0],
                   dptr[1] - aptr[1], dptr[2] - aptr[2], 0.0, afrc, dfrc,
                   sysUV, hcom->scee, hcom->scnb, qform);
#endif
#else
  AttenuatePairNrg(tp, Cfrc, Hfrc, hcom->a, hcom->d, dptr[0] - aptr[0],
                   dptr[1] - aptr[1], dptr[2] - aptr[2], sysUV, hcom->scee,
                   hcom->scnb, qform);
#endif
  }
}

//-----------------------------------------------------------------------------
// CmapFrc: computes the interaction of five atoms due to a CMAP interaction.
//          The CMAP is treated with the following notation:
//
//                B--C
//               /    \
//              A      D--E
//
//          Here, the two dihedral angles A-B-C-D and B-C-D-E are computed per
//          standard methods (see above in DiheFrc), but the energy value and
//          the magnitudes of its derivatives are taken from the bicubic spline
//          interpolation of the CMAP surface.
//-----------------------------------------------------------------------------
#if NEEDFORCE == 1
#if NEEDENERGY == 1
void CMAPCALC(double *aptr, double *bptr, double *cptr, double *dptr,
	      double *eptr, double *afrc, double *bfrc, double *cfrc,
	      double *dfrc, double *efrc, cmapcomm *mcom, prmtop *tp,
	      Energy *sysUV, int qform)
#else
void CMAPCALC(double *aptr, double *bptr, double *cptr, double *dptr,
	      double *eptr, double *afrc, double *bfrc, double *cfrc,
	      double *dfrc, double *efrc, cmapcomm *mcom, prmtop *tp,
	      int qform)
#endif
#else
void CMAPCALC(double *aptr, double *bptr, double *cptr, double *dptr,
	      double *eptr, cmapcomm *mcom, prmtop *tp, Energy *sysUV,
	      int qform)
#endif
{
  int i, cmapX, cmapY, cmapXp, cmapYp;
  double oneOverRBC, oneOverRCD, oneOverRUBA, oneOverRUBC, oneOverRUCD;
  double oneOverRUDE, dotBABC, dotCDBC, dotBCCD, dotDECD;
  double dot, phifrac, psifrac, cx, cy, cz;
  double phi, dotphi, cosphi, sinphi, psi, dotpsi, cospsi, sinpsi;
  double nrg00, nrg01, nrg10, nrg11, phi00, phi01, phi10, phi11;
  double psi00, psi01, psi10, psi11, dpp00, dpp01, dpp10, dpp11;
  double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23;
  double a30, a31, a32, a33, dPhi, dPsi, Ecmap, rad2degCoeff;
  double ba[3], bc[3], cd[3], de[3], ubc[3], ucd[3], v[3], w[3];
  double upba[3], upbc[3], upcd[3], upde[3];
  double upabc[3], upbcd[3], upbcd1[3], upcde[3];
  dmat *csurf, *cdphi, *cdpsi, *cdpp2;
  
  // Compute displacements
  for (i = 0; i < 3; i++) {
    ba[i] = aptr[i] - bptr[i];
    bc[i] = cptr[i] - bptr[i];
    cd[i] = dptr[i] - cptr[i];
    de[i] = eptr[i] - dptr[i];
  }

  // Information on the CMAP
  csurf = &tp->MParam[mcom->idx].esrf;
  cdphi = &tp->MParam[mcom->idx].dphi;
  cdpsi = &tp->MParam[mcom->idx].dpsi;
  cdpp2 = &tp->MParam[mcom->idx].d2pp;
  rad2degCoeff = 0.5 * csurf->row / PI; 
  
  // Calculate the first dihedral, and quantities relevant for its gradient
  oneOverRBC   = 1.0 / sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
  ubc[0]       = bc[0] * oneOverRBC;
  ubc[1]       = bc[1] * oneOverRBC;
  ubc[2]       = bc[2] * oneOverRBC;
  dotBABC      = ba[0]*ubc[0] + ba[1]*ubc[1] + ba[2]*ubc[2];
  upba[0]      = ba[0] - dotBABC * ubc[0];
  upba[1]      = ba[1] - dotBABC * ubc[1];
  upba[2]      = ba[2] - dotBABC * ubc[2];
  dotBABC     *= oneOverRBC;
  oneOverRUBA  = 1.0 / sqrt(upba[0]*upba[0] + upba[1]*upba[1] +
			    upba[2]*upba[2]);
  upba[0]     *= oneOverRUBA;
  upba[1]     *= oneOverRUBA;
  upba[2]     *= oneOverRUBA;
  dotCDBC      = cd[0]*ubc[0] + cd[1]*ubc[1] + cd[2]*ubc[2];
  upcd[0]      = cd[0] - dotCDBC * ubc[0];
  upcd[1]      = cd[1] - dotCDBC * ubc[1];
  upcd[2]      = cd[2] - dotCDBC * ubc[2];
  dotCDBC     *= oneOverRBC;
  oneOverRUCD  = 1.0 / sqrt(upcd[0]*upcd[0] + upcd[1]*upcd[1] +
			    upcd[2]*upcd[2]);
  upcd[0]     *= oneOverRUCD;
  upcd[1]     *= oneOverRUCD;
  upcd[2]     *= oneOverRUCD;
  dot          = upba[0]*upcd[0] + upba[1]*upcd[1] + upba[2]*upcd[2];
  cosphi       = MIN(MAX(dot, -1.0), 1.0);
  cx           = upba[1]*upcd[2] - upba[2]*upcd[1];
  cy           = upba[2]*upcd[0] - upba[0]*upcd[2];
  cz           = upba[0]*upcd[1] - upba[1]*upcd[0];
  dot          = cx*ubc[0] + cy*ubc[1] + cz*ubc[2];
  sinphi       = MIN(MAX(dot, -1.0), 1.0);
  phi          = acos(cosphi) * (sinphi >= 0.0 ? 1.0 : -1.0);
  upabc[0]     = (ubc[2]*upba[1] - ubc[1]*upba[2]) * oneOverRUBA;
  upabc[1]     = (ubc[0]*upba[2] - ubc[2]*upba[0]) * oneOverRUBA;
  upabc[2]     = (ubc[1]*upba[0] - ubc[0]*upba[1]) * oneOverRUBA;
  upbcd[0]     = (ubc[1]*upcd[2] - ubc[2]*upcd[1]) * oneOverRUCD;
  upbcd[1]     = (ubc[2]*upcd[0] - ubc[0]*upcd[2]) * oneOverRUCD;
  upbcd[2]     = (ubc[0]*upcd[1] - ubc[1]*upcd[0]) * oneOverRUCD;

  // Calculate the second dihedral
  oneOverRCD   = 1.0 / sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
  ucd[0]       = cd[0] * oneOverRCD;
  ucd[1]       = cd[1] * oneOverRCD;
  ucd[2]       = cd[2] * oneOverRCD;
  dotBCCD      = bc[0]*ucd[0] + bc[1]*ucd[1] + bc[2]*ucd[2];
  upbc[0]      = -bc[0] + dotBCCD*ucd[0];
  upbc[1]      = -bc[1] + dotBCCD*ucd[1];
  upbc[2]      = -bc[2] + dotBCCD*ucd[2];
  dotBCCD     *= oneOverRCD;
  oneOverRUBC  = 1.0 / sqrt(upbc[0]*upbc[0] + upbc[1]*upbc[1] +
			    upbc[2]*upbc[2]);
  upbc[0]     *= oneOverRUBC;
  upbc[1]     *= oneOverRUBC;
  upbc[2]     *= oneOverRUBC;
  dotDECD      = de[0]*ucd[0] + de[1]*ucd[1] + de[2]*ucd[2];
  upde[0]      = de[0] - dotDECD*ucd[0];
  upde[1]      = de[1] - dotDECD*ucd[1];
  upde[2]      = de[2] - dotDECD*ucd[2];
  dotDECD     *= oneOverRCD;
  oneOverRUDE  = 1.0 / sqrt(upde[0]*upde[0] + upde[1]*upde[1] +
			    upde[2]*upde[2]);
  upde[0]     *= oneOverRUDE;
  upde[1]     *= oneOverRUDE;
  upde[2]     *= oneOverRUDE;
  dot          = upbc[0]*upde[0] + upbc[1]*upde[1] + upbc[2]*upde[2];
  cospsi       = MIN(MAX(dot, -1.0), 1.0);
  cx           = upbc[1]*upde[2] - upbc[2]*upde[1];
  cy           = upbc[2]*upde[0] - upbc[0]*upde[2];
  cz           = upbc[0]*upde[1] - upbc[1]*upde[0];
  dot          = cx*ucd[0] + cy*ucd[1] + cz*ucd[2];
  sinpsi       = MIN(MAX(dot, -1.0), 1.0);
  psi          = acos(cospsi) * (sinpsi >= 0.0 ? 1.0 : -1.0);

  // Compute the cross terms
  upbcd1[0]    = (upbc[1]*ucd[2] - upbc[2]*ucd[1]) * oneOverRUBC;
  upbcd1[1]    = (upbc[2]*ucd[0] - upbc[0]*ucd[2]) * oneOverRUBC;
  upbcd1[2]    = (upbc[0]*ucd[1] - upbc[1]*ucd[0]) * oneOverRUBC;
  upcde[0]     = (ucd[1]*upde[2] - ucd[2]*upde[1]) * oneOverRUDE;
  upcde[1]     = (ucd[2]*upde[0] - ucd[0]*upde[2]) * oneOverRUDE;
  upcde[2]     = (ucd[0]*upde[1] - ucd[1]*upde[0]) * oneOverRUDE;
  phi         += PI;
  psi         += PI;
  cmapX        = phi * (csurf->row / (2.0 * PI));
  cmapY        = psi * (csurf->row / (2.0 * PI));
  phifrac      = (phi - cmapX*(2.0 * PI / csurf->row)) * csurf->row / (2.0*PI);
  psifrac      = (psi - cmapY*(2.0 * PI / csurf->row)) * csurf->row / (2.0*PI);

  // Pluck energy, dphi, dpsi, and dphi/dpsi grid values from the CMAP
  cmapXp = cmapX + 1;
  cmapYp = cmapY + 1;
  cmapXp -= (cmapXp == csurf->row) * csurf->row;
  cmapYp -= (cmapYp == csurf->col) * csurf->col;
  nrg00 = csurf->map[cmapX][cmapY];
  nrg01 = csurf->map[cmapX][cmapYp];
  nrg10 = csurf->map[cmapXp][cmapY];
  nrg11 = csurf->map[cmapXp][cmapYp];
  phi00 = cdphi->map[cmapX][cmapY];
  phi01 = cdphi->map[cmapX][cmapYp];
  phi10 = cdphi->map[cmapXp][cmapY];
  phi11 = cdphi->map[cmapXp][cmapYp];
  psi00 = cdpsi->map[cmapX][cmapY];
  psi01 = cdpsi->map[cmapX][cmapYp];
  psi10 = cdpsi->map[cmapXp][cmapY];
  psi11 = cdpsi->map[cmapXp][cmapYp];
  dpp00 = cdpp2->map[cmapX][cmapY];
  dpp01 = cdpp2->map[cmapX][cmapYp];
  dpp10 = cdpp2->map[cmapXp][cmapY];
  dpp11 = cdpp2->map[cmapXp][cmapYp];
  
  // Apply the energy values
  a00        =                   nrg00;
  a10        =                   phi00;
  a20        = (-3.0 * nrg00) + ( 3.0 * nrg10) -
               ( 2.0 * phi00) -                   phi10;
  a30        = ( 2.0 * nrg00) - ( 2.0 * nrg10) +
                       phi00  +                   phi10;
  a01        =                   psi00;
  a11        =                   dpp00;
  a21        = (-3.0 * psi00) + ( 3.0 * psi10) -
               ( 2.0 * dpp00) -                   dpp10;
  a31        = ( 2.0 * psi00) - ( 2.0 * psi10) +
                       dpp00  +                   dpp10;
  a02        = (-3.0 * nrg00) + ( 3.0 * nrg01) -
               ( 2.0 * psi00) -                   psi01;
  a12        = (-3.0 * phi00) + ( 3.0 * phi01) -
               ( 2.0 * dpp00) -                   dpp01;
  a22        = ( 9.0 * (nrg00 - nrg10 - nrg01 + nrg11) ) +
               ( 6.0 * phi00) + ( 3.0 * phi10) -
               ( 6.0 * phi01) - ( 3.0 * phi11) +
               ( 6.0 * psi00) - ( 6.0 * psi10) +
               ( 3.0 * psi01) - ( 3.0 * psi11) +
               ( 4.0 * dpp00) + ( 2.0 * dpp10) +
               ( 2.0 * dpp01) +                   dpp11;
  a32        = (-6.0 * (nrg00 - nrg10 - nrg01 + nrg11) ) +
               (-3.0 * (phi00 + phi10 - phi01 - phi11) ) +
               (-4.0 * psi00) + ( 4.0 * psi10) -
                       dpp01  -                   dpp11  +
               (-2.0 * (dpp00 + dpp10 + psi01 - psi11) );
  a03        = ( 2.0 * nrg00) - ( 2.0 * nrg01) +
                       psi00  +                   psi01;
  a13        = ( 2.0 * phi00) - ( 2.0 * phi01) +
                       dpp00  +                   dpp01;
  a23        = (-6.0 * (nrg00 - nrg10 - nrg01 + nrg11) ) +
               (-2.0 * (dpp00 + phi10 + dpp01 - phi11) ) +
               (-3.0 * (psi00 - psi10 + psi01 - psi11) ) +
               (-4.0 * phi00) - dpp10 +
               ( 4.0 * phi01) - dpp11;
  a33        = ( 4.0 * (nrg00 - nrg10 - nrg01 + nrg11) ) +
               ( 2.0 * (phi00 + phi10 - phi01 - phi11) ) +
               ( 2.0 * (psi00 - psi10 + psi01 - psi11) ) +
                        dpp00 +                   dpp10  +
                        dpp01 +                   dpp11;
  dPhi       = (3.0*a33*phifrac + 2.0*a23)*phifrac + a13;
  dPhi       = (3.0*a32*phifrac + 2.0*a22)*phifrac + a12 + dPhi*psifrac;
  dPhi       = (3.0*a31*phifrac + 2.0*a21)*phifrac + a11 + dPhi*psifrac;
  dPhi       = (3.0*a30*phifrac + 2.0*a20)*phifrac + a10 + dPhi*psifrac;
  dPsi       = (3.0*a33*psifrac + 2.0*a32)*psifrac + a31;
  dPsi       = (3.0*a23*psifrac + 2.0*a22)*psifrac + a21 + dPsi*phifrac;
  dPsi       = (3.0*a13*psifrac + 2.0*a12)*psifrac + a11 + dPsi*phifrac;
  dPsi       = (3.0*a03*psifrac + 2.0*a02)*psifrac + a01 + dPsi*phifrac;
  dPhi      *= rad2degCoeff;
  dPsi      *= rad2degCoeff;
  upabc[0]  *= dPhi;
  upabc[1]  *= dPhi;
  upabc[2]  *= dPhi;
  upbcd[0]  *= dPhi;
  upbcd[1]  *= dPhi;
  upbcd[2]  *= dPhi;
  upbcd1[0] *= dPsi;
  upbcd1[1] *= dPsi;
  upbcd1[2] *= dPsi;
  upcde[0]  *= dPsi;
  upcde[1]  *= dPsi;
  upcde[2]  *= dPsi;

  // Calculate gradients and accumulate forces
#if NEEDFORCE == 1
  for (i = 0; i < 3; i++) {
    v[i] =  dotBABC*upabc[i]  + dotCDBC*upbcd[i];
    w[i] = -dotBCCD*upbcd1[i] + dotDECD*upcde[i];
    afrc[i] -= upabc[i];
    bfrc[i] += upabc[i] - v[i] - upbcd1[i];
    cfrc[i] += v[i] + upbcd[i] - w[i] + upbcd1[i];
    dfrc[i] += upcde[i] + w[i] - upbcd[i];
    efrc[i] -= upcde[i];
  }
#endif

  // Calculate energy
#if NEEDENERGY == 1
  Ecmap = ((a33*psifrac + a32)*psifrac + a31)*psifrac + a30;
  Ecmap = ((a23*psifrac + a22)*psifrac + a21)*psifrac + a20 + phifrac*Ecmap;
  Ecmap = ((a13*psifrac + a12)*psifrac + a11)*psifrac + a10 + phifrac*Ecmap;
  Ecmap = ((a03*psifrac + a02)*psifrac + a01)*psifrac + a00 + phifrac*Ecmap;
  sysUV->cmap += Ecmap;
#endif
}

#undef BONDCALC
#undef ANGLCALC
#undef DIHECALC
#undef CMAPCALC
#undef ATTNPFRC
