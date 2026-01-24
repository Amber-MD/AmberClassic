//-----------------------------------------------------------------------------
// VectorBond: compute the energy and forces due to a bonds in one system for
//             a series of configurations.
//
// Arguments:
//   tp:       system topology
//   crdpop:   matrix of coordinates, giving X1, Y1, Z1, X2, ..., Zn on
//             successive rows
//   frcpop:   matrix of forces, being accumulated before this function is
//             called and accepting contributions from this function
//   skipintr: flag to have the function return after computing just the
//             bond lengths (set to 0 to have the energies computed)
//   sysUV:    vector of energy accumulators for each conformation
//   atm[ab]:  the atoms to which this bond a-b applies
//   r:        array to store the computed values of the bond length
//   Keq:      stiffness of the harmonic bond (or restraint)
//   Leq:      equilibrium length of the harmonic bond (a vector of numbers so
//             that different configurations can have different ideal lengths
//             in the case that this function is executing a set of restraints
//             instead of a bond in the actual topology)
//   Ftop:     the highest magnitude of the force that can be applied
//   compfrc:  flag to have forces computed if set to 1 (otherwise only the
//             energy will be computed)
//   nactive:  the number of configurations actively being manipulated
//-----------------------------------------------------------------------------
#if DOENERGY == 1
static void VectorBond(dmat *crdpop, int skipintr, Energy* sysUV, int atma,
		       int atmb, double* r, double Keq, double Kpulleq,
		       double Kpresseq, double* Leq, double Lpulleq,
		       double Lpresseq, int nactive)
#elif DOFORCE == 1
static void VectorBondF(dmat *crdpop, dmat *frcpop, Energy* sysUV, int atma,
			int atmb, double* r, double Keq, double Kpulleq,
			double Kpresseq, double* Leq, double Lpulleq,
			double Lpresseq, int nactive)
#elif DORESTRAINT == 1
static void VectorBondR(dmat *crdpop, dmat *frcpop, Energy* sysUV, int atma,
			int atmb, double* r, double Keq, double* Leq,
			double fbhw, double quadwin, int nactive)
#elif DORESTE == 1
static void VectorBondRU(dmat *crdpop, Energy* sysUV, int atma, int atmb,
			 double* r, double Keq, double* Leq, double fbhw,
			 double quadwin, int nactive)
#endif
{
  int i;
  double dx, dy, dz, r2, dl, fmag;
#if DOFORCE == 1 || DOENERGY == 1
  double dlpr, dlpu;
#endif
  double *xcrda, *xcrdb, *ycrda, *ycrdb, *zcrda, *zcrdb;
  double *xfrca, *xfrcb, *yfrca, *yfrcb, *zfrca, *zfrcb;

  // Loop over all active configurations
  xcrda = crdpop->map[3*atma];
  xcrdb = crdpop->map[3*atmb];
  ycrda = crdpop->map[3*atma+1];
  ycrdb = crdpop->map[3*atmb+1];
  zcrda = crdpop->map[3*atma+2];
  zcrdb = crdpop->map[3*atmb+2];

  // Set pointers for forces if they are needed
#if DOFORCE == 1 || DORESTRAINT == 1
  xfrca = frcpop->map[3*atma];
  xfrcb = frcpop->map[3*atmb];
  yfrca = frcpop->map[3*atma+1];
  yfrcb = frcpop->map[3*atmb+1];
  zfrca = frcpop->map[3*atma+2];
  zfrcb = frcpop->map[3*atmb+2];
#endif

  // For restraints, pre-compute the slope of the line that
  // the restraint will become after the force flatlines
#if DORESTRAINT == 1 || DORESTE == 1
  const double Ftop = 2.0*Keq*quadwin;
  const double quadU = Keq*quadwin*quadwin;
#endif

  // Loop over all active conformations
  for (i = 0; i < nactive; i++) {
    dx = xcrdb[i]-xcrda[i];
    dy = ycrdb[i]-ycrda[i];
    dz = zcrdb[i]-zcrda[i];
    r2 = dx*dx + dy*dy + dz*dz;
    r[i] = sqrt(r2);
#if DOENERGY == 1
    if (skipintr == 1) {
      continue;
    }
#endif
    dl = (Leq[i] >= 0.0) ? Leq[i] - r[i] : r[i] + Leq[i];
    
    // For flat-bottom harmonic restraints, we consider
    // whether dl is within the flat bottom region.
    // Always adjust dl so that it reflects the deviation
    // in the region where the potential is nonzero.  If
    // this is, instead,  a standard bond, then there are
    // restraint-like components to consider there as well.
#if DORESTRAINT == 1 || DORESTE == 1
    if (fabs(dl) < fbhw) {
      dl = 0.0;
    }
    else if (dl < 0.0) {
      dl += fbhw;
    }
    else {
      dl -= fbhw;
    }
#else
    dlpu = (r[i] > Lpulleq)  ? Lpulleq  - r[i] : 0.0;
    dlpr = (r[i] < Lpresseq) ? Lpresseq - r[i] : 0.0;
#endif

    // Simple harmonic potentials imply simple energies,
    // flat-bottom potentials have piecewise energies.
#if DOENERGY == 1 || DOFORCE == 1    
    sysUV[i].bond += Keq*dl*dl + Kpulleq*dlpu*dlpu + Kpresseq*dlpr*dlpr;
#endif
#if DOFORCE == 1
    fmag = -2.0 * (Keq*dl + Kpulleq*dlpu + Kpresseq*dlpr) / r[i];
#endif
#if DORESTRAINT == 1 || DORESTE == 1
    // Restraint energies get tallied in the reciprocal space
    // electrostatic energy, as this quantity is otherwise unused.
    if (fabs(dl) < quadwin) {
      sysUV[i].relec += Keq*dl*dl;
#if DORESTRAINT == 1
      fmag = -2.0*Keq*dl/r[i];
#endif
    }
    else {
      sysUV[i].relec += quadU + Ftop*(fabs(dl) - quadwin);
#if DORESTRAINT == 1
      if (dl < 0.0) {
	fmag = Ftop;
      }
      else {
	fmag = -Ftop;
      }
#endif
    }
#endif
#if DOFORCE == 1 || DORESTRAINT == 1
    xfrca[i] += fmag*dx;
    xfrcb[i] -= fmag*dx;
    yfrca[i] += fmag*dy;
    yfrcb[i] -= fmag*dy;
    zfrca[i] += fmag*dz;
    zfrcb[i] -= fmag*dz;
#endif
  }
}

//-----------------------------------------------------------------------------
// VectorAngle: analogous to VectorBond, computes the angle interactions of a
//              system with many configurations.
//
// Arguments:
//   tp:       system topology
//   crdpop:   matrix of coordinates, giving X1, Y1, Z1, X2, ..., Zn on
//             successive rows
//   frcpop:   matrix of forces, being accumulated before this function is
//             called and accepting contributions from this function
//   skipintr: flag to have the function return after computing just the
//             angles (set to 0 to have the energies computed)
//   sysUV:    vector of energy accumulators for each conformation
//   atm[abc]: the atoms to which this angle a-b-c applies
//   theta:    array to store the computed values of the angle
//   Keq:      stiffness of the harmonic angle (or restraint)
//   THeq:     equilibrium value of the harmonic angle (a vector of numbers so
//             that different configurations can have different ideal lengths
//             in the case that this function is executing a set of restraints
//             instead of a bond in the actual topology)
//   Ftop:     the highest magnitude of the force that can be applied
//   compfrc:  flag to have forces computed if set to 1 (otherwise only the
//             energy will be computed)
//   nactive:  the number of configurations actively being manipulated
//-----------------------------------------------------------------------------
#if DOENERGY == 1
static void VectorAngle(dmat *crdpop, int skipintr, Energy* sysUV, int atma,
			int atmb, int atmc, double* theta, double Keq,
			double* THeq, int nactive)
#elif DOFORCE == 1
static void VectorAngleF(dmat *crdpop, dmat *frcpop, Energy* sysUV, int atma,
			 int atmb, int atmc, double* theta, double Keq,
			 double* THeq, int nactive)
#elif DORESTRAINT == 1
static void VectorAngleR(dmat *crdpop, dmat *frcpop, Energy* sysUV, int atma,
			 int atmb, int atmc, double* theta, double Keq,
			 double* THeq, double fbhw, double quadwin,
			 int nactive)
#elif DORESTE == 1
static void VectorAngleRU(dmat *crdpop, Energy* sysUV, int atma, int atmb,
			  int atmc, double* theta, double Keq, double* THeq,
			  double fbhw, double quadwin, int nactive)
#endif
{
  int i;
  double costheta, dtheta, Ftop2, normfac;
  double mgba, mgbc, dA, mbabc, sqba, sqbc, invbabc;
  double acx, acy, acz, bax, bay, baz, bcx, bcy, bcz;
  double adfx, adfy, adfz, cdfx, cdfy, cdfz, adf2, cdf2, bdf2;
  double *xcrda, *xcrdb, *xcrdc, *ycrda, *ycrdb, *ycrdc;
  double *zcrda, *zcrdb, *zcrdc;
  double *xfrca, *xfrcb, *xfrcc, *yfrca, *yfrcb, *yfrcc;
  double *zfrca, *zfrcb, *zfrcc;
  anglcomm *thisac;

  // Compute the maximum force magnitude for convenience
#if DORESTRAINT == 1 || DORESTE == 1
  const double Ftop = 2.0*Keq*quadwin;
  const double quadU = Keq*quadwin*quadwin;
#endif

  // Loop over all active configurations
  xcrda = crdpop->map[3*atma];
  xcrdb = crdpop->map[3*atmb];
  xcrdc = crdpop->map[3*atmc];
  ycrda = crdpop->map[3*atma+1];
  ycrdb = crdpop->map[3*atmb+1];
  ycrdc = crdpop->map[3*atmc+1];
  zcrda = crdpop->map[3*atma+2];
  zcrdb = crdpop->map[3*atmb+2];
  zcrdc = crdpop->map[3*atmc+2];

  // Set pointers for forces if needed
#if DOFORCE == 1 || DORESTRAINT == 1
  xfrca = frcpop->map[3*atma];
  xfrcb = frcpop->map[3*atmb];
  xfrcc = frcpop->map[3*atmc];
  yfrca = frcpop->map[3*atma+1];
  yfrcb = frcpop->map[3*atmb+1];
  yfrcc = frcpop->map[3*atmc+1];
  zfrca = frcpop->map[3*atma+2];
  zfrcb = frcpop->map[3*atmb+2];
  zfrcc = frcpop->map[3*atmc+2];
#endif

  // Loop over all active conformations
  for (i = 0; i < nactive; i++) {
    bax = xcrda[i] - xcrdb[i];
    bay = ycrda[i] - ycrdb[i];
    baz = zcrda[i] - zcrdb[i];
    bcx = xcrdc[i] - xcrdb[i];
    bcy = ycrdc[i] - ycrdb[i];
    bcz = zcrdc[i] - zcrdb[i];
    acx = xcrdc[i] - xcrda[i];
    acy = ycrdc[i] - ycrda[i];
    acz = zcrdc[i] - zcrda[i];
    mgba = bax*bax + bay*bay + baz*baz;
    mgbc = bcx*bcx + bcy*bcy + bcz*bcz;
    invbabc = 1.0/sqrt(mgba*mgbc);
    costheta = (bax*bcx + bay*bcy + baz*bcz) * invbabc;
    costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
    theta[i] = acos(costheta);
#if DOENERGY == 1
    if (skipintr == 1) {
      continue;
    }
#endif
    dtheta = theta[i] - THeq[i];

    // For flat-bottom harmonic restraints, we consider
    // whether dtheta is within the flat bottom region.
    // Always adjust dtheta so that it reflects the deviation
    // in the region where the potential is nonzero.
#if DORESTRAINT == 1 || DORESTE == 1
    if (fabs(dtheta) < fbhw) {
      dtheta = 0.0;
    }
    else if (dtheta < 0.0) {
      dtheta += fbhw;
    }
    else {
      dtheta -= fbhw;
    }
#endif

    // Standard harmonic potentials are simple to compute,
    // but the flat-bottom well demands a bit more care.
#if DOENERGY == 1 || DOFORCE == 1
    sysUV[i].angl += Keq*dtheta*dtheta;
#endif
#if DOFORCE == 1
    dA = -2.0*Keq*dtheta / sqrt(1.0 - costheta*costheta);
#endif
#if DORESTRAINT == 1 || DORESTE == 1
    // Restraint energies get tallied in the reciprocal space
    // electrostatic energy, as this quantity is otherwise unused.
    if (fabs(dtheta) < quadwin) {
      sysUV[i].relec += Keq*dtheta*dtheta;
#if DORESTRAINT == 1
      dA = -2.0*Keq*dtheta / sqrt(1.0 - costheta*costheta);
#endif
    }
    else {
      sysUV[i].relec += quadU + Ftop*(fabs(dtheta) - quadwin);
#if DORESTRAINT == 1
      if (dtheta < 0.0) {
        dA = Ftop / sqrt(1.0 - costheta*costheta);
      }
      else {
        dA = -Ftop / sqrt(1.0 - costheta*costheta);
      }
#endif
    }
#endif
#if DOFORCE == 1 || DORESTRAINT == 1
    sqba = dA/mgba;
    sqbc = dA/mgbc;
    mbabc = dA * invbabc;
    adfx = bcx*mbabc - costheta*bax*sqba;
    cdfx = bax*mbabc - costheta*bcx*sqbc;
    adfy = bcy*mbabc - costheta*bay*sqba;
    cdfy = bay*mbabc - costheta*bcy*sqbc;
    adfz = bcz*mbabc - costheta*baz*sqba;
    cdfz = baz*mbabc - costheta*bcz*sqbc;
    xfrca[i] -= adfx;
    xfrcb[i] += adfx + cdfx;
    xfrcc[i] -= cdfx;
    yfrca[i] -= adfy;
    yfrcb[i] += adfy + cdfy;
    yfrcc[i] -= cdfy;
    zfrca[i] -= adfz;
    zfrcb[i] += adfz + cdfz;
    zfrcc[i] -= cdfz;
#endif
  }
}

//-----------------------------------------------------------------------------
// VectorTorsion: this computes the torsion angles of a system with many
//                configurations and returns them if requested.
//
// Arguments:
//   crdpop:    matrix of coordinates, giving X1, Y1, Z1, X2, ..., Zn on
//              successive rows
//   frcpop:    matrix of forces, being accumulated before this function is
//              called and accepting contributions from this function
//   atm[abcd]: the atoms to which this torsion a-b-c-d applies
//   theta:     array to store the computed values of the dihedral angle
//   skipintr:  flag to have the function return after computing just the
//              dihedal angles (set to 0 to have the energies computed)
//   nterm:     the number of Fourier terms in this dihedral
//   Kamp:      amplitudes of all Fourier series terms
//   Nperiod:   periodicities of all Fourier series terms
//   Phase:     phase angles of all conformations for all Fourier series terms
//              (each row is the set of phase angles for one configuration
//              and all terms)
//   compfrc:   flag to have forces computed if set to 1 (otherwise only the
//              energy will be computed)
//   nactive:   the number of configurations actively being manipulated
//-----------------------------------------------------------------------------
#if DOENERGY == 1
static void VectorTorsion(dmat *crdpop, int atma, int atmb, int atmc, int atmd,
			  double* theta, int skipintr, Energy* sysUV,
			  int nterm, double* Kamp, double* Nperiod,
			  dmat *Phase, int nactive)
#elif DOFORCE == 1
static void VectorTorsionF(dmat *crdpop, dmat *frcpop, int atma, int atmb,
			   int atmc, int atmd, double* theta, Energy* sysUV,
			   int nterm, double* Kamp, double* Nperiod,
			   dmat *Phase, int nactive)
#elif DORESTRAINT == 1
static void VectorTorsionR(dmat *crdpop, dmat *frcpop, int atma, int atmb,
			   int atmc, int atmd, double* theta, Energy* sysUV,
			   int nterm, double* Kamp, dmat *Phase, double fbhw,
			   double quadwin, int nactive)
#elif DORESTE == 1
static void VectorTorsionRU(dmat *crdpop, int atma, int atmb, int atmc,
			    int atmd, double* theta, Energy* sysUV, int nterm,
			    double* Kamp, dmat *Phase, double fbhw,
			    double quadwin, int nactive)
#endif
{
  int i, j;
  double costheta, sangle, normfac;
  double fmag, fa, fb1, fc1, fb2, fc2, fd, isinb2, isinc2;
  double mgab, mgbc, mgcd, invab, invbc, invcd, invabc, invbcd, cosb, cosc;
  double Ftop2, afmag, bfmag, cfmag, dfmag;
  double ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  double *phstmp;
  double *xcrda, *xcrdb, *xcrdc, *xcrdd, *ycrda, *ycrdb, *ycrdc, *ycrdd;
  double *zcrda, *zcrdb, *zcrdc, *zcrdd;
  double *xfrca, *xfrcb, *xfrcc, *xfrcd, *yfrca, *yfrcb, *yfrcc, *yfrcd;
  double *zfrca, *zfrcb, *zfrcc, *zfrcd;

  // Compute the maximum force magnitude for convenience
#if DORESTRAINT == 1 || DORESTE == 1
  double dtheta;
  const double Ftop = 2.0*Kamp[0]*quadwin;
  const double quadU = Kamp[0]*quadwin*quadwin;
#endif

  // Loop over all conformations
  xcrda = crdpop->map[3*atma];
  xcrdb = crdpop->map[3*atmb];
  xcrdc = crdpop->map[3*atmc];
  xcrdd = crdpop->map[3*atmd];
  ycrda = crdpop->map[3*atma+1];
  ycrdb = crdpop->map[3*atmb+1];
  ycrdc = crdpop->map[3*atmc+1];
  ycrdd = crdpop->map[3*atmd+1];
  zcrda = crdpop->map[3*atma+2];
  zcrdb = crdpop->map[3*atmb+2];
  zcrdc = crdpop->map[3*atmc+2];
  zcrdd = crdpop->map[3*atmd+2];

  // Set pointers to force arrays if needed
#if DOFORCE == 1 || DORESTRAINT == 1
  xfrca = frcpop->map[3*atma];
  xfrcb = frcpop->map[3*atmb];
  xfrcc = frcpop->map[3*atmc];
  xfrcd = frcpop->map[3*atmd];
  yfrca = frcpop->map[3*atma+1];
  yfrcb = frcpop->map[3*atmb+1];
  yfrcc = frcpop->map[3*atmc+1];
  yfrcd = frcpop->map[3*atmd+1];
  zfrca = frcpop->map[3*atma+2];
  zfrcb = frcpop->map[3*atmb+2];
  zfrcc = frcpop->map[3*atmc+2];
  zfrcd = frcpop->map[3*atmd+2];
#endif

  // Loop over all interactions
  for (i = 0; i < nactive; i++) {
    ab[0] = xcrdb[i] - xcrda[i];
    ab[1] = ycrdb[i] - ycrda[i];
    ab[2] = zcrdb[i] - zcrda[i];
    bc[0] = xcrdc[i] - xcrdb[i];
    bc[1] = ycrdc[i] - ycrdb[i];
    bc[2] = zcrdc[i] - zcrdb[i];
    cd[0] = xcrdd[i] - xcrdc[i];
    cd[1] = ycrdd[i] - ycrdc[i];
    cd[2] = zcrdd[i] - zcrdc[i];
    CrossP(ab, bc, crabbc);
    CrossP(bc, cd, crbccd);

    // Check the angles ABC and BCD: this routine is meant to be used in
    // extreme circumstances which are not commonly encountered in MD, so
    // it is good to check that the dihedral has not linearized.  Bail out
    // if that has happened.
    if (crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2] <
	1.0e-4 ||
	crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2] <
	1.0e-4) {
      continue;
    }
    
    costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
    costheta /=
      sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
           (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
    CrossP(crabbc, crbccd, scr);
    costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
    if (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) {
      theta[i] = acos(costheta);
    }
    else {
      theta[i] = -acos(costheta);
    }

    // This can only happen if the energy alone is being computed;
    // it shortens the calculation even further, to just the angles.
#if DOENERGY == 1
    if (skipintr == 1) {
      continue;
    }
#endif

    // If we're still here, energy and possibly the force must be
    // accumulated.  Either needs a pointer into the array of phases.
    phstmp = Phase->map[i];
#if DOENERGY == 1
    for (j = 0; j < nterm; j++) {
      sangle = Nperiod[j]*theta[i] - phstmp[j];
      sysUV[i].dihe += Kamp[j] * (1.0 + cos(sangle));
    }
#endif
#if DOFORCE == 1
    fmag = 0.0;
    for (j = 0; j < nterm; j++) {
      sangle = Nperiod[j]*theta[i] - phstmp[j];
      fmag += Kamp[j] * Nperiod[j] * sin(sangle);
      sysUV[i].dihe += Kamp[j] * (1.0 + cos(sangle));
    }
#endif

    // For restraints, there is only one term and the energy function
    // is no longer a Fourier series but instead another harmonic well.
#if DORESTRAINT == 1 || DORESTE == 1
    dtheta = theta[i] - phstmp[0];
    if (fabs(dtheta) < fbhw) {
      dtheta = 0.0;
    }
    else if (dtheta < 0.0) {
      dtheta += fbhw;
    }
    else {
      dtheta -= fbhw;
    }

    // Restraint energies get tallied in the reciprocal space
    // electrostatic energy, as this quantity is otherwise unused.
    if (fabs(dtheta) < quadwin) {
      sysUV[i].relec += Kamp[0]*dtheta*dtheta;
#if DORESTRAINT == 1
      fmag = -2.0*Kamp[0]*dtheta;
#endif
    }
    else {
      sysUV[i].relec += quadU + Ftop*(fabs(dtheta) - quadwin);
#if DORESTRAINT == 1
      if (dtheta < 0.0) {
	fmag = Ftop;
      }
      else {
        fmag = -Ftop;
      }
#endif
    }
#endif
#if DOFORCE == 1 || DORESTRAINT == 1
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
    isinb2 *= fmag;
    isinc2 *= fmag;
    invabc = invab*invbc;
    invbcd = invbc*invcd;
    for (j = 0; j < 3; j++) {
      crabbc[j] *= invabc;
      crbccd[j] *= invbcd;
    }

    // Transform the dihedral forces to cartesian coordinates
    fa = -invab * isinb2;
    fb1 = (mgbc - mgab*cosb) * invabc * isinb2;
    fb2 = cosc * invbc * isinc2;
    fc1 = (mgbc - mgcd*cosc) * invbcd * isinc2;
    fc2 = cosb * invbc * isinb2;
    fd = -invcd * isinc2;

    // Apply the dihedral forces
    xfrca[i] += crabbc[0] * fa;
    xfrcb[i] += fb1 * crabbc[0] - fb2 * crbccd[0];
    xfrcc[i] += -fc1 * crbccd[0] + fc2 * crabbc[0];
    xfrcd[i] += -fd * crbccd[0];
    yfrca[i] += crabbc[1] * fa;
    yfrcb[i] += fb1 * crabbc[1] - fb2 * crbccd[1];
    yfrcc[i] += -fc1 * crbccd[1] + fc2 * crabbc[1];
    yfrcd[i] += -fd * crbccd[1];
    zfrca[i] += crabbc[2] * fa;
    zfrcb[i] += fb1 * crabbc[2] - fb2 * crbccd[2];
    zfrcc[i] += -fc1 * crbccd[2] + fc2 * crabbc[2];
    zfrcd[i] += -fd * crbccd[2];
#endif
  }
}

//-----------------------------------------------------------------------------
// VectorCMAP: this computes the CMAP terms of a system with many
//             configurations.
//
// Arguments:
//   crdpop:      matrix of coordinates, giving X1, Y1, Z1, X2, ..., Zn on
//                successive rows
//   frcpop:      matrix of forces, being accumulated before this function is
//                called and accepting contributions from this function
//   atm[abcde]:  the atoms to which this CMAP a-b-c-d-e applies
//-----------------------------------------------------------------------------
#if !defined(DORESTRAINT) && !defined(DORESTE)
#if DOENERGY == 1
static void VectorCMAP(dmat *crdpop, int atma, int atmb, int atmc, int atmd,
                       int atme, Energy* sysUV, cmapcomm *mcom, prmtop *tp,
		       int nactive)
#elif DOFORCE == 1
static void VectorCMAPF(dmat *crdpop, dmat *frcpop, int atma, int atmb,
                        int atmc, int atmd, int atme, Energy *sysUV,
			cmapcomm *mcom, prmtop *tp, int nactive)
#endif
{
  int i, j;
  double aloc[3], bloc[3], cloc[3], dloc[3], eloc[3];
  double afrc[3], bfrc[3], cfrc[3], dfrc[3], efrc[3];

  // Loop over all conformations 
  for (i = 0; i < nactive; i++) {

    // This is not a large chunk of computation, so falling back on the
    // serial implementation just saves a lot of code duplication.
    for (j = 0; j < 3; j++) {
      aloc[j] = crdpop->map[3*atma + j][i];
      bloc[j] = crdpop->map[3*atmb + j][i];
      cloc[j] = crdpop->map[3*atmc + j][i];
      dloc[j] = crdpop->map[3*atmd + j][i];
      eloc[j] = crdpop->map[3*atme + j][i];
      afrc[j] = 0.0;
      bfrc[j] = 0.0;
      cfrc[j] = 0.0;
      dfrc[j] = 0.0;
      efrc[j] = 0.0;
    }
#if DOFORCE == 1
    CmapFrcNrg(aloc, bloc, cloc, dloc, eloc, afrc, bfrc, cfrc, dfrc, efrc,
               mcom, tp, &sysUV[i], 0);
    for (j = 0; j < 3; j++) {
      frcpop->map[3*atma + j][i] += afrc[j];
      frcpop->map[3*atmb + j][i] += bfrc[j];
      frcpop->map[3*atmc + j][i] += cfrc[j];
      frcpop->map[3*atmd + j][i] += dfrc[j];
      frcpop->map[3*atme + j][i] += efrc[j];
    }
#endif
#if DOENERGY == 1
    CmapNrg(aloc, bloc, cloc, dloc, eloc, mcom, tp, &sysUV[i], 0);
#endif
  }
}
#endif
