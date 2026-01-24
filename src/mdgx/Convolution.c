//-----------------------------------------------------------------------------
// ConvQBC: this convolutes the charge mesh Q with the reciprocal space pair
//          potential BC.  The transformation is done "in place" and so the   
//          charge mesh is overwritten with the electrostatic potential U.  If
//          the "nrg" variant of this function is used, the reciprocal space
//          energy will be returned, and if the "nrgvir" variant of this
//          function is used the reciprocal space energy and virial will be
//          returned.
//
// Arguments:                                                           
//   rcinp:  reciprocal space calculation control data                  
//   crd:    coordinates (used here for cell dimensions and transform   
//           matrices)                                                  
//   Q:      the charge mesh (it becomes the potential mesh after this  
//           function works on it)                                      
//   PPk:    a kit for computing the reciprocal space pair potential    
//-----------------------------------------------------------------------------
#if NEEDENERGY == 0
  #define CONVQBCFUNC ConvQBC
#else
  #if NEEDVIRIAL == 0
    #define CONVQBCFUNC ConvQBCnrg
  #else
    #define CONVQBCFUNC ConvQBCnrgvir
  #endif
#endif

#if NEEDENERGY == 0
void CONVQBCFUNC(reccon *rcinp, coord *crd, dbook *Q, bckit *PPk,
		 execon *etimers)
#else
void CONVQBCFUNC(reccon *rcinp, coord *crd, dbook *Q, bckit *PPk,
		 Energy *sysUV, execon *etimers)
#endif
{
  int i, j, k, ng2l, kmin;
  int ng[3];
  double dx, dy, dz, dx2, dxy2, dxyz2, SPx, SPxy, S, pivol, tpS;
  double mmx1, mmy1, mmz1, mmx2, mmy2, mmz2, mmx3, mmy3, mmz3, mm2;
  double e2axyz, ea2xyz, e2abxyz, bcfac;
  double *dtmp, *utmp;
  double *tBx, *tBy, *tBz, *tEx, *tEy, *tEz, *tmx, *tmy, *tmz;
  double **dtm2p;

  // Shortcuts
  ng[0] = Q->pag;
  ng[1] = Q->row;
  ng[2] = Q->col;
  ng2l = ng[2]/2+1;
  utmp = crd->U.data;
  S = rcinp->S;
  tBx = PPk->Bx;
  tBy = PPk->By;
  tBz = PPk->Bz;
  tEx = PPk->Ex;
  tEy = PPk->Ey;
  tEz = PPk->Ez;
  tmx = PPk->mx;
  tmy = PPk->my;
  tmz = PPk->mz;

#if NEEDENERGY == 1
  double relec = 0.0;

#if NEEDVIRIAL == 1
  double mhati, mhatj, mhatk, vterm;
  double fac2 = 8.0 * PI * PI * rcinp->S * rcinp->S;
  double vir00 = 0.0;
  double vir11 = 0.0;
  double vir22 = 0.0;
  double vir10 = 0.0;
  double vir20 = 0.0;
  double vir21 = 0.0;
#endif
#endif

  // Forward transform to create fQ
  fftw_execute(PPk->forwplan);
  etimers->nbFFT += mdgxStopTimer(etimers);

  // Neutralize the system
  Q->data[0] = 0.0;
  Q->data[1] = 0.0;

  // Create the BC mesh
  if (crd->isortho == 1) {

    // If this unit cell is orthorhombic, the computation is easy
    for (i = 0; i < ng[0]; i++) {
      dx2 = tmx[i]*tmx[i];
      dtm2p = Q->map[i];
      SPx = tBx[i]*tEx[i];
#if NEEDVIRIAL == 1
      mhati = tmx[i];
#endif
      for (j = 0; j < ng[1]; j++) {
        dxy2 = dx2 + tmy[j]*tmy[j];
        SPxy = SPx*tBy[j]*tEy[j];
        dtmp = dtm2p[j];
        kmin = (j == 0 && i == 0) ? 1 : 0;
#if NEEDVIRIAL == 1
	mhatj = tmy[j];
#endif
        for (k = kmin; k < ng2l; k++) {
          dxyz2 = dxy2 + tmz[k]*tmz[k];
          bcfac = SPxy*tBz[k]*tEz[k]/dxyz2;
#if NEEDENERGY == 1
  #if NEEDVIRIAL == 1
	  mhatk = tmz[k];
  #endif
	  int twok = 2*k;
	  double qreal = dtmp[twok];
	  double qimag = dtmp[twok+1];
	  double qfac = qreal*qreal + qimag*qimag;
	  dtmp[twok] = qreal*bcfac;
	  dtmp[twok+1] = qimag*bcfac;
  #if NEEDVIRIAL == 1
	  vterm = (fac2 + 2.0/dxyz2) * bcfac * qfac;
  #endif
	  if (k == 0) {
	    qfac *= 0.5;
  #if NEEDVIRIAL == 1
	    vterm *= 0.5;
  #endif
	  }
	  relec += qfac*bcfac;
  #if NEEDVIRIAL == 1
	  vir00 += vterm*mhati*mhati - bcfac*qfac;
	  vir11 += vterm*mhatj*mhatj - bcfac*qfac;
	  vir22 += vterm*mhatk*mhatk - bcfac*qfac;
	  vir10 += vterm*mhati*mhatj;
	  vir20 += vterm*mhati*mhatk;
	  vir21 += vterm*mhatj*mhatk;
  #endif
#else
          dtmp[2*k] *= bcfac;
          dtmp[2*k+1] *= bcfac;
#endif
        }
      }
    }
  }
  else {

    // The unit cell is not orthorhombic
    pivol = BIOQ/(PI*crd->invU.data[0]*crd->invU.data[4]*crd->invU.data[8]);
    tpS = 2.0*PI*S;
    tpS *= tpS;

    // Pre-compute some exponential factors
    dz = tpS*(utmp[6]*utmp[6] + utmp[7]*utmp[7] + utmp[8]*utmp[8]);
    for (i = 0; i < ng2l; i++) {
      tEz[i] = exp(-dz*tmz[i]*tmz[i]);
    }
    for (i = 0; i < ng[0]; i++) {

      // We multiply the [dx dy dz] vector by the transpose of U
      dx = tmx[i];
      mmx1 = utmp[0]*dx;
      mmy1 = utmp[1]*dx;
      mmz1 = utmp[2]*dx;
      SPx = pivol*tBx[i];
      dtm2p = Q->map[i];
      for (j = 0; j < ng[1]; j++) {
        dy = tmy[j];
        SPxy = SPx*tBy[j];
        mmx2 = mmx1 + utmp[3]*dy;
        mmy2 = mmy1 + utmp[4]*dy;
        mmz2 = mmz1 + utmp[5]*dy;
        dtmp = dtm2p[j];

        // Pre-computing these two exponentials relieves me of
        // computing them in the inner loop.                  
        ea2xyz = exp(-tpS*(mmx2*mmx2 + mmy2*mmy2 + mmz2*mmz2));
        e2abxyz = 1.0;
        e2axyz = exp(-2.0*tpS*(mmx2*utmp[6]+mmy2*utmp[7]+mmz2*utmp[8]));
        kmin = (j == 0 && i == 0) ? 1 : 0;
        for (k = kmin; k < ng2l; k++) {
          dz = tmz[k];
          mmx3 = mmx2 + utmp[6]*dz;
          mmy3 = mmy2 + utmp[7]*dz;
          mmz3 = mmz2 + utmp[8]*dz;
          mm2 = mmx3*mmx3 + mmy3*mmy3 + mmz3*mmz3;
          bcfac = SPxy*tBz[k] * (ea2xyz*e2abxyz) * tEz[k] / mm2;
#if NEEDENERGY == 1
          int twok = 2*k;
          double qreal = dtmp[twok];
          double qimag = dtmp[twok+1];
          double qfac = qreal*qreal + qimag*qimag;
          relec = (k > 0) ? relec + qfac*bcfac : relec + 0.5*qfac*bcfac;
          dtmp[twok] = qreal*bcfac;
          dtmp[twok+1] = qimag*bcfac;
#else
          dtmp[2*k] *= bcfac;
          dtmp[2*k+1] *= bcfac;
#endif
          e2abxyz *= e2axyz;
        }
      }
    }
  }
  etimers->nbCnv += mdgxStopTimer(etimers);

  // Perform the "backward," complex-to-real transform to create U
  fftw_execute(PPk->backplan);
  etimers->nbFFT += mdgxStopTimer(etimers);

#if NEEDENERGY == 1
  sysUV->relec = relec + PPk->SelfEcorr;
#if NEEDVIRIAL == 1
  sysUV->Vir[0] += vir00;
  sysUV->Vir[1] += vir10;
  sysUV->Vir[2] += vir20;
  sysUV->Vir[3] += vir10;
  sysUV->Vir[4] += vir11;
  sysUV->Vir[5] += vir21;
  sysUV->Vir[6] += vir20;
  sysUV->Vir[7] += vir21;
  sysUV->Vir[8] += vir22;
#endif

#endif
}

#undef CONVQBCFUNC
