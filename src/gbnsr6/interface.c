#include "copyright_c.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "prm.h"
#include "interface.h"

// Author: Mengjuei Hsieh, University of California Irvine

void prepb_read_(INT_T*,INT_T*,REAL_T*,REAL_T*,INT_T*);
void pb_read_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
	      int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
	      int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
	      int*,
	REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, 
	REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*,
	REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*,
	REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*,
	REAL_T*, REAL_T*,REAL_T*);
void pb_init_(int*,INT_T*,INT_T*,INT_T*,INT_T*,INT_T*,
        INT_T*,INT_T*,INT_T*,INT_T*,INT_T*,
	INT_T*,INT_T*,INT_T*,INT_T*,
	INT_T*,STRING_T*,STRING_T*,STRING_T*,REAL_T*,
	REAL_T*);
void mypb_force_(INT_T*,INT_T*,INT_T*,
	INT_T*,INT_T*,INT_T*,INT_T*,INT_T*,
	REAL_T*,REAL_T*,REAL_T*,REAL_T*,REAL_T*,REAL_T*,
	REAL_T*,REAL_T*,REAL_T*,REAL_T*);
//void clockstop_();
void pb_free_();

void pboptinit(PBOPTSTRUCT_T *myoption) {
       myoption->epsin         = 1.0;
       myoption->epsout        = 80.0;
       myoption->epsmemb       = 1.0;
       myoption->istrng        = 0.0;
       myoption->pbtemp        = 300.0;
       myoption->dprob         = 1.4;
       myoption->iprob         = 2.0;
       myoption->accept        = 0.001;
       myoption->fillratio     = 2.0;
       myoption->space         = 0.5;
       myoption->arcres        = 0.25;//need to be changed to something consistent
       myoption->cutres        = 99.0;
       myoption->cutfd         = 5.0;
       myoption->cutnb         = 0.0;
       myoption->sprob         = 0.557;
       myoption->vprob         = 1.30;
       myoption->rhow_effect   = 1.129;
       myoption->cavity_surften= 0.0378;
       myoption->cavity_offset = -0.5692;
       myoption->cutsa         = 9.0;
       myoption->fmiccg        = 0.90;
       myoption->ivalence      = 1.0;
       myoption->laccept       = 0.1;
       myoption->wsor          = 1.9;
       myoption->lwsor         = 1.95;
//     myoption->radinc        = myoption->dprob*0.5;
       myoption->radinc        = 0.8;
       myoption->expthresh     = 0.2;
       myoption->offx          = 0.0;
       myoption->offy          = 0.0;
       myoption->offz          = 0.0;
       myoption->sepbuf        = 4.0;
                               
       myoption->smoothopt     = 1;
       myoption->radiopt       = 1;
       myoption->npbopt        = 0;
       myoption->solvopt       = 1;
       myoption->maxitn        = 100;
       myoption->nbuffer       = 0;
       myoption->nfocus        = 2;
       myoption->fscale        = 8;
       myoption->npbgrid       = 1;
       myoption->dbfopt        = 1;
       myoption->bcopt         = 5;
       myoption->scalec        = 0;
       myoption->eneopt        = 2;
       myoption->frcopt        = 0;
       myoption->nsnbr         = 1;
       myoption->nsnba         = 1;
       myoption->phiout        = 0;
       myoption->phiform       = 0;
       myoption->npbverb       = 0;
       myoption->npopt         = 2;
       myoption->decompopt     = 2;
       myoption->use_rmin      = 1;
       myoption->use_sav       = 1;
       myoption->maxsph        = 400;
       myoption->maxarc        = 512;
       myoption->ndofd         = 1;
       myoption->ndosas        = 1;
       myoption->mpopt         = 0;
       myoption->lmax          = 80;
       myoption->inp           = 2;
       myoption->pbprint       = 0;
       myoption->maxarcdot     = 1500;
}

REAL_T epbsa(PBOPTSTRUCT_T *opt, PARMSTRUCT_T *prm, REAL_T *x,
	     REAL_T *grad, REAL_T *evdw, REAL_T *eelt, REAL_T *esurf,
	     REAL_T *edisp, int iteration, int freevectors){
	int i, ifcap, *dummy;
	REAL_T *f;
        REAL_T e_PBSA=0;
	ifcap = 0;
	//checking passed options

	if (freevectors==1 && opt!=NULL) {
           //clockstop_();
           pb_free_();
	   return (e_PBSA);
	} else if ( freevectors==1 ){
	   return (e_PBSA);
	} else if ( iteration < 1 ) {
	  //printf("%d\n",opt->ipb);
	  prepb_read_(&opt->ipb,&opt->inp,prm->Cn1,prm->Cn2,&prm->Nttyp);
	  pb_read_(
                &opt->smoothopt, &opt->radiopt, &opt->npbopt, &opt->solvopt,
		&opt->maxitn, &opt->nbuffer, &opt->nfocus, &opt->fscale,
		&opt->npbgrid, &opt->dbfopt, &opt->bcopt, &opt->scalec,
		&opt->eneopt, &opt->frcopt, &opt->nsnbr, &opt->nsnba,
		&opt->phiout,
		&opt->phiform, &opt->npbverb, &opt->npopt, &opt->decompopt,
		&opt->use_rmin, &opt->use_sav, &opt->maxsph, &opt->maxarc,
		&opt->ndofd, &opt->ndosas, &opt->mpopt, &opt->lmax,
		&opt->maxarcdot,
		&opt->pbprint,
		&opt->epsin, &opt->epsout, &opt->epsmemb, &opt->istrng, &opt->pbtemp,
		&opt->dprob, &opt->iprob, &opt->accept, &opt->fillratio,
		&opt->space, &opt->arcres, &opt->cutres, &opt->cutfd,
		&opt->cutnb, &opt->sprob, &opt->vprob, &opt->rhow_effect,
		&opt->cavity_surften, &opt->cavity_offset, &opt->cutsa,
		&opt->fmiccg, &opt->ivalence, &opt->laccept, &opt->wsor,
		&opt->lwsor, &opt->radinc, &opt->expthresh, &opt->offx,
		&opt->offy, &opt->offz, &opt->sepbuf);
          dummy = (int *) malloc(sizeof(int)*prm->Natom);
          pb_init_(&ifcap,&prm->Natom,&prm->Nres,&prm->Ntypes,&prm->Nbonh,
	     &prm->Nbona,prm->Ipres,prm->Iac,prm->Cno,prm->Iblo,prm->ExclAt,
	     prm->BondHAt1,prm->BondHAt2,prm->BondAt1,prm->BondAt2,
	     dummy,prm->ResNames,prm->AtomNames,prm->AtomSym,prm->Charges,
	     prm->Rborn);
	  free(dummy);
	}

	f = (REAL_T *) malloc(sizeof(REAL_T)*3*prm->Natom);

	mypb_force_(&prm->Natom,&prm->Nres,&prm->Ntypes,
	     prm->Ipres,prm->Iac,prm->Cno,prm->ExclAt,&iteration,
	     prm->Cn1,prm->Cn2,prm->Charges,x,f,&e_PBSA,
	     evdw,eelt,esurf,edisp);

	for (i=0;i<3*prm->Natom;i++){
		grad[i]=grad[i]-f[i];
	}

	free(f);
	return (e_PBSA);
}
/* namelist /pb/ epsin, epsout, smoothopt, istrng, pbtemp,     &
      radiopt, dprob, iprob, npbopt, solvopt, accept, maxitn,  &
      fillratio, space, nbuffer, nfocus, fscale, npbgrid,      &
      arcres,dbfopt,bcopt,scalec,eneopt,frcopt,cutres,cutfd,   &
      cutnb, nsnbr, nsnba,phiout, phiform, npbverb, npopt,     &
      decompopt, use_rmin, sprob, vprob, rhow_effect, use_sav, &
      cavity_surften, cavity_offset, maxsph, maxarc,           &
      cutsa, ndofd, ndosas, fmiccg, ivalence, laccept, wsor,   &
      lwsor, pbkappa, radinc, expthresh, offx, offy, offz,     &
      sepbuf, mpopt, lmax, saopt, intopt, ligandmask, buffer,  &
      xmin, ymin, zmin, xmax, ymax, zmax, isurfchg,            &
      ngrdblkx, ngrdblky, ngrdblkz, saltout, stern,            &
      xmblk, ymblk, zmblk, triopt, sasopt, maxtri, maxarcdot
 */

