/* double
 * epsin, epsout, istrng, pbtemp, dprob, iprob, accept,
 * fillratio, space, arcres, cutres, cutfd, cutnb, sprob,
 * vprob, rhow_effect, cavity_surften, cavity_offset, cutsa,
 * fmiccg, ivalence, laccept, wsor, lwsor, pbkappa, radinc,
 * expthresh, offx, offy, offz, sepbuf, buffer, stern
 * int
 * smoothopt, radiopt, npbopt, solvopt, maxitn, nbuffer,
 * nfocus, fscale, npbgrid, dbfopt, bcopt, scalec, eneopt,
 * frcopt, nsnbr, nsnba, phiout, phiform, npbverb, npopt, decompopt,
 * use_rmin, use_sav, maxsph, maxarc, ndofd, ndosas, mpopt,
 * lmax, saopt, intopt, isurfchg, saltout, triopt, sasopt, maxtri,
 * maxarcdot
 */
/*
 * LIGAND / MULTIBLOCK, not for NAB
 * double xmin, ymin, zmin, xmax, ymax, zmax
 * int    ngrdblkx, ngrdblky, ngrdblkz, xmblk, ymblk, zmblk
 * char * ligandmask
 */
typedef struct {
    REAL_T epsin, epsout, epsmemb, istrng, pbtemp, dprob, iprob, accept,
           fillratio, space, arcres, cutres, cutfd, cutnb, sprob,
           vprob, rhow_effect, cavity_surften, cavity_offset, cutsa,
           fmiccg, ivalence, laccept, wsor, lwsor, pbkappa, radinc,
           expthresh, offx, offy, offz, sepbuf, buffer, stern;
    int    smoothopt, radiopt, npbopt, solvopt, maxitn, nbuffer,
           nfocus, fscale, npbgrid, dbfopt, bcopt, scalec, eneopt,
           frcopt, nsnbr, nsnba, phiout, phiform, npbverb, npopt, decompopt,
           use_rmin, use_sav, maxsph, maxarc, ndofd, ndosas, mpopt,
           lmax, saopt, intopt, isurfchg, saltout, triopt, sasopt,
           maxtri, maxarcdot,
           inp, pbprint, ipb;
} PBOPTSTRUCT_T;

void pboptinit();

//PBOPTSTRUCT_T* pbopt;

REAL_T epbsa(PBOPTSTRUCT_T *opt, PARMSTRUCT_T *prm, REAL_T *x,
	REAL_T *grad, REAL_T *evdw, REAL_T *eelt, REAL_T *esurf,
	REAL_T *edisp, int iteration, int freevectors);

