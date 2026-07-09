#include <stdbool.h>
#ifdef BINTRAJ

#include        "basics.h"
#include        "defaults.h"
#include        "classes.h"
#include        "residue.h"
#include        "unitio.h"
#include        "netcdf.h"
#include        "avl.h"
#include        "sort.h"
#include        "cmap.h"
#include <assert.h>

#define LABLEN 4
#define MAXDESCLEN MAXSTRINGLENGTH

/* ================================================================
   unitio_netcdf.c - NetCDF prmtop writer for LEAP.

   Design:
   - NC_NETCDF4 mode makes nc_redef() cheap; each FLAG section is
     self-contained: redef → def_var → annotate → enddef → write.
   - No POINTERS or SOLVENT_POINTERS - sizes implicit in dimensions.
   - Torsion/CMAP atom indices stored 1-based, no AMBERINDEX multiply.
     Sign flags preserved in atom3/atom4 for torsions.
     Reader removes the /3 but sign extraction is unchanged.
   - CMAP_INDEX raw indices stored (no /3+1); reader applies at runtime.

   Reader-side changes required in sander/pmemd:
   - Torsion indices: remove AMBERINDEX inverse (the /3).
     Sign extraction for bCalc14/bProper unchanged.
   - CMAP_INDEX: remove /3+1 transform.

   Function order (forward-reference safe):
     1. Macros
     2. gbrad_for_element()
     3. gbscreen_for_element()
     4. write_mass_radii_screen()
     5. write_polar_section()
     6. write_torsions()
     7. SaveAmberParmCMAPNetcdf()
     8. UnitIOSaveAmberParmNetcdf()  -- in separate file unitio_netcdf_parm.c
     9. zUnitIOSaveAmberParmFormat() -- fork point
   ================================================================ */

#include <netcdf.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* ----------------------------------------------------------------
   Original macros from unitio.h
   RESTRAINTLOOP calls FortranWriteDouble so cannot be used in the
   NetCDF path; use NC_RESTRAINTLOOP instead.
   bPERT_* are clean booleans, used as-is.
   ---------------------------------------------------------------- */
/* NC_RESTRAINTLOOP: NetCDF equivalent of RESTRAINTLOOP.
   Collects values into buf[] at offset, preserves iParmIndex
   side effect identically to the original.                              */
#define NC_RESTRAINTLOOP(uUnit, type, field, buf, offset)                \
    do {                                                                 \
        int _ii, _iiMax, _jj = 0;                                       \
        SAVERESTRAINTt *_sr;                                             \
        if ((_iiMax = iVarArrayElementCount((uUnit)->vaRestraints))) {   \
            _sr = PVAI((uUnit)->vaRestraints, SAVERESTRAINTt, 0);        \
            for (_ii = 0; _ii < _iiMax; _ii++, _sr++) {                 \
                if (_sr->iType == (type)) {                              \
                    (buf)[(offset) + _jj] = _sr->field;                 \
                    _sr->iParmIndex = (offset) + _jj;                   \
                    _jj++;                                               \
                }                                                        \
            }                                                            \
        }                                                                \
    } while(0)

#define bPERT_BOND(bp,a1,a2) \
    (bp && (bAtomFlagsSet(a1,ATOMPERTURB) || bAtomFlagsSet(a2,ATOMPERTURB)))

#define bPERT_ANGLE(bp,a1,a2,a3) \
    (bp && (bAtomFlagsSet(a1,ATOMPERTURB) || bAtomFlagsSet(a2,ATOMPERTURB) \
         || bAtomFlagsSet(a3,ATOMPERTURB)))

#define bPERT_TORSION(bp,a1,a2,a3,a4) \
    (bp && (bAtomFlagsSet(a1,ATOMPERTURB) || bAtomFlagsSet(a2,ATOMPERTURB) \
         || bAtomFlagsSet(a3,ATOMPERTURB) || bAtomFlagsSet(a4,ATOMPERTURB)))


/* ----------------------------------------------------------------
   NC_CHECK: error handler for NetCDF API calls.
   ---------------------------------------------------------------- */
#define NC_CHECK(call) \
    do { int _e = (call); if (_e != NC_NOERR) { \
        fprintf(stderr, "NetCDF error %s: %s (%s:%d)\n", #call, nc_strerror(_e), __FILE__,__LINE__); \
        return _e; \
    }} while(0)


/* ================================================================
   3. write_mass_radii_screen()
   Merged single ParmSetAtom pass for MASS, RADII, SCREEN.
   ================================================================ */
static int write_mass_radii_screen(int ncid, UNITt *uUnit, int dimid_atoms)
{
    int n = iVarArrayElementCount(uUnit->vaAtoms);
    double *mass   = malloc(n * sizeof(double));
    double *radii  = malloc(n * sizeof(double));
    double *screen = malloc(n * sizeof(double));
    char sType[MAXTYPELEN], sDesc[PARMDESCRIPTIONLEN];
    double dMass, dPolar, dEpsilon, dRStar, dEpsilon14, dRStar14, dScreenF;
    int iElement, iHybridization;
    int vid_mass, vid_radii, vid_screen;

    SAVEATOMt *sa = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
    for (int i = 0; i < n; i++) {
        int iIdx = iParmSetFindAtom(uUnit->psParameters, sa[i].sType);
        ParmSetAtom(uUnit->psParameters, iIdx, sType, &dMass, &dPolar,
                    &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
                    &iElement, &iHybridization, sDesc);
        mass[i]   = dMass;
        radii[i]  = dGBRadiusForAtom(&sa[i], iElement, dMass, (i+1==n) );
        screen[i] = dGBScreenForElement(iElement);
    }

    NC_CHECK(nc_redef(ncid));
    NC_CHECK(nc_def_var(ncid, "MASS", NC_DOUBLE, 1, &dimid_atoms, &vid_mass));
    NC_CHECK(nc_put_att_text(ncid, vid_mass, "units", 5, "g/mol"));

    NC_CHECK(nc_def_var(ncid, "RADII", NC_DOUBLE, 1, &dimid_atoms, &vid_radii));
    NC_CHECK(nc_put_att_text(ncid, vid_radii, "units", 8, "angstrom"));
    NC_CHECK(nc_put_att_text(ncid, vid_radii, "long_name", 9, "GB radius"));

    NC_CHECK(nc_def_var(ncid, "SCREEN", NC_DOUBLE, 1, &dimid_atoms, &vid_screen));
    NC_CHECK(nc_put_att_text(ncid, vid_screen, "long_name", 20, "GB screening factor"));
    NC_CHECK(nc_enddef(ncid));

    NC_CHECK(nc_put_var_double(ncid, vid_mass, mass));
    NC_CHECK(nc_put_var_double(ncid, vid_radii, radii));
    NC_CHECK(nc_put_var_double(ncid, vid_screen, screen));

    free(mass); free(radii); free(screen);
    return NC_NOERR;
}

/* ================================================================
   4. write_polar_section()
   Merged single ParmSetAtom pass for POLARIZABILITY,
   DIPOLE_DAMP_FACTOR, and optionally PERT_POLARIZABILITY.
   ================================================================ */
static int write_polar_section(int ncid, UNITt *uUnit,
                               int dimid_atoms, bool bPert)
{
    int i, n = iVarArrayElementCount(uUnit->vaAtoms);
    int iCount = 0, iCountPerturbed = 0, iPertTot = 0;
    double *polar = malloc(n * sizeof(double));
    double *damp  = malloc(n * sizeof(double));
    double *pper  = bPert ? malloc(n * sizeof(double)) : NULL;
    char sType[MAXTYPELEN], sDesc[MAXDESCLEN];
    double dMass, dPolar, dEpsilon, dRStar, dEpsilon14, dRStar14, dScreenF;
    int iElement, iHybridization;
    int vid;
    SAVEATOMt *a = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);

    for (i = 0; i < n; i++, a++) {
        int iIdx = iParmSetFindAtom(uUnit->psParameters, a->sType);
        ParmSetAtom(uUnit->psParameters, iIdx, sType, &dMass, &dPolar,
                    &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
                    &iElement, &iHybridization, sDesc);
        if (dPolar == -1.0) { dPolar = 0.0; iCount++; }
        polar[i] = dPolar;
        damp[i]  = (dScreenF == 0.0) ? GDefaults.dDipoleDampFactor : dScreenF;
        if (bPert) {
            bool bTmp = bAtomPerturbed(a->aAtom);
            if (bTmp) {
                int iPIdx = iParmSetFindAtom(uUnit->psParameters, a->sPertType);
                ParmSetAtom(uUnit->psParameters, iPIdx, sType, &dMass, &dPolar,
                            &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
                            &iElement, &iHybridization, sDesc);
                if (dPolar == -1.0) { dPolar = 0.0; iCountPerturbed++; }
                iPertTot++;
            }
            pper[i] = dPolar;
        }
    }
    if (iCount > 0)
        VP0("Total atoms with default polarization=0.0: %d of %d\n", iCount, n);

    NC_CHECK(nc_def_var(ncid, "POLARIZABILITY", NC_DOUBLE, 1, &dimid_atoms, &vid));
    NC_CHECK(nc_put_att_text(ncid, vid, "units", 10, "angstrom^3"));
    NC_CHECK(nc_enddef(ncid));
    NC_CHECK(nc_put_var_double(ncid, vid, polar));

    NC_CHECK(nc_def_var(ncid, "DIPOLE_DAMP_FACTOR", NC_DOUBLE, 1, &dimid_atoms, &vid));
    NC_CHECK(nc_enddef(ncid));
    NC_CHECK(nc_put_var_double(ncid, vid, damp));

    if (bPert) {
        if (iCountPerturbed > 0)
            VP0("Total pert atoms with default polarization=0.0: %d of %d\n",
                iCountPerturbed, iPertTot);
        NC_CHECK(nc_def_var(ncid, "PERT_POLARIZABILITY", NC_DOUBLE, 1, &dimid_atoms, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 10, "angstrom^3"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_double(ncid, vid, pper));
        free(pper);
    }
    free(polar); free(damp);
    return NC_NOERR;
}

/* ================================================================
   5. write_torsions()
   SWAP preserved, sign flags preserved in atom3/atom4,
   1-based indices stored directly (no AMBERINDEX multiply).
   ================================================================ */
static int write_torsions(int ncid, UNITt *uUnit, bool bPert, int dimid_atom4)
{
    int nh = 0, nnh = 0;
    int dimid_dih_h, dimid_dih_nh;
    int vid_dih_h_atoms, vid_dih_h_parm;
    int vid_dih_nh_atoms, vid_dih_nh_parm;
    int dims2[2];
    int nT = iVarArrayElementCount(uUnit->vaTorsions);
    int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION);
    int (*ah)[4]  = malloc(4*(nT+nr)*sizeof(int));
    int *pih  = malloc(nT*sizeof(int));
    int (*anh)[4] = malloc(4*nT*sizeof(int));
    int *pinh = malloc(nT*sizeof(int));

    for (int i = 0; i < nT; i++) {
        SAVETORSIONt *t = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
        ATOMt *aA = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom1-1)->aAtom;
        ATOMt *aB = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom2-1)->aAtom;
        ATOMt *aC = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom3-1)->aAtom;
        ATOMt *aD = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom4-1)->aAtom;
        if (bPERT_TORSION(bPert, aA, aB, aC, aD)) continue;
        if (t->iAtom3 == 1 || t->iAtom4 == 1) {
            int iTemp;
            MESSAGE("Had to turn torsion around to avoid K,L == 0\n");
            SWAP(t->iAtom1, t->iAtom4, iTemp);
            SWAP(t->iAtom2, t->iAtom3, iTemp);
        }
        int iCalc14 = t->bCalc14 ? 1 : -1;
        int iProper  = t->bProper  ? 1 : -1;
        if (iAtomElement(aA)==HYDROGEN || iAtomElement(aB)==HYDROGEN ||
                     iAtomElement(aC)==HYDROGEN || iAtomElement(aD)==HYDROGEN) {
            ah[nh][0]=t->iAtom1; ah[nh][1]=t->iAtom2;
            ah[nh][2]=t->iAtom3*iCalc14; ah[nh][3]=t->iAtom4*iProper;
            pih[nh]=t->iParmIndex; nh++;
        } else {
            anh[nnh][0]=t->iAtom1; anh[nnh][1]=t->iAtom2;
            anh[nnh][2]=t->iAtom3*iCalc14; anh[nnh][3]=t->iAtom4*iProper;
            pinh[nnh]=t->iParmIndex; nnh++;
        }
    }
    if (nr) {
        int iMax = iVarArrayElementCount(uUnit->vaRestraints);
        SAVERESTRAINTt *srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        for (int i = 0; i < iMax; i++, srPRestraint++) {
            if (srPRestraint->iType == RESTRAINTTORSION) {
                if ((srPRestraint->iAtom3 == 1) || (srPRestraint->iAtom4 == 1)) {
                    int iTemp;
                    SWAP(srPRestraint->iAtom1, srPRestraint->iAtom4, iTemp);
                    SWAP(srPRestraint->iAtom2, srPRestraint->iAtom3, iTemp);
                }
                anh[nnh][0]=srPRestraint->iAtom1;
                anh[nnh][1]=srPRestraint->iAtom2;
                anh[nnh][2]=srPRestraint->iAtom3;
                anh[nnh][3]=srPRestraint->iAtom4;
                pinh[nnh]=srPRestraint->iParmIndex;
                nnh++;
            }
        }
    }

    NC_CHECK(nc_redef(ncid));
    NC_CHECK(nc_def_dim(ncid, "NPHIH", nh,  &dimid_dih_h));
    NC_CHECK(nc_def_dim(ncid, "MPHIA", nnh, &dimid_dih_nh));

    dims2[0]=dimid_dih_h; dims2[1]=dimid_atom4;
    NC_CHECK(nc_def_var(ncid, "DIHEDRALS_INC_HYDROGEN_ATOMS", NC_INT, 2, dims2, &vid_dih_h_atoms));
    NC_CHECK(nc_put_att_text(ncid, vid_dih_h_atoms, "long_name", 20, "1-based atom indices"));
    NC_CHECK(nc_put_att_text(ncid, vid_dih_h_atoms, "reference_dimension", 5, "NATOM"));
    NC_CHECK(nc_put_att_text(ncid, vid_dih_h_atoms, "note", 46,
                             "atom3,atom4 sign encodes bCalc14,bProper flags"));

    NC_CHECK(nc_def_var(ncid, "DIHEDRALS_INC_HYDROGEN_PARM", NC_INT, 1, &dimid_dih_h, &vid_dih_h_parm));
    NC_CHECK(nc_put_att_text(ncid, vid_dih_h_parm, "long_name", 18, "1-based parm index"));
    NC_CHECK(nc_put_att_text(ncid, vid_dih_h_parm, "reference",
                             59, "DIHEDRAL_FORCE_CONSTANT DIHEDRAL_PERIODICITY DIHEDRAL_PHASE"));

    dims2[0]=dimid_dih_nh; dims2[1]=dimid_atom4;
    NC_CHECK(nc_def_var(ncid, "DIHEDRALS_WITHOUT_HYDROGEN_ATOMS", NC_INT, 2, dims2, &vid_dih_nh_atoms));
    NC_CHECK(nc_put_att_text(ncid, vid_dih_nh_atoms, "long_name", 20, "1-based atom indices"));
    NC_CHECK(nc_put_att_text(ncid, vid_dih_nh_atoms, "reference_dimension", 5, "NATOM"));
    NC_CHECK(nc_put_att_text(ncid, vid_dih_nh_atoms, "note", 47,
                             "atom3/atom4 sign encodes bCalc14/bProper flags"));

    NC_CHECK(nc_def_var(ncid, "DIHEDRALS_WITHOUT_HYDROGEN_PARM", NC_INT, 1, &dimid_dih_nh, &vid_dih_nh_parm));
    NC_CHECK(nc_put_att_text(ncid, vid_dih_nh_parm, "long_name", 18, "1-based parm index"));
    NC_CHECK(nc_put_att_text(ncid, vid_dih_nh_parm, "reference",
                             59, "DIHEDRAL_FORCE_CONSTANT DIHEDRAL_PERIODICITY DIHEDRAL_PHASE"));
    NC_CHECK(nc_enddef(ncid));

    NC_CHECK(nc_put_var_int(ncid, vid_dih_h_atoms,  (int*)ah));
    NC_CHECK(nc_put_var_int(ncid, vid_dih_h_parm,   pih));
    NC_CHECK(nc_put_var_int(ncid, vid_dih_nh_atoms, (int*)anh));
    NC_CHECK(nc_put_var_int(ncid, vid_dih_nh_parm,  pinh));

    free(ah); free(pih);
    free(anh); free(pinh);
    return NC_NOERR;
}

/* ================================================================
   6. SaveAmberParmCMAPNetcdf() -- see unitio_cmap.c
   ================================================================ */

/* ================================================================
   7. UnitIOSaveAmberParmNetcdf()
   ================================================================ */
int UnitIOSaveAmberParmNetcdf(const char *fname, UNITt *uUnit,
                              bool bPert, bool bPolar,
                              VARARRAY vaExcludedAtoms,
                              VARARRAY vaExcludedCount,
                              VARARRAY vaNBIndexMatrix,
                              VARARRAY vaNBParameters,
                              VARARRAY vaNBIndex)
{
    int ncid;
    int dimid_atoms, dimid_res, dimid_name4;
    int dimid_bp;   /* NUMBND - shared by bond parm and augmentation sections */
    int dimid_ap;   /* NUMANG - shared by angle parm and CHARMM UB sections */
    int dimid_atom2, dimid_atom3, dimid_atom4;
    int n_atoms    = iVarArrayElementCount(uUnit->vaAtoms);
    int n_residues = iVarArrayElementCount(uUnit->vaResidues);
    char s1[8], s2[8], s3[8], s4[8], sd[80];
    double dKb, dR0, dKpull, dRpull0, dKpress, dRpress0;
    double dKt, dT0, dTkub, dRkub;
    double dKp, dP0, dScEE, dScNB;
    double dC, dD;
    int vid;

    NC_CHECK(nc_create(fname, NC_CLOBBER|NC_64BIT_OFFSET|NC_NOFILL, &ncid));

    /* Global attributes - AMBER NetCDF convention */
    {
        time_t tp; char sDate[64];
        time(&tp);
        strftime(sDate, sizeof(sDate), "%Y-%m-%dT%H:%M:%S", localtime(&tp));
        NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "date", strlen(sDate), sDate));
        NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "date_conventions", 8, "ISO 8601"));
    }
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "title",
                             strlen(sContainerName(uUnit)), sContainerName(uUnit)));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "application",        5, "AMBER"));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "applicationVersion", 2, "26"));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "program",            4, "leap"));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "programVersion",     3, "1.0"));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "Conventions",       11, "AMBERPRMTOP"));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "ConventionVersion",   3, "1.0"));

    /* Shared dimensions reused across many sections */
    NC_CHECK(nc_def_dim(ncid, "NTOTAT",      n_atoms,    &dimid_atoms));
    NC_CHECK(nc_def_dim(ncid, "NTOTRS",      n_residues, &dimid_res));
    NC_CHECK(nc_def_dim(ncid, "LABLEN",      LABLEN,     &dimid_name4));
    NC_CHECK(nc_def_dim(ncid, "atom_pair",   2,          &dimid_atom2));
    NC_CHECK(nc_def_dim(ncid, "atom_triple", 3,          &dimid_atom3));
    NC_CHECK(nc_def_dim(ncid, "atom_quad",   4,          &dimid_atom4));
    NC_CHECK(nc_enddef(ncid));

    if (iVarArrayElementCount(uUnit->vaRestraints)) {
        VP0(" Restraints:  Bond %d  Angle %d  Torsion %d\n",
            iUnitRestraintTypeCount(uUnit, RESTRAINTBOND),
            iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE),
            iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION));
    } else VP0(" (no restraints)\n");

    /* cell dimensions */
    if (bUnitUseBox(uUnit)) {
        int vid_a, vid_b, vid_c, vid_alpha, vid_beta, vid_gamma;
        double dAlpha = uUnit->dAlpha / DEGTORAD;
        double dBeta  = uUnit->dBeta  / DEGTORAD;
        double dGamma = uUnit->dGamma / DEGTORAD;
        double dX, dY, dZ;
        UnitGetBox(uUnit, &dX, &dY, &dZ);

        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "cell_length_a",    NC_DOUBLE, 0, NULL, &vid_a));
        NC_CHECK(nc_put_att_text(ncid, vid_a,     "units", 8, "angstrom"));
        NC_CHECK(nc_def_var(ncid, "cell_length_b",    NC_DOUBLE, 0, NULL, &vid_b));
        NC_CHECK(nc_put_att_text(ncid, vid_b,     "units", 8, "angstrom"));
        NC_CHECK(nc_def_var(ncid, "cell_length_c",    NC_DOUBLE, 0, NULL, &vid_c));
        NC_CHECK(nc_put_att_text(ncid, vid_c,     "units", 8, "angstrom"));
        NC_CHECK(nc_def_var(ncid, "cell_angle_alpha", NC_DOUBLE, 0, NULL, &vid_alpha));
        NC_CHECK(nc_put_att_text(ncid, vid_alpha, "units", 6, "degree"));
        NC_CHECK(nc_def_var(ncid, "cell_angle_beta",  NC_DOUBLE, 0, NULL, &vid_beta));
        NC_CHECK(nc_put_att_text(ncid, vid_beta,  "units", 6, "degree"));
        NC_CHECK(nc_def_var(ncid, "cell_angle_gamma", NC_DOUBLE, 0, NULL, &vid_gamma));
        NC_CHECK(nc_put_att_text(ncid, vid_gamma, "units", 6, "degree"));
        NC_CHECK(nc_enddef(ncid));

        NC_CHECK(nc_put_var_double(ncid, vid_a,     &dX));
        NC_CHECK(nc_put_var_double(ncid, vid_b,     &dY));
        NC_CHECK(nc_put_var_double(ncid, vid_c,     &dZ));
        NC_CHECK(nc_put_var_double(ncid, vid_alpha, &dAlpha));
        NC_CHECK(nc_put_var_double(ncid, vid_beta,  &dBeta));
        NC_CHECK(nc_put_var_double(ncid, vid_gamma, &dGamma));
    }

    /* -3- atom names */
    {
        int n = iVarArrayElementCount((uUnit)->vaAtoms);
        char *b = malloc(n * LABLEN);
        int vid;
        int dims2[2] = {dimid_atoms, dimid_name4};
        memset(b, ' ', n * LABLEN);
        for (int i = 0; i < n; i++) {
            const char *s = PVAI((uUnit)->vaAtoms,SAVEATOMt,i)->sName;
            size_t l = strlen(s);
            if (l > (size_t)LABLEN) s += l - LABLEN;
            memcpy(b+i*LABLEN, s, l < (size_t)LABLEN ? l : LABLEN);
        }
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "ATOM_NAME", NC_CHAR, 2, dims2, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 9, "Atom name"));
        NC_CHECK(nc_put_att_text(ncid, vid, "pdb_field",     4, "name"));
        if (GDefaults.bCIFReadAuth)
            NC_CHECK(nc_put_att_text(ncid, vid, "mmcif_field",    22, "atom_site.auth_atom_id"));
        else NC_CHECK(nc_put_att_text(ncid, vid, "mmcif_field",  23, "atom_site.label_atom_id"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_text(ncid, vid, b));
        free(b);
    }

/* ================================================================
   1. gbrad_for_element()
   ================================================================ */

    /* -4- charges (raw double, no AMBER unit scaling on write) */
    {
        int n = n_atoms;
        double *buf = malloc(n * sizeof(double));
        for (int i = 0; i < n; i++)
            buf[i] = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->dCharge;
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "CHARGE", NC_DOUBLE, 1, &dimid_atoms, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 18, "elementary charge"));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 51,
                                 "multiply by 18.2223 for AMBER internal charge units"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_double(ncid, vid, buf));
        free(buf);
    }

    /* -4b- atomic numbers */
    {
        int *buf = malloc(n_atoms * sizeof(int));
        for (int i = 0; i < n_atoms; i++)
            buf[i] = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iElement;
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "ATOMIC_NUMBER", NC_INT, 1, &dimid_atoms, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 13, "atomic number"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, buf));
        free(buf);
    }

    /* -5-, RADII, SCREEN: merged single ParmSetAtom pass */
    NC_CHECK(write_mass_radii_screen(ncid, uUnit, dimid_atoms));

    /* -6- atom type index */
    {
        int *buf = malloc(n_atoms * sizeof(int));
        for (int i = 0; i < n_atoms; i++) {
            int iAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iTypeIndex - 1;
            buf[i] = *PVAI(vaNBIndex, int, iAtom) + 1;
        }
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "ATOM_TYPE_INDEX", NC_INT, 1, &dimid_atoms, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name",          18, "1-based type index"));
        NC_CHECK(nc_put_att_text(ncid, vid, "reference_dimension", 7, "NTYPES2"));
        NC_CHECK(nc_put_att_text(ncid, vid, "reference",          20, "NONBONDED_PARM_INDEX"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, buf));
        free(buf);
    }

    /* -7- number of excluded atoms per atom */
    {
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "NUMBER_EXCLUDED_ATOMS", NC_INT, 1, &dimid_atoms, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 28,
                                 "excluded atom count per atom"));
        NC_CHECK(nc_put_att_text(ncid, vid, "encoding",  21,
                                 "cumulative_sum_offset"));
        NC_CHECK(nc_put_att_text(ncid, vid, "reference", 19,
                                 "EXCLUDED_ATOMS_LIST"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, PVAI(vaExcludedCount, int, 0)));
    }

    /* -8- nonbonded parm index matrix */
    int nttyp = iVarArrayElementCount(vaNBParameters);
    int ntypes = (sqrt(8*nttyp+1)-1)/2;
    int dimid_nb;
    {
        int ntypes2 = iVarArrayElementCount(vaNBIndexMatrix);
        assert(ntypes*ntypes == ntypes2);
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_dim(ncid, "NTYPES", ntypes, &dimid_nb));
        int dims2[2] = {dimid_nb, dimid_nb};
        NC_CHECK(nc_def_var(ncid, "NONBONDED_PARM_INDEX", NC_INT, 2, dims2, &vid));
        const char *long_name="1-based LJ parm echelon index: max(i,j)*(max(i,j)+1)/2+min(i,j)+1";
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", strlen(long_name),long_name));
        NC_CHECK(nc_put_att_text(ncid, vid, "reference", 39,
                                 "LENNARD_JONES_ACOEF LENNARD_JONES_BCOEF"));
        const char *note="It is probably faster to do the math than to use a lookup table";
        NC_CHECK(nc_put_att_text(ncid, vid, "note", strlen(note), note));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, PVAI(vaNBIndexMatrix, int, 0)));
    }

    /* -9- residue labels */
    {
        char *buf = malloc(n_residues * LABLEN);
        char *buf_type = malloc(n_residues);
        int vid_type;
        int iUnknownTypes = 0;
        const char *type_info = "residue type: w=water, p=protein, n=nucleic, s=saccharide,"
                                " i=ion, l=ligand, o=other, ?=undefined";
        int dims2[2] = { dimid_res, dimid_name4};
        memset(buf, ' ', n_residues * LABLEN);
        SAVERESIDUEt *srPRes = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        for (int i = 0; i < n_residues; i++, srPRes++) {
            const char *s = srPRes->sName;
            size_t l = strlen(s);
            if (l > 3) s += l - 3;
            memcpy(buf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
            buf_type[i] = srPRes->sResidueType[0];
            if (buf_type[i] == 0) buf_type[i]='?';
            if (buf_type[i] == '?') iUnknownTypes++;
        }
        if (iUnknownTypes>0)
            VPWARN("%d residues have undefined type\n",iUnknownTypes);
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "RESIDUE_LABEL", NC_CHAR, 2, dims2, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 12, "residue name"));
        NC_CHECK(nc_put_att_text(ncid, vid, "pdb_field",     7, "resName"));
        if (GDefaults.bCIFReadAuth)
            NC_CHECK(nc_put_att_text(ncid, vid, "mmcif_field",    22, "atom_site.auth_comp_id"));
        else NC_CHECK(nc_put_att_text(ncid, vid, "mmcif_field",  23, "atom_site.label_comp_id"));
        NC_CHECK(nc_def_var(ncid, "RESIDUE_TYPE", NC_CHAR, 1, dims2, &vid_type));
        NC_CHECK(nc_put_att_text(ncid, vid_type, "long_name", strlen(type_info), type_info));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_text(ncid, vid, buf));
        NC_CHECK(nc_put_var_text(ncid, vid_type, buf_type));
        free(buf);
    }





    /* -10- residue pointer */
    {
        int *buf = malloc(n_residues * sizeof(int));
        for (int i = 0; i < n_residues; i++)
            buf[i] = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iAtomStartIndex;
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "RESIDUE_POINTER", NC_INT, 1, &dimid_res, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 36,
                                 "1-based first atom index per residue"));
        NC_CHECK(nc_put_att_text(ncid, vid, "reference_dimension", 5, "NATOM"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, buf));
        free(buf);
    }

    /* -11,-12A- bond force constants and equil (+ restraints) */
    {
        int nb = iParmSetTotalBondParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
        int n  = nb + nr;
        double *kb = malloc(n*sizeof(double)), *r0 = malloc(n*sizeof(double));
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_dim(ncid, "NUMBND", n, &dimid_bp));
        NC_CHECK(nc_def_var(ncid, "BOND_FORCE_CONSTANT", NC_DOUBLE, 1, &dimid_bp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 19, "kcal/mol/angstrom^2"));
        NC_CHECK(nc_enddef(ncid));
        for (int i = 0; i < nb; i++) {
            ParmSetBond(uUnit->psParameters, i, s1, s2, &dKb, &dR0,
                        &dKpull, &dRpull0, &dKpress, &dRpress0, sd);
            kb[i] = dKb; r0[i] = dR0;
        }
        NC_RESTRAINTLOOP(uUnit, RESTRAINTBOND, dKx, kb, nb);
        NC_CHECK(nc_put_var_double(ncid, vid, kb));

        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "BOND_EQUIL_VALUE", NC_DOUBLE, 1, &dimid_bp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 8, "angstrom"));
        NC_CHECK(nc_enddef(ncid));
        NC_RESTRAINTLOOP(uUnit, RESTRAINTBOND, dX0, r0, nb);
        NC_CHECK(nc_put_var_double(ncid, vid, r0));
        free(kb); free(r0);
    }

    /* -12B..E- bond augmentation, same NUMBND dimension as above */
    if (BondAugmentationFound(uUnit) == 1) {
        int nb = iParmSetTotalBondParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
        int n  = nb + nr;
        double *kpull   = malloc(n*sizeof(double));
        double *rpull0  = malloc(n*sizeof(double));
        double *kpress  = malloc(n*sizeof(double));
        double *rpress0 = malloc(n*sizeof(double));
        for (int i = 0; i < nb; i++) {
            ParmSetBond(uUnit->psParameters, i, s1, s2, &dKb, &dR0,
                        &dKpull, &dRpull0, &dKpress, &dRpress0, sd);
            kpull[i]=dKpull; rpull0[i]=dRpull0;
            kpress[i]=dKpress; rpress0[i]=dRpress0;
        }
        for (int i = nb; i < n; i++) {
            kpull[i]=0.0; rpull0[i]=100.0;
            kpress[i]=0.0; rpress0[i]=-100.0;
        }
        NC_CHECK(nc_def_var(ncid, "BOND_STIFFNESS_PULL_ADJ", NC_DOUBLE, 1, &dimid_bp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 19, "kcal/mol/angstrom^2"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_double(ncid, vid, kpull));

        NC_CHECK(nc_def_var(ncid, "BOND_EQUIL_PULL_ADJ", NC_DOUBLE, 1, &dimid_bp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 8, "angstrom"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_double(ncid, vid, rpull0));

        NC_CHECK(nc_def_var(ncid, "BOND_STIFFNESS_PRESS_ADJ", NC_DOUBLE, 1, &dimid_bp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 19, "kcal/mol/angstrom^2"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_double(ncid, vid, kpress));

        NC_CHECK(nc_def_var(ncid, "BOND_EQUIL_PRESS_ADJ", NC_DOUBLE, 1, &dimid_bp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 8, "angstrom"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_double(ncid, vid, rpress0));
        free(kpull); free(rpull0); free(kpress); free(rpress0);
    }

    /* -13,-14- angle force constants and equil (+ restraints) */
    {
        int na = iParmSetTotalAngleParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
        int n  = na + nr;
        int vid_force, vid_equil;
        double *kt = malloc(n*sizeof(double)), *t0 = malloc(n*sizeof(double));
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_dim(ncid, "NUMANG", n, &dimid_ap));
        NC_CHECK(nc_def_var(ncid, "ANGLE_FORCE_CONSTANT", NC_DOUBLE, 1, &dimid_ap, &vid_force));
        NC_CHECK(nc_put_att_text(ncid, vid_force, "units", 17, "kcal/mol/radian^2"));
        NC_CHECK(nc_def_var(ncid, "ANGLE_EQUIL_VALUE", NC_DOUBLE, 1, &dimid_ap, &vid_equil));
        NC_CHECK(nc_put_att_text(ncid, vid_equil, "units", 6, "degrees"));
        NC_CHECK(nc_enddef(ncid));

        for (int i = 0; i < na; i++) {
            ParmSetAngle(uUnit->psParameters, i, s1, s2, s3, &dKt, &dT0,
                         &dTkub, &dRkub, sd);
            kt[i] = dKt; t0[i] = dT0 / DEGTORAD;
        }
        NC_RESTRAINTLOOP(uUnit, RESTRAINTANGLE, dKx, kt, na);
        NC_CHECK(nc_put_var_double(ncid, vid_force, kt));

        NC_RESTRAINTLOOP(uUnit, RESTRAINTANGLE, dX0, t0, na);
        NC_CHECK(nc_put_var_double(ncid, vid_equil, t0));
        free(kt); free(t0);
    }

    /* -15,-16,-17,-17B,-17C- torsion/improper params (+ restraints) */
    {
        int iN, dimid_tp;
        int nt  = iParmSetTotalTorsionParms(uUnit->psParameters);
        int ni  = iParmSetTotalImproperParms(uUnit->psParameters);
        int nr  = iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION);
        int n   = nt + ni + nr;
        double *kp   = malloc(n*sizeof(double));
        int *per  = malloc(n*sizeof(int));
        double *p0   = malloc(n*sizeof(double));
        double *scee = malloc(n*sizeof(double));
        double *scnb = malloc(n*sizeof(double));
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_dim(ncid, "NPTRA", n, &dimid_tp));
        NC_CHECK(nc_def_var(ncid, "DIHEDRAL_FORCE_CONSTANT", NC_DOUBLE, 1, &dimid_tp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 8, "kcal/mol"));
        NC_CHECK(nc_enddef(ncid));
        for (int i = 0; i < nt; i++) {
            ParmSetTorsion(uUnit->psParameters, i, s1, s2, s3, s4,
                           &iN, &dKp, &dP0, &dScEE, &dScNB, sd);
            kp[i]=dKp; per[i]=iN; p0[i]=dP0/DEGTORAD;
            scee[i]=(dScEE<0.0)?GDefaults.dSceeScaleFactor:dScEE;
            scnb[i]=(dScNB<0.0)?GDefaults.dScnbScaleFactor:dScNB;
        }
        for (int i = 0; i < ni; i++) {
            ParmSetImproper(uUnit->psParameters, i, s1, s2, s3, s4,
                            &iN, &dKp, &dP0, sd);
            kp[nt+i]=dKp; per[nt+i]=iN; p0[nt+i]=dP0/DEGTORAD;
            scee[nt+i]=0.0; scnb[nt+i]=0.0;
        }
        NC_RESTRAINTLOOP(uUnit, RESTRAINTTORSION, dKx, kp,  nt+ni);
        NC_RESTRAINTLOOP(uUnit, RESTRAINTTORSION, dX0, p0,  nt+ni);
        NC_RESTRAINTLOOP(uUnit, RESTRAINTTORSION, dN,  per, nt+ni);
        NC_CHECK(nc_put_var_double(ncid, vid, kp));

        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "DIHEDRAL_PERIODICITY", NC_INT, 1, &dimid_tp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 47,
                "torsion Fourier periodicity n: cos(n*phi-phase)"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, per));

        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "DIHEDRAL_PHASE", NC_DOUBLE, 1, &dimid_tp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 6, "degree"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_double(ncid, vid, p0));

        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "SCEE_SCALE_FACTOR", NC_DOUBLE, 1, &dimid_tp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 37, "1-4 electrostatic energy scale factor"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_double(ncid, vid, scee));

        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "SCNB_SCALE_FACTOR", NC_DOUBLE, 1, &dimid_tp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 27, "1-4 VdW energy scale factor"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_double(ncid, vid, scnb));
        free(kp); free(per); free(p0); free(scee); free(scnb);
    }

    int dimid_nttyp;
    /* -19,-20- LJ coefficients */
    {
        int vid_acoef, vid_bcoef;
        double *A = malloc(nttyp*sizeof(double)), *B = malloc(nttyp*sizeof(double));
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_dim(ncid, "NTTYP", nttyp, &dimid_nttyp));
        NC_CHECK(nc_def_var(ncid, "LENNARD_JONES_ACOEF", NC_DOUBLE, 1, &dimid_nttyp, &vid_acoef));
        NC_CHECK(nc_put_att_text(ncid, vid_acoef, "units", 20, "kcal/mol*angstrom^12"));
        NC_CHECK(nc_put_att_text(ncid, vid_acoef, "long_name", 48,
                                 "LJ parm A coefficients, size NTYPES*(NTYPES+1)/2"));
        NC_CHECK(nc_def_var(ncid, "LENNARD_JONES_BCOEF", NC_DOUBLE, 1, &dimid_nttyp, &vid_bcoef));
        NC_CHECK(nc_put_att_text(ncid, vid_bcoef, "units", 19, "kcal/mol*angstrom^6"));
        NC_CHECK(nc_put_att_text(ncid, vid_acoef, "long_name", 48,
                                 "LJ parm B coefficients, size NTYPES*(NTYPES+1)/2"));
        NC_CHECK(nc_enddef(ncid));
        for (int i = 0; i < nttyp; i++) {
            A[i] = PVAI(vaNBParameters, NONBONDACt, i)->dA;
            B[i] = PVAI(vaNBParameters, NONBONDACt, i)->dB;
        }
        NC_CHECK(nc_put_var_double(ncid, vid_acoef, A));
        NC_CHECK(nc_put_var_double(ncid, vid_bcoef, B));
        free(A); free(B);
    }

    /* -21,-22- bonds inc/excl hydrogen */
    {
        int nb = iVarArrayElementCount(uUnit->vaBonds);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
        int (*ah)[2]  = malloc(2*nb*sizeof(int));
        int *pih      = malloc(nb*sizeof(int));
        int (*anh)[2] = malloc(2*(nb+nr)*sizeof(int));
        int *pinh     = malloc((nb+nr)*sizeof(int));
        int nh = 0, nnh = 0;
        for (int i = 0; i < nb; i++) {
            SAVEBONDt *sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, i);
            ATOMt *aA = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1-1)->aAtom;
            ATOMt *aD = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2-1)->aAtom;
            if (bPERT_BOND(bPert, aA, aD)) continue;
            if (iAtomElement(aA)==HYDROGEN || iAtomElement(aD)==HYDROGEN) {
                ah[nh][0]=sbPBond->iAtom1; ah[nh][1]=sbPBond->iAtom2;
                pih[nh]=sbPBond->iParmIndex; nh++;
            } else {
                anh[nnh][0]=sbPBond->iAtom1; anh[nnh][1]=sbPBond->iAtom2;
                pinh[nnh]=sbPBond->iParmIndex; nnh++;
            }
        }
        if (iVarArrayElementCount(uUnit->vaRestraints)) {
            SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
            for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++)
                if (sr->iType == RESTRAINTBOND) {
                    anh[nnh][0]=sr->iAtom1; anh[nnh][1]=sr->iAtom2;
                    pinh[nnh]=sr->iParmIndex; nnh++;
                }
        }
        int dimid_bh, dimid_bnh;
        int vid_bh_atoms, vid_bh_parm, vid_bnh_atoms, vid_bnh_parm;
        int dims2[2];

        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_dim(ncid, "NBONH", nh,  &dimid_bh));
        NC_CHECK(nc_def_dim(ncid, "MBONA", nnh, &dimid_bnh));

        dims2[0]=dimid_bh; dims2[1]=dimid_atom2;
        NC_CHECK(nc_def_var(ncid, "BONDS_INC_HYDROGEN_ATOMS", NC_INT, 2, dims2, &vid_bh_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_bh_atoms, "long_name",          20, "1-based atom indices"));
        NC_CHECK(nc_put_att_text(ncid, vid_bh_atoms, "reference_dimension", 5, "NATOM"));

        NC_CHECK(nc_def_var(ncid, "BONDS_INC_HYDROGEN_PARM", NC_INT, 1, &dimid_bh, &vid_bh_parm));
        NC_CHECK(nc_put_att_text(ncid, vid_bh_parm, "long_name", 18, "1-based parm index"));
        NC_CHECK(nc_put_att_text(ncid, vid_bh_parm, "reference_dimenson", 6, "NUMBND"));
        NC_CHECK(nc_put_att_text(ncid, vid_bh_parm, "reference", 36,
                                 "BOND_FORCE_CONSTANT BOND_EQUIL_VALUE"));

        dims2[0]=dimid_bnh; dims2[1]=dimid_atom2;
        NC_CHECK(nc_def_var(ncid, "BONDS_WITHOUT_HYDROGEN_ATOMS", NC_INT, 2, dims2, &vid_bnh_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_bnh_atoms, "long_name",          20, "1-based atom indices"));
        NC_CHECK(nc_put_att_text(ncid, vid_bnh_atoms, "reference_dimension", 5, "NATOM"));

        NC_CHECK(nc_def_var(ncid, "BONDS_WITHOUT_HYDROGEN_PARM", NC_INT, 1, &dimid_bnh, &vid_bnh_parm));
        NC_CHECK(nc_put_att_text(ncid, vid_bnh_parm, "long_name", 18, "1-based parm index"));
        NC_CHECK(nc_put_att_text(ncid, vid_bnh_parm, "reference_dimenson", 6, "NUMBND"));
        NC_CHECK(nc_put_att_text(ncid, vid_bnh_parm, "reference", 36,
                                 "BOND_FORCE_CONSTANT BOND_EQUIL_VALUE"));
        NC_CHECK(nc_enddef(ncid));

        NC_CHECK(nc_put_var_int(ncid, vid_bh_atoms,  (int*)ah));
        NC_CHECK(nc_put_var_int(ncid, vid_bh_parm,   pih));
        NC_CHECK(nc_put_var_int(ncid, vid_bnh_atoms, (int*)anh));
        NC_CHECK(nc_put_var_int(ncid, vid_bnh_parm,  pinh));
        free(ah); free(pih);
        free(anh); free(pinh);
    }

    /* -23,-24- angles inc/excl hydrogen */
    {
        int na = iVarArrayElementCount(uUnit->vaAngles);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
        int (*ah)[3]  = malloc(3*na*sizeof(int));
        int (*anh)[3] = malloc(3*(na+nr)*sizeof(int));
        int *pih      = malloc(na*sizeof(int));
        int *pinh     = malloc((na+nr)*sizeof(int));
        int nh = 0, nnh = 0;
        for (int i = 0; i < na; i++) {
            SAVEANGLEt *saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
            ATOMt *aA = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1-1)->aAtom;
            ATOMt *aB = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2-1)->aAtom;
            ATOMt *aD = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3-1)->aAtom;
            if (bPERT_ANGLE(bPert, aA, aB, aD)) continue;
            if (iAtomElement(aA)==HYDROGEN||iAtomElement(aB)==HYDROGEN||
                iAtomElement(aD)==HYDROGEN) {
                ah[nh][0]=saPAngle->iAtom1; ah[nh][1]=saPAngle->iAtom2;
                ah[nh][2]=saPAngle->iAtom3; pih[nh]=saPAngle->iParmIndex; nh++;
            } else {
                anh[nnh][0]=saPAngle->iAtom1; anh[nnh][1]=saPAngle->iAtom2;
                anh[nnh][2]=saPAngle->iAtom3; pinh[nnh]=saPAngle->iParmIndex; nnh++;
            }
        }
        if (iVarArrayElementCount(uUnit->vaRestraints)) {
            SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
            for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++)
                if (sr->iType == RESTRAINTANGLE) {
                    anh[nnh][0]=sr->iAtom1; anh[nnh][1]=sr->iAtom2;
                    anh[nnh][2]=sr->iAtom3; pinh[nnh]=sr->iParmIndex; nnh++;
                }
        }
        int dimid_ah, dimid_anh;
        int vid_ah_atoms, vid_ah_parm, vid_anh_atoms, vid_anh_parm;
        int dims2[2];

        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_dim(ncid, "NTHETH", nh,  &dimid_ah));
        NC_CHECK(nc_def_dim(ncid, "MTHETA", nnh, &dimid_anh));

        dims2[0]=dimid_ah; dims2[1]=dimid_atom3;
        NC_CHECK(nc_def_var(ncid, "ANGLES_INC_HYDROGEN_ATOMS", NC_INT, 2, dims2, &vid_ah_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_ah_atoms, "long_name",          20, "1-based atom indices"));
        NC_CHECK(nc_put_att_text(ncid, vid_ah_atoms, "reference_dimension", 5, "NATOM"));

        NC_CHECK(nc_def_var(ncid, "ANGLES_INC_HYDROGEN_PARM", NC_INT, 1, &dimid_ah, &vid_ah_parm));
        NC_CHECK(nc_put_att_text(ncid, vid_ah_parm, "long_name", 18, "1-based parm index"));
        NC_CHECK(nc_put_att_text(ncid, vid_ah_parm, "reference_dimenson", 6, "NUMANG"));
        NC_CHECK(nc_put_att_text(ncid, vid_ah_parm, "reference", 38,
                                 "ANGLE_FORCE_CONSTANT ANGLE_EQUIL_VALUE"));

        dims2[0]=dimid_anh; dims2[1]=dimid_atom3;
        NC_CHECK(nc_def_var(ncid, "ANGLES_WITHOUT_HYDROGEN_ATOMS", NC_INT, 2, dims2, &vid_anh_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_anh_atoms, "long_name",          20, "1-based atom indices"));
        NC_CHECK(nc_put_att_text(ncid, vid_anh_atoms, "reference_dimension", 5, "NATOM"));

        NC_CHECK(nc_def_var(ncid, "ANGLES_WITHOUT_HYDROGEN_PARM", NC_INT, 1, &dimid_anh, &vid_anh_parm));
        NC_CHECK(nc_put_att_text(ncid, vid_anh_parm, "long_name", 18, "1-based parm index"));
        NC_CHECK(nc_put_att_text(ncid, vid_anh_parm, "reference_dimenson", 6, "NUMANG"));
        NC_CHECK(nc_put_att_text(ncid, vid_anh_parm, "reference", 38,
                                 "ANGLE_FORCE_CONSTANT ANGLE_EQUIL_VALUE"));
        NC_CHECK(nc_enddef(ncid));

        NC_CHECK(nc_put_var_int(ncid, vid_ah_atoms,  (int*)ah));
        NC_CHECK(nc_put_var_int(ncid, vid_ah_parm,   pih));
        NC_CHECK(nc_put_var_int(ncid, vid_anh_atoms, (int*)anh));
        NC_CHECK(nc_put_var_int(ncid, vid_anh_parm,  pinh));
        free(ah); free(pih);
        free(anh); free(pinh);
    }

    /* -25,-26- torsions */
    NC_CHECK(write_torsions(ncid, uUnit, bPert, dimid_atom4));

    /* -27- excluded atoms list */
    {
        int n = iVarArrayElementCount(vaExcludedAtoms);
        int dimid_excl;
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_dim(ncid, "NEXT", n, &dimid_excl));
        NC_CHECK(nc_def_var(ncid, "EXCLUDED_ATOMS_LIST", NC_INT, 1, &dimid_excl, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name",          19, "excluded atom list"));
        NC_CHECK(nc_put_att_text(ncid, vid, "reference_dimension", 5, "NATOM"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, PVAI(vaExcludedAtoms, int, 0)));
    }

    /* -28,-29,-30- hbond (skip if none) */
    {
        int n = iParmSetTotalHBondParms(uUnit->psParameters);
        if (n > 0) {
            int dimid_hb, vid_hba, vid_hbb;
            double *A   = malloc(n*sizeof(double));
            double *B   = malloc(n*sizeof(double));

            NC_CHECK(nc_redef(ncid));
            NC_CHECK(nc_def_dim(ncid, "NPHB", n, &dimid_hb));
            NC_CHECK(nc_def_var(ncid, "HBOND_ACOEF", NC_DOUBLE, 1, &dimid_hb, &vid_hba));
            NC_CHECK(nc_put_att_text(ncid, vid_hba, "units", 20, "kcal/mol*angstrom^12"));

            NC_CHECK(nc_def_var(ncid, "HBOND_BCOEF", NC_DOUBLE, 1, &dimid_hb, &vid_hbb));
            NC_CHECK(nc_put_att_text(ncid, vid_hbb, "units", 20, "kcal/mol*angstrom^10"));

            //NC_CHECK(nc_def_var(ncid, "HBCUT", NC_DOUBLE, 1, &dimid_hb, &vid_hbc));
            //NC_CHECK(nc_put_att_text(ncid, vid_hbc, "long_name", 21, "reserved, always zero"));

            NC_CHECK(nc_enddef(ncid));

            for (int i = 0; i < n; i++) {
                ParmSetHBond(uUnit->psParameters, i, s1, s2, &dC, &dD, sd);
                A[i] = dC; B[i] = dD;
            }
            NC_CHECK(nc_put_var_double(ncid, vid_hba, A));
            NC_CHECK(nc_put_var_double(ncid, vid_hbb, B));
            //NC_CHECK(nc_put_var_double(ncid, vid, cut));
            free(A); free(B);
        }
    }

    /* -31- amber atom type strings */
    {
        char *buf = malloc(n_atoms * LABLEN);
        memset(buf, ' ', n_atoms * LABLEN);
        for (int i = 0; i < n_atoms; i++) {
            const char *s = sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom);
            size_t l = strlen(s);
            if (l > 4) s += l - 4;
            memcpy(buf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
        }
        int dims2[2] = {dimid_atoms, dimid_name4};
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "AMBER_ATOM_TYPE", NC_CHAR, 2, dims2, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 26, "force field atom type name"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_text(ncid, vid, buf));
        free(buf);
    }

    /* -32- tree chain classification */
    {
        char *buf = malloc(n_atoms * LABLEN);
        for (int i = 0; i < n_atoms; i++) {
            ATOMt *aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom;
            const char *tc = "BLA ";
            double dT = dAtomTemp(aAtom);
            if      (dT==(double)'M') tc="M   ";
            else if (dT==(double)'E') tc="E   ";
            else if (dT==(double)'S') tc="S   ";
            else if (dT==(double)'B') tc="B   ";
            else if (dT==(double)'3') tc="3   ";
            else if (dT==(double)'4') tc="4   ";
            else if (dT==(double)'5') tc="5   ";
            else if (dT==(double)'6') tc="6   ";
            else if (dT==(double)'X') tc="X   ";
            memcpy(buf + i*LABLEN, tc, LABLEN);
        }
        int dims2[2] = { dimid_atoms, dimid_name4};
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "TREE_CHAIN_CLASSIFICATION", NC_CHAR, 2, dims2, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 52, "mainchain tree class: M,E,S,B,3,4,5,6,X; BLA=unknown"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_text(ncid, vid, buf));
        free(buf);
    }

    /* -35A..C- solvent/box info (conditional) */
    if (bUnitUseBox(uUnit)) {
        int iFirstSolvRes = n_residues;
        for (int i = 0; i < n_residues; i++) {
            if (PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sResidueType[0] == RESTYPESOLVENT) {
                iFirstSolvRes = i; break;
            }
        }
        UnitIOFindAndCountMolecules(uUnit);
        int nMol = iVarArrayElementCount(uUnit->vaAtomsPerMolecule);
        int iFS1 = uUnit->iFirstSolvent + 1;

        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "IPTRES", NC_INT, 0, NULL, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name",          46,
                                 "1-based residue index of first solvent residue"));
        NC_CHECK(nc_put_att_text(ncid, vid, "reference_dimension", 4, "NRES"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, &iFirstSolvRes));

        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "NSPSOL", NC_INT, 0, NULL, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name",          48,
                                 "1-based molecule index of first solvent molecule"));
        NC_CHECK(nc_put_att_text(ncid, vid, "reference_dimension", 4, "NSPM"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, &iFS1));

        {
            int dimid_mol;
            NC_CHECK(nc_redef(ncid));
            NC_CHECK(nc_def_dim(ncid, "NSPM", nMol, &dimid_mol));
            NC_CHECK(nc_def_var(ncid, "ATOMS_PER_MOLECULE", NC_INT, 1, &dimid_mol, &vid));
            NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 30, "number of atoms for molecule N"));
            NC_CHECK(nc_enddef(ncid));
            NC_CHECK(nc_put_var_int(ncid, vid,
                                    PVAI(uUnit->vaAtomsPerMolecule, int, 0)));
        }
    }

    /* -35D,-35E- cap info (conditional) */
    if (bUnitUseSolventCap(uUnit)) {
        int ires = uUnit->iCapTempInt;
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "CAP_INFO_ATOM", NC_INT, 0, NULL, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name",          45,
                             "Index of the last atom not in the Solvent Cap"));
        NC_CHECK(nc_put_att_text(ncid, vid, "reference_dimension", 5, "NATOM"));
        NC_CHECK(nc_enddef(ncid));

        double dX, dY, dZ, dR;
        UnitGetSolventCap(uUnit, &dX, &dY, &dZ, &dR);
        double cap[3] = { dX, dY, dZ };
        int dimid_vec, vid_r, vid_xyz;
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "CAP_INFO_RADIUS", NC_DOUBLE, 0, NULL, &vid_r));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name",          21,
                             "Radius of Solvent Cap"));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 8, "angstrom"));
        NC_CHECK(nc_def_dim(ncid, "vector", 3, &dimid_vec));
        NC_CHECK(nc_def_var(ncid, "CAP_INFO_ORIGIN", NC_DOUBLE, 1, &dimid_vec, &vid_xyz));
        NC_CHECK(nc_put_att_text(ncid, vid_xyz, "units", 8, "angstrom"));
        NC_CHECK(nc_put_att_text(ncid, vid_r, "long_name", 21,
                             "Origin of Solvent Cap"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, &ires));
        NC_CHECK(nc_put_var_double(ncid, vid_r, &dR));
        NC_CHECK(nc_put_var_double(ncid, vid_xyz, cap));
    }

    /* GB radius set - global attribute */
    {
        STRING sTemp;
        const char *rset="Unknown radius set";
        if (GDefaults.iGBparm >= 0 && GDefaults.iGBparm <= 8) {
            sprintf(sTemp,"%s (%s)",PBRadii_optionDesc[GDefaults.iGBparm],
                      PBRadii_options[GDefaults.iGBparm]);
            rset = sTemp;
        }
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "RADIUS_SET",
                                 strlen(rset), rset));
        NC_CHECK(nc_enddef(ncid));
    }

    /* IPOL - true scalar */
    {
        int val = GDefaults.iIPOL;
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "IPOL", NC_INT, 0, NULL, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 28, "Polarizable force field flag"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, &val));
    }

    /* -36A..M- perturbation (conditional) */
    if (bPert) {
        /* -36A- perturbed bond atoms */
        {
            int nb = iVarArrayElementCount(uUnit->vaBonds);
            int (*a)[2] = malloc(2*nb*sizeof(int));
            int n = 0;
            for (int pass = 0; pass < 2; pass++) {
                SAVEBONDt *sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
                for (int i = 0; i < nb; i++, sbPBond++) {
                    bool bnd = (sbPBond->fFlags &BOUNDARY) != 0;
                    if ((sbPBond->fFlags &PERTURBED) && (pass==0 ? !bnd : bnd))
                        { a[n][0]=sbPBond->iAtom1; a[n][1]=sbPBond->iAtom2; n++; }
                }
            }
            if (n > 0) {
                int dimid_pb;
                int dims2[2];
                NC_CHECK(nc_redef(ncid));
                NC_CHECK(nc_def_dim(ncid, "NBPER", n, &dimid_pb));
                dims2[0]=dimid_pb; dims2[1]=dimid_atom2;
                NC_CHECK(nc_def_var(ncid, "PERT_BOND_ATOMS", NC_INT, 2, dims2, &vid));
                NC_CHECK(nc_put_att_text(ncid, vid, "reference_dimension", 5, "NATOM"));
                NC_CHECK(nc_enddef(ncid));
                NC_CHECK(nc_put_var_int(ncid, vid, (int*)a));
            }
            free(a);
        }

        /* -36B- perturbed bond params */
        {
            int nb = iVarArrayElementCount(uUnit->vaBonds);
            int *p0 = malloc(nb*sizeof(int)), *p1 = malloc(nb*sizeof(int));
            int n = 0;
            for (int pass = 0; pass < 2; pass++) {
                SAVEBONDt *sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
                for (int i = 0; i < nb; i++, sbPBond++) {
                    bool bnd = (sbPBond->fFlags &BOUNDARY) != 0;
                    if ((sbPBond->fFlags &PERTURBED) && (pass==0 ? !bnd : bnd)) {
                        p0[n]=sbPBond->iParmIndex;
                        p1[n]=sbPBond->iPertParmIndex;
                        n++;
                    }
                }
            }
            if (n > 0) {
                int dimid_pb;
                NC_CHECK(nc_redef(ncid));
                NC_CHECK(nc_def_dim(ncid, "MBPER", n, &dimid_pb));
                NC_CHECK(nc_def_var(ncid, "PERT_BOND_PARAMS_L0", NC_INT, 1, &dimid_pb, &vid));
                NC_CHECK(nc_put_att_text(ncid, vid, "reference_dimenson", 6, "NUMBND"));
                NC_CHECK(nc_put_att_text(ncid, vid, "reference", 36,
                                         "BOND_FORCE_CONSTANT BOND_EQUIL_VALUE"));
                NC_CHECK(nc_enddef(ncid));
                NC_CHECK(nc_put_var_int(ncid, vid, p0));
                NC_CHECK(nc_def_var(ncid, "PERT_BOND_PARAMS_L1", NC_INT, 1, &dimid_pb, &vid));
                NC_CHECK(nc_put_att_text(ncid, vid, "reference_dimenson", 6, "NUMBND"));
                NC_CHECK(nc_put_att_text(ncid, vid, "reference", 36,
                                         "BOND_FORCE_CONSTANT BOND_EQUIL_VALUE"));
                NC_CHECK(nc_enddef(ncid));
                NC_CHECK(nc_put_var_int(ncid, vid, p1));
            }
            free(p0); free(p1);
        }

        /* -36G- perturbed residue names */
        {
            char *buf = malloc(n_residues * LABLEN);
            memset(buf, ' ', n_residues * LABLEN);
            for (int i = 0; i < n_residues; i++) {
                const char *s = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sName;
                size_t l = strlen(s);
                if (l > 3) s += l - 3;
                memcpy(buf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
            }
            int dims2[2] = { dimid_res, dimid_name4};
            NC_CHECK(nc_def_var(ncid, "PERT_RESIDUE_NAME", NC_CHAR, 2, dims2, &vid));
            NC_CHECK(nc_enddef(ncid));
            NC_CHECK(nc_put_var_text(ncid, vid, buf));
            free(buf);
        }

        /* -36H- perturbed atom names */
        {
            char *buf = malloc(n_atoms * LABLEN);
            memset(buf, ' ', n_atoms * LABLEN);
            for (int i = 0; i < n_atoms; i++) {
                const char *s = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sPertName;
                if (strlen(s) == 0)
                    s = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sName;
                size_t l = strlen(s);
                if (l > 4) s += l - 4;
                memcpy(buf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
            }
            int dims2[2] = { dimid_atoms, dimid_name4};
            NC_CHECK(nc_def_var(ncid, "PERT_ATOM_NAME", NC_CHAR, 2, dims2, &vid));
            NC_CHECK(nc_enddef(ncid));
            NC_CHECK(nc_put_var_text(ncid, vid, buf));
            free(buf);
        }

        /* -36I- perturbed atom symbols */
        {
            char *buf = malloc(n_atoms * LABLEN);
            memset(buf, ' ', n_atoms * LABLEN);
            for (int i = 0; i < n_atoms; i++) {
                const char *s = sAtomPertType(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom);
                if (strlen(s) == 0)
                    s = sAtomType(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom);
                size_t l = strlen(s);
                if (l > 4) s += l - 4;
                memcpy(buf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
            }
            int dims2[2] = { dimid_atoms, dimid_name4};
            NC_CHECK(nc_def_var(ncid, "PERT_ATOM_SYMBOL", NC_CHAR, 2, dims2, &vid));
            NC_CHECK(nc_enddef(ncid));
            NC_CHECK(nc_put_var_text(ncid, vid, buf));
            free(buf);
        }

        /* -36K- IAPER */
        {
            int *buf = malloc(n_atoms * sizeof(int));
            for (int i = 0; i < n_atoms; i++)
                buf[i] = bAtomPerturbed(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom) ? 1 : 0;
            NC_CHECK(nc_def_var(ncid, "IAPER", NC_INT, 1, &dimid_atoms, &vid));
            NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 22, "1 if atom is perturbed"));
            NC_CHECK(nc_enddef(ncid));
            NC_CHECK(nc_put_var_int(ncid, vid, buf));
            free(buf);
        }

        /* -36L- perturbed atom type index */
        {
            int *buf = malloc(n_atoms * sizeof(int));
            for (int i = 0; i < n_atoms; i++) {
                int iAtom = PVAI(uUnit->vaAtoms,SAVEATOMt,i)->iPertTypeIndex - 1;
                buf[i] = *PVAI(vaNBIndex, int, iAtom) + 1;
            }
            NC_CHECK(nc_def_var(ncid, "PERT_ATOM_TYPE_INDEX", NC_INT, 1, &dimid_atoms, &vid));
            NC_CHECK(nc_put_att_text(ncid, vid, "long_name",          18, "1-based type index"));
            NC_CHECK(nc_put_att_text(ncid, vid, "reference_dimension", 6, "NTYPES"));
            NC_CHECK(nc_enddef(ncid));
            NC_CHECK(nc_put_var_int(ncid, vid, buf));
            free(buf);
        }

        /* -36M- perturbed charges */
        {
            double *buf = malloc(n_atoms * sizeof(double));
            for (int i = 0; i < n_atoms; i++) {
                ATOMt *aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom;
                if (GDefaults.iGibbs)
                    buf[i] = bAtomPerturbed(aAtom) ?
                             dAtomPertCharge(aAtom) : dAtomCharge(aAtom);
                else
                    buf[i] = dAtomCharge(aAtom) + dAtomPertCharge(aAtom);
            }
            NC_CHECK(nc_def_var(ncid, "PERT_CHARGE", NC_DOUBLE, 1, &dimid_atoms, &vid));
            NC_CHECK(nc_put_att_text(ncid, vid, "units",    18, "elementary_charge"));
            NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 51,
                                 "multiply by 18.2223 for AMBER internal charge units"));
            NC_CHECK(nc_enddef(ncid));
            NC_CHECK(nc_put_var_double(ncid, vid, buf));
            free(buf);
        }
    }

    /* polar */
    if (bPolar)
        NC_CHECK(write_polar_section(ncid, uUnit, dimid_atoms, bPert));


    /* CHARMM 1-4 terms */
    if (GDefaults.iCharmm) {
        //int n = iVarArrayElementCount(vaNBParameters);
        int nubtypes = iParmSetTotalAngleParms(uUnit->psParameters);
        int dimid_nubtypes;
        double *A14  = malloc(nttyp*sizeof(double)), *B14 = malloc(nttyp*sizeof(double));
        double *tkub = malloc(nubtypes*sizeof(double)), *rkub = malloc(nubtypes*sizeof(double));
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "LENNARD_JONES_14_ACOEF", NC_DOUBLE, 1, &dimid_nttyp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 20, "kcal/mol*angstrom^12"));
        NC_CHECK(nc_enddef(ncid));
        for (int i = 0; i < nttyp; i++) {
            A14[i] = PVAI(vaNBParameters, NONBONDACt, i)->dA14;
            B14[i] = PVAI(vaNBParameters, NONBONDACt, i)->dB14;
        }
        NC_CHECK(nc_put_var_double(ncid, vid, A14));

        NC_CHECK(nc_def_var(ncid, "LENNARD_JONES_14_BCOEF", NC_DOUBLE, 1, &dimid_nttyp, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 19, "kcal/mol*angstrom^6"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_double(ncid, vid, B14));
        free(A14); free(B14);

        // CHARMM_UREY_BRADLEY_COUNT = nub, nubtypes
        //int nub, nubtypes, dimid_nub, dimid_nubtypes;
        //NC_CHECK(nc_def_dim(ncid, "CHARMM_NUB", nub, &dimid_nub));
        NC_CHECK(nc_def_dim(ncid, "CHARMM_NUBTYPES", nubtypes, &dimid_nubtypes));
        for (int i = 0; i < nubtypes; i++) {
            ParmSetAngle(uUnit->psParameters, i, s1, s2, s3, &dKt, &dT0, &dTkub, &dRkub, sd);
            tkub[i] = dTkub; rkub[i] = dRkub;
        }
        NC_CHECK(nc_def_var(ncid, "UREY_BRADLEY_FORCE_CONSTANT", NC_DOUBLE, 1, &dimid_nubtypes, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 19, "kcal/mol/angstrom^2"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_double(ncid, vid, tkub));

        NC_CHECK(nc_def_var(ncid, "UREY_BRADLEY_EQUIL_VALUE", NC_DOUBLE, 1, &dimid_nubtypes, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "units", 8, "angstrom"));
        NC_CHECK(nc_enddef(ncid));

        NC_CHECK(nc_put_var_double(ncid, vid, rkub));
        free(tkub); free(rkub);
    }

    /* PDB chain info */
    if (GDefaults.bPdbKeepChainId) {
        int *buf = malloc(n_residues * sizeof(int));
        for (int i = 0; i < n_residues; i++)
            buf[i] = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iPdbResSeq;
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "RESIDUE_NUMBER", NC_INT, 1, &dimid_res, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 10, "PDB resSeq"));
        NC_CHECK(nc_put_att_text(ncid, vid, "pdb_field",     6, "resSeq"));
        if (GDefaults.bCIFReadAuth)
            NC_CHECK(nc_put_att_text(ncid, vid, "mmcif_field",    21, "atom_site.auth_seq_id"));
        else NC_CHECK(nc_put_att_text(ncid, vid, "mmcif_field",  22, "atom_site.label_seq_id"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_int(ncid, vid, buf));
        free(buf);

        char *cbuf = malloc(n_residues * LABLEN);
        memset(cbuf, ' ', n_residues * LABLEN);
        for (int i = 0; i < n_residues; i++) {
            const char *s = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sChainId;
            size_t l = strlen(s);
            if (l > 2) s += l - 2;
            memcpy(cbuf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
        }
        int dims2[2] = { dimid_res, dimid_name4};
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "RESIDUE_CHAINID", NC_CHAR, 2, dims2, &vid));
        NC_CHECK(nc_put_att_text(ncid, vid, "long_name", 10, "PDB resSeq"));
        NC_CHECK(nc_put_att_text(ncid, vid, "pdb_field",     7, "chainID"));
        if (GDefaults.bCIFReadAuth)
            NC_CHECK(nc_put_att_text(ncid, vid, "mmcif_field",    22, "atom_site.auth_comp_id"));
        else NC_CHECK(nc_put_att_text(ncid, vid, "mmcif_field",  23, "atom_site.label_comp_id"));
        NC_CHECK(nc_enddef(ncid));
        NC_CHECK(nc_put_var_text(ncid, vid, cbuf));
        free(cbuf);
    }

    /* CMAP -- pending data structure review */
    NC_CHECK(SaveAmberParmCMAPNetcdf(uUnit, ncid));

    NC_CHECK(nc_close(ncid));
    printf("Successfully saved NetCDF PRMTOP file \"%s\"\n", fname);
    return NC_NOERR;
}

/* ================================================================
   UnitIOSaveAmberCoordNetcdf
   Author: Robin Betz (2013)
   Writes coordinates in UNIT to a NetCDF coordinate file.
   ================================================================ */
void
UnitIOSaveAmberCoordNetcdf(UNIT uUnit, char *filename)
{
    int ncid;
    int did_spatial, did_atom;
    int vid_spatial, vid_coord;
    int did_cell_spatial, did_cell_angular, did_label;
    int vid_cell_spatial, vid_cell_angular, vid_cell_length, vid_cell_angle,
        vid_time;
    int dimensionID[NC_MAX_VAR_DIMS];

    int iAtomCount = iVarArrayElementCount(uUnit->vaAtoms);
    printf("There are %i atoms\n", iAtomCount);

    if (nc_create(filename, NC_64BIT_OFFSET, &ncid) != NC_NOERR)
        VPFATALEXIT("%s: Error creating file\n", filename);

    if (nc_def_dim(ncid, "spatial", 3, &did_spatial) != NC_NOERR)
        VPFATALEXIT("%s: Error defining spatial dimension\n", filename);

    dimensionID[0] = did_spatial;
    if (nc_def_var(ncid, "spatial", NC_CHAR, 1, dimensionID, &vid_spatial) != NC_NOERR)
        VPFATALEXIT("%s: Error defining spatial variable\n", filename);

    if (nc_def_dim(ncid, "atom", iAtomCount, &did_atom) != NC_NOERR)
        VPFATALEXIT("%s: Error defining atom dimension\n", filename);

    if (nc_def_var(ncid, "time", NC_DOUBLE, 0, dimensionID, &vid_time) != NC_NOERR)
        VPFATALEXIT("%s: Error defining time variable\n", filename);
    if (nc_put_att_text(ncid, vid_time, "units", 10, "picosecond") != NC_NOERR)
        VPFATALEXIT("%s: Error setting time units to picosecond\n", filename);

    dimensionID[0] = did_atom;
    dimensionID[1] = did_spatial;
    if (nc_def_var(ncid, "coordinates", NC_DOUBLE, 2, dimensionID, &vid_coord) != NC_NOERR)
        VPFATALEXIT("%s: Error defining coordinate variable\n", filename);
    if (nc_put_att_text(ncid, vid_coord, "units", 8, "angstrom") != NC_NOERR)
        VPFATALEXIT("%s: Error setting coordinate units to angstrom\n", filename);

    if (bUnitUseBox(uUnit) == TRUE) {
        printf("Using the unit box\n");
        if (nc_def_dim(ncid, "cell_spatial", 3, &did_cell_spatial) != NC_NOERR)
            VPFATALEXIT("%s: Error defining cell spatial dimension\n", filename);
        dimensionID[0] = did_cell_spatial;
        if (nc_def_var(ncid, "cell_spatial", NC_CHAR, 1, dimensionID, &vid_cell_spatial) != NC_NOERR)
            VPFATALEXIT("%s: Error defining cell spatial variable\n", filename);

        if (nc_def_dim(ncid, "label", 5, &did_label) != NC_NOERR)
            VPFATALEXIT("%s: Error defining label dimension\n", filename);
        if (nc_def_dim(ncid, "cell_angular", 3, &did_cell_angular) != NC_NOERR)
            VPFATALEXIT("%s: Error defining cell angular dimension\n", filename);
        dimensionID[0] = did_cell_angular;
        dimensionID[1] = did_label;
        if (nc_def_var(ncid, "cell_angular", NC_CHAR, 2, dimensionID, &vid_cell_angular) != NC_NOERR)
            VPFATALEXIT("%s: Error defining cell angular variable\n", filename);

        dimensionID[0] = did_cell_spatial;
        if (nc_def_var(ncid, "cell_lengths", NC_DOUBLE, 1, dimensionID, &vid_cell_length) != NC_NOERR)
            VPFATALEXIT("%s: Error defining cell lengths\n", filename);
        if (nc_put_att_text(ncid, vid_cell_length, "units", 8, "angstrom") != NC_NOERR)
            VPFATALEXIT("%s: Error setting cell length units to angstrom\n", filename);

        dimensionID[0] = did_cell_angular;
        if (nc_def_var(ncid, "cell_angles", NC_DOUBLE, 1, dimensionID, &vid_cell_angle) != NC_NOERR)
            VPFATALEXIT("%s: Error defining cell angles variable\n", filename);
        if (nc_put_att_text(ncid, vid_cell_angle, "units", 6, "degree") != NC_NOERR)
            VPFATALEXIT("%s: Error setting cell angle units to degree\n", filename);
    }

    if (nc_put_att_text(ncid, NC_GLOBAL, "title",
                        strlen(sContainerName(uUnit)), sContainerName(uUnit)) != NC_NOERR)
        VPFATALEXIT("%s: Error writing title\n", filename);
    if (nc_put_att_text(ncid, NC_GLOBAL, "application",      5, "AMBER") != NC_NOERR)
        VPFATALEXIT("%s: Error writing application string\n", filename);
    if (nc_put_att_text(ncid, NC_GLOBAL, "program",          4, "leap") != NC_NOERR)
        VPFATALEXIT("%s: Error writing program string\n", filename);
    if (nc_put_att_text(ncid, NC_GLOBAL, "programVersion",   3, "1.0") != NC_NOERR)
        VPFATALEXIT("%s: Error writing program version string\n", filename);
    if (nc_put_att_text(ncid, NC_GLOBAL, "Conventions",     12, "AMBERRESTART") != NC_NOERR)
        VPFATALEXIT("%s: Error writing conventions\n", filename);
    if (nc_put_att_text(ncid, NC_GLOBAL, "ConventionVersion", 3, "1.0") != NC_NOERR)
        VPFATALEXIT("%s: Error writing conventions version\n", filename);

    if (nc_set_fill(ncid, NC_NOFILL, dimensionID) != NC_NOERR)
        VPFATALEXIT("%s: Error setting fill mode\n", filename);
    if (nc_enddef(ncid) != NC_NOERR)
        VPFATALEXIT("%s: NetCDF error on ending definitions\n", filename);

    size_t start[2]={0}, count[2]={0};
    count[0] = 3;
    if (nc_put_vara_text(ncid, vid_spatial, start, count, "xyz") != NC_NOERR)
        VPFATALEXIT("%s: Error writing spatial labels\n", filename);

    if (bUnitUseBox(uUnit) == TRUE) {
        count[1] = 1;
        if (nc_put_vara_text(ncid, vid_cell_spatial, start, count, "abc") != NC_NOERR)
            VPFATALEXIT("%s: Error writing cell spatial labels\n", filename);
        count[1] = 5;
        if (nc_put_vara_text(ncid, vid_cell_angular, start, count,
                             "alpha" "beta " "gamma") != NC_NOERR)
            VPFATALEXIT("%s: Error writing cell angular labels\n", filename);
    }

    double *data;
    MALLOC(data, double *, 3 * iAtomCount * sizeof(double));
    int counter = 0;
    VECTOR vPos;

    double time = 0.0;
    if (nc_put_var_double(ncid, vid_time, &time) != NC_NOERR)
        VPFATALEXIT("%s: Error writing start time\n", filename);

    double dX, dY, dZ, dX2, dY2, dZ2;
    if (bUnitUseBox(uUnit) == TRUE && GDefaults.nocenter == 0) {
        UnitGetBox(uUnit, &dX, &dY, &dZ);
        dX2 = dX * 0.5; dY2 = dY * 0.5; dZ2 = dZ * 0.5;
    } else {
        dX2 = dY2 = dZ2 = 0.0;
    }

    for (int i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        vPos = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->vPos;
        data[counter++] = dVX(&vPos) + dX2;
        data[counter++] = dVY(&vPos) + dY2;
        data[counter++] = dVZ(&vPos) + dZ2;
    }

    start[0] = start[1] = 0;
    count[0] = iAtomCount; count[1] = 3;
    if (nc_put_vara_double(ncid, vid_coord, start, count, data) != NC_NOERR) {
        FREE(data);
        VPFATALEXIT("%s: Error writing coordinate data\n", filename);
    }

    if (bUnitUseBox(uUnit) == TRUE) {
        count[0] = 3; count[1] = 0;
        double lengths[3] = { dX, dY, dZ };
        if (nc_put_vara_double(ncid, vid_cell_length, start, count, lengths) != NC_NOERR) {
            FREE(data);
            VPFATALEXIT("%s: Error writing cell lengths\n", filename);
        }
        lengths[0] = lengths[1] = lengths[2] = dUnitBeta(uUnit) / DEGTORAD;
        if (nc_put_vara_double(ncid, vid_cell_angle, start, count, lengths) != NC_NOERR) {
            FREE(data);
            VPFATALEXIT("%s: Error writing cell angles\n", filename);
        }
    }

    if (nc_close(ncid) != NC_NOERR) {
        FREE(data);
        VPFATALEXIT("%s: Error closing file\n", filename);
    }
    FREE(data);
    printf("Successfully saved NetCDF inpcrd file \"%s\"\n", filename);
}
#else

void UnitIOSaveAmberCoordNetcdf(UNIT uUnit, char *filename)
{
    VPFATALEXIT("Built without NETCDF support. Rebuild with -DBINTRAJ\n");
}

int UnitIOSaveAmberParmNetcdf(const char *fname, UNITt *uUnit,
                              bool bPert, bool bPolar,
                              VARARRAY vaExcludedAtoms,
                              VARARRAY vaExcludedCount,
                              VARARRAY vaNBIndexMatrix,
                              VARARRAY vaNBParameters,
                              VARARRAY vaNBIndex)
{
    VPFATALEXIT("Built without NETCDF support. Rebuild with -DBINTRAJ\n");
    return 0;
}

#endif
