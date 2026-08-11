
#include        "basics.h"
#include        "classes.h"
#include        "avl.h"
#include        "defaults.h"
#include        "fortran.h"
#include        "mathop.h"
#include        "sort.h"
#include        "cmap.h"
#include        "unitio.h"

#ifdef BINTRAJ
#  include      "netcdf.h"

/* ----------------------------------------------------------------
   NC_CHECK: error handler for NetCDF API calls.
   ---------------------------------------------------------------- */
#define NC_CHECK(call) \
    do { int _e = (call); if (_e != NC_NOERR) { \
        fprintf(stderr, "NetCDF error %s: %s (%s:%d)\n", #call, nc_strerror(_e), __FILE__,__LINE__); \
        return _e; \
    }} while(0)

#endif

/*
 *  UnitIOBuildCMAPTables
 *
 *  Build uUnit->vaCMAPs (VarArray of SAVECMAPt) from the torsion list.
 *  Follows SAVEBONDt pattern: search local parmset first, then parmlib,
 *  copy params into uUnit->psParameters, store 1-based iParmIndex.
 *
 *  SAVECMAPt (replaces PHIPSI):
 *    typedef struct {
 *        int iAtom[5];    raw 1-based atom indices
 *        int iParmIndex;  1-based index into uUnit->psParameters CMAPs
 *    } SAVECMAPt;
 *
 *  Key design points:
 *  - Two prospect lists (stdptPhi/stdptPsi) replace single stdpt0[]
 *  - Multi-term torsion skip via memcmp on iAtom1..4
 *  - Inner loop cutoff: direct iAtom integer comparison
 *  - CMAP is supplemental — not found is normal, not an error
 *  - iParmSetFindCMAP on unit parmset handles deduplication
 *  - ParmSetGetCMAP with bCopy=true for transfer to local parmset
 *
 *  Note: CMAP = Correction map
 */
static int dest_cmap_index[256];
/*
 *  UnitIOBuildCMAPTables
 *
 *  Build uUnit->vaCMAP (VarArray of SAVECMAPt) from the torsion list.
 *
 *  Algorithm:
 *  1. Loop over all torsions, skip multi-term duplicates.
 *  2. Screen each torsion as a phi candidate (atoms 1-4) using
 *     bParmSetCMAPHasTorsion — atom names only, fast pre-filter.
 *  3. For each phi match, enumerate atom5 candidates directly from
 *     the bond graph: neighbors of atom4, excluding atom3.
 *  4. Do full 5-atom iParmSetFindCMAP for each (phi, atom5) pair.
 *  5. If found, add to unit parmset and store SAVECMAPt entry.
 *
 *  No prospect lists needed — bond graph replaces psi torsion search.
 *  CMAP is supplemental — not found is normal, not an error.
 */
void UnitIOBuildCMAPTables(UNIT uUnit, PARMLIB plLib)
{
    int i, k, iNumDIH;
    PARMSET psTemp;
    for (i=0;i<256;i++) dest_cmap_index[i]=-1;

    /* ── destroy any previous CMAP table ── */
    if (uUnit->vaCMAPs) {
        VP0("Rebuilding CMAP parameters.\n");
        VarArrayDestroy(&uUnit->vaCMAPs);
    } else
        VP0("Building CMAP parameters.\n");
    uUnit->vaCMAPs = vaVarArrayCreate(sizeof(SAVECMAPt));

    if (!GDefaults.bCMAP) return;

    iNumDIH = iVarArrayElementCount(uUnit->vaTorsions);
    if (iNumDIH <= 0) return;

    /* base pointers — one PVAI each, direct C array indexing throughout */
    SAVEATOMt    *saAtoms    = PVAI(uUnit->vaAtoms,    SAVEATOMt,    0);
    SAVETORSIONt *saTorsions = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);

    for (i = 0; i < iNumDIH; i++) {
        SAVETORSIONt *saTor = &saTorsions[i];

        /* skip multi-term: same 4 atoms as previous torsion */
        if (i > 0 && !memcmp(&saTorsions[i-1].iAtom1,
                              &saTor->iAtom1,
                              4 * sizeof(int)))
            continue;

        /* atom pointers for phi torsion (atoms 1-4) */
        SAVEATOMt *sa1 = &saAtoms[saTor->iAtom1 - 1];
        SAVEATOMt *sa2 = &saAtoms[saTor->iAtom2 - 1];
        SAVEATOMt *sa3 = &saAtoms[saTor->iAtom3 - 1];
        SAVEATOMt *sa4 = &saAtoms[saTor->iAtom4 - 1];

        /* an5[0..3] filled here, [4] filled per bond-graph neighbor below */
        const char *rn5[5], *an5[5];
        int absres5[5], ri5[5], iAtom5[5];

        iAtom5[0] = saTor->iAtom1;
        iAtom5[1] = saTor->iAtom2;
        iAtom5[2] = saTor->iAtom3;
        iAtom5[3] = saTor->iAtom4;

        an5[0] = sAtomName(sa1->aAtom);
        an5[1] = sAtomName(sa2->aAtom);
        an5[2] = sAtomName(sa3->aAtom);
        an5[3] = sAtomName(sa4->aAtom);

        rn5[0] = sContainerName(cContainerWithin(sa1->aAtom));
        rn5[1] = sContainerName(cContainerWithin(sa2->aAtom));
        rn5[2] = sContainerName(cContainerWithin(sa3->aAtom));
        rn5[3] = sContainerName(cContainerWithin(sa4->aAtom));

        absres5[0] = iContainerTempInt(cContainerWithin(sa1->aAtom));
        absres5[1] = iContainerTempInt(cContainerWithin(sa2->aAtom));
        absres5[2] = iContainerTempInt(cContainerWithin(sa3->aAtom));
        absres5[3] = iContainerTempInt(cContainerWithin(sa4->aAtom));

        /* ── screen as phi: an5[0..3], fast ── */
        bool bPhi = false;
        ParmLibParmSetLoop(plLib);
        while (bParmLibNextParmSet(plLib, &psTemp)) {
            if (bParmSetCMAPHasPhi(psTemp, an5, rn5)) { bPhi = true; break; }
        }
        if (!bPhi) continue;

        /* ── enumerate atom5 candidates from bond graph ──
           neighbors of atom4, excluding atom3               */
        int nBonded = iAtomCoordination(sa4->aAtom);
        for (int nb = 0; nb < nBonded; nb++) {
            ATOMt *aAtom5 = aAtomBondedNeighbor(sa4->aAtom, nb);
            if (aAtom5 == sa3->aAtom) continue;  /* exclude atom3 */

            iAtom5[4] = iContainerTempInt(aAtom5);  /* 1-based atom index */
            SAVEATOMt *sa5 = &saAtoms[iAtom5[4] - 1];
            an5[4] = sAtomName(sa5->aAtom);
            rn5[4] = sContainerName(cContainerWithin(sa5->aAtom));
            absres5[4] = iContainerTempInt(cContainerWithin(sa5->aAtom));

            for (k = 0; k < 5; k++)
                ri5[k] = absres5[k] - absres5[0];

            /* ── lookup: unit parmset first, then parmlib ── */
            int iIndex = iParmSetFindCMAP(uUnit->psParameters, rn5, an5, ri5);
            if (iIndex == PARM_NOT_FOUND) {
                int iTemp = PARM_NOT_FOUND;
                PARMLIB_LOOP(plLib, psTemp,
                    (iTemp = iParmSetFindCMAP(psTemp, rn5, an5, ri5)));

                /* CMAP is supplemental — not found anywhere is normal so no error. */
                if (iTemp == PARM_NOT_FOUND) continue;

                /* fetch deep copy (bCopy=true) and add to unit parmset */
                CMAPt cmap;
                ParmSetCMAP(psTemp, iTemp, &cmap, true);
                iIndex = iParmSetAddCMAP(uUnit->psParameters, &cmap);
                dest_cmap_index[iTemp] = iIndex; // compaitibility hack
            }

            /* ── store SAVECMAPt entry ── */
            SAVECMAPt entry;
            for (k = 0; k < 5; k++) entry.iAtom[k] = iAtom5[k];
            entry.iParmIndex = iIndex + 1;   /* 1-based */
            VarArrayAdd(uUnit->vaCMAPs, (GENP)&entry);
        }
    }
}

void SaveAmberParmCMAP(UNIT uUnit, FILE * fOut)
{
    STRING sComment;
    int CMAP_TYPES = iParmSetTotalCMAPParms(uUnit->psParameters);
    int CMAP_KINDS = iVarArrayElementCount(uUnit->vaCMAPs);
    if (!(CMAP_TYPES && CMAP_KINDS)) return;


    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CMAP_COUNT");
    if (GDefaults.dPrmtopFormat>1.0) {
        FortranWriteString("%COMMENT SIZE:2; MEMBERS:CMAP_KINDS,CMAP_TYPES;");
        FortranWriteString("%COMMENT DESC:CMAP_KINDS=phi+psi set count, CMAP_TYPES=Energy map count;");
    }
    FortranWriteString("%FORMAT(2I8)");
    FortranFormat(10, "%8d");
    FortranWriteInt(CMAP_KINDS);
    FortranWriteInt(CMAP_TYPES);
    FortranEndLine();

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CMAP_RESOLUTION");
    if (GDefaults.dPrmtopFormat>1.0) {
        sprintf(sComment, "%%COMMENT DIMENSION:CMAP_TYPES; SIZE:%d;", CMAP_TYPES);
        FortranWriteString(sComment);
    }
    FortranWriteString("%FORMAT(20I4)");
    FortranFormat(20, "%4d");
    CMAP cmaps = MALLOC(sizeof(CMAPt)*CMAP_TYPES);
    for (int i = 0; i < CMAP_TYPES; i++) {
        ParmSetCMAP(uUnit->psParameters, i, &cmaps[i], false);
        FortranWriteInt(cmaps[i].resolution);
    }
    FortranEndLine();

    int cmap_renum[256]={0};
    if (GDefaults.bOrigCMAPOrder) { // FIXME: remove this compatibility feature hack
        int iIndex=0;
        for (int irenum=0;irenum<256;irenum++) if (dest_cmap_index[irenum]>=0) {
            int i = dest_cmap_index[irenum];
            cmap_renum[i]=++iIndex;
            FortranFormat(1, "%-80s");
            sprintf(sComment, "%%FLAG CMAP_PARAMETER_%02d",iIndex);
            FortranWriteString(sComment);
            int msize = cmaps[i].resolution * cmaps[i].resolution;
            if (GDefaults.dPrmtopFormat>1.0) {
                sprintf(sComment, "%%COMMENT DIMENSION:CMAP_RESOLUTION(%d),"
                       "CMAP_RESOLUTION(%d); SIZE:%d;",iIndex,iIndex,msize);
                FortranWriteString(sComment);
                FortranWriteString("%COMMENT UNITS:kcal/mol; DESC:CMAP energy grid: rows=phi, cols=psi;");
                sprintf(sComment, "%%COMMENT TITLE:%s;", cmaps[i].title);
            } else sprintf(sComment, "%%COMMENT  %s", cmaps[i].title);
            FortranWriteString(sComment);
            FortranWriteString("%FORMAT(8F9.5)");
            FortranFormat(8, "%9.5lf");
            for (int j = 0; j < msize; j ++) FortranWriteDouble(cmaps[i].map[j]);
            FortranEndLine();
        }
    } else {
    for (int i = 0; i < CMAP_TYPES; i++) {
        FortranFormat(1, "%-80s");
        sprintf(sComment, "%%FLAG CMAP_PARAMETER_%02d",i+1);
        FortranWriteString(sComment);
        int msize = cmaps[i].resolution * cmaps[i].resolution;
        if (GDefaults.dPrmtopFormat>1.0) {
            sprintf(sComment, "%%COMMENT DIMENSION:CMAP_RESOLUTION(%d),"
                   "CMAP_RESOLUTION(%d); SIZE:%d;",i+1,i+1,msize);
            FortranWriteString(sComment);
            FortranWriteString("%COMMENT UNITS:kcal/mol; DESC:CMAP energy grid: rows=phi, cols=psi;");
            sprintf(sComment, "%%COMMENT TITLE:%s;", cmaps[i].title);
        } else sprintf(sComment, "%%COMMENT %s", cmaps[i].title);
        FortranWriteString(sComment);
        FortranWriteString("%FORMAT(8F9.5)");
        FortranFormat(8, "%9.5lf");
        for (int j = 0; j < msize; j ++) FortranWriteDouble(cmaps[i].map[j]);
        FortranEndLine();
    }
    }

    FREE(cmaps);

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CMAP_INDEX");
    if (GDefaults.dPrmtopFormat>1.0) {
        sprintf(sComment, "%%COMMENT DIMENSION:CMAP_KINDS; STRIDE:6; SIZE:%d;", CMAP_KINDS*6);
        FortranWriteString(sComment);
        FortranWriteString("%COMMENT MEMBERS:ATOM_I,ATOM_J,ATOM_K,ATOM_L,ATOM_M,PARM_INDEX;");
        FortranWriteString("%COMMENT DESC:INDEX=1-based index into CMAP_RESOLUTION_<nn> CMAP_PARAMETER_<nn>");
    }
    FortranWriteString("%FORMAT(6I8)");
    FortranFormat(6, "%8d");
    for (int i = 0; i < CMAP_KINDS; i++) {
        SAVECMAPt *saveCMAP = PVAI(uUnit->vaCMAPs, SAVECMAPt, i);
        for (int j=0;j<5;j++) FortranWriteInt(saveCMAP->iAtom[j]);
      if (GDefaults.bOrigCMAPOrder) // FIXME: remove this compatibility feature hack
        FortranWriteInt(cmap_renum[saveCMAP->iParmIndex-1]);
      else
        FortranWriteInt(saveCMAP->iParmIndex);
    }
    FortranEndLine();
}


