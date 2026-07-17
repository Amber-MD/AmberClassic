
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

// Data written to CMAP_INDEX, size=CMAP_KINDS (var 'mk' in old code)
typedef struct {
    int iAtom[5];
    int iParmIndex;
} SAVECMAPt;

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
 */
static int dest_cmap_index[256];
void UnitIOBuildCMAPTables(UNIT uUnit, PARMLIB plLib)
{
    int i, j, k, iNumDIH;
    int nPhi = 0, nPsi = 0;
    PARMSET psTemp;
    for (i=0;i<256;i++) dest_cmap_index[i]=-1;

    /* ── destroy any previous CMAP table ── */
    if (uUnit->vaCMAPs) {
        VP0("Rebuilding CMAP parameters.\n");
        VarArrayDestroy(&uUnit->vaCMAPs);
    } else
        VP0("Building CMAP parameters.\n");
    uUnit->vaCMAPs = vaVarArrayCreate(sizeof(SAVECMAPt));

    if (!GDefaults.iCMAP) return;

    iNumDIH = iVarArrayElementCount(uUnit->vaTorsions);
    if (iNumDIH <= 0) return;

    /* base pointers — one PVAI each, direct C array indexing throughout */
    SAVEATOMt    *saAtoms    = PVAI(uUnit->vaAtoms,    SAVEATOMt,    0);
    SAVETORSIONt *saTorsions = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);

    SAVETORSIONt **stdptPhi = malloc(iNumDIH * sizeof(SAVETORSIONt*));
    SAVETORSIONt **stdptPsi = malloc(iNumDIH * sizeof(SAVETORSIONt*));

    /* ── pre-filter: build phi and psi prospect lists ──
       Skip multi-term torsions (same 4 atoms as previous entry).
       bParmSetCMAPHasTorsion checks atom names only — fast pre-filter.
       Only searches parmlib — no local cache yet at this stage.       */
    for (i = 0; i < iNumDIH; i++) {
        SAVETORSIONt *t = &saTorsions[i];

        /* skip multi-term: same 4 atoms as previous torsion */
        if (i > 0 && !memcmp(&saTorsions[i-1].iAtom1,
                              &t->iAtom1,
                              4 * sizeof(int)))
            continue;

        const char *an[4];
        an[0] = sAtomName(saAtoms[t->iAtom1-1].aAtom);
        an[1] = sAtomName(saAtoms[t->iAtom2-1].aAtom);
        an[2] = sAtomName(saAtoms[t->iAtom3-1].aAtom);
        an[3] = sAtomName(saAtoms[t->iAtom4-1].aAtom);

        bool bPhi = FALSE, bPsi = FALSE;
        ParmLibParmSetLoop(plLib);
        while (bParmLibNextParmSet(plLib, &psTemp)) {
            bool bP1 = FALSE, bP2 = FALSE;
            bParmSetCMAPHasTorsion(psTemp, an, &bP1, &bP2);
            bPhi |= bP1; bPsi |= bP2;
            if (bPhi && bPsi) break;
        }

        if (bPhi) stdptPhi[nPhi++] = t;
        if (bPsi) stdptPsi[nPsi++] = t;
    }

    /* ── pair-match: phi (outer) × psi (inner) ── */
    for (i = 0; i < nPhi; i++) {
        SAVETORSIONt *stPhi = stdptPhi[i];

        /* build phi 4-atom data — reused across all psi candidates */
        const char *rn5[5], *an5[5];
        int absres5[5], ri5[5], iAtom[5];

        iAtom[0] = stPhi->iAtom1;
        iAtom[1] = stPhi->iAtom2;
        iAtom[2] = stPhi->iAtom3;
        iAtom[3] = stPhi->iAtom4;
        for (k = 0; k < 4; k++) {
            SAVEATOMt *sa = &saAtoms[iAtom[k] - 1];
            an5[k]     = sAtomName(sa->aAtom);
            rn5[k]     = sContainerName(cContainerWithin(sa->aAtom));
            absres5[k] = iContainerTempInt(cContainerWithin(sa->aAtom));
        }

        for (j = 0; j < nPsi; j++) {
            SAVETORSIONt *stPsi = stdptPsi[j];

            /* ── central bond cutoff: psi iAtom1,2,3 == phi iAtom2,3,4 ── */
            if (stPsi->iAtom1 != iAtom[1] ||
                stPsi->iAtom2 != iAtom[2] ||
                stPsi->iAtom3 != iAtom[3])
                continue;

            /* ── fill atom5: psi's iAtom4 ── */
            iAtom[4] = stPsi->iAtom4;
            SAVEATOMt *sb = &saAtoms[iAtom[4] - 1];
            an5[4]     = sAtomName(sb->aAtom);
            rn5[4]     = sContainerName(cContainerWithin(sb->aAtom));
            absres5[4] = iContainerTempInt(cContainerWithin(sb->aAtom));

            for (k = 0; k < 5; k++)
                ri5[k] = absres5[k] - absres5[0];

            int iIndex = iParmSetFindCMAP(uUnit->psParameters,
                                          rn5, an5, ri5);
            if (iIndex == PARM_NOT_FOUND) {
                int iTemp = PARM_NOT_FOUND;
                PARMLIB_LOOP(plLib, psTemp,
                    (iTemp = iParmSetFindCMAP(psTemp, rn5, an5, ri5)));

                // CMAP is supplemental — not found anywhere is normal so no error.
                if (iTemp == PARM_NOT_FOUND) continue;

                /* fetch deep copy and add to unit parmset */
                CMAPt cmap;
                ParmSetCMAP(psTemp, iTemp, &cmap, TRUE);
                
                iIndex = iParmSetAddCMAP(uUnit->psParameters, &cmap);
                dest_cmap_index[iTemp] = iIndex;
            }

            /* ── store SAVECMAPt entry ── */
            SAVECMAPt entry;
            for (k = 0; k < 5; k++) entry.iAtom[k] = iAtom[k];
            entry.iParmIndex = iIndex + 1;   /* 1-based, same as bonds */
            VarArrayAdd(uUnit->vaCMAPs, (GENP)&entry);
        }
    }

    free(stdptPhi);
    free(stdptPsi);
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
        ParmSetCMAP(uUnit->psParameters, i, &cmaps[i], FALSE);
        FortranWriteInt(cmaps[i].resolution);
    }
    FortranEndLine();

    int cmap_renum[256]={0};
    if (GDefaults.orig_cmap_order) { // FIXME: remove this compatibility feature hack
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
      if (GDefaults.orig_cmap_order) // FIXME: remove this compatibility feature hack
        FortranWriteInt(cmap_renum[saveCMAP->iParmIndex-1]);
      else
        FortranWriteInt(saveCMAP->iParmIndex);
    }
    FortranEndLine();
}

#ifdef BINTRAJ
/* ================================================================
   6. SaveAmberParmCMAPNetcdf()
   NetCDF equivalent of SaveAmberParmCMAP().
   Differences from Fortran path:
   - CMAP_COUNT: 2-element NC_INT array, not formatted string
   - CMAP_PARAMETER_nn: 2D NC_DOUBLE grid [res][res]
   - CMAP_INDEX: raw atom indices, (originally direct index with CMAP_INDEX)
   ================================================================ */
/* ================================================================
   SaveAmberParmCMAPNetcdf()

   Dimensions:
     CMAP_KINDS  — number of unique phi/psi interaction pairs (mk)
     CMAP_TYPES  — number of unique map parameter sets (maptypes)
     CMAP_RESOLUTION_nn  — scalar int, grid size (res x res) for maptype nn

   Index arrays:
     CMAP_INDEX_ATOMS(CMAP_KINDS, cmap_atom_quint)
       5 raw 1-based atom indices defining two adjacent torsions.
       Atoms 1-4 = phi torsion, atoms 2-5 = psi torsion.
       Shared central bond atoms 2-3 defines the phi/psi junction.
       No /3+1 AMBERINDEX transform — reader applies at runtime.
     CMAP_INDEX_MAP(CMAP_KINDS)
       1-based index into /cmap/map_nn/ subgroup.
   Parameter arrays:
     CMAP_PARAMETER_nn(CMAP_RESOLUTION_nn,CMAP_RESOLUTION_nn) — energy surface in kcal/mol
   ================================================================ */
int SaveAmberParmCMAPNetcdf(UNIT uUnit, int ncid)
{
    int CMAP_TYPES = iParmSetTotalCMAPParms(uUnit->psParameters);
    int CMAP_KINDS = iVarArrayElementCount(uUnit->vaCMAPs); // = NumPhiPsi
    if (!(CMAP_TYPES && CMAP_KINDS)) return NC_NOERR;
    if (!GDefaults.iCMAP) return NC_NOERR;

    /* ── write CMAP_KINDS and CMAP_TYPES scalars ── */
    {

    /* ── CMAP_INDEX_ATOMS and CMAP_INDEX_MAP ── */

        int dimid_types, dimid_kinds, dimid_quint;
        int vid_atoms, vid_map;
        int dims2[2];
        int *atoms_buf = malloc(CMAP_KINDS * 5 * sizeof(int));
        int *map_buf   = malloc(CMAP_KINDS * sizeof(int));

        NC_CHECK(nc_redef(ncid));

        NC_CHECK(nc_def_dim(ncid, "CMAP_TYPES", CMAP_TYPES, &dimid_types));
        NC_CHECK(nc_def_dim(ncid, "CMAP_KINDS", CMAP_KINDS, &dimid_kinds));
        NC_CHECK(nc_def_dim(ncid, "atom_quint", 5, &dimid_quint));

        dims2[0] = dimid_kinds; dims2[1] = dimid_quint;
        NC_CHECK(nc_def_var(ncid, "CMAP_INDEX_ATOMS", NC_INT, 2, dims2, &vid_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_atoms, "long_name", 50,
                                 "5 raw 1-based atom indices: atoms 1-4=phi, 2-5=psi"));
        NC_CHECK(nc_put_att_text(ncid, vid_atoms, "reference_dimension", 6, "NTOTAT"));

        NC_CHECK(nc_def_var(ncid, "CMAP_INDEX_MAP", NC_INT, 1, &dimid_kinds, &vid_map));
        NC_CHECK(nc_put_att_text(ncid, vid_map, "long_name", 79,
                "1-based index into CMAP_RESOLUTION_<nn> dimension and CMAP_PARAMETER_<nn> data"));
        NC_CHECK(nc_put_att_text(ncid, vid_atoms, "reference_dimenson", 10, "CMAP_TYPES"));
        NC_CHECK(nc_enddef(ncid));

        for (int i=0, idx=0; i < CMAP_KINDS; i++, idx+=5) {
            SAVECMAPt *saveCMAP = PVAI(uUnit->vaCMAPs, SAVECMAPt, i);
            for (int j=0;j<5;j++) atoms_buf[idx+j]=saveCMAP->iAtom[j];
            map_buf[i]=saveCMAP->iParmIndex;
        }
        NC_CHECK(nc_put_var_int(ncid, vid_atoms, atoms_buf));
        NC_CHECK(nc_put_var_int(ncid, vid_map,   map_buf));
        free(atoms_buf); free(map_buf);
    }

    /* ── subgroups — one per active map type ──
       Each contains coordinate variables phi and psi (angle values
       in degrees) plus the GRID energy surface.                     */
    for (int i = 0; i < CMAP_TYPES; i++) {
        int dimid_res, dims2[2];
        int vid_grid;
        CMAPt cmap;
        ParmSetCMAP(uUnit->psParameters, i, &cmap, FALSE);

        NC_CHECK(nc_redef(ncid));

        /* RESOLUTION_nn dimension shared by both phi and psi axes. */
        STRING sResolution;
        sprintf(sResolution,"CMAP_RESOLUTION_%02d",i+1);
        NC_CHECK(nc_def_dim(ncid, sResolution, cmap.resolution, &dimid_res));

        /* GRID(RESOLUTION, RESOLUTION) — same dim for both axes */
        dims2[0] = dims2[1] = dimid_res;
        STRING sParam;
        sprintf(sParam,"CMAP_PARAMETER_%02d",i+1);
        NC_CHECK(nc_def_var(ncid, sParam, NC_DOUBLE, 2, dims2, &vid_grid));
        NC_CHECK(nc_put_att_text(ncid, vid_grid, "units",      8, "kcal/mol"));
        NC_CHECK(nc_put_att_text(ncid, vid_grid, "long_name",  36,
                                 "CMAP energy grid: rows=phi, cols=psi"));
        NC_CHECK(nc_put_att_text(ncid, vid_grid, "coordinates", 18, sResolution));
        NC_CHECK(nc_put_att_text(ncid, vid_grid, "title", strlen(cmap.title), cmap.title));
        NC_CHECK(nc_enddef(ncid));

        NC_CHECK(nc_put_var_double(ncid, vid_grid, cmap.map));
    }
    return NC_NOERR;
}
#endif
