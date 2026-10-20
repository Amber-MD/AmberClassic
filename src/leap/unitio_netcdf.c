#define AMBER_VERSION ""
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
#include        <assert.h>
#include        <netcdf.h>
#include        <stdarg.h>

#define LABLEN 4
#define MAXDESCLEN MAXSTRINGLENGTH

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

/* ----------------------------------------------------------------
   NC_CHECK: error handler for NetCDF API calls.
   ---------------------------------------------------------------- */
#define NC_CHECK(call) \
    do { int _e = (call); if (_e != NC_NOERR) { \
        VPFATAL("NetCDF error %s: %s (%s:%d)\n", #call, nc_strerror(_e), __FILE__,__LINE__); \
        return; \
    }} while(0)

/* convenience: put a text attribute using strlen() */
#define NC_ATT(ncid, vid, name, val) \
    do { const char *_s = (val); \
         NC_CHECK(nc_put_att_text(ncid, vid, name, strlen(_s), _s)); \
    } while(0)

/* ----------------------------------------------------------------
   nc_defvar: defines a variable and attaches its attributes.
   amber_name/units/long_name are given guaranteed positional slots
   since almost every variable sets some subset of these three; pass
   NULL for any that don't apply. Any remaining, less-common
   attributes (note, reference, index_dimension, pdb_field,
   mmcif_field, flag_values, flag_meanings, encoding, layout, ...)
   follow as NULL-terminated (name, value) pairs. A NULL value in the
   varargs list (not a NULL name) silently skips that attribute, so
   callers can pass optional fields without an `if` guard.

   Internal netCDF calls are checked manually (rc returned, not
   NC_CHECK'd) because NC_CHECK's error message stringifies the
   *literal* call expression -- inside a shared helper that would
   always print the same generic "nc_def_var(ncid, name, type, ...)"
   text rather than which actual variable failed. Callers wrap the
   whole nc_defvar() call in NC_CHECK instead, so the reported
   expression is the meaningful one (e.g. the real variable name).   */
static int
nc_defvar(int ncid, int *pvid, const char *name, nc_type type,
          int ndims, const int *dimids,
          const char *amber_name, const char *units, const char *long_name,
          ...)
{
    int rc = nc_def_var(ncid, name, type, ndims, dimids, pvid);
    if (rc != NC_NOERR) return rc;

    if (amber_name) {
        rc = nc_put_att_text(ncid, *pvid, "amber_name", strlen(amber_name), amber_name);
        if (rc != NC_NOERR) return rc;
    }
    if (units) {
        rc = nc_put_att_text(ncid, *pvid, "units", strlen(units), units);
        if (rc != NC_NOERR) return rc;
    }
    if (long_name) {
        rc = nc_put_att_text(ncid, *pvid, "long_name", strlen(long_name), long_name);
        if (rc != NC_NOERR) return rc;
    }

    va_list ap;
    va_start(ap, long_name);
    const char *attname;
    while ((attname = va_arg(ap, const char *)) != NULL) {
        const char *attval = va_arg(ap, const char *);
        if (attval) {
            rc = nc_put_att_text(ncid, *pvid, attname, strlen(attval), attval);
            if (rc != NC_NOERR) { va_end(ap); return rc; }
        }
    }
    va_end(ap);
    return NC_NOERR;
}

/* NC_PARAM: for a scalar NC_INT variable whose value is already known
   at definition time (e.g. a dimension length, or a fixed alias flag)
   but can't be written until after nc_enddef(). Defines the var, tags
   it "amber_name", and stashes (vid, val) in the param_vid/param_val
   arrays for the deferred write loop right after nc_enddef().         */
#define NC_PARAM(ncid, name, amber_name, val)                           \
    do {                                                                \
        int _vid;                                                       \
        NC_CHECK(nc_def_var(ncid, name, NC_INT, 0, NULL, &_vid));       \
        if (amber_name) NC_ATT(ncid, _vid, "amber_name", (amber_name)); \
        param_vid[nparam] = _vid; param_val[nparam++] = (val);          \
    } while(0)

/* NC_DIM: defines a dimension of the given length, plus a same-named
   scalar variable carrying "amber_name" (dimensions themselves carry
   no attributes in netCDF, so the amber_name tag is carried by this
   companion var, written via the same deferred param_vid/param_val
   mechanism as NC_PARAM -- classic/64-bit-offset netCDF is strictly
   define-then-data, so the companion var's value can't be set inline
   at nc_def_var time either).                                          */
#define NC_DIM(ncid, dimname, amber_name, len, pdimid)          \
    do {                                                        \
        NC_CHECK(nc_def_dim(ncid, (dimname), (len), (pdimid))); \
        NC_PARAM(ncid, (dimname), (amber_name), (len));         \
    } while(0)

void
UnitIOSaveAmberParmNetcdf(const char *fname, UNITt *uUnit,
                              bool bPert, bool bPolar, int iMaxAtoms,
                              VARARRAY vaExcludedAtoms,
                              VARARRAY vaExcludedCount,
                              VARARRAY vaNBIndexMatrix,
                              VARARRAY vaNBParameters,
                              VARARRAY vaNBIndex)
{
    int ncid;
    int dimid_atoms, vid_extrapoints, dimid_res;
    int dimid_bp, dimid_ap, dimid_atyp;
    int dimid_name4, dimid_chainid_len, dimid_icode_len, dimid_restype_len;
    int dimid_atom2, dimid_atom3, dimid_atom4, dimid_atom5;
    int n_atoms    = iVarArrayElementCount(uUnit->vaAtoms);
    int n_residues = iVarArrayElementCount(uUnit->vaResidues);
    char s1[8], s2[8];
    STRING sDesc, sFFType;
    double dKb, dR0, dKpull, dRpull0, dKpress, dRpress0;
    double dKt, dT0;
    double dKp, dP0, dScEE, dScNB;
    double dC, dD;
    int vid_name, vid_charge, vid_atomic_num, vid_itype;
    int param_vid[50], param_val[50], nparam=0;

    if (bPert) VPFATALEXIT("Pert atoms not supported\n");

    int iParmSets = iVarArrayElementCount(uUnit->vaParmSets);
    int title_length = 0, fname_length = 0;
    {
        int iPARM = 0;
        PARMSET psParm;
        for (int i = 0; i < iParmSets; i++) {
            psParm = *PVAI(uUnit->vaParmSets, PARMSET, i);
            int len = (int)strlen(psParm->sFname);
            if (len > fname_length) fname_length = len;
            len = (int)strlen(psParm->sTitle);
            if (len > title_length) title_length = len;
            if (!strncmp(psParm->sTitle, "PARM", 4)) iPARM = i;
        }
        psParm = *PVAI(uUnit->vaParmSets, PARMSET, iPARM);
        sprintf(sFFType, "AMBER " AMBER_VERSION " %.1024s", psParm->sTitle);
    }
    int chainid_len = 1;
    int icode_len = 1;
    int restype_len = 1;
    for (int i = 0; i < n_residues; i++) {
        const char *s = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sChainId;
        size_t l = strlen(s);
        if ((int)l > chainid_len) chainid_len = (int)l;
    }

    int CMAP_TYPES = iParmSetTotalCMAPParms(uUnit->psParameters);
    int CMAP_KINDS = iVarArrayElementCount(uUnit->vaCMAPs);

    /* ── counts for dimension definitions (no allocation) ── */
    int n_bonds  = iVarArrayElementCount(uUnit->vaBonds);
    int n_angles = iVarArrayElementCount(uUnit->vaAngles);
    int n_tors   = iVarArrayElementCount(uUnit->vaTorsions);
    int nr_bond  = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
    int nr_angle = iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
    int nr_tors  = iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION);
    int n_bond_total  = n_bonds  + nr_bond;
    int n_angle_total = n_angles + nr_angle;

    /* ── pre-build ALL torsion arrays in one pass ──
       nd = periodic interaction count, nh = harmonic interaction count.
       harm_map/dihe_map: stacked iParmIndex → split table index.

       dihedral_term numbering (replaces AMBER's sign-flip-on-atom3
       encoding of "calc14"): 0 = solo term (only Fourier term for this
       atom quad); for a multi-term set, ALL terms are numbered from 1
       (1, 2, 3, ...) -- including what would have been "term 0" -- so
       a reader can tell "this record belongs to a multi-term set" from
       a single row, without scanning neighbors. Harmonic terms are
       always solo and never take part in this numbering.

       Implementation: no lookahead/cache needed. Since harmonic terms
       are siphoned into a separate array, the previous PERIODIC record
       is always exactly dihe_term[nd-1] at the moment a repeat of the
       same 4 atoms is detected, so we can backpatch it in place.        */
    int n_dihe_total = 0, n_harm_total = 0;
    int nd_parm = 0, nh_parm = 0;
    bool bBondAug = bBondAugmentationFound(uUnit);

    int nt = iParmSetTotalTorsionParms(uUnit->psParameters);
    int ni = iParmSetTotalImproperParms(uUnit->psParameters);

    /* ── pass 1: classify parameters, build index maps ──
       nd_parm/nh_parm count only regular torsion/improper parm slots
       here; restraint-sourced parm slots (periodicity read directly off
       each SAVERESTRAINTt, not from the shared parm table) are appended
       past these counts once nr_tors_dihe/nr_tors_harm are known below. */
    int total_parms = nt + ni;  /* restraints share parm table */
    int *harm_map = CALLOC(sizeof(int) * (total_parms + 1));
    int *dihe_map = CALLOC(sizeof(int) * (total_parms + 1));
    {
        int iN;
        for (int i = 0; i < nt; i++) {
            ParmSetTorsion(uUnit->psParameters, i, NULL,NULL,NULL,NULL,
                           &iN, NULL, NULL, NULL, NULL, NULL);
            if (!iN) harm_map[i] = ++nh_parm;
            else     dihe_map[i] = ++nd_parm;
        }
        for (int i = 0; i < ni; i++) {
            ParmSetImproper(uUnit->psParameters, i, NULL,NULL,NULL,NULL,
                            &iN, NULL, NULL, NULL);
            if (!iN) harm_map[nt + i] = ++nh_parm;
            else     dihe_map[nt + i] = ++nd_parm;
        }
    }

    /* count how many torsion restraints are periodic vs harmonic, so the
       parm-value buffers below can be sized to hold both regular terms
       and restraint-sourced terms in one contiguous range each.         */
    int nr_tors_dihe = 0, nr_tors_harm = 0;
    if (nr_tors) {
        SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++) {
            if (sr->iType != RESTRAINTTORSION) continue;
            if (sr->dN == 0) nr_tors_harm++;
            else             nr_tors_dihe++;
        }
    }
    int nd_parm_total = nd_parm + nr_tors_dihe;
    int nh_parm_total = nh_parm + nr_tors_harm;

    /* ── pass 2: fill interaction and parameter arrays ── */
    int nmax = n_tors + nr_tors;
    int (*dihe_atoms)[4]  = MALLOC(nmax * sizeof(*dihe_atoms));
    int  *dihe_parm       = MALLOC(nmax * sizeof(int));
    char *dihe_proper     = MALLOC(nmax * sizeof(char));
    int  *dihe_term       = MALLOC(nmax * sizeof(int));
    double *dihe_force    = MALLOC(nd_parm_total * sizeof(double));
    int    *dihe_period   = MALLOC(nd_parm_total * sizeof(int));
    double *dihe_phase    = MALLOC(nd_parm_total * sizeof(double));
    double *dihe_scee     = MALLOC(nd_parm_total * sizeof(double));
    double *dihe_scnb     = MALLOC(nd_parm_total * sizeof(double));

    int (*harm_atoms)[4]  = MALLOC(nmax * sizeof(*harm_atoms));
    int  *harm_parm       = MALLOC(nmax * sizeof(int));
    char *harm_restraint  = MALLOC(nmax * sizeof(char));
    double *harm_force    = MALLOC(nh_parm_total * sizeof(double));
    double *harm_equil    = MALLOC(nh_parm_total * sizeof(double));
    int n_dihe_restraints = 0, n_harm_restraints = 0;

    /* fill parameter value arrays from parmset */
    {
        int iN;
        for (int i = 0; i < nt; i++) {
            ParmSetTorsion(uUnit->psParameters, i, NULL,NULL,NULL,NULL,
                           &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc);
            if (harm_map[i]) {
                int hi = harm_map[i] - 1;
                harm_force[hi] = dKp;
                harm_equil[hi] = dP0 * RADTODEG;
            } else {
                int di = dihe_map[i] - 1;
                dihe_force[di]  = dKp;
                dihe_period[di] = iN;
                dihe_phase[di]  = dP0 * RADTODEG;
                dihe_scee[di]   = (dScEE < 0.0) ? GDefaults.dSceeScaleFactor : dScEE;
                dihe_scnb[di]   = (dScNB < 0.0) ? GDefaults.dScnbScaleFactor : dScNB;
            }
        }
        for (int i = 0; i < ni; i++) {
            ParmSetImproper(uUnit->psParameters, i, NULL,NULL,NULL,NULL,
                            &iN, &dKp, &dP0, sDesc);
            if (harm_map[nt + i]) {
                int hi = harm_map[nt + i] - 1;
                harm_force[hi] = dKp;
                harm_equil[hi] = dP0 * RADTODEG;
            } else {
                int di = dihe_map[nt + i] - 1;
                dihe_force[di]  = dKp;
                dihe_period[di] = iN;
                dihe_phase[di]  = dP0 * RADTODEG;
                dihe_scee[di]   = 0.0;
                dihe_scnb[di]   = 0.0;
            }
        }
    }

    /* torsion restraint parm values: read directly off each
       SAVERESTRAINTt (dKx=force, dX0=phase, dN=periodicity -- restraints
       carry no SCEE/SCNB, written 0.0 and identified via *_restraint==1),
       appended past the regular nd_parm/nh_parm ranges in the same
       encounter order used by the interaction-fill loop above, so
       dihe_parm[nd]/harm_parm[nh] there line up with these slots here.  */
    if (nr_tors) {
        int iRestDihe = nd_parm, iRestHarm = nh_parm;
        SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++) {
            if (sr->iType != RESTRAINTTORSION) continue;
            if (sr->dN == 0) {
                harm_force[iRestHarm] = sr->dKx;
                harm_equil[iRestHarm] = sr->dX0 * RADTODEG;
                iRestHarm++;
            } else {
                dihe_force[iRestDihe]  = sr->dKx;
                dihe_period[iRestDihe] = sr->dN;
                dihe_phase[iRestDihe]  = sr->dX0 * RADTODEG;
                dihe_scee[iRestDihe]   = 0.0;
                dihe_scnb[iRestDihe]   = 0.0;
                iRestDihe++;
            }
        }
    }

    /* fill interaction arrays.
       Single pass per has-H/without-H bucket; dihedral_term is assigned
       and backpatched in place as each bucket is filled (see note above
       pass 2). Since same-atom-quad => same H-status by construction,
       every multi-term group stays entirely within one bucket, so the
       "previous periodic record is dihe_term[nd-1]" invariant holds
       across the has-H/without-H split.                                 */
    int dihe_nh = 0, dihe_nnh = 0;
    {
        int nd = 0, nh = 0;
        int iTerm;
        int iPrevAtom1, iPrevAtom2, iPrevAtom3, iPrevAtom4;
        bool has_h[n_tors];

        /* has-H bucket */
        iTerm = 0;
        iPrevAtom1 = iPrevAtom2 = iPrevAtom3 = iPrevAtom4 = -1;
        for (int i = 0; i < n_tors; i++) {
            SAVETORSIONt *t = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
            int pidx = t->iParmIndex - 1;  /* 0-based parm index */

            if (harm_map[pidx]) {
                harm_atoms[nh][0] = t->iAtom1;
                harm_atoms[nh][1] = t->iAtom2;
                harm_atoms[nh][2] = t->iAtom3;
                harm_atoms[nh][3] = t->iAtom4;
                harm_parm[nh]      = harm_map[pidx];
                harm_restraint[nh] = 0;
                nh++;
                continue;
            }

            ATOM aA = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom1 - 1)->aAtom;
            ATOM aB = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom2 - 1)->aAtom;
            ATOM aC = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom3 - 1)->aAtom;
            ATOM aD = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom4 - 1)->aAtom;
            has_h[i] = iAtomElement(aA) == HYDROGEN ||
                       iAtomElement(aB) == HYDROGEN ||
                       iAtomElement(aC) == HYDROGEN ||
                       iAtomElement(aD) == HYDROGEN;
            if (!has_h[i]) continue;

            bool bSameQuad = (t->iAtom1 == iPrevAtom1 && t->iAtom2 == iPrevAtom2 &&
                              t->iAtom3 == iPrevAtom3 && t->iAtom4 == iPrevAtom4);
            if (bSameQuad) {
                if (iTerm == 0) dihe_term[nd - 1] = 1;  /* backpatch solo -> first-of-multi */
                iTerm++;
            } else {
                iTerm = 0;
            }
            iPrevAtom1 = t->iAtom1; iPrevAtom2 = t->iAtom2;
            iPrevAtom3 = t->iAtom3; iPrevAtom4 = t->iAtom4;

            dihe_atoms[nd][0]  = t->iAtom1;
            dihe_atoms[nd][1]  = t->iAtom2;
            dihe_atoms[nd][2]  = t->iAtom3;
            dihe_atoms[nd][3]  = t->iAtom4;
            dihe_parm[nd]      = dihe_map[pidx];
            dihe_proper[nd]    = t->bProper ? 1 : 0;
            dihe_term[nd]      = iTerm;
            nd++;
            dihe_nh++;
        }

        /* without-H bucket */
        iTerm = 0;
        iPrevAtom1 = iPrevAtom2 = iPrevAtom3 = iPrevAtom4 = -1;
        for (int i = 0; i < n_tors; i++) {
            SAVETORSIONt *t = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
            int pidx = t->iParmIndex - 1;  /* 0-based parm index */
            if (harm_map[pidx] || has_h[i]) continue;

            bool bSameQuad = (t->iAtom1 == iPrevAtom1 && t->iAtom2 == iPrevAtom2 &&
                              t->iAtom3 == iPrevAtom3 && t->iAtom4 == iPrevAtom4);
            if (bSameQuad) {
                if (iTerm == 0) dihe_term[nd - 1] = 1;  /* backpatch solo -> first-of-multi */
                iTerm++;
            } else {
                iTerm = 0;
            }
            iPrevAtom1 = t->iAtom1; iPrevAtom2 = t->iAtom2;
            iPrevAtom3 = t->iAtom3; iPrevAtom4 = t->iAtom4;

            dihe_atoms[nd][0]  = t->iAtom1;
            dihe_atoms[nd][1]  = t->iAtom2;
            dihe_atoms[nd][2]  = t->iAtom3;
            dihe_atoms[nd][3]  = t->iAtom4;
            dihe_parm[nd]      = dihe_map[pidx];
            dihe_proper[nd]    = t->bProper ? 1 : 0;
            dihe_term[nd]      = iTerm;
            nd++;
            dihe_nnh++;
        }

        /* torsion restraints: always solo (restraints don't multi-term) */
        if (nr_tors) {
            SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
            for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++) {
                if (sr->iType != RESTRAINTTORSION) continue;
                int pidx = sr->iParmIndex - 1;
                if (harm_map[pidx]) {
                    harm_atoms[nh][0] = sr->iAtom1;
                    harm_atoms[nh][1] = sr->iAtom2;
                    harm_atoms[nh][2] = sr->iAtom3;
                    harm_atoms[nh][3] = sr->iAtom4;
                    harm_parm[nh]      = harm_map[pidx];
                    n_harm_restraints++;
                    nh++;
                } else {
                    dihe_atoms[nd][0]  = sr->iAtom1;
                    dihe_atoms[nd][1]  = sr->iAtom2;
                    dihe_atoms[nd][2]  = sr->iAtom3;
                    dihe_atoms[nd][3]  = sr->iAtom4;
                    dihe_parm[nd]      = dihe_map[pidx];
                    dihe_proper[nd]    = 1;
                    dihe_term[nd]      = 0;
                    n_dihe_restraints++;
                    nd++;
                }
            }
        }
        n_dihe_total = nd;
        n_harm_total = nh;
    }
    FREE(harm_map); FREE(dihe_map);

    START(tNetCDF);
    NC_CHECK(nc_create(fname, NC_CLOBBER|NC_64BIT_OFFSET|NC_NOFILL, &ncid));

    /* Global attributes */
    {
        time_t tp; char sDate[64];
        time(&tp);
        strftime(sDate, sizeof(sDate), "%Y-%m-%dT%H:%M:%S", localtime(&tp));
        NC_ATT(ncid, NC_GLOBAL, "date", sDate);
        NC_ATT(ncid, NC_GLOBAL, "date_conventions", "ISO 8601");
    }
    NC_ATT(ncid, NC_GLOBAL, "title", sContainerName(uUnit));
    NC_ATT(ncid, NC_GLOBAL, "Conventions", "MolecularDynamicsTopology");
    NC_ATT(ncid, NC_GLOBAL, "ConventionVersion", "1.0");
    NC_ATT(ncid, NC_GLOBAL, "application", "AMBER");
    NC_ATT(ncid, NC_GLOBAL, "applicationVersion", AMBER_VERSION);
    NC_ATT(ncid, NC_GLOBAL, "program", "leap");
    NC_ATT(ncid, NC_GLOBAL, "programVersion", LEAP_VERSION);
    NC_ATT(ncid, NC_GLOBAL, "force_field_type", sFFType);
    NC_ATT(ncid, NC_GLOBAL, "index_base", "1");
    NC_ATT(ncid, NC_GLOBAL, "index_base_note",
           "All atom, residue, molecule, and parameter-type indices in this "
           "file are 1-based unless otherwise noted. This does not apply to "
           "boolean flags (stored as byte 0/1) or ordinal counters such as "
           "dihedral_term, which are documented individually.");
    NC_ATT(ncid, NC_GLOBAL, "nonbond_lower_triangle_echelon_indexing","max(I,J)*(max(I,J)+1)/2 + min(I,J) + 1");
    NC_ATT(ncid, NC_GLOBAL, "nonbond_lower_triangle_echelon_variables", "I = atom_nonbond_type_index[atom_A], J = atom_nonbond_type_index[atom_B]");

    /* Shared dimensions */
    NC_DIM(ncid, "atom", "NATOMS", n_atoms, &dimid_atoms);
    NC_DIM(ncid, "residue", "NRES", n_residues, &dimid_res);
    NC_CHECK(nc_def_dim(ncid, "label_length", LABLEN, &dimid_name4));
    NC_CHECK(nc_def_dim(ncid, "chainid_length", chainid_len, &dimid_chainid_len));
    NC_CHECK(nc_def_dim(ncid, "icode_length", icode_len, &dimid_icode_len));
    NC_CHECK(nc_def_dim(ncid, "restype_length", restype_len, &dimid_restype_len));
    NC_CHECK(nc_def_dim(ncid, "atom_pair", 2, &dimid_atom2));
    NC_CHECK(nc_def_dim(ncid, "atom_triple", 3, &dimid_atom3));
    NC_CHECK(nc_def_dim(ncid, "atom_quad", 4, &dimid_atom4));
    if (CMAP_TYPES && CMAP_KINDS && GDefaults.bCMAP)
        NC_CHECK(nc_def_dim(ncid, "atom_quint", 5, &dimid_atom5));

    /* parameter set metadata */
    int dimid_fname_len, dimid_title_len, dimid_nparm;
    int vid_parm_file, vid_parm_title;
    int dimid_hash_length, vid_parm_hash, hash_length = 0; 
    for (int i = 0; i < iParmSets; i++) {
      PARMSET psParm = *PVAI(uUnit->vaParmSets, PARMSET, i);
      int len = strlen(psParm->cksum);
      if (len>hash_length) hash_length = len;
    }
    NC_CHECK(nc_def_dim(ncid, "parmset_title_length", title_length, &dimid_title_len));
    NC_CHECK(nc_def_dim(ncid, "parmset_fname_length", fname_length, &dimid_fname_len));
    NC_CHECK(nc_def_dim(ncid, "parmset_hash_length", hash_length, &dimid_hash_length));
    NC_CHECK(nc_def_dim(ncid, "parameter_set", iParmSets, &dimid_nparm));
    {
        int dims2[2] = {dimid_nparm, dimid_fname_len};
        NC_CHECK(nc_defvar(ncid, &vid_parm_file, "parm_file", NC_CHAR, 2, dims2,
                            NULL, NULL, NULL, NULL));
        dims2[1] = dimid_title_len;
        NC_CHECK(nc_defvar(ncid, &vid_parm_title, "parm_title", NC_CHAR, 2, dims2,
                            NULL, NULL, NULL, NULL));
        dims2[1] = dimid_hash_length;
        NC_CHECK(nc_defvar(ncid, &vid_parm_hash, "parm_hash", NC_CHAR, 2, dims2,
                            NULL, NULL,
                            "POSIX cksum(1) result for parm_file, formatted as \"<checksum>-<size>\"", NULL));
    }

    /* cell dimensions */
    int vid_cell_length, vid_cell_angle;
    if (bUnitUseBox(uUnit)) {
        int dim;
        NC_CHECK(nc_def_dim(ncid, "cell_spatial", 3, &dim));
        NC_CHECK(nc_defvar(ncid, &vid_cell_length, "cell_lengths", NC_DOUBLE, 1, &dim,
                            NULL, "angstrom", NULL,
                            "comment", "initial cell dimension is primarily for fixed dimension systems",
                            NULL));
        NC_CHECK(nc_def_dim(ncid, "cell_angular", 3, &dim));
        NC_CHECK(nc_defvar(ncid, &vid_cell_angle, "cell_angle", NC_DOUBLE, 1, &dim,
                            NULL, "degree", NULL, NULL));
    }

    /* atom properties */
    {
        int dims2[2] = {dimid_atoms, dimid_name4};
        NC_CHECK(nc_defvar(ncid, &vid_name, "atom_name", NC_CHAR, 2, dims2,
                            "ATOM_NAME", NULL, "atom name",
                            "pdb_name", "name",
                            "mmcif_name", GDefaults.bCIFReadAuth ? "atom_site.auth_atom_id"
                                                                   : "atom_site.label_atom_id",
                            NULL));
    }
    NC_CHECK(nc_defvar(ncid, &vid_charge, "atom_charge", NC_DOUBLE, 1, &dimid_atoms,
                        "CHARGE", "elementary charge", "partial atomic charge",
                        "comment", "multiply by 18.2223 for AMBER internal charge units",
                        NULL));
    NC_CHECK(nc_defvar(ncid, &vid_atomic_num, "atom_atomic_number", NC_INT, 1, &dimid_atoms,
                        "ATOMIC_NUMBER", NULL, "atomic number", NULL));
    NC_CHECK(nc_defvar(ncid, &vid_extrapoints, "extra_point", NC_INT, 0, NULL,
                        "NUMEXTRA", NULL, NULL, NULL));

    int vid_mass, vid_radii, vid_screen;
    NC_CHECK(nc_defvar(ncid, &vid_mass, "atom_mass", NC_DOUBLE, 1, &dimid_atoms,
                        "MASS", "g/mol", NULL, NULL));
    NC_CHECK(nc_defvar(ncid, &vid_radii, "atom_GB_radius", NC_DOUBLE, 1, &dimid_atoms,
                        "RADII", "angstrom", "GB radius", NULL));
    NC_CHECK(nc_defvar(ncid, &vid_screen, "atom_GB_screen", NC_DOUBLE, 1, &dimid_atoms,
                        "SCREEN", NULL, "GB screening factor", NULL));

    NC_CHECK(nc_defvar(ncid, &vid_itype, "atom_nonbond_type_index", NC_INT, 1, &dimid_atoms,
                        "ATOM_TYPE_INDEX", NULL, "nonbond type index",
                        "reference_equation", "nonbond_lower_triangle_echelon_indexing",
                        "index_dimension", "nonbond_atom_type",
                        NULL));

    /* excluded atom list.
       n_atom_nonbond_exclusions (amber_name NUMBER_EXCLUDED_ATOMS) is a
       per-atom COUNT, dimensioned over n_atoms.
       n_nonbond_excluded_atoms (amber_name NNB) is the dimension size of
       the flattened EXCLUDED_ATOMS_LIST array itself, i.e. the sum of
       all the per-atom counts; nonbond_excluded_atom is dimensioned
       over it and carries the amber_name NNB since NNB in the original
       format IS that array's declared size.                              */
    int vid_nbex_count, vid_nbex_list;
    NC_CHECK(nc_defvar(ncid, &vid_nbex_count, "atom_nonbond_exclusion", NC_INT, 1, &dimid_atoms,
                        "NUMBER_EXCLUDED_ATOMS", NULL, "excluded atom count per atom",
                        "index_variables", "nonbond_excluded_atom",
                        NULL));
    {
        int dimid_excl;
        int n = iVarArrayElementCount(vaExcludedAtoms);
        NC_DIM(ncid, "nonbond_excluded_atom", "NNB", n, &dimid_excl);
        NC_CHECK(nc_defvar(ncid, &vid_nbex_list, "nonbond_excluded_atom", NC_INT, 1, &dimid_excl,
                            "EXCLUDED_ATOMS_LIST", NULL, "excluded atom list",
                            "encoding", "cumulative_sum_offset by n_atom_nonbond_exclusions",
                            "index_dimension", "atom",
                            NULL));
    }

    /* nonbonded parm index + LJ coefficients */
    int vid_acoef, vid_bcoef, vid_14acoef, vid_14bcoef, vid_nbidx;
    int nttyp  = iVarArrayElementCount(vaNBParameters);
    {
        int dimid_nb, dimid_nttyp;
        int ntypes = (int)((sqrt(8.0*nttyp + 1.0) - 1.0) / 2.0);
        assert(ntypes * ntypes == iVarArrayElementCount(vaNBIndexMatrix));
        NC_CHECK(nc_def_dim(ncid, "nonbond_atom_type", ntypes, &dimid_nb));
        int dims2[2] = {dimid_nb, dimid_nb};
        NC_CHECK(nc_defvar(ncid, &vid_nbidx, "nonbond_pair_index", NC_INT, 2, dims2,
                            "ATOM_TYPE_INDEX", NULL, "nonbond parm index",
                            "comment", "sign encodes table: positive = LJ echelon index, negative = hbond index",
                            "reference_equation", "nonbond_lower_triangle_echelon_indexing",
                            "index_dimension", "nonbond_pair_types,-n_hbond_type",
                            NULL));
        NC_CHECK(nc_def_dim(ncid, "nonbond_pair_type", nttyp, &dimid_nttyp));
        NC_CHECK(nc_defvar(ncid, &vid_acoef, "lennard_jones_acoef", NC_DOUBLE, 1, &dimid_nttyp,
                            "LENNARD_JONES_ACOEF", "kcal/mol*angstrom^12",
                            "LJ parm A coefficients, size NTYPES*(NTYPES+1)/2",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_bcoef, "lennard_jones_bcoef", NC_DOUBLE, 1, &dimid_nttyp,
                            "LENNARD_JONES_BCOEF", "kcal/mol*angstrom^6",
                            "LJ parm B coefficients, size NTYPES*(NTYPES+1)/2",
                            NULL));
        if (GDefaults.bCharmm) {
            NC_CHECK(nc_defvar(ncid, &vid_14acoef, "lennard_jones_14_acoef", NC_DOUBLE, 1, &dimid_nttyp,
                                "LENNARD_JONES_14_ACOEF", "kcal/mol*angstrom^12", NULL, NULL));
            NC_CHECK(nc_defvar(ncid, &vid_14bcoef, "lennard_jones_14_bcoef", NC_DOUBLE, 1, &dimid_nttyp,
                                "LENNARD_JONES_14_BCOEF", "kcal/mol*angstrom^6", NULL, NULL));
        }
    }

    /* residue properties */
    int vid_resname, vid_res_type, vid_res_first_atom, vid_res_max_atoms;
    {
        int dims2[2] = {dimid_res, dimid_name4};
        NC_CHECK(nc_defvar(ncid, &vid_resname, "residue_label", NC_CHAR, 2, dims2,
                            "RESIDUE_LABEL", NULL, "residue name",
                            "pdb_name", "resName",
                            "mmcif_name", GDefaults.bCIFReadAuth ? "atom_site.auth_comp_id"
                                                                   : "atom_site.label_comp_id",
                            NULL));
        dims2[1] = dimid_restype_len;
        NC_CHECK(nc_defvar(ncid, &vid_res_type, "residue_type", NC_CHAR, 2, dims2,
                            "RESIDUE_TYPE", NULL,
                            "residue type: w=water, p=protein, n=nucleic, s=saccharide,"
                            " i=ion, l=ligand, o=other, ?=undefined",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_res_max_atoms, "residue_max_atoms", NC_INT, 0, NULL,
                            "NMXRS", NULL, "maximum number of atoms per residue", NULL));
    }
    NC_CHECK(nc_defvar(ncid, &vid_res_first_atom, "residue_first_atom", NC_INT, 1, &dimid_res,
                        "RESIDUE_POINTER", NULL, "first atom index per residue",
                        "index_dimension", "atom",
                        NULL));

    /* bond parameters */
    int vid_bond_force, vid_bond_eq;
    int vid_bond_stiffness_pull_adj, vid_bond_equil_pull_adj;
    int vid_bond_stiffness_press_adj, vid_bond_equil_press_adj;
    {
        int nb = iParmSetTotalBondParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
        NC_CHECK(nc_def_dim(ncid, "bond_type", nb + nr, &dimid_bp));
        NC_CHECK(nc_defvar(ncid, &vid_bond_force, "bond_force_constant", NC_DOUBLE, 1, &dimid_bp,
                            "BOND_FORCE_CONSTANT", "kcal/mol/angstrom^2", NULL, NULL));
        NC_CHECK(nc_defvar(ncid, &vid_bond_eq, "bond_equil_value", NC_DOUBLE, 1, &dimid_bp,
                            "BOND_EQUIL_VALUE", "angstrom", NULL, NULL));
        if (bBondAug) {
            NC_CHECK(nc_defvar(ncid, &vid_bond_stiffness_pull_adj, "bond_stiffness_pull_adj", NC_DOUBLE, 1, &dimid_bp,
                                "BOND_STIFFNESS_PULL_ADJ", "kcal/mol/angstrom^2", NULL, NULL));
            NC_CHECK(nc_defvar(ncid, &vid_bond_equil_pull_adj, "bond_equil_pull_adj", NC_DOUBLE, 1, &dimid_bp,
                                "BOND_EQUIL_PULL_ADJ", "angstrom", NULL, NULL));
            NC_CHECK(nc_defvar(ncid, &vid_bond_stiffness_press_adj, "bond_stiffness_press_adj", NC_DOUBLE, 1, &dimid_bp,
                                "BOND_STIFFNESS_PRESS_ADJ", "kcal/mol/angstrom^2", NULL, NULL));
            NC_CHECK(nc_defvar(ncid, &vid_bond_equil_press_adj, "bond_equil_press_adj", NC_DOUBLE, 1, &dimid_bp,
                                "BOND_EQUIL_PRESS_ADJ", "angstrom", NULL, NULL));
        }
    }

    /* angle parameters */
    int vid_angle_force, vid_angle_eq;
    {
        int na = iParmSetTotalAngleParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
        NC_CHECK(nc_def_dim(ncid, "angle_type", na + nr, &dimid_ap));
        NC_CHECK(nc_defvar(ncid, &vid_angle_force, "angle_force_constant", NC_DOUBLE, 1, &dimid_ap,
                            "ANGLE_FORCE_CONSTANT", "kcal/mol/radian^2", NULL, NULL));
        NC_CHECK(nc_defvar(ncid, &vid_angle_eq, "angle_equil_value", NC_DOUBLE, 1, &dimid_ap,
                            "ANGLE_EQUIL_VALUE", "degree", NULL,
                            "comment", "convert to radians for AMBER internal units", NULL));
    }

    /* CHARMM UB terms */
    VARARRAY ub_force = NULL, ub_equil = NULL;
    VARARRAY ub_atoms = NULL, ub_index = NULL;
    int dimid_UB = 0, dimid_UB_types = 0;
    int vid_UB13_atoms, vid_UB13_index, vid_UB13_force, vid_UB13_equil;
    if (GDefaults.bCharmm) {
        int na = iParmSetTotalAngleParms(uUnit->psParameters);
        int *ub_index_map = CALLOC(sizeof(int) * na);
        ub_force = vaVarArrayCreate(sizeof(double));
        ub_equil = vaVarArrayCreate(sizeof(double));
        for (int i = 0; i < na; i++) {
            double dTkub, dRkub;
            ParmSetAngle(uUnit->psParameters, i, NULL,NULL,NULL,NULL,NULL, &dTkub, &dRkub, NULL);
            if (!dTkub) continue;
            VarArrayAdd(ub_force, &dTkub);
            VarArrayAdd(ub_equil, &dRkub);
            ub_index_map[i] = iVarArrayElementCount(ub_force);
        }
        ub_atoms = vaVarArrayCreate(sizeof(int) * 2);
        ub_index = vaVarArrayCreate(sizeof(int));
        na = iVarArrayElementCount(uUnit->vaAngles);
        for (int i = 0; i < na; i++) {
            SAVEANGLEt *saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
            if (ub_index_map[saPAngle->iParmIndex]) {
                int pair[2] = {saPAngle->iAtom1, saPAngle->iAtom3};
                VarArrayAdd(ub_atoms, pair);
                VarArrayAdd(ub_index, &ub_index_map[saPAngle->iParmIndex]);
            }
        }
        FREE(ub_index_map);
        if (iVarArrayElementCount(ub_index)) {
            NC_CHECK(nc_def_dim(ncid, "urey_bradley", iVarArrayElementCount(ub_index), &dimid_UB));
            NC_CHECK(nc_def_dim(ncid, "urey_bradley_type", iVarArrayElementCount(ub_force), &dimid_UB_types));
            int dims2[2] = {dimid_UB, dimid_atom2};
            NC_CHECK(nc_defvar(ncid, &vid_UB13_atoms, "urey_bradley_atom", NC_INT, 2, dims2,
                                NULL, NULL, "angle atom indices (atoms 1 and 3)",
                                "index_dimension", "atom",
                                NULL));
            NC_CHECK(nc_defvar(ncid, &vid_UB13_index, "urey_bradley_type_index", NC_INT, 1, &dimid_UB,
                                NULL, NULL, "Urey-Bradley parameter index",
                                "index_dimension", "urey_bradley_type",
                                NULL));
            NC_CHECK(nc_defvar(ncid, &vid_UB13_force, "urey_bradley_force_constant", NC_DOUBLE, 1, &dimid_UB_types,
                                "CHARMM_UREY_BRADLEY_FORCE_CONSTANT", "kcal/mol/angstrom^2", NULL, NULL));
            NC_CHECK(nc_defvar(ncid, &vid_UB13_equil, "urey_bradley_equil_value", NC_DOUBLE, 1, &dimid_UB_types,
                                "CHARMM_UREY_BRADLEY_EQUIL_VALUE", "angstrom", NULL, NULL));
        }
    }

    /* torsion parameters */
    int vid_dihe_force, vid_dihe_period, vid_dihe_phase, vid_dihe_scee, vid_dihe_scnb;
    int vid_harm_dihe_force, vid_harm_dihe_equil;
    {
        int dimid_tp;
        NC_CHECK(nc_def_dim(ncid, "dihedral_type", nd_parm_total, &dimid_tp));
        NC_CHECK(nc_defvar(ncid, &vid_dihe_force, "dihedral_force_constant", NC_DOUBLE, 1, &dimid_tp,
                            "DIHEDRAL_FORCE_CONSTANT", "kcal/mol", NULL, NULL));
        NC_CHECK(nc_defvar(ncid, &vid_dihe_period, "dihedral_periodicity", NC_INT, 1, &dimid_tp,
                            "DIHEDRAL_PERIODICITY", NULL, "torsion Fourier periodicity",
                            "comment", "term is cos(n*phi - phase)",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_dihe_phase, "dihedral_phase", NC_DOUBLE, 1, &dimid_tp,
                            "DIHEDRAL_PHASE", "degree", NULL,
                            "comment", "convert to radians for AMBER internal units", NULL));
        NC_CHECK(nc_defvar(ncid, &vid_dihe_scee, "dihedral_scee_scale", NC_DOUBLE, 1, &dimid_tp,
                            "SCEE_SCALE_FACTOR", NULL, "1-4 electrostatic energy scale factor", NULL));
        NC_CHECK(nc_defvar(ncid, &vid_dihe_scnb, "dihedral_scnb_scale", NC_DOUBLE, 1, &dimid_tp,
                            "SCNB_SCALE_FACTOR", NULL, "1-4 VdW energy scale factor", NULL));
        if (nh_parm_total) {
            int dimid_tp_harm;
            NC_CHECK(nc_def_dim(ncid, "harmonic_dihedral_type", nh_parm_total, &dimid_tp_harm));
            NC_CHECK(nc_defvar(ncid, &vid_harm_dihe_force, "harmonic_dihedral_force_constant", NC_DOUBLE, 1, &dimid_tp_harm,
                                "CHARMM_IMPROPER_FORCE_CONSTANT", "kcal/mol/radian^2",
                                "CHARMM harmonic improper force constant",
                                NULL));
            NC_CHECK(nc_defvar(ncid, &vid_harm_dihe_equil, "harmonic_dihedral_equil_value", NC_DOUBLE, 1, &dimid_tp_harm,
                                "CHARMM_IMPROPER_EQUIL_VALUE", "degree",
                                "CHARMM harmonic improper equilibrium angle",
                                NULL));
        }
    }

    /* unified bond interactions.
       Splitting bond_atom from bond_type_index (rather than keeping
       AMBER's packed atom_i,atom_j,parm_index stride-3 layout) is not a
       NetCDF-driven departure: AMBER's own PERT_BOND_ATOMS/PERT_BOND_PARAMS
       already split atoms from parm index for perturbation bonds. This
       adopts that convention uniformly instead of keeping two
       inconsistent AMBER layouts (packed for regular bonds, split for
       perturbed ones).                                                  */
    int vid_bond_atom, vid_bond_parm, vid_bond_restraints = 0, vid_bonds_inc_h, vid_bonds_excl_h, vid_NBONA;
    {
        int dimid_bond;
        NC_CHECK(nc_def_dim(ncid, "bond", n_bond_total, &dimid_bond));
        NC_CHECK(nc_defvar(ncid, &vid_bonds_inc_h, "bonds_inc_hydrogen", NC_INT, 0, NULL,
                            "NBONH", NULL, "count of leading bond_atom entries containing hydrogen",
                            "index_dimension", "bond",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_bonds_excl_h, "bonds_without_hydrogen", NC_INT, 0, NULL,
                            "MBONA", NULL, "count of trailing bond_atom entries without hydrogen",
                            "index_dimension", "bond",
                            "comment", "bonds_inc_hydrogen + n_bonds_without_hydrogen + n_bond_restraints == n_bond",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_bond_restraints, "bond_restraint", NC_INT, 0, NULL,
                            NULL, NULL, "count of bond restraints",
                            "index_dimension", "bond",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_NBONA, "NBONA", NC_INT, 0, NULL,
                            NULL, NULL, "count of bonds without hydrogens plus bond restraints",
                            "index_dimension", "bond",
                            NULL));
        int dims2[2] = {dimid_bond, dimid_atom2};
        NC_CHECK(nc_defvar(ncid, &vid_bond_atom, "bond_atom", NC_INT, 2, dims2,
                            "PERT_BOND_ATOMS", NULL, "atom indices",
                            "index_dimension", "atom",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_bond_parm, "bond_type_index", NC_INT, 1, &dimid_bond,
                            "PERT_BOND_PARAMS", NULL, "parm index",
                            "index_dimension", "bond_type",
                            "index_variables", "bond_force_constant bond_equil_value",
                            NULL));
    }

    /* unified angle interactions */
    int vid_angle_atom, vid_angle_parm, vid_angle_restraints = 0, vid_angles_inc_h, vid_angles_excl_h, vid_NPHIA;
    {
        int dimid_angle;
        NC_CHECK(nc_def_dim(ncid, "angle", n_angle_total, &dimid_angle));
        NC_CHECK(nc_defvar(ncid, &vid_angles_inc_h, "angles_inc_hydrogen", NC_INT, 0, NULL,
                            "NTHETH", NULL, "count of leading angles containing hydrogen",
                            "index_dimension", "angle",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_angles_excl_h, "angles_without_hydrogen", NC_INT, 0, NULL,
                            "MTHETA", NULL, "count of angles without hydrogen",
                            "index_dimension", "angle",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_angle_restraints, "angle_restraint", NC_INT, 0, NULL,
                            NULL, NULL, "count of angle restraints",
                            "index_dimension", "angle",
                            "comment", "angles_inc_hydrogen + n_angles_without_hydrogen + n_angle_restraints == n_angle",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_NPHIA, "NPHIA", NC_INT, 0, NULL,
                            NULL, NULL, "count of angles without hydrogen plus angle restraints",
                            "index_dimension", "angle",
                            "comment", "angles_inc_hydrogen + n_angles_without_hydrogen + n_angle_restraints == n_angle",
                            NULL));
        int dims2[2] = {dimid_angle, dimid_atom3};
        NC_CHECK(nc_defvar(ncid, &vid_angle_atom, "angle_atom", NC_INT, 2, dims2,
                            NULL, NULL, "atom indices",
                            "index_dimension", "atom",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_angle_parm, "angle_type_index", NC_INT, 1, &dimid_angle,
                            NULL, NULL, "parm index",
                            "index_variables", "angle_force_constant angle_equil_value",
                            NULL));
        NC_ATT(ncid, vid_angle_parm, "index_variables", "angle_type");
    }

    /* unified periodic dihedral interactions */
    int vid_dihe_atom, vid_dihe_parm, vid_dihe_proper, vid_dihe_term, vid_dihe_restraints, vid_dihe_inc_h, vid_dihe_excl_h, vid_NTHETA;
    {
        int dimid_dihe;
        NC_CHECK(nc_def_dim(ncid, "dihedral", n_dihe_total, &dimid_dihe));
        NC_CHECK(nc_defvar(ncid, &vid_dihe_inc_h, "dihedrals_inc_hydrogen", NC_INT, 0, NULL,
                            "NPHIH", NULL, "count of leading dihedrals containing hydrogen",
                            "index_dimension", "dihedral",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_dihe_excl_h, "dihedrals_without_hydrogen", NC_INT, 0, NULL,
                            "MPHIA", NULL, "count of dihedrals without hydrogen",
                            "index_dimension", "dihedral",
                            "comment", "dihedrals_inc_hydrogen + n_dihedrals_without_hydrogen == n_dihedral",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_dihe_restraints, "dihedral_restraint", NC_INT, 0, NULL,
                            NULL, NULL, "count of trailing dihedral restraints",
                            "index_dimension", "dihedral",
                            "comment", "dihedrals_inc_hydrogen + n_dihedrals_without_hydrogen == n_dihedral",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_NTHETA, "NPHIA", NC_INT, 0, NULL,
                            NULL, NULL, "count of dihedral_atom entries without hydrogen",
                            "index_dimension", "dihedral",
                            "comment", "dihedrals_inc_hydrogen + n_dihedrals_without_hydrogen == n_dihedral",
                            NULL));
        int dims2[2] = {dimid_dihe, dimid_atom4};
        NC_CHECK(nc_defvar(ncid, &vid_dihe_atom, "dihedral_atom", NC_INT, 2, dims2,
                            NULL, NULL, "atom indices",
                            "index_dimension", "atom",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_dihe_parm, "dihedral_type_index", NC_INT, 1, &dimid_dihe,
                            NULL, NULL, "parm index",
                            "index_dimension", "dihedral_type",
                            "index_variables", "dihedral_force_constant dihedral_periodicity dihedral_phase",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_dihe_proper, "dihedral_proper", NC_BYTE, 1, &dimid_dihe,
                            NULL, NULL, "proper-torsion flag",
                            "comment", "0 = improper, 1 = proper",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_dihe_term, "dihedral_term", NC_INT, 1, &dimid_dihe,
                            NULL, NULL, "multi-term Fourier series ordinal (0=solo, else 1,2,3,... within a set)",
                            "comment", "calc14 = dihedral_term < 2 && dihedral_proper",
                            NULL));
    }

    /* unified harmonic dihedral interactions.
       Not H-sorted, unlike bonds/angles/periodic dihedrals: harmonic
       (CHARMM improper) terms are always solo, there is no legacy
       fixed-position reader depending on their ordering, and sorting
       would add bookkeeping for no compatibility benefit.               */
    int vid_harm_atom = 0, vid_harm_parm = 0, vid_harm_restraint = 0;
    if (n_harm_total > 0) {
        int dimid_harm;
        NC_CHECK(nc_def_dim(ncid, "harmonic_dihedral", n_harm_total, &dimid_harm));
        int dims2[2] = {dimid_harm, dimid_atom4};
        NC_CHECK(nc_defvar(ncid, &vid_harm_atom, "harmonic_dihedral_atom", NC_INT, 2, dims2,
                            NULL, NULL, "harmonic dihedral atom indices",
                            "index_dimension", "atom",
                            "comment", "typical use: CHARMM improper, first atom = center atom",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_harm_parm, "harmonic_dihedral_type_index", NC_INT, 1, &dimid_harm,
                            NULL, NULL, "harmonic dihedral parameter index",
                            "index_variables", "harmonic_dihedral_force_constant harmonic_dihedral_equil_value",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_harm_restraint, "harmonic_dihedral_restraint", NC_BYTE, 1, &dimid_harm,
                            NULL, NULL, "supplemental restraint flag",
                            "flag_meanings", "force_field_term restraint",
                            NULL));
         int flag_val[2]={0,1};
         NC_CHECK(nc_put_att_int(ncid, vid_harm_restraint, "flag_values", NC_INT, 2, flag_val));
    }

    /* NATYP */
    NC_CHECK(nc_def_dim(ncid, "amber_type", iParmSetTotalAtomParms(uUnit->psParameters), &dimid_atyp));

    /* hbond */
    int vid_hba = 0, vid_hbb = 0;
    {
        int n = iParmSetTotalHBondParms(uUnit->psParameters);
        if (n > 0) {
            int dimid_hb;
            NC_CHECK(nc_def_dim(ncid, "hbond_type", n, &dimid_hb));
            NC_CHECK(nc_defvar(ncid, &vid_hba, "hbond_acoef", NC_DOUBLE, 1, &dimid_hb,
                                "HBOND_ACOEF", "kcal/mol*angstrom^12", NULL, NULL));
            NC_CHECK(nc_defvar(ncid, &vid_hbb, "hbond_bcoef", NC_DOUBLE, 1, &dimid_hb,
                                "HBOND_BCOEF", "kcal/mol*angstrom^10", NULL, NULL));
        }
    }

    /* atom type and tree strings */
    int vid_atom_type, vid_atom_tree;
    {
        int dims2[2] = {dimid_atoms, dimid_name4};
        NC_CHECK(nc_defvar(ncid, &vid_atom_type, "atom_forcefield_type", NC_CHAR, 2, dims2,
                            "AMBER_ATOM_TYPE", NULL, "force field atom type name", NULL));
        NC_CHECK(nc_defvar(ncid, &vid_atom_tree, "atom_tree_chain_classification", NC_CHAR, 2, dims2,
                            "TREE_CHAIN_CLASSIFICATION", NULL,
                            "mainchain tree class: M,E,S,B,3,4,5,6,X; BLA=unknown",
                            NULL));
    }

    /* solvent/box */
    int vid_mol = 0, vid_fsr = 0, vid_fsm = 0, vid_res_mol = 0;
    int nMol = 0, iFS1 = 0, iFirstSolvRes = 0;
    if (bUnitUseBox(uUnit)) {
        iFirstSolvRes = n_residues;
        for (int i = 0; i < n_residues; i++) {
            if (bResidueFlagsSet(PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->rResidue,
                                 RESIDUEBULKSOLVENT)) {
                iFirstSolvRes = i; break;
            }
        }
        UnitIOFindAndCountMolecules(uUnit);
        nMol = iVarArrayElementCount(uUnit->vaAtomsPerMolecule);
        iFS1 = uUnit->iFirstSolvent + 1;
        int dimid_mol;
        NC_CHECK(nc_def_dim(ncid, "molecule", nMol, &dimid_mol));
        NC_CHECK(nc_defvar(ncid, &vid_mol, "atoms_per_molecule", NC_INT, 1, &dimid_mol,
                            "NSPM", NULL, "number of atoms in molecule N", NULL));
        NC_CHECK(nc_defvar(ncid, &vid_fsr, "first_solvent_residue", NC_INT, 0, NULL,
                            "IPTRES", NULL, "residue index of first solvent residue",
                            "index_dimension", "residue",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_fsm, "first_solvent_molecule", NC_INT, 0, NULL,
                            "NSPSOL", NULL, "molecule index of first solvent molecule",
                            "index_dimension", "molecule",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_res_mol, "residue_molecule", NC_INT, 1, &dimid_res,
                            "RESIDUE_MOLECULE", NULL, "molecule index",
                            "index_dimension", "molecule",
                            NULL));
    }

    /* cap info */
    int vid_cap_info_atom = 0, vid_cap_info_radius = 0, vid_cap_info_origin = 0;
    if (bUnitUseSolventCap(uUnit)) {
        int dimid_vec;
        NC_CHECK(nc_defvar(ncid, &vid_cap_info_atom, "cap_info_atom", NC_INT, 0, NULL,
                            "CAP_INFO", NULL, "index of the last atom not in the solvent cap",
                            "index_dimension", "atom",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_cap_info_radius, "cap_info_radius", NC_DOUBLE, 0, NULL,
                            "CAP_INFO2", "angstrom", "radius of solvent cap", NULL));
        NC_CHECK(nc_def_dim(ncid, "vector", 3, &dimid_vec));
        NC_CHECK(nc_defvar(ncid, &vid_cap_info_origin, "cap_info_origin", NC_DOUBLE, 1, &dimid_vec,
                            "CAP_INFO2", "angstrom", "origin of solvent cap", NULL));
    }

    /* GB radius set */
    {
        STRING sTemp;
        const char *rset = "Unknown radius set";
        if (GDefaults.iGBparm >= 0 && GDefaults.iGBparm <= GBPARM_MAX) {
            sprintf(sTemp, "%s (%s)", PBRadii_optionDesc[GDefaults.iGBparm],
                    PBRadii_options[GDefaults.iGBparm]);
            rset = sTemp;
        }
        NC_ATT(ncid, NC_GLOBAL, "radius_set", rset);
    }

    int vid_IPOL;
    NC_CHECK(nc_defvar(ncid, &vid_IPOL, "polarizability_type", NC_INT, 0, NULL,
                        "IPOL", NULL, "polarizable force field flag", NULL));

    int vid_polar = 0, vid_damp = 0;
    if (bPolar) {
        NC_CHECK(nc_defvar(ncid, &vid_polar, "polarizability", NC_DOUBLE, 1, &dimid_atoms,
                            "POLARIZABILITY", "angstrom^3", NULL, NULL));
        NC_CHECK(nc_defvar(ncid, &vid_damp, "dipole_damp_factor", NC_DOUBLE, 1, &dimid_atoms,
                            "DIPOLE_DAMP_FACTOR", NULL, NULL, NULL));
    }

    /* PDB chain info */
    int vid_resSeq = 0, vid_res_chain = 0, vid_res_icode = 0;
    bool bHasICodes = false;
    if (GDefaults.bPdbKeepChainId) {
        NC_CHECK(nc_defvar(ncid, &vid_resSeq, "residue_number", NC_INT, 1, &dimid_res,
                            "RESIDUE_NUMBER", NULL, "PDB resSeq",
                            "pdb_name", "resSeq",
                            "mmcif_name", GDefaults.bCIFReadAuth ? "atom_site.auth_seq_id"
                                                                   : "atom_site.label_seq_id",
                            NULL));
        int dims2[2] = {dimid_res, dimid_chainid_len};
        NC_CHECK(nc_defvar(ncid, &vid_res_chain, "residue_chainid", NC_CHAR, 2, dims2,
                            "RESIDUE_CHAINID", NULL, "PDB chainId",
                            "pdb_name", "chainId",
                            "mmcif_name", GDefaults.bCIFReadAuth ? "atom_site.auth_asym_id"
                                                                   : "atom_site.label_asym_id",
                            NULL));
        for (int i = 0; i < n_residues; i++) {
            if (PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sICode[0] != ' ') {
                bHasICodes = true; break;
            }
        }
        if (bHasICodes) {
            dims2[1] = dimid_icode_len;
            NC_CHECK(nc_defvar(ncid, &vid_res_icode, "residue_icode", NC_CHAR, 2, dims2,
                                "RESIDUE_ICODE", NULL, "PDB iCode",
                                "pdb_name", "iCode",
                                "mmcif_name", "atom_site.pdbx_PDB_ins_code",
                                NULL));
        }
    }

    /* CMAP */
    int vid_cmap_atoms = 0, vid_cmap_type_index = 0;
    int vid_cmap_resolution = 0;
    int vid_cmap_energy = 0, vid_cmap_title = 0;
    int cmap_energy_size = 0, cmap_title_len = 0;
    if (CMAP_TYPES && CMAP_KINDS && GDefaults.bCMAP) {
        for (int i = 0; i < CMAP_TYPES; i++) {
            CMAPt cmap;
            ParmSetCMAP(uUnit->psParameters, i, &cmap, false);
            int len = (int)strlen(cmap.title);
            if (len > cmap_title_len) cmap_title_len = len;
            cmap_energy_size += cmap.resolution * cmap.resolution;
        }
        int dimid_types, dimid_kinds, dimid_values, dimid_cmap_title_len;
        NC_CHECK(nc_def_dim(ncid, "cmap_phi_psi", CMAP_KINDS, &dimid_kinds));
        NC_CHECK(nc_def_dim(ncid, "cmap_type", CMAP_TYPES, &dimid_types));
        NC_CHECK(nc_def_dim(ncid, "cmap_parameter", cmap_energy_size, &dimid_values));
        NC_CHECK(nc_def_dim(ncid, "cmap_title_length", cmap_title_len, &dimid_cmap_title_len));
        int dims2[2] = {dimid_types, dimid_cmap_title_len};
        NC_CHECK(nc_defvar(ncid, &vid_cmap_title, "cmap_title", NC_CHAR, 2, dims2,
                            NULL, NULL, "descriptive title per CMAP type", NULL));
        NC_CHECK(nc_defvar(ncid, &vid_cmap_resolution, "cmap_resolution", NC_INT, 1, &dimid_types,
                            "CMAP_RESOLUTION", NULL, "grid resolution per map type",
                            "index_dimension", "cmap_type",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_cmap_energy, "cmap_parameter", NC_DOUBLE, 1, &dimid_values,
                            NULL, "kcal/mol", "Concatenated phi/psi torsion energy correction grids",
                            "row_axis", "phi",
                            "column_axis", "psi",
                            "axis_range", "[-180,+180) degrees",
                            "axis_units", "degrees",
                            "packing", "concatenated row-major grids; offset for type i is sum(cmap_resolution(1:i-1)**2)",
                            NULL));
        dims2[0] = dimid_kinds; dims2[1] = dimid_atom5;
        NC_CHECK(nc_defvar(ncid, &vid_cmap_atoms, "cmap_phi_psi_atom", NC_INT, 2, dims2,
                            "CMAP_INDEX", NULL, "5 atom indices: atoms 1-4=phi, 2-5=psi",
                            "index_dimension", "atom",
                            NULL));
        NC_CHECK(nc_defvar(ncid, &vid_cmap_type_index, "cmap_phi_psi_type_index", NC_INT, 1, &dimid_kinds,
                            "CMAP_INDEX", NULL, "index into cmap_title/cmap_resolution/cmap_energy",
                            "index_dimension", "cmap_type",
                            "index_variables", "cmap_title cmap_resolution cmap_energy",
                            NULL));
    }

    NC_CHECK(nc_enddef(ncid));

    for (int i=0;i<nparam;i++) {
      NC_CHECK(nc_put_var_int(ncid, param_vid[i], &param_val[i]));
    }

    /* write pre-built torsion arrays, free immediately */
    NC_CHECK(nc_put_var_int  (ncid, vid_dihe_inc_h, &dihe_nh));
    NC_CHECK(nc_put_var_int  (ncid, vid_dihe_excl_h, &dihe_nnh));
    dihe_nnh += n_dihe_restraints;
    NC_CHECK(nc_put_var_int  (ncid, vid_NPHIA, &dihe_nnh));
    NC_CHECK(nc_put_var_int  (ncid, vid_dihe_restraints, &n_dihe_restraints));
    NC_CHECK(nc_put_var_int  (ncid, vid_dihe_atom, (int*)dihe_atoms));
    NC_CHECK(nc_put_var_int  (ncid, vid_dihe_parm, dihe_parm));
    NC_CHECK(nc_put_var_schar(ncid, vid_dihe_proper, (signed char*)dihe_proper));
    NC_CHECK(nc_put_var_int  (ncid, vid_dihe_term, dihe_term));
    NC_CHECK(nc_put_var_double(ncid, vid_dihe_force, dihe_force));
    NC_CHECK(nc_put_var_int   (ncid, vid_dihe_period, dihe_period));
    NC_CHECK(nc_put_var_double(ncid, vid_dihe_phase, dihe_phase));
    NC_CHECK(nc_put_var_double(ncid, vid_dihe_scee, dihe_scee));
    NC_CHECK(nc_put_var_double(ncid, vid_dihe_scnb, dihe_scnb));
    FREE(dihe_atoms); FREE(dihe_parm); FREE(dihe_proper); FREE(dihe_term);
    FREE(dihe_force); FREE(dihe_period); FREE(dihe_phase); FREE(dihe_scee);
    FREE(dihe_scnb);

    if (n_harm_total > 0) {
        NC_CHECK(nc_put_var_int   (ncid, vid_harm_atom, (int*)harm_atoms));
        NC_CHECK(nc_put_var_int   (ncid, vid_harm_parm, harm_parm));
        NC_CHECK(nc_put_var_schar (ncid, vid_harm_restraint, (signed char*)harm_restraint));
        NC_CHECK(nc_put_var_double(ncid, vid_harm_dihe_force, harm_force));
        NC_CHECK(nc_put_var_double(ncid, vid_harm_dihe_equil, harm_equil));
    }
    FREE(harm_atoms); FREE(harm_parm); FREE(harm_restraint);
    FREE(harm_force); FREE(harm_equil);

    /* ── sequential temporary buffers for all remaining data ── */

    /* parameter set names */
    {
        PARMSET psParm;
        char *buf = MALLOC(iParmSets * fname_length);
        memset(buf, ' ', iParmSets * fname_length);
        for (int i = 0; i < iParmSets; i++) {
            psParm = *PVAI(uUnit->vaParmSets, PARMSET, i);
            memcpy(buf + i*fname_length, psParm->sFname, strlen(psParm->sFname));
        }
        NC_CHECK(nc_put_var_text(ncid, vid_parm_file, buf));
        FREE(buf);
        buf = MALLOC(iParmSets * title_length);
        memset(buf, ' ', iParmSets * title_length);
        for (int i = 0; i < iParmSets; i++) {
            psParm = *PVAI(uUnit->vaParmSets, PARMSET, i);
            memcpy(buf + i*title_length, psParm->sTitle, strlen(psParm->sTitle));
        }
        NC_CHECK(nc_put_var_text(ncid, vid_parm_title, buf));
        FREE(buf);
        buf = MALLOC(iParmSets * hash_length);
        memset(buf, ' ', iParmSets * hash_length);
        for (int i = 0; i < iParmSets; i++) {
            psParm = *PVAI(uUnit->vaParmSets, PARMSET, i);
            memcpy(buf + i*hash_length, psParm->cksum, strlen(psParm->cksum));
        }
        NC_CHECK(nc_put_var_text(ncid, vid_parm_hash, buf));
        FREE(buf);
    }

    /* cell dimensions */
    if (bUnitUseBox(uUnit)) {
        double vector[3] = {
            uUnit->dAlpha * RADTODEG,
            uUnit->dBeta  * RADTODEG,
            uUnit->dGamma * RADTODEG
        };
        NC_CHECK(nc_put_var_double(ncid, vid_cell_angle, vector));
        UnitGetBox(uUnit, &vector[0], &vector[1], &vector[2]);
        NC_CHECK(nc_put_var_double(ncid, vid_cell_length, vector));
    }

    /* atom names */
    {
        char *buf = MALLOC(n_atoms * LABLEN);
        memset(buf, ' ', n_atoms * LABLEN);
        for (int i = 0; i < n_atoms; i++) {
            const char *s = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sName;
            size_t l = strlen(s);
            if (l > LABLEN) s += l - LABLEN;
            memcpy(buf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
        }
        NC_CHECK(nc_put_var_text(ncid, vid_name, buf));
        FREE(buf);
    }

    /* charges */
    {
        double *buf = MALLOC(n_atoms * sizeof(double));
        for (int i = 0; i < n_atoms; i++)
            buf[i] = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->dCharge;
        NC_CHECK(nc_put_var_double(ncid, vid_charge, buf));
        FREE(buf);
    }

    /* atomic numbers */
    {
        int *buf = MALLOC(n_atoms * sizeof(int));
        for (int i = 0; i < n_atoms; i++)
            buf[i] = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iElement;
        NC_CHECK(nc_put_var_int(ncid, vid_atomic_num, buf));
        FREE(buf);
    }

    /* mass, GB radii, screen: merged single ParmSetAtom pass */
    {
        double *mass   = MALLOC(n_atoms * sizeof(double));
        double *radii  = MALLOC(n_atoms * sizeof(double));
        double *screen = MALLOC(n_atoms * sizeof(double));
        char sType[MAXTYPELEN];
        double dMass, dPolar, dEpsilon, dRStar, dEpsilon14, dRStar14, dScreenF;
        int iElement, iHybridization;
        SAVEATOMt *sa = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
        for (int i = 0; i < n_atoms; i++) {
            int iIdx = iParmSetFindAtom(uUnit->psParameters, sa[i].sType);
            ParmSetAtom(uUnit->psParameters, iIdx, sType, &dMass, &dPolar,
                        &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
                        &iElement, &iHybridization, sDesc);
            mass[i]   = dMass;
            radii[i]  = dGBRadiusForAtom(&sa[i], iElement, dMass, (i+1 == n_atoms));
            screen[i] = dGBScreenForElement(iElement);
        }
        NC_CHECK(nc_put_var_double(ncid, vid_mass, mass));
        NC_CHECK(nc_put_var_double(ncid, vid_radii, radii));
        NC_CHECK(nc_put_var_double(ncid, vid_screen, screen));
        FREE(mass); FREE(radii); FREE(screen);
    }

    /* atom type index */
    {
        int *buf = MALLOC(n_atoms * sizeof(int));
        for (int i = 0; i < n_atoms; i++) {
            int iAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iTypeIndex - 1;
            buf[i] = *PVAI(vaNBIndex, int, iAtom) + 1;
        }
        NC_CHECK(nc_put_var_int(ncid, vid_itype, buf));
        FREE(buf);
    }

    /* excluded atom counts and list */
    NC_CHECK(nc_put_var_int(ncid, vid_nbex_count, PVAI(vaExcludedCount, int, 0)));
    NC_CHECK(nc_put_var_int(ncid, vid_nbex_list, PVAI(vaExcludedAtoms, int, 0)));

    /* nonbonded parm index matrix */
    NC_CHECK(nc_put_var_int(ncid, vid_nbidx, PVAI(vaNBIndexMatrix, int, 0)));

    /* residue labels and type */
    {
        char *buf      = MALLOC(n_residues * LABLEN);
        char *buf_type = MALLOC(n_residues);
        int iUnknownTypes = 0;
        memset(buf, ' ', n_residues * LABLEN);
        SAVERESIDUEt *srPRes = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        for (int i = 0; i < n_residues; i++, srPRes++) {
            const char *s = srPRes->sName;
            size_t l = strlen(s);
            if (l > 3) s += l - 3;
            memcpy(buf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
            buf_type[i] = srPRes->sResidueType[0];
            if (!buf_type[i]) buf_type[i] = '?';
            if (buf_type[i] == '?') iUnknownTypes++;
        }
        if (iUnknownTypes > 0)
            VPWARN("%d residues have undefined type\n", iUnknownTypes);
        NC_CHECK(nc_put_var_text(ncid, vid_resname, buf));
        NC_CHECK(nc_put_var_text(ncid, vid_res_type, buf_type));
        NC_CHECK(nc_put_var_int(ncid, vid_res_max_atoms, &iMaxAtoms));
        FREE(buf); FREE(buf_type);
    }

    /* residue pointer */
    {
        int *buf = MALLOC(n_residues * sizeof(int));
        for (int i = 0; i < n_residues; i++)
            buf[i] = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iAtomStartIndex;
        NC_CHECK(nc_put_var_int(ncid, vid_res_first_atom, buf));
        FREE(buf);
    }

    /* bond parameters */
    {
        int nb = iParmSetTotalBondParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
        int n  = nb + nr;
        double *kb = MALLOC(n*sizeof(double)), *r0 = MALLOC(n*sizeof(double));
        for (int i = 0; i < nb; i++) {
            ParmSetBond(uUnit->psParameters, i, s1, s2, &dKb, &dR0,
                        &dKpull, &dRpull0, &dKpress, &dRpress0, sDesc);
            kb[i] = dKb; r0[i] = dR0;
        }
        NC_RESTRAINTLOOP(uUnit, RESTRAINTBOND, dKx, kb, nb);
        NC_RESTRAINTLOOP(uUnit, RESTRAINTBOND, dX0, r0, nb);
        NC_CHECK(nc_put_var_double(ncid, vid_bond_force, kb));
        NC_CHECK(nc_put_var_double(ncid, vid_bond_eq, r0));

        /* bond augmentation */
        if (bBondAug) {
            double *kpull  = kb, *rpull0  = r0;
            double *kpress = MALLOC(n*sizeof(double)), *rpress0 = MALLOC(n*sizeof(double));
            for (int i = 0; i < nb; i++) {
                ParmSetBond(uUnit->psParameters, i, s1, s2, &dKb, &dR0,
                            &dKpull, &dRpull0, &dKpress, &dRpress0, sDesc);
                kpull[i] = dKpull; rpull0[i] = dRpull0;
                kpress[i] = dKpress; rpress0[i] = dRpress0;
            }
            for (int i = nb; i < n; i++) {
                kpull[i] = 0.0; rpull0[i] = 100.0;
                kpress[i] = 0.0; rpress0[i] = -100.0;
            }
            NC_CHECK(nc_put_var_double(ncid, vid_bond_stiffness_pull_adj, kpull));
            NC_CHECK(nc_put_var_double(ncid, vid_bond_equil_pull_adj, rpull0));
            NC_CHECK(nc_put_var_double(ncid, vid_bond_stiffness_press_adj, kpress));
            NC_CHECK(nc_put_var_double(ncid, vid_bond_equil_press_adj, rpress0));
            FREE(kpress); FREE(rpress0);
        }
        FREE(kb); FREE(r0);
    }

    /* angle parameters */
    {
        int na = iParmSetTotalAngleParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
        int n  = na + nr;
        double *kt = MALLOC(n*sizeof(double)), *t0 = MALLOC(n*sizeof(double));
        for (int i = 0; i < na; i++) {
            ParmSetAngle(uUnit->psParameters, i, NULL,NULL,NULL, &dKt, &dT0,
                         NULL, NULL, NULL);
            kt[i] = dKt; t0[i] = dT0 * RADTODEG;
        }
        NC_RESTRAINTLOOP(uUnit, RESTRAINTANGLE, dKx, kt, na);
        NC_RESTRAINTLOOP(uUnit, RESTRAINTANGLE, dX0, t0, na);
        NC_CHECK(nc_put_var_double(ncid, vid_angle_force, kt));
        NC_CHECK(nc_put_var_double(ncid, vid_angle_eq, t0));
        FREE(kt); FREE(t0);
    }

    /* CHARMM UB terms */
    if (GDefaults.bCharmm && ub_atoms && iVarArrayElementCount(ub_atoms)) {
        NC_CHECK(nc_put_var_int   (ncid, vid_UB13_atoms, PVAI(ub_atoms, int, 0)));
        NC_CHECK(nc_put_var_int   (ncid, vid_UB13_index, PVAI(ub_index, int, 0)));
        NC_CHECK(nc_put_var_double(ncid, vid_UB13_force, PVAI(ub_force, double, 0)));
        NC_CHECK(nc_put_var_double(ncid, vid_UB13_equil, PVAI(ub_equil, double, 0)));
    }
    if (GDefaults.bCharmm && ub_atoms) {
        VarArrayDestroy(&ub_atoms); VarArrayDestroy(&ub_index);
        VarArrayDestroy(&ub_force); VarArrayDestroy(&ub_equil);
    }

    /* LJ coefficients */
    {
        double *A = MALLOC(nttyp*sizeof(double)), *B = MALLOC(nttyp*sizeof(double));
        for (int i = 0; i < nttyp; i++) {
            A[i] = PVAI(vaNBParameters, NONBONDACt, i)->dA;
            B[i] = PVAI(vaNBParameters, NONBONDACt, i)->dB;
        }
        NC_CHECK(nc_put_var_double(ncid, vid_acoef, A));
        NC_CHECK(nc_put_var_double(ncid, vid_bcoef, B));
        if (GDefaults.bCharmm) {
            for (int i = 0; i < nttyp; i++) {
                A[i] = PVAI(vaNBParameters, NONBONDACt, i)->dA14;
                B[i] = PVAI(vaNBParameters, NONBONDACt, i)->dB14;
            }
            NC_CHECK(nc_put_var_double(ncid, vid_14acoef, A));
            NC_CHECK(nc_put_var_double(ncid, vid_14bcoef, B));
        }
        FREE(A); FREE(B);
    }

    /* unified bond interactions — temporary buffer, write, free */
    {
        int (*buf)[2] = MALLOC(n_bond_total * sizeof(*buf));
        int  *parm    = MALLOC(n_bond_total * sizeof(int));
        char *restr   = MALLOC(n_bond_total * sizeof(char));
        int nb = 0, nh, nnh;
        for (int i = 0; i < n_bonds; i++) {
            SAVEBONDt *b = PVAI(uUnit->vaBonds, SAVEBONDt, i);
            ATOM aA = PVAI(uUnit->vaAtoms, SAVEATOMt, b->iAtom1 - 1)->aAtom;
            ATOM aB = PVAI(uUnit->vaAtoms, SAVEATOMt, b->iAtom2 - 1)->aAtom;
            if ( iAtomElement(aA) == HYDROGEN || iAtomElement(aB) == HYDROGEN ) {
                buf[nb][0] = b->iAtom1; buf[nb][1] = b->iAtom2;
                parm[nb] = b->iParmIndex; restr[nb] = 0; nb++;
            }
        }
        nh = nb;
        for (int i = 0; i < n_bonds; i++) {
            SAVEBONDt *b = PVAI(uUnit->vaBonds, SAVEBONDt, i);
            ATOM aA = PVAI(uUnit->vaAtoms, SAVEATOMt, b->iAtom1 - 1)->aAtom;
            ATOM aB = PVAI(uUnit->vaAtoms, SAVEATOMt, b->iAtom2 - 1)->aAtom;
            if ( iAtomElement(aA) != HYDROGEN && iAtomElement(aB) != HYDROGEN ) {
                buf[nb][0] = b->iAtom1; buf[nb][1] = b->iAtom2;
                parm[nb] = b->iParmIndex; restr[nb] = 0; nb++;
            }
        }
        nnh = nb - nh;
        if (nr_bond) {
            SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
            for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++) {
                if (sr->iType != RESTRAINTBOND) continue;
                buf[nb][0] = sr->iAtom1; buf[nb][1] = sr->iAtom2;
                parm[nb] = sr->iParmIndex; restr[nb] = 1; nb++;
            }
        }
        NC_CHECK(nc_put_var_int  (ncid, vid_bonds_inc_h, &nh));
        NC_CHECK(nc_put_var_int  (ncid, vid_bonds_excl_h, &nnh));
        nnh += nr_bond;
        NC_CHECK(nc_put_var_int  (ncid, vid_NBONA, &nnh));
        NC_CHECK(nc_put_var_int  (ncid, vid_bond_restraints, &nr_bond));
        NC_CHECK(nc_put_var_int  (ncid, vid_bond_atom, (int*)buf));
        NC_CHECK(nc_put_var_int  (ncid, vid_bond_parm, parm));
        FREE(buf); FREE(parm); FREE(restr);
    }

    /* unified angle interactions — temporary buffer, write, free */
    {
        int (*buf)[3] = MALLOC(n_angle_total * sizeof(*buf));
        int  *parm    = MALLOC(n_angle_total * sizeof(int));
        int na = 0, nh, nnh;
        for (int i = 0; i < n_angles; i++) {
            SAVEANGLEt *a = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
            ATOM aA = PVAI(uUnit->vaAtoms, SAVEATOMt, a->iAtom1 - 1)->aAtom;
            ATOM aB = PVAI(uUnit->vaAtoms, SAVEATOMt, a->iAtom2 - 1)->aAtom;
            ATOM aC = PVAI(uUnit->vaAtoms, SAVEATOMt, a->iAtom3 - 1)->aAtom;
            if ( iAtomElement(aA) == HYDROGEN || iAtomElement(aB) == HYDROGEN || iAtomElement(aC) == HYDROGEN ) {
                buf[na][0] = a->iAtom1; buf[na][1] = a->iAtom2; buf[na][2] = a->iAtom3;
                parm[na] = a->iParmIndex; na++;
            }
        }
        nh = na;
        for (int i = 0; i < n_angles; i++) {
            SAVEANGLEt *a = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
            ATOM aA = PVAI(uUnit->vaAtoms, SAVEATOMt, a->iAtom1 - 1)->aAtom;
            ATOM aB = PVAI(uUnit->vaAtoms, SAVEATOMt, a->iAtom2 - 1)->aAtom;
            ATOM aC = PVAI(uUnit->vaAtoms, SAVEATOMt, a->iAtom3 - 1)->aAtom;
            if ( iAtomElement(aA) != HYDROGEN && iAtomElement(aB) != HYDROGEN && iAtomElement(aC) != HYDROGEN ) {
                buf[na][0] = a->iAtom1; buf[na][1] = a->iAtom2; buf[na][2] = a->iAtom3;
                parm[na] = a->iParmIndex; na++;
            }
        }
        nnh = na - nh;
        if (nr_angle) {
            SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
            for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++) {
                if (sr->iType != RESTRAINTANGLE) continue;
                buf[na][0] = sr->iAtom1; buf[na][1] = sr->iAtom2; buf[na][2] = sr->iAtom3;
                parm[na] = sr->iParmIndex; na++;
            }
        }
        NC_CHECK(nc_put_var_int  (ncid, vid_angles_inc_h, &nh));
        NC_CHECK(nc_put_var_int  (ncid, vid_angles_excl_h, &nnh));
        nnh += nr_angle;
        NC_CHECK(nc_put_var_int  (ncid, vid_NPHIA, &nnh));
        NC_CHECK(nc_put_var_int  (ncid, vid_angle_restraints, &nr_angle));
        NC_CHECK(nc_put_var_int  (ncid, vid_angle_atom, (int*)buf));
        NC_CHECK(nc_put_var_int  (ncid, vid_angle_parm, parm));
        FREE(buf); FREE(parm);
    }

    /* hbond */
    {
        int n = iParmSetTotalHBondParms(uUnit->psParameters);
        if (n > 0) {
            double *A = MALLOC(n*sizeof(double)), *B = MALLOC(n*sizeof(double));
            for (int i = 0; i < n; i++) {
                ParmSetHBond(uUnit->psParameters, i, s1, s2, &dC, &dD, sDesc);
                A[i] = dC; B[i] = dD;
            }
            NC_CHECK(nc_put_var_double(ncid, vid_hba, A));
            NC_CHECK(nc_put_var_double(ncid, vid_hbb, B));
            FREE(A); FREE(B);
        }
    }

    /* atom type and tree strings — reuse buffer */
    {
        char *buf = MALLOC(n_atoms * LABLEN);
        memset(buf, ' ', n_atoms * LABLEN);
        int n_extra_points = 0;
        for (int i = 0; i < n_atoms; i++) {
            const char *s = sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom);
            if (!strcmp(s,"EP")) n_extra_points ++;
            size_t l = strlen(s);
            if (l > LABLEN) s += l - LABLEN;
            memcpy(buf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
        }
        NC_CHECK(nc_put_var_int(ncid, vid_extrapoints, &n_extra_points));
        NC_CHECK(nc_put_var_text(ncid, vid_atom_type, buf));
        memset(buf, ' ', n_atoms * LABLEN);
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
        NC_CHECK(nc_put_var_text(ncid, vid_atom_tree, buf));
        FREE(buf);
    }

    /* solvent/box */
    if (bUnitUseBox(uUnit)) {
        NC_CHECK(nc_put_var_int(ncid, vid_fsr, &iFirstSolvRes));
        NC_CHECK(nc_put_var_int(ncid, vid_fsm, &iFS1));
        NC_CHECK(nc_put_var_int(ncid, vid_mol, PVAI(uUnit->vaAtomsPerMolecule, int, 0)));
    }

    /* cap info */
    if (bUnitUseSolventCap(uUnit)) {
        int ires = uUnit->iCapTempInt;
        double dX, dY, dZ, dR;
        UnitGetSolventCap(uUnit, &dX, &dY, &dZ, &dR);
        double cap[3] = {dX, dY, dZ};
        NC_CHECK(nc_put_var_int   (ncid, vid_cap_info_atom, &ires));
        NC_CHECK(nc_put_var_double(ncid, vid_cap_info_radius, &dR));
        NC_CHECK(nc_put_var_double(ncid, vid_cap_info_origin, cap));
    }

    /* IPOL */
    { int val = GDefaults.iIPOL; NC_CHECK(nc_put_var_int(ncid, vid_IPOL, &val)); }

    /* polar */
    if (bPolar) {
        int iCount = 0;
        double *polar = MALLOC(n_atoms * sizeof(double));
        double *damp  = MALLOC(n_atoms * sizeof(double));
        char sType[MAXTYPELEN];
        double dMass, dPolar, dEpsilon, dRStar, dEpsilon14, dRStar14, dScreenF;
        int iElement, iHybridization;
        SAVEATOMt *a = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
        for (int i = 0; i < n_atoms; i++, a++) {
            int iIdx = iParmSetFindAtom(uUnit->psParameters, a->sType);
            ParmSetAtom(uUnit->psParameters, iIdx, sType, &dMass, &dPolar,
                        &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
                        &iElement, &iHybridization, sDesc);
            if (dPolar == -1.0) { dPolar = 0.0; iCount++; }
            polar[i] = dPolar;
            damp[i]  = (dScreenF == 0.0) ? GDefaults.dDipoleDampFactor : dScreenF;
        }
        if (iCount > 0)
            VP0("Total atoms with default polarization=0.0: %d of %d\n", iCount, n_atoms);
        NC_CHECK(nc_put_var_double(ncid, vid_polar, polar));
        NC_CHECK(nc_put_var_double(ncid, vid_damp, damp));
        FREE(polar); FREE(damp);
    }

    /* PDB chain info */
    if (GDefaults.bPdbKeepChainId) {
        int *buf = MALLOC(n_residues * sizeof(int));
        for (int i = 0; i < n_residues; i++)
            buf[i] = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iPdbResSeq;
        NC_CHECK(nc_put_var_int(ncid, vid_resSeq, buf));
        FREE(buf);
        char *cbuf = MALLOC(n_residues * chainid_len);
        memset(cbuf, ' ', n_residues * chainid_len);
        for (int i = 0; i < n_residues; i++) {
            char *chainid = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sChainId;
            for (int j = 0; chainid[j]; j++)
                cbuf[i*chainid_len + j] = chainid[j];
        }
        NC_CHECK(nc_put_var_text(ncid, vid_res_chain, cbuf));
        FREE(cbuf);
        if (bHasICodes) {
            char *sbuf = MALLOC(n_residues);
            for (int i = 0; i < n_residues; i++) {
                char iCode = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sICode[0];
                sbuf[i] = iCode ? iCode : ' ';
            }
            NC_CHECK(nc_put_var_text(ncid, vid_res_icode, sbuf));
            FREE(sbuf);
        }
    }

    /* CMAP */
    if (CMAP_TYPES && CMAP_KINDS && GDefaults.bCMAP) {
        int *atoms_buf = MALLOC(CMAP_KINDS * 5 * sizeof(int));
        int *index_buf = MALLOC(CMAP_KINDS * sizeof(int));
        for (int i = 0, idx = 0; i < CMAP_KINDS; i++, idx += 5) {
            SAVECMAPt *saveCMAP = PVAI(uUnit->vaCMAPs, SAVECMAPt, i);
            for (int j = 0; j < 5; j++) atoms_buf[idx+j] = saveCMAP->iAtom[j];
            index_buf[i] = saveCMAP->iParmIndex;
        }
        NC_CHECK(nc_put_var_int(ncid, vid_cmap_atoms, atoms_buf));
        NC_CHECK(nc_put_var_int(ncid, vid_cmap_type_index, index_buf));
        FREE(atoms_buf); FREE(index_buf);
        double *data_buf   = MALLOC(cmap_energy_size * sizeof(double));
        char   *title_buf  = MALLOC(cmap_title_len * CMAP_TYPES);
        int    *res_buf    = MALLOC(sizeof(int) * CMAP_TYPES);
        int data_offset = 0;
        memset(title_buf, ' ', cmap_title_len * CMAP_TYPES);
        for (int i = 0; i < CMAP_TYPES; i++) {
            CMAPt cmap;
            ParmSetCMAP(uUnit->psParameters, i, &cmap, false);
            memcpy(title_buf + i*cmap_title_len, cmap.title, strlen(cmap.title));
            int m = cmap.resolution * cmap.resolution;
            memcpy(data_buf + data_offset, cmap.map, m * sizeof(double));
            res_buf[i] = cmap.resolution;
            data_offset += m;
        }
        NC_CHECK(nc_put_var_text  (ncid, vid_cmap_title, title_buf));
        NC_CHECK(nc_put_var_int   (ncid, vid_cmap_resolution, res_buf));
        NC_CHECK(nc_put_var_double(ncid, vid_cmap_energy, data_buf));
        FREE(data_buf); FREE(title_buf); FREE(res_buf);
    }

    STOP(tNetCDF, "PRMTOP NetCDF prep+write");
    NC_CHECK(nc_close(ncid));
    VP0("Successfully saved NetCDF PRMTOP file \"%s\"\n", fname);
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
    int vid_cell_spatial, vid_cell_angular, vid_cell_length, vid_cell_angle, vid_time;
    int dimensionID[2];

    int iAtomCount = iVarArrayElementCount(uUnit->vaAtoms);
    VP0("There are %i atoms\n", iAtomCount);

    NC_CHECK(nc_create(filename, NC_64BIT_OFFSET, &ncid));
    NC_CHECK(nc_def_dim(ncid, "spatial", 3, &did_spatial));
    NC_CHECK(nc_def_var(ncid, "spatial", NC_CHAR, 1, &did_spatial, &vid_spatial));
    NC_CHECK(nc_def_dim(ncid, "atom", iAtomCount, &did_atom));
    NC_CHECK(nc_def_var(ncid, "time", NC_DOUBLE, 0, dimensionID, &vid_time));
    NC_CHECK(nc_put_att_text(ncid, vid_time, "units", 10, "picosecond"));

    dimensionID[0] = did_atom;
    dimensionID[1] = did_spatial;
    NC_CHECK(nc_def_var(ncid, "coordinates", NC_DOUBLE, 2, dimensionID, &vid_coord));
    NC_CHECK(nc_put_att_text(ncid, vid_coord, "units", 8, "angstrom"));

    if (bUnitUseBox(uUnit)) {
        VP0("Using the unit box\n");
        NC_CHECK(nc_def_dim(ncid, "cell_spatial", 3, &did_cell_spatial));
        NC_CHECK(nc_def_var(ncid, "cell_spatial", NC_CHAR, 1, &did_cell_spatial, &vid_cell_spatial));

        NC_CHECK(nc_def_dim(ncid, "label", 5, &did_label));
        NC_CHECK(nc_def_dim(ncid, "cell_angular", 3, &did_cell_angular));
        dimensionID[0] = did_cell_angular;
        dimensionID[1] = did_label;
        NC_CHECK(nc_def_var(ncid, "cell_angular", NC_CHAR, 2, dimensionID, &vid_cell_angular));

        NC_CHECK(nc_def_var(ncid, "cell_lengths", NC_DOUBLE, 1, &did_cell_spatial, &vid_cell_length));
        NC_CHECK(nc_put_att_text(ncid, vid_cell_length, "units", 8, "angstrom"));

        NC_CHECK(nc_def_var(ncid, "cell_angle", NC_DOUBLE, 1, &did_cell_angular, &vid_cell_angle));
        NC_CHECK(nc_put_att_text(ncid, vid_cell_angle, "units", 6, "degree"));
    }

    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "title",
                        strlen(sContainerName(uUnit)), sContainerName(uUnit)));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "application", 5, "AMBER"));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "program", 4, "leap"));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "programVersion", 3, LEAP_VERSION));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 12, "AMBERRESTART"));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "ConventionVersion", 3, "1.0"));

    NC_CHECK(nc_set_fill(ncid, NC_NOFILL, dimensionID));
    NC_CHECK(nc_enddef(ncid));

    NC_CHECK(nc_put_var_text(ncid, vid_spatial, "xyz"));

    if (bUnitUseBox(uUnit)) {
        NC_CHECK(nc_put_var_text(ncid, vid_cell_spatial, "abc"));
        NC_CHECK(nc_put_var_text(ncid, vid_cell_angular,
                             "alpha" "beta " "gamma"));
    }

    double *data;
    data = (double *)MALLOC(3 * iAtomCount * sizeof(double));
    int counter = 0;
    VECTOR vPos;

    double time = 0.0;
    NC_CHECK(nc_put_var_double(ncid, vid_time, &time));

    double dX, dY, dZ, dX2, dY2, dZ2;
    UnitGetBox(uUnit, &dX, &dY, &dZ);
    if (bUnitUseBox(uUnit) && GDefaults.bNoCenter == 0) {
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

    NC_CHECK(nc_put_var_double(ncid, vid_coord, data));

    if (bUnitUseBox(uUnit)) {
        double lengths[3] = { dX, dY, dZ };
        NC_CHECK(nc_put_var_double(ncid, vid_cell_length, lengths));
        lengths[0] = lengths[1] = lengths[2] = dUnitBeta(uUnit) * RADTODEG;
        NC_CHECK(nc_put_var_double(ncid, vid_cell_angle, lengths));
    }

    NC_CHECK(nc_close(ncid));
    FREE(data);
    VP0("Successfully saved NetCDF inpcrd file \"%s\"\n", filename);
}

#else

void UnitIOSaveAmberCoordNetcdf(UNIT uUnit, char *filename)
{
    VPFATALEXIT("Built without NetCDF support. Rebuild with -DBINTRAJ\n");
}

int UnitIOSaveAmberParmNetcdf(const char *fname, UNITt *uUnit,
                              bool bPert, bool bPolar,
                              VARARRAY vaExcludedAtoms,
                              VARARRAY vaExcludedCount,
                              VARARRAY vaNBIndexMatrix,
                              VARARRAY vaNBParameters,
                              VARARRAY vaNBIndex)
{
    VPFATALEXIT("Built without NetCDF support. Rebuild with -DBINTRAJ\n");
    return 0;
}

#endif
