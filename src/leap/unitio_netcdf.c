#include <stdbool.h>

#ifdef BINTRAJ

#include        "basics.h"
#include        "unit.h"
#include        "defaults.h"
#include        "classes.h"
#include        "residue.h"
#include        "unitio.h"
#include        "netcdf.h"
#include        "avl.h"
#include        "sort.h"
#include        "cmap.h"
#include        <assert.h>

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
   ================================================================ */

#include <netcdf.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* ----------------------------------------------------------------
   Original macros from unitio.h
   RESTRAINTLOOP calls FortranWriteDouble so cannot be used in the
   NetCDF path; use NC_RESTRAINTLOOP instead.
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

void
UnitIOSaveAmberParmNetcdf(const char *fname, UNITt *uUnit,
                              bool bPert, bool bPolar,
                              VARARRAY vaExcludedAtoms,
                              VARARRAY vaExcludedCount,
                              VARARRAY vaNBIndexMatrix,
                              VARARRAY vaNBParameters,
                              VARARRAY vaNBIndex)
{
    int ncid;
    int dimid_atoms, dimid_res, dimid_name4, dimid_chainid_len;
    int dimid_bp;   /* bond_type - shared by bond parm and augmentation */
    int dimid_ap;   /* angle_type - shared by angle parm and CHARMM UB */
    int dimid_atom2, dimid_atom3, dimid_atom4, dimid_atom5;
    int n_atoms    = iVarArrayElementCount(uUnit->vaAtoms);
    int n_residues = iVarArrayElementCount(uUnit->vaResidues);
    char s1[8], s2[8], s3[8], s4[8];
    STRING sDesc, sFFType;
    double dKb, dR0, dKpull, dRpull0, dKpress, dRpress0;
    double dKt, dT0;
    double dKp, dP0, dScEE, dScNB;
    double dC, dD;
    int vid_name, vid_charge, vid_atomic_num, vid_itype;

    if (bPert) VPFATALEXIT("Pert atoms not supported\n");

    int iParmSets = iVarArrayElementCount(uUnit->vaParmSets);
    int title_length = 0, fname_length = 0;
    {
        int iPARM = 0;
        PARMSET psLib;
        for (int i = 0; i < iParmSets; i++) {
            psLib = *PVAI(uUnit->vaParmSets, PARMSET, i);
            int len = (int)strlen(psLib->sFname);
            if (len > fname_length) fname_length = len;
            len = (int)strlen(psLib->sTitle);
            if (len > title_length) title_length = len;
            if (!strncmp(psLib->sTitle, "PARM", 4)) iPARM = i;
        }
        psLib = *PVAI(uUnit->vaParmSets, PARMSET, iPARM);
        sprintf(sFFType, "AMBER   %.1024s", psLib->sTitle);
    }
    int chainid_len=1;
    for (int i = 0; i < n_residues; i++) {
        const char *s = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sChainId;
        size_t l = strlen(s);
        if (l > chainid_len) chainid_len = l;
    }

    int CMAP_TYPES = iParmSetTotalCMAPParms(uUnit->psParameters);
    int CMAP_KINDS = iVarArrayElementCount(uUnit->vaCMAPs);

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
    NC_ATT(ncid, NC_GLOBAL, "title",              sContainerName(uUnit));
    NC_ATT(ncid, NC_GLOBAL, "application",        "AMBER");
    NC_ATT(ncid, NC_GLOBAL, "applicationVersion", "1.9");
    NC_ATT(ncid, NC_GLOBAL, "program",            "leap");
    NC_ATT(ncid, NC_GLOBAL, "programVersion",     LEAP_VERSION);
    NC_ATT(ncid, NC_GLOBAL, "Conventions",        "AMBERPRMTOP");
    NC_ATT(ncid, NC_GLOBAL, "ConventionVersion",  "1.0");
    NC_ATT(ncid, NC_GLOBAL, "force_field_type",     sFFType);

    /* Shared dimensions */
    NC_CHECK(nc_def_dim(ncid, "atom",         n_atoms,    &dimid_atoms));
    NC_CHECK(nc_def_dim(ncid, "residues",     n_residues, &dimid_res));
    NC_CHECK(nc_def_dim(ncid, "label_length", LABLEN,     &dimid_name4));
    NC_CHECK(nc_def_dim(ncid, "chainid_length", chainid_len, &dimid_chainid_len));
    NC_CHECK(nc_def_dim(ncid, "atom_pair",    2,          &dimid_atom2));
    NC_CHECK(nc_def_dim(ncid, "atom_triple",  3,          &dimid_atom3));
    NC_CHECK(nc_def_dim(ncid, "atom_quad",    4,          &dimid_atom4));
    if (CMAP_TYPES && CMAP_KINDS && GDefaults.bCMAP)
        NC_CHECK(nc_def_dim(ncid, "atom_quint",5,         &dimid_atom5));
    if (iVarArrayElementCount(uUnit->vaRestraints))
        VP0(" Restraints:  Bond %d  Angle %d  Torsion %d\n",
            iUnitRestraintTypeCount(uUnit, RESTRAINTBOND),
            iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE),
            iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION));
    else
        VP0(" (no restraints)\n");

    /* parameter set file/title metadata */
    int dimid_fname_len, dimid_title_len, dimid_nparm;
    int vid_parm_file, vid_parm_title;
    NC_CHECK(nc_def_dim(ncid, "parmset_title_length", title_length, &dimid_title_len));
    NC_CHECK(nc_def_dim(ncid, "parmset_fname_length", fname_length, &dimid_fname_len));
    NC_CHECK(nc_def_dim(ncid, "parameter_sets",       iParmSets,   &dimid_nparm));
    {
        int dims2[2] = {dimid_nparm, dimid_fname_len};
        NC_CHECK(nc_def_var(ncid, "parm_file",  NC_CHAR, 2, dims2, &vid_parm_file));
        dims2[1] = dimid_title_len;
        NC_CHECK(nc_def_var(ncid, "parm_title", NC_CHAR, 2, dims2, &vid_parm_title));
    }

    /* cell dimensions */
    int vid_a, vid_b, vid_c, vid_alpha, vid_beta, vid_gamma;
    if (bUnitUseBox(uUnit)) {
        NC_CHECK(nc_def_var(ncid, "cell_length_a",    NC_DOUBLE, 0, NULL, &vid_a));
        NC_ATT(ncid, vid_a,     "units", "angstrom");
        NC_CHECK(nc_def_var(ncid, "cell_length_b",    NC_DOUBLE, 0, NULL, &vid_b));
        NC_ATT(ncid, vid_b,     "units", "angstrom");
        NC_CHECK(nc_def_var(ncid, "cell_length_c",    NC_DOUBLE, 0, NULL, &vid_c));
        NC_ATT(ncid, vid_c,     "units", "angstrom");
        NC_CHECK(nc_def_var(ncid, "cell_angle_alpha", NC_DOUBLE, 0, NULL, &vid_alpha));
        NC_ATT(ncid, vid_alpha, "units", "degree");
        NC_CHECK(nc_def_var(ncid, "cell_angle_beta",  NC_DOUBLE, 0, NULL, &vid_beta));
        NC_ATT(ncid, vid_beta,  "units", "degree");
        NC_CHECK(nc_def_var(ncid, "cell_angle_gamma", NC_DOUBLE, 0, NULL, &vid_gamma));
        NC_ATT(ncid, vid_gamma, "units", "degree");
    }

    /* -3- atom names */
    {
        int dims2[2] = {dimid_atoms, dimid_name4};
        NC_CHECK(nc_def_var(ncid, "atom_name", NC_CHAR, 2, dims2, &vid_name));
        NC_ATT(ncid, vid_name, "long_name", "atom name");
        NC_ATT(ncid, vid_name, "pdb_field",  "name");
        NC_ATT(ncid, vid_name, "mmcif_field",
               GDefaults.bCIFReadAuth ? "atom_site.auth_atom_id"
                                      : "atom_site.label_atom_id");
    }

    /* -4- charges */
    NC_CHECK(nc_def_var(ncid, "atom_charge", NC_DOUBLE, 1, &dimid_atoms, &vid_charge));
    NC_ATT(ncid, vid_charge, "units",     "elementary charge");
    NC_ATT(ncid, vid_charge, "long_name", "multiply by 18.2223 for AMBER internal charge units");

    /* -4b- atomic numbers */
    NC_CHECK(nc_def_var(ncid, "atom_atomic_number", NC_INT, 1, &dimid_atoms, &vid_atomic_num));
    NC_ATT(ncid, vid_atomic_num, "long_name", "atomic number");

    /* -5- mass, GB radii, GB screen */
    int vid_mass, vid_radii, vid_screen;
    NC_CHECK(nc_def_var(ncid, "atom_mass",      NC_DOUBLE, 1, &dimid_atoms, &vid_mass));
    NC_ATT(ncid, vid_mass,   "units",     "g/mol");
    NC_CHECK(nc_def_var(ncid, "atom_GB_radius", NC_DOUBLE, 1, &dimid_atoms, &vid_radii));
    NC_ATT(ncid, vid_radii,  "units",     "angstrom");
    NC_ATT(ncid, vid_radii,  "long_name", "GB radius");
    NC_CHECK(nc_def_var(ncid, "atom_GB_screen", NC_DOUBLE, 1, &dimid_atoms, &vid_screen));
    NC_ATT(ncid, vid_screen, "long_name", "GB screening factor");

    /* -6- atom type index */
    NC_CHECK(nc_def_var(ncid, "atom_nonbond_type_index", NC_INT, 1, &dimid_atoms, &vid_itype));
    NC_ATT(ncid, vid_itype, "long_name",          "1-based type index");
    NC_ATT(ncid, vid_itype, "reference_dimension", "nonbond_pair_type");
    NC_ATT(ncid, vid_itype, "reference",           "nonbond_pair_index");

    /* -7- excluded atom count */
    int vid_excl;
    NC_CHECK(nc_def_var(ncid, "number_excluded_atoms", NC_INT, 1, &dimid_atoms, &vid_excl));
    NC_ATT(ncid, vid_excl, "long_name", "excluded atom count per atom");
    NC_ATT(ncid, vid_excl, "encoding",  "cumulative_sum_offset");
    NC_ATT(ncid, vid_excl, "reference", "excluded_atoms_list");

    /* -8- nonbonded parm index matrix + LJ coefficients */
    int vid_nbidx, vid_acoef, vid_bcoef, vid_14acoef, vid_14bcoef;
    int nttyp  = iVarArrayElementCount(vaNBParameters);
    int ntypes = (int)((sqrt(8.0*nttyp + 1.0) - 1.0) / 2.0);
    {
        int dimid_nb, dimid_nttyp;
        assert(ntypes * ntypes == iVarArrayElementCount(vaNBIndexMatrix));

        NC_CHECK(nc_def_dim(ncid, "nonbond_atom_type", ntypes, &dimid_nb));
        int dims2[2] = {dimid_nb, dimid_nb};
        NC_CHECK(nc_def_var(ncid, "nonbond_pair_index", NC_INT, 2, dims2, &vid_nbidx));
        NC_ATT(ncid, vid_nbidx, "long_name",
               "1-based nonbond parm echelon index: max(i,j)*(max(i,j)+1)/2+min(i,j)+1");
        NC_ATT(ncid, vid_nbidx, "reference_dimension", "nonbond_pair_type");

        NC_CHECK(nc_def_dim(ncid, "nonbond_pair_type", nttyp, &dimid_nttyp));
        NC_CHECK(nc_def_var(ncid, "nonbond_LJ_Acoef", NC_DOUBLE, 1, &dimid_nttyp, &vid_acoef));
        NC_ATT(ncid, vid_acoef, "units",     "kcal/mol*angstrom^12");
        NC_ATT(ncid, vid_acoef, "long_name", "LJ parm A coefficients, size NTYPES*(NTYPES+1)/2");
        NC_CHECK(nc_def_var(ncid, "nonbond_LJ_Bcoef", NC_DOUBLE, 1, &dimid_nttyp, &vid_bcoef));
        NC_ATT(ncid, vid_bcoef, "units",     "kcal/mol*angstrom^6");
        NC_ATT(ncid, vid_bcoef, "long_name", "LJ parm B coefficients, size NTYPES*(NTYPES+1)/2");

        if (GDefaults.bCharmm) {
            NC_CHECK(nc_def_var(ncid, "nonbond_LJ14_Acoef", NC_DOUBLE, 1, &dimid_nttyp, &vid_14acoef));
            NC_ATT(ncid, vid_14acoef, "units", "kcal/mol*angstrom^12");
            NC_CHECK(nc_def_var(ncid, "nonbond_LJ14_Bcoef", NC_DOUBLE, 1, &dimid_nttyp, &vid_14bcoef));
            NC_ATT(ncid, vid_14bcoef, "units", "kcal/mol*angstrom^6");
        }
    }

    /* -9- residue labels and type */
    int vid_resname, vid_res_type;
    {
        int dims2[2] = {dimid_res, dimid_name4};
        NC_CHECK(nc_def_var(ncid, "residue_label", NC_CHAR, 2, dims2, &vid_resname));
        NC_ATT(ncid, vid_resname, "long_name",  "residue name");
        NC_ATT(ncid, vid_resname, "pdb_field",  "resName");
        NC_ATT(ncid, vid_resname, "mmcif_field",
               GDefaults.bCIFReadAuth ? "atom_site.auth_comp_id"
                                      : "atom_site.label_comp_id");

        NC_CHECK(nc_def_var(ncid, "residue_type", NC_CHAR, 1, &dimid_res, &vid_res_type));
        NC_ATT(ncid, vid_res_type, "long_name",
               "residue type: w=water, p=protein, n=nucleic, s=saccharide,"
               " i=ion, l=ligand, o=other, ?=undefined");
    }

    /* -10- residue pointer */
    int vid_res_first_atom;
    NC_CHECK(nc_def_var(ncid, "residue_first_atom", NC_INT, 1, &dimid_res, &vid_res_first_atom));
    NC_ATT(ncid, vid_res_first_atom, "long_name",          "1-based first atom index per residue");
    NC_ATT(ncid, vid_res_first_atom, "reference_dimension", "atom");

    /* -11,-12A- bond force constants and equil (+ restraints) */
    int vid_bond_force, vid_bond_eq;
    {
        int nb = iParmSetTotalBondParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
        NC_CHECK(nc_def_dim(ncid, "bond_type", nb + nr, &dimid_bp));
        NC_CHECK(nc_def_var(ncid, "bond_force_constant", NC_DOUBLE, 1, &dimid_bp, &vid_bond_force));
        NC_ATT(ncid, vid_bond_force, "units", "kcal/mol/angstrom^2");
        NC_CHECK(nc_def_var(ncid, "bond_equil_value", NC_DOUBLE, 1, &dimid_bp, &vid_bond_eq));
        NC_ATT(ncid, vid_bond_eq, "units", "angstrom");
    }

    /* -12B..E- bond augmentation, same bond_type dimension */
    int vid_bond_stiffness_pull_adj, vid_bond_equil_pull_adj;
    int vid_bond_stiffness_press_adj, vid_bond_equil_press_adj;
    if (bBondAugmentationFound(uUnit)) {
        NC_CHECK(nc_def_var(ncid, "bond_stiffness_pull_adj",  NC_DOUBLE, 1, &dimid_bp, &vid_bond_stiffness_pull_adj));
        NC_ATT(ncid, vid_bond_stiffness_pull_adj,  "units", "kcal/mol/angstrom^2");
        NC_CHECK(nc_def_var(ncid, "bond_equil_pull_adj",      NC_DOUBLE, 1, &dimid_bp, &vid_bond_equil_pull_adj));
        NC_ATT(ncid, vid_bond_equil_pull_adj,      "units", "angstrom");
        NC_CHECK(nc_def_var(ncid, "bond_stiffness_press_adj", NC_DOUBLE, 1, &dimid_bp, &vid_bond_stiffness_press_adj));
        NC_ATT(ncid, vid_bond_stiffness_press_adj, "units", "kcal/mol/angstrom^2");
        NC_CHECK(nc_def_var(ncid, "bond_equil_press_adj",     NC_DOUBLE, 1, &dimid_bp, &vid_bond_equil_press_adj));
        NC_ATT(ncid, vid_bond_equil_press_adj,     "units", "angstrom");
    }

    /* -13,-14- angle force constants and equil (+ restraints) */
    int vid_angle_force, vid_angle_eq;
    {
        int na = iParmSetTotalAngleParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
        NC_CHECK(nc_def_dim(ncid, "angle_type", na + nr, &dimid_ap));
        NC_CHECK(nc_def_var(ncid, "angle_force_constant", NC_DOUBLE, 1, &dimid_ap, &vid_angle_force));
        NC_ATT(ncid, vid_angle_force, "units", "kcal/mol/radian^2");
        NC_CHECK(nc_def_var(ncid, "angle_equil_value", NC_DOUBLE, 1, &dimid_ap, &vid_angle_eq));
        NC_ATT(ncid, vid_angle_eq, "units", "degree");
    }

    /* CHARMM UB terms */
    VARARRAY ub_force = NULL, ub_equil = NULL;
    VARARRAY ub_atoms = NULL, ub_index = NULL;
    int dimid_UB, dimid_UB_types;
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

        if ( iVarArrayElementCount(ub_index) ) {
        NC_CHECK(nc_def_dim(ncid, "angle_UB13",      iVarArrayElementCount(ub_index), &dimid_UB));
        NC_CHECK(nc_def_dim(ncid, "angle_UB13_type", iVarArrayElementCount(ub_force), &dimid_UB_types));

        int dims2[2] = {dimid_UB, dimid_atom2};
        NC_CHECK(nc_def_var(ncid, "angle_UB13_atom",       NC_INT,    2, dims2,       &vid_UB13_atoms));
        NC_ATT(ncid, vid_UB13_atoms, "long_name",          "1-based atom indices (atoms 1 and 3)");
        NC_ATT(ncid, vid_UB13_atoms, "reference_dimension", "atom");
        NC_CHECK(nc_def_var(ncid, "angle_UB13_type_index", NC_INT,    1, &dimid_UB,   &vid_UB13_index));
        NC_ATT(ncid, vid_UB13_index, "long_name",  "1-based parm index");
        NC_ATT(ncid, vid_UB13_index, "reference",  "angle_UB13_force_constant angle_UB13_equil_value");
        NC_CHECK(nc_def_var(ncid, "angle_UB13_force_constant", NC_DOUBLE, 1, &dimid_UB_types, &vid_UB13_force));
        NC_ATT(ncid, vid_UB13_force, "units", "kcal/mol/angstrom^2");
        NC_CHECK(nc_def_var(ncid, "angle_UB13_equil_value",    NC_DOUBLE, 1, &dimid_UB_types, &vid_UB13_equil));
        NC_ATT(ncid, vid_UB13_equil, "units", "angstrom");
        }
    }

    /* -15..17C- torsion/improper params (+ restraints) */
    // proper+improper+restraint are all collected into one parameter table
    // but now we also subdivde the collection by harmonic form defined by iN==0
    // so it gets a bit convoluted
    int vid_dihe_period, vid_dihe_phase, vid_dihe_force, vid_dihe_scee, vid_dihe_scnb;
    int vid_harm_dihe_force, vid_harm_dihe_equil, n_dihe_parm = 0, n_harm_parm = 0;
    int *harm_map, *dihe_map;
    {
        int dimid_tp, dimid_tp_harm, iN;
        int nt = iParmSetTotalTorsionParms(uUnit->psParameters);
        int ni = iParmSetTotalImproperParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION);
        harm_map = CALLOC(sizeof(int)*(nt+ni+nr));
        dihe_map = CALLOC(sizeof(int)*(nt+ni+nr));
        int index = 0;
        for (int i = 0; i < nt; i++) {
            ParmSetTorsion(uUnit->psParameters, i, NULL, NULL, NULL, NULL,
                           &iN, NULL, NULL, NULL, NULL, NULL);
            if (!iN) harm_map[index] = ++n_harm_parm;
            else dihe_map[index] = ++n_dihe_parm;
            index++;
        }
        for (int i = 0; i < ni; i++) {
            ParmSetImproper(uUnit->psParameters, i, NULL, NULL, NULL, NULL,
                            &iN, NULL, NULL, NULL);
            if (!iN) harm_map[index] = ++n_harm_parm;
            else dihe_map[index] = ++n_dihe_parm;
            index++;
        }
        int nrt = iVarArrayElementCount((uUnit)->vaRestraints);
        for (int i=0; i< nrt; i++) {
            SAVERESTRAINTt *sr = PVAI((uUnit)->vaRestraints, SAVERESTRAINTt, i);
            if (sr->iType == RESTRAINTTORSION && !sr->iN) harm_map[index] = ++n_harm_parm;
            else dihe_map[index] = ++n_dihe_parm;
            index++;
        }

        NC_CHECK(nc_def_dim(ncid, "periodic_dihedral_type", n_dihe_parm, &dimid_tp));
        NC_CHECK(nc_def_var(ncid, "periodic_dihedral_force_constant", NC_DOUBLE, 1, &dimid_tp, &vid_dihe_force));
        NC_ATT(ncid, vid_dihe_force,  "units", "kcal/mol");
        NC_CHECK(nc_def_var(ncid, "periodic_dihedral_periodicity",    NC_INT,    1, &dimid_tp, &vid_dihe_period));
        NC_ATT(ncid, vid_dihe_period, "long_name", "torsion Fourier periodicity n: cos(n*phi-phase)");
        NC_CHECK(nc_def_var(ncid, "periodic_dihedral_phase",          NC_DOUBLE, 1, &dimid_tp, &vid_dihe_phase));
        NC_ATT(ncid, vid_dihe_phase,  "units", "degree");
        NC_CHECK(nc_def_var(ncid, "periodic_dihedral_scee_scale",     NC_DOUBLE, 1, &dimid_tp, &vid_dihe_scee));
        NC_ATT(ncid, vid_dihe_scee,   "long_name", "1-4 electrostatic energy scale factor");
        NC_CHECK(nc_def_var(ncid, "periodic_dihedral_scnb_scale",     NC_DOUBLE, 1, &dimid_tp, &vid_dihe_scnb));
        NC_ATT(ncid, vid_dihe_scnb,   "long_name", "1-4 VdW energy scale factor");

        if (n_harm_parm) {
        NC_CHECK(nc_def_dim(ncid, "harmonic_dihedral_type", n_harm_parm, &dimid_tp_harm));
        NC_CHECK(nc_def_var(ncid, "harmonic_dihedral_force_constant", NC_DOUBLE, 1, &dimid_tp_harm, &vid_harm_dihe_force));
        NC_ATT(ncid, vid_harm_dihe_force,  "units", "kcal/mol/rad^2");
        NC_CHECK(nc_def_var(ncid, "harmonic_dihedral_equil_value",    NC_DOUBLE, 1, &dimid_tp_harm, &vid_harm_dihe_equil));
        NC_ATT(ncid, vid_harm_dihe_equil,  "units", "degree");
        }
    }

    /* -21,-22- bonds inc/excl hydrogen
       Three dimensions:
         bonds_inc_h          — H bonds only (no restraints)
         bonds_excl_h         — non-H bonds only (no restraints)
         bonds_excl_h_r       — non-H bonds + restraints (write arrays use this) */
    int vid_bh_atoms, vid_bh_parm, vid_bnh_atoms, vid_bnh_parm;
    int (*bond_ah)[2], *bond_pih, (*bond_anh)[2], *bond_pinh;
    {
        int nb = iVarArrayElementCount(uUnit->vaBonds);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
        bond_ah   = MALLOC(2 * nb * sizeof(int));
        bond_pih  = MALLOC(nb * sizeof(int));
        bond_anh  = MALLOC(2 * (nb + nr) * sizeof(int));
        bond_pinh = MALLOC((nb + nr) * sizeof(int));
        int nh = 0, nnh = 0;
        for (int i = 0; i < nb; i++) {
            SAVEBONDt *sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, i);
            ATOMt *aA = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1-1)->aAtom;
            ATOMt *aD = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2-1)->aAtom;
            if (iAtomElement(aA) == HYDROGEN || iAtomElement(aD) == HYDROGEN) {
                bond_ah[nh][0] = sbPBond->iAtom1;
                bond_ah[nh][1] = sbPBond->iAtom2;
                bond_pih[nh]   = sbPBond->iParmIndex; nh++;
            } else {
                bond_anh[nnh][0] = sbPBond->iAtom1;
                bond_anh[nnh][1] = sbPBond->iAtom2;
                bond_pinh[nnh]   = sbPBond->iParmIndex; nnh++;
            }
        }
        if (iVarArrayElementCount(uUnit->vaRestraints)) {
            SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
            for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++)
                if (sr->iType == RESTRAINTBOND) {
                    bond_anh[nnh][0] = sr->iAtom1;
                    bond_anh[nnh][1] = sr->iAtom2;
                    bond_pinh[nnh]   = sr->iParmIndex; nnh++;
                }
        }
        int dimid_bh, dimid_bnh, dimid_bnh_r;
        NC_CHECK(nc_def_dim(ncid, "bond_inc_h",      nh,       &dimid_bh));
        NC_CHECK(nc_def_dim(ncid, "bond_excl_h",     nb - nh,  &dimid_bnh));
        NC_CHECK(nc_def_dim(ncid, "bonds_w_restraints",   nnh,      &dimid_bnh_r));

        int dims2[2] = {dimid_bh, dimid_atom2};
        NC_CHECK(nc_def_var(ncid, "bond_inc_hydrogen_atom",       NC_INT, 2, dims2,       &vid_bh_atoms));
        NC_ATT(ncid, vid_bh_atoms, "long_name",          "1-based atom indices");
        NC_ATT(ncid, vid_bh_atoms, "reference_dimension", "atom");
        NC_CHECK(nc_def_var(ncid, "bond_inc_hydrogen_type_index", NC_INT, 1, &dimid_bh,   &vid_bh_parm));
        NC_ATT(ncid, vid_bh_parm,  "long_name", "1-based parm index");
        NC_ATT(ncid, vid_bh_parm,  "reference", "bond_force_constant bond_equil_value");

        dims2[0] = dimid_bnh_r;
        NC_CHECK(nc_def_var(ncid, "bond_without_hydrogen_atom",       NC_INT, 2, dims2,       &vid_bnh_atoms));
        NC_ATT(ncid, vid_bnh_atoms, "long_name",          "1-based atom indices");
        NC_ATT(ncid, vid_bnh_atoms, "reference_dimension", "atom");
        NC_CHECK(nc_def_var(ncid, "bond_without_hydrogen_type_index", NC_INT, 1, &dimid_bnh_r, &vid_bnh_parm));
        NC_ATT(ncid, vid_bnh_parm,  "long_name", "1-based parm index");
        NC_ATT(ncid, vid_bnh_parm,  "reference", "bond_force_constant bond_equil_value");
    }

    /* -23,-24- angles inc/excl hydrogen */
    int vid_ah_atoms, vid_ah_parm, vid_anh_atoms, vid_anh_parm;
    int (*angle_ah)[3], *angle_pih, (*angle_anh)[3], *angle_pinh;
    {
        int na = iVarArrayElementCount(uUnit->vaAngles);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
        angle_ah   = MALLOC(3 * na * sizeof(int));
        angle_anh  = MALLOC(3 * (na + nr) * sizeof(int));
        angle_pih  = MALLOC(na * sizeof(int));
        angle_pinh = MALLOC((na + nr) * sizeof(int));
        int nh = 0, nnh = 0;
        for (int i = 0; i < na; i++) {
            SAVEANGLEt *saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
            ATOMt *aA = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1-1)->aAtom;
            ATOMt *aB = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2-1)->aAtom;
            ATOMt *aC = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3-1)->aAtom;
            if (iAtomElement(aA)==HYDROGEN || iAtomElement(aB)==HYDROGEN ||
                iAtomElement(aC)==HYDROGEN) {
                angle_ah[nh][0] = saPAngle->iAtom1;
                angle_ah[nh][1] = saPAngle->iAtom2;
                angle_ah[nh][2] = saPAngle->iAtom3;
                angle_pih[nh]   = saPAngle->iParmIndex; nh++;
            } else {
                angle_anh[nnh][0] = saPAngle->iAtom1;
                angle_anh[nnh][1] = saPAngle->iAtom2;
                angle_anh[nnh][2] = saPAngle->iAtom3;
                angle_pinh[nnh]   = saPAngle->iParmIndex; nnh++;
            }
        }
        if (iVarArrayElementCount(uUnit->vaRestraints)) {
            SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
            for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++)
                if (sr->iType == RESTRAINTANGLE) {
                    angle_anh[nnh][0] = sr->iAtom1;
                    angle_anh[nnh][1] = sr->iAtom2;
                    angle_anh[nnh][2] = sr->iAtom3;
                    angle_pinh[nnh]   = sr->iParmIndex; nnh++;
                }
        }
        int dimid_ah, dimid_anh, dimid_anh_r;
        NC_CHECK(nc_def_dim(ncid, "angle_inc_h",    nh,       &dimid_ah));
        NC_CHECK(nc_def_dim(ncid, "angle_excl_h",   na - nh,  &dimid_anh));
        NC_CHECK(nc_def_dim(ncid, "angles_w_restraints", nnh,      &dimid_anh_r));

        int dims2[2] = {dimid_ah, dimid_atom3};
        NC_CHECK(nc_def_var(ncid, "angle_inc_hydrogen_atom",       NC_INT, 2, dims2,        &vid_ah_atoms));
        NC_ATT(ncid, vid_ah_atoms, "long_name",          "1-based atom indices");
        NC_ATT(ncid, vid_ah_atoms, "reference_dimension", "atom");
        NC_CHECK(nc_def_var(ncid, "angle_inc_hydrogen_type_index", NC_INT, 1, &dimid_ah,    &vid_ah_parm));
        NC_ATT(ncid, vid_ah_parm,  "long_name", "1-based parm index");
        NC_ATT(ncid, vid_ah_parm,  "reference", "angle_force_constant angle_equil_value");

        dims2[0] = dimid_anh_r;
        NC_CHECK(nc_def_var(ncid, "angle_without_hydrogen_atom",       NC_INT, 2, dims2,        &vid_anh_atoms));
        NC_ATT(ncid, vid_anh_atoms, "long_name",          "1-based atom indices");
        NC_ATT(ncid, vid_anh_atoms, "reference_dimension", "atom");
        NC_CHECK(nc_def_var(ncid, "angle_without_hydrogen_type_index", NC_INT, 1, &dimid_anh_r, &vid_anh_parm));
        NC_ATT(ncid, vid_anh_parm,  "long_name", "1-based parm index");
        NC_ATT(ncid, vid_anh_parm,  "reference", "angle_force_constant angle_equil_value");
    }

    /* -25,-26- torsions */
    int vid_dih_h_atoms, vid_dih_h_parm, vid_dih_nh_atoms, vid_dih_nh_parm;
    int (*dihe_ah)[4], *dihe_pih, (*dihe_anh)[4], *dihe_pinh;
    int vid_harm_dih_atoms=0, vid_harm_dih_parm;
    int (*harm_a)[4], *harm_pi;
    {
        int nT = iVarArrayElementCount(uUnit->vaTorsions);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION);
        dihe_ah   = MALLOC(4 * (nT + nr) * sizeof(int));
        dihe_pih  = MALLOC(nT * sizeof(int));
        dihe_anh  = MALLOC(4 * (nT + nr) * sizeof(int));
        dihe_pinh = MALLOC((nT + nr) * sizeof(int));
        harm_a   = MALLOC(4 * (nT + nr) * sizeof(int));
        harm_pi  = MALLOC(nT * sizeof(int));
        int nh = 0, nnh = 0, nharm = 0;

        for (int i = 0; i < nT; i++) {
            SAVETORSIONt *t = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
            ATOMt *aA = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom1-1)->aAtom;
            ATOMt *aB = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom2-1)->aAtom;
            ATOMt *aC = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom3-1)->aAtom;
            ATOMt *aD = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom4-1)->aAtom;
            int iCalc14 = t->bCalc14 ? 1 : -1;
            int iProper  = t->bProper  ? 1 : -1;
            // XXX CHARM flag was swapping atom 1 <-> 3, which actually does not
            // matter because CHARM impropers are a different array! FIXME:review
            if (harm_map[t->iParmIndex]) {
                harm_a[nharm][0] = t->iAtom1;
                harm_a[nharm][1] = t->iAtom2;
                harm_a[nharm][2] = t->iAtom3 * iCalc14;
                harm_a[nharm][3] = t->iAtom4 * iProper;
                harm_pi[nharm]   = harm_map[t->iParmIndex];
                nharm++;
            } else if (iAtomElement(aA)==HYDROGEN || iAtomElement(aB)==HYDROGEN ||
                iAtomElement(aC)==HYDROGEN || iAtomElement(aD)==HYDROGEN) {
                dihe_ah[nh][0] = t->iAtom1;
                dihe_ah[nh][1] = t->iAtom2;
                dihe_ah[nh][2] = t->iAtom3 * iCalc14;
                dihe_ah[nh][3] = t->iAtom4 * iProper;
                dihe_pih[nh]   = dihe_map[t->iParmIndex];
                nh++;
            } else {
                dihe_anh[nnh][0] = t->iAtom1;
                dihe_anh[nnh][1] = t->iAtom2;
                dihe_anh[nnh][2] = t->iAtom3 * iCalc14;
                dihe_anh[nnh][3] = t->iAtom4 * iProper;
                dihe_pinh[nnh]   = dihe_map[t->iParmIndex];
                nnh++;
            }
        }
        if (nr) {
            SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
            for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++) {
                if (sr->iType != RESTRAINTTORSION) continue;
                if (harm_map[sr->iParmIndex]) {
                    harm_a[nharm][0] = sr->iAtom1;
                    harm_a[nharm][1] = sr->iAtom2;
                    harm_a[nharm][2] = sr->iAtom3;
                    harm_a[nharm][3] = sr->iAtom4;
                    harm_pi[nharm]   = harm_map[sr->iParmIndex];
                    nharm++;
                } else {
                    dihe_anh[nnh][0] = sr->iAtom1;
                    dihe_anh[nnh][1] = sr->iAtom2;
                    dihe_anh[nnh][2] = sr->iAtom3;
                    dihe_anh[nnh][3] = sr->iAtom4;
                    dihe_pinh[nnh]   = dihe_map[sr->iParmIndex];
                    nnh++;
                }
            }
        }
        FREE(harm_map); FREE(dihe_map);

        int dimid_dih_h, dimid_dih_nh, dimid_dih_nh_r, dimid_harm_dih;
        printf("sizes=%d,%d,%d,%d\n",nh,nT-nh,nnh,nharm);
        NC_CHECK(nc_def_dim(ncid, "periodic_dihedrals_inc_h",    nh,        &dimid_dih_h));
        NC_CHECK(nc_def_dim(ncid, "periodic_dihedrals_excl_h",   nT - nh,   &dimid_dih_nh));
        NC_CHECK(nc_def_dim(ncid, "periodic_dihedrals_w_restraints", nnh,       &dimid_dih_nh_r));

        int dims2[2] = {dimid_dih_h, dimid_atom4};
        NC_CHECK(nc_def_var(ncid, "periodic_dihedrals_inc_hydrogen_atom",       NC_INT, 2, dims2,          &vid_dih_h_atoms));
        NC_ATT(ncid, vid_dih_h_atoms,  "long_name", "1-based atom indices");
        NC_ATT(ncid, vid_dih_h_atoms,  "reference_dimension", "atom");
        NC_ATT(ncid, vid_dih_h_atoms,  "note", "atom3,atom4 sign encodes bCalc14,bProper flags");
        NC_CHECK(nc_def_var(ncid, "periodic_dihedrals_inc_hydrogen_type_index", NC_INT, 1, &dimid_dih_h,   &vid_dih_h_parm));
        NC_ATT(ncid, vid_dih_h_parm,   "long_name", "1-based parm index");
        NC_ATT(ncid, vid_dih_h_parm,   "reference", "periodic_dihedral_force_constant dihedral_periodicity dihedral_phase");

        dims2[0] = dimid_dih_nh_r;
        NC_CHECK(nc_def_var(ncid, "periodic_dihedrals_without_hydrogen_atom",       NC_INT, 2, dims2,          &vid_dih_nh_atoms));
        NC_ATT(ncid, vid_dih_nh_atoms, "long_name", "1-based atom indices");
        NC_ATT(ncid, vid_dih_nh_atoms, "reference_dimension", "atom");
        NC_ATT(ncid, vid_dih_nh_atoms, "note", "atom3,atom4 sign encodes bCalc14,bProper flags");
        NC_CHECK(nc_def_var(ncid, "periodic_dihedrals_without_hydrogen_type_index", NC_INT, 1, &dimid_dih_nh_r, &vid_dih_nh_parm));
        NC_ATT(ncid, vid_dih_nh_parm,  "long_name", "1-based parm index");
        NC_ATT(ncid, vid_dih_nh_parm,  "reference", "periodic_dihedral_force_constant dihedral_periodicity dihedral_phase");

        if (nharm) {
        NC_CHECK(nc_def_dim(ncid, "harmonic_dihedral",           nharm,     &dimid_harm_dih));
        dims2[0] = dimid_harm_dih;
        NC_CHECK(nc_def_var(ncid, "harmonic_dihedrals_atom", NC_INT, 2, dims2, &vid_harm_dih_atoms));
        NC_ATT(ncid, vid_harm_dih_atoms, "long_name", "1-based atom indices");
        NC_ATT(ncid, vid_harm_dih_atoms, "reference_dimension", "atom");
        NC_ATT(ncid, vid_harm_dih_atoms, "note", "atom3,atom4 sign encodes bCalc14,bProper flags");
        NC_CHECK(nc_def_var(ncid, "harmonic_dihedrals_type_index", NC_INT, 1, &dimid_harm_dih, &vid_harm_dih_parm));
        NC_ATT(ncid, vid_harm_dih_parm,  "long_name", "1-based parm index");
        NC_ATT(ncid, vid_harm_dih_parm,  "reference", "harmonic_dihedral_force_constant harmonic_dihedral_equil_value");
        }
    }

    /* -27- excluded atoms list */
    int vid_nbex_list;
    {
        int dimid_excl;
        NC_CHECK(nc_def_dim(ncid, "nonbond_excl_list_size",
                            iVarArrayElementCount(vaExcludedAtoms), &dimid_excl));
        NC_CHECK(nc_def_var(ncid, "excluded_atoms_list", NC_INT, 1, &dimid_excl, &vid_nbex_list));
        NC_ATT(ncid, vid_nbex_list, "long_name",          "excluded atom list");
        NC_ATT(ncid, vid_nbex_list, "reference_dimension", "atom");
    }

    /* -28,-29- hbond (skip if none) */
    int vid_hba, vid_hbb;
    {
        int n = iParmSetTotalHBondParms(uUnit->psParameters);
        if (n > 0) {
            int dimid_hb;
            NC_CHECK(nc_def_dim(ncid, "hbond_type", n, &dimid_hb));
            NC_CHECK(nc_def_var(ncid, "hbond_acoef", NC_DOUBLE, 1, &dimid_hb, &vid_hba));
            NC_ATT(ncid, vid_hba, "units", "kcal/mol*angstrom^12");
            NC_CHECK(nc_def_var(ncid, "hbond_bcoef", NC_DOUBLE, 1, &dimid_hb, &vid_hbb));
            NC_ATT(ncid, vid_hbb, "units", "kcal/mol*angstrom^10");
        }
    }

    /* -31- amber atom type strings */
    int vid_atom_type;
    {
        int dims2[2] = {dimid_atoms, dimid_name4};
        NC_CHECK(nc_def_var(ncid, "atom_amber_type", NC_CHAR, 2, dims2, &vid_atom_type));
        NC_ATT(ncid, vid_atom_type, "long_name", "force field atom type name");
    }

    /* -32- tree chain classification */
    int vid_atom_tree;
    {
        int dims2[2] = {dimid_atoms, dimid_name4};
        NC_CHECK(nc_def_var(ncid, "atom_tree_chain_classification", NC_CHAR, 2, dims2, &vid_atom_tree));
        NC_ATT(ncid, vid_atom_tree, "long_name", "mainchain tree class: M,E,S,B,3,4,5,6,X; BLA=unknown");
    }

    /* -35A..C- solvent/box info */
    int vid_mol, vid_fsr, vid_fsm;
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
        NC_CHECK(nc_def_var(ncid, "atoms_per_molecule",   NC_INT, 1, &dimid_mol, &vid_mol));
        NC_ATT(ncid, vid_mol, "long_name", "number of atoms in molecule N");
        NC_CHECK(nc_def_var(ncid, "first_solvent_residue", NC_INT, 0, NULL, &vid_fsr));
        NC_ATT(ncid, vid_fsr, "long_name",          "1-based residue index of first solvent residue");
        NC_ATT(ncid, vid_fsr, "reference_dimension", "residues");
        NC_CHECK(nc_def_var(ncid, "first_solvent_molecule", NC_INT, 0, NULL, &vid_fsm));
        NC_ATT(ncid, vid_fsm, "long_name",          "1-based molecule index of first solvent molecule");
        NC_ATT(ncid, vid_fsm, "reference_dimension", "molecules");
    }

    /* -35D,-35E- cap info */
    int vid_cap_info_atom, vid_cap_info_radius, vid_cap_info_origin;
    if (bUnitUseSolventCap(uUnit)) {
        int dimid_vec;
        NC_CHECK(nc_def_var(ncid, "cap_info_atom",   NC_INT,    0, NULL,      &vid_cap_info_atom));
        NC_ATT(ncid, vid_cap_info_atom,   "long_name",          "Index of the last atom not in the Solvent Cap");
        NC_ATT(ncid, vid_cap_info_atom,   "reference_dimension", "atom");
        NC_CHECK(nc_def_var(ncid, "cap_info_radius", NC_DOUBLE, 0, NULL,      &vid_cap_info_radius));
        NC_ATT(ncid, vid_cap_info_radius, "long_name", "Radius of Solvent Cap");
        NC_ATT(ncid, vid_cap_info_radius, "units",     "angstrom");
        NC_CHECK(nc_def_dim(ncid, "vector", 3, &dimid_vec));
        NC_CHECK(nc_def_var(ncid, "cap_info_origin", NC_DOUBLE, 1, &dimid_vec, &vid_cap_info_origin));
        NC_ATT(ncid, vid_cap_info_origin, "long_name", "Origin of Solvent Cap");
        NC_ATT(ncid, vid_cap_info_origin, "units",     "angstrom");
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

    /* IPOL scalar */
    int vid_IPOL;
    NC_CHECK(nc_def_var(ncid, "polarizability_type", NC_INT, 0, NULL, &vid_IPOL));
    NC_ATT(ncid, vid_IPOL, "long_name", "polarizable force field flag");

    /* polar */
    int vid_polar, vid_damp;
    if (bPolar) {
        NC_CHECK(nc_def_var(ncid, "polarizability",    NC_DOUBLE, 1, &dimid_atoms, &vid_polar));
        NC_ATT(ncid, vid_polar, "units", "angstrom^3");
        NC_CHECK(nc_def_var(ncid, "dipole_damp_factor", NC_DOUBLE, 1, &dimid_atoms, &vid_damp));
    }

    /* PDB chain info */
    int vid_resSeq, vid_res_chain, vid_res_icode;
    bool bHasICodes = false;
    if (GDefaults.bPdbKeepChainId) {
        NC_CHECK(nc_def_var(ncid, "residue_number", NC_INT, 1, &dimid_res, &vid_resSeq));
        NC_ATT(ncid, vid_resSeq, "long_name",  "PDB resSeq");
        NC_ATT(ncid, vid_resSeq, "pdb_field",  "resSeq");
        NC_ATT(ncid, vid_resSeq, "mmcif_field",
               GDefaults.bCIFReadAuth ? "atom_site.auth_seq_id"
                                      : "atom_site.label_seq_id");

        int dims2[2] = {dimid_res, dimid_chainid_len};
        NC_CHECK(nc_def_var(ncid, "residue_chainid", NC_CHAR, 2, dims2, &vid_res_chain));
        NC_ATT(ncid, vid_res_chain, "long_name",  "PDB chainId");
        NC_ATT(ncid, vid_res_chain, "pdb_field",  "chainId");
        NC_ATT(ncid, vid_res_chain, "mmcif_field",
               GDefaults.bCIFReadAuth ? "atom_site.auth_asym_id"
                                      : "atom_site.label_asym_id");

        for (int i = 0; i < n_residues; i++) {
            if (PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sICode[0] != ' ') {
                bHasICodes = true; break;
            }
        }
        if (bHasICodes) {
            NC_CHECK(nc_def_var(ncid, "residue_icode", NC_CHAR, 1, &dimid_res, &vid_res_icode));
            NC_ATT(ncid, vid_res_icode, "long_name",   "PDB iCode");
            NC_ATT(ncid, vid_res_icode, "pdb_field",   "iCode");
            NC_ATT(ncid, vid_res_icode, "mmcif_field", "atom_site.pdbx_PDB_ins_code");
        }
    }

    /* CMAP */
    int vid_cmap_atoms, vid_cmap_type_index;
    int vid_cmap_offset, vid_cmap_resolution, vid_cmap_energy, vid_cmap_title;
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
        NC_CHECK(nc_def_dim(ncid, "phi_psi_cmap",       CMAP_KINDS,       &dimid_kinds));
        NC_CHECK(nc_def_dim(ncid, "phi_psi_cmap_type",  CMAP_TYPES,       &dimid_types));
        NC_CHECK(nc_def_dim(ncid, "phi_psi_cmap_energy", cmap_energy_size, &dimid_values));
        NC_CHECK(nc_def_dim(ncid, "phi_psi_title_len",  cmap_title_len,   &dimid_cmap_title_len));

        int dims2[2] = {dimid_types, dimid_cmap_title_len};
        NC_CHECK(nc_def_var(ncid, "phi_psi_cmap_title",      NC_CHAR,   2, dims2,        &vid_cmap_title));
        NC_CHECK(nc_def_var(ncid, "phi_psi_cmap_offset",     NC_INT,    1, &dimid_types, &vid_cmap_offset));
        NC_ATT(ncid, vid_cmap_offset, "long_name",
               "offset into phi_psi_cmap_energy for this type; "
               "see also phi_psi_cmap_resolution");
        NC_CHECK(nc_def_var(ncid, "phi_psi_cmap_resolution", NC_INT,    1, &dimid_types, &vid_cmap_resolution));
        NC_ATT(ncid, vid_cmap_resolution, "long_name", "grid resolution per map type");

        NC_CHECK(nc_def_var(ncid, "phi_psi_cmap_energy",     NC_DOUBLE, 1, &dimid_values, &vid_cmap_energy));
        NC_ATT(ncid, vid_cmap_energy, "units",     "kcal/mol");
        NC_ATT(ncid, vid_cmap_energy, "long_name", "Concatenated phi/psi torsion energy correction grids");
        NC_ATT(ncid, vid_cmap_energy, "layout",    "row_major: rows=phi, cols=psi");
        NC_ATT(ncid, vid_cmap_energy, "axis_range","[-180,+180) degrees, periodic");
        NC_ATT(ncid, vid_cmap_energy, "navigation","use phi_psi_cmap_offset and phi_psi_cmap_resolution per type index");

        dims2[0] = dimid_kinds; dims2[1] = dimid_atom5;
        NC_CHECK(nc_def_var(ncid, "phi_psi_cmap_atom",       NC_INT, 2, dims2,        &vid_cmap_atoms));
        NC_ATT(ncid, vid_cmap_atoms, "long_name",          "5 raw 1-based atom indices: atoms 1-4=phi, 2-5=psi");
        NC_ATT(ncid, vid_cmap_atoms, "reference_dimension", "atom");
        NC_CHECK(nc_def_var(ncid, "phi_psi_cmap_type_index", NC_INT, 1, &dimid_kinds, &vid_cmap_type_index));
        NC_ATT(ncid, vid_cmap_type_index, "long_name", "1-based index into phi_psi_cmap_type arrays");
        NC_ATT(ncid, vid_cmap_type_index, "reference",
               "phi_psi_cmap_title phi_psi_cmap_offset phi_psi_cmap_resolution");
    }

    NC_CHECK(nc_enddef(ncid));

    /* ====================== write data ====================== */

    /* parameter set names */
    {
        PARMSET psLib;
        char *buf = MALLOC(iParmSets * fname_length);
        memset(buf, ' ', iParmSets * fname_length);
        for (int i = 0; i < iParmSets; i++) {
            psLib = *PVAI(uUnit->vaParmSets, PARMSET, i);
            memcpy(buf + i*fname_length, psLib->sFname, strlen(psLib->sFname));
        }
        NC_CHECK(nc_put_var_text(ncid, vid_parm_file, buf));
        FREE(buf);
        buf = MALLOC(iParmSets * title_length);
        memset(buf, ' ', iParmSets * title_length);
        for (int i = 0; i < iParmSets; i++) {
            psLib = *PVAI(uUnit->vaParmSets, PARMSET, i);
            memcpy(buf + i*title_length, psLib->sTitle, strlen(psLib->sTitle));
        }
        NC_CHECK(nc_put_var_text(ncid, vid_parm_title, buf));
        FREE(buf);
    }

    /* cell dimensions */
    if (bUnitUseBox(uUnit)) {
        double dX, dY, dZ;
        double dAlpha = uUnit->dAlpha * RADTODEG;
        double dBeta  = uUnit->dBeta  * RADTODEG;
        double dGamma = uUnit->dGamma * RADTODEG;
        UnitGetBox(uUnit, &dX, &dY, &dZ);
        NC_CHECK(nc_put_var_double(ncid, vid_a,     &dX));
        NC_CHECK(nc_put_var_double(ncid, vid_b,     &dY));
        NC_CHECK(nc_put_var_double(ncid, vid_c,     &dZ));
        NC_CHECK(nc_put_var_double(ncid, vid_alpha, &dAlpha));
        NC_CHECK(nc_put_var_double(ncid, vid_beta,  &dBeta));
        NC_CHECK(nc_put_var_double(ncid, vid_gamma, &dGamma));
    }

    /* -3- atom names */
    {
        int n = n_atoms;
        char *buf = MALLOC(n * LABLEN);
        memset(buf, ' ', n * LABLEN);
        for (int i = 0; i < n; i++) {
            const char *s = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sName;
            size_t l = strlen(s);
            if (l > LABLEN) s += l - LABLEN;
            memcpy(buf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
        }
        NC_CHECK(nc_put_var_text(ncid, vid_name, buf));
        FREE(buf);
    }

    /* -4- charges */
    {
        double *buf = MALLOC(n_atoms * sizeof(double));
        for (int i = 0; i < n_atoms; i++)
            buf[i] = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->dCharge;
        NC_CHECK(nc_put_var_double(ncid, vid_charge, buf));
        FREE(buf);
    }

    /* -4b- atomic numbers */
    {
        int *buf = MALLOC(n_atoms * sizeof(int));
        for (int i = 0; i < n_atoms; i++)
            buf[i] = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iElement;
        NC_CHECK(nc_put_var_int(ncid, vid_atomic_num, buf));
        FREE(buf);
    }

    /* -5- mass, GB radii, GB screen: merged single ParmSetAtom pass */
    {
        int n = n_atoms;
        double *mass   = MALLOC(n * sizeof(double));
        double *radii  = MALLOC(n * sizeof(double));
        double *screen = MALLOC(n * sizeof(double));
        char sType[MAXTYPELEN];
        double dMass, dPolar, dEpsilon, dRStar, dEpsilon14, dRStar14, dScreenF;
        int iElement, iHybridization;
        SAVEATOMt *sa = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
        for (int i = 0; i < n; i++) {
            int iIdx = iParmSetFindAtom(uUnit->psParameters, sa[i].sType);
            ParmSetAtom(uUnit->psParameters, iIdx, sType, &dMass, &dPolar,
                        &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
                        &iElement, &iHybridization, sDesc);
            mass[i]   = dMass;
            radii[i]  = dGBRadiusForAtom(&sa[i], iElement, dMass, (i+1 == n));
            screen[i] = dGBScreenForElement(iElement);
        }
        NC_CHECK(nc_put_var_double(ncid, vid_mass,   mass));
        NC_CHECK(nc_put_var_double(ncid, vid_radii,  radii));
        NC_CHECK(nc_put_var_double(ncid, vid_screen, screen));
        FREE(mass); FREE(radii); FREE(screen);
    }

    /* -6- atom type index */
    {
        int *buf = MALLOC(n_atoms * sizeof(int));
        for (int i = 0; i < n_atoms; i++) {
            int iAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iTypeIndex - 1;
            buf[i] = *PVAI(vaNBIndex, int, iAtom) + 1;
        }
        NC_CHECK(nc_put_var_int(ncid, vid_itype, buf));
        FREE(buf);
    }

    /* -7- excluded atom counts */
    NC_CHECK(nc_put_var_int(ncid, vid_excl, PVAI(vaExcludedCount, int, 0)));

    /* -8- nonbonded parm index matrix */
    NC_CHECK(nc_put_var_int(ncid, vid_nbidx, PVAI(vaNBIndexMatrix, int, 0)));

    /* -9- residue labels and type */
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
        NC_CHECK(nc_put_var_text(ncid, vid_resname,  buf));
        NC_CHECK(nc_put_var_text(ncid, vid_res_type, buf_type));
        FREE(buf); FREE(buf_type);
    }

    /* -10- residue pointer */
    {
        int *buf = MALLOC(n_residues * sizeof(int));
        for (int i = 0; i < n_residues; i++)
            buf[i] = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iAtomStartIndex;
        NC_CHECK(nc_put_var_int(ncid, vid_res_first_atom, buf));
        FREE(buf);
    }

    /* -11,-12A- bond force constants and equil */
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
        NC_CHECK(nc_put_var_double(ncid, vid_bond_eq,    r0));
        FREE(kb); FREE(r0);
    }

    /* -12B..E- bond augmentation */
    if (bBondAugmentationFound(uUnit)) {
        int nb = iParmSetTotalBondParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
        int n  = nb + nr;
        double *kpull   = MALLOC(n*sizeof(double));
        double *rpull0  = MALLOC(n*sizeof(double));
        double *kpress  = MALLOC(n*sizeof(double));
        double *rpress0 = MALLOC(n*sizeof(double));
        for (int i = 0; i < nb; i++) {
            ParmSetBond(uUnit->psParameters, i, s1, s2, &dKb, &dR0,
                        &dKpull, &dRpull0, &dKpress, &dRpress0, sDesc);
            kpull[i] = dKpull; rpull0[i] = dRpull0;
            kpress[i] = dKpress; rpress0[i] = dRpress0;
        }
        for (int i = nb; i < n; i++) {
            kpull[i] = 0.0;  rpull0[i]  = 100.0;
            kpress[i] = 0.0; rpress0[i] = -100.0;
        }
        NC_CHECK(nc_put_var_double(ncid, vid_bond_stiffness_pull_adj,  kpull));
        NC_CHECK(nc_put_var_double(ncid, vid_bond_equil_pull_adj,      rpull0));
        NC_CHECK(nc_put_var_double(ncid, vid_bond_stiffness_press_adj, kpress));
        NC_CHECK(nc_put_var_double(ncid, vid_bond_equil_press_adj,     rpress0));
        FREE(kpull); FREE(rpull0); FREE(kpress); FREE(rpress0);
    }

    /* -13,-14- angle force constants and equil */
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
        NC_CHECK(nc_put_var_double(ncid, vid_angle_eq,    t0));
        FREE(kt); FREE(t0);
    }

    /* CHARMM UB terms */
    if (GDefaults.bCharmm && ub_atoms) {
        if (iVarArrayElementCount(ub_atoms)) {
            NC_CHECK(nc_put_var_int   (ncid, vid_UB13_atoms, PVAI(ub_atoms, int,    0)));
            NC_CHECK(nc_put_var_int   (ncid, vid_UB13_index, PVAI(ub_index, int,    0)));
            NC_CHECK(nc_put_var_double(ncid, vid_UB13_force, PVAI(ub_force, double, 0)));
            NC_CHECK(nc_put_var_double(ncid, vid_UB13_equil, PVAI(ub_equil, double, 0)));
        }
        VarArrayDestroy(&ub_atoms);
        VarArrayDestroy(&ub_index);
        VarArrayDestroy(&ub_force);
        VarArrayDestroy(&ub_equil);
    }

    /* -15..17C- torsion/improper params */
    {
        int iN;
        int nt = iParmSetTotalTorsionParms(uUnit->psParameters);
        int ni = iParmSetTotalImproperParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION);
        int n  = nt + ni + nr;
        double *kp   = MALLOC(n * sizeof(double));
        int    *per  = MALLOC(n * sizeof(int));
        double *p0   = MALLOC(n * sizeof(double));
        double *scee = MALLOC(n * sizeof(double));
        double *scnb = MALLOC(n * sizeof(double));
        for (int i = 0; i < nt; i++) {
            ParmSetTorsion(uUnit->psParameters, i, s1, s2, s3, s4,
                           &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc);
            kp[i]  = dKp; per[i] = iN; p0[i] = dP0 * RADTODEG;
            scee[i] = (dScEE < 0.0) ? GDefaults.dSceeScaleFactor : dScEE;
            scnb[i] = (dScNB < 0.0) ? GDefaults.dScnbScaleFactor : dScNB;
        }
        for (int i = 0; i < ni; i++) {
            ParmSetImproper(uUnit->psParameters, i, s1, s2, s3, s4,
                            &iN, &dKp, &dP0, sDesc);
            kp[nt+i] = dKp; per[nt+i] = iN; p0[nt+i] = dP0 * RADTODEG;
            scee[nt+i] = 0.0; scnb[nt+i] = 0.0;
        }
        NC_RESTRAINTLOOP(uUnit, RESTRAINTTORSION, dKx, kp,  nt+ni);
        NC_RESTRAINTLOOP(uUnit, RESTRAINTTORSION, dX0, p0,  nt+ni);
        NC_RESTRAINTLOOP(uUnit, RESTRAINTTORSION, iN,  per, nt+ni);
        NC_CHECK(nc_put_var_double(ncid, vid_dihe_force,  kp));
        NC_CHECK(nc_put_var_int   (ncid, vid_dihe_period, per));
        NC_CHECK(nc_put_var_double(ncid, vid_dihe_phase,  p0));
        NC_CHECK(nc_put_var_double(ncid, vid_dihe_scee,   scee));
        NC_CHECK(nc_put_var_double(ncid, vid_dihe_scnb,   scnb));
        FREE(kp); FREE(per); FREE(p0); FREE(scee); FREE(scnb);
    }

    /* -19,-20- LJ coefficients */
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

    /* -21,-22- bonds */
    NC_CHECK(nc_put_var_int(ncid, vid_bh_atoms,  (int*)bond_ah));
    NC_CHECK(nc_put_var_int(ncid, vid_bh_parm,   bond_pih));
    FREE(bond_ah); FREE(bond_pih);
    NC_CHECK(nc_put_var_int(ncid, vid_bnh_atoms, (int*)bond_anh));
    NC_CHECK(nc_put_var_int(ncid, vid_bnh_parm,  bond_pinh));
    FREE(bond_anh); FREE(bond_pinh);

    /* -23,-24- angles */
    NC_CHECK(nc_put_var_int(ncid, vid_ah_atoms,  (int*)angle_ah));
    NC_CHECK(nc_put_var_int(ncid, vid_ah_parm,   angle_pih));
    FREE(angle_ah); FREE(angle_pih);
    NC_CHECK(nc_put_var_int(ncid, vid_anh_atoms, (int*)angle_anh));
    NC_CHECK(nc_put_var_int(ncid, vid_anh_parm,  angle_pinh));
    FREE(angle_anh); FREE(angle_pinh);

    /* -25,-26- torsions */
    NC_CHECK(nc_put_var_int(ncid, vid_dih_h_atoms,  (int*)dihe_ah));
    NC_CHECK(nc_put_var_int(ncid, vid_dih_h_parm,   dihe_pih));
    FREE(dihe_ah); FREE(dihe_pih);
    NC_CHECK(nc_put_var_int(ncid, vid_dih_nh_atoms, (int*)dihe_anh));
    NC_CHECK(nc_put_var_int(ncid, vid_dih_nh_parm,  dihe_pinh));
    FREE(dihe_anh); FREE(dihe_pinh);
    if (vid_harm_dih_atoms) {
    NC_CHECK(nc_put_var_int(ncid, vid_harm_dih_atoms, (int*)harm_a));
    NC_CHECK(nc_put_var_int(ncid, vid_harm_dih_parm,  harm_pi));
    }
    FREE(harm_a); FREE(harm_pi);

    /* -27- excluded atoms list */
    NC_CHECK(nc_put_var_int(ncid, vid_nbex_list, PVAI(vaExcludedAtoms, int, 0)));

    /* -28,-29- hbond */
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

    /* -31- amber atom type strings */
    {
        char *buf = MALLOC(n_atoms * LABLEN);
        memset(buf, ' ', n_atoms * LABLEN);
        for (int i = 0; i < n_atoms; i++) {
            const char *s = sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom);
            size_t l = strlen(s);
            if (l > LABLEN) s += l - LABLEN;
            memcpy(buf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
        }
        NC_CHECK(nc_put_var_text(ncid, vid_atom_type, buf));
        FREE(buf);
    }

    /* -32- tree chain classification */
    {
        char *buf = MALLOC(n_atoms * LABLEN);
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

    /* -35A..C- solvent/box */
    if (bUnitUseBox(uUnit)) {
        NC_CHECK(nc_put_var_int(ncid, vid_fsr, &iFirstSolvRes));
        NC_CHECK(nc_put_var_int(ncid, vid_fsm, &iFS1));
        NC_CHECK(nc_put_var_int(ncid, vid_mol, PVAI(uUnit->vaAtomsPerMolecule, int, 0)));
    }

    /* -35D,-35E- cap info */
    if (bUnitUseSolventCap(uUnit)) {
        int ires = uUnit->iCapTempInt;
        double dX, dY, dZ, dR;
        UnitGetSolventCap(uUnit, &dX, &dY, &dZ, &dR);
        double cap[3] = {dX, dY, dZ};
        NC_CHECK(nc_put_var_int   (ncid, vid_cap_info_atom,   &ires));
        NC_CHECK(nc_put_var_double(ncid, vid_cap_info_radius, &dR));
        NC_CHECK(nc_put_var_double(ncid, vid_cap_info_origin,  cap));
    }

    /* IPOL */
    {
        int val = GDefaults.iIPOL;
        NC_CHECK(nc_put_var_int(ncid, vid_IPOL, &val));
    }

    /* polar */
    if (bPolar) {
        int n = n_atoms;
        int iCount = 0;
        double *polar = MALLOC(n * sizeof(double));
        double *damp  = MALLOC(n * sizeof(double));
        char sType[MAXTYPELEN];
        double dMass, dPolar, dEpsilon, dRStar, dEpsilon14, dRStar14, dScreenF;
        int iElement, iHybridization;
        SAVEATOMt *a = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
        for (int i = 0; i < n; i++, a++) {
            int iIdx = iParmSetFindAtom(uUnit->psParameters, a->sType);
            ParmSetAtom(uUnit->psParameters, iIdx, sType, &dMass, &dPolar,
                        &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
                        &iElement, &iHybridization, sDesc);
            if (dPolar == -1.0) { dPolar = 0.0; iCount++; }
            polar[i] = dPolar;
            damp[i]  = (dScreenF == 0.0) ? GDefaults.dDipoleDampFactor : dScreenF;
        }
        if (iCount > 0)
            VP0("Total atoms with default polarization=0.0: %d of %d\n", iCount, n);
        NC_CHECK(nc_put_var_double(ncid, vid_polar, polar));
        NC_CHECK(nc_put_var_double(ncid, vid_damp,  damp));
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
            for (int j=0; chainid[j]; j++) cbuf[i*chainid_len+j] = chainid[j];
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
        NC_CHECK(nc_put_var_int(ncid, vid_cmap_atoms,      atoms_buf));
        NC_CHECK(nc_put_var_int(ncid, vid_cmap_type_index, index_buf));
        FREE(atoms_buf); FREE(index_buf);

        double *data_buf       = MALLOC(cmap_energy_size * sizeof(double));
        char   *title_buf      = MALLOC(cmap_title_len * CMAP_TYPES);
        int    *offset_buf     = MALLOC(sizeof(int) * CMAP_TYPES);
        int    *resolution_buf = MALLOC(sizeof(int) * CMAP_TYPES);
        int offset = 0;
        memset(title_buf, ' ', cmap_title_len * CMAP_TYPES);
        for (int i = 0; i < CMAP_TYPES; i++) {
            CMAPt cmap;
            ParmSetCMAP(uUnit->psParameters, i, &cmap, false);
            memcpy(title_buf + i*cmap_title_len, cmap.title, strlen(cmap.title));
            int m = cmap.resolution * cmap.resolution;
            memcpy(data_buf + offset, cmap.map, m * sizeof(double));
            offset_buf[i]     = offset;
            resolution_buf[i] = cmap.resolution;
            offset += m;
        }
        NC_CHECK(nc_put_var_text  (ncid, vid_cmap_title,      title_buf));
        NC_CHECK(nc_put_var_int   (ncid, vid_cmap_offset,     offset_buf));
        NC_CHECK(nc_put_var_int   (ncid, vid_cmap_resolution, resolution_buf));
        NC_CHECK(nc_put_var_double(ncid, vid_cmap_energy,     data_buf));
        FREE(data_buf); FREE(title_buf); FREE(offset_buf); FREE(resolution_buf);
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
    int vid_cell_spatial, vid_cell_angular, vid_cell_length, vid_cell_angle,
        vid_time;
    int dimensionID[NC_MAX_VAR_DIMS];

    int iAtomCount = iVarArrayElementCount(uUnit->vaAtoms);
    VP0("There are %i atoms\n", iAtomCount);

    NC_CHECK(nc_create(filename, NC_64BIT_OFFSET, &ncid));

    NC_CHECK(nc_def_dim(ncid, "spatial", 3, &did_spatial));

    dimensionID[0] = did_spatial;
    NC_CHECK(nc_def_var(ncid, "spatial", NC_CHAR, 1, dimensionID, &vid_spatial));

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
        dimensionID[0] = did_cell_spatial;
        NC_CHECK(nc_def_var(ncid, "cell_spatial", NC_CHAR, 1, dimensionID, &vid_cell_spatial));

        NC_CHECK(nc_def_dim(ncid, "label", 5, &did_label));
        NC_CHECK(nc_def_dim(ncid, "cell_angular", 3, &did_cell_angular));
        dimensionID[0] = did_cell_angular;
        dimensionID[1] = did_label;
        NC_CHECK(nc_def_var(ncid, "cell_angular", NC_CHAR, 2, dimensionID, &vid_cell_angular));

        dimensionID[0] = did_cell_spatial;
        NC_CHECK(nc_def_var(ncid, "cell_lengths", NC_DOUBLE, 1, dimensionID, &vid_cell_length));
        NC_CHECK(nc_put_att_text(ncid, vid_cell_length, "units", 8, "angstrom"));

        dimensionID[0] = did_cell_angular;
        NC_CHECK(nc_def_var(ncid, "cell_angle", NC_DOUBLE, 1, dimensionID, &vid_cell_angle));
        NC_CHECK(nc_put_att_text(ncid, vid_cell_angle, "units", 6, "degree"));
    }

    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "title",
                        strlen(sContainerName(uUnit)), sContainerName(uUnit)));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "application",      5, "AMBER"));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "program",          4, "leap"));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "programVersion",   3, LEAP_VERSION));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "Conventions",     12, "AMBERRESTART"));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "ConventionVersion", 3, "1.0"));

    NC_CHECK(nc_set_fill(ncid, NC_NOFILL, dimensionID));
    NC_CHECK(nc_enddef(ncid));

    size_t start[2]={0}, count[2]={0};
    count[0] = 3;
    NC_CHECK(nc_put_vara_text(ncid, vid_spatial, start, count, "xyz"));

    if (bUnitUseBox(uUnit)) {
        count[1] = 1;
        NC_CHECK(nc_put_vara_text(ncid, vid_cell_spatial, start, count, "abc"));
        count[1] = 5;
        NC_CHECK(nc_put_vara_text(ncid, vid_cell_angular, start, count,
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

    start[0] = start[1] = 0;
    count[0] = iAtomCount; count[1] = 3;
    NC_CHECK(nc_put_vara_double(ncid, vid_coord, start, count, data));

    if (bUnitUseBox(uUnit)) {
        count[0] = 3; count[1] = 0;
        double lengths[3] = { dX, dY, dZ };
        NC_CHECK(nc_put_vara_double(ncid, vid_cell_length, start, count, lengths));
        lengths[0] = lengths[1] = lengths[2] = dUnitBeta(uUnit) * RADTODEG;
        NC_CHECK(nc_put_vara_double(ncid, vid_cell_angle, start, count, lengths));
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
