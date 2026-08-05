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
    int dimid_atoms, dimid_res, dimid_name4;
    int dimid_bp;   /* NUMBND - shared by bond parm and augmentation sections */
    int dimid_ap;   /* NUMANG - shared by angle parm and CHARMM UB sections */
    int dimid_atom2, dimid_atom3, dimid_atom4;
    int n_atoms    = iVarArrayElementCount(uUnit->vaAtoms);
    int n_residues = iVarArrayElementCount(uUnit->vaResidues);
    char s1[8], s2[8], s3[8], s4[8];
    STRING sDesc;
    double dKb, dR0, dKpull, dRpull0, dKpress, dRpress0;
    double dKt, dT0, dTkub, dRkub;
    double dKp, dP0, dScEE, dScNB;
    double dC, dD;
    int vid_name, vid_charge, vid_atomic_num, vid_itype;

    if (bPert) VPFATALEXIT("Pert atoms not supported\n");

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
    NC_CHECK(nc_def_dim(ncid, "atoms",   n_atoms,    &dimid_atoms));
    NC_CHECK(nc_def_dim(ncid, "residues",n_residues, &dimid_res));
    NC_CHECK(nc_def_dim(ncid, "label_length",LABLEN,     &dimid_name4));
    NC_CHECK(nc_def_dim(ncid, "atom_pair",   2,          &dimid_atom2));
    NC_CHECK(nc_def_dim(ncid, "atom_triple", 3,          &dimid_atom3));
    NC_CHECK(nc_def_dim(ncid, "atom_quad",   4,          &dimid_atom4));

    if (iVarArrayElementCount(uUnit->vaRestraints)) {
        VP0(" Restraints:  Bond %d  Angle %d  Torsion %d\n",
            iUnitRestraintTypeCount(uUnit, RESTRAINTBOND),
            iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE),
            iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION));
    } else VP0(" (no restraints)\n");

    /* cell dimensions */
    int vid_a, vid_b, vid_c, vid_alpha, vid_beta, vid_gamma;
    if (bUnitUseBox(uUnit)) {
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
    }

    /* -3- atom names */
    {
        int dims2[2] = {dimid_atoms, dimid_name4};
        NC_CHECK(nc_def_var(ncid, "atom_name", NC_CHAR, 2, dims2, &vid_name));
        NC_CHECK(nc_put_att_text(ncid, vid_name, "long_name", 9, "atom name"));
        NC_CHECK(nc_put_att_text(ncid, vid_name, "pdb_field",     4, "name"));
        if (GDefaults.bCIFReadAuth)
            NC_CHECK(nc_put_att_text(ncid, vid_name, "mmcif_field",    22, "atom_site.auth_atom_id"));
        else NC_CHECK(nc_put_att_text(ncid, vid_name, "mmcif_field",  23, "atom_site.label_atom_id"));
    }

/* ================================================================
   1. gbrad_for_element()
   ================================================================ */

    /* -4- charges (raw double, no AMBER unit scaling on write) */
        NC_CHECK(nc_def_var(ncid, "atom_charge", NC_DOUBLE, 1, &dimid_atoms, &vid_charge));
        NC_CHECK(nc_put_att_text(ncid, vid_charge, "units", 18, "elementary charge"));
        NC_CHECK(nc_put_att_text(ncid, vid_charge, "long_name", 51,
                                 "multiply by 18.2223 for AMBER internal charge units"));

    /* -4b- atomic numbers */
        NC_CHECK(nc_def_var(ncid, "atom_atomic_number", NC_INT, 1, &dimid_atoms, &vid_atomic_num));
        NC_CHECK(nc_put_att_text(ncid, vid_atomic_num, "long_name", 13, "atomic number"));

    /* -5-, RADII, SCREEN: merged single ParmSetAtom pass */
    int vid_mass, vid_radii, vid_screen;
    NC_CHECK(nc_def_var(ncid, "mass", NC_DOUBLE, 1, &dimid_atoms, &vid_mass));
    NC_CHECK(nc_put_att_text(ncid, vid_mass, "units", 5, "g/mol"));

    NC_CHECK(nc_def_var(ncid, "radii", NC_DOUBLE, 1, &dimid_atoms, &vid_radii));
    NC_CHECK(nc_put_att_text(ncid, vid_radii, "units", 8, "angstrom"));
    NC_CHECK(nc_put_att_text(ncid, vid_radii, "long_name", 9, "GB radius"));

    NC_CHECK(nc_def_var(ncid, "screen", NC_DOUBLE, 1, &dimid_atoms, &vid_screen));
    NC_CHECK(nc_put_att_text(ncid, vid_screen, "long_name", 20, "GB screening factor"));

    /* -6- atom type index */
        NC_CHECK(nc_def_var(ncid, "atom_type_index", NC_INT, 1, &dimid_atoms, &vid_itype));
        NC_CHECK(nc_put_att_text(ncid, vid_itype, "long_name",          18, "1-based type index"));
        NC_CHECK(nc_put_att_text(ncid, vid_itype, "reference_dimension", 7, "NTYPES2"));
        NC_CHECK(nc_put_att_text(ncid, vid_itype, "reference",          20, "nonbonded_parm_index"));

    /* -7- number of excluded atoms per atom */
        int vid_excl;
        NC_CHECK(nc_def_var(ncid, "excluded_atoms", NC_INT, 1, &dimid_atoms, &vid_excl));
        NC_CHECK(nc_put_att_text(ncid, vid_excl, "long_name", 28, "excluded atom count per atom"));
        NC_CHECK(nc_put_att_text(ncid, vid_excl, "encoding",  21, "cumulative_sum_offset"));
        NC_CHECK(nc_put_att_text(ncid, vid_excl, "reference", 19, "excluded_atoms_list"));

    /* -8- nonbonded parm index matrix */
    int vid_nbidx, vid_acoef, vid_bcoef, vid_14acoef, vid_14bcoef;
    int nttyp = iVarArrayElementCount(vaNBParameters);
    int ntypes = (sqrt(8*nttyp+1)-1)/2;
    {
        int dimid_nb;
        int ntypes2 = iVarArrayElementCount(vaNBIndexMatrix);
        int dimid_nttyp;

        assert(ntypes*ntypes == ntypes2);
        NC_CHECK(nc_def_dim(ncid, "nonbond_atom_types", ntypes, &dimid_nb));
        int dims2[2] = {dimid_nb, dimid_nb};
        NC_CHECK(nc_def_var(ncid, "nonbond_pair_index", NC_INT, 2, dims2, &vid_nbidx));
        const char *long_name="1-based nonbond parm echelon index: max(i,j)*(max(i,j)+1)/2+min(i,j)+1";
        NC_CHECK(nc_put_att_text(ncid, vid_nbidx, "long_name", strlen(long_name),long_name));
        NC_CHECK(nc_put_att_text(ncid, vid_nbidx, "reference_dimension", 18, "nonbond_pair_types"));
        // -19-,-20-
        NC_CHECK(nc_def_dim(ncid, "nonbond_pair_types", nttyp, &dimid_nttyp)); // number of pair types
        NC_CHECK(nc_def_var(ncid, "nonbond_LJ_acoef", NC_DOUBLE, 1, &dimid_nttyp, &vid_acoef));
        NC_CHECK(nc_put_att_text(ncid, vid_acoef, "units", 20, "kcal/mol*angstrom^12"));
        NC_CHECK(nc_put_att_text(ncid, vid_acoef, "long_name", 48,
                                 "LJ parm A coefficients, size NTYPES*(NTYPES+1)/2"));
        NC_CHECK(nc_def_var(ncid, "nonbond_LJ_bcoef", NC_DOUBLE, 1, &dimid_nttyp, &vid_bcoef));
        NC_CHECK(nc_put_att_text(ncid, vid_bcoef, "units", 19, "kcal/mol*angstrom^6"));
        NC_CHECK(nc_put_att_text(ncid, vid_acoef, "long_name", 48,
                                 "LJ parm B coefficients, size NTYPES*(NTYPES+1)/2"));
        if (GDefaults.iCharmm) {
            NC_CHECK(nc_def_var(ncid, "nonbond_LJ14_acoef", NC_DOUBLE, 1, &dimid_nttyp, &vid_14acoef));
            NC_CHECK(nc_put_att_text(ncid, vid_14acoef, "units", 20, "kcal/mol*angstrom^12"));
            NC_CHECK(nc_def_var(ncid, "nonbond_LJ14_Bcoef", NC_DOUBLE, 1, &dimid_nttyp, &vid_14bcoef));
            NC_CHECK(nc_put_att_text(ncid, vid_14bcoef, "units", 19, "kcal/mol*angstrom^6"));
        }
    }

    /* -9- residue labels */
    int vid_resname, vid_res_type;
    {
        int dims2[2] = { dimid_res, dimid_name4};
        NC_CHECK(nc_def_var(ncid, "residue_label", NC_CHAR, 2, dims2, &vid_resname));
        NC_CHECK(nc_put_att_text(ncid, vid_resname, "long_name", 12, "residue name"));
        NC_CHECK(nc_put_att_text(ncid, vid_resname, "pdb_field",     7, "resName"));
        if (GDefaults.bCIFReadAuth)
            NC_CHECK(nc_put_att_text(ncid, vid_resname, "mmcif_field",    22, "atom_site.auth_comp_id"));
        else NC_CHECK(nc_put_att_text(ncid, vid_resname, "mmcif_field",  23, "atom_site.label_comp_id"));
        NC_CHECK(nc_def_var(ncid, "residue_type", NC_CHAR, 1, dims2, &vid_res_type));
        const char *type_info = "residue type: w=water, p=protein, n=nucleic, s=saccharide,"
                                " i=ion, l=ligand, o=other, ?=undefined";
        NC_CHECK(nc_put_att_text(ncid, vid_res_type, "long_name", strlen(type_info), type_info));
    }

    /* -10- residue pointer */
        int vid_res_first_atom;
        NC_CHECK(nc_def_var(ncid, "residue_first_atom", NC_INT, 1, &dimid_res, &vid_res_first_atom));
        NC_CHECK(nc_put_att_text(ncid, vid_res_first_atom, "long_name", 36,
                                 "1-based first atom index per residue"));
        NC_CHECK(nc_put_att_text(ncid, vid_res_first_atom, "reference_dimension", 5, "atoms"));

    /* -11,-12A- bond force constants and equil (+ restraints) */
    int vid_bond_force, vid_bond_eq;
    {
        int nb = iParmSetTotalBondParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
        int n  = nb + nr;
        NC_CHECK(nc_def_dim(ncid, "bond_types", n, &dimid_bp));
        NC_CHECK(nc_def_var(ncid, "bond_force_constant", NC_DOUBLE, 1, &dimid_bp, &vid_bond_force));
        NC_CHECK(nc_put_att_text(ncid, vid_bond_force, "units", 19, "kcal/mol/angstrom^2"));
        NC_CHECK(nc_def_var(ncid, "bond_equil_value", NC_DOUBLE, 1, &dimid_bp, &vid_bond_eq));
        NC_CHECK(nc_put_att_text(ncid, vid_bond_eq, "units", 8, "angstrom"));
    }

    /* -12B..E- bond augmentation, same NUMBND dimension as above */
    int vid_bond_stiffness_pull_adj, vid_bond_equil_pull_adj;
    int vid_bond_stiffness_press_adj, vid_bond_equil_press_adj;
    if (BondAugmentationFound(uUnit) == 1) {
        NC_CHECK(nc_def_var(ncid, "bond_stiffness_pull_adj", NC_DOUBLE, 1, &dimid_bp, &vid_bond_stiffness_pull_adj));
        NC_CHECK(nc_put_att_text(ncid, vid_bond_stiffness_pull_adj, "units", 19, "kcal/mol/angstrom^2"));
        NC_CHECK(nc_def_var(ncid, "bond_equil_pull_adj", NC_DOUBLE, 1, &dimid_bp, &vid_bond_equil_pull_adj));
        NC_CHECK(nc_put_att_text(ncid, vid_bond_equil_pull_adj, "units", 8, "angstrom"));
        NC_CHECK(nc_def_var(ncid, "bond_stiffness_press_adj", NC_DOUBLE, 1, &dimid_bp, &vid_bond_stiffness_press_adj));
        NC_CHECK(nc_put_att_text(ncid, vid_bond_stiffness_press_adj, "units", 19, "kcal/mol/angstrom^2"));
        NC_CHECK(nc_def_var(ncid, "bond_equil_press_adj", NC_DOUBLE, 1, &dimid_bp, &vid_bond_equil_press_adj));
        NC_CHECK(nc_put_att_text(ncid, vid_bond_equil_press_adj, "units", 8, "angstrom"));
    }

    /* -13,-14- angle force constants and equil (+ restraints) */
    int vid_angle_force, vid_angle_eq, vid_angle_UB_force, vid_angle_UB_equil;
    {
        int na = iParmSetTotalAngleParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
        int n  = na + nr;
        NC_CHECK(nc_def_dim(ncid, "angle_types", n, &dimid_ap));
        NC_CHECK(nc_def_var(ncid, "angle_force_constant", NC_DOUBLE, 1, &dimid_ap, &vid_angle_force));
        NC_CHECK(nc_put_att_text(ncid, vid_angle_force, "units", 17, "kcal/mol/radian^2"));
        NC_CHECK(nc_def_var(ncid, "angle_equil_value", NC_DOUBLE, 1, &dimid_ap, &vid_angle_eq));
        NC_CHECK(nc_put_att_text(ncid, vid_angle_eq, "units", 6, "degrees"));

        /* CHARMM UB terms */
        if (GDefaults.iCharmm) {
            NC_CHECK(nc_def_var(ncid, "angle_UB_force_constant", NC_DOUBLE, 1, &dimid_ap, &vid_angle_UB_force));
            NC_CHECK(nc_put_att_text(ncid, vid_angle_UB_force, "units", 19, "kcal/mol/angstrom^2"));
            NC_CHECK(nc_def_var(ncid, "angle_UB_equil_value", NC_DOUBLE, 1, &dimid_ap, &vid_angle_UB_equil));
            NC_CHECK(nc_put_att_text(ncid, vid_angle_UB_equil, "units", 8, "angstrom"));
        }
    }

    /* -15,-16,-17,-17B,-17C- torsion/improper params (+ restraints) */
    int vid_dihe_period, vid_dihe_phase, vid_dihe_force, vid_dihe_scee, vid_dihe_scnb;
    {
        int dimid_tp;
        int nt  = iParmSetTotalTorsionParms(uUnit->psParameters);
        int ni  = iParmSetTotalImproperParms(uUnit->psParameters);
        int nr  = iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION);
        int n   = nt + ni + nr;
        NC_CHECK(nc_def_dim(ncid, "dihedral_types", n, &dimid_tp));
        NC_CHECK(nc_def_var(ncid, "dihedral_force_constant", NC_DOUBLE, 1, &dimid_tp, &vid_dihe_force));
        NC_CHECK(nc_put_att_text(ncid, vid_dihe_force, "units", 8, "kcal/mol"));
        NC_CHECK(nc_def_var(ncid, "dihedral_periodicity", NC_INT, 1, &dimid_tp, &vid_dihe_period));
        NC_CHECK(nc_put_att_text(ncid, vid_dihe_period, "long_name", 47,
                "torsion Fourier periodicity n: cos(n*phi-phase)"));
        NC_CHECK(nc_def_var(ncid, "dihedral_phase", NC_DOUBLE, 1, &dimid_tp, &vid_dihe_phase));
        NC_CHECK(nc_put_att_text(ncid, vid_dihe_phase, "units", 6, "degree"));
        NC_CHECK(nc_def_var(ncid, "dihedral_scee_scale", NC_DOUBLE, 1, &dimid_tp, &vid_dihe_scee));
        NC_CHECK(nc_put_att_text(ncid, vid_dihe_scee, "long_name", 37, "1-4 electrostatic energy scale factor"));
        NC_CHECK(nc_def_var(ncid, "dihedral_scnb_scale", NC_DOUBLE, 1, &dimid_tp, &vid_dihe_scnb));
        NC_CHECK(nc_put_att_text(ncid, vid_dihe_scnb, "long_name", 27, "1-4 VdW energy scale factor"));
    }
    // 18=SOLTY; 19,20 = LJ A,B coeff (above)
    /* -21,-22- bonds inc/excl hydrogen */
    int vid_bh_atoms, vid_bh_parm, vid_bnh_atoms, vid_bnh_parm;
    int (*bond_ah)[2], *bond_pih, (*bond_anh)[2], *bond_pinh;
    {
        int nb = iVarArrayElementCount(uUnit->vaBonds);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
        int nh = 0, nnh = 0;
        bond_ah  = malloc(2*nb*sizeof(int));
        bond_pih  = malloc(nb*sizeof(int));
        bond_anh = malloc(2*(nb+nr)*sizeof(int));
        bond_pinh= malloc((nb+nr)*sizeof(int));
        for (int i = 0; i < nb; i++) {
            SAVEBONDt *sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, i);
            ATOMt *aA = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1-1)->aAtom;
            ATOMt *aD = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2-1)->aAtom;
            if (iAtomElement(aA)==HYDROGEN || iAtomElement(aD)==HYDROGEN) {
                bond_ah[nh][0]=sbPBond->iAtom1;
                bond_ah[nh][1]=sbPBond->iAtom2;
                bond_pih[nh]=sbPBond->iParmIndex; nh++;
            } else {
                bond_anh[nnh][0]=sbPBond->iAtom1;
                bond_anh[nnh][1]=sbPBond->iAtom2;
                bond_pinh[nnh]=sbPBond->iParmIndex; nnh++;
            }
        }
        if (iVarArrayElementCount(uUnit->vaRestraints)) {
            SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
            for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++)
                if (sr->iType == RESTRAINTBOND) {
                    bond_anh[nnh][0]=sr->iAtom1; bond_anh[nnh][1]=sr->iAtom2;
                    bond_pinh[nnh]=sr->iParmIndex; nnh++;
                }
        }
        int dimid_bh, dimid_bnh;
        NC_CHECK(nc_def_dim(ncid, "bonds_inc_h", nh,  &dimid_bh));
        NC_CHECK(nc_def_dim(ncid, "bonds_excl_h", nb-nh, &dimid_bnh));
        NC_CHECK(nc_def_dim(ncid, "bonds_w_restraints", nnh, &dimid_bnh));

        int dims2[2] = {dimid_bh, dimid_atom2};
        NC_CHECK(nc_def_var(ncid, "bonds_inc_hydrogen_atoms", NC_INT, 2, dims2, &vid_bh_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_bh_atoms, "long_name",          20, "1-based atom indices"));
        NC_CHECK(nc_put_att_text(ncid, vid_bh_atoms, "reference_dimension", 5, "atoms"));

        NC_CHECK(nc_def_var(ncid, "bonds_inc_hydrogen_parm", NC_INT, 1, &dimid_bh, &vid_bh_parm));
        NC_CHECK(nc_put_att_text(ncid, vid_bh_parm, "long_name", 18, "1-based parm index"));
        NC_CHECK(nc_put_att_text(ncid, vid_bh_parm, "reference_dimenson", 6, "bonds"));
        NC_CHECK(nc_put_att_text(ncid, vid_bh_parm, "reference", 36,
                                 "bond_force_constant bond_equil_value"));

        dims2[0]=dimid_bnh;
        NC_CHECK(nc_def_var(ncid, "bonds_without_hydrogen_atoms", NC_INT, 2, dims2, &vid_bnh_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_bnh_atoms, "long_name",          20, "1-based atom indices"));
        NC_CHECK(nc_put_att_text(ncid, vid_bnh_atoms, "reference_dimension", 5, "atoms"));

        NC_CHECK(nc_def_var(ncid, "bonds_without_hydrogen_parm", NC_INT, 1, &dimid_bnh, &vid_bnh_parm));
        NC_CHECK(nc_put_att_text(ncid, vid_bnh_parm, "long_name", 18, "1-based parm index"));
        NC_CHECK(nc_put_att_text(ncid, vid_bnh_parm, "reference_dimenson", 6, "bonds"));
        NC_CHECK(nc_put_att_text(ncid, vid_bnh_parm, "reference", 36,
                                 "bond_force_constant bond_equil_value"));

    }

    /* -23,-24- angles inc/excl hydrogen */
    int vid_ah_atoms, vid_ah_parm, vid_anh_atoms, vid_anh_parm;
    int (*angle_ah)[3], *angle_pih, (*angle_anh)[3], *angle_pinh;
    {
        int na = iVarArrayElementCount(uUnit->vaAngles);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
        angle_ah  = malloc(3*na*sizeof(int));
        angle_anh = malloc(3*(na+nr)*sizeof(int));
        angle_pih      = malloc(na*sizeof(int));
        angle_pinh     = malloc((na+nr)*sizeof(int));
        int nh = 0, nnh = 0;
        for (int i = 0; i < na; i++) {
            SAVEANGLEt *saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
            ATOMt *aA = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1-1)->aAtom;
            ATOMt *aB = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2-1)->aAtom;
            ATOMt *aD = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3-1)->aAtom;
            if (iAtomElement(aA)==HYDROGEN||iAtomElement(aB)==HYDROGEN||
                iAtomElement(aD)==HYDROGEN) {
                angle_ah[nh][0]=saPAngle->iAtom1; angle_ah[nh][1]=saPAngle->iAtom2;
                angle_ah[nh][2]=saPAngle->iAtom3; angle_pih[nh]=saPAngle->iParmIndex; nh++;
            } else {
                angle_anh[nnh][0]=saPAngle->iAtom1; angle_anh[nnh][1]=saPAngle->iAtom2;
                angle_anh[nnh][2]=saPAngle->iAtom3; angle_pinh[nnh]=saPAngle->iParmIndex; nnh++;
            }
        }
        if (iVarArrayElementCount(uUnit->vaRestraints)) {
            SAVERESTRAINTt *sr = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
            for (int i = 0; i < iVarArrayElementCount(uUnit->vaRestraints); i++, sr++)
                if (sr->iType == RESTRAINTANGLE) {
                    angle_anh[nnh][0]=sr->iAtom1; angle_anh[nnh][1]=sr->iAtom2;
                    angle_anh[nnh][2]=sr->iAtom3; angle_pinh[nnh]=sr->iParmIndex; nnh++;
                }
        }
        int dimid_ah, dimid_anh;
        int dims2[2] = { dimid_ah, dimid_atom3 };

        NC_CHECK(nc_def_dim(ncid, "angles_inc_h", nh,  &dimid_ah));
        NC_CHECK(nc_def_dim(ncid, "angles_excl_h", na-nh, &dimid_anh));
        NC_CHECK(nc_def_dim(ncid, "angles_w_restraints", nnh, &dimid_anh));

        NC_CHECK(nc_def_var(ncid, "angles_inc_hydrogen_atoms", NC_INT, 2, dims2, &vid_ah_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_ah_atoms, "long_name",          20, "1-based atom indices"));
        NC_CHECK(nc_put_att_text(ncid, vid_ah_atoms, "reference_dimension", 5, "atoms"));

        NC_CHECK(nc_def_var(ncid, "angles_inc_hydrogen_parm", NC_INT, 1, &dimid_ah, &vid_ah_parm));
        NC_CHECK(nc_put_att_text(ncid, vid_ah_parm, "long_name", 18, "1-based parm index"));
        NC_CHECK(nc_put_att_text(ncid, vid_ah_parm, "reference_dimenson", 6, "angles"));
        NC_CHECK(nc_put_att_text(ncid, vid_ah_parm, "reference", 38,
                                 "angle_force_constant angle_equil_value"));

        dims2[0]=dimid_anh;
        NC_CHECK(nc_def_var(ncid, "angles_without_hydrogen_atoms", NC_INT, 2, dims2, &vid_anh_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_anh_atoms, "long_name",          20, "1-based atom indices"));
        NC_CHECK(nc_put_att_text(ncid, vid_anh_atoms, "reference_dimension", 5, "atoms"));

        NC_CHECK(nc_def_var(ncid, "angles_without_hydrogen_parm", NC_INT, 1, &dimid_anh, &vid_anh_parm));
        NC_CHECK(nc_put_att_text(ncid, vid_anh_parm, "long_name", 18, "1-based parm index"));
        NC_CHECK(nc_put_att_text(ncid, vid_anh_parm, "reference_dimenson", 6, "angles"));
        NC_CHECK(nc_put_att_text(ncid, vid_anh_parm, "reference", 38,
                                 "angle_force_constant angle_equil_value"));
    }

    /* -25,-26- torsions */
    int vid_dih_h_atoms, vid_dih_h_parm;
    int vid_dih_nh_atoms, vid_dih_nh_parm;
    int (*dihe_ah)[4], *dihe_pih, (*dihe_anh)[4], *dihe_pinh;
    {
        int nh = 0, nnh = 0;
        int dimid_dih_h, dimid_dih_nh;
        int nT = iVarArrayElementCount(uUnit->vaTorsions);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION);
        dihe_ah  = malloc(4*(nT+nr)*sizeof(int));
        dihe_pih  = malloc(nT*sizeof(int));
        dihe_anh = malloc(4*nT*sizeof(int));
        dihe_pinh = malloc(nT*sizeof(int));

        for (int i = 0; i < nT; i++) {
            SAVETORSIONt *t = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
            ATOMt *aA = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom1-1)->aAtom;
            ATOMt *aB = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom2-1)->aAtom;
            ATOMt *aC = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom3-1)->aAtom;
            ATOMt *aD = PVAI(uUnit->vaAtoms, SAVEATOMt, t->iAtom4-1)->aAtom;
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
                dihe_ah[nh][0]=t->iAtom1; dihe_ah[nh][1]=t->iAtom2;
                dihe_ah[nh][2]=t->iAtom3*iCalc14; dihe_ah[nh][3]=t->iAtom4*iProper;
                dihe_pih[nh]=t->iParmIndex; nh++;
            } else {
                dihe_anh[nnh][0]=t->iAtom1; dihe_anh[nnh][1]=t->iAtom2;
                dihe_anh[nnh][2]=t->iAtom3*iCalc14; dihe_anh[nnh][3]=t->iAtom4*iProper;
                dihe_pinh[nnh]=t->iParmIndex; nnh++;
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
                    dihe_anh[nnh][0]=srPRestraint->iAtom1;
                    dihe_anh[nnh][1]=srPRestraint->iAtom2;
                    dihe_anh[nnh][2]=srPRestraint->iAtom3;
                    dihe_anh[nnh][3]=srPRestraint->iAtom4;
                    dihe_pinh[nnh]=srPRestraint->iParmIndex;
                    nnh++;
                }
            }
        }

        NC_CHECK(nc_def_dim(ncid, "dihedrals_inc_h", nh,  &dimid_dih_h));
        NC_CHECK(nc_def_dim(ncid, "dihedrals_excl_h", nT-nh, &dimid_dih_nh));
        NC_CHECK(nc_def_dim(ncid, "dihedrals_w_restraints", nnh, &dimid_dih_nh));

        int dims2[2]={ dimid_dih_h, dimid_atom4 };
        NC_CHECK(nc_def_var(ncid, "dihedrals_inc_hydrogen_atoms", NC_INT, 2, dims2, &vid_dih_h_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_dih_h_atoms, "long_name", 20, "1-based atom indices"));
        NC_CHECK(nc_put_att_text(ncid, vid_dih_h_atoms, "reference_dimension", 5, "atoms"));
        NC_CHECK(nc_put_att_text(ncid, vid_dih_h_atoms, "note", 46,
                                 "atom3,atom4 sign encodes bCalc14,bProper flags"));

        NC_CHECK(nc_def_var(ncid, "dihedrals_inc_hydrogen_parm", NC_INT, 1, &dimid_dih_h, &vid_dih_h_parm));
        NC_CHECK(nc_put_att_text(ncid, vid_dih_h_parm, "long_name", 18, "1-based parm index"));
        NC_CHECK(nc_put_att_text(ncid, vid_dih_h_parm, "reference",
                                 59, "dihedral_force_constant dihedral_periodicity dihedral_phase"));

        dims2[0]=dimid_dih_nh;
        NC_CHECK(nc_def_var(ncid, "dihedrals_without_hydrogen_atoms", NC_INT, 2, dims2, &vid_dih_nh_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_dih_nh_atoms, "long_name", 20, "1-based atom indices"));
        NC_CHECK(nc_put_att_text(ncid, vid_dih_nh_atoms, "reference_dimension", 5, "atoms"));
        NC_CHECK(nc_put_att_text(ncid, vid_dih_nh_atoms, "note", 47,
                                 "atom3/atom4 sign encodes bCalc14/bProper flags"));

        NC_CHECK(nc_def_var(ncid, "dihedrals_without_hydrogen_parm", NC_INT, 1, &dimid_dih_nh, &vid_dih_nh_parm));
        NC_CHECK(nc_put_att_text(ncid, vid_dih_nh_parm, "long_name", 18, "1-based parm index"));
        NC_CHECK(nc_put_att_text(ncid, vid_dih_nh_parm, "reference",
                                 59, "dihedral_force_constant dihedral_periodicity dihedral_phase"));
    }

    /* -27- excluded atoms list */
    int vid_nbex_list;
    {
        int n = iVarArrayElementCount(vaExcludedAtoms);
        int dimid_excl;
        NC_CHECK(nc_def_dim(ncid, "nonbond_excl_list_size", n, &dimid_excl));
        NC_CHECK(nc_def_var(ncid, "excluded_atoms_list", NC_INT, 1, &dimid_excl, &vid_nbex_list));
        NC_CHECK(nc_put_att_text(ncid, vid_nbex_list, "long_name",          19, "excluded atom list"));
        NC_CHECK(nc_put_att_text(ncid, vid_nbex_list, "reference_dimension", 5, "atoms"));
    }

    /* -28,-29,-30- hbond (skip if none) */
    int  vid_hba, vid_hbb;
    {
        int n = iParmSetTotalHBondParms(uUnit->psParameters);
        if (n > 0) {
            int dimid_hb;
            NC_CHECK(nc_def_dim(ncid, "hbond_types", n, &dimid_hb));
            NC_CHECK(nc_def_var(ncid, "hbond_acoef", NC_DOUBLE, 1, &dimid_hb, &vid_hba));
            NC_CHECK(nc_put_att_text(ncid, vid_hba, "units", 20, "kcal/mol*angstrom^12"));
            NC_CHECK(nc_def_var(ncid, "hbond_bcoef", NC_DOUBLE, 1, &dimid_hb, &vid_hbb));
            NC_CHECK(nc_put_att_text(ncid, vid_hbb, "units", 20, "kcal/mol*angstrom^10"));
        }
    }

    /* -31- amber atom type strings */
    int vid_atom_type;
    {
        int dims2[2] = {dimid_atoms, dimid_name4};
        NC_CHECK(nc_def_var(ncid, "atom_amber_type", NC_CHAR, 2, dims2, &vid_atom_type));
        NC_CHECK(nc_put_att_text(ncid, vid_atom_type, "long_name", 26, "force field atom type name"));
    }

    /* -32- tree chain classification */
    int vid_atom_tree;
    {
        int dims2[2] = { dimid_atoms, dimid_name4};
        NC_CHECK(nc_def_var(ncid, "atom_tree_chain_classification", NC_CHAR, 2, dims2, &vid_atom_tree));
        NC_CHECK(nc_put_att_text(ncid, vid_atom_tree, "long_name", 52, "mainchain tree class: M,E,S,B,3,4,5,6,X; BLA=unknown"));
    }

    /* -35A..C- solvent/box info (conditional) */
    int vid_mol, vid_fsr, vid_fsm;
    int nMol, iFS1, iFirstSolvRes;
    if (bUnitUseBox(uUnit)) {
        iFirstSolvRes = n_residues;
        for (int i = 0; i < n_residues; i++) {
            //if (PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sResidueType[0] == RESTYPESOLVENT) {
            if ( bResidueFlagsSet(PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->rResidue, RESIDUEBULKSOLVENT) ) {
                iFirstSolvRes = i; break;
            }
        }
        UnitIOFindAndCountMolecules(uUnit);
        nMol = iVarArrayElementCount(uUnit->vaAtomsPerMolecule);
        iFS1 = uUnit->iFirstSolvent + 1;
        int dimid_mol;

        NC_CHECK(nc_def_dim(ncid, "molecules", nMol, &dimid_mol));
        NC_CHECK(nc_def_var(ncid, "atoms_per_molecule", NC_INT, 1, &dimid_mol, &vid_mol));
        NC_CHECK(nc_put_att_text(ncid, vid_mol, "long_name", 30, "number of atoms in molecule N"));
        NC_CHECK(nc_def_var(ncid, "first_solvent_residue", NC_INT, 0, NULL, &vid_fsr));
        NC_CHECK(nc_put_att_text(ncid, vid_fsr, "long_name",          46,
                                 "1-based residue index of first solvent residue"));
        NC_CHECK(nc_put_att_text(ncid, vid_fsr, "reference_dimension", 4, "nres"));
        NC_CHECK(nc_def_var(ncid, "first_solvent_molecule", NC_INT, 0, NULL, &vid_fsm));
        NC_CHECK(nc_put_att_text(ncid, vid_fsm, "long_name",          48,
                                 "1-based molecule index of first solvent molecule"));
        NC_CHECK(nc_put_att_text(ncid, vid_fsm, "reference_dimension", 4, "molecules"));
    }

    /* -35D,-35E- cap info (conditional) */
    int vid_cap_info_atom, vid_cap_info_radius, vid_cap_info_origin;
    if (bUnitUseSolventCap(uUnit)) {
        int dimid_vec;
        NC_CHECK(nc_def_var(ncid, "cap_info_atom", NC_INT, 0, NULL, &vid_cap_info_atom));
        NC_CHECK(nc_put_att_text(ncid, vid_cap_info_atom, "long_name",          45,
                             "Index of the last atom not in the Solvent Cap"));
        NC_CHECK(nc_put_att_text(ncid, vid_cap_info_atom, "reference_dimension", 5, "atoms"));
        NC_CHECK(nc_def_var(ncid, "cap_info_radius", NC_DOUBLE, 0, NULL, &vid_cap_info_radius));
        NC_CHECK(nc_put_att_text(ncid, vid_cap_info_radius, "long_name",          21,
                             "Radius of Solvent Cap"));
        NC_CHECK(nc_put_att_text(ncid, vid_cap_info_radius, "units", 8, "angstrom"));
        NC_CHECK(nc_def_dim(ncid, "vector", 3, &dimid_vec));
        NC_CHECK(nc_def_var(ncid, "cap_info_origin", NC_DOUBLE, 1, &dimid_vec, &vid_cap_info_origin));
        NC_CHECK(nc_put_att_text(ncid, vid_cap_info_origin, "units", 8, "angstrom"));
        NC_CHECK(nc_put_att_text(ncid, vid_cap_info_radius, "long_name", 21,
                             "Origin of Solvent Cap"));
    }

    /* GB radius set - global attribute */
    {
        STRING sTemp;
        const char *rset="Unknown radius set";
        if (GDefaults.iGBparm >= 0 && GDefaults.iGBparm <= GBPARM_MAX) {
            sprintf(sTemp,"%s (%s)",PBRadii_optionDesc[GDefaults.iGBparm],
                      PBRadii_options[GDefaults.iGBparm]);
            rset = sTemp;
        }
        NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "radius_set", strlen(rset), rset));
    }

    /* IPOL - true scalar */
    int vid_IPOL;
    {
        NC_CHECK(nc_def_var(ncid, "polarizablity_type", NC_INT, 0, NULL, &vid_IPOL));
        NC_CHECK(nc_put_att_text(ncid, vid_IPOL, "long_name", 28, "polarizable force field flag"));
    }

    /* polar */
    int vid_polar, vid_damp;
    if (bPolar) {
        NC_CHECK(nc_def_var(ncid, "polarizability", NC_DOUBLE, 1, &dimid_atoms, &vid_polar));
        NC_CHECK(nc_put_att_text(ncid, vid_polar, "units", 10, "angstrom^3"));
        NC_CHECK(nc_def_var(ncid, "dipole_damp_factor", NC_DOUBLE, 1, &dimid_atoms, &vid_damp));
    }

    /* PDB chain info */
    int vid_resSeq, vid_res_chain;
    if (GDefaults.bPdbKeepChainId) {
        NC_CHECK(nc_def_var(ncid, "residue_number", NC_INT, 1, &dimid_res, &vid_resSeq));
        NC_CHECK(nc_put_att_text(ncid, vid_resSeq, "long_name", 10, "PDB resSeq"));
        NC_CHECK(nc_put_att_text(ncid, vid_resSeq, "pdb_field",     6, "resSeq"));
        if (GDefaults.bCIFReadAuth)
            NC_CHECK(nc_put_att_text(ncid, vid_resSeq, "mmcif_field",    21, "atom_site.auth_seq_id"));
        else NC_CHECK(nc_put_att_text(ncid, vid_resSeq, "mmcif_field",  22, "atom_site.label_seq_id"));

        int dims2[2] = { dimid_res, dimid_name4 };
        NC_CHECK(nc_def_var(ncid, "residue_chainid", NC_CHAR, 2, dims2, &vid_res_chain));
        NC_CHECK(nc_put_att_text(ncid, vid_res_chain, "long_name", 10, "PDB resSeq"));
        NC_CHECK(nc_put_att_text(ncid, vid_res_chain, "pdb_field",     7, "chainId"));
        if (GDefaults.bCIFReadAuth)
            NC_CHECK(nc_put_att_text(ncid, vid_res_chain, "mmcif_field",    22, "atom_site.auth_comp_id"));
        else NC_CHECK(nc_put_att_text(ncid, vid_res_chain, "mmcif_field",  23, "atom_site.label_comp_id"));
    }

    /* CMAP -- pending data structure review */
    int vid_phipsi_atoms, vid_phipsi_mapidx;
    int *vid_cmap_grid = NULL;
    int CMAP_TYPES = iParmSetTotalCMAPParms(uUnit->psParameters);
    int CMAP_KINDS = iVarArrayElementCount(uUnit->vaCMAPs); // = NumPhiPsi
    if (CMAP_TYPES && CMAP_KINDS && GDefaults.iCMAP) {

        /* ── write CMAP_KINDS and CMAP_TYPES scalars ── */
        /* ── CMAP_INDEX_ATOMS and CMAP_INDEX_MAP ── */

            int dimid_types, dimid_kinds, dimid_quint;
            NC_CHECK(nc_def_dim(ncid, "phi_psi_types", CMAP_TYPES, &dimid_types));
            NC_CHECK(nc_def_dim(ncid, "phi_psi", CMAP_KINDS, &dimid_kinds));
            NC_CHECK(nc_def_dim(ncid, "atom_quint", 5, &dimid_quint));
            int dims2[2] = { dimid_kinds, dimid_quint };
            NC_CHECK(nc_def_var(ncid, "phi_psi_atoms", NC_INT, 2, dims2, &vid_phipsi_atoms));
            NC_CHECK(nc_put_att_text(ncid, vid_phipsi_atoms, "long_name", 50,
                                     "5 raw 1-based atom indices: atoms 1-4=phi, 2-5=psi"));
            NC_CHECK(nc_put_att_text(ncid, vid_phipsi_atoms, "reference_dimension", 6, "atoms"));

            NC_CHECK(nc_def_var(ncid, "phi_psi_parm", NC_INT, 1, &dimid_kinds, &vid_phipsi_mapidx));
            NC_CHECK(nc_put_att_text(ncid, vid_phipsi_mapidx, "long_name", 79,
                    "1-based index into phi_psi_cmap_resolution_<nn> and phi_psi_cmap_<nn>"));
            NC_CHECK(nc_put_att_text(ncid, vid_phipsi_atoms, "reference_dimenson", 10, "phi_psi_types"));

        /* ── subgroups — one per active map type ──
           Each contains coordinate variables phi and psi (angle values
           in degrees) plus the GRID energy surface.                     */
        
        vid_cmap_grid = malloc(sizeof(int)*CMAP_TYPES);
        for (int i = 0; i < CMAP_TYPES; i++) {
            int dimid_cmap_res;
            CMAPt cmap;
            ParmSetCMAP(uUnit->psParameters, i, &cmap, FALSE);

            /* RESOLUTION_nn dimension shared by both phi and psi axes. */
            STRING sResolution;
            sprintf(sResolution,"phi_psi_cmap_resolution_%02d",i+1);
            NC_CHECK(nc_def_dim(ncid, sResolution, cmap.resolution, &dimid_cmap_res));

            /* GRID(RESOLUTION, RESOLUTION) — same dim for both axes */
            dims2[0] = dims2[1] = dimid_cmap_res;
            STRING sParam;
            sprintf(sParam,"phi_psi_cmap_%02d",i+1);
            NC_CHECK(nc_def_var(ncid, sParam, NC_DOUBLE, 2, dims2, &vid_cmap_grid[i]));
            NC_CHECK(nc_put_att_text(ncid, vid_cmap_grid[i], "units",      8, "kcal/mol"));
            NC_CHECK(nc_put_att_text(ncid, vid_cmap_grid[i], "long_name",  36,
                                     "CMAP energy grid: rows=phi, cols=psi"));
            NC_CHECK(nc_put_att_text(ncid, vid_cmap_grid[i], "coordinates", 26, sResolution));
            NC_CHECK(nc_put_att_text(ncid, vid_cmap_grid[i], "title", strlen(cmap.title), cmap.title));
        }
    }


    NC_CHECK(nc_enddef(ncid));
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /* cell dimensions */
    if (bUnitUseBox(uUnit)) {
        double dX, dY, dZ;
        double dAlpha = uUnit->dAlpha / DEGTORAD;
        double dBeta  = uUnit->dBeta  / DEGTORAD;
        double dGamma = uUnit->dGamma / DEGTORAD;
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
        int n = iVarArrayElementCount((uUnit)->vaAtoms);
        char *name = malloc(n * LABLEN);
        memset(name, ' ', n * LABLEN);
        for (int i = 0; i < n; i++) {
            const char *s = PVAI((uUnit)->vaAtoms,SAVEATOMt,i)->sName;
            size_t l = strlen(s);
            if (l > (size_t)LABLEN) s += l - LABLEN;
            memcpy(name+i*LABLEN, s, l < (size_t)LABLEN ? l : LABLEN);
        }
        NC_CHECK(nc_put_var_text(ncid, vid_name, name));
        free(name);
    }

    /* -4- charges (raw double, no AMBER unit scaling on write) */
    {
        double *buf = malloc(n_atoms * sizeof(double));
        for (int i = 0; i < n_atoms; i++)
            buf[i] = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->dCharge;
        NC_CHECK(nc_put_var_double(ncid, vid_charge, buf));
        free(buf);
    }

    /* -4b- atomic numbers */
    {
        int *buf = malloc(n_atoms * sizeof(int));
        for (int i = 0; i < n_atoms; i++)
            buf[i] = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iElement;
        NC_CHECK(nc_put_var_int(ncid, vid_atomic_num, buf));
        free(buf);
    }

    /* -5-, RADII, SCREEN: merged single ParmSetAtom pass */
    {
        int n = iVarArrayElementCount(uUnit->vaAtoms);
        double *mass   = malloc(n * sizeof(double));
        double *radii  = malloc(n * sizeof(double));
        double *screen = malloc(n * sizeof(double));
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
            radii[i]  = dGBRadiusForAtom(&sa[i], iElement, dMass, (i+1==n) );
            screen[i] = dGBScreenForElement(iElement);
        }
        NC_CHECK(nc_put_var_double(ncid, vid_mass, mass));
        NC_CHECK(nc_put_var_double(ncid, vid_radii, radii));
        NC_CHECK(nc_put_var_double(ncid, vid_screen, screen));
        free(mass); free(radii); free(screen);
    }

    /* -6- atom type index */
    {
        int *buf = malloc(n_atoms * sizeof(int));
        for (int i = 0; i < n_atoms; i++) {
            int iAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iTypeIndex - 1;
            buf[i] = *PVAI(vaNBIndex, int, iAtom) + 1;
        }
        NC_CHECK(nc_put_var_int(ncid, vid_itype, buf));
        free(buf);
    }

    /* -7- number of excluded atoms per atom */
    NC_CHECK(nc_put_var_int(ncid, vid_excl, PVAI(vaExcludedCount, int, 0)));

    /* -8- nonbonded parm index matrix */
    NC_CHECK(nc_put_var_int(ncid, vid_nbidx, PVAI(vaNBIndexMatrix, int, 0)));

    /* -9- residue labels */
    {
        char *buf = malloc(n_residues * LABLEN);
        char *buf_type = malloc(n_residues);
        int iUnknownTypes = 0;
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
        NC_CHECK(nc_put_var_text(ncid, vid_resname, buf));
        NC_CHECK(nc_put_var_text(ncid, vid_res_type, buf_type));
        free(buf_type); free(buf);
    }

    /* -10- residue pointer */
    {
        int *buf = malloc(n_residues * sizeof(int));
        for (int i = 0; i < n_residues; i++)
            buf[i] = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iAtomStartIndex;
        NC_CHECK(nc_put_var_int(ncid, vid_res_first_atom, buf));
        free(buf);
    }

    /* -11,-12A- bond force constants and equil (+ restraints) */
    {
        int nb = iParmSetTotalBondParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
        int n  = nb + nr;
        double *kb = malloc(n*sizeof(double)), *r0 = malloc(n*sizeof(double));
        for (int i = 0; i < nb; i++) {
            ParmSetBond(uUnit->psParameters, i, s1, s2, &dKb, &dR0,
                        &dKpull, &dRpull0, &dKpress, &dRpress0, sDesc);
            kb[i] = dKb; r0[i] = dR0;
        }
        NC_RESTRAINTLOOP(uUnit, RESTRAINTBOND, dKx, kb, nb);
        NC_CHECK(nc_put_var_double(ncid, vid_bond_force, kb));
        NC_RESTRAINTLOOP(uUnit, RESTRAINTBOND, dX0, r0, nb);
        NC_CHECK(nc_put_var_double(ncid, vid_bond_eq, r0));
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
                        &dKpull, &dRpull0, &dKpress, &dRpress0, sDesc);
            kpull[i]=dKpull; rpull0[i]=dRpull0;
            kpress[i]=dKpress; rpress0[i]=dRpress0;
        }
        for (int i = nb; i < n; i++) {
            kpull[i]=0.0; rpull0[i]=100.0;
            kpress[i]=0.0; rpress0[i]=-100.0;
        }
        NC_CHECK(nc_put_var_double(ncid, vid_bond_stiffness_pull_adj, kpull));
        NC_CHECK(nc_put_var_double(ncid, vid_bond_equil_pull_adj, rpull0));
        NC_CHECK(nc_put_var_double(ncid, vid_bond_stiffness_press_adj, kpress));
        NC_CHECK(nc_put_var_double(ncid, vid_bond_equil_press_adj, rpress0));
        free(kpull); free(rpull0); free(kpress); free(rpress0);
    }

    /* -13,-14- angle force constants and equil (+ restraints) */
    {
        int na = iParmSetTotalAngleParms(uUnit->psParameters);
        int nr = iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
        int n  = na + nr;
        double *kt = malloc(n*sizeof(double)), *t0 = malloc(n*sizeof(double));
        for (int i = 0; i < na; i++) {
            ParmSetAngle(uUnit->psParameters, i, s1, s2, s3, &dKt, &dT0,
                         &dTkub, &dRkub, sDesc);
            kt[i] = dKt; t0[i] = dT0 / DEGTORAD;
        }
        NC_RESTRAINTLOOP(uUnit, RESTRAINTANGLE, dKx, kt, na);
        NC_CHECK(nc_put_var_double(ncid, vid_angle_force, kt));

        NC_RESTRAINTLOOP(uUnit, RESTRAINTANGLE, dX0, t0, na);
        NC_CHECK(nc_put_var_double(ncid, vid_angle_eq, t0));
        free(kt); free(t0);

        /* CHARMM UB terms */
        if (GDefaults.iCharmm) {
            double *tkub = malloc(na*sizeof(double)), *rkub = malloc(na*sizeof(double));
            for (int i = 0; i < na; i++) {
                ParmSetAngle(uUnit->psParameters, i, s1, s2, s3, &dKt, &dT0, &dTkub, &dRkub, sDesc);
                tkub[i] = dTkub; rkub[i] = dRkub;
            }
            NC_CHECK(nc_put_var_double(ncid, vid_angle_UB_force, tkub));
            NC_CHECK(nc_put_var_double(ncid, vid_angle_UB_equil, rkub));
            free(tkub); free(rkub);
        }
    }

    /* -15,-16,-17,-17B,-17C- torsion/improper params (+ restraints) */
    {
        int iN;
        int nt  = iParmSetTotalTorsionParms(uUnit->psParameters);
        int ni  = iParmSetTotalImproperParms(uUnit->psParameters);
        int nr  = iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION);
        int n   = nt + ni + nr;
        double *kp   = malloc(n*sizeof(double));
        int *per  = malloc(n*sizeof(int));
        double *p0   = malloc(n*sizeof(double));
        double *scee = malloc(n*sizeof(double));
        double *scnb = malloc(n*sizeof(double));

        for (int i = 0; i < nt; i++) {
            ParmSetTorsion(uUnit->psParameters, i, s1, s2, s3, s4,
                           &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc);
            kp[i]=dKp; per[i]=iN; p0[i]=dP0/DEGTORAD;
            scee[i]=(dScEE<0.0)?GDefaults.dSceeScaleFactor:dScEE;
            scnb[i]=(dScNB<0.0)?GDefaults.dScnbScaleFactor:dScNB;
        }
        for (int i = 0; i < ni; i++) {
            ParmSetImproper(uUnit->psParameters, i, s1, s2, s3, s4,
                            &iN, &dKp, &dP0, sDesc);
            kp[nt+i]=dKp; per[nt+i]=iN; p0[nt+i]=dP0/DEGTORAD;
            scee[nt+i]=0.0; scnb[nt+i]=0.0;
        }
        NC_RESTRAINTLOOP(uUnit, RESTRAINTTORSION, dKx, kp,  nt+ni);
        NC_RESTRAINTLOOP(uUnit, RESTRAINTTORSION, dX0, p0,  nt+ni);
        NC_RESTRAINTLOOP(uUnit, RESTRAINTTORSION, dN,  per, nt+ni);
        NC_CHECK(nc_put_var_double(ncid, vid_dihe_force, kp));
        NC_CHECK(nc_put_var_int(ncid, vid_dihe_period, per));
        NC_CHECK(nc_put_var_double(ncid, vid_dihe_phase, p0));
        NC_CHECK(nc_put_var_double(ncid, vid_dihe_scee, scee));
        NC_CHECK(nc_put_var_double(ncid, vid_dihe_scnb, scnb));
        free(kp); free(per); free(p0); free(scee); free(scnb);
    }

    /* -19,-20- LJ coefficients */
    {
        double *A = malloc(nttyp*sizeof(double)), *B = malloc(nttyp*sizeof(double));
        for (int i = 0; i < nttyp; i++) {
            A[i] = PVAI(vaNBParameters, NONBONDACt, i)->dA;
            B[i] = PVAI(vaNBParameters, NONBONDACt, i)->dB;
        }
        NC_CHECK(nc_put_var_double(ncid, vid_acoef, A));
        NC_CHECK(nc_put_var_double(ncid, vid_bcoef, B));

        /* CHARMM 1-4 terms */
        if (GDefaults.iCharmm) {
            for (int i = 0; i < nttyp; i++) {
                A[i] = PVAI(vaNBParameters, NONBONDACt, i)->dA14;
                B[i] = PVAI(vaNBParameters, NONBONDACt, i)->dB14;
            }
            NC_CHECK(nc_put_var_double(ncid, vid_14acoef, A));
            NC_CHECK(nc_put_var_double(ncid, vid_14bcoef, B));
        }
        free(A); free(B);
    }


    /* -21,-22- bonds inc/excl hydrogen */
        NC_CHECK(nc_put_var_int(ncid, vid_bh_atoms,  (int*)bond_ah));
        NC_CHECK(nc_put_var_int(ncid, vid_bh_parm,   bond_pih));
        NC_CHECK(nc_put_var_int(ncid, vid_bnh_atoms, (int*)bond_anh));
        NC_CHECK(nc_put_var_int(ncid, vid_bnh_parm,  bond_pinh));
        free(bond_ah); free(bond_pih);
        free(bond_anh); free(bond_pinh);

    /* -23,-24- angles inc/excl hydrogen */
        NC_CHECK(nc_put_var_int(ncid, vid_ah_atoms,  (int*)angle_ah));
        NC_CHECK(nc_put_var_int(ncid, vid_ah_parm,   angle_pih));
        NC_CHECK(nc_put_var_int(ncid, vid_anh_atoms, (int*)angle_anh));
        NC_CHECK(nc_put_var_int(ncid, vid_anh_parm,  angle_pinh));
        free(angle_ah); free(angle_pih);
        free(angle_anh); free(angle_pinh);

    /* -25,-26- torsions */
        NC_CHECK(nc_put_var_int(ncid, vid_dih_h_atoms,  (int*)dihe_ah));
        NC_CHECK(nc_put_var_int(ncid, vid_dih_h_parm,   dihe_pih));
        NC_CHECK(nc_put_var_int(ncid, vid_dih_nh_atoms, (int*)dihe_anh));
        NC_CHECK(nc_put_var_int(ncid, vid_dih_nh_parm,  dihe_pinh));
        free(dihe_ah); free(dihe_pih); free(dihe_anh); free(dihe_pinh);

    /* -27- excluded atoms list */
        NC_CHECK(nc_put_var_int(ncid, vid_nbex_list, PVAI(vaExcludedAtoms, int, 0)));

    /* -28,-29,-30- hbond (skip if none) */
    {
        int n = iParmSetTotalHBondParms(uUnit->psParameters);
        if (n > 0) {
            double *A   = malloc(n*sizeof(double));
            double *B   = malloc(n*sizeof(double));
            for (int i = 0; i < n; i++) {
                ParmSetHBond(uUnit->psParameters, i, s1, s2, &dC, &dD, sDesc);
                A[i] = dC; B[i] = dD;
            }
            NC_CHECK(nc_put_var_double(ncid, vid_hba, A));
            NC_CHECK(nc_put_var_double(ncid, vid_hbb, B));
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
        NC_CHECK(nc_put_var_text(ncid, vid_atom_type, buf));
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
        NC_CHECK(nc_put_var_text(ncid, vid_atom_tree, buf));
        free(buf);
    }

    /* -35A..C- solvent/box info (conditional) */
    if (bUnitUseBox(uUnit)) {
        NC_CHECK(nc_put_var_int(ncid, vid_fsr, &iFirstSolvRes));
        NC_CHECK(nc_put_var_int(ncid, vid_fsm, &iFS1));
        NC_CHECK(nc_put_var_int(ncid, vid_mol, PVAI(uUnit->vaAtomsPerMolecule, int, 0)));
    }

    /* -35D,-35E- cap info (conditional) */
    if (bUnitUseSolventCap(uUnit)) {
        int ires = uUnit->iCapTempInt;
        double dX, dY, dZ, dR;
        UnitGetSolventCap(uUnit, &dX, &dY, &dZ, &dR);
        double cap[3] = { dX, dY, dZ };
        NC_CHECK(nc_put_var_int(ncid, vid_cap_info_atom, &ires));
        NC_CHECK(nc_put_var_double(ncid, vid_cap_info_radius, &dR));
        NC_CHECK(nc_put_var_double(ncid, vid_cap_info_origin, cap));
    }

    /* IPOL - true scalar */
    {
        int val = GDefaults.iIPOL;
        NC_CHECK(nc_put_var_int(ncid, vid_IPOL, &val));
    }

    /* polar */
    if (bPolar) {
        int i, n = iVarArrayElementCount(uUnit->vaAtoms);
        int iCount = 0, iCountPerturbed = 0, iPertTot = 0;
        double *polar = malloc(n * sizeof(double));
        double *damp  = malloc(n * sizeof(double));
        double *pper  = bPert ? malloc(n * sizeof(double)) : NULL;
        char sType[MAXTYPELEN];
        double dMass, dPolar, dEpsilon, dRStar, dEpsilon14, dRStar14, dScreenF;
        int iElement, iHybridization;
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
                if ( bAtomPerturbed(a->aAtom)) {
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
        NC_CHECK(nc_put_var_double(ncid, vid_polar, polar));
        NC_CHECK(nc_put_var_double(ncid, vid_damp, damp));
        free(polar); free(damp);
    }


    /* PDB chain info */
    if (GDefaults.bPdbKeepChainId) {
        int *buf = malloc(n_residues * sizeof(int));
        for (int i = 0; i < n_residues; i++)
            buf[i] = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iPdbResSeq;
        NC_CHECK(nc_put_var_int(ncid, vid_resSeq, buf));
        free(buf);

        char *cbuf = malloc(n_residues * LABLEN);
        memset(cbuf, ' ', n_residues * LABLEN);
        for (int i = 0; i < n_residues; i++) {
            const char *s = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sChainId;
            size_t l = strlen(s);
            if (l > 2) s += l - 2;
            memcpy(cbuf + i*LABLEN, s, l < LABLEN ? l : LABLEN);
        }
        NC_CHECK(nc_put_var_text(ncid, vid_res_chain, cbuf));
        free(cbuf);
    }

    /* CMAP -- pending data structure review */
    if (CMAP_TYPES && CMAP_KINDS && GDefaults.iCMAP) {
    /* ── write CMAP_KINDS and CMAP_TYPES scalars ── */
    /* ── CMAP_INDEX_ATOMS and CMAP_INDEX_MAP ── */

        int *atoms_buf = malloc(CMAP_KINDS * 5 * sizeof(int));
        int *map_buf   = malloc(CMAP_KINDS * sizeof(int));
        for (int i=0, idx=0; i < CMAP_KINDS; i++, idx+=5) {
            SAVECMAPt *saveCMAP = PVAI(uUnit->vaCMAPs, SAVECMAPt, i);
            for (int j=0;j<5;j++) atoms_buf[idx+j]=saveCMAP->iAtom[j];
            map_buf[i]=saveCMAP->iParmIndex;
        }
        NC_CHECK(nc_put_var_int(ncid, vid_phipsi_atoms, atoms_buf));
        NC_CHECK(nc_put_var_int(ncid, vid_phipsi_mapidx,   map_buf));
        free(atoms_buf); free(map_buf);

    /* ── subgroups — one per active map type ──
       Each contains coordinate variables phi and psi (angle values
       in degrees) plus the GRID energy surface.                     */
        for (int i = 0; i < CMAP_TYPES; i++) {
            CMAPt cmap;
            ParmSetCMAP(uUnit->psParameters, i, &cmap, FALSE);
            NC_CHECK(nc_put_var_double(ncid, vid_cmap_grid[i], cmap.map));
        }
        free(vid_cmap_grid);
    }


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
        VP0("Using the unit box\n");
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
    data = (double *)MALLOC(3 * iAtomCount * sizeof(double));
    int counter = 0;
    VECTOR vPos;

    double time = 0.0;
    if (nc_put_var_double(ncid, vid_time, &time) != NC_NOERR)
        VPFATALEXIT("%s: Error writing start time\n", filename);

    double dX, dY, dZ, dX2, dY2, dZ2;
    UnitGetBox(uUnit, &dX, &dY, &dZ);
    if (bUnitUseBox(uUnit) == TRUE && GDefaults.nocenter == 0) {
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
