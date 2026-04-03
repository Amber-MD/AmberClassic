#!/usr/bin/env ccp4-python
"""
supercell_collapse
==================
Convert P1 supercell structure factors (from gemmi sfcalc) into the
reciprocal-space ASU of the smaller primitive cell with a given space group.

Mathematical basis
------------------
Let the supercell be na x nb x nc copies of the primitive cell (all atoms
in P1, no crystallographic symmetry imposed).  The gemmi sfcalc output for
the ASU.pdb scenario -- where the same atoms are placed in the primitive cell
with space group SG -- applies the SG operators to every atom.  The resulting
structure factor at (H,K,L) is:

    F_SG(H,K,L) = Σ_i  exp(2πi H·t_i)  *  F_super( na*H'_i, nb*K'_i, nc*L'_i )

where the sum runs over all symmetry operators i of SG, and

    (H'_i, K'_i, L'_i) = R_direct_i^T * (H,K,L)

with R_direct_i^T the TRANSPOSE of the direct-space rotation matrix of
operator i (not the reciprocal-space apply_to_hkl action).

The phase factor  exp(2πi H·t_i)  in terms of gemmi's API:
    = exp( -i * op.phase_shift([H,K,L]) )          (gemmi stores -2π H·t)

If the required supercell index is not in the P1 half-sphere stored in the
MTZ, its Friedel mate is used:  F(-h,-k,-l) = conj(F(h,k,l)).

Usage
-----
  supercell_collapse  super.mtz  SPACEGROUP  super_mult=na,nb,nc  \\
                      [dmin=VALUE]  [outfile=collapsed.mtz]

All arguments are order-independent; keys are case-insensitive.

    super.mtz          P1 supercell MTZ from gemmi sfcalc (FC + PHIC columns)
    SPACEGROUP         target primitive-cell space group (e.g. P21, P212121)
    super_mult=na,nb,nc  supercell multipliers along a, b, c  (e.g. 2,2,2)
    dmin=VALUE         optional extra resolution cutoff in Angstroms
    outfile=FILE       output MTZ path (default: collapsed.mtz)

Dependencies: ccp4-python  (gemmi and numpy are bundled with CCP4 8+)
"""

import sys
import cmath
import numpy as np
import gemmi


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv):
    args = {
        'mtz':     None,
        'sg':      None,
        'na': 1,  'nb': 1,  'nc': 1,
        'dmin':    None,
        'outfile': 'collapsed.mtz',
    }
    for arg in argv[1:]:
        if '=' in arg:
            key, val = arg.split('=', 1)
            key = key.lower().strip()
            val = val.strip()
            if key in ('super_mult', 'mult', 'multipliers'):
                parts = [p.strip() for p in val.split(',')]
                if len(parts) != 3:
                    sys.exit("ERROR: super_mult must be na,nb,nc  (e.g. 2,2,2)")
                args['na'], args['nb'], args['nc'] = int(parts[0]), int(parts[1]), int(parts[2])
            elif key == 'dmin':
                args['dmin'] = float(val)
            elif key in ('outfile', 'output', 'out'):
                args['outfile'] = val
            elif key in ('sg', 'spacegroup', 'space_group', 'symm', 'symmetry'):
                args['sg'] = val
        elif arg.lower().endswith('.mtz') and args['mtz'] is None:
            args['mtz'] = arg
        elif args['sg'] is None:
            sg_test = gemmi.find_spacegroup_by_name(arg)
            if sg_test is not None:
                args['sg'] = arg
            # else ignore (could be a flag we don't understand)
    return args


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

def collapse(super_mtz_path, sg_name, na, nb, nc, dmin, out_path):

    # ------------------------------------------------------------------
    # 1. Load supercell MTZ
    # ------------------------------------------------------------------
    mtz = gemmi.read_mtz_file(super_mtz_path)
    col_labels = [col.label for col in mtz.columns]

    for label in ('H', 'K', 'L', 'FC', 'PHIC'):
        if label not in col_labels:
            sys.exit(f"ERROR: column '{label}' not found in {super_mtz_path}.\n"
                     f"  Available columns: {col_labels}")

    data    = np.array(mtz.array)
    h_i     = col_labels.index('H')
    k_i     = col_labels.index('K')
    l_i     = col_labels.index('L')
    fc_i    = col_labels.index('FC')
    phi_i   = col_labels.index('PHIC')

    hkl_sc  = data[:, [h_i, k_i, l_i]].astype(np.int32)
    F_sc    = (data[:, fc_i].astype(np.float64) *
               np.exp(1j * np.deg2rad(data[:, phi_i].astype(np.float64))))

    print(f"  Supercell reflections : {len(hkl_sc)}")
    print(f"  Supercell cell        : {mtz.cell}")
    print(f"  Supercell multipliers : {na} x {nb} x {nc}")

    # Build lookup dict; include Friedel mates implicitly on lookup.
    sc_dict = {}
    for (h, k, l), Fc in zip(hkl_sc.tolist(), F_sc.tolist()):
        sc_dict[(h, k, l)] = Fc

    def lookup(h, k, l):
        """Return F_super(h,k,l), using Friedel mate conj if needed."""
        v = sc_dict.get((h, k, l))
        if v is not None:
            return v
        v = sc_dict.get((-h, -k, -l))
        if v is not None:
            return v.conjugate()
        return None

    # ------------------------------------------------------------------
    # 2. Filter to divisible supercell indices  -> primitive-cell ASU candidates
    # ------------------------------------------------------------------
    mask = (
        (hkl_sc[:, 0] % na == 0) &
        (hkl_sc[:, 1] % nb == 0) &
        (hkl_sc[:, 2] % nc == 0)
    )
    hkl_prim_candidates = hkl_sc[mask].copy()
    hkl_prim_candidates[:, 0] //= na
    hkl_prim_candidates[:, 1] //= nb
    hkl_prim_candidates[:, 2] //= nc

    print(f"  After divisibility filter : {mask.sum()}")

    # Primitive cell and space group
    sc_cell = mtz.cell
    prim_cell = gemmi.UnitCell(
        sc_cell.a / na, sc_cell.b / nb, sc_cell.c / nc,
        sc_cell.alpha, sc_cell.beta, sc_cell.gamma
    )
    sg = gemmi.find_spacegroup_by_name(sg_name)
    if sg is None:
        sys.exit(f"ERROR: unknown space group '{sg_name}'")

    print(f"  Primitive cell        : {prim_cell}")
    print(f"  Space group           : {sg.xhm()}  ({len(list(sg.operations()))} operators)")

    # Apply optional resolution cutoff
    if dmin is not None:
        keep = np.array([
            prim_cell.calculate_d(hkl.tolist()) >= dmin
            for hkl in hkl_prim_candidates
        ])
        hkl_prim_candidates = hkl_prim_candidates[keep]
        print(f"  After dmin={dmin} A cutoff : {keep.sum()}")

    # Select only ASU reflections (identity path puts them in P1 half-sphere,
    # which is a superset of the target SG ASU).
    asu_obj = gemmi.ReciprocalAsu(sg)
    asu_mask = np.array([
        asu_obj.is_in(hkl.tolist()) for hkl in hkl_prim_candidates
    ])
    hkl_asu = hkl_prim_candidates[asu_mask]
    print(f"  ASU reflections       : {len(hkl_asu)}")

    # ------------------------------------------------------------------
    # 3. Compute F_SG(H,K,L) for each ASU reflection
    #
    #    F_SG(H) = Σ_i  exp(2πi H·t_i)  *  F_super( na*R_direct_i^T * H )
    #
    #    In gemmi's API:
    #      exp(2πi H·t)  =  exp( -i * op.phase_shift([H,K,L]) )
    #      R_direct^T * H uses TRANSPOSED indexing of op.rot:
    #        Hr = (rot[0][0]*H + rot[1][0]*K + rot[2][0]*L) // DEN
    #        Kr = (rot[0][1]*H + rot[1][1]*K + rot[2][1]*L) // DEN
    #        Lr = (rot[0][2]*H + rot[1][2]*K + rot[2][2]*L) // DEN
    # ------------------------------------------------------------------
    print("  Computing structure factors ...")
    ops = list(sg.operations())
    den = gemmi.Op.DEN

    F_out = []
    hkl_good = []
    n_missing = 0
    n_skipped = 0

    for hkl in hkl_asu.tolist():
        H, K, L = int(hkl[0]), int(hkl[1]), int(hkl[2])
        F_sum = 0j
        any_missing = False

        for op in ops:
            rot = op.rot
            # R_direct^T * (H,K,L)  -- note transposed [col][row] indexing
            Hr = (rot[0][0]*H + rot[1][0]*K + rot[2][0]*L) // den
            Kr = (rot[0][1]*H + rot[1][1]*K + rot[2][1]*L) // den
            Lr = (rot[0][2]*H + rot[1][2]*K + rot[2][2]*L) // den

            F_sc_val = lookup(na*Hr, nb*Kr, nc*Lr)
            if F_sc_val is None:
                n_missing += 1
                any_missing = True
                continue

            # exp(2πi H·t) = exp( -i * op.phase_shift([H,K,L]) )
            phase = -op.phase_shift([H, K, L])
            F_sum += F_sc_val * cmath.exp(1j * phase)

        if any_missing:
            n_skipped += 1
        else:
            F_out.append(F_sum)
            hkl_good.append(hkl)

    if n_missing:
        print(f"  WARNING: {n_missing} operator/reflection lookups failed "
              f"(supercell MTZ missing both the index and its Friedel mate)")
        print(f"  WARNING: {n_skipped} reflections skipped due to missing lookups "
              f"(typically resolution-boundary reflections)")
    hkl_asu = np.array(hkl_good, dtype=np.int32) if hkl_good else np.zeros((0, 3), dtype=np.int32)

    F_out = np.array(F_out, dtype=np.complex128)

    # ------------------------------------------------------------------
    # 4. Write output MTZ
    # ------------------------------------------------------------------
    out_mtz = gemmi.Mtz(with_base=False)
    out_mtz.spacegroup = sg
    out_mtz.cell = prim_cell

    base_ds = out_mtz.add_dataset("HKL_base")
    base_ds.wavelength = 0.0
    data_ds = out_mtz.add_dataset("SUPERCELL_COLLAPSE")
    data_ds.wavelength = next(
        (ds.wavelength for ds in mtz.datasets if ds.wavelength > 0), 1.0
    )
    ds_id = data_ds.id

    for label, ctype in [('H', 'H'), ('K', 'H'), ('L', 'H'),
                          ('FC', 'F'), ('PHIC', 'P')]:
        col = out_mtz.add_column(label, ctype)
        col.dataset_id = 0 if ctype == 'H' else ds_id

    n = len(hkl_asu)
    out_data = np.column_stack([
        hkl_asu[:, 0].astype(np.float32),
        hkl_asu[:, 1].astype(np.float32),
        hkl_asu[:, 2].astype(np.float32),
        np.abs(F_out).astype(np.float32),
        np.degrees(np.angle(F_out)).astype(np.float32),
    ])
    out_mtz.set_data(out_data)
    out_mtz.write_to_file(out_path)

    return n, prim_cell, sg


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    args = parse_args(sys.argv)

    if args['mtz'] is None or args['sg'] is None:
        print(__doc__)
        sys.exit("ERROR: must specify both a supercell MTZ and a space group")

    na, nb, nc = args['na'], args['nb'], args['nc']
    print(f"supercell_collapse")
    print(f"  Input MTZ   : {args['mtz']}")
    print(f"  Space group : {args['sg']}")
    print(f"  Multipliers : {na} x {nb} x {nc}")
    if args['dmin']:
        print(f"  dmin cutoff : {args['dmin']} A")
    print(f"  Output MTZ  : {args['outfile']}")
    print()

    n, cell, sg = collapse(
        args['mtz'], args['sg'],
        na, nb, nc,
        args['dmin'],
        args['outfile'],
    )

    print()
    print(f"Done.")
    print(f"  Output MTZ  : {args['outfile']}")
    print(f"  Reflections : {n}")
    print(f"  Space group : {sg.xhm()}")
    print(f"  Cell        : a={cell.a:.4f}  b={cell.b:.4f}  c={cell.c:.4f}  "
          f"alpha={cell.alpha:.3f}  beta={cell.beta:.3f}  gamma={cell.gamma:.3f}")


if __name__ == '__main__':
    main()
