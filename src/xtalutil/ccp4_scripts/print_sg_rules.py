#!/usr/bin/env ccp4-python
"""
print_sg_rules.py
=================
Print a human-readable table of the exact per-operator combination rules
used by supercell_collapse to fold a P1 supercell MTZ into a target space
group ASU.

For each symmetry operator i the table shows:
  • the operator in xyz notation (as gemmi reports it)
  • the re-indexed primitive-cell Miller index  (H'_i, K'_i, L'_i)
    derived from R_direct_i^T · (H,K,L)
  • the translational part t_i
  • the phase factor exp(2πi H·t_i) expressed symbolically, including
    the parity condition (if any) that determines its sign

When a supercell multiplier na,nb,nc is given the P1 lookup index
(na·H'_i, nb·K'_i, nc·L'_i) is also shown.

Usage
-----
  ccp4-python print_sg_rules.py  SPACEGROUP  [super_mult=na,nb,nc]

Examples
--------
  ccp4-python print_sg_rules.py  P212121
  ccp4-python print_sg_rules.py  P21       super_mult=1,2,1
  ccp4-python print_sg_rules.py  C2221     super_mult=2,2,2
"""

import sys
from fractions import Fraction
from math import gcd
import gemmi

DEN = gemmi.Op.DEN   # 24


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def format_hkl_prime(rot, labels=('H', 'K', 'L')):
    """
    Return a string like '-H, K, -L' representing R_direct^T · (H,K,L).

    R_direct^T means: column col of rot gives the coefficients for
    output index col, i.e. row col of R^T.
    """
    parts = []
    for col in range(3):
        coeffs = [rot[row][col] // DEN for row in range(3)]
        terms = []
        for c, lbl in zip(coeffs, labels):
            if c == 0:
                continue
            elif c == 1:
                terms.append(lbl)
            elif c == -1:
                terms.append(f'-{lbl}')
            else:
                terms.append(f'{c}{lbl}')
        if not terms:
            parts.append('0')
        else:
            s = terms[0]
            for t in terms[1:]:
                s += (t if t.startswith('-') else '+' + t)
            parts.append(s)
    return ', '.join(parts)


def format_translation(trn):
    """Return t as a string like '(1/2, 0, 1/2)'."""
    strs = []
    for v in trn:
        if v == 0:
            strs.append('0')
        else:
            g = gcd(abs(v), DEN)
            n, d = v // g, DEN // g
            strs.append(f'{n}/{d}')
    return '(' + ', '.join(strs) + ')'


def phase_factor_symbolic(trn):
    """
    Symbolically simplify exp(2πi · (t_a·H + t_b·K + t_c·L)).

    Returns (exponent_str, condition_str_or_None, factor_str).

    For the common crystallographic case where all t_j ∈ {0, 1/2}:
      the phase is (-1)^(sum of selected indices), which is ±1 depending
      on the parity of that sum.

    For t_j ∈ {0, 1/4, 1/2, 3/4} (less common):
      the phase cycles through {1, i, -1, -i}.

    For general t_j (e.g. rhombohedral, hexagonal):
      a symbolic form is returned without a parity table.
    """
    fracs = [Fraction(v, DEN) for v in trn]
    labels = ('H', 'K', 'L')

    # --- Case 1: all t_j in {0, 1/2}  -->  ±1 factor ---
    if all(f.denominator in (1, 2) for f in fracs):
        int_coeffs = [int(f * 2) for f in fracs]   # 0 or 1
        parity_lbls = [lbl for c, lbl in zip(int_coeffs, labels) if c % 2 != 0]
        if not parity_lbls:
            return 'exp(0) = +1', None, '+1'
        p = '+'.join(parity_lbls)
        cond = (f'({p}) % 2 == 0  ->  +1\n'
                f'({p}) % 2 != 0  ->  -1')
        return f'exp(πi·({p}))', cond, f'(-1)^({p})'

    # --- Case 2: all t_j in {0, 1/4, 1/2, 3/4}  -->  ±1 or ±i factor ---
    if all(f.denominator in (1, 2, 4) for f in fracs):
        qcoeffs = [int(f * 4) % 4 for f in fracs]  # 0,1,2,3
        qterms = []
        for c, lbl in zip(qcoeffs, labels):
            if c == 0:
                continue
            elif c == 1:
                qterms.append(lbl)
            elif c == 2:
                qterms.append(f'2{lbl}')
            elif c == 3:
                qterms.append(f'3{lbl}')
        if not qterms:
            return 'exp(0) = +1', None, '+1'
        q = '+'.join(qterms)
        rows = [f'({q}) % 4 == {n}  ->  {v}'
                for n, v in enumerate(['+1', '+i', '-1', '-i'])]
        cond = '\n'.join(rows)
        return f'exp(πi/2·({q}))', cond, f'i^({q})'

    # --- Case 3: general (hexagonal, rhombohedral, etc.) ---
    # Build the argument in units of 2π
    sym_terms = []
    for f, lbl in zip(fracs, labels):
        if f == 0:
            continue
        if f.denominator == 1:
            sym_terms.append(f'{f.numerator}·{lbl}')
        else:
            sym_terms.append(f'({f.numerator}/{f.denominator})·{lbl}')
    if not sym_terms:
        return 'exp(0) = +1', None, '+1'
    arg = ' + '.join(sym_terms)
    return f'exp(2πi·({arg}))', '(phase varies continuously; evaluate per HKL)', '(general)'


# ---------------------------------------------------------------------------
# Main table printer
# ---------------------------------------------------------------------------

def print_rules(sg_name, na=1, nb=1, nc=1):
    sg = gemmi.find_spacegroup_by_name(sg_name)
    if sg is None:
        sys.exit(f"ERROR: unknown space group '{sg_name}'")

    ops = list(sg.operations())
    mult_str = (f'{na} × {nb} × {nc}' if (na, nb, nc) != (1, 1, 1)
                else '1 × 1 × 1  (no supercell expansion)')

    print()
    print('=' * 72)
    print(f'  Space group : {sg.xhm()}  (#{sg.number})')
    print(f'  Supercell   : {mult_str}')
    print(f'  Operators   : {len(ops)}')
    print('=' * 72)
    print()
    print('  Output structure factor formula:')
    print()
    print('    F_out(H,K,L) = Σ_i  phase_i(H,K,L)  ×  F_P1(h_i, k_i, l_i)')
    print()
    print('  where:')
    print('    (H\'_i, K\'_i, L\'_i) = R_direct_i^T · (H, K, L)  [transposed rotation]')
    if (na, nb, nc) != (1, 1, 1):
        print(f'    (h_i, k_i, l_i)    = ({na}·H\'_i, {nb}·K\'_i, {nc}·L\'_i)  [supercell lookup index]')
    else:
        print('    (h_i, k_i, l_i)    = (H\'_i, K\'_i, L\'_i)')
    print('    phase_i            = exp(2πi · H·t_i)  [H·t_i = H·ta + K·tb + L·tc]')
    print()

    mults = (na, nb, nc)

    for idx, op in enumerate(ops):
        rot = op.rot
        trn = op.tran

        xyz = op.triplet()
        hkl_prime = format_hkl_prime(rot)
        t_str = format_translation(trn)
        exponent, condition, factor = phase_factor_symbolic(trn)

        # Build P1 lookup string
        raw_parts = hkl_prime.split(', ')
        scaled = []
        for m, rp in zip(mults, raw_parts):
            if m == 1:
                scaled.append(rp)
            else:
                needs_parens = ('+' in rp or
                                ('-' in rp and not rp.startswith('-')))
                scaled.append(f'{m}·({rp})' if needs_parens else f'{m}·{rp}')
        lookup_str = '(' + ', '.join(scaled) + ')'

        print(f'  Operator {idx+1}: {xyz}')
        print(f'    Re-indexed (H\',K\',L\') = {hkl_prime}')
        print(f'    Translation t          = {t_str}')
        print(f'    Phase factor           = {factor}')
        print(f'    Symbolic exponent      = {exponent}')
        if condition:
            for line in condition.strip().split('\n'):
                print(f'                             {line}')
        if (na, nb, nc) != (1, 1, 1):
            print(f'    P1 lookup index        = {lookup_str}')
        print()

    print('  ' + '-' * 70)
    print()
    print('  NOTES:')
    print('  • R_direct_i^T is the TRANSPOSE of the direct-space rotation matrix.')
    print('    This differs from gemmi\'s apply_to_hkl(), which uses the')
    print('    reciprocal-space (non-transposed) form.')
    print('  • If the P1 MTZ lacks a required index, its Friedel mate is used:')
    print('    F_P1(-h,-k,-l) = conj(F_P1(h,k,l))')
    print('  • gemmi stores phase_shift = -2π H·t, so the code uses')
    print('    exp(-i · op.phase_shift([H,K,L])) to obtain exp(+2πi H·t).')
    print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    args = sys.argv[1:]
    sg_name = None
    na, nb, nc = 1, 1, 1

    for arg in args:
        if '=' in arg:
            key, val = arg.split('=', 1)
            key = key.lower().strip()
            if key in ('super_mult', 'mult', 'multipliers'):
                parts = val.split(',')
                if len(parts) != 3:
                    sys.exit("ERROR: super_mult must be na,nb,nc")
                na, nb, nc = int(parts[0]), int(parts[1]), int(parts[2])
        else:
            sg_test = gemmi.find_spacegroup_by_name(arg)
            if sg_test is not None:
                sg_name = arg
            else:
                print(f"WARNING: could not parse '{arg}' as a space group name, ignoring.")

    if sg_name is None:
        print(__doc__)
        sys.exit("ERROR: must specify a space group name")

    print_rules(sg_name, na, nb, nc)


if __name__ == '__main__':
    main()
