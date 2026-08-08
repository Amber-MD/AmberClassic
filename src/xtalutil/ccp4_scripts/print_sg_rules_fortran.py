#!/usr/bin/env ccp4-python
"""
print_sg_rules.py
=================
Print a Fortran subroutine that folds a P1 supercell structure-factor array
into the ASU of a target space group, implementing exactly the combination
rules used by supercell_collapse.

The generated subroutine has the signature:

  SUBROUTINE COLLAPSE_<SGNAME>(F_P1, NA, NB, NC, IH, IK, IL, F_OUT)
    ! F_P1(ih,ik,il) : complex P1 supercell array (caller must size it)
    ! NA,NB,NC       : supercell multipliers
    ! IH,IK,IL       : primitive-cell ASU reflection to compute
    ! F_OUT          : output complex structure factor

For each symmetry operator the subroutine:
  1. Computes the re-indexed primitive-cell index (H',K',L') from
     R_direct^T * (IH,IK,IL)
  2. Scales by NA,NB,NC to get the supercell index
  3. Multiplies the looked-up F_P1 by the phase factor exp(2*pi*i * H.t)
     - for t in {0,1/2}:   an IF on the parity sum gives +/- REAL sign
     - for t in {0,1/4,...}: a SELECT CASE on mod(sum,4) gives +/-1/+/-i
     - general t:           EXP(CMPLX(0.0, phase)) computed numerically
  4. Accumulates into F_OUT

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
TWO_PI = '6.2831853071795864769'  # 2*pi as Fortran literal


# ---------------------------------------------------------------------------
# Helpers: index expressions
# ---------------------------------------------------------------------------

def hkl_prime_coeffs(rot):
    """
    Return the integer coefficients of R_direct^T * (IH,IK,IL) as a list of
    three lists: coeffs[out_idx] = [c_H, c_K, c_L].
    """
    result = []
    for col in range(3):
        result.append([rot[row][col] // DEN for row in range(3)])
    return result


def fortran_index_expr(coeffs, labels=('IH', 'IK', 'IL')):
    """
    Return a Fortran expression string for a linear combination of IH/IK/IL.
    E.g. coeffs=[-1,0,1] -> '-IH+IL'
    """
    terms = []
    for c, lbl in zip(coeffs, labels):
        if c == 0:
            continue
        elif c == 1:
            terms.append(lbl)
        elif c == -1:
            terms.append(f'-{lbl}')
        else:
            terms.append(f'{c}*{lbl}')
    if not terms:
        return '0'
    s = terms[0]
    for t in terms[1:]:
        s += (t if t.startswith('-') else '+' + t)
    return s


def fortran_lookup_expr(coeffs, mults, labels=('IH', 'IK', 'IL')):
    """
    Return the Fortran index expression for the supercell lookup:
    m * (R^T * index).  E.g. m=2, coeffs=[-1,0,1] -> '-2*IH+2*IL'
    """
    terms = []
    for c, m, lbl in zip(coeffs, mults, labels):
        net = c * m
        if net == 0:
            continue
        elif net == 1:
            terms.append(lbl)
        elif net == -1:
            terms.append(f'-{lbl}')
        else:
            terms.append(f'{net}*{lbl}')
    if not terms:
        return '0'
    s = terms[0]
    for t in terms[1:]:
        s += (t if t.startswith('-') else '+' + t)
    return s


# ---------------------------------------------------------------------------
# Helpers: phase factor code generation
# ---------------------------------------------------------------------------

def phase_fortran(trn, ih_expr, ik_expr, il_expr, op_idx):
    """
    Return a list of Fortran lines that compute PHASE_i (a COMPLEX(8) scalar)
    equal to exp(2*pi*i * (ta*IH + tb*IK + tc*IL)).

    ih_expr / ik_expr / il_expr are the already-computed re-indexed integer
    expressions (strings) to use in parity tests.

    op_idx is used to name any temporary integer variable (ISUM_i).
    """
    fracs = [Fraction(v, DEN) for v in trn]
    labels_expr = (ih_expr, ik_expr, il_expr)
    phase_var = f'PHASE_{op_idx}'

    # --- Case 1: all t_j in {0, 1/2}  -->  pure real ±1 ---
    if all(f.denominator in (1, 2) for f in fracs):
        int_coeffs = [int(f * 2) % 2 for f in fracs]   # 0 or 1
        parity_exprs = [expr for c, expr in zip(int_coeffs, labels_expr) if c != 0]

        if not parity_exprs:
            return [f'      {phase_var} = CMPLX(1.0D0, 0.0D0, KIND=8)']

        # Build parity sum expression
        isum = f'ISUM_{op_idx}'
        sum_expr = '+'.join(f'({e})' for e in parity_exprs)
        lines = [
            f'      {isum} = MOD({sum_expr}, 2)',
            f'      IF ({isum} .EQ. 0) THEN',
            f'        {phase_var} = CMPLX( 1.0D0, 0.0D0, KIND=8)',
            f'      ELSE',
            f'        {phase_var} = CMPLX(-1.0D0, 0.0D0, KIND=8)',
            f'      END IF',
        ]
        return lines

    # --- Case 2: all t_j in {0, 1/4, 1/2, 3/4}  -->  ±1 or ±i ---
    if all(f.denominator in (1, 2, 4) for f in fracs):
        qcoeffs = [int(f * 4) % 4 for f in fracs]
        q_exprs = []
        for c, expr in zip(qcoeffs, labels_expr):
            if c == 0:
                continue
            elif c == 1:
                q_exprs.append(f'({expr})')
            else:
                q_exprs.append(f'{c}*({expr})')
        if not q_exprs:
            return [f'      {phase_var} = CMPLX(1.0D0, 0.0D0, KIND=8)']

        isum = f'ISUM_{op_idx}'
        sum_expr = '+'.join(q_exprs)
        quarter_vals = [
            'CMPLX( 1.0D0,  0.0D0, KIND=8)',
            'CMPLX( 0.0D0,  1.0D0, KIND=8)',
            'CMPLX(-1.0D0,  0.0D0, KIND=8)',
            'CMPLX( 0.0D0, -1.0D0, KIND=8)',
        ]
        lines = [f'      {isum} = MOD({sum_expr}, 4)']
        # MOD can be negative in Fortran for negative arguments
        lines.append(f'      IF ({isum} .LT. 0) {isum} = {isum} + 4')
        lines.append(f'      SELECT CASE ({isum})')
        for n, val in enumerate(quarter_vals):
            lines.append(f'        CASE ({n})')
            lines.append(f'          {phase_var} = {val}')
        lines.append(f'      END SELECT')
        return lines

    # --- Case 3: general translation (hexagonal, rhombohedral, ...) ---
    # Compute phase angle = 2*pi * (ta*IH + tb*IK + tc*IL) numerically.
    phase_terms = []
    for f, expr in zip(fracs, labels_expr):
        if f == 0:
            continue
        if f.denominator == 1:
            phase_terms.append(f'{f.numerator}*DBLE({expr})')
        else:
            phase_terms.append(f'({f.numerator}.0D0/{f.denominator}.0D0)*DBLE({expr})')
    if not phase_terms:
        return [f'      {phase_var} = CMPLX(1.0D0, 0.0D0, KIND=8)']
    angle_expr = '+'.join(phase_terms)
    angle_var = f'ANGLE_{op_idx}'
    lines = [
        f'      {angle_var} = {TWO_PI} * ({angle_expr})',
        f'      {phase_var} = CMPLX(DCOS({angle_var}), DSIN({angle_var}), KIND=8)',
    ]
    return lines


# ---------------------------------------------------------------------------
# Main Fortran subroutine generator
# ---------------------------------------------------------------------------

def print_fortran(sg_name, na=1, nb=1, nc=1):
    sg = gemmi.find_spacegroup_by_name(sg_name)
    if sg is None:
        sys.exit(f"ERROR: unknown space group '{sg_name}'")

    ops = list(sg.operations())

    # Make a safe Fortran identifier from the space group name
    sg_tag = sg.xhm().replace(' ', '').replace('/', 'S').replace('-', 'M')
    sub_name = f'COLLAPSE_{sg_tag}'

    # Collect all temporary variable names needed
    needs_angle = []
    needs_isum  = []
    for idx, op in enumerate(ops):
        fracs = [Fraction(v, DEN) for v in op.tran]
        if not all(f == 0 for f in fracs):
            if all(f.denominator in (1, 2, 4) for f in fracs):
                needs_isum.append(idx + 1)
            else:
                needs_angle.append(idx + 1)
                needs_isum.append(idx + 1)  # not actually needed but harmless

    # -----------------------------------------------------------------------
    # Header comment
    # -----------------------------------------------------------------------
    L = []
    def emit(s=''):
        L.append(s)

    emit(f'! {"="*68}')
    emit(f'! Subroutine generated by print_sg_rules.py')
    emit(f'! Space group : {sg.xhm()}  (#{sg.number})')
    emit(f'! Supercell   : {na} x {nb} x {nc}')
    emit(f'! Operators   : {len(ops)}')
    emit(f'!')
    emit(f'! Folds a P1 supercell structure-factor array into a single')
    emit(f'! primitive-cell ASU reflection using the symmetry of {sg.xhm()}.')
    emit(f'!')
    emit(f'! Formula:')
    emit(f'!   F_OUT = SUM_i  phase_i(IH,IK,IL)  *  F_P1(ih_i, ik_i, il_i)')
    emit(f'! where (ih_i,ik_i,il_i) = (na*H\'_i, nb*K\'_i, nc*L\'_i)')
    emit(f'! and   (H\'_i,K\'_i,L\'_i) = R_direct_i^T * (IH,IK,IL)')
    emit(f'! {"="*68}')

    emit(f'      SUBROUTINE {sub_name}(F_P1, NA, NB, NC, IH, IK, IL, F_OUT)')
    emit(f'      IMPLICIT NONE')
    emit()
    emit(f'      ! Supercell multipliers (passed in; must match array dimensions)')
    emit(f'      INTEGER, INTENT(IN) :: NA, NB, NC')
    emit(f'      ! Primitive-cell ASU reflection to compute')
    emit(f'      INTEGER, INTENT(IN) :: IH, IK, IL')
    emit(f'      ! P1 supercell structure factor array F_P1(ih, ik, il)')
    emit(f'      ! Caller must dimension this to cover all required indices.')
    emit(f'      COMPLEX(8), INTENT(IN) :: F_P1(-NA*10:NA*10, -NB*10:NB*10, -NC*10:NC*10)')
    emit(f'      ! Output: accumulated structure factor for this ASU reflection')
    emit(f'      COMPLEX(8), INTENT(OUT) :: F_OUT')
    emit()

    # Declare temporaries
    emit(f'      ! --- local temporaries ---')
    for i, op in enumerate(ops):
        emit(f'      COMPLEX(8) :: PHASE_{i+1}')
    emit()
    # Integer temporaries for parity sums
    isum_needed = set()
    angle_needed = set()
    for i, op in enumerate(ops):
        fracs = [Fraction(v, DEN) for v in op.tran]
        if all(f == 0 for f in fracs):
            continue
        if all(f.denominator in (1, 2, 4) for f in fracs):
            isum_needed.add(i + 1)
        else:
            angle_needed.add(i + 1)
    if isum_needed:
        vars_ = ', '.join(f'ISUM_{i}' for i in sorted(isum_needed))
        emit(f'      INTEGER :: {vars_}')
    if angle_needed:
        vars_ = ', '.join(f'ANGLE_{i}' for i in sorted(angle_needed))
        emit(f'      REAL(8) :: {vars_}')
    emit()
    emit(f'      F_OUT = CMPLX(0.0D0, 0.0D0, KIND=8)')
    emit()

    # -----------------------------------------------------------------------
    # One block per operator
    # -----------------------------------------------------------------------
    for i, op in enumerate(ops):
        rot = op.rot
        trn = op.tran
        xyz = op.triplet()
        coeffs = hkl_prime_coeffs(rot)   # coeffs[0..2] = [cH,cK,cL] for H',K',L'

        ih_expr = fortran_index_expr(coeffs[0])
        ik_expr = fortran_index_expr(coeffs[1])
        il_expr = fortran_index_expr(coeffs[2])

        sc_ih = fortran_lookup_expr(coeffs[0], [na, nb, nc][0:1] + [0, 0],
                                    labels=('IH', 'IK', 'IL'))
        # Redo properly with per-axis multiplier
        sc_ih = fortran_lookup_expr(coeffs[0], (na, nb, nc))
        sc_ik = fortran_lookup_expr(coeffs[1], (na, nb, nc))
        sc_il = fortran_lookup_expr(coeffs[2], (na, nb, nc))

        emit(f'      ! --- Operator {i+1}: {xyz} ---')
        emit(f'      !     H\',K\',L\' = {ih_expr}, {ik_expr}, {il_expr}')
        emit(f'      !     supercell lookup: F_P1({sc_ih}, {sc_ik}, {sc_il})')

        phase_lines = phase_fortran(trn, ih_expr, ik_expr, il_expr, i + 1)
        for line in phase_lines:
            emit(line)

        emit(f'      F_OUT = F_OUT + PHASE_{i+1} * F_P1({sc_ih}, {sc_ik}, {sc_il})')
        emit()

    emit(f'      END SUBROUTINE {sub_name}')

    print('\n'.join(L))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


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

    print_fortran(sg_name, na, nb, nc)


if __name__ == '__main__':
    main()
