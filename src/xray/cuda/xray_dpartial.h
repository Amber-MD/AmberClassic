#ifndef AMBER_XRAY_DPARTIAL_H
#define AMBER_XRAY_DPARTIAL_H

#include "xray_common.h"

extern "C" {

void pmemd_xray_dpartial_init_gpu(
  int n_hkl,
  const int* hkl,
  const double* mss4,
  complex_double* f_calc,
  const double* abs_f_calc,
  int n_atom, int na, int nb, int nc, int sgn,
  const double* atom_b_factor,
  const double* atom_occupancy,
  const int* atom_scatter_type,
  int n_scatter_types,
  const double* atomic_scatter_factor
);

void pmemd_xray_dpartial_calc_d_target_d_frac(
  int n_atom, int na, int nb, int nc, int sgn,
  const double* frac,
  int n_hkl,
  const double* f_scale,
  const double* d_target_d_abs_f_calc,
  double* d_target_d_frac /* result variable */
);

void pmemd_xray_dpartial_finalize_gpu();
}

#endif //AMBER_XRAY_DPARTIAL_H
