#ifndef AMBER_XRAY_NON_BULK_H
#define AMBER_XRAY_NON_BULK_H

#include "xray_common.h"

extern "C" {

void pmemd_xray_non_bulk_init_gpu(
  int num_hkl,
  const int* hkl,
  complex_double* Fcalc,
  const double* mSS4,
  int n_atoms,
  int nb,
  const double* b_factor,
  const double* occupancy,
  int n_scatter_types,
  const int* scatter_type_index,
  const double* atomic_scatter_factor);

void pmemd_xray_non_bulk_calc_f_non_bulk_gpu(
  int n_atoms,
  const double* frac_xyz
);

void pmemd_xray_non_bulk_finalize_gpu();

}

#endif //AMBER_XRAY_NON_BULK_H
