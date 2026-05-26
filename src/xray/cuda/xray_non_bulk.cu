#include "xray_non_bulk.h"
#include "NonBulkCPU.h"
#include "NonBulkGPU.h"
#include <cassert>
#include <memory>
#include <chrono>
#include <iostream>
#include <iomanip>

namespace {

std::unique_ptr<xray::NonBulk>& non_bulk_instance(){
  static std::unique_ptr<xray::NonBulk> non_bulk;
  return non_bulk;
}

}

extern "C" void pmemd_xray_non_bulk_init_gpu(
  int n_hkl,
  const int* hkl,
  complex_double* f_non_bulk,
  const double* mSS4,
  int n_atom,
  int nb,
  const double* b_factor,
  const double* occupancy,
  int n_scatter_types,
  const int* scatter_type_index,
  const double* atomic_scatter_factor) {
  non_bulk_instance().reset(
    new xray::NonBulkGPU<xray::NonBulkKernelVersion::StraightForward, xray::KernelPrecision::CUDA_PRECISION>(
      n_hkl, hkl, f_non_bulk, mSS4, n_atom, nb, b_factor, occupancy, n_scatter_types, scatter_type_index, atomic_scatter_factor
   ));
}

extern "C" void pmemd_xray_non_bulk_calc_f_non_bulk_gpu(
  int n_atoms,
  const double* frac_xyz
) {
  assert(non_bulk_instance());
  non_bulk_instance()->calc_f_non_bulk(n_atoms, frac_xyz);
}


extern "C" void pmemd_xray_non_bulk_finalize_gpu() {
  assert(non_bulk_instance());
  non_bulk_instance().reset();
}
