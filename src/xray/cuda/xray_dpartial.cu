#include "xray_dpartial.h"
#include "DPartialCPU.h"
#include "DPartialGPU.h"
#include <memory>
#include <cassert>
#include <chrono>
#include <iomanip>

namespace {

std::unique_ptr<xray::DPartial>& dpartial_instance(){
  static std::unique_ptr<xray::DPartial> dpartial;
  return dpartial;
}

}

extern "C"
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
) {
  assert(!dpartial_instance());
  dpartial_instance().reset(
    new xray::DPartialGPU<xray::KernelPrecision::CUDA_PRECISION>(
      n_hkl,
      hkl,
      mss4,
      reinterpret_cast<std::complex<double>*>(f_calc),
      abs_f_calc,
      n_atom, na, nb, nc, sgn,
      atom_b_factor,
      atom_occupancy,
      atom_scatter_type,
      n_scatter_types,
      atomic_scatter_factor
    ));
}

extern "C"
void pmemd_xray_dpartial_calc_d_target_d_frac(
  int n_atom, int na, int nb, int nc, int sgn,
  const double* frac,
  int n_hkl,
  const double* f_scale,
  const double* d_target_d_abs_f_calc,
  double* d_target_d_frac /* result variable */
) {
  assert(dpartial_instance());
  dpartial_instance()->calc_d_target_d_frac(
    n_atom, na, nb, nc, sgn, frac, n_hkl, f_scale, d_target_d_abs_f_calc, d_target_d_frac
  );
}

extern "C"
void pmemd_xray_dpartial_finalize_gpu() {
  assert(dpartial_instance());
  dpartial_instance().reset();
}
