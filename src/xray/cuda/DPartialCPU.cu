#include <cassert>

#include "DPartialCPU.h"

void xray::DPartialCPU::calc_d_target_d_frac(
  int n_atom, int na, int nb, int nc, int sgn,
  const double* frac,
  int n_hkl,
  const double* f_scale,
  const double* d_target_d_abs_f_calc,
  double* d_target_d_frac) {

  assert(n_atom == m_n_atom);
  assert(n_hkl == m_n_hkl);

  std::fill(d_target_d_frac, d_target_d_frac + n_atom*3, 0.0);

  for (int i = 0; i < n_atom; ++i) {

    for (int i_hkl = 0; i_hkl < m_n_hkl; ++i_hkl) {
      auto f_abs = m_abs_f_calc[i_hkl];
      if (f_abs < 1e-3) {
        continue;
      }
      double phase = -(
        m_hkl[i_hkl * 3 + 0] * frac[i * 3 + 0] +
        m_hkl[i_hkl * 3 + 1] * frac[i * 3 + 1] +
        m_hkl[i_hkl * 3 + 2] * frac[i * 3 + 2]
      ) * 2 * M_PI;


      double scatter_factor = m_atomic_scatter_factor[(m_atom_scatter_type[i] - 1) * m_n_hkl + i_hkl];

      double f = scatter_factor * exp(m_mss4[i_hkl] * m_atom_b_factor[i]);

      std::complex<double> c_phase{cos(phase), sin(phase)};
      double tmp = 2 * M_PI * f * f_scale[i_hkl] *
        std::imag(c_phase * m_f_calc[i_hkl]) * d_target_d_abs_f_calc[i_hkl] / f_abs
        * m_atom_occupancy[i];

      d_target_d_frac[i * 3 + 0] += m_hkl[i_hkl * 3 + 0] * tmp;
      d_target_d_frac[i * 3 + 1] += m_hkl[i_hkl * 3 + 1] * tmp;
      d_target_d_frac[i * 3 + 2] += m_hkl[i_hkl * 3 + 2] * tmp;
    }
  }
}
