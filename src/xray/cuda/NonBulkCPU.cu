#include "NonBulkCPU.h"
#include <cassert>
#include <cstdio>


void xray::NonBulkCPU::calc_f_non_bulk(int n_atom, int na, int nb, 
     int nc, int sgn, const double* frac_xyz) {
  assert(n_atom == m_n_atom);
  for (int i_hkl = 0; i_hkl < m_n_hkl; ++i_hkl) {
    complex_double term{0, 0};
    for (int j_atom = 0; j_atom < n_atom; j_atom++) {
      double f = std::exp(m_mSS4[i_hkl] * m_b_factor[j_atom]) *
                 m_atomic_scatter_factor[(m_scatter_type_index[j_atom] - 1) * m_n_hkl + i_hkl];

      double angle = 2 * M_PI * (
        frac_xyz[j_atom * 3 + 0] * m_hkl[i_hkl * 3 + 0] +
        frac_xyz[j_atom * 3 + 1] * m_hkl[i_hkl * 3 + 1] +
        frac_xyz[j_atom * 3 + 2] * m_hkl[i_hkl * 3 + 2]
      );

      term += complex_double{f * std::cos(angle), f * std::sin(angle)} * m_occupancy[j_atom];
    }
    m_f_non_bulk[i_hkl] = term;
  }
}
