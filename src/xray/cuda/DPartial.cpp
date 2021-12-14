#include "xray/DPartial.h"


xray::DPartial::DPartial(int n_hkl, const int* hkl, const double* mss4, std::complex<double>* f_calc,
                         const double* abs_f_calc, int n_atom, const double* atom_b_factor,
                         const int* atom_scatter_type, int n_scatter_types, const double* atomic_scatter_factor) {
  m_n_hkl = n_hkl;
  m_hkl = hkl;
  m_mss4 = mss4;
  m_f_calc = f_calc;
  m_abs_f_calc = abs_f_calc;
  m_n_atom = n_atom;
  m_atom_b_factor = atom_b_factor;
  m_atom_scatter_type = atom_scatter_type;
  m_n_scatter_types = n_scatter_types;
  m_atomic_scatter_factor = atomic_scatter_factor;
}