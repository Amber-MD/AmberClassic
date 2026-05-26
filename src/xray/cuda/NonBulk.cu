#include "NonBulk.h"

xray::NonBulk::NonBulk(int n_hkl, const int* hkl, complex_double* f_non_bulk, const double* mSS4, int n_atom, int nb,
                       const double* b_factor, const double* occupancy, int n_scatter_types, const int* scatter_type_index,
                       const double* atomic_scatter_factor) {

  this->m_n_hkl = n_hkl;
  this->m_hkl = hkl;
  this->m_f_non_bulk = f_non_bulk;
  this->m_mSS4 = mSS4;

  this->m_n_atom = n_atom;
  this->m_nb = nb;
  this->m_b_factor = b_factor;
  this->m_occupancy = occupancy;
  this->m_n_scatter_types = n_scatter_types;
  this->m_scatter_type_index = scatter_type_index;
  this->m_atomic_scatter_factor = atomic_scatter_factor;

}
