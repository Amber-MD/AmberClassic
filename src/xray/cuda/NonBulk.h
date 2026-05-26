#ifndef AMBER_NON_BULK_H
#define AMBER_NON_BULK_H

#include "xray_common.h"

namespace xray {
  class NonBulk {
  public:
    NonBulk(int n_hkl,
            const int* hkl,
            complex_double* f_non_bulk,
            const double* mSS4,
            int n_atom,
            int nb,
            const double* b_factor,
            const double* occupancy,
            int n_scatter_types,
            const int* scatter_type_index,
            const double* atomic_scatter_factor);

    virtual ~NonBulk() = default;

    virtual void calc_f_non_bulk(
      int n_atoms,
      const double* frac_xyz
    ) = 0;

  protected:
    int m_n_hkl = 0;
    const int* m_hkl = nullptr;
    complex_double* m_f_non_bulk = nullptr;
    const double* m_mSS4 = nullptr;

    int m_n_atom = 0;
    int m_nb = 0;
    const double* m_b_factor = nullptr;
    const double* m_occupancy = nullptr;
    int m_n_scatter_types = 0;
    const int* m_scatter_type_index = nullptr;
    const double* m_atomic_scatter_factor = nullptr;
  };
}
#endif //AMBER_NON_BULK_H
