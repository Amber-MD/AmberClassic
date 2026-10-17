#ifndef AMBER_DPARTIAL_H
#define AMBER_DPARTIAL_H

#include <complex>

namespace xray {
  class DPartial {
  public:
    DPartial(
      int n_hkl,
      const int* hkl,
      const double* mss4,
      std::complex<double>* f_calc,
      const double* abs_f_calc,
      int n_atom, int na, int nb, int nc, int sgn,
      const double* atom_b_factor,
      const double* atom_occupancy,
      const int* atom_scatter_type,
      int n_scatter_types,
      const double* atomic_scatter_factor
    );

    virtual ~DPartial() = default;

    virtual void calc_d_target_d_frac(
      int n_atom, int na, int nb, int nc, int sgn,
      const double* frac,
      int n_hkl,
      const double* f_scale,
      const double* d_target_d_abs_f_calc,
      double* d_target_d_frac /* result variable */
    ) = 0;

  protected:
    int m_n_hkl;
    const int* m_hkl;
    const double* m_mss4;
    const std::complex<double>* m_f_calc;
    const double* m_abs_f_calc;
    int m_n_atom = 0;
    const double* m_atom_b_factor;
    const double* m_atom_occupancy;
    const int* m_atom_scatter_type;
    int m_n_scatter_types;
    const double* m_atomic_scatter_factor;
  };
}

#endif //AMBER_DPARTIAL_H
