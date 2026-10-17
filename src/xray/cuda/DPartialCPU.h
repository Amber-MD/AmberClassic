#ifndef AMBER_DPARTIAL_CPU_H
#define AMBER_DPARTIAL_CPU_H

#include "DPartial.h"

namespace xray {
  class DPartialCPU : public DPartial {
  public:
    using DPartial::DPartial;

    void calc_d_target_d_frac(
      int n_atom, int na, int nb, int nc, int sgn,
      const double* frac,
      int n_hkl,
      const double* f_scale,
      const double* d_target_d_abs_f_calc,
      double* d_target_d_frac /* result variable */
    ) override;

  };
}

#endif //AMBER_DPARTIAL_CPU_H
