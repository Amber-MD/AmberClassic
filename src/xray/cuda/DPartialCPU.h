#ifndef AMBER_DPARTIAL_CPU_H
#define AMBER_DPARTIAL_CPU_H

#include "DPartial.h"

namespace xray {
  class DPartialCPU : public DPartial {
  public:
    using DPartial::DPartial;

    void calc_d_target_d_frac(
      int n_atom,
      const double* frac,
      int n_hkl,
      const double* d_target_d_abs_f_calc,
      double* d_target_d_frac /* result variable */
    ) override;

  };
}

#endif //AMBER_DPARTIAL_CPU_H
