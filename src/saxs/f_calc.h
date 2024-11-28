#ifndef F_CALC_H
#define F_CALC_H

#include "saxsDS.h"

std::complex<double> form_factor (const std::vector<coordinate> &model,
                                  double q,
                                  const coordinate &Leb,
                                  double anom_f,
                                  bool flex,
                                  bool expli,
                                  const std::vector<coeff_f> &F_table);

std::complex<double> grid_factor (const dx_type &dx,
                                  int value_type,
                                  double atomic_factor,
                                  const std::vector<list_cutoff> &list,
                                  const coordinate &q_vector);

dx_type ex_elec (const std::vector<dx_type> &dx,
                 double anom_f,
                 bool expli,
                 const std::vector<coeff_f> &F_table);

#endif
