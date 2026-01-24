#ifndef ATM_F_H
#define ATM_F_H

#include <vector>
#include "saxsDS.h"

std::vector<coeff_f> assign_F_table ();

double atom_fact (const coeff_f &atom,
                  double q);

double count_e (const std::string &type,
                double anom_f);

double atom_rho (const std::string &type,
                 double dist_sqr,
                 const std::vector<coeff_f> &F_table);

double f_atm (const std::string &type_orig,
              double q_orig,
              size_t nHyd,
              bool expli,
              bool grid,
              const std::vector<coeff_f> &F_table);

#endif
