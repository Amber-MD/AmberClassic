#ifndef PDDF_H
#define PDDF_H

double smearing_e (const coordinate &atom,
                   const dx_type &dx,
                   int index,
                   size_t split,
                   const std::vector<coeff_f> &F_table);

void distributing_e (const coordinate &atom,
                     dx_type &dx,
                     double smear_cut,
                     double anom_f,
                     const std::vector<coeff_f> &F_table);

double calc_bulk_rho_e (const dx_type &dx,
                        const coordinate &max_bulk,
                        const coordinate &min_bulk);

std::vector<double> calc_pddf (const dx_type &dx1,
                               const std::vector<list_cutoff> &list1,
                               const dx_type &dx2,
                               const std::vector<list_cutoff> &list2,
                               double cutoff,
                               double dr);

#endif
