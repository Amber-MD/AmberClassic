#ifndef SAXS_MD_H
#define SAXS_MD_H

#include <string>
#include <vector>
#include <complex>
#include "saxsDS.h"


std::complex<double> mean_complex (const std::vector< std::complex<double> > &v);

double saxs_md_D11 (const std::vector< std::vector<coordinate> > &solu,
                    const std::vector<unsigned> &weight_solu,
                    const std::vector< std::vector<coordinate> > &solv,
                    const std::vector<unsigned> &weight_solv,
                    double q, double bulk_corr,
                    const coordinate &max_buffer_cut,
                    const coordinate &min_buffer_cut,
                    const coordinate &Leb,
                    double anom_f, bool isExpli,
                    const std::vector<coeff_f> &F_table);

double saxs_md_I (const std::vector< std::vector<coordinate> > &solu,
                  const std::vector<unsigned> &weight_solu,
                  const std::vector< std::vector<coordinate> > &solv,
                  const std::vector<unsigned> &weight_solv,
                  double q, double bulk_corr,
                  const coordinate &max_buffer_cut,
                  const coordinate &min_buffer_cut,
                  size_t rule, double anom_f, bool isExpli,
                  const std::vector<coeff_f> &F_table);

std::vector<double> saxs_md_calc (std::vector< std::vector<coordinate> > &solu_box,
                                  std::vector< std::vector<coordinate> > &solv_box,
                                  const std::vector<unsigned> &weight_solu,
                                  const std::vector<unsigned> &weight_solv,
                                  const std::vector<double> &q,
                                  double bulk_corr, double bulk_cutoff,
                                  double anom_f, double qcut, double dq,
                                  double dcutoff, bool expli, bool tight,
                                  unsigned ncpus);

#endif
