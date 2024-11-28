#ifndef SAXS_RISM_H
#define SAXS_RISM_H

void read_dx_dir (const std::string &dx_dir,
                  std::vector<dx_type> &dx,
                  double conc_wat, double conc_salt);

void generate_list (const std::vector<dx_type> &dx,
                    const std::vector<coordinate> &pdb_coord,
                    double cutoff,
                    bool isOff_cutoff, bool &isDecomp, bool isExclV,
                    std::vector< std::vector<list_cutoff> > &list_exclV,
                    std::vector< std::vector<list_cutoff> > &list_hyd,
                    std::vector<list_cutoff> &list_excl_coion,
                    bool isCoion, int &index_coion);

//bool index_sort (const list_cutoff &left,
//                 const list_cutoff &right);

bool dist_sort (const list_cutoff &left,
                const list_cutoff &right);

/*double increasing_center (const dx_type &dx,
                          double cutoff,
                          size_t Z,
                          double dr,
                          const std::vector<double> &rdf);*/

std::vector< std::complex<double> > saxs_rism_amplitude (const std::vector<dx_type> &dx,
                                                         bool isExcess,
                                                         const std::vector<coordinate> &pdb_coord,
                                                         const std::vector< std::vector<list_cutoff> > &list_hyd,
                                                         const std::vector< std::vector<list_cutoff> > &list_exclV,
                                                         const std::vector<list_cutoff> &list_excl_coion,
                                                         int index_coion, double q,
                                                         const coordinate &Leb,
                                                         bool isFlex, double anom_f, bool isExpli, bool isExclV,
                                                         const std::vector<coeff_f> &F_table);

std::vector<double> saxs_rism_I (const std::vector<dx_type> &dx,
                                 bool isExcess,
                                 const std::vector<coordinate> &pdb_coord,
                                 const std::vector< std::vector<list_cutoff> > &list_hyd,
                                 const std::vector< std::vector<list_cutoff> > &list_exclV,
                                 const std::vector<list_cutoff> &list_excl_coion,
                                 int index_coion, double q, size_t rule,
                                 bool isFlex, double anom_f, bool isExpli, bool isExclV,
                                 const std::vector<coeff_f> &F_table,
                                 std::vector<DoubPair> &dphase,
                                 std::vector<DoubPair> &averaged_phase);

std::vector< std::vector<double> > saxs_rism_calc (const std::vector<dx_type> &dx,
                                                   std::vector<coordinate> &pdb_coord,
                                                   const std::vector<double> &q,
                                                   double cutoff,
                                                   bool isTight, bool isFlex, double anom_f,
                                                   bool isExpli, bool isDecomp, bool isOff_cutoff,
                                                   bool isExclV, bool isCoion, int &index_coion, unsigned ncpus,
                                                   std::vector< std::vector<DoubPair> > &phase,
                                                   std::vector< std::vector<DoubPair> > &dphase);

#endif
