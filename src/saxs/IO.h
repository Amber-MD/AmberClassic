#ifndef IO_H
#define IO_H

void saxs_rism_usage ();
void saxs_md_usage ();
void pddf_rism_usage ();
void pddf_md_usage ();

void saxs_rism_output (const std::vector< std::vector<double> > &intensity,
                       const std::vector<dx_type> &dx,
                       const std::vector<double> &q,
                       std::string &exper,
                       std::string &dx_dir, std::string &pdb_file,
                       double conc_salt, double conc_wat, double cutoff,
                       double anom_f, bool isFlex, bool isOff_cutoff, bool isExpli,
                       bool isExclV, bool isTight, bool isDecomp, bool isPhase,
                       int index_coion,
                       const std::vector< std::vector<DoubPair> > &phase,
                       const std::vector< std::vector<DoubPair> > &dphase,
                       std::ofstream &OUTPUT);

void saxs_rism_output_intensity (const std::vector< std::vector<double> > &intensity,
                                 const std::vector<dx_type> &dx,
                                 const std::vector<double> &q,
                                 bool isExclV, bool isDecomp, int index_coion,
                                 std::ofstream &OUTPUT);

void saxs_rism_output_phase (const std::vector< std::vector<double> > &intensity,
                             const std::vector<dx_type> &dx,
                             const std::vector<double> &q,
                             bool isExclV, bool isDecomp, int index_coion,
                             const std::vector< std::vector<DoubPair> > &phase,
                             const std::vector< std::vector<DoubPair> > &dphase,
                             std::ofstream &OUTPUT);

void saxs_md_output (const std::vector<double> &intensity,
                     double dq,
                     std::string pdb_solu, std::string pdb_solv,
                     double dcutoff, double bulk_corr, double bulk_cutoff,
                     double anom_f, bool isExpli, bool isTight,
                     std::ofstream &OUTPUT);

void pddf_output (const std::vector<double> &pddf,
                  size_t type,
                  double dr,
                  std::ofstream &OUTPUT);

void read_RDF (const std::string &filename,
               double &dr,
               std::vector<double> &rdf);

void read_exp (const std::string &exper,
               std::vector<double> &q);

#endif
