#ifndef PDB_H
#define PDB_H

#include "saxsDS.h"

void read_pdb (const std::string &pdb_file,
               std::vector<coordinate> &pdb_coord);

void read_pdb_weight (const std::string &pdb_file,
                      std::vector< std::vector <coordinate> > &snapshot,
                      std::vector<unsigned> &weight);

void mergeH (std::vector<coordinate> &pdb_coord);

void maxmin_coord_pdb (const std::vector<coordinate> &pdb_coord,
                       coordinate &max,
                       coordinate &min);

void maxmin_coord_traj (const std::vector< std::vector<coordinate> > &snapshot,
                        coordinate &max,
                        coordinate &min);

bool check_solute (const coordinate &coord);

std::vector<coordinate> solute_coord (const std::vector<coordinate> &pdb_coord);

std::vector<coordinate> strip_buffer_atom (const std::vector<coordinate> &model,
                                           double dcutoff,
                                           const coordinate &max_solu,
                                           const coordinate &min_solu);

bool isSmear (const coordinate &atom,
              const coordinate &max_box,
              const coordinate &min_box);
#endif
