#ifndef VOL_DECOMP_H
#define VOL_DECOMP_H

std::vector<list_cutoff> dx_list_cutoff (const dx_type &dx,
                                         const std::vector<coordinate> &pdb_coord,
                                         double cutoff);

bool isEdge_point (const dx_type &dx,
                   const std::vector<bool> &bin_grid,
                   size_t x, size_t y, size_t z);

void empty_neighbor (std::vector<bool> &bin_grid,
                     const dx_type &dx,
                     size_t x, size_t y, size_t z, double cutoff);

std::vector<list_cutoff> excluded_volume_list (const dx_type &dx,
                                               const std::vector<coordinate> &pdb_coord);

void gen_exV_hyd_list (const std::vector<coordinate> &pdb_coord,
                       const dx_type &dx,
                       std::vector<list_cutoff> &list_exclV,
                       std::vector<list_cutoff> &list_hyd);

#endif
