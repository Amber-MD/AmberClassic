#ifndef DX_H
#define DX_H

void read_dx (std::string &dx_file,
              dx_type &dx);

void write_dx (const dx_type &dx,
               std::string &outfile);

bool check_dx (const std::vector<dx_type> &dx);

coordinate dx_1Dindex2coord (const dx_type &dx,
                             int index);

int dx_3Dindex_to_1Dindex (const dx_type &dx,
                           int x, int y, int z);

double grid2grid_dist (const dx_type &dx,
                       size_t index1, size_t index2);

void extract_dx (const dx_type &dx,
                 dx_type &dx_out,
                 coordinate &min,
                 coordinate &max);

dx_type scale_dx (const dx_type &dx,
                  double a, double b);

double integral_dx (const dx_type &dx);

dx_type gen_grid (const coordinate &max_coord,
                  const coordinate &min_coord,
                  const coordinate &res);

dx_type gen_coarse_dx (const dx_type &dx,
                       size_t merge_x, size_t merge_y, size_t merge_z);

void coarse_grid (const dx_type &fine_dx,
                  dx_type &coarse_dx,
                  size_t merge_x, size_t merge_y, size_t merge_z);
#endif
