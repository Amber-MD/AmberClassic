#include <vector>
#include <cmath>
#include "saxsDS.h"
#include "dx.h"
#include "pdb.h"
#include "const.h"

std::vector<list_cutoff> dx_list_cutoff (const dx_type &dx,
                                         const std::vector<coordinate> &pdb_coord,
                                         double cutoff) {
    // Set up a binary grid to encode which grid point lies within cutoff from the solute
    // The second argument stores the distance to the nearest solute atom
    std::vector< std::pair<bool, double> > binary_grid;
    binary_grid.resize (dx.value.size(), std::make_pair(0, INF));

    #pragma omp parallel for shared (binary_grid, cutoff)
    // Turn the values of every grid point within cutoff from any atom to 1
    for (size_t i = 0; i < pdb_coord.size(); i++) {
        double rcut = cutoff + pdb_coord[i].r;

        // Only consider local grids
        size_t minx = 0;
        double tmp = (pdb_coord[i].x - rcut - dx.origin[0]) / dx.delta[0];
        if (tmp > 0)
            minx = (size_t) ceil(tmp);

        size_t maxx = dx.ngrid[0] - 1;
        tmp = (pdb_coord[i].x + rcut - dx.origin[0]) / dx.delta[0];
        if (tmp < maxx)
            maxx = (size_t) floor(tmp);

        size_t miny = 0;
        tmp = (pdb_coord[i].y - rcut - dx.origin[1]) / dx.delta[1];
        if (tmp > 0)
            miny = (size_t) ceil(tmp);

        size_t maxy = dx.ngrid[1] - 1;
        tmp = (pdb_coord[i].y + rcut - dx.origin[1]) / dx.delta[1];
        if (tmp < maxy)
            maxy = (size_t) floor(tmp);

        size_t minz = 0;
        tmp = (pdb_coord[i].z - rcut - dx.origin[2]) / dx.delta[2];
        if (tmp > 0)
            minz = (size_t) ceil(tmp);

        size_t maxz = dx.ngrid[2] - 1;
        tmp = (pdb_coord[i].z + rcut - dx.origin[2]) / dx.delta[2];
        if (tmp < maxz)
            maxz = (size_t) floor(tmp);

        for (size_t x = minx; x <= maxx; x++)
            for (size_t y = miny; y <= maxy; y++)
                for (size_t z = minz; z <= maxz; z++) {
                    coordinate grid;
                    grid.x = x*dx.delta[0] + dx.origin[0];
                    grid.y = y*dx.delta[1] + dx.origin[1];
                    grid.z = z*dx.delta[2] + dx.origin[2];

                    double deltax = pdb_coord[i].x - grid.x;
                    double deltay = pdb_coord[i].y - grid.y;
                    double deltaz = pdb_coord[i].z - grid.z;
                    double dist = sqrt (deltax*deltax + deltay*deltay + deltaz*deltaz);
                    if (dist <= rcut) {
                        size_t index = dx_3Dindex_to_1Dindex (dx, x, y, z);
                        #pragma omp critical
                        {
                            binary_grid[index].first = 1;
                            if (dist < binary_grid[index].second)
                                binary_grid[index].second = dist;
                        }
                    }
                }
    }

    // Now, base on the binary grid, push the indexes and distances to the list vector
    std::vector<list_cutoff> result;
    for (size_t i = 0; i < binary_grid.size(); i++)
        if (binary_grid[i].first) {
            list_cutoff append_this;
            append_this.index = i;
            append_this.dist  = binary_grid[i].second;
            result.push_back (append_this);
        }

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////
// Empty (unpaint) grid points lying within cutoff from a particular points (x, y, z)
//////////////////////////////////////////////////////////////////////////////////////
void empty_neighbor (std::vector<bool> &bin_grid,
                     const dx_type &dx,
                     size_t x, size_t y, size_t z,
                     double cutoff) {

    size_t xmin = std::max ((int) ceil  (x - cutoff/dx.delta[0]), 0);
    size_t xmax = std::min ((int) floor (x + cutoff/dx.delta[0]), (int) dx.ngrid[0] - 1);
    size_t ymin = std::max ((int) ceil  (y - cutoff/dx.delta[1]), 0);
    size_t ymax = std::min ((int) floor (y + cutoff/dx.delta[1]), (int) dx.ngrid[1] - 1);
    size_t zmin = std::max ((int) ceil  (z - cutoff/dx.delta[2]), 0);
    size_t zmax = std::min ((int) floor (z + cutoff/dx.delta[2]), (int) dx.ngrid[2] - 1);

    for (size_t i = xmin; i <= xmax; i++)
        for (size_t j = ymin; j <= ymax; j++)
            for (size_t k = zmin; k <= zmax; k++) {
                size_t index = dx_3Dindex_to_1Dindex (dx, i, j, k);
                if (bin_grid[index]) {
                    double deltax = int (x - i) * dx.delta[0];
                    double deltay = int (y - j) * dx.delta[1];
                    double deltaz = int (z - k) * dx.delta[2];
                    double dist = deltax*deltax + deltay*deltay + deltaz*deltaz;
                    if (dist <= cutoff*cutoff)
                        bin_grid[index] = 0;
                }
            }
}

///////////////////////////////////////////////////////
bool isEdge_point (const dx_type &dx,
                   const std::vector<bool> &bin_grid,
                   size_t x, size_t y, size_t z) {

    size_t count = 0;
    if (x > 0)               count += bin_grid [dx_3Dindex_to_1Dindex (dx, x-1, y  , z  )];
    if (x < dx.ngrid[0] - 1) count += bin_grid [dx_3Dindex_to_1Dindex (dx, x+1, y  , z  )];
    if (y > 0)               count += bin_grid [dx_3Dindex_to_1Dindex (dx, x  , y-1, z  )];
    if (y < dx.ngrid[1] - 1) count += bin_grid [dx_3Dindex_to_1Dindex (dx, x  , y+1, z  )];
    if (z > 0)               count += bin_grid [dx_3Dindex_to_1Dindex (dx, x  , y  , z-1)];
    if (z < dx.ngrid[2] - 1) count += bin_grid [dx_3Dindex_to_1Dindex (dx, x  , y  , z+1)];

    if (count == 0) return 0;
    else return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Generate a list of indexes of grid points lying within the excluded volume of solute
// The idea is:
//      + Generate a volume within the solute solvent accessible surface, assuming r_wat = 1.4A
//      + Start unpainting every points lying within r_wat from any grid points that are not in the above volume
//      + The remainer will be the excluded volume
/////////////////////////////////////////////////////////////////////////////////////////
std::vector<list_cutoff> excluded_volume_list (const dx_type &dx,
                                               const std::vector<coordinate> &pdb_coord) {

    // First generate a list of all grid points lying within 1.4A from the solute vdw surface
    // This is the SASA
    const double solv_r = 1.4;
    std::vector<list_cutoff> prelist = dx_list_cutoff (dx, pdb_coord, solv_r);

    // Set up a binary grid
    std::vector<bool> sasa_grid;
    sasa_grid.resize (dx.value.size(), 0);
    for (size_t i = 0; i < prelist.size(); i++)
        sasa_grid [prelist[i].index] = 1;
    // Save a copy for excluded volume
    std::vector<bool> excl_vol_grid = sasa_grid;

    // Now at every edge point (i.e. grid point lying just outside the above volume),
    // zeroing all grid points within 1.4 A from it. The remainder will be the excluded volume

    // Again, only need to consider local grid points, within 2 grid points from the SASA
    coordinate max_pdb, min_pdb;
    maxmin_coord_pdb (pdb_coord, max_pdb, min_pdb);
    const double max_r = 2.0 + solv_r;   // Maximum radius of solute atom

    size_t minx = std::max ((int) ceil  ((min_pdb.x - max_r - dx.origin[0]) / dx.delta[0]) - 2, 0);
    size_t maxx = std::min ((int) floor ((max_pdb.x + max_r - dx.origin[0]) / dx.delta[0]) + 2, (int) dx.ngrid[0] - 1);
    size_t miny = std::max ((int) ceil  ((min_pdb.y - max_r - dx.origin[1]) / dx.delta[1]) - 2, 0);
    size_t maxy = std::min ((int) floor ((max_pdb.y + max_r - dx.origin[1]) / dx.delta[1]) + 2, (int) dx.ngrid[1] - 1);
    size_t minz = std::max ((int) ceil  ((min_pdb.z - max_r - dx.origin[2]) / dx.delta[2]) - 2, 0);
    size_t maxz = std::min ((int) floor ((max_pdb.z + max_r - dx.origin[2]) / dx.delta[2]) + 2, (int) dx.ngrid[2] - 1);

    #pragma omp parallel for shared (sasa_grid, excl_vol_grid) collapse(3)
    for (size_t x = minx; x <= maxx; x++)
        for (size_t y = miny; y <= maxy; y++)
            for (size_t z = minz; z <= maxz; z++) {
                size_t index = dx_3Dindex_to_1Dindex (dx, x, y, z);
                if ((not sasa_grid[index]) and (isEdge_point (dx, sasa_grid, x, y, z)))
                    #pragma omp critical
                    empty_neighbor (excl_vol_grid, dx, x, y, z, solv_r);
            }

    std::vector<list_cutoff> result;
    size_t j = 0;
    for (size_t i = 0; i < excl_vol_grid.size(); i++)
        if (excl_vol_grid[i]) {
            while (i != prelist[j].index)
                j++;
            result.push_back (prelist[j]);
        }

    return result;
}

//////////////////////////////////////
// Excluded volume and hydration list
/////////////////////////////////////
void gen_exV_hyd_list (const std::vector<coordinate> &pdb_coord,
                       const dx_type &dx,
                       std::vector<list_cutoff> &list_exclV,
                       std::vector<list_cutoff> &list_hyd) {

    list_exclV = excluded_volume_list (dx, pdb_coord);
//    list_exclV = dx_list_cutoff (dx, pdb_coord, 1.4);
    size_t index = 0;
    for (size_t i = 0; i < dx.value.size(); i++)
        if ((index >= list_exclV.size()) or (i != list_exclV[index].index)) {
            list_cutoff tmp;
            tmp.index = i;
            list_hyd.push_back (tmp);
        } else index++;
}
