#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <iomanip>

#include "saxsDS.h"
#include "dx.h"

// ********************************************
// ***            Read dx file              ***
// ********************************************
void read_dx (std::string &dx_file,
              dx_type &dx) {

    double tmp1, tmp2, tmp3;
    std::ifstream DXFILE (dx_file.c_str());
    if (DXFILE.is_open()) {
        std::string line;
        while (getline(DXFILE, line)) {
            std::istringstream iss(line);
            while (iss >> tmp1)
                dx.value.push_back(tmp1);
            if (line.find("object 1 class gridpositions counts") != std::string::npos) {
                line.replace(0,35," ");
                std::istringstream iss(line);
                dx.ngrid.resize(3);
                iss >> dx.ngrid[0] >> dx.ngrid[1] >> dx.ngrid[2];
            } else if (line.find("origin") != std::string::npos) {
                line.replace(0,6," ");
                std::istringstream iss(line);
                dx.origin.resize(3);
                iss >> dx.origin[0] >> dx.origin[1] >> dx.origin[2];
            } else if (line.find("delta") != std::string::npos) {
                line.replace(0,5," ");
                std::istringstream iss(line);
                dx.delta.resize(3);
                iss >> tmp1 >> tmp2 >> tmp3;
                if (tmp1 != 0)
                    dx.delta[0] = tmp1;
                else if (tmp2 != 0)
                    dx.delta[1] = tmp2;
                else if (tmp3 != 0)
                    dx.delta[2] = tmp3;
            }
        }
        DXFILE.close();
        if (dx.ngrid[0] * dx.ngrid[1] * dx.ngrid[2] != dx.value.size()) {
            std::cout << "Number of grid points in " << dx_file << " not matched!!\n";
            std::cout << "Number of grid points (x y z):\t" << dx.ngrid[0] << "\t" << dx.ngrid[1] << "\t" << dx.ngrid[2] << std::endl;
            std::cout << "Number of values:\t" << dx.value.size() << std::endl;
            exit (0);
        }
    } else {
        std::cerr << "Unable to open file " << dx_file << std::endl;
        exit (0);
    }
}

// ********************************************
// ***           Write dx file              ***
// ********************************************
void write_dx (const dx_type &dx,
               std::string &outfile) {

    std::ofstream OUT (outfile.c_str());
    if (OUT.is_open()) {
        OUT << "object 1 class gridpositions counts" << std::setw(8) << dx.ngrid[0] << std::setw(8) << dx.ngrid[1] << std::setw(8) << dx.ngrid[2] << std::endl;
        OUT << "origin " << std::setw(15) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.origin[0] <<\
            std::setw(15) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.origin[1] <<\
            std::setw(15) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.origin[2] << std::endl;
        OUT << "delta " << std::setw(16) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.delta[0] << " 0 0\n";
        OUT << "delta  0" << std::setw(16) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.delta[1] << " 0\n";
        OUT << "delta  0 0" << std::setw(16) << std::setprecision(8) << std::setiosflags(std::ios::fixed) << dx.delta[2] << std::endl;
        OUT << "object 2 class gridconnections counts" << std::setw(8) << dx.ngrid[0] << std::setw(8) << dx.ngrid[1] <<\
            std::setw(8) << dx.ngrid[2] << std::endl;
        OUT << "object 3 class array type double rank 0 items " << dx.ngrid[0]*dx.ngrid[1]*dx.ngrid[2] << " data follows";

        for (size_t i = 0; i < dx.value.size(); i++) {
            if (i % 3 == 0) OUT << std::endl;
            OUT << std::setw(16) << std::setprecision(5) << std::setiosflags(std::ios::fixed) << std::scientific << dx.value[i];
        }

        OUT << "\nobject \"Untitled\" call field";
        OUT.close();
    } else {
        std::cout << "Can't write to file " << outfile << std::endl;
        exit (0);
    }
}

/////////////////////////////////////////////
// Check to see dx grids are different or not
/////////////////////////////////////////////
bool check_dx (const std::vector<dx_type> &dx) {
    bool check = 0;
    for (size_t i = 1; i < dx.size(); i++) {
        if (dx[i].ngrid != dx[0].ngrid) {
            check = 1;
            break;
        }
        if (dx[i].origin != dx[0].origin) {
            check = 1;
            break;
        }
        if (dx[i].delta != dx[0].delta) {
            check = 1;
            break;
        }
    }
    return check;
}

// *****************************************************************
// *****    Get grid point coordinates from the 1d index       *****
// *****************************************************************
coordinate dx_1Dindex2coord (const dx_type &dx,
                             int index) {

    int x_index = floor (index/(dx.ngrid[1]*dx.ngrid[2]));
    int y_index = floor ((index - x_index * dx.ngrid[1] * dx.ngrid[2]) / dx.ngrid[2]);
    int z_index = index - x_index*dx.ngrid[1]*dx.ngrid[2] - y_index*dx.ngrid[2];

    coordinate coord;
    coord.x = dx.origin[0] + dx.delta[0] * x_index;
    coord.y = dx.origin[1] + dx.delta[1] * y_index;
    coord.z = dx.origin[2] + dx.delta[2] * z_index;

    return coord;
}

// ***************************************************
// *****    Return the 1D index from 3D index  *******
// ***************************************************
int dx_3Dindex_to_1Dindex (const dx_type &dx,
                           int x, int y, int z) {
    return x*dx.ngrid[1]*dx.ngrid[2] + y*dx.ngrid[2] + z;
}

//////////////////////////////////////////////////////////////////////////////
// Compute distance squared between two grid points, assuming orthogonal grids
//////////////////////////////////////////////////////////////////////////////
double grid2grid_dist (const dx_type &dx,
                       size_t index1, size_t index2) {

    coordinate coord1 = dx_1Dindex2coord (dx, index1);
    coordinate coord2 = dx_1Dindex2coord (dx, index2);
    double delx = coord1.x - coord2.x;
    double dely = coord1.y - coord2.y;
    double delz = coord1.z - coord2.z;
    return sqrt(delx*delx + dely*dely + delz*delz);
}

// ***************************************************************
// ***      Extract a smaller dx from the original dx          ***
// ***************************************************************
void extract_dx (const dx_type &dx,
                 dx_type &dx_out,
                 coordinate &min,
                 coordinate &max) {

    dx_out.delta = dx.delta;

    // Get the small and large indexes
    size_t imin = ceil ((min.x - dx.origin[0])/dx.delta[0]);
    size_t jmin = ceil ((min.y - dx.origin[1])/dx.delta[1]);
    size_t kmin = ceil ((min.z - dx.origin[2])/dx.delta[2]);

    dx_out.origin.resize (3);
    dx_out.origin[0] = dx.origin[0] + imin*dx.delta[0];
    dx_out.origin[1] = dx.origin[1] + jmin*dx.delta[1];
    dx_out.origin[2] = dx.origin[2] + kmin*dx.delta[2];

    size_t imax = floor ((max.x - dx.origin[0])/dx.delta[0]);
    size_t jmax = floor ((max.y - dx.origin[1])/dx.delta[1]);
    size_t kmax = floor ((max.z - dx.origin[2])/dx.delta[2]);

    dx_out.ngrid.resize (3);
    dx_out.ngrid[0] = imax - imin + 1;
    dx_out.ngrid[1] = jmax - jmin + 1;
    dx_out.ngrid[2] = kmax - kmin + 1;

    // Loop between the small and large indexes to get value for the extracted dx
    for (size_t i = imin; i <= imax; i++)
        for (size_t j = jmin; j <= jmax; j++)
            for (size_t k = kmin; k <= kmax; k++) {
                size_t index = (size_t) i*dx.ngrid[1]*dx.ngrid[2] + j*dx.ngrid[2] + k;
                dx_out.value.push_back(dx.value[index]);
            }

    // Sanity check
    if (dx_out.ngrid[0] * dx_out.ngrid[1] * dx_out.ngrid[2] != dx_out.value.size()) {
        std::cerr << "Number of grid point not matched\n";
        exit (0);
    }
}

// *****************************************************************
// ****         Adjust values of the dx: new = a*old + b         ***
// *****************************************************************
dx_type scale_dx (const dx_type &dx,
                  double a,
                  double b) {
    dx_type result;
    result.value.resize (dx.value.size());
    result.ngrid = dx.ngrid;
    result.origin = dx.origin;
    result.delta = dx.delta;

    for (size_t i = 0; i < dx.value.size(); i++)
        result.value[i] = a*dx.value[i] + b;

    return result;
}

// *************************************************************
// ****        Integrate all grid points of a dx       *********
// *************************************************************
double integral_dx (const dx_type &dx) {
    double result = 0;
    for (size_t i = 0; i < dx.value.size(); i++)
        result += dx.value[i];
    return result;
}

/////////////////////////////////////////////////////////////////////
//  Initiating a grid to enclose a region with a certain resolution
//      min and max are two extreme voxels
/////////////////////////////////////////////////////////////////////
dx_type gen_grid (const coordinate &max_coord,
                  const coordinate &min_coord,
                  const coordinate &res) {
    dx_type grid;
    grid.origin.resize(3);
    grid.origin[0] = min_coord.x;  grid.origin[1] = min_coord.y;  grid.origin[2] = min_coord.z;

    grid.ngrid.resize(3);
    grid.ngrid[0] = ceil ((max_coord.x - min_coord.x) / res.x);
    grid.ngrid[1] = ceil ((max_coord.y - min_coord.y) / res.y);
    grid.ngrid[2] = ceil ((max_coord.z - min_coord.z) / res.z);

    grid.delta.resize(3);
    grid.delta[0] = (max_coord.x - min_coord.x) / grid.ngrid[0];
    grid.delta[1] = (max_coord.y - min_coord.y) / grid.ngrid[1];
    grid.delta[2] = (max_coord.z - min_coord.z) / grid.ngrid[2];

    grid.conc = 1;
    grid.value.resize(grid.ngrid[0] * grid.ngrid[1] * grid.ngrid[2], 0);
    return grid;
}

//////////////////////////////
/// Initiating a coarser grid 
///////////////////////////////
dx_type gen_coarse_dx (const dx_type &dx,
                       size_t merge_x, size_t merge_y, size_t merge_z) {
    dx_type result;
    result.origin = dx.origin;
    result.type = dx.type;
    result.conc = dx.conc;

    result.delta.resize(3);
    result.delta[0] = dx.delta[0] * merge_x;
    result.delta[1] = dx.delta[1] * merge_y;
    result.delta[2] = dx.delta[2] * merge_z;

    result.ngrid.resize(3);
    result.ngrid[0] = floor (dx.ngrid[0]/merge_x);
    result.ngrid[1] = floor (dx.ngrid[1]/merge_y);
    result.ngrid[2] = floor (dx.ngrid[2]/merge_z);

    result.value.resize (result.ngrid[0]*result.ngrid[1]*result.ngrid[2], 0);
    return result;
}

/////////////////////////
// Coarse-graining grid
////////////////////////
void coarse_grid (const dx_type &fine_dx,
                  dx_type &coarse_dx,
                  size_t merge_x, size_t merge_y, size_t merge_z) {

    #pragma omp parallel for shared (coarse_dx, merge_x, merge_y, merge_z)
    for (size_t i = 0; i < coarse_dx.ngrid[0]*merge_x; i++)
        for (size_t j = 0; j < coarse_dx.ngrid[1]*merge_y; j++)
            for (size_t k = 0; k < coarse_dx.ngrid[2]*merge_z; k++) {
                size_t x = floor (i/merge_x);
                size_t y = floor (j/merge_y);
                size_t z = floor (k/merge_z);
                size_t index_old = dx_3Dindex_to_1Dindex (fine_dx, i, j, k);
                size_t index_new = dx_3Dindex_to_1Dindex (coarse_dx, x, y, z);
                coarse_dx.value[index_new] += fine_dx.value[index_old];
            }

    for (size_t i = 0; i < coarse_dx.value.size(); i++)
        coarse_dx.value[i] /= merge_x * merge_y * merge_z;
}
