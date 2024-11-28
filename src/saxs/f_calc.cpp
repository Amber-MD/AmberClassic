#include <complex>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "saxsDS.h"
#include "dx.h"
#include "atm_f.h"
#include "f_calc.h"
#include "const.h"

//////////////////////////////////////////////////
// Calculate the form factor of the pdb file A(q)
//////////////////////////////////////////////////
std::complex<double> form_factor (const std::vector<coordinate> &model,
                                  double q,
                                  const coordinate &Leb,
                                  double anom_f,
                                  bool flex,
                                  bool expli,
                                  const std::vector<coeff_f> &F_table) {

    std::complex<double> f (0, 0);
    for (size_t i = 0; i < model.size(); i++)
        if ((expli) or (model[i].type != "H")) { // Not account for H atoms for the implicit case
            double atomic_factor;
            atomic_factor = f_atm (model[i].type, q, model[i].nHyd, expli, 0, F_table);
            if ((model[i].type == "Rb+") or (model[i].type == "Cs+") or (model[i].type == "Br-"))
                atomic_factor += anom_f;
            if (flex)
                atomic_factor *= exp (-model[i].B_factor*q*q/(16*PI*PI));

            double qr = q * (Leb.x*model[i].x + Leb.y*model[i].y + Leb.z*model[i].z);
            f = f + atomic_factor * exp(std::complex<double> (0,1) * qr);
        }
    return f;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculate contribution of one particular region from one type of grid (water, counterions or electron).
// This region is specified by the list containing the indexes
// The boolean excess variable is used to specify whether the values are excess relatively to the bulk or not
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::complex<double> grid_factor (const dx_type &dx,
                                  int value_type,
                                  double atomic_factor,
                                  const std::vector<list_cutoff> &list,
                                  const coordinate &q_vector) {

    // 3D Fourier transform of the electron density in the unit volume
    // f = 8*rho*sin(q_x*a/2)*sin(q_y*a/2)*sin(q_z*a/2)/(q_x*q_y*q_z)
    double sinc_x, sinc_y, sinc_z;
    if (q_vector.x != 0)
        sinc_x = 2*sin(.5*q_vector.x*dx.delta[0]) / q_vector.x;
    else
        sinc_x = dx.delta[0];
    if (q_vector.y != 0)
        sinc_y = 2*sin(.5*q_vector.y*dx.delta[1]) / q_vector.y;
    else
        sinc_y = dx.delta[1];
    if (q_vector.z != 0)
        sinc_z = 2*sin(.5*q_vector.z*dx.delta[2]) / q_vector.z;
    else
        sinc_z = dx.delta[2];

    std::complex<double> sum (0, 0);

    for (size_t i = 0; i < list.size(); i++) {
        coordinate grid = dx_1Dindex2coord (dx, list[i].index);
        double value;

        if (value_type == 1)
            value = dx.value[list[i].index];
        else if (value_type == 0)
            value = dx.value[list[i].index] - 1;
        else if (value_type == -1)
            value = 1.;
        else {
            std::cout << "Invalid value_type   " << value_type << std::endl;
            exit (0);
        }

        double qr = q_vector.x*grid.x + q_vector.y*grid.y + q_vector.z*grid.z;
        sum = sum + value * exp(std::complex<double> (0, 1) * qr);
    }
    return atomic_factor * dx.conc * AVOGADRO * 1e-27 * sinc_x*sinc_y*sinc_z * sum;
}

////////////////////////////////////////////////////////////////////
// Map an excess electron map based on the atomic distribution grid
////////////////////////////////////////////////////////////////////
dx_type ex_elec (const std::vector<dx_type> &dx,
                 double anom_f,
                 bool expli,
                 const std::vector<coeff_f> &F_table) {

    dx_type result;
    result.ngrid = dx[0].ngrid;
    result.delta = dx[0].delta;
    result.origin = dx[0].origin;
    result.conc = 1;    // Concentration for each sites will be precalculated and stored
    result.type = "ex_e";
    result.value.resize(dx[0].value.size(), 0);

    for (size_t type = 0; type < dx.size(); type++)
        if ((expli) or (dx[type].type != "Hw")) {
            double Z;
            Z = f_atm (dx[type].type, 0, 0, expli, 1, F_table);
            if ((dx[type].type == "Rb+") or (dx[type].type == "Sr2+") or (dx[type].type == "Br-"))
                Z += anom_f;

            size_t nx = dx[type].ngrid[0];
            size_t ny = dx[type].ngrid[1];
            size_t nz = dx[type].ngrid[2];
            #pragma omp parallel for shared (result, type, Z) collapse(3)
            for (size_t x = 0; x < nx; x++)
                for (size_t y = 0; y < ny; y++)
                    for (size_t z = 0; z < nz; z++) {
                        size_t index_center = dx_3Dindex_to_1Dindex (dx[type], x, y, z);
                        double grid = (dx[type].value[index_center] - 1) * dx[type].conc * Z;
                        result.value[index_center] += grid;
                    }
        }
    return result;
}
