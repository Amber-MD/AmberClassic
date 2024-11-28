#ifndef SAXSDS_H
#define SAXSDS_H

#include <vector>
#include <string>

struct coeff_f {
    std::string atomtype;
    double a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6;
};

struct coordinate {
    // Atom type
    std::string type;
    double x, y, z;
    double r;
    double B_factor;
    // Number of hydrogen atoms attached
    size_t nHyd;
};

struct dx_type {
    std::string type;
    double conc;
    std::vector<double> value;
    std::vector<size_t> ngrid;
    std::vector<double> origin;
    std::vector<double> delta;
};

struct list_cutoff {
    size_t index;
    double dist;        // Distance to the nearest solute atom
};

typedef std::pair<double, double> DoubPair;

#endif
