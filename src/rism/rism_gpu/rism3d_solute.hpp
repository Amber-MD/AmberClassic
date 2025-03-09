#ifndef solute_struc_HPP
#define solute_struc_HPP
// #include "secondClass.hpp"
#include <iostream>
using namespace std;

namespace rism3d_c{

    struct solute_cpp{
      int numAtoms;
      double* mass;
      double* charge;
      double* origCharge;
      double* position;
      double* ljSigma;
      double* ljEpsilon;
      bool charged;
      double totalCharge;
      double* centerOfMass;
      double* translation;  
    };
    typedef struct solute_cpp solute;

}

// struct solute{
//   int numAtoms;
//   double* mass;
//   double* charge;
//   double* origCharge;
//   double* position;
//   double* ljSigma;
//   double* ljEpsilon;
//   bool charged;
//   double totalCharge;
//   double* centerOfMass;
//   double* translation;  
// };
// typedef struct solute solute;

#endif