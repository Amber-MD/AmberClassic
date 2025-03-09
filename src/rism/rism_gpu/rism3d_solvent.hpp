#ifndef solvent_struc_HPP
#define solvent_struc_HPP
// #include "secondClass.hpp"
#include <iostream>
using namespace std;

namespace rism3d_c{

    struct solvent_cpp{
      double temperature;
      double dielconst;
      double xappa;
      double xikt;
      double smear;
      double xikt_dT;
      int numAtomTypes;
      int numMolecules;
      int numRDFpoints;
      int* atomMultiplicity;
      int* numAtoms;
      string* atomName;
      // char[][] atomName2;
      double gridSpacingR;
      double gridSpacingK;
      double* waveNumbers;
      double* xvv;
      double* xvv_dT;
      double* charge;
      double* charge_sp;
      double* density;
      double* density_sp;
      double* ljSigma;
      double* ljEpsilon;
      double* coord;
      double* background_correction;
      double* delhv0;
      double* delhv0_dT;
      bool ionic;
      double xvv_version;
    };
    typedef struct solvent_cpp solvent;

}

// struct solvent{
//   double temperature;
//   double dielconst;
//   double xappa;
//   double xikt;
//   double smear;
//   double xikt_dT;
//   int numAtomTypes;
//   int numMolecules;
//   int numRDFpoints;
//   int* atomMultiplicity;
//   int* numAtoms;
//   string* atomName;
//   double gridSpacingR;
//   double gridSpacingK;
//   double* waveNumbers;
//   double* xvv;
//   double* xvv_dT;
//   double* charge;
//   double* charge_sp;
//   double* density;
//   double* density_sp;
//   double* ljSigma;
//   double* ljEpsilon;
//   double* coord;
//   double* background_correction;
//   double* delhv0;
//   double* delhv0_dT;
//   bool ionic;
//   double xvv_version;
// };
// typedef struct solvent solvent;

#endif