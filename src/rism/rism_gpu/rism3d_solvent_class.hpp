#ifndef solvent_class_HPP
#define solvent_class_HPP

#include <iostream>
#include "rism3d_solvent.hpp"
#include "rism3d_safemem.hpp"
#include "array_class.hpp"
#include "varTypes.hpp"
using namespace std;

namespace rism3d_c {

    class rism3d_solvent_class{
        private:
            rism3d_safemem* memalloc_p;

        public:
            double temperature;
            GPUtype dielconst;
            GPUtype xappa;
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
            
            // double* waveNumbers;
            array_class<GPUtype> waveNumbers;

            // double* xvv;
            array_class<GPUtype> xvv;

            double* xvv_dT;
            array_class<GPUtype> charge;
            array_class<double> charge_db;

            array_class<GPUtype> charge_sp;
            
            // double* density;
            array_class<GPUtype> density;
            array_class<double> density_db;
            
            double* density_sp;
            
            // double* ljSigma;
            array_class<GPUtype> ljSigma;

            // double* ljEpsilon;
            array_class<GPUtype> ljEpsilon;
            
            // double* coord;
            array_class<GPUtype> coord;
            
            double* background_correction;

            array_class<GPUtype> delhv0;
            
            double* delhv0_dT;
            bool ionic;
            double xvv_version;

            rism3d_solvent_class();
            rism3d_solvent_class(rism3d_safemem *memalloc, rism3d_c::solvent_cpp *solv_f);
            ~rism3d_solvent_class();

            void set_solvent(rism3d_c::solvent_cpp *solv_f);
    };
    
}

#endif