#ifndef solute_class_HPP
#define solute_class_HPP

#include <iostream>
#include "rism3d_solute.hpp"
// #include "rism3d_safemem.hpp"
#include "array_class.hpp"
#include "rism_util.hpp"
#include "varTypes.hpp"
using namespace std;

namespace rism3d_c {

    class rism3d_solute_class{
        private:
            rism3d_safemem* memalloc_p;

        public:
            //> Default origin
            GPUtype origin[3] = {0};

            //> Default translation
            GPUtype translation[3] = {0};

            int numAtoms;
            
            // double* mass;
            array_class<GPUtype> mass;
            
            // double* charge;
            array_class<GPUtype> charge;

            // double* origCharge;
            array_class<GPUtype> origCharge;
            
            // double* position;
            array_class<GPUtype> position;
            
            // double* ljSigma;
            array_class<GPUtype> ljSigma;

            // double* ljEpsilon;
            array_class<GPUtype> ljEpsilon;

            bool charged;
            
            // double totalCharge;
            GPUtype totalCharge;

            // double* centerOfMass;
            array_class<GPUtype> centerOfMass;

            rism3d_solute_class();
            rism3d_solute_class(rism3d_safemem *memalloc, rism3d_c::solute_cpp *solu_f);
            ~rism3d_solute_class();

            void set_solute(rism3d_c::solute_cpp *solu_f);
            void setCoord(double* solutePositions);
            void translate_solute(bool init_position = false);

    };

}

#endif