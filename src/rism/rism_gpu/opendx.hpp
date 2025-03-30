#ifndef OPENDX_HPP
#define OPENDX_HPP

#include <iostream>
#include "array_class.hpp"
#include "varTypes.hpp"

using namespace std;

namespace rism3d_c {

    class opendx{
        private:
            
        public:
            opendx();

            ~opendx();

            void write_dx(string *file, int nx, int ny, int nz, 
                          GPUtype hx, GPUtype hy, GPUtype hz, 
                          GPUtype translation[3], GPUtype *u);

            void write_swapped_dx(string *file, int nx, int ny, int nz, 
                          GPUtype hx, GPUtype hy, GPUtype hz, 
                          GPUtype translation[3], GPUtype *u);
    };

}

#endif // OPENDX_HPP