#ifndef RISM3DFFT_HPP
#define RISM3DFFT_HPP

#include <iostream>
#include "rism3d_grid.hpp"
#include "array_class.hpp"
#include "varTypes.hpp"
#include <cufft.h>

using namespace std;

namespace rism3d_c {

    class rism3d_fft{
        private:
            // aligned :: Needed for SIMD.
            //            FFTW_ALIGNED -   assume that the memory is 16-byte aligned.
            //            FFTW_UNALIGNED - assume that memory is not aligned.
            int aligned;

            // localtrans :: transpose data locally.  Not used yet.
            bool localtrans;

            //grid  :: grid object with array dimensions
            rism3d_grid* grid_p;
            
            //narray :: the number of arrays to transform
            int narray;
            
            // wrk :: work space to transform to and from Numerical Recipes memory layout
            // wrk_trans :: work space to do local transpose
            // array_class<double> wrk; // maybe use this version...
            // array_class<double> wrk_trans; // maybe use this version...
            double* wrk;
            double* wrk_trans;
            
        public:
            rism3d_fft(rism3d_grid* grid);

            ~rism3d_fft();

            cufftHandle plan;
            cufftHandle plan_inv;

            void setgrid(int ngr[3], GPUtype grdspc[3], int narray);

            void create_plans(int Nx, int Ny, int Nz);
            void destroy_plans();

    };

}

#endif // RISM3DFFT_HPP