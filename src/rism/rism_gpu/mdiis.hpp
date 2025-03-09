#ifndef MDIIS_HPP
#define MDIIS_HPP

#include <iostream>
#include "rism3d_grid.hpp"
#include "array_class.hpp"
#include "varTypes.hpp"
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include "rism3d_safemem.hpp"

using namespace std;

namespace rism3d_c {

    class mdiis{
        private:
            rism3d_safemem* memalloc_p;
            
            // total overlap matrix dimension
            int NVec;

            // current overlap matrix dimension: goes from 2 to NVec
            int current_NVec;
            
            GPUtype mdiis_del;

            // solution values
            GPUtype *xi;

            // residue values
            GPUtype *ri;

            // number of points in each solution array
            int np;

            // number of points considered in the real space
            int np_real;

            // handler for cublas
            cublasHandle_t handle;

            // handler for cusolver
            cusolverDnHandle_t handle_solv;

            // overlap matrix
            array_class<GPUtype> a;

            // b vector
            array_class<GPUtype> b;
            
        public:
            GPUtype residual;

            mdiis(int nvec, double del, rism3d_safemem *memalloc);

            ~mdiis();

            void resize(GPUtype *vecData, GPUtype *resVecData, int numAtomTypes, int numKx, int numKy, int numKz);

            void advance();

            int getWRKvec();

    };

}

#endif // MDIIS_HPP