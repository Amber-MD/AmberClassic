#include <iostream>
#include "rism3d_fft.hpp"
using namespace std;

namespace rism3d_c {

    rism3d_fft :: rism3d_fft(rism3d_grid *grid) {
        // cout << "FFT object created" << endl;

        grid_p = grid;
        
    }

    rism3d_fft :: ~rism3d_fft(){
        // cout << "FFT object destroyed" << endl;
    }

    void rism3d_fft :: setgrid(int ngr[3], GPUtype grdspc[3], int narray){
        int nrz = ngr[3];
        int nrzOff = 0;
        int nky = ngr[2];
        int nkyOff = 0;
        // old column major version:
        // int nkTotal = (ngr[1] + 2) * ngr[2] * ngr[3];

        // new row major version:
        int nkTotal = ngr[0] * ngr[1] * (ngr[2] + 2);

        int globalDimsK[3];
        int localDimsR[3];
        int localDimsK[3];
        int offsetR[3] = {0,0,0};
        int offsetK[3] = {0,0,0};

        // old column major version:
        // globalDimsK[0] = ngr[0]+2;
        // globalDimsK[1] = ngr[1];
        // globalDimsK[2] = ngr[2];

        // new row major version:
        globalDimsK[0] = ngr[0];
        globalDimsK[1] = ngr[1];
        globalDimsK[2] = ngr[2]+2;

        localDimsR[0] = ngr[0];
        localDimsR[1] = ngr[1];
        localDimsR[2] = ngr[2];

        // old column major version:
        // localDimsK[0] = ngr[0]+2;
        // localDimsK[1] = ngr[1];
        // localDimsK[2] = ngr[2];

        // new row major version:
        localDimsK[0] = ngr[0];
        localDimsK[1] = ngr[1];
        localDimsK[2] = ngr[2]+2;

        grid_p->resize(grdspc, ngr, globalDimsK, localDimsR, localDimsK, offsetR, offsetK);

        // Printing some parameters for checking
        // for(int i = 0; i < 3; i++){
        //     cout << "localDimsR[" << i << "] = " << grid_p->localDimsR[i] << endl;
        //     cout << "globalDimsR[" << i << "] = " << grid_p->globalDimsR[i] << endl;
        //     cout << "localDimsK[" << i << "] = " << grid_p->localDimsK[i] << endl;
        //     cout << "globalDimsK[" << i << "] = " << grid_p->globalDimsK[i] << endl;
        // }

    }

    void rism3d_fft :: create_plans(int Nx, int Ny, int Nz){
#if defined(RISMCUDA_DOUBLE)
        cufftPlan3d(&plan, Nx, Ny, Nz, CUFFT_D2Z);
        cufftPlan3d(&plan_inv, Nx, Ny, Nz, CUFFT_Z2D);
#else
        cufftPlan3d(&plan, Nx, Ny, Nz, CUFFT_R2C);
        cufftPlan3d(&plan_inv, Nx, Ny, Nz, CUFFT_C2R);
#endif // RISMCUDA_DOUBLE
    }

    void rism3d_fft :: destroy_plans(){
        cufftDestroy(plan);
        cufftDestroy(plan_inv);
    }  

}