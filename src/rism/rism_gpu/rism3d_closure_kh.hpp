#ifndef RISM3DCLOSUREKH_HPP
#define RISM3DCLOSUREKH_HPP

#include <iostream>
#include "rism3d_closure.hpp"
// #include <cub/cub.cuh> // requires C++14
// #include <thrust/device_vector.h> // requires C++14
// #include <thrust/reduce.h> // requires C++14
using namespace std;
#include "rism_util.hpp"

namespace rism3d_c {

    class kh : public rism3d_closure {     
        public:
            kh();
            kh(rism3d_solvent_class *solventclass, rism3d_solute_class *soluteclass, 
                           rism3d_grid *grid, rism3d_potential *pot, rism3d_safemem *memalloc);
            ~kh();

            void guv(GPUtype *uuv, GPUtype *guv, GPUtype *huv, GPUtype *cuv, int size) override;
            void cu_guv(GPUtype *uuv, GPUtype *guv, GPUtype *huv, GPUtype *cuv, int size);

            // array to store heaviside funtion calculated based on \Theta(-h(r))
            array_class<GPUtype> heaviside;

            // array to store excessChemicalPotential grid: already multiplied by density * voxelVolume;
            array_class<GPUtype> mu_grid;

            // array to store solvPotEne grid: already multiplied by density * voxelVolume;
            array_class<GPUtype> solvPotEne_grid;

            array_class<GPUtype> h_val;
            array_class<GPUtype> excessChemicalPotentialh2lr;
            array_class<GPUtype> excessChemicalPotentialhclr;

            // So far, we need the excChemPot_1d_db variable
            // to cast the result in double precision before
            // returning the result to amber_rism_interface.
            // This is because the wrapper functions for the
            // Fortran/C++ interface does not know about GPUtype
            // and I think that it would be too complicated
            // to mess with this at this moment: maybe in 
            // the future.
            GPUtype* excChemPot_1d;
            double* excChemPot_1d_db;

            GPUtype* solvPotEne_1d;
            double* solvPotEne_1d_db;

            // Function to calculate excessChemicalPotential grid points
            void cu_excessChemicalPotential_grid(GPUtype *heaviside, GPUtype *huv, GPUtype *cuv, GPUtype *mu_grid, 
                                                 int Nx, int Ny, int Nz, GPUtype density, GPUtype voxelVolume);

            void cu_excessChemicalPotential_lr_grid(GPUtype *heaviside, GPUtype *huv, GPUtype *cuv, GPUtype *mu_grid, 
                                                    int Nx, int Ny, int Nz, GPUtype density, GPUtype voxelVolume,
                                                    GPUtype solvent_charge, GPUtype* dcfLongRangeAsympR,
                                                    GPUtype solute_totalcharge, GPUtype solvent_charge_sp,
                                                    GPUtype* tcfLongRangeAsympR);

            // Wrap function to set heaviside on the GPU
            void cu_set_heaviside(GPUtype *heaviside, GPUtype *huv, int N);

            // Calculates the short range part of the excess chemical potential
            double* excessChemicalPotential(GPUtype *huv, GPUtype *cuv, int size);

            // Calculates the short range part of the excess chemical potential
            // with asymptotic correction in kT for each site.
            double* aexcessChemicalPotential(GPUtype *huv, GPUtype *cuv, int size, bool long_range_corr);

            // Calculates the solvent potential energy
            double* solvPotEne(GPUtype *guv, GPUtype *uuv, int size);

            // Function to calculate solvPotEne grid points
            void cu_solvPotEne_grid(GPUtype *guv, GPUtype *uuv, GPUtype *solvPotEne_grid, 
                                    int Nx, int Ny, int Nz,
                                    GPUtype density, GPUtype voxelVolume);

            array_class<GPUtype> cuv_int;
            double* cuv_int_dp;
            
            array_class<GPUtype> lj_corr;
            GPUtype partMolVol;
            double partMolVol_db;
            double partialMolarVolume(GPUtype* cuv);
            void DCFintegral(GPUtype* cuv);
            bool check_dbl_flt_max(GPUtype* ljCutoffs2, int solute_atoms, int solvent_atoms);
            void LJCorrection_DCF_int();

            GPUtype LJCorrection_TCF_int();

            double* get_cuv_int();

            array_class<GPUtype> kirkwoodBuff_1d;
            array_class<double> kirkwoodBuff_1d_db;
            double* kirkwoodBuff(GPUtype* huv, bool long_range_corr);

            void get_KB(GPUtype* huv, bool long_range_corr);

            array_class<GPUtype> excessPart_1d;
            array_class<double> excessPart_1d_db;

            array_class<GPUtype> excessPartBox_1d;
            array_class<double> excessPartBox_1d_db;

            double* excessParticles(GPUtype* huv);
            double* aexcessParticles(GPUtype* huv);

            GPUtype cu_akirkwoodBuff(GPUtype* huv, GPUtype solvent_charge_sp, GPUtype* tcfLongRangeAsympR, int Nx, int Ny, int Nz);




    };

}

#endif // RISM3DCLOSUREKH_HPP