#ifndef RISM3DPOT_HPP
#define RISM3DPOT_HPP

#include <iostream>
#include <float.h>
#include "rism3d_grid.hpp"
#include "rism3d_safemem.hpp"
#include "rism3d_solute.hpp"
#include "rism3d_solute_class.hpp"
#include "rism3d_solvent_class.hpp"
#include "rism3d_solvent.hpp"
#include "array_class.hpp"
#include "rism_util.hpp"
#include "varTypes.hpp"
#include <cstring>

#include <fstream>
#include <iomanip>

using namespace std;

namespace rism3d_c {

    class rism3d_potential{
        private:
            rism3d_grid* grid_p;
            
            // rism3d_c::solute *solu;
            // rism3d_c::solvent *solv;

            rism3d_safemem* memalloc_p;
            rism3d_solute_class* soluteclass_p;
            rism3d_solvent_class* solventclass_p;

            double factor;
            int LJpower;
            double ljTolerance;
        public:
            // Constructor of pot object:
            // Still in progress. New arguments will be added as required
            // It is receiving solute and solvent direct from fortran because C++ objects were not
            // defined, yet. If we can call this constructor inside rism3d constructor, it will
            // allow us to copy solu_f/solv_f to rism3d members first. Then, we will no longer need
            // passing things from fortran.
            rism3d_potential(rism3d_safemem *memalloc, rism3d_grid *grid,
                                        rism3d_solute_class* soluteclass,
                                        rism3d_solvent_class* solventclass,
                                        double cut, char* periodic_pot, 
                                        bool treeDCF, bool treeTCF, bool treeCoulomb,
                                        GPUtype chargeSmear);

            ~rism3d_potential();

// #if RISMCUDA_DOUBLE
//             GPUtype PI = M_PI;
// #else
//             GPUtype PI = static_cast<float>(M_PI);
// #endif // RISMCUDA_DOUBLE

            GPUtype PI = std::acos(-1.0);

#if defined(RISMCUDA_DOUBLE)
            GPUtype MAX_VAL = DBL_MAX;
#else
            GPUtype MAX_VAL = FLT_MAX;
#endif //RISMCUDA_DOUBLE

            char* periodicPotential;

            bool treeDCF;
            bool treeTCF;
            bool treeCoulomb;

            array_class<GPUtype> uuv;

            double cutoff;
            double cutoff2;
            array_class<GPUtype> ljCutoffs2;

            array_class<GPUtype> ljSigmaUV;
            array_class<GPUtype> ljEpsilonUV;
            array_class<GPUtype> ljAUV;
            array_class<GPUtype> ljBUV;

            GPUtype cut2_chlk;
            GPUtype chargeSmear;

            array_class<GPUtype> dcfLongRangeAsympK;
            array_class<GPUtype> dcfLongRangeAsympR;
            array_class<double> dcfLongRangeAsympR_db;

            // Long-range part of Huv(k) at k = 0 (2, solv%natom)
            array_class<GPUtype> huvk0;

            array_class<GPUtype> tcfLongRangeAsympR;
            array_class<GPUtype> tcfLongRangeAsympK;
            array_class<double> tcfLongRangeAsympR_db;

            array_class<GPUtype> cached_h2;
            array_class<GPUtype> cached_hc;
            array_class<GPUtype> cached_hn;

            bool applyLJCorrection = false;

            bool periodic = false;

            GPUtype DCFTolerance;
            GPUtype DCFCoefficient;
            GPUtype smearSq;

            void long_range_asymptotics();
            void potential_calc();

            void setcut_ljdistance(double cut);
            void setcut_ljtolerance(double ljtolerance);
            double lennard_jones_cut_calc(int power, double tolerance);

            // lennard_jones_cut and lennard_jones_cut_deriv must be static to be
            // passed as function arguments
            // static double lennard_jones_cut(double r);
            // static double lennard_jones_cut_deriv(double r);

            void mixSoluteSolventLJParameters();

            void dcf_tcf_long_range_asymptotics(bool charged, bool periodic);

            void setcut_asympKTolerance(GPUtype asympKSpaceTolerance, GPUtype boxVolume);

            GPUtype asympck_cut_calc(GPUtype *soluteCharge, int solute_size, 
                                     GPUtype *solventCharge, int solvent_size, 
                                     GPUtype chargeSmear, GPUtype boxVolume, 
                                     GPUtype asympKSpaceTolerance);

            GPUtype asympck_cut(GPUtype k2);
            GPUtype asympck_cut_deriv(GPUtype k2);

            void dcf_long_range_asymptotics_R(int target_x_low_ind, int target_x_high_ind,
                                              int target_y_low_ind, int target_y_high_ind,
                                              int target_z_low_ind, int target_z_high_ind,
                                              GPUtype grid_xmin, GPUtype grid_ymin, GPUtype grid_zmin,
                                              GPUtype grid_spacing_x, GPUtype grid_spacing_y, GPUtype grid_spacing_z,
                                              int grid_dim_x, int grid_dim_y, int grid_dim_z,
                                              int solute_numAtoms, int solute_numAtoms_idx_start,  
                                              GPUtype *solute_position_x, GPUtype *solute_position_y, GPUtype *solute_position_z, GPUtype *solute_charge,
                                              GPUtype eta, GPUtype *dcf_long_range_asymptotics);

            void tcf_long_range_asymptotics_R(int target_x_low_ind, int target_x_high_ind,
                                              int target_y_low_ind, int target_y_high_ind,
                                              int target_z_low_ind, int target_z_high_ind,
                                              GPUtype grid_xmin, GPUtype grid_ymin, GPUtype grid_zmin,
                                              GPUtype grid_spacing_x, GPUtype grid_spacing_y, GPUtype grid_spacing_z,
                                              int grid_dim_x, int grid_dim_y, int grid_dim_z,
                                              int solute_numAtoms, int solute_numAtoms_idx_start,  
                                              GPUtype *solute_position_x, GPUtype *solute_position_y, GPUtype *solute_position_z, GPUtype *solute_charge,
                                              GPUtype xappa, GPUtype chargeSmear, GPUtype solvent_dielconst, GPUtype *tcf_long_range_asymptotics);

            // function to test the kernel results agains the serial implementation: I will leave this here for now
            void computeSumCosSin_serial(GPUtype* waveVectorX,
                                        GPUtype* waveVectorY,
                                        GPUtype* waveVectorZ,
                                        GPUtype* waveVectors2,
                                        GPUtype cut2_chlk,
                                        GPUtype* position,
                                        GPUtype* charge,
                                        GPUtype* dcfLongRangeAsympK,
                                        GPUtype asympk_const,
                                        GPUtype smear2_4,
                                        int numAtoms,
                                        int numWaveVectors_2,
                                        int start_ind);                                  

            GPUtype* sumcos_0;
            GPUtype* sumsin_0;
            void calc_dcf_tcf_LongRangeAsympK(GPUtype* waveVectorX,
                              GPUtype* waveVectorY,
                              GPUtype* waveVectorZ,
                              GPUtype* waveVectors2,
                              GPUtype cut2_chlk,
                              GPUtype* position,
                              GPUtype* charge,
                              GPUtype* dcfLongRangeAsympK,
                              GPUtype asympk_const,
                              GPUtype smear2_4,
                              int numAtoms,
                              int numWaveVectors_2,
                              int start_ind, 
                              bool ionic, 
                              GPUtype* tcfLongRangeAsympK,
                              GPUtype xappa2, 
                              GPUtype solvent_dielconst);
            
            void calc_sum_cos_sin_huvk0_serial(GPUtype waveVectorX_0,
                            GPUtype waveVectorY_0,
                            GPUtype waveVectorZ_0,
                            GPUtype* position,
                            GPUtype* charge,
                            int numAtoms,
                            GPUtype* sumcos_0,
                            GPUtype* sumsin_0);

            void calc_sum_cos_sin_huvk0(GPUtype waveVectorX_0,
                            GPUtype waveVectorY_0,
                            GPUtype waveVectorZ_0,
                            GPUtype* position,
                            GPUtype* charge,
                            int numAtoms,
                            GPUtype* sumcos_0,
                            GPUtype* sumsin_0);

            void int_h2_hc(GPUtype* h2, GPUtype* hc);

            void int_h(GPUtype* h);

    };

}

#endif // RISM3DGRID_HPP