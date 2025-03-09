#include <iostream>
#include "rism3d_potential.hpp"
#include <math.h>
using namespace std;

namespace rism3d_c {

    rism3d_potential :: rism3d_potential(rism3d_safemem *memalloc, rism3d_grid *grid,
                                        rism3d_solute_class* soluteclass, rism3d_solvent_class* solventclass, 
                                        double cut, char* periodic_pot, bool treeDCF, bool treeTCF, bool treeCoulomb,
                                        GPUtype chargeSmear) :
                                        periodicPotential(periodic_pot),
                                        treeDCF(treeDCF),
                                        treeTCF(treeTCF),
                                        treeCoulomb(treeCoulomb),
                                        uuv(array_class<GPUtype>::ROW_MAJOR),
                                        ljSigmaUV(array_class<GPUtype>::ROW_MAJOR),
                                        ljEpsilonUV(array_class<GPUtype>::ROW_MAJOR),
                                        ljAUV(array_class<GPUtype>::ROW_MAJOR),
                                        ljBUV(array_class<GPUtype>::ROW_MAJOR),
                                        chargeSmear(chargeSmear),
                                        dcfLongRangeAsympK(array_class<GPUtype>::ROW_MAJOR),
                                        dcfLongRangeAsympR(array_class<GPUtype>::ROW_MAJOR),
                                        dcfLongRangeAsympR_db(array_class<double>::ROW_MAJOR),
                                        huvk0(array_class<GPUtype>::ROW_MAJOR),
                                        tcfLongRangeAsympR(array_class<GPUtype>::ROW_MAJOR),
                                        tcfLongRangeAsympK(array_class<GPUtype>::ROW_MAJOR)
                                        {
        // cout << "Pot object created" << endl;

        grid_p = grid;
        memalloc_p = memalloc;
        soluteclass_p = soluteclass;
        solventclass_p = solventclass;

        uuv.set_memalloc(memalloc_p, true);

        ljCutoffs2.set_memalloc(memalloc_p);
        ljCutoffs2.alloc_mem(soluteclass_p->numAtoms, solventclass_p->numAtomTypes);
        
        ljSigmaUV.set_memalloc(memalloc_p);
        ljSigmaUV.alloc_mem(soluteclass_p->numAtoms, solventclass_p->numAtomTypes);

        ljEpsilonUV.set_memalloc(memalloc_p);
        ljEpsilonUV.alloc_mem(soluteclass_p->numAtoms, solventclass_p->numAtomTypes);
        
        ljAUV.set_memalloc(memalloc_p, true);
        ljAUV.alloc_mem(soluteclass_p->numAtoms, solventclass_p->numAtomTypes);

        ljBUV.set_memalloc(memalloc_p, true);
        ljBUV.alloc_mem(soluteclass_p->numAtoms, solventclass_p->numAtomTypes);

        dcfLongRangeAsympK.set_memalloc(memalloc_p, true);

        dcfLongRangeAsympR.set_memalloc(memalloc_p, true);

        // If using float number, we need to have a double variable to pass 
        // values back to Fortran
#if !RISMCUDA_DOUBLE
            dcfLongRangeAsympR_db.set_memalloc(memalloc_p, true);
            tcfLongRangeAsympR_db.set_memalloc(memalloc_p, true);
#endif // RISMCUDA_DOUBLE

        huvk0.set_memalloc(memalloc_p, true);

        tcfLongRangeAsympR.set_memalloc(memalloc_p, true);

        tcfLongRangeAsympK.set_memalloc(memalloc_p, true);

        cached_h2.set_memalloc(memalloc_p, true);
        cached_h2.alloc_mem(solventclass_p->numAtomTypes);

        cached_hc.set_memalloc(memalloc_p, true);
        cached_hc.alloc_mem(solventclass_p->numAtomTypes);

        cached_hn.set_memalloc(memalloc_p, true);
        cached_hn.alloc_mem(solventclass_p->numAtomTypes);

        // cached_h2/hc/hn are too small to worth initializing
        // them on the GPU: so we can initialize them on the CPU
        // and prefetch the values before sending them to a cuda
        // kernel
        for(int iv = 0; iv < solventclass_p->numAtomTypes; iv++){
#if RISMCUDA_DOUBLE
            cached_h2.m_data[iv] = DBL_MAX;
            cached_hc.m_data[iv] = DBL_MAX;
            cached_hn.m_data[iv] = DBL_MAX;
#else
            cached_h2.m_data[iv] = FLT_MAX;
            cached_hc.m_data[iv] = FLT_MAX;
            cached_hn.m_data[iv] = FLT_MAX;
#endif // RISMCUDA_DOUBLE
        }

        setcut_ljdistance(cut);
        mixSoluteSolventLJParameters();
    }

    rism3d_potential :: ~rism3d_potential(){
        // cout << "Pot object destroyed" << endl;
    }

    void rism3d_potential :: setcut_ljdistance(double cut){
        cutoff = min(sqrt(DBL_MAX), cut);
        cutoff2 = cutoff*cutoff;

        for(int i = 0; i < soluteclass_p->numAtoms; i++){
            for(int j = 0; j < solventclass_p->numAtomTypes; j++){
                ljCutoffs2(i, j) = cutoff2;
            }
        }

    }

    void rism3d_potential :: setcut_ljtolerance(double tolerance){
        if(tolerance == 0){
            for(int i = 0; i < soluteclass_p->numAtoms; i++){
                for(int j = 0; j < solventclass_p->numAtomTypes; j++){
#if defined(RISMCUDA_DOUBLE)
                    ljCutoffs2(i, j) = DBL_MAX;
#else
                    ljCutoffs2(i, j) = FLT_MAX;
#endif //RISMCUDA_DOUBLE
                }
            }
        }
        else{
            for(int i = 0; i < soluteclass_p->numAtoms; i++){
                for(int j = 0; j < solventclass_p->numAtomTypes; j++){
                    applyLJCorrection = true;
                    cout << "Still need to implement the case ljTolerance != 0!!!" << endl;
                    // Would be some call like this:
                    // ljCutoffs2(i, j) = lennard_jones_cut_calc(6, tolerance);
                    abort();
                }
            }
        }

    }

    // Functions to be implemented
    // double rism3d_potential :: lennard_jones_cut_calc(int power, double tolerance){
    //     double cutoff;
    //     LJpower = power;
    //     ljTolerance = tolerance;
    //     factor = 1.0;
    //     if(power == 6){
    //         factor = 2.0;
    //     }

    //     cutoff = root_newton(&rism3d_c::rism3d_potential::lennard_jones_cut);

    //     return cutoff;
    // }

    // double rism3d_potential :: lennard_jones_cut(double r){
    //     return pow(r,2);
    // }

    // double rism3d_potential :: lennard_jones_cut_deriv(double r){
    //     return r;
    // }

    void rism3d_potential :: long_range_asymptotics(){
        // cout << "In rism3d_pot class" << endl;
        // cout << "grid_p->localDimsR[0] = " << grid_p->localDimsR[0] << endl;
        // cout << grid_p->localDimsR[0] << endl;
        // cout << grid_p->localDimsR[1] << endl;
        // cout << grid_p->localDimsR[2] << endl;

        // tcfLongRangeAsympR = memalloc_p->allocReal(tcfLongRangeAsympR,grid_p->localDimsR[0]*grid_p->localDimsR[1]*grid_p->localDimsR[2]);
        // dcfLongRangeAsympR = memalloc_p->allocReal(tcfLongRangeAsympR,grid_p->localDimsR[0]*grid_p->localDimsR[1]*grid_p->localDimsR[2]);

        // int count = 0;
        
        // for(int i = 0; i < grid_p->localDimsR[0]*grid_p->localDimsR[1]*grid_p->localDimsR[2]; i++){
        //     tcfLongRangeAsympR[i] = count;
        //     dcfLongRangeAsympR[i] = count+10;
        //     count++;
        // }
    }

    void rism3d_potential :: mixSoluteSolventLJParameters(){
        for(int i = 0; i < soluteclass_p->numAtoms; i++){
            for(int j = 0; j < solventclass_p->numAtomTypes; j++){
                
                ljSigmaUV(i,j) = soluteclass_p->ljSigma(i) + solventclass_p->ljSigma(j);

                ljEpsilonUV(i,j) = sqrt(soluteclass_p->ljEpsilon(i)*solventclass_p->ljEpsilon(j));

                ljAUV(i, j) = ljEpsilonUV(i, j)*pow(ljSigmaUV(i,j), 12);
                ljBUV(i, j) = 2.0*ljEpsilonUV(i, j)*pow(ljSigmaUV(i,j), 6);
            }
        }

// This was for priting A and B coefficients: I will leave it here for now in case I want to check these
// in the future.
// #if RISMCUDA_DOUBLE
//         ofstream myfile1;
//         myfile1.open("/home/fcarvalho/rism3d.cuda.test/Auv_db.txt");
//         for(int i = 0; i < soluteclass_p->numAtoms; i++){
//             for(int j = 0; j < solventclass_p->numAtomTypes; j++){
//                 myfile1 << setprecision(17) << i << "," << j << "," << ljAUV(i,j) << "\n";
//             }
//         }
//         myfile1.close();

//         ofstream myfile2;
//         myfile2.open("/home/fcarvalho/rism3d.cuda.test/Buv_db.txt");
//         for(int i = 0; i < soluteclass_p->numAtoms; i++){
//             for(int j = 0; j < solventclass_p->numAtomTypes; j++){
//                 myfile2 << setprecision(17) << i << "," << j << "," << ljBUV(i,j) << "\n";
//             }
//         }
//         myfile2.close();
// #else
//         ofstream myfile1;
//         myfile1.open("/home/fcarvalho/rism3d.cuda.test/Auv.txt");
//         for(int i = 0; i < soluteclass_p->numAtoms; i++){
//             for(int j = 0; j < solventclass_p->numAtomTypes; j++){
//                 myfile1 << setprecision(17) << i << "," << j << "," << ljAUV(i,j) << "\n";
//             }
//         }
//         myfile1.close();

//         ofstream myfile2;
//         myfile2.open("/home/fcarvalho/rism3d.cuda.test/Buv.txt");
//         for(int i = 0; i < soluteclass_p->numAtoms; i++){
//             for(int j = 0; j < solventclass_p->numAtomTypes; j++){
//                 myfile2 << setprecision(17) << i << "," << j << "," << ljBUV(i,j) << "\n";
//             }
//         }
//         myfile2.close();
// #endif //RISMCUDA_DOUBLE

    }

    void  rism3d_potential :: dcf_tcf_long_range_asymptotics(bool charged, bool periodic){
        // Check if system is neutral or non-periodic
        if(charged == false && periodic == false){
            return;
        }

        // Making sure the grid size was set
        if(grid_p->waveVectorX.m_data == nullptr && grid_p->waveVectorY.m_data == nullptr && grid_p->waveVectorZ.m_data == nullptr){
            cout << "In dcf_tcf_long_range_asymptotics: grid size not set" << endl;
            abort();    
        }

        if(grid_p->totalLocalPointsK != dcfLongRangeAsympK.get_dims(0) || 
           grid_p->localDimsR[0]*grid_p->localDimsR[1]*grid_p->localDimsR[2] != dcfLongRangeAsympR.get_dims(0)){
            // Allocating memory for dcfLongRangeAsympK and dcfLongRangeAsympR
            dcfLongRangeAsympK.alloc_mem(grid_p->totalLocalPointsK);
            dcfLongRangeAsympR.alloc_mem(grid_p->localDimsR[0]*grid_p->localDimsR[1]*grid_p->localDimsR[2]);

            // checking if solving is ionic
            if(solventclass_p->ionic == true){
                // Allocating memory for tcfLongRangeAsympR and tcfLongRangeAsympK
                tcfLongRangeAsympR.alloc_mem(grid_p->localDimsR[0]*grid_p->localDimsR[1]*grid_p->localDimsR[2]);
                tcfLongRangeAsympK.alloc_mem(grid_p->totalLocalPointsK);
            }
        }

        // Allocate long range part of Huv(k) if not done already.
        if(huvk0.m_data == nullptr){
            // Need to check if it is needed to change the order of
            // dimentions here
            huvk0.alloc_mem(solventclass_p->numAtomTypes,2);
            // the huvk0_dT should be allocated here as well.
        }

        // Set some values for the long range asymptotics
        GPUtype smear2_4 = chargeSmear * chargeSmear / 4.0;
        GPUtype xappa2 = solventclass_p->xappa * solventclass_p->xappa;
        GPUtype asympk_const = 4 * PI / grid_p->boxVolume;

        // Initialize correlation function arrays.
        cudaMemset(dcfLongRangeAsympK.m_data, 0, grid_p->totalLocalPointsK*sizeof(GPUtype));
        cudaMemset(dcfLongRangeAsympR.m_data, 0, grid_p->localDimsR[0]*grid_p->localDimsR[1]*grid_p->localDimsR[2]*sizeof(GPUtype));

        if(solventclass_p->ionic == true){
            // Initialize correlation function arrays.
            cudaMemset(tcfLongRangeAsympR.m_data, 0, grid_p->localDimsR[0]*grid_p->localDimsR[1]*grid_p->localDimsR[2]*sizeof(GPUtype));
            cudaMemset(tcfLongRangeAsympK.m_data, 0, grid_p->totalLocalPointsK*sizeof(GPUtype));
        }

        if(periodic == false){
            // This is the real-space long range portion of the
            // solute-solvent Ewald potential!

            // Asymptotic values for the DCF in R-space.
            // This is eq. 29 of Kovalenko/Hirata 2000.
            if(treeDCF == false){
                //? Not sure about the reason for this offset in the z direction
                //  Offset for the tested cases are all 0's.
                GPUtype offset_z = grid_p->spacing[2] * grid_p->OffsetR[2];

                // Call the BareTree version of direct sum here :)
                dcf_long_range_asymptotics_R(0, grid_p->localDimsR[0] - 1, 
                                             0, grid_p->localDimsR[1] - 1,
                                             0, grid_p->localDimsR[2] - 1,
                                             0, 0, 0 + offset_z,
                                             grid_p->spacing[0], grid_p->spacing[1], grid_p->spacing[2],
                                             grid_p->localDimsR[0], grid_p->localDimsR[1], grid_p->localDimsR[2],
                                             soluteclass_p->numAtoms, 0,
                                             soluteclass_p->position.m_data, soluteclass_p->position.m_data + soluteclass_p->numAtoms, soluteclass_p->position.m_data + 2*soluteclass_p->numAtoms, soluteclass_p->charge.m_data,
                                             chargeSmear, dcfLongRangeAsympR.m_data);
                
//                 // Writing dcf_long_range to compare with fortran: I will keep this here for now
// #if RISMCUDA_DOUBLE
//                 ofstream output_file("/home/fcarvalho/rism3d.cuda.test.chg/dcf_long_rang_db.txt");
// #else
//                 ofstream output_file("/home/fcarvalho/rism3d.cuda.test.chg/dcf_long_rang_float.txt");
// #endif // RISMCUDA_DOUBLE

//                 if (output_file.is_open()) {
//                     output_file << std::scientific << std::setprecision(16);
//                     for(int i = 0; i < grid_p->localDimsR[0]*grid_p->localDimsR[1]*grid_p->localDimsR[2]; i++){
//                         output_file << dcfLongRangeAsympR.m_data[i] << endl;
//                     }
//                     output_file.close();
//                 }
                
                
            }
            
            // Asymptotic values for  TCF in R-space.
            if(treeTCF == false){
                if(solventclass_p->ionic == true){
                    //? Not sure about the reason for this offset in the z direction
                    //  Offset for the tested cases are all 0's.
                    GPUtype offset_z = grid_p->spacing[2] * grid_p->OffsetR[2];

                    tcf_long_range_asymptotics_R(0, grid_p->localDimsR[0] - 1, 
                                                 0, grid_p->localDimsR[1] - 1,
                                                 0, grid_p->localDimsR[2] - 1,
                                                 0, 0, 0 + offset_z,
                                                 grid_p->spacing[0], grid_p->spacing[1], grid_p->spacing[2],
                                                 grid_p->localDimsR[0], grid_p->localDimsR[1], grid_p->localDimsR[2],
                                                 soluteclass_p->numAtoms, 0,
                                                 soluteclass_p->position.m_data, soluteclass_p->position.m_data + soluteclass_p->numAtoms, soluteclass_p->position.m_data + 2*soluteclass_p->numAtoms, soluteclass_p->charge.m_data,
                                                 solventclass_p->xappa, chargeSmear, solventclass_p->dielconst, tcfLongRangeAsympR.m_data);

//                 // Writing tcf_long_range to compare with fortran: I will keep this here for now
// #if RISMCUDA_DOUBLE
//                 ofstream output_file("/home/fcarvalho/rism3d.cuda.ionic.test.chg/tcf_long_rang_db.txt");
// #else
//                 ofstream output_file("/home/fcarvalho/rism3d.cuda.ionic.test.chg/tcf_long_rang_float.txt");
// #endif // RISMCUDA_DOUBLE

//                 if (output_file.is_open()) {
//                     output_file << std::scientific << std::setprecision(16);
//                     for(int i = 0; i < grid_p->localDimsR[0]*grid_p->localDimsR[1]*grid_p->localDimsR[2]; i++){
//                         output_file << tcfLongRangeAsympR.m_data[i] << endl;
//                     }
//                     output_file.close();
//                 }
                }
            }

            if(treeDCF == true || treeTCF == true){
                cout << "Treecode not available, yet." << endl;
                cout << "Use flag --notreeDCF and --notreeTCF to run direct sum version." << endl;
                abort();
            }

        }
        
        if(periodic == false || strcmp(periodicPotential, "ewald") == 0){
            int ig0;
            if(grid_p->OffsetR[2] == 0){
                ig0 = 1;
            }
            else{
                ig0 = 0;
            }

            /* in the Fortran version, we have the following check:
               waveVectors2(ig) > this%cut2_chlk
               Here, since the dcfLongRangeAsympK values have been already set to
               zero, we just need to do the calculations when 
               waveVectors2(ig) < this%cut2_chlk
               in the cuda kernel
            */

            if(periodic == true){
                cout << "Periodic boundary not supported, yet." << endl;
                abort();
            }
            else{
                calc_dcf_tcf_LongRangeAsympK(grid_p->waveVectorX.m_data,
                             grid_p->waveVectorY.m_data,
                             grid_p->waveVectorZ.m_data,
                             grid_p->waveVectors2.m_data,
                             cut2_chlk,
                             soluteclass_p->position.m_data,
                             soluteclass_p->charge.m_data,
                             dcfLongRangeAsympK.m_data,
                             asympk_const, 
                             smear2_4,
                             soluteclass_p->numAtoms,
                             grid_p->totalLocalPointsK/2,
                             ig0, solventclass_p->ionic, tcfLongRangeAsympK.m_data,
                             xappa2, solventclass_p->dielconst);

//                 // Writing dcf_long_range_K to compare both cuda kernel and serial implementation: I will keep this here for now
// #if RISMCUDA_DOUBLE
//                 ofstream output_file("/home/fcarvalho/rism3d.cuda.test.chg/dcfLongRangeAsympK_db.txt");
// #else
//                 ofstream output_file("/home/fcarvalho/rism3d.cuda.test.chg/dcfLongRangeAsympK_float.txt");
// #endif // RISMCUDA_DOUBLE
//                 if (output_file.is_open()) {
//                     output_file << std::scientific << std::setprecision(16);
//                     for(int i = 0; i < grid_p->totalLocalPointsK; i++){
//                         output_file << dcfLongRangeAsympK.m_data[i] << endl;
//                     }
//                     output_file.close();
//                 }
            }

            // Getting the difference between long-range asymptotic functions
            // of the TCF at k = 0.
            // Here we are still working if huvk0 in row major order allocated as
            // huvk0.alloc_mem(2, solventclass_p->numAtomTypes);
            if(periodic == false){
                // Calculate this only if offset in the z direction is zero
                if(grid_p->OffsetR[2] == 0){
                    // Allocate memory for sumcos_0 and sumsin_0: these will be needed for calculating 
                    // huvk0. I am storing these values so we do not need to do another kernel launch
                    // to calculate them, since these are calculated while running calc_dcfLongRangeAsympK
                    cudaMallocManaged((void**)&sumcos_0, 1*sizeof(GPUtype));
                    cudaMallocManaged((void**)&sumsin_0, 1*sizeof(GPUtype));

                    // Here I am calculating the sumcos and sumsin values in a cuda kernel
                    calc_sum_cos_sin_huvk0(grid_p->waveVectorX.m_data[0],
                             grid_p->waveVectorY.m_data[0],
                             grid_p->waveVectorZ.m_data[0],
                             soluteclass_p->position.m_data,
                             soluteclass_p->charge.m_data,
                             soluteclass_p->numAtoms,
                             sumcos_0, sumsin_0);
                    
                    // Calculating huvk0 is just a loop over the solvent atom types,
                    // so it seems that having a cuda kernel for it would not be optimal.
                    // For now I will calculate it on the host side, but later we can make it better
                    // Maybe calculating the long range DCF and huvk0 in the same kernel would be a good
                    // option
                    for(int i = 0; i < solventclass_p->numAtomTypes; i++){
                        // Calculate huvk0
                        // Real
                        huvk0.m_data[i*2 + 0] = solventclass_p->delhv0.m_data[i] * (*sumcos_0) / grid_p->boxVolume;

                        // Imaginary
                        huvk0.m_data[i*2 + 1] = solventclass_p->delhv0.m_data[i] * (*sumsin_0) / grid_p->boxVolume;

                        //Need to do the same for huvk0_dT
                    }

                    cudaFree(sumcos_0);
                    cudaFree(sumsin_0);
                }
            }
            else{
                cout << "Periodic boundary not supported, yet." << endl;
                abort();
            } 
        }
    }

    void rism3d_potential :: int_h(GPUtype* h){
        GPUtype solute_total_charge;

        solute_total_charge = cu_reduce_sum(soluteclass_p->charge.m_data, soluteclass_p->numAtoms);

        for(int iv = 0; iv < solventclass_p->numAtomTypes; iv++){
            if(solventclass_p->charge_sp.m_data[iv] != 0){
                h[iv] = -4 * PI / (solventclass_p->dielconst * solventclass_p->xappa * solventclass_p->xappa) * solventclass_p->charge_sp.m_data[iv] * solute_total_charge;
            }
            else{
                h[iv] = 0;
            }
        }

    }

    void rism3d_potential :: int_h2_hc(GPUtype* h2, GPUtype* hc){
        // If we have calculated h2 and hc before:
        if(cached_h2.m_data[0] != MAX_VAL && cached_hc.m_data[0] != MAX_VAL){
            for(int iv = 0; iv < solventclass_p->numAtomTypes; iv++){
                h2[iv] = cached_h2.m_data[iv];
                hc[iv] = cached_hc.m_data[iv];
            }
            return;
        }

        // Checking charges:
        if(soluteclass_p->charged == false || solventclass_p->ionic == false){
            for(int iv = 0; iv < solventclass_p->numAtomTypes; iv++){
                h2[iv] = 0.0;
                hc[iv] = 0.0;
            }
            return;
        }

        int Nintmx = 200;
        GPUtype* argum = new GPUtype[Nintmx];
        GPUtype* weight = new GPUtype[Nintmx];

        // Initialize Gauss-Legendre integration.
        gaussquad_legendre(0, 1, argum, weight, Nintmx);
        GPUtype sumhc = 0;
        GPUtype sumh2 = 0;

        int n0 = 0;
        int nf = Nintmx;

        // Implementing serial version first.
        // Part of the following code could be translated into a
        // cuda Kernel: will be the next step after finishing thermodynamic
        // calculations
        for(int ik = n0; ik < nf; ik++){
            GPUtype k = argum[ik] / (1 - argum[ik]);

            // Bessel part
            GPUtype sumb = 0;

            for(int i = 1; i < soluteclass_p->numAtoms; i++){
                for(int j = 0; j < i; j++){
                    GPUtype dx = soluteclass_p->position.m_data[i] - soluteclass_p->position.m_data[j];
                    GPUtype dy = soluteclass_p->position.m_data[i + soluteclass_p->numAtoms] - soluteclass_p->position.m_data[j + soluteclass_p->numAtoms];
                    GPUtype dz = soluteclass_p->position.m_data[i + 2*soluteclass_p->numAtoms] - soluteclass_p->position.m_data[j + 2*soluteclass_p->numAtoms];

                    GPUtype r2 = dx * dx + dy * dy + dz * dz;

                    GPUtype xarg = solventclass_p->xappa * k * sqrt(r2);

                    if(xarg == 0){
                        sumb = sumb + 1 * soluteclass_p->charge.m_data[i] * soluteclass_p->charge.m_data[j];
                    }
                    else{
                        sumb = sumb + sin(xarg) / xarg * soluteclass_p->charge.m_data[i] * soluteclass_p->charge.m_data[j];
                    }
                }
            }

            sumb = sumb * 2;

            for(int i = 0; i < soluteclass_p->numAtoms; i++){
                sumb = sumb + 1 * soluteclass_p->charge.m_data[i] * soluteclass_p->charge.m_data[i];
            }

            // End of Bessel part.

            GPUtype x2arg = (k * solventclass_p->xappa * chargeSmear) * (k * solventclass_p->xappa * chargeSmear);
            GPUtype k2 = k*k;

            GPUtype denom = argum[ik]*argum[ik] + (1 - argum[ik])*(1 - argum[ik]);

            sumhc = sumhc - exp(-x2arg / 2) / denom * sumb * weight[ik];
            
            sumh2 = sumh2 + exp(-x2arg / 2) * argum[ik]*argum[ik] / (denom*denom) * sumb * weight[ik];
        }

        // End of ik loop

        // Site xappa.
        for(int iv = 0; iv < solventclass_p->numAtomTypes; iv++){
            h2[iv] = 8 * PI / solventclass_p->dielconst * solventclass_p->charge_sp.m_data[iv] * solventclass_p->charge_sp.m_data[iv];
            hc[iv] = 8 * PI / solventclass_p->dielconst * solventclass_p->charge.m_data[iv] * solventclass_p->charge_sp.m_data[iv];
        }

        // Divided by total xappa.
        for(int iv = 0; iv < solventclass_p->numAtomTypes; iv++){
            if(solventclass_p->xappa != 0){
                h2[iv] = (h2[iv] / solventclass_p->xappa) * sumh2 / (PI * solventclass_p->dielconst);
                hc[iv] = (-hc[iv] / solventclass_p->xappa) * sumhc  / PI;
            }
            else{
                h2[iv] = 0;
                hc[iv] = 0;
            }
        }

        cached_h2.m_data = h2;
        cached_hc.m_data = hc;

        delete[] argum;
        delete[] weight;
    }

    void rism3d_potential :: setcut_asympKTolerance(GPUtype asympKSpaceTolerance, GPUtype boxVolume){
        if(asympKSpaceTolerance == 0){
#if defined(RISMCUDA_DOUBLE)
            cut2_chlk = DBL_MAX;
#else
            cut2_chlk = FLT_MAX;
#endif //RISMCUDA_DOUBLE
        }
        else{
            cut2_chlk = asympck_cut_calc(soluteclass_p->charge.m_data, soluteclass_p->numAtoms, 
                                         solventclass_p->charge.m_data, solventclass_p->numAtomTypes, 
                                         chargeSmear, boxVolume, 
                                         asympKSpaceTolerance);
            if (std::isnan(cut2_chlk)) {
                cout << "Could not converge k-space asymptotics cutoff. 'Try using a smaller error tolerance or no cutoff." << endl;
            }
        }
    }

    GPUtype rism3d_potential :: asympck_cut_calc(GPUtype *soluteCharge, int solute_size,
                                                 GPUtype *solventCharge, int solvent_size,
                                                 GPUtype chargeSmear, GPUtype boxVolume, 
                                                 GPUtype asympKSpaceTolerance){
        DCFTolerance = asympKSpaceTolerance;
        DCFCoefficient = get_maxval_abs(soluteCharge, solute_size) * 4*PI*sqrt(2.0)/boxVolume * get_maxval_abs(solventCharge, solvent_size);
        smearSq = chargeSmear*chargeSmear;

        GPUtype cutoff = root_newton([this](GPUtype k2) { return this->asympck_cut(k2); },
                                     [this](GPUtype k2) { return this->asympck_cut_deriv(k2); }, 
                                     1.0, 1e-16);

        return cutoff;
    }

    GPUtype rism3d_potential :: asympck_cut(GPUtype k2){
        GPUtype asympck_cut = DCFCoefficient/k2*exp(-0.25*k2*smearSq) - DCFTolerance;
        return asympck_cut;
    }

    GPUtype rism3d_potential :: asympck_cut_deriv(GPUtype k2){
        GPUtype asympck_cut_deriv = -DCFCoefficient/(k2*k2)*exp(-0.25*k2*smearSq)*(1+0.25*smearSq*k2);
        return asympck_cut_deriv;
    }

}