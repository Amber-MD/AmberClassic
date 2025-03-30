#include <iostream>
#include "rism3d_closure_kh.hpp"
#include <sstream>
#include <iostream>
using namespace std;

namespace rism3d_c {

    kh :: kh(){
        // cout << "Default object from rism3d_closure_kh is being created!" << endl;
    }

    kh :: kh(rism3d_solvent_class *solventclass, rism3d_solute_class *soluteclass, rism3d_grid *grid, rism3d_potential *pot, rism3d_safemem *memalloc){
        // cout << "Object from rism3d_closure_kh is being created!" << endl;
        solv = solventclass;
        solu = soluteclass;
        grid_p = grid;
        memalloc_p = memalloc;
        pot_p = pot;

        heaviside.set_memalloc(memalloc, true);
        heaviside.alloc_mem(solv->numAtomTypes, grid_p->globalDimsK[0], grid_p->globalDimsK[1], grid_p->globalDimsK[2]);

        mu_grid.set_memalloc(memalloc, true);
        mu_grid.alloc_mem(solv->numAtomTypes, grid_p->globalDimsR[0], grid_p->globalDimsR[1], grid_p->globalDimsR[2]);

        solvPotEne_grid.set_memalloc(memalloc, true);
        solvPotEne_grid.alloc_mem(solv->numAtomTypes, grid_p->globalDimsR[0], grid_p->globalDimsR[1], grid_p->globalDimsR[2]);

        cuv_int.set_memalloc(memalloc, true);
        cuv_int.alloc_mem(solv->numAtomTypes);

        lj_corr.set_memalloc(memalloc, true);
        lj_corr.alloc_mem(solv->numAtomTypes);

        excChemPot_1d = new GPUtype[solv->numAtomTypes];
        excChemPot_1d_db = new double[solv->numAtomTypes];

        solvPotEne_1d = new GPUtype[solv->numAtomTypes];
        solvPotEne_1d_db = new double[solv->numAtomTypes];

        cuv_int_dp = new double[solv->numAtomTypes];

        kirkwoodBuff_1d.set_memalloc(memalloc);
        // kirkwoodBuff_1d.alloc_mem(solv->numAtomTypes);

        kirkwoodBuff_1d_db.set_memalloc(memalloc);
        // kirkwoodBuff_1d_db.alloc_mem(solv->numAtomTypes);

        excessPart_1d.set_memalloc(memalloc);
        excessPart_1d.alloc_mem(solv->numAtomTypes);

        excessPart_1d_db.set_memalloc(memalloc);
        excessPart_1d_db.alloc_mem(solv->numAtomTypes);

        excessPartBox_1d.set_memalloc(memalloc);
        excessPartBox_1d.alloc_mem(solv->numAtomTypes);

        excessPartBox_1d_db.set_memalloc(memalloc);
        excessPartBox_1d_db.alloc_mem(solv->numAtomTypes);

        h_val.set_memalloc(memalloc);

        excessChemicalPotentialh2lr.set_memalloc(memalloc);
        excessChemicalPotentialhclr.set_memalloc(memalloc);

    }

    kh :: ~kh(){
        // cout << "Object from rism3d_closure_kh is being deleted!" << endl;
        delete[] excChemPot_1d;
        delete[] excChemPot_1d_db;
        delete[] solvPotEne_1d;
        delete[] solvPotEne_1d_db;
        delete[] cuv_int_dp;

    }

    void kh :: guv(GPUtype *uuv, GPUtype *guv, GPUtype *huv, GPUtype *cuv, int size){
        cu_guv(uuv, guv, huv, cuv, size);
    }    

    double* kh :: excessChemicalPotential(GPUtype *huv, GPUtype *cuv, int size){
        int size_real = grid_p->globalDimsR[0] * grid_p->globalDimsR[1] * grid_p->globalDimsR[2];
        
        cu_set_heaviside(heaviside.m_data, huv, solv->numAtomTypes * size);
        cudaDeviceSynchronize();

        for(int i = 0; i < solv->numAtomTypes; i++){
            cu_excessChemicalPotential_grid(heaviside.m_data + i * size, 
                            huv + i * size, 
                            cuv + i * size,
                            mu_grid.m_data + i * size_real,
                            grid_p->globalDimsR[0], 
                            grid_p->globalDimsR[1], 
                            grid_p->globalDimsR[2],
                            solv->density.m_data[i],
                            grid_p->voxelVolume);      

            excChemPot_1d[i] = 0;
        }
        cudaDeviceSynchronize();

        for(int i = 0; i < solv->numAtomTypes; i++){
            excChemPot_1d[i] = cu_reduce_sum(mu_grid.m_data + i * size_real, size_real);
        }
        cudaDeviceSynchronize();


#if defined(RISMCUDA_DOUBLE)
        return excChemPot_1d;
#else
        for(int i = 0; i < solv->numAtomTypes; i++){
            excChemPot_1d_db[i] = static_cast<double>(excChemPot_1d[i]);
        }
        return excChemPot_1d_db;
#endif //RISMCUDA_DOUBLE
    }

    double* kh :: aexcessChemicalPotential(GPUtype *huv, GPUtype *cuv, int size, bool long_range_corr){
        /*
        I have looked at the fortran version and for the case where the solute.charged is true,
        but solvent.ionic is false, the long range corrections for the excess chemical potential 
        are zeroed out, because huvlr = 0, excessChemicalPotentialh2lr = 0 
        and excessChemicalPotentialhclr(iv) = 0.
        */
        if(long_range_corr == true and solv->ionic == true and solu->charged == true){
            cu_set_heaviside(heaviside.m_data, huv, solv->numAtomTypes * size);
            cudaDeviceSynchronize();

            int size_real = grid_p->globalDimsR[0] * grid_p->globalDimsR[1] * grid_p->globalDimsR[2];

            // Working with excessChemicalPotentialh2lr and excessChemicalPotentialhclr as class members
            // avoided the issue of getting two different values for the excess chemical potential.
            // Since these arrays has the size of solvent atom types, it will not take too much memory
            // and should be fine leaving as it is.
            if(excessChemicalPotentialh2lr.m_data == nullptr || excessChemicalPotentialhclr.m_data == nullptr){
                excessChemicalPotentialh2lr.alloc_mem(solv->numAtomTypes);
                excessChemicalPotentialhclr.alloc_mem(solv->numAtomTypes);

                // Long range part.
                pot_p->int_h2_hc(excessChemicalPotentialh2lr.m_data, excessChemicalPotentialhclr.m_data);
            }

            for(int iv = 0; iv < solv->numAtomTypes; iv++){
                if(solu->totalCharge * solv->charge_sp.m_data[iv] > 0){
                    excessChemicalPotentialh2lr.m_data[iv] = 0.0;
                }
            }

            for(int iv = 0; iv < solv->numAtomTypes; iv++){
                // Getting excess potential grid with long range correction
                cu_excessChemicalPotential_lr_grid(heaviside.m_data + iv * size, 
                                                   huv + iv * size, 
                                                   cuv + iv * size,
                                                   mu_grid.m_data + iv * size_real,
                                                   grid_p->globalDimsR[0], 
                                                   grid_p->globalDimsR[1], 
                                                   grid_p->globalDimsR[2],
                                                   solv->density.m_data[iv],
                                                   grid_p->voxelVolume,
                                                   solv->charge.m_data[iv],
                                                   pot_p->dcfLongRangeAsympR.m_data,
                                                   solu->totalCharge,
                                                   solv->charge_sp.m_data[iv],
                                                   pot_p->tcfLongRangeAsympR.m_data);

                // initializing excess chemical potential array to return
                // at the end of this function
                excChemPot_1d[iv] = 0;
            }
            cudaDeviceSynchronize();

            for(int iv = 0; iv < solv->numAtomTypes; iv++){
                excChemPot_1d[iv] = cu_reduce_sum(mu_grid.m_data + iv * size_real, size_real) + solv->density.m_data[iv] * (excessChemicalPotentialh2lr.m_data[iv] - excessChemicalPotentialhclr.m_data[iv])/2;
            }
            cudaDeviceSynchronize();

#if defined(RISMCUDA_DOUBLE)
            return excChemPot_1d;
#else
            for(int i = 0; i < solv->numAtomTypes; i++){
                excChemPot_1d_db[i] = static_cast<double>(excChemPot_1d[i]);
            }
            return excChemPot_1d_db;
#endif //RISMCUDA_DOUBLE
        }
        else{
            return excessChemicalPotential(huv, cuv, size);
        }
        
    }

    double* kh :: solvPotEne(GPUtype *guv, GPUtype *uuv, int size){
        
        int size_real = grid_p->globalDimsR[0] * grid_p->globalDimsR[1] * grid_p->globalDimsR[2];

        for(int i = 0; i < solv->numAtomTypes; i++){
            cu_solvPotEne_grid(guv + i * size,
                               uuv + i * size,
                               solvPotEne_grid.m_data + i * size_real,
                               grid_p->globalDimsR[0], 
                               grid_p->globalDimsR[1], 
                               grid_p->globalDimsR[2],
                               solv->density.m_data[i],
                               grid_p->voxelVolume
                               );
        }
        cudaDeviceSynchronize();

        for(int i = 0; i < solv->numAtomTypes; i++){
            solvPotEne_1d[i] = cu_reduce_sum(solvPotEne_grid.m_data + i * size_real, size_real);
        }
        cudaDeviceSynchronize();

#if defined(RISMCUDA_DOUBLE)
            return solvPotEne_1d;
#else
            for(int i = 0; i < solv->numAtomTypes; i++){
                solvPotEne_1d_db[i] = static_cast<double>(solvPotEne_1d[i]);
            }
            return solvPotEne_1d_db;
#endif //RISMCUDA_DOUBLE
    }

    double kh :: partialMolarVolume(GPUtype* cuv){
        partMolVol = 0;
        DCFintegral(cuv);
#if defined(RISMCUDA_DOUBLE)
        for(int i = 0; i < solv->numAtomTypes; i++){
            partMolVol = partMolVol + solv->density.m_data[i] * cuv_int.m_data[i];
        }

        partMolVol = solv->xikt * (1.0 - partMolVol);

        return partMolVol;
#else
        for(int i = 0; i < solv->numAtomTypes; i++){
            partMolVol = partMolVol + solv->density.m_data[i] * cuv_int.m_data[i];
        }

        partMolVol = solv->xikt * (1.0 - partMolVol);

        partMolVol_db = static_cast<double>(partMolVol);

        return partMolVol_db;
#endif //RISMCUDA_DOUBLE
    }

    void kh :: DCFintegral(GPUtype* cuv){
        int size = grid_p->globalDimsK[0] * grid_p->globalDimsK[1] * grid_p->globalDimsK[2];

        LJCorrection_DCF_int();

        for(int i = 0; i < solv->numAtomTypes; i++){
            cuv_int.m_data[i] = cu_reduce_sum(cuv + i * size, size) * grid_p->voxelVolume + lj_corr.m_data[i];
        }
    }

    bool kh :: check_dbl_flt_max(GPUtype* ljCutoffs2, int solute_atoms, int solvent_atoms){
        for(int i = 0; i < solute_atoms * solvent_atoms; i++){
#if defined(RISMCUDA_DOUBLE)
            if (pot_p->ljCutoffs2.m_data[i] != DBL_MAX) {
                return false;  // Found a value that's not DBL_MAX
            }
            return true;
#else
            if (pot_p->ljCutoffs2.m_data[i] != FLT_MAX) {
                return false;  // Found a value that's not DBL_MAX
            }
            return true;
#endif //RISMCUDA_DOUBLE
        }
    }

    void kh :: LJCorrection_DCF_int(){
        if(check_dbl_flt_max(pot_p->ljCutoffs2.m_data, solu->numAtoms, solv->numAtomTypes) || pot_p->applyLJCorrection == false){
            for(int i = 0; i < solv->numAtomTypes; i++){
                lj_corr.m_data[i] = 0.0;
            }
        }
        else{
            cout << "in LJCorrection_DCF_int: LJ cutoff is not supported, yet." << endl;
            abort();
        }
    }

    double* kh :: get_cuv_int(){
#if defined(RISMCUDA_DOUBLE)
        return cuv_int.m_data;
#else
        for(int i = 0; i < solv->numAtomTypes; i++){
            cuv_int_dp[i] = static_cast<double>(cuv_int.m_data[i]);
        }
        return cuv_int_dp;
#endif //RISMCUDA_DOUBLE
    }
    
    void kh :: get_KB(GPUtype* huv, bool long_range_corr){
        if(long_range_corr == true and solv->ionic == true and solu->charged == true){
            int Nx = grid_p->localDimsR[0];
            int Ny = grid_p->localDimsR[1];
            int Nz = grid_p->localDimsR[2];

            int size = Nx * Ny * (Nz + 2);
            
            GPUtype corr = LJCorrection_TCF_int();

            if(h_val.m_data == NULL){
                h_val.alloc_mem(solv->numAtomTypes);
                pot_p->int_h(h_val.m_data);
            }

            if(kirkwoodBuff_1d.m_data == nullptr){
                kirkwoodBuff_1d.alloc_mem(solv->numAtomTypes);
                kirkwoodBuff_1d_db.alloc_mem(solv->numAtomTypes);

                for(int iv = 0; iv <  solv->numAtomTypes; iv++){
                    kirkwoodBuff_1d.m_data[iv] = cu_akirkwoodBuff(huv + iv * size, solv->charge_sp.m_data[iv], pot_p->tcfLongRangeAsympR.m_data, Nx, Ny, Nz) * grid_p->voxelVolume + corr + h_val.m_data[iv];
                }

            }
        }
        else{
            int size = grid_p->globalDimsK[0] * grid_p->globalDimsK[1] * grid_p->globalDimsK[2];
            int size_real = grid_p->globalDimsR[0] * grid_p->globalDimsR[1] * grid_p->globalDimsR[2];

            // This commneted block is for checking padded values. I will keep it here for now
            // cout << "Checking padded values: " << endl;
            // for(int i = 0; i < solv->numAtomTypes*size; i++){
            //     huv[i] = huv[i] - 1.0;
            // }
            // for(int idx = 0; idx < 10; idx++){
            //     for(int idy = 0; idy < 10; idy++){
            //         for(int idz = 0; idz < 2; idz++){
            //             int padded_idx = idx * grid_p->globalDimsR[1] * (grid_p->globalDimsR[2]+2) + idy * (grid_p->globalDimsR[2]+2) + grid_p->globalDimsR[2] + idz;
            //             cout << huv[padded_idx] << endl;    
            //         }
            //     }
            // }

            GPUtype corr = LJCorrection_TCF_int();

            if(kirkwoodBuff_1d.m_data == nullptr){
                kirkwoodBuff_1d.alloc_mem(solv->numAtomTypes);
                kirkwoodBuff_1d_db.alloc_mem(solv->numAtomTypes);

                for(int i = 0; i < solv->numAtomTypes; i++){
                    kirkwoodBuff_1d.m_data[i] = cu_reduce_sum(huv + i * size, size) * grid_p->voxelVolume + corr;    
                }
            }
        }
    }

    double* kh :: kirkwoodBuff(GPUtype* huv, bool long_range_corr){
        get_KB(huv,long_range_corr);
#if defined(RISMCUDA_DOUBLE)
        return kirkwoodBuff_1d.m_data;
#else
        for(int i = 0; i < solv->numAtomTypes; i++){
            kirkwoodBuff_1d_db.m_data[i] = static_cast<double>(kirkwoodBuff_1d.m_data[i]);
        }
        return kirkwoodBuff_1d_db.m_data;
#endif //RISMCUDA_DOUBLE       
    }

    GPUtype kh :: LJCorrection_TCF_int(){
        GPUtype correction = 0.0;

        if(check_dbl_flt_max(pot_p->ljCutoffs2.m_data, solu->numAtoms, solv->numAtomTypes) || pot_p->applyLJCorrection == false){
            return correction;
        }

        LJCorrection_DCF_int();

        cout << "Lennard-Jones cutoff not supported, yet" << endl;
        abort();
    }


    /*
    For excessParticles and aexcessParticles functions below I am working with two different objects: excessPart_1d
    and excessPartBox_1d. That is because on the C++ side, the function is returning a pointer. Thus, if we use the same 
    object for the cases where long_range_correct is true and false, as in these lines on amber_rism_interface.F90:
        rismthermo%excessParticlesBox => rism_3d%excessParticles(.false.)
        rismthermo%excessParticles => rism_3d%excessParticles(.true.)
    we will be pointing to the same address and the results for both variables will be the same.

    It may be possible to proceed in a different way by dealing with Shroud again. However,
    usually the number of solvent atom types is usually small and having these additional objects should
    not impact the memory usage.
    */
    double* kh :: excessParticles(GPUtype* huv){
        int size = grid_p->globalDimsK[0] * grid_p->globalDimsK[1] * grid_p->globalDimsK[2];

        for(int i = 0; i < solv->numAtomTypes; i++){
            excessPartBox_1d.m_data[i] = cu_reduce_sum(huv + i * size, size) * grid_p->voxelVolume * solv->density.m_data[i];
            // excessPart_1d.m_data[i] = KB[i] * solv->density.m_data[i];

        }

#if defined(RISMCUDA_DOUBLE)
        return excessPartBox_1d.m_data;
#else
        for(int i = 0; i < solv->numAtomTypes; i++){
            excessPartBox_1d_db.m_data[i] = static_cast<double>(excessPartBox_1d.m_data[i]);
        }
        return excessPartBox_1d_db.m_data;
#endif //RISMCUDA_DOUBLE         
    }

    double* kh :: aexcessParticles(GPUtype* huv){
        get_KB(huv,true);

        for(int iv = 0; iv < solv->numAtomTypes; iv++){
            excessPart_1d.m_data[iv] = kirkwoodBuff_1d.m_data[iv] * solv->density.m_data[iv];
        }

#if defined(RISMCUDA_DOUBLE)
        return excessPart_1d.m_data;
#else
        for(int i = 0; i < solv->numAtomTypes; i++){
            excessPart_1d_db.m_data[i] = static_cast<double>(excessPart_1d.m_data[i]);
        }
        return excessPart_1d_db.m_data;
#endif //RISMCUDA_DOUBLE         
    }

}