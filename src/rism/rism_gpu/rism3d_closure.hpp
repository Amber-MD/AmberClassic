#ifndef RISM3DCLOSURE_HPP
#define RISM3DCLOSURE_HPP

#include <iostream>
#include "rism3d_solvent_class.hpp"
#include "rism3d_solute_class.hpp"
#include "rism3d_grid.hpp"
#include "rism3d_potential.hpp"
#include "rism3d_safemem.hpp"
#include "varTypes.hpp"
using namespace std;

namespace rism3d_c {

    class rism3d_closure{
        protected:
            rism3d_solvent_class *solv;
            rism3d_solute_class *solu;
            rism3d_potential *pot_p;
            rism3d_grid *grid_p;
            rism3d_safemem *memalloc_p;
            
        public:
            rism3d_closure();
            rism3d_closure(rism3d_solvent_class *solventclass, rism3d_solute_class *soluteclass, 
                           rism3d_grid *grid, rism3d_potential *pot, rism3d_safemem *memalloc);
            virtual ~rism3d_closure();

            virtual void guv(GPUtype *uuv, GPUtype *guv, GPUtype *huv, GPUtype *cuv, int size);

            // int excChemPot_size;
            // GPUtype* excChemPot_1d;
            virtual double* excessChemicalPotential(GPUtype *huv, GPUtype *cuv, int size);
            virtual double* aexcessChemicalPotential(GPUtype *huv, GPUtype *cuv, int size, bool long_range_corr);

            // int excChemPotGF_size;
            double* excChemPotGF_1d;
            double* rism3d_closure_aexcessChemicalPotentialGF();

            // int solvPotEne_size;
            double* solvPotEne_1d;
            virtual double* solvPotEne(GPUtype *guv, GPUtype *uuv, int size);

            // int excessPart_size;
            double* excessPart_1d;
            virtual double* excessParticles(GPUtype* huv);

            virtual double* aexcessParticles(GPUtype* huv);

            // double partMolVol;
            virtual double partialMolarVolume(GPUtype* cuv);

            // int kirkwoodBuff_size;
            GPUtype* kirkwoodBuff_1d;
            double* kirkwoodBuff_1d_db;
            virtual double* kirkwoodBuff(GPUtype* huv, bool long_range_corr);

            // I am leaving "rism3d_closure_DCFintegral" commented out here, 
            // but right now the DCF integral is being carried out inside 
            // partialMolarVolume by calling DCFintegral function implemented 
            // in the child closure class.
            // We may want a DCFintegral function to be called from outside the 
            // closure class later, but right now it is not necessary
            // int DCFintegral_size;
            double* DCFintegral_1d;
            double* rism3d_closure_DCFintegral();
            
            // int solvationEnergy_size;
            double* solvationEnergy_1d;
            double* rism3d_closure_solvationEnergy();
            
            // int solvationEnergyGF_size;
            double* solvationEnergyGF_1d;
            double* rism3d_closure_solvationEnergyGF();

            double* solvationEnergyGF_site_map_1d;
            double* rism3d_closure_solvationEnergyGF_site_map(rism3d_grid* grid,
                                                    double* huv, double* huv_dT, double* cuv, double* cuv_dT);

            // int excessParticles_dT_size;
            double* excessParticles_dT_1d;
            double* rism3d_closure_excessParticles_dT();

            // int kirkwoodBuff_dT_size;
            double* kirkwoodBuff_dT_1d;
            double* rism3d_closure_kirkwoodBuff_dT();

            // int DCFintegral_dT_size;
            double* DCFintegral_dT_1d;
            double* rism3d_closure_DCFintegral_dT();

            virtual double* get_cuv_int();
    };

}

#endif // RISM3DCLOSURE_HPP