#include <iostream>
#include "rism3d_closure.hpp"
#include <sstream>
#include <iostream>
using namespace std;

namespace rism3d_c {

    rism3d_closure :: rism3d_closure(){
        // cout << "Default object from rism3d_closure is being created!" << endl;
    }

    rism3d_closure :: rism3d_closure(rism3d_solvent_class *solventclass, rism3d_solute_class *soluteclass, rism3d_grid *grid, rism3d_potential *pot, rism3d_safemem *memalloc){
        // cout << "Object from rism3d_closure is being created!" << endl;
        solv = solventclass;
        solu = soluteclass;
        grid_p = grid;
        memalloc_p = memalloc;
        pot_p = pot;
        
        kirkwoodBuff_1d = new GPUtype[solv->numAtomTypes];
        kirkwoodBuff_1d_db = new double[solv->numAtomTypes];
    }

    rism3d_closure :: ~rism3d_closure(){
        // cout << "Object from rism3d_closure is being deleted!" << endl;
        delete[] kirkwoodBuff_1d;
        delete[] kirkwoodBuff_1d_db;
    }

    void rism3d_closure :: guv(GPUtype *uuv, GPUtype *guv, GPUtype *huv, GPUtype *cuv, int size){
        cout << "guv in closure parent class" << endl;
        cout << "This should not have been called from here" << endl;
        abort();
    }  

    double* rism3d_closure :: excessChemicalPotential(GPUtype *huv, GPUtype *cuv, int size){
        cout << "excessChemicalPotential in closure parent class" << endl;
        cout << "This should not have been called from here" << endl;
        abort();
    }

    double* rism3d_closure :: aexcessChemicalPotential(GPUtype *huv, GPUtype *cuv, int size, bool long_range_corr){
        cout << "aexcessChemicalPotential in closure parent class" << endl;
        cout << "This should not have been called from here" << endl;
        abort();
    }

    double* rism3d_closure :: rism3d_closure_aexcessChemicalPotentialGF(){
        excChemPotGF_1d = new double[solv->numAtomTypes];
        for(int i = 0; i < solv->numAtomTypes; i++){
            excChemPotGF_1d[i] = 12.34 + i;
        }
        return excChemPotGF_1d;
    }

    double* rism3d_closure :: solvPotEne(GPUtype *guv, GPUtype *uuv, int size){
        cout << "solvPotEne in closure parent class" << endl;
        cout << "This should not have been called from here" << endl;
        abort();
    }

    double rism3d_closure :: partialMolarVolume(GPUtype* cuv){
        cout << "partialMolarVolume in closure parent class" << endl;
        cout << "This should not have been called from here" << endl;
        abort();
    }

    double* rism3d_closure :: excessParticles(GPUtype* huv){
        cout << "excessParticles in closure parent class" << endl;
        cout << "This should not have been called from here" << endl;
        abort();
    }

    double* rism3d_closure :: aexcessParticles(GPUtype* huv){
        cout << "aexcessParticles in closure parent class" << endl;
        cout << "This should not have been called from here" << endl;
        abort();
    }

    double* rism3d_closure :: kirkwoodBuff(GPUtype* huv, bool long_range_corr){
        cout << "kirkwoodBuff in closure parent class" << endl;
        cout << "This should not have been called from here" << endl;
        abort();
    }

    double* rism3d_closure :: rism3d_closure_DCFintegral(){
        DCFintegral_1d = new double[solv->numAtomTypes];
        for(int i = 0; i < solv->numAtomTypes; i++){
            DCFintegral_1d[i] = 3.13 + i;
        }
        return DCFintegral_1d;
    }

    double* rism3d_closure :: rism3d_closure_solvationEnergy(){
        solvationEnergy_1d = new double[solv->numAtomTypes];
        for(int i = 0; i < solv->numAtomTypes; i++){
            solvationEnergy_1d[i] = 1.31 + i;
        }
        return solvationEnergy_1d;
    }

    double* rism3d_closure :: rism3d_closure_solvationEnergyGF(){
        solvationEnergyGF_1d = new double[solv->numAtomTypes];
        for(int i = 0; i < solv->numAtomTypes; i++){
            solvationEnergyGF_1d[i] = 4.1 + i;
        }
        return solvationEnergyGF_1d;
    }

    double* rism3d_closure :: rism3d_closure_solvationEnergyGF_site_map(rism3d_grid* grid,
                                                            double* huv, double* huv_dT, double* cuv, double* cuv_dT){
        solvationEnergyGF_site_map_1d = new double[grid->localDimsR[0]*grid->localDimsR[1]*grid->localDimsR[2]*solv->numAtomTypes];
        int count = 616;
        for(int i = 0; i < grid->localDimsR[0]; i++){
            for(int j = 0; j < grid->localDimsR[1]; j++){
                for(int k = 0; k < grid->localDimsR[2]; k++){
                    for(int q = 0; q < solv->numAtomTypes; q++){
                        solvationEnergyGF_site_map_1d[i + grid->localDimsR[0]*j + 
                        grid->localDimsR[0]*grid->localDimsR[1]*k + grid->localDimsR[0]*grid->localDimsR[1]*grid->localDimsR[2]*q] = count;
                        count++;
                    }
                }
            }
        }
        return solvationEnergyGF_site_map_1d;
    }

    double* rism3d_closure :: rism3d_closure_excessParticles_dT(){
        excessParticles_dT_1d = new double[solv->numAtomTypes];
        for(int i = 0; i < solv->numAtomTypes; i++){
            excessParticles_dT_1d[i] = 4.32 + i;
        }
        return excessParticles_dT_1d;
    }

    double* rism3d_closure :: rism3d_closure_kirkwoodBuff_dT(){
        kirkwoodBuff_dT_1d = new double[solv->numAtomTypes];
        for(int i = 0; i < solv->numAtomTypes; i++){
            kirkwoodBuff_dT_1d[i] = 65.4 + i;
        }
        return kirkwoodBuff_dT_1d;
    }

    double* rism3d_closure :: rism3d_closure_DCFintegral_dT(){
        DCFintegral_dT_1d = new double[solv->numAtomTypes];
        for(int i = 0; i < solv->numAtomTypes; i++){
            DCFintegral_dT_1d[i] = 45.6 + i;
        }
        return DCFintegral_dT_1d;
    }

    double* rism3d_closure :: get_cuv_int(){
        cout << "get_cuv_int in closure parent class" << endl;
        cout << "This should not have been called from here" << endl;
        abort();
    }

}