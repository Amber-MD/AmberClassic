#include <iostream>
#include "rism3d_solvent_class.hpp"
using namespace std;

namespace rism3d_c {

    rism3d_solvent_class :: rism3d_solvent_class() : 
    xvv(array_class<GPUtype>::ROW_MAJOR),
    coord(array_class<GPUtype>::COLUMN_MAJOR)
    {
        // cout << "solvent class created!" << endl;
    }

    rism3d_solvent_class :: rism3d_solvent_class(rism3d_safemem *memalloc, rism3d_c::solvent_cpp *solv_f) : 
    // xvv must be ROW major
    xvv(array_class<GPUtype>::ROW_MAJOR),
    coord(array_class<GPUtype>::COLUMN_MAJOR)
    {
        // cout << "solvent class created!" << endl;

        memalloc_p = memalloc;

        set_solvent(solv_f);

    }

    rism3d_solvent_class :: ~rism3d_solvent_class(){
        // cout << "solvent class delete!" << endl;
    }

    void rism3d_solvent_class :: set_solvent(rism3d_c::solvent_cpp *solv_f){
        temperature = solv_f->temperature;
        dielconst = solv_f->dielconst;
        xappa = solv_f->xappa;
        xikt = solv_f->xikt;
        smear = solv_f->smear;
        xikt_dT = solv_f->xikt_dT;
        numAtomTypes = solv_f->numAtomTypes;
        numMolecules = solv_f->numMolecules;
        numRDFpoints = solv_f->numRDFpoints;
        gridSpacingR = solv_f->gridSpacingR;
        gridSpacingK = solv_f->gridSpacingK;
        ionic = solv_f->ionic;
        xvv_version = solv_f->xvv_version;

        atomMultiplicity = memalloc_p->allocInt(atomMultiplicity, numAtomTypes);
        for(int i = 0; i < numAtomTypes; i++){
            atomMultiplicity[i] = solv_f->atomMultiplicity[i];
        }

        numAtoms = memalloc_p->allocInt(numAtoms, numMolecules);
        for(int i = 0; i < numMolecules; i++){
            numAtoms[i] = solv_f->numAtoms[i];
        }
        
        // atomName is allocated here, but is set in a specific function
        atomName = new string[numAtomTypes];

        waveNumbers.set_memalloc(memalloc_p, true);
        waveNumbers.alloc_mem(numRDFpoints);
        for(int i = 0; i < numRDFpoints; i++){
            waveNumbers(i) = solv_f->waveNumbers[i];
        }


        // xvv = new double[numRDFpoints*numAtomTypes*numAtomTypes];
        // xvv = memalloc_p->allocReal(xvv, numRDFpoints, numAtomTypes, numAtomTypes);
        // xvv.alloc_mem(numRDFpoints, numAtomTypes, numAtomTypes);

        // Allocate data with order (numAtomTypes, numAtomTypes, numRDFpoints)
        // Pay attention to the xvv indexing ( (k,j,i) ) to make sure the data is
        // being copied correctly from Fortran's solv_f->xvv.
        // Data comes from Fortran as column major and we are copying it to C++
        // side as row major.
        xvv.set_memalloc(memalloc_p, true);
        xvv.alloc_mem(numAtomTypes, numAtomTypes, numRDFpoints);
        for(int i = 0; i < numRDFpoints; i++){
            for(int j = 0; j < numAtomTypes; j++){
                for(int k = 0; k < numAtomTypes; k++){
                    // old version for when xvv was column major on C++ side
                    // xvv(i,j,k) = solv_f->xvv[i + numRDFpoints*j + numRDFpoints*numAtomTypes*k];

                    // getting xvv as row major on C++ side
                    xvv(k,j,i) = solv_f->xvv[i + numRDFpoints*j + numRDFpoints*numAtomTypes*k];
                }
            }
        }

        // Note that we are not using xvv_dT so far, so this should be changed to get it
        // in the correct (row major) order
        xvv_dT = memalloc_p->allocReal(xvv_dT, numRDFpoints, numAtomTypes, numAtomTypes);
        for(int i = 0; i < numRDFpoints; i++){
            for(int j = 0; j < numAtomTypes; j++){
                for(int k = 0; k < numAtomTypes; k++){
                    xvv_dT[i + numRDFpoints*j + numRDFpoints*numAtomTypes*k] = 
                        solv_f->xvv_dT[i + numRDFpoints*j + numRDFpoints*numAtomTypes*k];
                }
            }
        }


        // charge = memalloc_p->allocReal(charge, numAtomTypes);
        charge.set_memalloc(memalloc_p, true);
        charge.alloc_mem(numAtomTypes);

        charge_db.set_memalloc(memalloc_p, true);
        charge_db.alloc_mem(numAtomTypes);

        for(int i = 0; i < numAtomTypes; i++){
            charge(i) = solv_f->charge[i];
        }
        
        charge_sp.set_memalloc(memalloc_p, true);
        charge_sp.alloc_mem(numAtomTypes);
        for(int i = 0; i < numAtomTypes; i++){
            charge_sp(i) = solv_f->charge_sp[i];
        }

        // Density is coming from Fortran with double precision,
        // so here we copy it into a specific variable and then
        // convert it to float if necessary.
        // We could probably free density_db below
        density.set_memalloc(memalloc_p);
        density.alloc_mem(numAtomTypes);

        density_db.set_memalloc(memalloc_p);
        density_db.alloc_mem(numAtomTypes);

        for(int i = 0; i < numAtomTypes; i++){
            density_db.m_data[i] = solv_f->density[i];
// Converting from double to float as needed, because density
// comes through Shroud interface as double
#if defined(RISMCUDA_DOUBLE)
            density.m_data[i] = solv_f->density[i];
#else
            density.m_data[i] = (float)solv_f->density[i];
#endif //RISMCUDA_DOUBLE
        }
        
        density_sp = memalloc_p->allocReal(density_sp, numMolecules);
        for(int i = 0; i < numMolecules; i++){
            density_sp[i] = solv_f->density_sp[i];
        }
        
        ljSigma.set_memalloc(memalloc_p);
        ljSigma.alloc_mem(numAtomTypes);
        for(int i = 0; i < numAtomTypes; i++){
            ljSigma(i) = solv_f->ljSigma[i];
        }
        
        ljEpsilon.set_memalloc(memalloc_p);
        ljEpsilon.alloc_mem(numAtomTypes);
        for(int i = 0; i < numAtomTypes; i++){
            ljEpsilon(i) = solv_f->ljEpsilon[i];
        }

        int max_atomMultiplicity = *max_element(atomMultiplicity, atomMultiplicity + numAtomTypes);
        int max_numAtoms = *max_element(numAtoms, numAtoms + numMolecules);

        coord.set_memalloc(memalloc_p);
        coord.alloc_mem(3, max_atomMultiplicity, max_numAtoms, numMolecules);
        
        // Here the coord variable is being handled with column major ordering
        // This is a residue from the beginning of my implementation. Changing this
        // to row major ordering will require looking at the cuda kernels that uses
        // this variable to make appropriate changes.
        // Since it is working properly so far, that will be done in the future.
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < max_atomMultiplicity; j++){
                for(int k = 0; k < max_numAtoms; k++){
                    for(int q = 0; q < numMolecules; q++){
                        coord(i,j,k,q) = solv_f->coord[i + 3*j + 3*max_atomMultiplicity*k + 3*max_atomMultiplicity*max_numAtoms*q];
                    }
                }
            }
        }

        background_correction = memalloc_p->allocReal(background_correction, numAtomTypes);
        for(int i = 0; i < numAtomTypes; i++){
            background_correction[i] = solv_f->background_correction[i];
        }

        delhv0.set_memalloc(memalloc_p,true);
        delhv0.alloc_mem(numAtomTypes);
        for(int i = 0; i < numAtomTypes; i++){
            delhv0(i) = solv_f->delhv0[i];
        }

        delhv0_dT = memalloc_p->allocReal(delhv0_dT, numAtomTypes);
        for(int i = 0; i < numAtomTypes; i++){
            delhv0_dT[i] = solv_f->delhv0_dT[i];
        }

    }
    
}