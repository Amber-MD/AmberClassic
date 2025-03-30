#include <iostream>
#include "rism3d_solute_class.hpp"
using namespace std;

namespace rism3d_c {

    rism3d_solute_class :: rism3d_solute_class() : 
                           position(array_class<GPUtype>::ROW_MAJOR)
    {
        // cout << "Solute class created!" << endl;
    }

    rism3d_solute_class :: rism3d_solute_class(rism3d_safemem *memalloc, rism3d_c::solute_cpp *solu_f) : 
                           position(array_class<GPUtype>::ROW_MAJOR)
    {
        // cout << "Solute class created!" << endl;

        memalloc_p = memalloc;

        set_solute(solu_f);



    }

    rism3d_solute_class :: ~rism3d_solute_class(){
        // cout << "Solute class deleted!" << endl;
    }

    void rism3d_solute_class :: set_solute(rism3d_c::solute_cpp *solu_f){
        numAtoms = solu_f->numAtoms;
        charged = solu_f->charged;
        totalCharge = solu_f->totalCharge;
        
        mass.set_memalloc(memalloc_p);
        mass.alloc_mem(numAtoms);
        for(int i = 0; i < numAtoms; i++){
            mass(i) = solu_f->mass[i];
        }

        charge.set_memalloc(memalloc_p, true);
        charge.alloc_mem(numAtoms);
        for(int i = 0; i < numAtoms; i++){
            charge(i) = solu_f->charge[i];
        }

        // origCharge holds the orignal values in charge() when _unsetcharges()
        // is used.  Otherwise it is not used. [sqrt(kT A)]
        //// origCharge = new double[numAtoms];
        // origCharge.set_memalloc(memalloc_p);
        // origCharge.alloc_mem(numAtoms);
        // for(int i = 0; i < numAtoms; i++){
        //     // origCharge[i] = solu_f->origCharge[i];
        //     solu_f->origCharge[i] = solu_f->charge[i];
        // }

        ljSigma.set_memalloc(memalloc_p);
        ljSigma.alloc_mem(numAtoms);

        for(int i = 0; i < numAtoms; i++){
            ljSigma(i) = solu_f->ljSigma[i];
        }

        ljEpsilon.set_memalloc(memalloc_p);
        ljEpsilon.alloc_mem(numAtoms);
        for(int i = 0; i < numAtoms; i++){
            ljEpsilon(i) = solu_f->ljEpsilon[i];
        }

        centerOfMass.set_memalloc(memalloc_p);
        centerOfMass.alloc_mem(3);
        for(int i = 0; i < 3; i++){
            centerOfMass(i) = solu_f->centerOfMass[i];
        }

        for(int i = 0; i < 3; i++){
            translation[i] = solu_f->translation[i];
        }

        // The coordinates received here are null. So memory is allocated and coordinates
        // set in a later call
        position.set_memalloc(memalloc_p, true);
        position.alloc_mem(3,numAtoms);
    }

    void rism3d_solute_class :: translate_solute(bool init_position){
        GPUtype to_translate[3];
        if(init_position == true){
            for(int i = 0; i < 3; i++){
                to_translate[i] = -translation[i];
            }
        }
        else{
            for(int i = 0; i < 3; i++){
                to_translate[i] = translation[i];
            }
        }
        translate(&position, numAtoms, to_translate);
    }

    void rism3d_solute_class :: setCoord(double* solutePositions){
        // if(charged == true){
#if defined(CUDA_NOSWAP)
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < numAtoms; j++){
                    position(i,j) = solutePositions[i + j*3];
                }
            }
#else
            for(int j = 0; j < numAtoms; j++){
                position(0,j) = solutePositions[2 + j*3];
                position(1,j) = solutePositions[1 + j*3];
                position(2,j) = solutePositions[0 + j*3];
            }
#endif
        // } else{
        //     for(int i = 0; i < 3; i++){
        //         for(int j = 0; j < numAtoms; j++){
        //             position(i,j) = solutePositions[i + j*3];
        //         }
        //     }
        // }

// #if RISMCUDA_DOUBLE
//             ofstream myfile2;
//             myfile2.open("/home/fcarvalho/rism3d.cuda.test/coords_db.txt");
//             for(int i = 0; i < numAtoms; i++){
//                 myfile2 << setprecision(17) << position(0,i) << ", " << position(1,i) << ", " << position(2,i) << "\n";
//             }
//             myfile2.close();
// #else
//             ofstream myfile2;
//             myfile2.open("/home/fcarvalho/rism3d.cuda.test/coords.txt");
//             for(int i = 0; i < numAtoms; i++){
//                 myfile2 << setprecision(17) << position(0,i) << ", " << position(1,i) << ", " << position(2,i) << "\n";
//             }
//             myfile2.close();
// #endif //RISMCUDA_DOUBLE
    }

}