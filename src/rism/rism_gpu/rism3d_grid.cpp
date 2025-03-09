#include <iostream>
#include "rism3d_grid.hpp"
using namespace std;

namespace rism3d_c {

    rism3d_grid :: rism3d_grid(rism3d_safemem *memalloc){
        // cout << "Grid object created" << endl;

        memalloc_p = memalloc;

        unitCellUnitVectorsR.set_memalloc(memalloc);
        unitCellUnitVectorsR.alloc_mem(3,3);

        unitCellVectorsR.set_memalloc(memalloc);
        unitCellVectorsR.alloc_mem(3,3);

        voxelVectorsR.set_memalloc(memalloc);
        voxelVectorsR.alloc_mem(3,3);

        voxelVectorsK.set_memalloc(memalloc);
        voxelVectorsK.alloc_mem(3,3);

        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                // unitCellVectorsK is not yet being used in the cuda version.
                // However, amber_rism_interface has a function to get these values,
                // so I am just assigning values here to have something to be passed back to
                // Fortran.
                unitCellVectorsK[i+3*j] = i+3*j + 119;
            }
        }
    }

    rism3d_grid :: ~rism3d_grid(){
        // cout << "Grid object destroyed" << endl;
    }

    void rism3d_grid :: set_spacing(double *grdspc){
        for(int i = 0; i < 3; i++){
            spacing[i] = grdspc[i];
        }
    }

    // This function is calculating unit cell parameters for periodic calculations.
    // I have started this implementation by calculating unitCellUnitVectorsR, but it
    // still needs to calculate unitCellUnitVectorsK as indicated below
    void rism3d_grid :: setUnitCellDimensions(double *unitCellDimensions, bool periodic_val){
        periodic = periodic_val;
        for(int i = 0; i < 3; i++){
            unitCellLengths[i] = unitCellDimensions[i];
            unitCellAngles[i] = unitCellDimensions[i+3] * M_PI / 180;

        // setting unitCellUnitVectorsR
        unitCellUnitVectorsR(0,0) = 1;
        unitCellUnitVectorsR(0,1) = 0;
        unitCellUnitVectorsR(0,2) = 0;

        unitCellUnitVectorsR(1,0) = cos(unitCellAngles[2]);
        unitCellUnitVectorsR(1,1) = sin(unitCellAngles[2]);
        unitCellUnitVectorsR(1,2) = 0;

        unitCellUnitVectorsR(2,0) = cos(unitCellAngles[1]);
        unitCellUnitVectorsR(2,1) = (cos(unitCellAngles[0]) - 
        unitCellUnitVectorsR(2, 0)*unitCellUnitVectorsR(1, 0)) / unitCellUnitVectorsR(1, 1);
        unitCellUnitVectorsR(2,2) = sqrt(1.0 - pow(unitCellUnitVectorsR(2, 0),2) - 
        pow(unitCellUnitVectorsR(2, 1),2));

        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                if(unitCellUnitVectorsR(i,j) < 1e-10){
                    unitCellUnitVectorsR(i,j) = 0;
                }
            }
        }

        // setting unitCellVectorsR
        for(int j = 0; j < 3; j++){
            unitCellVectorsR(0,j) = unitCellLengths[0]*unitCellUnitVectorsR(0,j);
            unitCellVectorsR(1,j) = unitCellLengths[1]*unitCellUnitVectorsR(1,j);
            unitCellVectorsR(2,j) = unitCellLengths[2]*unitCellUnitVectorsR(2,j);

            voxelVectorsR(0, j) = spacing[0] * unitCellUnitVectorsR(0, j);
            voxelVectorsR(1, j) = spacing[1] * unitCellUnitVectorsR(1, j);
            voxelVectorsR(2, j) = spacing[2] * unitCellUnitVectorsR(2, j);
        }

        // Still need to implement this part from Fortran code to get unitCellUnitVectorsK
        /*  u23 = cross(this%unitCellVectorsR(2, :), this%unitCellVectorsR(3, :))
            u31 = cross(this%unitCellVectorsR(3, :), this%unitCellVectorsR(1, :))
            u12 = cross(this%unitCellVectorsR(1, :), this%unitCellVectorsR(2, :))
            unitCellVolume = dot_product(this%unitCellVectorsR(1, :), u23)
            this%unitCellVectorsK(1, :) = u23(:) / unitCellVolume
            this%unitCellVectorsK(2, :) = u31(:) / unitCellVolume
            this%unitCellVectorsK(3, :) = u12(:) / unitCellVolume  */
        
        }
    }

    void rism3d_grid :: resize(GPUtype gridSpacing[3], int gDimsR[3], int gDimsK[3], 
                               int lDimsR[3], int lDimsK[3], int off_setR[3], int off_setK[3]){
        double unitCellDimensions[6];

        spacing[0] = gridSpacing[0];
        spacing[1] = gridSpacing[1];
        spacing[2] = gridSpacing[2];

        globalDimsR[0] = gDimsR[0];
        globalDimsR[1] = gDimsR[1];
        globalDimsR[2] = gDimsR[2];

        globalDimsK[0] = gDimsK[0];
        globalDimsK[1] = gDimsK[1];
        globalDimsK[2] = gDimsK[2];

        localDimsR[0] = lDimsR[0];
        localDimsR[1] = lDimsR[1];
        localDimsR[2] = lDimsR[2];

        localDimsK[0] = lDimsK[0];
        localDimsK[1] = lDimsK[1];
        localDimsK[2] = lDimsK[2];

        OffsetR[0] = off_setR[0];
        OffsetR[1] = off_setR[1];
        OffsetR[2] = off_setR[2];

        totalLocalPointsR = localDimsR[0]*localDimsR[1]*localDimsR[2];
        totalLocalPointsK = localDimsK[0]*localDimsK[1]*localDimsK[2];
        totalGlobalPointsR = globalDimsR[0]*globalDimsR[1]*globalDimsR[2];

        for(int i = 0; i < 3; i++){
            voxelVectorsR(i,0) = spacing[i]*unitCellUnitVectorsR(i,0);
            voxelVectorsR(i,1) = spacing[i]*unitCellUnitVectorsR(i,1);
            voxelVectorsR(i,2) = spacing[i]*unitCellUnitVectorsR(i,2);
            for(int j = 0; j < 3; j++){
                if(voxelVectorsR(i,j) != 0){
                    voxelVectorsK(i,j) = 1/voxelVectorsR(i,j);
                }
                else{
                    voxelVectorsK(i,j) = 0;
                }
            }
            boxLength[i] = globalDimsR[i]*spacing[i];
        }

        if(periodic == false){
            for(int i = 0; i < 6; i++){
                if(i < 3){
                    unitCellDimensions[i] = boxLength[i];
                }
                else{
                    unitCellDimensions[i] = 90.0;
                }
            }
            setUnitCellDimensions(unitCellDimensions,periodic);
        }

        if(periodic == true){
            cout << "Periodic systems not supported, yet." << endl;
            cout << "Message printed from 'resize' function at grid object." << endl;
            abort();
        }
        else{
            boxVolume = boxLength[0]*boxLength[1]*boxLength[2];
            voxelVolume = spacing[0]*spacing[1]*spacing[2];
        }

        // I am not sure what this is used for. However, I did not needed this so far.
        // It will be implemented when needed...
        /*
        width_a = abs(dot_product(this%unitCellVectorsR(1,:), &
         cross(this%unitCellVectorsR(2,:), this%unitCellVectorsR(3,:)))) &
         / abs(magnitude(cross(this%unitCellVectorsR(2,:), this%unitCellVectorsR(3,:))))
        width_b = abs(dot_product(this%unitCellVectorsR(2,:), &
                cross(this%unitCellVectorsR(3,:), this%unitCellVectorsR(1,:)))) &
                / abs(magnitude(cross(this%unitCellVectorsR(3,:), this%unitCellVectorsR(1,:))))
        width_c = abs(dot_product(this%unitCellVectorsR(3,:), &
                cross(this%unitCellVectorsR(1,:), this%unitCellVectorsR(2,:)))) &
                / abs(magnitude(cross(this%unitCellVectorsR(1,:), this%unitCellVectorsR(2,:))))
        this%inscribedSphereRadius = 0.5 * min(width_a, width_b, width_c)
        */

        setup_wavevector();

        // WRITE FILE TO CHECK WAVEVECTORS VALUES AGAINS FORTRAN
        // I will keep this here for now just in case we need reference on how to
        // get these files in the future
// #if RISMCUDA_DOUBLE
//         ofstream myfile1;
//         myfile1.open("/home/fcarvalho/rism3d.cuda.test/waveVector_db.txt");
//         for(int i = 0; i < totalLocalPointsK/2; i++){
//             myfile1 << setprecision(17) << waveVectorX(i) << "," << waveVectorY(i) << "," << waveVectorZ(i) << "\n";
//         }
//         myfile1.close();

//         ofstream myfile2;
//         myfile2.open("/home/fcarvalho/rism3d.cuda.test/waveNumbers_db.txt");
//         for(size_t i = 0; i < waveNumbers.get_dim1(); i++){
//             myfile2 << setprecision(17) << waveNumbers(i) << "\n";
//         }
//         myfile2.close();

//         ofstream myfile3;
//         myfile3.open("/home/fcarvalho/rism3d.cuda.test/waveVectors2_db.txt");
//         for(int i = 0; i < totalLocalPointsK/2; i++){
//             myfile3 << setprecision(17) << waveVectors2(i) << "\n";
//         }
//         myfile3.close();
// #else
//         ofstream myfile1;
//         myfile1.open("/home/fcarvalho/rism3d.cuda.test/waveVector.txt");
//         for(int i = 0; i < totalLocalPointsK/2; i++){
//             myfile1 << setprecision(17) << waveVectorX(i) << "," << waveVectorY(i) << "," << waveVectorZ(i) << "\n";
//         }
//         myfile1.close();

//         ofstream myfile2;
//         myfile2.open("/home/fcarvalho/rism3d.cuda.test/waveNumbers.txt");
//         for(size_t i = 0; i < waveNumbers.get_dim1(); i++){
//             myfile2 << setprecision(17) << waveNumbers(i) << "\n";
//         }
//         myfile2.close();

//         ofstream myfile3;
//         myfile3.open("/home/fcarvalho/rism3d.cuda.test/waveVectors2.txt");
//         for(int i = 0; i < totalLocalPointsK/2; i++){
//             myfile3 << setprecision(17) << waveVectors2(i) << "\n";
//         }
//         myfile3.close();
// #endif //RISMCUDA_DOUBLE

    }

}