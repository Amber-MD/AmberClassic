#ifndef RISM3DGRID_HPP
#define RISM3DGRID_HPP

#include <iostream>
#include <fstream>
#include <cuda.h> 
#include <cuda_runtime.h>
#include "rism3d_safemem.hpp"
#include <math.h>
#include "array_class.hpp"
#include "rism_util.hpp"
#include "varTypes.hpp"
using namespace std;

namespace rism3d_c {

    class rism3d_grid{
        private:
            rism3d_safemem* memalloc_p;
            
        public:
            //>  Grid spacing for all of the grids. [A]
            GPUtype spacing[3];

            //> Box size for 3D-RISM.  For PBC calculations, these should
            GPUtype boxLength[3];

            //> Volume of the box we are doing the calculation. [A^3]
            GPUtype boxVolume;

            //> Unit cell coordinate vectors of length unity, in Cartesian
            //  coordinates of the box grid.  These are derived from the unit
            //  cell interior angles (axis, vector). [A]
            array_class<double> unitCellUnitVectorsR; // It needs to be changed from <double> to <GPUtype> in the future

            //> Unit cell interior angles. [radians]
            double unitCellAngles[3];

            //> Unit cell coordinate vectors with lengths of the respective
            //  unit cell side, in Cartesian coordinates of the box grid.
            //  These are derived from the unit cell side lengths and interior
            //  angles (axis, vector). [A]
            array_class<double> unitCellVectorsR; // It needs to be changed from <double> to <GPUtype> in the future
            double unitCellVectorsK[3*3];

            //> Unit cell side lengths. [A]
            double unitCellLengths[3];

            //> Volume of a grid cell (the smallest volumetric grid element). [A^3]
            GPUtype voxelVolume;

            //> Unit cell coordinate vectors with length of the grid spacing
            //  in each respective dimension. These are derived from the grid
            //  spacing and unit cell interior angles (axis, vector). [A]
            array_class<double> voxelVectorsR; // It needs to be changed from <double> to <GPUtype> in the future
            array_class<double> voxelVectorsK; // It needs to be changed from <double> to <GPUtype> in the future

            //> Absolute value of the wave vector, |k|, also known as the wave
            //  number (this%waveNumberArraySize long).
            array_class<GPUtype> waveNumbers;

            //> Wave vector at every k-space grid point (vector, k-space grid points).
            //  By the time I did this implementation, it was easier to work of the x,y,z components
            //  instead of a single array.
            array_class<GPUtype> waveVectorX;
            array_class<GPUtype> waveVectorY;
            array_class<GPUtype> waveVectorZ;

            //> Self dot product of the wave vector at every k-space grid point.
            array_class<GPUtype> waveVectors2;

            //> Unit cell lattice vectors in frequency space.
            double inscribedSphereRadius;

            //> Map wave vector index to wave number index in ascending wave
            //  number value.
            array_class<int> waveVectorWaveNumberMap;

            //> True if grid has a defined periodic unit cell.
            bool periodic;

            //> Total number of MPI local r-space grid points.
            int totalLocalPointsR;

            //> Total number of MPI local k-space grid points.
            int totalLocalPointsK;

            //> Total number of global r-space grid points.
            int totalGlobalPointsR;

            //> Number of indices in array of wave vector magnitudes.
            long int waveNumberArraySize;

            //> Number of MPI local r-space grid spaces in each dimension.
            int localDimsR[3];

            //> Number of MPI local k-space grid spaces in each dimension.
            int localDimsK[3];

            //> Number of global r-space grid spaces in each dimension.
            int globalDimsR[3];

            //> Number of global k-space grid spaces in each dimension.
            int globalDimsK[3];

            //> r-space grid offset from (0, 0, 0) for MPI calculations.
            int OffsetR[3];

            //> k-space grid offset from (0, 0, 0) for MPI calculations.
            int OffsetK[3];

            int mpirank;
            int mpisize;
            int mpicomm;

            rism3d_grid(rism3d_safemem *memalloc);
            ~rism3d_grid();

            void set_spacing(double *grdspc);

            void resize(GPUtype gridSpacing[3], int gDimsR[3], int gDimsK[3], 
                        int lDimsR[3], int lDimsK[3], int off_setR[3], int off_setK[3]);

            void setup_wavevector();

            void setUnitCellDimensions(double *unitCellDimensions, bool periodic_val);

    };

}

#endif // RISM3DGRID_HPP