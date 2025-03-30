#include <stdio.h>
#include <cuda.h> 
#include <cuda_runtime.h>
#include "rism3d_grid.hpp"
using namespace std;

namespace rism3d_c {

    /////////////// KERNELS ///////////////

    __global__ void k_setup_wavevector(int n, int rx, int ry, int rz,
                                       GPUtype *wavevectorX, GPUtype *wavevectorY, GPUtype *wavevectorZ,
                                       GPUtype fatbLx, GPUtype fatbLy, GPUtype fatbLz,
                                       GPUtype *waveVectors2){
        int lgz, lgy, lgx, igk, igk2;
        GPUtype waveX, waveY, waveZ;
        int indX = blockIdx.x * blockDim.x + threadIdx.x;
        int indY = blockIdx.y * blockDim.y + threadIdx.y;
        int indZ = blockIdx.z * blockDim.z + threadIdx.z;

        // old column major version
        // if(indX >= rx/2 || indY >= ry || indZ >= rz){
        //     return;
        // }

        // new row major version
        // switching from rx/2,ry,rz to rx,ry,rz/2 is exchanging the values 
        // for WavevectorX and wavevectorZ while comparing to Fortran version.
        // This means we are considering different subvolumes in Fortran and
        // C++ implementations while working in the reciprocal space. So, results
        // of the direct Fourier transforms cannot be directly compared.
        // The wavenumber values does not change, though
        if(indX >= rx || indY >= ry || indZ >= rz/2){
            return;
        }

        // old column major version
        // igk = indX + (indY + indZ * ry) * rx / 2;

        // new row major version
        igk = indZ + (indY + indX * ry) * rz / 2;
    
        lgz = ((indZ + (rz / 2 - 1)) % rz) - (rz / 2 - 1);
        lgy = ((indY + (ry / 2 - 1)) % ry) - (ry / 2 - 1);
        lgx = ((indX + (rx / 2 - 1)) % rx) - (rx / 2 - 1);

        waveX = fatbLx*lgx;
        waveY = fatbLy*lgy;
        waveZ = fatbLz*lgz;

        wavevectorX[igk] = waveX;
        wavevectorY[igk] = waveY;
        wavevectorZ[igk] = waveZ;

        waveVectors2[igk] = waveX*waveX + waveY*waveY + waveZ*waveZ;

        // Nyquist frequencies

        // old column major version
        // igk2 = indY + (indZ + rz * rx/2) * ry;

        // new row major version
        igk2 = indY + (indX + rx * rz/2) * ry;

        // old column major version
        // GPUtype waveX2;
        // int lgx2;
        // lgx2 = ((rx / 2 + (rx / 2 - 1)) % rx) - (rx / 2 - 1);
        // waveX2 = fatbLx*lgx2;

        // new row major version
        GPUtype waveZ2;
        int lgz2;
        lgz2 = ((rz / 2 + (rz / 2 - 1)) % rz) - (rz / 2 - 1);
        waveZ2 = fatbLz*lgz2;

        wavevectorX[igk2] = waveX;
        wavevectorY[igk2] = waveY;
        wavevectorZ[igk2] = waveZ2;

        // old column major version
        // waveVectors2[igk2] = waveX2*waveX2 + waveY*waveY + waveZ*waveZ;

        // new row major version
        waveVectors2[igk2] = waveX*waveX + waveY*waveY + waveZ2*waveZ2;

    }

    /////////////// KERNEL WRAPPER FUNCTIONS ///////////////

    void rism3d_grid :: setup_wavevector(){
        GPUtype pi;
        GPUtype factor[3];

        // Defining and allocating memory for waveVectorToWaveVector2Map
        array_class<int> waveVectorToWaveVector2Map;
        waveVectorToWaveVector2Map.set_memalloc(memalloc_p);
        waveVectorToWaveVector2Map.alloc_mem(totalLocalPointsK/2);

        // allocate memory for waveVectorWaveNumberMap
        waveVectorWaveNumberMap.set_memalloc(memalloc_p, true);
        waveVectorWaveNumberMap.alloc_mem(totalLocalPointsK/2);

#if RISMCUDA_DOUBLE
        pi = M_PI;
#else
        pi = (float)M_PI;
#endif // RISMCUDA_DOUBLE

        waveVectorX.set_memalloc(memalloc_p, true);
        waveVectorX.alloc_mem(totalLocalPointsK/2);
        cudaMemPrefetchAsync(waveVectorX.m_data, (totalLocalPointsK/2)*sizeof(GPUtype), 0);

        waveVectorY.set_memalloc(memalloc_p, true);
        waveVectorY.alloc_mem(totalLocalPointsK/2);
        cudaMemPrefetchAsync(waveVectorY.m_data, (totalLocalPointsK/2)*sizeof(GPUtype), 0);

        waveVectorZ.set_memalloc(memalloc_p, true);
        waveVectorZ.alloc_mem(totalLocalPointsK/2);
        cudaMemPrefetchAsync(waveVectorZ.m_data, (totalLocalPointsK/2)*sizeof(GPUtype), 0);

        waveVectors2.set_memalloc(memalloc_p, true);
        waveVectors2.alloc_mem(totalLocalPointsK/2);
        cudaMemPrefetchAsync(waveVectors2.m_data, (totalLocalPointsK/2)*sizeof(GPUtype), 0);
        cudaDeviceSynchronize();
        cudaMemset(waveVectors2.m_data, 0, totalLocalPointsK/2 * sizeof(GPUtype));

        // ?? Do we need a synchronization step here?

        // The number of blocks in each thread dimention may need to be adjusted for optimization
        int n_blocksx = 8;
        int n_blocksy = 8;
        int n_blocksz = 8;
        dim3 block(n_blocksx,n_blocksy,n_blocksz);
        dim3 grid((globalDimsR[0] + n_blocksx-1) / n_blocksx, (globalDimsR[1] + n_blocksy-1) / n_blocksy, (globalDimsR[2] + n_blocksz-1) / n_blocksz);

        for(int i = 0; i < 3; i++){
            factor[i] = (2.0 * pi) / boxLength[i];
        }

        k_setup_wavevector<<<grid,block>>>(totalLocalPointsK/2,
                                           globalDimsR[0],globalDimsR[1],globalDimsR[2],
                                           waveVectorX.m_data, waveVectorY.m_data, waveVectorZ.m_data,
                                           factor[0],factor[1],factor[2], waveVectors2.m_data);
        cudaDeviceSynchronize();

        // Getting waveVectors2 data back to CPU
        cudaMemPrefetchAsync(waveVectors2.m_data, (totalLocalPointsK/2)*sizeof(GPUtype), cudaCpuDeviceId);

        // waveVectorToWaveVector2Map is an index array that will retrieve waveVectors2 in 
        // ascending order
        // Ideally, this indexArray function would be implemented as a CUDA kernel. However,
        // I think it is ok to have this as a serial implementation because maybe the speed 
        // up would not be significant. Also, the siftdown function may not be trivial to
        // parallelize
        indexArray(waveVectors2.m_data, waveVectorToWaveVector2Map.m_data, totalLocalPointsK/2);

        // print test for waveVecotrToWaveVector2Map
        // Comparing to Fortran results may have some different indexes because of repeated
        // values (e.g. 0.026127333114358848) happening more than once. Thus, some indexes 
        // refering to the same value can appear with switched positions in Fortran.
        // This print statement shows that using waveVectorToWaveVector2Map ensures
        // waveVectors2 will be get in ascending order
        // for(int i = 0; i < 20; i++){
        //     cout << "waveVectors2[" << " waveVectorToWaveVector2Map[" << i << "]] = " << waveVectors2(waveVectorToWaveVector2Map(i)) << endl;
        // }
        // for(int i = (totalLocalPointsK/2)-10; i < totalLocalPointsK/2; i++){
        //     cout << "waveVectors2[" << " waveVectorToWaveVector2Map[" << i << "]] = " << waveVectors2(waveVectorToWaveVector2Map(i)) << endl;
        // }

        waveNumbers.set_memalloc(memalloc_p, true);

        // Ideally, this sort_unique function would be implemented as a CUDA kernel,
        // since it is working with waveNumbers and waveVectorWaveNumberMap.
        // However, I think we can left it for later because:
        //     1 - It may not be that trivial to parallelize it
        //     2 - I do not think we would get a significant speed up
        //     3 - I think the slow down of sending data from CPU to GPU in this case
        //         would not be significant, because would be done only once.
        sort_unique(&waveVectors2, totalLocalPointsK/2, &waveNumbers, &waveVectorToWaveVector2Map, &waveVectorWaveNumberMap);

        // print test for waveVectorWaveNumberMap
        // Here we can see that waveNumbers(waveVectorWaveNumberMap(i)) will  return the correct value of
        // sqrt(waveVectors2(i)). Thus, waveVectorWaveNumberMap(i) maps the indexing correctly.
        // for(int i = 0; i < 20; i++){
        //     cout << "waveVectors2[" << i << "] = " << waveVectors2(i) << " || waveNumber = " << waveNumbers(waveVectorWaveNumberMap(i)) << endl;
        // }
        // for(int i = (totalLocalPointsK/2)-10; i < totalLocalPointsK/2; i++){
        //     cout << "waveVectors2[" << i << "] = " << waveVectors2(i) << " || waveNumber = " << waveNumbers(waveVectorWaveNumberMap(i)) << endl;
        // }

        waveNumberArraySize = waveNumbers.get_dim1();
        
    }

}