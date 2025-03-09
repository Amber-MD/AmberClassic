#include <stdio.h>
#include <cuda.h> 
#include <cuda_runtime.h>
#include "rism_util.hpp"
using namespace std;

namespace rism3d_c {

    __device__ void polynomialInterpolation(GPUtype *x0, GPUtype* y0, int n, GPUtype x, GPUtype *y){
        GPUtype *p = new GPUtype[n];
        // Copy values from y0 into p
        for(int j = 0; j < n; j++){
            p[j] = y0[j];
        }
        
        for(int m = 0; m < n - 2; m++){
            for(int i = 0; i < n-m-1; i++){
                p[i] = (x - x0[i + m + 1]) * p[i] - (x - x0[i]) * p[i+1];
                p[i] = p[i] / (x0[i] - x0[i+m + 1]);
            }
        }
        
        *y = ((x - x0[n-1]) * p[0] - (x - x0[0]) * p[1]) / (x0[0] - x0[n-1]);

        GPUtype error = (p[0] + p[1] - 2.0f * (*y)) / 2.0f;
        if(error > 1e-4){
            printf("Error = %f. Greater than 1e-4! \n", error);
        }
        delete p;
    
    }

    __global__ void polynomialInterpolation_kernel(int maxPointsToInterp, int numPoints, int numPointsToInterp,
                                            GPUtype *x0, GPUtype *y0, GPUtype *x, GPUtype *y, int iv1, int iv2, int Nyz){
        // numPoints == numRDFpoints
        // numPointsToInterp == waveNumberArraySize
        // x0 == solventWaveNumbers
        // y0 == xvv
        // x == gridWaveNumbers
        // y == xvva
        // Nyz = number of solvent atom types
        

        int idx = threadIdx.x + blockIdx.x * blockDim.x;

        if(idx < numPointsToInterp){
            int igk1 = 1;

            for(int igk = 0; igk < numPoints - maxPointsToInterp + 1; igk++){
                igk1 = igk;
                if(x0[igk1 + maxPointsToInterp/2] > x[idx]){
                    break;
                }
            }

            polynomialInterpolation(x0 + igk1, y0 + igk1, maxPointsToInterp, x[idx], &y[iv1*Nyz*numPointsToInterp + iv2*numPointsToInterp + idx]);

        }

    }

    __global__ void memcpy_kernel(GPUtype *src, GPUtype *dst, int size){
        int idx = threadIdx.x + blockIdx.x * blockDim.x;
        
        if(idx < size){
            dst[idx] = src[idx];
        }
    }

    __global__ void reduction_kernel(GPUtype* input, GPUtype* output, int N) {
        extern __shared__ GPUtype shared_data[];

        int idx = threadIdx.x + blockIdx.x * blockDim.x; // global index
        int tid = threadIdx.x; // thread's local id within the block (ranging from 0 to blockDim.x - 1).

        // Load elements into shared memory
        if (idx < N) {
            shared_data[tid] = input[idx];
        } else {
            shared_data[tid] = 0.0f;
        }

        __syncthreads();  // Synchronize threads in the block

        // Perform parallel reduction in shared memory:
        // This step will sum pair of values while reduicing the
        // the number of elements to be computed. Ex.: for the first block
        // containing 8 values
        // arr = [0 1 2 3 4 5 6 7]
        // for stride = 1:
        //     tid = 0, 2, 4, 6
        //     arr = [1 1 5 3 9 5 13 7]
        // for stride = 2:
        //     tid = 0, 4
        //     arr = [6 1 5 3 22 5 13 7]
        // for stride = 4
        //     tid = 0
        //     arr = [28 1 5 3 22 5 13 7]
        for (int stride = 1; stride < blockDim.x; stride *= 2) {
            if (tid % (2 * stride) == 0) {
                shared_data[tid] += shared_data[tid + stride];
            }
            __syncthreads();  // Ensure all threads have updated shared memory
        }

        // Write the result of the reduction in the first thread of the block
        // The value of the first thread of each block (which stores the total
        // sum for that block) will be asigned to the output array at the block 
        // id position
        if (tid == 0) {
            output[blockIdx.x] = shared_data[0];
        }
    }

    void test(GPUtype *xvv){
        for(int i = 0; i < 100; i++){
            cout << xvv[i] << endl;
        }
    }

    void cu_memcpy(GPUtype *src, GPUtype *dst, int size){
        int num_blocks = (size + 255) / 256;
        int num_threads = 256;

        memcpy_kernel<<<num_blocks, num_threads>>>(src, dst, size);
        cudaDeviceSynchronize();
    }

    void cu_polinomialInterpolation(int maxPointsToInterp, int waveNumberArraySize, 
                                    int numRDFpoints, int numAtomTypes,
                                    GPUtype *solventWaveNumbers, 
                                    GPUtype *xvv, GPUtype *gridWaveNumbers,
                                    GPUtype *xvva, int iv1, int iv2, int Nyz){
        

        ////////////////////////////////////////////////////////////////////////////////////////
        int num_blocks = (waveNumberArraySize+255)/256;
        int num_threads = 256;

        // test(xvv);

        // int device = -1;
        // cudaGetDevice(&device);

        // cudaMemPrefetchAsync(solventWaveNumbers, numRDFpoints*sizeof(GPUtype),device,NULL);
        // cudaMemPrefetchAsync(xvv, numRDFpoints*numAtomTypes*numAtomTypes*sizeof(GPUtype),device,NULL);
        // cudaMemPrefetchAsync(gridWaveNumbers, waveNumberArraySize*sizeof(GPUtype),device,NULL);

        // // Here I think we can just reserve the ammount of memory for each pair of solvent atom types,
        // // since we do not need memory beyond that for each kernel call and each kernel call will start at
        // // the correct contigous memory address.
        // cudaMemPrefetchAsync(xvva, waveNumberArraySize*sizeof(GPUtype),device,NULL);

        polynomialInterpolation_kernel<<<num_blocks,num_threads>>>(maxPointsToInterp, numRDFpoints, 
                                                                   waveNumberArraySize, 
                                                                   solventWaveNumbers, 
                                                                   xvv, 
                                                                   gridWaveNumbers, 
                                                                   xvva, iv1, iv2, Nyz);
        
    }

    GPUtype cu_reduce_sum(GPUtype* arr, int N){
        int block_size = 256; // number of threads in each block
        int grid_size = (N + block_size - 1) / block_size; // total number of blocks to be used
        int shared_memory_size = block_size * sizeof(GPUtype);

        GPUtype *output;
        cudaMallocManaged(&output, (N + 255) / 256 * sizeof(GPUtype));  // One per block

        // Launch the kernel
        // This will ensure that there will be lauched grid_size blocks,
        // each containing 256 threads and using enough shared memory to 
        // store shared_memory_size float values
        reduction_kernel<<<grid_size, block_size, shared_memory_size>>>(arr, output, N);

        cudaDeviceSynchronize();  // Wait for the kernel to finish

        // Final reduction on the host
        GPUtype final_sum = 0.0;
        for (int i = 0; i < grid_size; i++) {
            final_sum += output[i];
        }

        // Cleanup
        cudaFree(output);

        return final_sum;
    }

}