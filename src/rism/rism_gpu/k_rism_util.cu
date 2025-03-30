#include <stdio.h>
#include <cuda.h> 
#include <cuda_runtime.h>
#include "rism_util.hpp"
using namespace std;

#if defined(__CUDA_ARCH__ ) && __CUDA_ARCH__ < 600
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

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

    // Templates for cu_memcpy. These allow us copy memory to variables with different precision

    template <typename SrcType, typename DstType>
    __global__ void memcpy_kernel(const SrcType *src, DstType *dst, int size){
        int idx = threadIdx.x + blockIdx.x * blockDim.x;
        
        if(idx < size){
            dst[idx] = src[idx];
        }
    }

    template __global__ void memcpy_kernel<float, float>(const float *src, float *dst, int size);
    template __global__ void memcpy_kernel<double, double>(const double *src, double *dst, int size);
    template __global__ void memcpy_kernel<double, float>(const double *src, float *dst, int size);
    template __global__ void memcpy_kernel<float, double>(const float *src, double *dst, int size);

    template <typename SrcType, typename DstType>
    void cu_memcpy(const SrcType *src, DstType *dst, int size){
        int num_blocks = (size + 255) / 256;
        int num_threads = 256;

        memcpy_kernel<SrcType,DstType><<<num_blocks, num_threads>>>(src, dst, size);
        cudaDeviceSynchronize();
    }

    void cu_memcpy(const float *src, float *dst, int size) {
        cu_memcpy<float, float>(src, dst, size);
    }
    
    void cu_memcpy(const double *src, double *dst, int size) {
        cu_memcpy<double, double>(src, dst, size);
    }
    
    void cu_memcpy(const double *src, float *dst, int size) {
        cu_memcpy<double, float>(src, dst, size);
    }
    
    void cu_memcpy(const float *src, double *dst, int size) {
        cu_memcpy<float, double>(src, dst, size);
    }
    
    __global__ void reduction_kernel(GPUtype* input, GPUReduceAccumType* output, int N) {
        extern __shared__ GPUReduceAccumType shared_data[];

        int idx = threadIdx.x + blockIdx.x * blockDim.x; // global index
        int tid = threadIdx.x; // thread's local id within the block (ranging from 0 to blockDim.x - 1).

        // Load elements into shared memory
        if (idx < N) {
            shared_data[tid] = static_cast<GPUReduceAccumType>(input[idx]);
        } else {
            shared_data[tid] = 0.0;//f;
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
            // output[blockIdx.x] = static_cast<GPUtype>(shared_data[0]);
            output[blockIdx.x] = shared_data[0];
        }
    }

    void test(GPUtype *xvv){
        for(int i = 0; i < 100; i++){
            cout << xvv[i] << endl;
        }
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

        GPUReduceAccumType *output;
        cudaMallocManaged(&output, (N + 255) / 256 * sizeof(GPUtype));  // One per block
        cudaMemset(output, 0, (N + 255) / 256 * sizeof(GPUtype));

        // Launch the kernel
        // This will ensure that there will be lauched grid_size blocks,
        // each containing 256 threads and using enough shared memory to 
        // store shared_memory_size float values
        reduction_kernel<<<grid_size, block_size, shared_memory_size>>>(arr, output, N);

        cudaDeviceSynchronize();  // Wait for the kernel to finish

        // Final reduction on the host
        GPUReduceAccumType final_sum = 0.0;
        for (int i = 0; i < grid_size; i++) {
            final_sum += output[i];
        }

        // Cleanup
        cudaFree(output);

        return static_cast<GPUtype>(final_sum);
    }

    // Helper inline functions for FMA in device code.
    __device__ inline float myFMA(float a, float b, float c) {
        return __fmaf_rn(a, b, c);  // FMA in single precision
    }

    __device__ inline double myFMA(double a, double b, double c) {
        return __fma_rn(a, b, c);   // FMA in double precision
    }


    // Kernel to compute dot product with mixed precision.
    // Input arrays x and y are in single precision, but the accumulation is either float or double.
    // It may be possible to use a template this but the shared is problematic.
    __global__ void SDDotKernel(const float* x, const float* y, int n, float* result) {
        extern __shared__ float share_floats[];
        int tid = threadIdx.x;
        int idx = blockIdx.x * blockDim.x + tid;
        float sum = 0.0;
        
        // Each thread processes multiple elements if necessary.
        while (idx < n) {
            // sum += (float)x[idx] * (float)y[idx];
            sum = myFMA(x[idx], y[idx], sum);
            idx += blockDim.x * gridDim.x;
        }
        share_floats[tid] = sum;
        __syncthreads();
        
        // Reduction in shared memory.
        for (int s = blockDim.x / 2; s > 0; s >>= 1) {
            if (tid < s) {
                share_floats[tid] += share_floats[tid + s];
            }
            __syncthreads();
        }
        
        // The first thread in each block atomically adds the block's sum to the global result.
        if (tid == 0) {
            atomicAdd(result, share_floats[0]);
        }
    }

    __global__ void SDDotKernel(const float* x, const float* y, int n, double* result) {
        extern __shared__ double share_doubles[];
        int tid = threadIdx.x;
        int idx = blockIdx.x * blockDim.x + tid;
        double sum = 0.0;
        
        // Each thread processes multiple elements if necessary.
        while (idx < n) {
            // sum += (double)x[idx] * (double)y[idx];
            sum = myFMA((double)x[idx], (double)y[idx], sum);
            idx += blockDim.x * gridDim.x;
        }
        share_doubles[tid] = sum;
        __syncthreads();
        
        // Reduction in shared memory.
        for (int s = blockDim.x / 2; s > 0; s >>= 1) {
            if (tid < s) {
                share_doubles[tid] += share_doubles[tid + s];
            }
            __syncthreads();
        }
        
        // The first thread in each block atomically adds the block's sum to the global result.
        if (tid == 0) {
            atomicAdd(result, share_doubles[0]);
        }
    }

    // Mixed precision dot product: single precision inputs, double precision result.
    template <typename SrcType, typename DstType>
    void cu_SDDot(int n, const SrcType* d_x, const SrcType* d_y, DstType* d_result) {
        int blockSize = 256;
        int gridSize = (n + blockSize - 1) / blockSize;
        
        // Initialize the result to 0.
        cudaMemset(d_result, 0, sizeof(DstType));
        
        // Launch the kernel. Allocate shared memory for block reduction.
        SDDotKernel<<<gridSize, blockSize, blockSize * sizeof(DstType)>>>(d_x, d_y, n, d_result);
        
        // Synchronize to ensure the kernel has completed.
        cudaDeviceSynchronize();
    }
    void cu_SDDot(int n, const float* d_x, const float* d_y, float* d_result) {
        cu_SDDot<float, float>(n, d_x, d_y, d_result);
    }
    void cu_SDDot(int n, const float* d_x, const float* d_y, double* d_result) {
        cu_SDDot<float, double>(n, d_x, d_y, d_result);
    }
    // void cu_SDDot(int n, const double* d_x, const double* d_y, double* d_result) {
    //     cu_SDDot<double, double>(n, d_x, d_y, d_result);
    // }
}