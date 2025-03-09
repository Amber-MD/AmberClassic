#include <stdio.h>
#include <cuda.h> 
#include <cuda_runtime.h>
#include "rism3d_closure_kh.hpp"
#include <cmath>
using namespace std;

namespace rism3d_c {

    /////////////// KERNELS ///////////////
    __global__ void k_set_heaviside(GPUtype *heaviside, GPUtype *huv, int N) {
        int idx = blockIdx.x * blockDim.x + threadIdx.x;

        if (idx < N) {
            // Assign 0 if huv[idx] > 0, otherwise assign 1
            heaviside[idx] = (huv[idx] > 0) ? 0 : 1;
        }
    }

    __global__ void k_guv(GPUtype *uuv, GPUtype *guv, GPUtype *huv, GPUtype *cuv, int size){

        int idx = blockIdx.x * blockDim.x + threadIdx.x;
        GPUtype expoent;

        if(idx < size){
            expoent = -uuv[idx] + huv[idx] - cuv[idx];
            if(expoent >= 0){
                guv[idx] = 1 + expoent;
            }
            else{
                guv[idx] = exp(expoent);
            }
        }
    }

    __global__ void k_excessChemicalPotential_grid(GPUtype *heaviside, GPUtype *huv, GPUtype *cuv, 
                                                   GPUtype *mu_grid, int Nx, int Ny, int Nz,
                                                   GPUtype density, GPUtype voxelVolume){
        int gid = blockIdx.x * blockDim.x + threadIdx.x;

        if(gid < Nx * Ny * Nz){
            int idx = gid / (Ny * Nz);
            int idy = (gid / Nz) % Ny;
            int idz = gid % Nz;

            int id = idx * Ny * (Nz + 2) + idy * (Nz + 2) + idz;

            mu_grid[gid] = (0.5 * huv[id] * huv[id] * heaviside[id] - 0.5 * huv[id] * cuv[id] - cuv[id]) * density * voxelVolume;
        }
    }

    __global__ void k_excessChemicalPotential_lr_grid(GPUtype *heaviside, GPUtype *huv, GPUtype *cuv, 
                                                      GPUtype *mu_grid, int Nx, int Ny, int Nz,
                                                      GPUtype density, GPUtype voxelVolume,
                                                      GPUtype solvent_charge, GPUtype* dcfLongRangeAsympR,
                                                      GPUtype solute_totalcharge, GPUtype solvent_charge_sp,
                                                      GPUtype* tcfLongRangeAsympR){
        int gid = blockIdx.x * blockDim.x + threadIdx.x;

        if(gid < Nx * Ny * Nz){
            int idx = gid / (Ny * Nz);
            int idy = (gid / Nz) % Ny;
            int idz = gid % Nz;

            int id = idx * Ny * (Nz + 2) + idy * (Nz + 2) + idz;
            int id_nonpadded = idx * Ny * Nz + idy * Nz + idz;

            // mu_grid[gid] = ((0.5 * huv[id] * huv[id] * heaviside[id] - 0.5 * huv[id] * cuv[id] - cuv[id]) ) * density * voxelVolume;

            GPUtype cuvlr = solvent_charge * dcfLongRangeAsympR[id_nonpadded];
            GPUtype huvlr = solvent_charge_sp * tcfLongRangeAsympR[id_nonpadded];
            
            if(solute_totalcharge * solvent_charge_sp <= 0.0){
                mu_grid[gid] = ((0.5 * huv[id] * huv[id] * heaviside[id] - 0.5 * huv[id] * cuv[id] - cuv[id]) + 0.5 * cuvlr * huvlr) * density * voxelVolume;
            }
            else{
                mu_grid[gid] = ((0.5 * huv[id] * huv[id] * heaviside[id] - 0.5 * huv[id] * cuv[id] - cuv[id]) + 0.5 * huvlr * (cuvlr - huvlr)) * density * voxelVolume;
            }
            
        }
    }

    __global__ void k_solvPotEne_grid(GPUtype *guv, GPUtype *uuv, GPUtype *solvPotEne_grid, 
                                      int Nx, int Ny, int Nz,
                                      GPUtype density, GPUtype voxelVolume){
        int gid = blockIdx.x * blockDim.x + threadIdx.x;

        if(gid < Nx * Ny * Nz){
            int idx = gid / (Ny * Nz);
            int idy = (gid / Nz) % Ny;
            int idz = gid % Nz;

            int id = idx * Ny * (Nz + 2) + idy * (Nz + 2) + idz;

            solvPotEne_grid[gid] = (guv[id] * uuv[id]) * density * voxelVolume;
        }
    }

    __global__ void k_akirkwoodBuff(GPUtype* huv, GPUtype solvent_charge_sp, 
                                    GPUtype* tcfLongRangeAsympR, 
                                    int Nx, int Ny, int Nz,
                                    GPUtype* block_kb) {

        extern __shared__ GPUtype shared_data[];

        int idx = blockIdx.x * blockDim.x + threadIdx.x; 
        int idy = blockIdx.y * blockDim.y + threadIdx.y; 
        int idz = blockIdx.z * blockDim.z + threadIdx.z; 

        int tid = threadIdx.x * (blockDim.y * blockDim.z) + threadIdx.y * blockDim.z + threadIdx.z; // Flattened thread ID in the block

        int global_index = (idx * Ny * Nz) + (idy * Nz) + idz; // non-padded indexing
        int global_padded_index = (idx * Ny * (Nz + 2)) + (idy * (Nz + 2)) + idz; // padded indexing

        // Load elements into shared memory
        if (idx < Nx && idy < Ny && idz < Nz) {
            shared_data[global_index] = huv[global_padded_index] - solvent_charge_sp * tcfLongRangeAsympR[global_index];
        } else {
            shared_data[tid] = 0.0;
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

        for (int stride = 1; stride < (blockDim.x * blockDim.y * blockDim.z); stride *= 2) {
            if (tid % (2 * stride) == 0 && (tid + stride) < (blockDim.x * blockDim.y * blockDim.z)) {
                shared_data[tid] += shared_data[tid + stride];
            }
            __syncthreads();
        }

        // Write the result of the reduction in the first thread of the block
        // The value of the first thread of each block (which stores the total
        // sum for that block) will be asigned to the output array at the block 
        // id position
        if (tid == 0) {
            int block_id = (blockIdx.x * gridDim.y * gridDim.z) + (blockIdx.y * gridDim.z) + blockIdx.z; // ✅ Row-major block index
            block_kb[block_id] = shared_data[0];
        }
    }

    __global__ void k_akirkwoodBuff2(GPUtype* huv, GPUtype solvent_charge_sp, 
                                     GPUtype* tcfLongRangeAsympR, 
                                     int Nx, int Ny, int Nz, GPUtype* block_kb) {

        extern __shared__ GPUtype shared_data[];

        int gid = threadIdx.x + blockIdx.x * blockDim.x; // global index
        int tid = threadIdx.x; // thread's local id within the block (ranging from 0 to blockDim.x - 1).

        // Load elements into shared memory
        if (gid < Nx * Ny * Nz) {
            int idx = gid / (Ny * Nz);
            int idy = (gid / Nz) % Ny;
            int idz = gid % Nz;

            shared_data[tid] = huv[idx * Ny * (Nz + 2) + idy * (Nz + 2) + idz] - solvent_charge_sp * tcfLongRangeAsympR[gid];
        } else {
            shared_data[tid] = 0.0;
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
            block_kb[blockIdx.x] = shared_data[0];
        }
    }

    /////////////// KERNEL WRAPPER FUNCTIONS ///////////////

    void kh :: cu_guv(GPUtype *uuv, GPUtype *guv, GPUtype *huv, GPUtype *cuv, int size){
        int num_blocks = (size + 255) / 256;
        int num_threads = 256;

        k_guv<<<num_blocks, num_threads>>>(uuv, guv, huv, cuv, size);
        cudaDeviceSynchronize();
    }

    void kh :: cu_set_heaviside(GPUtype *heaviside, GPUtype *huv, int N){
        int num_threads = 256;
        int num_blocks = (N + num_threads - 1) / num_threads;

        k_set_heaviside<<<num_blocks, num_threads>>>(heaviside, huv, N);
    }

    void kh :: cu_excessChemicalPotential_grid(GPUtype *heaviside, GPUtype *huv, GPUtype *cuv, 
                                               GPUtype *mu_grid, int Nx, int Ny, int Nz,
                                               GPUtype density, GPUtype voxelVolume){
        int num_threads = 256;
        int num_blocks = ((Nx * Ny * Nz) + num_threads - 1) / num_threads;

        k_excessChemicalPotential_grid<<<num_blocks, num_threads>>>(heaviside, huv, cuv, 
                                                                    mu_grid, Nx, Ny, Nz,
                                                                    density, voxelVolume);
    }

    void kh :: cu_excessChemicalPotential_lr_grid(GPUtype *heaviside, GPUtype *huv, GPUtype *cuv, 
                                                  GPUtype *mu_grid, int Nx, int Ny, int Nz,
                                                  GPUtype density, GPUtype voxelVolume,
                                                  GPUtype solvent_charge, GPUtype* dcfLongRangeAsympR,
                                                  GPUtype solute_totalcharge, GPUtype solvent_charge_sp,
                                                  GPUtype* tcfLongRangeAsympR){
        int num_threads = 256;
        int num_blocks = ((Nx * Ny * Nz) + num_threads - 1) / num_threads;

        k_excessChemicalPotential_lr_grid<<<num_blocks, num_threads>>>(heaviside, huv, cuv, 
                                                                       mu_grid, Nx, Ny, Nz,
                                                                       density, voxelVolume,
                                                                       solvent_charge, dcfLongRangeAsympR,
                                                                       solute_totalcharge, solvent_charge_sp,
                                                                       tcfLongRangeAsympR);
    }

    void kh :: cu_solvPotEne_grid(GPUtype *guv, GPUtype *uuv, GPUtype *solvPotEne_grid, 
                                  int Nx, int Ny, int Nz,
                                  GPUtype density, GPUtype voxelVolume){
        int num_threads = 256;
        int num_blocks = ((Nx * Ny * Nz) + num_threads - 1) / num_threads;

        k_solvPotEne_grid<<<num_blocks, num_threads>>>(guv, uuv, solvPotEne_grid,
                                                       Nx, Ny, Nz,
                                                       density, voxelVolume);
    }

    GPUtype kh :: cu_akirkwoodBuff(GPUtype* huv, GPUtype solvent_charge_sp, 
                                   GPUtype* tcfLongRangeAsympR, 
                                   int Nx, int Ny, int Nz){

        int N = Nx * Ny * Nz;

        int block_size = 256; // number of threads in each block
        int grid_size = (N + block_size - 1) / block_size; // total number of blocks to be used
        int shared_memory_size = block_size * sizeof(GPUtype);

        GPUtype *block_kb;
        cudaMallocManaged(&block_kb, (N + 255) / 256 * sizeof(GPUtype));

        k_akirkwoodBuff2<<<grid_size, block_size, shared_memory_size>>>(huv, solvent_charge_sp, 
                                                                        tcfLongRangeAsympR, 
                                                                        Nx, Ny, Nz, block_kb);

        cudaDeviceSynchronize();

        GPUtype akb_iv = 0;
        for (int i = 0; i < grid_size; i++) {
            akb_iv += block_kb[i];
        }

        cudaFree(block_kb);
        return akb_iv;
    }

}