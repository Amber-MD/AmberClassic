#include <stdio.h>
#include <cuda.h> 
#include <cuda_runtime.h>
#include "rism3d.hpp"
using namespace std;

namespace rism3d_c{

    // Copy padded data into temporary array
    __global__ void getPadded(GPUtype *data, GPUtype *paddings, int Nx, int Ny, int Nz){
        int gid = blockIdx.x * blockDim.x + threadIdx.x;

        if(gid < Nx * Ny * 2){
            int idx = gid / (Ny * 2);   // X index
            int idy = (gid / 2) % Ny;   // Y index
            int idz = gid % 2;          // Z index

            // Checking indexes
            // printf("idx (%d), idy (%d), idz (%d)\n", idx, idy, idz);

            int padding_idx = gid;
            int guv_padded_idx = idx * Ny * (Nz+2) + idy * (Nz+2) + Nz + idz;
            paddings[padding_idx] = data[guv_padded_idx];
        }

    }

    // Copy data trasnlating its positions into temporary array
    __global__ void translate_data_nr(GPUtype *data, GPUtype *data_temp, int Nx, int Ny, int Nz){
        int gid = blockIdx.x * blockDim.x + threadIdx.x; // x-axis

        if(gid < Nx * Ny * Nz){
            int idx = gid / (Ny * Nz);
            int idy = (gid / Nz) % Ny;
            int idz = gid % Nz;

            // Checking indexes
            // printf("idx (%d), idy (%d), idz (%d)\n", idx, idy, idz);

            int fftw_id = idx * Ny * (Nz + 2) + idy * (Nz + 2) + idz;

            int nr_id = gid;

            data_temp[nr_id] = data[fftw_id];
        }
    }

    __global__ void translate_data_fftw(GPUtype *data, GPUtype *data_temp, int Nx, int Ny, int Nz){
        int gid = blockIdx.x * blockDim.x + threadIdx.x; // x-axis

        if(gid < Nx * Ny * Nz){
            int idx = gid / (Ny * Nz);
            int idy = (gid / Nz) % Ny;
            int idz = gid % Nz;

            // Checking indexes
            // printf("idx (%d), idy (%d), idz (%d)\n", idx, idy, idz);

            int fftw_id = idx * Ny * (Nz + 2) + idy * (Nz + 2) + idz;

            int nr_id = gid;

            data_temp[fftw_id] = data[nr_id];
        }
    }

    __global__ void cpy_padded(GPUtype *data_temp, GPUtype *paddings, int Nx, int Ny, int Nz){
        int gid = blockIdx.x * blockDim.x + threadIdx.x; // x-axis

        if(gid < Nx * Ny * 2){
            int idx = gid / (Ny * 2);
            int idy = (gid / 2) % Ny;
            int idz = gid % 2;

            // Checking indexes
            // printf("idx (%d), idy (%d), idz (%d)\n", idx, idy, idz);

            int padding_idx = gid;
            int guv_padded_idx = idx * Ny * (Nz+2) + idy * (Nz+2) + Nz + idz;

            data_temp[guv_padded_idx] = paddings[padding_idx];
        }
    }

    __global__ void set_dcf_longrange(GPUtype* cuv, GPUtype solvent_charge, GPUtype* dcfLongRangeAsympR,
                                      int Nx, int Ny, int Nz){
        
        int igx = blockIdx.x * blockDim.x + threadIdx.x; 
        int igy = blockIdx.y * blockDim.y + threadIdx.y; 
        int igz = blockIdx.z * blockDim.z + threadIdx.z; 
        
        if (igx < Nx && igy < Ny && igz < Nz){
            int ig = igx * Ny * Nz + igy * Nz + igz;
            int padded_ig = igx * Ny * (Nz + 2) + igy * (Nz + 2) + igz;

            cuv[padded_ig] = solvent_charge * dcfLongRangeAsympR[ig];
        }

    }

    __global__ void subtract_dcf_longrange(GPUtype* guv,  GPUtype* cuv, 
                                           GPUtype solvent_charge, GPUtype* dcfLongRangeAsympR,
                                           int Nx, int Ny, int Nz){
        
        int igx = blockIdx.x * blockDim.x + threadIdx.x; 
        int igy = blockIdx.y * blockDim.y + threadIdx.y; 
        int igz = blockIdx.z * blockDim.z + threadIdx.z; 
        
        if (igx < Nx && igy < Ny && igz < Nz){
            int ig = igx * Ny * Nz + igy * Nz + igz;
            int padded_ig = igx * Ny * (Nz + 2) + igy * (Nz + 2) + igz;

            guv[padded_ig] = cuv[padded_ig] - solvent_charge * dcfLongRangeAsympR[ig];
        }

    }

    __global__ void add_dcf_longrange(GPUtype* guv, GPUtype solvent_charge, 
                                      GPUtype* dcfLongRangeAsympK,
                                      int N){
        
        int ig = blockIdx.x * blockDim.x + threadIdx.x; 
        
        if (ig < N){
            guv[ig] = guv[ig] - solvent_charge * dcfLongRangeAsympK[ig];
        }

    }

    __global__ void subtract_tcf_longrange(GPUtype* huv, GPUtype solvent_charge_sp, 
                                           GPUtype* tcfLongRangeAsympK,
                                           int N){
        
        // As in the fortran code, we are skipping the real and imaginary part for
        // k = 0
        int ig = blockIdx.x * blockDim.x + threadIdx.x + 2; 
        
        if (ig < N){
            huv[ig] = huv[ig] + solvent_charge_sp * tcfLongRangeAsympK[ig];
        }

    }

    void rism3d :: set_dcf_longrange_cu(GPUtype* cuv, GPUtype* solvent_charge, GPUtype* dcfLongRangeAsympR){
        int Nx = grid.globalDimsR[0];
        int Ny = grid.globalDimsR[1];
        int Nz = grid.globalDimsR[2];

        dim3 blockDim(8, 8, 8);
        dim3 gridDim(
            (Nx + blockDim.x - 1) / blockDim.x,
            (Ny + blockDim.y - 1) / blockDim.y,
            (Nz + blockDim.z - 1) / blockDim.z
        );

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            set_dcf_longrange<<<gridDim, blockDim>>>(cuv + iv * Nx * Ny * (Nz + 2), 
                                                     solvent_charge[iv], dcfLongRangeAsympR,
                                                     Nx, Ny, Nz);
        }
        cudaDeviceSynchronize();

    }

    void rism3d :: subtract_dcf_longrange_cu(GPUtype* guv, GPUtype* cuv, GPUtype* solvent_charge, GPUtype* dcfLongRangeAsympR){
        int Nx = grid.globalDimsR[0];
        int Ny = grid.globalDimsR[1];
        int Nz = grid.globalDimsR[2];

        dim3 blockDim(8, 8, 8);
        dim3 gridDim(
            (Nx + blockDim.x - 1) / blockDim.x,
            (Ny + blockDim.y - 1) / blockDim.y,
            (Nz + blockDim.z - 1) / blockDim.z
        );

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            subtract_dcf_longrange<<<gridDim, blockDim>>>(guv + iv * Nx * Ny * (Nz + 2),
                                                          cuv + iv * Nx * Ny * (Nz + 2), 
                                                          solvent_charge[iv], dcfLongRangeAsympR,
                                                          Nx, Ny, Nz);
        }
        cudaDeviceSynchronize();

    }

    __global__ void add_tcf_longrange(GPUtype* huv,  GPUtype solvent_charge_sp, 
                                      GPUtype* tcfLongRangeAsympR,
                                      int Nx, int Ny, int Nz){
        
        int igx = blockIdx.x * blockDim.x + threadIdx.x; 
        int igy = blockIdx.y * blockDim.y + threadIdx.y; 
        int igz = blockIdx.z * blockDim.z + threadIdx.z; 
        
        if (igx < Nx && igy < Ny && igz < Nz){
            int ig = igx * Ny * Nz + igy * Nz + igz;
            int padded_ig = igx * Ny * (Nz + 2) + igy * (Nz + 2) + igz;

            huv[padded_ig] = huv[padded_ig] + solvent_charge_sp * tcfLongRangeAsympR[ig];
        }

    }

    void rism3d :: add_dcf_longrange_cu(GPUtype* guv, GPUtype* solvent_charge, GPUtype* dcfLongRangeAsympK){
        int totalLocalPointsK = grid.totalLocalPointsK;

        int num_threads = 256;
        int num_blocks = (totalLocalPointsK + num_threads-1) / num_threads;

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            add_dcf_longrange<<<num_blocks, num_threads>>>(guv + iv * totalLocalPointsK,
                                                           solvent_charge[iv], dcfLongRangeAsympK,
                                                           totalLocalPointsK);
        }
        cudaDeviceSynchronize();

    }

    void rism3d :: subtract_tcf_longrange_cu(GPUtype* huv, GPUtype* solvent_charge_sp, GPUtype* tcfLongRangeAsympK){
        int totalLocalPointsK = grid.totalLocalPointsK;

        int num_threads = 256;
        int num_blocks = (totalLocalPointsK + num_threads-1) / num_threads;

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            subtract_tcf_longrange<<<num_blocks, num_threads>>>(huv + iv * totalLocalPointsK,
                                                                solvent_charge_sp[iv], tcfLongRangeAsympK,
                                                                totalLocalPointsK);
        }
        cudaDeviceSynchronize();

    }

    void rism3d :: add_tcf_longrange_cu(GPUtype* huv, GPUtype* solvent_charge_sp, GPUtype* tcfLongRangeAsympR){
        int Nx = grid.globalDimsR[0];
        int Ny = grid.globalDimsR[1];
        int Nz = grid.globalDimsR[2];

        dim3 blockDim(8, 8, 8);
        dim3 gridDim(
            (Nx + blockDim.x - 1) / blockDim.x,
            (Ny + blockDim.y - 1) / blockDim.y,
            (Nz + blockDim.z - 1) / blockDim.z
        );

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            add_tcf_longrange<<<gridDim, blockDim>>>(huv + iv * Nx * Ny * (Nz + 2),
                                                     solvent_charge_sp[iv], tcfLongRangeAsympR,
                                                     Nx, Ny, Nz);
        }
        cudaDeviceSynchronize();

    }

    void rism3d :: convert2nr_cu(){
        // create array to store non-padded data
        array_class<GPUtype> guv_temp(array_class<GPUtype>::ROW_MAJOR);
        guv_temp.set_memalloc(&memalloc, true);
        guv_temp.alloc_mem(solventclass.numAtomTypes, grid.globalDimsR[0], grid.globalDimsR[1], grid.globalDimsR[2]);

        // create array store extra padding results
        array_class<GPUtype> paddings(array_class<GPUtype>::ROW_MAJOR);
        paddings.set_memalloc(&memalloc, true);
        paddings.alloc_mem(solventclass.numAtomTypes, grid.globalDimsR[0], grid.globalDimsR[1], 2);

        // Dimensions
        int Nx = grid.globalDimsR[0];
        int Ny = grid.globalDimsR[1];
        int Nz = grid.globalDimsR[2];
        int Nz_pad = Nz + 2;

        // Define the number of threads and blocks to be used
        int num_threads = 256;
        int num_blocks = ((Nx * Ny * 2) + (num_threads-1)) / num_threads;

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            getPadded<<<num_blocks, num_threads>>>(guv.m_data + iv * Nx * Ny * Nz_pad, 
                                                   paddings.m_data + iv * Nx * Ny * 2, 
                                                   Nx, Ny, Nz);   
            
            // Checking errors while launching kernel
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                cout << "getPadded kernel launch failed: " << cudaGetErrorString(err) << endl;
                abort();
            }
        }

        // Redefine number of blocks to be used in next step
        num_blocks = ((Nx * Ny * Nz) + (num_threads-1)) / num_threads;

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            translate_data_nr<<<num_blocks, num_threads>>>(guv.m_data + iv * Nx * Ny * Nz_pad, 
                                                           guv_temp.m_data + iv * Nx * Ny * Nz, 
                                                           Nx, Ny, Nz);

            // Checking errors while launching kernel
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                cout << "translate_data_nr kernel launch failed: " << cudaGetErrorString(err) << endl;
                abort();
            }
        }

        // Ensure all data is in temporary arrays before copying them back to data array
        cudaDeviceSynchronize();

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            // Here we can use cudaMemcpyAsync to copy everything simultaneously
            cudaMemcpyAsync(guv.m_data + iv * Nx * Ny * (Nz + 2), guv_temp.m_data + iv * Nx * Ny * Nz, Nx * Ny * Nz * sizeof(GPUtype), cudaMemcpyDeviceToDevice);
            cudaMemcpyAsync(guv.m_data + iv * Nx * Ny * (Nz + 2) + Nx * Ny * Nz, paddings.m_data + iv * 2 * Nx * Ny, 2 * Nx * Ny * sizeof(GPUtype), cudaMemcpyDeviceToDevice);
        }

        // Ensure all data was copied back to data array before proceeding with calculations
        cudaDeviceSynchronize();
        
    }

    void rism3d :: convert2fftw_cu(){
        // create array to store non-padded data
        array_class<GPUtype> huv_temp(array_class<GPUtype>::ROW_MAJOR);
        huv_temp.set_memalloc(&memalloc, true);
        huv_temp.alloc_mem(solventclass.numAtomTypes, grid.globalDimsR[0], grid.globalDimsR[1], grid.globalDimsR[2] + 2);

        // create array store extra padding results
        array_class<GPUtype> paddings(array_class<GPUtype>::ROW_MAJOR);
        paddings.set_memalloc(&memalloc, true);
        paddings.alloc_mem(solventclass.numAtomTypes, grid.globalDimsR[0], grid.globalDimsR[1], 2);

        // Dimensions
        int Nx = grid.globalDimsR[0];
        int Ny = grid.globalDimsR[1];
        int Nz = grid.globalDimsR[2];
        // int Nz_pad = Nz + 2;

        // Define the number of threads and blocks to be used
        int num_threads = 256;
        int num_blocks = ((Nx * Ny * Nz) + (num_threads-1)) / num_threads;

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            // copying padded data into paddings
            cudaMemcpyAsync(paddings.m_data + iv * 2 * Nx * Ny, huv.m_data + iv * Nx * Ny * (Nz + 2) + Nx * Ny * Nz, 2 * Nx * Ny * sizeof(GPUtype), cudaMemcpyDeviceToDevice);

            // translating data to be in fftw format
            translate_data_fftw<<<num_blocks,num_threads>>>(huv.m_data + iv * Nx * Ny * (Nz + 2), 
                                                            huv_temp.m_data + iv * Nx * Ny * (Nz + 2), 
                                                            Nx, Ny, Nz);
            
            // checking for errors while lauching kernel
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                cout << "translate_data_fftw kernel launch failed: " << cudaGetErrorString(err) << endl;
                abort();
            }
        }

        // Ensure data is translated and padded data is copied into temporary array before proceeding
        cudaDeviceSynchronize();

        // copying data from temp array into actual data array
        cudaMemcpy(huv.m_data, huv_temp.m_data, solventclass.numAtomTypes * Nx * Ny * (Nz+2) * sizeof(GPUtype), cudaMemcpyDeviceToDevice);
        
        // redefine number of blocks to be used
        num_blocks = ((Nx * Ny * 2) + (num_threads-1)) / num_threads;

        // copying padded data back into data array
        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            cpy_padded<<<num_blocks, num_threads>>>(huv.m_data + iv * Nx * Ny * (Nz + 2), 
                                                    paddings.m_data + iv * Nx * Ny * 2, 
                                                    Nx, Ny, Nz);   

            // checking for errors while lauching kernel
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                cout << "cpy_padded kernel launch failed: " << cudaGetErrorString(err) << endl;
                abort();
            }
        }

        // Ensure all data is back into data array before proceeding with calculations
        cudaDeviceSynchronize();

    }

    __global__ void set_padding2zero_cu(GPUtype *data, int Nx, int Ny, int Nz){
        int gid = blockIdx.x * blockDim.x + threadIdx.x;

        if(gid < Nx * Ny * 2){
            int idx = gid / (Ny * 2);   // X index
            int idy = (gid / 2) % Ny;   // Y index
            int idz = gid % 2;          // Z index

            int guv_padded_idx = idx * Ny * (Nz+2) + idy * (Nz+2) + Nz + idz;
            data[guv_padded_idx] = 0;
        }

    }

    void rism3d :: set_padding2zero(GPUtype* data){
        int num_threads = 256;
        int num_blocks = ((grid.globalDimsR[0] * grid.globalDimsR[1] * 2) + (num_threads-1)) / num_threads;

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            set_padding2zero_cu<<<num_blocks, num_threads>>>(data + iv * grid.globalDimsR[0] * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2), 
                                                   grid.globalDimsR[0], grid.globalDimsR[1], grid.globalDimsR[2]);   
            
            // Checking errors while launching kernel
            cudaError_t err = cudaGetLastError();
            if (err != cudaSuccess) {
                cout << "getPadded kernel launch failed: " << cudaGetErrorString(err) << endl;
                abort();
            }
        }

        // Ensure all padded values are set to zero before proceeding
        cudaDeviceSynchronize();

    }

    __global__ void get_residue_cu(GPUtype *residue_array, GPUtype *new_data, GPUtype *old_data, 
                                   int nat, int Nx, int Ny, int Nz, int Nz_pad){
        int gid = threadIdx.x + blockDim.x * blockIdx.x;

        if(gid < nat * Nx * Ny * Nz){
            int iv = gid / (Nx * Ny * Nz);
            int idx = (gid / (Ny * Nz)) % Nx;
            int idy = (gid / Nz) % Ny;
            int idz = gid % Nz;

            int id = iv * Nx * Ny * Nz_pad + idx * Ny * Nz_pad + idy * Nz_pad + idz;

            residue_array[id] = (new_data[id] - 1) - old_data[id];

        }

    }

    void rism3d :: get_residue(){
        int num_threads = 256;
        int num_blocks = ((solventclass.numAtomTypes * grid.globalDimsR[0] * grid.globalDimsR[1] * grid.globalDimsR[2]) + (num_threads-1)) / num_threads;

        get_residue_cu<<<num_blocks, num_threads>>>(cuvres.m_data, guv.m_data, huv.m_data,
                                                    solventclass.numAtomTypes, grid.globalDimsR[0],
                                                    grid.globalDimsR[1], grid.globalDimsR[2],
                                                    grid.globalDimsR[2] + 2);
        cudaDeviceSynchronize();

    }

    __global__ void get_h_k_cu(GPUtype *huv, GPUtype *guv, GPUtype *xvva, int *waveVectorWaveNumberMap, 
                               int totalLocalPointsK){
        int gid = threadIdx.x + blockIdx.x * blockDim.x;

        if(gid < totalLocalPointsK){
            int iga = waveVectorWaveNumberMap[gid/2];

            huv[gid] = huv[gid] + guv[gid] * xvva[iga];
        }

    }

    void rism3d :: get_h_k(){
        cudaMemset(huv.m_data, 0, solventclass.numAtomTypes*grid.totalLocalPointsK*sizeof(GPUtype));
        
        int num_threads = 256;
        int num_blocks = (grid.totalLocalPointsK + (num_threads-1)) / num_threads;

        for(int iv2 = 0; iv2 < solventclass.numAtomTypes; iv2++){
            for(int iv1 = 0; iv1 < solventclass.numAtomTypes; iv1++){
                get_h_k_cu<<<num_blocks, num_threads>>>(huv.m_data + iv1 * grid.totalLocalPointsK,
                                                        guv.m_data + iv2 * grid.totalLocalPointsK,
                                                        xvva.m_data + iv1 * solventclass.numAtomTypes * grid.waveNumberArraySize + iv2 * grid.waveNumberArraySize,
                                                        grid.waveVectorWaveNumberMap.m_data,
                                                        grid.totalLocalPointsK
                                                        );
            }
            cudaDeviceSynchronize();
        }
        // Do we need this second sybchronization or the one before would handle all synchronization?
        // I guess it would handle every synchronization needed
        // cudaDeviceSynchronize();    
    }

    ///////// Conversion between numerical recipes and fftw layout: CPU version

    void rism3d :: convert2nr(){
        ////////// CHANGING TO NR SCHEME
        // create array to point extra padding results
        array_class<GPUtype> paddings;
        paddings.set_memalloc(&memalloc);
        paddings.alloc_mem(grid.globalDimsR[0], grid.globalDimsR[1], 2);

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            // copy padded values into paddings array
            for(int ix = 0; ix < grid.globalDimsR[0]; ix++){
                for(int iy = 0; iy < grid.globalDimsR[1]; iy++){
                    for(int iz = 0; iz < 2; iz++){
                        paddings.m_data[ix*grid.globalDimsR[1]*2 + iy*2 + iz] = guv.m_data[iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
                                                                        ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
                                                                        iy*(grid.globalDimsR[2] + 2) + grid.globalDimsR[2] + iz];
                    }
                }
            }

            // move non-padded data to be stored contiguously
            for(int ix = 0; ix < grid.globalDimsR[0]; ix++){
                for(int iy = 0; iy < grid.globalDimsR[1]; iy++){
                    
                    int fftw_id = iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
                                    ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
                                    iy*(grid.globalDimsR[2] + 2);

                    // The numerical recipes layout is:
                    // "a contiguous piece of nx*ny*nz memory, with the Nyquist
                    // frequency data (2*ny*nz) starting at the end of this block"
                    // obs: this is true for each atom type. Thus, the numerical recipe index (nr_id)
                    // should still stride considering the paddings when looping over the atom type 
                    // index (iv) to make sure the xyz block starts at the correct positions for
                    // each atom type
                    int nr_id = iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
                                ix*grid.globalDimsR[1]*grid.globalDimsR[2] +
                                iy*grid.globalDimsR[2];
                    
                    for(int iz = 0; iz < grid.globalDimsR[2]; iz++){
                        guv.m_data[nr_id + iz] = guv.m_data[fftw_id + iz];
                    }
                }
            }
        
            // Calculate index where padded values should start
            int end_id = iv*(grid.globalDimsR[0] * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2)) + 
                            grid.globalDimsR[0]*grid.globalDimsR[1]*grid.globalDimsR[2];
            
            // Add padded values at the end of each chunk of memory containing data for each solvent index
            for(int i = 0; i < grid.globalDimsR[0]*grid.globalDimsR[1]*2; i++){
                // cout << "end_id = " << end_id << endl;
                guv.m_data[end_id + i] = paddings.m_data[i];
            }
        }
    }

    void rism3d :: convert2fftw(){

        ////////// CHANGING TO NR SCHEME
        // create array to point extra padding results
        array_class<GPUtype> paddings;
        paddings.set_memalloc(&memalloc);
        paddings.alloc_mem(grid.globalDimsR[0], grid.globalDimsR[1], 2);

        array_class<GPUtype> huv_new;
        huv_new.set_memalloc(&memalloc);
        huv_new.alloc_mem(solventclass.numAtomTypes, grid.globalDimsR[0], grid.globalDimsR[1], grid.globalDimsR[2] + 2);

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            for(int ix = 0; ix < grid.globalDimsR[0]; ix++){
                for(int iy = 0; iy < grid.globalDimsR[1]; iy++){
                    for (int iz = 0; iz < grid.globalDimsR[2]; iz++) {
                        huv_new.m_data[iv * grid.globalDimsR[0] * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2) + 
                                       ix * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2) + 
                                       iy * (grid.globalDimsR[2] + 2) + iz] = 
                                       huv.m_data[iv * grid.globalDimsR[0] * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2) + 
                                                  ix * grid.globalDimsR[1] * grid.globalDimsR[2] + iy * grid.globalDimsR[2] + iz];
                    }
                }
            }

            // Calculate index where padded values starts
            int end_id = iv*(grid.globalDimsR[0] * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2)) + 
                            grid.globalDimsR[0]*grid.globalDimsR[1]*grid.globalDimsR[2];
            
            // Copy padded values into padding array
            for(int i = 0; i < grid.globalDimsR[0]*grid.globalDimsR[1]*2; i++){
                paddings.m_data[i] = huv.m_data[end_id + i];
            }
            
            // Place padded values at the correct places
            for(int ix = 0; ix < grid.globalDimsR[0]; ix++){
                for(int iy = 0; iy < grid.globalDimsR[1]; iy++){
                    for(int iz = 0; iz < 2; iz++){
                        huv_new.m_data[iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
                                       ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
                                       iy*(grid.globalDimsR[2] + 2) + grid.globalDimsR[2] + iz] = 
                                       paddings.m_data[ix*grid.globalDimsR[1]*2 + iy*2 + iz];
                    }
                }
            }
        }

        int size = solventclass.numAtomTypes * grid.globalDimsR[0] * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2);
        copy(huv_new.m_data, huv_new.m_data + size, huv.m_data);
    }



}