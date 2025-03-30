#include <stdio.h>
#include <cuda.h> 
#include <cuda_runtime.h>
#include "rism3d_potential.hpp"
using namespace std;

namespace rism3d_c {

    // Define PI as a constant on device
    __device__ const GPUtype PI = 3.14159265358979323846;

#if RISMCUDA_DOUBLE
    __device__ const GPUtype max_value = DBL_MAX;
#else
    __device__ const GPUtype max_value = FLT_MAX;
#endif // RISMCUDA_DOUBLE

    // Kernel to calculate Lennard-Jones potential
    __global__ void k_potential_calc(GPUtype *lj, GPUtype *A, GPUtype *B, 
                                     int solvent_numAtomTypes, int solute_numAtoms, 
                                     int x_dim, int y_dim, int z_dim, 
                                     GPUtype *pos,
                                     GPUtype gridspc_x, GPUtype gridspc_y, GPUtype gridspc_z){
        GPUtype rx,ry,rz;
        GPUtype dz2,dy2,dx2;
        GPUtype rs2i, rs6i;
        GPUtype r2;

        GPUPotAccumType potential = 0.0;

        // Use long long int to handle large systems
        long long int idx = threadIdx.x + (long long int)blockIdx.x * blockDim.x;

        if (idx < (long long int)solvent_numAtomTypes * x_dim * y_dim * z_dim) {
            int i = idx / (x_dim * y_dim * z_dim); // solvent numAtomTypes index
            int m = (idx / (y_dim * z_dim)) % x_dim; // x index
            int l = (idx / z_dim) % y_dim; // y index
            int k = idx % z_dim; // z index
            for(int solu_n = 0; solu_n < solute_numAtoms; solu_n++){

                int j = solu_n; // solute numAtoms idex
    
                // Checking out-of-bound access
                // int mem_idx = i*x_dim*y_dim*(z_dim + 2) + m*y_dim*(z_dim + 2) + l*(z_dim + 2) + k;
                // int max_size = solvent_numAtomTypes * x_dim * y_dim * (z_dim + 2);
                // if (mem_idx >= max_size) {
                //     printf("Thread %lld attempted out-of-bounds access: %d >= %d\n", idx, mem_idx, max_size);
                // }

                rz = k*gridspc_z;
                // If pos is column major
                // dz2 = (rz - pos[2 + j*3])*(rz - pos[2 + j*3]);
                // If pos is row major
                dz2 = (rz - pos[2*solute_numAtoms + j])*(rz - pos[2*solute_numAtoms + j]);
                
                ry = l*gridspc_y;
                // If pos is column major
                // dy2 = (ry - pos[1 + j*3])*(ry - pos[1 + j*3]);
                // If pos is row major
                dy2 = (ry - pos[1*solute_numAtoms + j])*(ry - pos[1*solute_numAtoms + j]);

                rx = m*gridspc_x;
                // If pos is column major
                // dx2 = (rx - pos[0 + j*3])*(rx - pos[0 + j*3]);
                // If pos is row major
                dx2 = (rx - pos[0*solute_numAtoms + j])*(rx - pos[0*solute_numAtoms + j]);

                r2 = dx2 + dy2 + dz2;

                rs2i = (GPUtype)1.0/r2;
                rs6i = rs2i*rs2i*rs2i;

                potential += rs6i*(rs6i*A[j*solvent_numAtomTypes+i] - B[j*solvent_numAtomTypes+i]);

            }
            // old column major order
            // lj[m + l*x_dim + k*x_dim*y_dim + i*x_dim*y_dim*z_dim] += rs6i*(rs6i*A[j*solvent_numAtomTypes+i] - B[j*solvent_numAtomTypes+i]);

            // add 2 to z_dim in order to account the extra padding not used
            lj[i*x_dim*y_dim*(z_dim + 2) + m*y_dim*(z_dim + 2) + l*(z_dim + 2) + k] += static_cast<GPUtype>(potential);

            if(lj[i*x_dim*y_dim*(z_dim + 2) + m*y_dim*(z_dim + 2) + l*(z_dim + 2) + k] > max_value || isnan(lj[i*x_dim*y_dim*(z_dim + 2) + m*y_dim*(z_dim + 2) + l*(z_dim + 2) + k])){
                lj[i*x_dim*y_dim*(z_dim + 2) + m*y_dim*(z_dim + 2) + l*(z_dim + 2) + k] = sqrt(max_value);
            }

        }
    }

     // Kernel to calculate Coulomb potential
     __global__ void k_coulomb_potential_calc(int target_x_low_ind,  int target_x_high_ind,
        int target_y_low_ind,  int target_y_high_ind,
        int target_z_low_ind,  int target_z_high_ind,
        GPUtype grid_xmin,      GPUtype grid_ymin,      GPUtype grid_zmin,
        GPUtype grid_spacing_x,       GPUtype grid_spacing_y,       GPUtype grid_spacing_z,
        int grid_dim_x,   int grid_dim_y,   int grid_dim_z,
        int solute_numAtoms, int solute_numAtoms_idx_start,
        const GPUtype *solute_position_x, const GPUtype *solute_position_y, const GPUtype *solute_position_z, const GPUtype *solute_charge,
        GPUtype *potential, GPUtype solvent_charge){

        // Remeber: our potential array has (two) extra paddings in the z dimentions
        // Thus, we need to account this to get the correct indexes
        int target_yz_dim = grid_dim_y * (grid_dim_z + 2);

        // Compute 3D thread and block indices
        int ix = blockIdx.x * blockDim.x + threadIdx.x + target_x_low_ind;
        int iy = blockIdx.y * blockDim.y + threadIdx.y + target_y_low_ind;
        int iz = blockIdx.z * blockDim.z + threadIdx.z + target_z_low_ind;

        if (ix > target_x_high_ind || iy > target_y_high_ind || iz > target_z_high_ind) return;

        int ii = (ix * target_yz_dim) + (iy * (grid_dim_z + 2)) + iz;
        GPUPotAccumType temporary_potential = 0.0;

        GPUtype tx = grid_xmin + (ix - target_x_low_ind) * grid_spacing_x;
        GPUtype ty = grid_ymin + (iy - target_y_low_ind) * grid_spacing_y;
        GPUtype tz = grid_zmin + (iz - target_z_low_ind) * grid_spacing_z;

        for (int j = 0; j < solute_numAtoms; j++) {
            int jj = solute_numAtoms_idx_start + j;

            GPUtype dx = tx - solute_position_x[jj];
            GPUtype dy = ty - solute_position_y[jj];
            GPUtype dz = tz - solute_position_z[jj];
            GPUtype r  = sqrt(dx*dx + dy*dy + dz*dz);

            if (r > 0) {
            temporary_potential += static_cast<GPUPotAccumType>(solute_charge[jj] / r);
            }
        }

        // Use atomicAdd to safely update the shared potential array
        // atomicAdd(&potential[ii], temporary_potential);
        // No atomicAdd needed, since each thread handles a unique ii
        potential[ii] = static_cast<GPUtype>(temporary_potential * solvent_charge);

        if(potential[ii] > max_value || isnan(potential[ii])){
            potential[ii] = sqrt(max_value);
        }

    }

    void rism3d_potential :: potential_calc(){
        // TO DO: Check if grid size was changed
        
        // Allocating memory using globalDimsK to account the extra padding necessary for huv array
        uuv.alloc_mem(solventclass_p->numAtomTypes, grid_p->globalDimsK[0], grid_p->globalDimsK[1], grid_p->globalDimsK[2]);
        cudaMemset(uuv.m_data, 0, grid_p->globalDimsK[0]*grid_p->globalDimsK[1]*grid_p->globalDimsK[2]*solventclass_p->numAtomTypes*sizeof(GPUtype));

        // Prefetching arrays that will be used on the device functions
        int device = -1;
        cudaGetDevice(&device);
        cudaMemPrefetchAsync(ljAUV.m_data, solventclass_p->numAtomTypes * soluteclass_p->numAtoms * sizeof(GPUtype), device, NULL);
        cudaMemPrefetchAsync(ljBUV.m_data, solventclass_p->numAtomTypes * soluteclass_p->numAtoms * sizeof(GPUtype), device, NULL);
        cudaMemPrefetchAsync(soluteclass_p->position.m_data, 3 * soluteclass_p->numAtoms * sizeof(GPUtype), device, NULL);
        cudaMemPrefetchAsync(uuv.m_data, solventclass_p->numAtomTypes*grid_p->globalDimsK[0]*grid_p->globalDimsK[1]*grid_p->globalDimsK[2] * sizeof(GPUtype), device, NULL);

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            cout << "Probable error with cudaMemPrefetchAsync" << cudaGetErrorString(err) << endl;
            abort();
        }

        if(soluteclass_p->charged){
            if(periodic == true){
                cout << "Sorry: periodic systems not supported, yet!" << endl;
                abort();
            }
            else{
                if(treeCoulomb == true){
                    cout << "tree code version not available, yet." << endl;
                    cout << "Run using --notreeCoulomb flag" << endl;
                    abort();
                }
                else{
                    dim3 blockDim(8, 8, 8);
                    dim3 gridDim((grid_p->localDimsR[0] + blockDim.x) / blockDim.x,
                                 (grid_p->localDimsR[1] + blockDim.y) / blockDim.y,
                                 (grid_p->localDimsR[2] + blockDim.z) / blockDim.z);


                    // Iterate over solvent atom types to run calculations simultaneously
                    // and synchronize only after the for loop
                    for(int iv = 0; iv < solventclass_p->numAtomTypes; iv++){
                        k_coulomb_potential_calc<<<gridDim, blockDim>>>(0, grid_p->localDimsR[0] - 1,
                                                                        0, grid_p->localDimsR[1] - 1,
                                                                        0, grid_p->localDimsR[2] - 1,
                                                                        0, 0, 0,
                                                                        grid_p->spacing[0], grid_p->spacing[1], grid_p->spacing[2],
                                                                        grid_p->localDimsR[0], grid_p->localDimsR[1], grid_p->localDimsR[2],
                                                                        soluteclass_p->numAtoms, 0,
                                                                        soluteclass_p->position.m_data, 
                                                                        soluteclass_p->position.m_data + soluteclass_p->numAtoms, 
                                                                        soluteclass_p->position.m_data + 2*soluteclass_p->numAtoms, 
                                                                        soluteclass_p->charge.m_data,
                                                                        uuv.m_data + iv * grid_p->localDimsR[0] * grid_p->localDimsR[1] * (grid_p->localDimsR[2] + 2),
                                                                        solventclass_p->charge.m_data[iv]);
                    }
                    cudaError_t err2 = cudaGetLastError();
                    if (err2 != cudaSuccess) {
                        cout << "k_coulomb_potential_calc kernel launch failed: " << cudaGetErrorString(err2) << endl;
                        abort();
                    }
                    cudaDeviceSynchronize();

                    // Saving values to compare with fortran version: i am keeping this for now
// #if RISMCUDA_DOUBLE
//                     ofstream file1("/home/fcarvalho/rism3d.cuda.test.chg/coulomb_1_db.txt");
//                     ofstream file2("/home/fcarvalho/rism3d.cuda.test.chg/coulomb_2_db.txt");
// #else
//                     ofstream file1("/home/fcarvalho/rism3d.cuda.test.chg/coulomb_1_float.txt");
//                     ofstream file2("/home/fcarvalho/rism3d.cuda.test.chg/coulomb_2_float.txt");
// #endif // RISMCUDA_DOUBLE
//                     file1 << std::scientific << std::setprecision(16);
//                     file2 << std::scientific << std::setprecision(16);
//                     for(int ix = 0; ix < grid_p->localDimsR[0]; ix++){
//                         for(int iy = 0; iy < grid_p->localDimsR[1]; iy++){
//                             for(int iz = 0; iz < grid_p->localDimsR[2]; iz++){
//                                 file1 << uuv.m_data[0 * grid_p->localDimsR[0] * grid_p->localDimsR[1] * (grid_p->localDimsR[2] + 2) + 
//                                                     ix * grid_p->localDimsR[1] * (grid_p->localDimsR[2] + 2) + 
//                                                     iy * (grid_p->localDimsR[2] + 2) + 
//                                                     iz] << endl;
//                                 file2 << uuv.m_data[1 * grid_p->localDimsR[0] * grid_p->localDimsR[1] * (grid_p->localDimsR[2] + 2) + 
//                                                     ix * grid_p->localDimsR[1] * (grid_p->localDimsR[2] + 2) + 
//                                                     iy * (grid_p->localDimsR[2] + 2) + 
//                                                     iz] << endl;
//                             }
//                         }
//                     }
//                     file1.close();
//                     file2.close();
                }
            }
        }

        if(periodic == true){
            cout << "Sorry: periodic systems not supported, yet!" << endl;
            abort();
        } else{
            // using long long type to ensure we will not get overflow while computing num_blocks
            long long num_blocks = (static_cast<long long>(solventclass_p->numAtomTypes) * static_cast<long long>(soluteclass_p->numAtoms) * 
                                    static_cast<long long>(grid_p->globalDimsR[0]) * static_cast<long long>(grid_p->globalDimsR[1]) * 
                                    static_cast<long long>(grid_p->globalDimsR[2]) + 255) / 256;
            int num_threads = 256;

            // Checking GPU kernel properties (uncomment below to see values)
            // cudaDeviceProp prop;
            // cudaGetDeviceProperties(&prop, 0);
            // cout << prop.maxGridSize[0] << " " << prop.maxGridSize[1] << " " << prop.maxGridSize[2] << endl;

            k_potential_calc<<<num_blocks, num_threads>>>(uuv.m_data, ljAUV.m_data, ljBUV.m_data, 
                                                          solventclass_p->numAtomTypes, soluteclass_p->numAtoms, 
                                                          grid_p->globalDimsR[0], grid_p->globalDimsR[1], grid_p->globalDimsR[2], 
                                                          soluteclass_p->position.m_data,
                                                          grid_p->spacing[0], grid_p->spacing[1], grid_p->spacing[2]);
            
            cudaDeviceSynchronize();
            cudaError_t err3 = cudaGetLastError();
            if (err3 != cudaSuccess) {
                cout << "k_potential_calc kernel launch failed: " << cudaGetErrorString(err3) << endl;
                abort();
            }

        }

        // Saving values for full potential to compare with Fortran version: i am keeping this for now
// #if RISMCUDA_DOUBLE
//         ofstream file1("/home/fcarvalho/rism3d.cuda.test.chg/full_pot_1_db.txt");
//         ofstream file2("/home/fcarvalho/rism3d.cuda.test.chg/full_pot_2_db.txt");
// #else
//         ofstream file1("/home/fcarvalho/rism3d.cuda.test.chg/full_pot_1_float.txt");
//         ofstream file2("/home/fcarvalho/rism3d.cuda.test.chg/full_pot_2_float.txt");
// #endif // RISMCUDA_DOUBLE
//         file1 << std::scientific << std::setprecision(16);
//         file2 << std::scientific << std::setprecision(16);
//         for(int ix = 0; ix < grid_p->localDimsR[0]; ix++){
//             for(int iy = 0; iy < grid_p->localDimsR[1]; iy++){
//                 for(int iz = 0; iz < grid_p->localDimsR[2]; iz++){
//                     file1 << uuv.m_data[0 * grid_p->localDimsR[0] * grid_p->localDimsR[1] * (grid_p->localDimsR[2] + 2) + 
//                                         ix * grid_p->localDimsR[1] * (grid_p->localDimsR[2] + 2) + 
//                                         iy * (grid_p->localDimsR[2] + 2) + 
//                                         iz] << endl;
//                     file2 << uuv.m_data[1 * grid_p->localDimsR[0] * grid_p->localDimsR[1] * (grid_p->localDimsR[2] + 2) + 
//                                         ix * grid_p->localDimsR[1] * (grid_p->localDimsR[2] + 2) + 
//                                         iy * (grid_p->localDimsR[2] + 2) + 
//                                         iz] << endl;
//                 }
//             }
//         }
//         file1.close();
//         file2.close();
        
    }

    __global__ void K_dcf_long_range_asymptotics_R(int target_x_low_ind,    int target_x_high_ind,
                                                   int target_y_low_ind,    int target_y_high_ind,
                                                   int target_z_low_ind,    int target_z_high_ind,
                                                   GPUtype grid_xmin,      GPUtype grid_ymin,      GPUtype grid_zmin,
                                                   GPUtype grid_spacing_x,       GPUtype grid_spacing_y,       GPUtype grid_spacing_z,
                                                   int grid_dim_x,   int grid_dim_y,   int grid_dim_z,
                                                   int solute_numAtoms, int solute_numAtoms_idx_start,
                                                   const GPUtype *solute_position_x, const GPUtype *solute_position_y, const GPUtype *solute_position_z, const GPUtype *solute_charge,
                                                   GPUtype eta, GPUtype *dcf_long_range_asymptotics)
    {
        // Compute global indices for this thread
        int ix = blockIdx.x * blockDim.x + threadIdx.x + target_x_low_ind;
        int iy = blockIdx.y * blockDim.y + threadIdx.y + target_y_low_ind;
        int iz = blockIdx.z * blockDim.z + threadIdx.z + target_z_low_ind;

        // Check bounds
        if (ix > target_x_high_ind || iy > target_y_high_ind || iz > target_z_high_ind) {
            return;
        }

        int target_yz_dim = grid_dim_y * grid_dim_z;
        int ii = (ix * target_yz_dim) + (iy * grid_dim_z) + iz;

        GPUtype tx = grid_xmin + (ix - target_x_low_ind) * grid_spacing_x;
        GPUtype ty = grid_ymin + (iy - target_y_low_ind) * grid_spacing_y;
        GPUtype tz = grid_zmin + (iz - target_z_low_ind) * grid_spacing_z;

        GPUPotAccumType temporary_dcf_long_range_asymptotics = 0.0;

        // Each thread processes all source points for one (ix, iy, iz)
        for (int j = 0; j < solute_numAtoms; j++) {
            int jj = solute_numAtoms_idx_start + j;
            GPUtype dx = tx - solute_position_x[jj];
            GPUtype dy = ty - solute_position_y[jj];
            GPUtype dz = tz - solute_position_z[jj];
            GPUtype r  = sqrt(dx*dx + dy*dy + dz*dz);

            if (r > 0) {
                temporary_dcf_long_range_asymptotics -= static_cast<GPUPotAccumType>(solute_charge[jj] * erf(r / eta) / r);
            }
            else{
                temporary_dcf_long_range_asymptotics -= static_cast<GPUPotAccumType>(solute_charge[jj] / (sqrt(PI) * eta) * 2.0);
            }
        }

        // No atomicAdd needed, since each thread handles a unique ii
        dcf_long_range_asymptotics[ii] = static_cast<GPUtype>(temporary_dcf_long_range_asymptotics);
    }

    void rism3d_potential :: dcf_long_range_asymptotics_R(int target_x_low_ind, int target_x_high_ind,
                                                          int target_y_low_ind, int target_y_high_ind,
                                                          int target_z_low_ind, int target_z_high_ind,
                                                          GPUtype grid_xmin, GPUtype grid_ymin, GPUtype grid_zmin,
                                                          GPUtype grid_spacing_x, GPUtype grid_spacing_y, GPUtype grid_spacing_z,
                                                          int grid_dim_x, int grid_dim_y, int grid_dim_z,
                                                          int solute_numAtoms, int solute_numAtoms_idx_start,  
                                                          GPUtype *solute_position_x, GPUtype *solute_position_y, GPUtype *solute_position_z, GPUtype *solute_charge,
                                                          GPUtype eta, GPUtype *dcf_long_range_asymptotics){
    
        // Compute grid and block dimensions
        int x_range = target_x_high_ind - target_x_low_ind + 1;
        int y_range = target_y_high_ind - target_y_low_ind + 1;
        int z_range = target_z_high_ind - target_z_low_ind + 1;

        // Choose block dimensions (tune as necessary)
        dim3 blockDim(8, 8, 8);
        dim3 gridDim((x_range + blockDim.x - 1) / blockDim.x,
                     (y_range + blockDim.y - 1) / blockDim.y,
                     (z_range + blockDim.z - 1) / blockDim.z);

        // Choose a CUDA stream if gpu_async_stream_id corresponds to a valid cudaStream_t.
        // Here we assume you have a pre-created cudaStream_t array or handle mapping to gpu_async_stream_id.
        // If not using streams, just use the default stream (0).
        cudaStream_t stream = 0; // Replace with appropriate stream handle if available

        // Launch the CUDA kernel
        K_dcf_long_range_asymptotics_R<<<gridDim, blockDim, 0, stream>>>(target_x_low_ind, target_x_high_ind,
                                                                         target_y_low_ind, target_y_high_ind,
                                                                         target_z_low_ind, target_z_high_ind,
                                                                         grid_xmin, grid_ymin, grid_zmin,
                                                                         grid_spacing_x, grid_spacing_y, grid_spacing_z,
                                                                         grid_dim_x, grid_dim_y, grid_dim_z,
                                                                         solute_numAtoms, solute_numAtoms_idx_start,
                                                                         solute_position_x, solute_position_y, solute_position_z, solute_charge,
                                                                         eta, dcf_long_range_asymptotics);

        // Wait for the kernel to finish here
        // cudaStreamSynchronize(stream);
        cudaDeviceSynchronize();
    }

    __global__ void K_tcf_long_range_asymptotics_R(int target_x_low_ind, int target_x_high_ind,
                                                   int target_y_low_ind, int target_y_high_ind,
                                                   int target_z_low_ind, int target_z_high_ind,
                                                   GPUtype grid_xmin, GPUtype grid_ymin, GPUtype grid_zmin,
                                                   GPUtype grid_spacing_x, GPUtype grid_spacing_y, GPUtype grid_spacing_z,
                                                   int grid_dim_x, int grid_dim_y, int grid_dim_z,
                                                   int solute_numAtoms, int solute_numAtoms_idx_start,  
                                                   GPUtype *solute_position_x, GPUtype *solute_position_y, GPUtype *solute_position_z, GPUtype *solute_charge,
                                                   GPUtype kap, GPUtype eta, GPUtype solvent_dielconst, GPUtype *tcf_long_range_asymptotics){
        
        // Calculate 3D grid indices
        int ix = blockIdx.x * blockDim.x + threadIdx.x + target_x_low_ind;
        int iy = blockIdx.y * blockDim.y + threadIdx.y + target_y_low_ind;
        int iz = blockIdx.z * blockDim.z + threadIdx.z + target_z_low_ind;

        if (ix > target_x_high_ind || iy > target_y_high_ind || iz > target_z_high_ind)
            return; // Out-of-bounds check

        int target_yz_dim = grid_dim_y * grid_dim_z;
        GPUtype kap_eta_2 = kap * eta / 2.0;

        // Compute target index in potential array
        int ii = (ix * target_yz_dim) + (iy * grid_dim_z) + iz;

        // Compute target coordinates
        GPUtype tx = grid_xmin + (ix - target_x_low_ind) * grid_spacing_x;
        GPUtype ty = grid_ymin + (iy - target_y_low_ind) * grid_spacing_y;
        GPUtype tz = grid_zmin + (iz - target_z_low_ind) * grid_spacing_z;

        GPUPotAccumType temporary_tcf_long_range_asymptotics = 0.0;

        // Loop over sources
        for (int j = 0; j < solute_numAtoms; j++) {
            int jj = solute_numAtoms_idx_start + j;
            GPUtype dx = tx - solute_position_x[jj];
            GPUtype dy = ty - solute_position_y[jj];
            GPUtype dz = tz - solute_position_z[jj];
            GPUtype r  = sqrt(dx * dx + dy * dy + dz * dz);

            GPUtype kap_r = kap * r;
            GPUtype r_eta = r / eta;

            if (r > 0) {
                temporary_tcf_long_range_asymptotics += static_cast<GPUPotAccumType>(- solute_charge[jj] / r
                                                        * (exp(-kap_r) * erfc(kap_eta_2 - r_eta)
                                                        -  exp( kap_r) * erfc(kap_eta_2 + r_eta)) / 2);
            }
            else{
                temporary_tcf_long_range_asymptotics += static_cast<GPUPotAccumType>(- solute_charge[jj]
                                                        * (2 / (sqrt(PI) * eta)
                                                        -  exp( kap_eta_2 * kap_eta_2) * kap * erfc(kap_eta_2)) / exp( kap_eta_2 * kap_eta_2));
                }
        }

        // Atomic update to prevent race conditions
        // atmicAdd(&potential[ii], temporary_tcf_long_range_asymptotics);
        // No atomicAdd needed, since each thread handles a unique ii
        tcf_long_range_asymptotics[ii] = static_cast<GPUtype>(temporary_tcf_long_range_asymptotics * exp(kap_eta_2 * kap_eta_2) / solvent_dielconst);
    }

    void rism3d_potential :: tcf_long_range_asymptotics_R(int target_x_low_ind, int target_x_high_ind,
                                                          int target_y_low_ind, int target_y_high_ind,
                                                          int target_z_low_ind, int target_z_high_ind,
                                                          GPUtype grid_xmin, GPUtype grid_ymin, GPUtype grid_zmin,
                                                          GPUtype grid_spacing_x, GPUtype grid_spacing_y, GPUtype grid_spacing_z,
                                                          int grid_dim_x, int grid_dim_y, int grid_dim_z,
                                                          int solute_numAtoms, int solute_numAtoms_idx_start,  
                                                          GPUtype *solute_position_x, GPUtype *solute_position_y, GPUtype *solute_position_z, GPUtype *solute_charge,
                                                          GPUtype xappa, GPUtype chargeSmear, GPUtype solvent_dielconst, GPUtype *tcf_long_range_asymptotics){

        // Compute grid and block dimensions
        int x_range = target_x_high_ind - target_x_low_ind + 1;
        int y_range = target_y_high_ind - target_y_low_ind + 1;
        int z_range = target_z_high_ind - target_z_low_ind + 1;

        // Choose block dimensions (tune as necessary)
        dim3 blockDim(8, 8, 8);
        dim3 gridDim((x_range + blockDim.x - 1) / blockDim.x,
                     (y_range + blockDim.y - 1) / blockDim.y,
                     (z_range + blockDim.z - 1) / blockDim.z);

        // Choose a CUDA stream if gpu_async_stream_id corresponds to a valid cudaStream_t.
        // Here we assume you have a pre-created cudaStream_t array or handle mapping to gpu_async_stream_id.
        // If not using streams, just use the default stream (0).
        cudaStream_t stream = 0; // Replace with appropriate stream handle if available

        // Launch the CUDA kernel
        K_tcf_long_range_asymptotics_R<<<gridDim, blockDim, 0, stream>>>(target_x_low_ind, target_x_high_ind,
                                                                         target_y_low_ind, target_y_high_ind,
                                                                         target_z_low_ind, target_z_high_ind,
                                                                         grid_xmin, grid_ymin, grid_zmin,
                                                                         grid_spacing_x, grid_spacing_y, grid_spacing_z,
                                                                         grid_dim_x, grid_dim_y, grid_dim_z,
                                                                         solute_numAtoms, solute_numAtoms_idx_start,
                                                                         solute_position_x, solute_position_y, solute_position_z, solute_charge,
                                                                         xappa, chargeSmear, solvent_dielconst, tcf_long_range_asymptotics);

        // Wait for the kernel to finish here
        // cudaStreamSynchronize(stream);
        cudaDeviceSynchronize();

    }

    __global__ void K_calc_dcf_tcf_LongRangeAsympK(GPUtype* waveVectorX, // [numWaveVectors]
                                     GPUtype* waveVectorY, // [numWaveVectors]
                                     GPUtype* waveVectorZ, // [numWaveVectors]
                                     GPUtype* waveVectors2,
                                     GPUtype cut2_chlk,
                                     GPUtype* position,    // [3 * numAtoms] (x, y, z for all atoms)
                                     GPUtype* charge,      // [numAtoms] (charges of the solute atoms)
                                     GPUtype* dcfLongRangeAsympK,
                                     GPUtype asympk_const,
                                     GPUtype smear2_4,
                                     int numAtoms,
                                     int numWaveVectors_2,
                                     int start_ind,
                                     bool ionic, 
                                     GPUtype* tcfLongRangeAsympK,
                                     GPUtype xappa2, 
                                     GPUtype solvent_dielconst)
    {
        int ig = blockIdx.x * blockDim.x + threadIdx.x + start_ind;
        if (ig >= numWaveVectors_2) return; // Out of bounds check

        if(waveVectors2[ig] > cut2_chlk){
            return;
        }
        else{
            // Extract wave vector components for this thread
            GPUtype kx = waveVectorX[ig];
            GPUtype ky = waveVectorY[ig];
            GPUtype kz = waveVectorZ[ig];

            GPUPotAccumType sumCos = 0;
            GPUPotAccumType sumSin = 0;

            for (int iu = 0; iu < numAtoms; ++iu) {
                // Compute indices for x, y, z components of the atom
                GPUtype x = position[iu];
                GPUtype y = position[iu + numAtoms];
                GPUtype z = position[iu + 2 * numAtoms];

                // Compute the phase: dot product of wave vector and atom position
                GPUtype phase = kx * x + ky * y + kz * z;

                // Accumulate cosine and sine contributions
                sumCos = sumCos + static_cast<GPUPotAccumType>(charge[iu] * cos(phase));

                // Because we are using half of z axis for the in-place FFTW 
                // our CUDA version needs to use the 
                // complex conjugate of exp(i*dot(k,R))
                // (i.e. using a minus sign for the sine function)
                // to get the correct result.
                // There is a phase shift compared to Fortran version,
                // which takes half of x axis for the FFTW.
                sumSin = sumSin - static_cast<GPUPotAccumType>(charge[iu] * sin(phase));
            }

            GPUtype uc1g = asympk_const * exp(-smear2_4 * waveVectors2[ig]);
            GPUtype uc1gc = uc1g / waveVectors2[ig];

            dcfLongRangeAsympK[2*ig] = static_cast<GPUtype>(uc1gc * sumCos);
            dcfLongRangeAsympK[2*ig+1] = static_cast<GPUtype>(uc1gc * sumSin);

            if(ionic == true){
                GPUtype uc1gh = uc1g / ((waveVectors2[ig] + xappa2) * solvent_dielconst);
                tcfLongRangeAsympK[2*ig] = static_cast<GPUtype>(uc1gh * sumCos);
                tcfLongRangeAsympK[2*ig+1] = static_cast<GPUtype>(uc1gh * sumSin);
            }


        }
    }

    // Serial version for debugging: I will keep it here for now
    void rism3d_potential :: computeSumCosSin_serial(GPUtype* waveVectorX,
                                          GPUtype* waveVectorY,
                                          GPUtype* waveVectorZ,
                                          GPUtype* waveVectors2,
                                          GPUtype cut2_chlk,
                                          GPUtype* position,
                                          GPUtype* charge,
                                          GPUtype* dcfLongRangeAsympK,
                                          GPUtype asympk_const,
                                          GPUtype smear2_4,
                                          int numAtoms,
                                          int numWaveVectors_2,
                                          int start_ind){

        int index = 0;
        for(int i = start_ind; i < numWaveVectors_2; i++){
            if(waveVectors2[i] < cut2_chlk){
                GPUtype sumCos = 0;
                GPUtype sumSin = 0;
                for(int j = 0; j < numAtoms; j++){
                    GPUtype phase = position[j]*waveVectorX[i] + position[j + numAtoms]*waveVectorY[i] + position[j + 2*numAtoms]*waveVectorZ[i];
                    sumCos = sumCos + charge[j] * cos(phase);
                    sumSin = sumSin + charge[j] * sin(phase);
                }
                if(index < 4){
                    cout << sumCos << endl;
                    cout << sumSin << endl;
                    index = index + 1;
                }
                GPUtype uc1g = asympk_const * exp(-smear2_4 * waveVectors2[i]);
                GPUtype uc1gc = uc1g / waveVectors2[i];
                dcfLongRangeAsympK[2*i] = uc1gc * sumCos;
                dcfLongRangeAsympK[2*i+1] = uc1gc * sumSin;
            }
        }

    }

    void rism3d_potential :: calc_dcf_tcf_LongRangeAsympK(GPUtype* waveVectorX,
                                          GPUtype* waveVectorY,
                                          GPUtype* waveVectorZ,
                                          GPUtype* waveVectors2,
                                          GPUtype cut2_chlk,
                                          GPUtype* position,
                                          GPUtype* charge,
                                          GPUtype* dcfLongRangeAsympK,
                                          GPUtype asympk_const,
                                          GPUtype smear2_4,
                                          int numAtoms,
                                          int numWaveVectors_2,
                                          int start_ind, 
                                          bool ionic, GPUtype* tcfLongRangeAsympK,
                                          GPUtype xappa2, GPUtype solvent_dielconst){

        int num_thread = 256;
        int num_blocks = (numWaveVectors_2 + num_thread - 1) / num_thread;

        // Launch the CUDA kernel
        K_calc_dcf_tcf_LongRangeAsympK<<<num_blocks, num_thread>>>(waveVectorX, waveVectorY, waveVectorZ, 
                                                     waveVectors2, cut2_chlk,
                                                     position, charge, 
                                                     dcfLongRangeAsympK,
                                                     asympk_const, smear2_4, 
                                                     numAtoms, numWaveVectors_2,
                                                     start_ind, ionic,
                                                     tcfLongRangeAsympK, xappa2,
                                                     solvent_dielconst);

        // Synchronize to ensure kernel completion
        cudaDeviceSynchronize();

        // Serial version for debugging: I will keep it here for now
        // computeSumCosSin_serial(waveVectorX, waveVectorY, waveVectorZ, 
        //                         waveVectors2, cut2_chlk,
        //                         position, charge, 
        //                         dcfLongRangeAsympK, 
        //                         asympk_const, smear2_4,
        //                         numAtoms, numWaveVectors_2,
        //                         start_ind);

    }

    __global__ void K_calc_sum_cos_sin_huvk0_partial(GPUtype waveVectorX,
                                                    GPUtype waveVectorY,
                                                    GPUtype waveVectorZ,
                                                    GPUtype *positions,    
                                                    GPUtype *charges,      
                                                    int numAtoms,
                                                    GPUtype *blockCos,   // Partial sums (cosines)
                                                    GPUtype *blockSin) { // Partial sums (sines)
        extern __shared__ GPUPotAccumType sharedMemory[]; // Shared memory for partial sums
        GPUPotAccumType *sharedCos = sharedMemory;        // First half for cosines
        GPUPotAccumType *sharedSin = sharedMemory + blockDim.x; // Second half for sines

        int tid = threadIdx.x;
        int idx = blockIdx.x * blockDim.x + tid;

        // Initialize shared memory
        sharedCos[tid] = 0.0;
        sharedSin[tid] = 0.0;

        // Each thread processes a subset of atoms
        for (int i = idx; i < numAtoms; i += blockDim.x * gridDim.x) {
            // GPUtype x = positions[3 * i];
            // GPUtype y = positions[3 * i + 1];
            // GPUtype z = positions[3 * i + 2];
            GPUtype x = positions[i];
            GPUtype y = positions[numAtoms + i];
            GPUtype z = positions[2 * numAtoms + i];
            GPUtype charge = charges[i];

            // Calculate phase
            GPUtype phase = waveVectorX * x + waveVectorY * y + waveVectorZ * z;

            // Accumulate cos and sin contributions
            sharedCos[tid] += static_cast<GPUPotAccumType>(charge * cos(phase));
            sharedSin[tid] -= static_cast<GPUPotAccumType>(charge * sin(phase));
        }

        __syncthreads();

        // Reduce partial sums within the block
        for (int s = blockDim.x / 2; s > 0; s /= 2) {
            if (tid < s) {
                sharedCos[tid] += sharedCos[tid + s];
                sharedSin[tid] += sharedSin[tid + s];
            }
            __syncthreads();
        }

        // Store the block's result in global memory
        if (tid == 0) {
            blockCos[blockIdx.x] = static_cast<GPUtype>(sharedCos[0]);
            blockSin[blockIdx.x] = static_cast<GPUtype>(sharedSin[0]);
        }
    }

    void rism3d_potential :: calc_sum_cos_sin_huvk0_serial(GPUtype waveVectorX_0,
                            GPUtype waveVectorY_0,
                            GPUtype waveVectorZ_0,
                            GPUtype* position,
                            GPUtype* charge,
                            int numAtoms,
                            GPUtype* sumcos_0,
                            GPUtype* sumsin_0){
        *sumcos_0 = 0.0;
        *sumsin_0 = 0.0;
        // int count = 0;
        for(int j = 0; j < numAtoms; j++){
            GPUtype phase = position[j]*waveVectorX_0 + position[j + numAtoms]*waveVectorY_0 + position[j + 2*numAtoms]*waveVectorZ_0;
            *sumcos_0 += charge[j] * cos(phase);
            *sumsin_0 += charge[j] * sin(phase);

            // // Printing a few values to check
            // if(count < 5){
            //     count += 1;
            //     cout << "==================================" << endl;
            //     cout << "index = " << j << endl;
            //     cout << "(x,y,z) = " << position[j] << " " << position[j + numAtoms] << " " << position[j + 2*numAtoms] << endl;
            //     cout << "(kx,ky,kz) = " << waveVectorX_0 << " " << waveVectorY_0 << " " << waveVectorZ_0 << endl;
            //     cout << "charge = " << charge[j] << endl;
            //     cout << "phase = " << phase << endl;
            //     cout << "charge[j] * cos(phase) = " << charge[j] * cos(phase) << endl;
            //     cout << "charge[j] * sin(phase) = " << charge[j] * sin(phase) << endl;
            //     cout << "==================================" << endl;
            // }
        }
    }

    void rism3d_potential :: calc_sum_cos_sin_huvk0(GPUtype waveVectorX_0,
                            GPUtype waveVectorY_0,
                            GPUtype waveVectorZ_0,
                            GPUtype* position,
                            GPUtype* charge,
                            int numAtoms,
                            GPUtype* sumcos_0,
                            GPUtype* sumsin_0){

        // Calling serial version to compare
        // calc_sum_cos_sin_huvk0_serial(waveVectorX_0, waveVectorY_0, waveVectorZ_0,
        //                               position, charge, numAtoms,
        //                               sumcos_0, sumsin_0);

        int numThreads = 256;
        int numBlocks = (numAtoms + numThreads - 1) / numThreads;

        int sharedMemSize = 2 * numThreads * sizeof(GPUPotAccumType);

        GPUtype *blockCos, *blockSin;
        cudaMallocManaged((void**)&blockCos, numBlocks * sizeof(GPUtype));
        cudaMallocManaged((void**)&blockSin, numBlocks * sizeof(GPUtype));

        *sumcos_0 = 0.0;
        *sumsin_0 = 0.0;
        K_calc_sum_cos_sin_huvk0_partial<<<numBlocks, numThreads, sharedMemSize>>>(waveVectorX_0, waveVectorY_0, waveVectorZ_0, 
                                                                                   position, charge, numAtoms, blockCos, blockSin);

        cudaDeviceSynchronize();

        // Accumulate partial sums on the host
        for (int i = 0; i < numBlocks; ++i) {
            *sumcos_0 += blockCos[i];
            *sumsin_0 += blockSin[i];
        }

        // Free unified memory
        cudaFree(blockCos);
        cudaFree(blockSin);
    }

}