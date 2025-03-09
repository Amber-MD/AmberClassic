#include <iostream>
#include "mdiis.hpp"
using namespace std;

namespace rism3d_c {

    mdiis :: mdiis(int nvec, double del, rism3d_safemem *memalloc) :
    a(array_class<GPUtype>::ROW_MAJOR),
    b()
    {
        memalloc_p = memalloc;
        // cout << "MDIIS object created" << endl;
        NVec = nvec;
        mdiis_del = del;

        current_NVec = 1;

        cublasCreate(&handle);
        cusolverDnCreate(&handle_solv);

        a.set_memalloc(memalloc_p, true);
        a.alloc_mem(NVec+1,NVec+1);

        b.set_memalloc(memalloc_p, true);
        b.alloc_mem(NVec+1);

        cudaMemset(a.m_data, 0, (NVec + 1) * (NVec + 1) * sizeof(GPUtype));
        cudaMemset(b.m_data, 0, (NVec + 1) * sizeof(GPUtype));

        // The following lines will initialize a and b data on the cpu side, so
        // it will be needed to intilized on the device
        // or do a cudaMemPrefetch to ensure further optimization

        for(int i = 1; i < NVec+1; i++){
            a(0,i) = -1;
            a(i,0) = -1;
        }

        b(0) = -1;

        
    }

    mdiis :: ~mdiis(){
        cublasDestroy_v2(handle);
        cusolverDnDestroy(handle_solv);
        // cout << "MDIIS object destroyed" << endl;
    }

    void mdiis :: resize(GPUtype *vecData, GPUtype *resVecData, int numAtomTypes, int numKx, int numKy, int numKz){
        // Number of total points in the arrays
        np = numAtomTypes * numKx * numKy * numKz;

        // Number of total points in the real space: necessary for calculating the residue
        np_real = numAtomTypes * numKx * numKy * (numKz-2);

        // Array containing the solutions
        xi = vecData;

        // Array containing the residuals
        ri = resVecData;
    }

    void mdiis :: advance(){
        // Allocating memory for R*
        GPUtype *ri_star;
        cudaMallocManaged(&ri_star, np * sizeof(GPUtype));
        cudaMemset(ri_star, 0, np * sizeof(GPUtype));

        // Create temporary workspace for b
        GPUtype *temp_b;
        GPUtype *temp_a;
        GPUtype *new_sol;
        cudaMallocManaged(&temp_b, (NVec + 1) * sizeof(GPUtype));
        cudaMallocManaged(&temp_a, (NVec + 1) * (NVec + 1) * sizeof(GPUtype));
        cudaMallocManaged(&new_sol, np * sizeof(GPUtype));

        // Copy original b to temp_b
        cudaMemcpy(temp_b, b.m_data, (NVec + 1) * sizeof(GPUtype), cudaMemcpyDeviceToDevice);
        
        // Calculating overlaping matrix values 
        for(int i = 0; i < current_NVec; i++){
                for(int j = 0; j < current_NVec; j++){
#if RISMCUDA_DOUBLE
                    cublasDdot(handle, np, &ri[i*np], 1, &ri[j*np], 1, &a.m_data[(i+1) * (NVec + 1) + (j+1)]);
#else
                    cublasSdot(handle, np, &ri[i*np], 1, &ri[j*np], 1, &a.m_data[(i+1) * (NVec + 1) + (j+1)]);
#endif // RISMCUDA_DOUBLE
                }
        }

        // get mean square value of new residual, considering the number of points in
        // the real space
        residual = sqrt(a(current_NVec, current_NVec)/np_real);

        // Synchronize device before proceeding
        cudaDeviceSynchronize();

        // Copy a data into temp_a
        cudaMemcpy(temp_a, a.m_data, (NVec + 1) * (NVec + 1) * sizeof(GPUtype), cudaMemcpyDeviceToDevice);

        // Prepare the necessary variables for LU factorization
        int *d_ipiv; // Pivot indices
        int *d_info; // Info for the factorization status
        int lwork = 0; // Workspace size for the factorization

        // Allocate memory for pivot array and info
        cudaMallocManaged(&d_ipiv, (NVec + 1) * sizeof(int));
        cudaMallocManaged(&d_info, sizeof(int));

        // Solving Ax = b using LU decomposition and storing solution into a temp array
        GPUtype *d_work;
#if RISMCUDA_DOUBLE
        // Query the workspace size
        cusolverDnDgetrf_bufferSize(handle_solv, (current_NVec + 1), (current_NVec + 1), temp_a, (NVec + 1), &lwork);
        
        cudaMallocManaged(&d_work, lwork * sizeof(GPUtype));
        
        // Perform LU factorization
        cusolverDnDgetrf(handle_solv, (current_NVec + 1), (current_NVec + 1), temp_a, (NVec + 1), d_work, d_ipiv, d_info);
        cudaDeviceSynchronize();

        // Solve the system Ax = b
        cusolverDnDgetrs(handle_solv, CUBLAS_OP_T, (current_NVec + 1), 1, temp_a, (NVec + 1), d_ipiv, temp_b, (NVec + 1), d_info);
        cudaDeviceSynchronize();
#else
        // Query the workspace size
        cusolverDnSgetrf_bufferSize(handle_solv, (current_NVec + 1), (current_NVec + 1), temp_a, (NVec + 1), &lwork);
        
        cudaMallocManaged(&d_work, lwork * sizeof(GPUtype));

        // Perform LU factorization
        cusolverDnSgetrf(handle_solv, (current_NVec + 1), (current_NVec + 1), temp_a, (NVec + 1), d_work, d_ipiv, d_info);
        cudaDeviceSynchronize();

        // Solve the system Ax = b
        cusolverDnSgetrs(handle_solv, CUBLAS_OP_T, (current_NVec + 1), 1, temp_a, (NVec + 1), d_ipiv, temp_b, (NVec + 1), d_info);
        cudaDeviceSynchronize();
#endif // RISMCUDA_DOUBLE

        // Just checking results:
        // For NVec = 1:
        // A = | 0  -1  |
        //     | -1 S11 |
        // x = | lamb |
        //     |  c   |
        // b = | -1 |
        //     |  0 |
        // So, c = 1 and lamb = S11.
        // cout << "a(1,1) = " << a(1,1) << endl;
        // cout << "lambda = " << temp_b[0] << endl;
        // cout << "c1 = " << temp_b[1] << endl;
        // cout << "c2 = " << temp_b[2] << endl;
        // cout << "c3 = " << temp_b[3] << endl;

        // Caclulate new solution as ri_new = ri_new + mdiis_del * ri_star
        // first calculating ri_star values
        // then calculating the new solution as a linear combination of previous solutions
        // and finally adding mdiis_del * ri_star values to the new solution
        GPUtype del = mdiis_del;

        // After solving Ax = b, b stores the solutions to the problem.
        // cublasDaxpy calculates y = alpha * x + y, 
        // with alpha a constant, x and y being two arrays.
        // Thus, the following two for loops are calculating 
        // f = sum_i c_i * y_i, with y_i/f being ri/ri_star or xi/new_sol
        for(int i = 0; i < current_NVec; i++){
#if RISMCUDA_DOUBLE
            cublasDaxpy(handle, np, &temp_b[i+1], &ri[i*np], 1, ri_star, 1);
            cudaDeviceSynchronize();
#else
            cublasSaxpy(handle, np, &temp_b[i+1], &ri[i*np], 1, ri_star, 1);
            cudaDeviceSynchronize();
#endif // RISMCUDA_DOUBLE
        }

        for(int i = 0; i < current_NVec; i++){
#if RISMCUDA_DOUBLE
            cublasDaxpy(handle, np, &temp_b[i+1], &xi[i*np], 1, new_sol, 1);
            cudaDeviceSynchronize();
#else
            cublasSaxpy(handle, np, &temp_b[i+1], &xi[i*np], 1, new_sol, 1);
            cudaDeviceSynchronize();
#endif // RISMCUDA_DOUBLE
        }


        // The final solution in the MDIIS method is calculated as 
        // f_final = new_sol + del * ri_star
        // The code below is performing this calculation
#if RISMCUDA_DOUBLE
        cublasDaxpy(handle, np, &del, ri_star, 1, new_sol, 1);
        cudaDeviceSynchronize();
#else
        cublasSaxpy(handle, np, &del, ri_star, 1, new_sol, 1);
        cudaDeviceSynchronize();
#endif // RISMCUDA_DOUBLE

        // Move data starting from the second array to the beginning to free space for new solutions
        if(current_NVec == NVec){
            cudaMemcpy(xi, xi + np, (NVec - 1)*np*sizeof(GPUtype), cudaMemcpyDeviceToDevice);
            cudaMemcpy(ri, ri + np, (NVec - 1)*np*sizeof(GPUtype), cudaMemcpyDeviceToDevice);
            cudaDeviceSynchronize();
        }

        // Increment the current overlap matrix size
        if(current_NVec < NVec){
            current_NVec = current_NVec + 1;
        }

        // Need to copy solution to the next array to be used in the RISM/Picard[?] iteration
        if(current_NVec > 1){
            cudaMemcpy(xi + (current_NVec - 1)*np, new_sol, np*sizeof(GPUtype), cudaMemcpyDeviceToDevice);
            cudaDeviceSynchronize();
        }

        // free space
        cudaFree(temp_b);
        cudaFree(temp_a);
        cudaFree(d_ipiv);
        cudaFree(d_info);
        cudaFree(d_work);
        cudaFree(ri_star);
        cudaFree(new_sol);
    }

    int mdiis :: getWRKvec(){
        return current_NVec;
    }

}