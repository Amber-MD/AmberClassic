#ifndef RISMUTIL_HPP
#define RISMUTIL_HPP

#include <iostream>
#include <functional>
#include <cuda.h> 
#include "array_class.hpp"
#include "varTypes.hpp"
#include <functional>
#include <cmath>
using namespace std;

namespace rism3d_c{
    
    GPUtype root_newton(std::function<GPUtype(GPUtype)> func,
                        std::function<GPUtype(GPUtype)> deriv, 
                        GPUtype guess, GPUtype tol);

    // void findCenterOfMass(array_class<double> position, double* centerOfMass, array_class<double> mass, int numAtoms);
    void findCenterOfMass(array_class<GPUtype> *position, array_class<GPUtype> *centerOfMass, array_class<GPUtype> *mass, int numAtoms);

    void translate(array_class<GPUtype> *position, int numAtoms, GPUtype trans[3]);

    double get_maxval(array_class<GPUtype> *position, int numAtoms, int idx);
    GPUtype get_maxval_abs(GPUtype *array, int size);
    double get_minval(array_class<GPUtype> *position, int numAtoms, int idx);

    bool isFactorable(int number, int* factor, int size);

    int merge_int(int T_val, int F_val, bool mask);

    int count_unique(float* arr, int size);
    void sort_unique(array_class<GPUtype> *arr, int size, array_class<GPUtype> *sorted,
                    array_class<int> *wv_wv2_map, array_class<int> *wv_wn_map);

    void cu_polinomialInterpolation(int maxPointsToInterp, int waveNumberArraySize, 
                                    int numRDFpoints, int numAtomTypes, 
                                    GPUtype *solventWaveNumbers, 
                                    GPUtype *xvv, GPUtype *gridWaveNumbers,
                                    GPUtype *xvva, int idy, int idz, int Nyz);

    void cu_memcpy(const float *src, float *dst, int size);
    void cu_memcpy(const double *src, double *dst, int size);
    void cu_memcpy(const double *src, float *dst, int size);
    void cu_memcpy(const float *src, double *dst, int size);
                                                            
    void indexArray(const GPUtype* val, int* ptr, int n);

    void siftDown(const GPUtype* val, int* ptr, int root, int n);

    GPUtype cu_reduce_sum(GPUtype* arr, int N);

    void gaussquad_legendre(GPUtype a, GPUtype b, GPUtype* x, GPUtype* weights, int n);
    void legendre(GPUtype x, GPUtype& y, GPUtype& dy, int n);
    void cu_cpyFloatToDouble(const float* src, double* dst, int n);
    void cu_cpyDoubleToFloat(const double* src, float* dst, int n);

    void cu_SDDot(int n, const float* d_x, const float* d_y, float* d_result);
    void cu_SDDot(int n, const float* d_x, const float* d_y, double* d_result);
    // void cu_SDDot(int n, const double* d_x, const double* d_y, double* d_result);


}

#endif // RISMUTIL_HPP