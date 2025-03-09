#ifndef ARRAYCLASS_HPP
#define ARRAYCLASS_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cuda.h> 
#include <cuda_runtime.h>
#include "rism3d_safemem.hpp"
#include <stdexcept>

using namespace std;
using namespace rism3d_c;

////////////////////////////////////// Declaration //////////////////////////////////////

inline cudaError_t checkCuda(cudaError_t result)
{
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    assert(result == cudaSuccess);
    printf("Aborting...\n");
    abort();
  }
  return result;
}

template <typename T> 
class array_class {
public:
    // Define variables to store data as column or row major
    enum StorageOrder { ROW_MAJOR, COLUMN_MAJOR };

    // Store the pointer to a 1D array which contains all data
    T* m_data = nullptr;
    
    // Initializes the object:
    // IN:
    //    None, a, b, c, d, e: dimensions for the nD array (n from 1 to 5)
    //    obs: If no argument is given in constructor, one must call alloc_mem with desired arguments
    array_class(StorageOrder storage_order = COLUMN_MAJOR);

    void setStorageOrder(StorageOrder new_order);

    array_class(size_t a, rism3d_safemem *memalloc);
    array_class(size_t a, size_t b, rism3d_safemem *memalloc);
    array_class(size_t a, size_t b, size_t c, rism3d_safemem *memalloc);
    array_class(size_t a, size_t b, size_t c, size_t d, rism3d_safemem *memalloc);
    array_class(size_t a, size_t b, size_t c, size_t d, size_t e, rism3d_safemem *memalloc);

    ~array_class();

    void set_memalloc(rism3d_safemem *memalloc, bool uniMem_v = false);

    void alloc_mem(size_t a);
    void alloc_mem(size_t a, size_t b);
    void alloc_mem(size_t a, size_t b, size_t c);
    void alloc_mem(size_t a, size_t b, size_t c, size_t d);
    void alloc_mem(size_t a, size_t b, size_t c, size_t d, size_t e);

    // Counts the number of bytes for a given array of type T
    // IN:
    //    n: the total size of the array
    // OUT:
    //    the total number of bytes stored in memory
    int count_bytes(int n);

    // Returns the size of 1st, 2nd, ..., 5th dimention
    int get_dims(int id);
    size_t get_dim1();
    size_t get_dim2();
    size_t get_dim3();
    size_t get_dim4();
    size_t get_dim5();

    // Defining operators to retrive specific array index 
    // in a multidimentional representation.
    // IN:
    //    array position (e.g. (2,3,1) for a 3D array)
    // OUT:
    //    the value at the requested position
    inline T& operator()(size_t a);
    inline T& operator()(size_t a, size_t b);
    inline T& operator()(size_t a, size_t b, size_t c);
    inline T& operator()(size_t a, size_t b, size_t c, size_t d);
    inline T& operator()(size_t a, size_t b, size_t c, size_t d, size_t e);

private:
    // Define type of memory allocation
    bool uniMem;

    // Store the total number of bytes used
    int nbytes;

    // Store to array dimention sizes (default set to 0)
    size_t m_a = 0, m_b = 0, m_c = 0, m_d = 0, m_e = 0;

    // Storage order variable
    StorageOrder order;

    // Pointer to rism3d_safemem type object 
    // (necessary for keep track of memory allocation and deallocation)
    rism3d_safemem *memalloc_ptr;

    void call_exception();
};

////////////////////////////////////// Definition //////////////////////////////////////

template<typename T> 
array_class<T> :: array_class(StorageOrder storage_order) : order(storage_order){ }

template<typename T>
void array_class<T> :: setStorageOrder(StorageOrder new_order){
    int x_dim = (m_a == 0) ? (m_a + 1) : m_a;
    int y_dim = (m_b == 0) ? (m_b + 1) : m_b;
    int z_dim = (m_c == 0) ? (m_c + 1) : m_c;
    int u_dim = (m_d == 0) ? (m_d + 1) : m_d;
    int v_dim = (m_e == 0) ? (m_e + 1) : m_e;

    if(new_order != order){
        T *new_data = new T[x_dim * y_dim * z_dim * u_dim * v_dim];

        try{
            for (int i = 0; i < x_dim; ++i) {
                for (int j = 0; j < y_dim; ++j) {
                    for (int k = 0; k < z_dim; k++){
                        for(int p = 0; p < u_dim; p++){
                            for(int o = 0; o < v_dim; o++){
                                if (order == ROW_MAJOR) {
                                    new_data[i + j * x_dim + k * x_dim * y_dim + p * x_dim * y_dim * z_dim + o * x_dim * y_dim * z_dim * u_dim] = m_data[i * y_dim * z_dim * u_dim * v_dim + j * z_dim * u_dim * v_dim + k * u_dim * v_dim + o * v_dim + p];
                                } else {
                                    new_data[i * y_dim * z_dim * u_dim * v_dim + j * z_dim * u_dim * v_dim + k * u_dim * v_dim + o * v_dim + p] = m_data[i + j * x_dim + k * x_dim * y_dim + p * x_dim * y_dim * z_dim + o * x_dim * y_dim * z_dim * u_dim];
                                }
                            }
                        }
                    }
                }
            }
        } catch(...) {
            delete[] new_data;
            throw runtime_error("Failed to copy data during storage order change.");
            abort();
        }
        
        delete[] m_data;
        m_data = new_data;
        order = new_order;
    }
}

template<typename T> 
array_class<T> :: array_class(size_t a, rism3d_safemem *memalloc)
    : m_a(a), m_data(new T[a]){
    nbytes = count_bytes(a);
    if(is_same<T, double>::value){
        memalloc->memadd_r(nbytes);
    }
    if(is_same<T, int>::value){
        memalloc->memadd_i(nbytes);
    }
    memalloc_ptr = memalloc;
}

template<typename T> 
array_class<T> :: array_class(size_t a, size_t b, rism3d_safemem *memalloc)
    : m_a(a), m_b(b), m_data(new T[a * b]){
    nbytes = count_bytes(a*b);
    if(is_same<T, double>::value){
        memalloc->memadd_r(nbytes);
    }
    if(is_same<T, int>::value){
        memalloc->memadd_i(nbytes);
    }
    memalloc_ptr = memalloc;
}

template<typename T> 
array_class<T> :: array_class(size_t a, size_t b, size_t c, rism3d_safemem *memalloc)
    : m_a(a), m_b(b), m_c(c), m_data(new T[a * b * c]){
    nbytes = count_bytes(a*b*c);
    if(is_same<T, double>::value){
        memalloc->memadd_r(nbytes);
    }
    if(is_same<T, int>::value){
        memalloc->memadd_i(nbytes);
    }
    memalloc_ptr = memalloc;
}

template<typename T> 
array_class<T> :: array_class(size_t a, size_t b, size_t c, size_t d, rism3d_safemem *memalloc)
    : m_a(a), m_b(b), m_c(c), m_d(d), m_data(new T[a * b * c * d]){
    nbytes = count_bytes(a*b*c*d);
    if(is_same<T, double>::value){
        memalloc->memadd_r(nbytes);
    }
    if(is_same<T, int>::value){
        memalloc->memadd_i(nbytes);
    }
    memalloc_ptr = memalloc;        
}

template<typename T> 
array_class<T> :: array_class(size_t a, size_t b, size_t c, size_t d, size_t e, rism3d_safemem *memalloc)
    : m_a(a), m_b(b), m_c(c), m_d(d), m_e(e), m_data(new T[a * b * c * d * e]){
    nbytes = count_bytes(a*b*c*d*e);
    if(is_same<T, double>::value){
        memalloc->memadd_r(nbytes);
    }
    if(is_same<T, int>::value){
        memalloc->memadd_i(nbytes);
    }
    memalloc_ptr = memalloc;
}

template<typename T> 
array_class<T> :: ~array_class(){
    if(is_same<T, double>::value){
        memalloc_ptr->memremov_r(nbytes);
    }
    if(is_same<T, int>::value){
        memalloc_ptr->memremov_i(nbytes);
    }
    
    if(uniMem == false){
        delete[] m_data;
    }
    else{
        checkCuda(cudaFree(m_data));
    }
}

template<typename T>
void array_class<T> :: set_memalloc(rism3d_safemem *memalloc, bool uniMem_v){
    memalloc_ptr = memalloc;
    uniMem = uniMem_v;
}

template<typename T>
void array_class<T> :: alloc_mem(size_t a){
    m_a = a;
    if(uniMem == false){
        m_data = new T[a];
    }
    else{
        checkCuda(cudaMallocManaged((void**)&m_data, a*sizeof(T)));
    }

    nbytes = count_bytes(a);
    if(is_same<T, double>::value){
        memalloc_ptr->memadd_r(nbytes);
    }
    if(is_same<T, int>::value){
        memalloc_ptr->memadd_i(nbytes);
    }
}

template<typename T>
void array_class<T> :: alloc_mem(size_t a, size_t b){
    m_a = a;
    m_b = b;
    if(uniMem == false){
        m_data = new T[a*b];
    }
    else{
        checkCuda(cudaMallocManaged((void**)&m_data, a*b*sizeof(T)));
    }
    
    nbytes = count_bytes(a*b);
    if(is_same<T, double>::value){
        memalloc_ptr->memadd_r(nbytes);
    }
    if(is_same<T, int>::value){
        memalloc_ptr->memadd_i(nbytes);
    }
}

template<typename T>
void array_class<T> :: alloc_mem(size_t a, size_t b, size_t c){
    m_a = a;
    m_b = b;
    m_c = c;
    if(uniMem == false){
        m_data = new T[a*b*c];
    }
    else{
        checkCuda(cudaMallocManaged((void**)&m_data, a*b*c*sizeof(T)));
    }
    
    nbytes = count_bytes(a*b*c);
    if(is_same<T, double>::value){
        memalloc_ptr->memadd_r(nbytes);
    }
    if(is_same<T, int>::value){
        memalloc_ptr->memadd_i(nbytes);
    }
}

template<typename T>
void array_class<T> :: alloc_mem(size_t a, size_t b, size_t c, size_t d){
    m_a = a;
    m_b = b;
    m_c = c;
    m_d = d;
    if(uniMem == false){
        m_data = new T[a*b*c*d];
    }
    else{
        checkCuda(cudaMallocManaged((void**)&m_data, a*b*c*d*sizeof(T)));
    }

    nbytes = count_bytes(a*b*c*d);
    if(is_same<T, double>::value){
        memalloc_ptr->memadd_r(nbytes);
    }
    if(is_same<T, int>::value){
        memalloc_ptr->memadd_i(nbytes);
    }
}

template<typename T>
void array_class<T> :: alloc_mem(size_t a, size_t b, size_t c, size_t d, size_t e){
    m_a = a;
    m_b = b;
    m_c = c;
    m_d = d;
    m_e = e;
    if(uniMem == false){
        m_data = new T[a*b*c*d*e];
    }
    else{
        checkCuda(cudaMallocManaged((void**)&m_data, a*b*c*d*e*sizeof(T)));
    }

    nbytes = count_bytes(a*b*c*d*e);
    if(is_same<T, double>::value){
        memalloc_ptr->memadd_r(nbytes);
    }
    if(is_same<T, int>::value){
        memalloc_ptr->memadd_i(nbytes);
    }
}

template<typename T>
int array_class<T> :: count_bytes(int n){
    int nbytes = n*sizeof(T);
    return nbytes;
}

template<typename T>
int array_class<T> :: get_dims(int id){
    switch(id){
        case 0:
            return m_a;
        case 1:
            return m_b;
        case 2:
            return m_c;
        case 3:
            return m_d;
        case 4:
            return m_e;
        default:
            throw std::invalid_argument("Invalid value. Arrays with more than 5 dimensions not supported, yet.");
            abort();
    }
}

template<typename T>
size_t array_class<T> :: get_dim1(){
    return m_a;
}

template<typename T>
size_t array_class<T> :: get_dim2(){
    return m_b;
}

template<typename T>
size_t array_class<T> :: get_dim3(){
    return m_c;
}

template<typename T>
size_t array_class<T> :: get_dim4(){
    return m_d;
}

template<typename T>
size_t array_class<T> :: get_dim5(){
    return m_e;
}

template<typename T> 
inline T& array_class<T> :: operator()(size_t a){
    return m_data[a];
}

template<typename T> 
inline T& array_class<T> :: operator()(size_t a, size_t b){
    if(order == ROW_MAJOR){
        return m_data[a * m_b + b ];
    }
    else if(order == COLUMN_MAJOR){
        return m_data[a + b * m_a];
    }
    else{
        throw std::invalid_argument("Invalid order. Should be: ROW_Major or COLUMN_MAJOR");
    }
    
}

template<typename T> 
inline T& array_class<T> :: operator()(size_t a, size_t b, size_t c){
    if(order == ROW_MAJOR){
        return m_data[a * m_b * m_c + b * m_c + c];
    }
    else if(order == COLUMN_MAJOR){
        return m_data[a + b * m_a + c * m_a * m_b];
    }
    else{
        throw std::invalid_argument("Invalid order. Should be: ROW_Major or COLUMN_MAJOR");
        abort();
    }
}

template<typename T> 
inline T& array_class<T> :: operator()(size_t a, size_t b, size_t c, size_t d){
    if(order == ROW_MAJOR){
        return m_data[a * m_b * m_c * m_d + b * m_c * m_d + c * m_d + d ];
    }
    else if(order == COLUMN_MAJOR){
        return m_data[a + b * m_a + c * m_a * m_b + d * m_a * m_b * m_c];
    }
    else{
        throw std::invalid_argument("Invalid order. Should be: ROW_Major or COLUMN_MAJOR");
        abort();
    }
}

template<typename T> 
inline T& array_class<T> :: operator()(size_t a, size_t b, size_t c, size_t d, size_t e){
    if(order == ROW_MAJOR){
        return m_data[a * m_b * m_c * m_d * m_e + b * m_c * m_d * m_e + c * m_d * m_e + d * m_e + e];
    }
    else if(order == COLUMN_MAJOR){
        return m_data[a + b * m_a + c * m_a * m_b + d * m_a * m_b * m_c + e * m_a * m_b * m_c * m_d];
    }
    else{
        throw std::invalid_argument("Invalid order. Should be: ROW_Major or COLUMN_MAJOR");
        abort();
    }
}

#endif //ARRAYCLASS