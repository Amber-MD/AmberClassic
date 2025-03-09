#ifndef SAFEMEM_HPP
#define SAFEMEM_HPP

#include <iostream>
#include <bits/stdc++.h>
#include <cstring>
using namespace std;

namespace rism3d_c {

    class rism3d_safemem{
        private:
            
        public:

        rism3d_safemem();
        ~rism3d_safemem();
    
        int int_count = 0;
        int real_count = 0;
        int logical_count = 0;
        int string_count = 0;
        int total_count = 0;

        int maxint_count = 0;
        int maxreal_count = 0;
        int maxlogical_count = 0;
        int maxstring_count = 0;
        int max_count = 0;

        string *allocString(int N, string *stringlist);
        int* allocInt(int* data, int dim1);
        double* allocReal(double* data, int dim1);
        double* allocReal(double* data, int dim1, int dim2);
        double* allocReal(double* data, int dim1, int dim2, int dim3);
        double* allocReal(double* data, int dim1, int dim2, int dim3, int dim4);
        double* allocReal(double* data, int dim1, int dim2, int dim3, int dim4, int dim5);

        int memSize_i(int n);
        int memSize_c(string* str, int n_strs);
        int memSize_r(int n_strs);

        void memadd_i(int nbytes);
        void memadd_r(int nbytes);
        void memadd_s(int nbytes);


        // Add for array_class
        void memremov_i(int nbytes);
        void memremov_r(int nbytes);
    };

}

#endif // SAFEMEM_HPP