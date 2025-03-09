#include <iostream>
#include "rism3d_safemem.hpp"
using namespace std;

namespace rism3d_c {

    rism3d_safemem :: rism3d_safemem(){
        // cout << "Object from rism3d_safemem is being created!" << endl;
        // defining 1d_excessChemPot, so it can be retrieved from excessChemicalPotential(bool o_lr)
    }

    rism3d_safemem :: ~rism3d_safemem(){
        // cout << "Object from rism3d_safemem is being deleted!" << endl;
    }

    // It is possible to template the functions below!
    // I think we would only need to have some conditionals
    // for memadd_x to ensure we are adding the amount of
    // memory in the correct variable depending on its type.
    // However, here we also have overloaded functions... could it work?
    // Maybe using predefined variables!

    int* rism3d_safemem :: allocInt(int* data, int dim1){
        int nbytes;

        data = new int[dim1];

        nbytes = memSize_i(dim1);
        memadd_i(nbytes);

        return data;
    }

    string* rism3d_safemem :: allocString(int N, string* stringlist){
        int nbytes;

        // cout << "Allocate and set stringlist: N = " << N << endl;
        stringlist = new string[N];

        nbytes = memSize_c(stringlist, N);
        memadd_s(nbytes);

        return stringlist;
    }

    double* rism3d_safemem :: allocReal(double* data, int dim1){
        int nbytes;

        data = new double[dim1];

        nbytes = memSize_r(dim1);
        memadd_r(nbytes);

        return data;
    }

    double* rism3d_safemem :: allocReal(double* data, int dim1, int dim2){
        int nbytes;

        data = new double[dim1*dim2];

        nbytes = memSize_r(dim1*dim2);
        memadd_r(nbytes);

        return data;
    }

    double* rism3d_safemem :: allocReal(double* data, int dim1, int dim2, int dim3){
        int nbytes;

        data = new double[dim1*dim2*dim3];

        nbytes = memSize_r(dim1*dim2*dim3);
        memadd_r(nbytes);

        return data;
    }

    double* rism3d_safemem :: allocReal(double* data, int dim1, int dim2, int dim3, int dim4){
        int nbytes;

        cout << "Print dims" << endl;
        cout << dim1 << " "<< dim2 << " " << dim3 << " " <<dim4 << endl;
        data = new double[dim1*dim2*dim3*dim4];

        nbytes = memSize_r(dim1*dim2*dim3*dim4);
        memadd_r(nbytes);

        return data;
    }

    double* rism3d_safemem :: allocReal(double* data, int dim1, int dim2, int dim3, int dim4, int dim5){
        int nbytes;

        cout << "Print dims" << endl;
        cout << dim1 << " "<< dim2 << " " << dim3 << " " << dim4 << " " << dim5 << endl;
        data = new double[dim1*dim2*dim3*dim4*dim5];

        nbytes = memSize_r(dim1*dim2*dim3*dim4*dim5);
        memadd_r(nbytes);

        return data;
    }

    int rism3d_safemem :: memSize_i(int n){
        int nbytes;

        nbytes = n*sizeof(int);
        return nbytes;
    }

    int rism3d_safemem :: memSize_r(int n){
        int nbytes;

        // nbytes = sizeof(arr) + n*sizeof(double);
        nbytes = n*sizeof(double);
        return nbytes;
    }

    int rism3d_safemem :: memSize_c(string* str, int n_strs){
        int nbytes;

        // nbytes = sizeof(str);
        nbytes = 0;
        for(int i = 0; i < n_strs; i++){
            // nbytes = nbytes + sizeof(str[i]) + (str[i].capacity()+1)*sizeof(char);
            nbytes = (str[i].capacity()+1)*sizeof(char);
        }

        return nbytes;
    }


    void rism3d_safemem :: memadd_i(int nbytes){
        int_count = int_count + nbytes;
        total_count = total_count + nbytes;
        
        if(int_count > maxint_count){
            maxint_count = int_count;
        }
        if(total_count > max_count){
            max_count = total_count;
        }

    }

    void rism3d_safemem :: memadd_r(int nbytes){
        real_count = real_count + nbytes;
        total_count = total_count + nbytes;
        
        if(real_count > maxreal_count){
            maxreal_count = real_count;
        }
        if(total_count > max_count){
            max_count = total_count;
        }

    }

    void rism3d_safemem :: memadd_s(int nbytes){
        string_count = string_count + nbytes;
        total_count = total_count + nbytes;
        
        if(string_count > maxstring_count){
            maxstring_count = string_count;
        }
        if(total_count > max_count){
            max_count = total_count;
        }
    }

    // Add for array_class
    void rism3d_safemem :: memremov_r(int nbytes){
        real_count = real_count - nbytes;
        total_count = total_count - nbytes;
    }

    void rism3d_safemem :: memremov_i(int nbytes){
        int_count = int_count - nbytes;
        total_count = total_count - nbytes;
    }

}