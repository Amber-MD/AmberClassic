#include "rism_util.hpp"

namespace rism3d_c{
    
    GPUtype root_newton(std::function<GPUtype(GPUtype)> func,
                        std::function<GPUtype(GPUtype)> deriv, 
                        GPUtype guess, GPUtype tol){
        
        GPUtype root0 = guess;
        GPUtype root = root0 - func(root0)/deriv(root0);

        for(int i = 0; i < 1000; i++){
            if(abs((root-root0)/root) < tol){
                break;
            }
            root0 = root;
            root = root0 - func(root0)/deriv(root0);    
        }
        return root;
    }

    void findCenterOfMass(array_class<GPUtype> *position, array_class<GPUtype> *centerOfMass, array_class<GPUtype> *mass, int numAtoms){
        GPUtype total_mass = 0;
        for(int j = 0; j < numAtoms; j++){
            total_mass = total_mass + mass->m_data[j];
        }

        for(int i = 0; i < 3; i++){
            centerOfMass->m_data[i] = 0;
            for(int j = 0; j < numAtoms; j++){
                // If position was row major
                centerOfMass->m_data[i] = centerOfMass->m_data[i] + mass->m_data[j]*position->m_data[i*numAtoms + j];

                // If position is column major, we need to use idexing accordingly
                // centerOfMass->m_data[i] = centerOfMass->m_data[i] + mass->m_data[j]*position->m_data[i + j*3];
            }
            centerOfMass->m_data[i] = centerOfMass->m_data[i]/total_mass;
        }

    }

    void translate(array_class<GPUtype> *position, int numAtoms, GPUtype trans[3]){
        for(int i = 0; i < numAtoms; i++){
            // If position on C++ side is row major, we need to use
            position->m_data[0*numAtoms + i] = position->m_data[0*numAtoms + i] + trans[0];
            position->m_data[1*numAtoms + i] = position->m_data[1*numAtoms + i] + trans[1];
            position->m_data[2*numAtoms + i] = position->m_data[2*numAtoms + i] + trans[2];

            // Position comes from Fortran as column major, so if we use the same
            // ordering, the indexing must be done accordingly
            // position->m_data[0 + i*3] = position->m_data[0 + i*3] + trans[0];
            // position->m_data[1 + i*3] = position->m_data[1 + i*3] + trans[1];
            // position->m_data[2 + i*3] = position->m_data[2 + i*3] + trans[2];

        }
    }

    double get_maxval(array_class<GPUtype> *position, int numAtoms, int idx){
        double max_val = 0;

        for(int i = 0; i < numAtoms; i++){
            // If position was row major:
            if(max_val < position->m_data[idx*numAtoms + i]){
                max_val = position->m_data[idx*numAtoms + i];
            }

            // Because position is column major, we need to use the indexing accordingly
            // if(max_val < position->m_data[idx + i*3]){
            //     max_val = position->m_data[idx + i*3];
            // }
        }

        // cout << "max_val = " << max_val << endl;
        return max_val;
    }

    GPUtype get_maxval_abs(GPUtype *array, int size){
        GPUtype max_val = 0;

        for(int i = 0; i < size; i++){
            if(max_val < abs(array[i])){
                max_val = abs(array[i]);
            }
        }

        // cout << "max_val = " << max_val << endl;
        return max_val;
    }


    double get_minval(array_class<GPUtype> *position, int numAtoms, int idx){
        double min_val = 10000;

        for(int i = 0; i < numAtoms; i++){
            // If position was row major:
            if(min_val > position->m_data[idx*numAtoms + i]){
                min_val = position->m_data[idx*numAtoms + i];
            }

            // Because position is column major, we need to use the indexing accordingly
            // if(min_val > position->m_data[idx + i*3]){
            //     min_val = position->m_data[idx + i*3];
            // }
        }

        // cout << "min_val = " << min_val << endl;
        return min_val;
    }

    bool isFactorable(int number, int* factor, int size){
        int numold;
        int num;

        num = number;

        for(int i = 0; i < size; i++){
            numold = num + 1;
            while(numold > num){
                numold = num;
                if(num%factor[i] == 0){
                    num = num / factor[i];
                }
            }    
        }
        if(num == 1){
            return true;
        }
        else{
            return false;
        }    
    }

    int merge_int(int T_val, int F_val, bool mask){
        if(mask == true){
            return T_val;
        }
        else if(mask == false){
            return F_val;
        }
        else{
            cout << "Something wrong with merge_int function! Check boolean argument." << endl;
            abort();
        }
    }

    int count_unique(float* arr, int size){
        sort(arr, arr + size); 

        // Traverse the sorted array and get number of unique values
        int count = 0;
        for (int i = 0; i < size; i++) { 
            if (i == 0 || arr[i] != arr[i - 1]) 
                count = count + 1;
        }

        return count;
    }

    void sort_unique(array_class<GPUtype> *arr, int size, array_class<GPUtype> *sorted,
                    array_class<int> *wv_wv2_map, array_class<int> *wv_wn_map){
        // define temporary array
        GPUtype* temp = new GPUtype[size];

        // copy arr (wavenumber in this case) data into temp
        memcpy(temp, arr->m_data, sizeof(GPUtype)*size);

        // sort in ascending order
        sort(temp, temp + size); 

        // Traverse the sorted array and get number of unique values
        int count = 0;
        int igs;
        for (int i = 0; i < size; i++) { 
            // get index such that arr is in ascending order
            igs = wv_wv2_map->m_data[i];

            // Check if next value in temp is larger than the previous one and 
            // add count if true
            if (i == 0 || temp[i] > temp[i - 1] + std::numeric_limits<GPUtype>::epsilon()) {
                count = count + 1;
            }

            // create an array of indexes that maps the wavevectors2 array 
            // to wavenumbers array (which is stored in ascending order):
            //
            // wavenumbers(wv_wn_map(i)) = wavevectors2(i)
            wv_wn_map->m_data[igs] = count-1;
        }
        
        sorted->alloc_mem(count);
        
        int j = 0;
        for (int i = 0; i < size; i++) { 
            if (i == 0 || temp[i] > temp[i - 1] + std::numeric_limits<GPUtype>::epsilon()){ 
                sorted->m_data[j] = sqrt(temp[i]);
                j = j + 1;
            }
        }

        delete temp;
    }

    // Function to sort indices based on values using heapsort
    // Fills PTR with indices of VAL to give the values of VAL in accending order.
    // I.e.  VAL(PTR) will be in accending order.  The reproduces the functionality
    // of INDEXX from Numerical Recipes.
    // IN:
    //    val :: array of values
    //    ptr :: will be an array of VAL indices.  I.e. it is a 'pointer' to the VAL
    //           elements
    //    n   :: length of the arrays
    void indexArray(const GPUtype* val, int* ptr, int n) {
        // Initialize ptr with indices
        for (int i = 0; i < n; ++i) {
            ptr[i] = i;
        }

        // Build the heap
        for (int i = (n / 2) - 1; i >= 0; --i) {
            siftDown(val, ptr, i, n);
        }

        // Sort the heap
        for (int i = n - 1; i >= 0; --i) {
            // Move the largest element to the end
            std::swap(ptr[0], ptr[i]);

            // Restore the heap property
            siftDown(val, ptr, 0, i);
        }
    }

    // Function to maintain the heap property:
    // Used in the indexArray heapsort.  Brings the largest value between I and M
    // to index I while ensuring that values at 2*I and 2*I+1 are less than the value
    // at I.
    // IN:
    //    val :: array of values
    //    ptr :: will be an array of VAL indices.  I.e. it is a 'pointer' to the VAL
    //           elements
    //    i   :: the start element to use for the arrays
    //    m   :: the end element to use for the arrays
    //    n   :: length of the arrays
    void siftDown(const GPUtype* val, int* ptr, int root, int n) {
        int child = 2 * root + 1; // Left child
        int temp;

        while (child < n) {
            // Check if there is a right child and it is greater than the left child
            if (child + 1 < n && val[ptr[child + 1]] > val[ptr[child]]) {
                child++;
            }

            // If the child is greater than the root, swap them
            if (val[ptr[child]] > val[ptr[root]]) {
                temp = ptr[child];
                ptr[child] = ptr[root];
                ptr[root] = temp;

                // Move down the tree
                root = child;
                child = 2 * root + 1;
            } else {
                break;
            }
        }
    }

    void gaussquad_legendre(GPUtype a, GPUtype b, GPUtype* x, GPUtype* weights, int n){
        int nroot = (n + 1) / 2; // Roots are symmetric; only half need computation
        const GPUtype PI = std::acos(-1.0);

#if RISMCUDA_DOUBLE
        const GPUtype eps = 1e-14;
#else
        const GPUtype eps = 1e-6;
#endif // RISMCUDA_DOUBLE

        for (int iroot = 0; iroot < nroot; ++iroot) {
            // Initial guess for the root
            GPUtype root0 = std::cos(PI * (iroot + 0.75) / (n + 0.5));
            GPUtype root = 0.0, p = 0.0, dp = 0.0;

            // Newton-Raphson iteration
            while (true) {
                legendre(root0, p, dp, n);
                root = root0 - p / dp;
                if (std::abs(root - root0) < eps) break;
                root0 = root;
            }

            // Scale roots and compute weights
            x[iroot] = a + (1.0 - root) * (b - a) / 2.0;
            x[n - iroot - 1] = b - (1.0 - root) * (b - a) / 2.0;

            weights[iroot] = (b - a) / ((1.0 - root * root) * dp * dp);
            weights[n - iroot - 1] = weights[iroot];
        }
    }

    void legendre(GPUtype x, GPUtype& y, GPUtype& dy, int n){
        GPUtype y0 = 1.0, y1 = x;
        y = y0;
        dy = 0.0;

        if (n == 0) return;
        y = y1;
        dy = 1.0;

        if (n == 1) return;

        for (int i = 2; i <= n; ++i) {
            GPUtype yn = ((2.0 * i - 1.0) * x * y1 - (i - 1.0) * y0) / i;
            dy = (x * yn - y1) * i / (x * x - 1.0);
            y0 = y1;
            y1 = yn;
            y = yn;
        }
    }

}