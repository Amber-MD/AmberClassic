#include <iostream>
#include "rism3d.hpp"
#include "rism_util.hpp"
#include <sstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

namespace rism3d_c {

    void rism3d :: set_atomname(char** names){
        // cout << "setting atom names..." << endl;
        for(int i=0; i < solventclass.numAtomTypes; i++){
            solventclass.atomName[i] = names[i];
            // cout << "====================================================" << endl;
            // cout << solventclass.atomName[i] << endl;
            // cout << "====================================================" << endl;
        }

    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// BEGIN... FUNCTIONS CALLED FROM AMBER_RISM_INTERFACE.F90
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    // rism object constructor:
    rism3d::rism3d(solute_cpp *solu_f, solvent_cpp *solv_f, int centering, int ncuvsteps,
        char** closure, double cut, int mdiis_nvec, double mdiis_del, int mdiis_method,
        double mdiis_restart, bool treeDCF, bool treeTCF, bool treeCoulomb, 
        double treeDCFMAC, double treeTCFMAC, double treeCoulombMAC,
        int treeDCFOrder, int treeTCFOrder, int treeCoulombOrder, 
        int treeDCFN0, int treeTCFN0, int treeCoulombN0,
        double asympKSpaceTolerance, double ljTolerance, double chargeSmear,
        double o_buffer, double* o_grdspc, int o_mpicomm,
        char* o_periodic, double* o_unitCellDimensions, 
        double o_biasPotential) : 
        memalloc(),
        soluteclass(&memalloc, solu_f),
        solventclass(&memalloc, solv_f),
        grid(&memalloc),
        fft(&grid),
        pot(&memalloc, &grid, &soluteclass, &solventclass, cut, o_periodic, treeDCF, treeTCF, treeCoulomb, chargeSmear),
        mdiisclass(mdiis_nvec, mdiis_del, &memalloc),
        dxclass(),
        closure(nullptr),
        centering(centering),
        ncuvsteps(ncuvsteps),
        asympKSpaceTolerance(asympKSpaceTolerance),
        ljTolerance(ljTolerance),
        NVec(mdiis_nvec),
        mdiis_method(mdiis_method),
        deloz(mdiis_del),
        mdiis_restart(mdiis_restart),
        xvva(array_class<GPUtype>::ROW_MAJOR),
        cuv(array_class<GPUtype>::ROW_MAJOR),
        cuvWRK(array_class<GPUtype>::ROW_MAJOR),
        cuvres(array_class<GPUtype>::ROW_MAJOR),
        cuvresWRK(array_class<GPUtype>::ROW_MAJOR),
        guv(array_class<GPUtype>::ROW_MAJOR)
        {
        // Checking restrictions for running this first cuda version
        if(strcmp(o_periodic, "") != 0){
            periodic = true;
            cout << "Periodic boundary conditions not supported, yet." << endl;
            abort();
        }
        if(o_buffer <= 0){
            cout << "Buffer <= 0 is not supported, yet." << endl;
            abort();
        }
        if(ljTolerance != 0){
            cout << "ljTolerance != 0 is not supported, yet." << endl;
            abort();
        }

        nsolution = new int;
        nsolutionChg = new int;
        *nsolutionChg = 0;
        nsolutionNoChg = new int;
        *nsolutionNoChg = 0;
        nsolution = nsolutionChg;

        // create cuBLAS handle
        status = cublasCreate(&handle);
        if (status != CUBLAS_STATUS_SUCCESS) {
            std::cerr << "cuBLAS handle creation failed!" << std::endl;
            abort();
        }

        setbox_variable(&grid, o_buffer, o_grdspc);
        
        setcut(cut);

        grid.setUnitCellDimensions(o_unitCellDimensions, periodic);

    }


    // defining second version of the rism object constructor for the case of fixed box size:
    rism3d::rism3d(solute_cpp *solu_f, solvent_cpp *solv_f, int centering, int ncuvsteps,
        char** closure, double cut, int mdiis_nvec, double mdiis_del, int mdiis_method,
        double mdiis_restart, bool treeDCF, bool treeTCF, bool treeCoulomb, 
        double treeDCFMAC, double treeTCFMAC, double treeCoulombMAC,
        int treeDCFOrder, int treeTCFOrder, int treeCoulombOrder, 
        int treeDCFN0, int treeTCFN0, int treeCoulombN0,
        double asympKSpaceTolerance, double ljTolerance, double chargeSmear,
        double* o_boxlen, int* o_ng3, int o_mpicomm, 
        char* o_periodic, double* o_unitCellDimensions, 
        double o_biasPotential) :
        memalloc(), 
        soluteclass(&memalloc, solu_f),
        grid(&memalloc),
        fft(&grid),
        pot(&memalloc, &grid, &soluteclass, &solventclass, cut, o_periodic, treeDCF, treeTCF, treeCoulomb, chargeSmear),
        mdiisclass(mdiis_nvec, mdiis_del, &memalloc),
        dxclass(),
        // closure(&solventclass, &soluteclass, &grid, &pot, &memalloc),
        closure(nullptr),
        centering(centering),
        ncuvsteps(ncuvsteps),
        asympKSpaceTolerance(asympKSpaceTolerance),
        ljTolerance(ljTolerance),
        NVec(mdiis_nvec),
        mdiis_method(mdiis_method),
        deloz(mdiis_del),
        mdiis_restart(mdiis_restart),
        xvva(array_class<GPUtype>::ROW_MAJOR),
        cuv(array_class<GPUtype>::ROW_MAJOR),
        cuvWRK(array_class<GPUtype>::ROW_MAJOR),
        cuvres(array_class<GPUtype>::ROW_MAJOR),
        cuvresWRK(array_class<GPUtype>::ROW_MAJOR),
        guv(array_class<GPUtype>::ROW_MAJOR)
        {
        
        // This implementation was not clean or tested, yet.
        cout << "Fixed box size not supported/tested at this moment" << endl;
        abort();


        cout << "Object rism3d initialized!" << endl;

        NVec = mdiis_nvec;

        // Receive solv struct from F90, allocate memory and copy data:
        // set_solvent(solv);
        // set_solute(solu);

        // cout << "Setting solvent and solute" << endl;

        // solv = memalloc.set_solvent(solv_f,solv);
        // solu = memalloc.set_solute(solu_f,solu);
        // solv = set_solvent(solv_f,solv);
        // solu = set_solute(solu_f,solu);
        ///////////////////////////////////////////////////////////////

        nsolution = new int;
        nsolutionChg = new int;
        nsolutionNoChg = new int;

        // create cuBLAS handle
        status = cublasCreate(&handle);
        if (status != CUBLAS_STATUS_SUCCESS) {
            std::cerr << "cuBLAS handle creation failed!" << std::endl;
            abort();
        }

        // guv = memalloc.allocReal(guv,grid.totalLocalPointsK,solventclass.numAtomTypes);
        // guv_dT = memalloc.allocReal(guv,grid.totalLocalPointsK,solventclass.numAtomTypes);
        // huv = memalloc.allocReal(guv,grid.totalLocalPointsK,solventclass.numAtomTypes);
        // cuv = memalloc.allocReal(cuv,grid.localDimsR[0],grid.localDimsR[1],grid.localDimsR[2],solventclass.numAtomTypes);
        // cuvWRK = memalloc.allocReal(cuvWRK,grid.localDimsR[0],grid.localDimsR[1],grid.localDimsR[2],solventclass.numAtomTypes,NVec);
        // electronmap = memalloc.allocReal(electronmap,grid.localDimsR[0],grid.localDimsR[1],grid.localDimsR[2]);

        // cuvWRK.set_memalloc(&memalloc);
        // cuvWRK.alloc_mem(grid.localDimsR[0],grid.localDimsR[1],grid.localDimsR[2],solventclass.numAtomTypes);

        *nsolution = 1;
        periodic = false;

        //////////// print some values...
        cout << "centering = " << centering << endl;
        cout << "ncuvsteps = " << this->ncuvsteps << endl;
        cout << "mdiis_restart = " << mdiis_restart << endl;
        cout << "Closure names received to initialization are: " << endl;
        for(int i = 0; i < 3; i++){
            cout << closure[i] << endl;
        }
        
        setbox_fixed(o_boxlen, o_ng3);
        cout << "boxlen and num_grid_pts: " << endl;
        for(int i = 0; i < 3; i++){
            cout << "boxlen[" << i << "] = " << fixedBoxDimensionsR[i] << endl;
            cout << "num_grid_pts[" << i << "] = " << fixedNumGridPoints[i] << endl;
            // cout << "boxlen[" << i << "] = " << o_boxlen[i] << endl;
            // cout << "num_grid_pts[" << i << "] = " << o_ng3[i] << endl;
        }

        // for(int i = 0; i < 3; i++){
        //     cout << "o_grdspc[" << i << "] = " << o_grdspc[i] << endl;
        // }
        cout << "o_periodic = " << o_periodic << endl;
        cout << "solvent numAtomTypes = " << solventclass.numAtomTypes << endl;
        /////////////////////////////////////////////////////////////// 
        
        // fftw_planner = FFT_ESTIMATE
        fft_aligned = true;
        fftw_localtrans = true;

        verbose = 0;
        rism_failure = false;
        // buffer = 12.0;
        varbox = true;
        deloz = 0.7;
        periodicPotential = "testing...";

        ////////////////// CREATING VALUES JUST TO TEST
        // Populating 1D array to match a NxM matrix:
        // i = row indx
        // j = column indx
        // 1D_array[i + N*j]

        // CALLING FUNCTIONS FROM POTENTIAL CLASS HERE TO TEST ON FORTRAN SIDE
        pot.long_range_asymptotics();
        // cout << "Calling Pot method potential_calc" << endl;
        // pot.potential_calc(solventclass.numAtomTypes);

        // int count = 0;
        // for(int i = 0; i < grid.totalLocalPointsK; i++){
        //     for(int j = 0; j < solventclass.numAtomTypes; j++){
        //         guv[i + grid.totalLocalPointsK*j] = count;
        //         count++;
        //     }
        // }
        // count = 20;
        // for(int i = 0; i < grid.totalLocalPointsK; i++){
        //     for(int j = 0; j < solventclass.numAtomTypes; j++){
        //         guv_dT[solventclass.numAtomTypes*i + j] = count;
        //         count--;
        //     }
        // }

        // count = 121;
        // for(int i = 0; i < grid.totalLocalPointsK; i++){
        //     for(int j = 0; j < solventclass.numAtomTypes; j++){
        //         huv[i + grid.totalLocalPointsK*j] = count;
        //         count++;
        //     }
        // }

        // count = 0;
        // for(int i = 0; i < grid.localDimsR[0]; i++){
        //     for(int j = 0; j < grid.localDimsR[1]; j++){
        //         for(int k = 0; k < grid.localDimsR[2]; k++){
        //             electronmap[i + grid.localDimsR[0]*j + grid.localDimsR[0]*grid.localDimsR[1]*k] = count;
        //             count++;
        //         }
        //     }
        // }

        // count = 666;
        // for(int i = 0; i < grid.localDimsR[0]; i++){
        //     for(int j = 0; j < grid.localDimsR[1]; j++){
        //         for(int k = 0; k < grid.localDimsR[2]; k++){
        //             for(int q = 0; q < solventclass.numAtomTypes; q++){
        //                 cuvWRK.m_data[i + grid.localDimsR[0]*j + grid.localDimsR[0]*grid.localDimsR[1]*k + grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]*q] = count;
        //                 count++;
        //             }
        //         }
        //     }
        // }

    }

    // defining destructor
    rism3d::~rism3d(){
        cout << "Object deleted!" << endl;
        // Destroy the handle
        cublasDestroy(handle);
    }
    ////////////////////////////////////////////////////////////////////////////////////////

    void rism3d :: setcut(double cut){
        pot.setcut_ljdistance(cut);

    }

    void rism3d :: createelectrondensitymap(int iv, double* electronRDF, double electronRDFGridSpacing, int totalSolventElectrons, double density, double* electronmap_ptr){
        cout << "Sorry: createElectronDensityMap not implemented, yet." << endl;
        electronmap_ptr = electronmap.m_data;
    }

    // Sets verbosity of output.
    // IN:
    //   this :: rism3d object
    //   verbosity :: 0 - no output
    //                1 - memory allocation and steps for convergence
    //                2 - 1 + convergence progress      
    // MY COMMENT: IN .F90 PUBLIC LIST
    void rism3d :: setverbosity(int verbosity){
        verbose = verbosity;
    }

    // Set parent for this timer
    // IN:
    //   this : rism3d object
    //   parent : parent timer object      
    // The timer structure is defined in the rism_timer.hpp file.
    // So, it should not be hard to add the time tracking in this 
    // code. However, I am giving priority to get it working first.
    void rism3d :: settimerparent(timer_cpp *timer){
        // cout << "Sorry: rism3d_setTimerParent not implemented, yet." << endl;
    }

    // Sets solute coordinates.
    // IN:
    //   this :: rism3d object
    //   ratu :: coordinates   
    void rism3d :: setcoord(double* solutePositions){
        soluteclass.setCoord(solutePositions);
        
        // checking X,Y,Z of first two atoms:
        // ATOM      1 HH31 ACE     1       2.000   1.000  -0.000  
        // ATOM      2  CH3 ACE     1       2.000   2.090   0.000 
        // cout << "In cpp: checking solute positions:" << endl;
        // for(int i = 0; i < 3; i++){
        //     for(int j = 0; j < soluteclass.numAtoms; j++){
        //         // cout << "posi[" << i << ", " << j << "] = " << soluteclass.position[i+3*j] << endl;
        //         cout << "   position[" << i << ", " << j << "] = " << soluteclass.position(i,j) << endl;
        //     }
        // }

    }

    bool rism3d :: calculatesolution(int ksave, int kshow, int maxSteps, bool failure, double* tolerance, int tol_size){
        bool success = true;
        GPUtype asympKSpaceTolerance_local;
        double ljTolerance_local;

        // cout << "Sorry: Still implementing calculatesolution..." << endl;
        // cout << "tol_size = " << tol_size << endl;
        // cout << "cl_size = " << cl_size << endl;
        // cout << "   Print received values:" << endl;
        // cout << "   ksave = " << ksave << endl;
        // cout << "   kshow = " << kshow << endl;
        // cout << "   maxSteps = " << maxSteps << endl;
        // cout << "   failure = " << failure << endl;
        // cout << "   tolerances: " << tolerance[0] << ", " << tolerance[1] << ", " << tolerance[2] << endl;

        // 0) Quick check that the tolerance list is of the correct length.
        if(tol_size != cl_size){
            cout << "tol_size (" << tol_size << ") is not equal to number of closures (" << cl_size << ")." << endl;
            abort();
        }
        
        // 1) Reorient solute along the principal axis and resize the grids
        // if necessary

        // Step 1 is not implemented here in the Fortan version as well.

        // 2a) Set tolerances for LJ potential and long range asymptotics
        // asympKSpaceTolerance
        // Case 1: use the  asympKSpaceTolerance as is
        // Case 2: if asympKSpaceTolerance < 0, select the asympKSpaceTolerance < tolerance
        // Case 3: if asympKSpaceTolerance > 0, use the asympKSpaceTolerance        
        if(asympKSpaceTolerance < 0.0){
            asympKSpaceTolerance_local = tolerance[tol_size-1]/10;
        }
        else{
            asympKSpaceTolerance_local = asympKSpaceTolerance;
        }

        // non-periodic case
        if(periodic == false){
            //>> Should we specify buffer > 0, instead of != 0? To avoid someone providing negative values
            // LJTolerance
            // Case 1: if ljtolerance==0, then compute the LJ interation
            //     without cutoff
            // Case 2: if ljtolerance<0 and buffer!=0, fit the LJ cutoff to fit
            //     inside the solvent box so a correction may be applied
            // Case 3: if ljtolerance<0 and buffer==0, select the ljtolerance <
            //     tolerance but >0 and set buffer large enough to fit the cutoff
            // Case 4: if ljtolerance>0 and buffer==0, set buffer large enough
            //     to fit the cutoff
            // Case 5: if ljtolerance > 0 and buffer != 0, use ljtolerance as
            //     is.  Cutoff correction may not be applied.

            if(ljTolerance < 0.0 and buffer == 0){
                ljTolerance_local = tolerance[tol_size-1]/10;
                if(verbose > 0){
                    cout << "|Setting Lennard-Jones tolerance to " << ljTolerance_local << endl;
                }
            }
            else{
                ljTolerance_local = ljTolerance;
                if(ljTolerance > 0){
                    if(verbose > 0){
                        cout << "|Using provided Lennard-Jones tolerance : " << ljTolerance_local << endl;
                    }
                }
                else if(ljTolerance == 0){
                    if(verbose > 0){
                        cout << "|No Lennard-Jones cutoff" << endl;
                    }
                    cout << "WARNING> No LJ tolerance or cutoff correction used. For more" << endl;
                    cout << "accurate calculations, increase the tolerance, box" << endl;
                    cout << "dimensions, or use buffer=0" << endl;
                }
            }

            pot.setcut_ljtolerance(ljTolerance);

        }
        
        // k-space asymptotics cutoff depends on the box size and is called
        // in resizeBox()
            
        // 2b) Get the minimum box size for this frame.
        if(periodic == true or varbox == true or *nsolution == 0){
            resizeBox();
        }

        // By default, center solute in box for the infinite dilution case,
        // but not for the periodic case since the origin helps define
        // rotation axes.
        
        centerSolute(soluteclass.translation);

        if(verbose > 0){
            cout << "||Setting solvation box to" << endl;
            
            // cout << "|grid size: " << grid.globalDimsR[0] << " X " << grid.globalDimsR[1] << " X " << grid.globalDimsR[2] << endl;
            std::cout << "|grid size:" 
                      << std::setw(11) << grid.globalDimsR[0] << " X" 
                      << std::setw(11) << grid.globalDimsR[1] << " X" 
                      << std::setw(11) << grid.globalDimsR[2] << std::endl;

            // cout << "|box size [A]: " << grid.boxLength[0] << " X " << grid.boxLength[1] << " X " << grid.boxLength[2] << endl;
            std::cout << "|box size [A]:" 
                      << std::fixed << std::setprecision(3)
                      << std::setw(12) << grid.boxLength[0] << " X" 
                      << std::setw(11) << grid.boxLength[1] << " X" 
                      << std::setw(11) << grid.boxLength[2] << std::endl;
            
            // cout << "|grid spacing [A]: " << grid.spacing[0] << " X " << grid.spacing[1] << " X " << grid.spacing[2] << endl;
            std::cout << "|grid spacing [A]:"
                      << std::fixed << std::setprecision(3)
                      << std::setw(11) << grid.spacing[0] << " X"
                      << std::setw(11) << grid.spacing[1] << " X"
                      << std::setw(11) << grid.spacing[2] << std::endl;
            
            if(periodic){
                cout << "Periodic box not supported, yet" << endl;
            }
            
            effective_buffer();
            // cout << "|effective buffer [A]: " << effect_buffer[0] << " X " << effect_buffer[1] << " X " << effect_buffer[2] << endl;
            std::cout << "|effective buffer [A]:"
                      << std::fixed << std::setprecision(3)
                      << std::setw(10) << effect_buffer[0] << " X"
                      << std::setw(11) << effect_buffer[1] << " X"
                      << std::setw(11) << effect_buffer[2] << std::endl;
        } 

        // fit the LJ cutoff inside the solvent box

        if(buffer != 0 and ljTolerance < 0 and periodic == false){
            cout << "Still need to implement the case if: buffer != 0 and ljTolerance < 0 and periodic == false" << endl;
            abort();
        }
        if(soluteclass.charged == true){
            pot.setcut_asympKTolerance(asympKSpaceTolerance_local, grid.boxVolume);
        }

        // 2) Caclualte long range aymptotics of the direct and total 
        // correlation functions about the solute
        pot.dcf_tcf_long_range_asymptotics(soluteclass.charged, periodic);

        // 3) Calculate Lennard-Jones potential about the solute
        // to do: add electrostatic potential
        pot.potential_calc();

        // 4) Propagate previously saved solute-solvent DCF solutions to
        // create an initial guess for this solution
        guessDCF();

        // 5) Calculate 3D-RISM solution using MDIIS
        // If the user provided a list of closures, use it only for
        // the first solution (nsolution == 0) or solution propagation 
        // is turned off (ncuvsteps == 0). Ohterwise, the current closure 
        // will be the last one in the list.
        // For our first implementation, we will have ncuvsteps and nsolution
        // set to zero.
        if(*nsolution == 0 || ncuvsteps == 0){
            for(int iclosure = 0; iclosure < cl_size; iclosure++){
                cout << "|Switching to " << closurelist[iclosure] << " closure" << endl;
                setclosure(closurelist[iclosure]);
                solve3DRISM(maxSteps, tolerance[iclosure]);
            }
        }

        if(rism_failure == true){
            success = false;
        }

        return success;
    }

    void rism3d :: effective_buffer(){
        for(int i = 0; i < 3; i++){
            effect_buffer[i] = min(grid.boxLength[i] - get_maxval(&soluteclass.position, soluteclass.numAtoms, i), 
                               get_minval(&soluteclass.position, soluteclass.numAtoms, i));
        }
        
    }

    void rism3d :: resizeBox(){
        int ngr[3];
        int id;
        double boxlen[3];
        int primes[4] = {2, 3, 5, 7};
        int *primes_ptr;
        bool cuv_dimension_changed;
        int ngr_size;

        primes_ptr = primes;

        if(periodic == true){
            cout << "Sorry: periodic systems not implemented, yet." << endl;
            cout << "Terminating execution..." << endl;
            abort();
        }
        else if(varbox == true){
            // Get minimum box size defined by the buffer
            getMinBoxLength(boxlen);
            // cout << "boxlen = " << boxlen[0] << " " << boxlen[1] << " " << boxlen[2] << endl;

            ngr[0] = ceil(boxlen[0]/grid.spacing[0]);
            ngr[1] = ceil(boxlen[1]/grid.spacing[1]);
            ngr[2] = ceil(boxlen[2]/grid.spacing[2]);

            boxlen[0] = ngr[0]*grid.spacing[0];
            boxlen[1] = ngr[1]*grid.spacing[1];
            boxlen[2] = ngr[2]*grid.spacing[2];

            // cout << "ngr = " << ngr[0] << " " << ngr[1] << " " << ngr[2] << endl;
            // cout << "boxlen = " << boxlen[0] << " " << boxlen[1] << " " << boxlen[2] << endl;

            // Determine if the number of grid points in each dimension
            // product only has prime factors 2, 3, 5 or 7. If not
            // increment the number of points (in that dimension) until
            // this is true.

            // Make sure that each dimension is divisible by 2 and that
            // the y- and z-dimensions are divisible by mpisize if
            // mpisize > 1.

            ngr[0] = ngr[0] + ngr[0]%2;
            ngr[1] = ngr[1] + ngr[1]%2;
            ngr[2] = ngr[2] + ngr[2]%2;

            if(mpisize > 1){
                cout << "Sorry: MPI not supported, yet." << endl;
                cout << "Terminating execution..." << endl;
                abort();
            }

            // cout << "Is factorable? " << isFactorable(15,primes_ptr,4) << endl;
            for(int i = 0; i < 3; i++){
                while(isFactorable(ngr[i],primes_ptr,4) == false){
                    ngr[i] = ngr[i] + 2;
                }
            }

            boxlen[0] = ngr[0]*grid.spacing[0];
            boxlen[1] = ngr[1]*grid.spacing[1];
            boxlen[2] = ngr[2]*grid.spacing[2];

            // cout << "ngr = " << ngr[0] << " " << ngr[1] << " " << ngr[2] << endl;
            // cout << "boxlen = " << boxlen[0] << " " << boxlen[1] << " " << boxlen[2] << endl;

            cuv_dimension_changed = false;

            for(int i = 0; i < 3; i++){
                ngr_size = merge_int(ngr[i], ngr[i]/mpisize, i != 3);
                if(cuv.m_data == nullptr || cuv.get_dims(i) != ngr_size || 
                   oldcuvChg.get_dims(i) != ngr_size ||
                   oldcuvNoChg.get_dims(i) != ngr_size){

                    cuv_dimension_changed = true;
                    break;
                }
            }
            if(cuv_dimension_changed == true){
                // cout << "Calling reallocateBox!" << endl;
                reallocateBox(ngr, grid.spacing);
            }

        }
        else{
            cout << "Sorry: fixed box size not implemented, yet." << endl;
            cout << "Terminating execution..." << endl;
            abort();
        }
    }

    void rism3d :: getMinBoxLength(double boxlen[3]){
        int iu, iv, id;
        double maxdist[3];

        // Printing some values to compare position before and after centering the solute:
        // cout << "position atom 1 before: " << soluteclass.position(0,0) << " " << soluteclass.position(1,0) << " "  << soluteclass.position(2,0) << endl;
        centerSolute(soluteclass.translation, soluteclass.origin);
        // cout << "grid.origin = " << soluteclass.origin[0] << " " << soluteclass.origin[1] << " " << soluteclass.origin[2] << endl;
        // cout << "grid.translation = " << soluteclass.translation[0] << " " << soluteclass.translation[1] << " " << soluteclass.translation[2] << endl;
        // cout << "position atom 1 after: " << soluteclass.position(0,0) << " " << soluteclass.position(1,0) << " "  << soluteclass.position(2,0) << endl;

        if(buffer > 0){
            if(verbose > 0){
                cout << "|Setting solvent box from buffer." << endl;
            }
            for(int i = 0; i < 3; i++){
                maxdist[i] = max(get_maxval(&soluteclass.position, soluteclass.numAtoms, i) + buffer,
                                 abs(get_minval(&soluteclass.position, soluteclass.numAtoms, i) - buffer));
                boxlen[i] = 2*maxdist[i];
            }

        }
        else if(buffer == 0){
            cout << "Sorry: buffer = 0 not supported, yet." << endl;
            cout << "Terminating execution..." << endl;
            abort();
        }
        else if(buffer < 0){
            cout << "Buffer must be greater or equal to 0" << endl;
            cout << "Terminating execution..." << endl;
            abort();
        }

        // Return solute to its initial position
        soluteclass.translate_solute(true);
    }

    // void rism3d :: centerSolute(GPUtype translation_o[3], GPUtype origin_o[3], 
    //                             bool calc_translation, bool default_orig){
    void rism3d :: centerSolute(GPUtype *translation_o, GPUtype *origin_o){
        GPUtype origin[3];
        double gridpoints[3];
        array_class<GPUtype> mass;
        double projection[3];
        double unitCellCenter[3] = {0,0,0};

        mass.set_memalloc(&memalloc);
        mass.alloc_mem(soluteclass.numAtoms);

        if(abs(centering) % 2 == 1){
            memcpy(mass.m_data, soluteclass.mass.m_data, soluteclass.numAtoms * sizeof(GPUtype));
        }
        else if(abs(centering) % 2 == 0){
            for(int i = 0; i < soluteclass.numAtoms; i++){
                mass(i) = 1;
            }

        }

        if(centering > 0 or (centering < 0 and nsolution == 0)){
            findCenterOfMass(&soluteclass.position, &soluteclass.centerOfMass, &mass, soluteclass.numAtoms);
            // I will leave this for now just in case we want to check center of mass for a second test system
// #if RISMCUDA_DOUBLE
//             ofstream myfile5;
//             myfile5.open("/home/fcarvalho/rism3d.cuda.test/centerofmass_db.txt");
//             for(int i = 0; i < 3; i++){
//                 myfile5 << setprecision(17) << soluteclass.centerOfMass(i) << "\n";
//             }
//             myfile5.close();
// #else
//             ofstream myfile5;
//             myfile5.open("/home/fcarvalho/rism3d.cuda.test/centerofmass.txt");
//             for(int i = 0; i < 3; i++){
//                 myfile5 << setprecision(17) << soluteclass.centerOfMass(i) << "\n";
//             }
//             myfile5.close();
// #endif // RISMCUDA_DOUBLE
            // cout << "Calculated center of mass: " << soluteclass.centerOfMass(0) << " " << soluteclass.centerOfMass(1) << " " << soluteclass.centerOfMass(2) << endl;
        }
        if(centering >= 3 and centering <= 4){
            cout << "Sorry: centering = 3 or = 4 not implemented, yet." << endl;
            abort();
        }

        if(centering != 0){
            unitCellCenter[0] = 0;
            unitCellCenter[1] = 0;
            unitCellCenter[2] = 0;
            for(int j = 0; j < 3; j++){
                unitCellCenter[0] = unitCellCenter[0] + grid.unitCellVectorsR(0,j);
                unitCellCenter[1] = unitCellCenter[1] + grid.unitCellVectorsR(1,j);
                unitCellCenter[2] = unitCellCenter[2] + grid.unitCellVectorsR(2,j);
            }
            unitCellCenter[0] = unitCellCenter[0]/2;
            unitCellCenter[1] = unitCellCenter[1]/2;
            unitCellCenter[2] = unitCellCenter[2]/2;

            if(origin_o == NULL){
                // cout << "Using unitCellCenters..." << endl;
                origin[0] = unitCellCenter[0];
                origin[1] = unitCellCenter[1];
                origin[2] = unitCellCenter[2];
            }
            else{
                // cout << "Using default origin..." << endl;
                origin[0] = origin_o[0];
                origin[1] = origin_o[1];
                origin[2] = origin_o[2];
            }
            
            if(translation_o != NULL){
                translation_o[0] = origin[0] - soluteclass.centerOfMass(0);
                translation_o[1] = origin[1] - soluteclass.centerOfMass(1);
                translation_o[2] = origin[2] - soluteclass.centerOfMass(2);
            }
            
            // cout << "====unitCellCenter: " << unitCellCenter[0] << "  " << unitCellCenter[1] << "  " << unitCellCenter[2] << endl;
            // cout << "====translation_o: " << translation_o[0] << "  " << translation_o[1] << "  " << translation_o[2] << endl;
            // cout << "====COM: " << soluteclass.centerOfMass(0) << "  " << soluteclass.centerOfMass(1) << "  " << soluteclass.centerOfMass(2) << endl;
            // cout << "=====origin2: " << origin[0] << "," << origin[1] << "," << origin[2] << endl;

            // translate solute's center of mass to the origin
            soluteclass.translate_solute();
            
        }
    }

    void rism3d :: reallocateBox(int ngr[3], GPUtype spacing[3]){
        fft.setgrid(ngr, spacing, solventclass.numAtomTypes);
        // cout << "Printing unitCellUnitVectors and ::::" << endl;
        // for(int i = 0; i < 3; i++){
        //     for(int j = 0; j < 3; j++){
        //         cout << "unitCellUnitVectorsR(" << i << "," << j << ") = " << grid.unitCellUnitVectorsR(i,j) << endl;
        //         cout << "unitCellVectorsR(" << i << "," << j << ") = " << grid.unitCellVectorsR(i,j) << endl;
        //     }
        // }

        // Moving some memory allocation here, until it is put at proper location
        // electronmap = memalloc.allocReal(electronmap,grid.localDimsR[0],grid.localDimsR[1],grid.localDimsR[2]);
        guv.set_memalloc(&memalloc, true);
        guv_dT.set_memalloc(&memalloc);
        huv.set_memalloc(&memalloc, true);
        
        guv.alloc_mem(solventclass.numAtomTypes, grid.globalDimsK[0], grid.globalDimsK[1], grid.globalDimsK[2]);
        guv_dT.alloc_mem(solventclass.numAtomTypes, grid.globalDimsK[0], grid.globalDimsK[1], grid.globalDimsK[2]);
        huv.alloc_mem(solventclass.numAtomTypes, grid.globalDimsK[0], grid.globalDimsK[1], grid.globalDimsK[2]);

        xvva.set_memalloc(&memalloc, true);
        xvva.alloc_mem(solventclass.numAtomTypes, solventclass.numAtomTypes, grid.waveNumberArraySize);

        cudaMemset(guv.m_data, 0, solventclass.numAtomTypes * grid.globalDimsK[0] * grid.globalDimsK[1] * grid.globalDimsK[2]);
        cudaMemset(xvva.m_data, 0, solventclass.numAtomTypes * solventclass.numAtomTypes * grid.waveNumberArraySize);
        
        cuvWRK.set_memalloc(&memalloc, true);
        cuvWRK.alloc_mem(NVec, solventclass.numAtomTypes, grid.globalDimsK[0], grid.globalDimsK[1], grid.globalDimsK[2]);

        cuvresWRK.set_memalloc(&memalloc, true);
        cuvresWRK.alloc_mem(NVec, solventclass.numAtomTypes, grid.globalDimsK[0], grid.globalDimsK[1], grid.globalDimsK[2]);       
        
        // Set residue arrays to zero
        cudaMemset(cuvresWRK.m_data, 0, NVec * solventclass.numAtomTypes * grid.globalDimsK[0] * grid.globalDimsK[1] * grid.globalDimsK[2] * sizeof(GPUtype));

        // setting some important data in MDIIS class
        mdiisclass.resize(cuvWRK.m_data, cuvresWRK.m_data, solventclass.numAtomTypes, grid.globalDimsK[0], grid.globalDimsK[1], grid.globalDimsK[2]);

        // Writing xvv to a file: I will keep this for now
// #if RISMCUDA_DOUBLE
//         ofstream myfile1;
//         myfile1.open("/home/fcarvalho/rism3d.cuda.test/xvv_db.txt");
//         for(int i = 0; i < solventclass.numAtomTypes; i++){
//             for(int j = 0; j < solventclass.numAtomTypes; j++){
//                 for(int k = 0; k < solventclass.numRDFpoints; k++){
//                     myfile1 << setprecision(17) << i << "," << j << "," << k << "," << solventclass.xvv(i,j,k) << "\n";
//                 }
//             }
//         }
//         myfile1.close();

//         // myfile1.open("/home/fcarvalho/rism3d.cuda.test/xvv_db_new.txt");
//         // for(int i = 0; i < solventclass.numAtomTypes; i++){
//         //     for(int j = 0; j < solventclass.numAtomTypes; j++){
//         //         for(int k = 0; k < solventclass.numRDFpoints; k++){
//         //             int idxvv = i*solventclass.numAtomTypes*solventclass.numRDFpoints + j*solventclass.numRDFpoints + k;
//         //             myfile1 << setprecision(17) << idxvv << " " << solventclass.xvv.m_data[idxvv] << "\n";
//         //         }
//         //     }
//         // }
//         // myfile1.close();
// #else
//         ofstream myfile1;
//         myfile1.open("/home/fcarvalho/rism3d.cuda.test/xvv.txt");
//         for(int i = 0; i < solventclass.numAtomTypes; i++){
//             for(int j = 0; j < solventclass.numAtomTypes; j++){
//                 for(int k = 0; k < solventclass.numRDFpoints; k++){
//                     myfile1 << setprecision(17) << i << "," << j << "," << k << "," << solventclass.xvv(i,j,k) << "\n";
//                 }
//             }
//         }
//         myfile1.close();
// #endif //RISMCUDA_DOUBLE

    // Interpolating susceptibility
    interpolateSolventSusceptibility();
        
    }

    void rism3d :: interpolateSolventSusceptibility(){
        int maxPointsToInterp = 5;
        int numPointsToInterp = grid.waveNumberArraySize;
        int numRDFpoints = solventclass.numRDFpoints;

        int device = -1;
        cudaGetDevice(&device);

        cudaMemPrefetchAsync(solventclass.waveNumbers.m_data, numRDFpoints*sizeof(GPUtype),device,NULL);
        cudaMemPrefetchAsync(solventclass.xvv.m_data, numRDFpoints*solventclass.numAtomTypes*solventclass.numAtomTypes*sizeof(GPUtype),device,NULL);
        cudaMemPrefetchAsync(grid.waveNumbers.m_data, grid.waveNumberArraySize*sizeof(GPUtype),device,NULL);

        // Here I think we can just reserve the ammount of memory for each pair of solvent atom types,
        // since we do not need memory beyond that for each kernel call and each kernel call will start at
        // the correct contigous memory address.
        cudaMemPrefetchAsync(xvva.m_data, grid.waveNumberArraySize*solventclass.numAtomTypes*solventclass.numAtomTypes*sizeof(GPUtype),device,NULL);

        for(int i = 0; i < solventclass.numAtomTypes; i++){
            for(int j = 0; j < solventclass.numAtomTypes; j++){
                // Old column major ordering
                // cu_polinomialInterpolation(maxPointsToInterp, numPointsToInterp, 
                //                    numRDFpoints, solventclass.numAtomTypes, 
                //                    solventclass.waveNumbers.m_data, 
                //                    solventclass.xvv.m_data + i*numRDFpoints + j*numRDFpoints*solventclass.numAtomTypes, 
                //                    grid.waveNumbers.m_data, 
                //                    xvva.m_data, i, j, solventclass.numAtomTypes);

                // We have changed xvv to be row major, this is the new version
                cu_polinomialInterpolation(maxPointsToInterp, numPointsToInterp, 
                                   numRDFpoints, solventclass.numAtomTypes, 
                                   solventclass.waveNumbers.m_data, 
                                   solventclass.xvv.m_data + i*numRDFpoints*solventclass.numAtomTypes + j*numRDFpoints, 
                                   grid.waveNumbers.m_data, 
                                   xvva.m_data, i, j, solventclass.numAtomTypes);
            }
        }

        cudaDeviceSynchronize();

        // Writing interpolated xvv (xvva) to a file: I will keep this for now
// #if RISMCUDA_DOUBLE
//         ofstream myfile1;
//         myfile1.open("/home/fcarvalho/rism3d.cuda.test/xvva_db.txt");
//         // for(int i = 0; i < grid.waveNumberArraySize*solventclass.numAtomTypes*solventclass.numAtomTypes; i++){
//         //     myfile1 << setprecision(17) << i << "," << xvva.m_data[i] << "\n";
//         // }
//         for(int i = 0; i < solventclass.numAtomTypes; i++){
//             for(int j = 0; j < solventclass.numAtomTypes; j++){
//                 for(int k = 0; k < grid.waveNumberArraySize; k++){
//                     // myfile1 << setprecision(17) << i << ", " << j << ", " << k << " :: " << xvva(i,j,k) << "\n";
//                     int idwn = i*solventclass.numAtomTypes*grid.waveNumberArraySize + j*grid.waveNumberArraySize + k;
//                     myfile1 << setprecision(17) << idwn << " " << xvva.m_data[idwn] << "\n";
//                 }
//             }
//         }
//         myfile1.close();
// #else
//         ofstream myfile1;
//         myfile1.open("/home/fcarvalho/rism3d.cuda.test/xvva.txt");
//         // for(int i = 0; i < grid.waveNumberArraySize*solventclass.numAtomTypes*solventclass.numAtomTypes; i++){
//         //     myfile1 << setprecision(17) << i << "," << xvva(i) << "\n";
//         // }
//         for(int i = 0; i < solventclass.numAtomTypes; i++){
//             for(int j = 0; j < solventclass.numAtomTypes; j++){
//                 for(int k = 0; k < grid.waveNumberArraySize; k++){
//                     myfile1 << setprecision(17) << i << ", " << j << ", " << k << " :: " << xvva(i,j,k) << "\n";
//                 }
//             }
//         }
//         myfile1.close();
// #endif //RISMCUDA_DOUBLE
    }

    // We need to copy data over for the case of ncuvsteps > 0. Since it is not the case right now,
    // I will just skip this implementation for now
    void rism3d :: guessDCF(){
        // cout << "Nothing is being doing in guessDCF now, since we are using ncuvsteps = 0." << endl;
    }

    void rism3d :: setclosure(string type){
        if(closure){
            delete closure;
        }

        size_t ispse = type.find("PSE", 0, 3);

        if(type == "KH"){
            // cout << "| Setting " << type << endl;
            rism3d_closure* closure_child = new kh(&solventclass, &soluteclass, &grid, &pot, &memalloc);
            closure = closure_child;
        }
        else if(ispse != string::npos){
            cout << "| Setting pse-" << type[3] << endl;
            rism3d_closure* closure_child = new rism3d_closure_psen(&solventclass, &soluteclass, &grid, &pot, &memalloc);
            closure = closure_child;
        }
        else if(type == "HNC"){
            cout << "| setting hnc" << endl;
            rism3d_closure* closure_child = new rism3d_closure_hnc(&solventclass, &soluteclass, &grid, &pot, &memalloc);
            closure = closure_child;
        }
        else{
            cout << "amber_rism_interface should check the closure names" << endl;
            cout << "If this is being displayed, there is something wrong" << endl;
            cout << "or it supports other closures (besides HNC, PSE-n and KH)" << endl;
            cout << "that it is not supported in the c++ version" << endl;
            cout << "Not a valid closure!" << endl;
            abort();
        }


    }

    void rism3d :: solve3DRISM(int maxSteps, double tolerance){
        auto start = high_resolution_clock::now();

        if(*nsolution == 0 || ncuvsteps == 0){
            cudaMemset(cuvWRK.m_data, 0, NVec*solventclass.numAtomTypes*grid.globalDimsK[0]*grid.globalDimsK[1]*grid.globalDimsK[2]*sizeof(GPUtype));

            if(soluteclass.charged == true and periodic == false){
                // The first thing in the next call is to subtract these added values. Thus,
                // initial guess is still zero
                set_dcf_longrange_cu(cuvWRK.m_data, solventclass.charge.m_data, pot.dcfLongRangeAsympR.m_data);
            }
        }

        // Getting just the first work vector initially
        cuv.m_data = cuvWRK.m_data;
        cuvres.m_data = cuvresWRK.m_data;

        // Saving values to compare with fortran version: i am keeping this for now
// #if RISMCUDA_DOUBLE
//         ofstream file1("/home/fcarvalho/rism3d.cuda.test.chg/cuv_ini_1_db.txt");
//         ofstream file2("/home/fcarvalho/rism3d.cuda.test.chg/cuv_ini_2_db.txt");
// #else
//         ofstream file1("/home/fcarvalho/rism3d.cuda.test.chg/cuv_ini_1_float.txt");
//         ofstream file2("/home/fcarvalho/rism3d.cuda.test.chg/cuv_ini_2_float.txt");
// #endif // RISMCUDA_DOUBLE
//         file1 << std::scientific << std::setprecision(16);
//         file2 << std::scientific << std::setprecision(16);
//         for(int ix = 0; ix < grid.globalDimsR[0]; ix++){
//             for(int iy = 0; iy < grid.globalDimsR[1]; iy++){
//                 for(int iz = 0; iz < grid.globalDimsR[2]; iz++){
//                     file1 << cuv.m_data[0 * grid.globalDimsR[0] * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2) + 
//                                         ix * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2) + 
//                                         iy * (grid.globalDimsR[2] + 2) + 
//                                         iz] << endl;
//                     file2 << cuv.m_data[1 * grid.globalDimsR[0] * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2) + 
//                                         ix * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2) + 
//                                         iy * (grid.globalDimsR[2] + 2) + 
//                                         iz] << endl;
//                 }
//             }
//         }
//         file1.close();
//         file2.close();

        // prefetch data used for solving rism that have not been sent to GPU memory yet
        int device = -1;
        cudaGetDevice(&device);

        // DO I NEED TO PREFETCH CUV AND CUVRES? THEY ARE POINTERS TO DATA ALREADY IN
        // THE GPU MEMORY
        cudaMemPrefetchAsync(cuvresWRK.m_data, NVec*solventclass.numAtomTypes*grid.globalDimsK[0]*grid.globalDimsK[1]*grid.globalDimsK[2]*sizeof(GPUtype), device, NULL);
        cudaMemPrefetchAsync(cuvWRK.m_data, NVec*solventclass.numAtomTypes*grid.globalDimsK[0]*grid.globalDimsK[1]*grid.globalDimsK[2]*sizeof(GPUtype), device, NULL);
        cudaMemPrefetchAsync(cuvres.m_data, solventclass.numAtomTypes*grid.globalDimsK[0]*grid.globalDimsK[1]*grid.globalDimsK[2]*sizeof(GPUtype), device, NULL);
        cudaMemPrefetchAsync(cuv.m_data, solventclass.numAtomTypes*grid.globalDimsK[0]*grid.globalDimsK[1]*grid.globalDimsK[2]*sizeof(GPUtype), device, NULL);
        cudaMemPrefetchAsync(guv.m_data, solventclass.numAtomTypes*grid.globalDimsK[0]*grid.globalDimsK[1]*grid.globalDimsK[2]*sizeof(GPUtype), device, NULL);

        converged = false;
        for(int i = 0; i < maxSteps; i++){
            single3DRISMsolution(tolerance);
            
            if(verbose == 2){
                std::cout << " Step=" << std::setw(5) << i + 1 << "     Resid= ";
                if (std::abs(mdiisclass.residual) >= 1e-1) {
                    std::cout << std::fixed << std::setprecision(3) << mdiisclass.residual;
                    std::cout << "         IS=" << std::setw(3) << mdiisclass.getIS() << std::endl;
                } else {
                    std::cout << std::scientific << std::setprecision(3) << std::uppercase << mdiisclass.residual;
                    std::cout << "     IS=" << std::setw(3) << mdiisclass.getIS() << std::endl;
                }
            }

            if(mdiisclass.residual < tolerance){
                converged = true;
                cout << "|RXRISM converged in" << std::setw(6) << i + 1 << " steps" << endl;
                break;
            }
        }

        if(converged == true){
            rism_failure = false;
        }
        else{
            rism_failure = true;
            if(isnan(mdiisclass.residual)){
                cout << "RXRISM: NaN residual detected, aborting." << endl;
                abort();
            }
            else{
                cout << "RXRISM: reached limit number of relaxation steps: " << maxSteps << endl;
                abort();
            }
        }

        // Because we are doing the Fourier transform in place, and on C++ side the padded
        // values of huv, cuv and guv are not all at the end of the array, 
        // setting the paddings to zero here is a safety measure to not sum 
        // non-zero values while carrying out the integrals
        set_padding2zero(guv.m_data);
        set_padding2zero(huv.m_data);
        set_padding2zero(cuv.m_data);

        auto stop = high_resolution_clock::now();
        std::chrono::duration<double> duration = stop - start;

        // Comment line below to hide time to solve 3D-RISM equation. This will not
        // be necessary anymore after implementing the time tracking.
        if(verbose == 2){
            cout << "Time for solving 3D-RISM equations: " << duration.count() << endl;
        }
    }

    // Function: single3DRISMsolution
    // This function will use cuda libraries and custom kernels
    // to get one iteration of MDIIS method done.
    // Things to be done:
    //     - All data manipulation must be done in the GPU to avoid sending data back and forth between
    //       CPU and GPU;
    //     - So far it is solving MDIIS only for NVec = 1. So, it must be generalized.
    void rism3d :: single3DRISMsolution(double tolerance){
        // --------------------------------------------------------------
        // Subtract short-range part from Cuv(r) (if not periodic);
        // Cuv(r) is then loaded into the guv array.
        // >> It is not handling charged species, yet. Thus, it will only
        // >> load cuv data into guv array.
        // >> To do: add long range subtraction and turn cu_memcpy into
        // >> cu_subtractLongRange
        // --------------------------------------------------------------
        if(soluteclass.charged == true && periodic == false){
            subtract_dcf_longrange_cu(guv.m_data, cuv.m_data, solventclass.charge.m_data, pot.dcfLongRangeAsympR.m_data);

            // testing with serial code
            // int Nx = grid.globalDimsR[0];
            // int Ny = grid.globalDimsR[1];
            // int Nz = grid.globalDimsR[2];

            // for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            //     for(int ix = 0; ix < Nx; ix++){
            //         for(int iy = 0; iy < Ny; iy++){
            //             for(int iz = 0; iz < Nz; iz++){
            //                 guv.m_data[iv*Nx*Ny*(Nz+2) + ix*Ny*(Nz+2) + iy*(Nz+2) + iz] = cuv.m_data[iv*Nx*Ny*(Nz+2) + ix*Ny*(Nz+2) + iy*(Nz+2) + iz] - 
            //                                                                               solventclass.charge.m_data[iv] * pot.dcfLongRangeAsympR.m_data[ix*Ny*(Nz) + iy*(Nz) + iz];
            //             }
            //         }

            //     }
            // }
        }
        else{
            cu_memcpy(cuv.m_data, guv.m_data, solventclass.numAtomTypes * grid.totalLocalPointsK);
        }

        // Ensure that guv paddings are zero
        // It seems to be not necessary: I've ran with this commented and
        // results were the same
        // Old serial version:
        // for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
        //     for(int ix = 0; ix < grid.globalDimsR[0]; ix++){
        //         for(int iy = 0; iy < grid.globalDimsR[1]; iy++){
        //             for(int iz = 0; iz < 2; iz++){
        //                 guv.m_data[iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
        //                           ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
        //                           iy*(grid.globalDimsR[2] + 2) + grid.globalDimsR[2] + iz] = 0;
        //             }
        //         }
        //     }
        // } 
        // CUDA version:
        // set_padding2zero();    
        set_padding2zero(guv.m_data);

        // --------------------------------------------------------------
        // [Short-range part of] Cuv(r) FFT>K.
        // --------------------------------------------------------------
        fft.create_plans(grid.globalDimsR[0], grid.globalDimsR[1], grid.globalDimsR[2]);

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            GPUtype *data2transform = (guv.m_data + iv * grid.globalDimsK[0]*grid.globalDimsK[1]*grid.globalDimsK[2]);
#if defined(RISMCUDA_DOUBLE)
            cufftExecD2Z(fft.plan, reinterpret_cast<double*>(data2transform), reinterpret_cast<cufftDoubleComplex*>(data2transform));
#else
            cufftExecR2C(fft.plan, reinterpret_cast<float*>(data2transform), reinterpret_cast<cufftComplex*>(data2transform));
#endif

            cudaDeviceSynchronize();
        }
    
        
        // Normalize data and DO NOT take the complex conjugate
        // Old serial version:
        // for(int i = 0; i < solventclass.numAtomTypes * grid.totalLocalPointsK; i++){
        //     guv.m_data[i] = guv.m_data[i] / (grid.globalDimsR[0] * grid.globalDimsR[1] * grid.globalDimsR[2]);
        // }

        // New CUDA version:
        double normFactor = 1.0 / (grid.globalDimsR[0] * grid.globalDimsR[1] * grid.globalDimsR[2]);
        GPUtype normalization = static_cast<GPUtype>(normFactor);
#if RISMCUDA_DOUBLE
        cublasDscal(handle, solventclass.numAtomTypes * grid.totalLocalPointsK, &normalization, guv.m_data, 1);
#else
        cublasSscal(handle, solventclass.numAtomTypes * grid.totalLocalPointsK, &normalization, guv.m_data, 1);
#endif //RISMCUDA_DOUBLE
        cudaDeviceSynchronize();

        // Convert to numerical recipes storage scheme
        // Old serial version:
        // convert2nr();

        // New CUDA version
        convert2nr_cu();

        // --------------------------------------------------------------
        // Add long-range part to Cuv(k) in K-space.
        // --------------------------------------------------------------
        if(soluteclass.charged == true && periodic == false){
            add_dcf_longrange_cu(guv.m_data, solventclass.charge.m_data, pot.dcfLongRangeAsympK.m_data);
        }

        // --------------------------------------------------------------
        // Huv(k) by RISM.
        // --------------------------------------------------------------

        // Old serial version
        // int iga;
        // for(int iv1 = 0; iv1 < solventclass.numAtomTypes; iv1++){
        //     for(int ig1 = 0; ig1 < grid.totalLocalPointsK; ig1++){
        //         huv.m_data[iv1 * grid.totalLocalPointsK + ig1] = 0;
        //         iga = grid.waveVectorWaveNumberMap.m_data[ig1/2];
        //         for(int iv2 = 0; iv2 < solventclass.numAtomTypes; iv2++){
        //             huv.m_data[iv1 * grid.totalLocalPointsK + ig1] = huv.m_data[iv1 * grid.totalLocalPointsK + ig1] + 
        //                                                              guv.m_data[iv2 * grid.totalLocalPointsK + ig1] * 
        //                                                              xvva.m_data[iv1 * solventclass.numAtomTypes * grid.waveNumberArraySize + 
        //                                                                          iv2 * grid.waveNumberArraySize + iga];
        //         }
        //     }
        // }

        // New CUDA version
        get_h_k();
        
        // I have moved the function call to get back to FFTW format
        // inside the if statment below. It should be called after subtracting the
        // long range asymp, because the long range calculations assume numerical recipes
        // formatting.
        // ! Remember to include it for the periodic case when implementing it!
        if(!(periodic)){
            // --------------------------------------------------------------
            // Add long-range part of Huv(k) at k = 0, which was estimated by
            // long-range part of Cuv(k) at k = 0.
            // --------------------------------------------------------------
            if(soluteclass.charged == true){
                //! As in potential class, while handling huv0k: there is just a few datapoints here that maybe it is better not
                //  to lauch a cuda kernel for this.
                for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
                    // Real
                    huv.m_data[iv * grid.totalLocalPointsK + 0] = huv.m_data[iv * grid.totalLocalPointsK + 0] + pot.huvk0.m_data[iv*2 + 0];

                    // Imaginary
                    huv.m_data[iv * grid.totalLocalPointsK + 1] = huv.m_data[iv * grid.totalLocalPointsK + 1] + pot.huvk0.m_data[iv*2 + 1];
                }
            }

            // --------------------------------------------------------------
            // Subtract long-range part from huv in K-space.
            // --------------------------------------------------------------
            if(solventclass.ionic == true && soluteclass.charged == true){
                subtract_tcf_longrange_cu(huv.m_data, solventclass.charge_sp.m_data, pot.tcfLongRangeAsympK.m_data);
            }

            // Convert back to FFTW storage scheme
            // Old serial version:
            // convert2fftw();

            // New CUDA version:
            convert2fftw_cu();

            
        }else{
            // ---------------------------------------------------------------
            // Remove the background charge effect for periodic calculations
            // ---------------------------------------------------------------
            cout << "Periodic systems not supported, yet." << endl;
            abort();
            // something like this:
            //     if (this%mpirank == 0 .and. this%solute%charged) then
            //         this%huv(:2, :) = this%huv(:2, :) + this%potential%huvk0(:, :)
            //     endif
        }

        // --------------------------------------------------------------
        // Short-range part of Huv(k) FFT>R.
        // --------------------------------------------------------------

        for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
            GPUtype *data2transform_inv = (huv.m_data + iv * grid.totalLocalPointsK);
#if defined(RISMCUDA_DOUBLE)
            cufftExecZ2D(fft.plan_inv, reinterpret_cast<cufftDoubleComplex*>(data2transform_inv), reinterpret_cast<double*>(data2transform_inv));
#else
            cufftExecC2R(fft.plan_inv, reinterpret_cast<cufftComplex*>(data2transform_inv), reinterpret_cast<float*>(data2transform_inv));
#endif
        }
        cudaDeviceSynchronize();

        // --------------------------------------------------------------
        // To do: Add long-range part to huv in R-space.
        // --------------------------------------------------------------
        if(solventclass.ionic == true && periodic == false && soluteclass.charged == true){
            add_tcf_longrange_cu(huv.m_data, solventclass.charge_sp.m_data, pot.tcfLongRangeAsympR.m_data);
        }

//         // print huv in real space for testing: I will keep this here for now
// #if RISMCUDA_DOUBLE
//         ofstream file3("/home/fcarvalho/rism3d.cuda.ionic.test.chg/huv_r_3_db.txt");
//         ofstream file4("/home/fcarvalho/rism3d.cuda.ionic.test.chg/huv_r_4_db.txt");
// #else
//         ofstream file3("/home/fcarvalho/rism3d.cuda.ionic.test.chg/huv_r_3_float.txt");
//         ofstream file4("/home/fcarvalho/rism3d.cuda.ionic.test.chg/huv_r_4_float.txt");
// #endif // RISMCUDA_DOUBLE
//         file3 << std::scientific << std::setprecision(16);
//         file4 << std::scientific << std::setprecision(16);
//         for(int ix = 0; ix < grid.globalDimsR[0]; ix++){
//             for(int iy = 0; iy < grid.globalDimsR[1]; iy++){
//                 for(int iz = 0; iz < grid.globalDimsR[2]; iz++){
//                     file3 << ix << ", " << iy << ", " << iz << ", " << huv.m_data[2 * grid.globalDimsR[0] * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2) + 
//                                         ix * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2) + 
//                                         iy * (grid.globalDimsR[2] + 2) + 
//                                         iz] << endl;
//                     file4 << ix << ", " << iy << ", " << iz << ", " << huv.m_data[3 * grid.globalDimsR[0] * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2) + 
//                                         ix * grid.globalDimsR[1] * (grid.globalDimsR[2] + 2) + 
//                                         iy * (grid.globalDimsR[2] + 2) + 
//                                         iz] << endl;
//                 }
//             }
//         }
//         file3.close();
//         file4.close();
//         abort();

        // --------------------------------------------------------------
        // Solve the closure for the RDF.
        // -------------------------------------------------------------

        closure->guv(pot.uuv.m_data, guv.m_data, huv.m_data, cuv.m_data, solventclass.numAtomTypes * grid.totalLocalPointsK);

        // --------------------------------------------------------------
        // Calculate TCF residual for use in estimating DCF residual.
        // For now I will leave it using all the indexes explicitly,
        // just to make sure everything is being added correctly.
        // Later I will got back to the first implementation**
        // --------------------------------------------------------------

        // Calculate residual vector and update cuv using mdiis_vec = 0.5
        // Old serial version:
        // for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
        //     for(int ix = 0; ix < grid.globalDimsR[0]; ix++){
        //         for(int iy = 0; iy < grid.globalDimsR[1]; iy++){
        //             for(int iz = 0; iz < grid.globalDimsR[2]; iz++){
        //                 cuvres.m_data[iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
        //                               ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
        //                               iy*(grid.globalDimsR[2] + 2) + iz] = 
        //                               (guv.m_data[iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
        //                               ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
        //                               iy*(grid.globalDimsR[2] + 2) + iz] - 1) - 
        //                               huv.m_data[iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
        //                               ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
        //                               iy*(grid.globalDimsR[2] + 2) + iz];

        //                 // res = res + abs(cuvres.m_data[iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
        //                 //               ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
        //                 //               iy*(grid.globalDimsR[2] + 2) + iz]);

        //                 // cuv.m_data[iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
        //                 //           ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
        //                 //           iy*(grid.globalDimsR[2] + 2) + iz] = 
        //                 //           cuv.m_data[iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
        //                 //           ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
        //                 //           iy*(grid.globalDimsR[2] + 2) + iz] + 
        //                 //           0.5 * cuvres.m_data[iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
        //                 //           ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
        //                 //           iy*(grid.globalDimsR[2] + 2) + iz];
        //             }
        //         }
        //     }
        // }

        // New CUDA version
        get_residue();
        
        // --------------------------------------------------------------
        // MDIIS
        // --------------------------------------------------------------
        mdiisclass.advance(tolerance);

        // Getting next work vector
        cuv.m_data = cuvWRK.m_data + (mdiisclass.getWRKvec() - 1) * solventclass.numAtomTypes*grid.globalDimsK[0]*grid.globalDimsK[1]*grid.globalDimsK[2];
        cuvres.m_data = cuvresWRK.m_data + (mdiisclass.getWRKvec() - 1) * solventclass.numAtomTypes*grid.globalDimsK[0]*grid.globalDimsK[1]*grid.globalDimsK[2];

        // I am keeping this write part here for now in case we want to check things in the near future
        // ofstream myfile1;
// #if RISMCUDA_DOUBLE
//         myfile1.open("/home/fcarvalho/rism3d.cuda.test/cuv_afterAdvance.txt");
//         for(int iv = 0; iv < solventclass.numAtomTypes; iv++){
//             for(int ix = 0; ix < grid.globalDimsR[0]; ix++){
//                 for(int iy = 0; iy < grid.globalDimsR[1]; iy++){
//                     for(int iz = 0; iz < grid.globalDimsR[2]; iz++){
//                         myfile1 << setprecision(17) << cuv.m_data[iv*grid.globalDimsR[0]*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) + 
//                                       ix*grid.globalDimsR[1]*(grid.globalDimsR[2] + 2) +
//                                       iy*(grid.globalDimsR[2] + 2) + iz] << "\n";
//                     }
//                 }
//             }
//         }
//         // for(int i = 0; i < solventclass.numAtomTypes * grid.totalLocalPointsK; i++){
//         //     myfile1 << setprecision(17) << i << " :: " << guv.m_data[i] << "\n";
//         // }
//         myfile1.close();
// #endif //RISMCUDA_DOUBLE

        fft.destroy_plans();
    }

    void rism3d :: force(double* ff){
        cout << "Sorry: force not implemented, yet." << endl;
    }

    double rism3d :: excesschemicalpotential_tot(bool o_lr){
        double tot = 0;
        int size = grid.localDimsK[0] * grid.localDimsK[1] * grid.localDimsK[2];
        double *excesschempot = closure->aexcessChemicalPotential(huv.m_data, cuv.m_data, size, o_lr);
        for(int i = 0; i < solventclass.numAtomTypes; i++){
            tot = tot + excesschempot[i];
        }

        return tot;
    }

    double* rism3d :: excesschemicalpotential_tot_map(int* dim1, int* dim2, int* dim3){
        cout << "Sorry: excessChemicalPotential_tot_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        // double* excessChemicalPotential_t_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]];
        double* excessChemicalPotential_t_m;
        excessChemicalPotential_t_m = memalloc.allocReal(excessChemicalPotential_t_m, grid.localDimsR[0], grid.localDimsR[1], grid.localDimsR[2]);
        return excessChemicalPotential_t_m;
    }

    double* rism3d :: excesschemicalpotentialgf_tot_map(int* dim1, int* dim2, int* dim3){
        cout << "Sorry: excessChemicalPotentialGF_tot_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        double* excessChemicalPotentialGF_t_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]];
        return excessChemicalPotentialGF_t_m;
    }

    double* rism3d :: excesschemicalpotentialpcplus_tot_map(int* dim1, int* dim2, int* dim3){
        cout << "Sorry: excessChemicalPotentialPCPLUS_tot_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        double* excessChemicalPotentialPCPLUS_t_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]];
        return excessChemicalPotentialPCPLUS_t_m;
    }

    double* rism3d :: excesschemicalpotentialuc_tot_map(int* dim1, int* dim2, int* dim3, double* coeff){
        cout << "Sorry: excessChemicalPotentialUC_tot_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        double* excessChemicalPotentialUC_t_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]];
        return excessChemicalPotentialUC_t_m;
    }

    double* rism3d :: solventpotene_tot_map(int* dim1, int* dim2, int* dim3){
        cout << "Sorry: solventPotEne_tot_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        double* solventPotEne_t_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]];
        return solventPotEne_t_m;
    }

    double* rism3d :: solvationenergy_tot_map(int* dim1, int* dim2, int* dim3){
        cout << "Sorry: solvationEnergy_tot_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        double* solvationEnergy_t_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]];
        return solvationEnergy_t_m;
    }

    double* rism3d :: solvationenergygf_tot_map(int* dim1, int* dim2, int* dim3){
        cout << "Sorry: solvationEnergyGF_tot_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        double* solvationEnergyGF_t_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]];
        return solvationEnergyGF_t_m;
    }

    double* rism3d :: solvationenergypcplus_tot_map(int* dim1, int* dim2, int* dim3){
        cout << "Sorry: solvationEnergyPCPLUS_tot_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        double* solvationEnergyPCPLUS_t_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]];
        return solvationEnergyPCPLUS_t_m;
    }

    double* rism3d :: solvationenergyuc_tot_map(int* dim1, int* dim2, int* dim3, double* coeff){
        cout << "Sorry: solvationEnergyUC_tot_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        double* solvationEnergyUC_t_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]];
        return solvationEnergyUC_t_m;
    }

    double* rism3d :: excesschemicalpotential_site_map(int* dim1, int* dim2, int* dim3, int* dim4){
        cout << "Sorry: excessChemicalPotential_site_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        *dim4 = solventclass.numAtomTypes;
        double* excessChemicalPotential_s_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*
                                                         grid.localDimsR[2]*solventclass.numAtomTypes];
        return excessChemicalPotential_s_m;
    }

    double* rism3d :: excesschemicalpotentialgf_site_map(int* dim1, int* dim2, int* dim3, int* dim4){
        cout << "Sorry: excessChemicalPotentialGF_site_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        *dim4 = solventclass.numAtomTypes;
        double* excessChemicalPotentialgf_s_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*
                                                         grid.localDimsR[2]*solventclass.numAtomTypes];
        return excessChemicalPotentialgf_s_m;
    }

    double* rism3d :: solventpotene_site_map(int* dim1, int* dim2, int* dim3, int* dim4){
        cout << "Sorry: solventPotEne_site_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        *dim4 = solventclass.numAtomTypes;
        double* solventPotEne_s_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*
                                                         grid.localDimsR[2]*solventclass.numAtomTypes];
        return solventPotEne_s_m;
    }

    double* rism3d :: solvationenergy_site_map(int* dim1, int* dim2, int* dim3, int* dim4){
        cout << "Sorry: solvationEnergy_site_map not implemented, yet." << endl;
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        *dim4 = solventclass.numAtomTypes;
        double* solvationEnergy_s_m = new double[grid.localDimsR[0]*grid.localDimsR[1]*
                                                         grid.localDimsR[2]*solventclass.numAtomTypes];
        return solvationEnergy_s_m;
    }

    double* rism3d :: excesschemicalpotential(int* len, bool o_lr){
        int size = grid.localDimsK[0] * grid.localDimsK[1] * grid.localDimsK[2];
        *len = solventclass.numAtomTypes;
        return closure->aexcessChemicalPotential(huv.m_data, cuv.m_data, size, o_lr);
    }

    double* rism3d :: excesschemicalpotentialgf(int* len, bool o_lr){
        cout << "Sorry: excessChemicalPotentialGF not implemented, yet." << endl;
        cout << "   o_lr = " << o_lr << endl;
        double* test;
        test = closure->rism3d_closure_aexcessChemicalPotentialGF();
        *len = solventclass.numAtomTypes;
        return test;
    }

    double rism3d :: excesschemicalpotentialpcplus(bool o_lr){
        cout << "Sorry: excessChemicalPotentialPCPLUS not implemented, yet." << endl;
        cout << "   o_lr = " << o_lr << endl;
        // MUST CALL FUNCTION FROM RISM3D_CLOSURE CLASS
        return 1234.123;
    }

    double rism3d :: excesschemicalpotentialuc(double* coeff, bool o_lr){
        cout << "Sorry: excessChemicalPotentialUC not implemented, yet." << endl;
        cout << "   o_lr = " << o_lr << endl;
        return 12.21;
    }

    double* rism3d :: solventpotene(int* len){
        int size = grid.localDimsK[0] * grid.localDimsK[1] * grid.localDimsK[2];
        // double* solvPotEne;
        // solvPotEne = closure->solvPotEne();
        *len = solventclass.numAtomTypes;
        return closure->solvPotEne(guv.m_data, pot.uuv.m_data, size);
    }

    double rism3d :: partialmolarvolume(){
        // cout << "Sorry: still implementing partialmolarvolume." << endl;
        return closure->partialMolarVolume(cuv.m_data);
    }

    double* rism3d :: excessparticles(int* len, bool o_lr){
        *len = solventclass.numAtomTypes;
        if(o_lr == true and solventclass.ionic == true){
            return closure->aexcessParticles(huv.m_data);
        }
        else{
            return closure->excessParticles(huv.m_data);
        }
    }

    // std::vector<double> rism3d :: excessparticles(bool o_lr){
    //     cout << "Sorry: still implementing excessParticles." << endl;
    //     double* test;
    //     test = closure->rism3d_closure_excessPart();
    //     vector<double> test_arr(solventclass.numAtomTypes);
    //     for(int i = 0; i < solventclass.numAtomTypes; i++){
    //         test_arr[i] = test[i];
    //     }
    //     return test_arr;
    // }

    // double* rism3d :: kirkwoodbuff(int* len, bool o_lr){
    //     cout << "Sorry: kirkwoodBuff not implemented, yet." << endl;
    //     cout << "   o_lr = " << o_lr << endl;
    //     double* test;
    //     test = closure->rism3d_closure_kirkwoodBuff();
    //     *len = solventclass.numAtomTypes;
    //     return test;
    // }
    // std::vector<double> rism3d :: kirkwoodbuff(bool o_lr){
    double* rism3d :: kirkwoodbuff(int* len, bool o_lr){
        *len = solventclass.numAtomTypes;
        return closure->kirkwoodBuff(huv.m_data, o_lr);
    }

    // double* rism3d :: dcfintegral(int* len){
    //     cout << "Sorry: DCFintegral not implemented, yet." << endl;
    //     double* test;
    //     test = closure->rism3d_closure_DCFintegral();
    //     *len = solventclass.numAtomTypes;
    //     return test;
    // }
    // std::vector<double> rism3d :: dcfintegral(){
    double* rism3d :: dcfintegral(int* len){
        *len = solventclass.numAtomTypes;
        // cout << "Sorry: DCFintegral not implemented, yet." << endl;
        // double* test = new double[solventclass.numAtomTypes];
        // // test = closure->rism3d_closure_DCFintegral();
        // for(int i = 0; i < solventclass.numAtomTypes; i++){
        //     test[i] = i + 666;
        // }
        return closure->get_cuv_int();
    }

    bool rism3d :: cancalc_dt(){
        cout << "Sorry: canCalc_DT not implemented, yet." << endl;
        return false;
    }

    bool rism3d :: calculatesolution_dt(int kshow, int maxSteps, bool failure, double tolerance){
        cout << "Sorry: calculateSolution_dT not implemented, yet." << endl;
        cout << "   kshow = " << kshow << endl;
        cout << "   maxSteps = " << maxSteps << endl;
        cout << "   rism_failure = " << rism_failure << endl;
        cout << "   tolerance = " << tolerance << endl;
        return false;
    }

    double* rism3d :: solvationenergy(int* len, bool o_lr){
        cout << "Sorry: solvationEnergy not implemented, yet." << endl;
        cout << "   o_lr = " << o_lr << endl;
        double* test;
        test = closure->rism3d_closure_solvationEnergy();
        *len = solventclass.numAtomTypes;
        return test;
    }

    double* rism3d :: solvationenergygf(int* len, bool o_lr){
        cout << "Sorry: solvationEnergyGF not implemented, yet." << endl;
        cout << "   o_lr = " << o_lr << endl;
        double* test;
        test = closure->rism3d_closure_solvationEnergyGF();
        *len = solventclass.numAtomTypes;
        return test;
    }

    double* rism3d :: solvationenergygf_site_map(int* dim1, int* dim2, int* dim3, int* dim4){
        // *dim1 = grid.localDimsR[0];
        // *dim2 = grid.localDimsR[1];
        // *dim3 = grid.localDimsR[2];
        // *dim4 = solventclass.numAtomTypes;
        // closure->rism3d_closure_solvationEnergyGF_site_map(&grid,huv,huv_dT,cuv->m_data,cuv_dT);
        // return closure->solvationEnergyGF_site_map_1d;
        return NULL;
    }

    double* rism3d :: excessparticles_dt(int* len, bool o_lr){
        cout << "Sorry: excessParticles_dT not implemented, yet." << endl;
        cout << "   o_lr = " << o_lr << endl;
        double* test;
        test = closure->rism3d_closure_excessParticles_dT();
        *len = solventclass.numAtomTypes;
        return test;
    }

    double* rism3d :: kirkwoodbuff_dt(int* len, bool o_lr){
        cout << "Sorry: kirkwoodBuff_dT not implemented, yet." << endl;
        cout << "   o_lr = " << o_lr << endl;
        double* test;
        test = closure->rism3d_closure_kirkwoodBuff_dT();
        *len = solventclass.numAtomTypes;
        return test;
    }

    double* rism3d :: dcfintegral_dt(int* len){
        cout << "Sorry: DCFintegral_dT not implemented, yet." << endl;
        double* test;
        test = closure->rism3d_closure_DCFintegral_dT();
        *len = solventclass.numAtomTypes;
        return test;
    }

    double rism3d :: solvationenergypcplus(bool o_lr){
        cout << "Sorry: solvationEnergyPCPLUS not implemented, yet." << endl;
        cout << "   o_lr = " << o_lr << endl;
        return 76.43;
    }

    double rism3d :: solvationenergyuc(double* coeff, bool o_lr){
        cout << "Sorry: solvationEnergyUC not implemented, yet." << endl;
        cout << "   o_lr = " << o_lr << endl;
        return 21.12;
    }

    double rism3d :: partialmolarvolume_dt(){
        cout << "Sorry: partialMolarVolume_dT not implemented, yet." << endl;
        return 11.112;
    }

    void rism3d :: unsetcharges(){
        cout << "Sorry: unsetCharges not implemented, yet." << endl;
    }

    void rism3d :: resetcharges(){
        cout << "Sorry: resetCharges not implemented, yet." << endl;
    }

    void rism3d :: map_site_to_site_flat(double* thermo_map_flat, int center_site){
        cout << "Sorry: map_site_to_site_flat not implemented, yet." << endl;
    }

    // 3D array version of rism3d_map_site_to_site
    void rism3d :: map_site_to_site_3D(double* thermo_map, int center_site){
        cout << "Sorry: map_site_to_site_3D not implemented, yet." << endl;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// END... FUNCTIONS CALLED FROM AMBER_RISM_INTERFACE.F90
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// BEGIN... ATTRIBUTES CALLED FROM AMBER_RISM_INTERFACE.F90
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    void rism3d :: set_nsolution(int nsol){
        *nsolution = nsol;
    }

    int rism3d :: get_nsolution(){
        return *nsolution;
    }

    double* rism3d :: get_guv(int *dim1, int *dim2, int arg1, int arg2){
        *dim1 = 1;
        *dim2 = 1;
        arg1 = arg1 - 1;
        arg2 = arg2 - 1;

#if RISMCUDA_DOUBLE
        return &guv.m_data[arg1 + grid.totalLocalPointsK*arg2];
#else
        return (double*)&guv.m_data[arg1 + grid.totalLocalPointsK*arg2];
#endif // RISMCUDA_DOUBLE
        // return &guv.m_data[arg1 + grid.totalLocalPointsK*arg2];
        // return NULL;
    }
    double* rism3d :: get_guv(int *dim1, int *dim2, int arg){
        *dim1 = grid.totalLocalPointsK;
        *dim2 = 1;
        arg = arg - 1;

        double *temp_arr;
        temp_arr = memalloc.allocReal(temp_arr, grid.totalLocalPointsK);

        for(int i = 0; i < grid.totalLocalPointsK; i++){
            temp_arr[i] = guv.m_data[grid.totalLocalPointsK*arg + i];
        }

        return temp_arr;
        // return NULL;
    }

    double* rism3d :: get_huv(int *dim1, int *dim2, int arg1, int arg2){
        *dim1 = 1;
        *dim2 = 1;
        arg1 = arg1 - 1;
        arg2 = arg2 - 1;
#if RISMCUDA_DOUBLE
        return &huv.m_data[arg1 + grid.totalLocalPointsK*arg2];
#else
        return (double*)&huv.m_data[arg1 + grid.totalLocalPointsK*arg2];
#endif // RISMCUDA_DOUBLE
        // return &huv.m_data[arg1 + grid.totalLocalPointsK*arg2];
        // return NULL;
    }
    double* rism3d :: get_huv(int *dim1, int *dim2, int arg){
        *dim1 = grid.totalLocalPointsK;
        *dim2 = 1;
        arg = arg - 1;

        double *temp_arr;
        temp_arr = memalloc.allocReal(temp_arr, grid.totalLocalPointsK);

        for(int i = 0; i < grid.totalLocalPointsK; i++){
            temp_arr[i] = huv.m_data[grid.totalLocalPointsK*arg + i];
        }

        return temp_arr;
        // return NULL;
    }

    double rism3d :: get_cuv(int arg1, int arg2, int arg3, int arg4){
        arg1 = arg1 - 1;
        arg2 = arg2 - 1;
        arg3 = arg3 - 1;
        arg4 = arg4 - 1;

        return cuv.m_data[arg1 + grid.localDimsR[0]*arg2 + grid.localDimsR[0]*grid.localDimsR[1]*arg3 + grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]*arg4];
    }
    double* rism3d :: get_cuv(int *dim1, int *dim2, int *dim3, int *dim4, int arg){
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        *dim4 = 1;
        arg = arg - 1;

        double *temp_arr = new double [grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]];
        for(int i = 0; i < grid.localDimsR[0]; i++){
            for(int j = 0; j < grid.localDimsR[1]; j++){
                for(int k = 0; k < grid.localDimsR[2]; k++){
                    temp_arr[i + grid.localDimsR[0]*j + grid.localDimsR[0]*grid.localDimsR[1]*k] = cuv.m_data[i + grid.localDimsR[0]*j + grid.localDimsR[0]*grid.localDimsR[1]*k + grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]*arg];
                }
            }
        }

        return temp_arr;
    }

    // double* rism3d :: get_guv_dt(int *dim1, int *dim2){
    //     *dim1 = grid.totalLocalPointsK;
    //     *dim2 = solventclass.numAtomTypes;
    //     return guv_dT;
    // }
    double* rism3d :: get_guv_dt(int *dim, int arg){
        *dim = grid.totalLocalPointsK;
        arg = arg - 1;
        double *temp_arr = new double [grid.totalLocalPointsK];
        for(int i = 0; i < grid.totalLocalPointsK; i++){
            temp_arr[i] = guv_dT.m_data[solventclass.numAtomTypes*i + arg];
        }
        return temp_arr;
        // return NULL;
    }
    double* rism3d :: get_guv_dt(int *dim, int arg1, int arg2){
        *dim = 1;
        arg1 = arg1-1;
        arg2 = arg2-1;
        return &guv_dT.m_data[grid.totalLocalPointsK*arg1 + arg2];
        // return NULL;
    }

    double* rism3d :: get_electronmap(int *dim1, int *dim2, int *dim3){
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        return electronmap.m_data;
    }

    bool rism3d :: get_periodic(){
        return periodic;
    }

    bool rism3d :: get_failure(){
        return rism_failure;
    }

    /////////////////// From grid class
    double* rism3d :: get_spacing(int* len){
        *len = 3;
#if RISMCUDA_DOUBLE
        return grid.spacing;
#else
        return (double*)grid.spacing;
#endif // RISMCUDA_DOUBLE
    }

    double* rism3d :: get_voxelvectorsr(int* len1, int* len2){
        *len1 = 3;
        *len2 = 3;
        return grid.voxelVectorsR.m_data;
    }

    double* rism3d :: get_unitcellangles(int* len){
        *len = 3;
        return grid.unitCellAngles;
    }

    double* rism3d :: get_unitcellvectorsk(int* len1, int* len2){
        *len1 = 3;
        *len2 = 3;
        return grid.unitCellVectorsK;
    }

    double* rism3d :: get_unitcellvectorsr(int* len1, int* len2){
        *len1 = 3;
        *len2 = 3;
        return grid.unitCellVectorsR.m_data;
    }

    double* rism3d :: get_boxlength(int* len){
        *len = 3;
#if RISMCUDA_DOUBLE
        return grid.boxLength;
#else
        return (double*)grid.boxLength;
#endif // RISMCUDA_DOUBLE
    }

    int rism3d :: get_totallocalpointsr(){
        return grid.totalLocalPointsR;
    }

    double rism3d :: get_voxelvolume(){
        return grid.voxelVolume;
    }

    int* rism3d :: get_localdimsr(int* len){
        *len = 3;
        return grid.localDimsR;
    }

    int rism3d :: get_localdimsr(int indx){
        if (indx < 1 || indx > 3){
            cout << "Index must be > 0 and < 4" << endl;
            return grid.localDimsR[0];
        } else { 
            return grid.localDimsR[indx - 1];
        }
    }

    int* rism3d :: get_globaldimsr(int* len){
        *len = 3;
        return grid.globalDimsR;
    }

    int rism3d :: get_globaldimsr(int indx){
        if (indx < 1 || indx > 3){
            cout << "Index must be > 0 and < 4" << endl;
            return grid.globalDimsR[0];
        } else {
            return grid.globalDimsR[indx - 1];
        }
    }

    ////////////////// From solvent struct
    // void rism3d :: get_atomname(char** names, int* names_len, int str_len){
    //     *names_len = solventclass.numAtomTypes;
    //     char *buf = new char[solventclass.numAtomTypes*str_len];
    //     for(int i = 0; i < solventclass.numAtomTypes; i++){
    //         if(i == 0){
    //             strcpy(buf,const_cast<char*>(solventclass.atomName[i].c_str()));
    //         }
    //         else{
    //             strcat(buf,const_cast<char*>(solventclass.atomName[i].c_str()));
    //         }
    //     }
    //     *names = buf;
    // }
    void rism3d :: get_atomname(string** names, int* names_len){
        *names_len = solventclass.numAtomTypes;
        *names = solventclass.atomName;
    }

    double rism3d :: get_xikt_dt(){
        return solventclass.xikt_dT;
    }

    double rism3d :: get_temperature(){
        return solventclass.temperature;
    }
    // rism_3d%solvent%atomName // this one is pending...

    int rism3d :: get_numatomtypes(){
        return solventclass.numAtomTypes;
    }

    double* rism3d :: get_charge(int *len){
        *len = solventclass.numAtomTypes;
#if RISMCUDA_DOUBLE
        return solventclass.charge.m_data;;
#else
        for (int i = 0; i < *len; ++i) {
            solventclass.charge_db.m_data[i] = static_cast<double>(solventclass.charge.m_data[i]);
        }
        return solventclass.charge_db.m_data;
#endif // RISMCUDA_DOUBLE
    }

    double rism3d :: get_charge(int indx){
        return solventclass.charge.m_data[indx - 1];
    }

    double* rism3d :: get_density(int *len){
        *len = solventclass.numAtomTypes;
        // return solventclass.density;
#if defined(RISMCUDA_DOUBLE)
        return solventclass.density.m_data;
#else
        return solventclass.density_db.m_data;
#endif //RISMCUDA_DOUBLE
    }
    double rism3d :: get_density(int indx){
        // return solventclass.density[indx-1];
#if defined(RISMCUDA_DOUBLE)
        return solventclass.density.m_data[indx-1];
#else
        return solventclass.density_db.m_data[indx-1];
#endif //RISMCUDA_DOUBLE
    }

    bool rism3d :: get_ionic(){
        return solventclass.ionic;
    }

    // bool periodic;

    ////////////////// From solvent struct

    double* rism3d :: get_centerofmass(int* len){
        *len = 3;
#if RISMCUDA_DOUBLE
        return soluteclass.centerOfMass.m_data;
#else
        return (double*)soluteclass.centerOfMass.m_data;
#endif //RISMCUDA_DOUBLE
        // return soluteclass.centerOfMass;
    }

    double* rism3d :: get_translation(int* len){
        *len = 3;
#if RISMCUDA_DOUBLE
        return soluteclass.translation;
#else
        return (double*)soluteclass.translation;
#endif //RISMCUDA_DOUBLE
    }

    bool rism3d :: get_charged(){
        return soluteclass.charged;
    }

    int rism3d :: get_numatoms(){
        return soluteclass.numAtoms;
    }

    double* rism3d :: get_mass(int* dim1){
        *dim1 = soluteclass.numAtoms;
#if RISMCUDA_DOUBLE
        return soluteclass.mass.m_data;
#else
        return (double*)soluteclass.mass.m_data;
#endif //RISMCUDA_DOUBLE
    }

    ////////////////// From Pot class

    double* rism3d :: get_tcflongrangeasympr(int *dim1){
        *dim1 = grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2];
#if RISMCUDA_DOUBLE
        return pot.tcfLongRangeAsympR.m_data;
#else
        // Allocate memmory for pot.dcfLongRangeAsympR_db only if necessy
        pot.tcfLongRangeAsympR_db.alloc_mem(grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]);
        for(int i = 0; i < grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]; i++){
            pot.tcfLongRangeAsympR_db.m_data[i] = static_cast<double>(pot.tcfLongRangeAsympR.m_data[i]);
        }
        return pot.tcfLongRangeAsympR_db.m_data;
#endif //RISMCUDA_DOUBLE
    }

    double* rism3d :: get_dcflongrangeasympr(int *dim1){
        *dim1 = grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2];
#if RISMCUDA_DOUBLE
        return pot.dcfLongRangeAsympR.m_data;
#else
        // Allocate memmory for pot.dcfLongRangeAsympR_db only if necessy
        pot.dcfLongRangeAsympR_db.alloc_mem(grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]);
        for(int i = 0; i < grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]; i++){
            pot.dcfLongRangeAsympR_db.m_data[i] = static_cast<double>(pot.dcfLongRangeAsympR.m_data[i]);
        }
        return pot.dcfLongRangeAsympR_db.m_data;
#endif //RISMCUDA_DOUBLE
    }

    double rism3d :: get_dcflongrangeasympr(int indx){
#if RISMCUDA_DOUBLE
        return pot.dcfLongRangeAsympR.m_data[indx-1];
#else
        return static_cast<double>(pot.dcfLongRangeAsympR.m_data[indx-1]);
#endif //RISMCUDA_DOUBLE
    }

    double* rism3d :: get_uuv(int *dim1, int *dim2, int *dim3, int *dim4, int indx){
        *dim1 = grid.localDimsR[0];
        *dim2 = grid.localDimsR[1];
        *dim3 = grid.localDimsR[2];
        *dim4 = 1;

        double *temp_arr = new double [grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]];
        indx = indx - 1;
        for(int i = 0; i < grid.localDimsR[0]; i++){
            for(int j = 0; j < grid.localDimsR[1]; j++){
                for(int k = 0; k < grid.localDimsR[2]; k++){
                    temp_arr[i + grid.localDimsR[0]*j + grid.localDimsR[0]*grid.localDimsR[1]*k] = pot.uuv.m_data[i + grid.localDimsR[0]*j + grid.localDimsR[0]*grid.localDimsR[1]*k + grid.localDimsR[0]*grid.localDimsR[1]*grid.localDimsR[2]*indx];
                }
            }
        }

        return temp_arr;
    }

    void rism3d :: selftest(){
        cout << "Sorry: selftest not implemented, yet" << endl;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// END... ATTRIBUTES CALLED FROM AMBER_RISM_INTERFACE.F90
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    void rism3d :: setbox_variable(rism3d_grid *grid, double o_buffer, double* o_grdspc){
        grid->set_spacing(o_grdspc);
        buffer = o_buffer;
        varbox = true;
    }

    void rism3d :: setbox_fixed(double* o_boxlen, int* o_ng3){
        cout << "Setting boxlen and ng3" << endl;
        for(int i = 0; i < 3; i++){
            fixedNumGridPoints[i] = o_ng3[i];
            fixedBoxDimensionsR[i] = o_boxlen[i];
        }
    }

    // SetClosureList: will get list of closures to be used from Fortran and asign to closurelist
    void rism3d :: set_closurelist(char **names, int *N)
    {
        cl_size = *N;

        // Checking number of closures supported: right now, it is working only with KH
        if(cl_size > 1){
            cout << "Only KH closure is supported at this moment" << endl;
            abort();
        }

        closurelist = memalloc.allocString(cl_size, closurelist);
        for(int i = 0; i < cl_size; i++){
            closurelist[i] = names[i];
        }
    }

    // solvent rism3d :: set_solvent(solvent *arg, solvent solv){
    //     solv = *arg;
    //     // ALLOCATE MEMORY TO ALL POINTERS
    //     // AND INCLUDE ALL COPY LOOPS

    //     // solv.atomMultiplicity = new int[solv.numAtomTypes];
    //     solv.atomMultiplicity = memalloc.allocInt(solv.atomMultiplicity, solv.numAtomTypes);
    //     for(int i = 0; i < solv.numAtomTypes; i++){
    //         solv.atomMultiplicity[i] = arg->atomMultiplicity[i];
    //     }

    //     // solv.numAtoms = new int[solv.numMolecules];
    //     solv.numAtoms = memalloc.allocInt(solv.numAtoms, solv.numMolecules);
    //     for(int i = 0; i < solv.numMolecules; i++){
    //         solv.numAtoms[i] = arg->numAtoms[i];
    //     }
        
    //     //> This one is set in a specif function!!!
    //     solv.atomName = new string[solv.numAtomTypes];

    //     // solv.waveNumbers = new double[solv.numRDFpoints];
    //     solv.waveNumbers = memalloc.allocReal(solv.waveNumbers, solv.numRDFpoints);
    //     for(int i = 0; i < solv.numRDFpoints; i++){
    //         solv.waveNumbers[i] = arg->waveNumbers[i];
    //     }


    //     // solv.xvv = new double[solv.numRDFpoints*solv.numAtomTypes*solv.numAtomTypes];
    //     solv.xvv = memalloc.allocReal(solv.xvv, solv.numRDFpoints, solv.numAtomTypes, solv.numAtomTypes);
    //     for(int i = 0; i < solv.numRDFpoints; i++){
    //         for(int j = 0; j < solv.numAtomTypes; j++){
    //             for(int k = 0; k < solv.numAtomTypes; k++){
    //                 solv.xvv[i + solv.numRDFpoints*j + solv.numRDFpoints*solv.numAtomTypes*k] = 
    //                     arg->xvv[i + solv.numRDFpoints*j + solv.numRDFpoints*solv.numAtomTypes*k];
    //             }
    //         }
    //     }

    //     // solv.xvv_dT = new double[solv.numRDFpoints*solv.numAtomTypes*solv.numAtomTypes];
    //     solv.xvv_dT = memalloc.allocReal(solv.xvv_dT, solv.numRDFpoints, solv.numAtomTypes, solv.numAtomTypes);
    //     for(int i = 0; i < solv.numRDFpoints; i++){
    //         for(int j = 0; j < solv.numAtomTypes; j++){
    //             for(int k = 0; k < solv.numAtomTypes; k++){
    //                 solv.xvv_dT[i + solv.numRDFpoints*j + solv.numRDFpoints*solv.numAtomTypes*k] = 
    //                     arg->xvv_dT[i + solv.numRDFpoints*j + solv.numRDFpoints*solv.numAtomTypes*k];
    //             }
    //         }
    //     }


    //     // solv.charge = new double[solv.numAtomTypes];
    //     solv.charge = memalloc.allocReal(solv.charge, solv.numAtomTypes);
    //     for(int i = 0; i < solv.numAtomTypes; i++){
    //         solv.charge[i] = arg->charge[i];
    //     }
        
    //     // solv.charge_sp = new double[solv.numAtomTypes];
    //     solv.charge_sp = memalloc.allocReal(solv.charge_sp, solv.numAtomTypes);
    //     for(int i = 0; i < solv.numAtomTypes; i++){
    //         solv.charge_sp[i] = arg->charge_sp[i];
    //     }

    //     // solv.density = new double[solv.numAtomTypes];
    //     solv.density = memalloc.allocReal(solv.density, solv.numAtomTypes);
    //     for(int i = 0; i < solv.numAtomTypes; i++){
    //         solv.density[i] = arg->density[i];
    //     }
        
    //     // solv.density_sp = new double[solv.numMolecules];
    //     solv.density_sp = memalloc.allocReal(solv.density_sp, solv.numMolecules);
    //     for(int i = 0; i < solv.numMolecules; i++){
    //         solv.density_sp[i] = arg->density_sp[i];
    //     }
        
    //     // solv.ljSigma = new double[solv.numAtomTypes];
    //     solv.ljSigma = memalloc.allocReal(solv.ljSigma, solv.numAtomTypes);
    //     for(int i = 0; i < solv.numAtomTypes; i++){
    //         solv.ljSigma[i] = arg->ljSigma[i];
    //     }
        
    //     // solv.ljEpsilon = new double[solv.numAtomTypes];
    //     solv.ljEpsilon = memalloc.allocReal(solv.ljEpsilon, solv.numAtomTypes);
    //     for(int i = 0; i < solv.numAtomTypes; i++){
    //         solv.ljEpsilon[i] = arg->ljEpsilon[i];
    //     }

    //     int max_atomMultiplicity = *max_element(solv.atomMultiplicity, solv.atomMultiplicity + solv.numAtomTypes);
    //     int max_numAtoms = *max_element(solv.numAtoms, solv.numAtoms + solv.numMolecules);
    //     // solv.coord = new double[3*max_atomMultiplicity*max_numAtoms*solv.numMolecules];
    //     solv.coord = memalloc.allocReal(solv.coord, 3, max_atomMultiplicity, max_numAtoms, solv.numMolecules);
    //     for(int i = 0; i < 3; i++){
    //         for(int j = 0; j < max_atomMultiplicity; j++){
    //             for(int k = 0; k < max_numAtoms; k++){
    //                 for(int q = 0; q < solv.numMolecules; q++){
    //                     solv.coord[i + 3*j + 3*max_atomMultiplicity*k + 3*max_atomMultiplicity*max_numAtoms*q] = 
    //                         arg->coord[i + 3*j + 3*max_atomMultiplicity*k + 3*max_atomMultiplicity*max_numAtoms*q];
    //                 }
    //             }
    //         }
    //     }

    //     // solv.background_correction = new double[solv.numAtomTypes];
    //     solv.background_correction = memalloc.allocReal(solv.background_correction, solv.numAtomTypes);
    //     for(int i = 0; i < solv.numAtomTypes; i++){
    //         solv.background_correction[i] = arg->background_correction[i];
    //     }

    //     // solv.delhv0 = new double[solv.numAtomTypes];
    //     solv.delhv0 = memalloc.allocReal(solv.delhv0, solv.numAtomTypes);
    //     for(int i = 0; i < solv.numAtomTypes; i++){
    //         solv.delhv0[i] = arg->delhv0[i];
    //     }

    //     // cout << "Accuse Segmentation fault" << endl;
    //     // solv.delhv0_dT = new double[solv.numAtomTypes];
    //     solv.delhv0_dT = memalloc.allocReal(solv.delhv0_dT, solv.numAtomTypes);
    //     for(int i = 0; i < solv.numAtomTypes; i++){
    //         solv.delhv0_dT[i] = arg->delhv0_dT[i];
    //     }

    //     cout << "Solvent set" << endl;
    //     return solv;
    // }

    // solute rism3d :: set_solute(solute *arg, solute solu){
    //     solu = *arg;
    //     // NEED TO ALLOCATE MEMORY TO ALL POINTERS
    //     // NEED TO INCLUDE ALL COPY LOOPS

    //     // solu.mass = new double[solu.numAtoms];
    //     solu.mass = memalloc.allocReal(solu.mass, solu.numAtoms);
    //     for(int i = 0; i < solu.numAtoms; i++){
    //         solu.mass[i] = arg->mass[i];
    //     }

    //     // solu.charge = new double[solu.numAtoms];
    //     solu.charge = memalloc.allocReal(solu.charge, solu.numAtoms);
    //     for(int i = 0; i < solu.numAtoms; i++){
    //         solu.charge[i] = arg->charge[i];
    //     }

    //     // solu.origCharge = new double[solu.numAtoms];
    //     // for(int i = 0; i < solu.numAtoms; i++){
    //     //     solu.origCharge[i] = arg->origCharge[i];
    //     // }

    //     // solu.ljSigma = new double[solu.numAtoms];
    //     solu.ljSigma = memalloc.allocReal(solu.ljSigma, solu.numAtoms);
    //     for(int i = 0; i < solu.numAtoms; i++){
    //         solu.ljSigma[i] = arg->ljSigma[i];
    //     }

    //     // solu.ljEpsilon = new double[solu.numAtoms];
    //     solu.ljEpsilon = memalloc.allocReal(solu.ljEpsilon, solu.numAtoms);
    //     for(int i = 0; i < solu.numAtoms; i++){
    //         solu.ljEpsilon[i] = arg->ljEpsilon[i];
    //     }

    //     // solu.centerOfMass = new double[3];
    //     solu.centerOfMass = memalloc.allocReal(solu.centerOfMass, 3);
    //     for(int i = 0; i < 3; i++){
    //         solu.centerOfMass[i] = arg->centerOfMass[i];
    //     }

    //     // solu.translation = new double[3];
    //     solu.translation = memalloc.allocReal(solu.translation, 3);
    //     for(int i = 0; i < 3; i++){
    //         solu.translation[i] = arg->translation[i];
    //     }

    //     // solu.position = new double[3*solu.numAtoms];
    //     solu.position = memalloc.allocReal(solu.position, 3, solu.numAtoms);
    //     for(int i = 0; i < 3; i++){
    //         for(int j = 0; j < solu.numAtoms; j++){
    //             solu.position[i+3*j] = arg->position[i+3*j];
    //         }
    //     }

    //     cout << "Solute set" << endl;
    //     return solu;
    // }

    void rism3d :: PrintTests(){
        // cout << "================= " << "PrintTests " << "=================" << endl;
        // cout << "buffer = " << buffer << endl;
        // cout << "Closures: " << endl;
        // // for(int i = 0; i < cl_size; i++){
        // //     cout << closurelist[i] << endl;
        // // }
        // cout << centering << endl;
        // cout << varbox << endl;
        // cout << periodic << endl;
        // cout << periodicPotential << endl;

        // cout << "Testing structures..." << endl;
        // cout << "   Solvent" << endl;
        // cout << "numAtomTypes = " << solventclass.numAtomTypes << ", T = " << solventclass.temperature << endl;
        // cout << "densities = " << solventclass.density[0] << " and " << solventclass.density[1] << endl;
        // cout << "charges = " << solventclass.charge[0] << " and " << solventclass.charge[1] << endl;
        // cout << "xikt_dT = " << solventclass.xikt_dT << endl;
        // cout << "Atom names: " << endl;
        // // for(int i = 0; i < solventclass.numAtomTypes; i++){
        // //     cout << solventclass.atomName[i] << "..." << endl;
        // // }
        // cout << "Atom multiplicity: " << solventclass.atomMultiplicity[0] << " and " << solventclass.atomMultiplicity[1] << endl;
        
        // int max_atomMultiplicity = *max_element(solventclass.atomMultiplicity, solventclass.atomMultiplicity + solventclass.numAtomTypes);
        // int max_numAtoms = *max_element(solventclass.numAtoms, solventclass.numAtoms + solventclass.numMolecules);
        // cout << "max_atomMultiplicity, max_numAtoms and solventclass.numMolecules = "  << max_atomMultiplicity << " and " << max_numAtoms << " and " << solventclass.numMolecules << endl;
        // cout << "Coord: " << endl;
        // for(int i = 0; i < 3; i++){
        //     for(int j = 0; j < max_atomMultiplicity; j++){
        //         for(int k = 0; k < max_numAtoms; k++){
        //             for(int q = 0; q < solventclass.numMolecules; q++){
        //                 cout << "coord[" << i << ", " << j << ", " << k << ", " << q << "] = " << solventclass.coord[i + 3*j + 3*max_atomMultiplicity*k + 3*max_atomMultiplicity*max_numAtoms*q] << endl;
        //             }
        //         }
        //     }
        // }

        // cout << "xvv test: " << endl;
        // for(int i = 0; i < solventclass.numRDFpoints; i++){
        //     for(int j = 0; j < solventclass.numAtomTypes; j++){
        //         for(int k = 0; k < solventclass.numAtomTypes; k++){
        //             cout << "xvv[" << i << ", " << j << ", " << k << "] = " << solventclass.xvv[i + solventclass.numRDFpoints*j + solventclass.numRDFpoints*solventclass.numAtomTypes*k] << endl;
        //         }
        //     }
        // }

        // // char* part;
        // // memcpy(part, solventclass.atomName + 4 /* Offset */, 8 /* Length */);
        // // part[8] = 0; /* Add terminator */
        // // cout << part << "..." << endl; //", " << solventclass.atomName[1] <<

        // cout << "   Solute" << endl;
        // cout << soluteclass.numAtoms << ", " << soluteclass.mass[0] << " and " << soluteclass.mass[1] << endl;
        // for(int i = 0; i < 3; i++){
        //     for(int j = 0; j < soluteclass.numAtoms; j++){
        //         // cout << "posi[" << i << ", " << j << "] = " << soluteclass.position[i+3*j] << endl;
        //         cout << "posi[" << i << ", " << j << "] = " << soluteclass.position(i,j) << endl;
        //     }
        // }

        // cout << "Guv:" << endl;
        // for(int i = 0; i < grid.totalLocalPointsK; i++){
        //     for(int j = 0; j < solventclass.numAtomTypes; j++){
        //         cout << "guv[" << i << ", " << j << "] = " << guv[i + grid.totalLocalPointsK*j] << endl;
        //         cout << "   Global idx: " << i + grid.totalLocalPointsK*j << endl;
        //     }
        // }

        // cout << "In cpp: Voxel volume: " << grid.voxelVolume << endl;

        // // CALLING FUNCTIONS FROM POTENTIAL CLASS HERE TO TEST ON FORTRAN SIDE
        // pot.long_range_asymptotics();
        // cout << "Calling Pot method potential_calc" << endl;
        // pot.potential_calc(solventclass.numAtomTypes);
        
        // cout << "==================" << "===========" << "=================" << endl;
    }

    void rism3d :: test_access(){
        cout << "Testing access to centering: " << centering << endl;
        cout << "Testing access to pot.cutoff: "  << pot.cutoff << endl;
        cout << "Testing access to pot.cutoff2: "  << pot.cutoff2 << endl;
        cout << "Testing access to pot.ljCutoffs2[numAtoms-1,numAtomTypes-1]: "  << pot.ljCutoffs2(soluteclass.numAtoms -1, solventclass.numAtomTypes-1) << endl;
    }

    ////////////////////////////////////////////// FUNCTIONS FOR WRITING VOLUMETRIC DATA

    // So far, these functions will be called for each atom type
    void rism3d :: opendx_write_cpp(int atomType, string *file){
        // Uncomment these lines below to check the atomType index and file name:
        // cout << "atomType number = " << atomType << endl;
        // cout << "file name: " << *file << endl;
        // if(soluteclass.charged == true){
#if defined(CUDA_NOSWAP)
            dxclass.write_dx(file, grid.localDimsR[0], grid.localDimsR[1], grid.localDimsR[2], 
                             grid.spacing[0], grid.spacing[1], grid.spacing[2], soluteclass.translation,
                             guv.m_data + atomType * grid.localDimsR[0] * grid.localDimsR[1] * (grid.localDimsR[2] + 2));
#else
            dxclass.write_swapped_dx(file, grid.localDimsR[0], grid.localDimsR[1], grid.localDimsR[2], 
                             grid.spacing[0], grid.spacing[1], grid.spacing[2], soluteclass.translation,
                             guv.m_data + atomType * grid.localDimsR[0] * grid.localDimsR[1] * (grid.localDimsR[2] + 2));
#endif
        //     } else{
        //     dxclass.write_dx(file, grid.localDimsR[0], grid.localDimsR[1], grid.localDimsR[2], 
        //                      grid.spacing[0], grid.spacing[1], grid.spacing[2], soluteclass.translation,
        //                      guv.m_data + atomType * grid.localDimsR[0] * grid.localDimsR[1] * (grid.localDimsR[2] + 2));
        // }
    }

    void rism3d :: mrc_map_write_cpp(int atomType, string *file){
        cout << "Sorry: mrc_map_write_cpp not implemented, yet. No file will be generated" << endl;
    }

    void rism3d :: xyzv_write_cpp(int atomType, string *file){
        cout << "Sorry: xyzv_write_cpp not implemented, yet. No file will be generated" << endl;
    }

}