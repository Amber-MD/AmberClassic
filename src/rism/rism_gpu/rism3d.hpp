#ifndef CLASSES_HPP
#define CLASSES_HPP

#include <iostream>
#include "rism3d_closure.hpp"
#include "rism3d_closure_kh.hpp"
#include "rism3d_closure_hnc.hpp"
#include "rism3d_closure_psen.hpp"
#include "rism3d_potential.hpp"
#include "rism3d_grid.hpp"
#include "rism3d_fft.hpp"
#include "rism3d_solvent_class.hpp"
#include "rism3d_solute_class.hpp"
#include "rism3d_solute.hpp"
#include "rism3d_solvent.hpp"
#include "rism_timer.hpp"
#include "rism3d_safemem.hpp"
#include "mdiis.hpp"
#include "opendx.hpp"
#include <bits/stdc++.h>
#include <cstring>
#include <vector>
#include <cufft.h>
// #include "array_class.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////                           //////////////////////////////////////
///////////////////////////////////////   PRE-DEFINED VARIABLES   //////////////////////////////////////
///////////////////////////////////////                           //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

static double pd_grdspc[3] = {123, 234, 345};
static double pd_boxlen[3] = {0, 0, 0};
static int pd_ng3[3] = {0, 0, 0};
static char pd_periodic[] = " ";
static double pd_unitCellDimensions[6] = {0, 0, 0, 0, 0, 0};

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////                          //////////////////////////////////////
///////////////////////////////////////     DEFINE NAMESPACE     //////////////////////////////////////
///////////////////////////////////////                          //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

namespace rism3d_c {
    static timer_cpp parent;
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////                      //////////////////////////////////////
    ///////////////////////////////////////     DEFINE CLASS     //////////////////////////////////////
    ///////////////////////////////////////                      //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////// 

    class rism3d
    {

      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////             /////////////////////////////////////////////
      //////////////////////////////////////////////   PRIVATE   /////////////////////////////////////////////
      //////////////////////////////////////////////             /////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////

    private:
      // Create class objects needed here
      rism3d_safemem memalloc;
      rism3d_solute_class soluteclass;
      rism3d_solvent_class solventclass;
      rism3d_grid grid;
      rism3d_fft fft;
      rism3d_potential pot;
      mdiis mdiisclass;
      opendx dxclass;
      // rism3d_closure closure;
      rism3d_closure *closure;

      ///////////////////////////////////////////////////////////////////////////////
      //////// VARIABLES DEFINED ONLY ON C++ SIDE (NOT PASSED TO .YAML FILE) ////////
      ///////////////////////////////////////////////////////////////////////////////

      ///////////////////////////////////////////////////////////////////////////////
      ///////////////////////////// VARIABLE ATTRIBUTES /////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////

      // cuBLAS handle and status
      cublasHandle_t handle;
      cublasStatus_t status;

      // Create the class atribute to store set of closures to be used
      string *closurelist;
      int cl_size; // number of closures to be used

      int fftw_planner; // = FFT_ESTIMATE
      // .true.  - use aligned memory and to enable SIMD;
      // .false. - don't use aligned memory
      bool fft_aligned = true;
      // Transpose site number and spatial data locally before and after FFT.
      bool fftw_localtrans = true; 

      // Output verbosity.  Useful for debugging.
      // 0 - no ouput
      // 1 - memory allocation and steps for convergence
      // 2 - 1 + convergence progress
      int verbose = 0;

      // A failure has occurred
      bool rism_failure = false;

      // Tracking convergence progress
      bool converged;

      // This is a bit ugly and there may be a better solution.  We
      // need to keep track of the number of solutions for both charged
      // and un-charged solutes.  When we change between the solutes we
      // set the nsolutions pointer to the appropriate variable.
      // However, the 'target' attribute is not allowed in type
      // definitions so these variables have to be pointers and we have
      // to allocate memory for them.
      // int* nsolution;
      int* nsolutionChg;
      int* nsolutionNoChg;

      // Center the solute in the solvation box.
      // 0 - off
      // 1 - center of mass
      // 2 - center of geometry
      // 3 - center of mass shifted to the nearest grid point
      // 4 - center of geometry shifted to the nearest grid point
      // For negative numbers the centering translation is only
      // calculated the for the first solution and used for subsequent
      // calculations.  This allows the solute to drift in the box.
      int centering;

      // Number of past direct correlation function time step saves.
      int ncuvsteps;

      // Buffer distance to the edge of the box for the solvent. [A]
      double buffer;
      double effect_buffer[3];

      // Fixed box size for 3D-RISM.
      double fixedBoxDimensionsR[3];
      // Number of Cartesian grid points in each dimension for a fixed box size.
      int fixedNumGridPoints[3];

      // Variable box size.
      bool varbox;

      // long-range asymptotics k-space cut off tolerance.  Only grid
      // points that have an approximate value greater than this will be computed.
      // -1 - cutoff selected from calculation tolerance
      // 0  - no cutoff applied
      // >0 - used as the tolerance.  Should be a small value.  E.g., 1e-7      
      double asympKSpaceTolerance;

      // Lennard-Jones potential tolerance
      // -1 - cutoff selected from calculation tolerance
      // 0  - no cutoff applied
      // >0 - used as the tolerance.  Should be a small value.  E.g., 1e-7      
      double ljTolerance;

      // Number of vectors used for MDIIS (consequently, the number of
      // copies of CUV we need to keep for MDIIS).      
      int NVec;
      // MDIIS implementation to use.
      int mdiis_method;

      // 'Step size' for MDIIS.
      double deloz;
      // Restart threshold factor. Ratio of the current residual to the
      // minimum residual in the basis that causes a restart.
      double mdiis_restart;

      // MPI Support !! I do not think we will need these ones...
      int mpirank = 0;
      int mpicomm = 0;
      int mpisize = 1;

      // !! LARGE ARRAYS !!
      // 
      //  all arrays are declared as pointers to ensure we can reallocate them as necessary
      // 
      //
      //  xvva       :: solvent chi interpolated for our grid size
      //  guv        :: solvent distribution function
      //  huv        :: guv - 1
      //  cuv        :: solvent direct correlation function and points to the
      //               current active solution in cuvWRK
      //  cuvres     :: residual value for cuv calculation and points to the
      //               current active solution in cuvresWRK.
      //  cuvWRK     :: Working Cuv memory.  Holds Cuv from previous iterations.
      //  cuvresWRK  :: Working Cuvres memory.  Holds Cuvres from previous iterations.
      //  oldcuv     :: previous solutions of cuv. Points to oldcuvChg
      //               or oldcuvNoChg depending on the charge state of the
      //               calculation.
      //  oldcuvChg  :: previous solutions for the standard charged system
      //  oldcuvNoChg :: previous solutions for the chargeless system.  This is only allocated
      //                 if _unsetCharges() is called
      //  electronMap :: smeared solvent electron density map.
      double *oldcuv;

      array_class<GPUtype> xvva;
      array_class<GPUtype> cuv;

      array_class<GPUtype> cuvWRK;
      array_class<GPUtype> oldcuvChg;
      array_class<GPUtype> oldcuvNoChg;

      array_class<GPUtype> cuvres;
      array_class<GPUtype> cuvresWRK;

      array_class<GPUtype> guv;
      array_class<GPUtype> huv;
      array_class<double> electronmap;

      // double *guv, *huv;
      // double *electronmap; // do we need this one?

      // Temperature derivative memory
      // cuvk        :: k-space Cuv solution from 3D-RISM
      //               solution. NOTE: we should consider using Huv or
      //               Guv memory instead.  However, it has to be
      //               checked first that it is not used for any thermodynamics calculations
      // xvva_dT     :: solvent dT chi interpolated for our grid size
      // guv_dT      :: solvent dT distribution function
      // huv_dT      :: guvdT
      // cuv_dT      :: solvent dT direct correlation function and points to the
      //              current active solution in cuvWRK
      // cuvres_dT   :: residual value for dT cuv calculation and points to the
      //              current active solution in cuvresWRK.
      // cuvWRK_dT   :: Working dT Cuv memory.  Holds Cuv from previous iterations.
      // cuvresWRK_dT:: Working dT Cuvres memory.  Holds Cuvres from previous iterations.
      double *xvva_dT, *cuv_dT, *cuvres_dT, *cuvWRK_dT, *cuvresWRK_dT;
      double *cuvk, *huv_dT;
      array_class<double> guv_dT;

      // Lengths and interior angles of the unit cell. For aperiodic
      // systems, the interior angles are always 90 degrees.      
      double unitCellDimensions[6];

      // The abbreviated label of the periodic potential function used
      // for periodic calculations. See rism3d_potential for valid values.  
      string periodicPotential;

      ///////////////////////////////////////////////////////////////////////////////
      ///////////////////////////// FUNCTION ATTRIBUTES /////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////

      // void set_solvent(solvent *arg);   
      // void set_solute(solute *arg);   

      // Using the current orientation of the solute, define the minimum box size
      // and resize all associated grids.
      //
      // Non-periodic : Either a variable or fixed box size maybe specified
      //     (this%varbox). If fixed size, then the number of gridpoints provided by
      //     the user are directly applied but check to ensure they meet the minimum
      //     requirements of the MPI run, if specified. If variable size the buffer
      //     is either (a) taken from the user or (b), if buffer==0, determined from
      //     the ljTolerance cutoff. Grid spacing is taken from the user.
      //
      // Periodic : the unit cell is used to determine the number of grid points.
      //
      // @param[in] this rism3d object
      // MY COMMENT: IN .F90 PRIVATE LIST
      void resizeBox(); 

      // Get the minimum lengths of the solvent box using the buffer size and/or
      // ljTolerance.  If buffer >0, use that directly. If buffer = 0, use LJ
      // cutoffs in potential to ensure that the cutoffs are completely in the
      // solvent box. buffer <0 is an error.
      //
      // OUT:
      //     (_REAL_(3)) x, y, z dimensions
      // MY COMMENT: IN .F90 PRIVATE LIST
      void getMinBoxLength(double boxlen[3]); // return 1D array of size = 3

      // Calculate the effective minimum buffer distance between the
      // solute and the solvent box edge.
      // @param[in] this rism3d object.
      void effective_buffer(); // return 1D array of size = 3

      // Using the current box size and resize all associated grids and
      // variables.
      // @param[in,out] this rism3d object.
      // @param[in] ngr Number of grid points along each axis.
      // @param[in] grdspc Grid spacing along each axis.  
      // MY COMMENT: IN .F90 PRIVATE LIST    
      void reallocateBox(int ngr[3], GPUtype spacing[3]);

      // Using the current box size and resize all the dT associated grids
      // and variables.
      // IN:
      //   this :: rism3d object
      // MY COMMENT: IN .F90 PRIVATE LIST
      void reallocateBox_dT();

      // Prints the maximum amount of memory allocated at any one time so
      // far in the run.
      void rism3d_max_memory();

      // Interpolate the solvent-solvent susceptibility, solved on the
      // 1D-RISM grid, to the 3D-RISM grid.
      // @param[in,out] this rism3d object.
      // @param[in] xvv 1D-RISM Xvv or Xvv_dT data.
      // @param[out] xvva Interpolated result.
      // MY COMMENT: IN .F90 PRIVATE LIST
      void interpolateSolventSusceptibility();

      // Center the solute in the solvent box.
      // @param[in,out] this rism3d object.
      // @param[in] origin_o (optional) position to center on
      // @param[out] translation_o (optional) the amount that the solute is moved by
      // MY COMMENT: IN .F90 PRIVATE LIST
      // void centerSolute(GPUtype translation_o[3], GPUtype origin_o[3], 
      //                   bool calc_translation = true, bool default_orig = true);
      void centerSolute(GPUtype *translation_o = NULL, GPUtype *origin_o = NULL);

      //////////////// 
      // Subroutines to find the iterative 3D-RISM solution.
      ///////////////

      // Main driver for the 3D-RISM solver.
      // Makes an initial guess of the direct correlation function and
      // then solve the RISM and closure relations until either the
      // solution converges or the maximum of steps is reached.
      // @param[in,out] this rism3d object.
      // @param[in] ksave Save itermediate results every ksave interations
      //  (0 means no saves).
      // @param[in] kshow Print parameter for relaxation steps every kshow
      //  iteration (0 means no saves).
      // @param[in] maxSteps Maximum number of rism relaxation steps.
      // @param[in] tolerance Tolerance in.
      // MY COMMENT: IN .F90 PRIVATE LIST
      void solve3DRISM(int maxSteps, double tolerance);

      // One relaxation step for the UV-RISM equation with the HNC closure,
      // Guv(r) = exp(-Uuv(r) + Tuv(r) - DelHv0) + DelHv0
      // Cuv(r) = Guv(r) - 1 - Tvv(r)
      // Huv(k) = Cuv(k) * (Wvv(k) + Density * Hvv(k))
      // TuvRes(r) = Huv(r) - Guv(r) - 1
      // @param[in,out] this A rism3d object.
      // @param[in,out] residual ???
      // @param[in,out] converged Returns true if the solution has converged.
      // @param[in] tolerance Target residual tolerance for convergence.
      // MY COMMENT: IN .F90 PRIVATE LIST
      void single3DRISMsolution();

      void convert2nr();
      void convert2nr_cu();

      void convert2fftw();
      void convert2fftw_cu();

      void set_padding2zero(GPUtype* data);
      void get_residue();
      void get_h_k();

      // Set DCF long range values for charged particles
      void set_dcf_longrange_cu(GPUtype* cuv, GPUtype* solvent_charge, GPUtype* dcfLongRangeAsympR);

      // Function to subtract DCF long range in the real space
      void subtract_dcf_longrange_cu(GPUtype* guv, GPUtype* cuv,  GPUtype* solvent_charge, GPUtype* dcfLongRangeAsympR);

      // Function to add DCF long range in the reciprocal space
      void add_dcf_longrange_cu(GPUtype* guv, GPUtype* solvent_charge, GPUtype* dcfLongRangeAsympK);

      // Function to subtract TCF long range in reciprocal space
      void subtract_tcf_longrange_cu(GPUtype* huv, GPUtype* solvent_charge_sp, GPUtype* tcfLongRangeAsympK);

      // Function to add TCF long range in real space
      void add_tcf_longrange_cu(GPUtype* huv, GPUtype* solvent_charge_sp, GPUtype* tcfLongRangeAsympR);

      // Main driver for the 3D-RISM ΔT solver.
      // IN:
      //   this :: rism3d object
      //   kshow  :: print parameter for relaxation steps every kshow iterations
      //   maxSteps :: maximum number of rism relaxation steps
      //   mdiis_method :: MDIIS implementation to use
      // MY COMMENT: IN .F90 PRIVATE LIST
      void solve3DRISM_dT();

      // One relaxation step for the UV-RISM equation with the HNC
      // closure,
      // Guv(r) = exp(-this%potential%uuv(r) + Tuv(r) - DelHv0) + DelHv0
      // Cuv(r) = Guv(r) - 1 - Tvv(r)
      // Huv(k) = Cuv(k) * (Wvv(k) + Density * Hvv(k))
      // TuvRes(r) = Huv(r) - Guv(r) - 1
      // IN:
      //  this :: rism3d object
      //  residual ::
      //  converged ::
      //  mdiis  :: MDIIS object to accelerate convergence
      // MY COMMENT: IN .F90 PRIVATE LIST
      void single3DRISMsolution_dT();

      ///////////////////////// PROPAGATE PREVIOUS SOLUTIONS

      // Calculates a new initial guess for CUV based on the final solutions
      // from previous timesteps.  The maximum number of previous time
      // steps to use is provided by the user in ncuvsteps.  However, if
      // there are not enough previous timesteps only nsolution previous
      // timesteps will be used.
      //
      // See section 2.3.1 and eqs. 8-13 of doi:10.1021/ct900460m for
      // details.
      // @param[in] this rism3d object.
      // MY COMMENT: IN .F90 PRIVATE LIST
      void guessDCF();

      // Updates the values in the this%oldcuv queue.  The oldest value (the
      // ncuvstep index) is pushed out, the remainder of the data is
      // shifted and the newest solution is placed in the first index.
      // @param[in] this rism3d object.
      // MY COMMENT: IN .F90 PRIVATE LIST
      void updateDCFguessHistory();

      // Convert the user supplied a,b,a1,b1 Universal Correction
      // coefficients to a0,b0,a1,b1 coefficients used in
      // rism3d_closure_c.
      // IN:
      //   this :: rism3d object with computed solution
      //   coeff: coefficients for correction.  For the original correction
      //          (a, b) = coeff(1:2).
      //          Extra coefficients are a1 and b1 from Johnson et al. 2016.
      // OUT:
      //   array of length 4 containing a0,b0,a1,b1
      double* UC_temperature_coeff(); // return 1D array of size = 4

      // TO DO: add some description here...
      // returns wvv(:,:,:)
      double* rism3d_intramolecular();

      // Sets the parameters for a variable solvation box.
      // IN:
      //  this :: rism3d object
      //  buffer :: shortest distance between solute and solvent box boundary
      //  grdspc :: linear grid spacing for the solvent box in each dimension
      // MY COMMENT: IN .F90 PUBLIC LIST
      void setbox_variable(rism3d_grid *grid, double o_buffer, double* o_grdspc);

      // Sets the parameters for a fixed solvation box.
      // IN:
      //  this :: rism3d object
      //  boxlen :: solvent box size in each dimension in Angstroms
      //  ng3    :: number of grid points in each dimension      
      // MY COMMENT: IN .F90 PUBLIC LIST
      void setbox_fixed(double* o_boxlen, int* o_ng3);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////            //////////////////////////////////////////////
      //////////////////////////////////////////////   PUBLIC   //////////////////////////////////////////////
      //////////////////////////////////////////////            //////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////

    public:
      // solute solu;
      // solvent solv;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////// FUNCTION ATTRIBUTES /////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////

      // void get_atomnames(char** names, int numAtomType);
      // void set_atomnames(char** names, int numAtomType, int str_len);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      //// BEGIN... FUNCTIONS CALLED FROM AMBER_RISM_INTERFACE.F90
      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      // Defining constructor and destructor:
      // To do: add more information here...
      // MY COMMENT: IN .F90 PUBLIC LIST
      // I THINK WE CAN REMOCE THE DEFINED TYPES ARGUMENTS, 
      // SINCE WE CAN CREATE THE OBJECT AS A CLASS ATTRIBUTE 
      // AND ACCESS IT (AS DONE IN TUTORIAL EXAMPLE)
      // rism3d(solute_cpp *solu, solvent_cpp *solv, int centering, int ncuvsteps,
      //   char** closure, double cut, int mdiis_nvec, double mdiis_del, int mdiis_method,
      //   double mdiis_restart, bool treeDCF, bool treeTCF, bool treeCoulomb, 
      //   double treeDCFMAC, double treeTCFMAC, double treeCoulombMAC,
      //   int treeDCFOrder, int treeTCFOrder, int treeCoulombOrder, 
      //   int treeDCFN0, int treeTCFN0, int treeCoulombN0,
      //   double asympKSpaceTolerance, double ljTolerance, double chargeSmear,
      //   double o_buffer = 0, double* o_grdspc = pd_grdspc, int o_mpicomm = 0,
      //   char* o_periodic = pd_periodic, double* o_unitCellDimensions = pd_unitCellDimensions, 
      //   double o_biasPotential = 0); // constructor = rism3d_new

      rism3d(solute_cpp *solu_f, solvent_cpp *solv_f, int centering, int ncuvsteps,
        char** closure, double cut, int mdiis_nvec, double mdiis_del, int mdiis_method,
        double mdiis_restart, bool treeDCF, bool treeTCF, bool treeCoulomb, 
        double treeDCFMAC, double treeTCFMAC, double treeCoulombMAC,
        int treeDCFOrder, int treeTCFOrder, int treeCoulombOrder, 
        int treeDCFN0, int treeTCFN0, int treeCoulombN0,
        double asympKSpaceTolerance, double ljTolerance, double chargeSmear,
        double o_buffer, double* o_grdspc, int o_mpicomm,
        char* o_periodic, double* o_unitCellDimensions, 
        double o_biasPotential); // constructor = rism3d_new

      // rism3d(solute_cpp *solu, solvent_cpp *solv, int centering, int ncuvsteps,
      //   char** closure, double cut, int mdiis_nvec, double mdiis_del, int mdiis_method,
      //   double mdiis_restart, bool treeDCF, bool treeTCF, bool treeCoulomb, 
      //   double treeDCFMAC, double treeTCFMAC, double treeCoulombMAC,
      //   int treeDCFOrder, int treeTCFOrder, int treeCoulombOrder, 
      //   int treeDCFN0, int treeTCFN0, int treeCoulombN0,
      //   double asympKSpaceTolerance, double ljTolerance, double chargeSmear,
      //   double o_boxlen = 0, int* o_ng3 = pd_ng3, int o_mpicomm = 0, 
      //   char* o_periodic = pd_periodic, double* o_unitCellDimensions = pd_unitCellDimensions, 
      //   double o_biasPotential = 0); // constructor2 = rism3d_new

      rism3d(solute_cpp *solu_f, solvent_cpp *solv_f, int centering, int ncuvsteps,
        char** closure, double cut, int mdiis_nvec, double mdiis_del, int mdiis_method,
        double mdiis_restart, bool treeDCF, bool treeTCF, bool treeCoulomb, 
        double treeDCFMAC, double treeTCFMAC, double treeCoulombMAC,
        int treeDCFOrder, int treeTCFOrder, int treeCoulombOrder, 
        int treeDCFN0, int treeTCFN0, int treeCoulombN0,
        double asympKSpaceTolerance, double ljTolerance, double chargeSmear,
        double* o_boxlen, int* o_ng3, int o_mpicomm, 
        char* o_periodic, double* o_unitCellDimensions, 
        double o_biasPotential); // constructor2 = rism3d_new
        // add or remove?: rism3d_solute solute, rism3d_solvent solvent, 

      ~rism3d(); // destructor = rism3d_destroy

      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      //// BEGIN... ATTRIBUTES CALLED FROM AMBER_RISM_INTERFACE.F90
      ////////////////////////////////////////////////////////////////////////////////////////////////////////

      /////////////////// From grid class
      double* get_spacing(int* len);
      double* get_voxelvectorsr(int* len1, int* len2);
      double* get_unitcellangles(int* len);
      double* get_unitcellvectorsk(int* len1, int* len2);
      double* get_unitcellvectorsr(int* len1, int* len2);
      double* get_boxlength(int* len);
      int get_totallocalpointsr();
      double get_voxelvolume();
      int* get_localdimsr(int* len);
      int get_localdimsr(int indx);
      int* get_globaldimsr(int* len);
      int get_globaldimsr(int indx);

      ////////////////// From solvent struct
      double get_xikt_dt();
      double get_temperature();
      void set_atomname(char** names);
      void get_atomname(string** names, int* name_len);
      int get_numatomtypes();
      double* get_charge(int* len);
      double get_charge(int indx);
      double* get_density(int* len);
      double get_density(int indx);
      bool get_ionic();

      ////////////////// From solute struct
      double* get_centerofmass(int* len);
      double* get_translation(int* len);
      bool get_charged();
      int get_numatoms();
      double* get_mass(int* dim1);

      ////////////////// From rism3d class
      double* get_guv(int *dim1, int *dim2, int arg1, int arg2);
      double* get_guv(int *dim1, int *dim2, int arg);
      
      double* get_huv(int *dim1, int *dim2, int arg1, int arg2);
      double* get_huv(int *dim1, int *dim2, int arg);
      
      double get_cuv(int arg1, int arg2, int arg3, int arg4);
      double* get_cuv(int *dim1, int *dim2, int *dim3, int *dim4, int arg);
      
      // double* get_guv_dt(int *dim1, int *dim2); // add to yaml add default index value here to get value or array
      double* get_guv_dt(int *dim, int arg);
      double* get_guv_dt(int *dim, int arg1, int arg2);
      
      double* get_electronmap(int *dim1, int *dim2, int *dim3);
      int get_nsolution();
      void set_nsolution(int nsol);
      bool get_periodic(); 
      bool get_failure();

      ////////////////// From potential class
      double* get_tcflongrangeasympr(int *dim1);
      double* get_dcflongrangeasympr(int *dim1);
      double get_dcflongrangeasympr(int indx);
      double* get_uuv(int *dim1, int *dim2, int *dim3, int *dim4, int indx);

      // This is a bit ugly and there may be a better solution.  We
      // need to keep track of the number of solutions for both charged
      // and un-charged solutes.  When we change between the solutes we
      // set the nsolutions pointer to the appropriate variable.
      // However, the 'target' attribute is not allowed in type
      // definitions so these variables have to be pointers and we have
      // to allocate memory for them.
      int* nsolution;
      // int* nsolutionChg;
      // int* nsolutionNoChg;

      // If true, a periodic 3D-RISM calculation is performed. This
      // primarily differs from infinite dilution 3D-RISM by using
      // Ewald sum potential in place of Coulombic potential and
      // invoking the minimum image convention while calculating both
      // the Ewald sum and Lennard-Jones potentials.
      bool periodic = false;

      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      //// END... ATTRIBUTES CALLED FROM AMBER_RISM_INTERFACE.F90
      ////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Sets verbosity of output.
      // IN:
      //   this :: rism3d object
      //   verbosity :: 0 - no output
      //                1 - memory allocation and steps for convergence
      //                2 - 1 + convergence progress      
      // MY COMMENT: IN .F90 PUBLIC LIST
      void setverbosity(int verbosity);

      // Set parent for this timer
      // IN:
      //   this : rism3d object
      //   parent : parent timer object      
      void settimerparent(timer_cpp *timer);

      // Sets solute coordinates.
      // IN:
      //   this :: rism3d object
      //   ratu :: coordinates   
      void setcoord(double* solutePositions);

      // Sets the closure list and sets the current closure to the first one
      // in the list.  When there is no previous solution to work from, the
      // solver will use each closure in the list in turn. By choosing the
      // list to increase in order, it makes it possible to converge
      // otherwise difficult closures. Only the last closure is used for
      // thermodynamic output.
      // IN:
      //   this :: rism3d object
      //   closure :: array of closure types (see closure enumeration).
      void set_closurelist(char **names, int *N);

      // Calculates the full 3D-RISM solvent distribution.  This is required to
      // calculate thermodynamic quantities.
      // @param[in,out] this rism3d object.
      // @param[in,out] ksave Save intermediate results every ksave
      //            interations (0 means no saves).
      // @param[in] kshow Print parameter for relaxation steps every kshow
      //            iteration (0 means no print).
      // @param[in] maxSteps Maximum number of rism relaxation steps.
      // @param[in] tolerance Convergence tolerances. There should be one
      //          tolerance per closure in the closure list. 
      // MY COMMENT: IN .F90 PUBLIC LIST     
      bool calculatesolution(int ksave, int kshow, int maxSteps, bool failure, double* tolerance, int tol_size);

      // Calculates the forces on the solute contributed by the solvent according
      // to 3D-RISM.  Just a wrapper for rism3d_closure_force().
      // IN:
      //   this :: rism3d object with computed solution
      //   ff   :: 3D-RISM forces
      // MY COMMENT: IN .F90 PUBLIC LIST
      void force(double* ff);

      // Returns a 3D map of the excess chemical potential in kT.
      // Integrating this map gives the total excess chemical potential as
      // returned from rism3d_excessChemicalPotential_tot().  Memory is allocated into a
      // pointer and must be freed by the calling function. For MPI, only
      // the grid points local to this process are allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the excess chemical potential contributions
      // SIDEEFFECTS:
      //   memory is allocated for the grid
      double* excesschemicalpotential_tot_map(int* dim1, int* dim2, int* dim3); //

      // Returns a 3D map of the excess chemical potential in kT.
      // Integrating this map gives the site excess chemical potential
      // contribution as returned from rism3d_excessChemicalPotential_tot().  Memory is
      // allocated into a pointer and must be freed by the calling
      // function. For MPI, only the grid points local to this process are
      // allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the excess chemical potential contributions (nx, ny, nz, nsite)
      // SIDEEFFECTS:
      //   memory is allocated for the grid
      double* excesschemicalpotential_site_map(int* dim1, int* dim2, int* dim3, int* dim4); //

      // Returns a 3D map of the excess chemical potential in kT from the GF
      // approximation.  Integrating this map gives the total excess
      // chemical potential as returned from rism3d_excessChemicalPotential_tot().  Memory is
      // allocated into a pointer and must be freed by the calling
      // function. For MPI, only the grid points local to this process are
      // allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the excess chemical potential contributions
      // SIDEEFFECTS:
      //   memory is allocated for the grid
      double* excesschemicalpotentialgf_tot_map(int* dim1, int* dim2, int* dim3); //

      // Returns a 3D map of the excess chemical potential in kT using the
      // GF approximation.  Integrating this map gives the site excess
      // chemical potential contribution as returned from
      // rism3d_excessChemicalPotential_tot().  Memory is allocated into a pointer and must
      // be freed by the calling function. For MPI, only the grid points
      // local to this process are allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the excess chemical potential contributions (nx, ny, nz, nsite)
      // SIDEEFFECTS:
      //   memory is allocated for the grid
      double* excesschemicalpotentialgf_site_map(int* dim1, int* dim2, int* dim3, int* dim4); //

      // Returns a 3D map of the excess chemical potential in kT from the
      // PCPLUS.  Integrating this map gives the total excess chemical
      // potential as returned from rism3d_excessChemicalPotentialPCPLUS_tot().
      // Memory is allocated into a pointer and must be freed by the
      // calling function. For MPI, only the grid points local to this
      // process are allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the excess chemical potential contributions
      // SIDEEFFECTS:
      //   memory is allocated for the grid
      double* excesschemicalpotentialpcplus_tot_map(int* dim1, int* dim2, int* dim3); //

      // Returns a 3D map of the excess chemical potential in kT from the
      // Universal Correction correction.  Integrating this map gives the
      // total excess chemical potential as returned from
      // rism3d_excessChemicalPotential_tot().
      //
      // Memory is allocated into a pointer and must be freed by the
      // calling function. For MPI, only the grid points local to this
      // process are allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      //   coeff: coefficients for correction.  For the original correction
      //          (a, b) = coeff(1:2).
      //          Extra coefficients are a1 and b1 from Johnson et al. 2016.
      // OUT:
      //   a 3D-grid of the excess chemical potential contributions
      // SIDEEFFECTS:
      //   memory is allocated for the grid
      double* excesschemicalpotentialuc_tot_map(int* dim1, int* dim2, int* dim3, double* coeff);

      // Calculate the total excess chemical potential of solvation
      // IN:
      //   this :: rism3d object with computed solution
      //   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
      // OUT:
      //    total excess chemical potential of solvation   
      // MY COMMENT: IN .F90 PUBLIC LIST   
      double excesschemicalpotential_tot(bool o_lr = true);      

      // Calculate the excess chemical potential for each solvent species
      // IN:
      //   this :: rism3d object with computed solution
      //   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
      // OUT:
      //    excess chemical potential of solvation for each solvent species      
      // MY COMMENT: IN .F90 PUBLIC LIST
      double* excesschemicalpotential(int* len, bool o_lr = true);

      // Calculate the excess chemical potential of solvation for each solvent species
      // with the Gaussian fluctuation correction
      // IN:
      //   this :: rism3d object with computed solution
      //   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
      // OUT:
      //    Gaussian fluctuation excess chemical potential of solvation for each
      //    solvent species      
      // MY COMMENT: IN .F90 PUBLIC LIST
      double* excesschemicalpotentialgf(int* len, bool o_lr = true);

      // Calculate the total excess chemical potential of solvation with the
      // PC+/3D-RISM Correction
      // IN:
      //    this :: rism3d object with computed solution
      //    o_lr   :: (optional) (default=.true.) Apply asymptotic long range correction
      // OUT:
      //     total Universal Correction excess chemical potential of solvation      
      double excesschemicalpotentialpcplus(bool o_lr = true); //

      // Calculate the total excess chemical potential of solvation with the
      // Palmer et al. Universal Correction with optional temperature dependence.
      //
      // Uses the total molecular density of the solvent.  For pure water,
      // this gives the original Palmer correction.  There is no definition
      // for mixed solvents but this seems reasonable.
      //
      // IN:
      //   this :: rism3d object with computed solution
      //   coeff: coefficients for correction.  For the original correction
      //      (a, b) = coeff(1:2). Extra coefficients are a1, b1 from
      //      Johson et al. 2016
      //   o_lr :: (optional) (default = .true.) Apply asymptotic long
      //       range correction
      // OUT:
      //    total Universal Correction excess chemical potential of solvation
      double excesschemicalpotentialuc(double* coeff, bool o_lr = true);

      // Returns a 3D map of the solvent-solute potential energy in kT.
      // Integrating this map gives the site solvent-solute potential energy
      // contribution as returned from rism3d_closure_solvPotEne(, .false.).
      // Memory is allocated into a pointer and must be freed by the calling
      // function. For MPI, only the grid points local to this process are
      // allocated and calculated.
      // IN:
      //    this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the solvent-solute energy contributions (nx, ny, nz, nsite)
      // SIDEEFFECTS:
      //   memory is allocated for the grid      
      double* solventpotene_tot_map(int* dim1, int* dim2, int* dim3); // returns a 3D array

      // Returns a 3D map of the solvent-solute potential energy in kT.
      // Integrating this map gives the site solvent-solute potential energy
      // contribution as returned from rism3d_closure_solvPotEne(, .false.).
      // Memory is allocated into a pointer and must be freed by the calling
      // function. For MPI, only the grid points local to this process are
      // allocated and calculated.
      // IN:
      //    this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the solvent-solute energy contributions (nx, ny, nz, nsite)
      // SIDEEFFECTS:
      //   memory is allocated for the grid      
      double* solventpotene_site_map(int* dim1, int* dim2, int* dim3, int* dim4);

      // Calculate the solvation interaction energy: de = density sum g*u for
      // each solvent site.  I.e., the direct intection potential energy of
      // solute and solvent and not the total solvation energy (see solvationEnergy).
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //    the contribution of each solvent site to the total solvation interaction
      //    energy [kT]      
      double* solventpotene(int* len); // returns an 1D array with size = numAtomTypes

      // Calculating the partial molar volume of solute.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   partial molar volume
      double partialmolarvolume();

      // Calculating excess number of each solvent type associated with
      // the solute.
      // @param[in,out] this rism3d object with computed solution.
      // @param[in] o_lr (optional) (default = .true.)
      //                 Apply asymptotic long range correction.
      // @return Excess number of each solvent type associated with the solute.
      double* excessparticles(int* len, bool o_lr = true); // return 1D array of size = numAtomType
      // std::vector<double> excessparticles(bool o_lr = true);

      // Calculate the Kirkwood-Buff integral for the solute. This is the
      // all space integral of huv.
      //
      // J. G. Kirkwood; F. P. Buff. J. Chem. Phys. 1951, 19, 774-777
      // IN:
      //    this :: rism3d object with computed solution
      //    o_lr :: (optional) (default = .true.) Apply asymptotic long range
      //            correction
      // OUT:
      //    Kirkwood-Buff integeral for each solvent site
      double* kirkwoodbuff(int* len, bool o_lr = true); // return 1D array of size = numAtomType
      // std::vector<double> kirkwoodbuff(bool o_lr = true);

      // Calculates the direct correlation function integral for the solute. This is the
      // all space integral of cuv.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //    DCF integeral for each solvent site
      double* dcfintegral(int* len); // return 1D array of size = numAtomType

      // std::vector<double> dcfintegral(); // return 1D array of size = numAtomType

      // Check if we can perform a temperature derivative calculation
      // (i.e. all the necessary information is available and the closure
      // supports it).
      // IN:
      //   this : rism3d object
      // OUT:
      //   .true. if we can, .false. if we can't
      bool cancalc_dt();

      // Calculates the full 3D-RISM solvent distribution.  This is required to
      // calculate solvation energies and entropies.
      // IN:
      //   this :: rism3d object
      //   ksave  :: save itermediate results every ksave interations (0 means no saves)
      //   kshow  :: print parameter for relaxation steps every kshow iteration (0 means no saves)
      //   maxSteps :: maximum number of rism relaxation steps
      //   tolerance    :: convergence tolerance
      bool calculatesolution_dt(int kshow, int maxSteps, bool failure, double tolerance);

      // Calculate the solvation energy dE = dmu + dTS which includes both
      // solute-solvent interaction energy and the change in solvent
      // self-energy due to rearranging around the solvent.
      // IN:
      //   this :: rism3d object with computed solution
      //   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
      // OUT:
      //    solvation energy
      double* solvationenergy(int* len, bool o_lr = true); // return a 1D array with size = numAtomTypes

      // Calculate the solvation energy dE = dmu + dTS the Gaussian
      // fluctuation correction, which includes both solute-solvent
      // interaction energy and the change in solvent self-energy due to
      // rearranging around the solvent.
      // IN:
      //   this :: rism3d object with computed solution
      //   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
      // OUT:
      //    solvation energy
      double* solvationenergygf(int* len, bool o_lr = true); // return a 1D array with size = numAtomTypes

      // Calculating temperature derivative (T * d/dT) of the excess number of each
      // solvent type associated with the solute.
      // IN:
      //   this :: rism3d object with computed solution
      //   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
      // OUT:
      //    excess number temperature derivative of each solvent type
      //    associated with the solute
      double* excessparticles_dt(int* len, bool o_lr = true); // return 1D array of size = numAtomType

      // Calculate the Kirkwood-Buff temperature derivative (T*d/dT) integral for the
      // solute. This is the all space integral of huv_dT.
      //
      // J. G. Kirkwood; F. P. Buff. J. Chem. Phys. 1951, 19, 774-777
      // IN:
      //    this :: rism3d object with computed solution
      //    o_lr :: (optional) (default = .true.) Apply asymptotic long range
      //            correction
      // OUT:
      //    Kirkwood-Buff integeral temperature derivative for each solvent site
      double* kirkwoodbuff_dt(int* len, bool o_lr = true); // return 1D array of size = numAtomType

      // Calculates the direct correlation function integral temperature
      // derivative (T*d/dT) for the solute. This is the all space integral of cuv.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //    DCF integeral temperature derivative for each solvent site
      double* dcfintegral_dt(int* len); // return 1D array of size = numAtomType

      // Calculate the total solvation energy with the PC+/3D-RISM
      // Correction.
      // 
      // IN:
      //    this :: rism3d object with computed solution
      //    o_lr   :: (optional) (default=.true.) Apply asymptotic long range correction
      // OUT:
      //     total PC+/3D-RISM Correction excess chemical potential of
      //     solvation
      double solvationenergypcplus(bool o_lr = true);

      // Calculate the total solvation energy with the Palmer et
      // al. Universal Correction with optional temperature dependence.
      // 
      // Uses the total molecular density of the solvent.  For pure water,
      // this gives the original Palmer correction.  There is no
      // definition for mixed solvents but this seems reasonable.
      //
      // @param[in] this rism3d object with computed solution.
      // @param[in] coeff coefficients for correction.  For the original
      //     correction (a, b) = coeff(1:2). Extra coefficients are a1, b1
      //     from Johnson et al. 2016
      // @param[in] o_lr (optional) (default = .true.)
      //                 Apply asymptotic long range correction.
      // @return Total Universal Correction excess chemical potential of solvation.
      double solvationenergyuc(double* coeff, bool o_lr = true);

      // Returns a 3D map of the solvation energy in kT.
      // Integrating this map gives the total solvation energy as
      // returned from rism3d_excessChemicalPotential_tot(, .false.).  Memory is allocated into a
      // pointer and must be freed by the calling function. For MPI, only
      // the grid points local to this process are allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the solvation energy contributions
      // SIDEEFFECTS:
      //   memory is allocated for the grid     
      double* solvationenergy_tot_map(int* dim1, int* dim2, int* dim3);

      // Returns a 3D map of the solvation energy in kT.
      // Integrating this map gives the site solvation energy contribution as
      // returned from rism3d_excessChemicalPotential_tot(, .false.).  Memory is allocated into a
      // pointer and must be freed by the calling function. For MPI, only
      // the grid points local to this process are allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the solvation energy contributions (nx, ny, nz, nsite)
      // SIDEEFFECTS:
      //   memory is allocated for the grid      
      double* solvationenergy_site_map(int* dim1, int* dim2, int* dim3, int* dim4);

      // Returns a 3D map of the solvation energy in kT from the GF
      // approximation.  Integrating this map gives the total solvation
      // energy as returned from rism3d_excessChemicalPotential_tot(, .false.).  Memory is
      // allocated into a pointer and must be freed by the calling
      // function. For MPI, only the grid points local to this process are
      // allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the solvation energy contributions
      // SIDEEFFECTS:
      //   memory is allocated for the grid      
      double* solvationenergygf_tot_map(int* dim1, int* dim2, int* dim3);

      // Returns a 3D map of the solvation energy in kT from the GF approximation.
      // Integrating this map gives the site solvation energy contribution as
      // returned from rism3d_excessChemicalPotential_tot(, .false.).  Memory is allocated into a
      // pointer and must be freed by the calling function. For MPI, only
      // the grid points local to this process are allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the solvation energy contributions (nx, ny, nz, nsite)
      // SIDEEFFECTS:
      //   memory is allocated for the grid      
      double* solvationenergygf_site_map(int* dim1, int* dim2, int* dim3, int* dim4);

      // Returns a 3D map of the solvation energy in kT from the PCPLUS.
      // Integrating this map gives the total solvation energy as returned
      // from rism3d_excessChemicalPotential_tot(, .false.).  Memory is
      // allocated into a pointer and must be freed by the calling
      // function. For MPI, only the grid points local to this process are
      // allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the solvation energy contributions
      // SIDEEFFECTS:
      //   memory is allocated for the grid      
      double* solvationenergypcplus_tot_map(int* dim1, int* dim2, int* dim3);

      // Returns a 3D map of the solvation energy in kT from the Universal
      // Correction with optional temperature dependence.
      //
      // Integrating this map gives the total solvation energy as returned
      // from rism3d_excessChemicalPotential_tot(, .false.).  Memory is
      // allocated into a pointer and must be freed by the calling
      // function. For MPI, only the grid points local to this process are
      // allocated and calculated.
      // IN:
      //   this :: rism3d object with computed solution
      //   coeff: coefficients for correction.  For the original
      //       correction (a, b) = coeff(1:2). Extra coefficients are a1,
      //       b1 from Johnson et al. 2016
      // OUT:
      //  a 3D-grid of the solvation energy contributions
      // SIDEEFFECTS:
      //   memory is allocated for the grid
      double* solvationenergyuc_tot_map(int* dim1, int* dim2, int* dim3, double* coeff);

      // Calculating the partial molar volume temperature derivative
      // (T * d/dT) of solute.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   partial molar volume temperature derivative      
      double partialmolarvolume_dt();

      // Sets all solute partial charges to zero, resets MDIIS and wipes out
      // working memory.
      // IN:
      //   this :: rism3d object         
      void unsetcharges();

      // Sets all solute partial charges to to their original
      // values. (Undoes rism3d_unsetCharge().)
      // IN:
      //   this :: rism3d object
      void resetcharges();

      // THE FOLLOWING TWO FUNCTIONS ARE DEFINED UNDER "PRIVATE" COMMENT ON RISM3D_C.F90, BUT ARE CALLED FROM AMB_RISM_INTER

      // Transfer a density map of one site of a solvent molecule on
      // to another. E.g., this can be used to localize the hydrogen
      // energy distribution of water onto the oxygen site to which it is
      // bonded. By combining these transferred densities with that of a
      // central site, a molecule distribution can be constructed.
      //
      // The transfer is for local contributions only and makes use of the
      // intramolecular correlation function.  This is approximate and
      // will not produce meaningful results for densities with
      // significant contributions from the solute core region, such as
      // the direct correlation function. Note that the excluded volume is
      // zeroed for all sites, including the center site.  This prevents
      // non-local contributions where there are no sites being broadcast
      // onto physical densities.
      //
      // E.g., for the excess chemical potential, integrating over the
      // resulting distributions will not produce the same result as for
      // integrating the original site distributions. The result will be
      // lower value, much lower than that given by the Universal
      // Correction, which reduces but does not eliminate the positive
      // contribution from the excluded volume.
      //
      // For each non-center site, the following transform is applied
      //
      //  g_c(r) \times ( w_c,nc(r) * dens_nc(r) )
      //
      // where 'c' is the central site, 'nc' is the non-central site,
      // g_c(r) is the pair distribution function of the central site,
      // w_c,nc is the intramolecular correlation function, dens_nc(r) is
      // the density function of interest, \times is element-wise
      // multiplication and * is convolution.
      //
      // TODO:
      // * figure out molecule from the center_site
      // * CLI interface
      //
      // @param[inout] this rism3d object.
      // @param[inout] thermo_map thermodynamic density map of all
      //     sites. Sites attached to the central site are
      //     modified. Summing over all sites of the species will provide
      //     the molecular distribution.
      // @param[in] center_site (not implemented) the solvent site to use
      //     as the molecular center
      void map_site_to_site_flat(double* thermo_map_flat, int center_site);

      // 3D array version of rism3d_map_site_to_site
      void map_site_to_site_3D(double* thermo_map, int center_site);

      // Perform some internal tests for correctness. Some tests require a
      // converged standard [and temperature derivative] solution. A
      // warning message will report when tests cannot be run.
      void selftest();

      // Create an electron density map from a 3D solute-solvent RDF by
      // smearing the 3D RDF with a 1D solvent electron denisty RDF.
      // TODO: Currently only water oxygen is supported for electron
      // smearing.
      // @param[in] this rism3d object.
      // @param[in] electronRDF Solvent 1D electron density map.
      // @param[in] totalSolventElectrons Total electrons in solvent. Ex.: for water, Z = 10.
      // @param[out] electronMap Resulting smeared electron density map.
      void createelectrondensitymap(int iv, double* electronRDF, double electronRDFGridSpacing,
        int totalSolventElectrons, double density, double* electronmap_ptr);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      //// END... FUNCTIONS CALLED FROM AMBER_RISM_INTERFACE.F90
      ////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Just a test to print some members initialized in the constructor
      void PrintTests();

      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////// FUNCTIONS/SUBROUTINES FROM RISM3D_C.F90, BUT NOT
      ////////////////////// CALLED FROM AMBER_RISM_INTERFACE.F90
      ////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Check if we can calculate molecular reconstructions
      // IN:
      //   this : rism3d object
      // OUT:
      //   .true. if we can, .false. if we can't
      bool rism3d_canCalc_molReconstruct();

      // Sets the closure type.
      // IN:
      //   this :: rism3d object
      //   closure :: closure type (see closure enumeration).
      // MY COMMENT: IN .F90 PUBLIC LIST
      void setclosure(string type);

      // Sets the cut off distance for periodic potential and force calculations.
      // IN:
      //   this :: rism3d object
      //   cut     :: distance cutoff for potential and force calculations
      // MY COMMENT: IN .F90 PUBLIC LIST
      void setcut(double cut);

      // Sets MDIIS parameters
      // IN:
      //   this :: rism3d object!
      //   nvec :: number of MDIIS vectors (previous iterations) to keep
      //   del :: scaling factor (step size) applied to estimated gradient (residual)
      //   method :: which implementation of the algorithm
      //   restart :: restart threshold factor. Ratio of the current residual to the
      //              minimum residual in the basis that causes a restart
      // MY COMMENT: IN .F90 PUBLIC LIST
      void rism3d_setmdiis(); 

      // Calculate the total excess chemical potential of solvation with the Gaussian
      // fluctuation correction
      // IN:
      //   this :: rism3d object with computed solution
      //   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
      // OUT:
      //    total Gaussian fluctuation excess chemical potential of solvation
      // MY COMMENT: IN .F90 PUBLIC LIST
      double rism3d_excessChemicalPotentialGF_tot();
      
      // Calculate the total solvation interaction energy: de = density sum g*u for
      // each solvent site.  I.e., the direct intection potential energy of
      // solute and solvent and not the total solvation energy (see solvationEnergy).
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //    the total solvent-solute potential energy [kT]      
      double rism3d_solventPotEne_tot();

      // Calculate the total solvation energy dE = dmu + dTS which includes both
      // solute-solvent interaction energy and the change in solvent
      // self-energy due to rearranging around the solvent.
      // IN:
      //   this :: rism3d object with computed solution
      //   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
      // OUT:
      //    total solvation energy
      double rism3d_solvationEnergy_tot();

      // Calculate the total solvation energy dE = dmu + dTSthe Gaussian
      // fluctuation correction, which includes both solute-solvent
      // interaction energy and the change in solvent self-energy due to
      // rearranging around the solvent.
      // IN:
      //   this :: rism3d object with computed solution
      //   o_lr   :: (optional) (default = .true.) Apply asymptotic long range correction
      // OUT:
      //    total solvation energy
      double rism3d_solvationEnergyGF_tot();

      // Returns a 3D map of the partial molar volume temperature derivative.
      // Integrating this map gives the total PMV as returned from
      // rism3d_partialMolarVolume_tot().  Memory is allocated into a
      // pointer and must be freed by the calling function. For MPI, only
      // the grid points local to this process are allocated and
      // calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the partial molar volume contributions
      // SIDEEFFECTS:
      //   memory is allocated for the grid      
      double* rism3d_partialMolarVolume_tot_map();

      // Returns a 3D map of the partial molar volume .
      // Integrating this map gives the total PMV as returned from
      // rism3d_partialMolarVolume_tot().  Memory is allocated into a
      // pointer and must be freed by the calling function. For MPI, only
      // the grid points local to this process are allocated and
      // calculated.
      // IN:
      //   this :: rism3d object with computed solution
      // OUT:
      //   a 3D-grid of the partial molar volume contributions
      // SIDEEFFECTS:
      //   memory is allocated for the grid
      double* rism3d_partialMolarVolume_dT_tot_map();

      // Write to log file whether or not the test passes within the
      // prescribed error. This is a common utility function for
      // selftest().
      //
      // @param[in] this           rism3d object
      // @param[in] grid           array to integrate.  Assumes gridspacing of
      //                           this%grid%voxelVolume
      // @param[in] reference      pre-integrated reference value
      // @param[in] errorTolerance allowable error
      // @param[in] description    description of the tested quantity
      //
      // @returns .true. if passed.
      bool rism3d_integrateCompare();

      // Copy solvent defined type from fortran to a c++ struct
      // IN:
      //    *arg: pointer to fortran solvent defined type
      //    solv: empty solvent structure
      // OUT:
      //    filled solvent structure
      solvent set_solvent(solvent *arg, solvent solv);   

      // Copy solute defined type from fortran to a c++ struct
      // IN:
      //    *arg: pointer to fortran solute defined type
      //    solu: solute structure
      // OUT:
      //    filled solute structure
      solute set_solute(solute *arg, solute solu);

      // Void function to test access to class members:
      void test_access();

      //////////////////////////////////////////////////// FUNCTIONS FOR WRITING VOLUMETRIC DATA

      void opendx_write_cpp(int atomType, string *file);
      void mrc_map_write_cpp(int atomType, string *file);
      void xyzv_write_cpp(int atomType, string *file);

    };

}

#endif // CLASSES_HPP