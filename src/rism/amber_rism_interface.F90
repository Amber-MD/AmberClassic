!<compile=optimized>

#include "../include/dprec.fh"
#include "../include/assert.fh"

!> Book-keeping object for handling thermodynamic properties
!! calculated in amber_rism_interface.F90.  This stores the basic
!! properties, performs MPI reduction for distributed calculations and
!! can do basic global operations.
!!
!! Unset properties are flagged with a HUGE(1d0) value. This is
!! preserved through all operations. E.g., if partialMolarVolume is
!! HUGE(1D0) on one process, after reduction it is HUGE(1d0) on the
!! master process.
module rismthermo_c
  use safemem
  use rism_report_c

  !> Derived type to store thermodynamic results.  To facilitate MPI
  !! communication, many values are stored in a common array.
  !! Pointers, however, are used to access the data.
  type rismthermo_t

     ! Thermodynamic values that use mpi_buffer for memory.

     !> Excess chemical potential.
     _REAL_, pointer :: excessChemicalPotential(:) => NULL()
     !> Solvent potential energy
     _REAL_, pointer :: solventPotentialEnergy(:) => NULL()
     !> Solvation energy.
     _REAL_, pointer :: solvationEnergy(:) => NULL()
     !> Partial molar volume
     _REAL_, pointer :: partialMolarVolume => NULL()
     !> Total number of particles.
     _REAL_, pointer :: totalParticlesBox(:) => NULL()
     !> Excess number of particles.
     _REAL_, pointer :: excessParticlesBox(:) => NULL()
     !> Kirkwood-Buff integral (a.k.a. total correlation function
     !! integral).
     _REAL_, pointer :: kirkwoodBuff(:) => NULL()
     !> Direct correlation function integral.
     _REAL_, pointer :: DCFintegral(:) => NULL()

     integer :: mpirank = 0, mpisize = 1, mpicomm = 0

     !> All the results are stored in a single array so
     !! the necessary reductions are done in a single
     !! communication.  This is also used for the serial
     !! calculation for simplicity.
     !! Order and size:
     !!         excessChemicalPotential(solvent%numAtomTypes), 
     !!         solventPotentialEnergy(solvent%numAtomTypes),
     !!         solvationEnergy(solvent%numAtomTypes),
     !!         excessNum(solvent%numAtomTypes), 
     !!         kirkwoodBuff(solvent%numAtomTypes),
     !!         dcfi(solvent%numAtomTypes),
     !!         partialMolarVolume(1)
     _REAL_, private, pointer :: mpi_buffer(:) => NULL()
  end type rismthermo_t
  
contains

  !> Initializes and allocates memory for rismthermo objects.
  !! @param[in,out] this The rismthermo object.
  !! @param[in] nsite Number of solvent sites.
  !! @param[in] o_mpicomm  MPI communicator (optional).
  !! @sideeffects Allocates memory and sets value to huge.
  subroutine rismthermo_new(this, nsite, o_mpicomm)
#ifdef MPI
    use mpi
#endif
    implicit none
    type(rismthermo_t), intent(inout) :: this
    integer, intent(in) :: nsite
    integer, optional, intent(in) :: o_mpicomm
#ifdef MPI
    integer ::err
    if (present(o_mpicomm)) &
         this%mpicomm = o_mpicomm
    if (this%mpicomm == MPI_COMM_NULL) &
         call rism_report_error("RISMTHERMO: received NULL MPI communicator")
    call mpi_comm_rank(this%mpicomm, this%mpirank, err)
    if (err /= 0) call rism_report_error&
         ("(a, i8)", "RISMTHERMO: could not get MPI rank for communicator ", this%mpicomm)
    call mpi_comm_size(this%mpicomm, this%mpisize, err)
    if (err /= 0) call rism_report_error&
         ("(a, i8)", "RISMTHERMO interface: could not get MPI size for communicator ", this%mpicomm)
#else
    this%mpicomm=0
    this%mpisize=1
    this%mpirank=0
#endif /*MPI*/

    ! Setup memory space.
    this%mpi_buffer => safemem_realloc(this%mpi_buffer, 7*nsite + 2)
    ! Initialize for case that some values are not calculated.
    call rismthermo_reset(this)
    this%excessChemicalPotential => this%mpi_buffer(1:nsite)
    this%solventPotentialEnergy => this%mpi_buffer(nsite + 1:2*nsite)
    this%solvationEnergy => this%mpi_buffer(2*nsite + 1:3*nsite)
    this%excessParticlesBox => this%mpi_buffer(3*nsite + 1:4*nsite)
    this%kirkwoodBuff => this%mpi_buffer(4*nsite + 1:5*nsite)
    this%DCFintegral => this%mpi_buffer(5*nsite + 1:6*nsite)
    this%totalParticlesBox => this%mpi_buffer(6*nsite + 1:7*nsite)
    this%partialMolarVolume => this%mpi_buffer(7*nsite + 1)
  end subroutine rismthermo_new


!!!Resets the value of all thermodynamics to huge.  This is the flag
!!!indicating no value has been calculated
!!!IN:
!!! this : the rismthermo object
!!!SIDEEFFECTS:
!!! sets all values to huge
  subroutine rismthermo_reset(this)
    implicit none
    type(rismthermo_t), intent(inout) :: this
    this%mpi_buffer = huge(1d0)
  end subroutine rismthermo_reset

!!!Does an MPI reduction on all thermodynamic values
!!!IN:
!!! this : the rismthermo object
!!!SIDEEFFECTS:
!!! The 0 node has the total for each values.  If this is a distruted
!!! calculation, MPI is used to do this.  For serial runs this does
!!! nothing.
  subroutine rismthermo_mpireduce(this)
#ifdef MPI
    use mpi
#endif
    implicit none
    type(rismthermo_t), intent(inout) :: this
#ifdef MPI
    _REAL_ :: buffer_copy(ubound(this%mpi_buffer, 1))
    integer :: err
    buffer_copy = this%mpi_buffer
    if (this%mpirank == 0) then
       call mpi_reduce(MPI_IN_PLACE, this%mpi_buffer, &
            ubound(this%mpi_buffer, 1), MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, this%mpicomm, err)
    else
       call mpi_reduce(this%mpi_buffer, this%mpi_buffer, &
            ubound(this%mpi_buffer, 1), MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, this%mpicomm, err)
    end if
    where(buffer_copy == huge(1d0))
       this%mpi_buffer = huge(1d0)
    end where
#endif /*MPI*/
  end subroutine rismthermo_mpireduce

  
!!!Deallocates and nullifies the rismthermo objects
!!!IN:
!!! this : the rismthermo object
!!!SIDEEFFECTS:
!!! deallocates memory and nullifies all pointers
  subroutine rismthermo_destroy(this)
    implicit none
    type(rismthermo_t), intent(inout) :: this
    if (safemem_dealloc(this%mpi_buffer) /= 0 ) &
         call rism_report_error("Dealloc failed in rism_thermo")
    nullify(this%excessChemicalPotential)
    nullify(this%solventPotentialEnergy)
    nullify(this%solvationEnergy)
    nullify(this%partialMolarVolume)
    nullify(this%excessParticlesBox)
    nullify(this%totalParticlesBox)
    nullify(this%kirkwoodBuff)
    nullify(this%DCFintegral)
    this%mpicomm=0
    this%mpisize=1
    this%mpirank=0
  end subroutine rismthermo_destroy
end module rismthermo_c

!> Module to hold global data and preserve namespace.
!!
!! General interface between 3D-RISM and SANDER/SFF in the Amber suite.  Except
!! where noted, all data and subroutines may be called from either SFF or SANDER.
!!
!! To make this file interoperable with C while still maintaining the Fortran 95
!! standard, we must not put subroutines or functions inside of modules. However,
!! some variables must be available at the global scope either to be accessed
!! from outside routines (e.g. print statements) or be shared between local
!! routines (e.g. instances of derived types and parameters for the run).  This
!! is accomplished with two modules. AMBER_RISM_INTERFACE is always compiled and
!! contains RISM specific variables.  SANDER_RISM_INTERFACE is created for
!! SANDER but not for SFF. In the Fortran world (SANDER), the
!! SANDER_RISM_INTERFACE module provides access to subroutines and functions
!! with all the benefits of a module. In the C world, there is no module and
!! functions/subroutines may be directly called.
!!
!! Setup and initialization (serial and MPI):
!! The method of setting up and initializing is somewhat flexible and complicated
!! in order to maintain the correct output initialization and output sequence
!! of SANDER.  SANDER reads all of the input files and then prints a summary
!! from the master node.  RISM must do the same, so RISM_SETPARAM must be called
!! from the master node in SANDER.  However, it is safe to call it from all
!! nodes at the same time with the caveat that irism is defined and the same on
!! all nodes; this done in SFF.  Initializing the calculation
!! must be done in parallel so RISM_INIT must be call from all processes.
!!
!! Note that it is always safe to call RISM_SETPARAM and RISM_INIT as
!! long as irism is define on the master node in the relevent data
!! structure.
!!
!! To summarize:
!!
!! In SANDER:
!!
!!  if (master) then
!!    call rism_setparam(mdin, &
!!         commsander, &
!!         igb, numAtoms, ntypes, x(L15:L15 + numAtoms-1), &
!!         x(LMASS:LMASS + numAtoms-1), cn1, cn2, &
!!         ix(i04:i04 + ntypes**2-1), ix(i06:i06 + numAtoms-1))
!!  endif
!!  call rism_init(commsander)
!!
!! MPI and other subroutines and functions: All public functions are
!! safe to call from all nodes and usually must be.  Only
!! RISM_THERMO_PRINT and RISM_SETPARAM must exclusively be called by
!! the master.
module amber_rism_interface
  use rism3d_c
  use rism3d_solvent_c
  use rism3d_solute_c
  use rism_report_c
  use safemem
  use rismthermo_c

  !> Parameter derived type for storing calculation parameters. 
  type rismprm_t
     sequence
     !> Cutoff for rism calculations (separate from SANDER non-bond).
     _REAL_ :: solvcut
     !> Grid spacing for all of the grids.
     _REAL_ :: grdspc(3)
     !> 'Step size' for MDIIS.
     _REAL_ :: mdiis_del
     !> Restart threshold factor. Ratio of the current residual to the
     !! minimum residual in the basis that causes a restart.
     _REAL_ :: mdiis_restart
     !> Charge smearing parameter for Ewald,
     !! typically eta in the literature
     _REAL_ :: chargeSmear
     !> For backwards compatibility, we still need to read this.
     integer :: closureOrder
     !> Number of grid points in each dimension.
     integer :: ng3(3)
     !> Use 3D-RISM.
     integer :: rism
     !> Number of vectors used for MDIIS (consequently, the number of
     !! copies of CUV we need to keep for MDIIS).
     integer :: mdiis_nvec
     !> Mdiis implementation to use.
     integer :: mdiis_method
     !> Maximum number of rism relaxation steps.
     integer :: maxstep
     !> Number of past cuv time steps saves.
     integer :: npropagate
     !> 0 - Do nothing.
     !! 1 - Redistribute forces to get zero total force.
     integer :: zerofrc
     !> If 0, the 3D-RISM solution is calculated but the forces are not.
     integer :: apply_rism_force
     !> Size of rism multiple timestep.
     integer :: rismnrespa

     !> Save itermediate results every saveprogress interations (0
     !! means no saves).
     integer :: saveprogress
     integer :: ntwrism
     !> 0 - no ouput.
     !! 1 - memory allocation and steps for convergence.
     !! 2 - 1 + convergence progress.
     integer :: verbose
     !> Print parameter for relaxation steps every progress iteration
     !! (0 means no saves).
     integer :: progress

     !> Calculate and print out thermodynamics.  This is primarily used
     !! by sander but also serves as padding for alignment for NAB.
     integer :: write_thermo

     !> BPR: Make sure the number of INTEGERs is not odd, to shut up a
     !! compiler whinge about misalignment.
     !! Note: this should be commented if ever there are grounds for a
     !! new real INTEGER.
     integer :: padding
  end type rismprm_t
  
  !> Possible RISM calculation types for MD.
  !! RISM_NONE :: no RISM calculation
  !! RISM_FULL :: full RISM solution
  !! RISM_INTERP :: interpolation
  integer, parameter :: RISM_NONE=0, RISM_FULL=1, RISM_INTERP=2

  type(rismprm_t), save :: rismprm
  type(rismthermo_t), save :: rismthermo
  type(rismthermo_t), save :: rismthermo_pol
  type(rismthermo_t), save :: rismthermo_apol
  type(rism3d), save :: rism_3d
  type(rism3d_solvent), save :: solvent
  type(rism3d_solute), save :: solute
  _REAL_ :: centerOfMass(3)
  ! Read from prmtop file during parameter processing and passed to
  ! child MPI processes.
  _REAL_ :: unitCellDimensions(6)

  integer :: outunit

  !> List of closures to use in order.  Only the last closure is used
  !! for thermodynamic output.  This can be used to progressively
  !! increase the order of the closure to aid convergence.  Closure
  !! types: KH, HNC or PSEn where n is an integer. This is initialized
  !! to a length of 10 in defaults() to allow it to be used with the
  !! namelist in sander.
  character(len=8), pointer :: closurelist(:) => NULL()
  !> Name of potential used for periodic calculations.
  character(len=8) :: periodicPotential = 'pme'
  !> Residual tolerance for the solution of each closure in the
  !! list. On input this can be of length one, two or
  !! size(closurelist). If length one, use this value for the final
  !! closure and the default for all others (see sanity_check()). If
  !! length two, use the last value for the last closure and the first
  !! value for all other closures. Otherwise, match each value to each
  !! closure.
  _REAL_, pointer :: tolerancelist(:) => NULL()
  !> Default number of closures to start with in the list. Namelists
  !! don't support pointers in many ways so we need to initialize the
  !! closure and tolerance lists to some reasonable size.
  integer, parameter :: nclosuredefault = 10

  !I/O file names:
  !xvvfile       : (input) site-site solvent susceptibility from RISM1D (.xvv)
  !guvfile       : (output) pair distribution function.  Volumetric file.
  !huvfile       : (output) total correlation function.  Volumetric file.
  !cuvfile       : (output) direct correlation function.  Volumetric file.
  !uuvfile       : (output) solute-solvent potential.  Volumetric file.
  !quvfile       : (output) charge density distribution.  Volumetric file.
  !chgDistFile   : (output) charge distribution. Volumetric file.
  !excessChemicalPotentialfile    : (output) excess chemical potential map [kcal/mol/A^3]. Volumetric file.
  !solvationEnergyfile   : (output) solvation energy map [kcal/mol/A^3]. Volumetric file.
  !entropyfile   : (output) solvent entroy (-TS) map [kcal/mol/A^3]. Volumetric file.
  !solventPotentialEnergyfile     : (output) solvent-solute potential energy map [kcal/mol/A^3]. Volumetric file.
  !volfmt        : either 'mrc', 'dx', or 'xyzv'
  character(len=256) :: xvvfile='', guvfile='', huvfile='', cuvfile='', &
       uuvfile='', quvFile='', chgDistFile='', &
       excessChemicalPotentialfile='', solvationEnergyfile='', entropyfile='', &
       solventPotentialEnergyfile='', volfmt='mrc', crdFile=''

  integer :: mpirank = 0, mpisize = 1, mpicomm = 0

  !working memory for rism_force() so it is not reallocated every time
  !ff :: forces
  _REAL_, pointer :: ff(:, :) => NULL()

  private :: rism_mpi_bcast

end module amber_rism_interface

module sander_rism_interface
  use amber_rism_interface
  implicit none
contains

  !> Sets all input parameters for 3D-RISM.  This _must_ be called by the head
  !! node and may be called by all nodes.  If all nodes call this
  !! subroutine they must all agree on the value of rismprm%rism (SANDER) 
  !!
  !! SANDER prerequisites:
  !!   - Names of 3D-RISM specific I/O files (Xvv, Guv, etc.) are specified on the
  !!    command line and are read by mdfil.F90 into the variables in
  !!    AMBER_RISM_INTERFACE.
  !!   - igb should be set to 6 (vacuum electrostatics).
  !! SANDER IN:
  !! @param[in] mdin Name of the mdin file that SANDER name lists are read from.
  !!
  !! @param[in] userdata rismprm_t C struct equivalent with use options.
  !! @param[in] ntol Number of tolerances in the array.
  !! @param[in] tol Array of tolerances read in from the user.
  !! @param[in] closurelen Length of the closurechar strings.
  !! @param[in] nclosure Number of closures read in.
  !! @param[in] closurechar An array of nclosure strings of closurelen characters.
  !! @param[in] xvvlen Length of xvvchar array.
  !! @param[in] xvvchar Character array for Xvv input file name.
  !! @param[in] guvlen Length of guvchar array.
  !! @param[in] guvchar Character array for Guv output file name.
  !! @param[in] huvlen Length of huvchar array.
  !! @param[in] huvchar Character array for Huv output file name.
  !! @param[in] cuvlen Length of cuvchar array.
  !! @param[in] cuvchar Character array for Cuv output file name.
  !! @param[in] uuvlen Length of uuvchar array.
  !! @param[in] uuvchar Character array for Uuv output file name.
  !! @param[in] quvlen Length of quvchar array.
  !! @param[in] quvchar Character array for Quv output file name.
  !! @param[in] chgDistlen Length of chgDistchar array.
  !! @param[in] chgDistchar Character array for charge distribution
  !!   output file name.
  !! @param[in] excessChemicalPotentiallen Length of excessChemicalPotentialchar array.
  !! @param[in] excessChemicalPotentialchar Character array for excessChemicalPotential output file name.
  !! @param[in] solvationEnergylen Length of solvationEnergychar array.
  !! @param[in] solvationEnergychar Character array for solvationEnergy output file name.
  !! @param[in] entropylen Length of entropychar array.
  !! @param[in] entropychar Character array for entropy output file name.
  !! @param[in] solventPotentialEnergylen Length of solventPotentialEnergychar array.
  !! @param[in] solventPotentialEnergychar Character array for solute-solvent potential
  !!   energy output file name.
  !! @param[in] volfmtlen Length of volfmtchar array.
  !! @param[in] volfmtchar Character array for the format type for
  !!   volumetric data.
  !! @param[in] periodiclen Length of periodicchar array.
  !! @param[in] periodicchar Character array for the abbreviated label
  !!   of the periodic potential.
  !!
  !! IN:
  !! @param[in] comm  MPI communicator.
  !! @param[in] numAtoms Number of solute atoms.
  !! @param[in] numTypes Number of atom solute types.
  !! @param[in] charge Solute atom partial charges in Amber units.
  !! @param[in] mass Solute atom masses [AU].
  !! @param[in] ljA Lennard-Jones A parameter for each solute atom type pair.
  !! @param[in] ljB Lennard-Jones B parameter for each solute atom type pair.
  !! @param[in] atomTypeIndex Solute atom type index.
  !! @param[in] nonbondedParmIndex Solute nonbonded parameter index.
  subroutine rism_setparam( &
       mdin, &
       comm, &
       numAtoms, numTypes, &
       charge, mass, ljA, ljB, &
       atomTypeIndex, nonbondedParmIndex)
    use amber_rism_interface
    use binrestart, only : readUnitCellDimensionsFromCrd
#ifdef OPENMP
    use constants_rism, only : KB, omp_num_threads
#else
    use constants_rism, only : KB
#endif
#ifdef MPI
    use mpi
#endif
    implicit none
    integer, intent(in) :: numAtoms, numTypes, atomTypeIndex(numTypes**2), nonbondedParmIndex(numAtoms)
    _REAL_, intent(in) :: charge(numAtoms), mass(numAtoms), ljA(numTypes*(numTypes + 1)/2), &
         ljB(numTypes*(numTypes + 1)/2)
    character(*), intent(in) :: mdin
    integer :: mdin_unit=55
    integer, intent(in) :: comm
    character(len=16) :: whtspc
    integer :: i, stat, err
    integer :: iclosure, iclosurechar, un, ier
    logical :: op
    character(len=5) omp_threads

    write(whtspc, '(a16)')" "

#ifdef MPI
    mpicomm = comm
    if (mpicomm == MPI_COMM_NULL) &
         call rism_report_error("RISM3D interface: received NULL MPI communicator")
    call mpi_comm_rank(mpicomm, mpirank, err)
    if (err /= 0) call rism_report_error &
         ("(a, i8)", "RISM3D interface: could not get MPI rank for communicator ", mpicomm)
    call mpi_comm_size(mpicomm, mpisize, err)
    if (err /= 0) call rism_report_error &
         ("(a, i8)", "RISM3D interface: could not get MPI size for communicator ", mpicomm)
#endif /*MPI*/

    ! If this is not a RISM run, we're done.
    if (rismprm%rism == 0) return

    ! Rank 0 only.
    if (mpirank /= 0) return

#ifndef API
    call defaults()
    outunit = rism_report_getMUnit()
    inquire(file=mdin, opened=op, number=un)
    if (op) mdin_unit=un
    open(unit=mdin_unit, file=mdin, status='OLD', form='FORMATTED', iostat=stat)
    if (stat/= 0) then
       call rism_report_error('(a, i4)', "opening "//trim(mdin)//"failed. IOSTAT=", stat)
    end if
    call read_namelist(mdin_unit)
    if (.not.op) close(unit=mdin_unit)
#endif

    ! Initialize 3D-RISM solute and solvent.

#ifdef OPENMP
    ier = fftw_init_threads()
    if( ier == 0 ) then
       write(0,*) 'failure in fftw_plan_with_nthreads'
       stop 1
    else
       write(6,'(a,i2,a)') '| calling fftw_plan_with_nthreads(', &
          omp_num_threads,')'
    end if
    call fftw_plan_with_nthreads(omp_num_threads)
#endif

    call rism3d_solvent_new(solvent, xvvfile)    
    call rism3d_solute_new_sander(solute, numAtoms, numTypes, atomTypeIndex, &
         nonbondedParmIndex, charge, ljA, ljB, mass, solvent%temperature)

    call readUnitCellDimensionsFromCrd(crdFile, unitCellDimensions)

    call sanity_check()
    
#ifdef API
    if (rismprm%rism >= 1 .and. rismprm%verbose > 0) then
#else
    if (rismprm%rism >= 1) then
#endif
       write(outunit, '(a)') "3D-RISM:"
       if (rismprm%rism < 1) then
          write(outunit, '(5x, a, i10)') 'irism   =', rismprm%rism
       else if (rismprm%rism == 1) then
          write(outunit, '(5x, 3(a10, "=", 100a10))') &
               'closure'//whtspc, closurelist
          write(outunit, '(5x, 2(a10, "=", f10.5))') &
               'solvcut'//whtspc, rismprm%solvcut
          write(outunit, '(5x, a10, "=", 3(f10.5, 1x))') &
               'grd_spc'//whtspc, rismprm%grdspc
          write(outunit, '(5x, a10, "=", 3(i10, 1x))') &
               'ng3'//whtspc, rismprm%ng3
          write(outunit, '(5x, a10, "=", 1p, 100e10.2)')  &
               'tolerance'//whtspc, tolerancelist
          write(outunit, '(5x, a10, "=", f10.5, a10, "=", i10)')  &
               'mdiis_del'//whtspc, rismprm%mdiis_del, &
               ', mdiis_nvec'//whtspc, rismprm%mdiis_nvec
          write(outunit, '(5x, a10, "=", i10, a10, "=", 1p, e10.2)') &
               'mdiis_method'//whtspc, rismprm%mdiis_method, &
               ', mdiis_restart'//whtspc, rismprm%mdiis_restart
          write(outunit, '(5x, 3(a10, "=", i10))') &
               'maxstep'//whtspc, rismprm%maxstep, &
               ', npropagate'//whtspc, rismprm%npropagate
          write(outunit, '(5x, a10, "=", i10)') &
               'zerofrc'//whtspc, rismprm%zerofrc
          write(outunit, '(5x, a10, "=", i10)') &
               'apply_rism_force'//whtspc, rismprm%apply_rism_force
          write(outunit, '(5x, a10, "=", i10)') &
               'rismnrespa'//whtspc, rismprm%rismnrespa
          write(outunit, '(5x, a15, "=", i5, a10, "=  ", a8)') &
               'write_thermo'//whtspc, rismprm%write_thermo, &
               ', volfmt'//whtspc, volfmt
          write(outunit, '(5x, 3(a15, "=", i5))') &
               'saveprogress'//whtspc, rismprm%saveprogress, &
               ', ntwrism'//whtspc, rismprm%ntwrism, &
               ', verbose'//whtspc, rismprm%verbose
          write(outunit, '(5x, 3(a10, "=", i10))') &
               'progress'//whtspc, rismprm%progress
          write(outunit, '(5x, a14, "=", f6.3)') &
               'chargeSmear'//whtspc, rismprm%chargeSmear
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       end if
       call flush(outunit)
    end if
    return

  end subroutine rism_setparam

  !> Performs all of the initialization required for 3D-RISM
  !! calculations. All parameters should have already been set with
  !! RISM_SETPARAM. For MPI calculations this _must_ be called by all
  !! processes. Only the head node irism value is used.
  !! @param[in] comm MPI communicator.
  subroutine rism_init(comm)
    use amber_rism_interface
    use safemem
#ifdef MPI
    use mpi
#endif
    implicit none
    integer, intent(in) :: comm

    integer :: err

    call rism_report_setMUnit(6)
    call rism_report_setWUnit(6)
    call rism_report_setEUnit(6)

#ifdef MPI
    mpicomm = comm
    if (mpicomm == MPI_COMM_NULL) &
         call rism_report_error("RISM3D interface: received NULL MPI communicator")
    call mpi_comm_rank(mpicomm, mpirank, err)
    if (err /= 0) call rism_report_error&
         ("(a, i8)", "RISM3D interface: could not get MPI rank for communicator ", mpicomm)
    call mpi_comm_size(mpicomm, mpisize, err)
    if (err /= 0) call rism_report_error&
         ("(a, i8)", "RISM3D interface: could not get MPI size for communicator ", mpicomm)
    call rism_mpi_bcast(mpirank, mpisize, mpicomm)
#endif /*MPI*/

    ! STOP HERE IF THIS IS NOT A RISM RUN.
    if (rismprm%rism == 0) return

    ! 3D-RISM may have already been initialized. In the absence of a
    ! subroutine to set all of these parameters individually, we
    ! destroy the original instance and re-initialize. Since this is a
    ! rare event, it should not add to the expense of the calculation
    ! in any practical way.
    call rism3d_destroy(rism_3d)
    if (rismprm%ng3(1) == -1) then
       call rism3d_new(rism_3d, solute, solvent, rismprm%npropagate, &
          closurelist, rismprm%solvcut, &
          rismprm%mdiis_nvec, rismprm%mdiis_del, rismprm%mdiis_method, &
          rismprm%mdiis_restart, &
          rismprm%chargeSmear, &
          o_grdspc=rismprm%grdspc, o_mpicomm=mpicomm, &
          o_periodic=periodicPotential,&
          o_unitCellDimensions=unitCellDimensions)
    else
       call rism3d_new(rism_3d, solute, solvent, rismprm%npropagate, &
          closurelist, rismprm%solvcut, &
          rismprm%mdiis_nvec, rismprm%mdiis_del, rismprm%mdiis_method, &
          rismprm%mdiis_restart, &
          rismprm%chargeSmear, &
          o_ng3=rismprm%ng3, o_mpicomm=mpicomm, &
          o_periodic=periodicPotential, &
          o_unitCellDimensions=unitCellDimensions)
    end if
    call rism3d_setverbosity(rism_3d, rismprm%verbose)

    call rismthermo_new(rismthermo, rism_3d%solvent%numAtomTypes, mpicomm)

    ! Allocate working memory.
    ff => safemem_realloc(ff, 3, rism_3d%solute%numAtoms)

    ! Free up a bit of memory.
    call rism3d_solvent_destroy(solvent)
    call rism3d_solute_destroy(solute)

    call flush(outunit)

  end subroutine rism_init


  !> Returns the RISM calculation type for this step.  I.e., no
  !! calculation (RISM_NONE), full RISM solution (RISM_FULL),
  !! interpolation (RISM_INTERP).
  !! @param[in] irespa Respa iteration.
  function rism_calc_type(irespa) result(calc_type)
    use amber_rism_interface

    implicit none
    integer, intent(in) :: irespa
    integer :: calc_type
    calc_type = RISM_FULL
    ! Test for RESPA.
    if (mod(irespa, rismprm%rismnrespa) /= 0) then
       calc_type = RISM_NONE
       ! write(6,'(a,i4)') '|  setting RISM_NONE for irespa = ', irespa
    end if
  end function rism_calc_type

  
  !> Driver routine for 3D-RISM.
  !! Calculates the 3D-RISM solution, energy and forces.  Determines
  !! if an interpolation step can be used.  If so, the interpolation
  !! code is called, otherwise a full 3D-RISM calculation is
  !! performed.  Energies are only valid for full 3D-RISM
  !! calculations.
  !! @param[in] atomPositions_md Atom positions for solute.
  !! @param[in,out] frc Pre-3D-RISM force.  The 3D-RISM forces are added to this.
  !! @param[in] epol Polarization energy, excess chemical potential.
  !! @param[in] irespa Respa iteration.
  !! @param[in] imin Sander imin value.  If not 0, full thermodynamics
  !!     are calculated immediately after a full RISM solution. This
  !!     is to work around the imin=1,5 MPI implementation in sander,
  !!     which calls force() in a infinite loop on the non-master
  !!     processes.
  subroutine rism_force(atomPositions_md, frc, epol, irespa, imin)
    use amber_rism_interface
    use constants_rism, only : KB
    use rism3d_c, only : rism3d_calculateSolution
    use rism_util, only: corr_drift
#ifdef MPI
    use mpi
#endif
    implicit none
#include "def_time.h"


    integer, intent(in) :: irespa
    _REAL_, intent(in) :: atomPositions_md(3, rism_3d%solute%numAtoms)
    _REAL_, intent(inout) :: frc(3, rism_3d%solute%numAtoms)
    _REAL_, intent(out) :: epol
    integer, intent(in) :: imin
    ! Solvation energy and entropy.
    _REAL_ :: epol_e, epol_ts

    !iclosure :: counter for closures
    integer :: iclosure
    integer :: i, iatu
    integer :: err
    _REAL_ :: mpi_temp, partialMolarVolume

   epol = 0
   ff=0

    if (rism_calc_type(irespa) == RISM_NONE) then
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!No forces this step. DO NOTHING!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
       if (rismprm%verbose >= 2) call rism_report_message("|FULL RISM!!!")
!!!!!!!!!!!!!!!!!!!!!!!!
!!!Full RISM SOLUTION!!!
!!!!!!!!!!!!!!!!!!!!!!!!
       call rism3d_setCoord(rism_3d, atomPositions_md)
       call rism3d_calculateSolution(rism_3d, rismprm%saveprogress, &
            rismprm%progress, rismprm%maxstep, tolerancelist, &
            rismprm%ng3, rismprm%verbose)
       if(imin /= 0) then
          call rism_solvdist_thermo_calc(.false., 0)
       end if
       ! Ugly, ugly hack.  SANDER runs force on the first frame twice.
       ! This messes up the solution propagation.  Here we set the
       ! solution counter back one to ignore one of the duplicate
       ! solutions.
       if (irespa == 1) then
          rism_3d%nsolution = 1
       end if
       if (rismprm%apply_rism_force == 1) then
          call rism3d_force(rism_3d, ff)
          if (rismprm%zerofrc == 1) then
#ifdef MPI
             call corr_drift(ff, rism_3d%solute%mass, rism_3d%solute%numAtoms, &
                  mpirank, mpisize, mpicomm)
#else
             call corr_drift(ff, rism_3d%solute%mass, rism_3d%solute%numAtoms)
#endif /*MPI*/
          end if
          ! Convert to [kcal/mol/A].
          ff = ff * KB * rism_3d%solvent%temperature
       end if

       ! Get the excess chemical potential.
       call timer_start(TIME_EXCESSCHEMICALPOTENTIAL)
       epol = rism3d_excessChemicalPotential_tot(rism_3d)*KB*rism_3d%solvent%temperature
       call timer_stop(TIME_EXCESSCHEMICALPOTENTIAL)

       if (rismprm%rismnrespa >1) then
          if (rismprm%verbose>=2) call rism_report_message("|IMPULSE FORCE!!!")
          ff=rismprm%rismnrespa*ff
       end if

    end if

    if (rismprm%apply_rism_force==1) frc = frc + ff

    call flush(outunit)
  end subroutine rism_force

  
  !> Calculates and stores thermodynamics for 3D-RISM, optionally
  !! printing distribution files.
  !! 
  !! Since performing polar decomposition destroys solvent distributions
  !! for both standard solutions and entropic decompositions,
  !! distributions must be output here.  This is also necessary since
  !! the temperature derivative solution is only calculated here.
  !! 
  !! @param[in] writedist .true. to write distributions.
  !! @param[in] step Time step number.
  !! @sideeffects The first time through, memory is allocated.  This
  !!    must be destroyed later.
  subroutine rism_solvdist_thermo_calc(writedist, step)
    use amber_rism_interface
    use constants_rism, only : kb, COULOMB_CONST_E
#ifdef MPI
    use mpi
#endif
    implicit none
    logical*4, intent(in) :: writedist
    integer, intent(in) :: step
    integer :: err
    call rismthermo_reset(rismthermo)
    if (associated(rismthermo_pol%excessChemicalPotential)) &
         call rismthermo_reset(rismthermo_pol)
    if (associated(rismthermo_apol%excessChemicalPotential)) &
         call rismthermo_reset(rismthermo_apol)

    ! Calculate thermodynamics.
    rismthermo%excessChemicalPotential = &
         rism3d_excessChemicalPotential(rism_3d)* KB * rism_3d%solvent%temperature
    rismthermo%solventPotentialEnergy = rism3d_solventPotEne(rism_3d) * KB * rism_3d%solvent%temperature
    rismthermo%partialMolarVolume = rism3d_partialMolarVolume(rism_3d)
    rismthermo%excessParticlesBox = rism3d_excessParticles(rism_3d)
    rismthermo%totalParticlesBox = rismthermo%excessParticlesBox &
         + rism_3d%grid%voxelVolume &
         * rism_3d%grid%totalLocalPointsR * rism_3d%solvent%density
         ! this is the local volume for the MPI process.
    rismthermo%kirkwoodBuff = rism3d_kirkwoodbuff(rism_3d)
    rismthermo%DCFintegral = rism3d_DCFintegral(rism_3d)

    ! Output distributions.
    if (writedist .and. step >= 0) &
         call rism_writeVolumetricData(rism_3d, step)

    call rismthermo_mpireduce(rismthermo)
  end subroutine rism_solvdist_thermo_calc

  
  !> The format for NAB programs will consist of one line for each
  !! catgory of data: energy, volume and excess.  RISM-specific values
  !! are prefixed with "rism_".  Values at least partly calculated
  !! outside RISM have no prefix.
  !! @param[in] description If true, output a description of the table
  !!                        that will be printed but do not print out
  !!                        any values.
  !! @param[in] soluPot Total and component potential energy of the
  !!                    solute. Expected order is: total, LJ, elec,
  !!                    bond, angle, dihedral, H-bond, LJ-14, elec-14,
  !!                    restraints, 3D-RISM total excess chemical
  !!                    potential.
  subroutine rism_thermo_print(description, soluPot)
    use amber_rism_interface
    use constants_rism, only : kb, COULOMB_CONST_E, PI
    implicit none
    logical*4, intent(in) :: description
    _REAL_, intent(in) :: soluPot(25)
    
    integer :: iv, iu, ig, il, err
    integer :: ix, iy, iz

    if (description) then
       if (mpirank==0) then

          ! Thermodynamics key.

          ! Free energy-type properties.
          call thermo_print_descr_line('solutePotentialEnergy', '[kcal/mol]', 'Total', '', &
               (/'LJ        ', 'Coulomb   ', 'Bond      ', 'Angle     ', 'Dihedral  ', &
               'H-Bond    ', 'LJ-14     ', 'Coulomb-14', 'Restraints', '3D-RISM   '/), 10)
          call thermo_print_descr_line('rism_excessChemicalPotential', '[kcal/mol]', 'Total', 'ExcChemPot_', &
               rism_3d%solvent%atomName, rism_3d%solvent%numAtomTypes)
          call thermo_print_descr_line('rism_solventPotentialEnergy', '[kcal/mol]', 'Total', 'UV_potential_', &
               rism_3d%solvent%atomName, rism_3d%solvent%numAtomTypes)

          ! Thermodynamic properties not related to free energy.
          call thermo_print_descr_line('rism_partialMolarVolume', '[A^3]', 'PMV', '', &
               rism_3d%solvent%atomName, 0)
          call thermo_print_descr_line('rism_totalParticlesBox', '[#]', '', 'TotalNum_', &
               rism_3d%solvent%atomName, rism_3d%solvent%numAtomTypes)
          call thermo_print_descr_line('rism_totalChargeBox', '[e]', 'Total', 'TotalChg_', &
               rism_3d%solvent%atomName, rism_3d%solvent%numAtomTypes)
          call thermo_print_descr_line('rism_excessParticlesBox', '[#]', '', 'ExNum_', &
               rism_3d%solvent%atomName, rism_3d%solvent%numAtomTypes)
          call thermo_print_descr_line('rism_excessChargeBox', '[e]', 'Total', 'ExChg_', &
               rism_3d%solvent%atomName, rism_3d%solvent%numAtomTypes)
          call thermo_print_descr_line('rism_KirkwoodBuff', '[A^3]', '', 'KB_', &
               rism_3d%solvent%atomName, rism_3d%solvent%numAtomTypes)
          call thermo_print_descr_line('rism_DCFintegral', '[A^3]', '', 'DCFI_', &
               rism_3d%solvent%atomName, rism_3d%solvent%numAtomTypes)

       end if
    else
       ! DATA: free energy-based properties.
       if (mpirank==0 .and. associated(rismthermo%excessChemicalPotential)) then

          ! Cast the SANDER array into the order of the NAB array.
          call thermo_print_results_line('solutePotentialEnergy', soluPot(1), (/soluPot(2), &
               soluPot(3), soluPot(5), soluPot(6), soluPot(7), soluPot(13), &
               soluPot(8), soluPot(9), soluPot(10), soluPot(24)/), 10)
          call thermo_print_results_line('rism_excessChemicalPotential', sum(rismthermo%excessChemicalPotential), &
               rismthermo%excessChemicalPotential, rism_3d%solvent%numAtomTypes)
          call thermo_print_results_line('rism_solventPotentialEnergy', &
               sum(rismthermo%solventPotentialEnergy), &
               rismthermo%solventPotentialEnergy, &
               rism_3d%solvent%numAtomTypes)

          ! Non-free energy-based properties.
          call thermo_print_results_line('rism_partialMolarVolume', &
               rismthermo%partialMolarVolume, (/1d0/), 0)
          call thermo_print_results_line('rism_totalParticlesBox', &
               HUGE(1d0), rismthermo%totalParticlesBox, &
               rism_3d%solvent%numAtomTypes)
          call thermo_print_results_line('rism_totalChargeBox', &
               sum(rismthermo%totalParticlesBox * rism_3d%solvent%charge) &
               * sqrt((KB * rism_3d%solvent%temperature) / COULOMB_CONST_E), &
               rismthermo%totalParticlesBox * rism_3d%solvent%charge &
               * sqrt((KB * rism_3d%solvent%temperature) / COULOMB_CONST_E), &
               rism_3d%solvent%numAtomTypes)
          call thermo_print_results_line('rism_excessParticlesBox', &
               HUGE(1d0), rismthermo%excessParticlesBox, rism_3d%solvent%numAtomTypes)
          call thermo_print_results_line('rism_excessChargeBox', &
               sum(rismthermo%excessParticlesBox * rism_3d%solvent%charge) &
               * sqrt((KB * rism_3d%solvent%temperature) / COULOMB_CONST_E), &
               rismthermo%excessParticlesBox * rism_3d%solvent%charge &
               * sqrt((KB * rism_3d%solvent%temperature) / COULOMB_CONST_E), &
               rism_3d%solvent%numAtomTypes)
          call thermo_print_results_line('rism_KirkwoodBuff', HUGE(1d0), &
               rismthermo%kirkwoodBuff, rism_3d%solvent%numAtomTypes)
          call thermo_print_results_line('rism_DCFintegral', HUGE(1d0), &
               rismthermo%DCFintegral, rism_3d%solvent%numAtomTypes)

       end if
    end if
    call flush(outunit)
  end subroutine rism_thermo_print
  
  !> Prints out a description line for thermodynamics output.
  !! @param[in] category Name of the thermodynamic quantity.
  !! @param[in] units    Units of output.
  !! @param[in] total    String to identify the system wide amount.
  !! @param[in] prefix   Prefix attached to each item name.
  !! @param[in] item     Array of decomposition names.  E.g., site names.
  !! @param[in] nitem    Number of items.  Needed for the NAB case
  !!                     where we don't use a module.
  !! @sideeffect Writes to outunit.
  subroutine thermo_print_descr_line(category, units, total, prefix, item, nitem)
    use amber_rism_interface
    implicit none
    integer, intent(in) :: nitem
    character(len=*), intent(in) :: category, units, total, prefix, item(nitem)
    ! Format for category (calculation type, e.g. excess chemical
    ! potential).
    character(len=64) :: catFmt = "(a40)"
    ! Format for category (calculation type, e.g. excess chemical
    ! potential) with comment bar and units.
    character(len=32) :: catbarFMT = "('|', a40, ' ', a10)"
    ! Format for column headings.
    character(len=32) :: titleFmt = "(a21)"
    ! Long string of whitespace that can be used to effect a
    ! left-justified string.  Otherwise strings are right-justified.
    ! Simply concatenate this to the end of the string you wish
    ! left-justified.
    character(len=64) :: whtspc
    integer :: i
    write(whtspc, '(a64)')" "
    write(outunit, catbarFmt, advance='no') trim(category)//whtspc, trim(units)//whtspc
    write(outunit, titleFmt, advance='no') trim(Total)
    do i=1, nitem
       write(outunit, titleFmt, advance='no') trim(prefix)//trim(item(i))
    end do
    write(outunit, '(a)')
  end subroutine thermo_print_descr_line

  
  !> Prints out a results line for thermodynamics output. Printing of
  !! the total is supressed if the value HUGE() is passed in.  Print of
  !! everything is supressed the total is HUGE() and nitem==0 or if any
  !! values in item are HUGE(). This can be used to supress output for
  !! undefined values at runtime.
  !! @param[in] category Name of the thermodynamic quantity.
  !! @param[in] total If applicable, total value for the system.  If not, use
  !!                  HUGE(1d0).
  !! @param[in] item Array of decomposition values.  E.g., site
  !!                 contributions. If any of these values are huge the entire
  !!                 output line is supressed.
  !! @param[in] nitem Number of items. Needed for the NAB case where we don't
  !!                  use a module
  !! @sideeffect Writes to outunit.
  subroutine thermo_print_results_line(category, total, item, nitem)
    use amber_rism_interface
    implicit none
    character(len=*), intent(in) :: category
    integer, intent(in) :: nitem
    _REAL_, intent(in) :: total, item(nitem)
    ! Format for category (calculation type, e.g. excess chemical
    ! potential).
    character(len=64) :: catFmt = "(a30)"
    ! Format for a string the same width as valfmt.
    character(len=32) :: strFmt = "(a16)"
    ! Format for floating point values.
    character(len=32) :: valFmt = '(1p, 2x, e14.6e3)'
    ! Long string of whitespace that can be used to effect a
    ! left-justified string.  Otherwise strings are right-justified.
    ! Simply concatenate this to the end of the string you wish
    ! left-justified.
    character(len=64) :: whtspc
    integer :: i
    write(whtspc, '(a64)')" "

    if (any(item == HUGE(1d0)) .or. (total == HUGE(1d0) .and. nitem == 0)) return
    write(outunit, catFmt, advance='no') trim(category) // whtspc
    if (total /= HUGE(1d0)) then
       write(outunit, valFmt, advance='no') total
    else
       write(outunit, strFmt, advance='no') ""
    end if
    do i = 1, nitem
       write(outunit, valFmt, advance='no') item(i)
    end do
    write(outunit, '(a)')
  end subroutine thermo_print_results_line

  
  !> Finalizes all of the 3D-RISM objects and frees memory.
  subroutine rism_finalize()
    use amber_rism_interface
    use safemem
    use rism3d_solvent_c
    use rism3d_solute_c
    implicit none
    integer :: err
    integer*8 :: memstats(10)
    if (rismprm%rism == 1) then
       call rism3d_destroy(rism_3d)
       call rismthermo_destroy(rismthermo)
       call rismthermo_destroy(rismthermo_pol)
       call rismthermo_destroy(rismthermo_apol)
       call rism3d_solvent_destroy(solvent)
       call rism3d_solute_destroy(solute)
       if (safemem_dealloc(ff)/= 0) &
            call rism_report_error("Deallocation in Amber-RISM interface failed")
       if (safemem_dealloc(closurelist) /= 0) &
            call rism_report_error("Deallocation in Amber-RISM interface failed")

       if (safemem_dealloc(tolerancelist) /= 0) &
            call rism_report_error("Deallocation in Amber-RISM interface failed")
       call rism_max_memory()
    end if
  end subroutine rism_finalize

  
  !> Prints the maximum amount of memory allocated at any one time so
  !! far in the run.
  subroutine rism_max_memory()
    use amber_rism_interface
    use safemem
#ifdef MPI
    use mpi
#endif
    implicit none
    integer*8 :: memstats(10), tmemstats(10)
    integer :: err, irank
    outunit = rism_report_getMUnit()
    memstats = memStatus()
#ifdef MPI
    if (mpirank==0) then
       call MPI_REDUCE(MPI_IN_PLACE, memstats, ubound(memstats, 1), MPI_INTEGER8, &
            MPI_SUM, 0, mpicomm, err)
    else
       call MPI_REDUCE(memstats, memstats, ubound(memstats, 1), MPI_INTEGER8, &
            MPI_SUM, 0, mpicomm, err)
    end if
    if (err/= 0) call rism_report_warn("RISM_MAX_MEMORY: MPI_REDUCE failed.")
#endif
    if (mpirank==0) then
       write(outunit, '(a)')
       write(outunit, '(a)') "|3D-RISM memory allocation summary"
       write(outunit, '(a)') "|Type          Maximum        Current   "
       write(outunit, '(a, 2(f12.5, a))') "|Integer  ", &
            dble(memstats(6))/BYTES_PER_GB, " GB", &
            dble(memstats(1))/BYTES_PER_GB, " GB"
       write(outunit, '(a, 2(f12.5, a))') "|Real     ", &
            dble(memstats(7))/BYTES_PER_GB, " GB", &
            dble(memstats(2))/BYTES_PER_GB, " GB"
       write(outunit, '(a, 2(f12.5, a))') "|Logical  ", &
            dble(memstats(8))/BYTES_PER_GB, " GB", &
            dble(memstats(3))/BYTES_PER_GB, " GB"
       write(outunit, '(a, 2(f12.5, a))') "|Character", &
            dble(memstats(9))/BYTES_PER_GB, " GB", &
            dble(memstats(4))/BYTES_PER_GB, " GB"
       write(outunit, '(a)') "|---------------------------------------"
       write(outunit, '(a, 2(f12.5, a))') "|Total    ", &
            dble(memstats(10))/BYTES_PER_GB, " GB", &
            dble(memstats(5))/BYTES_PER_GB, " GB"
    end if
  end subroutine rism_max_memory

  
!!! I/O: performs RISM related I/O for files that only deal with RISM data

  !> Outputs volumetric data, such as solvent and electric potential
  !! distributions, to their respective files. Each distribution
  !! is written in a separate file with the step number before the
  !! suffix.  File names are taken from guvfile, huvfile, cuvfile, and
  !! similar variables.
  !! @param[in,out] this 3D-RISM object.
  !! @param[in] step Step number used as a suffix.
  subroutine rism_writeVolumetricData(this, step)
    use constants_rism, only : COULOMB_CONST_E, KB, PI
    use amber_rism_interface
    use rism_io
    use rism3d_mrc
    use rism3d_opendx
    use rism3d_xyzv
    use safemem
    use rism_util, only : checksum
#ifdef MPI
    use mpi
#endif
    implicit none
    
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: step

    integer :: iv, ivv, nsolv, igx, igy, igz, ios, i, j, k, ig
    character(len=16) :: cstep
    character(len=64) :: suffix
    character(len=6) :: extension
  
#ifdef MPI
    integer :: err
#endif /*MPI*/

    procedure (writeVolumeInterface), pointer :: writeVolume => NULL()

    if (volfmt .eq. 'mrc') then
       extension = '.mrc'
       writeVolume => rism3d_mrc_map_write
    else if (volfmt .eq. 'dx') then
       extension = '.dx'
       writeVolume => rism3d_opendx_write
    else if (volfmt .eq. 'xyzv') then
       extension = '.xyzv'
       writeVolume => rism3d_xyzv_write
    else
       call rism_report_warn("Volume format "//trim(volfmt) &
            //" is not one of the currently handled volumetric data formats." &
            //" No volume data will be written.")
       return
    end if

    call rism_writePdfTcfDcf(this, writeVolume, step, extension)
    call rism_writeThermo(this, writeVolume, step, extension)
    
  end subroutine rism_writeVolumetricData

  !> Outputs PDF, TCF, and DCF. Each distribution
  !! is written in a separate file with the step number before the
  !! suffix. 
  !! @param[in,out] this 3D-RISM object.
  !! @param[in] step Step number used as a suffix.
  subroutine rism_writePdfTcfDcf(this, writeVolume, step, extension)
    use constants_rism, only : COULOMB_CONST_E, KB, PI
    use amber_rism_interface
    use rism_io
    use rism3d_mrc
    use rism3d_opendx
    use rism3d_xyzv
    use safemem
    implicit none

    type(rism3d), intent(inout) :: this
    character(len=*), intent(in) :: extension
    integer, intent(in) :: step
    procedure (writeVolumeInterface) :: writeVolume 
    
    character(len=64) :: suffix
    character(len=16) :: cstep

    _REAL_, pointer :: work(:, :, :) => NULL()
    integer :: i, j, k, ig, iv

    
    if (len_trim(guvfile) /= 0 .or. len_trim(huvfile) /= 0) then
       work => safemem_realloc(work, this%grid%globalDimsR(1), this%grid%globalDimsR(2), this%grid%localDimsR(3))
    end if
    do iv = 1, this%solvent%numAtomTypes
       write(cstep, '(i16)') step
       cstep = adjustl(cstep)
       suffix = '.'//trim(rism_3d%solvent%atomName(iv))//'.'//trim(cstep)
       suffix = trim(suffix)//extension
       
       ! Solute-solvent RDF.
       if (len_trim(guvfile) /= 0)  then
#  if defined(MPI)
          do k = 1, this%grid%localDimsR(3)
             do j = 1, this%grid%globalDimsR(2)
                do i = 1, this%grid%globalDimsR(1)
                   work(i, j, k) = &
                        this%guv(i + (j - 1) * (this%grid%globalDimsR(1) + 2) &
                        + (k - 1) * (this%grid%globalDimsR(1) + 2) &
                        * this%grid%globalDimsR(2), iv)
                end do
             end do
          end do
          call writeVolume(trim(guvfile)//suffix, work, this%grid, &
               this%solute, mpirank, mpisize, mpicomm)
#  else
          call writeVolume(trim(guvfile)//suffix, this%guv(:, iv), &
               this%grid, this%solute)

#  endif /*defined(MPI)*/
       endif

       ! Solute-solvent TCF.
       if (len_trim(huvfile) /= 0)  then
#  if defined(MPI)
          do k = 1, this%grid%localDimsR(3)
             do j = 1, this%grid%globalDimsR(2)
                do i = 1, this%grid%globalDimsR(1)
                   work(i, j, k) = &
                        this%huv(i + (j - 1) * (this%grid%globalDimsR(1) + 2) &
                        + (k - 1) * (this%grid%globalDimsR(1) + 2) &
                        * this%grid%globalDimsR(2), iv)
                end do
             end do
          end do
          call writeVolume(trim(huvfile)//suffix, work, this%grid, &
               this%solute, mpirank, mpisize, mpicomm)
#  else
          call writeVolume(trim(huvfile)//suffix, this%huv(:, iv), this%grid, this%solute)
#  endif /*defined(MPI)*/
       endif

       ! Solute-solvent DCF.
       if (len_trim(cuvfile) /= 0)  then
          do k = 1, this%grid%localDimsR(3)
             do j = 1, this%grid%globalDimsR(2)
                do i = 1, this%grid%globalDimsR(1)
                   ig = i + (j - 1) * this%grid%localDimsR(1) + &
                        (k - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                   work(i, j, k) = this%cuv(i,j,k, iv)
                end do
             end do
          end do
          call writeVolume(trim(cuvfile)//suffix, work, &
               this%grid, this%solute, mpirank, mpisize, mpicomm)
          ! call writeVolume(trim(cuvfile)//suffix, this%cuv(:, :, :, iv), &
          !      this%grid, this%solute, mpirank, mpisize, mpicomm)

       endif
    end do
    if (safemem_dealloc(work)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate WORK")
  end subroutine rism_writePdfTcfDcf

  !> Outputs thermodynamics distributions. Each distribution
  !! is written in a separate file with the step number before the
  !! suffix. 
  !! @param[in,out] this 3D-RISM object.
  !! @param[in] step Step number used as a suffix.
  subroutine rism_writeThermo(this, writeVolume, step, extension)
    use constants_rism, only : COULOMB_CONST_E, KB, PI
    use amber_rism_interface
    use rism_io
    use rism3d_mrc
    use rism3d_opendx
    use rism3d_xyzv
    use safemem
    implicit none

    type(rism3d), intent(inout) :: this
    character(len=*), intent(in) :: extension
    integer, intent(in) :: step
    procedure (writeVolumeInterface) :: writeVolume 
    
    character(len=64) :: suffix
    character(len=16) :: cstep

    integer :: i, j, k, ig, iv

    _REAL_, pointer :: work(:, :, :) => NULL()
    _REAL_, pointer :: excessChemicalPotential_map(:, :, :) => NULL(), &
         solvationEnergy_map(:, :, :) => NULL(), &
         solventPotentialEnergy_map(:, :, :) => NULL()
    _REAL_, pointer :: excessChemicalPotential_V_map(:, :, :, :) => NULL(), &
         solvationEnergy_V_map(:, :, :, :) => NULL(), &
         solventPotentialEnergy_V_map(:, :, :, :) => NULL(),&
         solventEntropy_V_map(:, :, :, :) => NULL()

    ! Outputting charge and thermodynamic distributions.
    if (len_trim(quvfile) /= 0 .or. len_trim(chgDistFile) /= 0 .or. &
         len_trim(excessChemicalPotentialfile) /= 0 .or. &
         len_trim(solvationEnergyfile)/= 0 .or. &
         len_trim(entropyfile) /= 0 .or. &
         len_trim(solventPotentialEnergyfile) /= 0 ) then
       write(cstep, '(i16)') step
       cstep = adjustl(cstep)
       suffix = '.'//trim(cstep)
       suffix = trim(suffix)//extension
       work => safemem_realloc(work, this%grid%globalDimsR(1), this%grid%globalDimsR(2), this%grid%localDimsR(3))
       work = 0
       ! Sum the contributions from each solvent type at each grid
       ! point and convert units to [e/A^3].
       do iv = 1, this%solvent%numAtomTypes
#  if defined(MPI)
          do k = 1, this%grid%localDimsR(3)
             do j = 1, this%grid%globalDimsR(2)
                do i = 1, this%grid%globalDimsR(1)
                   work(i, j, k) = work(i, j, k) &
                        + this%guv(i + (j - 1) * (this%grid%globalDimsR(1) + 2) &
                        & + (k - 1) * (this%grid%globalDimsR(1) + 2) * this%grid%globalDimsR(2), iv) &
                        * sqrt((KB *this%solvent%temperature) / COULOMB_CONST_E) &
                        * this%solvent%charge(iv) * this%solvent%density(iv)
                end do
             end do
          end do
#  else
          call DAXPY(this%grid%totalLocalPointsR, sqrt((KB *this%solvent%temperature)/COULOMB_CONST_E) &
               *this%solvent%charge(iv)*this%solvent%density(iv), this%guv(1, iv), 1, work, 1)
#  endif /*defined(MPI)*/
       end do
       if (len_trim(quvfile) /= 0) then
          call writeVolume(trim(quvfile)//suffix, work, this%grid, &
               this%solute, mpirank, mpisize, mpicomm)
       end if
       if (len_trim(chgDistFile) /= 0)  then
          call DSCAL(this%grid%totalLocalPointsR, this%grid%voxelVolume, work, 1)
          call writeVolume(trim(chgDistfile)//suffix, work, this%grid, &
               this%solute, mpirank, mpisize, mpicomm)
       end if
       
    endif
  contains
    subroutine writeThermo(wv, root, suffix, data, data_sites)
      implicit none
      ! we are passing writeVolume to work around what appears to be a
      ! PGI compiler bug.  writeThermo has scope to see writeVolume so
      ! we should be able to use it directly.  Doing so works with GNU
      ! and Intel but not PGI
      procedure(writeVolumeInterface) :: wv
      character(len=*), intent(in) :: root, suffix
      _REAL_, intent(in) :: data(:,:,:), data_sites(:,:,:,:)
      call wv(trim(root)//trim(suffix), data, this%grid, &
           this%solute, mpirank, mpisize, mpicomm)
      do iv = 1, this%solvent%numAtomTypes
         call wv(trim(root)//'.'//trim(rism_3d%solvent%atomName(iv))//trim(suffix), &
              data_sites(:, :, :, iv), this%grid, &
              this%solute, mpirank, mpisize, mpicomm)
      end do
    end subroutine writeThermo
   end subroutine rism_writeThermo
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MPI
  !> Broadcasts initialization information about the system from the
  !! master to all the other processes.
  !!
  !! @param[in] mrank MPI process rank
  !! @param[in] msize Number of MPI processes
  !! @param[in] mcomm MPI communicator
  !!
  !! @sideeffects sets MPI parmeters and distributes run information
  !!              to all processes
  subroutine rism_mpi_bcast(mrank, msize, mcomm)
    use amber_rism_interface
    use mpi
    implicit none
    ! add the 'm' to avoid clashing with intrinsics, like size()
    integer, intent(in) :: mrank !< MPI process rank.
    integer, intent(in) :: msize !< Number of MPI processes.
    integer, intent(in) :: mcomm !< MPI communicator index.
    integer :: err
    integer :: nclosure
    ! Private subroutine so there should be no timer.

    mpirank = mrank
    mpisize = msize
    mpicomm = mcomm

    ! Could be done by creating and passing an MPI derived type, but
    ! this is done once so it is not worth the effort.
    call mpi_bcast(rismprm%rism, 1, mpi_integer, 0, mpicomm, err)
    if (rismprm%rism == 1) then
       ! Broadcast the entire rismprm object. Not everything actually needs
       ! to be transferred, but there are not that many exceptions.
       call mpi_bcast(rismprm%solvcut, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast SOLVCUT")
       call mpi_bcast(rismprm%grdspc, 3, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast GRDSPC")
       call mpi_bcast(rismprm%mdiis_del, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast MDIIS_DEL")
       call mpi_bcast(rismprm%mdiis_restart, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast MDIIS_RESTART")
       call mpi_bcast(rismprm%chargeSmear, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast chargeSmear")

       call mpi_bcast(rismprm%ng3, 3, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast NG3")
       call mpi_bcast(rismprm%maxstep, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast MAXSTEP")
       call mpi_bcast(rismprm%npropagate, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast NPROPAGATE")
       call mpi_bcast(rismprm%zerofrc, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast ZEROFRC")
       call mpi_bcast(rismprm%apply_rism_force, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast APPLY_RISM_FORCE")
       call mpi_bcast(rismprm%rismnrespa, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast RISMNRESPA")
       call mpi_bcast(rismprm%saveprogress, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast SAVEPROGRESS")
       call mpi_bcast(rismprm%ntwrism, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast NTWRISM")
       call mpi_bcast(rismprm%verbose, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast VERBOSE")
       call mpi_bcast(rismprm%progress, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast PROGRESS")       

       if (mpirank==0) &
            nclosure=ubound(closurelist, 1)
       call mpi_bcast(nclosure, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast NCLOSURE")
       if (mpirank/= 0) then
          closurelist=>safemem_realloc(closurelist, len(closurelist), nclosure)
          tolerancelist=>safemem_realloc(tolerancelist, nclosure)
       end if
       call mpi_bcast(closurelist, len(closurelist)*nclosure, mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast CLOSURE")
       call mpi_bcast(tolerancelist, nclosure, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast TOLERANCE")

       call mpi_bcast(rismprm%write_thermo, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast WRITE_THERMO")

       ! I/O
       ! Special output files that all nodes write to.
       call mpi_bcast(guvfile, len(guvfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast GUVFILE")
       call mpi_bcast(huvfile, len(huvfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast HUVFILE")
       call mpi_bcast(cuvfile, len(cuvfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast CUVFILE")
       call mpi_bcast(uuvfile, len(uuvfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast UUVFILE")
       call mpi_bcast(quvfile, len(quvfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast QUVFILE")
       call mpi_bcast(chgdistfile, len(chgdistfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast CHGDISTFILE")
       call mpi_bcast(excessChemicalPotentialfile, len(excessChemicalPotentialfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast EXCESSCHEMICALPOTENTIALFILE")
       call mpi_bcast(solvationEnergyfile, len(solvationEnergyfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast SOLVATIONENERGYFILE")
       call mpi_bcast(entropyfile, len(entropyfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast ENTROPYFILE")
       call mpi_bcast(solventPotentialEnergyfile, len(solventPotentialEnergyfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast POTUVFILE")
       call mpi_bcast(volfmt, len(volfmt), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast VOLFMT")
       call mpi_bcast(periodicPotential, len(periodicPotential), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast PERIODICPOTENTIAL")
       call mpi_bcast(unitCellDimensions, 6, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast UNITCELLDIMENSIONS")
    end if
  end subroutine rism_mpi_bcast
#endif /*MPI*/

  !> Sets default values for 3D-RISM paramters.
  subroutine defaults()
    use amber_rism_interface
    use constants_rism, only: NO_INPUT_VALUE, NO_INPUT_VALUE_FLOAT
    implicit none

    closurelist => safemem_realloc(closurelist, len(closurelist), nclosuredefault)
    closurelist(2:)           = ''
    closurelist(1)            = 'KH'
    rismprm%closureOrder      = 1
    periodicPotential         = 'pme'

    !solvation box
#ifdef API
    rismprm%solvcut         = 8.d0
    rismprm%grdspc          = 0.8d0
#else
    rismprm%solvcut         = 9.d0
    rismprm%grdspc          = 0.5d0
#endif
    rismprm%ng3             = -1

    !convergence
    tolerancelist => safemem_realloc(tolerancelist, nclosuredefault)
    tolerancelist             = HUGE(1d0)
    tolerancelist(1)          = 1d-7
    rismprm%mdiis_del         = 0.4d0
    rismprm%mdiis_nvec        = 5
    rismprm%mdiis_method      = 2
    rismprm%mdiis_restart     = 10d0
    rismprm%maxstep           = 10000
    rismprm%npropagate        = 5

    !imin = 1 (minimization)
    rismprm%zerofrc          = 1

    !imin = 5 (trajectory analysis)
    rismprm%apply_rism_force = 1

    !imin = 0 (MD)
    rismprm%rismnrespa       = 1

    !output
    rismprm%saveprogress     = 0
    rismprm%ntwrism          = -1
#ifdef API
    rismprm%verbose          = -100
#else
    rismprm%verbose          = 0
#endif
    rismprm%progress         = 1
    volfmt                   = 'mrc'

    !charge smear
    rismprm%chargeSmear = 1d0
    rismprm%write_thermo=1
  end subroutine defaults
  
  !> Reads the RISM namelist from the input file.
  subroutine read_namelist(mdin_unit)
    use amber_rism_interface
    implicit none
    !mdin_unit :: unit number for mdin file
    character(len=8) :: closure(nclosuredefault)
    _REAL_ :: tolerance(nclosuredefault)
    integer :: closureOrder
    integer, intent(in) :: mdin_unit
    character(len=8) :: periodic
    _REAL_ :: solvcut
    _REAL_ :: grdspc(3)
    integer ::  ng3(3)
    _REAL_ :: mdiis_del
    integer :: mdiis_nvec
    integer :: mdiis_method
    _REAL_ :: mdiis_restart
    integer :: maxstep
    integer :: npropagate
    integer :: zerofrc
    integer :: apply_rism_force
    integer :: rismnrespa
    integer :: saveprogress
    integer :: ntwrism
    integer :: verbose
    integer :: progress
    _REAL_ :: chargeSmear
    integer :: write_thermo
    namelist /rism/ &
         closure, closureOrder, periodic, &
         grdspc, solvcut, ng3, &
         tolerance, mdiis_del, mdiis_nvec, mdiis_method, &
         mdiis_restart, maxstep, npropagate, zerofrc, &
         apply_rism_force, &
         rismnrespa, chargeSmear, write_thermo, &
         saveprogress, ntwrism, verbose, progress, volfmt
    
    call flush(0)

    closure = closurelist
    tolerance = tolerancelist
    closureOrder = rismprm%closureOrder
    periodic = periodicPotential
    solvcut = rismprm%solvcut
    grdspc= rismprm%grdspc
    ng3 = rismprm%ng3
    mdiis_del = rismprm%mdiis_del
    mdiis_nvec = rismprm%mdiis_nvec
    mdiis_method = rismprm%mdiis_method
    mdiis_restart = rismprm%mdiis_restart
    maxstep = rismprm%maxstep
    npropagate = rismprm%npropagate
    zerofrc = rismprm%zerofrc
    apply_rism_force = rismprm%apply_rism_force
    rismnrespa = rismprm%rismnrespa
    saveprogress = rismprm%saveprogress
    ntwrism = rismprm%ntwrism
    verbose= rismprm%verbose
    progress = rismprm%progress
    chargeSmear = rismprm%chargeSmear
    write_thermo = rismprm%write_thermo

!!$  !resize tolerance to the size of closure
!!$  tolerancelist => safemem_realloc(tolerancelist, size(closurelist))
!!$  tolerancelist(2:)=HUGE(1d0)

    rewind(mdin_unit)
    read(mdin_unit, nml=rism)

    closurelist = closure
    tolerancelist = tolerance
    rismprm%closureOrder = closureOrder
    periodicPotential = periodic
    ! Solvation box.
    rismprm%grdspc=grdspc
    rismprm%solvcut=solvcut
    rismprm%ng3=ng3
    ! Convergence.
    rismprm%mdiis_del=mdiis_del
    rismprm%mdiis_nvec=mdiis_nvec
    rismprm%mdiis_method=mdiis_method
    rismprm%mdiis_restart=mdiis_restart
    rismprm%maxstep=maxstep
    rismprm%npropagate=npropagate
    ! Minimization.
    rismprm%zerofrc=zerofrc
    ! imin=5
    rismprm%apply_rism_force=apply_rism_force
    ! md
    rismprm%rismnrespa=rismnrespa
    ! Output.
    rismprm%saveprogress=saveprogress
    rismprm%ntwrism=ntwrism
    rismprm%verbose=verbose
    rismprm%progress=progress
    rismprm%chargeSmear = chargeSmear
    rismprm%write_thermo = write_thermo

    ! Set the RISM cutoff if not set by the user.
    if (rismprm%solvcut < 0) then
       rismprm%solvcut = 9.0
    end if

  end subroutine read_namelist

  !> Checks user input to ensure that it is not completely crazy.
  subroutine sanity_check()
    use amber_rism_interface
    use rism_util, only : caseup
    use array_util, only : array_index
    implicit none
    character(len=32) :: fmt
    !iclosure :: counter for closures
    integer :: iclosure

    ! Ensure that an apropriate file format has been chosen for
    ! volumetric output.
    if (.not. (volfmt .eq. "mrc" .or. volfmt .eq. "dx" .or. volfmt .eq. "xyzv")) then
       call rism_report_error("Only 'mrc', 'dx', and 'xyzv' volumetric data formats are supported")
    end if

    ! Resize closure list to the appropriate size.
    if (len_trim(closurelist(size(closurelist))) == 0) then
       closurelist => safemem_realloc(closurelist, len(closurelist), array_index(closurelist, '')-1)
    end if

    ! Resize and setup the tolerance list.
    ! 1) Get rid of extraneous values.
    if (array_index(tolerancelist, HUGE(1d0))>0) then
       tolerancelist=>safemem_realloc(tolerancelist, &
            array_index(tolerancelist, HUGE(1d0))-1)
    end if

    ! 2) If there is only one closure, use only the last tolerance.
    if (size(closurelist) == 1) then
       tolerancelist(1) = tolerancelist(size(tolerancelist))
       tolerancelist=>safemem_realloc(tolerancelist, 1)
       ! 3) If there is one tolerance, the default for intermediate closures is 1.
    else if (size(tolerancelist) == 1) then
       tolerancelist=>safemem_realloc(tolerancelist, size(closurelist))
       tolerancelist(size(tolerancelist)) = tolerancelist(1)
       tolerancelist(1:size(tolerancelist)-1) = 1d0
       ! 4) If there are two tolerances, the first is for intermediate closures.
    else if (size(tolerancelist) == 2) then
       tolerancelist=>safemem_realloc(tolerancelist, size(closurelist))
       tolerancelist(size(tolerancelist)) = tolerancelist(2)
       tolerancelist(1:size(tolerancelist)-1) = tolerancelist(1)
       ! 5) Otherwise, there should be the same number of tolerances and closures.
    else if (size(tolerancelist) /= size(closurelist)) then
       call rism_report_error("number of tolerances must be 1, 2 or the number of closures.")
    end if

    ! If a closure number is given, map it to a name.
    do iclosure = 1, ubound(closurelist, 1)
       if (trim(closurelist(iclosure)) .eq. "0") then
          closurelist(iclosure) = "HNC"
       else if (trim(closurelist(iclosure)) .eq. "1") then
          closurelist(iclosure) = "KH"
       else if (trim(closurelist(iclosure)) .eq. "2") then
          closurelist(iclosure) = "PSEN"
       end if
       ! If the old method of indicating the PSE order has been used
       ! (closureOrder) then reformat to the new method.
       call caseup(closurelist(iclosure))
       if (trim(closurelist(iclosure)).eq."PSEN" .or. trim(closurelist(iclosure)).eq."PSE") then
          ! Check if 'closureOrder' is used with a list of closures.
          if (iclosure > 1) &
               call rism_report_error("'closureOrder' is depricated and not compatible "&
               //"with closure lists")
          write(fmt, '(a, i4, a)') '(a, i', int(log10(dble(rismprm%closureOrder))) + 1, ')'
          write(closurelist, fmt) "PSE", rismprm%closureOrder
       end if
    end do

    ! charge smear
    if (rismprm%chargeSmear .lt. 0d0) then
       call rism_report_error('(a,g)',"'chargeSmear' must be >= 0. Got",rismprm%chargeSmear)
    end if
  end subroutine sanity_check

  !> Writes the contents of a C string (array of chars) to a Fortran
  !! string.  Will not write past the end of the Fortran string.
  !! @param[in] fstr Fortran string to write to.
  !! @param[in] cstr C string to write from.
  !! @param[in] nchar Number of chars in cstr (not including null).
  subroutine cstr2fstr(fstr, cstr, nchar)
    implicit none
    character(len=*), intent(out) :: fstr
    integer, intent(in) :: nchar
    integer(kind=1), intent(in) :: cstr(nchar + 1)
    integer :: i
    fstr = ""
    do i = 1, min(nchar, len(fstr))
       if (cstr(i) == 0) exit
       fstr = trim(fstr)//char(cstr(i), 1)
    end do
  end subroutine cstr2fstr


!!!Calculates the least common multiple of integers a and b.
!!!IN:
!!!  a : integer
!!!  b : integer
!!!OUT:
!!! least common multiple
  function mylcm(a, b)
    use rism_util, only: lcm
    implicit none
    integer :: mylcm, a, b
    mylcm = lcm(a, b)
  end function mylcm

#if 0
! Initializes a struct with RISM options to all default values
!
! Parameters
! ----------
! inp : type(rismprm_t)
!     struct of RISM input options that will be filled by this subroutine
subroutine rism_sander_input(inp)

   implicit none
   type(rismprm_t), intent(out) :: inp
   call defaults

end subroutine rism_sander_input
#endif

end module sander_rism_interface

#ifdef OPENMP
  subroutine set_omp_num_threads_rism()
    use constants_rism, only: omp_num_threads
    implicit none
    character(len=5) :: omp_threads
    integer :: ier

    call get_environment_variable('OMP_NUM_THREADS', omp_threads, status=ier)
    if( ier .ne. 1 ) read( omp_threads, * ) omp_num_threads
#ifndef API
    write(6,'(a,i3,a)') '| Running RISM OpenMP with ',omp_num_threads,' threads'
#endif
  end subroutine set_omp_num_threads_rism
#endif

subroutine rism_defaults()
   use sander_rism_interface, only: defaults
   implicit none
   call defaults()
   return
end subroutine rism_defaults

#if 0
!  Keep outside of the module, to avoid name mangling
#ifdef API
  subroutine rism_setparam2( solvcut, grdspc, verbose )
     use amber_rism_interface
     use constants_rism, only: NO_INPUT_VALUE_FLOAT, NO_INPUT_VALUE
     implicit none
     double precision, intent(in):: solvcut, grdspc
     integer, intent(in) :: verbose

     if( solvcut <= 0.d0 ) then
        rismprm%solvcut = 8.d0  !default
     else
        rismprm%solvcut = solvcut
     endif
     if( grdspc <= 0.d0 ) then
        rismprm%grdspc(:) = 0.5d0  !default
     else
        rismprm%grdspc(:) = grdspc
     endif
     if( verbose == -100 ) then
        rismprm%verbose = -1  !default
     else
        rismprm%verbose = verbose
     endif

     return
  end subroutine rism_setparam2
#endif
#endif
