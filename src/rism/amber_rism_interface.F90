!<compile=optimized>

! The 3D-RISM-KH software found here is copyright (c) 2010-2012 by
! Andriy Kovalenko, Tyler Luchko and David A. Case.
!
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
!
! You should have received a copy of the GNU General Public License in the
! ../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
! Users of the 3D-RISM capability found here are requested to acknowledge
! use of the software in reports and publications.  Such acknowledgement
! should include the following citations:
!
! 1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
! ibid. 112:10391-10417 (2000).
!
! 2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by
! F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.
!
! 3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
! and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010).
!
! 4) I. Omelyan and A. Kovalenko, Mol. Simul. 39:25-48 (2013).

#include "../include/dprec.fh"

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
  use rism_timer_c

  !> Derived type to store thermodynamic results.  To facilitate MPI
  !! communication, many values are stored in a common array.
  !! Pointers, however, are used to access the data.
  type rismthermo_t

     ! Thermodynamic values that use mpi_buffer for memory.

     !> Excess chemical potential.
     _REAL_, pointer :: excessChemicalPotential(:) => NULL()
     !> Gaussian fluctuation excess chemical potential.
     _REAL_, pointer :: excessChemicalPotentialGF(:) => NULL()
     !> Initial State Correction excess chemical potential
     _REAL_, pointer :: excessChemicalPotentialPCPLUS => NULL()
     !> Universal Correction excess chemical potential.
     _REAL_, pointer :: excessChemicalPotentialUC => NULL()
     !> Solute-solvent interaction energy.
     _REAL_, pointer :: solventPotentialEnergy(:) => NULL()
     !> Solvation energy.
     _REAL_, pointer :: solvationEnergy(:) => NULL()
     !> Gaussian fluctuation solvation energy.
     _REAL_, pointer :: solvationEnergyGF(:) => NULL()
     !> Initial State Correction solvation energy
     _REAL_, pointer :: solvationEnergyPCPLUS => NULL()
     !> Palmer correction solvation energy.
     _REAL_, pointer :: solvationEnergyUC => NULL()
     !> Partial molar volume.
     _REAL_, pointer :: partialMolarVolume => NULL()
     !> Partial molar volume temperature derivative.
     _REAL_, pointer :: partialMolarVolume_dT => NULL()
     !> Total number of particles.
     _REAL_, pointer :: totalParticlesBox(:) => NULL()
     !> Excess number of particles.
     _REAL_, pointer :: excessParticlesBox(:) => NULL()
     !> Excess number of particles after asymptotic TCF correction.
     _REAL_, pointer :: excessParticles(:) => NULL()
     !> Temperature derivative of the excess number of particles.
     _REAL_, pointer :: excessParticles_dT(:) => NULL()
     !> Kirkwood-Buff integral (a.k.a. total correlation function
     !! integral).
     _REAL_, pointer :: kirkwoodBuff(:) => NULL()
     !> Direct correlation function integral.
     _REAL_, pointer :: DCFintegral(:) => NULL()
     !> Kirkwood-Buff integral temperature derivative.
     _REAL_, pointer :: kirkwoodBuff_dT(:) => NULL()
     !> Temperature derivative of the direct correlation function integral.
     _REAL_, pointer :: DCFintegral_dT(:) => NULL()

     !> TODO: Document me.
     integer :: mpirank = 0, mpisize = 1, mpicomm = 0

     !> All the results are stored in a single array so
     !! the necessary reductions are done in a single
     !! communication.  This is also used for the serial
     !! calculation for simplicity.
     !! Order and size:
     !!         excessChemicalPotential(solvent%numAtomTypes), excessChemicalPotentialGF(solvent%numAtomTypes),
     !!         solventPotentialEnergy(solvent%numAtomTypes),
     !!         solvationEnergy(solvent%numAtomTypes),solvationEnergyGF(solvent%numAtomTypes),
     !!         excessNum(solvent%numAtomTypes), excessNum_dT(solvent%numAtomTypes),
     !!         kirkwoodBuff(solvent%numAtomTypes),
     !!         kirkwoodBuff_dT(solvent%numAtomTypes), dcfi(solvent%numAtomTypes),
     !!         DCFintegral_dT(solvent%numAtomTypes),
     !!         excessChemicalPotentialPCPLUS(1), solvationEnergyPCPLUS(1),
     !!         excessChemicalPotentialUC(1), solvationEnergyUC(1),
     !!         partialMolarVolume(1), partialMolarVolume_dT(1)
     _REAL_, private, pointer :: mpi_buffer(:) => NULL()
#if ! defined( USE_MPI_IN_PLACE) && defined( MPI )
     _REAL_, private, pointer :: tmpi_buffer(:) => NULL()
#endif
  end type rismthermo_t
  
contains

  !> Initializes and allocates memory for rismthermo objects.
  !! @param[in,out] this The rismthermo object.
  !! @param[in] nsite Number of solvent sites.
  !! @param[in] o_mpicomm  MPI communicator (optional).
  !! @sideeffects Allocates memory and sets value to huge.
  subroutine rismthermo_new(this, nsite, o_mpicomm)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
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
    this%mpi_buffer => safemem_realloc(this%mpi_buffer, 13*nsite + 6)
    ! Initialize for case that some values are not calculated.
    call rismthermo_reset(this)
#if !defined( USE_MPI_IN_PLACE) && defined( MPI )
    this%tmpi_buffer => &
         safemem_realloc(this%tmpi_buffer, ubound(this%mpi_buffer, 1))
#endif
    this%excessChemicalPotential => this%mpi_buffer(1:nsite)
    this%excessChemicalPotentialGF => this%mpi_buffer(nsite + 1:2*nsite)
    this%solventPotentialEnergy => this%mpi_buffer(2*nsite + 1:3*nsite)
    this%solvationEnergy => this%mpi_buffer(3*nsite + 1:4*nsite)
    this%solvationEnergyGF => this%mpi_buffer(4*nsite + 1:5*nsite)
    this%excessParticlesBox => this%mpi_buffer(5*nsite + 1:6*nsite)
    this%excessParticles_dT => this%mpi_buffer(6*nsite + 1:7*nsite)
    this%kirkwoodBuff => this%mpi_buffer(7*nsite + 1:8*nsite)
    this%kirkwoodBuff_dT => this%mpi_buffer(8*nsite + 1:9*nsite)
    this%DCFintegral => this%mpi_buffer(9*nsite + 1:10*nsite)
    this%DCFintegral_dT => this%mpi_buffer(10*nsite + 1:11*nsite)
    this%totalParticlesBox => this%mpi_buffer(11*nsite + 1:12*nsite)
    this%excessParticles => this%mpi_buffer(12*nsite + 1:13*nsite)
    this%excessChemicalPotentialPCPLUS => this%mpi_buffer(13*nsite + 1)
    this%solvationEnergyPCPLUS => this%mpi_buffer(13*nsite + 2)
    this%excessChemicalPotentialUC => this%mpi_buffer(13*nsite + 3)
    this%solvationEnergyUC => this%mpi_buffer(13*nsite + 4)
    this%partialMolarVolume => this%mpi_buffer(13*nsite + 5)
    this%partialMolarVolume_dT => this%mpi_buffer(13*nsite + 6)
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

  
!!!Subtracts all of the values from rismthermo object C from B and
!!!puts the result in A
!!!
!!!  A = B - C
!!!
!!!Any values that are HUGE in either B or C will be HUGE in A.
!!!
!!!This is useful for calculating polar contributions from full and
!!!apolar calculations, as an example.
!!!
!!!IN:
!!! A : the rismthermo object in which to place the result
!!! B : the rismthermo object substracted from
!!! C : the rismthermo object subtracted
  subroutine rismthermo_sub(A, B, C)
    implicit none
    type(rismthermo_t), intent(inout) :: A
    type(rismthermo_t), intent(in) :: B, C
    where(B%mpi_buffer /= HUGE(1d0) .and. C%mpi_buffer /= HUGE(1D0))
       A%mpi_buffer = B%mpi_buffer - C%mpi_buffer
    elsewhere
       A%mpi_buffer = huge(1d0)
    end where
  end subroutine rismthermo_sub

  
!!!Does an MPI reduction on all thermodynamic values
!!!IN:
!!! this : the rismthermo object
!!!SIDEEFFECTS:
!!! The 0 node has the total for each values.  If this is a distruted
!!! calculation, MPI is used to do this.  For serial runs this does
!!! nothing.
  subroutine rismthermo_mpireduce(this)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
    type(rismthermo_t), intent(inout) :: this
#ifdef MPI
    _REAL_ :: buffer_copy(ubound(this%mpi_buffer, 1))
    integer :: err
    buffer_copy = this%mpi_buffer
#  ifdef USE_MPI_IN_PLACE
    if (this%mpirank == 0) then
       call mpi_reduce(MPI_IN_PLACE, this%mpi_buffer, &
            ubound(this%mpi_buffer, 1), MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, this%mpicomm, err)
    else
       call mpi_reduce(this%mpi_buffer, this%mpi_buffer, &
            ubound(this%mpi_buffer, 1), MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, this%mpicomm, err)
    end if
#  else
    call mpi_reduce(this%mpi_buffer, this%tmpi_buffer, &
         ubound(this%mpi_buffer, 1), MPI_DOUBLE_PRECISION, &
         MPI_SUM, 0, this%mpicomm, err)
    this%mpi_buffer = this%tmpi_buffer
#  endif /*USE_MPI_IN_PLACE*/
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
#if !defined( USE_MPI_IN_PLACE) && defined( MPI )
    if (safemem_dealloc(this%tmpi_buffer) /= 0 ) &
         call rism_report_error("Dealloc failed in rism_thermo")
#endif
    nullify(this%excessChemicalPotential)
    nullify(this%excessChemicalPotentialGF)
    nullify(this%excessChemicalPotentialPCPLUS)
    nullify(this%excessChemicalPotentialUC)
    nullify(this%solventPotentialEnergy)
    nullify(this%solvationEnergy)
    nullify(this%solvationEnergyGF)
    nullify(this%solvationEnergyPCPLUS)
    nullify(this%solvationEnergyUC)
    nullify(this%partialMolarVolume)
    nullify(this%partialMolarVolume_dT)
    nullify(this%excessParticlesBox)
    nullify(this%totalParticlesBox)
    nullify(this%excessParticles)
    nullify(this%excessParticles_dT)
    nullify(this%kirkwoodBuff)
    nullify(this%kirkwoodBuff_dT)
    nullify(this%DCFintegral)
    nullify(this%DCFintegral_dT)
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
!! In SFF:
!!  rism_setparam_( &rismData, &xvvlen, xvvfile,
!!                  &guvlen, guvfile, &huvlen, huvfile, &cuvlen, cuvfile,
!!                  &uuvlen, uuvfile, &asymplen, asympfile,
!!                  &quvlen, quvfile, &chgdistlen, chgdistfile,
!!                  &volfmtlen, volfmt,
!!                  &comm,
!!                  &gb, &(prm->NumAtoms), &(prm->Ntypes), prm->Charges,
!!                  prm->Masses, prm->Cn1, prm->Cn2, prm->Iac, prm->Cno);
!!  rism_init_(&comm);
!!
!! MPI and other subroutines and functions: All public functions are
!! safe to call from all nodes and usually must be.  Only
!! RISM_THERMO_PRINT and RISM_SETPARAM must exclusively be called by
!! the master.
module amber_rism_interface
#ifndef CUDA
  use rism3d_c
#else
  use rism3d_rism3d_c_mod
#endif /*CUDA*/
  use rism3d_solvent_c
  use rism3d_solute_c
#ifdef RISM_CRDINTERP
  use fce_c
#endif /*RISM_CRDINTERP*/
  use rism_report_c
  use rism_timer_c
  use safemem
  use rismthermo_c

  !> Parameter derived type for storing calculation parameters. This can be use
  !! to transfer parameters from a C program.
  !! This must match the RismData struct in 'sff.h'.
  type rismprm_t
     sequence
     !> Cutoff for rism calculations (separate from SANDER non-bond).
     _REAL_ :: solvcut
     !> Buffer distance to the edge of the box for the solvent.
     _REAL_ :: buffer
     !> Grid spacing for all of the grids.
     _REAL_ :: grdspc(3)
     !> Box size for 3d-rism. For PBC calculations, these should generally be equal.
     _REAL_ :: solvbox(3)
     !> 'Step size' for MDIIS.
     _REAL_ :: mdiis_del
     !> Restart threshold factor. Ratio of the current residual to the
     !! minimum residual in the basis that causes a restart.
     _REAL_ :: mdiis_restart
     !> Fce cutoff distance.
     _REAL_ :: fcecut
     !> ???
     _REAL_ :: fceenormsw
     !> Coefficients for the temperature dependent universal correction.
     !! a,b,a1,b1
     _REAL_ :: uccoeff(4)
     !> Uniform bias applied to the electrostatic potential.
     _REAL_ :: biasPotential

     !> treecode multipole acceptance criterion parameter for DCF
     _REAL_ :: treeDCFMAC
     !> treecode multipole acceptance criterion parameter for TCF
     _REAL_ :: treeTCFMAC
     !> treecode multipole acceptance criterion parameter for Coulomb
     _REAL_ :: treeCoulombMAC

     !> long-range asymptotics k-space cut off tolerance.  Only grid
     !! points that have an approximate value greater than this will be computed.
     !! -1 - cutoff selected from calculation tolerance
     !! 0  - no cutoff applied
     !! >0 - used as the tolerance.  Should be a small value.  E.g., 1e-7
     _REAL_ :: asympKSpaceTolerance
     !> Tolerance for LJ potential calculation
     !! -1 - cutoff selected from calculation tolerance
     !! 0  - no cutoff applied
     !! >0 - used as the tolerance.  Should be a small value.  E.g., 1e-7
     _REAL_ :: ljtolerance
     !> Charge smearing parameter for long-range asymptotics and Ewald,
     !! typically eta in the literature
     _REAL_ :: chargeSmear
     
     !> For backwards compatibility, we still need to read this.
     integer :: closureOrder
     !> Number of grid points in each dimension.
     integer :: ng3(3)
     !> Use 3D-RISM.
     integer :: rism
     !> Use long range asymptotic corrections for thermodynamic calculations.
     logical*4 :: asympcorr
     !> Number of vectors used for MDIIS (consequently, the number of
     !! copies of CUV we need to keep for MDIIS).
     integer :: mdiis_nvec
     !> Mdiis implementation to use.
     integer :: mdiis_method
     !> Maximum number of rism relaxation steps.
     integer :: maxstep
     !> Number of past cuv time steps saves.
     integer :: npropagate
     !> Center the solute in the solvation box.
     !!          0 - off
     !!          1 - center of mass
     !!          2 - center of geometry
     !!          3 - center of mass shifted to the nearest grid point
     !!          4 - center of geometry shifted to the nearest grid point
     !! For negative numbers the centering translation is only
     !! calculated for the first solution and used for subsequent
     !! calculations.  This allows the solute to drift in the box.
     integer :: centering
     !> 0 - Do nothing.
     !! 1 - Redistribute forces to get zero total force.
     integer :: zerofrc
     !> If 0, the 3D-RISM solution is calculated but the resulting forces are not.
     integer :: apply_rism_force
     !> Do a polar/apolar decomposition of the chemical potential.
     integer :: polarDecomp
     !> Do a energy/entropy decomposition of the chemical potential.
     integer :: entropicDecomp
     !> Include Gaussian fluctuation functional in output
     integer :: gfCorrection
     !> Include PC+/3D-RISM correction in output
     integer :: pcplusCorrection
     !> Size of rism multiple timestep.
     integer :: rismnrespa
     !> Fce MTS stride length.
     integer :: fcestride
     !> Number of allFCE basis vectors.
     integer :: fcenbasis
     !> Number of leading FCE basis vectors.
     integer :: fcenbase
     !> Fce coordinate basis type.
     integer :: fcecrd
     !> ???
     integer :: fceweigh
     !> Type of the force extrapolation.
     integer :: fcetrans
     !> Sorting of coordinate-force pairs.
     integer :: fcesort

     ! New FCE parameters, for GSFE (Generalized Solvent Force Extrapolation)

     !fceifreq :: updating frequency of the extended to basic mapping list
     integer :: fceifreq
     !fcentfrcor :: net force correction due to individual extrapolation
     integer :: fcentfrcor
     !fcewrite :: write FCE basis set at each full 3D-RISM solution (1=yes)
     integer :: fcewrite
     !fcewrite :: read in FCE basis set from files 'fcecoord.rst' & 'fceforce.rst' (1=yes)
     integer :: fceread
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

     !> treecode order of Taylor approximation for DCF
     integer :: treeDCFOrder
     !> treecode order of Taylor approximation for TCF
     integer :: treeTCFOrder
     !> treecode order of Taylor approximation for Coulomb
     integer :: treeCoulombOrder

     !> maximum leaf size of tree for DCF
     integer :: treeDCFN0
     !> maximum leaf size of tree for TCF
     integer :: treeTCFN0
     !> maximum leaf size of tree for Coulomb
     integer :: treeCoulombN0

     !> perform internal consistency test. Done after output.
     logical*4 :: selftest

     !> do treecode DCF
     logical*4 :: treeDCF
     !> do treecode TCF
     logical*4 :: treeTCF
     !> do treecode Coulomb
     logical*4 :: treeCoulomb

     !> output molecular_reconstruction
     logical*4 :: molReconstruct

     !> fail if solution does not converge
     logical*4 :: stoponfail
     
     !> BPR: Make sure the number of INTEGERs is not odd, to shut up a
     !! compiler whinge about misalignment.
     !! Note: this should be commented if ever there are grounds for a
     !! new real INTEGER.
      ! integer :: padding
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
#ifdef RISM_CRDINTERP  
  type(fce), save :: fce_o
#endif /*RISM_CRDINTERP*/
  type(rism3d_solvent), save :: solvent
  type(rism3d_solute), save :: solute
#ifdef CUDA
  type(solvent_cpp), save :: solv_cpp
  type(solute_cpp), save :: solu_cpp
  type(timer_cpp), save :: timer_c
#endif /*CUDA*/
  type(rism_timer), save :: timer
  type(rism_timer), save :: timer_write
  type(rism_timer), save :: timer_init
  integer :: pa_orient, rmsd_orient
  _REAL_ :: centerOfMass(3)
  ! Read from prmtop file during parameter processing and passed to
  ! child MPI processes.
  _REAL_ :: unitCellDimensions(6)

  integer :: outunit

  !> If true, calculate Universal Correction excess chemical
  !! potential and, if possible, solvation energy and entropy
  logical :: canCalculateUC

  !> If true, calculate partial molar volume temperature, Universal
  !! Correction and PCPLUS solvation ! temperature derivative. I.e., the
  !! solvation energy and entropy
  logical :: canCalculatePartialMolarVolume_dT
  
  !> List of closures to use in order.  Only the last closure is used
  !! for thermodynamic output.  This can be used to progressively
  !! increase the order of the closure to aid convergence.  Closure
  !! types: KH, HNC or PSEn where n is an integer. This is initialized
  !! to a length of 10 in defaults() to allow it to be used with the
  !! namelist in sander.
  character(len=8), pointer :: closurelist(:) => NULL()
  !> Name of potential used for periodic calculations.
  character(len=8) :: periodicPotential = ''
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
  !asympfile     : (output) long range asymptotics function.  Volumetric file.
  !quvfile       : (output) charge density distribution.  Volumetric file.
  !chgDistFile   : (output) charge distribution. Volumetric file.
  !excessChemicalPotentialfile    : (output) excess chemical potential map [kcal/mol/A^3]. Volumetric file.
  !solvationEnergyfile   : (output) solvation energy map [kcal/mol/A^3]. Volumetric file.
  !entropyfile   : (output) solvent entroy (-TS) map [kcal/mol/A^3]. Volumetric file.
  !excessChemicalPotentialGFfile  : (output) Gaussian fluctuation excess chemical potential map [kcal/mol/A^3]. Volumetric file.
  !solvationEnergyGFfile : (output) Gaussian fluctuation solvation energy map [kcal/mol/A^3]. Volumetric file.
  !entropyGFfile : (output) Gaussian fluctuation solvent entroy (-TS) map [kcal/mol/A^3]. Volumetric file.
  !excessChemicalPotentialPCPLUSfile  : (output) PCPLUS excess chemical potential map [kcal/mol/A^3]. Volumetric file.
  !solvationEnergyPCPLUSfile : (output) PCPLUS solvation energy map [kcal/mol/A^3]. Volumetric file.
  !entropyPCPLUSfile : (output) PCPLUS solvent entroy (-TS) map [kcal/mol/A^3]. Volumetric file.
  !excessChemicalPotentialUCfile  : (output) Universal Correction excess chemical potential map [kcal/mol/A^3]. Volumetric file.
  !solvationEnergyUCfile : (output) Universal Correction solvation energy map [kcal/mol/A^3]. Volumetric file.
  !entropyUCfile : (output) Universal Correction solvent entroy (-TS) map [kcal/mol/A^3]. Volumetric file.
  !solventPotentialEnergyfile     : (output) solvent-solute potential energy map [kcal/mol/A^3]. Volumetric file.
  !electronMapFile : (output) solvent electron density map. Volumetric file.
  !periodicPotential : Specify periodic potential used for periodic calculations.
  !                    Either 'ewald' or 'pme' or 'pmekh'.
  !volfmt        : either 'mrc', 'ccp4', 'dx', or 'xyzv'
  character(len=256) :: xvvfile='', guvfile='', huvfile='', cuvfile='', &
       uuvfile='', asympfile='', quvFile='', chgDistFile='', &
       excessChemicalPotentialfile='', solvationEnergyfile='', entropyfile='', &
       excessChemicalPotentialGFfile='', solvationEnergyGFfile='', entropyGFfile='', &
       excessChemicalPotentialPCPLUSfile='', solvationEnergyPCPLUSfile='', entropyPCPLUSfile='', &
       excessChemicalPotentialUCfile='', solvationEnergyUCfile='', entropyUCfile='', &
       solventPotentialEnergyfile='', electronMapFile='', &
       volfmt='mrc', crdFile=''

  integer :: mpirank = 0, mpisize = 1, mpicomm = 0

  !working memory for rism_force() so it is not reallocated every time
  !ff :: forces
  !atomPositions_fce :: coordinates for FCE
  _REAL_, pointer :: ff(:, :) => NULL()
#if defined(RISM_CRDINTERP)
  _REAL_, pointer :: atomPositions_fce(:, :) => NULL()
#endif /*RISM_CRDINTERP*/

  private :: rism_mpi_bcast

  

end module amber_rism_interface

#ifdef SANDER
module sander_rism_interface
  use amber_rism_interface
  implicit none
contains
#endif


  !> Sets all input parameters for 3D-RISM.  This _must_ be called by the head
  !! node and may be called by all nodes.  If all nodes call this
  !! subroutine they must all agree on the value of rismprm%rism
  !! (SANDER) or userData%rism(SFF).
  !!
  !! SANDER prerequisites:
  !!   - Names of 3D-RISM specific I/O files (Xvv, Guv, etc.) are specified on the
  !!    command line and are read by mdfil.F90 into the variables in
  !!    AMBER_RISM_INTERFACE.
  !!   - igb should be set to 6 (vacuum electrostatics).
  !! SANDER IN:
  !! @param[in] mdin Name of the mdin file that SANDER name lists are read from.
  !!
  !! SFF prerequisites:
  !!   - All user options (including file names) are read from mm_options.l.
  !!    Non-string options are read into a C struct equivalent to rismprm_t.
  !!    String options are supplied as char* array, integer length pairs.
  !!   - gb should be set to 0 (vacuum electrostatics).
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
  !! @param[in] asymplen Length of asympchar array.
  !! @param[in] asympchar Character array for asymptotics output file name.
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
  !! @param[in] excessChemicalPotentialGFlen Length of excessChemicalPotentialGFchar array.
  !! @param[in] excessChemicalPotentialGFchar Character array for Gaussian fluctuation
  !!   excessChemicalPotential output file name.
  !! @param[in] solvationEnergyGFlen Length of solvationEnergyGFchar array.
  !! @param[in] solvationEnergyGFchar Character array for Gaussian fluctuation
  !!   solvationEnergy output file name.
  !! @param[in] entropyGFlen Length of entropyGFchar array.
  !! @param[in] entropyGFchar Character array for Gaussian fluctuation
  !!   entropy output file name.
  !! @param[in] excessChemicalPotentialPCPLUSlen Length of excessChemicalPotentialPCPLUSchar array.
  !! @param[in] excessChemicalPotentialPCPLUSchar Character array for PCPLUS
  !!   excessChemicalPotential output file name.
  !! @param[in] solvationEnergyPCPLUSlen Length of solvationEnergyPCPLUSchar array.
  !! @param[in] solvationEnergyPCPLUSchar Character array for PCPLUS
  !!   solvationEnergy output file name.
  !! @param[in] entropyPCPLUSlen Length of entropyPCPLUSchar array.
  !! @param[in] entropyPCPLUSchar Character array for PCPLUS entropy output file name.
  !! @param[in] excessChemicalPotentialUClen Length of excessChemicalPotentialUCchar array.
  !! @param[in] excessChemicalPotentialUCchar Character array for Universal Correction
  !!   excessChemicalPotential output file name.
  !! @param[in] solvationEnergyUClen Length of solvationEnergyUCchar array.
  !! @param[in] solvationEnergyUCchar Character array for Universal Correction
  !!   solvationEnergy output file name.
  !! @param[in] entropyUClen Length of entropyUCchar array.
  !! @param[in] entropyUCchar Character array for Universal Correction
  !!   entropy output file name.
  !! @param[in] solventPotentialEnergylen Length of solventPotentialEnergychar array.
  !! @param[in] solventPotentialEnergychar Character array for solute-solvent potential
  !!   energy output file name.
  !! @param[in] electronMaplen Length of electronMapchar array.
  !! @param[in] electronMapchar Character array for the smeared
  !!   electron density map.
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
#ifdef SANDER
       mdin, &
#else /* not SANDER */
       userData, ntol, tol, &
       closurelen, nclosure, closurechar, &
       xvvlen, xvvchar, &
       guvlen, guvchar, huvlen, huvchar, cuvlen, cuvchar, &
       uuvlen, uuvchar, asymplen, asympchar, quvlen, quvchar, &    
       chgDistlen, chgDistchar, &
       excessChemicalPotentiallen, excessChemicalPotentialchar, &
       solvationEnergylen, solvationEnergychar, entropylen, entropychar, &
       excessChemicalPotentialGFlen, excessChemicalPotentialGFchar, &
       solvationEnergyGFlen, solvationEnergyGFchar, entropyGFlen, entropyGFchar, &
       excessChemicalPotentialPCPLUSlen, excessChemicalPotentialPCPLUSchar, &
       solvationEnergyPCPLUSlen, solvationEnergyPCPLUSchar, entropyPCPLUSlen, entropyPCPLUSchar, &
       excessChemicalPotentialUClen, excessChemicalPotentialUCchar, &
       solvationEnergyUClen, solvationEnergyUCchar, entropyUClen, entropyUCchar, &
       solventPotentialEnergylen, solventPotentialEnergychar, &
       electronMaplen, electronMapchar, &
       volfmtlen, volfmtchar, &
       periodiclen, periodicchar, &
       rstlen, rstchar, &
#endif /*SANDER*/
       comm, &
       numAtoms, numTypes, &
       charge, mass, ljA, ljB, &
       atomTypeIndex, nonbondedParmIndex)
    use amber_rism_interface
    use rism_io, only : readUnitCellDimensionsFromCrd
    use constants, only : KB
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
    integer, intent(in) :: numAtoms, numTypes, atomTypeIndex(numTypes**2), nonbondedParmIndex(numAtoms)
    _REAL_, intent(in) :: charge(numAtoms), mass(numAtoms), ljA(numTypes*(numTypes + 1)/2), &
         ljB(numTypes*(numTypes + 1)/2)
#ifdef SANDER
    !  character(*), intent(in) :: xvvfile, mdin
    character(*), intent(in) :: mdin
    integer :: mdin_unit=55
#else /*SANDER*/
    type(rismprm_t), intent(in) :: userData;
    integer, intent(in) :: ntol
    _REAL_, intent(in)  :: tol(ntol)
    integer, intent(in) :: closurelen, nclosure, &
         xvvlen, guvlen, huvlen, cuvlen, uuvlen, &
         asymplen, quvlen, chgDistLen, &
         excessChemicalPotentiallen, solvationEnergylen, entropylen, &
         excessChemicalPotentialGFlen, solvationEnergyGFlen, entropyGFlen, &
         excessChemicalPotentialPCPLUSlen, solvationEnergyPCPLUSlen, entropyPCPLUSlen, &
         excessChemicalPotentialUClen, solvationEnergyUClen, entropyUClen, &
         solventPotentialEnergylen, &
         electronMaplen, &
         volfmtlen, periodiclen, rstlen
    integer(kind=1), intent(in) :: closurechar(closurelen*nclosure), &
         xvvchar(xvvlen + 1), guvchar(guvlen + 1), &
         huvchar(huvlen + 1), cuvchar(cuvlen + 1), uuvchar(uuvlen + 1), &
         asympchar(asymplen + 1), quvchar(quvlen + 1), chgDistchar(chgDistlen + 1), &
         excessChemicalPotentialchar(excessChemicalPotentiallen + 1), &
         solvationEnergychar(solvationEnergylen + 1), entropychar(entropylen + 1), &
         excessChemicalPotentialGFchar(excessChemicalPotentiallen + 1), &
         solvationEnergyGFchar(solvationEnergylen + 1), entropyGFchar(entropylen + 1), &
         excessChemicalPotentialPCPLUSchar(excessChemicalPotentiallen + 1), &
         solvationEnergyPCPLUSchar(solvationEnergylen + 1), entropyPCPLUSchar(entropylen + 1), &
         excessChemicalPotentialUCchar(excessChemicalPotentiallen + 1), &
         solvationEnergyUCchar(solvationEnergylen + 1), entropyUCchar(entropylen + 1), &
         solventPotentialEnergychar(solventPotentialEnergylen + 1), &
         electronMapchar(electronMaplen + 1), &
         volfmtchar(volfmtlen + 1), periodicchar(periodiclen + 1), &
         rstchar(rstlen + 1)
    !  character(xvvlen):: xvvfile
#endif /*SANDER*/
    integer, intent(in) :: comm
    character(len=16) :: whtspc
    integer :: i, stat, err
    !iclosure :: counter for closures
    !iclosurechar :: current index in the closurechar array from sff
    integer :: iclosure, iclosurechar
    logical :: op
    integer ::  un

    ! In case this is not the first time init has been called (i.e.,
    ! multiple runs) destroy timers and re-create them.
    call rism_timer_destroy(timer_write)
    call rism_timer_destroy(timer_init)
    call rism_timer_destroy(timer)
    call rism_timer_new(timer, "3D-RISM Total")
    call rism_timer_new(timer_init, "3D-RISM initialization", timer)
    call rism_timer_start(timer_init)
    call rism_timer_new(timer_write, "3D-RISM Output", timer)

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
#ifdef SANDER
    if (rismprm%rism == 0) then
       call rism_timer_stop(timer_init)
       return
    end if
#else /*SANDER*/
    if (userData%rism == 0) then
       call rism_timer_stop(timer_init)
       return
    end if
#endif /*SANDER*/

    ! Rank 0 only.
    if (mpirank /= 0) then
       call rism_timer_stop(timer_init)
       return
    end if


    outunit = rism_report_getMUnit()
    call defaults()

#ifdef SANDER
    inquire(file=mdin, opened=op, number=un)
    if (op) mdin_unit=un
    open(unit=mdin_unit, file=mdin, status='OLD', form='FORMATTED', iostat=stat)
    if (stat/= 0) then
       call rism_report_error('(a, i4)', "opening "//trim(mdin)//"failed. IOSTAT=", stat)
    end if
    call read_namelist(mdin_unit)
    if (.not.op) close(unit=mdin_unit)
#else /*.not.SANDER*/
    if (ntol /= 0) then
       tolerancelist => safemem_realloc(tolerancelist, ntol)
       tolerancelist = tol(1:ntol)
    end if
    call update_param(userData)
    iclosurechar = 1
    closurelist => safemem_realloc(closurelist, len(closurelist), max(nclosure, 1))
    do iclosure = 1, nclosure
       ! The default in SFF is to pass an empty string, but, since
       ! this would overwrite our defaults, ignore it.
       if (closurechar((iclosure - 1) * closurelen + 1) == 0) exit
       call cstr2fstr(closurelist(iclosure), closurechar((iclosure-1)*closurelen + 1), closurelen)
    end do
    call cstr2fstr(xvvfile, xvvchar, xvvlen)
    call cstr2fstr(guvfile, guvchar, guvlen)
    call cstr2fstr(huvfile, huvchar, huvlen)
    call cstr2fstr(cuvfile, cuvchar, cuvlen)
    call cstr2fstr(uuvfile, uuvchar, uuvlen)
    call cstr2fstr(asympfile, asympchar, asymplen)
    call cstr2fstr(quvfile, quvchar, quvlen)
    call cstr2fstr(chgDistfile, chgDistchar, chgDistlen)
    call cstr2fstr(excessChemicalPotentialfile, excessChemicalPotentialchar, excessChemicalPotentiallen)
    call cstr2fstr(solvationEnergyfile, solvationEnergychar, solvationEnergylen)
    call cstr2fstr(entropyfile, entropychar, entropylen)
    call cstr2fstr(excessChemicalPotentialGFfile, excessChemicalPotentialGFchar, excessChemicalPotentialGFlen)
    call cstr2fstr(solvationEnergyGFfile, solvationEnergyGFchar, solvationEnergyGFlen)
    call cstr2fstr(entropyGFfile, entropyGFchar, entropyGFlen)
    call cstr2fstr(excessChemicalPotentialPCPLUSfile, excessChemicalPotentialPCPLUSchar, excessChemicalPotentialPCPLUSlen)
    call cstr2fstr(solvationEnergyPCPLUSfile, solvationEnergyPCPLUSchar, solvationEnergyPCPLUSlen)
    call cstr2fstr(entropyPCPLUSfile, entropyPCPLUSchar, entropyPCPLUSlen)
    call cstr2fstr(excessChemicalPotentialUCfile, excessChemicalPotentialUCchar, excessChemicalPotentialUClen)
    call cstr2fstr(solvationEnergyUCfile, solvationEnergyUCchar, solvationEnergyUClen)
    call cstr2fstr(entropyUCfile, entropyUCchar, entropyUClen)
    call cstr2fstr(solventPotentialEnergyfile, solventPotentialEnergychar, solventPotentialEnergylen)
    call cstr2fstr(electronMapFile, electronMapchar, electronMaplen)
    if (volfmtlen > 0) &
         call cstr2fstr(volfmt, volfmtchar, volfmtlen)
    call cstr2fstr(periodicPotential, periodicchar, periodiclen)
    call cstr2fstr(crdFile, rstchar, rstlen)
#endif /*SANDER*/

    ! Initialize 3D-RISM solute and solvent.
    call rism3d_solvent_new(solvent, xvvfile)    
    call rism3d_solute_new_sander(solute, numAtoms, numTypes, atomTypeIndex, &
         nonbondedParmIndex, charge, ljA, ljB, mass, solvent%temperature)

    if (periodicPotential /= '') then
       call readUnitCellDimensionsFromCrd(crdFile, unitCellDimensions)
    end if

    call sanity_check()

#ifdef CUDA
    call cpy_solute(solu_cpp,solute)
    call cpy_solvent(solv_cpp,solvent)
#endif /*CUDA*/
    
#ifdef SANDER
    if (rismprm%rism >= 1) then
       write(outunit, '(a)') "3D-RISM:"
       if (rismprm%rism < 1) then
          write(outunit, '(5x, a, i10)') 'irism   =', rismprm%rism
       else if (rismprm%rism == 1) then
          write(outunit, '(5x, 3(a10, "=", 100a10))') &
               'closure'//whtspc, closurelist
          write(outunit, '(5x, a10, "= ",1p, 4(e12.5, 1x))') &
               'uccoeff'//whtspc, rismprm%uccoeff
          write(outunit, '(5x, 3(a10, "="f10.5))') &
               'solvcut'//whtspc, rismprm%solvcut, &
               ', buffer'//whtspc, rismprm%buffer
          write(outunit, '(5x, a10, "=", 3(f10.5, 1x))') &
               'grd_spc'//whtspc, rismprm%grdspc
          write(outunit, '(5x, a10, "=", 3(i10, 1x))') &
               'ng3'//whtspc, rismprm%ng3
          write(outunit, '(5x, a10, "=", 3(f10.5, 1x))') &
               'solvbox'//whtspc, rismprm%solvbox
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
          write(outunit, '(5x, 3(a10, "=", i10))') &
               'centering'//whtspc, rismprm%centering, &
               ', zerofrc'//whtspc, rismprm%zerofrc
          write(outunit, '(5x, a10, "=", i10, a11, "=", l10)') &
               'apply_rism_force'//whtspc, rismprm%apply_rism_force, &
               ', asympcorr'//whtspc, rismprm%asympcorr
          write(outunit, '(5x, 2(a10, "=", i10), a10, "=", f10.5)') &
               'rismnrespa'//whtspc, rismprm%rismnrespa, &
               ', fcestride'//whtspc, rismprm%fcestride, &
               ', fcecut'//whtspc, rismprm%fcecut
          write(outunit, '(5x, 3(a10, "=", i10))') &
               'fcenbasis'//whtspc, rismprm%fcenbasis, &
               ', fcenbase'//whtspc, rismprm%fcenbase, &
               ', fcecrd'//whtspc, rismprm%fcecrd
          write(outunit,'(3(a15,"=",i10))') '|     fceweigh ', rismprm%fceweigh, &
               ', fcetrans      ', rismprm%fcetrans,   ', fcesort    ', rismprm%fcesort
          write(outunit,'(a15,"=",i10,a12,"=",d10.2,a12,"=",i10)') '|     fceifreq ', &
             rismprm%fceifreq,      ', fceenormsw', rismprm%fceenormsw, &
             ', fcentfrcor', rismprm%fcentfrcor
          write(outunit,'(a15,"=",i5, a20,"=",i5)') '|     fcewrite ', rismprm%fcewrite, &
               ', fceread  ', rismprm%fceread
          write(outunit, '(5x, 2(a20, "=", i10))') &
               'polarDecomp'//whtspc, rismprm%polardecomp, &
               ', entropicDecomp'//whtspc, rismprm%entropicDecomp
          write(outunit, '(5x, 2(a20, "=", i10))') &
               'gfCorrection'//whtspc, rismprm%gfCorrection, &
               ', pcplusCorrection'//whtspc, rismprm%pcplusCorrection
          write(outunit, '(5x, 2(a20, "= ", a8))') &
               'periodic'//whtspc, periodicPotential
          write(outunit, '(5x, 1(a10, "=", i10), a10, "=  ", a8)') &
               'write_thermo'//whtspc, rismprm%write_thermo, &
               ', volfmt'//whtspc, volfmt
          write(outunit, '(5x, 3(a10, "=", i10))') &
               'saveprogress'//whtspc, rismprm%saveprogress, &
               ', ntwrism'//whtspc, rismprm%ntwrism, &
               ', verbose'//whtspc, rismprm%verbose
          write(outunit, '(5x, a15, "=", l5)') 'stoponfail'//whtspc, rismprm%stoponfail
          write(outunit, '(5x, a15, "=", l5)') 'molReconstruct'//whtspc, rismprm%molReconstruct
          write(outunit, '(5x, 3(a10, "=", i10))') &
               'progress'//whtspc, rismprm%progress
          write(outunit, '(5x, 2(a15, "=", l5), 1(a19 "=", l5))') &
               'treeDCF'//whtspc, rismprm%treeDCF, &
               ', treeTCF'//whtspc, rismprm%treeTCF, &
               ', treeCoulomb'//whtspc, rismprm%treeCoulomb
          write(outunit, '(5x, 2(a15, "=", f6.3), 1(a19, "=", f6.3))') &
               'treeDCFMAC'//whtspc, rismprm%treeDCFMAC, &
               ', treeTCFMAC'//whtspc, rismprm%treeTCFMAC, &
               ', treeCoulombMAC'//whtspc, rismprm%treeCoulombMAC
          write(outunit, '(5x, 2(a15, "=", i5), 1(a19, "=", i5))') &
               'treeDCFOrder'//whtspc, rismprm%treeDCFOrder, &
               ', treeTCFOrder'//whtspc, rismprm%treeTCFOrder, &
               ', treeCoulombOrder'//whtspc, rismprm%treeCoulombOrder
          write(outunit, '(5x, 2(a15, "=", i5), 1(a19, "=", i5))') &
               'treeDCFN0'//whtspc, rismprm%treeDCFN0, &
               ', treeTCFN0'//whtspc, rismprm%treeTCFN0, &
               ', treeCoulombN0'//whtspc, rismprm%treeCoulombN0
          write(outunit, '(5x, a21, "=", es10.2, a14, "=", f6.3)') &
               'asympKSpaceTolerance'//whtspc, rismprm%asympKSpaceTolerance, &
               ', chargeSmear'//whtspc, rismprm%chargeSmear
          write(outunit, '(5x, a10, "=", es10.2)') &
               'ljTolerance'//whtspc, rismprm%ljTolerance
          write(outunit, '(5x, a10, "=", f10.5)') &
               'biasPotential'//whtspc, rismprm%biasPotential
          !!!! Extrapolation method warnings
          if(rismprm%fcetrans/=0.and.rismprm%fcetrans/=1.and.rismprm%fcetrans/= &
               2.and.rismprm%fcetrans/=3.and.rismprm%fcetrans/=4.and.rismprm%fcetrans/=5 & 
            .and.rismprm%fcetrans/=6) then
             write(6,'(/,a)') 'trans must be equal to 0, 1, 2, 3, 4, 5, or 6'
             stop
          end if
          if(rismprm%fcenbasis < rismprm%fcenbase) then
             write(6,'(/,a)') 'nbasis must be equal or larger than nbase'
             stop
          end if
          if(rismprm%fceenormsw < 0 .and. rismprm%fcetrans == 2) then
             write(6,'(/,a)') 'enormsw must be positive if fcetrans = 2'
             stop
          end if
          if(rismprm%fceenormsw < 0 .and. rismprm%fcetrans == 3) then
             write(6,'(/,a)') 'enormsw must be positive if fcetrans = 3'
             stop
          end if
          if(rismprm%fceenormsw < 0 .and. rismprm%fcetrans == 5) then
             write(6,'(/,a)') 'enormsw must be positive if fcetrans = 5'
             stop
          end if
          if(rismprm%fceenormsw < 0 .and. rismprm%fcetrans == 6) then
             write(6,'(/,a)') 'enormsw must be positive if fcetrans = 6'
             stop
          end if
          if(rismprm%fceifreq <= 0 .and. rismprm%fcetrans == 2) then
             write(6,'(/,a)') 'ifreq must be larger than 0 if fcetrans = 2'
             stop
          end if
          if(rismprm%fceifreq <= 0 .and. rismprm%fcetrans == 3) then
             write(6,'(/,a)') 'ifreq must be larger than 0 if fcetrans = 3'
             stop
          end if
          if(rismprm%fceifreq <= 0 .and. rismprm%fcetrans == 5) then
             write(6,'(/,a)') 'ifreq must be larger than 0 if fcetrans = 5'
             stop
          end if
          if(rismprm%fceifreq <= 0 .and. rismprm%fcetrans == 6) then
             write(6,'(/,a)') 'ifreq must be larger than 0 if fcetrans = 6'
             stop
          end if
          if(rismprm%fcenbasis > rismprm%fcenbase .and. rismprm%fcetrans == 4) then
             write(6,'(/,a)') 'nbasis must be equal to nbase if fcetrans = 4'
             stop
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       end if
       call flush(outunit)
    end if
#endif /*SANDER*/

    call rism_timer_stop(timer_init)

  end subroutine rism_setparam


  !! Obtain the periodic potential early in sander input processing.
  !! TODO: This is a hacky solution to ensure igb is set correctly in
  !! mdread1.F90 which is before rism parameters are normally are set
  !! by sander. rism_setparam ends up doing the same parmeter file
  !! parsing later. Either all rism processing could be done when this
  !! function is called, allowing this to be merged into
  !! rism_setparam, or else only the parameter file could be parsed
  !! here, making rism_setparam not need to do so later.
  subroutine rism_getPeriodicPotential(mdin, periodicPotentialLabel)
    use amber_rism_interface
    implicit none
    character(*), intent(in) :: mdin
    character(len=8), intent(out) :: periodicPotentialLabel
    ! character(len=8) :: periodic = ''
    integer :: un
    integer :: stat
    logical :: op

    character(len=8) :: closure(nclosuredefault)
    _REAL_ :: tolerance(nclosuredefault)
    integer :: closureOrder
    integer :: mdin_unit
    logical :: asympCorr
    integer :: entropicDecomp
    integer :: polarDecomp
    integer :: gfcorrection
    integer :: pcpluscorrection
    character(len=8) :: periodic = ''
    _REAL_ :: uccoeff(size(rismprm%uccoeff))
    _REAL_ :: biasPotential
    _REAL_ :: solvcut
    _REAL_ :: buffer
    _REAL_ :: grdspc(3)
    integer ::  ng3(3)
    _REAL_ :: solvbox(3)
    _REAL_ :: mdiis_del
    integer :: mdiis_nvec
    integer :: mdiis_method
    _REAL_ :: mdiis_restart
    integer :: maxstep
    logical :: stoponfail
    integer :: npropagate
    integer :: centering
    integer :: zerofrc
    integer :: apply_rism_force
    integer :: rismnrespa
    integer :: fcestride
    _REAL_ :: fcecut
    integer ::  fcenbasis
    integer ::  fcenbase
    integer ::  fcecrd
    integer ::  fceweigh
    integer ::  fcetrans
    integer ::  fcesort
    integer ::  fceifreq
    integer ::  fcentfrcor
    integer ::  fcewrite
    integer ::  fceread
    _REAL_  ::  fceenormsw
    integer :: saveprogress
    integer :: ntwrism
    integer :: verbose
    integer :: progress
    logical :: treeDCF
    logical :: treeTCF
    logical :: treeCoulomb
    _REAL_ :: treeDCFMAC
    _REAL_ :: treeTCFMAC
    _REAL_ :: treeCoulombMAC
    integer :: treeDCFOrder
    integer :: treeTCFOrder
    integer :: treeCoulombOrder
    integer :: treeDCFN0
    integer :: treeTCFN0
    integer :: treeCoulombN0
    _REAL_ :: asympKSpaceTolerance
    _REAL_ :: ljTolerance
    _REAL_ :: chargeSmear
    logical :: molReconstruct
#ifdef SANDER
    integer :: write_thermo
#endif
    namelist /rism/ &
         ! closure
         closure, closureOrder, uccoeff, &
         entropicDecomp, polarDecomp, &
         gfCorrection, pcpluscorrection,&
         biasPotential, &
         ! thermodynamics
         asympCorr, periodic, &
         ! solvation box
         buffer, grdspc, solvcut, &
         ng3, solvbox, &
         ! convergence
         tolerance, mdiis_del, mdiis_nvec, mdiis_method, &
         mdiis_restart, maxstep, npropagate, stoponfail, &
         ! minimization
         centering, zerofrc, &
         ! imin=5
         apply_rism_force, pa_orient, rmsd_orient, &
         ! md
         rismnrespa, &
         ! treecode and direct sum
         treeDCF, treeTCF, treeCoulomb, &
         treeDCFMAC, treeTCFMAC, treeCoulombMAC, &
         treeDCFOrder, treeTCFOrder, treeCoulombOrder, &
         treeDCFN0, treeTCFN0, treeCoulombN0, &
         ! long-range asymptotics parameters
         asympKSpaceTolerance, chargeSmear, &
         ! LJ potential
         ljTolerance, &
         ! molecular reconstruction
         molReconstruct, &
#ifdef RISM_CRDINTERP
         fcestride, fcecut, fcenbasis, fcenbase, fcecrd, &
         fceweigh, fcetrans, fcesort, fceifreq, fcentfrcor, &
         fcewrite, fceread, fceenormsw, &
#endif /*RISM_CRDINTERP*/
         !output
#ifdef SANDER
         write_thermo, &
#endif
         saveprogress, ntwrism, verbose, progress, volfmt

    inquire(file=mdin, opened=op, number=un)
    if (op) mdin_unit=un
    open(unit=mdin_unit, file=mdin, status='OLD', form='FORMATTED', iostat=stat)
    if (stat /= 0) then
       call rism_report_error('(a, i4)', "opening "//trim(mdin)//"failed. IOSTAT=", stat)
    end if

    
    rewind(mdin_unit)
    read(mdin_unit, nml=rism)
    
    if (.not. op) close(unit=mdin_unit)

    periodicPotentialLabel = periodic
    
  end subroutine rism_getPeriodicPotential

  
  !> Performs all of the initialization required for 3D-RISM
  !! calculations. All parameters should have already been set with
  !! RISM_SETPARAM. For MPI calculations this _must_ be called by all
  !! processes. Only the head node irism value is used.
  !! @param[in] comm MPI communicator.
  subroutine rism_init(comm)
    use amber_rism_interface
    use safemem
#ifdef RISM_CRDINTERP
    use fce_c
#endif /*RISM_CRDINTERP*/
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
    integer, intent(in) :: comm

    integer :: err

    
    ! Ensure that timers have been created incase RISM_SETPARAM was
    ! only called by the master node.
    if (trim(timer%name) .ne. "3D-RISM Total") &
         call rism_timer_new(timer, "3D-RISM Total")
    if (trim(timer_init%name) .ne. "3D-RISM initialization") &
         call rism_timer_new(timer_init, "3D-RISM initialization", timer)
    call rism_timer_start(timer_init)
    if (trim(timer_write%name) .ne. "3D-RISM Output") &
         call rism_timer_new(timer_write, "3D-RISM Output", timer)

#ifdef SANDER
    call rism_report_setMUnit(6)
    call rism_report_setWUnit(6)
    call rism_report_setEUnit(6)
#endif /*SANDER*/

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
    if (rismprm%rism == 0) then
       call rism_timer_stop(timer_init)
       return
    end if

#ifdef RISM_CRDINTERP
    call fce_new(fce_o, solute%numAtoms, rismprm%fcenbasis, rismprm%fcenbase, &
         rismprm%fcecrd, rismprm%fceweigh, rismprm%fcetrans, &
         rismprm%fcesort, rismprm%fceifreq,rismprm%fcentfrcor, & 
         rismprm%fceenormsw, rismprm%fcecut, mpicomm)

    ! Read in saved extrapolation basis
    if(rismprm%fceread == 1) then
       call fce_readbasis(fce_o)
    endif
#endif /*RISM_CRDINTERP*/

    ! 3D-RISM may have already been initialized. In the absence of a
    ! subroutine to set all of these parameters individually, we
    ! destroy the original instance and re-initialize. Since this is a
    ! rare event, it should not add to the expense of the calculation
    ! in any practical way.
#ifndef CUDA
    call rism3d_destroy(rism_3d)
#else
    call rism_3d%delete
#endif /*CUDA*/

    if (rismprm%buffer >= 0) then
!     print*,"print ljSigma(1) and (6)", solute%ljSigma(1),solute%ljSigma(6)
!     print*,"solute numAtoms", solute%numAtoms
!     print*,"charge ", solute%charge
#ifndef CUDA
       call rism3d_new(rism_3d, solute, solvent, rismprm%centering, rismprm%npropagate, &
            closurelist, rismprm%solvcut, &
            rismprm%mdiis_nvec, rismprm%mdiis_del, rismprm%mdiis_method, rismprm%mdiis_restart, &
            rismprm%treeDCF, rismprm%treeTCF, rismprm%treeCoulomb, &
            rismprm%treeDCFMAC, rismprm%treeTCFMAC, rismprm%treeCoulombMAC, &
            rismprm%treeDCFOrder, rismprm%treeTCFOrder, rismprm%treeCoulombOrder, &
            rismprm%treeDCFN0, rismprm%treeTCFN0, rismprm%treeCoulombN0, &
            rismprm%asympKSpaceTolerance, rismprm%ljTolerance, rismprm%chargeSmear, &
            o_buffer=rismprm%buffer, o_grdspc=rismprm%grdspc, o_mpicomm=mpicomm, &
            o_periodic=periodicPotential, o_unitCellDimensions=unitCellDimensions, &
            o_biasPotential=rismprm%biasPotential)
#else 
       rism_3d = rism3d_new_0(solu_cpp, solv_cpp, rismprm%centering, rismprm%npropagate, &
            closurelist, rismprm%solvcut, rismprm%mdiis_nvec, rismprm%mdiis_del, rismprm%mdiis_method, &
            rismprm%mdiis_restart, rismprm%treeDCF, rismprm%treeTCF, rismprm%treeCoulomb, rismprm%treeDCFMAC, &
            rismprm%treeTCFMAC, rismprm%treeCoulombMAC, rismprm%treeDCFOrder, rismprm%treeTCFOrder, &
            rismprm%treeCoulombOrder, rismprm%treeDCFN0, rismprm%treeTCFN0, rismprm%treeCoulombN0, &
            rismprm%asympKSpaceTolerance, rismprm%ljTolerance, rismprm%chargeSmear, rismprm%buffer, &
            rismprm%grdspc, mpicomm, periodicPotential, unitCellDimensions, &
            rismprm%biasPotential)

       call rism3d_set_atomname(rism_3d, solvent%atomName)
       call rism3d_set_closurelist(rism_3d, closurelist, size(closurelist))
#endif /*CUDA*/
    else
#ifndef CUDA
       call rism3d_new(rism_3d, solute, solvent, rismprm%centering, rismprm%npropagate, &
            closurelist, rismprm%solvcut, &
            rismprm%mdiis_nvec, rismprm%mdiis_del, rismprm%mdiis_method, rismprm%mdiis_restart, &
            rismprm%treeDCF, rismprm%treeTCF, rismprm%treeCoulomb, &
            rismprm%treeDCFMAC, rismprm%treeTCFMAC, rismprm%treeCoulombMAC, &
            rismprm%treeDCFOrder, rismprm%treeTCFOrder, rismprm%treeCoulombOrder, &
            rismprm%treeDCFN0, rismprm%treeTCFN0, rismprm%treeCoulombN0, &
            rismprm%asympKSpaceTolerance, rismprm%ljTolerance, rismprm%chargeSmear, &
            o_boxlen=rismprm%solvbox, o_ng3=rismprm%ng3, o_mpicomm=mpicomm, &
            o_periodic=periodicPotential, o_unitCellDimensions=unitCellDimensions,&
            o_biasPotential=rismprm%biasPotential)
#else
       rism_3d = rism3d_new_1(solu_cpp, solv_cpp, rismprm%centering, rismprm%npropagate, &
            closurelist, rismprm%solvcut, &
            rismprm%mdiis_nvec, rismprm%mdiis_del, rismprm%mdiis_method, rismprm%mdiis_restart, &
            rismprm%treeDCF, rismprm%treeTCF, rismprm%treeCoulomb, &
            rismprm%treeDCFMAC, rismprm%treeTCFMAC, rismprm%treeCoulombMAC, &
            rismprm%treeDCFOrder, rismprm%treeTCFOrder, rismprm%treeCoulombOrder, &
            rismprm%treeDCFN0, rismprm%treeTCFN0, rismprm%treeCoulombN0, &
            rismprm%asympKSpaceTolerance, rismprm%ljTolerance, rismprm%chargeSmear, &
            rismprm%solvbox, rismprm%ng3, mpicomm, &
            periodicPotential, unitCellDimensions,&
            rismprm%biasPotential)

       call rism3d_set_atomname(rism_3d, solvent%atomName)
       call rism3d_set_closurelist(rism_3d, closurelist, size(closurelist))
#endif /*CUDA*/
    end if
    call rism3d_setverbosity(rism_3d, rismprm%verbose)
#ifndef CUDA
    call rism3d_setTimerParent(rism_3d, timer)
#else
    call rism3d_setTimerParent(rism_3d, timer_c)
#endif /*CUDA*/

    call rismthermo_new(rismthermo, rism3d_get_numAtomTypes(rism_3d), mpicomm)
    if (rismprm%polarDecomp == 1) then
       call rismthermo_new(rismthermo_pol, rism3d_get_numAtomTypes(rism_3d), mpicomm)
       call rismthermo_new(rismthermo_apol, rism3d_get_numAtomTypes(rism_3d), mpicomm)
    end if

    ! Allocate working memory.
    ff => safemem_realloc(ff, 3, rism3d_get_numAtoms(rism_3d))

#if defined(RISM_CRDINTERP)
    if (rismprm%fcestride > 0) &
         atomPositions_fce => safemem_realloc(atomPositions_fce, 3, rism3d_get_numAtoms(rism_3d))
#endif /* RISM_CRDINTERP */

    call temperatureDerivativeCheck()


    ! Free up a bit of memory.
    call rism3d_solvent_destroy(solvent)
    call rism3d_solute_destroy(solute)

    call flush(outunit)
    call rism_timer_stop(timer_init)

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
    ! REQUIRE( rismprm%rismnrespa /= 0 )
    if (rismprm%rismnrespa == 0) then
       call rism_report_error("Programmer error; rismprm%rismnrespa is 0 !")
    end if
    if (mod(irespa, rismprm%rismnrespa) /= 0) then
       calc_type = RISM_NONE
#ifdef RISM_CRDINTERP
    else if (rismprm%fcestride > 0 .and. fce_o%nsample >= fce_o%nbasis) then
       if (mod(irespa, rismprm%rismnrespa*rismprm%fcestride) /= 0) &
          calc_type = RISM_INTERP
#endif /*RISM_CRDINTERP*/
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
    use constants, only : KB
#ifndef CUDA    
    use rism3d_c, only : rism3d_calculateSolution
#else
    use rism3d_rism3d_c_mod
#endif /*CUDA*/
    use rism_util, only : corr_drift, alignorient, translate, findCenterOfMass, rotationalVelocity
    implicit none
#include "def_time.h"

#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/

    integer, intent(in) :: irespa
    _REAL_, intent(in) :: atomPositions_md(3, rism3d_get_numAtoms(rism_3d))
    _REAL_, intent(inout) :: frc(3, rism3d_get_numAtoms(rism_3d))
    _REAL_, intent(out) :: epol
    integer, intent(in) :: imin
    ! Solvation energy and entropy.
    _REAL_ :: epol_e, epol_ts

    !iclosure :: counter for closures
    integer :: iclosure
    integer :: i, iatu
    integer :: err
    _REAL_ :: mpi_temp, partialMolarVolume

#if !defined(SANDER)
    integer, external :: rism_calc_type
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Adding new variables

  _REAL_ :: sff(3, rism3d_get_numAtoms(rism_3d))
  _REAL_ :: ffm(3, rism3d_get_numAtoms(rism_3d))
  _REAL_ :: forcnetr

   DOUBLE PRECISION deviat

   integer jrespa,fsestride,ijrespe
   integer :: idirom=0, mmidirom=0,llidirom=0,iupdate=0

   logical :: converged
   
   save jrespa,fsestride,ijrespe
   save idirom,mmidirom,llidirom,iupdate

! Initialization of some new variables
   
   if (irespa == 0) then
      jrespa=0
   else if (irespa == 1) then
      jrespa=1
      fsestride=1
      ijrespe=0
   else
      jrespa=jrespa+1
   end if

   call rism_timer_start(timer)
   epol = 0
   ff=0
   
    !
    !Test for interpolation, minimum samples and if this is an interpolation step
    !
    if (rism_calc_type(jrespa) == RISM_NONE) then
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!No forces this steps. DO NOTHING!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(RISM_CRDINTERP)
    else if (rism_calc_type(jrespa) == RISM_INTERP) then
     !!!!!!!!!!!!!!!!
     !!!FCE forces!!!
     !!!!!!!!!!!!!!!!
       if (rismprm%verbose>=2) call rism_report_message("|LINEAR PROJECTION PREDICT!!!")

       call timer_start(TIME_REORIENT)
       atomPositions_fce = atomPositions_md
       call orient(rism_3d%solute, atomPositions_fce, rism3d_get_nsolution(rism_3d))
       call timer_stop(TIME_REORIENT)
       !linproj predict
       call timer_start(TIME_CRDINTERP)
     if (rismprm%apply_rism_force==1) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use untransformed, three versions of the transformed extrapolation, or
! the old non-normalized individual scheme

     ! old methods and associated objects (sff, sratu_fce, ratu_fce)

        ! if (fce_o%trans==0) call fce_forcea(fce_o,ff,atomPositions_fce)

        ! if (fce_o%trans==1) call fce_forceb(fce_o,ff,atomPositions_fce,sff,sratu_fce)

        ! if (fce_o%trans==2) call fce_forcebm(fce_o,ff,atomPositions_fce,sff,sratu_fce)

        ! if (fce_o%trans==3) call fce_forcec(fce_o,ff,atomPositions_fce,sff,sratu_fce)

        ! if (fce_o%trans==4) call fce_force(fce_o,ff,atomPositions_fce)

     if(fce_o%trans==0) call fce_forcea(fce_o,ff,atomPositions_fce)
        
     if(fce_o%trans==1) call fce_forceb(fce_o,ff,atomPositions_fce)

     if(fce_o%trans==2) call fce_forcebm(fce_o,ff,atomPositions_fce,iupdate,idirom)                ! ASFE

     if(fce_o%trans==3) call fce_forcebm(fce_o,ff,atomPositions_fce,iupdate,idirom)                ! ASFE
        
     if(fce_o%trans==4) call fce_force(fce_o,ff,atomPositions_fce)

     if(fce_o%trans==5) call fce_forcesa(fce_o,ff,atomPositions_fce,forcnetr,iupdate,idirom)       ! GSFE

     if(fce_o%trans==6) call fce_forcesan(fce_o,ff,atomPositions_fce,forcnetr,iupdate,idirom)      ! GSFE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     end if

       call timer_stop(TIME_CRDINTERP)

!!$     call timer_start(TIME_REORIENT)
!!$     call unorient(rism_3d%solu, atomPositions_md)
!!$     call timer_stop(TIME_REORIENT)

       if (rismprm%zerofrc==1) then
#ifdef MPI
          call corr_drift(ff, rism3d_get_mass(rism_3d), rism3d_get_numAtoms(rism_3d), &
               mpirank, mpisize, mpicomm)
#else
          call corr_drift(ff, rism3d_get_mass(rism_3d), rism3d_get_numAtoms(rism_3d))
#endif /*MPI*/
       end if
       if (rismprm%fcestride > 1) then
          if (rismprm%verbose >= 2) call rism_report_message("|IMPULSE FORCE - INTERP!!!")
          ff=rismprm%rismnrespa*ff
       end if

#endif /*RISM_CRDINTERP*/
    else
       if (rismprm%verbose >= 2) call rism_report_message("|FULL RISM!!!")
!!!!!!!!!!!!!!!!!!!!!!!!
!!!Full RISM SOLUTION!!!
!!!!!!!!!!!!!!!!!!!!!!!!
       call rism3d_setCoord(rism_3d, atomPositions_md)

#ifndef CUDA
       converged = rism3d_calculateSolution(rism_3d, rismprm%saveprogress, rismprm%progress, &
            rismprm%maxstep, rismprm%stoponfail, tolerancelist)
#else
       converged = rism3d_calculateSolution(rism_3d, rismprm%saveprogress, rismprm%progress, &
            rismprm%maxstep, rismprm%stoponfail, tolerancelist, ubound(tolerancelist,1))
#endif /*CUDA*/

       if(imin /= 0) then
          call rism_solvdist_thermo_calc(.false., 0)
       end if
#ifdef SANDER
       ! Ugly, ugly hack.  SANDER runs force on the first frame twice.
       ! This messes up the solution propagation.  Here we set the
       ! solution counter back one to ignore one of the duplicate
       ! solutions.
       if (irespa == 1) then
          call rism3d_set_nsolution(rism_3d,1)
       end if
#else /*SANDER*/
#endif /*SANDER*/
       if (rismprm%apply_rism_force == 1) then
          call rism3d_force(rism_3d, ff)
          if (rismprm%zerofrc == 1) then
#ifdef MPI
             call corr_drift(ff, rism3d_get_mass(rism_3d), rism3d_get_numAtoms(rism_3d), &
                  mpirank, mpisize, mpicomm)
#else
             call corr_drift(ff, rism3d_get_mass(rism_3d), rism3d_get_numAtoms(rism_3d))
#endif /*MPI*/
          end if
          ! Convert to [kcal/mol/A].
          ff = ff * KB * rism3d_get_temperature(rism_3d)
       end if

       ! Get the excess chemical potential.
       call timer_start(TIME_EXCESSCHEMICALPOTENTIAL)
       epol = rism_3d%excessChemicalPotential_tot(rismprm%asympCorr)*KB*rism3d_get_temperature(rism_3d)
#if 0
       ! TODO: figure out when we need to see a detail like this:
       write(6,'(a,f15.5)') 'excessChemicalPotential: ', epol
#endif
       call timer_stop(TIME_EXCESSCHEMICALPOTENTIAL)

#ifdef RISM_CRDINTERP
       if (rismprm%fcestride >0 .and. rismprm%apply_rism_force==1) then
          call timer_start(TIME_REORIENT)
          atomPositions_fce = atomPositions_md
          call orient(rism_3d%solute, atomPositions_fce, rism3d_get_nsolution(rism_3d))
          call timer_stop(TIME_REORIENT)
          call timer_start(TIME_SAVECRDINTERP)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate the accuracy of the extrapolation 

        if (fce_o%nsample >= fce_o%nbase .and. rismprm%verbose >= 1) then

        !   call fce_estimate(fce_o,ff,atomPositions_fce,sff,sratu_fce,ffm,deviat)

          call fce_estimate(fce_o,ff,atomPositions_fce,ffm,deviat,forcnetr,iupdate,idirom)

! Print the deviations using the root processor

        if (mpirank==0) then

             if(fce_o%trans==5.or.fce_o%trans==6) then
                 write(6,256) 'Error of the extrapolation:',0.5d0*deviat*100.d0, &
                      ' %    |  Net force accuracy:',forcnetr*100.d0,' %'
256              format(a31,f8.3,a28,f8.3,a2)
              else
                 write(6,257) 'Error of the extrapolation:',0.5d0*deviat*100.d0,' %'
257              format(a31,f8.3,a2)
              end if
              
              if((fce_o%trans==2.or.fce_o%trans==3.or.fce_o%trans==5.or.fce_o%trans==6).and.fce_o%ifreq>1) then
                 ! Calculating the averaged number of steps for the frequency regime
                 mmidirom=mmidirom+idirom
                 llidirom=llidirom+1
                 write(6,265) '|Averaged inverse frequency after updating:', &
                      dble(mmidirom)/llidirom
265              format(a42,f11.3)
              end if
        end if

        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update the basic vectors

        call fce_update(fce_o, ff, atomPositions_fce)

        iupdate=0

        ! Print the current number of the basic vectors and the size of the outer step
        if(mpirank == 0 .and. rismprm%verbose >= 1) then
           if(fce_o%nsample<fce_o%nbasis-1.and.fsestride<=rismprm%fcestride) then
              write(6,275) '|Number of samples:',fce_o%nsample,'  /  ', &
                   'Size of the outer time step (in dt):',fsestride*rismprm%rismnrespa
275           format(a18,i5,a5,a36,i5)
           end if
           if(fce_o%nsample == fce_o%nbasis-1) then
              write(6,276) '|Number of samples:',fce_o%nsample
            !  write(6,276) '|Number of samples:',fce_o%nsample+1
276           format(a18,i5)
           end if
           if(fsestride==rismprm%fcestride-1) then
              write(6,277) '|Size of the outer time step (in dt):',fsestride*rismprm%rismnrespa
         !     write(6,277) '|Size of the outer time step (in dt):',(fsestride+1)*rismprm%rismnrespa
277           format(a36,i5)
           end if
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   On the very beginning the exact 3D-RISM forces are applied every   !
!   rismnrespa*dt time step by fcenbase times. Then the extrapolation  !
!   starts and the 3D-RISM forces are calculated with a decreasing     !
!   frequency every 2*rismnrespa*dt, 3*rismnrespa*dt, and so on up     !
!   to fcestride*rismnrespa*dt steps. In other words, after each       !
!   outer time interval passed, it increases by 1 until achieves       !
!   the value fcestride*rismnrespa*dt.                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (fce_o%nsample >= fce_o%nbase) then

      if (ijrespe==1) then
         fsestride=rismprm%fcestride
      else
         if (fsestride < rismprm%fcestride) then

! Before the last step when the outer interval is enlarged it is
! additionally adjusted to make the total number of steps to be
! exactly multiple of fcestride*rismnrespa

            if (fsestride == rismprm%fcestride-1) then
               fsestride=rismprm%fcestride- &
                    mod(irespa,rismprm%rismnrespa*rismprm%fcestride)/rismprm%rismnrespa

               ijrespe=1

            else

               fsestride=fsestride+1

            end if
         end if
      end if

      jrespa=0

   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Perform transformation of basic force-coordinate pairs after updating

   if((fce_o%trans==1.or.fce_o%trans==2.or.fce_o%trans==3).and.fce_o%nsample>=fce_o%nbase) then
        call fce_transformi(fce_o)
   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          call timer_stop(TIME_SAVECRDINTERP)
!!$        call timer_start(TIME_REORIENT)
!!$        call unorient(rism_3d%solu, atomPositions_fce)
!!$        call timer_stop(TIME_REORIENT)
       
    endif
#endif /*RISM_CRDINTERP*/

       !         if (rismnrespa >1 .and. (.not. interpcrd>0 .or. .not. nsample >= fcenbasis)) then
       if (rismprm%rismnrespa >1) then
          if (rismprm%verbose>=2) call rism_report_message("|IMPULSE FORCE!!!")
          ff=rismprm%rismnrespa*ff
       end if

    end if

    if (rismprm%apply_rism_force==1) frc = frc + ff

#if defined(RISM_CRDINTERP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! If specified with 'fcewrite', write extrapolation basis
    
    if (mpirank == 0 .and. rismprm%fcewrite > 0) then
       if (mod(irespa,rismprm%fcewrite) == 0) call fce_wrtbasis(fce_o, jrespa)
    endif
#endif

    call flush(outunit)
    call rism_timer_stop(timer)
  end subroutine rism_force

  
  !> Calculates and stores thermodynamics for 3D-RISM, optionally
  !! printing distribution files.  Decomposition is performed based on
  !! the global parameters 'polarDecomp' and 'entropicDecomp'.  If
  !! polarDecomp==.true. then repeat the most recent calculation with solute
  !! charges turned off.  Report the usual quanities but also decompose
  !! the chemical potential into polar and non-polar terms. May be
  !! combined with entropicDecomp.  If entropicDecomp ==
  !! .true. calculate the temperature derivative of the most recent
  !! calculation.  Report the energetic and entropic components of the
  !! excess chemical potential.  May be combined with polarDecomp.
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
    use constants, only : kb, COULOMB_CONST_E
    implicit none
#ifdef MPI
    include "mpif.h"
#endif /*MPI*/
    logical*4, intent(in) :: writedist
    integer, intent(in) :: step
    integer :: err
    logical :: converged
    call rism_timer_start(timer_write)
    call rismthermo_reset(rismthermo)
    if (associated(rismthermo_pol%excessChemicalPotential)) &
         call rismthermo_reset(rismthermo_pol)
    if (associated(rismthermo_apol%excessChemicalPotential)) &
         call rismthermo_reset(rismthermo_apol)

    ! Calculate thermodynamics.
    call rism_timer_stop(timer_write)
    rismthermo%excessChemicalPotential = &
         rism_3d%excessChemicalPotential(rismprm%asympCorr)* KB * rism3d_get_temperature(rism_3d)
    if (rismprm%gfCorrection == 1) then
       rismthermo%excessChemicalPotentialGF = &
            rism_3d%excessChemicalPotentialGF(rismprm%asympCorr)* KB * rism3d_get_temperature(rism_3d)
    end if
    if (rismprm%pcplusCorrection == 1) then
       rismthermo%excessChemicalPotentialPCPLUS = &
            rism_3d%excessChemicalPotentialPCPLUS(rismprm%asympCorr)*KB*rism3d_get_temperature(rism_3d)
    end if
    if (canCalculateUC) then
       rismthermo%excessChemicalPotentialUC = &
            rism_3d%excessChemicalPotentialUC(rismprm%Uccoeff, rismprm%asympCorr) &
            * KB * rism3d_get_temperature(rism_3d)
    end if
    rismthermo%solventPotentialEnergy = rism3d_solventPotEne(rism_3d) * KB * rism3d_get_temperature(rism_3d)
    rismthermo%partialMolarVolume = rism3d_partialMolarVolume(rism_3d)
#ifndef CUDA
    rismthermo%excessParticlesBox = rism_3d%excessParticles(.false.)
#else
    rismthermo%excessParticlesBox => rism_3d%excessParticles(.false.)
#endif /*CUDA*/
    rismthermo%totalParticlesBox = rismthermo%excessParticlesBox &
         + rism3d_get_voxelVolume(rism_3d) &
         * rism3d_get_totalLocalPointsR(rism_3d) * rism_3d%get_density()
         ! this is the local volume for the MPI process.
    if (.not. rism3d_get_periodic(rism_3d)) then
#ifndef CUDA
       rismthermo%excessParticles = rism_3d%excessParticles(.true.)
#else
       rismthermo%excessParticles => rism_3d%excessParticles(.true.)
#endif /*CUDA*/
    end if
#ifndef CUDA
    rismthermo%kirkwoodBuff = rism_3d%kirkwoodbuff(rismprm%asympCorr)
    rismthermo%DCFintegral = rism3d_DCFintegral(rism_3d)
#else
    rismthermo%kirkwoodBuff => rism_3d%kirkwoodbuff(rismprm%asympCorr)
    rismthermo%DCFintegral => rism3d_DCFintegral(rism_3d)
#endif /*CUDA*/
    if (rismprm%entropicDecomp == 1 .and. rism3d_canCalc_DT(rism_3d)) then
       converged =  rism3d_calculateSolution_dT(rism_3d, rismprm%progress, &
            rismprm%maxstep, rismprm%stoponfail, tolerancelist(size(tolerancelist)))
       rismthermo%solvationEnergy = &
            rism_3d%solvationEnergy(rismprm%asympCorr) * KB * rism3d_get_temperature(rism_3d)
       if (rismprm%gfCorrection == 1) then
          rismthermo%solvationEnergyGF = &
               rism_3d%solvationEnergyGF(rismprm%asympCorr) * KB * rism3d_get_temperature(rism_3d)
       end if
       rismthermo%excessParticles_dT = rism_3d%excessParticles_dT(rismprm%asympCorr) &
            / rism3d_get_temperature(rism_3d)
       rismthermo%kirkwoodBuff_dT = &
            rism_3d%kirkwoodbuff_dT(rismprm%asympCorr) &
            / rism3d_get_temperature(rism_3d)
       rismthermo%DCFintegral_dT = rism3d_DCFintegral_dT(rism_3d) &
            / rism3d_get_temperature(rism_3d)
       ! PMV_dT based quantities
       if (canCalculatePartialMolarVolume_dT) then
          rismthermo%solvationEnergyPCPLUS = &
               rism_3d%solvationEnergyPCPLUS(rismprm%asympCorr)*KB*rism3d_get_temperature(rism_3d)
          if (canCalculateUC) then
             rismthermo%solvationEnergyUC = &
                  rism_3d%solvationEnergyUC(rismprm%Uccoeff, rismprm%asympCorr)&
                  * KB * rism3d_get_temperature(rism_3d)
          end if
          rismthermo%partialMolarVolume_dT = rism3d_partialMolarVolume_dT(rism_3d) &
               / rism3d_get_temperature(rism_3d)
       end if
    end if

    ! Output distributions.
    if (writedist .and. step >= 0) &
         call rism_writeVolumetricData(rism_3d, step)
    call rism_timer_start(timer_write)

    call rismthermo_mpireduce(rismthermo)

    ! By doing the second RISM calculation here output from the RISM
    ! routines does not break up the thermodynamic output below.
    if (rismprm%polarDecomp == 1) then
       ! Redo calculation with charges off.
       call rism_timer_stop(timer_write)
       call rism3d_unsetCharges(rism_3d)

#ifndef CUDA
       converged = rism3d_calculateSolution(rism_3d, rismprm%saveprogress, rismprm%progress, &
            rismprm%maxstep, rismprm%stoponfail, tolerancelist)
#else
       converged = rism3d_calculateSolution(rism_3d, rismprm%saveprogress, rismprm%progress, &
            rismprm%maxstep, rismprm%stoponfail, tolerancelist, ubound(tolerancelist,1))
#endif /*CUDA*/

       ! Calculate thermodynamics.
       rismthermo_apol%excessChemicalPotential = &
            rism_3d%excessChemicalPotential(rismprm%asympCorr)*KB*rism3d_get_temperature(rism_3d)
       if (rismprm%gfCorrection == 1) then
          rismthermo_apol%excessChemicalPotentialGF = &
               rism_3d%excessChemicalPotentialGF(rismprm%asympCorr)*KB*rism3d_get_temperature(rism_3d)
       end if
       if (rismprm%pcplusCorrection == 1) then
          rismthermo_apol%excessChemicalPotentialPCPLUS = &
               rism_3d%excessChemicalPotentialPCPLUS(rismprm%asympCorr)*KB*rism3d_get_temperature(rism_3d)
       end if
       if (canCalculateUC) then 
            rismthermo_apol%excessChemicalPotentialUC = &
            rism_3d%excessChemicalPotentialUC(rismprm%Uccoeff, rismprm%asympCorr) &
            *KB*rism3d_get_temperature(rism_3d)
       end if
       rismthermo_apol%solventPotentialEnergy = rism3d_solventPotEne(rism_3d)*KB*rism3d_get_temperature(rism_3d)
       rismthermo_apol%partialMolarVolume = rism3d_partialMolarVolume(rism_3d)
#ifndef CUDA
       rismthermo_apol%excessParticles = rism_3d%excessParticles(rismprm%asympCorr)
       rismthermo_apol%DCFintegral = rism3d_DCFintegral(rism_3d)
       rismthermo_apol%kirkwoodBuff = rism_3d%kirkwoodbuff(rismprm%asympCorr)
#else
       rismthermo_apol%excessParticles => rism_3d%excessParticles(rismprm%asympCorr)
       rismthermo_apol%DCFintegral => rism3d_DCFintegral(rism_3d)
       rismthermo_apol%kirkwoodBuff => rism_3d%kirkwoodbuff(rismprm%asympCorr)
#endif /*CUDA*/
       if (rismprm%entropicDecomp == 1 .and. rism3d_canCalc_DT(rism_3d)) then
          converged = rism3d_calculateSolution_dT(rism_3d, rismprm%progress, &
               rismprm%maxstep, rismprm%stoponfail, tolerancelist(size(tolerancelist)))
          rismthermo_apol%solvationEnergy = &
               rism_3d%solvationEnergy(rismprm%asympCorr)*KB*rism3d_get_temperature(rism_3d)
          if (rismprm%gfCorrection == 1) then
             rismthermo_apol%solvationEnergyGF = &
                  rism_3d%solvationEnergyGF(rismprm%asympCorr)*KB*rism3d_get_temperature(rism_3d)
          end if
          rismthermo_apol%excessParticles_dT = rism_3d%excessParticles_dT(rismprm%asympCorr) &
               /rism3d_get_temperature(rism_3d)
          rismthermo_apol%kirkwoodBuff_dT = &
               rism_3d%kirkwoodbuff_dT(rismprm%asympCorr) &
               /rism3d_get_temperature(rism_3d)
          rismthermo_apol%DCFintegral_dT = rism3d_DCFintegral_dT(rism_3d) &
               /rism3d_get_temperature(rism_3d)
          if (canCalculatePartialMolarVolume_dT) then
             rismthermo_apol%solvationEnergyPCPLUS = &
                  rism_3d%solvationEnergyPCPLUS(rismprm%asympCorr)*KB*rism3d_get_temperature(rism_3d)
             if (canCalculateUC) then
                rismthermo_apol%solvationEnergyUC = &
                     rism_3d%solvationEnergyUC(rismprm%Uccoeff, rismprm%asympCorr)&
                     *KB*rism3d_get_temperature(rism_3d)
             end if
             rismthermo_apol%partialMolarVolume_dT= &
                  rism3d_partialMolarVolume_dT(rism_3d)/rism3d_get_temperature(rism_3d)
          end if
       end if
       call rism3d_resetCharges(rism_3d)
       call rism_timer_start(timer_write)
       call rismthermo_mpireduce(rismthermo_apol)

       ! Calculate the polar contribution.
       call rismthermo_sub(rismthermo_pol, rismthermo, rismthermo_apol)
    end if
    if (rismprm%selftest) then
       call rism_timer_stop(timer_write)
       call rism3d_selftest(rism_3d)
       call rism_timer_start(timer_write)
    end if
    call rism_timer_stop(timer_write)
  end subroutine rism_solvdist_thermo_calc

  
  !> The format for NAB programs will consist of one line for each
  !! category of data: energy, volume and excess.  RISM-specific values
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
    use constants, only : kb, COULOMB_CONST_E, PI
    use rism_report_c
    implicit none
    logical*4, intent(in) :: description
#ifdef SANDER
    _REAL_, intent(in) :: soluPot(25)
#else
    _REAL_, intent(in) :: soluPot(11)
#endif
    
    integer :: iv, iu, ig, il, err
    integer :: ix, iy, iz

    character(len = 4) :: atomName(rism3d_get_numAtomTypes(rism_3d))
    call rism3d_get_atomName(rism_3d, atomName)

    call rism_timer_start(timer_write)
    if(rism3d_get_failure(rism_3d)) then
       call rism_report_warn('SKIPPING THERMODYNAMIC OUTPUT DUE TO CONVERGENCE FAILURE.')
       return
    end if

!     call rism_timer_start(timer_write)
    
    if (description) then
       if (mpirank==0) then

          ! Thermodynamics key.

          ! Free energy-type properties.
          call thermo_print_descr_line('solutePotentialEnergy', '[kcal/mol]', 'Total', '', &
               (/'LJ        ', 'Coulomb   ', 'Bond      ', 'Angle     ', 'Dihedral  ', &
               'H-Bond    ', 'LJ-14     ', 'Coulomb-14', 'Restraints', '3D-RISM   '/), 10)
          call thermo_print_descr_line('rism_excessChemicalPotential', '[kcal/mol]', 'Total', 'ExcChemPot_', &
               atomName, rism3d_get_numAtomTypes(rism_3d))
          if (rismprm%gfCorrection == 1) then
             call thermo_print_descr_line('rism_excessChemicalPotentialGF', '[kcal/mol]', &
                  'Total', 'ExcChemPot_GF_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
          end if
          if (rismprm%pcplusCorrection == 1) then                       
             call thermo_print_descr_line('rism_excessChemicalPotentialPCPLUS','[kcal/mol]',&
                  'Total','ExcChemPot_PCPLUS_',&
                  atomName,0)
          end if
          if (canCalculateUC) then
             call thermo_print_descr_line('rism_excessChemicalPotentialUC', &
                  '[kcal/mol]', 'Total', 'ExcChemPot_UC_', &
                  atomName, 0)
          end if
          call thermo_print_descr_line('rism_solventPotentialEnergy', '[kcal/mol]', 'Total', 'UV_potential_', &
               atomName, rism3d_get_numAtomTypes(rism_3d))

          ! Thermodynamic properties not related to free energy.
          call thermo_print_descr_line('rism_partialMolarVolume', '[A^3]', 'PMV', '', &
               atomName, 0)
          call thermo_print_descr_line('rism_totalParticlesBox', '[#]', '', 'TotalNum_', &
               atomName, rism3d_get_numAtomTypes(rism_3d))
          call thermo_print_descr_line('rism_totalChargeBox', '[e]', 'Total', 'TotalChg_', &
               atomName, rism3d_get_numAtomTypes(rism_3d))
          call thermo_print_descr_line('rism_excessParticlesBox', '[#]', '', 'ExNum_', &
               atomName, rism3d_get_numAtomTypes(rism_3d))
          call thermo_print_descr_line('rism_excessChargeBox', '[e]', 'Total', 'ExChg_', &
               atomName, rism3d_get_numAtomTypes(rism_3d))
          if (.not. rism3d_get_periodic(rism_3d)) then
             call thermo_print_descr_line('rism_excessParticles', '[#]', '', 'ExNum_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_excessCharge', '[e]', 'Total', 'ExChg_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
          end if
          call thermo_print_descr_line('rism_KirkwoodBuff', '[A^3]', '', 'KB_', &
               atomName, rism3d_get_numAtomTypes(rism_3d))
          call thermo_print_descr_line('rism_DCFintegral', '[A^3]', '', 'DCFI_', &
               atomName, rism3d_get_numAtomTypes(rism_3d))

          ! Entropic decomposition.
          if (rismprm%entropicDecomp == 1 .and. rism3d_canCalc_DT(rism_3d)) then
             call thermo_print_descr_line('rism_-TS', '[kcal/mol]', 'Total', '-TS_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             if (rismprm%gfCorrection == 1) then             
                call thermo_print_descr_line('rism_-TS_GF', '[kcal/mol]', 'Total', '-TS_GF_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
             end if
             if (canCalculatePartialMolarVolume_dT) then
                if (rismprm%pcplusCorrection == 1) then
                   call thermo_print_descr_line('rism_-TS_PCPLUS','[kcal/mol]','Total','-TS_PCPLUS_',&
                        atomName,0)
                end if
                if (canCalculateUC) then
                   call thermo_print_descr_line('rism_-TS_UC', '[kcal/mol]', 'Total', '', &
                        atomName, 0)
                end if
             end if
             call thermo_print_descr_line('rism_solvationEnergy', '[kcal/mol]', 'Total', 'SolvEne_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             if (rismprm%gfCorrection == 1) then
                call thermo_print_descr_line('rism_solvationEnergyGF', '[kcal/mol]', 'Total', 'SolEne_GF_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
             end if
             if (canCalculatePartialMolarVolume_dT) then
                if (rismprm%pcplusCorrection == 1) then
                   call thermo_print_descr_line('rism_solvationEnergyPCPLUS','[kcal/mol]',&
                        'Total','SolvEne_PCPLUS_',&
                        atomName,0)
                end if
                if (canCalculateUC) then 
                   call thermo_print_descr_line('rism_solvationEnergyUC','[kcal/mol]','Total','SolvEne_UC_',&
                        atomName,0)
                end if
                call thermo_print_descr_line('rism_partialMolarVolume_dT', '[A^3/K]', 'PMV_dT', '', &
                     atomName, 0)
             end if
             call thermo_print_descr_line('rism_excessParticles_dT', '[#/K]', '', 'ExNum_dT_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_excessCharge_dT', '[e/K]', 'Total', 'ExChg_dT_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_KirkwoodBuff_dT', '[A^3/K]', '', 'KB_dT_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_DCFintegral_dT', '[A^3/K]', '', 'DCFI_dT_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
          end if

          ! Polar decomposition.
          if (rismprm%polarDecomp==1) then
             call thermo_print_descr_line('rism_polarExcessChemicalPotential', '[kcal/mol]', 'Total', 'polar_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_apolarExcessChemicalPotential', '[kcal/mol]', 'Total', 'apolar_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             if (rismprm%gfCorrection == 1) then
                call thermo_print_descr_line('rism_polarExcessChemicalPotentialGF', '[kcal/mol]', &
                     'Total', 'polar_GF_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_descr_line('rism_apolarExcessChemicalPotentialGF', '[kcal/mol]', &
                     'Total', 'apolar_GF_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
             end if
             if (rismprm%pcplusCorrection == 1) then
             call thermo_print_descr_line('rism_polarExcessChemicalPotentialPCPLUS','[kcal/mol]','Total','polar_PCPLUS_',&
                  atomName,0)
             call thermo_print_descr_line('rism_apolarExcessChemicalPotentialPCPLUS','[kcal/mol]','Total','apolar_PCPLUS_',&
                  atomName,0)
             end if
             if (canCalculateUC) then
                call thermo_print_descr_line('rism_polarExcessChemicalPotentialUC', &
                     '[kcal/mol]', 'Total', '', &
                     atomName, 0)
                call thermo_print_descr_line('rism_apolarExcessChemicalPotentialUC', &
                     '[kcal/mol]', 'Total', '', &
                     atomName, 0)
             end if

             ! Solute-solvent potential energy.
             call thermo_print_descr_line('rism_polarSolventPotentialEnergy', '[kcal/mol]', 'Total', 'polar_UV_pot_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_apolarSolventPotentialEnergy', '[kcal/mol]', 'Total', 'apolar_UV_pot_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_polarPartialMolarVolume', '[A^3]', 'polar_PMV', '', &
                  atomName, 0)
             call thermo_print_descr_line('rism_apolarPartialMolarVolume', '[A^3]', 'apolar_PMV', '', &
                  atomName, 0)
             call thermo_print_descr_line('rism_polarExcessParticles', '[#]', '', 'polar_ExNum_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_apolarExcessParticles', '[#]', '', 'apolar_ExNum_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_polarExcessCharge', '[e]', 'Total', 'polar_ExChg_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_apolarExcessCharge', '[e]', 'Total', 'apolar_ExChg_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_polarKirkwoodBuff', '[A^3]', '', 'polar_KB_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_apolarKirkwoodBuff', '[A^3]', '', 'apolar_KB_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_polarDCFintegral', '[A^3]', '', 'polar_DCFI_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_descr_line('rism_apolarDCFintegral', '[A^3]', '', 'apolar_DCFI_', &
                  atomName, rism3d_get_numAtomTypes(rism_3d))

             ! Entropic decomposition w/ polar decomposition.
             if (rismprm%entropicDecomp == 1 .and. rism3d_canCalc_DT(rism_3d)) then
                call thermo_print_descr_line('rism_polar-TS', '[kcal/mol]', 'Total', 'polar_-TS_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_descr_line('rism_apolar-TS', '[kcal/mol]', 'Total', 'apolar_-TS_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                if (rismprm%gfCorrection == 1) then
                   call thermo_print_descr_line('rism_polar-TS_GF', '[kcal/mol]', 'Total', 'polar_-TS_GF_', &
                        atomName, rism3d_get_numAtomTypes(rism_3d))
                   call thermo_print_descr_line('rism_apolar-TS_GF', '[kcal/mol]', 'Total', 'apolar_-TS_GF_', &
                        atomName, rism3d_get_numAtomTypes(rism_3d))
                end if
                if (canCalculatePartialMolarVolume_dT) then 
                   if (rismprm%pcplusCorrection == 1) then
                      call thermo_print_descr_line('rism_polar-TS_PCPLUS','[kcal/mol]','Total','polar_-TS_PCPLUS_',&
                           atomName,0)
                      call thermo_print_descr_line('rism_apolar-TS_PCPLUS','[kcal/mol]','Total','apolar_-TS_PCPLUS_',&
                           atomName,0)
                   end if
                   if (canCalculateUC) then
                      call thermo_print_descr_line('rism_polar-TS_UC', '[kcal/mol]', 'Total', '', &
                           atomName, 0)
                      call thermo_print_descr_line('rism_apolar-TS_UC', '[kcal/mol]', 'Total', '', &
                           atomName, 0)
                   end if
                end if
                call thermo_print_descr_line('rism_polarSolvationEnergy', '[kcal/mol]', 'Total', 'polar_solv_ene_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_descr_line('rism_apolarSolvationEnergy', '[kcal/mol]', 'Total', 'apolar_solv_ene_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                if (rismprm%gfCorrection == 1) then
                   call thermo_print_descr_line('rism_polarSolvationEnergyGF', '[kcal/mol]', &
                        'Total', 'polar_solv_ene_GF_', &
                        atomName, rism3d_get_numAtomTypes(rism_3d))
                   call thermo_print_descr_line('rism_apolarSolvationEnergyGF', '[kcal/mol]', &
                        'Total', 'apolar_solv_ene_GF_', &
                        atomName, rism3d_get_numAtomTypes(rism_3d))
                end if
                if (canCalculatePartialMolarVolume_dT) then
                   if (rismprm%pcplusCorrection == 1) then
                      call thermo_print_descr_line('rism_polarSolvationEnergyPCPLUS',&
                           '[kcal/mol]','Total','polar_solv_ene_PCPLUS_',&
                           atomName,0)
                      call thermo_print_descr_line('rism_apolarSolvationEnergyPCPLUS',&
                           '[kcal/mol]','Total','apolar_solv_ene_PCPLUS_',&
                           atomName,0)
                   end if
                   if (canCalculateUC) then
                      call thermo_print_descr_line('rism_polarSolvationEnergyUC',&
                           '[kcal/mol]','Total','polar_solv_ene_UC_',&
                           atomName,0)
                      call thermo_print_descr_line('rism_apolarSolvationEnergyUC',&
                           '[kcal/mol]','Total','apolar_solv_ene_UC_',&
                           atomName,0)
                   end if
                end if
                ! Solute-solvent potential energy.
                call thermo_print_descr_line('rism_polarSolventPotentialEnergy_dT', '[kcal/mol]', 'Total', 'polar_UV_pot_dT_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_descr_line('rism_apolarSolventPotentialEnergy_dT', '[kcal/mol]', 'Total', 'apolar_UV_pot_dT_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                if (canCalculatePartialMolarVolume_dT) then
                   call thermo_print_descr_line('rism_polarPartialMolarVolume_dT', &
                        '[A^3]', 'polar_partialMolarVolume_dT', '', &
                        atomName, 0)
                   call thermo_print_descr_line('rism_apolarPartialMolarVolume_dT', &
                        '[A^3]', 'apolar_partialMolarVolume_dT', '', &
                        atomName, 0)
                end if
                call thermo_print_descr_line('rism_polarExcessParticles_dT', '[#]', '', 'polar_ExNum_dT_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_descr_line('rism_apolarExcessParticles_dT', '[#]', '', 'apolar_ExNum_dT_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_descr_line('rism_polarExcessCharge_dT', '[e]', 'Total', 'polar_ExChg_dT_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_descr_line('rism_apolarExcessCharge_dT', '[e]', 'Total', 'apolar_ExChg_dT_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_descr_line('rism_polarKirkwoodBuff_dT', '[A^3]', '', 'polar_KB_dT_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_descr_line('rism_apolarKirkwoodBuff_dT', '[A^3]', '', 'apolar_KB_dT_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_descr_line('rism_polarDCFintegral_dT', '[A^3]', '', 'polar_DCFI_dT_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_descr_line('rism_apolarDCFintegral_dT', '[A^3]', '', 'apolar_DCFI_dT_', &
                     atomName, rism3d_get_numAtomTypes(rism_3d))
             end if
          end if
       end if
    else
       ! DATA: free energy-based properties.
       if (mpirank==0 .and. associated(rismthermo%excessChemicalPotential)) then

#ifdef SANDER
          ! Cast the SANDER array into the order of the NAB array.
          call thermo_print_results_line('solutePotentialEnergy', soluPot(1), (/soluPot(2), &
               soluPot(3), soluPot(5), soluPot(6), soluPot(7), soluPot(13), &
               soluPot(8), soluPot(9), soluPot(10), soluPot(24)/), 10)
#else
          call thermo_print_results_line('solutePotentialEnergy', soluPot(1), soluPot(2:), 10)
#endif
          call thermo_print_results_line('rism_excessChemicalPotential', sum(rismthermo%excessChemicalPotential), &
               rismthermo%excessChemicalPotential, rism3d_get_numAtomTypes(rism_3d))
          if (rismprm%gfCorrection == 1) then
             call thermo_print_results_line('rism_excessChemicalPotentialGF', sum(rismthermo%excessChemicalPotentialGF), &
                  rismthermo%excessChemicalPotentialGF, rism3d_get_numAtomTypes(rism_3d))
          end if
          if (rismprm%pcplusCorrection == 1) then
             call thermo_print_results_line('rism_excessChemicalPotentialPCPLUS',rismthermo%excessChemicalPotentialPCPLUS,&
                  (/1d0/),0)
          end if
          if (canCalculateUC) then
             call thermo_print_results_line('rism_excessChemicalPotentialUC', rismthermo%excessChemicalPotentialUC, &
                  (/1d0/), 0)
          end if
          call thermo_print_results_line('rism_solventPotentialEnergy', &
               sum(rismthermo%solventPotentialEnergy), &
               rismthermo%solventPotentialEnergy, &
               rism3d_get_numAtomTypes(rism_3d))

          ! Non-free energy-based properties.
          call thermo_print_results_line('rism_partialMolarVolume', &
               rismthermo%partialMolarVolume, (/1d0/), 0)
          call thermo_print_results_line('rism_totalParticlesBox', &
               HUGE(1d0), rismthermo%totalParticlesBox, &
               rism3d_get_numAtomTypes(rism_3d))
          call thermo_print_results_line('rism_totalChargeBox', &
               sum(rismthermo%totalParticlesBox * rism_3d%get_charge()) &
               * sqrt((KB * rism3d_get_temperature(rism_3d)) / COULOMB_CONST_E), &
               rismthermo%totalParticlesBox * rism_3d%get_charge() &
               * sqrt((KB * rism3d_get_temperature(rism_3d)) / COULOMB_CONST_E), &
               rism3d_get_numAtomTypes(rism_3d))
          call thermo_print_results_line('rism_excessParticlesBox', &
               HUGE(1d0), rismthermo%excessParticlesBox, rism3d_get_numAtomTypes(rism_3d))
          call thermo_print_results_line('rism_excessChargeBox', &
               sum(rismthermo%excessParticlesBox * rism_3d%get_charge()) &
               * sqrt((KB * rism3d_get_temperature(rism_3d)) / COULOMB_CONST_E), &
               rismthermo%excessParticlesBox * rism_3d%get_charge() &
               * sqrt((KB * rism3d_get_temperature(rism_3d)) / COULOMB_CONST_E), &
               rism3d_get_numAtomTypes(rism_3d))
          if (.not. rism3d_get_periodic(rism_3d)) then
             call thermo_print_results_line('rism_excessParticles', &
                  HUGE(1d0), rismthermo%excessParticles, &
                  rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_excessCharge', &
                  sum(rismthermo%excessParticles * rism_3d%get_charge()) &
                  * sqrt((KB * rism3d_get_temperature(rism_3d)) / COULOMB_CONST_E), &
                  rismthermo%excessParticles * rism_3d%get_charge() &
                  * sqrt((KB * rism3d_get_temperature(rism_3d)) / COULOMB_CONST_E), &
                  rism3d_get_numAtomTypes(rism_3d))
          end if
          call thermo_print_results_line('rism_KirkwoodBuff', HUGE(1d0), &
               rismthermo%kirkwoodBuff, rism3d_get_numAtomTypes(rism_3d))
          call thermo_print_results_line('rism_DCFintegral', HUGE(1d0), &
               rismthermo%DCFintegral, rism3d_get_numAtomTypes(rism_3d))

          ! Entropic decomposition.
          if (rismprm%entropicDecomp == 1 .and. rism3d_canCalc_DT(rism_3d)) then
             call thermo_print_results_line('rism_-TS', &
                  -sum(rismthermo%solvationEnergy - rismthermo%excessChemicalPotential), &
                  -1d0*(rismthermo%solvationEnergy - rismthermo%excessChemicalPotential), rism3d_get_numAtomTypes(rism_3d))
             if (rismprm%gfCorrection == 1) then
                call thermo_print_results_line('rism_-TS_GF', &
                     -sum(rismthermo%solvationEnergyGF - rismthermo%excessChemicalPotentialGF), &
                     -1d0*(rismthermo%solvationEnergyGF - rismthermo%excessChemicalPotentialGF), &
                     rism3d_get_numAtomTypes(rism_3d))
             end if
             ! special case. If the solvation energy is HUGE() the
             ! expression below will not be and will not be caught by
             ! thermo_print_results_line
             if (canCalculatePartialMolarVolume_dT) then
                if (rismprm%pcplusCorrection == 1) then
                   call thermo_print_results_line('rism_-TS_PCPLUS', &
                        -(rismthermo%solvationEnergyPCPLUS - rismthermo%excessChemicalPotentialPCPLUS),&
                        (/-1d0*(rismthermo%solvationEnergyPCPLUS - rismthermo%excessChemicalPotentialPCPLUS)/),0)
                end if
                if (canCalculateUC) then
                   call thermo_print_results_line('rism_-TS_UC', &
                        -1d0*(rismthermo%solvationEnergyUC - rismthermo%excessChemicalPotentialUC), &
                        (/rismthermo%solvationEnergyUC/), 0)
                end if
             end if
             call thermo_print_results_line('rism_solvationEnergy', &
                  sum(rismthermo%solvationEnergy), rismthermo%solvationEnergy, rism3d_get_numAtomTypes(rism_3d))
             if (rismprm%gfCorrection == 1) then
                call thermo_print_results_line('rism_solvationEnergyGF', &
                     sum(rismthermo%solvationEnergyGF), rismthermo%solvationEnergyGF, rism3d_get_numAtomTypes(rism_3d))
             end if
             if (rismprm%pcplusCorrection == 1) then
                call thermo_print_results_line('rism_solvationEnergyPCPLUS', &
                     rismthermo%solvationEnergyPCPLUS,(/rismthermo%solvationEnergyPCPLUS/),0)
             end if
             if (canCalculateUC) then
                call thermo_print_results_line('rism_solvationEnergyUC', &
                     rismthermo%solvationEnergyUC, (/rismthermo%solvationEnergyUC/), 0)
             end if

             ! Non-free energy-based properties.
             call thermo_print_results_line('rism_partialMolarVolume_dT', &
                  rismthermo%partialMolarVolume_dT, (/rismthermo%partialMolarVolume_dT/), 0)
             call thermo_print_results_line('rism_excessParticles_dT', &
                  HUGE(1d0), rismthermo%excessParticles_dT, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_excessCharge_dT', &
                  sum(rismthermo%excessParticles_dT*rism_3d%get_charge()) &
                  *sqrt((KB *rism3d_get_temperature(rism_3d))/COULOMB_CONST_E), &
                  rismthermo%excessParticles_dT*rism_3d%get_charge()&
                  *sqrt((KB *rism3d_get_temperature(rism_3d))/COULOMB_CONST_E), &
                  rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_KirkwoodBuff_dT', HUGE(1d0), &
                  rismthermo%kirkwoodBuff_dT, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_DCFintegral_dT', HUGE(1d0), &
                  rismthermo%DCFintegral_dT, rism3d_get_numAtomTypes(rism_3d))
          end if

          ! Polar decomposition.
          if (rismprm%polarDecomp == 1) then
             call thermo_print_results_line('rism_polarExcessChemicalPotential', &
                  sum(rismthermo_pol%excessChemicalPotential), &
                  rismthermo_pol%excessChemicalPotential, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_apolarExcessChemicalPotential', &
                  sum(rismthermo_apol%excessChemicalPotential), &
                  rismthermo_apol%excessChemicalPotential, rism3d_get_numAtomTypes(rism_3d))
             if (rismprm%gfCorrection == 1) then
                call thermo_print_results_line('rism_polarExcessChemicalPotentialGF', &
                     sum(rismthermo_pol%excessChemicalPotentialGF), &
                     rismthermo_pol%excessChemicalPotentialGF, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_results_line('rism_apolarExcessChemicalPotentialGF', &
                     sum(rismthermo_apol%excessChemicalPotentialGF), &
                     rismthermo_apol%excessChemicalPotentialGF, rism3d_get_numAtomTypes(rism_3d))
             end if
             if (rismprm%pcplusCorrection == 1) then
                call thermo_print_results_line('rism_polarExcessChemicalPotentialPCPLUS',&
                     rismthermo_pol%excessChemicalPotentialPCPLUS,&
                     (/rismthermo_pol%excessChemicalPotentialPCPLUS/),0)
                call thermo_print_results_line('rism_apolarExcessChemicalPotentialPCPLUS',&
                     rismthermo_apol%excessChemicalPotentialPCPLUS,&
                     (/rismthermo_apol%excessChemicalPotentialPCPLUS/),0)
             end if
             if (canCalculateUC) then
                call thermo_print_results_line('rism_polarExcessChemicalPotentialUC', rismthermo_pol%excessChemicalPotentialUC, &
                     (/1d0/), 0)
                call thermo_print_results_line('rism_apolarExcessChemicalPotentialUC', rismthermo_apol%excessChemicalPotentialUC, &
                     (/1d0/), 0)
             end if
             call thermo_print_results_line('rism_polarSolventPotentialEnergy', &
                  sum(rismthermo_pol%solventPotentialEnergy), rismthermo_pol%solventPotentialEnergy, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_apolarSolventPotentialEnergy', &
                  sum(rismthermo_apol%solventPotentialEnergy), rismthermo_apol%solventPotentialEnergy, rism3d_get_numAtomTypes(rism_3d))

             ! Non-free energy-based properties.
             call thermo_print_results_line('rism_polarPartialMolarVolume', &
                  rismthermo_pol%partialMolarVolume, (/1d0/), 0)
             call thermo_print_results_line('rism_apolarPartialMolarVolume', &
                  rismthermo_apol%partialMolarVolume, (/1d0/), 0)
             call thermo_print_results_line('rism_polarExcessParticles', &
                  HUGE(1d0), rismthermo_pol%excessParticles, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_apolarExcessParticles', &
                  HUGE(1d0), rismthermo_apol%excessParticles, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_polarExcessCharge', &
                  sum(rismthermo_pol%excessParticles*rism_3d%get_charge()) &
                  *sqrt((KB *rism3d_get_temperature(rism_3d))/COULOMB_CONST_E), &
                  rismthermo_pol%excessParticles*rism_3d%get_charge()&
                  *sqrt((KB *rism3d_get_temperature(rism_3d))/COULOMB_CONST_E), &
                  rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_apolarExcessCharge', &
                  sum(rismthermo_apol%excessParticles*rism_3d%get_charge()) &
                  *sqrt((KB *rism3d_get_temperature(rism_3d))/COULOMB_CONST_E), &
                  rismthermo_apol%excessParticles*rism_3d%get_charge()&
                  *sqrt((KB *rism3d_get_temperature(rism_3d))/COULOMB_CONST_E), &
                  rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_polarKirkwoodBuff', HUGE(1d0), &
                  rismthermo_pol%kirkwoodBuff, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_apolarKirkwoodBuff', HUGE(1d0), &
                  rismthermo_apol%kirkwoodBuff, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_polarDCFintegral', HUGE(1d0), &
                  rismthermo_pol%DCFintegral, rism3d_get_numAtomTypes(rism_3d))
             call thermo_print_results_line('rism_apolarDCFintegral', HUGE(1d0), &
                  rismthermo_apol%DCFintegral, rism3d_get_numAtomTypes(rism_3d))
             ! Entropic decomposition w/ polar decomposition.
             if (rismprm%entropicDecomp == 1 .and. rism3d_canCalc_DT(rism_3d)) then
                call thermo_print_results_line('rism_polar-TS', &
                     -sum(rismthermo_pol%solvationEnergy - rismthermo_pol%excessChemicalPotential), &
                     -1d0*(rismthermo_pol%solvationEnergy - rismthermo_pol%excessChemicalPotential), rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_results_line('rism_apolar-TS', &
                     -sum(rismthermo_apol%solvationEnergy - rismthermo_apol%excessChemicalPotential), &
                     -1d0*(rismthermo_apol%solvationEnergy - rismthermo_apol%excessChemicalPotential), rism3d_get_numAtomTypes(rism_3d))
                if (rismprm%gfCorrection == 1) then
                   call thermo_print_results_line('rism_polar-TS_GF', &
                        -sum(rismthermo_pol%solvationEnergyGF - rismthermo_pol%excessChemicalPotentialGF), &
                        -1d0*(rismthermo_pol%solvationEnergyGF - rismthermo_pol%excessChemicalPotentialGF), &
                        rism3d_get_numAtomTypes(rism_3d))
                   call thermo_print_results_line('rism_apolar-TS_GF', &
                        -sum(rismthermo_apol%solvationEnergyGF - rismthermo_apol%excessChemicalPotentialGF), &
                        -1d0*(rismthermo_apol%solvationEnergyGF - rismthermo_apol%excessChemicalPotentialGF), &
                        rism3d_get_numAtomTypes(rism_3d))
                end if
                ! special case. If the solvation energy is HUGE() the
                ! expression below will not be and will not be caught by
                ! thermo_print_results_line
                if (canCalculatePartialMolarVolume_dT) then
                   if (rismprm%pcplusCorrection == 1) then
                      call thermo_print_results_line('rism_polar-TS_PCPLUS', &
                           -(rismthermo_pol%solvationEnergyPCPLUS - rismthermo_pol%excessChemicalPotentialPCPLUS),&
                           (/-1d0*(rismthermo_pol%solvationEnergyPCPLUS - rismthermo_pol%excessChemicalPotentialPCPLUS)/),0)
                      call thermo_print_results_line('rism_apolar-TS_PCPLUS', &
                           -(rismthermo_apol%solvationEnergyPCPLUS - rismthermo_apol%excessChemicalPotentialPCPLUS),&
                           (/-1d0*(rismthermo_apol%solvationEnergyPCPLUS - rismthermo_apol%excessChemicalPotentialPCPLUS)/),0)
                   end if
                   if (canCalculateUC) then
                      call thermo_print_results_line('rism_polar-TS_UC', &
                           -1d0*(rismthermo_pol%solvationEnergyUC - rismthermo_pol%excessChemicalPotentialUC), &
                           (/rismthermo_pol%solvationEnergyUC/), 0)
                      call thermo_print_results_line('rism_apolar-TS_UC', &
                           -1d0*(rismthermo_apol%solvationEnergyUC - rismthermo_apol%excessChemicalPotentialUC), &
                           (/rismthermo_apol%solvationEnergyUC/), 0)
                   end if
                end if
                call thermo_print_results_line('rism_polarSolvationEnergy', &
                     sum(rismthermo_pol%solvationEnergy), rismthermo_pol%solvationEnergy, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_results_line('rism_apolarSolvationEnergy', &
                     sum(rismthermo_apol%solvationEnergy), rismthermo_apol%solvationEnergy, rism3d_get_numAtomTypes(rism_3d))
                if (rismprm%gfCorrection == 1) then
                   call thermo_print_results_line('rism_polarSolvationEnergyGF', &
                        sum(rismthermo_pol%solvationEnergyGF), rismthermo_pol%solvationEnergyGF, rism3d_get_numAtomTypes(rism_3d))
                   call thermo_print_results_line('rism_apolarSolvationEnergyGF', &
                        sum(rismthermo_apol%solvationEnergyGF), rismthermo_apol%solvationEnergyGF, rism3d_get_numAtomTypes(rism_3d))
                end if
                if (rismprm%pcplusCorrection == 1) then
                   call thermo_print_results_line('rism_polarSolvationEnergyPCPLUS', &
                        rismthermo_pol%solvationEnergyPCPLUS,(/rismthermo_pol%solvationEnergyPCPLUS/),0)
                   call thermo_print_results_line('rism_apolarSolvationEnergyPCPLUS', &
                        rismthermo_apol%solvationEnergyPCPLUS,(/rismthermo_apol%solvationEnergyPCPLUS/),0)
                end if
                if (canCalculateUC) then
                   call thermo_print_results_line('rism_polarSolvationEnergyUC', &
                        rismthermo_pol%solvationEnergyUC, (/rismthermo_pol%solvationEnergyUC/), 0)
                   call thermo_print_results_line('rism_apolarSolvationEnergyUC', &
                        rismthermo_apol%solvationEnergyUC, (/rismthermo_apol%solvationEnergyUC/), 0)
                end if
                ! Non-free energy-based properties.
                call thermo_print_results_line('rism_polarPartialMolarVolume_dT', &
                     rismthermo_pol%partialMolarVolume_dT, (/rismthermo_pol%partialMolarVolume_dT/), 0)
                call thermo_print_results_line('rism_apolarPartialMolarVolume_dT', &
                     rismthermo_apol%partialMolarVolume_dT, (/rismthermo_apol%partialMolarVolume_dT/), 0)
                call thermo_print_results_line('rism_polarExcessParticles_dT', &
                     HUGE(1d0), rismthermo_pol%excessParticles_dT, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_results_line('rism_apolarExcessParticles_dT', &
                     HUGE(1d0), rismthermo_apol%excessParticles_dT, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_results_line('rism_polarExcessCharge_dT', &
                     sum(rismthermo_pol%excessParticles_dT*rism_3d%get_charge()) &
                     *sqrt((KB *rism3d_get_temperature(rism_3d))/COULOMB_CONST_E), &
                     rismthermo_pol%excessParticles_dT*rism_3d%get_charge()&
                     *sqrt((KB *rism3d_get_temperature(rism_3d))/COULOMB_CONST_E), &
                     rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_results_line('rism_apolarExcessCharge_dT', &
                     sum(rismthermo_apol%excessParticles_dT*rism_3d%get_charge()) &
                     *sqrt((KB *rism3d_get_temperature(rism_3d))/COULOMB_CONST_E), &
                     rismthermo_apol%excessParticles_dT*rism_3d%get_charge()&
                     *sqrt((KB *rism3d_get_temperature(rism_3d))/COULOMB_CONST_E), &
                     rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_results_line('rism_polarKirkwoodBuff_dT', HUGE(1d0), &
                     rismthermo_pol%kirkwoodBuff_dT, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_results_line('rism_apolarKirkwoodBuff_dT', HUGE(1d0), &
                     rismthermo_apol%kirkwoodBuff_dT, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_results_line('rism_polarDCFintegral_dT', HUGE(1d0), &
                     rismthermo_pol%DCFintegral_dT, rism3d_get_numAtomTypes(rism_3d))
                call thermo_print_results_line('rism_apolarDCFintegral_dT', HUGE(1d0), &
                     rismthermo_apol%DCFintegral_dT, rism3d_get_numAtomTypes(rism_3d))
             end if
          end if
       end if
    end if
    call flush(outunit)
    call rism_timer_stop(timer_write)
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
    character(len=*), intent(in) :: category, units, total, prefix, item(nitem)
    integer, intent(in) :: nitem
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
    _REAL_, intent(in) :: total, item(nitem)
    integer, intent(in) :: nitem
    ! Format for category (calculation type, e.g. excess chemical
    ! potential).
    character(len=64) :: catFmt = "(a40)"
    ! Format for a string the same width as valfmt.
    character(len=32) :: strFmt = "(a26)"
    ! Format for floating point values.
    character(len=32) :: valFmt = '(1p, 2x, e24.16e3)'
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

  
  !> Prints out the heirarchical timer summary.
  subroutine rism_printTimer()
    use amber_rism_interface
    implicit none
    call rism_timer_start(timer)
    call rism_timer_summary(timer, '|', outunit, mpicomm)
    call rism_timer_stop(timer)
  end subroutine rism_printTimer

  
  !> Finalizes all of the 3D-RISM objects and frees memory.
  subroutine rism_finalize()
    use amber_rism_interface
#ifdef RISM_CRDINTERP
    use fce_c, only : fce_destroy
#endif /*RISM_CRDINTERP*/
    use safemem
    use rism3d_solvent_c
    use rism3d_solute_c
    implicit none
    integer :: err
    integer(kind=int64) :: memstats(10)
    call rism_timer_destroy(timer_write)
    call rism_timer_destroy(timer_init)
    call rism_timer_destroy(timer)
    if (rismprm%rism == 1) then
#ifndef CUDA
       call rism3d_destroy(rism_3d)
#else
       print*, "In .F90: deleting C++ rism_3d object"
       call rism_3d%delete
#endif /*CUDA*/
#ifdef RISM_CRDINTERP
       call fce_destroy(fce_o)
#endif /*RISM_CRDINTERP*/
       call rismthermo_destroy(rismthermo)
       call rismthermo_destroy(rismthermo_pol)
       call rismthermo_destroy(rismthermo_apol)
       call rism3d_solvent_destroy(solvent)
       call rism3d_solute_destroy(solute)
       if (safemem_dealloc(ff)/= 0) &
            call rism_report_error("Deallocation in Amber-RISM interface failed")
#if defined(RISM_CRDINTERP)
       if (safemem_dealloc(atomPositions_fce)/= 0) &
            call rism_report_error("Deallocation in Amber-RISM interface failed")
#endif /*RISM_CRDINTERP*/

       if (safemem_dealloc(closurelist) /= 0) &
            call rism_report_error("Deallocation in Amber-RISM interface failed")

       if (safemem_dealloc(tolerancelist) /= 0) &
            call rism_report_error("Deallocation in Amber-RISM interface failed")
    end if
    call rism_max_memory()
  end subroutine rism_finalize

  
  !> Prints the maximum amount of memory allocated at any one time so
  !! far in the run.
  subroutine rism_max_memory()
    use amber_rism_interface
    use safemem
    implicit none
#ifdef MPI
    include "mpif.h"
#endif /*MPI*/
    integer(kind=int64) :: memstats(10), tmemstats(10)
    integer :: err, irank
    outunit = rism_report_getMUnit()
    memstats = memStatus()
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
    if (mpirank==0) then
       call MPI_REDUCE(MPI_IN_PLACE, memstats, ubound(memstats, 1), MPI_INTEGER8, &
            MPI_SUM, 0, mpicomm, err)
    else
       call MPI_REDUCE(memstats, memstats, ubound(memstats, 1), MPI_INTEGER8, &
            MPI_SUM, 0, mpicomm, err)
    end if
#  else /*USE_MPI_IN_PLACE*/
    call MPI_REDUCE(memstats, tmemstats, ubound(memstats, 1), MPI_INTEGER8, &
         MPI_SUM, 0, mpicomm, err)
    memstats = tmemstats
#  endif /*USE_MPI_IN_PLACE*/
    if (err/= 0) call rism_report_warn("RISM_MAX_MEMORY: MPI_REDUCE failed.")
#endif
    if (mpirank==0) then
       write(outunit, '(a)')
       write(outunit, '(a)') "|3D-RISM memory allocation summary"
#ifdef SANDER
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
#else
       write(outunit, '(a)') "|Type          Maximum"
       write(outunit, '(a, f12.5, a)') "|Integer  ", &
            dble(memstats(6))/BYTES_PER_GB, " GB"
       write(outunit, '(a, f12.5, a)') "|Real     ", &
            dble(memstats(7))/BYTES_PER_GB, " GB"
       write(outunit, '(a, f12.5, a)') "|Logical  ", &
            dble(memstats(8))/BYTES_PER_GB, " GB"
       write(outunit, '(a, f12.5, a)') "|Character", &
            dble(memstats(9))/BYTES_PER_GB, " GB"
       write(outunit, '(a)') "|------------------------"
       write(outunit, '(a, f12.5, a)') "|Total    ", &
            dble(memstats(10))/BYTES_PER_GB, " GB"
#endif /* SANDER */
    end if
  end subroutine rism_max_memory

  
!!! I/O: performs RISM related I/O for files that only deal with RISM data

  
  !> Provides access to writeVolumetricData for non-Fortran code.
  !! @param[in] step Step number used as a suffix.
  subroutine rism_writeVolumetricDataC(step)
    use amber_rism_interface
    implicit none
    integer, intent(in) :: step
    call rism_writeVolumetricData(rism_3d, step)
  end subroutine rism_writeVolumetricDataC

  !> Outputs volumetric data, such as solvent and electric potential
  !! distributions, to their respective files. Each distribution
  !! is written in a separate file with the step number before the
  !! suffix.  File names are taken from guvfile, huvfile, cuvfile, and
  !! similar variables.
  !! @param[in,out] this 3D-RISM object.
  !! @param[in] step Step number used as a suffix.
  subroutine rism_writeVolumetricData(this, step)
    use constants, only : COULOMB_CONST_E, KB, PI
    use amber_rism_interface
    use rism_io
    use rism3d_mrc
    use rism3d_opendx
    use rism3d_xyzv
    use safemem
    use rism_util, only : checksum
    implicit none
    
    
#if defined(MPI)
    include 'mpif.h'
#endif
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: step

    integer :: iv, ivv, nsolv, igx, igy, igz, ios, i, j, k, ig
    character(len=16) :: cstep
    character(len=64) :: suffix
    character(len=6) :: extension
  
#ifdef MPI
    integer :: err
#endif /*MPI*/

#ifndef CUDA
    procedure (writeVolumeInterface), pointer :: writeVolume => NULL()
#else
    procedure (writeVolumeInterface_cpp), pointer :: writeVolume => NULL()
#endif /*CUDA*/

    character(len = 4) :: atomName(rism3d_get_numAtomTypes(rism_3d))
    call rism3d_get_atomName(rism_3d, atomName)

    if(rism3d_get_failure(rism_3d)) then
       call rism_report_warn('SKIPPING VOLUMETRIC OUTPUT DUE TO CONVERGENCE FAILURE.')
       return
    end if
    
    call rism_timer_start(timer_write)

#ifdef CUDA
     if(volfmt .ne. 'dx') then
          call rism_report_warn("Only dx format is supported for the GPU version at the moment." &
                              //" No volume data will be written.")
          call rism_timer_stop(timer_write)
          return
     end if
#endif /*CUDA*/

    if (volfmt .eq. 'mrc') then
       extension = '.mrc'
#ifndef CUDA
       writeVolume => rism3d_mrc_map_write
#else
       writeVolume => rism3d_mrc_map_write_cpp
#endif /*CUDA*/

    else if (volfmt .eq. 'ccp4') then
       extension = '.ccp4'
       ! MRC is compatible with CCP4. However, readers may ignore some
       ! of the MRC data, such as the origin offset.
#ifndef CUDA
       writeVolume => rism3d_mrc_map_write
#else
       writeVolume => rism3d_mrc_map_write_cpp
#endif /*CUDA*/

    else if (volfmt .eq. 'dx') then
       extension = '.dx'
#ifndef CUDA
       writeVolume => rism3d_opendx_write
#else
       writeVolume => rism3d_opendx_write_cpp
#endif /*CUDA*/

    else if (volfmt .eq. 'xyzv') then
       extension = '.xyzv'
#ifndef CUDA
       writeVolume => rism3d_xyzv_write
#else
       writeVolume => rism3d_xyzv_write_cpp
#endif /*CUDA*/
    else
       call rism_report_warn("Volume format "//trim(volfmt) &
            //" is not one of the currently handled volumetric data formats." &
            //" No volume data will be written.")
       return
    end if

#ifndef CUDA
    call rism_writePdfTcfDcf(this, writeVolume, step, extension)
    call rism_writePotentialAsymptotics(this, writeVolume, step, extension)
    call rism_writeThermodynamicsDistributions(this, writeVolume, step, extension)
#ifndef MPI
    call rism_writeSmearElectronMap(this, writeVolume, step, extension)
#endif /*MPI*/
#else 
    call rism_writePdfTcfDcf(this, writeVolume, step, extension)
#endif /*CUDA*/

   call rism_timer_stop(timer_write)
  end subroutine rism_writeVolumetricData

  !> Outputs PDF, TCF, and DCF. Each distribution
  !! is written in a separate file with the step number before the
  !! suffix. 
  !! @param[in,out] this 3D-RISM object.
  !! @param[in] step Step number used as a suffix.
  subroutine rism_writePdfTcfDcf(this, writeVolume, step, extension)
    use constants, only : COULOMB_CONST_E, KB, PI
    use amber_rism_interface
    use rism_io
    use rism3d_mrc
    use rism3d_opendx
    use rism3d_xyzv
#ifndef CUDA
    use rism3d_c, only: rism3d_map_site_to_site
#else
    use rism3d_rism3d_c_mod
#endif /*CUDA*/
    use safemem
    implicit none

    type(rism3d), intent(inout) :: this
    character(len=*), intent(in) :: extension
    integer, intent(in) :: step
#ifndef CUDA
    procedure (writeVolumeInterface) :: writeVolume
#else
    procedure (writeVolumeInterface_cpp) :: writeVolume
#endif /*CUDA*/
    
    character(len=64) :: suffix
    character(len=16) :: cstep

    _REAL_, pointer :: work(:, :, :) => NULL()
    integer :: i, j, k, ig, iv

    character(len = 4) :: atomName(rism3d_get_numAtomTypes(rism_3d))
    call rism3d_get_atomName(rism_3d, atomName)
    
    if (len_trim(guvfile) /= 0 .or. len_trim(huvfile) /= 0) then
       work => safemem_realloc(work, this%get_globalDimsR(1), this%get_globalDimsR(2), this%get_localDimsR(3))
    end if
    do iv = 1, rism3d_get_numAtomTypes(this)
       write(cstep, '(i16)') step
       cstep = adjustl(cstep)
       suffix = '.'//trim(atomName(iv))//'.'//trim(cstep)
       suffix = trim(suffix)//extension
       
       ! Solute-solvent RDF.
       if (len_trim(guvfile) /= 0)  then

#ifndef CUDA
#  if defined(MPI)
          do k = 1, this%get_localDimsR(3)
             do j = 1, this%get_globalDimsR(2)
                do i = 1, this%get_globalDimsR(1)
                   work(i, j, k) = &
                        this%get_guv(i + (j - 1) * (this%get_globalDimsR(1) + 2) &
                        + (k - 1) * (this%get_globalDimsR(1) + 2) &
                        * this%get_globalDimsR(2), iv)
                end do
             end do
          end do
          call writeVolume(trim(guvfile)//suffix, work, &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
               mpirank, mpisize, mpicomm)
          if (rismprm%entropicDecomp == 1 .and. rism3d_canCalc_DT(rism_3d)) then
             do k = 1, this%get_localDimsR(3)
                do j = 1, this%get_globalDimsR(2)
                   do i = 1, this%get_globalDimsR(1)
                      work(i, j, k) = this%get_guv_dT(i + (j - 1) * (this%get_globalDimsR(1) + 2) &
                           & + (k - 1) * (this%get_globalDimsR(1) + 2) * this%get_globalDimsR(2), iv)
                   end do
                end do
             end do
             call writeVolume(trim(guvfile)//"_dT"//suffix, work, &
                  rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
                  rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
                  rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
                  rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
                  mpirank, mpisize, mpicomm)
          end if
#  else
          call writeVolume(trim(guvfile)//suffix, this%get_guv(iv), &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this))
          if (rismprm%entropicDecomp == 1 .and. rism3d_canCalc_DT(rism_3d)) then
             call writeVolume(trim(guvfile)//"_dT"//suffix, this%get_guv_dT(iv), &
                  rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
                  rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
                  rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
                  rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this))
          end if

#  endif /*defined(MPI)*/
#else /*CUDA*/
          call writeVolume(this, iv - 1, trim(guvfile)//suffix)
#endif /*CUDA*/
       endif

       ! Solute-solvent TCF.
       if (len_trim(huvfile) /= 0)  then
#ifndef CUDA
#  if defined(MPI)
          do k = 1, this%get_localDimsR(3)
             do j = 1, this%get_globalDimsR(2)
                do i = 1, this%get_globalDimsR(1)
                   work(i, j, k) = &
                        this%get_huv(i + (j - 1) * (this%get_globalDimsR(1) + 2) &
                        + (k - 1) * (this%get_globalDimsR(1) + 2) &
                        * this%get_globalDimsR(2), iv)
                end do
             end do
          end do
          call writeVolume(trim(huvfile)//suffix, work, &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
               mpirank, mpisize, mpicomm)
#  else
          call writeVolume(trim(huvfile)//suffix, this%get_huv(iv), &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this))
#  endif /*defined(MPI)*/
#else /*CUDA*/
          call writeVolume(this, iv - 1, trim(huvfile)//suffix)
#endif /*CUDA*/
       endif

       ! Solute-solvent DCF.
       if (len_trim(cuvfile) /= 0)  then
#ifndef CUDA
          do k = 1, this%get_localDimsR(3)
             do j = 1, this%get_globalDimsR(2)
                do i = 1, this%get_globalDimsR(1)
                   ig = i + (j - 1) * this%get_localDimsR(1) + &
                        (k - 1) * this%get_localDimsR(2) * this%get_localDimsR(1)
                   work(i, j, k) = &
                        this%get_cuv(i,j,k, iv)&
                        -this%get_charge(iv)*this%get_dcfLongRangeAsympR(ig)
                end do
             end do
          end do
          call writeVolume(trim(cuvfile)//suffix, work, &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
               mpirank, mpisize, mpicomm)
          ! call writeVolume(trim(cuvfile)//suffix, this%cuv(:, :, :, iv), &
          !      this%grid, this%solute, mpirank, mpisize, mpicomm)

          if (rismprm%entropicDecomp == 1 .and. rism3d_canCalc_DT(rism_3d)) then
             call writeVolume(trim(cuvfile)//"_dT"//suffix, this%get_cuv(iv), &
                  rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
                  rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
                  rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
                  rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
                  mpirank, mpisize, mpicomm)
          end if
#else /*CUDA*/
          call writeVolume(this, iv - 1, trim(cuvfile)//suffix)
#endif /*CUDA*/
       endif
    end do
    if (safemem_dealloc(work)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate WORK")
  end subroutine rism_writePdfTcfDcf

#ifndef MPI
  !> Outputs smeared solvent electron map.. Each distribution
  !! is written in a separate file with the step number before the
  !! suffix. 
  !! @param[in,out] this 3D-RISM object.
  !! @param[in] step Step number used as a suffix.
  subroutine rism_writeSmearElectronMap(this, writeVolume, step, extension)
    use constants, only : COULOMB_CONST_E, KB, PI
    use amber_rism_interface
    use rism_io
    use rism3d_mrc
    use rism3d_opendx
    use rism3d_xyzv
#ifndef CUDA
    use rism3d_c, only: createElectronDensityMap
#else
    use rism3d_rism3d_c_mod
#endif /*CUDA*/
    use rism_io, only : readRDF1D
    use safemem
    implicit none

    type(rism3d), intent(inout) :: this
    character(len=*), intent(in) :: extension
    integer, intent(in) :: step
    procedure (writeVolumeInterface) :: writeVolume 
    
    character(len=64) :: suffix
    character(len=16) :: cstep

    integer :: i, j, k, iv, elec_tot
    _REAL_, pointer :: electronRDF(:) => NULL()
    _REAL_ :: electronRDFGridSpacing
    character(len=1024) :: amberhome, rdfFilename
    
    character(len = 4) :: atomName(rism3d_get_numAtomTypes(rism_3d))
    _REAL_, pointer :: electronMap(:,:,:)

    electronMap => rism3d_get_electronMap(this)
    call rism3d_get_atomName(rism_3d, atomName)

    do iv = 1, rism3d_get_numAtomTypes(this)
       write(cstep, '(i16)') step
       cstep = adjustl(cstep)
       suffix = '.'//trim(atomName(iv))//'.'//trim(cstep)
       suffix = trim(suffix)//extension
       
       ! Smeared solvent electron map.
       if (len_trim(electronMapFile) /= 0) then

          ! dac strategy: do this for all values of iv, skipping H,
          !   assumed to be in water for now.  This will potentially
          !   create more files than we need, but they could be summed
          !   in a separate script later to create a full map.

          if( atomName(iv)(1:1) .eq. 'H' ) cycle

          call getenv('AMBERHOME', amberhome)
          rdfFilename = trim(amberhome) // '/dat/rism3d/electron_rdf/' & 
               // trim(atomName(iv)) // '.rdf'
          call readRDF1D(trim(rdfFilename), elec_tot,  &
               electronRDF, electronRDFGridSpacing)

#ifndef CUDA
          call createElectronDensityMap(this, iv, electronRDF, &
               electronRDFGridSpacing, elec_tot, &
               this%get_density(iv),  electronMap)
#else
          call rism3d_createElectronDensityMap(this, iv, electronRDF, &
               electronRDFGridSpacing, elec_tot, &
               this%get_density(iv),  electronMap)
#endif /*CUDA*/
          call writeVolume(trim(electronMapFile)//suffix, &
               electronMap, &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
               mpirank, mpisize, mpicomm)
          if (safemem_dealloc(electronRDF) /= 0) then
             call rism_report_error( &
                "rism_writeVolumetricData: Deallocate of electron RDF failed.")
          end if
       end if

    end do


    if (safemem_dealloc(electronRDF)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate ELECTRON_RDF")
  end subroutine rism_writeSmearElectronMap
#endif    


  !> Outputs PDF, TCF, and DCF. Each distribution
  !! is written in a separate file with the step number before the
  !! suffix. 
  !! @param[in,out] this 3D-RISM object.
  !! @param[in] step Step number used as a suffix.
  subroutine rism_writePotentialAsymptotics(this, writeVolume, step, extension)
    use constants, only : COULOMB_CONST_E, KB, PI
    use amber_rism_interface
    use rism_io
    use rism3d_mrc
    use rism3d_opendx
    use rism3d_xyzv
#ifndef CUDA
    use rism3d_c, only: rism3d_map_site_to_site
#else
    use rism3d_rism3d_c_mod
#endif /*CUDA*/
    use safemem
    implicit none

    type(rism3d), intent(inout) :: this
    integer, intent(in) :: step
    character(len=*), intent(in) :: extension
    procedure (writeVolumeInterface) :: writeVolume 

    character(len=64) :: suffix
    character(len=16) :: cstep
    integer :: i, j, k, iv

    logical :: warned_h, warned_all

    character(len = 4) :: atomName(rism3d_get_numAtomTypes(rism_3d))
    call rism3d_get_atomName(rism_3d, atomName)

    warned_h = .false.
    warned_all = .false.
    
    do iv = 1, rism3d_get_numAtomTypes(this)
       ! Solute-solvent electric potential.
       if (len_trim(uuvfile) /= 0)  then
          write(cstep, '(i16)') step
          cstep = adjustl(cstep)
          suffix = '.'//trim(atomName(iv))//'.'//trim(cstep)
          suffix = trim(suffix)//extension

          call writeVolume(trim(uuvfile)//suffix, rism3d_get_uuv(this,iv), &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
               mpirank, mpisize, mpicomm)
       endif
       ! Asymptotics files.  There are four but it is a simple loop.
       if (len_trim(asympfile) /= 0 .and. rism3d_get_charged(this)) then
          write(cstep, '(i16)') step
          cstep = adjustl(cstep)
          suffix = '.'//trim(cstep)
          suffix = trim(suffix)//extension
          ! h(r)
          if (rism3d_get_ionic(this)) then
             call writeVolume(trim(asympfile)//"hr"//suffix, &
                  & rism3d_get_tcfLongRangeAsympR(this), &
                  rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
                  rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
                  rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
                  rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
                  mpirank, mpisize, mpicomm)
          else if(.not. warned_h) then
             call rism_report_warn("Not writing long-range asymptotics of h(r); not used for non-ionic solvent.")
             warned_h = .true.
          end if
          ! c(r)
          call writeVolume(trim(asympfile)//"cr"//suffix, &
               & this%get_dcfLongRangeAsympR(), &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
               mpirank, mpisize, mpicomm)

       else if (.not.rism3d_get_charged(this) .and. .not. warned_all) then
          call rism_report_warn("Not writing long-range asymptotics; not used for uncharged solute.")
          warned_all = .true.
       endif
    end do
  end subroutine rism_writePotentialAsymptotics

  !> Outputs thermodynamics distributions. Each distribution
  !! is written in a separate file with the step number before the
  !! suffix. 
  !! @param[in,out] this 3D-RISM object.
  !! @param[in] step Step number used as a suffix.
  subroutine rism_writeThermodynamicsDistributions(this, writeVolume, step, extension)
    use constants, only : COULOMB_CONST_E, KB, PI
    use amber_rism_interface
    use rism_io
    use rism3d_mrc
    use rism3d_opendx
    use rism3d_xyzv
#ifndef CUDA
    use rism3d_c, only: rism3d_map_site_to_site
#else
    use rism3d_rism3d_c_mod
#endif /*CUDA*/
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
         len_trim(solventPotentialEnergyfile)/= 0 .or.  &
         len_trim(excessChemicalPotentialGFfile) /= 0 .or. &
         len_trim(solvationEnergyGFfile)/= 0 .or. &
         len_trim(entropyGFfile) /= 0 .or. &
         len_trim(excessChemicalPotentialPCPLUSfile) /= 0 .or. &
         len_trim(solvationEnergyPCPLUSfile)/= 0 .or. &
         len_trim(entropyPCPLUSfile) /= 0 .or. &
         len_trim(excessChemicalPotentialUCfile) /= 0 .or. &
         len_trim(solvationEnergyUCfile)/= 0 .or. &
         len_trim(entropyUCfile) /= 0)  then
       write(cstep, '(i16)') step
       cstep = adjustl(cstep)
       suffix = '.'//trim(cstep)
       suffix = trim(suffix)//extension
       work => safemem_realloc(work, this%get_globalDimsR(1), this%get_globalDimsR(2), this%get_localDimsR(3))
       work = 0
       ! Sum the contributions from each solvent type at each grid
       ! point and convert units to [e/A^3].
       do iv = 1, rism3d_get_numAtomTypes(this)
#  if defined(MPI)
          do k = 1, this%get_localDimsR(3)
             do j = 1, this%get_globalDimsR(2)
                do i = 1, this%get_globalDimsR(1)
                   work(i, j, k) = work(i, j, k) &
                        + this%get_guv(i + (j - 1) * (this%get_globalDimsR(1) + 2) &
                        & + (k - 1) * (this%get_globalDimsR(1) + 2) * this%get_globalDimsR(2), iv) &
                        * sqrt((KB *rism3d_get_temperature(this)) / COULOMB_CONST_E) &
                        * this%get_charge(iv) * this%get_density(iv)
                end do
             end do
          end do
#  else
          call DAXPY(rism3d_get_totalLocalPointsR(this), sqrt((KB *rism3d_get_temperature(this))/COULOMB_CONST_E) &
               *this%get_charge(iv)*this%get_density(iv), this%get_guv(1, iv), 1, work, 1)
#  endif /*defined(MPI)*/
       end do
       if (len_trim(quvfile) /= 0) then
          call writeVolume(trim(quvfile)//suffix, work, &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
               mpirank, mpisize, mpicomm)
       end if
       if (len_trim(chgDistFile) /= 0)  then
          call DSCAL(rism3d_get_totalLocalPointsR(this), rism3d_get_voxelVolume(this), work, 1)
          call writeVolume(trim(chgDistfile)//suffix, work, &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
               mpirank, mpisize, mpicomm)
       end if

       
       ! Free energy, energy and entropy.
       if (len_trim(excessChemicalPotentialfile) /= 0 .or. len_trim(entropyfile) /= 0) then
          ! Get the solvation free energy map and convert to kcal/mol.
          call rism_timer_stop(timer_write)
          if (safemem_dealloc(excessChemicalPotential_map)/= 0) &
               call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate EXCESSCHEMICALPOTENTIAL_MAP")
          excessChemicalPotential_map => rism3d_excessChemicalPotential_tot_map(this)
          excessChemicalPotential_map = excessChemicalPotential_map * KB * rism3d_get_temperature(rism_3d)
          if (safemem_dealloc(excessChemicalPotential_V_map)/= 0) &
               call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate EXCESSCHEMICALPOTENTIAL_V_MAP")
          excessChemicalPotential_V_map => rism3d_excessChemicalPotential_site_map(this)
          excessChemicalPotential_V_map = excessChemicalPotential_V_map * KB * rism3d_get_temperature(rism_3d)
          call rism_timer_start(timer_write)
       end if
       if (len_trim(solvationEnergyfile) /= 0 .or. len_trim(entropyfile) /= 0) then
          ! Get the solvation energy map and convert to kcal/mol.
          call rism_timer_stop(timer_write)
          !        solvationEnergy_map => safemem_realloc(solvationEnergy_map, &
          !             this%get_globalDimsR(1), this%get_globalDimsR(2), this%get_localDimsR(3))
          if (safemem_dealloc(solventPotentialEnergy_map)/= 0) &
               call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVENTPOTENTIALENERGY_MAP")
          solvationEnergy_map => rism3d_solvationEnergy_tot_map(this)
          solvationEnergy_map = solvationEnergy_map * KB * rism3d_get_temperature(rism_3d)
          if (safemem_dealloc(solvationEnergy_V_map)/= 0) &
               call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVATIONENERGY_V_MAP")
          solvationEnergy_V_map => rism3d_solvationEnergy_site_map(this)
          solvationEnergy_V_map = solvationEnergy_V_map * KB * rism3d_get_temperature(rism_3d)
          call rism_timer_start(timer_write)
       end if
       if (len_trim(solventPotentialEnergyfile) /= 0) then
          ! Get the solvation energy map and convert to kcal/mol.
          call rism_timer_stop(timer_write)
          if (safemem_dealloc(solventPotentialEnergy_map)/= 0) &
               call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVENTPOTENTIALENERGY_MAP")
          solventPotentialEnergy_map => rism3d_solventPotEne_tot_map(this)
          solventPotentialEnergy_map = solventPotentialEnergy_map * KB * rism3d_get_temperature(rism_3d)
          if (safemem_dealloc(solventPotentialEnergy_V_map)/= 0) &
               call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVENTPOTENTIALENERGY_V_MAP")
          solventPotentialEnergy_V_map => rism3d_solventPotEne_site_map(this)
          solventPotentialEnergy_V_map = solventPotentialEnergy_V_map * KB * rism3d_get_temperature(rism_3d)
          call rism_timer_start(timer_write)
       end if
       ! Outputting excess chemical potential map.
       if (len_trim(excessChemicalPotentialfile) /= 0) then
          call writeThermo(writeVolume,excessChemicalPotentialfile, suffix, excessChemicalPotential_map, &
            excessChemicalPotential_V_map)
       end if
       ! Outputting solvation energy map.
       if (len_trim(solvationEnergyfile) /= 0 .and. rismprm%entropicDecomp == 1 .and. &
            rism3d_canCalc_DT(rism_3d)) then
          call writeThermo(writeVolume,solvationEnergyfile, suffix, solvationEnergy_map,&
               solvationEnergy_V_map)
       end if
       ! Outputting solvent-solute energy map.
       if (len_trim(solventPotentialEnergyfile) /= 0) then
          call writeThermo(writeVolume,solventPotentialEnergyfile, suffix, solventPotentialEnergy_map, &
               solventPotentialEnergy_V_map)
       end if
       ! Outputting entropy map.
       if (len_trim(entropyfile) /= 0 .and. rismprm%entropicDecomp == 1 .and. &
            rism3d_canCalc_DT(rism_3d)) then
          call DAXPY(rism3d_get_totalLocalPointsR(this), -1d0, solvationEnergy_map, 1, excessChemicalPotential_map, 1)
          call DAXPY(rism3d_get_totalLocalPointsR(this)*rism3d_get_numAtomTypes(this), -1d0, &
               solvationEnergy_V_map, 1, excessChemicalPotential_V_map, 1)
          call writeThermo(writeVolume,entropyfile, suffix, excessChemicalPotential_map,&
               excessChemicalPotential_V_map)
       end if


       ! PCPLUS free energy, energy and entropy.
       if (len_trim(excessChemicalPotentialPCPLUSfile) /= 0 .or. len_trim(entropyPCPLUSfile) /= 0) then
          ! Get the solvation free energy map and convert to kcal/mol.
          call rism_timer_stop(timer_write)
          if (safemem_dealloc(excessChemicalPotential_map)/= 0) &
               call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate EXCESSCHEMICALPOTENTIAL_MAP")
          excessChemicalPotential_map => rism3d_excessChemicalPotentialPCPLUS_tot_map(this)
          excessChemicalPotential_map = excessChemicalPotential_map*KB*rism3d_get_temperature(rism_3d)
          call rism_timer_start(timer_write)
       end if
       ! Outputting excess chemical potential map.
       if (len_trim(excessChemicalPotentialPCPLUSfile) /= 0) then
          call writeVolume(trim(excessChemicalPotentialPCPLUSfile)//suffix, excessChemicalPotential_map, &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
               mpirank, mpisize, mpicomm)
       end if
       if (rismprm%entropicDecomp == 1 .and.  rism3d_canCalc_DT(rism_3d) &
            .and. rismthermo%partialMolarVolume_dT /= huge(1d0)) then
          if (len_trim(solvationEnergyPCPLUSfile) /= 0 .or. len_trim(entropyPCPLUSfile) /= 0) then
             ! Get the solvation energy map and convert to kcal/mol.
             call rism_timer_stop(timer_write)
             if (safemem_dealloc(solvationEnergy_map)/= 0) &
                  call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVATIONENERGY_MAP")
             solvationEnergy_map => rism3d_solvationEnergyPCPLUS_tot_map(this)
             solvationEnergy_map = solvationEnergy_map*KB*rism3d_get_temperature(rism_3d)
             call rism_timer_start(timer_write)
          end if
          ! Outputting solvation energy map.
          if (len_trim(solvationEnergyPCPLUSfile) /= 0) then
             call writeVolume(trim(solvationEnergyPCPLUSfile)//suffix, solvationEnergy_map, &
                  rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
                  rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
                  rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
                  rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
                  mpirank, mpisize, mpicomm)
          end if
          ! Outputting entropy map.
          if (len_trim(entropyPCPLUSfile) /= 0) then
             call DAXPY(rism3d_get_totalLocalPointsR(this), -1d0, solvationEnergy_map, 1, excessChemicalPotential_map, 1)
             call DAXPY(rism3d_get_totalLocalPointsR(this)*rism3d_get_numAtomTypes(this), -1d0, &
                  solvationEnergy_V_map, 1, excessChemicalPotential_V_map, 1)
             call writeVolume(trim(entropyPCPLUSfile)//suffix, excessChemicalPotential_map, &
                  rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
                  rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
                  rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
                  rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
                  mpirank, mpisize, mpicomm)
          end if
       end if

       ! Universal Correction free energy, energy and entropy.
       if ((len_trim(excessChemicalPotentialUCfile) /= 0 &
            .or. len_trim(entropyUCfile) /= 0 &
            .or. len_trim(solvationEnergyUCfile) /= 0) &
            .and. .not. canCalculateUC) then
          call rism_report_warn("UC coefficients not set. Cannot output UC grid files.")
       else
          if (len_trim(excessChemicalPotentialUCfile) /= 0 .or. len_trim(entropyUCfile) /= 0) then
             ! Get the solvation free energy map and convert to kcal/mol.
             call rism_timer_stop(timer_write)
             if (safemem_dealloc(excessChemicalPotential_map)/= 0) &
                  call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate EXCESSCHEMICALPOTENTIAL_MAP")
             excessChemicalPotential_map => rism3d_excessChemicalPotentialUC_tot_map(this, rismprm%uccoeff)
             excessChemicalPotential_map = excessChemicalPotential_map*KB*rism3d_get_temperature(rism_3d)
             call rism_timer_start(timer_write)
          end if
          ! Outputting excess chemical potential map.
          if (len_trim(excessChemicalPotentialUCfile) /= 0) then
             call writeVolume(trim(excessChemicalPotentialUCfile)//suffix, excessChemicalPotential_map, &
                  rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
                  rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
                  rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
                  rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
                  mpirank, mpisize, mpicomm)
          end if
          if (rismprm%entropicDecomp == 1 .and.  rism3d_canCalc_DT(rism_3d) &
               .and. rismthermo%partialMolarVolume_dT /= huge(1d0)) then
             if (len_trim(solvationEnergyUCfile) /= 0 .or. len_trim(entropyUCfile) /= 0) then
                ! Get the solvation energy map and convert to kcal/mol.
                call rism_timer_stop(timer_write)
                if (safemem_dealloc(solvationEnergy_map)/= 0) &
                     call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVATIONENERGY_MAP")
                solvationEnergy_map => rism3d_solvationEnergyUC_tot_map(this, rismprm%uccoeff)
                solvationEnergy_map = solvationEnergy_map*KB*rism3d_get_temperature(rism_3d)
                call rism_timer_start(timer_write)
             end if
             ! Outputting solvation energy map.
             if (len_trim(solvationEnergyUCfile) /= 0) then
                call writeVolume(trim(solvationEnergyUCfile)//suffix, solvationEnergy_map, &
                     rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
                     rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
                     rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
                     rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
                     mpirank, mpisize, mpicomm)
             end if
             ! Outputting entropy map.
             if (len_trim(entropyUCfile) /= 0) then
                call DAXPY(rism3d_get_totalLocalPointsR(this), -1d0, solvationEnergy_map, 1, excessChemicalPotential_map, 1)
                call DAXPY(rism3d_get_totalLocalPointsR(this)*rism3d_get_numAtomTypes(this), -1d0, &
                     solvationEnergy_V_map, 1, excessChemicalPotential_V_map, 1)
                call writeVolume(trim(entropyUCfile)//suffix, excessChemicalPotential_map, &
                     rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
                     rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
                     rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
                     rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
                     mpirank, mpisize, mpicomm)
             end if
          end if
       end if

       ! Gaussian fluctuation correction free energy, energy and entropy
       if (len_trim(excessChemicalPotentialGFfile) /= 0 .or. len_trim(entropyGFfile) /= 0) then
          ! Get the solvation free energy map and convert to kcal/mol.
          call rism_timer_stop(timer_write)
          if (safemem_dealloc(solventPotentialEnergy_map)/= 0) &
               call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVENTPOTENTIALENERGY_MAP")
          excessChemicalPotential_map => rism3d_excessChemicalPotentialGF_tot_map(this)
          excessChemicalPotential_map = excessChemicalPotential_map*KB*rism3d_get_temperature(rism_3d)
          if (safemem_dealloc(excessChemicalPotential_V_map)/= 0) &
               call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate EXCESSCHEMICALPOTENTIAL_V_MAP")
          excessChemicalPotential_V_map => rism3d_excessChemicalPotentialGF_site_map(this)
          excessChemicalPotential_V_map = excessChemicalPotential_V_map*KB*rism3d_get_temperature(rism_3d)
          call rism_timer_start(timer_write)
       end if
       if (len_trim(solvationEnergyGFfile) /= 0 .or. len_trim(entropyGFfile) /= 0) then
          ! Get the solvation energy map and convert to kcal/mol.
          call rism_timer_stop(timer_write)
          !        solvationEnergy_map => safemem_realloc(solvationEnergy_map, &
          !             this%get_globalDimsR(1), this%get_globalDimsR(2), this%get_localDimsR(3))
          if (safemem_dealloc(solvationEnergy_map)/= 0) &
               call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVATIONENERGY_MAP")
          solvationEnergy_map => rism3d_solvationEnergyGF_tot_map(this)
          solvationEnergy_map = solvationEnergy_map*KB*rism3d_get_temperature(rism_3d)
          if (safemem_dealloc(solvationEnergy_V_map)/= 0) &
               call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVATIONENERGY_V_MAP")
          solvationEnergy_V_map => rism3d_solvationEnergyGF_site_map(this)
          solvationEnergy_V_map = solvationEnergy_V_map*KB*rism3d_get_temperature(rism_3d)
          call rism_timer_start(timer_write)
       end if
       ! Outputting excess chemical potential map.
       if (len_trim(excessChemicalPotentialGFfile) /= 0) then
          call writeThermo(writeVolume,excessChemicalPotentialGFfile, suffix, excessChemicalPotential_map, &
               excessChemicalPotential_V_map)
       end if
       ! Outputting solvation energy map.
       if (len_trim(solvationEnergyGFfile) /= 0 .and. rismprm%entropicDecomp == 1 .and. &
            rism3d_canCalc_DT(rism_3d)) then
          call writeThermo(writeVolume,solvationEnergyGFfile, suffix, solvationEnergy_map, &
               solvationEnergy_V_map)
       end if
       ! Outputting entropy map.
       if (len_trim(entropyGFfile) /= 0 .and. rismprm%entropicDecomp == 1 .and. &
            rism3d_canCalc_DT(rism_3d)) then
          call DAXPY(rism3d_get_totalLocalPointsR(this), -1d0, solvationEnergy_map, 1, excessChemicalPotential_map, 1)
          call DAXPY(rism3d_get_totalLocalPointsR(this)*rism3d_get_numAtomTypes(this), -1d0, &
               solvationEnergy_V_map, 1, excessChemicalPotential_V_map, 1)
          call writeThermo(writeVolume,entropyGFfile, suffix, excessChemicalPotential_map, &
               excessChemicalPotential_V_map)
       end if
       if (rismprm%molReconstruct) then
          call rism_writeMolReconstruct(this, writeVolume, suffix)
       end if
    endif
    if (safemem_dealloc(excessChemicalPotential_map)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate EXCESSCHEMICALPOTENTIAL_MAP")
    if (safemem_dealloc(solvationEnergy_map)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVATIONENERGY_MAP")
    if (safemem_dealloc(solventPotentialEnergy_map)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVENTPOTENTIALENERGY_MAP")
    if (safemem_dealloc(excessChemicalPotential_V_map)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate EXCESSCHEMICALPOTENTIAL_V_MAP")
    if (safemem_dealloc(solvationEnergy_V_map)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVATIONENERGY_V_MAP")
    if (safemem_dealloc(solventEntropy_V_map)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVENTENTROPY_V_MAP")
    if (safemem_dealloc(solventPotentialEnergy_V_map)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVENTPOTENTIALENERGY_V_MAP")
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

      character(len = 4) :: atomName(rism3d_get_numAtomTypes(rism_3d))
      call rism3d_get_atomName(rism_3d, atomName)

      call wv(trim(root)//trim(suffix), data, &
           rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
           rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
           rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
           rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
           mpirank, mpisize, mpicomm)
      do iv = 1, rism3d_get_numAtomTypes(this)
         call wv(trim(root)//'.'//trim(atomName(iv))//trim(suffix), &
              data_sites(:, :, :, iv), &
              rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
              rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
              rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
              rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
              mpirank, mpisize, mpicomm)
      end do
    end subroutine writeThermo
   end subroutine rism_writeThermodynamicsDistributions
  
  !> Outputs molecular reconstruction of volumetric data. Each distribution
  !! is written in a separate file with the step number before the
  !! suffix. 
  !! @param[in,out] this 3D-RISM object.
  !! @param[in] step Step number used as a suffix.
  subroutine rism_writeMolReconstruct(this, writeVolume, suffix)
    use constants, only : COULOMB_CONST_E, KB, PI
    use amber_rism_interface
    use rism_io
    use rism3d_mrc
    use rism3d_opendx
    use rism3d_xyzv
#ifndef CUDA
    use rism3d_c, only: rism3d_map_site_to_site
#else
    use rism3d_rism3d_c_mod
#endif /*CUDA*/
    use safemem
    implicit none

    type(rism3d), intent(inout) :: this
    character(len=64), intent(in) :: suffix
    procedure (writeVolumeInterface) :: writeVolume 
    _REAL_, pointer::excessChemicalPotential_V_map(:, :, :, :) => NULL(), &
         solvationEnergy_V_map(:, :, :, :) => NULL(), &
         solventPotentialEnergy_V_map(:, :, :, :) => NULL(),&
         solventEntropy_V_map(:, :, :, :) => NULL()
    integer :: center_site
    integer :: iv
    
    !
    ! excess chemical potential
    !
    center_site = 1
    if (len(trim(excessChemicalPotentialfile)) /= 0 .or. len(trim(entropyfile)) /=0) then
       call rism_timer_stop(timer_write)
       excessChemicalPotential_V_map => rism3d_excessChemicalPotential_site_map(this)
       call rism_timer_start(timer_write)
       call writeReconstruct(excessChemicalPotentialfile, suffix, &
            excessChemicalPotential_V_map)
    end if
    !
    ! Solvation energy
    !
    if (len(trim(solvationEnergyfile)) /= 0 .or. len(trim(entropyfile)) /=0 &
         .and. rismprm%entropicDecomp == 1 .and. rism3d_canCalc_DT(rism_3d)) then

       call rism_timer_stop(timer_write)
       solvationEnergy_V_map => rism3d_solvationEnergy_site_map(this)
       call rism_timer_start(timer_write)
       call writeReconstruct(solvationEnergyFile, suffix, &
            solvationEnergy_V_map)
    end if
    !
    ! solute-solvent potential energy
    !
    if (len(trim(solventPotentialEnergyfile)) /= 0 ) then

       call rism_timer_stop(timer_write)
       solventPotentialEnergy_V_map => rism3d_solventPotEne_site_map(this)
       call rism_timer_start(timer_write)
       call writeReconstruct(solventPotentialEnergyFile, suffix, &
            solventPotentialEnergy_V_map)
    end if
    !
    ! Solvation Entropy
    !
    if (len(trim(entropyfile)) /=0 .and. rismprm%entropicDecomp == 1 .and. &
         rism3d_canCalc_DT(rism_3d)) then
       solventEntropy_V_map => safemem_realloc(solventEntropy_V_map, &
            this%get_globalDimsR(1), this%get_globalDimsR(2), this%get_localDimsR(3), &
            rism3d_get_numAtomTypes(this))
       call DCOPY(rism3d_get_totalLocalPointsR(this)*rism3d_get_numAtomTypes(this),&
            excessChemicalPotential_V_map,1,solventEntropy_V_map,1)
       
       call DAXPY(rism3d_get_totalLocalPointsR(this)*rism3d_get_numAtomTypes(this), -1d0, &
            solvationEnergy_V_map, 1, solventEntropy_V_map, 1)
       call writeVolume(trim(entropyfile)//'.mol'//suffix, &
            solvententropy_V_map(:, :, :, center_site), &
            rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
            rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
            rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
            rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
            mpirank, mpisize, mpicomm)
    end if
    if (safemem_dealloc(excessChemicalPotential_V_map)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate EXCESSCHEMICALPOTENTIAL_V_MAP")
    if (safemem_dealloc(solvationEnergy_V_map)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVATIONENERGY_V_MAP")
    if (safemem_dealloc(solventEntropy_V_map)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVENTENTROPY_V_MAP")
    if (safemem_dealloc(solventPotentialEnergy_V_map)/= 0) &
         call rism_report_error("RISM_WRITESOLVEDIST: failed to deallocate SOLVENTPOTENTIALENERGY_V_MAP")
  contains
    subroutine writeReconstruct(root, suffix, data_sites)
      implicit none
      character(len=*), intent(in) :: root, suffix
      _REAL_, pointer, intent(inout) :: data_sites(:,:,:,:)
      call rism3d_map_site_to_site(rism_3d, data_sites, center_site)
       data_sites = data_sites * KB * rism3d_get_temperature(rism_3d)
       do iv = 2,2
          data_sites(:, :, :, center_site) = &
               data_sites(:, :, :, center_site)&
               + data_sites(:, :, :, iv)
       end do
       if (len(trim(root)) /= 0 ) then
          call writeVolume(trim(root)//'.mol'//trim(suffix), &
               data_sites(:, :, :, center_site), &
               rism3d_get_periodic(this), rism3d_get_boxLength(this), this%get_localDimsR(), &
               rism3d_get_unitCellVectorsK(this), rism3d_get_unitCellVectorsR(this), this%get_globalDimsR(), &
               rism3d_get_unitCellAngles(this), rism3d_get_voxelVectorsR(this), rism3d_get_spacing(this), &
               rism3d_get_totalLocalPointsR(this), rism3d_get_translation(this), rism3d_get_centerOfMass(this), &
               mpirank, mpisize, mpicomm)
       end if
    end subroutine writeReconstruct
  end subroutine rism_writeMolReconstruct

  
!!!! INTERPOLATION RESTART FILE I/O
#if 0
#if defined(RISM_CRDINTERP)
  !> Writes the interpolation restart file.
  !! @param[in] this rism3d object.
  !! @param[in] atomPositions The position of each solute atom for each step.
  !! @param[in] frc The force of each solute atom for each step.
  !! @param[in] nstep Number of steps.
  !! FIXME: modify this to use the NetCDF format
  subroutine rismRestartWrite(this, atomPositions, frc, nstep)
    use amber_rism_interface
    use rism_util, only : freeUnit
    implicit none
#include "files.h"
    type(rism3d) :: this
    _REAL_, intent(in) :: atomPositions(3, this%solute%numAtoms, nstep), frc(3, this%solute%numAtoms, nstep)
    integer, intent(in) :: nstep
    integer :: istep
    integer :: unit
    !this is performed by the master node only:
    if (mpirank == 0) then
#ifdef RISM_DEBUG
       write(outunit, *) 'entering rismrestrtwrit'
       call flush(6)
#endif /*RISM_DEBUG*/
       unit = freeUnit()
       if (len_trim(rismcrdfil) /= 0)  then
          open(unit=unit, file=rismcrdfil, status='new', form='FORMATTED', iostat=stat)
          write(unit, '(i8)') nstep
          do istep = 1, nstep
             call corpac(atomPositions(1:3, 1:this%solute%numAtoms, istep), 1, this%solute%numAtoms*3, unit, .true.)
          end do
          close(unit)
       end if
       if (len_trim(rismfrcfil) /= 0)  then
          open(unit=unit, file=rismfrcfil, status='new', form='FORMATTED', iostat=stat)
          write(unit, '(i8)') nstep
          do istep = 1, nstep
             call corpac(frc(1:3, 1:this%solute%numAtoms, istep), 1, this%solute%numAtoms*3, unit, .true.)
          end do
          close(unit)
       end if
#ifdef RISM_DEBUG
       write(outunit, *) 'done rismrestrtwrit'
       call flush(6)
#endif /*RISM_DEBUG*/
    end if
  end subroutine rismRestartWrite

  
  !! Reads in the interpolation restart file
  !!IN:
  !!   atomPositions    :: the position of each atom for each step read in
  !!   frc     :: the force of each atom for each step read in
  !!   this%solute%numAtoms    :: number of solute atoms
  !!   nstep   :: maximum number of steps to read
  !!   nsample :: will hold the total number of steps read in
  !!TODO: Modify this to use the NetCDF format.
  subroutine rismRestartRead(this, atomPositions, frc, nstep, nsample)
    use amber_rism_interface
    use rism_util, only : freeUnit
    implicit none
#include "files.h"
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d) :: this
    _REAL_, intent(out) :: atomPositions(3, this%solute%numAtoms, nstep), frc(3, this%solute%numAtoms, nstep)
    integer, intent(in) :: nstep
    integer, intent(out) :: nsample
    integer :: istep, iatu, csteps, fsteps, err
    integer :: unit
    !this is performed by the master node only:
    if (mpirank == 0) then
#ifdef RISM_DEBUG
       write(outunit, *) "Reading interpolation restart files..."
       write(outunit, *) 'entering rismrestrtread'
       call flush(6)
#endif /*RISM_DEBUG*/
       unit=freeeUnit()
       if (len_trim(rismcrdrstfil) /= 0)  then
          open(unit=unit, file=rismcrdfil, status='old', form='FORMATTED', iostat=stat)
          read(unit, '(i12)') csteps
          !discard the extra samples at the beginning
#ifdef RISM_DEBUG
          write(outunit, *) csteps, nstep, csteps-nstep
          call flush(6)
#endif /*RISM_DEBUG*/
          do istep = 1, csteps-nstep
             !            read(unit, '(10f8.3)') (atomPositions(1:3, iatu, 1), iatu=1, this%solute%numAtoms)
             read(unit, '(10f21.16)') (atomPositions(1:3, iatu, 1), iatu=1, this%solute%numAtoms)
          end do
          do istep = 1, min(csteps, nstep)
             !            read(unit, '(10f8.3)') (atomPositions(1:3, iatu, istep), iatu=1, this%solute%numAtoms)
             read(unit, '(10f21.16)') (atomPositions(1:3, iatu, istep), iatu=1, this%solute%numAtoms)
          end do
          close(unit)
       endif
       if (len_trim(rismfrcrstfil) /= 0)  then
          open(unit=unit, file=rismfrcfil, status='old', form='FORMATTED', iostat=stat)
          read(unit, '(i12)') fsteps
          !discard the extra samples at the beginning
          do istep = 1, fsteps-nstep
             !            read(unit, '(10f8.3)') (frc(1:3, iatu, 1), iatu=1, this%solute%numAtoms)
             read(unit, '(10f21.16)') (frc(1:3, iatu, 1), iatu=1, this%solute%numAtoms)
          end do
          do istep = 1, min(fsteps, nstep)
             !            read(unit, '(10f8.3)') (frc(1:3, iatu, istep), iatu=1, this%solute%numAtoms)
             read(unit, '(10f21.16)') (frc(1:3, iatu, istep), iatu=1, this%solute%numAtoms)
          end do
          close(unit)
          if (csteps /= fsteps) then
             call rism_report_error('RISM interpolation restart files have different numbers of entries.')
          end if
          nsample = min(csteps, nstep)
       end if
#ifdef RISM_DEBUG
       write(outunit, *) nsample, " samples read"
       write(outunit, *) 'done rismrestrtread'
       call flush(6)
#endif /*RISM_DEBUG*/
    end if
    !
    !broadcast results to the other processes: nsamples, atomPositions and frc
    !
!!$#if defined(MPI)
!!$      CALL MPI_BCAST(nsample, 1, MPI_INTEGER, 0, mpicomm, err)
!!$      CALL MPI_BCAST(atomPositions, 3*this%solute%numAtoms*nstep, MPI_DOUBLE_PRECISION, 0, mpicomm, err)
!!$      CALL MPI_BCAST(frc, 3*this%solute%numAtoms*nstep, MPI_DOUBLE_PRECISION, 0, mpicomm, err)
!!$#endif /*defined(MPI)*/

  end subroutine rismRestartRead
#endif /*RISM_CRDINTERP*/
#endif  /* 0 */

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
    implicit none
    include 'mpif.h'
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
       call mpi_bcast(rismprm%buffer, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast BUFFER")
       call mpi_bcast(rismprm%grdspc, 3, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast GRDSPC")
       call mpi_bcast(rismprm%solvbox, 3, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast SOLVBOX")
       call mpi_bcast(rismprm%mdiis_del, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast MDIIS_DEL")
       call mpi_bcast(rismprm%mdiis_restart, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast MDIIS_RESTART")
       call mpi_bcast(rismprm%fcecut, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCECUT")
       call mpi_bcast(rismprm%uccoeff, size(rismprm%uccoeff), mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCECUT")
       call mpi_bcast(rismprm%treeDCFMAC, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeDCFMAC")
       call mpi_bcast(rismprm%treeTCFMAC, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeTCFMAC")
       call mpi_bcast(rismprm%treeCoulombMAC, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeCoulombMAC")
       call mpi_bcast(rismprm%asympKSpaceTolerance, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast asympKSpaceTolerance")
       call mpi_bcast(rismprm%ljTolerance, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast ljTolerance")
       call mpi_bcast(rismprm%chargeSmear, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast chargeSmear")


       call mpi_bcast(rismprm%uccoeff, size(rismprm%uccoeff), mpi_double_precision, 0, mpicomm, err)

       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast uccoeff")
       call mpi_bcast(rismprm%ng3, 3, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast NG3")
       call mpi_bcast(rismprm%polarDecomp, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast POLARDECOMP")
       call mpi_bcast(rismprm%entropicDecomp, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast ENTROPICDECOMP")
       call mpi_bcast(rismprm%gfCorrection, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast GFCORRECTION")
       call mpi_bcast(rismprm%pcplusCorrection, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast pcplusCorrection")
       call mpi_bcast(rismprm%biasPotential, 1, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast BIASPOTENTIAL")
       call mpi_bcast(rismprm%asympCorr, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast ASYMPCORR")
       call mpi_bcast(rismprm%maxstep, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast MAXSTEP")
       call mpi_bcast(rismprm%npropagate, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast NPROPAGATE")
       call mpi_bcast(rismprm%centering, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast CENTERING")
       call mpi_bcast(rismprm%zerofrc, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast ZEROFRC")
       call mpi_bcast(rismprm%apply_rism_force, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast APPLY_RISM_FORCE")
       call mpi_bcast(rismprm%rismnrespa, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast RISMNRESPA")
       call mpi_bcast(rismprm%fcestride, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCESTRIDE")
       call mpi_bcast(rismprm%fcenbasis, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCENBASIS")
       call mpi_bcast(rismprm%fcenbase,1,mpi_integer,0,mpicomm,err)
       if (err /=0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCENBASE")
       call mpi_bcast(rismprm%fcecrd, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCECRD")
       call mpi_bcast(rismprm%fceweigh,1,mpi_integer,0,mpicomm,err)
       if (err /=0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCEWEIGH")
       call mpi_bcast(rismprm%fcetrans,1,mpi_integer,0,mpicomm,err)
       if (err /=0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCETRANS")
       call mpi_bcast(rismprm%fcesort,1,mpi_integer,0,mpicomm,err)
       if (err /=0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCESORT")
       call mpi_bcast(rismprm%fceifreq,1,mpi_integer,0,mpicomm,err)
       if(err /=0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCEIFREQ")
       call mpi_bcast(rismprm%fcentfrcor,1,mpi_integer,0,mpicomm,err)
       if(err /=0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCENTFRCOR")
       call mpi_bcast(rismprm%fceenormsw,1,mpi_double_precision,0,mpicomm,err)
       if (err /=0) call rism_report_error&
            ("RISM3D interface: could not broadcast FCENORMSW")
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
       call mpi_bcast(rismprm%treeDCFOrder, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeDCFOrder")
       call mpi_bcast(rismprm%treeTCFOrder, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeTCFOrder")
       call mpi_bcast(rismprm%treeCoulombOrder, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeCoulombOrder")
       call mpi_bcast(rismprm%treeDCFN0, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeDCFN0")
       call mpi_bcast(rismprm%treeTCFN0, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeTCFN0")
       call mpi_bcast(rismprm%treeCoulombN0, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeCoulombN0")

       call mpi_bcast(rismprm%selftest, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast SELFTEST")
       call mpi_bcast(rismprm%treeDCF, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeDCF")
       call mpi_bcast(rismprm%treeTCF, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeTCF")
       call mpi_bcast(rismprm%treeCoulomb, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast treeCoulomb")

       if (mpirank==0) &
            nclosure=ubound(closurelist, 1)
       call mpi_bcast(nclosure, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast PROGRESS")
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

       call mpi_bcast(rismprm%stoponfail, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast STOPONFAIL")
       
#ifdef SANDER
       call mpi_bcast(rismprm%write_thermo, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast WRITE_THERMO")
#endif

       ! These are not being used currently.
       call mpi_bcast(pa_orient, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast PA_ORIENT")
       call mpi_bcast(rmsd_orient, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast RMSD_ORIENT")

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
       call mpi_bcast(asympfile, len(asympfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast ASYMPFILE")
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
       call mpi_bcast(excessChemicalPotentialGFfile, len(excessChemicalPotentialGFfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast EXCESSCHEMICALPOTENTIALGFFILE")
       call mpi_bcast(solvationEnergyGFfile, len(solvationEnergyGFfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast SOLVATIONENERGYGFFILE")
       call mpi_bcast(entropyGFfile, len(entropyGFfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast ENTROPYGFFILE")
       call mpi_bcast(excessChemicalPotentialPCPLUSfile, len(excessChemicalPotentialPCPLUSfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast EXCESSCHEMICALPOTENTIALPCPLUSFILE")
       call mpi_bcast(solvationEnergyPCPLUSfile, len(solvationEnergyPCPLUSfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast SOLVATIONENERGYPCPLUSFILE")
       call mpi_bcast(entropyPCPLUSfile, len(entropyPCPLUSfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast ENTROPYPCPLUSFILE")
       call mpi_bcast(excessChemicalPotentialUCfile, len(excessChemicalPotentialUCfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast EXCESSCHEMICALPOTENTIALUCFILE")
       call mpi_bcast(solvationEnergyUCfile, len(solvationEnergyUCfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast SOLVATIONENERGYUCFILE")
       call mpi_bcast(entropyUCfile, len(entropyUCfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast ENTROPYUCFILE")
       call mpi_bcast(solventPotentialEnergyfile, len(solventPotentialEnergyfile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast POTUVFILE")
#ifndef MPI
       call mpi_bcast(electronMapFile, len(electronMapFile), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast ELECTRONMAPFILE")
#endif
       call mpi_bcast(rismprm%molReconstruct, 1, mpi_integer, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast MOLRECONSTRUCT")
       call mpi_bcast(volfmt, len(volfmt), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast VOLFMT")
       call mpi_bcast(periodicPotential, len(periodicPotential), mpi_character, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast PERIODIC")
       call mpi_bcast(unitCellDimensions, 6, mpi_double_precision, 0, mpicomm, err)
       if (err /= 0) call rism_report_error&
            ("RISM3D interface: could not broadcast UNITCELLDIMENSIONS")
    end if
  end subroutine rism_mpi_bcast
#endif /*MPI*/

  !> Sets default values for 3D-RISM paramters.
  subroutine defaults()
    use amber_rism_interface
    implicit none

    closurelist => safemem_realloc(closurelist, len(closurelist), nclosuredefault)
    closurelist(2:)           = ''
    closurelist(1)            = 'KH'
    rismprm%closureOrder      = 1
    rismprm%asympCorr         = .true.
    rismprm%uccoeff           = 0d0
    rismprm%biasPotential     = 0
    rismprm%polarDecomp       = 0
    rismprm%entropicDecomp    = 0
    rismprm%gfCorrection     = 0
    rismprm%pcplusCorrection = 0
    periodicPotential         = ''

    !solvation box
    rismprm%solvcut         = -1
    rismprm%buffer          = 14d0
    rismprm%grdspc          = 0.5d0
    rismprm%ng3             = -1
    rismprm%solvbox         = -1d0

    !convergence
    tolerancelist => safemem_realloc(tolerancelist, nclosuredefault)
    tolerancelist             = HUGE(1d0)
    tolerancelist(1)          = 1d-5
    rismprm%mdiis_del         = 0.7d0
    rismprm%mdiis_nvec        = 5
    rismprm%mdiis_method      = 2
    rismprm%mdiis_restart     = 10d0
    rismprm%maxstep           = 10000
    rismprm%npropagate        = 5
    rismprm%stoponfail        = .true.
    
    !imin = 1 (minimization)
    rismprm%centering        = -9999
    rismprm%zerofrc          = 1

    !imin = 5 (trajectory analysis)
    rismprm%apply_rism_force = 1
    pa_orient        = 0
    rmsd_orient      = 0

    !imin = 0 (MD)
    rismprm%rismnrespa       = 1
#ifdef RISM_CRDINTERP
    rismprm%fcestride        = 0
    rismprm%fcecut           = 9999d0
    rismprm%fcenbasis        = 20
    rismprm%fcenbase         = 20
    rismprm%fcecrd           = 0
    rismprm%fceweigh         = 0
    rismprm%fcetrans         = 0
    rismprm%fcesort          = 0
    rismprm%fceifreq         = 1
    rismprm%fcentfrcor       = 0
    rismprm%fcewrite         = 0
    rismprm%fceread          = 0
    rismprm%fceenormsw       = 0.d0
#endif

    !output
    rismprm%saveprogress     = 0
    rismprm%ntwrism          = -1
    rismprm%verbose          = 0
    rismprm%progress         = 1
    volfmt                   = 'mrc'
    rismprm%selftest         = .false.

    ! molecular reconstruction
    rismprm%molReconstruct = .false.
    
    !treecode/direct sum
    rismprm%treeDCF     = .true.
    rismprm%treeTCF     = .true.
    rismprm%treeCoulomb = .false.
    rismprm%treeDCFMAC     = 0.1d0
    rismprm%treeTCFMAC     = 0.1d0
    rismprm%treeCoulombMAC = 0.1d0
    rismprm%treeDCFOrder     = 2
    rismprm%treeTCFOrder     = 2
    rismprm%treeCoulombOrder = 2
    rismprm%treeDCFN0     = 500
    rismprm%treeTCFN0     = 500
    rismprm%treeCoulombN0 = 500
    
    !long-range asymptotics k-space tolerance
    rismprm%asympKSpaceTolerance = -1
    !Lennard-Jones tolerance
    rismprm%ljTolerance = -1
    !charge smear
    rismprm%chargeSmear = 1d0
#ifdef SANDER
    rismprm%write_thermo=1
#endif
  end subroutine defaults

  
  !> Transfers RISM settings from a rismprm_t type (user) to the
  !! rismprm parameter set for the 3D-RISM calculation.  For each
  !! possible setting the new value is used IFF the value > -9999.
  subroutine update_param(user)
    use amber_rism_interface
    implicit none
    type(rismprm_t), intent(in) :: user
    rismprm%rism = user%rism
    rismprm%molReconstruct = user%molReconstruct
    ! asympCorr, periodic(?), direct sums and selftest are logical, so we can't perform the test
    rismprm%asympCorr = user%asympCorr
    rismprm%selftest = user%selftest
    rismprm%treeDCF = user%treeDCF
    rismprm%treeTCF = user%treeTCF
    rismprm%treeCoulomb = user%treeCoulomb
    rismprm%stoponfail = user%stoponfail
    
    if (user%closureOrder > -9999) rismprm%closureOrder = user%closureOrder
    if (user%uccoeff(1) > -9999) rismprm%uccoeff = user%uccoeff
    if (user%biasPotential /= 0) rismprm%biasPotential = user%biasPotential
    if (user%polarDecomp > -9999) rismprm%polarDecomp = user%polarDecomp
    if (user%entropicDecomp > -9999) rismprm%entropicDecomp = user%entropicDecomp
    if (user%gfCorrection > -9999) rismprm%gfCorrection = user%gfCorrection
    if (user%pcplusCorrection > -9999) rismprm%pcplusCorrection = user%pcplusCorrection
    ! We need to know if this is a periodic or open boundary to set the default centering.
    if (user%centering == -9999) then
       if (periodicPotential /= '') then
          ! For periodic calculations, the default is 0 (no centering)
          rismprm%centering = 0
       else
          ! For open-boundary calculations, the default is 1 (CoM centering)
          rismprm%centering = 1
       end if
    else
       ! use the user preferences
       rismprm%centering = user%centering
    end if

    if (user%buffer > -9999) rismprm%buffer = user%buffer
    if (user%solvcut > 0) then
       rismprm%solvcut = user%solvcut
    else
       rismprm%solvcut = rismprm%buffer
    end if

    if (user%grdspc(1) > -9999) rismprm%grdspc = user%grdspc
    if (user%ng3(1) > -9999) rismprm%ng3 = user%ng3
    if (user%solvbox(1) > -9999) rismprm%solvbox = user%solvbox
    if (user%mdiis_del > -9999) rismprm%mdiis_del = user%mdiis_del
    if (user%mdiis_nvec > -9999) rismprm%mdiis_nvec = user%mdiis_nvec
    if (user%mdiis_method > -9999) rismprm%mdiis_method = user%mdiis_method
    if (user%mdiis_restart > -9999) rismprm%mdiis_restart = user%mdiis_restart
    if (user%maxstep > -9999) rismprm%maxstep = user%maxstep
    if (user%npropagate > -9999) rismprm%npropagate = user%npropagate
    if (user%zerofrc > -9999) rismprm%zerofrc = user%zerofrc
    if (user%apply_rism_force > -9999) rismprm%apply_rism_force = user%apply_rism_force
    if (user%rismnrespa > -9999) rismprm%rismnrespa = user%rismnrespa
    if (user%fcestride > -9999) rismprm%fcestride = user%fcestride
    if (user%fcecut > -9999) rismprm%fcecut = user%fcecut
    if (user%fcenbasis > -9999) rismprm%fcenbasis = user%fcenbasis
    if (user%fcenbase > -9999) rismprm%fcenbase = user%fcenbase
    if (user%fcecrd > -9999) rismprm%fcecrd = user%fcecrd
    if (user%fceweigh > -9999) rismprm%fceweigh = user%fceweigh
    if (user%fcetrans > -9999) rismprm%fcetrans = user%fcetrans
    if (user%fcesort > -9999) rismprm%fcesort = user%fcesort
    if (user%fceifreq > -9999) rismprm%fceifreq = user%fceifreq
    if (user%fcentfrcor > -9999) rismprm%fcentfrcor = user%fcentfrcor
    if (user%fcewrite > -9999) rismprm%fcewrite = user%fcewrite
    if (user%fceread > -9999) rismprm%fceread = user%fceread
    if (user%fceenormsw > -9999) rismprm%fceenormsw = user%fceenormsw
    if (user%treeDCFMAC > -9999) rismprm%treeDCFMAC = user%treeDCFMAC
    if (user%treeTCFMAC > -9999) rismprm%treeTCFMAC = user%treeTCFMAC
    if (user%treeCoulombMAC > -9999) rismprm%treeCoulombMAC = user%treeCoulombMAC
    if (user%treeDCFOrder > -9999) rismprm%treeDCFOrder = user%treeDCFOrder
    if (user%treeTCFOrder > -9999) rismprm%treeTCFOrder = user%treeTCFOrder
    if (user%treeCoulombOrder > -9999) rismprm%treeCoulombOrder = user%treeCoulombOrder
    if (user%treeDCFN0 > -9999) rismprm%treeDCFN0 = user%treeDCFN0
    if (user%treeTCFN0 > -9999) rismprm%treeTCFN0 = user%treeTCFN0
    if (user%treeCoulombN0 > -9999) rismprm%treeCoulombN0 = user%treeCoulombN0
    if (user%asympKSpaceTolerance > -9999) rismprm%asympKSpaceTolerance = user%asympKSpaceTolerance
    if (user%ljTolerance > -9999) rismprm%ljTolerance = user%ljTolerance
    if (user%chargeSmear > -9999) rismprm%chargeSmear = user%chargeSmear
    if (user%saveprogress > -9999) rismprm%saveprogress = user%saveprogress
    if (user%ntwrism > -9999) rismprm%ntwrism = user%ntwrism
    if (user%verbose > -9999) rismprm%verbose = user%verbose
    if (user%progress > -9999) rismprm%progress = user%progress
  end subroutine update_param

  
  !> Reads the RISM namelist from the input file.
  subroutine read_namelist(mdin_unit)
    use amber_rism_interface
    implicit none
    !mdin_unit :: unit number for mdin file
    character(len=8) :: closure(nclosuredefault)
    _REAL_ :: tolerance(nclosuredefault)
    integer :: closureOrder
    integer, intent(in) :: mdin_unit
    logical :: asympCorr
    integer :: entropicDecomp
    integer :: polarDecomp
    integer :: gfCorrection
    integer :: pcplusCorrection
    character(len=8) :: periodic
    _REAL_ :: uccoeff(size(rismprm%uccoeff))
    _REAL_ :: biasPotential
    _REAL_ :: solvcut
    _REAL_ :: buffer
    _REAL_ :: grdspc(3)
    integer ::  ng3(3)
    _REAL_ :: solvbox(3)
    _REAL_ :: mdiis_del
    integer :: mdiis_nvec
    integer :: mdiis_method
    _REAL_ :: mdiis_restart
    integer :: maxstep
    integer :: npropagate
    logical :: stoponfail
    integer :: centering
    integer :: zerofrc
    integer :: apply_rism_force
    integer :: rismnrespa
    integer :: fcestride
    _REAL_ :: fcecut
    integer ::  fcenbasis
    integer ::  fcenbase
    integer ::  fcecrd
    integer ::  fceweigh
    integer ::  fcetrans
    integer ::  fcesort
    integer ::  fceifreq
    integer ::  fcentfrcor
    integer ::  fcewrite
    integer ::  fceread
    _REAL_  ::  fceenormsw
    integer :: saveprogress
    logical :: molReconstruct
    integer :: ntwrism
    integer :: verbose
    integer :: progress
    logical :: selftest
    logical :: treeDCF
    logical :: treeTCF
    logical :: treeCoulomb
    _REAL_ :: treeDCFMAC
    _REAL_ :: treeTCFMAC
    _REAL_ :: treeCoulombMAC
    integer :: treeDCFOrder
    integer :: treeTCFOrder
    integer :: treeCoulombOrder
    integer :: treeDCFN0
    integer :: treeTCFN0
    integer :: treeCoulombN0
    _REAL_ :: asympKSpaceTolerance
    _REAL_ :: ljTolerance
    _REAL_ :: chargeSmear
#ifdef SANDER
    integer :: write_thermo
#endif
    namelist /rism/ &
         ! closure
         closure, closureOrder, uccoeff, &
         entropicDecomp, polarDecomp, &
         gfCorrection, pcplusCorrection, &
         biasPotential, &
         ! thermodynamics
         asympCorr, periodic, &
         ! solvation box
         buffer, grdspc, solvcut, &
         ng3, solvbox, &
         ! convergence
         tolerance, mdiis_del, mdiis_nvec, mdiis_method, &
         mdiis_restart, maxstep, npropagate, stoponfail, &
         ! minimization
         centering, zerofrc, &
         ! imin=5
         apply_rism_force, pa_orient, rmsd_orient, &
         ! md
         rismnrespa, &
         ! treecode and direct sum
         treeDCF, treeTCF, treeCoulomb, &
         treeDCFMAC, treeTCFMAC, treeCoulombMAC, &
         treeDCFOrder, treeTCFOrder, treeCoulombOrder, &
         treeDCFN0, treeTCFN0, treeCoulombN0, &
         asympKSpaceTolerance, chargeSmear, &
         ljTolerance, &
         molReconstruct, &
#ifdef RISM_CRDINTERP
         fcestride, fcecut, fcenbasis, fcenbase, fcecrd, &
         fceweigh,fcetrans,fcesort,fceifreq,fcentfrcor, &
         fceenormsw, fcewrite, fceread, &
#endif /*RISM_CRDINTERP*/
         !output
#ifdef SANDER
         write_thermo, &
#endif
         saveprogress, ntwrism, verbose, progress, volfmt, selftest
    
    call flush(0)

    closure = closurelist
    tolerance = tolerancelist
    closureOrder = rismprm%closureOrder
    asympCorr = rismprm%asympCorr
    entropicDecomp = rismprm%entropicDecomp
    polarDecomp = rismprm%polarDecomp
    gfCorrection = rismprm%gfCorrection
    pcplusCorrection = rismprm%pcplusCorrection
    periodic = periodicPotential
    uccoeff = rismprm%uccoeff
    biasPotential = rismprm%biasPotential
    solvcut = rismprm%solvcut
    buffer = rismprm%buffer
    grdspc= rismprm%grdspc
    ng3 = rismprm%ng3
    solvbox = rismprm%solvbox
    mdiis_del = rismprm%mdiis_del
    mdiis_nvec = rismprm%mdiis_nvec
    mdiis_method = rismprm%mdiis_method
    mdiis_restart = rismprm%mdiis_restart
    maxstep = rismprm%maxstep
    npropagate = rismprm%npropagate
    stoponfail = rismprm%stoponfail
    centering = rismprm%centering
    zerofrc = rismprm%zerofrc
    apply_rism_force = rismprm%apply_rism_force
    rismnrespa = rismprm%rismnrespa
    fcestride = rismprm%fcestride
    fcecut = rismprm%fcecut
    fcenbasis = rismprm%fcenbasis
    fcenbase = rismprm%fcenbase
    fcecrd = rismprm%fcecrd
    fceweigh = rismprm%fceweigh
    fcetrans = rismprm%fcetrans
    fcesort  = rismprm%fcesort
    fceifreq = rismprm%fceifreq
    fcentfrcor = rismprm%fcentfrcor
    fceenormsw = rismprm%fceenormsw
    fcewrite = rismprm%fcewrite
    fceread = rismprm%fceread
    saveprogress = rismprm%saveprogress
    molreconstruct = rismprm%molReconstruct
    ntwrism = rismprm%ntwrism
    verbose= rismprm%verbose
    progress = rismprm%progress
    selftest = rismprm%selftest
    treeDCF     = rismprm%treeDCF    
    treeTCF     = rismprm%treeTCF    
    treeCoulomb = rismprm%treeCoulomb
    treeDCFMAC     = rismprm%treeDCFMAC    
    treeTCFMAC     = rismprm%treeTCFMAC     
    treeCoulombMAC = rismprm%treeCoulombMAC
    treeDCFOrder     = rismprm%treeDCFOrder    
    treeTCFOrder     = rismprm%treeTCFOrder    
    treeCoulombOrder = rismprm%treeCoulombOrder
    treeDCFN0     = rismprm%treeDCFN0    
    treeTCFN0     = rismprm%treeTCFN0    
    treeCoulombN0 = rismprm%treeCoulombN0
    asympKSpaceTolerance = rismprm%asympKSpaceTolerance
    ljTolerance = rismprm%ljTolerance
    chargeSmear = rismprm%chargeSmear
#ifdef SANDER
    write_thermo = rismprm%write_thermo
#endif

!!$  !resize tolerance to the size of closure
!!$  tolerancelist => safemem_realloc(tolerancelist, size(closurelist))
!!$  tolerancelist(2:)=HUGE(1d0)

    rewind(mdin_unit)
    read(mdin_unit, nml=rism)

    closurelist = closure
    tolerancelist = tolerance
    rismprm%closureOrder = closureOrder
    rismprm%asympCorr=asympCorr
    rismprm%entropicDecomp = entropicDecomp
    rismprm%polarDecomp = polarDecomp
    rismprm%gfCorrection = gfCorrection
    rismprm%pcplusCorrection = pcplusCorrection
    periodicPotential = periodic
    rismprm%uccoeff = uccoeff
    rismprm%biasPotential = biasPotential
    ! Solvation box.
    rismprm%buffer=buffer
    rismprm%grdspc=grdspc
    rismprm%solvcut=solvcut
    rismprm%ng3=ng3
    rismprm%solvbox=solvbox
    ! Convergence.
    rismprm%mdiis_del=mdiis_del
    rismprm%mdiis_nvec=mdiis_nvec
    rismprm%mdiis_method=mdiis_method
    rismprm%mdiis_restart=mdiis_restart
    rismprm%maxstep=maxstep
    rismprm%npropagate=npropagate
    rismprm%stoponfail = stoponfail
    ! Minimization.
    ! We need to know if this is a periodic or open boundary to set the default centering.
    if (centering == -9999) then
       if (periodicPotential /= '') then
          ! For periodic calculations, the default is 0 (no centering)
          rismprm%centering = 0
       else
          ! For open-boundary calculations, the default is 1 (CoM centering)
          rismprm%centering = 1
       end if
    else
       ! use the user preferences
       rismprm%centering = centering
    end if
    rismprm%zerofrc=zerofrc
    ! imin=5
    rismprm%apply_rism_force=apply_rism_force
    pa_orient=pa_orient
    rmsd_orient=rmsd_orient
    ! md
    rismprm%rismnrespa=rismnrespa
#ifdef RISM_CRDINTERP
    rismprm%fcestride=fcestride
    rismprm%fcecut=fcecut
    rismprm%fcenbasis=fcenbasis
    rismprm%fcenbase=fcenbase
    rismprm%fcecrd=fcecrd
    rismprm%fceweigh=fceweigh
    rismprm%fcetrans=fcetrans
    rismprm%fcesort=fcesort
    rismprm%fceifreq=fceifreq
    rismprm%fcentfrcor=fcentfrcor
    rismprm%fceenormsw=fceenormsw
    rismprm%fcewrite=fcewrite
    rismprm%fceread=fceread
#endif /*RISM_CRDINTERP*/
    ! Output.
    rismprm%saveprogress=saveprogress
    rismprm%molReconstruct = molReconstruct
    rismprm%ntwrism=ntwrism
    rismprm%verbose=verbose
    rismprm%progress=progress
    rismprm%selftest=selftest
    rismprm%treeDCF     = treeDCF    
    rismprm%treeTCF     = treeTCF    
    rismprm%treeCoulomb = treeCoulomb
    rismprm%treeDCFMAC     = treeDCFMAC    
    rismprm%treeTCFMAC     = treeTCFMAC
    rismprm%treeCoulombMAC = treeCoulombMAC
    rismprm%treeTCFOrder     = treeTCFOrder    
    rismprm%treeDCFOrder     = treeDCFOrder    
    rismprm%treeCoulombOrder = treeCoulombOrder
    rismprm%treeTCFN0     = treeTCFN0    
    rismprm%treeDCFN0     = treeDCFN0    
    rismprm%treeCoulombN0 = treeCoulombN0
    rismprm%asympKSpaceTolerance = asympKSpaceTolerance
    rismprm%ljTolerance = ljTolerance
    rismprm%chargeSmear = chargeSmear
#ifdef SANDER
    rismprm%write_thermo = write_thermo
#endif

    ! Set the RISM cutoff if not set by the user.
    if (rismprm%solvcut < 0) then
       rismprm%solvcut = rismprm%buffer
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

    if (periodicPotential /= '') then
       if (periodicPotential /= "ewald" & 
             .and. periodicPotential /= "ewaldn"  &
             .and. periodicPotential /= "pme"  &
             .and. periodicPotential /= "pmen" ) then
          call rism_report_error( &
          "Only 'ewald', 'ewaldn', 'pme' and 'pmen' periodic potentials are supported.")
       end if 
       ! Do not apply asymptotic corrections for periodic 3D-RISM.
       rismprm%asympcorr = .false.
       !if (rismprm%apply_rism_force /= 0 ) call rism_report_error( &
       !   "RISM forces are not yet correct for periodic systems")
    end if

    ! Ensure that solvbox and ng3 have been set if buffer < 0.
    if (rismprm%buffer < 0) then
       if (minval(rismprm%ng3) < 0) &
            call rism_report_error('if buffer < 0 ng3 (grid size) must be set.')
       if (minval(rismprm%solvbox) < 0) &
            call rism_report_error('if buffer < 0 solvbox must be set.')
       ! Ensure that the cutoff is set to a reasonable value. This can
       ! happen if buffer is < 0.
       if (rismprm%solvcut < 0 .and. periodicPotential /= '') then
          call rism_report_error('solvcut must be >= 0.')
       end if
    elseif (rismprm%buffer == 0d0 .and. rismprm%ljTolerance == 0d0) then
       call rism_report_error('buffer and ljTolerance cannot both be 0.')
    end if

    ! Ensure that an appropriate file format has been chosen for
    ! volumetric output.
    if (.not. (volfmt .eq. "ccp4" .or. volfmt .eq. "mrc" .or. volfmt .eq. "dx" .or. volfmt .eq. "xyzv")) then
       call rism_report_error("Only 'mrc', 'ccp4', 'dx', and 'xyzv' volumetric data formats are supported")
    end if

    ! Make sure that a valid centering method is used and doesn't
    ! conflict with other options.
    if (rismprm%centering > 4 .or. rismprm%centering < -4) &
         call rism_report_error("Centering must be between -4 and 4")
    if ((rismprm%centering > 2 .or. rismprm%centering < -2) .and. rismprm%fcestride > 0) &
         call rism_report_error("CENTERING must be between -2 and 2 when FCESTRIDE > 0")

    ! We can only run temperature derivative if the Xvv file has this information.
    if (rismprm%entropicDecomp == 1 .and. .not.rism3d_solvent_canCalc_dT(solvent)) &
         call rism_report_warn("cannot perform temperature derivative. No 1D-RISM temperature derivative data found")

    ! Solvation energy and entropy maps only make sense for entropicDecomp runs.
    if (rismprm%entropicDecomp /=1 .and. &
         (len_trim(solvationEnergyfile) /= 0 .or. len_trim(entropyfile) /= 0)) then
       call rism_report_error(&
            "Solvation energy and entropy maps require energy/entropy decomposition")
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
               call rism_report_error("'closureOrder' is deprecated and not compatible "&
               //"with closure lists")
          write(fmt, '(a, i4, a)') '(a, i', int(log10(dble(rismprm%closureOrder))) + 1, ')'
          write(closurelist, fmt) "PSE", rismprm%closureOrder
       end if
    end do

    ! check treecode parameters
    if (rismprm%treeDCFMAC .lt.  0.) then
       call rism_report_error('(a,e16.8e3)',"'treeDCFMAC' must >= 0. Got ",rismprm%treeDCFMAC)
    end if
    if (rismprm%treeTCFMAC .lt.  0.) then
       call rism_report_error('(a,e16.8e3)',"'treeTCFMAC' must >= 0. Got ",rismprm%treeTCFMAC)
    end if
    if (rismprm%treeCoulombMAC .lt.  0.) then
       call rism_report_error('(a,e16.8e3)',"'treeCoulombMAC' must >= 0. Got",rismprm%treeCoulombMAC)
    end if
    if (rismprm%treeDCFOrder .lt.  0.) then
       call rism_report_error('(a,i8)',"'treeDCFOrder' must >= 0. Got ",rismprm%treeDCFOrder)
    end if
    if (rismprm%treeTCFOrder .lt.  0.) then
       call rism_report_error('(a,i8)',"'treeTCFOrder' must >= 0. Got ",rismprm%treeTCFOrder)
    end if
    if (rismprm%treeCoulombOrder .lt.  0.) then
       call rism_report_error('(a,i8)',"'treeCoulombOrder' must >= 0. Got",rismprm%treeCoulombOrder)
    end if
    if (rismprm%treeDCFN0 .lt.  1.) then
       call rism_report_error('(a,i8)',"'treeDCFN0' must >= 1. Got ",rismprm%treeDCFN0)
    end if
    if (rismprm%treeTCFN0 .lt.  1.) then
       call rism_report_error('(a,i8)',"'treeTCFN0' must >= 1. Got ",rismprm%treeTCFN0)
    end if
    if (rismprm%treeCoulombN0 .lt.  1.) then
       call rism_report_error('(a,i8)',"'treeCoulombN0' must >= 1. Got",rismprm%treeCoulombN0)
    end if

    ! long-range asymptotics error tolerance
    if (rismprm%asympKSpaceTolerance .lt. 0d0 .and. rismprm%asympKSpaceTolerance .ne. -1d0) then
       call rism_report_error('(a,e10.2)',"'asympKSpaceTolerance' must be >= 0 or -1. Got",rismprm%asympKSpaceTolerance)
    end if

    ! Lennard-Jones potential error tolerance
    if (rismprm%ljTolerance .lt. 0d0 .and. rismprm%ljTolerance .ne. -1d0) then
       call rism_report_error('(a,e10.2)',"'ljTolerance' must be >= 0 or -1. Got",rismprm%ljTolerance)
    end if

    ! charge smear
    if (rismprm%chargeSmear .lt. 0d0) then
       call rism_report_error('(a,g)',"'chargeSmear' must be >= 0. Got",rismprm%chargeSmear)
    end if
#ifdef MPI
    if (len_trim(electronMapFile)>0) then
       call rism_report_error('Electron density maps are not supported for MPI calculations')
    end if
#endif    
  end subroutine sanity_check

  !> Determines if several temperature derivative based quantities
  !! (PMV_dT, PCPLUS, UCT, NgBT, etc) can be calculated from the given
  !! input and change to internal units
  !! @sideeffects
  !!   sets canCalculateUC and canCalculatePartialMolarVolume_dT.
  subroutine temperatureDerivativeCheck()
    use amber_rism_interface
    use constants, only : kb
    implicit none
    ! check if the user has supplied UC coefficients
    if (all(rismprm%uccoeff == 0d0)) then
       canCalculateUC = .false.
    else
       canCalculateUC = .true.
       ! Unit conversion for user supplied parameters.
       rismprm%uccoeff = rismprm%uccoeff / (KB * rism3d_get_temperature(rism_3d))
    end if

    ! check if we can calculate temperature derivatives for PMV,
    ! needed for PCPLUS and UC
    if (rism3d_get_xikt_dt(rism_3d) == huge(1d0)) then
       canCalculatePartialMolarVolume_dT = .false.
    else
       canCalculatePartialMolarVolume_dT = .true.
    end if
       
    ! Warn the user that some values cannot be calculate.
    if (rismprm%entropicDecomp == 1 .and. rism3d_canCalc_DT(rism_3d) &
         .and. .not. canCalculatePartialMolarVolume_dT) then
       call rism_report_warn("Input Xvv file version < 1.001. "&
            //"Cannot calculate partialMolarVolume_dT, or solvation energy or entropy of UC and PCPLUS.")
    end if
  end subroutine temperatureDerivativeCheck
  
  !> Modifies the position of the solute using the centering method
  !! requested. (This used to also rotate the solute in the box but
  !! this has been disabled.)
  !! If centering is /= 0 we want to move the solute to the center of
  !! the solvent box. However, if centering < 0 we only figure out the
  !! displacement required the _first_ time we see the solute. Thus,
  !! for centering <= 0, the solute's CM can move relative to the
  !! grid.
  !! @param[in] solu solute object.
  !! @param[in,out] atomPositions The x, y, z position of each solute atom.
  !!    This is modified to place the solute in the center of the box.
  !! @param[in] nsolution Number of times a RISM solution has been
  !!    calculated. Methods < 0 only are used if nsolution == 0.
  subroutine orient(mySolute, atomPositions, nsolution)
    use amber_rism_interface
    use rism_util, only : findCenterOfMass, translate
    use constants, only : PI
    implicit none
    type(rism3d_solute) :: mySolute
    _REAL_, intent(inout) :: atomPositions(3, mySolute%numAtoms)
    integer, intent(in) :: nsolution
    _REAL_ :: weight(mySolute%numAtoms)
    
#ifdef RISM_DEBUG
    write(outunit, *) "RISM REORIENT"
    call flush(6)
#endif /*RISM_DEBUG*/

    if (abs(rismprm%centering) == 1) then
       weight = mySolute%mass
    else if (abs(rismprm%centering) == 2) then
       weight = 1
    end if
    if (rismprm%centering > 0 .or. (rismprm%centering < 0 .and. nsolution == 0)) then
       call findCenterOfMass(atomPositions, centerOfMass, weight, mySolute%numAtoms)
    end if
    if (rismprm%centering /= 0) then
       call translate(atomPositions, mySolute%numAtoms, -1 * centerOfMass)
    end if
  end subroutine orient


  !! Return the solute to its original position.
  !! IN:
  !!    solute :: solute object
  !!    atomPositions  :: the x, y, z position of each solute atom.  This is modified
  !!             to return the solute to its original MD position and orientation
  subroutine unorient(mySolute, atomPositions)
    use amber_rism_interface
    use rism_util, only: findCenterOfMass, translate
    implicit none
    type(rism3d_solute) :: mySolute
    _REAL_, intent(inout) :: atomPositions(3, mySolute%numAtoms)

    _REAL_ :: cm(3), weight(mySolute%numAtoms)
    !counters
    integer :: iu, id

    ! We don't know if additional translations have been performed on
    ! the mySolute.  1) The first step is to calculate the current COM
    ! and translate this to the origin.  2) Then the system and forces
    ! are rotated according to qback.  3) Finally, using centerOfMass,
    ! we translate the system back to its original MD COM.

#ifdef RISM_DEBUG
    write(outunit, *)"RISM UNORIENT", centerOfMass
    ! write(outunit, *) qback
    call flush(6)
#endif /*RISM_DEBUG*/
    if (abs(rismprm%centering)==1) then
       weight = mySolute%mass
       call findCenterOfMass(atomPositions, cm, weight, mySolute%numAtoms)
       call translate(atomPositions, mySolute%numAtoms, -1*cm)
       call findCenterOfMass(atomPositions, cm, weight, mySolute%numAtoms)
    else if (abs(rismprm%centering)==2) then
       weight=1
       call findCenterOfMass(atomPositions, cm, weight, mySolute%numAtoms)
       call translate(atomPositions, mySolute%numAtoms, -1*cm)
       call findCenterOfMass(atomPositions, cm, weight, mySolute%numAtoms)
    end if
    if (rismprm%centering /= 0) then
       call translate(atomPositions, mySolute%numAtoms, centerOfMass)
    end if

  end subroutine unorient


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


#ifdef SANDER
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
#else /*SANDER*/
!!!stubs for SANDER timers for non-SANDER executables
  subroutine timer_start( label )
    integer label
  end subroutine timer_start

  subroutine timer_stop( label )
    integer label
  end subroutine timer_stop
#endif /*SANDER*/

#ifdef SANDER
end module sander_rism_interface
#endif
