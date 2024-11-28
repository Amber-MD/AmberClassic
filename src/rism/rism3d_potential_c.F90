!<compile=optimized>

#include "../include/dprec.fh"

!> Electrostatic potential class for 3D-RISM. Used to calculate/store
!! quantities that are potential dependent and do not change while
!! converging a solution calculation.
!!
!! Pointers to solute and solvent objects are maintained.  So, if
!! values change in these objects, they are automatically used in the
!! potential calculation.
!!
!! This class is generally MPI agnostic.  That is, MPI is only an
!! issue for setting the grid size where both information about the
!! size of the local slab and the total grid must be supplied.  There
!! is no MPI communication within the class.

module rism3d_potential_c
  use safemem
  use rism3d_solute_c
  use rism3d_solvent_c
  use rism3d_grid_c
  use rism3d_fft_c ! Only used for periodic code.
#include "def_time.h"
  type rism3d_potential
     !> Only called by periodic code.
     type(rism3d_fft), pointer :: fft => NULL()

     !> Cutoff for RISM potential calculations.
     _REAL_ :: cutoff
     !> Cutoff**2 for RISM LJ potential calculations.
     _REAL_ :: cutoff2
     _REAL_, pointer :: ljCutoffs2(:,:) => NULL()

     !> Pointer to grid object.
     type(rism3d_grid), pointer :: grid => NULL()

     !> Pointer to solute object.
     type(rism3d_solute), pointer :: solute => NULL()
     !> Pointer to solvent object.
     type(rism3d_solvent), pointer :: solvent => NULL()

     !> Potential energy of the solvent about the solute.  This is
     !! recalculated for each solution but we want to reserve the
     !! memory and may want to change it if the box dimensions
     !! change. The fourth dimension is per solvent atom. [kT]
     _REAL_, pointer :: uuv(:,:,:,:) => NULL()

     !> The long range portion of the Ewald potential in k-space.
     !! Every two items are the real and imaginary components of the
     !! complex exponential. Thus even though only half the box points
     !! have their potential calculated, the array must have
     !! dimensions of the total number of k-space grid points. This
     !! also allows it to store the real space output which is for
     !! every grid point.
     _REAL_, pointer :: uuv1d(:,:) => NULL()

     !> Solute-solvent sigma interaction matrix.  Calculated once and
     !! used in later calls. [A]
     _REAL_, pointer :: ljSigmaUV(:,:) => NULL()
     !> Solute-solvent epsilon interaction matrix.  Calculated once and
     !! used in later calls. [kT]
     _REAL_, pointer :: ljEpsilonUV(:,:) => NULL()

     !> Long-range part of Huv(k) at k = 0 (2, solv%natom).
     _REAL_, pointer :: huvk0(:,:) => NULL()

     !> If true, a periodic 3D-RISM calculation is performed. This
     !! mostly differs from simple 3D-RISM by using Ewald sum potential
     !! in place of Coulombic potential and invoking the minimum image
     !! convention while calculating potentials.
     logical :: periodic = .false.

     !> Specifies the periodic potential to use.
     !! Current valid values include:
     !!   'pme'   = Particle Mesh Ewald potential
     character(len=256) :: periodicPotential = ''

     !> Charge smearing parameter for long-range
     !! asymtotics and Ewald, typically eta in the literature
     _REAL_ :: chargeSmear

     !> indicates if it is approriate to apply the LJ truncation
     !! correction for thermodynamic values
     logical :: applyLJCorrection = .false.
  end type rism3d_potential
  
  private mixSoluteSolventLJParameters
  private uvCoulombicPotential
  private uvLennardJonesPotentialWithCutoff 
  private uvLJrEwaldMinImage, uvPMErecip
contains


  !> Constructor.
  !! @param[in,out] this potential object
  !! @param[in] grid grid object.  A pointer to this will be retained.
  !! @param[in] solv solvent object.  A pointer to this will be retained.
  !! @param[in] solu solute object.  A pointer to this will be retained.
  !! @param[in] cut Cutoff.
  !! @param[in] fft Fast Fourier Transform object.
  !! @param[in] periodic True when calculating potentials for periodic solute.
  !! @param[in] chargeSmear :: Charge smearing parameter for long-range
  !!       asymtotics and Ewald, typically eta in the literature
  subroutine rism3d_potential_new(this, grid, solv, solu, cut, fft, &
       periodicPotential, chargeSmear)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    type(rism3d_grid), target, intent(in) :: grid
    type(rism3d_solute), target, intent(in) :: solu
    type(rism3d_solvent), target, intent(in) :: solv
    _REAL_, intent(in):: cut
    type(rism3d_fft), target, intent(in) :: fft
    character(len=*), intent(in) :: periodicPotential
    _REAL_, intent(in) :: chargeSmear
    this%grid => grid
    this%solvent => solv
    this%solute => solu
    this%fft => fft
    this%periodicPotential = periodicPotential
    if (this%periodicPotential /= '') then
       this%periodic = .true.
    end if
    this%chargeSmear = chargeSmear
    this%ljCutoffs2 => safemem_realloc(this%ljCutoffs2, this%solute%numAtoms, this%solvent%numAtomTypes, .false.)
    this%ljSigmaUV => safemem_realloc(this%ljSigmaUV, this%solute%numAtoms, this%solvent%numAtomTypes, .false.)
    this%ljEpsilonUV => safemem_realloc(this%ljEpsilonUV, this%solvent%numAtomTypes,this%solute%numAtoms, .false.)
    call rism3d_potential_setCut_ljdistance(this, cut)
    call mixSoluteSolventLJParameters(this)
    
  end subroutine rism3d_potential_new

  !> Directly a distance cut off for potential and force.
  !! @param[in,out] this potential object.
  !! @param[in] cut Distance cutoff for potential and force calculations.
  subroutine rism3d_potential_setcut_ljdistance(this, cut)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(in) :: cut

    ! Assign potential cutoff.
    ! Ensure that the square won't overflow.
    this%cutoff = min(sqrt(huge(1d0)), cut)
    this%cutoff2 = this%cutoff**2

    ! non-periodic code now uses this variable
    this%ljCutoffs2 = this%cutoff**2

  end subroutine rism3d_potential_setcut_ljdistance
  
  !> Calculates the potential on the grid.
  subroutine rism3d_potential_calc(this)
    use rism_util, only : checksum
    use rism3d_opendx, only : rism3d_opendx_write
#ifdef MPI
    use mpi
#endif
    implicit none
    type(rism3d_potential), intent(inout) :: this !< potential object.

    integer :: id

    integer :: ix, iy, iz, ierr
    character(len=30) :: filename
    integer :: iu 
    logical, save :: first = .true.

    ! Ensure a grid size has been set.
    if (.not. associated(this%grid%waveVectors)) then
       call rism_report_error("rism3d_potential_calc: grid size not set")
       stop
    end if
    ! Check if the grid size has changed.
    if (ubound(this%uuv, 1) /= this%grid%localDimsR(1) .or. &
         ubound(this%uuv, 2) /= this%grid%localDimsR(2) .or. &
         ubound(this%uuv, 3) /= this%grid%localDimsR(3) .or. &
         .not. associated(this%uuv)) then
       this%uuv => safemem_realloc(this%uuv, &
            this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3),&
            this%solvent%numAtomTypes, .false.)
       this%uuv1d => safemem_realloc(this%uuv1d, this%grid%totalLocalPointsK, &
            this%solvent%numAtomTypes, o_preserve = .false., o_aligned = .true.)
    end if
     
    this%uuv = 0

    call timer_start(TIME_UCOULU)
    if (this%solute%charged) call uvPMErecip(this, this%uuv)
    call timer_stop(TIME_UCOULU)

    call timer_start(TIME_ULJUV)
    call uvLJrEwaldMinImage(this, this%uuv)
    call timer_stop(TIME_ULJUV)

  end subroutine rism3d_potential_calc

  !> Frees all memory and resets values.
  !! @param[in,out] this potential object.
  subroutine rism3d_potential_destroy(this)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    nullify(this%grid)
    nullify(this%solvent)
    nullify(this%solute)
    this%cutoff2 = 0
    if (safemem_dealloc(this%uuv) /= 0) &
         call rism_report_error("Uuv deallocation failed")
    if (safemem_dealloc(this%uuv1d,o_aligned=.true.) /= 0) &
         call rism_report_error("Uuv1d deallocation failed")
    if (safemem_dealloc(this%ljSigmaUV) /= 0) &
         call rism_report_error("LjSigmaUV deallocation failed")
    if (safemem_dealloc(this%ljEpsilonUV) /= 0) &
         call rism_report_error("EPSuv deallocation failed")
    if (safemem_dealloc(this%ljCutoffs2) /= 0) &
         call rism_report_error("ljCutoffs2 deallocation failed")
    if (safemem_dealloc(this%huvk0) /= 0) &
         call rism_report_error("huvk0 deallocation failed")
  end subroutine rism3d_potential_destroy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                               PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !> Mix per-site Lennard-Jones solute and solvent parameters to obtain LJ
  !! solute-solvent interaction parameters.
  !! @param[in,out] this potential object.
  subroutine mixSoluteSolventLJParameters(this)
    use safemem
    implicit none
    type(rism3d_potential), intent(inout) :: this
    integer :: iv, iu
    !GMG> where to write, read lj matrix
    integer, parameter :: out_unit=20
    integer, parameter :: sout_unit=22
    integer, parameter :: in_unit=21
    logical :: exist
    integer :: status
    _REAL_ :: sm,em
    !-------------------------------

    ! Read LJ-solute (uu) corrections in case ljsolute.mods.txt file is
    ! present.
    !--------------------------------------------------------------------
    inquire(file="ljsolute.mods.txt", exist=exist)
    if (exist) then
        open(in_unit, file="ljsolute.mods.txt", status="old", action="read")
        do
          read (in_unit,'(i6,2e16.8)',IOSTAT=status) iu,sm,em
          if (status /= 0) exit
          ! write (6,'(a5,i6,2e16.8)') "ljsolute-mods>",iu,sm,em
          ! call
          ! rism_report_message("(a20,i6,2e16.8)","ljsolute-mods>",iu,sm,em)
          this%solute%ljSigma(iu)   = sm
          this%solute%ljEpsilon(iu) = em
        enddo
        call rism_report_message("|-> Read LJ-solute modifications from ljsolute.mods.txt.")

       ! Write the LJ-solute (uu) parameters.
       !--------------------------------------------------------------------
       open (unit=sout_unit,file="ljsolute.orig.txt",action="write",&
             status="replace")
       do iu = 1, this%solute%numAtoms
          write (sout_unit,'(i6,2e16.8)') iu,this%solute%ljSigma(iu),&
                                             this%solute%ljEpsilon(iu)
       end do
       close(sout_unit)
    endif

    ! Compute the LJ-matrix (uv).
    ! ----------------------------------------------------------
    do iv = 1, this%solvent%numAtomTypes
       do iu = 1, this%solute%numAtoms
          this%ljSigmaUV(iu, iv) = this%solute%ljSigma(iu) + &
                                   this%solvent%ljSigma(iv)
          this%ljEpsilonUV(iv, iu) = sqrt(this%solute%ljEpsilon(iu) * &
                                          this%solvent%ljEpsilon(iv))
       end do
    end do

    ! Read LJ-matrix (uv) corrections in case ljmatrix.mods.txt file is present.
    !--------------------------------------------------------------------
    !check file exist to trigger the LJ-matrix update
    inquire(file="ljmatrix.mods.txt", exist=exist)
    !if true, than read until eof has been reached
    if (exist) then
       open(in_unit, file="ljmatrix.mods.txt", status="old", action="read")
       do
         read (in_unit,'(2i6,2e16.8)',IOSTAT=status) iu,iv,sm,em
         if (status /= 0) exit
         !write (6,'(a5,2i6,2e16.8)') "mods>", iu,iv,sm,em
         !call rism_report_message("(a20,2i6,2e16.8)", "lj-matrix-mods>",
         !iu,iv,sm,em)
         this%ljSigmaUV(iu,iv) = sm
         this%ljEpsilonUV(iv,iu) = em
       enddo
       call rism_report_message("|-> Read LJ-matrix modifications from ljmatrix.mods.txt.")


       ! Write the LJ-matrix (uv) to ljmatrix.orig.txt.
       ! ----------------------------------------------------------
       open (unit=out_unit,file="ljmatrix.orig.txt",action="write",&
             status="replace")
       do iv = 1, this%solvent%numAtomTypes
          do iu = 1, this%solute%numAtoms
             write (out_unit,'(2i6,2e16.8)') iu,iv,this%ljSigmaUV(iu,iv), &
                                                   this%ljEpsilonUV(iv,iu)
          end do
       end do
       close (out_unit)
    endif
    flush(6)
  end subroutine mixSoluteSolventLJParameters

  !> Tabulate the solute-solvent 12-6 Lennard-Jones potential, and
  !! the short-range Ewald electrostatic potential,  in the
  !! box subject to the minimum image convention.
  !! @param[in] this potential object
  !! @param[in,out] ulj grid to add potential to
  subroutine uvLJrEwaldMinImage(this, ulj)
    use constants_rism, only : omp_num_threads
#ifdef MPI
    use mpi
#endif
    implicit none
    type(rism3d_potential), intent(in) :: this
    _REAL_, intent(inout) :: ulj(:,:,:,:)
    ! Grid, solute, slovent, and dimension indices.
    integer :: igx, igy, igz, iu, iv
    ! Grid point position.
    _REAL_ :: rx(3), ry(3), rz(3), gridp(3)
    ! Distance of grid point from solute.
    _REAL_ :: sd2
    ! Base term in LJ equation (ratio of sigma and distance).
    _REAL_ :: ljBaseTerm(this%solvent%numAtomTypes)
    _REAL_ :: solutePosition(3)
    _REAL_ :: sd, sr
    _REAL_ :: sigma(this%solvent%numAtomTypes,this%solute%numAtoms), beta

    beta = 1.d0/this%chargeSmear
    do iu = 1, this%solute%numAtoms
       do iv = 1, this%solvent%numAtomTypes
          sigma(iv,iu) = 1.d0/this%ljSigmaUV(iu, iv)**2
       end do
    end do


!$omp parallel do private (rx,ry,rz,solutePosition,sd2,sd,sr,ljBaseTerm, &
!$omp&   igx,igy,igz,iu) num_threads(omp_num_threads)

    do igz = 1, this%grid%localDimsR(3)
       rz = (igz - 1 + this%grid%offsetR(3)) * this%grid%voxelVectorsR(3, :)
       do igy = 1, this%grid%localDimsR(2)
          ry = (igy - 1) * this%grid%voxelVectorsR(2, :)
          do igx = 1, this%grid%localDimsR(1)
             rx = (igx - 1) * this%grid%voxelVectorsR(1, :)

             do iu = 1, this%solute%numAtoms

                solutePosition = rx + ry + rz - this%solute%position(:, iu)
                solutePosition = minimumImage(this, solutePosition)

                sd2 = max(4d-6, dot_product(solutePosition, solutePosition))
                if (sd2 < this%cutoff2) then

                   sd = sqrt(sd2)
                   sr = this%solute%charge(iu) * erfc(sd * beta) / sd
                   ljBaseTerm(:) = 1.d0 / (sd2 * sigma(:,iu))**3

                   ulj(igx,igy,igz,:) = ulj(igx,igy,igz,:) &
                           + this%ljEpsilonUV(:,iu) &
                             * ljBaseTerm(:) * (ljBaseTerm(:) - 2.d0) &
                           + sr * this%solvent%charge(:)
                end if
             end do
          end do
       end do
    end do
!$omp end parallel do
  end subroutine uvLJrEwaldMinImage

  !> Long-range portion of the Particle Mesh Ewald (PME) electric potential.
  subroutine uvPMErecip(this, ucu)
    use, intrinsic :: iso_c_binding
    use bspline
    use constants_rism, only : pi
    use FFTW3
    use rism3d_opendx, only : rism3d_opendx_write
    use rism_util, only: r2c_pointer
#ifdef MPI
    use mpi
#endif
    implicit none
    type(rism3d_potential), intent(inout) :: this
    _REAL_, intent(inout) :: ucu(:,:,:,:)

    ! Ewald charge smear parameter.
    ! _REAL_, parameter :: smear = this%chargeSmear
    ! Order of b-spline interpolation.
    integer, parameter :: splineOrder = 6

    !_REAL_, parameter :: zeta = this%chargeSmear * this%chargeSmear
    !_REAL_, parameter :: t = -0.25 / zeta

    _REAL_ :: smear, zeta, t

    logical, parameter :: DEBUG = .false.

    integer :: i, j, k
    integer :: ix, iy, iz, ixyz
    integer :: igx, igy, igz
    integer :: igk, igk_fort, igk_cpp, igk_fort_global, igk_cpp_global
    integer :: id
    integer :: iu, iv

    integer :: gridPoints(splineOrder, 3)

    _REAL_, pointer :: kxi(:,:) => NULL()
    _REAL_, pointer :: kyi(:,:) => NULL()
    _REAL_, pointer :: kzi(:,:) => NULL()
    _REAL_ :: k2
    _REAL_ :: waveVector(3)
    integer :: lgx, lgy, lgz

    _REAL_ :: byz, bxyz

    _REAL_ :: weights(splineOrder, 3)
    _REAL_, pointer :: bsplineFourierCoeffX(:) => NULL()
    _REAL_, pointer :: bsplineFourierCoeffY(:) => NULL()
    _REAL_, pointer :: bsplineFourierCoeffZ(:) => NULL()
    _REAL_, pointer :: gaussianFourierCoeff(:) => NULL()
    _REAL_, pointer :: kernel(:) => NULL()

    type(C_PTR) :: uuv1d_r_cptr, uuv1d_c_cptr, uuv1d_final_cptr
    real(C_DOUBLE), pointer :: uuv1d_r(:,:,:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: uuv1d_c(:,:,:) => NULL()
    real(C_DOUBLE), pointer :: uuv1d_final(:,:,:) => NULL()
    integer(C_INTPTR_T) :: localPtsK
    real(C_DOUBLE), pointer :: outr_1d(:) => NULL()
    integer(C_INTPTR_T) :: L, M, N
    integer(C_INTPTR_T) :: local_N, local_k_offset

    type(C_PTR) :: planfwd = C_NULL_PTR, planbwd = C_NULL_PTR

    _REAL_ :: reciprocalPos(3)
    _REAL_ :: chargeCorrection
    ! integer :: numGridPointsK

    integer :: gridDimX_k

    character(len=120) :: filename, suffix
    integer :: ierr

    smear = this%chargeSmear
    zeta = this%chargeSmear * this%chargeSmear
    t = -0.25 * zeta

    L = this%grid%localDimsR(1)
    M = this%grid%localDimsR(2)
    N = this%grid%globalDimsR(3)

!    call timer_start(TIME_UCOULULR)

    gridDimX_k = this%grid%localDimsR(1) / 2 + 1
    kxi => safemem_realloc(kxi,  gridDimX_k, 3, .false.)
    kyi => safemem_realloc(kyi,  this%grid%localDimsR(2), 3, .false.)
    kzi => safemem_realloc(kzi,  this%grid%localDimsR(3), 3, .false.)

    bsplineFourierCoeffX => safemem_realloc(bsplineFourierCoeffX, this%grid%globalDimsR(1), .false.)
    bsplineFourierCoeffY => safemem_realloc(bsplineFourierCoeffY, this%grid%globalDimsR(2), .false.)
    bsplineFourierCoeffZ => safemem_realloc(bsplineFourierCoeffZ, this%grid%globalDimsR(3), .false.)

    gaussianFourierCoeff => safemem_realloc(gaussianFourierCoeff, &
         gridDimX_k * this%grid%localDimsR(2) * this%grid%localDimsR(3), .false.)

#ifdef MPI
    ! localPtsK = fftw_mpi_local_size_3d_transposed( &
    localPtsK = fftw_mpi_local_size_3d(N, M, L / 2 + 1, &
         this%grid%mpicomm, local_N, local_k_offset)
         ! localPtsK_y, localOffsetK_y)
#else
    localPtsK = gridDimX_k * this%grid%localDimsR(2) * this%grid%localDimsR(3)
    local_N = this%grid%localDimsR(3)
    local_k_offset = 0
#endif

    uuv1d_r_cptr = fftw_alloc_real(2 * localPtsK)
    uuv1d_c_cptr = fftw_alloc_complex(localPtsK)
#ifdef MPI
    call c_f_pointer(uuv1d_r_cptr, uuv1d_r, [2*(L/2+1), M, local_N])
    call c_f_pointer(uuv1d_r_cptr, outr_1d, [2*(L/2+1) * M * local_N])
    call c_f_pointer(uuv1d_c_cptr, uuv1d_c, [L/2+1, M, local_N])
    if (this%grid%mpirank == 0) then
       uuv1d_final_cptr = fftw_alloc_real(L * M * N)
       call c_f_pointer(uuv1d_final_cptr, uuv1d_final, [L, M, N])
    end if
#else
    call c_f_pointer(uuv1d_r_cptr, uuv1d_r, [L, M, local_N])
    call c_f_pointer(uuv1d_r_cptr, outr_1d, [L * M * local_N])
    call c_f_pointer(uuv1d_c_cptr, uuv1d_c, [L/2+1, M, local_N])
#endif

    kernel => safemem_realloc(kernel, &
         gridDimX_k * this%grid%localDimsR(2) * this%grid%localDimsR(3), .false.)

#if 0
    if (any(mod(this%grid%localDimsR, 2) > 0)) then
       !TODO: call RISM error functions
#ifdef MPI
       if( this%grid%mpirank == 0 ) &
#endif
       write(6,'(a,a)')  "| PME implementation prefers an even", &
            " number of grid points on each axis."
    end if
#endif

   ! Angular wave numbers.
   !FIXME: This needs to be removed and X should be the halved axis.
   !FIXME: Scrutinize grid dims in loops (K v. R, local vs. global).
   do igk = 0, gridDimX_k - 1
      kxi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(1,:) * igk
   end do
   do igk = 0, this%grid%localDimsR(2) - 1
      kyi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(2,:) &
           * merge(igk, -(this%grid%localDimsR(2) - igk), &
           igk < this%grid%localDimsR(2) / 2 + 1)
   end do
   do igk = 0, this%grid%localDimsR(3) - 1
      kzi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(3,:) &
           * merge(igk + this%grid%offsetR(3), &
           -(this%grid%globalDimsR(3) - (igk + this%grid%offsetR(3))), &
           (igk + this%grid%offsetR(3)) < this%grid%globalDimsR(3) / 2 + 1)
   end do

    ! Compute the discrete Fourier transform coefficients of the b-spline.
    call cardinal_bspline(merge(0d0, 0.5d0, mod(splineOrder, 2) == 0), &
         splineOrder, weights(:,1))
    !TODO: For triclinic case, b-spline may have axial interdependence
    ! and thus so will its Fourier coefficients.

    !FIXME: Make these functions MPI aware by passing in global +
    ! local dims and offset. For now, just compute all coeffs on every
    ! node and use the grid offset to index.
    call cardinal_bspline_Fourier_coefficients( &
         this%grid%globalDimsR(1), splineOrder, &
         weights(:,1), bsplineFourierCoeffX, .true.)
    call cardinal_bspline_Fourier_coefficients( &
         this%grid%globalDimsR(2), splineOrder, &
         weights(:,1), bsplineFourierCoeffY, .false.)
    call cardinal_bspline_Fourier_coefficients( &
         this%grid%globalDimsR(3), splineOrder, &
         weights(:,1), bsplineFourierCoeffZ, .false.)

    !TODO: Much of this code can be linearized like the C++ counterpart.

    ! Reciprocal space kernel with b-spline discrete Fourier transform
    ! correction.
    !TODO: If b-spline Fourier coefficients are stored in a linearized
    ! array, then this could be simplified to a single loop.
    ixyz = 1
    do iz = 0, this%grid%localDimsR(3) - 1
       do iy = 0, this%grid%localDimsR(2) - 1
          byz = bsplineFourierCoeffY(iy + 1) * bsplineFourierCoeffZ(iz + 1 + this%grid%offsetR(3))
          do ix = 0, gridDimX_k - 1

!#if defined(MPI)
!             igk = 1 + ix + &
!                  iy * (this%grid%localDimsK(1) / 2) + &
!                  iz * this%grid%localDimsK(2) * (this%grid%localDimsK(1) / 2)
!#else
!
!             igk = 1 + ix + iy * (this%grid%localDimsR(1) / 2 ) +&
!                    iz * this%grid%localDimsR(2) * (this%grid%localDimsR(1) / 2 )
!             if (ix .eq. gridDimX_k - 1) then
!                 igk = 1 + iy + iz * this%grid%globalDimsR(2) &
!               + this%grid%globalDimsR(1) / 2 * this%grid%globalDimsR(2) * this%grid%globalDimsR(3)
!             end if
!#endif
             bxyz = bsplineFourierCoeffX(ix + 1) * byz
             waveVector = kxi(ix+1,:) + kyi(iy+1,:) + kzi(iz+1,:)

             k2 = dot_product(waveVector, waveVector)
             !k2 = this%grid%waveVectors2(igk)
!            if ( abs(k2-this%grid%waveVectors2(igk)) .gt. 1e-8   ) then
!               write (6,*) "-|: ", iz,iy,ix,ixyz,igk,k2-this%grid%waveVectors2(igk), this%grid%mpirank
!            end if
             kernel(ixyz) = (4 * pi / k2) * exp(k2 * t) / bxyz
             ixyz = ixyz + 1
          end do
       end do
    end do

    if (this%grid%offsetR(3) == 0) then
       ! Remove the k = 0 term (tinfoil boundary conditions).
       kernel(1) = 0d0
    end if

    ! Allocate the FFT.
    !TODO: Consider switching to in-place, transposed in/out FFT.
    ! Speed vs. memory usage tradeoff.
#ifdef MPI
    planfwd = fftw_mpi_plan_dft_r2c_3d(N, M, L, &
         uuv1d_r, uuv1d_c, this%grid%mpicomm, &
         !TODO: Enable FFTW_MPI_TRANSPOSED_OUT for speedup.
         ! ior(ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT), FFT_ALIGNED))
         ! ior(FFTW_ESTIMATE, FFTW_MPI_TRANSPOSED_OUT))
         FFTW_ESTIMATE)

    planbwd = fftw_mpi_plan_dft_c2r_3d(N, M, L, &
         uuv1d_c, uuv1d_r, this%grid%mpicomm, &
         ! ior(ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN), FFT_ALIGNED))
         FFTW_ESTIMATE)
#else
    planfwd = fftw_plan_dft_r2c_3d( &
         this%grid%globalDimsR(3), &
         this%grid%globalDimsR(2), &
         this%grid%globalDimsR(1), &
         uuv1d_r, uuv1d_c, FFTW_ESTIMATE)
    planbwd = fftw_plan_dft_c2r_3d( &
         this%grid%globalDimsR(3), &
         this%grid%globalDimsR(2), &
         this%grid%globalDimsR(1), &
         uuv1d_c, uuv1d_r, FFTW_ESTIMATE)
#endif

    ! Initialize grids.
    ucu = 0
    uuv1d_r = 0
    uuv1d_c = 0
    this%uuv1d = 0

    ! Spread the Gaussian charges to the grid using b-splines.
    do iu = 1, this%solute%numAtoms
#ifdef MPI
       !TODO: Perform a check earlier for whether the atom
       ! can contribute to this portion of the grid. If not, skip it.
#endif
#if 1
       ! Convert Cartesian position to reciprocal space by projecting
       ! to reciprocal unit cell vectors.
       do id = 1, 3
          reciprocalPos(id) = dot_product(this%solute%position(:,iu), &
               this%grid%unitCellVectorsK(id, :))
       end do

       ! Set the box length to unity since reciprocal space already
       ! divides positions by the box length.
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(1), 1d0, &
            this%grid%globalDimsR(1), &
            splineOrder, gridPoints(:,1), weights(:,1))
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(2), 1d0, &
            this%grid%globalDimsR(2), &
            splineOrder, gridPoints(:,2), weights(:,2))
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(3), 1d0, &
            this%grid%globalDimsR(3), &
            splineOrder, gridPoints(:,3), weights(:,3))
#else
       ! Assume orthorhombic unit cell.
       call cardinal_bspline_periodic_box( &
            this%solute%position(1,iu), this%grid%boxLength(1), &
            this%grid%globalDimsR(1), &
            splineOrder, weights(:,1), gridPoints(:,1))
       call cardinal_bspline_periodic_box( &
            this%solute%position(2,iu), this%grid%boxLength(2), &
            this%grid%globalDimsR(2), &
            splineOrder, weights(:,2), gridPoints(:,2))
       call cardinal_bspline_periodic_box( &
            this%solute%position(3,iu), this%grid%boxLength(3), &
            this%grid%globalDimsR(3), &
            splineOrder, weights(:,3), gridPoints(:,3))
#endif

       do k = 1, splineOrder
#ifdef MPI
          ! Check if grid point is within local grid.
          if (gridPoints(k,3) + 1 < this%grid%offsetR(3) + 1 &
               .or. gridPoints(k,3) + 1 > this%grid%offsetR(3) + this%grid%localDimsR(3)) then
             cycle
          end if
#endif
          do j = 1, splineOrder
             byz = weights(j,2) * weights(k,3)
             do i = 1, splineOrder
                uuv1d_r(1 + gridPoints(i,1), 1 + gridPoints(j,2), &
                     1 + gridPoints(k,3) - this%grid%offsetR(3)) = &
                     uuv1d_r(1 + gridPoints(i,1), 1 + gridPoints(j,2), &
                     1 + gridPoints(k,3) - this%grid%offsetR(3)) &
                     + weights(i,1) * byz * this%solute%charge(iu)
             end do
          end do
       end do
    end do


    ! Convert the charge density into reciprocal space.

#ifdef MPI
    call fftw_mpi_execute_dft_r2c(planfwd, uuv1d_r, uuv1d_c)
#else
    call fftw_execute_dft_r2c(planfwd, uuv1d_r, uuv1d_c)
#endif

    ! Convolute with the kernel.
    do igz = 0, local_N - 1
       do igy = 0, M - 1
          do igx = 0, (L / 2 + 1) - 1
             igk_fort = 1 + igx + (igy + igz * M) * (L / 2 + 1)
             uuv1d_c(igx + 1, igy + 1, igz + 1) &
                  = uuv1d_c(igx + 1, igy + 1, igz + 1) * kernel(igk_fort)
          end do
       end do
    end do

    ! Evaluate the recip-space potential at the grid points.
    ! outr = 0
    uuv1d_r = 0
#ifdef MPI
    call fftw_mpi_execute_dft_c2r(planbwd, uuv1d_c, uuv1d_r)
#else
    call fftw_execute_dft_c2r(planbwd, uuv1d_c, uuv1d_r)
#endif

    chargeCorrection = - pi * (this%solute%totalCharge/this%grid%boxVolume) &
            * zeta
    do igz = 0, local_N - 1
       do igy = 0, M - 1
          do igx = 0, L - 1
             ucu(igx + 1, igy + 1, igz + 1, :) = uuv1d_r(igx + 1, igy + 1, igz + 1) &
                  / this%grid%boxVolume + chargeCorrection
          end do
       end do
    end do

!    call timer_stop(TIME_UCOULULR)

    ! Deallocate the FFT plans.
    call fftw_destroy_plan(planfwd)
    call fftw_destroy_plan(planbwd)

    ! Calculate electrostatic potential energy on each grid point.
    ! No 1/2 term is required since solute affects solvent, but
    ! solvent does not affect solute, hence no double counting.
    do iv = this%solvent%numAtomTypes, 1, -1
          ucu(:,:,:,iv) = ucu(:,:,:,1) * this%solvent%charge(iv)
    end do

    !TODO: Still causes crashes in rare circumstances. Need to debug
    ! more thoroughly.
    if (safemem_dealloc(kxi) /= 0) then
       call rism_report_error("uvPMErecip: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(kyi) /= 0) then
       call rism_report_error("uvPMErecip: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(kzi) /= 0) then
       call rism_report_error("uvPMErecip: Failed to deallocate arrays.")
    end if
    ! if (safemem_dealloc(k2i) /= 0) then
    !    call rism_report_error("uvPMErecip: Failed to deallocate arrays.")
    ! end if
    if (safemem_dealloc(bsplineFourierCoeffX) /= 0) then
       call rism_report_error("uvPMErecip: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(bsplineFourierCoeffY) /= 0) then
       call rism_report_error("uvPMErecip: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(bsplineFourierCoeffZ) /= 0) then
       call rism_report_error("uvPMErecip: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(gaussianFourierCoeff) /= 0) then
       call rism_report_error("uvPMErecip: Failed to deallocate arrays.")
    end if
    call fftw_free(uuv1d_r_cptr)
    call fftw_free(uuv1d_c_cptr)
#ifdef MPI
    if (this%grid%mpirank == 0) then
       call fftw_free(uuv1d_final_cptr)
    end if
#endif
    if (safemem_dealloc(kernel) /= 0) then
       call rism_report_error("uvPMErecip: Failed to deallocate arrays.")
    end if
  end subroutine uvPMErecip

  !> Applying minimum-image convention to find distance from grid
  !! point to nearest solute atom image, which may be in an adjacent
  !! cell. Closest image is based on solute unit cell dimensions.
  !! The closest image is whichever one is less than half a
  !! unit cell away, i.e. minimum image convention.
  !! This is only useful for periodic solute and requires a defined
  !! unit cell.
  !! Note that the minimum image convention implicitly implies a
  !! close-range interaction cutoff at half the unit cell length,
  !! which may or may not be physically realistic depending on the
  !! system.
  !! DAC note: might be worth having a special routine for
  !!   orthogonal boxes
  function minimumImage(this, position)
   implicit none
    _REAL_ :: minimumImage(3)
    type(rism3d_potential), intent(in) :: this !< potential object.
    _REAL_, intent(in) :: position(3) !< Position vector.
    !! The vector origin is usually the grid point a calculation is
    !! performed at.

    integer :: id

    ! Applying minimal-image convention to find distance from grid
    ! point to nearest solute atom image. Closest image is based on
    ! solute unit cell dimensions.

    _REAL_ :: f(3)
    do id = 1, 3
       ! 1. Transform to fractional coordinates (i.e., reciprocal space).
       f(id) = dot_product(position, this%grid%unitCellVectorsK(id, :))
       ! 2. Round to the nearest whole unit and subtract from
       ! coordinate to obtain the minimum image.
       f(id) = f(id) - anint(f(id))
    end do
    ! f = f - anint(f)
    do id = 1, 3
       ! 3. Transform back to Cartesian coordinates using the
       ! transpose of the matrix of Cartesian unit cell vectors.
       minimumImage(id) = dot_product(f, this%grid%unitCellVectorsR(:, id))
    end do
  end function minimumImage

end module rism3d_potential_c

