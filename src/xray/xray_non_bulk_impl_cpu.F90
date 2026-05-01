#include "../include/assert.fh"

module xray_non_bulk_impl_cpu_module
  
  use xray_pure_utils, only : real_kind
  use xray_contracts_module
  use xray_non_bulk_data_module
  use xray_interface2_data_module, only : spacegroup_number
  
  implicit none
  
  private
  
  public :: init
  public :: finalize
  public :: calc_f_non_bulk
  public :: get_f_non_bulk

  real(real_kind), parameter :: M_TWOPI = 2 * 3.1415926535897932384626433832795d0 ! FIXME: use constants module

contains
  
  subroutine init(hkl_, mSS4_, b_factor_, scatter_type_index_, occupancy_)
    implicit none
    integer, intent(in), target :: hkl_(:, :)
    real(real_kind), intent(in), target :: mSS4_(:)
    real(real_kind), intent(in) :: b_factor_(:)
    integer, intent(in) :: scatter_type_index_(:)
    real(real_kind), intent(in) :: occupancy_(:)
    
    ASSERT(size(hkl_, 1) == 3)
    ASSERT(size(hkl_, 2) == size(mSS4_))
    ASSERT(size(b_factor_) == size(scatter_type_index_))
    ASSERT(size(b_factor_) == size(occupancy_))
    
    hkl => hkl_
    mSS4 => mSS4_
    b_factor = b_factor_
    scatter_type_index = scatter_type_index_
    
    occupancy = occupancy_
    
    allocate(F_non_bulk(size(mSS4_)))
    allocate(f(size(b_factor_)))
    allocate(angle(size(b_factor_)))
  
  end subroutine init
  
  subroutine finalize()
    implicit none
    hkl => null()
    mSS4 => null()
    xyz => null()
    if(allocated(b_factor)) deallocate(b_factor)
    if(allocated(scatter_type_index)) deallocate(scatter_type_index)
    if (allocated(occupancy)) deallocate(occupancy)
    if (allocated(F_non_bulk)) deallocate(F_non_bulk)
    if (allocated(f)) deallocate(f)
    if (allocated(angle)) deallocate(angle)
  end subroutine finalize
  
  !-------------------------------------------------------------------
  ! Caution: Some literature uses S to represent S^2
  !
  ! (for d*, see p. 93 of Glusker, Lewis, Rossi)
  ! S == d* = sqrt(sum(HKL * orth_to_frac)^2) = sqrt(-4*mSS4)
  ! mSS4 = -S*S/4
  ! mSS4 is more relevant to the formulas used than S.
  
  function get_f_non_bulk() result(result)
    complex(real_kind) :: result(size(F_non_bulk))
    result(:) = F_non_bulk(:)
  end function get_f_non_bulk
  
  subroutine calc_f_non_bulk(frac)
    use xray_atomic_scatter_factor_module, only : atomic_scatter_factor
    use constants_xray, only: xray_num_threads
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    ! locals
    integer :: ihkl, hkls(3)
    complex(real_kind) :: fcalcs, fhkl
    real(real_kind) :: delta_phi
    double precision :: time0, time1

    call wallclock( time0 )

    ASSERT(associated(hkl))
    ASSERT(associated(mSS4))
    ASSERT(allocated(atomic_scatter_factor))
    ASSERT(allocated(b_factor))
    ASSERT(allocated(scatter_type_index))
    
    ASSERT(size(frac, 1) == 3)
    ASSERT(size(frac, 2) == size(b_factor))
    ASSERT(size(frac, 2) == size(scatter_type_index))
    ASSERT(size(hkl, 2) == size(atomic_scatter_factor, 1))

    !$omp parallel do private(ihkl,f,angle,hkls,fcalcs)  num_threads(xray_num_threads)
    do ihkl = 1, size(hkl, 2)
      
      ! Fhkl = SUM( fj * exp(M_TWOPI * i * (h * xj + k * yj + l * zj)) ),
      !      j = 1,num_selected_atoms
      ! where:
      !    The sum is versus j, over all selected atoms
      !    fj is the atomic scatter for atom j:   atomic_scatter_factor(j)
      !    h,k,l are from the hkl list for index ihkl:   hkl(1:3,ihkl)
      !    x,y,z are coordinates for atom j:   xyz(1:3,j)
      !        xyz(:) may be a reduced list.
      !
      ! Rather than using a complex exponential where the real part is
      ! always zero, this is optimized to calculate sin and cosine parts,
      ! then convert to a complex number
      ! after the A and B components are summed over all selected atoms.
      ! This can be written as:
      !
      ! Ahkl = SUM( fj * cos(M_TWOPI * (h * xj + k * yj + l * zj)) ),
      ! Bhkl = SUM( fj * sin(M_TWOPI * (h * xj + k * yj + l * zj)) ),
      !    j = 1,num_selected_atoms
      
!     f(:) = exp(mSS4(ihkl) * b_factor(:)) * occupancy(:)
      f(:) = exp(mSS4(ihkl) * b_factor(:)) &
          * atomic_scatter_factor(ihkl, scatter_type_index(:))

      ! original hkl for space group symmetrization:
      angle(:) = M_TWOPI*(hkl(1, ihkl)*frac(1,:) + &
                          hkl(2, ihkl)*frac(2,:) + &
                          hkl(3, ihkl)*frac(3,:))
      F_non_bulk(ihkl) = cmplx(sum(f(:) * cos(angle(:))), &
          sum(f(:) * sin(angle(:))), real_kind)

      if( spacegroup_number .eq. 19 ) then

         ! set #2:  -h,-k,l
         hkls(1) = -hkl(1,ihkl)
         hkls(2) = -hkl(2,ihkl)
         hkls(3) =  hkl(3,ihkl)
         angle(:) = M_TWOPI*(hkls(1)*frac(1,:) + &
                             hkls(2)*frac(2,:) + &
                             hkls(3)*frac(3,:))
         fcalcs = cmplx(sum(f(:) * cos(angle(:))), &
             sum(f(:) * sin(angle(:))), real_kind)
         if( mod(hkls(1)/na + hkls(3)/nc, 2) .ne. 0 ) fcalcs = -fcalcs
         F_non_bulk(ihkl) = F_non_bulk(ihkl) + fcalcs

         ! set #3:  -h,k,-l
         hkls(1) = -hkl(1,ihkl)
         hkls(2) =  hkl(2,ihkl)
         hkls(3) = -hkl(3,ihkl)
         angle(:) = M_TWOPI*(hkls(1)*frac(1,:) + &
                             hkls(2)*frac(2,:) + &
                             hkls(3)*frac(3,:))
         fcalcs = cmplx(sum(f(:) * cos(angle(:))), &
             sum(f(:) * sin(angle(:))), real_kind)
         if( mod(hkls(2)/nb + hkls(3)/nc, 2) .ne. 0 ) fcalcs = -fcalcs
         F_non_bulk(ihkl) = F_non_bulk(ihkl) + fcalcs
    
         ! set #4:  h,-k,-l
         hkls(1) =  hkl(1,ihkl)
         hkls(2) = -hkl(2,ihkl)
         hkls(3) = -hkl(3,ihkl)
         angle(:) = M_TWOPI*(hkls(1)*frac(1,:) + &
                             hkls(2)*frac(2,:) + &
                             hkls(3)*frac(3,:))
         fcalcs = cmplx(sum(f(:) * cos(angle(:))), &
             sum(f(:) * sin(angle(:))), real_kind)
         if( mod(hkls(1)/na + hkls(2)/nb, 2) .ne. 0 ) fcalcs = -fcalcs
         F_non_bulk(ihkl) = F_non_bulk(ihkl) + fcalcs

      else if ( spacegroup_number .eq. 4 ) then

         ! set #2:  -h,k,-l
         hkls(1) = -hkl(1,ihkl)
         hkls(2) =  hkl(2,ihkl)
         hkls(3) = -hkl(3,ihkl)
         angle(:) = M_TWOPI*(hkls(1)*frac(1,:) + &
                             hkls(2)*frac(2,:) + &
                             hkls(3)*frac(3,:))
         fcalcs = cmplx(sum(f(:) * cos(angle(:))), &
             sum(f(:) * sin(angle(:))), real_kind)
         if( mod(hkls(2)/nb, 2) .ne. 0 ) fcalcs = -fcalcs
         F_non_bulk(ihkl) = F_non_bulk(ihkl) + fcalcs

      end if

    end do
    !$omp end parallel do

    call wallclock( time1 )
    ! write(0,'(a,f8.3)') 'ihkl time: ', time1 - time0
  
  end subroutine calc_f_non_bulk


end module xray_non_bulk_impl_cpu_module
