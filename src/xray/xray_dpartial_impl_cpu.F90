#include "../include/assert.fh"

module xray_dpartial_impl_cpu_module
  
  use xray_contracts_module
  use xray_dpartial_data_module
  use xray_pure_utils, only : real_kind
  use constants_xray, only : xray_num_threads
  use xray_interface2_data_module, only : spacegroup_number
  use xray_non_bulk_data_module, only : na, nb, nc
  
  implicit none
  private
  
  public :: calc_partial_d_target_d_frac
  public :: calc_partial_d_vls_d_frac
  public :: finalize
  public :: init

contains
  
  function calc_partial_d_target_d_frac(frac, f_scale, &
        d_target_d_abs_Fcalc) result(d_target_d_frac)
    use xray_atomic_scatter_factor_module, only : atomic_scatter_factor
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    real(real_kind), intent(in) :: f_scale(:)
    real(real_kind), intent(in) :: d_target_d_abs_Fcalc(:)
    real(real_kind) :: d_target_d_frac(3, size(frac, 2))

    ! locals:
    real(real_kind) :: hkl_v(3)
    real(real_kind) :: phase, fa, delta_phi
    complex(real_kind) :: f, fhkl
    integer :: i, ihkl, hkls(3)
    real(real_kind) :: time0, time1

    real(real_kind), parameter :: M_TWOPI = 2 * 3.1415926535897932384626433832795d0 ! FIXME: use constants module
    call wallclock( time0 )
    
    ASSERT(size(frac, 1) == 3)
    ASSERT(size(frac, 2) == size(atom_b_factor))
    ASSERT(size(frac, 2) == size(atom_scatter_type))
    ASSERT(size(f_scale) == size(hkl, 2))
    ASSERT(size(d_target_d_abs_Fcalc) == size(hkl, 2))
    
    ASSERT(all(abs_Fcalc >= 0))
    ASSERT(all(mSS4 <= 0))
    
    d_target_d_frac = 0

!$omp parallel do private(i,ihkl,hkl_v,hkls,phase,f,fa) num_threads(xray_num_threads)
    do i = 1, size(frac, 2)
      do ihkl = 1, size(hkl, 2)
        
        if (abs_Fcalc(ihkl) < 1e-3) then
          ! Note: when Fcalc is approximately zero the phase is undefined,
          ! so no force can be determined even if the energy is high. (Similar
          ! to a linear bond angle.)
          cycle
        end if
        
        hkl_v(:) = hkl(:, ihkl) * M_TWOPI
        
        !  fa should be the same for all symmetry mates
        fa = atomic_scatter_factor(ihkl, atom_scatter_type(i)) &
            * exp(mSS4(ihkl) * atom_b_factor(i)) &
            * f_scale(ihkl) * d_target_d_abs_Fcalc(ihkl) / abs_Fcalc(ihkl)
        
        ! original hkl for symmetrization:
        phase = hkl_v(1)*frac(1,i) + hkl_v(2)*frac(2,i) + hkl_v(3)*frac(3,i)
        f = fa * cmplx(cos(phase), sin(phase), real_kind)
        d_target_d_frac(:, i) = d_target_d_frac(:, i) + hkl_v(:) &
            * (real(f)*aimag(Fcalc(ihkl)) - aimag(f)*real(Fcalc(ihkl))) 

        if( spacegroup_number .eq. 19 ) then

           ! set #2:  -h,-k,l
           hkls(1) = -hkl(1,ihkl)
           hkls(2) = -hkl(2,ihkl)
           hkls(3) =  hkl(3,ihkl)
           hkl_v(:) = hkls(:) * M_TWOPI

           phase = hkl_v(1)*frac(1,i) + hkl_v(2)*frac(2,i) + hkl_v(3)*frac(3,i)
           f = fa * cmplx(cos(phase), sin(phase), real_kind)
           if( mod(hkls(1)/na + hkls(3)/nc, 2) .ne. 0 ) f = -f
           d_target_d_frac(:, i) = d_target_d_frac(:, i) + hkl_v(:) &
               * (real(f)*aimag(Fcalc(ihkl)) - aimag(f)*real(Fcalc(ihkl))) 

           ! set #3:  -h,k,-l
           hkls(1) = -hkl(1,ihkl)
           hkls(2) =  hkl(2,ihkl)
           hkls(3) = -hkl(3,ihkl)
           hkl_v(:) = hkls(:) * M_TWOPI
           phase = hkl_v(1)*frac(1,i) + hkl_v(2)*frac(2,i) + hkl_v(3)*frac(3,i)
           f = fa * cmplx(cos(phase), sin(phase), real_kind)
           if( mod(hkls(2)/nb + hkls(3)/nc, 2) .ne. 0 ) f = -f
           d_target_d_frac(:, i) = d_target_d_frac(:, i) + hkl_v(:) &
               * (real(f)*aimag(Fcalc(ihkl)) - aimag(f)*real(Fcalc(ihkl))) 

           ! set #4:   h,-k,-l
           hkls(1) =  hkl(1,ihkl)
           hkls(2) = -hkl(2,ihkl)
           hkls(3) = -hkl(3,ihkl)
           hkl_v(:) = hkls(:) * M_TWOPI
           phase = hkl_v(1)*frac(1,i) + hkl_v(2)*frac(2,i) + hkl_v(3)*frac(3,i)
           f = fa * cmplx(cos(phase), sin(phase), real_kind)
           if( mod(hkls(1)/na + hkls(2)/nb, 2) .ne. 0 ) f = -f
           d_target_d_frac(:, i) = d_target_d_frac(:, i) + hkl_v(:) &
               * (real(f)*aimag(Fcalc(ihkl)) - aimag(f)*real(Fcalc(ihkl))) 

        else if( spacegroup_number .eq. 4 ) then
           ! set #2:   -h,k,-l
           hkls(1) = -hkl(1,ihkl)
           hkls(2) =  hkl(2,ihkl)
           hkls(3) = -hkl(3,ihkl)
           hkl_v(:) = hkls(:) * M_TWOPI
           phase = hkl_v(1)*frac(1,i) + hkl_v(2)*frac(2,i) + hkl_v(3)*frac(3,i)
           f = fa * cmplx(cos(phase), sin(phase), real_kind)
           if( mod(hkls(2)/nb, 2) .ne. 0 ) f = -f
           d_target_d_frac(:, i) = d_target_d_frac(:, i) + hkl_v(:) &
               * (real(f)*aimag(Fcalc(ihkl)) - aimag(f)*real(Fcalc(ihkl))) 

        end if

      end do
    end do
!$omp end parallel do

    call wallclock( time1 )
    ! write(0,'(a,f8.3)') 'dhkl time: ', time1 - time0
  
  end function calc_partial_d_target_d_frac
  
  function calc_partial_d_vls_d_frac(frac, f_scale) &
         result(d_target_d_frac)
    use xray_atomic_scatter_factor_module, only : atomic_scatter_factor
    use xray_target_vector_least_squares_data_module, only: derivc
    use xray_interface2_data_module, only : n_work
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    real(real_kind), intent(in) :: f_scale(:)
    real(real_kind) :: d_target_d_frac(3, size(frac, 2))
    real(real_kind) :: hkl_v(3)
    real(real_kind) :: f, phase
    integer :: ihkl, i
    
    real(real_kind), parameter :: M_TWOPI = 2 * 3.1415926535897932384626433832795d0 ! FIXME: use constants module

    ASSERT(size(frac, 1) == 3)
    ASSERT(size(frac, 2) == size(atom_b_factor))
    ASSERT(size(frac, 2) == size(atom_scatter_type))
    ASSERT(size(f_scale) == size(hkl, 2))
    
    ASSERT(all(abs_Fcalc >= 0))
    ASSERT(all(mSS4 <= 0))
    
    d_target_d_frac = 0

!$omp parallel do private(i,ihkl,hkl_v,phase,f) num_threads(xray_num_threads)
    do i = 1, size(frac, 2)
      do ihkl = 1, n_work

        if (abs_Fcalc(ihkl) < 1e-3) then
          ! Note: when Fcalc is approximately zero the phase is undefined,
          ! so no force can be determined even if the energy is high. (Similar
          ! to a linear bond angle.)
          cycle
        end if
        
        hkl_v = hkl(:, ihkl) * M_TWOPI
        
        phase = -sum(hkl_v * frac(:, i))
        f = atomic_scatter_factor(ihkl, atom_scatter_type(i)) &
            * exp(mSS4(ihkl) * atom_b_factor(i)) &
            * ( sin(phase) * real(derivc(ihkl)) + cos(phase) * aimag(derivc(ihkl)) )
        
        d_target_d_frac(:, i) = d_target_d_frac(:, i) &
            + atom_occupancy(i) * f_scale(ihkl) * hkl_v(:) * f

      end do
    end do
!$omp end parallel do
  
  end function calc_partial_d_vls_d_frac
  
  
  subroutine init(hkl_, mss4_, Fcalc_, abs_Fcalc_, atom_b_factor_,  &
        atom_occupancy_, atom_scatter_type_)
    implicit none
    integer, target, intent(in) :: hkl_(:, :)
    real(real_kind), target, intent(in) :: mSS4_(:)
    complex(real_kind), target, intent(in) :: Fcalc_(:)
    real(real_kind), target, intent(in) :: abs_Fcalc_(:)
    real(real_kind), intent(in) :: atom_b_factor_(:)
    real(real_kind), intent(in) :: atom_occupancy_(:)
    integer, intent(in) :: atom_scatter_type_(:)
    
    ASSERT(size(hkl_, 1) == 3)
    ASSERT(size(mSS4_) == size(hkl_, 2))
    ASSERT(size(abs_Fcalc_) == size(hkl_, 2))
    ASSERT(size(Fcalc_) == size(hkl_, 2))
    
    ASSERT(size(atom_scatter_type_) == size(atom_b_factor_))
    
    hkl => hkl_
    mSS4 => mss4_
    Fcalc => Fcalc_
    abs_Fcalc => abs_Fcalc_
    atom_b_factor = atom_b_factor_
    atom_occupancy = atom_occupancy_
    atom_scatter_type = atom_scatter_type_
  
  end subroutine init
  
  subroutine finalize()
    hkl => null()
    mSS4 => null()
    Fcalc => null()
    abs_Fcalc => null()
    !  dac: These are also deallocated in xray_interface2_data.F90
    !  deallocate(atom_b_factor)
    ! deallocate(atom_scatter_type)
  end subroutine finalize

end module xray_dpartial_impl_cpu_module
