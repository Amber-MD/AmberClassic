#include "../include/assert.fh"

module xray_dpartial_impl_gpu_module
  
  use xray_contracts_module
  use xray_dpartial_data_module
  use xray_pure_utils, only : real_kind
  use xray_non_bulk_data_module, only : na, nb, nc
  
  implicit none
  private
  
  public :: calc_partial_d_target_d_frac
  public :: calc_partial_d_vls_d_frac
  public :: finalize
  public :: init
  
  interface
    
    subroutine pmemd_xray_dpartial_init_gpu(&
        & n_hkl, &
        & hkl, &
        & mss4, &
        & f_calc, &
        & abs_f_calc, &
        & n_atom, nb, &
        & atom_b_factor, &
        & atom_occupancy, &
        & atom_scatter_type, &
        & n_scatter_types, &
        & atomic_scatter_factor &
        ) bind(C)
      use iso_c_binding
      implicit none
      
      integer(c_int), value :: n_hkl
      integer(c_int), target, intent(in) :: hkl(3, n_hkl)
      real(c_double), target, intent(in) :: mss4(n_hkl)
      complex(c_double_complex), intent(in) :: f_calc(n_hkl)
      real(c_double), target, intent(in) :: abs_f_calc(n_hkl)
      integer(c_int), value :: n_atom
      integer(c_int), value :: nb
      real(c_double), target, intent(in) :: atom_b_factor(n_atom)
      real(c_double), target, intent(in) :: atom_occupancy(n_atom)
      integer(c_int), target, intent(in) :: atom_scatter_type(n_atom)
      integer(c_int), value :: n_scatter_types
      real(c_double), target, intent(in) :: atomic_scatter_factor(n_hkl, n_scatter_types)
    
    end subroutine pmemd_xray_dpartial_init_gpu
    
    subroutine pmemd_xray_dpartial_calc_d_target_d_frac(&
        & n_atom, nb, &
        & frac, &
        & n_hkl, &
        & f_scale, &
        & d_target_d_abs_f_calc, &
        & d_target_d_frac &
        ) bind(C)
      use iso_c_binding
      implicit none
      integer(c_int), value :: n_atom, nb
      real(c_double), intent(in) :: frac(3, n_atom)
      integer(c_int), value :: n_hkl
      real(c_double), intent(in) :: f_scale(n_hkl)
      real(c_double), intent(in) :: d_target_d_abs_f_calc(3, n_hkl)
      real(c_double), intent(inout) :: d_target_d_frac(3, n_atom)
    
    end subroutine pmemd_xray_dpartial_calc_d_target_d_frac
    
    subroutine pmemd_xray_dpartial_finalize_gpu() bind(C)
      implicit none
    end subroutine pmemd_xray_dpartial_finalize_gpu
  end interface


contains

  ! stub for un-written GPU vls implementation:
  function calc_partial_d_vls_d_frac() result(d_vls_d_frac)
  implicit none
  real(real_kind) :: d_vls_d_frac      
  d_vls_d_frac = 0.    ! to avoid a compiler warning
  end function calc_partial_d_vls_d_frac
  
  function calc_partial_d_target_d_frac(frac, f_scale, &
         d_target_d_abs_Fcalc) result(d_target_d_frac)
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    real(real_kind), intent(in) :: f_scale(:)
    real(real_kind), intent(in) :: d_target_d_abs_Fcalc(:)
    real(real_kind) :: d_target_d_frac(3, size(frac, 2))
    
    ASSERT(size(frac, 1) == 3)
    ASSERT(size(frac, 2) == size(atom_b_factor))
    ASSERT(size(frac, 2) == size(atom_scatter_type))
    ASSERT(size(f_scale) == size(hkl, 2))
    ASSERT(size(d_target_d_abs_Fcalc) == size(hkl, 2))
    
    ASSERT(all(abs_Fcalc >= 0))
    ASSERT(all(mSS4 <= 0))
    
    call pmemd_xray_dpartial_calc_d_target_d_frac(&
        & size(frac, 2), nb, &
        & frac, &
        & size(hkl, 2), &
        & f_scale, &
        & d_target_d_abs_Fcalc, &
        & d_target_d_frac &
        &)
  
  end function calc_partial_d_target_d_frac
  
  
  subroutine init(hkl_, mss4_, Fcalc_, abs_Fcalc_, atom_b_factor_,  &
        atom_occupancy_, atom_scatter_type_)
    use xray_dpartial_impl_cpu_module, only : cpu_init => init
    implicit none
    integer, target, intent(in) :: hkl_(:, :)
    real(real_kind), target, intent(in) :: mSS4_(:)
    complex(real_kind), target, intent(in) :: Fcalc_(:)
    real(real_kind), target, intent(in) :: abs_Fcalc_(:)
    real(real_kind), intent(in) :: atom_b_factor_(:)
    real(real_kind), intent(in) :: atom_occupancy_(:)
    integer, intent(in) :: atom_scatter_type_(:)
    
    call cpu_init(hkl_, mss4_, Fcalc_, abs_Fcalc_, atom_b_factor_, &
         atom_occupancy_, atom_scatter_type_)
    call gpu_init()
  
  end subroutine init
  
  subroutine finalize()
    use xray_dpartial_impl_cpu_module, only : cpu_finalize => finalize
    implicit none
    
    call gpu_finalize()
    call cpu_finalize()
  end subroutine finalize
  
  subroutine gpu_init()
    ! use xray_dpartial_data_module
    use xray_atomic_scatter_factor_module, only : atomic_scatter_factor
    implicit none
    
    ASSERT(size(atomic_scatter_factor, 1) == size(hkl, 2))
    
    call pmemd_xray_dpartial_init_gpu(&
        & size(hkl, 2), &
        & hkl, &
        & mSS4, &
        & Fcalc, &
        & abs_Fcalc, &
        & size(atom_b_factor), &
        & nb, &
        & atom_b_factor, &
        & atom_occupancy, &
        & atom_scatter_type, &
        & size(atomic_scatter_factor, 2), &
        & atomic_scatter_factor &
        &)
  end subroutine gpu_init
  
  subroutine gpu_finalize()
    implicit none
    call pmemd_xray_dpartial_finalize_gpu()
  end subroutine gpu_finalize

end module xray_dpartial_impl_gpu_module
