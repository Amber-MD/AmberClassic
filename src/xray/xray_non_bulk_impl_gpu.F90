#include "../include/assert.fh"

module xray_non_bulk_impl_gpu_module
  
  use xray_contracts_module
  use xray_pure_utils, only : real_kind
  use xray_non_bulk_impl_cpu_module, only : get_f_non_bulk
  
  implicit none
  
  private
  
  public :: init
  public :: finalize
  public :: calc_f_non_bulk
  public :: get_f_non_bulk
  
  interface
    subroutine pmemd_xray_non_bulk_init_gpu(n_hkl, hkl, f_non_bulk, mSS4, n_atom, &
        & nb, b_factor, occupancy, &
        & n_scatter_types, scatter_type_index, atomic_scatter_factor) bind(C)
      use iso_c_binding
      implicit none
      integer(c_int), value :: n_hkl
      integer(c_int), intent(in) :: hkl(3, n_hkl)
      complex(c_double_complex), intent(out) :: f_non_bulk(n_hkl)
      real(c_double), intent(in) :: mSS4(n_hkl)
      integer(c_int), value :: n_atom
      integer(c_int), value :: nb
      real(c_double), intent(in) :: b_factor(n_atom)
      real(c_double), intent(in) :: occupancy(n_atom)
      integer(c_int), value :: n_scatter_types
      integer(c_int), intent(in) :: scatter_type_index(n_atom)
      real(c_double), intent(in) :: atomic_scatter_factor(n_hkl, n_scatter_types)
    
    end subroutine pmemd_xray_non_bulk_init_gpu
  end interface
  
  interface
    subroutine pmemd_xray_non_bulk_calc_f_non_bulk_gpu(n_atom, nb, frac_xyz) bind(C)
      use iso_c_binding
      implicit none
      integer(c_int), value :: n_atom
      integer(c_int), value :: nb
      real(c_double), intent(in) :: frac_xyz(3, n_atom)
    end subroutine pmemd_xray_non_bulk_calc_f_non_bulk_gpu
  end interface
  
  interface
    subroutine pmemd_xray_non_bulk_finalize_gpu() bind(C)
    end subroutine pmemd_xray_non_bulk_finalize_gpu
  end interface

contains
  
  subroutine init(hkl_, mSS4_, b_factor_, scatter_type_index_, occupancy_)
    use xray_non_bulk_impl_cpu_module, only : cpu_init => init
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
    
    call cpu_init(hkl_, mSS4_, b_factor_, scatter_type_index_, occupancy_)
    call gpu_init()
  end subroutine init
  
  
  subroutine finalize()
    use xray_non_bulk_impl_cpu_module, only : cpu_finalize => finalize
    implicit none
    call gpu_finalize()
    call cpu_finalize()
  end subroutine finalize
  
  subroutine calc_f_non_bulk(frac)
    use xray_non_bulk_data_module
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    
    ASSERT(size(frac, 1) == 3)
    call pmemd_xray_non_bulk_calc_f_non_bulk_gpu(size(frac, 2), nb, frac)
  end subroutine calc_f_non_bulk
  
  
  subroutine gpu_init()
    use xray_non_bulk_data_module
    use xray_atomic_scatter_factor_module, only : atomic_scatter_factor
    implicit none

    ASSERT(allocated(atomic_scatter_factor))
    
    call pmemd_xray_non_bulk_init_gpu( &
        & size(hkl, 2), hkl, f_non_bulk, mSS4, size(b_factor), nb, &
        & b_factor, &
        & occupancy, &
        & size(atomic_scatter_factor, 2), scatter_type_index, atomic_scatter_factor &
    )
  end subroutine gpu_init
  
  subroutine gpu_finalize()
    implicit none
    call pmemd_xray_non_bulk_finalize_gpu()
  end subroutine gpu_finalize

end module xray_non_bulk_impl_gpu_module
