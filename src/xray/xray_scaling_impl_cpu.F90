#include "../include/assert.fh"

module xray_scaling_impl_cpu_module
  
  use xray_pure_utils, only : real_kind, is_sorted
  use xray_contracts_module
  use xray_scaling_data_module

  implicit none

  private

  public :: init
  public :: optimize_scale_factors
  public :: combine
  public :: get_f_scale
  public :: rescale
  public :: finalize


  public :: update_k_iso_exp                  ! Made public only for unit tests
  public :: update_k_aniso                    ! Made public only for unit tests
  public :: update_k_bulk_k_iso               ! Made public only for unit tests
  public :: update_k_bulk_k_iso_via_cubic_eq  ! Made public only for unit tests

contains
  
  ! - - - - - - - - - - - - - - - !
  ! Public functions/subroutines  !
  ! - - - - - - - - - - - - - - - !
  
  subroutine init(resolution, reciprocal_norms, num_work_flags, hkl, max_resolution_bins, n_reflections_in_worst_resolution_bin, min_bin_size)
   
    implicit none
    real(real_kind), intent(in) :: resolution(:)
    real(real_kind), intent(in) :: reciprocal_norms(3)
    integer, intent(in) :: num_work_flags
    integer, intent(in) :: hkl(3, size(resolution))
  
    ! Default values set to match defaults of set.log_binning in cctbx
    integer, intent(in), optional :: max_resolution_bins
    integer, intent(in), optional :: n_reflections_in_worst_resolution_bin
    integer, intent(in), optional :: min_bin_size

    ! Actual values of parameters
    integer :: max_resolution_bins_value
    integer :: n_reflections_in_worst_resolution_bin_value
    integer :: min_bin_size_value

    if (present(max_resolution_bins)) then
      max_resolution_bins_value = max_resolution_bins
    else
      max_resolution_bins_value = 30
    end if
  
    if (present(n_reflections_in_worst_resolution_bin)) then
      n_reflections_in_worst_resolution_bin_value = n_reflections_in_worst_resolution_bin
    else
      n_reflections_in_worst_resolution_bin_value = 100
    end if
  
    if (present(min_bin_size)) then
      min_bin_size_value = min_bin_size
    else
      min_bin_size_value = 50
    end if
    
    call init_impl(resolution, reciprocal_norms, num_work_flags, hkl, max_resolution_bins_value, n_reflections_in_worst_resolution_bin_value, min_bin_size_value)
    
  end subroutine init
  
  subroutine optimize_scale_factors(absFobs, Fprot, Fbulk, neg_s_norm2_div4, hkl)
    use xray_pure_utils, only : calc_r_factor
    implicit none
    
    real(real_kind), intent(in) :: absFobs(:) ! Magnitude of experimental structure factors
    complex(real_kind), intent(in) :: Fprot(size(absFobs)) ! Unscaled protein (non-bulk) structure factors
    complex(real_kind), intent(in) :: Fbulk(size(absFobs)) ! Unscaled bulk structure factors
    real(real_kind), intent(in) :: neg_s_norm2_div4(size(absFobs))
    integer, intent(in) :: hkl(3, size(absFobs))
    
    real(real_kind) :: prev_r_work, r_work
    real(real_kind) :: k_bulk_bin(n_resolution_bins), &
        k_iso_bin(n_resolution_bins)
    logical :: updated
    integer :: i
    
    ASSERT(size(Fprot) == size(absFobs))
    ASSERT(size(Fbulk) == size(absFobs))
    ASSERT(size(neg_s_norm2_div4) == size(absFobs))
    ASSERT(size(hkl, 2) == size(absFobs))
    ASSERT(size(k_iso) == size(absFobs))
    ASSERT(size(k_iso_exp) == size(absFobs))
    ASSERT(size(k_aniso) == size(absFobs))
    ASSERT(size(k_bulk) == size(absFobs))
    ASSERT(all(k_iso >= 0))
    ASSERT(all(k_iso_exp >= 0))
    ASSERT(all(k_aniso >= 0))
    !ASSERT(all(k_bulk >= 0))
    !ASSERT(all(k_bulk <= 1))

    k_bulk = 0
    k_iso = 1
    k_iso_exp = 1
    k_aniso = 1
    
    r_work = calc_r_factor(absFobs(:n_work), abs(k_iso(:n_work) * k_iso_exp(:n_work) * k_aniso(:n_work) * (Fprot(:n_work) + k_bulk(:n_work) * Fbulk(:n_work))))

    do i = 1, 2
      updated = update_k_iso_exp(absFobs, k_iso, abs(Fprot + k_bulk * Fbulk), neg_s_norm2_div4, r_work)
      updated = update_k_bulk_k_iso(absFobs, Fprot, Fbulk, neg_s_norm2_div4, k_bulk_bin, k_iso_bin, r_work)
      updated = update_k_iso_exp(absFobs, k_iso, abs(Fprot + k_bulk * Fbulk), neg_s_norm2_div4, r_work)
    end do
    updated = update_k_aniso(absFobs, k_iso * k_iso_exp * abs(Fprot + k_bulk * Fbulk), hkl(1,:), hkl(2,:), hkl(3,:), r_work)
    
    prev_r_work = r_work
    
    do i = 2, 20 ! Skip one cycle to keep it closer to initial sf.f90 code
      updated = update_k_bulk_k_iso_via_cubic_eq(absFobs, Fprot, Fbulk, neg_s_norm2_div4, k_bulk_bin, k_iso_bin, r_work)
      updated = update_k_aniso(absFobs, k_iso * k_iso_exp * abs(Fprot + k_bulk * Fbulk), hkl(1,:), hkl(2,:), hkl(3,:), r_work)

      if (r_work > prev_r_work - 1e-4) then
        exit
      end if
      
      prev_r_work = r_work
    end do
    
    ASSERT(all(k_iso >= 0))
    ASSERT(all(k_iso_exp >= 0))
    ASSERT(all(k_aniso >= 0))
    ! ASSERT(all(k_bulk >= 0))
    ! ASSERT(all(k_bulk <= 1))
    
  end subroutine optimize_scale_factors
  
  function combine (Fprot, Fbulk) result(result)
    use xray_interface2_data_module, only : has_Fuser
    complex(real_kind), intent(in) :: Fprot(:) ! Unscaled protein (non-bulk) structure factors
    complex(real_kind), intent(in) :: Fbulk(size(Fprot)) ! Unscaled bulk structure factors
    complex(real_kind) :: result(size(Fprot)) ! Final Fcalc structure factors
    
    if( has_Fuser > 0 ) then
       result = Fprot + Fbulk  ! don't scale Fuser relative to Fprot
       write(6,*) 'combining Fprot and Fbulk'
    else
       result = Fprot + k_bulk * Fbulk
    endif
  
  end function combine
  
  function rescale(Fcalc) result(result)
    complex(real_kind), intent(in)  :: Fcalc(:) ! Unscaled Fcalc (non-bulk) structure factors
    complex(real_kind) :: result(size(Fcalc))
    
    result = (k_iso * k_iso_exp * k_aniso) * Fcalc
  end function rescale
  
  subroutine finalize()
  
    if (allocated(hkl_scale_resolution_bin)) deallocate(hkl_scale_resolution_bin)
    if (allocated(work_scale_resolution_bin_start)) deallocate(work_scale_resolution_bin_start)
    if (allocated(work_scale_resolution_bin_size)) deallocate(work_scale_resolution_bin_size)
    if (allocated(free_scale_resolution_bin_start)) deallocate(free_scale_resolution_bin_start)
    if (allocated(free_scale_resolution_bin_size)) deallocate(free_scale_resolution_bin_size)
 
    if (allocated(k_bulk)) deallocate(k_bulk)
    if (allocated(k_iso)) deallocate(k_iso)
    if (allocated(k_iso_exp)) deallocate(k_iso_exp)
    if (allocated(k_aniso)) deallocate(k_aniso)
    if (allocated(MVcryst_base)) deallocate(MVcryst_base)
    
  end subroutine finalize
  
  ! - - - - - - - - - - - - - - - !
  ! Private functions/subroutines !
  ! - - - - - - - - - - - - - - - !
  
  function calc_k_aniso(Uaniso, h, k, l) result(result)
    integer, intent(in) :: h(:), k(:), l(:)  !< Miller indices
    real(real_kind), intent(in) :: Uaniso(6)
    real(real_kind) :: result(size(h))

    result = max(0.0_real_kind, exp( &
        Uaniso(1) * h**2 + Uaniso(2) * k**2 + Uaniso(3) * l**2 + &
            Uaniso(4) * h * k * 2 + Uaniso(5) * h * l * 2 + Uaniso(6) * k * l * 2))

  end function calc_k_aniso

  function calc_k_aniso_from_V(Vaniso, V_base) result(result)
    real(real_kind), intent(in) :: V_base(:, :)  !< MVcryst_base array
    real(real_kind), intent(in) :: Vaniso(12)
    real(real_kind) :: result(size(V_base(1, :)))
    integer :: i

    result = 1.0_real_kind
    do i = 1, 12
      result = result + Vaniso(i) * V_base(i, :)
    end do

  end function calc_k_aniso_from_V
  
  function calc_k_iso_bin(absFobs, absFcalc) result(result)
    use xray_pure_utils, only: calc_k_overall
    implicit none
    real(real_kind), intent(in) :: absFobs(:)
    real(real_kind), intent(in) :: absFcalc(size(absFobs))
    real(real_kind) :: result(n_resolution_bins)
    
    integer :: i, s, e
    real(real_kind) :: bin_k_iso
    
    do i = 1, n_resolution_bins
      s = work_scale_resolution_bin_start(i)
      e = s + work_scale_resolution_bin_size(i) - 1
      
      result(i) = calc_k_overall(absFobs(s:e), absFcalc(s:e))
    end do
  
  end function calc_k_iso_bin
  
  pure function calc_k_iso_exp(neg_s_norm2_div4, factor, decay) result(result)
    real(real_kind), intent(in) :: neg_s_norm2_div4(:)
    
    real(real_kind), intent(in) :: factor, decay
    real(real_kind) :: result(size(neg_s_norm2_div4))
    
    result = factor * exp(neg_s_norm2_div4 * decay)
  end function
  
  function calc_Uaniso(absFobs, absFcalc, h, k, l) result(Uaniso)
    implicit none
    real(real_kind), intent(in) :: absFobs(:)
    real(real_kind), intent(in) :: absFcalc(size(absFobs))
    integer, dimension(size(absFobs)), intent(in) :: h, k, l  !< Miller indices
    real(real_kind) :: Uaniso(6)
    real(real_kind) :: b(6)
    real(real_kind) :: b_vector_base(size(absFobs))
    
    ASSERT(size(h) == size(k))
    ASSERT(size(h) == size(l))
    ASSERT(size(h) == size(absFobs))
    ASSERT(size(h) == size(absFcalc))
    
    b_vector_base = log(absFobs / absFcalc) / (size(absFobs) ** 2)
    
    b = (/&
            sum(b_vector_base * h ** 2), &
            sum(b_vector_base * k ** 2), &
            sum(b_vector_base * l ** 2), &
            sum(b_vector_base * h * k) * 2, &
            sum(b_vector_base * h * l) * 2, &
            sum(b_vector_base * k * l) * 2 &
        /)
    
    Uaniso = matmul(MUcryst_inv, b) ! in Phenix it is u_star  multiplied by (-2 *pi *pi)
  
  end function calc_Uaniso

  function calc_Vaniso(absFobs, absFcalc, V_base) result(Vaniso)
    implicit none
    real(real_kind), intent(in) :: absFobs(:)
    real(real_kind), intent(in) :: absFcalc(size(absFobs))
    real(real_kind), intent(in) :: V_base(12, size(absFobs))
    real(real_kind) :: Vaniso(12)
    real(real_kind) :: b(12)
    real(real_kind) :: MVcryst_inv(12, 12)
    real(real_kind) :: b_vector_base(size(absFobs))
    integer :: i

    b_vector_base = (absFobs - absFcalc) * absFcalc / (size(absFobs) ** 2)

    b = 0.0_real_kind
    do i = 1, 12
      b(i) = sum(b_vector_base * V_base(i, :))
    end do

    MVcryst_inv = init_MVcryst_inv(absFcalc, V_base)

    Vaniso = matmul(MVcryst_inv, b)

  end function calc_Vaniso
  
  function count_high_resolution_reflexes(resolution) result(result)
    
    implicit none
    
    real(real_kind), intent(in) :: resolution(:)
    real(real_kind) :: max_resolution
    integer :: result
    integer :: i

    ! Preconditions:
    ASSERT(size(resolution) > 1)
    ASSERT(is_sorted(resolution))
    
    result = 0
    max_resolution = min(4.0_real_kind, resolution(1) + 1.0_real_kind) ! Best of 4A and highest resolution + 1A
    
    do i = 1, size(resolution)
      if (resolution(i) < max_resolution) then
        result = i
      else
        return
      end if
    end do
  
  end function count_high_resolution_reflexes

  subroutine grid_search_k_iso_exp(absFobs, absFcalc, neg_s_norm2_div4, b_lower_limit, factor_kb, decay_kb)
    use xray_pure_utils, only: calc_unscaled_r_factor, calc_k_overall

    real(real_kind), intent(in) :: absFobs(:)
    real(real_kind), intent(in) :: absFcalc(size(absFobs))
    real(real_kind), intent(in) :: neg_s_norm2_div4(size(absFobs))
    real(real_kind), intent(in) :: b_lower_limit

    real(real_kind), intent(out) :: factor_kb
    real(real_kind), intent(out) :: decay_kb

    real(real_kind) :: factor_test, decay_test
    real(real_kind) :: r_kb, r_test
    integer :: i

    factor_kb = 1
    decay_kb = 0

    r_kb = 1.0e10
    do i = 0, 200
      decay_test = b_lower_limit + i * 1.0
      factor_test = calc_k_overall(absFobs, exp(decay_test * neg_s_norm2_div4) * absFcalc)
      r_test = calc_unscaled_r_factor(absFobs, calc_k_iso_exp(neg_s_norm2_div4, factor_test, decay_test) * absFcalc)
      if (r_test < r_kb) then
        factor_kb = factor_test
        decay_kb = decay_test
        r_kb = r_test
      end if
    end do
  end subroutine grid_search_k_iso_exp

  subroutine grid_search_k_bulk(absFobs, Fprot, Fbulk, k_bulk_grid, k_bulk, r, k_overall, k_upper_limit)
    use xray_pure_utils, only: calc_r_factor, calc_unscaled_r_factor, calc_k_overall
    implicit none
    
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(in) :: Fprot(size(absFobs))
    complex(real_kind), intent(in) :: Fbulk(size(absFobs))
    real(real_kind), intent(in) :: k_bulk_grid(:)
    real(real_kind), intent(inout) :: k_bulk
    real(real_kind), intent(inout) :: r
    real(real_kind), optional, intent(inout) :: k_overall
    real(real_kind), optional, intent(in) :: k_upper_limit
    
    real(real_kind) :: k, tmp_r, tmp_k_overall
    real(real_kind) :: tmp_k_bulk_grid(size(k_bulk_grid) + 1)
    logical :: valid_k_condition
    integer :: j
    
    tmp_k_overall = 1
    tmp_k_bulk_grid(1) = 0.0_real_kind
    tmp_k_bulk_grid(2:size(k_bulk_grid) + 1) = k_bulk_grid
    
    do j = 1, size(tmp_k_bulk_grid)
      k = tmp_k_bulk_grid(j)
      if (present(k_upper_limit)) then
        valid_k_condition = (0.0_real_kind <= k .and. k <= k_upper_limit)
      else
        valid_k_condition = 0.0_real_kind <= k
      end if
      if (valid_k_condition) then
        if (present(k_overall)) then
          tmp_k_overall = calc_k_overall(absFobs, abs(Fprot + k * Fbulk))
          tmp_r = calc_unscaled_r_factor(absFobs, abs(Fprot + k * Fbulk) * tmp_k_overall)
        else
          tmp_r = calc_r_factor(absFobs, abs(Fprot + k * Fbulk))
        end if
        if (tmp_r < r) then
          r = tmp_r
          k_bulk = k
          if (present(k_overall)) then
            k_overall = tmp_k_overall
          end if
        end if
      end if
    end do
  
  end subroutine grid_search_k_bulk

  subroutine init_impl(resolution, reciprocal_norms, num_work_flags, hkl, max_resolution_bins, n_reflections_in_worst_resolution_bin, min_bin_size)
 
    use xray_pure_utils, only : create_logspace_resolution_bins, assign_resolution_bin_indices, set_start_size_from_bin_index, sorted
    
    implicit none
    real(real_kind), intent(in) :: resolution(:)
    real(real_kind), intent(in) :: reciprocal_norms(3)
    integer, intent(in) :: num_work_flags
    integer, intent(in) :: hkl(3, size(resolution))
    
    integer, intent(in):: max_resolution_bins
    integer, intent(in):: n_reflections_in_worst_resolution_bin
    integer, intent(in):: min_bin_size
    
    real(real_kind) :: sorted_resolution(size(resolution))
    
    integer :: bin_start(max_resolution_bins)
    integer :: bin_size(max_resolution_bins)

    ! Preconditions
    ASSERT(is_sorted(resolution(:num_work_flags)))
    ASSERT(is_sorted(resolution(num_work_flags+1:)))
    ASSERT(size(hkl, 2) == size(resolution))
    ASSERT(n_reflections_in_worst_resolution_bin > 0)
!    ASSERT(min_bin_size >= n_reflections_in_worst_resolution_bin)
    
    n_work = num_work_flags
    sorted_resolution = sorted(resolution)
    
    call create_logspace_resolution_bins(&
        sorted_resolution, &
        n_reflections_in_worst_resolution_bin, &
        min_bin_size, &
        bin_start, bin_size, n_resolution_bins &
        )
    
    ! Allocate & copy to globals
    allocate(work_scale_resolution_bin_start(n_resolution_bins))
    allocate(work_scale_resolution_bin_size(n_resolution_bins))
    allocate(free_scale_resolution_bin_start(n_resolution_bins))
    allocate(free_scale_resolution_bin_size(n_resolution_bins))
    
    allocate(hkl_scale_resolution_bin(size(resolution)))
    
    ! Count number of high-resolution reflexes
    n_work_high_resolution = count_high_resolution_reflexes(resolution(:n_work))
    
    ! Assign resolution bin to work flags
    call assign_resolution_bin_indices(&
        resolution(:n_work), &
        sorted_resolution(bin_start(:n_resolution_bins)), &
        hkl_scale_resolution_bin(:n_work) &
        )

    call set_start_size_from_bin_index( &
        hkl_scale_resolution_bin(:n_work), &
        n_resolution_bins, &
        work_scale_resolution_bin_start, &
        work_scale_resolution_bin_size &
        )
    
    ! Assign resolution bin to free flags
    call assign_resolution_bin_indices(&
        resolution(n_work + 1:), &
        sorted_resolution(bin_start(:n_resolution_bins)), &
        hkl_scale_resolution_bin(n_work + 1:) &
        )
    
    call set_start_size_from_bin_index( &
        hkl_scale_resolution_bin(n_work + 1:), &
        n_resolution_bins, &
        free_scale_resolution_bin_start, &
        free_scale_resolution_bin_size &
        )
    
    ! Free flag indices has `n_work` offset in common [WORK|FREE] arrays
    free_scale_resolution_bin_start = free_scale_resolution_bin_start + n_work
    
    ! Allocate scaling arrays
    allocate(&
        k_bulk(size(resolution)), &
        k_iso(size(resolution)), &
        k_iso_exp(size(resolution)), &
        k_aniso(size(resolution))   &
        )
    
    ! Initialize scaling arrays
    k_bulk = 0
    k_iso = 1
    k_iso_exp = 1
    k_aniso = 1
    
    call init_MUcryst_inv(&
        hkl(1, :n_work), &
        hkl(2, :n_work), &
        hkl(3, :n_work) &
        )

    allocate(MVcryst_base(12, size(resolution)))

    call init_MVcryst_base(&
        hkl(1, :), &
        hkl(2, :), &
        hkl(3, :), &
        resolution(:) * 2.0_real_kind, &
        reciprocal_norms &
        )
  
  end subroutine init_impl
  
  subroutine init_MUcryst_inv(h, k, l)
    use xray_pure_utils, only: inverse
    implicit none
    integer, intent(in) :: h(:) !< Miller indices
    integer, intent(in) :: k(:) !< Miller indices
    integer, intent(in) :: l(:) !< Miller indices
    
    real(real_kind) :: Ucryst(6, 6)

    ! Preconditions
    ASSERT(size(h) == size(k))
    ASSERT(size(h) == size(l))
  
    Ucryst(1, 1) = sum(1.0_real_kind * h**2 * h**2)
    Ucryst(1, 2) = sum(1.0_real_kind * k**2 * h**2)
    Ucryst(1, 3) = sum(1.0_real_kind * l**2 * h**2)
    Ucryst(1, 4) = sum(1.0_real_kind * h * k * h**2)
    Ucryst(1, 5) = sum(1.0_real_kind * h * l * h**2)
    Ucryst(1, 6) = sum(1.0_real_kind * k * l * h**2)
    
    Ucryst(2, 1) = sum(1.0_real_kind * h**2 * k**2)
    Ucryst(2, 2) = sum(1.0_real_kind * k**2 * k**2)
    Ucryst(2, 3) = sum(1.0_real_kind * l**2 * k**2)
    Ucryst(2, 4) = sum(1.0_real_kind * h * k * k**2)
    Ucryst(2, 5) = sum(1.0_real_kind * h * l * k**2)
    Ucryst(2, 6) = sum(1.0_real_kind * k * l * k**2)
    
    Ucryst(3, 1) = sum(1.0_real_kind * h**2 * l**2)
    Ucryst(3, 2) = sum(1.0_real_kind * k**2 * l**2)
    Ucryst(3, 3) = sum(1.0_real_kind * l**2 * l**2)
    Ucryst(3, 4) = sum(1.0_real_kind * h * k * l**2)
    Ucryst(3, 5) = sum(1.0_real_kind * h * l * l**2)
    Ucryst(3, 6) = sum(1.0_real_kind * k * l * l**2)
    
    Ucryst(4, 1) = sum(1.0_real_kind * h**2 * h * k)
    Ucryst(4, 2) = sum(1.0_real_kind * k**2 * h * k)
    Ucryst(4, 3) = sum(1.0_real_kind * l**2 * h * k)
    Ucryst(4, 4) = sum(1.0_real_kind * h * k * h * k)
    Ucryst(4, 5) = sum(1.0_real_kind * h * l * h * k)
    Ucryst(4, 6) = sum(1.0_real_kind * k * l * h * k)
    
    Ucryst(5, 1) = sum(1.0_real_kind * h**2 * h * l)
    Ucryst(5, 2) = sum(1.0_real_kind * k**2 * h * l)
    Ucryst(5, 3) = sum(1.0_real_kind * l**2 * h * l)
    Ucryst(5, 4) = sum(1.0_real_kind * h * k * h * l)
    Ucryst(5, 5) = sum(1.0_real_kind * h * l * h * l)
    Ucryst(5, 6) = sum(1.0_real_kind * k * l * h * l)
    
    Ucryst(6, 1) = sum(1.0_real_kind * h**2 * k * l)
    Ucryst(6, 2) = sum(1.0_real_kind * k**2 * k * l)
    Ucryst(6, 3) = sum(1.0_real_kind * l**2 * k * l)
    Ucryst(6, 4) = sum(1.0_real_kind * h * k * k * l)
    Ucryst(6, 5) = sum(1.0_real_kind * h * l * k * l)
    Ucryst(6, 6) = sum(1.0_real_kind * k * l * k * l)
    
    Ucryst(4:6, :) = Ucryst(4:6, :) * 2
    Ucryst(:, 4:6) = Ucryst(:, 4:6) * 2
    Ucryst = Ucryst / size(h) ** 2
    
    MUcryst_inv = inverse(Ucryst)
  
  end subroutine init_MUcryst_inv

  subroutine init_MVcryst_base(h, k, l, resolution, reciprocal_norms)
    implicit none
    integer, intent(in) :: h(:) !< Miller indices
    integer, intent(in) :: k(:) !< Miller indices
    integer, intent(in) :: l(:) !< Miller indices
    real(real_kind), intent(in) :: resolution(:)
    real(real_kind), intent(in) :: reciprocal_norms(3)

    ! Preconditions
    call check_precondition(size(h) == size(k))
    call check_precondition(size(h) == size(l))
    call check_precondition(size(h) == size(resolution))

    MVcryst_base(1, :) = h ** 2 * reciprocal_norms(1) **2
    MVcryst_base(2, :) = h ** 2 * resolution ** 2 * reciprocal_norms(1) **2
    MVcryst_base(3, :) = k ** 2 * reciprocal_norms(2) **2
    MVcryst_base(4, :) = k ** 2 * resolution ** 2 * reciprocal_norms(2) **2
    MVcryst_base(5, :) = l ** 2 * reciprocal_norms(3) **2
    MVcryst_base(6, :) = l ** 2 * resolution ** 2 * reciprocal_norms(3) **2
    MVcryst_base(7, :) = k * l * reciprocal_norms(2) * reciprocal_norms(3)
    MVcryst_base(8, :) = k * l * resolution ** 2 * reciprocal_norms(2) * reciprocal_norms(3)
    MVcryst_base(9, :) = h * l * reciprocal_norms(1) * reciprocal_norms(3)
    MVcryst_base(10, :) = h * l * resolution ** 2 * reciprocal_norms(1) * reciprocal_norms(3)
    MVcryst_base(11, :) = h * k * reciprocal_norms(1) * reciprocal_norms(2)
    MVcryst_base(12, :) = h * k * resolution ** 2 * reciprocal_norms(1) * reciprocal_norms(2)
    MVcryst_base(7:12, :) = MVcryst_base(7:12, :) * 2.0_real_kind
  end subroutine init_MVcryst_base

  function init_MVcryst_inv(absFcalc, V_base) result(MVcryst_inv)
    use xray_pure_utils, only: inverse
    implicit none
    real(real_kind), intent(in) :: absFcalc(:), V_base(:, :)

    real(real_kind) :: Vcryst(12, 12), MVcryst_inv(12, 12)
    integer :: i, j

    ! Preconditions
    call check_precondition(size(V_base(1, :)) == size(absFcalc))
    call check_precondition(size(V_base(:, 1)) == 12)

    do i = 1, 12
      do j = i, 12
        Vcryst(i, j) = sum(absFcalc(:) ** 2 * V_base(i, :) * V_base(j, :))
        if (i /= j) then
          Vcryst(j, i) = Vcryst(i, j)
        end if
      end do
    end do

!    Vcryst(7:12, :) = Vcryst(7:12, :) * 2
!    Vcryst(:, 7:12) = Vcryst(:, 7:12) * 2
    Vcryst = Vcryst / size(absFcalc) ** 2

    MVcryst_inv = inverse(Vcryst)

  end function init_MVcryst_inv

  subroutine populate_k_bulk_linear_interpolation(absFobs, Fprot, Fbulk, k_bulk_bin, neg_s_norm2_div4, test_k_bulk)
    use xray_pure_utils, only: calc_r_factor, linear_interpolation
    implicit none
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(in) :: Fprot(size(absFobs))
    complex(real_kind), intent(in) :: Fbulk(size(absFobs))
    real(real_kind), intent(in) :: k_bulk_bin(n_resolution_bins)
    real(real_kind), intent(in) :: neg_s_norm2_div4(size(absFobs))
    real(real_kind), intent(out) :: test_k_bulk(:)
    integer :: i, s, e
    real(real_kind) :: k, b
    real(real_kind) :: s2_max, s2_min
    real(real_kind) :: r_linear, r_const
    
    ! Preconditions
    ASSERT(size(absFobs) == size(Fprot))
    ASSERT(size(absFobs) == size(Fprot))
    ASSERT(size(absFobs) == size(neg_s_norm2_div4))
    ASSERT(size(absFobs) == size(test_k_bulk))
    
    do i = n_resolution_bins, 2, -1
      s = work_scale_resolution_bin_start(i)
      e = s + work_scale_resolution_bin_size(i) - 1
      
      s2_max = abs(neg_s_norm2_div4(s))
      s2_min = abs(neg_s_norm2_div4(e))
      
      ASSERT(s2_max > s2_min)
      
      call linear_interpolation(s2_max, s2_min, k_bulk_bin(i - 1), k_bulk_bin(i), k, b)
      test_k_bulk(s:e) = max(k * abs(neg_s_norm2_div4(s:e)) + b, 0.0_real_kind)
      
      r_linear = calc_r_factor(absFobs(s:e), k_iso_exp(s:e) * k_aniso(s:e) * abs(Fprot(s:e) + test_k_bulk(s:e) * Fbulk(s:e)))
      r_const = calc_r_factor(absFobs(s:e), k_iso_exp(s:e) * k_aniso(s:e) * abs(Fprot(s:e) + k_bulk_bin(i) * Fbulk(s:e)))
      
      if (r_const < r_linear) then
        test_k_bulk(s:e) = k_bulk_bin(i)
        ! Update k_bulk for free flags
        s = free_scale_resolution_bin_start(i)
        e = s + free_scale_resolution_bin_size(i) - 1
        test_k_bulk(s:e) = k_bulk_bin(i)
      else
        ! TODO: update k_iso to match updated linear k_bulk
        ! test_k_mask(s:e) is already set
        ! Update k_bulk for free flags
        s = free_scale_resolution_bin_start(i)
        e = s + free_scale_resolution_bin_size(i) - 1
        test_k_bulk(s:e) = max(0.0d0, k * abs(neg_s_norm2_div4(s:e)) + b)
      end if
    end do
    
    ! Assign k_bulk without interpolation in last resolution bin
    ! work set
    s = work_scale_resolution_bin_start(1)
    e = s + work_scale_resolution_bin_size(1) - 1
    test_k_bulk(s:e) = k_bulk_bin(1)
    ! free set
    s = free_scale_resolution_bin_start(1)
    e = s + free_scale_resolution_bin_size(1) - 1
    test_k_bulk(s:e) = k_bulk_bin(1)
    
    ASSERT(all(test_k_bulk >=0))
  
  end subroutine populate_k_bulk_linear_interpolation
  
  
  function update_k_aniso(absFobs, absFcalc, h, k, l, r_factor) result(updated)
    use xray_pure_utils, only : calc_r_factor
    implicit none
    
    real(real_kind), intent(in) :: absFobs(:)
    real(real_kind), intent(in) :: absFcalc(size(absFobs))
    real(real_kind) :: k_aniso_V(size(absFobs))
    integer, dimension(size(absFobs)), intent(in) :: h, k, l  !< Miller indices
    real(real_kind), intent(inout) :: r_factor
    real(real_kind) :: Uaniso(6), Vaniso(12)
    real(real_kind) :: r_U, r_V
    logical :: updated
    
    updated = .FALSE.
    
    Uaniso = calc_Uaniso(absFobs(:n_work), absFcalc(:n_work), h(:n_work), k(:n_work), l(:n_work))
    Vaniso = calc_Vaniso(absFobs(:n_work), absFcalc(:n_work), MVcryst_base(:, :n_work))
    r_U = calc_r_factor(absFobs(:n_work), calc_k_aniso(Uaniso, h(:n_work), k(:n_work), l(:n_work)) * absFcalc(:n_work))
    k_aniso_V = calc_k_aniso_from_V(Vaniso, MVcryst_base)
    r_V = calc_r_factor(absFobs(:n_work), k_aniso_V(:n_work) * absFcalc(:n_work))
    if (r_U < r_factor) then
      r_factor = r_U
      k_aniso = calc_k_aniso(Uaniso, h, k, l)
      updated = .TRUE.
    end if
    if (r_V < r_factor .and. all(k_aniso_V > 0)) then
      r_factor = r_V
      k_aniso = k_aniso_V
      updated = .TRUE.
    end if
    
    ASSERT(all(k_aniso >= 0))
    
  end function update_k_aniso
  
  function update_k_bulk_k_iso(absFobs, Fprot, Fbulk, neg_s_norm2_div4, k_bulk_in_bin, k_iso_in_bin, r_work)  result(updated)
    use xray_pure_utils, only: calc_k_overall, calc_unscaled_r_factor, smooth_k_bulk
    implicit none
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(in) :: Fprot(size(absFobs))
    complex(real_kind), intent(in) :: Fbulk(size(absFobs))
    real(real_kind), intent(in) :: neg_s_norm2_div4(size(absFobs))
    real(real_kind), intent(out) :: k_bulk_in_bin(n_resolution_bins)
    real(real_kind), intent(out) :: k_iso_in_bin(n_resolution_bins)
    real(real_kind), intent(inout) :: r_work
    logical :: updated
    
    integer, parameter :: grid_size = 101
    real(real_kind), parameter :: min_k_bulk = 0.0
    real(real_kind), parameter :: max_k_bulk = 1.0
    
    real(real_kind) :: k_iso_bin, k_bulk_bin, k_overall, r_bin
    real(real_kind) :: best_k_iso_bin, best_k_bulk_bin, best_r_bin
    real(real_kind) :: k_bulk_in_bin_smooth(n_resolution_bins)
    real(real_kind) :: absFcalc(size(absFobs))
    real(real_kind) :: test_k_bulk(size(absFobs))
    real(real_kind) :: r
    integer :: i, j, s, e
    
    updated = .FALSE.
    
    do i = 1, n_resolution_bins
      s = work_scale_resolution_bin_start(i)
      e = s + work_scale_resolution_bin_size(i) - 1
      
      best_k_bulk_bin = min_k_bulk
      absFcalc(s:e) = k_aniso(s:e) * k_iso_exp(s:e) * abs(Fprot(s:e) + min_k_bulk * Fbulk(s:e)) ! Note: k_iso is excluded
      best_k_iso_bin = calc_k_overall(absFobs(s:e), absFcalc(s:e))
      best_r_bin = calc_unscaled_r_factor(absFobs(s:e), absFcalc(s:e) * best_k_iso_bin)
      
      do j = 2, grid_size
        k_bulk_bin = min_k_bulk + (max_k_bulk - min_k_bulk) / (grid_size - 1) * (j - 1)

        absFcalc(s:e) = k_aniso(s:e) * k_iso_exp(s:e) * abs(Fprot(s:e) + k_bulk_bin * Fbulk(s:e)) ! Note: k_iso is excluded
        k_iso_bin = calc_k_overall(absFobs(s:e), absFcalc(s:e))
        r_bin = calc_unscaled_r_factor(absFobs(s:e), absFcalc(s:e) * k_iso_bin)
        
        if (r_bin < best_r_bin) then
          best_k_bulk_bin = k_bulk_bin
          best_k_iso_bin = k_iso_bin
          best_r_bin = r_bin
        end if
      end do
      
      k_bulk_in_bin(i) = best_k_bulk_bin
      k_iso_in_bin(i) = best_k_iso_bin
    end do

    k_bulk_in_bin_smooth = smooth_k_bulk(k_bulk_in_bin, 1 / sqrt(abs(neg_s_norm2_div4(work_scale_resolution_bin_start) * 4)) )
    call populate_k_bulk_linear_interpolation(absFobs, Fprot, Fbulk, k_bulk_in_bin_smooth, neg_s_norm2_div4, test_k_bulk)
    ! k_iso_in_bin = calc_k_iso_bin(absFobs, &
    !         k_aniso * k_iso_exp * abs(Fprot + test_k_bulk * Fbulk))
    absFcalc(:n_work) = k_iso_in_bin(hkl_scale_resolution_bin(:n_work)) * &
            k_aniso(:n_work) * k_iso_exp(:n_work) * &
            abs(Fprot(:n_work) + test_k_bulk(:n_work) * Fbulk(:n_work))
    k_overall = calc_k_overall(absFobs(:n_work), absFcalc(:n_work))
    absFcalc(:n_work) = k_overall * absFcalc(:n_work)
    k_iso_in_bin = k_overall * k_iso_in_bin

    r = calc_unscaled_r_factor(absFobs(:n_work), absFcalc(:n_work))
    
    if (r < r_work) then
      r_work = r
      k_iso = k_iso_in_bin(hkl_scale_resolution_bin)
      k_bulk = test_k_bulk
      updated = .TRUE.
    end if

    ASSERT(all(k_bulk >= 0))
    ASSERT(all(k_iso >= 0))
    
  end function update_k_bulk_k_iso
  
  function update_k_bulk_k_iso_via_cubic_eq(absFobs, Fprot, Fbulk, neg_s_norm2_div4, k_bulk_bin, k_iso_bin, r_work) result(updated)
    use xray_pure_utils, only: calc_k_overall, calc_r_factor, calc_k_bulk_cubic, linspace, smooth_k_bulk
    implicit none
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(in) :: Fprot(size(absFobs))
    complex(real_kind), intent(in) :: Fbulk(size(absFobs))
    real(real_kind), intent(in) :: neg_s_norm2_div4(size(absFobs))
    real(real_kind), intent(inout) :: k_bulk_bin(n_resolution_bins)
    real(real_kind), intent(inout) :: k_iso_bin(n_resolution_bins)
    real(real_kind), intent(inout) :: r_work
    
    real(real_kind) :: tmp_scale(n_work)
    complex(real_kind) :: tmp_f_calc(n_work)
    
    real(real_kind) :: test_k_bulk_bin(n_resolution_bins), &
                       smooth_test_k_bulk_bin(n_resolution_bins), &
                       test_k_iso_bin(n_resolution_bins)
    
    real(real_kind) :: test_k_bulk(size(absFobs))
    
    real(real_kind), parameter :: shift = 0.05_real_kind
    integer, parameter :: grid_size = 11
    real(real_kind) :: scale_k1, min_f_mask
    real(real_kind) :: k_bulk_cubic(3)
    real(real_kind) :: k_overall
    
    real(real_kind) :: r_best, k_best, test_r

    integer :: i, s, e
    logical :: updated
    
    updated = .FALSE.
    
    ! Grid size must be odd in order to include `shift` exactly
    ASSERT(mod(grid_size, 2) == 1) 
    
    ! Estimate scale in "high resolution bin"
    scale_k1 = calc_k_overall( &
        absFobs(:n_work_high_resolution), &
        abs(k_iso_exp(:n_work_high_resolution) * k_aniso(:n_work_high_resolution) * &
            (Fprot(:n_work_high_resolution) + k_bulk(:n_work_high_resolution) * Fbulk(:n_work_high_resolution))) &
        )
    
    !write(*, *) 'k1', scale_k1
    tmp_scale = k_aniso(1:n_work) * scale_k1 * k_iso_exp(1:n_work)
    
    do i = 1, n_resolution_bins
      s = work_scale_resolution_bin_start(i)
      e = s + work_scale_resolution_bin_size(i) - 1
      min_f_mask = minval(tmp_scale(s:e) * abs(Fbulk(s:e)))
      
      k_best = k_bulk_bin(i)
      k_overall = k_iso_bin(i)

      tmp_f_calc(s:e) = tmp_scale(s:e) * (Fprot(s:e) + k_best * Fbulk(s:e))
      r_best = calc_r_factor(absFobs(s:e), abs(tmp_f_calc(s:e)))

      if (min_f_mask > 1.0e-9) then
        k_bulk_cubic = calc_k_bulk_cubic(absFobs(s:e), tmp_scale(s:e) * Fprot(s:e), tmp_scale(s:e) * Fbulk(s:e))
        ! Note: select the best k_bulk WITHOUT k_overall (just as in cctbx)
        call grid_search_k_bulk( &
            absFobs(s:e), tmp_scale(s:e) * Fprot(s:e), tmp_scale(s:e) * Fbulk(s:e), &
            k_bulk_cubic, k_best, r_best &
        )
      end if

      k_best = min(k_best, 1.0_real_kind)
      r_best = calc_r_factor(absFobs(s:e), abs(tmp_scale(s:e) * (Fprot(s:e) + k_best * Fbulk(s:e))))

      ! fine-sample k_mask around minimum of LS to fall into minimum of R
      ! Note: grid search WITH k_overall
      call grid_search_k_bulk( &
          absFobs(s:e), tmp_scale(s:e) * Fprot(s:e), tmp_scale(s:e) * Fbulk(s:e), &
          linspace(k_best - shift, k_best + shift, grid_size, endpoint=.TRUE.), &
          k_best, r_best, k_overall, 1.0_real_kind &
      )
      
      test_k_iso_bin(i) = k_overall
      test_k_bulk_bin(i) = k_best
    end do

    smooth_test_k_bulk_bin =  smooth_k_bulk(test_k_bulk_bin, 1 / sqrt(abs(neg_s_norm2_div4(work_scale_resolution_bin_start) * 4)))
    call populate_k_bulk_linear_interpolation(absFobs, Fprot, Fbulk, smooth_test_k_bulk_bin, neg_s_norm2_div4, test_k_bulk)
    test_k_iso_bin = calc_k_iso_bin(absFobs, &
        k_aniso * k_iso_exp * abs(Fprot + test_k_bulk * Fbulk))
    test_r = calc_r_factor(absFobs(:n_work), test_k_iso_bin(hkl_scale_resolution_bin(:n_work)) * k_aniso(:n_work) * k_iso_exp(:n_work) * &
        abs(Fprot(:n_work) + test_k_bulk(:n_work) * Fbulk(:n_work)) &
    )
    if (test_r < r_work) then
      k_iso = test_k_iso_bin(hkl_scale_resolution_bin)
      k_bulk = test_k_bulk
      r_work = test_r

      ! TODO: investigate whether this makes fit better, differs from cctbx
      ! k_bulk_bin = test_k_bulk_bin
      ! k_iso_bin = test_k_iso_bin
      ! END TODO
    end if

    ASSERT(all(k_bulk >= 0))
    ASSERT(all(k_iso >= 0))
    
  end function update_k_bulk_k_iso_via_cubic_eq
  
  function update_k_iso_exp(absFobs, test_k_iso, absFcalc, neg_s_norm2_div4, r_factor) result(updated)
    use xray_pure_utils, only : exponential_fit_1d_analytical, calc_r_factor
    implicit none
    
    real(real_kind), intent(inout) :: r_factor
    real(real_kind), intent(in) :: absFobs(:)
    real(real_kind), intent(in) :: test_k_iso(size(absFobs))
    real(real_kind), intent(in) :: absFcalc(size(absFobs))
    real(real_kind), intent(in) :: neg_s_norm2_div4(size(absFobs))
    
    ! https://github.com/cctbx/cctbx_project/blob/f476ac7af391357632757bf698f269333d1756e3/mmtbx/bulk_solvent/scaler.py#L225
    real(real_kind), parameter :: b_lower_limit = -100
    
    real(real_kind) :: factor, decay, factor_kb, decay_kb
    real(real_kind) :: r, r_kb
    logical :: updated
    
    updated = .FALSE.
    
    call exponential_fit_1d_analytical(absFobs(:n_work), absFcalc(:n_work), abs(neg_s_norm2_div4(:n_work)), factor, decay)

    if (decay < b_lower_limit) then
      return
    end if

    call grid_search_k_iso_exp(absFobs(:n_work), absFcalc(:n_work), neg_s_norm2_div4(:n_work), b_lower_limit, factor_kb, decay_kb)

    r = calc_r_factor(absFobs(:n_work), calc_k_iso_exp(neg_s_norm2_div4(:n_work), factor, decay) * test_k_iso(:n_work) * absFcalc(:n_work))
    r_kb = calc_r_factor(absFobs(:n_work), calc_k_iso_exp(neg_s_norm2_div4(:n_work), factor_kb, decay_kb) * test_k_iso(:n_work) * absFcalc(:n_work))

    if (r_kb <= r) then
      r = r_kb
      factor = factor_kb
      decay = decay_kb
    end if

    if (r < r_factor) then
      r_factor = r
      k_iso_exp = calc_k_iso_exp(neg_s_norm2_div4, factor, decay)
      updated = .TRUE.
    end if

    ASSERT(all(k_iso_exp >= 0))
  
  end function update_k_iso_exp
  
  function get_f_scale(n_hkl) result(result)
    implicit none
    integer, intent(in) :: n_hkl
    real(real_kind) :: result(n_hkl)
    
    ASSERT(size(k_iso) == n_hkl)
    ASSERT(size(k_iso_exp) == n_hkl)
    ASSERT(size(k_aniso) == n_hkl)
    
    result = k_iso * k_iso_exp * k_aniso
  end function get_f_scale
  
end module xray_scaling_impl_cpu_module
