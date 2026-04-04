#include "../include/assert.fh"

module xray_interface2_data_module
  
  use xray_contracts_module
  use xray_unit_cell_module
  use xray_pure_utils, only : real_kind
  use xray_non_bulk_data_module, only : na, nb, nc
  implicit none
  
  public
  
  !! All input arrays are reordered into [work, free] set,
  !! i.e. "work" reflexes alaways precede "free" ones.
  !! All reflexes within "work" and "free" partitions are sorted by resolution
  !! from best (lower values) to worst (higher values)
  !! Initial reflex oreder is kept in hkl_io_order

  !! Reflex data
  integer, save :: n_hkl
  integer, save :: n_work
  integer, allocatable, save :: hkl(:, :)           ! size = (3, n_hkl)
  integer, allocatable, save :: hkl_io_order(:)  ! Initial input order of reflexes ; size = (n_hkl)
  integer, allocatable, save :: new_order(:)     ! Sorted  order of reflexes ; size = (n_hkl)
  complex(real_kind), allocatable, save :: Fobs(:)  ! size = (n_hkl)
  complex(real_kind), allocatable, save :: Fcalc(:) ! size = (n_hkl)
  real(real_kind), allocatable, save :: sigma_Fobs(:)  ! size = (n_hkl)
  real(real_kind), allocatable, save :: resolution(:)
  real(real_kind), allocatable, save :: penalty(:)  ! per-reflection penalty
  integer :: spacegroup_number

  !! Atomic data
  integer, save :: n_atom
  type(unit_cell_t), save :: unit_cell
  logical, allocatable, save :: atom_is_not_bulk(:)
  integer, allocatable, save :: non_bulk_atom_indices(:)
  real(real_kind), allocatable, save :: atom_b_factor(:)
  integer, allocatable, save :: atom_scatter_type(:)
  real(real_kind), allocatable, save :: atom_occupancy(:)
  real(real_kind), allocatable, save :: scatter_coefficients(:, :, :) ! Fourier coefficients (2,n_scatter_coeffs,n_scatter_types)

  !! info about bulk solvent model
  integer, save :: model_id

  !! options for ls target
  real(real_kind), save :: ls_r3, ls_r4

contains
  
  subroutine init(input_hkl, input_Fobs, input_sigma_Fobs, input_work_flag, &
      input_unit_cell, input_scatter_coefficients, &
      input_atom_b_factor, input_atom_occupancy, input_atom_scatter_type, &
      input_atom_selection, r3, r4 , input_spacegroup_number, iix, iiy, iiz )
    use xray_pure_utils, only: index_partition, index_sort, calc_resolution, pack_index
    
    implicit none
    
    integer, intent(in) :: input_hkl(:, :)
    complex(real_kind), intent(in) :: input_Fobs(:)
    real(real_kind), intent(in) :: input_sigma_Fobs(:)
    logical, intent(in) :: input_work_flag(:)
    class(unit_cell_t), intent(in) :: input_unit_cell
    real(real_kind), intent(in) :: input_scatter_coefficients(:,:,:) ! Fourier coefficients (2,n_scatter_coeffs,n_scatter_types)
    real(real_kind), intent(in) :: input_atom_b_factor(:)
    real(real_kind), intent(in) :: input_atom_occupancy(:)
    integer, intent(in) :: input_atom_scatter_type(:)
    logical, intent(in) :: input_atom_selection(:)
    real(real_kind), intent(in) :: r3, r4
    integer, intent(in) :: input_spacegroup_number, iix, iiy, iiz
    
    ! locals
    integer :: j
    
    
    ASSERT(size(input_hkl, 1) == 3)
    ASSERT(size(input_hkl, 2) == size(input_Fobs))
    ASSERT(size(input_hkl, 2) == size(input_sigma_Fobs))
    ASSERT(size(input_hkl, 2) == size(input_work_flag))
    ASSERT(size(input_atom_b_factor) == size(input_atom_occupancy))
    ASSERT(size(input_atom_b_factor) == size(input_atom_scatter_type))
    ASSERT(size(input_scatter_coefficients, 1) == 2)
    ASSERT(minval(input_atom_b_factor, input_atom_selection) >= 0)
    ASSERT(all(input_atom_occupancy <= 1.0))
    ASSERT(all(input_atom_occupancy >= 0.0))
    ASSERT(minval(input_atom_scatter_type) >= 1)
    ASSERT(maxval(input_atom_scatter_type) <= size(input_scatter_coefficients, 3))
    
    n_hkl = size(input_hkl, 2)
    n_atom = size(input_atom_occupancy)
    n_work = count(input_work_flag)
    
    ASSERT(n_work > 0) ! Number of 'work' reflections must be >=1
    ASSERT(n_hkl - n_work > 0) ! Number of 'free' reflections must be >=1
    
    unit_cell = input_unit_cell
  
    resolution = 1 / sqrt(unit_cell%get_s2(input_hkl))
    allocate(penalty(n_hkl))

    allocate(new_order(n_hkl))
    allocate(hkl_io_order(n_hkl))
    
    do j = 1, n_hkl
      new_order(j) = j
      hkl_io_order(j) = j
    end do
    
    ! Partition and sort reflexes by resolution
    call index_partition(input_work_flag, new_order)
    call index_sort(resolution, new_order(:n_work))
    call index_sort(resolution, new_order(n_work + 1:))

    ! (note: following lines invoke an automatic allocation of the new arrays;
    !   cf. Section 6.7 of Metcalf, Modern Fortran Explained)
    hkl = input_hkl(:, new_order)
    Fobs = input_Fobs(new_order)
    sigma_Fobs = input_sigma_Fobs(new_order)
    resolution = resolution(new_order)
    hkl_io_order(new_order) = hkl_io_order
    
    allocate(Fcalc(n_hkl))

    atom_is_not_bulk = input_atom_selection
    non_bulk_atom_indices = pack_index(input_atom_selection)
    
    atom_b_factor = input_atom_b_factor
    atom_scatter_type = input_atom_scatter_type
    atom_occupancy = input_atom_occupancy
    scatter_coefficients = input_scatter_coefficients

    ASSERT(size(atom_b_factor) == n_atom)
    ASSERT(size(atom_occupancy) == n_atom)
    ASSERT(size(atom_scatter_type) == n_atom)
    ASSERT(size(atom_is_not_bulk) == n_atom)

    ASSERT(size(hkl, 1) == 3)
    ASSERT(size(hkl, 2) == n_hkl)
    ASSERT(size(Fobs) == n_hkl)
    ASSERT(size(sigma_Fobs) == n_hkl)
    ASSERT(size(hkl_io_order) == n_hkl)

    ls_r3 = r3
    ls_r4 = r4
    spacegroup_number = input_spacegroup_number
    na = iix
    nb = iiy
    nc = iiz
    
  end subroutine init
  
  subroutine finalize()
  
    if(allocated(hkl)) deallocate(hkl)
    if(allocated(hkl_io_order)) deallocate(hkl_io_order)
    if(allocated(Fobs)) deallocate(Fobs)
    if(allocated(Fcalc)) deallocate(Fcalc)
    
    if(allocated(non_bulk_atom_indices)) deallocate(non_bulk_atom_indices)
    if(allocated(atom_is_not_bulk)) deallocate(atom_is_not_bulk)
    if(allocated(atom_scatter_type)) deallocate(atom_scatter_type)
    if(allocated(atom_occupancy)) deallocate(atom_occupancy)
    if(allocated(atom_b_factor)) deallocate(atom_b_factor)
    if(allocated(scatter_coefficients)) deallocate(scatter_coefficients)
    
  end subroutine finalize
end module xray_interface2_data_module
