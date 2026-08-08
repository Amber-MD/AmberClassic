module xray_non_bulk_data_module
  
  use xray_pure_utils, only: real_kind
  
  implicit none
  
  public
  
  integer, pointer :: hkl(:, :)
  real(real_kind), pointer :: mSS4(:)
  real(real_kind), pointer :: xyz(:, :)
  real(real_kind), allocatable :: b_factor(:)
  integer, allocatable :: scatter_type_index(:)
  real(real_kind), allocatable :: occupancy(:)
  
  ! automatic
  complex(real_kind), allocatable :: F_non_bulk(:)
  real(real_kind), allocatable :: f(:)
  real(real_kind), allocatable :: angle(:)

  integer na, nb, nc
  
end module xray_non_bulk_data_module
