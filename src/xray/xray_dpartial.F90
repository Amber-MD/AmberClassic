module xray_dpartial_module

  use xray_dpartial_impl_cpu_module

  implicit none
  private
  
  public :: calc_partial_d_target_d_frac
  public :: calc_partial_d_vls_d_frac
  public :: finalize
  public :: init
  
end module xray_dpartial_module
