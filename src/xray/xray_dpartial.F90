module xray_dpartial_module

#ifdef CUDA
  use xray_dpartial_impl_gpu_module
#else
  use xray_dpartial_impl_cpu_module
#endif

  implicit none
  private
  
  public :: calc_partial_d_target_d_frac
  public :: calc_partial_d_vls_d_frac
  public :: finalize
  public :: init
  
end module xray_dpartial_module
