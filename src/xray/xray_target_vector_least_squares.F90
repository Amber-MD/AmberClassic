module xray_target_vector_least_squares_module

!   #ifdef CUDA
!     use xray_target_vector_least_squares_impl_gpu_module
!   #else
    use xray_target_vector_least_squares_impl_cpu_module
!   #endif

    ! Public module interface
    public :: init
    public :: calc_partial_d_target_d_absFcalc
    public :: finalize

end module xray_target_vector_least_squares_module
