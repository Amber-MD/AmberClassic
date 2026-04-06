module xray_fft_interface_module
  
  implicit none
  private
  
  public :: dft_i2c_3d

contains
  
  subroutine dft_i2c_3d(dim, data_in, data_out)
  
    use iso_fortran_env
    implicit none
#ifdef MKL
    include 'fftw/fftw3.f'
#else
    include 'fftw3.f'
#endif
    
    integer, intent(in) :: dim(3)
    integer, allocatable, intent(in) :: data_in(:)
    complex(real64), allocatable, intent(inout) :: data_out(:)
    
    double precision :: mask_bs_grid_3d(dim(3), dim(2), dim(1))
    double complex ::   mask_bs_grid_3d_fft(dim(3)/2 + 1, dim(2), dim(1))
    integer(int64) :: plan_forward(1)
    
    mask_bs_grid_3d = reshape(data_in, [dim(3) , dim(2), dim(1)])
    call dfftw_plan_dft_r2c_3d(plan_forward, dim(3), &
        dim(2), dim(1), mask_bs_grid_3d, &
        mask_bs_grid_3d_fft, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan_forward, mask_bs_grid_3d, &
        mask_bs_grid_3d_fft)
    data_out = reshape(mask_bs_grid_3d_fft, &
        [(dim(3) / 2 + 1) * dim(2) * dim(1)])
    call dfftw_destroy_plan(plan_forward)

  end subroutine dft_i2c_3d
  
end module xray_fft_interface_module
