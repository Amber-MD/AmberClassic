!> FFTW module using Fortran 2003 interfaces.
module FFTW3
  use, intrinsic :: iso_c_binding
#  ifdef MPI
     include 'fftw3-mpi.f03'
#  else
     include 'fftw3.f03'
#  endif
end module FFTW3
