!> FFTW module using Fortran 2003 interfaces.
module FFTW3
  use, intrinsic :: iso_c_binding
#  ifdef MPI
#    ifdef MKL
        include 'fftw/fftw3-mpi.f03'
#    else
        include 'fftw/fftw3-mpi.f03'
#    endif

#  else
#    ifdef MKL
        include 'fftw/fftw3.f03'
#    else
        include 'fftw3.f03'
#    endif

#  endif
end module FFTW3
