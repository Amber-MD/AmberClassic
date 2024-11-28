!<compile=optimized>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!A simple class to handle communicating messages, warnings and errors to the
!!!user.  Provides a consistent look for the user.  Three types of reports are
!!!available:
!!!o message: Simple message to the user.
!!!o warning: A message with the prefix "WARNING>".  The program will return 
!!!           unless stop_on_warn is set to true.
!!!o error  : A message with the prefix "ERROR>".  The program will stop 
!!!           immediately 
!!!
!!!The primary limitation is that the 'stop' function cannot be specified. Some
!!!programs have a prefered method for calling stop and this does not work here.
!!!
!!!For MPI code, any process may call any procedure; however, only the master 
!!!process will write anything.  All processes will call stop instructions, if
!!!applicable.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "../include/dprec.fh"
  module rism_report_c
    !munit : fortran unit number for message output
    !wunit : fortran unit number for warning output
    !eunit : fortran unit number for error output
    integer, private ::  munit=6, wunit=6, eunit=0

    !stopOnWarn : .true. - program will stop on warnings
    !            .false. - program will not stop on warnings    
    logical, private :: stopOnWarn = .false.

    !mpirank :: process id
    !mpicomm :: communicator
    !mpisize :: number of processes
#ifdef MPI
    integer, private :: mpirank=0, mpicomm=0, mpisize=1
#endif

    !Fortran does not support variable number arguments so we have bunch of 
    !subroutines to cover all the required instances.  Add as needed.

    interface rism_report_message
       module procedure message, message_i, message_is, message_isi,&
            message_isis,&
            message_isisi, message_isisisi, message_isrsi, &
            message_r, message_rs, message_rsr, message_rsrsr, &
            message_rsrsrsr, message_ra
    end interface rism_report_message

    interface rism_report_warn
       module procedure warn, warn_r, warn_rsr
    end interface rism_report_warn

    interface rism_report_error
       module procedure error, error_i, error_ia, error_is, error_isi, error_isisisi,&
            error_r, error_rs, error_rsr, error_rsrsr, error_rsrsrsr, error_ra, error_l_s_ia_s_ia
    end interface rism_report_error

    interface mwrite
       module procedure mwrite_, mwrite_i, mwrite_ia, mwrite_is, mwrite_isi, mwrite_isis,&
            mwrite_isisi, mwrite_isisisi,&
            mwrite_isrsi, mwrite_r, mwrite_rs, mwrite_rsr, mwrite_rsrsr, mwrite_rsrsrsr,&
            mwrite_ra, mwrite_l_s_ia_s_ia
    end interface mwrite
  contains

#ifdef MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets up MPI.  This should be called for any MPI code that uses this class
!!!IN:
!!!   comm : MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism_report_MPI(comm)
    use mpi
    implicit none
    integer, intent(in) :: comm
    integer :: err

    mpicomm = comm
    if(comm == MPI_COMM_NULL)&
         call rism_report_error("RISM3D_REPORT: received NULL MPI communicator")
    call mpi_comm_rank(comm,mpirank,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_REPORT: could not get MPI rank for communicator")
    call mpi_comm_size(comm,mpisize,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_REPORT: could not get MPI size for communicator")
    end subroutine rism_report_MPI
#endif /*MPI*/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets the unit number for output messages
!!!IN:
!!!   unit : unit number to write messages to
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism_report_setMUnit(unit)
      implicit none
      integer, intent(in) :: unit
      munit = unit
    end subroutine rism_report_setMUnit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets the unit number for warning messages
!!!IN:
!!!   unit : unit number to write warnings to
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism_report_setWUnit(unit)
      implicit none
      integer, intent(in) :: unit
      wunit = unit
    end subroutine rism_report_setWUnit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets the unit number for error messages
!!!IN:
!!!   unit : unit number to write error messages to
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism_report_setEUnit(unit)
      implicit none
      integer, intent(in) :: unit
      eunit = unit
    end subroutine rism_report_setEUnit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the unit number for output messages
!!!OUT:
!!!   unit number to write messages to
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism_report_getMUnit() result(unit)
      implicit none
      integer :: unit
      unit = munit
    end function rism_report_getMUnit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the unit number for warning messages
!!!OUT:
!!!   unit number to write warnings to
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism_report_getWUnit() result(unit)
      implicit none
      integer :: unit
      unit = wunit
    end function rism_report_getWUnit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the unit number for error messages
!!!OUT:
!!!   unit number to write errors to
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism_report_getEUnit() result(unit)
      implicit none
      integer :: unit
      unit = eunit
    end function rism_report_getEUnit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets whether warnings will stop the program.
!!!IN:
!!!   halt : .true. - program will stop on warnings
!!!          .false. - program will not stop on warnings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism_report_stopOnWarning(halt)
      implicit none
      logical, intent(in) :: halt
      stopOnWarn = halt
    end subroutine rism_report_stopOnWarning

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Flushes all the I/O units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism_report_flush()
      implicit none
      call flush(munit)
      call flush(wunit)
      call flush(eunit)
    end subroutine rism_report_flush

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   string : message
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message(string)
      implicit none
      character(len=*), intent(in) :: string
      call mwrite(munit, string)
    end subroutine message

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_i(form,string,i1)
      implicit none
      character(len=*), intent(in) :: form, string
      integer, intent(in) :: i1
      call mwrite(munit,form, string,i1)
    end subroutine message_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_is(form,string,i1,s2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      integer, intent(in) :: i1
      call mwrite(munit,form, string,i1,s2)
    end subroutine message_is

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   i2     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_isi(form,string,i1,s2,i2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      integer, intent(in) :: i1,i2
      call mwrite(munit,form, string,i1,s2,i2)
    end subroutine message_isi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   i2     : integer
!!!   s3     : string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_isis(form,string,i1,s2,i2,s3)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3
      integer, intent(in) :: i1,i2
      call mwrite(munit,form, string,i1,s2,i2,s3)
    end subroutine message_isis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   i2     : integer
!!!   s3     : string
!!!   i3     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_isisi(form,string,i1,s2,i2,s3,i3)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3
      integer, intent(in) :: i1,i2,i3
      call mwrite(munit,form, string,i1,s2,i2,s3,i3)
    end subroutine message_isisi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   i2     : integer
!!!   s3     : string
!!!   i3     : integer
!!!   s4     : string
!!!   i4     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_isisisi(form,string,i1,s2,i2,s3,i3,s4,i4)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3,s4
      integer, intent(in) :: i1,i2,i3,i4
      call mwrite(munit,form, string,i1,s2,i2,s3,i3,s4,i4)
    end subroutine message_isisisi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   r1     : real
!!!   s3     : string
!!!   i2     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_isrsi(form,string,i1,s2,r1,s3,i2)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3
      _REAL_, intent(in) :: r1
      integer, intent(in) :: i1,i2
      call mwrite(munit,form, string,i1,s2,r1,s3,i2)
    end subroutine message_isrsi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_r(form,string,r1)
      implicit none
      character(len=*), intent(in) :: form, string
      _REAL_, intent(in) :: r1
      call mwrite(munit, form, string, r1)
    end subroutine message_r
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_rs(form,string,r1,s2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      _REAL_, intent(in) :: r1
      call mwrite(munit,form, string,r1,s2)
    end subroutine message_rs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!   r2     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_rsr(form,string,r1,s2,r2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      _REAL_, intent(in) :: r1, r2
      call mwrite(munit,form, string,r1,s2,r2)
    end subroutine message_rsr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!   r2     : real
!!!   s3     : string
!!!   r3     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_rsrsr(form,string,r1,s2,r2,s3,r3)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3
      _REAL_, intent(in) :: r1, r2,r3
      call mwrite(munit,form, string,r1,s2,r2,s3,r3)
    end subroutine message_rsrsr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!   r2     : real
!!!   s3     : string
!!!   r3     : real
!!!   s4     : string
!!!   r4     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_rsrsrsr(form,string,r1,s2,r2,s3,r3,s4,r4)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3,s4
      _REAL_, intent(in) :: r1, r2,r3,r4
      call mwrite(munit,form, string,r1,s2,r2,s3,r3,s4,r4)
    end subroutine message_rsrsrsr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   ra1    : real array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine message_ra(form,string,ra1)
      implicit none
      character(len=*), intent(in) :: form, string
      _REAL_, intent(in) :: ra1(:)
      call mwrite(munit,form, string,ra1)
    end subroutine message_ra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   string : message
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_(unit,string)
      implicit none
      character(len=*), intent(in) :: string
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,'(a)') string
    end subroutine mwrite_

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_i(unit,form,string,i1)
      implicit none
      character(len=*), intent(in) :: form, string
      integer, intent(in) :: i1
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,i1
    end subroutine mwrite_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   ia1     : integer array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_ia(unit,form,string,ia1)
      implicit none
      character(len=*), intent(in) :: form, string
      integer, intent(in) :: ia1(:)
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,ia1
    end subroutine mwrite_ia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_is(unit,form,string,i1,s2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      integer, intent(in) :: i1
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,i1,s2
    end subroutine mwrite_is

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   i2     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_isi(unit,form,string,i1,s2,i2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      integer, intent(in) :: i1,i2
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,i1,s2,i2
    end subroutine mwrite_isi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   i2     : integer
!!!   s3     : string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_isis(unit,form,string,i1,s2,i2,s3)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3
      integer, intent(in) :: i1,i2
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,i1,s2,i2,s3
    end subroutine mwrite_isis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   i2     : integer
!!!   s3     : string
!!!   i3     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_isisi(unit,form,string,i1,s2,i2,s3,i3)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3
      integer, intent(in) :: i1,i2,i3
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,i1,s2,i2,s3,i3
    end subroutine mwrite_isisi
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   i2     : integer
!!!   s3     : string
!!!   i3     : integer
!!!   s4     : string
!!!   i4     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_isisisi(unit,form,string,i1,s2,i2,s3,i3,s4,i4)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3,s4
      integer, intent(in) :: i1,i2,i3,i4
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,i1,s2,i2,s3,i3,s4,i4
    end subroutine mwrite_isisisi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   r1     : real
!!!   s3     : string
!!!   i2     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_isrsi(unit,form,string,i1,s2,r1,s3,i2)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3
      _REAL_, intent(in) :: r1
      integer, intent(in) :: i1,i2
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,i1,s2,r1,s3,i2
    end subroutine mwrite_isrsi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!   r2     : real
!!!   s3     : string
!!!   r3     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_rsrsr(unit,form,string,r1,s2,r2,s3,r3)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3
      _REAL_, intent(in) :: r1, r2,r3
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,r1,s2,r2,s3,r3
    end subroutine mwrite_rsrsr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!   r2     : real
!!!   s3     : string
!!!   r3     : real
!!!   s4     : string
!!!   r4     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_rsrsrsr(unit,form,string,r1,s2,r2,s3,r3,s4,r4)
      implicit none
      character(len=*), intent(in) :: form, string, s2,s3,s4
      _REAL_, intent(in) :: r1, r2,r3,r4
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,r1,s2,r2,s3,r3,s4,r4
    end subroutine mwrite_rsrsrsr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_rs(unit,form,string,r1,s2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      _REAL_, intent(in) :: r1
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,r1,s2
    end subroutine mwrite_rs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!   r2     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_rsr(unit,form,string,r1,s2,r2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      _REAL_, intent(in) :: r1, r2
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,r1,s2,r2
    end subroutine mwrite_rsr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_r(unit,form,string,r1)
      implicit none
      character(len=*), intent(in) :: form, string
      _REAL_, intent(in) :: r1
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,r1
    end subroutine mwrite_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   ra1    : real array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_ra(unit,form,string,ra1)
      implicit none
      character(len=*), intent(in) :: form, string
      _REAL_, intent(in) :: ra1(:)
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string,ra1
    end subroutine mwrite_ra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a message to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   l      : logical
!!!   s2     : string
!!!   ia1    : integer array
!!!   s3     : string
!!!   ia2    : integer array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mwrite_l_s_ia_s_ia(unit, form, string, l, s2, ia1, s3, ia2)
      implicit none
      character(len=*), intent(in) :: form, string, s2, s3
      logical, intent(in) :: l
      integer, intent(in) :: ia1(:), ia2(:)
      integer, intent(in) :: unit
#ifdef MPI
      if(mpirank == 0)&
#endif /*MPI*/
           write(unit,form) string, l, s2, ia1, s3, ia2
    end subroutine mwrite_l_s_ia_s_ia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a warning to the user
!!!IN:
!!!   string : warning
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine warn(string)
      implicit none
      character(len=*), intent(in) :: string
      call mwrite(wunit,"WARNING> "//string)
      if(stopOnWarn) call halt()
    end subroutine warn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a warning to the user
!!!IN:
!!!   form : print format
!!!   string : warning
!!!   r1     : real
!!!   s2     : warning
!!!   r2     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine warn_rsr(form,string,r1,s2,r2)
      implicit none
      character(len=*), intent(in) :: form,string, s2
      _REAL_, intent(in) :: r1, r2
      call mwrite(wunit,form,"WARNING> "//string,r1,s2,r2)
      if(stopOnWarn) call halt()
    end subroutine warn_rsr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out a warning to the user
!!!IN:
!!!   form : print format
!!!   string : warning
!!!   r1     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine warn_r(form,string,r1)
      implicit none
      character(len=*), intent(in) :: form,string
      _REAL_, intent(in) :: r1
      call mwrite(wunit,form,"WARNING> "//string,r1)
      if(stopOnWarn) call halt()
    end subroutine warn_r
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   string : error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error(string)
      implicit none
      character(len=*), intent(in) :: string
      call mwrite_(eunit,"ERROR> "//string)
      call halt()
    end subroutine error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_i(form,string,i1)
      implicit none
      character(len=*), intent(in) :: form, string
      integer, intent(in) :: i1
      call mwrite_i(eunit,form,"ERROR> "//string,i1)
      call halt()
    end subroutine error_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   ia1     : integer array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_ia(form,string,ia1)
      implicit none
      character(len=*), intent(in) :: form, string
      integer, intent(in) :: ia1(:)
      call mwrite_ia(eunit,form,"ERROR> "//string,ia1)
      call halt()
    end subroutine error_ia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_is(form,string,i1,s2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      integer, intent(in) :: i1
      call mwrite_is(eunit,form,"ERROR> "//string,i1,s2)
      call halt()
    end subroutine error_is

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   i2     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_isi(form,string,i1,s2,i2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      integer, intent(in) :: i1, i2
      call mwrite_isi(eunit,form,"ERROR> "//string,i1,s2,i2)
      call halt()
    end subroutine error_isi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   i1     : integer
!!!   s2     : string
!!!   i2     : integer
!!!   s3     : string
!!!   i3     : integer
!!!   s4     : string
!!!   i4     : integer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_isisisi(form,string,i1,s2,i2,s3,i3,s4,i4)
      implicit none
      character(len=*), intent(in) :: form, string, s2, s3, s4
      integer, intent(in) :: i1, i2, i3, i4
      call mwrite_isisisi(eunit,form,"ERROR> "//string,i1,s2,i2,s3,i3,s4,i4)
      call halt()
    end subroutine error_isisisi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_r(form,string,r1)
      implicit none
      character(len=*), intent(in) :: form, string
      _REAL_, intent(in) :: r1
      call mwrite_r(eunit,form,"ERROR> "//string,r1)
      call halt()
    end subroutine error_r
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_rs(form,string,r1,s2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      _REAL_, intent(in) :: r1
      call mwrite_rs(eunit,form,"ERROR> "//string,r1,s2)
      call halt()
    end subroutine error_rs
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!   r2     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_rsr(form,string,r1,s2,r2)
      implicit none
      character(len=*), intent(in) :: form, string, s2
      _REAL_, intent(in) :: r1, r2
      call mwrite_rsr(eunit,form,"ERROR> "//string,r1,s2,r2)
      call halt()
    end subroutine error_rsr
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!   r2     : real
!!!   s3     : string
!!!   r3     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_rsrsr(form,string,r1,s2,r2,s3,r3)
      implicit none
      character(len=*), intent(in) :: form, string, s2, s3
      _REAL_, intent(in) :: r1, r2, r3
      call mwrite_rsrsr(eunit,form,"ERROR> "//string,r1,s2,r2,s3,r3)
      stop 1
    end subroutine error_rsrsr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   r1     : real
!!!   s2     : string
!!!   r2     : real
!!!   s3     : string
!!!   r3     : real
!!!   s4     : string
!!!   r4     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_rsrsrsr(form,string,r1,s2,r2,s3,r3,s4,r4)
      implicit none
      character(len=*), intent(in) :: form, string, s2, s3,s4
      _REAL_, intent(in) :: r1, r2, r3,r4
      call mwrite_rsrsrsr(eunit,form,"ERROR> "//string,r1,s2,r2,s3,r3,s4,r4)
      stop 1
    end subroutine error_rsrsrsr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   ra1    : real array
!!!   s2     : string
!!!   r2     : real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_ra(form,string,ra1)
      implicit none
      character(len=*), intent(in) :: form, string
      _REAL_, intent(in) :: ra1(:)
      call mwrite_ra(eunit,form,"ERROR> "//string,ra1)
      call halt()
    end subroutine error_ra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Write out an error to the user
!!!IN:
!!!   form : print format
!!!   string : error
!!!   l      : logical
!!!   s2     : string
!!!   ia1    : integer array
!!!   s3     : string
!!!   ia2    : integer array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine error_l_s_ia_s_ia(form, string, l, s2, ia1, s3, ia2)
      implicit none
      character(len=*), intent(in) :: form, string, s2, s3
      logical, intent(in) :: l
      integer, intent(in) :: ia1(:), ia2(:)
      call mwrite_l_s_ia_s_ia(eunit, form, "ERROR> "//string, l, s2, ia1, s3, ia2)
      call halt()
    end subroutine error_l_s_ia_s_ia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Halts the program with an optional error code
!!!IN:
!!!   o_code: (optional) integer error code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine halt(o_code)
      implicit none
      integer, optional, intent(in) :: o_code
      integer :: code
      ! integer, pointer :: p=>NULL()
      code =1
      if(present(o_code)) code = o_code
      stop 1
      !using the below statement instead of 'stop' will trigger a 
      !segfault, which can be useful if tracebacks are enabled through 
      !the compiler
      !p = code + p
    end subroutine halt

  end module rism_report_c
