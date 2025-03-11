!<compile=optimized>

#include "../include/dprec.fh"

!> Functions and subroutines to safely allocate, reallocate (data preserving)
!! and deallocate pointers.
!!
!! Typical usage example:
!!
!!   _REAL_, pointer :: p(:)=>NULL()
!!   call safemem_realloc(p, N)
!!
!!   ! Somewhere else, not necessarily in the same scope ('p' may be a
!!   ! different pointer handle here), resize the array.
!!   call safemem_realloc(p, M)
!!
!!   ! Somewhere else, not necessarily in the same scope ('p' may be a
!!   ! different pointer handle here), deallocate the memory.
!!   if (safemem_dealloc(p) /= 0) then
!!     ! Deal with error.
!!   end if
!!
!!   ! At any point you can get a summary of the current and maximum
!!   ! allocated memory in bytes.
!!   integer(kind=int64) :: memstats(10)
!!   memstats = memStatus()
!!   write(unit,'(a)') "Type         Current         Maximum"
!!   write(unit,'(a,i12,a,f12.5,a)') "Integer  ",memstats(1)," B ",&
!!        dble(memstats(6))/BYTES_PER_GB," GB"
!!   write(unit,'(a,i12,a,f12.5,a)') "Real     ",memstats(2)," B ",&
!!        dble(memstats(7))/BYTES_PER_GB," GB"
!!   write(unit,'(a,i12,a,f12.5,a)') "Logical  ",memstats(3)," B ",&
!!        dble(memstats(8))/BYTES_PER_GB," GB"
!!   write(unit,'(a,i12,a,f12.5,a)') "Character",memstats(4)," B ",&
!!        dble(memstats(9))/BYTES_PER_GB," GB"
!!   write(unit,'(a)') "---------------------------------------"
!!   write(unit,'(a,i12,a,f12.5,a)') "Total    ",memstats(5)," B ",&
!!        dble(memstats(10))/BYTES_PER_GB," GB"
!!
!! Notes on proper usage:
!! - safemem_realloc preserves data by default.  Since Fortran has no
!!   built in facility for this, it is generally slow.  If you don't
!!   need this feature, use 'preserve=.false.'.
!! - Uninitiallized pointers should be set to NULL.  Either use
!!   'p=>NULL()' or 'nullify(p)'.
!! - As long as you are careful to set unitialized data to NULL, you
!!   can mix these routine calls with your own allocated/deallocate
!!   statements on any pointer.  However, the memory summary will be
!!   invalid.
!! - If compiled and linked with FFTW, 16-byte aligned memory, useful
!!   for SIMD operations can be allocated by using
!!   'o_aligned=.true.'. If you allocate aligned memory, you _must_
!!   deallocate it with the aligned option.  It is not possible for
!!   these routines to determine how the memory was allocated.
!! - The memory summary only tracks heap memory allocated with these
!!   routines.
!! - The memory used can be compared against the resident memory.
!!   pmap -x <PID> works well here.  Note that resident memory is 
!!   (usually) memory that has been allocated and accessed and that 
!!   memory from untracked parts of the code will also be include in 
!!   the resident memory.

module safemem
  use rism_report_c
  use iso_fortran_env, only: int64
#ifndef CUDA
  use FFTW3
#else
  use iso_c_binding
#endif /*CUDA*/

  implicit none
  public
  type memTracker
     private
     ! int     :: Current number of bytes of integer memory allocated
     ! real    :: Current number of bytes of real memory allocated
     ! logical :: Current number of bytes of logical memory allocated
     ! char    :: Current number of bytes of character memory allocated
     ! total   :: Current number of bytes of total memory allocated
     integer(kind=int64) :: int = 0, real = 0, logical = 0, char = 0, total = 0
     ! maxint     :: Current number of bytes of integer memory allocated
     ! maxreal    :: Current number of bytes of real memory allocated
     ! maxlogical :: Current number of bytes of logical memory allocated
     ! maxchar    :: Current number of bytes of character memory allocated
     ! max        :: Current number of bytes of total memory allocated
     integer(kind=int64) :: maxint = 0, maxreal = 0, maxlogical = 0, maxchar = 0, max = 0
  end type memTracker

  ! BYTES_PER_GIGABYTES :: used to convert between bytes and GB
  integer(kind=int64), parameter :: BYTES_PER_GB = 1024**3
  ! BYTES_PER_MEGABYTES :: used to convert between bytes and MB
  integer(kind=int64), parameter :: BYTES_PER_MB = 1024**2
  ! BYTES_PER_KILOBYTES :: used to convert between bytes and KB
  integer(kind=int64), parameter :: BYTES_PER_KB = 1024**1

  type(memTracker), private, save :: totalMem

  !> Reallocate memory of a specific pointer during the simulation. This can
  !! handle up to 5 dimensions and the optiontional last argument can indicate
  !! whether or not to preserve the contents of the array.  Preservation is the
  !! default.
  !! @param[in] p original pointer
  !! @param[in] n1 size of the first dimension
  !! @param[in] n2 size of the second dimension (optional)
  !! @param[in] n3 size of the third dimension (optional)
  !! @param[in] n4 size of the fourth dimension (optional)
  !! @param[in] n5 size of the fifth dimension (optional)
  !! @param[in] preserve boolean indicating whether or not to preserve the contents of
  !!                     the array (optional)
  interface safemem_realloc
     module procedure safemem_realloc_1d_real_int32
     module procedure safemem_realloc_1d_real
     module procedure safemem_realloc_2d_real_int32
     module procedure safemem_realloc_2d_real
     module procedure safemem_realloc_3d_real_int32
     module procedure safemem_realloc_3d_real
     module procedure safemem_realloc_4d_real_int32
     module procedure safemem_realloc_4d_real
     module procedure safemem_realloc_5d_real_int32
     module procedure safemem_realloc_5d_real
     module procedure safemem_realloc_5d_real_interp_int32
     module procedure safemem_realloc_5d_real_interp

     module procedure safemem_realloc_1d_integer_int32
     module procedure safemem_realloc_1d_integer
     module procedure safemem_realloc_1d_int64
     module procedure safemem_realloc_2d_integer_int32
     module procedure safemem_realloc_2d_integer

     module procedure safemem_realloc_1d_character
     module procedure safemem_realloc_1d_character_int32
     module procedure safemem_realloc_2d_character
     module procedure safemem_realloc_2d_character_int32

     module procedure safemem_realloc_1d_logical
     module procedure safemem_realloc_1d_logical_int32
  end interface safemem_realloc


  !> Check if the pointer is allocated before deallocating.  Returns any error code.
  !! @param[in] p pointer to be deallocated
  !! @param[out] ?? error code
  interface safemem_dealloc
     module procedure safemem_dealloc_pointer_1d_real
     module procedure safemem_dealloc_pointer_2d_real
     module procedure safemem_dealloc_pointer_3d_real
     module procedure safemem_dealloc_pointer_4d_real
     module procedure safemem_dealloc_pointer_5d_real
     module procedure safemem_dealloc_pointer_1d_int
     module procedure safemem_dealloc_pointer_2d_int
     module procedure safemem_dealloc_pointer_3d_int
     module procedure safemem_dealloc_pointer_4d_int
     module procedure safemem_dealloc_pointer_5d_int
     module procedure safemem_dealloc_pointer_1d_int64
     module procedure safemem_dealloc_pointer_1d_character
     module procedure safemem_dealloc_pointer_2d_character
     module procedure safemem_dealloc_pointer_1d_logical
  end interface safemem_dealloc


  !> Calculates the size in bytes of an array of arbitrary dimensions.
  !! @param[in] array array to get the size of
  !! @param[in] N     number of elements
  interface arraySize
     module procedure arraySize_1d_real
     module procedure arraySize_2d_real
     module procedure arraySize_3d_real
     module procedure arraySize_4d_real
     module procedure arraySize_5d_real
     module procedure arraySize_1d_int
     module procedure arraySize_2d_int
     module procedure arraySize_3d_int
     module procedure arraySize_4d_int
     module procedure arraySize_5d_int
     module procedure arraySize_1d_int64
     module procedure arraySize_1d_character
     module procedure arraySize_2d_character
     module procedure arraySize_1d_logical
     ! module procedure arraySize_int
     ! module procedure arraySize_character
     ! module procedure arraySize_logical
  end interface arraySize
   
  !> Returns an array of current and maximum memory allocated in bytes so
  !! far.  No distinction is made between different
  !! kinds of the same type:
  !! (1) : Current integer memory
  !! (2) : Current real memory
  !! (3) : Current logical memory
  !! (4) : Current char memory
  !! (5) : Current logical memory
  !! (6) : Max integer memory
  !! (7) : Max real memory
  !! (8) : Max logical memory
  !! (9) : Max char memory
  !! (10): Max logical memory
  interface memStatus
     module procedure memStatus_obj, memStatus_global
  end interface memStatus


  public safemem_realloc, safemem_dealloc, memStatus


  private memadd_i, memadd_r, memadd_c, memadd_l


contains


  !> Given a value x,y,z and arrays of coordinates x1, x2, x3, this function locates the index of
  !! x1,2,3 that is closest to x, y, z
  !! IN:
  !!    x1 :: values in dimension 1 (ascending order)
  !!    x2 :: values in dimension 2 (ascending order)
  !!    x3 :: values in dimension 3 (ascending order)
  !!    n1 :: number of values
  !!    n2 :: number of values
  !!    n3 :: number of values
  !!    x  :: x value
  !!    y  :: y value
  !!    z  :: z value
  !! OUT:
  !!     three
  function nearestIndex(x1,x2,x3,n1,n2,n3,x,y,z)
    implicit none
    integer, intent(in) :: n1, n2, n3
    _REAL_, intent(in) :: x1(n1),x2(n2),x3(n3),x,y,z
    integer :: nearestIndex(3)
    integer :: i
    nearestIndex = 1
    do i = 2, n1
       if (abs(x1(i)-x) < abs(x1(nearestindex(1))-x)) then
          nearestIndex(1) = i
       end if
    end do
    do i = 2, n2
       if (abs(x2(i)-y) < abs(x2(nearestindex(2))-y)) then
          nearestIndex(2) = i
       end if
    end do
    do i = 2, n3
       if (abs(x3(i)-z) < abs(x3(nearestindex(3))-z)) then
          nearestIndex(3) = i
       end if
    end do
    return
  end function nearestIndex


  !> Reallocate memory of a specific pointer during the simulation. This can
  !! handle 1 dimension and preserves the contents of the array.
  !! IN:
  !!    p        :: original pointer
  !!    n1       :: size of the first dimension
  !!    o_preserve :: (optional) boolean indicating whether or not to
  !!                  preserve the contents of the array
  !!    o_aligned :: (optional) Use the FFTW library to allocate 16-byte
  !!                 aligned memory necessary to SIMD operations
  function safemem_realloc_1d_real_int32(p, n1, o_preserve, o_aligned)
    implicit none
    _REAL_, pointer, dimension(:) :: safemem_realloc_1d_real_int32
    _REAL_, pointer, dimension(:):: p
    logical, optional, intent(in) :: o_preserve,o_aligned
    integer, intent(in) :: n1
    safemem_realloc_1d_real_int32 => safemem_realloc_1d_real(p, int(n1,kind=int64), o_preserve, o_aligned)    
  end function safemem_realloc_1d_real_int32
  
  function safemem_realloc_1d_real(p, n1, o_preserve, o_aligned)
    implicit none
    _REAL_, pointer, dimension(:) :: safemem_realloc_1d_real
    _REAL_, pointer, dimension(:):: p
    logical, optional, intent(in) :: o_preserve,o_aligned
    integer(kind=int64), intent(in) :: n1
    logical :: preserve,aligned
    integer(kind=int64) :: i,nold1
    integer :: err

    type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
    write(6,*) "safemem_realloc_1d_real", associated(p),n1,preserve
    call flush(6)
#endif /*RISM_DEBUG*/
    aligned = .false.
    if (present(o_aligned)) aligned = o_aligned

    preserve = .true.
    if (present(o_preserve)) preserve = o_preserve

    if (.not.preserve) then
       if (safemem_dealloc(p,aligned) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 1D REAL")
    end if
    !get new memory
    if (aligned) then
#ifndef CUDA
       err = 0
       cptr = fftw_alloc_real(int(n1, C_SIZE_T))
       if (.not.c_associated(cptr)) then
          err = 1
       else
          call c_f_pointer(cptr,safemem_realloc_1d_real,[n1])
       end if
#else
       call rism_report_error('1: Cannot allocate aligned memory when using CUDA. Check you Fortran calls to safemem.')
#endif /*CUDA*/
    else
       allocate(safemem_realloc_1d_real(n1), STAT=err)
    end if
    
    if (err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 1D REAL")
    call memadd_r(totalmem,arraySize(safemem_realloc_1d_real))

    !transfer data to the new memory
    if (associated(p) .and. preserve .and. n1 >0) then
       nold1 = min(ubound(p,1), n1)
       call dcopy(nold1, p(1), 1, &
            safemem_realloc_1d_real(1),1)
       safemem_realloc_1d_real(1:nold1) = p(1:nold1)
       !free old memory
       if (safemem_dealloc(p,aligned) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 1D REAL")
    end if
  end function safemem_realloc_1d_real


  !> Reallocate memory of a specific pointer during the simulation. This can
  !! handle 2 dimensions and preserves the contents of the array.
  !! IN:
  !!    p     :: original pointer
  !!    n1    :: size of the first dimension
  !!    n2    :: size of the second dimension
  !!    o_preserve :: (optional) boolean indicating whether or not to
  !!                  preserve the contents of the array
  !!    o_aligned :: (optional) Use the FFTW library to allocate 16-byte
  !!                 aligned memory necessary to SIMD operations
  function safemem_realloc_2d_real_int32(p,n1, n2,o_preserve,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:) :: safemem_realloc_2d_real_int32
    _REAL_, pointer, dimension(:,:) :: p
    logical,optional, intent(in) :: o_preserve,o_aligned
    integer, intent(in) :: n1,n2
    safemem_realloc_2d_real_int32 => safemem_realloc_2d_real(p, int(n1, kind=int64), int(n2, kind=int64), o_preserve, o_aligned)
  end function safemem_realloc_2d_real_int32

  function safemem_realloc_2d_real(p,n1, n2,o_preserve,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:) :: safemem_realloc_2d_real
    _REAL_, pointer, dimension(:,:) :: p
    logical,optional, intent(in) :: o_preserve,o_aligned
    integer(kind=int64), intent(in) :: n1,n2
    logical :: preserve, aligned
    integer(kind=int64) :: nold1, nold2
    integer :: err
    integer(kind=int64) :: i2
    type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
    write(6,*) "safemem_realloc_2d_real", associated(p),n1,n2,preserve
    call flush(6)
#endif /*RISM_DEBUG*/

    aligned = .false.
    if (present(o_aligned)) aligned = o_aligned

    preserve=.true.
    if (present(o_preserve)) preserve = o_preserve

    if (.not.preserve) then
       if (safemem_dealloc(p,aligned) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 2D REAL 1")
    end if
    if (aligned) then
#ifndef CUDA
       err = 0
       cptr = fftw_alloc_real(int(n1*n2, C_SIZE_T))
       if (.not.c_associated(cptr)) then
          err = 1
       else
          call c_f_pointer(cptr,safemem_realloc_2d_real,[n1,n2])
       end if
#else
        call backtrace
        call rism_report_error('2: Cannot allocate aligned memory when using CUDA. Check you Fortran calls to safemem.')
#endif /*CUDA*/
    else
       allocate(safemem_realloc_2d_real(1:n1,1:n2), STAT=err)
    end if
    if (err /= 0)&
         call rism_report_error("(a,i4)","reallocation_pointer allocate error: 2D REAL 2",err)
    call memadd_r(totalmem,arraySize(safemem_realloc_2d_real))
    if (associated(p) .and.preserve .and. n1*n2 > 0) then
       safemem_realloc_2d_real = 0d0
       nold1 = min(ubound(p,1), n1)
       nold2 = min(ubound(p,2), n2)
       do i2=1,nold2
          call dcopy(nold1, p(1, i2), 1, &
               safemem_realloc_2d_real(1, i2),1)
       end do
       if (safemem_dealloc(p,aligned) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 2D REAL 3")
    end if
  end function safemem_realloc_2d_real

  !
  !! reallocate memory of a specific pointer during the simulation. This can
  !! handle 3 dimension and can preserve the contents of the array.
  !! IN:
  !!    p     :: original pointer
  !!    n1    :: size of the first dimension
  !!    n2    :: size of the second dimension
  !!    n3    :: size of the third dimension
  !!    o_preserve :: (optional) boolean indicating whether or not to
  !!                  preserve the contents of the array
  !!    o_aligned :: (optional) Use the FFTW library to allocate 16-byte
  !!                 aligned memory necessary to SIMD operations
  !
  function safemem_realloc_3d_real_int32(p, n1, n2, n3, o_preserve,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:,:) :: safemem_realloc_3d_real_int32
    _REAL_, pointer, dimension(:,:,:):: p
    logical, optional, intent(in) :: o_preserve,o_aligned
    integer, intent(in) :: n1, n2, n3
    safemem_realloc_3d_real_int32 => safemem_realloc_3d_real(p, &
         int(n1, kind=int64), int(n2, kind=int64), int(n3, kind=int64), &
         o_preserve,o_aligned)
  end function safemem_realloc_3d_real_int32

  function safemem_realloc_3d_real(p, n1, n2, n3, o_preserve,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:,:) :: safemem_realloc_3d_real
    _REAL_, pointer, dimension(:,:,:):: p
    logical, optional, intent(in) :: o_preserve,o_aligned
    integer(kind=int64), intent(in) :: n1, n2, n3
    logical :: preserve,aligned
    integer(kind=int64) :: nold1, nold2, nold3
    integer :: err
    integer(kind=int64) :: i2, i3
    type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
    write(6,*) "safemem_realloc_3d_real", associated(p),n1,n2,n3,preserve
    call flush(6)
#endif /*RISM_DEBUG*/
    aligned = .false.
    if (present(o_aligned)) aligned = o_aligned

    preserve = .true.
    if (present(o_preserve)) preserve = o_preserve

    if (.not.preserve) then
       if (safemem_dealloc(p,aligned) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 3D REAL")
    end if
    if (aligned) then
#ifndef CUDA    
       err = 0
       cptr = fftw_alloc_real(int(n1*n2*n3, C_SIZE_T))
       if (.not.c_associated(cptr)) then
          err = 1
       else
          call c_f_pointer(cptr,safemem_realloc_3d_real,[n1,n2,n3])
       end if
#else
        call rism_report_error('3: Cannot allocate aligned memory when using CUDA. Check you Fortran calls to safemem.')
#endif /*CUDA*/
    else
       allocate(safemem_realloc_3d_real(1:n1,1:n2,1:n3), STAT=err)
    end if
    if (err /= 0)&
         call rism_report_error("(a,1x,i4)", "reallocation_pointer allocate error: 3D REAL: ",err)
    call memadd_r(totalmem,arraySize(safemem_realloc_3d_real))

    if (associated(p) .and. preserve .and. n1*n2*n3 > 0) then
       safemem_realloc_3d_real = 0d0
       nold1 = min(ubound(p,1), n1)
       nold2 = min(ubound(p,2), n2)
       nold3 = min(ubound(p,3), n3)
       do i3=1,nold3
          do i2=1,nold2
             call dcopy(nold1, p(1, i2, i3), 1, &
                  safemem_realloc_3d_real(1, i2, i3),1)
          end do
       end do
       if (safemem_dealloc(p,aligned) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 3D REAL")
    end if
  end function SAFEMEM_REALLOC_3D_REAL

  !
  !! reallocate memory of a specific pointer during the simulation. This can
  !! handle 4 dimension and can preserve the contents of the array.
  !! IN:
  !!    p     :: original pointer
  !!    n1    :: size of the first dimension
  !!    n2    :: size of the second dimension
  !!    n3    :: size of the third dimension
  !!    n4    :: size of the fourth dimension
  !!    o_preserve :: (optional) boolean indicating whether or not to
  !!                  preserve the contents of the array
  !!    o_center :: (optional) when preserving the contents of the
  !!                array, center the data in the array instead of
  !!                preserving the value at each index
  !!    o_aligned :: (optional) Use the FFTW library to allocate 16-byte
  !!                 aligned memory necessary to SIMD operations
  !
  function safemem_realloc_4d_real_int32(p, n1, n2, n3, n4, o_preserve, o_center,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:,:,:) :: safemem_realloc_4d_real_int32
    _REAL_, pointer, dimension(:,:,:,:):: p
    logical,optional, intent(in) :: o_preserve, o_center,o_aligned
    integer, intent(in) :: n1, n2, n3, n4
    safemem_realloc_4d_real_int32 => safemem_realloc_4d_real(p, &
         int(n1,kind=int64), int(n2,kind=int64), int(n3,kind=int64), int(n4,kind=int64),&
         o_preserve, o_center,o_aligned)
  end function safemem_realloc_4d_real_int32
  
  function safemem_realloc_4d_real(p, n1, n2, n3, n4, o_preserve, o_center,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:,:,:) :: safemem_realloc_4d_real
    _REAL_, pointer, dimension(:,:,:,:):: p
    logical,optional, intent(in) :: o_preserve, o_center,o_aligned
    integer(kind=int64), intent(in) :: n1, n2, n3, n4
    logical :: preserve, center,aligned
    integer(kind=int64) :: nold1, nold2, nold3, nold4, offset(4)
    integer :: err
    integer(kind=int64) :: i2,i3,i4
    type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
    write(0,*) "safemem_realloc_4d_real", associated(p),n1,n2,n3,n4,preserve
    call flush(6)
#endif /*RISM_DEBUG*/
    aligned = .false.
    if (present(o_aligned)) aligned = o_aligned

    preserve = .true.
    if (present(o_preserve)) preserve=o_preserve

    center = .false.
    if (present(o_center)) center=o_center

    if (.not.preserve) then
       if (safemem_dealloc(p,aligned) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 4D REAL")
    end if
    if (aligned) then
#ifndef CUDA
       err = 0
       cptr = fftw_alloc_real(int(n1*n2*n3*n4, C_SIZE_T))
       if (.not.c_associated(cptr)) then
          err = 1
       else
          call c_f_pointer(cptr,safemem_realloc_4d_real,[n1,n2,n3,n4])
       end if
#else
       call rism_report_error('4: Cannot allocate aligned memory when using CUDA. Check you Fortran calls to safemem.')
#endif /*CUDA*/
    else
       allocate(safemem_realloc_4d_real(n1,n2,n3,n4), STAT=err)
    end if
    if (err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 4D REAL")
    call memadd_r(totalmem,arraySize(safemem_realloc_4d_real))

    if (associated(p) .and. preserve .and. n1*n2*n3*n4 > 0) then
       safemem_realloc_4d_real = 0d0
       nold1 = min(ubound(p,1), n1)
       nold2 = min(ubound(p,2), n2)
       nold3 = min(ubound(p,3), n3)
       nold4 = min(ubound(p,4), n4)
       if (center) then
          offset(1) = (ubound(p,1)- n1)/2
          offset(2) = (ubound(p,2)- n2)/2
          offset(3) = (ubound(p,3)- n3)/2
          offset(4) = (ubound(p,4)- n4)/2
       else
          offset=0
       end if
       do i4=1,nold4
          do i3=1,nold3
             do i2=1,nold2
                call dcopy(nold1, &
                     p(max(offset(1),0)+1, &
                     max(offset(2),0)+i2, &
                     max(offset(3),0)+i3, &
                     max(offset(4),0)+i4),1, &
                     safemem_realloc_4d_real(&
                     abs(min(offset(1),0))+1,&
                     abs(min(offset(2),0))+i2,&
                     abs(min(offset(3),0))+i3,&
                     abs(min(offset(4),0))+i4),1)
             end do
          end do
       end do
       if (safemem_dealloc(p,aligned) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 4D REAL")
    end if
  end function SAFEMEM_REALLOC_4D_REAL

  !
  !! reallocate memory of a specific pointer during the simulation. This can
  !! handle 5 dimensions. The original contents of the array are 3D polynomial
  !! interpolated to fill the new array.  That is, only the first three
  !! dimensions are interpolated
  !! IN:
  !!    p  :: original pointer
  !!    n1 :: size of the first dimension
  !!    n2 :: size of the second dimension
  !!    n3 :: size of the third dimension
  !!    n4 :: size of the fourth dimension
  !!    n5 :: size of the fifth dimension
  !!    xi1   :: input linear values of the first dimension grid points
  !!    xi2   :: input linear values of the second dimension grid points
  !!    xi2   :: input linear values of the third dimension grid points
  !!    xo1   :: input linear values of the first dimension grid points
  !!    xo2   :: input linear values of the second dimension grid points
  !!    xo2   :: input linear values of the third dimension grid points
  !!    o_aligned :: (optional) Use the FFTW library to allocate 16-byte
  !!                 aligned memory necessary to SIMD operations
  !
  function safemem_realloc_5d_real_interp_int32(p, n1, n2, n3, n4, n5,&
       xi1,xi2,xi3,xo1,xo2,xo3,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:,:,:,:) :: safemem_realloc_5d_real_interp_int32
    _REAL_, pointer, dimension(:,:,:,:,:):: p
    _REAL_, dimension(:), intent(in) :: xi1,xi2,xi3,xo1,xo2,xo3
    logical ,optional, intent(in) :: o_aligned
    integer, intent(in) :: n1, n2, n3, n4, n5
    safemem_realloc_5d_real_interp_int32 => safemem_realloc_5d_real_interp(p, &
         int(n1,kind=int64), int(n2,kind=int64), int(n3,kind=int64), int(n4,kind=int64), int(n5,kind=int64), &
         xi1,xi2,xi3,xo1,xo2,xo3,o_aligned)
  end function safemem_realloc_5d_real_interp_int32
  
  function safemem_realloc_5d_real_interp(p, n1, n2, n3, n4, n5,&
       xi1,xi2,xi3,xo1,xo2,xo3,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:,:,:,:) :: safemem_realloc_5d_real_interp
    _REAL_, pointer, dimension(:,:,:,:,:):: p
    _REAL_, dimension(:), intent(in) :: xi1,xi2,xi3,xo1,xo2,xo3
    logical ,optional, intent(in) :: o_aligned
    integer(kind=int64), intent(in) :: n1, n2, n3, n4, n5
    logical :: aligned
    _REAL_ :: dy
    integer(kind=int64) :: nold1, nold2,nold3, nold4,nold5,order =3
    integer :: err    
    integer(kind=int64) :: i1,i2,i3, i4,i5,index(3), id
    _REAL_ :: x(0:1,0:1,0:1),r,s,t
    type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
    write(6,*) "safemem_realloc_5d_real_interp", associated(p),n1,n2,n3,n4,n5
    call flush(6)
#endif /*RISM_DEBUG*/
    aligned = .false.
    if (present(o_aligned)) aligned = o_aligned

    nold1 = min(Ubound(p,1), n1)
    nold2 = min(ubound(p,2), n2)
    nold3 = min(ubound(p,3), n3)
    nold4 = min(ubound(p,4), n4)
    nold5 = min(ubound(p,5), n5)
    if (aligned) then
#ifndef CUDA
       err = 0
       cptr = fftw_alloc_real(int(n1*n2*n3*n4*n5, C_SIZE_T))
       if (.not.c_associated(cptr)) then
          err = 1
       else
          call c_f_pointer(cptr,safemem_realloc_5d_real_interp,[n1,n2,n3,n4,n5])
       end if
#else
       call rism_report_error('5: Cannot allocate aligned memory when using CUDA. Check you Fortran calls to safemem.')
#endif /*CUDA*/
    else
       allocate(safemem_realloc_5d_real_interp(n1,n2,n3,n4,n5), STAT=err)
    end if
    if (err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 5D REAL Interpolate")
    call memadd_r(totalmem,arraySize(safemem_realloc_5d_real_interp))
    if (.not. associated(p)) return
    nold1 = min(Ubound(p,1), n1)
    nold2 = min(ubound(p,2), n2)
    nold3 = min(ubound(p,3), n3)
    nold4 = min(ubound(p,4), n4)
    nold5 = min(ubound(p,5), n5)

    do i5 = 1 , n5
       do i4 = 1 , n4
          do i3 = 1 , n3
             do i2 = 1 , n2
                do i1 = 1 , n1
                   !find nearest original point
                   index = nearestindex(xi1,xi2,xi3,&
                        ubound(p,1),ubound(p,2),ubound(p,3),&
                        xo1(i1),xo2(i2),xo3(i3))
                   if (index(1) == 1 .or. index(2) == 1 .or. index(3) == 1 .or.&
                        index(1) == ubound(p,1) .or. index(2) == ubound(p,2) .or. index(3) == ubound(p,3)) then
                      safemem_realloc_5d_real_interp(i1,i2,i3,i4,i5) = p(index(1),index(2),index(3),i4,i5)
                   else
                      safemem_realloc_5d_real_interp(i1,i2,i3,i4,i5) = p(index(1),index(2),index(3),i4,i5)
                   end if
                end do
             end do
          end do
       end do
    end do
    if (safemem_dealloc(p,aligned) /= 0)&
         call rism_report_error("reallocation_pointer deallocate error: 5D REAL INTERPOLATE")
  end function safemem_realloc_5d_real_interp


  !> reallocate memory of a specific pointer during the simulation. This can
  !! handle 5 dimension and can preserve the contents of the array.
  !! IN:
  !!    p     :: original pointer
  !!    n1    :: size of the first dimension
  !!    n2    :: size of the second dimension
  !!    n3    :: size of the third dimension
  !!    n4    :: size of the fourth dimension
  !!    n5    :: size of the fifth dimension
  !!    o_preserve :: boolean indicating whether or not to preserve the contents of
  !!            the array (optional)
  !!    o_center :: when preserving the contents of the array, center the
  !!              data in the array instead of preserving the value at each
  !!              index
  !!    o_aligned :: (optional) Use the FFTW library to allocate 16-byte
  !!                 aligned memory necessary to SIMD operations
  function safemem_realloc_5d_real_int32(p, n1, n2, n3, n4, n5, o_preserve, o_center,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:,:,:,:) :: safemem_realloc_5d_real_int32
    _REAL_, pointer, dimension(:,:,:,:,:):: p
    logical, optional, intent(in) :: o_preserve, o_center,o_aligned
    integer, intent(in) :: n1, n2, n3, n4, n5
    safemem_realloc_5d_real_int32 => safemem_realloc_5d_real(p, &
         int(n1, kind=int64), int(n2, kind=int64), int(n3, kind=int64), int(n4, kind=int64), int(n5, kind=int64), &
         o_preserve, o_center,o_aligned)
  end function safemem_realloc_5d_real_int32

  function safemem_realloc_5d_real(p, n1, n2, n3, n4, n5, o_preserve, o_center,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:,:,:,:) :: safemem_realloc_5d_real
    _REAL_, pointer, dimension(:,:,:,:,:):: p
    logical, optional, intent(in) :: o_preserve, o_center,o_aligned
    integer(kind=int64), intent(in) :: n1, n2, n3, n4, n5
    logical :: preserve,center,aligned
    integer(kind=int64) :: nold1, nold2, nold3, nold4, nold5, offset(5)
    integer :: err
    integer(kind=int64) :: i1,i2,i3,i4,i5
    type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
    write(6,*) "safemem_realloc_5d_real", associated(p),n1,n2,n3,n4,n5,preserve
    call flush(6)
#endif /*RISM_DEBUG*/
    aligned = .false.
    if (present(o_aligned)) aligned = o_aligned

    preserve = .true.
    if (present(o_preserve)) preserve = o_preserve

    center=.false.
    if (present(o_center)) center = o_center

    if (.not.preserve) then
       if (safemem_dealloc(p,aligned) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 5D REAL")
    end if
    if (aligned) then
#ifndef CUDA
       err = 0
       cptr = fftw_alloc_real(int(n1*n2*n3*n4*n5, C_SIZE_T))
       if (.not.c_associated(cptr)) then
          err = 1
       else
          call c_f_pointer(cptr,safemem_realloc_5d_real,[n1,n2,n3,n4,n5])
       end if
#else
       call rism_report_error('6: Cannot allocate aligned memory when using CUDA. Check you Fortran calls to safemem.')
#endif /*CUDA*/
    else
       allocate(safemem_realloc_5d_real(n1,n2,n3,n4,n5), STAT=err)
    end if
    if (err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 5D REAL")
    call memadd_r(totalmem,arraySize(safemem_realloc_5d_real))
    if (associated(p) .and. preserve .and. n1*n2*n3*n4*n5 > 0) then
       safemem_realloc_5d_real = 0d0
       nold1 = min(ubound(p,1), n1)
       nold2 = min(ubound(p,2), n2)
       nold3 = min(ubound(p,3), n3)
       nold4 = min(ubound(p,4), n4)
       nold5 = min(ubound(p,5), n5)
       if (center) then
          offset(1) = (ubound(p,1)- n1)/2
          offset(2) = (ubound(p,2)- n2)/2
          offset(3) = (ubound(p,3)- n3)/2
          offset(4) = (ubound(p,4)- n4)/2
          offset(5) = (ubound(p,5)- n5)/2
       else
          offset = 0
       end if
       do i5=1,nold5
          do i4=1,nold4
             do i3=1,nold3
                do i2=1,nold2
                   call dcopy(nold1, &
                        p(max(offset(1),0)+1, &
                        max(offset(2),0)+i2, &
                        max(offset(3),0)+i3, &
                        max(offset(4),0)+i4, &
                        max(offset(5),0)+i5),1,&
                        safemem_realloc_5d_real(&
                        abs(min(offset(1),0))+1,&
                        abs(min(offset(2),0))+i2,&
                        abs(min(offset(3),0))+i3,&
                        abs(min(offset(4),0))+i4,&
                        abs(min(offset(5),0))+i5),1)
                end do
             end do
          end do
       end do
       if (safemem_dealloc(p,aligned) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 5D REAL")
    endif
  end function safemem_realloc_5d_real


  !> reallocate memory of a specific pointer during the simulation. This can
  !! handle 1 dimension and preserves the contents of the array.
  !! IN:
  !!    p     :: original pointer
  !!    n1    :: size of the first dimension
  !!    o_preserve :: boolean indicating whether or not to preserve the contents of
  !!            the array (optional)
  function safemem_realloc_1d_integer_int32(p, n1, o_preserve)
    implicit none
    integer, pointer, dimension(:) :: safemem_realloc_1d_integer_int32
    integer, pointer, dimension(:) :: p
    logical, optional, intent(in) :: o_preserve
    integer, intent(in) :: n1
    safemem_realloc_1d_integer_int32 => safemem_realloc_1d_integer(p, int(n1, kind=int64), o_preserve)
  end function safemem_realloc_1d_integer_int32
  
  function safemem_realloc_1d_integer(p, n1, o_preserve)
    implicit none
    integer, pointer, dimension(:) :: safemem_realloc_1d_integer
    integer, pointer, dimension(:) :: p
    logical, optional, intent(in) :: o_preserve
    integer(kind=int64), intent(in) :: n1
    logical :: preserve
    integer(kind=int64) :: nold1
    integer :: err    
    integer(kind=int64) :: i1

#ifdef RISM_DEBUG
    write(6,*) "safemem_realloc_1d_integer", associated(p),n1,preserve
    call flush(6)
#endif /*RISM_DEBUG*/
    preserve = .true.
    if (present(o_preserve)) preserve = o_preserve

    if (.not.preserve) then
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 1D integer")
    end if
    allocate(safemem_realloc_1d_integer(1:n1), STAT=err)
    call memadd_i(totalmem,arraySize(safemem_realloc_1d_integer))
    if (err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 1D integer")
    if (associated(p) .and. preserve .and. n1 > 0) then
       nold1 = min(ubound(p,1), n1)
       !explicit loop to prevent stack overflow with intel compilers
       do i1=1,nold1
          safemem_realloc_1d_integer(i1) = p(i1)
       end do
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 1D integer")
    end if
  end function SAFEMEM_REALLOC_1D_integer

  function safemem_realloc_1d_int64(p, n1, o_preserve)
    implicit none
    integer(kind=int64), pointer, dimension(:) :: safemem_realloc_1d_int64
    integer(kind=int64), pointer, dimension(:) :: p
    logical, optional, intent(in) :: o_preserve
    integer(kind=int64), intent(in) :: n1
    logical :: preserve
    integer(kind=int64) :: nold1
    integer :: err    
    integer(kind=int64) :: i1

#ifdef RISM_DEBUG
    write(6,*) "safemem_realloc_1d_int64", associated(p),n1,preserve
    call flush(6)
#endif /*RISM_DEBUG*/
    preserve = .true.
    if (present(o_preserve)) preserve = o_preserve

    if (.not.preserve) then
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 1D int64")
    end if
    allocate(safemem_realloc_1d_int64(1:n1), STAT=err)
    call memadd_i(totalmem,arraySize(safemem_realloc_1d_int64))
    if (err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 1D int64")
    if (associated(p) .and. preserve .and. n1 > 0) then
       nold1 = min(ubound(p,1), n1)
       !explicit loop to prevent stack overflow with intel compilers
       do i1=1,nold1
          safemem_realloc_1d_int64(i1) = p(i1)
       end do
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 1D int64")
    end if
  end function safemem_realloc_1d_int64

  !
  !! reallocate memory of a specific pointer during the simulation. This can
  !! handle 2 dimension and preserves the contents of the array.
  !! IN:
  !!    p     :: original pointer
  !!    n1    :: size of the first dimension
  !!    n2    :: size of the second dimension
  !!    o_preserve :: boolean indicating whether or not to preserve the contents of
  !!            the array (optional)
  !
  function safemem_realloc_2d_integer_int32(p, n1, n2,o_preserve)
    implicit none
    integer, pointer, dimension(:,:) :: safemem_realloc_2d_integer_int32
    integer, pointer, dimension(:,:) :: p
    logical, optional, intent(in) :: o_preserve
    integer, intent(in) :: n1,n2
    logical :: preserve
    safemem_realloc_2d_integer_int32 => safemem_realloc_2d_integer(p, int(n1, kind=int64), int(n2, kind=int64), o_preserve)
  end function safemem_realloc_2d_integer_int32

  function safemem_realloc_2d_integer(p, n1, n2,o_preserve)
    implicit none
    integer, pointer, dimension(:,:) :: safemem_realloc_2d_integer
    integer, pointer, dimension(:,:) :: p
    logical, optional, intent(in) :: o_preserve
    integer(kind=int64), intent(in) :: n1,n2
    logical :: preserve
    integer(kind=int64) :: nold1,nold2
    integer :: err    
    integer(kind=int64) :: i1, i2

#ifdef RISM_DEBUG
    write(6,*) "safemem_realloc_2d_integer", associated(p),n1,n2,preserve
    call flush(6)
#endif /*RISM_DEBUG*/
    preserve=.true.
    if (present(o_preserve)) preserve = o_preserve

    if (.not.preserve) then
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 2D integer")
    end if
    allocate(safemem_realloc_2d_integer(1:n1,1:n2), STAT=err)
    call memadd_i(totalmem,arraySize(safemem_realloc_2d_integer))
    if (err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 2D integer")
    if (associated(p) .and. preserve .and. n1*n2 > 0) then
       nold1 = min(ubound(p,1), n1)
       nold2 = min(ubound(p,2), n2)
       !explicit loop to prevent stack overflow with intel compilers
       do i1=1,nold1
          do i2=1,nold2
             safemem_realloc_2d_integer(i1,i2) = p(i1,i2)
          end do
       end do
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 2D integer")
    end if
  end function SAFEMEM_REALLOC_2D_integer

  !
  !! reallocate memory of a specific pointer during the simulation. This can
  !! handle 1 dimension and preserves the contents of the array.
  !! IN:
  !!    p     :: original pointer
  !!    length :: character length for each element
  !!    n1    :: size of the first dimension
  !!    o_preserve :: boolean indicating whether or not to preserve the contents of
  !!            the array (optional)
  !
  function safemem_realloc_1d_character_int32(p, length, n1,o_preserve)
    implicit none
    integer, intent(in) :: length
    character(len=length), pointer, dimension(:) :: safemem_realloc_1d_character_int32
    character(len=length), pointer, dimension(:) :: p
    logical, optional, intent(in) :: o_preserve
    integer, intent(in) :: n1
    safemem_realloc_1d_character_int32 => safemem_realloc_1d_character(p, length, int(n1, kind=int64), o_preserve)
  end function safemem_realloc_1d_character_int32
  
  function safemem_realloc_1d_character(p, length, n1,o_preserve)
    implicit none
    integer, intent(in) :: length
    character(len=length), pointer, dimension(:) :: safemem_realloc_1d_character
    character(len=length), pointer, dimension(:) :: p
    logical, optional, intent(in) :: o_preserve
    integer(kind=int64), intent(in) :: n1
    logical :: preserve
    integer(kind=int64) :: nold1
    integer :: err
    integer(kind=int64) :: i1

#ifdef RISM_DEBUG
    write(6,*) "safemem_realloc_1d_character", associated(p),n1,o_preserve
    call flush(6)
#endif /*RISM_DEBUG*/
    preserve = .true.
    if (present(o_preserve)) preserve = o_preserve

    if (.not.preserve) then
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 1D CHARACTER")
    end if
    allocate(safemem_realloc_1d_character(1:n1), STAT=err)
    call memadd_c(totalmem,arraySize(safemem_realloc_1d_character))
    if (err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 1D CHARACTER")
    if (associated(p) .and. preserve .and. n1 > 0) then
       nold1 = min(ubound(p,1), n1)
       !explicit loop to prevent stack overflow with intel compilers
       do i1=1,nold1
          safemem_realloc_1d_character(i1) = trim(p(i1))
       end do
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 1D CHARACTER")
    end if
  end function SAFEMEM_REALLOC_1D_CHARACTER

  !
  !! reallocate memory of a specific pointer during the simulation. This can
  !! handle 2 dimensions and preserves the contents of the array.
  !! IN:
  !!    p     :: original pointer
  !!    length :: character length for each element
  !!    n1    :: size of the first dimension
  !!    n2    :: size of the second dimension
  !!    o_preserve :: boolean indicating whether or not to o_preserve the contents of
  !!            the array (optional)
  !
  function safemem_realloc_2d_character_int32(p, length, n1, n2,o_preserve)
    implicit none
    integer, intent(in) :: length
    character(len=length), pointer, dimension(:,:) :: safemem_realloc_2d_character_int32
    character(len=length), pointer, dimension(:,:) :: p
    logical, optional, intent(in) :: o_preserve
    integer, intent(in) :: n1, n2
    safemem_realloc_2d_character_int32 => safemem_realloc_2d_character(p, length, int(n1, kind=int64), int(n2, kind=int64), o_preserve)
  end function safemem_realloc_2d_character_int32
  function safemem_realloc_2d_character(p, length, n1, n2,o_preserve)
    implicit none
    integer, intent(in) :: length
    character(len=length), pointer, dimension(:,:) :: safemem_realloc_2d_character
    character(len=length), pointer, dimension(:,:) :: p
    logical, optional, intent(in) :: o_preserve
    integer(kind=int64), intent(in) :: n1, n2
    logical :: preserve
    integer(kind=int64) :: nold1, nold2
    integer :: err
    integer(kind=int64) :: i1, i2

#ifdef RISM_DEBUG
    write(6,*) "safemem_realloc_2d_character", associated(p),n1,n2,o_preserve
    call flush(6)
#endif /*RISM_DEBUG*/
    preserve = .true.
    if (present(o_preserve)) preserve = o_preserve

    if (.not.preserve) then
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 2D CHARACTER")
    end if
    allocate(safemem_realloc_2d_character(1:n1,1:n2), STAT=err)
    call memadd_c(totalmem,arraySize(safemem_realloc_2d_character))
    if (err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 1D CHARACTER")
    if (associated(p) .and. preserve .and. n1*n2 > 0) then
       nold1 = min(ubound(p,1), n1)
       nold2 = min(ubound(p,2), n2)
       !explicit loop to prevent stack overflow with intel compilers
       do i1=1,nold1
          do i2=1,nold2
             safemem_realloc_2d_character(i1,i2) = trim(p(i1,i2))
          end do
       end do
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 2D CHARACTER")
    end if
  end function SAFEMEM_REALLOC_2D_CHARACTER


  !! reallocate memory of a specific pointer during the simulation. This can
  !! handle 1 dimension and preserves the contents of the array.
  !! IN:
  !!    p     :: original pointer
  !!    n1    :: size of the first dimension
  !!    o_preserve :: boolean indicating whether or not to preserve the contents of
  !!            the array (optional)

  function safemem_realloc_1d_logical_int32(p, n1,o_preserve)
    implicit none
    logical, pointer, dimension(:) :: safemem_realloc_1d_logical_int32
    logical, pointer, dimension(:) :: p
    logical, optional, intent(in) :: o_preserve
    integer, intent(in) :: n1
    safemem_realloc_1d_logical_int32 => safemem_realloc_1d_logical(p, int(n1, kind=int64) ,o_preserve)
  end function safemem_realloc_1d_logical_int32
  
  function safemem_realloc_1d_logical(p, n1,o_preserve)
    implicit none
    logical, pointer, dimension(:) :: safemem_realloc_1d_logical
    logical, pointer, dimension(:) :: p
    logical, optional, intent(in) :: o_preserve
    integer(kind=int64), intent(in) :: n1
    logical :: preserve
    integer(kind=int64) :: nold1
    integer :: err
    integer(kind=int64) :: i1

#ifdef RISM_DEBUG
    write(6,*) "safemem_realloc_1d_logical", associated(p),n1,preserve
    call flush(6)
#endif /*RISM_DEBUG*/
    preserve = .true.
    if (present(o_preserve)) preserve = o_preserve

    if (.not.preserve) then
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 1D LOGICAL")
    end if
    allocate(safemem_realloc_1d_logical(1:n1), STAT=err)
    call memadd_l(totalmem,arraySize(safemem_realloc_1d_logical))
    if (err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 1D LOGICAL")
    if (associated(p) .and. preserve .and. n1 > 0) then
       nold1 = min(ubound(p,1), n1)
       !explicit loop to prevent stack overflow with intel compilers
       do i1=1,nold1
          safemem_realloc_1d_logical(i1) = p(i1)
       end do
       if (safemem_dealloc(p) /= 0)&
            call rism_report_error("reallocation_pointer deallocate error: 1D LOGICAL")
    end if
  end function safemem_realloc_1d_logical

  !!
  !! checks if the pointer is allocated before deallocating.  Returns any error code
  !! IN:
  !!    p  :: original pointer
  !!    o_aligned :: (optional) Use the FFTW library to free 16-byte
  !!                 aligned memory.  It is not possible to tell how
  !!                 memory was originally allocated
  !!
  function safemem_dealloc_pointer_1d_real(p,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:) :: p
    logical ,optional, intent(in) :: o_aligned
    logical :: aligned
    type(c_ptr) :: cptr
    integer :: safemem_dealloc_pointer_1d_real,temp
    temp=0

    aligned = .false.
    if (present(o_aligned)) aligned = o_aligned

    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_1d_real)
    if (associated(p)) then
       call memadd_r(totalmem,-arraySize(p))
       if (aligned) then
#ifndef CUDA
          cptr = c_loc(p(1))
          call fftw_free(cptr)
          nullify(p)
          temp=0
#else
          call rism_report_error('How do we got here???')
#endif /*CUDA*/
       else
          deallocate(p,STAT=temp)
       end if
    end if
    safemem_dealloc_pointer_1d_real=temp
  end function safemem_dealloc_pointer_1d_real

  !!
  !! checks if the pointer is allocated before deallocating.  Returns any error code
  !! IN:
  !!    p  :: original pointer
  !!    o_aligned :: (optional) Use the FFTW library to free 16-byte
  !!                 aligned memory.  It is not possible to tell how
  !!                 memory was originally allocated
  !!
  function safemem_dealloc_pointer_2d_real(p,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:) :: p
    logical ,optional, intent(in) :: o_aligned
    logical :: aligned
    type(c_ptr) :: cptr
    integer :: safemem_dealloc_pointer_2d_real,temp
    temp=0

    aligned = .false.
    if (present(o_aligned)) aligned = o_aligned

    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_2d_real)
    if (associated(p)) then
       call memadd_r(totalmem,-arraySize(p))
       if (aligned) then
#ifndef CUDA
          cptr = c_loc(p(1,1))
          call fftw_free(cptr)
          nullify(p)
#else
          call rism_report_error('How do we got here???')
#endif /*CUDA*/
       else
          deallocate(p,STAT=temp)
       end if
    end if
    safemem_dealloc_pointer_2d_real=temp
  end function safemem_dealloc_pointer_2d_real

  !!
  !! checks if the pointer is allocated before deallocating.  Returns any error code
  !! IN:
  !!    p  :: original pointer
  !!    o_aligned :: (optional) Use the FFTW library to free 16-byte
  !!                 aligned memory.  It is not possible to tell how
  !!                 memory was originally allocated
  !!
  function safemem_dealloc_pointer_3d_real(p,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:,:) :: p
    logical ,optional, intent(in) :: o_aligned
    logical :: aligned
    type(c_ptr) :: cptr
    integer :: safemem_dealloc_pointer_3d_real,temp
    temp=0

    aligned = .false.
    if (present(o_aligned)) aligned = o_aligned

    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_3d_real)
    if (associated(p)) then
       call memadd_r(totalmem,-arraySize(p))
       if (aligned) then
#ifndef CUDA
          cptr = c_loc(p(1,1,1))
          call fftw_free(cptr)
          nullify(p)
#else
          call rism_report_error('How do we got here???')
#endif /*CUDA*/
       else
          deallocate(p,STAT=temp)
       end if
    end if
    safemem_dealloc_pointer_3d_real=temp
  end function safemem_dealloc_pointer_3d_real

  !!
  !! checks if the pointer is allocated before deallocating.  Returns any error code
  !! IN:
  !!    p  :: original pointer
  !!    o_aligned :: (optional) Use the FFTW library to free 16-byte
  !!                 aligned memory.  It is not possible to tell how
  !!                 memory was originally allocated
  !!
  function safemem_dealloc_pointer_4d_real(p,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:,:,:) :: p
    logical ,optional, intent(in) :: o_aligned
    logical :: aligned
    type(c_ptr) :: cptr
    integer :: safemem_dealloc_pointer_4d_real,temp
    temp=0

    aligned = .false.
    if (present(o_aligned)) aligned = o_aligned

    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_4d_real)
    if (associated(p)) then
       call memadd_r(totalmem,-arraySize(p))
       if (aligned) then
#ifndef CUDA
          cptr = c_loc(p(1,1,1,1))
          call fftw_free(cptr)
          nullify(p)
#else
          call rism_report_error('How do we got here???')
#endif /*CUDA*/
       else
          deallocate(p,STAT=temp)
       end if
    end if
    safemem_dealloc_pointer_4d_real=temp
  end function safemem_dealloc_pointer_4d_real

  !!
  !! checks if the pointer is allocated before deallocating.  Returns any error code
  !! IN:
  !!    p  :: original pointer
  !!    o_aligned :: (optional) Use the FFTW library to free 16-byte
  !!                 aligned memory.  It is not possible to tell how
  !!                 memory was originally allocated
  !!
  function safemem_dealloc_pointer_5d_real(p,o_aligned)
    implicit none
    _REAL_, pointer, dimension(:,:,:,:,:) :: p
    logical ,optional, intent(in) :: o_aligned
    logical :: aligned
    type(c_ptr) :: cptr
    integer :: safemem_dealloc_pointer_5d_real,temp
    temp=0

    aligned = .false.
    if (present(o_aligned)) aligned = o_aligned

    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_5d_real)
    if (associated(p)) then
       call memadd_r(totalmem,-arraySize(p))
       if (aligned) then
#ifndef CUDA
          cptr = c_loc(p(1,1,1,1,1))
          call fftw_free(cptr)
          nullify(p)
#else
          call rism_report_error('How do we got here???')
#endif /*CUDA*/
       else
          deallocate(p,STAT=temp)
       end if
    end if
    safemem_dealloc_pointer_5d_real=temp
  end function safemem_dealloc_pointer_5d_real

  
  !> Checks if the pointer is allocated before deallocating.  Returns
  !! any error code.
  function safemem_dealloc_pointer_1d_int(p)
    implicit none
    integer, pointer, dimension(:) :: p
    integer :: safemem_dealloc_pointer_1d_int,temp
    temp=0
    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_1d_int)
    if (associated(p)) then
       call memadd_i(totalmem,-arraySize(p))
       deallocate(p,STAT=temp)
    end if
    safemem_dealloc_pointer_1d_int=temp
  end function safemem_dealloc_pointer_1d_int
  
  function safemem_dealloc_pointer_1d_int64(p)
    implicit none
    integer(kind = int64), pointer, dimension(:) :: p
    integer :: safemem_dealloc_pointer_1d_int64, temp
    temp=0
    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_1d_int64)
    if (associated(p)) then
       call memadd_i(totalmem, -arraySize(p))
       deallocate(p,STAT=temp)
    end if
    safemem_dealloc_pointer_1d_int64=temp
  end function safemem_dealloc_pointer_1d_int64

  
  !> Checks if the pointer is allocated before deallocating.  Returns
  !! any error code.
  function safemem_dealloc_pointer_2d_int(p)
    implicit none
    integer, pointer, dimension(:,:) :: p
    integer :: safemem_dealloc_pointer_2d_int,temp
    temp=0
    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_2d_int)
    if (associated(p)) then
       call memadd_i(totalmem,-arraySize(p))
       deallocate(p,STAT=temp)
    end if
    safemem_dealloc_pointer_2d_int=temp
  end function safemem_dealloc_pointer_2d_int

  
  !> Checks if the pointer is allocated before deallocating.  Returns
  !! any error code.
  function safemem_dealloc_pointer_3d_int(p)
    implicit none
    integer, pointer, dimension(:,:,:) :: p
    integer :: safemem_dealloc_pointer_3d_int,temp
    temp=0
    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_3d_int)
    if (associated(p)) then
       call memadd_i(totalmem,-arraySize(p))
       deallocate(p,STAT=temp)
    end if
    safemem_dealloc_pointer_3d_int=temp
  end function safemem_dealloc_pointer_3d_int

  
  !> Checks if the pointer is allocated before deallocating.  Returns
  !! any error code.
  function safemem_dealloc_pointer_4d_int(p)
    implicit none
    integer, pointer, dimension(:,:,:,:) :: p
    integer :: safemem_dealloc_pointer_4d_int,temp
    temp=0
    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_4d_int)
    if (associated(p)) then
       call memadd_i(totalmem,-arraySize(p))
       deallocate(p,STAT=temp)
    end if
    safemem_dealloc_pointer_4d_int=temp
  end function safemem_dealloc_pointer_4d_int

  
  !> Checks if the pointer is allocated before deallocating.  Returns
  !! any error code.
  function safemem_dealloc_pointer_5d_int(p)
    implicit none
    integer, pointer, dimension(:,:,:,:,:) :: p
    integer :: safemem_dealloc_pointer_5d_int,temp
    temp=0
    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_5d_int)
    if (associated(p)) then
       call memadd_i(totalmem,-arraySize(p))
       deallocate(p,STAT=temp)
    end if
    safemem_dealloc_pointer_5d_int=temp
  end function safemem_dealloc_pointer_5d_int


  !> Checks if the pointer is allocated before deallocating.  Returns
  !! any error code.
  function safemem_dealloc_pointer_1d_character(p)
    implicit none
    character(len=*), pointer, dimension(:) :: p
    integer :: safemem_dealloc_pointer_1d_character,temp
    temp=0
    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_1d_character)
    if (associated(p)) then
       call memadd_c(totalmem,-arraySize(p))
       deallocate(p,STAT=temp)
    end if
    safemem_dealloc_pointer_1d_character = temp
  end function safemem_dealloc_pointer_1d_character

  
  !> Checks if the pointer is allocated before deallocating.  Returns
  !! any error code.
  function safemem_dealloc_pointer_2d_character(p)
    implicit none
    character(len=*), pointer, dimension(:,:) :: p
    integer :: safemem_dealloc_pointer_2d_character,temp
    temp=0
    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_2d_character)
    if (associated(p)) then
       call memadd_c(totalmem,-arraySize(p))
       deallocate(p,STAT=temp)
    end if
    safemem_dealloc_pointer_2d_character = temp
  end function safemem_dealloc_pointer_2d_character

  
  !> Check if the pointer is allocated before deallocating.  Returns any error code.
  function safemem_dealloc_pointer_1d_logical(p)
    implicit none
    logical, pointer, dimension(:) :: p
    integer :: safemem_dealloc_pointer_1d_logical,temp
    temp=0
    !work around for old gfortran bug
    !  if (associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_1d_logical)
    if (associated(p)) then
       call memadd_l(totalmem,-arraySize(p))
       deallocate(p,STAT=temp)
    end if
    safemem_dealloc_pointer_1d_logical = temp
  end function safemem_dealloc_pointer_1d_logical


  !> Returns an array of current and maximum memory allocated in bytes so
  !! far for a specific memTracker.  No distinction is made between different
  !! kinds of the same type:
  !! (1) : Current integer memory
  !! (2) : Current real memory
  !! (3) : Current logical memory
  !! (4) : Current char memory
  !! (5) : Current total memory
  !! (6) : Max integer memory
  !! (7) : Max real memory
  !! (8) : Max logical memory
  !! (9) : Max char memory
  !! (10): Max total memory
  !! IN:
  !!    this : memTracker object
  function memStatus_obj(this) result(mem)
    implicit none
    type(memTracker), intent(in) :: this
    integer(kind=int64) :: mem(10)
    mem = (/this%int, this%real, this%logical, this%char, this%total,&
         this%maxint, this%maxreal, this%maxlogical, this%maxchar, this%max/)
  end function memStatus_obj


  !> Returns an array of global current and maximum memory allocated in bytes so
  !! far.  No distinction is made between different kinds of the same type:
  !! (1) : Current integer memory
  !! (2) : Current real memory
  !! (3) : Current logical memory
  !! (4) : Current char memory
  !! (5) : Current logical memory
  !! (6) : Max integer memory
  !! (7) : Max real memory
  !! (8) : Max logical memory
  !! (9) : Max char memory
  !! (10): Max logical memory
  function memStatus_global() result(mem)
    implicit none
    integer(kind=int64) :: mem(10)
    mem = memStatus_obj(totalMem)
  end function memStatus_global

#define BITS_PER_BYTE 8

! If available, we use STORAGE_SIZE(), which is the 2008 standard.
! This returns the number of bits for *one* element. KIND() is not
! guaranteed to give the number of bytes per element, but before
! gfortran 4.6 it did and STORAGE_SIZE() was not available on these
! versions.
#if defined(__GFORTRAN__) && __GNUC__ == 4 && __GNUC_MINOR__ <= 6
  ! bits per element for gfortran < 4.6
#  define STORAGE_SIZE BITS_PER_BYTE*kind
#else
  ! bits per element for Fortran 2008 standard
#  define STORAGE_SIZE storage_size
#endif
  function arraySize_1d_real(array) result (nbytes)
    implicit none
    _REAL_, intent(in) :: array(:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE
  end function arraySize_1d_real
    
  function arraySize_2d_real(array) result (nbytes)
    implicit none
    _REAL_, intent(in) :: array(:,:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE
  end function arraySize_2d_real
    
  function arraySize_3d_real(array) result (nbytes)
    implicit none
    _REAL_, intent(in) :: array(:,:,:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE
  end function arraySize_3d_real
    
  function arraySize_4d_real(array) result (nbytes)
    implicit none
    _REAL_, intent(in) :: array(:,:,:,:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE
  end function arraySize_4d_real
    
  function arraySize_5d_real(array) result (nbytes)
    implicit none
    _REAL_, intent(in) :: array(:,:,:,:,:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE
  end function arraySize_5d_real
    
  function arraySize_1d_int(array) result (nbytes)
    implicit none
    integer, intent(in) :: array(:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE

  end function arraySize_1d_int

  function arraySize_1d_int64(array) result (nbytes)
    implicit none
    integer(kind=int64), intent(in) :: array(:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE

  end function arraySize_1d_int64
  
  function arraySize_2d_int(array) result (nbytes)
    implicit none
    integer, intent(in) :: array(:,:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE
  end function arraySize_2d_int
    
  function arraySize_3d_int(array) result (nbytes)
    implicit none
    integer, intent(in) :: array(:,:,:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE
  end function arraySize_3d_int
    
  function arraySize_4d_int(array) result (nbytes)
    implicit none
    integer, intent(in) :: array(:,:,:,:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE
  end function arraySize_4d_int
    
  function arraySize_5d_int(array) result (nbytes)
    implicit none
    integer, intent(in) :: array(:,:,:,:,:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE
  end function arraySize_5d_int

  function arraySize_1d_character(array) result (nbytes)
    implicit none
    character(len=*), intent(in) :: array(:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE

  end function arraySize_1d_character
    
  function arraySize_2d_character(array) result (nbytes)
    implicit none
    character(len=*), intent(in) :: array(:,:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE
  end function arraySize_2d_character
    
  function arraySize_1d_logical(array) result (nbytes)
    implicit none
    logical, intent(in) :: array(:)
    integer(kind=int64) :: nbytes
    nbytes = size(array,kind=int64)*STORAGE_SIZE(array)/BITS_PER_BYTE

  end function arraySize_1d_logical
      

  !> Adds the size in bytes of the array and increments the integer counter.  Use
  !! a negative number to subtract.
  !! IN:
  !!    this  : memtracker object
  !!    nbyte : number of bytes
  subroutine memadd_i(this,nbyte)
    implicit none
    type(memTracker),intent(inout) :: this
    integer(kind=int64), intent(in) :: nbyte
    this%int = this%int+nbyte
    this%total = this%total+nbyte
    if (this%int > this%maxint) this%maxint = this%int
    if (this%total > this%max) this%max = this%total
  end subroutine memadd_i


  !> Adds the size in bytes of the array and increments the real counter.  Use
  !! a negative number to subtract.
  !! IN:
  !!    this  : memtracker object
  !!    nbyte : number of bytes
  subroutine memadd_r(this,nbyte)
    implicit none
    type(memTracker),intent(inout) :: this
    integer(kind=int64), intent(in) :: nbyte
    this%real = this%real+nbyte
    this%total = this%total+nbyte
    if (this%real > this%maxreal) this%maxreal = this%real
    if (this%total > this%max) this%max = this%total
  end subroutine memadd_r


  !> Adds the size in bytes of the array and increments the character counter.  Use
  !! a negative number to subtract.
  !! IN:
  !!    this  : memtracker object
  !!    nbyte : number of bytes
  subroutine memadd_c(this,nbyte)
    implicit none
    type(memTracker),intent(inout) :: this
    integer(kind=int64), intent(in) :: nbyte
    this%char = this%char+nbyte
    this%total = this%total+nbyte
    if (this%char > this%maxchar) this%maxchar = this%char
    if (this%total > this%max) this%max = this%total
  end subroutine memadd_c

  !> Adds the size in bytes of the array and increments the logical counter.  Use
  !! a negative number to subtract.
  !! IN:
  !!    this  : memtracker object
  !!    nbyte : number of bytes
  subroutine memadd_l(this,nbyte)
    implicit none
    type(memTracker),intent(inout) :: this
    integer(kind=int64), intent(in) :: nbyte
    this%logical = this%logical+nbyte
    this%total = this%total+nbyte
    if (this%logical > this%maxlogical) this%maxlogical = this%logical
    if (this%total > this%max) this%max = this%total
  end subroutine memadd_l
end module safemem
