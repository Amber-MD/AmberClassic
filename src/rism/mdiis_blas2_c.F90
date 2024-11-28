!<compile=optimized>

#include "../include/dprec.fh"

!> Original MDIIS implementation by Andriy Kovalenko.
!! Parallelization by Sergey Gusarov.
!! Further modifications and objectification by Tyler Luchko.
!! New memory model and BLASification Tyler Luchko December, 2009
!! Memory swaps removed Tyler Luchko January, 2011
!! 
!! To use this implementation, the user allocates two large chucks of
!! contiguous memory; one for data and one for the residual that
!! should be minimized. Regardless of the number of dimensions of the
!! actual data, the each memory chuck should be (ndata,nvector) where
!! nvector is the number of data sets from previous iterations.  On
!! initialization, the first vector in each array is active and should
!! be updated.  After subsequent mdiis_advance calls, the active
!! vector will be updated and the vector number can be obtained from
!! mdiis_getVectorNumber.
!! 
!! The usage pattern allows MDIIS to avoid costly memory swaps.
!! Ideally, the memory should be allocated in this module and a
!! pointer to the active array returned to the user.  However, the
!! user data may be of any dimension and pointer rank remapping is not
!! available until Fortran 2003, which is not generally support
!! (especially this feature).
module mdiis_blas2_c
  use rism_report_c
  implicit none
  type mdiis_blas2
     private
     !istep  :: number of steps since start/restart
     !np     :: number of data points (first dimension length)
     !nvec    :: max number of vectors
     integer ::  istep = 0, np=0, nvec=0
     !overlap :: overlap matrix
     _REAL_,dimension(:,:),pointer :: overlap=>NULL()
     !ratio of the current residual to the minimum residual found that causes a restart
     _REAL_ ::  restart
     integer :: mpirank=0, mpisize=1, mpicomm=0
     !vecMap :: maps vectors order onto ri and xi.  E.g. vecMap(1) gives the index of xi or 
     !          ri that is most recently updated
     integer, pointer :: vecMap(:)=>NULL()
     !delta :: coefficient for residual gradient.  Should be between 0 and 1.
     !tol   :: target residual
     _REAL_ :: delta, tol
     !pointers to working data
     !xi    :: array of vector data.  np X nvec
     !ri    :: array of residual data.  np X nvec
     _REAL_,pointer :: xi(:,:)=>NULL(), ri(:,:)=>NULL()

  end type mdiis_blas2

  interface mdiis_blas2_new
     module procedure mdiis_blas2_new_serial, mdiis_blas2_new_mpi
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Create new MDIIS object.  Provide convergence parameters and
!!!working memory.  Note that the working memory must me nvec times
!!!larger than the data array.  Further, the number of data points
!!!cannot change between iterations.  xi(:,1) and ri(:,1) are the active vectors.
!!!IN:
!!!   this  :: mdiis object
!!!   delta :: coefficient for residual gradient.  Should be between 0 and 1.
!!!   tol   :: target residual
!!!   xi    :: array of vector data.  np X nvec
!!!   ri    :: array of residual data.  np X nvec
!!!   np    :: number of data points (first dimension length)
!!!   nvec   :: length of DIIS vectors
!!!   restart :: restart threshold factor. Ratio of the current residual to the 
!!!              minimum residual in the basis that causes a restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_blas2_new_serial (this,delta,tol,restart)
    implicit none
    type(mdiis_blas2),intent(out) :: this
    _REAL_,intent(in) :: delta, tol, restart

    this%delta = delta
    this%tol = tol
    this%restart = restart
    this%istep=0
  end subroutine mdiis_blas2_new_serial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Create new MDIIS object.  Provide convergence parameters, working
!!!memory and MPI parameters.  Note that the working memory must me
!!!nvec times larger than the data array.  Further, the number of data
!!!points cannot change between iterations.  xi(:,1) and ri(:,1) are
!!!the active vectors.
!!!IN:
!!!   this  :: mdiis object
!!!   delta :: coefficient for residual gradient.  Should be between 0 and 1.
!!!   tol   :: target residual
!!!   xi    :: array of vector data.  np X nvec
!!!   ri    :: array of residual data.  np X nvec
!!!   np    :: number of data points (first dimension length)
!!!   nvec   :: length of DIIS vectors
!!!   restart :: restart threshold factor. Ratio of the current residual to the 
!!!              minimum residual in the basis that causes a restart
!!!   rank  :: MPI process rank
!!!   size  :: Number of processes
!!!   comm  :: MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_blas2_new_mpi(this,delta,tol,restart,&
       rank,size,comm)
    implicit none
    type(mdiis_blas2),intent(out) :: this
    _REAL_,intent(in) :: delta, tol, restart
    integer, intent(in) ::  rank, size, comm
    call mdiis_blas2_new_serial(this,delta,tol,restart)
    this%mpirank = rank
    this%mpisize = size
    this%mpicomm = comm
  end subroutine mdiis_blas2_new_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Destroy MDIIS object. All internal variables and pointers are
!!!reset. Working memory remains intact memory is swapped such that
!!!xi(:,1) and ri(:,1) become the active vectors.
!!!IN:
!!!   this :: mdiis object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_blas2_destroy(this)
    use safemem
    implicit none
    type(mdiis_blas2),intent(inout) :: this
    call mdiis_blas2_reset(this)
    this%delta=0
    this%tol=0
    this%np=0
    this%nvec=0
    if(safemem_dealloc(this%overlap)/=0) call rism_report_error("MDIIS_BLAS2: failed to deallocate OVERLAP")
    if(safemem_dealloc(this%vecMap)/=0) call rism_report_error("MDIIS_BLAS2: failed to deallocate VECRANK")
    nullify(this%xi)
    nullify(this%ri)
  end subroutine mdiis_blas2_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reset MDIIS object. Progress is reset to a single basis
!!!vector. Working memory remains intact memory and is
!!!swapped such that xi(:,1) and ri(:,1) become the active vectors.
!!!IN:
!!!   this :: mdiis object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_blas2_reset(this)
    use safemem
    implicit none
    type(mdiis_blas2),intent(inout) :: this
    this%overlap => safemem_realloc(this%overlap,this%nvec,this%nvec,.false.)
    this%overlap = 0
    !transfer the current working vector to the first index
    if(associated(this%vecMap))then
       if(this%vecMap(1) /=1)then
          call DCOPY(this%np,this%xi(1,this%vecMap(1)),1,this%xi(1,1),1)
          call DCOPY(this%np,this%ri(1,this%vecMap(1)),1,this%ri(1,1),1)
       end if
    end if
    this%istep=0
    this%vecMap => safemem_realloc(this%vecMap,this%nvec,.false.)
    this%vecMap = 0
    if(this%nvec>0)&
         this%vecMap(1) = 1
  end subroutine mdiis_blas2_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Resizes and resets MDIIS object. Progress is reset to a single basis
!!!vector. Working memory remains intact memory and is
!!!swapped such that xi(:,1) and ri(:,1) become the active vectors.
!!!IN:
!!!   this :: mdiis object
!!!   xi    :: array of vector data.  np X nvec
!!!   ri    :: array of residual data.  np X nvec
!!!   np    :: number of data points (first dimension length)
!!!   nvec   :: length of DIIS vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_blas2_resize(this,xi,ri,np,nvec)
    implicit none
    type(mdiis_blas2),intent(inout) :: this
    integer,intent(in) ::  np
    integer, intent(in) :: nvec
    _REAL_,target, intent(in) :: xi(np,nvec), ri(np,nvec)

    !transfer the current working vector to the first index
    this%np = np
    this%nvec = nvec
    this%xi => xi
    this%ri => ri
    call mdiis_blas2_reset(this)
  end subroutine mdiis_blas2_resize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the vector length currently used by the MDIIS object
!!!IN:
!!!   this :: mdiis object
!!!OUT:
!!!   The number of vectors currently used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mdiis_blas2_getCurrentNvec (this) result(nvec)
    implicit none
    type(mdiis_blas2),intent(in) :: this
    integer :: nvec
    nvec = min(this%istep,this%nvec)!count(this%vecMap>0)
  end function mdiis_blas2_getCurrentNvec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the index number of the working vector
!!!IN: 
!!!   this :: mdiis object
!!!OUT:
!!!    the number of the vector to read/write
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mdiis_blas2_getWorkVector (this) result(index)
    implicit none
    type(mdiis_blas2),intent(in) :: this
    integer :: index
    index = this%vecMap(1)
  end function mdiis_blas2_getWorkVector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!One MDIIS iteration. Query mdiis_getVectorNumer() to get the new active vector.
!!!IN:
!!!   this  :: mdiis object
!!!   rms1  :: residual at the end of the step
!!!   conver:: have we converged
!!!   tolerance :: tolerance for this calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  mdiis_blas2_advance (this, rms1,conver,tolerance_o)
#ifdef MPI
    use mpi
#endif
    implicit none 
#include "def_time.h" 
    
    type(mdiis_blas2), intent(inout) :: this
    logical, intent(out) ::  conver
    _REAL_, intent(out) ::  rms1
    _REAL_, optional, intent(in) ::  tolerance_o

    _REAL_ :: tolerance

    ! indx :: for LAPACK routines
    integer ::  indx(0:this%nvec)
    _REAL_ :: toverlap(1:this%nvec)

    ! aij :: overlap matrix for LAPACK
    ! bi  :: linear coefficients from LAPACK
    _REAL_, allocatable :: aij(:,:), bi(:,:)

    ! nvecWRK :: number of vectors with data
    integer :: nvecWRK

    ! counter variables
    integer :: is2, ivec
    integer :: err = 0
    tolerance = this%tol
    if (present(tolerance_o)) tolerance = tolerance_o
#ifdef RISM_DEBUG
    write(0, * )"MDIIS start", this%delta, rms1, tolerance, this%np, this%nvec, conver
#endif

    !.................. initialize counters and switches ...................
    nvecWRK = count(this%vecMap > 0)
    !.................. increment step and MDIIS counters ..................
    this%istep = this%istep + 1

    allocate(aij(0:this%nvec, 0:this%nvec), bi(0:this%nvec,1))

    !............. calculate diagonal overlap of new residual ..............
    !                              and
    !........... calculate nondiagonal overlaps of new residual ............

    !This is dot product of the new residual vector with itself and
    !all of the older residuals.  Results from previous iterations are
    !also present in the overlap array.  New results are first stored
    !in the vecMap(1) row of the overlap matrix. These are then copied
    !to the vecMap(1) column.  Note that the overlap array is a
    !symmetric matrix.
    call timer_start(TIME_MDIIS_DATA)
#if defined(MPI)      
    call DGEMV ('T', this%np, nvecWRK, 1d0, this%ri, this%np, &
         this%ri(1, this%vecMap(1)), 1, 0d0, toverlap, 1)
    CALL MPI_AllREDUCE(toverlap(1:nvecWRK), this%overlap(this%vecMap(1), 1:nvecWRK), nvecWRK, &
         MPI_DOUBLE_PRECISION, MPI_SUM, this%mpicomm, err)
    if (err /= 0) call rism_report_error &
         ("MDIIS_BLAS2_UPDATE: could not reduce OVERLAP")
#else
    call DGEMV ('T', this%np, nvecWRK, 1d0, this%ri(1,1), this%np, &
         this%ri(1, this%vecMap(1)), 1,0d0, this%overlap(this%vecMap(1), 1), this%nvec)
#endif /* defined(MPI) */


    this%overlap(:nvecWRK, this%vecmap(1)) = this%overlap(this%vecmap(1), :nvecWRK)
    call timer_stop(TIME_MDIIS_DATA)

    !................ get mean square value of new residual ................
    rms1 = sqrt( this%overlap(this%vecMap(1), this%vecMap(1)) / (this%np * this%mpisize))
    !.............. check mean square residual for tolerance ...............
    if (rms1 <= tolerance)  then
       conver = .true.
       return
    else
       conver = .false.
    end if

    call checkrestart(this, rms1, nvecWRK)

    !We now solve for the non - trivial coefficients that minimizes the
    !linear combination of residuals.  I.e., we get the eigenvalues of
    !the aij matrix (that contains the overlap matrix as a submatrix).

    ! Note that the form is different than in the paper.  Here we have
    !
    !       |  0  -1  ... -1  |
    ! aij = | -1  S11 ... S1M |
    !       |  :   :       :  |
    !       | -1  SM1 ... SMM |
    !
    !      | lambda |
    ! X  = |   C1   |
    !      |   :    |
    !      |   CM   |
    !
    !      | -1 |
    ! bi = |  0 |
    !      |  : |
    !      |  0 |

    call timer_start(TIME_MDIIS_LAPACK)
    !......................... load DIIS matrices ..........................
    aij(0, 0) = 0d0
    bi = 0
    bi(0, 1) = -1d0
    aij(1:nvecWRK, 0) = -1d0 
    aij(0, 1:nvecWRK) = -1d0
    do is2 = 1, nvecWRK
       call DCOPY(nvecWRK, this%overlap(1, is2), 1,aij(1, is2), 1)
    end do

    call DGESV(nvecWRK + 1, 1,aij, this%nvec + 1, indx, bi, this%nvec + 1, err)
    if (err < 0) then
       err = err * ( - 1)
       call rism_report_error("Linear equation solver failed.")
    end if

    !......... get DIIS minimum, MDIIS correction, and next point ..........
    !update vecMap so we know which will be the new active vector.  If
    !this vector has data, it should the be the oldest data and will
    !be overwritten
    !the equivalent <where> statement gives an error in Valgrind for
    !the first instance but this works fine. Go figure.
    do ivec = 1, nvecWRK
       if (this%vecMap(ivec) == this%nvec)&
            this%vecMap(ivec) = 0
    end do
    this%vecMap(1:min(nvecWRK + 1, this%nvec)) = this%vecMap(1:min(nvecWRK + 1, this%nvec)) +1
    nvecWRK = count(this%vecMap > 0)

    !DIIS part.  New vector is a linear combition of all the previous
    !vectors.  This must be done in two parts (before and after
    !vecMap(1) vector) to avoid memory aliasing.
    if (this%vecMap(1) == 1 ) then
       !if vecMap(1) == 1, DGEMV thinks there is no work to be done
       !and the scalar product is also skipped.
       call DSCAL(this%np, bi(this%vecMap(1), 1), this%xi(1, this%vecMap(1)), 1)
    else
       call DGEMV ('N', this%np, this%vecMap(1) - 1, 1d0, this%xi(1,1), this%np, &
            bi(1, 1), 1,bi(this%vecMap(1), 1), this%xi(1, this%vecMap(1)), 1)
    end if
    if (nvecWRK - this%vecMap(1) > 0)&
         call DGEMV ('N', this%np, nvecWRK - this%vecMap(1), &
         1d0, this%xi(1, this%vecMap(1) + 1), this%np, &
         bi(this%vecMap(1) + 1, 1), 1,1d0, this%xi(1, this%vecMap(1)), 1)

    ! calculate the projected residual for the modified update.  It is
    ! expected that the external program will overwrite these values.
    if (this%vecMap(1) == 1 ) then
       !if vecMap(1) == 1, DGEMV thinks there is no work to be done
       !and the scalar product is also skipped.
       call DSCAL(this%np, bi(this%vecMap(1), 1), this%ri(1, this%vecMap(1)), 1)
    else
       call DGEMV ('N', this%np, this%vecMap(1) - 1, 1d0, this%ri(1,1), this%np, &
            bi(1, 1), 1,bi(this%vecMap(1), 1), this%ri(1, this%vecMap(1)), 1)
    end if
    if (nvecWRK - this%vecMap(1) > 0)&
         call DGEMV ('N', this%np, nvecWRK - this%vecMap(1), &
         1d0, this%ri(1, this%vecMap(1) + 1), this%np, &
         bi(this%vecMap(1) + 1, 1), 1,1d0, this%ri(1, this%vecMap(1)), 1)
    
    !Modified part.  Add the predicted residual to multiplied by the
    !step size to the DIIS result
    ! call DGEMV ('N', this%np, nvecWRK, this%delta, this%ri, this%np,&
    !      bi(1,1),1,1d0,this%xi(1,this%vecMap(1)),1)
    call DAXPY (this%np, this%delta, this%ri(1, this%vecMap(1)), 1, &
         this%xi(1,this%vecMap(1)), 1)
    call timer_stop(TIME_MDIIS_LAPACK)
    deallocate(aij, bi)

  end subroutine mdiis_blas2_advance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Checks if a restart is require (rms1 > this%restart * smallest rms in memory).
!!!If so, the lowest rms vector is chosen and all others are discarded.
!!!IN:
!!!   this  :: mdiis object
!!!   rms1  :: residual at the end of the step
!!!   nvecWRK  :: number of vectors with data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine checkRestart (this,rms1,nvecWRK)
    implicit none
    type(mdiis_blas2), intent(inout) :: this
    _REAL_, intent(in) :: rms1
    integer, intent(inout) :: nvecWRK

    _REAL_ ::  rmsmin, rmsmi2
    integer :: isirst, isimin, isimi2, is

    !............... get two residuals, minimal in DIIS set ................
    rmsmin = this%overlap(1, 1)
    isimin = 1
    rmsmi2 = rmsmin
    isimi2 = isimin
    do is = 2, nvecWRK
       if (this%overlap(is, is) < rmsmin)  then
          rmsmi2 = rmsmin
          isimi2 = isimin
          rmsmin = this%overlap(is, is)
          isimin = is
       else if (this%overlap(is, is) < rmsmi2)  then
          rmsmi2 = this%overlap(is, is)
          isimi2 = is
       end if
    end do
    rmsmin = sqrt( rmsmin / this%np / this%mpisize)
    rmsmi2 = sqrt( rmsmi2 / this%np / this%mpisize)

    !................ if convergence is poor, restart MDIIS ................
    if (nvecWRK > 1 .AND. rms1 > this%restart * rmsmin)  then
       !.......... choose restarting vector so as to prevent cycling .........

       !if nvecWRK = 1 then we have problems
       !if nvecWRK = 2 ditch the old vector and just use the new vector
       !otherwise we use the minimal residue vector
       
       if (nvecWRK == 1) then
       else if (nvecWRK == 2) then
          isirst = 1
       else
          isirst = isimin
       end if

       !................... restore vector to restart from ...................
       if (isirst /= 1)  then
          call timer_start(TIME_MDIIS_DATA)
          call DCOPY(this%np, this%ri(1, isirst), 1, this%ri(1,1), 1)
          call DCOPY(this%np, this%xi(1, isirst), 1, this%xi(1,1), 1)
          call timer_stop(TIME_MDIIS_DATA)
          this%overlap(1, 1) = this%overlap(isirst, isirst)
       end if

       !.................... reset MDIIS vectors counters .....................
       nvecWRK = 1
       this%vecMap = 0
       this%vecMap(1) = 1
       this%istep = 1
    end if
  end subroutine checkRestart
end module mdiis_blas2_c
