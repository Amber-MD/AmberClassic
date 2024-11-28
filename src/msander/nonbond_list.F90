!<compile=optimized>
#include "../include/assert.fh"
#include "../include/dprec.fh"

module nblist

_REAL_, dimension(:,:), allocatable, save :: imagcrds, fraccrd, savfrac, &
                                             dfrac, savcrd
integer, dimension(:),  allocatable, save :: indatg, atmcell, numatg, indoff, &
                                             bckptr, nlogrid, nhigrid, &
                                             my_grids, lstmask, nummask, &
                                             maskptr, iwa, iwh, iws, &
                                             numvdw, numhbnd, numsc, atmlist, &
                                             mygrdlist, numimg
integer, dimension(:), allocatable, save :: indexlo, indexhi
integer, dimension(:), allocatable, save :: nvdwcls
integer, dimension(:,:),  allocatable, save :: nghbptr, nghtran
integer, dimension(:),    allocatable, save :: itran
integer, save :: nucgrd1_0, nucgrd2_0, nucgrd3_0
integer, dimension(1:7,1:10), save :: xtran
integer, save :: myindexlo, myindexhi, nblist_allint, nblist_allreal
integer, private, allocatable, dimension(:), save :: exclude

logical :: first_list_flag
integer, save :: num_calls_nblist=0

_REAL_, dimension(1:3,1:18), save :: tranvec

! Contents of ew_unitcell.h
integer, parameter :: BC_EWUCR=54, BC_EWUCI=3
!#define BC_EWUCR 54
!#define BC_EWUCI 3

_REAL_ :: a, b, c, alpha, beta, gamma, volume
_REAL_ :: ucell, recip, dirlng, reclng, sphere, cutoffnb
_REAL_ :: olducell, oldrecip
_REAL_ :: skinnb, cutlist, nbfilter
common/unitcell/ &
      ucell(3,3), recip(3,3), dirlng(3), reclng(3),   &! 24
      olducell(3,3), oldrecip(3,3),                   &! 42
      a, b, c, alpha, beta, gamma, volume,            &! 49
      sphere, cutoffnb, skinnb, cutlist, nbfilter      ! 54

integer :: nbflag, nbtell, steps_since_list_build
common/ew_upd/nbflag, nbtell, steps_since_list_build
! End contents of ew_unitcell.h

! Contents of ew_localnb.h
!
! Unit cell grid (NUCGRD) of mapped unit cell coords for preimaging.  Here, a
! unit cell defines one of the (equally sized) partitions of the larger
! simulation cell. nucgrd1, nucgrd2, and nucgrd3 are the numbers of unit cells
! along the first, second, and third box vectors.
!
! The neighborhood of each subcell is obtained by considering cells within +-
! nghb1, nghb2, and nghb3 of one another along the first, second, and third
! box vectors, respectively.

! The distance between parallel faces of a subcell is then reclng(1)/nucgrd1,
! reclng(2)/nucgrd2, or reclng(3)/nucgrd3, where reclng was determined in the
! ew_box.F90 module as the distance between parallel faces of the simulation
! cell (this can account for non-orthogonal box angles).  Therefore, the short
! range cutoff is the minimum of nghb1*reclng(1)/nucgrd1,
! nghb2*reclng(2)/nucgrd2 and nghb3*reclng(3)/nucgrd3.  mxatcell is the maximum
! number of atoms per subcell.  Note that nghb1 < nucgrd1, nghb2 < nucgrd2, and
! nghb3 < nucgrd3 are required.  imagptr(i) is the subcell number in the imaged
! grid corresponding to subcell i in the unit cell grid.  nghbptr is an array
! from which the neighbor atoms in the image grid can be rapidly retrieved.

! neighbor cell #
integer nghb
parameter (nghb = 3)
integer, parameter :: BC_DIRPARS=13
integer :: numnptrs, nucgrd1, nucgrd2, nucgrd3, nucgmax, nucgrd
integer :: nghb1, nghb2, nghb3
integer :: maxnblst, maxnptrs, maximage, mxlstmsk

common/dirpars/ &
      numnptrs, nucgrd1, nucgrd2, &
      nucgrd3, nghb1, nghb2, nghb3, nucgmax, nucgrd, &
      maxnblst, maxnptrs, maximage, mxlstmsk

integer, parameter :: BC_NB_GINDEX = 2
integer gindexlo, gindexhi

common/nb_gindex/gindexlo, gindexhi
! End contents of ew_localnb.h

contains

!------------------------------------------------------------------------------
! nblist_allocate: allocate memory to store a non-bonded pairlist (a Verlet
!                  list).  This routine will report an estimate of the Verlet
!                  list's memory usage.
!
! Arguments:
!   natom:       the number of atoms in the system
!   ntypes:      the number of atom types in the topology
!   num_direct:  number of CPU cores allotted to direct space calculations
!   numtasks:    number of tasks in the direct space calculation that can be
!                split amongst the CPU cores
!------------------------------------------------------------------------------
subroutine nblist_allocate(natom, ntypes, num_direct, numtasks)

  implicit none

  integer, intent(in) :: natom,ntypes, num_direct, numtasks
  integer :: allocate_error, allint, allreal

  allreal = 0   ! Report back how many reals have been allocated
  allint = 0    ! Report back how many integers have been allocated

  allocate(imagcrds(1:3,1:natom), stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + 3*natom

  allocate(fraccrd(1:3,1:natom), savfrac(1:3,1:natom), dfrac(1:3,1:natom), &
           savcrd(1:3,1:natom), stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allreal = allreal + 12*natom

  allocate (indatg(1:natom), atmcell(1:natom), atmlist(1:natom+natom/10), &
            stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + 3*natom + natom/10

  allocate(nlogrid(1:nucgmax), nhigrid(1:nucgmax), my_grids(1:nucgmax), &
           indoff(1:nucgmax), numimg(1:nucgmax), numatg(1:nucgmax), &
           stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + 6*nucgmax

  allocate(nghbptr(0:maxnptrs,1:nucgmax), bckptr(1:natom), &
           stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint+natom + (maxnptrs+1)*nucgmax

  allocate(nummask(1:natom), maskptr(1:natom), lstmask(1:mxlstmsk), &
           stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + 2*natom + mxlstmsk

  allocate(iwa(1:natom), iwh(1:natom), numvdw(1:natom), numhbnd(1:natom), &
           stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + 4*natom

#ifdef MPI /* SOFT CORE */
  allocate(iws(1:natom), numsc(1:natom), stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + 2*natom
#endif

  allocate(itran(1:natom+natom/10), stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + natom + natom/10

  allocate(indexlo(0:numtasks-1), indexhi(0:numtasks-1), stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + 2*numtasks

  ! Private arrays
  allocate(exclude(1:natom), mygrdlist(1:natom), stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + 2*natom
   
#ifdef MPI
  allocate(nghtran(1:maxnptrs, 1+nucgmax/num_direct), stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + maxnptrs*(nucgmax/num_direct + 1)
#else
  REQUIRE(num_direct == 1 .or. .true. ) ! hack to prevent warnings
  allocate (nghtran(1:maxnptrs, 1:nucgmax+1), stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + maxnptrs*(nucgmax + 1)
#endif

  allocate (nvdwcls(1:ntypes), stat=allocate_error)
  REQUIRE(allocate_error == 0)
  allint = allint + ntypes
  nblist_allint  = allint
  nblist_allreal = allreal

  return
end subroutine nblist_allocate

!------------------------------------------------------------------------------
! nblist_deallocate: de-allocate memory for the non-bonded Verlet list.
!------------------------------------------------------------------------------
subroutine nblist_deallocate()
  implicit none

  integer :: allocate_error

  if (.not. allocated(imagcrds) ) return ! None of this memory has been
                                         !  allocated, just return

  deallocate(imagcrds, stat=allocate_error)
  REQUIRE(allocate_error == 0)

  deallocate(fraccrd, savfrac, dfrac, savcrd, stat=allocate_error)
  REQUIRE(allocate_error == 0)

  deallocate (indatg,atmcell,atmlist, stat=allocate_error)
  REQUIRE(allocate_error == 0)

  deallocate(nlogrid, nhigrid, my_grids, indoff, numimg, numatg, &
             stat=allocate_error)
  REQUIRE(allocate_error == 0)

  deallocate(nghbptr, bckptr, stat=allocate_error)
  REQUIRE(allocate_error == 0)

  deallocate(nummask, maskptr, lstmask, stat=allocate_error)
  REQUIRE(allocate_error == 0)

  deallocate(iwa, iwh, numvdw, numhbnd, stat=allocate_error)
  REQUIRE(allocate_error == 0)

#ifdef MPI /* SOFT CORE */
  deallocate (iws, numsc, stat=allocate_error)
  REQUIRE(allocate_error == 0)
#endif

  deallocate (itran, stat=allocate_error)
  REQUIRE(allocate_error == 0)

  deallocate(indexlo, indexhi, stat=allocate_error)
  REQUIRE(allocate_error == 0)

  ! Private arrays
  deallocate(exclude, mygrdlist, stat=allocate_error)
  REQUIRE(allocate_error == 0)
   
  deallocate(nghtran, stat=allocate_error)
  REQUIRE(allocate_error == 0)

  deallocate(nvdwcls, stat=allocate_error)
  REQUIRE(allocate_error == 0)

  return

end subroutine nblist_deallocate

!------------------------------------------------------------------------------
! nonbond_list: handles set-up and error checking for calling of get_nb_list,
!               which creates the nonbond list.
!
! Arguments:
!   crd:            coordinates of all atoms in the system
!   iac:            atom types of each atom in the system, which serve as
!                   indices into the two-dimensional array ico (see below)
!   ico:            used to look up hydrogen bonding interactions, an array of
!                   ntypes x ntypes integers that gives indices into tables of
!                   van-der Waals and hydrogen bonding coefficients.  In
!                   practice this gets unrolled into a one-dimensional array.
!   iblo:           the number of atoms masked by i, as listed in the array
!                   inb.
!   inb:            indices of other atoms masked by each atom in the system.
!                   These indices are stored back-to-back.
!   ntypes:         the number of atom types to consider
!   natom:          number of atoms in the system
!   x:              total dynamic _REAL_ memory allocated by sander
!   ix:             total dynamic integer memory allocated by sander
!   ipairs:         indices of non-bonded atom pairs
!   ntnb:           
!   ibelly:         indices of atoms in a belly
!   belly:          flag to indicate that there is a belly in effect
!   newbalance:
!   qsetup:         flag to have QM/MM simulation setup take place, also
!                   indicates that this is the first call to the subroutine
!   do_list_update: flag to have the non-bonded list reconstructed
!------------------------------------------------------------------------------
subroutine nonbond_list(crd, iac, ico, iblo, inb, ntypes, natom, x, ix, &
                        ipairs, ntnb, ibelly, belly, newbalance, qsetup, &
                        do_list_update)
#ifdef MPI
   use mpi
#endif
   implicit none

   integer :: natom, ntnb, newbalance
   integer :: ipairs(*), ix(*), ibelly(*)
   _REAL_  :: crd(3,natom)
   integer :: iac(natom), ico(*), iblo(*), inb(*), ntypes
   _REAL_  :: x(*)
   logical qsetup

   integer :: ngrd1, ngrd2, ngrd3, sizgrdprs
   integer last_numlist, isiz
   logical belly, trial, balance
#  include "extra.h"
#  include "ew_cntrl.h"
#  include "ew_mpole.h"
   
#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
   integer listdiff, listdiffmax
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
   integer tmplist(0:MPI_MAX_PROCESSORS), alllist(0:MPI_MAX_PROCESSORS)
#else   /* not parallel needs numtasks and mytaskid */
   integer numtasks, mytaskid
#endif
   
!#  include "ew_tran.h"
#  include "def_time.h"
   
   integer ier
   integer, dimension(:), allocatable :: nnghbptr 
   integer, dimension(:), allocatable :: nghtranptr
   integer, dimension(:), allocatable :: gridpairs   
   integer ifail,listtot,listtotall
#ifdef MPI
   integer inddel, i, ierr
#endif
   logical, intent(out) :: do_list_update
   save trial, balance
   save last_numlist
#ifdef MPI
   save listdiffmax
#endif
   
   ! Code starts here
#ifndef MPI
   mytaskid=0
   numtasks=1
#endif


  ! First time setup stuff: if qsetup is true then this is the first
  ! time through so do the setup stuff.  The qsetup flag will not be
  ! set to false until after the if(.not. qsetup ) block.
  if ( qsetup ) then
    call ew_startup(natom, iblo, inb, x, ix)
    nucgrd1_0 = -1
    nucgrd2_0 = -1
    nucgrd3_0 = -1
    call save_crds(natom, crd)
    balance = .false.
    newbalance = 0
    if (periodic == 0) then
      balance = (numtasks > 1)
    end if
    if (balance) then
      newbalance = 2
    end if
    last_numlist = 0
    steps_since_list_build = 0
  end if

  ! Done with first-time setup stuff ?
  trial = (newbalance > 0)
  if (.not. qsetup) then
     
    ! The skin check is done in runmd for spatial decomposition.
    ! do_list_update is passed in as optional arg.  A skin (buffer)
    ! check will always be done to see if new list is required,
    ! unless this is the first pass or if the nocutoff flag is set.
    ! In these cases such a list is undoubtedly required or its
    ! contents are easily seen to be all pairs in the system.
    if (nocutoff) then
      return
    end if

    !   if nbflag=0, use old logic; i.e. only update if ntnb .ne. 0
    !   if nbflag=1, the skin check determines it
    if (nbflag == 0) then
      if (nbtell .ne. 0) then ! list update info requested
        if(master)write(6,'(1x,A,I2)') 'OLD LIST LOGIC; ntnb = ',ntnb
      end if
      if (ntnb == 0) then
        return
      end if
    else
      call check_skin(crd,do_list_update)
      if (do_list_update) then
        call save_crds(natom,crd)
      else
        return
      end if
    end if
  end if
  qsetup = .false.
  ! End contingency for (.not. qsetup)

  ! Start the list build.  First, do the list grid setup.
  num_calls_nblist = num_calls_nblist + 1
  call timer_start(TIME_BLDLST)
  if (master .and. nbtell .ne. 0) then
 
    ! List update info requested
    write(6, '(1x,A,I7,A)') 'Building list: ', &
          steps_since_list_build, ' steps since previous list build.'
  end if
  steps_since_list_build = 0

#ifdef MPI
  if (i_do_direct) then
    call mpi_comm_size(direct_comm, numtasks, ierr)
    call mpi_comm_rank(direct_comm, mytaskid, ierr)
#endif
  ! Get nb grid information.  This is necessary in the case where the
  ! system dimensions are changing: nonperiodic, or constant pressure.
  call map_coords(crd,natom,recip)
  call save_frac_crds(natom)
  call setup_grids(periodic, nogrdptrs, verbose)
  if (nucgrd1*nucgrd2*nucgrd3 > nucgmax) then
    call sander_bomb('nonbond_list', &
                     ' volume of ucell too big, too many subcells', &
                     ' list grid memory needs to be reallocated, &
                     &restart sander')
  end if

  ! Contingency for nonperiodic systems:
  if (periodic == 0) then

    ! If system has too few cells for the pointer method to be efficient, use
    ! the no-grid-pointer system where all forward cells are checked for pair
    ! atoms rather than just the forward cells on the grid-pointer-list.
    nogrdptrs = nogrdptrs .or. (nucgrd1 <= 2*nghb + 2 .or. &
                                nucgrd2 <= 2*nghb + 2 .or. &
                                nucgrd3 <= 2*nghb + 2)
    if (verbose >= 3) then
      write(6,*) " List using nogrdptrs ", nogrdptrs
    end if

    ! Make the subcells about 3 A in size now so that there
    ! are more subcells, easier to load balance.
    if (nogrdptrs) then
      if (dirlng(1)/nucgrd1 > 3) then
        ngrd1 = int(dirlng(1)/3)
        ngrd2 = int(dirlng(2)/3)
        ngrd3 = int(dirlng(3)/3)
        if (verbose >= 1) then
          write(6,'("|    New Grids set up for nogrdptrs ")')
          write(6, '(5X,a,/,5X,i9,1X,i9,1X,i9)') &
                'Number of grids per unit cell in x,y,z:', ngrd1, ngrd2, ngrd3
        end if
        nucgrd1 = ngrd1
        nucgrd2 = ngrd2
        nucgrd3 = ngrd3
      end if
    end if
  end if
  ! End of contingency for nonperiodic grid

  call fill_tranvec()
  nucgrd = nucgrd1*nucgrd2*nucgrd3

#ifdef MPI
    ! For setup in a parallel run, an initial guess at the distribution
    ! of subcells must be made, thus "trial" is TRUE.  For a periodic system,
    ! the balance is assumed to be good for evenly dividing the subcells,
    ! meaning the variable "balance" is FALSE.
      if (trial .or. .not. balance) then
        inddel = (nucgrd-1) / numtasks + 1
        if (inddel == 0) then
          inddel = 1
        end if
        myindexlo = 1 + (mytaskid)*inddel
        myindexhi = myindexlo + inddel - 1
        if (mytaskid == numtasks-1) then
          myindexhi = nucgrd
        end if
        last_numlist = 0
      end if
#else
  ! For non parallel runs, do all the subcells
  myindexlo = 1
  myindexhi = nucgrd
#endif

  ! Assign atoms to cells (atmcell) and make cell atom list (indatg)
  call grid_ucell(natom, periodic)
  if (periodic == 0) then

    ! Nonperiodic systems:
    isiz = max(2, (myindexhi-myindexlo+1)*(maxnptrs+1))
    if (nogrdptrs) then
      isiz = 1
    end if
    allocate(nnghbptr(isiz), stat=ier)
    REQUIRE(ier == 0)
    isiz = max(2, (myindexhi-myindexlo+1)*(maxnptrs))
    if (nogrdptrs) then
      isiz = 1
    end if
    allocate(nghtranptr(isiz), stat=ier)
    REQUIRE(ier == 0)
    call grid_pointers(nnghbptr, nghtranptr, periodic, nogrdptrs)
  else

    ! Periodic systems: only call grid_pointers if the unit cell has changed
    ! so much that nucgrd[123] have changed.  Otherwise the grid pointers do
    ! not change.
    if (nucgrd1 .ne. nucgrd1_0 .or. nucgrd2 .ne. nucgrd2_0 .or. &
        nucgrd3 .ne. nucgrd3_0) then
      nucgrd1_0 = nucgrd1
      nucgrd2_0 = nucgrd2
      nucgrd3_0 = nucgrd3
      call grid_pointers(nghbptr, nghtran, periodic, nogrdptrs)
    end if
  end if
  ! End branch for periodicity considerations
      
  call grid_image(natom, verbose, periodic)
  if (balance) then
    sizgrdprs = nucgrd*2
  else
    sizgrdprs=2
  end if

  allocate(gridpairs(sizgrdprs), stat=ier)
  REQUIRE(ier == 0)
  if (periodic == 0) then
    call get_nb_list(iac, ico, ntypes, ifail, listtot, natom, ipairs, &
                     nnghbptr, nghtranptr, tranvec, belly, ibelly, balance, &
                     gridpairs, periodic, nogrdptrs)
  else
    call get_nb_list(iac, ico, ntypes, ifail, listtot, natom, ipairs, &
                     nghbptr, nghtran, tranvec, belly, ibelly, balance, &
                     gridpairs, periodic, nogrdptrs)
  end if
  if (ifail == 1) then
    write(6, '(5x,a,i10)') 'SIZE OF NONBOND LIST = ', listtot
    call sander_bomb('nonbond_list','Non bond list overflow!', &
                     'check MAXPR in locmem.F90')
  end if

#ifdef MPI
    ! Attempt load balancing in parallel simulations
    if (trial) then
      ASSERT(periodic == 0)
      call mpi_allreduce(gridpairs, gridpairs(1+nucgrd), nucgrd, &
                         mpi_integer, mpi_sum, commsander, ierr)
      call fix_grid_balance(gridpairs(1+nucgrd), nucgrd, numtasks, mytaskid, &
                            listdiffmax)
      deallocate(gridpairs, nghtranptr, nnghbptr)
      isiz = max(2, (myindexhi-myindexlo+1)*(maxnptrs+1))
      if (nogrdptrs) then
        isiz = 1
      end if
      allocate(nnghbptr(isiz), stat=ier)
      REQUIRE(ier == 0)
      isiz = max(2, (myindexhi-myindexlo+1)*(maxnptrs))
      if (nogrdptrs) then
        isiz = 1
      end if
      allocate(nghtranptr(isiz), stat=ier)
      REQUIRE(ier == 0)
      allocate(gridpairs(2*nucgrd), stat=ier)
      REQUIRE(ier == 0)
      call grid_pointers(nnghbptr, nghtranptr, periodic, nogrdptrs)
      call grid_image(natom, verbose, periodic)
      call get_nb_list(iac, ico, ntypes, ifail, listtot, natom, ipairs, &
                       nnghbptr, nghtranptr, tranvec, belly, ibelly, balance, &
                       gridpairs, periodic, nogrdptrs)
      trial = .false.
      last_numlist = listtot
      newbalance = 0
    else if (balance) then

      ! If this was not a rebalance step, then check whether the list has
      ! gotten more than a tolerance away from the last balance. If so,
      ! trigger a new balance by setting the variables trial and newbalance.
      ! newbalance will need to be communicated to all processors before the
      ! next list build.
      listdiff = iabs(last_numlist-listtot)
      if (listdiff > listdiffmax) then
        newbalance = 2
      end if
    end if 
    ! End contingency for attempting load balancing

    if ((nbtell > 1) .or. (balance .and. (verbose > 0))) then
      if (master) then
        write(6,105)
      end if
      if (master) then
        write(6,106)
      end if
      do i = 0, numtasks-1
        tmplist(i) = 0
      end do
      tmplist(mytaskid) = listtot
      call mpi_allreduce(tmplist(0), alllist(0), numtasks, mpi_integer, &
                         mpi_sum, commsander, ierr)
      if (master) then
        write(6,110) (i, alllist(i), i=0, numtasks-1)
        do i = 1, numtasks-1
          tmplist(0) = tmplist(0) + alllist(i)
        end do
        write(6,103) tmplist(0)
      end if
      110 format("|       ", 2i12)
      103 format("|            Total: ", i12)
      105 format("| ", '----------------List Breakdown----------------')
      106 format("|       ", 'list processor   listtot')
    end if
#else
  if (master) then
    if (verbose > 0) then
      write(6,*) 'listtot = ', listtot
    end if
  end if
#endif /* MPI */
      
  deallocate(gridpairs)
  if (periodic == 0) then
    deallocate(nghtranptr, nnghbptr)
  end if
  listtotall = listtot
#ifdef MPI
    if (first_list_flag .or. (nbtell >= 1)) then
      call mpi_allreduce(listtot, listtotall, 1, mpi_integer, &
                         mpi_sum, commsander, ierr)
    end if
    call mpi_comm_rank(commsander,mytaskid,ierr)
    call mpi_comm_size(commsander,numtasks,ierr)
  end if
  call mpi_barrier(commsander,ierr)
  ! End contingency for i_do_direct, started circa line 415.
  ! The inclusion of an extra conditional that only exists in
  ! the parallel code gave rise to some funky indentation.
#endif /* MPI */

#ifndef API
  if ((first_list_flag .or. (nbtell == 1)) .and. master) then
    write(6, '(a,i10)') '| Local SIZE OF NONBOND LIST = ', listtot
    write(6, '(a,i10)') '| TOTAL SIZE OF NONBOND LIST = ', listtotall
  end if
  call flush(6)
#endif /* API */
  first_list_flag = .false.
  call timer_stop(TIME_BLDLST)

  return

end subroutine nonbond_list 

!------------------------------------------------------------------------------
! get_nb_list: this routine lays out the short range part of the non-bonded
!              calculation.  This is configured for parallelism, using a
!              gridding routine, or geometric hashing with renumbering of atoms
!              to speed list generation and promote locality of reference,
!              important for cache memory usage.  The pre-imaging is to avoid
!              expensive coordinate transforms in the minimum image pair
!              calculations.  Pointers have been set up between subcells to
!              speed this process.
!
! Arguments:
!   iac:
!   ico:
!   ntypes:
!   ifail:
!   listtot:
!   natom:
!   ipairs:
!   nghbptr:
!   nghtran:
!   tranvec:
!   belly:
!   ibelly:
!   balance:
!   gridpairs:
!   periodic:
!   nogrdptrs: flag to indicate that there are no grid pointers in effect
!------------------------------------------------------------------------------
subroutine get_nb_list(iac, ico, ntypes, ifail, listtot, natom, ipairs, &
                       nghbptr, nghtran, tranvec, belly, ibelly, balance, &
                       gridpairs, periodic, nogrdptrs)

#ifdef MPI
   use mpi
#endif
   implicit none
   
#  include "extra.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#else   /* not parallel needs numtasks and mytaskid */
   integer numtasks,mytaskid
#endif
   
  ! Formal arguments
  integer iac(*), ico(*), ntypes, ifail
  integer listtot, natom
  integer ipairs(maxnblst)
  integer nghbptr(0:maxnptrs,*)
  integer nghtran(maxnptrs,*)
  _REAL_  tranvec(*)
  logical belly
  integer ibelly(*)
  logical, intent(in) :: balance  ! unused in serial code
  integer gridpairs(nucgrd,2)     ! unused in serial code
  integer periodic
  logical nogrdptrs

  ! Local variables
  integer numlist
  integer index, index0, index1
  integer j, jj, j0, k, ncell_lo, ncell_hi, m, m1, m2
  integer i, kk, numpack
  integer nstart, tranyz, xtindex, jtran, indi, nstop, numindex
  _REAL_ xk, yk, zk
  _REAL_ cutoffsq

#ifdef MPI
  integer np0
#endif
   
#ifndef MPI
  numtasks = 1
  mytaskid = 0
#endif
   
  mygrdlist(1:natom) = 0
  ifail = 0
  cutoffsq = cutlist*cutlist 
  exclude(1:natom) = 0  

  numpack = 0
# ifdef MPI
  if (balance) then
    np0 = 0
    do i = 1, nucgrd
      gridpairs(i,1) = 0
    end do
  end if
# endif
  do index = myindexlo, myindexhi
    index0 = index - myindexlo + 1
    if (numimg(index) > 0) then
      ncell_lo = nlogrid(index)
      ncell_hi = nhigrid(index)
      numlist = 0

      ! Contingency for the pointer-free approach
      if (nogrdptrs) then
        if (periodic /= 0) then
          call sander_bomb('get_nb_list', 'Cannot run nogrdptrs with ntb=1', &
                           'turn off nogrdptrs ')
        end if
        do j0 = index, nucgrd
          m1 = nlogrid(j0)
          m2 = nhigrid(j0)
          if (m2 >= m1) then
            do m = m1, m2
              numlist = numlist+1
              atmlist(numlist) = m
            end do
          end if
        end do
        if (numlist > 0) then
          do k = ncell_lo, ncell_hi
            kk = k - ncell_lo + 1
            i = bckptr(k)
            xk = imagcrds(1,k)
            yk = imagcrds(2,k)
            zk = imagcrds(3,k)
            mygrdlist(k) = 1
            call pack_nb_nogrdptrs(kk, i, xk, yk, zk, imagcrds, cutoffsq, &
                                   numlist, numpack, iac, ico, ntypes, &
                                   ipairs, ifail, belly, ibelly)
            if (ifail == 1) then
              listtot = maxnblst
              return
            end if
          end do
        end if
      else

        ! Engage the approach based on grid pointers.  Get the list of
        ! atoms in the neighborhood of cell(index0).
        numindex = numnptrs
        if (periodic == 0) then
          numindex = nghbptr(0,index0)
        end if
        do j0 = 1, numindex
          index1 = nghbptr(j0,index0)
          if (j0 == 1) then
            nstart = 0
          else
            nstart = -nghb1
          end if
          nstop = nghb1
          if (periodic == 0) then
            indi = mod(index1-1, nucgrd1) + 1
            if (indi > nucgrd1-nghb1) then
              nstop = nucgrd1 - indi
            end if
            if (indi <= nghb1) then
              nstart = 1 - indi
            end if
            if (j0 == 1) then
              nstart = 0
            end if
          end if
          xtindex = ishft(nghtran(j0,index0), -8)
          REQUIRE(xtindex <= 10 .and. xtindex >= 0)
          tranyz = nghtran(j0,index0) - ishft(xtindex, 8)
          do j = nstart, nstop
            jtran = tranyz + xtran(j-nstart+1,xtindex)
            jj = index1 + j - xtran(j-nstart+1,xtindex)*nucgrd1
            m1 = nlogrid(jj)
            m2 = nhigrid(jj)
            if (m2 >= m1) then
              do m = m1, m2
                numlist = numlist + 1
                atmlist(numlist) = m
                itran(numlist) = jtran
              end do
            end if
          end do
        end do
        ! End loop over grid pointers

        if (numlist > 0) then
          do k = ncell_lo, ncell_hi
            kk = k - ncell_lo + 1
            i = bckptr(k)
            xk = imagcrds(1,k)
            yk = imagcrds(2,k)
            zk = imagcrds(3,k)
            mygrdlist(k) = 1
            call pack_nb_list(kk, i, xk, yk, zk, imagcrds, cutoffsq, &
                              numlist, numpack, iac, ico, ntypes, ipairs, &
                              ifail, tranvec, belly, ibelly)
            if (ifail == 1) then
              listtot = maxnblst
              return
            end if
          end do
        end if
      end if
      ! End branch for whether to use the grid pointers approach

    end if
    ! Contingency for non-zero number of (numimg(index) > 0)
#  ifdef MPI
    if (balance) then
      gridpairs(index,1) = numpack - np0
      np0 = numpack
    end if
#  endif /* MPI */
  end do
  ! End loop over applicable indices (for this thread, or all
  ! such indices in the case of a serial run)

  listtot = numpack

  return

end subroutine get_nb_list 

!------------------------------------------------------------------------------
! grid_image: this subroutine grids the imaged atoms in the unit cell plus
!             neighbor cells while building the array of sorted image atoms.
!             The sorting should improve locality of reference.
!
! Arguments:
!   natom:    the number of atoms in the system
!   verbose:  verbosity level
!   periodic: integer that functions more like a logical, indicates whether the
!             system is periodic
!------------------------------------------------------------------------------
subroutine grid_image(natom, verbose, periodic)

  use constants, only : zero, half

  implicit none
  integer, intent(in) ::  natom, verbose, periodic
  _REAL_  f1, f2, f3, shft
  integer index, i, j, n
  integer i0, i1, i2, numimage
#include "extra.h"

  shft = half
  if (periodic == 0) then
    shft = zero
  end if
  numimage = 0
  do index = 1, nucgrd
    if (my_grids(index) == 1) then
      n = numatg(index)
      if (numimage + n > natom) then
        write(6,*)'natom = ',natom
        write(6,*)'numimage = ',numimage
        call sander_bomb('grid_image', 'num image atoms exceeds natom!!', &
                         '<ew_direct.f>grid_image() ')
      end if
      i0 = indoff(index)
      i1 = i0 + 1
      i2 = i0 + n
      nlogrid(index) = numimage + 1
      nhigrid(index) = numimage + n
      numimg(index) = n
      do i = i1, i2
        j = indatg(i)
        f1 = fraccrd(1,j) + shft
        f2 = fraccrd(2,j) + shft
        f3 = fraccrd(3,j) + shft
        numimage = numimage + 1
        imagcrds(1,numimage) = f1*ucell(1,1) + f2*ucell(1,2) + f3*ucell(1,3)
        imagcrds(2,numimage) = f1*ucell(2,1) + f2*ucell(2,2) + f3*ucell(2,3)
        imagcrds(3,numimage) = f1*ucell(3,1) + f2*ucell(3,2) + f3*ucell(3,3)
        bckptr(numimage) = j
      end do
    end if  ! (my_grids(index) == 1)
  end do  !  index = 1,nucgrd
  if (master) then
    if (verbose == 1) then
      write(6,*) 'grid_image: numimage = ', numimage
    end if
  end if

  return

end subroutine grid_image 

!------------------------------------------------------------------------------
! grid_pointers: list nearby cell line centers that need to be searched for
!                creating the non-bonded list.
!
!  Example:   2d grid of cells, X is cell of interest.  Here we are using 3
!             cells to the cutoff, so need 3 cells distant in each direction
!             We only point to the center of each string of 7 cells in the
!             x direction, and the nb list routine knows to check this cell
!             and three on each side of it in the +x and -x directions.
!             This routine will list the center of cell strings:
!
!         A                          B
!    ..................    ..................   and three more like layer B
!    ..................    ..................
!    .....ooo#ooo......    .....ooo#ooo......
!    .....ooo#ooo......    .....ooo#ooo......
!    .....ooo#ooo......    .....ooo#ooo......
!    ........Xooo......    .....ooo#ooo......
!    ..................    .....ooo#ooo......
!    ..................    .....ooo#ooo......
!    ..................    .....ooo#ooo......
!    ..................    ..................
!
!      cell X and its       Next layer ahead
!       neighbors o
!      center cell #
!    (ahead only search)
!
!    The pointer list contains the identity of all the # cells and the X cell
!
! Arguments:
!   nghbptr:  
!   nghtran:   
!   periodic:  flag to indicate that the system is periodic (integer works like
!              a logical)
!   nogrdptrs: flag to indicate that grid pointers are not in use
!------------------------------------------------------------------------------
subroutine grid_pointers(nghbptr, nghtran, periodic, nogrdptrs)

  implicit none
  integer nghbptr(0:maxnptrs,*), nghtran(maxnptrs,*)
  integer index0
  integer periodic

  integer i1, i2, i3, j2, j3, i2grd, i3grd, j3grd, j2grd, jj2
  integer j2strt, j2stop, jj2strt, jj2stop, j3stop
  integer index, num
  integer k3off
  integer xtindex_1, xtindex_2
  integer sizgrd12
  integer index0_min, index0_max

  logical nogrdptrs

  index0_min = maxnptrs
  index0_max = 0
  sizgrd12 = nucgrd1*nucgrd2
  my_grids(:) = 0
  num = 0

  ! Get forward pointers for unit cell grid.  Neighbor cells must be later
  ! in the subcell list to get unique atom - imageatom pairing

  ! Non-Periodic cell list pointers
  if (periodic == 0) then
    if (nogrdptrs) then
      do i1 = myindexlo, nucgrd1*nucgrd2*nucgrd3
        my_grids(i1) = 1
      end do
    else
      do i3 = 1,nucgrd3
        i3grd = sizgrd12*(i3-1)
        do i2 = 1,nucgrd2
          i2grd = i3grd + nucgrd1*(i2-1)
          do i1 = 1, nucgrd1
            index = i2grd+i1
            if (index < myindexlo .or. index > myindexhi) then
              goto 180
            end if
            my_grids(index) = 1
            index0 = index - myindexlo + 1
            index0_max = max(index0_max, index0)
            index0_min = min(index0_min, index0)

            ! Set up xtran index for this cell.  UNSTABLE until fixed:
            ! atoms must not drift to the edges.  Otherwise this
            ! simplified xtindex will not work.
            xtindex_1 = 1
            xtindex_2 = 1
            j2strt = -nghb2
            if (i2 + j2strt < 1) then
              j2strt = 1 - i2
            end if
            j2stop = nghb2
            if (i2 + j2stop > nucgrd2) then
              j2stop = nucgrd2 - i2
            end if
            jj2strt = -nghb1
            if (i1 + jj2strt < 1) then
              jj2strt = 1 - i1
            end if
            jj2stop = nghb1
            if (i1 + jj2stop > nucgrd1) then
              jj2stop = nucgrd1 - i1
            end if

            ! First, do the plane that this cell is in.  The
            ! first row of neighbors starts with the cell itself (Cell0).
            num = 1
            nghbptr(num, index0) = index

            ! No y or z translate, this cell HAS to be in uc
            nghtran(num,index0) = 5+ishft(1,8)
            do jj2 = 1, jj2stop
              my_grids(index+jj2) = 1
            end do

            ! Next rows of neighbors are the nghb2 rows of 2*nghb1+1 cells
            ! above.  Use the center cell of the row as the pointer (with
            ! the same x position as Cell0).  This way the pointer cell
            ! must be in the uc wrt the x direction.
            do j2 = 1, j2stop
              j2grd = index + j2*nucgrd1
              if (j2 + i2 <= nucgrd2) then
                num = num + 1
                nghtran(num,index0) = 5 + ishft(xtindex_2,8)
                nghbptr(num,index0) = j2grd
                do jj2 = jj2strt, jj2stop
                  my_grids(j2grd + jj2 - (nucgrd1 * &
                                          xtran(jj2+1+nghb1,xtindex_2))) = 1
                end do
              end if
            end do

            ! Now there are nghb3 planes of (2*nghb1+1)(2*nghb2+1) cells
            ! ahead that are good neighbors.
            j3stop = nghb3
            if (i3 + j3stop > nucgrd3) then
              j3stop = nucgrd3 - i3
            end if
            do j3 = 1, j3stop
              j3grd = j3 * sizgrd12
              do j2 = j2strt, j2stop
                num = num + 1
                j2grd = index + j3grd + j2*nucgrd1
                nghtran(num,index0) = 5 + ishft(xtindex_2,8)
                nghbptr(num,index0) = j2grd
                do jj2 = jj2strt, jj2stop
                  my_grids(j2grd+jj2) = 1
                end do
              end do
            end do
            nghbptr(0,index0) = num
            180 continue
          end do  !  i1 = 1,nucgrd1
        end do  !  i2 = 1,nucgrd2
      end do  !  i3 = 1,nucgrd3
    end if
    ! End contingency for no grid pointers

    numnptrs = maxnptrs
    return
  end if
  ! End contingency for non-periodic systems

  do i3 = 1,nucgrd3
    i3grd = sizgrd12*(i3 - 1)
    do i2 = 1,nucgrd2
      i2grd = i3grd + nucgrd1*(i2 - 1)
      do i1 = 1,nucgrd1
        index = i2grd + i1
        if ((index >= myindexlo) .and. (index <= myindexhi)) then
          my_grids(index) = 1
          index0 = index - myindexlo + 1
          index0_max = max(index0_max, index0)
          index0_min = min(index0_min, index0)

          ! Set up xtran index for this cell
          if (i1 + nghb1 > nucgrd1) then
            xtindex_1 = 1 + i1 + nghb1 - nucgrd1
            xtindex_2 = xtindex_1 + 2*nghb1
          else if(i1 - nghb1 < 1)then
            xtindex_1 = 1
            xtindex_2 = 2*nghb1 + 2 - i1
          else
            xtindex_1 = 1
            xtindex_2 = 1
          end if

          ! First do the plane that this cell is in.  The
          ! first row of neighbors starts with the cell itself (Cell0).
          num = 1
          nghbptr(num,index0) = index

          ! No y or z translate, this cell HAS to be in uc
          nghtran(num,index0) = 5 + ishft(xtindex_1,8)
          do jj2 = 1, nghb1
            my_grids(index + jj2 - (nucgrd1 * &
                                    xtran(jj2+1,xtindex_1))) = 1
          end do

          ! Next rows of neighbors are the nghb2 rows of 2*nghb1+1 cells
          ! above.  Use the center cell of the row as the pointer (with
          ! the same x position as Cell0).  This way the pointer cell
          ! must be in the uc wrt the x direction.
          do j2 = 1, nghb2
            num = num + 1
            j2grd = index + j2*nucgrd1
            if (j2+i2 <= nucgrd2) then
              nghtran(num,index0) = 5 + ishft(xtindex_2,8)
            else
              j2grd = j2grd - sizgrd12
              nghtran(num,index0) = 8 + ishft(xtindex_2,8)
            end if
            nghbptr(num,index0) = j2grd
            do jj2 = -nghb1, nghb1
              my_grids(j2grd + jj2 - (nucgrd1 * &
                                      xtran(jj2+1+nghb1,xtindex_2))) = 1
            end do
          end do

          ! Now there are nghb3 planes of (2*nghb1+1)(2*nghb2+1) cells
          ! ahead that are good neighbors.
          do j3 = 1, nghb3
            j3grd = j3 * sizgrd12
            if (j3 + i3 > nucgrd3) then
              j3grd = j3grd - sizgrd12*nucgrd3
              k3off = 9
            else
              k3off=0
            end if
            do j2 = -nghb2, nghb2
              num = num + 1
              j2grd = index + j3grd + j2*nucgrd1
              if (j2 + i2 > nucgrd2) then
                nghtran(num,index0) = 8 + k3off + ishft(xtindex_2,8)
                j2grd = j2grd - sizgrd12
              else if (j2 + i2 < 1) then
                nghtran(num,index0) = 2 + k3off + ishft(xtindex_2,8)
                j2grd = j2grd + sizgrd12
              else
                nghtran(num,index0) = 5 + k3off + ishft(xtindex_2,8)
              end if
              nghbptr(num,index0) = j2grd
              do jj2 = -nghb1, nghb1
                my_grids(j2grd + jj2 - (nucgrd1 * &
                                        xtran(jj2+1+nghb1,xtindex_2))) = 1
              end do
            end do
          end do
        endif
      end do
    end do
  end do
  ! End triply nested loop, scanning over the third, second, and first
  ! dimensions of the cell grid as per Fortran order.

  numnptrs = num
  return

end subroutine grid_pointers 

!------------------------------------------------------------------------------
! setup_grids: this routine checks on the necessary resources for the unit cell
!              and image cell grids used for short range particle pair
!              calculations.  It is assumed that unit cell setup has already
!              occurred.  This routine then checks the short range cutoff.  The
!              unit cell will be split into NUCGRD1 x NUCGRD2 x NUCGRD3
!              geometrically similar subcells of size dirlng(1)/NUCGRD1 by
!              dirlng(2)/NUCGRD2  by  dirlng(3)/NUCGRD3.  The short range
!              interactions will involve pairs in the subcell neighborhood of
!              +/- NGHB1 by +/- NGHB2 by +/- NGHB3 subcells, about any given
!              central subcell.  The distances between parallel faces of the
!              simulation cell are, respectively reclng(1), reclng(2) and
!              reclng(3).  Each subcell neighborhood is guaranteed to contain
!              all points within the minimum of (NGHB1/NUCGRD1)*reclng(1),
!              (NGHB2/NUCGRD2)*reclng(2),and (NGHB3/NUCGRD3)*1.0/reclng(3).
!              This minimum is taken to be the short range cutoff.
!
! Arguments:
!   periodic:  flag to indicate that the system is periodic (integer acts as a
!              logical)
!   nogrdptrs: flag to indicate that the grid-based list construction approach
!              is NOT in effect
!   verbose:   verbosity level
!------------------------------------------------------------------------------
subroutine setup_grids(periodic, nogrdptrs, verbose)
  implicit none
#include "extra.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#else
  integer mytaskid, numtasks
  parameter (mytaskid=0,numtasks=1)
#endif
  _REAL_ sizmaxhb
  integer periodic, verbose
  integer nghb0
  integer ngrd1, ngrd2, ngrd3
  logical :: nogrdptrs

  _REAL_ dc1, dc2, dc3, cut

  ! Maximum allowed bond length for Hydrogen to heavy atoms
  parameter(sizmaxhb = 1.34d0)

  nghb0 = nghb
  nghb1 = nghb0
  nghb2 = nghb0
  nghb3 = nghb0
  dc1 = cutlist / nghb1
  dc2 = cutlist / nghb2
  dc3 = cutlist / nghb3
  nucgrd1 = max(1, int(reclng(1)/dc1))
  nucgrd2 = max(1, int(reclng(2)/dc2))
  nucgrd3 = max(1, int(reclng(3)/dc3))

  if (periodic == 0) then

    ! If system has too few cells for the pointer method
    ! to be efficient, use the no-grid-pointer system
    ! where all forward cells are checked for pair atoms
    ! rather than just the forward cells on the grid-pointer-list
    nogrdptrs = nogrdptrs .or. ((nucgrd1 <= 2*nghb+2) .or. &
                                (nucgrd2 <= 2*nghb+2) .or. &
                                (nucgrd3 <= 2*nghb+2))
    if (verbose >= 3) then
      write(6,*) " List using nogrdptrs ", nogrdptrs
    end if

    ! Make the subcells about 3 A in size now so that there
    ! are more subcells, easier to load balance.
    if (nogrdptrs) then
      if (dirlng(1)/nucgrd1 > 3) then
        ngrd1 = int(dirlng(1)/3)
        ngrd2 = int(dirlng(2)/3)
        ngrd3 = int(dirlng(3)/3)
        if (verbose >= 1) then
          write(6,'("|    New Grids set up for nogrdptrs ")')
          write(6, '(5X,a,/,5X,i9,1X,i9,1X,i9)') &
                'Number of grids per unit cell x,y,z:', &
                ngrd1, ngrd2, ngrd3
        end if
        nucgrd1 = ngrd1
        nucgrd2 = ngrd2
        nucgrd3 = ngrd3
      end if
    end if
  end if
  ! Check the short range cutoff:
   
  dc1 = reclng(1)/nucgrd1
  dc2 = reclng(2)/nucgrd2
  dc3 = reclng(3)/nucgrd3
  cut = nghb1*dc1
  if (nghb2*dc2 < cut) then
    cut = nghb2*dc2
  end if
  if (nghb3*dc3 < cut) then
    cut = nghb3*dc3
  end if
  if (nogrdptrs) then
    cut = cutlist
  end if
#ifdef MPI
  if (master) then
#endif
    if (verbose >= 1) then
      write(6, '(5X,a,/,5X,i9,1X,i9,1X,i9)') &
            'Number of grids per unit cell in each dimension:', &
            nucgrd1, nucgrd2, nucgrd3
      write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
            'Unit cell edge lengths in each dimension:', &
            dirlng(1), dirlng(2), dirlng(3)
      write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
            'Distance between parallel faces of unit cell:', &
            reclng(1), reclng(2), reclng(3)
      write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
            'Distance between faces of short range grid subcells:', &
            dc1, dc2, dc3
      write(6, '(5X,a,F9.3)') &
            'Resulting cutoff from subcell neighborhoods is ', cut
    end if
#ifdef MPI
  end if
#endif
  if (cut < cutlist) then
    call sander_bomb('setup_grids', &
                     'Resulting cutoff is too small for your lower limit', ' ')
  end if

  return

end subroutine setup_grids 

!------------------------------------------------------------------------------
! setup_grid_sizes: 
!
! Arguments:
!   periodic:  flag to indicate that the system is periodic (integer acts as a
!              logical)
!   nogrdptrs: flag to indicate that the grid-based list construction approach
!              is NOT in effect
!   verbose:   verbosity level
!------------------------------------------------------------------------------
subroutine setup_grid_sizes(periodic, nogrdptrs, verbose)

  implicit none
#include "extra.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#endif
  _REAL_ sizmaxhb
  integer periodic, verbose
  integer nghb0
  integer ngrd1, ngrd2, ngrd3
  logical nogrdptrs
  _REAL_ dc1, dc2, dc3, cut

  ! Maximum allowed bond length for Hydrogen to heavy atoms
  parameter(sizmaxhb = 1.34d0)

  nghb0 = nghb
  nghb1 = nghb0
  nghb2 = nghb0
  nghb3 = nghb0
  dc1 = (cutlist + sizmaxhb) / nghb1
  dc2 = (cutlist + sizmaxhb) / nghb2
  dc3 = (cutlist + sizmaxhb) / nghb3
  nucgrd1 = max(1, int(reclng(1)/dc1))
  nucgrd2 = max(1, int(reclng(2)/dc2))
  nucgrd3 = max(1, int(reclng(3)/dc3))

  if (periodic == 0) then

    ! If system has too few cells for the pointer method
    ! to be efficient, use the no-grid-pointer system
    ! where all forward cells are checked for pair atoms
    ! rather than just the forward cells on the grid-pointer-list
    nogrdptrs = nogrdptrs .or. ((nucgrd1 <= 2*nghb+2) .or. &
                                (nucgrd1 <= 2*nghb+2) .or. &
                                (nucgrd1 <= 2*nghb+2))
    if (verbose >= 3) then
      write(6,*) " List using nogrdptrs ", nogrdptrs
    end if

    ! Make the subcells about 3 A in size now so that there
    ! are more subcells, easier to load balance.
    if (nogrdptrs) then
      if (dirlng(1)/nucgrd1 > 3) then
        ngrd1 = int(dirlng(1)/3)
        ngrd2 = int(dirlng(2)/3)
        ngrd3 = int(dirlng(3)/3)
        if (verbose >= 1) then
          write(6,'("|    New Grids set up for nogrdptrs ")')
          write(6, '(5X,a,/,5X,i9,1X,i9,1X,i9)') &
                'Number of grids per unit cell x,y,z:', ngrd1, ngrd2, ngrd3
        end if
        nucgrd1 = ngrd1
        nucgrd2 = ngrd2
        nucgrd3 = ngrd3
      end if
    end if
  end if
   
  ! Check the short range cutoff:
  dc1 = reclng(1) / nucgrd1
  dc2 = reclng(2) / nucgrd2
  dc3 = reclng(3) / nucgrd3
  cut = nghb1 * dc1
  if (nghb2*dc2 < cut) then
    cut = nghb2*dc2
  end if
  if (nghb3*dc3 < cut) then
    cut = nghb3*dc3
  end if
  if (nogrdptrs) then
    cut = cutlist
  end if
#ifdef MPI
  if (master) then
#endif
    if (verbose >= 1) then
      write(6, '(5X,a,/,5X,i9,1X,i9,1X,i9)') &
            'Number of grids per unit cell in each dimension:', &
            nucgrd1, nucgrd2, nucgrd3
      write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
            'Unit cell edge lengths in each dimension:', &
            dirlng(1), dirlng(2), dirlng(3)
      write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
            'Distance between parallel faces of unit cell:', &
            reclng(1), reclng(2), reclng(3)
      write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
            'Distance between faces of short range grid subcells:', &
            dc1, dc2, dc3
      write(6, '(5X,a,F9.3)') &
            'Resulting cutoff from subcell neighborhoods is ', cut
    end if
#ifdef MPI
  end if
#endif
  if ( cut < cutlist )then
    call sander_bomb('setup_grids', &
                     'Resulting cutoff is too small for your lower limit', ' ')
  end if

  return

end subroutine setup_grid_sizes

!------------------------------------------------------------------------------
! grid_ucell: this routine grids the mapped atoms in the unit cell
!             into the NUCGRD1 x NUCGRD2 x NUCGRD3 subcells according
!             to their fractional coordinates.
!
! Arguments:
!   natom:     the number of atoms in the system
!   periodic:  flag to indicate whether the system is periodic
!------------------------------------------------------------------------------
subroutine grid_ucell(natom, periodic)

  use constants, only : half

  implicit none
  integer natom, periodic
  integer i, j, i1, i2, i3, index

  _REAL_ shft

  shft = half
  if (periodic == 0) then
    shft = 0.0
  end if

  do index = 1, nucgrd
    numatg(index) = 0
  end do

  ! Find out which ucgrd subcell each atom is in.
  ! atmcell(i) is the ucgrd subcell that contains atom i.
  ! numatg(index) is the number of atoms in the index ucgrd.
  do i = 1, natom
    i1 = int((fraccrd(1,i) + shft) * dble(nucgrd1)) + 1
    i2 = int((fraccrd(2,i) + shft) * dble(nucgrd2)) + 1
    i3 = int((fraccrd(3,i) + shft) * dble(nucgrd3)) + 1
    index = nucgrd1*nucgrd2*(i3-1)+nucgrd1*(i2-1)+i1
    atmcell(i) = index
    numatg(index) = numatg(index) + 1
  end do

  ! Find the offset of the starting atoms for each ucgrd
  ! subcell.  Zero the numatg()entries in the process.
  indoff(1) = 0
  do i = 2, nucgrd
    indoff(i) = indoff(i-1) + numatg(i-1)
    numatg(i-1) = 0
  end do
  numatg(nucgrd) = 0

  ! Fill indatg() as a list of atoms in each subcell such that
  ! the list of atoms in subcell 1 are at the beginning,
  ! subcell 2 list is right after that (starting at indoff(2)+1)
  do i = 1, natom
    index = atmcell(i)
    numatg(index) = numatg(index) + 1
    j = numatg(index) + indoff(index)
    indatg(j) = i
  end do

  return

end subroutine grid_ucell 

!------------------------------------------------------------------------------
! pack_nb_list: pack the non-bonded list with likely interacting atom pairs.
!
! Arguments:
!   kk:       the atom for which the list is being packed, the kth atom from
!             a series within the packed array of atoms populating each grid
!             cell
!   i:        the current atom to build a pairlist for
!   xk:
!   yk:       imaged coordinates of the atom in x, y, and z
!   zk:
!   imagcrds: imaged coordinates of all atoms in the system, rearranged as they
!             appear in each cell
!   cutoffsq:
!   numlist:
!   numpack:
!   iac:
!   ico:
!   ntypes:
!   ipairs:
!   ifail:
!   tranvec:
!   belly:
!   ibelly:
!------------------------------------------------------------------------------
subroutine pack_nb_list(kk, i, xk, yk, zk, imagcrds, cutoffsq, numlist, &
                        numpack, iac, ico, ntypes, ipairs, ifail, tranvec, &
                        belly, ibelly)
#ifdef LES
  use les_data, only : cnum, lestyp
#endif
#ifdef MPI /* SOFT CORE */
  use softcore, only : nsc, ifsc
#endif
  use constants, only : zero

  implicit none

  integer numpack
  integer kk, i, numlist, iac(*), ico(*), ntypes, ifail 
  _REAL_ xk, yk, zk, imagcrds(3,*), cutoffsq
  integer ipairs(*), ibelly(*)
  logical belly, deadi, deadik
  integer num

  !Local variables for QMMM
  integer qm_temp_count2

#ifdef MPI
#  include "parallel.h"
#endif

  _REAL_ dx, dy, dz, r2
  integer iaci, ic, index, k, lpr, lps, lhb, m, n, npr

  integer jtran
  _REAL_ x_tran, y_tran, z_tran, tranvec(3,*)

  num = 0
  m = maskptr(i)
  lpr = 0
  lps = 0
  lhb = 0
  iaci = ntypes*(iac(i) - 1)
   
  do n = 1, nummask(i)
    k = lstmask(m+n)
    exclude(k) = i
  end do

  deadi = .false.
  deadi = ((ibelly(i) == 0) .and. belly)

#ifdef MPI /* SOFT CORE */
  softcore_on: if (ifsc .ne. 0) then

    ! Softcore potential in use, check which list each atom goes into.
    ! The first thing to know is whether i is a softcore atom.
    check_softcore: if (nsc(i) == 0) then
      do m = kk+1,numlist            
        jtran = itran(m)
        x_tran = tranvec(1,itran(m))
        y_tran = tranvec(2,itran(m))
        z_tran = tranvec(3,itran(m))
        n = atmlist(m)
        k = bckptr(n)
        deadik = ((ibelly(k) == 0) .and. deadi)
        if ((exclude(k) /= i) .and. .not. deadik) then
          dx = imagcrds(1,n) - xk + x_tran
          dy = imagcrds(2,n) - yk + y_tran
          dz = imagcrds(3,n) - zk + z_tran
          r2 = dx*dx + dy*dy + dz*dz
          if (r2 < cutoffsq) then
#  ifdef TEST_CLOSEATOMS
            if (r2 < .5) then
              write(6,190) i, k, r2
            end if
190         format("<pack_nb_list> Atoms close:", 2i8, "  r**2 ", f10.6)
#  endif
            mygrdlist(n) = ior(mygrdlist(n), 1)
            num = num + 1
            index = iaci + iac(k)
            ic = ico(index)
            if (ic > 0) then
              if (nsc(k) == 0) then
                lpr = lpr+1
                iwa(lpr) = ior(n, ishft(itran(m), 27)) ! bitwise optimization
              else
                lps = lps+1
                iws(lps) = ior(n, ishft(itran(m), 27)) ! bitwise optimization
              end if
            else
              if (nsc(k) == 0) then
                lhb = lhb+1
                iwh(lhb) = ior(n, ishft(itran(m), 27)) ! bitwise optimization
              else
                lps = lps+1
                iws(lps) = ior(n, ishft(itran(m), 27)) ! bitwise optimization
              end if
            end if
          end if
        end if
      end do  !  m = kk+1,numlist
    else check_softcore ! atom i IS a softcore atom
      do m = kk+1, numlist            
        jtran = itran(m)
        x_tran = tranvec(1,itran(m))
        y_tran = tranvec(2,itran(m))
        z_tran = tranvec(3,itran(m))
        n = atmlist(m)
        k = bckptr(n)
        deadik = ((ibelly(k) == 0) .and. deadi)
        if ((exclude(k) /= i) .and. .not. deadik) then
          dx = imagcrds(1,n) - xk + x_tran
          dy = imagcrds(2,n) - yk + y_tran
          dz = imagcrds(3,n) - zk + z_tran
          r2 = dx*dx + dy*dy + dz*dz
          if (r2 < cutoffsq) then
            mygrdlist(n) = ior(mygrdlist(n), 1)
            num = num + 1
            index = iaci + iac(k)
            ic = ico(index)
            lps = lps + 1
            iws(lps) = ior(n, ishft(itran(m), 27)) ! bitwise optimization
          end if
        end if
      end do
    end if check_softcore
  else softcore_on
# endif /* MPI, SOFTCORE */

  ! Softcore potential not on, build list the normal way
  do m = kk+1, numlist
    jtran = itran(m)
    x_tran = tranvec(1,itran(m))
    y_tran = tranvec(2,itran(m))
    z_tran = tranvec(3,itran(m))
    n = atmlist(m)
    k = bckptr(n)
    deadik = ((ibelly(k) == 0) .and. deadi)
    if ((exclude(k) /= i) .and. .not. deadik) then
      dx = imagcrds(1,n) - xk + x_tran
      dy = imagcrds(2,n) - yk + y_tran
      dz = imagcrds(3,n) - zk + z_tran
      r2 = dx*dx + dy*dy + dz*dz
      if (r2 < cutoffsq) then
#  ifdef TEST_CLOSEATOMS
        if (r2 < .5) then
          write(6,190) i, k, r2
        end if
190     format("<pack_nb_list> Atoms close:", 2i8, "  r**2 ", f10.6)
#  endif
        mygrdlist(n) = ior(mygrdlist(n), 1)
        num = num + 1
        index = iaci + iac(k)
        ic = ico(index)
        if (ic >= 0) then
          lpr = lpr+1
          iwa(lpr) = ior(n, ishft(itran(m), 27)) ! bitwise optimization
        else
          lhb = lhb+1
          iwh(lhb) = ior(n, ishft(itran(m), 27)) ! bitwise optimization
        end if
      end if
    end if
  end do
  ! End loop for normal list build

# ifdef MPI /* SOFT CORE */
  ! Having this in pre-processor directives changes
  ! the indentation of the loop above.
  end if softcore_on
# endif 
  if (num + numpack > maxnblst) then
#ifdef MPI
    write(6, '(/,a,i12,i12,a,i12,a,i4)') ' * NB pairs ', num, numpack, &
          ' exceeds capacity (', maxnblst, ')', mytaskid
#else
    write(6, '(/,a,i12,i12,a,i12,a)') ' * NB pairs ', num, numpack, &
           ' exceeds capacity (', maxnblst, ')'
#endif
    ifail=1
    return
  end if

  ! Now put all pairs into iwa
  do k = 1, lhb
    iwa(lpr+k) = iwh(k)
  end do
  npr = lpr + lhb
#ifdef MPI /* SOFT CORE */
  iwa(npr+1:npr+lps) = iws(1:lps)
  npr = npr + lps
  numsc(i) = lps
#endif

  numvdw(i) = lpr
  numhbnd(i) = lhb
  do k = 1, npr
    numpack = numpack + 1
    ipairs(numpack) = iwa(k)
  end do

  return

end subroutine pack_nb_list 

!------------------------------------------------------------------------------
! pack_nb_nogrdptrs: pack the non-bonded list without grid pointers.
!
! Arguments:
!------------------------------------------------------------------------------
subroutine pack_nb_nogrdptrs(kk, i, xk, yk, zk, imagcrds, cutoffsq, numlist, &
                             numpack, iac, ico, ntypes, ipairs, ifail, belly, &
                             ibelly)

#ifdef MPI /* SOFT CORE */
  use softcore, only : nsc, ifsc
#endif

  implicit none
  integer numpack
  integer kk, i, numlist, iac(*), ico(*), ntypes, ifail 
  _REAL_ xk, yk, zk, imagcrds(3,*), cutoffsq
  integer ipairs(*), ibelly(*)
  logical belly, deadi, deadik
  integer jtran, num

  ! Local variables for QMMM
  integer qm_temp_count2

#ifdef MPI
#  include "parallel.h"
#endif

  _REAL_ dx, dy, dz, r2
  integer iaci, ic, index, k, lpr, lps, lhb, m, n, npr

  num = 0
  m = maskptr(i)
  lpr = 0
  lps = 0
  lhb = 0
  iaci = ntypes*(iac(i)-1)
  do n = 1, nummask(i)
    k = lstmask(m+n)
    exclude(k) = i
  end do

  deadi = .false.
  deadi = ((ibelly(i) == 0) .and. belly)
#ifdef MPI /* SOFT CORE */
  softcore_on: if (ifsc /= 0) then
    check_softcore: if (nsc(i) == 0) then ! atom i is NOT a softcore atom
      do m = kk+1,numlist
        jtran = 5
        n = atmlist(m)
        k = bckptr(n)
        deadik = ((ibelly(k) == 0) .and. deadi)
        if ((exclude(k) /= i) .and. .not. deadik) then
          dx = imagcrds(1,n) - xk
          dy = imagcrds(2,n) - yk
          dz = imagcrds(3,n) - zk
          r2 = dx*dx + dy*dy + dz*dz
          if (r2 < cutoffsq) then
# ifdef TEST_CLOSEATOMS
            if (r2 < .5) then
              write(6,190) i, k, r2
            end if
190         format("<pack_nb_list> Atoms close:", 2i8, "  r**2 ", f10.6)
# endif
            if (mygrdlist(n) == 0) then
              mygrdlist(n) = 1
            end if
            num = num + 1
            index = iaci + iac(k)
            ic = ico(index)
            if (ic >= 0) then
              if (nsc(k) == 0) then
                lpr = lpr + 1
 
                ! Bitwise optimization
                iwa(lpr) = ior(n, ishft(5, 27))
              else
                lps = lps + 1
 
                ! Bitwise optimization
                iws(lps) = ior(n, ishft(itran(m), 27))
              end if
            else
              lhb = lhb + 1

              ! bitwise optimization
              iwh(lhb) = ior(n, ishft(5, 27))
            end if
          end if
        end if
      end do  !  m = kk+1,numlist
    else check_softcore ! atom i IS a softcore atom
      do m = kk+1, numlist
        jtran = 5
        n = atmlist(m)
        k = bckptr(n)
        deadik = ((ibelly(k) == 0) .and. deadi)
        if ((exclude(k) /= i) .and. .not. deadik) then
          dx = imagcrds(1,n) - xk
          dy = imagcrds(2,n) - yk
          dz = imagcrds(3,n) - zk
          r2 = dx*dx + dy*dy + dz*dz
          if (r2 < cutoffsq) then
            if (mygrdlist(n) == 0) then
              mygrdlist(n) = 1
            end if
            num = num + 1
            index = iaci + iac(k)
            ic = ico(index)
            if (ic >= 0) then
              lps = lps + 1
              iws(lps) = ior(n, ishft(5, 27)) ! bitwise optimization
            else ! unlikely case atom is softcore and 10-12 flagged
              lhb = lhb + 1
              iwh(lhb) = ior(n, ishft(5,27)) ! bitwise optimization
            end if
          end if
        end if
      end do  !  m = kk+1,numlist
    end if check_softcore
  else softcore_on ! softcore is not on, build list the usual way
#endif /* SOFT CORE  */
    do m = kk+1, numlist
      jtran = 5
      n = atmlist(m)
      k = bckptr(n)
      deadik = ((ibelly(k) == 0) .and. deadi)
      if ((exclude(k) /= i) .and. .not. deadik) then
        dx = imagcrds(1,n) - xk
        dy = imagcrds(2,n) - yk
        dz = imagcrds(3,n) - zk
        r2 = dx*dx + dy*dy + dz*dz
        if (r2 < cutoffsq) then
# ifdef TEST_CLOSEATOMS
          if (r2 < .5) then
            write(6,190) i, k, r2
          end if
190       format("<pack_nb_list> Atoms close:", 2i8, "  r**2 ", f10.6)
# endif
          if (mygrdlist(n) == 0) then
            mygrdlist(n) = 1
          end if
          num = num + 1
          index = iaci + iac(k)
          ic = ico(index)
          if (ic >= 0) then
            lpr = lpr+1

            ! Bitwise optimization
            iwa(lpr) = ior(n, ishft(5, 27))
          else
            lhb = lhb+1
            iwh(lhb) = ior(n, ishft(5, 27)) ! bitwise optimization
          end if
        end if
      end if
    end do  !  m = kk+1,numlist
#ifdef MPI /* SOFT CORE */
  end if softcore_on
#endif

  if (num + numpack > maxnblst) then
#ifdef MPI
    write(6, '(/,a,i12,i12,a,i12,a,i4)') ' * NB pairs ', num, numpack, &
          ' exceeds capacity (', maxnblst, ')', mytaskid
#else
    write(6, '(/,a,i12,i12,a,i12,a)') ' * NB pairs ', num, numpack, &
          ' exceeds capacity (', maxnblst, ')'
#endif
    ifail = 1
    return
  end if

  ! Now put all pairs into iwa
  do k = 1, lhb
    iwa(lpr+k) = iwh(k)
  end do
  npr = lpr + lhb

#ifdef MPI /* SOFT CORE */
  iwa(npr+1:npr+lps) = iws(1:lps)
  npr = npr + lps
  numsc(i) = lps
#endif
  numvdw(i) = lpr
  numhbnd(i) = lhb
  do k = 1, npr
    numpack = numpack + 1
    ipairs(numpack) = iwa(k)
  end do
  return
end subroutine pack_nb_nogrdptrs 

!------------------------------------------------------------------------------
! fix_grid_balance: re-balance the grid by divvying cells among processors.
!
! Arguments:
!   gridpairs:
!   nucgrd:       three element array holding the dimensions of the cell grid
!   numtaks:      the number of multi-processor tasks to divide amongst the
!                 processors
!   mytaskid:     task index that the current processor sees
!   listdiffmax:  sends back the threshold at which the grid gets rebalanced
!------------------------------------------------------------------------------
subroutine fix_grid_balance(gridpairs, nucgrd, numtasks, mytaskid, listdiffmax)

  implicit none

  integer gridpairs(*)
  integer nucgrd, listtot, numtasks, mytaskid, listdiffmax
  integer i, n, lsum, i0, ishare

  n = 0
  lsum = 0
  myindexlo = 0
  myindexhi = 0
  listtot = 0
  do i = 1, nucgrd
    listtot = listtot + gridpairs(i)
  end do
  i0 = 1
  ishare = listtot / numtasks
   
  ! Set trigger for new balance as 10 percent of list size  
  listdiffmax = max(1000, ishare/10)
   
  ! Give cells to each processor .lt. mytaskid until
  ! they have ishare, then give the next
  ! cells to this pe till it has its share, then return   
  do i = 1, nucgrd
    lsum = lsum + gridpairs(i)
    if (lsum >= ishare) then
      if (n == mytaskid) then
        myindexlo = i0
        myindexhi = i
        return
      end if
      lsum = 0
      i0 = i + 1
      n = n + 1
    end if
  end do
  if (n == mytaskid) then
    myindexlo = i0
    myindexhi = nucgrd
    return
  end if
  myindexlo = nucgrd

  return

end subroutine fix_grid_balance 

!------------------------------------------------------------------------------
! adjust_imagcrds: needed in case coordinates are bing wrapped in the box.
!                  The code can be run without wrapping (since it is easier to
!                  analyze the results) but it should work equally well,
!                  either way.
!
! Arguments:
!   crd:      the coordinates of all atoms
!   natom:    the number of atoms in the system
!------------------------------------------------------------------------------
subroutine adjust_imagcrds(crd, natom)

   use constants, only : half
   implicit none
   integer, intent(in) :: natom
   _REAL_, dimension(3,natom), intent(in) :: crd

   integer i, j
   _REAL_ f1, f2, f3
   _REAL_ anint
#  include "parallel.h"
#  include "ew_cntrl.h"
#  include "box.h"
  ! For nonperiodic, there is no imaging so frac coords are not
  ! used and no wrapping is done, but the imgcrds need to be
  ! in the right order.
  if (periodic == 0) then
    do i = 1, natom
      if (mygrdlist(i) == 1) then
        j = bckptr(i)
        imagcrds(1,i) = crd(1,j)+xbox0
        imagcrds(2,i) = crd(2,j)+ybox0
        imagcrds(3,i) = crd(3,j)+zbox0
      end if
    end do
    return
  end if
   
  !   ---- Periodic systems
  !        find change since last step:
  i = 0
  do while(i < natom)
    i = i + 1
    if (mygrdlist(i) == 1) then
      j = bckptr(i)
      dfrac(1,j) = fraccrd(1,j) - savfrac(1,j)
      dfrac(2,j) = fraccrd(2,j) - savfrac(2,j)
      dfrac(3,j) = fraccrd(3,j) - savfrac(3,j)
      dfrac(1,j) = dfrac(1,j) - anint(dfrac(1,j))
      dfrac(2,j) = dfrac(2,j) - anint(dfrac(2,j))
      dfrac(3,j) = dfrac(3,j) - anint(dfrac(3,j))
      savfrac(1,j) = savfrac(1,j) + dfrac(1,j)
      savfrac(2,j) = savfrac(2,j) + dfrac(2,j)
      savfrac(3,j) = savfrac(3,j) + dfrac(3,j)
      f1 = savfrac(1,j)+half
      f2 = savfrac(2,j)+half
      f3 = savfrac(3,j)+half
      if (ifbox == 1) then
        imagcrds(1,i) = f1*ucell(1,1)
        imagcrds(2,i) = f2*ucell(2,2)
        imagcrds(3,i) = f3*ucell(3,3)
      else
        imagcrds(1,i) = f1*ucell(1,1) + f2*ucell(1,2) + f3*ucell(1,3)
        imagcrds(2,i) = f1*ucell(2,1) + f2*ucell(2,2) + f3*ucell(2,3)
        imagcrds(3,i) = f1*ucell(3,1) + f2*ucell(3,2) + f3*ucell(3,3)
      end if
    end if
  end do
   
  return

end subroutine adjust_imagcrds 

!------------------------------------------------------------------------------
! map_coords: this routine takes the cartesian coordinates of the atoms and
!             maps them to fractional coords.  Everything is done atom by atom,
!             i.e. translations are not done by residues or by molecules.  The
!             fractional coordinates for atoms are obtained from the dot
!             products with the reciprocal lattice vectors.
!
! Arguments:
!   crd:      the coordinates of all atoms
!   natom:    the number of atoms in the system
!   recip:    the reciprocal space lattice vectors (transformation matrix to
!             take coordinates into box space)
!------------------------------------------------------------------------------
subroutine map_coords(crd,natom,recip)
   
  use constants, only : half, zero, one

  implicit none
#include "ew_cntrl.h"
#include "extra.h"
#include "box.h"

  integer natom
  _REAL_ crd(3,natom)
  _REAL_ recip(3,3)
  integer i
  _REAL_ anint, fracmax, fracmin
   
  fracmax = half
  fracmin = half
  if (periodic == 0) then
    do i = 1, natom
      fraccrd(1,i) = (crd(1,i)+xbox0)*recip(1,1)
      fraccrd(2,i) = (crd(2,i)+ybox0)*recip(2,2)
      fraccrd(3,i) = (crd(3,i)+zbox0)*recip(3,3)
      fracmax = max(fracmax, fraccrd(1,i), fraccrd(2,i), fraccrd(3,i))
      fracmin = min(fracmin, fraccrd(1,i), fraccrd(2,i), fraccrd(3,i))
    end do
    if (fracmax > one .or. fracmin < zero) then
      write(6,*) "Frac coord min, max: ", fracmin, fracmax
      write(6,*) "The system has extended beyond "
      write(6,*) "    the extent of the virtual box."
      write(6,*) "Restarting sander will recalculate"
      write(6,*) "   a new virtual box with 30 Angstroms"
      write(6,*) "   extra on each side, if there is a"
      write(6,*) "   restart file for this configuration."
      call sander_bomb('Routine: map_coords (ew_force.f)', &
                       'Atom out of bounds. If a restart has been written,', &
                       'restarting should resolve the error')
    end if
    return
  end if

  ! Get fractional coordinates: ifbox == 1 indicates
  ! an orthonormal simulation box.
  if (ifbox == 1) then
    do i = 1, natom
      fraccrd(1,i) = crd(1,i) * recip(1,1)
      fraccrd(2,i) = crd(2,i) * recip(2,2)
      fraccrd(3,i) = crd(3,i) * recip(3,3)
    end do
  else
    do i = 1, natom
      fraccrd(1,i) = crd(1,i)*recip(1,1) + crd(2,i)*recip(2,1) + &
                     crd(3,i)*recip(3,1)
      fraccrd(2,i) = crd(1,i)*recip(1,2) + crd(2,i)*recip(2,2) + &
                     crd(3,i)*recip(3,2)
      fraccrd(3,i) = crd(1,i)*recip(1,3) + crd(2,i)*recip(2,3) + &
                     crd(3,i)*recip(3,3)
    end do
  end if

  ! Check if system has gone out of box for nonperiodic
  ! simulations with finite cutoff.
  if (periodic == 0) then
    if (.not. nocutoff) then
      boxbad = .false.
      do i = 1, natom
        if (anint(fraccrd(1,i)-half) /= zero .or. &
            anint(fraccrd(2,i)-half) /= zero .or. &
            anint(fraccrd(3,i)-half) /= zero) then
          boxbad = .true.
        end if
      end do
    end if
    if (boxbad) then
      write(6,*) "**********BOX IS BAD****************"
    end if
  else

    !      --- map them inside, if periodic:
    do i = 1, natom
      fraccrd(1,i) = fraccrd(1,i) - anint(fraccrd(1,i))
      fraccrd(2,i) = fraccrd(2,i) - anint(fraccrd(2,i))
      fraccrd(3,i) = fraccrd(3,i) - anint(fraccrd(3,i))
    end do
  end if

  return

end subroutine map_coords 

!------------------------------------------------------------------------------
! save_frac_crds: needed in order to wrap coordinates inside the box.  Note,
!                 td actually runs ewald without wrapping since it is easier to
!                 analyze but it should work equally well, either way.
!
! Arguments:
!   natom:     the number of atoms in the system
!------------------------------------------------------------------------------
subroutine save_frac_crds(natom)

  implicit none

  integer natom
  ! These are actually two-dimensional (3,natom), but to enable
  ! vectorization on IA32 SSE platforms they are treated as
  ! one-dimensional; this may also improve software pipelining !

  integer i
  do i = 1, natom
    savfrac(1,i) = fraccrd(1,i)
    savfrac(2,i) = fraccrd(2,i)
    savfrac(3,i) = fraccrd(3,i)
  end do

  return

end subroutine save_frac_crds 

!------------------------------------------------------------------------------
! save_crds: save the atomic coordinates that were used to build the list.
!            This is needed for the skin test for buffered pairlists, where
!            the current coordinates are compared with the saved coordinates
!            relative to the skin criterion; see check_skin.
!
! Arguments:
!   natom:     the number of atoms in the system
!   crd:       the coordinates of all atoms in the system
!------------------------------------------------------------------------------
subroutine save_crds(natom, crd)

  implicit none

  integer natom
  _REAL_ crd(3,natom)

  !     These are actually two-dimensional (3,natom), but to enable
  !     vectorization on IA32 SSE platforms they are treated as
  !     one-dimensional; this may also improve software pipelining.
  integer i
  do i = 1,natom
    savcrd(1,i) = crd(1,i)
    savcrd(2,i) = crd(2,i)
    savcrd(3,i) = crd(3,i)
  end do

  return

end subroutine save_crds 

!------------------------------------------------------------------------------
! check_skin: Check if any atom has moved more than half the skin distance,
!             which is half the nbskin added to the cutoff in generating the
!             verlet list; in which case a list build is flagged.  An obvious
!             parallel implementation, in which each processor checks its atoms
!             and communicates the results to all other processors, produced on
!             Linux clusters large losses with large numbers of processors.
!             The separate sections based on nbtell exist for improved
!             performance; computing the list update info has a small
!             performance cost.
!
! Arguments:
!   crd:             current coordinates of all atoms (input only)
!   do_list_update:  true if a list update is needed, (output)
!------------------------------------------------------------------------------
subroutine check_skin(crd,do_list_update)

  use constants, only : zero, fourth
  implicit none
#include "extra.h"
#include "../include/memory.h"

  _REAL_, intent(in) :: crd(3,natom)
  logical, intent(out) :: do_list_update

  integer first_atom
  integer last_atom
  integer i
  integer nmoved_atoms     ! total number of atoms triggering a list update
  _REAL_ dx, dy, dz, dis2
  _REAL_ maxdis2

  steps_since_list_build = steps_since_list_build + 1
  first_atom = 1
  last_atom  = natom
  maxdis2    = zero
  if (nbtell == 0) then

    ! List update info not requested
    do i = first_atom, last_atom
      dx = crd(1,i) - savcrd(1,i)
      dy = crd(2,i) - savcrd(2,i)
      dz = crd(3,i) - savcrd(3,i)
      dis2 = dx**2 + dy**2 + dz**2
      maxdis2 = max(dis2,maxdis2)
    end do
    do_list_update = maxdis2 > fourth * skinnb * skinnb
  else

    ! List update info requested
    nmoved_atoms = 0
    do i = first_atom, last_atom
      dx = crd(1,i) - savcrd(1,i)
      dy = crd(2,i) - savcrd(2,i)
      dz = crd(3,i) - savcrd(3,i)
      dis2 = dx**2 + dy**2 + dz**2
      if (dis2 > fourth*skinnb*skinnb) then
        maxdis2 = max(dis2,maxdis2)
        nmoved_atoms = nmoved_atoms + 1
      end if
    end do
    do_list_update = (nmoved_atoms > 0)
    if (do_list_update) then
      if (master) then
        write(6, '(1x,A,I7,/,1x,A,F8.4,I7)') &
             'List Build Triggered: Number of atoms triggering = ', &
             nmoved_atoms, ' Maximum distance moved = ', sqrt(maxdis2)
      end if
    end if
  end if

  return

end subroutine check_skin

!------------------------------------------------------------------------------
! fill_tranvec: 
!------------------------------------------------------------------------------
subroutine fill_tranvec()

  implicit none
  integer iv, i0, i1, i2, i3
   
  iv=0
  do i3 = 0, 1
    do i2 = -1, 1
      do i1 = -1, 1
        iv = iv + 1
        do i0 = 1, 3
          tranvec(i0,iv) = i1*ucell(i0,1) + i2*ucell(i0,2) + i3*ucell(i0,3)
        end do
      end do
    end do
  end do

  return

end subroutine fill_tranvec 

!------------------------------------------------------------------------------
! fill_xtran: neighbor cells of a cell of interest (call it A) that are only
!             "forward" of that cell reside in the plane of the cell,
!             and three planes "forward" in the z direction will fill the
!             cell translation array in order to construct a neighbor list.
!
!             A = cell for which we want to get neighbor list
!             x = forward neighbor cell within 3 cells
!             o = same as x except this cell has same x index as A
!
!            i3          i3+1            i3+2         i3+3
!        ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!        ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!        ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!             Axxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!                      ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!                      ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!                      ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!
!
!
! A cell and its neighbors span the x direction over 7 cells (3 on each side
! and the cell iself).  The xtran array contains a 0, 1, or -1 for whether the
! neighbor cell has x index outside the unit cell and needs to be translated
! along the x uc vector positive one cell, negative, or not at all.  There are
! 10 cases of neighbor cells in x direction
!
! - xtran(*,1) (0000000) All have x indices within the unit cell
!
! - cases 2, 3, and 4 are special for the row of neighbors containing the cell
!   A itself (see arrow). This is the only case where neighbors with x index
!   are not included in the list search since those cells are "before" the cell
!   of interest, and this is a "look forward only" method.  So, for this row of
!   cells, only 4 xtran values are needed: cell A, and the three cells to
!   right.  The cases represent whether the neighbors extend out of the unit
!   cell by one, two, or three cells. Entry 1 is for cell A and must be 0
!   since it must be in the unit cell. (The last 4 entries are ignored for
!   this set):
!          (*,2) (0001000)
!          (*,3) (0011000)
!          (*,4) (0111000)
!
! - cases 5, 6, and 7 are same as 2,3,4 except that there are 7 cells in all
!   other rows:
!          (*,5) (0000001)
!          (*,6) (0000011)
!          (*,7) (0000111)
!
! - cases 8, 9, and 10 are for neighbors that extend to the left out of the UC.
!          (*,8) (-1000000)
!          (*,9) (-1-100000)
!          (*,10)(-1-1-10000)
!------------------------------------------------------------------------------
subroutine fill_xtran()

  implicit none

  integer i, j

  ! Hard Wired for nghb1=3 thus the hard dimensions (7,10)
  REQUIRE(nghb1 == 3)

  ! Most cells will not be translated, so set xtran to 0
  ! for all possibilities.  Then, fill in the 1 and -1
  ! entries for neighbor cells that are beyond uc edges.

  ! Case 1: zero out entire array
  do i = 1, 2*nghb1+1
    do j = 1, nghb1*3+1
      xtran(i,j) = 0
    end do
  end do

  ! Cases 2, 3, and 4
  do i = 0, nghb1-1
    do j = nghb1+1-i, nghb1+1
      xtran(j,i+2)=1
     end do
  end do

  ! Cases 5, 6, and 7
  do i = 1, nghb1
    do j = 1, i
      xtran(j,nghb1+1+i) = -1
    end do
  end do

  ! Cases 8, 9, and 10
  do i = 0, nghb1-1
    do j = 0,i
      xtran(2*nghb1+1-j,2*nghb1+2+i)=1
    end do
  end do

  return

end subroutine fill_xtran 

end module nblist
