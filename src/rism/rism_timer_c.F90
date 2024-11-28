! <compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2011-2012 by 
!Tyler Luchko and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
!ibid. 112:10391-10417 (2000).   
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 

#include "../include/dprec.fh"

!> Timer module.  Used by but not limited to RISM.
!! - Each timer has parent and an unlimited number of children.
!! - If this timer is the top most timer, its parent is NULL.
!! - The total time for a timer includes the time from its children.
!! - To avoid double counting,
!!   + Children start their parent, if it is not running.
!!   + Children stop their parent, if they started it.
!!   + The timer must be stopped/started in order to start/stop.
!!   + No other branch may be running. I.e., a parent may only have one child running.
!! - When a timer is destroyed and it is has a parent, the timer is duplicated
!!  and stored with its parent.  These copies are called 'displaced'. This preserves 
!!  it for timer summaries.  All duplicates are destroyed when the top level timer 
!!  is destroyed.
module rism_timer_c
  use rism_report_c
  implicit none
  type rism_timer_p
     type(rism_timer), pointer :: p=>NULL()
  end type rism_timer_p
  integer, private, parameter :: CLEN=32
  type rism_timer
     !name :: text name for this timer
     character(len=CLEN) :: name=''
     !running :: is this timer running?
     !startedParent:: Did this timer start its own parent?
     !displaced:: true if the instantiating object has destroyed this class.  
     !           This should be true to safely destroy this object.
     logical :: running, startedParent, displaced
     !sublevels : number of levels of child timers below this one.
     integer :: subLevels
     !total :: total accumulated time in seconds.
     _REAL_ :: total
     !timeStamp :: Time stamp for the last time this timer was
     !             started is seconds.  
     _REAL_ :: timeStamp
     !parent :: pointer to parent timer
     type(rism_timer_p) :: parent
     !child :: list of child timers
     !displacedChild :: clones of timers that have been displaced.  Kept around
     !                   for timer summaries
     type(rism_timer_p),pointer :: child(:)=>NULL(), displacedChild(:)=>NULL()
     !childID :: index of this timer in its parent's child() array
     integer :: childID
  end type rism_timer
public rism_timer_new, rism_timer_destroy, rism_timer_start, rism_timer_stop, &
     rism_timer_getTime, rism_timer_getTotalTime
private setChild, summary, destroy, destroyTree, displacedChild, replace
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Creates a new timer with a given name and a parent.
!!!IN:
!!!   this : timer object
!!!   name : text identifier for this timer
!!!   parent : (optional) parent timer.  If this the new timer will be the top
!!!            timer, do not use this argument
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism_timer_new(this,name, parent)
    implicit none
    type(rism_timer), intent(inout) :: this
    type(rism_timer), optional, target, intent(inout) :: parent
    character(len=*) :: name
    this%name = trim(name)
    this%running = .false.
    this%startedParent = .false.
    this%displaced = .false.
    this%total = 0
    this%timestamp = 0
    this%childID = 0
    this%subLevels = 0
    if(present(parent))then
       call rism_timer_setParent(this,parent)
    end if
    if(associated(this%child))then
       deallocate(this%child)
    end if
    allocate(this%child(0)) !no children until this is identified as a parent
    if(associated(this%displacedChild)) then
       deallocate(this%displacedChild)
    end if
    allocate(this%displacedChild(0))
  end subroutine rism_timer_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Destroys a timer. iff it has a parent timer a copy is maintained.  When the 
!!!root timer is destroyed all of the copies are destroyed.  If non-copy 
!!!children still exist, it is an error. This means that all timers exist 
!!!when a timer summary is requested.  However, the timer tree only grows during 
!!!a run.  Further more, if the root timer is destroyed before timed objects are 
!!!destroyed then null pointers exist.
!!!IN:
!!!   this : timer object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism_timer_destroy(this)
    implicit none
    type(rism_timer), intent(inout) :: this
    this%displaced = .true.
    if(.not.associated(this%parent%p))then
       call destroy(this)
    else
       !MAKE A CLONE AND STORE IT LOCALLY.  DESTROY THE ORIGINAL.  WE CAN'T COUNT 
       !ON THE ORIGINAL MEMORY BEING MAINTIAINED.
       call displacedChild(this%parent%p,this)
    end if
  end subroutine rism_timer_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!sets the parent timer
!!!IN:
!!!   this   : timer object
!!!   parent : parent timer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism_timer_setParent(this,parent)
    implicit none
    type(rism_timer), target, intent(inout) :: this, parent
    this%parent%p=>parent
    this%childID = setChild(parent,this)
  end subroutine rism_timer_setParent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Starts a timer object.  Object must not be running. No other timer
!!!branches may be running.
!!!IN:
!!!   this : timer object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine rism_timer_start(this)
    implicit none
    type(rism_timer), intent(inout) :: this
    integer :: count, count_rate, count_max
    character(len=10) :: date, time, zone
    integer :: values(8), ichild

    if(this%running)&
         call rism_report_error("timer already running. Cannot start: "&
         //trim(this%name))
    !If the parent is running, start it
    if(associated(this%parent%p))then
       if(.not. this%parent%p%running)then
          call rism_timer_start(this%parent%p)
          this%startedParent = .true.
       elseif(this%parent%p%running)then
          !check if another sibling is already running.  Since this
          !recursively climbs the tree this checks if another branch is
          !running
          do ichild=1,ubound(this%parent%p%child,1)
             if(this%parent%p%child(ichild)%p%running)&
                  call rism_report_error("another timer branch ("&
                  //trim(this%parent%p%child(ichild)%p%name)&
                  //") already running. Cannot start: "&
                  //trim(this%name)//". Common parent: "&
                  //trim(this%parent%p%name)//".")
          end do
       end if
    end if
    !get time stamp
    call cpu_time(this%timestamp)
    this%running = .true.
  end subroutine rism_timer_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Stops a timer object.  Object must be running. Parent must
!!!be running.
!!!IN:
!!!   this : timer object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine rism_timer_stop(this)
    implicit none
    type(rism_timer), intent(inout) :: this
    _REAL_ :: endtimestamp
    integer :: ichild
    if(.not. this%running)&
         call rism_report_error("timer not running. Cannot stop: "&
         //trim(this%name))
    this%running = .false.
    !no children should be running, even if they started this timer
    do ichild=1, ubound(this%child,1)
       if(this%child(ichild)%p%running)&
            call rism_report_error("Cannot stop timer '"//trim(this%name)&
            //"'. Child still running: "//trim(this%child(ichild)%p%name))
    end do
    !get time stamp
    call cpu_time(endtimestamp)
    !increment total time
    this%total= this%total + (endtimestamp - this%timestamp)
    !stop parent
    if(this%startedParent)then
       call rism_timer_stop(this%parent%p)
       this%startedParent = .false.
    end if
  end subroutine rism_timer_stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the total time is milliseconds this timer object has been running.
!!!IN:
!!!   this : timer object
!!!OUT:
!!!   total accumulated run time so far
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism_timer_getTotalTime(this) result(time)
    implicit none
    type(rism_timer), intent(inout) :: this
    _REAL_ :: time
    if(this%running)then
       !get time stamp
       call cpu_time(time)
       !increment total time
       this%total= this%total + (time - this%timestamp)
       !update the time stamp
       this%timestamp = time
    end if
    time = this%total
  end function rism_timer_getTotalTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns an array for the total time is milliseconds for all of the children 
!!!of this timer object and the time not spent in childer timers (other).
!!!IN:
!!!   this : timer object
!!!OUT:
!!!   One index for each of the children plus the last one for the remaining
!!!   time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism_timer_getTime(this) result(time)
    implicit none
    type(rism_timer), intent(inout) :: this
    integer :: nchild
    _REAL_ :: time(ubound(this%child,1)+1)
    _REAL_ :: total
    integer :: ichild
    nchild=ubound(this%child,1)
    total = rism_timer_getTotalTime(this)
    do ichild=1,nchild+1
       time(ichild) = rism_timer_getTotalTime(this%child(ichild)%p)
    end do
    time(nchild+1) = total - sum(time(1:nchild))
  end function rism_timer_getTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Writes out a recurrsive summary of the time taken for this timer and its 
!!!children 
!!!IN:
!!!   this : timer object
!!!   o_comment : (optional) Line comment character to place in front of 
!!!               the line
!!!   o_unit    : (optional) output unit number
!!!   o_mpicomm : (optional) MPI communicator
!!!Sideeffect:
!!!   writes to standard out unless a unit number is provided
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism_timer_summary(this,o_comment,o_unit,o_mpicomm)
    implicit none
#ifdef MPI
  include 'mpif.h'
#endif /*MPI*/  
    type(rism_timer), intent(inout) :: this
    character(len=*), optional, intent(in) :: o_comment
    integer, optional, intent(in) :: o_unit, o_mpicomm
    character(len=4) :: comment
    integer :: unit
    integer :: mpirank, mpicomm, err
    comment=''
    if(present(o_comment)) comment = o_comment
    unit = rism_report_getMUnit()
    if(present(o_unit)) unit = o_unit

    mpirank = 0
#ifdef MPI
    mpicomm =  MPI_COMM_NULL
    if(present(o_mpicomm)) mpicomm = o_mpicomm
    if(mpicomm /= MPI_COMM_NULL)&
         call mpi_comm_rank(mpicomm,mpirank,err)
    if(err /=0) call rism_report_error&
         ("(a,i8)","RISM_TIMER: could not get MPI rank for communicator ",mpicomm)
#else
    mpicomm = 0
#endif /*MPI*/  
    call summary(this,0,trim(comment), 0.01d0, unit, mpirank, mpicomm)
    if(mpirank == 0)then
       write(unit,'(a)')
    end if
  end subroutine rism_timer_summary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Does the actual destruction of the timer.
!!!IN:
!!!   this : timer object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine destroy(this)
    implicit none
    type(rism_timer), intent(inout) :: this
    integer :: ichild
    if(.not.this%displaced) &
         call rism_report_error("TIMER: attempting to destroy active timer: "//trim(this%name))
    this%name = ''
    this%running = .false.
    this%total = 0
    this%timestamp = 0
    nullify(this%parent%p)
    this%subLevels = 0
    
    if (associated(this%child)) then
       do ichild=1, ubound(this%child,1)
          ! Remove the parent of each child.
          nullify(this%child(ichild)%p%parent%p)
          ! Falsify 'startedParent' so we don't have the case where
          ! the child tries to stop a non-existant parent.
          this%child(ichild)%p%startedParent=.false.
       end do
       deallocate(this%child)
       nullify(this%child)
    end if
    ! Displaced children have already had a destroy request. It is now
    ! safe to get rid of them completely.
    if (associated(this%displacedChild)) then
       do ichild = 1, ubound(this%displacedChild, 1)
          call destroy(this%displacedChild(ichild)%p)
          deallocate(this%displacedChild(ichild)%p)
       end do
       deallocate(this%displacedChild)
       nullify(this%displacedChild)
    end if
  end subroutine destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Destroyes the timer tree from this timer down.  Recursively calls this method for
!!!all children
!!!IN:
!!!   this : timer object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine destroyTree(this)
    implicit none
    type(rism_timer), intent(inout) :: this
    integer :: ichild
    do ichild=1, ubound(this%child,1)
       call destroyTree(this%child(ichild)%p)
    end do
    call destroy(this)
  end subroutine destroyTree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!When a child timer becomes displaced (it is destroyed by its host object) it 
!!!needs to be preserved.  This subroutine stores a clone of the child locally
!!!for future use (timer summaries) and destroys the original.
!!!IN:
!!!   this : timer object
!!!   child : child timer being displaced
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine displacedChild(this, child)
    implicit none
    type(rism_timer), intent(inout) :: this
    type(rism_timer),target, intent(inout) :: child
    !displaced :: the copy of child that will be persisent in memory
    type(rism_timer),pointer :: displaced
    !tempChild :: duplicate of this%child
    type(rism_timer_p):: tempChild(ubound(this%Child,1))
    !tempDisplaced :: duplicate of this%displacedChild
    type(rism_timer_p):: tempDisplace(ubound(this%displacedChild,1))
    integer :: ichild
    
    !grow displacedChild
    tempDisplace = this%displacedChild
    if(associated(this%displacedChild))&
         deallocate(this%displacedChild)
    allocate(this%displacedChild(ubound(tempDisplace,1)+1))
    this%displacedChild(1:ubound(tempDisplace,1)) = tempDisplace

    !shrink child array
    tempChild=this%child
    if(associated(this%Child))&
         deallocate(this%Child)
    allocate(this%child(ubound(tempChild,1)-1))
    if(child%childID > 1)&
         this%Child(1:child%childID-1) = tempChild(1:child%childID-1)
    if(child%childID<=ubound(this%child,1))&
         this%Child(child%childID:) = tempChild(child%childID+1:)
    !Update childIDs
    do ichild = child%childID, ubound(this%child,1)
       this%child(ichild)%p%childID=ichild
    end do

    !point the final index of displacedChild at child and update ID
    this%displacedChild(ubound(this%displacedChild,1))%p => child
    child%childID = ubound(this%displacedChild,1)
    !make a copy of child and destroy it since the original may be reused
    allocate(displaced)
    call replace(child, displaced)
  end subroutine displacedChild

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!replaces THIS with CLONE, which is otherwise a copy.  THIS is destroyed.
!!!IN:
!!!   this : timer object
!!!   clone : copy of this
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine replace(this,clone)
    implicit none
    type(rism_timer), intent(inout) :: this
    type(rism_timer), target, intent(inout) :: clone
    integer :: ichild

    !copy vital statistics
    clone%name = trim(this%name)    
    clone%running = this%running
    clone%startedParent = this%startedParent
    clone%displaced = this%displaced
    clone%total = this%total
    clone%timestamp = this%timestamp
    clone%childID = this%childID
    !update parent and change pointer in the appropriate list
    !This is the actual replacement
    if(associated(this%parent%p))then
       clone%parent%p => this%parent%p
       if(.not.this%displaced)then
          clone%parent%p%child(clone%childID)%p => clone
       else
          clone%parent%p%displacedChild(clone%childID)%p => clone
       end if
    end if
    !copy children and update their parent
    allocate(clone%child(ubound(this%child,1)))
    clone%child = this%child
    do ichild = 1, ubound(clone%child,1)
       clone%child(ichild)%p%parent%p => clone
    end do
    !copy displaced and update their parent
    allocate(clone%displacedChild(ubound(this%displacedChild,1)))
    clone%displacedChild = this%displacedChild
    do ichild = 1, ubound(clone%displacedChild,1)
       clone%displacedChild(ichild)%p%parent%p => clone
    end do
    !deallocating the displacedChild array prevents destroy() from
    !deleting all of the copies.
    deallocate(this%displacedChild)
    call destroy(this)

  end subroutine replace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Writes out a recurrsive summary of the time taken for this timer and its 
!!!children 
!!!IN:
!!!   this : timer object
!!!   level : indentation level. level=0 indicates this is the top
!!!           level and a header is printed
!!!   comment : Line comment character to place in front of the line
!!!   cutoff : percentage of parents total time required to be reported
!!!   unit   : unit to output to
!!!   mpirank : rank of MPI process
!!!   mpicomm : MPI communicator
!!!Sideeffect:
!!!   writes to standard out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  recursive subroutine summary(this,level,comment,cutoff, unit, mpirank, mpicomm,o_offset)
    implicit none
#ifdef MPI
  include 'mpif.h'
#endif /*MPI*/  
    type(rism_timer), intent(inout) :: this
    integer, intent(in) :: level
    character(len=*), intent(in) :: comment
    _REAL_, intent(in) :: cutoff
    integer, intent(in) :: unit, mpirank, mpicomm
    integer, optional, intent(in) :: o_offset
    integer :: offset
    integer :: ichild, i
#ifdef MPI
    integer :: mpisize, err
    !In MPI we report extra information.
    !(1) : average across all processes
    !(2) : minimum across all processes
    !(3) : maximum across all processes
    _REAL_ :: otherTime(3), totalTime(3), parentTime(3),&
         mpi_buffer(3)
#else
    _REAL_ :: otherTime(1), totalTime(1), parentTime(1)
#endif /*MPI*/
    character(len=CLEN) :: name
    character(len=40) :: namefmt, timefmt, colheadfmt, titlefmt
    character(len=2) :: tab
    character(len=128) :: teststr
    tab(:)=' '
    offset = this%subLevels*len(tab)
    if(present(o_offset)) offset = o_offset
    write(namefmt,'(a,i3,a)') '(a',CLEN,')'
    if(offset >=0)then
       write(timefmt,'(a,i3,a)') "(",offset,"x,f10.3,3x,'(',f6.2,'%)')"
    else
       timefmt = "(f10.3,3x,'(',f6.2,'%)')"
    end if
    if(level == 0 .and. mpirank == 0)then
       write(unit,'(a)')
       !write out title in the space over the timer names and max indentation
       write(titlefmt,'(a,i3,a)') '(a',CLEN+offset,')'
       write(unit,titlefmt,advance='no') trim(comment)//"CPU Timer summary for "//trim(this%name)//&
            repeat(" ",CLEN+offset)
#ifdef MPI
!!$       write(unit,'(a)',advance='no') repeat(" ",offset)
       !create a sample string to get the length of the numeric output
       write(teststr,timefmt) 1d0, 1d0
       !use this to create a format.  This will truncate extra white
       !space so we won't worry about rounding below
       write(colheadfmt,'(a,i3,a)') '(a',len_trim(teststr),')'
       !write the headers
       write(unit,colheadfmt,advance='no') repeat(" ",len_trim(teststr)/2-len("Average")/2)//&
            "Average"//repeat(" ",len_trim(teststr)/2-len("Average")/2)
       write(unit,colheadfmt,advance='no') repeat(" ",len_trim(teststr)/2-len("Minimum")/2)//&
            "Minimum"//repeat(" ",len_trim(teststr)/2-len("Minimum")/2)
       write(unit,colheadfmt,advance='no') repeat(" ",len_trim(teststr)/2-len("Maximum")/2)//&
            "Maximum"//repeat(" ",len_trim(teststr)/2-len("Maximum")/2)
#endif /*MPI*/
       write(unit,'(a)')
    endif

    !Each process collect their repsective times
    totalTime = rism_timer_getTotalTime(this)
    otherTime = totalTime
    if(associated(this%parent%p))then
       parentTime = rism_timer_getTotalTime(this%parent%p)
    else
       parentTime = totalTime
    end if
    do ichild = 1, ubound(this%child,1)
       otherTime = otherTime - rism_timer_getTotalTime(this%child(ichild)%p)
       call summary(this%child(ichild)%p, level+1,comment, cutoff, unit,mpirank,mpicomm,offset)
    end do
    do ichild = 1, ubound(this%displacedChild,1)
       otherTime = otherTime - rism_timer_getTotalTime(this%displacedChild(ichild)%p)
       call summary(this%displacedChild(ichild)%p, level+1,comment, cutoff, unit,mpirank,mpicomm,offset)
    end do

#ifdef MPI
     call mpi_comm_size(mpicomm,mpisize,err)
     if(err /=0) call rism_report_error&
          ("(a,i8)","RISM_TIMER: could not get MPI size for communicator ",mpicomm)
    !get average (first index), min (second index) and max (third index) for each catagory
     call MPI_REDUCE((/otherTime(1),totalTime(1),parentTime(1)/), mpi_buffer, 3,MPI_DOUBLE_PRECISION,&
          MPI_SUM,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("(a,i8)","RISM_TIMER: MPI_REDUCE failed ",err)
     mpi_buffer = mpi_buffer/mpisize
     otherTime(1) = mpi_buffer(1)
     totalTime(1) = mpi_buffer(2)
     parentTime(1) = mpi_buffer(3)
     
     call MPI_REDUCE((/otherTime(2),totalTime(2),parentTime(2)/), mpi_buffer, 3,MPI_DOUBLE_PRECISION,&
          MPI_MIN,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("(a,i8)","RISM_TIMER: MPI_REDUCE failed ",err)
     otherTime(2) = mpi_buffer(1)
     totalTime(2) = mpi_buffer(2)
     parentTime(2) = mpi_buffer(3)

     call MPI_REDUCE((/otherTime(3),totalTime(3),parentTime(3)/), mpi_buffer, 3,MPI_DOUBLE_PRECISION,&
          MPI_MAX,0,mpicomm,err)
     if(err /=0) call rism_report_error&
          ("(a,i8)","RISM_TIMER: MPI_REDUCE failed ",err)
     otherTime(3) = mpi_buffer(1)
     totalTime(3) = mpi_buffer(2)
     parentTime(3) = mpi_buffer(3)
#endif /*MPI*/    
    if(mpirank == 0)then
       if(ubound(this%child,1) > 0 .and. otherTime(1)/totalTime(1)*100d0 >= cutoff )then
          write(unit,'(a)',advance='no') comment//repeat(tab,level+1)
          call flush(unit)
          write(unit,namefmt,advance='no') 'Other'//repeat(' ',CLEN-len_trim('Other'))
          do i=1, size(totalTime)
             write(unit,timefmt,advance='no')     otherTime(i),otherTime(i)/totalTime(i)*100d0       
          end do
          write(unit,'(a)')
          call flush(unit)
       end if
       if(totalTime(1)/parentTime(1)*100d0 >= cutoff)then
          write(unit,'(a)',advance='no') comment//repeat(tab,level)
          call flush(unit)
          write(unit,namefmt,advance='no') this%name//repeat(' ',CLEN-len_trim(this%name))
          do i=1, size(totalTime)
             write(unit,timefmt,advance='no') totalTime(i), totalTime(i)/parentTime(i)*100d0
          end do
          write(unit,'(a)')
          call flush(unit)
       end if
    end if
  end subroutine summary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!adds a child timer
!!!IN:
!!!   this   : timer object
!!!   child  : child timer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function setChild(this,child) result(childID)
    implicit none
    type(rism_timer), target, intent(inout) :: this
    type(rism_timer), target, intent(in) :: child
    type(rism_timer_p),pointer :: temp(:)=>NULL();
    integer :: childID
    type(rism_timer),pointer :: p=>NULL(),c=>NULL()
    if(.not.associated(this%child))then
       childID = 1
       allocate(this%child(1))
    else
       childID = ubound(this%child,1)+1
       allocate(temp(childID - 1))
       temp = this%child
       deallocate(this%child)
       allocate(this%child(childID))
       this%child(:childID-1) = temp
       deallocate(temp)
    end if
    !update the number of sublevels in parent timers.
    !The number of sublevels for a parent should be one greater than
    !the max sublevels of all its children. So, if the new number of
    !sublevels for the child is >= that of the parent, increment the
    !parents sublevels.  Otherwise we are done.
    this%child(childID)%p => child
    p=>this
    c=>child
    !walk up the direct parents in the tree and stop at the top
    do while(associated(p))
       !Only update if this really adds a sublevel (no other child has more sublevels.)
       if(p%subLevels <= c%subLevels)then
          p%subLevels = c%subLevels+1
          c => p
          p => p%parent%p
       else
          exit
       end if
    end do
  end function setChild
end module rism_timer_c
