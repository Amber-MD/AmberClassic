#include "copyright.h"
#include "../include/dprec.fh"
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ The GBNSR6 main program
program gbnsr6main
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !  Saeed Izadi,
   !  Boris Aguilar,
   !  Abhishek Mukhopadhyay,
   !  Alexey Onufriev,
   !  Lab for Theoretical and Computational Molecular Biophysics,
   !  Virginia Tech, Blacksburg
   !
   !  Call gbnsr6 to perform calculations.
   !
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

#  include "files.h"
#  include "extra.h"

   numgroup = 1
   call gbnsr6file()
   
   master   = .true.

   call gbnsr6()

   if ( master ) call mexit(6,0); call mexit(0,0)

end program gbnsr6main

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ The PBSA driver
subroutine gbnsr6()

!  use parms ! Saeed: replaced with parms.h 
   use memory_module !Saeed: replaced with memory.h
   use pbtimer_module
   use decomp, only : allocate_int_decomp, allocate_real_decomp, &
                      deallocate_int_decomp, deallocate_real_decomp
   implicit none

#  include "files.h"
#  include "bc_memory.h"
#  include "box.h"
#  include "md.h"
#  include "parms.h"
#  include "timer.h"
#  include "extra.h"

   logical                          erstop
   integer                          nr
   integer                          i_stack(1)
   integer                          atm1, atm2, ntrns, rotopt, mrot, mrotx, msph
   integer                          nrotx, nsph
   integer, allocatable          :: ix(:), ipairs(:)
   character(len=8 )                initial_date, setup_end_date, final_date
   character(len=10)                initial_time, setup_end_time, final_time
   character(len=4), allocatable :: ih(:)
   _REAL_                           ene(51)
   _REAL_                           carrms
   _REAL_                           r_stack(1)
   _REAL_                           delta
   _REAL_,  allocatable          :: x(:)

   ! LOCAL
   logical                          belly
   integer                          ier, ncalls
   !integer                          iatm

   ier = 0
   
   ! Initialize the cpu timer. Needed for machines where returned cpu times
   ! are relative.

   call date_and_time( initial_date, initial_time )
   call pbtimer_init()
    
   erstop = .false.

   ! Only the master node (only node when single-process)
   ! performs the initial setup and reading/writing

   call pbtimer_start(PBTIME_TOTAL)
   if ( master ) then

      !        --- first, initial reads to determine memry sizes:

      call pbtimer_start(PBTIME_READ)
      call mdread1()
!     call openparm(parm//CHAR(0))

      call myopen(8,parm,'O','F','R')
      call rdparm1(8)
      call gb_read(natom)

      !        --- now, we can allocate memory:

      call locmem
   end if

   !        --- dynamic memory allocation:

   REQUIRE( lastr>0 .and. lasti>0 .and. lastpr>0 .and. lasth>0 )
   allocate( x      (lastr), stat = ier ); REQUIRE(ier==0)
   allocate( ix     (lasti), stat = ier ); REQUIRE(ier==0)
   allocate( ipairs(lastpr), stat = ier ); REQUIRE(ier==0)
   allocate( ih     (lasth), stat = ier ); REQUIRE(ier==0)

   if( idecomp > 0 ) then
      call allocate_int_decomp(natom, nres)
   else
      call allocate_int_decomp(1, 1)
   endif

   lastrst = 1
   lastist = 1
   r_stack(1) = 0.0d0
   i_stack(1) = 0

   if ( master ) then
      write(6,'(/,a,5x,a)') '|','Memory Use     Allocated'
      write(6,'(a,5x,a,i14)') '|', 'Real      ', lastr
      write(6,'(a,5x,a,i14)') '|', 'Hollerith ', lasth
      write(6,'(a,5x,a,i14)') '|', 'Integer   ', lasti
      write(6,'(a,5x,a,i14)') '|', 'Max Pairs ', lastpr
      write(6,'(a,5x,a,i14)') '|', 'Max Rstack', lastrst
      write(6,'(a,5x,a,i14)') '|', 'Max Istack', lastist
      write(6,'(a,5x,a,i14,a)') '|', '  Total   ', &
           (8*(lastr+lastrst) + 4*(lasth+lasti+lastpr+lastist))/1024, ' kbytes'

      !        --- second reads to finalize

!     call rdparm2(x,ix,ih)
      call rdparm2(x,ix,ih,ipairs,8,i_stack) !saeed: reading crd/top
      call mdread2(x,ix,ih)
   end if

!  --- alloc memory for decomp module that needs info from mdread2
   if( idecomp == 1 .or. idecomp == 2 ) then
      call allocate_real_decomp(nres)
   else if( idecomp == 3 .or. idecomp == 4 ) then
      call allocate_real_decomp(npdec*npdec)
   end if

   ! EVALUATE SOME CONSTANTS FROM MDREAD SETTINGS

   nr = nrp
   belly = ibelly > 0

   if ( master ) then
      ! READ COORDINATES AND VELOCITIES

      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,t)
      if ( belly ) call bellyf()
      call pbtimer_stop(PBTIME_READ)
      
      ! Compatibility print-out please ignore

      write(6,'(" Number of triangulated 3-point waters found: ",i8)') 0

      ! OPEN THE DATA DUMPING FILES AND POSITION IT DEPENDING ON THE TYPE OF RUN

      ! call open_dump_files

      call amflsh(6)

      ! end of master process setup
   end if  ! (master)
   
   call date_and_time( setup_end_date, setup_end_time )

   ! Use the debugf namelist to activate
   call debug_frc(x,ix,ih,ipairs,x(lcrd),x(lforce),cn1,cn2)

   ! Now do the dynamics or minimization.

   if ( master ) write(6,'(/80(1H-)/''   3.  RESULTS'',/80(1H-)/)')

   ! Input flag imin determines the type of calculation: MD, minimization, ...

   select case ( imin )
   case ( 0 )
      ! Dynamics:

      call pbtimer_start(PBTIME_RUNMD)
      call runmd(x,ix,ih,ipairs, &
              x(lcrd),x(lwinv),x(lmass),x(lforce), &
              x(lvel),x(lvel2),x(l45),x(lcrdr), &
              x(l50),x(l95),ix(i70),x(l75),erstop,r_stack,i_stack)
      call pbtimer_stop(PBTIME_RUNMD)

      if (master) call amflsh(6)

   case ( 1 )
      ! Minimization:

      ! Input flag ntmin determines the method of minimization
      select case ( ntmin )
      case ( 0, 1, 2 )
         call runmin(x,ix,ih,ipairs,x(lcrd),x(lforce),x(lvel), &
                 ix(iibh),ix(ijbh),x(l50),x(lwinv),ix(lasti), &
                 x(l95),ene,r_stack,i_stack, carrms)   !saeed: up to here
!everything is needed 
      case default
         ! invalid ntmin
         ! ntmin input validation occurs in mdread.f
         if (master) write(6,'(/2x,a,i3,a)') 'Error: Invalid NTMIN (',ntmin,').'
         call mexit(6,0)
      end select

      if (master) call amflsh(6)

   case ( 5 )
      ! trajene option not supported
      ! imin input validation is in mdread.f
      if (master) write (6,'(a)') 'Error: Post-processing of trajectory not supported'
      call mexit(6,0)
   case ( 6 )
      ! force check
      atm1 = 1
      atm2 = 2
      ntrns = 1
      rotopt = 1 ! 1, spherical; otherwise, random
      nrotx = 3
      nsph =  40 
      mrot = 100
      mrotx = 100 
      msph = 1000
      delta = 1.0d-4
     ! call chkfrc(x,ix,ih,ipairs,x(lcrd),x(lforce),ene, &
     !         r_stack,i_stack,atm1,atm2,ntrns,rotopt,nrotx,nsph,mrot,mrotx,msph,delta)
   case default
      ! invalid imin
      ! imin input validation is in mdread.f
      if (master) write(6,'(/2x,a,i3,a)') 'Error: Invalid IMIN (',imin,').'
      call mexit(6,0)
   end select

!  do iatm = 1, natom
!     atm1 = lcrd+3*(iatm-1)
!     write(58,'(a,3(f15.10),2x)') ih(m04+iatm-1), x(atm1:atm1+2)
!  enddo

   !     -- calc time spent running vs setup

   call pbtimer_stop(PBTIME_TOTAL)
   call date_and_time( final_date, final_time )

   if ( master ) then

   !call pbtimer_summary()

      !     --- write out final times

      write(6,'(12(a))') '|           Job began  at ', initial_time(1:2), &
           ':', initial_time(3:4), ':', initial_time(5:10), '  on ',&
           initial_date(5:6), '/', initial_date(7:8), '/', initial_date(1:4)
      write(6,'(12(a))') '|           Setup done at ', setup_end_time(1:2),  &
           ':', setup_end_time(3:4), ':', setup_end_time(5:10), '  on ', &
           setup_end_date(5:6), '/',setup_end_date(7:8),'/',setup_end_date(1:4)
      write(6,'(12(a))') '|           Run   done at ', final_time(1:2),  &
           ':', final_time(3:4), ':', final_time(5:10), '  on ', &
           final_date(5:6), '/', final_date(7:8), '/', final_date(1:4)
      call nwallclock( ncalls )
      write(6, '(''|'',5x,''wallclock() was called'',I8,'' times'')') ncalls
   end if

   !     --- dynamic memory deallocation:

   call pb_free()
   deallocate(     ih, stat = ier ); REQUIRE(ier==0)
   deallocate( ipairs, stat = ier ); REQUIRE(ier==0)
   deallocate(     ix, stat = ier ); REQUIRE(ier==0)
   deallocate(      x, stat = ier ); REQUIRE(ier==0)


end subroutine gbnsr6

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ The gbnsr6 file handler
subroutine gbnsr6file

   implicit none

#  include "files.h"

   character(len=80) argbuf
   integer iarg ! index of the current argument
   integer pb_iargc ! wrapper to intrinsic that returns the index of the last argument
                    ! from either the command line or a string
   integer last_arg_index ! index of the last argument

   !     --- default file names ---

   mdin   = 'mdin'
   mdout  = 'mdout'
   inpcrd = 'inpcrd'
   parm   = 'prmtop'
   restrt = 'restrt'
   refc   = 'refc'
   mdvel  = 'mdvel'
   mden   = 'mden'
   mdcrd  = 'mdcrd'
   mdinfo = 'mdinfo'
   vecs   = 'vecs'
   freqe  = 'dummy'
   rstdip = 'rstdip'
   inpdip = 'inpdip'
   mddip  = 'mddip'
   radii  = 'radii'
   cpin   = 'cpin'
   cpout  = 'cpout'
   cprestrt = 'cprestrt'
   newtop = ' ' !B.Aguilar for mdar6
   if (numgroup == 1) groups = ' '

   !     --- default status of output: New

   owrite = 'N'

   !     --- get command line arguments ---

   iarg = 0
   last_arg_index = pb_iargc()
   do while (iarg < last_arg_index)
      iarg = iarg + 1

      call getpb_arg(iarg,argbuf)
      select case (argbuf)
      case (' ')
         continue
      case ('-O')
         owrite = 'R' !      status of output: Replace
      case ('-i')
         iarg = iarg + 1
         call getpb_arg(iarg,mdin)
      case ('-o')
         iarg = iarg + 1
         call getpb_arg(iarg,mdout)
      case ('-p')
         iarg = iarg + 1
         call getpb_arg(iarg,parm)
      case ('-c')
         iarg = iarg + 1
         call getpb_arg(iarg,inpcrd)
      case ('-md')
         iarg = iarg + 1
         call getpb_arg(iarg,newtop)
      case ('-vec')
         iarg = iarg + 1
         call getpb_arg(iarg,vecs)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-radii')
         iarg = iarg + 1
         call getpb_arg(iarg,radii)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-f')
         iarg = iarg + 1
         call getpb_arg(iarg,freqe)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-r')
         iarg = iarg + 1
         call getpb_arg(iarg,restrt)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-ref', '-z')
         iarg = iarg + 1
         call getpb_arg(iarg,refc)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-e')
         iarg = iarg + 1
         call getpb_arg(iarg,mden)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-v')
         iarg = iarg + 1
         call getpb_arg(iarg,mdvel)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-x', '-t')
         iarg = iarg + 1
         call getpb_arg(iarg,mdcrd)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-inf')
         iarg = iarg + 1
         call getpb_arg(iarg,mdinfo)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-idip')
         iarg = iarg + 1
         call getpb_arg(iarg,inpdip)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-rdip')
         iarg = iarg + 1
         call getpb_arg(iarg,rstdip)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-mdip')
         iarg = iarg + 1
         call getpb_arg(iarg,mddip)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-cpin')
         iarg = iarg + 1
         call getpb_arg(iarg,cpin)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-cpout')
         iarg = iarg + 1
         call getpb_arg(iarg,cpout)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case ('-cprestrt')
         iarg = iarg + 1
         call getpb_arg(iarg,cprestrt)
         write(6,'(a)') "Unsupported command arguments, exit.";call mexit(6,0)
      case default
         write(6,'(/,5x,a,a)') 'ERROR: Unknown argument: ',argbuf
         write(6,9000)
         call mexit(6, 1)
      end select 
   end do  !  while (iarg < last_arg_index)
 

   9000 format(/,5x, &
         'usage: gbnsr6  [-O] -i mdin -o mdout -p prmtop -c inpcrd ', &
         '-r restrt',/19x,'[-ref refc -x mdcrd -v mdvel -e mden ', &
         '-idip inpdip -rdip rstdip -mdip mddip ', &
         '-inf mdinfo -radii radii]' &
         , /, 'Consult the manual for additional options.')

end subroutine gbnsr6file

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ RESERVED FOR FUTURE MPI EXTENSION
integer function pb_iargc()
   implicit none
   integer iargc
   pb_iargc = iargc()

end function pb_iargc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ RESERVED FOR FUTURE MPI EXTENSION
subroutine getpb_arg(iarg, arg)
   implicit none
   integer iarg
   character(len=*) arg
   ! Intrinsic getarg requires a 4 byte integer argument; 
   ! this guards the argument for builds with default 8 byte integers.
   call getarg(int(iarg,4), arg)
 
end subroutine getpb_arg

!  Rewritten by: Meng-Juei Hsieh
subroutine mexit(filenum, exitstatus)
   implicit none
   integer filenum
   integer exitstatus

   if (filenum > 0) then  ! close this unit if greater than zero
      close(unit=filenum)
   endif
   call exit(exitstatus)
end subroutine mexit
