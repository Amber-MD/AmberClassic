#include "../include/dprec.fh"
#include "../include/assert.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read initial coordinates (and velocities) from restart file
#ifdef MPI
subroutine getcor(nr,x,v,f,ntx,box,irest,tt,temp0,writeflag,solvph,solve,stagid)
#else
subroutine getcor(nr,x,v,f,ntx,box,irest,tt,writeflag)
#endif
   
   !     --- reads initial coords, vel, and box for MD.
   
   !     Rev A mods: amopen calls, EOF protection, XC() no longer
   !                 read from restart file.  improved error message.
   !                 output new box lengths.
   
   !     EWALD: read a modified form of the coordinate/restart file
   
   !      IMPLICIT _REAL_ (A-H,O-Z)

   ! writeflag was introduced to avoid false failures in test.sander
   ! after atommask() dependency on getcor() was fixed by calling getcor()
   ! twice: once before atommask() calls and second time in its usual place.
   use constants, only : NO_INPUT_VALUE_FLOAT
   use AmberNetcdf_mod, only: NC_checkRestart
   use binrestart, only  : read_nc_restart, read_nc_remd_dimension, &
                           read_nc_remd_types
   use file_io_dat
   use sander_lib, only  : check_inpcrd_overflow, get_num_tokens
#ifdef MPI
   use remd, only : rem
   use sgld, only : isgld
#endif

   implicit none
#  include "ew_cntrl.h"
   integer lun,nr,ntx,irest,stagid
   _REAL_ x(*),v(*),f(*),box(3),tt
#ifdef MPI
   _REAL_ temp0,solvph,solve
   integer, dimension(:), allocatable :: lremd_types
#endif
   _REAL_, dimension(:), allocatable :: lremd_values ! for REMD restart
   
   logical form, writeflag,isRstValid,netcdf_input
   _REAL_ fvar(7),a,b,c,alpha,beta,gamma,input_time
   integer i,ivar,ifld(7)
   character(len=4) ihol(1)
   integer natom,nr3,ier,ibead
   character(len=MAX_LINE_BUF_LEN) line_test
   integer :: line_num = 0
   integer :: nwords,remd_values_dim,alloc_stat
  
   isRstValid = .false.
   netcdf_input = NC_checkRestart(inpcrd)
  
#ifdef MPI
   ! Setting the size of the local_remd_values array
   if (rem.ne.-1) then
     remd_values_dim = 1
   ! If doing Multi-D REMD and NetCDF restart
   else if (netcdf_input) then
     ! Getting remd_values_dim from restart file
     call read_nc_remd_dimension(inpcrd,title1,remd_values_dim)
     ! Is the remd_dimension info present in the restart file?
     if (remd_values_dim.ne.-1) then
       ! Allocating local_remd_types
       allocate(lremd_types(remd_values_dim), &
                stat = alloc_stat)
       REQUIRE(alloc_stat==0)
       ! Getting lremd_types from restart file
       call read_nc_remd_types(inpcrd,title1,lremd_types,remd_values_dim)
       ! Is the remd_types info present in the restart file?
       if (lremd_types(1).ne.-1) then
         isRstValid = .true.
       else
         remd_values_dim = 1
       end if
    else
      remd_values_dim = 1
    end if
   end if
#else
   remd_values_dim = 1
#endif
  
   ! Allocating local_remd_values
   allocate(lremd_values(remd_values_dim), &
            stat = alloc_stat)
   REQUIRE(alloc_stat==0)

   lremd_values(:) = 0.0d0
   nr3 = 3*nr
   lun = INPCRD_UNIT !FIXME: use lun=-1 for new amopen()

   if (writeflag) write(6,9108)
   
   !     ----- OPEN THE COORDINATE FILE -----

   ! Netcdf Restart
   if (netcdf_input) then
      call read_nc_restart(inpcrd,title1,ntx,nr,x,v,lremd_values,remd_values_dim,input_time)
      ! ntx=1 No Velocity Read, set to 0.0 
      if (ntx == 1) v(1:nr3) = 0.d0
#ifdef MPI
      ! in REM, overwrite temp0, solvph and/or solve if it's present in inpcrd/restrt
      ! DAN ROE: May want to just use temps in input file, so only 
      !          do this if a restart is requested.
      if(rem == 1 .and. irest==1) then
        if (isgld > 0) then
          if (lremd_values(1) /= NO_INPUT_VALUE_FLOAT) then
            write(6, '(a)') '| Correcting stagid using value from inpcrd file'
            stagid = lremd_values(1)
          else
            stagid = 0
          end if
        else if (lremd_values(1) /= NO_INPUT_VALUE_FLOAT) then
          write(6, '(a)') '| Overwriting temp0 from mdin with temp0 from &
                          &netcdf inpcrd file'
          temp0 = lremd_values(1)
        end if
      else if(rem == 4 .and. lremd_values(1) /= NO_INPUT_VALUE_FLOAT .and. irest==1) then
         write(6, '(a)') '| Overwriting solvph from mdin with solvph from &
                          &netcdf inpcrd file'
         solvph = lremd_values(1)
      else if(rem == 5 .and. lremd_values(1) /= NO_INPUT_VALUE_FLOAT .and. irest==1) then
         write(6, '(a)') '| Overwriting solve from mdin with solve from &
                          &netcdf inpcrd file'
         solve = lremd_values(1)
      ! Multi-D REMD
      else if(rem == -1 .and. irest==1 .and. isRstValid) then
          do i = 1, remd_values_dim
          if(isgld > 0) then
             if (lremd_types(i) == 1 .and. lremd_values(i) /= NO_INPUT_VALUE_FLOAT) then
               write(6, '(a)') '| Correcting stagid using value from inpcrd file'
               stagid = lremd_values(1)
             else
               stagid = 0
             end if

          else if (lremd_types(i) == 1 .and. lremd_values(i) /= NO_INPUT_VALUE_FLOAT) then
             write(6, '(a)') '| Overwriting temp0 from mdin with temp0 from &
                              &netcdf inpcrd file'
             temp0 = lremd_values(i)
           else if (lremd_types(i) == 4 .and. lremd_values(i) /= NO_INPUT_VALUE_FLOAT) then
             write(6, '(a)') '| Overwriting solvph from mdin with solvph from &
                              &netcdf inpcrd file'
             solvph = lremd_values(i)
           else if (lremd_types(i) == 5 .and. lremd_values(i) /= NO_INPUT_VALUE_FLOAT) then
             write(6, '(a)') '| Overwriting solve from mdin with solve from &
                              &netcdf inpcrd file'
             solve = lremd_values(i)
           end if
         end do
      end if
#endif
      if (writeflag) then
         write(6,9008) title1
         write(6,9009) input_time
      end if
      ! If restarting, set the time from input_time
      if (irest == 1) tt = input_time
      return
   endif

   ! Standard formatted or unformatted restart
   form = (ntx == 1 .or. ntx == 5 .or. ntx == 7)
   
   if (form) then
      !        subr amopen(lun,fname,fstat,fform,facc)
      call amopen(lun,inpcrd,'O','F','R')
      read(lun,9008) title1
      
      read(lun,'(a80)') line_test
      if( line_test(6:6) == ' ' ) then ! this is an old, i5 file
        read(line_test,9010) natom,input_time,lremd_values(1)
      elseif( line_test(7:7) == ' ' ) then ! sander 7/8/9/10 large system format...
        read(line_test,9011) natom,input_time,lremd_values(1)
      elseif( line_test(8:8) == ' ' ) then ! Sander 11 - 1 mil+ format
        read(line_test,9012) natom,input_time,lremd_values(1)
      else                   ! assume amber 11 VERY large system format. 10 mil+
        read(line_test,9013,err=666) natom,input_time,lremd_values(1)
      end if
      ! See how many words were in the first line. If there were 3, that means
      ! one of them was the temperature (or pH). If that was the case, allow pH
      ! REMD simulations to set solvph to zero.
      call get_num_tokens(line_test, nwords)
      
      9010 format(i5,2e15.7)
      9011 format(i6,2e15.7)
      9012 format(i7,2e15.7)
      9013 format(i8,2e15.7)

#ifdef MPI
      ! in REM, overwrite temp0 if it's present in inpcrd/restrt
      ! DAN ROE: May want to just use temps in input file, so only 
      !          do this if a restart is requested.
      if(rem == 1 .and. irest==1) then
        if(isgld > 0) then
          if (lremd_values(1) > 0) then
            write(6, '(a)') '| Correcting stagid using value from inpcrd file'
            stagid = lremd_values(1)
          else
            stagid = 0
          end if
        else if (lremd_values(1) > 0) then
          write(6, '(a)') '| Overwriting temp0 from mdin with temp0 from &
                          &netcdf inpcrd file'
          temp0 = lremd_values(1)
        end if
      else if(rem == 4 .and. nwords == 3 .and. irest==1) then
         write(6, '(a)') '| Overwriting solvph from mdin with solvph from &
                          &netcdf inpcrd file'
         solvph = lremd_values(1)
      else if(rem == 5 .and. nwords == 3 .and. irest==1) then
         write(6, '(a)') '| Overwriting solve from mdin with solve from &
                          &netcdf inpcrd file'
         solve = lremd_values(1)
      end if
#endif
      ! If restarting, set the time from input_time
      if (irest == 1) tt = input_time
      
      if(natom == nr) then
         read(lun,9028,end=667,err=668) (x(i),i=1,natom*3)
      else
         write(6,9118)
         call mexit(6, 1)
      end if
      
      if(ntx == 1) then
         do i = 1,nr3
            v(i) = 0.d0
         end do
         if (writeflag) then
            write(6,9008) title1
            write(6,9009) input_time
         end if
         close(lun, iostat=ier)
         return
      end if
      
      !     ----- LOAD THE VELOCITY -----
      
      read(lun,9028,end=669,err=670) (v(i),i=1,nr3)
      if (writeflag) then
         write(6,9008) title1
         write(6,9009) input_time
      end if
      close(lun, iostat=ier)
      return
      
      !     ----- BINARY READING -----
      
   else
      call amopen(lun,inpcrd,'O','U','R')
      if(ntx == 2)then
         read(lun) title1
         read(lun) natom
         if(natom /= nr) then
            write(6,9118)
            call mexit(6, 1)
         end if
         read(lun,end=1000,err=1000) (x(i),i = 1,nr3)
         do i = 1,nr3
            v(i) = 0.d0
         end do
         write(6,9008) title1
         write(6,9009) input_time
         close(lun, iostat=ier)
         return
      end if
      
      if(ntx == 3)then
         read(lun) title1
         read(lun) natom
         if(natom /= nr) then
            write(6,9118)
            call mexit(6, 1)
         end if
         read(lun,end=1000,err=1000) (x(i),i = 1,nr3)
         read(lun) (f(i),i = 1,nr3)
         write(6,9008) title1
         write(6,9009) input_time
         close(lun, iostat=ier)
         return
      end if
      
      read(lun) title1
      read(lun) natom,input_time
      if(natom /= nr) then
         write(6,9118)
         call mexit(6, 1)
      end if
      if (irest == 1) tt = input_time
      read(lun,end=1000,err=1000) (x(i),i = 1,nr3)
      read(lun,end=1010,err=1010) (v(i),i = 1,nr3)
      if(ntx < 6) then
         write(6,9008) title1
         write(6,9009) input_time
         close(lun, iostat=ier)
         return
      end if
      
      !     READ(lun,end=1020,err=1020) a,b,c,alpha,beta,gamma
      do i=1,6
         ifld(i)=3
         fvar(i)=0.0d0
      end do
      ifld(7)=0
      ihol(1)=' '
      call rfree(ifld,ihol,ivar,fvar,lun,6)
      if((fvar(4) /= 0).or.(fvar(5) /= 0).or.(fvar(6) /= 0)) then
         alpha=fvar(4)
         beta=fvar(5)
         gamma=fvar(6)
      else
         !     -- defaults:
         alpha=90.0d0
         beta=90.0d0
         gamma=90.0d0
      end if
      a=fvar(1)
      b=fvar(2)
      c=fvar(3)
      box(1) = a
      box(2) = b
      box(3) = c
      if ( alpha < 1.d0 ) then
         write(6,'(/,a)') 'EWALD: BAD BOX PARAMETERS in inpcrd!'
         call mexit(6, 1)
      end if
      call fill_ucell(a,b,c,alpha,beta,gamma)
      write(6,9129) a,b,c,alpha,beta,gamma
   end if
   
   write(6,9008) title1
   write(6,9009) input_time
   close(lun, iostat=ier)
   return

666 continue
   ! We hit here if we couldn't find NATOM on the 2nd line -- it's not an inpcrd
   write(6, '(a)')  'ERROR: I could not find the number of atoms or the time on'
   write(6, '(3a)') '       the second line of your inpcrd file [', trim(inpcrd), &
                    ']. Bad INPCRD file!'
   call mexit(6, 1)

667 continue
   ! We hit here if the coordinate section was truncated
   write(6, '(a)') 'ERROR: I could not find enough coordinates in ', trim(inpcrd)
   call mexit(6, 1)

669 continue
   ! We hit here if the velocity section was truncated
   write(6, '(a)') 'ERROR: I could not find enough velocities in ', trim(inpcrd)
   call mexit(6, 1)

668 continue
   ! We hit here if we had an error reading floating points into the crd array

   ! Calculate which line the error occurred on. The first 2 lines are title,
   ! and info (NATOM, etc.). i-1 is the last successfully-read float, and there
   ! are 6 floats on each line, starting on line 3
   line_num = 3 + (i-1) / 6
   go to 671

670 continue
   ! We hit here if we had an error reading floating points into the vel array

   ! Calculate which line the error occurred on. Same as above for coordinates,
   ! but adjusted for the fact that all coordinates are present. There are 6 #s
   ! (2 atoms) on each line, with an odd atom getting its own line at the end
   line_num = natom / 2 + mod(natom, 2) + 3 + (i-1) / 6
   go to 671

671 continue
   write(6, '(2a)') 'ERROR: Problem reading coordinates or velocities from ', &
                    trim(inpcrd)
   write(6, '()')
   write(6, '(a,i5,a)') 'I could not understand line ', line_num, ' :'

! I want to print the offending line. However, there's no convenient way to
! do this in Fortran. So I will rewind the whole file, read ic+2 lines, then
! print the ic+3'th line to unit 6
   rewind(INPCRD_UNIT)
   do i = 1, line_num
      read(INPCRD_UNIT, '(a80)') line_test
   end do

   write(6, '(a)') line_test
   write(6, '()')

   ! Not periodic -- that is caught in peek_ewald_inpcrd in ewald_setup.f
   ! before here
   call check_inpcrd_overflow(line_test, .false.)
   close (INPCRD_UNIT)
   call mexit(6, 1)

1000 continue
   write(6,'(a,a)') 'FATAL: Could not read coords from ',trim(inpcrd)
   call mexit(6, 1)
   1010 continue
   write(6,'(a,a)') 'FATAL: Could not read velocities from ',trim(inpcrd)
   call mexit(6, 1)
!  1020 continue
!  write(6,'(a,a)') 'FATAL: Could not read BOX from ',trim(inpcrd)
!  call mexit(6, 1)
   
   9008 format(a80)
   9009 format(t2,'begin time read from input coords =', &
         f10.3,' ps'/)
   9028 format(6f12.7)
   
   9108 format &
         (/80("-")/,'   3.  ATOMIC COORDINATES AND VELOCITIES',/80("-")/)
   9118 format(/2x,'FATAL: NATOM mismatch in coord and ', &
         'topology files')
   9129 format(t2,'NEW EWALD BOX PARAMETERS from inpcrd file:', &
         /5x,'A     =',f10.5,'  B    =',f10.5,'  C     =',f10.5,/, &
         /5x,'ALPHA =',f10.5,'  BETA =',f10.5,'  GAMMA =',f10.5,/)
end subroutine getcor 
