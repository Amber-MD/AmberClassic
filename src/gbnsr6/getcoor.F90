#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Amber restart file reader
subroutine getcor(nr,x,v,f,ntx,box,mytime)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Rewritten by:
   ! Meng-Juei Hsieh, Luo Research Group, UC-Irvine
   !
   ! This is a Amber restart file reader that preserves the behaviors of getcor
   ! for compatibility.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

#  include "files.h"

   ! Passed variable
   integer filenum,nr,ntx
   _REAL_ x(*),v(*),f(*),box(3),mytime
   
   ! Local variables
   logical isascii
   _REAL_ fvar(7),a,b,c,alpha,beta,gamma
   integer i,ivar,ifld(7),ihol(1)
   integer natom,nr3,ier
   character(len=256) testbuffer
   
   nr3 = 3*nr
   filenum = 9
   
   write(6,9108)
   
   isascii = (ntx == 1 .or. ntx == 5 .or. ntx == 7 )
   
   ! open the restart file
   if(isascii)then
      call myopen(filenum,inpcrd,'O','F','R')
      read(filenum,9008) title1
      read(filenum,'(a80)') testbuffer
      if( testbuffer(6:6) == ' ' ) then ! this is an old, i5 file
         read(testbuffer,9010) natom,mytime
      else                             ! assume a new, i6 file
         read(testbuffer,9011) natom,mytime
      end if
      if(natom /= nr) then
         write(6,9118)
         close(filenum, iostat=ier)
         call mexit(6, 1)
      end if
      read(filenum,9028) (x(i),i=1,nr3)
   else
      call myopen(filenum,inpcrd,'O','U','R')
      read(filenum) title1
   endif

   ! read the file according to the option speficied by ntx
   select case(ntx)
      case (1)! ASCII
         ! X is read formatted with no initial velocity information (default)
         do i = 1,nr3
            v(i) = 0.d0
         end do
         write(6,9008) title1
         write(6,9009) mytime
         close(filenum, iostat=ier)
         return
      case (2)! BINARY
         ! X is read unformatted with no initial velocity information
         read(filenum) natom
         if(natom /= nr) then
            write(6,9118)
            close(filenum, iostat=ier)
            call mexit(6, 1)
         end if
         read(filenum,end=1000,err=1000) (x(i),i = 1,nr3)
         do i = 1,nr3
            v(i) = 0.d0
         end do
         write(6,9008) title1
         write(6,9009) mytime
         close(filenum, iostat=ier)
         return
      case (3)! BINARY
         ! X and F are read unformatted.
         read(filenum) natom
         if(natom /= nr) then
            write(6,9118)
            close(filenum, iostat=ier)
            call mexit(6, 1)
         end if
         read(filenum,end=1000,err=1000) (x(i),i = 1,nr3)
         read(filenum) (f(i),i = 1,nr3)
         write(6,9008) title1
         write(6,9009) mytime
         close(filenum, iostat=ier)
         return
      case (4)! BINARY
         ! X and V are read unformatted.
         read(filenum) natom,mytime
         if(natom /= nr) then
            write(6,9118)
            close(filenum, iostat=ier)
            call mexit(6, 1)
         end if
         read(filenum,end=1000,err=1000) (x(i),i = 1,nr3)
         read(filenum,end=1010,err=1010) (v(i),i = 1,nr3)
         write(6,9008) title1
         write(6,9009) mytime
         close(filenum, iostat=ier)
         return
      case (5,7)! ASCII
         ! X and V are read formatted; box information will be read if ntb>0.
         ! The velocity information will only be used if irest=1.
         read(filenum,9028,end=1010) (v(i),i=1,nr3)
         write(6,9008) title1
         write(6,9009) mytime
         close(filenum, iostat=ier)
         return
      case (6) ! BINARY
         ! X, V and BOX(1..3) are read unformatted; in other respects, 
         ! this is the same as option "5".
         read(filenum) natom,mytime
         if(natom /= nr) then
            write(6,9118)
            call mexit(6, 1)
         end if
         read(filenum,end=1000,err=1000) (x(i),i = 1,nr3)
         read(filenum,end=1010,err=1010) (v(i),i = 1,nr3)
         do i=1,6
!           ifld(i)=3
            fvar(i)=0.0d0
         enddo
!        ifld(7)=0
!        ihol(1)=0
!        call rfree(ifld,ihol,ivar,fvar,filenum,6)
         read(filenum,end=1020,err=1020) (fvar(i),i = 1,6)
         if((fvar(4) /= 0).or.(fvar(5) /= 0).or.(fvar(6) /= 0)) then
            alpha=fvar(4)
            beta=fvar(5)
            gamma=fvar(6)
         else
            ! defaults:
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
         !call fill_ucell(a,b,c,alpha,beta,gamma)
         write(6,9129) a,b,c,alpha,beta,gamma
         write(6,9008) title1
         write(6,9009) mytime
         close(filenum, iostat=ier)
         return
      case default
         close(filenum, iostat=ier)
         call mexit(6,1)
   end select
   close(filenum, iostat=ier)

   1000 continue
   write(6,'(a,a)') 'FATAL: Could not read coords from ',inpcrd
   call mexit(6, 1)
   1010 continue
   write(6,'(a,a)') 'FATAL: Could not read velocities from ',inpcrd
   call mexit(6, 1)
   1020 continue
   write(6,'(a,a)') 'FATAL: Could not read BOX from ',inpcrd
   call mexit(6, 1)
   
   9008 format(a80)
   9009 format(t2,'begin time read from input coords =', &
         f10.3,' ps'/)
   9010 format(i5,e15.7)
   9011 format(i6,e15.7)
   9028 format(6f12.7)

   9108 format &
         (/80(1h-)/,'   2.  ATOMIC COORDINATES AND VELOCITIES',/80(1h-)/)
   9118 format(/2x,'FATAL: NATOM mismatch in coord and ', &
         'topology files')
!  9128 format(t2,'NEW BOX DIMENSIONS from inpcrd file:', &
!        /5x,'X =',f10.5,'  Y =',f10.5,'  Z =',f10.5,/)
   9129 format(t2,'NEW EWALD BOX PARAMETERS from inpcrd file:', &
         /5x,'A     =',f10.5,'  B    =',f10.5,'  C     =',f10.5,/, &
         /5x,'ALPHA =',f10.5,'  BETA =',f10.5,'  GAMMA =',f10.5,/)


end subroutine getcor 
