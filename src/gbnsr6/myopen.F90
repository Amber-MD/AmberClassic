#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ A simple file open routine that mimicks AMBER open
subroutine myopen(ifilenum,fname,fstat,fform,facc)

   implicit none

   character(len=* ) fname !fixme to have safer buffer
   character(len=1 ) fstat
   character(len=1 ) fform
   character(len=1 ) facc

   character(len=7 ) mystat
   character(len=11) myform
   integer ierror,ifilenum

   select case (fstat)
      case ('N')
         mystat = 'NEW'
      case ('O')
         mystat = 'OLD'
      case ('R')
         mystat = 'REPLACE'
      case ('U')
         mystat = 'UNKNOWN'
      case default
         write(6,'(a)') 'myopen: wrong fstat argument'; call mexit(6,1)
   end select
   select case (fform)
      case ('F')
         myform = 'FORMATTED'
      case ('U')
         myform = 'UNFORMATTED'
      case default
         write(6,'(a)') 'myopen: wrong fform argument'; call mexit(6,1)
   end select
   select case (facc)
      case ('A')
         open(unit=ifilenum,file=fname,status=mystat,form=myform, &
            position='APPEND',iostat=ierror)
      case ('R')
         open(unit=ifilenum,file=fname,status=mystat,form=myform, &
            action='READ',iostat=ierror)
      case ('W')
         open(unit=ifilenum,file=fname,status=mystat,form=myform, &
            action='READWRITE',iostat=ierror)
      case default
         write(6,'(a)') 'myopen: wrong facc argument'; call mexit(6,1)
   end select
   if (ierror /= 0) then
      if (ifilenum /= 6) then
         write(6,'(/,2x,a,i4,a,a)') 'Unit ', ifilenum, &
               ' Error on OPEN: ',fname
         close(unit=6)
      end if
      write(0,'(/,2x,a,i4,a,a)') 'Unit ', ifilenum, &
            ' Error on OPEN: ',fname
      call mexit(6, 1)
   end if
   rewind(ifilenum)!this is presumably wrong on APPEND(?)


end subroutine
