!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine amopen here]
subroutine amopen(lun,fname,fstat,fform,facc)

   !  When this is converted to Fortran 90
   !  these codes can be replaced with module types, eg,
   !       Character(*), public, parameter :: unknown = 'unknown'
   !  facc is used for file appending

   implicit none

   !     INPUT:

   integer lun
   !        ... logical unit number
   character(len=*) fname
   !        ... file name 
   character(len=1) fstat
   !        ... status code: 'N', 'O', 'R', 'U' = new, old, replace, unknown
   character(len=1) fform
   !        ... format code: 'F', 'U' = formatted, unformatted
   character(len=1) facc
   !        ... access code: 'R', 'W', 'A' = read, read/write, append

   !     INTERNAL:

   character(len=7) stat
   !        ... status keyword
   character(len=11) kform
   !        ... form keyword
   character(len=11) pos
   !        ... position keyword
   integer ios 
   !        ... i/o status variable
   character(len=40) errstr

   if (fstat == 'N') then
      stat = 'NEW'
   else if (fstat == 'O') then
      stat = 'OLD'
   else if (fstat == 'R') then
      stat = 'REPLACE'
   else if (fstat == 'U') then
      stat = 'UNKNOWN'
   else
      write(6,'(/,2x,a,i4)') &
            'amopen: bogus fstat, unit ', lun
      call mexit(6, 1)
   end if

   if (fform == 'U') then
      kform = 'UNFORMATTED'
   else if (fform == 'F') then
      kform = 'FORMATTED'
   else
      write(6,'(/,2x,a,i4)') &
            'amopen: bogus fform, unit', lun
      call mexit(6, 1)
   end if

   if (facc == 'A') then
      pos = "APPEND"
   else
      pos = "ASIS" ! default f90
   end if

   if (facc == 'R') then
      errstr = "is missing or unreadable"
   else
      errstr = "exists (use -O to overwrite)"
   end if

   open(unit=lun,file=fname,status=stat,form=kform,iostat=ios,position=pos)

   if (ios /= 0) then
      if (lun == 6) then
         write(0,'(/,2x,a,i4,a,a,a,a)') 'Error opening unit ', &
               lun, ': File "', trim(fname), '" ', errstr
      else
         write(6,'(/,2x,a,i4,a,a,a,a)') 'Error opening unit ', &
               lun, ': File "', trim(fname), '" ', errstr
         close(unit=6)
         write(0,'(/,2x,a,i4,a,a,a,a)') 'Error opening unit ', &
               lun, ': File "', trim(fname), '" ', errstr
      end if
      call mexit(6, 1)
   end if
   if (pos /= "APPEND") rewind(lun)
   return
end subroutine amopen 
