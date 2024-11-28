#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rfree here]
subroutine rfree(ifld,ihol,ivar,fvar,in,iout)
   
   implicit none
   
   !     Author:  George Seibel
   
   !     This is a free format reader for mixed Hollerith and numeric data.
   !     This code is a complete re-write of the old Amber rfree, and is
   !     now machine independent ANSI fortran 77.
   !     Rev 02-May-89: changed return on EOF back to stop. (Edit no longer
   !                    needs this bogus feature.)
   !     Rev 14-Mar-89: add else to elseif check on ifld()
   !     Rev 01-Mar-89: initialize ierr to 0
   !     Rev 22-Feb-89: fixed bug in ifld() interpretation
   !     Rev 20-Feb-89: changed stop on EOF to return
   !     Rev 20-Jan-92: made ch2int() ebcdic-capable (BR)
   
   !     PARAMETERS:
   
   integer lenbuf
   !        ... length of character buffer for input line
   parameter (lenbuf=132)
   
   !     INPUT: Local variables
   
   integer ifld(20)  !XP: only rgroup called with dimension 20
   !        ... code for field types to be read:  1 = Hollerith (4 byte int)
   !            2 = integer   3 = float  0 = no more fields
   integer in, iout
   !        ... input and output logical unit numbers
   
   !     OUTPUT:
   
   character(len=4) ihol(*)
   !        ... extracted Hollerith data (4 byte max for Amber)
   integer ivar(20)
   !        ... extracted integer data
   _REAL_  fvar(20)
   !        ... extracted floating pt data
   
   !     INTERNAL:
   
   character(len=lenbuf) buf, token
   !        ... input line buffer,  temp for undecoded tokens
   character(len=4) blank
   !        ... just 4 bytes of space
   integer nvar
   !        ... number of variables to read
   integer ntoken, nint, nhol, nflt
   !        ... counters for tokens, int, hol, and real variables
   logical inword
   !        ... true if tokenizer is in a word
   integer ipt, i
   !        ... pointer into char buffer, loop index
   integer ibeg, iend,lenstr
   !        ... buf indices = beginning and end of current token
   integer ival, ierr
   !        ... temp for integer values, error return for ch2int()
   
   
   ierr = 0
   nvar = 0
   ntoken = 0
   nint = 0
   nhol = 0
   nflt = 0
   blank = ' '
   ibeg = 1
   iend = lenbuf
   inword = .false.
   
   !     --- initialize the output arrays ---
   
   do 100 i = 1, 80
      if (ifld(i) <= 0) goto 110
      nvar = nvar+1
      if (ifld(i) == 1) then
         nhol = nhol + 1
         read(blank,'(a4)') ihol(nhol)
      else if (ifld(i) == 2) then
         nint = nint + 1
         ivar(nint) = 0
      else if (ifld(i) == 3) then
         nflt = nflt + 1
         fvar(nflt) = 0.0d0
      else
         write(iout,'(5x,a)') 'rfree: bogus ifld()'
         call mexit(iout, 1)
      end if
   100 continue
   
   110 continue
   
   !     --- read entire line into character buffer ---
   
   read(in,'(a)',end=1000) buf
   
   !     --- tokenize buf using any whitespace as delimitter ---
   
   nint = 0
   nhol = 0
   nflt = 0
   do 200 ipt = 1, lenbuf
      if (ntoken >= nvar) return
      if (.not. inword) then
         !               --- look for start of word = non-whitespace --
         if (buf(ipt:ipt) > ' ') then
            inword = .true.
            ibeg = ipt
         end if
      else
         !               --- look for end of word = whitespace or end of buf ---
         if (buf(ipt:ipt) <= ' ' .or. ipt >= lenbuf) then
            inword = .false.
            ntoken = ntoken + 1
            iend = ipt
            token = buf(ibeg:iend)
            lenstr = iend - ibeg
            
            !                    --- decode according to ifld() ---
            
            if (ifld(ntoken) == 1) then
               !                         --- Hollerith was a Great Man ---
               nhol = nhol + 1
               read(token,'(a4)') ihol(nhol)
            else if (ifld(ntoken) == 2) then
               !                         --- Integer Field ---
               nint = nint + 1
               call ch2int(lenstr,token,ival,ierr)
               if (ierr /= 0) goto 900
               ivar(nint) = ival
            else if (ifld(ntoken) == 3) then
               !                         --- Floating Point field ---
               nflt = nflt + 1
               if (index(token,'.') > 0) then
                  !                              --- if decimal pt, use internal read ---
                  read(token,'(f20.0)',err=900) fvar(nflt)
               else
                  !                              --- no decimal, use char to int routine ---
                  call ch2int(lenstr,token,ival,ierr)
                  if (ierr /= 0) goto 900
                  fvar(nflt) = ival
               end if
            end if
         end if  ! (buf(ipt:ipt) <= ' ' .or. ipt >= lenbuf)
      end if  ! (.not. inword)
   200 continue
   return
   900 continue
   
   !     --- token could not be decoded ---
   
   write(iout,'(/5x,a,i3,i3,a,/,a)')'rfree: Error decoding variable', &
         ntoken, ifld(ntoken), ' from:', buf(1:iend)
      write(iout,'(/5x,a)')'this indicates that your input contains ',&
   ' incorrect information'
      write(iout,'(/5x,a,i3,a)') 'field ',ntoken,' was supposed to ',&
   ' have a (1=character, 2=integer, 3=decimal) value'
	
   call mexit(iout, 1)
   1000 continue
   
   !     --- hit EOF ---
   
   write(iout,'(/5x,a,i3)') 'rfree: End of file on unit ', in
   call mexit(iout, 1)
end subroutine rfree 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ch2int here]
subroutine ch2int(lenstr, string,ival,ierr)
   
   !     converts character representations of cardinal numbers to integer
   !     Author:  George Seibel
   !     Rev 01-Mar-89: initialize ierr to 0
   !     Initial Rev: 19-Dec-87
   !     implicit none
   
   !     INPUT:
   implicit none
   integer lenstr
   !        ... length of string
   character string*(*)
   !        ... character representation of legitimate integer
   
   !     OUTPUT:
   
   integer ival
   !        ... the integer result
   integer ierr
   !        ... returned as one if any decoding error, zero if ok
   
   !     INTERNAL:
   
   integer i, j, num, ifirst, last, idec
   logical isneg
   
   ierr = 0
   
   !     --- look for minus sign ---
   
   isneg = (index(string,'-') > 0)
   
   !     --- find first and last numeric character ---
   
   last = lenstr
   do 200 i = 1, lenstr
      if (string(i:i) >= '0'.and. string(i:i) <= '9') then
         ifirst = i
         do 100 j = ifirst, lenstr
            if (string(i:i) < '0'.or.string(i:i) > '9') then
               last = j - 1
               goto 300
            end if
         100 continue
         goto 300
      end if
   200 continue
   
   !     --- no numerics found - error return ---
   
   ierr = 1
   return
   
   !     --- crunch the number ---
   
   300 continue
   num = 0
   idec = 0
   do 400 i = last, ifirst, -1
      num = num + (ichar(string(i:i)) - ichar('0')) * 10**idec
      idec = idec + 1
   400 continue
   if (isneg) num = -num
   ival = num
   return
end subroutine ch2int 
!-----------------------------------------------------------------------
