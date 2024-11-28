#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Amber runmin place holder, only one force call and printing
subroutine runmin(xx,ix,ih,ipairs,x,fg,w,ib,jb,conp, &
                  winv,igrp,skips,ene,r_stack, i_stack,carrms)

   use memory_module

   implicit none

   ! global variables

#  include "pb_constants.h"
#  include "md.h"
#  include "box.h"
#  include "files.h"
#  include "pb_md.h"
#  include "extra.h"

   ! passed variables

   _REAL_   xx(*)
   integer  ix(*), ipairs(*)
   character(len=4) ih(*)
   _REAL_   x(*),fg(*),w(*)
   integer  ib(*),jb(*)
   _REAL_   conp(*),winv(*)
   integer  igrp(*)
   logical  skips(*)
   _REAL_   ene(51)
   _REAL_   r_stack(*)
   integer  i_stack(*)
   _REAL_   carrms

   ! external Functions

   _REAL_ ddot

   ! local variables

   _REAL_ vir(4)
!  _REAL_ dxst, betax, ddspln, dfpr, dxsth
   _REAL_ f, fmin, fnq!, fold, fch, finit, gama, gamden
!  _REAL_ ginit, gmin, gnew, gspln, gsqrd, sbound, step
   _REAL_ sum!, stepch, stmin, work
!  _REAL_ ecopy(51),vircopy(3)
   _REAL_ rms,fndfp!,fdmax,swork
   logical skip,newstr,steep,belly
   integer n,nr,nct,ndfp,nstcyc!,i,nitp,ier
   integer mstcyc,linmin!,iterrs
   integer irsdx,irsdg,iginit,ixopt,igopt,iterc,n_force_calls
!  integer iterfm,iretry,nfbeg,nfopt
!  integer imes
   integer ier
!  integer iprint

   logical do_list_update

   ! Only used by subroutine force when PSANDER is defined.
   ! By default for minimization list updating is controlled the old
   ! way, ie, nbflag is 0 and updating occurs every nsnb steps;
   ! ntnb is the actual variable that is tested by nonbond_list.

   data do_list_update / .true. /

   integer maxlin,mxfcon,kstcyc
   parameter ( maxlin = 10 )  ! maximum number of line searches ?
   parameter ( mxfcon =  4 )  ! maximum force con ?
   parameter ( kstcyc =  4 )  ! number of starting cycles ?
   _REAL_  crits, dxstm, dfpred
   parameter ( dxstm  = TEN_TO_MINUS5 )  ! ?
   parameter ( crits  = TEN_TO_MINUS6 )  ! ?
   parameter ( dfpred = ONE )            ! ? in kcal/mol

   !     ----- EVALUATE SOME CONSTANTS -----

   if (imin /= 5 .and. master) call myopen(MDINFO_UNIT,mdinfo,'U','F','W')
   fmin = 0.0d0
   nr = nrp
   n = 3*nr
   belly = ibelly > 0
   ier = 0
   nct = 0
   if (ntc == 2) nct = nbonh
   if (ntc == 3) nct = nbonh + nbona
   ndfp = n-nct
   if(belly) ndfp = 3*natbel-nct
   ntnb = 1
   fndfp = ndfp
   fnq = sqrt(fndfp)
   rms = 0.0d0
   skip = .false.
   newstr = .false.
   ! determine the number of steepest descent steps
   steep = .false.
   nstcyc = 0
   mstcyc = kstcyc
   if(ntmin == 2) mstcyc = maxcyc
   if(ntmin == 1) mstcyc = ncyc
   if(ntmin > 0) steep = .true.

!  fold = 0.0d0
!  dxst = dx0
   linmin = 0

   !     ----- PARTITION THE WORKING ARRAY -----

   irsdx = n
   irsdg = irsdx+n
   iginit = irsdg+n
   ixopt = iginit+n
   igopt = ixopt+n

   !     ----- SET SOME PARAMETERS TO BEGIN THE CALCULATION -----

   iterc = 0
   n_force_calls = 0
!  iterfm = iterc

   n_force_calls = n_force_calls + 1
   ntnb = 1
   steep = .true.

   !     ----- CALCULATE THE FORCE AND ENERGY -----

   ! APPLY SHAKE TO CONSTRAIN BONDS IF NECESSARY -----
   if(ntc /= 1) then
      write(6,'(a)') 'SHAKE OPTIONS ARE DISABLED'
      call mexit(6,1)
   end if

   ! reset pb-related flags

   if(master)then
      if ( igb == 10 .or. ipb >= 1 ) then
         if ( mod(n_force_calls,npbgrid) == 0 .and. n_force_calls /= maxcyc ) pbgrid = .true.
         if ( mod(n_force_calls,ntpr) == 0 .or. n_force_calls ==maxcyc ) pbprint = .true.
         if ( mod(n_force_calls,nsnbr) == 0 .and. n_force_calls /=maxcyc ) ntnbr = 1
         if ( mod(n_force_calls,nsnba) == 0 .and. n_force_calls /=maxcyc ) ntnba = 1
         npbstep = n_force_calls
      end if
   endif

   call force(xx,ix,ih,ipairs,x,fg,ene,vir,r_stack,i_stack, &
           xx(l96), xx(l97), xx(l98),do_list_update)

   f = ene(1)
   ntnb = 0
   sum = ddot(n,fg,1,fg,1)
   rms = sqrt(sum)/fnq

   !     ----- PRINT THE INTERMEDIATE RESULTS -----
   !           ih(m04) = atom names

   call report_min_progress( n_force_calls, rms, fg, ene, ih(m04) )  

   if (master) then
      write(6,'(//,a)') '  Maximum number of minimization cycles reached.'
      close(unit=MDINFO_UNIT)
   end if

   !     ----- WRITE THE FINAL RESULTS -----

   call report_min_results( n_force_calls, rms, x, fg, ene, ih(m04), xx, ix, ih )  
   carrms = rms


end subroutine runmin

