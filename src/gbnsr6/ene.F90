! <compile=optimized>
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ bond energy place holder
subroutine bond(eb)

   implicit none

   ! passed variables
   _REAL_ eb
   
   eb = 0.d0
   

end subroutine bond 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ angle energy place holder
subroutine angl(eba)

   implicit none

   ! passed variables
   _REAL_ eba

   eba = 0.d0


end subroutine angl
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dihedral energy place holder
!+ it will compute only 1-4ELE B.Aguilar
subroutine ephi(nphiin, ip, jp, kp, lp, icp, cg, x,  ix,ep,enbp,eelp)
!subroutine ephi(ep,enbp,eelp)

   
   use variable_module, only : eps0, epsin

   implicit none
#  include "parms.h"

   ! passed variables
   integer ix(*),ip(*),jp(*),kp(*),lp(*),icp(*)
   _REAL_ cg(*)
   _REAL_ x(*)
   integer nphiin
   _REAL_  ep,enbp,eelp

   ! local variables
   logical skip
   _REAL_ scee0
   _REAL_ intdieli
   integer max190,maxlen
   parameter (max190=190)
   _REAL_ xij(max190),yij(max190),zij(max190)
   _REAL_ ct(max190),cphi(max190),sphi(max190)
   integer lenc,ia1,ia2,ibig,isml,ii,jj,kk,ll
   integer ic,inc, ic0
   integer i3,j3,k3,l3,k3t,l3t
   integer piece,start,end,newnb
   integer jn,istc,ist,nphi
   _REAL_ xa,ya,za,f1,f2,r1,r2,r6,r12,dfn, g
   _REAL_ cgi, cgj, crfac

   enbp = 0.d0
   eelp = 0.d0
   ep   = 0.d0


   start=1
   end=nphiin

   nphi = end
   ist = start -1

   intdieli =eps0/epsin

   4200 continue
  
   maxlen = max190
   skip = (ist+maxlen) > nphi
   if(skip) maxlen = nphi-ist
   if(maxlen <= 0) goto 820

   do jn = 1,maxlen
      i3 = ip(jn+ist)
      l3t = lp(jn+ist)
      l3 = iabs(l3t)
      xij(jn) = x(i3+1)-x(l3+1)
      yij(jn) = x(i3+2)-x(l3+2)
      zij(jn) = x(i3+3)-x(l3+3)
   end do

   do jn = 1,maxlen
        ct(jn) = xij(jn)*xij(jn)+yij(jn)*yij(jn)+zij(jn)*zij(jn)
   end do

   cphi(1:maxlen) = 0.0d0
   sphi(1:maxlen) = 0.0d0 

   do jn = 1,maxlen
        !Check if we should do this 1-4 interaction or not.
        k3t = kp(jn+ist)
        l3t = lp(jn+ist)
        if (k3t < 0 .or. l3t < 0) cycle

        ic0 = icp(jn+ist)
      ! scnb0 = one_scnb(ic0)
        scee0 = one_scee(ic0)
!        scee0 = 1.0d0 / 1.2d0
        i3 = ip(jn+ist)
        l3 = iabs(l3t)
        ii = (i3+3)/3
        jj = (l3+3)/3

        !             ----- CALCULATE THE 14-EEL ENERGY -----
        r2 = 1.0d0/ct(jn)
        r1 = sqrt(r2)

        cgi = cg(ii)
        cgj = cg(jj)

        crfac = r1*intdieli
        g = cgi*cgj*crfac
        sphi(jn) = g*scee0

   end do 
  
   do jn = 1,maxlen
      eelp = eelp+sphi(jn)
   end do

   ist = ist + maxlen
   820 if(.not.skip) goto 4200


end subroutine ephi
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ capwat restraining force place holder
subroutine capwat()

   implicit none


end subroutine capwat
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ position restraining force place holder
subroutine xconst(econ)

   implicit none

   ! passed variables
   _REAL_ econ

   econ = 0.d0


end subroutine xconst
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ belly restraining force place holder
subroutine bellyf()

   implicit none


end subroutine bellyf
