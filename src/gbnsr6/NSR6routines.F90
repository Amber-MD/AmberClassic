! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"
#include "timer.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine calc_NSR6(acrd,xm,ym,zm,xmymzm,nbndx,nbndy,nbndz,iepsavx,iepsavy,&
                    iepsavz,gox,goy,goz,h,fedgex,fedgey,fedgez,sasopt,natom)

   use solvent_accessibility, only : dprob, radi, arccrd, narcdot, ntri
   use genborn, only : reff, onereff !, B 
   implicit none

#  include "pb_constants.h"
#  include "md.h"

   ! Passed variables
   _REAL_ acrd(3,*),secx(3,nbndx),secy(3,nbndy),secz(3,nbndz)
   _REAL_  dsx(nbndx),dsy(nbndy),dsz(nbndz)
   _REAL_  rpx(3,nbndx),rpy(3,nbndy),rpz(3,nbndz)
   integer xm,ym,zm,xmymzm,nbndx,nbndy,nbndz
   !integer iepsav(4,xmymzm), iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   integer iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   _REAL_ gox,goy,goz,h
   _REAL_ fedgex(xmymzm), fedgey(xmymzm), fedgez(xmymzm)
   integer sasopt

   ! Local variables

   integer i, j, k, iatm, ip, ii
   integer cnt_dot, rnt_dot, dim_dot, tri_dot
   _REAL_ x(3), crd(3)
   _REAL_ rn(1:3), rsphere, dr, r1, r2, r3, h2, hh
   _REAL_ ds1, total_s1, ds2, total_s2
   _REAL_ ds, total_s , cnt_s, rnt_s, dim_s, tri_s
   _REAL_ dss, tss, rx, ry, rz, ess, e0, e1

  ! added by SI
   integer natom
   real*8::sec(3) ! Intersection of grid and surface
   real*8::rb(3),rp(3),rb_dot_rp,rb2,rb6,drb,total_rb
   real*8::atom(3), rinv, rp1,rcut 
   _REAL_ mytime
   call wallclock(mytime)

!   open(100,file='testfile')
!   open(101,file='pp2new.rinv')
 !   write(100,*) natom
 !   write(100,*) mytime
   ! mjhsieh: warning eliminator
   rsphere = -1d0
   h2 = h*h
   hh = HALF*h
   rcut = 7.0d0
   onereff = 0.0
!   print*, 'crdcrd',' zone'

   cnt_dot = 0; rnt_dot = 0; dim_dot = 0; tri_dot = 0
   total_s = ZERO; cnt_s = ZERO; rnt_s = ZERO
   dim_s = ZERO; tri_s = ZERO
   dss = ZERO; tss = ZERO; ess = ZERO; e0 = ZERO; e1 = ZERO
   ds1 = ZERO; ds2 = ZERO; total_s1 = ZERO; total_s2 = ZERO
   secx = ZERO; secy = ZERO; secz = ZERO
   dsx = ZERO; dsy = ZERO; dsz = ZERO
   rpx = ZERO; rpy = ZERO; rpz = ZERO

   do ip = 1, nbndx
      i = iepsavx(1,ip); j = iepsavx(2,ip); k = iepsavx(3,ip); iatm = iepsavx(4,ip)
      crd(1) = gox + h*i + hh; crd(2) = goy + h*j; crd(3) = goz + h*k
      secx(1,ip) = gox + h*i + fedgex(ip)*h; secx(2,ip) = goy + h*j; secx(3,ip) = goz + h*k
      if ( iatm == 0 ) then
         write(6,'(a)') 'PBMD FATAL ERROR: cannot find owner of boundary grid points'
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
         rpx(1,ip) = secx(1,ip) - x(1)
         rpx(2,ip) = secx(2,ip) - x(2)
         rpx(3,ip) = secx(3,ip) - x(3)
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
         rpx(1,ip) = x(1) - secx(1,ip)
         rpx(2,ip) = x(2) - secx(2,ip)
         rpx(3,ip) = x(3) - secx(3,ip)
      end if
      dr = abs(rn(1))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(1)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      dsx(ip) = ds
      total_s = total_s + ds
  enddo

   do ip = 1, nbndy
      i = iepsavy(1,ip); j = iepsavy(2,ip); k = iepsavy(3,ip); iatm = iepsavy(4,ip)
      crd(1) = gox + h*i; crd(2) = goy + h*j + hh; crd(3) = goz + h*k
      secy(1,ip) = gox + h*i; secy(2,ip) = goy + h*j+ fedgey(ip)*h ; secy(3,ip) = goz + h*k
      if ( iatm == 0 ) then
         write(6,'(a)') 'PBMD FATAL ERROR: cannot find owner of boundary grid points'
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
         rpy(1,ip) = secy(1,ip) - x(1)
         rpy(2,ip) = secy(2,ip) - x(2)
         rpy(3,ip) = secy(3,ip) - x(3)
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
         rpy(1,ip) = x(1) - secy(1,ip)
         rpy(2,ip) = x(2) - secy(2,ip)
         rpy(3,ip) = x(3) - secy(3,ip)
      end if
      dr = abs(rn(2))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(1)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      dsy(ip) = ds
      total_s = total_s + ds
  enddo

  do ip = 1, nbndz
      i = iepsavz(1,ip); j = iepsavz(2,ip); k = iepsavz(3,ip); iatm = iepsavz(4,ip)
      crd(1) = gox + h*i; crd(2) = goy + h*j; crd(3) = goz + h*k + hh
      secz(1,ip) = gox + h*i; secz(2,ip) = goy + h*j ; secz(3,ip) = goz + h*k + fedgez(ip)*h
      if ( iatm == 0 ) then
         write(6,'(a)') 'PBMD FATAL ERROR: cannot find owner of boundary grid points'
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
         rpz(1,ip) = secz(1,ip) - x(1)
         rpz(2,ip) = secz(2,ip) - x(2)
         rpz(3,ip) = secz(3,ip) - x(3)
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
         rpz(1,ip) = x(1) - secz(1,ip)
         rpz(2,ip) = x(2) - secz(2,ip)
         rpz(3,ip) = x(3) - secz(3,ip)
      end if
      dr = abs(rn(3))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(1)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      dsz(ip) = ds
      total_s = total_s + ds
  enddo


  do ii=1,natom
   cnt_dot = 0; rnt_dot = 0; dim_dot = 0; tri_dot = 0
   total_s = ZERO; cnt_s = ZERO; rnt_s = ZERO
   dim_s = ZERO; tri_s = ZERO
   dss = ZERO; tss = ZERO; ess = ZERO; e0 = ZERO; e1 = ZERO
   ds1 = ZERO; ds2 = ZERO; total_s1 = ZERO; total_s2 = ZERO
   atom(1:3) = acrd(1:3,ii)
   total_rb = ZERO; drb = ZERO
   do ip = 1, nbndx
      rb(1) = secx(1,ip)-atom(1)
      rb(2) = secx(2,ip)-atom(2)
      rb(3) = secx(3,ip)-atom(3)         
      rb2 = rb(1)*rb(1) + rb(2)*rb(2) + rb(3)*rb(3)
     ! if (rb2>rcut*rcut) cycle
      rb_dot_rp = rb(1)*rpx(1,ip) + rb(2)*rpx(2,ip) + rb(3)*rpx(3,ip)
            rb6 = rb2*rb2*rb2
             rp1 = sqrt(rpx(1,ip)*rpx(1,ip) + rpx(2,ip)*rpx(2,ip) + rpx(3,ip)*rpx(3,ip))
            drb = rb_dot_rp/(rb6*rp1) * dsx(ip)
       total_rb = total_rb + drb 
   end do !nbndx
   do ip = 1, nbndy
      rb(1) = secy(1,ip)-atom(1)
      rb(2) = secy(2,ip)-atom(2)
      rb(3) = secy(3,ip)-atom(3)
      rb2 = rb(1)*rb(1) + rb(2)*rb(2) + rb(3)*rb(3)
     ! if (rb2>rcut*rcut) cycle
      rb_dot_rp = rb(1)*rpy(1,ip) + rb(2)*rpy(2,ip) + rb(3)*rpy(3,ip)
            rb6 = rb2*rb2*rb2
             rp1 = sqrt(rpy(1,ip)*rpy(1,ip) + rpy(2,ip)*rpy(2,ip) + rpy(3,ip)*rpy(3,ip))

            drb = rb_dot_rp/(rb6*rp1) * dsy(ip)
       total_rb = total_rb + drb 
   end do !nbndy
   do ip = 1, nbndz
      rb(1) = secz(1,ip)-atom(1)
      rb(2) = secz(2,ip)-atom(2)
      rb(3) = secz(3,ip)-atom(3)
      rb2 = rb(1)*rb(1) + rb(2)*rb(2) + rb(3)*rb(3)
     ! if (rb2>rcut*rcut) cycle
      rb_dot_rp = rb(1)*rpz(1,ip) + rb(2)*rpz(2,ip) + rb(3)*rpz(3,ip)
            rb6 = rb2*rb2*rb2
             rp1 = sqrt(rpz(1,ip)*rpz(1,ip) + rpz(2,ip)*rpz(2,ip) + rpz(3,ip)*rpz(3,ip))
            drb = rb_dot_rp/(rb6*rp1) * dsz(ip)
       total_rb = total_rb + drb 
   end do !nbndz

   total_s = total_s*h2
   total_rb = total_rb * h2 / (FOURPI)
   rinv = (total_rb) ** (1./3.)
    onereff(ii) = rinv
   enddo ! ii : over n   
   call wallclock(mytime)
  ! write(100,*) mytime
end subroutine calc_NSR6
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine calc_chunk_inv(acrd,xm,ym,zm,xmymzm,nbndx,nbndy,nbndz,iepsavx,iepsavy,&
                    iepsavz,gox,goy,goz,h,fedgex,fedgey,fedgez,sasopt,natom, atmid, rinvi )

   use solvent_accessibility, only : dprob, radi, arccrd, narcdot, ntri
   !use genborn, only : reff, onereff !, B 
   implicit none

#  include "pb_constants.h"
#  include "md.h"

   ! Passed variables
   _REAL_ acrd(3,*),secx(3,nbndx),secy(3,nbndy),secz(3,nbndz)
   _REAL_  dsx(nbndx),dsy(nbndy),dsz(nbndz)
   _REAL_  rpx(3,nbndx),rpy(3,nbndy),rpz(3,nbndz)
   integer xm,ym,zm,xmymzm,nbndx,nbndy,nbndz
   !integer iepsav(4,xmymzm), iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   integer iepsavx(4,xmymzm), iepsavy(4,xmymzm), iepsavz(4,xmymzm)
   _REAL_ gox,goy,goz,h
   _REAL_ fedgex(xmymzm), fedgey(xmymzm), fedgez(xmymzm)
   integer sasopt
   integer atmid
   _REAL_ rinvi

   ! Local variables

   integer i, j, k, iatm, ip, ii
   integer cnt_dot, rnt_dot, dim_dot, tri_dot
   _REAL_ x(3), crd(3)
   _REAL_ rn(1:3), rsphere, dr, r1, r2, r3, h2, hh
   _REAL_ ds1, total_s1, ds2, total_s2
   _REAL_ ds, total_s , cnt_s, rnt_s, dim_s, tri_s
   _REAL_ dss, tss, rx, ry, rz, ess, e0, e1

  ! added by SI
   integer natom
   real*8::sec(3) ! Intersection of grid and surface
   real*8::rb(3),rp(3),rb_dot_rp,rb2,rb6,drb,total_rb
   real*8::atom(3), rinv, rp1,rcut
   _REAL_ mytime
   call wallclock(mytime)

!   open(100,file='testfile')
!   open(101,file='pp2new.rinv')
!    write(100,*) natom
!    write(100,*) mytime
   ! mjhsieh: warning eliminator
   rsphere = -1d0
   h2 = h*h
   hh = HALF*h
   rcut = 7.0d0
   rinvi = 0.0
!   print*, 'crdcrd',' zone'

   cnt_dot = 0; rnt_dot = 0; dim_dot = 0; tri_dot = 0
   total_s = ZERO; cnt_s = ZERO; rnt_s = ZERO
   dim_s = ZERO; tri_s = ZERO
   dss = ZERO; tss = ZERO; ess = ZERO; e0 = ZERO; e1 = ZERO
   ds1 = ZERO; ds2 = ZERO; total_s1 = ZERO; total_s2 = ZERO
   secx = ZERO; secy = ZERO; secz = ZERO
   dsx = ZERO; dsy = ZERO; dsz = ZERO
   rpx = ZERO; rpy = ZERO; rpz = ZERO

   do ip = 1, nbndx
      i = iepsavx(1,ip); j = iepsavx(2,ip); k = iepsavx(3,ip); iatm = iepsavx(4,ip)
      crd(1) = gox + h*i + hh; crd(2) = goy + h*j; crd(3) = goz + h*k
      secx(1,ip) = gox + h*i + fedgex(ip)*h; secx(2,ip) = goy + h*j; secx(3,ip) = goz + h*k
      if ( iatm == 0 ) then
         write(6,'(a)') 'PBMD FATAL ERROR: cannot find owner of boundary grid points'
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
         rpx(1,ip) = secx(1,ip) - x(1)
         rpx(2,ip) = secx(2,ip) - x(2)
         rpx(3,ip) = secx(3,ip) - x(3)
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
         rpx(1,ip) = x(1) - secx(1,ip)
         rpx(2,ip) = x(2) - secx(2,ip)
         rpx(3,ip) = x(3) - secx(3,ip)
      end if
      dr = abs(rn(1))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(1)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      dsx(ip) = ds
      total_s = total_s + ds
   enddo

   do ip = 1, nbndy
      i = iepsavy(1,ip); j = iepsavy(2,ip); k = iepsavy(3,ip); iatm = iepsavy(4,ip)
      crd(1) = gox + h*i; crd(2) = goy + h*j + hh; crd(3) = goz + h*k
      secy(1,ip) = gox + h*i; secy(2,ip) = goy + h*j+ fedgey(ip)*h ; secy(3,ip) = goz + h*k
      if ( iatm == 0 ) then
         write(6,'(a)') 'PBMD FATAL ERROR: cannot find owner of boundary grid points'
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
         rpy(1,ip) = secy(1,ip) - x(1)
         rpy(2,ip) = secy(2,ip) - x(2)
         rpy(3,ip) = secy(3,ip) - x(3)
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
         rpy(1,ip) = x(1) - secy(1,ip)
         rpy(2,ip) = x(2) - secy(2,ip)
         rpy(3,ip) = x(3) - secy(3,ip)
      end if
      dr = abs(rn(2))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(1)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      dsy(ip) = ds
      total_s = total_s + ds
   enddo

   do ip = 1, nbndz
      i = iepsavz(1,ip); j = iepsavz(2,ip); k = iepsavz(3,ip); iatm = iepsavz(4,ip)
      crd(1) = gox + h*i; crd(2) = goy + h*j; crd(3) = goz + h*k + hh
      secz(1,ip) = gox + h*i; secz(2,ip) = goy + h*j ; secz(3,ip) = goz + h*k + fedgez(ip)*h
      if ( iatm == 0 ) then
         write(6,'(a)') 'PBMD FATAL ERROR: cannot find owner of boundary grid points'
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
         rn = crd - x
         rpz(1,ip) = secz(1,ip) - x(1)
         rpz(2,ip) = secz(2,ip) - x(2)
         rpz(3,ip) = secz(3,ip) - x(3)
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn =  x - crd
         rpz(1,ip) = x(1) - secz(1,ip)
         rpz(2,ip) = x(2) - secz(2,ip)
         rpz(3,ip) = x(3) - secz(3,ip)
      end if
      dr = abs(rn(3))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1
      ds1 = dr*rsphere**2*r3
      ! ds = dr/rsphere
!     ds2 = dr*rsphere**2*r3*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(1)**2*r2))
!     total_s1 = total_s1 + ds1
!     total_s2 = total_s2 + ds2
      ds = ds1
      dsz(ip) = ds
      total_s = total_s + ds
   enddo

   ii=atmid ! this is just for one atom
   cnt_dot = 0; rnt_dot = 0; dim_dot = 0; tri_dot = 0
   total_s = ZERO; cnt_s = ZERO; rnt_s = ZERO
   dim_s = ZERO; tri_s = ZERO
   dss = ZERO; tss = ZERO; ess = ZERO; e0 = ZERO; e1 = ZERO
   ds1 = ZERO; ds2 = ZERO; total_s1 = ZERO; total_s2 = ZERO
   atom(1:3) = acrd(1:3,ii)
   total_rb = ZERO; drb = ZERO
   do ip = 1, nbndx
      rb(1) = secx(1,ip)-atom(1)
      rb(2) = secx(2,ip)-atom(2)
      rb(3) = secx(3,ip)-atom(3)
      rb2 = rb(1)*rb(1) + rb(2)*rb(2) + rb(3)*rb(3)
   !   if (rb2>rcut*rcut) cycle
      rb_dot_rp = rb(1)*rpx(1,ip) + rb(2)*rpx(2,ip) + rb(3)*rpx(3,ip)
            rb6 = rb2*rb2*rb2
             rp1 = sqrt(rpx(1,ip)*rpx(1,ip) + rpx(2,ip)*rpx(2,ip) + rpx(3,ip)*rpx(3,ip))
            drb = rb_dot_rp/(rb6*rp1) * dsx(ip)
       total_rb = total_rb + drb
   end do !nbndx
   do ip = 1, nbndy
      rb(1) = secy(1,ip)-atom(1)
      rb(2) = secy(2,ip)-atom(2)
      rb(3) = secy(3,ip)-atom(3)
      rb2 = rb(1)*rb(1) + rb(2)*rb(2) + rb(3)*rb(3)
   !   if (rb2>rcut*rcut) cycle
      rb_dot_rp = rb(1)*rpy(1,ip) + rb(2)*rpy(2,ip) + rb(3)*rpy(3,ip)
            rb6 = rb2*rb2*rb2
             rp1 = sqrt(rpy(1,ip)*rpy(1,ip) + rpy(2,ip)*rpy(2,ip) + rpy(3,ip)*rpy(3,ip))

            drb = rb_dot_rp/(rb6*rp1) * dsy(ip)
       total_rb = total_rb + drb
   end do !nbndy
   do ip = 1, nbndz
      rb(1) = secz(1,ip)-atom(1)
      rb(2) = secz(2,ip)-atom(2)
      rb(3) = secz(3,ip)-atom(3)
      rb2 = rb(1)*rb(1) + rb(2)*rb(2) + rb(3)*rb(3)
   !   if (rb2>rcut*rcut) cycle
      rb_dot_rp = rb(1)*rpz(1,ip) + rb(2)*rpz(2,ip) + rb(3)*rpz(3,ip)
            rb6 = rb2*rb2*rb2
             rp1 = sqrt(rpz(1,ip)*rpz(1,ip) + rpz(2,ip)*rpz(2,ip) + rpz(3,ip)*rpz(3,ip))
            drb = rb_dot_rp/(rb6*rp1) * dsz(ip)
       total_rb = total_rb + drb
   end do !nbndz

   total_s = total_s*h2
   total_rb = total_rb * h2 / (FOURPI)
   !rinvi = (total_rb) ** (1./3.)
   rinvi = total_rb


end subroutine calc_chunk_inv





