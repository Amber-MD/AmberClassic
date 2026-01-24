! <compile=optimized>
#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine load_debug here]
subroutine load_debug(nf)
   implicit none
   integer nf
#  include "flocntrl.h"
#  include "debug.h"
#  include "md.h"

   integer ifind,j
   namelist /debugf/do_dir,do_rec,do_adj,do_self,do_bond, &
         do_angle,do_ephi,doxconst,do_cap,do_14, &
         do_debugf,neglgdel,zerochg,zerovdw,zerodip, &
         atomn,nranatm, &
         ranseed,chkvir,dumpfrc,rmsfrc,do_tgt,&
         do_pbdir,do_pbnp,do_pbfd
         

   ! default flow control all force routines turned on
   do_dir = 1
   do_rec = 1
   do_adj = 1
   do_self = 1
   do_bond = 1
   do_angle = 1
   do_ephi = 1
   doxconst = 1
   do_cap = 1
   do_14 = 1
   do_tgt = 1
   do_pbdir = 1
   do_pbnp = 1
   do_pbfd = 1

   ! debug off by default
   do_debugf = 0
   neglgdel = 5
   zerochg = 0
   zerovdw = 0
   zerodip = 0
   nranatm = 0
   ranseed = 71277
   chkvir = 0
   dumpfrc = 0
   rmsfrc = 0
   do j = 1,natomn
      atomn(j) = 0
   end do
   rewind(5)
   call mynmlsrc('debugf',nf,ifind)
   if ( ifind /= 0)then
      read(5,nml=debugf,err=190)
   end if

   if ( do_debugf == 0 )then
      do_dir = 1
      do_rec = 1
      do_adj = 1
      do_self = 1
      do_bond = 1
      do_angle = 1
      do_ephi = 1
      doxconst = 1
      do_tgt = 1
      do_cap = 1
      do_pbdir = 1
      do_pbnp = 1
      do_pbfd = 1

      neglgdel = 5
      zerochg = 0
      zerovdw = 0
      zerodip = 0
      nranatm = 0
      ranseed = 71277
      chkvir = 0
      dumpfrc = 0
      rmsfrc = 0
   end if
   return
   190 write(6,'(a)')'Error in &debugf namelist'
   call mexit(6,1)
end subroutine load_debug 
!----------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine debug_frc here]
subroutine debug_frc(xx,ix,ih,ipairs,x,f, &
      cn1,cn2,r_stack,i_stack)
   
   use memory_module

   implicit none
#  include "flocntrl.h"
#  include "debug.h"
#  include "extra.h"
#  include "md.h"
   _REAL_ xx(*),x(*),f(3,*)
   _REAL_ cn1(*),cn2(*),r_stack(*)
   integer ipairs(*),ix(*),i_stack(*)
   character(len=4) ih(*)
   _REAL_ vir(4),ener(51),time,onefac(3),rms
   integer iout7,nstep,nitp,nits,j,k,atomnum
   _REAL_ apfrc(3),del,y
   _REAL_ apavir,apmvir,exavir,exmvir,exvir
   _REAL_ cpfrc1(3,natomn), &
         cpfrc2(3,natomn)
   _REAL_ ene(52),dudv
   _REAL_ intvir(3,3),duda(3,3),apduda(3,3)
   _REAL_ mduda(3,3),mapduda(3,3)
   _REAL_ mvten(3,3),avten(3,3), diff
   integer ranatm(natomn),type

   if ( do_debugf == 0 )return

   ! first call force for analytic result
   if ( master )then
      write(6,'(a)')'DEBUG FORCE!; calling force routine'
   end if
   call get_analfrc(xx,ix,ih,ipairs,x,f, &
         vir,ener,r_stack,i_stack)
   if ( master )then
      write(6,'(a)')'DEBUG FORCE!; back from force routine'
   end if
   do j = 1,51
      ene(j) = ener(j)
   end do

   iout7 = 0
   onefac(1) = 1.d0
   onefac(2) = 1.d0
   onefac(3) = 1.d0
   nstep = 0
   nitp = 0
   nits = 0
   ene(2) = 0.d0
   ene(3) = 0.d0
   ene(1) = ener(23)
   ene(22) = vir(1) + vir(2) + vir(3)
   if ( master )then
      call prntmd(nstep,nitp,nits,t,ene,onefac,iout7,.false.)
   end if

   ! get the delta for numerical force, virial calcs
   del = 1.d0
   do j = 1,neglgdel
      del = del/1.d1
   end do
   ! save copies of analytic forces for user defined and random atoms
   do j = 1,natomn
      atomnum = atomn(j)
      if ( atomnum == 0 )goto 50
      cpfrc1(1,j) = f(1,atomnum)
      cpfrc1(2,j) = f(2,atomnum)
      cpfrc1(3,j) = f(3,atomnum)
   end do
   50 continue
   if ( nranatm > natomn )then
      if ( master) write(6,166)natomn
      166 format(1x,'MAX NUM Random atoms: ',i6)
      call mexit(6,1)
   end if
   if ( nranatm > 0 )then
!     call amrset(ranseed)
      do j = 1,nranatm
!        call amrand(y)
         y = (y*4657+2897)/1523
         y = y - int(y)
         atomnum = y * natom
         ranatm(j) = atomnum
         cpfrc2(1,j) = f(1,atomnum)
         cpfrc2(2,j) = f(2,atomnum)
         cpfrc2(3,j) = f(3,atomnum)
      end do
   end if
   ! now do user defined atoms
   if ( atomn(1) > 0 )then
      if ( master )write(6,168)
      if ( master ) &
            write(6,'(a)')'----------------------------------------------'
   end if
   do j = 1,natomn
      atomnum = atomn(j)
      if ( atomnum == 0 )goto 100
      call get_numfrc(xx,ix,ih,ipairs,x,f, &
            vir,ener,atomnum,del,apfrc,r_stack,i_stack)
      call rmsdiff(apfrc,cpfrc1(1,j),rms)
      if ( master )then
         write(6,60)atomnum
         do k = 1,3
            diff=apfrc(k)-cpfrc1(k,j)
            write(6,66)k,apfrc(k),cpfrc1(k,j),diff
         end do
         write(6,70)rms
      end if
   end do
   100 continue
   if ( master ) &
         write(6,'(a)')'--------------------------------------------'
   if ( nranatm > 0 )then
      if ( master )write(6,167)
      if ( master ) &
            write(6,'(a)')'--------------------------------------------'
      do j = 1,nranatm
         atomnum = ranatm(j)
         call get_numfrc(xx,ix,ih,ipairs,x,f, &
               vir,ener,atomnum,del,apfrc,r_stack,i_stack)
         call rmsdiff(apfrc,cpfrc2(1,j),rms)
         if ( master )then
            write(6,60)atomnum
            do k = 1,3
               diff = apfrc(k)-cpfrc2(k,j)
               write(6,66)k,apfrc(k),cpfrc2(k,j),diff
            end do
            write(6,70)rms
         end if
      end do
      if ( master ) &
            write(6,'(a)')'--------------------------------------------'
   end if
   if ( master )then
      call mexit(6,0)
   else
      call mexit(0,0)
   end if
   60 format(9x,'NUMERICAL, ANALYTICAL FORCES (diff) from atom ',i6)
   66 format(1x,i6,f14.8,1x,f14.8,1x,f14.8)
   67 format(1x,3(1x,f14.8))
   70 format(1x,'RMS force error = ',e10.3)
   167 format(1x,'Checking numerical force for random atoms')
   168 format(1x,'Checking numerical force for user chosen atoms')
end subroutine debug_frc 
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine crd_mod here]
subroutine crd_mod(crd,atomn,j,del,save)
   implicit none
   _REAL_ crd(3,*),del,save
   integer atomn,j

   save = crd(j,atomn)
   crd(j,atomn) = save + del
   return
end subroutine crd_mod 
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine crd_rst here]
subroutine crd_rst(crd,atomn,j,del,save)
   implicit none
   _REAL_ crd(3,*),del,save
   integer atomn,j

   crd(j,atomn) = save
   return
end subroutine crd_rst 
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_analfrc here]
subroutine get_analfrc(xx,ix,ih,ipairs,x,f, &
      vir,ene,r_stack,i_stack)

   use memory_module

   implicit none
   character(kind=1,len=11) :: routine="get_analfrc"
#  include "md.h"
   _REAL_ xx(*),x(*),f(3,*), &
         ene(*),vir(*),r_stack(*)
   integer ipairs(*),ix(*),i_stack(*)
   character(len=4) ih(*)
   _REAL_ ekcmt(4)
   integer ltmp
   logical, parameter :: do_list_update=.false.

   call fix_xr(x,natom,nspm,ix(i70),xx(l75), &
         ekcmt,xx(l45),xx(lvel),xx(lmass))
   call force(xx,ix,ih,ipairs,x,f,ene(23),vir,r_stack,i_stack, &
         xx(l96), xx(l97), xx(l98), do_list_update)
   return
end subroutine get_analfrc 
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fix_xr here]
subroutine fix_xr(x,natom,nspm,nsp,tma,ekcmt,xr,v,amass)
   implicit none
   _REAL_ x(3,*),tma(*),ekcmt(*), &
         xr(3,*),v(*),amass(*)
   integer natom,nspm,nsp(*),i

   do i = 1,natom
      xr(1,i) = x(1,i)
      xr(2,i) = x(2,i)
      xr(3,i) = x(3,i)
   end do
   return
end subroutine fix_xr 
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_numfrc here]
subroutine get_numfrc(xx,ix,ih,ipairs,x,f, &
      vir,ene,atomn,del,apfrc,r_stack,i_stack)
   
   use memory_module

   implicit none
#  include "md.h"

   _REAL_ xx(*),x(*),f(*), &
         ene(*),vir(*),del,apfrc(3),r_stack(*)
   integer ipairs(*),ix(*),atomn,i_stack(*)
   character(len=4) ih(*)
   _REAL_  enep,enem,save,delta
   integer j
   logical, parameter :: do_list_update=.false.

   do j = 1,3
      delta = del
      call crd_mod(x,atomn,j,delta,save)
      call force(xx,ix,ih,ipairs,x,f,ene(23),vir,r_stack,i_stack, &
            xx(l96), xx(l97), xx(l98), do_list_update)
      enep = ene(23)
      call crd_rst(x,atomn,j,delta,save)
      delta = -del
      call crd_mod(x,atomn,j,delta,save)
      call force(xx,ix,ih,ipairs,x,f,ene(23),vir,r_stack,i_stack, &
            xx(l96), xx(l97), xx(l98), do_list_update)
      enem = ene(23)
      call crd_rst(x,atomn,j,delta,save)
      ! force is negative of gradient
      apfrc(j) = -(enep - enem)/(2.d0*del)
   end do
   return
end subroutine get_numfrc 
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rmsdiff here]
subroutine rmsdiff(apfrc,frc,rms)
   implicit none
   _REAL_ apfrc(3),frc(3)
   _REAL_ num,den,small,rms
   small = 1.d-6
   den = frc(1)**2 + frc(2)**2 + frc(3)**2 + small
   num = (apfrc(1)-frc(1))**2 + (apfrc(2)-frc(2))**2 + &
         (apfrc(3)-frc(3))**2
   rms = sqrt(num/den)
   return
end subroutine rmsdiff 
!----------------------------------------------------
