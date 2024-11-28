! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Amber force field force driver/interface routine 
subroutine force(xx,ix,ih,ipairs,x,f,ener,vir,r_stack,i_stack, &
      fs, rborn, reff ,do_list_update )

   use genborn,           only : egb
!   use poisson_boltzmann, only : pb_force
   use dispersion_cavity, only : npopt, np_force
   use pbtimer_module

   implicit none
   logical do_list_update
   integer   ipairs(*)
   _REAL_ xx(*)
   integer   ix(*)
   character(len=4) ih(*)
   _REAL_ fs(*),rborn(*),reff(*)!,dvdl

#  include "pb_constants.h"
#  include "timer.h"
#  include "md.h"
#  include "box.h"
#  include "memory.h"
#  include "parms.h"
#  include "files.h"
#  include "flocntrl.h"
#  include "pb_md.h"

   logical belly!,nocrst
   _REAL_  r_stack(*)
   integer i_stack(*)

   _REAL_  enmr(3),entr!,devang(4),devdis(4),devtor(4)

   _REAL_  x(*),f(*),ene(30),vir(*)
   _REAL_  ener(*) ! offsets in this ener array = offsets in runmd ener - 22
   save ene

   integer i,m!,nttyp,npair,nhb
!  _REAL_  virvsene,evdw,eelt,e3bod,epol,esurf,edisp
   _REAL_  evdw,eelt,e3bod,epol,esurf,edisp
   _REAL_  epolar,aveper,aveind,avetot!,dipiter,dipole_temp,emtot
!  integer l_r2x,l_rjx,l_tmp1,l_tmp2,l_tmp3,l_tmp4,l_tmp5
!  integer l_tmp6,l_tmp7,l_tmp8,l_jj,l_skipv, l_kvls,l_jvls,l_psi
!  integer l_da,l_sumd
!  integer newbalance
!  save newbalance

   call pbtimer_start(PBTIME_FORCE)
   ene(:) = ZERO

   belly = ibelly > 0
   !nocrst = .false.
   !nttyp = ntypes*(ntypes+1)/2

   ! ZERO OUT THE ENERGIES AND FORCES

   aveper=0.d0
   aveind=0.d0
   avetot=0.d0
   !dipiter=0.d0
   !dvdl=0.d0
   !dipole_temp=0.d0
   do i=1,3
      enmr(i) = 0.d0
   end do
   do i=1,4
      vir(i) = 0.d0
   end do
   !virvsene = 0.d0
   do i=1,3*natom
      f(i) = 0.d0
   end do

   epolar = 0.d0
   e3bod = 0.d0

   ! part I: bonded terms

   call pbtimer_start(PBTIME_BOND)

   ! bonds with H

   if( ntf < 2 ) then
      call bond(ene(6))
   end if

   ! bonds without H

   if( ntf < 3 ) then
      call bond(ene(7))
   end if

   ! angles with H

   if( ntf < 4 ) then
      call angl(ene(8))
   end if

   ! angles without H

   if( ntf < 5 ) then
      call angl(ene(9))
   end if

   ! dihedrals with H

   if( ntf < 6 ) then
      !call ephi(ene(10),ene(11),ene(12))
      call ephi(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48),xx(l15),x,ix,ene(10),ene(11),ene(12))
   end if

   ! dihedrals without H

   if( ntf < 7 ) then
      !call ephi(ene(13),ene(14),ene(15))
      call ephi(nphia,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58),xx(l15),x,ix,ene(13),ene(14),ene(15))
   end if
   call pbtimer_stop(PBTIME_BOND)

   ! part II: restraining terms

   ! positional restraints or targeted MD

   if( natc > 0 ) then
      if (ntr == 1) then
         call xconst(entr)
         ene(20) = entr
      end if
   end if

   ! cap water restraining forces

   if(ifcap > 0) call capwat()

   ! noesy volume penalty energy

   ene(22) = 0.0

   ! part III: implicit solvent nonbonded treatments

   call pbtimer_start(PBTIME_NONBON)

   esurf = 0.d0

   ! gb options

   if( ipb >= 1 ) then

      !call egb(epol,eelt,evdw,esurf)

      ene(2) = 0.0d0
      ene(3) = 0.0d0
      ene(4) = 0.0d0
      ene(23) = 0.0d0

   end if  ! ( ipb >= 1 )
   
   ! pb options
   ! all nonpolar interactions are done in np_force:

   if( ipb >= 1 ) then

      call pbtimer_start(PBTIME_PBFORCE)
      !call pb_force(natom,nres,ntypes,npdec,ix(i02),ix(i04),ix(i06),ix(i10), &
      !              cn1,cn2,xx(l15),x,f,evdw,eelt,epol)
 
      call egb(natom,nres,ntypes,npdec,ix(i02),ix(i04),ix(i06),ix(i10), &
                    cn1,cn2,xx(l15),x,f,evdw,eelt,epol, nbonh, nbonh+nbona, &
                    ix(iibh), ix(ijbh), ih(m06)  )
       

      if ( pbgrid ) pbgrid = .false.
      if ( pbinit ) pbinit = .false.
      ene(2) = evdw
      ene(3) = eelt
      ene(4) = epol
      call pbtimer_stop(PBTIME_PBFORCE)
   end if  ! ( ipb >= 1 )

   if ( inp /= 0 ) then
      call pbtimer_start(PBTIME_NPFORCE)
      esurf = 0.0d0; edisp = 0.0d0
      call np_force(natom,nres,ntypes,ix(i02),ix(i04),ix(i06),&
                    cn1,cn2,x,f,esurf,edisp)
      if ( pbprint ) pbprint = .false.
      ene(23) = esurf
      ene(24) = edisp
      call pbtimer_stop(PBTIME_NPFORCE)
   end if

   call pbtimer_stop(PBTIME_NONBON)

   ! part IV: summary of energy components for printing
   !
   !    ene(1):    total energy
   !    ene(2):    van der Waals
   !    ene(3):    electrostatic energy
   !    ene(4):    10-12 (hb) energy, or GB/PB energy when igb.gt.0
   !    ene(5):    bond energy
   !    ene(6):    angle energy
   !    ene(7):    torsion angle energy
   !    ene(8):    1-4 nonbonds
   !    ene(9):    1-4 electrostatics
   !    ene(10):   constraint energy
   !    ene(11-19):  used a scratch, but not needed further below
   !    ene(20):   position constraint energy
   !    ene(21):   charging free energy result
   !    ene(22):   noe volume penalty
   !    ene(23):   surface-area dependent solvation energy or cavity energy
   !    ene(24):   surface-area dependent dispersion energy

   do m = 2,15
      ene(1) = ene(1) + ene(m)
   end do
   ene(1) = ene(1) + epolar + e3bod + ene(23) + ene(24)

   ene(5) = ene(6)+ene(7)
   ene(6) = ene(8)+ene(9)
   ene(7) = ene(10)+ene(13)
   ene(8) = ene(11)+ene(14)
   ene(9) = ene(12)+ene(15)
   ene(10) = ene(17)+ene(20)+enmr(1)+enmr(2)+enmr(3)
   ene(1) = ene(1)+ene(10)

   ! transfer the energy array to external usage for printing

   ener(1:10) = ene(1:10)
   ener(11) = epolar
   ener(12) = aveper
   ener(13) = aveind
   ener(14) = avetot
   ener(15) = ene(23)
   ener(16) = e3bod
   ener(17) = ene(21)
   ener(18) = ene(24)

   ! part VI: if belly is on, set belly force to zero

   if (belly) call bellyf()

   call pbtimer_stop(PBTIME_FORCE)


end subroutine force
