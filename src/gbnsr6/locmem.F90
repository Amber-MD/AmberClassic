#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Amber memory interface routine
subroutine locmem
 


  use memory_module

   implicit none

   ! included variables

 
#  include "box.h"
!#  include "memory.h"
#  include "md.h"
#  include "dynph.h"

   ! local variables

!  integer ntbond,ntangl,ntdih,m7,istartr, &
!        istarti, iendr,iendi, istomp,ida_max
   integer ntbond,ntangl,ntdih,m7
!  integer r_ptr,i_ptr,h_ptr,maxpr
   integer maxpr
!  _REAL_ maxpr_float,natom_float,n2_float
   
   ! assign standard partition lengths
   
   maxdup = 2000
   ntbond = nbonh  + nbona 
   ntangl = ntheth + ntheta 
   ntdih  = nphih  + nphia + 2*maxdup
   m7     = nphih  + maxdup
   
   ! set pointers for real arrays
   
   l15 = 1                      ! CG: PARTIAL CHARGES FOR ATOMS
   lwinv = l15 + natom          ! AMASS: ATOMIC MASSES
   lcrd = lwinv + natom         ! C: COORDINATES
   lforce = lcrd + 3*natom      ! F: FORCE
   lvel = lforce + 3*natom + 40 ! V: VELOCITY for MD
   if (imin == 0) then
      lvel2 = lvel + 3*natom    ! VOLD: OLD VELOCITY for MD
      l45 = lvel2  + 3*natom    ! XR: Coords rel. to COM of each mol
   else
      lvel2 = lvel + 6*(3*natom)
      l45 = lvel2 
   end if
   l50 = l45 + 3*natom  ! CONP
   lcrdr = l50 + 0      ! XC (disabled)
   lmass = lcrdr + ntbond  ! Mass
   l75 = lmass + natom  ! TMA
   l95 = l75 + natom    ! Real scratch in shake (2*ntbond) 
   l96 = l95 + 2*ntbond ! FS/GB
   l97 = l96 + natom    ! RBORN/GB
   l98 = l97 + natom    ! REFF/GB
   lastr = l98 + natom
   
   ! set pointers for hollerith arrays
    
   m02 = 1              ! LBRES
   m04 = m02 + nres + 1 ! IGRAPH
   m06 = m04 + natom    ! ISYMBL
   m08 = m06 + natom    ! ITREE
   m12 = m08 + natom    ! N14
   m16 = m12 + natom    ! IARX 
   lasth = m16 + 2*natom

   ! set points for integer arrays
   
   i01 = 1                   ! Marker of principle atom for imaging
   i02 = i01 + natom         ! IPRES
   i04 = i02 + nres + 1      ! IAC
   i06 = i04 + natom         ! NO
   i08 = i06 + ntypes*ntypes ! IBLO
   i10 = i08 + natom         ! INB

   iibh = i10 + 2*nnb        ! IBH
   
   !     ----- BOND ARRAYS -----
   
   ijbh  = iibh  + ntbond    ! JBH
   iicbh = ijbh  + ntbond    ! ICBH
   iiba  = iibh  + nbonh     ! IBA
   ijba  = ijbh  + nbonh     ! JBA
   iicba = iicbh + nbonh     ! ICBA
   i24   = iicbh  + ntbond   ! ITH
   
   !     ----- ANGLE ARRAYS -----
   
   i26 = i24 + ntangl        ! JTH
   i28 = i26 + ntangl        ! KTH
   i30 = i28 + ntangl        ! ICTH
   i32 = i24 + ntheth        ! ITA
   i34 = i26 + ntheth        ! JTA
   i36 = i28 + ntheth        ! KTA
   i38 = i30 + ntheth        ! ICTA
   i40 = i30 + ntangl        ! IPH
   
   !     ----- DIHEDRAL ARRAYS -----
   
   i42 = i40 + ntdih         ! JPH
   i44 = i42 + ntdih         ! KPH
   i46 = i44 + ntdih         ! LPH
   i48 = i46 + ntdih         ! ICPH
   i50 = i40 + m7            ! IPA
   i52 = i42 + m7            ! JPA
   i54 = i44 + m7            ! KPA
   i56 = i46 + m7            ! LPA
   i58 = i48 + m7            ! ICPA
    
   ! misc integer arrays
    
!  CQ
   ibellygp = i58 + natom         ! IBELLYGP
   i66 = ibellygp + natom         ! IROTAT

   i70 = i66 + natom + 1     ! NSP
   lasti = i70
    
   maxpr = 1
   lastpr = maxpr
    
    
end subroutine locmem
