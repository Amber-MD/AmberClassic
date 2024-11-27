 !<compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

module parms

_REAL_, allocatable, dimension(:) :: rk         ! bond force constant
_REAL_, allocatable, dimension(:) :: req        ! bond equilibrium value
_REAL_, allocatable, dimension(:) :: tk         ! angle force constant
_REAL_, allocatable, dimension(:) :: teq        ! angle equilibrium value
_REAL_, allocatable, dimension(:) :: pk         ! torsion barrier height
_REAL_, allocatable, dimension(:) :: pn         ! torsion periodicity
_REAL_, allocatable, dimension(:) :: phase      ! torsion phase
_REAL_, allocatable, dimension(:) :: cn1        ! LJ A coefficient
_REAL_, allocatable, dimension(:) :: cn2        ! LJ B coefficient
_REAL_, allocatable, dimension(:) :: one_scnb   ! 1 / LJ scaling factor
_REAL_, allocatable, dimension(:) :: one_scee   ! 1 / EEL scaling factor
_REAL_, allocatable, dimension(:) :: solty      ! SOLTY parm section (unused)
_REAL_, allocatable, dimension(:) :: gamc       ! cos(phase)*pk
_REAL_, allocatable, dimension(:) :: gams       ! sin(phase)*pk
_REAL_, allocatable, dimension(:) :: asol       ! 10-12 A coefficient
_REAL_, allocatable, dimension(:) :: bsol       ! 10-12 B coefficient
_REAL_, allocatable, dimension(:) :: hbcut      ! Hydrogen bond cutoff

_REAL_, allocatable, dimension(:) :: cn3        !  v
_REAL_, allocatable, dimension(:) :: cn4        ! mjhsieh: for another vdwmodel
_REAL_, allocatable, dimension(:) :: cn5        !  ^

integer, allocatable, dimension(:) :: ipn       ! integer torsion periodicities

! NPHB is the number of h-bond parameters. NIMPRP is the number of
! improper torsional parameters (NPTRA-NIMPRP is the number of regular
! torsional parameters).

! nttyp = ntypes*(ntypes+1)/2  (number of LJ type pairs)

integer :: numbnd
integer :: numang
integer :: nptra
integer :: nphb
integer :: nimprp
integer :: nttyp
integer :: orig_numbnd ! necessary for QM/MM simulations with QM SHAKE

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Allocate space for all of the parm arrays
subroutine allocate_parms

   use memory_module, only : natyp

   implicit none

   integer :: ierror

   ! Do my memory allocation steps
   allocate(rk(numbnd), req(numbnd),             & ! bond params
            tk(numang), teq(numang),             & ! angle params
            pk(nptra), pn(nptra), phase(nptra),  & ! dihedral params
            one_scnb(nptra), one_scee(nptra),    & ! 1-4 scaling factors
            gams(nptra), gamc(nptra),            & ! dihedral 'work' arrays
            cn1(nttyp), cn2(nttyp), cn3(nttyp),  & ! vdw arrays
            cn4(nttyp), cn5(nttyp),              & ! more vdw arrays
            solty(natyp),                        & ! unused ??
            asol(nphb), bsol(nphb), hbcut(nphb), & ! h-bond arrays
            stat=ierror)

   if (ierror /= 0) then
      write(6,'(a)') 'ERROR allocating parameter arrays!'
      call mexit(6, 1)
   end if

   allocate(ipn(nptra), stat=ierror)

   if (ierror /= 0) then
      write(6,'(a)') 'ERROR allocating ipn array!'
      call mexit(6, 1)
   end if

   return

end subroutine allocate_parms

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Allocate space for all of the parm arrays
subroutine add_qmmm_bonds(new_rk, new_req)

   implicit none
   
   ! Formal arguments
   _REAL_, dimension(*), intent(in) :: new_rk      ! New force constants
   _REAL_, dimension(*), intent(in) :: new_req     ! New equilibrium values

   ! Local arguments
   _REAL_, dimension(orig_numbnd) :: rk_holder
   _REAL_, dimension(orig_numbnd) :: req_holder

   integer :: ierror, i

   ! We only have something to do if numbnd > orig_numbnd
   if (numbnd <= orig_numbnd) return

   ! back up rk and req
   rk_holder(1:orig_numbnd) = rk(1:orig_numbnd)
   req_holder(1:orig_numbnd) = req(1:orig_numbnd)

   ! Deallocate our rk/req
   deallocate(rk, req)

   ! Reallocate with the larger numbnd
   allocate(rk(numbnd), req(numbnd), stat=ierror)
   if (ierror /= 0) then
      write(6, '(a)') 'ERROR extending rk/req in add_qmmm_bonds (parms.F90)'
      call mexit(6, 1)
   end if

   ! Restore our original rk, req and then extend it with our new values
   rk(1:orig_numbnd) = rk_holder(1:orig_numbnd)
   req(1:orig_numbnd) = req_holder(1:orig_numbnd)

   ! Now add in our additional rk and req values
   do i = 1, numbnd - orig_numbnd
      rk(orig_numbnd + i) = new_rk(i)
      req(orig_numbnd + i) = new_req(i)
   end do

   ! Now update orig_numbnd (allows us to use add_qmmm_bonds multiple times)
   orig_numbnd = numbnd

   return

end subroutine add_qmmm_bonds

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ deallocates the various data structures
subroutine clean_parms

   implicit none

   if (allocated(rk))      deallocate(rk)
   if (allocated(req))     deallocate(req)
   if (allocated(tk))      deallocate(tk)
   if (allocated(teq))     deallocate(teq)
   if (allocated(pk))      deallocate(pk)
   if (allocated(pn))      deallocate(pn)
   if (allocated(phase))   deallocate(phase)
   if (allocated(cn1))     deallocate(cn1)
   if (allocated(cn2))     deallocate(cn2)
   if (allocated(cn3))     deallocate(cn3)
   if (allocated(cn4))     deallocate(cn4)
   if (allocated(cn5))     deallocate(cn5)
   if (allocated(solty))   deallocate(solty)
   if (allocated(one_scnb))deallocate(one_scnb)
   if (allocated(one_scee))deallocate(one_scee)
   if (allocated(gamc))    deallocate(gamc)
   if (allocated(gams))    deallocate(gams)
   if (allocated(asol))    deallocate(asol)
   if (allocated(bsol))    deallocate(bsol)
   if (allocated(hbcut))   deallocate(hbcut)
   if (allocated(ipn))     deallocate(ipn)

end subroutine clean_parms


end module parms
