 !<compile=optimized>
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
_REAL_, allocatable, dimension(:) :: cn6        ! 12-6-4 LJ C coefficient
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
            cn4(nttyp), cn5(nttyp), cn6(nttyp),  & ! more vdw arrays
            solty(natyp),                        & ! unused ??
            asol(nphb), bsol(nphb), hbcut(nphb), & ! h-bond arrays
            stat=ierror)

   if (ierror /= 0) then
      write(6,*) 'ERROR allocating parameter arrays!'
      call mexit(6, 1)
   end if

   allocate(ipn(nptra), stat=ierror)

   if (ierror /= 0) then
      write(6,*) 'ERROR allocating ipn array!'
      call mexit(6, 1)
   end if

   return

end subroutine allocate_parms

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
   if (allocated(cn6))     deallocate(cn6)
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

#ifdef MPI
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ broadcast parm arrays
subroutine bcast_parms

   use memory_module, only: natyp
   use mpi

   implicit none

#  include "extra.h"
#  include "parallel.h"

   integer ierror

   ! Broadcast the pointers
   call mpi_bcast(numbnd, 1, MPI_INTEGER, 0, commsander, ierror)
   call mpi_bcast(numang, 1, MPI_INTEGER, 0, commsander, ierror)
   call mpi_bcast(nptra, 1, MPI_INTEGER, 0, commsander, ierror)
   call mpi_bcast(nphb, 1, MPI_INTEGER, 0, commsander, ierror)
   call mpi_bcast(nimprp, 1, MPI_INTEGER, 0, commsander, ierror)
   call mpi_bcast(nttyp, 1, MPI_INTEGER, 0, commsander, ierror)

   ! Our slave nodes have not yet allocated the arrays, so do that here
   if (.not. master) then
      call clean_parms() 
      call allocate_parms()
   end if    

   ! Broadcast the real arrays
   call mpi_bcast(rk, numbnd, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(req, numbnd, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(tk, numang, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(teq, numang, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(pk, nptra, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(pn, nptra, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(phase, nptra, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(one_scnb, nptra, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(one_scee, nptra, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(solty, natyp, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(gamc, nptra, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(gams, nptra, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(asol, nphb, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(bsol, nphb, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(hbcut, nphb, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(cn1, nttyp, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(cn2, nttyp, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(cn3, nttyp, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(cn4, nttyp, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(cn5, nttyp, MPI_DOUBLE_PRECISION, 0, commsander, ierror)
   call mpi_bcast(cn6, nttyp, MPI_DOUBLE_PRECISION, 0, commsander, ierror)

   ! Broadcast the integer array
   call mpi_bcast(ipn, nptra, MPI_INTEGER, 0, commsander, ierror)
   
   return

end subroutine bcast_parms
#endif /* MPI */

end module parms
