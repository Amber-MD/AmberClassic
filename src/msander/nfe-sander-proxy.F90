#include "nfe-utils.h"
#include "nfe-config.h"

!
! barf [ba:rf]  2. "He suggested using FORTRAN, and everybody barfed."
!    - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition)
!

module nfe_sander_proxy

use file_io_dat

implicit none

private

public :: sander_mdout_name
public :: sander_mdin_name
public :: sander_mdin_unit

public :: multisander_rem
public :: multisander_initremd
public :: multisander_numgroup

public :: multisander_numwatkeep

public :: terminate
public :: is_master

public :: proxy_finalize

public :: sander_imin
public :: sander_natoms
public :: sander_mdtime
public :: sander_tempsg
public :: sander_sgft
public :: sander_sgff
public :: sander_psgldg
public :: sander_temp0
public :: sander_timestep
public :: sander_init
public :: sander_nstlim
public :: sander_ntp
public :: sander_ntb
public :: sander_nsolut 

#ifdef MPI
public :: set_sander_temp0
#endif

public :: sander_atom_name
public :: remember_atom_names

character(len = 4), private, pointer, save :: atom_names(:) => null()

public :: flush_UNIT

#ifndef NFE_DISABLE_ASSERT
private :: afailed
#endif /* NFE_DISABLE_ASSERT */

public :: remember_rem
public :: remember_initremd

public :: nfe_prt

integer, public,  save :: infe = 0
integer, private, save :: saved_rem = -3212341
#ifdef MPI
logical, private, save :: saved_initremd = .true.
#else
logical, private, save :: saved_initremd = .false.
#endif /* MPI */

! different from PMEMD, we set the default value of
! nfe_real_mdstep to false since force() is called 
! much more often in sander
logical, public,  save :: nfe_real_mdstep = .false.

! NFE potential energy contribution, added by F. Pan
type nfe_pot_ene_rec
    sequence
    double precision    :: smd
    double precision    :: pmd
    double precision    :: abmd
    double precision    :: bbmd
    double precision    :: stsm
    double precision    :: total  ! restraint energy in NFE module
end type nfe_pot_ene_rec

integer, parameter, public    :: nfe_pot_ene_rec_size = 6

type(nfe_pot_ene_rec), parameter, public      :: null_nfe_pot_ene_rec = &
    nfe_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)

type(nfe_pot_ene_rec), public, save :: nfe_pot_ene = null_nfe_pot_ene_rec

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

subroutine terminate()
   implicit none
   call mexit(6, 1)
end subroutine terminate

!-----------------------------------------------------------------------------

character(len=MAX_FN_LEN) function sander_mdin_name()
   use file_io_dat, only : mdin
   implicit none
   sander_mdin_name = mdin
end function sander_mdin_name

!-----------------------------------------------------------------------------

character(len=MAX_FN_LEN) function sander_mdout_name()
   use file_io_dat, only : mdout
   implicit none
   sander_mdout_name = mdout
end function sander_mdout_name

!-----------------------------------------------------------------------------

pure integer function sander_mdin_unit()
   implicit none
   sander_mdin_unit = 5
end function sander_mdin_unit

!-----------------------------------------------------------------------------

pure integer function multisander_numgroup()
   use file_io_dat, only : numgroup
   implicit none
   multisander_numgroup = numgroup
end function multisander_numgroup

!-----------------------------------------------------------------------------

NFE_PURE_EXCEPT_ASSERT integer function multisander_rem()
   implicit none
   ! nfe_assert(saved_rem.ne.-3212341)
   multisander_rem = saved_rem
end function multisander_rem

!-----------------------------------------------------------------------------

subroutine remember_rem(r)
   implicit none
   integer, intent(in) :: r
   saved_rem = r
end subroutine remember_rem

!-----------------------------------------------------------------------------

subroutine remember_initremd(i)
   implicit none
   logical, intent(in) :: i
   saved_initremd = i
end subroutine remember_initremd

!-----------------------------------------------------------------------------

pure logical function multisander_initremd()
   implicit none
   multisander_initremd = saved_initremd
end function multisander_initremd

!-----------------------------------------------------------------------------

pure integer function multisander_numwatkeep()
   implicit none
#  include "../include/md.h"
   multisander_numwatkeep = numwatkeep
end function multisander_numwatkeep

!-----------------------------------------------------------------------------

!
! for use in nfe_assert() & co [where unneeded indirection is acceptable]
!

logical function is_master()

   implicit none

#ifdef MPI
#  include "nfe-mpi.h"
   nfe_assert(commsander /= mpi_comm_null)
   is_master = (sanderrank == 0)
#else
   is_master = .true.
#endif /* MPI */

end function is_master

!-----------------------------------------------------------------------------

subroutine proxy_finalize()
   implicit none

   if (associated(atom_names)) &
      deallocate(atom_names)
end subroutine proxy_finalize

!-----------------------------------------------------------------------------

pure integer function sander_imin()
   implicit none
#include "../include/md.h"
   sander_imin = imin
end function sander_imin

!-----------------------------------------------------------------------------
pure integer function sander_nsolut()
   implicit none
#include "../include/md.h"
#include "../include/memory.h"
   sander_nsolut = natom-(nres-ibgwat+1)*4
end function sander_nsolut
!-----------------------------------------------------------------------------

pure integer function sander_natoms()
   implicit none
#include "../include/memory.h"
   sander_natoms = natom
end function sander_natoms

!-----------------------------------------------------------------------------

pure NFE_REAL function sander_mdtime()
   implicit none
#include "../include/md.h"
   sander_mdtime = t
end function sander_mdtime

!-----------------------------------------------------------------------------

pure NFE_REAL function sander_tempsg()
   use sgld, only:tempsg
   implicit none
#include "../include/md.h"
   sander_tempsg = tempsg
end function sander_tempsg

!-----------------------------------------------------------------------------

pure NFE_REAL function sander_sgft()
   use sgld, only:sgft
   implicit none
#include "../include/md.h"
   sander_sgft = sgft
end function sander_sgft

!-----------------------------------------------------------------------------

pure NFE_REAL function sander_sgff()
   use sgld, only:sgff
   implicit none
#include "../include/md.h"
   sander_sgff = sgff
end function sander_sgff

!-----------------------------------------------------------------------------

pure NFE_REAL function sander_psgldg()
   use sgld, only:psgldg
   implicit none
#include "../include/md.h"
   sander_psgldg = psgldg
end function sander_psgldg

!-----------------------------------------------------------------------------

pure NFE_REAL function sander_temp0()
   implicit none
#include "../include/md.h"
   sander_temp0 = temp0
end function sander_temp0

!-----------------------------------------------------------------------------

#ifdef MPI
subroutine set_sander_temp0(new_temp0)
   implicit none
   NFE_REAL, intent(in) :: new_temp0
#include "../include/md.h"
   temp0 = new_temp0
end subroutine set_sander_temp0
#endif /* MPI */

!-----------------------------------------------------------------------------

pure NFE_REAL function sander_timestep()
   implicit none
#include "../include/md.h"
   sander_timestep = dt
end function sander_timestep

!-----------------------------------------------------------------------------

pure integer function sander_init()
   implicit none
#include "../include/md.h"
   sander_init = init
end function sander_init

!-----------------------------------------------------------------------------

pure integer function sander_nstlim()
   implicit none
#include "../include/md.h"
   sander_nstlim = nstlim
end function sander_nstlim

!-----------------------------------------------------------------------------

pure integer function sander_ntp()
   implicit none
#include "../include/md.h"
   sander_ntp = ntp
end function sander_ntp

!-----------------------------------------------------------------------------

pure integer function sander_ntb()
   implicit none
#include "box.h"
   sander_ntb = ntb
end function sander_ntb

!-----------------------------------------------------------------------------

character(len = 4) function sander_atom_name(n)

   implicit none

   integer, intent(in) :: n

#  include "../include/memory.h"

   nfe_assert(n > 0)
   nfe_assert(n <= sander_natoms())
   nfe_assert(associated(atom_names))

   sander_atom_name = atom_names(n)

end function sander_atom_name

!-----------------------------------------------------------------------------

subroutine remember_atom_names(ih)

   use nfe_constants

   implicit none

   character(len = 4), intent(in) :: ih(*)

#include "../include/memory.h"

   integer :: n, error

   if (associated(atom_names)) &
      deallocate(atom_names)

   nfe_assert(natom > 0)

   allocate(atom_names(natom), stat = error)
   if (error /= 0) then
      write (unit = ERR_UNIT, fmt = '(a,a)') &
         NFE_ERROR, 'out of memory in remember_atom_names()'
      call terminate()
   end if

   do n = 1, natom
      atom_names(n) = ih(m04 + n - 1)
   end do

end subroutine remember_atom_names

!-----------------------------------------------------------------------------

subroutine flush_UNIT(lun)

   implicit none

   integer, intent(in) :: lun

   call flush(lun)

end subroutine flush_UNIT

!-----------------------------------------------------------------------------

#ifndef NFE_DISABLE_ASSERT
subroutine afailed(filename, lineno)

   use nfe_constants, only : ERR_UNIT

   implicit none

   character(len = *), intent(in) :: filename
   integer,            intent(in) :: lineno

   write (unit = ERR_UNIT, fmt = '(/a,a,a,i3,a/)') &
      NFE_ERROR, filename, ':', lineno, ': nfe_assert() failed'

   call terminate()

end subroutine afailed
#endif /* NFE_DISABLE_ASSERT */

!-----------------------------------------------------------------------------

subroutine nfe_prt(i)
  implicit none
  integer, intent(in) :: i

  if (infe .eq. 0) return

  write(i, 100) nfe_pot_ene%smd, nfe_pot_ene%pmd, nfe_pot_ene%abmd
  write(i, 101) nfe_pot_ene%bbmd, nfe_pot_ene%stsm
  write(i, 102)
  return

100 format(' NFE restraints:    SMD  :',f9.3,4x,'PMD  : ',f9.3,4x, &
          'ABMD : ',f9.3)
101 format(20x,'BBMD :',f9.3,4x,'STSM : ',f9.3)
102 format(79('='))
end subroutine nfe_prt

end module nfe_sander_proxy
