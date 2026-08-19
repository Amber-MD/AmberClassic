!<compile=optimized>

#include "nfe-utils.h"
#include "nfe-config.h"

!
! barf [ba:rf]  2. "He suggested using FORTRAN, and everybody barfed."
!    - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition)
!

module nfe_sander_hooks

use nfe_constants, only : SL => STRING_LENGTH

implicit none

private

!-----------------------------------------------------------------------------
!                          T H E    H O O K S
!-----------------------------------------------------------------------------

! remd.f
! old REMD subroutines, not used now --------
#ifdef MPI
public :: on_delta
public :: on_exchange
#endif /* MPI */
! -------------------------------------------

! multisander.f
public :: on_multisander_exit

!-----------------------------------------------------------------------------

! sander.f
public :: on_sander_init
public :: on_sander_exit

! force.f
public :: on_force

! mdwrit.f
public :: on_mdwrit

! mdread.f
! old REMD subroutines, not used now
#ifdef MPI
public :: on_mdread1
#endif /* MPI */

! runmd.f
#ifdef NFE_ENABLE_BBMD
public :: on_mdstep
#endif /* NFE_ENABLE_BBMD */

!-----------------------------------------------------------------------------
! old REMD subroutines and parms, not used now 
#ifdef MPI
character(SL), private, save :: initial_mdin_name
NFE_REAL, private, save :: initial_mdin_temp0
NFE_REAL, private, save :: initial_mdin_tempsg,initial_mdin_sgft,initial_mdin_sgff,initial_mdin_psgldg

private :: rem_preinit
#endif /* MPI */

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

#ifdef MPI
subroutine on_delta(o_masterrank, need_U_xx, U_mm, U_mo, U_om, U_oo)

   use nfe_constants, only : ZERO

   use nfe_pmd_hooks, only : pmd_on_delta => on_delta
   use nfe_abmd_hooks, only : abmd_on_delta => on_delta

   implicit none

   integer, intent(in) :: o_masterrank
   logical, intent(in) :: need_U_xx

   NFE_REAL, intent(out) :: U_mm, U_mo, U_om, U_oo

   U_mm = ZERO
   U_mo = ZERO
   U_om = ZERO
   U_oo = ZERO

   call pmd_on_delta(o_masterrank, need_U_xx, U_mm, U_mo, U_om, U_oo)
   call abmd_on_delta(o_masterrank, need_U_xx, U_mm, U_mo, U_om, U_oo)

end subroutine on_delta

!-----------------------------------------------------------------------------

subroutine on_exchange(o_masterrank)

   use nfe_pmd_hooks, only : pmd_on_exchange => on_exchange
   use nfe_abmd_hooks, only : abmd_on_exchange => on_exchange

   implicit none

   integer, intent(in) :: o_masterrank

   call pmd_on_exchange(o_masterrank)
   call abmd_on_exchange(o_masterrank)

end subroutine on_exchange
#endif /* MPI */

!-----------------------------------------------------------------------------

subroutine on_multisander_exit()

   use nfe_sander_proxy, only : proxy_finalize
   use nfe_pmd_hooks, only : pmd_on_multisander_exit => on_multisander_exit
   use nfe_abmd_hooks, only : abmd_on_multisander_exit => on_multisander_exit
   use nfe_stsm_hooks, only : stsm_on_multisander_exit => on_multisander_exit

   implicit none

#  include "nfe-mpi.h"

   call pmd_on_multisander_exit()
   call abmd_on_multisander_exit()
   call stsm_on_multisander_exit()

   NFE_MASTER_ONLY_BEGIN
   call proxy_finalize()
   NFE_MASTER_ONLY_END

end subroutine on_multisander_exit

!-----------------------------------------------------------------------------

!
! 'on_sander_init()' is called at the point when
! MPI is initialized and MDIN/PRMTOP/INPCRD are loaded
!

subroutine on_sander_init(ih, amass, acrds, rem)

   use nfe_constants
   use nfe_sander_proxy

   use nfe_smd_hooks, only : smd_on_sander_init => on_sander_init
   use nfe_pmd_hooks, only : pmd_on_sander_init => on_sander_init
   use nfe_abmd_hooks, only : abmd_on_sander_init => on_sander_init
   use nfe_stsm_hooks, only : stsm_on_sander_init => on_sander_init
#ifdef NFE_ENABLE_BBMD
   use nfe_bbmd_hooks, only : bbmd_on_sander_init => on_sander_init
#endif /* NFE_ENABLE_BBMD */

   implicit none

   character(len = 4), intent(in) :: ih(*)

   NFE_REAL, intent(in) :: amass(*)
   NFE_REAL, intent(in) :: acrds(*)

   integer, intent(in) :: rem
   integer, save :: my_unit = 5

#ifdef MPI
#  include "nfe-mpi.h"
   logical, save :: first_call = .true.

   if (first_call) then
#endif /* MPI */
      call remember_rem(rem)

      NFE_MASTER_ONLY_BEGIN
      call remember_atom_names(ih)
#ifdef MPI
      NFE_MASTER_ONLY_END

      first_call = .false.
   endif ! first_call
#endif /* MPI */

#ifdef MPI
   if (multisander_rem().eq.0) &
      call smd_on_sander_init(my_unit, amass, acrds)
#else
   call smd_on_sander_init(my_unit, amass, acrds)
#endif /* MPI */

   call pmd_on_sander_init(my_unit, amass)
   call abmd_on_sander_init(my_unit, amass)

#ifdef NFE_ENABLE_BBMD
   if (multisander_rem().eq.0.and.multisander_numgroup().gt.1) &
      call bbmd_on_sander_init(my_unit, amass)
#endif /* NFE_ENABLE_BBMD */

#ifdef MPI
   if (multisander_rem().eq.0) &
      call stsm_on_sander_init(my_unit, amass)
#else
   call stsm_on_sander_init(my_unit, amass)
#endif /* MPI */

end subroutine on_sander_init

!-----------------------------------------------------------------------------

subroutine on_sander_exit()

   use nfe_smd_hooks, only : smd_on_sander_exit => on_sander_exit
   use nfe_pmd_hooks, only : pmd_on_sander_exit => on_sander_exit
   use nfe_abmd_hooks, only : abmd_on_sander_exit => on_sander_exit
   use nfe_stsm_hooks, only : stsm_on_sander_exit => on_sander_exit
#ifdef NFE_ENABLE_BBMD
   use nfe_bbmd_hooks, only : bbmd_on_sander_exit => on_sander_exit
#endif /* NFE_ENABLE_BBMD */

   implicit none

   call smd_on_sander_exit()
   call pmd_on_sander_exit()
   call abmd_on_sander_exit()
   call stsm_on_sander_exit()

#ifdef NFE_ENABLE_BBMD
   call bbmd_on_sander_exit()
#endif /* NFE_ENABLE_BBMD */

end subroutine on_sander_exit

!-----------------------------------------------------------------------------

subroutine on_force(x, f, pot)

! Modified by M Moradi
   use nfe_constants, only : ZERO
! Moradi end
   use nfe_sander_proxy

   use nfe_smd_hooks, only : smd_on_force => on_force
   use nfe_pmd_hooks, only : pmd_on_force => on_force
   use nfe_abmd_hooks, only : abmd_on_force => on_force
   use nfe_stsm_hooks, only : stsm_on_force => on_force
#ifdef NFE_ENABLE_BBMD
   use nfe_bbmd_hooks, only : bbmd_on_force => on_force
#endif /* NFE_ENABLE_BBMD */

   implicit none

   NFE_REAL, intent(in) :: x(*)

   NFE_REAL, intent(inout) :: f(*)
   NFE_REAL, intent(inout) :: pot

! Modified by M Moradi
! for driven ABMD
   NFE_REAL :: wdriven = ZERO
   NFE_REAL :: udriven = ZERO

   nfe_pot_ene = null_nfe_pot_ene_rec

   call smd_on_force(x, f, wdriven, udriven, nfe_pot_ene%smd)
   call pmd_on_force(x, f, nfe_pot_ene%pmd)
   call abmd_on_force(x, f, wdriven, udriven, nfe_pot_ene%abmd)
   call stsm_on_force(x, f, nfe_pot_ene%stsm)
#ifdef NFE_ENABLE_BBMD
   call bbmd_on_force(x, f, wdriven, udriven, nfe_pot_ene%bbmd)
#endif /* NFE_ENABLE_BBMD */
! Moradi end
   
   nfe_pot_ene%total = nfe_pot_ene%smd + nfe_pot_ene%pmd + &
                       nfe_pot_ene%abmd + nfe_pot_ene%bbmd + &
                       nfe_pot_ene%stsm
   pot = nfe_pot_ene%total
 
   nfe_real_mdstep = .false. ! reset nfe_real_mdstep to false

end subroutine on_force

!-----------------------------------------------------------------------------

subroutine on_mdwrit()

   use nfe_abmd_hooks, only : abmd_on_mdwrit => on_mdwrit
#ifdef NFE_ENABLE_BBMD
   use nfe_bbmd_hooks, only : bbmd_on_mdwrit => on_mdwrit
#endif /* NFE_ENABLE_BBMD */

   implicit none

   call abmd_on_mdwrit()
#ifdef NFE_ENABLE_BBMD
   call bbmd_on_mdwrit()
#endif /* NFE_ENABLE_BBMD */

end subroutine on_mdwrit

!-----------------------------------------------------------------------------

! Modified by M Moradi
! for selection algorithm (runmd.F90 was modified accordingly as well)
#ifdef NFE_ENABLE_BBMD
subroutine on_mdstep(eptot, x, v, ekmh)

   use nfe_bbmd_hooks, only : bbmd_on_mdstep => on_mdstep
   use nfe_abmd_hooks, only : abmd_on_mdstep => on_mdstep

   implicit none

   NFE_REAL, intent(in) :: eptot
   NFE_REAL, intent(inout) :: x(*)
   NFE_REAL, intent(inout) :: v(*)
   NFE_REAL, intent(inout) :: ekmh ! self-explaining, scientific variable

   call bbmd_on_mdstep(eptot, v, ekmh)
   call abmd_on_mdstep(x, v)

end subroutine on_mdstep
#endif /* NFE_ENABLE_BBMD */
! Moradi end

#ifdef MPI

! ><><><><><><><><><><><><><><><>< R E M D ><><><><><><><><><><><><><><><><><><

!
! remembers MDIN/temp0 [needed for proper REMD restart]
!

subroutine on_mdread1()

   use nfe_sander_proxy

   implicit none

   initial_mdin_name  = sander_mdin_name()
   initial_mdin_temp0 = sander_temp0()
   initial_mdin_tempsg = sander_tempsg()
   initial_mdin_sgft = sander_sgft()
   initial_mdin_sgff = sander_sgff()
   initial_mdin_psgldg = sander_psgldg()

end subroutine on_mdread1

!
! gets proper MDIN filename
!

subroutine rem_preinit(my_mdin, my_idx)

   use nfe_utils
   use nfe_constants
   use nfe_sander_proxy
   use sgld, only: isgld,sorttempsg,tempsglookup

   implicit none

   character(len = SL), intent(out) :: my_mdin
   integer, intent(out) :: my_idx

#  include "nfe-mpi.h"

   NFE_REAL, allocatable :: all_initial_temp0(:)
   NFE_REAL, allocatable :: all_current_temp0(:)

   NFE_REAL, allocatable :: all_initial_tempsg(:)
   NFE_REAL, allocatable :: all_current_tempsg(:)
   NFE_REAL, allocatable :: all_initial_sgft(:)
   NFE_REAL, allocatable :: all_current_sgft(:)
   NFE_REAL, allocatable :: all_initial_sgff(:)
   NFE_REAL, allocatable :: all_current_sgff(:)
   NFE_REAL, allocatable :: all_initial_psgldg(:)
   NFE_REAL, allocatable :: all_current_psgldg(:)

   NFE_REAL, parameter :: TINY = 0.00010000000000000000D0 ! NFE_TO_REAL(0.0001)

   integer :: error, i, j, src_rank, dst_rank

   nfe_assert(multisander_rem().ne.0)
   nfe_assert(multisander_numgroup().gt.1)
   nfe_assert(commmaster.ne.mpi_comm_null)

   ! get all the temperatures

   allocate(all_initial_temp0(multisander_numgroup()), &
            all_current_temp0(multisander_numgroup()), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   call mpi_allgather(initial_mdin_temp0, 1, MPI_DOUBLE_PRECISION, &
                      all_initial_temp0, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
   nfe_assert(error.eq.0)

   call mpi_allgather(sander_temp0(), 1, MPI_DOUBLE_PRECISION, &
                      all_current_temp0, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
   nfe_assert(error.eq.0)

   if (isgld > 0) then
      allocate(all_initial_tempsg(multisander_numgroup()), &
            all_current_tempsg(multisander_numgroup()), &
            all_initial_sgft(multisander_numgroup()), &
            all_current_sgft(multisander_numgroup()), &
            all_initial_sgff(multisander_numgroup()), &
            all_current_sgff(multisander_numgroup()), &
            all_initial_psgldg(multisander_numgroup()), &
            all_current_psgldg(multisander_numgroup()), stat = error)
      if (error.ne.0) NFE_OUT_OF_MEMORY

      call mpi_allgather(initial_mdin_tempsg, 1, MPI_DOUBLE_PRECISION, &
                      all_initial_tempsg, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_allgather(sander_tempsg(), 1, MPI_DOUBLE_PRECISION, &
                      all_current_tempsg, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      nfe_assert(error.eq.0)
      call mpi_allgather(initial_mdin_sgft, 1, MPI_DOUBLE_PRECISION, &
                      all_initial_sgft, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_allgather(sander_sgft(), 1, MPI_DOUBLE_PRECISION, &
                      all_current_sgft, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      nfe_assert(error.eq.0)
      call mpi_allgather(initial_mdin_sgff, 1, MPI_DOUBLE_PRECISION, &
                      all_initial_sgff, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_allgather(sander_sgff(), 1, MPI_DOUBLE_PRECISION, &
                      all_current_sgff, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      nfe_assert(error.eq.0)
      call mpi_allgather(initial_mdin_psgldg, 1, MPI_DOUBLE_PRECISION, &
                      all_initial_psgldg, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      nfe_assert(error.eq.0)

      call mpi_allgather(sander_psgldg(), 1, MPI_DOUBLE_PRECISION, &
                      all_current_psgldg, 1, MPI_DOUBLE_PRECISION, &
                      commmaster, error)
      nfe_assert(error.eq.0)
   endif
   ! src_rank -- rank of the master that has MDIN <--> sander_temp0()
   ! dst_rank -- rank of the master that needs MDIN <--> mdin_temp0

   src_rank = mpi_proc_null
   dst_rank = mpi_proc_null

   do i = 1, multisander_numgroup()
      if (isgld > 0) then
         do j = i + 1, multisander_numgroup()
            if (abs(all_initial_tempsg(i) - all_initial_tempsg(j)).lt.TINY) &
               call fatal('same tempsg in different replicas')
         end do
         if (abs(all_initial_tempsg(i) - sander_tempsg()).lt.TINY) &
            src_rank = i - 1
         if (abs(all_current_tempsg(i) - initial_mdin_temp0).lt.TINY) &
            dst_rank = i - 1
      else
         do j = i + 1, multisander_numgroup()
            if (abs(all_initial_temp0(i) - all_initial_temp0(j)).lt.TINY) &
               call fatal('same temp0 in different replicas')
         end do
         if (abs(all_initial_temp0(i) - sander_temp0()).lt.TINY) &
            src_rank = i - 1
         if (abs(all_current_temp0(i) - initial_mdin_temp0).lt.TINY) &
            dst_rank = i - 1
      end if
   end do

   if (src_rank.eq.mpi_proc_null) then
      write (unit = ERR_UNIT, fmt = '(/a,a,f8.3/)') NFE_ERROR, &
         'could not find MDIN with temp0 = ', sander_temp0()
      call terminate()
   end if

   if (dst_rank.eq.mpi_proc_null) then
      write (unit = ERR_UNIT, fmt = '(/a,a,f8.3/)') NFE_ERROR, &
         'could not find replica that needs MDIN with temp0 = ', &
         initial_mdin_temp0
      call terminate()
   end if

   call mpi_sendrecv(initial_mdin_name, SL, MPI_CHARACTER, dst_rank, 3, &
                     my_mdin, SL, MPI_CHARACTER, src_rank, 3, &
                     commmaster, MPI_STATUS_IGNORE, error)
   nfe_assert(error.eq.0)

   if (isgld > 0) then
      ! Sort temperatures
      call sorttempsg(multisander_numgroup(),all_current_temp0,all_current_tempsg,all_current_sgft,all_current_sgff,all_current_psgldg)
      ! Determine this replca's ID
      my_idx=tempsglookup(multisander_numgroup(),sander_tempsg(), &
                   all_current_tempsg)
      deallocate(all_current_tempsg, all_initial_tempsg,all_current_sgft, all_initial_sgft, &
                     all_current_sgff, all_initial_sgff,all_current_psgldg, all_initial_psgldg)
   else
   ! (bubble) sort the temperatures
   do i = 1, multisander_numgroup()
      do j = i + 1, multisander_numgroup()
         if (all_current_temp0(j).lt.all_current_temp0(i)) &
            call swap(all_current_temp0(i), all_current_temp0(j))
      end do
   end do

   my_idx = 0 ! my index in sorted T-table
   do i = 1, multisander_numgroup()
      if (abs(sander_temp0() - all_current_temp0(i)).lt.TINY) then
         my_idx = i
         exit
      end if
   end do
   endif
   
   nfe_assert(my_idx.gt.0)

   deallocate(all_current_temp0, all_initial_temp0)

end subroutine rem_preinit

! ><><><><><><><><><><><><><><><><><><><><+><><><><><><><><><><><><><><><><><><

#endif /* MPI */

end module nfe_sander_hooks
