!<compile=optimized>
#include "../include/dprec.fh"
#include "../include/assert.fh"

!*******************************************************************************
!
! Module: barostats_mod
!
! Description: Monte-Carlo barostat
!              
!*******************************************************************************

module barostats

   use random

   implicit none

   private

   _REAL_, save :: dvmax = 0.02d0

   integer, save :: total_mcbar_attempts = 0
   integer, save :: total_mcbar_successes = 0

   integer, save      :: mcbar_attempts = 0
   integer, save      :: mcbar_successes = 0
   integer, parameter :: dvmax_interval = 10

   type(rand_gen_state), save :: mcbar_gen

! Variable explanations
!
! Public variables
!   dvmax           : Size of trial dV move
!   total_mcbar_attempts  : # of trial moves we've attempted with MC barostat
!   total_mcbar_successes : # of trial moves we've accepted with MC barostat
!
!   mcbar_attempts  : # of trial moves; used to adjust size of volume move
!   mcbar_successes : # of successes; used to adjust size of volume move
!   dvmax_interval  : # of exchange attempts to do before deciding whether or
!                     not to change the size of the volume move

   public mcbar_trial, mcbar_setup, mcbar_summary

contains

!*******************************************************************************
!
! Subroutine:  mcbar_setup
!
! Description: Sets up the random number generator for the MC barostat
!              
!*******************************************************************************

#ifdef MPI
! Make sure all threads on the same replica have the exact same random # stream
subroutine mcbar_setup(ig)

   use mpi
   implicit none
#include "parallel.h"
#include "extra.h"
   integer, intent(in) :: ig
   integer :: local_ig, ierr

   local_ig = ig + 10

   call mpi_bcast(local_ig, 1, mpi_integer, 0, commsander, ierr)

   call amrset_gen(mcbar_gen, local_ig)

   return

end subroutine mcbar_setup
#else
! Non-MPI case below
subroutine mcbar_setup(ig)

   implicit none

   integer, intent(in) :: ig

   call amrset_gen(mcbar_gen, ig + 10)

   return

end subroutine mcbar_setup
#endif

!*******************************************************************************
!
! Subroutine:  mcbar_trial
!
! Description: Perform the trial move for the MC barostat
!              
!*******************************************************************************

subroutine mcbar_trial(xx, ix, ih, ipairs, x, xc, f, vir, fs, rborn, reff, &
                onereff, qsetup, do_list_update, nstep, nsp, amass, nmropt)

   use constants, only : AVOGADRO, KB, THIRD, HALF, JPKC, TEN_TO_MINUS25
   use memory_module, only : natom, lcrdr
   use nblist, only : ucell, volume
   use state, only : state_rec

   implicit none

   ! need nmropt, nspm, pres0, and ntp
#  include "../include/md.h"
#  include "extra.h"

! Passed parameters
  
   _REAL_, dimension(*)        :: xx          ! Global real array
   _REAL_, dimension(*)        :: x           ! Position array
   _REAL_, dimension(*)        :: xc          ! Restraint position array
   _REAL_, dimension(*)        :: f           ! Force array
   _REAL_, dimension(*)        :: amass       ! Atomic mass array
   _REAL_, dimension(*)        :: fs          ! GB Screen (unused here)
   _REAL_, dimension(*)        :: rborn       ! GB radii (unused here)
   _REAL_, dimension(*)        :: reff        ! effective GB radii (unused here)
   _REAL_, dimension(*)        :: onereff     ! 1 / eff GB radii (unused here)
   _REAL_, dimension(4)        :: vir         ! Virial array (unused here)

   integer, dimension(*)       :: ix          ! Global integer array
   integer, dimension(*)       :: ipairs      ! Pairlist array
   integer, dimension(*)       :: nsp         ! Molecule count array
   integer                     :: nstep       ! Step number we are on
   integer                     :: nmropt      ! from nmr.h

   character(len=4), dimension(*) :: ih       ! Global Hollerith array

   logical, intent(inout) :: qsetup
   logical, intent(out)   :: do_list_update

! Local parameters

   _REAL_, dimension(3) :: dv
   _REAL_, dimension(3) :: rmu

   _REAL_ :: randval
   _REAL_ :: pv_work
   _REAL_ :: delta_area
   _REAL_ :: orig_vol
   _REAL_ :: expfac
   _REAL_ :: nbeta

   integer :: aniso_dim

   type(state_rec)  :: orig_ener   ! Energy of the current state
   type(state_rec)  :: new_ener    ! Energy of the proposed state
   ! converts dyne.A^2/cm to kcal/mol (2.390057e-27 * 6.022e23)
   _REAL_, parameter :: TENSION_CONV = 0.0014393264316443592

   ! This is another mcbar attempt
   total_mcbar_attempts = total_mcbar_attempts + 1
   mcbar_attempts = mcbar_attempts + 1

   dv(:) = 0.d0
   rmu(:) = 0.d0

   ! Back up the original volume
   orig_vol = volume

   ! Get the original energy
   call force(xx, ix, ih, ipairs, x, f, orig_ener, vir,  &
      fs, rborn, reff, onereff, qsetup, do_list_update, nstep)

   ! Decrement the NMR opt counter since this force call should not 'count'
   if (nmropt .gt. 0) call nmrdcp

   ! Get the dV move we plan on doing
   if (ntp .eq. 1) then
      ! Isotropic
      call amrand_gen(mcbar_gen, randval)
      dv(1) = (randval - HALF) * dvmax
      dv(2) = dv(1)
      dv(3) = dv(1)
   else if (ntp .eq. 2) then
    ! Anisotropic -- pick one dimension to change(if baroscalingdir==0)
      if (baroscalingdir .eq. 0) then
    ! 'baroscalingdir' -- controlling the direction for anisotropic
    ! pressure scaling. Choose one direction(x, y or z), the box only
    ! scales along the chosen direction, its lengths along the other two
    ! directions are fixed. Default value 0, randomly pick x, y or
    ! z-direction to apply box size change. 'baroscalingdir = 1'-->
    ! only x-direction is allowed to scale; '= 2'--> y-direction;
    ! '= 3'--> z-direction. --Yeyue Xiong
         call amrand_gen(mcbar_gen, randval)
         aniso_dim = int(randval * 3.d0 * 0.99999999d0) + 1
         call amrand_gen(mcbar_gen, randval)
         dv(:) = 0.d0
         dv(aniso_dim) = (randval - HALF) * dvmax
      else
         select case (baroscalingdir)
            case(1)
               aniso_dim = 1
               call amrand_gen(mcbar_gen, randval)
               dv(:) = 0.d0
               dv(aniso_dim) = (randval - HALF) * dvmax
            case(2)
               aniso_dim = 2
               call amrand_gen(mcbar_gen, randval)
               dv(:) = 0.d0
               dv(aniso_dim) = (randval - HALF) * dvmax
            case(3)
               aniso_dim = 3
               call amrand_gen(mcbar_gen, randval)
               dv(:) = 0.d0
               dv(aniso_dim) = (randval - HALF) * dvmax
         end select
      end if
   end if

   rmu(1) = (1.d0 + dv(1)) ** THIRD
   rmu(2) = (1.d0 + dv(2)) ** THIRD
   rmu(3) = (1.d0 + dv(3)) ** THIRD

   call scale_system_volume(rmu, natom, x, xc, amass, nsp)

   ! Get the proposed state energy
   call force(xx, ix, ih, ipairs, x, f, new_ener, vir,  &
      fs, rborn, reff, onereff, qsetup, do_list_update, nstep)

   ! Decrement the NMR opt counter since this force call should not 'count'
   if (nmropt .gt. 0) call nmrdcp

   ! p*dV (6.02204d-2 / 4184.0d0 converts from bar*A^3/particle to kcal/mol)
   pv_work = pres0 * (volume - orig_vol) * AVOGADRO * TEN_TO_MINUS25 / JPKC

   nbeta = - 1 / (temp0 * KB)

   expfac = new_ener%pot%tot - orig_ener%pot%tot + pv_work + &
            nspm * log(rmu(1) * rmu(2) * rmu(3)) / nbeta

   call amrand_gen(mcbar_gen, randval)

   ! Monte carlo decision

   ! DEBUG
!  if (master) &
!     write(6, '(a)', advance='NO') '| Attempting MC barostat volume change: '

   if (randval .lt. exp(nbeta * expfac)) then
      ! DEBUG
!     if(master) write(6, '(a)') 'Accepted'
      mcbar_successes = mcbar_successes + 1
      total_mcbar_successes = total_mcbar_successes + 1
   else
      ! Reject -- rescale everything and put back original forces
      rmu(1) = 1.d0 / rmu(1)
      rmu(2) = 1.d0 / rmu(2)
      rmu(3) = 1.d0 / rmu(3)

      call scale_system_volume(rmu, natom, x, xc, amass, nsp)

      ! DEBUG
!     if(master) write(6, '(a)') 'Rejected'
   end if

   if (mod(mcbar_attempts, dvmax_interval) .eq. 0) then

      ! If our success fraction is too large or too small, adjust dvmax

      if (mcbar_successes .ge. 0.75d0 * mcbar_attempts) then

         dvmax = dvmax * 1.1d0
         mcbar_attempts = 0
         mcbar_successes = 0
         if (master) &
            write(6, '(a)') '| MC Barostat: Increasing size of volume moves'

      else if (mcbar_successes .le. 0.25d0 * mcbar_attempts) then

         dvmax = dvmax / 1.1d0
         mcbar_attempts = 0
         mcbar_successes = 0
         if (master) &
            write(6, '(a)') '| MC Barostat: Decreasing size of volume moves'

      end if

   end if

   return
   
end subroutine mcbar_trial

!*******************************************************************************
!
! Subroutine:  scale_system_volume
!
! Description: Scales the system volume
!              
!*******************************************************************************

subroutine scale_system_volume(rmu, natom, crd, rcrd, amass, nsp)
   
   use nblist, only : volume, oldrecip, ucell, fill_tranvec
#ifdef MPI
   use softcore, only : sc_pscale, ifsc
#endif

#ifdef MPI
   use mpi
#endif
   implicit none

#  include "box.h"
#  include "../include/md.h"
#ifdef MPI
#  include "extra.h"
#  include "parallel.h"
#endif

   ! The following variables are used from the headers:
   !     box         Periodic box vectors    (box.h)
   !     nspm        Number of molecules     (md.h) 
   !     nspcal      Do molecular scaling?   (md.h)
   !     master      Is this process 0?      (extra.h)
   !     commsander  Sander comm             (parallel.h)

! Formal arguments
  
   integer, intent(in)               :: natom   ! number of atoms
   integer, dimension(*), intent(in) :: nsp     ! Molecule list

   _REAL_, intent(in), dimension(3)    :: rmu   ! Box scaling factors
   _REAL_, intent(inout), dimension(*) :: crd   ! Atomic coordinates
   _REAL_, intent(inout), dimension(*) :: rcrd  ! Restraint coordinates
   _REAL_, intent(in), dimension(*)    :: amass ! Atomic mass array

! Local variables
   
#ifdef MPI
   integer :: ierr
#endif

   box(1:3) = box(1:3)*rmu(1:3)
   
   !    WARNING!!   This is not correct for non-orthogonal boxes if
   !    NTP > 1 (i.e. non-isotropic scaling).  Currently general cell
   !    updates which allow cell angles to change are not implemented.
   !    The viral tensor computed for ewald is the general Nose Klein,
   !    however the cell response needs a more general treatment.
   
   call redo_ucell(rmu)
   ! keep tranvec up to date, rather than recomputing each MD step.
   call fill_tranvec()  ! tranvec is dependent on only ucell

#ifdef MPI /* SOFT CORE */
   ! if softcore potentials and the dual topology approach are used
   ! C.O.M. scaling has to be changed to account for different masses 
   ! of the same molecule in V0 and V1. This is quite inefficient and is
   ! therefore done in a separate routine in softcore.f
   ! only both masters actually do the computation for ifsc==1
   ! the scaled coordinates are then broadcast to the nodes
   if (icfe /= 0 .and. ifsc == 1) then
      if (master) then
         call sc_pscale(crd,amass,nspm,nsp,oldrecip,ucell)
      end if
      call mpi_bcast(crd,natom*3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   else
#endif
      call ew_pscale(natom,crd,amass,nspm,nsp,npscal)
#ifdef MPI /* SOFT CORE */
   end if
#endif
   if (ntr > 0 .and. nrc > 0) &
      call ew_pscale(natom,rcrd,amass,nspm,nsp,npscal)

   ! Unlike pmemd, the skin check and nonbond list update is all handled inside
   ! the call to force, so we can skip the skincheck here.

   return
  
end subroutine scale_system_volume

!*******************************************************************************
!
! Subroutine: mcbar_summary
!
! Description: Print out a summary of the MC barostat statistics
!
!*******************************************************************************

subroutine mcbar_summary()

   implicit none

#  include "extra.h"
   _REAL_ :: success_percent

   if (.not. master) return

   success_percent = dble(total_mcbar_successes) / &
                     dble(total_mcbar_attempts) * 100.d0

   write(6, '(("| MC Barostat: ",i10," volume changes attempted."))') &
      total_mcbar_attempts
  
   write(6, '(("| MC Barostat: ",i10," changes successful (",1f6.2,"%)"))') &
      total_mcbar_successes, success_percent
   write(6, '(t2,78("-"),/)')

   return

end subroutine mcbar_summary

end module barostats
