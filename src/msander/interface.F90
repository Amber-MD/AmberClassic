! Do not optimize right now -- segfaults with ifort

! This file contains the code necessary to initialize the sander state without
! reading command-line flags or input files, as well as computing energies and
! forces via an API.
#include "../include/assert.fh"

#ifdef LES
module sanderles_api
#else
module sander_api
#endif

   use constants, only: MAX_QUANTUM_ATOMS
   use file_io_dat, only: MAX_FN_LEN
   use prmtop_type, only: prmtop_struct, destroy_prmtop_struct, read_prmtop_file
   use qmmm_module, only: qmmm_input_options
   use state, only: potential_energy_rec

! qmmm_input options is a data structure defined in qmmm_module, shown below:
!
! type qmmm_input_options
!    ! Allow a way to input options programmatically through an API rather than
!    ! requiring an input file
!    sequence
!    integer :: iqmatoms(MAX_QUANTUM_ATOMS), qmgb, lnk_atomic_no, &
!               ndiis_matrices, ndiis_attempts, lnk_method, qmcharge, &
!               corecharge, buffercharge, spin, qmqmdx, verbosity, &
!               printcharges, printdipole, print_eigenvalues, peptide_corr, &
!               itrmax, printbondorders, qmshake, qmmmrij_incore, &
!               qmqm_erep_incore, pseudo_diag, qm_ewald, qm_pme, kmaxqx, &
!               kmaxqy, kmaxqz, ksqmaxq, kappa, qmmm_int, adjust_q, tight_p_conv, &
!               diag_routine, density_predict, fock_predict, vsolv, &
!               dftb_maxiter, dftb_disper, dftb_chg, abfqmmm, hot_spot, &
!               qmmm_switch, core_iqmatoms(MAX_QUANTUM_ATOMS), &
!               buffer_iqmatoms(MAX_QUANTUM_ATOMS)
!    _REAL_ :: qmcut, lnk_dis, scfconv, errconv, dftb_telec, dftb_telec_step, &
!              fockp_d1, fockp_d2, fockp_d3, fockp_d4, damp, vshift, &
!              pseudo_diag_criteria, min_heavy_mass, r_switch_hi, r_switch_lo
!    character(len=8192) :: qmmask, coremask, buffermask, centermask
!    character(len=256) :: dftb_3rd_order, dftb_slko_path
!    character(len=12) :: qm_theory
! end type qmmm_input_options

! potenetial_energy_rec is a data structure defined in state.F90, shown below:

! type potential_energy_rec
!   sequence
!   _REAL_  :: tot
!   _REAL_  :: vdw
!   _REAL_  :: elec
!   _REAL_  :: gb
!   _REAL_  :: bond
!   _REAL_  :: angle
!   _REAL_  :: dihedral
!   _REAL_  :: vdw_14
!   _REAL_  :: elec_14
!   _REAL_  :: constraint
!   _REAL_  :: polar
!   _REAL_  :: hbond
!   _REAL_  :: surf
!   _REAL_  :: scf
!   _REAL_  :: disp
!   _REAL_  :: dvdl
!   _REAL_  :: angle_ub
!   _REAL_  :: imp
!   _REAL_  :: cmap
!   _REAL_  :: emap
!   _REAL_  :: les
!   _REAL_  :: noe
!   _REAL_  :: pb
!   _REAL_  :: rism
!   _REAL_  :: ct
!   _REAL_  :: amd_boost
! end type potential_energy_rec

   ! Datatypes
   type :: sander_input
      sequence
      ! Floats (input parameters)
      double precision :: extdiel
      double precision :: intdiel
      double precision :: rgbmax
      double precision :: saltcon
      double precision :: cut
      double precision :: dielc
      double precision :: rdt
      double precision :: fswitch
      double precision :: restraint_wt
      double precision :: ala
      double precision :: arg
      double precision :: asn
      double precision :: asp
      double precision :: cys
      double precision :: gln
      double precision :: glu
      double precision :: gly
      double precision :: his
      double precision :: hip
      double precision :: ile
      double precision :: leu
      double precision :: lys
      double precision :: met
      double precision :: phe
      double precision :: pro
      double precision :: ser
      double precision :: thr
      double precision :: triptophan
      double precision :: tyr
      double precision :: valine
      double precision :: bb      

      ! Integers (toggle options)
      integer :: igb
      integer :: alpb
      integer :: gbsa
      integer :: lj1264
      integer :: plj1264 !New2021
      integer :: ipb
      integer :: inp
      integer :: vdwmeth
      integer :: ew_type
      integer :: ntb
      integer :: ifqnt
      integer :: jfastw
      integer :: ntf
      integer :: ntc
      integer :: ntr
      integer :: ibelly
      integer :: mask_from_ref ! Antoine Marion

      ! Character string arrays
      character(len=256) :: restraintmask
      character(len=256) :: bellymask
      character(len=256) :: refc
   end type sander_input

   private

   logical, save :: is_setup_ = .false.

   integer, dimension(:), allocatable :: ipairs

   public :: gas_sander_input, pme_sander_input, sander_setup, energy_forces, &
             set_positions, set_box, sander_input, qmmm_input_options, &
             potential_energy_rec, sander_cleanup, MAX_FN_LEN, get_box, &
             qm_sander_input, sander_natom, prmtop_struct, read_inpcrd_file, &
             get_inpcrd_natom, destroy_prmtop_struct, read_prmtop_file, &
#ifdef NO_ALLOCATABLES_IN_TYPE
             is_setup, get_positions
#else
             sander_setup2, is_setup, get_positions
#endif

contains

! Returns .true. if sander is set up and .false. otherwise
logical function is_setup()

   implicit none

   is_setup = is_setup_

   return

end function is_setup

! Sets up input options for a gas-phase (i.e., no periodic boundaries) system
!
! Parameters
! ----------
! inp : type(sander_input)
!     struct of input options that will be filled by this subroutine
! gb : integer (optional)
!     If GB model is desired, igb will be set to this.  If not present or 0, a
!     gas phase calculation is set up. If set to 10, a PB calculation is set up
!
! If gb is given, the given GB model is set up. Otherwise, the system is set up
! for a gas-phase simulation
subroutine gas_sander_input(inp, gb)

   use constants, only: NO_INPUT_VALUE

   implicit none

   type(sander_input), intent(out) :: inp
   integer, optional, intent(in)   :: gb

   inp%ntb = 0
   inp%lj1264 = 0
   inp%plj1264 = 0 !New2021
   inp%alpb = 0
   inp%vdwmeth = 0
   inp%ew_type = 0
   inp%inp = 0
   inp%ipb = 0
   inp%igb = 6
   inp%gbsa = 0
   inp%jfastw = 0
   inp%ifqnt = 0
   inp%extdiel = 1.d0
   inp%intdiel = 1.d0
   inp%rgbmax = 25.d0
   inp%saltcon = 0.d0
   inp%cut = 1000.d0
   inp%dielc = 1.d0
   inp%rdt = 0.d0
   inp%ntc = 1
   inp%ntf = 1
   inp%fswitch = -1.d0
   inp%ntr = 0
   inp%ibelly = 0
   inp%mask_from_ref = 0 ! Antoine Marion
   inp%restraint_wt = 0.d0
   inp%restraintmask = ''
   inp%bellymask = ''
   inp%refc = ''
   inp%ala = 0.d0
   inp%arg = 0.d0
   inp%asn = 0.d0
   inp%asp = 0.d0
   inp%cys = 0.d0
   inp%gln = 0.d0
   inp%glu = 0.d0
   inp%gly = 0.d0
   inp%his = 0.d0
   inp%hip = 0.d0
   inp%ile = 0.d0
   inp%leu = 0.d0
   inp%lys = 0.d0
   inp%met = 0.d0
   inp%phe = 0.d0
   inp%pro = 0.d0
   inp%ser = 0.d0
   inp%thr = 0.d0
   inp%triptophan = 0.d0
   inp%tyr = 0.d0
   inp%valine = 0.d0
   inp%bb = 0.d0

   if (present(gb)) then
      if (gb /= 0 .and. gb /= 1 .and. gb /= 2 .and. gb /= 5 .and. &
          gb /= 6 .and. gb /= 7 .and. gb /= 8 .and. gb /= 10) then
         write(0,*) 'Illegal gb model. Setting to vacuum'
      else
         inp%igb = gb
         inp%extdiel = 78.5
      end if
      ! Vacuum
      if (gb == 0) inp%igb = 6
      ! PB
      if (gb == 10) then
         inp%ipb = 1
         inp%inp = 1
      end if
   end if

end subroutine gas_sander_input

! Sets up input options for a periodic system treated with PME
!
! Parameters
! ----------
! inp : type(sander_input)
!     struct of input options that will be filled by this subroutine
subroutine pme_sander_input(inp)

   use constants, only: NO_INPUT_VALUE

   implicit none

   type(sander_input), intent(out) :: inp

   inp%ntb = 1
   inp%igb = 0
   inp%alpb = 0
   inp%lj1264 = 0
   inp%plj1264 = 0 !New2021
   inp%ipb = 0
   inp%inp = 0
   inp%vdwmeth = 1
   inp%ew_type = 0
   inp%gbsa = 0
   inp%ifqnt = 0
   inp%jfastw = 0
   inp%extdiel = 1.d0
   inp%intdiel = 1.d0
   inp%cut = 8.d0
   inp%rgbmax = 25.d0
   inp%saltcon = 0.d0
   inp%dielc = 1.d0
   inp%rdt = 0.d0
   inp%ntc = 1
   inp%ntf = 1
   inp%fswitch = -1.d0
   inp%ntr = 0
   inp%ibelly = 0
   inp%mask_from_ref = 0 ! Antoine Marion
   inp%restraint_wt = 0.d0
   inp%restraintmask = ''
   inp%bellymask = ''
   inp%refc = ''
   inp%ala = 0.d0
   inp%arg = 0.d0
   inp%asn = 0.d0
   inp%asp = 0.d0
   inp%cys = 0.d0
   inp%gln = 0.d0
   inp%glu = 0.d0
   inp%gly = 0.d0
   inp%his = 0.d0
   inp%hip = 0.d0
   inp%ile = 0.d0
   inp%leu = 0.d0
   inp%lys = 0.d0
   inp%met = 0.d0
   inp%phe = 0.d0
   inp%pro = 0.d0
   inp%ser = 0.d0
   inp%thr = 0.d0
   inp%triptophan = 0.d0
   inp%tyr = 0.d0
   inp%valine = 0.d0
   inp%bb = 0.d0

end subroutine pme_sander_input

! Initializes a struct with QM options to all default values
!
! Parameters
! ----------
! inp : type(qmmm_input_options)
!     struct of QM input options that will be filled by this subroutine
subroutine qm_sander_input(inp)

   use qmmm_module, only: default_qmmm_input_options
   implicit none

   type(qmmm_input_options), intent(out) :: inp

   call default_qmmm_input_options(inp)

end subroutine qm_sander_input

! Initializes the major data structures needed to evaluate energies and forces
!
! Parameters
! ----------
! prmname : string
!     Name of the topology file
! coordinates : double precision(3*natom)
!     Starting coordinates
! inbox : double precision(6)
!     Box dimensions
! input_options : type(sander_input)
!     struct of input options used to set up the calculation
! qmmm_options : type(qmmm_input_options) [optional]
!     struct of input options used to set up the QM/MM part of the calculation
! ierr : integer
!     Set to 0 for success or 1 if the setup failed
subroutine sander_setup(prmname, coordinates, inbox, input_options, qmmm_options, ierr)

#define USE_PRMTOP_FILE 1
#include "interface_setup.F90"
#undef USE_PRMTOP_FILE

end subroutine sander_setup

#ifndef NO_ALLOCATABLES_IN_TYPE
! Initializes the major data structures needed to evaluate energies and forces
!
! Parameters
! ----------
! parmdata : prmtop_struct
!     Struct with all of the prmtop data stored in it
! coordinates : double precision(3*natom)
!     Starting coordinates
! inbox : double precision(6)
!     Box dimensions
! input_options : type(sander_input)
!     struct of input options used to set up the calculation
! qmmm_options : type(qmmm_input_options) [optional]
!     struct of input options used to set up the QM/MM part of the calculation
! ierr : integer
!     Set to 0 for success or 1 if the setup failed
subroutine sander_setup2(parmdata, coordinates, inbox, input_options, qmmm_options, ierr)

#include "interface_setup.F90"

end subroutine sander_setup2
#endif /* NO_ALLOCATABLES_IN_TYPE */

#undef rem

subroutine api_mdread1(input_options, ierr)
#define API 1
#include "mdread1.F90"
#undef API
end subroutine api_mdread1

subroutine api_mdread2(x, ix, ih, ierr)
#define API 1
#include "mdread2.F90"
#undef API
end subroutine api_mdread2

! Sets the atomic positions
subroutine set_positions(positions)

   use memory_module, only: natom, lcrd, x

   implicit none

   double precision, intent(in), dimension(natom*3) :: positions

   integer :: i, start

   if (.not. is_setup_) return ! Not sure how to warn...

   start = lcrd - 1
   do i = 1, natom * 3
      x(start + i) = positions(i)
   end do

end subroutine set_positions

! Sets the periodic box vectors
subroutine set_box(a, b, c, alpha, beta, gamma)

    implicit none

    double precision, intent(in) :: a, b, c, alpha, beta, gamma

    if (.not. is_setup_) return ! Not sure how to warn...

    call fill_ucell(a, b, c, alpha, beta, gamma)

end subroutine set_box

subroutine get_box(a_, b_, c_, alpha_, beta_, gamma_)

    use nblist, only : a, b, c, alpha, beta, gamma

    implicit none

    double precision, intent(out) :: a_, b_, c_, alpha_, beta_, gamma_

    if (.not. is_setup_) return ! not sure how to warn

    a_ = a
    b_ = b
    c_ = c
    alpha_ = alpha
    beta_ = beta
    gamma_ = gamma

end subroutine get_box

! Returns the currently active positions
subroutine get_positions(positions)

    use memory_module, only: natom, lcrd, x

    implicit none

    double precision, intent(out), dimension(natom*3) :: positions
    integer :: start, i

    positions(:) = 0.d0
    if(.not. is_setup_) return

    start = lcrd - 1
    do i = 1, natom * 3
        positions(i) = x(start + i)
    end do

end subroutine get_positions

! Computes the energies and forces with the current positions
subroutine energy_forces(energy, forces)

   use memory_module
   use state

   implicit none

   ! Input parameters
   type(potential_energy_rec), intent(out) :: energy
   double precision, dimension(natom*3), intent(out) :: forces

   ! Private variables
   type(state_rec) :: ener
   logical :: qsetup, do_list_update
   integer :: nstep, i, start

   ! Initialization
   energy = null_potential_energy_rec
   forces(:) = 0.d0
   qsetup = .true.
   do_list_update = .true.
   nstep = 1

   if (.not. is_setup_) return ! Not sure how to warn...

   call force(x,ix,ih,ipairs,x(lcrd),x(lforce),ener,ener%vir, &
              x(l96),x(l97),x(l98),x(l99), qsetup, &
              do_list_update,nstep)

   start = lforce - 1
   do i = 1, natom*3
      forces(i) = x(start+i)
   end do
   energy = ener%pot

   return

end subroutine energy_forces

subroutine sander_cleanup()

   use qmmm_adaptive_module, only: adaptive_reset
   use amoeba_runmd, only: AM_RUNMD_reset
   use charmm_mod, only : charmm_active, charmm_deallocate_arrays
   use decomp, only : deallocate_int_decomp, deallocate_real_decomp
#ifdef LES
   use genbornles, only: deallocate_gb
#else
   use genborn, only: deallocate_gb
#endif
   use memory_module, only: memory_free
   use molecule, only : deallocate_molecule
   use nblist, only: nblist_deallocate
   use parms, only: clean_parms
   use qmmm_module, only: deallocate_qmmm, qmmm_nml, qmmm_struct, qmmm_vsolv, &
                          qm2_params, qmewald
   use stack, only: deallocate_stacks
!   use xray_interface_module, only : xray_fini

   implicit none
#include "../include/md.h"
#include "ew_cntrl.h"

   if (.not. is_setup_) return ! Not sure how to warn...

   call memory_free
   call clean_parms
   if (igb /= 0 .and. igb /= 10 .and. ipb == 0) call deallocate_gb
   if (idecomp == 1 .or. idecomp == 2) then
      call deallocate_int_decomp
      call deallocate_real_decomp
   end if
   call deallocate_stacks
   call nblist_deallocate
!   call xray_fini
   if (allocated(ipairs)) deallocate(ipairs)
   ! periodic does not get reset when ntb==0 inside sander, so do it here
   periodic = 0
   ! QM/MM cleanup
   if (qmmm_nml%ifqnt) then
      call deallocate_qmmm(qmmm_nml, qmmm_struct, qmmm_vsolv, qm2_params)
      call get_qm2_forces_reset
      qmmm_nml%ifqnt = .false.
      ! Initialize all of the first call flags
      qmewald%ewald_startup = .true.
      qmmm_struct%qm_mm_first_call = .true.
      qmmm_struct%fock_first_call = .true.
      qmmm_struct%fock2_2atm_first_call  = .true.
      qmmm_struct%qm2_allocate_e_repul_first_call = .true.
      qmmm_struct%qm2_calc_rij_eqns_first_call = .true.
      qmmm_struct%qm2_scf_first_call = .true.
      qmmm_struct%zero_link_charges_first_call = .true.
      qmmm_struct%adj_mm_link_pair_crd_first_call = .true.
   end if
   ! CHARMM stuff
   if (charmm_active) call charmm_deallocate_arrays()
   charmm_active = .false.
   ! Reset all of the "first"
   call AM_RUNMD_reset
   call mdeng_reset
   call mdwrit_reset
   call minwrit_reset
   call pcshift_reset
   call adaptive_reset
   call shake_reset
   call ewald_self_reset
   call handle_induced_reset
   call nxtsec_reset
   call deallocate_molecule

   is_setup_ = .false.

end subroutine sander_cleanup

! Returns the number of atoms defined in the current system
subroutine sander_natom(n)

   implicit none

#include "../include/memory.h"

   integer, intent(out) :: n

   if (.not. is_setup_) then
      n = 0
   else
      n = natom
   end if

   return

end subroutine sander_natom

! Parses an inpcrd file and fills the coordinate and box arrays with the values
! in the inpcrd file. This routine assumes that the caller already knows how
! many atoms are present and has pre-allocated enough space in the coordinate
! array to hold all of the coordinates. The box variable should have 6 elements.
! If no box is present, then all box lengths and angles will be set to 0. If
! there was any error reading the inpcrd file, ierr will be set to 1.
! Otherwise, it will be 0 on success
subroutine read_inpcrd_file(filename, coordinates, box, ierr)

#ifdef BINTRAJ
   use AmberNetcdf_mod
   use binrestart
#endif
   use file_io_dat, only : INPCRD_UNIT

   implicit none

   ! Passed variables

   character(len=*), intent(in)  :: filename
   double precision, intent(out) :: coordinates(*)
   double precision, intent(out) :: box(6)
   integer, intent(out)          :: ierr

   ! Local variables

   integer           :: ncid
   integer           :: alloc_failed
   integer           :: natom
   integer           :: i
   character(len=80) :: title
   double precision  :: time
   double precision, dimension(1) :: remd_values
   double precision, allocatable  :: velocities(:)

   box(:) = 0.d0
   ierr = 0

#ifdef BINTRAJ
   ! Make sure the NetCDF module doesn't scream at us...
   verbose_netcdf = .false.

   if (NC_checkRestart(filename)) then
      ! See if there are box coordinates
      if (NC_openRead(filename, ncid)) then
         go to 666 ! error label
      end if
      if (NC_readRestartBox(ncid,box(1),box(2),box(3),box(4),box(5),box(6))) &
         box(:) = 0.d0 ! There _is_ no box...
      call NC_close(ncid)
      call read_nc_restart(filename, title, 1, natom, coordinates, velocities, &
                           remd_values, 1, time)
   else
#endif
      ! Try reading the raw file
      open(unit=INPCRD_UNIT, file=filename, status='OLD', form='FORMATTED', &
           iostat=alloc_failed)
      if (alloc_failed /= 0) then
         ierr = 1
         return
      end if

      read(INPCRD_UNIT, '(a80)') title
      read(INPCRD_UNIT, *, err=666) natom
      allocate(velocities(natom*3), stat=alloc_failed)
      if (alloc_failed == 1) then
         coordinates(1:natom*3) = 0.d0
         box(:) = 0.d0
         ierr = 1
         return
      end if
      read(INPCRD_UNIT, '(6f12.7)', end=666, err=666) &
            (coordinates(i), i=1,natom*3)
      ! Now get the next 6 numbers. If these are the last numbers in the file,
      ! then we know it is really the set of box dimensions. If there are more
      ! numbers but not enough, it's a corrupt file, otherwise, read all of the
      ! velocities
      read(INPCRD_UNIT, '(6f12.7)', end=667, err=666) &
            (velocities(i), i=1,6)
      read(INPCRD_UNIT, '(6f12.7)', end=668, err=666) &
            (velocities(i), i=7,12)
      read(INPCRD_UNIT, '(6f12.7)', end=666, err=666) &
            (velocities(i), i=13,natom*3)
      ! Now see if we can read our box
      read(INPCRD_UNIT, '(6f12.7)', end=667, err=666) &
            (box(i), i=1,6)
      deallocate(velocities)
      close(INPCRD_UNIT)
#ifdef BINTRAJ
   end if
#endif

   return

666 continue
   ierr = 1
   close(INPCRD_UNIT)
   if (allocated(velocities)) deallocate(velocities)
   return

667 continue
   ! We come here if we hit EOF trying to read the box, so the box must not
   ! exist. Just deallocate our velocities and leave
   box(:) = 0.d0
   close(INPCRD_UNIT)
   if (allocated(velocities)) deallocate(velocities)
   return

668 continue
   ! We come here if we read in 6 velocities but could not read any more. That
   ! means that the velocities were really the box. Copy them over to the box
   ! and bail out
   close(INPCRD_UNIT)
   if (allocated(velocities)) then
      box(:) = velocities(1:6)
      deallocate(velocities)
   else
      ierr = 1  ! This should never happen
   end if
   return

end subroutine read_inpcrd_file

! This subroutine extracts the number of atoms defined in an inpcrd file. This
! can be used to protect against buffer overruns. natom is set to -1 upon errors
subroutine get_inpcrd_natom(filename, natom)

#ifdef BINTRAJ
   use AmberNetcdf_mod
   use binrestart
#endif
   use file_io_dat, only : INPCRD_UNIT

   implicit none

   ! Passed variables

   character(len=*), intent(in) :: filename
   integer, intent(out)         :: natom

   ! Local variables

   integer           :: ncid
   integer           :: alloc_failed
   character(len=80) :: title
   double precision  :: time
   integer           :: id1, id2, id3

#ifdef BINTRAJ
   ! Make sure the NetCDF module doesn't scream at us...
   verbose_netcdf = .false.

   if (NC_checkRestart(filename)) then
      ! Get the number of atoms
      if (NC_openRead(filename, ncid)) then
         natom = -1
         return
      end if
      if (NC_setupRestart(ncid, title, natom, id1, id2, id3, time)) &
         natom = -1
      call NC_close(ncid)
      return
   else
#endif
      ! Try reading the raw file
      open(unit=INPCRD_UNIT, file=filename, status='OLD', form='FORMATTED', &
           iostat=alloc_failed)
      if (alloc_failed /= 0) then
         natom = -1
         return
      end if

      read(INPCRD_UNIT, '(a80)', end=666) title
      read(INPCRD_UNIT, *, err=666, end=666) natom
      close(INPCRD_UNIT)
#ifdef BINTRAJ
   end if
#endif

   return

666 continue
   natom = -1
   close(INPCRD_UNIT)
   return

end subroutine get_inpcrd_natom

#ifdef LES
end module sanderles_api
#else
end module sander_api
#endif

! Below here we provide the functions _outside_ the module so that the names do
! not become name-mangled (and therefore compiler-dependent). These are just
! wrappers around the module subroutines

! Sets up input options for a gas-phase (i.e., no periodic boundaries) system
!
! Parameters
! ----------
! inp : type(sander_input)
!     struct of input options that will be filled by this subroutine
! gb : integer (optional)
!     If GB model is desired, igb will be set to this.  If not present or 0, a
!     gas phase calculation is set up. If set to 10, a PB calculation is set up
!
! If gb is given, the given GB model is set up. Otherwise, the system is set up
! for a gas-phase simulation

#ifdef LES
#  define SANDER_API_MOD sanderles_api
#else
#  define SANDER_API_MOD sander_api
#endif

subroutine ext_gas_sander_input(inp, gb)

   use SANDER_API_MOD, only : mod_func => gas_sander_input, sander_input

   implicit none

   type(sander_input), intent(out) :: inp
   integer, optional, intent(in)   :: gb

   if (present(gb)) then
      call mod_func(inp, gb)
   else
      call mod_func(inp)
   end if

end subroutine ext_gas_sander_input

! Sets up input options for a periodic system treated with PME
!
! Parameters
! ----------
! inp : type(sander_input)
!     struct of input options that will be filled by this subroutine
subroutine ext_pme_sander_input(inp)

   use SANDER_API_MOD, only : mod_func => pme_sander_input, sander_input

   implicit none

   type(sander_input), intent(out) :: inp

   call mod_func(inp)

end subroutine ext_pme_sander_input

! Initializes a struct with QM options to all default values
!
! Parameters
! ----------
! inp : type(qmmm_input_options)
!     struct of QM input options that will be filled by this subroutine
subroutine ext_qm_sander_input(inp)

   use SANDER_API_MOD, only : mod_func => qm_sander_input, qmmm_input_options

   implicit none

   type(qmmm_input_options), intent(out) :: inp

   call mod_func(inp)

end subroutine ext_qm_sander_input

! Initializes the major data structures needed to evaluate energies and forces
!
! Parameters
! ----------
! prmname : string
!     Name of the topology file
! coordinates : double precision(3*natom)
!     The starting coordinates
! box : double precision(6)
!     The box dimensions
! input_options : type(sander_input)
!     struct of input options used to set up the calculation
! qmmm_options : type(qmmm_input_options) [optional]
!     struct of input options used to set up the QM/MM part of the calculation
! ierr : integer
!     Set to 0 for success or 1 if the setup failed
subroutine ext_sander_setup(prmname, coordinates, box, input_options, qmmm_options, ierr)

   use SANDER_API_MOD, only : mod_func => sander_setup, sander_input, &
                              qmmm_input_options, MAX_FN_LEN

   implicit none

   ! Input parameters
   character(len=MAX_FN_LEN) :: prmname
   double precision, dimension(6), intent(in) :: box
   double precision, dimension(*), intent(in) :: coordinates
   type(sander_input) :: input_options
   type(qmmm_input_options), optional :: qmmm_options
   integer, intent(out) :: ierr

   call mod_func(prmname, coordinates, box, input_options, qmmm_options, ierr)

   return

end subroutine ext_sander_setup

#ifndef NO_ALLOCATABLES_IN_TYPE
! Initializes the major data structures needed to evaluate energies and forces
!
! Parameters
! ----------
! parmdata : prmtop_struct
!     Name of the topology file
! coordinates : double precision(3*natom)
!     The starting coordinates
! box : double precision(6)
!     The box dimensions
! input_options : type(sander_input)
!     struct of input options used to set up the calculation
! qmmm_options : type(qmmm_input_options) [optional]
!     struct of input options used to set up the QM/MM part of the calculation
! ierr : integer
!     Set to 0 for success or 1 if the setup failed
subroutine ext_sander_setup2(parmdata, coordinates, box, input_options, qmmm_options, ierr)

   use SANDER_API_MOD, only : mod_func => sander_setup2, sander_input, &
                              qmmm_input_options, prmtop_struct

   implicit none

   ! Input parameters
   type(prmtop_struct), intent(in) :: parmdata
   double precision, dimension(6), intent(in) :: box
   double precision, dimension(*), intent(in) :: coordinates
   type(sander_input) :: input_options
   type(qmmm_input_options), optional :: qmmm_options
   integer, intent(out) :: ierr

   call mod_func(parmdata, coordinates, box, input_options, qmmm_options, ierr)

   return

end subroutine ext_sander_setup2

#else

subroutine ext_sander_setup2(parmdata, coordinates, box, input_options, qmmm_options, ierr)

   use SANDER_API_MOD, only : sander_input, qmmm_input_options, prmtop_struct

   implicit none

   ! Input parameters
   type(prmtop_struct), intent(in) :: parmdata
   double precision, dimension(6), intent(in) :: box
   double precision, dimension(*), intent(in) :: coordinates
   type(sander_input) :: input_options
   type(qmmm_input_options), optional :: qmmm_options
   integer, intent(out) :: ierr

   write(0,*) 'Compiler does not support allocatables in structs. Recompile'
   write(0,*) 'with a more modern compiler'

   return
end subroutine ext_sander_setup2

#endif /* NO_ALLOCATABLES_IN_TYPE */

! Sets the atomic positions
subroutine ext_set_positions(positions)

   use memory_module, only: natom
   use SANDER_API_MOD, only : mod_func => set_positions

   implicit none

   double precision, intent(in), dimension(natom*3) :: positions

   call mod_func(positions)

   return

end subroutine ext_set_positions

! Gets the periodic box dimensions
subroutine ext_get_box(a, b, c, alpha, beta, gamma)

    use SANDER_API_MOD, only : mod_func => get_box

    implicit none

    double precision, intent(out) :: a, b, c, alpha, beta, gamma

    call mod_func(a, b, c, alpha, beta, gamma)

    return

end subroutine ext_get_box

! Sets the periodic box vectors
subroutine ext_set_box(a, b, c, alpha, beta, gamma)

   use SANDER_API_MOD, only : mod_func => set_box

   implicit none

   double precision, intent(in) :: a, b, c, alpha, beta, gamma

   call mod_func(a, b, c, alpha, beta, gamma)

   return

end subroutine ext_set_box

! Gets the currently active positions
subroutine ext_get_positions(positions)

    use memory_module, only: natom
    use SANDER_API_MOD, only: mod_func => get_positions

    implicit none

    double precision, intent(out), dimension(natom*3) :: positions

    call mod_func(positions)

    return

end subroutine ext_get_positions

! Computes the energies and forces with the current positions
subroutine ext_energy_forces(energy, forces)

   use SANDER_API_MOD, only : potential_energy_rec, mod_func => energy_forces
   use memory_module, only : natom

   implicit none

   ! Input parameters
   type(potential_energy_rec), intent(out) :: energy
   double precision, dimension(natom*3), intent(out) :: forces

   call mod_func(energy, forces)

   return

end subroutine ext_energy_forces

subroutine ext_sander_cleanup()

   use SANDER_API_MOD, only : mod_func => sander_cleanup

   implicit none

   call mod_func

   return

end subroutine ext_sander_cleanup

! Returns the number of atoms defined in the current system
subroutine ext_sander_natom(n)

    use SANDER_API_MOD, only : mod_func => sander_natom

    implicit none

    integer, intent(out) :: n

    call mod_func(n)

    return

end subroutine ext_sander_natom

! Parses the input coordinate file and properly determines whether velocities
! and/or box are present
subroutine ext_read_inpcrd_file(filename, coordinates, box, ierr)

   use SANDER_API_MOD, only : mod_func => read_inpcrd_file, MAX_FN_LEN

   implicit none

   character(len=MAX_FN_LEN), intent(in) :: filename
   double precision, intent(out) :: coordinates(*)
   double precision, intent(out) :: box(6)
   integer, intent(out)          :: ierr

   call mod_func(filename, coordinates, box, ierr)

   return

end subroutine ext_read_inpcrd_file

! Gets the number of atoms from an inpcrd file
subroutine ext_get_inpcrd_natom(filename, natom)

   use SANDER_API_MOD, only : mod_func => get_inpcrd_natom, MAX_FN_LEN

   implicit none

   character(len=MAX_FN_LEN), intent(in) :: filename
   integer, intent(out) :: natom

   call mod_func(filename, natom)

   return

end subroutine ext_get_inpcrd_natom

! Read a topology file and fill a prmtop_struct with it
subroutine ext_read_prmtop_file(filename, parm, ierr)

    use SANDER_API_MOD, only: mod_func => read_prmtop_file, prmtop_struct, &
                              MAX_FN_LEN

    implicit none

    character(len=MAX_FN_LEN), intent(in) :: filename
    type(prmtop_struct)  :: parm
    integer, intent(out) :: ierr

    ierr = 0

    call mod_func(filename, parm, ierr)

end subroutine ext_read_prmtop_file

! Sets setup to 0 if sander is not set up and 1 if it is
subroutine ext_is_setup(setup)

   use SANDER_API_MOD, only: is_setup

   implicit none

   integer, intent(out) :: setup

   setup = 0
   if (is_setup()) setup = 1

   return

end subroutine ext_is_setup
