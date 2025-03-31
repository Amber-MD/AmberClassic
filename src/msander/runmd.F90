!<compile=optimized>
#include "../include/dprec.fh"
#include "../include/assert.fh"
#include "nfe-config.h"

!  In vim, :set foldmethod=marker, then zo,zc, to open/close single folds
!          zO,zC to open close folds recursively
!          zM closes all folds; zR opens all folds

!------------------------------------------------------------------------------
! runmd: main driver routine for molecular dynamics.  This is over 2300 lines
!        long and incorporates virtually every method that sander offers.
!        Prior to calling this routine, the sander subroutine itself will take
!        user input, read system specifications from coordinates and topology
!        files, and lay out a series of arrays to hold the information that
!        runmd then operates with.
!
! Units:
!   Runmd operates in kcal/mol units for energies, amu for masses, and
!   Angstroms for distances.  To convert the input time parameters from
!   picoseconds to internal units, multiply by 10.0*sqrt(4.184) = 20.455.
!
! Arguments:
!   xx:        global real array. See locmem.F90 for structure/pointers.
!   ix:        global integer array. See locmem.F90 for structure/pointers.
!   ih:        global Hollerith array holding atom names, residues names, atom
!              types, and other information (see the indexing fully described
!              in locmem.F90)
!   ipairs:    pairlist of nonbonded interactions
!   x:         global position array *
!   winv:      array with inverse masses *
!   amass:     mass array *
!   f:         force array, used to hold old coordinates temporarily, too
!   v:         velocity array
!   vold:      old velocity array, from the previous step
!   xc:        array of reals, matching the size of x itself, used for scratch
!              space in various subroutine calls
!   conp:      bond parameters for SHAKE
!   skip:      logical skip array for SHAKE (and QM/MM too, I think)
!   nsp:       submolecule index array (?)
!   erstop:    should we stop in error (?)
!   qsetup:    Flag to activate setup of multiple components, .false. on
!              first call
!------------------------------------------------------------------------------

module runmd_module

!------------------------------------------------------------------------------
! modules used:  {{{

  use state

#if !defined(DISABLE_NFE) && defined(NFE_ENABLE_BBMD)
  use nfe_sander_hooks, only : nfe_on_mdstep => on_mdstep
#endif
  use nfe_sander_proxy, only : infe, nfe_real_mdstep, nfe_prt

  use commandline_module, only: cpein_specified
  use md_scheme, only: thermostat_step, ntt
  use sander_rism_interface, only: rismprm, RISM_NONE, RISM_FULL, &
                                   RISM_INTERP, rism_calc_type, &
                                   rism_solvdist_thermo_calc, mylcm

  use file_io_dat
  use constants, only: third, ten_to_minus3
  use stack

  use fastwt
  use bintraj, only: end_binary_frame
  use nblist,only: fill_tranvec,volume,oldrecip,ucell
  use md_scheme, only: gamma_ln

  use sgld, only: isgld, sgenergy, sgldw, sgmdw

  use sgld, only: isgsta,isgend,sg_fix_degree_count

#ifdef MPI
  use remd, only: rem, mdloop, remd_ekmh, repnum, stagid, my_remd_data, &
                  hybrid_remd_ene, next_rem_method, remd_types, replica_indexes, &
                  my_pressure, my_volume
  use softcore, only: ifsc, sc_dvdl, sc_tot_dvdl, sc_tot_dvdl_partner, &
                      sc_dvdl_ee, sc_tot_dvdl_ee, sc_tot_dvdl_partner_ee, &
                      extra_atoms, mix_temp_scaling, sc_pscale, &
                      adj_dvdl_stat, sc_mix_velocities, sc_nomix_frc, &
                      sc_sync_x, sc_print_energies, calc_softcore_ekin, &
                      sc_ener, sc_ener_ave, sc_ener_rms, sc_lngdyn, &
                      sc_ener_tmp, sc_ener_tmp2, sc_ener_old, sc_ener_old2, &
                      sc_mix_position, sc_print_dvdl_values, &
                      sc_degrees_o_freedom, dynlmb, sc_change_clambda, &
                      ti_ene_cnt, sc_compare
  use mbar, only : ifmbar, bar_intervall, calc_mbar_energies, &
                   bar_collect_cont, do_mbar
  use mpi
#endif /* MPI */

  use emap, only:temap,emap_move
  use barostats, only : mcbar_trial, mcbar_summary

  ! use memory_module, only: mass
  use random

#ifdef MPI
  use sgld, only : trxsgld
#endif

  ! SINR (Stochastic Isokinetic Nose-Hoover RESPA integrator)
  use sinr_t

  ! }}}

  ! Local variables  {{{
  !  factt       : degree-of-freedom correction factor for temperature scaling
  !  nr          : local copy of nrp, number of atoms
  !  nr3         : 3 * nr, used for runtime efficiency
  !
  ! Common memory variables
  !  nrp         : number of atoms, adjusted for LES copies

!------------------------------------------------------------------------------
! module variables:

  character(kind=1,len=5) :: routine="runmd"
#ifdef MPI
#  include "parallel.h"
  _REAL_ mpitmp(8) !Use for temporary packing of mpi messages.
  integer ist(MPI_STATUS_SIZE), partner
#else
  ! mdloop and REM is always 0 in serial
  integer, parameter :: mdloop = 0, rem = 0, remd_types(1) = 0, replica_indexes(1) = 0
  integer, parameter :: numtasks = 1
  integer, parameter :: mytaskid = 0
#endif

  ! The following variables are needed since nstep and nstlim behave
  ! differently in a REMD run.  In certain places where output is required,
  ! total_nstep and total_nstlim take the place of nstep and nstlim. This
  ! allows replica runs to output less often than every exchange.  They are
  ! the absolute step # of the REMD or MD simulation.
  integer total_nstep, total_nstlim

  ! Stochastic Isokinetic Nose-Hoover RESPA integrator (SINR)
  type(sinr) :: sinrdata ! Variables for SINR
  _REAL_ gammai, sinrtaui

#include "../include/md.h"
#include "box.h"
#include "nmr.h"
#include "tgtmd.h"
#include "multitmd.h"
#include "../include/memory.h"
#include "extra.h"
#include "ew_frc.h"
#include "ew_cntrl.h"
#include "ew_mpole.h"
#include "def_time.h"
#include "extra_pts.h"

  _REAL_ sgsta_rndfp, sgend_rndfp, ignore_solvent
  _REAL_ sysx, sysy, sysz, sysrange(3,2)
  logical mv_flag

  integer iatom, ierr
  ! Workaround for gfortran optimization bug that retains support for
  ! gfortran 4.2 and lower
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
  integer, volatile :: m
#else
  integer :: m
#endif

  _REAL_  Ekin2_tot,tmp
  integer :: idim
  _REAL_ :: exp1, exp2

  logical ivscm
  logical resetvelo
  _REAL_ etot_save,ekpbs

  logical do_list_update
  logical belly, lout, loutfm, erstop, vlim, onstep
  ! Fortran does not guarantee short circuit logical expressions:
  type(state_rec) :: ener   ! energy values per time step
  type(state_rec) :: enert  ! energy values tallied over the time steps
  type(state_rec) :: enert2 ! energy values squared tallied over the time steps
  type(state_rec) :: enert_old, enert2_old
  type(state_rec) :: enert_tmp, enert2_tmp
  type(state_rec) :: edvdl
  type(state_rec) :: edvdl_r

#ifdef MPI
  type(state_rec) :: ecopy
  _REAL_ :: clfac, tmpvir(3,3)
#endif
  _REAL_ rmu(3), fac(3), onefac(3), etot_start
  _REAL_ tspan, atempdrop, fln, scaltp
  _REAL_ vel, vel2, vcmx, vcmy, vcmz, vmax
  _REAL_ winf, aamass, rterm, ekmh, ekph, wfac, rsd
  _REAL_ fit, fiti, fit2

  ! Variables to control a Langevin dynamics simulation
  _REAL_ ekins0
  _REAL_ dtx, dtxinv, dt5, factt, ekin0, ekinp0, dttp
  _REAL_ rndf, rndfs, rndfp, boltz2, pconv, tempsu
  _REAL_ xcm(3), acm(3), ocm(3), vcm(3), ekcm, ekrot
  _REAL_ emtmd

  ! Variables and parameters for constant surface tension:
  ! ten_conv converts dyne/cm to bar angstroms
  _REAL_, parameter :: ten_conv = 100.0d0
  _REAL_  :: pres0x
  _REAL_  :: pres0y
  _REAL_  :: pres0z
  _REAL_  :: gamma_ten_int
  _REAL_  :: press_tan_ave

  integer idumar(4)
  integer l_temp
  integer i, j, im, i3, nitp, nits, iskip_start, iskip_end
  integer nstep, nrep, nrek, iend, istart3, iend3
  integer nrx, nr, nr3, ntcmt, izero, istart
  logical ixdump, ivdump, itdump, ifdump
  logical irismdump

  integer nvalid, nvalidi
  _REAL_ eke

  _REAL_, allocatable, dimension(:) :: frcti

  _REAL_ small
  data small/1.0d-7/

  !--- VARIABLES FOR DIPOLE PRINTING ---
  integer prndipngrp
  integer prndipfind
  character(len=4) prndiptest

  _REAL_,parameter :: pressure_constant = 6.85695d+4

  ! variables used in middle scheme
  _REAL_, allocatable :: xold(:)       ! for shake
  ! Ek(t) used to be written to mdout
  ! Ek(t+dt/2) is now written to mdout
  ! DAC: need to figure out exactly how the above matches the code
  character(len=*),parameter :: file_ek = "eke.dat"
  _REAL_ :: ekhf, ekhf2

  _REAL_ :: kinetic_E_save(2)

  ! PLUMED related variables
  _REAL_ :: plumed_box(3,3), plumed_virial(3,3), plumed_kbt
  integer :: plumed_version, plumed_stopflag
  _REAL_ :: plumed_energyUnits, plumed_timeUnits, plumed_lengthUnits
  _REAL_ :: plumed_chargeUnits
  ! }}}

private
public :: runmd

contains

subroutine runmd(xx, ix, ih, ipairs, x, winv, amass, f, v, vold, xc, &
                 conp, skip, nsp, erstop, qsetup)

  implicit none
  integer, intent(in) ::   ipairs(:), ix(*), nsp(*)
  _REAL_, intent(inout) ::  xx(*)
  character(len=4), intent(in) :: ih(*)
  _REAL_, intent(inout) ::  x(:), winv(:), amass(:), f(:), v(:), vold(:), &
                            xc(:), conp(:)
  logical, intent(inout) ::  erstop, qsetup
  _REAL_, intent(in) ::  skip(:)   ! N.B.: skip is really a logical variable,
                                   ! but we are just really passing an opaque
                                   ! pointer here
  ! Used to save averaged velocity information for on-step output
  _REAL_ onstepvel(3*natom)
  ! Used to save coordinates at time=t for on-step output
  _REAL_ onstepcrd(3*natom)
  ! scratch for ndvptx
  integer iscrth(natom)


!------------------------------------------------------------------------------
!  execution/initialization begins here:

  call initialize_runmd(x,ix,v)
  onstepvel(:) = 0.d0
  onstepcrd(:) = 0.d0

  if (init == 3 .or. nstlim == 0) then
  ! Make a first dynamics step. {{{
  ! init = 3: general startup if not continuing a previous run

!----------------------------------------------------------------------------

     !  Correct the input velocities for SHAKE:
     if (ntc /= 1) then
        ! RATTLE-V, correct velocities
        call rattlev(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
        winv,conp,skip,x,v,nitp,belly,ix(iifstwt),ix(noshake))
        ! use SETTLE to deal with water model
        call quick3v(x, v, ix(iifstwr), natom, nres, ix(i02))
     end if  

    ! Calculate the force.  Set irespa to get full
    ! energies calculated on step "0":
    irespa = 0
    iprint = 1
    call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, xx(l96), xx(l97), &
               xx(l98), xx(l99), qsetup, do_list_update, nstep)
    call modwt_reset()  ! since TEMP0 might have been updated

    ! This force call does not count as a "step". CALL NMRDCP to decrement
    ! local NMR step counter and MTMDUNSTEP to decrease the local MTMD step
    ! counter
    call nmrdcp
    call mtmdunstep

    ! PLUMED force is added in this routine.
    plumed_stopflag=0
#ifdef PLUMED
    if (plumed == 1) then
#     include "Plumed_force.inc"
    end if
#endif

#ifdef MPI
    onstep = .false.
    call thermodynamic_integration(f)
#endif
    call modwt_reset()  ! since TEMP0 might have been updated
    irespa = 1
    if (ntp > 0) then
      ener%volume = volume
      ener%density = tmass / (0.602204d0*volume)
      ener%cmt(4) = 0.d0
      ener%vir(4) = 0.d0
      ener%pres(4) = 0.d0
      do m = 1,3
        ener%cmt(m)  = ener%cmt(m) * 0.5d0
        ener%cmt(4)  = ener%cmt(4) + ener%cmt(m)
        ener%vir(4)  = ener%vir(4) + ener%vir(m)
        ener%pres(m) = (pconv+pconv) * (ener%cmt(m)-ener%vir(m)) / volume
        ener%pres(4) = ener%pres(4) + ener%pres(m)
      end do
      ener%pres(4) = ener%pres(4) / 3.d0
    end if
    ntnb = 0
    i3 = 0
    tempsu = 0.0d0

!----------------------------------------------------------------------------
    ! update the velocities; only real atoms are handled in this loop:
    !  TODO:  fdist has been called, but I think each MPI process only
    !         knows the forces assigned to its atoms(?)
    do j = 1,nrp
      winf = winv(j) * dt5
      aamass = amass(j)
      do m = 1,3
        i3 = i3+1
        rterm = v(i3)*v(i3) * aamass
        tempsu = tempsu + rterm
        v(i3) = v(i3) - f(i3) * winf
        if (vlim) v(i3) = sign(min(abs(v(i3)), vlimit), v(i3))
      end do
    end do

!----------------------------------------------------------------------------
    ! Now do velocity updates for the extra variables
    !       only final node has the correct forces for these extra variables
    do im = 1, iscale
       tempsu = tempsu + scalm * v(nr3+im)*v(nr3+im)
       if( mytaskid == numtasks -1 )  &
          v(nr3+im) = v(nr3+im) - f(nr3+im) * dt5 / scalm
    end do
    ener%kin%solt = tempsu * 0.5d0

#ifdef MPI /* SOFT CORE */
    if (ifsc /= 0) then
      call calc_softcore_ekin(amass, v, v, istart, iend)
      sc_ener(13) = sc_ener(6) + sc_ener(12)
    end if
#endif

!----------------------------------------------------------------------------
    ! middle scheme for constrained MD
    if (ntc /= 1) then
       ! RATTLE-V, correct velocities
       call rattlev(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
       winv,conp,skip,x,v,nitp,belly,ix(iifstwt),ix(noshake))
       ! use SETTLE to deal with water model
       call quick3v(x, v, ix(iifstwr), natom, nres, ix(i02))
    end if  
    ener%kin%tot = ener%kin%solt
    ener%tot = ener%kin%tot + ener%pot%tot

  ! This ends a HUGE conditional branch in which init == 3, general startup
  ! when not continuing a previous dynamics run.  }}}
  end if

!------------------------------------------------------------------------------
  ! Continue init=3 or start init=4: {{{
  ! What follows applies to the case of init = 4: continuation of a previous
  ! trajectory.  This will also done for init = 3 (that code is simply run as
  ! startup and effectively performs a dynamics step to prime the pump).
  !
  ! Note: if the last printed energy from the previous trajectory was at time
  ! "t", then the restrt file has velocities at time t + 0.5dt and
  ! coordinates at time t + dt.
  ekmh = 0.0d0

  ! kinetic energy calculation:
  i3 = 0
  do j = 1, nrp
    aamass = amass(j)
    do m = 1, 3
      i3 = i3+1
      rterm = v(i3) * v(i3) * aamass
      ekmh = ekmh + rterm
    end do
  end do

#ifdef MPI /* SOFT CORE */
  if (ifsc /= 0) then
    call calc_softcore_ekin(amass, v, v, istart, iend)
    sc_ener(13) = sc_ener(6) + sc_ener(12)
  end if
#endif /* MPI */

   do im = 1, iscale
      ekmh = ekmh + scalm*v(nr3+im)*v(nr3+im)
   end do
   ekmh = ekmh * 0.5d0

  vold(1:nr3+iscale) = v(1:nr3+iscale)

  ! }}}

!------------------------------------------------------------------------------
  if (init .ne. 4 .or. nstlim == 0) then
  ! More setup {{{

    ! Print the initial energies and temperatures

    if (rismprm%rism == 1 .and. rismprm%write_thermo==1 .and. &
        nstep <= 0 .and. facc .ne. 'A') then
      if (rism_calc_type(0) == RISM_FULL) then
        if (nstlim == 0) then
          call rism_solvdist_thermo_calc(.true., 0)
        else
          call rism_solvdist_thermo_calc(.false., 0)
        end if
      end if
    end if

    if (nstep <= 0 .and. master .and. facc .ne. 'A') then
      if (isgld > 0) call sgenergy(ener)
      rewind(7)

      call prntmd(nstep, t, ener, onefac, 7, .false.)
#ifdef MPI /* SOFT CORE */
      if (ifsc .ne. 0) call sc_print_energies(6, sc_ener)
      if (ifsc .ne. 0) call sc_print_energies(7, sc_ener)
#endif
      if (nmropt > 0) call nmrptx(6)
      if (infe == 1) call nfe_prt(6)
      call flush(7)
    end if

!----------------------------------------------------------------------------
    ! Clean exit for a zero-step run:
    if (nstlim == 0) then
#ifdef MPI
      call xdist(x, xx(lfrctmp), 3*natom+iscale)
      call xdist(v, xx(lfrctmp), 3*natom+iscale)
#endif
      if (master) then
        call mdwrit(nstep, nr, ntxo, ntb, x, v, t, temp0,solvph,solve)
        if (ntwx>0) call corpac(x, 1, nrx, MDCRD_UNIT, loutfm)
        if (ntwv>0) call corpac(v, 1, nrx, MDVEL_UNIT, loutfm)
        if (ntwe>0) call mdeng(15, nstep, t, ener, onefac, ntp)
      end if
      return
    end if
    init = 4
  ! End of contingencies primarily related to init not equal to 4 }}}
  end if

!------------------------------------------------------------------------------

#ifdef MPI
  ! If this is a replica run and we are on exchange > 1, restore the
  ! old ekmh value since it was reset after we left runmd last time.
  if (rem /= 0 .and. mdloop >= 1) ekmh = remd_ekmh
#endif

  if (ntt == 10) then
     !   initialize the SINR integrator: {{{
    ! convert from ps time units to internal times:
    gammai = gamma_ln/20.455d0
    sinrtaui = sinrtau*20.455d0

#ifdef MPI
    call sinr_init(natom, nkija, dtx, boltz2, temp0, gammai, sinrtaui, &
                   sinrdata, commsander)
#else
    call sinr_init(natom, nkija, dtx, boltz2, temp0, gammai, sinrtaui, &
                   sinrdata)
#endif
    if (irest == 1) then
      call sinr_read_vels(v, sinrtaui, sinrdata)
    else
      call init_sinr_vels(v, amass, sinrdata)
    endif

#ifdef MPI
    if (numtasks > 1) then
      call sinr_mpi_init(v,numtasks,sinrdata)
    end if
#endif
    call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, &
               xx(l96), xx(l97), xx(l98), xx(l99), qsetup, &
               do_list_update,nstep)
    call iLndt(v, amass, istart, iend, sinrdata)
    call iLvdt(v, f, amass, istart, iend, sinrdata)
    call iLudt(x, v, istart, iend, sinrdata)
    ! }}}
  endif


!------------------------------------------------------------------------------
  ! The main loop for performing dynamics steps: at this point, the
  ! coordinates are a half-step "ahead" of the velocities; the variable
  ! EKMH holds the kinetic energy at these "-1/2" velocities, which are
  ! stored in the array vold.
  260 continue
  onstep = mod(irespa,nrespa) == 0
  iprint = 0
  if (nstep == 0 .or. nstep+1 == nstlim) iprint = 1
  if (rem .eq. 0 .or. mdloop .gt. 0) nfe_real_mdstep = .True.
#ifdef MPI
  ! Set do_mbar for the force contributions
  if (ifmbar /= 0) then
    do_mbar = .false.
    if (mod(nstep+1, bar_intervall) == 0) do_mbar = .true.
  end if
#endif

!------------------------------------------------------------------------------
  ! Step 1a: initial trial move for MC barostat:
  if (ntp > 0 .and. barostat == 2 .and. mod(total_nstep+1, mcbarint) == 0) &
    call mcbar_trial(xx, ix, ih, ipairs, x, xc, f, ener%vir, xx(l96), &
                     xx(l97), xx(l98), xx(l99), qsetup, do_list_update, &
                     nstep, nsp, amass, nmropt)

!------------------------------------------------------------------------------
  ! Step 1b: get the forces for the system's current coordinates
  !     [This is where the force() routine mainly gets called:]
  call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, xx(l96), xx(l97), &
             xx(l98), xx(l99), qsetup, do_list_update, nstep)
#ifdef PLUMED
  if (plumed == 1) then
#     include "Plumed_force.inc"
  end if
#endif

#ifdef MPI
   call thermodynamic_integration(f)
#endif

!------------------------------------------------------------------------------
  ! Reset quantities depending on TEMP0 and TAUTP  {{{
  ! (which may have been changed by MODWT during FORCE call).
  ekinp0 = fac(2)*temp0
  ekins0 = fac(3)*temp0
  ekin0 = fac(1)*temp0
  ! }}}

  ! Constant pressure conditions: {{{
  if (ntp > 0) then
    ener%volume = volume
    ener%density = tmass / (0.602204d0*volume)
    ener%cmt(4) = 0.d0
    ener%vir(4) = 0.d0
    ener%pres(4) = 0.d0
    do m = 1,3
      ener%cmt(m)  = ener%cmt(m)*0.5d0
      ener%cmt(4)  = ener%cmt(4) + ener%cmt(m)
      ener%vir(4)  = ener%vir(4) + ener%vir(m)
      ener%pres(m) = (pconv + pconv) * (ener%cmt(m) - ener%vir(m)) / volume
      ener%pres(4) = ener%pres(4) + ener%pres(m)
    end do
    ener%pres(4) = ener%pres(4) / 3.d0
  end if
  ! End contingency for constant pressure conditions }}}

!------------------------------------------------------------------------------
#ifdef MPI
  ! Replica Exchange Molecular Dynamics {{{
  ! if rem /= 0 and mdloop == 0, this is
  ! the first sander call and we don't want to actually do any MD or change the
  ! initial coordinates.  Exit here since we only wanted to get the potential
  ! energy for the first subrem exchange probability calc.
  if (rem /= 0 .and. mdloop == 0) then
#  ifdef VERBOSE_REMD
    if (master) &
      write (6,'(a,i3)') '| REMD: Exiting runmd after getting initial energies &
                          &for replica', repnum
#  endif /* VERBOSE_REMD */

    ! This diverts all the way to the end of the runmd loop
    goto 480
  endif
  ! End contingency for the first sander call in REMD
  ! (rem is not 0, and mdloop == 0)

  ! }}}
#endif /* MPI */

!------------------------------------------------------------------------------
  ! Step 2: Do the velocity update: {{{
  ! Step 2a: apply quenched MD if needed.  This is useful in NEB>0
  call timer_start(TIME_VERLET)

  if (vv == 1) call quench(f, v)

  if (isgld > 0) then
    call sgldw(natom, istart, iend, dtx, temp0, ener, amass, winv, &
               x, f, v)
  else if (ntt == 10) then

  ! Stochastic Isokinetic Nose-Hoover RESPA (SINR) integrator)
    call iLvdt(v, f, amass, istart, iend, sinrdata)
    call iLndt(v, amass, istart, iend, sinrdata)
    call iLndt(v, amass, istart, iend, sinrdata)
    call iLvdt(v, f, amass, istart, iend, sinrdata)
    call iLudt(x, v, istart, iend, sinrdata)
    if (mod(nstep+1, ntpr*nrespa) == 0) then
      call sinr_temp(v, amass, istart3, iend3, sinrdata)
    end if
    if (mod(nstep+1, ntwr) == 0 .or. nstep+1 == nstlim) then
      call sinr_write_vels(v, nstep, istart3, iend3, sinrtaui, sinrdata)
    end if

  else
  ! leap-frog middle scheme
  ! the 1st step for updating p
     i3 = 3*(istart - 1)
     do j = istart, iend
       wfac = winv(j) * dtx
       v(i3+1) = v(i3+1) + f(i3+1)*wfac
       v(i3+2) = v(i3+2) + f(i3+2)*wfac
       v(i3+3) = v(i3+3) + f(i3+3)*wfac
       i3 = i3+3
     end do

     ! for constrained MD
     if (ntc /= 1) then
       ! RATTLE-V, correct velocities
       call rattlev(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
       winv,conp,skip,x,v,nitp,belly,ix(iifstwt),ix(noshake))

       ! use SETTLE to deal with water model
       call quick3v(x, v, ix(iifstwr), natom, nres, ix(i02))
     end if

  end if

  ! Consider vlimit
  if (vlim) then
    vmax = 0.0d0
    do i = istart3, iend3
      vmax = max(vmax, abs(v(i)))
      v(i) = sign(min(abs(v(i)), vlimit), v(i))
    end do

    ! Only violations on the master node are actually reported
    ! to avoid both MPI communication and non-master writes.
    if (vmax > vlimit) then
      if (master) &
        write(6,'(a,i6,a,f10.4)') 'vlimit exceeded for step ', nstep, &
              '; vmax = ', vmax
    end if
  end if

  !  Update velocities for the "extra" variables
  if( mytaskid == numtasks - 1 ) then
     do im = 1, iscale
        v(nr3+im) = (v(nr3+im) + f(nr3+im)*dtx/scalm)
     end do
  endif
  ! }}}

!------------------------------------------------------------------------------
  ! Update EMAP rigid domains
  if (temap) call emap_move()

  ! Determine if we need to save the on-step coordinates
  if (isynctraj.ne.0) then
    ! Decision to write trajectory coords
    if (ntwx > 0 .and. (mod(total_nstep+1, ntwx) == 0)) then
      do i3 = istart3, iend3
        onstepcrd(i3) = x(i3)
      end do
    end if
  endif

!------------------------------------------------------------------------------
  ! Step 3: update the positions, and apply the thermostat,  {{{
    ! the step for updating x-T-x
  if (ntt == 10) then
    do i3 = istart3, iend3
      f(i3) = x(i3)
    end do
  else   ! default middle scheme code follows:
    do i3 = istart3, iend3
       f(i3) = x(i3)  ! temporarily(?) stores the starting coordinates
       x(i3) = x(i3) + v(i3)*dt5
    end do

    ! DAC note: at this point, v is at the same point as x: are these the
    !   velocities and positions we want to write if ntwv=-1?
    !   also: is this not the place to compute the kinetic energy?

    ! for random number, controled by ig
    if (no_ntt3_sync == 1) then
      iskip_start = 0
      iskip_end = 0
    else
      iskip_start = 3*(istart-1)
      iskip_end = 3*(nr-iend)
    endif
    !  note: following is a no-op if ntt=0
    call thermostat_step(v, winv, dtx, istart, iend, iskip_start,iskip_end, &
          iscale)
    do i3 = istart3, iend3
       x(i3) = x(i3) + v(i3)*dt5   ! completes the full-step update of x
    end do
  endif

  call timer_stop(TIME_VERLET)
  ! }}}

!------------------------------------------------------------------------------
  if (ntc .ne. 1) then
  ! Step 4a: if shake is being used, update the positions {{{
    call timer_start(TIME_SHAKE)
    xold(istart3:iend3) = x(istart3:iend3)
    call shake(nrp, nbonh, nbona, 0, ix(iibh), ix(ijbh), ix(ibellygp), &
               winv, conp, skip, f, x, nitp, belly, ix(iifstwt), &
               ix(noshake))
    call quick3(f, x, ix(iifstwr), natom, nres, ix(i02))
    if (nitp == 0) then
      erstop = .true.
      goto 480
    end if

    ! Need to synchronize coordinates after shake for TI
#ifdef MPI
    if (icfe .ne. 0) then
      call timer_barrier( commsander )
      call timer_stop_start(TIME_SHAKE,TIME_DISTCRD)
      if (.not. mpi_orig .and. numtasks > 1) &
        call xdist(x, xx(lfrctmp), nr3+iscale)

      ! In dual-topology this is done within softcore.f
      if (ifsc .ne. 1) then
        if (master) call mpi_bcast(x, nr3, MPI_DOUBLE_PRECISION, &
                         0, commmaster, ierr)
      else
        if (master) call sc_sync_x(x, nr3)
      end if
      if (numtasks > 1) &
        call mpi_bcast(x, nr3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
      call timer_stop_start(TIME_DISTCRD, TIME_SHAKE)
    end if
#endif  /* MPI */

!------------------------------------------------------------------------------
    ! Step 4b: Now fix the velocities and calculate KE.
    ! Re-estimate the velocities from differences in positions.
    !  DAC: is this still needed?  Might not rattlev do the job?
    v(istart3:iend3) = v(istart3:iend3) &
        + (x(istart3:iend3) - xold(istart3:iend3))*dtxinv

    ! RATTLE-V, correct velocities
    call rattlev(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
      winv,conp,skip,x,v,nitp,belly,ix(iifstwt),ix(noshake))

    ! use SETTLE to deal with water model
    call quick3v(x, v, ix(iifstwr), natom, nres, ix(i02))

    call timer_stop(TIME_SHAKE)
  ! }}}
  end if
  call timer_start(TIME_VERLET)

!------------------------------------------------------------------------------
  if (onstep) then
  ! Step 4c: get the KE, averaging:  {{{
    eke = 0.d0
    ekph = 0.d0
    ekpbs = 0.d0

    ! LF-Middle: use velocity v(t) for KE calculation {{{
      i3 = 3*(istart-1)
      do j = istart, iend
        aamass = amass(j)
        do m = 1, 3
          i3 = i3 + 1
          eke = eke + aamass*0.25d0*(v(i3) + vold(i3))**2
          ! Try pseudo KE from Eq. 4.7b of Pastor, Brooks & Szabo,
          ! Mol. Phys. 65, 1409-1419 (1988):
          ekpbs = ekpbs + aamass*v(i3)*vold(i3)
          ekph = ekph + aamass*v(i3)**2
        end do
      end do

    ! Sum up the partial kinetic energies:
#ifdef MPI
    if (.not. mpi_orig .and. numtasks > 1) then
      mpitmp(1) = eke
      mpitmp(2) = ekph
      mpitmp(3) = ekpbs
      call mpi_allreduce(MPI_IN_PLACE, mpitmp, 3, MPI_DOUBLE_PRECISION, &
                         mpi_sum, commsander, ierr)
      eke = mpitmp(1)
      ekph = mpitmp(2)
      ekpbs = mpitmp(3)
    end if
    ! }}}

    ! Calculate Ekin of the softcore part of the system
    if (ifsc .ne. 0) then
      call calc_softcore_ekin(amass, v, vold, istart, iend)
      sc_ener(13) = sc_ener(6) + sc_ener(12)
    end if
#endif

    ! All processors handle the "extra" variables:
    do im = 1, iscale
      eke = eke + scalm*0.25d0*(v(nr3+im) + vold(nr3+im))**2
      ekpbs = ekpbs + scalm*v(nr3+im)*vold(nr3+im)
      ekph = ekph + scalm*v(nr3+im)**2
    end do
    eke = eke * 0.5d0
    ekph = ekph * 0.5d0
    ekpbs = ekpbs * 0.5d0
  ! End contingency for onstep; end of step 4c
  ! }}}
  end if

!------------------------------------------------------------------------------
  ! Step 5: several tasks related to dumping of trajectory information  {{{
  !
  ! Determine if trajectory, velocity, or restart writing is
  ! imminent, or if the center of mass motion will be removed.
  ! These requires distribution of velocities or dipoles to
  ! all processes by the subroutine xdist in parallel runs.
  !
  ! Modified so that when running REMD, writing can occur less often
  ! than exchanges (e.g. ntwx > nstlim).  Two new variables, total_nstep
  ! and total_nstlim were added.  For non-REMD runs,
  ! total_nstep = nstep + 1 and total_nstlim = nstlim as before.
  !
  ! For REMD runs, total_nstep = (mdloop-1)*nstlim + nstep + 1, where
  ! mdloop is the current exchange - this is the current
  ! replica exchange MD step.  total_nstlim = numexchg*nstlim, which
  ! is the maximum number of REMD steps.
  total_nstep = nstep + 1
  total_nstlim = nstlim

#ifdef MPI
  if (rem .ne. 0) then
    total_nstep = (mdloop - 1) * nstlim + nstep + 1
    total_nstlim = nstlim * numexchg
  endif
#endif

  ! Decision to write trajectory coords
  itdump = .false.
  if (ntwx > 0) itdump = (mod(total_nstep, ntwx) == 0)
  ! Decision to write velocities
  ivdump = .false.
  if (ntwv > 0) ivdump = (mod(total_nstep,ntwv) == 0)
  ! Determine whether we want to save the on-step velocities.
  if ((isynctraj.ne.0).and.ivdump) then
    i3 = (3*(istart-1)) + 1
    do j = istart, iend
      onstepvel(i3  ) = (v(i3  ) + vold(i3  )) / 2.d0
      onstepvel(i3+1) = (v(i3+1) + vold(i3+1)) / 2.d0
      onstepvel(i3+2) = (v(i3+2) + vold(i3+2)) / 2.d0
      i3 = i3 + 3
    enddo
  endif

  ! Decision to write forces
  ifdump = .false.
  if (ntwf > 0) ifdump = (mod(total_nstep,ntwf) == 0)

  ! Decision to write a restart file, or the final restart file
  ixdump = .false.
  if ( ntwr /= 0) ixdump = mod(total_nstep, ntwr ) == 0
  if (total_nstep >= total_nstlim) ixdump = .true.

  ! Decision to remove velocity of the system center of mass
  ivscm  = .false.
  if (nscm > 0) ivscm = mod(total_nstep,nscm) == 0

  ! Combined coordinate and velocity file writing
  if (ntwv == -1 .and. itdump) ivdump = .true.

  ! Write RISM files this step?
  irismdump = .false.
  if (rismprm%rism == 1) then
    if (rismprm%ntwrism > 0) then
      irismdump = (mod(nstep+1,rismprm%ntwrism) == 0)
      if (nstep + 1 >= nstlim) irismdump = .true.
    end if
  end if

#ifdef MPI
!------------------------------------------------------------------------------
  ! Distribute the coordinates, dipoles, and velocities as necessary
  call timer_stop_start(TIME_VERLET,TIME_DISTCRD)
  call timer_barrier(commsander)
  if (.not. mpi_orig .and. numtasks > 1) then 
     call xdist(x, xx(lfrctmp), nr3+iscale)
     if (isynctraj.ne.0.and.itdump) call xdist(onstepcrd, xx(lfrctmp), natom)
  endif

  ! DAC/knut change: force the coordinates to be the same on both masters.
  ! For certain compilers, addition may not be strictly commutative, so
  ! the forces on group 0 may be different by roundoff from the forces on
  ! group 1.  This can lead to divergent trajectories.  The interval at
  ! which they are resynchronized is hard-wired here to 20, which seems to
  ! work fine in our tests.
  !
  ! jwk change: coordinates are synchronized when shake is enabled above
  if (icfe .ne. 0 .and. mod(nstep+1, 20) == 0 .and. ntc == 1) then

    ! In dual-topology this is done within softcore.f
    if (ifsc .ne. 1) then
      if (master) &
        call mpi_bcast(x, nr3, MPI_DOUBLE_PRECISION, 0, commmaster, ierr)
    else
      if (master) then

        ! First, check if coordinates have desynced,
        ! then do the same for velocities
        call sc_compare(x, nr3, 'CRD')
        if (numtasks == 1) call sc_compare(v,nr3,'VEL')

        ! Resync coordinates and velocities
        call sc_sync_x(x,nr3)
      end if
    end if
    if (numtasks > 1) &
      call mpi_bcast(x, nr3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
  end if
  call timer_stop(TIME_DISTCRD)
#endif  /* MPI */

!------------------------------------------------------------------------------
  ! Fix lone pair positions
  if (numextra > 0) call local_to_global(x, xx, ix)

#ifdef MPI
  if (.not. mpi_orig .and. numtasks > 1) then
    call timer_start(TIME_DISTCRD)

!------------------------------------------------------------------------------
    ! Here we provide every processor a full copy of the velocities
    ! for removal of center of mass motion, or for archiving.
    ! (Note: this is actually over-kill: for example, only the master
    ! node really needs the velocities for archiving.  But the extra
    ! overhead of doing it this way is probably small in most cases.)
    if (ivdump .or. ivscm .or. ixdump) then
      call xdist(v, xx(lfrctmp), nr3+iscale)
      call xdist(onstepvel, xx(lfrctmp), natom)
    endif
    call timer_stop(TIME_DISTCRD)
  end if
  call timer_start(TIME_VERLET)
  ! This is the end of major broadcasting operations for parallel runs
#endif  /* MPI */
  ! }}}

!------------------------------------------------------------------------------
  ! Step 6: zero COM velocity if requested.  {{{
  ! This is used for preventing
  ! the "block of ice flying thru space" phenomenon (caused by Berendsen
  ! thermocoupling, or by coarse Ewald approximations).  It also
  ! prevents accumulations of rotational momentum in vacuum simulations.
  if (ivscm) then
    if ( nsnb /= 0 ) then
      if (mod(nstep,nsnb) == 0) ntnb = 1
    end if
    if (ifbox == 0) then
      if (ntt > 0) then

        ! Get current center of the system
        call get_position(nr, x, vcmx, vcmy, vcmz, sysrange, 0)
#ifdef MPI /* SOFT CORE */
        if (ifsc == 1) call sc_mix_position(vcmx, vcmy, vcmz, clambda)
#endif /* MPI */
        ! Center the system to the original center
        call re_position(nr, ntr, x, xc, vcmx, vcmy, vcmz, sysx, sysy, &
                         sysz, sysrange, mv_flag, 0)

      else
        ! Non-periodic simulation: remove both translation and rotation.
        ! Back the coords up 1/2 step, so that they correspond to the
        ! velocities; temporarily store in the F() array:
        f(1:nr3) = x(1:nr3) - v(1:nr3)*dt5

        ! Now compute the com motion, remove it, and recompute (just
        ! to check that it is really gone.....)
        call cenmas(nr, f, v, amass, ekcm, xcm, vcm, acm, ekrot, ocm, 4)
        call stopcm(nr, f, v, xcm, vcm, ocm, .true.)
        call cenmas(nr, f, v, amass, ekcm, xcm, vcm, acm, ekrot, ocm, 4)
      end if
    end if
    ! End of contingency for isolated, non-periodic systems (ifbox == 0)
  end if
  ! End of contingency for zeroing system center of mass velocity 

!------------------------------------------------------------------------------
  !  Also zero out the non-moving velocities if a belly is active:
  if (belly) call bellyf(nr, ix(ibellygp), v)

!------------------------------------------------------------------------------
  ! Put current velocities into vold
  vold(istart3:iend3) = v(istart3:iend3)
  ! }}}

!------------------------------------------------------------------------------
  ! Step 7: miscellaneous, get total energy {{{

  ! Pastor, Brooks, Szabo conserved quantity
  ! for harmonic oscillator: Eq. 4.7b of Mol.
  ! Phys. 65:1409-1419, 1988
  ener%kin%solv = ekpbs + ener%pot%tot
  if( ntt == 3 ) then
     ener%kin%solt = ekph  ! seems to be what Jian Li really wants for LD
  else
     ener%kin%solt = eke   ! original sander, shows energy conservation
  endif
  ener%kin%tot  = ener%kin%solt

  ! If velocities were reset, the KE is not accurate; fudge it
  ! here to keep the same total energy as on the previous step.
  ! Note that this only affects printout and averages for Etot
  ! and KE -- it has no effect on the trajectory, or on any
  ! averages of potential energy terms.
  if (resetvelo) ener%kin%tot = etot_save - ener%pot%tot

!------------------------------------------------------------------------------
  ! Total energy is sum of KE + PE:
  ener%tot = ener%kin%tot + ener%pot%tot
  etot_save = ener%tot
  ! }}}

!------------------------------------------------------------------------------
  ! Step 8: update the step counter and the integration time: {{{
  nstep = nstep + 1
  t = t + dt

  ! Full energies are only calculated every nrespa steps.
  ! nvalid is the number of steps where all energies are calculated.
  if (onstep) then
    nvalid = nvalid + 1

    ! Update all elements of these sequence types
    enert  = enert + ener
    enert2 = enert2 + (ener*ener)
    ekhf = ekhf + ekph
    ekhf2 = ekhf2 + ekph*ekph
#ifdef MPI
    if (ifsc .ne. 0) then
      sc_ener_ave(1:ti_ene_cnt) = sc_ener_ave(1:ti_ene_cnt) + &
                                  sc_ener(1:ti_ene_cnt)
      sc_ener_rms(1:ti_ene_cnt) = sc_ener_rms(1:ti_ene_cnt) + &
                                  sc_ener(1:ti_ene_cnt)**2
    end if
#endif /* MPI */
    if (nvalid == 1) etot_start = ener%tot

    kinetic_E_save(2) = kinetic_E_save(1)
    kinetic_E_save(1) = ener%kin%tot
  end if
  ! End contingency to calculate energies when on a reportable step

  ! Added for rbornstat
  ! !FIX: TL - do we need to put in rismnrespa here?
  if (mod(irespa,nrespai) == 0 .or. irespa < 2) nvalidi = nvalidi + 1
  ntnb = 0
  if ( nsnb /= 0 ) then
    if (mod(nstep,nsnb) == 0) ntnb = 1
  end if

  ! Since nstep has been incremented, total_nstep is now equal to
  ! (mdloop-1)*nstlim+nstep for REMD and nstep for MD.
  lout = .false.
  if (ntpr > 0) lout = (mod(total_nstep,ntpr) == 0 .and. onstep)
  irespa = irespa + 1

  ! }}}

!------------------------------------------------------------------------------
  ! Step 9: output from this step if required: {{{
  !    RISM dumping: {{{
  ! Some 3D-RISM files require all processes to participate in output
  ! due to the distributed memory.  RISM archive:
  if (rismprm%rism == 1) then
    ! Combined thermodynamics and distribution output.
    ! Execute if we need to do either.
    if (irismdump .or. (rism_calc_type(nstep) == RISM_FULL .and. &
        rismprm%write_thermo == 1 .and. lout)) &
      call rism_solvdist_thermo_calc(irismdump, nstep)
  end if
   ! }}}

  ! Begin writing output on the master process.
  if (master) then

    !    Restart file writing {{{
    if (ixdump) then

      ! NOTE - This assumes that if numextra > 0, then velocities are
      !        found in the array v...
      if (numextra > 0) call zero_extra_pnts_vec(v, ix)

!------------------------------------------------------------------------------
      nr = nrp
      call mdwrit(nstep, nr, ntxo, ntb, x, v, t, temp0,solvph,solve)
!------------------------------------------------------------------------------

    end if
    ! End decision process for restart file writing (ixdump flag) }}}

    !    Coordinate archive: {{{
    !   for formatted writes and replica exchange, write out a header line.
    if (itdump) then
#ifdef MPI
      ! Write out current replica#, exchange#, step#, and mytargettemp
      ! If mdloop==0 this is a normal md run (since REMD never calls
      ! corpac when mdloop==0) and we don't want the REMD header.
      ! total_nstep is set in step 5.
      if (mdloop > 0 .and. loutfm) then
        if (trxsgld) then
          write (MDCRD_UNIT,'(a,4(1x,i8))') "RXSGLD ", repnum, mdloop, &
                total_nstep, stagid
        else
          write (MDCRD_UNIT,'(a,3(1x,i8),1x,f8.3)') "REMD ", repnum, mdloop, &
                total_nstep, my_remd_data%mytargettemp
        end if
      end if
#endif /* MPI */
      ! Dump coordinates:
      if (isynctraj.ne.0) then
        call corpac(onstepcrd, 1, nrx, MDCRD_UNIT, loutfm)
      else
        call corpac(x,1,nrx,MDCRD_UNIT,loutfm)
      endif
      if (ntb > 0) call corpac(box, 1, 3, MDCRD_UNIT, loutfm)

    end if
    ! End contingency for trajectory writing on the master process (itdump) }}}

    !    Velocity archive: {{{
    if (ivdump) then

      ! NOTE - This assumes that if numextra > 0, then velocities are
      !        found in the array v...
      if (numextra > 0) call zero_extra_pnts_vec(v, ix)
#ifdef MPI
      ! Write out current replica#, exchange#, step#, and mytargettemp
      ! If mdloop==0 this is a normal md run (since REMD never calls corpac
      ! when mdloop==0) and we don't want the REMD header.
      if (mdloop > 0 .and. loutfm) then
        if (trxsgld) then
          write (MDVEL_UNIT,'(a,4(1x,i8))') "RXSGLD ", repnum, mdloop, &
                total_nstep, stagid
        else
          write (MDVEL_UNIT,'(a,3(1x,i8),1x,f8.3)') "REMD ", repnum, mdloop, &
                total_nstep, my_remd_data%mytargettemp
        end if
      end if
#endif /* MPI */

      if (isynctraj.eq.0) then
         call corpac(v, 1, nrx, MDVEL_UNIT, loutfm)
      else
         call corpac(onstepvel, 1, nrx, MDVEL_UNIT, loutfm)
      endif
    end if
    ! }}}

!------------------------------------------------------------------------------
    ! Energy archive: (total_nstep set in Step 5.)
    if ( ntwe > 0 ) then
      if (mod(total_nstep,ntwe) == 0 .and. onstep) &
        call mdeng(15,nstep,t,ener,onefac,ntp)
    end if

    if (ioutfm > 0) then
      if (itdump) call end_binary_frame(MDCRD_UNIT)
      if (ivdump .and. ntwv > 0) call end_binary_frame(MDVEL_UNIT)
      if (ifdump .and. ntwf > 0) call end_binary_frame(MDFRC_UNIT)
    end if

!------------------------------------------------------------------------------
    !    General printed output:  {{{
    if (lout) then
      if (facc .ne. 'A') rewind(7)

      call prntmd(total_nstep, t, ener, onefac, 7, .false.)

#ifdef MPI /* SOFT CORE */
      if (ifsc .ne. 0) call sc_print_energies(6, sc_ener)
      if (ifsc .ne. 0) call sc_print_energies(7, sc_ener)
#endif

!------------------------------------------------------------------------------
      if (nmropt > 0) call nmrptx(6)
      if (infe > 0) call nfe_prt(6)
      if (itgtmd == 2) then
        emtmd = 0.0d0
        call mtmdcall(emtmd, xx(lmtmd01), ix(imtmd02), x, f, ih(m04), &
                      ih(m02), ix(i02), ih(m06), xx(lmass), natom, &
                      nres, 'PRNT')
      end if
      call flush(7)
    end if
    ! end of giant "if (lout)" contingency related to data output }}}

!------------------------------------------------------------------------------
    !    Output running averages: {{{
    ! total_nstep = Total nstep REMD/MD, set in step 5
    if (ntave > 0) then
      if (mod(total_nstep,ntave) == 0 .and. onstep) then
        write(6, 542)
        if (rismprm%rism == 1) then
          tspan = ntave / mylcm(nrespa, rismprm%rismnrespa)
        else
          tspan = ntave / nrespa
        end if

        ! Update all elements of these sequence types
        enert_tmp  = enert - enert_old
        enert2_tmp = enert2 - enert2_old
        enert_old  = enert
        enert2_old = enert2
        enert_tmp  = enert_tmp/tspan
        enert2_tmp = enert2_tmp/tspan - enert_tmp*enert_tmp
        call zero_neg_values_state(enert2_tmp)
        enert2_tmp = sqrt(enert2_tmp)
#ifdef MPI
        if (ifsc .ne. 0) then
          do m = 1,ti_ene_cnt
            sc_ener_tmp(m) = sc_ener_ave(m) - sc_ener_old(m)
            sc_ener_tmp2(m) = sc_ener_rms(m) - sc_ener_old2(m)
            sc_ener_old(m) = sc_ener_ave(m)
            sc_ener_old2(m) = sc_ener_rms(m)
            sc_ener_tmp(m) = sc_ener_tmp(m) / tspan
            sc_ener_tmp2(m) = sc_ener_tmp2(m)/tspan - sc_ener_tmp(m)**2
            if (sc_ener_tmp2(m) < 0.0d0) then
              sc_ener_tmp2(m) = 0.0d0
            end if
            sc_ener_tmp2(m) = sqrt(sc_ener_tmp2(m))
          end do
        end if
#endif
        if (rismprm%rism == 1) then
          write(6, 540) ntave / mylcm(nrespa, rismprm%rismnrespa)
        else
          write(6, 540) ntave/nrespa
        end if
        call prntmd(total_nstep, t, enert_tmp, onefac, 0, .false.)
#ifdef MPI
        if (ifsc .ne. 0) call sc_print_energies(6, sc_ener_tmp)
#endif /* MPI */
        write(6, 550)
        call prntmd(total_nstep, t, enert2_tmp, onefac, 0, .true.)
#ifdef MPI /* SOFT CORE */
        if (ifsc .ne. 0) call sc_print_energies(6, sc_ener_tmp2)
#endif /* MPI */
        if (icfe > 0) then
          if (rismprm%rism == 1) then
            write (6, 541) ntave / mylcm(nrespa, rismprm%rismnrespa)
          else
            write (6, 541) ntave/nrespa
          end if
          edvdl_r = edvdl_r/tspan
          edvdl_r%pot%dvdl = enert_tmp%pot%dvdl  ! fix for DV/DL output
          edvdl_r%virvsene = 0.d0 ! virvsene should not but included here
          call prntmd(total_nstep, t, edvdl_r, onefac, 0, .false.)
          edvdl_r = null_state_rec
        end if
        write(6,542)
      end if
    end if
    ! End contingency to output running averages (ntave)      }}}
  end if
  ! End output work on the master process

  ! Soft core output:  {{{
#ifdef MPI /* SOFT CORE */
  if (ntave > 0 .and. icfe > 0 .and. dynlmb > 0) then
    if (mod(nstep,ntave) == 0 .and. onstep) then

!------------------------------------------------------------------------------
      ! For runs with dynamically changing lambda, raise lambda here
      ! and flush all buffers for the next averages
      clambda = clambda + dynlmb
      call sc_change_clambda(clambda)
      if (master) then
        sc_ener(1:ti_ene_cnt) = 0.0d0
        sc_ener_ave(1:ti_ene_cnt) = 0.0d0
        sc_ener_rms(1:ti_ene_cnt) = 0.0d0
        sc_ener_old(1:ti_ene_cnt) = 0.0d0
        sc_ener_old2(1:ti_ene_cnt) = 0.0d0
        enert = null_state_rec
        enert2 = null_state_rec
        enert_old = null_state_rec
        enert2_old = null_state_rec
        write (6, *)
        write (6, '(a,f12.4,a,f12.4)') &
              'Dynamically changing lambda: Increased clambda by ', &
              dynlmb, ' to ', clambda
        write (6,*)
      end if
    end if
  end if
#endif
  ! }}}

  ! }}}
!------------------------------------------------------------------------------
  call timer_stop(TIME_VERLET)
  ! Miscellaneous stuff at the end of each step: {{{
#if !defined(DISABLE_NFE) && defined(NFE_ENABLE_BBMD)
  if (infe == 1) then
    call xdist(x, xx(lfrctmp), nr3+iscale)
    call xdist(v, xx(lfrctmp), nr3+iscale)
    call nfe_on_mdstep(ener%pot%tot, x, v, ekmh)
    vold(istart3:iend3) = v(istart3:iend3)
  end if
#endif /* DISABLE_NFE is NOT defined, but NFE_ENABLE_BBMD is */

  if (plumed .ne. 0 .and. plumed_stopflag .ne. 0) goto 480
  ! }}}

  ! This is where we actually cycle back to the next MD step
  if (nstep < nstlim) goto 260

  480 continue

#ifdef MPI
!------------------------------------------------------------------------------
  ! Replica Exchange MD post-dynamics work {{{
  if (next_rem_method == 1) then
    remd_ekmh = ekmh

    ! Hybrid REMD
    if (numwatkeep >= 0) then
      ! This is a hybrid REMD run. Get energy of
      ! stripped system for next exchange.
      call hybrid_remd_ene(xx, ix, ih, ipairs, qsetup, numwatkeep, hybridgb, &
                           igb, ntr, nspm, t, temp0, ntb, cut, ener, &
                           do_list_update, nstep, onefac)
    else

      ! The positions are currently one step ahead of the energy ener%pot%tot,
      ! since force was called prior to the position propagation. Thus, call
      ! force one more time to update ener%pot%tot to reflect the current
      ! coordinates.
      call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, xx(l96), xx(l97), &
                 xx(l98), xx(l99), qsetup, do_list_update, nstep)
    end if
    ! End branch for additional energy computations in Hybrid REMD

    ! Set myeptot, mytemp, and mytargettemp
    my_remd_data%mytemp = ener%kin%tot * onefac(1)
    my_remd_data%myeptot = ener%pot%tot
    my_pressure = pres0
    my_volume = ener%volume

    my_remd_data%mytargettemp = temp0
#  ifdef VERBOSE_REMD
    if (master) write(6, '(a,f15.4,2(a,f6.2))') "| REMD: myEptot= ", &
           my_remd_data%myeptot, " myTargetTemp= ", &
           my_remd_data%mytargettemp, " mytemp= ", my_remd_data%mytemp
#  endif /* VERBOSE_REMD */
  else if (next_rem_method == 3) then
    remd_ekmh = ekmh
    if (mdloop > 0) then
      my_remd_data%mytemp = ener%kin%tot * onefac(1)
    end if
    my_remd_data%mytargettemp = temp0

    ! This call to force will bring all energies up-to-date
    call force(xx, ix, ih, ipairs, x, f, ener, ener%vir, xx(l96), xx(l97), &
               xx(l98), xx(l99), qsetup, do_list_update, nstep)
    my_remd_data%myeptot = ener%pot%tot
    my_pressure = pres0
    my_volume = ener%volume
    ! Call nmrdcp to decrement the NMR counter, since this should not count as
    ! a real step (JMS 2/12). This is OK, since the counter got incremented at
    ! the very end of nmrcal, so we haven't already printed an unwanted value.
    if (nmropt /= 0) call nmrdcp

    ! Call xdist such that master has all the velocities
    call xdist(v, xx(lfrctmp), 3*natom+iscale)
  else if (next_rem_method == 4 .or. next_rem_method == 5) then
    remd_ekmh = ekmh
  endif
#endif /* MPI */
  ! End Replica Exchange MD post-dynamics work }}}

  ! Stochastic Isokinetic Nose-Hoover RESPA (SINR) integrator clean-up
  if (ntt==10) call sinr_cleanup(sinrdata)

  ! Print averages {{{
#ifdef MPI
  ! Turn off avg. for REMD. and explicit solvent CpHMD, since it's not
  ! accumulated correctly in that case for each compiler
  if (master .and. rem == 0) then
#else
  if (master) then
#endif /*MPI*/
    tspan = nvalid
    if (nvalid > 0) then

      ! Update all elements of these sequence types
      enert  = enert/tspan
      enert2 = enert2/tspan - enert*enert
      call zero_neg_values_state(enert2)
      enert2 = sqrt(enert2)
      edvdl = edvdl/tspan

      if (rem == 0 .and. master) then
         ekhf = ekhf / tspan
         ekhf2 = ekhf2/tspan - ekhf*ekhf
         ekhf2 = sqrt(ekhf2)
      endif
#ifdef MPI
      if (ifsc .ne. 0) then
        do m = 1, ti_ene_cnt
          sc_ener_ave(m) = sc_ener_ave(m)/tspan
          sc_ener_rms(m) = sc_ener_rms(m)/tspan - sc_ener_ave(m)**2
          if (sc_ener_rms(m) < 0.0d0) then
            sc_ener_rms(m) = 0.0d0
          end if
          sc_ener_rms(m) = sqrt(sc_ener_rms(m))
        end do
      end if
#endif
      write(6, 540) nvalid
      call prntmd(total_nstep, t, enert, onefac, 0, .false.)
#ifdef MPI /* SOFT CORE */
      if (ifsc .ne. 0) call sc_print_energies(6, sc_ener_ave)
#endif
      if (nmropt > 0) call nmrptx(6)
      write(6, 550)
      call prntmd(total_nstep, t, enert2, onefac, 0, .true.)
#ifdef MPI
      if (ifsc .ne. 0) call sc_print_energies(6, sc_ener_rms)
      if (ifsc .ne. 0) call sc_print_dvdl_values()
      if (icfe > 0) then
        write(6, 541) nvalid
        edvdl%pot%dvdl = enert%pot%dvdl  ! fix for DV/DL output
        edvdl%virvsene = 0.d0 ! virvsene should not but included here
        call prntmd(total_nstep, t, edvdl, onefac, 0, .false.)
      end if
#endif /* MPI */
      if (nmropt >= 1) then
        write(6, 500)
        if (iredir(7) .ne. 0) then
          call pcshift(-1, x, f)
        end if
        
        call ndvptx(x, f, ih(m04), ih(m02), ix(i02), nres, iscrth, &
                    natom, xx(lwinv), xx(lnmr01), ix(inmr02), 6)
      end if

!------------------------------------------------------------------------------
      ! Print Born radii statistics
      if (rbornstat == 1 .and. (igb .ne. 0 .or. ipb .ne. 0)) then

        ! Born radii stats collected every nrespai step not nrespa step
        tspan = nvalidi
        write(6, 580) nstep
        write(6, 590)
        do m = 1, natom
          xx(l188-1+m) = xx(l188-1+m)/tspan
          xx(l189-1+m) = xx(l189-1+m)/tspan - xx(l188-1+m)*xx(l188-1+m)
          xx(l189-1+m) = sqrt(xx(l189-1+m))
          write(6, 600) m, xx(l186-1+m), xx(l187-1+m), xx(l188-1+m), &
                        xx(l189-1+m)
        end do
      end if

      enert%kin%tot = enert%kin%tot*onefac(1)
      enert2%kin%tot = enert2%kin%tot*onefac(1)
      enert%kin%solt = enert%kin%solt*onefac(2)
      enert2%kin%solt = enert2%kin%solt*onefac(2)
      enert%kin%solv = enert%kin%solv*onefac(3)
      enert2%kin%solv = enert2%kin%solv*onefac(3)
      temp = enert%kin%tot
    end if
    ! End contingency for nvalid > 0, signifying that
    ! all energies must be calculated
    if (ntp > 0 .and. barostat == 2) call mcbar_summary
  end if
  ! End of contingency for work on the master process;
  ! this started about 120 lines ago and the way it started
  ! depends on whether MPI is part of the compilation. }}}

  ! deallocates: {{{
  if( ntc>1 ) then
    deallocate( xold, stat=ierr )
    REQUIRE( ierr == 0 )
  end if
  if (icfe .ne. 0) then
    deallocate( frcti, stat = ierr )
    REQUIRE( ierr == 0 )
  end if
  if (plumed .ne. 0) call plumed_f_gfinalize()
  ! }}}

  ! format statements: {{{

  500 format(/,' NMR restraints on final step:'/)
  540 format(/5x,' A V E R A G E S   O V E R ',i7,' S T E P S',/)
  541 format(/5x,' DV/DL, AVERAGES OVER ',i7,' STEPS',/)
  542 format('|',79('='))
  550 format(/5x,' R M S  F L U C T U A T I O N S',/)
  580 format('STATISTICS OF EFFECTIVE BORN RADII OVER ',i7,' STEPS')
  590 format('ATOMNUM     MAX RAD     MIN RAD     AVE RAD     FLUCT')
  600 format(i4,2x,4f12.4)
  ! }}}
  return
end subroutine runmd

subroutine initialize_runmd(x,ix,v)

  implicit none
  integer, intent(in) :: ix(*)
  _REAL_, intent(in) :: x(*), v(*)

  ! Initialize some variables {{{
#ifdef MPI
  if (master) then
    ! In Replica Exchange Molecular Dynamics (REMD), runmd will be called many
    ! times, so we dont want to open files every time.  For normal md, mdloop
    ! will just be 0.
    if (mdloop .eq. 0) then
      call amopen(7, mdinfo, 'U', 'F', facc)
    end if
  end if

#else
  call amopen(7, mdinfo, 'U', 'F', 'W')
#endif

  vlim = (vlimit > small)
  ntcmt = 0
  izero = 0
  belly = (ibelly > 0)
  lout = .true.
  loutfm = (ioutfm <= 0)
  nr = nrp
  nr3 = 3*nr
  ekmh = 0.d0

  etot_save = 0.d0
  pres0x = 0.d0
  pres0y = 0.d0
  pres0z = 0.d0
  gamma_ten_int = 0.d0
  dttp = 0.d0
  ekph = 0.d0
  ekpbs = 0.d0
  eke = 0.d0

  do_list_update = .false.
#ifdef MPI
  if (mpi_orig) then
    istart = 1
    iend = natom
  else
    istart = iparpt(mytaskid) + 1
    iend = iparpt(mytaskid+1)
  end if
#else
  istart = 1
  iend = nr
#endif
  istart3 = 3*istart -2
  iend3 = 3*iend
  if( mytaskid == numtasks -1 )  iend3 = iend3 + iscale

  if( ntc>1 ) allocate (xold(3*natom+iscale))
#ifdef MPI
  if (icfe /= 0) then
    allocate(frcti(nr3 + 3*extra_atoms), stat=ierr)
    REQUIRE(ierr == 0)
  end if
#endif

  ! If ntwprt.NE.0, only print the atoms up to this value
  nrx = nr3
  if (ntwprt > 0) nrx = ntwprt*3

  ! Cleanup the velocity if belly run
  if (belly) call bellyf(nr,ix(ibellygp),v)
  ! Determine system degrees of freedom (for T scaling, reporting)
#   include "degcnt.inc"

!   }}}

!------------------------------------------------------------------------------
!    Pressure/temp units  {{{
  ! Begin unit conversion.  pconv eventually becomes a factor to convert
  ! pressure in kcal/mole-A^3 to bar.
  boltz2 = 8.31441d-3 * 0.5d0
  pconv = 1.6604345d+04
  boltz2 = boltz2/4.184d0
  dtx = dt*20.455d+00
  dtxinv = 1.0d0 / dtx
  dt5 = dtx * 0.5d0
  pconv = pconv*4.184d0

  ! fac() are #deg freedom * kboltz / 2.  Multiply by T to get the expected
  ! kinetic energy.  fac(1) is for the entire system.
  fac(1) = boltz2*rndf
  fac(2) = boltz2*rndfp
  if (rndfp < 0.1d0) fac(2) = 1.d-6

  fac(3) = boltz2*rndfs
  if (rndfs < 0.1d0) fac(3) = 1.d-6
  onefac(1) = 1.0d0 / fac(1)
  onefac(2) = 1.0d0 / fac(2)
  onefac(3) = 1.0d0 / fac(3)
  factt = rndf/(rndf+ndfmin)

  ! These are "desired" kinetic energies based on the number of
  ! degrees of freedom and target temperature.  They will be used
  ! for calculating the velocity scaling factor
  ekinp0 = fac(2)*temp0
  ekins0 = fac(3)*temp0
  ekin0  = fac(1)*temp0
  ! }}}

!------------------------------------------------------------------------------
  !    Langevin dynamics setup  {{{
  if (nscm > 0) then
     if (ifbox == 0) then
        call get_position(nr, x, sysx, sysy, sysz, sysrange, 0)
#ifdef MPI
        ! Soft core position mixing
        if (ifsc == 1) call sc_mix_position(sysx, sysy, sysz, clambda)
#endif /* MPI */
     end if ! ifbox==0: non-periodic
  end if ! nscm is enabled
  ! }}}

  !   General initialization:  {{{
  nrek = 4
  nrep = 15

  nvalid = 0
  nvalidi = 0
  nstep = 0
  total_nstep = 0
#ifdef MPI
  ! For REMD, total_nstep is the number of steps * the number of exchanges
  ! we've already attempted
  if (rem .ne. 0) total_nstep = (mdloop - 1) * nstlim
#endif /* MPI */
  fit = 0.d0
  fiti = 0.d0
  fit2 = 0.d0

!------------------------------------------------------------------------------
  ! Zero all elements of these sequence types
  ener       = null_state_rec
  enert      = null_state_rec
  enert2     = null_state_rec
  enert_old  = null_state_rec
  enert2_old = null_state_rec
  edvdl      = null_state_rec
  edvdl_r    = null_state_rec

!------------------------------------------------------------------------------

  ener%kin%pres_scale_solt = 1.d0
  ener%kin%pres_scale_solv = 1.d0
  ener%box(1:3) = box(1:3)
  ener%cmt(1:4) = 0.d0
  nitp = 0
  nits = 0
  ekhf = 0.0d0
  ekhf2 = 0.0d0

!------------------------------------------------------------------------------
#ifdef PLUMED
  ! PLUMED initialization.  PLUMED is an open-source plugin that
  ! confers the functionality of a number of enhanced sampling methods.
  if (plumed == 1) then
#   include "Plumed_init.inc"
  endif
#endif
  ! }}}

end subroutine initialize_runmd

#ifdef MPI
subroutine thermodynamic_integration(f)

   implicit none
   _REAL_, intent(inout) :: f(nr3)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  ! If softcore potentials are used, collect their dvdl contributions:
  if (ifsc .ne. 0) then
    call mpi_reduce(sc_dvdl, sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, commsander, ierr)
    sc_dvdl=0.0d0 ! zero for next step
    call mpi_reduce(sc_dvdl_ee, sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, commsander, ierr)
    sc_dvdl_ee = 0.0d0 ! zero for next step
    call mpi_reduce(sc_ener, sc_ener_tmp, ti_ene_cnt, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, 0, commsander, ierr)
    sc_ener(1:ti_ene_cnt) = sc_ener_tmp(1:ti_ene_cnt)
  end if
  if (ifsc == 2) then

    ! If this is a perturb to nothing run, scale forces and calculate dvdl
    call sc_nomix_frc(f,nr3,ener)
    if (numtasks > 1) then
      call mpi_bcast(f, nr3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
      call mpi_bcast(ener, state_rec_len, MPI_DOUBLE_PRECISION, 0, &
                     commsander, ierr)
    end if
  end if

!------------------------------------------------------------------------------
  ! Multi-state Bennet Acceptance Ratio upkeep
  if (ifmbar .ne. 0 .and. do_mbar) call bar_collect_cont()

!------------------------------------------------------------------------------
  if (icfe .ne. 0) then
  ! Free energies using thermodynamic integration:

    ! First, partners exchange forces and energies:
    if (master) then
      partner = ieor(masterrank, 1)
      call mpi_sendrecv(f, nr3, MPI_DOUBLE_PRECISION, partner, 5, frcti, &
                        nr3 + 3*extra_atoms, MPI_DOUBLE_PRECISION, partner, &
                        5, commmaster, ist, ierr )
      call mpi_sendrecv(ener, state_rec_len, MPI_DOUBLE_PRECISION, partner, &
                        5, ecopy, state_rec_len, MPI_DOUBLE_PRECISION, &
                        partner, 5, commmaster, ist, ierr)

      ! Exchange sc-dvdl contributions between masters:
      call mpi_sendrecv(sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                        sc_tot_dvdl_partner, 1, MPI_DOUBLE_PRECISION, &
                        partner, 5, commmaster, ist, ierr)
      call mpi_sendrecv(sc_tot_dvdl_ee, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                        sc_tot_dvdl_partner_ee, 1, MPI_DOUBLE_PRECISION, &
                        partner, 5, commmaster, ist, ierr )

      ! Collect statistics for free energy calculations
      if (onstep) then
        if (masterrank == 0) then
          if (klambda == 1) then
            edvdl = edvdl - ener + ecopy
            edvdl_r = edvdl_r - ener + ecopy
          else
            clfac = klambda*(1.d0 - clambda)**(klambda-1)
            edvdl = edvdl - (ener - ecopy)*clfac
            edvdl_r = edvdl_r - (ener - ecopy)*clfac
          end if
        else
          if (klambda == 1) then
            edvdl = edvdl + ener - ecopy
            edvdl_r = edvdl_r + ener - ecopy
          else
            clfac = klambda*(1.d0 - clambda)**(klambda-1)
            edvdl = edvdl + (ener - ecopy)*clfac
            edvdl_r = edvdl_r + (ener - ecopy)*clfac
          end if
        end if

        ! This includes the sc-dvdl contribution into the vdw-part
        ! and potential energy parts of the dvdl-statistics
        if (ifsc == 1) call adj_dvdl_stat(edvdl, edvdl_r)
      end if

      ! Do energy collection for MBAR FEP runs
      if (ifmbar .ne. 0 .and. do_mbar) &
        call calc_mbar_energies(ener%pot%tot, ecopy%pot%tot)

      if (masterrank == 0) then
        call mix_frcti(frcti, ecopy, f, ener, nr3, clambda, klambda)
      else
        call mix_frcti(f, ener, frcti, ecopy, nr3, clambda, klambda)
      endif
    endif

    if (numtasks > 1) then
      call mpi_bcast(f, nr3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
      call mpi_bcast(ener, state_rec_len, MPI_DOUBLE_PRECISION, 0, &
                     commsander, ierr)
    end if
  ! End contingency for free energies by Thermodynamic Integration }}}
  end if

end subroutine thermodynamic_integration
#endif

subroutine modwt_reset()

   implicit none
  ! Reset quantities depending on TEMP0, which may have been changed
  ! by  MODWT during FORCE call.
  ekinp0 = fac(2)*temp0

  ekins0 = fac(3)*temp0
  ekin0 = fac(1)*temp0

end subroutine modwt_reset

end module runmd_module
