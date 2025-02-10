! <compile=optimized>
#include "../include/dprec.fh"
#include "../include/assert.fh"
#include "nfe-config.h"

!------------------------------------------------------------------------------
! force: The main driver routine to compute energies and forces in sander.
!        This is a simple nexus of many different modules, each with their own
!        contributions to how the atoms will move next.
!
! Arguments:
!   xx:             global real array (holding e.g. coordinates at position
!                   l15 etc., see locmem.F90)
!   ix:             global array holding all integer arrays (see locmem.F90)
!   ih:             Hollerith array containing atom names, types,
!                   residue names, and more
!   ipairs:         ?? Global pairlist ?? -- needs explanation according to JMS
!   x:              Coordinates of all atoms
!   f:              Forces on all atoms
!   ener:           Energy with components (see the state_rec type in the
!                   state.F90 module)
!   vir:            Virial (four element _REAL_ vector)
!   fs:             
!   rborn:          The Generalized Born radii (they will be recalculated
!                   within this subroutine)
!   reff:           The effective Born radii (recalculated in this subroutine)
!   onereff:        The inverse effective Born radii, 1.0/reff 
!   qsetup:         Flag to activate setup of multiple components, .false. on
!                   first call
!   do_list_update: flag to have the non-bonded list updated (returns .TRUE. or
!                   .FALSE. to the calling subroutine following a call to
!                   nonbond_list)
!   nstep:          The step number
!------------------------------------------------------------------------------
subroutine force(xx, ix, ih, ipairs, x, f, ener, vir, fs, rborn, reff, &
                 onereff, qsetup, do_list_update, nstep)

#if !defined(DISABLE_NFE)
  use nfe_sander_hooks, only: nfe_on_force => on_force
  use nfe_sander_proxy, only: infe
#endif /* DISABLE_NFE */
  use file_io_dat
#ifdef LES
  use genbornles
#  ifdef MPI
  use les_data, only: elesa, elesb, elesd, elesp
#  endif /* MPI */
#else
  use genborn
#endif /* LES */
#ifdef MPI
  use softcore, only: sc_ener
#endif /* MPI */
#if defined(LES) && defined(MPI)
  use remd, only: rem ! wasn't used for LES above
#endif /* LES && MPI */
  use sander_rism_interface, only: rismprm, rism_force
  use stack
  use constants, only: zero, one
  use relax_mat
  use ew_recip
  use parms, only: cn1, cn2, cn6, asol, bsol, pk, rk, tk, numbnd, numang, &
                   nptra, nphb, nimprp, cn3, cn4, cn5 ! for another vdw model
  use nblist, only: nonbond_list, a, b, c, alpha, beta, gamma
#ifdef DSSP
  use dssp, only: fdssp, edssp, idssp
#endif /* DSSP */

  use emap, only: temap, emapforce

  use linear_response, only: ilrt, ee_linear_response, energy_m0, energy_w0, &
                             energy_vdw0, cn1_lrt, cn2_lrt, crg_m0, crg_w0, &
                             do_lrt, f_scratch, lrt_solute_sasa
  use xray_interface_module, only: xray_get_derivative
  use xray_globals_module, only: atom_bfactor, xray_energy, xray_active

  ! CHARMM Force Field Support
  use charmm_mod, only: charmm_active, charmm_calc_impropers, &
                        charmm_calc_cmap, charmm_calc_urey_bradley,&
                        charmm_dump_gold, &
                        do_charmm_dump_gold
  use ff11_mod, only: cmap_active, calc_cmap
  use state
  use les_data, only: temp0les
  use music_module, only : music_force
#ifdef MPI
   use mpi
#endif

  implicit none
  
  integer, intent(in) :: nstep

#ifdef MPI
  integer ierr
#endif
  integer   ipairs(*)
  _REAL_ xx(*)
  integer   ix(*)
  character(len=4) ih(*)
  _REAL_ fs(*), rborn(*), reff(*), dvdl, xray_e
  _REAL_, intent(out) :: onereff(*)
#include "def_time.h"
#include "ew_frc.h"
#include "ew_cntrl.h"
#include "extra_pts.h"
#include "parallel.h"

#ifdef MPI
#  include "ew_parallel.h"
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
  integer gb_rad_mpistart, j3, j, i3
  _REAL_ :: temp_amd_totdih
#endif /* MPI */

  logical belly
#include "../include/md.h"
#include "box.h"
#include "nmr.h"
#include "../include/memory.h"
#include "extra.h"
#include "tgtmd.h"
#include "multitmd.h"
#include "flocntrl.h"
  integer istart, iend
  _REAL_ evdwex, eelex
  _REAL_ enemap

  logical, intent(inout) :: qsetup
  logical, intent(out) :: do_list_update

  _REAL_  enmr(6), devdis(4), devang(4), devtor(4), devpln(4), devplpt(4), &
          devgendis(4), entr, ecap, enfe
  _REAL_, target, intent(in) :: x(3*natom)
  _REAL_, target, intent(out) :: f(3*natom)
  _REAL_  vir(4)
  type(state_rec)  ener

  ! Local
  _REAL_                     :: ene(30)    !Used locally ONLY
  type(potential_energy_rec) :: pot        !Used locally ONLY
  logical, save :: first=.true.

  _REAL_, pointer :: x3(:,:)     ! to pass to xray_get_derivative()
  _REAL_, pointer :: f3(:,:)     ! to pass to xray_get_derivative()

#ifndef LES
  _REAL_ escf
#endif /* LES */

  integer i
  _REAL_  virvsene, eelt, epol, esurf, edisp
  _REAL_ erism

  ! Charge transfer
  _REAL_ ect

  _REAL_ epolar, aveper, aveind, avetot, emtot, dipiter, dipole_temp
  integer, save :: newbalance
  integer, save :: xray_nstep = 0
   
  ! Aceelerated MD variables
  _REAL_ amd_totdih

  ! MuSiC
  _REAL_ :: music_vdisp, music_vang, music_vgauss, music_spohr89

  x3(1:3,1:natom) => x(1:3*natom)   ! for xray_get_derivative
  f3(1:3,1:natom) => f(1:3*natom)   ! for xray_get_derivative

  ect = 0.0


  call timer_start(TIME_FORCE)
  ene(:) = ZERO 
  call zero_pot_energy(pot)

  belly = ibelly > 0

#ifdef MPI
  if (mpi_orig) then

    ! Check to see if we are done yet in mpi_orig case (tec3).
    ! This is done by monitoring the status of an integer notdone.
    ! If notdone .eq. 1 then we keep going.  notdone is set to zero
    ! when we no longer want to call force().  This perhaps is not the
    ! most efficient means to implement this check...
    call mpi_bcast(notdone,1,mpi_integer,0,commsander,ierr)
    if (notdone /= 1) return

    ! Send copies of xyz coords, setbox common block, vir array
    ! and NTNB value to all nodes from master with a broadcast.
    if (numtasks > 1) then
      call mpi_bcast(box, BC_BOXR, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
      call mpi_bcast(ntb, BC_BOXI, mpi_integer, 0, commsander, ierr)
      call mpi_bcast(vir, 3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
      call mpi_bcast(xx(lcrd), 3*natom, MPI_DOUBLE_PRECISION, 0, &
                     commsander, ierr)
      call mpi_bcast(ntnb, 1, mpi_integer, 0, commsander, ierr)
    end if
  end if
  istart = iparpt(mytaskid) + 1
  iend = iparpt(mytaskid+1)
#else
  istart = 1
  iend = natom
#endif /* MPI */

  ! Zero out the energies and forces
  enoe = 0.d0
  aveper = 0.d0
  aveind = 0.d0
  avetot = 0.d0
  dipiter = 0.d0
  dvdl = 0.d0
  dipole_temp = 0.d0
  enmr(1:6) = 0.d0
  enfe = 0.d0
  vir(1:4) = 0.d0
  virvsene = 0.d0
  f(1:3*natom+iscale) = 0.d0

  if (igb == 0 .and. ipb == 0 .and. iyammp == 0) then

    ! For GB: do all nonbondeds together below
    call timer_start(TIME_NONBON)
    call timer_start(TIME_LIST)
    call nonbond_list(x, ix(i04), ix(i06), ix(i08), ix(i10), ntypes, &
                      natom, xx, ix, ipairs, ntnb, ix(ibellygp), &
                      belly, newbalance, qsetup, do_list_update)
    call timer_stop(TIME_LIST)
    call timer_stop(TIME_NONBON)
  end if

#if !defined(DISABLE_NFE)
  if (infe == 1) then
    call nfe_on_force(x, f, enfe)
  end if
#endif

  ! Do weight changes, if requested.  Updated 9/2007 by Matthew Seetin
  ! to enable plane-point and plane-plane restraints.
  if (nmropt > 0) then
    call nmrcal(x, f, ih(m04), ih(m02), ix(i02), xx(lwinv), enmr, devdis, &
                devang, devtor, devplpt, devpln, devgendis, temp0, tautp, &
                cut, xx(lnmr01), ix(inmr02), xx(l95), 31, 6, rk, tk, pk, cn1, &
                cn2, asol, bsol, xx(l15), numbnd, numang, nptra-nimprp, &
                nimprp, nphb, natom, natom, ntypes, nres, rad, wel, radhb, &
                welhb, rwell, tgtrmsd, temp0les, -1, 'WEIT')
  end if
  epolar = 0.d0

  ! EGB: if Generalized Born is in effect (there is a GB solvent,
  ! not just electrostatics in a vacuum or some Poisson-Boltzmann
  ! solvent), then we need to calculate the GB radii for this
  ! structure.  If we are doing GB in the context of QM (that is,
  ! qm_gb = 2) then we need to calculate the GB radii before
  ! calling qm_mm.
  if (igb > 0 .and. igb /= 6 .and. igb /= 10 .and. ipb == 0 .and. &
      (irespa < 2 .or. mod(irespa,nrespai) == 0)) then
#ifdef MPI
    gb_rad_mpistart = mytaskid+1
#endif
    call timer_start(TIME_EGB)
    call timer_start(TIME_GBRAD1)

#ifdef LES
#  ifdef MPI
    call egb_calc_radii(igb, natom, x, fs, reff, onereff, fsmax, rgbmax, &
                        rborn, offset, rbornstat, xx(l188), xx(l189), &
                        xx(l186), xx(l187), gbneckscale, ncopy, rdt, &
                        gbalpha, gbbeta, gbgamma, gb_rad_mpistart)
#  else
    call egb_calc_radii(igb, natom, x, fs, reff, onereff, fsmax, rgbmax, &
                        rborn, offset, rbornstat, xx(l188), xx(l189), &
                        xx(l186), xx(l187), gbneckscale, ncopy, rdt, &
                        gbalpha, gbbeta, gbgamma)
#  endif
#else
#  ifdef MPI
    call egb_calc_radii(igb, natom, x, fs, reff, onereff, fsmax, rgbmax, &
                        rborn, offset, rbornstat, xx(l188), xx(l189), &
                        xx(l186), xx(l187), gbneckscale, ncopy, rdt, &
                        xx(l2402), xx(l2403), xx(l2404), gb_rad_mpistart)
#  else
    call egb_calc_radii(igb, natom, x, fs, reff, onereff, fsmax, rgbmax, &
                        rborn, offset, rbornstat, xx(l188), xx(l189), &
                        xx(l186), xx(l187), gbneckscale, ncopy, rdt, &
                        xx(l2402), xx(l2403), xx(l2404))
#  endif
#endif
    call timer_stop(TIME_GBRAD1)
    call timer_stop(TIME_EGB)
  end if
  ! End EGB
 
  ! Calculate the non-bonded contributions
  call timer_start(TIME_NONBON)

  if (igb == 0 .and. ipb == 0 .and. iyammp == 0) then

    ! (for GB: do all nonbondeds together below)
    call timer_start(TIME_EWALD)
    if (ilrt /= 0) then

      ! Modifications for computing interaction energy
      ! according to the Linear Response Theory, LIE module
      if (do_lrt) then

        ! call with molecule charges set to zero
        call ewald_force(x, natom, ix(i04), ix(i06), crg_m0, cn1, cn2, &
                         cn6, energy_m0, epolar, f_scratch, xx, ix, &
                         ipairs, virvsene, xx(lpol), &
                         xx(lpol2), .false. , cn3, cn4, cn5)

        ! call with water charges set to zero
        call ewald_force(x, natom, ix(i04), ix(i06), crg_w0, cn1, cn2, &
                         cn6, energy_w0, epolar, f_scratch, xx, ix, &
                         ipairs, virvsene, xx(lpol), &
                         xx(lpol2), .false. , cn3, cn4, cn5)
        ! call with full charges but no vdw interaction
        ! between solute and solvent
        call ewald_force(x, natom, ix(i04), ix(i06), xx(l15), cn1_lrt, &
                         cn2_lrt, cn6, eelt, epolar, f_scratch, xx, ix, &
                         ipairs, virvsene, xx(lpol), &
                         xx(lpol2), .false. , cn3, cn4, cn5)
        energy_vdw0 = evdw
        call lrt_solute_sasa(x,natom, xx(l165))
      end if

      ! call normal_ewald force this will overwrite everything 
      ! computed above except energy_m0 and energy_w0
      call ewald_force(x, natom, ix(i04), ix(i06), xx(l15), cn1, cn2, &
                       cn6, eelt, epolar, f, xx, ix, ipairs, &
                       virvsene, xx(lpol), &
                       xx(lpol2), .false. , cn3, cn4, cn5)
      energy_vdw0 = evdw - energy_vdw0

      ! count call to ltr, maybe calculate Eee and print it
      call ee_linear_response(eelt, master)
    else ! just call ewald_force normally
      call ewald_force(x, natom, ix(i04), ix(i06), xx(l15), cn1, cn2, &
                       cn6, eelt, epolar, f, xx, ix, ipairs, &
                       virvsene, xx(lpol), &
                       xx(lpol2), .false. , cn3, cn4, cn5)
    end if ! ilrt /= 0

    call timer_stop(TIME_EWALD)
#ifdef MPI
      if (mytaskid == 0) then
#endif
        pot%vdw   = evdw
        pot%elec  = eelt
        pot%hbond = ehb  !whereis ehb?
#ifdef MPI
      else
        ! energies have already been reduced to the master
        ! node in ewald_force, so here we zero out elements
        ! on non-master nodes:
        pot%vdw      = 0.d0
        pot%elec     = 0.d0
        pot%hbond    = 0.d0
      end if
#endif
   end if  ! ( igb == 0 .and. ipb == 0 .and. iyammp == 0 )

   ! End of non-bonded computations
   call timer_stop(TIME_NONBON)

   ! Calculate other contributions to the forces

  ! Bonds with H
  call timer_start(TIME_BOND)

#ifdef MPI /* SOFT CORE */
  ! zero only once, sc bond energy is sum of H and non-H terms
  sc_ener(1) = 0.0d0
#endif
  if (ntf < 2) then
    ebdev = 0.d0
    call bond(nbonh, ix(iibh), ix(ijbh), ix(iicbh), x, f, ene(6))
    pot%bond = pot%bond + ene(6)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesb
    endif
#  endif
#endif
  end if

  ! Bonds without H
  if (ntf < 3) then
    call bond(nbona+nbper, ix(iiba), ix(ijba), ix(iicba), x, f, ene(7))
    pot%bond = pot%bond + ene(7)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesb
    endif
#  endif
#endif
    if (nbonh+nbona > 0) then
      ebdev = sqrt(ebdev/(nbonh+nbona))
    end if
  end if

  ! Angles with H
  if (ntf < 4) then

#ifdef MPI /* SOFT CORE */
    ! zero only once, sc bond energy is sum of H and non-H terms
    sc_ener(2) = 0.0d0
#endif
    eadev = 0.d0

    call angl(ntheth, ix(i24), ix(i26), ix(i28), ix(i30), x, f, ene(8))
    pot%angle = pot%angle + ene(8)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesa
    endif
#  endif
#endif
  end if

  ! Angles without H
  if (ntf < 5) then
    call angl(ntheta+ngper, ix(i32), ix(i34), ix(i36), ix(i38), x, f, &
              ene(9)) 
    pot%angle = pot%angle + ene(9)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesa
    end if
#  endif
#endif
    if (ntheth+ntheta > 0) then
      eadev = 57.296*sqrt( eadev/(ntheth+ntheta) )
    end if
  end if

  ! Dihedrals with H
  if (ntf < 6) then

#ifdef MPI /* SOFT CORE */
    ! zero only once, sc bond energy is sum of H and non-H terms
    sc_ener(3) = 0.0d0
#endif
    call ephi(nphih, ix(i40), ix(i42), ix(i44), ix(i46), ix(i48), xx(l15), &
              ix(i04), x, f, dvdl, ene(10), ene(11), ene(12), xx(l190))

    ! Combine contributions from dihedrals with H, including
    ! the torsions themselves and 1:4 vdW / elec interactions.
    pot%dihedral = pot%dihedral + ene(10)
    pot%vdw_14   = pot%vdw_14   + ene(11)
    pot%elec_14  = pot%elec_14  + ene(12)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesd
    endif
#  endif
#endif
  end if

  ! Dihedrals without H
  if (ntf < 7) then
    call ephi(nphia+ndper, ix(i50), ix(i52), ix(i54), ix(i56), ix(i58), &
              xx(l15), ix(i04), x, f, dvdl, ene(13), ene(14), ene(15), &
              xx(l190))

    ! Combine contributions from dihedrals without H, including
    ! the torsions themselves and 1:4 vdW / elec interactions.
    pot%dihedral = pot%dihedral + ene(13)
    pot%vdw_14   = pot%vdw_14   + ene(14)
    pot%elec_14  = pot%elec_14  + ene(15)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
       pot%les = pot%les + elesd
    endif
#  endif
#endif

    ! Do CHARMM impropers, if the CHARMM force field is in use
    ! Note: CHARMM does not distinguish between impropers with and 
    ! without hydrogen.  Hence, it is not possible to strictly
    ! conform to the ntf options of sander. Here CHARMM impropers
    ! are calculated as long as ntf < 7 - so CHARMM impropers are
    ! essentially considered to be in the same set as dihedrals
    ! NOT involving hydrogen.
    if (charmm_active) then
      call charmm_calc_urey_bradley(x, pot%angle_ub, f)
      call charmm_calc_impropers(x, pot%imp, f)
      call charmm_calc_cmap(x, pot%cmap, f)
    end if
    ! End CHARMM impropers computation

    if (cmap_active) then
      call calc_cmap(x,pot%cmap,f)
    end if
  end if
  ! End computations for ntf < 7

  ! End of valence contribution computations
  call timer_stop(TIME_BOND)

  ! Step 1/2 to save emap + entr + ecap forces

  ! Calculate the position constraint energy
#ifdef MPI /* SOFT CORE */
  ! Zero all restraint/constraint energies
  sc_ener(14:19) = 0.d0
#endif
  if (natc > 0 .and. ntr==1) then   ! ntr=1 (positional restraints)
    call xconst(natc, entr, ix(icnstrgp), x, f, xx(lcrdr), xx(l60))
    pot%constraint = entr
  end if
  if (itgtmd==1 .and. (nattgtfit > 0 .or. nattgtrms > 0)) then

    ! Calculate rmsd for targeted md (or minimization) if requested.
    ! All nodes do rms fit, could just be master then broadcast.
    ! All nodes need all coordinates for this.
    call rmsfit(xx(lcrdr), x, xx(lmass), ix(itgtfitgp), ix(itgtrmsgp), &
                rmsdvalue, nattgtrms, nattgtfit, rmsok)
    if (.not. rmsok) then
      if (master) then
        write (6,*) 'Fatal error: Error calculating RMSD!'
      end if
      call mexit(6, 1)
    end if

    call xtgtmd(entr, ix(itgtrmsgp), x, f, xx(lcrdr), xx(lmass), tgtrmsd, &
                tgtmdfrc, rmsdvalue, nattgtrms)
    pot%constraint = entr
  else if (itgtmd == 2) then
    call mtmdcall(entr, xx(lmtmd01), ix(imtmd02), x, f, ih(m04), ih(m02), &
                  ix(i02), ih(m06), xx(lmass), natom, nres, 'CALC')
    pot%constraint = entr
  end if

  if (ifcap == 1 .or. ifcap == 2) then
    call capwat(natom, x, f, ecap)
    pot%constraint = pot%constraint + ecap
  else if (ifcap == 3) then
    write(6,*) 'No energy expression for spherical boundary known yet'
    call mexit(6,1)
  else if(ifcap == 4) then
    write(6,*) 'No energy expression for orthorhombic boundary known yet'
    call mexit(6,1)
    !call orth(natom,ix(ibellygp),x,f,eorth)
    !ene(20) = ene(20) + eorth
  end if

  ! No energy expression for ifcap == 5 given because only
  !    one step of minimization is allowed with this.

  if (igb == 0 .and. ipb == 0 .and. iyammp == 0) then
    ener%virvsene = virvsene
    ener%diprms = diprms
    ener%dipiter = dipiter
    ener%dipole_temp = dipole_temp
  end if

  ! Get the noesy volume penalty energy
  pot%noe = 0.d0
  if (iredir(4) /= 0) then
    call timer_start(TIME_NOE)
    call noecalc(x, f, xx, ix)
    call timer_stop(TIME_NOE)
  end if
  ! Do we need a pot%noe here?  mjw TODO

  ! When any form of Generalized Born is active but Poisson-Boltzman is
  ! not, all nonbonded interactions are done in the subroutine egb:
  esurf = 0.d0
  if (igb /= 0 .and. igb /= 10 .and. ipb == 0) then
    call timer_start(TIME_EGB)
    call egb(x, f, rborn, fs, reff, onereff, xx(l15), ix(i04), ix(i06), &
             ix(i08), ix(i10), xx(l190), cut, ntypes, natom, natbel, epol, &
             eelt, evdw, esurf, dvdl, xx(l165), ix(i82), xx(l170), xx(l175), &
             xx(l180), xx(l185), ncopy &
#ifndef LES
             , xx(l2402),xx(l2403),xx(l2404) &
#endif
             )
    pot%vdw  = evdw
    pot%elec = eelt
    pot%gb   = epol
    pot%surf = esurf
    pot%dvdl = dvdl
    call timer_stop(TIME_EGB)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesp
    endif
#  endif
#endif
  end if  ! ( igb /= 0 .and. igb /= 10 .and. ipb == 0 )
  ! End handoff of nonbonded computations to subroutine egb

  ! Force and energy computations by the Reference Interaction Site Model
  if (rismprm%rism == 1) then
    call timer_start(TIME_RISM)
    call rism_force(x, f, erism, irespa, imin)
    pot%rism = erism
    call timer_stop(TIME_RISM)
  endif

  if (master) then
    !  These parts of the NMR energies are not parallelized, so only
    !  are done on the master node:
    eshf = 0.d0
    epcshf = 0.d0
    ealign = 0.d0
    ecsa = 0.d0
    if (iredir(5) /= 0) call cshf(natom,x,f)
    if (iredir(7) /= 0) call pcshift(natom,x,f)
    if (iredir(9) /= 0) call csa1(natom,x,f)
    if (iredir(8) /= 0) call align1(natom,x,f,xx(lmass))
    if (temap) then
      call emapforce(natom, enemap, xx(lmass), x, f)
      pot%emap = enemap
    end if

  end if

  ! MuSiC - GAL17 force field
  call music_force(ipairs, music_vdisp, music_vang, music_vgauss, music_spohr89)

  ! Built-in X-ray target function and gradient
  xray_energy = 0.d0
  if( xray_active .and. master) then
     call xray_get_derivative(x3,f3,xray_nstep,xray_e)
     xray_nstep = xray_nstep + 1
  endif

#ifdef MPI
  call timer_start(TIME_COLLFRC)
  call timer_barrier( commsander )
  !     add force, ene, vir, copies from all nodes
  !            also add up newbalance for nonperiodic.

  ! Remember to work on the local instance of the
  ! potential energy array, i.e. pot and NOT the global one,
  ! i.e. ener%pot
  call fdist(f,xx(lfrctmp),pot,vir,newbalance,3*natom+iscale)
  call timer_stop(TIME_COLLFRC)
#endif

  ! ---- at this point, the parallel part of the force calculation is
  !      finished, and the forces have been distributed to their needed
  !      locations.  All forces below here are computed redundantly on
  !      all processors, and added into the force vector.  Hence, below
  !      is the place to put any component of the force calculation that
  !      has not (yet) been parallelized.

  ! Calculate the NMR restraint energy contributions, if requested.
  ! (Even though this is not parallelized, it needs to be run on all
  ! threads, since this code is needed for weight changes as well as
  ! for NMR restraint energy analysis.  The whole thing could stand a
  ! major re-write....)
  if (nmropt > 0) then
    call nmrcal(x, f, ih(m04), ih(m02), ix(i02), xx(lwinv), enmr, devdis, &
                devang, devtor, devplpt, devpln, devgendis, temp0, tautp, &
                cut, xx(lnmr01), ix(inmr02), xx(l95), 31, 6, rk, tk, pk, cn1, &
                cn2, asol, bsol, xx(l15), numbnd, numang, nptra-nimprp, &
                nimprp, nphb, natom, natom, ntypes, nres, rad, wel, radhb, &
                welhb, rwell, tgtrmsd, temp0les, -1, 'CALC')
  end if
#ifdef MPI
  call mpi_reduce(enoe, pot%noe, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, &
                  commsander, ierr)
  enoe = pot%noe ! so all processors now have the full enoe value
#else
  pot%noe = enoe
#endif

#ifdef DSSP
  if (idssp > 0) then
    call fdssp(natom, x, f, edssp)
  else
    edssp = 0.d0
  end if
#endif

    ! Calculate the total energy and group the components
#ifndef LES
  if (igb == 0 .and. ipb == 0) then
    pot%vdw_14   = pot%vdw_14   + enb14
    pot%elec_14  = pot%elec_14  + ee14
  endif
#endif
  pot%constraint = pot%constraint + eshf + epcshf + pot%noe + &
       sum(enmr(1:6)) + ealign + ecsa + pot%emap + xray_energy + enfe
#ifdef DSSP
  pot%constraint = pot%constraint + edssp
#endif
  pot%polar = epolar
  pot%tot   = pot%vdw + pot%elec + pot%gb + pot%pb + pot%bond + pot%angle + &
              pot%dihedral + pot%vdw_14 + pot%elec_14 + pot%hbond + &
              pot%constraint + pot%rism + pot%ct
  pot%tot = pot%tot + pot%polar + pot%surf + pot%scf + pot%disp

  !Charmm related
  pot%tot = pot%tot + pot%angle_ub + pot%imp + pot%cmap 

  ! MuSiC - GAL17 force field
  pot%tot = pot%tot + music_vdisp + music_vang + music_vgauss + music_spohr89

  ! The handover
  ener%pot = pot
  ener%aveper = aveper
  ener%aveind = aveind
  ener%avetot = avetot
   
  ! If freezemol has been set, zero out all of the forces for
  ! the real atoms. (It is no longer necessary to set ibelly.)
  if (ifreeze > 0) then
    do i = 1, 3*natom
      f(i) = 0.d0
    end do
  end if

  ! If a bellymask is being used, set the belly atom forces to zero.
  if (belly) call bellyf(natom,ix(ibellygp),f)

  ! Dump forces in CHARMM format if appropriate and requested
  if (charmm_active) then
    if (do_charmm_dump_gold == 1) then
      call charmm_dump_gold(f, natom, ener)
    endif 
  end if

  ! End force computations and exit
  call timer_stop(TIME_FORCE)
  return

end subroutine force


