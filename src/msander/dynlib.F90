! <compile=optimized>
#include "../include/dprec.fh"
#include "../include/assert.fh"

!------------------------------------------------------------------------------
! corpac: Emit (write out) coordinates or velocities, r(istart:n), to unit nf.
!
! Arguments:
!   r:      the coordinates or velocities to write out
!   istart: the index of the data at which to start (i.e. "start at atom 6")
!   n:      the length of the coordinate or velocity data
!   nf:     the number of the output file unit
!   loutfm: the logical translation of ioutfm as seen in mdin, describing the
!           output file format.  In ioutfm, 0 = .crd and 1 = netcdf, but loutfm
!           is set to .true. if ioutfm <= 0 (.true. for formatted writes) and
!           false otherwise (binary writes).
!------------------------------------------------------------------------------
subroutine corpac(r, istart, n, nf, loutfm)
  use bintraj

  implicit none
  integer, intent(in) :: istart, n, nf
  _REAL_, intent(in) ::  r(n)
  logical, intent(in) ::  loutfm  ! true for formatted output

  integer i,j
  logical three_digit_fractional_part
  _REAL_  rmax,rmin
  parameter (rmax = 9999.99d0)
  parameter (rmin = -999.99d0)

  ! Bail out if the start position is beyond the length of the array
  if (istart > n) then
    return
  end if


  if (.not. loutfm) then

    ! Unformatted writes
    call write_binary_traj(r, istart, n, nf)
  else

    ! Formatted writes
    three_digit_fractional_part = .true.
    do i = istart, n
      if (r(i) > rmax .or. r(i) < rmin) then
        three_digit_fractional_part = .false.
        exit
      end if
    end do

    if (three_digit_fractional_part) then
      write(nf, '(10f8.3)') (r(j), j = istart, n)
    else
      write(nf, '(10f8.2)') (r(j), j = istart, n)
    end if
  end if

  return
end subroutine corpac

!------------------------------------------------------------------------------
! open_dummp_files: open the coordinate, velocity, energy and constant pH
!                   and constant Redox potential output files.
!------------------------------------------------------------------------------
subroutine open_dump_files
  use bintraj, only: open_binary_files
  use file_io_dat

  implicit none
#  include "../include/md.h"
#  include "extra.h"

  ! Only the master process does file I/O
  if (master) then
    if (ioutfm <= 0) then

      ! Formatted dumping
      if (ntwx > 0) then
        call amopen(MDCRD_UNIT, mdcrd, 'U', 'F', facc)
        if (facc /= 'A') then
          write(MDCRD_UNIT, 1000) title
        end if
      end if
      if (ntwv > 0) then
        call amopen(MDVEL_UNIT, mdvel, owrite, 'F', 'W')
        write(MDVEL_UNIT, 1000) title
      end if
      if (ntwf > 0) then
        call amopen(MDFRC_UNIT,mdfrc,owrite,'F','W')
        if (facc /= 'A') then
          write(MDFRC_UNIT, 1000) title
        end if
      end if
    else
      call open_binary_files
    end if  ! (ioutfm <= 0)
    if (ntwe > 0) then
      call amopen(MDEN_UNIT, mden, owrite, 'F', 'W')
    end if
  end if
  ! End contingency for work on the master process

  1000 format(a80)

  return
end subroutine open_dump_files

!------------------------------------------------------------------------------
! close_dump_files: close the coordinate, velocity, energy, constant pH
!                   and constant Redox potential output files.
!------------------------------------------------------------------------------
subroutine close_dump_files
  use bintraj, only: close_binary_files
  use file_io_dat

  implicit none
#  include "extra.h"
#  include "../include/md.h"

  if (master) then
    if (ioutfm > 0) then
      call close_binary_files
    else
      if (ntwx > 0) close(MDCRD_UNIT)
      if (ntwv > 0) close(MDVEL_UNIT)
      if (ntwf > 0) close(MDFRC_UNIT)
    end if
    if (ntwe > 0) close(MDEN_UNIT)
    if (imin == 5) close(INPTRAJ_UNIT)
    if (ntpr > 0) close(7)
  end if

  return
end subroutine close_dump_files

!------------------------------------------------------------------------------
! prntmd: energy output for md, in human-readable form.
!
! Arguments:
!   nstep:     the current time step
!   time:      the simulation time in ps
!   ener:      struct storing various energy components (see state.F90)
!   onefac:    the inverse of the degrees of freedom factors (a three-element
!              array of _REAL_ numbers).  See runmd.F90 and amoebarunmd.F90
!              to understand how each element is set.
!   iout7:     the output file unit for mdinfo
!   rms:       flag to have RMS statistics for each quantity printed (disabled
!              in Locally Enhanced Sampling, LES)
!------------------------------------------------------------------------------
subroutine prntmd(nstep, time, ener, onefac, iout7, rms)
  use file_io_dat

#if defined( MPI )
  use remd, only : rem, repnum, mdloop, stagid, &
                   remd_dimension, group_num, replica_indexes
#endif /* MPI */
#ifdef LES
  use les_data, only : temp0les
#endif
  use sander_rism_interface, only: rismprm, RISM_FULL, &
                                   rism_calc_type, &
                                   rism_thermo_print

  use xray_interface_module, only: xray_write_md_state
  use xray_globals_module, only: xray_active
  use state
  use charmm_mod, only: charmm_active
  use sgld, only: isgld,sglabel,sgfti,sgffi
  use ff11_mod, only: cmap_active
  use emap,only: temap
  use md_scheme, only: ntt

  implicit none

  ! Passed in
  integer, intent(in) :: nstep, iout7
  _REAL_,  intent(in) :: onefac(3), time
  type(state_rec), intent(in) :: ener
  logical, intent(in) :: rms ! true if this is a print of RMS fluctuations

  ! Local
  integer :: total_nstlim, total_nstep
  _REAL_ :: enb14, edihed, eangle, ebond, ehbond, eel, enonb, epot, egb, epb
  _REAL_ :: virt
  _REAL_ :: dvdl, esurf, epol, econst, eel14, ekcmt
  _REAL_ :: edisp, escf, enemap
  _REAL_ :: press, densit, volume
  _REAL_ :: ektot, etot
#ifdef LES
  _REAL_ :: eksolt, eksolv
#else
  _REAL_ :: rms_pbs
#endif
  _REAL_ :: erism
  _REAL_ :: pot_array(potential_energy_rec_len)
  _REAL_ :: ect
  _REAL_ :: amd_boost
#ifdef MPI
  integer :: n
#endif /* MPI */

#  include "../include/md.h"
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#  include "tgtmd.h"

  ! The logical argument rms is always passed in, but only used
  ! if LES is NOT defined (which is most of the builds).  However,
  ! when LES is defined, this creates compiler warnings, which I
  ! would rather silence so that it is easier to see if something
  ! is truly wrong.  This code will make use of rms so that the
  ! compiler will not complain.  etot will be immediately redefined
  ! with the same value after this, so it has no effect.
  if (rms) then
    etot = ener%tot
  end if

  ! Define various terms
  etot   = ener%tot
  ektot  = ener%kin%tot

  temp = ener%kin%tot*onefac(1)
#ifdef LES
  eksolt = ener%kin%solt*onefac(2)

  ! Modified for Locally Enhanced Sampling, ener(4) is
  ! LES KE now, not solvent KE, so it should be reported
  ! along with temperature for LES region.
  if (temp0les < 0.d0) then
    eksolv = 0.0d0
  else
    eksolv = ener%kin%solv*onefac(3)
  end if

#else
  rms_pbs = ener%kin%solv
#endif

  volume  = ener%volume
  densit  = ener%density

  press   = ener%pres(4)
  ekcmt   = ener%cmt(4)
  virt    = ener%vir(4)

  epot    = ener%pot%tot
  enonb   = ener%pot%vdw
  eel     = ener%pot%elec
  ehbond  = ener%pot%hbond
  egb     = ener%pot%gb
  epb     = ener%pot%pb
  ebond   = ener%pot%bond
  eangle  = ener%pot%angle
  edihed  = ener%pot%dihedral
  enb14   = ener%pot%vdw_14
  eel14   = ener%pot%elec_14
  econst  = ener%pot%constraint
  epol    = ener%pot%polar
  esurf   = ener%pot%surf
  dvdl    = ener%pot%dvdl
  escf    = ener%pot%scf
  edisp   = ener%pot%disp
  enemap   = ener%pot%emap
  amd_boost = ener%pot%amd_boost
  erism   = ener%pot%rism
  ect     = ener%pot%ct

  write(6, 9018) nstep,time,temp,press
  write(6, 9028) etot,ektot,epot
  write(6, 9038) ebond,eangle,edihed
  if (charmm_active) then
    write(6, 9039) ener%pot%angle_ub, ener%pot%imp, ener%pot%cmap
  end if
  write(6, 9048) enb14, eel14, enonb

  if (igb == 0 .and. ipb == 0 .and. rismprm%rism == 0) then
    write(6, 9058) eel, ehbond, econst
  else if (igb == 10 .or. ipb /= 0) then
    write(6, 9060) eel, epb, econst
  else if (rismprm%rism == 1) then
    write(6, 9061) eel, erism, econst
  else
    write(6, 9059) eel, egb, econst
  end if


  if (gbsa > 0) then
    write(6, 9077) esurf
  end if
  if (igb == 10 .or. ipb /= 0) then
    write(6, 9074) esurf, edisp
  end if
  if (icfe > 0) then
    write(6, 9100) dvdl
  end if

#ifndef LES
  if (rms .and. ntt == 0) then
    write(6, 9075) rms_pbs
  end if
#endif /* LES is not defined */
  if (ntp > 0.0 .or. volume /= 0.0) then
    write(6, 9078) ekcmt, virt, volume
  end if
#ifdef LES

  ! If LES T bath is used, write temperature for LES and non-LES
  if (temp0les >= 0.d0) then
    write(6, 9067) eksolt, eksolv
  end if
#endif /* LES */

  if (cmap_active .and. epol /= 0.0) then
    write(6, 9068) epol, ener%pot%cmap
  else
    if (cmap_active) then
      write(6, 9069) ener%pot%cmap
    end if
    if (epol /= 0.0) then
      write(6, 9070) epol
    end if
  end if

  if (ntp > 0.0 .or. volume /= 0.0) then
    write(6, 9079) densit
  end if

  if (itgtmd == 1) then
    write (6,'(a,f8.3)') "Current RMSD from reference: ", rmsdvalue
    write (6,'(a,f8.3)') "Current target RMSD:         ", tgtrmsd
  end if

  ! Printout SGLD guiding information
  if (isgld > 0) then
    write(6, 1005) sglabel,ener%sgld%sgscale, ener%sgld%templf, ener%sgld%temphf, &
                   ener%sgld%epotlf, ener%sgld%epothf, ener%sgld%epotllf, &
                   ener%sgld%sgwt
  endif
  if (xray_active) call xray_write_md_state(6)

  if (econst /= 0.0) then
    write(6, 9076) epot-econst
  end if

#ifdef MPI
  ! Print current REMD info (replica#, temp0, excgh#) only for iout7 > 0
  ! (Not for average/rms)
  if (rem /= 0 .and. rem /= 4 .and. rem /= 5 .and. rem /= -1 .and. iout7 > 0) then
    if (isgld > 0) then
      write (6, 9064) temp0, sgfti, sgffi, stagid, repnum, mdloop
    else
      write (6, 9065) temp0, repnum, mdloop
    endif
  else if (rem == 4 .and. iout7 > 0) then
    write (6, 9066) solvph, repnum, mdloop
  else if (rem == 5 .and. iout7 > 0) then
    write (6, 9071) solve, repnum, mdloop
  else if (rem == -1 .and. iout7 > 0) then
    write(6, 967) remd_dimension, mdloop
    write(6, 968) (group_num(n), n=1, remd_dimension)
    write(6, 969) (replica_indexes(n), n=1, remd_dimension)
  endif
#endif
  ! wxw: EMAP energy
  if (temap ) then
    write (6, 9062) enemap
  endif
  ! write(6, 8088)

  ! Flush i/o buffer
  call flush(6)
  if (iout7 == 0) return

  ! Output the info file if requested
  write(7, 9018) nstep, time, temp, press
  write(7, 9028) etot, ektot, epot
  write(7, 9038) ebond, eangle, edihed
  ! CHARMM specific energy terms
  if (charmm_active) then
    write(7, 9039) ener%pot%angle_ub, ener%pot%imp, ener%pot%cmap
  end if

  write(7, 9048) enb14, eel14, enonb
  if (igb == 0 .and. ipb == 0 .and. rismprm%rism == 0) then
    write(7, 9058) eel, ehbond, econst
  else if ( igb == 10 .or. ipb /= 0) then
    write(7, 9060) eel, epb, econst
  else if (rismprm%rism == 1) then
    write(7, 9061) eel, erism, econst
  else
    write(7, 9059) eel, egb, econst
  end if
  if (xray_active) call xray_write_md_state(7)

#ifdef MPI

  ! Print current REMD info (replica#, temp0, excgh#) only for iout7 > 0
  ! (Not for average/rms)
  if (rem /= 0 .and. rem /= 4 .and. rem /= 5 .and. rem /= -1 .and. iout7 > 0) then
    if (isgld > 0) then
      write(7, 9064) temp0, sgfti, sgffi, stagid, repnum, mdloop
    else
      write(7, 9065) temp0, repnum, mdloop
    endif
  else if (rem == 4) then
    write (7, 9066) solvph, repnum, mdloop
  else if (rem == 5) then
    write (7, 9071) solve, repnum, mdloop
  else if (rem == -1 .and. iout7 > 0) then
    write(7, 967) remd_dimension, mdloop
    write(7, 968) (group_num(n), n=1, remd_dimension)
    write(7, 969) (replica_indexes(n), n=1, remd_dimension)
  endif
#endif /* MPI */
  ! wxw: EMAP energy
  if (temap .and. iout7 > 0) then
    write (7, 9062) enemap
  endif

   if (gbsa > 0) then
     write(7, 9077) esurf
   end if
   if (igb == 10 .or. ipb /= 0) then
     write(7, 9074) esurf, edisp
   end if
   if (econst /= 0.0) then
     write(7, 9076) epot-econst
   end if
   if (icfe > 0) then
     write(7, 9100) dvdl
   end if
#ifndef LES
   if (rms .and. ntt == 0) then
     write(6, 9075) rms_pbs
   end if
#endif
   if (ntp > 0.0 .or. volume /= 0.0) then
     write(7, 9078) ekcmt, virt, volume
   end if
#ifdef LES

  ! If LES T bath is used, write temperature for LES and non-LES
  if (temp0les >= 0.d0) then
    write(7, 9067) eksolt, eksolv
  end if
#endif

  if (cmap_active .and. epol /= 0.0) then
    write(7,9068) epol, ener%pot%cmap
  else
    if (cmap_active) then
      write(7, 9069) ener%pot%cmap
    end if
    if (epol /= 0.0) then
      write(7,9070) epol
    end if
  endif
  if (ntp > 0.0 .or. volume /= 0.0) then
    write(7, 9079) densit
  end if

  ! Printout SGLD guiding information
  if (isgld > 0) then
    write(7, 1005) sglabel,ener%sgld%sgscale, ener%sgld%templf, ener%sgld%temphf, &
                   ener%sgld%epotlf, ener%sgld%epothf, ener%sgld%epotllf, &
                   ener%sgld%sgwt
  endif
  call nmrptx(7)

  if (rismprm%rism == 1 .and. rismprm%write_thermo == 1) then
    if (rism_calc_type(nstep) == RISM_FULL) then
      call rism_thermo_print(.false., transfer(ener%pot, pot_array))
    end if
  end if

#ifndef NO_DETAILED_TIMINGS
  !Print Timing estimates to mdinfo.
#  ifdef MPI
  if (rem .ne. 0) then
    total_nstlim = nstlim * numexchg
    total_nstep = nstlim * (mdloop-1) + nstep
  else
    total_nstlim = nstlim
    total_nstep = nstep
  end if
#  else
  total_nstlim = nstlim
  total_nstep = nstep
#  endif
#  ifndef API
  if (total_nstep /= 0) then
    call print_ongoing_time_summary(total_nstlim, nstep, dt, 7)
  end if
#  endif
#endif

  8088 format(t2,78('-'),/)
  9018 format(/1x, 'NSTEP =',i9,3x,'TIME(PS) =',f12.3,2x, &
         'TEMP(K) =',f9.2,2x,'PRESS =',f8.1)
  9028 format (1x,'Etot   = ',f14.4,2x,'EKtot   = ',f14.4,2x, &
         'EPtot      = ',f14.4)
  9038 format (1x,'BOND   = ',f14.4,2x,'ANGLE   = ',f14.4,2x, &
         'DIHED      = ',f14.4)
  9039 format(1x, 'UB     = ', f14.4, 2x, 'IMP     = ', f14.4, 2x, &
         'CMAP       = ', f14.4)
  9048 format (1x,'1-4 NB = ',f14.4,2x,'1-4 EEL = ',f14.4,2x, &
         'VDWAALS    = ',f14.4)
  9058 format (1x,'EELEC  = ',f14.4,2x,'EHBOND  = ',f14.4,2x, &
         'RESTRAINT  = ',f14.4)
  9059 format (1x,'EELEC  = ',f14.4,2x,'EGB     = ',f14.4,2x, &
         'RESTRAINT  = ',f14.4)
  9060 format (1x,'EELEC  = ',f14.4,2x,'EPB     = ',f14.4,2x, &
         'RESTRAINT  = ',f14.4)
  9061 format (1x,'EELEC  = ',f14.4,2x,'ERISM   = ',f14.4,2x, &
         'RESTRAINT  = ',f14.4)
  9062 format (1x,'EMAP   = ',f14.4)
  9180 format (1x, 'EAMD_BOOST  = ', f14.4)
#ifdef MPI
  ! EVB formats
   777 format( 1x, 'V_TOT  = ', f14.4, 2x, 'V_EVB   = ', f14.4, 2x &
              , 'V_MAP      = ', f14.4 )
   888 format( 1x, 'V_TOT  = ', f14.4, 2x, 'V_EVB   = ', f14.4, 2x &
              , 'V_UMB      = ', f14.4 )
   999 format( 1x, A, (4(2x,f14.4)) )
  ! DAN ROE: Added Temp, Rep, Exchange
  ! Xiongwu: add sgft, tempsg for RXSGLD
   9064 format (1x,'TEMP= ',f6.1,1x,'SGFT= ',f7.4, 1X,'SGFF= ',f7.4, 1X, &
   'STAGE= ',i4,1x,'REPNUM= ',i4,1x,'EXCH= ',i6)
  9065 format (1x,'TEMP0  = ',f14.4,2x,'REPNUM  = ',i14,2x, &
         'EXCHANGE#  = ',i14)
  9066 format (1x,'SOLVPH = ',f14.4,2x,'REPNUM  = ',i14,2x, &
         'EXCHANGE#  = ',i14)
  9071 format (1x,'SOLVE  = ',f14.6,2x,'REPNUM  = ',i14,2x, &
         'EXCHANGE#  = ',i14)
   967 format(1x,'REMD_DIMENSION = ',i14,2x,'EXCHANGE#  = ',i14)
   968 format(1x,'MULTI-D REMD. GROUP INDEXES:   ',12i5)
   969 format(1x,'MULTI-D REMD. REPLICA INDEXES: ', 12i5)
#endif
  9072 format(52x, 'SURFTEN    = ', f14.4)
  9074 format (1x,'ECAVITY= ',f14.4,2x,'EDISPER = ',f14.4)
#ifndef LES
  9075 format ('|E(PBS) = ',f14.4)
#endif /* LES not defined */
  9076 format (1x,'EAMBER (non-restraint)  = ',f14.4)
  9077 format (1x,'ESURF= ',f14.4)
  9078 format (1x,'EKCMT  = ',f14.4,2x,'VIRIAL  = ',f14.4,2x, &
         'VOLUME     = ',f14.4)
  9079 format (52x,'Density    = ',f14.4)
  9080 format (1x,'PM3ESCF= ',f14.4)
  9081 format (1x,'AM1ESCF= ',f14.4)
  9981 format (1x,'AM1DESCF= ',f14.4)
  9082 format (1x,'MNDOESCF= ',f13.4)
  9982 format (1x,'MNDODESCF= ',f13.4)
  9083 format (1x,'PDDGPM3-ESCF= ',f12.4)
  9084 format (1x,'PDDGMNDO-ESCF= ',f11.4)
  9085 format (1x,'PM3CARB1-ESCF= ',f11.4)
  9086 format (1x,'DFTBESCF= ',f13.4)
  9087 format (1x,'RM1ESCF= ',f14.4)
  9088 format (1x,'PDDGPM3_08-ESCF= ',f12.4)
  9089 format (1x,'PM6ESCF= ',f14.4)
  9090 format (1x,'PM3MMXESCF= ',f14.4)
  9091 format (1x,'PM3ZNBESCF= ',f14.4)
  9092 format (1x,'EXTERNESCF= ',f14.4)
  9093 format (1x,'PM3MAISESCF= ',f14.4)
  9096 format (1x,'PM3MMX2ESCF= ',f14.4)
  9099 format (1x,'ECRG   = ',f14.4)

#ifdef LES
  ! LES and non-LES temperatures (no solvent/solute)
  9067 format (1x,'T_non-LES =',f10.4,2x, 'T_LES     = ',f10.4)
#endif
  9068 format(1x, 'EPOLZ  = ',f14.4,2x,'CMAP    = ',f14.4)
  9069 format(1x, 'CMAP   = ',f14.4)
  9070 format (1x,'EPOLZ  = ',f14.4)
  9190 format(1x,'Dipole convergence: rms = ',e10.3,' iters = ',f6.2)
  9191 format(1x,'Dipole convergence: rms = ',e10.3, &
         ' temperature = ',f6.2)
  9100 format (1x,'DV/DL  = ',f14.4)
  9188 format (1x,'Ewald error estimate: ', e12.4)
  1005 format(1x,A6,F9.4,F8.2,F8.2,X,F12.2,F12.2,F12.2,F10.4)

  return
end subroutine prntmd

!------------------------------------------------------------------------------
! setvel: assign velocities from a Maxwellian distribution for initialization
!         or Andersen thermocoupling.
!
! Arguments:
!   nr:      number of atoms
!   v:       velocity array
!   winv:    inverse mass array
!   tempi:   temperature to randomize velocities to
!   iscale:  nmr restraint-related variable
!   scalm:   nmr restraint-related variable
!------------------------------------------------------------------------------
subroutine setvel(nr, v, winv, tempi, iscale, scalm)
   use random, only: gauss
   implicit none
#  include "../include/memory.h"

  ! Passed variables
  _REAL_, intent(out) :: v(*)
  _REAL_, intent(in) :: winv(*), tempi, scalm
  integer, intent(in) :: nr, iscale

  ! Local variables
  !  boltz : boltzmann's constant in appropriate units
  !  sd    : Maxwell-Boltzmann factor from which to derive temperature
  integer :: i, j, nr3
  _REAL_ :: boltz, sd

  nr3 = 3 * nr

  ! Assign velocities from a Maxwellian distribution
  if (tempi < 1.d-6) then
    v(1:nr3+iscale) = 0.d0
    return
  end if
  boltz = (8.31441d-3) * tempi/4.184d0
  i = 0
  do j = 1, nr
    sd = sqrt(boltz*winv(j))
    call gauss(0.d0, sd, v(i+1))
    call gauss(0.d0, sd, v(i+2))
    call gauss(0.d0, sd, v(i+3))
    i = i+3
  end do
  if (iscale > 0) then
    sd = sqrt(boltz/scalm)
    do j = 1, iscale
      call gauss(0.d0, sd, v(i+j))
    end do
  end if

  return
end subroutine setvel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Energy output for mden file, in human-readable form.
subroutine mdeng(nf,nstep,time,ener,onefac,ntp)
   use state
   use constants, only : zero
   implicit none
#  include "box.h"
   integer nf,nstep,ntp,i
   _REAL_ onefac(3),time
   type(state_rec), intent(in) :: ener
   logical first
   character(len=16) labs(42)
   save first
   save labs
   data first/.true./
   data labs/'Nsteps  ','time(ps)  ','Etot  ','EKinetic  ', &
         'Temp  ','T_solute ','T_solv  ','Pres_scal_solu ', &
         'Pres_scal_solv ','BoxX  ','BoxY  ','BoxZ  ', &
         'volume  ','pres_X  ','pres_Y  ','pres_Z  ', &
         'Pressure ','EKCoM_x ','EKCoM_y ','EKCoM_z', &
         'EKComTot ','VIRIAL_x ','VIRIAL_y ','VIRIAL_z ', &
         'VIRIAL_tot ','E_pot  ','E_vdw  ','E_el  ', &
         'E_hbon  ','E_bon  ','E_angle  ','E_dih  ', &
         'E_14vdw  ','E_14el  ','E_const  ','E_pol  ', &
         'AV_permMoment ','AV_indMoment ','AV_totMoment ', &
         'Density', 'dV/dlambda', 'surften'/

   !     ----- DEFINE VARIOUS TERMS -----

   if (first) then
      !       -- up to Ekinetic:
      write(nf,1) 'L0 ', (labs(i),i=1,4)
      !       -- up to Pres_scal_solu:
      write(nf,1) 'L1 ', (labs(i),i=5,8)
      !       -- up to boxZ:
      write(nf,1) 'L2 ', (labs(i),i=9,12)
      !       -- up to pres_Z:
      write(nf,1) 'L3 ', (labs(i),i=13,16)
      !       -- up to EKCoM_z:
      write(nf,1) 'L4 ', (labs(i),i=17,20)
      !       -- up to VIRIAL_z:
      write(nf,1) 'L5 ', (labs(i),i=21,24)
      !       -- up to E_el:
      write(nf,1) 'L6 ', (labs(i),i=25,28)
      !       -- up to E_dih:
      write(nf,1) 'L7 ', (labs(i),i=29,32)
      !       -- up to E_pol:
      write(nf,1) 'L8 ', (labs(i),i=33,36)
      !       -- up to Density:
      write(nf,1) 'L9 ', (labs(i),i=37,41)
      1 format(a,10(1x,a))
      first = .false.
   end if

   !     ----- write values for this step -----

   !     -- up to Ekinetic:
   write(nf, 2) 'L0 ', nstep, time, ener%tot, ener%kin%tot
   !     -- up to Pres_scal_solu:
   write(nf, 3) 'L1 ', ener%kin%tot*onefac(1), ener%kin%solt*onefac(2), &
         ener%kin%solv*onefac(3), ener%kin%pres_scale_solt
   !     -- up to boxZ:
   write(nf, 3) 'L2 ', ener%kin%pres_scale_solv, ener%box(1), ener%box(2), ener%box(3)
   !     -- up to pres_Z:
   write(nf, 3) 'L3 ', ener%volume, ener%pres(1), ener%pres(2), ener%pres(3)
   !     -- up to EKCoM_z:
   write(nf, 3) 'L4 ', ener%pres(4), ener%cmt(1), ener%cmt(2), ener%cmt(3)
   !     -- up to VIRIAL_z:
   if( ntp > 0 ) then
      write(nf, 3) 'L5 ', ener%cmt(4), ener%vir(1), ener%vir(2), ener%vir(3)
   else
      write(nf, 3) 'L5 ', ener%cmt(4), ZERO, ZERO, ZERO
   end if
   !     -- up to E_el:
   write(nf, 3) 'L6 ', ener%vir(4), ener%pot%tot, ener%pot%vdw, ener%pot%elec
   !     -- up to E_dih:
   write(nf, 3) 'L7 ', ener%pot%hbond, ener%pot%bond, ener%pot%angle, ener%pot%dihedral
   !      -- up to E_pol:
   write(nf, 3) 'L8 ', ener%pot%vdw_14, ener%pot%elec_14, ener%pot%constraint, &
                                                   ener%pot%polar
   !     -- up to dV/dlambda:
   write(nf, 3) 'L9 ', ener%aveper, ener%aveind, ener%avetot, &
                                    ener%density, ener%pot%dvdl

   2 format(a, i8, 20(2x,e16.9))
   3 format(a, 20(e16.9,2x))
   return

    entry mdeng_reset
    first = .true.

end subroutine mdeng

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute the Center of Mass motion.
subroutine cenmas(natom,x,v,amass,ekcm,xcm,vcm,acm,ekrot,ocm,icm)
   implicit none
   integer:: i, icm, j, lh, m, mh, n, natom
   _REAL_ :: aamass, acm, amass, comvel, crit, d, ekcm, ekrot, ocm, &
        tcm, tmass, tmassinv, v, vcm, x, x1, x2, x3, xcm, &
        xx, xy, xz, yy, yz, zz
#  include "extra.h"

   !     ----- ROUTINE TO CALCULATE THE TRANSLATIONAL AND ROTATIONAL
   !           KINETIC ENERGIES AND VELOCITIES -----

   !     icm=0   return, do nothing
   !     icm=1   just return com coords in XCM
   !     icm=2   also calculate energy of com motion
   !     icm=3   also calcualte the angular momemtum
   !     icm=4   also calculate the inertia tensor
   !     icm=5   same as icm=4, but does not write anything to file - used by neb

   dimension x(*),v(*),amass(*),xcm(*),vcm(*),acm(*),ocm(*)
   dimension tcm(3,3),lh(3),mh(3)

   data crit/1.0d-06/

   if (icm == 0) return

   !     ----- CALCULATE THE CENTER OF MASS COORDINATES -----

   xcm(1) = 0.d0
   xcm(2) = 0.d0
   xcm(3) = 0.d0

   i = 0
   tmass = 0.d0
   do j = 1,natom
      aamass = amass(j)
      tmass = tmass + aamass
      xcm(1) = xcm(1) + x(i+1)*aamass
      xcm(2) = xcm(2) + x(i+2)*aamass
      xcm(3) = xcm(3) + x(i+3)*aamass
      i = i + 3
   end do
   tmassinv = 1.d0/tmass
   xcm(1) = xcm(1) * tmassinv
   xcm(2) = xcm(2) * tmassinv
   xcm(3) = xcm(3) * tmassinv

   if (iabs(icm) < 2) return

   !     ----- CALCULATE THE VELOCITY AND THE TRANSLATIONAL
   !           KINETIC ENERGY OF THE CENTRE OF MASS -----

   ekcm = 0.d0
   vcm(1) = 0.0d0
   vcm(2) = 0.0d0
   vcm(3) = 0.0d0

   i = 0
   do j = 1,natom
      aamass = amass(j)
      vcm(1) = vcm(1) + v(i+1)*aamass
      vcm(2) = vcm(2) + v(i+2)*aamass
      vcm(3) = vcm(3) + v(i+3)*aamass
      i = i + 3
   end do
   do i=1,3
      vcm(i) = vcm(i) * tmassinv
      ekcm = ekcm + vcm(i)*vcm(i)
   end do
   ekcm = ekcm * tmass * 0.5d0
   comvel = sqrt(vcm(1)*vcm(1)+vcm(2)*vcm(2)+vcm(3)*vcm(3))

   if(iabs(icm) < 3) return

   !     ----- CALCULATE THE ANGULAR MOMENTUM ABOUT THE
   !           CENTER OF MASS ----

   acm(1) = 0.0d0
   acm(2) = 0.0d0
   acm(3) = 0.0d0

   i = 0
   do j = 1,natom
      aamass = amass(j)
      acm(1) = acm(1) + (x(i+2)*v(i+3)-x(i+3)*v(i+2)) * aamass
      acm(2) = acm(2) + (x(i+3)*v(i+1)-x(i+1)*v(i+3)) * aamass
      acm(3) = acm(3) + (x(i+1)*v(i+2)-x(i+2)*v(i+1)) * aamass
      i = i+3
   end do

   acm(1) = acm(1) - (xcm(2)*vcm(3)-xcm(3)*vcm(2)) * tmass
   acm(2) = acm(2) - (xcm(3)*vcm(1)-xcm(1)*vcm(3)) * tmass
   acm(3) = acm(3) - (xcm(1)*vcm(2)-xcm(2)*vcm(1)) * tmass

   if (iabs(icm) < 4) return

   !     ----- CALCULATE THE INERTIA TENSOR -----

   xx = 0.d0
   xy = 0.d0
   xz = 0.d0
   yy = 0.d0
   yz = 0.d0
   zz = 0.d0

   i = 0
   do j = 1,natom
      x1 = x(i+1)-xcm(1)
      x2 = x(i+2)-xcm(2)
      x3 = x(i+3)-xcm(3)
      aamass = amass(j)
      xx = xx+x1*x1*aamass
      xy = xy+x1*x2*aamass
      xz = xz+x1*x3*aamass
      yy = yy+x2*x2*aamass
      yz = yz+x2*x3*aamass
      zz = zz+x3*x3*aamass
      i = i+3
   end do
   tcm(1,1) = yy+zz
   tcm(2,1) = -xy
   tcm(3,1) = -xz
   tcm(1,2) = -xy
   tcm(2,2) = xx+zz
   tcm(3,2) = -yz
   tcm(1,3) = -xz
   tcm(2,3) = -yz
   tcm(3,3) = xx+yy

   !     ----- INVERT THE INERTIA TENSOR -----

   call matinv(tcm,3,d,lh,mh)
   if(abs(d) <= crit) then
      write(6,307)
      307 format(/5x,'%CENMAS-F-INERTIA_TENSOR,  determinant', &
            ' is zero ... stop')
      call mexit(6,1)
   end if

   !     ----- CALCULATE THE ANGULAR VELOCITY ABOUT THE CENTER OF
   !           MASS AND THE ROTATIONAL KINETIC ENERGY -----

   ekrot = 0.d0
   do m = 1,3
      ocm(m) = 0.d0
      do n = 1,3
        ocm(m) = ocm(m)+tcm(m,n)*acm(n)
      end do
      ekrot = ekrot+ocm(m)*acm(m)
   end do
   ekrot = ekrot * 0.5d0

   if (icm < 0) return
   if (master.AND.icm.NE.5) then
      write(6,'(/3x,a,f11.4,3x,a,f11.4,3x,a,f12.6)') 'KE Trans =', &
            ekcm, 'KE Rot =', ekrot,'C.O.M. Vel =',comvel
      call flush(6)
   end if
   return
end subroutine cenmas

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Remove Center of Mass transrotational motion.
subroutine stopcm(nr,x,v,xcm,vcm,ocm, printmsg)

   implicit none
   integer nr
   _REAL_  x(*),v(*),xcm(*),vcm(*),ocm(*)
   logical printmsg !if true then a message is printed stating translational
                    !and rotational motion removed. Else no message is printed
                    !Useful for neb to stop too many messages being printed to screen.

#ifdef MPI
#  include "parallel.h"
#endif
#  include "extra.h"

   integer i, j, m
   _REAL_  x1, x2, x3

   !     ----- STOP THE CENTER OF MASS TRANSLATION -----

   i = 0
   do j = 1,nr
      do m = 1,3
         i = i+1
         v(i) = v(i)-vcm(m)
      end do
   end do

   !     ----- STOP THE ROTATION ABOUT THE CENTER OF MASS -----

   i = 0
   do j = 1,nr
      x1 = x(i+1)-xcm(1)
      x2 = x(i+2)-xcm(2)
      x3 = x(i+3)-xcm(3)
      v(i+1) = v(i+1)-ocm(2)*x3+ocm(3)*x2
      v(i+2) = v(i+2)-ocm(3)*x1+ocm(1)*x3
      v(i+3) = v(i+3)-ocm(1)*x2+ocm(2)*x1
      i = i+3
   end do
   if (master .AND. printmsg) write(6,9008)
   9008 format(/3x,'Translational and rotational motion removed')
   return
end subroutine stopcm
