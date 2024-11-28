! <compile=optimized>

! The 3D-RISM-KH software found here is copyright (c) 2010-2012 by
! Andriy Kovalenko, Tyler Luchko, Takeshi Yamazaki and David A. Case.
! 
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
! 
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
! 
! You should have received a copy of the GNU General Public License in the
! ../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
! 
! Users of the 3D-RISM capability found here are requested to acknowledge
! use of the software in reports and publications.  Such acknowledgement
! should include the following citations:
! 
! 1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999);
! ibid. 112:10391-10417 (2000).
! 
! 2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by
! F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.
! 
! 3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
! and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010).

#include "../include/dprec.fh"

!> TODO: Describe module.
!!
module rism1d_m
  use rism1d_c
  use rism_report_c
  use rism_timer_c
  implicit none
  integer, parameter :: CLEN = 256
  character(len = CLEN) :: fileroot
  type(rism1d), save :: rism
  type(rism_timer), save :: timer, ioTimer, inputTimer, outputTimer

  ! THEORY
  !> 1D-RISM theory to use.  May be DRISM or XRISM.
  character(len = CLEN) :: theory
  !> Closure to use. May be KH, PSE, HNC, PY or MV0.
  character(len = CLEN) :: closure
  !> Order parameter for closures like PSE.
  integer :: closure_order

  ! GRID SIZE
  !> Grid spacing [A].
  _REAL_ :: dr
  !> Number of grid points.
  integer :: nr

  ! OUTPUT
  !> List of output files to produce. See output() for details.
  character(len = CLEN) :: outlist
  
  !> Max r-space index to output.
  integer :: nrout
  !> Max k-space index to output.
  integer :: nkout
  !> Display iteration progress every 'progress' iterations.
  integer :: progress
  !> Write a restart file every ksave iterations.
  integer :: ksave
  !> 1 - run the default self-test.
  !! 0 - (default) do not run the self-test.
  !! -1 - run the full self-test, including properties we
  !!      don't expect to pass.
  integer :: selftest
  !> 1 - (default) output the Pettitt-Rossky form of exchem.
  !! 0 - do not output the Pettitt-Rossky form of exchem.
  integer exchem_pr
  !> 1 - output the Singer-Chandler form of exchem.
  !! 0 - (default) do not output the Singer-Chandler form of exchem.
  integer exchem_sc
  !> 1 - output the Schmeer-Maurer form of exchem.
  !! 0 - (default) do not output the Schmeer-Maurer form of exchem.
  integer exchem_sm

  !> Option to specify the use of extra precision in key pieces of code.
  integer :: extra_precision
  !> Perform a temperature derivative to get energy / entropy decomposition.
  integer :: entropicDecomp

  !> Max r-space distance to output.
  _REAL_ :: rout
  !> Max k-space wavenumber to output.
  _REAL_ :: kout

  ! SOLUTION CONVERGENCE
  !> Number of MDIIS vectors.
  integer :: mdiis_nvec
  !> Maximum number of iterations.
  integer :: maxstep
  !> MDIIS step size.
  _REAL_ :: mdiis_del
  !> Maximum tolerance for solution.
  _REAL_ :: tolerance
  !> Restart threshold factor. Ratio of the current residual to the
  !! minimum residual in the basis that causes a restart.
  _REAL_ :: mdiis_restart

  ! SOLVENT DESCRIPTION
  !> Temperature of the solvent [K].
  _REAL_ :: temperature
  !> Dielectric constant of the solvent.
  _REAL_ :: dieps

  ! ELECTROSTATICS
  !> Smear parameter for long range electrostatics [A].
  _REAL_ ::  smear
  !> Coefficient for DRISM (see Perkyns and Pettitt. J. Chem. Phys. 97, 1992, Eq. 34, variable 'a').
  _REAL_ :: adbcor


  ! DEPRECATED
  !> Synonym for closure.
  character(len = CLEN) :: closur
  !> Synonym for outlist.
  character(len = CLEN) :: outlst
  

  !> Synonym for rout.
  integer :: routup
  !> Synonym for kout.
  integer :: toutup
  !> Synonym for progress.
  integer :: kshow
  !> Synonym for mdiis_nvec.
  integer :: nis
  !> Synonym for maxstep.
  integer :: maxste
  !> Synonym for mdiis_del.
  _REAL_ :: delvv
  !> Synonym for tolerance.
  _REAL_ :: tolvv
  !> Synonym for temperature.
  _REAL_ :: temper

contains

  !> Reads user input and initialize rism1d object.
  subroutine initialize()
    use rism_util, only : freeUnit
    use array_util, only : array_index
    use constants, only : avogadro
    use solvMDL_c
    use rism1d_nub_c, only : rism1d_nub_initialize
    implicit none
#include "../xblas/f77/blas_namedconstants.fh"

    character(len = clen) :: inpfile
    integer :: unit, i, iostat
    type(solvMDL) :: mdl

    !nsp :: number of solvent species
    !ncoeff :: number of closure coefficients
    integer :: nsp, ncoeff
    integer, parameter :: ncoeff_buffer = 30
    !closure_coeff_p :: coefficients for closures.  This is the
    !pointer that will be passed but a static array used for
    !the namelist
    _REAL_, pointer :: closure_coeff_p(:) => NULL()
    _REAL_ :: closure_coeff(ncoeff_buffer)
    !density :: density of the solvent species
    !chg_scale :: scale all charges in MDL by this amount
    _REAL_ :: density, chg_scale = 1d0
    !model :: MDL file
    character(len = CLEN) :: model, units

    namelist /parameters/ outlist, theory, closure, nr, dr, rout, kout, &
         mdiis_nvec, mdiis_del, mdiis_restart, tolerance, ksave, progress, &
         maxstep, smear, adbcor, temperature, exchem_pr, exchem_sc, exchem_sm, selftest, &
         dieps, nsp, closure_order, closure_coeff, extra_precision, entropicDecomp, &
         closur, outlst, routup, toutup, kshow, nis, maxste, delvv, tolvv, temper
    namelist /species/ density, model, units, chg_scale

    call rism_timer_new(timer, "Total")
    call rism_timer_start(timer)
    call rism_timer_new(ioTimer, "I/O")
    call rism_timer_setParent(ioTimer, timer)
    call rism_timer_new(inputTimer, "Input")
    call rism_timer_setParent(inputTimer, ioTimer)
    call rism_timer_new(outputTimer, "Output")
    call rism_timer_setParent(outputTimer, ioTimer)
    call defaults(closure_coeff)

    call rism_timer_start(inputTimer)
    call  getarg (1, fileroot)
    inpfile = trim(fileroot) // ".inp"

    call rism_report_message('reading input data file: ' // trim(inpfile))
    unit = freeUnit()
    open(unit = unit, file = trim(inpfile), status = 'old', iostat = iostat)
    if (iostat /=  0) &
         call rism_report_error('(a, i4)', "Could not open " // trim(inpfile) // ":", iostat)

    read (unit, parameters)
    call rism_timer_stop(inputTimer)

    call sanity_check()
    !remove buffer values from coefficient array and transfer to a
    !properly sized array
    ncoeff = array_index(closure_coeff, huge(1d0))-1
    if (ncoeff<0) ncoeff = ubound(closure_coeff, 1)
    closure_coeff_p => safemem_realloc(closure_coeff_p, ncoeff)
    if (ncoeff >0)&
         closure_coeff_p = closure_coeff(1:ncoeff)
    call rism1d_new(rism, theory, closure, closure_coeff_p, &
         temperature, dieps, smear, adbcor, nr, dr, mdiis_nvec, mdiis_del, &
         mdiis_restart, trim(fileroot) // '.sav', extra_precision)
    if (rout == 0) then
       nrout = nr
    else
       nrout = min(nr, nint(rout / dr) + 1)
    end if
    if (kout == 0) then
       nkout = nr
    else
       nkout = min(nr, int(kout / rism%pot%dk))
    end if

    !---read paramters for each species of solvent---!
    do i = 1, nsp
       call rism_timer_start(inputTimer)
       !default units
       units = 'M'
       read (unit, species)
       call rism_timer_stop(inputTimer)
       call solvMDL_new(mdl, model)
       !the precision of the density has implications for the
       !precision of calculating Ak in single3DRISMsolution().  Different
       !compilers (e.g. GNU and Intel) have subtle differences in
       !the unit conversion that cause small numerical discrepancies
       !in the final solution. Using XBLAS removes these small
       !differences when using different compilers on the same
       !hardware.
       select case (units)
       case ('M')
          if (extra_precision < 1) then
             density = density*avogadro*1d-27
          else
             call BLAS_DAXPBY_X(1, 0d0, 0d0, 1, avogadro*1d-27, density, 1, BLAS_PREC_EXTRA)
          end if
       case ('mM')
          if (extra_precision < 1) then
             density = density*avogadro*1d-30
          else
             call BLAS_DAXPBY_X(1, 0d0, 0d0, 1, avogadro*1d-30, density, 1, BLAS_PREC_EXTRA)
          end if
       case ('1/A^3')
       case ('g/cm^3')
          if (extra_precision < 1) then
             density = density / sum(mdl%mass*mdl%multiplicity)*avogadro*1d-24
          else
             call BLAS_DAXPBY_X(1, 0d0, 0d0, 1, avogadro*1d-24 / dble(sum(mdl%mass*mdl%multiplicity)), &
                  density, 1, BLAS_PREC_EXTRA)
          end if
       case ('kg/m^3')
          if (extra_precision < 1) then
             density = density / sum(mdl%mass*mdl%multiplicity)*avogadro*1d-27
          else
             call BLAS_DAXPBY_X(1, 0d0, 0d0, 1, avogadro*1d-27 / dble(sum(mdl%mass*mdl%multiplicity)), &
                  density, 1, BLAS_PREC_EXTRA)
          end if
       case default
          call rism_report_error("'"//trim(units)//"' are not valid density units")
       end select
       mdl%charge = mdl%charge*chg_scale
       call rism1d_addSpecies(rism, mdl, density)

       call solvMDL_destroy(mdl)
    end do

    ! close the input file
    close(unit)

    ! read Non-bond modification
    open(unit = unit, file = trim(fileroot)//'.nbfix', status = 'old', iostat = iostat)
    if (iostat == 0) then
       close(unit)
       call rism1d_readNBFix(rism, trim(fileroot)//'.nbfix')
    end if

    ! initialize some rism_closure objects
    if (rism%closure%type .eq. "NUB")&
         call rism1d_nub_initialize(rism%closure%nub, trim(fileroot)//'.nub', &
         rism%pot)

    if (safemem_dealloc(closure_coeff_p) /= 0)&
         call rism_report_error("rism1d: failed to deallocate closure_coeff_p")
  end subroutine initialize


  !> Set default user parameters.
  subroutine defaults(closure_coeff)
    implicit none
    _REAL_, intent(out) :: closure_coeff(:)
    theory = 'DRISM'
    closure = 'KH'
    outlist = ''
    extra_precision = 1
    entropicDecomp = 1
    nr = 2**14
    dr = 0.025d0
    rout = 0
    kout = 0
    mdiis_nvec = 20
    mdiis_del = 0.3d0
    mdiis_restart = 10d0
    tolerance = 1d-12
    ksave = -2
    progress = 1
    maxstep = 10000
    smear = 1d0
    adbcor = 0.5d0
    temperature = 298.15d0
    ! Must be set by user if theory = 'drism'.
    dieps = -1d0             
    closure_order = 3
    ! Place holder values
    closure_coeff = huge(0d0)        
    exchem_pr = 1
    exchem_sc = 0
    exchem_sm = 0
    selftest = 0
    ! Deprecated variable names.
    closur = ''
    outlst = ''
    routup = -1
    toutup = -1
    kshow = -1
    nis   = -1
    maxste = -1
    delvv = -1d0
    tolvv = -1d0
    temper = -1d0
  end subroutine defaults


  !> Check user parameters.
  subroutine sanity_check()
    use rism_util, only : caseup
    implicit none
    character(len = 32) :: fmt
    !set everything to upper case
    call  caseup (outlist)
    call  caseup (outlst)
    call  caseup (theory)
    call  caseup (closure)
    call  caseup (closur)

    !handle deprecated parameter names
    if (closur .ne. '') then
       call rism_report_warn("'closur' is deprecated")
       if (closure .ne. 'KH') call rism_report_error("Both 'closure' and 'closur' defined")
       closure = closur
    end if
    if (outlst .ne. '') then
       call rism_report_warn("'outlst' is deprecated")
       if (outlist .ne. '') call rism_report_error("Both 'outlist' and 'outlst' defined")
       outlist = outlst
    end if
    if (routup /= -1) then
       call rism_report_warn("'routup' is deprecated")
       if (rout /=  0) call rism_report_error("Both 'rout' and 'routup' defined")
       rout = routup
    end if
    if (toutup /= -1) then
       call rism_report_warn("'toutup' is deprecated")
       if (kout /=  0) call rism_report_error("Both 'kout' and 'toutup' defined")
       kout = toutup
    end if
    if (kshow /= -1) then
       call rism_report_warn("'kshow' is deprecated")
       if (progress /=  1) call rism_report_error("Both 'progress' and 'kshow' defined")
       progress = kshow
    end if
    if (nis /= -1) then
       call rism_report_warn("'nis' is deprecated")
       if (mdiis_nvec /=  20) call rism_report_error("Both 'mdiis_nvec' and 'nis' defined")
       mdiis_nvec = nis
    end if
    if (maxste /= -1) then
       call rism_report_warn("'maxste' is deprecated")
       if (maxstep /=  10000) call rism_report_error("Both 'maxstep' and 'maxste' defined")
       maxstep = maxste
    end if
    if (delvv /= -1d0) then
       call rism_report_warn("'delvv' is deprecated")
       if (mdiis_del /=  0.3d0) call rism_report_error("Both 'mdiis_del' and 'delvv' defined")
       mdiis_del = delvv
    end if
    if (tolvv /= -1d0) then
       call rism_report_warn("'tolvv' is deprecated")
       if (tolerance /=  1d-12) call rism_report_error("Both 'tolerance' and 'tolvv' defined")
       tolerance = tolvv
    end if
    if (temper /= -1d0) then
       call rism_report_warn("'temper' is deprecated")
       if (temperature /=  298.15d0) call rism_report_error("Both 'temperature' and 'temper' defined")
       temperature = temper
    end if

    !outlist: hmmm... leave this for later...

    if (rout < 0) &
         call rism_report_error("ROUT must be >= 0")
    if (kout < 0) &
         call rism_report_error("KOUT must be >= 0")

    !Check the dielectric constant
    if (theory .eq. "DRISM" .and. dieps < 0) &
         call rism_report_error("DIEPS must be set for DRISM")

    !optional input method for closures
    if (trim(closure).eq."PSEN" .or. trim(closure).eq."PSE") then
       write(fmt, '(a, i4, a)') '(a, i', int(log10(dble(closure_order))) + 1, ')'
       write(closure, fmt) "PSE", closure_order
    end if
  end subroutine sanity_check


!!!Invokes the rism1d solver

  subroutine solve()
    implicit none
    logical :: converged
    converged = rism1d_solve(rism, ksave, progress, maxstep, tolerance)
    if (.not. converged) then
       call rism_report_error('(a, i5)', '1D-RISM failed to converge.  Reached steps limit Maxstep = ', maxstep)
    end if

    !Temperature derivatives
    if (closure .ne. "KH" .and. closure .ne. "HNC" .and. index(closure, "PSE") /= 1) then
       call rism_report_message("Temperature derivatives not supported for " // closure)
       entropicDecomp = 0
       return
    end if
    if (entropicDecomp == 1) then
       converged = rism1d_dt_solve(rism, ksave, kshow, maxstep, tolerance)
       if (.not. converged) then
          call rism_report_error('(a, i5)', '1D-RISM DT failed to converge.  Reached steps limit MaxStep = ', maxstep)
       end if
    end if
  end subroutine solve


!!!Cleanup memory

  subroutine finalize()
    implicit none
    integer*8 :: memstats(10)
    integer :: unit
    unit = rism_report_getMUnit()
    call rism_timer_summary(timer)
    call rism1d_destroy(rism)
    call rism_timer_destroy(inputTimer)
    call rism_timer_destroy(outputTimer)
    call rism_timer_destroy(ioTimer)
    call rism_timer_destroy(timer)
    memstats = memStatus()
    write(unit, '(a)') "1D-RISM memory allocation summary"
    write(unit, '(a)') "Type         Current         Maximum"
    write(unit, '(a, i12, a,f12.5, a)') "Integer  ", memstats(1), " B ", &
         dble(memstats(6)) / BYTES_PER_GB, " GB"
    write(unit, '(a, i12, a,f12.5, a)') "Real     ", memstats(2), " B ", &
         dble(memstats(7)) / BYTES_PER_GB, " GB"
    write(unit, '(a, i12, a,f12.5, a)') "Logical  ", memstats(3), " B ", &
         dble(memstats(8)) / BYTES_PER_GB, " GB"
    write(unit, '(a, i12, a,f12.5, a)') "Character", memstats(4), " B ", &
         dble(memstats(9)) / BYTES_PER_GB, " GB"
    write(unit, '(a)') "---------------------------------------"
    write(unit, '(a, i12, a,f12.5, a)') "Total    ", memstats(5), " B ", &
         dble(memstats(10)) / BYTES_PER_GB, " GB"

  end subroutine finalize


!!!Calculate and output the thermodynamic properties the user has requested

  subroutine output()
    use constants, only : KB, COULOMB_CONST_E
    implicit none
    _REAL_, pointer :: exnvv(:, :,:) => NULL(), nvv(:, :,:) => NULL(), &
         bvv(:, :) => NULL()
    _REAL_ :: qEx(rism%pot%nv, rism%pot%nv)
    integer :: iv
    ! NOTE: We should be able to combine writeUVV and writeRunning Coordination
    ! into something more general.
    if (index(outlist, 'U') /=  0)  then
       call writeVV(rism%pot%uvv, 'uvv', 'POTENTIAL ENERGY [kT]')
    end if
    if (index(outlist, 'X') /=  0)  then
       call writeXVV()
    end if
    if (index(outlist, 'G') /=  0)  then
       call writeVV(rism%gvv, 'gvv', 'PAIR DISTRIBUTION')
       if (entropicDecomp == 1)&
            call writeVV(rism%gvv_dT, 'gvv_dT', 'DT PAIR DISTRIBUTION')
    end if
    if (index(outlist, 'H') /=  0)  then
       call writeVV(rism%hvv, 'hvv', 'TOTAL CORRELATION', k_space = .true.)
       if (entropicDecomp == 1)&
            call writeVV(rism%hvv_dT, 'hvv_dT', 'DT TOTAL CORRELATION', k_space = .true.)
    end if
    if (index(outlist, 'C') /=  0)  then
       call writeVV(rism%cvv, 'cvv', 'DIRECT CORRELATION')
       if (entropicDecomp == 1)&
            call writeVV(rism%cvv_dT, 'cvv_dT', 'DT DIRECT CORRELATION')
    end if
    if (index(outlist, 'B') /=  0)  then
       bvv => rism1d_bvv(rism)
       call writeVV(bvv, 'bvv', 'BRIDGE FUNCTION')
       if (safemem_dealloc(bvv) /= 0)&
            call rism_report_error("OUTPUT: dealloc of Bvv failed")
    end if
    if (index(outlist, 'T') /=  0)  then
       ! call writeTD()
       call writeTherm()
    end if
    if (index(outlist, 'E') /=  0)  then
       exnvv => rism1d_getRunExNumber(rism)
       call writeRunningCoordination(exnvv, 'exnvv', "RUNNING EXCESS COORDINATION NUMBER")
       if (safemem_dealloc(exnvv) /= 0) &
            call rism_report_error("Deallocating memory for running excess number")
       call writeTotalExcess(rism1d_getExNumber(rism), 'n00', &
            "Total excess coordination number")
    end if
    if (index(outlist, 'Q') /=  0)  then
       qEx = rism1d_getExNumber(rism)
       do iv =1, rism%pot%nv
          qEx(iv, :) = qEx(iv, :)*rism%pot%qv*rism%pot%mtv*sqrt(KB*rism%pot%temperature / COULOMB_CONST_E)
       end do
       call writeTotalExcess(qEx, 'q00', &
            "Total excess coordinated charge [e]", .true.)
    end if
    if (index(outlist, 'N') /=  0)  then
       nvv => rism1d_getRunNumber(rism)
       call writeRunningCoordination(nvv, 'nvv', "RUNNING COORDINATION NUMBER")
       if (safemem_dealloc(nvv) /= 0)&
            call rism_report_error("Deallocating memory for running number")
    end if
    if (index(outlist, 'S') /=  0)  then
       call writeSVV()
    end if
    if (selftest == 1) then
       call rism1d_selftest(rism, trim(fileroot)//'.self.test')
    else if (selftest == -1) then
       call rism1d_selftest(rism, trim(fileroot)//'.self.test', o_all = .true.)
    end if
  end subroutine output


!!!Writes the site-site suceptibility to file

  subroutine writeXVV()
    use rism_util
    use constants, only : KB, COULOMB_CONST_E, JPKC
    implicit none
    character(len = clen) :: file, metafmt, intfmt, real3fmt, real5fmt, textfmt
    integer :: unit, iostat
    character(len = 8) :: date
    character(len = 10) :: time
    character(len = 5) :: zone
    integer :: values(8)
    integer :: iv, ir, iv1, iv2, isp, imlt
    _REAL_ :: delhv0(rism%pot%nv), delhv0_dT(rism%pot%nv), background_correction(rism%pot%nv), &
         kappa, compressibility, compressibility_dT
    _REAL_, pointer :: xvv(:, :, :) => NULL(), xvv_dT(:, :, :) => NULL(), coord(:) => NULL()

    write(metafmt, '(a)') '(A)'
    write(intfmt, '(a)') '(10I8)'
    write(real3fmt, '(a)') '(1P3E24.16)'
    write(real5fmt, '(a)') '(1P5E24.16)'
    write(textfmt, '(a)') '(20A4)'

    background_correction = rism1d_get3DRISMbackground(rism)
    delhv0 = rism1d_getDelHvLimit(rism)
    delhv0_dT = rism1d_getDelHvLimit_DT(rism)
    kappa = rism1d_getInvDebyeLen(rism)
    compressibility = rism1d_getCompressibility(rism)
    if (entropicDecomp == 1) then
       compressibility_dT = rism1d_getCompressibility_dT(rism)
    end if
    xvv => rism1d_getSusceptibility(rism)
    xvv_dT => rism1d_getSusceptibility_DT(rism)
    call rism_timer_start(outputTimer)

    !count number of atoms
    iv = 0
    do isp = 1, rism%POT%NSP
       do iv1 = 1, rism%pot%nat(isp)
          do imlt = 1, rism%pot%mta(iv1, isp)
             iv = iv + 1
          end do
       end do
    end do
    !allocate memory
    coord => safemem_realloc(coord, iv*3)
    !copy over coordinates
    iv = 0
    do isp = 1, rism%POT%NSP
       do iv1 = 1, rism%pot%nat(isp)
          do imlt = 1, rism%pot%mta(iv1, isp)
             coord(3*iv + 1:3*(iv + 1)) = rism%pot%rma(:, imlt, iv1, isp)
             iv = iv + 1
          end do
       end do
    end do

    unit = freeUnit()
    file = trim(fileroot)//'.xvv'
    call rism_report_message('outputting Xvv(K) to file: '//trim(file))
    open(unit, file = file, status = 'replace', iostat = iostat)
    if (iostat /= 0)&
         call rism_report_error("failed to open :"//trim(file))

    call date_and_time(date, time, zone, values)
    write (unit, '(12a)') "%VERSION  VERSION_STAMP = V0001.001  DATE = ", &
         date(5:6), ":", date(7:8), ":", date(3:4), " ", time(1:2), ":", time(3:4), ":", time(5:6)
    write (unit, metafmt) "%COMMENT NR, NV, NSP"
    write (unit, metafmt) "%FLAG POINTERS"
    write (unit, metafmt) "%FORMAT"//trim(intfmt)
    write (unit, intfmt) rism%POT%NR, rism%POT%NV, rism%POT%NSP
    write (unit, metafmt) "%COMMENT TEMPERATURE [K], DIEPS, KAPPA [1/A], COMPRESSIBILITY [A^3], DR [A], SMEAR [A]"
    write (unit, metafmt) "%FLAG THERMO"
    write (unit, metafmt) "%FORMAT"//trim(real5fmt)
    write (unit, real5fmt) rism%POT%TEMPERATURE, rism%pot%dielconst, kappa, &
         compressibility, rism%POT%DR, rism%POT%SMEAR
    if (entropicDecomp == 1) then
       write (unit, metafmt) "%COMMENT COMPRESSIBILITY_DT [A^3/K]"
       write (unit, metafmt) "%FLAG THERMO_DT"
       write (unit, metafmt) "%FORMAT"//trim(real5fmt)
       write (unit, real5fmt) compressibility_dT
    end if
    write (unit, metafmt) "%FLAG ATOM_NAME"
    write (unit, metafmt) "%FORMAT"//trim(textfmt)
    write (unit, textfmt) rism%pot%namev
    write (unit, metafmt) "%COMMENT SITE MULTIPLICITY"
    write (unit, metafmt) "%FLAG MTV"
    write (unit, metafmt) "%FORMAT"//trim(intfmt)
    write (unit, intfmt) rism%pot%mtv
    write (unit, metafmt) "%COMMENT NUMBER OF SITES FOR EACH SPECIES"
    write (unit, metafmt) "%FLAG NVSP"
    write (unit, metafmt) "%FORMAT"//trim(intfmt)
    write (unit, intfmt) rism%pot%nat
    write (unit, metafmt) "%COMMENT MASS [g/mol]"
    write (unit, metafmt) "%FLAG MASS"
    write (unit, metafmt) "%FORMAT"//trim(real5fmt)
    write (unit, real5fmt) rism%pot%mass
    write (unit, metafmt) "%COMMENT RHOV [A^{-3}]"
    write (unit, metafmt) "%FLAG RHOV"
    write (unit, metafmt) "%FORMAT"//trim(real5fmt)
    write (unit, real5fmt) rism%pot%densityv
    write (unit, metafmt) "%COMMENT RHOSP [A^{-3}]"
    write (unit, metafmt) "%FLAG RHOSP"
    write (unit, metafmt) "%FORMAT"//trim(real5fmt)
    write (unit, real5fmt) rism%pot%densitysp
    write (unit, metafmt) "%COMMENT QV [sqrt(kT A)]"
    write (unit, metafmt) "%FLAG QV"
    write (unit, metafmt) "%FORMAT"//trim(real5fmt)
    write (unit, real5fmt) rism%pot%qv
    write (unit, metafmt) "%COMMENT QSPV [sqrt(kT A)]"
    write (unit, metafmt) "%FLAG QSPV"
    write (unit, metafmt) "%FORMAT"//trim(real5fmt)
    write (unit, real5fmt) rism%pot%qspv
    write (unit, metafmt) "%COMMENT EPSV [kT]"
    write (unit, metafmt) "%FLAG EPSV"
    write (unit, metafmt) "%FORMAT"//trim(real5fmt)
    write (unit, real5fmt) rism%pot%ljEpsilonVV
    write (unit, metafmt) "%COMMENT RMIN2V [A]"
    write (unit, metafmt) "%FLAG RMIN2V"
    write (unit, metafmt) "%FORMAT"//trim(real5fmt)
    write (unit, real5fmt) rism%pot%ljRminVV
    write (unit, metafmt) "%COMMENT 3D-RISM h(k=0) BACKGROUND CHARGE CORRECTION [sqrt(kT A)]"
    write (unit, metafmt) "%FLAG BACKGROUND_CORRECTION"
    write (unit, metafmt) "%FORMAT"//trim(real5fmt)
    write (unit, real5fmt) background_correction
    write (unit, metafmt) "%COMMENT DELHV0 [sqrt(kT A)]"
    write (unit, metafmt) "%FLAG DELHV0"
    write (unit, metafmt) "%FORMAT"//trim(real5fmt)
    write (unit, real5fmt) delhv0
    if (entropicDecomp == 1) then
       write (unit, metafmt) "%COMMENT DELHV0_DT [sqrt(kT A)]"
       write (unit, metafmt) "%FLAG DELHV0_DT"
       write (unit, metafmt) "%FORMAT"//trim(real5fmt)
       write (unit, real5fmt) (delhv0_dT)
    end if
    write (unit, metafmt) "%COMMENT COORD [A]"
    write (unit, metafmt) "%FLAG COORD"
    write (unit, metafmt) "%FORMAT"//trim(real3fmt)
    write (unit, real3fmt) coord
    write (unit, metafmt) "%COMMENT COLUMN MAJOR NR X NAT X NAT"
    write (unit, metafmt) "%FLAG XVV"
    write (unit, metafmt) "%FORMAT"//trim(real5fmt)
    write (unit, real5fmt)  xvv
    if (entropicDecomp == 1) then
       write (unit, metafmt) "%COMMENT COLUMN MAJOR NR X NAT X NAT"
       write (unit, metafmt) "%FLAG XVV_DT"
       write (unit, metafmt) "%FORMAT"//trim(real5fmt)
       write (unit, real5fmt)  xvv_dT
    end if
    close (unit)
    close(unit)
    if (safemem_dealloc(coord) /= 0)&
         call rism_report_error("WRITEXVV: : failed to deallocate COORD")
    call rism_timer_stop(outputTimer)
  end subroutine writeXVV


!!!Writes the r-space site-site correlation function to file
!!!IN:
!!!   vv    : correlation function (nr, nvv)
!!!   suffix: suffix for the file
!!!   descript : description of the column
!!!   k_space : (optional) if true, data is in reciprocal space and the frequency
!!!            will be printed out.  if false (default), data is in real space
!!!            and the separation will be printed out

  subroutine writeVV(vv, suffix, descript, k_space)
    use rism_util
    implicit none
    _REAL_, intent(in) :: vv(:, :)
    character(len = *), intent(in) :: suffix, descript
    logical, optional, intent(in) :: k_space
    logical :: kspace
    character(len = clen) :: file, usuffix
    integer :: unit, iostat
    character(len = 64) :: fmt
    integer :: ir, ivv
    call rism_timer_start(outputTimer)
    kspace = .false.
    if (present(k_space)) kspace = k_space
    unit = freeUnit()
    file = trim(fileroot)//'.'//suffix
    usuffix = suffix
    call caseup(usuffix(1:1))
    call rism_report_message('outputting '//trim(usuffix)//'(R) to file: '//trim(file))
    open(unit, file = file, status = 'replace', iostat = iostat)
    if (iostat /= 0)&
         call rism_report_error("failed to open :"//trim(file))
    if (kspace) then
       call writeVVHeader(unit, &
            "#RISM1D ATOM-ATOM INTERACTIONS: "//trim(descript)//" VS. FREQUENCY [1/A]", &
            17, k_space = kspace)
    else
       call writeVVHeader(unit, &
            "#RISM1D ATOM-ATOM INTERACTIONS: "//trim(descript)//" VS. SEPARATION [A]", 17)
    end if
    write(fmt, '(a, i8, a)') "(1p, 1x, ", rism%pot%nvv + 1, "(1x, E16.8E3))"
    !      write(fmt, '(a, i8, a)') "(1p, 1x, ", rism%pot%nvv + 1, "(1x, E24.16E3))"
    if (kspace) then
       do ir = 1, nkout
          write (unit, fmt)  (ir-1)*rism%pot%dk, rmExPrec(vv(ir, :))
       end do
    else
       do ir = 1, nrout
          write (unit, fmt)  (ir-1)*rism%pot%dr, rmExPrec(vv(ir, :))
       end do
    end if
    close(unit)
    call rism_timer_stop(outputTimer)
  end subroutine writeVV

  
  !> Writes bulk thermodynamic quantities to file.
  subroutine writeTherm()
    use rism_util, only : freeUnit
    use constants, only : boltzmann, KB
    implicit none
    _REAL_ :: pmv(rism%pot%nsp)
    !we need exchem for each type form of the equation (PR, SM or SC)
    _REAL_ :: exchemPR(rism%pot%nv), exchemspPR(rism%pot%nsp), &
         exchemSC(rism%pot%nv), exchemspSC(rism%pot%nsp), &
         exchemSM(rism%pot%nv), exchemspSM(rism%pot%nsp)
    _REAL_ :: solvene(rism%pot%nv), solvenesp(rism%pot%nsp)
    _REAL_ :: compressibility, compressibility_dT, pressureFE, freeEnergy
    integer :: isp, iv, iat
    character(len = clen) :: file
    integer :: unit, iostat
    !descriptFmt : Format for category (calculation type, e.g. excess chemical potential)
    !varFmt      : Format for variable name
    !unitFmt     : Format for units
    !titleFmt    : Format for column headings
    !valFmt      : Format for floating point values
    character(len = 16) :: descriptFmt = "(a34)", varFMT = "(a12)", unitFmt = "(a15)", &
         titleFmt = "(a17)", valFmt = '(1p,e16.8e3,1x)'
    !whtspc : long string of whitespace that can be used to effect a
    !         left-justified string.  Otherwise strings are right-justified.
    !         Simply concatenate this to the end of the string you wish
    !         left-justified.
    character(len = 34) :: whtspc


    !pre-calculate some values

    exchemPR = huge(1d0)
    exchemspPR = huge(1d0)
    exchemSC = huge(1d0)
    exchemspSC = huge(1d0)
    exchemSM = huge(1d0)
    exchemspSM = huge(1d0)

    if (exchem_PR == 1) then
       exchemPR = rism1d_getExChem(rism, 'PR')
       exchemspPR = 0d0
    end if
    if (exchem_SC == 1) then
       exchemSC = rism1d_getExChem(rism, 'SC')
       exchemspSC = 0d0
    end if
    if (exchem_SM == 1) then
       exchemSM = rism1d_getExChem(rism, 'SM')
       exchemspSM = 0d0
    end if

    if (entropicDecomp == 1) then
       solvene = rism1d_getSolvene(rism)
    else
       solvene = huge(1d0)
    end if
    solvenesp = 0
    iv = 0
    do isp = 1, rism%pot%nsp
       do iat = 1, rism%pot%nat(isp)
          iv = iv + 1
          if (exchemPR(iv) /= huge(1d0))&
               exchemspPR(isp) = exchemspPR(isp) + exchemPR(iv)
          if (exchemSC(iv) /= huge(1d0))&
               exchemspSC(isp) = exchemspSC(isp) + exchemSC(iv)
          if (exchemSM(iv) /= huge(1d0))&
               exchemspSM(isp) = exchemspSM(isp) + exchemSM(iv)
          solvenesp(isp) = solvenesp(isp) + solvene(iv)
       end do
    end do
    pmv = rism1d_getPMV(rism)
    compressibility = rism1d_getCompressibility(rism)
    if (entropicDecomp == 1) then
       compressibility_dT = rism1d_getCompressibility_dT(rism)&
            /rism%POT%TEMPERATURE
    end if
    pressureFE = rism1d_getPressureFE(rism)
    freeEnergy = rism1d_getFreeEnergy(rism)

    ! Open the file.
    write(whtspc, descriptFmt) " "
    unit = freeUnit()
    file = trim(fileroot)//'.therm'
    call rism_report_message('outputting thermodynamics to file: '//trim(file))
    open(unit, file = file, status = 'replace', iostat = iostat)

    ! GLOBAL PROPERTIES
    write(unit, descriptFmt) "#Global properties"//whtspc
    write(unit, descriptFmt, advance = 'no') "#Description"//whtspc
    write(unit, varFmt, advance = 'no') "Variable"//whtspc
    write(unit, unitFmt, advance = 'no') "Units"//whtspc
    write(unit, titleFmt) "Value"//whtspc

    ! COMPRESSIBILITY
    write(unit, descriptFmt, advance = 'no') "Compressibility"//whtspc
    write(unit, varFmt, advance = 'no') "xi"//whtspc
    write(unit, unitFmt, advance = 'no') "[1/MPa]"//whtspc
    write(unit, valFmt) compressibility&
         *1d-24 / (boltzmann*rism%pot%temperature)

    ! COMPRESSIBILITY dT
    if (entropicDecomp == 1) then
       write(unit, descriptFmt, advance = 'no') "Compressibility_dT"//whtspc
       write(unit, varFmt, advance = 'no') "xi_dT"//whtspc
       write(unit, unitFmt, advance = 'no') "[1/MPa/K]"//whtspc
       write(unit, valFmt) compressibility_dT&
            *1d-24 / (boltzmann*rism%pot%temperature)
    end if

    ! PRESSURE - FREE ENERGY
    if (pressureFE /= HUGE(1d0)) then
       write(unit, descriptFmt, advance = 'no') "Pressure_(Free_Energy)"//whtspc
       write(unit, varFmt, advance = 'no') "Pfe"//whtspc
       write(unit, unitFmt, advance = 'no') "[MPa]"//whtspc
       write(unit, valFmt) pressureFE&
            *1.d24 * boltzmann*rism%pot%temperature
    end if

    ! PRESSURE - VIRIAL
    ! Only output the virial pressure for pure atomic fluids where
    ! pressureFE cannot be computed.
    if (pressureFE == HUGE(1d0) .and. sum(rism%pot%mtv) == rism%pot%nv .and. &
         sum(rism%pot%nat) == rism%pot%nsp) then
       write(unit, descriptFmt, advance = 'no') "Pressure_(Virial)"//whtspc
       write(unit, varFmt, advance = 'no') "Pvir"//whtspc
       write(unit, unitFmt, advance = 'no') "[MPa]"//whtspc
       write(unit, valFmt) rism1d_getPressureVirial(rism)&
            *1.d24 * boltzmann*rism%pot%temperature
    end if

    ! TOTAL FREE ENERGY
    if (freeEnergy /= HUGE(1d0)) then
       write(unit, descriptFmt, advance = 'no') "Excess_free_energy"//whtspc
       write(unit, varFmt, advance = 'no') "FE"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol/A^3]"//whtspc
       write(unit, valFmt) freeEnergy&
            *KB*rism%pot%temperature
    end if

    !Species
    write(unit, '(a)')
    write(unit, descriptFmt) "#Species properties"//whtspc
    write(unit, descriptFmt, advance = 'no') "#Description"//whtspc
    write(unit, varFmt, advance = 'no') "Variable"//whtspc
    write(unit, unitFmt, advance = 'no') "Units"//whtspc
    do isp = 1, rism%pot%nsp
       write(unit, titleFmt, advance = 'no') rism%pot%namesp(isp)//whtspc
    end do
    write(unit, '(a)')

    ! EXCESS CHEMICAL POTENTIAL
    ! Test site value as species value is a sum.
    if (exchemPR(1) /= HUGE(1d0)) then
       write(unit, descriptFmt, advance = 'no') "Excess_chemical_potential_PR"//whtspc
       write(unit, varFmt, advance = 'no') "EXCHEMsp_PR"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do isp = 1, rism%pot%nsp
          write(unit, valFmt, advance = 'no') exchemspPR(isp)*KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if
    if (exchemSC(1) /= HUGE(1d0)) then
       write(unit, descriptFmt, advance = 'no') "Excess_chemical_potential_SC"//whtspc
       write(unit, varFmt, advance = 'no') "EXCHEMsp_SC"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do isp = 1, rism%pot%nsp
          write(unit, valFmt, advance = 'no') exchemspSC(isp)*KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if
    if (exchemSC(1) /= HUGE(1d0)) then
       write(unit, descriptFmt, advance = 'no') "Excess_chemical_potential_SM"//whtspc
       write(unit, varFmt, advance = 'no') "EXCHEMsp_SM"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do isp = 1, rism%pot%nsp
          write(unit, valFmt, advance = 'no') exchemspSM(isp)*KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if

    ! SOLVATION ENERGY
    ! Test site value as species value is a sum.
    if (solvene(1) /= HUGE(1d0) .and. entropicDecomp == 1) then
       write(unit, descriptFmt, advance = 'no') "Solvation_energy"//whtspc
       write(unit, varFmt, advance = 'no') "ESOLVsp"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do isp = 1, rism%pot%nsp
          write(unit, valFmt, advance = 'no') solveneSP(isp)*KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if

    ! SOLVATION ENTROPY
    ! Test site value as species value is a sum.
    if (solvene(1) /= HUGE(1d0) .and. exchemPR(1) /= HUGE(1d0) .and. entropicDecomp == 1) then
       write(unit, descriptFmt, advance = 'no') "-Temperature*solvation_entropy_PR"//whtspc
       write(unit, varFmt, advance = 'no') "-TSsp_PR"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do isp = 1, rism%pot%nsp
          write(unit, valFmt, advance = 'no') -(solveneSP(isp)-exchemspPR(isp))&
               *KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if
    if (solvene(1) /= HUGE(1d0) .and. exchemSC(1) /= HUGE(1d0) .and. entropicDecomp == 1) then
       write(unit, descriptFmt, advance = 'no') "-Temperature*solvation_entropy_SC"//whtspc
       write(unit, varFmt, advance = 'no') "-TSsp_SC"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do isp = 1, rism%pot%nsp
          write(unit, valFmt, advance = 'no') -(solveneSP(isp)-exchemspSC(isp))&
               *KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if
    if (solvene(1) /= HUGE(1d0) .and. exchemSM(1) /= HUGE(1d0) .and. entropicDecomp == 1) then
       write(unit, descriptFmt, advance = 'no') "-Temperature*solvation_entropy_SM"//whtspc
       write(unit, varFmt, advance = 'no') "-TSsp_SM"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do isp = 1, rism%pot%nsp
          write(unit, valFmt, advance = 'no') -(solveneSP(isp)-exchemspSM(isp))&
               *KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if

    ! PARTIAL MOLAR VOLUME
    write(unit, descriptFmt, advance = 'no') "Partial_molar_volume"//whtspc
    write(unit, varFmt, advance = 'no') "PMV"//whtspc
    write(unit, unitFmt, advance = 'no') "[A^-3]"//whtspc
    do isp = 1, rism%pot%nsp
       write(unit, valFmt, advance = 'no') pmv(isp)
    end do
    write(unit, '(a)')

    ! Site.
    write(unit, '(a)')
    write(unit, descriptFmt) "#Site properties"//whtspc
    write(unit, descriptFmt, advance = 'no') "#Description"//whtspc
    write(unit, varFmt, advance = 'no') "Variable"//whtspc
    write(unit, unitFmt, advance = 'no') "Units"//whtspc
    do iv = 1, rism%pot%nv
       write(unit, titleFmt, advance = 'no') rism%pot%namev(iv)//whtspc
    end do
    write(unit, '(a)')


    ! EXCESS CHEMICAL POTENTIAL
    if (exchemPR(1) /= HUGE(1d0)) then
       write(unit, descriptFmt, advance = 'no') "Excess_chemical_potential_PR"//whtspc
       write(unit, varFmt, advance = 'no') "EXCHEMv_PR"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do iv = 1, rism%pot%nv
          write(unit, valFmt, advance = 'no') exchemPR(iv)*KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if
    if (exchemSC(1) /= HUGE(1d0)) then
       write(unit, descriptFmt, advance = 'no') "Excess_chemical_potential_SC"//whtspc
       write(unit, varFmt, advance = 'no') "EXCHEMv_SC"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do iv = 1, rism%pot%nv
          write(unit, valFmt, advance = 'no') exchemSC(iv)*KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if
    if (exchemSM(1) /= HUGE(1d0)) then
       write(unit, descriptFmt, advance = 'no') "Excess_chemical_potential_SM"//whtspc
       write(unit, varFmt, advance = 'no') "EXCHEMv_SM"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do iv = 1, rism%pot%nv
          write(unit, valFmt, advance = 'no') exchemSM(iv)*KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if

    ! SOLVATION ENERGY
    if (solvene(1) /= HUGE(1d0) .and. entropicDecomp == 1) then
       write(unit, descriptFmt, advance = 'no') "Solvation_energy"//whtspc
       write(unit, varFmt, advance = 'no') "ESOLVv"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do iv = 1, rism%pot%nv
          write(unit, valFmt, advance = 'no') solvene(iv)*KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if

    ! SOLVATION ENTROPY
    if (solvene(1) /= HUGE(1d0) .and. exchemPR(1) /= HUGE(1d0) .and. entropicDecomp == 1) then
       write(unit, descriptFmt, advance = 'no') "-Temperature*solvation_entropy_PR"//whtspc
       write(unit, varFmt, advance = 'no') "-TSv_PR"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do iv = 1, rism%pot%nv
          write(unit, valFmt, advance = 'no') -(solvene(iv)-exchemPR(iv))&
               *KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if
    if (solvene(1) /= HUGE(1d0) .and. exchemSC(1) /= HUGE(1d0) .and. entropicDecomp == 1) then
       write(unit, descriptFmt, advance = 'no') "-Temperature*solvation_entropy_SC"//whtspc
       write(unit, varFmt, advance = 'no') "-TSv_SC"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do iv = 1, rism%pot%nv
          write(unit, valFmt, advance = 'no') -(solvene(iv)-exchemSC(iv))&
               *KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if
    if (solvene(1) /= HUGE(1d0) .and. exchemSM(1) /= HUGE(1d0) .and. entropicDecomp == 1) then
       write(unit, descriptFmt, advance = 'no') "-Temperature*solvation_entropy_SM"//whtspc
       write(unit, varFmt, advance = 'no') "-TSv_SM"//whtspc
       write(unit, unitFmt, advance = 'no') "[kcal/mol]"//whtspc
       do iv = 1, rism%pot%nv
          write(unit, valFmt, advance = 'no') -(solvene(iv)-exchemSM(iv))&
               *KB*rism%pot%temperature
       end do
       write(unit, '(a)')
    end if

    if (iostat /= 0)&
         call rism_report_error("failed to open :"//trim(file))
    close(unit)
  end subroutine writeTherm


  !> Writes bulk thermodynamic quantities to file.
  subroutine writeTD()
    use rism_util, only : freeUnit
    use constants, only : boltzmann, KB
    implicit none
    _REAL_ :: pmv(rism%pot%nsp), exchem(rism%pot%nv), exchemsp(rism%pot%nsp), &
         exchemIon(rism%pot%nv), exchemIonsp(rism%pot%nsp), &
         compressibility, pressureFE, pressureVirial, FE
    _REAL_ :: solvene(rism%pot%nv), solvenesp(rism%pot%nsp)
    integer :: isp, iv, iat
    character(len = clen) :: file
    integer :: unit, iostat
    !field : Fortran likes to right-justify if the text string is too short
    !        for the specifier.  'field' is used as a work around
    character(len = 40) :: field
    pmv = rism1d_getPMV(rism)
    exchem = rism1d_getExChem(rism)
    exchemIon = rism1d_getExChemIon(rism)
    compressibility = rism1d_getCompressibility(rism)&
         *1d-20 / (boltzmann*rism%pot%temperature)
    pressureFE = rism1d_getPressureFE(rism)&
         *1.d24 * boltzmann*rism%pot%temperature
    pressureVirial = rism1d_getPressureVirial(rism)&
         *1.d24 * boltzmann*rism%pot%temperature
    FE = rism1d_getFreeEnergy(rism)&
         *KB*rism%pot%temperature
    call rism_timer_start(outputTimer)
    unit = freeUnit()
    file = trim(fileroot)//'.td'
    call rism_report_message('outputting thermodynamics to file: '//trim(file))
    open(unit, file = file, status = 'replace', iostat = iostat)
    if (iostat /= 0)&
         call rism_report_error("failed to open :"//trim(file))

    field = "!total compressibility    "
    write(unit, '(a40, a)') field, "Xi  [10-4/MPa]"
    field = "!pressure                 "
    write(unit, '(a40, a)') field, "P   [MPa]"
    field = "!excess free energy       "
    write(unit, '(a40, a)') field, "xA  [kcal/mol]"
    write(unit, '(a)') "&THERMO"
    write (unit, '(A, 1p, (g16.8, :, ", "))')  "XITOT = ", compressibility

    write (unit, '(A, 1p, (g16.8, :, ", "))')  "FE_PRESURE = ", pressureFE

    write (unit, '(A, 1p, (g16.8, :, ", "))')  "VIRIAL_PRESURE = ", pressureVirial

    write(unit, '(A, 1p, (g16.8, :, ", "))')   "XFE = ", FE
    write(unit, '(a)') "/"
    exchemsp = 0
    exchemIonsp = 0
    iv = 0
    do isp = 1, rism%pot%nsp
       do iat = 1, rism%pot%nat(isp)
          iv = iv + 1
          exchemsp(isp) = exchemsp(isp) + exchem(iv)
          exchemIonsp(isp) = exchemIonsp(isp) + exchemIon(iv)
       end do
    end do
    field = "!Molecular Species Contributions"
    write(unit, '(a40)') field
    field = "!partial molar volume"
    write(unit, '(a40, a)') field, "PMV [1/A^3]"
    field = "!excess chemical potential"
    write(unit, '(a40, a)') field, "xmu [kcal/mol]"
    write(unit, '(a)') "&PMV"
    write(unit, '(A, 1p, (g16.8, :, ", "))') "PVM = ", pmv
    write(unit, '(A, 1p, (g16.8, :, ", "))') "XMUSP = ", exchemsp*KB*rism%pot%temperature
    write(unit, '(A, 1p, (g16.8, :, ", "))') "XMUISP = ", exchemIonsp*KB*rism%pot%temperature
    write(unit, '(a)') "/"
    field = "!Site Contributions"
    write(unit, '(a40)') field
    field = "!excess chemical potential"
    write(unit, '(a40, a)') field, "xmu [kcal/mol]"
    write(unit, '(a)') "&XMU"
    write(unit, '(A, 1p, (g16.8, :, ", "))') "XMUV = ", exchem*KB*rism%pot%temperature
    write(unit, '(A, 1p, (g16.8, :, ", "))') "XMUIV = ", exchemIon*KB*rism%pot%temperature
    write(unit, '(a)') "/"

    solvene = rism1d_getSolvene(rism)
    solvenesp = 0
    iv = 0
    do isp = 1, rism%pot%nsp
       do iat = 1, rism%pot%nat(isp)
          iv = iv + 1
          solvenesp(isp) = solvenesp(isp) + solvene(iv)
       end do
    end do

    field = "!Analytical Temperature Derivative"
    write(unit, '(a40)') field
    field = "!Molecular Species Contributions"
    write(unit, '(a40)') field
    !
    field = "!excess chemical potential"
    write(unit, '(a40, a)') field, "xmu [kcal/mol]"
    field = "!solvation energy"
    write(unit, '(a40, a)') field, "esolv [kcal/mol]"
    field = "!solvation entropy"
    write(unit, '(a40, a)') field, "-ts [kcal/mol]"
    write(unit, '(a)') "&DECOMP"
    write(unit, '(A, 1p, (g16.8, :, ", "))') "XMUSP = ", exchemsp*KB*rism%pot%temperature
    write(unit, '(A, 1p, (g16.8, :, ", "))') "ESOLVSP = ", solvenesp*KB*rism%pot%temperature
    write(unit, '(A, 1p, (g16.8, :, ", "))') "-TSSP = ", (exchemsp-solvenesp)*KB*rism%pot%temperature
    write(unit, '(a)') "/"

    field = "!Site Contributions"
    write(unit, '(a40)') field
    field = "!excess chemical potential"
    write(unit, '(a40, a)') field, "xmu [kcal/mol]"
    field = "!solvation energy"
    write(unit, '(a40, a)') field, "esolv [kcal/mol]"
    field = "!solvation entropy"
    write(unit, '(a40, a)') field, "-ts [kcal/mol]"
    write(unit, '(a)') "&XMU"
    write(unit, '(A, 1p, (g16.8, :, ", "))') "XMUV = ", exchem*KB*rism%pot%temperature
    write(unit, '(a)') "&ESOLV"
    write(unit, '(A, 1p, (g16.8, :, ", "))') "ESOLVV = ", solvene*KB*rism%pot%temperature
    write(unit, '(a)') "&-TSV"
    write(unit, '(A, 1p, (g16.8, :, ", "))') "-TSV = ", (exchem-solvene)*KB*rism%pot%temperature
    write(unit, '(a)') "/"

    close(unit)
    call rism_timer_stop(outputTimer)
  end subroutine writeTD


  !> Writes the coordination number to file.  Use Gvv for coordination
  !! number and Hvv for excess coordination number.
  !! @param[in] vv Gvv or Hvv
  !! @param[in] suffix suffix for the end of the file
  !! @param[in] descript description of the column
  subroutine writeRunningCoordination(vv, suffix, descript)
    use rism_util
    use constants, only : pi
    implicit none
    _REAL_, intent(in) :: vv(:, :, :)
    character(len = *), intent(in) :: suffix, descript
    character(len = clen) :: file
    integer :: unit, iostat
    integer :: ivv, iv1, iv2, ir
    call rism_timer_start(outputTimer)
    unit = freeUnit()
    file = trim(fileroot)//'.'//suffix
    call rism_report_message('outputting running coordination number to file: '//trim(file))
    open(unit, file = file, status = 'replace', iostat = iostat)
    if (iostat /= 0)&
         call rism_report_error("failed to open :"//trim(file))
    call writeVVHeader(unit, &
         "#RISM1D ATOM-ATOM INTERACTIONS: "//trim(descript)//" A:B (# of B around A) VS. SEPARATION [A]", 17, .true.)

    do ir = 1, nrout
       write (unit, '(1p, 5050(1x, e16.8e3))')  (ir-1)*rism%pot%dr, rmExPrec(vv(ir, :, :))
    end do

    close(unit)
    call rism_timer_stop(outputTimer)
  end subroutine writeRunningCoordination


  !> Writes the total excess coordination number and site-site
  !! structure factor to file.
  !! @param[in] excess Excess quantity of sepcies in dimension 1 about
  !!  species in dimension 2 (nv, nv).
  !! @param[in] suffix Suffix for the file.
  !! @param[in] descript Description of the column.
  !! @param[in] total (optional) Write out the sum of each row.
  subroutine writeTotalExcess(excess, suffix, descript, total)
    use rism_util, only : freeUnit
    implicit none
    _REAL_, intent(in) :: excess(:, :)
    character(len = *), intent(in) :: suffix, descript
    logical, optional, intent(in) :: total
    logical :: optTotal
    character(len = clen) :: file
    integer :: unit, iostat
    integer :: ivv, iv1, iv2, ir
    character(len = CLEN) :: colTitle_fmt, colData_fmt, advance
    call rism_timer_start(outputTimer)
    unit = freeUnit()

    advance = 'yes'
    optTotal = .false.
    if (present(total)) optTotal = total
    if (optTotal) advance = 'no'

    file = trim(fileroot)//'.'//suffix
    call rism_report_message('outputting Nvv(R = infty) to file: '//trim(file))
    open(unit, file = file, status = 'replace', iostat = iostat)
    if (iostat /= 0)&
         call rism_report_error("failed to open :"//trim(file))

    write(colTitle_fmt, '(a, i4, a)') '(4x, ', rism%pot%nv, '(7x, a4, 6x))'
    write(colData_fmt, '(a, i4, a)') '(a4, 1p, ', rism%pot%nv, '(1x, e16.8e3))'

    !we will assume that atom names are four characters long form simplicity
    write(unit, '(a)') "#"//descript//" of column site about row site"
    write(unit, colTitle_fmt, advance = advance) rism%pot%namev
    if (optTotal) then
       write(unit, '(a)') "   Total charge  "
    end if
    do iv1 = 1, rism%pot%nv
       write (unit, colData_fmt, advance = advance)  rism%pot%namev(iv1), excess(iv1, :)
       if (optTotal) then
          write(unit, '(1p, 2(1x, e16.8e3))') sum(excess(iv1, :))
       end if
    end do
    close (unit)
    call rism_timer_stop(outputTimer)
  end subroutine writeTotalExcess


  !> Writes the total excess coordination number and site-site structure factor to
  !! file.
  subroutine writeSVV()
    implicit none
    _REAL_, pointer :: wrkvv(:, :) => NULL()
    !.......................... outputting Svv(k) ..........................
    wrkvv => rism1d_getStructFactor(rism)
    call writeVV(wrkvv, 'svv', 'STRUCTURE FACTOR', k_space = .true.)
    if (safemem_dealloc(wrkvv) /= 0)&
         call rism_report_error("memory deallocation failed in writeSVVN00")
  end subroutine writeSVV


  !> Writes out a header for VV files to describe the contents.  Columns are
  !! labelled by species and site numbers. All columns should be the same width to
  !! acheive the correct centering of columns.
  !! @param[in] unit : Unit to write to.
  !! @param[in] head : One line description of the file.
  !! @param[in] width : Column width in characters.
  !! @param[in] full : (optional) If true do the full NV*NV header. If false do the
  !!  site-site header only.  Default: .false.
  !! @param[in] k_space (optional) If true, data is in reciprocal space and the frequency
  !!  will be printed out.  If false (default), data is in real space
  !!  and the separation will be printed out.
  subroutine writeVVHeader(unit, head, width, full, k_space)
    implicit none
    integer, intent(in) :: unit, width
    character(len = *), intent(in) :: head
    logical, optional, intent(in) :: full, k_space
    logical :: all, kspace
    integer, pointer :: interactions(:, :) => NULL()
    integer :: ivv, iv1, iv2, isp, isp2, iat, iat2
    !sepbuffer : number of buffer spaces before and after 'separation', the first column
    !pairbuffer : number of buffer spaces before and after each atom pair column
    !pairlength : number of spaces required for the species-site
    !numbers.  e.g. S1A3:S11A100 would be (1, 1, 2, 3) for the
    !          character lengths of "1", "3", "11", "100"
    integer :: sepbuffer(2), pairbuffer(2), pairlength(4)
    !form : we will write the format to this string
    !form1col : for the first column in the head.  This gets reused
    !xlabel : SEPARATION or FREQUENCY depending on the values of kspace
    character(len = 256) :: form, form1col, xlabel
    all = .false.
    if (present(full)) all = full
    kspace = .false.
    if (present(k_space)) kspace = k_space

    xlabel = "SEPARATION"
    if (kspace) xlabel = "FREQUENCY"

    write (unit, '(a)') head
    write(unit, '(a)') "#S=SPECIES, A=ATOM"

    !get the correct order for the labels.  We assume that the internal RISM1D
    !ordering is used
    ivv = 0
    if (all) then
       interactions => safemem_realloc(interactions, 4, rism%pot%nv**2)
       do isp = 1, rism%pot%nsp
          do iat = 1, rism%pot%nat(isp)
             do isp2 = 1, rism%pot%nsp
                do iat2 = 1, rism%pot%nat(isp2)
                   ivv = ivv + 1
                   interactions(:, ivv) = (/isp, iat, isp2, iat2/)
                end do
             end do
          end do
       end do
    else
       interactions => safemem_realloc(interactions, 4, rism%pot%nvv)
       do isp = 1, rism%pot%nsp
          do iat = 1, rism%pot%nat(isp)
             do isp2 = 1, isp
                do iat2 = 1, rism%pot%nat(isp2)
                   if (isp2 .eq. isp .and. iat2 .gt. iat) then
                      exit
                   else
                      ivv = ivv + 1
                      interactions(:, ivv) = (/isp, iat, isp2, iat2/)
                   end if
                end do
             end do
          end do
       end do
    end if

    !first column of the species-site number row
    !get the number of blanks needed to fill the column width
    sepbuffer(1) = width - len("#"//trim(xlabel))
    !use half minus one for the right side of the column
    sepbuffer(2) = sepbuffer(1) - sepbuffer(1) / 2 - 1
    !use half plus one for the left side of the column
    sepbuffer(1) = sepbuffer(1) / 2 + 1

    write(form1col, '(a, i2, a, i2, a)') '("#", ', sepbuffer(1), &
         'x, "'//trim(xlabel)//'", ', sepbuffer(2), "x)"
    write(unit, form1col, advance = 'no')
    do ivv = 1, ubound(interactions, 2)
       !do the same for the atom pair columns.  Assume only one character
       !for the species and site numbers
       pairlength = floor(log10(dble(interactions(:, ivv)))) + 1
       pairbuffer(1) = width - len("SA:SA") - sum(pairlength)
       pairbuffer(2) = pairbuffer(1) - pairbuffer(1) / 2 -1
       pairbuffer(1) = pairbuffer(1) / 2 +1

       !write the format to string
       write(form, '(a, i2, a, i2, a, i2, a, i2, a, i2, a, i2, a)') "(1000(:, ", &
            pairbuffer(1), &
            'x, "S", i', pairlength(1), ', "A", i', pairlength(2), &
            ', ":S", i', pairlength(3), ', "A", i', pairlength(4), &
            ', ', pairbuffer(2), 'x))'
       !write column headers
       write(unit, form, advance = 'no') interactions(:, ivv)
    end do
    write(unit, *)

    !do it again with atom names
    write(unit, form1col, advance = 'no')
    ivv = 0
    if (all) then
       do iv1 = 1, rism%pot%nv
          do iv2 = 1, rism%pot%nv
             pairbuffer(1) = width - len(trim(rism%pot%namev(iv1))//":"//trim(rism%pot%namev(iv2)))
             pairbuffer(2) = pairbuffer(1) - pairbuffer(1) / 2 - 1
             pairbuffer(1) = pairbuffer(1) / 2 + 1
             write(form, '(a, i2, a, i2, a)') '(', pairbuffer(1), 'x, a, ', pairbuffer(2), 'x)'
             write(unit, form, advance = 'no')trim(rism%pot%namev(iv1))//":"//trim(rism%pot%namev(iv2))
          end do
       end do
    else
       do iv1 = 1, rism%pot%nv
          do iv2 = 1, iv1
             pairbuffer(1) = width - len(trim(rism%pot%namev(iv1))//":"//trim(rism%pot%namev(iv2)))
             pairbuffer(2) = pairbuffer(1) - pairbuffer(1) / 2 - 1
             pairbuffer(1) = pairbuffer(1) / 2 + 1
             write(form, '(a, i2, a, i2, a)') '(', pairbuffer(1), 'x, a, ', pairbuffer(2), 'x)'
             write(unit, form, advance = 'no') trim(rism%pot%namev(iv1))//":"//trim(rism%pot%namev(iv2))
          end do
       end do
    end if
    write(unit, *)
    if (safemem_dealloc(interactions) /= 0)&
         call rism_report_error("writeVVHeader: could not deallocate interactions")
  end subroutine writeVVHeader

end module rism1d_m


program rism1d_program
  use rism1d_m
  implicit none
  call initialize()
  call solve()
  call output()
  call finalize()
end program rism1d_program


subroutine timer_start(label)
  integer label
end subroutine timer_start


subroutine timer_stop(label)
  integer label
end subroutine timer_stop
