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

!> Site-site correlations and thermodynamics of a molecular mixture.
!! o Site-site PDFs, TCFs, DCFs, and susceptibility and their
!!   temperature derivatives.
!! o Compressibility.
!! o Site-site RISM equation and its temperature derivative.
!! o HNC, KH, PSE-n, MV0, PY closures
!! o HNC, KH and PSE-n closure temperature derivatives
!! o temperature derivatives expressed as T*d/dT. This is used throughout.
!! o Solving for the Direct Correlation Functions.
!! o Using the linearly spaced 1D-FFT.
!! o Calculating the residual as r*(G_clos(r)-1-H_rism(r)).
!! o Compact triangle ordering of site-site vectors.
!! o Modified direct inversion in the iterative subspace (MDIIS).
!! o Coulomb + 12-6 Lennard-Jones site-site potentials.
!! o units:  energy       [kT]
!!           distances    [A]           (Angstroms)
!!           site charges [sqrt(kT A)]
!!           temperature  [K]
!!           density      [#/A^3]
!! o To convert [e] to [sqrt(kT A)] * sqrt(COULOMB_CONST_E / (KB * temperature)).
module rism1d_c
  use safemem
  use rism1d_potential_c
  use rism1d_closure_c
  use rism_report_c
  use rism_timer_c
  use mdiis_c
  implicit none
  integer,private ,parameter :: charlen = 8, maxep0=4
  type rism1d
     !closureID : string identifier for the closure
     character(len=charlen) :: closureID=''
     !savefile  : save file to allow restarts
     character(len=256) :: savefile='rism.sav'

     !pot : potential object
     type(rism1d_potential) :: pot

     !closure : Closure object.  Current options are 'KH', 'HNC', 'PY', VM0'
     type(rism1d_closure) :: closure

     !extra_precision :: controls the precision in key parts of the algorithm.
     !                   0 - no extra precision
     !                   1 - XBLAS DGEMM for Ak and Bk in single3DRISMsolution()
     integer :: extra_precision=1
     
     !Mdiis_nvec : number of MDIIS vectors
     integer :: Mdiis_nvec=0

     !mdiis_del       : MDIIS step size
     !mdiis_restart :: restart threshold factor. Ratio of the current residual to the 
     !                 minimum residual in the basis that causes a restart
     _REAL_ :: mdiis_del, mdiis_restart
     
     integer :: mdiis_method =2
    
     !gvv   : pair correlation function
     !hvv   : total correlation function
     !cvvWRK  : last Mdiis_nvec steps of cvv
     !cvvresWRK: last Mdiis_nvec steps of cvvres
     _REAL_, pointer :: &
          gvv(:,:)=>NULL(), hvv(:,:)=>NULL(), &
          cvvWRK(:,:,:)=>NULL(), cvvresWRK(:,:,:)=>NULL()

     !gvv   : temperature derivative of pair correlation function
     !hvv   : temperature derivative of total correlation function
     !cvvWRK_dT  : last mdiis_nvec steps of cvv_dT
     !cvvresWRK_dT: last mdiis_nvec steps of cvv_dTres
     _REAL_, pointer :: &
          gvv_dT(:,:)=>NULL(), hvv_dT(:,:)=>NULL(), &
          cvvWRK_dT(:,:,:)=>NULL(), cvvresWRK_dT(:,:,:)=>NULL()

     !convenience pointers
     !cvv   : direct correlation function.  Points to cvvWRK(:,:,1)
     !cvvres: residual direct correlation function (residual between iterations).
     !        Note that cvvres is also used as workspace in single3DRISMsolution* and doesn't 
     !        necessarily contain residuals after a call to solve3DRISM*. Points to 
     !        cuvresWRK(:,:,1)
     _REAL_,pointer :: cvv(:,:)=>NULL(), cvvres(:,:)=>NULL()

     !convenience pointers
     !cvv_dT   : temperature derivative of direct correlation function.  Points to cvvWRK_dT(:,:,1)
     !cvv_dTres: residual temperature derivative direct correlation function (residual between iterations)
     !        Points to cvvresWRK_dT(:,:,1)
     _REAL_,pointer :: cvv_dT(:,:)=>NULL(), cvv_dTres(:,:)=>NULL()

     !Computing the k-space contribution to the free energy and the pressure 
     !along the free energy route can be done at the same for the same cost
     !as a single calculation.  We store the values here.  These are reset each
     !solution to indicate that they need to be recalculated

     !TIMERS.  Subtimers only account for computation.  We ignore setup etc.
     !timer :: timer for this class.  Activated for all public routines
     !potTimer :: specifically times potential calculation
     !fftTimer :: specifically times FFT calculation
     !solveTimer :: specifically times rism1d_solve calculation
     !solve3DRISMTimer :: specifically times solve3DRISM calculation
     !single3DRISMsolutionTimer :: specifically times single3DRISMsolution calculation
     !fft_dTTimer :: specifically times FFT calculation for temperature
     !              derivative calculation
     !solve_dTTimer :: specifically times rism1d_solve calculation for
     !                temperature derivative calculation
     !solve3DRISM_dTTimer :: specifically times solve3DRISM calculation for
     !                 temperature derivative calculation
     !single3DRISMsolution_dTTimer :: specifically times single3DRISMsolution calculation for
     !                 temperature derivative calculation
     !thermoTimer :: specifically times thermodynamics calculations
     !selftestTimer :: specifically times self-test calculations
     type(rism_timer) :: timer, potTimer, &
          fftTimer, solventTimer, solve3DRISMTimer, single3DRISMsolutionTimer, &
          fft_dTTimer, solvent_dTTimer, solve3DRISM_dTTimer, single3DRISMsolution_dTTimer, &
          thermoTimer, selftestTimer
     
  end type rism1d

  interface rism1d_new
     module procedure new_
  end interface rism1d_new

  public rism1d_new, rism1d_destroy, rism1d_solve, rism1d_addSpecies,&
       rism1d_getInvDebyeLen, rism1d_getCompressibility , rism1d_getDelHvLimit,& 
       rism1d_getExNumber, rism1d_getStructFactor, rism1d_getRunNumber, &
       rism1d_getRunExNumber, rism1d_getPressureVirial, rism1d_getPressureFE, rism1d_getFreeEnergy, &
       rism1d_getPMV, rism1d_getExChem, rism1d_getSolvene, rism1d_readNBFix, rism1d_readNBFixtsf
  private new_, solve3DRISM, single3DRISMsolution, sanity_check
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initializes the in rism1d object
!!!IN:
!!!   this : rism1d object
!!!   theory : 'xrism' or 'drism'
!!!   closure : 'kh','hnc', 'py', 'vm0', 'psen', 'v*', 'nub'
!!!   coeff  : (optional) coefficients for the selected closure
!!!   vstarB : (optional) 'B' parameter for Verlet (V*) closure
!!!   temperature : system temperature in [K]
!!!   dielconst : dielconst for DRISM calculations
!!!   smear  : electrostatic smear parameter
!!!   adbcor : ???  
!!!   nr     : number of grid points
!!!   dr     : grid spacing
!!!   Mdiis_nvec    : number of MDIIS vectors
!!!   mdiis_del    : MDIIS step size
!!!   mdiis_restart :: restart threshold factor. Ratio of the current residual to the 
!!!                    minimum residual in the basis that causes a restart
!!!   savefile : name of the restart file
!!!   extra_precision : use extra_precision in certain key parts of the code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine new_(this, theory, closure, coeff, temperature, dielconst, smear, adbcor, &
       nr, dr, Mdiis_nvec,mdiis_del,mdiis_restart,savefile,extra_precision)
    use constants, only : pi
    use rism_util, only : caseup
    implicit none
    type(rism1d), intent(inout) :: this
    character(len=*), intent(in) :: theory, closure
    _REAL_, optional, intent(in) :: coeff(:)
    character(len=*), optional, intent(in) :: savefile
    integer, intent(in) :: nr, Mdiis_nvec, extra_precision
    _REAL_, intent(in) :: temperature, dielconst, smear, adbcor, dr, mdiis_del, mdiis_restart
    integer :: factored_nr
    _REAL_  :: l_vstarA, l_vstarB
    l_vstarA=0
    l_vstarB=0
    call rism_timer_new(this%timer, "1D-RISM Total")
    call rism_timer_start(this%timer)
    call rism_timer_new(this%thermoTimer, "Thermodynamics",this%timer)
    call rism_timer_new(this%solventTimer, "Solve 1D-RISM",this%timer)
    call rism_timer_new(this%potTimer, "Potential",this%solventTimer)
    call rism_timer_new(this%solve3DRISMTimer, "RXRISM",this%solventTimer)
    call rism_timer_new(this%single3DRISMsolutionTimer, "R1RISM",this%solve3DRISMTimer)
    call rism_timer_new(this%fftTimer, "FFT",this%single3DRISMsolutionTimer)
    call rism_timer_new(this%solvent_dTTimer, "Solve 1D-RISM dT",this%timer)
    call rism_timer_new(this%solve3DRISM_dTTimer, "RXRISM dT",this%solvent_dTTimer)
    call rism_timer_new(this%single3DRISMsolution_dTTimer, "R1RISM dT",this%solve3DRISM_dTTimer)
    call rism_timer_new(this%fft_dTTimer, "FFT dT",this%single3DRISMsolution_dTTimer)
    call rism_timer_new(this%selftestTimer, "Self-test",this%timer)

    call caseup(closure)
    this%closureID = closure
    this%Mdiis_nvec = Mdiis_nvec
    this%mdiis_del = mdiis_del
    this%mdiis_restart = mdiis_restart
    this%extra_precision = extra_precision

    call sanity_check(this)

    !THIS%POT MUST BE INITIALIZED FIRST
    call rism1d_potential_new(this%pot,theory,temperature,dielconst,smear,adbcor,nr,dr,&
         this%extra_precision)

    if(present(coeff)) then
       call rism1d_closure_new(this%closure,this%closureID,this%pot,&
            coeff=coeff)
    else
       call rism1d_closure_new(this%closure,this%closureID,this%pot)
    end if

    if(present(savefile))then
       this%savefile = trim(savefile)
    end if

    call rism_timer_stop(this%timer)
  end subroutine new_

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Set parent for this timer
!!!IN:
!!!   this : rism1d object
!!!   parent : parent timer object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_setTimerParent(this, parent)
    implicit none
    type(rism1d), intent(inout) :: this
    type(rism_timer), intent(inout) :: parent
    call rism_timer_start(this%timer)
    call rism_timer_setParent(this%timer,parent)
    call rism_timer_stop(this%timer)
  end subroutine rism1d_setTimerParent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Add solvent species.  Allocates necessary memory to accomodate the additional
!!!species
!!!IN:
!!!   this : rism1d object
!!!   mdl  : solvMDL object
!!!   density : number density (A^{-3})
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_addSpecies(this,mdl,density)
    use solvMDL_c
    implicit none
    type(rism1d), intent(inout) :: this
    type(solvMDL), intent(in) :: mdl
    _REAL_, intent(in) :: density
    integer*8 :: memstats(10)
    call rism_timer_start(this%timer)

    call rism1d_potential_addSpecies(this%pot,mdl,density)
    !data arrays in the original 1D-RISM code had an offset of 0 in the first index (like C).
    !The offset is 1 in the current code (like most Fortran).
    this%gvv => safemem_realloc(this%gvv,this%pot%nr,this%pot%nvv)
    this%hvv => safemem_realloc(this%hvv,this%pot%nr,this%pot%nvv)
    this%cvvWRK => safemem_realloc(this%cvvWRK,this%pot%nr,this%pot%nvv,this%Mdiis_nvec)
    this%cvvresWRK => safemem_realloc(this%cvvresWRK,this%pot%nr,this%pot%nvv,this%Mdiis_nvec)
    this%cvv => this%cvvWRK(:,:,1)
    this%cvvres => this%cvvresWRK(:,:,1)


    this%gvv_dT => safemem_realloc(this%gvv_dT,this%pot%nr,this%pot%nvv)
    this%hvv_dT => safemem_realloc(this%hvv_dT,this%pot%nr,this%pot%nvv)
    this%cvvWRK_dT => safemem_realloc(this%cvvWRK_dT,this%pot%nr,this%pot%nvv,this%mdiis_nvec)
    this%cvvresWRK_dT => safemem_realloc(this%cvvresWRK_dT,this%pot%nr,this%pot%nvv,this%mdiis_nvec)

    this%cvv_dT => this%cvvWRK_dT(:,:,1)
    this%cvv_dTres => this%cvvresWRK_dT(:,:,1)
    call rism_timer_stop(this%timer)
  end subroutine rism1d_addSpecies

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Solves the 1D-RISM system of equations.
!!!IN:
!!!   this : rism1d object
!!!   maxstep : max number of iterations
!!!   ksave  : frequency of restart file saves
!!!   progress  : frequency of progress updates
!!!   tolerance    : tolerance
!!!OUT:
!!!    .true. if converged.  .false. otherwise.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_solve(this, ksave, progress, maxstep, tolerance) result(converged)
    implicit none
    type(rism1d), intent(inout) :: this
    integer, intent(in) :: maxstep, ksave, progress
    _REAL_, intent(in) :: tolerance
    logical :: converged
    _REAL_ :: charge
    integer :: err
    call rism_timer_start(this%solventTimer)


    !................ for MV0 closure, solving RISM at Q=0 .................
    if (rism1d_closure_isCharged(this%closure)) then
       charge = 1d0
    else
       charge = 0d0
    end if
    call rism_timer_start(this%potTimer)
    call rism1d_potential_calc(this%pot,charge)
    if (rism1d_closure_type(this%closure) .eq. "MV0") &
         call rism1d_mv0_calcUvv(this%closure%mv0, this%pot%nr, this%pot%dr, &
         this%pot%qv, this%pot%ljEpsilonVV, this%pot%ljRminVV, this%pot%smear)
    call rism_timer_stop(this%potTimer)

    !........................ solving RISM by MDIIS ........................
    converged =  solve3DRISM(this,ksave,progress,maxstep,tolerance)

    call rism_timer_stop(this%solventTimer)
  end function rism1d_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Solves the 1D-RISM temperature derivative system of equations.
!!!IN:
!!!   this : rism1d object
!!!   maxstep : max number of iterations
!!!   ksave  : frequency of restart file saves
!!!   progress  : frequency of progress updates
!!!   tolerance    : tolerance
!!!OUT:
!!!    .true. if converged.  .false. otherwise.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_dt_solve(this, ksave, progress, maxstep, tolerance) result(converged)
    implicit none
    type(rism1d), intent(inout) :: this
    integer, intent(in) :: maxstep, ksave, progress
    _REAL_, intent(in) :: tolerance
    logical :: converged
    _REAL_ :: charge
    integer :: err
    call rism_timer_start(this%solvent_dTTimer)

    !........................ solving RISM DT by MDIIS ........................
    converged =  solve3DRISM_dT(this,ksave,progress,maxstep,tolerance)
    call rism_timer_stop(this%solvent_dTTimer)

  end function rism1d_dt_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Check internal consistency by calculating a number of properties
!!!and writing the results to file.  This is a self-test for the
!!!particular solution (w/ or w/o temperature derivatives) already
!!!solved.  Explicit pass/fails are not generated as this depends on
!!!the numeric precision request for the solution.  Results are
!!!written to the given file name. By default, only tests that are
!!!expected to pass are performed.  E.g., the simple virial expression
!!!for pressure is not correct for multi-site species so pressure is
!!!not tested for such cases.  To run all tests, regardless of
!!!applicability, use 'all=.true.'.
!!!IN:
!!!   this : rism1d object
!!!   filename : name of file to write results to
!!!   o_all: (optional) if .true., run all test calculations regardless
!!!          of applicability
!!!SIDEFFECTS:
!!!   Creates/overwrites the given file name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_selftest(this,filename,o_all)
    use rism_util, only : freeunit
    implicit none
    type(rism1d), intent(inout) :: this
    character(len=*) , intent(in) :: filename
    logical, optional, intent(in) :: o_all
    logical :: all
    character(len=4) :: char_nsite
    integer :: unit, iostat
    _REAL_ :: qEx(this%pot%nv,this%pot%nv)
    _REAL_ :: exchemSC(this%pot%nv),exchemPR(this%pot%nv),exchemSM(this%pot%nv)
    _REAL_ :: FE, pressureFE, pressureVirial
    integer :: iv

    all = .false.
    if(present(o_all)) all = o_all

    ! total excess charge from Kirkwood-Buff theory
    qEx = rism1d_getExNumber(this)

    !only calculate exchem for XRISM or if there are no site charges
    !in the system (DRISM == XRISM in this case).
    !If the closure does not support these calculations, HUGE() will be returned
    if(all .or. this%pot%theory .eq. "XRISM" .or. sum(abs(this%pot%qv)) == 0d0) then
       exchemPR = rism1d_getExChem(this,'PR')
       exchemSC = rism1d_getExChem(this,'SC')
    else 
       exchemPR = HUGE(1d0)
       exchemSC = HUGE(1d0)
    end if
    exchemSM = rism1d_getExChem(this,'SM')

    fe = rism1d_closure_getFreeEnergy(this%closure,this%gvv,this%cvv)

    pressureFE=rism1d_closure_getPressureFE(this%closure,this%gvv,this%cvv)
    !virial only works for point particles. Only do the calculation if
    !each species has one type and each site has a multiplicity of 1
    if(all .or. &
       (sum(this%pot%mtv) == this%pot%nv .and. &
       sum(this%pot%nat) == this%pot%nsp))then
       pressureVirial=rism1d_closure_getPressureVirial(this%closure,this%gvv)
    else
       pressureVirial = HUGE(1d0)
    end if

    call rism_timer_start(this%selftestTimer)

    write(char_nsite,'(i4)') this%pot%nv

    unit = freeUnit()
    open(unit,file=filename,status='replace',iostat=iostat)
    if(iostat/=0)&
         call rism_report_error("failed to open :"//trim(filename))
    
    !1) charge neutrality - input vs. total excess
    write(unit,'(a)') 'NET CHARGE NEUTRALITY [sqrt(kT A)]'
    write(unit,'(1p,a,e24.16)') '  ZERO Input from MDL        ',sum(this%pot%qv*this%pot%densityv)
    !get excess number of all other sites about each site
    !convert to charge and get the sum for each site (excess charge around each site)
    do iv =1, this%pot%nv
       qEx(iv,1) = sum(qEx(iv,:)*this%pot%qv*this%pot%mtv)
    end do
    !print out the sum of the excess charges for each site.
    write(unit,'(1p,a,e24.16)') '  ZERO Sum of excess charges ',sum(qEx(:,1))

    write(unit,*)

    !2) exchem - Singer-Chandler vs. Pettitt-Rossky vs. Schmeer-Maurer
    if(exchemSM(1) /= huge(1d0))then
       write(unit,'(a)') 'EXCESS CHEMICAL POTENTIAL [kT]'
       write(unit,'(a,'//trim(char_nsite)//'(a25))')       '  Equation             ',&
            this%pot%namev
       if(exchemPR(1) /= huge(1d0))&
            write(unit,'(1p,a,'//trim(char_nsite)//'(e24.16))') '  Pettitt-Rossky (PR)  ',exchemPR
       if(exchemSC(1) /= huge(1d0))&
            write(unit,'(1p,a,'//trim(char_nsite)//'(e24.16))') '  Singer-Chandler (SC) ',exchemSC
       write(unit,'(1p,a,'//trim(char_nsite)//'(e24.16))') '  Schmeer-Maurer (SM)  ',exchemSM
       if(exchemSC(1) /= huge(1d0) .or.exchemPR(1) /= huge(1d0) )&
            write(unit,'(a)') 'Relative difference [kT]'
       if(exchemPR(1) /= huge(1d0))then
          write(unit,'(1p,a,'//trim(char_nsite)//'(e24.16))') '  ZERO ExChem SC to PR ',&
               (exchemSC-exchemPR)/exchemPR
          write(unit,'(1p,a,'//trim(char_nsite)//'(e24.16))') '  ZERO ExChem SM to PR ',&
               (exchemSM-exchemPR)/exchemPR
       end if
       if(exchemSC(1) /= huge(1d0))&
            write(unit,'(1p,a,'//trim(char_nsite)//'(e24.16))') '  ZERO ExChem SM to SC ',&
            (exchemSM-exchemSC)/exchemSC
       write(unit,*)
    end if

    !3) Free energy vs. sum of exchem
    if(exchemSM(1) /= huge(1d0))then
       write(unit,'(a)') 'TOTAL EXCESS FREE ENERGY PER UNIT VOLUME [kT/A^3]'
       write(unit,'(1p,a,e24.16)') '  Free energy                          ',FE
       if(exchemPR(1) /= huge(1d0))&
            write(unit,'(1p,a,e24.16)') '  Sum of ExChem (PR) - Excess pressure ',&
            sum(exchemPR * this%pot%densityv/this%pot%mtv)-pressureFE+sum(this%pot%densitysp)
       if(exchemSC(1) /= huge(1d0))&
            write(unit,'(1p,a,e24.16)') '  Sum of ExChem (SC) - Excess pressure ',&
            sum(exchemSC * this%pot%densityv/this%pot%mtv)-pressureFE+sum(this%pot%densitysp)
       write(unit,'(1p,a,e24.16)') '  Sum of ExChem (SM) - Excess pressure ',&
            sum(exchemSM * this%pot%densityv/this%pot%mtv)-pressureFE+sum(this%pot%densitysp)
       write(unit,'(a)') 'Relative difference [kT/A^3]'
       if(exchemPR(1) /= huge(1d0))&
            write(unit,'(1p,a,e24.16)') '  ZERO ExChem (PR) - ExP to Free Energy',&
            (FE-(sum(exchemPR * this%pot%densityv/this%pot%mtv)-pressureFE+sum(this%pot%densitysp)))/FE
       if(exchemSC(1) /= huge(1d0))&
            write(unit,'(1p,a,e24.16)') '  ZERO ExChem (SC) - ExP to Free Energy',&
            (FE-(sum(exchemSC * this%pot%densityv/this%pot%mtv)-pressureFE+sum(this%pot%densitysp)))/FE
       write(unit,'(1p,a,e24.16)') '  ZERO ExChem (SM) - ExP to Free Energy',&
            (FE-(sum(exchemSM * this%pot%densityv/this%pot%mtv)-pressureFE+sum(this%pot%densitysp)))/FE
       write(unit,*)
    end if

    !4) pressure - virial vs. free energy
    !only applies to XRISM right now
    if(pressureFE /= HUGE(1D0) .and. pressureVirial /= HUGE(1d0))then
       write(unit,'(a)') 'Pressure [kT/A^3]'
       write(unit,'(1p,a,e24.16)') '  Pressure (free energy) ',pressureFE
       write(unit,'(1p,a,e24.16)') '  Pressure (virial)      ',pressureVirial
       write(unit,'(a)') 'Relative difference [kT/A^3]'
       write(unit,'(1p,a,e24.16)') '  ZERO Pressure          ',&
            (pressureFE-pressureVirial)/pressureFE
    end if

    !other ideas
    !a) dielectric constant - calculated vs. predicted(XRISM)/requested(DRISM)

    close(unit)
    call rism_timer_stop(this%selftestTimer)
  end subroutine rism1d_selftest

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the inverse Debye length (1/A) of the solvent
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getInvDebyeLen(this) result(kappa)
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: kappa
    call rism_timer_start(this%thermoTimer)
    kappa=rism1d_potential_getInvDebyeLen(this%pot)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getInvDebyeLen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the compressibility of the solvent [A^3/kT]
!!!IN:
!!!   this : rism1d object
!!!OUT:
!!!   isothermal compressibility temperature derivative
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getCompressibility(this) result(xikt)
    use constants, only : pi
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: xikt
    call rism_timer_start(this%thermoTimer)
    xikt = rism1d_closure_getCompressibility(this%closure,this%cvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getCompressibility

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the compressibility temperature derivative (T*d/dT) of the
!!!solvent [A^3]
!!!IN:
!!!   this : rism1d object
!!!OUT:
!!!   isothermal compressibility temperature derivative
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getCompressibility_dT(this) result(xikt_dT)
    use constants, only : pi
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: xikt_dT
    call rism_timer_start(this%thermoTimer)
    xikt_dT = rism1d_closure_getCompressibility_dT(this%closure,this%cvv,this%cvv_dT)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getCompressibility_dT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Returns the the extrapolated value for
!!! DelHv=-Lim_k->0 ( Sum_v1 Qv1*Xv1v2(k)4pi/k^2 )
!!! This is the background charge contribution for periodic and
!!! non-periodic RISM without Q_solute, Q (Eq. 20 from Kovalenko &
!!! Hirata, JCP 2000).
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_get3DRISMbackground(this) result(background)
    use constants, only : pi
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: background(this%pot%nv)
    call rism_timer_start(this%thermoTimer)
    background = rism1d_closure_get3DRISMBackground(this%closure,this%hvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_get3DRISMbackground

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Returns the the extrapolated value for
!!! DelHv=-Lim_k->0 ( Sum_v1 Qv1*Xv1v2(k)4pi/k^2 - hlkv0 )
!!! This is the entire h(k=0) for non-periodic RISM without Q_solute,
!!! combining \Delta Q (Eq. 20 from Kovalenko & Hirata, JCP 2000) and
!!! the long-range asymptotics value at k=0 (Eq. 11 from Kaminski et
!!! al., JPCA, 2010 and Eq. 10a from Genheden et al. JPCB, 2010).
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getDelHvLimit(this) result(delhv0)
    use constants, only : pi
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: delhv0(this%pot%nv)
    call rism_timer_start(this%thermoTimer)
    delhv0 = rism1d_closure_getDelHvLimit(this%closure,this%hvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getDelHvLimit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the the extrapolated value for the temperature derivative
!!!(T*d/dT) of DelHv0.
!!!DelHv_dT=-Lim_k->0 ( Sum_v1 Qv1*Xv1v2_dT(k)4pi/k^2 - hlkv0 )
!!!This is used by 3D-RISM long range asymptotics
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getDelHvLimit_DT(this) result(delhv0_dT)
    use constants, only : pi
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: delhv0_dT(this%pot%nv)
    delhv0_dT = rism1d_closure_getDelHvLimit_DT(this%closure,this%hvv,this%hvv_dT)
  end function rism1d_getDelHvLimit_DT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a pointer to the site-site susceptibility.
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getSusceptibility(this) result(xvv)
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_, pointer :: xvv(:,:,:)
    call rism_timer_start(this%thermoTimer)
    call rism1d_closure_calcXvv(this%closure, this%hvv)
    xvv => this%closure%xvv
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getSusceptibility

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a pointer to the temperature derivative of site-site susceptibility.
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getSusceptibility_DT(this) result(xvv_dT)
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_, pointer :: xvv_dT(:,:,:)
    call rism1d_closure_calcXvv_DT(this%closure,this%hvv_dT)
    xvv_dT => this%closure%xvv_dT
  end function rism1d_getSusceptibility_DT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns an NV X NV array of total excess coordination numbers excluding
!!!multiplicity
!!!IN:
!!!   this : rism1d object
!!!OUT:
!!!    nv*nv array (iv1,iv2) of the excess of site iv2 around iv1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getExNumber(this) result(exvv)
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: exvv(this%pot%nv,this%pot%nv)
    call rism_timer_start(this%thermoTimer)
    exvv = rism1d_closure_getExNumber(this%closure,this%hvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getExNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a pointer to an NR X NVV array of structure factors.  This memory 
!!!must be freed (preferably with safemem_dealloc) as it is not freed locally or
!!!after the object instance is destroyed.
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getStructFactor(this) result(svv)
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_, pointer :: svv(:,:)
    call rism_timer_start(this%thermoTimer)
    svv => rism1d_closure_getStructFactor(this%closure,this%hvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getStructFactor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a pointer to an NR X NVV array of the running site-site excess number.
!!!This is the excess number of a site within a given radius. The memory for 
!!!this pointer must be freed (preferably with safemem_dealloc) as it is not 
!!!freed locally or after the object instance is destroyed.
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getRunExNumber(this) result(exnvv)
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_, pointer :: exnvv(:,:,:)
    nullify(exnvv)
    call rism_timer_start(this%thermoTimer)
    exnvv => rism1d_closure_getRunExNumber(this%closure,this%gvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getRunExNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a pointer to an NR X NVV array of the running site-site number.
!!!This is the number of a site within a given radius. The memory for 
!!!this pointer must be freed (preferably with safemem_dealloc) as it is not 
!!!freed locally or after the object instance is destroyed.
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getRunNumber(this) result(nvv)
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_, pointer :: nvv(:,:,:)
    nullify(nvv)
    call rism_timer_start(this%thermoTimer)
    nvv => rism1d_closure_getRunNumber(this%closure,this%gvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getRunNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the pressure in [kT / A^3] of the system using the virial 
!!!path.  To convert to Pacals, for example, multiply by 
!!!1.d27 * kb * temperature
!!!where kb is Boltzmann's constant [j/K] and temperature is in [K]
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getPressureVirial(this) result(pressure)
    use safemem
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: pressure
    call rism_timer_start(this%thermoTimer)
    pressure = rism1d_closure_getPressureVirial(this%closure,this%gvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getPressureVirial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the pressure in [kT / A^3] of the system using the free energy 
!!!path.  To convert to Pacals, for example, multiply by 
!!!1.d27 * kb * temperature
!!!where kb is Boltzmann's constant [j/K] and temperature is in [K]
!!!IN:
!!!   this : rism1d object
!!!OUT:
!!!  Total pressure in kT / A^3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getPressureFE(this) result(pressure)
    use safemem
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: pressure
    call rism_timer_start(this%thermoTimer)
    pressure = rism1d_closure_getPressureFE(this%closure,this%gvv,this%cvv)
    call rism_timer_stop(this%thermoTimer)
   end function rism1d_getPressureFE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the freeEnergy in kT
!!!IN:
!!!   this : rism1d object
!!!OUT:
!!!   the excess free energy of the system in kT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getFreeEnergy(this) result(fe)
    use safemem
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: fe
    call rism_timer_start(this%thermoTimer)
    fe = rism1d_closure_getFreeEnergy(this%closure,this%gvv,this%cvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getFreeEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the partial molar volume of each species in A^3
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getPMV(this) result(pmv)
    use constants, only : pi
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: pmv(this%pot%nsp)
    call rism_timer_start(this%thermoTimer)
    pmv = rism1d_closure_getPMV(this%closure,this%cvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getPMV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site.
!!!IN:
!!!   this : rism1d object
!!!   o_form : specific form to use for the k-space contribution: 
!!!            'PR' Pettit-Rossky (default), 
!!!            'SC' Singer-Chandler, 'SM' Schmeer-Maurer.
!!!            S. J. Singer; D. Chandler. Molecular Physics, 1985, 55,
!!!            3, 621—625
!!!            B. M. Pettitt; P. J. Rossky. J. Chem. Phys. 1986, 15,
!!!            5836-5844
!!!            G. Schmeer and A. Maurer. Phys. Chem. Chem. Phys.,
!!!            2010, 12, 2407–2417
!!!OUT:
!!!   the excess chemical potential for each site in kT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getExChem(this, o_form) result(exchem)
    use safemem
    use rism_util, only : caseup
    implicit none
    type(rism1d), intent(inout) :: this
    character(len=*), optional, intent(in):: o_form
    character(len=2) form
    _REAL_ :: exchem(this%pot%nv)
    call rism_timer_start(this%thermoTimer)
    form = 'PR'
    if(present(o_form))then
       form = o_form
       call caseup(form)
    end if
    if(form .eq. 'PR')then
       exchem = rism1d_closure_getExChem_PR(this%closure,this%gvv,this%cvv)
    elseif(form .eq. 'SM')then
       exchem = rism1d_closure_getExChem_SM(this%closure,this%gvv,this%cvv)
    elseif(form .eq. 'SC')then
       exchem = rism1d_closure_getExChem_SC(this%closure,this%gvv,this%hvv,this%cvv)
    else
       call rism_report_error("'"//form//"' is not a valid excess chemical potential formula")
    end if
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getExChem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the ionic(?) excess chemical potential in kT for each site.
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getExChemIon(this) result(exchem)
    use safemem
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: exchem(this%pot%nv)
    call rism_timer_start(this%thermoTimer)
    exchem = rism1d_closure_getExChemIon(this%closure,this%gvv,this%cvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_getExChemIon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess solvation energy in kT for each site. This is
!!!the interaction of the solvent with the solute and the change in
!!!the solvent self interaction. Requires solution to both regular and
!!!temperature derivative correlation functions
!!!IN:
!!!   this : rism1d object
!!!OUT:
!!!   excess solvation energy in kT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_getSolvene(this) result(solvene)
    use safemem
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_ :: solvene(this%pot%nv)
    solvene = rism1d_closure_getSolvene(this%closure,this%gvv,this%pot%uvv,&
            this%cvv,this%gvv_dT,this%cvv_dT)
  end function rism1d_getSolvene

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Bvv (bridge function) for the closure used.The memory for 
!!!this pointer must be freed (preferably with safemem_dealloc) as it is not 
!!!freed locally or after the object instance is destroyed.
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_bvv(this) result(bvv)
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_, pointer :: bvv(:,:)
    call rism_timer_start(this%thermoTimer)
    bvv=>rism1d_closure_bvv(this%closure,this%gvv, this%cvv)
    call rism_timer_stop(this%thermoTimer)
  end function rism1d_bvv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory.  Makes object pristine
!!!IN:
!!!   this : rism1d object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_destroy(this)
    use safemem
    implicit none
    type(rism1d), intent(inout) :: this
    call rism_timer_start(this%timer)
    call rism1d_closure_destroy(this%closure)
    call rism1d_potential_destroy(this%pot)
    this%Mdiis_nvec=0
    if(safemem_dealloc(this%gvv) /=0) call rism_report_error("RISM1D: failed to deallocate GVV")
    if(safemem_dealloc(this%hvv) /=0) call rism_report_error("RISM1D: failed to deallocate HVV")
    if(safemem_dealloc(this%cvvWRK) /=0) call rism_report_error("RISM1D: failed to deallocate CVVWRK")
    if(safemem_dealloc(this%cvvresWRK) /=0) call rism_report_error("RISM1D: failed to deallocate CVVRESWRK")
    if(safemem_dealloc(this%gvv_dT) /=0) call rism_report_error("RISM1D: failed to deallocate GVVDT")
    if(safemem_dealloc(this%hvv_dT) /=0) call rism_report_error("RISM1D: failed to deallocate HVVDT")
    if(safemem_dealloc(this%cvvWRK_dT) /=0) call rism_report_error("RISM1D: failed to deallocate CVVDTWRK")
    if(safemem_dealloc(this%cvvresWRK_dT) /=0) call rism_report_error("RISM1D: failed to deallocate CVVDTRESWRK")
    nullify(this%cvv)
    nullify(this%cvvres)
    call rism_timer_destroy(this%potTimer)
    call rism_timer_destroy(this%fftTimer)
    call rism_timer_destroy(this%solventTimer)
    call rism_timer_destroy(this%thermoTimer)
    call rism_timer_destroy(this%solve3DRISMTimer)
    call rism_timer_destroy(this%single3DRISMsolutionTimer)
    call rism_timer_destroy(this%solvent_dTTimer)
    call rism_timer_destroy(this%solve3DRISM_dTTimer)
    call rism_timer_destroy(this%single3DRISMsolution_dTTimer)
    call rism_timer_destroy(this%fft_dTTimer)
    call rism_timer_destroy(this%selftestTimer)
    call rism_timer_stop(this%timer)
    call rism_timer_destroy(this%timer)
  end subroutine rism1d_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                               PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Check user parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sanity_check(this)
    implicit none
    type(rism1d), intent(inout) :: this
    
    !perhaps this can be moved to MDIIS code...

    if(this%mdiis_del <= 0)&
         call rism_report_error("MDIIS_MDIIS_DEL must be > 0")
    if(this%mdiis_nvec <=0)&
         call rism_report_error("MDIIS_mdiis_nvec must be > 0")

  end subroutine sanity_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Relaxing the RISM equation
!!!IN:
!!!   this : rism1d object
!!!   ksave : save restart frequency
!!!   progress : write progress frequency
!!!   maxstep : max number of steps
!!!   tolerance : target tolerance
!!!OUT:
!!!    .true. if converged.  .false. otherwise.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function  solve3DRISM (this,ksave,progress,maxstep,tolerance) result (converged)
    use rism_util, only : polynomialInterpolation_progressive, freeUnit
    use constants, only : pi
#ifdef __PGI
    use ieee_arithmetic, only : ieee_is_nan
#endif
    implicit none
    type(rism1d), intent(inout) :: this
    integer, intent(in) ::  maxstep,ksave,progress
    _REAL_, intent(in) :: tolerance
    logical :: converged
    integer :: unit
    logical ::  saved, start
    integer ::  nro, nvvo, ir,ivv, istep, kis
    _REAL_ ::  residual, r, err0

    !.................... k=0 extrapolation array size .....................
    _REAL_ ::  ep0(maxep0)

    _REAL_ :: k
    type(mdiis) :: mdiis_o

    call rism_timer_start(this%solve3DRISMTimer)

    call mdiis_new(mdiis_o,2,&
         this%mdiis_del,0d0, this%mdiis_restart)
    call mdiis_setData(mdiis_o,this%cvvWRK, this%cvvresWRK,&
         this%pot%nr*this%pot%nvv, this%Mdiis_nvec)
    this%cvv=>this%cvvWRK(:,:,mdiis_getWorkVector(mdiis_o))
    this%cvvres=>this%cvvresWRK(:,:,mdiis_getWorkVector(mdiis_o))

    if(ksave < -1) &
         call rism_report_error("KSAVE must be >= -1") 
    if(progress < 0) &
         call rism_report_error("PROGRESS must be >= 0") 
    if(maxstep <= 0) &
         call rism_report_error("MAXSTEP must be > 0") 
    if(tolerance <= 0) &
         call rism_report_error("TOLERANCE must be > 0") 

    converged = .false.
    unit = freeUnit()

    !.................. checking extrapolation grid size ...................
    if (maxep0 > this%pot%nr) &
       call rism_report_error('(a,i4,a,i4)','RXRISM: k->0 array size MaxEp0=',&
            maxep0,' > Nr=',this%pot%nr)

    !......................... reading saved r*Cvv .........................
    inquire (file=this%savefile,exist=saved)
    if (saved)  then
       call rism_report_message('reading saved Cvv file: '//trim(this%savefile))
       open (unit,file=this%savefile,form='unformatted',status='old')
       read (unit)  nro, nvvo
       if (nro /= this%pot%nr .OR. nvvo /= this%pot%nvv)  then
          close (unit)
          call rism_report_error('(a,i6,a,i4,a,i6,a,i4)', &
               'RXRISM: dimensions discrepancy with saved data' &
               //'               Nr         Nvv  ' &
               //'  actual    ',this%pot%nr,'       ',this%pot%nvv, &
               '   saved    ',nro,'       ',nvvo)
       endif
       read (unit)  this%cvv
       close (unit)

       !....................... initial guess of r*Cvv ........................
    else
       do ivv=1,this%pot%nvv
          do ir=2,this%pot%nr
             this%cvv(ir,ivv) = - this%pot%ulrvv(ir,ivv)
          enddo
       enddo

    endif

    !............................ relaxing RISM ............................
    call rism_report_message('relaxing RISM:')
    converged = .false.
    start = .true.

    do istep=1,maxstep

       ! The first element is never set but is read by the closure and
       ! MDIIS. This ensures that is it initialized properly. It
       ! doesn't change the result but does allow Valgrind to run
       ! cleanly.
       this%cvv(1,:)=0
       this%cvvres(1,:)=0
       this%gvv(1,:)=0
       this%hvv(1,:)=0
       !..................... one relaxation step of RISM .....................
       call  single3DRISMsolution (this,residual,tolerance, start,converged,mdiis_o)
       !................. terminate if residual is NaN ........................
       ! isnan is not standard, but GNU doesn't support Fortran 2003's
       ! ieee_arithmetic until version 5. Intel 15 supports both; use
       ! the PGI define to enable the most support for old compilers.
#ifdef __PGI
       if (ieee_is_nan(residual)) &
#else
       if (isnan(residual)) &
#endif
          call rism_report_error( 'RXRISM: residual is NaN' )

       !................. screen outputting relaxation steps ..................
       if (progress /= 0)  then
          if (converged .or. mod(istep,progress) == 0)  then
             call rism_report_message('(a,i4,a,1pe24.16,a,i3)',"step=",istep,&
                  "     Res=",residual,"     MDIIS=",getCurrentNVec(mdiis_o))
             call flush(6)
          endif
       endif

       !............ saving intermediate and last relaxation steps ............
       if (ksave == 0) then
          ! do nothing; this guards against division by 0 in remainder function
          ! because Fortran does not guarantee short circuit logical evaluation.
       else if (converged .OR. (ksave > 0 .AND. mod(istep,ksave) == 0))  then
          call rism_report_message('saving Cvv to file: '//trim(this%savefile)//' ...')
          open(unit,file=this%savefile,form='unformatted',status='unknown')
          write (unit)  this%pot%nr, this%pot%nvv
          write (unit)  this%cvv
          close (unit)
          call rism_report_message('done.')
       endif

       if(converged .and. this%closureID .eq. "MV0" &
            .and. .not. rism1d_closure_isCharged(this%closure))then
          converged = .false.
          call rism1d_closure_useCharged(this%closure,.true.)
          call rism_timer_stop(this%solve3DRISMTimer)
          call rism_timer_start(this%potTimer)
          call rism1d_potential_calc(this%pot,1d0)
          call rism_timer_stop(this%potTimer)
          call rism_timer_start(this%solve3DRISMTimer)
          call mdiis_reset(mdiis_o)
          this%cvv=>this%cvvWRK(:,:,mdiis_getWorkVector(mdiis_o))
          this%cvvres=>this%cvvresWRK(:,:,mdiis_getWorkVector(mdiis_o))
          do ivv=1,this%pot%nvv
             do ir=2,this%pot%nr
                this%cvv(ir,ivv) = - this%pot%ulrvv(ir,ivv)
             enddo
          enddo
       end if

       !............... exiting relaxation loop on convergence ................
       if (converged)  exit
    enddo
    if(.not. converged) &
       call rism_report_error('(a,i5)','RXRISM: reached steps limit Maxstep=',maxstep)

    !....................... eliminating R from Cvv ........................
    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          r = (ir-1)*this%pot%dr
          this%cvv(ir,ivv) = this%cvv(ir,ivv) / r
       enddo
    enddo

    !................. extrapolating Cvv(r=0) and Gvv(r=0) .................
    do ir=1,maxep0
       r = ir*this%pot%dr
       ep0(ir) = r
    enddo

    do ivv=1,this%pot%nvv
       call  polynomialInterpolation_progressive (ep0,this%cvv(2:1+maxep0,ivv),maxep0, 0.d0,this%cvv(1,ivv), err0)
    enddo

    do ivv=1,this%pot%nvv
       call  polynomialInterpolation_progressive (ep0,this%gvv(2:1+maxep0,ivv),maxep0, 0.d0,this%gvv(1,ivv), err0)
    enddo

    !....................... eliminating K from Hvv ........................
    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          k = (ir-1)*this%pot%dk
          this%hvv(ir,ivv) = this%hvv(ir,ivv) / k
       enddo
    enddo

    !....................... extrapolating Hvv(k=0) ........................
    do ir=1,maxep0
       k = ir*this%pot%dk
       ep0(ir) = k
    enddo

    do ivv=1,this%pot%nvv
       call  polynomialInterpolation_progressive (ep0,this%hvv(2:1+maxep0,ivv),maxep0, 0.d0,this%hvv(1,ivv), err0)
    enddo

    call mdiis_destroy(mdiis_o)
    call rism_timer_stop(this%solve3DRISMTimer)
  end function solve3DRISM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! One relaxation step for the direct correlation function
!!! in the RISM equation,
!!!   k*Hvv(k) = (k-(Wvv+Zvv)^(tr)*k*Cvv*Density_v)^(-1)
!!!              * k*(Wvv+Zvv)^(tr)*k*Cvv*(Wvv+Zvv) + k*Zvv
!!! where '^(tr)' means 'transpose', with the selected closure.
!!!OUT:
!!!    .true. if converged.  .false. otherwise.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  single3DRISMsolution (this,residual,tolerance,start,converged,mdiis_o)
    implicit none
#include "../xblas/f77/blas_namedconstants.fh"    
    type(rism1d), intent(inout) :: this
    _REAL_, intent(out) :: residual
    _REAL_, intent(in) :: tolerance
    logical, intent(out) ::  start, converged
    type(mdiis), intent(inout) :: mdiis_o

    integer :: err

    integer ::  ir, ivv, iv1,iv2,iv3
    _REAL_ ::  dind, r,t, arg, tvv, tvvr, al

    integer ::  indx(this%pot%nv)
    _REAL_ ::  ak(this%pot%nv,this%pot%nv), bk(this%pot%nv,this%pot%nv), wck(this%pot%nv,this%pot%nv)

    call rism_timer_start(this%single3DRISMsolutionTimer)
    !............................. FFT>K r*Cvv .............................
    call DCOPY(this%pot%nvv*this%pot%nr,this%cvv,1,this%cvvres,1)
    call DAXPY(this%pot%nvv*this%pot%nr,1d0,this%pot%ulrvv,1,this%cvvres,1)

    call rism_timer_start(this%fftTimer)
    do ivv=1,this%pot%nvv
       call  sinfti (this%cvvres(2:,ivv),this%pot%nr-1, this%pot%dr, +1)
    enddo
    call rism_timer_stop(this%fftTimer)

    call DAXPY(this%pot%nvv*this%pot%nr,-1d0,this%pot%ulkvv,1,this%cvvres,1)

    !................... enumerating K-space gridpoints ....................
    do ir=2,this%pot%nr
       t = (ir-1)*this%pot%dk

       !............. loading and symmetrizing k*Cvv at current K .............
       ivv = 0
       do iv2=1,this%pot%nv
          do iv1=1,iv2
             ivv = ivv + 1
             bk(iv1,iv2) = this%cvvres(ir,ivv)
             bk(iv2,iv1) = bk(iv1,iv2)
          enddo
       enddo

       !.................... getting (Wvv+Zvv)^(tr)*k*Cvv .....................
       if(this%extra_precision == 0)then
          call DGEMM('T','N',this%pot%nv,this%pot%nv,this%pot%nv,1d0,&
               this%pot%wzvv(ir,:,:),this%pot%nv,&
               bk,this%pot%nv,&
               0d0,wck,this%pot%nv)
       else
          call BLAS_DGEMM_X(BLAS_TRANS,BLAS_NO_TRANS,this%pot%nv,this%pot%nv,this%pot%nv,1d0,&
               this%pot%wzvv(ir,:,:),this%pot%nv,&
               bk,this%pot%nv,&
               0d0,wck,this%pot%nv,&
               BLAS_PREC_EXTRA)
       end if

       !............. getting Bk=k*(Wvv+Zvv)^(tr)*k*Cvv*(Wvv+Zvv) .............
       if(this%extra_precision == 0)then
          call DGEMM('N','N',this%pot%nv,this%pot%nv,this%pot%nv,t,&
               wck,this%pot%nv,&
               this%pot%wzvv(ir,:,:),this%pot%nv,&
               0d0,bk,this%pot%nv)
       else
          call BLAS_DGEMM_X(BLAS_NO_TRANS,BLAS_NO_TRANS,this%pot%nv,this%pot%nv,this%pot%nv,t,&
               wck,this%pot%nv,&
               this%pot%wzvv(ir,:,:),this%pot%nv,&
               0d0,bk,this%pot%nv,&
               BLAS_PREC_EXTRA)
       end if

       !................. getting (Wvv+Zvv)^(tr)*k*Cvv*Density_v ..................
       do iv2=1,this%pot%nv
          call DSCAL(this%pot%nv,this%pot%densityv(iv2),wck(:,iv2),1)
       enddo

       !.............. getting Ak=(k-(Wvv+Zvv)^(tr)*k*Cvv*Density_v) ..............
       ak=0
       do iv1=1,this%pot%nv
          ak(iv1,iv1) = t
       enddo
       call DAXPY(this%pot%nv**2, -1d0,wck,1,ak,1)

       !.................... calculating k*Hvv=Ak^(-1)*Bk .....................
       !Invert Ak using LU decomposition and multiply by Bk
       call DGESV(this%pot%nv,this%pot%nv,ak,this%pot%nv,indx,bk,this%pot%nv,err)
       if(err < 0)then
          err = err*(-1)
          call rism_report_error("Linear equation solver failed.")
       endif

       !......................... getting k*Hvv+k*Zvv .........................
       call DAXPY(this%pot%nv**2,1d0,this%pot%zkvv(ir,:,:),1,bk,1)

       !........................ unloading k*Hvv+k*Zvv ........................
       ivv = 0
       do iv2=1,this%pot%nv
          do iv1=1,iv2
             ivv = ivv + 1
             this%hvv(ir,ivv) = bk(iv1,iv2)
          enddo
       enddo
    enddo
    
    !............................. FFT>R k*Hvv .............................
    call DCOPY(this%pot%nvv*this%pot%nr,this%hvv,1,this%cvvres,1)
    call DAXPY(this%pot%nvv*this%pot%nr,-1d0,this%pot%hlkvv,1,this%cvvres,1)
    
    call rism_timer_start(this%fftTimer)
    do ivv=1,this%pot%nvv
       call  sinfti (this%cvvres(2:,ivv),this%pot%nr-1, this%pot%dr, -1)
    enddo
    call rism_timer_stop(this%fftTimer)

    call DAXPY(this%pot%nvv*this%pot%nr,1d0,this%pot%hlrvv,1,&
         this%cvvres,1)

    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          r = (ir-1)*this%pot%dr
          this%cvvres(ir,ivv) = this%cvvres(ir,ivv)/r
          this%cvv(ir,ivv) = this%cvv(ir,ivv)/r
       end do
    end do
    call rism1d_closure_gvv(this%closure,this%gvv, this%cvvres, this%cvv)
    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          r = (ir-1)*this%pot%dr
          this%cvvres(ir,ivv) = (this%gvv(ir,ivv)-1d0 - this%cvvres(ir,ivv))*r
          this%cvv(ir,ivv) = this%cvv(ir,ivv)*r
       end do
    end do
    !........................ performing MDIIS step ........................

    call mdiis_advance(mdiis_o,residual,converged,tolerance)
    this%cvv=>this%cvvWRK(:,:,mdiis_getWorkVector(mdiis_o))
    this%cvvres=>this%cvvresWRK(:,:,mdiis_getWorkVector(mdiis_o))
    call rism_timer_stop(this%single3DRISMsolutionTimer)
  end subroutine single3DRISMsolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Relaxing the RISM DT equation
!!!IN:
!!!   this : rism1d object
!!!   ksave : save restart frequency
!!!   progress : write progress frequency
!!!   maxstep : max number of steps
!!!   tolerance : target tolerance
!!!OUT:
!!!    .true. if converged.  .false. otherwise.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function  solve3DRISM_dT (this,ksave,progress,maxstep,tolerance) result (converged)
    use mdiis_c
    use rism_util, only : polynomialInterpolation_progressive, freeUnit
    use constants, only : pi
    implicit none
    type(rism1d), intent(inout) :: this
    integer, intent(in) ::  maxstep,ksave,progress
    _REAL_, intent(in) :: tolerance
    logical :: converged
    type(mdiis), save :: mdiis_o
    integer :: unit
    logical ::  saved, start
    integer ::  nro, nvvo, ir,ivv, istep, kis
    _REAL_ ::  residual, r, err0

    !.................... k=0 extrapolation array size .....................
    _REAL_ ::  ep0(maxep0)

    _REAL_ :: k

    !xvv :: site-site solvent suseptibility from standard calculation
    _REAL_, pointer :: xvv(:,:,:)
    !mdiis_method :: implementation of MDIIS to use
    integer :: mdiis_method =2

    ! Fortran does not guarantee short circuit logical expressions:
    logical is_remainder0

    call rism_timer_start(this%solve3DRISM_dTTimer)
!    call rism_report_message('retrieving Xvv ...')
    call rism1d_closure_calcXvv(this%closure,this%hvv)
    xvv => this%closure%xvv



    converged = .false.
    call mdiis_new(mdiis_o,mdiis_method,&
         this%mdiis_del,tolerance, this%mdiis_restart)

    call mdiis_setData(mdiis_o,this%cvvWRK_dT,this%cvvresWRK_dT,&
         this%pot%nr*this%pot%nvv, this%Mdiis_nvec)
    unit = freeUnit()

    !.................. checking extrapolation grid size ...................
    if (maxep0 > this%pot%nr) &
       call rism_report_error('(a,i4,a,i4)','RXRISMDT: k->0 array size MaxEp0=',&
            maxep0,' > Nr=',this%pot%nr)

    !.... FFT>K r*Cvv
    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          r = (ir-1)*this%pot%dr
          this%cvvres(ir,ivv) = this%cvv(ir,ivv)*r + this%pot%ulrvv(ir,ivv)
       enddo
    enddo
    do ivv=1,this%pot%nvv
       call  sinfti (this%cvvres(2:,ivv),this%pot%nr-1, this%pot%dr, +1)
    enddo
    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          this%cvvres(ir,ivv) = this%cvvres(ir,ivv) - this%pot%ulkvv(ir,ivv)
       enddo
    enddo

    !... initial guess of r*dCvv/dt
    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          this%cvv_dT(ir,ivv) = this%pot%ulrvv(ir,ivv)
       enddo
    enddo

    !............................ relaxing RISM ............................
    call rism_report_message('relaxing RISM DT:')
    converged = .false.
    start = .true.


    do istep=1,maxstep

       !the first element is never set but is read by the closure and MDIIS.  This
       !ensures that is it initialized properly.  It doesn't change the result but
       !does allow Valgrind to run cleanly
       this%cvv_dT(1,:)=0
       this%cvv_dTres(1,:)=0
       !..................... one relaxation step of RISM .....................
       call  single3DRISMsolution_dT (this,xvv,residual,tolerance, start,converged, mdiis_o)
       !................. screen outputting relaxation steps ..................
       if (progress /= 0)  then
          is_remainder0 = .false.
          if ( ksave > 0 ) is_remainder0 = mod(istep,ksave) == 0
          if (converged .OR. mod(istep,progress) == 0 .OR. is_remainder0)  then
             call rism_report_message('(a,i4,a,1pe24.16,a,i3)',"step=",istep,&
                  "     Res=",residual,"     MDIIS=",getCurrentNVec(mdiis_o))
             call flush(6)
          endif
       endif

       !............... exiting relaxation loop on convergence ................
       if (converged)  exit
    enddo
    if(.not. converged) &
       call rism_report_error('(a,i5)','RXRISMDT: reached steps limit Maxstep=',maxstep)

    !....................... eliminating R from dCvv/dt ........................
    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          r = (ir-1)*this%pot%dr
          this%cvv_dT(ir,ivv) = this%cvv_dT(ir,ivv) / r
       enddo
    enddo

    !................. extrapolating dCvv/dt(r=0) and dGvv/dt(r=0) .................
    do ir=1,maxep0
       r = ir*this%pot%dr
       ep0(ir) = r
    enddo

    do ivv=1,this%pot%nvv
       call  polynomialInterpolation_progressive (ep0,this%cvv_dT(2:1+maxep0,ivv),maxep0, 0.d0,this%cvv_dT(1,ivv), err0)
    enddo

    do ivv=1,this%pot%nvv
       call  polynomialInterpolation_progressive (ep0,this%gvv_dT(2:1+maxep0,ivv),maxep0, 0.d0,this%gvv_dT(1,ivv), err0)
    enddo

    !....................... eliminating K from dHvv/dt ........................
    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          k = (ir-1)*this%pot%dk
          this%hvv_dT(ir,ivv) = this%hvv_dT(ir,ivv) / k
       enddo
    enddo

    !....................... extrapolating dHvv(k=0)/dt ........................
    do ir=1,maxep0
       k = ir*this%pot%dk
       ep0(ir) = k
    enddo

    do ivv=1,this%pot%nvv
       call  polynomialInterpolation_progressive (ep0,this%hvv_dT(2:1+maxep0,ivv),maxep0, 0.d0,this%hvv_dT(1,ivv), err0)
    enddo

    call mdiis_destroy(mdiis_o)
    call rism_timer_stop(this%solve3DRISM_dTTimer)
  end function solve3DRISM_dT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! One relaxation step for the direct correlation function        *
!!! in the RISM equation with the HNC or PLHNC closure,          *
!!! k*Hvv(k) = (k-(Wvv+Zvv)^(tr)*k*Cvv*Density_v)^(-1)            *
!!! * k*(Wvv+Zvv)^(tr)*k*Cvv*(Wvv+Zvv) + k*Zvv ,   *
!!! Gvv(r) = exp(-b*Uvv + Hvv - Cvv)  or                      *
!!! Gvv(r) = exp(X), X<0                                      *
!!! = 1 + X,  X>0                                      *
!!! X = - b*Uvv + Hvv - Cvv                                 *
!!!OUT:
!!!    .true. if converged.  .false. otherwise.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  single3DRISMsolution_dT (this,xvv,residual,tolerance,start,converged,mdiis_o)
    use mdiis_c
    implicit none
    type(rism1d), intent(inout) :: this
    _REAL_, intent(out) :: residual
    _REAL_, intent(in) :: tolerance
    logical, intent(out) ::  start, converged
    type(mdiis), intent(inout) :: mdiis_o

    integer :: err

    integer ::  ir, ivv, iv1,iv2,iv3
    _REAL_ ::  dind, r,t, arg, tvv, tvvr, al

    integer ::  indx(this%pot%nv)
    _REAL_ ::  ak(this%pot%nv,this%pot%nv), bk(this%pot%nv,this%pot%nv), wck(this%pot%nv,this%pot%nv)

    _REAL_, intent(in) :: xvv(:,:,:)
    _REAL_ :: ck(this%pot%nv,this%pot%nv), dzk(this%pot%nv,this%pot%nv), pdz(this%pot%nv,this%pot%nv)
    _REAL_ :: wck_dT(this%pot%nv,this%pot%nv), dck(this%pot%nv,this%pot%nv)

    call rism_timer_start(this%single3DRISMsolution_dTTimer)
    !............................. FFT>K r*dCvv/dt .............................
    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          this%cvv_dTres(ir,ivv) = this%cvv_dT(ir,ivv) - this%pot%ulrvv(ir,ivv)
       enddo
    enddo

    call rism_timer_start(this%fft_dTTimer)
    do ivv=1,this%pot%nvv
       call  sinfti (this%cvv_dTres(2:,ivv),this%pot%nr-1, this%pot%dr, +1)
    enddo
    call rism_timer_stop(this%fft_dTTimer)

    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          this%cvv_dTres(ir,ivv) = this%cvv_dTres(ir,ivv) + this%pot%ulkvv(ir,ivv)
       enddo
    enddo

    !................... enumerating K-space gridpoints ....................
    do ir=2,this%pot%nr
       t = (ir-1)*this%pot%dk

       !............. loading and symmetrizing k*Cvv and k*dCvv/dt at current K .............
       ivv = 0
       do iv2=1,this%pot%nv
          do iv1=1,iv2
             ivv = ivv + 1
             bk(iv1,iv2) = this%cvvres(ir,ivv)
             bk(iv2,iv1) = bk(iv1,iv2)
             ck(iv1,iv2) = this%cvv_dTres(ir,ivv)
             ck(iv2,iv1) = ck(iv1,iv2)
          enddo
       enddo

       !.... getting d[Wvv+Zvv]/dt*k and dZvv/dt 
       if(this%pot%pcdiel == 0.d0)then
         dzk = 0.d0
         pdz = 0.d0
       else
         do iv2=1,this%pot%nv
            do iv1=1,this%pot%nv
               dzk(iv1,iv2) = (this%pot%pcdiel+3.d0)/this%pot%pcdiel *  this%pot%zkvv(ir,iv1,iv2)
               pdz(iv1,iv2) = this%pot%densityv(iv1) * dzk(iv1,iv2) / t
            enddo
         enddo
       endif

       !.................... getting (Wvv+Zvv)^(tr)*k*Cvv .....................
       do iv2=1,this%pot%nv
          do iv1=1,this%pot%nv
             wck(iv1,iv2) = 0.d0
             do iv3=1,this%pot%nv
                wck(iv1,iv2) = wck(iv1,iv2) + this%pot%wzvv(ir,iv3,iv1)*bk(iv3,iv2)
             enddo
          enddo
       enddo

       !.................... getting (Wvv+Zvv)^(tr)*k*dCvv/dt .....................
       do iv2=1,this%pot%nv
          do iv1=1,this%pot%nv
             wck_dT(iv1,iv2) = 0.d0
             do iv3=1,this%pot%nv
                wck_dT(iv1,iv2) = wck_dT(iv1,iv2) + this%pot%wzvv(ir,iv3,iv1)*ck(iv3,iv2)
             enddo
          enddo
       enddo

       !.................... getting (d[Wvv+Zvv]/dt)^(tr)*k*Cvv .....................
       do iv2=1,this%pot%nv
          do iv1=1,this%pot%nv
             dck(iv1,iv2) = 0.d0
             do iv3=1,this%pot%nv
                dck(iv1,iv2) = dck(iv1,iv2) + pdz(iv3,iv1) * bk(iv3,iv2)
             enddo
          enddo
       enddo

       !..... getting Bk = ( [(d[Wvv+Zvv]/dt)^(tr)*k*Cvv + (Wvv+Zvv)^(tr)*k*dCvv/dt]*Xvv + dZvv/dt ) * k
       do iv2=1,this%pot%nv
          do iv1=1,this%pot%nv
             bk(iv1,iv2) = 0.d0
             do iv3=1,this%pot%nv
                bk(iv1,iv2) = bk(iv1,iv2) + (dck(iv1,iv3)+wck_dT(iv1,iv3))*xvv(ir,iv3,iv2)
             enddo
             bk(iv1,iv2) = t*(bk(iv1,iv2) + dzk(iv1,iv2))
          enddo
       enddo

       !................. getting (Wvv+Zvv)^(tr)*k*Cvv*Density_v ..................
       do iv2=1,this%pot%nv
          do iv1=1,this%pot%nv
             wck(iv1,iv2) = wck(iv1,iv2)*this%pot%densityv(iv2)
          enddo
       enddo

       !.............. getting Ak=(k-(Wvv+Zvv)^(tr)*k*Cvv*Density_v) ..............
       do iv2=1,this%pot%nv
          do iv1=1,this%pot%nv
             ak(iv1,iv2) = - wck(iv1,iv2)
          enddo
       enddo

       do iv1=1,this%pot%nv
          ak(iv1,iv1) = t + ak(iv1,iv1)
       enddo

       !.................... calculating k*dHvv/dt=Ak^(-1)*Bk .....................
       call DGETRF(this%pot%nv,this%pot%nv,ak,this%pot%nv,indx,err)
       if(err > 0)then
          call rism_report_error("LU-factorization failed.  U = 0")
       elseif(err<0)then
          err = err*(-1)
          call rism_report_error("LU-factorization failed.")
       endif
       call DGETRS('N', this%pot%nv, this%pot%nv,ak,this%pot%nv,indx,bk,this%pot%nv,err)
       if(err < 0)then
          err = err*(-1)
          call rism_report_error("Linear equation solver failed.")
       endif

       !........................ unloading k*dHvv/dt ........................
       ivv = 0
       do iv2=1,this%pot%nv
          do iv1=1,iv2
             ivv = ivv + 1
             this%hvv_dT(ir,ivv) = bk(iv1,iv2)
          enddo
       enddo

    enddo

    !............................. FFT>R k*dHvv/dt .............................
    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          this%cvv_dTres(ir,ivv) = this%hvv_dT(ir,ivv) + this%pot%hlkvv(ir,ivv)
       enddo
    enddo
    call rism_timer_start(this%fft_dTTimer)
    do ivv=1,this%pot%nvv
       call  sinfti (this%cvv_dTres(2:,ivv),this%pot%nr-1, this%pot%dr, -1)
    enddo
    call rism_timer_stop(this%fft_dTTimer)

    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          this%cvv_dTres(ir,ivv) = this%cvv_dTres(ir,ivv) - this%pot%hlrvv(ir,ivv)
       enddo
    enddo

    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          r = (ir-1)*this%pot%dr
          this%cvv_dTres(ir,ivv) = this%cvv_dTres(ir,ivv)/r
          this%cvv_dT(ir,ivv) = this%cvv_dT(ir,ivv)/r
       end do
    end do
    call rism1d_closure_gvv_dT(this%closure,this%gvv_dT, this%gvv, this%cvv, this%cvv_dTres, this%cvv_dT)
    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          r = (ir-1)*this%pot%dr
          this%cvv_dTres(ir,ivv) = (this%gvv_dT(ir,ivv) - this%cvv_dTres(ir,ivv))*r
          this%cvv_dT(ir,ivv) = this%cvv_dT(ir,ivv)*r
       end do
    end do
    !........................ performing MDIIS step ........................

    call mdiis_advance(mdiis_o,residual,converged)
    this%cvv_dT=>this%cvvWRK_dT(:,:,mdiis_getWorkVector(mdiis_o))
    this%cvv_dTres=>this%cvvresWRK_dT(:,:,mdiis_getWorkVector(mdiis_o)) 
    call rism_timer_stop(this%single3DRISMsolution_dTTimer)
 end subroutine single3DRISMsolution_dT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Read modification for nonbond file
!!! this: rism1d object
!!! nbfix: the filename of nonbond_fix file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_readNBFix(this,nbfix)
    implicit none 
    type(rism1d), intent(inout) :: this
    character(len=*), intent(in) :: nbfix

    call rism1d_potential_initNBFix(this%pot,nbfix)
  end subroutine rism1d_readNBFix



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Read modification for nonbond file 12-6-4
!!! this: rism1d object
!!! nbfix: the filename of nonbond_fix file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_readNBFixtsf(this,nbfix)
    implicit none
    type(rism1d), intent(inout) :: this
    character(len=*), intent(in) :: nbfix

    call rism1d_potential_initNBFixtsf(this%pot,nbfix)
  end subroutine rism1d_readNBFixtsf
end module rism1d_c
