#include "../include/dprec.fh"
#include "nfe-config.h"
#include "../include/assert.fh"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Open input files and read cntrl namelist.
subroutine mdread1()
#undef API
#ifdef API
#  define FATAL_ERROR ierr = 1; return
#else
#  define FATAL_ERROR call mexit(6, 1)
#endif /* API */

   use file_io_dat
   use lmod_driver, only : read_lmod_namelist
   use constants, only : RETIRED_INPUT_OPTION, zero, one, two, three, seven, &
                         eight, NO_INPUT_VALUE_FLOAT, NO_INPUT_VALUE
   use md_scheme, only: ntt, gamma_ln
   use les_data, only : temp0les
   use stack, only: lastist,lastrst
   use nmr, only: echoin
   use sgld, only : isgld, isgsta,isgend,nsgsize, &
                    tsgavg,sgft,sgff,sgfg,tsgavp
   use emap,only: temap,gammamap
#ifdef DSSP
   use dssp, only: idssp
#endif /* DSSP */

#ifndef API
   use xray_globals_module, only: xray_active
#endif /* API */
#ifdef MPI /* SOFT CORE */
   use softcore, only : scalpha,scbeta,ifsc,scmask,logdvdl,dvdl_norest,dynlmb, &
                        sceeorder, tishake
   use mbar, only : ifmbar, bar_intervall, bar_l_min, bar_l_max, bar_l_incr
   use remd, only  : rem
#endif /* MPI */
   ! Parameter for LIE module
   use linear_response, only: ilrt, lrt_interval, lrtmask
#  ifndef API
   use sander_rism_interface, only: xvvfile, guvfile, huvfile, cuvfile,&
        uuvfile, quvFile, chgDistFile,  &
        excessChemicalPotentialfile, solvationEnergyfile, entropyfile, &
        solventPotentialEnergyfile
#  endif /* API */
   use sander_rism_interface, only: rismprm
   use nfe_sander_proxy, only: infe
   implicit none
#  include "box.h"
#  include "def_time.h"
#  include "ew_cntrl.h"
#  include "ew_pme_recip.h"
#  include "../include/md.h"
#  include "../include/memory.h"
#  include "nmr.h"
#  include "tgtmd.h"
#  include "multitmd.h"
#  include "ew_erfc_spline.h"

   character(len=4) watdef(4),watnam,owtnm,hwtnm1,hwtnm2
   _REAL_      dele
   integer     ierr  ! internal error reporting for sander API
   integer     imcdo
   integer     itotst
   integer     inerr
   logical     mdin_cntrl, mdin_lmod, mdin_qmmm  ! true if namelists are in mdin
   integer     mxgrp
   integer     iemap
   _REAL_      dtemp  ! retired
   _REAL_      dxm  ! retired
   _REAL_      heat  ! retired
   _REAL_      timlim ! retired
   character(11), parameter :: error_hdr =      '| ERROR:   '

! pmemd functionality that is not supported in sander:
   character(len=256) :: scmask1, scmask2, timask1, timask2

#ifdef API
   ! Input options passed in to the setup API routine
   type(sander_input), intent(in) :: input_options
#else
   integer     ifind
   character(len=8) date
   character(len=10) time
   character(len=512) :: char_tmp_512
#endif /* API */

   integer irism, stat
   character(len=256) line
   character(len=8) periodicPotential

!  N.B.: If you make changes to this namelist, you also need to make
!        corresponding changes in ./sander.h and in
!         ../pysander/sander/src/{pysandermodule,pysandermoduletypes}.c

   namelist /cntrl/ irest,ibelly, &
         ntx,ntxo,ntcx,ig,tempi, &
         ntb,temp0,tautp, &
         ntp,pres0,comp,barostat,mcbarint, &
         nscm,nstlim,t,dt, &
         ntc,ntcc,nconp,tol,ntf,ntn,nsnb, &
         cut,dielc, &
         ntpr,ntwx,ntwv,ntwe,ntwf,ntave,ntpp,ioutfm, &
         ntr,nrc,ntrx,taur,nmropt, &
         ivcap,cutcap,xcap,ycap,zcap,fcap, &
         xlorth,ylorth,zlorth,xorth,yorth,zorth,forth, &
         imin,drms,dele,dx0, &
         pencut,ipnlty,iscale,scalm,noeskp, &
         maxcyc,ncyc,ntmin,vlimit, &
         mxsub,jfastw,watnam,owtnm,hwtnm1,hwtnm2, iesp, &
         skmin, skmax, vv,vfac, tmode, &
         ntt, gamma_ln, &
         iemap,gammamap, &
         isgld,isgsta,isgend,nsgsize,tsgavg,sgft,sgff,sgfg,tsgavp,&
         jar, isynctraj, &
         numexchg, repcrd, numwatkeep, hybridgb, reservoir_exchange_step, &
         ntwprt,tausw, &
         ntwr,iyammp,imcdo, &
         plumed,plumedfile, &
         igb,alpb,Arad,rgbmax,saltcon,offset,gbsa,vrand, &
         surften,nrespa,nrespai,extdiel,intdiel, &
         cut_inner,icfe,clambda,klambda, rbornstat,lastrst,lastist,  &
         itgtmd,tgtrmsd,tgtmdfrc,tgtfitmask,tgtrmsmask, dec_verbose, &
         temp0les,restraintmask,restraint_wt,bellymask, &
         noshakemask,crgmask, &
         mask_from_ref, &
         rdt, &
         ievb, profile_mpi, &
         ipb, inp, nkija, idistr, sinrtau, &
         gbneckscale, &
         gbalphaH,gbbetaH,gbgammaH, &
         gbalphaC,gbbetaC,gbgammaC, &
         gbalphaN,gbbetaN,gbgammaN, &
         gbalphaOS,gbbetaOS,gbgammaOS, &
         gbalphaP,gbbetaP,gbgammaP, &
         Sh,Sc,Sn,So,Ss,Sp, &
         lj1264, fswitch, &
         infe, baroscalingdir, &
#ifdef MPI /* SOFT CORE */
         scalpha, scbeta, ifsc, scmask, logdvdl, dvdl_norest, dynlmb, &
         sceeorder, &
         ifmbar, bar_intervall, bar_l_min, bar_l_max, bar_l_incr, tishake, &
#endif
         ilrt, lrt_interval, lrtmask, &
#ifdef DSSP
         idssp, &
#endif
         irism,&
         vdwmodel, & ! mjhsieh - the model used for van der Waals
         ! retired:
         dtemp, dxm, heat, timlim, &
         ! pmemd functionality that is not supported in sander:
         scmask1, scmask2, timask1, timask2

   ! Define default water residue name and the names of water oxygen & hydrogens

   data watdef/'WAT ','O   ','H1  ','H2  '/

   !     ----- READ THE CONTROL DATA AND OPEN DIFFERENT FILES -----

#ifndef API /* If this is NOT the API */
   if (mdout /= "stdout" ) &
         call amopen(6,mdout,owrite,'F','W')
   call amopen(5,mdin,'O','F','R')
   write(6,9308)
   call date_and_time( DATE=date, TIME=time )
   write(6,'(12(a))') '| Run on ', date(5:6), '/', date(7:8), '/',  &
        date(1:4), ' at ', time(1:2), ':', time(3:4), ':', time(5:6)

   ! Write the path of the current executable and working directory
   call get_command_argument(0, char_tmp_512)
   write(6,'(/,a,a)') '|   Executable path: ', trim(char_tmp_512)
   call getcwd(char_tmp_512)
   write(6,'(a,a)') '| Working directory: ', trim(char_tmp_512)
! Write the hostname if we can get it from environment variable
! Note: get_environment_variable is part of the F2003 standard but seems
!       to be supported by GNU, Intel, IBM and Portland (2010+) compilers
   call get_environment_variable("HOSTNAME", char_tmp_512, inerr)
   if (inerr .eq. 0) then
     write(6,'(a,a,/)') '|          Hostname: Unknown'
   else
     write(6,'(a,a,/)') '|          Hostname: ', trim(char_tmp_512)
   end if

   if (owrite /= 'N') write(6, '(2x,a)') '[-O]verwriting output'

   ! Echo the file assignments to the user:

   write(6,9700) 'MDIN'   ,mdin(1:70)  , 'MDOUT' ,mdout(1:70) , &
         'INPCRD' ,inpcrd(1:70), 'PARM'  ,parm(1:70)  , &
         'RESTRT',restrt(1:70) , 'REFC'  ,refc(1:70)  , &
         'MDVEL' ,mdvel(1:70)  , 'MDFRC' ,mdfrc(1:70) , &
         'MDEN'   ,mden(1:70)  , &
         'MDCRD' ,mdcrd(1:70)  , 'MDINFO' ,mdinfo(1:70), &
         'MTMD'  ,mtmd(1:70)   , 'INPDIP', inpdip(1:70), &
         'RSTDIP', rstdip(1:70), 'INPTRAJ', inptraj(1:70)
#  ifdef MPI
   write(6,9702) 'REMLOG',     trim(remlog), &
                 'REMTYPE',    trim(remtype), &
                 'REMSTRIP',   trim(remstripcoord), &
                 'SAVEENE',    trim(saveenefile), &
                 'CLUSTERINF', trim(clusterinfofile), &
                 'RESERVOIR',  trim(reservoirname), &
                 'REMDDIM',    trim(remd_dimension_file)
#  endif
   if (len_trim(xvvfile) > 0) &
        write(6,9701) 'Xvv', trim(xvvfile)
   if (len_trim(guvfile) > 0) &
        write(6,9701) 'Guv', trim(Guvfile)
   if (len_trim(huvfile) > 0) &
        write(6,9701) 'Huv', trim(Huvfile)
   if (len_trim(cuvfile) > 0) &
        write(6,9701) 'Cuv', trim(Cuvfile)
   if (len_trim(uuvfile) > 0) &
        write(6,9701) 'Uuv', trim(Uuvfile)
   if (len_trim(quvfile) > 0) &
        write(6,9701) 'Quv', trim(Quvfile)
   if (len_trim(chgDistfile) > 0) &
        write(6,9701) 'ChgDist', trim(chgDistfile)
   if (len_trim(excessChemicalPotentialfile) > 0) &
        write(6,9701) 'ExChem', trim(excessChemicalPotentialfile)
   if (len_trim(solvationEnergyfile) > 0) &
        write(6,9701) 'SolvEne', trim(solvationEnergyfile)
   if (len_trim(entropyfile) > 0) &
        write(6,9701) 'Entropy', trim(entropyfile)
   if (len_trim(solventPotentialEnergyfile) > 0) &
        write(6,9701) 'PotUV', trim(solventPotentialEnergyfile)

   ! Echo the input file to the user:
   call echoin(5,6)
   !     ----- READ DATA CHARACTERIZING THE MD-RUN -----
   read(5,'(a80)') title
   !       ----read input in namelist format, first setting up defaults
#endif /* ifndef API */

   ierr = 0

   dtemp = RETIRED_INPUT_OPTION
   dxm   = RETIRED_INPUT_OPTION
   heat  = RETIRED_INPUT_OPTION
   timlim = RETIRED_INPUT_OPTION
! pmemd functionality that is not supported in sander:
   scmask1 = ''
   scmask2 = ''
   timask1 = ''
   timask2 = ''

   irest = 0
   ibelly = 0
   mask_from_ref = 0 ! Antoine Marion
   ipol = RETIRED_INPUT_OPTION
   iesp = 0
   ntx = 1
   ntxo = NO_INPUT_VALUE
   ig = -1
   tempi = ZERO
   ntb = NO_INPUT_VALUE
   temp0 = 300.0d0
   ntt = 3    ! was ithermostat=1
   gamma_ln = 5.d0  ! was therm_par = 5.0d0
   sinrtau = 0.05
   plumed = 0
   plumedfile = 'plumed.dat'
#ifdef LES
   ! alternate temp for LES copies, if negative then use single bath
   ! single bath not the same as 2 baths with same target T
   temp0les = -ONE
#endif
   rdt = 0

   tautp = ONE
   ntp = 0
   fswitch = -1.d0
   barostat = 2
   mcbarint = 100
   pres0 = ONE
   comp = 44.6d0
   npscal = 1
   nscm = 1000
   nstlim = 1
   t = ZERO
   dt = 0.001d0
   ntc = 1
   tol = 0.00001
   ntf = 1
   nsnb = NO_INPUT_VALUE
   cut =  NO_INPUT_VALUE_FLOAT
   dielc = ONE
   ntpr = 50
!RCW - Amber 16 change ntwr to equal nstlim
!0 here is to detect if user set it or not.
   ntwr = 0
   ntwx = 0
   ntwv = 0
   ntwf = 0
   ntwe = 0
   ipb = 0
   inp = 2
   irism = 0
   ntave = 0
   ioutfm = 1
   ntr = 0
   ntrx = 1
   ivcap = 0
   natcap = 0
   fcap = 1.5d0
   cutcap = 0.0d0
   xcap = 0.0d0
   ycap = 0.0d0
   zcap = 0.0d0
   forth = 1.5d0
   xlorth = -1.0d0
   ylorth = -1.0d0
   zlorth = -1.0d0
   xorth = 47114711.0d0
   yorth = 47114711.0d0
   zorth = 47114711.0d0
   numexchg = 0
   repcrd   = 1
   lj1264 = -1

   profile_mpi = 0 !whether to write profile_mpi timing file - default = 0 (NO).

   ! number of waters to keep for hybrid model,
   ! numwatkeep: the number of closest
   ! waters to keep. close is defined as close to non-water.
   ! for simulations with ions, ions should be stripped too
   ! or at least ignored in the "closest" calculation. this
   ! is not currently done.

   ! if it stays at -1 then we keep all waters
   ! 0 would mean to strip them all

    numwatkeep=-1

   ! hybridgb: gb model to use with hybrid REMD.
   hybridgb=0

   ! Reservoir Exchange Frequency
   reservoir_exchange_step=2

   ! carlos targeted MD, like ntr

   itgtmd=0
   tgtrmsd=0.
   tgtmdfrc=0.
   tgtfitmask=''
   tgtrmsmask=''

   ! Flag for writing on-step velocities and coordinates (instead of
   ! +half step velocities and +1 step coords which is the default).
   isynctraj=0

   pencut = 0.1d0
   taumet = 0.0001d0
   omega = 500.0d0
   ipnlty = 1
   scalm = 100.0d0
   iscale = 0
   noeskp = 1
   nmropt = 0
   jar = 0
   tausw = 0.1d0
   imin = 0
   isftrp = 0
   rwell = ONE
   maxcyc = 1
   ncyc = 10
   ntmin = 3
   dx0 = 0.01d0
   drms = 1.0d-4
   vlimit = 20.0d0
   mxsub = 1
   jfastw = 0
   watnam = '    '
   owtnm =  '    '
   hwtnm1 = '    '
   hwtnm2 = '    '
   ntwprt = 0
   igb = 0
   alpb = 0
   Arad = 15.0d0
   rgbmax = 25.d0
   saltcon = ZERO

   !  default offset depends on igb value, and users need to
   !  be able to modify it, so we need to set a dummy value. if it's still the
   !  dummy after we read the namelist, we set the default based on igb. if not,
   !  we leave it at what the user set.
   !  best solution would be to create a GB namelist.
   offset = -999999.d0
   gbneckscale = -999999.d0

   iyammp = 0
   imcdo = -1
   gbsa = 0
   vrand=1000
   surften = 0.005d0
   nrespa = 1
   nrespai = 1
   irespa = 1
   extdiel = 78.5d0
   intdiel = ONE
   gbgamma = ZERO
   gbbeta = ZERO
   gbalpha = ONE

   ! set default parameters for igb = 8
   ! NOTE THAT NONE OF THESE ARE USED UNLESS IGB=8, SO USERS SHOULD NOT EVEN SET
   ! THEM
   gbalphaH = 0.788440d0
   gbbetaH = 0.798699d0
   gbgammaH = 0.437334d0
   gbalphaC = 0.733756d0
   gbbetaC = 0.506378d0
   gbgammaC = 0.205844d0
   gbalphaN = 0.503364d0
   gbbetaN = 0.316828d0
   gbgammaN = 0.192915d0
   gbalphaOS = 0.867814d0
   gbbetaOS = 0.876635d0
   gbgammaOS = 0.387882d0
   gbalphaP   = 0.418365d0 !!Use gbneck2nu parms for P
   gbbetaP    = 0.290054d0
   gbgammaP   = 0.1064245d0
   !scaling parameters below will only be used for igb=8.
   ! the actual code does not use these variables, it uses X(l96)
   ! if igb=8, we will use these to set the X(l96) array.
   Sh = 1.425952d0
   Sc = 1.058554d0
   Sn = 0.733599d0
   So = 1.061039d0
   Ss = -0.703469d0
   Sp = 0.5d0          !P parameters are not optimized for protein
   ! update gbneck2nu pars
   ! using offset and gb_neckscale parameters from GB8-protein for gbneck2nu
   ! Scaling factors
   ! name of variables are different from sander's igb8
   ! we use below names in pmemd
   ! gbneck2nu
   screen_hnu = 1.696538d0
   screen_cnu = 1.268902d0
   screen_nnu = 1.4259728d0
   screen_onu = 0.1840098d0
   screen_pnu = 1.5450597d0
   !alpha, beta, gamma for each atome element
   gb_alpha_hnu = 0.537050d0
   gb_beta_hnu = 0.362861d0
   gb_gamma_hnu = 0.116704d0
   gb_alpha_cnu = 0.331670d0
   gb_beta_cnu = 0.196842d0
   gb_gamma_cnu = 0.093422d0
   gb_alpha_nnu = 0.686311d0
   gb_beta_nnu = 0.463189d0
   gb_gamma_nnu = 0.138722d0
   gb_alpha_osnu = 0.606344d0
   gb_beta_osnu = 0.463006d0
   gb_gamma_osnu = 0.142262d0
   gb_alpha_pnu = 0.418365d0
   gb_beta_pnu = 0.290054d0
   gb_gamma_pnu = 0.1064245d0
   ! End gbneck2nu

   iconstreff = 0
   cut_inner = EIGHT
   icfe = 0
   clambda = ZERO
   klambda = 1
   ievb = 0
   rbornstat = 0
   ! added a flag to control output of BDC/SDC synonymous with mmpbsa.py's
   ! version of the same variable.
   dec_verbose = 3
   lastrst = 1
   lastist = 1
   restraintmask=''
   restraint_wt = ZERO
   bellymask=''
   noshakemask=''
   crgmask=''

   skmin = 50 !used by neb calculation
   skmax = 100 !used by neb calculation
   vv = 0 !velocity verlet -- off if vv/=1
   vfac = 0 !velocity verlet scaling factor, 0 by default
   tmode = 1 !default tangent mode for NEB calculation

   iemap=0     ! no emap constraint
   gammamap=1     ! default friction constant for map motion, 1/ps
   isgld = 0   ! no self-guiding
   isgsta=1    ! Begining index of SGLD range
   isgend=0    ! Ending index of SGLD range
   nsgsize=1       !   1--self atom only,
                   !   2--bond and angle atoms,
                   !   3--bond, angle, and dihedral atoms
   tsgavg=0.2d0    !  Local averaging time of SGLD simulation
   sgft=0.0d0      !  Guiding factor of SGLD simulation
   sgff=0.0d0      !  Guiding factor of SGLD simulation
   sgfg=0.0d0      !  Guiding factor of SGLD simulation
   tsgavp=2.0d0    !  Convergency time of SGLD simulation

   !     Check to see if "cntrl" namelist has been defined.
   mdin_cntrl=.false.
   mdin_qmmm = .false.
   mdin_ewald=.false.
   mdin_pb=.false.
   mdin_lmod=.false.
   mdin_amoeba=.false.
#ifdef MPI /* SOFT CORE */
   scalpha=0.5
   scbeta=12.0
   sceeorder=2
   ifsc=0
   logdvdl=0
   dvdl_norest=0
   dynlmb=0.0
   ifmbar=0
   bar_intervall=100
   bar_l_min=0.1
   bar_l_max=0.9
   bar_l_incr=0.1
   tishake = 0
#endif
   ilrt = 0
   lrt_interval = 50
   lrtmask=''
#ifdef DSSP
   idssp = 0
#endif

#ifdef API
   igb = input_options%igb
   alpb = input_options%alpb
   gbsa = input_options%gbsa
   lj1264 = input_options%lj1264
   ipb = input_options%ipb
   inp = input_options%inp
   vdwmeth = input_options%vdwmeth
   extdiel = input_options%extdiel
   intdiel = input_options%intdiel
   rgbmax = input_options%rgbmax
   saltcon = input_options%saltcon
   cut = input_options%cut
   dielc = input_options%dielc
   jfastw = input_options%jfastw
   ntf = input_options%ntf
   ntc = input_options%ntc
   fswitch = input_options%fswitch
   ntr = input_options%ntr
   ibelly = input_options%ibelly
   restraint_wt = input_options%restraint_wt
   restraintmask = input_options%restraintmask
   bellymask = input_options%bellymask
   refc = input_options%refc
   mask_from_ref = input_options%mask_from_ref ! Antoine Marion
#  ifdef LES
   rdt = input_options%rdt
#  endif
#else /* NOT the API */
   call nmlsrc('cntrl',5,ifind)
   if (ifind /= 0) mdin_cntrl=.true.
   rewind(5)

   call nmlsrc('ewald',5,ifind)
   if (ifind /= 0) mdin_ewald=.true.
   rewind 5

   call nmlsrc('pb',5,ifind)
   if (ifind /= 0) mdin_pb=.true.
   rewind 5

   call nmlsrc('qmmm', 5, ifind)
   if (ifind /= 0) mdin_qmmm = .true.
   rewind 5

   call nmlsrc('lmod',5,ifind)
   if (ifind /= 0) mdin_lmod=.true.
   rewind 5

   call nmlsrc('amoeba',5,ifind)
   if (ifind /= 0) mdin_amoeba=.true.
   rewind 5

   call nmlsrc('xray',5,ifind)
   xray_active = (ifind /= 0)
   rewind 5

   if ( mdin_cntrl ) then
      read(5,nml=cntrl,iostat=stat)
      if (stat /= 0) then
        backspace(5)
        read(5, fmt='(A)') line
        write(6, '(A)') 'Invalid line in &cntrl namelist '//': '//trim(line)
        call mexit(6,1)
      end if

   else
      write(6, '(1x,a,/)') 'Could not find cntrl namelist'
      FATAL_ERROR
   end if
#endif /* ifdef API */

   ! Set ntwr to nstlim if the user did not specify anything.
   ! Note nstlim can actually be set to zero but zero is not permitted for ntwr
   ! so if it is still zero set it to 1.
   if ( nstlim .eq. 0 ) then
      ntwr = 1
   else
      if ( ntwr .eq. 0 ) then
         if ( numexchg .eq. 0 ) then
            ntwr = nstlim
         else
            ntwr = nstlim * numexchg   ! for replica exchange
         end if
      end if
   end if

   if ( ntwr .eq. 0 ) ntwr = nstlim
   if ( ntwr .eq. 0 ) ntwr = 1

   if ( igb == 10 .and. ipb == 0 ) ipb = 2
   if ( igb == 0  .and. ipb /= 0 ) igb = 10

   ! middle scheme or SINR are now the only thermostats:
   if (ntt.ne.0 .and. ntt.ne.2 .and. ntt.ne.3 .and. ntt.ne.10) then
      write(6,'(1x,a,/)') &
         'Middle scheme/SINR: ntt is only available for 0,2,3 and 10'
      FATAL_ERROR
   end if
   if (gamma_ln < 0d0) then
      write(6,'(1x,a,/)') 'Middle scheme: gamma_ln MUST be non-negative'
      FATAL_ERROR
   end if
   if (sinrtau < 0.025d0) then
      write(6,'(1x,a,/)') 'sinrtau must be >= 0.025 when ntt = 10'
      FATAL_ERROR
   end if

   ! Now that we've read the input file, set up the defaults for variables
   ! whose values depend on other input values (ntb, cut)
   if (ntb == NO_INPUT_VALUE) then
      if (ntp > 0) then
         ntb = 2
      else if (igb > 0) then
         ntb = 0
      else
         ntb = 1
      end if
   end if

   if (ntxo == NO_INPUT_VALUE) then
      ntxo = 2
   end if

   if (cut == NO_INPUT_VALUE_FLOAT) then
      if (igb == 0) then
         cut = EIGHT
      else
         cut = 9999.d0
      end if
   end if

   if (nsnb == NO_INPUT_VALUE) then
      nsnb = 25   ! old default from sander
      if (imin>0) nsnb =1   ! safer if there are big initial forces
   endif

   ! Force igb=6 to get vacuum electrostatics or igb=0 for periodic
   ! boundary conditions. This must be done ASAP to ensure SANDER's
   ! electrostatics are initialized properly.

   rismprm%rism=irism

   if (irism /= 0) then
      periodicPotential = 'pme'

#   ifndef API
      write(6,'(a)') "|periodic 3D-RISM Forcing igb=0"
#   endif
      igb = 0
   end if

   if ( mdin_lmod ) then
      rewind 5
      call read_lmod_namelist()
   end if

   !--------------------------------------------------------------------
   !     --- vars have been read ---
   !--------------------------------------------------------------------

#ifndef API
   write(6,9309)
#endif

   ! emit warnings for retired cntrl namelist variables

   if ( dtemp /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/)') 'Warning: dtemp has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if
   if ( dxm /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/)') 'Warning: dxm has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
            ! '  The step length will be unlimited.'
   end if
   if ( heat /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/)') 'Warning: heat has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if
   if ( timlim /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/)') 'Warning: timlim has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if

   ! emit warnings for pmemd functionality that is not supported in sander:

   if ( scmask1 /= '' ) then
      write(6,'(/,a,/,a,/,a/)') 'Warning: scmask1 is a pmemd variable.', &
            '  It is not supported in sander.', &
            '  Check the Softcore Potentials section in the manual.'
   end if
   if ( scmask2 /= '' ) then
      write(6,'(/,a,/,a,/,a/)') 'Warning: scmask2 is a pmemd variable.', &
            '  It is not supported in sander.', &
            '  Check the Softcore Potentials section in the manual.'
   end if
   if ( timask1 /= '' ) then
      write(6,'(/,a,/,a,/,a/)') 'Warning: timask1 is a pmemd variable.', &
            '  It is not supported in sander.', &
            '  Check the Thermodynamic Integration section in the manual.'
   end if
   if ( timask2 /= '' ) then
      write(6,'(/,a,/,a,/,a/)') 'Warning: timask2 is a pmemd variable.', &
            '  It is not supported in sander.', &
            '  Check the Thermodynamic Integration section in the manual.'
   end if

! baroscalingdir valid options
   if (baroscalingdir > 0) then
      if (baroscalingdir < 0 .or. baroscalingdir > 3) then
         write(6,'(/2x,a)') &
         'Invalid baroscalingdir value. baroscalingdir must be between 0 and 3'
         FATAL_ERROR
      end if
      if (ntb /= 2) then
         write(6,'(/2x,a)') &
         'ntb invalid. ntb must be 2 for directional pressure scaling.'
         FATAL_ERROR
      end if
      if (ntp /= 2) then
         write(6,'(/2x,a)') &
         'ntp invalid. ntp must be 2 for directional pressure scaling.'
         FATAL_ERROR
      end if
   end if

! MC Barostat valid options. Some of these may work, but disable them until they
! are fully tested.

   if (ntp > 0 .and. barostat == 2) then
      inerr = 0
      if (icfe /= 0) then
         write(6, '(/2x,a)') 'TI is not compatible with the MC Barostat'
         inerr = 1
      end if
#ifdef LES
      write(6, '(/2x,a)') 'LES is not compatible with the MC Barostat'
      inerr = 1
#endif
      ! Any others? Hopefully most or all of the above can be made compatible.
      if (inerr == 1) then
         FATAL_ERROR
      end if
   end if

#ifndef API
   call printflags()
#endif

   !--------------------------------------------------------------------
   ! If user has requested ewald electrostatics, read some more input
   !--------------------------------------------------------------------

#ifdef API
   if( ipb == 0 ) call load_ewald_info(ntp)
#else
   if( ipb == 0 ) call load_ewald_info(inpcrd,ntp)
#endif

   temap=iemap>0
   ishake = 0
   if (ntc > 1) ishake = 1

   !--------------------------------------------------------------------
   ! Set up some parameters for GB simulations:
   !--------------------------------------------------------------------
   ! update offset = 0.09d0 for igb /= 8
   ! I add this step because I want to use different offset value as default value
   ! for igb = 8

   if ( igb == 8 ) then
     if (offset == -999999.d0) then
        offset = 0.195141d0  !set to default for igb=8
     end if
     if (gbneckscale == -999999.d0) then
        gbneckscale = 0.826836d0
     end if
   else
      ! not igb=8, use old defaults
      if (offset == -999999.d0) then
         offset = 0.09d0
      end if
      if (gbneckscale == -999999.d0) then
         gbneckscale = 0.361825d0
      end if
   endif

   if( igb == 2 .or. hybridgb == 2 ) then
      !       --- use our best guesses for Onufriev/Case GB  (GB^OBC I)

      gbgamma = 2.90912499999d0  ! (the "99999" to force roundoff on print)
      gbbeta = ZERO
      gbalpha = 0.8d0
   end if

   if( igb == 5 .or. hybridgb == 5 ) then

      !       --- use our second best guesses for Onufriev/Case GB (GB^OBC II)

      gbgamma = 4.850d0
      gbbeta = 0.8d0
      gbalpha = ONE
   end if

   if( igb == 7 ) then

      !       --- use parameters for Mongan et al. CFA GBNECK

      gbgamma = 2.50798245d0
      gbbeta = 1.90792938d0
      gbalpha = 1.09511284d0
   end if

   ! -------------------------------------------------------------------
   ! If the user has requested NMR restraints, do a cursory read of the
   ! restraints file(s) now to determine the amount of memory necessary
   ! for these restraints:
   ! -------------------------------------------------------------------

   if (jar == 1 ) nmropt = 1
   intreq = 0
   irlreq = 0
   if (nmropt > 0) then
      mxgrp = 0
      itotst = 1

      ! Set ITOTST to 0 if IMIN equals 1 (i.e. if minimization, not dynamics)
      ! This will cause any "time-averaged" requests to be over-ridden.

      if (imin == 1) then
         itotst = 0
      end if
      !         CALL AMOPEN(31,NMR,'O','F','R')
      call restlx(5,itotst,mxgrp,dt,6,ierr)
      !         CLOSE(31)
   end if

   ! Set the definition of the water molecule. The default definition is in
   ! WATDEF(4).

   read(watdef(1),'(A4)') iwtnm
   read(watdef(2),'(A4)') iowtnm
   read(watdef(3),'(A4)') ihwtnm(1)
   read(watdef(4),'(A4)') ihwtnm(2)
   if (watnam /= '    ') read(watnam,'(A4)') iwtnm
   if (owtnm /= '    ') read(owtnm, '(A4)') iowtnm
   if (hwtnm1 /= '    ') read(hwtnm1,'(A4)') ihwtnm(1)
   if (hwtnm2 /= '    ') read(hwtnm2,'(A4)') ihwtnm(2)

   return

#ifndef API

   ! --- input file polar opts read err trapping:

   9308 format(/10x,55('-'),/10x, &
         'MSANDER                                      2021', &
         /10x,55('-')/)
   9309 format(/80('-')/'   1.  RESOURCE   USE: ',/80('-')/)
   9700 format(/,'File Assignments:',/,15('|',a6,': ',a,/))
   9701 format('|',a6,': ',a)
#  ifdef MPI
   9702 format(7('|',a10,': ',a,/))
#  endif /* MPI */
#endif /* API */

#undef FATAL_ERROR
end subroutine mdread1 

!======================================================================
!          MDREAD2
!======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize to defaults and print the inputable variables.
subroutine mdread2(x,ix,ih)
#undef API
#ifdef API
#  define FATAL_ERROR ierr = 1; return
#  define DELAYED_ERROR ierr = 1; return
#else
#  define FATAL_ERROR call mexit(6, 1)
#  define DELAYED_ERROR inerr = 1
#endif /* API */

   use lmod_driver, only : LMOD_NTMIN_LMOD, LMOD_NTMIN_XMIN, write_lmod_namelist
   use findmask
#ifdef LES
   use genbornles, only: isnucat ! gbneck2nu: check if atom belongs to nuc or protein
#else
   use genborn, only: isnucat ! gbneck2nu: check if atom belongs to nuc or protein
#endif
   use constants, only : ZERO, ONE, TWO
   use parms, only: req
   use md_scheme, only: ntt, gamma_ln
   use nblist, only: a,b,c,alpha,beta,gamma,nbflag,skinnb,sphere,nbtell,cutoffnb
   use file_io_dat
   use sander_lib, only: upper
#ifdef LES
   use les_data, only : lestyp, lesfac, nlesty, temp0les
#endif
#ifndef API
   use emap, only : temap, emap_options
   use linear_response, only : lrt_interval
#else
#endif /* API */
#ifdef MPI
   use softcore, only : ifsc, scalpha, scbeta, dvdl_norest, &
                        sceeorder, logdvdl, dynlmb
   use mbar, only : ifmbar, bar_intervall, bar_l_min, bar_l_max, bar_l_incr
! REMD
   use remd, only : rem, rremd
   use sgld, only : isgld ! for RXSGLD
   use mpi
#endif /* MPI */
   use linear_response, only: ilrt, lrtmask
   use nfe_sander_proxy, only: infe
   use commandline_module, only: cpein_specified

   implicit none
   _REAL_ x(*)
#  include "../include/memory.h"
   integer ix(lasti)
   character(len=4) ih(*)
   integer ierr, nbond
   integer atom1,atom2
   logical belly,konst
   character(len=2) atype
   integer atomicnumber, hybridization
   integer ngrp,inerr,nr,ir,i,mxresat,j
   integer noshakegp( natom ), natnos
   integer ier
   logical errFlag
   _REAL_ emtmd, wallc
#ifndef LES
   logical newstyle
#endif /* LES */
#ifndef API
   integer ntmp
#endif /* API */

#ifdef MPI
   !     =========================== AMBER/MPI ===========================
#  ifdef MPI_DOUBLE_PRECISION
#     undef MPI_DOUBLE_PRECISION
#  endif
#  include "parallel.h"
   integer ist(MPI_STATUS_SIZE), partner, nbonh_c, num_noshake_c
   integer nquant_c, noshake_overlap_c
   integer crggp( natom )
   _REAL_ val
   _REAL_, allocatable :: repvals(:)
   !     ========================= END AMBER/MPI =========================
#endif /* MPI */
#  include "../include/md.h"
#  include "box.h"
#  include "nmr.h"
#  include "extra_pts.h"
#  include "ew_cntrl.h"
#  include "ew_pme_recip.h"
#  include "ew_erfc_spline.h"
#  include "ew_mpole.h"
#  include "ew_legal.h"
#  include "def_time.h"
#  include "tgtmd.h"
#  include "multitmd.h"
   character(11), parameter :: error_hdr =      '| ERROR:   '

   ! -------------------------------------------------------------------
   !      --- set up resat array, containing string identifying
   !          residue for each atom
   ! -------------------------------------------------------------------
   ierr = 0
   mxresat = min( natom, matom )
   ir = 0
   do i=1,mxresat
      if (i >= ix(ir+i02)) ir=ir+1
      write(resat(i),'(a4,1x,a4,i4)') ih(m04+i-1), &
            ih(ir+m02-1),ir
      !                                    ---null terminator:
      resat(i)(14:14) = char(0)
   end do
   close(unit=8)

   ! -------------------------------------------------------------------
   !     ----- SET THE DEFAULT VALUES FOR SOME VARIABLES -----
   ! -------------------------------------------------------------------

   nrp = natom

#ifndef API
   if (ifbox == 1) write(6, '(/5x,''BOX TYPE: RECTILINEAR'')')
   if (ifbox == 2) write(6, '(/5x,''BOX TYPE: TRUNCATED OCTAHEDRON'')')
   if (ifbox == 3) write(6, '(/5x,''BOX TYPE: GENERAL'')')
#endif

   if (ntr.eq.1 .or. ibelly.ne.0) nscm = 0

   nsolut =  nrp
   if ( nscm > 0 .and. ntb == 0 ) then
      ndfmin = 6   ! both translation and rotation com motion removed
      if (nsolut == 1) ndfmin = 3
      if (nsolut == 2) ndfmin = 5
   else if ( nscm > 0 ) then
      ndfmin = 3    ! just translation com will be removed
   else
      ndfmin = 0
   end if
   if(nscm <= 0) nscm = 0
   if (ntt == 3) ndfmin = 0 ! No COM motion removal for Langevin

   init = 3
   if (irest > 0) init = 4
   if (dielc <= ZERO ) dielc = ONE
   if (tautp <= ZERO ) tautp = 0.2d0

   !     ----- RESET THE CAP IF NEEDED -----

   ! ivcap == 0: Cap will be in effect if it is in the prmtop file (ifcap = 1)

   if(ivcap == 1) then
      ! Cap will be in effect even if not in prmtop file
      !   requires additional information in sander.in file as in the case of ivcap == 3, 4, or 5
      ifcap = 2
   else if(ivcap == 2) then
      ! Inactivate cap
      ifcap = 0
   else if(ivcap == 3) then
      ! Sphere -> not yet implemented
      ifcap = 3
   else if(ivcap == 4) then
      ! Orthorhombus
      ifcap = 4
   else if(ivcap == 5) then
      ! Shell of waters around solute
      ifcap = 5
   end if

   !Support for random seed using time of day in microsecs
   if( ig==-1 ) then
     !Turn on NO_NTT3_SYNC when ig=-1. This means the code no
     !longer synchronized the random numbers between streams when
     !running in parallel giving better scaling.
     no_ntt3_sync = 1
     !  old: call microsec(ig)
     call wallclock( wallc ) ! GNU fortran yields wallc to milliseconds
     wallc = modulo( 1.d3*wallc, 1.d6) ! should give six digits, positive
     ig = wallc
#ifdef MPI
     write (6, '(a,i8,a)') "| Note: ig = -1. Setting random seed to ", ig ," based on wallclock"
     write (6, '(a)') "|     and disabling the synchronization of random &
                               &numbers between tasks"
     write (6, '(a)') "|     to improve performance."
#else
#  ifndef API
     write (6, '(a,i8,a)') "| Note: ig = -1. Setting random seed to ", ig ," based on wallclock &
                                &time in microseconds."
#  endif /* API */
#endif
   end if

#ifdef MPI
   ! For some runs using multisander where the coordinates of the 
   ! two 'replicas' need to be identical, for example TI, it is 
   ! critical that the random number stream is synchronized between 
   ! all replicas.  Only use the IG value from worldrank=0. Ok to 
   ! broadcast between the various sander masters since
   ! they all call mdread2.
   if ( icfe > 0 ) then
      ! no_ntt3_sync currently does not work with softcore TI simulations
      ! see sc_lngdyn in softcore.F90
      if ( ifsc > 0 ) no_ntt3_sync = 0
      call mpi_bcast(ig, 1, MPI_INTEGER, 0, commmaster, ierr)
   end if
#endif

   ! -------------------------------------------------------------------
   !     ----- PRINT DATA CHARACTERIZING THE RUN -----
   ! -------------------------------------------------------------------

   nr = nrp
#ifndef API
   write(6,9328)
   write(6,9008) title
   write(6,'(/a)') 'General flags:'
   write(6,'(5x,2(a,i8))') 'imin    =',imin,', nmropt  =',nmropt
#endif

   ! Error Checking for REMD
#ifdef MPI
   if (rem/=0) then
      ! Make sure that the number of replicas is even so
      !  that they all have partners in the exchange step
      if (mod(numgroup,2).ne.0) then
         write (6,'(a)') "==================================="
         write (6,'(a)') "REMD requires an even # of replicas"
         write (6,'(a)') "==================================="
         FATAL_ERROR
      endif

      write(6,'(/a)') 'Replica exchange'
      write(6,'(5x,4(a,i8))') 'numexchg=',numexchg,', rem=',rem

      ! REPCRD option temporarily disabled
      if (repcrd == 0) write(6,'(a)') &
         "REMD WARNING: repcrd disabled. Only replica &
         &trajectories/output can be written."

      ! Check for correct number of exchanges
      if (numexchg <= 0) then
         write(6,'(a)') "REMD ERROR: numexchg must be > 0, "
         FATAL_ERROR
      endif

      ! RXSGLD
      if (isgld > 0 .and. rem < 0) then
         write(6, '(a)') 'Multi-D REMD and replica-exchange SGLD are not &
                         &supported yet!'
         FATAL_ERROR
      end if

      ! Hybrid GB
      if (numwatkeep >= 0) then
         write(6,'(5x,4(a,i8))') 'numwatkeep=',numwatkeep,', hybridgb=',hybridgb
         ! Check that user specified GB model for hybrid REMD
         if (hybridgb /= 2 .and. hybridgb /= 5 .and. hybridgb /= 1) then
            write(6,'(a)') "HYBRID REMD ERROR: hybridgb must be 1, 2, or 5."
            FATAL_ERROR
         endif
      else
         !Check that user did not specify GB model if no hybrid run.
         if (hybridgb /= 0) then
            write(6,'(a)') &
            "HYBRID REMD ERROR: numwatkeep must be >= 0 when hybridgb is set."
            FATAL_ERROR
         endif
      endif
      ! RREMD
      if (rremd>0) then
         write(6,'(5x,4(a,i8))') "rremd=",rremd, ", reservoir_exchange_step=", &
                 & reservoir_exchange_step
         if (mod(reservoir_exchange_step, 2) .ne. 0) then
            write(6,'(a)') "reservoir_exchange_step must be a multiple of 2"
            FATAL_ERROR
         end if
      endif

      ! M-REMD (rem < 0) requires netcdf output.
      if (rem < 0 .and. ioutfm .ne. 1) then
         write(6,'(a)') "ERROR: Multi-D REMD (rem < 0) requires NetCDF &
                        &trajectories (ioutfm=1)"
         FATAL_ERROR
      endif

#  ifdef LES
      ! DAN ROE: Temporarily disable LES REMD until it is verified with new
      !           REMD code
      if (rem==2) then
         write (6,*) "******* LES REM (rem==2) temporarily disabled. Stop. *******"
         FATAL_ERROR
      endif

      if (rem==2 .and. igb/=1) then
         write (6,*) ' partial REM (rem=2) only works with igb=1'
         FATAL_ERROR
      endif
#  else
      if (rem==2) then
         write(6,*) '*******  For rem ==2, partial REM'
         write(6,*) 'use sander.LES with topology created by addles'
         FATAL_ERROR
      endif
#  endif

   endif ! rem>0

#else /* _NO_ MPI */
   ! Check if user set numexchg with no MPI
   if (numexchg>0) write(6,'(a)') &
      "WARNING: numexchg > 0 - for REMD run please recompile sander for &
      &parallel runs."

   ! Check if user set numwatkeep or hybridgb with no MPI - not sensible.
   if (numwatkeep>=0) write(6,'(a)') &
      "WARNING: numwatkeep >= 0 - for hybrid REMD run please recompile &
      &sander for parallel runs."

   if (hybridgb>0) write(6,'(a)') &
      "WARNING: hybridgb > 0 - for hybrid REMD run please recompile &
      &sander for parallel runs."
#endif /* MPI */
   ! End error checking for REMD

#ifndef API
   write(6,'(/a)') 'Nature and format of input:'
   write(6,'(5x,4(a,i8))') 'ntx     =',ntx,', irest   =',irest, &
         ', ntrx    =',ntrx

   write(6,'(/a)') 'Nature and format of output:'
   write(6,'(5x,4(a,i8))') 'ntxo    =',ntxo,', ntpr    =',ntpr, &
         ', ntrx    =',ntrx,', ntwr    =',ntwr
   write(6,'(5x,5(a,i8))') 'ntwx    =',ntwx, &
         ', ntwv    =',ntwv,', ntwe    =',ntwe
   if (isynctraj > 0) &
     write(6, '(5x, a,i8)') 'isynctraj    =',isynctraj
   write(6,'(5x,2(a,i8),a,i7)') 'ioutfm  =',ioutfm, &
         ', ntwprt  =',ntwprt, &
         ', rbornstat=',rbornstat
   if (ntwf > 0) &
      write(6,'(5x, a,i8)') 'ntwf    =',ntwf
   write(6,'(/a)') 'Potential function:'
   write(6,'(5x,5(a,i8))') 'ntf     =',ntf,', ntb     =',ntb, &
         ', igb     =',igb,', nsnb    =',nsnb
   write(6,'(5x,3(a,i8))') 'ipol    =',ipol,', gbsa    =',gbsa, &
         ', iesp    =',iesp
   write(6,'(5x,3(a,f10.5))') 'dielc   =',dielc, &
         ', cut     =',cut,', intdiel =',intdiel
   if (lj1264 /= 0) &
      write(6, '(5x,a,i8)') 'lj1264  =',lj1264

   if ( igb /= 0 .and. igb /= 6 .and. igb /= 10 .and. ipb == 0 .and. igb /= 8) &
                                   then
      write(6,'(5x,3(a,f10.5))') 'saltcon =',saltcon, &
            ', offset  =',offset,', gbalpha= ',gbalpha
      write(6,'(5x,3(a,f10.5))') 'gbbeta  =',gbbeta, &
            ', gbgamma =',gbgamma,', surften =',surften
      write(6,'(5x,3(a,f10.5))') 'rdt     =',rdt, ', rgbmax  =',rgbmax, &
            '  extdiel =',extdiel
      write(6,'(5x,3(a,i8))') 'alpb  = ',alpb
   end if

   ! print output for igb=8
   if ( igb == 8 ) then
      write(6,'(5x,3(a,f10.5))') 'saltcon =',saltcon, &
            ', offset  =',offset,', surften =',surften
      write(6,'(5x,3(a,f10.5))') 'rdt     =',rdt, ', rgbmax  =',rgbmax, &
            '  extdiel =',extdiel
      write(6,'(5x,3(a,i8))') 'alpb  = ',alpb
      write(6,'(5x,3(a,f10.5))') 'gbalphaH  =',gbalphaH, &
            ', gbbetaH   =',gbbetaH,',  gbgammaH  = ',gbgammaH
      write(6,'(5x,3(a,f10.5))') 'gbalphaC  =',gbalphaC, &
            ', gbbetaC   =',gbbetaC,',  gbgammaC  = ',gbgammaC
      write(6,'(5x,3(a,f10.5))') 'gbalphaN  =',gbalphaN, &
            ', gbbetaN   =',gbbetaN,',  gbgammaN  = ',gbgammaN
      write(6,'(5x,3(a,f10.5))') 'gbalphaOS =',gbalphaOS, &
            ', gbbetaOS  =',gbbetaOS,',  gbgammaOS = ',gbgammaOS
      write(6,'(5x,3(a,f10.5))') 'gbalphaP  =',gbalphaP, &
            ', gbbetaP   =',gbbetaP,',  gbgammaP  = ',gbgammaP
      ! gbneck2nu pars
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_hnu  =', gb_alpha_hnu, &
        ', gb_beta_hnu   =', gb_beta_hnu, ',  gb_gamma_hnu  = ', gb_gamma_hnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_cnu  =', gb_alpha_cnu, &
        ', gb_beta_cnu   =', gb_beta_cnu, ',  gb_gamma_cnu  = ', gb_gamma_cnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_nnu  =', gb_alpha_nnu, &
        ', gb_beta_nnu   =', gb_beta_nnu, ',  gb_gamma_nnu  = ', gb_gamma_nnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_onu  =', gb_alpha_osnu, &
        ', gb_beta_onu   =', gb_beta_osnu, ',  gb_gamma_onu  = ', gb_gamma_osnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_pnu  =', gb_alpha_pnu, &
        ', gb_beta_pnu   =', gb_beta_pnu, ',  gb_gamma_pnu  = ', gb_gamma_pnu
   end if


   if( alpb /= 0 ) then
      write(6,'(5x,3(a,f10.5))') 'Arad =', Arad
   end if

   write(6,'(/a)') 'Frozen or restrained atoms:'
   write(6,'(5x,4(a,i8))') 'ibelly  =',ibelly,', ntr     =',ntr
   if( ntr == 1 ) write(6,'(5x,a,f10.5)') 'restraint_wt =', restraint_wt

   if( imin /= 0 ) then

      write(6,'(/a)') 'Energy minimization:'
      ! print inputable variables applicable to all minimization methods.
      write(6,'(5x,4(a,i8))') 'maxcyc  =',maxcyc,', ncyc    =',ncyc, &
            ', ntmin   =',ntmin
      write(6,'(5x,2(a,f10.5))') 'dx0     =',dx0, ', drms    =',drms

      ! Input flag ntmin determines the method of minimization
      select case ( ntmin )
      case ( LMOD_NTMIN_XMIN, LMOD_NTMIN_LMOD )
         call write_lmod_namelist( )
      case default
         ! invalid ntmin
         write(6,'(/2x,a,i3,a)') 'Error: Invalid NTMIN (',ntmin,').'
         stop
      end select

   else

      write(6,'(/a)') 'Molecular dynamics:'
      write(6,'(5x,4(a,i10))') 'nstlim  =',nstlim,', nscm    =',nscm, &
            ', nrespa  =',nrespa
      write(6,'(5x,3(a,f10.5))') 't       =',t, &
            ', dt      =',dt,', vlimit  =',vlimit

      if( ntt == 0 .and. tempi > 0.0d0 .and. irest == 0 ) then
         write(6,'(/a)') 'Initial temperature generation:'
         write(6,'(5x,a,i8)') 'ig      =',ig
         write(6,'(5x,a,f10.5)') 'tempi   =',tempi
      else if( ntt == 3) then
         write(6,'(/a)') 'Langevin dynamics temperature regulation:'
         write(6,'(5x,a,i8)') 'ig      =',ig
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, &
               ', tempi   =',tempi,', gamma_ln=', gamma_ln
      else if( ntt == 2 ) then
         write(6,'(/a)') 'Anderson (strong collision) temperature regulation:'
         write(6,'(5x,2(a,i8))') 'ig      =',ig, ', vrand   =',vrand
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, ', tempi   =',tempi
      end if

      if( ntp /= 0 ) then
         write(6,'(/a)') 'Pressure regulation:'
         write(6,'(5x,4(a,i8))') 'ntp     =',ntp
         write(6,'(5x,a,f10.5)') 'pres0   =',pres0
         if (barostat == 2) then
            write(6, '(5x,a)') 'Monte-Carlo Barostat:'
            write(6, '(5x,a,i8)') 'mcbarint  =', mcbarint
         end if
      end if

      if (baroscalingdir /= 0) then
         write(6,'(/a)') 'Directional pressure scaling:'
         write(6,'(5x,a,i8)') 'baroscalingdir  =', baroscalingdir
      end if

   end if

   if( ntc /= 1 ) then
      write(6,'(/a)') 'SHAKE:'
      write(6,'(5x,4(a,i8))') 'ntc     =',ntc,', jfastw  =',jfastw
      write(6,'(5x,3(a,f10.5))') 'tol     =',tol
   end if

   if( ifcap == 1 .or. ifcap == 2 .or. ifcap == 3 ) then
      write(6,'(/a)') 'Water cap:'
      write(6,'(5x,2(a,i8))') 'ivcap   =',ivcap,', natcap  =',natcap
      write(6,'(5x,2(a,f10.5))') 'fcap    =',fcap, ', cutcap  =',cutcap
      write(6,'(5x,3(a,f10.5))') 'xcap    =',xcap, ', ycap    =',ycap,    &
                                 ', zcap    =',zcap
   else if( ifcap == 4 ) then
      write(6,'(/a)') 'Orthorhombus:'
      write(6,'(5x,1(a,i8))')    'ivcap   =',ivcap
      write(6,'(5x,1(a,f10.5))') 'forth   =',forth
      write(6,'(5x,3(a,f10.5))') 'xlorth  =',xlorth,', ylorth  =',ylorth, &
                                 ', zlorth  =',zlorth
      write(6,'(5x,3(a,f10.5))') 'xorth   =',xorth, ', yorth   =',yorth,  &
                                 ', zorth   =',zorth
   else if( ifcap == 5 ) then
      write(6,'(/a)') 'Water shell:'
      write(6,'(5x,(a,i8,a,f10.5))') 'ivcap   =',ivcap,', cutcap  =',cutcap
   endif

   if( nmropt > 0 ) then
      write(6,'(/a)') 'NMR refinement options:'
      write(6,'(5x,4(a,i8))')'iscale  =',iscale,', noeskp  =',noeskp, &
            ', ipnlty  =',ipnlty,', mxsub   =',mxsub
      write(6,'(5x,3(a,f10.5))') 'scalm   =',scalm, &
            ', pencut  =',pencut,', tausw   =',tausw
   end if

   if( numextra > 0 ) then
      write(6,'(/a)') 'Extra-points options:'
      write(6,'(5x,4(a,i8))') 'frameon =',frameon, &
            ', chngmask=',chngmask
   end if

#  ifdef MPI /* SOFT CORE */
   if( icfe /= 0 .or. ifsc/=0) then
      write(6,'(/a)') 'Free energy options:'
      write(6,'(5x,3(a,i8))')       'icfe    =', icfe   , ', ifsc    =', ifsc,    ', klambda =', klambda
      write(6,'(5x,3(a,f8.4))')     'clambda =', clambda, ', scalpha =', scalpha, ', scbeta  =', scbeta
      write(6,'(5x,2(a,i8))')       'sceeorder =', sceeorder, ' dvdl_norest =', dvdl_norest
      write(6,'(5x,a,f8.4,a,i8)')   'dynlmb =', dynlmb, ' logdvdl =', logdvdl
   end if
   if ( ifmbar /= 0 ) then
      write (6,'(/a)') 'FEP MBAR options:'
      write(6,'(5x,a,i8,a,i8)') 'ifmbar  =', ifmbar, ',  bar_intervall = ', bar_intervall
      write(6,'(5x,3(a,f6.4))') 'bar_l_min =', bar_l_min, ',  bar_l_max =', bar_l_max, &
                                ',  bar_l_incr =', bar_l_incr
   end if
#  endif /* MPI */

   if (ilrt /= 0) then
      write (6,*)
      write (6,'(a,i4,a,i4)') ' Linear Response Theory: ilrt =', ilrt, ' lrt_interval =', lrt_interval
      write (6,*)
   end if

   if( itgtmd == 1 ) then
      write(6,'(/a)') 'Targeted molecular dynamics:'
      write(6,'(5x,3(a,f10.5))') 'tgtrmsd =',tgtrmsd, &
            ', tgtmdfrc=',tgtmdfrc
   end if

   if( ntb > 0 ) then
      write(6,'(/a)') 'Ewald parameters:'
      ew_type = 0   ! no longer an option
      write(6,'(5x,4(a,i8))') 'verbose =',verbose, &
            ', ew_type =',ew_type,', nbflag  =',nbflag, &
            ', use_pme =',use_pme
      write(6,'(5x,4(a,i8))') 'vdwmeth =',vdwmeth, &
            ', eedmeth =',eedmeth,', netfrc  =',netfrc
      write(6, 9002) a, b, c
      write(6, 9003) alpha, beta, gamma
      write(6, 9004) nfft1, nfft2, nfft3
      write(6, 9006) cutoffnb, dsum_tol
      write(6, 9007) ew_coeff
      write(6, 9005) order
      9002 format (5x,'Box X =',f9.3,3x,'Box Y =',f9.3,3x,'Box Z =',f9.3)
      9003 format (5x,'Alpha =',f9.3,3x,'Beta  =',f9.3,3x,'Gamma =',f9.3)
      9004 format (5x,'NFFT1 =',i5  ,7x,'NFFT2 =',i5  ,7x,'NFFT3 =',i5)
      9005 format (5x,'Interpolation order =',i5)
      9006 format (5x,'Cutoff=',f9.3,3x,'Tol   =',e10.3)
      9007 format (5x,'Ewald Coefficient =',f9.5)
   end if

#ifdef MPI
  ! Prints an error if REMD is to be performed with different cell sizes across
  ! the different replicas
  if (rem .ne. 0 .and. ntp .eq. 0 .and. ntb .ne. 0) then
    allocate(repvals(numgroup), stat=ierr)
    REQUIRE(ierr == 0)
    do i = 1, 6
       repvals(:) = 0.d0
       if (i .eq. 1) then
         val = a
       else if (i .eq. 2) then
         val = b
       else if (i .eq. 3) then
         val = c
       else if (i .eq. 4) then
         val = alpha
       else if (i .eq. 5) then
         val = beta
       else if (i .eq. 6) then
         val = gamma
       end if
       call mpi_allgather(val, 1, mpi_double_precision, &
                          repvals, 1, mpi_double_precision, &
                          commmaster, ierr)
       do j = 2, numgroup
         if (abs(repvals(j)-repvals(1)) .gt. 1.0d-5) then
           write(6,'(a,a,a)') error_hdr, "The box sizes do not match for all REMD replicas."
           write(6,'(11x,a)') "This can happen, for example, if you are using different input coordinate files"
           write(6,'(11x,a)') "for the different replicas and these structures come from a NPT run."
           write(6,'(11x,a)') "Make sure the box lengths and angles are the same for all replicas in all input"
           write(6,'(11x,a)') "coordinate files."
           FATAL_ERROR
         end if
       end do
    end do
  end if
#endif

! ----EMAP Options-----
   if(temap)call emap_options(5)
! ---------------------

#  ifdef MPI
! --- MPI TIMING OPTIONS ---
      write(6, '(/a)') '| MPI Timing options:'
      write(6, '("|",5x," profile_mpi = ",i8)') profile_mpi
! Sanity check for profile_mpi
      call int_legal_range('profile_mpi',profile_mpi,0,1)
! --------------------------
#  endif /* MPI */
#endif /* API */

   cut = cut*cut
   cut_inner = cut_inner*cut_inner

   !------------------------------------------------------------------------
   ! If user has requested generalized born electrostatics, set up variables
   !------------------------------------------------------------------------

   if( igb == 0 .and. gbsa > 0 ) then
      write(0,*) 'GB/SA calculation is performed only when igb>0'
      FATAL_ERROR
   end if
   if( gbsa == 2 .and. &
       ((imin == 0 .and. nstlim > 1) .or. &
        (imin == 1 .and. maxcyc > 1)) ) then
      write(0,*) 'GBSA=2 only works for single point energy calc'
      FATAL_ERROR
   end if
   if( gbsa == 3 ) then
      write(0,*) 'GBSA=3 only works for pmemd, not sander'
      FATAL_ERROR
   end if
   if ( igb /= 0 .and. igb /= 10 .and. ipb == 0 ) then
#if defined (LES) && !defined (API)
      write(6,*) 'igb=1,5,7 are working with LES, no SA term included'
#endif
      ! igb7 uses special S_x screening params.
      ! overwrite the tinker values read from the prmtop
      if (igb == 7) then
         do i=1,natom
            if(ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if

            if (atomicnumber .eq. 6) then
               x(l96+i-1) = 4.84353823306d-1
            else if (atomicnumber .eq. 1) then
               x(l96+i-1) = 1.09085413633d0
            else if (atomicnumber .eq. 7) then
               x(l96+i-1) = 7.00147318409d-1
            else if (atomicnumber .eq. 8) then
               x(l96+i-1) = 1.06557401132d0
            else if (atomicnumber .eq. 16) then
               x(l96+i-1) = 6.02256336067d-1
            else if (atomicnumber .eq. 15) then
               x(l96+i-1) = 5d-1
            else
               x(l96+i-1) = 5d-1
            end if
         end do
      end if

      ! changing S_x screening params for igb = 8
      ! overwrite the tinker values read from the prmtop

      if (igb == 8) then
         do i=1,natom
            if (ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if

            call isnucat(nucat,i,nres,60,ix(i02:i02+nres-1),ih(m02:m02+nres-1))
            ! check if atom belong to nucleic or protein residue
            !60 = the total number of different nucleic acid residue type in AMBER
            !Need to update this number if adding new residue type
            !update isnucat function in gb_ene.F90 too
            if (nucat == 1) then
                ! If atom belongs to nuc, use gbneck2nu pars
                if (atomicnumber .eq. 6) then
                  x(l96+i-1) = screen_cnu
                else if (atomicnumber .eq. 1) then
                  x(l96+i-1) = screen_hnu
                else if (atomicnumber .eq. 7) then
                  x(l96+i-1) = screen_nnu
                else if (atomicnumber .eq. 8) then
                  x(l96+i-1) = screen_onu
                else if (atomicnumber .eq. 15) then
                  x(l96+i-1) = screen_pnu
                else
                  x(l96+i-1) = 5.d-1 ! not optimized
                end if
            else
               ! protein
               if (atomicnumber .eq. 6) then
                  x(l96+i-1) = Sc
               else if (atomicnumber .eq. 1) then
                  x(l96+i-1) = Sh
               else if (atomicnumber .eq. 7) then
                  x(l96+i-1) = Sn
               else if (atomicnumber .eq. 8) then
                  x(l96+i-1) = So
               else if (atomicnumber .eq. 16) then
                  x(l96+i-1) = Ss
               else if (atomicnumber .eq. 15) then
                  x(l96+i-1) = Sp ! We still don't have an optimized Sp parameter
                                  ! if nuc, use gbneck2nu pars
               else
                  !for atom type Cl,Br,...
                  !These parameters are also not optimized.
                  x(l96+i-1) = 5d-1
               end if
            end if !check isnucat
         end do !loop i=1, natom
      end if ! igb = 8

      ! Set up for igb == 2, 5, 7, 8
      ! Put gb parameters in arrays
      if ( igb == 2 .or. igb == 5 .or. igb == 7 .or. &
           hybridgb == 2 .or. hybridgb == 5) then
        do i=1,natom
            x(l2402+i-1) = gbalpha
            x(l2403+i-1) = gbbeta
            x(l2404+i-1) = gbgamma
        end do
      end if

      ! IGB = 8, update gb_alpha, gb_beta, gb_gamma
       if ( igb == 8 ) then
           do i=1,natom
              if (ix(i100) .eq. 1) then
                 atomicnumber = ix(i100+i)
              else
                 call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                    atomicnumber, errFlag)
#ifdef LES
                 if(errFlag) then
                    call get_atomic_number(ih(m04+i-1), &
                       x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                       atomicnumber, errFlag)
                 end if
#endif
              end if

          ! check if atom belongs to protein or nucleic acid
          call isnucat(nucat,i,nres,60,ix(i02:i02+nres-1),ih(m02:m02+nres-1))
          if (nucat == 1) then
             !if atom belong to nucleic part, use nuc pars
             if (atomicnumber .eq. 1) then
               x(l2402+i-1) = gb_alpha_hnu
               x(l2403+i-1) = gb_beta_hnu
               x(l2404+i-1) = gb_gamma_hnu
             else if (atomicnumber .eq. 6) then
               x(l2402+i-1) = gb_alpha_cnu
               x(l2403+i-1) = gb_beta_cnu
               x(l2404+i-1) = gb_gamma_cnu
             else if (atomicnumber .eq. 7) then
               x(l2402+i-1) = gb_alpha_nnu
               x(l2403+i-1) = gb_beta_nnu
               x(l2404+i-1) = gb_gamma_nnu
             else if (atomicnumber .eq. 8) then
               x(l2402+i-1) = gb_alpha_osnu
               x(l2403+i-1) = gb_beta_osnu
               x(l2404+i-1) = gb_gamma_osnu
             else if (atomicnumber .eq. 16) then
               x(l2402+i-1) = gb_alpha_osnu
               x(l2403+i-1) = gb_beta_osnu
               x(l2404+i-1) = gb_gamma_osnu
             else if (atomicnumber .eq. 15) then
               x(l2402+i-1) = gb_alpha_pnu
               x(l2403+i-1) = gb_beta_pnu
               x(l2404+i-1) = gb_gamma_pnu
             else
               ! Use GB^OBC (II) (igb = 5) set for other atoms
               x(l2402+i-1) = 1.0d0
               x(l2403+i-1) = 0.8d0
               x(l2404+i-1) = 4.851d0
             end if ! assign nucleic acid pars

          else
             !not nucleic part, use protein pars
             if (atomicnumber .eq. 1) then
                     x(l2402+i-1) = gbalphaH
                     x(l2403+i-1) = gbbetaH
                     x(l2404+i-1) = gbgammaH
             else if (atomicnumber .eq. 6) then
                     x(l2402+i-1) = gbalphaC
                     x(l2403+i-1) = gbbetaC
                     x(l2404+i-1) = gbgammaC
             else if (atomicnumber .eq. 7) then
                     x(l2402+i-1) = gbalphaN
                     x(l2403+i-1) = gbbetaN
                     x(l2404+i-1) = gbgammaN
             else if (atomicnumber .eq. 8) then
                     x(l2402+i-1) = gbalphaOS
                     x(l2403+i-1) = gbbetaOS
                     x(l2404+i-1) = gbgammaOS
             else if (atomicnumber .eq. 16) then
                     x(l2402+i-1) = gbalphaOS
                     x(l2403+i-1) = gbbetaOS
                     x(l2404+i-1) = gbgammaOS
             else if (atomicnumber .eq. 15) then
                     x(l2402+i-1) = gbalphaP
                     x(l2403+i-1) = gbbetaP
                     x(l2404+i-1) = gbgammaP
             else
                     !use GBOBC set for other atom types
                     x(l2402+i-1) = 1.0d0
                     x(l2403+i-1) = 0.8d0
                     x(l2404+i-1) = 4.851d0
             end if ! assign protein pars
          end if !check nucat
        end do !loop i=1,natom
       end if ! end IGB = 8 section

     !       put fs(i)*(rborn(i) - offset) into the "fs" array
      fsmax = 0.d0
      do i=1,natom
         x(l96-1+i) = x(l96-1+i)*( x(l97-1+i) - offset )
         fsmax = max( fsmax, x(l96-1+i) )
         if (rbornstat == 1) then
            x(l186-1+i) = 0.d0
            x(l187-1+i) = 999.d0
            x(l188-1+i) = 0.d0
            x(l189-1+i) = 0.d0
         end if
      end do

      !     ---------------------------------------------------------------------
      !       ---get Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
      !         T = 298.15, epsext=78.5,

      kappa = sqrt( 0.10806d0 * saltcon )

      !       ---scale kappa by 0.73 to account(?) for lack of ion exclusions:

      kappa = 0.73d0* kappa

      !   ---------------------------------------------------------------------

      if ( gbsa == 1 ) then

         !     --- assign parameters for calculating SASA according to the
         !         LCPO method ---

         do i=1,natom
            ix(i80+i-1)=0
         end do

         !         --- get the number of bonded neighbors for each atom:

         do i=1,nbona
            atom1=ix(iiba+i-1)/3+1
            atom2=ix(ijba+i-1)/3+1
            ix(i80+atom1-1)=ix(i80+atom1-1)+1
            ix(i80+atom2-1)=ix(i80+atom2-1)+1
         end do

         do i=1,natom
            ix(i80-i)=ix(i80+i-1)
         end do

         do i=1,nbonh
            atom1=ix(iibh+i-1)/3+1
            atom2=ix(ijbh+i-1)/3+1
            ix(i80-atom1)=ix(i80-atom1)+1
            ix(i80-atom2)=ix(i80-atom2)+1
         end do


         !         --- construct parameters for SA calculation; note that the
         !             radii stored in L165 are augmented by 1.4 Ang.

         do i=1,natom
            write(atype,'(a2)') ih(m06+i-1)
            call upper(atype)
            if (ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if
            hybridization = ix(i80-i)
            nbond=ix(i80+i-1)
            if (atomicnumber .eq. 6) then
               if (hybridization .eq. 4) then
                  if (nbond == 1) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.77887d0
                     x(l175-1+i) = -0.28063d0
                     x(l180-1+i) = -0.0012968d0
                     x(l185-1+i) = 0.00039328d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.56482d0
                     x(l175-1+i) = -0.19608d0
                     x(l180-1+i) = -0.0010219d0
                     x(l185-1+i) = 0.0002658d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.23348d0
                     x(l175-1+i) = -0.072627d0
                     x(l180-1+i) = -0.00020079d0
                     x(l185-1+i) = 0.00007967d0
                  else if (nbond == 4) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.00000d0
                     x(l175-1+i) = 0.00000d0
                     x(l180-1+i) = 0.00000d0
                     x(l185-1+i) = 0.00000d0
                  else
                     write(6,*) 'Unusual nbond for CT:', i, nbond, &
                        ' Using default carbon LCPO parameters'
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.77887d0
                     x(l175-1+i) = -0.28063d0
                     x(l180-1+i) = -0.0012968d0
                     x(l185-1+i) = 0.00039328d0
                  end if
               else
                  if (nbond == 2) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.51245d0
                     x(l175-1+i) = -0.15966d0
                     x(l180-1+i) = -0.00019781d0
                     x(l185-1+i) = 0.00016392d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.070344d0
                     x(l175-1+i) = -0.019015d0
                     x(l180-1+i) = -0.000022009d0
                     x(l185-1+i) = 0.000016875d0
                  else
                     write(6,*) 'Unusual nbond for C :', i, nbond, &
                       ' Using default carbon LCPO parameters'
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.77887d0
                     x(l175-1+i) = -0.28063d0
                     x(l180-1+i) = -0.0012968d0
                     x(l185-1+i) = 0.00039328d0
                  end if
               end if
            else if (atomicnumber .eq. 8) then
               if (atype == 'O ') then
                  x(l165-1+i) = 1.60d0 + 1.4d0
                  x(l170-1+i) = 0.68563d0
                  x(l175-1+i) = -0.1868d0
                  x(l180-1+i) = -0.00135573d0
                  x(l185-1+i) = 0.00023743d0
               else if (atype == 'O2') then
                  x(l165-1+i) = 1.60d0 + 1.4d0
                  x(l170-1+i) = 0.88857d0
                  x(l175-1+i) = -0.33421d0
                  x(l180-1+i) = -0.0018683d0
                  x(l185-1+i) = 0.00049372d0
               else
                  if (nbond == 1) then
                     x(l165-1+i) = 1.60d0 + 1.4d0
                     x(l170-1+i) = 0.77914d0
                     x(l175-1+i) = -0.25262d0
                     x(l180-1+i) = -0.0016056d0
                     x(l185-1+i) = 0.00035071d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.60d0 + 1.4d0
                     x(l170-1+i) = 0.49392d0
                     x(l175-1+i) = -0.16038d0
                     x(l180-1+i) = -0.00015512d0
                     x(l185-1+i) = 0.00016453d0
                  else
                     write(6,*) 'Unusual nbond for O:', i, nbond, &
                        ' Using default oxygen LCPO parameters'
                     x(l165-1+i) = 1.60d0 + 1.4d0
                     x(l170-1+i) = 0.77914d0
                     x(l175-1+i) = -0.25262d0
                     x(l180-1+i) = -0.0016056d0
                     x(l185-1+i) = 0.00035071d0
                  end if
               end if
            else if(atomicnumber .eq. 7) then
               if (atype == 'N3') then
                  if (nbond == 1) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.078602d0
                     x(l175-1+i) = -0.29198d0
                     x(l180-1+i) = -0.0006537d0
                     x(l185-1+i) = 0.00036247d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.22599d0
                     x(l175-1+i) = -0.036648d0
                     x(l180-1+i) = -0.0012297d0
                     x(l185-1+i) = 0.000080038d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.051481d0
                     x(l175-1+i) = -0.012603d0
                     x(l180-1+i) = -0.00032006d0
                     x(l185-1+i) = 0.000024774d0
                  else
                     write(6,*) 'Unusual nbond for N3:', i, nbond, &
                        ' Using default nitrogen LCPO parameters'
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.078602d0
                     x(l175-1+i) = -0.29198d0
                     x(l180-1+i) = -0.0006537d0
                     x(l185-1+i) = 0.00036247d0
                  end if
               else
                  if (nbond == 1) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.73511d0
                     x(l175-1+i) = -0.22116d0
                     x(l180-1+i) = -0.00089148d0
                     x(l185-1+i) = 0.0002523d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.41102d0
                     x(l175-1+i) = -0.12254d0
                     x(l180-1+i) = -0.000075448d0
                     x(l185-1+i) = 0.00011804d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.062577d0
                     x(l175-1+i) = -0.017874d0
                     x(l180-1+i) = -0.00008312d0
                     x(l185-1+i) = 0.000019849d0
                  else
                     write(6,*) 'Unusual nbond for N:', i, nbond, &
                        ' Using default nitrogen LCPO parameters'
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.078602d0
                     x(l175-1+i) = -0.29198d0
                     x(l180-1+i) = -0.0006537d0
                     x(l185-1+i) = 0.00036247d0
                  end if
               end if
            else if(atomicnumber .eq. 16) then
               if (atype == 'SH') then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.7722d0
                  x(l175-1+i) = -0.26393d0
                  x(l180-1+i) = 0.0010629d0
                  x(l185-1+i) = 0.0002179d0
               else
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.54581d0
                  x(l175-1+i) = -0.19477d0
                  x(l180-1+i) = -0.0012873d0
                  x(l185-1+i) = 0.00029247d0
               end if
            else if (atomicnumber .eq. 15) then
               if (nbond == 3) then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.3865d0
                  x(l175-1+i) = -0.18249d0
                  x(l180-1+i) = -0.0036598d0
                  x(l185-1+i) = 0.0004264d0
               else if (nbond == 4) then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.03873d0
                  x(l175-1+i) = -0.0089339d0
                  x(l180-1+i) = 0.0000083582d0
                  x(l185-1+i) = 0.0000030381d0
               else
                  write(6,*) 'Unusual nbond for P:', i, nbond, &
                     ' Using default phosphorus LCPO parameters'
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.3865d0
                  x(l175-1+i) = -0.18249d0
                  x(l180-1+i) = -0.0036598d0
                  x(l185-1+i) = 0.0004264d0
               end if
            else if (atype(1:1) == 'Z') then
               x(l165-1+i) = 0.00000d0 + 1.4d0
               x(l170-1+i) = 0.00000d0
               x(l175-1+i) = 0.00000d0
               x(l180-1+i) = 0.00000d0
               x(l185-1+i) = 0.00000d0
            else if (atomicnumber .eq. 1) then
               x(l165-1+i) = 0.00000d0 + 1.4d0
               x(l170-1+i) = 0.00000d0
               x(l175-1+i) = 0.00000d0
               x(l180-1+i) = 0.00000d0
               x(l185-1+i) = 0.00000d0
            else if (atype == 'MG') then
               !  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.45A: Aqvist 1992
               x(l165-1+i) = 1.18d0 + 1.4d0
               !  The following values were taken from O.sp3 with two bonded
               !  neighbors -> O has the smallest van der Waals radius
               ! compared to all other elements which had been parametrized
               x(l170-1+i) = 0.49392d0
               x(l175-1+i) = -0.16038d0
               x(l180-1+i) = -0.00015512d0
               x(l185-1+i) = 0.00016453d0
            else if (atype == 'F') then
               x(l165-1+i) = 1.47d0 + 1.4d0
               x(l170-1+i) = 0.68563d0
               x(l175-1+i) = -0.1868d0
               x(l180-1+i) = -0.00135573d0
               x(l185-1+i) = 0.00023743d0
            else
               ! write( 0,* ) 'bad atom type: ',atype
               ! call mexit( 6,1 )
               x(l165-1+i) = 1.70 + 1.4;
               x(l170-1+i) = 0.51245;
               x(l175-1+i) = -0.15966;
               x(l180-1+i) = -0.00019781;
               x(l185-1+i) = 0.00016392;
               write(6,'(a,a)') 'Using carbon SA parms for atom type', atype
            end if
         end do  !  i=1,natom
         !
      else if ( gbsa == 2 ) then

         !     --- assign parameters for calculating SASA according to the
         !         ICOSA method; the radii are augmented by 1.4 A ---

         do i=1,natom
            write(atype,'(a2)') ih(m06+i-1)
            if(ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if
            if (atomicnumber .eq. 7) then
               x(L165-1+i) = 1.55d0 + 1.4d0
            else if (atomicnumber .eq. 6) then
               x(L165-1+i) = 1.70d0 + 1.4d0
            else if (atomicnumber .eq. 1) then
               x(L165-1+i) = 1.20d0 + 1.4d0
            else if (atomicnumber .eq. 8) then
               x(L165-1+i) = 1.50d0 + 1.4d0
            else if (atomicnumber .eq. 15) then
               x(L165-1+i) = 1.80d0 + 1.4d0
            else if (atomicnumber .eq. 16) then
               x(L165-1+i) = 1.80d0 + 1.4d0
            else if (atomicnumber .eq. 17) then
               ! Cl radius
               x(L165-1+i) = 1.70d0 + 1.4d0
            else if (atomicnumber .eq. 9) then
               ! F radius
               x(L165-1+i) = 1.50d0 + 1.4d0
            else if (atomicnumber .eq. 35) then
               ! Br radius
               ! Bondi, J. Phys. Chem. 1964, 68, 441.
               x(L165-1+i) = 1.85d0 + 1.4d0
            else if (atomicnumber .eq. 20) then
               ! Ca radius
               ! Calculated from Aqvist, J. Phys. Chem. 1990, 94, 8021.
               x(L165-1+i) = 1.33d0 + 1.4d0
            else if (atomicnumber .eq. 11) then
               ! Na radius
               ! Calculated from Aqvist, J. Phys. Chem. 1990, 94, 8021.
               x(L165-1+i) = 1.87d0 + 1.4d0
            else if (atomicnumber .eq. 30) then
               ! Zn radius
               ! Hoops, Anderson, Merz, JACS 1991, 113, 8262.
               x(L165-1+i) = 1.10d0 + 1.4d0
            else if (atomicnumber .eq. 12) then
               !  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.45A: Aqvist 1992
               x(L165-1+i) = 1.18d0 + 1.4d0
            else
               write( 0,* ) 'bad atom type: ',atype
               FATAL_ERROR
            end if

            ! dummy LCPO values:
            x(L170-1+i) = 0.0d0
            x(L175-1+i) = 0.0d0
            x(L180-1+i) = 0.0d0
            x(L185-1+i) = 0.0d0
            !  write(6,*) i,' ',atype,x(L165-1+i)
         end do  !  i=1,natom

      end if ! ( gbsa == 1 )

   end if  ! ( igb /= 0 .and. igb /= 10 .and. ipb == 0 )

   !-----------------------------------
   ! If a LRT calculation is requested,
   ! setup the icosa-SASA parameters
   ! this code is copied from above
   !-----------------------------------

   if ( ilrt /= 0 ) then
      !     --- assign parameters for calculating SASA according to the
      !         ICOSA method; the radii are augmented by 1.4 A ---
      do i=1,natom
         write(atype,'(a2)') ih(m06+i-1)
         if(ix(i100) .eq. 1) then
            atomicnumber = ix(i100+i)
         else
            call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
               atomicnumber, errFlag)
#ifdef LES
            if(errFlag) then
               call get_atomic_number(ih(m04+i-1), &
                  x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                  atomicnumber, errFlag)
            end if
#endif
         end if
         if (atomicnumber .eq. 7) then
            x(L165-1+i) = 1.55d0 + 1.4d0
         else if (atomicnumber .eq. 6) then
            x(L165-1+i) = 1.70d0 + 1.4d0
         else if (atomicnumber .eq. 1 .or. &
              ! added for lone pairs
              atype == 'EP') then
            x(L165-1+i) = 1.20d0 + 1.4d0
         else if (atomicnumber .eq. 8) then
            x(L165-1+i) = 1.50d0 + 1.4d0
         else if (atomicnumber .eq. 15) then
            x(L165-1+i) = 1.80d0 + 1.4d0
         else if (atomicnumber .eq. 16) then
            x(L165-1+i) = 1.80d0 + 1.4d0
         else if (atomicnumber .eq. 12) then
            !             Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
            !             Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
            !             Mg radius = 1.45A: Aqvist 1992
            x(L165-1+i) = 1.18d0 + 1.4d0
         else
            write( 0,* ) 'bad atom type: ',atype,' cannot perform SASA calculation'
            FATAL_ERROR
         end if  !  atype(1:1) == 'N'
         x(L170-1+i) = 0.0d0
         x(L175-1+i) = 0.0d0
         x(L180-1+i) = 0.0d0
         x(L185-1+i) = 0.0d0
         !write(6,*) i,' ',atype,x(L165-1+i)
      end do  !  i=1,natom

   end if ! ( ilrt /= 0 )

   ! -------------------------------------------------------------------
   !
   ! -- add check to see if the space in nmr.h is likely to be
   !     too small for this run:
   !     [Note: this check does *not* indicate if MXTAU and MXP are
   !      too small.  There is no easy way to ensure this, since
   !      the experimental intensities are read in a namelist
   !      command: if too many intensities are input, the read
   !      statment may cause a coredump before returning control
   !      to the main program.  Be careful.  sigh....]

   if (natom > matom .and. nmropt > 1) then
      write(6,*) 'WARNING: MATOM in nmr.h is smaller than the ', &
            natom,' atoms in this molecule.'
      write(6,*) 'Printout of NMR violations may be compromised.'
   end if

   ! -------------------------------------------------------------------
   !     --- checks on bogus data ---
   ! -------------------------------------------------------------------

   inerr = 0

   if( icfe < 0 .or. icfe > 1 ) then
      write(6,*) 'icfe must be 0 or 1 (icfe=2 is no longer supported)'
      DELAYED_ERROR
   end if
   if( icfe /= 0 .and. numgroup /= 2 ) then
      write(6,*) 'numgroup must be 2 if icfe is set'
      DELAYED_ERROR
   end if
   if (ievb>0) then
#ifdef MPI
  if( numgroup /= 2 ) then
     write(6,*) 'numgroup must be 2 if ievb is set'
     DELAYED_ERROR
  end if
#else
      write(6,'(/2x,a)') 'Setting ievb>0 requires compilation with MPI'
      DELAYED_ERROR
#endif
   end if
   if( (igb > 0 .or. ipb /= 0) .and. numextra > 0) then
      if (igb /= 6) then
         write(6,'(a)') 'Cannot use igb>0 (except igb=6) with extra-point force fields'
         DELAYED_ERROR
      end if
   end if

   if (igb /= 0 .and. igb /= 1 .and. igb /= 2 .and. igb /= 5 &
         .and. igb /= 6 .and. igb /= 7 .and. igb /= 8 .and. igb /= 10) then
      write(6,'(/2x,a,i3,a)') 'IGB (',igb,') must be 0,1,2,5,6,7,8 or 10.'
      DELAYED_ERROR
   end if
   if (alpb /= 0 .and. alpb /= 1 )  then
      write(6,'(/2x,a,i3,a)') 'ALPB (',alpb,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (alpb /= 0 .and. igb /= 1 .and. igb /= 2 .and. igb /= 5 .and. igb /=7 )  then
      write(6,'(/2x,a,i3,a)') 'IGB (',igb,') must be 1,2,5, or 7 if ALPB > 0.'
      DELAYED_ERROR
   end if

   if (jar < 0 .or. jar > 1) then
      write(6,'(/2x,a,i3,a)') 'JAR (',jar,') must be 0 or 1'
      DELAYED_ERROR
   end if

#ifdef LES
   if( igb /= 0 .and. igb /= 1 .and. igb /= 5 .and. igb /=7 ) then
      write(6,'(/,a)') 'Error: LES is only compatible with IGB > 0,1,5,7'
      DELAYED_ERROR
   end if
   if( alpb /= 0) then
      write(6,'(/,a)') 'Error: LES is not compatible with ALPB'
      DELAYED_ERROR
   end if
   if( gbsa > 0 ) then
      write(6,'(/,a)') 'Error: LES is not compatible with GBSA > 0'
      DELAYED_ERROR
   end if
   if( ipol > 0 ) then
      write(6,'(/,a)') 'Error: LES is not compatible with IPOL > 0'
      DELAYED_ERROR
   end if
   if (temp0les >= 0.d0 .and. iscale > 0 ) then
      write (6,'(/,a)') 'Error: iscale cannot be used with temp0les'
      DELAYED_ERROR
   end if
#else
   rdt = 0
#endif
   if (irest /= 0 .and. irest /= 1) then
      write(6,'(/2x,a,i3,a)') 'IREST (',irest,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (ibelly /= 0 .and. ibelly /= 1) then
      write(6,'(/2x,a,i3,a)') 'IBELLY (',ibelly,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (imin < 0) then
      write(6,'(/2x,a,i3,a)') 'IMIN (',imin,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntxo .eq. 0) then
      write(6, '(/2x,a,a)') 'Old style binary restart files (ntxo=0) are no longer supported.'
      DELAYED_ERROR
   end if
   if (ntxo /= 1 .and. ntxo /= 2) then
      write(6, '(/2x,a,i3,a)') 'NTXO (',ntxo,') must be 1 or 2.'
      DELAYED_ERROR
   end if

   if (imin == 5) then
      if (ifbox /= 0 .and. ntb == 2) then
         write(6,'(/2x,a)') 'WARNING: IMIN=5 with changing periodic boundaries (NTB=2) can result in'
         write(6,'(/2x,a)') '         odd energies being calculated. Use with caution.'
      endif
#ifdef MPI
      if (sandersize > 1 .and. ntb == 2) then
         write(6,'(/2x,a)') 'ERROR: IMIN=5 and NTB=2 cannot be run with multiple processors.'
         DELAYED_ERROR
      endif
#endif
   end if


   if (fswitch > 0) then
      if (fswitch >= cut) then
         write(6, '(/2x,a)') 'fswitch must be less than the cutoff'
         DELAYED_ERROR
      end if
      if (igb /= 0 .or. ntb == 0) then
         write(6, '(/2x,a)') 'fswitch is only compatible with periodic simulations'
         DELAYED_ERROR
      end if
      if (lj1264 > 0) then
         write(6, '(/2x,a)') 'fswitch is incompatible with the 12-6-4 L-J model'
         DELAYED_ERROR
      end if
   end if

   if (iscale > mxvar) then
      write(6,9501) iscale,mxvar
      9501 format('ERROR: ISCALE (',i5,') exceeds MXVAR (',i5, &
            '). See nmr.h')
      DELAYED_ERROR
   end if
   if (ntx < 1 .or. ntx > 7) then
      write(6,'(/2x,a,i3,a)') 'NTX (',ntx,') must be in 1..7'
      DELAYED_ERROR
   end if

   if (ntb /= 0 .and. ntb /= 1 .and. ntb /= 2) then
      write(6,'(/2x,a,i3,a)') 'NTB (',ntb,') must be 0, 1 or 2.'
      DELAYED_ERROR
   end if

   if (infe /= 0 .and. infe /= 1) then
      write(6,'(/2x,a,i3,a)') 'IFE (',infe,') must be 0 or 1.'
      DELAYED_ERROR
   end if

   if (ntp /= 0 .and. ntp /= 1 .and. ntp /= 2 .and. ntp /= 3) then
      write(6,'(/2x,a,i3,a)') 'NTP (',ntp,') must be 0, 1, 2, or 3.'
      DELAYED_ERROR
   end if
   if (npscal < 0 .or. npscal > 1) then
      write(6,'(/2x,a,i3,a)') 'NPSCAL (',npscal,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (ntp > 0) then
      if (barostat /= 2) then
         write(6, '(/2x,a,i3,a)') 'BAROSTAT (', barostat, ') must be 2'
         DELAYED_ERROR
      end if
      if (barostat == 2) then
         if (mcbarint <= 0) then
            write(6, '(/2x,a,i3,a)') 'MCBARINT (',mcbarint,') must be positive'
            DELAYED_ERROR
         end if
         if (mcbarint >= nstlim) then
            write(6, '(a)') 'WARNING: mcbarint is greater than the number of &
                     &steps. This is effectively constant volume.'
         end if
      end if
   end if

   if (ntc < 1 .or. ntc > 4) then
      write(6,'(/2x,a,i3,a)') 'NTC (',ntc,') must be 1,2,3 or 4.'
      DELAYED_ERROR
   end if
   if (jfastw < 0 .or. jfastw > 4) then
      write(6,'(/2x,a,i3,a)') 'JFASTW (',jfastw,') must be 0->4.'
      DELAYED_ERROR
   end if

   if (ntf < 1 .or. ntf > 8) then
      write(6,'(/2x,a,i3,a)') 'NTF (',ntf,') must be in 1..8.'
      DELAYED_ERROR
   end if

   if (ioutfm /= 0 .and. ioutfm /= 1) then
      write(6,'(/2x,a,i3,a)') 'IOUTFM (',ioutfm,') must be 0 or 1.'
      DELAYED_ERROR
   end if

   if (ntpr < 0) then
      write(6,'(/2x,a,i3,a)') 'NTPR (',ntpr,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntwx < 0) then
      write(6,'(/2x,a,i3,a)') 'NTWX (',ntwx,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntwv < -1) then
      write(6,'(/2x,a,i3,a)') 'NTWV (',ntwv,') must be >= -1.'
      DELAYED_ERROR
   end if
   if (ntwf < -1) then
      write(6, '(/2x,a,i3,a)') 'NTWF (',ntwf,') must be >= -1.'
      DELAYED_ERROR
   end if
   if (ntwv == -1 .and. ioutfm /= 1) then
      write (6, '(/2x,a)') 'IOUTFM must be 1 for NTWV == -1.'
      DELAYED_ERROR
   end if
   if (ntwf == -1 .and. ioutfm /= 1) then
      write(6, '(/2x,a)') 'IOUTFM must be 1 for NTWF == -1.'
      DELAYED_ERROR
   end if
   if (ntwv == -1 .and. ntwx == 0) then
      write (6, '(/2x,a)') 'NTWX must be > 0 for NTWV == -1.'
      DELAYED_ERROR
   end if
   if (ntwe < 0) then
      write(6,'(/2x,a,i3,a)') 'NTWE (',ntwe,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntave < 0) then
      write(6,'(/2x,a,i3,a)') 'NTAVE (',ntave,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntr /= 0 .and. ntr /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTR (',ntr,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (ntrx /= 0 .and. ntrx /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTRX (',ntrx,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (nmropt < 0 .or. nmropt > 2) then
      write(6,'(/2x,a,i3,a)') 'NMROPT (',nmropt,') must be in 0..2.'
      DELAYED_ERROR
   end if

   ! check settings related to ivcap

   if(ivcap == 3 .or. ivcap == 4) then
      write(6,'(/2x,a)') 'IVCAP == 3 and IVCAP == 4 currently not implemented'
      DELAYED_ERROR
   endif
   if (ivcap < 0 .and. ivcap > 5) then
      write(6,'(/2x,a)') 'IVCAP must be 0 ... 5'
      DELAYED_ERROR
   end if
   if ((ivcap == 1 .or. ivcap == 5) .and. igb /= 10 .and. ipb == 0) then
      write(6,'(/2x,a)') 'IVCAP == 1,5 only works with Poisson Boltzmann (igb=10 or ipb/=0)'
      DELAYED_ERROR
   end if
   if((ivcap == 1 .or. ivcap == 3 .or. ivcap == 5 ) .and. cutcap <= 0.0d0) then
      write(6,'(/2x,a)') 'For IVCAP == 1,3, or 5, cutcap must be > 0.0'
      DELAYED_ERROR
   endif
   if (ivcap == 4 .and. &
         (xlorth < ZERO .or. ylorth < ZERO .or. zlorth < ZERO .or. &
!  give magic numbers a name  srb aug 2007 !
          xorth > 47114710.0d0 .or. &
!  give magic numbers a name !
          yorth > 47114710.0d0 .or. &
!  give magic numbers a name !
          zorth > 47114710.0d0)) then
      write(6,'(/2x,a)') &
      'For IVCAP == 4, xlorth, ylorth, zlorth, xorth, yorth, zorth must be set'
      DELAYED_ERROR
   end if
   if ((ivcap == 3 .or. ivcap == 4) .and. ibelly == 0) then
      write(6,'(/2x,a,a)') &
         'For IVCAP == 3 or 4, ibelly must be 1 and all atoms', &
         '  not in the spherical or orthorhombic region must be set NOT moving'
      DELAYED_ERROR
   end if
   if (ivcap == 5 .and. (imin /= 1 .or. maxcyc > 1)) then
      write(6,'(/2x,a,a)') &
         'IVCAP == 5 only works for single-point energy calculation'
      DELAYED_ERROR
   end if

   ! check if ifbox variable from prmtop file matches actual angles:

   if ( ipb == 0 .and. ntb /= 0 ) then
      if ( ifbox == 1 ) then
         if ( abs(alpha - 90.0d0) > 1.d-5 .or. &
           abs(beta  - 90.0d0) > 1.d-5 .or. &
           abs(gamma - 90.0d0) > 1.d-5 ) then
           ifbox =3
#ifndef API
           write(6,'(a)') '     Setting ifbox to 3 for non-orthogonal unit cell'
#endif
         end if
      end if

      if ( ifbox == 2 ) then
         if ( abs(alpha - 109.4712190d0) > 1.d-5 .or. &
              abs(beta  - 109.4712190d0) > 1.d-5 .or. &
              abs(gamma - 109.4712190d0) > 1.d-5 ) then
              write(6,'(/2x,a)') &
              'Error: ifbox=2 in prmtop but angles are not correct'
              DELAYED_ERROR
         end if
      end if
   end if

   ! checks for targeted MD
   if (itgtmd /= 0 .and. itgtmd > 2) then
      write(6,'(/2x,a,i3,a)') 'ITGTMD (',itgtmd,') must be 0, 1 or 2.'
      DELAYED_ERROR
   end if
   if (itgtmd == 1 .and. ntr == 1) then
      if (len_trim(tgtfitmask) > 0 .or. len_trim(tgtrmsmask) <= 0) then
         write(6,'(/2x,a)') 'ITGTMD: tgtrmsmask (and not tgtfitmask) ' //  &
                            'should be specified if NTR=1'
         DELAYED_ERROR
      end if
   end if

   !     -- consistency checking

   if (imin > 0.and.nrespa > 1)  then
      write(6,'(/2x,a)') 'For minimization, set nrespa,nrespai=1'
      DELAYED_ERROR
   end if
   if (ntp > 0 .and. nrespa > 1) then
      write(6,'(/2x,a)') 'nrespa must be 1 if ntp>0'
      DELAYED_ERROR
   end if
   if  (ntx < 4.and.init /= 3)  then
      write(6,'(/2x,a)') 'NTX / IREST inconsistency'
      DELAYED_ERROR
   end if
   if (ntb == 2 .and. ntp == 0) then
      write(6,'(/2x,a)') 'NTB set but no NTP option (must be 1 or 2)'
      DELAYED_ERROR
   end if
   if (ntp /= 0 .and. ntb /= 2) then
      write(6,'(/,a,a)')' NTP > 0 but not constant pressure P.B.C.', &
            ' (NTB = 2) must be used'
      DELAYED_ERROR
   end if
   if (ntb /= 0 .and. ifbox == 0 .and. ntp /= 0) then
      write(6,'(/,a)') ' (NTB /= 0 and NTP /= 0) but IFBOX == 0'
      write(6,'(/,a)') ' This combination is not supported'
      DELAYED_ERROR
   end if
   if (ntb /= 0 .and. &
         ( box(1) < 1.d0  .or. &
         box(2) < 1.d0  .or. &
         box(3) < 1.d0 ) ) then
      write(6,'(/,a,3f10.3)') ' BOX is too small: ',box(1),box(2),box(3)
      DELAYED_ERROR
   else if (ntb /= 0 .and. &
         (sqrt(cut) >= box(1)*0.5d0 .or. &
         sqrt(cut) >= box(2)*0.5d0 .or. &
         sqrt(cut) >= box(3)*0.5d0) ) then
      write(6,'(/,a)') ' CUT must be < half smallest box dimension'
      DELAYED_ERROR
   end if
#if 0
   if (ntb /= 0 .and. (igb > 0 .or. ipb /= 0) ) then
      write(6,'(/,a)') ' igb>0 is only compatible with ntb=0'
      DELAYED_ERROR
   end if
#endif
   if ( ntb == 0 .and. sqrt(cut) < 8.05 .and. igb /= 10 .and. ipb == 0 .and. &
        igb /= 6 ) then
      write(6,'(/,a,f8.2)') ' unreasonably small cut for non-periodic run: ', &
         sqrt(cut)
      DELAYED_ERROR
   end if
   if ( rgbmax < 5.d0*fsmax ) then
      write(6,'(/,a,f8.2)') ' rgbmax must be at least ', 5.d0*fsmax
      DELAYED_ERROR
   end if
   if (icfe /= 0 .and. ibelly /= 0 ) then
      write(6,'(/,a)') ' ibelly cannot be used with icfe'
      DELAYED_ERROR
   end if
#ifdef MPI /* SOFT CORE */
   if (icfe /= 0 .and. dvdl_norest /= 0 ) then
      write(6,'(/,a)') 'dvdl_norest must == 0!'
      write(6,'(/,a)') 'The dvdl_norest option is deprecated.'
      write(6,'(/,a)') 'Restraint energies are seperated, &
        &and do not contribute to dvdl.'
      DELAYED_ERROR
   end if
#endif
   ! Modification done by Ilyas Yildirim
   if (icfe == 1 .and. (klambda < 1 .or. klambda > 6)) then
     write(6,'(/,a)') ' klambda must be between 1 and 6'
     DELAYED_ERROR
   end if
   ! End of modification done by Ilyas Yildirim

   if (clambda < 0.d0 .or. clambda > 1.d0 ) then
      write(6,'(/,a)') ' clambda must be between 0 and 1'
      DELAYED_ERROR
   end if

#ifdef MPI /* SOFT CORE */
   if (ifsc /= 0) then
      if (ifsc == 2) then
         write (6,'(/,a)') 'The ifsc=2 option is no longer supported. &
               &Internal energies of the soft core region.'
         write (6,'(a)') 'are now handled implicitly and setting ifsc=2 &
               &is no longer needed'
         DELAYED_ERROR
      end if
      if (icfe /= 1 .and. ifsc==1) then
         write (6,'(/,a)') ' Softcore potential requires a standard TI run, set icfe to 1'
         DELAYED_ERROR
      end if
      if ( igb > 0 .or. ipb /= 0 ) then
         write (6,'(/,a)') ' Softcore potential is incompatible with GB (for now)'
         DELAYED_ERROR
      end if
      if ( ntf > 1 ) then
         write (6,'(/,a)') ' Softcore potentials require ntf=1 because SHAKE &
               &constraints on some bonds might be removed'
         DELAYED_ERROR
      end if
      if (clambda > 0.995 .or. clambda < 0.005) then
         write (6,'(/,a)') ' Softcore potentials cannot be used with clambda < 0.005 or > 0.995'
         DELAYED_ERROR
      end if
      if (klambda /= 1) then
         write (6,'(/,a)') ' Softcore potential requires linear mixing, set klambda to 1'
         DELAYED_ERROR
      end if
   end if
   ! Following is temporary (dec 2014), until/unless problem gets fixed
   if ( logdvdl /= 0 .and. ifsc == 0 ) then
      write (6,'(/,a)') ' Cannot request logdvdl unless ifsc>0'
      DELAYED_ERROR
   end if
   if (ifmbar /= 0) then
      if (icfe /= 1) then
         write (6,'(/,a)') ' MBAR requires a standard TI run, set icfe to 1'
         DELAYED_ERROR
      end if
      if (ifsc /=0 .and. (bar_l_max > 0.995 .or. bar_l_min < 0.005) ) then
         write (6,'(/,a)') ' Softcore potentials cannot be used with &
               &bar_l_min < 0.005 or bar_l_max > 0.995'
         DELAYED_ERROR
      end if
      if (klambda /= 1) then
         write (6,'(/,a)') ' MBAR requires linear mixing, set klambda to 1'
         DELAYED_ERROR
      end if
   end if
#endif

   if (ilrt /= 0 .and. lrtmask == '') then
      write (6,'(a)') 'Linear Response Theory activated, but lrtmask is not set'
      DELAYED_ERROR
   end if

   !-----------------------------------------------------
   !     ----sanity checks for Ewald
   !-----------------------------------------------------

   if( igb == 0 .and. ipb == 0 ) then
      call float_legal_range('skinnb: (nonbond list skin) ', &
            skinnb,skinlo,skinhi)

      !  --- Will check on sanity of settings after coords are read in
      !      and the extent of the system is determined.

      if(periodic == 1)then
         call float_legal_range('skinnb+cutoffnb: (nonbond list cut) ', &
               skinnb+cutoffnb,zero,sphere)
      end if
      if (ntb==0 .and. use_pme/=0) then
         write(6,'(/,a)') &
         'Using PME with a non-periodic simulation does not make sense. Set either ntb>0 of use_pme=0.'
         DELAYED_ERROR
      end if
      call float_legal_range('a: (unit cell size) ',a,boxlo,boxhi)
      call float_legal_range('b: (unit cell size) ',b,boxlo,boxhi)
      call float_legal_range('c: (unit cell size) ',c,boxlo,boxhi)
      call float_legal_range('alpha: (unit cell angle) ', &
            alpha,anglo,anghi)
      call float_legal_range('beta: (unit cell angle) ', &
            beta,anglo,anghi)
      call float_legal_range('gamma: (unit cell angle) ', &
            gamma,anglo,anghi)
      call int_legal_range('order: (interpolation order) ', &
            order,orderlo,orderhi)
      call opt_legal_range('verbose: ',verbose,0,4)
      call opt_legal_range('netfrc: ',netfrc,0,1)
      call opt_legal_range('nbflag: ',nbflag,0,1)
      call opt_legal_range('nbtell: ',nbtell,0,2)
      call opt_legal_range('vdwmeth: ',vdwmeth,0,2)
      call opt_legal_range('eedmeth: ',eedmeth,1,6)
      call opt_legal_range('ee_type: ',ee_type,1,2)
      call opt_legal_range('maxiter: ',maxiter,1,50)
      call opt_legal_range('indmeth: ',indmeth,0,3)
      call opt_legal_range('fix_quad: ',fix_quad,0,1)
      call float_legal_range('eedtbdns: (erfc table density) ', &
            eedtbdns,denslo,denshi)
   end if  ! ( igb == 0 .and. ipb == 0 )

   ! ---WARNINGS:

   if ( ibelly == 1 .and. igb == 0 .and. ipb == 0 .and. ntb /= 0 ) then
#ifndef API
      write(6,'(/,a,/,a,/,a)') &
            'Warning: Although EWALD will work with belly', &
            '(for equilibration), it is not strictly correct!'
#endif
   end if

   if (inerr == 1) then
      write(6, '(/,a)') ' *** input error(s)'
      FATAL_ERROR
   end if

   ! Load the restrained atoms (ntr=1) or the belly atoms (ibelly=1)
   ! or atoms for targeted md (itgtmd=1). Selections are read from
   ! &cntrl variables or, if these are not defined, it falls back to
   ! the old group input format.

   if(mtmd /= 'mtmd') then
     itgtmd = 2
     ntr = 0
     emtmd = 0.0d0
   end if
   konst = ntr > 0
   dotgtmd = itgtmd > 0
   belly = .false.
   natc = 0
   ngrp = 0
   natbel = 0
   nattgtfit = 0  ! number of atoms for tgtmd fitting (=overlap)
   nattgtrms = 0  ! number of atoms for tgtmd rmsd calculation
   nrc = 0
   if(konst.or.dotgtmd) then

      ! inserted here to fix the bug that coords are not available
      ! yet when distance based selection (<,>) is requested
#ifdef LES
#  ifdef MPI
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0les,.FALSE.,solvph,solve)
#  else
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,.FALSE.)
#  endif
#else
#  ifndef API
#  ifdef MPI
         call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0, &
            .FALSE.,solvph,solve,0)
#  else
         call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,.FALSE.)
#  endif
#  endif /* API */
#endif /* LES */

      if(itgtmd == 2) then
        call mtmdcall(emtmd,x(lmtmd01),ix(imtmd02),x(lcrd),x(lforce),ih(m04),ih(m02),ix(i02),&
                    ih(m06),x(lmass),natom,nres,'READ')
      else
#ifndef API
        if (konst) write(6,9408)
        if (dotgtmd) write(6,9409)
#endif
        call rdrest(natom,ntrx,refc,x(lcrdr))
        !close(10)

        ! VH - tgtmd change: preferably call atommask() instead of rgroup()
        if (konst) then
          if( len_trim(restraintmask) <= 0 ) then
              call rgroup(natom,natc,nres,ngrp,ix(i02),ih(m02),ih(m04), &
                          ih(m06),ih(m08),ix(icnstrgp), &
                          npdec,x(l60),konst,dotgtmd,belly,5,.true.)
          else
              if (mask_from_ref > 0) then ! Antoine Marion : base mask on reference coordinates
                call atommask( natom, nres, 0, ih(m04), ih(m06), &
                  ix(i02), ih(m02), x(lcrdr), restraintmask, ix(icnstrgp) )
              else
                call atommask( natom, nres, 0, ih(m04), ih(m06), &
                  ix(i02), ih(m02), x(lcrd), restraintmask, ix(icnstrgp) )
              endif

              ! for now, emulate the "GATHER ALL THE CONSTRAINED ATOMS TOGETHER"
              ! section of rgroup(); later, the various masks should be done
              ! differently, i.e. without the "gather", as in the following:
              !     x(l60:l60+natom-1) = restraint_wt
              !     natc = sum(ix(icnstrgp:icnstrgp+natom-1))

              natc = 0
              do i=1,natom
                if( ix(icnstrgp-1+i) <= 0 ) cycle
                natc = natc + 1
                ix(icnstrgp-1+natc) = i
                x(l60-1+natc) = restraint_wt
              end do
#ifndef API
              write(6,'(a,a,a,i7,a)') '     Mask ', &
              restraintmask(1:len_trim(restraintmask)), ' matches ',natc,' atoms'
#endif
          end if
        end if
        nrc = natc

        if (itgtmd == 1) then
          if (len_trim(tgtfitmask) <= 0 .and. len_trim(tgtrmsmask) <= 0) then
              ! the following if-endif can be deleted when we stop
              ! supporting rgroup()
              if (konst) then
                ! cannot do both ntr and tgtmd together using old group format
                write(6,'(/2x,a)') 'NTR must be 0 for targeted MD (TGTMD=1)'
                call mexit(6,1)
              else  ! the following only for backward compatibility
                call rgroup(natom,natc,nres,ngrp,ix(i02),ih(m02),ih(m04), &
                            ih(m06),ih(m08),ix(icnstrgp), &
                            npdec,x(l60),konst,dotgtmd,belly,5,.true.)
                ! tgtmd atoms are now stored in nattgt, igroup -> icnstrgp
                nattgtfit = natc
                nattgtrms = natc
                do i=1,nattgtfit
                    ix(itgtfitgp-1+i) = ix(icnstrgp-1+i)
                    ix(itgtrmsgp-1+i) = ix(icnstrgp-1+i)
                end do
              end if
          else
              if (ntr == 0) then  ! read tgtfitmask only if ntr=1
                ! read in atom group for tgtmd fitting (=overlap region)
                call atommask( natom, nres, 0, ih(m04), ih(m06), &
                    ix(i02), ih(m02), x(lcrd), tgtfitmask, ix(itgtfitgp) )
                ! see comments above (for ntr) for the following reduction cycle
                nattgtfit = 0
                do i=1,natom
                  if( ix(itgtfitgp-1+i) <= 0 ) cycle
                  nattgtfit = nattgtfit + 1
                  ix(itgtfitgp-1+nattgtfit) = i
                end do
                write(6,'(a,a,a,i7,a)')  &
                '     Mask "', tgtfitmask(1:len_trim(tgtfitmask)-1),  &
                '" matches ',nattgtfit,' atoms'
              end if
              ! read in atom group for tgtmd rmsd calculation
              call atommask( natom, nres, 0, ih(m04), ih(m06), &
                ix(i02), ih(m02), x(lcrd), tgtrmsmask, ix(itgtrmsgp) )
              nattgtrms = 0
              do i=1,natom
                if( ix(itgtrmsgp-1+i) <= 0 ) cycle
                nattgtrms = nattgtrms + 1
                ix(itgtrmsgp-1+nattgtrms) = i
              end do
              write(6,'(a,a,a,i7,a)')  &
              '     Mask "', tgtrmsmask(1:len_trim(tgtrmsmask)-1),  &
              '" matches ',nattgtrms,' atoms'
          end if
        end if

      end if ! (itgtmd == 2)

   end if  ! (konst.or.dotgtmd)

   ! dotgtmd may be false here even if doing tgtmd
   ! this is so belly info is read properly? following existing KONST code

   dotgtmd=.false.
   konst = .false.
   belly = ibelly > 0
   ngrp = 0
   if(belly) then
      ! inserted here to fix the bug that coords are not available
      ! yet when distance based selection (<,>) is requested
#ifdef LES
#  ifdef MPI
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0les,.FALSE.,solvph,solve)
#  else
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,.FALSE.)
#  endif
#else
#  ifdef MPI
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0,.FALSE.,solvph,solve, 0)
#  else
#    ifndef API
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,.FALSE.)
#    endif
#  endif
#endif /* LES */
#ifndef API
      write(6,9418)
#endif /* API */
      if( len_trim(bellymask) <= 0 ) then
         call rgroup(natom,natbel,nres,ngrp,ix(i02),ih(m02),ih(m04),ih(m06), &
                     ih(m08),ix(ibellygp),npdec, &
                     x(l60),konst,dotgtmd,belly,5,.true.)
      else
         if (mask_from_ref > 0) then ! Antoine Marion : base mask on reference coordinates
           call rdrest(natom,ntrx,refc,x(lcrdr))
           call atommask( natom, nres, 0, ih(m04), ih(m06), &
              ix(i02), ih(m02), x(lcrdr), bellymask, ix(ibellygp) )
         else
           call atommask( natom, nres, 0, ih(m04), ih(m06), &
              ix(i02), ih(m02), x(lcrd), bellymask, ix(ibellygp) )
         endif
         natbel = sum(ix(ibellygp:ibellygp+natom-1))
#ifndef API
         write(6,'(a,a,a,i7,a)') '     Mask ', &
            bellymask(1:len_trim(bellymask)), ' matches ',natbel,' atoms'
#endif
      end if
   end if
   call setvar(ix,belly)

   !  see if the user has input a noshakemask string, and process it:
   natnos = 0
   if( len_trim(noshakemask) > 0 ) then
      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), noshakemask, noshakegp )
      natnos = sum(noshakegp(1:natom))
      write(6,*)
      write(6,'(a,a,a,i7,a)') 'Noshake mask ', &
         noshakemask(1:len_trim(noshakemask)), ' matches ',natnos,' atoms'
      call setnoshake(ix,noshakegp,ntc,num_noshake)
      if( ntf > 1 ) then
         write(6,'(a)') '   Setting ntf to 1'
         ntf = 1
      end if
   end if

#ifdef MPI /* SOFT CORE */
   ! lower charges if a crgmask is set
   if ( len_trim(crgmask) > 0 ) then
      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), crgmask, crggp )
      write(6,'(a,a,a,i7,a)') 'Zero-Charge Mask ',crgmask(1:len_trim(crgmask)), ' matches ',sum(crggp(1:natom)),' atoms'
      call remove_charges(crggp, natom, x(l15))
   end if
#endif

   konst = .false.
   belly = .false.

   if( ibelly > 0 .and. (igb > 0 .or. ipb /= 0) ) then

      !          ---here, the only allowable belly has just the first
      !             NATBEL atoms in the moving part.  Check to see that this
      !             requirement is satisfied:

      do i=natbel+1,natom
         if( ix(ibellygp+i-1) /= 0 ) then
            write(6,*) 'When igb>0, the moving part must be at the'
            write(6,*) '   start of the molecule.  This does not seem'
            write(6,*) '   to be the case here.'
            write(6,*) 'natbel,i,igroup(i) = ' &
                  ,natbel,i,ix(ibellygp+i-1)
            call mexit(6,1)
         end if
      end do
   end if


   !     ----- CALCULATE THE SQUARE OF THE BOND PARAMETERS FOR SHAKE
   !           THE PARAMETERS ARE PUT SEQUENTIALLY IN THE ARRAY CONP -----

   do i=1,nbonh + nbona + nbper
      j = ix(iicbh+i-1)
      x(l50+i-1) = req(j)**2
   end do

#ifdef MPI
      if( icfe /= 0 ) then

      !  use the masses of the prmtop file for the first group for both groups:
      !  [only the master nodes communicate here, since non-master nodes
      !   have not yet allocated space]
      ! This leads to problems for dual topology runs, and is therefore skipped
      ! if ifsc is set to one, the masses from both prmtop files are used
         if (ifsc == 0) then
            call mpi_bcast(x(lmass),natom,MPI_DOUBLE_PRECISION,0,commmaster,ierr)
            call mpi_bcast(x(lwinv),natom,MPI_DOUBLE_PRECISION,0,commmaster,ierr)
         end if
         tmass = sum(x(lmass:lmass+natom-1))
         tmassinv = 1.d0/tmass

      !  next, do a minimal sanity check that the SHAKE parameters are
      !  consistent on the two processors:

         ! For Softcore this might be allowed
         ! Put a better check here later
         if( ntc == 2 .and. ifsc == 0) then
            partner = ieor(masterrank,1)
            call mpi_sendrecv( nbonh, 1, MPI_INTEGER, partner, 5, &
                               nbonh_c, 1, MPI_INTEGER, partner, 5, &
                               commmaster, ist, ierr )
            call mpi_sendrecv( num_noshake, 1, MPI_INTEGER, partner, 5, &
                               num_noshake_c, 1, MPI_INTEGER, partner, 5, &
                               commmaster, ist, ierr )
            if( nbonh - num_noshake /= nbonh_c - num_noshake_c ) then
               write(6,*) 'SHAKE lists are not compatible in the two groups!'
               call mexit(6,1)
            end if
         else if( ntc == 3 ) then
            write(6,*) 'ntc = 3 is not compatible with icfe>0'
            call mexit(6,1)
         end if

      end if
#endif

      if( igb == 0 .and. ipb == 0 ) &
         call init_extra_pts( &
         ix(iibh),ix(ijbh),ix(iicbh), &
         ix(iiba),ix(ijba),ix(iicba), &
         ix(i24),ix(i26),ix(i28),ix(i30), &
         ix(i32),ix(i34),ix(i36),ix(i38), &
         ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
         ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
         ih(m06),ix,x,ix(i08),ix(i10), &
         nspm,ix(i70),tmass,tmassinv,x(lmass),x(lwinv),req)

   !  DEBUG input; force checking
#ifdef API
   call load_debug
#else
   call load_debug(5)
#endif

   return
   ! -------------------------------------------------------------------------
   ! Standard format statements:

#ifndef API
   9328 format(/80('-')/,'   2.  CONTROL  DATA  FOR  THE  RUN',/80('-')/)
   9408 format(/4x,'LOADING THE CONSTRAINED ATOMS AS GROUPS',/)
   9409 format(/4x,'LOADING THE TARGETED MD ATOMS AS GROUPS',/)
   9418 format(/4x,'LOADING THE BELLY ATOMS AS GROUPS',/)
   9428 format(/4x,'LOADING THE DECOMP ATOMS AS GROUPS',/)
   9008 format(a80)
#endif

#undef FATAL_ERROR
#undef DELAYED_ERROR
end subroutine mdread2 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit defined preprocessor names, ie, flags.
#include "mdreadutils.h"
