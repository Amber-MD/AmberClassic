#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Open input files and read cntrl namelist.
subroutine mdread1()
   
   use memory_module
   implicit none

#  include "pb_constants.h"
#  include "files.h"
#  include "box.h"
#  include "md.h"
#  include "parms.h"

   integer     ifind
   logical     mdin_cntrl  ! true if this namelist exists in mdin
   character(len=8) date
   character(len=10) time
   _REAL_ rnum
   
   namelist /cntrl/ ntx,    ipb,    inp,   igb,            &
                    imin,   ntf,    ntb,   dielc,  cut,    & !compatibility
                    nsnb,   maxcyc, ntmin,                 & !compatibility
                    ivcap,  cutcap, xcap,  ycap,   zcap,   & !compatibility
                    idecomp, ig

   if (mdout /= "stdout" ) call myopen(6,mdout,owrite,'F','W')

   call myopen(5,mdin,'O','F','R')

   write(6,9308)

   call date_and_time( DATE=date, TIME=time )

   write(6,'(12(a))') '| Run on ', date(5:6), '/', date(7:8), '/',  &
        date(1:4), ' at ', time(1:2), ':', time(3:4), ':', time(5:6)

   if (owrite /= 'N') write(6, '(2x,a)') '[-O]verwriting output'
   
   ! Echo the file assignments to the user:
   
   write(6,9700) 'MDIN'   ,mdin(1:70)  , 'MDOUT' ,mdout(1:70) , &
                 'INPCRD' ,inpcrd(1:70), 'PARM'  ,parm(1:70)
   
   ! Echo the input file to the user:
   
   call myechoin(5,6)
   
   !     ----- READ DATA CHARACTERIZING THE MD-RUN -----
   
   read(5,'(20a4)') title
   
   ! read input in namelist format, first setting up defaults
   
   imin = 1
   ntx = 1
   ipb = 2
   inp = 2
   igb = 0

   ! the following are to be retired, but included for backward compatibility

   irest = 0
   ibelly = 0
   idecomp= 0
   ntxo = 1
   ig = 71277
!  tempi = ZERO
   ntt = 0
!  temp0 = 300.0d0

   tautp = ONE
   ntp = 0
   pres0 = ONE
   comp = 44.6d0
   taup = ONE
   npscal = 1
   nscm = 1000
   nstlim = 1
   t = ZERO
   dt = 0.001d0
   ntc = 1
   tol = 0.00001
   ntpr = 50
   ntwr = 500
   ntwx = 0
   ntwv = 0
   ntwe = 0
   ntave = 0
   ioutfm = 0
   ntr = 0
   ntrx = 1
   fcap = 1.5d0
   
   isftrp = 0
   rwell = ONE
   dx0 = 0.01d0
   drms = 1.0d-4
   vlimit = 20.0d0
   mxsub = 1
   ntwprt = 0
!   plevel = 1
   rgbmax = 25.d0
   iyammp = 0
   vrand=1000
   iwrap = 0
   nrespa = 1
   nrespai = 1
   irespa = 1
   gamma_ln = ZERO
   iconstreff = 0
   cut_inner = EIGHT
   icfe = 0
   clambda = ZERO
   klambda = 1
   rbornstat = 0
   lastrst = 2000000
   lastist = 2000000
   restraintmask=''
   restraint_wt = ZERO
   bellymask=''
   saltcon= ZERO
   surften= 0.005d0
   offset = 0.09d0
   intdiel= ONE
   extdiel= 78.5d0
   gbsa   = 0
   ncyc   = 10

   icnstph = 0
   solvph = SEVEN
   ntcnstph = 10
   
   ! check what namelists are defined
   
   mdin_cntrl=.false.
   mdin_pb=.false.
   call mynmlsrc('cntrl',5,ifind)
   if (ifind /= 0) mdin_cntrl=.true.
   call mynmlsrc('pb',5,ifind)
   if (ifind /= 0) mdin_pb=.true.

   ! read cntrl namelist

   rewind 5
   if ( mdin_cntrl ) then
      read(5,nml=cntrl)
   else
      write(6, '(1x,a,/)') 'Error: Could not find cntrl namelist'
      call mexit(6,1)
   end if
   ! disabled cntrl flag
   ntf   =1
   ntb   =1
   dielc  = ONE
   cut    = EIGHT
   nsnb   = 25
   !scnb   = TWO
   !scee   = 1.2d0
   maxcyc = 1
   ntmin  = 1
   ivcap  = 0
   xcap   = 0
   ycap   = 0
   zcap   = 0

   if ( ig == -1 ) then
      call wallclock(rnum)
      ig = nint(rnum)
   end if
 
   ! read pb namelist
 
   if (ipb <= 0) then
      write(6,'(x,a)') "Error: ipb can only be > 0 in PBSA"
      call mexit(6,1)
   else
      call pb_read
   endif

   write(6,9309)

   call printflags()
  

   9308 format(/10x,55('-'),/10x, &
         'GBNSR6 2015                       ', &
         /10x,55('-')/)
   9309 format(/80('-')/'   1.  RESOURCE   USE: ',/80('-')/)
   9700 format(/,'File Assignments:',/,12('|',a6,': ',a,/))

end subroutine mdread1 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize to defaults and print the inputable variables.
subroutine mdread2(x,ix,ih)

   use decomp, only: jgroup, index, irespw

   implicit none

   _REAL_ x(*)
   integer ix(*)
   character(len=4) ih(*)
   integer inerr
   integer nmropt
   integer ntmp, ngrp

#  include "pb_constants.h"
#  include "files.h"
#  include "parms.h"
#  include "memory.h"
#  include "box.h"
#  include "md.h"
   
   close(unit=8)
   
   !     ----- SET THE DEFAULT VALUES FOR SOME VARIABLES -----
   
   nrp = natom
   
   if (ifbox == 1) write(6, '(/5x,''BOX TYPE: RECTILINEAR'')')
   if (ifbox == 2) write(6, '(/5x,''BOX TYPE: TRUNCATED OCTAHEDRON'')')
   if (ifbox == 3) write(6, '(/5x,''BOX TYPE: GENERAL'')')
   
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
   if (ibelly > 0) then   ! No COM Motion Removal, ever.
      nscm = 9999999
      ndfmin = 0
   end if
   if(nscm <= 0) nscm = 9999999
   init = 3
   if (irest > 0) then
!     init = 4
      write(6,'(a)') "We are sorry, but irest have to be 0"
      call mexit(6,0)
   endif
   !if (scnb == ZERO ) scnb = TWO
   if (dielc <= ZERO ) dielc = ONE
   if (tautp <= ZERO ) tautp = 0.2d0
   if (taup <= ZERO ) taup = 0.2d0
   
   ! RESET THE CAP IF NEEDED
   
   if(ivcap == 2) ifcap = 0
   
   ! PRINT DATA CHARACTERIZING THE RUN
   
   nmropt = 0
   iesp = 0
   ipol = 0
  ! write(6,9328)
  ! write(6,9008) title
  ! write(6,'(/a)') 'General flags:'
  ! write(6,'(5x,2(a,i8))') 'imin    =',imin,', nmropt  =',nmropt

   !write(6,'(/a)') 'Nature and format of input:'
   !write(6,'(5x,4(a,i8))') 'ntx     =',ntx,', irest   =',irest, &
   !      ', ntrx    =',ntrx

   !write(6,'(/a)') 'Nature and format of output:'
   !write(6,'(5x,4(a,i8))') 'ntxo    =',ntxo,', ntpr    =',ntpr, &
   !      ', ntrx    =',ntrx,', ntwr    =',ntwr
   !write(6,'(5x,4(a,i8))') 'iwrap   =',iwrap,', ntwx    =',ntwx, &
   !      ', ntwv    =',ntwv,', ntwe    =',ntwe
   !write(6,'(5x,3(a,i8),a,i7)') 'ioutfm  =',ioutfm, &
   !      ', ntwprt  =',ntwprt, &
   !      ', idecomp =',idecomp,', rbornstat=',rbornstat

   !write(6,'(/a)') 'Potential function:'
   !write(6,'(5x,5(a,i8))') 'ntf     =',ntf,', ntb     =',ntb, &
   !      ', ipb     =',ipb,', nsnb    =',nsnb
   !write(6,'(5x,4(a,i8))') 'ipol    =',ipol,', gbsa    =',gbsa, &
   !      ', iesp    =',iesp
   !write(6,'(5x,3(a,f10.5))') 'dielc   =',dielc, &
   !      ', cut     =',cut,', intdiel =',intdiel
   
   !write(6,'(5x,3(a,f10.5))') 'scnb    =',scnb, &
   !      ', scee    =',scee
   
  ! write(6,'(/a)') 'Frozen or restrained atoms:'
  ! write(6,'(5x,4(a,i8))') 'ibelly  =',ibelly,', ntr     =',ntr

  ! if( imin /= 0 ) then
  !    write(6,'(/a)') 'Energy minimization:'
      ! print inputable variables applicable to all minimization methods.
   !   write(6,'(5x,4(a,i8))') 'maxcyc  =',maxcyc,', ncyc    =',ncyc, &
   !         ', ntmin   =',ntmin
   !   write(6,'(5x,2(a,f10.5))') 'dx0     =',dx0, ', drms    =',drms

!     ! Input flag ntmin determines the method of minimization
!     select case ( ntmin )
!     case ( 0, 1, 2 )
!        ! no specific output
!     case default
!        ! invalid ntmin
!        write(6,'(a,i3)') ' ERROR: Invalid NTMIN.', ntmin
!        call mexit(6, 1)
!     end select

!   else
!      write(6,'(/a)') 'Molecular dynamics:'
!      write(6,'(5x,4(a,i8))') 'nstlim  =',nstlim,', nscm    =',nscm, &
!            ', nrespa  =',nrespa
!      write(6,'(5x,3(a,f10.5))') 't       =',t, &
!            ', dt      =',dt,', vlimit  =',vlimit
!   end if

   cut = cut*cut
   cut_inner = cut_inner*cut_inner
   
   ! If user has requested Poisson-Boltzmann electrostatics, set up variables
    
   call pb_init(ifcap,natom,nres,ntypes,nbonh,nbona,ix(i02),ix(i04),ix(i06),ix(i08),ix(i10),&
                ix(iibh),ix(ijbh),ix(iiba),ix(ijba),ix(ibellygp),ih(m02),ih(m04),ih(m06),x(l15),x(l97))


   ntmp = 0; ngrp = 0 
   if(idecomp > 0) then
      write(6,9428)
      call rgroup(natom,ntmp,nres,ngrp,ix(i02),ih(m02), &
            ih(m04),ih(m06),ih(m08),ix(ibellygp), &
            jgroup,index,irespw,npdec, &
            x(l60),x(lcrdr),.false.,.false.,.false.,idecomp,5,.true.)
   end if

   ! checking of not supported options
 
   if( ifcap /= 0 ) then
      write(6, '(a,i3)') ' ERROR: cap water not supported.', ifcap
      call mexit(6, 1)
   end if

   if( igb /= 10 .and. ipb < 1 ) then
      write(6, '(a,i3)') ' ERROR: Non PB potential not supported', igb
      call mexit(6, 1)
   end if

   !     --- checks on bogus data ---
   
   inerr = 0
   
   if (imin /= 0 .and. imin /= 1 .and. imin /= 6) then
      write(6,'(/2x,a,i3,a)') 'IMIN (',imin,') must be 0 or 1 or 6.'
      inerr = 1
   end if
   if (ntx < 1 .or. ntx > 7) then
      write(6,'(/2x,a,i3,a)') 'NTX (',ntx,') must be in 1..7'
      inerr = 1
   end if
   
   !        ---WARNINGS:
   
   if (inerr == 1) then
      write(6, '(/,a)') ' *** input error(s)'
      call mexit(6,1)
   end if

   !  DEBUG input; force checking
   call load_debug(5)
   
   ! Standard format statements:
   
   9328 format(/80('-')/,'   2.  CONTROL  DATA  FOR  THE  RUN',/80('-')/)
   9428 format(/4x,'LOADING THE DECOMP ATOMS AS GROUPS',/)
   9008 format(a80)


end subroutine mdread2 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit defined preprocessor names, ie, flags.
subroutine printflags()

   implicit none
   integer     max_line_length
   parameter ( max_line_length = 80 )

   character(len=max_line_length) line  ! output string of active flags
   integer n                            ! len(line)
   
   line = '| Flags:'
   n = 8
   
   call printflags2(' ',1,n,line,.true.)
   return
end subroutine printflags 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Primitive pre-Fortran90 implementation of printflags.
subroutine printflags2(flag,flag_len,line_len,line,last)

   implicit none
   integer     max_line_length
   parameter ( max_line_length = 80 )

   character(*) flag                ! flag name with blank prefix, intent(in)
   integer flag_len                 ! len(flag), intent(in)
   integer line_len                 ! len(line), intent(inout)
   character(len=max_line_length) line ! intent(inout)
   logical last                     ! is this the last flag ?, intent(in)

   if (line_len + flag_len > max_line_length) then
      write( 6,'(a)') line
      ! begin another line
      line = '| Flags:'
      line_len=8
   end if
   line=line(1:line_len) // flag(1:flag_len)
   line_len=line_len+flag_len
   if(last)write( 6,'(a)') line
   return
end subroutine printflags2 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of a float; abort on illegal values.
subroutine float_legal_range(string,param,lo,hi)
   implicit none
   _REAL_ param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'Ewald PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',e12.5)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',e12.5,' Upper limit: ',e12.5)
   63 format(1x,'Check ew_legal.h')
   return
end subroutine float_legal_range 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of an integer; abort on illegal values.
subroutine int_legal_range(string,param,lo,hi)
   implicit none
   integer param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',i8)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',i8,' Upper limit: ',i8)
   63 format(1x,'The limits may be adjustable; search in the .h files ')
   return
end subroutine int_legal_range 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of an integer option; abort on illegal values.
subroutine opt_legal_range(string,param,lo,hi)
   implicit none
   integer param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'Ewald OPTION CHECKING: ')
   60 format(1x,'option ',a,' has value ',i5)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',i5,' Upper limit: ',i5)
   63 format(1x,'Check the manual')
   return
end subroutine opt_legal_range 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBMD module parsing and initialization
subroutine pb_read
   
   ! Module variables
   
   use variable_module, only : epsin, epsout, epsmemb, smoothopt,   &
         istrng, pbtemp, npbopt, solvopt, accept, maxitn,    &
         fillratio, nbuffer, nfocus, fscale, dbfopt, bcopt,  &
         buffer, xmin, xmax, ymin, ymax, zmin, zmax,         &
         eneopt, frcopt, cutres, cutfd, cutnb, phiform,      &
         cutsa, fmiccg, ivalence, laccept, wsor, lwsor,      &
         pbkappa, offx, offy, offz, sepbuf, mpopt, lmax,     &
         saopt, intopt, ligand, ligandmask, outphi, scalerf, &
         srsas, level, savxm, savym, savzm, savxmym,         &
         savxmymzm, savh, savbcopt, lastp, pbgamma_int,      &
         pbgamma_ext, eps0, fioni, pbkb, isurfchg,           &
         ngrdblkx, ngrdblky, ngrdblkz, xmblk, ymblk, zmblk,  &
         multiblock, saltout, outsalt, stern, sasopt,        &
         membraneopt, mthick, mctrdz, outlvlset, outmlvlset, &
         poretype,poreradius,augtoltype, augctf, augtol,     &
         fmiccg2  
         
   use dispersion_cavity
   use solvent_accessibility
   implicit none
   
   ! Common variables
   
#  include "pb_def.h"
#  include "md.h"
#  include "files.h"
#  include "pb_md.h"
   
     
   ! Local variables
       
   integer npbverb, l, phiout, scalec
   _REAL_ space
    
   ! begin code
 
   namelist /pb/ epsin, epsout, epsmemb, smoothopt, istrng, pbtemp,     &
      radiopt, dprob, iprob, npbopt, solvopt, accept, maxitn,  &
      fillratio, space, nbuffer, nfocus, fscale, npbgrid,      &
      arcres,dbfopt,bcopt,scalec,eneopt,frcopt,       cutfd,   &
      cutnb,        nsnba,phiout, phiform, npbverb, npopt,     &
      decompopt, use_rmin, sprob, vprob, rhow_effect, use_sav, &
      cavity_surften, cavity_offset, maxsph, maxarc,           &
      cutsa, ndofd, ndosas, fmiccg, ivalence, laccept, wsor,   &
      lwsor, pbkappa, radinc, expthresh, offx, offy, offz,     &
      sepbuf, mpopt, lmax, saopt, intopt, ligandmask, buffer,  &
      xmin, ymin, zmin, xmax, ymax, zmax, isurfchg,            &
      ngrdblkx, ngrdblky, ngrdblkz, saltout, stern,            &
      xmblk, ymblk, zmblk, triopt, sasopt, maxtri, maxarcdot,  &
      membraneopt, mthick, mctrdz, outlvlset, outmlvlset,      &
      poretype, poreradius, augtoltype, augtol, augctf,        &
      fmiccg2

   ! initialize with default values
    
   outphi = .false.
   !levelset output file writting controls
   outlvlset = .false.
   outmlvlset = .false.
   phiout = 0
   phiform = 0
    
   !WMBS - membrane options
   membraneopt = 0
   mthick = 20.0d0
   mctrdz = 0.0d0
   poretype = 0
   poreradius = 4.0
   !!!!

   !WMBS - Aug IIM options
   augtoltype=1
   augctf=0.00
   augtol=1.0d-5
   !!!!

   epsin  = 1.0d0
   epsout = 80.0d0
   epsmemb= 1.0d0
   istrng = 0.0d0
   ivalence = 1.0d0
   pbtemp = 300.0d0
    
   scalec = 0
   scalerf = .false.
    
   srsas = .true.
   smoothopt = 1
   radiopt = 1
   radiscale = 1.0d0
   npopt = 2
   decompopt = 2
   use_rmin = 1
   use_sav = 1
   rhow_effect = 1.129d0
   maxsph   = 400
   maxarcdot= 1500
   maxarc   = 512
   maxtri = 10
   nbuffer = 0
   dprob  = 1.40d0
   sprob  = 0.557d0
   vprob  = 1.30d0
   iprob  = 2.00d0
!  optimized radinc = dprob*0.57
!  fixed to be backward compatible
   radinc = 0.8d0
   expthresh = 0.2d0
   arcres = 0.25d0
   cavity_surften = 0.0378d0
   cavity_offset = -0.5692d0
 
   nfocus = 2
   fscale = 8
   level = 1

   bcopt = 5
   space = 0.5d0
   savxm(1) = 0
   savym(1) = 0
   savzm(1) = 0
   savxmym(1) = 0
   savxmymzm(1) = 0
   savh(1) = 0
   savbcopt(1) = 0
   offx = 0.0d0
   offy = 0.0d0
   offz = 0.0d0
   fillratio= 2.0d0
   xmin = 0.0d0
   xmax = 0.0d0
   ymin = 0.0d0
   ymax = 0.0d0
   zmin = 0.0d0
   zmax = 0.0d0
   buffer = 16*0.5 ! this needs to be updated for NAB
   saopt = 0 ! this needs to be updated for NAB
   
   npbopt = 0
   solvopt = 1
   fmiccg = 0.90d0
   fmiccg2= 0.90d0
   accept = 1.0d-3
   laccept = 0.1d0
   wsor = 1.9d0
   lwsor = 1.95d0
   maxitn = 100
   isurfchg = 0
  
   pbverbose = .false. 
   pbprint = .true. 
   pbgrid = .true.
   pbinit = .true.
   npbverb = 0
   npbgrid = 1
   ndofd = 1
   dofd = 1
   donpsa = .true.
   ndosas = 1
   nsaslag = 100
   nsnbr = 1
   nsnba = 1
   ntnbr = 1
   ntnba = 1
   cutres = 99.0d0
   cutfd = 5.0d0
   cutnb = 0.0d0
   cutsa = 9.0d0
   lastp = 0
   pbgamma_int = 1.0
   pbgamma_ext = 65.0
   sepbuf = 4.0d0
   mpopt = 0
   lmax = 80
   dbfopt = -1
   eneopt = -1
   frcopt = 0
   intopt = 1 ! this needs to be updated for NAB
   triopt = 1
   sasopt = 0

   ! Focusing with partial 2nd level region

   ligandmask=''
   outsalt = .false.
   saltout = 0
   stern = 0.0d0
   ngrdblkx = 0
   ngrdblky = 0
   ngrdblkz = 0
   xmblk = 0
   ymblk = 0
   zmblk = 0

   ! reading parameters
     
   if ( mdin_pb ) then
      rewind 5
      read(5, nml=pb) 
   end if
     
   ! check and post processing of options
     
   !WMBS - check membrane options compatibility
   if ( membraneopt > 0 ) then
      if ( membraneopt == 1 ) then
         if ( .not.(bcopt==10) ) then
            write(6,'(/a)') 'PB Warning(): membrane setup only suggested for use with bcopt=10 (periodic)'
         end if
         if (.not. (sasopt==2)) then
            write(6,'(/a)') 'PB Info in pb_read(): membrane only compatible with sasopt=2, sasopt set to 2'
            sasopt=2
         end if
         if ( (.not. (solvopt==7)) .and. (.not. (solvopt==1 .or. solvopt==3))  ) then
            write(6,'(/a)') 'PB Info in pb_read(): membrane requires periodic solver -',&
                       ' solvopt 3 (pcg) or solvopt 7 (fft IIM), using default (piccg) solvopt=1'
            solvopt=1
            if (.not. (ipb==2)) then
               write(6,'(/a)') 'PB Warning in pb_read(): PCG solver only suggested for use with ipb==2'
            end if
         else if (solvopt == 7) then
            if (.not. (ipb==5)) then
               write(6,'(/a)') 'PB Info in pb_read(): FFT solver requires Augmented Variable IIM to function for non-vacuum systems'
               ipb=5
            end if
         else
            if (.not. (ipb==2)) then
               write(6,'(/a)') 'PB Warning in pb_read(): PCG solver not suggested for use except with ipb==2'
            end if
         end if
      else
         write(6,'(/a)') 'PB Bomb in pb_read(): unrecognized membrane option'
         call mexit(6,1)
      end if
   end if

   if ( npbverb > 0 ) pbverbose = .true.

   dofd  = ndofd
   dosas  = ndosas

   epsin  = epsin*eps0
   epsout = epsout*eps0
   epsmemb= epsmemb*eps0
   istrng = fioni * istrng
   pbkappa  = SQRT( 2.0d0 * istrng / (epsout * pbkb * pbtemp) )

   ! check arcres
   if ( arcres > space * 0.5d0 ) then
      arcres = space * 0.5d0
      write(6, '(a,f14.8)') 'PB Warning in pb_read(): arcres was too big and is now set to', arcres
   end if
   
   ! check surface
   if ( ipb == 1 .and. sasopt == 2 ) then
      write(6, '(a)') 'PB Bomb in pb_read(): the smooth Gaussian surface must be built with a level set funtion (ipb > 1)'
      call mexit(6, 1)
   end if
    
   ! set up solvent probes for different solvation components
   ! if end user requests different probes for different solvation
   ! components, no need to set up different probes here.

!  if ( dprob == 0.0d0 ) dprob = sprob

   ! set up solvent accessible arc resolutions and limits
   ! Mengjuei: if maxarcdot is default, automatically set it up.

   if (maxarcdot == 1500) maxarcdot = 1500*nint(0.5d0/arcres)

   ! check dprob and iprob

   if ( dprob > iprob ) then
      if ( sasopt == 1 ) then
         write(6, '(a)') 'PB Bomb in pb_read(): ion probe cannot be smaller than solvent prob when the SAS is used (sasopt = 1)'
         call mexit(6, 1)
      else
         write(6, '(a)') 'PB Warning in pb_read(): ion probe is smaller than solvent prob'
      end if
   end if
   if ( dprob <= space .and. ipb == 1 ) then
      write(6, '(a)') 'PB Warning in pb_read(): setting grid spacing larger than solvent probe'
      write(6, '(a)') 'may cause numerical instability if ipb=1'
   end if
   if ( dprob <= space .and. ipb == 2 ) then
      write(6, '(a)') 'PB Bomb in pb_read(): solvent probe cannot be smaller than grid spacing if ipb=2'
      call mexit(6, 1)
   end if

   ! check pb options

   ! Paranoid Check by Mengjuei

   if ( ipb < 1 ) then
      write(6, '(a)') 'PB Bomb in pb_read(): ipb can only be larger than or equal to 1.'
      call mexit(6, 1)
   else if (igb == 10 .and. ipb == 0) then
      ipb=2
   else if (igb /= 0 .and. igb /= 10) then
!     write(6, *) 'Error: igb can only be 10 in PBSA'
!     write(6, *) 'PB Info in pb_read(): igb has been overwritten by ipb'
      igb=10
   else if ( igb == 10 .and. ipb /= 2 ) then
!     write(6, *) 'PB Info in pb_read(): igb has been overwritten by ipb'
   end if

   if ( ipb == 4 ) then
      if ( bcopt < 4 ) then
         bcopt = 6
         write(6, '(a)') 'PB Info in pb_read(): bcopt has been reset to 6 with ipb=4'
      end if
      if ( nfocus /= 1 ) then
         nfocus = 1
         write(6, '(a)') 'PB Info in pb_read(): nfocus has been reset to 1 with ipb=4'
      end if
!     if ( frcopt == 1 ) then 
!        write(6, *) 'PB Info in pb_read(): frcopt=1 conflicts with ipb=4' 
!        call mexit(6, 1)
!     end if
      if ( istrng /= 0 ) then
         write(6, '(a)') 'PB Info in pb_read(): ipb=4 cannot deal with salt solution' 
         call mexit(6, 1)
      end if
   end if

   ! check focusing options

   if ( nfocus > 1 .and. fscale < 2 ) then
      write(6, '(a)') 'PB Info in pb_read(): nfocus reset to 1 due small fscale.'
      nfocus = 1
   endif
   
   ! check nonpolar options

   if ( inp /= npopt ) then
!     write(6, *) 'PB Info in pb_read(): npopt has been overwritten with inp.'
      npopt = inp
   endif
!  if ( inp == 2 .and. radiopt /= 1 ) then
!     write(6, *) 'PB Warning in pb_read(): radiopt has been reset to 1 because inp=2.'
!     radiopt = 1
!  end if
!   if ( inp == 1 .and. use_sav == 1 ) then
!      write(6, '(a)') &
!         'PB Warning in pb_read(): the cavity energy should be correlated to &
!         & the SASA when inp=1, so use_sav is reset to 0.'
!      use_sav = 0
!   end if
   if ( inp /= 2 .and. use_rmin == 1 ) then
!     write(6, *) 'PB Warning in pb_read(): use van der Waals sigma values when inp/=2, &
!                & so use_rmin is reset to 0.'
      use_rmin = 0
   end if
!   if ( inp == 1 .and. sprob == 0.557d0 ) then
!      write(6, '(a)') &
!         'PB Warning in pb_read(): sprob=0.557 is optimized for inp=2 and &
!         & should not be used with inp=1. It has been reset to 1.4.'
!      sprob = 1.4d0 
!   end if
!   if ( inp == 1 .and. cavity_surften == 0.0378d0 ) then
!      write(6, '(a)') &
!         'PB Warning in pb_read(): cavity_surften=0.0378 is optimized for inp=2 &
!         & and should not be used with inp=1. It has been reset to 0.04356.'
!      cavity_surften =  0.0072d0
!   end if
!   if ( inp == 1 .and. cavity_offset == -0.5692d0 ) then
!      write(6, '(a)') &
!         'PB Warning in pb_read(): cavity_offset=-0.5692 is optimized for inp=2 &
!         & and should not be used with inp=1. It has been reset to -1.008.'
!      cavity_offset =  -1.008d0
!   end if
!  if ( radiopt == 0 ) then
!     donpsa = .false.
!  end if

   ! check force options

   if ( eneopt == -1 ) then
      if ( dbfopt == 0 ) then 
         eneopt = 1
      else if ( dbfopt == 1 .or. dbfopt == -1 ) then
         eneopt = 2
      else
         write(6, '(a)') 'PB Info in pb_read(): only dbfopt = 0 or 1 are supported'
         write(6, '(a)') 'PB Info in pb_read(): dbfopt is replaced by eneopt'
         call mexit(6, 1)
      end if
   else
      if ( dbfopt == 0 .and. eneopt == 1 ) then 
!        write(6, *) 'PB Info in pb_read(): dbfopt has been overwritten by eneopt'
      else if ( dbfopt == 1 .and. eneopt == 2 ) then 
!        write(6, *) 'PB Info in pb_read(): dbfopt has been overwritten by eneopt'
      else if ( dbfopt /= -1 ) then
         write(6, '(a)') 'PB Info in pb_read(): dbfopt is ignored when eneopt is set'
      end if
   end if
   if ( eneopt < 1 .or. eneopt > 3 ) then
      write(6, '(a)') 'PB Bomb in pb_read(): only eneopt= 1-3 are supported'
      call mexit(6, 1)
   end if
   if ( frcopt < 0 .or. frcopt > 4 ) then
      write(6, *) 'PB Bomb in pb_read(): only frcopt= 0-4 are supported'
      call mexit(6, 1)
   end if
   if ( eneopt == 1 .and. ( frcopt >= 2 ) ) then
      write(6, '(a,2i2)') 'PB Bomb in pb_read(): combination of eneopt and frcopt is not supported', eneopt, frcopt
      call mexit(6, 1)
   end if
   if ( eneopt == 2 .and. frcopt == 1 ) then
!     write(6, *) 'PB Bomb in pb_read(): combination of eneopt and frcopt is not supported', eneopt, frcopt
!     call mexit(6, 1)
   end if
   if ( ( frcopt >= 2 .and. frcopt <= 4 ) .and. bcopt == 5 ) then
      bcopt = 6
      write(6, '(a)') 'PB Info in pb_read(): bcopt has been reset from 5 to 6 for frcopt= 2 - 4'
   end if
   if ( frcopt > 0 .and. smoothopt == 0 ) then
!     smoothopt = 1
!     write(6, *) 'PB Info in pb_read(): smoothopt has been reset to 1 for force compuation'
   end if

   if ( ivcap /= 0 .and. (eneopt /= 1 .or. frcopt /= 0) ) then
!     write(6,*) "cap water should only be run with eneopt=1 and frcopt=0"
!     call mexit(6,1)
   end if

   if ( saopt /= 0 .and. abs(saopt) /= 1 .and. abs(saopt) /= 2 ) then
      saopt = 0
      write(6, '(a)') 'PB Info in pb_read(): saopt has been reset to 0 for surface area compuation'
   end if
    
   ! check numerical solution options

   if ( solvopt == 7 ) then
      continue
   else if ( npbopt == 0 .and. solvopt > 5 .and. (.not. solvopt == 8) ) then
      write(6, '(a)') 'PB Bomb in pb_read(): solvopt>5 cannot be used to solve linear PB equation'
      call mexit(6, 1)
   endif
   if ( ipb /= 4 .and. solvopt == 2 .and. bcopt == 6 ) then
      write(6, '(a)') 'PB Bomb in pb_read(): bcopt=6 cannot be used with solvopt=2'
      call mexit(6, 1)
   end if
   if ( npbopt == 1 .and. bcopt == 10 ) then
      write(6, '(a)') 'PB Bomb in pb_read(): bcopt=10 can be used only with npbopt=0'
      call mexit(6, 1)
   end if
   if ( solvopt /= 3 .and. solvopt /=8 .and. solvopt /=1 .and. bcopt == 10 ) then
      write(6, '(a)') 'PB Bomb in pb_read(): bcopt=10 can be used only with solvopt=1,3, or 8'
      call mexit(6, 1)
   end if
   if ( npbopt == 1 .and. solvopt > 6 ) then
      write(6, '(a)') 'PB Bomb in pb_read(): nonsupported solvopt (>6) for e nonlinear PB equation'
      call mexit(6, 1)
   endif
   if ( npbopt == 1 .and. eneopt == 2 ) then
      eneopt = 1
      write(6, '(a)') 'PB Info in pb_read(): eneopt has been reset to be 1 for nonlinear PB equation'
   end if
   if ( npbopt == 1 .and. frcopt > 0 ) then
      write(6, '(a)') 'PB Bomb in pb_read(): force computatoin is not supported for nonlinear PB equation'
      call mexit(6, 1)
   end if

   ! check cutoff options
    
   if ( cutfd > cutsa ) then
      cutsa = cutfd
   end if
   if ( cutnb /= 0 .and. cutfd > cutnb ) then
      cutnb = cutfd
      write(6, '(a)') 'PB Info in pb_read(): cutnb has been reset to be equal to cutfd'
   end if
   if ( max(cutnb,cutsa,cutfd) > cutres ) then
      write(6, '(a,3f8.4)') 'PB Bomb in pb_read(): cutnb/cutfd must be <= cutres', cutnb, cutfd, cutres
      call mexit(6, 1)
   end if
   if ( cutnb /= 0 .and. bcopt > 5 .and. bcopt < 10 ) then
      cutnb = 0
      write(6, '(a,2f8.4)') 'PB Info in pb_read(): cutnb has been reset to 0 with 5 < bcopt < 10', cutnb, bcopt
   end if
   if ( cutnb == 0 .and. eneopt == 1 .and. bcopt < 6 ) then
      write(6, '(a,2f8.4)') 'PB Bomb in pb_read(): cutnb=0 cannot be used with eneopt=1', cutnb, eneopt
      call mexit(6, 1)
   end if
   cutfd = cutfd**2
   cutsa = cutsa**2
   cutnb = cutnb**2
   cutres = cutres**2
    
   ! set buffer zone between the fine FD grid boundary and the solute surface:
    
   if ( nbuffer == 0 ) then
      if ( istrng == 0.0d0 ) then
         nbuffer = int(2.0d0*dprob/space)+1
      else
         if ( dprob >= iprob ) then
            nbuffer = int(2.0d0*dprob/space)+1
         else
            nbuffer = int(2.0d0*iprob/space)+1
         end if
      end if
      if ( nbuffer >= fscale ) then
         nbuffer = 2*nbuffer+1
      else
         nbuffer = 2*fscale+1
      end if
   end if
!if ( npbverb > 0 .or. pbverbose ) then
!   print *, ' dprob =', dprob,'    sprob =',sprob
!   print *, 'eneopt =',eneopt,'   dbfopt =',dbfopt
!   print *, 'frcopt =',frcopt,'smoothopt =',smoothopt
!   print *, 'npbopt =',npbopt,'  solvopt =',solvopt
!   print *, ' bcopt =', bcopt,'  nbuffer =',nbuffer
!   print *, ' cutfd =', cutfd,'    cutsa =',cutsa
!   print *, ' cutnb =', cutnb,'   cutres =',cutres
!end if
    
   ! set flag to scale induced surface charges:
    
   if ( scalec == 1) scalerf = .true. 
    
   ! set saved grid options

   savxm(nfocus) = 0
   savym(nfocus) = 0
   savzm(nfocus) = 0
   savxmym(nfocus) = 0
   savxmymzm(nfocus) = 0
   savh(nfocus) = 0
   savbcopt(nfocus) = 0

   do l = 1, nfocus
      savbcopt(l) = bcopt
   end do
   savh(nfocus) = space
   do l = nfocus - 1, 1, -1
      savh(l) = savh(l+1)*fscale
   end do

   do l = 1, nfocus
      savxm(l) = 1
      savym(l) = 1
      savzm(l) = 1
      savxmym(l) = 1
      savxmymzm(l) = 1
   end do

   if ( nfocus > 2 ) then !for sanity
      write(6,'(a)') "Focusing with more than two levels is not implemented"
      call mexit(6,1)
   end if

   ! set MULTIBLOK focusing option
   multiblock = .false.
   if ( ngrdblkx < 1 .or. ngrdblky < 1 .or. ngrdblkz < 1 ) then
      if ( xmblk > 0 .or.    ymblk > 0 .or.    zmblk > 0 ) then
         write(6,'(a)') "pb_read: [x,y,z]mblk require non-zero ngrdblk[x,y,z]."
         write(6,'(a)') "pb_read: [x,y,z]mblk are reset to zero."
         xmblk = 1; ymblk = 1; zmblk = 1
      end if 
   else if ( mod(ngrdblkx-1,fscale) /= 0 .or. &
             mod(ngrdblky-1,fscale) /= 0 .or. &
             mod(ngrdblkz-1,fscale) /= 0       ) then
         write(6,'(a)') "pb_read: wrong ngrdblk[x,y,z] - fscale combination"
         call mexit(6,1)
   else if ( nfocus /= 2 ) then
      write(6,'(a)') "pb_read: Multiblock focusing need nfocus = 2."
      call mexit(6,1)
   else if ( frcopt /= 0 ) then
      write(6, '(a)') 'pb_read: Multiblock is incompatible with frcopt /= 0.'
      call mexit(6,1)
   else if ( eneopt == 2 ) then
      write(6, '(a)') 'pb_read: Multiblock is incompatible with eneopt == 2.'
      call mexit(6,1)
   else if ( ivcap /= 0 ) then
      write(6, '(a)') 'pb_read: Multiblock is incompatible with cap water.'
      call mexit(6,1)
   else
      multiblock = .true.
      ! nbuffer is an odd number
      buffer  = max(savh(nfocus)*(nbuffer-1)/2,buffer)
      !nbuffer = max(nbuffer,nint(buffer/savh(nfocus)*2))
   end if

   ! set up ligand focusing option
   
   ligand = .false.
   if ( len_trim(ligandmask) > 0 ) ligand = .true.
   if ( xmax-xmin > 0 .and. &
        ymax-ymin > 0 .and. &
        zmax-zmin > 0 ) ligand = .true.

   if (ligand) then
      if ( frcopt /= 0 ) then
         write(6, '(a)') 'pb_read: Ligandmask is not compatible with frcopt /= 0.'
         call mexit(6,1)
      else if ( eneopt == 2 ) then
         write(6, '(a)') 'pb_read: Ligandmask is not compatible with eneopt == 2.'
         call mexit(6,1)
      else if ( nfocus < 2 ) then
         write(6, '(a)') 'pb_read: You cannot use ligandmask without focusing.'
         call mexit(6,1)
      else if ( multiblock ) then
         write(6, '(a)') 'pb_read: You cannot use multiblock with ligandmask.'
         call mexit(6,1)
      end if
      !buffer = max(space*nbuffer/2,buffer)
      !nbuffer = max(nbuffer,nint(buffer/space*2))
   end if

   ! set phimap output options when requested

   if ( phiout == 1 ) then
      outphi = .true.
      radiopt = 2
      donpsa = .false.
      npopt = 0
   end if

   ! set saltmap output

   if ( saltout == 1 ) outsalt = .true.

end subroutine pb_read 

subroutine myresmask ( maskstr, masklength, resmask, nres )
   ! Mengjuei Hsieh
   implicit none
   character(*) maskstr
   integer      resmask(*)
   integer      masklength, nres

   ! Local
   integer      strpos
   integer      i, range_start, range_end
   logical      assigned
   logical      res_mode, num_mode, rng_mode

   i=0
   res_mode=.false.
   num_mode=.false.
   rng_mode=.false.
   assigned=.false.
   ! mjhsieh: warning elminator
   range_start = -1
   range_end = -1

   do strpos = 1, masklength
      if ( maskstr(strpos:strpos)==' ' ) cycle
      if ( maskstr(strpos:strpos)==':' .and. res_mode ) then
         print *,"Error parsing residue mask: dupe commas "; call mexit(6,1)
      else if ( maskstr(strpos:strpos) == ':' ) then
         res_mode = .true.
         cycle
      else if ( .not. res_mode ) then
         print *,"Error parsing residue mask: lack commas "; call mexit(6,1)
      else if ( maskstr(strpos:strpos) == '@' ) then
         print *,"mjhsieh: unsupported mask syntax "; call mexit(6,1)
      else if ( maskstr(strpos:strpos) >= 'A' .and. &
                maskstr(strpos:strpos) <= 'z' ) then
         print *,"mjhsieh: unsupported mask syntax "; call mexit(6,1)
      else if ( maskstr(strpos:strpos) >= '0' .and. &
                maskstr(strpos:strpos) <= '9' ) then
         num_mode = .true.
         assigned = .false.
         i = i*10 + ichar(maskstr(strpos:strpos))-ichar('0')
         if ( i > nres ) then 
            print *,"mjhsieh: mask out of range"; call mexit(6,1)
         end if 
      else if ( maskstr(strpos:strpos) == '-' ) then
         if ( rng_mode ) then
            print *,"mjhsieh: mask syntax error"
         else if ( num_mode ) then
            rng_mode = .true.
            num_mode = .false.
            range_start = i; i = 0
         else
            print *,"mjhsieh: mask syntax error"
         end if
      else if ( maskstr(strpos:strpos) == ',' ) then
         if ( rng_mode ) then
            range_end = i; i = 0
            resmask(range_start:range_end)=1
            rng_mode = .false.
         else if ( num_mode ) then
            resmask(i:i) = 1; i = 0
            assigned=.true.
         else
            print *,"mjhsieh; mask syntax error"
         end if
      else
         print *,"mjhsieh: unsupported mask syntax "; call mexit(6,1)
      end if
   end do
   if ( assigned ) return

   if ( .not. rng_mode ) range_start=i
   range_end=i
   if ( range_start > range_end ) then
      print *,"reverse number order in mask, are you sure?"; call mexit(6,1)
   end if
   resmask(range_start:range_end)=1
   !print *, range_start,"-",range_end
end subroutine myresmask

subroutine myechoin (ilun,iout)
   ! Rewritten by: Meng-Juei Hsieh
   ! A rewrite of echoin
   implicit none
   integer ilun, iout

   logical boolinside, boolhastag
   character(len=80) mylnbuf
   boolinside=.false.
   boolhastag=.false.

   do while (.true.)
      read(ilun,'(80a)',end=1949) mylnbuf
      if (boolinside) then
         if (mylnbuf(1:14)=='!! BEGIN Amber') then
            write(iout,"('Error: Nested interface tags')")
            call mexit(iout,1)
         elseif (mylnbuf(1:14)=='!! END   Amber') then
            boolinside=.false.
         else
            write(iout,'(80a)') mylnbuf(1:79)
         endif
      elseif (mylnbuf(1:14)=='!! BEGIN Amber') then
         write(iout,'()')
         write(iout,"(' The Interface script used to generate the input file:')")
         write(iout,'()')
         boolinside=.true.
         boolhastag=.true.
      elseif (mylnbuf(1:14)=='!! END   Amber') then
         write(iout,"('Error: unclosed interface tag')")
         call mexit(iout,1)
      endif
   enddo
   
   1949 rewind(ilun)
   if (boolhastag) write(iout,"(79('-'))")
   write(iout,'()')
   write(iout,"(' Here is the input file:')")
   write(iout,'()')

   boolinside=.false.
   do while (.true.)
      read(ilun,'(80a)',end=9562) mylnbuf
      if (boolinside) then
         if (mylnbuf(1:14)=='!! END   Amber') then
            boolinside=.false.
         endif
      elseif (mylnbuf(1:14)=='!! BEGIN Amber') then
         boolinside=.true.
      else
         write(iout,'(80a)') mylnbuf(1:79)
      endif
   enddo

   9562 rewind(ilun)

   return
end subroutine myechoin

subroutine mynmlsrc(srchkey,ilun,ifound)
   ! Author: Meng-Juei Hsieh
   implicit none
   character(len=*) srchkey
   integer ilun

   character(len=80) realkey
   character(len=80) mylnbuf
   integer ifound, ilen, ipos

   ilen = len(srchkey)
   realkey(2:ilen+1)=srchkey(1:ilen)
   ifound = 0
   ipos = 0

   do while (ifound == 0)
      read(ilun,'(80a)',end=824,err=9528) mylnbuf
      realkey(1:1)='&'
      ipos=index(mylnbuf, realkey(1:ilen+1), .false.)
      if (ipos>0) then
         ifound = 1
      else
         realkey(1:1)='$'
         ipos=index(mylnbuf, realkey(1:ilen+1), .false.)
         if (ipos>0) ifound = 1
      endif
   enddo

   824 continue
   !write(6,*) realkey(1:ilen+1)

   if (ifound == 1 ) then
      backspace(ilun)
   else if (ifound == 0) then
      rewind(ilun)
   else
      write(6,'(a)') "mynmlsrc Error: exception"
      call mexit(6,1)
   endif
   return

   9528 rewind(ilun)
   return
end subroutine mynmlsrc
