! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Amber runmd place holder, only one step force call and printing
subroutine runmd(xx,ix,ih,ipairs,x,winv,amass,f, &
                 v,vold,xr,xc,conp,skip,nsp,tma,erstop,r_stack,i_stack)
   
   use memory_module

   implicit none

! Passed variables
   integer   ipairs(*)
   _REAL_ xx(*)
   integer   ix(*)
   character(len=4) ih(*)
   _REAL_ r_stack(*)
   integer i_stack(*)
 
#  include "files.h"
#  include "md.h"
#  include "box.h"
#  include "extra.h"
#  include "timer.h"
#  include "pb_constants.h"
#  include "pb_md.h"

!  character(len=6)fnam

!  logical resetvelo
!  integer nshak
!  _REAL_ etot_save,ekpbs,ekgs,eold3,eold4
   _REAL_ ekpbs
   
   logical do_list_update
!  logical skip(*),belly,lout,loutfm,erstop,vlim,onstep
   logical skip(*),belly,lout,erstop
   _REAL_ x(*),winv(*),amass(*),f(*),v(*),vold(*), &
         xr(*),xc(*),conp(*),vol
   _REAL_ enert(51),enert2(51),ener(51),vir(4),ekcmt(4)
!  _REAL_ enert_old(51),enert2_old(51),ecopy(51),edvdl(51), &
!        enert_tmp(51),enert2_tmp(51),etot_start,edvdl_r(51)
   _REAL_                                        edvdl(51)
   _REAL_ pres(4),fac(3)!,vircopy(3),clfac,rmu(3)
   _REAL_ tma(*)

   _REAL_ tspan,scaltp!,atempdrop,fln
!  _REAL_ vel,vel2,vcmx,vcmy,vcmz,vmax,vx,vy,vz
!  _REAL_ ekmh,ekph,wfac,winf,aamass,ekpht,ekav,rsd,rterm
   _REAL_ ekmh,ekph
!  _REAL_ fit,fiti,fit2
!  _REAL_ gammai,c_implic,c_explic,c_ave,sdfac,ekins0
!  _REAL_ dtx,dtxinv,dt5,factt,ekin0,ekinp0,dtcp,dttp
   _REAL_ dtx
!  _REAL_ rndf,rndfs,rndfp,onet,boltz2,pconv,ibelsv,tempsu
   _REAL_ rndf,rndfs,rndfp,boltz2,pconv
!  _REAL_ ekrot,ekcm,acm(3),ocm(3),vcm(3),xcm(3)

   integer nsp(*)
!  integer idumar(4)
!  integer l_temp
   integer i,m,nitp,nits!,i3,im,j
!  integer nstep,nrep,nrek,nren,iend,istart3,iend3
   integer nstep,nrek,nren
!  integer nrx,nr,nr3,ntcmt,izero,istart
   integer nrx,nr,nr3,izero
!  logical ixdump,ivarch,itdump
   
   equivalence (scaltp,ener(5)),(vol,ener(10))
   equivalence (pres(1),ener(11)),(ekcmt(1),ener(15))
   equivalence (vir(1),ener(19))
   integer nvalid
   _REAL_ eke!,eket
!  _REAL_ extent

!  _REAL_ xcen,ycen,zcen,extents(3,2),centertest
   integer ier
   
!  _REAL_ small
!  data small/1.0d-7/
   data nren/51/
   
   !--- VARIABLES FOR DIPOLE PRINTING ---
!  integer prndipngrp
!  integer prndipfind
!  character(len=4) prndiptest
   !--- END VARIABLES FOR DIPOLE PRINTING ---

   !  Runmd operates in kcal/mol units for energy, amu for masses,
   !     and angstoms for distances.  To convert the input time parameters
   !     from picoseconds to internal units, multiply by 20.455
   !     (which is 10.0*sqrt(4.184)).
   
   !==========================================================================
   
   !     ----- INITIALIZE SOME VARIABLES -----
   
   if( master ) call myopen(MDINFO_UNIT,mdinfo,'U','F','W')
!  vlim = vlimit > small
!  ntcmt = 0
   izero = 0
   belly = ibelly > 0
   lout = .true.
!  loutfm = ioutfm <= 0
   nr = nrp
   nr3 = 3*nr
   ekmh = 0.d0
!  istart = 1
!  iend = natom
!  istart3 = 3*istart -2
!  iend3 = 3*iend

   ! If NTWPRT.NE.0, only print the atoms up to this value

   nrx  = nr3
   if (ntwprt > 0) nrx = ntwprt*3
   
   ! Cleanup the velocity if belly run

   if(belly) call bellyf()
   
   ! Determine system degrees of freedom (for T scaling, reporting)
   
   ! RNDFP = # degrees of freedom for solute
   ! RNDFS = # degrees of freedom for solvent
   ! RNDF = total number of degrees of freedom.
   
   rndfp = 1
   rndfs = 0
   rndf = 1
   
!  onet = 1.d0/3.d0
   
   boltz2 = 8.31441d-3 * 0.5d0
   pconv = 1.6604345d+04  ! factor to convert the pressure kcal/mole to bar
   
   !     ---convert to kcal/mol units
   
   boltz2 = boltz2/4.184d0   ! k-sub-B/2
   dtx = dt*20.455d+00
!  dtxinv = 1.0d0 / dtx
!  dt5 = dtx * 0.5d0
   pconv = pconv*4.184d0
   
   ! FAC() are #deg freedom * kboltz / 2
   ! multiply by T to get expected kinetic energy
   ! FAC(1) is for total system
   
   fac(1) = boltz2*rndf
   fac(2) = boltz2*rndfp
   if(rndfp < 0.1d0) fac(2) = 1.d-6

   fac(3) = boltz2*rndfs
   if(rndfs < 0.1d0) fac(3) = 1.d-6
!  factt = rndf/(rndf+ndfmin)
   
   ! these are "desired" kinetic energies based on
   ! # degrees freedom and target temperature
   ! they will be used for calculating the velocity scaling factor
   
!  ekin0  = fac(1)*temp0
!  ekinp0 = fac(2)*temp0
!  ekins0 = fac(3)*temp0

   ! LD setup:
   
!  gammai = gamma_ln/20.455d0
!  c_implic = 1.d0/(1.d0+gammai*dt5)
!  c_explic = 1.d0 - gammai*dt5
!  c_ave    = 1.d0+gammai*dt5
   !sdfac = sqrt( 4.d0*gammai*boltz2*temp0/dtx )

   nrek = 4
!  nrep = 15
   
   nvalid = 0
   nstep = 0
!  fit = 0.d0
!  fiti = 0.d0
!  fit2 = 0.d0

   do i = 1,nren
      ener(i) = 0.0d0
      enert(i) = 0.0d0
      enert2(i) = 0.0d0
!     enert_old(i) = 0.d0
!     enert2_old(i) = 0.d0
      edvdl(i) = 0.d0
!     edvdl_r(i) = 0.d0
   end do
   
   ener(5) = 1.d0
   ener(6) = 1.d0
   do m = 1,3
      ener(m+6) = box(m)
   end do
   
   nitp = 0
   nits = 0
   
   !=======================================================================
   !     ----- MAKE A FIRST DYNAMICS STEP -----
   !=======================================================================

   !  init = 3:  general startup if not continuing a previous run
   
   if( init == 3 ) then
      
      ! ----- CALCULATE THE FORCE -----

      npbstep = nstep + 1

      call force(xx,ix,ih,ipairs,x,f,ener(23),vir, &
              r_stack,i_stack, xx(l96), xx(l97), xx(l98), do_list_update)
      
      ener(1) = ener(2)+ener(23)

      if(ntt == 1) then
         ekmh = max(ener(3),fac(1)*10.d0)
      end if

      ! PRINT THE INITIAL ENERGIES AND TEMPERATURES
      
      if (nstep <= 0 .and. master) then
          rewind(MDINFO_UNIT)
          call prntmd(nstep,nitp,nits,t,ener,fac,MDINFO_UNIT,.false.)
          call amflsh(MDINFO_UNIT)
      end if
      
   end if  ! ( init == 3 )
   
   !-------------------------------------------------------------------------
   ! init = 4: continuation of a previous trajectory. this code is also done for init=3
   !
   ! Note: if the last printed energy from the previous trajectory was
   !       at time "t", then the restrt file has velocities at time
   !       t + 0.5dt, and coordinates at time t + dt
   !-------------------------------------------------------------------------
   
   ekmh = 0.0d0
   
   if (nstlim == 0) return
   init = 4
   
   !---------------------------------------------------------------
   !  ---Step 1a: do some setup for pressure calculations:
   !---------------------------------------------------------------
   
   !--------------------------------------------------------------
   !  ---Step 1b: Get the forces for the current coordinates:
   !--------------------------------------------------------------
   
   npbstep = nstep + 1

   call force(xx,ix,ih,ipairs,x,f,ener(23),vir, &
           r_stack,i_stack,xx(l96),xx(l97),xx(l98),do_list_update)

   ! Reset quantities depending on TEMP0 and TAUTP (which may have been
   ! changed by MODWT during FORCE call).
   
   !  Pressure coupling:
   
   !----------------------------------------------------------------
   !  ---Step 1c: do randomization of velocities, if needed:
   !----------------------------------------------------------------

   !-----------------------------------------------------
   !  ---Step 2: Do the velocity update:
   !-----------------------------------------------------
   
   !-------------------------------------------------------------------
   !   Step 3: update the positions, putting the "old" positions into F:
   !-------------------------------------------------------------------
   
   if (ntc /= 1) then
      write(6,'(a)') 'SHAKE OPTIONS ARE DISABLED'
      call mexit(6,1)
   end if

!  if( ntt == 1 .or. onstep ) then
      
      !-----------------------------------------------------------------
      !   Step 4c: get the KE, either for printing or for Berendsen:
      !-----------------------------------------------------------------
      
      eke = 0.d0
      ekph = 0.d0
      ekpbs = 0.d0
      
!  end if  ! ( ntt == 1 .or. onstep; end of step 4c )
 
   !-----------------------------------------------------------------
   !   Step 5:  several tasks related to dumping of trajectory information
   !-----------------------------------------------------------------

   !  --- determine if restart writing is imminent and
   !      requires xdist of v and dipole information in parallel runs:
 
!  ivarch = .false.
!  ixdump = .false.
!  itdump = .false.
 
   !-------------------------------------------------------------------
   !   Step 6: zero COM velocity if requested; used for preventing
   !   ewald "block of ice flying thru space" phenomenon, or accumulation
   !   of rotational momentum in vacuum simulations
   !-------------------------------------------------------------------
 
   ntnb = 1
   
   !  Also zero out the non-moving velocities if a belly is active:

   if (belly) call bellyf()
   
   !-------------------------------------------------------------------
   !  Step 7: scale coordinates if constant pressure run:
   !-------------------------------------------------------------------
   
   ener(4) = ekpbs + ener(23)
   ener(3) = eke
   ener(2) = ener(3)
   
   ener(1) = ener(2)+ener(23)
!  etot_save = ener(1)
   
   !-------------------------------------------------------------------
   !  Step 8:  update the step counter and the integration time:
   !-------------------------------------------------------------------
   
   nstep = nstep+1
   t = t+dt
   nvalid = nvalid+1
   
   !     ---full energies are only calculated every nrespa steps
   !     nvalid is the number of steps where all energies are calculated
   
   ntnb = 1
   lout = .true.
     
   ! reset pb-related flags

   if ( ipb >= 1 ) then
      if ( mod(nstep+1,npbgrid) == 0 .and. nstep+1 /= nstlim ) pbgrid = .true.
      if ( mod(nstep+1,ntpr) == 0 .or. nstep+1 == nstlim ) pbprint = .true.
      if ( mod(nstep+1,nsnbr) == 0 .and. nstep+1 /= nstlim ) ntnbr = 1
      if ( mod(nstep+1,nsnba) == 0 .and. nstep+1 /= nstlim ) ntnba = 1
   end if

   !-------------------------------------------------------------------
   !  Step 9:  output from this step if required:
   !-------------------------------------------------------------------
   
   !     ...only the master needs to do the output
   
   if (master) then
      if (lout) then
         rewind(MDINFO_UNIT)
         call prntmd(nstep,nitp,nits,t,ener,fac,MDINFO_UNIT,.false.)
         call amflsh(MDINFO_UNIT)
      end if
   end if  ! (master)
   
   !=======================================================================
   !     ----- PRINT AVERAGES -----
   !=======================================================================
   
   if (master) then
      tspan = nvalid
      if (nvalid > 0) then
         do m = 1,nren
            enert(m) = enert(m)/tspan
            enert2(m) = enert2(m)/tspan - enert(m)*enert(m)
            if(enert2(m) < 0.d0) enert2(m) = 0.d0
            enert2(m) =  sqrt(enert2(m))
            edvdl(m) = edvdl(m)/tspan
         end do
         
         write(6,540) nvalid
         call prntmd(nstep,izero,izero,t,enert,fac,0,.false.)
         write(6,550)
         call prntmd(nstep,izero,izero,t,enert2,fac,0,.true.)
         
         do m = 2,nrek
            enert(m) = enert(m)/fac(m-1)
            enert2(m) = enert2(m)/fac(m-1)
         end do
         temp = enert(2)
      end if  ! (nvalid > 0)
      
      close(unit=MDINFO_UNIT)
   end if  ! (master)
   
   540 format(/5x,' A V E R A G E S   O V E R ',i7,' S T E P S',/)
!  542 format('|',79('='))
   550 format(/5x,' R M S  F L U C T U A T I O N S',/)
!  580 format('STATISTICS OF EFFECTIVE BORN RADII OVER ',i7,' STEPS')
!  590 format('ATOMNUM     MAX RAD     MIN RAD     AVE RAD     FLUCT')
!  600 format(i4,2x,4f12.4)


end subroutine runmd 
