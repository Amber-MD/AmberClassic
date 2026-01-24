#include "copyright.h"
#include "../include/dprec.fh"

subroutine gb_read(natom) 

  
   
   use genborn, only : gb_alpb, gb_arad,gb_chunk, gb_B, gb_depth, gb_Rs, gb_tau, gb_ROH, gb_cha, gb_dgij
   use variable_module 
   use solvent_accessibility
   use dispersion_cavity

   implicit none

#  include "pb_def.h"
#  include "md.h"
#  include "pb_md.h"

   !character(len=4) ih(*)
   integer natom 

   ! Local variables
   integer npbverb, l, phiout, scalec
   _REAL_ space
   _REAL_ B
   _REAL_ Rs
   _REAL_ tau
   _REAL_ ROH
   integer chagb
   integer depth
   integer ifind
   logical mdin_gb
   integer dgij

   namelist /gb/ epsin, epsout, istrng,              &
      dprob, space, arcres,                   &
      alpb,  depth, B, rbornstat, depth,      &
      cavity_surften , Rs , tau , ROH, chagb, dgij, radiopt
  
   ! initialized the default gb parameters
   B = 0.028d0
   Rs = 0.52d0
   tau = 1.47d0
   ROH = 0.58588d0
   alpb = 1
   chagb = 0   ! 1 = CHAGB, 0 = GB_canonical 
   depth = 2
   rbornstat = 0
   dgij = 0
   ipb = 2 ! Saeed: to test  
 
   ! initialized the default pb parameters 
   outphi = .false.
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

   !! Physical parameters
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
   radiopt = 0  ! this option is to read radii from topology
   radiscale = 1.0d0
   npopt = 2
   decompopt = 2
   use_rmin = 1
   use_sav = 0
   rhow_effect = 1.129d0
   maxsph   = 400
   maxarcdot= 1500
   maxarc   = 512
   maxtri = 10
   nbuffer = 0
   dprob  = 1.40d0
   sprob  = dprob
   vprob  = 1.30d0
   iprob  = 2.00d0

   radinc = 0.8d0
   expthresh = 0.2d0
   arcres = 0.20d0
   cavity_surften = 0.005d0
   cavity_offset = 0.0d0

   nfocus = 1
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
   fillratio= 1.5d0
   xmin = 0.0d0
   xmax = 0.0d0
   ymin = 0.0d0
   ymax = 0.0d0
   zmin = 0.0d0
   zmax = 0.0d0
   buffer = 16*0.5 ! this needs to be updated for NAB
   saopt = 1 ! this needs to be updated for NAB

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
   !donpsa = .true.
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

   !Focusing with partial 2nd level region
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

   mdin_gb=.false.
   call mynmlsrc('gb',5,ifind)
   if (ifind /= 0) mdin_gb=.true. 

   if ( mdin_gb ) then
      rewind 5
      read(5, nml=gb)
   end if
   
   ! updating genborn module
   gb_B = B
   gb_alpb = alpb
   gb_arad = 30  
   gb_depth = depth 
   gb_Rs = Rs
   gb_tau = tau
   gb_ROH = ROH
   gb_cha = chagb
   gb_dgij = dgij
   !This block was copied from gb_read.F90


   if ( npbverb > 0 ) pbverbose = .true.
!if ( npbverb > 0 .or. pbverbose ) then
!   print *, ' ndofd =', ndofd,'ndosas =',ndosas
!   print *, ' epsin =', epsin,'epsout =',epsout
!   print *, 'istrng =',istrng,' fioni =', fioni
!   print *, 'pbtemp =',pbtemp
!end if

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
      write(6, '(a,f12.4)') 'GB Warning in gb_read(): arcres was too big and is now set to', arcres
   end if

   ! check surface
   if ( ipb == 1 .and. sasopt == 2 ) then
      write(6, '(a)') 'GB Bomb in gb_read(): the smooth Gaussian surface must be built with a level set funtion (ipb > 1)'
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
         write(6, '(a)') 'GB Bomb in gb_read(): ion probe cannot be smaller than solvent prob when the SAS is used (sasopt = 1)'
         call mexit(6, 1)
      else
         write(6, '(a)') 'GB Warning in gb_read(): ion probe is smaller than solvent prob'
      end if
   end if
   if ( dprob <= space .and. ipb == 1 ) then
      write(6, '(a)') 'GB Warning in gb_read(): setting grid spacing larger than solvent probe'
      write(6, '(a)') 'may cause numerical instability if ipb=1'
   end if
   if ( dprob <= space .and. ipb == 2 ) then
      write(6, '(a)') 'GB Bomb in gb_read(): solvent probe cannot be smaller than grid spacing if ipb=2'
      call mexit(6, 1)
   end if



   ! check pb options

   ! Paranoid Check by Mengjuei

   if ( ipb < 1 ) then
      write(6, *) 'ipb =', ipb
      write(6, '(a)') 'GB Bomb in gb_read(): ipb can only be larger than or equal to 1.'
      call mexit(6, 1)
   else if (igb == 10 .and. ipb == 0) then
      ipb=2
   else if (igb /= 0 .and. igb /= 10) then
!     write(6, *) 'Error: igb can only be 10 in GBSA'
!     write(6, *) 'GB Info in gb_read(): igb has been overwritten by ipb'
      igb=10
   else if ( igb == 10 .and. ipb /= 2 ) then
!     write(6, *) 'GB Info in gb_read(): igb has been overwritten by ipb'
   end if

   if ( ipb == 4 ) then
      if ( bcopt < 4 ) then
         bcopt = 6
         write(6,'(a)') 'GB Info in gb_read(): bcopt has been reset to 6 with ipb=4'
      end if
      if ( nfocus /= 1 ) then
         nfocus = 1
         write(6,'(a)') 'GB Info in gb_read(): nfocus has been reset to 1 with ipb=4'
      end if
!     if ( frcopt == 1 ) then 
!        write(6, *) 'GB Info in gb_read(): frcopt=1 conflicts with ipb=4' 
!        call mexit(6, 1)
!     end if
      if ( istrng /= 0 ) then
         write(6,'(a)') 'GB Info in gb_read(): ipb=4 cannot deal with salt solution' 
         call mexit(6, 1)
      end if
   end if

   ! check focusing options

   if ( nfocus > 1 .and. fscale < 2 ) then
      write(6,'(a)') 'GB Info in gb_read(): nfocus reset to 1 due small fscale.'
      nfocus = 1
   endif
   
   ! check nonpolar options

   if ( inp /= npopt ) then
!     write(6, *) 'GB Info in gb_read(): npopt has been overwritten with inp.'
      npopt = inp
   endif
!  if ( inp == 2 .and. radiopt /= 1 ) then
!     write(6, *) 'GB Warning in gb_read(): radiopt has been reset to 1 because inp=2.'
!     radiopt = 1
!  end if
   if ( inp == 1 .and. use_sav == 1 ) then
      write(6, '(a)') &
         'GB Warning in gb_read(): the cavity energy should be correlated to &
         & the SASA when inp=1, so use_sav is reset to 0.'
      use_sav = 0
   end if
   if ( inp /= 2 .and. use_rmin == 1 ) then
!     write(6, *) 'GB Warning in gb_read(): use van der Waals sigma values when inp/=2, &
!                & so use_rmin is reset to 0.'
      use_rmin = 0
   end if
   if ( inp == 1 .and. sprob == 0.557d0 ) then
      write(6, '(a)') &
         'GB Warning in gb_read(): sprob=0.557 is optimized for inp=2 and &
         & should not be used with inp=1. It has been reset to 1.4.'
      sprob = 1.4d0 
   end if
   if ( inp == 1 .and. cavity_surften == 0.0378d0 ) then
      write(6, '(a)') &
         'GB Warning in gb_read(): cavity_surften=0.0378 is optimized for inp=2 &
         & and should not be used with inp=1. It has been reset to 0.04356.'
      cavity_surften =  0.0072d0
   end if
   if ( inp == 1 .and. cavity_offset == -0.5692d0 ) then
      write(6, '(a)') &
         'GB Warning in gb_read(): cavity_offset=-0.5692 is optimized for inp=2 &
         & and should not be used with inp=1. It has been reset to -1.008.'
      cavity_offset =  -1.008d0
   end if
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
         write(6, '(a)') 'GB Info in gb_read(): only dbfopt = 0 or 1 are supported'
         write(6, '(a)') 'GB Info in gb_read(): dbfopt is replaced by eneopt'
         call mexit(6, 1)
      end if
   else
      if ( dbfopt == 0 .and. eneopt == 1 ) then 
!        write(6, *) 'GB Info in gb_read(): dbfopt has been overwritten by eneopt'
      else if ( dbfopt == 1 .and. eneopt == 2 ) then 
!        write(6, *) 'GB Info in gb_read(): dbfopt has been overwritten by eneopt'
      else if ( dbfopt /= -1 ) then
         write(6, *) 'GB Info in gb_read(): dbfopt is ignored when eneopt is set'
      end if
   end if
   if ( eneopt < 1 .or. eneopt > 3 ) then
      write(6, '(a)') 'GB Bomb in gb_read(): only eneopt= 1-3 are supported'
      call mexit(6, 1)
   end if
   if ( frcopt < 0 .or. frcopt > 4 ) then
      write(6, '(a)') 'GB Bomb in gb_read(): only frcopt= 0-4 are supported'
      call mexit(6, 1)
   end if
   if ( eneopt == 1 .and. ( frcopt >= 2 ) ) then
      write(6, '(a)') 'GB Bomb in gb_read(): combination of eneopt and frcopt is not supported', eneopt, frcopt
      call mexit(6, 1)
   end if
   if ( eneopt == 2 .and. frcopt == 1 ) then
!     write(6, *) 'GB Bomb in gb_read(): combination of eneopt and frcopt is not supported', eneopt, frcopt
!     call mexit(6, 1)
   end if
   if ( ( frcopt >= 2 .and. frcopt <= 4 ) .and. bcopt == 5 ) then
      bcopt = 6
      write(6, '(a)') 'GB Info in gb_read(): bcopt has been reset from 5 to 6 for frcopt= 2 - 4'
   end if
   if ( frcopt > 0 .and. smoothopt == 0 ) then
!     smoothopt = 1
!     write(6, *) 'GB Info in gb_read(): smoothopt has been reset to 1 for force compuation'
   end if

   if ( ivcap /= 0 .and. (eneopt /= 1 .or. frcopt /= 0) ) then
!     write(6,*) "cap water should only be run with eneopt=1 and frcopt=0"
!     call mexit(6,1)
   end if

   if ( saopt /= 0 .and. abs(saopt) /= 1 .and. abs(saopt) /= 2 ) then
      saopt = 0
      write(6, '(a)') 'GB Info in gb_read(): saopt has been reset to 0 for surface area compuation'
   end if
    
   ! check numerical solution options

   if ( solvopt == 7 ) then
      continue
   else if ( npbopt == 0 .and. solvopt > 5 .and. (.not. solvopt == 8) ) then
      write(6, '(a)') 'GB Bomb in gb_read(): solvopt>5 cannot be used to solve linear GB equation'
      call mexit(6, 1)
   endif
   if ( ipb /= 4 .and. solvopt == 2 .and. bcopt == 6 ) then
      write(6, '(a)') 'GB Bomb in gb_read(): bcopt=6 cannot be used with solvopt=2'
      call mexit(6, 1)
   end if
   if ( npbopt == 1 .and. bcopt == 10 ) then
      write(6, '(a)') 'GB Bomb in gb_read(): bcopt=10 can be used only with npbopt=0'
      call mexit(6, 1)
   end if
   if ( solvopt /= 3 .and. solvopt /=8 .and. solvopt /=1 .and. bcopt == 10 ) then
      write(6, '(a)') 'GB Bomb in gb_read(): bcopt=10 can be used only with solvopt=1,3, or 8'
      call mexit(6, 1)
   end if
   if ( npbopt == 1 .and. solvopt > 6 ) then
      write(6, '(a)') 'GB Bomb in gb_read(): nonsupported solvopt (>6) for e nonlinear GB equation'
      call mexit(6, 1)
   endif
   if ( npbopt == 1 .and. eneopt == 2 ) then
      eneopt = 1
      write(6, '(a)') 'GB Info in gb_read(): eneopt has been reset to be 1 for nonlinear GB equation'
   end if
   if ( npbopt == 1 .and. frcopt > 0 ) then
      write(6, '(a)') 'GB Bomb in gb_read(): force computatoin is not supported for nonlinear GB equation'
      call mexit(6, 1)
   end if

   ! check cutoff options
    
   if ( cutfd > cutsa ) then
      cutsa = cutfd
   end if
   if ( cutnb /= 0 .and. cutfd > cutnb ) then
      cutnb = cutfd
      write(6, '(a)') 'GB Info in gb_read(): cutnb has been reset to be equal to cutfd'
   end if
   if ( max(cutnb,cutsa,cutfd) > cutres ) then
      write(6, '(a)') 'GB Bomb in gb_read(): cutnb/cutfd must be <= cutres', cutnb, cutfd, cutres
      call mexit(6, 1)
   end if
   if ( cutnb /= 0 .and. bcopt > 5 .and. bcopt < 10 ) then
      cutnb = 0
      write(6, '(a)') 'GB Info in gb_read(): cutnb has been reset to 0 with 5 < bcopt < 10', cutnb, bcopt
   end if
   if ( cutnb == 0 .and. eneopt == 1 .and. bcopt < 6 ) then
      write(6, '(a)') 'GB Bomb in gb_read(): cutnb=0 cannot be used with eneopt=1', cutnb, eneopt
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
         write(6,'(a)') "gb_read: [x,y,z]mblk require non-zero ngrdblk[x,y,z]."
         write(6,'(a)') "gb_read: [x,y,z]mblk are reset to zero."
         xmblk = 1; ymblk = 1; zmblk = 1
      end if 
   else if ( mod(ngrdblkx-1,fscale) /= 0 .or. &
             mod(ngrdblky-1,fscale) /= 0 .or. &
             mod(ngrdblkz-1,fscale) /= 0       ) then
         write(6,'(a)') "gb_read: wrong ngrdblk[x,y,z] - fscale combination"
         call mexit(6,1)
   else if ( nfocus /= 2 ) then
      write(6,'(a)') "gb_read: Multiblock focusing need nfocus = 2."
      call mexit(6,1)
   else if ( frcopt /= 0 ) then
      write(6, '(a)') 'gb_read: Multiblock is incompatible with frcopt /= 0.'
      call mexit(6,1)
   else if ( eneopt == 2 ) then
      write(6, '(a)') 'gb_read: Multiblock is incompatible with eneopt == 2.'
      call mexit(6,1)
   else if ( ivcap /= 0 ) then
      write(6, '(a)') 'gb_read: Multiblock is incompatible with cap water.'
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
         write(6, '(a)') 'gb_read: Ligandmask is not compatible with frcopt /= 0.'
         call mexit(6,1)
      else if ( eneopt == 2 ) then
         write(6, '(a)') 'gb_read: Ligandmask is not compatible with eneopt == 2.'
         call mexit(6,1)
      else if ( nfocus < 2 ) then
         write(6, '(a)') 'gb_read: You cannot use ligandmask without focusing.'
         call mexit(6,1)
      else if ( multiblock ) then
         write(6, '(a)') 'gb_read: You cannot use multiblock with ligandmask.'
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

 ! For CHAGB the probe is reduced by Rs for invariant atom-water dist.
   if ( gb_cha == 1 ) dprob = dprob - gb_Rs 

  ! write(6,*)   "READING GB " , B
  ! if (  (gb_B  0.0d0) .and. ( ntom < )   )
   if (gb_cha == 0) then
   if (  (gb_B > 0.0d0)  .and.  (natom < 50) ) then
     write(6, '(a)') 'GB Warning in  gb_read(): B=0 is recommended for molecules with more than 50 atoms'
   endif 
   endif
end subroutine gb_read


