#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBMD cleanup routine
subroutine pb_free

   use variable_module
   use solvent_accessibility

   implicit none
#  include "md.h"
!#ifdef SANDER
!#  include "../../../src/sander/box.h"
!#else
#  include "box.h"
!#endif

   integer alloc_err(64)

   alloc_err(1:64) = 0

   deallocate(   icrd, stat = alloc_err(1 ) )
   deallocate( grdcrg, stat = alloc_err(2 ) )
   deallocate(qgrdcrg, stat = alloc_err(3 ) )
   deallocate(   acrg, stat = alloc_err(4 ) )
   deallocate(   gcrg, stat = alloc_err(5 ) )
   deallocate(  nshrt, stat = alloc_err(6 ) )
   deallocate(    nex, STAT = alloc_err(7 ) )
   deallocate(    iex, STAT = alloc_err(8 ) )
    
   if ( ifcap == 0 .or. (ifcap >= 3 .and. ifcap <= 5) ) then

   deallocate(   acrd, stat = alloc_err(9 ) )
   deallocate(   gcrd, stat = alloc_err(10) )
    
   deallocate(  mdsig, stat = alloc_err(11) )
   deallocate(   rmin, stat = alloc_err(12) )
   deallocate(   radi, stat = alloc_err(13) )
   deallocate(  radip, stat = alloc_err(14) )
   deallocate( radip2, stat = alloc_err(15) )
   deallocate( radip3, stat = alloc_err(16) )
   deallocate( nzratm, stat = alloc_err(17) )
   deallocate(   nmax, stat = alloc_err(18) )
   deallocate(   nexp, stat = alloc_err(19) )
   deallocate(sumnmax, stat = alloc_err(20) )
   deallocate(sumnexp, stat = alloc_err(21) )
   deallocate( avnmax, stat = alloc_err(22) )
   deallocate( avnexp, stat = alloc_err(23) )
   deallocate(   scrd, stat = alloc_err(24) )
 
   else
 
   deallocate(    acrd, stat = alloc_err(9 ) )
   deallocate(    gcrd, stat = alloc_err(10) )
   deallocate(    radi, stat = alloc_err(11) )
   deallocate(  radip3, stat = alloc_err(12) )
   deallocate(  nzratm, stat = alloc_err(13) )

   end if
    
   deallocate(  iar1pb, stat = alloc_err(25) )
   deallocate( iprshrt, stat = alloc_err(26) )
   !deallocate( iprlong, stat = alloc_err(27) )
   deallocate(   cn1pb, stat = alloc_err(28) )
   deallocate(   cn2pb, stat = alloc_err(29) )
   deallocate(   cn3pb, stat = alloc_err(30) )

   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
        alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
        alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
        alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
        alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30) /= 0 ) then
      write(6, '(a,30i6)') 'PB Bomb in pb_free(): Deallocation aborted', alloc_err(1:30)
      call mexit(6, 1)
   end if

   deallocate(    phi, stat = alloc_err(1 ) )
   deallocate(  chgrd, stat = alloc_err(2 ) )
   if(ipb /= 4 .and. ipb /= 5) deallocate(   epsx, stat = alloc_err(4 ) )
   if(ipb /= 4 .and. ipb /= 5) deallocate(   epsy, stat = alloc_err(5 ) )
   if(ipb /= 4 .and. ipb /= 5) deallocate(   epsz, stat = alloc_err(6 ) )
   if(ipb /= 4 .and. ipb /= 5) deallocate(saltgrd, stat = alloc_err(7 ) )
   deallocate(     bv, stat = alloc_err(3 ) )

   deallocate(  insas, stat = alloc_err(8 ) )
   deallocate( atmsas, stat = alloc_err(9 ) )
   deallocate( lvlset, stat = alloc_err(10) )
   deallocate(     zv, stat = alloc_err(11) )

   deallocate(   cphi, stat = alloc_err(13) )
   if(ipb /= 4 .and. ipb /= 5) deallocate( fedgex, stat = alloc_err(14) )
   if(ipb /= 4 .and. ipb /= 5) deallocate( fedgey, stat = alloc_err(15) )
   if(ipb /= 4 .and. ipb /= 5) deallocate( fedgez, stat = alloc_err(16) )
   deallocate( iepsav, stat = alloc_err(17) )
   if(ipb /= 4 .and. ipb /= 5) deallocate(iepsavx, stat = alloc_err(18) )
   if(ipb /= 4 .and. ipb /= 5) deallocate(iepsavy, stat = alloc_err(19) )
   if(ipb /= 4 .and. ipb /= 5) deallocate(iepsavz, stat = alloc_err(20) )

   deallocate(     xs, stat = alloc_err(30) )
   deallocate(outflag, stat = alloc_err(31) )
   deallocate(outflagorig, stat = alloc_err(32) )
   deallocate( mapout, stat = alloc_err(33) )

   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
        alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
        alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
        alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
        alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30)+&
        alloc_err(31)+alloc_err(32)+alloc_err(33) /= 0 ) then
      write(6,  '(a,33i6)') 'PB Bomb in pb_free(): Deallocation aborted', alloc_err(1:33)
      call mexit(6, 1)
   end if

   if ( ligand .or. multiblock ) then
      deallocate(liveflag, stat = alloc_err(1) )
      deallocate(realflag, stat = alloc_err(2) )
      if ( SUM(alloc_err(1:2)) /= 0 ) then
         write(6,  '(a,2i6)') 'PB Bomb in pb_free(): Deallocation aborted', alloc_err(1:2)
         call mexit(6, 1)
      end if
   endif

end subroutine pb_free
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBMD initialization routine
subroutine pb_init(ifcap,natom,nres,ntypes,nbonh,nbona,ipres,iac,ico,numex,natex,ibh,jbh,iba,jba,ibel,lbres,igraph,isymbl,cg,rin)
     
   ! Module variables
     
   !use genborn ! this to clean GB variables
   use genborn, only : gb_cha, gb_Rs
   use variable_module
   use solvent_accessibility
   use dispersion_cavity

   implicit none
     
   ! Common variables
   ! Currently _REAL_ is not defined
     
!#ifdef SANDER
!#  include "../../../src/sander/md.h"
!#else
#  include "md.h"
!#endif
#  include "parms.h"
#  include "pb_md.h"
     
   ! Passed variables
      
   integer ifcap,natom, nres, ntypes, nbonh, nbona
   integer ipres(*), iac(*), ico(*), numex(*), natex(*), ibh(*), jbh(*), iba(*), jba(*), ibel(*)
   character (len=4) :: lbres(*), igraph(*), isymbl(*)
   _REAL_ cg(natom), rin(natom)
     
   ! Local variables
     
   integer ires, iatm, jatm, maxmax, ic, i, j, jp, idum
   integer alloc_err(64)
   _REAL_ maxnba_l
   _REAL_ rinchk
   _REAL_ ucrgh(natom), ucrga(natom)
   character (len=4) :: residue, resid(natom)
   integer focusresidue(nres)
    
   ! begin of code
     
   alloc_err(1:64) = 0

   ! allocate topology informations
     
   allocate(   icrd( 3,  natom), stat = alloc_err(1 ) )
   allocate( grdcrg( 3,8*natom), stat = alloc_err(2 ) )
   allocate(qgrdcrg(   8*natom), stat = alloc_err(3 ) )
   allocate(   acrg(     natom), stat = alloc_err(4 ) )
   allocate(   gcrg( 8,  natom), stat = alloc_err(5 ) )
   allocate(  nshrt( 0:  natom), stat = alloc_err(6 ) )
   allocate(    nex(     natom), STAT = alloc_err(7 ) )
   allocate(    iex(64,  natom), STAT = alloc_err(8 ) )
    
   if ( ifcap == 0 .or. (ifcap >= 3 .and. ifcap <= 5) ) then

   allocate(   acrd( 3,  natom), stat = alloc_err(9 ) )
   allocate(   gcrd( 3,  natom), stat = alloc_err(10) )
    
   ! allocate sas informations
    
   allocate(  mdsig(  natom  ), stat = alloc_err(11) )
   allocate(   rmin(  natom  ), stat = alloc_err(11) )
   allocate(   radi(  natom  ), stat = alloc_err(12) )
   allocate(  radip(  natom  ), stat = alloc_err(13) )
   allocate( radip2(  natom  ), stat = alloc_err(14) )
   allocate( radip3(  natom  ), stat = alloc_err(15) )
   allocate( nzratm(  natom  ), stat = alloc_err(16) )
   allocate(   nmax(  natom  ), stat = alloc_err(17) )
   allocate(   nexp(  natom  ), stat = alloc_err(18) )
   allocate(sumnmax(  natom  ), stat = alloc_err(19) )
   allocate(sumnexp(  natom  ), stat = alloc_err(20) )
   allocate( avnmax(  natom  ), stat = alloc_err(21) )
   allocate( avnexp(  natom  ), stat = alloc_err(22) )
   allocate(   scrd(3,maxsph ), stat = alloc_err(23) )
 
   else
 
   allocate(   acrd( 3, 0:natom), stat = alloc_err(9 ) )
   allocate(   gcrd( 3, 0:natom), stat = alloc_err(10) )

   allocate(   radi(    0: 0   ), stat = alloc_err(12) )
   allocate( radip3(    1: 1   ), stat = alloc_err(15) )
   allocate( nzratm(    1: 1   ), stat = alloc_err(16) )

   end if
    
   ! allocate pb nblists
    
   maxnba_l = dble(natom) * ( sqrt( max(cutnb,cutsa,cutfd) ) )**3 / 3.0d0
   if ( natom >= 65536 ) then
      write(6,'(a)') "PB Warnning: natom**2 exceeds integer limit (2147483647)."
      maxmax = 2147483647
   else
      maxmax = ceiling(dble(natom)/2*dble(natom))
   end if
   if ( maxnba_l > maxmax ) then
      maxnba = maxmax
   else
      maxnba = int(maxnba_l)
   end if

   allocate( iar1pb (4,0:natom), stat = alloc_err(24) )
   allocate( iprshrt(  maxnba ), stat = alloc_err(25) )
   allocate( cn1pb  (  maxnba ), stat = alloc_err(27) )
   allocate( cn2pb  (  maxnba ), stat = alloc_err(28) )
   allocate( cn3pb  (  maxnba ), stat = alloc_err(29) )

   ! allocate ibelly and outflag

   if ( ifcap == 3 .or. ifcap == 4 ) then
     allocate( ibelly(natom), stat = alloc_err(30) )
   end if
   allocate( outflag(  natom  ), stat = alloc_err(31) )
   allocate( outflagorig(  natom  ), stat = alloc_err(32) )
   allocate(  mapout(  natom  ), stat = alloc_err(33) )

   if ( alloc_err( 1)+alloc_err( 2)+alloc_err( 3)+alloc_err( 4)+alloc_err( 5)+&
        alloc_err( 6)+alloc_err( 7)+alloc_err( 8)+alloc_err( 9)+alloc_err(10)+&
        alloc_err(11)+alloc_err(12)+alloc_err(13)+alloc_err(14)+alloc_err(15)+&
        alloc_err(16)+alloc_err(17)+alloc_err(18)+alloc_err(19)+alloc_err(20)+&
        alloc_err(21)+alloc_err(22)+alloc_err(23)+alloc_err(24)+alloc_err(25)+&
        alloc_err(26)+alloc_err(27)+alloc_err(28)+alloc_err(29)+alloc_err(30)+&
        alloc_err(31)+alloc_err(32)+alloc_err(33) /= 0 ) then
      write(6,  '(a,31i6)') 'PB Bomb in pb_init(): Allocation aborted', alloc_err(1:31)
      call mexit(6, 1)
   end if 
   if ( pbverbose ) then
      write(6,  '()')
      write(6,  '(a)') '======== Implicit Solvent Initialization ========'
      write(6,  '()')
      write(6,'(5x,a,2i9)') 'Max Nonbonded Pairs:', maxnba, maxmax
      write(6,  '()')
   endif

   ! getting some topology info into pb setup ...

   do ires = 1, nres 
      write(residue,'(a4)') lbres(ires)
      do iatm = ipres(ires), ipres(ires+1) - 1
         resid(iatm) = residue
      enddo
   enddo

   ! atomic and total charge in electron for printing only
       
   totcrgp = ZERO
   totcrgn = ZERO
   do iatm = 1, natom
       
      acrg(iatm) = cg(iatm)*INV_AMBER_ELECTROSTATIC
      if ( acrg(iatm) > ZERO) then
         totcrgp = totcrgp + acrg(iatm)
      else
         totcrgn = totcrgn + acrg(iatm)
      end if
       
      ! get info about belly atoms for GBSP

      if(ifcap == 3 .or. ifcap == 4) then
        ibelly(iatm) = ibel(iatm)
      end if

   end do
   totcrg = totcrgp + totcrgn

   ! for pure implicit solvent, set up radii arrays

   if ( ifcap == 0 .or. (ifcap >= 3 .and. ifcap <= 5) ) then

      ! set up group charges for cavity radii analysis

      ucrgh(1:natom) = acrg(1:natom)
      do idum = 1, nbonh
         iatm = ibh(idum)/3 + 1
         jatm = jbh(idum)/3 + 1
         if (isymbl(iatm)(1:1) == 'H' ) then
            ucrgh(jatm) = ucrgh(jatm) + acrg(iatm)
         else
            ucrgh(iatm) = ucrgh(iatm) + acrg(jatm)
         endif
      enddo
                                                                                                          
      ucrga(1:natom) = ucrgh(1:natom)
      do idum = 1, nbona
         iatm = iba(idum)/3 + 1
         jatm = jba(idum)/3 + 1
         ucrga(iatm) = ucrga(iatm) + ucrgh(jatm)
         ucrga(jatm) = ucrga(jatm) + ucrgh(iatm)
      enddo
                                                                                                       
      do idum = 1, natom
         if (isymbl(idum)(1:1) == 'H' ) then
            ucrga(idum) = 0.0d0
            ucrgh(idum) = 0.0d0
         endif
      enddo

      ! van der Waals sigma radii for nonpolar solvation
       
      do iatm = 1, natom
         ic = ico(ntypes*(iac(iatm)-1) + iac(iatm))
         if (cn2(ic) /= ZERO) then
            mdsig(iatm) = (cn1(ic)/cn2(ic))**(SIXTH)/2 ! this is sigma
            rmin(iatm) = mdsig(iatm)*(2.0d0**(SIXTH)) ! this is Rmin 
         else
            mdsig(iatm) = ZERO
            rmin(iatm) = ZERO
         endif
      end do

      ! cavity radii for polar solvation if not passing in from driver
      if (( radiopt == 0 ).or.( radiopt==1 )) then
          if (gb_cha==0) then
              rinchk = ZERO
              do iatm = 1, natom
                rinchk = rinchk+rin(iatm)
              end do
              if (rinchk == ZERO) then
                write(6,'(a)') 'PB Bomb in pb_init(): Requested radi to be read in, but found none'
                call mexit(6,1)
              end if
              radi = rin ! for pb
              mdsig = rin ! for np
          else if  ( gb_cha==1 ) then
              if ( radiopt == 0 ) then 
                 call cha_rad( isymbl, natom, gb_Rs,radi,rin)
                 mdsig = radi
              elseif (radiopt == 1) then
                 rin = rin + gb_Rs
                 mdsig = rin
                 radi = rin
              endif
         endif  
      else
         write(6,'(a,i3)') 'PB Bomb in pb_init(): Unknown radi assigment option', radiopt
         call mexit(6,1)
      end if
      if ( pbverbose ) then
         write(6,*) ' no. of atoms processed in PB initialization:', natom
         write(6, '(a5,3a6,5a10)') 'NUM','RESI','NAME','TYPE',&
         'CHARGE', 'ATM CRG/H', 'GRP CRG', 'PB RADI', 'NP RADI'
         do iatm = 1, natom
            if ( use_rmin == 0 ) write(6, '(i5,3a6,5f10.6)') iatm,resid(iatm),igraph(iatm),isymbl(iatm),&
            real(acrg(iatm)),real(ucrgh(iatm)), real(ucrga(iatm)),&
            real(radi(iatm)), real(mdsig(iatm))
            if ( use_rmin == 1 ) write(6, '(i5,3a6,5f10.6)') iatm,resid(iatm),igraph(iatm),isymbl(iatm),&
            real(acrg(iatm)),real(ucrgh(iatm)), real(ucrga(iatm)),&
            real(radi(iatm)), real(rmin(iatm))
         end do
         write(6,'()')
         write(6,'(a,3f14.4)') ' total system charges (+/-) for PB', totcrg, totcrgp, totcrgn
         write(6,'(a,f14.4,a,f14.4)') ' cavity_surften =', cavity_surften, ' cavity_offset =', cavity_offset
         write(6, '()')
      end if

      ! initialization for sas surface

      if ( srsas ) then

         call sa_sphere(maxsph, scrd)
         if ( pbverbose ) write(6,'(2x,a,i6)') 'SAS Surface: surface dots generated: ', maxsph

         ! nmax and nexp accumulators

         sumnmax(1:natom) = ZERO
         sumnexp(1:natom) = ZERO

      ! initialization for for vdw surface

      else
         if ( pbverbose ) write(6, '(a)') ' VDW Surface: setting up working radii'
         radip3(1:natom) = radi(1:natom)
      end if

   else
      acrd(1:3,0:0) = ZERO; gcrd(1:3,0:0) = ZERO
   end if ! if ( ifcap == 0 .or. (ifcap >= 3 .and. ifcap <= 5) ) then

   ! assinging atom-based pointers to exclusion list

   nshrt(0) = 0
   do i = 1, natom
      nshrt(i) = nshrt(i-1) + numex(i)
   end do
   nex = 0
   do i = 1, natom-1
      do jp = nshrt(i-1) + 1, nshrt(i)
         j = natex(jp)
         if (j == 0)  cycle
         nex(i) = nex(i) + 1
         iex(nex(i),i) = j
         nex(j) = nex(j) + 1
         iex(nex(j),j) = i
      end do
   end do

   ! set up green 3-d array for pb_fdcoulomb

   call pb_green

   ! set up variables if using ligand or multiblock
   if ( ligand .or. multiblock ) then

      allocate( liveflag(natom), stat = alloc_err(1) )
      allocate( realflag(natom), stat = alloc_err(2) )
      if ( SUM(alloc_err(1:2)) /= 0 ) then
         write(6, '(a,2i6)') 'PB Bomb in pb_init(): Allocation aborted', alloc_err(1:2)
         call mexit(6,1)
      end if
     
      ! Mengjuei -
      ! Basically liveflag HERE is a temporary switch storing the list of
      ! atoms labelled ligandmask. However if the geormetries of ligand
      ! box were specified, we don't need liveflag that early.
      ! Later liveflag will become the list of atom inside the block of
      ! rectangular bounding box.
     
      if ( multiblock ) then
         continue!noop
      else if ( xmax - xmin > ZERO .and. &
                ymax - ymin > ZERO .and. &
                zmax - zmin > ZERO       ) then
         ! If the geometries of bounding box are correctly specified,
         ! skip the ligandmask
         continue
      else if ( ligand .and. len_trim(ligandmask) > 0 ) then
         focusresidue = 0
             liveflag = 0
         call myresmask ( ligandmask, LEN_TRIM(ligandmask), focusresidue, nres )
         j = 0
         do i = 1, nres
            if ( focusresidue(i) == 0 ) cycle
            j = j + 1
            liveflag(ipres(i):ipres(i+1)-1)=1
         enddo
         write(6,'(a,a,a,i5)') ' Focusing Mask ', &
              ligandmask(1:LEN_TRIM(ligandmask)), ' matches residue', j
         xmax = ZERO; ymax = ZERO; zmax = ZERO
         xmin = ZERO; ymin = ZERO; zmin = ZERO
      else
         write(6,'(a)') 'wrong ligand box'
         call mexit(6,1)
      endif

   endif

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Setup table of precomputed FD Green's function.
subroutine pb_green

   implicit none

   ! Common variables

   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
   
   _REAL_ rdsqr(0:20, 0:20, 0:20)
   common /blk_rdsqr/ rdsqr
   
   ! Local variables
   
   integer i, j, k
   _REAL_ green_data

   do i = 0, 20
      do j = 0, i
         do k = 0, j
            green_data = green(i, j, k)
            green(i, k, j) = green_data
            green(j, i, k) = green_data
            green(k, i, j) = green_data
            green(j, k, i) = green_data
            green(k, j, i) = green_data
         end do
      end do
   end do
   do i = 0, 20
      do j = 0, 20
         do k = 0, 20
            if ( i == 0 .and. j == 0 .and. k == 0 ) then
               green_data = 0.0d0
            else
               green_data = 1.0d0/sqrt( dble(i**2 + j**2 + k**2) )
            end if
            rdsqr(i, j, k) = green_data
         end do
      end do
   end do

end subroutine pb_green

end subroutine pb_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Data for precomputed FD Green's function
block data blkgreen

   implicit none
   
   ! Common variables
   
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
   
   data green( 0, 0, 0) / 0.3175911535616087d+01 /

   data green( 1, 0, 0) / 0.1081516433229559d+01 /
   data green( 1, 1, 0) / 0.6935560104573442d+00 /
   data green( 1, 1, 1) / 0.5476217521280413d+00 /

   data green( 2, 0, 0) / 0.5389630219348477d+00 /
   data green( 2, 1, 0) / 0.4515298460159590d+00 /
   data green( 2, 1, 1) / 0.4016874937936074d+00 /
   data green( 2, 2, 0) / 0.3525516965096933d+00 /
   data green( 2, 2, 1) / 0.3287203646281152d+00 /
   data green( 2, 2, 2) / 0.2846454604239108d+00 /

   data green( 3, 0, 0) / 0.3461423143146174d+00 /
   data green( 3, 1, 0) / 0.3207333596066411d+00 /
   data green( 3, 1, 1) / 0.3020027893481455d+00 /
   data green( 3, 2, 0) / 0.2774048788858750d+00 /
   data green( 3, 2, 1) / 0.2658750216219361d+00 /
   data green( 3, 2, 2) / 0.2405705562200611d+00 /
   data green( 3, 3, 0) / 0.2351770067410159d+00 /
   data green( 3, 3, 1) / 0.2283341440179216d+00 /
   data green( 3, 3, 2) / 0.2117699730470565d+00 /
   data green( 3, 3, 3) / 0.1912514753585858d+00 /

   data green( 4, 0, 0) / 0.2549574255267407d+00 /
   data green( 4, 1, 0) / 0.2453175397274443d+00 /
   data green( 4, 1, 1) / 0.2371124798390072d+00 /
   data green( 4, 2, 0) / 0.2242171672109520d+00 /
   data green( 4, 2, 1) / 0.2182173966320647d+00 /
   data green( 4, 2, 2) / 0.2034878875608144d+00 /
   data green( 4, 3, 0) / 0.1997919973195232d+00 /
   data green( 4, 3, 1) / 0.1956539205386707d+00 /
   data green( 4, 3, 2) / 0.1849465532311655d+00 /
   data green( 4, 3, 3) / 0.1707329776707052d+00 /
   data green( 4, 4, 0) / 0.1764982221836034d+00 /
   data green( 4, 4, 1) / 0.1736810643669396d+00 /
   data green( 4, 4, 2) / 0.1661064732588335d+00 /
   data green( 4, 4, 3) / 0.1555791639459903d+00 /
   data green( 4, 4, 4) / 0.1438328713882703d+00 /

   data green( 5, 0, 0) / 0.2023320789662321d+00 /
   data green( 5, 1, 0) / 0.1977723261939040d+00 /
   data green( 5, 1, 1) / 0.1936022159068010d+00 /
   data green( 5, 2, 0) / 0.1863537924857621d+00 /
   data green( 5, 2, 1) / 0.1829579023549173d+00 /
   data green( 5, 2, 2) / 0.1740288693193717d+00 /
   data green( 5, 3, 0) / 0.1715517461569988d+00 /
   data green( 5, 3, 1) / 0.1689523673466971d+00 /
   data green( 5, 3, 2) / 0.1619280863139817d+00 /
   data green( 5, 3, 3) / 0.1520949562005345d+00 /
   data green( 5, 4, 0) / 0.1560216045727995d+00 /
   data green( 5, 4, 1) / 0.1540869242343024d+00 /
   data green( 5, 4, 2) / 0.1487427520543840d+00 /
   data green( 5, 4, 3) / 0.1410348417220078d+00 /
   data green( 5, 4, 4) / 0.1320865779400301d+00 /
   data green( 5, 5, 0) / 0.1412654538129695d+00 /
   data green( 5, 5, 1) / 0.1398369570685601d+00 /
   data green( 5, 5, 2) / 0.1358239205529447d+00 /
   data green( 5, 5, 3) / 0.1298867252110265d+00 /
   data green( 5, 5, 4) / 0.1227931068178285d+00 /
   data green( 5, 5, 5) / 0.1152121187708051d+00 /

   data green( 6, 0, 0) / 0.1679457485086538d+00 /
   data green( 6, 1, 0) / 0.1654261095272136d+00 /
   data green( 6, 1, 1) / 0.1630403628732706d+00 /
   data green( 6, 2, 0) / 0.1586657175626259d+00 /
   data green( 6, 2, 1) / 0.1565927725874027d+00 /
   data green( 6, 2, 2) / 0.1509133479959037d+00 /
   data green( 6, 3, 0) / 0.1492383485996070d+00 /
   data green( 6, 3, 1) / 0.1475356222270095d+00 /
   data green( 6, 3, 2) / 0.1428030248125766d+00 /
   data green( 6, 3, 3) / 0.1359109018716310d+00 /
   data green( 6, 4, 0) / 0.1386403575195265d+00 /
   data green( 6, 4, 1) / 0.1372868018013602d+00 /
   data green( 6, 4, 2) / 0.1334762662655059d+00 /
   data green( 6, 4, 3) / 0.1278188743773593d+00 /
   data green( 6, 4, 4) / 0.1210307023559824d+00 /
   data green( 6, 5, 0) / 0.1279377996163346d+00 /
   data green( 6, 5, 1) / 0.1268792625248693d+00 /
   data green( 6, 5, 2) / 0.1238671701342395d+00 /
   data green( 6, 5, 3) / 0.1193168193951752d+00 /
   data green( 6, 5, 4) / 0.1137433220221452d+00 /
   data green( 6, 5, 5) / 0.1076311304178983d+00 /
   data green( 6, 6, 0) / 0.1177570181416862d+00 /
   data green( 6, 6, 1) / 0.1169336535622994d+00 /
   data green( 6, 6, 2) / 0.1145701568434098d+00 /
   data green( 6, 6, 3) / 0.1109457930352897d+00 /
   data green( 6, 6, 4) / 0.1064235923531429d+00 /
   data green( 6, 6, 5) / 0.1013649659819519d+00 /
   data green( 6, 6, 6) / 0.9607596998651532d-01 /

   data green( 7, 0, 0) / 0.1436379703676282d+00 /
   data green( 7, 1, 0) / 0.1420921430389773d+00 /
   data green( 7, 1, 1) / 0.1406021945384286d+00 /
   data green( 7, 2, 0) / 0.1377905085064394d+00 /
   data green( 7, 2, 1) / 0.1364436836006640d+00 /
   data green( 7, 2, 2) / 0.1326596236361408d+00 /
   data green( 7, 3, 0) / 0.1315010244137252d+00 /
   data green( 7, 3, 1) / 0.1303404244861499d+00 /
   data green( 7, 3, 2) / 0.1270539222330616d+00 /
   data green( 7, 3, 3) / 0.1221266551245613d+00 /
   data green( 7, 4, 0) / 0.1240707811356643d+00 /
   data green( 7, 4, 1) / 0.1231023802935189d+00 /
   data green( 7, 4, 2) / 0.1203389703389623d+00 /
   data green( 7, 4, 3) / 0.1161437166739061d+00 /
   data green( 7, 4, 4) / 0.1109732503425302d+00 /
   data green( 7, 5, 0) / 0.1162054481591678d+00 /
   data green( 7, 5, 1) / 0.1154131802988344d+00 /
   data green( 7, 5, 2) / 0.1131365927016754d+00 /
   data green( 7, 5, 3) / 0.1096390267904409d+00 /
   data green( 7, 5, 4) / 0.1052645749390114d+00 /
   data green( 7, 5, 5) / 0.1003580913834594d+00 /
   data green( 7, 6, 0) / 0.1083996007720677d+00 /
   data green( 7, 6, 1) / 0.1077581123377295d+00 /
   data green( 7, 6, 2) / 0.1059035863647062d+00 /
   data green( 7, 6, 3) / 0.1030236836257126d+00 /
   data green( 7, 6, 4) / 0.9937207589038648d-01 /
   data green( 7, 6, 5) / 0.9521398174874596d-01 /
   data green( 7, 6, 6) / 0.9078697558798884d-01 /
   data green( 7, 7, 0) / 0.1009545776025604d+00 /
   data green( 7, 7, 1) / 0.1004370965875413d+00 /
   data green( 7, 7, 2) / 0.9893330927694677d-01 /
   data green( 7, 7, 3) / 0.9657601000730712d-01 /
   data green( 7, 7, 4) / 0.9355002327258850d-01 /
   data green( 7, 7, 5) / 0.9005596172135705d-01 /
   data green( 7, 7, 6) / 0.8628155587586626d-01 /
   data green( 7, 7, 7) / 0.8238481603995024d-01 /

   data green( 8, 0, 0) / 0.1255135043441952d+00 /
   data green( 8, 1, 0) / 0.1244938756234683d+00 /
   data green( 8, 1, 1) / 0.1235011515169599d+00 /
   data green( 8, 2, 0) / 0.1215968001009854d+00 /
   data green( 8, 2, 1) / 0.1206765773688242d+00 /
   data green( 8, 2, 2) / 0.1180491816994469d+00 /
   data green( 8, 3, 0) / 0.1172256514690563d+00 /
   data green( 8, 3, 1) / 0.1164059153275179d+00 /
   data green( 8, 3, 2) / 0.1140548340090682d+00 /
   data green( 8, 3, 3) / 0.1104537583185242d+00 /
   data green( 8, 4, 0) / 0.1118731000057510d+00 /
   data green( 8, 4, 1) / 0.1111641239933093d+00 /
   data green( 8, 4, 2) / 0.1091209423292397d+00 /
   data green( 8, 4, 3) / 0.1059655290946281d+00 /
   data green( 8, 4, 4) / 0.1019922051725589d+00 /
   data green( 8, 5, 0) / 0.1059981467610767d+00 /
   data green( 8, 5, 1) / 0.1053972859411583d+00 /
   data green( 8, 5, 2) / 0.1036576323080125d+00 /
   data green( 8, 5, 3) / 0.1009487738681923d+00 /
   data green( 8, 5, 4) / 0.9750168346024086d-01 /
   data green( 8, 5, 5) / 0.9356030135870605d-01 /
   data green( 8, 6, 0) / 0.9996433592021670d-01 /
   data green( 8, 6, 1) / 0.9946155575783373d-01 /
   data green( 8, 6, 2) / 0.9799966070635584d-01 /
   data green( 8, 6, 3) / 0.9570560330503135d-01 /
   data green( 8, 6, 4) / 0.9275661615274926d-01 /
   data green( 8, 6, 5) / 0.8934581858423703d-01 /
   data green( 8, 6, 6) / 0.8565480527199545d-01 /
   data green( 8, 7, 0) / 0.9402703640079557d-01 /
   data green( 8, 7, 1) / 0.9360923311145095d-01 /
   data green( 8, 7, 2) / 0.9238978732931179d-01 /
   data green( 8, 7, 3) / 0.9046268273899294d-01 /
   data green( 8, 7, 4) / 0.8796200130961226d-01 /
   data green( 8, 7, 5) / 0.8503812765429647d-01 /
   data green( 8, 7, 6) / 0.8183730074804783d-01 /
   data green( 8, 7, 7) / 0.7848807686205882d-01 /
   data green( 8, 8, 0) / 0.8834709960853802d-01 /
   data green( 8, 8, 1) / 0.8800081609344552d-01 /
   data green( 8, 8, 2) / 0.8698673290013459d-01 /
   data green( 8, 8, 3) / 0.8537410642727687d-01 /
   data green( 8, 8, 4) / 0.8326364507844906d-01 /
   data green( 8, 8, 5) / 0.8077123734867642d-01 /
   data green( 8, 8, 6) / 0.7801307831495344d-01 /
   data green( 8, 8, 7) / 0.7509490959573568d-01 /
   data green( 8, 8, 8) / 0.7210598883606888d-01 /

   data green( 9, 0, 0) / 0.1114675533037391d+00 /
   data green( 9, 1, 0) / 0.1107585061476348d+00 /
   data green( 9, 1, 1) / 0.1100638094274859d+00 /
   data green( 9, 2, 0) / 0.1087176090578064d+00 /
   data green( 9, 2, 1) / 0.1080627359638833d+00 /
   data green( 9, 2, 2) / 0.1061726441411663d+00 /
   data green( 9, 3, 0) / 0.1055711554782572d+00 /
   data green( 9, 3, 1) / 0.1049738798870201d+00 /
   data green( 9, 3, 2) / 0.1032452821957364d+00 /
   data green( 9, 3, 3) / 0.1005551629045936d+00 /
   data green( 9, 4, 0) / 0.1016157734602918d+00 /
   data green( 9, 4, 1) / 0.1010851206318975d+00 /
   data green( 9, 4, 2) / 0.9954456947893656d-01 /
   data green( 9, 4, 3) / 0.9713375096738032d-01 /
   data green( 9, 4, 4) / 0.9404557003280233d-01 /
   data green( 9, 5, 0) / 0.9715142575127113d-01 /
   data green( 9, 5, 1) / 0.9668907145426274d-01 /
   data green( 9, 5, 2) / 0.9534254334095228d-01 /
   data green( 9, 5, 3) / 0.9322316123408030d-01 /
   data green( 9, 5, 4) / 0.9048763199919399d-01 /
   data green( 9, 5, 5) / 0.8730870895555073d-01 /
   data green( 9, 6, 0) / 0.9243812005955071d-01 /
   data green( 9, 6, 1) / 0.9204070685710865d-01 /
   data green( 9, 6, 2) / 0.9087979193006089d-01 /
   data green( 9, 6, 3) / 0.8904223038675743d-01 /
   data green( 9, 6, 4) / 0.8665250254084390d-01 /
   data green( 9, 6, 5) / 0.8385108509327174d-01 /
   data green( 9, 6, 6) / 0.8077561652110507d-01 /
   data green( 9, 7, 0) / 0.8767773610633496d-01 /
   data green( 9, 7, 1) / 0.8733910756827704d-01 /
   data green( 9, 7, 2) / 0.8634710578292538d-01 /
   data green( 9, 7, 3) / 0.8476858114697289d-01 /
   data green( 9, 7, 4) / 0.8270091342713728d-01 /
   data green( 9, 7, 5) / 0.8025644363570066d-01 /
   data green( 9, 7, 6) / 0.7754815339335515d-01 /
   data green( 9, 7, 7) / 0.7467922519730323d-01 /
   data green( 9, 8, 0) / 0.8301344674702572d-01 /
   data green( 9, 8, 1) / 0.8272629913286701d-01 /
   data green( 9, 8, 2) / 0.8188294949837628d-01 /
   data green( 9, 8, 3) / 0.8053444819367288d-01 /
   data green( 9, 8, 4) / 0.7875626139571232d-01 /
   data green( 9, 8, 5) / 0.7663722110817615d-01 /
   data green( 9, 8, 6) / 0.7426886866755487d-01 /
   data green( 9, 8, 7) / 0.7173711721029552d-01 /
   data green( 9, 8, 8) / 0.6911706841658720d-01 /
   data green( 9, 9, 0) / 0.7853816503164655d-01 /
   data green( 9, 9, 1) / 0.7829513028986744d-01 /
   data green( 9, 9, 2) / 0.7757970492213170d-01 /
   data green( 9, 9, 3) / 0.7643074722956934d-01 /
   data green( 9, 9, 4) / 0.7490645943831237d-01 /
   data green( 9, 9, 5) / 0.7307660435016609d-01 /
   data green( 9, 9, 6) / 0.7101468041178351d-01 /
   data green( 9, 9, 7) / 0.6879146254640564d-01 /
   data green( 9, 9, 8) / 0.6647054573580437d-01 /
   data green( 9, 9, 9) / 0.6410594713334371d-01 /

   data green(10, 0, 0) / 0.1002577912116933d+00 /
   data green(10, 1, 0) / 0.9974437853334095d-01 /
   data green(10, 1, 1) / 0.9923922073815769d-01 /
   data green(10, 2, 0) / 0.9825372155413650d-01 /
   data green(10, 2, 1) / 0.9777189449728949d-01 /
   data green(10, 2, 2) / 0.9637064580629341d-01 /
   data green(10, 3, 0) / 0.9592013925418678d-01 /
   data green(10, 3, 1) / 0.9547306881115049d-01 /
   data green(10, 3, 2) / 0.9417060949950148d-01 /
   data green(10, 3, 3) / 0.9211915732388094d-01 /
   data green(10, 4, 0) / 0.9292870877604903d-01 /
   data green(10, 4, 1) / 0.9252330606117132d-01 /
   data green(10, 4, 2) / 0.9133977791419053d-01 /
   data green(10, 4, 3) / 0.8946851694471304d-01 /
   data green(10, 4, 4) / 0.8703843787554651d-01 /
   data green(10, 5, 0) / 0.8947837360825321d-01 /
   data green(10, 5, 1) / 0.8911734567323766d-01 /
   data green(10, 5, 2) / 0.8806103347227387d-01 /
   data green(10, 5, 3) / 0.8638403618452394d-01 /
   data green(10, 5, 4) / 0.8419416553647674d-01 /
   data green(10, 5, 5) / 0.8161451740901197d-01 /
   data green(10, 6, 0) / 0.8575380141340755d-01 /
   data green(10, 6, 1) / 0.8543660996393478d-01 /
   data green(10, 6, 2) / 0.8450650433993552d-01 /
   data green(10, 6, 3) / 0.8302373791444409d-01 /
   data green(10, 6, 4) / 0.8107653910052529d-01 /
   data green(10, 6, 5) / 0.7876742203076122d-01 /
   data green(10, 6, 6) / 0.7620040733214006d-01 /
   data green(10, 7, 0) / 0.8190959917492843d-01 /
   data green(10, 7, 1) / 0.8163356160340368d-01 /
   data green(10, 7, 2) / 0.8082241804381261d-01 /
   data green(10, 7, 3) / 0.7952410551005101d-01 /
   data green(10, 7, 4) / 0.7780969180796847d-01 /
   data green(10, 7, 5) / 0.7576315982585451d-01 /
   data green(10, 7, 6) / 0.7347146592495993d-01 /
   data green(10, 7, 7) / 0.7101674376460161d-01 /
   data green(10, 8, 0) / 0.7806508211085719d-01 /
   data green(10, 8, 1) / 0.7782634159391273d-01 /
   data green(10, 8, 2) / 0.7712340542901212d-01 /
   data green(10, 8, 3) / 0.7599404356180570d-01 /
   data green(10, 8, 4) / 0.7449488865948030d-01 /
   data green(10, 8, 5) / 0.7269391793199149d-01 /
   data green(10, 8, 6) / 0.7066292377435307d-01 /
   data green(10, 8, 7) / 0.6847117504139806d-01 /
   data green(10, 8, 8) / 0.6618109564712377d-01 /
   data green(10, 9, 0) / 0.7430591829056880d-01 /
   data green(10, 9, 1) / 0.7410015733590297d-01 /
   data green(10, 9, 2) / 0.7349322668675839d-01 /
   data green(10, 9, 3) / 0.7251471189729476d-01 /
   data green(10, 9, 4) / 0.7120945300914271d-01 /
   data green(10, 9, 5) / 0.6963202182589029d-01 /
   data green(10, 9, 6) / 0.6784114522187117d-01 /
   data green(10, 9, 7) / 0.6589465567126211d-01 /
   data green(10, 9, 8) / 0.6384585224195553d-01 /
   data green(10, 9, 9) / 0.6174134887775983d-01 /
   data green(10,10, 0) / 0.7068920367524822d-01 /
   data green(10,10, 1) / 0.7051211509511597d-01 /
   data green(10,10, 2) / 0.6998890177534903d-01 /
   data green(10,10, 3) / 0.6914268527305115d-01 /
   data green(10,10, 4) / 0.6800883100276685d-01 /
   data green(10,10, 5) / 0.6663099443635093d-01 /
   data green(10,10, 6) / 0.6505687404676844d-01 /
   data green(10,10, 7) / 0.6333437786552837d-01 /
   data green(10,10, 8) / 0.6150862452163437d-01 /
   data green(10,10, 9) / 0.5961995993771525d-01 /
   data green(10,10,10) / 0.5770290160588142d-01 /

   data green(11, 0, 0) / 0.9110168046434086d-01 /
   data green(11, 1, 0) / 0.9071780901336175d-01 /
   data green(11, 1, 1) / 0.9033897123132434d-01 /
   data green(11, 2, 0) / 0.8959641334512164d-01 /
   data green(11, 2, 1) / 0.8923197433615886d-01 /
   data green(11, 2, 2) / 0.8816622137768942d-01 /
   data green(11, 3, 0) / 0.8782111225409407d-01 /
   data green(11, 3, 1) / 0.8747858387098068d-01 /
   data green(11, 3, 2) / 0.8647572516935614d-01 /
   data green(11, 3, 3) / 0.8488152431551316d-01 /
   data green(11, 4, 0) / 0.8551135644152588d-01 /
   data green(11, 4, 1) / 0.8519581285733283d-01 /
   data green(11, 4, 2) / 0.8427063105410360d-01 /
   data green(11, 4, 3) / 0.8279594800076019d-01 /
   data green(11, 4, 4) / 0.8085968878575372d-01 /
   data green(11, 5, 0) / 0.8280161135253408d-01 /
   data green(11, 5, 1) / 0.8251568011609525d-01 /
   data green(11, 5, 2) / 0.8167599625304525d-01 /
   data green(11, 5, 3) / 0.8033360087301226d-01 /
   data green(11, 5, 4) / 0.7856390481651128d-01 /
   data green(11, 5, 5) / 0.7645522053526901d-01 /
   data green(11, 6, 0) / 0.7982349651519470d-01 /
   data green(11, 6, 1) / 0.7956773670882406d-01 /
   data green(11, 6, 2) / 0.7881543356002289d-01 /
   data green(11, 6, 3) / 0.7760901231430582d-01 /
   data green(11, 6, 4) / 0.7601171240344266d-01 /
   data green(11, 6, 5) / 0.7409882242897636d-01 /
   data green(11, 6, 6) / 0.7194903612834178d-01 /
   data green(11, 7, 0) / 0.7669385142896425d-01 /
   data green(11, 7, 1) / 0.7646729557509413d-01 /
   data green(11, 7, 2) / 0.7579982415468997d-01 /
   data green(11, 7, 3) / 0.7472616108426046d-01 /
   data green(11, 7, 4) / 0.7329854895221447d-01 /
   data green(11, 7, 5) / 0.7158000929949575d-01 /
   data green(11, 7, 6) / 0.6963746184484755d-01 /
   data green(11, 7, 7) / 0.6753593972476768d-01 /
   data green(11, 8, 0) / 0.7350884699773956d-01 /
   data green(11, 8, 1) / 0.7330954273829536d-01 /
   data green(11, 8, 2) / 0.7272142607912221d-01 /
   data green(11, 8, 3) / 0.7177270358992764d-01 /
   data green(11, 8, 4) / 0.7050594767835185d-01 /
   data green(11, 8, 5) / 0.6897329874034243d-01 /
   data green(11, 8, 6) / 0.6723095574112561d-01 /
   data green(11, 8, 7) / 0.6533450739352768d-01 /
   data green(11, 8, 8) / 0.6333543002203640d-01 /
   data green(11, 9, 0) / 0.7034274444422710d-01 /
   data green(11, 9, 1) / 0.7016821199487679d-01 /
   data green(11, 9, 2) / 0.6965247774671877d-01 /
   data green(11, 9, 3) / 0.6881812686676836d-01 /
   data green(11, 9, 4) / 0.6769973408104581d-01 /
   data green(11, 9, 5) / 0.6634002301146850d-01 /
   data green(11, 9, 6) / 0.6478574078950290d-01 /
   data green(11, 9, 7) / 0.6308391147473533d-01 /
   data green(11, 9, 8) / 0.6127891097843042d-01 /
   data green(11, 9, 9) / 0.5941049821204494d-01 /
   data green(11,10, 0) / 0.6724957721046959d-01 /
   data green(11,10, 1) / 0.6709713546389921d-01 /
   data green(11,10, 2) / 0.6664607738923635d-01 /
   data green(11,10, 3) / 0.6591447785283995d-01 /
   data green(11,10, 4) / 0.6493021401024528d-01 /
   data green(11,10, 5) / 0.6372810888526606d-01 /
   data green(11,10, 6) / 0.6234678864803328d-01 /
   data green(11,10, 7) / 0.6082572759401585d-01 /
   data green(11,10, 8) / 0.5920284093129349d-01 /
   data green(11,10, 9) / 0.5751276763382698d-01 /
   data green(11,10,10) / 0.5578584437122584d-01 /
   data green(11,11, 0) / 0.6426621872093206d-01 /
   data green(11,11, 1) / 0.6413321438931355d-01 /
   data green(11,11, 2) / 0.6373918858005193d-01 /
   data green(11,11, 3) / 0.6309858023634662d-01 /
   data green(11,11, 4) / 0.6223381115297828d-01 /
   data green(11,11, 5) / 0.6117317353634055d-01 /
   data green(11,11, 6) / 0.5994843804384092d-01 /
   data green(11,11, 7) / 0.5859255646394460d-01 /
   data green(11,11, 8) / 0.5713771112635687d-01 /
   data green(11,11, 9) / 0.5561383423715520d-01 /
   data green(11,11,10) / 0.5404766636162724d-01 /
   data green(11,11,11) / 0.5246225811515642d-01 /

   data green(12, 0, 0) / 0.8348105551187303d-01 /
   data green(12, 1, 0) / 0.8318644007573157d-01 /
   data green(12, 1, 1) / 0.8289503946321708d-01 /
   data green(12, 2, 0) / 0.8232188867859169d-01 /
   data green(12, 2, 1) / 0.8203976357893702d-01 /
   data green(12, 2, 2) / 0.8121127206628913d-01 /
   data green(12, 3, 0) / 0.8094159670018782d-01 /
   data green(12, 3, 1) / 0.8067380953205813d-01 /
   data green(12, 3, 2) / 0.7988678238266739d-01 /
   data green(12, 3, 3) / 0.7862664221116486d-01 /
   data green(12, 4, 0) / 0.7912508082373361d-01 /
   data green(12, 4, 1) / 0.7887532038124223d-01 /
   data green(12, 4, 2) / 0.7814052333646752d-01 /
   data green(12, 4, 3) / 0.7696172849284801d-01 /
   data green(12, 4, 4) / 0.7540012125756697d-01 /
   data green(12, 5, 0) / 0.7696508124722977d-01 /
   data green(12, 5, 1) / 0.7673557789782000d-01 /
   data green(12, 5, 2) / 0.7605959809871861d-01 /
   data green(12, 5, 3) / 0.7497277718386208d-01 /
   data green(12, 5, 4) / 0.7352861983168400d-01 /
   data green(12, 5, 5) / 0.7179149079871779d-01 /
   data green(12, 6, 0) / 0.7455624179495302d-01 /
   data green(12, 6, 1) / 0.7434790514590302d-01 /
   data green(12, 6, 2) / 0.7373352172905660d-01 /
   data green(12, 6, 3) / 0.7274343090725080d-01 /
   data green(12, 6, 4) / 0.7142351511344813d-01 /
   data green(12, 6, 5) / 0.6982953516262391d-01 /
   data green(12, 6, 6) / 0.6802124435664933d-01 /
   data green(12, 7, 0) / 0.7198657625585808d-01 /
   data green(12, 7, 1) / 0.7179925722029065d-01 /
   data green(12, 7, 2) / 0.7124617598124484d-01 /
   data green(12, 7, 3) / 0.7035276687354203d-01 /
   data green(12, 7, 4) / 0.6915777579243522d-01 /
   data green(12, 7, 5) / 0.6770876518667349d-01 /
   data green(12, 7, 6) / 0.6605734858473514d-01 /
   data green(12, 7, 7) / 0.6425496074263920d-01 /
   data green(12, 8, 0) / 0.6933230991524528d-01 /
   data green(12, 8, 1) / 0.6916510150934985d-01 /
   data green(12, 8, 2) / 0.6867080424528094d-01 /
   data green(12, 8, 3) / 0.6787049670629464d-01 /
   data green(12, 8, 4) / 0.6679651152816783d-01 /
   data green(12, 8, 5) / 0.6548893633145812d-01 /
   data green(12, 8, 6) / 0.6399180585521055d-01 /
   data green(12, 8, 7) / 0.6234963341352993d-01 /
   data green(12, 8, 8) / 0.6060464826331222d-01 /
   data green(12, 9, 0) / 0.6665570130238527d-01 /
   data green(12, 9, 1) / 0.6650721515618704d-01 /
   data green(12, 9, 2) / 0.6606775879817525d-01 /
   data green(12, 9, 3) / 0.6535470380718020d-01 /
   data green(12, 9, 4) / 0.6439465414595712d-01 /
   data green(12, 9, 5) / 0.6322123372119758d-01 /
   data green(12, 9, 6) / 0.6187161185218396d-01 /
   data green(12, 9, 7) / 0.6038392502567606d-01 /
   data green(12, 9, 8) / 0.5879492104517095d-01 /
   data green(12, 9, 9) / 0.5713828397746269d-01 /
   data green(12,10, 0) / 0.6400502526129169d-01 /
   data green(12,10, 1) / 0.6387361729901000d-01 /
   data green(12,10, 2) / 0.6348428440102762d-01 /
   data green(12,10, 3) / 0.6285118378330068d-01 /
   data green(12,10, 4) / 0.6199631718152376d-01 /
   data green(12,10, 5) / 0.6094746014574964d-01 /
   data green(12,10, 6) / 0.5973584209308586d-01 /
   data green(12,10, 7) / 0.5839388844520180d-01 /
   data green(12,10, 8) / 0.5695330474238792d-01 /
   data green(12,10, 9) / 0.5544362901680640d-01 /
   data green(12,10,10) / 0.5389128954425403d-01 /
   data green(12,11, 0) / 0.6141586482286285d-01 /
   data green(12,11, 1) / 0.6129980328751115d-01 /
   data green(12,11, 2) / 0.6095559091004567d-01 /
   data green(12,11, 3) / 0.6039476279480380d-01 /
   data green(12,11, 4) / 0.5963534456532297d-01 /
   data green(12,11, 5) / 0.5870028682792838d-01 /
   data green(12,11, 6) / 0.5761565979527127d-01 /
   data green(12,11, 7) / 0.5640887125213167d-01 /
   data green(12,11, 8) / 0.5510709620621658d-01 /
   data green(12,11, 9) / 0.5373604601398634d-01 /
   data green(12,11,10) / 0.5231910960458832d-01 /
   data green(12,11,11) / 0.5087684983588978d-01 /
   data green(12,12, 0) / 0.5891302695087103d-01 /
   data green(12,12, 1) / 0.5881060485306964d-01 /
   data green(12,12, 2) / 0.5850656139122329d-01 /
   data green(12,12, 3) / 0.5801028134948139d-01 /
   data green(12,12, 4) / 0.5733650734563099d-01 /
   data green(12,12, 5) / 0.5650415323776656d-01 /
   data green(12,12, 6) / 0.5553491399554743d-01 /
   data green(12,12, 7) / 0.5445186071589581d-01 /
   data green(12,12, 8) / 0.5327817184342189d-01 /
   data green(12,12, 9) / 0.5203608945849587d-01 /
   data green(12,12,10) / 0.5074617404203012d-01 /
   data green(12,12,11) / 0.4942680345415775d-01 /
   data green(12,12,12) / 0.4809393985782028d-01 /

   data green(13, 0, 0) / 0.7703889286051638d-01 /
   data green(13, 1, 0) / 0.7680780784877862d-01 /
   data green(13, 1, 1) / 0.7657885648161755d-01 /
   data green(13, 2, 0) / 0.7612735480939857d-01 /
   data green(13, 2, 1) / 0.7590459662952359d-01 /
   data green(13, 2, 2) / 0.7524833172745232d-01 /
   data green(13, 3, 0) / 0.7503388018421212d-01 /
   data green(13, 3, 1) / 0.7482080992672581d-01 /
   data green(13, 3, 2) / 0.7419272062131753d-01 /
   data green(13, 3, 3) / 0.7318130773651561d-01 /
   data green(13, 4, 0) / 0.7358180949861125d-01 /
   data green(13, 4, 1) / 0.7338111821262081d-01 /
   data green(13, 4, 2) / 0.7278907741513559d-01 /
   data green(13, 4, 3) / 0.7183435369978027d-01 /
   data green(13, 4, 4) / 0.7056038791500062d-01 /
   data green(13, 5, 0) / 0.7183639760675302d-01 /
   data green(13, 5, 1) / 0.7164988357711846d-01 /
   data green(13, 5, 2) / 0.7109919274902818d-01 /
   data green(13, 5, 3) / 0.7020968495076238d-01 /
   data green(13, 5, 4) / 0.6901997231711041d-01 /
   data green(13, 5, 5) / 0.6757725655646277d-01 /
   data green(13, 6, 0) / 0.6986648648824506d-01 /
   data green(13, 6, 1) / 0.6969509509805502d-01 /
   data green(13, 6, 2) / 0.6918858725960286d-01 /
   data green(13, 6, 3) / 0.6836899083728129d-01 /
   data green(13, 6, 4) / 0.6727002313769816d-01 /
   data green(13, 6, 5) / 0.6593336931168965d-01 /
   data green(13, 6, 6) / 0.6440466342840299d-01 /
   data green(13, 7, 0) / 0.6773853969851734d-01 /
   data green(13, 7, 1) / 0.6758248943681053d-01 /
   data green(13, 7, 2) / 0.6712087922523301d-01 /
   data green(13, 7, 3) / 0.6637255979700815d-01 /
   data green(13, 7, 4) / 0.6536654354185854d-01 /
   data green(13, 7, 5) / 0.6413898353206067d-01 /
   data green(13, 7, 6) / 0.6272985285955927d-01 /
   data green(13, 7, 7) / 0.6117985897116904d-01 /
   data green(13, 8, 0) / 0.6551253322149389d-01 /
   data green(13, 8, 1) / 0.6537147941471408d-01 /
   data green(13, 8, 2) / 0.6495383450432213d-01 /
   data green(13, 8, 3) / 0.6427553935673967d-01 /
   data green(13, 8, 4) / 0.6336126542297392d-01 /
   data green(13, 8, 5) / 0.6224200292294987d-01 /
   data green(13, 8, 6) / 0.6095235287013211d-01 /
   data green(13, 8, 7) / 0.5952795421599506d-01 /
   data green(13, 8, 8) / 0.5800335834597368d-01 /
   data green(13, 9, 0) / 0.6323969853930737d-01 /
   data green(13, 9, 1) / 0.6311289983475619d-01 /
   data green(13, 9, 2) / 0.6273711613693157d-01 /
   data green(13, 9, 3) / 0.6212570788110104d-01 /
   data green(13, 9, 4) / 0.6129947758851307d-01 /
   data green(13, 9, 5) / 0.6028472092863377d-01 /
   data green(13, 9, 6) / 0.5911112149888466d-01 /
   data green(13, 9, 7) / 0.5780958664855207d-01 /
   data green(13, 9, 8) / 0.5641045008074607d-01 /
   data green(13, 9, 9) / 0.5494210921724269d-01 /
   data green(13,10, 0) / 0.6096177336742466d-01 /
   data green(13,10, 1) / 0.6084824131830983d-01 /
   data green(13,10, 2) / 0.6051147740455500d-01 /
   data green(13,10, 3) / 0.5996260987778013d-01 /
   data green(13,10, 4) / 0.5921904200534028d-01 /
   data green(13,10, 5) / 0.5830297184260553d-01 /
   data green(13,10, 6) / 0.5723964663155539d-01 /
   data green(13,10, 7) / 0.5605566269769010d-01 /
   data green(13,10, 8) / 0.5477745135858935d-01 /
   data green(13,10, 9) / 0.5343008136421920d-01 /
   data green(13,10,10) / 0.5203642093270375d-01 /
   data green(13,11, 0) / 0.5871131147958240d-01 /
   data green(13,11, 1) / 0.5860992763757320d-01 /
   data green(13,11, 2) / 0.5830894563701332d-01 /
   data green(13,11, 3) / 0.5781759540620831d-01 /
   data green(13,11, 4) / 0.5715038369964721d-01 /
   data green(13,11, 5) / 0.5632592909675871d-01 /
   data green(13,11, 6) / 0.5536560665383577d-01 /
   data green(13,11, 7) / 0.5429216180765198d-01 /
   data green(13,11, 8) / 0.5312847357179801d-01 /
   data green(13,11, 9) / 0.5189651806270230d-01 /
   data green(13,11,10) / 0.5061662988344395d-01 /
   data green(13,11,11) / 0.4930701504170588d-01 /
   data green(13,12, 0) / 0.5651261115736957d-01 /
   data green(13,12, 1) / 0.5642221695898542d-01 /
   data green(13,12, 2) / 0.5615365050733514d-01 /
   data green(13,12, 3) / 0.5571454670488019d-01 /
   data green(13,12, 4) / 0.5511696029537030d-01 /
   data green(13,12, 5) / 0.5437654736202797d-01 /
   data green(13,12, 6) / 0.5351107456354735d-01 /
   data green(13,12, 7) / 0.5254016839197281d-01 /
   data green(13,12, 8) / 0.5148344032567817d-01 /
   data green(13,12, 9) / 0.5036005053568500d-01 /
   data green(13,12,10) / 0.4918796525939784d-01 /
   data green(13,12,11) / 0.4798350354093404d-01 /
   data green(13,12,12) / 0.4676107633011698d-01 /
   data green(13,13, 0) / 0.5438293854794247d-01 /
   data green(13,13, 1) / 0.5430239580873501d-01 /
   data green(13,13, 2) / 0.5406292945926400d-01 /
   data green(13,13, 3) / 0.5367083569063271d-01 /
   data green(13,13, 4) / 0.5313612954572275d-01 /
   data green(13,13, 5) / 0.5247180487975907d-01 /
   data green(13,13, 6) / 0.5169301929977736d-01 /
   data green(13,13, 7) / 0.5081621825075060d-01 /
   data green(13,13, 8) / 0.4985829813304905d-01 /
   data green(13,13, 9) / 0.4883588613494076d-01 /
   data green(13,13,10) / 0.4776477385374071d-01 /
   data green(13,13,11) / 0.4665949931494588d-01 /
   data green(13,13,12) / 0.4553309191119101d-01 /
   data green(13,13,13) / 0.4439694584548196d-01 /

   data green(14, 0, 0) / 0.7152106737902052d-01 /
   data green(14, 1, 0) / 0.7133644846943507d-01 /
   data green(14, 1, 1) / 0.7115329116241083d-01 /
   data green(14, 2, 0) / 0.7079135807103833d-01 /
   data green(14, 2, 1) / 0.7061246305834240d-01 /
   data green(14, 2, 2) / 0.7008408174835440d-01 /
   data green(14, 3, 0) / 0.6991090029050440d-01 /
   data green(14, 3, 1) / 0.6973873586981673d-01 /
   data green(14, 3, 2) / 0.6923001537559118d-01 /
   data green(14, 3, 3) / 0.6840705362340681d-01 /
   data green(14, 4, 0) / 0.6873326165206863d-01 /
   data green(14, 4, 1) / 0.6856980798794006d-01 /
   data green(14, 4, 2) / 0.6808655775175572d-01 /
   data green(14, 4, 3) / 0.6730396207928187d-01 /
   data green(14, 4, 4) / 0.6625339181590958d-01 /
   data green(14, 5, 0) / 0.6730524087683480d-01 /
   data green(14, 5, 1) / 0.6715191811703045d-01 /
   data green(14, 5, 2) / 0.6669832368069810d-01 /
   data green(14, 5, 3) / 0.6596283314900651d-01 /
   data green(14, 5, 4) / 0.6497373989289887d-01 /
   data green(14, 5, 5) / 0.6376632297558286d-01 /
   data green(14, 6, 0) / 0.6567754974272230d-01 /
   data green(14, 6, 1) / 0.6553521859292928d-01 /
   data green(14, 6, 2) / 0.6511384427578483d-01 /
   data green(14, 6, 3) / 0.6442965547815595d-01 /
   data green(14, 6, 4) / 0.6350774816081910d-01 /
   data green(14, 6, 5) / 0.6237959694822497d-01 /
   data green(14, 6, 6) / 0.6108029045912914d-01 /
   data green(14, 7, 0) / 0.6390066305816181d-01 /
   data green(14, 7, 1) / 0.6376968542677226d-01 /
   data green(14, 7, 2) / 0.6338163439737193d-01 /
   data green(14, 7, 3) / 0.6275064014342555d-01 /
   data green(14, 7, 4) / 0.6189865252359852d-01 /
   data green(14, 7, 5) / 0.6085336848297798d-01 /
   data green(14, 7, 6) / 0.5964591141622549d-01 /
   data green(14, 7, 7) / 0.5830857838473431d-01 /
   data green(14, 8, 0) / 0.6202169259292453d-01 /
   data green(14, 8, 1) / 0.6190201825790708d-01 /
   data green(14, 8, 2) / 0.6154719053702444d-01 /
   data green(14, 8, 3) / 0.6096936559276550d-01 /
   data green(14, 8, 4) / 0.6018752220423795d-01 /
   data green(14, 8, 5) / 0.5922576056140686d-01 /
   data green(14, 8, 6) / 0.5811136947640624d-01 /
   data green(14, 8, 7) / 0.5687293808133353d-01 /
   data green(14, 8, 8) / 0.5553870562483351d-01 /
   data green(14, 9, 0) / 0.6008238144574942d-01 /
   data green(14, 9, 1) / 0.5997364878817465d-01 /
   data green(14, 9, 2) / 0.5965100821437892d-01 /
   data green(14, 9, 3) / 0.5912486557138306d-01 /
   data green(14, 9, 4) / 0.5841144320285861d-01 /
   data green(14, 9, 5) / 0.5753152253822041d-01 /
   data green(14, 9, 6) / 0.5650881589803534d-01 /
   data green(14, 9, 7) / 0.5536839257610165d-01 /
   data green(14, 9, 8) / 0.5413527324432756d-01 /
   data green(14, 9, 9) / 0.5283330745149469d-01 /
   data green(14,10, 0) / 0.5811812243149051d-01 /
   data green(14,10, 1) / 0.5801975368794781d-01 /
   data green(14,10, 2) / 0.5772766207525820d-01 /
   data green(14,10, 3) / 0.5725063650001787d-01 /
   data green(14,10, 4) / 0.5660250107278653d-01 /
   data green(14,10, 5) / 0.5580103112068265d-01 /
   data green(14,10, 6) / 0.5486667282456567d-01 /
   data green(14,10, 7) / 0.5382124152985159d-01 /
   data green(14,10, 8) / 0.5268673548953817d-01 /
   data green(14,10, 9) / 0.5148435924257055d-01 /
   data green(14,10,10) / 0.5023380395107076d-01 /
   data green(14,11, 0) / 0.5615776412464538d-01 /
   data green(14,11, 1) / 0.5606904720271656d-01 /
   data green(14,11, 2) / 0.5580543338970727d-01 /
   data green(14,11, 3) / 0.5537432781376447d-01 /
   data green(14,11, 4) / 0.5478742910652269d-01 /
   data green(14,11, 5) / 0.5405986640236535d-01 /
   data green(14,11, 6) / 0.5320916584769497d-01 /
   data green(14,11, 7) / 0.5225419710896263d-01 /
   data green(14,11, 8) / 0.5121417395083095d-01 /
   data green(14,11, 9) / 0.5010782843387915d-01 /
   data green(14,11,10) / 0.4895275493344969d-01 /
   data green(14,11,11) / 0.4776497198258449d-01 /
   data green(14,12, 0) / 0.5422395630267672d-01 /
   data green(14,12, 1) / 0.5414411200811974d-01 /
   data green(14,12, 2) / 0.5390670532987885d-01 /
   data green(14,12, 3) / 0.5351795786443523d-01 /
   data green(14,12, 4) / 0.5298773145086994d-01 /
   data green(14,12, 5) / 0.5232885052202225d-01 /
   data green(14,12, 6) / 0.5155627763397595d-01 /
   data green(14,12, 7) / 0.5068625338986604d-01 /
   data green(14,12, 8) / 0.4973547852248339d-01 /
   data green(14,12, 9) / 0.4872040509921060d-01 /
   data green(14,12,10) / 0.4765666252074589d-01 /
   data green(14,12,11) / 0.4655866142412146d-01 /
   data green(14,12,12) / 0.4543932724603427d-01 /
   data green(14,13, 0) / 0.5233380868254357d-01 /
   data green(14,13, 1) / 0.5226203885963645d-01 /
   data green(14,13, 2) / 0.5204850859807088d-01 /
   data green(14,13, 3) / 0.5169843299377453d-01 /
   data green(14,13, 4) / 0.5122010797380921d-01 /
   data green(14,13, 5) / 0.5062437620559734d-01 /
   data green(14,13, 6) / 0.4992397215999869d-01 /
   data green(14,13, 7) / 0.4913282817274011d-01 /
   data green(14,13, 8) / 0.4826540039071987d-01 /
   data green(14,13, 9) / 0.4733607223940629d-01 /
   data green(14,13,10) / 0.4635866283639944d-01 /
   data green(14,13,11) / 0.4534606137208064d-01 /
   data green(14,13,12) / 0.4430997696147362d-01 /
   data green(14,13,13) / 0.4326079968904013d-01 /
   data green(14,14, 0) / 0.5049968959520108d-01 /
   data green(14,14, 1) / 0.5043521193844083d-01 /
   data green(14,14, 2) / 0.5024326658770582d-01 /
   data green(14,14, 3) / 0.4992822183372374d-01 /
   data green(14,14, 4) / 0.4949704607535031d-01 /
   data green(14,14, 5) / 0.4895889776851988d-01 /
   data green(14,14, 6) / 0.4832459056585090d-01 /
   data green(14,14, 7) / 0.4760604922941294d-01 /
   data green(14,14, 8) / 0.4681575164406249d-01 /
   data green(14,14, 9) / 0.4596622851054324d-01 /
   data green(14,14,10) / 0.4506964905799747d-01 /
   data green(14,14,11) / 0.4413748707906946d-01 /
   data green(14,14,12) / 0.4318029388013415d-01 /
   data green(14,14,13) / 0.4220755165047038d-01 /
   data green(14,14,14) / 0.4122760615822885d-01 /

   data green(15, 0, 0) / 0.6674171875620340d-01 /
   data green(15, 1, 0) / 0.6659187626741821d-01 /
   data green(15, 1, 1) / 0.6644306187808230d-01 /
   data green(15, 2, 0) / 0.6614851872448553d-01 /
   data green(15, 2, 1) / 0.6600271652459495d-01 /
   data green(15, 2, 2) / 0.6557119780529386d-01 /
   data green(15, 3, 0) / 0.6542943043154954d-01 /
   data green(15, 3, 1) / 0.6528841873268112d-01 /
   data green(15, 3, 2) / 0.6487094337922872d-01 /
   data green(15, 3, 3) / 0.6419305893649080d-01 /
   data green(15, 4, 0) / 0.6446200317440251d-01 /
   data green(15, 4, 1) / 0.6432725556938332d-01 /
   data green(15, 4, 2) / 0.6392815733194791d-01 /
   data green(15, 4, 3) / 0.6327958205488594d-01 /
   data green(15, 4, 4) / 0.6240458491748328d-01 /
   data green(15, 5, 0) / 0.6328040403480169d-01 /
   data green(15, 5, 1) / 0.6315303360801033d-01 /
   data green(15, 5, 2) / 0.6277559682964562d-01 /
   data green(15, 5, 3) / 0.6216162923380718d-01 /
   data green(15, 5, 4) / 0.6133217276781156d-01 /
   data green(15, 5, 5) / 0.6031381745761460d-01 /
   data green(15, 6, 0) / 0.6192247000281843d-01 /
   data green(15, 6, 1) / 0.6180321891405282d-01 /
   data green(15, 6, 2) / 0.6144964590958335d-01 /
   data green(15, 6, 3) / 0.6087387607540391d-01 /
   data green(15, 6, 4) / 0.6009482044067303d-01 /
   data green(15, 6, 5) / 0.5913648671195317d-01 /
   data green(15, 6, 6) / 0.5802606125039635d-01 /
   data green(15, 7, 0) / 0.6042682501769198d-01 /
   data green(15, 7, 1) / 0.6031608892723563d-01 /
   data green(15, 7, 2) / 0.5998756755294807d-01 /
   data green(15, 7, 3) / 0.5945197784349463d-01 /
   data green(15, 7, 4) / 0.5872609139627624d-01 /
   data green(15, 7, 5) / 0.5783130821831069d-01 /
   data green(15, 7, 6) / 0.5679200983354748d-01 /
   data green(15, 7, 7) / 0.5563391259654712d-01 /
   data green(15, 8, 0) / 0.5883054056663740d-01 /
   data green(15, 8, 1) / 0.5872841553574899d-01 /
   data green(15, 8, 2) / 0.5842525730505976d-01 /
   data green(15, 8, 3) / 0.5793042949678512d-01 /
   data green(15, 8, 4) / 0.5725864842363673d-01 /
   data green(15, 8, 5) / 0.5642877666766111d-01 /
   data green(15, 8, 6) / 0.5546243506118998d-01 /
   data green(15, 8, 7) / 0.5438262937840914d-01 /
   data green(15, 8, 8) / 0.5321244867937117d-01 /
   data green(15, 9, 0) / 0.5716747818605419d-01 /
   data green(15, 9, 1) / 0.5707382135237229d-01 /
   data green(15, 9, 2) / 0.5679563190406549d-01 /
   data green(15, 9, 3) / 0.5634101957273299d-01 /
   data green(15, 9, 4) / 0.5572276898064677d-01 /
   data green(15, 9, 5) / 0.5495736571903830d-01 /
   data green(15, 9, 6) / 0.5406381760529760d-01 /
   data green(15, 9, 7) / 0.5306250075899824d-01 /
   data green(15, 9, 8) / 0.5197405216562960d-01 /
   data green(15, 9, 9) / 0.5081846716524512d-01 /
   data green(15,10, 0) / 0.5546730897034690d-01 /
   data green(15,10, 1) / 0.5538179897216711d-01 /
   data green(15,10, 2) / 0.5512765596265851d-01 /
   data green(15,10, 3) / 0.5471185459646648d-01 /
   data green(15,10, 4) / 0.5414542573182933d-01 /
   data green(15,10, 5) / 0.5344265264292526d-01 /
   data green(15,10, 6) / 0.5262013628417619d-01 /
   data green(15,10, 7) / 0.5169578736931021d-01 /
   data green(15,10, 8) / 0.5068791149732016d-01 /
   data green(15,10, 9) / 0.4961440920393869d-01 /
   data green(15,10,10) / 0.4849217111094843d-01 /
   data green(15,11, 0) / 0.5375509994314404d-01 /
   data green(15,11, 1) / 0.5367729297734811d-01 /
   data green(15,11, 2) / 0.5344591211035581d-01 /
   data green(15,11, 3) / 0.5306691385409971d-01 /
   data green(15,11, 4) / 0.5254976929865351d-01 /
   data green(15,11, 5) / 0.5190679051068153d-01 /
   data green(15,11, 6) / 0.5115237592932710d-01 /
   data green(15,11, 7) / 0.5030217909766109d-01 /
   data green(15,11, 8) / 0.4937232974689406d-01 /
   data green(15,11, 9) / 0.4837875535715296d-01 /
   data green(15,11,10) / 0.4733663250489944d-01 /
   data green(15,11,11) / 0.4625998666427240d-01 /
   data green(15,12, 0) / 0.5205132976719224d-01 /
   data green(15,12, 1) / 0.5198070713377539d-01 /
   data green(15,12, 2) / 0.5177057177370504d-01 /
   data green(15,12, 3) / 0.5142600326328514d-01 /
   data green(15,12, 4) / 0.5095508469594236d-01 /
   data green(15,12, 5) / 0.5036839000487217d-01 /
   data green(15,12, 6) / 0.4967834876824720d-01 /
   data green(15,12, 7) / 0.4889857276229880d-01 /
   data green(15,12, 8) / 0.4804320466212800d-01 /
   data green(15,12, 9) / 0.4712633518823360d-01 /
   data green(15,12,10) / 0.4616152497708758d-01 /
   data green(15,12,11) / 0.4516144063839250d-01 /
   data green(15,12,12) / 0.4413761040741662d-01 /
   data green(15,13, 0) / 0.5037218983506634d-01 /
   data green(15,13, 1) / 0.5030819577830918d-01 /
   data green(15,13, 2) / 0.5011768321366715d-01 /
   data green(15,13, 3) / 0.4980496727161907d-01 /
   data green(15,13, 4) / 0.4937693044528065d-01 /
   data green(15,13, 5) / 0.4884262412773155d-01 /
   data green(15,13, 6) / 0.4821274121040875d-01 /
   data green(15,13, 7) / 0.4749907501888514d-01 /
   data green(15,13, 8) / 0.4671397395469189d-01 /
   data green(15,13, 9) / 0.4586985003808540d-01 /
   data green(15,13,10) / 0.4497875876463045d-01 /
   data green(15,13,11) / 0.4405208065595908d-01 /
   data green(15,13,12) / 0.4310028697684704d-01 /
   data green(15,13,13) / 0.4213279564511804d-01 /
   data green(15,14, 0) / 0.4873004822033072d-01 /
   data green(15,14, 1) / 0.4867211887171267d-01 /
   data green(15,14, 2) / 0.4849957494990365d-01 /
   data green(15,14, 3) / 0.4821607432849154d-01 /
   data green(15,14, 4) / 0.4782747401865364d-01 /
   data green(15,14, 5) / 0.4734149736357721d-01 /
   data green(15,14, 6) / 0.4676732632738748d-01 /
   data green(15,14, 7) / 0.4611514877118640d-01 /
   data green(15,14, 8) / 0.4539571341238746d-01 /
   data green(15,14, 9) / 0.4461991250424643d-01 /
   data green(15,14,10) / 0.4379842510874424d-01 /
   data green(15,14,11) / 0.4294142817619453d-01 /
   data green(15,14,12) / 0.4205838565540346d-01 /
   data green(15,14,13) / 0.4115790515090587d-01 /
   data green(15,14,14) / 0.4024766066184958d-01 /
   data green(15,15, 0) / 0.4713398637027776d-01 /
   data green(15,15, 1) / 0.4708156970795920d-01 /
   data green(15,15, 2) / 0.4692537241822766d-01 /
   data green(15,15, 3) / 0.4666849437149211d-01 /
   data green(15,15, 4) / 0.4631591110315855d-01 /
   data green(15,15, 5) / 0.4587421042808160d-01 /
   data green(15,15, 6) / 0.4535126072662866d-01 /
   data green(15,15, 7) / 0.4475584974911247d-01 /
   data green(15,15, 8) / 0.4409731117925764d-01 /
   data green(15,15, 9) / 0.4338517668094200d-01 /
   data green(15,15,10) / 0.4262887742715928d-01 /
   data green(15,15,11) / 0.4183748511320042d-01 /
   data green(15,15,12) / 0.4101951111721654d-01 /
   data green(15,15,13) / 0.4018278410841121d-01 /
   data green(15,15,14) / 0.3933436638804425d-01 /
   data green(15,15,15) / 0.3848050739283826d-01 /

   data green(16, 0, 0) / 0.6256173785534200d-01 /
   data green(16, 1, 0) / 0.6243844639174706d-01 /
   data green(16, 1, 1) / 0.6231589398967027d-01 /
   data green(16, 2, 0) / 0.6207301304091344d-01 /
   data green(16, 2, 1) / 0.6195263888283642d-01 /
   data green(16, 2, 2) / 0.6159578784374864d-01 /
   data green(16, 3, 0) / 0.6147831824980046d-01 /
   data green(16, 3, 1) / 0.6136147417263093d-01 /
   data green(16, 3, 2) / 0.6101481140651049d-01 /
   data green(16, 3, 3) / 0.6045025840164824d-01 /
   data green(16, 4, 0) / 0.6067441210999405d-01 /
   data green(16, 4, 1) / 0.6056211184627168d-01 /
   data green(16, 4, 2) / 0.6022900533937064d-01 /
   data green(16, 4, 3) / 0.5968611484000638d-01 /
   data green(16, 4, 4) / 0.5895060206860541d-01 /
   data green(16, 5, 0) / 0.5968664048187906d-01 /
   data green(16, 5, 1) / 0.5957980867166215d-01 /
   data green(16, 5, 2) / 0.5926279033867567d-01 /
   data green(16, 5, 3) / 0.5874571452491548d-01 /
   data green(16, 5, 4) / 0.5804444603921494d-01 /
   data green(16, 5, 5) / 0.5717926111522229d-01 /
   data green(16, 6, 0) / 0.5854360263478289d-01 /
   data green(16, 6, 1) / 0.5844285666000566d-01 /
   data green(16, 6, 2) / 0.5814377143531852d-01 /
   data green(16, 6, 3) / 0.5765554620191340d-01 /
   data green(16, 6, 4) / 0.5699253953850877d-01 /
   data green(16, 6, 5) / 0.5617332413449223d-01 /
   data green(16, 6, 6) / 0.5521909022886302d-01 /
   data green(16, 7, 0) / 0.5727509978117892d-01 /
   data green(16, 7, 1) / 0.5718082130192768d-01 /
   data green(16, 7, 2) / 0.5690079674166779d-01 /
   data green(16, 7, 3) / 0.5644324417841674d-01 /
   data green(16, 7, 4) / 0.5582113875623978d-01 /
   data green(16, 7, 5) / 0.5505113784397305d-01 /
   data green(16, 7, 6) / 0.5415242088301828d-01 /
   data green(16, 7, 7) / 0.5314562222955732d-01 /
   data green(16, 8, 0) / 0.5591041657697739d-01 /
   data green(16, 8, 1) / 0.5582276675261970d-01 /
   data green(16, 8, 2) / 0.5556230863531540d-01 /
   data green(16, 8, 3) / 0.5513629980544873d-01 /
   data green(16, 8, 4) / 0.5455629991021001d-01 /
   data green(16, 8, 5) / 0.5383710873660453d-01 /
   data green(16, 8, 6) / 0.5299600423868570d-01 /
   data green(16, 8, 7) / 0.5205154089892618d-01 /
   data green(16, 8, 8) / 0.5102262395054150d-01 /
   data green(16, 9, 0) / 0.5447699555932852d-01 /
   data green(16, 9, 1) / 0.5439595611034077d-01 /
   data green(16, 9, 2) / 0.5415500622113444d-01 /
   data green(16, 9, 3) / 0.5376057128533996d-01 /
   data green(16, 9, 4) / 0.5322273383505084d-01 /
   data green(16, 9, 5) / 0.5255467073550780d-01 /
   data green(16, 9, 6) / 0.5177165040662200d-01 /
   data green(16, 9, 7) / 0.5089033163751801d-01 /
   data green(16, 9, 8) / 0.4992770978282948d-01 /
   data green(16, 9, 9) / 0.4890056655482078d-01 /
   data green(16,10, 0) / 0.5299955439372710d-01 /
   data green(16,10, 1) / 0.5292496046331108d-01 /
   data green(16,10, 2) / 0.5270306376166042d-01 /
   data green(16,10, 3) / 0.5233948868884210d-01 /
   data green(16,10, 4) / 0.5184298191078303d-01 /
   data green(16,10, 5) / 0.5122516421019803d-01 /
   data green(16,10, 6) / 0.5049951630795567d-01 /
   data green(16,10, 7) / 0.4968076685985674d-01 /
   data green(16,10, 8) / 0.4878414966205977d-01 /
   data green(16,10, 9) / 0.4782479327295563d-01 /
   data green(16,10,10) / 0.4681709294461325d-01 /
   data green(16,11, 0) / 0.5149961051240483d-01 /
   data green(16,11, 1) / 0.5143119600163065d-01 /
   data green(16,11, 2) / 0.5122762579619197d-01 /
   data green(16,11, 3) / 0.5089361960634859d-01 /
   data green(16,11, 4) / 0.5043696592747502d-01 /
   data green(16,11, 5) / 0.4986771945218328d-01 /
   data green(16,11, 6) / 0.4919764038077680d-01 /
   data green(16,11, 7) / 0.4843981579138766d-01 /
   data green(16,11, 8) / 0.4760775523686294d-01 /
   data green(16,11, 9) / 0.4671499362988930d-01 /
   data green(16,11,10) / 0.4577460017821876d-01 /
   data green(16,11,11) / 0.4479880759614868d-01 /
   data green(16,12, 0) / 0.4999531828354619d-01 /
   data green(16,12, 1) / 0.4993274040589635d-01 /
   data green(16,12, 2) / 0.4974642365558785d-01 /
   data green(16,12, 3) / 0.4944051496950049d-01 /
   data green(16,12, 4) / 0.4902167971325357d-01 /
   data green(16,12, 5) / 0.4849864188848311d-01 /
   data green(16,12, 6) / 0.4788174118665058d-01 /
   data green(16,12, 7) / 0.4718238159583703d-01 /
   data green(16,12, 8) / 0.4641253189972664d-01 /
   data green(16,12, 9) / 0.4558427478044640d-01 /
   data green(16,12,10) / 0.4470932087589460d-01 /
   data green(16,12,11) / 0.4379877774565288d-01 /
   data green(16,12,12) / 0.4286288663955817d-01 /
   data green(16,13, 0) / 0.4850156047567791d-01 /
   data green(16,13, 1) / 0.4844443684197702d-01 /
   data green(16,13, 2) / 0.4827427728529055d-01 /
   data green(16,13, 3) / 0.4799466606651963d-01 /
   data green(16,13, 4) / 0.4761132876340764d-01 /
   data green(16,13, 5) / 0.4713180490596692d-01 /
   data green(16,13, 6) / 0.4656510534579313d-01 /
   data green(16,13, 7) / 0.4592118865082891d-01 /
   data green(16,13, 8) / 0.4521060783404478d-01 /
   data green(16,13, 9) / 0.4444404600888210d-01 /
   data green(16,13,10) / 0.4363200909176084d-01 /
   data green(16,13,11) / 0.4278451130381724d-01 /
   data green(16,13,12) / 0.4191086791816503d-01 /
   data green(16,13,13) / 0.4101958930345760d-01 /
   data green(16,14, 0) / 0.4703018608507870d-01 /
   data green(16,14, 1) / 0.4697811254212711d-01 /
   data green(16,14, 2) / 0.4682294346210170d-01 /
   data green(16,14, 3) / 0.4656770855738484d-01 /
   data green(16,14, 4) / 0.4621738290485727d-01 /
   data green(16,14, 5) / 0.4577846677038516d-01 /
   data green(16,14, 6) / 0.4525871988411324d-01 /
   data green(16,14, 7) / 0.4466688601074147d-01 /
   data green(16,14, 8) / 0.4401217707200717d-01 /
   data green(16,14, 9) / 0.4330407962672791d-01 /
   data green(16,14,10) / 0.4255192171875387d-01 /
   data green(16,14,11) / 0.4176470641331672d-01 /
   data green(16,14,12) / 0.4095089051355488d-01 /
   data green(16,14,13) / 0.4011825293460711d-01 /
   data green(16,14,14) / 0.3927381459573205d-01 /
   data green(16,15, 0) / 0.4559034080406678d-01 /
   data green(16,15, 1) / 0.4554291046935628d-01 /
   data green(16,15, 2) / 0.4540150900068592d-01 /
   data green(16,15, 3) / 0.4516872817858282d-01 /
   data green(16,15, 4) / 0.4484890654432647d-01 /
   data green(16,15, 5) / 0.4444755660324370d-01 /
   data green(16,15, 6) / 0.4397142013050519d-01 /
   data green(16,15, 7) / 0.4342811437260578d-01 /
   data green(16,15, 8) / 0.4282570604354332d-01 /
   data green(16,15, 9) / 0.4217254014370911d-01 /
   data green(16,15,10) / 0.4147687226494933d-01 /
   data green(16,15,11) / 0.4074682944127523d-01 /
   data green(16,15,12) / 0.3999001530951307d-01 /
   data green(16,15,13) / 0.3921351027966460d-01 /
   data green(16,15,14) / 0.3842379284910723d-01 /
   data green(16,15,15) / 0.3762664807465570d-01 /
   data green(16,16, 0) / 0.4418883756063733d-01 /
   data green(16,16, 1) / 0.4414565200885717d-01 /
   data green(16,16, 2) / 0.4401685589181507d-01 /
   data green(16,16, 3) / 0.4380470141208303d-01 /
   data green(16,16, 4) / 0.4351280552344674d-01 /
   data green(16,16, 5) / 0.4314600966900181d-01 /
   data green(16,16, 6) / 0.4271013877248645d-01 /
   data green(16,16, 7) / 0.4221176390703201d-01 /
   data green(16,16, 8) / 0.4165795894672368d-01 /
   data green(16,16, 9) / 0.4105604547294645d-01 /
   data green(16,16,10) / 0.4041336165945691d-01 /
   data green(16,16,11) / 0.3973714818851525d-01 /
   data green(16,16,12) / 0.3903424476096395d-01 /
   data green(16,16,13) / 0.3831107516538472d-01 /
   data green(16,16,14) / 0.3757360649754232d-01 /
   data green(16,16,15) / 0.3682717433535779d-01 /
   data green(16,16,16) / 0.3607655607088319d-01 /
   data green(17, 0, 0) / 0.5887493353388418d-01 /
   data green(17, 1, 0) / 0.5877226104072109d-01 /
   data green(17, 1, 1) / 0.5867013568763937d-01 /
   data green(17, 2, 0) / 0.5846752810364193d-01 /
   data green(17, 2, 1) / 0.5836699854695342d-01 /
   data green(17, 2, 2) / 0.5806859127542865d-01 /
   data green(17, 3, 0) / 0.5797021626458378d-01 /
   data green(17, 3, 1) / 0.5787227807355474d-01 /
   data green(17, 3, 2) / 0.5758147047665457d-01 /
   data green(17, 3, 3) / 0.5710660228451652d-01 /
   data green(17, 4, 0) / 0.5729527820732241d-01 /
   data green(17, 4, 1) / 0.5720076444608333d-01 /
   data green(17, 4, 2) / 0.5692006421652476d-01 /
   data green(17, 4, 3) / 0.5646145230555712d-01 /
   data green(17, 4, 4) / 0.5583797030553406d-01 /
   data green(17, 5, 0) / 0.5646181387001448d-01 /
   data green(17, 5, 1) / 0.5637141448056780d-01 /
   data green(17, 5, 2) / 0.5610287530693918d-01 /
   data green(17, 5, 3) / 0.5566379001382651d-01 /
   data green(17, 5, 4) / 0.5506636902089202d-01 /
   data green(17, 5, 5) / 0.5432622827887320d-01 /
   data green(17, 6, 0) / 0.5549169293687296d-01 /
   data green(17, 6, 1) / 0.5540591908181451d-01 /
   data green(17, 6, 2) / 0.5515100409178197d-01 /
   data green(17, 6, 3) / 0.5473399283452457d-01 /
   data green(17, 6, 4) / 0.5416601154013759d-01 /
   data green(17, 6, 5) / 0.5346143868309417d-01 /
   data green(17, 6, 6) / 0.5263699003776753d-01 /
   data green(17, 7, 0) / 0.5440811393786864d-01 /
   data green(17, 7, 1) / 0.5432731011183933d-01 /
   data green(17, 7, 2) / 0.5408707803184717d-01 /
   data green(17, 7, 3) / 0.5369379699307785d-01 /
   data green(17, 7, 4) / 0.5315751440066769d-01 /
   data green(17, 7, 5) / 0.5249140399071511d-01 /
   data green(17, 7, 6) / 0.5171068048732165d-01 /
   data green(17, 7, 7) / 0.5083191182490546d-01 /
   data green(17, 8, 0) / 0.5323433132226412d-01 /
   data green(17, 8, 1) / 0.5315868319020285d-01 /
   data green(17, 8, 2) / 0.5293368909872921d-01 /
   data green(17, 8, 3) / 0.5256506601450348d-01 /
   data green(17, 8, 4) / 0.5206183659795096d-01 /
   data green(17, 8, 5) / 0.5143582530128132d-01 /
   data green(17, 8, 6) / 0.5070089593107328d-01 /
   data green(17, 8, 7) / 0.4987202653867786d-01 /
   data green(17, 8, 8) / 0.4896479191389495d-01 /
   data green(17, 9, 0) / 0.5199261309253627d-01 /
   data green(17, 9, 1) / 0.5192216711364216d-01 /
   data green(17, 9, 2) / 0.5171255750827315d-01 /
   data green(17, 9, 3) / 0.5136884531018861d-01 /
   data green(17, 9, 4) / 0.5089911606230395d-01 /
   data green(17, 9, 5) / 0.5031390073357141d-01 /
   data green(17, 9, 6) / 0.4962559731470560d-01 /
   data green(17, 9, 7) / 0.4884763030747716d-01 /
   data green(17, 9, 8) / 0.4799454894596799d-01 /
   data green(17, 9, 9) / 0.4707995355480811d-01 /
   data green(17,10, 0) / 0.5070349335486637d-01 /
   data green(17,10, 1) / 0.5063818302354651d-01 /
   data green(17,10, 2) / 0.5044377037107320d-01 /
   data green(17,10, 3) / 0.5012473487807578d-01 /
   data green(17,10, 4) / 0.4968823496561772d-01 /
   data green(17,10, 5) / 0.4914349839237318d-01 /
   data green(17,10, 6) / 0.4850172140260486d-01 /
   data green(17,10, 7) / 0.4777497376996495d-01 /
   data green(17,10, 8) / 0.4697597082927561d-01 /
   data green(17,10, 9) / 0.4611749201508532d-01 /
   data green(17,10,10) / 0.4521186575827368d-01 /
   data green(17,11, 0) / 0.4938530141002867d-01 /
   data green(17,11, 1) / 0.4932497198848633d-01 /
   data green(17,11, 2) / 0.4914532026757512d-01 /
   data green(17,11, 3) / 0.4885027916224602d-01 /
   data green(17,11, 4) / 0.4844601044378896d-01 /
   data green(17,11, 5) / 0.4794093907460038d-01 /
   data green(17,11, 6) / 0.4734469816219383d-01 /
   data green(17,11, 7) / 0.4666816498819426d-01 /
   data green(17,11, 8) / 0.4592271032413171d-01 /
   data green(17,11, 9) / 0.4511982713153132d-01 /
   data green(17,11,10) / 0.4427071950637440d-01 /
   data green(17,11,11) / 0.4338608119433198d-01 /
   data green(17,12, 0) / 0.4805392731762422d-01 /
   data green(17,12, 1) / 0.4799836308571572d-01 /
   data green(17,12, 2) / 0.4783282819587404d-01 /
   data green(17,12, 3) / 0.4756069866536667d-01 /
   data green(17,12, 4) / 0.4718756415730453d-01 /
   data green(17,12, 5) / 0.4672054402247768d-01 /
   data green(17,12, 6) / 0.4616829490843453d-01 /
   data green(17,12, 7) / 0.4554040841149645d-01 /
   data green(17,12, 8) / 0.4484700708067738d-01 /
   data green(17,12, 9) / 0.4409840205136780d-01 /
   data green(17,12,10) / 0.4330473071615870d-01 /
   data green(17,12,11) / 0.4247571819889255d-01 /
   data green(17,12,12) / 0.4162040424971019d-01 /
   data green(17,13, 0) / 0.4672279502570380d-01 /
   data green(17,13, 1) / 0.4667173044232976d-01 /
   data green(17,13, 2) / 0.4651954315419128d-01 /
   data green(17,13, 3) / 0.4626921409766704d-01 /
   data green(17,13, 4) / 0.4592547518394455d-01 /
   data green(17,13, 5) / 0.4549470203837398d-01 /
   data green(17,13, 6) / 0.4498441486762699d-01 /
   data green(17,13, 7) / 0.4440308010579855d-01 /
   data green(17,13, 8) / 0.4375969682307911d-01 /
   data green(17,13, 9) / 0.4306346361016566d-01 /
   data green(17,13,10) / 0.4232349839097180d-01 /
   data green(17,13,11) / 0.4154860757859177d-01 /
   data green(17,13,12) / 0.4074704535407191d-01 /
   data green(17,13,13) / 0.3992651571646749d-01 /
   data green(17,14, 0) / 0.4540294021332417d-01 /
   data green(17,14, 1) / 0.4535608857278361d-01 /
   data green(17,14, 2) / 0.4521640388167961d-01 /
   data green(17,14, 3) / 0.4498647041313176d-01 /
   data green(17,14, 4) / 0.4467041979012532d-01 /
   data green(17,14, 5) / 0.4427376587555599d-01 /
   data green(17,14, 6) / 0.4380312285573023d-01 /
   data green(17,14, 7) / 0.4326607725043573d-01 /
   data green(17,14, 8) / 0.4267011503076935d-01 /
   data green(17,14, 9) / 0.4202391159334237d-01 /
   data green(17,14,10) / 0.4133542222461956d-01 /
   data green(17,14,11) / 0.4061265905451390d-01 /
   data green(17,14,12) / 0.3986316800096901d-01 /
   data green(17,14,13) / 0.3909452798202293d-01 /
   data green(17,14,14) / 0.3831122569940225d-01 /
   data green(17,15, 0) / 0.4410321569033952d-01 /
   data green(17,15, 1) / 0.4406027811665782d-01 /
   data green(17,15, 2) / 0.4393222363733416d-01 /
   data green(17,15, 3) / 0.4372128475333481d-01 /
   data green(17,15, 4) / 0.4343100890052237d-01 /
   data green(17,15, 5) / 0.4306623667591733d-01 /
   data green(17,15, 6) / 0.4263276021834217d-01 /
   data green(17,15, 7) / 0.4213705761468117d-01 /
   data green(17,15, 8) / 0.4158614238615563d-01 /
   data green(17,15, 9) / 0.4098727524690763d-01 /
   data green(17,15,10) / 0.4034776767784514d-01 /
   data green(17,15,11) / 0.3967475976571049d-01 /
   data green(17,15,12) / 0.3897507206847028d-01 /
   data green(17,15,13) / 0.3825513494823103d-01 /
   data green(17,15,14) / 0.3752081328271912d-01 /
   data green(17,15,15) / 0.3677744566329784d-01 /
   data green(17,16, 0) / 0.4283052065609781d-01 /
   data green(17,16, 1) / 0.4279119777941231d-01 /
   data green(17,16, 2) / 0.4267388256869234d-01 /
   data green(17,16, 3) / 0.4248049329759331d-01 /
   data green(17,16, 4) / 0.4221415575166544d-01 /
   data green(17,16, 5) / 0.4187900910842762d-01 /
   data green(17,16, 6) / 0.4148010128387013d-01 /
   data green(17,16, 7) / 0.4102312937711716d-01 /
   data green(17,16, 8) / 0.4051426661919115d-01 /
   data green(17,16, 9) / 0.3995993823001194d-01 /
   data green(17,16,10) / 0.3936663924520244d-01 /
   data green(17,16,11) / 0.3874080590713080d-01 /
   data green(17,16,12) / 0.3808858643477965d-01 /
   data green(17,16,13) / 0.3741580027282519d-01 /
   data green(17,16,14) / 0.3672789943251679d-01 /
   data green(17,16,15) / 0.3602979433095416d-01 /
   data green(17,16,16) / 0.3532593790517709d-01 /
   data green(17,17, 0) / 0.4159006209256801d-01 /
   data green(17,17, 1) / 0.4155406018697424d-01 /
   data green(17,17, 2) / 0.4144661866078738d-01 /
   data green(17,17, 3) / 0.4126939534253214d-01 /
   data green(17,17, 4) / 0.4102508994270987d-01 /
   data green(17,17, 5) / 0.4071730343901121d-01 /
   data green(17,17, 6) / 0.4035043068719424d-01 /
   data green(17,17, 7) / 0.3992944339086500d-01 /
   data green(17,17, 8) / 0.3945975286993446d-01 /
   data green(17,17, 9) / 0.3894706971456777d-01 /
   data green(17,17,10) / 0.3839715881851206d-01 /
   data green(17,17,11) / 0.3781575950497194d-01 /
   data green(17,17,12) / 0.3720841689190447d-01 /
   data green(17,17,13) / 0.3658044007026569d-01 /
   data green(17,17,14) / 0.3593679533093446d-01 /
   data green(17,17,15) / 0.3528202046199543d-01 /
   data green(17,17,16) / 0.3462027113800623d-01 /
   data green(17,17,17) / 0.3395524878883717d-01 /

   data green(18, 0, 0) / 0.5559881044038090d-01 /
   data green(18, 1, 0) / 0.5551239783537756d-01 /
   data green(18, 1, 1) / 0.5542640139406329d-01 /
   data green(18, 2, 0) / 0.5525562275966436d-01 /
   data green(18, 2, 1) / 0.5517082896385630d-01 /
   data green(18, 2, 2) / 0.5491880411949875d-01 /
   data green(18, 3, 0) / 0.5483562665466853d-01 /
   data green(18, 3, 1) / 0.5475277441927072d-01 /
   data green(18, 3, 2) / 0.5450650013975773d-01 /
   data green(18, 3, 3) / 0.5410346084434892d-01 /
   data green(18, 4, 0) / 0.5426370122299168d-01 /
   data green(18, 4, 1) / 0.5418344445648763d-01 /
   data green(18, 4, 2) / 0.5394483535830598d-01 /
   data green(18, 4, 3) / 0.5355423293266306d-01 /
   data green(18, 4, 4) / 0.5302152029728531d-01 /
   data green(18, 5, 0) / 0.5355444172502194d-01 /
   data green(18, 5, 1) / 0.5347733347075403d-01 /
   data green(18, 5, 2) / 0.5324800435614933d-01 /
   data green(18, 5, 3) / 0.5287239552775495d-01 /
   data green(18, 5, 4) / 0.5235980974810366d-01 /
   data green(18, 5, 5) / 0.5172248312126407d-01 /
   data green(18, 6, 0) / 0.5272478421375259d-01 /
   data green(18, 6, 1) / 0.5265123073319172d-01 /
   data green(18, 6, 2) / 0.5243241937505021d-01 /
   data green(18, 6, 3) / 0.5207385762736174d-01 /
   data green(18, 6, 4) / 0.5158416219448582d-01 /
   data green(18, 6, 5) / 0.5097467028443662d-01 /
   data green(18, 6, 6) / 0.5025860755858098d-01 /
   data green(18, 7, 0) / 0.5179293694775466d-01 /
   data green(18, 7, 1) / 0.5172324790197263d-01 /
   data green(18, 7, 2) / 0.5151589208191134d-01 /
   data green(18, 7, 3) / 0.5117582062599070d-01 /
   data green(18, 7, 4) / 0.5071099857021139d-01 /
   data green(18, 7, 5) / 0.5013182652509036d-01 /
   data green(18, 7, 6) / 0.4945047507932946d-01 /
   data green(18, 7, 7) / 0.4868037470804959d-01 /
   data green(18, 8, 0) / 0.5077747718396287d-01 /
   data green(18, 8, 1) / 0.5071183771802008d-01 /
   data green(18, 8, 2) / 0.5051644688389025d-01 /
   data green(18, 8, 3) / 0.5019582330021627d-01 /
   data green(18, 8, 4) / 0.4975716951277947d-01 /
   data green(18, 8, 5) / 0.4920989705775394d-01 /
   data green(18, 8, 6) / 0.4856520345829524d-01 /
   data green(18, 8, 7) / 0.4783527942550182d-01 /
   data green(18, 8, 8) / 0.4703296794010305d-01 /
   data green(18, 9, 0) / 0.4969652385771705d-01 /
   data green(18, 9, 1) / 0.4963500968753524d-01 /
   data green(18, 9, 2) / 0.4945185687290201d-01 /
   data green(18, 9, 3) / 0.4915106634242440d-01 /
   data green(18, 9, 4) / 0.4873918952468697d-01 /
   data green(18, 9, 5) / 0.4822468103987832d-01 /
   data green(18, 9, 6) / 0.4761764194758322d-01 /
   data green(18, 9, 7) / 0.4692919210294014d-01 /
   data green(18, 9, 8) / 0.4617106840883924d-01 /
   data green(18, 9, 9) / 0.4535502848547590d-01 /
   data green(18,10, 0) / 0.4856712395319882d-01 /
   data green(18,10, 1) / 0.4850972997899163d-01 /
   data green(18,10, 2) / 0.4833875428645069d-01 /
   data green(18,10, 3) / 0.4805782647581010d-01 /
   data green(18,10, 4) / 0.4767275968949198d-01 /
   data green(18,10, 5) / 0.4719112636541165d-01 /
   data green(18,10, 6) / 0.4662196665594107d-01 /
   data green(18,10, 7) / 0.4597541541279396d-01 /
   data green(18,10, 8) / 0.4526201354865622d-01 /
   data green(18,10, 9) / 0.4449256965094449d-01 /
   data green(18,10,10) / 0.4367761465062792d-01 /
   data green(18,11, 0) / 0.4740483003571147d-01 /
   data green(18,11, 1) / 0.4735147337390667d-01 /
   data green(18,11, 2) / 0.4719248994793394d-01 /
   data green(18,11, 3) / 0.4693107422428211d-01 /
   data green(18,11, 4) / 0.4657239169403432d-01 /
   data green(18,11, 5) / 0.4612311203206099d-01 /
   data green(18,11, 6) / 0.4559146811412607d-01 /
   data green(18,11, 7) / 0.4498639443497940d-01 /
   data green(18,11, 8) / 0.4431751464690460d-01 /
   data green(18,11, 9) / 0.4359457578081609d-01 /
   data green(18,11,10) / 0.4282722282962292d-01 /
   data green(18,11,11) / 0.4202469429485825d-01 /
   data green(18,12, 0) / 0.4622343524631051d-01 /
   data green(18,12, 1) / 0.4617398050160367d-01 /
   data green(18,12, 2) / 0.4602656855687063d-01 /
   data green(18,12, 3) / 0.4578403726420281d-01 /
   data green(18,12, 4) / 0.4545088846600396d-01 /
   data green(18,12, 5) / 0.4503314566459500d-01 /
   data green(18,12, 6) / 0.4453799536394557d-01 /
   data green(18,12, 7) / 0.4397350663222339d-01 /
   data green(18,12, 8) / 0.4334827633485979d-01 /
   data green(18,12, 9) / 0.4267109349054533d-01 /
   data green(18,12,10) / 0.4195075944771441d-01 /
   data green(18,12,11) / 0.4119571396804787d-01 /
   data green(18,12,12) / 0.4041396773031416d-01 /
   data green(18,13, 0) / 0.4503488348780513d-01 /
   data green(18,13, 1) / 0.4498915477261323d-01 /
   data green(18,13, 2) / 0.4485281630779927d-01 /
   data green(18,13, 3) / 0.4462834202999478d-01 /
   data green(18,13, 4) / 0.4431967755687356d-01 /
   data green(18,13, 5) / 0.4393217496703627d-01 /
   data green(18,13, 6) / 0.4347218624599026d-01 /
   data green(18,13, 7) / 0.4294686425610071d-01 /
   data green(18,13, 8) / 0.4236391011185976d-01 /
   data green(18,13, 9) / 0.4173122603988050d-01 /
   data green(18,13,10) / 0.4105678361695903d-01 /
   data green(18,13,11) / 0.4034819040544022d-01 /
   data green(18,13,12) / 0.3961288350944592d-01 /
   data green(18,13,13) / 0.3885764661217822d-01 /
   data green(18,14, 0) / 0.4384926674740239d-01 /
   data green(18,14, 1) / 0.4380706336706453d-01 /
   data green(18,14, 2) / 0.4368118740717580d-01 /
   data green(18,14, 3) / 0.4347379928619551d-01 /
   data green(18,14, 4) / 0.4318840133285651d-01 /
   data green(18,14, 5) / 0.4282963185319920d-01 /
   data green(18,14, 6) / 0.4240315039681451d-01 /
   data green(18,14, 7) / 0.4191529945222977d-01 /
   data green(18,14, 8) / 0.4137284857841832d-01 /
   data green(18,14, 9) / 0.4078301621832606d-01 /
   data green(18,14,10) / 0.4015281668785910d-01 /
   data green(18,14,11) / 0.3948931152651212d-01 /
   data green(18,14,12) / 0.3879919685014722d-01 /
   data green(18,14,13) / 0.3808873572324255d-01 /
   data green(18,14,14) / 0.3736374873808740d-01 /
   data green(18,15, 0) / 0.4267493288689973d-01 /
   data green(18,15, 1) / 0.4263603488197849d-01 /
   data green(18,15, 2) / 0.4251997907426849d-01 /
   data green(18,15, 3) / 0.4232866310578393d-01 /
   data green(18,15, 4) / 0.4206514415076208d-01 /
   data green(18,15, 5) / 0.4173349966322484d-01 /
   data green(18,15, 6) / 0.4133866816238355d-01 /
   data green(18,15, 7) / 0.4088628341909525d-01 /
   data green(18,15, 8) / 0.4038240674535827d-01 /
   data green(18,15, 9) / 0.3983337328013360d-01 /
   data green(18,15,10) / 0.3924560541274894d-01 /
   data green(18,15,11) / 0.3862540662742558d-01 /
   data green(18,15,12) / 0.3797887367158370d-01 /
   data green(18,15,13) / 0.3731177517656965d-01 /
   data green(18,15,14) / 0.3662946155913213d-01 /
   data green(18,15,15) / 0.3593683038083195d-01 /
   data green(18,16, 0) / 0.4151861311816773d-01 /
   data green(18,16, 1) / 0.4148279599825051d-01 /
   data green(18,16, 2) / 0.4137590490814477d-01 /
   data green(18,16, 3) / 0.4119957562230880d-01 /
   data green(18,16, 4) / 0.4095649773039808d-01 /
   data green(18,16, 5) / 0.4065023471178629d-01 /
   data green(18,16, 6) / 0.4028514230088729d-01 /
   data green(18,16, 7) / 0.3986616092077289d-01 /
   data green(18,16, 8) / 0.3939867456434226d-01 /
   data green(18,16, 9) / 0.3888832745196970d-01 /
   data green(18,16,10) / 0.3834088005545602d-01 /
   data green(18,16,11) / 0.3776194969322479d-01 /
   data green(18,16,12) / 0.3715715658237222d-01 /
   data green(18,16,13) / 0.3653172595857034d-01 /
   data green(18,16,14) / 0.3589058956420365d-01 /
   data green(18,16,15) / 0.3523828368230370d-01 /
   data green(18,16,16) / 0.3457893824564939d-01 /
   data green(18,17, 0) / 0.4038560539778671d-01 /
   data green(18,17, 1) / 0.4035264350183123d-01 /
   data green(18,17, 2) / 0.4025424177733911d-01 /
   data green(18,17, 3) / 0.4009183464728781d-01 /
   data green(18,17, 4) / 0.3986776378749193d-01 /
   data green(18,17, 5) / 0.3958514609734598d-01 /
   data green(18,17, 6) / 0.3924779708184409d-01 /
   data green(18,17, 7) / 0.3886006952136728d-01 /
   data green(18,17, 8) / 0.3842674471295594d-01 /
   data green(18,17, 9) / 0.3795281331143387d-01 /
   data green(18,17,10) / 0.3744341078094923d-01 /
   data green(18,17,11) / 0.3690366495071960d-01 /
   data green(18,17,12) / 0.3633856674209219d-01 /
   data green(18,17,13) / 0.3575291896902526d-01 /
   data green(18,17,14) / 0.3515124639010405d-01 /
   data green(18,17,15) / 0.3453773683192799d-01 /
   data green(18,17,16) / 0.3391624149942888d-01 /
   data green(18,17,17) / 0.3329022528150614d-01 /
   data green(18,18, 0) / 0.3927995384830721d-01 /
   data green(18,18, 1) / 0.3924962682471758d-01 /
   data green(18,18, 2) / 0.3915909633156443d-01 /
   data green(18,18, 3) / 0.3900951982471178d-01 /
   data green(18,18, 4) / 0.3880304230685258d-01 /
   data green(18,18, 5) / 0.3854235002495100d-01 /
   data green(18,18, 6) / 0.3823076854020312d-01 /
   data green(18,18, 7) / 0.3787217365150290d-01 /
   data green(18,18, 8) / 0.3747073898955337d-01 /
   data green(18,18, 9) / 0.3703092123482025d-01 /
   data green(18,18,10) / 0.3655729878348245d-01 /
   data green(18,18,11) / 0.3605446665151452d-01 /
   data green(18,18,12) / 0.3552694092275994d-01 /
   data green(18,18,13) / 0.3497907049641263d-01 /
   data green(18,18,14) / 0.3441501856430099d-01 /
   data green(18,18,15) / 0.3383861205307182d-01 /
   data green(18,18,16) / 0.3325345650540754d-01 /
   data green(18,18,17) / 0.3266278259570786d-01 /
   data green(18,18,18) / 0.3206951241064816d-01 /

   data green(19, 0, 0) / 0.5266832259107634d-01 /
   data green(19, 1, 0) / 0.5259491421484098d-01 /
   data green(19, 1, 1) / 0.5252181141154750d-01 /
   data green(19, 2, 0) / 0.5237654127111076d-01 /
   data green(19, 2, 1) / 0.5230435953373048d-01 /
   data green(19, 2, 2) / 0.5208961514449295d-01 /
   data green(19, 3, 0) / 0.5201866958439184d-01 /
   data green(19, 3, 1) / 0.5194797715484865d-01 /
   data green(19, 3, 2) / 0.5173765465795430d-01 /
   data green(19, 3, 3) / 0.5139277601490361d-01 /
   data green(19, 4, 0) / 0.5152996172125972d-01 /
   data green(19, 4, 1) / 0.5146126110138179d-01 /
   data green(19, 4, 2) / 0.5125682893671041d-01 /
   data green(19, 4, 3) / 0.5092155726774003d-01 /
   data green(19, 4, 4) / 0.5046317116345187d-01 /
   data green(19, 5, 0) / 0.5092171416390186d-01 /
   data green(19, 5, 1) / 0.5085543925991835d-01 /
   data green(19, 5, 2) / 0.5065819554639128d-01 /
   data green(19, 5, 3) / 0.5033456857898344d-01 /
   data green(19, 5, 4) / 0.4989188520044252d-01 /
   data green(19, 5, 5) / 0.4933974472612889d-01 /
   data green(19, 6, 0) / 0.5020716654741713d-01 /
   data green(19, 6, 1) / 0.5014367147385186d-01 /
   data green(19, 6, 2) / 0.4995464468839594d-01 /
   data green(19, 6, 3) / 0.4964435518044589d-01 /
   data green(19, 6, 4) / 0.4921963960748020d-01 /
   data green(19, 6, 5) / 0.4868943687848241d-01 /
   data green(19, 6, 6) / 0.4806444043188928d-01 /
   data green(19, 7, 0) / 0.4940075306968572d-01 /
   data green(19, 7, 1) / 0.4934028972920873d-01 /
   data green(19, 7, 2) / 0.4916024895432470d-01 /
   data green(19, 7, 3) / 0.4886456902620757d-01 /
   data green(19, 7, 4) / 0.4845953841162737d-01 /
   data green(19, 7, 5) / 0.4795347996662469d-01 /
   data green(19, 7, 6) / 0.4735619234062338d-01 /
   data green(19, 7, 7) / 0.4667855830606328d-01 /
   data green(19, 8, 0) / 0.4851740100846755d-01 /
   data green(19, 8, 1) / 0.4846014603934912d-01 /
   data green(19, 8, 2) / 0.4828959922915432d-01 /
   data green(19, 8, 3) / 0.4800938518211196d-01 /
   data green(19, 8, 4) / 0.4762523652983048d-01 /
   data green(19, 8, 5) / 0.4714476180226954d-01 /
   data green(19, 8, 6) / 0.4657701267234445d-01 /
   data green(19, 8, 7) / 0.4593200746352408d-01 /
   data green(19, 8, 8) / 0.4522030240394082d-01 /
   data green(19, 9, 0) / 0.4757190884105121d-01 /
   data green(19, 9, 1) / 0.4751795462401835d-01 /
   data green(19, 9, 2) / 0.4735720559680244d-01 /
   data green(19, 9, 3) / 0.4709295264614891d-01 /
   data green(19, 9, 4) / 0.4673032098668063d-01 /
   data green(19, 9, 5) / 0.4627633427950465d-01 /
   data green(19, 9, 6) / 0.4573915142388196d-01 /
   data green(19, 9, 7) / 0.4512803463558859d-01 /
   data green(19, 9, 8) / 0.4445265405858352d-01 /
   data green(19, 9, 9) / 0.4372294528620843d-01 /
   data green(19,10, 0) / 0.4657844516406008d-01 /
   data green(19,10, 1) / 0.4652781551983176d-01 /
   data green(19,10, 2) / 0.4637692397445072d-01 /
   data green(19,10, 3) / 0.4612872110179154d-01 /
   data green(19,10, 4) / 0.4578788987203487d-01 /
   data green(19,10, 5) / 0.4536068423956008d-01 /
   data green(19,10, 6) / 0.4485453071553357d-01 /
   data green(19,10, 7) / 0.4427791104284514d-01 /
   data green(19,10, 8) / 0.4363956602104892d-01 /
   data green(19,10, 9) / 0.4294895515199568d-01 /
   data green(19,10,10) / 0.4221423864193425d-01 /
   data green(19,11, 0) / 0.4555017246376767d-01 /
   data green(19,11, 1) / 0.4550283530708807d-01 /
   data green(19,11, 2) / 0.4536172483782476d-01 /
   data green(19,11, 3) / 0.4512944673555686d-01 /
   data green(19,11, 4) / 0.4481023839276159d-01 /
   data green(19,11, 5) / 0.4440968711098475d-01 /
   data green(19,11, 6) / 0.4393450858472005d-01 /
   data green(19,11, 7) / 0.4339231573376540d-01 /
   data green(19,11, 8) / 0.4279113022968248d-01 /
   data green(19,11, 9) / 0.4213922762251332d-01 /
   data green(19,11,10) / 0.4144495396003144d-01 /
   data green(19,11,11) / 0.4071630472646075d-01 /
   data green(19,12, 0) / 0.4449900839847521d-01 /
   data green(19,12, 1) / 0.4445488517798557d-01 /
   data green(19,12, 2) / 0.4432330611585193d-01 /
   data green(19,12, 3) / 0.4410664735747223d-01 /
   data green(19,12, 4) / 0.4380855162163242d-01 /
   data green(19,12, 5) / 0.4343410794754284d-01 /
   data green(19,12, 6) / 0.4298938355303930d-01 /
   data green(19,12, 7) / 0.4248111062616419d-01 /
   data green(19,12, 8) / 0.4191658717626175d-01 /
   data green(19,12, 9) / 0.4130337700512753d-01 /
   data green(19,12,10) / 0.4064903950840235d-01 /
   data green(19,12,11) / 0.3996093975891628d-01 /
   data green(19,12,12) / 0.3924611316472242d-01 /
   data green(19,13, 0) / 0.4343548564985747d-01 /
   data green(19,13, 1) / 0.4339445870439769d-01 /
   data green(19,13, 2) / 0.4327207523670244d-01 /
   data green(19,13, 3) / 0.4307041246832503d-01 /
   data green(19,13, 4) / 0.4279280779371831d-01 /
   data green(19,13, 5) / 0.4244368582609741d-01 /
   data green(19,13, 6) / 0.4202847402494354d-01 /
   data green(19,13, 7) / 0.4155323437898869d-01 /
   data green(19,13, 8) / 0.4102451331790256d-01 /
   data green(19,13, 9) / 0.4044916040641412d-01 /
   data green(19,13,10) / 0.3983404691600849d-01 /
   data green(19,13,11) / 0.3918589148727772d-01 /
   data green(19,13,12) / 0.3851123818506239d-01 /
   data green(19,13,13) / 0.3781615625989714d-01 /
   data green(19,14, 0) / 0.4236872966790754d-01 /
   data green(19,14, 1) / 0.4233064751481610d-01 /
   data green(19,14, 2) / 0.4221705447194933d-01 /
   data green(19,14, 3) / 0.4202976060216406d-01 /
   data green(19,14, 4) / 0.4177172080936307d-01 /
   data green(19,14, 5) / 0.4144689028187904d-01 /
   data green(19,14, 6) / 0.4106003445862264d-01 /
   data green(19,14, 7) / 0.4061660376605981d-01 /
   data green(19,14, 8) / 0.4012246043312818d-01 /
   data green(19,14, 9) / 0.3958380940815934d-01 /
   data green(19,14,10) / 0.3900680287445832d-01 /
   data green(19,14,11) / 0.3839763737515740d-01 /
   data green(19,14,12) / 0.3776225014439741d-01 /
   data green(19,14,13) / 0.3710626322907910d-01 /
   data green(19,14,14) / 0.3643491837726490d-01 /
   data green(19,15, 0) / 0.4130643091292196d-01 /
   data green(19,15, 1) / 0.4127115670238986d-01 /
   data green(19,15, 2) / 0.4116588941357448d-01 /
   data green(19,15, 3) / 0.4099220534395089d-01 /
   data green(19,15, 4) / 0.4075273694866115d-01 /
   data green(19,15, 5) / 0.4045097708752012d-01 /
   data green(19,15, 6) / 0.4009115561081315d-01 /
   data green(19,15, 7) / 0.3967811911403152d-01 /
   data green(19,15, 8) / 0.3921709660198509d-01 /
   data green(19,15, 9) / 0.3871367213811946d-01 /
   data green(19,15,10) / 0.3817341531833542d-01 /
   data green(19,15,11) / 0.3760192188923826d-01 /
   data green(19,15,12) / 0.3700463168160356d-01 /
   data green(19,15,13) / 0.3638671464532157d-01 /
   data green(19,15,14) / 0.3575302474495262d-01 /
   data green(19,15,15) / 0.3510804455738082d-01 /
   data green(19,16, 0) / 0.4025502838463846d-01 /
   data green(19,16, 1) / 0.4022238295630504d-01 /
   data green(19,16, 2) / 0.4012492423543754d-01 /
   data green(19,16, 3) / 0.3996407217619601d-01 /
   data green(19,16, 4) / 0.3974211206792933d-01 /
   data green(19,16, 5) / 0.3946213556682676d-01 /
   data green(19,16, 6) / 0.3912787496101781d-01 /
   data green(19,16, 7) / 0.3874366238307352d-01 /
   data green(19,16, 8) / 0.3831414400324829d-01 /
   data green(19,16, 9) / 0.3784430105776799d-01 /
   data green(19,16,10) / 0.3733918083659542d-01 /
   data green(19,16,11) / 0.3680384535975049d-01 /
   data green(19,16,12) / 0.3624323422463845d-01 /
   data green(19,16,13) / 0.3566209786569426d-01 /
   data green(19,16,14) / 0.3506490501636884d-01 /
   data green(19,16,15) / 0.3445581935113427d-01 /
   data green(19,16,16) / 0.3383863944261860d-01 /
   data green(19,17, 0) / 0.3921971014722767d-01 /
   data green(19,17, 1) / 0.3918952895475911d-01 /
   data green(19,17, 2) / 0.3909939162776711d-01 /
   data green(19,17, 3) / 0.3895051983504493d-01 /
   data green(19,17, 4) / 0.3874496539645023d-01 /
   data green(19,17, 5) / 0.3848542416553417d-01 /
   data green(19,17, 6) / 0.3817521572884985d-01 /
   data green(19,17, 7) / 0.3781815446919554d-01 /
   data green(19,17, 8) / 0.3741840767162204d-01 /
   data green(19,17, 9) / 0.3697797254949323d-01 /
   data green(19,17,10) / 0.3650868343984736d-01 /
   data green(19,17,11) / 0.3600782492238904d-01 /
   data green(19,17,12) / 0.3548230118330761d-01 /
   data green(19,17,13) / 0.3493645764644598d-01 /
   data green(19,17,14) / 0.3437441971494789d-01 /
   data green(19,17,15) / 0.3380002037582742d-01 /
   data green(19,17,16) / 0.3321681993373137d-01 /
   data green(19,17,17) / 0.3262805605118717d-01 /
   data green(19,18, 0) / 0.3820462799260847d-01 /
   data green(19,18, 1) / 0.3817672550555491d-01 /
   data green(19,18, 2) / 0.3809338487900765d-01 /
   data green(19,18, 3) / 0.3795569232774473d-01 /
   data green(19,18, 4) / 0.3776543449210863d-01 /
   data green(19,18, 5) / 0.3752497054617601d-01 /
   data green(19,18, 6) / 0.3723727278237021d-01 /
   data green(19,18, 7) / 0.3690568702804106d-01 /
   data green(19,18, 8) / 0.3653392223363360d-01 /
   data green(19,18, 9) / 0.3612593219182070d-01 /
   data green(19,18,10) / 0.3568579023487571d-01 /
   data green(19,18,11) / 0.3521760986997396d-01 /
   data green(19,18,12) / 0.3472547847357045d-01 /
   data green(19,18,13) / 0.3421332657378441d-01 /
   data green(19,18,14) / 0.3368494121744788d-01 /
   data green(19,18,15) / 0.3314387509937145d-01 /
   data green(19,18,16) / 0.3259343200314951d-01 /
   data green(19,18,17) / 0.3203663703797065d-01 /
   data green(19,18,18) / 0.3147624220909137d-01 /
   data green(19,19, 0) / 0.3721295036200426d-01 /
   data green(19,19, 1) / 0.3718716318741537d-01 /
   data green(19,19, 2) / 0.3711018983789283d-01 /
   data green(19,19, 3) / 0.3698281245220978d-01 /
   data green(19,19, 4) / 0.3680673859429633d-01 /
   data green(19,19, 5) / 0.3658406402949283d-01 /
   data green(19,19, 6) / 0.3631734094277288d-01 /
   data green(19,19, 7) / 0.3600951818423553d-01 /
   data green(19,19, 8) / 0.3566396021399994d-01 /
   data green(19,19, 9) / 0.3528413680251755d-01 /
   data green(19,19,10) / 0.3487371493720784d-01 /
   data green(19,19,11) / 0.3443639993880904d-01 /
   data green(19,19,12) / 0.3397586723418568d-01 /
   data green(19,19,13) / 0.3349571681473099d-01 /
   data green(19,19,14) / 0.3299941083948185d-01 /
   data green(19,19,15) / 0.3249019659806088d-01 /
   data green(19,19,16) / 0.3197117995205517d-01 /
   data green(19,19,17) / 0.3144515151521055d-01 /
   data green(19,19,18) / 0.3091470713645940d-01 /
   data green(19,19,19) / 0.3038217784602753d-01 /

   data green(20, 0, 0) / 0.5003147782279721d-01 /
   data green(20, 1, 0) / 0.4996857965445155d-01 /
   data green(20, 1, 1) / 0.4990592247536536d-01 /
   data green(20, 2, 0) / 0.4978132851137523d-01 /
   data green(20, 2, 1) / 0.4971937925735290d-01 /
   data green(20, 2, 2) / 0.4953494008792210d-01 /
   data green(20, 3, 0) / 0.4947394504959887d-01 /
   data green(20, 3, 1) / 0.4941314671138683d-01 /
   data green(20, 3, 2) / 0.4923212330593416d-01 /
   data green(20, 3, 3) / 0.4893484925909421d-01 /
   data green(20, 4, 0) / 0.4905315836873730d-01 /
   data green(20, 4, 1) / 0.4899391684095938d-01 /
   data green(20, 4, 2) / 0.4881749294172963d-01 /
   data green(20, 4, 3) / 0.4852768919278038d-01 /
   data green(20, 4, 4) / 0.4813059090258144d-01 /
   data green(20, 5, 0) / 0.4852782635136227d-01 /
   data green(20, 5, 1) / 0.4847048599323722d-01 /
   data green(20, 5, 2) / 0.4829966476323132d-01 /
   data green(20, 5, 3) / 0.4801906131872140d-01 /
   data green(20, 5, 4) / 0.4763432615989711d-01 /
   data green(20, 5, 5) / 0.4715328238449140d-01 /
   data green(20, 6, 0) / 0.4790840534122267d-01 /
   data green(20, 6, 1) / 0.4785325062817174d-01 /
   data green(20, 6, 2) / 0.4768874689471643d-01 /
   data green(20, 6, 3) / 0.4741885631163341d-01 /
   data green(20, 6, 4) / 0.4704843283701984d-01 /
   data green(20, 6, 5) / 0.4658482952280266d-01 /
   data green(20, 6, 6) / 0.4603659701058205d-01 /
   data green(20, 7, 0) / 0.4720642965826296d-01 /
   data green(20, 7, 1) / 0.4715368498112246d-01 /
   data green(20, 7, 2) / 0.4699651017203263d-01 /
   data green(20, 7, 3) / 0.4673803322381157d-01 /
   data green(20, 7, 4) / 0.4638334985815126d-01 /
   data green(20, 7, 5) / 0.4593908580919895d-01 /
   data green(20, 7, 6) / 0.4541321976167775d-01 /
   data green(20, 7, 7) / 0.4481464929706800d-01 /
   data green(20, 8, 0) / 0.4643397221710913d-01 /
   data green(20, 8, 1) / 0.4638378786640200d-01 /
   data green(20, 8, 2) / 0.4623421007879004d-01 /
   data green(20, 8, 3) / 0.4598814666422973d-01 /
   data green(20, 8, 4) / 0.4565023020834109d-01 /
   data green(20, 8, 5) / 0.4522666656359477d-01 /
   data green(20, 8, 6) / 0.4472474834059118d-01 /
   data green(20, 8, 7) / 0.4415278101997594d-01 /
   data green(20, 8, 8) / 0.4351952287220272d-01 /
   data green(20, 9, 0) / 0.4560317119147676d-01 /
   data green(20, 9, 1) / 0.4555564911061586d-01 /
   data green(20, 9, 2) / 0.4541397233902599d-01 /
   data green(20, 9, 3) / 0.4518079199214250d-01 /
   data green(20, 9, 4) / 0.4486034621674326d-01 /
   data green(20, 9, 5) / 0.4445829053634860d-01 /
   data green(20, 9, 6) / 0.4398139912901106d-01 /
   data green(20, 9, 7) / 0.4343725517858835d-01 /
   data green(20, 9, 8) / 0.4283399236086775d-01 /
   data green(20, 9, 9) / 0.4217997689855560d-01 /
   data green(20,10, 0) / 0.4472582547950323d-01 /
   data green(20,10, 1) / 0.4468100461301278d-01 /
   data green(20,10, 2) / 0.4454735393196586d-01 /
   data green(20,10, 3) / 0.4432727663044406d-01 /
   data green(20,10, 4) / 0.4402463372808131d-01 /
   data green(20,10, 5) / 0.4364452957916899d-01 /
   data green(20,10, 6) / 0.4319320182886031d-01 /
   data green(20,10, 7) / 0.4267755656728302d-01 /
   data green(20,10, 8) / 0.4210513323283776d-01 /
   data green(20,10, 9) / 0.4148351651003902d-01 /
   data green(20,10,10) / 0.4082054116507249d-01 /
   data green(20,11, 0) / 0.4381307546249753d-01 /
   data green(20,11, 1) / 0.4377095371101127d-01 /
   data green(20,11, 2) / 0.4364531879091178d-01 /
   data green(20,11, 3) / 0.4343833937353846d-01 /
   data green(20,11, 4) / 0.4315352063348728d-01 /
   data green(20,11, 5) / 0.4279542870015960d-01 /
   data green(20,11, 6) / 0.4236976369881548d-01 /
   data green(20,11, 7) / 0.4188282822856045d-01 /
   data green(20,11, 8) / 0.4134144164340638d-01 /
   data green(20,11, 9) / 0.4075269368545795d-01 /
   data green(20,11,10) / 0.4012369279918613d-01 /
   data green(20,11,11) / 0.3946140600066386d-01 /
   data green(20,12, 0) / 0.4287518153769498d-01 /
   data green(20,12, 1) / 0.4283567480856715d-01 /
   data green(20,12, 2) / 0.4271799653034785d-01 /
   data green(20,12, 3) / 0.4252390298167427d-01 /
   data green(20,12, 4) / 0.4225663470799169d-01 /
   data green(20,12, 5) / 0.4192033167626557d-01 /
   data green(20,12, 6) / 0.4152008220526816d-01 /
   data green(20,12, 7) / 0.4106162536119000d-01 /
   data green(20,12, 8) / 0.4055116066255539d-01 /
   data green(20,12, 9) / 0.3999514241291746d-01 /
   data green(20,12,10) / 0.3940011953799236d-01 /
   data green(20,12,11) / 0.3877254677151808d-01 /
   data green(20,12,12) / 0.3811863129385605d-01 /
   data green(20,13, 0) / 0.4192133669216255d-01 /
   data green(20,13, 1) / 0.4188450917347730d-01 /
   data green(20,13, 2) / 0.4177444701991690d-01 /
   data green(20,13, 3) / 0.4159292964307606d-01 /
   data green(20,13, 4) / 0.4134277675682798d-01 /
   data green(20,13, 5) / 0.4102760307620792d-01 /
   data green(20,13, 6) / 0.4065236349316348d-01 /
   data green(20,13, 7) / 0.4022180987765270d-01 /
   data green(20,13, 8) / 0.3974173667864837d-01 /
   data green(20,13, 9) / 0.3921799805143671d-01 /
   data green(20,13,10) / 0.3865653184023669d-01 /
   data green(20,13,11) / 0.3806328468245283d-01 /
   data green(20,13,12) / 0.3744413328067516d-01 /
   data green(20,13,13) / 0.3680422099751433d-01 /
   data green(20,14, 0) / 0.4095983579882056d-01 /
   data green(20,14, 1) / 0.4092543862131844d-01 /
   data green(20,14, 2) / 0.4082276289584826d-01 /
   data green(20,14, 3) / 0.4065336782003948d-01 /
   data green(20,14, 4) / 0.4041973677170600d-01 /
   data green(20,14, 5) / 0.4012518625800629d-01 /
   data green(20,14, 6) / 0.3977392197957901d-01 /
   data green(20,14, 7) / 0.3937048701491352d-01 /
   data green(20,14, 8) / 0.3891995260739173d-01 /
   data green(20,14, 9) / 0.3842768894409253d-01 /
   data green(20,14,10) / 0.3789911150927097d-01 /
   data green(20,14,11) / 0.3733963496263833d-01 /
   data green(20,14,12) / 0.3675451566370086d-01 /
   data green(20,14,13) / 0.3614879874174875d-01 /
   data green(20,14,14) / 0.3552720322599993d-01 /
   data green(20,15, 0) / 0.3999758949423679d-01 /
   data green(20,15, 1) / 0.3996556510143659d-01 /
   data green(20,15, 2) / 0.3986995957713585d-01 /
   data green(20,15, 3) / 0.3971212318823773d-01 /
   data green(20,15, 4) / 0.3949428519624955d-01 /
   data green(20,15, 5) / 0.3921944343829193d-01 /
   data green(20,15, 6) / 0.3889126891956751d-01 /
   data green(20,15, 7) / 0.3851387758554787d-01 /
   data green(20,15, 8) / 0.3809185179989451d-01 /
   data green(20,15, 9) / 0.3763001978464734d-01 /
   data green(20,15,10) / 0.3713329300515899d-01 /
   data green(20,15,11) / 0.3660660411591864d-01 /
   data green(20,15,12) / 0.3605479618631158d-01 /
   data green(20,15,13) / 0.3548250152044416d-01 /
   data green(20,15,14) / 0.3489411070849359d-01 /
   data green(20,15,15) / 0.3429369651071714d-01 /
   data green(20,16, 0) / 0.3904064107594118d-01 /
   data green(20,16, 1) / 0.3901086315005439d-01 /
   data green(20,16, 2) / 0.3892192474230244d-01 /
   data green(20,16, 3) / 0.3877507997659928d-01 /
   data green(20,16, 4) / 0.3857226834189545d-01 /
   data green(20,16, 5) / 0.3831613143862127d-01 /
   data green(20,16, 6) / 0.3800996336589046d-01 /
   data green(20,16, 7) / 0.3765745811863105d-01 /
   data green(20,16, 8) / 0.3726268225634342d-01 /
   data green(20,16, 9) / 0.3683007903856923d-01 /
   data green(20,16,10) / 0.3636399275887506d-01 /
   data green(20,16,11) / 0.3586896240633952d-01 /
   data green(20,16,12) / 0.3534938852627099d-01 /
   data green(20,16,13) / 0.3480954650034342d-01 /
   data green(20,16,14) / 0.3425348848066810d-01 /
   data green(20,16,15) / 0.3368502359590180d-01 /
   data green(20,16,16) / 0.3310763879633136d-01 /
   data green(20,17, 0) / 0.3809398272830940d-01 /
   data green(20,17, 1) / 0.3806632104778947d-01 /
   data green(20,17, 2) / 0.3798369806461166d-01 /
   data green(20,17, 3) / 0.3784718961162297d-01 /
   data green(20,17, 4) / 0.3765851969143245d-01 /
   data green(20,17, 5) / 0.3742010292771095d-01 /
   data green(20,17, 6) / 0.3713477761555953d-01 /
   data green(20,17, 7) / 0.3680588493584639d-01 /
   data green(20,17, 8) / 0.3643708857966097d-01 /
   data green(20,17, 9) / 0.3603228912714095d-01 /
   data green(20,17,10) / 0.3559550588748400d-01 /
   data green(20,17,11) / 0.3513081237666212d-01 /
   data green(20,17,12) / 0.3464224128988944d-01 /
   data green(20,17,13) / 0.3413370195484366d-01 /
   data green(20,17,14) / 0.3360892975252311d-01 /
   data green(20,17,15) / 0.3307143961047686d-01 /
   data green(20,17,16) / 0.3252453510631315d-01 /
   data green(20,17,17) / 0.3197119789894768d-01 /
   data green(20,18, 0) / 0.3716169895451119d-01 /
   data green(20,18, 1) / 0.3713602039660262d-01 /
   data green(20,18, 2) / 0.3705930419318157d-01 /
   data green(20,18, 3) / 0.3693249226034289d-01 /
   data green(20,18, 4) / 0.3675715278262096d-01 /
   data green(20,18, 5) / 0.3653536830426567d-01 /
   data green(20,18, 6) / 0.3626967751379730d-01 /
   data green(20,18, 7) / 0.3596307078844869d-01 /
   data green(20,18, 8) / 0.3561880703236430d-01 /
   data green(20,18, 9) / 0.3524041516736324d-01 /
   data green(20,18,10) / 0.3483149613796071d-01 /
   data green(20,18,11) / 0.3439575084062178d-01 /
   data green(20,18,12) / 0.3393685752837700d-01 /
   data green(20,18,13) / 0.3345828814020074d-01 /
   data green(20,18,14) / 0.3296361447729214d-01 /
   data green(20,18,15) / 0.3245603972135603d-01 /
   data green(20,18,16) / 0.3193861576343361d-01 /
   data green(20,18,17) / 0.3141416095142756d-01 /
   data green(20,18,18) / 0.3088525276086289d-01 /
   data green(20,19, 0) / 0.3624705523409605d-01 /
   data green(20,19, 1) / 0.3622322755717346d-01 /
   data green(20,19, 2) / 0.3615202421731307d-01 /
   data green(20,19, 3) / 0.3603429045111626d-01 /
   data green(20,19, 4) / 0.3587137851187712d-01 /
   data green(20,19, 5) / 0.3566517035636952d-01 /
   data green(20,19, 6) / 0.3541790556293278d-01 /
   data green(20,19, 7) / 0.3513226077115970d-01 /
   data green(20,19, 8) / 0.3481111387823038d-01 /
   data green(20,19, 9) / 0.3445764063397681d-01 /
   data green(20,19,10) / 0.3407508698676612d-01 /
   data green(20,19,11) / 0.3366678808399975d-01 /
   data green(20,19,12) / 0.3323608916713310d-01 /
   data green(20,19,13) / 0.3278619867827280d-01 /
   data green(20,19,14) / 0.3232031894378669d-01 /
   data green(20,19,15) / 0.3184145395069676d-01 /
   data green(20,19,16) / 0.3135242257876048d-01 /
   data green(20,19,17) / 0.3085587359135744d-01 /
   data green(20,19,18) / 0.3035421417126546d-01 /
   data green(20,19,19) / 0.2984964849460441d-01 /
   data green(20,20, 0) / 0.3535259625594454d-01 /
   data green(20,20, 1) / 0.3533049032276067d-01 /
   data green(20,20, 2) / 0.3526438293584833d-01 /
   data green(20,20, 3) / 0.3515515101923312d-01 /
   data green(20,20, 4) / 0.3500382602864701d-01 /
   data green(20,20, 5) / 0.3481213994405296d-01 /
   data green(20,20, 6) / 0.3458210578320652d-01 /
   data green(20,20, 7) / 0.3431606785956563d-01 /
   data green(20,20, 8) / 0.3401660733136788d-01 /
   data green(20,20, 9) / 0.3368657119914283d-01 /
   data green(20,20,10) / 0.3332887598261000d-01 /
   data green(20,20,11) / 0.3294651975468069d-01 /
   data green(20,20,12) / 0.3254255301965124d-01 /
   data green(20,20,13) / 0.3211989650126146d-01 /
   data green(20,20,14) / 0.3168146640461020d-01 /
   data green(20,20,15) / 0.3123003434764223d-01 /
   data green(20,20,16) / 0.3076825488954206d-01 /
   data green(20,20,17) / 0.3029851145532374d-01 /
   data green(20,20,18) / 0.2982312486803485d-01 /
   data green(20,20,19) / 0.2934416045219572d-01 /
   data green(20,20,20) / 0.2886350273699321d-01 /

end block data blkgreen
