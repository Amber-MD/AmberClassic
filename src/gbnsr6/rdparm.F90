#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ A lite version of Amber rdparm1
subroutine rdparm1(nf)
   
   implicit none
   
#  include "md.h"
#  include "memory.h"
#  include "parms.h"
#  include "files.h"
#  include "box.h"
   integer nf
   integer iok!,nspsol
   integer nhparm,nttyp!,idum
   integer mbper,mgper,mdper,mbona,mtheta,mphia ! read but ignored
   integer numextra
  ! integer ifpert,nbper,ngper,ndper! local, not going anywhere
   character(len=80) fmt
   character(len=80) fmtin,ifmt,afmt,rfmt,type
   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   
   !     ----- READ THE MOLECULAR TOPOLOGY -----
   
   !nspsol = 0
   
   !     ----- FORMATTED INPUT -----
   
   fmtin = '(A80)'
   type = 'TITLE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmtin) title
   
   fmtin = ifmt
   type = 'POINTERS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
         numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper, &
         mbper,mgper,mdper,ifbox,nmxrs,ifcap,numextra,ncopy
   if (nparm == 1) then
      write(6,'(a)') ' *** THIS VERSION WILL NOT ACCEPT TOPOLOGY FILES'
      write(6,'(a)') '     THAT WERE CREATED BY ADDLES, WITH NPARM=1'
      write(6,'(a)') '     USE A VERSION COMPILED WITH -DLES '
      call mexit(6,1)
   end if
   write(6,8118) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia,numbnd, &
         numang,nptra,natyp,nphb,ifbox,nmxrs,ifcap,numextra,ncopy
   8118 format(t2, &
         'NATOM  = ',i7,' NTYPES = ',i7,' NBONH = ',i7,' MBONA  = ',i7, &
         /' NTHETH = ',i7,' MTHETA = ',i7,' NPHIH = ',i7,' MPHIA  = ',i7, &
         /' NHPARM = ',i7,' NPARM  = ',i7,' NNB   = ',i7,' NRES   = ',i7, &
         /' NBONA  = ',i7,' NTHETA = ',i7,' NPHIA = ',i7,' NUMBND = ',i7, &
         /' NUMANG = ',i7,' NPTRA  = ',i7,' NATYP = ',i7,' NPHB   = ',i7, &
         /' IFBOX  = ',i7,' NMXRS  = ',i7,' IFCAP = ',i7,' NEXTRA = ',i7 &
         ,/' NCOPY  = ',i7/)
   
   ! make sure we do not exceed memory limits in commons
   nttyp = ntypes*(ntypes+1)/2
   if (numbnd > 5000 .or. numang > 900 .or. nptra > 1200 .or. &
         nphb > 200 .or. natyp > 60 .or. nttyp > 1830) then
      write(6,'(/,5x,a)') 'rdparm: a parameter array overflowed'
      write(6,'(/,5x,a)') '       (e.g. the table of dihedral params)'
      call mexit(6, 1)
   end if
   
   if(nbona /= mbona .or. ntheta /= mtheta .or. nphia /= mphia) then
      write(6,'(a)') 'Sander no longer allows constraints in prmtop'
      write(6,'(a)') '...must have nbona=mbona, ntheta=mtheta, nphi=mphi'
      call mexit(6,1)
   end if
   
   return
end subroutine rdparm1 

!=======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ A lite version of Amber rdparm2
subroutine rdparm2(x,ix,ih,ipairs,nf,i_stack)
   implicit none
   integer nf
   _REAL_ x(*)
   integer ix(*),ipairs(*),i_stack(*)
   character(len=4) ih(*)
   
#  include "md.h"
#  include "memory.h"
#  include "parms.h"
#  include "files.h"
#  include "box.h"

   ! Local variables
   integer nttyp,ntype,i,iok
   integer j,jj,k,l,n,nn
   integer iptres,nspsol,natsm!,idum,ip14
   _REAL_ dumd,oldbeta,duma,dumb,dumc

   character(len=80) fmt
   character(len=80) fmtin,ifmt,afmt,rfmt,type
   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   nttyp = ntypes*(ntypes+1)/2
   ntype = ntypes*ntypes
   
   !     ----- READ THE SYMBOLS AND THE CHARGES AND THE MASSES -----
   
   fmtin = afmt
   type = 'ATOM_NAME'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ih(m04+i-1),i = 1,natom)
   
   fmtin = rfmt
   type = 'CHARGE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(l15+i-1),i = 1,natom)
   
   fmtin = rfmt
   type = 'MASS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(lwinv+i-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'ATOM_TYPE_INDEX'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i04+i-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'NUMBER_EXCLUDED_ATOMS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i08-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'NONBONDED_PARM_INDEX'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i06-1),i = 1,ntype)
   
   fmtin = afmt
   type = 'RESIDUE_LABEL'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ih(i+m02-1),i=1,nres)
   
   fmtin = ifmt
   type = 'RESIDUE_POINTER'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i02-1),i=1,nres)
   ix(i02+nres) = natom+1
   
   !     ----- READ THE PARAMETERS -----
   
   fmtin = rfmt
   type = 'BOND_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (rk(i),    i = 1,numbnd)
   
   fmtin = rfmt
   type = 'BOND_EQUIL_VALUE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (req(i),   i = 1,numbnd)
   
   fmtin = rfmt
   type = 'ANGLE_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (tk(i),    i = 1,numang)
   
   fmtin = rfmt
   type = 'ANGLE_EQUIL_VALUE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (teq(i),   i = 1,numang)
   
   fmtin = rfmt
   type = 'DIHEDRAL_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (pk(i),    i = 1,nptra)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PERIODICITY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (pn(i),    i = 1,nptra)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PHASE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (phase(i), i = 1,nptra)

   !   ----- READ VARIABLE SCEE -----
   ! --- SCEE ---
   fmtin = rfmt
   type = 'SCEE_SCALE_FACTOR'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   if (iok == 0) then
      !we found the SCEE scale factor data so read it in, there should be
      !a value for each dihedral type. Note while there is one for each
      !dihedral type not all are actually used in the calculation since
      !1-4's are only done over unique dihedral terms. If multiple dihedrals
      !for a given set of atom types exist, with different pn's for example
      !then the last value of SCEE/SCNB should be used.
      write (6,'(a)') ''
      write (6,'(a)') '| Note: 1-4 EEL scale factors are being read from the topology file.'

      read(nf,fmt) (one_scee(i), i = 1,nptra)
      one_scee(1:nptra) = 1.0d0  / one_scee(1:nptra)
   else
      !We will use default scee of 1.2
      write (6,'(a)') ''
      write (6,'(a)') '| Note: 1-4 EEL scale factors were NOT found in the topology file.'
      write (6,'(a)') '|       Using default value of 1.2.'
      do i=1,nptra
        one_scee(i)=1.0d0/1.2d0
      end do
   end if

   
   fmtin = rfmt
   type = 'SOLTY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (solty(i), i = 1,natyp)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_ACOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (cn1(i),   i = 1,nttyp)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_BCOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (cn2(i),   i = 1,nttyp)
   
   !     ----- READ THE BONDING INFORMATIONS -----
   
   fmtin = ifmt
   type = 'BONDS_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+iibh-1),ix(i+ijbh-1),ix(i+iicbh-1), &
         i = 1,nbonh)
   
   fmtin = ifmt
   type = 'BONDS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt)(ix(i+iiba-1),ix(i+ijba-1),ix(i+iicba-1),i = 1,nbona)
   
   fmtin = ifmt
   type = 'ANGLES_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i24-1),ix(i+i26-1),ix(i+i28-1),ix(i+i30-1), &
         i = 1,ntheth)
   
   fmtin = ifmt
   type = 'ANGLES_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i32-1),ix(i+i34-1),ix(i+i36-1),ix(i+i38-1), &
         i = 1,ntheta)
   
   fmtin = ifmt
   type = 'DIHEDRALS_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i40-1),ix(i+i42-1),ix(i+i44-1),ix(i+i46-1), &
         ix(i+i48-1),i = 1,nphih)
   
   fmtin = ifmt
   type = 'DIHEDRALS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i50-1),ix(i+i52-1),ix(i+i54-1),ix(i+i56-1), &
         ix(i+i58-1),i = 1,nphia)
   
   fmtin = ifmt
   type = 'EXCLUDED_ATOMS_LIST'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i10-1),i=1,nnb)
   
   !     ----- READ THE H-BOND PARAMETERS -----
   
   fmtin = rfmt
   type = 'HBOND_ACOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (asol(i),i=1,nphb)
   
#ifndef HAS_10_12
   do i=1,nphb
      if( asol(i) /= 0.d0 ) then
         write(6,'(a)') 'Found a non-zero 10-12 coefficient, but source', &
               ' was not compiled with -DHAS_10_12.'
         write(6,'(a)') 'If you are using a pre-1994 force field, you', &
               ' will need to re-compile with this flag.'
         call mexit(6,1)
      end if
   end do
#endif
   
   fmtin = rfmt
   type = 'HBOND_BCOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (bsol(i),i=1,nphb)
   
   fmtin = rfmt
   type = 'HBCUT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (hbcut(i),i=1,nphb)
   
   !     ----- READ ISYMBL,ITREE,JOIN AND IROTAT ARRAYS -----
   
   fmtin = afmt
   type = 'AMBER_ATOM_TYPE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ih(i+m06-1),i=1,natom)
   
!  fmtin = afmt
!  type = 'TREE_CHAIN_CLASSIFICATION'
!  call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
!  read(nf,fmt) (ih(i+m08-1),i=1,natom)
   
!  fmtin = ifmt
!  type = 'JOIN_ARRAY'
!  call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
!  read(nf,fmt) (ix(i+i64-1),i=1,natom)
   
!  fmtin = ifmt
!  type = 'IROTAT'
!  call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
!  read(nf,fmt) (ix(i+i66-1),i=1,natom)
   
   !     ----- READ THE BOUNDARY CONDITION STUFF -----
   
   nspm = 1
   ix(i70) = natom
   if (ifbox > 0) then
      
      fmtin = ifmt
      type = 'SOLVENT_POINTERS'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) iptres,nspm,nspsol
      
      fmtin = ifmt
      type = 'ATOMS_PER_MOLECULE'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ix(i+i70-1),i=1,nspm)
      
      fmtin = rfmt
      type = 'BOX_DIMENSIONS'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) oldbeta,duma,dumb,dumc
      
      !       ---(above values are still read, for backward compatibility, but
      !           ignored. Box info must come from the coord. file or from
      !           the &ewald namelist of the input file)
      
      if( ipb >= 1  .or.  ntb == 0 )then
         box(1)=0.0d0
         box(2)=0.0d0
         box(3)=0.0d0
      end if
      
   end if  ! (ifbox > 0)
   
   !     ----- LOAD THE CAP INFORMATION IF NEEDED -----
   
   if(ifcap > 0) then
      fmtin = '(I6)'
      type = 'CAP_INFO'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) natcap
      
      fmtin = '(4E16.8)'
      type = 'CAP_INFO2'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) cutcap,xcap,ycap,zcap
   else
      natcap = 0
      cutcap = 0.0d0
      xcap = 0.0d0
      ycap = 0.0d0
      zcap = 0.0d0
   end if
   
   if( ipb >= 1 .and. iok == -1 ) then
      write(6,'(a)') 'GB calculations now require a new-style prmtop file'
      write(6,'(a)') 'GB calculations now require a new-style prmtop file'
      call mexit(6,1)
   end if
   
   if ( ipb >= 1 ) then
      fmtin = rfmt
      type = 'RADII'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (x(l97+i-1),i=1,natom)
      type = 'SCREEN'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (x(l96+i-1),i=1,natom)
   end if
   
   !     ----- CALCULATE INVERSE, TOTAL MASSES -----
   
   !       -- save the masses for removal of mass weighted velocity,
   !          leaving the inverse masses in the legacy, Lwinv area
   
   
   tmass = 0.0d0
   !     -- index over molecules
   j = l75-1
   jj = i70-1
   !     -- index over mass->invmass
   k = lwinv-1
   !     -- index over saved mass
   l = lmass-1
   do n = 1,nspm
      j = j + 1
      jj = jj + 1
      x(j) = 0.0d0
      natsm = ix(jj)
      do nn = 1,natsm
         k = k+1
         l = l+1
         
         !         -- sum molecule
         
         x(j) = x(j) + x(k)
         
         !         -- save mass in "new" Lmass area
         
         x(l) = x(k)
         
         !         -- make inverse in "old" Lwinv area
         
         x(k) = 1.0d0 / x(k)
      end do
      tmass = tmass + x(j)
   end do
   tmassinv = 1.0d0 / tmass
   
   !     ----- SCALE THE CHARGES IF DIELC.GT.1.0E0 -----
   
   if (dielc /= 1.0e0) then
      dumd = sqrt(dielc)
      do i = 1,natom
         x(i+l15-1) = x(i+l15-1)/dumd
      end do
   end if
   
   !     ----- INVERT THE HBCUT ARRAY -----
   
   do i = 1,nphb
      if(hbcut(i) <= 0.001e0) hbcut(i) = 1.0d-10
      hbcut(i) = 1.0e0/hbcut(i)
   end do
   
   ! We no longer support pre-LEaP parmtop file anymore, therefore
   ! following 3 calls are disabled.
   !     ----- duplicate dihedral pointers for vector ephi -----
   !call dihdup(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48),pn)
   !call dihdup(nphia,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58),pn)
   
   !     --- pre-calculate some parameters for vector ephi ---
   !call dihpar(nptra,pk,pn,phase,gamc,gams,ipn,fmn)
   
#ifdef CHARMM
   write(6,'(a)') 'CHARMM is not supported by PBSA'
   call mexit(6,1)
   !    ---read in   1-4 parameters at end of parm file:
   !read(nf,9128) (cn114(i),   i = 1,nttyp)
   !read(nf,9128) (cn214(i),   i = 1,nttyp)
   
   !    --- read in   urey-bradley parameters:
   !read(nf,9128) (rkub(i),i=1,numang)
   !read(nf,9128) (rub(i),i=1,numang)
#endif
   
   return
   !9108 format(20a4)
   !9128 format(5e16.8)
   !9138 format(80A1)
end subroutine rdparm2 
