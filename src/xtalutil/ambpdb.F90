!   Simple filter program to convert AMBER coordinate files into PDB files
!   (and to do other similar conversions).

program ambpdb

!  implicit none

   character(len=80) :: arg,prmtop,title,ititl
   character(len=8) :: hbenec
   character(len=6) :: arg1,arg2
   character(len=4), dimension(:), allocatable :: igraph, lbres, isymbl
   character(len=1), dimension(:), allocatable :: ftype
   integer :: ntf, ioffset, iarg, indx, iargc, natom, nres
   integer :: maxatom, maxres, maxfhb ! --- required for -first option ---
   integer :: natomlen, nbondlen, nreslen
   integer :: nbond, ier, j, nhb
   logical :: alttit,center,aatm,bres,ter,bin,hasradii,remediate
   logical :: ext_pdb_data
   
   double precision, dimension(:), allocatable :: coords, fhbene, chg
   real, dimension(:), allocatable :: radius, atom_occupancy, atom_bfactor
   double precision :: hbene

   integer, dimension(:), allocatable :: ix, residue_number
   character(len=1), dimension(:), allocatable :: residue_chainid, &
      atom_altloc, residue_icode
   character(len=4), dimension(:), allocatable :: atom_element

   integer, dimension(:), allocatable :: ipres,lastat,ib,jb,fhybrid, &
                                         fhbdon,fhbh,fhbacc,fhbunit,itf,jtf

   ! ------ check argument options

   arg1 = '-pdb'
   prmtop = 'prmtop'
   alttit = .false.
   hbene = 1.0E10
   center = .false.
   aatm = .false.
   bres = .false.
   ter = .true.
   bin = .false.
   remediate = .true.
   ext_pdb_data = .false.
   ioffset = 0
   iarg = 0
   indx = iargc()
   
   if (indx .ne. 0) then
      do
         iarg = iarg + 1
         call getarg(iarg,arg)
         
         if (arg .eq. '-p') then
            iarg = iarg + 1
            call getarg(iarg,prmtop)
         else if (arg .eq. '-tit') then
            iarg = iarg + 1
            call getarg(iarg,title)
            alttit = .true.
         else if (arg .eq. '-ene') then
            iarg = iarg + 1
            call getarg(iarg,hbenec)
            read(hbenec,'(f8.0)' ) hbene                
         else if (arg .eq. '-sas') then
            arg1 = '-sas'
         else if (arg .eq. '-pqr') then
            arg1 = '-pqr'
         else if (arg .eq. '-bnd') then
            arg1 = '-bnd'
         else if (arg .eq. '-atm') then
            arg1 = '-atm'
         else if (arg .eq. '-first') then
            arg1 = '-first'
         else if (arg .eq. '-mol2') then
            arg1 = '-mol2'
         else if (arg .eq. '-ctr') then
            center = .true.
         else if (arg .eq. '-aatm') then
            aatm = .true.
         else if (arg .eq. '-bres') then
            bres = .true.
         else if (arg .eq. '-noter') then
            ter = .false.
         else if (arg .eq. '-noremediate') then
            remediate = .false.
         else if (arg .eq. '-ext') then
            ext_pdb_data = .true.
         else if (arg .eq. '-bin') then
            bin = .true.
         else if (arg .eq. '-offset') then
            iarg = iarg + 1
            call getarg(iarg,arg2)
            read(arg2,'(i5)' ) ioffset
         else if (arg .eq. '-h' .or. arg .eq. '-help' &
                  .or. arg .eq. '--help') then
            call usage()
         else
            write(6,'(/,5x,a,a)') 'unknown flag: ',arg
            stop
         end if

         if (iarg .ge. indx) exit
      end do
   end if
   
   ! ----- OPEN THE PARM FILE AND LOAD THE NECESSARY STUFF -----
   
   call amopen(10,prmtop,'O','F','R')
   call getnam0(natom,nres,nbond,10)
   rewind(unit=10)
   
   ! ------ Allocate memory: ------
   natomlen = natom
   nbondlen = nbond
   nreslen = nres
   if (arg1 .eq. '-first') then
      maxatom = 2*natom
      maxres = 2*natom
      maxfhb = natom
      natomlen = maxatom
      nbondlen = maxatom
      nreslen = maxres
   end if

   allocate( coords(3*natomlen), ix(3*natomlen), igraph(natomlen), &
            ipres(nreslen+1), lbres(nreslen), lastat(nreslen), &
            ib(nbondlen), jb(nbondlen), &
            chg(natomlen), ftype(natomlen), fhybrid(natomlen), &
            fhbdon(natomlen), fhbh(natomlen), fhbacc(natomlen), &
            fhbene(natomlen), fhbunit(natomlen), &
            itf(natomlen), jtf(natomlen), &
            radius(natomlen), &
            residue_number(nreslen), &
            residue_chainid(nreslen), atom_altloc(natomlen), &
            atom_element(natomlen), residue_icode(nreslen), &
            isymbl(natomlen), atom_occupancy(natomlen), &
            atom_bfactor(natomlen), stat = ier)
  
   if (ier /= 0) then
      write(0,*) 'memory allocation failure'
      call mexit(6,1)
   end if

   call getnam(natom,nres,igraph,ipres,lbres,ititl,10, &
               coords,ix,ib,jb,nbond,chg,radius,lastat,ter,hasradii, &
               residue_number, residue_chainid, atom_altloc, &
               atom_element, residue_icode, ext_pdb_data, isymbl, &
               atom_occupancy, atom_bfactor)
   close(unit=10)

   ! ----- READ THE COORDINATE FILE -----

   call getcor(bin,natom,coords,5)

   if (arg1 .eq. '-first') then      

      ! ----- Prepare for FIRST output -----
      call getftype(natom,nres,igraph,ipres,lbres,ib,jb,nbond,ftype)
      
      call gethybrid(nres,igraph,ipres,ib,jb,nbond,coords,fhybrid)
      
      call findtf(maxatom,nbond,ib,jb,fhybrid,igraph, &
                  ntf,itf,jtf)
      call corbondl(maxatom,natom,igraph,nres,ipres,lbres, &
                    nbond,ib,jb,coords)
      call findhbond(maxfhb,natom,igraph,ipres,lbres, &
                     ib,jb,nbond,coords,ftype,fhybrid, &
                     nhb,fhbdon,fhbh,fhbacc,fhbene,fhbunit)
      !call addsugarpsatoms(maxatom,maxres,natom,nres, &
      !                    igraph,ipres,lbres,ib,jb,nbond,coords, &
      !                    ftype,fhybrid)
      call addtether(maxfhb,maxatom,maxres,natom,nres, &
                     igraph,ipres,lbres,ib,jb,nbond,coords, &
                     ftype,fhybrid, &
                     nhb,fhbdon,fhbh,fhbacc,fhbene)
   end if

   ! ----- OUTPUT THE PDB, atm, pqr, first, bnd, or pluto  FILE -----

   if (     arg1 .eq. '-pdb' .or. arg1 .eq. '-atm' & 
       .or. arg1 .eq. '-pqr' .or. arg1 .eq. '-sas' &
       .or. arg1 .eq. '-first' ) then
      call genpdb(natom,nres,coords,igraph,ipres,lbres,ititl,6,arg1, &
                  alttit,title,center,chg,aatm,bres,ioffset,lastat, &
                  ftype,radius,hasradii,remediate, &
                  residue_number, residue_chainid, atom_altloc, &
                  atom_element, residue_icode, ext_pdb_data, &
                  atom_occupancy, atom_bfactor)
      if (arg1 .eq. '-first') then
          call writebond(6,nbond,ib,jb)
          call writetf(6,ntf,itf,jtf)
          call writehb(6,nhb,ftype,fhybrid,fhbdon,fhbh, &
                      fhbacc, fhbene, igraph, hbene, fhbunit)
      end if

   else if (arg1.eq.'-bnd') then
      do j=1,nbond
         write(6,'(2i5)') ib(j)/3+1,jb(j)/3+1
      end do

   else if (arg1.eq.'-plu') then
      call pluto(coords,natom,igraph,title)

   else if (arg1 .eq. '-mol2') then
      call genmol2(natom, nbond, nres, igraph, isymbl, coords, lbres, &
                   ipres, chg, ib, jb, ititl )
   else
      write(6,*) 'Needs -pdb, -pqr, -mol2, -sas, -atm, -bnd, -plu, or &
                 &-first flag'
      stop
   end if

   deallocate( coords, ix, igraph, ipres, lbres, lastat, ib, jb, chg, &
              ftype, fhybrid, fhbdon,fhbh, fhbacc, fhbene, itf, jtf, &
              stat = ier)
   if (ier /= 0 ) then
      write(0,*) 'memory deallocation failure'
      call mexit(6,1)
   end if
   call mexit(6,0)

end program ambpdb 

!=====================================================================

subroutine getnam0(natom,nres,mbona,nf)

   implicit none

   ! Arguments
   integer, intent(in) :: nf
   integer, intent(out) :: natom, nres, mbona

   ! Local variables
   integer :: i, iok, ntypes, nbonh, ntheth, nphih, &
              mtheta, mphia, nhparm, nparm, nnb
   character(len=80) :: fmt,fmtin,ifmt,afmt,rfmt,type,ititl
   
   ifmt = '(12i6)'
   afmt = '(20a4)'
   rfmt = '(5e16.8)'

   ! ----- Read the molecular topology -----
   fmtin = afmt
   type = 'TITLE'
   call nxtsec(nf,0,0,fmtin,type,fmt,iok)
   read(nf,'(a)') ititl

   fmtin = ifmt
   type = 'POINTERS'
   call nxtsec(nf,6,0,fmtin,type,fmt,iok)
   read(nf,fmt) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
                nhparm,nparm,nnb,nres

   ! ---- fix for allocating bond arrays ------
   mbona = mbona+nbonh

   return

end subroutine getnam0

!=====================================================================

subroutine getnam(natom,nres,igraph,ipres,lbres,ititl,nf, &
                  coords,ix,ib,jb,nbond,chg,radius,lastat,ter,hasradii, &
                  residue_number,residue_chainid,atom_altloc, &
                  atom_element,residue_icode,ext_pdb_data,isymbl, &
                  atom_occupancy, atom_bfactor)

   implicit none

   ! Arguments and returned values
   logical, intent(in) :: ter
   logical, intent(out) :: hasradii
   logical, intent(inout) :: ext_pdb_data
   integer, intent(in) :: nf
   integer, intent(out) :: nbond, ipres(*), residue_number(nres), &
                           lastat(*)
   integer, intent(inout) :: natom, nres, ix(*), ib(*), jb(*)
   double precision, intent(out) :: chg(*) 
   real, intent(out)             :: radius(*)
   double precision, intent(out) :: coords(*)
   character(len=4), intent(out) :: igraph(*), lbres(*), &
                                    atom_element(natom), isymbl(*)
   character(len=1), intent(out) :: residue_chainid(nres), atom_altloc(natom), &
                                    residue_icode(nres)
   real, intent(out) :: atom_occupancy(natom), atom_bfactor(natom)
   
   ! Local variables
   character(len=80) :: fmt,fmtin,ifmt,afmt,rfmt,type,ititl
   real :: xdummy
   integer :: i, iok, ifpert, ifbox, ifcap, imol, ires, ibond, &
              i1, i2, j1, j2, i13, i23, j13, j23, &
              mbona, mtheta, mphia, mbper, mgper, mdper, &
              ntypes, nbonh, ntheth, nphih, nparm, nnb, &
              nbona, ntheta, nphia, numbnd, numang, nptra, &
              natyp, nphb, nbper, ngper, ndper, nxmrs, &
              nhparm, ntype, nttyp, &
              idummy, jdummy, kdummy, ldummy, mdummy, &
              iptres, nspm, nspsol
   
   ifmt = '(12i6)'
   afmt = '(20a4)'
   rfmt = '(5e16.8)'
   hasradii = .false.

   ! ----- Read the molecular topology -----

   fmtin = afmt
   type = 'TITLE'
   call nxtsec(nf,0,0,fmtin,type,fmt,iok)
   read(nf,'(a)') ititl

   fmtin = ifmt
   type = 'POINTERS'
   call nxtsec(nf,6,0,fmtin,type,fmt,iok)
   read(nf,fmt) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
              nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
              numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper, &
              ndper,mbper,mgper,mdper,ifbox,nxmrs,ifcap

   ntype = ntypes*ntypes

   ! ----- READ THE SYMBOLS AND THE CHARGES AND THE MASSES -----

   fmtin = afmt
   type = 'ATOM_NAME'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (igraph(i),i = 1,natom)

   fmtin = rfmt
   type = 'CHARGE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (chg(i),i = 1,natom)
   do i=1,natom
      chg(i) = chg(i)/18.2223
   end do

   if (iok .eq. -1) then   ! this is an old-style prmtop file

      read(nf,40) (coords(i), i=1,natom)
      read(nf,30) (ix(i), i=1,natom)
      read(nf,30) (ix(i), i=1,natom)
      read(nf,30) (ix(i), i=1,ntype)
      read(nf,20) (lbres(i), i=1,nres)
      read(nf,30) (ipres(i), i=1,nres)
      ipres(nres+1) = natom+1
      20   format(20a4)
      30   format(12i6)
      40   format(5e16.8)

      ! ----- READ THE PARAMETERS -----

      read(nf,40) (coords(i), i=1,numbnd)
      read(nf,40) (coords(i), i=1,numbnd)
      read(nf,40) (coords(i), i=1,numang)
      read(nf,40) (coords(i), i=1,numang)
      read(nf,40) (coords(i), i=1,nptra)
      read(nf,40) (coords(i), i=1,nptra)
      read(nf,40) (coords(i), i=1,nptra)
      read(nf,40) (coords(i), i=1,natyp)

      nttyp = ntypes*(ntypes+1)/2

      read(nf,40) (coords(i),   i = 1,nttyp)
      read(nf,40) (coords(i),   i = 1,nttyp)

      ! ----- READ THE BONDING INFORMATION -----

      read(nf,30) (ib(i),jb(i),ix(i), i=1,nbonh)
      read(nf,30) (ib(i),jb(i),ix(i), i=nbonh+1,nbonh+nbona)
      nbond = nbonh + nbona

   else    !   this is new-style prmtop file:

      fmtin = afmt
      type = 'RESIDUE_LABEL'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (lbres(i),i=1,nres)

      fmtin = ifmt
      type = 'RESIDUE_POINTER'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ipres(i),i=1,nres)
      ipres(nres+1) = natom+1

      fmtin = ifmt
      type = 'BONDS_INC_HYDROGEN'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ib(i),jb(i),ix(i),i = 1,nbonh)

      fmtin = ifmt
      type = 'BONDS_WITHOUT_HYDROGEN'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ib(i),jb(i),ix(i),i = nbonh+1,nbonh+nbona)
      nbond = nbonh + nbona

      fmtin = rfmt
      type = 'RADII'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      if (iok .eq. -2) then
        write (0,*)
        write (0,*) 'WARNING: RADII FLAG NOT FOUND IN NEW PRMTOP'
        hasradii = .false.
      else  
        read(nf,fmt) (radius(i), i=1,natom)
        hasradii = .true.
      endif
      
      fmtin = afmt
      type = 'AMBER_ATOM_TYPE'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (isymbl(i), i=1,natom)

      !  defaults, if not using extended pdb data:
      residue_chainid = ' '
      residue_icode = ' '
      atom_altloc = ' '
      atom_occupancy = 1.0
      atom_bfactor = 0.0
      atom_element = '    '
      
      if (ext_pdb_data) then

         call nxtsec(nf,6,1,'*','RESIDUE_NUMBER',fmt,iok)
         if (iok == 0) then
            read(nf,fmt) residue_number
         else
            write(*,'(A)') &
                  'PRMTOP file did not have extended PDB data.'
            ext_pdb_data = .false.
            go to 10
         end if

         call nxtsec(nf,6,0,'*','RESIDUE_CHAINID',fmt,iok)
         if (iok == 0 ) read(nf,fmt) residue_chainid
         call nxtsec(nf,6,0,'*','ATOM_ELEMENT',fmt,iok)
         if (iok == 0 ) read(nf,fmt) atom_element

         call nxtsec(nf,6,1,'*','RESIDUE_ICODE',fmt,iok)
         if (iok == 0) read(nf,fmt) residue_icode
         call nxtsec(nf,6,1,'*','ATOM_ALTLOC',fmt,iok)
         if (iok == 0) read(nf,fmt) atom_altloc

         call nxtsec(nf,6,1,'*','ATOM_OCCUPANCY',fmt,iok)
         if (iok == 0) read(nf,fmt) atom_occupancy

         call nxtsec(nf,6,1,'*','ATOM_BFACTOR',fmt,iok)
         if (iok == 0) read(nf,fmt) atom_bfactor

      end if
   end if
  10  continue

#if 0
   if (ifbox .gt. 0) then

      ! --- prmtop file will have information needed to generate where the
      !   TER cards in the output should be:

      if (iok .eq. -1) then     !  old-style prmtop:
         read(nf,30) (idummy,jdummy,kdummy,ldummy,i=1,ntheth)
         read(nf,30) (idummy,jdummy,kdummy,ldummy,i=1,ntheta)
         read(nf,30) (idummy,jdummy,kdummy,ldummy,mdummy,i=1,nphih)
         read(nf,30) (idummy,jdummy,kdummy,ldummy,mdummy,i=1,nphia)
         read(nf,30) (idummy, i=1,nnb)
         read(nf,40) (xdummy,i=1,nphb)
         read(nf,40) (xdummy,i=1,nphb)
         read(nf,40) (xdummy,i=1,nphb)
         read(nf,20) (idummy,i=1,natom)
         read(nf,20) (idummy,i=1,natom)
         read(nf,30) (idummy,i=1,natom)
         read(nf,30) (idummy,i=1,natom)
      end if

      ! ----get molecule info to put TER cards in output:

      lastat(1) = natom
      fmtin = ifmt
      type = 'SOLVENT_POINTERS'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) iptres,nspm,nspsol

      fmtin = ifmt
      type = 'ATOMS_PER_MOLECULE'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (lastat(i), i=1,nspm)
      do i=2,nspm
         lastat(i) = lastat(i-1) + lastat(i)
      end do
      
   else if (ter) then
#endif

   if (ter) then
      
      ! --- check bonding file to find end of molecules: assume that if there
      !     are no bonds between two adjacent residues, there should be a TER
      !     card:
      
      imol = 1
      do ires=2,nres
         i1 = ipres(ires-1)
         i2 = ipres(ires) - 1
         j1 = ipres(ires)
         j2 = ipres(ires+1) - 1
         i13 = 3*(i1-1)
         i23 = 3*(i2-1)
         j13 = 3*(j1-1)
         j23 = 3*(j2-1)
         do ibond=1,nbond
            if(    ib(ibond).ge.i13 .and. ib(ibond).le.i23 &
                      .and. jb(ibond).ge.j13 .and. jb(ibond).le.j23 ) goto 50
            if(    jb(ibond).ge.i13 .and. jb(ibond).le.i23 &
                      .and. ib(ibond).ge.j13 .and. ib(ibond).le.j23 ) goto 50
         end do
         ! --- here if there is no bond between the two residues:
         lastat(imol) = i2
         imol = imol + 1
         50 continue
      end do
      lastat(imol) = natom
      
   else
      
      ! ----  user has requested no TER cards:
      lastat(1) = 0
      
   end if
   
   return
   
end subroutine getnam

!=====================================================================

subroutine genpdb(natom,nres,coords,igraph,ipres,lbres,ititl,nf,arg1, &
                  alttit,title,center,chg,aatm,bres,ioffset,lastat, &
                  ftype,radius,hasradii,remediate, &
                  residue_number, residue_chainid, atom_altloc, &
                  atom_element, residue_icode, ext_pdb_data, &
                  atom_occupancy, atom_bfactor)
   
   implicit none

   ! Arguments and returning values
   logical, intent(in) :: alttit, center, aatm, bres, hasradii, &
                          remediate, ext_pdb_data
   integer, intent(in) :: natom, nres, nf, ipres(*), ioffset, &
                          lastat(*), residue_number(nres)
   real, intent(in) :: radius(*)
   double precision, intent(in) :: chg(*)
   double precision, intent(inout) :: coords(*)
   character(len=1), intent(in) :: ftype(*)
   character(len=4), intent(in) :: igraph(*), atom_element(natom)
   character(len=1), intent(in) :: residue_chainid(nres), atom_altloc(natom), &
                                   residue_icode(nres)
   character(len=4), intent(inout) :: lbres(*)
   character(len=6), intent(in) :: arg1
   character(len=40),intent(in) :: title
   character(len=80),intent(in) :: ititl
   real, intent(in) :: atom_occupancy(natom), atom_bfactor(natom)

   ! Local variables   
   integer :: i, iat, imol, isrn, itype, &
              j, j1, j2, kc, jpat,j_print, &
              k, k_print, m, nn
   real :: elrad(15)   
   double precision :: xs, ys, zs
   character(len=1) :: type
   character(len=2) :: element(0:94), resnam
   character(len=4) :: code, atnam, tmpnam
   character(len=23) :: occb
   

   ! from Sitkoff, Sharp & Honig, PARSE values, J. Phys. Chem. 98:1978(1994)
   !elrad=(/2.0, 1.4, 1.0, 1.5, 1.0, 0.0, 1.7, 1.85, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/)

   ! Bondi radii, from Huheey, Inorganic Chemistry: Principles of Structure
   ! and Reactivity, 2nd ed. p. 232.
   elrad=(/1.85, 1.5, 1.0, 1.55, 1.0, 0.0, 1.7, 1.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.2/)

   ! dummy occupation and b-factor string:
   occb='  1.00  0.00           '

   ! element names:
      element(0) = '  '
      element(1) = ' H'
      element(2) = 'He'
      element(3) = 'Li'
      element(4) = 'Be'
      element(5) = ' B'
      element(6) = ' C'
      element(7) = ' N'
      element(8) = ' O'
      element(9) = ' F'
      element(10) = 'Ne'
      element(11) = 'Na'
      element(12) = 'Mg'
      element(13) = 'Al'
      element(14) = 'Si'
      element(15) = ' P'
      element(16) = ' S'
      element(17) = 'Cl'
      element(18) = 'Ar'
      element(19) = ' K'
      element(20) = 'Ca'
      element(21) = 'Sc'
      element(22) = 'Ti'
      element(23) = ' V'
      element(24) = 'Cr'
      element(25) = 'Mn'
      element(26) = 'Fe'
      element(27) = 'Co'
      element(28) = 'Ni'
      element(29) = 'Cu'
      element(30) = 'Zn'
      element(31) = 'Ga'
      element(32) = 'Ge'
      element(33) = 'As'
      element(34) = 'Se'
      element(35) = 'Br'
      element(36) = 'Kr'
      element(37) = 'Rb'
      element(38) = 'Sr'
      element(39) = ' Y'
      element(40) = 'Zr'
      element(41) = 'Nb'
      element(42) = 'Mo'
      element(43) = 'Tc'
      element(44) = 'Ru'
      element(45) = 'Rh'
      element(46) = 'Pd'
      element(47) = 'Ag'
      element(48) = 'Cd'
      element(49) = 'In'
      element(50) = 'Sn'
      element(51) = 'Sb'
      element(52) = 'Te'
      element(53) = ' I'
      element(54) = 'Xe'
      element(55) = 'Cs'
      element(56) = 'Ba'
      element(57) = 'La'
      element(58) = 'Ce'
      element(59) = 'Pr'
      element(60) = 'Nd'
      element(61) = 'Pm'
      element(62) = 'Sm'
      element(63) = 'Eu'
      element(64) = 'Gd'
      element(65) = 'Tb'
      element(66) = 'Dy'
      element(67) = 'Ho'
      element(68) = 'Er'
      element(69) = 'Tm'
      element(70) = 'Yb'
      element(71) = 'Lu'
      element(72) = 'Hf'
      element(73) = 'Ta'
      element(74) = ' W'
      element(75) = 'Re'
      element(76) = 'Os'
      element(77) = 'Ir'
      element(78) = 'Pt'
      element(79) = 'Au'
      element(80) = 'Hg'
      element(81) = 'Tl'
      element(82) = 'Pb'
      element(83) = 'Bi'
      element(84) = 'Po'
      element(85) = 'At'
      element(86) = 'Rn'
      element(87) = 'Fr'
      element(88) = 'Ra'
      element(89) = 'Ac'
      element(90) = 'Th'
      element(91) = 'Pa'
      element(92) = ' U'
      element(93) = 'Np'
      element(94) = 'Pu'

   ! ---isrn =2 makes this a "normal" atom for ms
   isrn = 2
   nn = 1
   if (alttit) then
      write(nf,100) title
   else
      write(nf,100) ititl(1:40)
   end if

   ! -- if -ctr was requested, need to center molecule at origin

   if (center) then
      xs = 0.0
      ys = 0.0
      zs = 0.0
      k = 0
      do iat=1,natom
         xs = xs + coords(k+1)
         ys = ys + coords(k+2)
         zs = zs + coords(k+3)
         k = k + 3
      end do
      xs = xs /natom
      ys = ys /natom
      zs = zs /natom
      k = 0
      do iat=1,natom
         coords(k+1) = coords(k+1) - xs
         coords(k+2) = coords(k+2) - ys
         coords(k+3) = coords(k+3) - zs
         k = k + 3
      end do
   end if
   
   imol = 1
   do j=1,nres
      j1 = ipres(j)
      j2 = ipres(j+1)-1
      
      if (bres) then
         ! ---convert protein residue names back to more like Brookhaven format:

         if (lbres(j).eq.'HID ' .or. lbres(j).eq.'HIE ' .or. &
            lbres(j).eq.'HIP '.or. lbres(j).eq.'HIC ') lbres(j) = 'HIS '
         if (lbres(j).eq.'CYX ') lbres(j) = 'CYS '
         if (lbres(j).eq.'CYM ') lbres(j) = 'CYS '
         if (lbres(j).eq.'MEM ') lbres(j) = 'MET '
         if (lbres(j).eq.'ASH ') lbres(j) = 'ASP '
         if (lbres(j).eq.'GLH ') lbres(j) = 'GLU '

         ! ---also for nucleic acid names:

         if( lbres(j).eq.'G   ' ) lbres(j) = '  G '
         if( lbres(j).eq.'DG  ' ) lbres(j) = ' DG '
         if( lbres(j).eq.'C   ' ) lbres(j) = '  C '
         if( lbres(j).eq.'DC  ' ) lbres(j) = ' DC '
         if( lbres(j).eq.'A   ' ) lbres(j) = '  A '
         if( lbres(j).eq.'DA  ' ) lbres(j) = ' DA '
         if( lbres(j).eq.'U   ' ) lbres(j) = '  U '
         if( lbres(j).eq.'DT  ' ) lbres(j) = ' DT '
         if( lbres(j).eq.'G5  ' ) lbres(j) = '  G '
         if( lbres(j).eq.'C5  ' ) lbres(j) = '  C '
         if( lbres(j).eq.'A5  ' ) lbres(j) = '  A '
         if( lbres(j).eq.'U5  ' ) lbres(j) = '  U '
         if( lbres(j).eq.'G3  ' ) lbres(j) = '  G '
         if( lbres(j).eq.'C3  ' ) lbres(j) = '  C '
         if( lbres(j).eq.'A3  ' ) lbres(j) = '  A '
         if( lbres(j).eq.'U3  ' ) lbres(j) = '  U '
      end if
      
      do k=j1,j2
         k_print = mod(k,100000)
         if (aatm) then
            write(atnam,'(a4)') igraph(k)
         else
            
            ! ---convert atom names to closely resemble those used by the
            ! wwPDB in its version 3 (aka "remdiated") files:
            
            ! ---First, assume that there are no two-character element names
            ! (like Fe or Ca or Na).  Then, according to Brookhaven rules,
            ! column 13 will be blank, and the name will be left-justified
            ! starting in column 14.  UNLESS, the name is four characters
            ! long!  In that case, don't use the first blank.
            
            resnam = lbres(j)(1:3)
            write(tmpnam,'(a4)') igraph(k)
            if (tmpnam(4:4) .eq. ' ') then 
               atnam(1:1) = ' '
               atnam(2:4) = tmpnam(1:3)
            else
               atnam(1:4) = tmpnam(1:4)
            endif
            
            ! --- Now, special case out the two-character element names:
            if (atnam(1:4) .eq. ' Na+' .or. atnam(1:4) .eq. ' NA+' .or. &
                atnam(1:3) .eq. ' Fe' .or. atnam(1:3) .eq. ' FE' .or. &
                atnam(1:3) .eq. ' Cl' .or. atnam(1:3) .eq. ' CL' .or. &
                atnam(1:3) .eq. ' Zn' .or. atnam(1:3) .eq. ' ZN' .or. &
                atnam(1:4) .eq. ' Li+' .or. atnam(1:4) .eq. ' LI+' .or. &
                atnam(1:4) .eq. ' Ca+' .or. atnam(1:4) .eq. ' CA+' .or. &
                atnam(1:4) .eq. ' Mg+' .or. atnam(1:4) .eq. ' MG+' .or. &
                atnam(1:4) .eq. ' Br-' .or. atnam(1:4) .eq. ' BR-' ) then
               atnam(1:1) = atnam(2:2)
               atnam(2:2) = atnam(3:3)
               atnam(3:3) = atnam(4:4)
               atnam(4:4) = ' '
            end if
         end if

         kc = 3*k-3
         if( ext_pdb_data ) then
            j_print = residue_number(j)
         else
            j_print = mod(j+ioffset,10000)
         endif
         if (arg1 .eq. '-pdb') then
            if (ext_pdb_data) then
               write(nf, &
                  '(A4,I7,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A4)')  &
                  'ATOM',k_print,atnam,atom_altloc(j),lbres(j), &
                  residue_chainid(j),residue_number(j), &
                  residue_icode(j)(1:1),coords(kc+1:kc+3), &
                  atom_occupancy(k), atom_bfactor(k), atom_element(k)
            else
               write(nf, &
                  '(A4,I7,1X,A4,1X,A3,2X,I4,4X,3F8.3,2F6.2,10X)') 'ATOM', &
                  k_print,atnam,lbres(j),j_print,(coords(kc+m),m=1,3), &
                  1.0,0.0
            end if
         else if (arg1 .eq. '-atm' .or. arg1 .eq. '-pqr'  &
                                   .or. arg1 .eq. '-sas') then
            code = igraph(k)
            type = code(1:1)
            if (type .eq. 'H') then
               itype = 15
            else if (type .eq. 'C') then
               itype = 7
            else if (type .eq. 'N') then
               itype = 4
            else if (type .eq. 'O') then
               itype = 2
            else if (type .eq. 'P') then
               itype = 1
            else if (type .eq. 'S') then
               itype = 8
            else if (type .eq. 'L') then
               itype = 16
            else if (type .eq. 'Z') then
               itype = 5
            else
               itype = 3
            end if
            if (arg1 .eq. '-atm') then
               write(nf,50) (coords(kc+m),m=1,3), &
                            itype,isrn,nn,lbres(j),j+ioffset,igraph(k)
            else if (arg1 .eq. '-pqr') then
               if (hasradii) then
                  if (ext_pdb_data) then
                     write(nf, &
                  '(A4,I7,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F8.4,6X,A2)')  &
                        'ATOM',k_print,atnam,atom_altloc(j),lbres(j), &
                        residue_chainid(j),residue_number(j), &
                        residue_icode(j),coords(kc+1:kc+3),chg(k), &
                        radius(k),atom_element(k)
                  else
                     write(nf,81) k_print,atnam,lbres(j),j_print, &
                       (coords(kc+m),m=1,3),chg(k),radius(k)
                  endif
               else
                  if (ext_pdb_data) then
                     write(nf, &
                  '(A4,I7,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F8.4,6X,A2)')  &
                        'ATOM',k_print,atnam,atom_altloc(j),lbres(j), &
                        residue_chainid(j),residue_number(j), &
                        residue_icode(j),coords(kc+1:kc+3),chg(k), &
                        radius(k),atom_element(k)
                  else
                     write(nf,81) k_print,atnam,lbres(j),j_print, &
                       (coords(kc+m),m=1,3),chg(k),elrad(itype)
                  endif
               end if
            else if (arg1.eq.'-sas') then
               write(nf,81) k_print,atnam,lbres(j),j_print, &
                    (coords(kc+m),m=1,3),chg(k),elrad(itype)+1.4
            end if
         else if (arg1 .eq. '-first') then
            if (atnam .eq. ' X' .and. lbres(j) .eq. 'BMH') then
               write(nf,67) k_print,atnam,lbres(j),j_print - jpat, &
                            (coords(kc+m),m=1,3), 0.00, 0.00, &
                            j_print - jpat, ftype(k)
            else
               write(nf,65) k_print,atnam,lbres(j),j_print, &
                            (coords(kc+m),m=1,3),0.00, 0.00, j_print, ftype(k)
               jpat = j_print
            end if
         end if
         if (k .eq. lastat(imol)) then
            write(nf,120)
            imol = imol + 1
         end if
      end do
   end do

   write(nf,130)

   return

    50 format(3f10.3,3i5,5x,a4,i5,1x,a4)
    60 format('ATOM',2X,I5,1X,A4,1X,A4,1X, I4,4X,3F8.3,A23,A3)
    61 format('ATOM',2X,I5,1X,A4,1X,A4,1X, I4,4X,3F8.3)
    65 format('ATOM',2X,I5,1X,A4,1X,A4,1X, I4,4X,3F8.3,2F6.2,I5,1X,A1)
    67 format('HETATM', I5,1X,A4,1X,A4,'G',I4,4X,3F8.3,2F6.2,I5,1X,A1)
    81 format('ATOM',2X,I5,1X,A4,1X,A4,1X, I4,4X,3F8.3,f8.4,f8.3)
    90 format('REMARK  ',15a4)
   100 format('REMARK  ',a40)
   120 format('TER   ')
   130 format('END   ')
  
end subroutine genpdb

!=====================================================================

subroutine getcor(bin,natom,coords,nf)

   implicit none

   ! Arguments and return values
   logical, intent(in) :: bin
   integer, intent(in) :: natom, nf
   double precision, intent(inout) :: coords(*)

   ! Local variables   
   integer :: i, nat3, matom
   character(len=80) :: line,ititl

   nat3 = 3*natom

   if (.not. bin) then
      read(nf,'(a)') ititl
      read(nf,'(a)') line
      if (line(6:6) .eq. ' ') then
         read(line,'(i5)') matom
      else
         read(line,'(i6)') matom
      end if
      if (matom .ne. natom) then
         write(0,20) natom,matom
         stop
      end if
      read(nf,50) (coords(i),i=1,nat3)
      return
   else
      read(nf) ititl
      read(nf) matom
      if (matom .ne. natom) then
         write(0,20) natom,matom
         stop
      end if
      read(nf,end=15,err=15) (coords(i), i=1,nat3)
      return
   end if

   15 continue
   write(0,'(a,a)') 'COULD NOT READ COORDINATES FROM RESTRT FILE'
   stop

   20 format(/2x,'ATOMS DO NOT MATCH BETWEEN PARM AND MIN FILES',2i8)
   50 format(6f12.7)

end subroutine getcor

!=====================================================================

subroutine pluto(coords,natom,igraph,title)

   implicit none

   integer, intent(in) :: natom
   double precision, intent(in) :: coords(*)
   character(len=4), intent(in) :: igraph(*)
   character(len=40), intent(in) :: title

   integer :: j, k
   character(len=1) :: blnk, c1
   character(len=4) :: atnam

   ! ------ program to write pluto input file -----

   write(6,120)
   blnk = ' '
   k = 0
   do j=1,natom
      atnam = igraph(j)
      c1 = atnam(1:1)
      if (j.le.9) then
         write(6,130) c1,j,coords(k+1),coords(k+2),coords(k+3)
      else if (j.le.99) then
         write(6,140) c1,j,coords(k+1),coords(k+2),coords(k+3)
      else
         write(6,150) c1,j,coords(k+1),coords(k+2),coords(k+3)
      end if
      k = k + 3
   end do
   write(6,30)
   write(6,20) title
   write(6,50)
   write(6,90)
   write(6,100)
   write(6,110)
   write(6,60)
   write(6,70)
   write(6,80)
   return

    20 format('TITLE',1x,a40)
    30 format('*')
    50 format('JOIN RADII C 0.8 N 0.8 H 0.4 O 0.74 P 1.1 S 1.1')
    60 format('STEREO')
    70 format('VIEW YORIGIN')
    80 format('PLOT')
    90 format('SOLID')
   100 format('RADII ATOMS C 0.3 N 0.4 O 0.4 H 0.1')
   110 format('RADII BONDS 0.05 8 TAPER 5')
   120 format('DATA',1x,a40)
   130 format(a1,i1,t9,3f8.3)
   140 format(a1,i2,t9,3f8.3)
   150 format(a1,i3,t9,3f8.3)

end subroutine pluto

!=====================================================================

subroutine usage()
   
   implicit none
   
   write(*,'(A)') 'Usage:', &
                  'ambpdb [OPTION]... < restrt > out.pdb', &
                  '', &
                  'Options:', &
                  ' -p PRMTOP     Define PRMTOP filename (default:"prmtop").', &
                  ' -tit TITLE    Write a REMARK record containing TITLE.', &
                  '                   (default: use prmtop title)', &
                  ' -aatm         Left-justified Amber atom names.', &
                  ' -bres         Brookhaven Residue names (HIE->HIS, etc.).', &
                  ' -ctr          Center molecule on (0,0,0).', &
                  ' -noter        Do not write TER records.', &
                  ' -ext          Use PRMTOP extended PDB info, if present.', &
                  ' -ene FLOAT    Define H-bond energy cutoff for FIRST.', &
                  ' -bin          The coordinate file is in binary form.', &
                  ' -offset INT   Add offset to residue numbers.', &
                  '', &
                  'Options for alternate output format (give only one of these):', &
                  ' -pqr          PQR (MEAD) format with charges and radii.', &
                  ' -sas          PQR with 1.4 added to atom radii.', &
                  ' -mol2         TRIPOS MOL2 format.', &
                  ' -bnd          list bonds from the PRMTOP.', &
                  ' -atm          Mike Connolly surface/volume format.', &
                  ' -first        Add REMARKs for input to FIRST.'

   call mexit(6,0)

end subroutine usage

subroutine genmol2(natom, nbond, nres, igraph, isymbl, coords, lbres, &
                   ipres, charge, ib, jb, ititl )
   ! This subroutine prints a mol2 file to stdout
   implicit none

   ! Passed arguments
   integer, intent(in)           :: natom, nbond, nres, ib(*), jb(*), ipres(*)
   character (len=4), intent(in) :: igraph(*), isymbl(*), lbres(*)
   character (len=80), intent(in) :: ititl
   double precision, intent(in)  :: charge(*), coords(*)

   ! Local variables
   integer     :: i, residue_container(natom)
! Variable descriptions
   !
   ! Passed Variables
   !  natom             : number of atoms in system
   !  nbond             : number of bonds in system
   !  nres              : number of residues in system
   !  igraph            : array of atom names
   !  isymbl            : array of atom types
   !  coords            : cartesian coordinates of each atom
   !  lbres             : residue labels
   !  ipres             : residue pointers
   !  charge            : partial charge array
   !  ib                : 1st atom in each bond array
   !  jb                : 2nd atom in each bond array
   !  ititl             : title of system from prmtop
   !
   !
   !  i                 : loop counter
   !  residue_container : lists which residue each atom belongs to

   ! First call the subroutine to assign each atom to a residue
   call assign_res(natom, nres, ipres, residue_container)

   ! MOLECULE section
   write(6,10) 'MOLECULE'
   write(6,'(a)') ititl
   write(6,30) natom, nbond, nres, 0, 1
   write(6,40) 'SMALL'
   write(6,40) 'USER_CHARGES'

   ! ATOM section
   write(6,10) 'ATOM'
   do i = 1, natom
      write(6,50) i, igraph(i), coords(3*i-2), coords(3*i-1), coords(3*i), &
                  isymbl(i), residue_container(i), lbres(residue_container(i)), &
                  charge(i), '****'
   end do

   ! BOND section
   write(6,10) 'BOND'
   do i = 1, nbond
      write(6,60) i, ib(i)/3+1, jb(i)/3+1, 1
   end do

   ! SUBSTRUCTURE section
   write(6,10) 'SUBSTRUCTURE'
   do i = 1, nres
      write(6,70) i, lbres(i), ipres(i), '****', 0, '****', '****'
   end do

   return

   10 format('@<TRIPOS>',a)
   20 format(20a4)
   30 format(i6,1x,i7,1x,i6,1x,i3,1x,i3)
   40 format(a)
   50 format(i6,2x,1a4,1x,3f10.3,2x,a4,1x,i6,1x,a4,1x,f9.4,1x,a4)
   60 format(i6,1x,2i7,2x,i1)
   70 format(i6,1x,a4,1x,i7,1x,a4,5x,i1,1x,a4,1x,a4)

end subroutine genmol2

subroutine assign_res(natom, nres, ipres, residue_container)
   ! This subroutine assigns each atom to a specific residue
   implicit none

   ! Passed arguments
   integer, intent(in)  :: natom, nres, ipres(*)
   integer, intent(out) :: residue_container(*)

   ! Local variables
   integer  :: i, j, curres, nextres

   ! Variable descriptions
   !
   ! Passed Variables
   !  natom             : number of atoms in system
   !  nres              : number of residues in system
   !  ipres             : array of "first-atom" residue pointers
   !  residue_container : lists which residue each atom belongs to
   !
   ! Local Variables
   !  i, j              : loop counters
   !  curres            : first atom of the current residue
   !  nextres           : first atom of the next residue

   ! Loop through all residue pointers except the last one
   do i = 1, nres - 1
      curres  = ipres(i  )
      nextres = ipres(i+1)

      do j = curres, nextres - 1
         residue_container(j) = i
      end do
   end do

   ! Now set the last residue
   do j = ipres(nres), natom
      residue_container(j) = nres
   end do

end subroutine assign_res
