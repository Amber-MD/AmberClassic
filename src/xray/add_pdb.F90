program add_pdb
!-------------------------------------------------------------------------------
! This program adds data entries to the PRMTOP file from a matching PDB file.
! It is very simple, and requires all residues to be in the same order.
!
! The following property arrays are added:
!
!  %FLAG ATOM_OCCUPANCY
!  %COMMENT Atom occupancies read from the PDB file
!  %FORMAT(10F8.2)

!  %FLAG ATOM_BFACTOR
!  %COMMENT B-factor read from PDB file; DIMENSION(NATOM)
!  %FORMAT(10F8.2)
!
!  %FLAG ATOM_NUMBER   ! not written in this version
!  %COMMENT Atom serial number from the PDB file
!  %FORMAT(10I8)

!  %FLAG RESIDUE_NUMBER
!  %COMMENT Residue number (resSeq) read from PDB file; DIMENSION(NRES)
!  %FORMAT(20I4)
!
!  %FLAG RESIDUE_CHAINID
!  %COMMENT Residue chain ID (chainId) read from PDB file; DIMENSION(NRES)
!  %FORMAT(80A1)
!
!  %FLAG RESIDUE_ICODE
!  %COMMENT Residue insertion code (iCode) read from PDB file; DIMENSION(NRES)
!  If the PDB file does not use iCodes (most do not), then the following comment
!  is added to help keep the new data mostly self-documenting:
!  %COMMENT Residue insertion code (iCode) not present in PDB file
!  %COMMENT If all iCodes are blank, this array is not stored in the PRMTOP file.
!  %FORMAT(80A1)
!
!  %FLAG ATOM_ELEMENT
!  %COMMENT Atom element and formal charge, as read from PDB file; DIMENSION(NATOM)
!  %FORMAT(20A4)
!
!  The ELEMENT field is the 2-letter element name, right justified, followed
!  by 2 characters to define the formal charge for ions. For example:
!  ' H  ' is hydrogen, 'HE  ' is helium, and 'FE+2' is an iron ion.
!
!  %FLAG ATOM_ALTLOC
!  %COMMENT Alternate location code; DIMENSION(NATOM)
!  %FORMAT(80A1)
!
!-------------------------------------------------------------------------------
   implicit none
   integer, parameter :: real_kind=8, rk=real_kind
   integer, parameter :: path_max=256
   integer, parameter :: stdin=5, stdout=6, stderr=0
   
   integer, parameter :: pdb_lun=10, in_lun=11, out_lun=12, crd_lun=13
   
   character(len=path_max) :: &
         prmtop_infile='', &
         prmtop_outfile='', &
         pdb_infile='', &
         crd_infile=''

   integer, parameter :: num_elements = 111
   character(len=2), parameter :: element_name(num_elements) = (/ &
       ' H', 'HE', 'LI', 'BE', ' B', ' C', ' N', ' O', ' F', 'NE', &
       'NA', 'MG', 'AL', 'SI', ' P', ' S', 'CL', 'AR', ' K', 'CA', &
       'SC', 'TI', ' V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', &
       'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', ' Y', 'ZR', &
       'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', &
       'SB', 'TE', ' I', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', &
       'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', &
       'LU', 'HF', 'TA', ' W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', &
       'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', &
       'PA', ' U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', &
       'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', &
       'RG' /)

   integer :: i, ierr
   character(len=32) :: fmt
   integer :: pointers(32)
   integer :: natom, nres
   equivalence(natom,pointers(1))
   equivalence(nres,pointers(12))
   
   character(len=4) :: name,segID
   character(len=3) :: resName,prev_resName
   character(len=4) :: element
   character(len=1) :: altLoc,chainID,iCode,prev_chainID,prev_iCode,prev_altLoc
   character(len=5) :: serial
   integer :: resSeq,prev_resSeq
   real :: xyz(3),occupancy,tempFactor
   real, allocatable :: coor(:,:)
   real :: d
   logical :: warned_missing_elements = .false.
   
   ! Existing PRMTOP data:
   character(len=4), allocatable :: atom_name(:)
   character(len=4), allocatable :: residue_label(:)
   integer, allocatable :: residue_pointer(:)
   ! New PRMTOP data:
   character(len=1), allocatable :: residue_chainid(:), residue_icode(:), &
                                    atom_altloc(:)
   character(len=4), allocatable :: atom_element(:)
   character(len=5), allocatable :: atom_number(:)
   integer, allocatable :: residue_number(:)
   real, allocatable :: atom_bfactor(:), atom_occupancy(:)
   integer :: pdb_nres, pdb_natom
   logical :: guess_all
   
   ! Current PRMTOP files are limited to 80-character lines.
   ! This allows a bit more for future versions, or long comments.
   character(len=128) :: buf
   
   character(len=32) :: arg
   integer :: argc, iarg
   ! begin
   guess_all = .false.
   argc = command_argument_count()
   if (argc==0) call usage()
   iarg = 0
   do while(iarg<argc)
      iarg=iarg+1
      call getarg(iarg,arg)
      select case(trim(arg))
      case('-h')
         call usage()
      case('-guess')
         guess_all = .true.
      case('-i')
         call get_next_arg(prmtop_infile)
      case('-o')
         call get_next_arg(prmtop_outfile)
      case('-p')
         call get_next_arg(pdb_infile)
      case('-c')
         call get_next_arg(crd_infile)
      case default
         write(*,*) 'Unrecognize option: ',trim(arg)
         call usage()
      end select
   end do
   
   if (prmtop_infile == '' .or. prmtop_outfile == '' &
         .or. pdb_infile == '') then
      write(*,*) &
            'Error: you must define prmtop in and out files, and a PDB file.'
      call usage()
   end if
   
   open(unit=pdb_lun,file=pdb_infile,status='OLD', &
         action='READ',form='FORMATTED')
   open(unit=in_lun,file=prmtop_infile,status='OLD', &
         action='READ',form='FORMATTED')
   open(unit=out_lun,file=prmtop_outfile,status='UNKNOWN', &
         action='WRITE',form='FORMATTED')
   
   call nxtsec(in_lun,STDOUT,0,'*','POINTERS',fmt,ierr)
   read(in_lun,fmt) pointers
   
   if (crd_infile /= '') then
      open(unit=crd_lun,file=crd_infile,status='OLD', &
            action='READ',form='FORMATTED')
      read(crd_lun,'(A)') buf
      read(crd_lun,*) i
      if (i /= natom) then
         write(*,*) 'INPCRD file has mismatched atom count'
         stop 666
      end if
      allocate(coor(3,natom))
      read(crd_lun,'(6F12.7)') coor
   end if
   
   allocate(atom_name(natom),residue_label(nres),residue_pointer(nres+1), &
         residue_chainid(nres), residue_icode(nres), residue_number(nres), &
         atom_element(natom),atom_bfactor(natom),atom_occupancy(natom), &
         atom_number(natom),atom_altloc(natom) )
   
   call nxtsec(in_lun,STDOUT,0,'*','ATOM_NAME',fmt,ierr)
   read(in_lun,fmt) atom_name
   
   call nxtsec(in_lun,STDOUT,0,'*','RESIDUE_LABEL',fmt,ierr)
   read(in_lun,fmt) residue_label
   
   call nxtsec(in_lun,STDOUT,0,'*','RESIDUE_POINTER',fmt,ierr)
   read(in_lun,fmt) residue_pointer(1:nres)
   residue_pointer(nres+1)=natom+1
   
   residue_chainid(:) = ' '
   residue_icode(:) = ' '
   residue_number(:) = 0
   atom_element(:) = '????'
   atom_bfactor(:) = 15.0
   atom_occupancy(:) = 1.0
   atom_number(:) = '    1'
   atom_altloc(:) = ' '
   
   prev_iCode='*'
   prev_resSeq=HUGE(prev_resSeq)
   prev_chainID='*'
   prev_resName='***'
   pdb_nres=0
   pdb_natom=0
   
   do
      read(pdb_lun,'(A)',end=99) buf
      if (buf(1:6)=='END   ') exit
      if (buf(1:6)=='ATOM  ' .or. buf(1:6)=='HETATM') then
         !  write(0,*) buf(1:70)
         read(buf,'(6X,A5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,2A4)') &
               serial,name,altLoc,resName,chainID,resSeq,iCode, &
               xyz,occupancy,tempFactor,segID,element
         pdb_natom=pdb_natom+1
         ! If any of these properties change, a new residue has begun.
         if (prev_resSeq /= resSeq .or. prev_chainID /= chainID &
               .or. prev_iCode /= iCode) then
            prev_resSeq = resSeq
            prev_iCode = iCode
            prev_chainID = chainID
            prev_resName = resName
            pdb_nres=pdb_nres+1
            if (pdb_nres > nres) then
               write(*,*)'PDB file has too many residues!'
               exit
            end if
            if (residue_label(pdb_nres)(1:3) /= resName) then
               if (residue_label(pdb_nres) /= 'WAT') then
                  write(*,*) 'Warning: PDB resName "',resName, &
                       '" /= PRMTOP label "',residue_label(pdb_nres)(1:3),'"'
               end if
            end if
            residue_chainid(pdb_nres) = chainID
            residue_icode(pdb_nres) = iCode
            residue_number(pdb_nres) = resSeq
         end if
         i = find_atom(pdb_nres,adjustl(name),altLoc)
         if (i<0) then
            write(*,*) 'Atom not found: ',trim(buf)
            stop
         end if
         if (crd_infile /= '') then
            d = sqrt(sum(coor(:,i)-xyz)**2)
            if (d > 0.3) then
               write(*,*)'Atom is ',d,' from PDB coordinate'
               write(*,*) coor(:,i)
               write(*,*) trim(buf)
               stop
            end if
         end if
         
         if (element(1:1) == char(0)) stop 'BAD ELEMENT'
         ! This element determination requires that names are element-aligned
         ! unless file contains explicit element+charge data.
         if (element == '    ') then
            if (guess_all) then
               element = guess_element(atom_name(i))
            else
               if (.not.warned_missing_elements) then
                  warned_missing_elements = .true.

                  write(*,'(A)') &
           'Warning: This PDB file does not contain atomic element data', &
           'for all atoms. Element determination will assume that atom', &
           'labels are element-aligned: 1-letter element codes must', &
           'contain a leading space or a wrapped digit.'

               end if
               if (name(1:1) >= '0' .and. name(1:1) <= '9') then
                  element = ' '//name(2:2)//'  '
               else
                  element = name(1:2)//'  '
               end if
            end if
         end if
#if 0
         if (atom_element(i) /= '????') then
            write(*,*) 'Duplicate atom: ',trim(buf)
            stop
         end if
#endif
         atom_element(i) = element
         atom_altloc(i) = altLoc
         atom_bfactor(i) = tempFactor
         atom_number(i) = serial
         atom_occupancy(i) = occupancy
      end if
   end do

   ! For all atoms not in the PDB file, guess the element:
99 continue
   do i=1,natom
      if (atom_element(i) == '????') then
         atom_element(i) = guess_element(atom_name(i))
      end if
   end do
   
   write(stdout,*) 'PRMTOP: num. atoms=',natom,', num. residues=',nres
   write(stdout,*) 'PDB:    num. atoms=',pdb_natom,', num. residues=',pdb_nres
   
   !---------------------------------------------------------------------------
   ! Copy infile content to outfile (crude with Fortran)
   rewind(in_lun)
   do
      read(in_lun,'(A)',end=1) buf
      write(out_lun,'(A)') buf(1:80)
   end do
   1 continue
   close(in_lun)
   
   !---------------------------------------------------------------------------
   ! Append new data to outfile
   
   fmt='(15I5)'
   write(out_lun,'(A)') &
         '%FLAG RESIDUE_NUMBER', &
         '%COMMENT Residue number (resSeq) read from PDB file; DIMENSION(NRES)'
   if (.not. any(residue_icode /= ' ')) then
       write(out_lun,'(A)') &
          '%COMMENT Residue insertion code (iCode) not present in PDB file', &
          '%COMMENT If present, look for: %FLAG RESIDUE_ICODE'
   endif
   write(out_lun,'(A)') '%FORMAT'//trim(fmt)
   write(out_lun,fmt) residue_number
   
   fmt='(80A1)'
   write(out_lun,'(A)') &
         '%FLAG RESIDUE_CHAINID', &
         '%COMMENT Residue chain ID (chainId) read from PDB file; DIMENSION(NRES)', &
         '%FORMAT'//fmt
   write(out_lun,fmt) residue_chainid
   
   if (any(residue_icode /= ' ')) then
      fmt='(80A1)'
      write(out_lun,'(A)') &
         '%FLAG RESIDUE_ICODE', &
         '%COMMENT Residue insertion code (iCode) read from PDB file; DIMENSION(NRES)', &
         '%FORMAT'//fmt
      write(out_lun,fmt) residue_icode
   end if
   
   fmt='(20A4)'
   write(out_lun,'(A)') &
         '%FLAG ATOM_ELEMENT', &
         '%COMMENT Atom element and formal charge, read from PDB file; DIMENSION(NATOM)', &
         '%FORMAT'//fmt
   write(out_lun,fmt) atom_element
   
   fmt='(80A1)'
   write(out_lun,'(A)') &
         '%FLAG ATOM_ALTLOC', &
         '%COMMENT Atom alternate location, read from PDB file; DIMENSION(NATOM)', &
         '%FORMAT'//fmt
   write(out_lun,fmt) atom_altloc
   
   fmt='(10F8.2)'
   write(out_lun,'(A)') &
         '%FLAG ATOM_OCCUPANCY', &
         '%COMMENT occupancies read from PDB file; DIMENSION(NATOM)', &
         '%FORMAT'//fmt
   write(out_lun,fmt) atom_occupancy
   
   fmt='(10F8.2)'
   write(out_lun,'(A)') &
         '%FLAG ATOM_BFACTOR', &
         '%COMMENT B-factor read from PDB file; DIMENSION(NATOM)', &
         '%FORMAT'//fmt
   write(out_lun,fmt) atom_bfactor
   
#if 0
   fmt='(15A5)'
   write(out_lun,'(A)') &
         '%FLAG ATOM_NUMBER', &
         '%COMMENT atom serial number read from PDB file; DIMENSION(NATOM)', &
         '%FORMAT'//fmt
   write(out_lun,fmt) atom_number
#endif
   
   close(out_lun)
   
contains
   subroutine get_next_arg(string)
      
      implicit none
      
      character(len=*), intent(out) :: string
      if (iarg >= argc) then
         write(*,*)'Error: option "',trim(arg),'" requires an argument.'
         call usage()
      end if
      iarg = iarg + 1
      call getarg(iarg,string)
   end subroutine get_next_arg

   subroutine usage()
      
      implicit none
      
      write(*,'(A)') &
            'add_pdb -i prmtop -p pdb -o prmtop [-guess]', &
            '   prmtop  : amber topology', &
            '   pdb     : matching PDB file', &
            '  -guess   : Guess atomic elements when absent from the PDB file.', &
            '             (default assumes proper element-aligned names)', &
            '', &
            'Residues not in the PDB file are assigned chainID="*" and resSeq=0'
      stop
   end subroutine usage

   function find_atom(res_index,name,altLoc) result(atom_index)
      
      implicit none
      
      integer :: atom_index
      integer, intent(in) :: res_index
      character(len=4), intent(in) :: name
      character(len=1), intent(in) :: altLoc
      ! local
      integer :: i
      ! begin
      do i = residue_pointer(res_index), residue_pointer(res_index+1)-1
         if (atom_name(i)==name) then
            if(altLoc==' ' .or. altLoc=='A') then
               atom_index = i
            elseif(altLoc=='B') then 
               atom_index = i+1
            elseif(altLoc=='C') then 
               atom_index = i+2
            elseif(altLoc=='D') then 
               atom_index = i+3
            elseif(altLoc=='E') then 
               atom_index = i+4
            else
               write(6,*) 'Bad altLoc: ', res_index, name, altLoc
               stop
            end if
            return
         end if
      end do
      atom_index=-1
   end function find_atom

   ! TODO: estimate formal charges?
   function guess_element(name) result(element)
      
      implicit none
      
      character(len=4), intent(in) :: name
      character(len=4) :: element
      ! locals
      character(len=4) :: lname
      ! begin
      lname = adjustl(name)

      ! Leading non-alpha is a wrapped 5th character; insert a leading space.
      if (lname(1:1) < 'A' .or. lname(1:1) > 'Z') then
         element = ' '//lname(2:2)//'  '

      ! If second character is lower case, it is a 2-letter element.
      ! (Not standard, but it would make things a lot easier if it was.)
      else if (lname(1:1) >= 'a' .and. lname(1:1) <= 'z') then
         element = lname(1:2)//'  '

      ! If guessing, prefer the 1-letter element if it is valid.
      else if (index("BCFHIKNOPSUVWY",lname(1:1))>0) then
         element = ' '//lname(1:1)
!        if (lname(2:2)>='A' .and. lname(2:2)<='Z') then
!           if (any(element_name == lname(1:2))) then
!              write(*,*) &
!                 'Warning: guessing atom "',name,'" as element "',element,'"'
!           end if
!        end if

      !Otherwise it must be a 2-letter element.
      else
         element = lname(1:2)
      end if
   end function guess_element
end program add_pdb
