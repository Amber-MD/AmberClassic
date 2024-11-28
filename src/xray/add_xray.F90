program add_xray
!-------------------------------------------------------------------------------
! This program adds data entries to the PRMTOP file needed for X-ray refinement.
! The current implementation is a bit incomplete (i.e. no log output).
! Only P1 symmetry will be supported for now, so the symmetry info is not really
! needed.
!
! Currently, only single-letter elements are supported, because the input PRMTOP
! files do not include atomic element information.
!
! Data appended to the PRMTOP:
!  %FLAG XRAY_ATOM_SCATTER_TYPE_INDEX    ! dimension(NATOM); index into
!                                          XRAY_SCATTER_COEFFICIENTS
!  %FLAG XRAY_NUM_SCATTER_TYPES          ! dimension(1); scalar value defining
!                                          the scatter-factor table size
!  %FLAG XRAY_ATOM_SCATTER_TYPE          ! dimension(NATOM); element-name for
!                                          this scatter type
!  %FLAG XRAY_SCATTER_COEFFICIENTS       ! dimension(2,5,XRAY_NUM_SCATTER_TYPES)
!
!  %FLAG XRAY_NUM_SYMMOPS                ! dimension(1); scalar dimension
!                                          for XRAY_SYMMOPS
!  %FLAG XRAY_SYMMOPS                    ! dimension(3,4,XRAY_NUM_SYMMOPS)
!                                          Each symmop is a 3x3 translation and
!                                          a 3x1 translation, as a 3x4 matrix.
!-------------------------------------------------------------------------------
   implicit none
   integer, parameter :: real_kind=8, rk=real_kind
   integer, parameter :: PATH_MAX=256
   integer, parameter :: stdin=5, stdout=6, stderr=0
   
   integer, parameter :: scatter_ncoeffs = 5 ! 4 gaussians + 1 constant.
   integer, parameter :: infile_lun=10, scatter_lun=11, &
         symmop_lun=12, outfile_lun=13
   
   integer, parameter :: MAX_SYMMOPS=96
   
   character(len=PATH_MAX) :: &
         symmop_filename, &
         scatter_filename, &
         prmtop_infile, &
         prmtop_outfile
   
   character(len=16) :: spacegroup_name='P1'
   
   integer, parameter :: MAX_SCATTER_TYPES = 32
   integer :: num_scatter_types = 0
   character(len=4) :: scatter_types(MAX_SCATTER_TYPES)
   
   integer :: i, j, iok, num_atoms, num_symmops
   integer, allocatable :: scatter_type_index(:)
   character(len=4), allocatable :: atom_name(:), atom_element(:)
   real(real_kind), allocatable :: scatter_coefficients(:,:,:)
   real(real_kind) :: symmop(3,4,MAX_SYMMOPS)
   character(len=32) :: fmt
   ! Current PRMTOP files are limited to 80-character lines.
   ! This allows a bit more for future versions, or long comments.
   character(len=128) :: buf
   
   character(len=32) :: arg
   character(len=PATH_MAX) :: ambhome
   integer :: argc, iarg, amblength
   integer :: iargc
   logical :: electronsf = .false.

   ! begin

   call get_environment_variable( 'AMBERCLASSICHOME', ambhome, amblength )
   if( amblength > 0 ) then
     symmop_filename = ambhome(1:amblength) // '/dat/xray/symop.lib'
   else
     symmop_filename = ' '
   endif
   if( amblength > 0 ) then
      if (electronsf) then
         scatter_filename = ambhome(1:amblength) // &
            '/dat/xray/atomsf_electron.lib'
      else
         scatter_filename = ambhome(1:amblength) // '/dat/xray/atomsf.lib'
      endif
   else
      scatter_filename = ' '
   endif

   argc = command_argument_count()
   if (argc==0) call usage()
   iarg = 0
   do while(iarg<argc)
      iarg=iarg+1
      call getarg(iarg,arg)
      select case(trim(arg))
      case('-h')
         call usage()
      case('-i')
         call get_next_arg(prmtop_infile)
      case('-o')
         call get_next_arg(prmtop_outfile)
      case('-symm')
         call get_next_arg(symmop_filename)
      case('-sg')
         call get_next_arg(spacegroup_name)
      case('-electronsf')
         electronsf = .true.
      case default
         write(*,*) 'Unrecognize option: ',trim(arg)
         call usage()
      end select
   end do
   
   call get_environment_variable( 'AMBERCLASSICHOME', ambhome, amblength )
   if( amblength .le. 0 ) then
      write(0,*) 'Error: AMBERCLASSICHOME is not set!'
      stop
   endif

   symmop_filename = ambhome(1:amblength) // '/dat/xray/symop.lib'
   if (electronsf) then
      scatter_filename = ambhome(1:amblength) // '/dat/xray/atomsf_electron.lib'
   else
      scatter_filename = ambhome(1:amblength) // '/dat/xray/atomsf.lib'
   endif

   open(unit=infile_lun,file=prmtop_infile,status='OLD', &
         action='READ',form='FORMATTED')
   open(unit=outfile_lun,file=prmtop_outfile,status='UNKNOWN', &
         action='WRITE',form='FORMATTED')
   open(unit=scatter_lun,file=scatter_filename,status='OLD', &
         action='READ',form='FORMATTED')
   
   call nxtsec(infile_lun,STDOUT,0,'*','POINTERS',fmt,iok)
   read(infile_lun,fmt) num_atoms
   write(stdout,*) 'Number of atoms = ',num_atoms
   
   allocate( &
         atom_name(num_atoms), atom_element(num_atoms), &
         scatter_type_index(num_atoms))
   
   call nxtsec(infile_lun,STDOUT,0,'*','ATOM_NAME',fmt,iok)
   read(infile_lun,fmt) atom_name

   call nxtsec(infile_lun,STDOUT,0,'*','ATOM_ELEMENT',fmt,iok)
   read(infile_lun,fmt) atom_element

   ! scattering tables use element names like Fe, not FE; convert
   !   second character to lower case if it is upper case to begin with:
   do i = 1,num_atoms
      if( atom_element(i)(1:1) .ne. ' ' .and.  &
          iachar(atom_element(i)(2:2)) .le. 90 ) &
          atom_element(i)(2:2) = achar(iachar(atom_element(i)(2:2)) + 32)
   end do

   ! fix up ion elements to match those in X-ray scattering tables:
   if( .not. electronsf ) then
      do i = 1,num_atoms
         if( atom_element(i) == 'Ca  ' ) atom_element(i) = 'Ca+2'
         if( atom_element(i) == 'Cd  ' ) atom_element(i) = 'Cd+2'
         if( atom_element(i) == 'Mg  ' ) atom_element(i) = 'Mg+2'
         if( atom_element(i) == 'Na  ' ) atom_element(i) = 'Na+1'
         if( atom_element(i) == ' K  ' ) atom_element(i) = ' K+1'
         if( atom_element(i) == 'Cl  ' ) atom_element(i) = 'Cl-1'
         if( atom_element(i) == 'Fe  ' ) atom_element(i) = 'Fe+2'
         if( atom_element(i) == 'Zn  ' ) atom_element(i) = 'Zn+2'
         if( atom_element(i) == 'Ep  ' ) atom_element(i) = ' H  '
      end do
   end if

   num_scatter_types=1
   scatter_types(1)=atom_element(1)
   do i=2,num_atoms
      if (any(scatter_types(1:num_scatter_types)==atom_element(i))) cycle
      if (num_scatter_types == MAX_SCATTER_TYPES) then
         stop 'ERROR: Too many scatter types (bug?)'
      end if
      num_scatter_types = num_scatter_types + 1
      scatter_types(num_scatter_types) = atom_element(i)
   end do

   scatter_type_index(:) = 0
   do i=1,num_atoms
      do j=1,num_scatter_types
         if (scatter_types(j)==atom_element(i)) then
            scatter_type_index(i) = j
            exit
         end if
      end do
   end do
   
   allocate(scatter_coefficients(2,scatter_ncoeffs,num_scatter_types))
   call read_ccp4_atomsf(scatter_lun,num_scatter_types, &
         scatter_types,scatter_coefficients)
   close(scatter_lun)
   
   if (symmop_filename /= ' ') then
      open(unit=symmop_lun,file=symmop_filename,status='OLD', &
            action='READ',form='FORMATTED')
      call read_ccp4_symlib(symmop_lun,trim(spacegroup_name),num_symmops,symmop)
      close(symmop_lun)
   end if
   
   !---------------------------------------------------------------------------
   ! Copy infile content to outfile (crude with Fortran)
   rewind(infile_lun)
   do
      read(infile_lun,'(A)',end=1) buf
      write(outfile_lun,'(A)') buf(1:80)
   end do
   1 continue
   close(infile_lun)
   
   !---------------------------------------------------------------------------
   ! Append new data to outfile
   
   fmt='(I4)'
   write(outfile_lun,'(A)') '%FLAG XRAY_NUM_SCATTER_TYPES'
   write(outfile_lun,'(A)') '%COMMENT INTEGER, DIMENSION(1)'
   write(outfile_lun,'(2A)') '%FORMAT',fmt
   write(outfile_lun,fmt) num_scatter_types
   
   fmt='(20I4)'
   write(outfile_lun,'(A)') '%FLAG XRAY_ATOM_SCATTER_TYPE_INDEX'
   write(outfile_lun,'(A)') '%COMMENT REAL, DIMENSION(NUM_ATOMS)'
   write(outfile_lun,'(2A)') '%FORMAT',fmt
   write(outfile_lun,fmt) scatter_type_index

   fmt='(20A4)'
   write(outfile_lun,'(A)') '%FLAG XRAY_SCATTER_TYPE'
   write(outfile_lun,'(A)') '%COMMENT REAL, DIMENSION(XRAY_NUM_SCATTER_TYPES)'
   write(outfile_lun,'(2A)') '%FORMAT',fmt
   write(outfile_lun,fmt) scatter_types(1:num_scatter_types)
   
   fmt='(5F12.6)'
   write(outfile_lun,'(A)') '%FLAG XRAY_SCATTER_COEFFICIENTS'
   write(outfile_lun,'(A,I1,A)') '%COMMENT REAL, DIMENSION(2,', &
         scatter_ncoeffs,',XRAY_NUM_SCATTER_TYPES)'
   write(outfile_lun,'(2A)') '%FORMAT',fmt
   write(outfile_lun,fmt) scatter_coefficients
   
#if 0
   fmt='(I4)'
   write(outfile_lun,'(A)') '%FLAG XRAY_NUM_SYMMOPS'
   write(outfile_lun,'(A)') '%COMMENT INTEGER, DIMENSION(1)'
   write(outfile_lun,'(2A)') '%FORMAT',fmt
   write(outfile_lun,fmt) num_symmops
   
   fmt='(3F12.8)'
   write(outfile_lun,'(A)') '%FLAG XRAY_SYMMOPS'
   write(outfile_lun,'(A)') '%COMMENT REAL, DIMENSION(3,4,XRAY_NUM_SYMMOPS)'
   write(outfile_lun,'(A)') &
         '%COMMENT Each symmop is a 3x3 rotation matrix and a 3x1 translation'
   write(outfile_lun,'(2A)') '%FORMAT',fmt
   write(outfile_lun,fmt) symmop(:,:,1:num_symmops)
#endif
   
   close(outfile_lun)
   
contains
   subroutine get_next_arg(string)
      character(len=*), intent(out) :: string
      if (iarg >= argc) then
         write(*,*)'Error: option "',trim(arg),'" requires an argument.'
         call usage()
      end if
      iarg = iarg + 1
      call getarg(iarg,string)
   end subroutine get_next_arg
   subroutine usage()
      write(*,'(A)') &
            'add_xray -i prmtop -o prmtop -symm symlib -sf atomsf -sg sgname', &
            '   prmtop  : amber topology, input and output files', &
            '   symlib  : symmetry operator file, CCP4 format (symop.lib)', &
            '   atomsf  : scatter factor file, CCP4 format (atomsf.lib)', &
            '   sgname  : spacegroup name (P1)', &
            ''
      stop
   end subroutine usage
   ! --------------------------------------------------------------------------
   ! Lookup and calculation of standard gaussian-sum atomic structure factors.
   !
   ! COEFF(2,num_scatter_coeffs) is a merged array of gaussian A,B coefficients.
   ! The lst coefficent is the constant term. With B=0, A is the same as the
   ! constant C because exp(0)==1
   !
   ! Standard coefficient tables use 4 gaussians plus a constant term.
   !
   !  [ A1, B1 ]
   !  [ A2, B2 ]
   !   ...  ...
   !  [ An, Bn ]
   !  [ C , 0  ]
   subroutine read_ccp4_atomsf(lun,num_types,names,coeffs)
      implicit none
      integer, intent(in) :: lun
      integer, intent(in) :: num_types
      character(len=4), intent(in) :: names(num_types)
      real(real_kind), intent(out) :: coeffs(2,scatter_ncoeffs,num_types)
      
      integer :: ios
      integer :: n
      integer :: nelec
      character(len=1) :: ngauss
      character(len=4) :: atom_type
      real(real_kind) :: coeff(2,scatter_ncoeffs)
      real(real_kind) :: fp,fpp
      character(len=90) :: line
      
      coeffs(2,scatter_ncoeffs,:) = -1.0_rk ! undefined value flag
      do
#if 0   /* later turn this into a program option */
         read(lun,'(a)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:1) == '#') cycle
         read(line, * ) atom_type, nelec, coeff(1,:), coeff(2,:)
#else  /* xray form factors are read here: */
         read(lun,'(A4,1X,A1)',iostat=ios) atom_type, ngauss
         if (ios /= 0) exit
         if (atom_type(1:1) == ' ' .or. atom_type(1:2) == 'AD'  &
               .or. ngauss == '2') cycle
         ! The two skipped 8X fields are actually integers for atomic weight
         ! and number of electrons. (Currently not used here.)
         read (lun,'(2(2X,8X),2X,F14.6)',err=1) coeff(1,5) ! constant term
         coeff(2,5) = 0.0_rk ! constant term has zero exponent
         read (lun,'(4(2X,F14.6))',err=1) coeff(1,1:4)
         read (lun,'(4(2X,F14.6))',err=1) coeff(2,1:4)
         read (lun,'(2(2X,F14.6))',err=1) fp,fpp ! wavelength-dependent f' and f''
#endif
         do n = 1,num_types
            if (adjustl(names(n)) == atom_type) then
               coeffs(:,:,n) = coeff
               exit
            end if
         end do
         
      end do
      close(lun)
      
      do n = 1,num_types
         if (coeffs(2,scatter_ncoeffs,n) == -1.0_rk) then
            coeffs(:,:,n) = 0.0_rk
            write(*,*) 'ERROR: No scatter for type(',names(n),')'
            call exit(1)
         end if
      end do
      return
   1  continue
      write(*,*) 'ERROR while parsing atomsf.lib'
      call exit(1)
   end subroutine read_ccp4_atomsf
   
   subroutine read_ccp4_symlib(symlib,spacegroup_name,num_symmops,symmop)
      implicit none
      integer, intent(in) :: symlib
      character(len=*), intent(in) :: spacegroup_name
      integer, intent(out) :: num_symmops
      real(real_kind) :: symmop(3,4,MAX_SYMMOPS)
      
      character(len=80) :: line
      integer :: ios
      
      integer :: i, i1, i2, op, p, axis
      real(real_kind) :: sign
      
      integer :: spacegroup_number
      
      ! These hold the data parsed for each spacegroup.
      ! Example:    4 2 2 P21 PG2 MONOCLINIC 'P 1 21 1'
      integer :: sg_number               ! 4
      integer :: nsymm                   ! 2
      integer :: nprim_symm              ! 2
      character(len=7) :: sg_shortname   ! P21
      character(len=6) :: pointgroup     ! PG2
      character(len=12) :: xtal_system   ! MONOCLINIC
      character(len=16) :: sg_longname   ! P 1 21 1
      character(len=16) :: sg_longname_condensed   ! P1211
      
      ! begin
      ! Allow for the spacegroup_name to be given as a number.
      read(spacegroup_name,*,iostat=ios) spacegroup_number
      if (ios /= 0) spacegroup_number=0
      
      ios=0
      do while (ios==0)
         read (symlib,'(A80)',iostat=ios) line
         if (ios /= 0) then
            write(stderr,*) &
                  'Cannot find symmetry operators for spacegroup "', &
                  spacegroup_name,'" !'
            call exit(1)
         end if
         if (line(1:1)>='1' .and. line(1:1)<='9') then
            read(line,*,iostat=ios) sg_number, nsymm, nprim_symm, sg_shortname, &
                  pointgroup, xtal_system, sg_longname
            if (ios==0) then
               ! Convert all-caps xtal_system to capitalized.
               do i=2,Len_Trim(xtal_system)
                  xtal_system(i:i)=Char(ibset(IChar(xtal_system(i:i)),5))
               end do
               
               ! Remove spaces from sg_longname, saved in sg_longname_condensed
               i1=1
               do i2=1,Len_Trim(sg_longname)
                  if (sg_longname(i2:i2)/=' ') then
                     sg_longname_condensed(i1:i1)=sg_longname(i2:i2)
                     i1=i1+1
                  end if
               end do
               
               if (sg_number==spacegroup_number &
                     .or. sg_longname==spacegroup_name &
                     .or. sg_shortname==spacegroup_name &
                     .or. sg_longname_condensed==spacegroup_name) then
                 exit
               end if
            end if
         end if
      end do
      
      num_symmops = nsymm
      
      ! Parse symmetry operator strings of the form: -X,1/2+Y,1/2+Z
      ! This is easier than it sounds. Each comma-delimited part represents
      ! one row of a 3x3 rotation matrix, and one element of a 3x1 translation
      ! vector. An X, Y or Z puts -1 or 1 in the corresponding column of the
      ! current row. Fractions are evaulated and added to the translation
      ! vector. Here, symmops are stored as a 3x4 matrix, with the last row
      ! being the translation vector.
      do op=1,nsymm
         read(symlib,'(A)') line
         p=0
         do axis = 1,3
            p=p+1
            sign=1.0
            do while (line(p:p)/=',' .and. p<=Len_Trim(line))
               select case(line(p:p))

               ! Check for a sign, addition or subtraction.
               ! A '+' has no effect. Here, proper syntax is not enforced, so
               ! a missing operator is treated as an implicit '+'.
               case ('-')
                  sign=-1.0
               case ('+',' ')
                  continue

               ! If a fraction is found, add it to the translation vector.
               ! Only these fractions are valid: 1/2, 3/4, 1/4, 1/3, 2/3, 1/6, 5/6
               case ('1','2','3','5')
                  if (line(p+1:p+1)/='/' .or. index('2346',line(p+1:p+2))==0) then
                     write(stderr,'(A)') &
                           'ERROR parsing CCP4 SYMOP library: unexpected digit.', &
                           'All digits must be part of a fraction: [123]/[2346]'
                     call exit(1)
                  end if
                  symmop(axis,4,op)=sign * real(IChar(line(p:p))-IChar('0')) &
                        / real(IChar(line(p+2:p+2))-IChar('0'))
                  p=p+2
                  sign=1.0

               ! For X,Y,Z, set the matching part of the current matrix row.
               case ('x','X')
                  symmop(1,axis,op)=sign
                  sign=1.0
               case ('y','Y')
                  symmop(2,axis,op)=sign
                  sign=1.0
               case ('z','Z')
                  symmop(3,axis,op)=sign
                  sign=1.0

               case default
                  write(stderr,'(3A)') &
                        'Error parsing CCP4 SymLib expression: "',trim(line),'"'
                  call exit(1)
               end select
               p=p+1
            end do ! p  (character position)
         end do ! axis=1,3
      end do ! op=1,nsymm
   end subroutine read_ccp4_symlib
   
end program add_xray
