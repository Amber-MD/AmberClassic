#include "../include/assert.fh"
! This is a module to encapsulate all of the AMBER common blocks and
! other global junk into a module namespace.
module xray_common_module
use file_io_dat
end module xray_common_module

module xray_cpu_module
   use xray_globals_module
   use xray_contracts_module
   use xray_bulk_mask_data_module, only: k_sol, b_sol
   use xray_interface_pre_init_data, only: num_scatter_types => n_scatter_types
   implicit none
   private

   public :: finalize
   public :: init
   public :: xray_get_derivative
   public :: xray_read_mdin
   public :: xray_read_parm
   public :: xray_write_md_state
   public :: xray_write_min_state
   public :: xray_write_options
   public :: xray_write_pdb
   public :: xray_write_fmtz

   namelist /xray/ &
         pdb_infile, pdb_outfile, &
         fave_outfile, fmtz_outfile, &
         pdb_read_coordinates, &
         pdb_use_segid, &
         pdb_wrap_names, &
         spacegroup_name, spacegroup_number, &
         reflection_infile, &
         xray_weight, &
         target, &
         solvent_mask_probe_radius, &
         solvent_mask_adjustment, &
         ntwsf, &
         sf_outfile, &
         atom_selection_mask, &
         k_sol, b_sol, ls_r3, ls_r4, ixp, iyp, izp, &
         mask_update_period, scale_update_period, &
         ml_update_period, bulk_solvent_model,  &
         timeavg
   
   !-------------------------------------------------------------------
contains

   subroutine xray_read_mdin(mdin_lun)
      implicit none
      integer, intent(in) :: mdin_lun
#include "nmr.h"

      character(len=512) :: line
      integer :: stat
      
      call xray_init_globals()
      if (.not.xray_active) then
        return
      end if

      rewind(mdin_lun)
      read(unit=mdin_lun,nml=xray,iostat=stat)

      if (stat /= 0) then
        backspace(mdin_lun)
        read(mdin_lun, fmt='(A)') line
        write(stdout, '(A)') 'Invalid line in &xray namelist '//': '//trim(line)
        call mexit(stdout,1)
      end if

      ! Check input
      if (ntwsf < 0) then
        write(stdout, '(A)') 'ntwsf must be >= 0'
        call mexit(stdout,1)
      end if
      
      ! initialize the weight-change value, in case other weight changes
      ! might be used:
      wxray = xray_weight
      
      !write(unit=6,nml=xray)
   end subroutine xray_read_mdin

   subroutine xray_write_options()
      use xray_bulk_mask_data_module, only: k_sol, b_sol
      implicit none

      write(stdout,'(/,A)') 'X-ray Refinement Parameters:'
      write(stdout,'(5X,2A)') 'PDB InFile: ',trim(pdb_infile)
      if( pdb_outfile /= '' ) &
         write(stdout,'(5X,2A)') 'PDB OutFile:',trim(pdb_outfile)
      if( fave_outfile /= '' ) &
         write(stdout,'(5X,2A)') 'FCALC_AVE OutFile:',trim(fave_outfile)
      if( fmtz_outfile /= '' ) &
         write(stdout,'(5X,2A)') 'FMTZ OutFile:',trim(fmtz_outfile)
      write(stdout,'(5X,A,L1)') 'PDB Read Coordinates: ',pdb_read_coordinates
      write(stdout,'(5X,A,L1)') 'PDB Use SegID: ',pdb_use_segid
      write(stdout,'(5X,A,L1)') 'PDB Wrap Names: ',pdb_wrap_names
      write(stdout,'(5X,3A,i3)') 'Spacegroup: ',trim(spacegroup_name), &
          ' number ',spacegroup_number
      write(stdout,'(5X,2A)') 'Reflection InFile: ',trim(reflection_infile)
      write(stdout,'(5X,A,E10.3)') 'X-ray weight: ', xray_weight
      if( target(1:2) .eq. 'ls' ) then
         write(stdout,'(5X,A,A4,2F8.3)') 'Use target: ',target, ls_r3, ls_r4
      else
         write(stdout,'(5X,A,A4)') 'Use target: ',target
      endif

      write(stdout,'(5X,A,I5)') 'Scale update Interval: ',scale_update_period
      ! write(stdout,'(5X,A,F8.3)') 'Solvent mask probe radius: ',solvent_mask_probe_radius
      ! write(stdout,'(5X,A,F8.3)') 'Solvent mask expand: ',solvent_mask_expand
      ! write(stdout,'(5X,2A)') 'Solvent Mask OutFile:',trim(solvent_mask_outfile)
      ! write(stdout,'(5X,2A)') 'Solvent Mask Reflection OutFile:',trim(solvent_mask_reflection_outfile)
      write(stdout,'(5X,A,I5)') 'Solvent Mask Update Interval: ',mask_update_period
      write(stdout,'(5X,2(A,F8.3))') 'Solvent scale:',k_sol,', B-factor:', b_sol
      ! write(stdout,'(5X,2(A,F8.3))') 'B-Factor Min:',bfactor_min,', Max: ',bfactor_max
      ! write(stdout,'(5X,A,I4)') 'B-factor Refinement Interval: ',bfactor_refinement_interval
      write(stdout,'(5X,2A)') 'Atom Selection Mask: ',trim(atom_selection_mask)
      return
   end subroutine xray_write_options

   ! Read X-ray data from the PRMTOP file, and also save pointers to global
   ! PRMTOP data.
   subroutine xray_read_parm(prmtop_lun, out_lun)
      use memory_module, only: natom, nres
      use xray_interface_pre_init_data, only: n_scatter_coeffs, scatter_coefficients
      implicit none
      integer, intent(in) :: prmtop_lun, out_lun
      ! local
      character(len=32) :: fmt
      integer :: ierr

      ! if (pdb_outfile /= '') then
         num_atoms = natom
         num_residues = nres

         allocate(atom_bfactor(natom), atom_occupancy(natom), &
            atom_selection(natom), residue_chainid(nres), residue_icode(nres), &
            atom_element(natom), atom_altloc(natom), residue_number(nres))

         call nxtsec_reset()
         call nxtsec(prmtop_lun,STDOUT,0,'(20I4)','RESIDUE_NUMBER',fmt,ierr)
         read(prmtop_lun,fmt) residue_number
         call nxtsec(prmtop_lun,STDOUT,0,'(20A4)','RESIDUE_CHAINID',fmt,ierr)
         read(prmtop_lun,fmt) residue_chainid

         call nxtsec(prmtop_lun,STDOUT,1,'*','RESIDUE_ICODE',fmt,ierr)
         if (ierr==0) then
            read(prmtop_lun,fmt) residue_icode
         else
            residue_icode=' '
         end if

         call nxtsec(prmtop_lun,STDOUT,0,'(20A4)','ATOM_ELEMENT',fmt,ierr)
         read(prmtop_lun,fmt) atom_element

         call nxtsec(prmtop_lun,STDOUT,1,'*','ATOM_ALTLOC',fmt,ierr)
         if (ierr==0) then
            read(prmtop_lun,fmt) atom_altloc
         else
            atom_altloc=' '
         end if

      ! end if

      if (reflection_infile == '') return

      call nxtsec(prmtop_lun,out_lun,0,'(I4)','XRAY_NUM_SCATTER_TYPES',fmt,ierr)
      if (fmt=='*') then
         write(stdout,'(A)') &
            'ERROR: XRAY_NUM_SCATTER_TYPES not found in PRMTOP file.'
         call mexit(stdout,1)
      end if
      read(prmtop_lun,fmt) num_scatter_types

      allocate(atom_scatter_type(natom),  &
            scatter_coefficients(2, n_scatter_coeffs, num_scatter_types))

      call nxtsec(prmtop_lun,out_lun,0,'(20I4)','XRAY_ATOM_SCATTER_TYPE_INDEX',fmt,ierr)
      read(prmtop_lun,fmt) atom_scatter_type
      call nxtsec(prmtop_lun,out_lun,0,'(F12.6)','XRAY_SCATTER_COEFFICIENTS',fmt,ierr)
      read(prmtop_lun,fmt) scatter_coefficients
      call nxtsec(prmtop_lun,out_lun,1,'*','XRAY_SYMMETRY_TYPE',fmt,ierr)
      if (ierr==-2) then
         write(STDOUT,*) &
               'XRAY_SYMMETRY_TYPE not found in PRMTOP file; assuming P1'
         num_symmops = 1
         au_type = 1
      else
         stop 'ONLY P1 SUPPORTED FOR NOW'
         read(prmtop_lun,fmt) num_symmops, spacegroup_number, au_type, system
         call nxtsec(prmtop_lun,out_lun,1,'*', &
               'XRAY_SYMMETRY_OPERATORS',fmt,ierr)
         ! ...
      end if
   end subroutine xray_read_parm

   subroutine xray_read_pdb(filename)
      use memory_module, only: residue_label,atom_name,coordinate
      implicit none
      character(len=*), intent(in) :: filename
      ! locals
      character(len=4) :: name,segID,element
      character(len=3) :: resName
      character(len=1) :: altLoc,chainID,iCode
      integer :: resSeq
      real(real_kind) :: xyz(3),occupancy,tempFactor
      character(len=80) :: line
      integer :: unit, iostat, iatom, ires, i, j, ndup, nmiss
      real(real_kind), parameter :: MISSING = -999.0_rk_
      ! begin
      atom_occupancy(:)=MISSING
      call amopen(allocate_lun(unit),filename,'O','F','R')
      ndup=0
      iatom=1
      ires=1
      do
         read(unit,'(A)',iostat=iostat) line
         if (iostat/=0) exit
         if (line(1:6)=='END   ') exit
         if (line(1:6)=='ATOM  ' .or. line(1:6)=='HETATM') then
            read(line,'(12X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,2A4)') &
                  name,altLoc,resName,chainID,resSeq,iCode, &
                  xyz,occupancy,tempFactor,segID,element
            i = find_atom(name,resName,chainID,resSeq,iCode,altLoc)
            if (i<0) then
               write(stdout,'(A)') 'Atom not found:'
               write(stdout,'(A)') trim(line)
               stop
            ! else
            !    write(stdout,'(a,a)') 'Found: ', line(1:54)
            end if
            if (atom_occupancy(i) >= 0) then
               ndup=ndup+1
               if (ndup<10) then
                  write(stdout,'(3(A,1X),A,I4,A)') 'PDB: Duplicate ATOM:', &
                        name,resName,chainID(1:1),resSeq,iCode(1:1)
               end if
            end if
            if (pdb_read_coordinates) coordinate(1:3,i) = xyz
            atom_bfactor(i) = tempFactor
            atom_occupancy(i) = occupancy
         end if
      end do
      nmiss = count(atom_occupancy==MISSING)
      if (nmiss>0) then
         write(stdout,'(A,I7,A)') 'PDB: missing data for ',nmiss,' atoms.'
         j=0
         do i=1,num_atoms
            if (atom_occupancy(i)==MISSING) then
               atom_occupancy(i)=0
               j=j+1
               if (j<=10) then
                  ! would be better to print more info here, but it's 
                  ! awkward, since find_atom() only returns the atom number.
                  ! Let's hope this is rare, and provides enough info....
                  write(stdout,'(A,I5,1x,A)') &
                     'PDB: Missing ATOM:', i, atom_name(i)
               end if
            end if
         end do
      end if
      if (nmiss==0 .and. ndup==0) then
         write(stdout,'(A)') 'PDB: All atoms read successfully.'
      end if
      close(unit)
      return
   end subroutine xray_read_pdb

   function find_atom(name,resName,chainID,resSeq,iCode,altLoc) result(atom_serial)
      use memory_module, only: residue_pointer,residue_label,atom_name
      implicit none
      integer :: atom_serial
      character(len=4), intent(in) :: name
      character(len=3), intent(in) :: resName
      character(len=1), intent(in) :: chainID, iCode, altLoc
      integer, intent(in) :: resSeq
      ! locals
      character(len=4) :: lname
      integer, save :: ires = 1
      integer :: i,j
      lname = adjustl(name)
      ! first find the matching residue; no need to match resname:
      do i=1,num_residues
         if (        resSeq==mod(residue_number(ires),10000)  &
               .and. chainID==residue_chainid(ires) &
               .and. iCode==residue_icode(ires)) then
            ! then find the matching atom name:
            do j = residue_pointer(ires),residue_pointer(ires+1)-1
               if (lname==atom_name(j) .and. altLoc==atom_altLoc(j)) then
                  atom_serial = j
                  return
               end if
            end do
            ! Continue searching, just in case there is a residue
            ! that has been split into two parts.
         end if
         ires = ires + 1
      end do
      atom_serial = -1  ! not-found indicator
      return
   end function find_atom

   subroutine xray_write_pdb(filename)
      use xray_common_module, only: title, title1
      use memory_module, only: &
           residue_pointer,residue_label,atom_name,coordinate
      implicit none
      character(len=*), intent(in) :: filename
      ! locals
      integer :: unit, iatom, ires, ierr
      integer :: first1, last1, first2, last2, ibond
      integer :: iatom_p, ires_p
      character(len=4) :: name
      character(len=8) :: date
      character(len=10) :: time
      ! character(len=1) :: altloc
      ! character(len=4) :: segid
      character(len=3) :: resName
      logical          :: isStandardRes
      character(len=3), parameter :: standard_pdb_residues(28) = (/ &
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE", &
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL", &
            " DG"," DA"," DT"," DC","  G","  A","  U","  C" /)
      ! Amber modres types: "CYX","HID","HIE","HIP",
      character(len=*), parameter :: pdbfmt_MODRES = &
            '("MODRES",1X,A4,1X,A3,1X,A1,1X,I4,A1,1X,A3,2X,A41)'
      ! GMS: Fix for pgf90 compiler
      character(len=4) :: this_residue_chainid

      call amopen(allocate_lun(unit),filename,'U','F','R')
      call date_and_time(date,time)
      if (title/='') write(unit,'(2A)') 'REMARK  ', title
      if (title1/='') write(unit,'(2A)') 'REMARK  ', title1
      write(unit,'(12A)') 'REMARK  Written by MSANDER, ', &
            date(1:4),'.',date(5:6),'.',date(7:8),'  ', &
            time(1:2),':',time(3:4),':',time(5:6)

      ! Actually '(6A,3F9.3A9,3F7.2,1X,A11,I4)', with last value = Z
      write(unit,'(A6,3F9.3,3F7.2,1X,A11)') &
            'CRYST1', unit_cell%as_array(), spacegroup_name

      do ires = 1,num_residues
         do iatom = residue_pointer(ires), residue_pointer(ires+1)-1
            ! ***NOTE***
            ! This code only adds a leading space to give element-alignment
            ! where possible. It is impossible to follow the PDB version 3
            ! "remediated" format correctly, because it has no alignment rules.
            ! Instead, it assumes you have a complete database of all known
            ! residues, and any other residue names are a fatal error.
            name = atom_name(iatom)
            if (atom_element(iatom)(1:1)==' ' &
                  .and. name(1:1) == atom_element(iatom)(2:2)) then
               if (len_trim(name) < 4 .or. pdb_wrap_names) then
                  name = name(4:4)//name(1:3)
               end if
            end if
            resName=residue_label(ires)(1:3)
            resName=adjustr(resName)
            ! GMS: Fix for pgf90 compiler
            this_residue_chainid = residue_chainid(ires)
            ! DRR: PGI does not seem to like any() inside merge() intrinsic.
            isStandardRes = any(resName==standard_pdb_residues)
            ! don't overflow atom or residue numbers:
            iatom_p = mod( iatom, 100000 )
            ires_p = mod( residue_number(ires), 10000 )
            write(unit,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2)')&
                  merge('ATOM  ', 'HETATM', isStandardRes), &
                  iatom_p,name,atom_altloc(iatom)(1:1), &
                  resName,residue_chainid(ires)(1:1), &
                  ires_p,residue_icode(ires)(1:1), &
                  coordinate(1:3,iatom), &
                  atom_occupancy(iatom), &
                  atom_bfactor(iatom), &
                  atom_element(iatom)(3:4)
         end do
      end do
      write(unit,'(A)') 'END'
      close(unit)
      return
   end subroutine xray_write_pdb

   subroutine xray_write_fmtz(filename)

   use xray_globals_module
   use xray_interface2_data_module, only:  Fcalc, Fobs, hkl, resolution, &
       sigma_Fobs, new_order, penalty
   use xray_target_module, only : target_function_id
   implicit none
   character(len=*), intent(in) :: filename
   integer rfree(num_hkl)
   real(real_kind) :: phicalc, zero=0.0
   complex(real_kind) :: Fcalcav
   integer :: i

   rfree = test_flag(new_order)

   open(20,file=trim(fmtz_outfile),action='write')
   if( target_function_id == 1 ) then
      write(20,'(15a)') 'h',achar(9),'k',achar(9),'l',achar(9), &
         'resolution', achar(9), 'fobsr',achar(9),'fobsi',achar(9), &
         'fcalcr',achar(9),'fcalci'
      write(20,'(15a)') '4N',achar(9),'4N',achar(9),'4N',achar(9), &
         '15N', achar(9), '15N',achar(9), '15N',achar(9),'15N', achar(9),'15N'
      do i=1,num_hkl
         write(20,&
          '(i4,a,i4,a,i4,a,f8.3,a,f12.3,a,f12.3,a,f12.3,a,f12.3)') &
          hkl(1,i), achar(9),hkl(2,i), achar(9), hkl(3,i), achar(9), &
          resolution(i), achar(9), real(Fobs(i)), achar(9), &
          aimag(Fobs(i)), achar(9), real(Fcalc(i)), achar(9), aimag(Fcalc(i))
      end do
   else
      write(20,'(19a)') 'h',achar(9),'k',achar(9),'l',achar(9), &
         'resolution', achar(9), 'fobs',achar(9),'sigfobs',achar(9), &
         'fcalc',achar(9),'phicalc',achar(9),'fcalcr', &
         achar(9),'fcalci'
      write(20,'(19a)') '4N',achar(9),'4N',achar(9),'4N',achar(9), &
         '15N', achar(9), '15N',achar(9), '15N',achar(9),'15N', &
         achar(9),'15N',achar(9),'15N',achar(9),'15N'
      do i=1,num_hkl
#if 0
       if( mod(i,4) == 1 ) then
         Fcalcav = Fcalc(i)
       else if( mod(i,4) == 2 .or. mod(i,4) == 3 ) then
         Fcalcav = Fcalcav + Fcalc(i)
       else
         Fcalcav = (Fcalcav + Fcalc(i)) / 4.0
       endif
#endif
         phicalc = atan2( aimag(Fcalc(i)), real(Fcalc(i)) ) * 57.2957795d0
         write(20,&
         '(i4,a,i4,a,i4,a,f8.3,a,f12.3,a,f12.3,a,f12.3,a,f12.3,a,f12.3,a,f12.4)') &
          hkl(1,i), achar(9),hkl(2,i), achar(9), hkl(3,i), achar(9), &
          resolution(i), achar(9), abs(Fobs(i)), achar(9), &
          sigma_Fobs(i), achar(9), abs(Fcalc(i)), achar(9), phicalc, &
          achar(9), real(Fcalc(i)),achar(9),aimag(Fcalc(i))
#if 0
         if( mod(i,4) == 0 ) write(20,&
         '(i4,a,i4,a,i4,a,f8.3,a,f12.3,a,f12.3,a,f12.3,a,f12.3,a,i4,a,f12.4)') &
          hkl(1,i), achar(9),hkl(2,i), achar(9), hkl(3,i), achar(9), &
          resolution(i), achar(9), abs(Fobs(i)), achar(9), &
          sigma_Fobs(i), achar(9), abs(Fcalcav), achar(9), phicalc, &
          achar(9), rfree(i),achar(9),zero
#endif
      end do
   end if
   close(20)
   end subroutine xray_write_fmtz
  

   subroutine init()

      use xray_common_module, only: inpcrd
      use AmberNetcdf_mod, only: NC_checkRestart
      use binrestart, only: read_nc_restart_box
      use xray_pure_utils, only: pack_index
      use xray_unit_cell_module, only: unit_cell_t
      use findmask, only: atommask
      use memory_module, only: natom,nres,ih,m02,m04,m06,ix,i02,x,lcrd,i100
      use xray_interface_pre_init_data, only: scatter_coefficients
      use xray_interface2_module, only: init_interface2 => init
      ! use xray_debug_dump_module, only: xray_dump => dump  ! FIXME: remove this line in release
      use constants, only : DEG_TO_RAD
      implicit none
      ! local
      integer :: hkl_lun, i, j, k, NAT_for_mask1
      real(real_kind) :: resolution, fabs_solvent, phi_solvent
      complex(real_kind), allocatable, dimension(:) :: Fobs
      real(real_kind) :: phi,a,b,c,alpha,beta,gamma
      real(real_kind) :: abs_Fuser, phase_Fuser
      integer :: has_Fuser, alloc_status

      if (pdb_infile /= '') call xray_read_pdb(trim(pdb_infile))
      if (reflection_infile == '') xray_active = .false.

      ! get the values for ucell:
      if ( NC_checkRestart(inpcrd) ) then
        write(6,'(a,a)') ' getting box info from netcdf file ',trim(inpcrd)
        call read_nc_restart_box(inpcrd,a,b,c,alpha,beta,gamma)
      else
        write(6,'(a,a)') ' getting box info from bottom of ',trim(inpcrd)
        call peek_ewald_inpcrd(inpcrd,a,b,c,alpha,beta,gamma)
      endif

      write(stdout,'(A,3F9.3,3F7.2)') &
            'XRAY: UNIT CELL= ',a, b, c, alpha, beta, gamma
      call unit_cell%init(a, b, c, alpha, beta, gamma)

      !--------------------------------------------------------------
      ! Read reflection data
      call amopen(allocate_lun(hkl_lun),reflection_infile,'O','F','R')
      read(hkl_lun,*,end=1,err=2) num_hkl, has_Fuser

      allocate(hkl_index(3,num_hkl),abs_Fobs(num_hkl),sigFobs(num_hkl), &
            & test_flag(num_hkl))

      if (fave_outfile /= '') then
         write(stdout,'(A)') 'fave_outfile is not yet implemented'
         call mexit(stdout, 1)
      endif

      !  each line contains h,k,l and two reals
      !  if target == "ls" or "ml", these are Fobs, sigFobs (for diffraction)
      !  if target == "vls",  these are Fobs, phiFobs (for cryoEM)

      if (has_Fuser > 0 ) then
         allocate( Fuser(num_hkl) )
         do i = 1,num_hkl
            read(hkl_lun,*,end=1,err=2) &
               hkl_index(1:3,i),abs_Fobs(i),sigFobs(i),test_flag(i), &
               abs_Fuser, phase_Fuser
            Fuser(i) = abs_Fuser*cmplx( cos( DEG_TO_RAD*phase_Fuser) , &
                                        sin( DEG_TO_RAD*phase_Fuser), rk_ )
            test_flag(i) = min(test_flag(i),1)
         end do
      else
         do i = 1,num_hkl
            read(hkl_lun,*,end=1,err=2) &
               hkl_index(1:3,i),abs_Fobs(i),sigFobs(i),test_flag(i)
            test_flag(i) = min(test_flag(i),1)
         end do
      end if
      write(6,'(a,i8,a,i2)') '| found ', num_hkl, &
         ' reflections, with has_Fuser = ', has_Fuser
      if( has_Fuser > 0 ) then
         bulk_solvent_model = 'afonine-2013'
         write(6,'(a)') '|     setting bulk_solvent_model to afonine-2013'
      end if

      if (atom_selection_mask/='') then
         call atommask(natom=natom,nres=nres,prnlev=0, &
               igraph=ih(m04),isymbl=ih(m06),ipres=ix(i02), &
               lbres=ih(m02),crd=x(lcrd), &
               maskstr=atom_selection_mask,mask=atom_selection)
         NAT_for_mask1 = sum(atom_selection)
         write(6,'(a,i6,a,a)') 'Found ',NAT_for_mask1, &
              ' atoms in ', atom_selection_mask
         !  also ignore any atoms with zero occupancy:
         do i=1,natom
            if( atom_occupancy(i) == 0._rk_) atom_selection(i) = 0
         end do
         NAT_for_mask = sum(atom_selection)
         if( NAT_for_mask1 /= NAT_for_mask ) &
            write(6,'(a,i4,a)') 'Removing ',NAT_for_mask1 - NAT_for_mask, &
              ' additional atoms with zero occupancy'
      end if

      ! set up complex Fobs(:), if vector target is requested
      allocate(Fobs(num_hkl),stat=alloc_status)
      REQUIRE(alloc_status==0)
      if( target(1:3) == 'vls' ) then
         do i = 1,num_hkl
            !  sigFobs() here is assumed to be really phi(), in degrees
            phi = sigFobs(i) * DEG_TO_RAD
            Fobs(i) = cmplx( abs_Fobs(i)*cos(phi), abs_Fobs(i)*sin(phi), rk_ )
         end do
      else
         Fobs = abs_Fobs
      endif
      
      call init_interface2( &
         & target, bulk_solvent_model, &
         & hkl_index, Fobs, sigFobs, test_flag==1, &
         & unit_cell, scatter_coefficients, &
         & atom_bfactor, atom_occupancy, atom_scatter_type, &
         & atom_selection==1, ix(i100+1:i100+natom), &
         & mask_update_period, scale_update_period, &
         & ml_update_period, k_sol, b_sol, &
         & solvent_mask_adjustment, solvent_mask_probe_radius, ls_r3, ls_r4, &
         & spacegroup_number, ixp, iyp, izp, timeavg )
      
      ! should be able to do some deallocations here:
      deallocate(hkl_index,Fobs,sigFobs, &
           atom_scatter_type, stat=alloc_status)
      if( alloc_status .ne. 0 ) then
         write(6,*) 'error in deallocation after init_interface2()'
         call mexit(6,1)
      end if

      return
      1 continue
      write(stdout,'(A)') 'End-of-file reading HKL file.'
      call mexit(stdout,1)
      2 continue
      write(stdout,'(A)') 'Error reading HKL file.'
      call mexit(stdout,1)
   end subroutine init

   subroutine xray_init_globals()
      pdb_infile = ''
      pdb_outfile = ''
      fave_outfile = ''
      fmtz_outfile = ''
      sf_outfile = ''
      pdb_read_coordinates = .false.
      pdb_use_segid = .false.
      pdb_wrap_names = .false.
      target = 'ls  '
      spacegroup_name = 'P 1'
      spacegroup_number = 1
      reflection_infile = ''
      xray_weight = 1.0
      solvent_mask_probe_radius = 0.9
      solvent_mask_adjustment = 1.1
      solvent_mask_reflection_outfile = ''
      solvent_mask_outfile = ''
      atom_selection_mask = '!@H='
      mask_update_period = 100
      scale_update_period = 100
      ml_update_period = 100
      bulk_solvent_model = 'none'
      timeavg = .false.
      return
   end subroutine xray_init_globals

   ! Write X-ray output files and deallocate. Bond info is included
   ! here only to check for places to insert TER in PDB output files.
   subroutine finalize()
      use xray_interface2_module, only: finalize2 => finalize
      implicit none
      ! local
      integer :: i
      real(real_kind) :: phi

      if (.not.xray_active) return

      call finalize2()
      
#if 0
      deallocate(atom_bfactor,atom_occupancy,atom_scatter_type, &
            atom_selection,residue_chainid,residue_icode, &
            atom_element,atom_altloc,residue_number, &
            hkl_index,abs_Fobs,sigFobs,test_flag)
#endif

   end subroutine finalize

   ! gets xray_energy and derivatives
   subroutine xray_get_derivative(xyz, force, current_step, xray_e)
      use xray_interface2_module, only: calc_force2 => calc_force, get_r_factors
      implicit none
#include "../include/md.h"
#include "def_time.h"
#include "nmr.h"
      real(real_kind), intent(in) :: xyz(:, :)
      real(real_kind), intent(out) :: force(:, :)
      integer, intent(in) :: current_step
      real(real_kind), intent(out) :: xray_e

      call check_precondition(size(xyz, 1) == 3)
      call check_precondition(size(xyz, 2) == size(force, 2))
      call check_precondition(size(force, 1) == 3)
!      call check_precondition(size(force, 2) == n_atom)

      if (.not. xray_active) then
         xray_e = 0
         xray_energy = 0
         return
      end if

      call timer_start(TIME_XRAY)
      if(nmropt.gt.0) xray_weight = wxray

      call calc_force2(xyz, current_step, xray_weight, force, xray_e, Fuser)
      xray_energy = xray_e
      call get_r_factors(r_work, r_free)
      call timer_stop(TIME_XRAY)

   end subroutine xray_get_derivative

   subroutine xray_write_md_state(unit)
      integer, intent(in) :: unit
      write(unit,'(3(1x,A,f14.4))') &
        'Exray  = ', xray_energy, ' Rwork   = ', r_work, ' Rfree      = ', r_free
   end subroutine xray_write_md_state

   subroutine xray_write_min_state(unit)
      integer, intent(in) :: unit
      write(unit,'(3(1x,A,f13.4))') &
        'Exray   = ', xray_energy, ' Rwork   = ', r_work, ' Rfree      = ', r_free
   end subroutine xray_write_min_state


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Allocate the next unused logical UNIT number from a predefined range.
   ! Intended to be used with amopen, which aborts on error, because there
   ! is no corresponding deallocate function.
   function allocate_lun(lun_return) result(lun)
      integer :: lun
      integer, intent(out), optional :: lun_return
      integer, parameter :: FILE_UNIT_FIRST=201, FILE_UNIT_LAST=250
      integer, save :: unit = FILE_UNIT_FIRST-1
      ! locals
      integer :: i
      logical :: opened
      ! begin
      lun = FILE_UNIT_FIRST
      do i = FILE_UNIT_FIRST, FILE_UNIT_LAST
         unit = unit + 1
         if (unit > FILE_UNIT_LAST) unit = FILE_UNIT_FIRST
         ! This assumes that an allocated unit always results in an opened state.
         inquire(unit=unit,opened=opened)
         if (.not. opened) then
            lun = unit
            if (present(lun_return)) lun_return = unit
            return
         end if
      end do
      write(stdout,'(A)') &
              'ERROR in ALLOCATE_LUN(): ran out of available LUNs!'
      call mexit(stdout,2)
   end function allocate_lun
   
end module xray_cpu_module
