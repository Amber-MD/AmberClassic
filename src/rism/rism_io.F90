!<compile=optimized>

! The 3D-RISM-KH software found here is copyright (c) 2010-2012 by
! Andriy Kovalenko, Tyler Luchko and David A. Case.
!
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
!
! You should have received a copy of the GNU General Public License in the
! ../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
! Users of the 3D-RISM capability found here are requested to acknowledge
! use of the software in reports and publications.  Such acknowledgement
! should include the following citations:
!
! 1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999);
! ibid. 112:10391-10417 (2000).
!
! 2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by
! F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.
!
! 3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
! and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010).

#include "../include/dprec.fh"

!> Minimimal support for text format OpenDX rectangular grid output.
!! There is parallel support for writting but this is done through
!! communicating with the master node so it is very slow.
module rism_io
  use safemem
  use rism_report_c
  implicit none

  !> Abstract interace for functions that write volumetric data
  !! @param[in] file Name of output file
  !! @param[in] data 3D array of _REAL_ data to write
  !! @param[in] grid rism3d_grid object for data
  !! @param[in] solute rism3d_solute object for data
  !! @param[in] o_rank (optional) MPI rank
  !! @param[in] o_nproc (optional) Number of MPI proceses
  !! @param[in] o_comm (optional) MPI communicator
  abstract interface
     subroutine writeVolumeInterface(file, data, &
     periodic, boxLength, localDimsR, unitCellVectorsK, unitCellVectorsR, globalDimsR, &
     unitCellAngles, voxelVectorsR, spacing, totalLocalPointsR, &
     translation, centerOfMass, &
     o_rank, o_nproc, o_comm)
      !  use rism3d_grid_c
      !  use rism3d_solute_c
       character(len=*), intent(in) :: file
      !  type(rism3d_grid), intent(in) :: grid
      !  _REAL_, target, intent(in) :: data(grid%localDimsR(1), grid%localDimsR(2), grid%localDimsR(3))
      !  type(rism3d_solute), intent(in) :: solute

       integer :: globalDimsR(3)
       integer :: localDimsR(3)
       _REAL_, target, intent(in) :: data(localDimsR(1), localDimsR(2), localDimsR(3)) !, centerOfMass(3)
       _REAL_ :: boxLength(3)
       logical :: periodic
       _REAL_ :: voxelVectorsR(3, 3)
       _REAL_ :: unitCellAngles(3)
       _REAL_ :: unitCellVectorsR(3, 3)
       _REAL_ :: unitCellVectorsK(3, 3)
       _REAL_ :: spacing(3)
       integer :: totalLocalPointsR
       _REAL_ :: translation(3)
       _REAL_ :: centerOfMass(3)

       integer, optional :: o_rank, o_nproc, o_comm
     end subroutine writeVolumeInterface

#ifdef CUDA
     subroutine writeVolumeInterface_cpp(this, atomType, file)
       use iso_c_binding, only : C_INT
       use rism3d_rism3d_c_mod
       class(rism3d) :: this
       integer(C_INT), value, intent(IN) :: atomType
       character(len=*), intent(IN) :: file
     end subroutine writeVolumeInterface_cpp
#endif /*CUDA*/

  end interface
  
contains

  !> Read a 1D RDF from a space separated variable file.
  !! The file structure consists of commented lines beginning with
  !! '#', and lines containing individual value pairs of grid position
  !! and its RDF value. Note that grid positions must be consecutive
  !! to simplify grid spacing consistency checking.
  !! @param[in] filename Name of file to be read.
  !! @param[out] rdf 1D radial distribution function.
  !! @param[out] gridSpacing Distance between grid points.
  subroutine readRDF1D(filename, elec_tot, rdf, gridSpacing)
    use rism_util, only : freeUnit
    use safemem
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: elec_tot
    _REAL_, pointer :: rdf(:)
    _REAL_, intent(out) :: gridSpacing

    integer :: unit, iostat
    integer :: id
    character(len=1024) :: buffer
    _REAL_ :: rdfPoint(2)
    _REAL_ :: currentGridPosition, prevGridPosition

    rdf => safemem_realloc(rdf, 400)
    
    unit = freeUnit()
    open(unit=unit, file=filename, status='old', iostat=iostat)
    if (iostat /= 0) then
       call rism_report_error("(a,i4)", "Could not open "//trim(filename)//":", iostat)
    end if

    ! read total number of electrons in the top line
    call nextline(unit, buffer)
    read(buffer, *) elec_tot

    ! silence warning.  Values is defined first time through the loop
    prevGridPosition = 0
    id = 1
    do
       call nextline(unit, buffer)

       ! End of file.
       if (buffer == '') exit
       
       read(buffer, *) rdfPoint

       currentGridPosition = rdfPoint(1)
       rdf(id) = rdfPoint(2)
       
       if (id > 2 .and. &
            ((currentGridPosition - prevGridPosition - gridSpacing) .gt. 1E-10)) then
          print *, "id: ", id
          call rism_report_error("(a,f12.10,a,f12.10,a,f10.2)", &
               "Inconsistent 1D electron RDF grid spacing: ", &
                gridSpacing, " then ", currentGridPosition - prevGridPosition, &
                " at grid position ", currentGridPosition)
       end if
       gridSpacing = currentGridPosition - prevGridPosition
       prevGridPosition = currentGridPosition

       id = id + 1

       ! Enlarge rdf if necessary.
       if (size(rdf) < id) then
          rdf => safemem_realloc(rdf, id * 2, .true.)
       end if
    end do

    ! Shrink rdf to its true size.
    rdf => safemem_realloc(rdf, id - 1, .true.)
    
    close(unit)
  end subroutine readRDF1D
  
  
  !> Read unit cell dimensions from a crd / rst file.  Abort if box info
  !! is not found.
  subroutine readUnitCellDimensionsFromCrd(file, unitCellDimensions)
    use rism_util, only : freeUnit
    implicit none

    character(len=*), intent(in) :: file
    _REAL_, intent(out) :: unitCellDimensions(6)
    integer :: unit, iostat
    character(len=1024) :: buffer

    integer :: numAtoms
    
    unit = freeUnit()
    
    open(unit = unit, file = trim(file), status = 'old', iostat = iostat)
    if (iostat /= 0) then
       call rism_report_error( &
            "(a, i4)", "readUnitCellDimensionsFromCrd: could not open " // trim(file) // ":", iostat)
    end if

    ! ! Skip the molecule name.
    ! call nextline(unit, buffer)

    ! ! Read the atom count.
    ! call nextline(unit, buffer)
    ! read(line, '(i8, e15.7)', iostat = iostat) numAtoms
    ! if (iostat /= 0) then
    !    call rism_report_error("Atom count is missing from second line of crd/rst file.")
    ! else if (numAtoms <= 2) then
    !    call rism_report_error( &
    !         "Crd/rst file must have more than two atoms to deduce unit cell information.")
    ! end if

    ! The unit cell dimensions should be the last non-empty line in
    ! the file.
    call fseek(unit, 0, 2)
    !TODO: Handle case of empty last lines prior to end of file.
    !TODO: Handle case where unit cell dimensions are not provided.
    ! buffer = ""
    ! do while (trim(buffer(1:1)) == "")
    backspace(unit, iostat = iostat)
    if (iostat /= 0) then
       call rism_report_error( &
            "readUnitCellDimensionsFromCrd: could seek last line of crd/rst file.")
    end if
    read(unit, '(a)', iostat = iostat) buffer
    if (iostat /= 0) then
       call rism_report_error( &
            "readUnitCellDimensionsFromCrd: could not read last line of crd/rst file.")
    end if
    ! end do
    

    ! Convert dimensions string to floats.

    !FIXME: This conversion creates errors in the last floating point
    ! decimal place, which is a greater precision than the string it
    ! is obtained from.
    ! read(buffer, '(6d12.7)', iostat = iostat) unitCellDimensions
    read(buffer, *, iostat = iostat) unitCellDimensions
    if (iostat /= 0) then
       call rism_report_error( &
            "readUnitCellDimensionsFromCrd: last line of crd/rst file does not " &
            // "appear to be unit cell dimensions.")
    end if

    !TODO: Handle multiple unit cell dimensions formats.
    !     if (ic == justcrd + 1 .or. ic == vel + 1) then
    !        write(6,'(a)') '| peek_ewald_inpcrd: Box info found'
    !        a = x1
    !        b = x2
    !        c = x3
    !        if (x4 > 0.d0 .and. x5 == 0.d0 .and. x6 == 0.d0) then
    !           ! Only has beta.
    !           alpha = 90.d0
    !           beta = x4
    !           gamma = 90.d0
    !        else if (x4 == 0.d0 .and. x5 == 0.d0 .and. x6 == 0.d0) then
    !           ! No angles in input: assume they are all 90:
    !           alpha = 90.d0
    !           beta  = 90.d0
    !           gamma = 90.d0
    !        else
    !           ! Found the angles.
    !           alpha = x4
    !           beta = x5
    !           gamma = x6
    !        end if
    !        return
    !     end if
    
    close(unit)

  end subroutine readUnitCellDimensionsFromCrd

  
  !> Reads in the next non-comment line from the file.  A comment is
  !! defined as starting with a '#'.
  !! @param[in] unit Open unit.
  !! @param[out] buffer A character buffer to read into.
  subroutine nextline(unit, buffer)
    implicit none
    integer, intent(in) :: unit
    character(len=*), intent(out) :: buffer
    integer :: iostat
    do
       read(unit,'(a)',iostat=iostat) buffer
       if (iostat /= 0) then
          buffer=""
          exit
       end if
       if (buffer(1:1) /= "#") exit
    end do
  end subroutine nextline
  
end module rism_io
