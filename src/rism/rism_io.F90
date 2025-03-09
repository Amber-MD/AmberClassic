!<compile=optimized>

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
     subroutine writeVolumeInterface(file, data, grid, solute, o_rank, o_nproc, o_comm)
       use rism3d_grid_c
       use rism3d_solute_c
       character(len=*), intent(in) :: file
       type(rism3d_grid), intent(in) :: grid
       _REAL_, target, intent(in) :: data(grid%localDimsR(1), grid%localDimsR(2), grid%localDimsR(3))
       type(rism3d_solute), intent(in) :: solute
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
