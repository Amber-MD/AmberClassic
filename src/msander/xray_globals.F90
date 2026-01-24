! <compile=optimized>
#include "../include/dprec.fh"
module xray_globals_module
   use file_io_dat, only : MAX_FN_LEN
   use xray_unit_cell_module, only: unit_cell_t
   implicit none
   public

   _REAL_, private :: r_dummy
   integer, parameter :: real_kind = kind(r_dummy)
   integer, parameter :: rk_ = real_kind

   ! HP Fortran uses unit 7 for stderr. Most others use 0.
   ! F2003 gives the actual values in intrinsic module ISO_FORTRAN_ENV.
   ! SANDER does not use stderr, so it is assigned to stdout here,
   ! which is really mdout.
   integer, parameter :: STDERR=6, STDIN=5, STDOUT=6

   logical, save :: xray_active = .false.

   !-------------------------------------------------------------------
   ! NameList Input Parameters
   integer, parameter :: REFL_LABEL_MAXLEN=32

   character(len=MAX_FN_LEN), save :: pdb_infile, pdb_outfile, sf_outfile, &
                                      fave_outfile, fmtz_outfile
   integer, save :: n_fcalc_ave = 0
   integer, save :: ntwsf = 0
   
   ! If true, PDB coordinates will overwrite INPCRD coordinates.
   ! NOTE: the cell still comes from the INPCRD!
   logical, save :: pdb_read_coordinates

   ! If TRUE, write the full 4-character ChainID to the SegID field.
   logical, save :: pdb_use_segid, pdb_wrap_names

   ! use time-averaging?
   logical, save :: timeavg

   ! Xray energy target function, see xray_target.F90 for details
   character(len=4), save :: target
   ! variable parameters for the least-squares target:
   real(real_kind), save :: ls_r3 = 0.0, ls_r4 = 0.5

   ! Standard condensed spacegroup name, or integer spacegroup number
   character(len=16), save :: spacegroup_name

   ! Filename for reflection input file.
   character(len=MAX_FN_LEN), save :: reflection_infile

   ! Weight term for X-ray force:
   real(real_kind), save :: xray_weight

   !> Increment to be added to atomic radii of the atoms selected
   !  by atom_selection_mask as a part of the algorithm to build bulk mask
   real(real_kind), save :: solvent_mask_adjustment

   !> The radius of solvent probe to apply as a part of the algorithm
   !  to build bulk solvent mask (shrinks non-bulk volume)
   real(real_kind), save :: solvent_mask_probe_radius

   ! Output file for bulk-solvent reflections (Fbulk) and mask
   character(len=MAX_FN_LEN), save :: solvent_mask_reflection_outfile
   character(len=MAX_FN_LEN), save :: solvent_mask_outfile

   character(len=256) :: atom_selection_mask

   real(real_kind), save :: ihkl_duration=0._rk_, dhkl_duration=0._rk_

   !----------------------------------------------------------------------------
   ! GLOBALS:

   !----------------------------------------------------------------------------
   ! Atom data:
   integer, save :: num_atoms, num_residues, NAT_for_mask
   real(real_kind), allocatable, save :: atom_bfactor(:), atom_occupancy(:)
   integer, allocatable, save :: atom_scatter_type(:)
   integer, allocatable, save :: atom_selection(:)

   ! Residue and atom data that may become SANDER globals:
   character(len=1), allocatable, save :: residue_chainid(:), &
                                          residue_icode(:), atom_altloc(:)
   character(len=4), allocatable, save :: atom_element(:)
   integer, allocatable, save :: residue_number(:)

   !----------------------------------------------------------------------------
   ! Reflection data:

   integer, save :: num_hkl
   integer, save :: num_free_flags
   integer, save :: num_work_flags
   
   integer, allocatable, target, save :: hkl_index(:,:) ! (3,num_hkl)

   real(real_kind), allocatable, save :: abs_Fobs(:), sigFobs(:)
   complex(real_kind), allocatable, save :: Fuser(:)
   integer, allocatable, save :: test_flag(:)  ! 0 -- "free set" ; 1 -- "work set"
   integer, save :: scale_update_period, &
          ml_update_period, mask_update_period

   !----------------------------------------------------------------------------
   ! Symmetry and transformations:
   type(unit_cell_t) :: unit_cell
   integer, parameter :: MAX_SYMMOPS = 16 ! actually 96
   real(real_kind), save :: symmop(3,4,MAX_SYMMOPS), symmop_inv(3,4,MAX_SYMMOPS)
   integer, save :: num_symmops

   integer, save :: spacegroup_number
   integer, save :: au_type ! Laue code index
   integer, save :: system  ! i.e. SYMM_TRICLINIC

   !  number of unit-cells in supercell (x,y,z):
   integer, save :: ixp=1, iyp=1, izp=1

   real(real_kind), save :: xray_energy, r_work, r_free

   ! bulk solvent model
   character(len=16), save :: bulk_solvent_model

end module xray_globals_module
