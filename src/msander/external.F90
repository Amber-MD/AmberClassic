!*******************************************************************************
!
! Module: external_module
!
! Description: This module houses the capabilities of calling external libraries
!              to compute energies and forces, instead of using the force field ones.
!
!              Adapted here by Vinicius Wilian D. Cruzeiro and Delaram Ghoreishi
!
!*******************************************************************************

#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

module external_module

  use UtilitiesModule, only: Upcase
#ifdef QUICK
  use quick_api_module, only : setQuickJob, getQuickEnergyGradients, deleteQuickJob
#endif
#ifdef TORCHANI
  use torchani, only: &
      torchani_energy_force, &
      torchani_energy_force_pbc, &
      torchani_energy_force_from_external_neighbors, &
      convert_sander_neighborlist_to_ani_fmt, &
      torchani_init_model
#endif
  implicit none

  private

  character(256), public :: extprog, json, keywords, outfprefix
  character(80), public :: method, basis
  integer, public :: charge, mult

  ! The active namelist:

  private       :: extpot

  namelist /extpot/      extprog, json, outfprefix, method, basis, charge, mult, keywords

  public :: &
    external_init, &
    gb_external, &
    pme_external, &
    pme_external_from_neighborlist, &
    pme_external_needs_neighborlist, &
    pme_external_dummyiblo, &
    external_cleanup

  ! Read by force, if .true.,
  ! it dispatches "pme_external_from_neighborlist" instead of "pme_external".
  ! Currently used by Torchani-Amber
  logical, protected :: pme_external_needs_neighborlist = .false.  ! Read-only
  ! Dummymask is only init if pme_external_needs_neighborlist == .true.
  ! It is passed to pme_list to generate a pairlist that includes all atom pairs
  integer, allocatable, protected :: pme_external_dummyiblo(:)  ! Read-only

#ifdef TORCHANI
  ! NOTE: TorchANI does not rely on the "extpot" namelist, it uses "ani" instead
  ! "ani" namelist vars
  character(256) :: model_type
  logical :: use_cuda_device
  logical :: use_double_precision
  logical :: use_amber_neighborlist
  logical :: use_cuaev
  ! "ani" namelist vars (advanced and dev-only)
  integer :: model_index
  integer :: cuda_device_index
  namelist /ani/&
      use_cuda_device,&
      use_double_precision,&
      use_amber_neighborlist,&
      use_cuaev,&
      cuda_device_index,&
      model_index,&
      model_type
#endif  /* TORCHANI */

  contains

  subroutine external_init(ix,ih,xx)

    use memory_module,   only : natom, m04, nres, m02, i02, i100, lmass
    use qmmm_module,     only : get_atomic_number

    IMPLICIT NONE

    integer :: i, ifind, cnt, ierr
    _REAL_ :: coord(3*natom)
    character(len=5), allocatable :: monomers(:)
    integer, allocatable :: nats(:)
    character(len=5), allocatable :: at_name(:)
    integer :: nmon, val
    integer :: ix(*)
    character(len=4) :: ih(*)
    _REAL_ :: xx(*)
    integer :: atomicnumber(natom)
    logical :: isvalid, errFlag
#ifdef TORCHANI
    integer :: iostat
    iostat = 0
#endif
#ifdef QUICK
    character(256) :: keys
#endif

!    write (6,*) " -> Printing atom names"
!    do i=1,natom
!      write (6,*) "  * ", ih(m04+i-1)
!    end do
!    write (6,*) " -> Printing residue names"
!    do i=1,nres
!      write (6,*) "  * ", ih(m02+i-1)
!    end do
!    write (6,*) " -> Printing residue pointers"
!    do i=1,nres
!      write (6,*) "  * ", ix(i02+i-1)
!    end do

    extprog = ''
    json = ''
    outfprefix = 'quick'
    method = ''
    basis = ''
    charge = 0
    mult = 1
    keywords = ''

    ! Read input in namelist format:

    rewind(5)                  ! Insurance against maintenance mods.

    call nmlsrc('extpot', 5, ifind)

    if (ifind .ne. 0) then        ! Namelist found. Read it:
      read(5, nml = extpot)
    else                          ! No namelist was present,
      write(6, '(a)') 'ERROR: Could not find the namelist extpot in the mdin file'
      call mexit(6, 1)
    end if

    if (upcase(extprog) .eq. 'KMMD') then
        call kmmd_init(trim(json)//CHAR(0))
        return
    end if
    
    
    do i=1, natom
      if(ix(i100) .eq. 1) then
        atomicnumber(i) = ix(i100+i)
      else
        call get_atomic_number(ih(m04+i-1), xx(lmass+i-1), atomicnumber(i), errFlag)
      end if
    end do

    ! For MBX
    if (Upcase(extprog) .eq. 'MBX') then
#ifdef MBX
      nmon = nres

      allocate(nats(nmon),monomers(nmon),at_name(natom))

      cnt = 1
      do i=1,nmon
         isvalid = .True.
         if (ih(m02+i-1) == "WAT") then
            val = 3
            monomers(i) = "h2o"
            at_name(cnt+0) = "O"
            if (atomicnumber(cnt+0) .ne. 8) isvalid = .False.
            at_name(cnt+1) = "H"
            if (atomicnumber(cnt+1) .ne. 1) isvalid = .False.
            at_name(cnt+2) = "H"
            if (atomicnumber(cnt+2) .ne. 1) isvalid = .False.
            cnt = cnt + val
         else if (ih(m02+i-1) == "N2O") then
            val = 7
            monomers(i) = "n2o5"
            at_name(cnt+0) = "O"
            if (atomicnumber(cnt+0) .ne. 8) isvalid = .False.
            at_name(cnt+1) = "N"
            if (atomicnumber(cnt+1) .ne. 7) isvalid = .False.
            at_name(cnt+2) = "N"
            if (atomicnumber(cnt+2) .ne. 7) isvalid = .False.
            at_name(cnt+3) = "O"
            if (atomicnumber(cnt+3) .ne. 8) isvalid = .False.
            at_name(cnt+4) = "O"
            if (atomicnumber(cnt+4) .ne. 8) isvalid = .False.
            at_name(cnt+5) = "O"
            if (atomicnumber(cnt+5) .ne. 8) isvalid = .False.
            at_name(cnt+6) = "O"
            if (atomicnumber(cnt+6) .ne. 8) isvalid = .False.
            cnt = cnt + val
         else
            write(6, '(a,a,a)') 'ERROR: The residue ',ih(m02+i-1),' is not recognized by MBX!'
            call mexit(6, 1)
         end if
         if (val == ix(i02+i)-ix(i02+i-1)) then
            nats(i) = val
         else
            write(6, '(a,a,a)') 'ERROR: The number of atoms in residue ',ih(m02+i-1),' does not match the expected by MBX!'
            call mexit(6, 1)
         end if
         if (.not. isvalid) then
            write(6, '(a,a,a)') 'ERROR: The order or type of the atoms in residue ',ih(m02+i-1),&
                                    ' does not match the expected by MBX!'
            call mexit(6, 1)
         end if
      end do

      do i =1, natom
        at_name(i)=trim(at_name(i))//CHAR(0)
      end do
      do i=1,nmon
        monomers(i) = trim(monomers(i))//CHAR(0)
      end do

      if (json /= '') then
        call initialize_system(coord, nats, at_name, monomers, nmon, trim(json)//CHAR(0))
      else
        call initialize_system(coord, nats, at_name, monomers, nmon)
      end if
#endif
#ifdef TORCHANI
    elseif (upcase(extprog) == 'TORCHANI') then
        ! Init the "ani" namelist vars
        model_type = 'ani1x'
        use_cuaev = .false.
        use_cuda_device = .true.
        use_double_precision = .true.
        use_amber_neighborlist = .false.
        model_index = -1
        cuda_device_index = 0

        ! Find the "ani" namelist in mdin, and read it to override defaults
        rewind(5)
        read(5, nml=ani, iostat=iostat)
        if (iostat /= 0) then
            call sander_bomb(&
                'external_init (external)',&
                '"ani" namelist read error',&
                'Error reading "ani" namelist. Please check your mdin input'&
            )
        end if

        if (use_amber_neighborlist) then
            pme_external_needs_neighborlist = .true.
            allocate(pme_external_dummyiblo(natom))
            pme_external_dummyiblo = 0
        endif

        call torchani_init_model( &
            atomicnumber, &
            cuda_device_index, &
            model_type, &
            model_index, &
            use_double_precision, &
            use_cuda_device, &
            use_cuaev &
        )
#endif  /* TORCHANI */
    ! For Quick
    else if (Upcase(extprog) .eq. 'QUICK') then
#ifdef QUICK
      ! Check if input flags are ok
      if (keywords .eq. '' .and. method .eq. '') then
        write(6, '(a,a)') 'ERROR: Please specify the keywords flag or method and basis',&
                              ' flags for QUICK in the extpot namelist!'
        call mexit(6, 1)
      end if
      if (keywords .eq. '' .and. method .ne. '' .and. basis .eq. '') then
        write(6, '(a,a)') 'ERROR: Please specify the basis set for QUICK using the',&
                              ' basis flag in the extpot namelist!'
        call mexit(6, 1)
      end if

      ! Constructing the keywords flag, if needed
      if (keywords .eq. '') then
        if (Upcase(method) .eq. 'HF') then
          write(keys, '(4a,2(i0,a))') trim(method), ' BASIS=', trim(basis),&
                                      ' CHARGE=', charge, ' MULT=', mult,' GRADIENT'
        else
          write(keys, '(5a,2(i0,a))') 'DFT ', trim(method), ' BASIS=', trim(basis),&
                                      ' CHARGE=', charge, ' MULT=', mult,' GRADIENT'
        end if
      else
        write(keys, '(a)') trim(keywords)
      end if
      call setQuickJob(outfprefix, keys, natom, atomicnumber, .true., ierr)
      if ( ierr /= 0 ) then
        write(6, '(a)') 'ERROR: setting up Quick in setQuickJob at external.F90'
        call mexit(6, 1)
      end if
#endif
    else
      write(6, '(a,a,a)') 'ERROR: External program ',trim(extprog),&
             ' is not valid! Please set a valid value in the extprog flag'
      call mexit(6, 1)
    end if

  end subroutine

  subroutine gb_external(crd, frc, pot_ene)

    use memory_module,   only : natom
    use constants, only :  EV_TO_KCAL, AU_TO_EV, BOHRS_TO_A

    IMPLICIT NONE

    _REAL_   ::  crd(*)
    _REAL_   ::  frc(*)
    _REAL_   ::  pot_ene
    integer  ::  i, i3
    _REAL_   ::  coord(3*natom), grads(3*natom)
#ifdef QUICK
    integer  ::  ierr
    _REAL_, dimension(3,natom) :: crd_quick, frc_quick
    _REAL_, dimension(:,:), pointer :: xc_coord       => null()
    _REAL_, dimension(:,:), pointer :: ptchgGrad      => null()
#endif

    ! For MBX
    if (Upcase(extprog) .eq. 'MBX') then
#ifdef MBX
      do i = 1, natom
        i3=3*(i-1)
        coord(3*(i-1)+1) = crd(i3+1)
        coord(3*(i-1)+2) = crd(i3+2)
        coord(3*(i-1)+3) = crd(i3+3)
      end do

      call get_energy_g(coord, natom, pot_ene, grads)

      do i = 1, natom
        i3=3*(i-1)
        frc(i3+1) = frc(i3+1) - grads(3*(i-1)+1)
        frc(i3+2) = frc(i3+2) - grads(3*(i-1)+2)
        frc(i3+3) = frc(i3+3) - grads(3*(i-1)+3)
      end do
#endif
#ifdef TORCHANI
    elseif (upcase(extprog) .eq. "TORCHANI") then
        ! Delegate non-PBC force calculation to TorchANI-Amber
        ! NOTE: Unclear whether zeroing the frc array and pot_ene is needed or not
        ! In principle Sander routines should zero frc before calling TorchANI-Amber,
        ! If that is so, this is a no-op
        pot_ene = 0.0d0
        do i = 1, natom
          frc(3 * (i - 1) + 1) = 0.0d0
          frc(3 * (i - 1) + 2) = 0.0d0
          frc(3 * (i - 1) + 3) = 0.0d0
        end do
        ! The crd and frc arrays underlying buffers already have the correct memory layout
        ! for TorchANI-Amber, namely [x0 y0 z0 x1 y1 z1 x2 y2 z2...]
        ! The numel is 3 * natom for both
        call torchani_energy_force(&
            natom, &
            crd, &
            frc, &
            pot_ene &
        )
#endif /* TORCHANI */
    ! For Quick
    else if (Upcase(extprog) .eq. 'QUICK') then
#ifdef QUICK
      do i = 1, natom
        crd_quick(1,i) = crd(3*(i-1)+1)
        crd_quick(2,i) = crd(3*(i-1)+2)
        crd_quick(3,i) = crd(3*(i-1)+3)
        frc_quick(1,i) = 0.0
        frc_quick(2,i) = 0.0
        frc_quick(3,i) = 0.0
      end do
      call getQuickEnergyGradients(crd_quick, 0, xc_coord, pot_ene, frc_quick, ptchgGrad, ierr)
      if ( ierr /= 0 ) then
        write(6, '(a)') 'ERROR: getting energy and gradient with Quick in getQuickEnergyGradients at external.F90'
        call mexit(6, 1)
      end if
      do i = 1, natom
        frc(3*(i-1)+1) = - frc_quick(1,i) * AU_TO_EV * EV_TO_KCAL / BOHRS_TO_A
        frc(3*(i-1)+2) = - frc_quick(2,i) * AU_TO_EV * EV_TO_KCAL / BOHRS_TO_A
        frc(3*(i-1)+3) = - frc_quick(3,i) * AU_TO_EV * EV_TO_KCAL / BOHRS_TO_A
      end do
      pot_ene = pot_ene * AU_TO_EV * EV_TO_KCAL
#endif
    else if (upcase(extprog) .eq. 'KMMD') then
      do i = 1, natom
        coord(3*(i-1)+1) = crd(3*(i-1)+1)
        coord(3*(i-1)+2) = crd(3*(i-1)+2)
        coord(3*(i-1)+3) = crd(3*(i-1)+3)
        grads(3*(i-1)+1) = 0.0
        grads(3*(i-1)+2) = 0.0
        grads(3*(i-1)+3) = 0.0
      end do

      call kmmd_frccalc(coord, grads, pot_ene)

      do i = 1, natom
        frc(3*(i-1)+1) = frc(3*(i-1)+1) + grads(3*(i-1)+1)
        frc(3*(i-1)+2) = frc(3*(i-1)+2) + grads(3*(i-1)+2)
        frc(3*(i-1)+3) = frc(3*(i-1)+3) + grads(3*(i-1)+3)
      end do
    end if

  end subroutine

  subroutine pme_external(crd, frc, pot_ene)

    use memory_module, only : natom
    use nblist, only: ucell

    IMPLICIT NONE

    double precision           ::  crd(*)
    double precision           ::  frc(*)
    double precision           ::  pot_ene
    integer                    ::  i, i3
    double precision           ::  coord(3*natom), grads(3*natom), box(9)

#ifdef TORCHANI
    if (upcase(extprog) .eq. 'TORCHANI') then
        ! Delegate PBC force calculation to TorchANI-Amber
        ! NOTE: Unclear whether zeroing the frc array and pot_ene is needed or not
        ! In principle Sander routines should zero frc before calling TorchANI-Amber,
        ! If that is so, this is a no-op
        pot_ene = 0.0d0
        do i = 1, natom
          frc(3 * (i - 1) + 1) = 0.0d0
          frc(3 * (i - 1) + 2) = 0.0d0
          frc(3 * (i - 1) + 3) = 0.0d0
        end do
        ! The crd and frc arrays underlying buffers already have the correct memory layout
        ! for TorchANI-Amber, namely [x0 y0 z0 x1 y1 z1 x2 y2 z2...]
        ! The numel is 3 * natom for both
        call torchani_energy_force_pbc( &
            natom, &
            crd, &
            ucell, &
            frc, &
            pot_ene &
        )
    endif
#endif  /* TORCHANI */

    ! For MBX
    if (Upcase(extprog) .eq. 'MBX') then
#ifdef MBX
      do i = 1, natom
        i3=3*(i-1)
        coord(3*(i-1)+1) = crd(i3+1)
        coord(3*(i-1)+2) = crd(i3+2)
        coord(3*(i-1)+3) = crd(i3+3)
      end do

      box(1) = ucell(1,1)
      box(2) = ucell(2,1)
      box(3) = ucell(3,1)
      box(4) = ucell(1,2)
      box(5) = ucell(2,2)
      box(6) = ucell(3,2)
      box(7) = ucell(1,3)
      box(8) = ucell(2,3)
      box(9) = ucell(3,3)

      call get_energy_pbc_g(coord, natom, box, pot_ene, grads)

      do i = 1, natom
        i3=3*(i-1)
        frc(i3+1) = frc(i3+1) - grads(3*(i-1)+1)
        frc(i3+2) = frc(i3+2) - grads(3*(i-1)+2)
        frc(i3+3) = frc(i3+3) - grads(3*(i-1)+3)
      end do
#endif
    end if
    
    if (upcase(extprog) .eq. 'KMMD') then
      do i = 1, natom
        i3=3*(i-1)
        coord(3*(i-1)+1) = crd(i3+i)
        coord(3*(i-1)+2) = crd(i3+i+1)
        coord(3*(i-1)+3) = crd(i3+i+2)
        grads(3*(i-1)+1) = 0.0
        grads(3*(i-1)+2) = 0.0
        grads(3*(i-1)+3) = 0.0
      end do

      call kmmd_frccalc(coord, grads, pot_ene)

      do i = 1, natom
        i3=3*(i-1)
        frc(i3+1) = frc(i3+1) + grads(3*(i-1)+1)
        frc(i3+2) = frc(i3+2) + grads(3*(i-1)+2)
        frc(i3+3) = frc(i3+3) + grads(3*(i-1)+3)
      end do
    end if

  end subroutine

  subroutine external_cleanup()

    IMPLICIT NONE

    integer :: ierr

    ! For Quick
    if (Upcase(extprog) .eq. 'QUICK') then
#ifdef QUICK
       call deleteQuickJob(ierr)
       if ( ierr /= 0 ) then
         write(6, '(a)') 'ERROR: ending Quick in deleteQuickJob at external.F90'
         call mexit(6, 1)
       end if
#endif
    end if

  end subroutine

! Calculate the energy and force from an external potential, using sander's neighborlist
! At the point this function is called, the sander_neighborlist needs to be built
! using a dummy mask that includes all pairs of atoms
! NOTE: The coordinates passed to this function are image-ordered, and thus don't match
! with the atomic numbers passed to external_init, re-ordering is necessary if this
! is important for the external potential
subroutine pme_external_from_neighborlist( &
    sander_neighborlist, &
    crd, &
    frc, &
    ene_pot &
)
    ! numvdw and numbnd take atom-idx, not img-idx
    ! backptr takes img-idx and converts it to an atom-idx
    use nblist, only: tranvec, bckptr, numvdw, numhbnd, imagcrds, recip, &
        adjust_imagcrds, map_coords
    use memory_module, only : natom
    integer, intent(in) :: sander_neighborlist(*)
    double precision, intent(in) :: crd(*)
    double precision, intent(out) :: frc(*)
    double precision, intent(out) :: ene_pot
#ifdef TORCHANI
    integer :: j
    double precision, allocatable :: atom_ordered_coords(:,:)
    integer, allocatable :: ani_neighborlist(:, :)
    double precision, allocatable :: ani_shifts(:, :)
#endif

    call map_coords(crd, natom, recip)
    call adjust_imagcrds(crd, natom)

#ifdef TORCHANI
    allocate(atom_ordered_coords(3, natom))

    ! Both the "atom_ordered_coords" and the "image coords" are mapped to the central cell
    ! But the "atom_ordered_coords" have the same order os the "atomic nums" passed to external_init
    ! This could be re-written in such a way that only the atomic nums are
    ! flipped once and cached, but 1 copy of the coords is not a bottleneck for ANI,
    ! so that optimization is not really worth it
    do j = 1, natom
        atom_ordered_coords(:, bckptr(j)) = imagcrds(:, j)
    end do

    if (sum(numhbnd) /= 0) then
        call sander_bomb(&
            'pme_external_from_neighborlist (external)',&
            'torchani does not support 10-12 H-bond interactions',&
            'Error, trying to run torchani with a 10-12 H-bond force field'&
        )
    endif

    if (upcase(extprog) == "TORCHANI") then
        ! sander_neighborlist holds 1-idxing img_idxs (into "img coords"), it
        ! has shape (max-num-pairs,). Only the first (num-pairs,) are useful.
        ! The order is relevant in this list, and corresponds to the img-order

        ! The obtained ani_neighborlist holds 0-idxing atom_idxs.
        ! ani_neighborlist and ani_shifts are fshape (num-pairs, 2) and (3, num-pairs). The
        ! order of the pairs is not relevant, but it is consistent for both
        call convert_sander_neighborlist_to_ani_fmt( &
            bckptr, &
            numvdw, &
            tranvec, &
            sander_neighborlist, &
            ani_neighborlist, &
            ani_shifts &
        )

        call torchani_energy_force_from_external_neighbors( &
            natom, &
            size(ani_neighborlist, dim=1), &
            atom_ordered_coords, &
            ani_neighborlist, &
            ani_shifts, &
            frc, &
            ene_pot &
        )
    end if
#endif
end subroutine
end module external_module
