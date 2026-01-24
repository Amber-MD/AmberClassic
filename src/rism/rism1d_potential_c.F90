! <compile=optimized>

! The 3D-RISM-KH software found here is copyright (c) 2011-2012 by
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

!> Potential class for 1D-RISM.  Used to calculate/store quantities that are
!! potential dependent and do not change during the calculation.
module rism1d_potential_c
  use safemem
  implicit none
  integer, private, parameter :: charlen = 8

 type rism1d_potential
     !> Type of RISM theory.  Current options are 'rism' (XRISM) and
     !! 'drism' (dielectrically conistent RISM).
     character(len = charlen) :: theory = ''

     !> Controls the precision in key parts of the algorithm.
     !! 0 - no extra precision.
     !! 1 - extra precision for wzvv and zkvv.
     integer :: extra_precision

     !> r-space grid spacing.
     _REAL_ :: dr = 0
     !> k-space grid spacing.
     _REAL_ :: dk = 0

     !> 1D grid size.
     integer :: nr = 0

     !> (input) Solvent temperature.
     _REAL_ :: temperature
     !> (input) Dielectric constant.
     _REAL_ :: dielconst
     !> (input) Charge smearing parameter for long range asymptotics.
     _REAL_ :: smear
     !> (input) Parameter use in calculating Zvv(k).
     _REAL_ :: adbcor

     !> densitydm * cdiel for rism_dT.
     _REAL_ :: pcdiel

     !> Number of molecular species.
     integer :: nsp = 0
     !> Number of sites (atoms).
     integer :: nv = 0
     !> Number of site-site interactions (nv * (nv + 1) / 2).
     integer :: nvv = 0

     !> Name of the species.
     character(len=80), pointer :: namesp(:) => NULL()

     !> Atom site names.
     character(len=4), pointer :: namev(:) => NULL()

     !> Number of atom types for each species.
     integer, pointer :: nat(:) => NULL()
     !> Multiplicity for each atom in each type.
     integer, pointer :: mta(:,:) => NULL()
     !> Multiplicity for each atom by site.
     integer, pointer :: mtv(:) => NULL()
     !> Triangle packing pointer.
     integer, pointer :: jvv(:,:) => NULL()
     !> Modified non-bond pair.
     integer, pointer :: modnbpair(:) => NULL()
     integer, pointer :: modnbpairtsf(:) => NULL()

     !> (output) Inverse Debye length. [1/A]
     _REAL_ :: kappa = HUGE(1d0)

     !> Mass of each atom type. [g/mol]
     _REAL_, pointer :: mass(:) => NULL()
     !> Number density of each molecular species. [A^{-3}]
     _REAL_, pointer :: densitysp(:) => NULL()
     !> Number density of each atomic site. [A^{-3}]
     _REAL_, pointer :: densityv(:) => NULL()
     !> Charge of each molecular species. [sqrt(kT A)]
     _REAL_, pointer :: qsp(:) => NULL()
     !> Charge of each atomic site in each molecular species. [sqrt(kT A)]
     _REAL_, pointer :: qat(:,:) => NULL()
     !> Charge of each atomic site. [sqrt(kT A)]
     _REAL_, pointer :: qv(:) => NULL()
     !> Charge of the parent molecular species for each site. [sqrt(kT A)]
     _REAL_, pointer :: qspv(:) => NULL()
     !> LJ epsilon for each site in. [kT]
     _REAL_, pointer :: ljEpsilonVV(:) => NULL()
     !> LJ sigma* (rmin) for each site in. [A]
     _REAL_, pointer :: ljRminVV(:) => NULL()
     !> Positions of all atoms in each species.
     _REAL_, pointer :: rma(:,:,:,:) => NULL()
     !> wlmvv
     _REAL_, pointer :: wlmvv(:,:,:) => NULL()
     !> Omega_VV. Intramolecular distance matrix.
     _REAL_, pointer :: wvv(:,:,:) => NULL()
     !> wvv + densityv * zkvv/k
     _REAL_, pointer :: wzvv(:,:,:) => NULL()
     !> The dielectric correction for DRISM in reciprocal space
     !> (commonly written as \chi in the literature) multiplied by k.
     _REAL_, pointer :: zkvv(:,:,:) => NULL()
     !> uvv
     _REAL_, pointer :: uvv(:,:) => NULL()
     !> r-space long range smeared Coulomb site-site potential.
     _REAL_, pointer :: ulrvv(:,:) => NULL()
     !> k-space long range smeared Coulombsite-site potential.
     _REAL_, pointer :: ulkvv(:,:) => NULL()
     !> r-space long range smeared potential contribution to (asymptotics) of hvv.
     _REAL_, pointer :: hlrvv(:,:) => NULL()
     !> k-space long range smeared potential contribution to (asymptotics) of hvv.
     _REAL_, pointer :: hlkvv(:,:) => NULL()
     !> r-space derivative of the potential.
     _REAL_, pointer :: duvv(:,:) => NULL()
     !> Modified Lennard-Jones parameters.
     _REAL_, pointer :: modnb(:,:) => NULL()
     _REAL_, pointer :: modnbtsf(:,:) => NULL()
  end type rism1d_potential

  private intramolecular_dist, sanity_check

contains

  !> Initializes new potential
  !! @param[in,out]  this Potential object.
  !! @param[in] theory 'XRISM' or 'DRISM'.
  !! @param[in] temperature Temperature in [K].
  !! @param[in] dielconst Dielectric constant of the solution.
  !! @param[in] smear Charge smear parameter.
  !! @param[in] adbcor DRISM parameter.
  !! @param[in] nr Number of grid points.
  !! @param[in] dr r-space grid spacing.
  !! @param[in] extra_precision If > 0, use extra precision in key locations.
  subroutine rism1d_potential_new(this, theory, temperature, &
       dielconst, smear, adbcor, nr, dr, extra_precision)
    use constants, only : pi
    use rism_util, only : caseup
    implicit none
    type(rism1d_potential), intent(inout) ::this
    character(len=*), intent(in) :: theory
    integer, intent(in) :: nr, extra_precision
    _REAL_, intent(in) :: temperature, dielconst, smear, adbcor, dr
    call caseup(this%theory)
    this%theory = trim(theory)
    this%temperature = temperature
    this%dielconst = dielconst
    this%smear = smear
    this%adbcor = adbcor
    this%nr = nr
    this%dr = dr
    this%dk = pi / (nr * dr)
    this%extra_precision = extra_precision

    ! Label kappa as 'undefined'.
    this%kappa = HUGE(1d0)

    call sanity_check(this)
  end subroutine rism1d_potential_new


  !> Add solvent species.  Allocates necessary memory to accomodate
  !! the additional species.
  !! @param[in,out] this potential object.
  !! @param[in] mdl solvMDL object.
  !! @param[in] density Number density [A^{-3}].
  subroutine rism1d_potential_addSpecies(this, mdl, density)
    use solvMDL_c
    implicit none
    type(rism1d_potential), intent(inout) :: this
    type(solvMDL), intent(in) :: mdl
    _REAL_, intent(in) :: density
    integer :: oldnv
    integer :: iv1, iv2, ivv, iv, iat, imlt, isp

    ! Increment scalar variables.
    oldnv = this%nv
    this%nsp = this%nsp + 1
    this%nv = this%nv + mdl%ntype
    this%nvv = this%nv * (this%nv + 1) / 2

    ! (Re)allocate arrays and fill with data.
    this%namesp => safemem_realloc(this%namesp, len(this%namesp), this%nsp)
    this%namesp(this%nsp) = mdl%title
    this%namev => safemem_realloc(this%namev, len(this%namev), this%nv)
    this%namev(oldnv + 1:this%nv) = mdl%atmname
    this%nat => safemem_realloc(this%nat, this%nsp)
    this%nat(this%nsp) = mdl%ntype
    this%mta => safemem_realloc(this%mta, max(ubound(this%mta, 1), mdl%ntype), this%nsp)
    this%mta(1:mdl%ntype, this%nsp) = mdl%multiplicity

    ! Zero out multiplicity for unused indicies in previously add
    ! species.
    do isp = 1, this%nsp
       this%mta(this%nat(isp) + 1:, isp)=0
    end do

    this%mtv => safemem_realloc(this%mtv, this%nv)
    this%mtv(oldnv + 1:this%nv) = mdl%multiplicity

    this%jvv => safemem_realloc(this%jvv, this%nv, this%nv)

    ! Setting triangle packing pointer.
    ivv = 0
    do iv2 = 1, this%nv
       do iv1 = 1, iv2
          ivv = ivv + 1
          this%jvv(iv1, iv2) = ivv
          this%jvv(iv2, iv1) = ivv
       end do
    end do

    this%mass => safemem_realloc(this%mass, this%nv)
    this%mass(oldnv + 1:this%nv) = mdl%mass

    this%densitysp => safemem_realloc(this%densitysp, this%nsp)
    this%densitysp(this%nsp) = density

    this%densityv => safemem_realloc(this%densityv, this%nv)
    this%densityv(oldnv + 1:this%nv) = mdl%multiplicity * density

    this%qsp => safemem_realloc(this%qsp, this%nsp)
    this%qsp(this%nsp) = sum(mdl%charge * mdl%multiplicity) / sqrt(this%temperature)

    this%qat => safemem_realloc(this%qat, ubound(this%mta, 1), this%nsp)
    this%qat(1:mdl%ntype, this%nsp) = mdl%charge / sqrt(this%temperature)

    this%qv => safemem_realloc(this%qv, this%nv)
    this%qv(oldnv + 1:this%nv) = mdl%charge / sqrt(this%temperature)

    this%qspv => safemem_realloc(this%qspv, this%nv)
    this%qspv(oldnv + 1:this%nv) = this%qsp(this%nsp)

    this%ljEpsilonVV => safemem_realloc(this%ljEpsilonVV, this%nv)
    this%ljEpsilonVV(oldnv + 1:this%nv) = mdl%ljEpsilon / this%temperature

    this%ljRminVV => safemem_realloc(this%ljRminVV, this%nv)
    this%ljRminVV(oldnv + 1:this%nv) = mdl%rmin / 2d0

    this%rma => safemem_realloc(this%rma, 3, maxval(this%mta), maxval(this%nat), this%nsp, .true., .false.)
    iv = 0
    do iat = 1, this%nat(this%nsp)
       do imlt = 1, this%mta(iat, this%nsp)
          iv = iv + 1
          this%rma(:, imlt, iat, this%nsp) = mdl%coord(:, iv)
       end do
    end do

    this%wlmvv => safemem_realloc(this%wlmvv, maxval(this%mta), this%nv, this%nv)
    call intramolecular_dist(this)

    ! Data arrays in the original 1D-RISM code had an offset of 0 in
    ! the first index (like C).  The offset is 1 in the current code
    ! (like most Fortran).
    this%wvv => safemem_realloc(this%wvv, this%nr, this%nv, this%nv)
    this%wzvv => safemem_realloc(this%wzvv, this%nr, this%nv, this%nv)
    this%zkvv => safemem_realloc(this%zkvv, this%nr, this%nv, this%nv)
    this%uvv => safemem_realloc(this%uvv, this%nr, this%nvv)
    this%ulrvv => safemem_realloc(this%ulrvv, this%nr, this%nvv)
    this%ulkvv => safemem_realloc(this%ulkvv, this%nr, this%nvv)
    this%hlrvv => safemem_realloc(this%hlrvv, this%nr, this%nvv)
    this%hlkvv => safemem_realloc(this%hlkvv, this%nr, this%nvv)
  end subroutine rism1d_potential_addSpecies


  !> Destroys a potential object
  !! @param[in,out] this Potential object.
  subroutine rism1d_potential_destroy(this)
    implicit none
    type(rism1d_potential), intent(inout) :: this
    this%theory = ''
    this%dielconst = 0
    this%smear = 0
    this%adbcor = 0
    this%nr = 0
    this%dr = 0
    this%dk = 0
    this%nsp = 0
    this%nv = 0
    this%nvv = 0
    this%kappa = HUGE(1d0)
    if (safemem_dealloc(this%nat) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate NAT")
    if (safemem_dealloc(this%namesp) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate NAMESP")
    if (safemem_dealloc(this%namev) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate NAMEV")
    if (safemem_dealloc(this%mta) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate MTA")
    if (safemem_dealloc(this%mtv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate MTV")
    if (safemem_dealloc(this%jvv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate JVV")
    if (safemem_dealloc(this%mass) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate MASS")
    if (safemem_dealloc(this%densitysp) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate DENSITYSP")
    if (safemem_dealloc(this%densityv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate DENSITYV")
    if (safemem_dealloc(this%qsp) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate QSP")
    if (safemem_dealloc(this%qat) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate QAT")
    if (safemem_dealloc(this%qv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate QV")
    if (safemem_dealloc(this%qspv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate QSPV")
    if (safemem_dealloc(this%ljEpsilonVV) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate EPSV")
    if (safemem_dealloc(this%ljRminVV) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate LJRMINVV")
    if (safemem_dealloc(this%rma) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate RMA")
    if (safemem_dealloc(this%wlmvv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate WLMVV")
    if (safemem_dealloc(this%wvv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate WVV")
    if (safemem_dealloc(this%wzvv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate WZVV")
    if (safemem_dealloc(this%zkvv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate XKVV")
    if (safemem_dealloc(this%uvv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate UVV")
    if (safemem_dealloc(this%ulrvv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate ULRVV")
    if (safemem_dealloc(this%ulkvv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate ULKVV")
    if (safemem_dealloc(this%hlrvv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate HLRVV")
    if (safemem_dealloc(this%hlkvv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate HLKVV")
    if (safemem_dealloc(this%duvv) /= 0) &
         call rism_report_error("POTENTIAL: failed to deallocate DUVV")
  end subroutine rism1d_potential_destroy


  !> Tabulates the intramolecular matrix
  !! and the site-site potential comprising
  !! the Coulomb and Lennard-Jones (LJ) terms,
  !! the latter for Composite Sphere (CS) sites
  !! as well as for usual point sites.
  !!  Units: charge[sqrt(kT A)], energy[kT], distance [A]
  subroutine  rism1d_potential_calc(this, charge)
    use constants, only : pi, KB
    use rism_util, only : spherical_bessel_j0, spherical_bessel_j1
    implicit none
#include "../xblas/f77/blas_namedconstants.fh"
    type(rism1d_potential), intent(inout) :: this
    _REAL_, intent(in) :: charge

    integer ::  id, imt, imt1, iat, iat1, iat2, isp, isp1, isp2, &
         iv, iv1, iv2, ivv
    integer :: ir

    logical :: found
    _REAL_ ::  dipmom, densitydm, densitydm2, ydiel, cdiel, &
         k, kl, krx, kry, krz, hck, wk, qvv, epsvv, ljRminVV, &
         r, rs, rs6, ri6, usr, usra
    _REAL_ :: alphavv, ri4
    _REAL_ ::  d0x(this%nv), d0y(this%nv), d1z(this%nv), &
         wkvv(this%nv, this%nv), erfc_test, err

    !............... shifting and orienting species dipoles ................
    call dipord(this)

    !............................ screen output ............................
    call rism_report_message('tabulating potentials and intramolecular functions ...')

    !............ calculating total Density and Density * D^2 of dipoles .............
    densitydm = 0.d0
    densitydm2 = 0.d0
    do isp = 1, this%nsp

       do id = 1, 3
          dipmom = 0.d0
          do iat = 1, this%nat(isp)
             do imt = 1, this%mta(iat, isp)
                dipmom = dipmom + charge * this%qat(iat, isp) * this%rma(id, imt, iat, isp)
             end do
          end do
          if (id /= 3 .and. abs(dipmom) > 1.d-8) then
             call rism_report_error('POTENT: dipole moment not in z-direction')
          end if
       end do

       if (abs(dipmom) > 1.d-8) then
          densitydm = densitydm + this%densitysp(isp)
          densitydm2 = densitydm2 + this%densitysp(isp) * dipmom**2
       end if

    end do
    !............ calculating dielectric asymptotics parameters ............
    ydiel = 4.d0 / 9.d0 * pi * densitydm2
    if (ydiel == 0.d0) then
       this%dielconst = 1.d0
       cdiel = 0.d0
    else if (this%theory == 'DRISM') then
       cdiel = ((this%dielconst - 1.d0) / ydiel - 3.d0) / densitydm
    else
       this%dielconst = 1.d0 + 3.d0 * ydiel
       cdiel = 0.d0
    end if


    !.... this%pcdiel

    this%pcdiel = densitydm * cdiel

    !.................. calculating inverse Debye length ...................
    this%kappa = 0.d0
    do isp = 1, this%nsp
       this%kappa = this%kappa + this%densitysp(isp) * charge * this%qsp(isp)**2
    end do
    this%kappa = sqrt(4.d0 * pi * this%kappa / this%dielconst)

    !................ calculating intramolecular distances .................
    iv2 = 0
    do isp2 = 1, this%nsp
       do iat2 = 1, this%nat(isp2)
          iv2 = iv2 + 1

          iv1 = 0
          do isp1 = 1, this%nsp
             do iat1 = 1, this%nat(isp1)
                iv1 = iv1 + 1

                do imt1 = 1, this%mtv(iv1)

                   if (isp1 /= isp2) then
                      this%wlmvv(imt1, iv1, iv2) = -1.d0
                   else
                      this%wlmvv(imt1, iv1, iv2) = sqrt(&
                           (this%rma(1, imt1, iat1, isp1) - this%rma(1, 1, iat2, isp2))**2 &
                           + (this%rma(2, imt1, iat1, isp1) - this%rma(2, 1, iat2, isp2))**2 &
                           + (this%rma(3, imt1, iat1, isp1) - this%rma(3, 1, iat2, isp2))**2)
                   end if
                end do
             end do
          end do
       end do
    end do



    ! Enumerating K-grid points.
    do ir = 1, this%nr
       k = (ir - 1) * this%dk
       hck = cdiel * exp(-(this%adbcor * k / 2.d0)**2)

       ! Calculating j0x, j0y, j1z.
       iv = 0
       do isp = 1, this%nsp
          do iat = 1, this%nat(isp)
             iv = iv + 1

             krx = k * this%rma(1, 1, iat, isp)
             if (krx == 0.d0) then
                d0x(iv) = 1.d0
             else
                d0x(iv) = spherical_bessel_j0(krx)
             end if

             kry = k * this%rma(2, 1, iat, isp)
             if (kry == 0.d0) then
                d0y(iv) = 1.d0
             else
                d0y(iv) = spherical_bessel_j0(kry)
             end if

             krz = k * this%rma(3, 1, iat, isp)
             if (krz == 0.d0) then
                d1z(iv) = 0.d0
             else
                d1z(iv) = spherical_bessel_j1(krz, err)
             end if

          end do
       end do

       ! Getting and reducing the intramolecular matrix Wvv(k).
       do iv2 = 1, this%nv
          do iv1 = 1, this%nv
             wkvv(iv1, iv2) = 0.d0
             do imt1 = 1, this%mtv(iv1)
                kl = k * this%wlmvv(imt1, iv1, iv2)
                if (this%wlmvv(imt1, iv1, iv2) < 0.d0) then
                   wk = 0.d0
                else if (kl == 0.d0) then
                   wk = 1.d0
                else
                   wk = spherical_bessel_j0(kl)
                end if
                wkvv(iv1, iv2) = wkvv(iv1, iv2) + wk
             end do
          end do
       end do

       ! Loading the intramolecular matrix Wvv(k).
       do iv2 = 1, this%nv
          do iv1 = 1, this%nv
             this%wvv(ir, iv1, iv2) = wkvv(iv1, iv2)
          end do
       end do

       !........................... getting Zvv(k) ............................
       !this should be symmetric.  Try an force it to be
       do iv2 = 1, this%nv
          do iv1 = 1, iv2
             this%zkvv(ir, iv1, iv2) = d0x(iv1) * d0y(iv1) * d1z(iv1) * hck &
                  * d0x(iv2) * d0y(iv2) * d1z(iv2)
             if (iv1 /= iv2) then
                this%zkvv(ir, iv2, iv1) = d0x(iv2) * d0y(iv2) * d1z(iv2) * hck &
                     * d0x(iv1) * d0y(iv1) * d1z(iv1)
             end if
          end do
       end do

       !.............. getting Wvv(k) + Density_v * Zvv(k) and k * Zvv(k) ...............
       do iv1 = 1, this%nv
          if (this%extra_precision == 0) then
             call DCOPY(this%nv, wkvv(iv1, :), 1, this%wzvv(ir, iv1, :), 1)
             call DAXPY(this%nv, this%densityv(iv1), this%zkvv(ir, iv1, :), 1, &
                  this%wzvv(ir, iv1, :), 1)
             call DSCAL(this%nv, k, this%zkvv(ir, iv1, :), 1)
          else
             call BLAS_DWAXPBY_X(this%nv, 1d0, wkvv(iv1, :), 1, this%densityv(iv1), this%zkvv(ir, iv1, :), 1, &
                  this%wzvv(ir, iv1, :), 1, BLAS_PREC_EXTRA)
             call BLAS_DAXPBY_X(this%nv, 0d0, wkvv(iv1, :), 1, k, this%zkvv(ir, iv1, :), 1, &
                  BLAS_PREC_EXTRA)
          end if
       end do
!!$       write(0,*) 'zkvv', ir, this%zkvv(ir, :,:)
    end do
!!$    write(0, '(a, 2(e24.16, 1x))') 'wkvv', wkvv(:, :)
!!$    stop
    ivv = 0
    do iv2 = 1, this%nv
       do iv1 = 1, iv2
          ivv = ivv + 1
          !ensure all elements are initialized
          this%uvv(1, ivv) = 0
          !.......................... Coulomb potential ..........................
          !.......................... 12-6 LJ potential ..........................
          qvv = charge * this%qv(iv1) * this%qv(iv2)
          found = .false.
          if (associated(this%modnb)) then
             do id = 1, size(this%modnbpair)
                if (this%modnbpair(id) == ivv) then
                   ljRminVV = this%modnb(id, 1)
                   epsvv = this%modnb(id, 2)
                   !call rism_report_message("(a, i10)", 'NB_Fix pair ', ivv)
                   write (6,'(a,a,a)') "# LJ 12-6 NBFIX: " , this%namev(iv1), this%namev(iv2)
                   write (6,'(i8, E20.12, E20.12)') ivv , ljRminVV, epsvv*KB*this%temperature
                   found = .true.
                   exit;
                end if
             end do
          end if
          if (.not. found) then
             ljRminVV = (this%ljRminVV(iv1) + this%ljRminVV(iv2))
             epsvv = sqrt(this%ljEpsilonVV(iv1) * this%ljEpsilonVV(iv2))
             if (associated(this%modnb)) then
                write (6,'(a, a, a)') "# LJ 12-6 Mixing Rule: " , this%namev(iv1), this%namev(iv2)
                write (6,'(i8, E20.12, E20.12)') ivv , ljRminVV, epsvv*KB*this%temperature
             end if
          end if
          do ir = 2, this%nr
             r = (ir - 1) * this%dr
             rs = r / ljRminVV
             rs = max(rs, 0.002d0)
             rs6 = rs**6
             ri6 = 1.d0 / rs6
             usr = (ri6 - 2.d0) * ri6
             usr = usr * epsvv
             this%uvv(ir, ivv) = qvv / r + usr
          end do
       end do
    end do


    if (associated(this%modnbtsf)) then
        call rism_report_message ("LJ 12-6-4 corrections.")
        ivv = 0
        call rism_report_message ("#NB_Fix pairs")
        do iv2 = 1, this%nv
           do iv1 = 1, iv2
              ivv = ivv + 1
              this%uvv(1, ivv) = 0
              qvv = charge * this%qv(iv1) * this%qv(iv2)
              found = .false.
              do id = 1, size(this%modnbpairtsf)
                if (this%modnbpairtsf(id) == ivv) then
                   ljRminVV = this%modnbtsf(id, 1)
                   epsvv = this%modnbtsf(id, 2)
                   alphavv = this%modnbtsf(id, 3)
                   write (6,'(a, a, a)') "# LJ 12-6-4 NBFIX: " , this%namev(iv1), this%namev(iv2)
                   write (6,'(i8, E20.12, E20.12, E20.12)') ivv , ljRminVV, epsvv*KB*this%temperature, alphavv
                   found = .true.
                   exit;
                end if
              end do
              if (found) then
                  do ir = 2, this%nr
                     r = (ir - 1) * this%dr
                     rs = r / ljRminVV
                     rs = max(rs, 0.002d0)
                     rs6 = rs**6
                     ri6 = 1.d0 / rs6
                     ri4 = rs**4
                     usr = (ri6 - 2.d0  - 2d0*r*r*alphavv) * ri6
                     usr = usr * epsvv
                     this%uvv(ir, ivv) = qvv / r + usr
                  end do
               else
                ljRminVV = (this%ljRminVV(iv1) + this%ljRminVV(iv2))
                epsvv = sqrt(this%ljEpsilonVV(iv1) * this%ljEpsilonVV(iv2))
                write (6,'(a,a,a)') "# LJ 12-6-4 Mixing Rule: " , this%namev(iv1), this%namev(iv2)
                write (6,'(i8, E20.12, E20.12, E20.12)') ivv , ljRminVV, epsvv*KB*this%temperature, 0d0
               endif
           end do
        end do
    endif

    !......... smeared Coulomb potential: r * Ulvv(r) and k * Ulvv(k) ..........
    ivv = 0
    do iv2 = 1, this%nv
       do iv1 = 1, iv2
          ivv = ivv + 1

          !ensure all elements are initialized
          this%ulrvv(1, ivv) = 0
          this%ulkvv(1, ivv) = 0
          qvv = charge * this%qv(iv1) * this%qv(iv2)
          do ir = 2, this%nr
             r = (ir - 1) * this%dr
             k = (ir - 1) * this%dk
             this%ulrvv(ir, ivv) = qvv * (1.d0 - erfc(r / this%smear))
             this%ulkvv(ir, ivv) = qvv * 4.d0 * pi * exp(-(this%smear * k / 2.d0)**2) / k
          end do

       end do
    end do

    !........ smeared Coulomb term of TCF: r * Hlvv(r) and k * Hlvv(k) .........
    ivv = 0
    do iv2 = 1, this%nv
       do iv1 = 1, iv2
          ivv = ivv + 1

          qvv = charge * this%qspv(iv1) * this%qspv(iv2)
          !ensure all elements are initialized
          this%hlrvv(1, ivv) = 0
          this%hlkvv(1, ivv) = 0
          do ir = 2, this%nr
             r = (ir - 1) * this%dr
             k = (ir - 1) * this%dk
             !for large grids with short Debye lengths the positive
             !exponent can overflow.  erfc is already zero at this
             !point.  So we test erfc to see if we can avoid the
             !exponent
             erfc_test = erfc(this%kappa * this%smear / 2.d0 + r / this%smear)
             if (erfc_test > sqrt(tiny(1d0))) then
                this%hlrvv(ir, ivv) = -qvv / this%dielconst &
                     * 0.5d0 * exp((this%smear * this%kappa / 2.d0)**2) &
                     * (exp(-this%kappa * r) * erfc(this%kappa * this%smear / 2.d0 - r / this%smear) &
                     - exp(this%kappa * r) * erfc_test)
             else
                this%hlrvv(ir, ivv) = -qvv / this%dielconst &
                     * 0.5d0 * exp((this%smear * this%kappa / 2.d0)**2) &
                     * (exp(-this%kappa * r) * erfc(this%kappa * this%smear / 2.d0 - r / this%smear))
             end if
             this%hlkvv(ir, ivv) = -qvv / this%dielconst &
                  * 4.d0 * pi * exp(-(0.5d0 * this%smear * k)**2) * k / (k**2 + this%kappa**2)
          end do
       end do
    end do

    return
  end subroutine rism1d_potential_calc


  !> Tabulates the r derivative of the site-site potential comprising
  !! the Coulomb and Lennard-Jones (LJ) terms.  The result is stored in duvv.
  !! This is only run once so the results may be reused.
  subroutine  rism1d_potential_derivative(this)
    implicit none
    type(rism1d_potential), intent(inout) :: this
    _REAL_ :: qvv, ljRminVV, epsvv, r, ri6, rs, rs6, usr
    integer :: ivv, iv1, iv2, id
    integer :: ir
    logical :: found

    if (.not.associated(this%duvv)) then
       this%duvv => safemem_realloc(this%duvv, this%nr, this%nvv)

       ivv = 0
       do iv2 = 1, this%nv
          do iv1 = 1, iv2
             ivv = ivv + 1
             !ensure all elements are initialized
             this%duvv(1, ivv) = 0
             !.......................... Coulomb potential ..........................
             !.......................... 12-6 LJ potential ..........................
             qvv = -this%qv(iv1) * this%qv(iv2)
             found = .false.
             if (associated(this%modnb)) then
                do id = 1, size(this%modnbpair)
                   if (this%modnbpair(id) == ivv) then
                      ljRminVV = this%modnb(id, 1)
                      epsvv = this%modnb(id, 2)
                      !                      call rism_report_message("(a, i10)", 'NB_Fix pair ', ivv)
                      found = .true.
                      exit;
                   end if
                end do
             end if
             if (.not. found) then
                ljRminVV = (this%ljRminVV(iv1) + this%ljRminVV(iv2))
                epsvv = sqrt(this%ljEpsilonVV(iv1) * this%ljEpsilonVV(iv2))
             end if
             do ir = 2, this%nr
                r = (ir - 1) * this%dr
                rs = r / ljRminVV
                rs = max(rs, 0.002d0)
                rs6 = rs**6
                ri6 = 1.d0 / rs6
                usr = (1d0 - ri6) * ri6
                usr = usr * epsvv * 12d0
                this%duvv(ir, ivv) = (qvv / r + usr) / r
             end do
          end do
       end do
    end if
  end subroutine rism1d_potential_derivative


  !> Returns the inverse Debye length (1/A) of the solvent. rism1d_potential_calc
  !! must be called first or huge(1d0) will be returned.
  !! @param[in,out] this rism1d object.
  function rism1d_potential_getInvDebyeLen(this) result(kappa)
    implicit none
    type(rism1d_potential), intent(inout) :: this
    _REAL_ :: kappa
    kappa = this%kappa
  end function rism1d_potential_getInvDebyeLen



!!!                               PRIVATE


  !> Check user parameters.
  subroutine sanity_check(this)
    implicit none
    type(rism1d_potential), intent(inout) :: this
    integer :: factored_nr
    !supported theory...
    if (.not.(this%theory.eq.'XRISM' .or. this%theory.eq.'DRISM')) then
       call rism_report_error("Only DRISM and XRISM theories are supported. Invalid theory: "//trim(this%theory))
    end if

    !check grid size
    if (this%nr <1) &
         call rism_report_error("NR must be > 1")
    !check the grid size prime factorization...
    factored_nr = this%nr
    do while(mod(factored_nr, 2) == 0)
       factored_nr = factored_nr / 2
    end do
    do while(mod(factored_nr, 3) == 0)
       factored_nr = factored_nr / 3
    end do
    do while(mod(factored_nr, 5) == 0)
       factored_nr = factored_nr / 5
    end do
    if (factored_nr /= 1)&
         call rism_report_error("NR must be a product of 2, 3 and 5 only")

    !grid spacing
    if (this%dr <=0) &
         call rism_report_error("DR must be > 0")

    !temperature
    if (this%temperature <=0) &
         call rism_report_error("TEMPERATURE must be > 0")

    !smear
    if (this%smear <= 0) &
         call rism_report_error("SMEAR must be > 0")

    !adbcor
    if (this%adbcor <= 0) &
         call rism_report_error("ADBCOR must be > 0")

    !dielconst
    if (this%dielconst <= 0 .and. this%theory .eq. 'DRISM') &
         call rism_report_error("dielectric constant must be > 0 for DRISM")

  end subroutine sanity_check


  !> Ordering the species dipole moments:
  !! Shifting the molecule center of |charge| to the coordinate origin
  !! and orienting the dipole moment along the OZ axis.
  subroutine  dipord(this)
    use quaternion
    use rism_util, only : principalAxes, orientToPrincipalAxes, cross
    implicit none
    type(rism1d_potential), intent(inout) :: this
    integer ::  imt, iat, isp, id, iv
    !rqcm : center of absolute charge vector
    !rdm : dipole moment vector
    !qcm : total absolute charge
    !dm : magnitude of dipole moment
    _REAL_ ::  rqcm(3), rdm(3), qcm, dm
    !quat : quaternion for rotations
    !angle : the angle to rotate about
    !rotvec : the vector to rotate about
    !checkvec : vector used to determine the sign of angle
    !zaxis : the zaxis
    _REAL_ :: quat(4), angle, rotvec(3), checkvec(3), &
         zaxis(3)=(/0d0, 0d0, 1d0/)!, yaxis(3)=(/0d0, 1d0, 0d0/), xaxis(3)=(/1d0, 0d0, 0d0/), &
    !pa : holds principal axes
    _REAL_ :: pa(3, 3)
    !ratu : temporary array for coordinates while aligning principal axes
    !mass : variable for mass in PA calculations.  Will be absolute
    !       charge in this case
    _REAL_, pointer :: ratu(:,:) => NULL(), mass(:) => NULL()

    do isp = 1, this%nsp

       !
       !shifting molecule to center of absolute charge
       !

       !total absolute charge
       qcm = 0.d0
       do iat = 1, this%nat(isp)
          do imt = 1, this%mta(iat, isp)
             qcm = qcm + abs(this%qat(iat, isp))
          end do
       end do
       if (qcm == 0.d0)  cycle

       !get center of |charge|
       do id = 1, 3
          rqcm(id) = 0.d0
       end do
       do iat = 1, this%nat(isp)
          do imt = 1, this%mta(iat, isp)
             do id = 1, 3
                rqcm(id) = rqcm(id) + abs(this%qat(iat, isp)) * this%rma(id, imt, iat, isp)
             end do
          end do
       end do
       do id = 1, 3
          rqcm(id) = rqcm(id) / qcm
       end do

       !translate
       do iat = 1, this%nat(isp)
          do imt = 1, this%mta(iat, isp)
             do id = 1, 3
                this%rma(id, imt, iat, isp) = this%rma(id, imt, iat, isp) - rqcm(id)
             end do
          end do
       end do

       !
       ! Orient dipole moment on z-axis.
       !

       ! Calculate dipole moment.
       do id = 1, 3
          rdm(id) = 0.d0
       end do
       do iat = 1, this%nat(isp)
          do imt = 1, this%mta(iat, isp)
             do id = 1, 3
                rdm(id) = rdm(id) + this%qat(iat, isp) * this%rma(id, imt, iat, isp)
             end do
          end do
       end do

       ! Magnitude of the dipole moment.
       dm = sqrt(sum(rdm**2))

       if (dm < 1.d-16)  cycle

       ! Rotate the dipole moment to the z-axis.

       ! For DRISM it is necessary to rotate the molecule such that the
       ! dipole is oriented along the z-axis.  There is also a
       ! dependence on the position of the sites in the xy-plane and it
       ! is not clear what the correct orientation should be, if there
       ! is one.  In light of this, we attempt to at least be
       ! consistent such that the same molecule always ends up with the
       ! same orientation.

       ! This is done by finding the vector perpendicular to the dipole
       ! and the z-axis and then rotating about this vector by the
       ! angle between the dipole and z-axis. The orientation in the
       ! XY-plane is then fixed by aligning the XY principal axis
       ! (using absolute charge instead of mass) on the x-axis.  It is
       ! possible to get rotations of 180 degrees, which does not
       ! change the final solution.

       ! Orient dipole.
       rotvec = cross(zaxis, rdm / dm)
       checkvec = cross(rotvec, rdm / dm)
       angle = sign(acos(dot_product(rdm / dm, zaxis)), dot_product(checkvec, zaxis))
       quat = quaternionFromEulerAxis(angle, rotvec)
       do iat = 1, this%nat(isp)
          do imt = 1, this%mta(iat, isp)
             call rotate_quat(this%rma(:, imt, iat, isp), quat)
          end do
       end do

       !orient XY-principal axis

       !1) trans coordinates to an contiguous array
       ratu=> safemem_realloc(ratu, 3, sum(this%mta(:, isp)))
       mass=> safemem_realloc(mass, sum(this%mta(:, isp)))
       iv = 0
       do iat = 1, this%nat(isp)
          do imt = 1, this%mta(iat, isp)
             iv = iv + 1
             ratu(:, iv) =this%rma(:, imt, iat, isp)
             mass(iv) = abs(this%qat(iat, isp))
          end do
       end do
       !avoid rotating the dipole from the proper alignment
       ratu(3, :) = 0d0

       !2) calculate PA and rotate
       call principalAxes(ratu, mass, pa)
       call orientToPrincipalAxes(ratu, pa, quat)

       !3) transfer coordinates back
       iv = 0
       do iat = 1, this%nat(isp)
          do imt = 1, this%mta(iat, isp)
             iv = iv + 1
             this%rma(1:2, imt, iat, isp) = ratu(1:2, iv)
          end do
       end do

       !Alternately, we first align the xy-dipole component along the
       !x-axis and the align the xz-component (there is no more
       !y-component) along the z-axis.  This was the original method
       !but could lead to arbitrary rotations in the xz-plane of 90
       !degrees depending on platform/initial orientation.  Such
       !rotations, however, do not change the final results
!!$       !rotate about z-axis to align xy-dipole component with the x-axis.
!!$       angle = -acos(rdm(1)/sqrt(sum(rdm**2)))
!!$       write(0,*) "angle 1", angle, sign(zaxis, rdm(2))
!!$       !The sign of the y-component determines the direction of the axis of rotation
!!$       quat = new_quat(angle, sign(zaxis, rdm(2)))
!!$       do iat = 1, this%nat(isp)
!!$          do imt = 1, this%mta(iat, isp)
!!$             call rotate_quat(this%rma(:, imt, iat, isp), quat)
!!$          end do
!!$       end do
!!$       do id = 1, 3
!!$          rdm(id) = 0.d0
!!$       end do
!!$       do iat = 1, this%nat(isp)
!!$          do imt = 1, this%mta(iat, isp)
!!$             do id = 1, 3
!!$                rdm(id) = rdm(id) + this%qat(iat, isp)*this%rma(id, imt, iat, isp)
!!$             end do
!!$                write(0,*) "COORD", imt, iat, isp, this%rma(:, imt, iat, isp)
!!$          end do
!!$       end do
!!$       write(0,*) "RDM x", rdm
!!$       !rotate about y-axis to align xz-dipole component with the
!!$       !z-axis
!!$       angle = -acos(rdm(3) / sqrt(sum(rdm**2)))
!!$       write(0,*) "angle 2", angle
!!$       !now the dipole is positive x and + / - z. So, the y unit
!!$       !vector is always the correct rotation axis
!!$       quat = new_quat(angle, yaxis)
!!$       do iat = 1, this%nat(isp)
!!$          do imt = 1, this%mta(iat, isp)
!!$             call rotate_quat(this%rma(:, imt, iat, isp), quat)
!!$          end do
!!$       end do
    end do
    if (safemem_dealloc(ratu) /= 0) then
       call rism_report_error("DIPORD: deallocation failed: RATU")
    end if
    if (safemem_dealloc(mass) /= 0) then
       call rism_report_error("DIPORD: deallocation failed: MASS")
    end if
  end subroutine dipord


  !> Calculate the intramolecular distance matrix.
  !! @param[in,out] this rism1d object
  subroutine intramolecular_dist(this)
    implicit none
    type(rism1d_potential), intent(inout) :: this
    integer :: iv1, iv2, iat1, iat2, imlt, isp1, isp2
    iv2 = 0
    do isp2 = 1, this%nsp
       do iat2 = 1, this%nat(isp2)
          iv2 = iv2 + 1

          iv1 = 0
          do isp1 = 1, this%nsp
             do iat1 = 1, this%nat(isp1)
                iv1 = iv1 + 1

                do imlt = 1, this%mtv(iv1)

                   if (isp1 /= isp2) then
                      this%wlmvv(imlt, iv1, iv2) = -1.d0
                   else
                      this%wlmvv(imlt, iv1, iv2) = sqrt(&
                           (this%rma(1, imlt, iat1, isp1) - this%rma(1, 1, iat2, isp2))**2 &
                           + (this%rma(2, imlt, iat1, isp1) - this%rma(2, 1, iat2, isp2))**2 &
                           + (this%rma(3, imlt, iat1, isp1) - this%rma(3, 1, iat2, isp2))**2)
                   end if
                end do
             end do
          end do
       end do
    end do
  end subroutine intramolecular_dist


  !> Initialize NB modification / parse the nbfix file
  !! @param[in,out] this rism1d object
  !! @param[in] nbfix the filename of nonbond_fix file
  subroutine rism1d_potential_initNBFix(this, nbfix)
    use constants, only: KB
    use rism_util, only : freeunit
    implicit none
    type(rism1d_potential), intent(inout) :: this
    character(len=*), intent(in) :: nbfix
    integer :: unit, iostat, nline
    integer, parameter :: CLEN = 256
    character(len = CLEN) :: line
    _REAL_ :: beta

    beta = 1d0 / (KB * this%temperature)

    ! Read the nonbond-fix input file.
    call rism_report_message('reading NBFIX: '//nbfix)
    unit = freeUnit()
    open(unit = unit, file = nbfix, status = 'old', iostat = iostat)
    if (iostat /= 0) &
         call rism_report_error('(a, i4)', "Could not open"//nbfix//":", iostat)

    ! Skip the comment on the head and count the number of lines.
    nline = 0
    do while (.true.)
       read(unit, *, iostat = iostat) line
       if (iostat /= 0) exit
       line = adjustl(line)
       if (line(1:1) == '#') cycle
       nline = nline + 1
    end do

    this%modnbpair => safemem_realloc(this%modnbpair, nline)
    this%modnb => safemem_realloc(this%modnb, nline, 2)
    ! Read the modified nonbond parameters.
    rewind(unit)
    nline = 1
    do while (.true.)
       read(unit, *, iostat = iostat) line
       if (iostat /= 0) exit
       line = adjustl(line)
       if (line(1:1) == '#') cycle
       backspace(unit)
       read(unit, *, iostat = iostat) this%modnbpair(nline), this%modnb(nline, 1:2)
       this%modnb(nline, 2) = this%modnb(nline, 2) * beta
       nline = nline + 1
    end do
    close(unit)
  end subroutine rism1d_potential_initNBFix




  !> Initialize NB modification / parse the nbfix file
  !! @param[in,out] this rism1d object
  !! @param[in] nbfix the filename of nonbond_fix file
  subroutine rism1d_potential_initNBFixtsf(this, nbfix)
    use constants, only: KB
    use rism_util, only : freeunit
    implicit none
    type(rism1d_potential), intent(inout) :: this
    character(len=*), intent(in) :: nbfix
    integer :: unit, iostat, nline
    integer, parameter :: CLEN = 256
    character(len = CLEN) :: line
    _REAL_ :: beta

    beta = 1d0 / (KB * this%temperature)

    ! Read the nonbond-fix input file.
    call rism_report_message('reading NBFIX: '//nbfix)
    unit = freeUnit()
    open(unit = unit, file = nbfix, status = 'old', iostat = iostat)
    if (iostat /= 0) &
         call rism_report_error('(a, i4)', "Could not open"//nbfix//":", iostat)

   ! Skip the comment on the head and count the number of lines.
    nline = 0
    do while (.true.)
       read(unit, *, iostat = iostat) line
       if (iostat /= 0) exit
       line = adjustl(line)
       if (line(1:1) == '#') cycle
       nline = nline + 1
    end do

    this%modnbpairtsf => safemem_realloc(this%modnbpairtsf, nline)
    this%modnbtsf => safemem_realloc(this%modnbtsf, nline, 3)
   ! Read the modified nonbond parameters.
    rewind(unit)
    nline = 1
    do while (.true.)
       read(unit, *, iostat = iostat) line
       if (iostat /= 0) exit
       line = adjustl(line)
       if (line(1:1) == '#') cycle
       backspace(unit)
       read(unit, *, iostat = iostat) this%modnbpairtsf(nline), this%modnbtsf(nline,1:3)
       this%modnbtsf(nline, 2) = this%modnbtsf(nline, 2) * beta
       nline = nline + 1
    end do
    close(unit)
  end subroutine rism1d_potential_initNBFixtsf


end module rism1d_potential_c
