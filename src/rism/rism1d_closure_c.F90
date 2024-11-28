!<compile=optimized>

! The 3D-RISM-KH software found here is copyright (c) 2010-2012 by
! Andriy Kovalenko, Tyler Luchko, Takeshi Yamazaki and David A. Case.
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


!> Closure super class for 1D-RISM.  Closure sub-classes (i.e.,
!! actual closure implementations) are registered here.  Subroutine
!! calls then call the appropriate subroutine of the
!! subclass. interface.  This is an explicit implementation of class
!! inheritance. See V. K. Decyk, C. D. Norton, B. K. Szymanski.  How
!! to express C++ concepts in Fortran 90. Scientific Programming. 6,
!! 363-390 (1997).
!!
!! Some closure independent properties are calculated within this class. Uvv is
!! the site-site potential.
module rism1d_closure_c
  ! Add new closure modules here.
  use rism1d_kh_c
  use rism1d_hnc_c
  use rism1d_py_c
  use rism1d_mv0_c
  use rism1d_psen_c
  use rism1d_nub_c
  use rism1d_polyt_c
  use rism1d_devb_c
  use rism1d_potential_c
  use safemem
  implicit none

  integer, private, parameter :: maxep0 = 20
  type rism1d_closure
     ! r-space contribution to free energy and pressure (free energy route).
     _REAL_ :: FE_press_r = HUGE(1d0)
     ! k-space pressure contribution.
     _REAL_ :: pressk = HUGE(1d0)
     ! k-space free energy contribution.
     _REAL_ :: fek = HUGE(1d0)
     ! Site-site susceptibility.
     _REAL_, pointer :: xvv(:,:,:) => NULL()
     ! Site-site succeseptibility temperature derivative.
     _REAL_, pointer :: xvv_dT(:,:,:) => NULL()
     ! k-space representation of Cvv.
     _REAL_, pointer :: cvvk(:,:) => NULL()

     type(rism1d_potential), pointer :: pot => NULL()
     type(rism1d_kh), pointer :: kh => NULL()
     type(rism1d_hnc), pointer :: hnc => NULL()
     type(rism1d_py), pointer :: py => NULL()
     type(rism1d_mv0), pointer :: mv0 => NULL()
     type(rism1d_psen), pointer :: psen => NULL()
     type(rism1d_nub), pointer :: nub
     type(rism1d_polyt), pointer :: polyt
     type(rism1d_devb), pointer :: devb => NULL()
     character(len = 8) :: type
  end type rism1d_closure

contains

  !> Creates a new closure object of the requested type.
  !! @param[in,out] this The closure object..
  !! @param[in] type One of 'KH', 'HNC', 'MV0', 'PSEN', 'POLYT', 'V*', and 'NUB'
  !!  order of the  PSE-n closure.
  !! @param[in] pot rism1d_potential object.  Must be initialized.
  !! @param[in] coeff (optional) Coefficients for the selected closure.
  subroutine rism1d_closure_new(this, type, pot, coeff)
    implicit none
    type(rism1d_closure), intent(inout) :: this
    type(rism1d_potential), target, intent(in)  :: pot
    character(len = *), intent(in) :: type
    integer :: order, iostat
    _REAL_, optional, intent(in) :: coeff(:)

    ! Reset press(r/k) and fe(r/k) so they will be recalculated.
    this%FE_press_r = HUGE(1d0)
    this%pressk = HUGE(1d0)
    this%fek = HUGE(1d0)

    this%pot => pot

    this%type = trim(type)
    if (type .eq. "KH") then
       allocate(this%kh)
       call rism1d_kh_new(this%kh)
    else if (index(type,"PSE") == 1) then
       read(type(4:),*, iostat=iostat) order
       if (iostat /= 0) &
            call rism_report_error(trim(type)//" not a valid closure")
       allocate(this%psen)
       call rism1d_psen_new(this%psen,order)
    else if (type .eq. "POLYT") then
       if (.not. present(coeff)) &
            call rism_report_error("coefficient array must be defined for the Poly-T closure.")
       allocate(this%polyt)
       call rism1d_polyt_new(this%polyt,coeff)
    else if (type .eq. "HNC") then
       allocate(this%hnc)
       call rism1d_hnc_new(this%hnc)
    else if (type .eq. "PY") then
       allocate(this%py)
       call rism1d_py_new(this%py)
    else if (type .eq. "MV0") then
       allocate(this%mv0)
       call rism1d_mv0_new(this%mv0)
    else if (type .eq. "NUB") then
       allocate(this%nub)
       call rism1d_nub_new(this%nub)
    else if (type .eq. "DEVB") then
       allocate(this%devb)
       if (.not. present(coeff)) &
            call rism_report_error("coefficient array must be defined for the DEVB closure.")
       call rism1d_devb_new(this%devb,coeff)
    else
       call rism_report_error(trim(type)//" not a valid closure")
    end if
  end subroutine rism1d_closure_new


  !> Returns a identifier string for the closure type.
  !! @param[in] this The closure object.
  function rism1d_closure_type(this) result(type)
    implicit none
    type(rism1d_closure), intent(in) :: this
    character(len = 8) :: type
    type = this%type
  end function rism1d_closure_type


  !> Returns true if the full potential is being used.  This is always
  !! true except for MV0.
  !! @param[in] this The closure object.
  function rism1d_closure_isCharged(this) result(charged)
    implicit none
    type(rism1d_closure), intent(in) :: this
    logical :: charged
    if (associated(this%MV0)) then
       charged = rism1d_mv0_charged(this%mv0)
       return
    end if
    charged = .true.
  end function rism1d_closure_isCharged


  !! Set to true to use the full potential, false otherwise.  This
  !! only applies to MV0.
  !! @param[in] this The closure object.
  !! @param[in] charged True to use full potential.
  subroutine rism1d_closure_useCharged(this,charged)
    implicit none
    type(rism1d_closure), intent(inout) :: this
    logical, intent(in) :: charged
    if (associated(this%MV0)) then
       call rism1d_mv0_usecharged(this%mv0, charged)
    end if
  end subroutine rism1d_closure_useCharged


  !! Calculates Gvv from Uvv, Hvv, and Cvv using the associated closure
  !! IN:
  !!    this : The closure object.
  !!    gvv  : site-site pair correlation function
  !!    hvv  : site-site total correlation function
  !!    cvv  : site-site direct correlation function
  subroutine rism1d_closure_gvv(this,gvv, hvv, cvv)
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(out) :: gvv(:,:)
    _REAL_, intent(in) :: hvv(:,:),cvv(:,:)
    if (associated(this%KH)) then
       call rism1d_kh_gvv(this%kh,gvv,this%pot%uvv,hvv,cvv)
    else if (associated(this%PSEN)) then
       call rism1d_psen_gvv(this%psen,gvv,this%pot%uvv,hvv,cvv)!
    else if (associated(this%POLYT)) then
       call rism1d_polyt_gvv(this%polyt,gvv,this%pot%uvv,hvv,cvv)
    else if (associated(this%HNC)) then
       call rism1d_hnc_gvv(this%hnc,gvv,this%pot%uvv,hvv,cvv)
    else if (associated(this%PY)) then
       call rism1d_py_gvv(this%py,gvv,this%pot%uvv,hvv,cvv)
    else if (associated(this%MV0)) then
       call rism1d_mv0_gvv(this%mv0,gvv,this%pot%uvv,hvv,cvv)
    else if (associated(this%NUB)) then
       call rism1d_nub_gvv(this%nub,gvv,this%pot%uvv,hvv,cvv)
    else if (associated(this%DEVB)) then
       call rism1d_devb_gvv(this%devb,gvv,this%pot%uvv,hvv,cvv)
    end if
  end subroutine rism1d_closure_gvv


  !! Calculates temperature derivative Gvv from Gvv, Uvv, Hvv_dT, and Cvv_dT using
  !! the associated closure
  !! IN:
  !!    this   : The closure object.
  !!    gvv_dT  : site-site temperature derivative pair correlation function
  !!    gvv    : site-site pair correlation function
  !!    cvv  : site-site direct correlation function
  !!    hvv_dT  : site-site temperature derivative total correlation function
  !!    cvv_dT  : site-site temperature derivative direct correlation function

  subroutine rism1d_closure_gvv_dT(this,gvv_dT, gvv, cvv, hvv_dT, cvv_dT)
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(out) :: gvv_dT(:,:)
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:), hvv_dT(:,:),cvv_dT(:,:)
    if (associated(this%KH)) then
       call rism1d_kh_gvv_dT(this%kh,gvv_dT,this%pot%uvv,gvv,hvv_dT,cvv_dT)
    else if (associated(this%PSEN)) then
       call rism1d_psen_gvv_dT(this%psen,gvv_dT,this%pot%uvv,gvv,cvv,hvv_dT,cvv_dT)
    else if (associated(this%POLYT)) then
       !       call rism1d_polyt_gvv_dT(this%polyt,gvv_dT,this%pot%uvv,gvv,cvv,hvv_dT,cvv_dT)
    else if (associated(this%HNC)) then
       call rism1d_hnc_gvv_dT(this%hnc,gvv_dT,this%pot%uvv,gvv,hvv_dT,cvv_dT)
    else
       call rism_report_error("Temperature derivative not supported for "//this%type)
    end if
  end subroutine rism1d_closure_gvv_dT


  !! Calculates Bvv (bridge function) from Uvv, Gvv, and Cvv using the associated closure
  !! IN:
  !!    this : The closure object.
  !!    gvv  : site-site total correlation function
  !!    cvv  : site-site direct correlation function
  !! OUT:
  !!    bvv  : site-site birdge function

  function rism1d_closure_bvv(this, gvv, cvv) result(bvv)
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, pointer :: bvv(:,:)
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:)
    nullify(bvv)
    if (associated(this%KH)) then
       bvv=> rism1d_kh_bvv(this%kh,this%pot%uvv,gvv,cvv)
    else if (associated(this%PSEN)) then
       bvv=> rism1d_psen_bvv(this%psen,this%pot%uvv,gvv,cvv)
    else if (associated(this%POLYT)) then
       bvv=> rism1d_polyt_bvv(this%polyt,this%pot%uvv,gvv,cvv)
    else if (associated(this%HNC)) then
       bvv=> rism1d_hnc_bvv(this%hnc,this%pot%uvv,gvv,cvv)
    else if (associated(this%PY)) then
       bvv=> rism1d_py_bvv(this%py,this%pot%uvv,gvv,cvv)
    else if (associated(this%MV0)) then
       bvv=> rism1d_mv0_bvv(this%mv0,this%pot%uvv,gvv,cvv)
    else if (associated(this%NUB)) then
       bvv=> rism1d_nub_bvv(this%nub,this%pot%uvv,gvv,cvv)
    else if (associated(this%DEVB)) then
       bvv=> rism1d_devb_bvv(this%devb,this%pot%uvv,gvv,cvv)
    end if
  end function rism1d_closure_bvv


  !! Frees memory and resets object state
  !! IN:
  !!    this : The closure object.

  subroutine rism1d_closure_destroy(this)
    use safemem
    implicit none
    type(rism1d_closure), intent(inout) :: this
    !............... reset press(r/k) and fe(r/k) so they will be recalculated
    this%FE_press_r = HUGE(1d0)
    this%pressk = HUGE(1d0)
    this%fek = HUGE(1d0)
    if (associated(this%KH)) then
       call rism1d_kh_destroy(this%kh)
       deallocate(this%kh)
    end if
    if (associated(this%PSEN)) then
       call rism1d_psen_destroy(this%psen)
       deallocate(this%psen)
    end if
    if (associated(this%POLYT)) then
       call rism1d_polyt_destroy(this%polyt)
       deallocate(this%polyt)
    end if
    if (associated(this%HNC)) then
       call rism1d_hnc_destroy(this%hnc)
       deallocate(this%hnc)
    end if
    if (associated(this%PY)) then
       call rism1d_py_destroy(this%py)
       deallocate(this%py)
    end if
    if (associated(this%MV0)) then
       call rism1d_mv0_destroy(this%mv0)
       deallocate(this%mv0)
    end if
    if (associated(this%NUB)) then
       call rism1d_nub_destroy(this%nub)
       deallocate(this%nub)
    end if
    if (associated(this%DEVB)) then
       call rism1d_devb_destroy(this%devb)
       deallocate(this%devb)
    end if
    nullify(this%pot)
    if (safemem_dealloc(this%xvv) /= 0) then
       call rism_report_error("deallocating Xvv failed in closure")
    end if
    if (safemem_dealloc(this%xvv_dT) /= 0) then
       call rism_report_error("deallocating Xvv_dT failed in closure")
    end if
    if (safemem_dealloc(this%cvvk) /= 0) then
       call rism_report_error("deallocating Cvvk failed in closure")
    end if
  end subroutine rism1d_closure_destroy


  !! Returns the compressibility of the solvent [A^3/kT]
  !! IN:
  !!    this : rism1d closure object
  !!    cvv  : direct correlation function
  !! OUT:
  !!    isothermal compressibility temperature derivative

  function rism1d_closure_getCompressibility(this,cvv) result(xikt)
    use constants, only : pi
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: cvv(:,:)
    _REAL_ :: xikt
    _REAL_ :: r, ck0, ck0r
    integer :: ir, ivv, iv2, iv1, msym
    ck0 = 0.d0
    do ir=2,this%pot%nr
       r = (ir-1)*this%pot%dr
       ck0r = 0.d0
       ivv = 0
       do iv2=1,this%pot%nv
          do iv1=1,iv2
             ivv = ivv + 1
             if (iv1 == iv2)  then
                msym = 1
             else
                msym = 2
             end if
             ck0r = ck0r + msym*this%pot%densityv(iv1)*this%pot%densityv(iv2)*cvv(ir,ivv)
          end do
       end do
       ck0 = ck0 + r**2*ck0r
    end do
    ck0 = ck0 * 4.d0*pi*this%pot%dr
    xikt = 1.d0 / (sum(this%pot%densitysp)-ck0)
  end function rism1d_closure_getCompressibility


  !! Returns the compressibility temperature derivative of the solvent [A^3/K]
  !! IN:
  !!    this : rism1d closure object
  !!    cvv  : direct correlation function
  !!    cvv_dT  : direct correlation function temperature derivative
  !! OUT:
  !!    isothermal compressibility temperature derivative

  function rism1d_closure_getCompressibility_dT(this,cvv,cvv_dT) result(xikt_dT)
    use constants, only : pi
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: cvv(:,:),cvv_dT(:,:)
    _REAL_ :: xikt_dT
    _REAL_ :: xikt
    _REAL_ :: r, ck0, ck0r
    integer :: ir, ivv, iv2, iv1, msym

    xikt = rism1d_closure_getCompressibility(this,cvv)
    ck0 = 0.d0
    do ir=2,this%pot%nr
       r = (ir-1)*this%pot%dr
       ck0r = 0.d0
       ivv = 0
       do iv2=1,this%pot%nv
          do iv1=1,iv2
             ivv = ivv + 1
             if (iv1 == iv2)  then
                msym = 1
             else
                msym = 2
             end if
             ck0r = ck0r + msym*this%pot%densityv(iv1)*this%pot%densityv(iv2)*cvv_dT(ir,ivv)
          end do
       end do
       ck0 = ck0 + r**2*ck0r
    end do
    ck0 = ck0 * 4.d0*pi*this%pot%dr
    xikt_dT = (-xikt + xikt**2*ck0)
  end function rism1d_closure_getCompressibility_dT


  !! Returns the the extrapolated value for
  !! DelHv=-Lim_k->0 ( Sum_v1 Qv1*Xv1v2(k)4pi/k^2 )
  !! This is the background charge contribution for periodic and
  !! non-periodic RISM without Q_solute, Q (Eq. 20 from Kovalenko &
  !! Hirata, JCP 2000).
  !! IN:
  !!    this : rism1d closure object
  !!    hvvk  : k-space total correlation function
  function rism1d_closure_get3DRISMbackground(this,hvvk) result(background)
    use constants, only : pi
    use rism_util, only : polynomialInterpolation_progressive
    implicit none
#include "../xblas/f77/blas_namedconstants.fh"
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: hvvk(:, :)
    _REAL_ :: background(this%pot%nv)
    integer :: ivv, iv1, iv2, ir
    _REAL_ :: k, ep0(maxep0), cep0(0:maxep0), err0

    call rism1d_closure_calcXvv(this, hvvk)
    
    do ir = 2, maxep0 +1
       ep0(ir - 1) = (ir - 1) * this%pot%dk
    end do
    
    do iv2 = 1, this%pot%nv
       ! write(6,'(a,i4)') 'background extrapolation for solvent site ', iv2
       ! write(6,'(a,5f12.7)') '    k:         ', ep0(1:5)
       ! do iv1 = 1, this%pot%nv
       !     write(6,'(a,i3,a,5f12.7)') '    alpha: ',iv1,':', &
       !          -4d0*pi*this%pot%qv(iv1)*this%xvv(2:6,iv1,iv2)
       ! end do
       call DGEMV("N", maxep0 + 1, this%pot%nv, - 4d0 * pi, this%xvv(:maxep0 + 1, :, iv2), maxep0 + 1, &
            this%pot%qv, 1, 0d0, cep0, 1)
       ! write(6,'(a,5f12.7)') '    sum:       ',cep0(1:5)
       cep0(1:maxep0) = cep0(1:maxep0) / ep0**2
       ! write(6,'(a,5f12.7)') '    sum/k^2:   ',cep0(1:5)
       call polynomialInterpolation_progressive (ep0, cep0(1:maxep0), maxep0, 0.d0, background(iv2), err0)
       ! write(6,'(a,2f12.5)') '    Value/error at zero: ', background(iv2), err0
    end do
  end function rism1d_closure_get3DRISMbackground

  !! Returns the the extrapolated value for
  !! DelHv=-Lim_k->0 ( Sum_v1 Qv1*Xv1v2(k)4pi/k^2 - hlkv0 )
  !! This is the entire h(k=0) for non-periodic RISM without Q_solute,
  !! combining \Delta Q (Eq. 20 from Kovalenko & Hirata, JCP 2000) and
  !! the long-range asymptotics value at k=0 (Eq. 11 from Kaminski et
  !! al., JPCA, 2010 and Eq. 10a from Genheden et al. JPCB, 2010).
  !! IN:
  !!    this : rism1d closure object
  !!    hvvk  : k-space total correlation function
  function rism1d_closure_getDelHvLimit(this,hvvk) result(delhv0)
    use constants, only : pi
    use rism_util, only : polynomialInterpolation_progressive
    implicit none
#include "../xblas/f77/blas_namedconstants.fh"
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: hvvk(:, :)
    _REAL_ :: delhv0(this%pot%nv)
    ! hlkv0 : molecular version of hlkvv for 3D-RISM
    _REAL_ :: hlkv0(maxep0 + 1, this%pot%nv)
    integer :: ivv, iv1, iv2, ir
    _REAL_ :: k, ep0(maxep0), cep0(0:maxep0), err0

    ! Getting the long-range asymptotic function of H(k) for 3D-RISM.
    do iv2 = 1, this%pot%nv
       do ir = 2, maxep0 + 1
          k = (ir - 1) * this%pot%dk
          hlkv0(ir, iv2) = -this%pot%qspv(iv2) / this%pot%dielconst &
               * 4.d0 * pi * exp( - (0.5d0 * this%pot%smear * k)**2) &
               / (k**2 + this%pot%kappa**2)
       end do
    end do
    
    call rism1d_closure_calcXvv(this, hvvk)
    
    do ir = 2, maxep0 +1
       ep0(ir - 1) = (ir - 1) * this%pot%dk
    end do
    
    do iv2 = 1, this%pot%nv
       call DGEMV("N", maxep0 + 1, this%pot%nv, - 4d0 * pi, this%xvv(:maxep0 + 1, :, iv2), maxep0 + 1, &
            this%pot%qv, 1, 0d0, cep0, 1)
       cep0(1:maxep0) = cep0(1:maxep0) / ep0**2
       call DAXPY(maxep0, - 1d0, hlkv0(2:maxep0 + 1, iv2), 1, &
            cep0(1:maxep0), 1)
       call polynomialInterpolation_progressive (ep0, cep0(1:maxep0), maxep0, 0.d0, delhv0(iv2), err0)
    end do
  end function rism1d_closure_getDelHvLimit


  !! Returns the the extrapolated value for the temperature derivative
  !! (T*d/dT) of DelHv0.
  !! DelHv_dT=-Lim_k->0 ( Sum_v1 Qv1*Xv1v2_dT(k)4pi/k^2 - hlkv0 )
  !! This is used by 3D-RISM long range asymptotics
  !! IN:
  !!    this   : rism1d closure object.
  !!    hvvk    : k-space total correlation function
  !!    hvvk_dT : k-space temperature derivative total correlation function
  !! OUT:
  !!    temperature derivative delhv0

  function rism1d_closure_getDelHvLimit_DT(this,hvvk,hvvk_dT) result(delhv0_dT)
    use constants, only : pi
    use rism_util, only : polynomialInterpolation_progressive
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: hvvk(:,:), hvvk_dT(:,:)
    _REAL_ :: delhv0_dT(this%pot%nv)
    !hlkv0 : molecular version of hlkvv for 3D-RISM
    _REAL_ :: hlkv0(maxep0+1,this%pot%nv)
    integer :: ivv, iv1, iv2, ir
    _REAL_ :: k, ep0(maxep0), cep0(maxep0), err0

    !.....getting the long-range asymptotic function of H(k) for 3D-RISM
    call rism1d_closure_calcXvv_DT(this,hvvk_dT)
    do iv2=1,this%pot%nv
       do ir=2,maxep0+1
          k = (ir-1)*this%pot%dk
          ep0(ir-1) = k
          cep0(ir-1) = 0.d0
          do iv1=1,this%pot%nv
             cep0(ir-1) = cep0(ir-1) + this%pot%qv(iv1)*this%xvv_dT(ir,iv1,iv2)
          end do

          cep0(ir-1) = - 4.d0*pi*cep0(ir-1)/k**2
          cep0(ir-1) = cep0(ir-1)
       end do
       call  polynomialInterpolation_progressive (ep0,cep0,4, 0.d0,delhv0_dT(iv2), err0)
    end do
    delhv0_dT = delhv0_dT - rism1d_closure_getDelHvLimit(this,hvvk)
  end function rism1d_closure_getDelHvLimit_DT

  !! !! !! !

  !! If Xvv has not been calculated since our last solution, allocate
  !! memory and calculate it.
  !! @param[] this rism1d closure object..
  !! @param[] hvvk k-space total correlation function
  subroutine rism1d_closure_calcXvv(this, hvvk)
    implicit none
#include "../xblas/f77/blas_namedconstants.fh"
    type(rism1d_closure),intent(inout) :: this
    _REAL_, intent(in) :: hvvk(:,:)
    integer :: iv1, iv2, ivv, ir
    if (associated(this%xvv)) return
    this%xvv => safemem_realloc(this%xvv,this%pot%nr,this%pot%nv,this%pot%nv)
    this%xvv = 0
    !.................. getting Xvv(k)=Wvv(k)+DensityV*Hvv(k) ..................
    ivv = 0
    do iv2 = 1, this%pot%nv
       do iv1 = 1, iv2
          ivv = ivv + 1
!!$          call BLAS_DWAXPBY_X(this%pot%nr,1d0,this%pot%wvv(:,iv1,iv2),1,&
!!$               this%pot%densityv(iv1),hvvk(:,ivv),1,&
!!$               this%xvv(:,iv1,iv2),1,BLAS_PREC_EXTRA)
!!$          if (iv1 /= iv2) &
!!$               call BLAS_DWAXPBY_X(this%pot%nr,1d0,this%pot%wvv(:,iv2,iv1),1,&
!!$               this%pot%densityv(iv2),hvvk(:,ivv),1,&
!!$               this%xvv(:,iv2,iv1),1,BLAS_PREC_EXTRA)
          do ir = 1, this%pot%nr
             this%xvv(ir,iv1,iv2) = this%pot%wvv(ir,iv1,iv2) + this%pot%densityv(iv1) * hvvk(ir,ivv)
             if (iv1 /= iv2) &
                  this%xvv(ir,iv2,iv1) = this%pot%wvv(ir,iv2,iv1) + this%pot%densityv(iv2) * hvvk(ir,ivv)
          end do
       end do
    end do
  end subroutine rism1d_closure_calcXvv


  !! If Xvv_dT has not been calculated since our last solution, allocate memory and
  !! calculate it
  !! IN:
  !!    this    : rism1d closure object.
  !!    hvvk_dT : k-space temperature derivative total correlation function
  subroutine rism1d_closure_calcXvv_DT(this,hvvk_dT)
    implicit none
    type(rism1d_closure),intent(inout) :: this
    _REAL_, intent(in) :: hvvk_dT(:,:)
    integer :: iv1, iv2, ivv, ir
    if (associated(this%xvv_dT)) return
    this%xvv_dT => safemem_realloc(this%xvv_dT,this%pot%nr,this%pot%nv,this%pot%nv)
    this%xvv_dT=0
    !.................. getting Xvv(k)=Wvv(k)+DensityV*Hvv(k) ..................
    ivv = 0
    do iv2=1,this%pot%nv
       do iv1=1,iv2
          ivv = ivv + 1
          do ir=1,this%pot%nr
             this%xvv_dT(ir,iv1,iv2) = this%pot%densityv(iv1)*hvvk_dT(ir,ivv)
             if (iv1 /= iv2) &
                  this%xvv_dT(ir,iv2,iv1) = this%pot%densityv(iv2)*hvvk_dT(ir,ivv)
          end do
       end do
    end do
  end subroutine rism1d_closure_calcXvv_DT


  !! Returns an NV X NV array of total excess coordination numbers excluding
  !! multiplicity
  !! IN:
  !!    this : rism1d closure object.
  !!    hvv  : total correlation function

  function rism1d_closure_getExNumber(this,hvv) result(exvv)
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: hvv(:,:)
    _REAL_ :: exvv(this%pot%nv,this%pot%nv)
    integer :: ivv, iv1, iv2
    !......... excess number (excluding Mult) of sites 2 around 1 ..........
    ivv = 0
    do iv2=1,this%pot%nv
       do iv1=1,iv2
          ivv = ivv + 1
          exvv(iv1,iv2) = this%pot%densityv(iv2)/this%pot%mtv(iv2) * hvv(1,ivv)
          exvv(iv2,iv1) = this%pot%densityv(iv1)/this%pot%mtv(iv1) * hvv(1,ivv)
       end do
    end do
  end function rism1d_closure_getExNumber


  !! Returns a pointer to an NR X NVV array of structure factors.  This memory
  !! must be freed (preferably with safemem_dealloc) as it is not freed locally or
  !! after the object instance is destroyed.
  !! IN:
  !!    this : rism1d closure object.
  !!    hvv  : total correlation function

  function rism1d_closure_getStructFactor(this,hvv) result(svv)
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: hvv(:,:)
    _REAL_, pointer :: svv(:,:)
    _REAL_ :: densitytot
    integer :: ivv, iv1, iv2, ir
    nullify(svv)
    svv=>safemem_realloc(svv, this%pot%nr, this%pot%nvv)
    densitytot = sum(this%pot%densitysp)
    ivv = 0
    do iv2=1,this%pot%nv
       do iv1=1,iv2
          ivv = ivv + 1
          do ir=1,this%pot%nr
             svv(ir,ivv) = hvv(ir,ivv) &
                  * this%pot%densityv(iv1)/this%pot%mtv(iv1) * this%pot%densityv(iv2)/this%pot%mtv(iv2) &
                  / densitytot
          end do
       end do
    end do
  end function rism1d_closure_getStructFactor


  !! Returns a pointer to an NR X NVV array of the running site-site excess number.
  !! This is the excess number of a site within a given radius. The memory for
  !! this pointer must be freed (preferably with safemem_dealloc) as it is not
  !! freed locally or after the object instance is destroyed.
  !! IN:
  !!    this : rism1d closure object.
  !!    gvv  : pair distribution function

  function rism1d_closure_getRunExNumber(this,gvv) result(exnvv)
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:)
    _REAL_, pointer :: exnvv(:,:,:)
    nullify(exnvv)
    exnvv => exNumber(this,gvv,.true.)
  end function rism1d_closure_getRunExNumber


  !! Returns a pointer to an NR X NVV array of the running site-site number.
  !! This is the number of a site within a given radius. The memory for
  !! this pointer must be freed (preferably with safemem_dealloc) as it is not
  !! freed locally or after the object instance is destroyed.
  !! IN:
  !!    this : rism1d closure object.
  !!    gvv  : pair distribution function

  function rism1d_closure_getRunNumber(this,gvv) result(nvv)
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:)
    _REAL_, pointer :: nvv(:,:,:)
    nullify(nvv)
    nvv => exNumber(this,gvv,.false.)
  end function rism1d_closure_getRunNumber


  !! Calculates the pressure in [kT / A^3] of the system using the virial
  !! path.  To convert to Pacals, for example, multiply by
  !! 1.d27 * kb * temperature
  !! where kb is Boltzmann's constant [j/K] and temperature is in [K]
  !! IN:
  !!    this : rism1d closure object.
  !!    gvv  : pair distribution function

  function rism1d_closure_getPressureVirial(this,gvv) result(pressure)
    use constants, only :pi
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: gvv(:,:)
    _REAL_ :: pressure
    _REAL_ :: r
    integer :: ir, ivv, iv1, iv2, cnt

    !make sure the derivative of the potential is calculated
    call rism1d_potential_derivative(this%pot)

    !calculate the integral for each site
    pressure=0
    ivv=0
    do iv1 = 1, this%pot%nv
       do iv2 = 1, iv1
          ivv = ivv + 1
          do ir = 2, this%pot%nr
             r = (ir -1)*this%pot%dr
             !the full eqn. is a double sum over sites. so, if the
             !sites are different, we need to double count the
             !contribution.
             if (iv1 == iv2)  then
                cnt = 1
             else
                cnt = 2
             end if
             pressure = pressure + &
                  cnt*this%pot%densityv(iv1)*this%pot%densityv(iv2)*&
                  gvv(ir, ivv)*this%pot%duvv(ir,ivv)*r**3
          end do
       end do
    end do
    pressure = sum(this%pot%densitysp) - this%pot%dr*PI*2d0/3d0*pressure
  end function rism1d_closure_getPressureVirial


  !! Calculates the pressure in [kT / A^3] of the system using the free energy
  !! path.  To convert to Pacals, for example, multiply by
  !! 1.d27 * kb * temperature
  !! where kb is Boltzmann's constant [j/K] and temperature is in [K]
  !! IN:
  !!    this : rism1d closure object.
  !!    gvv  : pair distribution function
  !!    cvv  : direct correlation function
  !! OUT:
  !!     pressure from the free energy route if defined for the given
  !!     closure.  If not, returns HUGE

  function rism1d_closure_getPressureFE(this,gvv,cvv) result(pressure)
    use safemem
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:)
    _REAL_ :: pressure

    !r-space contribution. This is closure dependent
    call fe_press_r(this,gvv,cvv)

    !check that the closure supports the pressure calculation
    if (this%FE_press_r==HUGE(1d0)) then
       pressure=this%FE_press_r
       return
    end if

    !k-space contribution. This is not closure dependent
    call fe_press_k(this,cvv)

    !check that the k-space contribution was correctly calculated
    if (this%pressk == tiny(0d0)) then
       pressure = huge(1d0)
       return
    end if

    pressure = sum(this%pot%densitysp)+ this%pressk + this%FE_press_r
  end function rism1d_closure_getPressureFE


  !! Calculates the freeEnergy in kT
  !! IN:
  !!    this : rism1d object
  !!    gvv  : pair distribution function
  !!    cvv  : direct correlation function
  !! OUT:
  !!     free energy if defined for the given closure.  If not, returns HUGE

  function rism1d_closure_getFreeEnergy(this,gvv,cvv) result(fe)
    use safemem
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:)
    _REAL_ :: fe

    !r-space contribution. This is closure dependent
    call fe_press_r(this,gvv,cvv)

    !check that the closure supports the pressure calculation
    if (this%FE_press_r==HUGE(1d0)) then
       fe=this%FE_press_r
       return
    end if

    !k-space contribution. This is not closure dependent
    call fe_press_k(this,cvv)

    !check that the k-space contribution was correctly calculated
    if (this%fek == tiny(0d0)) then
       fe = 0
       return
    end if
    fe = (this%fek + this%FE_press_r)
  end function rism1d_closure_getFreeEnergy


  !! Calculates the partial molar volume of each species in A^3
  !! IN:
  !!    this : rism1d closure object.
  !!    cvv  : direct correlation function

  function rism1d_closure_getPMV(this,cvv) result(pmv)
    use constants, only : pi
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: cvv(:,:)
    _REAL_ :: pmv(this%pot%nsp)
    _REAL_ :: compressibility, r, cvvk0, cvvk0r
    integer :: ivv, iv1, iv2, iv21, iv22, ir, isp, iat
    compressibility = rism1d_closure_getCompressibility(this,cvv)
    iv22 = 0
    do isp=1,this%pot%nsp
       iv21 = iv22 + 1
       do iat=1,this%pot%nat(isp)
          iv22 = iv22 + 1
       end do
       !calculate Cvv(k=0)
       cvvk0 = 0.d0
       do ir=2,this%pot%nr
          r = (ir-1)*this%pot%dr
          cvvk0r = 0.d0
          do iv2=iv21,iv22
             do iv1=1,this%pot%nv
                ivv = this%pot%jvv(iv1,iv2)
                cvvk0r = cvvk0r + this%pot%densityv(iv1)*this%pot%mtv(iv2)*cvv(ir,ivv)
             end do
          end do
          cvvk0 = cvvk0 + r**2*cvvk0r
       end do
       cvvk0 = cvvk0 * 4.d0*pi*this%pot%dr
       pmv(isp) = compressibility*(1.d0-cvvk0)
    end do
  end function rism1d_closure_getPMV


  !! Calculates the excess chemical potential in kT for each site using
  !! the Pettitt-Rossky formula
  !! B. M. Pettitt; P. J. Rossky. J. Chem. Phys. 1986, 15, 5836-5844
  !! IN:
  !!    this : rism1d closure object.
  !!    gvv  : pair distribution function
  !!    cvv  : direct correlation function
  !! OUT:
  !!     excess chemical potential if defined for the given closure.  If not, returns HUGE

  function rism1d_closure_getExChem_PR(this,gvv,cvv) result(exchem)
    use safemem
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:)
    _REAL_ :: exchem(this%pot%nv)
    integer :: isp, iat, iv

    if (associated(this%KH)) then
       exchem = rism1d_kh_exChem_PR(this%kh,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else if (associated(this%PSEN)) then
       exchem = rism1d_psen_exChem_PR(this%psen,gvv,this%pot%uvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else if (associated(this%POLYT)) then
       !       exchem = rism1d_polyt_exChem(this%polyt,gvv,this%pot%uvv,cvv,this%pot%mtv,this%pot%jvv,&
       !            this%pot%densityv,densitytrgtv,this%pot%dr)
    else if (associated(this%HNC)) then
       exchem = rism1d_hnc_exChem_PR(this%hnc,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else if (associated(this%NUB)) then
       !       exchem = rism1d_nub_exChem(this%nub,gvv,cvv,this%pot%mtv,this%pot%jvv,&
       !            this%pot%densityv,densitytrgtv,this%pot%dr)
    else if (associated(this%DEVB)) then
       exchem = rism1d_devb_exChem_PR(this%devb,gvv,this%pot%uvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else
       exchem = HUGE(1d0)
    end if
  end function rism1d_closure_getExChem_PR


  !! Calculates the excess chemical potential in kT for each site using
  !! the Schmeer-Maurer formula
  !! G. Schmeer and A. Maurer. Phys. Chem. Chem. Phys., 2010, 12, 2407–2417
  !! IN:
  !!    this : rism1d closure object.
  !!    gvv  : pair distribution function
  !!    cvv  : direct correlation function
  !! OUT:
  !!     excess chemical potential if defined for the given closure.  If not, returns HUGE

  function rism1d_closure_getExChem_SM(this,gvv,cvv) result(exchem)
    use safemem
    use constants, only : pi
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:)
    _REAL_ :: exchem(this%pot%nv)
    _REAL_ :: exchemk(this%pot%nv)
    !wcvvk : wvv(k)*cvv(k) for a given k
    _REAL_ :: wcvvk(this%pot%nv,this%pot%nv)
    !ak : temp variable to hold the matrix to be inverted
    !bk : temp variable.  First hold wcvvk then hold ak^-1*wcvvk
    _REAL_ :: ak(this%pot%nv,this%pot%nv),bk(this%pot%nv,this%pot%nv)
    !k2 : k^2 for the spherical integral in k-space
    _REAL_ :: k2
    integer :: isp, iat1, iat2, iv, ir
    integer :: iv1, iv2, iv3, iu1, iu2
    integer ::  indx(this%pot%nv)
    integer :: err

    !get the closure dependent, r-space part
    if (associated(this%KH)) then
       exchem = rism1d_kh_exChem_r(this%kh,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else if (associated(this%PSEN)) then
       exchem = rism1d_psen_exChem_r(this%psen,gvv,this%pot%uvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else if (associated(this%HNC)) then
       exchem = rism1d_hnc_exChem_r(this%hnc,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else
       exchem = HUGE(1d0)
       return
    end if
    !get the k-space part that is common to HNC-like closures
    call calculate_cvvk(this,cvv)
    !this is the last two terms for Eqn. 21 of Schmeer-Maurer 2010.
    exchemk=0d0
    do ir = 2, this%pot%nr
       !get wvv(k)*cvv(k)
       wcvvk=0d0
       do iv2 = 1, this%pot%nv
          do iv1 = 1, this%pot%nv
             do iv3 = 1, this%pot%nv
                wcvvk(iv1,iv2) = wcvvk(iv1,iv2)+&
                     this%pot%wvv(ir,iv3,iv1)*this%cvvK(ir,this%pot%jvv(iv3,iv2))
             end do
          end do
       end do
       !place  wvv(k)*cvv(k)*density into ak
       do iv2 = 1, this%pot%nv
          do iv1 = 1, this%pot%nv
             ak(iv1,iv2) = -wcvvk(iv1,iv2)*this%pot%densityv(iv2)
          end do
       end do
       !get ak = 1-wvv(k)*cvv(k)*density
       do iv=1,this%pot%nv
          ak(iv,iv) = 1d0+ak(iv,iv)
       end do
       !set bk=wvv(k)*cvv(k)
       bk = wcvvk

       !Invert Ak using LU decomposition and multiply by Bk.  The result is in bk
       call DGESV(this%pot%nv,this%pot%nv,ak,this%pot%nv,indx,bk,this%pot%nv,err)
       if (err < 0) then
          err = err*(-1)
          call rism_report_error("Linear equation solver failed.")
       end if

       !sum the two terms into bk
       bk=wcvvk-bk

       !multiply by the derivative of the density w.r.t. the density
       !of the given site and take the trace.  This is the same as
       !assigning to site i b(i,i).  Also multiply by k^2 (spherical
       !integral)
       k2=((ir-1)*this%pot%dk)**2
       do iv=1,this%pot%nv
          exchemk(iv) = exchemk(iv) + bk(iv,iv)*this%pot%mtv(iv)*k2
       end do
    end do
    !Note that this term has a factor of 1/2 as well as the 4*PI from
    !spherical itegration and (2*PI)**-3 for the Fourier transform
    !normalization
    exchemk=exchemk*(2d0*PI)**(-2)*this%pot%dk
    exchem = exchem+exchemk
  end function rism1d_closure_getExChem_SM


  !! Calculates the excess chemical potential in kT for each site using
  !! the original Singer-Chandler formula
  !! IN:
  !!    this : rism1d closure object.
  !!    gvv  : pair distribution function
  !!    hvvk : k-space total correlation function
  !!    cvv  : direct correlation function
  !! OUT:
  !!     excess chemical potential if defined for the given closure.  If not, returns HUGE

  function rism1d_closure_getExChem_SC(this,gvv,hvvk,cvv) result(exchem)
    use safemem
    use constants, only : pi
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: gvv(:,:), hvvk(:,:), cvv(:,:)
    _REAL_ :: exchem(this%pot%nv)
    _REAL_ :: exchemk(this%pot%nv)
    integer :: isp, iat1, iat2, iv, ir
    integer :: iv1, iv2, iu1, iu2

    !get the closure dependent part
    if (associated(this%KH)) then
       exchem = rism1d_kh_exChem_r(this%kh,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else if (associated(this%PSEN)) then
       exchem = rism1d_psen_exChem_r(this%psen,gvv,this%pot%uvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else if (associated(this%HNC)) then
       exchem = rism1d_hnc_exChem_r(this%hnc,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else
       exchem = HUGE(1d0)
       return
    end if
    !get the k-space part that is common to HNC-like closures
    call rism1d_closure_calcXvv(this,hvvk)
    call calculate_cvvk(this,cvv)
    !this is the last term for Eqn. 2.8 of Singer-Chandler 1985.
    !Indexing variables iu* and iv* are for the solute and solvent
    !respectively.  Here we take all species to be part of the
    !solvent, even if identical to the solute. Note that this term in
    !2.8 leaves out a factor of omega that appears in Eqn. 2.7
    exchemk=0d0
    iu1=0
    !For the contribution from each solute site we need to iterate
    !over all of the other solute sites on that species.
    do isp=1,this%pot%nsp
       !iterate over each site in this species
       do iat1=1,this%pot%nat(isp)
          iu1=iu1+1
          !iterate over all of the sites in the species again.  First,
          !figure out the global site number for the first site in the
          !species
          if (isp>1) then
             iu2=sum(this%pot%nat(1:isp-1))
          else
             iu2=0
          end if
          do iat2=1,this%pot%nat(isp)
             iu2=iu2+1
             !iterate over the solvent
             do iv1=1,this%pot%nv
                do iv2=1,this%pot%nv
                   do ir=2,this%pot%nr
                      exchemk(iu1) = exchemk(iu1)&
                           - this%cvvK(ir,this%pot%jvv(iu1,iv1))&
                           * this%pot%mtv(iu2)&
                           * this%pot%wvv(ir,iu1,iu2)&
                           * this%pot%densityv(iv2)&
                           * this%xvv(ir,iv1,iv2)&
                           * this%cvvK(ir,this%pot%jvv(iv2,iu2))&
                           * ((ir-1)*this%pot%dk)**2
                   end do
                end do
             end do
          end do
       end do
    end do
    !Note that this term has a factor of 1/2 as well as the 4*PI from
    !spherical itegration and (2*PI)**-3 for the Fourier transform
    !normalization
    exchemk=exchemk*(2d0*PI)**(-2)*this%pot%dk
    exchem = exchem+exchemk
  end function rism1d_closure_getExChem_SC


  !! Calculates the ionic(?) excess chemical potential in kT for each site.
  !! IN:
  !!    this : rism1d closure object.
  !!    gvv  : pair distribution function
  !!    cvv  : direct correlation function
  !! OUT:
  !!     ionic excess chemical potential if defined for the given closure.  If not, returns HUGE

  function rism1d_closure_getExChemIon(this,gvv,cvv) result(exchem)
    use safemem
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:)
    _REAL_ :: exchem(this%pot%nv)
    integer :: isp, iat, iv

    if (associated(this%KH)) then
       exchem = rism1d_kh_exChemIon(this%kh,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else if (associated(this%PSEN)) then
       exchem = rism1d_psen_exChemIon(this%psen,gvv,this%pot%uvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else if (associated(this%POLYT)) then
       !       exchem = rism1d_polyt_exChemIon(this%polyt,gvv,this%pot%uvv,cvv,this%pot%mtv,this%pot%jvv,&
       !            this%pot%densityv,densitytrgtv,this%pot%dr)
    else if (associated(this%HNC)) then
       exchem = rism1d_hnc_exChemIon(this%hnc,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else
       exchem = HUGE(1d0)
    end if
  end function rism1d_closure_getExChemIon


  !! Calculates the solvation energy in kT for each site.
  !! IN:
  !!    this   : rism1d closure object.
  !!    gvv    : pair distribution function
  !!    uvv    : site-site potential
  !!    cvv    : direct correlation function
  !!    gvv_dT : temperature derivative pair distribution function
  !!    cvv_dT : temperature derivative direct correlation function
  !! OUT:
  !!     solvation energy if defined for the given closure.  If not, returns HUGE

  function rism1d_closure_getSolvene(this,gvv,uvv,cvv,gvv_dT,cvv_dT) result(solvene)
    use safemem
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:), uvv(:,:), cvv(:,:), gvv_dT(:,:), cvv_dT(:,:)
    _REAL_ :: solvene(this%pot%nv)
    integer :: isp, iat, iv

    if (associated(this%KH)) then
       solvene = rism1d_kh_Solvene(this%kh,gvv,cvv,gvv_dT,cvv_dT,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else if (associated(this%PSEN)) then
       solvene = rism1d_psen_Solvene(this%psen,gvv,uvv,cvv,gvv_dT,cvv_dT,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else if (associated(this%POLYT)) then
       !       solvene = rism1d_polyt_Solvene(this%polyt,gvv,uvv,cvv,gvv_dT,cvv_dT,this%pot%mtv,this%pot%jvv,&
       !            this%pot%densityv,densitytrgtv,this%pot%dr)
    else if (associated(this%HNC)) then
       solvene = rism1d_hnc_Solvene(this%hnc,gvv,cvv,gvv_dT,cvv_dT,this%pot%mtv,this%pot%jvv,&
            this%pot%densityv,this%pot%dr)
    else
       solvene = HUGE(1d0)
    end if
  end function rism1d_closure_getSolvene


  !! !! !! !! !! !! !! !! !! !! !! !! !PRIVATE!! !! !! !! !! !! !! !! !! !! !! !!


  !! Returns a pointer to an NR X NV X NV array of the running
  !! (excess) site-site number.  This is the (excess) number of a site
  !! within a given radius. The memory for this pointer must be freed
  !! (preferably with safemem_dealloc) as it is not freed locally or
  !! after the object instance is destroyed.
  !! IN:
  !!    this : rism1d object
  !!    gvv  : pair distribution function
  !!    excess : .true. for excess, .false. for total
  function exNumber(this,gvv,excess) result(exnvv)
    use constants, only : pi
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:)
    logical, intent(in) ::excess
    _REAL_, pointer :: exnvv(:,:,:)
    _REAL_ :: tempvv0, tempvv1
    integer :: ivv, iv1, iv2, ir
    nullify(exnvv)
    exnvv => safemem_realloc(exnvv, this%pot%nr, this%pot%nv, this%pot%nv)

    !the cummulative integral is done inplace using the trapazoidal rule
    !for either case of hvv or gvv.  This means using temporary variables to
    !hold values of hvv and gvv before they are overwritten
    ivv = 0
    do iv2=1,this%pot%nv
       do iv1=1,iv2
          ivv = ivv + 1
          !select hvv or gvv
          if (excess) then
             exnvv(:,iv1,iv2) = gvv(:,ivv)-1d0
          else
             exnvv(:,iv1,iv2) = gvv(:,ivv)
          end if
          !integrate
          tempvv0 = exnvv(1,iv1,iv2)
          exnvv(1,iv1,iv2) = 0.d0
          do ir=2,this%pot%nr
             tempvv1 = exnvv(ir,iv1,iv2)
             exnvv(ir,iv1,iv2) = exnvv(ir-1,iv1,iv2) &
                  + 2.d0*pi*this%pot%dr &
                  *(((ir-2)*this%pot%dr)**2*tempvv0 &
                  + ((ir-1)*this%pot%dr)**2*exnvv(ir,iv1,iv2) )
             tempvv0=tempvv1
          end do
       end do
    end do
    do iv2=1,this%pot%nv
       do iv1=1,iv2
          !          exnvv(:,iv2,iv1) = this%pot%densityv(iv1)/this%pot%mtv(iv1) * exnvv(:,iv1,iv2)
          exnvv(:,iv2,iv1) = this%pot%densityv(iv2) * exnvv(:,iv1,iv2)
          if (iv1 /= iv2)&
               !               exnvv(:,iv1,iv2) = this%pot%densityv(iv2)/this%pot%mtv(iv2) * exnvv(:,iv1,iv2)
               exnvv(:,iv1,iv2) = this%pot%densityv(iv1) * exnvv(:,iv1,iv2)
       end do
    end do
  end function exNumber


  !! Calculate and store the k-space representation of Cvv unless it
  !! already is stored
  !! IN:
  !!    this : rism1d closure object.
  !!    cvv  : direct correlation function
  !! SIDEEFFECT:
  !!    If this%Cvvk is NULL, memory is allocated and the k-space form
  !!    of Cvv is stored
  subroutine calculate_cvvk(this,cvv)
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: cvv(:,:)
    integer :: ivv, ir
    this%cvvk => safemem_realloc(this%cvvk, this%pot%nr, this%pot%nvv)
    !.............................. FFT>K Cvv ..............................
    do ivv=1,this%pot%nvv
       do ir=1,this%pot%nr
          this%cvvk(ir,ivv) = (ir-1)*this%pot%dr*cvv(ir,ivv)
       end do
    end do

    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          this%cvvk(ir,ivv) = this%cvvk(ir,ivv) + this%pot%ulrvv(ir,ivv)
       end do
    end do

    do ivv=1,this%pot%nvv
       call  sinfti (this%cvvk(2,ivv),this%pot%nr-1, this%pot%dr, +1)
    end do

    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          this%cvvk(ir,ivv) = this%cvvk(ir,ivv) - this%pot%ulkvv(ir,ivv)
       end do
    end do

    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          this%cvvk(ir,ivv) = this%cvvk(ir,ivv) / ((ir-1)*this%pot%dk)
       end do
    end do
  end subroutine calculate_cvvk


  !> Calculates the r-space contributions to the free energy and pressure
  !! along the free energy path.  This is done the first time the
  !! subroutine is called and the result is stored.
  !! @param[in,out] this rism1d closure object.
  !! @param[in] gvv Pair distribution function.
  !! @param[in] cvv Direct correlation function.
  !! @sideeffects
  !!    If this%FE_PRESS_r is not HUGE(1d0), the calculation is carried
  !!    out and the result stored in this%FE_press_r.
  subroutine fe_press_r(this,gvv,cvv)
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:)
    !r-space contribution. This is closure dependent
    if (this%FE_press_r /= HUGE(1d0)) then
       return
    end if
    if (associated(this%KH)) then
       this%FE_press_r = rism1d_kh_freeEnergyPressure_r(this%kh,gvv,cvv,&
            this%pot%mtv,this%pot%densityv,this%pot%dr)
    else if (associated(this%PSEN)) then
       this%FE_press_r = rism1d_psen_freeEnergyPressure_r(this%psen,gvv,&
            this%pot%uvv,cvv,this%pot%mtv,this%pot%densityv,this%pot%dr)
    else if (associated(this%HNC)) then
       this%FE_press_r = rism1d_hnc_freeEnergyPressure_r(this%hnc,gvv,cvv,&
            this%pot%mtv,this%pot%densityv,this%pot%dr)
    else if (associated(this%DEVB)) then
       this%FE_press_r = rism1d_devb_freeEnergyPressure_r(this%devb,gvv,&
            this%pot%uvv,cvv,this%pot%mtv,this%pot%densityv,this%pot%dr)
    end if
  end subroutine fe_press_r


  !> Calculates the k-space contributions to the free energy and pressure
  !! along the free energy path.  This is done the first time the
  !! subroutine is called and the result is stored.
  !! @param[in,out] this rism1d closure object.
  !! @param[in] cvv Direct correlation function.
  !! @sideeffects
  !! If this%fek and this%pressk are not HUGE(1d0), the calculation is
  !! carried out and the result stored in this%fek and this%pressk.
  subroutine fe_press_k(this,cvv)
    use constants, only : pi
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: cvv(:,:)
    _REAL_ :: fe_k, pr_k, det, wcvvkp(this%pot%nv,this%pot%nv), ak(this%pot%nv,this%pot%nv)
    integer :: indx(this%pot%nv), err, ir,iv, iv1, iv2, iv3
    _REAL_ :: k2
    if (this%fek /= HUGE(1d0)) then
       return
    end if
    call calculate_cvvk(this,cvv)
    fe_k=0d0
    pr_k=0d0
    do ir = 2, this%pot%nr
       !k^2 for spherical integration
       k2=((ir-1)*this%pot%dk)**2

       !get wvv(k)*cvv(k)
       wcvvkp=0d0
       do iv2 = 1, this%pot%nv
          do iv1 = 1, this%pot%nv
             do iv3 = 1, this%pot%nv
                wcvvkp(iv1,iv2) = wcvvkp(iv1,iv2)+&
                     this%pot%wvv(ir,iv3,iv1)*this%cvvK(ir,this%pot%jvv(iv3,iv2))
             end do
             wcvvkp(iv1,iv2) = wcvvkp(iv1,iv2)*this%pot%densityv(iv1)
          end do
       end do

       !Add Tr[wvv(k)*cvv(k)*p] to the free energy now.  wcvvkp is
       !modified below
       do iv1=1,this%pot%nv
          fe_k = fe_k+wcvvkp(iv1,iv1)*k2
       end do

       !place  wvv(k)*cvv(k)*density into ak
       do iv2 = 1, this%pot%nv
          do iv1 = 1, this%pot%nv
             ak(iv1,iv2) = -wcvvkp(iv1,iv2)
          end do
       end do

       !get ak = 1-wvv(k)*cvv(k)*density
       do iv=1,this%pot%nv
          ak(iv,iv) = 1d0+ak(iv,iv)
       end do

       ! LU decompostion of ak to calculate determinant
       call DGETRF(this%pot%nv,this%pot%nv,ak,this%pot%nv,indx,err)
       !DGETRF does not calculate the determinant of the permumation
       !maxtrix (as the Numerical Recipes subroutine does) so we
       !calculate it here IPIV is 'compressed' permutation matrix
       !where the index of the exchange row is record for each row
       !index.  So we count the number of rows that are not
       !exchanged with themselves and divide by two to get the total
       !number of exchanges.  The detminant is then
       !(-1)^number_of_exchanges
       det =1d0
       do iv1=1,this%pot%nv
          if (indx(iv1) /= iv1) then
             det = -det* ak(iv1,iv1)
          else
             det = det * ak(iv1,iv1)
          end if
       end do
       if (det <= 0.d0)  then
          call rism_report_error("fe__k: non-positive Det[1-Wvv*Cvv*Density]")
       end if

       !use the LU decomposition to invert the matrix and multiply it
       !by wcvvkp.  The result is stored in wcvvkp
       !............ calculating (1-Wvv*Cvv*Density)^(-1)*Wvv*Cvv*Density .............
       !        call DGETRS(TRANS, N, NRHS,A,LDA,IPIV,B,LDB,INFO)
       call DGETRS('N', this%pot%nv, this%pot%nv,ak,this%pot%nv,indx,wcvvkp,this%pot%nv,err)
       if (err < 0) then
          err = err*(-1)
          call rism_report_error("fe_press_k: Linear equation solver failed.")
       end if

       ! add in the remaining contributions to free energy and pressure
       fe_k = fe_k + log(det)*k2
       pr_k = pr_k + log(det)*k2
       do iv=1,this%pot%nv
          pr_k = pr_k + wcvvkp(iv,iv)*k2
       end do
    end do
    this%fek = fe_k*(2d0*pi)**(-2)*this%pot%dk
    this%pressk = -pr_k*(2d0*pi)**(-2)*this%pot%dk
  end subroutine fe_press_k

end module rism1d_closure_c
