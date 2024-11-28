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

!> Partial series expansion of order-n (PSE-n) closure
!!
!! Stefan M. Kast and Thomas Kloss. J. Chem. Phys. 129, 236101 (2008)
!!
!! Interpolates between the Kovalenko-Hirata (KH) and hypernetted chain equation
!! closures (HNC).  Order-1 gives KH and order-infinity gives HNC.
!!
!! This class provides PSE-n specific contributions for thermodynamic
!! relations.  This includes the radial distribution function, gvv;
!! its temperature derivative, gvv_dT; excess chemical potential; free
!! energy; and pressure (free energy route).
!!
!! Excess Chemical Potential
!! -------------------------
!! Under XRISM there are three equivalent formulations for the excess
!! chemical potential for the hypernetted-chain equation closure:
!!
!! S. J. Singer; D. Chandler. Molecular Physics, 1985, 55, 3, 621—625
!!
!! B. M. Pettitt; P. J. Rossky. J. Chem. Phys. 1986, 15, 5836-5844
!!
!! G. Schmeer and A. Maurer. Phys. Chem. Chem. Phys., 2010, 12, 2407–2417
!!
!! PSE-n has a similar form where the r-space contribution is closure
!! dependent and the k-space is not.  The Pettitt-Rossky equivalent
!! can be found in:
!!
!! Stefan M. Kast and Thomas Kloss. J. Chem. Phys. 129, 236101 (2008)
!!
!! For the total excess chemical potential, we use the simplest
!! formulation, Pettitt-Rossky/Kast-Kloss.  All three
!! forumlations have the same r-space form and this is calculated
!! separately in rism1d_psen_exchem_r().
!!
!! Free energy and pressure
!! ------------------------
!! Only the closure dependent, r-space component of the free energy
!! and pressure are provided here.  This contribution is the same for
!! both quantities and is calculated with a single function,
!! rism1d_psen_freeEnergyPressure_r().
module rism1d_psen_c
  use rism_report_c
  use safemem
  
  ! The PSE-n type.
  type rism1d_psen
     !> Order of the series expansion, >= 1.
     integer :: order = 0
     !> order + 1.
     integer :: order1
     ! (order+1)!
     _REAL_ :: order1fac
  end type rism1d_psen

  public rism1d_psen_new, rism1d_psen_destroy!, rism1d_psen_gvv
  
contains

  !> Initializes the PSEN closure
  !! @param[in,out] this PSEN object.
  subroutine rism1d_psen_new(this, order)
    implicit none
    type(rism1d_psen), intent(inout) :: this
    integer, intent(in) :: order
    integer :: i
    double precision :: orderfac

    if (order < 1) then
       call rism_report_error("PSE-n closure: order must be > 0.")
    end if
    this%order = order
    this%order1 = order+1
    this%order1fac = 1
    orderfac = 1

    ! Calculate the factorial of (order+1) for future reference.
    ! Check that we are not losing precision or overflowing the
    ! precision.
    do i = 2, this%order1
       orderfac = this%order1fac
       this%order1fac = this%order1fac * i
       ! Keep 16 significant digits.
       if (abs((this%order1fac / i) / orderfac - 1d0) > 1e-16) then
          call rism_report_error("(a,i4)",&
               "PSE-n: numerical precision exceeded.  Use a smaller order. Try ",i-1)
       end if
    end do
  end subroutine rism1d_psen_new


  !> Calculates Gvv from Uvv, Hvv, and Cvv using the PSEN closure.
  !! @param[] this The PSEN closure object.
  !! @param[] gvv Site-site pair correlation function.
  !! @param[] uvv Site-site potential.
  !! @param[] hvv Site-site total correlation function.
  !! @param[] cvv Site-site direct correlation function.
  subroutine rism1d_psen_gvv(this,gvv, uvv, hvv, cvv)
    implicit none
    type(rism1d_psen), intent(in) :: this
    _REAL_, intent(out) :: gvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),hvv(:,:),cvv(:,:)
    integer :: ivv, ir, i
    _REAL_ :: tvv, orderfac
    do ivv = 1, ubound(gvv,2)
       do ir = 1, ubound(gvv,1)
          tvv = -uvv(ir,ivv) +hvv(ir,ivv)-cvv(ir,ivv)
          if (tvv>=0d0)then
             gvv(ir,ivv) = 1d0
             orderfac = 1d0
             do i=1,this%order
                orderfac = orderfac*i
                gvv(ir,ivv) = gvv(ir,ivv) + (tvv**i)/orderfac
             end do
          else
             gvv(ir,ivv) = exp(tvv)
          end if
       end do
    end do
  end subroutine rism1d_psen_gvv


  !! Calculates temperature derivative Gvv from Gvv, Uvv, Hvv_dT, and Cvv_dT using
  !! the associated closure
  !! IN:
  !!    this   : the closure object
  !!    gvv_dT  : site-site temperature derivative pair correlation function
  !!    gvv    : site-site pair correlation function
  !!    cvv  : site-site direct correlation function
  !!    hvv_dT  : site-site temperature derivative total correlation function
  !!    cvv_dT  : site-site temperature derivative direct correlation function
  subroutine rism1d_psen_gvv_dT(this, gvv_dT, uvv, gvv, cvv, hvv_dT, cvv_dT)
    implicit none
    type(rism1d_psen), intent(in) :: this
    _REAL_, intent(out) :: gvv_dT(:,:)
    _REAL_, intent(in) :: uvv(:,:),gvv(:,:),cvv(:,:),hvv_dT(:,:),cvv_dT(:,:)
    integer :: ivv, ir, i
    _REAL_ :: tvv_dT, tvv, orderfac
    do ivv = 1, ubound(gvv_dT,2)
       do ir = 1, ubound(gvv_dT,1)
          tvv_dT = uvv(ir,ivv) +hvv_dT(ir,ivv)-cvv_dT(ir,ivv)
          tvv = -uvv(ir,ivv) +gvv(ir,ivv)-1d0-cvv(ir,ivv)
          if (gvv(ir,ivv)>=1d0)then
             gvv_dT(ir,ivv) = 1d0
             orderfac = 1d0
             do i=1,this%order-1
                orderfac = orderfac*i
                gvv_dT(ir,ivv) = gvv_dT(ir,ivv) + (tvv**i)/orderfac
             end do
             gvv_dT(ir,ivv) = gvv_dT(ir,ivv)*tvv_dT
          else
             gvv_dT(ir,ivv) = tvv_dT*gvv(ir,ivv)
          end if
       end do
    end do
  end subroutine rism1d_psen_gvv_dT



  !! Calculates Bvv (bridge function) from Uvv, Gvv, and Cvv using the PSEN closure
  !! IN:
  !!    this : the PSEN closure object
  !!    uvv  : site-site potential
  !!    gvv  : site-site total correlation function
  !!    cvv  : site-site direct correlation function
  !! OUT:
  !!    bvv  : site-site bridge function
  function rism1d_psen_bvv(this, uvv, gvv, cvv) result(bvv)
    implicit none
    type(rism1d_psen), intent(in) :: this
    _REAL_, pointer :: bvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),gvv(:,:),cvv(:,:)
    integer :: ivv, ir, i
    _REAL_ :: tvv, orderfac
    nullify(bvv)
    bvv =>safemem_realloc(bvv,ubound(gvv,1),ubound(gvv,2))
    do ivv = 1, ubound(bvv,2)
       do ir = 1, ubound(bvv,1)
          tvv = -uvv(ir,ivv) +gvv(ir,ivv)-1d0-cvv(ir,ivv)
          if (tvv>=0d0)then
             orderfac = 1d0
             bvv(ir,ivv) = 1d0
             do i=1,this%order
                orderfac = orderfac*i
                bvv(ir,ivv) = bvv(ir,ivv) + (tvv**i)/orderfac
             end do
             bvv(ir,ivv) = log(bvv(ir,ivv))
             bvv(ir,ivv) = bvv(ir,ivv)-tvv
          else
             bvv(ir,ivv) = 0
          end if
       end do
    end do
  end function rism1d_psen_bvv


  !! Calculate the excess chemical potential in kT using the
  !! Pettitt-Rossky formula
  !! IN:
  !!    this : the closure object
  !!    gvv  : site-site pair correlation function
  !!    uvv  : site-site potential
  !!    cvv  : site-site direct correlation function
  !!    mtv  : multiplicity of each site
  !!    densityv : number density of each site
  !!    dr   : r-space grid spacing
  function rism1d_psen_exchem_PR(this, gvv, uvv, cvv, mtv, jvv, densityv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_psen), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),uvv(:,:),cvv(:,:), densityv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))
    !hnc    : HNC contribution
    !tsvv   : t^* = u +h -c
    !series : contribution from the series expansion
    _REAL_ :: r, hnc,hvv,tsvv, exchemv,series
    integer :: ivv, iv1, iv2, ir,i
    do iv2=1,ubound(mtv,1) !nv
       exchem(iv2) = 0.d0
       do ir=2,ubound(gvv,1) !nr
          r = (ir-1)*dr
          exchemv = 0.d0
          do iv1=1,ubound(mtv,1) !nv
             ivv = jvv(iv1,iv2)
             hvv = gvv(ir,ivv) - 1.d0
             tsvv = -uvv(ir,ivv) + hvv -cvv(ir,ivv)
             if (tsvv >= 0.d0)  then
                series = tsvv**(this%order1)/this%order1fac
             else
                series = 0
             endif
             hnc = 0.5d0*hvv*hvv - cvv(ir,ivv)  - 0.5d0*hvv*cvv(ir,ivv)
             exchemv = exchemv + densityv(iv1)*mtv(iv2) &
                  *(hnc-series)
          enddo
          exchem(iv2) = exchem(iv2) + exchemv*r**2
       enddo
    enddo
    exchem = exchem * 4.d0*pi*dr
  end function rism1d_psen_exchem_PR


  !! Calculates the closure dependent, r-space component of excess chemical
  !! potential in kT.
  !!
  !!   \int 1/2*h^2*\theta(-h) - c dr
  !!
  !! This is the r-space part of the PSE-n version of the
  !! Singer-Chandler and Schmeer-Maurer forms and is the Pettitt-Rossky
  !! form without the 'hc/2' term.
  !! IN:
  !!    this : the closure object
  !!    gvv  : site-site pair correlation function
  !!    uvv  : site-site potential
  !!    cvv  : site-site direct correlation function
  !!    mtv  : multiplicity of each site
  !!    densityv : number density of each site
  !!    dr   : r-space grid spacing
  function rism1d_psen_exchem_r(this, gvv, uvv, cvv, mtv, jvv, densityv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_psen), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),uvv(:,:),cvv(:,:), densityv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))
    !hnc    : HNC contribution
    !tsvv   : t^* = u +h -c
    !series : contribution from the series expansion
    _REAL_ :: r, hnc,hvv,tsvv, exchemv,series
    integer :: ivv, iv1, iv2, ir,i
    do iv2=1,ubound(mtv,1) !nv
       exchem(iv2) = 0.d0
       do ir=2,ubound(gvv,1) !nr
          r = (ir-1)*dr
          exchemv = 0.d0
          do iv1=1,ubound(mtv,1) !nv
             ivv = jvv(iv1,iv2)
             hvv = gvv(ir,ivv) - 1.d0
             tsvv = -uvv(ir,ivv) + hvv -cvv(ir,ivv)
             if (tsvv >= 0.d0)  then
                series = tsvv**(this%order1)/this%order1fac
             else
                series = 0
             endif
             hnc = 0.5d0*hvv*hvv - cvv(ir,ivv)
             exchemv = exchemv + densityv(iv1)*mtv(iv2) &
                  *(hnc-series)
          enddo
          exchem(iv2) = exchem(iv2) + exchemv*r**2
       enddo
    enddo
    exchem = exchem * 4.d0*pi*dr
  end function rism1d_psen_exchem_r


  !! Calculate the ionic(?) excess chemical potential in kT.  This seems to take
  !! the contribution only from the last solvent type.  I'm not sure what the point is
  !! here.
  !! IN:
  !!    this : the closure object
  !!    gvv  : site-site pair correlation function
  !!    cvv  : site-site direct correlation function
  !!    mtv  : multiplicity of each site
  !!    densityv : number density of each site
  !!    dr   : r-space grid spacing
  function rism1d_psen_exchemIon(this, gvv, uvv, cvv, mtv, jvv, densityv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_psen), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),uvv(:,:),cvv(:,:), densityv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))
    !hnc    : HNC contribution
    !tsvv   : t^* = u +h -c
    !series : contribution from the series expansion
    _REAL_ :: r, hnc,hvv, tsvv, exchemv,series
    integer :: ivv, iv1, iv2, ir, i

    do iv2=1,ubound(mtv,1) !nv
       exchem(iv2) = 0.d0
       do ir=2,ubound(gvv,1) !nr
          r = (ir-1)*dr
          exchemv = 0.d0
          !          do iv1=1,ubound(mtv,1) !nv
          iv1=ubound(mtv,1)
          ivv = jvv(iv1,iv2)
          hvv = gvv(ir,ivv) - 1.d0
          tsvv = -uvv(ir,ivv) + hvv -cvv(ir,ivv)
          if (tsvv >= 0.d0)  then
             series = tsvv**(this%order1)/this%order1fac
          else
             series = 0
          endif
          hnc = 0.5d0*hvv*hvv - cvv(ir,ivv)  - 0.5d0*hvv*cvv(ir,ivv)
          exchemv = exchemv + densityv(iv2)*mtv(iv1) &
               *(hnc - series)
          !         enddo
          exchem(iv2) = exchem(iv2) + exchemv*r**2
       enddo
    enddo
    exchem = exchem * 4.d0*pi*dr
  end function rism1d_psen_exchemIon


  !! Calculate the solvation energy in kT
  !! IN:
  !!    this : the closure object
  !!    gvv  : site-site pair correlation function
  !!    uvv  : site-iste potential
  !!    cvv  : site-site direct correlation function
  !!    gvv_dT  : temperature derivative of site-site pair correlation function
  !!    cvv_dT  : temperature derivative of site-site direct correlation function
  !!    mtv  : multiplicity of each site
  !!    densityv : number density of each site
  !!    dr   : r-space grid spacing
  function rism1d_psen_solvene(this, gvv, uvv, cvv, gvv_dT, cvv_dT, mtv, jvv, densityv, dr) result(solvene)
    use constants, only : pi
    implicit none
    type(rism1d_psen), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),uvv(:,:), cvv(:,:), gvv_dT(:,:), cvv_dT(:,:), &
         densityv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: solvene(ubound(mtv,1))
    _REAL_ :: r, h2c,hvv, solvenev, hvv_dT
    integer :: ivv, iv1, iv2, ir
    _REAL_ :: orderfac, tsvv, tvv_dT, series
    orderfac = this%order1fac/this%order1

    do iv2=1,ubound(mtv,1) !nv
       solvene(iv2) = 0.d0
       do ir=2,ubound(gvv,1) !nr
          r = (ir-1)*dr
          solvenev = 0.d0
          do iv1=1,ubound(mtv,1) !nv
             ivv = jvv(iv1,iv2)
             hvv = gvv(ir,ivv) - 1.d0
             hvv_dT = gvv_dT(ir,ivv)
             if (hvv >= 0.d0)  then
                tsvv = -uvv(ir,ivv) + hvv -cvv(ir,ivv)
                tvv_dT = uvv(ir,ivv) + hvv_dT - cvv_dT(ir,ivv)
                series = tsvv**(this%order)/orderfac*tvv_dT
             else
                series = 0
             end if
             h2c = hvv*hvv_dT - cvv_dT(ir,ivv)
             solvenev = solvenev + densityv(iv1)*mtv(iv2) &
                  *(h2c - 0.5d0*hvv_dT*cvv(ir,ivv) - 0.5d0*hvv*cvv_dT(ir,ivv) -series)
          enddo
          solvene(iv2) = solvene(iv2) + solvenev*r**2
       enddo
    enddo
    solvene = -1.d0 * solvene * 4.d0*pi*dr
  end function rism1d_psen_solvene


  !! Calculate the r-space contribution to the pressure in internal units
  !! using the free energy route
  !! IN:
  !!    this : the closure object
  !!    gvv  : site-site pair correlation function
  !!    cvv  : site-site direct correlation function
  !!    mtv  : multiplicity of each site
  !!    densityv : number density of each site
  !!    dr   : r-space grid spacing
  function rism1d_psen_pressureFE(this, gvv, uvv, cvv, mtv, densityv, dr) result(pr)
    use constants, only : pi
    implicit none
    type(rism1d_psen), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),uvv(:,:),cvv(:,:),densityv(:),dr
    integer, intent(in) :: mtv(:)
    _REAL_ :: pr
    _REAL_ :: prv
    !hnc    : h^2*c/2 contribution
    !tsvv   : t^* = u +h -c
    !series : contribution from the series expansion
    _REAL_ :: r, hvv, h2c, tsvv, series
    integer :: ir,ivv,iv1,iv2, cnt, i
    !r-space
    pr = 0.d0
    do ir=2,ubound(gvv,1) !nr
       r = (ir-1)*dr
       prv = 0.d0
       ivv = 0
       do iv2=1,ubound(mtv,1) !nv
          do iv1=1,iv2
             ivv = ivv + 1
             !if the sites are different, we need to double count the contribution
             if (iv1 == iv2)  then
                cnt = 1
             else
                cnt = 2
             endif
             hvv = gvv(ir,ivv) - 1.d0
             tsvv = -uvv(ir,ivv) + hvv -cvv(ir,ivv)
             if (tsvv >= 0.d0)  then
                series = tsvv**(this%order1)/this%order1fac
             else
                series = 0
             endif
             h2c = 0.5d0*hvv*hvv - cvv(ir,ivv)
             prv = prv + cnt*mtv(iv1)*densityv(iv1)*densityv(iv2)*(h2c-series)
          enddo
       enddo
       pr = pr + prv*r**2
    enddo
    pr = pr * 2.d0*pi*dr
  end function rism1d_psen_pressureFE


  !! Calculate the r-space contribution to the free energy and pressure
  !! (free energy route) in internal units.  This is the same quanitity
  !! for both values.
  !! IN:
  !!    this : the closure object
  !!    gvv  : site-site pair correlation function
  !!    cvv  : site-site direct correlation function
  !!    mtv  : multiplicity of each site
  !!    densityv : number density of each site
  !!    dr   : r-space grid spacing
  function rism1d_psen_freeEnergyPressure_r(this, gvv, uvv, cvv, mtv, densityv, dr) result(fe)
    use constants, only : pi
    implicit none
    type(rism1d_psen), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),uvv(:,:),cvv(:,:),densityv(:),dr
    integer, intent(in) :: mtv(:)
    _REAL_ :: fe
    _REAL_ :: fev
    !h2c    : h^2*c/2 contribution
    !tsvv   : t^* = u +h -c
    !series : contribution from the series expansion
    _REAL_ :: r, hvv, h2c, tsvv, series
    integer :: ir,ivv,iv1,iv2, cnt
    !r-space
    fe = 0.d0
    do ir=2,ubound(gvv,1) !nr
       r = (ir-1)*dr
       fev = 0.d0
       ivv = 0
       do iv2=1,ubound(mtv,1) !nv
          do iv1=1,iv2
             ivv = ivv + 1
             !if the sites are different, we need to double count the contribution
             if (iv1 == iv2)  then
                cnt = 1
             else
                cnt = 2
             endif
             hvv = gvv(ir,ivv) - 1.d0
             tsvv = -uvv(ir,ivv) + hvv -cvv(ir,ivv)
             if (tsvv >= 0.d0)  then
                series = tsvv**(this%order1)/this%order1fac
             else
                series = 0
             endif
             h2c = 0.5d0*hvv*hvv - cvv(ir,ivv)
             fev = fev + cnt*densityv(iv1)*densityv(iv2)*(h2c-series)
          enddo
       enddo
       fe = fe + fev*r**2
    enddo
    fe = fe * 2.d0*pi*dr
  end function rism1d_psen_freeEnergyPressure_r


  !> Frees memory and resets the PSEN closure.
  !! @param[in,out] this PSEN object.
  subroutine rism1d_psen_destroy(this)
    implicit none
    type(rism1d_psen), intent(inout) :: this
    this%order = 0
  end subroutine rism1d_psen_destroy
  
end module rism1d_psen_c
