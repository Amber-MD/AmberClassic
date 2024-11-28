! <compile=optimized>

! The 3D-RISM-KH software found here is copyright (c) 2011-2012 by
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

!> Hypernetted-chain (HNC) equation closure class for 1D-RISM.  This class
!! provides HNC specific contributions for thermodynamic relations.
!! This includes the radial distribution function, gvv; its
!! temperature derivative, gvv_dT; excess chemical potential; free
!! energy; and pressure (free energy route).
!!
!! Excess Chemical Potential
!! -------------------------
!! Under XRISM there are three equivalent formulations for the excess
!! chemical potential:
!!
!! S. J. Singer; D. Chandler. Molecular Physics, 1985, 55, 3, 621—625
!!
!! B. M. Pettitt; P. J. Rossky. J. Chem. Phys. 1986, 15, 5836-5844
!!
!! G. Schmeer and A. Maurer. Phys. Chem. Chem. Phys., 2010, 12, 2407–2417
!!
!! For the total excess chemical potential, we use the simplest
!! formulation, Pettitt-Rossky.  All three forumlations have the same
!! R-space form and this is calculated separately in
!! rism1d_hnc_exchem_r().
!!
!! Free energy and pressure
!! ------------------------
!! Only the closure dependent, R-space component of the free energy
!! and pressure are provided here.  This contribution is the same for
!! both quantities and is calculated with a single function,
!! rism1d_hnc_freeEnergyPressure_r().
module rism1d_hnc_c
  use safemem
  
  !> The HNC type is actually empty.
  type rism1d_hnc
     logical hnc
  end type rism1d_hnc

  public rism1d_hnc_new, rism1d_hnc_destroy, rism1d_hnc_gvv, rism1d_hnc_exchem,&
       rism1d_hnc_exchemIon, rism1d_hnc_pressureFE, rism1d_hnc_freeEnergy, &
       rism1d_hnc_gvv_dT, rism1d_hnc_solvene

contains

  !> Initializes the HNC closure.
  !! @param[in,out] this HNC object.
  subroutine rism1d_hnc_new(this)
    implicit none
    type(rism1d_hnc), intent(inout) :: this
  end subroutine rism1d_hnc_new


  !> Calculates Gvv from Uvv, Hvv, and Cvv using the HNC closure.
  !! @param[in] this The HNC closure object.
  !! @param[out] gvv Site-site pair correlation function.
  !! @param[in] uvv Site-site potential.
  !! @param[in] hvv Site-site total correlation function.
  !! @param[in] cvv Site-site direct correlation function.
  subroutine rism1d_hnc_gvv(this,gvv, uvv, hvv, cvv)
    implicit none
    type(rism1d_hnc), intent(in) :: this
    _REAL_, intent(out) :: gvv(:,:)
    _REAL_, intent(in) :: uvv(:,:), hvv(:,:), cvv(:,:)
    integer :: ivv, ir
    _REAL_ :: tvv
    
    do ivv = 1, ubound(gvv, 2)
       do ir = 1, ubound(gvv, 1)
          tvv = -uvv(ir, ivv) + hvv(ir, ivv) - cvv(ir, ivv)
          gvv(ir, ivv) = exp(tvv)
       end do
    end do
  end subroutine rism1d_hnc_gvv


  !> Calculates temperature derivative Gvv from Gvv, Uvv, Hvv_dT, and Cvv_dT using
  !! the associated closure
  !! @param[in] this The closure object.
  !! @param[out] gvv_dT Site-site temperature derivative pair correlation function.
  !! @param[in] gvv Site-site pair correlation function.
  !! @param[in] hvv_dT Site-site temperature derivative total correlation function.
  !! @param[in] cvv_dT Site-site temperature derivative direct correlation function.
  subroutine rism1d_hnc_gvv_dT(this, gvv_dT, uvv, gvv, hvv_dT, cvv_dT)
    implicit none
    type(rism1d_hnc), intent(in) :: this
    _REAL_, intent(out) :: gvv_dT(:,:)
    _REAL_, intent(in) :: uvv(:,:),gvv(:,:),hvv_dT(:,:),cvv_dT(:,:)
    integer :: ivv, ir
    _REAL_ :: tvv_dT

    do ivv = 1, ubound(gvv_dT,2)
       do ir = 1, ubound(gvv_dT,1)
          tvv_dT = uvv(ir,ivv) +hvv_dT(ir,ivv)-cvv_dT(ir,ivv)
          gvv_dT(ir,ivv) = tvv_dT*gvv(ir,ivv)
       end do
    end do
  end subroutine rism1d_hnc_gvv_dT


  !> Calculates Bvv (bridge function) from Uvv, Hvv, and Cvv using the HNC closure
  !! @param[] this The HNC closure object.
  !! @param[] uvv Site-site potential.
  !! @param[] hvv Site-site total correlation function.
  !! @param[] cvv Site-site direct correlation function.
  !! @return Site-site bridge function.
  function rism1d_hnc_bvv(this, uvv, hvv, cvv) result(bvv)
    implicit none
    type(rism1d_hnc), intent(in) :: this
    _REAL_, pointer :: bvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),hvv(:,:),cvv(:,:)

    nullify(bvv)
    bvv => safemem_realloc(bvv, ubound(hvv,1), ubound(hvv,2))
    bvv=0
  end function rism1d_hnc_bvv


  !> Calculate the excess chemical potential in kT using the
  !! Pettitt-Rossky formula.
  !! @param[] this The closure object.
  !! @param[] gvv Site-site pair correlation function.
  !! @param[] cvv Site-site direct correlation function.
  !! @param[] mtv Multiplicity of each site.
  !! @param[] densityv Number density of each site.
  !! @param[] dr R-space grid spacing.
  !! @return Excess chemical potential in kT.
  function rism1d_hnc_exchem_PR(this, gvv, cvv, mtv, jvv, densityv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_hnc), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:), densityv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))
    _REAL_ :: r, h2c,hvv, exchemv!, temp(ubound(mtv,1)), tempv
    integer :: ivv, iv1, iv2, ir

    do iv2=1,ubound(mtv,1) !nv
       exchem(iv2) = 0.d0
       do ir=2,ubound(gvv,1) !nr
          r = (ir-1)*dr
          exchemv = 0.d0
          !          tempv=0d0
          do iv1=1,ubound(mtv,1) !nv
             ivv = jvv(iv1,iv2)
             hvv = gvv(ir,ivv) - 1.d0
             h2c = 0.5d0*hvv*hvv - cvv(ir,ivv)
             exchemv = exchemv + densityv(iv1)*mtv(iv2) &
                  *(h2c - 0.5d0*hvv*cvv(ir,ivv))
             !             tempv = tempv  -densitytrgtv(iv1)*mtv(iv2)* 0.5d0*hvv*cvv(ir,ivv)
          enddo
          exchem(iv2) = exchem(iv2) + exchemv*r**2
          !          temp(iv2) = temp(iv2) + tempv*r**2
       enddo
    enddo
    !    write(0,*) "HC", temp*4d0*pi*dr
    exchem = exchem * 4.d0*pi*dr
  end function rism1d_hnc_exchem_PR


  !> Calculates the closure dependent, R-space component of excess chemical
  !! potential in kT.
  !!
  !!   \int 1/2*h^2 - c dr
  !!
  !! This is the R-space part of the Singer-Chandler and Schmeer-Maurer
  !! forms and is the Pettitt-Rossky form without the 'hc/2' term.
  !! @param this The closure object.
  !! @param[] gvv Site-site pair correlation function.
  !! @param[] cvv Site-site direct correlation function.
  !! @param[] mtv Multiplicity of each site.
  !! @param[] densityv Number density of each site.
  !! @param[] dr R-space grid spacing.
  !! @return
  function rism1d_hnc_exchem_r(this, gvv, cvv, mtv, jvv, densityv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_hnc), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:), densityv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))
    _REAL_ :: r, h2c,hvv, exchemv
    integer :: ivv, iv1, iv2, iv3, ir

    do iv2=1,ubound(mtv,1) !nv
       exchem(iv2) = 0.d0
       do ir=2,ubound(gvv,1) !nr
          r = (ir-1)*dr
          exchemv = 0.d0
          do iv1=1,ubound(mtv,1) !nv
             ivv = jvv(iv1,iv2)
             hvv = gvv(ir,ivv) - 1.d0
             h2c = 0.5d0*hvv*hvv - cvv(ir,ivv)

             exchemv = exchemv + densityv(iv1)*mtv(iv2) &
                  *(h2c)
          enddo
          exchem(iv2) = exchem(iv2) + exchemv*r**2
       enddo
    enddo
    exchem = exchem * 4.d0*pi*dr
  end function rism1d_hnc_exchem_r


  !! Calculate the ionic(?) excess chemical potential in kT.  This seems to take
  !! the contribution only from the last solvent type.  I'm not sure what the point is
  !! here.
  !! IN:
  !!    this : the closure object
  !!    gvv  : site-site pair correlation function
  !!    cvv  : site-site direct correlation function
  !!    mtv  : multiplicity of each site
  !!    densityv : number density of each site
  !!    dr   : R-space grid spacing
  function rism1d_hnc_exchemIon(this, gvv, cvv, mtv, jvv, densityv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_hnc), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:), densityv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))
    _REAL_ :: r, h2c,hvv, exchemv
    integer :: ivv, iv1, iv2, ir

    do iv2=1,ubound(mtv,1) !nv
       exchem(iv2) = 0.d0
       do ir=2,ubound(gvv,1) !nr
          r = (ir-1)*dr
          exchemv = 0.d0
          !          do iv1=1,ubound(mtv,1) !nv
          iv1=ubound(mtv,1)
          ivv = jvv(iv1,iv2)
          hvv = gvv(ir,ivv) - 1.d0
          h2c = 0.5d0*hvv*hvv - cvv(ir,ivv)
          !                     exchemv = exchemv + densityv(iv1)*mtv(iv2) &
          exchemv = exchemv + densityv(iv2)*mtv(iv1) &
               *(h2c - 0.5d0*hvv*cvv(ir,ivv))
          !         enddo
          exchem(iv2) = exchem(iv2) + exchemv*r**2
       enddo
    enddo
    exchem = exchem * 4.d0*pi*dr
  end function rism1d_hnc_exchemIon


  !! Calculate the solvation energy in kT
  !! IN:
  !!    this   : the closure object
  !!    gvv    : site-site pair correlation function
  !!    cvv    : site-site direct correlation function
  !!    gvv_dT  : temperature derivative of site-site pair correlation function
  !!    cvv_dT  : temperature derivative of site-site direct correlation function
  !!    mtv    : multiplicity of each site
  !!    jvv    : triangle packing pointer
  !!    densityv   : number density of each site
  !!    dr     : R-space grid spacing
  function rism1d_hnc_solvene(this, gvv, cvv, gvv_dT, cvv_dT, mtv, jvv, densityv, dr) result(solvene)
    use constants, only : pi
    implicit none
    type(rism1d_hnc), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:), gvv_dT(:,:),cvv_dT(:,:), densityv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: solvene(ubound(mtv,1))
    _REAL_ :: r, h2c,hvv, solvenev, hvv_dT
    integer :: ivv, iv1, iv2, ir

    do iv2=1,ubound(mtv,1) !nv
       solvene(iv2) = 0.d0
       do ir=2,ubound(gvv,1) !nr
          r = (ir-1)*dr
          solvenev = 0.d0
          do iv1=1,ubound(mtv,1) !nv
             ivv = jvv(iv1,iv2)
             hvv = gvv(ir,ivv) - 1.d0
             hvv_dT = gvv_dT(ir,ivv)
             h2c = hvv*hvv_dT - cvv_dT(ir,ivv)
             solvenev = solvenev + densityv(iv1)*mtv(iv2) &
                  *(h2c - 0.5d0*hvv_dT*cvv(ir,ivv) - 0.5d0*hvv*cvv_dT(ir,ivv))
          enddo
          solvene(iv2) = solvene(iv2) + solvenev*r**2
       enddo
    enddo
    solvene = -1.d0 * solvene * 4.d0*pi*dr
  end function rism1d_hnc_solvene


  !> Calculate the R-space contribution to the free energy and pressure
  !! (free energy route) in internal units.  This is the same quanitity
  !! for both values.
  !! @param[in] this The closure object.
  !! @param[in] gvv  Site-site pair correlation function.
  !! @param[in] cvv  Site-site direct correlation function.
  !! @param[in] mtv  Multiplicity of each site.
  !! @param[in] densityv Number density of each site.
  !! @param[in] dr   R-space grid spacing.
  !! @return R-space contribution to the free energy and pressure
  !!  (free energy route) in internal units.
  function rism1d_hnc_freeEnergyPressure_r(this, gvv, cvv, mtv, densityv, dr) result(fe)
    use constants, only : pi
    implicit none
    type(rism1d_hnc), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:),densityv(:),dr
    integer, intent(in) :: mtv(:)
    _REAL_ :: fe
    _REAL_ :: fev
    _REAL_ :: r, hvv, h2c
    integer :: ir,ivv,iv1,iv2, cnt
    !R-space
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
             h2c = 0.5d0*hvv*hvv - cvv(ir,ivv)
             fev = fev + cnt*densityv(iv1)*densityv(iv2)*h2c
          enddo
       enddo
       fe = fe + fev*r**2
    enddo
    fe = fe * 2.d0*pi*dr
  end function rism1d_hnc_freeEnergyPressure_r


  !> Frees memory and resets the HNC closure
  !! @param[in,out] this HNC object.
  subroutine rism1d_hnc_destroy(this)
    implicit none
    type(rism1d_hnc), intent(inout) :: this

    !TODO: Free memory?
  end subroutine rism1d_hnc_destroy
end module rism1d_hnc_c
