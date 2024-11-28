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

!> Kovalenko-Hirata (KH) closure class for 1D-RISM.
!!
!! Kovalenko, A.; Hirata, F. J. Chem. Phys. 1999, 110, 10095–10112.
!!
!! This class provides KH specific contributions for thermodynamic
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
!! KH has a similar form where the R-space contribution is closure
!! dependent and the k-space is not.  The Pettitt-Rossky equivalent
!! can be found in:
!!
!! Kovalenko, A.; Hirata, F. J. Chem. Phys. 1999, 110, 10095–10112.
!!
!! For the total excess chemical potential, we use the simplest
!! formulation, Pettitt-Rossky/Kovalenko-Hirata.  All three
!! forumlations have the same R-space form and this is calculated
!! separately in rism1d_kh_exchem_r().
!!
!! Free energy and pressure
!! ------------------------
!! Only the closure dependent, R-space component of the free energy
!! and pressure are provided here.  This contribution is the same for
!! both quantities and is calculated with a single function,
!! rism1d_kh_freeEnergyPressure_r().
module rism1d_kh_c
  use safemem
  ! The KH type is actually empty.
  type rism1d_kh
     logical kh
  end type rism1d_kh

  public rism1d_kh_new, rism1d_kh_destroy, rism1d_kh_gvv, rism1d_kh_gvv_dT, rism1d_kh_solvene

contains

  !> Initializes the KH closure.
  !! @param[in] this KH object.
  subroutine rism1d_kh_new(this)
    implicit none
    type(rism1d_kh), intent(inout) :: this
  end subroutine rism1d_kh_new


  !> Calculates Gvv from Uvv, Hvv, and Cvv using the KH closure.
  !! @param[in] this The KH closure object.
  !! @param[out] gvv Site-site pair correlation function.
  !! @param[in] uvv Site-site potential.
  !! @param[in] hvv Site-site total correlation function.
  !! @param[in] cvv Site-site direct correlation function.
  subroutine rism1d_kh_gvv(this,gvv, uvv, hvv, cvv)
    implicit none
    type(rism1d_kh), intent(in) :: this
    _REAL_, intent(out) :: gvv(:,:)
    _REAL_, intent(in) :: uvv(:,:), hvv(:,:), cvv(:,:)
    integer :: ivv, ir
    _REAL_ :: tvv
    do ivv = 1, ubound(gvv, 2)
       do ir = 1, ubound(gvv, 1)
          tvv = -uvv(ir, ivv) + hvv(ir, ivv) - cvv(ir, ivv)
          if(tvv>=0d0)then
             gvv(ir, ivv) = 1d0 + tvv
          else
             gvv(ir, ivv) = exp(tvv)
          end if
       end do
    end do
  end subroutine rism1d_kh_gvv


  !> Calculates temperature derivative Gvv from Gvv, Uvv, Hvv_dT, and
  !! Cvv_dT using the associated closure.
  !! param[in] this The closure object.
  !! param[out] gvv_dT Site-site temperature derivative pair correlation function.
  !! param[in] gvv Site-site pair correlation function.
  !! param[in] hvv_dT Site-site temperature derivative total correlation function.
  !! param[in] cvv_dT Site-site temperature derivative direct correlation function.
  subroutine rism1d_kh_gvv_dT(this,gvv_dT, uvv, gvv, hvv_dT, cvv_dT)
    implicit none
    type(rism1d_kh), intent(in) :: this
    _REAL_, intent(out) :: gvv_dT(:,:)
    _REAL_, intent(in) :: uvv(:,:),gvv(:,:),hvv_dT(:,:),cvv_dT(:,:)
    integer :: ivv, ir
    _REAL_ :: tvv_dT
    do ivv = 1, ubound(gvv_dT,2)
       do ir = 1, ubound(gvv_dT,1)
          tvv_dT = uvv(ir,ivv) +hvv_dT(ir,ivv)-cvv_dT(ir,ivv)
          if(gvv(ir,ivv)>=1d0)then
             gvv_dT(ir,ivv) = tvv_dT
          else
             gvv_dT(ir,ivv) = tvv_dT*gvv(ir,ivv)
          end if
       end do
    end do
  end subroutine rism1d_kh_gvv_dT


  !! Calculates Bvv (bridge function) from Uvv, Hvv, and Cvv using the KH closure
  !! IN:
  !!    thisthe KH closure object
  !!    uvv site-site potential
  !!    gvv site-site pair correlation function
  !!    cvv site-site direct correlation function
  !! OUT:
  !!    bvv site-site bridge function

  function rism1d_kh_bvv(this, uvv, gvv, cvv) result(bvv)
    implicit none
    type(rism1d_kh), intent(in) :: this
    _REAL_, pointer :: bvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),gvv(:,:),cvv(:,:)
    integer :: ivv, ir
    _REAL_ :: tvv
    nullify(bvv)
    bvv =>safemem_realloc(bvv,ubound(gvv,1),ubound(gvv,2))
    do ivv = 1, ubound(bvv,2)
       do ir = 1, ubound(bvv,1)
          tvv = -uvv(ir,ivv) +gvv(ir,ivv)-1d0-cvv(ir,ivv)
          if(tvv>=0d0)then
             bvv(ir,ivv) = -tvv + log(1 + tvv)
          else
             bvv(ir,ivv) = 0
          end if
       end do
    end do
  end function rism1d_kh_bvv


  !! Calculate the excess chemical potential in kT using the
  !! Pettitt-Rossky formula
  !! IN:
  !!    thisthe closure object
  !!    gvv site-site pair correlation function
  !!    cvv site-site direct correlation function
  !!    mtv multiplicity of each site
  !!    densityvnumber density of each site
  !!    dr  R-space grid spacing
  function rism1d_kh_exchem_PR(this, gvv, cvv, mtv, jvv, densityv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_kh), intent(in) :: this
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
          do iv1=1,ubound(mtv,1) !nv
             ivv = jvv(iv1,iv2)
             hvv = gvv(ir,ivv) - 1.d0
             if (hvv >= 0.d0)  then
                h2c = - cvv(ir,ivv)
             else
                h2c = 0.5d0*hvv*hvv - cvv(ir,ivv)
             endif
             exchemv = exchemv + densityv(iv1)*mtv(iv2) &
                  *(h2c - 0.5d0*hvv*cvv(ir,ivv))
          enddo
          exchem(iv2) = exchem(iv2) + exchemv*r**2
       enddo
    enddo
    exchem = exchem * 4.d0*pi*dr
  end function rism1d_kh_exchem_PR


  !! Calculates the closure dependent, R-space component of excess chemical
  !! potential in kT.
  !!
  !!   \int 1/2*h^2*\theta(-h) - c dr
  !!
  !! This is the R-space part of the KH version of the Singer-Chandler
  !! and Schmeer-Maurer forms and is the Pettitt-Rossky form without the
  !! 'hc/2' term.
  !! IN:
  !!    this the closure object
  !!    gvv site-site pair correlation function
  !!    cvv site-site direct correlation function
  !!    mtv multiplicity of each site
  !!    densityv number density of each site
  !!    dr R-space grid spacing
  function rism1d_kh_exchem_r(this, gvv, cvv, mtv, jvv, densityv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_kh), intent(in) :: this
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
          do iv1=1,ubound(mtv,1) !nv
             ivv = jvv(iv1,iv2)
             hvv = gvv(ir,ivv) - 1.d0
             if (hvv >= 0.d0)  then
                h2c = - cvv(ir,ivv)
             else
                h2c = 0.5d0*hvv*hvv - cvv(ir,ivv)
             endif
             exchemv = exchemv + densityv(iv1)*mtv(iv2) &
                  *(h2c)
          enddo
          exchem(iv2) = exchem(iv2) + exchemv*r**2
       enddo
    enddo
    exchem = exchem * 4.d0*pi*dr
  end function rism1d_kh_exchem_r


  !! Calculate the ionic(?) excess chemical potential in kT.  This seems to take
  !! the contribution only from the last solvent type.  I'm not sure what the point is
  !! here.
  !! IN:
  !!    thisthe closure object
  !!    gvv site-site pair correlation function
  !!    cvv site-site direct correlation function
  !!    mtv multiplicity of each site
  !!    densityvnumber density of each site
  !!    dr  R-space grid spacing
  function rism1d_kh_exchemIon(this, gvv, cvv, mtv, jvv, densityv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_kh), intent(in) :: this
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
          if (hvv >= 0.d0)  then
             h2c = - cvv(ir,ivv)
          else
             h2c = 0.5d0*hvv*hvv - cvv(ir,ivv)
          endif
          exchemv = exchemv + densityv(iv2)*mtv(iv1) &
               *(h2c - 0.5d0*hvv*cvv(ir,ivv))
          !         enddo
          exchem(iv2) = exchem(iv2) + exchemv*r**2
       enddo
    enddo
    exchem = exchem * 4.d0*pi*dr
  end function rism1d_kh_exchemIon


  !! Calculate the solvation energy in kT
  !! IN:
  !!    thisthe closure object
  !!    gvv site-site pair correlation function
  !!    cvv site-site direct correlation function
  !!    gvv_dT temperature derivative of site-site pair correlation function
  !!    cvv_dT temperature derivative of site-site direct correlation function
  !!    mtv multiplicity of each site
  !!    densityvnumber density of each site
  !!    dr  R-space grid spacing
  function rism1d_kh_solvene(this, gvv, cvv, gvv_dT, cvv_dT, mtv, jvv, densityv, &
       dr) result(solvene)
    use constants, only : pi
    implicit none
    type(rism1d_kh), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:), gvv_dT(:,:), cvv_dT(:,:), densityv(:), dr
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
             if (hvv >= 0.d0)  then
                h2c = - cvv_dT(ir,ivv)
             else
                h2c = hvv*hvv_dT - cvv_dT(ir,ivv)
             endif
             solvenev = solvenev + densityv(iv1)*mtv(iv2) &
                  *(h2c - 0.5d0*hvv_dT*cvv(ir,ivv) - 0.5d0*hvv*cvv_dT(ir,ivv))
          enddo
          solvene(iv2) = solvene(iv2) + solvenev*r**2
       enddo
    enddo
    solvene = -1.d0 * solvene * 4.d0*pi*dr
  end function rism1d_kh_solvene


  !> Calculate the R-space contribution to the free energy and
  !! pressure (free energy route) in internal units.  This is the same
  !! quanitity for both values.
  !! @param[in] this The closure object.
  !! @param[in] gvv Site-site pair correlation function.
  !! @param[in] cvv Site-site direct correlation function.
  !! @param[in] mtv Multiplicity of each site.
  !! @param[in] densityv Number density of each site.
  !! @param[in] dr R-space grid spacing.
  !! @return R-space contribution to the free energy and
  !! pressure (free energy route) in internal units
  function rism1d_kh_freeEnergyPressure_r(this, gvv, cvv, mtv, densityv, dr) result(fe)
    use constants, only : pi
    implicit none
    type(rism1d_kh), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:), densityv(:), dr
    integer, intent(in) :: mtv(:)
    _REAL_ :: fe
    _REAL_ :: fev
    _REAL_ :: r, hvv, h2c
    integer :: ir, ivv, iv1, iv2, cnt
    ! R-space.
    fe = 0.d0
    do ir = 2, ubound(gvv, 1) !nr
       r = (ir - 1) * dr
       fev = 0.d0
       ivv = 0
       do iv2 = 1, ubound(mtv, 1) !nv
          do iv1 = 1, iv2
             ivv = ivv + 1
             ! If the sites are different, we need to double count the
             ! contribution.
             if (iv1 == iv2)  then
                cnt = 1
             else
                cnt = 2
             endif
             hvv = gvv(ir, ivv) - 1.d0
             if (hvv >= 0.d0)  then
                h2c = - cvv(ir, ivv)
             else
                h2c = 0.5d0 * hvv * hvv - cvv(ir,ivv)
             endif
             fev = fev + cnt * densityv(iv1) * densityv(iv2) * h2c
          enddo
       enddo
       fe = fe + fev * r**2
    enddo
    fe = fe * 2.d0 * pi * dr
  end function rism1d_kh_freeEnergyPressure_r


  !> Frees memory and resets the KH closure.
  !! @param[in,out] this KH object.
  subroutine rism1d_kh_destroy(this)
    implicit none
    type(rism1d_kh), intent(inout) :: this

    !TODO: Free memory?
  end subroutine rism1d_kh_destroy

end module rism1d_kh_c
