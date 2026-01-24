! <compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2011-2012 by 
!Andriy Kovalenko, Tyler Luchko and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
!ibid. 112:10391-10417 (2000).   
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 

#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Modified Verlet closure.
!!!M. Kinoshita, T. Imai, A. Kovalenko and F. Hirata. Chem. Phys. Lett. 2001, 348, 337-342.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module rism1d_mv0_c
  use safemem
  !the MV0 type
  type rism1d_mv0
     !uvvlr :: long range LJ potential using the WCA partition
     !bvv   :: bridge function
     _REAL_, pointer :: uvvlr(:,:)=>NULL(), bvv(:,:)=>NULL()
     !charged :: flag to indicate if we are using the full potential or just the LJ potential
     logical :: charged
  end type rism1d_mv0

  public rism1d_mv0_new, rism1d_mv0_destroy, rism1d_mv0_gvv, rism1d_mv0_exchem,&
       rism1d_mv0_exchemIon, rism1d_mv0_pressureFE, rism1d_mv0_freeEnergy, &
       rism1d_calcUvv, rism1d_mv0_isCharged, rism1d_mv0_useCharged

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initializes the MV0 closure
!!!IN:
!!!   this : MV0 object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_mv0_new(this)
    implicit none
    type(rism1d_mv0), intent(inout) :: this
    this%charged = .false.
  end subroutine rism1d_mv0_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns true if the potential with charges is being used.
!!!IN:
!!!   this : the closure object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_mv0_charged(this) result(charged)
    implicit none
    type(rism1d_mv0), intent(in) :: this
    logical :: charged
    charged = this%charged
  end function rism1d_mv0_charged

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Set to true to use the full potential, false otherwise.
!!!IN:
!!!   this : the closure object
!!!   charged : true to use full potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_mv0_useCharged(this,charged)
    implicit none
    type(rism1d_mv0), intent(inout) :: this
    logical,intent(in) :: charged
    this%charged =  charged
  end subroutine rism1d_mv0_useCharged

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Gvv from Uvv, Hvv, and Cvv using the MV0 closure
!!!IN:
!!!   this : the MV0 closure object
!!!   gvv  : site-site pair correlation function
!!!   uvv  : site-site potential
!!!   hvv  : site-site total correlation function
!!!   cvv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_mv0_gvv(this,gvv, uvv, hvv, cvv)
    implicit none
    type(rism1d_mv0), intent(inout) :: this
    _REAL_, intent(out) :: gvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),hvv(:,:),cvv(:,:)
    integer :: ivv, ir
    _REAL_ :: tvv, tvvr
    _REAL_, parameter :: alpha=0.8d0

    if(.not. this%charged) then
       do ivv = 1, ubound(gvv,2)
          do ir = 2, ubound(gvv,1)
             tvv = -uvv(ir,ivv) +hvv(ir,ivv)-cvv(ir,ivv)
             tvvr = -this%uvvlr(ir,ivv) +hvv(ir,ivv)-cvv(ir,ivv)
             this%bvv(ir,ivv) = -0.5d0*tvvr**2 / (1.d0 + alpha*max(tvvr,0.d0))
             gvv(ir,ivv) = exp(min(tvv+this%bvv(ir,ivv),10d0))
          end do
       end do
    else
       do ivv = 1, ubound(gvv,2)
          do ir = 2, ubound(gvv,1)
             tvv = -uvv(ir,ivv) +hvv(ir,ivv)-cvv(ir,ivv)
             gvv(ir,ivv) = exp(min(tvv+this%bvv(ir,ivv),10d0))      
          end do
       end do
    endif
  end subroutine rism1d_mv0_gvv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Bvv (bridge function) from Uvv, Gvv, and Cvv using the MV0 closure
!!!IN:
!!!   this : the MV0 closure object
!!!   uvv  : site-site potential
!!!   gvv  : site-site total correlation function
!!!   cvv  : site-site direct correlation function
!!!OUT:
!!!   bvv  : (pointer) site-site bridge function.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_mv0_bvv(this, uvv, gvv, cvv) result(bvv)
    implicit none
    type(rism1d_mv0), intent(in) :: this
    _REAL_,pointer :: bvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),gvv(:,:),cvv(:,:)
    nullify(bvv)
    bvv =>safemem_realloc(bvv,ubound(this%bvv,1),ubound(this%bvv,2))
    bvv = this%bvv
  end function rism1d_mv0_bvv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the excess chemical potential in kT
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site pair correlation function
!!!   cvv  : site-site direct correlation function
!!!   mtv  : multiplicity of each site
!!!   densityv : number density of each site
!!!   densitytrgt : (optional) The final, physical density for thermodynamics.  This
!!!             can be used as an effective correction for some closures that
!!!             over estimate the pressure (e.g. MV0 and MV0)
!!!   dr   : r-space grid spacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_mv0_exchem(this, gvv, cvv, mtv, jvv, densityv, densitytrgtv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_mv0), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:), densityv(:), densitytrgtv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))
    
    !No analytic expression for chemical potential so return NaN
    exchem = 0d0
    exchem = exchem/exchem

end function rism1d_mv0_exchem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the ionic(?) excess chemical potential in kT.  This seems to take
!!!the contribution only from the last solvent type.  I'm not sure what the point is
!!!here.
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site pair correlation function
!!!   cvv  : site-site direct correlation function
!!!   mtv  : multiplicity of each site
!!!   densityv : number density of each site
!!!   densitytrgt : (optional) The final, physical density for thermodynamics.  This
!!!             can be used as an effective correction for some closures that
!!!             over estimate the pressure (e.g. MV0 and MV0)
!!!   dr   : r-space grid spacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_mv0_exchemIon(this, gvv, cvv, mtv, jvv, densityv, densitytrgtv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_mv0), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:), densityv(:), densitytrgtv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))

    !No analytic expression for chemical potential so return NaN
    exchem = 0d0
    exchem = exchem/exchem
  end function rism1d_mv0_exchemIon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the r-space contribution to the pressure in internal units
!!!using the free energy route
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site pair correlation function
!!!   cvv  : site-site direct correlation function
!!!   mtv  : multiplicity of each site
!!!   densityv : number density of each site
!!!   dr   : r-space grid spacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_mv0_pressureFE(this, gvv, cvv, mtv, densityv, dr) result(pr)
    use constants, only : pi
    implicit none
    type(rism1d_mv0), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:),densityv(:),dr
    integer, intent(in) :: mtv(:)
    _REAL_ :: pr

    !No analytic expression for chemical potential so return NaN
    pr = 0d0
    pr = pr/pr
  end function rism1d_mv0_pressureFE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the r-space contribution to the free energy in internal units
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site pair correlation function
!!!   cvv  : site-site direct correlation function
!!!   mtv  : multiplicity of each site
!!!   densityv : number density of each site
!!!   dr   : r-space grid spacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_mv0_freeEnergy(this, gvv, cvv, mtv, densityv, dr) result(fe)
    use constants, only : pi
    implicit none
    type(rism1d_mv0), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:),densityv(:),dr
    integer, intent(in) :: mtv(:)
    _REAL_ :: fe

    !No analytic expression for chemical potential so return NaN
    fe = 0d0
    fe = fe/fe
  end function rism1d_mv0_freeEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets the MV0 closure
!!!IN:
!!!   this : MV0 object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_mv0_destroy(this)
    use safemem
    implicit none
    type(rism1d_mv0), intent(inout) :: this
    if(safemem_dealloc(this%uvvlr) /=0) then
       call rism_report_error("deallocating UvvLR in MV0 closure failed")
    end if
    if(safemem_dealloc(this%bvv) /=0) then
       call rism_report_error("deallocating Bvv in MV0 closure failed")
    end if
  end subroutine rism1d_mv0_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvent site-site potential
!!!IN:
!!!   this : the closure object
!!!   nr   : number of grid points for potential/correlation functions
!!!   qv   : charge of each site in e
!!!   epsv : LJ epsilon value for each site in kT
!!!   rminv: LJ r_min value for each site in A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_mv0_calcUvv(this, nr, dr, qv, epsv, rminv, smear)
    use safemem
    implicit none
    type(rism1d_mv0), intent(inout) :: this
    _REAL_, intent(in) :: dr, smear
    integer, intent(in) :: nr
    _REAL_, intent(in) :: qv(:),epsv(:),rminv(:)

    _REAL_:: epsvv, rminvv, qvv, r, ri6, rs, rs6, usr
    integer :: nv, nvv
    integer :: ivv, iv1, iv2, ir

    nv = ubound(qv,1)
    nvv = nv*(nv+1)/2
    this%uvvLR => safemem_realloc(this%uvvlr,nr,nvv,.false.)
    this%bvv => safemem_realloc(this%bvv,nr,nvv,.false.)
    ivv = 0
    do iv2=1,nv
       do iv1=1,iv2
          ivv = ivv + 1

          !.......................... 12-6 LJ potential ..........................
          rminvv = (rminv(iv1)+rminv(iv2))
          epsvv = sqrt(epsv(iv1)*epsv(iv2))
          do ir=1,nr
             r = (ir-1)*dr
             rs = r/rminvv
             rs = max( rs, 0.002d0*2**(1d0/6d0))
             rs6 = rs**6
             ri6 = 1.d0/rs6
             usr = (ri6-2.d0)*ri6
             usr = usr * epsvv
             if(rs6 < 1d0)then
                this%uvvlr(ir,ivv) = -epsvv
             else 
                this%uvvlr(ir,ivv) = usr
             end if
          enddo
       enddo
    enddo
  end subroutine rism1d_mv0_calcUvv
end module rism1d_mv0_c
