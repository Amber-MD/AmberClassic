! <compile=optimized>
#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! bridge for development
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module rism1d_devb_c
  use rism_report_c
  use safemem
  !the PSE-n type
  type rism1d_devb
     !coeff :: coefficients
     _REAL_, pointer :: coeff(:)=>NULL()
     integer :: order=0
  end type rism1d_devb

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initializes the DEVB closure
!!!IN:
!!!   this : DEVB object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_devb_new(this,coeff)
    implicit none
    type(rism1d_devb), intent(inout) :: this
    _REAL_, intent(in) :: coeff(:)

    this%coeff => safemem_realloc(this%coeff,ubound(coeff,1))
    this%coeff = coeff
  end subroutine rism1d_devb_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Gvv from Uvv, Hvv, and Cvv using the DEVB closure
!!!IN:
!!!   this : the DEVB closure object
!!!   gvv  : site-site pair correlation function
!!!   uvv  : site-site potential
!!!   hvv  : site-site total correlation function
!!!   cvv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_devb_gvv(this,gvv, uvv, hvv, cvv)
    implicit none
    type(rism1d_devb), intent(in) :: this
    _REAL_, intent(out) :: gvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),hvv(:,:),cvv(:,:)
    integer :: ivv, ir, i, order
    _REAL_ :: tvv, cuvv, orderfac

    do ivv = 1, ubound(gvv,2)
       do ir = 1, ubound(gvv,1)
          cuvv = uvv(ir,ivv) + cvv(ir,ivv) 
          tvv = -uvv(ir,ivv) +hvv(ir,ivv)-cvv(ir,ivv)
!          gvv(ir,ivv) = exp( this%coeff(ivv) * tvv )
          gvv(ir,ivv) =exp( tvv + cuvv*this%coeff(ivv) )
          !if ( tvv >= 0.d0 ) then 
          !   gvv(ir,ivv) = exp( tvv )
          !else
          !   gvv(ir,ivv) = exp( tvv + this%coeff(ivv) * tvv**2 / ( tvv - 1.0d0 ) )
          !end if
       end do
    end do
  end subroutine rism1d_devb_gvv

!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!Calculates dTGvv from Gvv, Uvv, dTHvv, and dTCvv using the associated 
!!$!!!closure
!!$!!!IN:
!!$!!!   this   : the closure object
!!$!!!   gvvdt  : site-site temperature derivative pair correlation function
!!$!!!   gvv    : site-site pair correlation function
!!$!!!   cvv  : site-site direct correlation function
!!$!!!   hvvdt  : site-site temperature derivative total correlation function
!!$!!!   cvvdt  : site-site temperature derivative direct correlation function
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  subroutine rism1d_devb_gvvdt(this, gvvdt, uvv, gvv, cvv, hvvdt, cvvdt)
!!$    implicit none
!!$    type(rism1d_devb), intent(in) :: thish:jjkk!!$    _REAL_, intent(out) :: gvvdt(:,:)
!!$    _REAL_, intent(in) :: uvv(:,:),gvv(:,:),cvv(:,:),hvvdt(:,:),cvvdt(:,:)
!!$    integer :: ivv, ir, i
!!$    _REAL_ :: tvvdt, tvv, orderfac
!!$    do ivv = 1, ubound(gvvdt,2)
!!$       do ir = 1, ubound(gvvdt,1)
!!$          tvvdt = uvv(ir,ivv) +hvvdt(ir,ivv)-cvvdt(ir,ivv)
!!$          tvv = -uvv(ir,ivv) +gvv(ir,ivv)-1d0-cvv(ir,ivv)
!!$          gvvdt(ir,ivv) = 1d0
!!$          orderfac = 1d0
!!$          do i=1,this%order-1
!!$             orderfac = orderfac*i
!!$             gvvdt(ir,ivv) = gvvdt(ir,ivv) + (tvv**i)/orderfac
!!$          end do
!!$          gvvdt(ir,ivv) = gvvdt(ir,ivv)*tvvdt
!!$       end do
!!$    end do
!!$  end subroutine rism1d_devb_gvvdt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Bvv (bridge function) from Uvv, Gvv, and Cvv using the DEVB closure
!!!IN:
!!!   this : the DEVB closure object
!!!   uvv  : site-site potential
!!!   gvv  : site-site total correlation function
!!!   cvv  : site-site direct correlation function
!!!OUT:
!!!   bvv  : site-site bridge function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_devb_bvv(this, uvv, gvv, cvv) result(bvv)
    implicit none
    type(rism1d_devb), intent(in) :: this
    _REAL_, pointer :: bvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),gvv(:,:),cvv(:,:)
    integer :: ivv, ir, i
    _REAL_ :: tvv
    nullify(bvv)
    bvv => safemem_realloc(bvv,ubound(gvv,1),ubound(gvv,2))
    do ivv = 1, ubound(bvv,2)
       do ir = 1, ubound(bvv,1)
          tvv = -uvv(ir,ivv) +gvv(ir,ivv)-1d0-cvv(ir,ivv)
          bvv(ir,ivv) = ( this%coeff(ivv) - 1.d0 ) * tvv
!          if ( tvv >= 0.d0 ) then
!             bvv(ir,ivv) = this%coeff(ivv) * tvv**2 / ( 1.d0 + tvv )
!          else
!             bvv(ir,ivv) = this%coeff(ivv) * tvv**2 / ( this%coeff(ivv*3) - tvv )**2
!          end if
       end do
    end do
  end function rism1d_devb_bvv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the excess chemical potential in kT using the
!!!Pettitt-Rossky formula
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site pair correlation function
!!!   uvv  : site-site potential
!!!   cvv  : site-site direct correlation function
!!!   mtv  : multiplicity of each site
!!!   densityv : number density of each site
!!!   dr   : r-space grid spacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_devb_exchem_PR(this, gvv, uvv, cvv, mtv, jvv, densityv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_devb), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),uvv(:,:),cvv(:,:), densityv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))
    !hnc    : HNC contribution
    !tsvv   : t^* = -u +h -c
    !bridge : contribution from the bridge expansion
    _REAL_ :: r, hnc,hvv,tsvv, exchemv,bridge
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
             bridge = hvv - 1.0d0 / this%coeff(ivv) * ( exp( this%coeff(ivv) * tsvv ) - 1.0d0 )
!             if ( tsvv >= 0.0d0 ) then
!                bridge = 0.0d0
!             else
!                bridge = hvv - 1.0d0 / ( this%coeff(ivv) + 1.0d0 ) * ( exp( ( this%coeff(ivv) + 1.0d0 ) * tsvv ) - 1.0d0 )
!             end if
             hnc = 0.5d0*hvv*hvv - cvv(ir,ivv)  - 0.5d0*hvv*cvv(ir,ivv)           
             exchemv = exchemv + densityv(iv1)*mtv(iv2) &
                  *(hnc+bridge)
          enddo
          exchem(iv2) = exchem(iv2) + exchemv*r**2
       enddo
!write(*,*) iv2, exchem(iv2)
    enddo
    exchem = exchem * 4.d0*pi*dr
  end function rism1d_devb_exchem_PR

!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!Calculate the solvation energy in kT
!!$!!!IN:
!!$!!!   this : the closure object
!!$!!!   gvv  : site-site pair correlation function
!!$!!!   uvv  : site-iste potential
!!$!!!   cvv  : site-site direct correlation function
!!$!!!   gvvdt  : temperature derivative of site-site pair correlation function
!!$!!!   cvvdt  : temperature derivative of site-site direct correlation function
!!$!!!   mtv  : multiplicity of each site
!!$!!!   densityv : number density of each site
!!$!!!   densitytrgt : (optional) The final, physical density for thermodynamics.  This
!!$!!!             can be used as an effective correction for some closures that
!!$!!!             over estimate the pressure (e.g. HNC and KH)
!!$!!!   dr   : r-space grid spacing
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  function rism1d_devb_solene(this, gvv, uvv, cvv, gvvdt, cvvdt, mtv, jvv, densityv, densitytrgtv, dr) result(solene)
!!$    use constants, only : pi
!!$    implicit none
!!$    type(rism1d_devb), intent(in) :: this
!!$    _REAL_, intent(in) :: gvv(:,:),uvv(:,:), cvv(:,:), gvvdt(:,:), cvvdt(:,:), &
!!$         densityv(:), densitytrgtv(:), dr
!!$    integer, intent(in) :: mtv(:),jvv(:,:)
!!$    _REAL_ :: solene(ubound(mtv,1))
!!$    _REAL_ :: r, h2c,hvv, solenev, hvvdt
!!$    integer :: ivv, iv1, iv2, ir
!!$    _REAL_ :: orderfac, tsvv, tvvdt, series
!!$    orderfac = this%order1fac/this%order1
!!$    
!!$    do iv2=1,ubound(mtv,1) !nv
!!$       solene(iv2) = 0.d0
!!$       do ir=2,ubound(gvv,1) !nr
!!$          r = (ir-1)*dr
!!$          solenev = 0.d0
!!$          do iv1=1,ubound(mtv,1) !nv
!!$             ivv = jvv(iv1,iv2)
!!$             hvv = gvv(ir,ivv) - 1.d0
!!$             hvvdt = gvvdt(ir,ivv)
!!$             if (hvv >= 0.d0)  then
!!$                tsvv = -uvv(ir,ivv) + hvv -cvv(ir,ivv)
!!$                tvvdt = uvv(ir,ivv) + hvvdt - cvvdt(ir,ivv)
!!$                series = tsvv**(this%order)/orderfac*tvvdt
!!$             else
!!$                series = 0
!!$             end if
!!$             h2c = hvv*hvvdt - cvvdt(ir,ivv)
!!$             solenev = solenev + densitytrgtv(iv1)*mtv(iv2) &
!!$                  *(h2c - 0.5d0*hvvdt*cvv(ir,ivv) - 0.5d0*hvv*cvvdt(ir,ivv) -series)
!!$          enddo
!!$          solene(iv2) = solene(iv2) + solenev*r**2
!!$       enddo
!!$    enddo
!!$    solene = -1.d0 * solene * 4.d0*pi*dr
!!$  end function rism1d_devb_solene
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!Calculate the r-space contribution to the pressure in internal units
!!$!!!using the free energy route
!!$!!!IN:
!!$!!!   this : the closure object
!!$!!!   gvv  : site-site pair correlation function
!!$!!!   cvv  : site-site direct correlation function
!!$!!!   mtv  : multiplicity of each site
!!$!!!   densityv : number density of each site
!!$!!!   dr   : r-space grid spacing
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  function rism1d_devb_pressureFE(this, gvv, uvv, cvv, mtv, densityv, dr) result(pr)
!!$    use constants, only : pi
!!$    implicit none
!!$    type(rism1d_devb), intent(in) :: this
!!$    _REAL_, intent(in) :: gvv(:,:),uvv(:,:),cvv(:,:),densityv(:),dr
!!$    integer, intent(in) :: mtv(:)
!!$    _REAL_ :: pr
!!$    _REAL_ :: prv
!!$    !hnc    : h^2*c/2 contribution
!!$    !tsvv   : t^* = u +h -c
!!$    !series : contribution from the series expansion
!!$    _REAL_ :: r, hvv, h2c, tsvv, series
!!$    integer :: ir,ivv,iv1,iv2, cnt, i
!!$    !r-space
!!$    pr = 0.d0
!!$    do ir=2,ubound(gvv,1) !nr
!!$       r = (ir-1)*dr
!!$       prv = 0.d0
!!$       ivv = 0
!!$       do iv2=1,ubound(mtv,1) !nv
!!$          do iv1=1,iv2
!!$             ivv = ivv + 1
!!$             !if the sites are different, we need to double count the contribution
!!$             if (iv1 == iv2)  then
!!$                cnt = 1
!!$             else
!!$                cnt = 2
!!$             endif
!!$             hvv = gvv(ir,ivv) - 1.d0
!!$             tsvv = -uvv(ir,ivv) + hvv -cvv(ir,ivv)
!!$             if (tsvv >= 0.d0)  then
!!$                series = tsvv**(this%order1)/this%order1fac
!!$             else
!!$                series = 0
!!$             endif
!!$             h2c = 0.5d0*hvv*hvv - cvv(ir,ivv)
!!$             prv = prv + cnt*mtv(iv1)*densityv(iv2)*(h2c-series)
!!$          enddo
!!$       enddo
!!$       pr = pr + prv*r**2
!!$    enddo
!!$    pr = pr * 2.d0*pi*dr
!!$  end function rism1d_devb_pressureFE
!!$
  function rism1d_devb_freeEnergyPressure_r(this, gvv, uvv, cvv, mtv, densityv, dr) result(fe)
    use constants, only : pi
    implicit none
    type(rism1d_devb), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),uvv(:,:),cvv(:,:),densityv(:),dr
    integer, intent(in) :: mtv(:)
    _REAL_ :: fe
    _REAL_ :: fev
    !h2c    : h^2*c/2 contribution
    !tsvv   : t^* = -u +h -c
    !bridge : contribution from the bridge correction
    _REAL_ :: r, hvv, h2c, tsvv, bridge
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
             bridge = hvv - 1.0d0 / this%coeff(ivv) * ( exp( this%coeff(ivv) * tsvv ) - 1.0d0 )
             !if ( tsvv >= 0.0d0 ) then
             !   bridge = 0.0d0
             !else
             !   bridge = hvv - 1.0d0 / ( this%coeff(ivv) + 1.0d0 ) * ( exp( ( this%coeff(ivv) + 1.0d0 ) * tsvv ) - 1.0d0 )
             !end if
    
             h2c = 0.5d0*hvv*hvv - cvv(ir,ivv)
             fev = fev + cnt*densityv(iv1)*densityv(iv2)*(h2c+bridge)
          enddo
       enddo
       fe = fe + fev*r**2
    enddo
    fe = fe * 2.d0*pi*dr
  end function rism1d_devb_freeEnergyPressure_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets the DEVB closure
!!!IN:
!!!   this : DEVB object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_devb_destroy(this)
    implicit none
    type(rism1d_devb), intent(inout) :: this
    if(safemem_dealloc(this%coeff)/=0)&
         call rism_report_error("DEVB closure: could not deallocate coeff")
  end subroutine rism1d_devb_destroy
end module rism1d_devb_c
