! <compile=optimized>
#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Numerical Bridge (NUB) closure
!!!
!!!Precalculated bridge function is provided as an input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module rism1d_nub_c
  use rism_report_c
  use rism1d_potential_c
  use safemem
  !the NUB type
  type rism1d_nub
     ! bridge for every pair
     _REAL_, pointer :: bvv(:,:)=>NULL()
  end type rism1d_nub

  public rism1d_nub_exchem

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initializes the NUB closure
!!!IN:
!!!   this : NUB object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_nub_new(this)
    implicit none
    type(rism1d_nub), intent(inout) :: this
  end subroutine rism1d_nub_new
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initialize NUB closure
!!!IN:
!!!   this : the NUB closure object
!!!   nubfile : .nub filename
!!!   pot : rism1d_potential object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_nub_initialize(this,nubfile,pot)
    use rism_util, only : freeUnit
    implicit none
    type(rism1d_nub), intent(inout) :: this
    character(len=*), intent(in) :: nubfile
    type(rism1d_potential), intent(in) :: pot
    integer :: unit, iostat
    integer, parameter :: linelen=20
    character(len=linelen) :: line
    integer :: ir, ivv

    ! initialize bvv
    this%bvv => safemem_realloc(this%bvv,pot%nr,pot%nvv)
    this%bvv = 0.0

    ! read the numerical bridge function
    call rism_report_message('reading numerical bridge file: '//nubfile)
    unit = freeUnit()
    open(unit=unit,file=nubfile,status='old',iostat=iostat)
    if(iostat /= 0) &
       call rism_report_error('(a,i4)',"Could not open "//nubfile//":",iostat)
    ! skip the leading lines of comments (starting with #)
    do while(.true.)
       read(unit,*,iostat=iostat) line
       if ( iostat /= 0 ) exit
       line = adjustl(line)
       if ( line(1:1) == '#' ) cycle
       backspace(unit)
       exit
    enddo
    ! format of the line
    write(line,'(i10)') pot%nvv
    line = '(18x,'//trim(adjustl(line))//'(x,E16.8E3))'
    ! read numbers
    do ir=1,pot%nr
       read(unit,line,iostat=iostat) (this%bvv(ir,ivv),ivv=1,pot%nvv)
       if ( iostat /= 0 ) exit
    enddo
    close(unit)
  end subroutine rism1d_nub_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Gvv from Uvv, Hvv, and Cvv using the NUB closure
!!!IN:
!!!   this : the NUB closure object
!!!   gvv  : site-site pair correlation function
!!!   uvv  : site-site potential
!!!   hvv  : site-site total correlation function
!!!   cvv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_nub_gvv(this,gvv, uvv, hvv, cvv)
    implicit none
    type(rism1d_nub), intent(in) :: this
    _REAL_, intent(out) :: gvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),hvv(:,:),cvv(:,:)
    integer :: ivv, ir, i
    _REAL_ :: tvv
    do ivv = 1, ubound(gvv,2)
       do ir = 1, ubound(gvv,1)
          tvv = -uvv(ir,ivv) +hvv(ir,ivv)-cvv(ir,ivv) +this%bvv(ir,ivv)
          gvv(ir,ivv) = exp(tvv)
       end do
    end do
  end subroutine rism1d_nub_gvv

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
!!$  subroutine rism1d_nub_gvvdt(this, gvvdt, uvv, gvv, cvv, hvvdt, cvvdt)
!!$    implicit none
!!$    type(rism1d_nub), intent(in) :: this
!!$    _REAL_, intent(out) :: gvvdt(:,:)
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
!!$  end subroutine rism1d_nub_gvvdt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Bvv (bridge function) from Uvv, Gvv, and Cvv using the NUB closure
!!!IN:
!!!   this : the NUB closure object
!!!   uvv  : site-site potential
!!!   gvv  : site-site total correlation function
!!!   cvv  : site-site direct correlation function
!!!OUT:
!!!   bvv  : site-site bridge function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_nub_bvv(this, uvv, gvv, cvv) result(bvv)
    implicit none
    type(rism1d_nub), intent(in) :: this
    _REAL_, pointer :: bvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),gvv(:,:),cvv(:,:)
    nullify(bvv)
    bvv =>safemem_realloc(bvv,ubound(gvv,1),ubound(gvv,2))
    bvv = this%bvv
  end function rism1d_nub_bvv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the excess chemical potential in kT
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site pair correlation function
!!!   uvv  : site-site potential
!!!   cvv  : site-site direct correlation function
!!!   mtv  : multiplicity of each site
!!!   densityv : number density of each site
!!!   densitytrgt : (optional) The final, physical density for thermodynamics.  This
!!!             can be used as an effective correction for some closures that
!!!             over estimate the pressure (e.g. HNC and NUB)
!!!   dr   : r-space grid spacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_nub_exchem(this, gvv, cvv, mtv, jvv, densityv, densitytrgtv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_nub), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:), densityv(:), densitytrgtv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))
    _REAL_ :: r, h2c, hvv, bvv, exchemv
    integer :: ivv, iv1, iv2, ir
    
    do iv2=1,ubound(mtv,1) !nv
       exchem(iv2) = 0.d0
       do ir=2,ubound(gvv,1) !nr
          r = (ir-1)*dr
          exchemv = 0.d0
          do iv1=1,ubound(mtv,1) !nv
             ivv = jvv(iv1,iv2)
             hvv = gvv(ir,ivv) - 1.d0
             bvv = this%bvv(ir,ivv)
             h2c = hvv*( 0.5d0*(hvv-cvv(ir,ivv)) + bvv ) - cvv(ir,ivv) + bvv
             exchemv = exchemv + densitytrgtv(iv1)*mtv(iv2)*h2c
          enddo
          exchem(iv2) = exchem(iv2) + exchemv*r**2
       enddo
    enddo
    exchem = exchem * 4.d0*pi*dr
  end function rism1d_nub_exchem

!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!Calculate the ionic(?) excess chemical potential in kT.  This seems to take
!!$!!!the contribution only from the last solvent type.  I'm not sure what the point is
!!$!!!here.
!!$!!!IN:
!!$!!!   this : the closure object
!!$!!!   gvv  : site-site pair correlation function
!!$!!!   cvv  : site-site direct correlation function
!!$!!!   mtv  : multiplicity of each site
!!$!!!   densityv : number density of each site
!!$!!!   densitytrgt : (optional) The final, physical density for thermodynamics.  This
!!$!!!             can be used as an effective correction for some closures that
!!$!!!             over estimate the pressure (e.g. HNC and NUB)
!!$!!!   dr   : r-space grid spacing
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  function rism1d_nub_exchemIon(this, gvv, uvv, cvv, mtv, jvv, densityv, densitytrgtv, dr) result(exchem)
!!$    use constants, only : pi
!!$    implicit none
!!$    type(rism1d_nub), intent(in) :: this
!!$    _REAL_, intent(in) :: gvv(:,:),uvv(:,:),cvv(:,:), densityv(:), densitytrgtv(:), dr
!!$    integer, intent(in) :: mtv(:),jvv(:,:)
!!$    _REAL_ :: exchem(ubound(mtv,1))
!!$    !hnc    : HNC contribution
!!$    !tsvv   : t^* = u +h -c
!!$    !series : contribution from the series expansion
!!$    _REAL_ :: r, hnc,hvv, tsvv, exchemv,series
!!$    integer :: ivv, iv1, iv2, ir, i
!!$    
!!$    do iv2=1,ubound(mtv,1) !nv
!!$       exchem(iv2) = 0.d0
!!$       do ir=2,ubound(gvv,1) !nr
!!$          r = (ir-1)*dr
!!$          exchemv = 0.d0
!!$!          do iv1=1,ubound(mtv,1) !nv
!!$          iv1=ubound(mtv,1)
!!$             ivv = jvv(iv1,iv2)
!!$             hvv = gvv(ir,ivv) - 1.d0
!!$             tsvv = -uvv(ir,ivv) + hvv -cvv(ir,ivv)
!!$             if (tsvv >= 0.d0)  then
!!$                series = tsvv**(this%order1)/this%order1fac
!!$             else
!!$                series = 0
!!$             endif
!!$             hnc = 0.5d0*hvv*hvv - cvv(ir,ivv)  - 0.5d0*hvv*cvv(ir,ivv)           
!!$             !                     exchemv = exchemv + densityv(iv1)*mtv(iv2) &
!!$             exchemv = exchemv + densitytrgtv(iv2)*mtv(iv1) &
!!$                  *(hnc - series)
!!$ !         enddo
!!$          exchem(iv2) = exchem(iv2) + exchemv*r**2
!!$       enddo
!!$    enddo
!!$    exchem = exchem * 4.d0*pi*dr
!!$  end function rism1d_nub_exchemIon
!!$
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
!!$  function rism1d_nub_solene(this, gvv, uvv, cvv, gvvdt, cvvdt, mtv, jvv, densityv, densitytrgtv, dr) result(solene)
!!$    use constants, only : pi
!!$    implicit none
!!$    type(rism1d_nub), intent(in) :: this
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
!!$  end function rism1d_nub_solene
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
!!$  function rism1d_nub_pressureFE(this, gvv, uvv, cvv, mtv, densityv, dr) result(pr)
!!$    use constants, only : pi
!!$    implicit none
!!$    type(rism1d_nub), intent(in) :: this
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
!!$  end function rism1d_nub_pressureFE
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!Calculate the r-space contribution to the free energy in internal units
!!$!!!IN:
!!$!!!   this : the closure object
!!$!!!   gvv  : site-site pair correlation function
!!$!!!   cvv  : site-site direct correlation function
!!$!!!   mtv  : multiplicity of each site
!!$!!!   densityv : number density of each site
!!$!!!   dr   : r-space grid spacing
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  function rism1d_nub_freeEnergy(this, gvv, uvv, cvv, mtv, densityv, dr) result(fe)
!!$    use constants, only : pi
!!$    implicit none
!!$    type(rism1d_nub), intent(in) :: this
!!$    _REAL_, intent(in) :: gvv(:,:),uvv(:,:),cvv(:,:),densityv(:),dr
!!$    integer, intent(in) :: mtv(:)
!!$    _REAL_ :: fe
!!$    _REAL_ :: fev
!!$    !h2c    : h^2*c/2 contribution
!!$    !tsvv   : t^* = u +h -c
!!$    !series : contribution from the series expansion
!!$    _REAL_ :: r, hvv, h2c, tsvv, series
!!$    integer :: ir,ivv,iv1,iv2, cnt
!!$    !r-space
!!$    fe = 0.d0
!!$    do ir=2,ubound(gvv,1) !nr
!!$       r = (ir-1)*dr
!!$       fev = 0.d0
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
!!$             fev = fev + cnt*densityv(iv1)*densityv(iv2)*(h2c-series)
!!$          enddo
!!$       enddo
!!$       fe = fe + fev*r**2
!!$    enddo
!!$    fe = -fe * 2.d0*pi*dr
!!$  end function rism1d_nub_freeEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets the NUB closure
!!!IN:
!!!   this : NUB object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_nub_destroy(this)
    implicit none
    type(rism1d_nub), intent(inout) :: this
    if(safemem_dealloc(this%bvv)/=0)&
         call rism_report_error("NUB closure: could not deallocate bvv")
  end subroutine rism1d_nub_destroy
end module rism1d_nub_c
