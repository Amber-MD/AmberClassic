! <compile=optimized> 
!  -*- mode: f90; coding: iso-8859-15; -*-

! DFTB3
!
! gam12_dftb3:
! Calculate gamma or gamma^h function and Gamma-function
! for an atom pair
!
! gam121_dftb3:
! Ditto, but gradient
!
! Author: Andreas W. Goetz
! Date  : September, 2016

#include "../include/dprec.fh"

subroutine gam12_dftb3(dist, uhi, uhj, uhdi, xhgamma, zeta, &
     gval, gder)

  implicit none

  _REAL_, intent(in) :: dist, uhi, uhj, uhdi, zeta
  logical, intent(in) :: xhgamma
  _REAL_, intent(out) :: gval, gder

  _REAL_, parameter :: small = 1.0d-04
  _REAL_ :: a, b, ar, g, s, a2, a3, a4, a6, b2, b3, b4, b6
  _REAL_ :: fab, fba, dfabda, dfbada, dgda, dsdu
  _REAL_ :: h, dhdu
  _REAL_ :: ear, ebr

  gval = 0.0d0
  gder = 0.0d0

  ! alpha = 16/5 U_a
  a = 3.2d0 * uhi
  b = 3.2d0 * uhj

  if ((a+b) < small) then
     return
  end if

  if (dist < small) then

     if (abs(a-b) < small) then
        gval = 0.15625d0*(a+b)  ! with a=b: gval = 5/16*a = U
        gder = 0.5d0            ! with a=b: gder = 1/2 (for easier summation that is wanted!)
     else
        gval = 0.5d0*( a*b/(a+b)+a*a*b*b/((a+b)**3) )
        gder = 1.6d0*( b/(a+b) -a*b/(a+b)**2 +2.0d0*a*b*b/(a+b)**3 &
             -3.0d0*a*a*b*b/(a+b)**4 )
     end if

  else

     ! 1/r-S
     if (abs(a-b) < small) then  ! same atom type
        ar   = (a+b)/2.0d0*dist
        ear  = exp(-ar)
        g    = (48.0d0+33.0d0*ar+9.0d0*ar*ar+ar*ar*ar)/(dist*48.0d0)
        s    = ear*g

        dgda = (33.0d0+18.0d0*ar+3.0d0*ar*ar)/48.0d0
        dsdu = 3.2d0*ear*(dgda-dist*g)
     else
        ear = exp(-a*dist)
        ebr = exp(-b*dist)
        a2   = a*a
        a3   = a2*a
        a4   = a2*a2
        a6   = a4*a2
        b2   = b*b
        b3   = b2*b
        b4   = b2*b2
        b6   = b4*b2
        fab  = a*b4/(2.0d0*(a2-b2)**2)-(b6-3.0d0*a2*b4)/((a2-b2)**3*dist)
        fba  = b*a4/(2.0d0*(b2-a2)**2)-(a6-3.0d0*b2*a4)/((b2-a2)**3*dist)
        s    = ear*fab + ebr*fba

        dfabda = -(b6+3.0d0*a2*b4)/(2.0d0*(a2-b2)**3) &
                         -12.0d0*a3*b4/(dist*(a2-b2)**4)
        dfbada = 2.0d0*b3*a3/((b2-a2)**3) &
                         +12.0*a3*b4/(dist*(b2-a2)**4)
        dsdu = 3.2d0*(ear*(dfabda-dist*fab)+ebr*dfbada)
     end if

     ! special case for hydrogen atoms (gamma^h)
     if (xhgamma) then
        h    = exp(-((a+b)*0.15625d0)**zeta*dist*dist)
        dhdu = -h*zeta*dist*dist*((a+b)*0.15625d0)**(zeta-1.0d0)*0.5d0
        gval = 1.0d0/dist - s*h
        gder = -(dsdu*h+s*dhdu)
     else
        gval = 1.0d0/dist - s
        gder = -dsdu
     end if

  end if

  gder = gder*uhdi

end subroutine gam12_dftb3


subroutine gam121_dftb3(r,ui,uj,udi,xhgamma,zeta,dcdr,dcdr3)

  implicit none

  _REAL_, intent(in) :: r, ui, uj, udi, zeta
  logical, intent(in) :: xhgamma
  _REAL_, intent(out) :: dcdr, dcdr3

  ! internal variables
  _REAL_, parameter ::  small = 1.0d-4
  _REAL_ :: r2,a,b,a2,a3,a4,b2,b3,b4
  _REAL_ :: z,z2,zr,g,dgdr,dsdr
  _REAL_ :: fab,fba,dfabdr,dfbadr
  _REAL_ :: dcdudr,dsdudr,dgdadr,dgda,dfabdadr,dfbadadr,dfabda,dfbada
  _REAL_ :: h,dhdu,dhdr,dhdudr,s,dsdu

  r2 = r*r

  ! alpha = 16/5 U_a
  a  = 3.2d0 * ui
  b  = 3.2d0 * uj

  if ( ( (a+b) < small) .or. (r < small) ) then

     dcdr  = 0.0d0
     dcdr3 = 0.0d0

  else ! here 1/r-s

     if (dabs(a-b) < 1.0d-5) then

        z    = 0.5d0*(a+b)
        z2   = z*z
        zr   = z*r
        g    = (48.0d0+33.0d0*zr+9.0d0*zr*zr+zr*zr*zr)/(48.0d0*r)
        dgdr = -1.0d0/r2+3.0d0*z2/16.0d0+z2*zr/24.0d0
        dsdr = dexp(-zr)*(dgdr-z*g)           
  
        dgda   = (33.0d0+18.0d0*zr+3.0d0*zr*zr)/48.0d0
        dgdadr = 0.375d0*z + 0.125d0*z2*r
        dsdudr = 3.2d0*dexp(-zr)*(g*(zr-1.0d0)-z*dgda+dgdadr-r*dgdr)
        if(xhgamma) then
           s    = dexp(-zr)*g
           dsdu = 3.2d0*dexp(-zr)*(dgda-r*g)
        endif
        
     else
        a2  = a*a
        a3  = a2*a
        a4  = a2*a2
        b2  = b*b
        b3  = b2*b
        b4  = b2*b2
        fab = a*b4/(2.0d0*(a2-b2)**2)-(b4*b2-3.0d0*a2*b4) &
                                          /((a2-b2)**3*r)
        fba = b*a4/(2.0d0*(b2-a2)**2)-(a4*a2-3.0d0*b2*a4) &
                                          /((b2-a2)**3*r)
        dfabdr = (b4*b2-3.0d0*a2*b4)/((a2-b2)**3*r2)
        dfbadr = (a4*a2-3.0d0*b2*a4)/((b2-a2)**3*r2)
        dsdr = dexp(-a*r)*(dfabdr-a*fab)+dexp(-b*r)*(dfbadr-b*fba)

        dfabda   = -(b2*b4+3.0d0*a2*b4)/(2.0d0*(a2-b2)**3) &
                     -12.0d0*a3*b4/(r*(a2-b2)**4)
        dfbada   = 2.0d0*b3*a3/((b2-a2)**3) &
                     +12.0d0*a3*b4/(r*(b2-a2)**4)
        dfabdadr =  12.0d0*a3*b4/(r2*(a2-b2)**4)
        dfbadadr = -12.0d0*a3*b4/(r2*(b2-a2)**4)
        dsdudr = 3.2d0*(dexp(-a*r)*(fab*(a*r-1.0d0)-a*dfabda &
                   +dfabdadr-r*dfabdr) +dexp(-b*r)*(dfbadadr-b*dfbada) )
        if(xhgamma) then
           s   =dexp(-a*r)*fab + dexp(-b*r)*fba
           dsdu=3.2d0*(dexp(-a*r)*(dfabda-r*fab)+dexp(-b*r)*dfbada)
        endif
        
     endif

     ! check for gamma^h (special case for hydrogen atoms)
     if(xhgamma) then

        h      = dexp(-((a+b)*0.15625d0)**zeta*r*r)
        dhdu   = -h*zeta*r*r*((a+b)*0.15625d0)**(zeta-1.0d0)*0.5d0
        dhdr   = -h*2.0d0*r*((a+b)*0.15625d0)**zeta
        dhdudr = h*zeta*r*((a+b)*0.15625d0)**(zeta-1.0d0)* &
                   (r*r*((a+b)*0.15625d0)**zeta-1.0d0)
        dcdr   = -1.0d0/r2-(dsdr*h+s*dhdr)
        dcdudr = -(dsdudr*h+dsdu*dhdr+dsdr*dhdu+s*dhdudr)
        dcdr3  = dcdudr * udi
        
     else

        dcdr   = -1.0d0/r2 - dsdr
        dcdudr = -dsdudr
        dcdr3  = dcdudr * udi

     endif
     
  endif ! end 1/r-s

  !AWGDEBUG
  !write(6,*) r,ui,uj,udi,xhgamma,zeta,dcdr,dcdr3
  
end subroutine gam121_dftb3
