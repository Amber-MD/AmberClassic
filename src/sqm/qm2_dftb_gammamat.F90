! <compile=optimized>
!-*- mode: f90; coding: iso-8859-15; -*-

#include "../include/dprec.fh"

subroutine gammamatrix(dftb3,natom,qm_coords,atomtype,ishydrogen, &
     uhubb,uhder,zeta,gammamat,gamma_der)
!===========================================================================
! Build lower triangular Gamma matrix containing short range terms
!
!  INPUT Parameter:
!   INTEGER natom          number of atoms
!   REAL*8 qm_coords(3,*)      position of atoms
!   REAL*8 u(*)          hubbard parameters
!                                                                                                                
!  OUTPUT:
!   REAL*8 gammamat(*,*) matrix containing the values of the ewlad potential
!                       in the upper triangular part
!   !!! NOTE THAT phi(ri - rj) = phi(rj - ri) !!!
!   !!! NOTE THAT shortrange(ri - rj) = shortrange(rj - ri) !!!
!
! DFTB3: Calculate gamma_der (Gamma matrix, derivative of gamma)
! Implementation by Andreas W Goetz (SDSC), August/September 2016
!
!============================================================================

  use constants, only: A_TO_BOHRS

  implicit none

!Passed in:
  logical, intent(in)  :: dftb3
  integer, intent(in)  :: natom              ! number of atoms
  _REAL_ , intent(in)  :: qm_coords(3,natom) ! position of atoms
  integer, intent(in)  :: atomtype(natom)    ! atom types
  logical, intent(in)  :: ishydrogen(*)      ! whether atom type is hydrogen
  _REAL_ , intent(in)  :: uhubb(*)           ! hubbard parameters U
  _REAL_ , intent(in)  :: uhder(*)           ! hubbard derivatives dU/dq
  _REAL_ , intent(in)  :: zeta               ! exponent for gamma^h
  _REAL_ , intent(out) :: gammamat(natom,natom)  ! matrix containing kernel of electron-electron interaction
  _REAL_ , intent(out) :: gamma_der(natom,natom)  ! DFTB3: Gamma = dgamma/dq

!Locals
  integer :: i, j, ati, atj
  logical :: xhgamma
  _REAL_  :: r(3)
  _REAL_  :: gval, gder, norm
  external GAM12


  if (dftb3) then
     
     ! code path for DFTB3
     do i=1,natom
        ati = atomtype(i)
        do j=1,natom
           atj = atomtype(j)
           
           ! Distance between the 2 atoms
           r(1:3)=(qm_coords(1:3,i)-qm_coords(1:3,j))*A_TO_BOHRS
           norm   = sqrt(r(1)**2+r(2)**2+r(3)**2)

           xhgamma = ishydrogen(ati) .or. ishydrogen(atj)

           gval = 0.0d0
           gder = 0.0d0
           
           ! get gamma/gamma^h and Gamma
           call gam12_dftb3(norm, uhubb(ati), uhubb(atj), &
                uhder(ati), xhgamma, zeta, &
                gval, gder)
           
           gammamat(i,j) = gval
           gamma_der(i,j) = gder

!           write(6,*) 'gamma(',i,',',j,')=', gammamat(i,j)
!           write(6,*) 'Gamma(',i,',',j,')=', gamma_der(i,j)
     
        end do
     end do

  else

     ! DFTB2
     do i=1,natom
        do j=1,i
     
           ! Distance between the 2 atoms
           r(1:3)=(qm_coords(1:3,i)-qm_coords(1:3,j))*A_TO_BOHRS
           norm   = sqrt(r(1)**2+r(2)**2+r(3)**2)
     
           gval = 0.0d0
           
           ! get value for gamma
           call GAM12(norm,uhubb(atomtype(i)),uhubb(atomtype(j)),gval)
           gammamat(i,j)=gval
!           write(6,*) 'gamma(',i,',',j,')=', gammamat(i,j)
     
        end do
     end do

  end if

end subroutine gammamatrix

