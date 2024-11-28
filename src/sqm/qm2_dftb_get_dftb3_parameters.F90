! <compile=optimized> 
!  -*- mode: f90; coding: iso-8859-15; -*-

! DFTB3
! Author: Andreas W. Goetz
! Date  : August, 2016

#include "../include/dprec.fh"

subroutine qm2_dftb_get_dftb3_parameters(silence, natom, ntyp, atyp, atomic_number)

  use qm2_dftb_module, only: uhder  ! Hubbard derivatives dU/dq
  use qm2_dftb_module, only: zeta   ! exponent for gamma^h function
  use ElementOrbitalIndex, only: NumberElements, ElementSymbol

  implicit none

  logical, intent(in) :: silence
  integer, intent(in) :: natom
  integer, intent(in) :: ntyp
  integer, dimension(*), intent(in) :: atyp
  integer, dimension(*), intent(in) :: atomic_number

  _REAL_, dimension(NumberElements) :: uhder_params
  _REAL_, parameter :: tiny = 1.0d-10
  character(2), dimension(ntyp) :: atyp_element
  integer :: i

  ! REFERENCES
  ! A = JCTC 9 (2013) 338
  ! B = JCTC 10 (2014) 1518
  ! C = JPCB 119 (2015) 1062
  ! D = JCTC 11 (2015) 332
  
  uhder_params = 0.0d0
  uhder_params(1)  = -0.1857d0   ! H : A, Table 1
  uhder_params(6)  = -0.1492d0   ! C : -"-
  uhder_params(7)  = -0.1535d0   ! N : -"-
  uhder_params(8)  = -0.1575d0   ! O : -"-
  uhder_params(9)  = -0.1623d0   ! F : D, Table 1
  uhder_params(11) = -0.0454d0   ! Na: D, Table 1
  uhder_params(12) = -0.02d0     ! Mg: C, Table 1
  uhder_params(15) = -0.14d0     ! P : B, Table 1
  uhder_params(16) = -0.11d0     ! S : -"-
  uhder_params(17) = -0.0697d0   ! Cl: D, Table 1
  uhder_params(19) = -0.0339d0   ! K : -"-
  uhder_params(20) = -0.0340d0   ! Ca: -"-
  uhder_params(30) = -0.03d0     ! Zn: C, Table 1
  uhder_params(35) = -0.0573d0   ! Br: D, Table 1
  uhder_params(53) = -0.0433d0   ! I : -"-

  zeta = 4.0d0 ! A, Table 1

  do i = 1, natom
     uhder(atyp(i)) = uhder_params(atomic_number(i))
     atyp_element(atyp(i)) = ElementSymbol(atomic_number(i))
     ! sanity check - quit if no Hubbard derivative exists for this element
     if ( dabs(uhder(atyp(i))) < tiny ) then
        call sander_bomb('qm2_dftb_get_dftb3_parameters',&
             'Hubbard derivative dU/dq not found.', &
             'Element '//trim(atyp_element(atyp(i)))//' not supported.')
     end if
  end do

  if (.not. silence) then
     write(6,*)
     write(6,'(a)') "QMMM: Hubbard Derivatives dU/dq:"
     do i = 1, ntyp
        write(6,'(a,a3,2x,f10.6)') "QMMM:  ", atyp_element(i), uhder(i)
     end do
     write(6,*)
     write(6,'(a,f10.6)') "QMMM: zeta = ", zeta
  end if

end subroutine qm2_dftb_get_dftb3_parameters
