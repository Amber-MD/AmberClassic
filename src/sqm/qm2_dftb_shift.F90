! <compile=optimized>
!-*- mode: f90; coding: iso-8859-15; -*-

#include "../include/dprec.fh"
#include "copyright.h"

! DFTB GB implementation by Gustavo Seabra (UFL) and Ross Walker (TSRI), 2005
!
! DFTB3 implementation by Andreas W Goetz (SDSC), 2016

subroutine hamilshift(dftb3, natom, atomtype, gammamat, gamma_der, scf_mchg, &
     shift, shift3, shift3A)

   !=======================================================
   ! get the hubbard contribution to the H matrix elements
   !=======================================================

   use qm2_dftb_module, only : mol, mcharge

   implicit none

   ! Passed in:
   logical, intent(in)  :: dftb3
   integer, intent(in)  :: natom               ! number of atoms
   integer, intent(in)  :: atomtype(natom)     ! 
                                             !              build up gammamat
   _REAL_ , intent(in) :: gammamat(natom,natom)
   _REAL_ , intent(in) :: gamma_der(natom,natom) !AWG DFTB3
   _REAL_ , intent(out) :: scf_mchg(natom)
   _REAL_ , intent(out) :: shift(natom)          ! array contains shifts for hamilton matrix elements
   _REAL_ , intent(out) :: shift3(natom)  ! AWG DFTB3
   _REAL_ , intent(out) :: shift3A(natom) ! AWG DFTB3
   
   ! Locals
   integer :: i,j
   _REAL_ :: tmpvalue, qdiff_i, qdiff_j

  
   scf_mchg(1:natom) = mcharge%qzero(atomtype(1:natom)) &
                                        - mol%qmat(1:natom)

   ! Calculate atomic hamilton shift (=sum over gamma * charges)
   do i=1,natom
      do j=1,natom
         ! gammamat is lower diagonal. All elements where j > i are zero.
         if (j > i) then
            tmpvalue = gammamat(j,i)
         else
            tmpvalue = gammamat(i,j)
         endif
         shift(i) = shift(i) - scf_mchg(j)*tmpvalue
!         write(6,*) 'shift  (',i,') =', shift(i)
      end do
   end do

   if (dftb3) then

      shift3(:) = 0.0d0
      shift3A(:) = 0.0d0
      do i = 1, natom
         qdiff_i = -scf_mchg(i)
         do j = 1, natom
            qdiff_j = -scf_mchg(j)
            shift3(i) = shift3(i) + qdiff_j*gamma_der(i,j)
            shift3A(i) = shift3A(i) + qdiff_j*qdiff_j*gamma_der(j,i)
         end do
         shift3(i) = shift3(i)*qdiff_i
!         write(6,*) 'shift3 (',i,') =', shift3(i)
!         write(6,*) 'shift3A(',i,') =', shift3A(i)
      end do

   end if

end subroutine hamilshift

