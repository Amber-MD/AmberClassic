#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit the final minimization report 
!-----------------------------------------------------------------------
!     --- REPORT_MIN_RESULTS ---
!-----------------------------------------------------------------------
! Find the maximum component of the gradient (grdmax),
! emit this and other gradient details (printe),
! optionally do pseudocontact shift constraints,
! optionally decompose energies for mm_pbsa,
! and emit nmr related information (nmrptx).

subroutine report_min_results( nstep, gradient_rms, coordinates, &
      forces, energies, igraph, xx, ix, ih )

   use state
   use file_io_dat
   implicit none

   integer, intent(in)             :: nstep
   _REAL_,  intent(in)             :: gradient_rms
   _REAL_,  intent(in)             :: coordinates(*)
   _REAL_,  intent(in)             :: forces(*)
   type(state_rec), intent(in)     :: energies
   character(len=4), intent(in)    :: igraph(*)    ! atom name map
   _REAL_,  intent(inout)          :: xx(*)        ! real dynamic memory
   integer, intent(inout)          :: ix(*)        ! integer dynamic memory
   character(len=4), intent(inout) :: ih(*)        ! hollerith dynamic memory

#  include "box.h"
#  include "extra.h"
#  include "../include/md.h"
#  include "../include/memory.h"
#  include "nmr.h"
#  include "tgtmd.h"
#  include "multitmd.h"

   character(len=4) :: atom_name_of_gmax
   integer          :: atom_number_of_gmax
   _REAL_           :: gradient_max
   _REAL_ emtmd

   if (master) then
      write(6, '(/ /20x,a,/)' ) 'FINAL RESULTS'
      call grdmax( forces, gradient_max, atom_number_of_gmax )
      atom_name_of_gmax = igraph(atom_number_of_gmax)
      if (imin /= 5) rewind(MDINFO_UNIT)
      call printe( nstep, gradient_rms, gradient_max, energies, &
             atom_number_of_gmax, atom_name_of_gmax )
      if (nmropt > 0) then
         if (iredir(7) /= 0) call pcshift(-1,coordinates,forces)
         call nmrptx(6)
         call nmrptx(MDINFO_UNIT)
         call ndvptx(coordinates,forces,ih(m04),ih(m02),ix(i02),nres, &
               xx(l95),natom,xx(lwinv),xx(lnmr01),ix(inmr02),6)
      end if
      if (itgtmd == 2) then
        emtmd = 0.0d0
        call mtmdcall(emtmd,xx(lmtmd01),ix(imtmd02),coordinates,forces,ih(m04),ih(m02),ix(i02),&
                                      ih(m06),xx(lmass),natom,nres,'PRNT')
      end if
      if (imin /= 5) call flush(MDINFO_UNIT)
   end if

   return
end subroutine report_min_results


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit a minimization progress report 
!-----------------------------------------------------------------------
!     --- REPORT_MIN_PROGRESS ---
!-----------------------------------------------------------------------
! Find the maximum component of the gradient (grdmax),
! emit this and other gradient details (printe),
! and emit nmr related information (nmrptx).

subroutine report_min_progress( nstep, gradient_rms, forces, energies, &
                                igraph, charge )

   use state
   use file_io_dat
   implicit none

   integer, intent(in)           :: nstep
   _REAL_,  intent(in)           :: gradient_rms
   _REAL_,  intent(in)           :: forces(*)
   type(state_rec), intent(in)   :: energies
   character(len=4), intent(in)  :: igraph(*)    ! atom name map
   _REAL_, intent(in)            :: charge(*)

#  include "extra.h"
#  include "../include/md.h"
#  include "nmr.h"

   character(len=4) :: atom_name_of_gmax
   integer          :: atom_number_of_gmax
   _REAL_           :: gradient_max

   if (master) then
      call grdmax( forces, gradient_max, atom_number_of_gmax )
      atom_name_of_gmax = igraph(atom_number_of_gmax)
      if (imin /= 5) rewind(MDINFO_UNIT)
      call printe( nstep, gradient_rms, gradient_max, energies, &
             atom_number_of_gmax, atom_name_of_gmax )
      if (nmropt > 0) then
         call nmrptx(6)
         call nmrptx(MDINFO_UNIT)
      end if
      if (imin /= 5) call flush(MDINFO_UNIT)
   end if

   return
end subroutine report_min_progress

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute the maximum gradient component and the corresponding atom
!-----------------------------------------------------------------------
!     --- GRDMAX ---
!-----------------------------------------------------------------------

subroutine grdmax( gradient, gradient_max, atom_number_of_gmax )

   use constants, only : zero
   implicit none
   _REAL_,  intent(in)  :: gradient(*)
   !     This is actually two-dimensional (3,natoms), but to enable
   !     vectorization on IA32 SSE platforms they are treated as
   !     one-dimensional; this may also improve software pipelining !

   _REAL_,  intent(out) :: gradient_max
   integer, intent(out) :: atom_number_of_gmax

#  include "../include/memory.h"

   integer :: i
   _REAL_  :: gi

   gradient_max = ZERO
   atom_number_of_gmax = 1
   do i = 1,3*natom
      gi = abs(gradient(i))
      if (gi > gradient_max) then
         gradient_max = gi
         atom_number_of_gmax = i
      end if
   end do
   atom_number_of_gmax = (atom_number_of_gmax - 1)/3 + 1

   return
end subroutine grdmax


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Print status of a minimization calculation.
!-----------------------------------------------------------------------
!     --- PRINTE ---
!-----------------------------------------------------------------------
! Output the step number, the gradient rms, the max gradient component,
! and the atom label and atom number of the max gradient component.

subroutine printe( nstep, gradient_rms, gradient_max, ene, &
      atom_number_of_gmax, atom_name_of_gmax )
   
   use sander_rism_interface, only : rismprm, RISM_NONE, RISM_FULL, RISM_INTERP,&
        rism_calc_type, rism_thermo_print
   use state ! Access to energy_rec
   use charmm_mod, only : charmm_active
   use emap,only : temap,scemap
   use ff11_mod, only : cmap_active
   use xray_interface_module, only: xray_write_min_state
   use xray_globals_module, only: xray_active

   implicit none
   integer, intent(in)          :: nstep
   _REAL_,  intent(in)          :: gradient_rms
   _REAL_,  intent(in)          :: gradient_max
   type(state_rec), intent(in)  :: ene
   integer, intent(in)          :: atom_number_of_gmax
   character(len=4), intent(in) :: atom_name_of_gmax

#  include "../include/md.h"
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#  include "tgtmd.h"

   _REAL_  epot,enonb,enele,ehbond,ebond,eangle,edihed,enb14,eel14,egb,epb
   _REAL_  econst,epolar,aveper,aveind,avetot,esurf,edisp,diprms,dipiter
   _REAL_  dipole_temp,escf,dvdl,enemap
   _REAL_ :: erism
   _REAL_ :: pot_array(potential_energy_rec_len)
   _REAL_ :: ect

   !     ----- Extract Energies. -----

   epot        = ene%pot%tot
   enonb       = ene%pot%vdw
   enele       = ene%pot%elec
   ehbond      = ene%pot%hbond
   egb         = ene%pot%gb 
   epb         = ene%pot%pb
   ebond       = ene%pot%bond
   eangle      = ene%pot%angle
   edihed      = ene%pot%dihedral
   enb14       = ene%pot%vdw_14
   eel14       = ene%pot%elec_14
   econst      = ene%pot%constraint
   epolar      = ene%pot%polar
   dvdl        = ene%pot%dvdl
   aveper      = ene%aveper
   aveind      = ene%aveind
   avetot      = ene%avetot
   esurf       = ene%pot%surf
   diprms      = ene%diprms
   dipiter     = ene%dipiter
   dipole_temp = ene%dipole_temp
   
   escf    = ene%pot%scf
   edisp   = ene%pot%disp
   enemap   = ene%pot%emap
   erism   = ene%pot%rism
   ect = ene%pot%ct

   write(6,9018)
   write(6,9028) nstep, epot, gradient_rms, gradient_max, &
         atom_name_of_gmax, atom_number_of_gmax
   write(6,9038) ebond,eangle,edihed
! CHARMM SPECIFIC ENERGY TERMS
   if (charmm_active) write(6,9039) ene%pot%angle_ub, &
                                    ene%pot%imp,      &
                                    ene%pot%cmap

   if( igb == 0 .and. ipb == 0 .and. rismprm%rism == 0) then
      write(6,9048) enonb,enele,ehbond
   else if ( igb == 10 .or. ipb /= 0 ) then
      write(6,9050) enonb,enele,epb
   else if(rismprm%rism == 1 )then
      write(6,9051) enonb,enele,erism
   else
      write(6,9049) enonb,enele,egb
   end if
   write(6,9058) enb14,eel14,econst

!  wxw: EMAP energy
   if (temap) then
      write (6,9062) enemap,scemap
   endif

   if( gbsa > 0 ) write(6,9077) esurf
   if (igb == 10 .or. ipb /= 0) write(6,9074) esurf,edisp
      if (cmap_active .and. ipol > 0 ) then
          write(6,9066) epolar, ene%pot%cmap
      else
          if (cmap_active) write(6,9067) ene%pot%cmap
          if (epolar /= 0.0) write(6,9068) epolar
      end if
   if (econst /= 0.0) write(6,9078) epot-econst
   if ( dvdl /= 0.d0) write(6,9100) dvdl
   if (itgtmd == 1) then
      write (6,'(a,f8.3)') "Current RMSD from reference: ",rmsdvalue
      write (6,'(a,f8.3)') "Current target RMSD:         ",tgtrmsd
   end if
   if (xray_active) call xray_write_min_state(6)

   call flush(6)

   !     ----- SEND IT TO THE INFO FILE -----

   if (imin /= 5) then
      write(7,9018)
      write(7,9028) nstep, epot, gradient_rms, gradient_max, &
            atom_name_of_gmax, atom_number_of_gmax
      write(7,9038) ebond,eangle,edihed
! CHARMM SPECIFIC ENERGY TERMS
   if (charmm_active) write(7,9039) ene%pot%angle_ub, &
                                    ene%pot%imp,      &
                                    ene%pot%cmap


      if( igb == 0 .and. ipb == 0 .and. rismprm%rism == 0) then
         write(7,9048) enonb,enele,ehbond
      else if ( igb == 10 .or. ipb /= 0 ) then
         write(7,9050) enonb,enele,epb
      else if ( rismprm%rism == 1 ) then
         write(7,9051) enonb,enele,erism
      else
         write(7,9049) enonb,enele,egb
      end if
      write(7,9058) enb14,eel14,econst

!  wxw: EMAP energy
   if (temap) then
      write (7,9062) enemap,scemap
   endif
   if (xray_active) call xray_write_min_state(7)

      if( gbsa > 0 ) write(7,9077) esurf
      if ( igb == 10 .or. ipb /= 0 ) write(7,9074) esurf,edisp
! FF11 CMAP SPECIFIC ENERGY TERMS
      if (cmap_active .and. epolar /= 0.0 ) then
          write(7,9066) epolar, ene%pot%cmap
      else
          if (cmap_active) write(7,9067) ene%pot%cmap
          if (epolar /= 0.0) write(7,9068) epolar
      end if
      if (econst /= 0.0) write(7,9078) epot-econst
      if ( dvdl /= 0.d0) write(7,9100) dvdl
   end if

   if(rismprm%rism==1 .and. rismprm%write_thermo==1)then
      if(rism_calc_type(nstep) == RISM_FULL)&
           call rism_thermo_print(.false.,transfer(ene%pot,pot_array))
   end if

   9018 format (/ /,3x,'NSTEP',7x,'ENERGY',10x,'RMS',12x,'GMAX',9x, &
         'NAME',4x,'NUMBER')
   9028 format(1x,i6,2x,3(2x,1pe13.4),5x,a4,2x,i7,/)
   9038 format (1x,'BOND    = ',f13.4,2x,'ANGLE   = ',f13.4,2x, &
         'DIHED      = ',f13.4)
   9039 format (1x, 'UB      = ', f13.4, 2x, 'IMP     = ', f13.4, 2x, &
            'CMAP       = ', f13.4)
   9048 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'HBOND      = ',f13.4)
   9049 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'EGB        = ',f13.4)
   9050 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'EPB        = ',f13.4)
   9051 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'ERISM      = ',f13.4)
   9058 format (1x,'1-4 VDW = ',f13.4,2x,'1-4 EEL = ',f13.4,2x, &
         'RESTRAINT  = ',f13.4)
   9062 format (1x,'EMAP   = ',f14.4,   '  EMSCORE = ',f8.4)
   9066 format (1x,'EPOLAR  = ',f13.4,2x,'CMAP    = ',f13.4)
   9067 format (1x,'CMAP    = ',f13.4)
   9068 format (1x,'EPOLAR  = ',f13.4)
   9074 format (1x,'ECAVITY = ',f13.4,2x,'EDISPER = ',f13.4)
   9077 format (1x,'ESURF   = ',f13.4)
   9078 format (1x,'EAMBER  = ',f13.4)
   9080 format (1x,'PM3ESCF =',f14.4)
   9081 format (1x,'AM1ESCF =',f14.4)
   9981 format (1x,'AM1DESCF =',f14.4)
   9082 format (1x,'MNDOESCF=',f14.4)
   9982 format (1x,'MNDODESCF=',f14.4)
   9083 format (1x,'PDDGPM3-ESCF=',f14.4)
   9084 format (1x,'PDDGMNDO-ESCF=',f14.4)
   9085 format (1x,'PM3CARB1-ESCF=',f14.4)
   9086 format (1x,'DFTBESCF=',f14.4)
   9087 format (1x,'RM1ESCF =',f14.4)
   9088 format (1x,'PDDGPM3_08-ESCF=',f14.4)
   9089 format (1x,'PM6ESCF =',f14.4)
   9090 format (1x,'PM3MMXESCF =',f14.4)
   9091 format (1x,'PM3ZNBESCF =',f14.4)
   9092 format (1x,'EXTERNESCF =',f14.4)
   9093 format (1x,'PM3MAISESCF =',f14.4)
   9099 format (1x,'ECRG    =',f14.4)
   9096 format (1x,'PM3MMX2ESCF =',f14.4)
   9190 format(1x,'Dipole convergence: rms = ',e10.3,' iters = ',f6.2)
   9191 format(1x,'Dipole convergence: rms = ',e10.3, &
         ' temperature = ',f6.2)
   9100 format (1x,'DV/DL  = ',f14.4)
   return
end subroutine printe
