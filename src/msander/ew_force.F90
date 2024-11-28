! <compile=optimized>
#include "../include/assert.fh"
#include "../include/dprec.fh"

!----------------------------------------------------------------------
!                   EWALD_FORCE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+
subroutine ewald_force(crd,numatoms,iac,ico,charge, &
      cn1,cn2,cn6,eelt,epol,frc,x,ix,ipairs, virvsene,pol,pol2,qm_pot_only, &
      cn3,cn4,cn5)
   use ew_recip
   use stack
   use nblist, only: recip, cutoffnb,nbfilter, volume, &
                     maxnblst,nucgrd1,nucgrd2,nucgrd3, &
                     adjust_imagcrds, map_coords, &
                     nvdwcls
   use constants, only : INV_AMBER_ELECTROSTATIC2, INV_AMBER_ELECTROSTATIC, &
                         zero, AMBER_ELECTROSTATIC
   use file_io_dat
#ifdef LES
   use les_data, only : les_vir, cnum, nlesadj, eeles
#else
   use charmm_mod, only : charmm_active, charmm_cn114, charmm_cn214
   use parms, only : one_scee, one_scnb
#endif /* LES */
#ifdef MPI
   use mpi
#endif

   implicit none
#  include "extra.h"
#  include "ew_cntrl.h"
#  include "def_time.h"
#  include "ew_time.h"
#  include "ew_pme_recip.h"
#  include "ew_mpole.h"
#  include "ew_erfc_spline.h"
#  include "../include/memory.h"
#  include "../include/md.h"
#  include "extra_pts.h"
#  include "box.h"
#ifndef LES
   ! Dummy vars that are not protected by LES ifdefs
   _REAL_ les_vir(3,3), eeles
#endif
#  include "ew_frc.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   integer ierr
#endif

   integer numatoms,iac(*),ico(*)
   _REAL_ crd(3,*),charge(*),cn1(*),cn2(*),cn6(*), &
          eelt,epol,frc(3,*),virvsene
   _REAL_ pol(*), pol2(*)
   _REAL_ cn3(*), cn4(*), cn5(*)
!!
   logical, intent(in) :: qm_pot_only ! Flag for whether we do a regular 
                                      ! full pme (gradients as well) or just 
                                      ! get the potential for qm scf procedure.
#ifdef MPI
   integer ndel
   _REAL_ rl_temp(BC_EW_COMM3)
#endif

   integer nstart,ntop
   _REAL_ x(*)
   integer ix(*),ipairs(*)
   integer mm_no,qm_no, lnk_no
   _REAL_ forcemod(3)
   _REAL_ boltz2
   ! Intra Gaussian polarization energy
   _REAL_ epolg
#ifdef LES
   integer icopy
#endif

   !     ARGUMENTS:
   !     CRD (INPUT) is the array holding atomic coords.
   !     NUMATOMS (INPUT) is the number of atoms
   !     IAC, ICO and NTYPES (INPUT) are used to look up vdw and hbond
   !       interactions.  ICO is an NTYPES by NTYPES array giving lookup
   !       into VDW and HBOND coefficients. IAC(i) is the atom type of
   !       atom i. The coefficients for 6-12 interactions for atoms i
   !       and j are indexed by ICO(iac(i),iac(j)).  In practice ICO
   !       is unrolled to a 1 dimensional array.  They are needed here
   !       to split nonbond list into vdw,hbonds (see pack_nb_list())
   !     CHARGE (INPUT) is the array of atomic partial charges.
   !     CN1 and CN2 (INPUT) are the VDW coefficients
   !     ASOL and BSOL (INPUT) are the HBOND coefficients
   !     EELT,EVDW,EHB (OUTPUT) are the energies computed and
   !     FRC (OUTPUT) is the force array
   !     X and IX are real and integer arrays which comprise the total
   !       "dynamic" memory  for amber; coords,forces,bond lists etc are
   !       accessed as offsets in them.
   
   integer i,commsander_mytaskid,commsander_numtasks, qm_temp_count

#ifdef MPI
   commsander_mytaskid = mytaskid
   commsander_numtasks = numtasks
#else
   commsander_mytaskid = 0
   commsander_numtasks = 1
#endif
 call timer_start(TIME_EWFSTRT)
   
   ! Clear all energies and virials plus our scratch force array
   
   eelt = 0.d0
   eer = 0.d0
   eed = 0.d0
   eea = 0.d0
   ees = 0.d0
   epol = 0.d0
   epold = 0.d0
   epola = 0.d0
   epols = 0.d0
   evdw = 0.d0
   ehb = 0.d0
   eedvir = 0.d0
   evdwr = 0.d0
   !     if( numextra.gt.0 ) then
   ee14 = 0.d0
   enb14 = 0.d0
   epol14 = 0.d0

   ! Intra Gaussian polarization energy
   epolg = 0.d0

   !     endif
   eeles=0.d0
   rec_vir(1:3,1:3) = 0.d0
   rec_vird(1:3,1:3) = 0.d0
   dir_vir(1:3,1:3) = 0.d0
   adj_vir(1:3,1:3) = 0.d0
   self_vir(1:3,1:3) = 0.d0
   les_vir(1:3,1:3) = 0.d0
   e14vir(1:3,1:3) = 0.d0
   framevir(1:3,1:3) = 0.d0
   virvsene = 0.d0

   frcx(1:3) = 0.d0
   
   call timer_stop(TIME_EWFSTRT)
   if( numextra > 0 ) call local_to_global(crd,x,ix)
   
   ! Set use_pme=2 to always run PME, e.g. when rismrespa is set as well)
   if( use_pme == 2 .or. (use_pme == 1 .and. mod(irespa,nrespa) == 0)) then
      
      !--------------------------------------------------------
      ! SELF ENERGY
      !--------------------------------------------------------
      
      call timer_start(TIME_SELF)
      if ( master ) call self(charge,numatoms,ees,ew_coeff,volume,self_vir)
      call timer_stop(TIME_SELF)

      !--------------------------------------------------------
      ! RECIPROCAL ENERGY
      !--------------------------------------------------------


      call timer_start(TIME_REC)
#ifdef MPI
        if(i_do_recip)then
           call mpi_comm_size(recip_comm,numtasks,ierr)
           call mpi_comm_rank(recip_comm,mytaskid,ierr)
           master = mytaskid.eq.0
#endif

           if ( ew_type == 0 )then
              call do_pme_recip(mpoltype,numatoms,crd,charge,frc,     &
                    X(linddip),x(lfield),x(lprefac1),x(lprefac2),    &
                    X(lprefac3),x(lfftable),qm_pot_only )
           else
              ASSERT( .false. )  ! ew_type input validation occurs in mdread.f
           end if
        
           ! Scaleup forces, field for respa
         
           if ( use_pme == 1 .and. nrespa > 1 ) &
               call respa_scale(numatoms,frc,nrespa)
#ifdef MPI
           numtasks = commsander_numtasks
           mytaskid = commsander_mytaskid
           master = mytaskid.eq.0
        end if         !loop for ido_recip
#endif
      call timer_stop(TIME_REC)
      
      !-------------------------------------------------------
      ! LONG RANGE DISPERSION CONTRIBUTION
      !-------------------------------------------------------

      call timer_start(TIME_SELF)
      if ( vdwmeth == 1)then
         if ( master ) then
               call vdw_correction(ico,ntypes,nvdwcls, &
               volume,evdwr,rec_vird,cn2,cutoffnb)
         end if
      end if
      call timer_stop(TIME_SELF)
      
   end if      ! (use_pme and respa check)

   !--------------------------------------------------------
   ! DIRECT PART OF EWALD PLUS VDW, HBOND
   !--------------------------------------------------------
   
   call timer_start(TIME_DIR)
#ifdef MPI
   if(i_do_direct)then
      call mpi_comm_size(direct_comm,numtasks,ierr)
      call mpi_comm_rank(direct_comm,mytaskid,ierr)
#endif

      !         Get fractional coordinates, and adjust the imaged cartesians
      !         for pre-imaged short range interactions

      call map_coords(crd,numatoms,recip)
      call adjust_imagcrds(crd,natom )

      call get_nb_energy(iac,ico,ntypes,charge, &
            cn1,cn2,cn6,frc,numatoms, ipairs, &
            ew_coeff,eedtbdns,x(leed_cub),x(leed_lin), &
            maxnblst,eed,evdw,ehb,dir_vir,eedvir, &
            nbfilter,ee_type,eedmeth,dxdr, &
            cn3, cn4, cn5, epold,x(linddip),x(lfield))

#ifdef MPI
       numtasks = commsander_numtasks
       mytaskid = commsander_mytaskid
   end if
#endif
   call timer_stop(TIME_DIR)

   !--------------------------------------------------------
   ! ADJUST ENERGIES,FORCES FOR MASKED OUT PAIRS
   !--------------------------------------------------------
   if( use_pme /= 0 ) then
      call timer_start(TIME_ADJ)

      call nb_adjust(charge,eea,crd, &
            ix(imask1),ix(imask2),numadjst, &
            ew_coeff,eedtbdns, &
            x(leed_cub),x(leed_lin),frc,numatoms,adj_vir, &
            ee_type,eedmeth)

      call timer_stop(TIME_ADJ)
   end if

#ifdef LES
   
   ! ADJUST ENERGIES,FORCES FOR LES INTRA-COPY PAIR SCALING FACTOR
   
   if (nlesadj > 0) then
      call timer_start(TIME_LESADJ)
      call nb_adjust_les(charge,eeles,crd, frc, &
            numatoms,les_vir,use_pme)
      call timer_stop(TIME_LESADJ)
   end if
  
#endif /* LES */
   
#ifndef LES
   !--------------------------------------------------------
   ! do the 1-4 here & remove them from ephi
   !--------------------------------------------------------

   if (charmm_active) then
     call get_14_cg(charge,crd,frc,iac, &
           charmm_cn114,charmm_cn214, ee14,enb14,one_scee,one_scnb, &
           e14vir,ix,commsander_mytaskid,commsander_numtasks,eedmeth)
   ! mjhsieh: note that ff11cmap doesn't need to touch this
   else
     call get_14_cg(charge,crd,frc,iac, &
           cn1,cn2, ee14,enb14,one_scee,one_scnb, &
           e14vir,ix,commsander_mytaskid,commsander_numtasks,eedmeth)
   end if

   ! Now transfer force and torque from extra points:
   if (numextra > 0) call orient_frc(crd,frc,framevir,ix)
#endif

   !--------------------------------------------------------
   ! ADJUST FORCES FOR net force (comes from recip work)
   !--------------------------------------------------------

   call timer_start(TIME_EWFRCADJ)
#ifdef MPI
   ndel=(numatoms+numtasks-1)/numtasks
   nstart = ndel*mytaskid+1
   ntop = nstart+ndel-1
   if(ntop > numatoms)ntop=numatoms
   call mpi_allreduce(MPI_IN_PLACE,frcx,3,MPI_DOUBLE_PRECISION, &
         mpi_sum,commsander,ierr)

#else

   !  ---- Put force back in frc array from temp. force array ----

   nstart=1
   ntop=numatoms
#endif /* MPI */

   frcx(1:3)=frcx(1:3)/dble(numatoms - numextra)
   if ( verbose > 0 )then
      if ( master )write(6,33)frcx(1),frcx(2),frcx(3)
      33 format(1x,'NET FORCE PER ATOM: ',3(1x,e12.4))
   end if
   !     ---do not remove net force if netfrc = 0; e.g. in minimization
   
   if ( netfrc > 0 )then
      do i = nstart,ntop

#ifdef LES
#else
         frc(1,i) = frc(1,i) - frcx(1)
         frc(2,i) = frc(2,i) - frcx(2)
         frc(3,i) = frc(3,i) - frcx(3)
#endif
      end do
#ifndef LES
      ! Re-zero any extra points forces that got whacked above...
      if (numextra > 0) call zero_extra_pnts_vec(frc, ix)
#endif
   end if
   call timer_stop_start(TIME_EWFRCADJ,TIME_EWVIRIAL)   

#ifdef MPI
   
   !     Accumulate the virials and energies
   call mpi_allreduce(MPI_IN_PLACE,eer,BC_EW_COMM3, &
         MPI_DOUBLE_PRECISION, mpi_sum,commsander,ierr)

#  ifdef LES
   call mpi_allreduce(MPI_IN_PLACE,eeles,10,MPI_DOUBLE_PRECISION, &
         mpi_sum,commsander,ierr)
#  endif
#endif /* MPI */
   
#ifdef LES
   eelt = ees + eer + eed + eea + eeles
#else
   eelt = ees + eer + eed + eea
#endif
   epol = epold + epola + epols 
   epol = epol + epol14
   
   evdw = evdw + evdwr
   
   call timer_stop(TIME_EWVIRIAL)  

   return
end subroutine ewald_force 

#ifdef LES
!--------------------------------------------------------------------
subroutine assign_copy_charge(numatoms,charge,chg_single,id)
   use les_data, only : cnum
   implicit none
   integer numatoms,id
#  include "../include/memory.h"
   _REAL_  charge(*), chg_single(*)
   integer i

   do i=1,numatoms
      if( cnum(i).eq.0 ) then
         chg_single(i) = charge(i)
      else if( cnum(i).eq.id ) then
         chg_single(i) = charge(i)*ncopy
      else 
         chg_single(i) = 0.0
      end if
    end do
end subroutine
#endif

!--------------------------------------------------------------------
!         DO_PME_RECIP
!--------------------------------------------------------------------
subroutine do_pme_recip(mpoltype,numatoms,crd,charge,frc,dipole,   &
      field,prefac1,prefac2,prefac3,fftable,qm_pot_only)
   use ew_recip
   use nblist, only: recip, volume
#ifdef MPI
   use mpi
#endif
   implicit none
#  include "../include/memory.h"

   integer mpoltype, numatoms
   _REAL_ crd(3,numatoms),charge(numatoms),frc(3,numatoms)
   _REAL_ dipole(3,numatoms),field(3,numatoms)
   logical,intent(in) :: qm_pot_only

#  include "ew_frc.h"
#  include "../include/md.h"
#  include "ew_pme_recip.h"

#ifdef MPI
#  include "parallel.h"
#endif
   ! OUTPUT
   !       eer:  ewald reciprocal or k-space  energy
   !       frc forces incremented by k-space sum
   !       virial:  virial due to k-space sum (valid for atomic scaling;
   !                rigid molecule virial needs a correction term not
   !                computed here
   
   ! HEAP STORAGE:  These arrays need to be preserved throughout simulation
   _REAL_ prefac1(*),prefac2(*),prefac3(*),fftable(*)

#ifdef LES
   integer icopy
   _REAL_ eer_sum
#endif

#ifdef LES      
   call do_pmesh_kspace(numatoms,crd,charge,frc, &
          prefac1,prefac2,prefac3,fftable,qm_pot_only)
          frcx(:) = frcx(:) * dble(nrespa) ! scale up for respa

#else /* LES */

     !  column_fft stuff is not working, November, 2021
     ! if(column_fft_flag)then
     !    call spatial_do_pmesh_kspace(numatoms,crd,charge,frc, &
     !         prefac1,prefac2,prefac3,qm_pot_only)
     ! else
        call do_pmesh_kspace(numatoms,crd,charge,frc, &
            prefac1,prefac2,prefac3,fftable,qm_pot_only)
     ! endif
     frcx(:) = frcx(:) * dble(nrespa) ! scale up for respa

#endif /* LES */
end subroutine do_pme_recip

!     --- NB_ADJUST ---
!     ...the part of ewald due to gaussian counterion about an atom
!     you are bonded to or otherwise for which you do not compute the
!     nonbond pair force. NECESSARY since you ARE computing this pair
!     in the reciprocal sum...

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nb_adjust here]
subroutine nb_adjust(charge,eea,crd, &
      mask1,mask2, &
      numadjst,ewaldcof,eedtbdns, &
      eed_cub,eed_lin,frc,numatoms, &
      adj_vir, ee_type,eedmeth)

   use constants, only : INVSQRTPI, third, half
   use file_io_dat
#ifdef LES
   use les_data, only : lfac, lestyp, lesfac, cnum, nlesty
   use nblist, only : cutoffnb
#endif
#ifdef MPI
   use mpi
#endif
   implicit none

   _REAL_ charge(*),eea,crd(3,*)
   integer numadjst,numatoms
   integer mask1(*),mask2(*)
   integer ee_type,eedmeth
   _REAL_ ewaldcof,eedtbdns,eed_cub(4,*),eed_lin(2,*)
   _REAL_ adj_vir(3,3),frc(3,numatoms)

#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
#endif

   integer numlo,numhi
#ifdef MPI
   integer numdel,numleft
#endif

   integer i,k,nx,ind
   _REAL_ r2,r,x,df,cgi,cgk,erfc,derfc
   _REAL_ dfx,dfy,dfz
   _REAL_ del,dx,xx,delr2inv
   _REAL_ delx,dely,delz
   _REAL_ vxx,vxy,vxz,vyy,vyz,vzz
   _REAL_ d0,d1
#include "box.h"
#include "../include/md.h"

#ifdef LES
   _REAL_ xbig, ecur
   _REAL_ xsmall,rsmall,esmall,ezero,efit, &
         fsmall,fzero,ffit,prefac
   
   ! for LES, we do NOT use the scale factor for nb_adjust.
   ! this is because we are correcting for the interactions calculated
   ! in the reciprocal sum that should not have been included.
   ! so, we want to subtract exactly what we calculated- the non-scaled value.
   ! we do not scale the direct sum part since we CANNOT scale
   ! the reciprocal component.
   ! for part PIMD, we do use the scale factor for nb_adjust.
   ! since we have calculate the reciprocal energy correctly
   ! by call ew_recip P times
#endif


#  include "flocntrl.h"
   ! FLOW CONTROL FLAG (debug, possible vacuum sims and future respa)
   if ( do_adj == 0 )return

   !
   !     initialize virial
   !
   vxx = 0.d0
   vxy = 0.d0
   vxz = 0.d0
   vyy = 0.d0
   vyz = 0.d0
   vzz = 0.d0
   !
   !     return if no terms to adjust!
   !
   if ( numadjst .eq. 0) return

   
   !     Get numlo, numhi for MPI and 1processor version
   
#ifdef MPI
   
   !     numdel is the amount of numadjst that each PE gets
   
   numdel = numadjst / numtasks
   if ( numdel == 0 )numdel = 1
   
   !--MFC-- numleft is the amout of numadjst left after numadjst is distributed
   !   with numdel to each PE. It is the number of PEs that
   !   need an extra numadjst.
   
   numleft = mod(numadjst,numtasks)
   
   !--MFC-- PE0 always gets an extra (inddel + 1) unless ind is evenly
   !   divisible by number of PEs
   
   if(mytaskid == 0 )then
      if( numleft == 0)then
         numlo = 1
         numhi = numdel
      else
         numlo = 1
         numhi = numdel + 1
      end if
      
      !--MFC-- Other PEs get either numdel or (numdel + 1) depending on
      !   whether all the numleft is used up.
      
   else
      if(mytaskid < numleft) then
         numlo = mytaskid*(numdel+1)+1
         numhi = (numdel+1) * (mytaskid+1)
      else
         numlo = numleft*(numdel+1) + &
               (mytaskid - numleft)*(numdel) +1
         numhi = (numlo - 1) + numdel
      end if
   end if
#else    /*   MPI   */
   numlo = 1
   numhi = numadjst
#endif
#ifdef LES
   
   !     Carlos/Tom D. modifications for LES
   !     this is to account for possible short and long distances
   !     between atoms that are excluded
   
   prefac = 2.d0 * INVSQRTPI
   xsmall = 0.01d0
   rsmall = xsmall / ewaldcof
   call erfcfun(xsmall,erfc)
   derfc = prefac*exp(-xsmall**2)
   esmall = (erfc - 1.d0)/rsmall
   ezero = -ewaldcof * prefac
   
   fsmall = (erfc - 1.d0)/rsmall + ewaldcof*derfc
   fsmall = fsmall/(rsmall*rsmall)

   fzero = 0.d0

#endif /* LES */
   vxx = 0.d0
   vxy = 0.d0
   vxz = 0.d0
   vyy = 0.d0
   vyz = 0.d0
   vzz = 0.d0
   eea = 0.d0
   del = 1.d0 / eedtbdns
   do nx = numlo,numhi
      i = mask1(nx)
      k = mask2(nx)


      cgi = charge(i)
      cgk = charge(k)
      delx = crd(1,k) - crd(1,i)
      dely = crd(2,k) - crd(2,i)
      delz = crd(3,k) - crd(3,i)
      r2 = delx*delx+dely*dely+delz*delz
      r=sqrt(r2)
#ifdef TEST_CLOSEATOMS
      if(r < .8) &
            write(6,'("<NB_adjust> Atoms close:",2i8,e10.3)')i,k,r
#endif
      x = ewaldcof*r
#ifdef LES
      xbig=cutoffnb*ewaldcof
      
      !       special cases below for LES
      !       xbig is for excluded atoms that are far apart in space
      !       xsmall is for copies that may be overlapping
      
      if ( x > xbig )then
         
         !         put 1/r2 here so it won't cause NaN in xsmall routine
         
         r2=1.d0/r2
         eea = eea - cgi*cgk/r
         df = -cgi*cgk*r2/r
         dfx = delx*df
         dfy = dely*df
         dfz = delz*df
      else if (x <= xsmall)then
         efit = (x/xsmall)*esmall + ((xsmall - x)/xsmall)*ezero
         ffit = (x/xsmall)*fsmall + ((xsmall - x)/xsmall)*fzero
         eea = eea + cgi*cgk*efit
         
         !         *r2 accounted for in fsmall calc above
         
         !         df = cgi*cgk*(ffit)*r2
         df = cgi*cgk*(ffit)
         dfx = delx*df
         dfy = dely*df
         dfz = delz*df
      else

#endif /* LES */
         
         ! similar code to that in short_ene; however the only valid option
         !         here is the erfc switch, and dxdr = ewaldcof.
         
         if ( eedmeth == 1 )then
            
            !         ---cubic spline on erfc derfc
            
            ind = int(eedtbdns*x) + 1
            dx = x - (ind-1)*del
            erfc = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
                  dx*(eed_cub(3,ind)+dx*eed_cub(4,ind)*third)*half)
            derfc = (eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
                  dx*eed_cub(4,ind)*half))
         else if ( eedmeth == 2 )then
            
            !         ---linear lookup on erfc, deriv
            
            xx = eedtbdns*x + 1
            ind = int(xx)
            dx = xx - ind
            erfc = (1.d0 - dx)*eed_lin(1,ind) + &
                  dx*eed_lin(1,ind+1)
            derfc = (1.d0 - dx)*eed_lin(2,ind) + &
                  dx*eed_lin(2,ind+1)
         else if ( eedmeth == 3 )then
            call get_ee_func(x,erfc,derfc,ee_type)
         end if
         delr2inv = 1.d0 / r2
         d0 = (erfc - 1.d0)/r
         d1 = (d0 - ewaldcof*derfc)*delr2inv
#ifdef LES
         eea = eea + cgi*cgk*d0
         df = cgi*cgk*d1
#else
         eea = eea + cgi*cgk*d0
         df = cgi*cgk*d1
#endif /* LES */
         dfx = delx*df
         dfy = dely*df
         dfz = delz*df
#ifdef LES
         
         !       endif for (x.gt.xbig) above
         
      end if  ! ( x > xbig )
#endif
      vxx = vxx - dfx*delx
      vxy = vxy - dfx*dely
      vxz = vxz - dfx*delz
      vyy = vyy - dfy*dely
      vyz = vyz - dfy*delz
      vzz = vzz - dfz*delz
      frc(1,k) = frc(1,k) + dfx
      frc(2,k) = frc(2,k) + dfy
      frc(3,k) = frc(3,k) + dfz
      frc(1,i) = frc(1,i) - dfx
      frc(2,i) = frc(2,i) - dfy
      frc(3,i) = frc(3,i) - dfz
   end do  !  nx = numlo,numhi
   adj_vir(1,1) = vxx
   adj_vir(1,2) = vxy
   adj_vir(2,1) = vxy
   adj_vir(1,3) = vxz
   adj_vir(3,1) = vxz
   adj_vir(2,2) = vyy
   adj_vir(2,3) = vyz
   adj_vir(3,2) = vyz
   adj_vir(3,3) = vzz
   return
end subroutine nb_adjust 
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nb_adjust_dipole here]
subroutine nb_adjust_dipole(charge,eea,crd, &
      mask1,mask2,numadjst,ewaldcof,eedtbdns, &
      eed_cub,eed_lin,frc,numatoms, &
      adj_vir,ee_type,eedmeth, &
      dipole,field,epola)
   use constants, only : third, half
#ifdef MPI
   use mpi
#endif
   implicit none

   _REAL_ charge(*),eea,crd(3,*)
   integer mask1(*),mask2(*),numadjst,numatoms
   integer ee_type,eedmeth
   _REAL_ ewaldcof,eedtbdns,eed_cub(4,*),eed_lin(2,*)
   _REAL_ adj_vir(3,3),frc(3,numatoms), &
         dipole(3,*),field(3,*),epola

#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
#endif

   integer numlo,numhi
#ifdef MPI
   integer numdel,numleft
#endif

   integer i,j,k,nx,ind
   _REAL_ r2,r,x,cgi,cgk,erfc,derfc
   _REAL_ dfx,dfy,dfz
   _REAL_ del,dx,xx,delr2inv
   _REAL_ delx,dely,delz
   _REAL_ d0,d1,d2,d3,fac,dotir,dotkr,dotik,dotp
   _REAL_ term,term0,term1,termi,termk,cgp
   _REAL_ dphii_dx,dphii_dy,dphii_dz, &
         dphik_dx,dphik_dy,dphik_dz

#  include "flocntrl.h"
   !     FLOW CONTROL FLAG (debug, possible vacuum sims and future respa)
   if ( do_adj == 0 )return
   
   !     Get numlo, numhi for MPI and 1processor version
   
#ifdef MPI
   
   !     numdel is the amount of numadjst that each PE gets
   
   numdel = numadjst / numtasks
   if ( numdel == 0 )numdel = 1

   !--MFC-- numleft is the amout of numadjst left after numadjst is distributed
   !   with numdel to each PE. It is the number of PEs that
   !   need an extra numadjst.
   
   numleft = mod(numadjst,numtasks)
   
   !--MFC-- PE0 always gets an extra (inddel + 1) unless ind is evenly
   !   divisible bu number of PEs
   
   if(mytaskid == 0 )then
      if( numleft == 0)then
         numlo = 1
         numhi = numdel
      else
         numlo = 1
         numhi = numdel + 1
      end if
      
      !--MFC-- Other PEs get either numdel or (numdel + 1) depending on
      !   whether all the numleft is used up.
      
   else
      if(mytaskid < numleft) then
         numlo = mytaskid*(numdel+1)+1
         numhi = (numdel+1) * (mytaskid+1)
      else
         numlo = numleft*(numdel+1) + &
               (mytaskid - numleft)*(numdel) +1
         numhi = (numlo - 1) + numdel
      end if
   end if
#else
   numlo = 1
   numhi = numadjst
#endif
   do j = 1,3
      do i = 1,3
         adj_vir(i,j) = 0.d0
      end do
   end do
   eea = 0.d0
   epola = 0.d0
   del = 1.d0 / eedtbdns
   fac = 2.d0*ewaldcof*ewaldcof
   do nx = numlo,numhi
      i = mask1(nx)
      k = mask2(nx)
      cgi = charge(i)
      cgk = charge(k)
      delx = crd(1,k) - crd(1,i)
      dely = crd(2,k) - crd(2,i)
      delz = crd(3,k) - crd(3,i)
      r2 = delx*delx+dely*dely+delz*delz
      r=sqrt(r2)
      x = ewaldcof*r
      
      !       similar code to that in short_ene; however the only valid option
      !       here is the erfc switch, and dxdr = ewaldcof.
      
      if ( eedmeth == 1 )then
         
         !         ---cubic spline on erfc derfc
         
         ind = int(eedtbdns*x) + 1
         dx = x - (ind-1)*del
         erfc = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
               dx*(eed_cub(3,ind)+dx*eed_cub(4,ind)*third)*half)
         derfc = (eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
               dx*eed_cub(4,ind)*half))
      else if ( eedmeth == 2 )then
         
         !         ---linear lookup on erfc, deriv
         
         xx = eedtbdns*x + 1
         ind = int(xx)
         dx = xx - ind
         erfc = (1.d0 - dx)*eed_lin(1,ind) + &
               dx*eed_lin(1,ind+1)
         derfc = (1.d0 - dx)*eed_lin(2,ind) + &
               dx*eed_lin(2,ind+1)
      else if ( eedmeth == 3 )then
         call get_ee_func(x,erfc,derfc,ee_type)
      end if
      delr2inv = 1.d0 / r2
      d0 = (erfc - 1.d0)/r
      d1 = (d0 - ewaldcof*derfc)*delr2inv
      d2 = (3.d0*d1 - fac*ewaldcof*derfc)*delr2inv
      d3 = (5.d0*d2 - fac*fac*ewaldcof*derfc)*delr2inv
      dotkr = dipole(1,k)*delx+dipole(2,k)*dely+dipole(3,k)*delz
      dotir = dipole(1,i)*delx+dipole(2,i)*dely+dipole(3,i)*delz
      dotp = dotkr*dotir
      cgp = cgi*cgk
      dotik = dipole(1,i)*dipole(1,k)+dipole(2,i)*dipole(2,k)+ &
            dipole(3,i)*dipole(3,k)
      term = cgk*dotir-cgi*dotkr+dotik
      term0 = cgp*d0 + term*d1 - dotp*d2
      term1 = cgp*d1 + term*d2 - dotp*d3
      termi = cgi*d1+dotir*d2
      termk = cgk*d1 - dotkr*d2
      eea = eea + cgp*d0
      epola = epola + term*d1 - dotp*d2
      dfx = term1*delx + termi*dipole(1,k) - termk*dipole(1,i)
      dfy = term1*dely + termi*dipole(2,k) - termk*dipole(2,i)
      dfz = term1*delz + termi*dipole(3,k) - termk*dipole(3,i)
      adj_vir(1,1) = adj_vir(1,1) - dfx*delx
      adj_vir(1,2) = adj_vir(1,2) - dfx*dely
      adj_vir(1,3) = adj_vir(1,3) - dfx*delz
      adj_vir(2,1) = adj_vir(2,1) - dfy*delx
      adj_vir(2,2) = adj_vir(2,2) - dfy*dely
      adj_vir(2,3) = adj_vir(2,3) - dfy*delz
      adj_vir(3,1) = adj_vir(3,1) - dfz*delx
      adj_vir(3,2) = adj_vir(3,2) - dfz*dely
      adj_vir(3,3) = adj_vir(3,3) - dfz*delz
      frc(1,k) = frc(1,k) + dfx
      frc(2,k) = frc(2,k) + dfy
      frc(3,k) = frc(3,k) + dfz
      frc(1,i) = frc(1,i) - dfx
      frc(2,i) = frc(2,i) - dfy
      frc(3,i) = frc(3,i) - dfz
      dphii_dx = termk*delx + d1*dipole(1,k)
      dphii_dy = termk*dely + d1*dipole(2,k)
      dphii_dz = termk*delz + d1*dipole(3,k)
      dphik_dx = -termi*delx + d1*dipole(1,i)
      dphik_dy = -termi*dely + d1*dipole(2,i)
      dphik_dz = -termi*delz + d1*dipole(3,i)
      field(1,i) = field(1,i) - dphii_dx
      field(2,i) = field(2,i) - dphii_dy
      field(3,i) = field(3,i) - dphii_dz
      field(1,k) = field(1,k) - dphik_dx
      field(2,k) = field(2,k) - dphik_dy
      field(3,k) = field(3,k) - dphik_dz
   end do  !  nx = numlo,numhi
   return
end subroutine nb_adjust_dipole 
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine self here]
subroutine self(cg,numatoms,ene,ewaldcof,volume,self_vir)
   use constants, only : PI, INVSQRTPI
   use linear_response, only: ilrt
   use file_io_dat
#ifdef LES
   use les_data, only : lfac, lesfac, lestyp, nlesty, cnum
#endif

   implicit none
#  include "../include/memory.h"
   integer numatoms
   _REAL_ cg(*),ene,ewaldcof,ee_plasma,volume,self_vir(3,3)
   _REAL_ sumq,sumq2
   integer i,j
   _REAL_ factor, d0, d1
   logical :: first = .true.
#  include "flocntrl.h"
#  include "../include/md.h"

   save sumq,sumq2
   data sumq,sumq2/0.d0,0.d0/
   
   if (do_self == 0 .and. ilrt == 0 ) return
   !     ---only compute sumq and sumq2 at beginning. They don't change
   ! if LIE module is in use, this is called repeatedly with different charge arrays, so recompute

   d0 = -ewaldcof*INVSQRTPI
   
   if (first) then
      sumq = 0.d0
      sumq2 = 0.d0
      factor = sqrt(pi/(ewaldcof*ewaldcof*volume))*sqrt(0.5d0)
      do i = 1,numatoms
         sumq = sumq + cg(i)
#ifdef LES
         sumq2 = sumq2 + cg(i)*cg(i)
#else
         sumq2 = sumq2 + cg(i)*cg(i)
#endif
      end do
   end if
   ene = sumq2*d0
   factor = pi / (ewaldcof*ewaldcof*volume)
   ee_plasma = -0.5d0*factor*sumq*sumq
   ! force related to ee_plasma
   ene = ene + ee_plasma
   do j = 1,3
      do i = 1,3
         self_vir(i,j) = 0.d0
      end do
      self_vir(j,j) = -ee_plasma
   end do

   return

   entry ewald_self_reset
   first = .true.

end subroutine self 
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine self_dipole here]
subroutine self_dipole(numatoms,epol,ewaldcof,efield,dipole, &
      iproc,nproc)
   use constants, only : SQRTPI
   implicit none
   integer numatoms,iproc,nproc
   _REAL_ efield(3,*),dipole(3,*),epol,ewaldcof

   _REAL_ fac_dip
   integer i
#  include "flocntrl.h"

   if ( do_self == 0 )return
   epol = 0.d0
   fac_dip = 4.d0*ewaldcof**3/(3.d0*SQRTPI)
   do i = 1,numatoms
      if ( iproc == mod(i,nproc) )then
         epol = epol - fac_dip* &
               (dipole(1,i)**2+dipole(2,i)**2+dipole(3,i)**2)
         efield(1,i) = efield(1,i) + fac_dip*dipole(1,i)
         efield(2,i) = efield(2,i) + fac_dip*dipole(2,i)
         efield(3,i) = efield(3,i) + fac_dip*dipole(3,i)
      end if
   end do
   epol = 0.5d0*epol
   return
end subroutine self_dipole 
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dip_field_corr here]
subroutine dip_field_corr(efield,dipole,rec_vir,numatoms)
   implicit none
   _REAL_ efield(3,*),dipole(3,*),rec_vir(3,3)
   integer numatoms !,iproc,nproc
   
   ! The virial needs to be adjusted when dipoles are present
   ! the structure factors contribute a term, which is equal to
   ! the tensor product of the dipole and the reciprocal sum electric field
   ! this routine computes this term
   
   integer i,j,n
   do n = 1,numatoms
      !      if ( iproc .eq. mod(n,nproc) )then
      do j = 1,3
         do i = 1,3
            rec_vir(i,j) = rec_vir(i,j) + efield(i,n)*dipole(j,n)
         end do
      end do
      !      endif
   end do

   return
end subroutine dip_field_corr 
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine respa_scale here]
subroutine respa_scale(numatoms,frc,nrespa)
   implicit none
   integer numatoms,nrespa
   _REAL_ frc(3,*)
   integer n
   do n = 1, numatoms
      frc(1,n) = nrespa*frc(1,n)
      frc(2,n) = nrespa*frc(2,n)
      frc(3,n) = nrespa*frc(3,n)
      !       field(1,n) = nrespa*field(1,n)
      !       field(2,n) = nrespa*field(2,n)
      !       field(3,n) = nrespa*field(3,n)
   end do
   return
end subroutine respa_scale 
!-------------------------------------------------------------------

#ifdef LES

!-------------------------------------------------------------------
!     --- NB_ADJUST_LES ---

!     Author: Carlos Simmerling

!     The electrostatic energies have not been corrected for the
!     intra-copy scaling factors. Here we go through all pairs
!     that need correcting (calculated at program start) and
!     correct the energy. The calculated energy was cgi*cgk/r,
!     and the proper LES energy is cgi*cgj/r*lfac where lfac
!     corresponds to the number of copies of that pair interaction.
!     The PME electrostatic energy was therefore underestimated by
!     the standard code.

!     This routine is based on the nb_adjust routine.



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nb_adjust_les here]
subroutine nb_adjust_les(charge,ene,crd, &
      frc,numatoms,adj_vir,use_pme)

   use les_data, only : eeles, lfac, lesfac, lestyp, nlesty, nlesadj, ileslst, &
                        jleslst
#ifdef MPI
   use mpi
#endif
   implicit none
   _REAL_ charge(*),ene,crd(3,*)
   integer numatoms,use_pme
   _REAL_ frc(3,numatoms)
   _REAL_ adj_vir(3,3)
   integer numlo,numhi


#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#  ifdef MPI_DOUBLE_PRECISION
#     undef MPI_DOUBLE_PRECISION
#  endif
   integer numleft,numdel
#endif /* MPI */
   !-------------------------------------------------------------------
   _REAL_ tmpreal

   integer i,k,nx
   _REAL_ r2,r,df,cgi,cgk
   _REAL_ dfx,dfy,dfz
   _REAL_ delx,dely,delz
   _REAL_ vxx,vxy,vxz,vyy,vyz,vzz
   
#  include "box.h"
#  include "flocntrl.h"

   !     FLOW CONTROL FLAG (debug, possible vacuum sims and future respa)
   if ( do_adj == 0 )return
   
   ! if we are not doing PME, then return because the only calculation
   ! was for the direct sum and therefore does not need correcting here
   ! since it was done at the time the direct sum was calculated
   
   if (use_pme == 0) return

#ifdef MPI
   
   ! numdel is the amount of nlesadj that each PE gets
   
   numdel = nlesadj / numtasks
   if ( numdel == 0 )numdel = 1
   
   !--MFC-- numleft is the amout of nlesadj left after nlesadj is distributed
   !   with numdel to each PE. It is the number of PEs that
   !   need an extra nlesadj.
   
   numleft = mod(nlesadj,numtasks)

   !--MFC-- PE0 always gets an extra (inddel + 1) unless ind is evenly
   !   divisible bu number of PEs
   
   if(mytaskid == 0 )then
      if( numleft == 0)then
         numlo = 1
         numhi = numdel
      else
         numlo = 1
         numhi = numdel + 1
      end if
      !--MFC-- Other PEs get either numdel or (numdel + 1) depending on
      !   whether all the numleft is used up.
      
   else
      if(mytaskid < numleft) then
         numlo = mytaskid*(numdel+1)+1
         numhi = (numdel+1) * (mytaskid+1)
      else
         numlo = numleft*(numdel+1) + &
               (mytaskid - numleft)*(numdel) +1
         numhi = (numlo - 1) + numdel
      end if
   end if
#else
   numlo = 1
   numhi = nlesadj
#endif
   vxx = 0.d0
   vxy = 0.d0
   vxz = 0.d0
   vyy = 0.d0
   vyz = 0.d0
   vzz = 0.d0
   eeles = 0.d0
   do nx = numlo,numhi
      i = ileslst(nx)
      k = jleslst(nx)
      cgi = charge(i)
      cgk = charge(k)
      lfac=lesfac((lestyp(i)-1)*nlesty+lestyp(k))
      
      !           should check here to make sure lfac is not 1.0! if it
      !           _is_ 1.0, then there is no need for the correction, and
      !           the list is incorrect.
      
      delx = crd(1,k) - crd(1,i)
      dely = crd(2,k) - crd(2,i)
      delz = crd(3,k) - crd(3,i)

      !       r is really 1/r
      
      r = 1.d0/sqrt(delx*delx+dely*dely+delz*delz)
      r2 = r*r
      
      !       (lfac-1) adds the correct term and subtracts the incorrect
      !       one already added
      
      !       routine made more efficient with a
      !       temp variable for cgi*cgk*(lfac - 1.d0)*r
      
      tmpreal=cgi*cgk*(lfac - 1.d0)*r
      
      eeles = eeles + tmpreal
      
      !       have to correct the force in the same way, adding
      !       correct force and subtracting incorrect one
      
      df = tmpreal*r2
      dfx = delx*df
      dfy = dely*df
      dfz = delz*df
      vxx = vxx - dfx*delx
      vxy = vxy - dfx*dely
      vxz = vxz - dfx*delz
      vyy = vyy - dfy*dely
      vyz = vyz - dfy*delz
      vzz = vzz - dfz*delz
      frc(1,k) = frc(1,k) + dfx
      frc(2,k) = frc(2,k) + dfy
      frc(3,k) = frc(3,k) + dfz
      frc(1,i) = frc(1,i) - dfx
      frc(2,i) = frc(2,i) - dfy
      frc(3,i) = frc(3,i) - dfz
   end do  !  nx = numlo,numhi
   ene = eeles
   adj_vir(1,1) = vxx
   adj_vir(1,2) = vxy
   adj_vir(2,1) = vxy
   adj_vir(1,3) = vxz
   adj_vir(3,1) = vxz
   adj_vir(2,2) = vyy
   adj_vir(2,3) = vyz
   adj_vir(3,2) = vyz
   adj_vir(3,3) = vzz
   return
end subroutine nb_adjust_les 
!-------------------------------------------------------------------
#endif /*LES*/
