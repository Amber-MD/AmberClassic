#ifdef MPI
#include "../include/dprec.fh"

!     The AMBER/MPI implementation and support routines were
!     originally and independently implemented and contributed
!     by James Vincent (JV) 7/94.  Modified by tec3, dac and JV.


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ broadcast data from master to all other nodes, at beginning
subroutine startup(xx,ix,ih)
   !************************************************************
   !     Send data needed by all nodes once at startup from master
   !     after master has read in all data
   !************************************************************

   use parms, only : bcast_parms
   use charmm_mod, only : mpi_bcast_charmm_params
   use ff11_mod, only : mpi_bcast_cmap_params
#ifdef LES
   use les_data, only : nlesty, lesfac, BC_LESR, BC_LESI
#endif
   use nblist, only : ucell,bc_ewucr,bc_ewuci,nbflag, &
                     BC_DIRPARS,numnptrs
   use file_io_dat
   use md_scheme, only: ntt, gamma_ln
! SOFT CORE
   use softcore, only : ifsc, scalpha, scbeta, scmask, dynlmb, &
                       sceeorder, tishake
! end SOFT CORE
   use mbar, only : ifmbar, bar_intervall, bar_l_min, bar_l_max, bar_l_incr
   use linear_response, only : ilrt, lrt_interval, lrtmask
! SGLD
   use sgld, only : isgld, isgsta,isgend,nsgsize, &
                    sgft,sgff,sgfg,tsgavg,tsgavp
! EMAP parameters
   use emap, only : temap,gammamap,nemap,nrigid
! bcast variables from mdfil.F90
   use commandline_module, only : commandline_bcast, cpein_specified

   use mpi
   implicit none
#  include "parallel.h"
#  include "ew_parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   integer ierr
#  include "extra.h"
#  include "../include/md.h"
#  include "../include/memory.h"
#  include "nmr.h"
#  include "box.h"
#  include "extra_pts.h"
#  include "ew_pme_recip.h"
#  include "ew_mpole.h"
#  include "ew_erfc_spline.h"
#  include "ew_cntrl.h"
#  include "flocntrl.h"
#  include "debug.h"
#  include "new_time.h"
#  include "tgtmd.h"

   _REAL_ xx(*)
   integer ix(*), ier
   character(len=4) ih(*)

   !     Send and receive common blocks from the master node:

   !  file_io_dat

   call mpi_bcast(ntpr, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(ntwr, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(ntwx, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(ntwv, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(ntwf, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(ntwe, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(ntpp, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(ioutfm, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(ntwprt, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(ntave, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(iredir, 9, MPI_INTEGER, 0, commsander, ierr)

   !  nmr.h:

   call mpi_bcast(nmropt,7,MPI_INTEGER,0,commsander,ierr)   ! /nmr1/
   call mpi_bcast(intreq,6,MPI_INTEGER,0,commsander,ierr)   ! /nmrstf/
   call mpi_bcast(wnoesy,6,MPI_DOUBLE_PRECISION,0,commsander,ierr) ! /wremar/
   call mpi_bcast(scalm,7,MPI_DOUBLE_PRECISION,0,commsander,ierr) ! /nmr1/
   call mpi_bcast(dobsu,BC_ALIGNR,MPI_DOUBLE_PRECISION,0,commsander,ierr)
                                                                  ! /align/
   call mpi_bcast(ndip, BC_ALIGNI,MPI_INTEGER,0,commsander,ierr)  ! /align/
   call mpi_bcast(nath,BC_METHYLI,MPI_INTEGER,0,commsander,ierr)  ! /methyli/
   call mpi_bcast(tau,BC_METHYLR,MPI_DOUBLE_PRECISION,0,commsander,ierr)
                                                                  ! /methylr/

   !  md.h:

   call mpi_bcast(nrp,BC_MDI,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(t,BC_MDR,MPI_DOUBLE_PRECISION,0,commsander,ierr)

   !  box.h:

   call mpi_bcast(ntb,BC_BOXI,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(box,BC_BOXR,MPI_DOUBLE_PRECISION,0,commsander,ierr)

   !  parms.f:
   call bcast_parms()

#ifdef LES
   !   les.h:
   call mpi_bcast(lesfac,BC_LESR,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(nlesty,BC_LESI,MPI_INTEGER,0,commsander,ierr)
#endif

   call mpi_bcast(ntt, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(gamma_ln, 1, MPI_DOUBLE_PRECISION, 0, commsander, ierr)

   ! carlos: targeted MD

   call mpi_bcast(itgtmd,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(tgtrmsd,2,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   ! end targeted md

   !  ew_pme_recip.h:

   call mpi_bcast(sizfftab,BC_PME_PARS_INT,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(dsum_tol,BC_PME_PARS_REAL,MPI_DOUBLE_PRECISION, &
         0,commsander,ier)

   ! ew_mpole.h

   call mpi_bcast(ifirst,BC_MULTPOLE,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(diptol,BC_INDDIPR,MPI_DOUBLE_PRECISION,0,commsander,ier)
   call mpi_bcast(maxiter,BC_INDDIPI,MPI_INTEGER,0,commsander,ier)

   !  erfc_spline.h

   call mpi_bcast(leed_cub,6,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(eedtbdns,2,MPI_DOUBLE_PRECISION,0,commsander,ier)



   !---- from nonbond_list.f module nblist -----------------------
   !     common/dirpars/
   call mpi_bcast(numnptrs,BC_DIRPARS,MPI_INTEGER,0,commsander,ier)
   ! was in  ew_unitcell.h now in module nblist
   call mpi_bcast(ucell,BC_EWUCR,MPI_DOUBLE_PRECISION,0,commsander,ier)
   call mpi_bcast(nbflag,BC_EWUCI,MPI_INTEGER,0,commsander,ier)


   !  ewcntrl.h

   call mpi_bcast(verbose,BC_EWCNTRL,MPI_INTEGER,0,commsander,ier)
   inocutoff=0
   inogrdptrs=0
   if(master)then
      if(nocutoff)inocutoff=1
      if(nogrdptrs)inogrdptrs=1
   end if
   call mpi_bcast(inogrdptrs,BC_EWCNTRL_NP,MPI_INTEGER, &
         0,commsander,ier)
   nogrdptrs=inogrdptrs == 1
   nocutoff=inocutoff == 1

   ! flocntrl.h

   call mpi_bcast(do_dir,BC_FLOCNTRL,MPI_INTEGER,0,commsander,ier)

   ! debug.h

   call mpi_bcast(do_debugf,BC_DEBUG,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(lscg,BC_DEB_HEAP,MPI_INTEGER,0,commsander,ier)

   ! timer info new_time.h

   call mpi_bcast(tpar_p,BC_TIME_PAR,MPI_INTEGER,0,commsander,ier)

   ! extra points

   call mpi_bcast(ifrtyp,BC_EXTRA_PT,MPI_INTEGER,0,commsander,ier)

   !  ew_parallel.h

   call mpi_bcast(indz,BC_SLABS,MPI_INTEGER,0,commsander,ier)


   ! sgld

   call mpi_bcast(isgld,1,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(isgsta,1,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(isgend,1,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(nsgsize,1,MPI_INTEGER,0,commsander,ier)
   call mpi_bcast(tsgavg,1,MPI_DOUBLE_PRECISION,0,commsander,ier)
   call mpi_bcast(tsgavp,1,MPI_DOUBLE_PRECISION,0,commsander,ier)
   call mpi_bcast(sgft,1,MPI_DOUBLE_PRECISION,0,commsander,ier)
   call mpi_bcast(sgff,1,MPI_DOUBLE_PRECISION,0,commsander,ier)
   call mpi_bcast(sgfg,1,MPI_DOUBLE_PRECISION,0,commsander,ier)

   !     IX,XX,IH

   call mpi_bcast(ix(1),lasti,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(ih(1),4*lasth,MPI_CHARACTER,0,commsander,ierr)
   call mpi_bcast(xx(1),lastr,MPI_DOUBLE_PRECISION,0,commsander,ierr)

! SOFT CORE
   ! Get all nodes informed about the SC parameters
   call mpi_bcast(ifsc,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(scalpha,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(scbeta,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(sceeorder,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(scmask,256,MPI_CHARACTER,0,commsander,ierr)
   call mpi_bcast(dynlmb,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(tishake,1,MPI_INTEGER,0,commsander,ierr)
! end SOFT CORE
   call mpi_bcast(ifmbar,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(bar_intervall,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(bar_l_min,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(bar_l_max,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(bar_l_incr,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)

   ! Parameters for LIE module
   call mpi_bcast(ilrt,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(lrtmask,256,MPI_CHARACTER,0,commsander,ierr)
   call mpi_bcast(lrt_interval,1,MPI_INTEGER,0,commsander,ierr)

! EMAP
   call mpi_bcast(temap,1,MPI_LOGICAL,0,commsander,ierr)
   call mpi_bcast(NEMAP,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(NRIGID,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(GAMMAMAP,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
! end EMAP

! broadcast commandline info
   call commandline_bcast(ierr)

   ! Charmm force field support - will simply return if charmm is not in use.
   call mpi_bcast_charmm_params(master)
   ! FF11 Cmap support - will simply return if ff11cmap is not in use.
   call mpi_bcast_cmap_params(master)

   call mpi_barrier(commsander,ierr)

   !   ---- divide atoms up among the processors, always splitting on
   !        residue boundaries:

   call setpar(nspm, ix(i70), ntp, ntc, ix(i02), xx(lmass))

   if ( no_ntt3_sync == 1 ) then
     !Here we are not synching the random number generator across threads for
     !NTT=3 but we don't want every thread to have the same random seed. So for
     !the moment just add mytask id to ig.
     ig = ig+mytaskid
   end if

   return
end subroutine startup
!----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fdist here]
subroutine fdist(f,forcetmp,pot,vir,newbalance,size)

   use state
   use mpi
   implicit none
#  include "../include/memory.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   integer ierr
#  include "../include/md.h"
#  include "extra.h"
#  include "nmr.h"

   !     Parameters:

   integer, intent(in) ::  size
   _REAL_, intent(inout)           :: f(size),vir(4)
   _REAL_                          :: forcetmp(size)
   type(potential_energy_rec), intent(inout)  :: pot
   type(potential_energy_rec)                 :: pot_tmp
   integer, intent(out) ::  newbalance

   !     Local:

   integer i,j

   if (numtasks == 1) return

   if( mpi_orig .or. ievb>0 .or. icfe>0 ) then

      !  ---Reduce the force array and energies back to the master node;
      !      (hence, the master will know all coordinates and all forces):

      if (master) then
        call mpi_reduce(MPI_IN_PLACE,f,size, &
            MPI_DOUBLE_PRECISION,mpi_sum,0,commsander,ierr)
        call mpi_reduce(MPI_IN_PLACE,pot,potential_energy_rec_len, &
            MPI_DOUBLE_PRECISION,mpi_sum,0,commsander,ierr)
      else
        call mpi_reduce(f,0,size, &
            MPI_DOUBLE_PRECISION,mpi_sum,0,commsander,ierr)
        call mpi_reduce(pot,0,potential_energy_rec_len, &
            MPI_DOUBLE_PRECISION,mpi_sum,0,commsander,ierr)
      end if

   else

      ! ---Add all copies of energy and put result back on ALL nodes:

      newbalance=0
      !if(f(j+32) > 0.d0)newbalance=1 !MJW TO FIX

      call mpi_allreduce(MPI_IN_PLACE,pot,potential_energy_rec_len, &
            MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
      newbalance=0
      !if(forcetmp(j+32) > 0.d0)newbalance=1  !TODO mjw

      if (init /= 3) then
         !  ---Do a distributed sum of the force array:
         call fsum(f,forcetmp)
      else

         !  Due to lack of parallelization in the initial parts
         !    of runmd in the init=3 case, the more efficient
         !    reduce_scatter needs to be replaced with an mpi_allreduce call;
         !    this is also required for NEB simulations so that all the forces
         !    are available for dot products and removing forces along path
         !    tangent.
         !    in the future we might consider a respa-type NEB to overcome this
         !    extra work
         !   CARLOS: ONLY MASTER NEEDS FORCES FOR NEB
         !   AWG: adaptive QM/MM needs the complete force arrays for force mixing

         call mpi_allreduce(MPI_IN_PLACE, f, 3*natom, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)

      end if

   end if  ! mpi

   return
end subroutine fdist

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fsum here]
subroutine fsum(f,tmp)

   !     equivalent to an MPI_REDUCE_SCATTER on f:  all processors contribute
   !       to f, and the appropriate part of the result winds up on each
   !       processor

   use mpi
   implicit none
   _REAL_ f(*),tmp(*)

#  include "parallel.h"
#  include "extra.h"
#  include "../include/memory.h"
#ifdef MPI_DOUBLE_PRECISION
#  undef MPI_DOUBLE_PRECISION
#endif
   integer ierr

   !Used for Binary Tree
   integer other,ncyclesm1,k,bit,cs,cr,ns,nr,istart,iend
   integer ist(mpi_status_size)

   if (numtasks <= 1) return

   ! If we have a power of two cpus do a binary tree:

   if (logtwo(numtasks)>=1) then  !We have a power of two cpus
     ncyclesm1 = logtwo(numtasks) - 1
     bit = ishft(numtasks,-1)

     do k = 0,ncyclesm1

        other=ieor(mytaskid,bit)

        !        send chunk:

        cs = ishft(other,-((ncyclesm1)-k))*bit
        ns = iparpt3(cs+bit)-iparpt3(cs)

        !        recv chunk:

        cr = ishft(mytaskid,-((ncyclesm1)-k))*bit
        istart = iparpt3(cr)
        iend = iparpt3(cr+bit)
        nr = iend-istart

        call mpi_sendrecv( &
              f(iparpt3(cs)+1),ns,MPI_DOUBLE_PRECISION,other,5, &
              tmp(iparpt3(cr)+1),nr,MPI_DOUBLE_PRECISION,other,5, &
              commsander, ist, ierr )
        f(istart+1:iend) = f(istart+1:iend) + tmp(istart+1:iend)

        bit = ishft(bit,-1)

     end do  !  k = 0,ncyclesm1

   else

     ! We don't have a power of two - do things the old fashioned way.
#ifdef RED_SCAT_INPLACE
     call mpi_reduce_scatter(f, f(iparpt3(mytaskid)+1), &
           rcvcnt3, MPI_DOUBLE_PRECISION, mpi_sum, &
           commsander, ierr)
#else
     ! Reduce scaller used to work with most MPI installations even if the send
     ! and receive buffers aliased each other but now it seems that newer
     ! mpi installations are checking this and quiting with an error. The
     ! standard does not allow this in MPI v1.0. For now we will make the 
     ! non-inplace version the default.
     call mpi_reduce_scatter(f, tmp(iparpt3(mytaskid)+1), &
           rcvcnt3, MPI_DOUBLE_PRECISION, mpi_sum, &
           commsander, ierr)
     f(iparpt3(mytaskid)+1:iparpt3(mytaskid+1)) = tmp(iparpt3(mytaskid)+1:iparpt3(mytaskid+1))
#endif
   end if !Power of two cpus.
   return
end subroutine fsum

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Distribute the coordinates to all processors.
subroutine xdist(x, tmp, size)

   use mpi
   implicit none

   integer, intent(in) :: size   ! will be 3*natom + iscale
   _REAL_ x(size), tmp(size)

#  include "parallel.h"
#    ifdef MPI_DOUBLE_PRECISION
#      undef MPI_DOUBLE_PRECISION
#    endif
   integer ierr

   !Used for Binary Tree
   integer other,ncyclesm1,k,bit,cs,cr,ns,nr
   integer ist(mpi_status_size)

   if (numtasks <= 1) return

   ! If we have a power of two cpus do a binary tree:

   if (logtwo(numtasks)>=1) then
     ncyclesm1 = logtwo(numtasks) - 1
     bit=1
     do k = 0,ncyclesm1
      other=ieor(mytaskid,bit)
      cs = ishft(mytaskid,-k)*bit
      cr = ishft(other,-k)*bit
      ns = iparpt3(cs+bit)-iparpt3(cs)
      nr = iparpt3(cr+bit)-iparpt3(cr)
      call mpi_sendrecv( &
            x(iparpt3(cs)+1),ns,MPI_DOUBLE_PRECISION,other,5, &
            x(iparpt3(cr)+1),nr,MPI_DOUBLE_PRECISION,other,5, &
            commsander, ist, ierr )
      bit = ishft(bit,1)
     end do

   else

     ! We don't have a power of two - do things the old fashioned way.
     call mpi_allgatherv( &
           x(iparpt3(mytaskid)+1),rcvcnt3(mytaskid), &
           MPI_DOUBLE_PRECISION,tmp,rcvcnt3,iparpt3, &
           MPI_DOUBLE_PRECISION,commsander, ierr)
     x(1:size) = tmp(1:size)
   endif
   return
end subroutine xdist

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fgblsum here]
subroutine fgblsum(x,tmp, size)
   implicit none
   integer, intent(in) :: size
   _REAL_ x(size),tmp(size)

   call fsum(x,tmp)
   call xdist(x, tmp, size)

   return
end subroutine fgblsum

#else
subroutine dummy_parallel()
end subroutine dummy_parallel
#endif  /* MPI  */
