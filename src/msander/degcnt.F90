! <compile=optimized>
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Count solute and solvent degrees of freedom.
subroutine degcnt(ibelly,nat,igrp,nsolut,nbonh,nbona,nbper, &
      ibh,jbh,iba,jba,ibp,jbp,ntc,natrst,nbrst, &
      narst,nprst,itorty,rndfp,rndfs &
      )

   ! Subroutine DEGgree CouNT
   
   ! This routine determines how many degrees of freedom can be ascribed
   ! to the "solute", and how many to the "solvent". Unlike any previous
   ! version of AMBER, it works correctly with the belly option, with
   ! separate solute/solvent coupling, etc.
   
   ! It then determines the number of degrees of freedom in the solute and in
   ! the solvent which are lost to SHAKE and TORCON constraints.
   
   ! Returned are three values:
   
   !    IBELSV = number of moving "solvent" atoms
   !    RNDFP = net number of degrees of freedom for the solute
   !    RNDFS = net number of degrees of freedom for the solvent
   
   ! INPUT:
   ! -----
   ! IBELLY > 0 if belly run is being performed.
   ! NAT is the number of atoms in the system.
   ! IGRP(I) is > 0 if atom I is part of moving belly.
   ! NSOLUT is the number of the last solute atom.
   ! NBONH, IBH(), JBH() are the NBONH pointers to bonds to hydrogens.
   ! NBONA, IBA(), JBA() are the NBONA pointers to bonds to non-hydrogens
   ! NBPER, IBP(), JBP() are the NBPER pointers to pert bonds.
   ! NTC is the shake flag (1=no shake; 2,4=shake on bonds to H; 3=shake on
   !       all bonds)
   ! NATRST(4,I) are the 2-4 atoms making up added constraint I.
   ! NBRST is the number of added bond constraints
   ! NARST is the number of added angle constraints
   ! NPRST is the number of added torsion constraints
   ! ITORTY(I) =2 if this is an added _con_straint.

   ! Local:   
   ! RSTSSL = number of degrees of freedom in the solute lost to SHAKE.
   ! RSTSSV = number of degrees of freedom in the solvent lost to SHAKE.
   !         (RSTSSL and RSTSSV omit degrees of freedom already lost in the
   !          belly case because the bond corresponds to two non-moving atoms)
   
   ! Author: David Pearlman
   ! Date: 5/91
   ! Modifications for QMMM by Ross Walker (5/2005)


#ifdef MPI /* SOFT CORE */
   use softcore, only : sc_dof_shaked, ifsc
   use softcore, only : nsc, tishake
#endif /* MPI */
   implicit none

!Passed In
   integer :: ibelly, nat, igrp(*), nsolut, nbonh, nbona, nbper
   integer :: ibh(*), jbh(*), iba(*), jba(*), ibp(*), jbp(*)
   integer :: ntc, natrst(4,*), nbrst, narst, nprst, itorty(*)
   _REAL_ :: rndfp, rndfs

!Local
   integer :: i, ibelsl, ibelsv, ib, jb, inum, j, im, m
   _REAL_ :: rstssl, rstssv, fract
   integer ::  rstssmm, rstssqm
   integer :: ibelsl_qm, ibelsl_mm
   
   ! count up degrees of freedom in the belly, if any.
   
   ibelsl = 0
   ibelsv = 0
   rstssmm = 0
   rstssqm = 0
   ibelsl_qm = 0
   ibelsl_mm = 0

   if (ibelly > 0) then
      do i = 1,nat
         if(igrp(i) > 0) then
            if (i <= nsolut) then
               ibelsl = ibelsl + 1
            else
               ibelsv = ibelsv + 1
            end if
         end if
      end do
   end if

   
   ! now, if shake is one, loop over the appropriate bond. Add up all bonds
   ! in the solute/solvent which cannot move because of shake.
   ! For a belly run, do not count bonds which already cannot move because both
   ! atoms are part of the non-moving section.
   ! If a bond is 1/2 in solute, 1/2 in solvent, assign 1/2 to both counters.
   
   rstssl = 0.0d0
   rstssv = 0.0d0
#ifdef MPI /* SOFT CORE */
   if (ifsc /= 0) sc_dof_shaked = 0
#endif
   
   ! bonds to H
   
   if (ntc >= 2) then
      do i = 1,nbonh
         ib = ibh(i)/3 + 1
         jb = jbh(i)/3 + 1
         if (ibelly <= 0) then
            if (ib <= nsolut .and. jb <= nsolut) then
               rstssl = rstssl + 1.0d0
            else if (ib > nsolut .and. jb > nsolut) then
               rstssv = rstssv + 1.0d0
            else
               rstssl = rstssl + 0.5d0
               rstssv = rstssv + 0.5d0
            end if
# ifdef MPI /* SOFT CORE */
            ! Check if this bond lies in the softcore region
            ! SHAKE is removed from bonds crossing into or out 
            ! of the SC part by setnoshake_sc() in set.f
            if (ifsc /= 0) then
               if (tishake == 1) then
                 if (nsc(ib) == 1 .and. nsc(jb) == 1) sc_dof_shaked = sc_dof_shaked + 1
               else
                 if (nsc(ib) == 1 .or. nsc(jb) == 1) sc_dof_shaked = sc_dof_shaked + 1
               end if
            end if
# endif
         else if (igrp(ib) > 0 .or. igrp(jb) > 0) then
            if (ib <= nsolut .and. jb <= nsolut) then
               rstssl = rstssl + 1.0d0
            else if (ib > nsolut .and. jb > nsolut) then
               rstssv = rstssv + 1.0d0
            else
               rstssl = rstssl + 0.5d0
               rstssv = rstssv + 0.5d0
            end if
         end if  ! (ibelly <= 0)
      end do
   end if  !  20 i = 1,nbonh
   
   ! bonds to heavy atoms
   
   if (ntc == 3) then
      do i = 1,nbona
         ib = iba(i)/3 + 1
         jb = jba(i)/3 + 1
         if (ibelly <= 0) then
            if (ib <= nsolut .and. jb <= nsolut) then
               rstssl = rstssl + 1.0d0
            else if (ib > nsolut .and. jb > nsolut) then
               rstssv = rstssv + 1.0d0
            else
               rstssl = rstssl + 0.5d0
               rstssv = rstssv + 0.5d0
            end if
# ifdef MPI /* SOFT CORE */
            ! Check if this bond lies in the softcore region
            ! SHAKE is removed from bonds crossing into or out 
            ! of the SC part by setnoshake_sc() in set.f if tishake == 1
            if (ifsc /= 0) then
               if (tishake == 1) then
                 if (nsc(ib) == 1 .and. nsc(jb) == 1) sc_dof_shaked = sc_dof_shaked + 1
               else
                 if (nsc(ib) == 1 .or. nsc(jb) == 1) sc_dof_shaked = sc_dof_shaked + 1
               end if
            end if
# endif
         else if (igrp(ib) > 0 .or. igrp(jb) > 0) then
            if (ib <= nsolut .and. jb <= nsolut) then
               rstssl = rstssl + 1.0d0
            else if (ib > nsolut .and. jb > nsolut) then
               rstssv = rstssv + 1.0d0
            else
               rstssl = rstssl + 0.5d0
               rstssv = rstssv + 0.5d0
            end if
         end if  ! (ibelly <= 0)
      end do

      ! bonds to pert atoms
      
      do i = 1,nbper
         ib = ibp(i)/3 + 1
         jb = jbp(i)/3 + 1
         if (ibelly <= 0) then
            if (ib <= nsolut .and. jb <= nsolut) then
               rstssl = rstssl + 1.0d0
            else if (ib > nsolut .and. jb > nsolut) then
               rstssv = rstssv + 1.0d0
            else
               rstssl = rstssl + 0.5d0
               rstssv = rstssv + 0.5d0
            end if
         else if (igrp(ib) > 0 .or. igrp(jb) > 0) then
            if (ib <= nsolut .and. jb <= nsolut) then
               rstssl = rstssl + 1.0d0
            else if (ib > nsolut .and. jb > nsolut) then
               rstssv = rstssv + 1.0d0
            else
               rstssl = rstssl + 0.5d0
               rstssv = rstssv + 0.5d0
            end if
         end if  ! (ibelly <= 0)
      end do ! i = 1,nbper
   end if
   
   ! Added constraints. We ascribe to the solvent and the solute the fraction
   ! of the lost constraint equal to the fraction of the internal from that
   ! source:
   
   do i = 1,nprst
      if (itorty(i) /= 2) cycle
      
      if (i > narst) then
         inum = 4
      else if (i > nbrst) then
         inum = 3
      else
         inum = 2
      end if
      fract = 1.0d0/dble(inum)
      
      if (ibelly <= 0) then
         do j = 1,inum
            im = natrst(j,i)/3 + 1
            if (im <= nsolut) rstssl = rstssl + fract
            if (im > nsolut) rstssv = rstssv + fract
         end do
      else
         do m = 1,inum
            if (igrp(natrst(m,i)) > 0) then
               do j = 1,inum
                  im = natrst(j,i)/3 + 1
                  if (im <= nsolut) rstssl = rstssl + fract
                  if (im > nsolut) rstssv = rstssv + fract
               end do
               exit
            end if
         end do
      end if
   end do
   
   ! Now determine RNDFP (net number of degrees of freedom for solute) and
   !               RNDFS (net number of degrees of freedom for solvent).
   
   if (ibelly <= 0) then
      rndfp = 3*nsolut - rstssl
      rndfs = 3*(nat-nsolut) - rstssv
   else
      rndfp = 3*ibelsl - rstssl
      rndfs = 3*ibelsv - rstssv
   end if
   
   return
end subroutine degcnt 
