subroutine getftype(natom,nres,igraph,ipres,lbres,ib,jb,nbond,ftype)
   
   ! Assign FIRST atom types
   ! 
   ! Results for proteins agree with FIRST results.
   ! Note that carboxylate groups, phosphate groups etc. in
   ! non-protein residues (except nucleic acids) are assigned
   ! acceptor/donor types, but no negative/positive types.
   ! This is in agreement with FIRST. 
   !
   ! Holger Gohlke - 06.12.2002
   !
   !
   ! RNA is explicitly considered. 
   ! Note that both terminal oxygens of the phosphate groups
   ! are assigned negative.
   ! For further details about the RNA parametrization see
   ! Biophys. J. 2008, DOI:10.1529/biophysj.107.113415.
   !
   ! Note that not all modified RNA nucleosides are explicitly
   ! considered so far!
   ! Considered are pseudouridine (PSU), the methylation of 2'-O
   ! position of the ribose sugar (OMA, OMG, OMC, OMU),
   ! 6-hydro-1-methyladenosine (1MA), and 3-methyluridine (UR3). 
   !
   ! DNA is considered according to the RNA derived parametrization.
   !
   ! Simone Fulle - 25.02.2008
   !
   ! help functions
   !   logical: is_negative, is_positive, is_donor_acceptor
   !   others: ass_H, gethybrid
   !
   ! Assigned atom types:
   !   E: negative, C: positive
   !   A: acceptor, D: donor, B: donor & acceptor
   !   V: polar H,  N: not defined
   
   implicit none
   
   character(len=4), intent(in) :: igraph(*), lbres(*)
   character(len=4) :: tmpname
   
   integer, intent(in) :: natom, nres, nbond, ipres(*), ib(*), jb(*)
   integer :: i, j1, j2, k
   
   character(len=1), intent(out) :: ftype(*)
   
   logical :: found, is_negative, is_positive, is_donor_acceptor
   
   ! --- Init ftype ---
   do i=1,natom
      ftype(i) = ' '
   end do
   
   ! --- Consider non-H atoms for the moment ---
   do i=1,nres
      j1 = ipres(i)
      j2 = ipres(i+1) - 1
      do k=j1,j2
         write(tmpname,'(a4)') igraph(k)
         if (tmpname(1:1) .ne. 'H') then
            found = .false.
            found = is_negative(lbres(i), tmpname, ftype(k))
            if (.not. found) &
                found = is_positive(k, nbond, ib, jb, igraph, &
                                    lbres(i), tmpname, ftype(k))
            if (.not. found) &
                found = is_donor_acceptor(k, nbond, ib, jb, igraph, &
                                          lbres(i), tmpname, ftype(k))
            if (.not. found) ftype(k) = 'N'
         end if 
      end do
   end do
   
   ! --- Assign H atoms now ---
   do i=1,nres
      j1 = ipres(i)
      j2 = ipres(i+1) - 1
      do k=j1,j2
         write(tmpname,'(a4)') igraph(k)
         if (tmpname(1:1) .eq. 'H') call ass_H(k, nbond, ib, jb, ftype)
      end do
   end do
   
   return
   
end subroutine getftype

!=====================================================================

logical function is_negative(lbres, aname, ftype)
   
   ! Find deprotonated atoms
   !
   ! Holger Gohlke
   !   06.12.2002
   
   implicit none
   
   character(len=4), intent(in) :: lbres, aname
   
   character(len=1), intent(inout) :: ftype
   
   is_negative = .false.
   
   if ((lbres .eq. 'GLU' .and. aname(1:2) .eq. 'OE') .or. &
       (lbres .eq. 'ASP' .and. aname(1:2) .eq. 'OD') .or. &
       (lbres .eq. 'CYM' .and. aname(1:2) .eq. 'SG') .or. &
       aname(1:3) .eq. 'O1P' .or. &
       aname(1:3) .eq. 'O2P') then
      ftype = 'E'
      is_negative = .true.
   end if
   
   return
   
end function is_negative

!=====================================================================

logical function is_positive(nat,nbond,ib,jb,igraph,lbres,aname,ftype)
   
   ! Find protonated atoms
   !
   ! Holger Gohlke
   !   06.12.2002
   
   implicit none
   
   character(len=4), intent(in) :: igraph(*), lbres, aname
   
   integer, intent(in) :: nat, nbond, ib(*), jb(*)
   integer :: i, i1, i2, ihcnt, ibcnt
   
   character(len=1), intent(inout) :: ftype
   
   is_positive = .false.
   if (aname(1:1) .eq. 'N') then
      if ((lbres .eq. 'LYS' .and. aname(1:2) .eq. 'NZ') .or. &
          (lbres .eq. 'ARG' .and. &
           (aname(1:2).eq.'NE' .or. aname(1:2) .eq. 'NH')) .or. &
          (lbres .eq. 'HID' .and. aname(1:2) .eq. 'ND') .or. &
          (lbres .eq. 'HIE' .and. aname(1:2) .eq. 'NE') .or. &
          (lbres .eq. 'HIP' .and. &
           (aname(1:2) .eq. 'ND' .or. aname(1:2) .eq. 'NE'))) then
         ftype = 'C'
         is_positive = .true.
      else if (lbres .eq. '1MA' .and. aname(1:2) .eq. 'N1') then
         ftype = 'C'
         is_positive = .true.
      else
         ihcnt = 0
         ibcnt = 0
         do i=1,nbond
            i1 = ib(i)/3+1
            i2 = jb(i)/3+1
            if (i1 .eq. nat) then
               ibcnt = ibcnt + 1
               if (igraph(i2)(1:1) .eq. 'H') ihcnt = ihcnt + 1
            else if (i2 .eq. nat) then
               ibcnt = ibcnt + 1
               if (igraph(i1)(1:1) .eq. 'H') ihcnt = ihcnt + 1
            end if
         end do
         if (ibcnt .eq. 4 .and. ihcnt .gt. 0) then
            ftype = 'C'
            is_positive = .true.
         end if
      end if
   end if
   
   return
   
end function is_positive

!=====================================================================

logical function is_donor_acceptor(nat,nbond,ib,jb,igraph,lbres, &
                                   aname,ftype)
   
   ! Find donor/acceptor atoms
   !
   ! Holger Gohlke
   !   06.12.2002
   
   implicit none
   
   character(len=4), intent(in) :: igraph(*), lbres, aname
      
   integer, intent(in) :: nat, nbond, ib(*), jb(*)
   integer :: i, i1, i2
   
   character(len=1), intent(inout) :: ftype
   
   is_donor_acceptor = .false.
   ! -----------------------------------------------------------      
   if (aname(1:1) .eq. 'O') then
      if ((lbres .eq. 'SER' .and. aname(1:2) .eq. 'OG') .or. & 
          (lbres .eq. 'THR' .and. aname(1:2) .eq. 'OG') .or. &
          (lbres .eq. 'TYR' .and. aname(1:2) .eq. 'OH') .or. &
          (lbres(1:1) .eq. 'R' .and. aname .eq. 'O2''') .or. &
          (lbres .eq. 'WAT' .or. &
           lbres .eq. 'HOH' .or. lbres .eq. 'H2O' .or. &
           lbres .eq. 'OH2' .or. &
           lbres .eq. 'DOD' .or. lbres .eq. 'D2O' .or. &
           lbres .eq. 'OD2')) then
         ftype = 'B'
         is_donor_acceptor = .true.
      ! --- modified nucleosides of RNA ---
      else if ((lbres .eq. '1MA' .and. aname .eq. 'O2''') .or. &
               (lbres .eq. 'PSU' .and. aname .eq. 'O2''') .or. &
               (lbres .eq. 'UR3' .and. aname .eq. 'O2''')) then
                ftype = 'B'
                is_donor_acceptor = .true.    
      else if (aname .ne. 'O') then
         do i=1,nbond
            i1 = ib(i)/3+1
            i2 = jb(i)/3+1
            if ((i1 .eq. nat .and. igraph(i2)(1:1) .eq. 'H') .or. &
                (i2 .eq. nat .and. igraph(i1)(1:1) .eq. 'H')) then
               if (lbres .eq. 'ASH' .or. lbres .eq. 'GLH') then
                  ftype = 'D'
               else
                  ftype = 'B'
               end if
               is_donor_acceptor = .true.
               exit
            end if
         end do
      end if
      ! --- e.g. for methylated 2'-O position of the ribose sugar---
      if (.not. is_donor_acceptor) then
         ftype = 'A'
         is_donor_acceptor = .true.
      end if                  
   else if(aname(1:1) .eq. 'N') then
      if (lbres .eq. 'LYN' .and. aname(1:2) .eq. 'NZ') then
         ftype = 'B'
         is_donor_acceptor = .true.
      else if (aname .eq. 'N') then
         ftype = 'D'
         is_donor_acceptor = .true.
      else if (lbres(1:1) .eq. 'R' .or. lbres(1:1) .eq. 'D' .or. &
               lbres .eq. 'PSU' .or. lbres(1:2) .eq. 'OM' .or. &
               lbres .eq. '1MA' .or. lbres .eq. 'UR3') then   
         if (lbres(2:2) .eq. 'A' .or. lbres .eq. 'OMA' .or. &
            lbres .eq. '1MA') then
            !if (aname .eq. 'N1' .or. aname .eq. 'N3' .or. &
            !    aname .eq. 'N7') then
            !   ftype = 'A'
            if (aname.eq.'N6') then
               ftype = 'B'
               is_donor_acceptor = .true.    
            else if (aname.eq.'N9') then
               ftype = 'N'
               is_donor_acceptor = .true.
            end if
         else if (lbres(2:2) .eq. 'G' .or. lbres .eq. 'OMG') then
            !if (aname.eq.'N3' .or. aname.eq.'N7') then
            !   ftype = 'A'
            if (aname .eq. 'N1') then
               ftype = 'D'
               is_donor_acceptor = .true.
            else if (aname .eq. 'N2') then
               ftype = 'B'
               is_donor_acceptor = .true.
            else if (aname .eq. 'N9') then
               ftype = 'N'
               is_donor_acceptor = .true.
            end if
         else if(lbres(2:2).eq.'C' .or. lbres.eq.'OMC') then
            !if (aname .eq. 'N3') then
            !   ftype = 'A'
            if (aname .eq. 'N4') then
               ftype = 'B'
               is_donor_acceptor = .true.
            else if (aname .eq. 'N1') then
               ftype = 'N'
               is_donor_acceptor = .true.
            end if
         else if (lbres(2:2) .eq. 'U' .or. lbres .eq. 'OMU' .or. &
                  lbres(2:2) .eq. 'T') then
            if (aname .eq. 'N3') then
               ftype = 'D'
               is_donor_acceptor = .true.
            else if (aname .eq. 'N1') then
               ftype = 'N'
               is_donor_acceptor = .true.
            end if
         else if (lbres .eq. 'PSU') then
            if (aname .eq. 'N1' .or. aname .eq. 'N3') then
               ftype = 'D'
               is_donor_acceptor = .true.
            end if
         else if (lbres .eq. 'UR3') then
            if (aname .eq. 'N1' .or. aname .eq. 'N3') then
               ftype = 'N'
               is_donor_acceptor = .true.
            end if 
         end if
      else          
         do i=1,nbond
            i1 = ib(i)/3+1
            i2 = jb(i)/3+1
            if ((i1 .eq. nat .and. igraph(i2)(1:1) .eq. 'H') .or. &
                (i2 .eq. nat .and. igraph(i1)(1:1) .eq. 'H')) then
               ftype = 'D'
               is_donor_acceptor = .true.
               exit
            end if
         end do
      end if
      if (.not. is_donor_acceptor) then
         ftype = 'A'
         is_donor_acceptor = .true.
      end if
   else if(aname(1:1).eq.'S') then
      do i=1,nbond
         i1 = ib(i)/3+1
         i2 = jb(i)/3+1
         if ((i1 .eq. nat .and. igraph(i2)(1:1) .eq. 'H') .or. &
             (i2 .eq. nat .and. igraph(i1)(1:1) .eq. 'H')) then
            ftype = 'D'
            is_donor_acceptor = .true.
            exit
         end if
      end do
      if (.not. is_donor_acceptor) then
         ftype = 'A'
         is_donor_acceptor = .true.
      end if
   end if
   
   return
   
end function is_donor_acceptor

!=====================================================================

subroutine ass_H(nat, nbond, ib, jb, ftype)
   
   ! Assign H atoms
   !
   ! Holger Gohlke
   !   06.12.2002
   
   implicit none
   
   integer, intent(in) :: nat, nbond, ib(*), jb(*)
   integer :: i, i1, i2
   
   character(len=1), intent(inout) :: ftype(*)
   
   do i=1,nbond
      i1 = ib(i)/3+1
      i2 = jb(i)/3+1
      if (i1 .eq. nat) then
         if(ftype(i2) .ne. 'N') then
            ftype(i1) = 'V'
         else
            ftype(i1) = 'N'
         end if
         exit
      else if (i2 .eq. nat) then
         if (ftype(i1) .ne. 'N') then
            ftype(i2) = 'V'
         else
            ftype(i2) = 'N'
         end if
         exit
      end if
   end do
   
   return
   
end subroutine ass_H

!=====================================================================

subroutine gethybrid(nres,igraph,ipres,ib,jb,nbond,c,fhybrid)
   
   ! Determine hybridization of atoms
   !
   ! Holger Gohlke
   !   06.12.2001
   
   implicit none
   
   double precision, intent(in) :: c(*)
   double precision :: ang, &
                       ck1, ck2, ck3, c11, c12, c13, &
                       c21, c22, c23, cn1, cn2
   
   character(len=4), intent(in) :: igraph(*)
   character(len=4) :: tmpname
   
   integer, intent(in) :: nbond, nres, ipres(*), ib(*), jb(*)
   integer, intent(inout) :: fhybrid(*)
   integer :: i, i1, i2, j, j1, j2, k, l1, l2, m1, m2, &
              ibcnt, ibond(4)
   
   do i=1,nres
      j1 = ipres(i)
      j2 = ipres(i+1) - 1
      do k=j1,j2
         write(tmpname,'(a4)') igraph(k)
         if (tmpname(1:1) .eq. 'H') then
            fhybrid(k) = 0;
         else if (tmpname(1:1) .eq. 'C') then
            ibcnt = 0
            do j=1,nbond
               i1 = ib(j)/3+1
               i2 = jb(j)/3+1
               if (i1 .eq. k .or. i2 .eq. k) ibcnt = ibcnt + 1
            end do
            if (ibcnt .eq. 4) then
               fhybrid(k) = 3
            else if (ibcnt .eq. 3) then
               fhybrid(k) = 2
            else if (ibcnt .eq. 2 .or. ibcnt .eq. 1) then
               fhybrid(k) = 1
            end if
         else if (tmpname(1:1) .eq. 'O' .or. tmpname(1:1) .eq. 'S') &
         then
            ibcnt = 0
            do j=1,nbond
               i1 = ib(j)/3+1
               i2 = jb(j)/3+1
               if (i1 .eq. k .or. i2 .eq. k) ibcnt = ibcnt + 1
            end do
            if (ibcnt .eq. 2) then
               fhybrid(k) = 3
            else if (ibcnt .eq. 1) then
               fhybrid(k) = 2
            end if
         else if (tmpname(1:1) .eq. 'P') then
            fhybrid(k) = 3
         else if (tmpname(1:1) .eq. 'N') then
            ibcnt = 0
            do j=1,nbond
               i1 = ib(j)/3+1
               i2 = jb(j)/3+1
               if (i1 .eq. k .or. i2 .eq. k) then
                  ibcnt = ibcnt + 1
                  if(i1 .eq. k) ibond(ibcnt) = i2
                  if(i2 .eq. k) ibond(ibcnt) = i1
               end if
            end do
            if (ibcnt .eq. 4) then
               fhybrid(k) = 3
            else if (ibcnt .eq. 3) then
               ang = 0.0d0
               ck1 = c(3*k-2)
               ck2 = c(3*k-1)
               ck3 = c(3*k  )
               do l1=1,2
                  m1 = ibond(l1)
                  c11 = c(3*m1-2)
                  c12 = c(3*m1-1)
                  c13 = c(3*m1  )
                  cn1 = sqrt((c11-ck1)*(c11-ck1) + &
                             (c12-ck2)*(c12-ck2) + &
                             (c13-ck3)*(c13-ck3))
                  do l2 = l1+1,3
                     m2 = ibond(l2)
                     c21 = c(3*m2-2)
                     c22 = c(3*m2-1)
                     c23 = c(3*m2  )
                     cn2 = sqrt((c21-ck1)*(c21-ck1) + &
                                (c22-ck2)*(c22-ck2) + &
                                (c23-ck3)*(c23-ck3))
                     ang = ang + ((c21-ck1)*(c11-ck1) + &
                                  (c22-ck2)*(c12-ck2) + &
                                  (c23-ck3)*(c13-ck3))/ &
                                  (cn1 * cn2)
                  end do
               end do
               ang = ang / 3.0d0
               ! 114.5 degree = 1.9984 rad <=> cos(1.9984) = -0.41469
               if (ang .lt. -0.41469) then
                  fhybrid(k) = 2
               else
                  fhybrid(k) = 3
               end if
            else if (ibcnt .eq. 2) then
               fhybrid(k) = 2
            else if (ibcnt .eq. 1) then
               fhybrid(k) = 1
            end if
         else
            !write(6,*) 'WARNING: Assigning hybrid=0 to atom ',k, &
            !           ' (',tmpname,')'
            fhybrid(k) = 0
         endif
      end do
   end do
   
   return
   
end subroutine gethybrid
