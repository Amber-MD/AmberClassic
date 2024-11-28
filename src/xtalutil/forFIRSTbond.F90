subroutine corbondl(maxatom,natom,igraph,nres,ipres,lbres, &
                    nbond,ib,jb,c)
   
   ! Corrects bond list:
   !   - removes bonds between H's of water
   !   - adds bonds between metal ion and surroundings
   !
   ! Holger Gohlke
   !   06.12.2001
   
   implicit none
   
   integer, intent(in) :: maxatom, natom, nres, ipres(*)
   integer, intent(inout) :: nbond, ib(*), jb(*)
   double precision, intent(in) :: c(*)
   character(len=4), intent(in) :: igraph(*), lbres(*)
   
   logical :: found
   integer :: i, i1, i2, j, k, k1, k2, l
   double precision, parameter :: dist = 2.30
   double precision, parameter :: dist2 = 5.29
   double precision :: x1, y1, z1, x2, y2, z2, xd, yd, zd, tmp
   character(len=4) :: lbr
   
   ! --- Remove bonds between H's in water ---
   
   do i=1,nres
      lbr = lbres(i)
      if (lbr.eq.'WAT' .or. &
          lbr.eq.'HOH' .or. lbr.eq.'H2O' .or. &
          lbr.eq.'OH2' .or. &
          lbr.eq.'DOD' .or. lbr.eq.'D2O' .or. &
          lbr.eq.'OD2') then
         i1 = ipres(i)
         i2 = ipres(i+1) - 1
         l = 0
         do j=i1,i2
            if (igraph(j)(1:1) .eq. 'H' .or. igraph(j)(1:1) .eq. 'D') then
               do k=1,nbond
                  k1 = ib(k)/3+1
                  k2 = jb(k)/3+1
                  if((k1 .eq. j .and. (igraph(k2)(1:1) .eq. 'H' .or. &
                                       igraph(k2)(1:1) .eq. 'D')) &
                            .or. &
                     (k2 .eq. j .and. (igraph(k1)(1:1) .eq. 'H' .or. &
                                       igraph(k1)(1:1) .eq. 'D'))) then
                     l = k
                     goto 10
                  end if
               end do
            end if
         end do
      10 continue
         if (l .gt. 0) then
            !write(6,*) 'Removing bond between ', ib(l)/3+1, jb(l)/3+1
            ib(l) = ib(nbond)
            jb(l) = jb(nbond)
            ib(nbond) = 0
            jb(nbond) = 0
            nbond = nbond - 1
         end if
      end if
   end do
   
   ! --- Add bonds between metal ions and surroundings ---
   
   do i=1,natom
      found = .false.
      do j=1,nbond
         i1 = ib(j)/3+1
         i2 = jb(j)/3+1
         if (i .eq. i1 .or. i .eq. i2) then
            found = .true.
            exit
         end if
      end do
      
      if (.not. found) then
         ! --- Atom w/o bonds found => assuming metal ion ---
         !write(6,*) 'Found metal ion ', i
         x1 = c(3*i-2)
         y1 = c(3*i-1)
         z1 = c(3*i  )
         do k=1,natom
            if (i .ne. k) then
               x2 = c(3*k-2)
               y2 = c(3*k-1)
               z2 = c(3*k  )
               xd = x2 - x1
               yd = y2 - y1
               zd = z2 - z1
               if (abs(xd) .lt. dist .or. &
                   abs(yd) .lt. dist .or. &
                   abs(zd) .lt. dist) then
                  tmp = xd*xd + yd*yd + zd*zd
                  if (tmp .lt. dist2) then
                     !write(6,*) '  Found surrounding ', k
                     nbond = nbond + 1
                     if (nbond .gt. maxatom) then
                        write(6,*) 'Too many bonds: ', nbond
                        stop
                     end if
                     ib(nbond) = (i-1)*3
                     jb(nbond) = (k-1)*3
                  end if
               end if
            end if
         end do
      end if
   end do
   
   return
   
end subroutine corbondl

!=====================================================================

subroutine findtf(maxatom,nbond,ib,jb,fhybrid,igraph,ntf,itf,jtf)
   
   ! Finds bonds which have to be tethered
   !
   ! Holger Gohlke
   !   06.12.2001
   
   implicit none
   
   character(len=4), intent(in) :: igraph(*)
   integer, intent(in) :: maxatom, nbond, ib(*), jb(*), &
                          fhybrid(*)
   integer, intent(inout) :: ntf, itf(*), jtf(*)
   
   integer :: fh1, fh2, i, i1, i2, j, j1, j2, n1, n2
   integer, parameter :: maxrsize = 6
   
   logical :: tether, isinring, inring1, inring2, webstyle
   
   ! --- TRUE: CZ-NH{1,2} in ARG are tethered,
   !     but not C(O)-NH2 bonds in ASN, GLN ---
   data webstyle/.true./
   
   ntf = 0
   do i=1,nbond
      i1 = ib(i)/3 + 1
      i2 = jb(i)/3 + 1
      if (fhybrid(i1) .eq. 2 .and. fhybrid(i2) .eq. 2) then
         tether = .true.
         ! --- Don't tether sp2-sp2 bonds with one terminal atom ---
         n1 = 0
         n2 = 0
         do j=1,nbond
            if (i .ne. j) then
               j1 = ib(j)/3 + 1
               j2 = jb(j)/3 + 1
               fh1 = fhybrid(j1)
               fh2 = fhybrid(j2)
               if (webstyle) then
                  if ((fh2 .gt. 0 .or. igraph(j1)(1:2) .eq. 'NH') .and. &
                      i1 .eq. j1 .or. &
                      (fh1 .gt. 0 .or. igraph(j2)(1:2) .eq. 'NH') .and. &
                      i1 .eq. j2) n1 = n1 + 1
                  if ((fh2 .gt. 0 .or. igraph(j1)(1:2) .eq. 'NH') .and. &
                      i2 .eq. j1 .or. &
                      (fh1 .gt. 0 .or. igraph(j2)(1:2) .eq. 'NH') .and. &
                      i2 .eq. j2) n2 = n2 + 1
               else
                  if (fh2 .gt. 0 .and. i1 .eq. j1 .or. &
                      fh1 .gt. 0 .and. i1 .eq. j2) n1 = n1 + 1
                  if (fh2 .gt. 0 .and. i2 .eq. j1 .or. &
                      fh1 .gt. 0 .and. i2 .eq. j2) n2 = n2 + 1
               end if
               if (n1 .gt. 0 .and. n2 .gt. 0) exit
            end if
         end do
         
         if (n1 .eq. 0 .or. n2 .eq. 0) tether = .false.
         ! --- Don't tether bonds in rings ---
         if (tether) then
            inring1 = isinring(maxrsize, 0, nbond, ib, jb, i1, 0, i1)
            inring2 = isinring(maxrsize, 0, nbond, ib, jb, i2, 0, i2)
            if (inring1 .and. inring2) tether = .false.
         end if
         ! --- Tether bond ---
         if (tether) then
            ntf = ntf + 1
            if (ntf .gt. maxatom) then
               write(6,*) 'Too many tether bonds: ', ntf
               stop
            end if
            itf(ntf) = (i1-1)*3
            jtf(ntf) = (i2-1)*3
         end if
      end if
   end do
   
   return
   
end subroutine findtf

!=====================================================================

recursive function isinring(maxrsize,rsize,nbond,ib,jb,at,from,curr) &
                   result(inring)
   
   ! Determines if atom at is in a ring of size maxrsize or smaller
   !
   ! Holger Gohlke
   !   06.12.2001
   
   implicit none
   
   integer, intent(in) :: maxrsize, rsize, nbond, ib(*), jb(*), &
                          at, from, curr
   
   integer :: i, i1, i2
   logical :: inring
   
   inring = .false.
   if (rsize .gt. 0 .and. curr .eq. at) then
      !write(6,*) 'Atom in ring ', at
      inring = .true.
   else if ((rsize+1) .le. maxrsize) then
      do i=1,nbond
         i1 = ib(i)/3+1
         i2 = jb(i)/3+1
         if (i1 .eq. curr .and. i2 .ne. from) then
            inring = inring .or. &
                     isinring(maxrsize,rsize+1, & 
                              nbond,ib,jb,at,curr,i2)
         else if (i2 .eq. curr .and. i1 .ne. from) then
            inring = inring .or. &
                     isinring(maxrsize,rsize+1, &
                              nbond,ib,jb,at,curr,i1)
         end if
         if (inring) exit
      end do
   end if
   
   return
   
end function isinring
