subroutine findhbond(maxfhb,natom,igraph,ipres,lbres, &
                     ib,jb,nbond,c,ftype,fhybrid, &
                     nhb,fhbdon,fhbh,fhbacc,fhbene,fhbunit)
   
   ! Identifies H-bonds
   !
   ! Holger Gohlke - 06.12.2001
   !
   ! H-bonds occuring in proteins (1), in nucleic acids (2), and
   ! between (1) can be distinguished in the last column of the outfile.
   ! Useful to apply different Ehb in nucleic acid/protein complexes.
   !
   ! Simone Fulle - 26.02.2008
   
   implicit none
   
   double precision, intent(in) :: c(*)
   double precision, intent(inout) :: fhbene(*)
   double precision :: d(3), r(3), dtmp, donacc, hydacc, theta, &
                       ene, hbene, &
                       x1, y1, z1, x2, y2, z2, x3, y3, z3
   double precision, parameter :: t = 1.396263  ! Minimum don-H...acc angle
   double precision, parameter :: hx = 1.30     ! (heuristic) Maximum value for X-H bond length
   
   integer, intent(in) :: maxfhb, natom, nbond, ipres(*), ib(*), jb(*), &
                          fhybrid(*)
   integer, intent(inout) :: fhbdon(*), fhbacc(*), fhbh(*), fhbunit(*)
   integer, intent(out) :: nhb
   integer :: i, i1, i2, j, j1, j2, k, l, &
              residi, residj, unit, fhbtype
   
   character(len=4), intent(in) :: igraph(*), lbres(*)
   
   character(len=1), intent(in) :: ftype(*)
   character(len=1) :: ft1, ft2
   
   logical :: ishb
   
   ! --- Explicitly considered modified RNA nucleosides ---
   !     'PSU', 'OMA', 'OMG', 'OMC', 'OMU', '1MA', 'UR3'
   
   ! --- Parameters as described in suppl. to PNAS 2002, 99, 3540 ---
   !     Standard  Sulfur  Salt
   
   ! --- Max. donor-acceptor distance ---
   data d/3.6,      4.0,    4.6/
   
   ! --- Max. hydrogen-acceptor distance ---
   data r/2.6,      3.0,    3.6/
   
   nhb = 0
   residi = 1
   residj = 2
   
   do i=1,natom-1
      if (i .eq. ipres(residi+1)) then
         residi = residi+1
      end if
      if (i .eq. ipres(residi+1)-1) then
         residj = residi+1
      else
         residj = residi
      end if
      ! --- allows to use different Ehb in nucleic acid and protein units ---
      if (lbres(residi)(1:1) .eq. 'R' .or. &
          lbres(residi)(1:1) .eq. 'D' .or. &
          lbres(residi) .eq. 'PSU' .or. &
          lbres(residi)(1:2) .eq. 'OM' .or. &
          lbres(residi) .eq. '1MA' .or. &
          lbres(residi) .eq. 'UR3') then
         unit = 2
      else
         unit = 1
      end if
      
      ft1 = ftype(i)
      
      if (ft1 .ne. 'N' .and. ft1 .ne. 'V') then
         x1=c(3*i-2)
         y1=c(3*i-1)
         z1=c(3*i  )
         
         do j=i+1,natom
            if (lbres(residj)(1:1) .ne. 'R' .and. &
                lbres(residj)(1:1) .ne. 'D' .and. &
                lbres(residi) .ne. 'PSU' .and. &
                lbres(residi)(1:2) .ne. 'OM' .and. &
                lbres(residi) .ne. '1MA' .and. &
                lbres(residi) .ne. 'UR3') then
               unit = 1
               ! --- HBs between nucleic acids and proteins are marked ---
               ! --- according to HBs in proteins ---
            end if
            
            ft2 = ftype(j)
            if (ft2 .ne. 'N' .and. ft2 .ne. 'V' .and. &
                (((ft1 .eq. 'E' .or. ft1 .eq. 'A' .or. ft1 .eq. 'B') .and. &
                  (ft2 .eq. 'C' .or. ft2 .eq. 'D' .or. ft2 .eq. 'B')) .or. &
                 ((ft2 .eq. 'E' .or. ft2 .eq. 'A' .or. ft2 .eq. 'B') .and. &
                  (ft1 .eq. 'C' .or. ft1 .eq. 'D' .or. ft1 .eq. 'B')))) then
               ! --- Is HB or salt bridge ---
               x2=c(3*j-2)
               y2=c(3*j-1)
               z2=c(3*j  )
               if (igraph(i)(1:1) .eq. 'S' .or. &
                   igraph(j)(1:1) .eq. 'S') then
                  ! --- Sulfur involved ---
                  fhbtype = 2
               else if ((ft1 .eq. 'E' .and. ft2 .eq. 'C') .or. &
                        (ft1 .eq. 'C' .and. ft2 .eq. 'E')) then
                  ! --- Salt bridge ---
                  fhbtype = 3
               else
                  ! --- Standard HB ---
                  fhbtype = 1
               end if
               
               dtmp = r(fhbtype) + hx
               if (abs(x1-x2).le.dtmp .and. &
                   abs(y1-y2).le.dtmp .and. &
                   abs(z1-z2).le.dtmp) then
                  donacc = sqrt((x1-x2)*(x1-x2) + &
                                (y1-y2)*(y1-y2) + &
                                (z1-z2)*(z1-z2))
                  ishb = .false.
                  if (ft1 .eq. 'C' .or. ft1 .eq. 'D' .or. ft1 .eq. 'B') then
                     do k=1,nbond
                        i1 = ib(k)/3+1
                        i2 = jb(k)/3+1
                        if (i1 .eq. i .and. igraph(i2)(1:1) .eq. 'H' .or. &
                            i2 .eq. i .and. igraph(i1)(1:1) .eq. 'H') then
                           if (i1 .eq. i .and. igraph(i2)(1:1) .eq. 'H') then
                              l = i2
                              x3 = c(3*i2-2)
                              y3 = c(3*i2-1)
                              z3 = c(3*i2  )
                           else if (i2 .eq. i .and. &
                                    igraph(i1)(1:1) .eq. 'H') then
                              l = i1
                              x3 = c(3*i1-2)
                              y3 = c(3*i1-1)
                              z3 = c(3*i1  )
                           end if
                           call hbgeom(x1,y1,z1, x3,y3,z3, x2,y2,z2, &
                                       hydacc, theta)
                           if ((donacc .le. d(fhbtype) .or. &
                                hydacc .le. r(fhbtype)) &
                                .and. theta .ge. t) then
                              nhb = nhb + 1
                              if (nhb > maxfhb) then
                                 write(6,*) 'Too many HB found: ',nhb
                                 stop
                              end if
                              ene = hbene(c, nbond, ib, jb, i, l, j, &
                                          fhbtype, fhybrid, donacc, theta)
                              !write(6,*) 'Hbond: ',i,'-',l,'...',j,'  ', &
                              !           donacc, hydacc, theta, ene
                              fhbdon(nhb) = i
                              fhbh(nhb) = l
                              fhbacc(nhb) = j
                              fhbene(nhb) = ene
                              fhbunit(nhb) = unit
                              ishb = .true.
                           end if
                        end if
                     end do
                  end if
                  
                  if (.not.ishb .and. &
                      (ft2 .eq. 'C' .or. ft2 .eq. 'D' .or. ft2 .eq. 'B')) then
                     do k=1,nbond
                        j1 = ib(k)/3+1
                        j2 = jb(k)/3+1
                        if (j1 .eq. j .and. igraph(j2)(1:1) .eq. 'H' .or. &
                            j2 .eq. j .and. igraph(j1)(1:1) .eq. 'H') then
                           if (j1 .eq. j .and. igraph(j2)(1:1) .eq. 'H') then
                              l = j2
                              x3 = c(3*j2-2)
                              y3 = c(3*j2-1)
                              z3 = c(3*j2  )
                           else if (j2 .eq. j .and. &
                                    igraph(j1)(1:1) .eq. 'H') then
                              l = j1
                              x3 = c(3*j1-2)
                              y3 = c(3*j1-1)
                              z3 = c(3*j1  )
                           end if
                           call hbgeom(x2,y2,z2, x3,y3,z3, x1,y1,z1, &
                                       hydacc, theta)
                           if ((donacc .le. d(fhbtype) .or. &
                                hydacc .le. r(fhbtype)) &
                                .and. theta .ge. t) then
                              nhb = nhb + 1
                              if (nhb > maxfhb) then
                                 write(6,*) 'Too many HB found: ',nhb
                                 stop
                              end if
                              ene = hbene(c, nbond, ib, jb, j, l, i, &
                                          fhbtype, fhybrid, donacc, theta)
                              !write(6,*) 'Hbond: ',j,'-',l,'...',i,'  ', &
                              !            donacc, hydacc, theta, ene
                              fhbdon(nhb) = j
                              fhbh(nhb) = l
                              fhbacc(nhb) = i
                              fhbene(nhb) = ene
                              fhbunit(nhb) = unit
                              ishb = .true.
                           end if
                        end if
                     end do
                  end if
               end if          
            end if
         end do
      end if
   end do
   
   return
   
end subroutine findhbond

!=====================================================================

subroutine hbgeom(xd,yd,zd,xh,yh,zh,xa,ya,za,hydacc,theta)
   
   ! Calculates H...acc distance and don...h...acc angle
   !
   ! Holger Gohlke
   !   06.12.2001
   
   implicit none
   
   double precision, intent(in) :: xd,yd,zd, xh,yh,zh, xa,ya,za
   double precision, intent(inout) :: hydacc, theta
   
   double precision :: hyddon
   
   hydacc = sqrt((xa-xh)*(xa-xh) + (ya-yh)*(ya-yh) + (za-zh)*(za-zh))
   hyddon = sqrt((xd-xh)*(xd-xh) + (yd-yh)*(yd-yh) + (zd-zh)*(zd-zh))
   theta = acos(((xa-xh)*(xd-xh) + (ya-yh)*(yd-yh) + (za-zh)*(zd-zh)) / &
                          (hydacc * hyddon))
   
   return
   
end subroutine hbgeom

!=====================================================================
double precision function hbene(c,nbond,ib,jb,don,h,acc, &
                                fhbtype,fhybrid,donacc,theta)
   
   ! Calcs HB energy according to Mayo energy function 
   !   (mod. for FIRST, s. Suppl. to PNAS 2002, 99, 3540)
   !
   ! Holger Gohlke
   !   06.12.2001
   
   implicit none
   
   double precision, intent(in) :: c(*), donacc, theta
   double precision, parameter :: pi = 3.141592654
   double precision, parameter :: v0 = 8.0
   double precision, parameter :: r0 = 2.8
   double precision, parameter :: vs = 10.0
   double precision, parameter :: rs = 3.2
   double precision, parameter :: a = 0.375
   double precision :: phi, phimax, gamma, &
                       xd,yd,zd, xh,yh,zh, xa,ya,za, &
                       xv1,yv1,zv1, xv2,yv2,zv2, &
                       f, tmp, tmp3, tmp6, tmp10, tmp12
   
   integer, intent(in) :: nbond, ib(*), jb(*), fhybrid(*), &
                          fhbtype, don, h, acc
   integer :: i, i1, i2, j, k(2), kcnt, basecnt
   
   if (fhbtype .eq. 3) then
      ! --- Salt bridge ---
      tmp = rs / (donacc + 0.375)
      tmp3 = tmp*tmp*tmp
      tmp10 = tmp3*tmp3*tmp3*tmp
      tmp12 = tmp10*tmp*tmp
      hbene = vs * (5.0*tmp12 - 6.0*tmp10)
   else
      ! --- Normal HB (including sulfur) ---
      xd = c(3*don-2)
      yd = c(3*don-1)
      zd = c(3*don  )
      xh = c(3*h  -2)
      yh = c(3*h  -1)
      zh = c(3*h    )
      xa = c(3*acc-2)
      ya = c(3*acc-1)
      za = c(3*acc  )
      ! --- Calc distance dependent part for energy ---
      tmp = r0 / donacc
      tmp3 = tmp*tmp*tmp
      tmp10 = tmp3*tmp3*tmp3*tmp
      tmp12 = tmp10*tmp*tmp
      hbene = v0 * (5.0*tmp12 - 6.0*tmp10)
      ! --- Find base attached to acc, calc phi ---
      basecnt = 0
      phimax = 0.0
      do i=1,nbond
         i1 = ib(i)/3+1
         i2 = jb(i)/3+1
         if (i1 .eq. acc .or. i2 .eq. acc) then
            if (i1 .eq. acc) then
               j = i2
            else if (i2 .eq. acc) then
               j = i1
            end if
            basecnt = basecnt + 1
            call hbgeom(xh,yh,zh, xa,ya,za, &
                        c(3*j-2),c(3*j-1),c(3*j  ), &
                        tmp, phi)
            if (phi .gt. phimax) phimax = phi
         end if
      end do
      
      if (fhybrid(don) .eq. 2 .and. fhybrid(acc) .eq. 2) then
         ! --- Find atom attached to donor (not H) ---
         kcnt = 0
         do i=1,nbond
            i1 = ib(i)/3+1
            i2 = jb(i)/3+1
            if (i1 .eq. don .and. i2 .ne. h) then
               kcnt = kcnt + 1
               k(kcnt) = i2
            else if (i2 .eq. don .and. i1 .ne. h) then
               kcnt = kcnt + 1
               k(kcnt) = i1
            end if
            if (kcnt .eq. 2) exit
         end do
         
         i1 = k(1)
         i2 = k(2)
         call vecprod(c(3*i1-2)-xh, c(3*i1-1)-yh, c(3*i1  )-zh, &
                      c(3*i2-2)-xh, c(3*i2-1)-yh, c(3*i2  )-zh, &
                      xv1,yv1,zv1)
         ! --- Find atom attached to acc (except base) or 
         !     else attached to base (except acc) ---
         kcnt = 0
         do i=1,nbond
            i1 = ib(i)/3+1
            i2 = jb(i)/3+1
            if (i1 .eq. j .and. i2 .ne. acc) then
               kcnt = kcnt + 1
               k(kcnt) = i2
            else if (i2 .eq. j .and. i1 .ne. acc) then
               kcnt = kcnt + 1
               k(kcnt) = i1
            end if
            if (kcnt .eq. 2) exit
         end do
         
         i1 = k(1)
         i2 = k(2) 
         call vecprod(c(3*i1-2)-xa, c(3*i1-1)-ya, c(3*i1  )-za, &
                      c(3*i2-2)-xa, c(3*i2-1)-ya, c(3*i2  )-za, &
                      xv2,yv2,zv2)

         ! --- Calc gamma ---
         call hbgeom(xv1,yv1,zv1, 0.d0,0.d0,0.d0, xv2,yv2,zv2, &
                     tmp, gamma)
         if (gamma < 0.5*pi) gamma = pi - gamma
      end if
      
      ! --- Calc angular dependency ---
      tmp = cos(theta)
      f = tmp*tmp
      if (fhybrid(don) .eq. 2 .and. fhybrid(acc) .eq. 3) f = f*f

      tmp = pi - theta
      tmp3 = tmp*tmp*tmp
      tmp6 = tmp3*tmp3
      if (fhybrid(don) .eq. 2 .and. fhybrid(acc) .eq. 3) then
         f = f * exp(-2.0*tmp6)
      else
         f = f * exp(-tmp6)
      end if
      
      if (fhybrid(don) .eq. 3 .and. fhybrid(acc) .eq. 3) then
         tmp = cos(phimax - 1.911135)
      else if (fhybrid(don) .eq. 3 .and. fhybrid(acc) .eq. 2) then
         tmp = cos(phimax)
      else if (fhybrid(don) .eq. 2 .and. fhybrid(acc) .eq. 2) then
         tmp = cos(max(phimax, gamma))
      else
         tmp = 1.0
      end if
      
      f = f * tmp*tmp
      ! --- Calc total energy ---
      !write(6,*) don, h, acc, donacc, theta, phimax, gamma, &
      !           hbene, f
      hbene = hbene * f
   end if
   
   return
   
end function hbene

!=====================================================================
subroutine vecprod(x1,y1,z1, x2,y2,z2, xo,yo,zo)
   
   ! Calcs vector product
   !
   ! Holger Gohlke
   !   06.12.2001
   
   implicit none
   
   double precision, intent(in) :: x1,y1,z1, x2,y2,z2
   double precision, intent(out) :: xo,yo,zo
   
   xo = y1*z2 - z1*y2
   yo = z1*x2 - x1*z2
   zo = x1*y2 - y1*x2
   
   return
   
end subroutine vecprod

