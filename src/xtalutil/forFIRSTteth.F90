subroutine addtether(maxfhb,maxatom,maxres,natom,nres, &
                     igraph,ipres,lbres,ib,jb,nbond,c, &
                     ftype,fhybrid, &
                     nhb,fhbdon,fhbh,fhbacc,fhbene)
   
   ! Adds tether atoms to simulate hydrophobic interactions.
   !
   ! Holger Gohlke - 06.12.2001
   !
   !
   ! Base stacking interactions between bases are limited to 1. 
   !   This avoids over-rigidification of nucleic acid structures.
   ! The threshold for hydrophobic interactions in proteins is set 
   !   to 0.25, and in nucleic acids to 0.15. 
   ! Note that the parametrization for nucleic acids has been
   !   tested so far only on RNA structures!
   !
   ! Simone Fulle - 26.02.2008
   
   implicit none
   
   double precision, intent(inout) :: c(*), fhbene(*)
   double precision :: vecfac, cutoff, &
                       thresh, threshprot, threshrna, &
                       vdw, vdw1, vdw2, sum, sum2, dist2, ene, &
                       x1, y1, z1, x2, y2, z2, xd, yd, zd, array, &
                       xcoord, ycoord, zcoord
   
   character(len=4), intent(inout) :: igraph(*), lbres(*)
   character(len=4) :: aname
   
   ! --- Either MAXTETH or MAXDOF must be >= 0 ---
   ! --- MAXTETH >= 0: fixed number of tethering atoms applied ---
   ! --- MAXDOF  >= 0: number of degrees of freedom REMOVED per tether from the system ---
   integer, parameter :: maxneigh = 14
   integer, parameter :: maxdepth = 3
   integer, parameter :: maxteth = 3
   integer, parameter :: maxdof = -1
   integer, intent(in) :: maxfhb, maxatom, maxres
   integer, intent(inout) :: natom, nres, ipres(*), nbond, nhb, &
                             ib(*), jb(*), fhybrid(*), fhbh(*), &
                             fhbdon(*), fhbacc(*)
   integer :: d, i, i2, ith, j, k, n, t, nat, locali, localj, &
              neigh(maxneigh, maxatom), neighnum(maxatom), &
              neighnumfst(maxatom), tethnum, residjth, &
              residi, residj, tethercut, count, residith, nrs
   
   character(len=1), intent(inout) :: ftype(*)
   
   logical :: go
   
   dimension vdw(2)
   dimension array((nres*2),5)
   
   ! --- Bondi radii for C, S ---
   data vdw/1.7, 1.8/
   
   ! --- Threshold: dist < vdw1 + vdw2 + thresh <=> hydrophobic contact ---
   data threshprot/0.25/
   data threshrna/0.15/
   
   ! --- Default energy value for tether ---
   data ene/-9.999999/
   
   ! --- number of tethers considered for base stacking = {1,2} ---
   data tethercut/1/
   
   do i=1,maxatom                                                  
      neighnum(i) = 0
   end do
   
   ! --- Check consistency of maxteth and maxdof
   if ((maxteth .ge. 0 .and. maxdof .ge. 0) .or. &
       (maxteth .lt. 0 .and. maxdof .lt. 0)) then
      write(6,*) 'Either MAXTETH or MAXDOF, but not both, must be >= 0', &
                  maxteth, maxdof
      stop
   end if
   
   ! --- Pre-calculate some values for fixed maxteth
   if (maxteth .ge. 0) then
      tethnum = maxteth
      ! --- Factor to calc coords of tether atoms ---
      vecfac = 1.0 / dble(tethnum + 1)
   end if
   
   ! --- Build list of bonded neighb up to maxdepth ---
   do i=1,natom
      if (igraph(i)(1:1) .eq. 'C' .or. igraph(i)(1:1) .eq. 'S') then
         neighnum(i) = 0
         neighnumfst(i) = 0
         call buildneigh(maxneigh, maxatom, maxdepth, &
                          i, 0, i, nbond, ib, jb, igraph, &
                          0, neigh, neighnum, neighnumfst)
      end if
   end do
   
   ! --- Find tether atoms. Tethers can not be between atoms 
   !     of bond distance <= maxdepth ---
   
   nat = natom
   nrs = nres
   residi = 1
   residj = 2
   residjth = 0
   residith = nres+1
   cutoff = 10000.0
   
   do i=1,(nrs*2)
      array(i,5) = cutoff
   enddo
   
   ! --- Loop over all atoms i ---
   do i=1,natom-1
      
      if (i .eq. ipres(residi+1)) then
         residi = residi+1
         residjth = 0
      end if
      
      if (i .eq. (ipres(residi+1)-1)) then
         residj = residi+1
      else 
         residj = residi
      end if
      
      ! --- set threshold for hc ---
      if (lbres(residi)(1:1) .eq. 'R' .or. &
          lbres(residi)(1:1) .eq. 'D' .or. &
          lbres(residi) .eq. 'PSU' .or. &
          lbres(residi)(1:2) .eq. 'OM' .or. &
          lbres(residi) .eq. '1MA' .or. &
          lbres(residi) .eq. 'UR3') then
         thresh = threshrna
      else
         thresh = threshprot 
      end if
      
      aname = igraph(i)
      if (aname(1:1) .eq. 'C' .or. aname(1:1) .eq. 'S') then
         
         if (aname(1:1) .eq. 'C') then
            vdw1 = vdw(1)
         else
            vdw1 = vdw(2)
         end if
         
         x1 = c(3*i-2)
         y1 = c(3*i-1)
         z1 = c(3*i  )
         
         ! --- Loop over all atoms j > i ---
         do j=i+1,natom
            
            if (j .eq. ipres(residj+1)) residj = residj+1
            
            ! ---  set threshold for hc between nucleic acid and protein units ---
            if (lbres(residj)(1:1) .ne. 'R' .and. &
                lbres(residj)(1:1) .ne. 'D' .and. &
                lbres(residi) .ne. 'PSU' .and. &
                lbres(residi)(1:2) .ne. 'OM' .and. &
                lbres(residi) .ne. '1MA' .and. &
                lbres(residi) .ne. 'UR3') thresh = threshprot
            
            ! --- for all carbon atoms in nucleic acid bases ---
            if ((igraph(i) .eq. 'C2' .or. igraph(i) .eq. 'C4' .or. &
                 igraph(i) .eq. 'C5' .or. igraph(i) .eq. 'C6' .or. &
                 igraph(i) .eq. 'C8') .and. &
                (igraph(j) .eq. 'C2' .or. igraph(j) .eq. 'C4' .or. &
                 igraph(j) .eq. 'C5' .or. igraph(j) .eq. 'C6' .or. &
                 igraph(j) .eq. 'C8')) then
               
               ! --- which do not belong to the same nucleotide---
               if (residi .ne. residj) then
                  aname = igraph(j)
                  go = .true.
                  n = neighnum(j)
                  ! --- Check that no bond distance <= MAXDEPTH ---
                  do k=1,n
                     if (neigh(k,j) .eq. i) then
                        go = .false.
                        exit
                     end if
                  end do
                  
                  if (go) then
                     vdw2 = vdw(1)
                     ! --- Check distance ---
                     sum = vdw1 + vdw2 + threshrna
                     sum2 = sum * sum
                     x2 = c(3*j-2)
                     y2 = c(3*j-1)
                     z2 = c(3*j)
                     xd = x2 - x1
                     yd = y2 - y1
                     zd = z2 - z1
                     
                     if (abs(xd).lt.sum .or. &
                         abs(yd).lt.sum .or. &
                         abs(zd).lt.sum) then
                        dist2 = xd*xd + yd*yd + zd*zd
                        if (dist2 .lt. sum2) then
                           
                           ! --- candidate for hc ---
                           do ith=1, residjth
                              if (residj .eq. array(ith*2-1,2)) then
                                 residith = ith
                                 exit
                              end if
                              ! --- new residue
                              residith = residjth+1 
                           end do
                           
                           ! --- save per base stacking 2 hc with smallest distance ---
                           if (residith .le. residjth) then
                              ! --- old residue j ---
                              if (dist2 .lt. array(residith*2-1,1)) then
                                 array(residith*2,1) = array(residith*2-1,1)
                                 array(residith*2,2) = array(residith*2-1,2)
                                 array(residith*2,3) = array(residith*2-1,3)
                                 array(residith*2,4) = array(residith*2-1,4)
                                 array(residith*2-1,1) = dist2
                                 array(residith*2-1,2) = residj
                                 array(residith*2-1,3) = i
                                 array(residith*2-1,4) = j
                              else if (dist2 .lt. array(residith*2,1)) then
                                 array(residith*2,1) = dist2
                                 array(residith*2,2) = residj
                                 array(residith*2,3) = i
                                 array(residith*2,4) = j
                              end if
                           else
                              ! --- new residue j ---
                              residjth=residjth+1
                              array(residjth*2-1,1)=dist2
                              array(residjth*2-1,2)=residj
                              array(residjth*2-1,3)=i
                              array(residjth*2-1,4)=j
                           end if  
                        end if
                     end if  
                  end if
               end if  
            ! --- end different residues ---
            
            ! --- within the sugar rings ---
            
            ! Ben Roberts says: I discovered code here that would identify atoms
            ! within the sugar rings. But it was an empty if. It appears
            ! that once upon a time, code existed to, and I quote,
            !    forbid tether between C1' and C5' of the same residue    
            !    required if 2 pseudoatoms are inserted into one sugar ring
            ! After it, came a catch-all "else". To keep the identification of
            ! sugar rings so tethers were ignored, but remove the empty if,
            ! I reversed the conditional.
            ! Original conditional:
            ! if (((igraph(i) .eq. 'C1''' .and. igraph(j) .eq. 'C5''') .or. &
            !      (igraph(i) .eq. 'C5''' .and. igraph(j) .eq. 'C1''')) .and. &
            !      residi .eq. residj) then
            !    apparently, do nothing
            ! else...insert tether for all other cases (see below)
            
            ! --- insert tether for all other cases ---
            else if (((igraph(i) .ne. 'C1''' .or. igraph(j) .ne. 'C5''') .and. &
                      (igraph(i) .ne. 'C5''' .or. igraph(j) .ne. 'C1''')) .or. &
                      residi .ne. residj) then
                     
               aname = igraph(j)
               go = .true.
               n = neighnum(j)
               
               ! --- Check that no bond distance <= MAXDEPTH ---
               do k=1,n
                  if (neigh(k,j).eq.i) then
                     go = .false.
                     exit
                  end if
               end do
               
               if (go .and. &
                   (aname(1:1) .eq. 'C' .or. aname(1:1) .eq. 'S')) then
                  if (aname(1:1) .eq. 'C') then
                     vdw2 = vdw(1)
                  else
                     vdw2 = vdw(2)
                  end if
                  
                  ! --- Check distance ---
                  sum = vdw1 + vdw2 + thresh
                  sum2= sum * sum
                  x2=c(3*j-2)
                  y2=c(3*j-1)
                  z2=c(3*j  )
                  xd = x2 - x1
                  yd = y2 - y1
                  zd = z2 - z1
                  if (abs(xd) .lt. sum .or. &
                      abs(yd) .lt. sum .or. &
                      abs(zd) .lt. sum) then
                     dist2 = xd*xd + yd*yd + zd*zd
                     if (dist2 .lt. sum2) then
                        
                        ! --- Add tether atoms ---
                        
                        ! --- Update residue information ---
                        nres = nres + 1
                        
                        if (nres .gt. maxres) then
                           write(6,*) 'Too many residues: ', nres
                           stop
                        end if
                        
                        lbres(nres) = 'BMH'
                        ipres(nres) = nat + 1
                        
                        ! --- Calc some values if MAXDOF >= 0 --- 
                        if (maxdof .ge. 0) then
                           tethnum = -maxdof + 1 + neighnumfst(i) + &
                                     neighnumfst(j)
                           vecfac = 1.0 / dble(tethnum + 1)
                        end if
                        
                        ! --- Update atom and bond information ---                  
                        xd = xd * vecfac
                        yd = yd * vecfac
                        zd = zd * vecfac
                        
                        do k=1,tethnum
                           
                           nat = nat + 1
                           
                           if (nat .gt. maxatom) then
                              write(6,*) 'Too many atoms: ', nat
                              stop
                           end if
                           
                           ! BPR: Deal with round-off error differences
                           ! between compilers - use a definite rounding
                           ! algorithm.
                           xcoord = x1 + (k*xd)
                           c(3*nat-2) = dble(idnint(xcoord*1000))/1000
                           ycoord = y1 + (k*yd)
                           c(3*nat-1) = dble(idnint(ycoord*1000))/1000
                           zcoord = z1 + (k*zd)
                           c(3*nat) = dble(idnint(zcoord*1000))/1000
                           
                           igraph(nat) = 'X'
                           ftype(nat) = 'N'
                           fhybrid(nat) = 0
                           nbond = nbond + 1
                           
                           if (nbond .gt. maxatom) then
                              write(6,*) 'Too many bonds: ', nbond
                              stop
                           end if
                           
                           if (k .eq. 1) then
                              ib(nbond) = (i-1)*3
                           else
                              ib(nbond) = (nat-2)*3
                           end if
                           
                           jb(nbond) = (nat-1)*3
                        end do
                        
                        ! --- Update H-bond information ---
                        nhb = nhb+1
                        if (nhb .gt. maxfhb) then
                           write(6,*) 'Too many H-bonds: ', nhb
                           stop
                        end if
                        fhbdon(nhb) = nat - 1
                        fhbh(nhb) = nat
                        fhbacc(nhb) = j
                        fhbene(nhb) = ene
                     end if
                  end if
               end if 
            end if  
         end do
         
         ! --- add hc between bases --
         if (((lbres(residi)(1:2).eq.'RC' .or. &
               lbres(residi)(1:2).eq.'RU' .or. &
               lbres(residi)(1:2).eq.'DC' .or. &
               lbres(residi)(1:2).eq.'DT' .or. &
               lbres(residi).eq.'PSU' .or. &
               lbres(residi).eq.'OMC' .or. &
               lbres(residi).eq.'OMU' .or. &
               lbres(residi).eq.'UR3') .and. &
               igraph(i).eq.'C2') .or. &
             ((lbres(residi)(1:2).eq.'RA' .or. &
               lbres(residi)(1:2).eq.'RG' .or. &
               lbres(residi)(1:2).eq.'DA' .or. &
               lbres(residi)(1:2).eq.'DG' .or. &
               lbres(residi).eq.'OMA' .or. &
               lbres(residi).eq.'OMG' .or. &
               lbres(residi).eq.'1MA') .and. &
              igraph(i).eq.'C4')) then
            
            do t=1,tethercut
               if (t .eq. 1) then
                  d = 1
               else
                  d = 0  
               end if
               
               do count=1,residjth
                  locali = array(count*2-d,3)
                  localj = array(count*2-d,4)
                  vdw1 = vdw(1)
                  vdw2 = vdw(1)
                  sum = vdw1 + vdw2 + threshrna
                  sum2 = sum * sum
                  x1 = c(3*locali-2)
                  y1 = c(3*locali-1)
                  z1 = c(3*locali  )
                  x2 = c(3*localj-2)
                  y2 = c(3*localj-1)
                  z2 = c(3*localj)
                  xd = x2 - x1
                  yd = y2 - y1
                  zd = z2 - z1
                  dist2 = array(count*2-d,1)
                  
                  ! --- Add tether atoms ---
                  if (dist2 .lt. cutoff .and. dist2 .lt. sum2) then
                     ! --- Update residue information ---
                     nres = nres + 1
                     if (nres .gt. maxres) then
                        write(6,*) 'Too many residues: ', nres
                        stop
                     end if
                     lbres(nres) = 'BMH'
                     ipres(nres) = nat + 1
                     
                     ! --- Calc some values if MAXDOF >= 0 ---
                     if (maxdof .ge. 0) then
                        tethnum = -maxdof + 1 + neighnumfst(locali) + &
                                  neighnumfst(localj)
                        vecfac = 1.0 / dble(tethnum + 1)
                     end if
                     
                     ! --- Update atom and bond information ---                  
                     xd = xd * vecfac
                     yd = yd * vecfac
                     zd = zd * vecfac
                     
                     do k=1,tethnum
                        nat = nat + 1
                        if (nat .gt. maxatom) then
                           write(6,*) 'Too many atoms: ', nat
                           stop
                        end if
                        
                        ! BPR: Deal with round-off error differences
                        ! between compilers - use a definite rounding
                        ! algorithm.
                        xcoord = x1 + (k*xd)
                        c(3*nat-2) = dble(idnint(xcoord*1000))/1000
                        ycoord = y1 + (k*yd)
                        c(3*nat-1) = dble(idnint(ycoord*1000))/1000
                        zcoord = z1 + (k*zd)
                        c(3*nat) = dble(idnint(zcoord*1000))/1000
                        
                        igraph(nat) = 'X'
                        ftype(nat) = 'N'
                        fhybrid(nat) = 0
                        nbond = nbond + 1
                        
                        if (nbond .gt. maxatom) then
                           write(6,*) 'Too many bonds: ', nbond
                           stop
                        end if
                        
                        if (k .eq. 1) then
                           ib(nbond) = (locali-1)*3
                        else
                           ib(nbond) = (nat-2)*3
                        end if
                        
                        jb(nbond) = (nat-1)*3
                     end do
                     
                     ! --- Update H-bond information ---
                     nhb = nhb + 1
                     if (nhb .gt. maxfhb) then
                        write(6,*) 'Too many H-bonds: ', nhb
                        stop
                     end if
                     
                     fhbdon(nhb) = nat - 1
                     fhbh(nhb) = nat
                     fhbacc(nhb) = localj
                     fhbene(nhb) = ene
                  end if  
               end do
            end do
            
            do i2=1,(nrs*2)
               array(i2,1) = cutoff
               array(i2,2) = 0
            end do
            
         end if         
      end if  
   end do
   
   natom = nat
   if((nres+1) .gt. maxres) then
      write(6,*) 'Too many residues: ', nres+1
      stop
   end if
   
   ipres(nres+1) = natom+1
   
   return
   
end subroutine addtether

!=====================================================================

recursive subroutine buildneigh(maxneigh, maxatom, maxdepth, &
                                at, from, curr, nbond, ib, jb, igraph, &
                                ndepth, neigh, neighnum, neighnumfst)
   
   ! Recursively find neighbors of iat
   !
   ! Holger Gohlke
   !   06.12.2001
   
   implicit none
   
   character(len=4), intent(in) :: igraph(*)
   
   integer, intent(in) :: maxneigh, maxatom, maxdepth, &
                          at, from, curr, nbond, ib(*), jb(*), ndepth
   integer, intent(inout) :: neigh(maxneigh, maxatom), &
                             neighnum(maxatom), neighnumfst(maxatom)
   integer :: to, currdepth, nn, i, i1, i2
   
   logical :: go
   
   currdepth = ndepth + 1
   if (currdepth .le. maxdepth) then
      do i=1,nbond
         i1 = ib(i)/3+1
         i2 = jb(i)/3+1
         go = .false.
         if (i1 .eq. curr .and. i2 .ne. from) then
            to = i2
            go = .true.
         else if (i2 .eq. curr .and. i1 .ne. from) then
            to = i1
            go = .true.
         end if
         
         if (go) then
            ! --- Store number of neighbors in first shell ---
            if (currdepth .eq. 1) then
               neighnumfst(at) = neighnumfst(at) + 1
            end if
            ! --- Store neighbor information ---
            if (igraph(to)(1:1) .eq. 'C' .or. &
                igraph(to)(1:1) .eq. 'S') then
               nn = neighnum(at)
               nn = nn + 1
               if (nn .gt. maxneigh) then
                  write(6,*) 'Too many neighbors: ', nn
                  stop
               end if
               neighnum(at) = nn
               neigh(nn, at) = to
               !write(6,*) ndepth+1, at, ' (', nn,') ', to
            end if
            call buildneigh(maxneigh, maxatom, maxdepth, &
                            at, curr, to, nbond, ib, jb, igraph, &
                            currdepth, neigh, neighnum, neighnumfst)
         end if
      end do
   end if
   
   return
   
end subroutine buildneigh

