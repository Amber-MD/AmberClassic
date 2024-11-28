subroutine writebond(nf,nbond,ib,jb)
   
   ! Writes bond list for FIRST input file
   !
   ! Holger Gohlke
   !   06.12.2001
   
   implicit none
   
   integer, intent(in) :: nf, nbond, ib(*), jb(*)
   integer :: i, i1, i2
   
   write(nf,100)
   write(nf,110)
   
   do i=1,nbond
      i1 = ib(i)/3 + 1
      i2 = jb(i)/3 + 1
      if (i1 .lt. i2) then
         write(nf,120) i1, i2
      else
         write(nf,120) i2, i1
      end if
   end do
   
   100 format('REMARK:L:----------------------------------------------', &
              '----------------------------------') 
   110 format('REMARK:cf','    so','    sf','     (atom-label pairs)') 
   120 format('REMARK:CF',2i6)
   
   return
   
end subroutine writebond

!=====================================================================
subroutine writetf(nf,ntf,itf,jtf)
   
   ! Writes tf list for FIRST input file
   !
   ! Holger Gohlke
   !   06.12.2001
   
   implicit none
   
   integer, intent(in) :: nf, ntf, itf(*), jtf(*)
   integer :: i, i1, i2
   
   write(nf,100)
   write(nf,110)
   
   do i=1,ntf
      i1 = itf(i)/3 + 1
      i2 = jtf(i)/3 + 1
      if (i1 .lt. i2) then
         write(nf,120) i1, i2
      else
         write(nf,120) i2, i1
      end if
   end do
   
   100 format('REMARK:L:----------------------------------------------', &
              '----------------------------------') 
   110 format('REMARK:tf','    so','    sf','     (atom-label pairs)') 
   120 format('REMARK:TF',2i6)
   
   return
   
end subroutine writetf

!=====================================================================
subroutine writehb(nf, nhb, ftype, fhybrid, fhbdon, fhbh, &
                   fhbacc, fhbene, igraph, hbene, fhbunit)
   
   ! Writes HB for FIRST input file
   !
   ! Holger Gohlke
   !   06.12.2001
   
   implicit none
   
   double precision, intent(in) :: fhbene(*)
   double precision :: ene, eneold, tmp, hbene, this_hbene
   
   character(len=4), intent(in) :: igraph(*)
   
   integer, intent(in) :: nf, nhb, fhybrid(*), fhbdon(*), fhbh(*), &
                          fhbacc(*), fhbunit(*)
   integer :: i, iold, icurr, iacc, idon, j, cnt
   
   character(len=1), intent(in) :: ftype(*)
   
   write(nf, 100)
   write(nf, 110)
   
   ene = fhbene(1)
   icurr = 1
   do i=2,nhb
      if (fhbene(i) .gt. ene) then
         ene = fhbene(i)
         icurr = i
      end if
   end do
   
   cnt = 0
   
   do i=1,nhb
      iacc = fhbacc(icurr)
      idon = fhbdon(icurr)
      if (fhbene(icurr) .le. hbene) then
         
         ! For printing purposes, assume that any number between
         ! -0.000005 and 0.000005 is zero. Thus, we stop complaints
         ! about differences between "0.00000" and "-0.00000".
         this_hbene = fhbene(icurr)
         if (this_hbene .lt. 0.000005 .and. this_hbene .gt. -0.000005) &
            this_hbene = 0.00000
         
         cnt = cnt + 1
         if (ftype(iacc) .eq. 'E' .and. ftype(idon) .eq. 'C' .or. &
             ftype(iacc) .eq. 'C' .and. ftype(idon) .eq. 'E') then
            write(nf, 130) cnt, this_hbene, fhbdon(icurr), &
               fhbh(icurr), fhbacc(icurr)
         else if(igraph(iacc)(1:1) .eq. 'X' .or. &
                 igraph(idon)(1:1) .eq. 'X') then
            write(nf, 140) cnt, this_hbene, fhbdon(icurr), &
                           fhbh(icurr), fhbacc(icurr)
         else
            write(nf, 120) cnt, this_hbene, fhbdon(icurr), &
                           fhbh(icurr), fhbacc(icurr), fhybrid(idon), &
                           fhybrid(iacc), fhbunit(icurr)
         end if
      end if
      
      eneold = ene
      iold = icurr
      ene = -1.0d10
      do j=1,nhb
         tmp = fhbene(j)
         if (tmp .eq. eneold .and. j .gt. iold) then
            icurr = j
            ene = tmp
            goto 10
         else if (tmp .lt. eneold .and. tmp .gt. ene) then
            icurr = j
            ene = tmp
         end if          
      end do
   10 continue
   end do   
   
   100 format('REMARK:L:----------------------------------------------', &
              '----------------------------------') 
   110 format('REMARK:hb','   ID','   E(Kcal/Mol)','  Donor', &
              '  Hydrogen',' Acceptor','  Description') 
   120 format('REMARK:HB',i5,1x,f12.5,3(1x,i7),4x, &
              'HB',1x,'Dsp',i1,1x,'Asp',i1,2x,'hbtype:',1x,i1)
   130 format('REMARK:HB',i5,1x,f12.5,3(1x,i7),4x, &
              'SB',1x,'no energy')
   140 format('REMARK:HB',i5,1x,f12.5,3(1x,i7),4x, &
              'PH',1x,'hydr',1x,'phob')
   
   return
   
end subroutine writehb

