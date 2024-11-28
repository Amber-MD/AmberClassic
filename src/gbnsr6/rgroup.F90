#include "copyright.h"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rgroup here]
subroutine rgroup(natom,natc,nres,ngrp,ipres,lbres,igraph,isymbl, &
      itree,igroup,jgroup,index,irespw,npdec, &
      weit,xc,konst,dotgtmd,belly,idecomp,nfu,writeout)
   
   implicit none
   
   integer:: i, i1, i2, idecomp, iend, ifld, igroup, igrp, iii, &
        index, ipres, irespw, isear, isrch, istart, itime, ivar, izero, j, &
        jfld, jgroup, jgrp, k, l, lfind, lsign, natc, natmg, natom, nf, &
        nfu, ngrp, npdec, nres
   _REAL_ :: fvar, weit, wt, xc

   !     Mods for Rev A by gls:
   !     - cpp selectable precision
   !     - wildcard functionality fixed. (see findp)
   !     - format change to print less garbage
   !     - unreferenced statement labels removed
   !     - changed title write from 12 to 19 A4. (not 20 to avoid wraps)
   !     - changed to read from unit nf -- dac 12/90
   
   logical konst,misc,belly,flres,frres,dotgtmd
   character(len=4) lbres,igraph,isymbl,itree

   ! Ross Walker modifications for dipole printing - flag to control writing
   ! of group data to output file
   ! if (writeout) then  - wrapper placed around all write(6,*) statements
   logical writeout

   !     ----- READ IN GROUPS EACH ATOM IS IN -----
   
   !           WILL KEEP READING IN GROUPS OF CARDS UNTIL A BLANK CARD IS
   !           READ
   !           ALL ATOMS IN THIS CARD GROUP WILL BE PUT IN THE SAME GROUP
   !           THE ONLY EXCEPTION IS A RES CARD WHICH STARTS WITH A -
   !           RES NO. THE RESIDUES IN THIS GROUP WILL ALL BE PUT IN INDIV.
   !           GROUPS ANY OTHER GROUPS IN THIS SECTION MUST ALSO START
   !           WITH A  I.E. THEY MUST ALSO BE INDIV. GROUPS
   !           ANY TIME A "FIND" CARD IS READ, ONLY THOSE ATOMS THAT MEET
   !           THE SPECIFICATIONS ON AT LEAST 1 OF THE FIND CARDS
   !           WILL BE INCLUDED IN THE GROUP
   !           THIS OPTION MAY BE ENDED BY READING IN "NFIND",TERMINATING
   !           THE INPUT FOR THE GROUPS, OR BY READING IN ANOTHER "FIND"
   
   integer, parameter :: max_finds=100
   character(len=4), dimension(1:max_finds) :: jgraph,jresnm,jsymbl,jtree
   common/propf/jgraph,jresnm,jsymbl,jtree,isrch
   character(len=4) title1(20)
   
   dimension ipres(*),igroup(*),jgroup(*),index(*),irespw(*), &
         weit(*),igrp(8),jgrp(8)
   dimension lbres(*),igraph(*),isymbl(*),itree(*),xc(*)
   character(len=4) ihol,katn,ifind,nfind,iiend,ksear,ires,ilres,irres,ktypg
   dimension ifld(20),jfld(20),ihol(20),ivar(20),fvar(20)
   character(len=4) igrap(8),isymb(8),iwild,jwild
   
   data ifld  / 2*0, 13*2 , 5*0 /
   data jfld  / 4*1, 16*0 /
   data katn  /'ATOM'/
   data ifind /'FIND'/
   data nfind /'NFIN'/
   data iiend /'END '/
   data ksear /'SEAR'/
   data ires  /'RES '/
   data ilres /'LRES'/
   data irres /'RRES'/
   data isear / 8/
   data igrap /'C   ','O   ','N   ','H   ','CA  ', &
         'HA  ','HA2 ','HA3 '/
   data isymb /'C   ','O   ','N   ','H   ','CT  ', &
         'H1  ','H1  ','H1  '/
   data iwild /'*   '/
   data jwild /'*   '/
   
   !     ----- RES CARD LISTS RESIDUE GROUPS IGRP(I) TO JGRP(I) ----
   
   !           IF 1ST VALUE IS NEGATIVE, THEN EVERY RESIDUE FROM IGRP(I)
   !           TO JGRP(I) IS CONSIDERED A SEPARATE GROUP
   !           IF 2ND VALUE   = 0, THEN SUBGROUP CONSISTS OF 1 RESIDUE
   !           ATOM CARD READS IN GROUPS OF ATOMS
   !           IF 2ND ATOM IN EACH PAIR  = 0 , THEN SUBGROUP CONSISTS OF
   !           1 ATOM RES AND ATOM CARDS MAY BE READ IN ANY ORDER
   !           END INPUT WITH AN "END " CARD
   
   !           ZERO NGRP BEFORE CALLING THIS ROUTINE
   !           ROUTINE WILL RETURN WITH THE NUMBER OF THE LAST GROUP READ
   

   ngrp = 0
   npdec = 0
   flres = .false.
   frres = .false.
   nf = nfu
   if (nf <= 0) nf=5
   
   !     ----- INITIALISE THE GROUP ARRAY -----
   
   do 100 i = 1,natom
      igroup(i) = 0
      if (konst) weit(i) = 0.0d0
   100 continue

   22 continue
   itime = 0
   lsign = 0
   izero = 0
   isrch = 0
   natmg = 0
   
   !       ----- READ DIFFERENT GROUPS -----
   
   read(nf,9208) (title1(k),k=1,19)
   if(title1(1) == iiend) then
      if(writeout) write(6, '(4x,a)') '----- END OF GROUP READ -----'
      goto 900
   end if
   ngrp = ngrp+1
   if(writeout) then
     write(6,'(4x,a,i5,a)') '----- READING GROUP ', ngrp, '; TITLE:'
     write(6,9218) (title1(k),k=1,19)
   end if
   
   !       ----- IF CONSTRAINED GROUPS READ THE WEIGHT FOR EACH GROUP -----
   !       ----- NOTE NO WEIGHT FOR TARGETED MD -----
   
   if (konst) then
      ifld(1) = 3
      ifld(2) = 0
      call rfree(ifld,ihol,ivar,fvar,nf,6)
      wt = fvar(1)
      if(writeout) write(6,9018) ngrp,wt
   end if


   10 continue
   !       ----- READ THE GROUP CARDS -----
   
   ifld(1) = 1
   ifld(2) = 2
   call rfree(ifld,ihol,ivar,fvar,nf,6)
   ktypg = ihol(1)
   k = 1
   do i = 1,7
      igrp(i) = ivar(k)
      jgrp(i) = ivar(k+1)
      k = k+2
   enddo
   if (ktypg == iiend) goto 16
   if (idecomp > 0 .and. &
         (ktypg /= ires .and. ktypg /= ilres .and. ktypg /= irres)) &
         goto 36
   if (idecomp > 0 .and. ktypg == ilres) flres = .true.
   if (idecomp > 0 .and. ktypg == irres) frres = .true.
   if (ktypg == nfind) then
      if(writeout) write(6,199)
      isrch = 0
      goto 10
   end if
   if (ktypg == ifind) then
      
      !         ----- FIND OPTION ... READ THE ATOM SPECIFICS -----
      
      if(writeout) write(6,200)
      do iii = 1,max_finds
         call rfree(jfld,ihol,ivar,fvar,nf,6)
         jgraph(iii) = ihol(1)
         jsymbl(iii) = ihol(2)
         jtree(iii) = ihol(3)
         jresnm(iii) = ihol(4)
         if(jgraph(iii) == ksear) exit
         if(writeout) write(6,202) jgraph(iii),jsymbl(iii),jtree(iii),jresnm(iii)
      enddo
      
      isrch = iii-1
      if(isrch > max_finds) then
        if(writeout) write(6,66) isrch
      end if
      !         ----- NOW READ IN RES AND ATOMS TO BE SEARCHED -----
      
      goto 10
   end if
   itime = itime+1
   
   if (idecomp > 0 .and. &
         (ktypg == ilres .or. ktypg == irres)) then
      
      !         ----- CHECK LRES or RRES CARD -----
      
      do 40 i = 1,7
         !           --- Sign does not play a role ---
         i1 = igrp(i)
         if (i1 == 0) goto 10
         if(i1 < 0) i1 = -i1
         if(i1 > nres) goto 36
         i2 = jgrp(i)
         if(i2 < 0) i2 = -i2
         if(i2 > nres) i2 = nres
         do 41 j = i1,i2
            if(idecomp > 2) then
               npdec = npdec + 1
               index(j) = npdec
               irespw(npdec) = j
            end if
            istart = ipres(j)
            iend = ipres(j+1)-1
            do 42 k = istart,iend
               !               --- Find backbone atoms ---
               do 50 l = 1,isear
                  jgraph(l) = igrap(l)
                  jresnm(l) = iwild
                  jsymbl(l) = isymb(l)
                  jtree(l)  = jwild
               50 continue
               isrch = isear
               call findp(k,j,lfind,nres,ipres,lbres,isymbl,itree,igraph)
               if (lfind == 1) then
                  if (ktypg == irres) then
                     jgroup(k) = nres + j
                  else !Ligand atom
                     jgroup(k) = -nres - j
                  end if
               else
                  !                 --- Store as sidechain atoms ---
                  if (ktypg == irres) then
                     jgroup(k) = j
                  else !Ligand atom
                     jgroup(k) = -j
                  end if
               end if
            42 continue
         41 continue
      40 continue
      ngrp = ngrp - 1
      goto 10
   end if
   
   if (ktypg /= katn) then
      
      !         ----- CHECK RES CARD -----
      
      !         ----- 1ST GROUP OF 1ST CARD MUST BE - IF ANY - NUMBERS ARE
      !               FOUND -----
      
      if(itime == 1.and.igrp(1) < 0) lsign = 1
      do 12 i = 1,7
         i1 = igrp(i)
         if (i1 == 0) goto 10
         i2 = jgrp(i)
         if(i2 > nres) i2 = nres
         if(i1 > 0.and.lsign == 1) goto 36
         if(i1 < 0.and.lsign == 0) goto 36
         if(i1 < 0) i1 = -i1
         if(i2 <= 0) i2 = i1
         if(lsign == 0) then
           if(writeout) write(6,14) ngrp,i1,i2
         end if
         do 13 j = i1,i2
            istart = ipres(j)
            iend = ipres(j+1)-1
            do 45 k = istart,iend
               if (isrch > 0) &
                  call findp(k,j,lfind,nres,ipres,lbres,isymbl,itree,igraph)
               if (isrch > 0.and.lfind == 0) goto 45
               igroup(k) = ngrp
               if(konst) weit(k) = wt
               natmg = natmg+1
            45 continue
            if(lsign == 1) then
              if(writeout) write(6,46) ngrp,j
            end if
            if(lsign == 1) ngrp = ngrp+1
         13 continue
      12 continue
      goto 10
   end if
   
   !       ----- ATOM TYPE CONSTRAINTS -----
   
   if(lsign == 1) goto 36
   if(writeout) write(6,51) ngrp
   do 17 i = 1,7
      i1 = igrp(i)
      if(i1 < 0) goto 36
      if(i1 == 0) goto 10
      i2 = jgrp(i)
      if(i2 > natom) i2 = natom
      if(i2 <= 0) i2 = i1
      if(writeout) write(6,52) i1,i2
      do 18 j = i1,i2
         if(isrch > 0) &
            call findp(j,izero,lfind,nres,ipres,lbres,isymbl,itree,igraph)
         if(isrch > 0.and.lfind == 0) goto 18
         natmg = natmg+1
         igroup(j) = ngrp
         if(konst) weit(j) = wt
      18 continue
   17 continue
   goto 10
   
   16 if(lsign == 1) ngrp = ngrp-1
   if(itime == 0) ngrp = ngrp-1
   !       IF(ISRCH.GT.0) WRITE(6,199)
   if(writeout) write(6,222) natmg
   goto 22
   
   36 continue
   if(writeout) write(6,127) ktypg,(igrp(i),jgrp(i),i = 1,7)
   goto 10
   
   !     ----- ALL GROUPS ARE READ RETURN -----
   
   900 continue
   if (konst.or.dotgtmd) then
      
      !       ----- GATHER ALL THE CONSTRAINED ATOMS TOGETHER -----
      
      natc = 0
      do 920 i = 1,natom
         if(igroup(i) <= 0) goto 920
         natc = natc+1
         igroup(natc) = i
         
         ! WEIT will not be used for targeted MD
         
         weit(natc) = weit(i)
      920 continue
      
      !       ----- do not PRINT THE HISTORY OF CONSTRAINED ATOMS -----
      ! if(writeout) write(6,9108)
      ! do i = 1,natc
      !    j = igroup(i)
      !    j3 = 3*j-3
      !    if(writeout) write(6,9118) i,j,igraph(j),isymbl(j),itree(j), &
      !          (xc(j3+k),k=1,3), weit(i)
      ! end do
      
   else if (idecomp > 0) then
      
      !       ----- Special treatment for energy decomposition -----
      
      if(.not.flres .and. .not.frres) then
         !         --- Assign all atoms to "Protein" ---
         !         --- Set all residues to be printed ---
         do 43 j = 1,nres
            if(idecomp > 2) then
               npdec = npdec + 1
               index(j) = npdec
               irespw(npdec) = j
            end if
            istart = ipres(j)
            iend = ipres(j+1)-1
            do 44 k = istart,iend
               igroup(k) = 1
               !              --- Find backbone atoms ---
               do 53 l = 1,isear
                  jgraph(l) = igrap(l)
                  jresnm(l) = iwild
                  jsymbl(l) = isymb(l)
                  jtree(l)  = jwild
               53 continue
               isrch = isear
               call findp(k,j,lfind,nres,ipres,lbres,isymbl,itree,igraph)
               if (lfind == 1) then
                  jgroup(k) = nres + j
               else
                  !                --- Store as sidechain atoms ---
                  jgroup(k) = j
               end if
            44 continue
         43 continue
      end if

      !       --- Print assignment of atoms ---
      do 47 j = 1,nres
         istart = ipres(j)
         iend = ipres(j+1)-1
         do 48 k = istart,iend
            if(writeout) write(6,'(a,i6,a,i6,a,2i6)') 'Atom ',k,' (',j,') : ',jgroup(k),igroup(k)
         48 continue
      47 continue
      
   else if (.not.belly) then
      
      !       ----- PUT THE ATOMS WHICH ARE NOT IN THE DEFINED GROUPS
      !             AS THE LAST GROUP -----
      
      misc = .false.
      do 820 i = 1,natom
         if(igroup(i) /= 0) goto 820
         misc = .true.
         igroup(i) = ngrp+1
      820 continue
      if(misc) ngrp = ngrp+1
      !       IF(MISC) WRITE(6,9308) NGRP
   end if
   
   ! for targeted MD, support only 1 group for now
   
   if (dotgtmd.and.ngrp > 1) then
      write (6,9200)
      write (6,9205)
      call mexit(6, 1)
      9200 format ("ERROR IN TARGETED MD GROUP INPUT (ITGTMD=1)")
      9205 format ("ONLY 1 GROUP ALLOWED")
   end if

   !  11 FORMAT(A4,1X,7(2I5))
   199 format(6x,'END OF ATOM SPECIFICATION',/)
   200 format(6x,'ALL ATOMS THAT MEET 1 OF THE FOLLOWING', &
         ' SPECIFICATIONS WILL BE INCLUDED IN GROUP BELOW',/)
   !  61 FORMAT(A4,A2,A1,A4)
   202 format(6x,'GRAPH NAME  = ',a4,2x,'SYMBOL  = ',a2,4x, &
         'TREE SYMBOL  = ',a1,5x,'RESIDUE TYPE  = ',a4,/)
   66 format(6x,'**** NUMBER OF FIND CARDS  = ',i5,2x, &
         'IS TOO BIG ******',/)
   14 format(' GRP',i5,' RES',i5,' TO ',i5)
   46 format(6x,'GROUP',i5,3x,'CONSISTS OF RESIDUE',i5,/)
   51 format(6x,'GROUP',i5,3x,'CONSISTS OF ATOMS -',/)
   52 format(35x,i5,2x,'TO',i5)
   222 format(6x,'Number of atoms in this group  = ',i5)
   127 format(6x,'***PROBLEMS WITH GROUP',a4,14i5,'*******',/)
   9018 format(/5x,'GROUP',i5,' HAS HARMONIC CONSTRAINTS',f12.5)
   9118 format(i5,i6,1x,a4,1x,2a4,3f10.4,f12.5)
   9108 format(/ /10x,'HISTORY OF CONSTRAINED ATOMS',/ /)
   9208 format(20a4)
   9218 format(1x,20a4)
   !9308 FORMAT(/5X,'THE GROUP ',I4, ' CONTAINS ALL ATOMS NOT DEFINED',
   !    +       ' AS GROUPS BY THE INPUT',/)
   return

   contains   

   subroutine findp(iatom,ires,iloc,nres,ipres,lbres,isymbl,itree,igraph)
   
   implicit none
   integer:: iatom, iloc, ipres, ires, n, nres

   character(len=4) lbres,isymbl,itree,igraph,iwild,jwild

   !     Rev A mod (G. Seibel, Apr 89)
   !     Changed iblank, jblank to iwild, jwild, to give the wildcard
   !     functionality promised by the doc.
   !     isymbl() dimensioned. (this was long-standing bug)
   
   dimension ipres(*),lbres(*),itree(*),igraph(*),isymbl(*)
   
   !     ----- CHECKS IF A GIVEN ATOM HAS CERTAIN CHARACTERISTICS -----
   
   iwild = '*   '
   jwild = '*   '
   iloc = 0
   if(ires == 0) call findrs(iatom,ires,nres,ipres)
   do 10 n = 1,isrch
      if((jresnm(i) /= iwild).and.(jresnm(i) /= lbres(ires))) goto 10
      if((jgraph(i) /= iwild).and.(jgraph(i) /= igraph(iatom))) goto 10
      if((jtree(i) /= jwild).and.(jtree(i) /= itree(iatom))) goto 10
      if((jsymbl(i) /= jwild).and.(jsymbl(i) /= isymbl(iatom))) goto 10
      iloc = 1
      goto 20
   10 continue
   20 continue
   return
   end subroutine findp 
!-----------------------------------------------------------------------
end subroutine rgroup 
!-----------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine findrs here]
subroutine findrs(numa,ires,nres,ipres)
   
   dimension ipres(*)
   
   if(numa <= 0) goto 11
   if(numa >= ipres(nres)) ires = nres
   if(numa >= ipres(nres)) return
   im = nres-1
   do 10 i = 1,im
      if(numa >= ipres(i).and.numa < ipres(i+1)) goto 12
   10 continue
   11 write(6,100) numa
   goto 200
   12 ires = i
   return
   200 continue
   100 format(/2x,'PROBLEMS FINDING RESIDUE OF ATOM ',i5)
   call mexit(6, 1)
end subroutine findrs 
