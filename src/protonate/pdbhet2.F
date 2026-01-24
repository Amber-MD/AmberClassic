c======================================================================c
c
      subroutine pdbhet2 (ird,ipr,iwr, natom,atmnam,resnam,cord
     &                                ,occ,temp)
      character*80 line
      character*4  atmnam(1)
      character*3  resnam, tag
      dimension    cord(3,1),occ(1),temp(1)
      equivalence  (line,tag)
c
c                                                ---  Read  ---
      if (ird.le.0) go to 6
      natom = 0
   51 read  (ird,501,end=6) line
  501 format (a80)
      if (tag.eq.'END') go to 6
      if (tag.eq.'HET' .or. tag.eq.'ATO') then
         natom = natom + 1
         read   (line,511) atmnam(natom),resnam,(cord(j,natom),j=1,3)
     &                    ,occ(natom),temp
  511    format (12x,a4,1x,a3,10x,3f8.0,2f6.0)
      endif
      go to 51
c
c                                                ---  Print  ---
    6 if (ipr.le.0) go to 8
      do 610 i=1,natom
         write  (ipr,601) i,atmnam(i),resnam,(cord(j,i),j=1,3),
     &                    occ(i),temp(i)
  601    format (' HETATM',i5,1x,a4,1x,a3,2x,'   1',4x,3f8.3,2f6.2)
  610 continue
c
c                                                ---  Write  ---
    8 if (iwr.le.0) return
      do 810 i=1,natom
         write  (iwr,801) i,atmnam(i),resnam,(cord(j,i),j=1,3),
     &                    occ(i),temp(i)
  801    format ('HETATM',i5,1x,a4,1x,a3,2x,'   1',4x,3f8.3,2f6.2)
  810 continue
      write  (iwr,809)
  809 format ('TER   '/'END   ')
      return
      end
