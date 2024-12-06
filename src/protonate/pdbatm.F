c
c======================================================================c
c
c Function that reads and writes PDB files.
c JPH added Feb 97 - also write PDBQ files - used by Autodock.
c
c
      subroutine pdbatm( maxatm, ird, ipr, iwr,
     .                   natom, atmnam, resnam, chid, resno, cord,
     .                   iterm )
c
      integer maxatm, ird, ipr, iwr, natom
      character*4  atmnam(maxatm)
      character*3  resnam(maxatm)
      character*1  chid(maxatm)
      integer      resno(maxatm), iterm(maxatm)
      real         cord(3,maxatm), chrg(1)
c
c  ---  Read  ---
c
      if ( ird .gt. 0 ) 
     .   call pdbin( ird, maxatm, natom, atmnam, resnam, chid,
     .               resno, cord, iterm )
c
c  ---  Print  ---
c
      if ( ipr .gt. 0 )
     .   call pdbout( ipr, 2, maxatm, natom, atmnam, resnam, chid,
     .                resno, cord, iterm, 0, 0, chrg )
c
c  ---  Print  ---
c
      if ( iwr .gt. 0 )
     .   call pdbout( iwr, 1, maxatm, natom, atmnam, resnam, chid,
     .                resno, cord, iterm, 0, 0, chrg )
c
      return
      end
c
c======================================================================cc
      subroutine pdbin( iun, maxatm, natom, atmnam, resnam, chid,
     .                   resno, cord, iterm )
c
c Read a PDB file
c
c
c
      integer iun, maxatm, natom, resno(maxatm), iterm(maxatm)
      character*4  atmnam(maxatm)
      character*3  resnam(maxatm), tag
      character*1  chid(maxatm)
      real         cord(3,maxatm)
      character*80 line
      equivalence  (line,tag)
c
      natom = 0
  100 continue
      read( iun, 200, end=500 ) line
  200 format( a80 )
      if ( tag .eq. 'END' ) go to 500
      if ( tag .eq. 'TER' ) iterm(natom) = 1
      if ( tag .eq. 'HET' .or. tag .eq. 'ATO' ) then
         natom = natom + 1
         if ( natom .gt. maxatm ) then
            write(0,*) 'pdbatm: maxatm exceeded'
            call mexit(0,1)
         end if
         iterm( natom ) = 0
         read( line, 300 ) atmnam(natom), resnam(natom), chid(natom),
     .                     resno(natom), cord(1,natom),
     .                     cord(2,natom), cord(3,natom)
  300    format( 12x, a4, 1x, a3, 1x, a1, i4, 4x, 3f8.0, 16x )
      endif
      go to 100
c
c
  500 continue
      return
      end
c
c======================================================================c
c
      subroutine pdbout( iun, ipos, maxatm, natom, atmnam, resnam, chid,
     .                   resno, cord, iterm, iq, iqr, chrg )
c
c Write a PDB or PDBQ file
c
c
c
      integer iun, ipos, maxatm, natom, iq, iqr, iend,
     .        resno(maxatm), iterm(maxatm), i
      character*4  atmnam(maxatm)
      character*3  resnam(maxatm)
      character*1  chid(maxatm)
      real         cord(3,maxatm), chrg(maxatm), q1
      character*80 line
c
      if (iun .gt. 0 ) then
         iend = ipos + 77
         if ( iq .ne. 0 ) iend = ipos + 75
         do 100 i = 1, natom
            line = ' '
            write( line(ipos:), 200 ) i, atmnam(i), resnam(i), chid(i),
     &         resno(i),cord(1,i), cord(2,i), cord(3,i), atmnam(i)(1:2)
c
c Add the charges for a PDBQ file
c
            if ( iq .ne. 0 ) then
c
c NB Remove conversion factor from PARM
c
               q1 = chrg(i) / 18.2223d0
c
c iqr describes the format for the charges
c     = 0  write charges in F6.3 format - truncating
c       1  write charges in F6.3 format - rounding
c      -1  write charges in F8.5 format
c
               if ( iqr .lt. 0 ) then
                  write( line(ipos+70:), 400 ) q1
               else
                  if ( iqr .gt. 0 ) q1 = q1 + .0005
                  write( line(ipos+70:), 500 ) q1
               end if
            end if
            write( iun, 1000 ) line(1:iend)
            if ( iterm(i) .eq. 1 ) then
               line(ipos:) = 'TER'
               write( iun, 1000 ) line(1:ipos+2)
            endif
  100    continue
         line(ipos:) = 'END'
         write( iun, 1000 ) line(1:ipos+2)
      end if
c
  200 format( 'ATOM  ', i5, 1x, a4, 1x, a3, 1x, a1, i4, 4x, 3f8.3, 
     &        '  1.00  0.00          ', a2 )
  400 format( f8.5 )
  500 format( f6.3 )
 1000 format( a )
c
      return
      end
