      program protonate
      implicit none
c
c -- Function:-
c     Add all protons to a protein structure according to formulae
c     specified in the PROTON_INFO file
c
c -- Created:-
c     This routine is (loosely) based on code originated by K. Cross.
c     The bulk of the present implementation was done by G.P. Gippert,
c     with some input by D.A. Case.
c
c -- Modified:-
c     Feb 1997  JPH  Added switches to specify the file names
c                    Combined UCSF and TSRI versions
c                    Fixed 'atfind' to check chain id as well as residue,
c                    and residue insertion codes at various places
c
#include "protonate.h"
c
      character*(NLINSZ) inline, line
      character*3 cdresnam(NR), cpro, creslast
      character*4 cmnatm(3), ctemp
      character*4 camide, camidet(3)
      character*1 ctemp1
      logical first, header, aminot, nterm(MAXAT), ter(MAXAT)
      integer minres, maxres, ksave, lname, ios,
     .        iat, icout, ic, iar, nscat, i, j, k, l, m, n,
     .        nrs, nhyd, lres
c
      real protp(3,4)
c
      data cmnatm / ' N  ', ' CA ', ' C  ' /
      data cpro  / 'PRO' /
      data first, aminot / .true., .true. /
      data nterm / MAXAT*.false. /
      data ter   / MAXAT*.false. /
c
c -- initialise the unit numbers and read arguments etc
c -- process output, printing options etc
c
      call start( line )
c
      minres = 0
      maxres = 999999
c
c  -- read PDB input file from stdin
c
   30 continue
      iat = 1
      do icout = 1, MAXAT
        read( nfilin, '(a80)', iostat=ios ) inline
        if( ios>0 ) go to 40
        if ( inline(1:3) .eq. 'END' ) go to 40
c
c -- (following line assumes that the alternate conformation flag
c     should be either blank or "A")
c
        if ( ( inline(1:4) .eq. 'ATOM' .or.
     .         inline(1:6) .eq. 'HETATM' ) .and.
     .       ( inline(17:17) .eq. ' ' .or.
     .         inline(17:17) .eq. 'A' ) ) then
          read( inline, 35 ) ckey(iat), cnamat(iat), cresnam(iat),
     .                       chainid(iat), ires(iat), crescod(iat),
     .                       (coord(ic,iat), ic=1,3)
   35     format( A6, 6X, A4, 1X, A3, 1X, A1, I4, A1, 3X, 3F8.3 )
c         if ( cresnam(iat) .eq. 'HIS' ) cresnam(iat) = 'HIE'
#if 0
c
c  -- left justify the input atom name, as Amber does internally:
c
          if( cnamat(iat)(1:1).eq.' ') then
            cnamat(iat)(1:1) = cnamat(iat)(2:2)
            cnamat(iat)(2:2) = cnamat(iat)(3:3)
            cnamat(iat)(3:3) = cnamat(iat)(4:4)
            cnamat(iat)(4:4) = ' '
          end if
c
          lname = 3
          if( cnamat(iat)(3:3).eq.' ' ) lname = 2
          if( cnamat(iat)(2:2).eq.' ' ) lname = 1
c
          if( cnamat(iat)(1:1).eq.'1' .or. cnamat(iat)(1:1).eq.'2' .or.
     .        cnamat(iat)(1:1).eq.'3' ) then
            ctemp1 = cnamat(iat)(1:1)
            do i=1,lname
              cnamat(iat)(i:i) = cnamat(iat)(i+1:i+1)
            end do
            cnamat(iat)(lname+1:lname+1) = ctemp1
          end if
          if( cnamat(iat)(3:3).eq.'*' ) cnamat(iat)(3:3) = ''''
#endif
c
c  -- remember positions of protons in input file for later checking:
c
          if ( cnamat(iat)(1:1).eq.'H' .or. 
     .         cnamat(iat)(1:2).eq.' H' ) then
            isoldp(iat) = .true.
            if( cresnam(iat).eq.'WAT' .or. cresnam(iat).eq.'HOH')
     .         isoldp(iat) = .false.
            matched(iat) = .false.
          else
            isoldp(iat) = .false.
#if 0
c
c           --- try to recognize heavy atoms either with or without
c               a blank in column 13:
c
            if( cnamat(iat)(1:1).ne.' ' ) then
              ctemp(1:1) = ' '
              ctemp(2:4) = cnamat(iat)(1:3)
              cnamat(iat) = ctemp
            end if
#endif
          end if
          if ( cnamat(iat)(1:2) .eq. ' N' .and. aminot ) then
            if ( iat .eq. 1 ) then
              nterm(iat) = .true.
            else if ( cresnam(iat-1) .ne. 'ACE' ) then
              nterm(iat) = .true.
            end if
            aminot = .false.
          end if
          iat = iat + 1
        else if ( inline(1:3) .eq. 'TER' ) then
          aminot = .true.
          ter(iat) = .true.
        else if ( inline(1:3) .eq. 'END' ) then
          go to 40
        endif
      end do
      call maxerr( iat )
c
   40 continue
      natom = iat - 1
      if ( natom .gt. MAXAT ) then
         call maxerr( iat )
      end if
      if ( natom .le. 0 ) go to 900
c
c -- read file containing information about protons in each
c    type of residue
c
      if ( .not. first ) go to 80
      first = .false.
      read( ndatin, 1000 ) camide, camidet
      do 60 iar = 1, NR
        read( ndatin, 1010, end=70 ) cdresnam(iar), nscat
        read( ndatin, 1020 ) (cdatname(iar,i), i=1,nscat)
        read( ndatin, 1030 ) nprot(iar)
        n = nprot(iar)
        if ( n .gt. NHSET ) then
           write( nerrout, * ) 'number of sets of protons is too great!'
           call mexit( 0, 1 )
        end if
        do 50 i = 1, n
           read( ndatin, 1040 ) iprot(iar,i), (indat(iar,i,j), j=1,4),
     .                          (cdhname(iar,i,k), k=1,3)
           m = indat(iar,i,1)
           cdatnamh(iar,i) = cdatname(iar,m)
   50   continue
   60 continue
      nrs = NR
      go to 80
   70 continue
      nrs = iar-1
 1000 format( 15X, 4(A4,1X) )
 1010 format( A3, I3 )
 1020 format( 20A4 )
 1030 format( 3I2 )
 1040 format( I2, 4I3, 3(1X,A4) )
c
c -- proceed through protein data file atom by atom.
c    attach protons to appropriate atoms according to information
c    in stdprt file.
c    write proton data on new file.
c
c -- hydros does the real work of finding the proton positions
c
   80 continue
      creslast = '   '
      icount   = 1
      ilast    = 1
      iresct   = 0
      iresp    = ires(1)
      chainp   = chainid(1)
      cinsert  = crescod(1)
c
      do 171 j = 1, natom
c
c -- skip old protons; we'll deal with them later:
c
        if ( isoldp(j) ) go to 170
        if ( ires(j) .lt. minres .or.
     .       ires(j) .gt. maxres ) go to 900
c
        if ( atout )
     .    call writeout( ckey(j), cnamat(j), cresnam(j), chainid(j),
     .                   ires(j), crescod(j), coord(1,j), coord(2,j),
     .                   coord(3,j) )
c
c -- Check for N atom
c
        if ( cnamat(j) .ne. cmnatm(1) ) then
c
c -- side-chain protons and ca protons
c
          if ( creslast .eq. cresnam(j) ) then
            call asch( ksave, j, ires(j), chainid(j), crescod(j) )
          else
c
c -- find the correct residue in the stdprt file
c
            do 100 k = 1, nrs
              if ( cresnam(j) .ne. cdresnam(k) ) go to 100
              creslast = cresnam(j)
              ksave = k
              call asch( ksave, j, ires(j), chainid(j), crescod(j) )
              go to 170
  100       continue
          end if
        else
c
c -- n-h protons on the backbone
c
          if ( cresnam(j) .eq. cpro ) go to 170
c
          do 120 l = 1, 3
            catom(l) = cmnatm(l)
  120     continue
c
          if ( nterm(j) ) then
            call atfind( 3, ires(j), chainid(j), crescod(j) )
            call hydros( 10, 3, nhyd, protp )
            do 130 n =1, nhyd
              ctemp = camidet(n)
              call writeout( ckey(j), ctemp, cresnam(j), chainid(j),
     .                       ires(j), crescod(j), protp(1,n),
     .                       protp(2,n), protp(3,n) )
  130       continue
          else
            lres = ires(j) - 1
            catom(1) = cmnatm(3)
            call atfind( 1, lres, chainid(j), crescod(j) )
            do 150 ic = 1, 3
              scoord(ic,3) = scoord(ic,1)
  150       continue
            catom(1) = cmnatm(1)
            call atfind( 2, ires(j), chainid(j), crescod(j) )
            call hydros( 8, 3, nhyd, protp )
            do 160 n = 1, nhyd
              ctemp = camide
              call writeout( ckey(j), ctemp, cresnam(j), chainid(j),
     .                       ires(j), crescod(j), protp(1,n),
     .                       protp(2,n), protp(3,n) )
  160       continue
          end if
        end if
  170   if( ter(j) ) write(nfilout,'(a3)') 'TER'
  171 continue
c
      icount = icount - 1
      if ( ndbgout .ne. 0 ) then
         iresct = iresct + 1
         write( ndbgout, '(3I7)' ) iresct, iresp, icount - ilast
      end if
c
c -- print out summary of unmatched protons to stderr:
c
      header = .true.
      do 500 i = 1, natom
        if ( isoldp(i) .and. .not. matched(i) ) then
          if ( header ) then
            write( nerrout, 510 ) line, line
            write( nerrout, 520 )
            header = .false.
          end if
          write( nerrout, 530 ) cnamat(i), cresnam(i), ires(i),
     .                          coord(1,i), coord(2,i), coord(3,i)
        end if
  500 continue
  510 format( A40, A40 )
  520 format( 'Here are the mystery protons from input file:' / )
  530 format( A4, 1X, A3, I4, 3F7.2 )
c
c --- loop back and read another pdb file (if present):
c
      write( nfilout, '(a3)' ) 'END'
      go to 30
  900 continue
      end
c
c
      subroutine start( line )
      implicit none
c
#include "protonate.h"
c
      character*(*) line
      character*(NLINSZ) arg
      integer iargc, iarg, nargs, lendf, len_trim
c
c     external iargc, getarg, len_trim
c
      line     = '----------------------------------------'
      keep     = .false.
      ipopt    = .false.
      mismch   = .false.
      lprint   = .false.
      atout    = .true.
      call getenv( 'AMBERHOME', datafile )
      lendf = len_trim( datafile )
      if (lendf.gt.0 ) then
        datafile = datafile(1:lendf) // '/dat/PROTON_INFO'
      else
        datafile = 'PROTON_INFO'
      end if
      infile   = '-'
      outfile  = '-'
      logfile  = '-'
      dbgfile  = '-'
c
      nerrout  =  0
      ndbgout  =  0
      nfilin   =  5
      nfilout  =  6
      nprmout  =  7
      nedtout  =  8
      nlnkout  =  9
      ndatin   = 10
c
#ifdef HITACHI
      iarg = 1
#else
      iarg = 0
#endif
      nargs = iargc()
      if (nargs.eq.iarg) go to 20
   10 continue
          iarg = iarg + 1
          call getarg(iarg,arg)
          if ( arg .eq. '-f' .or. arg .eq. '-force' ) then
            ipopt = .true.
          else if ( arg .eq. '-h' .or. arg .eq. '-hyd_only' ) then
            atout = .false.
          else if ( arg .eq. '-k' .or. arg .eq. '-keep' ) then
            keep = .true.
          else if ( arg .eq. '-m' .or. arg .eq. '-mismatch' ) then
            mismch = .true.
          else if ( arg .eq. '-p' .or. arg .eq. '-print' ) then
            lprint = .true.
          else if ( arg .eq. '-r' .or. arg .eq. '-renumber' ) then
            lrenum = .true.
          else if ( arg .eq. '-d' .or. arg .eq. '-data' ) then
            iarg = iarg + 1
            call getarg(iarg,datafile)
          else if ( arg .eq. '-i' .or. arg .eq. '-in' ) then
            iarg = iarg + 1
            call getarg(iarg,infile)
          else if ( arg .eq. '-o' .or. arg .eq. '-out' ) then
            iarg = iarg + 1
            call getarg(iarg,outfile)
          else if ( arg .eq. '-l' .or. arg .eq. '-log' ) then
            iarg = iarg + 1
            call getarg(iarg,logfile)
          else if ( arg .eq. '-debug' ) then
            iarg = iarg + 1
            call getarg(iarg,dbgfile)
          else
            write(nerrout,*) 'Bad argument: ',arg(1:len_trim(arg))
            write(nerrout,*)
     .         'Usage: protonate [-fhkmpr] [-d datafile]',
     .         ' [-i input-pdb-file] [-o output-pdb-file] [-l logfile]'
            write(nerrout,*)
     .         '   -f to force write of atoms found (debugging)'
            write(nerrout,*)
     .         '   -h to write ONLY hydrogens to output file'
            write(nerrout,*)
     .         '   -k to keep original coordinates of matched protons'
            write(nerrout,*)
     .         '   -m to list mismatched protons'
            write(nerrout,*)
     .         '   -p to print proton substitutions'
            write(nerrout,*)
     .         '   -r to renumber residues starting at 1'
            write(nerrout,*)
     .       '   -d to specify datafile',
     .       ' (default is $AMBERHOME/dat/PROTON_INFO)'
            write(nerrout,*)
     .         '   -i to specify input file (default is stdin)'
            write(nerrout,*)
     .         '   -o to specify output file (default is stdout)'
            write(nerrout,*)
     .         '   -l to specify logfile (default is stderr)'
            call mexit(0,1)
          end if
      if ( iarg .lt. nargs ) go to 10
   20 continue
      call    amopen( ndatin,  datafile, 'O', 'F', 'R' )
      if ( infile .ne. '-' ) then
         nfilin = 11
         call amopen( nfilin,  infile,   'O', 'F', 'R' )
      end if
      if ( outfile .ne. '-' ) then
         nfilout = 12
         call amopen( nfilout, outfile,  'U', 'F', 'W' )
      end if
      if ( logfile .ne. '-' ) then
         nerrout = 13
         call amopen( nerrout, logfile,  'U', 'F', 'W' )
      end if
      if ( dbgfile .ne. '-' ) then
         ndbgout = 14
         call amopen( ndbgout, dbgfile,  'U', 'F', 'W' )
      end if
c
      return
      end
c
      subroutine asch( k, j, kres, chain, resc )
      implicit none
c
c -- for a given atom belonging to a known residue, check whether it
c    has attached proton(s).  If so, call appropriate routines to
c    determine the position(s) of the proton(s)
c
#include "protonate.h"
c
      integer k, j, kres
      character*1 chain, resc
c
      character*4 atoml
      real protp(3,4)
      integer l, m, mm, n, nat, nhyd, npt, kind
c
      npt=nprot(k)
c
c!!      write(nerrout,'(3i5,1H:,a4,1H:,a3,1H:)') k,j,kres,
c!!     .                                         cnamat(j),cdatnamh(k,1)
      do 30 l = 1, npt
          if ( cdatnamh(k,l) .ne. cnamat(j) ) go to 30
          kind = iprot(k,l)
          nat = 3
          if ( kind .eq.  3 .or. kind .eq. 4 .or.
     .         kind .eq.  6 .or. kind .eq. 9 .or.
     .         kind .eq. 12 ) nat = 4
          do 10 m = 1, nat
             mm = indat(k,l,m)
             catom(m) = cdatname(k,mm)
   10     continue
          call atfind( nat, kres, chain, resc )
          call hydros( kind, nat, nhyd, protp )
          do 20 n = 1, nhyd
            atoml = cdhname(k,l,n)
            call writeout( ckey(j), atoml, cresnam(j), chainid(j),
     .                     kres, resc, protp(1,n), protp(2,n),
     .                     protp(3,n) )
   20     continue
   30 continue
      return
      end
c----------------------------------------------------------------------
      subroutine atfind( nat, npres, chain, resc )
      implicit none
c
c -- obtain coordinates from data file of specified set of non-proton atoms
c
c -- Modified:-
c     Feb 1997 JPH  Check for chain identity as well as residue numbers
c                   Also check for residue insertion codes
c
c
#include "protonate.h"
c
      integer nat, npres
      character*1 chain, resc
c
      integer mink, iat, k, ic
c
      mink = 1
      do 40 iat = 1, nat
         do 20 k = mink, natom
            if ( ires(k) .eq. npres .and.
     .           crescod(k) .eq. resc .and.
     .           chainid(k) .eq. chain .and.
     .           cnamat(k) .eq. catom(iat) ) then
               do 10 ic = 1, 3
                  scoord(ic,iat) = coord(ic,k)
   10          continue
               mink = k - 100
               if ( mink .lt. 1 ) mink = 1
               go to 30
            end if
   20    continue
c
   30    continue
         if ( ipopt )
     .      call writeout( ckey(iat), catom(iat), cresnam(iat),
     .                     chainid(iat), npres, resc, scoord(1,iat),
     .                     scoord(2,iat), scoord(3,iat) )
   40 continue
c
      return
      end
c----------------------------------------------------------------------
      subroutine hydros( ip, na, nhyd, protp )
      implicit none
c
#include "protonate.h"
c
c -- calculate proton coordinates given the proton type and the
c    coordinates of required non-proton atoms
c
      logical first
      integer ip, na, nhyd, ierror, ic
      real rad, c1, c2, c3, c4, c5, c6, ar5, bl, deta
      real protp(3,4), r1(3), r2(3), r3(3), r4(3), r5(3), r6(3)
      real dh(13), det(3)
      intrinsic sqrt, cos, abs
c
c -- following values taken from CRC handbook, 59th edition,
c      p. F-215:
c
      data dh /1.084, 1.073, 1.070, 1.084, 1.101, 1.084,
     .         0.96,  1.01,  1.01,  1.014, 1.073, 1.070, 1.330/
      data first / .true. /
c
      save first, dh, rad, c1, c2, c3, c4, c5, c6
c
      if ( first ) then
         first = .false.
         rad = 3.141593 / 180.0
         c1 = cos(54.73*rad)
         c2 = cos(35.27*rad)
         c3 = cos(30.0*rad)
         c4 = cos(60.0*rad)
         c5 = cos(70.53*rad)
         c6 = cos(19.47*rad)
      end if
c
      go to (10,10,90,120,140,160,180,10,120,140,10,90,192),ip
      call exit_error('error on computed goto #1 in subr hydros',
     .        ip,na,nhyd,protp)
c
c -- case 1    benzene type
c
   10 continue
      call unit(1,2,r1,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #1.100a, bond r1 has zero length.  ',
     .              ip,na,nhyd,protp)
      end if
      call unit(1,3,r2,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #1.100b, bond r2 has zero length.  ',
     .             ip,na,nhyd,protp)
      end if
      do 20 ic=1,3
         r5(ic)=r1(ic) + r2(ic)
   20 continue
      ar5=sqrt(r5(1)*r5(1)+r5(2)*r5(2)+r5(3)*r5(3))
      do 30 ic=1,3
         r5(ic)=r5(ic)/ar5
   30 continue
      if(ip.eq.2.or.ip.eq.11) go to 70
      if(ip.eq.8) go to 50
c
c -- benzene c-h
c
      do 40 ic=1,3
         protp(ic,1)=scoord(ic,1) + r5(ic)*dh(1)
   40 continue
      go to 210
c
c -- aromatic or backbone n-h
c
   50 continue
      do 60 ic=1,3
         protp(ic,1)=scoord(ic,1) +  r5(ic)*dh(8)
   60 continue
      go to 210
c
c -- case 2    ch2 type
c
   70 continue
      call mult(r6,r1,r2)
      bl=dh(2)
      if(ip.eq.11) bl=dh(11)
      do 80 ic=1,3
         protp(ic,1)=scoord(ic,1) + (r5(ic)*c1+r6(ic)*c2)*bl
         protp(ic,2)=scoord(ic,1) + (r5(ic)*c1-r6(ic)*c2)*bl
   80 continue
      go to 210
c
c -- case 3    tetrahedral c-h
c
   90 continue
      call unit(1,2,r1,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #3.900a, bond r1 has zero length.  ',
     .                ip,na,nhyd,protp)
      end if
      call unit(1,3,r2,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #3.900b, bond r2 has zero length.  ',
     .              ip,na,nhyd,protp)
      end if
      call unit(1,4,r3,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #3.900c, bond r3 has zero length.  ',
     .              ip,na,nhyd,protp)
      end if
      deta=r1(1)*(r2(2)*r3(3)-r2(3)*r3(2)) - r1(2)*(r2(1)*r3(3)
     .          -r2(3)*r3(1)) + r1(3)*(r2(1)*r3(2)-r2(2)*r3(1))
      if(abs(deta).lt.1.0e-6) go to 200
      det(1)=r2(2)*r3(3) - r2(3)*r3(2) - r1(2)*(r3(3)-r2(3))
     .           + r1(3)*(r3(2)-r2(2))
      det(2)=r1(1)*(r3(3)-r2(3)) - r2(1)*r3(3) + r2(3)*r3(1)
     .           + r1(3)*(r2(1)-r3(1))
      det(3)=r1(1)*(r2(2)-r3(2)) - r1(2)*(r2(1)-r3(1))
     .           + r2(1)*r3(2) - r2(2)*r3(1)
      do 100 ic=1,3
         r5(ic)=det(ic)/deta
  100 continue
      ar5=sqrt(r5(1)*r5(1)+r5(2)*r5(2)+r5(3)*r5(3))
      bl=dh(3)
      if(ip.eq.12) bl=dh(12)
      do 110 ic=1,3
         protp(ic,1)=scoord(ic,1) + r5(ic)*bl/ar5
  110 continue
      go to 210
c
c -- case 4    planar ethylene ch2 or nh2
c
  120 continue
      call unit(1,2,r1,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #4.1200a, bond r1 has zero length. ',
     .              ip,na,nhyd,protp)
      end if
      call unit(4,3,r2,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #4.1200b, bond r2 has zero length. ',
     .              ip,na,nhyd,protp)
      end if
      call mult(r3,r1,r2)
      call mult(r4,r1,r3)
      bl=dh(4)
      if(ip.eq.9) bl=dh(9)
      do 130 ic=1,3
         protp(ic,1)=scoord(ic,1) + (r1(ic)*c4 + r4(ic)*c3)*bl
         protp(ic,2)=scoord(ic,1) + (r1(ic)*c4 - r4(ic)*c3)*bl
  130 continue
      go to 210
c
c -- case 5    methyl protons or nh3+ assumed to be in staggered config
c
  140 continue
      call unit(1,2,r1,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #5.1400a, bond r1 has zero length. ',
     .              ip,na,nhyd,protp)
      end if
      call unit(3,2,r2,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #5.1400b, bond r2 has zero length. ',
     .              ip,na,nhyd,protp)
      end if
      call mult(r3,r2,r1)
      call mult(r4,r1,r3)
      bl=dh(5)
      if(ip.eq.10) bl=dh(10)
      do 150 ic=1,3
          protp(ic,1)=scoord(ic,1) +
     .                (r3(ic)*c3*c6 + r1(ic)*c5 + r4(ic)*c4*c6)*bl
          protp(ic,2)=scoord(ic,1) -
     .                (r3(ic)*c3*c6 - r1(ic)*c5 - r4(ic)*c4*c6)*bl
          protp(ic,3)=scoord(ic,1) + (r1(ic)*c5 - r4(ic)*c6)*bl
  150 continue
      go to 210
c
c -- case 6    heme methyl protons
c
  160 continue
      call unit(1,2,r1,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #6.1600a, bond r1 has zero length. ',
     .              ip,na,nhyd,protp)
      end if
      call unit(3,2,r2,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #6.1600b, bond r2 has zero length. ',
     .              ip,na,nhyd,protp)
      end if
      call unit(4,2,r3,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #6.1600c, bond r3 has zero length. ',
     .              ip,na,nhyd,protp)
      end if
      call mult(r4,r2,r3)
      call mult(r6,r4,r1)
      do 170 ic=1,3
         protp(ic,1)=scoord(ic,1) + (r4(ic)*c6+r1(ic)*c5)*dh(6)
         protp(ic,2)=scoord(ic,1) -
     .                   ((r4(ic)*c4-r6(ic)*c3)*c6 - r1(ic)*c5)*dh(6)
         protp(ic,3)=scoord(ic,1) -
     .                   ((r4(ic)*c4+r6(ic)*c3)*c6 - r1(ic)*c5)*dh(6)
  170 continue
      go to 210
c
c -- case 7    ch2oh type hydroxyl proton
c
  180 continue
      call unit(1,2,r1,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #7.1800a, bond r1 has zero length. ',
     .              ip,na,nhyd,protp)
      end if
      call unit(3,2,r2,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #7.1800b, bond r2 has zero length. ',
     .              ip,na,nhyd,protp)
      end if
      call mult(r3,r2,r1)
      call mult(r4,r1,r3)
      do 190 ic=1,3
         protp(ic,1)=scoord(ic,1) + (r1(ic)*c5 - r4(ic)*c6)*dh(7)
  190 continue
      go to 210
c
c -- case 13    ch2sh type proton
c
  192 continue
      call unit(1,2,r1,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #7.1800a, bond r1 has zero length. ',
     .              ip,na,nhyd,protp)
      end if
      call unit(3,2,r2,ierror)
      if (ierror .gt. 0) then
            call exit_error('case #7.1800b, bond r2 has zero length. ',
     .              ip,na,nhyd,protp)
      end if
      call mult(r3,r2,r1)
      call mult(r4,r1,r3)
      do 194 ic=1,3
         protp(ic,1)=scoord(ic,1) + (r1(ic)*c5 - r4(ic)*c6)*dh(13)
  194 continue
      go to 210
c
  200 continue
      write(nerrout,*)'eqns. for h position have no unique soln.',
     .                ' all four atoms are co-planar'
  210 continue
      go to (220,230,220,230,240,240,220,220,230,240,230,220,220),ip
      call exit_error('error in computed goto #2 of subr hydros',
     1              ip,na,nhyd,protp)
  220 continue
      nhyd=1
      go to 250
  230 continue
      nhyd=2
      go to 250
  240 continue
      nhyd=3
  250 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine exit_error(str, l1, l2, l3, rprotp)
      implicit none
#include "protonate.h"
      character*40 str
      integer l1, l2, l3, j
      real  rprotp(3,4)
 1000 format('ip= ',i5,' na= ',i5,' nhyd= ',i5,' protp(1-3,1-4)= ',
     .          12f6.2)
 1010 format('error:  ',a40)
      write(nerrout,1010)str
      write(nerrout,1000)l1,l2,l3,(rprotp(1,j),j=1,4),
     .                            (rprotp(2,j),j=1,4),
     .                            (rprotp(3,j),j=1,4)
      call mexit(0,1)
      end
c----------------------------------------------------------------------
      subroutine unit(j1, j2, r1, ierror)
      implicit none
c
#include "protonate.h"
c
c -- find the unit vector along the bond joining two atoms
c
      integer j1, j2, ierror, ic, k
      real r1(3), ar1
      intrinsic sqrt
c
      ierror = 0
      do ic=1,3
         r1(ic)=scoord(ic,j1) - scoord(ic,j2)
      end do
      ar1=sqrt(r1(1)*r1(1)+r1(2)*r1(2)+r1(3)*r1(3))
      if (ar1 .lt. 1e-20) then
          write(nerrout,*) 'error in finding unit vector:'
          write(nerrout,'(2i5,6f8.3)') j1,j2,
     .                     (scoord(k,j1),k=1,3),(scoord(k,j2),k=1,3)
          ierror=1
          return
      end if
      do ic=1,3
         r1(ic)=r1(ic)/ar1
      end do
      return
      end
c----------------------------------------------------------------------
      subroutine mult(r3,r1,r2)
      implicit none
c
c -- find the unit vector perpendicular to the two given vectors
c
      real r1(3),r2(3),r3(3), ar3
      integer ic
      intrinsic sqrt
c
      r3(1)=r1(2)*r2(3) - r2(2)*r1(3)
      r3(2)=-(r1(1)*r2(3) - r1(3)*r2(1))
      r3(3)=r1(1)*r2(2) - r2(1)*r1(2)
      ar3=sqrt(r3(1)*r3(1)+r3(2)*r3(2)+r3(3)*r3(3))
      do ic=1,3
         r3(ic)=r3(ic)/ar3
      end do
      return
      end
c
c
      subroutine writeout( bstr, atomname, resname, id, resnum, resc,
     .                     x, y, z )
      implicit none
c
#include "protonate.h"
c
c -- centralize the output procedure   gpg
c
      character*6 bstr
      character*4 atomname, printa
      character*3 resname
      character*1 mark, id, resc
      character*40 line
      logical header, longbond
      integer resnum, i, imatch
      real x, y, z, TOL, ONE, ZERO, dis, dish, disn
      parameter (TOL=2.0, ONE=1.00, ZERO=0.00)
      intrinsic sqrt
c
      data header /.true./
c
      line = '----------------------------------------' 
      printa = atomname
c
c
c -- check to see if this proton overlaps with one from input:
c
      if ( atomname(1:1).eq.'H' ) then
c 
c -- process hydrogens:  check for matches with this (built)
c    proton and any that may have been in the input file:
c
         dis = 999.
         imatch = 0
c
c -- find the closest proton in the same residue in the input file:
c
         do i = 1, natom
           if ( isoldp(i) .and.
     .          .not. matched(i) .and. 
     .          ires(i)    .eq. resnum .and.
     .          crescod(i)  .eq. resc .and.
     .          chainid(i) .eq. id ) then
             disn = ( ( x - coord(1,i) ) * ( x - coord(1,i) ) ) +
     .              ( ( y - coord(2,i) ) * ( y - coord(2,i) ) ) +
     .              ( ( z - coord(3,i) ) * ( z - coord(3,i) ) )
             disn = sqrt(disn)
             if ( disn .lt. dis ) then
               dis = disn
               imatch = i
             end if
           end if
         end do
         if ( imatch .eq. 0 ) go to 3000
c
c -- check if closest proton is bonded to same atom as proton
c    currently being considered.  Bonded (heavy) atom is in
c    the first three elements of scoord.
c
         i = imatch
         dish = ( (scoord(1,1)-coord(1,i))*(scoord(1,1)-coord(1,i)) ) +
     .          ( (scoord(2,1)-coord(2,i))*(scoord(2,1)-coord(2,i)) ) +
     .          ( (scoord(3,1)-coord(3,i))*(scoord(3,1)-coord(3,i)) )
         dish = sqrt(dish)
         longbond = dish .gt. 1.2
         if( resname.eq.'CYS' .and. atomname.eq.'HG  ' .and.
     .         dish.lt.1.95 ) longbond = .false.
c
         if ( mismch .and. longbond ) then
           write(nerrout,1010)
           write(nerrout,1020) i,atomname,resname,cnamat(i),dis,dish
           write(nerrout,1030) x, y, z
           write(nerrout,1030) coord(1,i),  coord(2,i),  coord(3,i)
           write(nerrout,1030) scoord(1,1), scoord(2,1), scoord(3,1)
 1010      format( 'proton mismatch:' )
 1020      format( I4, 1X, A4, 1X, A3, 1X, A4, 1X, 2F9.2 )
 1030      format( 3F9.2 )
         end if
c
c        ---convert atom names PDB version 3:
c
c        ---First, assume that there are no two-character element names
c           (like Fe or Ca or Na).  Then, according to Brookhaven rules,
c           column 13 will be blank, and the name will be left-justified
c           starting in column 14.  UNLESS, the name is four characters
c           long!  In that case, don't use the first blank.
c
 3000       if (atomname(4:4) .eq. ' ') then 
               printa(1:1) = ' '
               printa(2:4) = atomname(1:3)
            else
               printa(1:4) = atomname(1:4)
            endif
c
c       --- Special fixes where (old) Amber nucleic acid atom names differ from
c           version 3 pdb names:
c
            if( printa(1:4) .eq. 'H5''1' ) printa(1:4) = ' H5'''
            if( printa(1:4) .eq. 'H5''2' ) printa(1:4) = 'H5'''''
            if( printa(1:4) .eq. 'H2''1' ) printa(1:4) = ' H2'''
            if( printa(1:4) .eq. 'H2''2' ) printa(1:4) = 'H2'''''
            if( printa(1:4) .eq. ' H5T' ) printa(1:4) = 'HO5'''
            if( printa(1:4) .eq. ' H3T' ) printa(1:4) = 'HO3'''
            if( printa(1:4) .eq. 'HO''2' ) printa(1:4) = 'HO2'''
c
c       --- Now, special case out the two-character element names:
c
            if( printa(1:4).eq.' Na+' .or. printa(1:4).eq.' NA+' .or.
     +          printa(1:3).eq.' Fe' .or. printa(1:3).eq.' FE' .or.
     +          printa(1:3).eq.' Cl' .or. printa(1:3).eq.' CL' .or.
     +          printa(1:3).eq.' Zn' .or. printa(1:3).eq.' ZN' .or.
     +          printa(1:4).eq.' Li+' .or. printa(1:4).eq.' LI+' .or.
     +          printa(1:4).eq.' Ca+' .or. printa(1:4).eq.' CA+' .or.
     +          printa(1:4).eq.' Mg+' .or. printa(1:4).eq.' MG+' .or.
     +          printa(1:4).eq.' Br-' .or. printa(1:4).eq.' BR-' ) then
              printa(1:1) = printa(2:2)
              printa(2:2) = printa(3:3)
              printa(3:3) = printa(4:4)
              printa(4:4) = ' '
            end if
c
         if ( imatch .eq. 0 ) go to 5000
c
c -- check for matches, but add special code so that atom "H" in
c    residue 1 can never be matched:
c
         if ( (dis .lt. TOL .or.
     .         (resname.eq.'CYS' .and. atomname.eq.'HG  ' .and.
     .          dis.lt.2.5) )
     .        .and. .not.longbond 
     .        .and. ( cnamat(i) .ne. ' H  ' .or.
     .                ires(i) .ne. 1 ) ) then
c
c -- a "matching" proton has been found; mark it as matched;
c    substitute its coords for the model built coords if
c    -k has been requested.
c
           matched(i) = .true.
           if ( keep ) then
             x = coord(1,i)
             y = coord(2,i)
             z = coord(3,i)
           end if
           if ( lprint ) then
              mark = ' '
              if ( dis .gt. 0.5 ) mark = '*'
              if ( printa .ne. cnamat(i) ) then
                if ( header ) then
                  write(nerrout,2010)
                  write(nerrout,2020) line
                  write(nerrout,2030)
                  write(nerrout,2020) line
                  header = .false.
                end if
                write(nerrout,2040) printa, resname, resnum, x, y, z,
     .                              cnamat(i), cresnam(i), ires(i),
     .                              coord(1,i), coord(2,i), coord(3,i),
     .                              dis, mark
              end if
           end if
         end if
 2010    format( '                       Proton substitutions' )
 2020    format( A )
 2030    format( '         in output file             ',
     .           '         in input file' )
 2040    format( A4, 1X, A3, I4, 3F7.2, ' | ', A4, 1X, A3, I4, 3F7.2,
     .           F9.2, A1 )
c
      end if
c
c -- Write a TER record between chains
c
      if ( id .ne. chainp ) then
        write( nfilout, 3010 )
 3010   format( 'TER' )
      end if
c
c -- End of residue ???
c
      if ( resnum .ne. iresp .or.
     .     id     .ne. chainp .or.
     .     resc   .ne. cinsert ) then
         iresp = resnum
         chainp = id
         cinsert = resc
         iresct = iresct + 1
         if ( ndbgout .ne. 0 ) then
            write( ndbgout, 3020 ) iresct, iresp, icount - ilast
 3020       format( 3I7 )
            ilast = icount
         end if
      end if
 5000 continue
      if ( lrenum ) then
         write( nfilout, 6000 ) bstr, icount, printa, resname, id,
     .                          iresct+1, ' ', x, y, z
      else
         write( nfilout, 6000 ) bstr, icount, printa, resname, id,
     .                          resnum, resc, x, y, z
      endif
      icount = icount + 1
c
      return
 6000 format( A6, I5, 1X, A4, 1X, A3, 1X, A1, I4, A1, 3X, 3F8.3)
      end
c
      subroutine maxerr( iat )
      implicit none
      integer iat
#include "protonate.h"
c
      write(nerrout,*) iat, ' atoms read from file'
      write(nerrout,*) 'This version only dimensioned for ',
     .                 MAXAT,' atoms!'
      call mexit(0,1)
      return
      end
