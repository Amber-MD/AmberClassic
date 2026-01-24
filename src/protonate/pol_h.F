      program  Pol_H
      implicit none
c
c*************************************************************************
c  Program : Pol_H - reset positions of polar hydrogens
c  Author : Shuichi Miyamoto
c  Modified: D.A. Case
c  Modified: J.P.Holland 18/7/96 -- explicit declaration of all variables
c*************************************************************************
c
      integer      mcirc, maxatm, maxpol
      parameter    ( mcirc = 180, maxatm=10000, maxpol=200 )
      character*1024 infile, outfile
      character*80 title, arg
      character*4  atnam(maxatm), atnamo(maxatm), rsnamt(maxatm)
      character*3  rsnam(maxatm), rsnamo(maxatm)
      character*1  chid(maxatm)
      integer      resno(maxatm), nprmin, nfilin, nfilout
      logical      pdbq
      real         cord(3,maxatm), chrg(maxatm),
     .             chrgm(maxatm), crdcir(3,mcirc,maxpol), elecir(mcirc),
     .             torcir(mcirc,maxpol), emin(maxpol),
     .             cord1(3), cord2(3), cord3(3)
      integer      iterm(maxatm), kres(maxatm), ipatm(maxpol),
     .             lys(maxpol), i13(maxpol), kmin(maxpol), isort(maxpol)
      real         oh, chrg_ow, chrg_hw, dtor, pi, rad, rad2, weitwo,
     .             bl, ba, enesum, esum, r2, enemin
      integer      nargs, iarg, natom, nres, ifbox, mcycle,
     .             no, i, iargc, ifound, npol, j, k, l,
     .             icycle, is, iatm, ii, iis, jatm, iiatm, lend, kirmin,
     .             mcirc3, iqr
c     external     iargc, getargs
c
      data  oh /0.9572/ , chrg_ow /-15.1974/ , chrg_hw /7.59870/
      data  dtor /0.0174532/ , pi /3.14159/
      data  rad /3.2/
c
      rad2    = rad * rad
      infile  = '-'
      outfile = '-'
      nprmin  = 1
      nfilin  = 5
      nfilout = 6
      nargs   = iargc()
#ifdef HITACHI
      iarg = 1
#else
      iarg = 0
#endif
      if (nargs.eq.iarg) go to 20
   10 continue
          iarg = iarg + 1
          call getarg(iarg,arg)
          if ( arg .eq. '-i' ) then
            iarg = iarg + 1
            call getarg(iarg,infile)
          else if ( arg .eq. '-o' ) then
            iarg = iarg + 1
            call getarg(iarg,outfile)
          else
            call badarg( iarg, arg )
          end if
      if ( iarg .lt. nargs ) go to 10
c
   20 continue
      if ( infile .ne. '-' ) then
         nfilin = 31
         call amopen( nfilin, infile,   'O', 'F', 'R' )
      end if
      if ( outfile .ne. '-' ) then
         nfilout = 32
         call amopen( nfilout, outfile, 'U', 'F', 'W' )
      end if
c
      call pdbatm( maxatm, nfilin, 0, 0,
     .      natom, atnam, rsnam, chid, resno, cord, iterm, chrg )
c
c ---"5" seems like a good default value here:
c
      mcycle = 5
c
c   ----- Charge Setup -----
c
      do i=1,natom
        chrgm(i) = chrg(i)
      end do
c
c     ---cycle through atoms: set chrgm to zero for polar O's, N's, H's;
c
      do i=1,natom
c
         if ((rsnam(i).eq.'SER' .or. rsnam(i).eq.'THR') .and.
     &       (atnam(i)(2:3).eq.'OG' .or. atnam(i)(2:3) .eq. 'HG') 
     &   .or. (rsnam(i).eq.'TYR' .and. atnam(i)(2:3).eq.'OH') 
     &   .or. (rsnam(i).eq.'LYS' .and. atnam(i)(2:3).eq.'NZ') 
     &   .or. (rsnam(i).eq.'TYR' .and. atnam(i)(2:3).eq.'HH') 
     &   .or. (rsnam(i).eq.'LYS' .and. atnam(i)(2:3).eq.'HZ') )
     &        chrgm(i) = 0.
c
      end do
c
c    ----- Coord Setup -----
c
      npol = 0
c
c     --- again cycle through atoms
c           create circles of possible H positions for OH groups
c           (seems to have some hard-wired atom positions here)
c
      do i=1,natom
c
         if (rsnam(i).eq.'SER' .and. atnam(i)(2:3) .eq. 'OG') then
            npol = npol + 1
            ipatm(npol) = i + 1
            lys(npol) = 0
            i13(npol) = i - 3
            do j=1,3
               cord1(j) = cord(j,i-5)
               cord2(j) = cord(j,i-3)
               cord3(j) = cord(j,i)
            end do
            bl = 0.96
            ba = 108.5
            call circl2 (cord1,cord2,cord3,bl,ba,npol,maxpol, crdcir)
            do k=1,180
               torcir(k,npol) = 0.5 * ( 1. + cos (6.* float(k) * dtor) )
            end do
         endif
c
         if (rsnam(i).eq.'THR' .and. atnam(i)(2:3) .eq. 'OG') then
            npol = npol + 1
            ipatm(npol) = i + 1
            lys(npol) = 0
            i13(npol) = i - 6
            do j=1,3
               cord1(j) = cord(j,i-8)
               cord2(j) = cord(j,i-6)
               cord3(j) = cord(j,i)
            end do
            bl = 0.96
            ba = 108.5
            call circl2 (cord1,cord2,cord3,bl,ba,npol,maxpol, crdcir)
            do k=1,180
               torcir(k,npol) = 0.5 * ( 1. + cos (6.* float(k) * dtor) )
            end do
         endif
c
         if (rsnam(i).eq.'TYR' .and. atnam(i)(2:3) .eq. 'OH') then
            npol = npol + 1
            ipatm(npol) = i + 1
            lys(npol) = 0
            i13(npol) = i - 2
            do j=1,3
               cord1(j) = cord(j,i-3)
               cord2(j) = cord(j,i-1)
               cord3(j) = cord(j,i)
            end do
            bl = 0.96
            ba = 113.0
            call circl2 (cord1,cord2,cord3,bl,ba,npol,maxpol, crdcir)
            do k=1,180
               torcir(k,npol) = 1.8 * ( 1. +
     .                                   cos (4.* float(k)*dtor + pi) )
            end do
         endif
c
         if (rsnam(i).eq.'LYS' .and. atnam(i)(1:3) .eq. ' NZ') then
            npol = npol + 1
            ipatm(npol) = i + 1
            lys(npol) = 1
            i13(npol) = i - 3
            do j=1,3
               cord1(j) = cord(j,i-6)
               cord2(j) = cord(j,i-3)
               cord3(j) = cord(j,i)
            end do
            bl = 1.01
            ba = 109.5
            call circl2 (cord1,cord2,cord3,bl,ba,npol,maxpol, crdcir)
            do k=1,180
               torcir(k,npol) = 1.4 * ( 1. + cos (6.* float(k) * dtor) )
            end do
         endif
      end do
c
      do i=1,npol
         isort(i) = i
      end do
c
c ----- Main Cycle -----
c
      icycle = 1
  160 continue
      enesum = 0.
c
c --- Calc elec ---
c
      do i=1,npol
         is = isort(i)
         iatm = ipatm(is)
c
         do k=1,mcirc
            esum = 0.0
            do 170 ii=1,natom
               if (ii.eq.iatm .or. ii.eq.i13(is)) go to 170
               if ((ii.eq.iatm+1 .or. ii.eq.iatm+2) .and. lys(is).eq.1)
     &            go to 170
               r2 = ( crdcir(1,k,is) - cord(1,ii) ) ** 2
     &            + ( crdcir(2,k,is) - cord(2,ii) ) ** 2
     &            + ( crdcir(3,k,is) - cord(3,ii) ) ** 2
               esum = esum + chrg(iatm) * chrgm(ii) / sqrt(r2)
  170       continue
c
            if (icycle.ne.1 .and. i.ne.1 ) then
c
               do ii=1,i-1
                  iis = isort(ii)
                  iiatm = ipatm(iis)
                  lend = lys(iis) * 2
                  do l=-1,lend
                     jatm = iiatm + lend
                     r2 = ( crdcir(1,k,is) - cord(1,jatm) ) ** 2
     &                  + ( crdcir(2,k,is) - cord(2,jatm) ) ** 2
     &                  + ( crdcir(3,k,is) - cord(3,jatm) ) ** 2
                     esum = esum + chrg(iatm) * chrg(jatm) / sqrt(r2)
                  end do
               end do
c
            end if
            elecir(k) = esum
         end do
c
c --- Calc min ---
c
         if (lys(is).eq.0) then
            enemin = 999.e10
            do k=1,mcirc
               if (enemin .gt. elecir(k) + torcir(k,is) ) then
                   enemin = elecir(k) + torcir(k,is)
                   kirmin = k
               endif
            end do
         else
            enemin = 999.e10
            mcirc3 = mcirc / 3
            do k=1,mcirc3
               if (enemin .gt. elecir(k) + elecir(k+60) + elecir(k+120)
     &                                   + torcir(k,is) )  then
                   enemin = elecir(k) + elecir(k+60) + elecir(k+120)
     &                    + torcir(k,is)
                   kirmin = k
               endif
            end do
         endif
         emin(is) = enemin
         kmin(is) = kirmin
         enesum = enesum + emin(is)
c
c --- Put coord ---
c
         if (lys(is).eq.0) then
            cord(1,iatm) = crdcir(1,kmin(is),is)
            cord(2,iatm) = crdcir(2,kmin(is),is)
            cord(3,iatm) = crdcir(3,kmin(is),is)
         else
            do j=1,3
               cord(j,iatm)   = crdcir(j,kmin(is),is)
               cord(j,iatm+1) = crdcir(j,kmin(is)+60,is)
               cord(j,iatm+2) = crdcir(j,kmin(is)+120,is)  
            end do
         endif
c
      end do
c
c --- Sort emin ---
c
      call sortf_as (npol,emin, isort)
#ifdef DEBUG
c
         write  (9,270) icycle,enesum,(isort(i),i=1,npol)
  270    format (/' after',i2,' th cycle'/f8.2/10(9x,10i7/))
         write  (9,280) (emin(i),i=1,npol)
         write  (9,290) (kmin(i),i=1,npol)
  280    format (9x,10f7.2)
  290    format (9x,10i7)
c
#endif
      icycle = icycle + 1
      if (icycle.le.mcycle) go to 160
c
c ----- End of Main Cycle -----
c
c ----- Write Coords -----
c
      do i=1,npol
         iatm = ipatm(i)
         if (lys(i).eq.0) then
            cord(1,iatm) = crdcir(1,kmin(i),i)
            cord(2,iatm) = crdcir(2,kmin(i),i)
            cord(3,iatm) = crdcir(3,kmin(i),i)
         else
            do j=1,3
               cord(j,iatm)   = crdcir(j,kmin(i),i)
               cord(j,iatm+1) = crdcir(j,kmin(i)+60,i)
               cord(j,iatm+2) = crdcir(j,kmin(i)+120,i)  
            end do
         endif
      end do
c
      call pdbatm( maxatm, 0, 0, nfilout,
     .        natom, atnam, rsnam, chid, resno, cord, iterm, chrg )
      call mexit( 0, 0 )
      end
c
c======================================================================c
c
      subroutine circl2 (cord1,cord2,cord3,bl,ba,ipatm,maxpol,crdcir)
      implicit none
c
      integer ipatm, maxpol, i
      real cord1(3), cord2(3), cord3(3), crdcir(3,180,maxpol),
     .     tcord(3), bl, ba, ta
c
      do i=1,180
         ta = float(i) * 2.0
         call cortor(cord1,cord2,cord3,bl,ba,ta,tcord)
         crdcir(1,i,ipatm) = tcord(1)
         crdcir(2,i,ipatm) = tcord(2)
         crdcir(3,i,ipatm) = tcord(3)
      end do
      return
      end
c
c======================================================================c
c
      subroutine badarg( iarg, arg )
c
      integer iarg
      character*(*) arg
c
      write(0,*) 'Bad argument: ',iarg,arg
      write(0,*) 'Usage: pol_h ',
     .           '[-i input-pqr-file] ',
     .           '[-o output-pdb-file]'
      call mexit(0,1)
      end
