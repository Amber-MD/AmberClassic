      program     GuesWatH
c
c**********************************************************************c
c                                                                      c
c    Program : GuesWatH - create hydrogen positions of water           c
c    Author : Shuichi Miyamoto                                         c
c    Modified: D.A. Case                                               c
c    Date of last update : February 8, 1994                            c
c                                                                      c
c**********************************************************************c
c
      parameter   (maxsph=15000,maxcir=1000)
      parameter   (maxatm=100000,maxwat=3000,maxwa3=9000)
      character*80 title,prmtop,watpdb,arg
      character*4  atnam(maxatm),atnamw(maxatm),rsnamt(maxatm) 
      character*3  rsnam(maxatm),rsnamw(maxatm)
      character*1  chid(maxatm) 
      integer*4   resno(maxatm) 
      dimension   msph(3,maxsph), ncirc(maxsph),icirc(maxcir,maxsph)
      dimension   lp(3,2), lpp(2,maxcir,maxsph)
      dimension   cord(3,maxatm), kres(maxatm)
     &          , chrg(maxatm), elecp(maxsph,maxwat)
     &                        , elecw(maxsph,maxwat)
      dimension   occ(maxwat),temp(maxwat)
      dimension   cordo(3,maxwat),cordh(3,2,maxwat),cordw(3,maxwa3)
      dimension   dmin(maxwat),emin(maxwat),isort(maxwat)
      dimension   iterm(maxwa3)
      dimension   rh(2), hnew(3,2)
      dimension   vech1(3),vech2(3),vecx(3), unset(3,3)
c
      data  oh /0.9572/ , ohcos /0.586077/, hh2 /0.7568/
     &     ,chrg_ow /-15.1974/ , chrg_hw /7.59870/
c
c
      prmtop = 'prmtop'
      watpdb = 'watpdb'
      iddd = 0
      iord = 0
      nargs = iargc()
      iarg = 0
      if (nargs.eq.iarg) go to 20
   10 iarg = iarg + 1
          call getarg(iarg,arg)
          if (arg.eq.'-p') then
            iarg = iarg + 1
            call getarg(iarg,prmtop)
          else if (arg.eq.'-w') then
            iarg = iarg + 1
            call getarg(iarg,watpdb)
          else if (arg.eq.'-consd') then
            iddd = 1
          else if (arg.eq.'-elstat') then
            iord = 1
          else if (arg.eq.'-h' .or. arg.eq.'--help') then
            write(0,*) 'Usage: gwh [-p prmtop ] [-w watpdb ]',
     .       ' [-consd] [-elstat] < input-pdb-file > output-pdb-file'
      write(0,*) '  watpdb contains new waters to have hydrogens fixed'
      write(0,*) '  -consd: use a constant dielectric'
      write(0,*) '  -elstat: sort on electrostatics; otherwise sort on'
      write(0,*) '           minimum H-bond distance'
      write(0,*) 'Note: waters in input-pdb-file will not be modified'
            call mexit(0,1)
          else
            write(0,*) 'Bad argument: ',iarg,arg
            write(0,*) 'Usage: gwh [-p prmtop ] [-w watpdb ]',
     .       ' [-consd] [-elstat] < input-pdb-file > output-pdb-file'
      write(0,*) '  watpdb contains new waters to have hydrogens fixed'
      write(0,*) '  -consd: use a constant dielectric'
      write(0,*) '  -elstat: sort on electrostatics; otherwise sort on'
      write(0,*) '           minimum H-bond distance'
      write(0,*) 'Note: waters in input-pdb-file will not be modified'
            call mexit(0,1)
          end if
      if (iarg.lt.nargs) go to 10
   20 continue
      call amopen(1,prmtop,'O', 'F', 'R')
      call top_3 (1, title,natom,atnam,chrg,nres,rsnamt,kres,ifbox)
      close(1)
      call pdbatm (maxatm,5,0,0, natomp,atnam,
     .      rsnam,chid,resno,cord,iterm)
      call amopen(3,watpdb,'O','F','R')
      call pdbhet2 (3,0,0, no,atnamw,rsnamw,cordo,occ,temp)
      close(3)
c
      if (natom.ne.natomp) then
         write  (0,'(a,2i6)') 'diff in natom: ',natom,natomp
         call mexit(0,1)
      endif
c
      mcycle =5
c
c
c                                             ----- Set up of Sphere -----
      weig = 1.0
      nsph = 0
      do i=1,21
         do j=1,21
            do k=1,21
               ijk = (i-11)*(i-11) + (j-11)*(j-11) + (k-11)*(k-11) 
               if (ijk.gt.81 .and. ijk.le.100) then
                  nsph = nsph + 1
                  msph(1,nsph) = i - 11
                  msph(2,nsph) = j - 11
                  msph(3,nsph) = k - 11
               endif
            end do
         end do
      end do
c
      do i=1,nsph-1
         nhit = 0
         do j=i+1,nsph
            kx = (msph(1,i) - msph(1,j)) ** 2
            ky = (msph(2,i) - msph(2,j)) ** 2
            kz = (msph(3,i) - msph(3,j)) ** 2
            kxyz = kx + ky + kz
            if (kxyz.gt.213 .and. kxyz.le.243) then
               nhit = nhit + 1
               icirc(nhit,i) = j
               vech1(1) = float(msph(1,i)) / 10.
               vech1(2) = float(msph(2,i)) / 10.
               vech1(3) = float(msph(3,i)) / 10.
               vech2(1) = float(msph(1,j)) / 10.
               vech2(2) = float(msph(2,j)) / 10.
               vech2(3) = float(msph(3,j)) / 10.
               vecx(1) = vech1(1) + vech2(1)
               vecx(2) = vech1(2) + vech2(2)
               vecx(3) = vech1(3) + vech2(3)
               call uniset (1,vecx,vech2, unset)
               lp(1,1) = - unset(1,1) * 5.80 + unset(1,3) * 8.15
               lp(2,1) = - unset(2,1) * 5.80 + unset(2,3) * 8.15
               lp(3,1) = - unset(3,1) * 5.80 + unset(3,3) * 8.15
               lp(1,2) = - unset(1,1) * 5.80 - unset(1,3) * 8.15
               lp(2,2) = - unset(2,1) * 5.80 - unset(2,3) * 8.15
               lp(3,2) = - unset(3,1) * 5.80 - unset(3,3) * 8.15 
               do l=1,2
                  kmin = 99999
                  do k=1,nsph
                     kx = (lp(1,l) - msph(1,k)) ** 2
                     ky = (lp(2,l) - msph(2,k)) ** 2
                     kz = (lp(3,l) - msph(3,k)) ** 2  
                     if (kmin .gt. kx + ky + kz) then
                         kmin  =   kx + ky + kz
                        lpp(l,nhit,i) = k
                     endif
                  end do
               end do
            endif
         end do
         ncirc(i) = nhit
      end do
c
c                                               --- Calc Dist from Prot ---
      do io=1,no
         dismin = 999.e10
         do k=1,natom
            dist2 = 0.
            do j=1,3
               dist2 = dist2 + (cordo(j,io) - cord(j,k)) ** 2
            end do
            if (dismin.gt.dist2) dismin = dist2
         end do
         dmin(io) = dismin
      end do
c
c                                               --- Calc Elec from Prot ---
c
      do io=1,no
         enemin = 999.e10
         do i=1,nsph
            xs = float(msph(1,i)) / 10. + cordo(1,io)
            ys = float(msph(2,i)) / 10. + cordo(2,io)
            zs = float(msph(3,i)) / 10. + cordo(3,io)
            esum = 0.
            do k=1,natom
               r2 = ( ( xs - cord(1,k) ) * ( xs - cord(1,k) )
     &              + ( ys - cord(2,k) ) * ( ys - cord(2,k) )
     &              + ( zs - cord(3,k) ) * ( zs - cord(3,k) ) )
               if (iddd.eq.0) then
                  esum = esum + chrg(k) / r2
               else
                  esum = esum + chrg(k) / sqrt(r2)
               endif
            end do
            elecp(i,io) = esum
            if (enemin .gt. elecp(i,io) )  enemin = elecp(i,io)
         end do
         emin(io) = enemin
      end do
c
      if (iord.eq.0) then
         call sortf_as (no,dmin, isort)
      else
         call sortf_as (no,emin, isort)
      endif
c
c                                                 ----- Water Loop -----
      do l=1,mcycle
c
      etot = 0.
      rhh2 = 0.
c
      do io=1,no
         iso = isort(io)
         if (io.eq.1) go to 170
c                                                 --- Calc ElecW ---
         do i=1,nsph
            xs = float(msph(1,i)) / 10. + cordo(1,iso)
            ys = float(msph(2,i)) / 10. + cordo(2,iso)
            zs = float(msph(3,i)) / 10. + cordo(3,iso)
            esum = 0.
            do iio=1,no
               if (l.eq.1 .and. iio.eq.io) exit
               if (iio.eq.io) cycle
               iiso = isort(iio)
               if (k.eq.io) cycle
               r2 = ( xs - cordo(1,iiso) ) ** 2
     &            + ( ys - cordo(2,iiso) ) ** 2
     &            + ( zs - cordo(3,iiso) ) ** 2
               if (iddd.eq.0) then
                  esum = esum + chrg_ow / r2
               else
                  esum = esum + chrg_ow / sqrt(r2)
               endif
               do m=1,2
                  rh(m) = ( xs - cordh(1,m,iiso) ) ** 2
     &               + ( ys - cordh(2,m,iiso) ) ** 2
     &               + ( zs - cordh(3,m,iiso) ) ** 2
                  if (iddd.eq.0) then
                     esum = esum + chrg_hw / rh(m)
                  else
                     esum = esum + chrg_hw / sqrt(rh(m))
                  endif
               end do
            end do
c
            elecw(i,iso) = esum
         end do
c
c
  170    continue
c                                                 --- Calc Min ---
         emint = 999.e10
         do i=1,nsph
            if (ncirc(i).eq.0) cycle
            eminj = 999.e10
            do j=1,ncirc(i)
               elp = elecp(lpp(1,j,i),iso) + elecw(lpp(1,j,i),iso)
     &             + elecp(lpp(2,j,i),iso) + elecw(lpp(2,j,i),iso)
               if (eminj.gt.elecp(icirc(j,i),iso)
     &                    + elecw(icirc(j,i),iso) - elp*weig) then
                   eminj =  elecp(icirc(j,i),iso)
     &                    + elecw(icirc(j,i),iso) - elp*weig
                   nminj = icirc(j,i)
               endif
            end do
            if (emint .gt. elecp(i,iso) + elecw(i,iso) + eminj) then
                emint  =   elecp(i,iso) + elecw(i,iso) + eminj
                nmini = i
                nmink = nminj
            endif
         end do
c
         etot = etot + emint
c
         vech1(1) = float(msph(1,nmini)) / 10.
         vech1(2) = float(msph(2,nmini)) / 10.
         vech1(3) = float(msph(3,nmini)) / 10.
         vech2(1) = float(msph(1,nmink)) / 10.
         vech2(2) = float(msph(2,nmink)) / 10.
         vech2(3) = float(msph(3,nmink)) / 10.
         vecx(1) = vech1(1) + vech2(1)
         vecx(2) = vech1(2) + vech2(2)
         vecx(3) = vech1(3) + vech2(3)
         call uniset (1,vecx,vech2, unset)
         do m=1,3
            hnew(m,1) = cordo(m,iso)
     &                + ohcos * unset(m,1) + hh2 * unset(m,2)
            hnew(m,2) = cordo(m,iso)
     &                + ohcos * unset(m,1) - hh2 * unset(m,2)
            rhh2 = rhh2 + ( hnew(m,1) - cordh(m,1,iso) ) ** 2
     &                  + ( hnew(m,2) - cordh(m,2,iso) ) ** 2
            cordh(m,1,iso) = hnew(m,1)
            cordh(m,2,iso) = hnew(m,2)
         end do
      end do
c
c
      avet = etot / float(no)
      rms = sqrt ( 0.5 * rhh2 / float(no) )
c
      write  (0,220) avet,rms
  220 format ('REMARK   Av.E =',f8.3,' rms =',f8.3)
c
      end do
c
c                                                 ----- End of Water Loop -----
c
      iok = 0
      do 260 io=1,no
         iok = iok + 1
         iio  = iok * 3 - 2
         iih1 = iio + 1
         iih2 = iio + 2 
         do 240 m=1,3
            cordw(m,iio) = cordo(m,io)
            cordw(m,iih1) = cordh(m,1,io)
            cordw(m,iih2) = cordh(m,2,io)
  240    continue
         atnamw(iio)  = ' O  '
         atnamw(iih1) = ' H1 '
         atnamw(iih2) = ' H2 '
         do 250 m=1,3
            rsnamw(iio-1+m) = 'HOH'
            chid  (iio-1+m) = ' '
            resno (iio-1+m) = io + nres
  250    continue
         iterm(iio)  = 0
         iterm(iih1) = 0
         iterm(iih2) = 1
  260 continue
c
      natw = 3 * iok
      call pdbatm(maxatm,0,0,6,
     .      natw,atnamw,rsnamw,chid,resno,cordw,iterm)
c
      call mexit(6,0)
      end
