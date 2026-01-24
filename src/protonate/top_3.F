c======================================================================c
c
      subroutine top_3 (intop, title,natom,atnam,chrg,nres,resnam,kres
     &                        ,ifbox)
      character*80 title
      character*20 line
      character*4  atnam(*),resnam(*)
      dimension    kres(*),chrg(*)
c
      read   (intop,*)
      read   (intop,*)
      read   (intop,*)
      read   (intop,'(a80)') title
c
      read   (intop,*)
      read   (intop,*)
      read   (intop,'(2i8)') natom,ntypes
      read   (intop,'(8x,i8)') nres
      read   (intop,'(56x,i8)') ifbox
      read   (intop,*)
c
      read   (intop,*)
      read   (intop,*)
      read   (intop,'(20a4)') (atnam(i),i=1,natom)
c
      read   (intop,*)
      read   (intop,*)
      read   (intop,'(5e16.8)') (chrg(i),i=1,natom)
c
      do i=1,99999
         read   (intop,'(a20)') line
         if( line(1:19) .eq. '%FLAG RESIDUE_LABEL' ) exit
      end do
      read   (intop,'(a20)') line
      read   (intop,'(20a4)') (resnam(i),i=1,nres)
c
      read   (intop,*)
      read   (intop,*)
      read   (intop,'(10i8)') (kres(i),i=1,nres)
c
c
c      read (NF,9108) ITITL
c      read (NF,9118) NTOTAT,NTYPES,NBONH,NBONA,NTHETH,NTHETS,NPHIH,
c     +          NPHIA,JHPARM,JPARM,NEXT,NTOTRS,MBONA,MTHETS,MPHIA,
c     +          MUMBND,MUMANG,MPTRA,NATYP,NHB,IFPERT,NBPER,NGPER,
c     +          NDPER,MBPER,MGPER,MDPER,IFBOX,NMXRS,IFCAP
c      NTYPE = NTYPES*NTYPES
c         WRITE(NF,9108) (IGRAPH(I),I = 1,NTOTAT)
c         WRITE(NF,9128) (CHRG(I),  I = 1,NTOTAT)
c         WRITE(NF,9128) (AMASS(I), I = 1,NTOTAT)
c         WRITE(NF,9118) (IAC(I),   I = 1,NTOTAT)
c         WRITE(NF,9118) (NUMEX(I), I = 1,NTOTAT)
c         WRITE(NF,9118) (NNO(I),   I = 1,NTYPE)
c         WRITE(NF,9108) (LABRES(I),I = 1,NTOTRS)
c         WRITE(NF,9118) (IPRES(I), I = 1,NTOTRS)
c         WRITE(NF,9128) (RK(I),    I = 1,MUMBND)
c         WRITE(NF,9128) (REQ(I),   I = 1,MUMBND)
c         WRITE(NF,9128) (TK(I),    I = 1,MUMANG)
c         WRITE(NF,9128) (TEQ(I),   I = 1,MUMANG)
c         WRITE(NF,9128) (PK(I),    I = 1,MPTRA)
c         WRITE(NF,9128) (PN(I),    I = 1,MPTRA)
c         WRITE(NF,9128) (PHASE(I), I = 1,MPTRA)
c         WRITE(NF,9128) (SOLTY(I), I = 1,NATYP)
c         NTTYP = (NTYPES*(NTYPES+1))/2
c         WRITE(NF,9128) (CN1(I),   I = 1,NTTYP)
c         WRITE(NF,9128) (CN2(I),   I = 1,NTTYP)
c 9108 FORMAT(20A4)
c 9118 FORMAT(12I6)
c 9128 FORMAT(5E16.8)
      return
      end
