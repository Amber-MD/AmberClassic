        program resp
C
C       RESP   version 2.1     October 1994 Jim Caldwell
C       RESP   version 2.0     September 1992
C          Author: Christopher Bayly
C
C       ESPFIT version 1.0 modified by Ian Gould to run in conjunction
C              with gaussian 90.
C
C       ESPFIT version 1.0 (G80UCSF):
C
C                 U.CHANDRA SINGH AND P.A.KOLLMAN
C
C       All authors:
C
C                 DEPARTMENT OF PHARMACEUTICAL CHEMISTRY
C                 SCHOOL OF PHARMACY
C                 UNIVERSITY OF CALIFORNIA
C                 SAN FRANCISCO   CA 94143
C
C
C---------------------------------------------------------------------
C
C     THIS PROGRAM FITS THE QUANTUM MECHANICALLY CALCULATED
C     POTENTIAL AT MOLECULAR SURFACES USING AN ATOM-CENTERED
C     POINT CHARGE MODEL. THE MOLECULAR SURFACES ARE GENERATED
C     BEYOND VANDER WAAL SURFACE IN ORDER TO MINIMISE OTHER
C     CONTRIBUTIONS SUCH AS EXCHANGE REPULSION AND CHARGE TRANSFER
C
C---------------------------------------------------------------------
C
C     -1st-   TITLE       FORMAT(10A8)
C
C---------------------------------------------------------------------
C
c     -2nd-
C
C           OPTIONS FOR THE JOB begin with " &cntrl"
c                               end with   " &end"
c           note leading blanks !!!!!!!!!!
C
C
C        INOPT   =  0  ... NORMAL RUN
C                =  1  ... CYCLE THROUGH A LIST OF DIFFERENT qwt
c                                read from -w unit 
c
C       IOUTOPT  =  0   NORMAL RUN
C                =  1   write restart info of new esp etc to 
c                              unit -e (esout unit)  
C
C         IQOPT  =  0  ... use the q's which are read the -i unit 
C                =  1  ... RESET ALL INITIAL CHARGES TO ZERO
C                =  2  ... READ IN NEW INITIAL CHARGES FROM -q (qwt) 
C                =  3  ... READ IN NEW INITIAL CHARGES FROM  -q (qwt)
C                                  AND PERFORM AVERAGING OF THOSE NEW 
C                                  INITIAL CHARGES ACCORDING TO IVARY VALUES
C
C         ihfree =  0  ... ALL ATOMS ARE RESTRAINED
C                =  1  ... HYDROGENS NOT RESTRAINED
C
C      irstrnt   =  0  ... HARMONIC RESTRAINTS (old style)
C                =  1  ... HYPERBOLIC RESTRAINT TO CHARGE OF ZERO (default)
C                =  2  ... ONLY ANALYSIS OF INPUT CHARGES; NO
C                          CHARGE FITTING IS CARRIED OUT
c
c      iunits    =  0  ... atom coordinates in angstroms
c                =  1       "     "          "  bohrs
c 
c 
c          qwt   =  restraint weight if irstrnt = 1
c
c NOTE: ESP coordinates must always be in Bohrs
C
C--------------------------------------------------------------------------
C
C     -3rd- wtmol .... relative weight for the molecule if 
c                    multiple molecule fit (1.0 otherwise) 
C
C             FORMAT(F10.5)
C
C--------------------------------------------------------------------------
c
c     -4th- subtitle for molecule
c
C------------------------------------------------------
C--------------------------------------------------------------------------
C
c     -5th- CHARGE,IUNIQ ( THE NUMBER OF UNIQUE CENTERS for this molecule)
c
C           FORMAT(2I5)
C
C--------------------------------------------------------------------------
C     -6th-  ONE CARD FOR EACH UNIQUE CENTER
C
C     FORMAT(I5,i5)
C
C      Name, IVARY
C
C             NAME = ATOMIC number  
C
C             IVARY = CONTROL OF CHARGE VARIATION OF THIS CENTER
C                   =  0 CHARGE VARIED INDEPENDENTLY OF PREVIOUS CENTERS
C                   = -n CHARGE FROZEN AT "INITIAL CHARGE" VALUE
C                   =  n CHARGE FITTED TOGETHER WITH CENTER n
C
C-------------------------------------------------------------------------
C
C     -7th-  intra molecule charge constraints...  blank line if no constr
C
C             FORMAT(I5,F10.5)
C
C     ngrp = number of centers in the group associated with this
C            constraint (i.e. the number of centers to be read in)
C
C    grpchg(i) = charge to which the associated group of atoms
C               (given on the next card) is to be constrained
C
C
C     -7.1-  
C       imol,iatom (16I5) (repeat if more than 8  centers
C
C    the list (ngrp long) of the atom indices of those atoms to be
C    constrained to the charge specified on the previous card.
C
c          blank to end
C
C------------------------------------------------------------------------
c  -8th-
c       intermolecular charge constraints  (atoms must sum to the 
c                                           specified value)
c       same format as indvidual molecule constraints
c       blank to end
c
C------------------------------------------------------------------------
c
c -9th- 
c      Multiple molecule constraints....constrain atoms on i to be
c                                            the same as on j
c      NGRP (I5) number of constaints 
C      (imol,iatom) (16I5) (repeat card if more than 8 groups)
c
c      blank to end
c
C------------------------------------------------------------------------
c
c
C  Unit 3 (qin) input of replacement charges if requested 
c  iqopt = 2,3
c
c       (8f10.6)  (i = 1,iuniq)
c
c-------------------------------------------------------------
c  Unit 4 input if new weight factors if requested
c
c    (i5)  nqwt  number of new weights to cycle thru
c    (f10.5)  new weights (nqwt lines)
c
c--------------------------------------------------------------
C
C Unit 10 input of ESP's  (mandatory)
C
C      natoms,nesp (2i5)
C              X , Y , Z  .   FORMAT (17X,3E16.7)
C      QUPOT , X , Y , Z  .   FORMAT (1X,4E16.7)
C
C          QUPOT = THE QUANTUM MECHANICAL ELECTROSTATIC
C                  POTENTIAL ( A.U )
C
C          X,Y,Z = THE COORDINATE AT WHICH THE POTENTIAL
C                  IS CALCULATED ( A.U )
C
C      NOTE : THE PROGRAM G80UCSF WRITES IN THIS FORMAT BUT THE
C             OUTPUT OF G90 MUST BE TRANSLATED (PROGRAM BOHR).
C
C--------------------------------------------------------------
C
c
c    usage: resp -i input -o output -p punch -q qin -t qout \
c                -e espot -w qwts -s esout 
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer icycle
      character*8   TITLE,    keywd
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      common /files/ input,output,qin,qout,punch,espot,qwts,esout,
     .               owrite
      character*80 input,output,qin,qout,punch,espot,qwts,esout
      character owrite
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10), keywd( 4,4)
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq),a(maxq,maxq),b(maxq),
     & qwtval(maxq),iwttyp(maxq),iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
C
C
c get the file names
c
      do jn = 1,maxq
         q0(jn) = 0.0d0
      enddo
      call filein
      call amopen(5,input,'O','F','R')
C
c read the atomic centers and q0's (readin), then read the potential inf
c
      call readin
c
c if its a multiple molecule run (nmol>0), do mult. mol. input
c
      if( nmol .gt. 1 )call mult_mol
c
c center & reorient molecule once in preparation for dipole & quadrupole
c
      if( nmol .eq. 1 ) call reornt
c
c read in the qm esp, forming the matrices apot(awt) and bpot(bwt)
c
      call matpot
c
c process the input (freezing, equivalencing charges)
c
      call data_prep
c
c set up cycle control structure: if icycle .gt. 0, come back to stmt 10
c subroutine cycle and read new "qwt"
c
c  icycle is initially set in readin 
c  subsequently decremented in cycle
c
      icycle= 0
c
c call cycle to see if we are supposed to cycle. reset icycle as needed
c
   10 continue
      call cycle(icycle)
c
      if( irstrnt .eq. 2 ) then
c
c if irstrnt= 2 then we just want to compare esp's to q0's
c
        do k= 1, iuniq
           qcal(k) = q0(k)
        enddo
        qwt= 0.0d0
      else
c 
c do the charge fitting
c
        call charge_opt
      endif
c
c calculate residuals sum-of-squares (chi-square) for the esp's
c
      call evlchi( ssvpot, bpot, apot, qcal, iuniq, maxq, chipot)
c
c now calculate and print dipole and quadrupole moments
c
      if(nmol .eq. 1) call elec_mom
c
c now punch short summary of charges & important evaluation criteria
c
      call pun_sum
c
      if( icycle .ne. 0 ) go to 10
c
c now calculate and print sum-of-squares, sigma, and rms
c
      call pot_out
c
      if( ioutopt .eq. 1 ) then
c 
c write the coords, old esps, new esps and residuals
c
         call wrt_pot
      endif
c
      END
C-----------------------------------------------------------------------
      SUBROUTINE readin
c
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      character *80 card
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10),keywd( 4,4)
      character*8   TITLE,    keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
C-----------------------------------------------------------------------
      common /files/ input,output,qin,qout,punch,espot,qwts,esout,
     .               owrite
      character*80 input,output,qin,qout,punch,espot,qwts,esout
      character owrite
c
      namelist /cntrl/
     &  ich, INOPT, ioutopt, iuniq, nmol ,IQOPT, ihfree,  qwt,
     &  irstrnt,iunits
C
      UNITS = 1.D0/0.529177249d0
c
      do i= 1,maxlgr
        do j= 1,maxq
           lgrcnt(i,j)= 0
        enddo
      enddo
c
      nlgrng= 0
c
c
c start of molecule input
c
      READ(5,1000) (TITLE(I),I=1,10)
 1000 FORMAT(10A8)
      WRITE(6,1010) (TITLE(I),I=1,10)
 1010 FORMAT(/,t2,'-----------------------------------------------',
     $       /,t2,'     Restrained ESP Fit 2.3  Amber 4.1',
     $       /,t2,'-----------------------------------------------',
     $       /,t2,10A8, 
     $       /,t2,'-----------------------------------------------'/) 
c 
c read in charge, number of charge centers, and control parameters 
c
      ich = 0
      iuniq = 0
      nmol = 1
      IQOPT = 0
      irstrnt = 1
      ihfree = 1
      qwt = 0.0005d0
      iunits = 0
c
      do 7 icard=1,1000
        read(5,'(a80)',end=10) card
        if (card(3:7).eq.'cntrl') go to 9
    7 continue
    9 backspace (5)
      read(5,cntrl,end=10)
      go to 20
   10 continue
         write(6,'(''Sorry, you must use namelist input'')')
         stop
   20 continue
      WRITE(6,1030) INOPT,ioutopt,nmol,IQOPT,
     $ ihfree,irstrnt,iunits,qwt
 1030 FORMAT(
     $       /t2,'inopt       = ',I5,'   ioutopt     = ',I5,
     $       /t2,'nmol        = ',I5,'   iqopt       = ',I5,
     $       /t2,'ihfree      = ',I5,'   irstrnt     = ',I5,
     $       /t2,'iunits      = ',i5,'   qwt         = ',f12.8)
c
c if nmol > 1, this is a multiple molecule run
c
c and so we should leave this routine now because routine mult_mol
c is responsible for the rest of the multiple-molecule reading in.
c
      if( nmol .gt. 1 ) then
        return
      endif
c
c read in fitting weight for q0 and esp point weighting
c
      READ(5,'(f10.0)') wtmol(1) 
      write(6,'(t2,''wtmol(1) = '',f12.6)')wtmol(1)
      READ(5,1000) (TITLE(I),I=1,10)
      write(6,'(t2,''subtitle:'',/ /,10a8)')(TITLE(I),I=1,10)
      read(5,'(2i5)')ich,iuniq
      write(6,'(''ich = '',i3,''  iuniq = '',i3)')ich,iuniq
      ibeg(1) = 1
      iend(1) = iuniq
      wtmol(1) = 1.d0
c
c
c readin  this is a single-molecule run 
c  
c readin  read in nuclear positions crd(i), initial charges q0(i)
c readin  convert angstroms to bohrs
c
      DO 12 I=1,iuniq
        READ(5,'(2i5)') IZAN(I),IVARY(I)
        write(6,'(3i5)') i,IZAN(I),IVARY(I)
   12 CONTINUE
 9708 FORMAT(2I5)
 9709 FORMAT(2I5)
c 
c
c read in lagrange (charge) constraints
c
      call lagrange( ICH, ibeg(1),iend(1), nmol)
c
c
c replacement initial charges q0 from unit IUNTQ0 if IQOPT=2
c
      if( IQOPT .gt. 1) then
        call amopen(3,qin,'O','F','R')
        WRITE(6,'(t2,''new q0 values to be read'' ,i4)') iuniq
c
c
        READ(3,'(8f10.6)',end=55,err=55) (Q0(I),i=1,iuniq)
        CLOSE(3)
        go to 59
c
c now is a simple trap for when not enough replacement q0 are given.
c tactic: just keep going with old q0 (assume IQOPT=2 was not intended)
c
   55 WRITE(6,'(t2,a,/,t2,a)')
     .   ' not enough (possibly none) q0 are given in file',
     .   ' ESP.Q0, so the remaining old ones will be used.'
   59 continue
      endif
c
c end of "replacement initial charge" section.
c
c set initial charges to 0; done if IQOPT=1
c
      if( IQOPT .eq. 1 ) then
        WRITE(6,'(t2,''IQOPT=1, all q0 values will be set to 0''/)')
        DO I=1,iuniq
           Q0(I)= 0.0d0
        enddo
      endif
c
      WRITE (6,9088)
 9088 FORMAT(/,2X,38(2H--),/,2X,'   ATOM   ',20X,  'COORDINATES',25X,
     $                             'CHARGE',/,23X,'X',13X, 'Y',
     $                            16X,'Z',/,2X,38(2H--))
      WRITE(6,9128)
 9128 FORMAT(2X,38(2H--),/)
      WRITE (6,9168) ICH,iuniq,qwt
 9168 FORMAT(/t2,'Charge on the molecule(ich) =',I5,
     $       /t2,'Total number of atoms (iuniq)     =',I5,
     $       /t2,'Weight factor on initial charge restraints(qwt)='
     $           ,d16.5,/)
c
c charge constraint info
c
      write(6,'(/t2,'' there are'',i3,'' charge constraints:'')')nlgrng
      write(6,'(/)')
      do 165 i= 1, iuniq
        write(6,'(t2,i5,5x,6(20i3,/,10x),20i3)')
     $   i,(lgrcnt(j,i), j= 1,nlgrng)
  165 continue
c     write(6,'(/t2'' charge: '',7f10.3)') (grpchg(i), i=1,nlgrng)
c
      IF (iuniq .LE. maxq ) RETURN
      WRITE (6,'(t2,''Number of atoms exceeds program dimensions'')')
      stop
      END
C------------------------------------------------------------------------
      SUBROUTINE mult_mol 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10), keywd( 4,4)
      character*8   TITLE,     keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
C-----------------------------------------------------------------------
      integer itmp(maxq),imoll(maxq)
      common /files/ input,output,qin,qout,punch,espot,qwts,esout,
     .               owrite
      character*80 input,output,qin,qout,punch,espot,qwts,esout
      character owrite
      DATA ONE,ZERO/1.0D0,0.0D0/
c
c this routine reads in multiple molecule input.  In routine readin
c it has already read the control variable input for the run, namely:
c ICH,IUNIQ,INOPT,IQOPT,ihfree, irstrnt, nlgrng, nmol
c where nmol > 0 caused this routine to be called.
c The input form for the other molecules is that their entire control
c decks are appended to the initial 2 control lines just read in,
c each control deck separated by a blank line, and then comes the
c multiple-molecule specific input, which is
c  - equivalencing of centers between molecules in the series
c  - the lagrange constraints to be applied between molecules
c
c the control characters read in the individual job decks are ignored
c except for ICH, and icntrs.  The lagrange (charge) constraints
c contained in the individual-molecule inputs ARE included.
c
c NOTE: the following are NOT implemented:
c       - output QM/calculated esp's (ioutOPT=1)
c
 1000 FORMAT(10A8)
 1010 FORMAT(10A8)
 1030 FORMAT(/,t2,'Total charge (ich):',I3,
     *       /,t2,'Number of centers:',I3)
 9708 FORMAT(I5,6x,4F12.6,i5,i5,f10.5)
      UNITS = 1.D0/0.529177249d0
      imol  = 0
      iuniq = 0
      mmlgr= nlgrng
      nlgrng= 0
c
      WRITE(6,'(/,t2,a,i3,a)')
     .   '%RESP-I-MULT_MOL,  multiple-molecule run of ', 
     .   nmol, ' molecules'
c
      do imol= 1, nmol
c
c read in the molecule weight and control variables
c  - the lagrange constraints to be applied between molecules
c
         READ(5,'(f10.5)') wtmol(imol)
         WRITE(6,'(/,t2,a,i3,a,f10.3)')
     .      'Reading input for molecule ', imol,
     .      ' weight:', wtmol(imol)
c
         READ(5,1000) (TITLE(I),I=1,10)
         WRITE(6,1010) (TITLE(I),I=1,10)
c
c read in charge, number of charge centers, and control parameters
c
         READ(5,'(2i5)') ICH,icntrs
         WRITE(6,1030) ICH, icntrs
c
c read in fitting weight for q0 and esp point weighting
c
c now some book-keeping: IUNIQ is the global variable for the total
c number of centers over all molecules.  The first center of this
c mol therefore starts in IUNIQ+1 and goes to IUNIQ+icntrs.
c
         ibeg(imol)= iuniq+1
         iend(imol)= iuniq+icntrs
c
c trap for having too many centers
c
         if( iend(imol) .gt. maxq ) then
           write (6,'(t2,''ERROR: more than '',i5,'' centers'')') maxq
           stop
         endif
c
c Read in nuclear positions crd(i), initial charges q0(i), and ivary(i)
c Since IVARY(i) is supposed to correspond to a center-number in the
c same molecule, this has to be adjusted to IVARY(i)+ibeg(imol)-1
c convert angstroms to bohrs if necessary
C
         DO I= ibeg(imol),iend(imol)
            READ(5,'(2i5)') IZAN(I), IVARY(I)
            write(6,'(3i5)') i,IZAN(I), IVARY(I)
            if (IVARY(i) .gt. 0) IVARY(i)= IVARY(i)+ibeg(imol)-1
         enddo
c
c now reset IUNIQ to IUNIQ+icntrs
c
         iuniq= iend(imol)
c
c now read in the lagrange (charge) restraints for this molecule
c
         call lagrange(ICH, ibeg(imol), iend(imol), imol)
      enddo
c
c end of molecule input, now do other preparation stuff
c
c read past a blank line after the final molecule job deck and then
c read in inter-molecule lagrange (charge) constraints (mmlgr).
c The "-99" for the total charge tells lgrange to drop the total charge
c constraint
c
      call lagrange( -99, 1, iuniq, imol)
c
c  this is for different types of charge resetting according to IQOPT :
c
c  replacement initial charges q0 from unit 3 if IQOPT=1
c
      if( IQOPT .gt. 1 ) then
        call amopen(3,qin,'O','F','R')
        WRITE(6,'(t2,'' since IQOPT=1,'',i4,'' new q0 values'')') iuniq
        WRITE(6,'(t2,'' will be read in from file ESP.Q0 (unit 3)'')')
c
c now read in replacement charges
c
        READ(3,'(8f10.6)',end=55,err=55) (Q0(I),i=1,iuniq)
        CLOSE(3)
        goto 59
c
c now is a simple trap for when not enough replacement q0 are given.
c tactic: just keep going with old q0 (assume IQOPT=2 was not intended)
c
   55   WRITE(6,'(t2,a,/,t2,a)')
     .     'not enough (possibly none) q0 are given in file',
     .     'ESP.Q0 (unit 3), so the remaining old ones will be used.'
   59   continue
      endif
c 
c end of "replacement initial charge" section.
c begin section: setting initial charges to 0; done if IQOPT=1
c
      if( IQOPT .eq. 1 ) then
        WRITE(6,'(/t2,''Iqopt =1: all q0 values will be set to 0'')')
        DO I=1,iuniq
           Q0(I)= 0.0
        enddo
      endif
c
c now carry out the inter-molecule equivalencing.  This is done by
c
c First : read the cards saying how many centers will be read in in the
c         next card. a zero means we have finished input
c
c Second: read the first-occurrence-in-each-molecule of the centers to
c        be equivalenced.  
c
c        The specifcations MUST be in ascending order.
c
c        The expanding of the centers within each
c        molecule is based on the IVARY values for the individual mol.
c
c        if ivary for mol 2+ is zero it is replaced with the atom number
c        of mol 1.  
c
      write(6,'(t2,''--------------------------------'')')
      write(6,'(t2,''reading mult_mol constraint info'')')
      write(6,'(t2,''--------------------------------'')')
  600 read(5,'(i5)') ntmp1
      if( ntmp1 .gt. 0 ) then
c
        read(5,'(16i5)') (imoll(j),itmp( j), j= 1,ntmp1)
        write(6,'(16i5)') (imoll(j),itmp( j), j= 1,ntmp1)
        do i = 1,ntmp1
           icntrs = ibeg(imoll(i)) - 1
           itmp(i) = icntrs + itmp(i)
        enddo
        do i= 2,ntmp1
           IVARY( itmp(i))= itmp(1)
        enddo
        go to 600
      endif
      do i= 1,iuniq
        ntmp1= IVARY(i)
        if( ntmp1 .gt. 0 ) then
          ntmp2= IVARY(ntmp1)
          if( ntmp2 .gt. 0 ) then 
             IVARY(i)= ntmp2
          endif
        endif
      enddo
c
      WRITE (6,9088)
 9088 FORMAT(/,2X,10(2H--),/,5x,'Atom   Ivary',/,2X,10(2H--))
      icnt = 1
      jcnt = 1
      DO 1100 IAT = 1,iuniq
        WRITE (6,'(2i5)')IZAN(IAT),ivary(iat)
      jcnt = jcnt + 1
      if (jcnt.gt.iend(icnt))then
         write(6,'( )') 
         icnt = icnt + 1
      endif
 1100 CONTINUE
      WRITE(6,9128)
 9128 FORMAT(2X,38(2H--),/)
      WRITE (6,9168) iuniq,qwt
 9168 FORMAT(/t2,'Total number of atoms      =',I5,
     *       /t2,'Weight factor on initial charge restraints=',f10.6,/)

c charge constraint info
 
      write(6,'(/,t2,''There are'',i3,'' charge constraints'')')nlgrng
      return
      end
C-------------------------------------------------------------------
      SUBROUTINE lagrange( ncharge, ifirst, ilast, imol)
c
c read in and assign lagrange constraint pointers
c
c called from "readin" and "mult_mol"
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10),keywd( 4,4)
      character*8   TITLE,    keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
C-----------------------------------------------------------------------
      integer itmp(maxq),imoll(maxq)
c
c
c section: read in explicit Lagrange constraints on charge groups
c 
c 
c read in constraints: charge, center to which it applies
c
  10    continue
          read(5,'(i5,f10.5)') ntmp, gtemp
          if(ntmp.eq.0) go to 20
          nlgrng= nlgrng + 1
          if( (nlgrng) .gt. maxlgr) then
             write(6,'(t2,a,i3,/,t2,a,i3)')
     .          ' Too many charge-group constraints', nlgrng,
     .          ' Maximum allowed: ', maxlgr
             stop
          endif
          grpchg(nlgrng) = gtemp 
          read(5,'(16i5)') (imoll(j),itmp( j), j= 1,ntmp)
          write(6,'(16i5)') (imoll(j),itmp( j), j= 1,ntmp)
c
          do ii = 1,ntmp
             jmol = imoll(ii)
             icntrs = ibeg(jmol) - 1
             itmp(ii) = icntrs + itmp(ii)
          enddo
          do j= 1, ntmp
            if( itmp(j) .gt. 0 ) then
              itmp(j)= itmp(j) + ifirst -1
              lgrcnt(nlgrng, itmp(j))= 1
            elseif( itmp(j) .lt. 0 ) then
              itmp(j)= itmp(j) - ifirst +1
              lgrcnt(nlgrng, itmp(j))= -1
            endif
          enddo
        go to 10
   20 continue
c
c as long as ncharge is not -99, implement the "total charge" constraint
c
      if( ncharge .gt. -99 ) then
        nlgrng= nlgrng + 1
c
        if( nlgrng .gt. maxlgr) then
           write(6,'(t2,a,i3,/,t2,a,i3)')
     .          ' Too many charge-group constraints', nlgrng,
     .          ' Maximum allowed: ', maxlgr
           stop
        endif
        grpchg(nlgrng)= float( ncharge)
        do  j= ifirst, ilast
            lgrcnt(nlgrng,j)= 1
        enddo
      endif
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE matpot
c
c read in the electrostatic potential points used in the fitting,
c building up as we go the matrices for LU decomposition
c
c called from Main
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10), keywd( 4,4)
      character*8   TITLE,         keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
      common /files/ input,output,qin,qout,punch,espot,qwts,esout,
     .               owrite
      character*80 input,output,qin,qout,punch,espot,qwts,esout
      character owrite
c
      DO k=1,iuniq
        bpot(k)= 0.0d0
        bwt(k) = 0.0d0
        DO j=1,iuniq
           apot(j,k)= 0.0d0
           awt(j,k) = 0.0d0
        enddo
      enddo
c
      call amopen(10,espot,'O','F','R')
c
      ssvpot= 0.0d0
      vavrg = 0.0d0
      vavtmp = 0.0d0
      ioff = 1
c
      if( nmol .gt. 0 ) then
        inmol= nmol
      else
        inmol= 1
        ibeg(1)= 1
        iend(1)= iuniq
        wtmol(1)= 1.0d0
      endif
c
      do imol= 1, inmol
        read(10,'(2i5)') inat,nesp
        WRITE(6,'(/,t2,a,i3,/,t2,a,i5,/,t2,a,i5)')
     .     'Reading esp"s for molecule ',imol,
     .     'total number of atoms      = ',inat,
     .     'total number of esp points = ',NESP
        WRITE(6,'(/ /,a)')
     .     ' center     X       Y       Z '
        do i = 1,inat
            read(10,52)crd(1,ioff),crd(2,ioff),crd(3,ioff)
            write(6,152)i,crd(1,ioff),crd(2,ioff),crd(3,ioff)
            ioff = ioff + 1
        enddo
   52      format(17X,3e16.7) 
  152      format(1X,i4,3e16.7) 
c
c build up matrix elements Ajk according to (SUMi 1/Rik SUMj 1/Rij)
c
        DO i= 1,nesp
           read(10, 53, err=940, end=930) espi, xi, yi, zi
   53      format(1X,4e16.7) 
           wt= wtmol(imol)
           wt2= wt*wt
           vavtmp = vavtmp + wt*espi
           ssvpot = ssvpot + wt2*espi*espi
           vavrg  = vavrg  + vavtmp/float(NESP)
           DO k= ibeg(imol),iend(imol)
              rik = sqrt( (xi - CRD(1,k))**2 + (yi - CRD(2,k))**2 +
     $                                      (zi - CRD(3,k))**2)
              rik = 1.d0/rik
              bpot(k)  = bpot(k) +     espi*rik
              bwt(k)   = bwt(k)  + wt2*espi*rik
              apot(k,k)= apot(k,k) +     rik*rik 
              awt(k,k) = awt(k,k)  + wt2*rik*rik 
              DO j= k+1, iend(imol)
                 rij = sqrt((xi - CRD(1,j))**2+(yi - CRD(2,j))**2+
     $                                         (zi - CRD(3,j))**2)
                 rij = 1.d0/rij
                 apot(j,k)= apot(j,k) +     rij*rik
                 awt(j,k) = awt(j,k)  + wt2*rij*rik
              enddo
           enddo
        enddo
      enddo
c
c symmetrize the potenitial and weighted potential matrices
c
      DO k=1,iuniq
         DO j= k+1,iuniq
            awt(k,j)= awt(j,k)
            apot(k,j)= apot(j,k)
         enddo
      enddo
      close(10)
      write(6,900)ssvpot
  900 format(t2,'Initial ssvpot =',f10.3/)
      return
c
  930 write(6,'(5X,''premature end of potential file'')')
      stop
  940 write(6,'(5X,a)') 
     .    'Error in reading potential input file (1X,4e16.7)'
      stop
      end
c-----------------------------------------------------------------------
      SUBROUTINE data_prep
c
c  setup pointers for groups of charges  based on "ivary" info
c
c  called from Main
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10), keywd( 4,4)
      character*8   TITLE,     keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
C***********************************************************************
c
c
c begin section: set lists for combined and frozen charges
c
c IVARY(i) = 0, it is a new charge center to be fitted
c IVARY(i) =+n, it is a charge center to be fitted with center n
c                    (center n must be a previous center entered with 
c                    IVARY(n) = 0
c IVARY(i) =-n, it is a frozen charge center to be kept at q0(i)
c
c************************************************************************
      NAT = 0
      do i= 1,iuniq
        if(IVARY(i) .eq. 0) then
           NAT= NAT + 1
           iqcntr( i)= NAT
        elseif(IVARY(i) .gt. 0) then
           iqcntr(i)= iqcntr(IVARY(i))
c
           if( iqcntr( IVARY(i)) .gt. NAT ) then
              write(6,'(t2,''data_prep: charge vary input is screwy'')')
              stop
           endif
        else
           iqcntr( i)= -1
        endif
      enddo
c
      WRITE (6,'(/t2,''Number of unique UNfrozen centers='',i5)') NAT
      if( NAT .eq. 0 ) then 
            write(6,'(t2,''ALL charges are frozen!!!'')')
      endif
c
c finish off list with Lagrange constraints
c
      do 150 i= 1,nlgrng
        iqcntr(iuniq+i)= nat + i
  150 continue
c
c set NATPL1 to the total number of row elements 
c (charges to be independantly fit + constraints )
c in fitting matrix
c
      natpl1= nat + nlgrng
c
c     write(6,'(t2,''iqcntr: '')')
c     write(6,'(20i3)')(iqcntr(jn),jn=1,iuniq+nlgrng)
c
c done adding Lagrange constraints to elements list
c
c red in charges must now be averaged 
c a posteriori averaging of replacement charges according
c              to current IVARY charge-combining pointers
c
      if( IQOPT .eq. 3 ) then
        do i= 1,iuniq-1
          qcntrs= q0(i)
          tmpctr= 1.0d0
          do j= i+1,iuniq
            if( IVARY(j) .eq. i ) then
              qcntrs= qcntrs + q0(j)
              tmpctr= tmpctr + 1.0d0
            endif
          enddo
          if(tmpctr .gt. 0.99d0 )then
             qcntrs= qcntrs/tmpctr
             q0(i)= qcntrs
             do j= i+1,iuniq
                if( IVARY(j) .eq. i ) q0(j)= qcntrs
             enddo
          endif
        enddo
      endif
      RETURN
      END
C-------------------------------------------------------------------
      SUBROUTINE matbld
c
c called from "chgopt"
c 
c build up matrices for LU decomposition:
c
c   stage 1: copy weighted matrices awt and bwt to work arrays awork and bwork
c            (which are destroyed in the LU decomp & back subst)
c
c   stage 2: if charge restraints are to be included,
c            then modify awork and bwork appropriately
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10), keywd( 4,4)
      character*8   TITLE,           keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
C
      do k=1,iuniq
        b(k) = bwt(k)
        do J=1,iuniq
           a(j,k) = awt(j,k)
        enddo
      enddo
c
c fill in the final columns & rows of A with the Lagrange
c constraints which keep the charge on groups of atoms to a 
c constant 
c
c note index counters!
c
      do i=1,nlgrng
        b(iuniq+i)  = grpchg(i)
        do j= 1,iuniq+nlgrng
          a(iuniq+i,j)= float(lgrcnt(i,j))
          a(j,iuniq+i)= float(lgrcnt(i,j))
        enddo
      enddo
  450 continue
c
c add restraint to initial charge q0(i):
c
      call rstran
c
c build awork and bwork based on "combined and frozen centers" info:
c
c 1) frozen centers do not appear in the matrix of fitted charges
c 2) combined centers appear as one single charge center for fitting
c
c first, since we accumulate values, zero out awork & bwork up to natpl1
c (the independant + contraint number):
c
      do i= 1,natpl1
        bwork( i)   = 0.0d0
        awork( i, i)= 0.0d0
        do j= i+1,natpl1
          awork( j, i)= 0.0d0
          awork( i, j)= 0.0d0
        enddo
      enddo
c
c loop over all centers, building awork & bwork from A and
c B based on iqcntr: for each center, iqcntr(i) dictates which of
c the fitted charges it is and therefore where it goes in the matrices.
c If iqcntr(j) < 1, this center is a frozen charge and it is skipped as
c far as forming a row in awork, and its esp contribution is subtracted
c from bwork to take care of it's awork jth column-element for each i.
c
      do i= 1,iuniq+nlgrng
        icntr= iqcntr(i)
        if(icntr .ge. 1) then
c i is "active" 
          bwork(icntr)= bwork(icntr) + b(i)
          do j = 1,iuniq+nlgrng
            jcntr = iqcntr( j)
            if( jcntr .ge. 1 ) then
c j i "active"
              awork(icntr,jcntr)= awork(icntr,jcntr) + a(i,j)
            else
c j is  a frozen charge
              bwork(icntr)= bwork(icntr)       - q0(j)*a(i,j)
            endif
          enddo
        endif
      enddo
c
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE rstran
c 
c routine to assign the retraint weights 
c to the diagonal of A and to B
c 
c called from "matbld"
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10),keywd( 4,4)
      character*8   TITLE,    keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
C----------------------------------------------------------------------
c two kinds of restraint are available:
c
c a) a harmonic restraint to the initial charge.  Fine as long as there
c  aren't any large charges that SHOULD be large... these really feel a
c  strong force if they are restrained to a low value.
c
c b) a hyperbolic restraint to a charge of 0.  This gets asymptotic at
c  "large" values, so "large" charges aren't pulled down any stronger
c  than some (reasonable) limiting force.  This is a non-linear
c  weighting function, so the fit procedure is iterative.
c
c other options for restraints to initial charge q0(i):
c if requested, restrain the charges by modifying the sum-of-squares
c cost function derivative.  The scheme for doing this is as follows:
c
c if control variable ihfree > 0, let hydrogen charges float free
c                                   (i.e. reset their qwtval to 0.0).
c
c-----------------------------------------------------------------------
c
      do i = 1,iuniq
        qwtval(i)= qwt
        if(ihfree.gt.0.and. IZAN(i).eq.1) then
           qwtval(i)= 0.0D0
        endif
c
        if(irstrnt .eq. 0) then
          a(i,i)= a(i,i) + qwtval(i)
c
c q0 has the initial and/or frozen  charge
c
          b(i)  = b(i)   + qwtval(i)*q0(i)
        elseif (irstrnt .gt. 0 .and. qwtval(i) .gt. 0.1d-10) then
c
c use analytic gradient of the hyperbola 
c
c qcal has the current (calculated) charge
c
          qwtval(i) = qwt/sqrt(qcal(i)*qcal(i) + 0.01d0)
          a(i,i)= a(i,i) + qwtval(i)
        endif
      enddo
c
c  if all qwtval(i) are 0.0, no restraints so set irstrnt= -1
c
      ihit = 0
      do i= 1,iuniq
         if(qwtval(i) .gt. 0.1d-10 ) ihit = 1
      enddo
      if(ihit .eq. 0) irstrnt = -1
      return
      end
c------------------------------------------------------------------------
      SUBROUTINE cycle( icycle)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      integer icycle, nqwt
      character*8   TITLE,    keywd
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/  NAT,IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/ TITLE(10),keywd( 4,4)
      COMMON/ESPCOM/ apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     &               qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/   q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/ awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
C-----------------------------------------------------------------------
      save nqwt
      common /files/ input,output,qin,qout,punch,espot,qwts,esout,
     .               owrite
      character*80 input,output,qin,qout,punch,espot,qwts,esout
      character owrite
c
      DATA ZERO/0.0D0/
c
c
c if INOPT=1, reads values new of qwt from unit 4 and 
c writes a summary of the resulting fit, 
c as well as a list of the fit charges.
c
      if( INOPT .lt. 1 ) return
      if( icycle .eq. 0 ) then
         nqwt= 0
c
c first pass
c
         call amopen(4,qwts,'O','F','R')
         read(4,'(i5)') nqwt
         if(nqwt.eq.0) then
            write(6,'('' subroutine icycle: INPUT ERROR...'')')
            write(6,'(a,a)')
     .         ' INOPT=1 so qwt cycling expected, but',
     .         ' reading non-zero nqwt failed'
            return
         endif
   37    READ(4,'(f10.5)') qwt
         icycle= 1
         write(6,'(t2,''cycle   1: weighting factor= '',f10.4)')qwt
      else
c
c rest of the time
c
         READ(4,'(f10.5)') qwt
         icycle= icycle + 1
         write(6,'(/t2,''cycle'',i4,'': weighting factor= '',f10.4)')
     &             icycle,qwt
         if( icycle .ge. nqwt ) icycle= 0
      endif
      return
      end
c-----------------------------------------------------------------------------
      SUBROUTINE reornt
c
c translates molecule to center of mass and 
c reorients it along principal axes
c in preparation for dipole and quadrupole
c  moment calculation.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10), keywd( 4,4)
      character*8   TITLE,     keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
C-----------------------------------------------------------------------
      DIMENSION IATOM(maxq),WT(35)
C
C     ---- ATOMIC WEIGHT ARRAY FOR CENTER OF MASS ----
C
      DATA ZERO/0.0D0/,BOHR/0.52917725D0/
      DATA (WT(I),I=1,21) /
     *                     0.0000D0, 1.0079D0, 4.0026D0,
     *                     6.9410D0, 9.0122D0,10.8100D0,12.0110D0,
     *                    14.0067D0,15.9994D0,18.9984D0,20.1700D0,
     *                    22.9897D0,24.3050D0,26.9815D0,28.0860D0,
     *                    30.9738D0,32.0600D0,35.4530D0,39.9480D0,
     *                    39.0900D0,40.0800D0 /
C
C     ----- INITIALISE SOME VARIABLES -----
C
      XC = ZERO
      YC = ZERO
      ZC = ZERO
      DO 10 I = 1,iuniq
         IATOM(I) = IZAN(I)+1
   10 CONTINUE
C
C     ----- CALCULATE THE CENTER OF MASS -----
C
      CALL CMASS(XC,YC,ZC,CRD,IATOM,iuniq,WT,maxq)
      CMAS(1)= XC
      CMAS(2)= YC
      CMAS(3)= ZC
c convert center of mass from bohr to angstroms
      DO I = 1,3
         CMAS(I) = CMAS(I)*BOHR
      enddo
C
C     ----- MOVE THE ORIGIN TO CENTER OF MASS AND
C           FORM NEW COORDINATE -----
C
      CALL CMOVE(XC,YC,ZC,CRD,iuniq,CO,maxq)
C
C     ----- CALCULATE THE MOMENT OF INERTIA -----
C
      CALL MOMIN(CO,IATOM,iuniq,WT,maxq)
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE elec_mom
c
C ROUTINE TO CALCULATE THE DIPOLE AND QUADRUPOLE MOMENTS 
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10), keywd( 4,4)
      character*8   TITLE,     keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
C***********************************************************************
      ZERO=0.0D0
      BOHR=0.52917725D0
      debye=2.541765d0
C
C calculate dipole moment
C
      DIPOL(1) =ZERO
      DIPOL(2) =ZERO
      DIPOL(3) =ZERO
      DO I =1,iuniq
        DIPOL(1) =DIPOL(1) +QCAL(I)*CRD(1,I)
        DIPOL(2) =DIPOL(2) +QCAL(I)*CRD(2,I)
        DIPOL(3) =DIPOL(3) +QCAL(I)*CRD(3,I)
      enddo
      CONTINUE
      dipmom= dsqrt( DIPOL(1)*DIPOL(1) + DIPOL(2)*DIPOL(2) +
     *                                               DIPOL(3)*DIPOL(3))
 
C     ----- CALCULATE THE QUADRUPOLE MOMENT -----
 
      CALL QUADM(CRD,QCAL,iuniq,QUAD,maxq)
 
C convert dipoles from a.u. to debyes, and quadrupoles to debye*angstroms
 
      DO I = 1,3
         DIPOL(I) = DIPOL(I)*debye
         QUAD(I) = QUAD(I)*debye*BOHR
         QUAD(I+3) = QUAD(I+3)*debye*BOHR
      enddo
      dipmom= dipmom*debye
 
      RETURN
      END
C***********************************************************************
      SUBROUTINE CMASS(XC,YC,ZC,C,IATOM,NATOM,WT,maxq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
c  called from "reornt"
 
C      THIS SUBROUTINE CALCULATES THE CENTER OF MASS OF THE
C      MOLECULE.
 
      DIMENSION C(3,maxq),IATOM(maxq),WT(2)
 
      DATA ZERO/0.0D0/
 
      SUMX = ZERO
      SUMY = ZERO
      SUMZ = ZERO
      SUM  = ZERO
      DO 3 I =1,NATOM
        INDEX=IATOM(I)
        SUMX = SUMX + C(1,I) * WT(INDEX)
        SUMY = SUMY + C(2,I) * WT(INDEX)
        SUMZ = SUMZ + C(3,I) * WT(INDEX)
    3 SUM = SUM + WT(INDEX)
      XC = SUMX/SUM
      YC = SUMY/SUM
      ZC = SUMZ/SUM
      RETURN
      END
C***********************************************************************
      SUBROUTINE CMOVE(XC,YC,ZC,C,NATOM,CO,maxq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  
c called from "reornt"
 
C     THIS SUBROUTINE MOVES THE ORIGIN TO THE CENTER OF MASS.
 
      DIMENSION C(3,maxq),CO(3,maxq)
C
      DO 5 I=1,NATOM
      CO(1,I) = C(1,I) - XC
      CO(2,I) = C(2,I) - YC
      CO(3,I) = C(3,I) - ZC
    5 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE MOMIN(CO,IATOM,NATOM,WT,maxq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
c 
c called from  "reornt"
c
C     THIS SUBROUTINE CALCULATES THE MOMENTS OF INERTIA AND
C     THE PRINCIPAL AXES OF ROTATION OF THE MOLECULE.  IT THEN
C     REORIENTS THE MOLECULE ALONG THE PRINCIPAL AXES OF
C     ROTATION (WITH THE ORIGIN AT THE CENTER OF MASS) IN
C     PREPARATION FOR CALCULATION OF QUADRAPOLE MOMENT
C     COMPONENTS.
C
C
      DIMENSION D(3),AIN(3,3),S(3,3)
      DIMENSION CO(3,maxq),IATOM(maxq),WT(2)
C
      DATA ZERO /0.0D0/
C
      SXX=ZERO
      SYY=ZERO
      SZZ=ZERO
      SXY=ZERO
      SXZ=ZERO
      SYZ=ZERO
      DO 10 I=1,NATOM
        INDEX=IATOM(I)
        XX=CO(1,I)*CO(1,I)
        YY=CO(2,I)*CO(2,I)
        ZZ=CO(3,I)*CO(3,I)
        SXX=SXX + (WT(INDEX)*(YY + ZZ))
        SYY=SYY + (WT(INDEX)*(XX + ZZ))
        SZZ=SZZ + (WT(INDEX)*(XX + YY))
        SXY=SXY + (WT(INDEX)*CO(1,I)*CO(2,I))
        SXZ=SXZ + (WT(INDEX)*CO(1,I)*CO(3,I))
   10 SYZ=SYZ + (WT(INDEX)*CO(2,I)*CO(3,I))
      AIN(1,1)=SXX
      AIN(1,2)=-SXY
      AIN(1,3)=-SXZ
      AIN(2,1)=AIN(1,2)
      AIN(2,2)=SYY
      AIN(2,3)=-SYZ
      AIN(3,1)=AIN(1,3)
      AIN(3,2)=AIN(2,3)
      AIN(3,3)=SZZ
C
C     ----- CALCULATE PRINCIPAL AXES OF INERTIA -----
C     ----- EV ARE PRINCIPLE MOMENTS OF INERTIA -----
C
      CALL DIAGM(AIN,S)
C
C      ----- ARRANGE SMALLEST TO LARGEST -----
C
      DO 30 I=1,2
        I1=I + 1
        DO 30 K=I1,3
          IF(AIN(I,I)-AIN(K,K))30,30,20
   20     RET=AIN(I,I)
          AIN(I,I)=AIN(K,K)
          AIN(K,K)=RET
          DO 25 L=1,3
            D(L)=S(L,I)
            S(L,I)=S(L,K)
   25     S(L,K)=D(L)
   30 CONTINUE
      DO 35 I=1,NATOM
        XT=S(1,1)*CO(1,I) + S(2,1)*CO(2,I) + S(3,1)*CO(3,I)
        YT=S(1,2)*CO(1,I) + S(2,2)*CO(2,I) + S(3,2)*CO(3,I)
        ZT=S(1,3)*CO(1,I) + S(2,3)*CO(2,I) + S(3,3)*CO(3,I)
        CO(1,I)=XT
        CO(2,I)=YT
   35 CO(3,I)=ZT
C
C MOMENT OF INERTIA WITH THE PRINCIPAL AXES
C
      XPI=AIN(1,1)
      YPI=AIN(2,2)
      ZPI=AIN(3,3)
      RETURN
      END
C***********************************************************************
      SUBROUTINE DIAGM(AIN,S)
c
c  called from "momin"
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     THIS SUBROUTINE PERFORMS DIAGONALIZATION OF A ROTATION
C     MATRIX BY JACOBI ROTATION METHOD.  IT IS NECESSARY TO
C     OBTAIN THE MOMENTS OF INERTIA AND COORDINATES OF THE
C     PRINCIPAL ROTATIONAL AXES.  THIS SUBROUTINE IS FROM THE
C     DIAGONALIZATION SUBROUTINE IN MM2 WRITTEN BY M. MILLER.
C
      DIMENSION AIN(3,3),S(3,3)
C
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/
      DATA THREE,PT5,TENP8/3.0D0,0.5D0,1.0D+08/
C
      I=1
      K=1
      VI=ZERO
   10 BD=ZERO
   20 DO 80 I=1,3
   30    DO 70 K=1,3
   40       IF(I-K)60,50,60
   50       S(I,K)=ONE
            GO TO 70
   60       S(I,K)=ZERO
            VI=VI + (AIN(I,K)*AIN(I,K))
   70    CONTINUE
   80 CONTINUE
      VI=DSQRT(VI)
      VF=VI/TENP8
      SGM=THREE
      V=VI
      IF(VI)280,280,90
   90 V=V/SGM
  100 M=2
  110 L=1
  120 IF(DABS(AIN(L,M))-V)210,130,130
  130 BD=ONE
      ALM=-AIN(L,M)
      UM=PT5*(AIN(L,L)-AIN(M,M))
      OMG=ALM/(DSQRT((ALM*ALM)+(UM*UM)))
      IF(UM)140,150,150
  140 OMG=-OMG
  150 SN=OMG/(DSQRT(TWO*(ONE+(DSQRT(ONE-OMG*OMG)))))
      CS=DSQRT(ONE-SN*SN)
      I=1
  160 C1=AIN(I,L)
      C2=AIN(I,M)
      AIN(I,L)=C1*CS-C2*SN
      AIN(I,M)=C1*SN+C2*CS
      C1=S(I,L)
      C2=S(I,M)
      S(I,L)=C1*CS-C2*SN
      S(I,M)=C1*SN+C2*CS
      IF(I-3)170,180,170
  170 I=I+1
      GO TO 160
  180 C1=AIN(L,L)
      C2=AIN(M,M)
      C3=AIN(L,M)
      C4=AIN(M,L)
      AIN(L,L)=C1*CS-C4*SN
      AIN(M,M)=C2*CS+C3*SN
      AIN(L,M)=C3*CS-C2*SN
      AIN(M,L)=AIN(L,M)
      I=1
  190 AIN(L,I)=AIN(I,L)
      AIN(M,I)=AIN(I,M)
      IF(I-3)200,210,200
  200 I=I+1
      GO TO 190
  210 IF(L-M+1)220,230,220
  220 L=L+1
      GO TO 120
  230 IF(M-3)240,250,240
  240    M=M+1
         GO TO 110
  250    IF(BD-ONE)260,270,260
  260       IF(V-VF)280,280,90
  270          BD=ZERO
               GO TO 100
  280 CONTINUE
      RETURN
      END
C------------------------------------------------------------------------
      SUBROUTINE QUADM( C, q, NATOM, QUAD, maxq)
 
c  called from "elec_mom"
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
C   THIS SUBROUTINE CALCULATES THE COMPONENTS OF THE QUADROPOLE MOMENT 
 
      DIMENSION C(3,maxq), q(maxq), QUAD(6)
 
      QXX=0.0d0
      QYY=0.0d0
      QZZ=0.0d0
      QXY=0.0d0
      QXZ=0.0d0
      QYZ=0.0d0

      DO 7 I=1,NATOM
         X2 = C(1,I)*C(1,I)
         Y2 = C(2,I)*C(2,I)
         Z2 = C(3,I)*C(3,I)
         RQ = X2 + Y2 + Z2
            QXX = QXX + Q(I)*( 3.d0* (C(1,I)*C(1,I)) - RQ )
            QYY = QYY + Q(I)*( 3.d0* (C(2,I)*C(2,I)) - RQ )
            QZZ = QZZ + Q(I)*( 3.d0* (C(3,I)*C(3,I)) - RQ )
            QXY = QXY + Q(I)*( 3.d0* (C(1,I)*C(2,I)) - RQ )
            QXZ = QXZ + Q(I)*( 3.d0* (C(1,I)*C(3,I)) - RQ )
            QYZ = QYZ + Q(I)*( 3.d0* (C(2,I)*C(3,I)) - RQ )
    7 CONTINUE
C
C STORE IN QUAD
C
      QUAD(1) = QXX
      QUAD(2) = QYY
      QUAD(3) = QZZ
      QUAD(4) = QXY
      QUAD(5) = QXZ
      QUAD(6) = QYZ
C
      RETURN
      END
C***********************************************************************
      subroutine filein
c
      implicit double precision (a-h,o-z)
c
c     OUTPUT: (to common)
c
      common /files/ input,output,qin,qout,punch,espot,qwts,esout,
     .               owrite
      character*80 input,output,qin,qout,punch,espot,qwts,esout
      character owrite
c
      character *80 arg
c
      integer iarg, narg
      input = 'input'       ! -i
      output = 'output'     ! -o
      punch = 'punch'       ! -p
      qin= 'qin'            ! -q
      qout= 'qout'          ! -t 
      espot = 'espot'       ! -e
      qwts = 'qwts'         ! -w
      esout = 'esout'       ! -s
c
c     --- default for output files: 'N'ew
c
      owrite = 'N'
c
c     --- get com line arguments ---
c
      indx = iargc()
      iarg = 0
      if (indx.eq.iarg) goto 20
   10 continue
           iarg = iarg + 1
           call getarg(iarg,arg)
           if (arg .eq. '-O') then
                owrite = 'U'
           elseif (arg .eq. '-i') then
                iarg = iarg + 1
                call getarg(iarg,input)
           elseif (arg .eq. '-o') then
                iarg = iarg + 1
                call getarg(iarg,output)
           elseif (arg .eq. '-p') then
                iarg = iarg + 1
                call getarg(iarg,punch)
           elseif (arg .eq. '-q') then
                iarg = iarg + 1
                call getarg(iarg,qin)
           elseif (arg .eq. '-t') then
                iarg = iarg + 1
                call getarg(iarg,qout)
           elseif (arg .eq. '-e') then
                iarg = iarg + 1
                call getarg(iarg,espot)
           elseif (arg .eq. '-s') then
                iarg = iarg + 1
                call getarg(iarg,esout)
           elseif (arg .eq. '-w') then
                iarg = iarg + 1
                call getarg(iarg,qwts)
           else
                if (arg .eq. ' ') go to 20
                write(6,'(/,5x,a,a)') 'unknown flag: ',arg
                write(6,9000)
                stop
           endif
      if (iarg .lt. indx) go to 10
c
   20 continue
      narg = iarg - 1
      if (narg .lt. 2) then
           write(6,9000)
           stop
      endif
c
      if(output .ne. 'screen') then
          call amopen(6,output,owrite,'F','W')
      endif
      return
 9000 format(/,t2,
     $'usage: resp [-O] -i input -o output -p punch -q qin -t qout','
     $   -e espot -w qwts -s esout ')
      end
c------------------------------------------------------------------------
      subroutine evlchi(ssvpot,bpot,apot,qcal,iuniq,maxq,chipot)

      implicit double precision (a-h,o-z)
      dimension bpot(maxq), apot(maxq,maxq), 
     $         qcal(maxq)
c
c called from Main
c
c Evaluate chi-square for linear function  yclci= sum-j(paramj*termij),
c where j= number of terms, paramj is the coefficient to termij, and
c chi-square is the merit function: chi-square= sum-i((yi-yclci)**2),
c where i is the number of data points for which y is known.
c
c To avoid going through all i data points every time, the chi-square
c function is expanded into: chi-square= sum-i(yi**2 - 2*yi*yclci + yclci**2),
c and re-expressed as sum-i(yi**2) - 2*sum-i(yi*yclci) + sum-i(yclci**2).
c The first term is calculated once as the data is read in (sum-i(yi**2)== ssy).
c The second and third term depend upon the parameters: for each parameter
c paramj, sum-i(yi*termij) is the jth element in xprod, and sum-i(yclci**2) is
c a row in apot where element ajk= sum-i(termij*termik).  These elements
c are also built up once as the (possibly large) set of data is read in.
c
      cross= 0.0d0
      ssyclc= 0.0d0
      do j= 1,iuniq
        cross= cross + qcal(j)*bpot(j)
        do k= 1,iuniq
          ssyclc= ssyclc + qcal(j)*qcal(k)*apot(j,k)
        enddo
      enddo
c
      chipot = ssvpot - 2.0d0*cross + ssyclc
c
      return
      end
c------------------------------------------------------------------------
      subroutine pun_sum
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995,1997                 **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
c
c  called from Main
c
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10),keywd( 4,4)
      character*8   TITLE,    keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
C***********************************************************************
      common /files/ input,output,qin,qout,punch,espot,qwts,esout,
     .               owrite
      character*80 input,output,qin,qout,punch,espot,qwts,esout
      character owrite
c
c
      QCRTRN= SQRT(chipot/ssvpot)
      ssvtot = ssvpot
c
c calculate standard error of estimate and correlation coefficients
c
      SIGMA = SQRT(chipot/FLOAT(NESP))
c
c punch name, one-line summary, and charges:
c
      call amopen(7,punch,owrite,'F','W')
      WRITE(7,1110) (TITLE(I),I=1,10)
 1110 FORMAT(/,10A8)
      write(7, 1230) IQOPT, irstrnt, ihfree,qwt
 1230 format(/,'iqopt   irstrnt  ihfree     qwt'/,3(i3,4x),f12.6)
      write(7,1235)  QCRTRN, dipmom, (QUAD(I), I=1,3)
 1235 format(/t2,'rel.rms   dipole mom       Qxx      Qyy      Qzz',/
     $ t2,5f10.5)
      WRITE(7,1200)
 1200 FORMAT(/,10X,'Point charges before & after optimization',/,4X,
     .  'NO',3X,'At.No.',4x,'q0',11X,'q(opt)   IVARY  d(rstr)/dq ')
      icnt = 1
      jcnt = 1
      DO j=1,iuniq
      WRITE(7,1210) j,izan(j),q0(j), qcal(j), IVARY(j), qwtval(j)
      jcnt = jcnt + 1
      if (jcnt.gt.iend(icnt))then
         write(6,'( )') 
         icnt = icnt + 1
      endif
      enddo
 1210 format(2x,i4,2x,i4,1x,f10.6,5x,f10.6, i5, f12.6, f12.3)
c
      call amopen(19,qout,owrite,'F','W')
      write(19,'(8f10.6)') (Qcal(I),i=1,iuniq)
      close(unit=19)
c
      RETURN
      END
C---------------------------------------------------------------------
      SUBROUTINE POT_OUT
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995,1997                 **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
c
c  called from Main
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, iuniq,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10),keywd( 4,4)
      character*8   TITLE,    keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
      DATA ZERO/0.0D0/,AU2CAL/627.5095D0/
C
C
C     ---- PRINT THE OPTIMIZED CHARGES AND COORDINATES ----
C
      WRITE(6,1110) (TITLE(I),I=1,10)
c
c print the charges
c
      WRITE(6,1200)
      icnt = 1
      jcnt = 1
      chge = 0.0
      DO j=1,iuniq
         WRITE(6,1210) j,izan(j),q0(j), qcal(j), IVARY(j), qwtval(j)
         chge = chge + qcal(j)
         jcnt = jcnt + 1
         if (jcnt.gt.iend(icnt))then
            write(6,'( )') 
            icnt = icnt + 1
         endif
      enddo
      write(6,'(t2,''Sum over the calculated charges: '',f10.3)')chge
c
c
      QCRTRN= SQRT(chipot/ssvpot)
c
c calculate standard error of estimate
c
      SIGMA = SQRT(chipot/FLOAT(NESP))
c
c now write all this stuff out
c
      WRITE(6,1040) ssvpot
      WRITE(7,1040) ssvpot
      WRITE(6,1050) chipot 
      WRITE(7,1050) chipot 
      WRITE(6,1080) SIGMA
      WRITE(7,1080) SIGMA
      WRITE(6,1090) QCRTRN
      WRITE(7,1090) QCRTRN
C
C     ----- PRINT THE DIPOLE , QUADRUPOLE AND CENTER OF MASS ----
C
      if(nmol.eq.1) then
         WRITE(6,1340)
         write(6,1350) ( CMAS(I), I=1,3)
         WRITE(6,1360)
         write(6,1370) ( DIPOL(I), I=1,3)
         write(6,1375) dipmom
         write(7,1375) dipmom
         WRITE(6,1380)
         write(6,1390) ( QUAD(I), I=1,3)
         write(6,1395) ( QUAD(I), I=4,6)
      endif
c
 1000 FORMAT(/ /,10X,'Point charges after optimization',/ /,4X,
     *'NO',15X,'X',15X,'Y',15X,'Z',10X, 'Qoptim',/)
 1010 FORMAT(2X,I4,2X,4(5X,F10.6))
 1040 FORMAT(/,8X,'Statistics of the fitting:',
     *       /,2x,'The initial sum of squares (ssvpot)',t50,f15.3)
 1050 FORMAT(2x,'The residual sum of squares (chipot)',t50,f15.3)
 1080 FORMAT(2x,
     *'The std err of estimate (sqrt(chipot/N))',t50,f15.5)
 1090 FORMAT(2x, 'ESP relative RMS (SQRT(chipot/ssvpot))',t50,f15.5)
 1110 FORMAT(/,10A8)
 1130 FORMAT(/,t2,'Net charge on the system =',F10.7)
 1200 FORMAT(/,10X,
     *    'Point Charges Before & After Optimization',/ /,4X,'no.',
     *    2x,'At.no.',4x,'q(init)',7X,'q(opt)     ivary    d(rstr)/dq')
 1210 format(t2,2i4, 2(5x,f10.6), i7, f15.6, f12.3)
 1340 FORMAT(/,t2,'Center of Mass (Angst.):',/)
 1350 FORMAT(t2,'X  =',F10.5,5X,' Y  =',F10.5,5X,' Z  =',F10.5)
 1360 FORMAT(/,t2,'Dipole (Debye):',/)
 1370 FORMAT(t2,'X  =',F10.5,5X,' Y  =',F10.5,5X,' Z  =',F10.5)
 1375 FORMAT(/,t2,'Dipole Moment (Debye)=',F10.5)
 1380 FORMAT(/,t2,'Quadrupole (Debye*Angst.):',/)
 1390 FORMAT(t2,'Qxx =',F10.5,5X,'QYY =',F10.5,5X,'QZZ =',F10.5)
 1395 FORMAT(t2,'Qxy =',F10.5,5X,'QXZ =',F10.5,5X,'QYZ =',F10.5)
c
      return
      end
c-------------------------------------------------------------------------
      SUBROUTINE wrt_pot
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995,1997                 **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
c
c called from Main
c
c read in the electrostatic potential points used in the fitting,
c calculate esp using existing charges, and write out both esp's & residual
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, iuniq,NESP,natpl1, Ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10), keywd( 4,4)
      character*8   TITLE,     keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
      common /files/ input,output,qin,qout,punch,espot,qwts,esout,
     .               owrite
      character*80 input,output,qin,qout,punch,espot,qwts,esout
      character owrite
c
      DATA AU2CAL/627.5095d0/,BOHR/0.52917725d0/
c
c open the file containing the qm esp points & read in the no. of points
c
      call amopen(10, ESPOT,'O','F','R')
      rewind(10)
      call amopen(20, ESOUT,owrite,'F','W')
      read(10,'(2i5)') idum,nesp
      ssvkcl= ssvpot*au2cal*au2cal
      chikcl= chipot*au2cal*au2cal
      write( 20, '(i5,i6,4x,2f20.10)') nesp, izan(1), ssvkcl, chikcl
c
c build up matrix elements Ajk according to (SUMi 1/Rik SUMj 1/Rij)
c
      do i = 1,idum
         read(10,100, end=930)  xi, yi, zi
 100     format(17x,3e16.7)
      enddo
      do i= 1,nesp
         read(10,101, end=930) espqmi, xi, yi, zi
 101     format(1x,4e16.7)
         espclc= 0.0d0
         do k=1,iuniq
           Xik    = xi - CRD(1,k)
           Yik    = yi - CRD(2,k)
           Zik    = zi - CRD(3,k)
           espclc= espclc + qcal(k)/ SQRT( Xik*Xik + Yik*Yik + Zik*Zik)
         enddo
         vresid= espqmi - espclc
         xa= xi*BOHR
         ya= yi*BOHR
         za= zi*BOHR
         espclc= espclc*AU2CAL
         espqmi= espqmi*AU2CAL
         vresid= vresid*AU2CAL
         write( 20, '(3f10.5, 3f12.5)') xa,ya,za, espqmi,espclc,vresid
      enddo
c
      close(10)
      close(20)
c
      return
  930 write (6, *) ' unexpected eof in ', espot
      stop
      end
C-------------------------------------------------------------------
      SUBROUTINE charge_opt
c
c  driver for the charge determinization/optimizaton
c
c  called from Main
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "limits.h"
c      parameter (maxq   = 8000)
c      parameter (maxlgr = 800)
c      parameter (maxmol = 300)
c
      COMMON/IOSTUF/inopt,ioutopt,IQOPT,iunits
      COMMON/INFOA/NAT, IUNIQ,NESP,natpl1, ihfree,irstrnt
      COMMON/RUNLAB/TITLE(10),keywd( 4,4)
      character*8   TITLE,    keywd
      COMMON/ESPCOM/apot(maxq,maxq), bpot(maxq), grad(maxq),
     &               awt(maxq,maxq), bwt(maxq), ssvpot,chipot,vavrg
      COMMON/CALCUL/ QCAL(maxq), a(maxq,maxq), b(maxq),
     & qwtval(maxq), iwttyp(maxq), iqcntr(maxq)
      common/LAGRNG/ grpchg(maxlgr), lgrcnt(maxlgr,maxq), nlgrng
      COMMON/ORIG/q0(maxq),CRD(3,maxq),IVARY(maxq),IZAN(maxq),qwt,q0tot
      COMMON/worker/awork(maxq,maxq),bwork(maxq),scr1(maxq),iscr1(maxq)
      COMMON/propty/ CO(3,maxq), CMAS(3), DIPOL(3), dipmom, QUAD(6)
      COMMON/mltmol/ wtmol(maxmol), moleqv(4,maxmol), ibeg(maxmol),
     &               iend(maxmol), nmol
      dimension qold(maxq)
      integer istat
c
      gs= 0.0d0
      irsave= 0
      nitern= 0
      do i= 1,maxq
         qold(i)= 0.0d0
      enddo
c
c qtol & maxit are criteria for convergence & maximum iterations for
c the non-linear optimizations.
c
      qtol= 0.1d-5
      maxit=  24
c
c only on first pass through this subroutine (indicated by nitern= 0),
c if irstrnt > 0, transfer irstrnt to irsave and reset irstrnt to 0,
c in order to get an initial guess using a harmonic constraint.  This is
c done so restraint subroutine (rstran) will use a harmonic restraint.
c
      if( irstrnt .gt. 0 ) then
        irsave= irstrnt
        irstrnt= 0
        WRITE(6, '(/,t2,a)')'Non-linear optimization requested.'
      endif
c
c now go do a "harmonic restraint" run, restraint= qwt(qcal(i)-q0(i))**2
c -- loop to convergence
c
  100 continue
         call matbld
c
c        solve (Ax = b) where A and b are input, x is solution
c                     awork x = bwork
c
c        the solution "x" is returned in "b" (bwork)
c
c        -- condition the matrix diagonal to avoid DGETRF() detecting
c           singularity
c
         do jn = 1,natpl1
           if (abs(awork(jn,jn)).lt. 1.d-10) awork(jn,jn) = 1.d-10
         enddo
c
         call DGETRF( NATPL1,NATPL1,awork,maxq,iscr1,istat)
         if(istat.lt. 0) then
           write(6,'('' chgopt: LU decomp has an illegal value'')')
           stop
         elseif(istat.gt. 0) then
           write(6,'('' chgopt: LU decomp gave almost-singular U'')')
           stop
         endif
c
         call DGETRS( 'N',NATPL1,1,awork,maxq,iscr1,bwork,NATPL1,istat)
         if(istat.ne. 0) then
           write(6,'('' chgopt: LU backsubst has an illegal value'')')
           stop
         endif
c
c        -- copy solution vector "bwork" to 'calculated charges' vector 
c           qcal
c
         do k= 1, iuniq
           icntr= iqcntr(k)
           if( icntr .ge. 1 ) then
c
c            -- new charge
c
             qcal(k) = bwork( icntr)
           else
c
c            -- frozen charge
c
             qcal(k) = q0(k)
           endif
         enddo
c
c        -- a quick check from rstrn: if irstrnt is now negative, 
c           there are no restraints because no qwtval(i) > 0.1e-10, 
c           so reset irsave= 0
c
         if (irstrnt .lt. 0) then
           irsave=0
           WRITE(6,'(/,t2,a,/,t2,a)')
     .        'WARNING: Restraints were requested, but',
     .        '         the restraint weights were all zero'
         endif
c
c        -- we're finished if it's only a "harmonic restraint" run, 
c           but if it's a non-linear optimization (irsave>0)... 
c           we've only just begun (i.e. we have our initial guess)
c
         if( irsave .le. 0 ) then
           return
         else
c
c          -- it's a non-linear optimization: reset irstrnt (to now 
c             calculate the proper non-linear restraint derivatives 
c             in routine rstran)
c
           irstrnt= irsave
         endif
c
c        -- begin iterative optimization loop with comparison of 
c           old & new charges; calculate the convergence and replace 
c           the old charges with the new
c
         qchnge= 0.0d0
         xuniq = dble(iuniq)
         do i= 1,iuniq
            qdiff = qcal(i) - qold(i)
            qchnge = qchnge + (qdiff*qdiff)
            qold(i) = qcal(i)
         enddo
c
c        -- hp compiler bug requires that xuniq be calcd
c           before loop 11/94
c
         qchnge = sqrt(qchnge) / xuniq
         write(6,'(t2,''qchnge ='',g20.10)')qchnge
c
c        -- if this is less than qtol then we're done
c
         if (qchnge .lt. qtol .and. nitern .gt. 1) then
           write(6,'(/t2,''Convergence in'',i5,'' iterations'')')nitern
           return
         elseif( nitern .ge. maxit ) then
c
           write(6, 
     $        '(t2,''after '',i5,'' iterations, no convergence!'')') 
     $        maxit
           return
         endif
c
c loop again
c
         nitern= nitern + 1
      go to 100
      end
C-----------------------------------------------------------------------
      subroutine amopen(lun,fname,fstat,fform,facc)

c     INPUT:
c
      integer lun
c        ... logical unit number
      character*(*) fname
c        ... file name (not used in VAX/VMS implementation)
      character*1 fstat
c        ... status code: "N", "O", or "U" = new, old, unk.
      character*1 fform
c        ... format code: "U", "F" = unform., form.
      character*1 facc
c        ... access code: "R", "W", "A" = read, read/write, append
c
c     THIS IS UNIX VERSION
c     Author: George Seibel
c     Rev 13-Jun-90:  add rewind after open.
 
c     INTERNAL:
 
      character*7 stat
c        ... status keyword
      character*11 kform
c        ... form keyword
      integer ios
c        ... i/o status variable
 
      if (fstat .eq. 'N') then
           stat = 'NEW'
      elseif (fstat .eq. 'O') then
           stat = 'OLD'
      elseif (fstat .eq. 'U') then
           stat = 'UNKNOWN'
      else
           write(6,'(/,2x,a,i4)')
     $           'amopen: bogus fstat, unit ', lun
           stop
      endif
c
      if (fform .eq. 'U') then
           kform = 'UNFORMATTED'
      elseif (fform .eq. 'F') then
           kform = 'FORMATTED'
      else
           write(6,'(/,2x,a,i4)')
     $           'amopen: bogus fform, unit', lun
           stop
      endif
c
      open(unit=lun,file=fname,status=stat,form=kform,iostat=ios)
c
      if (ios .ne. 0) then
           if (lun .eq. 6) then
                write(0,'(/,2x,a,i4,a,a)') 'Unit ', lun, 
     $                                    ' Error on OPEN: ',fname
                close(unit=0)
           else
                write(6,'(/,2x,a,i4,a,a)') 'Unit ', lun, 
     $                                    ' Error on OPEN: ',fname
                close(unit=6)
                write(0,'(/,2x,a,i4,a,a)') 'Unit ', lun, 
     $                                    ' Error on OPEN: ',fname
                close(unit=0)
           endif
           stop
      endif
      rewind(lun)
      return
      end
