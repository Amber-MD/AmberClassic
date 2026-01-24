      SUBROUTINE CORTOR (CORD1,CORD2,CORD3,BL,ALPHA,THETA,CORDO)
      DIMENSION CORD1(3),CORD2(3),CORD3(3),CORDO(3),VEC12(3),VEC23(3),
     &          FCORD(3),ECORD(3),UNSETF(3,3)
C----------------------------------------------------------------------C
C                                                                      C
C                                                 Z                    C
C                           1                     !      O   Y         C
C                           C :                   !     /  /           C
C                               :         ALPHA   !   / /              C
C                                 :               ! //                 C
C                                   2 ----------> 3    -------> X      C
C                                   B             A                    C
C                                                                      C
C----------------------------------------------------------------------C
      DATA DTOR /0.01745329/
C
      DO 10 J=1,3
      VEC23(J)=CORD3(J)-CORD2(J)
      VEC12(J)=CORD2(J)-CORD1(J)
   10 CONTINUE
      FCORD(1)=-BL*COS(ALPHA*DTOR)
      FCORD(2)= BL*SIN(ALPHA*DTOR)*COS(THETA*DTOR)
      FCORD(3)= BL*SIN(ALPHA*DTOR)*SIN(THETA*DTOR)
      CALL UNISET (2,VEC23,VEC12,UNSETF)
      CALL FTOE1 (UNSETF,FCORD,ECORD)
      DO 20 J=1,3
      CORDO(J)=CORD3(J)+ECORD(J)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE FTOE1 (UNSETF,FCORD,ECORD)
      DIMENSION UNSETF(3,3),FCORD(3),ECORD(3)
C----------------------------------------------------------------------C
C          THIS SUB. CHANGES BASES  ( E --> F )                        C
C          THAT IS TO SAY, A COORDINATE ON F IS CHANGED TO             C
C                          A COORDINATE ON E .                         C
C          E IS A SET OF THREE FUNDAMENTAL VECTORS (1,0,0),(0,1,0),    C
C                                          (0,0,1)                     C
C----------------------------------------------------------------------C
      DO 10 J=1,3
      ECORD(J)=UNSETF(J,1)*FCORD(1)+UNSETF(J,2)*FCORD(2)
     &        +UNSETF(J,3)*FCORD(3)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE UNISET (IND,VECA,VECB,UNSET)
C
      DIMENSION VECA(3),VECB(3),UNSET(3,3),
     &          VECX(3),VECY(3),VECZ(3),UECX(3),UECY(3),UECZ(3)
C
      IPRINT=0
      JND = IND
      IF (IND.LT.0) THEN
         IPRINT=1
         JND=-IND
      ENDIF
      IF (JND.EQ.0 .OR. JND.GE.8) GO TO 999
C
      GO TO (1,2,3,4,5,6,7),JND
C
    1 CALL VECPRD (VECA,VECB,VECZ)
      CALL VECPRD (VECZ,VECA,VECY)
      GO TO 91
    2 CALL VECPRD (VECB,VECA,VECZ)
      CALL VECPRD (VECZ,VECA,VECY)
      GO TO 91
    3 CALL VECPRD (VECB,VECA,VECY)
      CALL VECPRD (VECA,VECY,VECZ)
      GO TO 91
C
    4 CALL VECPRD (VECB,VECA,VECZ)
      CALL VECPRD (VECA,VECZ,VECX)
      GO TO 94
    5 CALL VECPRD (VECA,VECB,VECX)
      CALL VECPRD (VECX,VECA,VECZ)
      GO TO 94
C
    6 CALL VECPRD (VECA,VECB,VECY)
      CALL VECPRD (VECY,VECA,VECX)
      GO TO 96
    7 CALL VECPRD (VECB,VECA,VECX)
      CALL VECPRD (VECA,VECX,VECY)
      GO TO 96
C
   91 CALL UNIVEC (VECA,UECX)
      CALL UNIVEC (VECY,UECY)
      CALL UNIVEC (VECZ,UECZ)
      GO TO 100
   94 CALL UNIVEC (VECX,UECX)
      CALL UNIVEC (VECA,UECY)
      CALL UNIVEC (VECZ,UECZ)
      GO TO 100
   96 CALL UNIVEC (VECX,UECX)
      CALL UNIVEC (VECY,UECY)
      CALL UNIVEC (VECA,UECZ)
C
  100 DO 10 I=1,3
      UNSET(I,1)=UECX(I)
      UNSET(I,2)=UECY(I)
      UNSET(I,3)=UECZ(I)
   10 CONTINUE
      IF (IPRINT.EQ.0) RETURN
C
      WRITE (6,601) IND,((UNSET(J,I),I=1,3),J=1,3)
  601 FORMAT (/1H ,' *** UNSET WILL BE PRINTED BY IND=',I3,' AT UNISET'
     &            ,' ***'
     &       //1H ,4X,'1(X)',4X,'2(Y)',4X,'3(Z)'
     &        /1H ,'X',3F8.3/1H ,'Y',3F8.3/1H ,'Z',3F8.3/)
      RETURN
C
  999 WRITE(6,699) IND
  699 FORMAT(/1H ,' *** ERROR IS DETECTED AT UNISET : IND=',I4,' ***'/)
      RETURN
      END
      SUBROUTINE UNIVEC (VEC,UEC)
      DIMENSION VEC(3),UEC(3)
      ALENGH=SQRT(VEC(1)**2+VEC(2)**2+VEC(3)**2)
      DO 10 J=1,3
      UEC(J)=VEC(J)/ALENGH
   10 CONTINUE
      RETURN
      END
      SUBROUTINE VECPRD (VECA,VECB,VECP)
      DIMENSION VECA(3),VECB(3),VECP(3)
      VECP(1)=VECA(2)*VECB(3)-VECA(3)*VECB(2)
      VECP(2)=VECA(3)*VECB(1)-VECA(1)*VECB(3)
      VECP(3)=VECA(1)*VECB(2)-VECA(2)*VECB(1)
      RETURN
      END
