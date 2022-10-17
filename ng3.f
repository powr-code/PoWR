      SUBROUTINE NG3 (ND,N,RNE,NCHARG,POPNEW,POP,POP1,POP2,NATOM,ABXYZ,
     $         NFIRST,NLAST,RADIUS,W,TNEW,T,TOLD1,TOLD2,NOTEMP,NGWEIGHT,
     $         NDONE, NTDONE)
C***********************************************************************
C***  NG'S ACCELERATION METHOD (K.C.NG 1974, CHEM. PHYS. 61, 2680):
C***  EXTRAPOLATION BY CONSIDERING THE LAST THREE ITERATIONS
C***                      (NG'S PAPER: LAST FOUR ITERATIONS!!!)
C***  === VERSION ===: EXTRAPOLATION OF LOGARITHMIC POP. NUMBERS
C***                   IF (.NOT.NOTEMP): ALSO EXTRAPOLATION OF TEMPERATURE
C***  "INTEGRATION" WEIGHTS:  
C***    LOGICAL SWITCH: NGWEIGHT
C***      .TRUE. = WEIGHTED WITH 1/VARIABLE (OPTION: NGWEIGHT)
C***     .FALSE. = UNWEIGHTED               (DEFAULT)
C***
C***  TESTED : INCREASING CORRECTIONS IN GENERAL (A>1.)
C***             ==> NO EXTRAP DONE AT ALL
C***           INCREASING CORRECTIONS ONLY FOR SOME POPNUMBERS
C***             ==> NO EXTRAP DONE FOR THIS POPNUMBER
C***********************************************************************
 
      DIMENSION POPNEW(ND,N),POP(ND,N),POP1(ND,N),POP2(ND,N)
      DIMENSION RNE(ND),NCHARG(ND),RADIUS(ND),W(ND)
      DIMENSION TNEW(ND),T(ND),TOLD1(ND),TOLD2(ND)
      DIMENSION ABXYZ(NATOM),NFIRST(NATOM),NLAST(NATOM)
      LOGICAL NOTEMP, NGWEIGHT
      COMMON / COMNEGT / NEGT

C***********************************************************************
C***  SET MINIMUM TEMPERATURE FOR EXTRAPOLATION OF T(R)
      TMIN=3000.
C***********************************************************************

C***  INITIALIZE COUNTER FOR NEGATIVE TEMPERATURE WARNING
      NEGT=0
      NDONE = 0
      NTDONE = 0

      DO 101 J=1,N
      DO 102 L=1,ND
        IF ( POP(L,J) .LT. 1.E-99 ) POP(L,J)=1.E-99
        IF ( POP1(L,J) .LT. 1.E-99 ) POP1(L,J)=1.E-99
        IF ( POP2(L,J) .LT. 1.E-99 ) POP2(L,J)=1.E-99
102   CONTINUE
101   CONTINUE

C***  INITIALIZATION OF "INTEGRATION" WEIGHTS (UNWEIGHTED)
      DO 2 L=1, ND
      W(L)=1.
    2 CONTINUE


***  LOOP FOR EACH LEVEL  ---------------------------------------------
      DO 10 J=1,N
      A1=0.0
      C1=0.0

C***  INITIALIZATION OF "INTEGRATION" WEIGHTS
      IF (NGWEIGHT) THEN
         DO 1 L=1,ND
         IF (POP(L,J) .GT. 0.) THEN
            W(L)=1./POP(L,J)
            ELSE
            W(L) = 1.
            ENDIF
    1    CONTINUE
         ENDIF

C***  CALCULATION OF THE EXTRAPOLATION CONSTANTS FOR NG'S ACCELERATION METHOD
      DO 5 L=1,ND
      WL=W(L)
      D01=ALOG10(POP(L,J)*POP2(L,J)/POP1(L,J)/POP1(L,J))
      DN=ALOG10(POP(L,J)/POP1(L,J))
      A1=A1+D01*D01*WL
      C1=C1+D01*DN*WL
    5 CONTINUE
      IF (A1 .EQ. 0.) THEN
         A=0.0
         ELSE
         A=C1/A1
         ENDIF
      DO 12 L=1,ND
        POPNEW(L,J) = POP(L,J)
   12 CONTINUE
C***  CONTROL IN GENERAL
      IF (A .LT. 1.) THEN
        DO 8 L=1,ND
C***      CONTROL FOR EACH POPNUMBER
          IF (ABS(POP2(L,J)-POP1(L,J)) .GE. ABS(POP1(L,J)-POP(L,J)))
     >    THEN
            NDONE = NDONE + 1
            POPNEW(L,J) = POP(L,J)*(POP1(L,J)/POP(L,J))**A
C***        CORRECTION NOT REATER THAN 10 TIMES OF LAST CORRECTION
            IF (ABS(POPNEW(L,J)-POP(L,J)) .GT. 
     >          10. *ABS(POP(L,J)-POP1(L,J))) THEN
              POPNEW(L,J) = POP(L,J) + 10.*(POP(L,J)-POP1(L,J))
            ENDIF
          ENDIF
    8   CONTINUE
      ENDIF
   10 CONTINUE
C***  ENDLOOP  ---------------------------------------------------------
 
C***  POPNUMBERS ARE SCALED TO ENSURE NUMBER CONSERVATION
C***  ELECTRON DENSITY IS UPDATED
      DO 20 L=1,ND
      RNE(L)=0.0
      DO 20 NA=1,NATOM
      NFIRNA=NFIRST(NA)
      NLANA=NLAST(NA)
      POPSUM=0.0
      DO 15 J=NFIRNA,NLANA
      POPSUM=POPSUM+POPNEW(L,J)
   15 CONTINUE
      POPSUM=POPSUM/ABXYZ(NA)
      DO 17 J=NFIRNA,NLANA
      POPNEW(L,J)=POPNEW(L,J)/POPSUM
      RNE(L)=RNE(L)+NCHARG(J)*POPNEW(L,J)
   17 CONTINUE
   20 CONTINUE
 
C***  IN CASE OF TEMPERATURE CORRECTIONS: ALSO EXTRAPOLATION OF T(R)
      IF (.NOT. NOTEMP) THEN

C***     INITIALIZATION OF "INTEGRATION" WEIGHTS
         DO 11 L=1,ND
           IF (NGWEIGHT) THEN
             W(L)=1. / T(L)
           ELSE
             W(L) = 1.
           ENDIF
   11    CONTINUE
         A1=0.0
         C1=0.0
C***     CALCULATION OF THE EXTRAPOLATION CONSTANTS
         DO 25 L=1,ND
         WL=W(L)
         DN=T(L)-TOLD1(L)
         D01=DN-(TOLD1(L)-TOLD2(L))
         A1=A1+D01*D01*WL
         C1=C1+D01*DN*WL
   25    CONTINUE
      IF (A1 .EQ. 0.) THEN
         A=0.0
         ELSE
         A=C1/A1
         ENDIF
C!!!      WRITE(*,*) 'TESTOUTPUT TEMPERATUR-EXTRAP'
C!!!      WRITE(*,*) 'A=',A
C***  NO TEMPERATURE EXTRAPOLATION IF CORRECTION INCREASES IN GENERAL (A >= 1.)
         IF (A .LT. 1.) THEN
           DO 28 L=1, ND
C***  CONTROL FOR EACH DEPTH=POINT
C!!!      WRITE(*,*) 'L=',L,' T=',T(L),
C!!!     >           ' TOLD1=',TOLD1(L),' TOLD2=',TOLD2(L)
             IF (ABS(TOLD1(L)-TOLD2(L)) .GE. ABS(T(L)-TOLD1(L))) THEN
C!!!      WRITE(*,*) 'DONE'
               NTDONE = NTDONE + 1
               TNEW(L) = (1.-A) * T(L) + A * TOLD1(L)
C***  CORRECTION NOT GREATER THAN 10 TIMES OF THE LAST CORRECTION
               IF (ABS(TNEW(L) - T(L)) .GT. 10.*ABS(T(L) - TOLD1(L)))
     >           THEN
C!!!      WRITE(*,*) 'CORRECTION BRAKED'
                 TNEW(L) = T(L) + 10. * (T(L) - TOLD1(L))
               ENDIF
               IF (TNEW(L) .LT. TMIN) THEN
                 TNEW(L) = TMIN
                 NEGT = NEGT + 1
               ENDIF
             ELSE
               TNEW(L) = T(L)
             ENDIF
   28     CONTINUE
        ELSE
          DO 29 L=1, ND
            TNEW(L) = T(L)
   29     CONTINUE
        ENDIF
      ELSE
C***  NO TEMPERATURE EXTRAPOLATION
        DO 38 L=1, ND
          TNEW(L) = T(L)
   38   CONTINUE
      ENDIF

      RETURN
      END
