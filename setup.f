      SUBROUTINE SETUP (L,A,B,C,W,JMAX,ND,NP,NPDIM,OPA,ETA,THOMSON,
     $          XLAM,Z,RADIUS,BCORE,DBDR,XIMINUS,ENTOT,K, IVERS)
C***********************************************************************
C***  ANGLE-DEPENDENT CONTINUOUS RADIATION TRANSFER IN SPHERICAL SYMMETRY: 
C***  SET UP THE COEFFICINT MATRICES FOR A FEAUTRIER SCHEME:
C***     A (DIAGONAL), B (FULL), C (DIAGONAL) AND W (VECTOR)
C***  CALLED FROM: SUBROUTINE ELIMIN
C***********************************************************************

      DIMENSION A(NPDIM),B(NPDIM,NPDIM),C(NPDIM),W(NPDIM)
      DIMENSION RADIUS(ND),OPA(ND),ETA(ND),THOMSON(ND)
      DIMENSION Z(ND,NP)
      LOGICAL PLOT

      JMAX=NP+1-L
      JMM=JMAX-1
 
C***  EVERY L = 1 ... ND
      X=OPA(L)
      G=-X*THOMSON(L)
      ETAL=ETA(L)
 
C***  MEAN INTENSITY INTEGRATION WEIGHTS FROM SUBROUTINE MOMENT0 (VEKTOR W)
      CALL MOMENT0 (ND,RADIUS,L,JMAX,Z,W,DUMMY,.TRUE.)
      DO 1 J=1,JMAX
      WJG=W(J)*G
      DO 1 JS=1,JMAX
    1 B(JS,J)=WJG
      DO 3 J=1,JMAX
    3 W(J)=ETAL
 
      IF(L.EQ.1) GOTO 9
      IF(L.EQ.ND) GOTO 10
 
C***  ALL NON-BOUNDARY POINTS  L= 2 ... ND-1
      XP=(X+OPA(L+1))/2.
      XM=(X+OPA(L-1))/2.
      DO 2 J=1,JMM
      ZLPLUS=Z(L+1,J)
      ZLJ=Z(L,J)
      ZLMIN=Z(L-1,J)
      DT=2./(ZLMIN-ZLPLUS)
      DTM=XM*(ZLMIN-ZLJ)
      DTP=XP*(ZLJ-ZLPLUS)
      A(J)=DT/DTM
      C(J)=DT/DTP
    2 B(J,J)=B(J,J)+A(J)+C(J)+X
 
C     LAST ROW OF BLOCK, J=JMAX
      ZLMIN=Z(L-1,JMAX)
      DT=ZLMIN*XM
      A(JMAX)=2.*X/DT/DT
      B(JMAX,JMAX)=B(JMAX,JMAX)+A(JMAX)+X
      RETURN
 
C***  OUTER BOUNDARY CONDITION (AT L=1) 
    9 CONTINUE
      
C***  DIFFERENT VERSIONS FOR IMINUS
      IF (IVERS .EQ. 0) THEN
         XIMINUS=0.

      ELSE IF (IVERS .EQ. 1) THEN
         TAU = RADIUS(1) * X
         SBOUND = ETAL / X
         XIMINUS = SBOUND * (1.-EXP(-TAU))
   
      ELSE IF (IVERS .EQ. 4) THEN
C***     TESTPLOT FACILITY: DELETE "C" + SET "K" !!!
C***     PROGRAM             WRCONT: K = CONTINUUM INDEX
C***               ETL, CMF, FORMAL: K = LINE INDEX
         PLOT=.FALSE.
C         PLOT = K .EQ. 13
         CALL SFIT (ND, OPA, ETA, ENTOT, XLAM, SBOUND, PLOT)
         TAU = RADIUS(1) * X
         IF (TAU .LT. 1.E-3) THEN
            XIMINUS = SBOUND * TAU / 3.
         ELSE
            EXPTB = EXP(-TAU)
            FAK1 = 2./TAU
            FAK2 = FAK1/TAU
            XIMINUS = SBOUND * (1.-FAK1+FAK2-FAK2*EXPTB)
         ENDIF

      ELSE
         PRINT *,'INVALID VERSION OF OUTER BOUNDARY CONDITION'
         STOP 'ERROR'
         ENDIF

      XP=(X+OPA(2))/2.

      DO 8 J=1,JMM
      ZLPLUS=Z(2,J)
      ZLJ=Z(1,J)
      DT=XP*(ZLJ-ZLPLUS)
C***  MODIFICATION FOR NONZERO INCIDENT RADIATION  FROM TRUNCATED LAYERS
      W(J)=ETAL + XIMINUS * X *2./DT
      C(J)=2.*X/DT/DT
      B(JMAX,J)=.0
    8 B(J,J)=B(J,J)+C(J)+X+2.*X/DT
      B(JMAX,JMAX)=X
      W(JMAX)=XIMINUS * X
      RETURN
 
C***  INNER BOUNDARY CONDITION    L = ND
   10 XM=(X+OPA(ND-1))/2.
      DO 14 J=1,JMM
      ZLMIN=Z(ND-1,J)
      ZLJ=Z(ND,J)
      DT=XM*(ZLMIN-ZLJ)
      A(J)=2.*X/DT/DT
      B(J,J)=B(J,J)+A(J)+X+2.*X/DT
      PLUSI=BCORE+DBDR*ZLJ/X
   14 W(J)=ETAL+PLUSI*2.*X/DT
      A(JMAX)=2.*X/DT/DT
      B(JMAX,JMAX)=B(JMAX,JMAX)+A(JMAX)+X
      W(JMAX)=ETAL
      RETURN
      END
