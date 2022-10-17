C******************************************************************************
      SUBROUTINE SFIT (ND, OPA, ETA, ENTOT, XLAM, SBOUND, PLOT)
C******************************************************************************
C***  THIS SUBROUTINE PREPARES SBOUND = SOURCE FUNCTION AT OUTER BOUNDARY,
C***  OBTAINED BY A SMOOTH EXTRAPOLATION. AT THE OUTERMOST MBOUND
C***  DEPTH POINTS, THE LINE SOURCE FUNCTION IS CONVERTED INTO RADIATION
C***  TEMPERATURES, WHICH ARE THEN FITTED BY A WEIGHTED LEAST-SQUARE FIT
C***  POLYNOMIAL. 
C***  INDEPENDETN VATIABLE: LOG ENTOT(L) - LOG ENTOT(L=1)
C***  WEIGHTS: AS APPROPRIATE FOR AN INTEGRATION (TRAPEZOIDAL RULE). THUS 
C***  THE CROWDED POINTS CLOSE TO THE BOUNDARY ARE NOT HEAVILY WEIGHTED.
C***
C***  THE APPLIED SUBROUTINES E02ADF AND E02AEF ARE FROM THE NAG FORTRAN LIBRARY
C***  THE FOLLOWING PARAMETERS MUST BE SUITABLE CHOOSEN:
C***     MBOUND = NUMBER OF DEPTH POINTS ACCOUNTED FOR IN THE EXTRAPOLATION
C***     KPLUS1 = DEGREE OF FIT POLYNOMIAL + 1 (E.G. QUADRATIC ... KPLUS1=3)
C******************************************************************************

      PARAMETER ( MBOUND = 30 )
      PARAMETER ( KPLUS1 =  4 )
C***  Note: Other degrees are not supported by SUBR. POLYFIT !

      DIMENSION OPA(ND), ETA(ND), ENTOT(ND)
      DIMENSION XFIT(MBOUND),YFIT(MBOUND),WFIT(MBOUND),SDEV(KPLUS1)
      DIMENSION AFIT(KPLUS1,KPLUS1), WORK1(3*MBOUND),WORK2(2*KPLUS1)
      DIMENSION A(KPLUS1,KPLUS1), B(KPLUS1)
      LOGICAL PLOT

      DIMENSION ATEST(KPLUS1,KPLUS1), BTEST(KPLUS1), DTEST(KPLUS1)
      DIMENSION SCRATCH(2*KPLUS1), XKOEFF(KPLUS1)

C***  Operating system:                    
      COMMON / COMOS / OPSYS
      CHARACTER*8 OPSYS

      IF (MBOUND .GT. ND) THEN
         WRITE (0,*) 'MBOUND EXCEEDS ND'
         WRITE (0,*) 'MBOUND, ND=', MBOUND, ND
         PRINT *, 'MBOUND EXCEEDS ND'
         STOP 'ERROR'
         ENDIF

C***  DEFINE THE TABLE OF (X,Y)-VALUES
      A1=ALOG10(ENTOT(1))
      DO 10 L=1,MBOUND
      XFIT(L)=ALOG10(ENTOT(L))-A1
      IF (OPA(L) .LE. .0) THEN
         YFIT(L)=.0
         ELSE
         S=ETA(L)/OPA(L)
         TRAD=TRADFUN(XLAM,S)
         YFIT(L)=TRAD
         ENDIF
   10    CONTINUE

C***  DEFINE WEIGHTS
      WFIT(1)=XFIT(2)-XFIT(1)
      DO 11 L=2,MBOUND-1
      WFIT(L)=XFIT(L+1)-XFIT(L-1)
   11 CONTINUE
      WFIT(MBOUND)=XFIT(MBOUND)-XFIT(MBOUND-1)

C***  Former CRAY branch deleted wrh  1-Jun-2021
        CALL POLYFIT (XFIT, YFIT, WFIT, MBOUND, KPLUS1,
     >                XKOEFF, A, B, SCRATCH, ATEST, BTEST, DTEST)

C***   As the X scale is defined such that X=0 for RMAX, 
C***      the value of the polynomial is simply the last coefficient
ccc        CALL HORNER (XFIT(1), TRADBOUND, KPLUS1, XKOEFF)
        TRADBOUND = XKOEFF(KPLUS1)

C***  If the fit polynomial has a minimum inside the interval, 
C***     the minimum value is taken instead of the value at x=0
C***  Note: This branch is also restricted to KPLUS1=4 !!
cccC***    First test: negative derivative at x=0
ccc        IF (XKOEFF(3) .GT. .0) GOTO 16
        IF (XKOEFF(KPLUS1) .LE. 0.) GOTO 16
ccc     Fix added: Skip for XKOEFF(1) == 0    (ansander, 06-Oct-2016)
        IF (ABS(XKOEFF(1)) <= 1.E-20) GOTO 16
        PHALB =    XKOEFF(2) / (3.*XKOEFF(1))
        Q     =    XKOEFF(3) / (3.*XKOEFF(1))
C***    Second test: Existence of (real) solution for f'=0
        WURZ = PHALB*PHALB - Q
        IF (WURZ .LT. .0) GOTO 16 
        WURZ = SQRT(WURZ)

C***    First solution of quadratic equation
        X1 = -PHALB + WURZ
C***    Test if X1 inside interval 0, X(MBOUND)
        IF (X1 .LT. .0) GOTO 15
        IF (X1 .GT. XFIT(MBOUND)) GOTO 15
C***    Test, if minimum, not maximum (f" > 0)
        F2STRICH = 2. * XKOEFF(2) + 6. * XKOEFF(1) * X1  
        IF (F2STRICH .LT. .0) GOTO 15
C***    X1 is minimum!        
        CALL HORNER (X1, TRADMIN, KPLUS1, XKOEFF)
        IF (TRADMIN .LT. TRADBOUND) TRADBOUND = TRADMIN

   15   CONTINUE
C***    Now the same for second solution of quadratic equation
        X1 = -phalb - WURZ
C***    Test if X1 inside interval 0, X(MBOUND)
        IF (X1 .LT. .0) GOTO 16
        IF (X1 .GT. XFIT(MBOUND)) GOTO 16
C***    Test, if minimum, not maximum (f" > 0)
        F2STRICH = 2. * XKOEFF(2) + 6. * XKOEFF(1) * X1  
        IF (F2STRICH .LT. .0) GOTO 16
C***    X1 is minimum!        
        CALL HORNER (X1, TRADMIN, KPLUS1, XKOEFF)
        IF (TRADMIN .LT. TRADBOUND) TRADBOUND = TRADMIN

   16   CONTINUE

      SBOUND=BNUE(XLAM,TRADBOUND)

C***  TEST PLOT, IF REQUESTED:
      IF (PLOT) CALL SFITPLO (ENTOT,ETA,OPA,ND,MBOUND,AFIT,
     >                        KPLUS1,XLAM,XKOEFF, TRADBOUND)

      RETURN
      END
