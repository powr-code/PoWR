      SUBROUTINE XEXTINC(EXTINC, XEV, COLDENS)
C*********************************************************************
C*** SEE MORRISON AND McCAMMON APJ 1983 SIDE 119-122
C***
C***  CALCULATING EXTINCTION IN X-RAY 0.03-10.KEV 
C*********************************************************************

      DIMENSION XKEVB(15)
      DIMENSION C0(14), C1(14), C2(14)

C*** ENERGY-INTERVALS     
      DATA XKEVB / 0.03 , 0.1  , 0.284, 0.4 ,  0.532, 
     $             0.707, 0.867, 1.303, 1.84,  2.471,
     $             3.21 , 4.038, 7.111, 8.331, 10.    /


C*** C0,C1,C2 COEFFICIENTS OF ANALYTIC FIT TO CROSSECTION
      DATA C0/17.3,34.6,78.1,71.4,95.5,308.9,120.6,141.3,202.7,342.7,
     $        352.2,433.9,629.,701.2/
      DATA C1/608.1,267.9,18.8,66.8,145.8,-380.6,169.3,146.8,104.7,
     $        18.7,18.7,-2.4,30.9,25.2/
      DATA C2/-2150.,-476.1,4.3,-51.4,-61.1,294.,-47.7,-31.5,-17.,.0,
     $        .0,0.75,.0,.0/
  
      NXKEVB=15
      XKEV = XEV * 0.001

C*** CROSSECTION XSIG PER HYDROGEN ATOM IN CM^2
C*** 1.36248E-2 KEV IS LYMAN-EDGE
      IF (XKEV .LT. 1.36248E-2) THEN
         XSIG = .0
         GOTO 10

      ELSE IF (XKEV .LT. XKEVB(1)) THEN
         XSIG = 1.08E-17 * (1.36248E-2 / XKEV )**2.7
         GOTO 10

      ELSE IF (XKEV .LT. XKEVB(NXKEVB)) THEN
        DO 1 K=1,NXKEVB-1
        IF ((XKEV .GE. XKEVB(K)) .AND. (XKEV .LT. XKEVB(K+1))) THEN
          XSIG = (((C2(K)*XKEV+C1(K))*XKEV)+C0(K))*XKEV**(-3.)*1.E-24
          GOTO 10
        ENDIF
    1   CONTINUE

      ELSE 
         XSIG = 0.
      ENDIF

   10 CONTINUE    

C***  EXTINCTION = EXP(-TAU)
C***  TAU = SIGMA * COLUMN_DENSITY
      EXTINC = EXP(-XSIG * COLDENS)

      RETURN
      END
