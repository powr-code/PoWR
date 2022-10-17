      FUNCTION TRADFUN (XLAMBDA,XJ)
C***********************************************************************
C***  RADIATION TEMPERAURE IN KELVIN, FROM XJ = J-NUE (CGS UNITS)
C***  AND XLAMBDA IN ANGSTROEM
C***  CONSTANTS :  C1 = H * C / K   (DIMENSION ANGSTROEM * KELVIN )
C***               C2 = 2 * H * C
C***  Version improved to prevent overflow for almost-zero XJ
C***        wrh 14-Mar-2005 17:15:26
C***********************************************************************

      DATA C1,C2 / 1.4388E8, 3.9724E+8 /
C***  Threshold for Taylor expansion of ALOG(1+X)
      DATA EPS / 1.E-14 /

C***  Zero Trad for negative J
cc      IF (XJ .LE. .0) THEN 
cc wegen Absturz:
      IF (XJ .LE. 1.E-100) THEN 
         TRADFUN=.0
         RETURN
      ENDIF

      W=1./XLAMBDA
      W3 = W * W * W

      IF (C2 * W3 / EPS .LT. XJ) THEN 
C***  Equivalent to:      IF (X .LT. 1.E-14) THEN
C***     SERIES EXPANSION: LN(1+X) = X
         TRADFUN = C1 * XLAMBDA * XLAMBDA * XJ / C2
      ELSE
         X = C2 * W3 / XJ
         TRADFUN=C1*W/ALOG(1.+X)
      ENDIF

      RETURN
      END
