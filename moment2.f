      SUBROUTINE MOMENT2 (R,JMAX,P,U,XK)
C***********************************************************************
C***  INTEGRATION OF THE 2. MOMENT XK OF THE RADIATION FIELD U
C***  FEAUTRIER-INTENSITY U(J), IMPACT PARAMETER MESH P(J)
C***  AND RADIUS POINT R ARE GIVEN.
C***  WEIGHTS ARE ACCORDING TO TRAPEZOIDAL RULE IN Z*Z*DZ, Z=SQRT(R*R-P*P)
C***********************************************************************

      DIMENSION P(JMAX),U(JMAX)

c     DON'T change this back, fp consistency needed!
c      RR=R*R
C***  FIRST STEP, IMPLYING P(1)=0
      Z=R
      ZQ=R*R
      PJ=P(2)
      ZNQ=R*R-PJ*PJ
      ZNQ=MAX(ZNQ,0.)
      ZNEXT=SQRT(ZNQ)
      W=Z*(3*ZQ-ZNQ)-ZNEXT*(ZQ+ZNQ)
      XK=W*U(1)
C***  MIDDLE STEPS
      DO 1 J=3,JMAX
      ZLAST=Z
      ZLQ=ZQ
      Z=ZNEXT
      ZQ=ZNQ
      PJ=P(J)
      ZNQ=R*R-PJ*PJ
      ZNQ=MAX(ZNQ,0.)
      ZNEXT=SQRT(ZNQ)
      W=Z*(ZLQ-ZNQ)+ZLAST*(ZLQ+ZQ)-ZNEXT*(ZQ+ZNQ)
    1 XK=XK+W*U(J-1)
C***  LAST STEP, IMPLYING P(JMAX)=R
      W=Z*ZQ
      XK=XK+W*U(JMAX)
c     old:
c      XK=XK/R/RR/12.
      XK=XK/R/R/R/12.

      RETURN
      END
