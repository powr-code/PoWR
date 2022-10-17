      FUNCTION ERF(X)
C***********************************************************************
C***  DIESE ROUTINE IST NEU ABGELOCHT UND NOCH NICHT GETESTET
C***  INTEGRAL OVER NORMALIZED GAUSS FUNCTION FROM -INFINITY TO X
C***********************************************************************

      REAL A(4)
      DATA A / .25482 9592 , -.28449 6736,  1.42141 3741,-1.453152027 /
      Z=ABS(X)
      T=1./(1.+.3275911*Z)
      ERF=1.06140 5429
      DO 1 I=1,4
      K=5-I
    1 ERF=A(K)+T*ERF
      ERF=T*ERF*EXP(-Z*Z)*.5
      IF (X .GT. .0) ERF=1.-ERF

      RETURN
      END
