
      FUNCTION FIERFC(X)
C**********************************************************************
C***  FIRST REPEATED INTEGRAL OF THE ERROR FUNCTION (COMPLEMENT)
C***  - RATIONAL APPROXIMATION FROM ABRAMOWITZ AND STEGUN, P. 299
C***  - NOT VALID FOR 
C**********************************************************************
      DIMENSION A(3)
      DATA A / 0.3480242, -0.0958798, 0.7478556 /
      DATA WPIINV / 0.564189583549 /

      IF ( X .LT. 0. ) THEN
        WRITE (0,*) 'FUNCTION FIERFC CALLED WITH NEGATIVE ARGUMENT'
        STOP 'ERROR'
      ENDIF
      T = 1. / (1. + 0.47047 * X)
      POLY = ((A(3) * T + A(2)) * T + A(1)) * T
      FIERFC = (WPIINV - X*POLY) * EXP(-X*X)

      RETURN
      END
