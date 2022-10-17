      FUNCTION ISAMAX (N, X, INC)
C***  Find index of vector element with highest absolute value

ccc   The following call works on intel compiler, but not with gfortran
ccc      ISAMAX = IDAMAX (N, X, INC)

C***   Hand-made version
      DIMENSION X(N)

      XMAX = ABS(X(1))
      IMAX = 1

      DO I=2*INC, N, INC
        IF (ABS(X(I)) .GT. XMAX) THEN
          XMAX = ABS(X(I))
          IMAX = I
        ENDIF
      ENDDO

      ISAMAX = IMAX

      RETURN
      END
