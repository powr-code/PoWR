      FUNCTION MY_ISAMAX (N, X, INC)

      MY_ISAMAX = ISAMAX (N, X, INC)

c      DIMENSION X(N)

c      XMAX = ABS(X(1))
c      IMAX = 1

c      DO I=2*INC, N, INC
c        IF (ABS(X(I)) .GT. XMAX) THEN
c          XMAX = ABS(X(I))
c          IMAX = I
c        ENDIF
c      ENDDO

c      ISAMAX = IMAX

      RETURN
      END
