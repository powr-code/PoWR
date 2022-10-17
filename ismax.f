      FUNCTION ISMAX (N, X, INC)

      DIMENSION X(N)

      XMAX = X(1)
      IMAX = 1

      DO I=2*INC, N, INC
        IF (X(I) .GT. XMAX) THEN
          XMAX = X(I)
          IMAX = I
        ENDIF
      ENDDO

      ISMAX = IMAX

      RETURN
      END
