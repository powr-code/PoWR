      FUNCTION SDOT (N, X, IX, Y, IY)

      DIMENSION X(N), Y(N)

      SUM = 0.

      DO I=0, N-1
        SUM = SUM + (X(IX+I) * Y(IY+I))
      ENDDO

      SDOT = SUM

      RETURN
      END
