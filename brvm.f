      SUBROUTINE BRVM (V, A, TV, N, NDIM)
C********************************************************************
C***  V = V * A
C********************************************************************
      DIMENSION  V(NDIM), A(NDIM,NDIM), TV(NDIM)


      DO 1 I=1,N
        TV(I) = 0.
1     CONTINUE

      DO 2 I=1,N
        DO 2 J=1,N
          TV(I) = TV(I) + V(J) * A(J,I) 
2     CONTINUE

      DO 4 I=1,N
        V(I) = TV(I)
4     CONTINUE

      RETURN
      END
