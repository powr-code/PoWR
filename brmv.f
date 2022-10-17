      SUBROUTINE BRMV (V, A, TV, N, NDIM)
C********************************************************************
C***  V = A * V
C********************************************************************

      DIMENSION A(NDIM,NDIM), V(NDIM), TV(NDIM)

      DO 1 I=1,N
        TV(I) = 0.
1     CONTINUE

      DO 2 I=1,N
        DO 2 J=1,N
          TV(I) = TV(I) + A(I,J) * V(J)
2     CONTINUE

      DO 3 I=1,N
        V(I) = TV(I)
3     CONTINUE

      RETURN
      END
