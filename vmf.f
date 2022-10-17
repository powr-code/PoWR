      SUBROUTINE VMF (V2,V1,A,N,NDIM)
C***********************************************************************
C***  MULTIPLICATION VECTOR = VECTOR * MATRIX (FULL)  --  V2 = V1 * A
C***********************************************************************
      DIMENSION V1(NDIM),V2(NDIM),A(NDIM,NDIM)

      DO 1 J=1,N
      SUM=.0
      DO 2 I=1,N
    2 SUM=SUM+V1(I)*A(I,J)
    1 V2(J)=SUM

      RETURN
      END
