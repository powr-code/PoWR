      SUBROUTINE BRVVDY (A,V1,V2,N,NDIM)
C*****************************************************************
C***  MULTIPLICATION DYADE = VECTOR * VECTOR ;
C***  THE RESULTING MATRIX IS IMMEDIATELY SUBSTRACTED FROM MATRIX A
C***  A := A - V1 * V2
C*****************************************************************

      DIMENSION V1(NDIM), V2(NDIM), A(NDIM,NDIM)

      DO 1 I=1,N
        DO 2 J=1,N
          A(I,J) = A(I,J) - V1(I) * V2(J)
    2   CONTINUE
    1 CONTINUE

      RETURN
      END
