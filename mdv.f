      SUBROUTINE MDV (A,W,N)
C***********************************************************************
C*** MATRIX A (DIAGONAL)  *  VEKTOR W
C***  ERGEBNIS-VEKTOR UEBERSCHREIBT  W
C***********************************************************************

      DIMENSION  A(N),W(N)

      DO 1 I=1,N
    1 W(I)=A(I)*W(I)

      RETURN
      END
