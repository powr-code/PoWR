      SUBROUTINE MDMV (A,B,JMAX,NP)
C***********************************************************************
C***  MATRIX (DIAGONAL)  A  *  MATRIX (VOLL)  B
C***  ERGEBNIS-MATRIX UEBERSCHREIBT  B
C***********************************************************************

      DIMENSION A(NP),B(NP,NP)

      DO 1 I=1,JMAX
      AI=A(I)
      DO 1 K=1,JMAX
    1 B(I,K)=B(I,K)*AI
      RETURN
      END
