      SUBROUTINE MSUB (A,B,JMAX,NP)
C***********************************************************************
C***  A := A - B
C***********************************************************************

      DIMENSION A(NP,NP),B(NP,NP)

      DO 1 K=1,JMAX
      DO 1 I=1,JMAX
    1 A(I,K)=A(I,K)-B(I,K)

      RETURN
      END
