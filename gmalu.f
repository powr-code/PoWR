      SUBROUTINE GMALU (GA,U,V,LMAX)
C***********************************************************************
C***  ALGEBRAIC ROUTINE CALLED FROM CMFRAY
C***********************************************************************
      DIMENSION GA (LMAX),U(LMAX),V(LMAX)

      LZ=LMAX-1
      DO 1 L=1,LZ
    1 V(L)=GA(L)*(U(L)-U(L+1))

      RETURN
      END
