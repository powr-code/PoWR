      SUBROUTINE CALC_S(S, ETAK, OPAK, ND)
C****************************************************************
C***  Calculation of the Source Function
C***    Called by COLI
C****************************************************************

      DIMENSION S(ND), ETAK(ND), OPAK(ND)

      DO L=1, ND
        S(L) = ETAK(L) / OPAK(L)
      ENDDO

      RETURN
      END
