      SUBROUTINE POPMIN_NULLING (ZERO_RATES, POPNUM, ND, N)
C***********************************************************************
C***  Sets all POPNUMS that were flagged by ZERO_RATES to zero 
C***  Reason: these level populations have been set to POPMIN in steal 
C***      but should be set to 0.0 in all radiative-transfer programs
C***      in oder to avoid any artificial contributions to the 
C***      emissivities and opacities 
C***  Called from: WRCONT, COMO, COLI, FORMAL  
C***********************************************************************

      INTEGER, INTENT(IN) :: ND, N

      REAL,    DIMENSION(ND, N) :: POPNUM
      LOGICAL, DIMENSION(N, ND) :: ZERO_RATES
    
      INTEGER :: L

      DO L=1, ND
        DO J=1, N
          IF (ZERO_RATES(J,L)) THEN
            POPNUM(L,J) = .0
          ENDIF
        ENDDO        
      ENDDO

      RETURN
      END
