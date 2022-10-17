      SUBROUTINE MACROCLUMP (OPA, POROLENGTH, CORRFAC_MACROCLUMP)
C*****************************************************************
C***  For a given opacity OPA and porosity length, this routine 
C***  returns the correction factor: OPA_EFF = OPA * CORRFAC_MACROCLUMP 
C***  due to clumps of the optical thickness TAUCLUMP = OPA * POROLENGTH
C*****************************************************************

      IF (POROLENGTH .EQ. .0) THEN
         CORRFAC_MACROCLUMP = 1.
         RETURN
      ENDIF
 
      TAUCLUMP = OPA * POROLENGTH
 
      IF (TAUCLUMP .GT. 1.E-10) THEN
         CORRFAC_MACROCLUMP = (1. - EXP(-TAUCLUMP)) / TAUCLUMP
      ELSE IF (TAUCLUMP .GT. .0) THEN
         CORRFAC_MACROCLUMP = 1. - TAUCLUMP * 0.5
      ELSE
         CORRFAC_MACROCLUMP = 1.
      ENDIF
 
      RETURN
      END
