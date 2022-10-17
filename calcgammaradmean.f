      SUBROUTINE CALCGAMMARADMEAN(ARAD, AGRAV, RADIUS, TAUROSS, 
     >                            ND, RCON, RI, GAMMARADMEAN)
C***********************************************************************
C**** Calculates the mean full Eddington Gamma (based on ARAD/AGRAV)
C****    in the hydrostatic part of a model. This is done by a
C****    weighted mean, starting from TAUROSS_MIN=0.1 (but inside RCON)
C****    and with weights exp(-TAUROSS) (full Rosseland-tau!)
C****
C**** Note: This routine does not limit GAMMARADMEAN to any maximum 
C****       value (This is intended and should be done afterwards!)
C****
C**** called from STEAL subroutines INITFCORR, CALCMASSFROMGEFF
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'

      INTEGER, INTENT(IN) :: ND
      
      REAL, DIMENSION(ND), INTENT(IN) :: TAUROSS, RADIUS
      REAL, DIMENSION(ND-1), INTENT(IN) :: ARAD, AGRAV, RI
      
      REAL, INTENT(IN) :: RCON
      REAL, INTENT(OUT) :: GAMMARADMEAN
      
      REAL :: TAUNORM, WTAU, RINT, TAUINT, RLAST, TAULAST,
     >        ARADINT, AGRAVINT, GAMMAINT
      INTEGER :: L, LCON      
      
      REAL, PARAMETER :: TAUROSS_MIN = 0.1

      GAMMARADMEAN = 0.
      TAUNORM = 0.
      LCON = ND-1
      IF (RCON > 0. .AND. RCON < RADIUS(1)) THEN
        DO L=ND-1, 1, -1
          TAUINT = 0.5 * (TAUROSS(L) + TAUROSS(L+1))
          IF (RADIUS(L) > RCON .OR. TAUINT < TAUROSS_MIN) THEN
            !add contribution of last, incomplete grid cell
            IF (RADIUS(L) > RCON) THEN
              RLAST = RCON
              CALL SPLINPOX(TAULAST, RCON, TAUROSS, RADIUS, ND)
            ELSEIF (TAUINT < TAUROSS_MIN) THEN
              CALL SPLINPOX(RLAST, TAUROSS_MIN, RADIUS, TAUROSS, ND)
              TAULAST = TAUROSS_MIN
            ENDIF
            RINT = 0.5 * (RLAST + RADIUS(LCON))
            TAUINT = 0.5 * (TAULAST + TAUROSS(LCON))
            WTAU = (TAUROSS(LCON) - TAULAST) * EXP(-TAUINT)
            
            CALL SPLINPOX(ARADINT, RINT, ARAD, RI, ND-1)
            CALL SPLINPOX(AGRAVINT, RINT, AGRAV, RI, ND-1)
            GAMMAINT = ARADINT / AGRAVINT
            GAMMARADMEAN = GAMMARADMEAN + GAMMAINT * WTAU
            TAUNORM = TAUNORM + WTAU
            EXIT            
          ENDIF
          LCON = L
          WTAU = (TAUROSS(L+1) - TAUROSS(L)) * EXP(-TAUINT)
          GAMMARADMEAN = GAMMARADMEAN + ARAD(L)/AGRAV(L) * WTAU
          TAUNORM = TAUNORM + WTAU
        ENDDO
      ENDIF
      IF (LCON /= ND-1) THEN
        GAMMARADMEAN = GAMMARADMEAN / TAUNORM
      ELSE
        !If there is no hydrostatic part, take innermost value for GAMMARADMEAN
        GAMMARADMEAN = ARAD(ND-1)/AGRAV(ND-1)
      ENDIF      

      RETURN
      END
