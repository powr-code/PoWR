      SUBROUTINE INITVELBETAPAR (RMAX, RCONin, VCON, 
     >                           bUseMaxForPar1)
C*******************************************************************************
C***  RAW INITIALIZATION OF THE VELOCITY-FIELD PARAMETERS
C***  unifies code parts that were previously in 
C***   INITVEL, VELTHIN and ENSURETAUMAX
C*******************************************************************************

      IMPLICIT NONE

      REAL RMAX, RCONin, VCON
      LOGICAL bUseMaxForPar1, bPAR1AtCon

      REAL :: VFINAL, VMAX, BETA, BETA2, BETA2FRACTION, VMIN,
     >        VPAR1, VPAR2, HSCALE, VPAR1_2, VPAR2_2, RCON,
     >        Q, S, T, RM

      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >            BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

       bPAR1AtCon = .NOT. bUseMaxForPar1

      RCON = RCONin

      VMAX = VFINAL * (1.-BETA2FRACTION)

      IF (BETA == .0) THEN
        Q = (Vcon/VMAX)**2
        VPAR2 = (RMAX**Q / RCON)**(1./(Q-1.))
        VPAR1 = VMAX / SQRT(ALOG(RMAX/VPAR2))
      ELSEIF (BETA < 0.) THEN
        Q = (VCON/VMAX)**(1./ABS(BETA))
        VPAR2 = (Q - 1.) / (Q/RMAX - 1./RCON)
        IF (bPAR1AtCon) THEN
          VPAR1 = Vcon / (1.-VPAR2/RCON)**(ABS(BETA))
        ELSE
          VPAR1 = VMAX / (1.-VPAR2/RMAX)**(ABS(BETA))         !INITVEL approach
        ENDIF
        IF (BETA2FRACTION > 0.) THEN
          VPAR2_2 = RCON
          VPAR1_2 = BETA2FRACTION * 
     >      VFINAL / (1.-VPAR2_2/RMAX)**(ABS(BETA2))     !note: 2BETA has always RMAX here
        ENDIF
      ELSE
        Q = (VCON/VMAX)**(1./BETA)
        S = (RMAX - 1. + RCON) / 2.
        T = (RMAX - Q * RCON) / (1.-Q) - RMAX * RCON
        VPAR2 = SQRT(S*S + T) - S
        IF (bPAR1AtCon) THEN
          RM = RCON
          IF (VPAR2 + RCON <= 1.) RM = 1. + ABS(VPAR2) + 1.E-10
          VPAR1 = Vcon/(1.-1./(VPAR2+RM))**BETA             !VELTHIN approach
        ELSE
          VPAR1 = VMAX / (1.-1./(VPAR2+RMAX))**BETA         !INITVEL approach
        ENDIF
        IF (BETA2FRACTION > 0.) THEN
          VPAR2_2 = 1. - RCON   
          VPAR1_2 = BETA2FRACTION * 
     >             VFINAL/(1.-1./(VPAR2_2+RMAX))**BETA2     !note: 2BETA has always RMAX here
        ENDIF
      ENDIF

      RETURN
      END
