      FUNCTION VELOBETA(R)
C***********************************************************************
C***  VELOCITY FIELD from analytic (beta) law
C***  BETA < 0: different fine parametrization of beta law with abs(beta)
C***  BETA = 0: SWITCH TO SQRT-LOG-LAW
C***
C***  called by WRVEL, DELTAGRTHIN, VELTHIN
C***********************************************************************

      IMPLICIT NONE

      REAL :: VELOBETA
      REAL, INTENT(IN) :: R
      
      REAL :: VFINAL, VMIN, BETA, VPAR1, VPAR2, RCON, HSCALE,
     >        BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2, RM 
      LOGICAL :: bSMOCO

      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >               BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

        IF (BETA == 0.) THEN
          VELOBETA = VPAR1 * SQRT(LOG(R/VPAR2))
        ELSEIF (BETA < 0.) THEN
          VELOBETA = VPAR1 * (1.-VPAR2/R)**(ABS(BETA))
          IF (BETA2FRACTION > .0) THEN 
             VELOBETA = VELOBETA + 
     >                VPAR1_2 * (1.-VPAR2_2/R)**(ABS(BETA2))
          ENDIF
        ELSE
C***     New parametrization of beta law, wrh  5-Mar-2001 
          RM = (VPAR2+R)
          IF (RM <= 1.) RM = 1. + 1.E-10
          VELOBETA = VPAR1 * (1.-1./RM)**BETA
          IF (BETA2FRACTION > .0) THEN 
             VELOBETA = VELOBETA + 
     >                VPAR1_2 * (1.-1./(VPAR2_2+R))**BETA2   
          ENDIF
        ENDIF

      RETURN
      
      END
