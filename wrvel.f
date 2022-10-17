      FUNCTION WRVEL(R)
C***********************************************************************
C***  VELOCITY FIELD
C***  NEWVELO = TRUE: ANALYTIC LAW;
C***  NEWVELO = FALSE: BELOW RCON: INTERPOLATION IN TABLE OLDVELO
C***            OVER OLDRADI
C***            (NEWVELO is called after VELTHIN routine has been passed)
C***  BETA .LE. 0: SWITCH TO SQRT-LOG-LAW
C***********************************************************************
      
      REAL, INTENT(IN) :: R
      INTEGER :: NDVAL

C***  COMMON/VELPAR/ TRANSFERS VELOCITY-FIELD PARAMETERS TO FUNCTION WRVEL
C***  The second line has additional parameters for the 2-beta-law
C***     -- wrh  6-Apr-2006 17:21:57
      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >     BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

C***  (OLD) RADIUS AND VELOCITY ARE TRANSFERRED TO FUNCTION WRVEL 
C***     BY SPECIAL COMMON BLOCKS:
      COMMON /COMRADI/ OLDRADI(2)
      COMMON /COMVELO/ NEWVELO, NDVAL, OLDVELO(2)
      LOGICAL NEWVELO
      
      REAL, EXTERNAL :: VELOBETA

      IF (R .LE. 1.) THEN
         WRVEL=VMIN
      ELSE IF (R. LT. RCON) THEN
         IF (NEWVELO) THEN
            WRVEL = (VMIN) * EXP(MIN(100.,(R-1.)/HSCALE))
         ELSE
           IF (R < OLDRADI(NDVAL)) THEN
             WRVEL = OLDVELO(NDVAL)
           ELSEIF (R > OLDRADI(1)) THEN
             WRVEL = OLDVELO(1)
           ELSE
             CALL SPLINPO(WRVEL, R, OLDVELO, OLDRADI, NDVAL)
           ENDIF
         ENDIF
      ELSE
C***    Calculate velocity via (double-)beta law      
        WRVEL = VELOBETA(R)
      ENDIF


      IF (WRVEL.GT.VFINAL) WRVEL=VFINAL

      RETURN

      END
