      FUNCTION STARKHOLTSMARK(XI, VDOP, VDOPDU_FINE, GRIEMPAR,
     >                        NLOC, MAXLAP, ZFINE, LRIP, RFINERIP, 
     >                        PR, QR, RADIUS, ND, PJPJ, bHYDROGEN)
C***********************************************************************      
C***  Called from: ZONEINT
C***  Profile function for pressure-broadened line opacity profiles
C***   that are broadened by the linear Stark effect
C****  
C***  GRIEMPAR = convertion factor for Griem's beta parameter 
C***              (has been multiplied with VDOP in STARKBROAD)
C***
C***  
C***  LRIP (1...ND) gives a start index for the interpolation over radius 
C***  If RFINERIP .NE. .0 it is assumed that interpolation weights for 
C***  radius are already known (function multiply called for different lines)
C***  HENCE DON'T FOGET TO SET RFINERIP = .0 befor calling this function
C***  at a new spatial point!
C***********************************************************************      

      INTEGER, INTENT(IN) :: ND

      REAL, DIMENSION(MAXLAP,ND) :: GRIEMPAR
      REAL, DIMENSION(ND) :: RADIUS
      
      REAL, INTENT(IN) :: VDOP, VDOPDU_FINE
      REAL, INTENT(INOUT) :: RFINERIP
      
      REAL :: BETA, BETADOP, FGRIEM, DLAMDOPREL
      
      LOGICAL :: bHYDROGEN

      REAL, PARAMETER :: CLIGHT  = 2.99792458E5    ! c in km / s
      
      
C***  Special interpolation in radius for efficiency:
C***  most likely, the radius index did not change
      IF (RFINERIP .EQ. .0) THEN
         RFINERIP = SQRT(PJPJ+ZFINE*ZFINE)
         RFINERIP = AMAX1 (RFINERIP, 1.0)
         RFINERIP = AMIN1 (RFINERIP, RADIUS(1))
  200    CONTINUE

         IF (RFINERIP .LE. RADIUS(LRIP)) THEN
            IF (RFINERIP .GE. RADIUS(LRIP+1)) THEN
               QR = (RADIUS(LRIP) - RFINERIP)/
     >              (RADIUS(LRIP) - RADIUS(LRIP+1))
               PR = 1. - QR
            ELSE
               LRIP = LRIP + 1
               GOTO 200
            ENDIF
         ELSE
            LRIP = LRIP - 1
            GOTO 200
         ENDIF
      ENDIF
C***  End of radius interpolation weights: PR, QR established

      FGRIEM = PR * GRIEMPAR(NLOC,LRIP) + QR * GRIEMPAR(NLOC,LRIP+1)
      
      DLAMDOPREL = VDOP / CLIGHT
      ALN = LOG(1. - DLAMDOPREL)
      BETA = FGRIEM * ABS( EXP(ALN*XI) - 1. )
      BETADOP = FGRIEM * DLAMDOPREL
      
      STARKHOLTSMARK = PHIHOLTSMARK(BETA, BETADOP, bHYDROGEN)

      RETURN 
      END
