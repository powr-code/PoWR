      FUNCTION STARKVOIGT (XI, AVOIGT, NLOC, MAXLAP, 
     >     ZFINE, LRIP, RFINERIP, PR, QR, RADIUS, ND, PJPJ)
C***********************************************************************      
C***  Called from: ZONEINT
C***  Profile function for pressure-broadened line opacity profiles
C***   that are represented by a Voigt function with depth-dependent
C***   Voigt parameter AVOIGT
C***  The profile is obtained by interpolation in AVOIGT and call
C***  of the VOIGTH function. 
C***  Table AVOIGT must have been prepared by SUBROUTINE STARKBROAD
C***  LRIP (1...ND) gives a start index for the interpolation over radius 
C***  If RFINERIP .NE. .0 it is assumed that interpolation weights for 
C***  radius are already known (function multiply called for different lines)
C***  HENCE DON'T FOGET TO SET RFINERIP = .0 befor calling this function
C***  at a new spatial point!
C***********************************************************************      

      DIMENSION AVOIGT(MAXLAP,ND)
      DIMENSION RADIUS(ND)

C***  Special interpolation in radius for efficiency:
C***  most likely, the radius index did not change
      IF (RFINERIP .EQ. .0) THEN
         RFINERIP = SQRT(PJPJ+ZFINE*ZFINE)
         RFINERIP = AMAX1 (RFINERIP, 1.0)
         RFINERIP = AMIN1 (RFINERIP, RADIUS(1))
  200    CONTINUE

ccc for tests only
cc         IF (LRIP .LT. 1 .OR. LRIP .GT. ND) THEN
cc            WRITE (0,*) 'RFINERIP=', RFINERIP
cc            STOP '*** INTERNAL ERROR IN STARKPROF'
cc         ENDIF


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

      AV   = PR * AVOIGT(NLOC,LRIP) + QR * AVOIGT(NLOC,LRIP+1)

      STARKVOIGT = VOIGTH(AV,XI)

      RETURN
      END
