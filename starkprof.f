      FUNCTION STARKPROF (XI, IPOINTERPHITAB, ZFINE, LRIP, RFINERIP, 
     >          PR, QR, PHITAB, NFDIMPHITAB, NLDIMPHITAB, RADIUS, ND, 
     >          NDDIM, PJPJ, DXMAX)
C***********************************************************************      
C***  Called from: ZONEINT
C***  Profile function for pressure-broadened line opacity profiles
C***  (used for H I and He II)
C***  The profile is obtained by interpolation in the table PHITAB
C***  Table PHITAB must have been prepared by SUBROUTINE STARKBROAD
C***  The index IPOINTERPHITAB is the pointer to the line index in that table  
C***  LRIP (1...ND) gives a start index for the interpolation over radius 
C***  If RFINERIP .NE. .0 it is assumed that interpolation weights for 
C***  radius are already known (function multiply called for different lines)
C***  HENCE DON'T FOGET TO SET RFINERIP = .0 befor calling this function
C***  at a new spatial point!
C***********************************************************************      

      DIMENSION PHITAB(-NFDIMPHITAB:NFDIMPHITAB, NDDIM, NLDIMPHITAB)
      DIMENSION RADIUS(ND)

C***  Prepare frequency interpolation
      XIK = XI / DXMAX
      K = INT(XIK)

C***  If frequency outside the bandwidth: PHI = 0.0
      IF (K .LE. -NFDIMPHITAB .OR. K .GE. NFDIMPHITAB) THEN
         STARKPROF = .0
         RETURN
      ENDIF

      Q = ABS( XIK - AINT(XIK) )
      IF (XI .LT. .0) THEN
         KK = K-1
      ELSE
         KK = K+1
      ENDIF

C***  Special interpolation in radius for efficiency:
C***  most likely, the radius index did not change
      IF (RFINERIP .EQ. .0) THEN
         RFINERIP = SQRT(PJPJ+ZFINE*ZFINE)
         RFINERIP = AMAX1 (RFINERIP, 1.0)
         RFINERIP = AMIN1 (RFINERIP, RADIUS(1))
  200    CONTINUE

ccc for tests only
ccc         IF (LRIP .LT. 1 .OR. LRIP .GT. ND) THEN
ccc            WRITE (0,*) 'RFINERIP=', RFINERIP
ccc            STOP '*** INTERNAL ERROR IN STARKPROF'
ccc         ENDIF


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

      IP = IPOINTERPHITAB
      PHITABK   = PR * PHITAB(K ,LRIP,IP) + QR * PHITAB(K ,LRIP+1,IP)
      PHITABKK  = PR * PHITAB(KK,LRIP,IP) + QR * PHITAB(KK,LRIP+1,IP)

      STARKPROF = (1.-Q) * PHITABK + Q * PHITABKK

      RETURN
      END
