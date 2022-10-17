      SUBROUTINE PREPMACROCLUMP (MACROCLUMPLINE, DENSCON, VELO, RADIUS, 
     >            TAUROSS, ND, POROLENGTH)
C******************************************************************** 
C***  Decodes the MACROCLUMP option from input, and creates a vector
C***  over depth index  POROLENGTH(L)  
C******************************************************************** 

      CHARACTER MACROCLUMPLINE*(*), ACTPAR*20, ACTPAR2*20
      DIMENSION DENSCON(ND), VELO(ND), RADIUS(ND), POROLENGTH(ND)
      DIMENSION TAUROSS(ND)

      DATA PI / 3.141592654 /

C***  Defaults
      CLUMP_SEP0 = 1.
      TAU1 = -1.
      TAU2 = -1.
      VELO1 = -1.
      VELO2 = -1.

      CALL SARGC (MACROCLUMPLINE, NPAR)
      DO IPAR=2, NPAR, 2
         CALL SARGV (MACROCLUMPLINE, IPAR, ACTPAR)
         IF (NPAR .LT. IPAR+1) GOTO 90
         CALL SARGV (MACROCLUMPLINE, IPAR+1, ACTPAR2)
         IF (ACTPAR .EQ. 'CLUMP_SEP') THEN
            READ (ACTPAR2, '(F10.0)', ERR=91) CLUMP_SEP0
         ELSE IF (ACTPAR .EQ. 'TAU1') THEN
            READ (ACTPAR2, '(F10.0)', ERR=91) TAU1
         ELSE IF (ACTPAR .EQ. 'TAU2') THEN
            READ (ACTPAR2, '(F10.0)', ERR=91) TAU2
         ELSE IF (ACTPAR .EQ. 'VELO1') THEN
            READ (ACTPAR2, '(F10.0)', ERR=91) VELO1
         ELSE IF (ACTPAR .EQ. 'VELO2') THEN
            READ (ACTPAR2, '(F10.0)', ERR=91) VELO2
         ELSE
            GOTO 92
         ENDIF      
      ENDDO

      IF ((TAU1.NE.-1. .OR.  TAU2.NE.-1.) .AND. 
     >   (VELO1.NE.-1. .OR. VELO2.NE.-1.)) GOTO 93


      DO L=1, ND
         POROLENGTH(L) = DENSCON(L)**(2./3.) * CLUMP_SEP0 
     >       * (VELO(L) / VELO(1) * RADIUS(L)**2)**(1./3.)
      ENDDO


C***  Optional suppression of macroclumping at large Rosseland depth
C***    or at low velocities. Soft switch between two boundaries  

C***  Soft switch with TAU-Roseland parameters
      IF (TAU1 .NE. -1. .OR.  TAU2 .NE. -1.) THEN
         IF (TAU2 .LT. .0) TAU2 = TAU1
         IF (TAU1 .LT. .0) TAU1 = TAU2

C***   Put TAU1, TAU2 in increasing order
         IF (TAU2 .LT. TAU1) THEN
            TEMP = TAU1
            TAU1 = TAU2
            TAU2 = TEMP
         ENDIF 

         DO L=1, ND
            IF (TAUROSS(L) .LT. TAU1) CYCLE 
            IF (TAUROSS(L) .GT. TAU2) THEN
               POROLENGTH(L) = .0 
            ELSE
               W = (TAU2 - TAUROSS(L)) / (TAU2 - TAU1) 
               Q = 0.5 - 0.5 * COS(PI*W)
               POROLENGTH(L) = POROLENGTH(L) * Q 
            ENDIF
         ENDDO

C***  Soft-switch with VELOCITY parameters
      ELSEIF (VELO1 .NE. -1. .OR.  VELO2 .NE. -1.) THEN
         IF (VELO2 .LT. .0) VELO2 = VELO1
         IF (VELO1 .LT. .0) VELO1 = VELO2

C***   Put VELO1, VELO2 in increasing order
         IF (VELO2 .LT. VELO1) THEN
            TEMP = VELO1
            VELO1 = VELO2
            VELO2 = TEMP
         ENDIF 

         DO L=1, ND
            IF (VELO(L) .GT. VELO2) CYCLE 
            IF (VELO(L) .LT. VELO1) THEN
               POROLENGTH(L) = .0 
            ELSE
               W = (VELO2 - VELO(L)) / (VELO2 - VELO1) 
               Q = 0.5 + 0.5 * COS(PI*W)
               POROLENGTH(L) = POROLENGTH(L) * Q 
            ENDIF
         ENDDO

      ENDIF

C*** Output
      WRITE (*,'(A)') MACROCLUMPLINE(:IDX(MACROCLUMPLINE))
      WRITE (*,'(A)') 'PARAMETERS:'
      WRITE (*,'(A, G10.2)') 'CLUMP_SEP = ', CLUMP_SEP0

      IF (TAU1 .GE. .0) THEN 
         CALL LIPO (V1, TAU1, VELO, TAUROSS, ND)  
         WRITE (*,'(A, G10.2, A, G10.2)') 
     >       ' TAU1 = ', TAU1, '  ---> VELO1= ', V1
      ENDIF

      IF (TAU2 .GE. .0) THEN
         CALL LIPO (V2, TAU2, VELO, TAUROSS, ND)  
         WRITE (*,'(A, G10.2, A, G10.2)') 
     >       ' TAU2 = ', TAU2, '  ---> VELO2= ', V2
      ENDIF

      IF (VELO1 .GE. .0) THEN
         CALL LIPO (T1, VELO1, TAUROSS, VELO, ND)  
         WRITE (*,'(A, G10.2, A, G10.2)') 
     >       ' VELO1= ', VELO1, '  ---> TAU1 = ', T1
      ENDIF

      IF (VELO2 .GE. .0) THEN
         CALL LIPO (T2, VELO2, TAUROSS, VELO, ND)  
         WRITE (*,'(A, G10.2, A, G10.2)') 
     >       ' VELO2= ', VELO2, '  ---> TAU2 = ', T2
      ENDIF

cc      WRITE (*,*) 
cc      DO L=1, ND
cc         write (*,'(I3, G10.2)') l, porolength(l)
cc      ENDDO

      MACROCLUMPLINE = 'DONE'

      RETURN

C**** ERROR BRANCHES *******************
   90 WRITE (0,*) 
     >    'FATAL ERROR: MACROCLUMP OPTION HAS TOO FEW PARAMETERS'
      GOTO 99

   91 WRITE (0,*) 'FATAL ERROR: CANNOT DECODE ',ACTPAR2(:IDX(ACTPAR2)),
     >       ' AS FLOATING POINT NUMBER'
      GOTO 99

   92 WRITE (0,*) 'FATAL ERROR: UNRECOGNIZED PARAMETER: ',
     >       ACTPAR(:IDX(ACTPAR))
      GOTO 99

   93 WRITE (0,*) 'FATAL ERROR: TAU and VELO range cannot ',
     >      'be used simultaneously!'
      GOTO 99

   99 WRITE (0,*) 'THE ERROR OCCURRED IN THE FOLLOWING LINE:'
      WRITE (0,*) MACROCLUMPLINE
      STOP 'ERROR DETECTED IN SUBR. PREPMACROCLUMP'

       END
