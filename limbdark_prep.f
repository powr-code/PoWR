      SUBROUTINE LIMBDARK_PREP (LIMB_LINE, XOBS0, DXOBS, NFOBS,
     >                          XLAM, ALN, WEIGHT_LIMBDARK, BLIMB)

C***************************************************************************
C*** This routine decodes the LIMBDARKENING options
C*** and prepares the vector WEIGHT_LIMBDARK for the
C*** integration/extraction of the intensities I_nu(p) 
C*** The options allow to extract for either
C*** (1) a specific wavelength, 
C*** (2) a wavelenght range, or 
C*** (3) a broad-band color
C*** Syntax:
C*** (1) LIMB LAM x.x [plotoptions]
C***     provides I_nu(p) at given wavelength x.x 
C*** (2) LIMB RANGE x.x y.y [plotoptions]
C***     provides I_nu(p) averaged over the given wavelength x.x to y.y 
C*** (3) LIMB BAND U|B|V|u|b|v|y [plotoptions]
C***     provides I_nu(p) averaged over the specified filter 
C************************************************************************** 

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: NFOBS 
      CHARACTER LIMB_LINE*(*), SELECTEDBAND*8 
      REAL, INTENT(IN) :: DXOBS, XOBS0, XLAM, ALN
      REAL WEIGHT_LIMBDARK(NFOBS)

C***  For the filter funcions:
      INTEGER, PARAMETER :: MAXNFILT = 200
      INTEGER NFILT
      REAL, DIMENSION(MAXNFILT) :: FLAM, FILT

      REAL XLAMTEST, FILTVAL

C***  LOCAL VARIABLES:
      CHARACTER ACTPAR*80
      REAL DUMMY, XLAM_LIMBDARK, XLAM1, XLAM2, XK, XK1, XK2, P, Q 
      REAL XOBS1, XLAMMIN, XLAMMAX, XBLUE, XRED, WEIGHTSUM, DNUEDX
      INTEGER :: I, J, L, NPAR, PLIPO, K, K1, K2, IDX
      LOGICAL BLIMB, BDEFINED
      
      BLIMB = .TRUE.

C***  frequency range is from XOBS0 to XOBS1
C***  X... are dimensionless frequencies
      XBLUE = -(XOBS0 - DXOBS)
      XRED  = -(XOBS0 - DXOBS*NFOBS)
C***  current LAMBDA  range
      XLAMMIN = XLAM * EXP(-XBLUE*ALN)
      XLAMMAX = XLAM * EXP(-XRED *ALN)

ccc      write (0,*) 'XBLUE=', XBLUE
ccc      write (0,*) 'XRED =', XRED
ccc      write (0,*) 'XLAMMIN=', XLAMMIN
ccc      write (0,*) 'XLAMMAX=', XLAMMAX


      DO K=1, NFOBS
         WEIGHT_LIMBDARK(K) = .0
      ENDDO

C***  Check: only one BAND or LAMBDA or RANGE can be specified !
      BDEFINED = .FALSE.

C***  Decode the line LIMB_LINE specified in FORMAL_CARDS
      CALL SARGC(LIMB_LINE, NPAR)
      IF ((NPAR .LT. 2) .OR. (NPAR .GT. 12)) GOTO 100

      DO I=2, NPAR
         CALL SARGV(LIMB_LINE, I, ACTPAR)

C*******************************************************************
         IF (ACTPAR(:3) .EQ. 'LAM') THEN
C                             ===
            IF (BDEFINED) GOTO 90
            BDEFINED = .TRUE.
            IF (I+1 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, I+1, ACTPAR)
            READ (ACTPAR, '(F12.0)', ERR=100) XLAM_LIMBDARK

C***        Skip if XLAM_LIMBDARK-XLAMMIN not covered:
            IF ((XLAM_LIMBDARK-XLAMMIN)*(XLAM_LIMBDARK-XLAMMAX) .GE. .0) 
     >         GOTO 200

C*          Find non-integer K-index of requested XLAM_LIMBDARK
            XK   = ( XOBS0 + (ALOG(XLAM)-ALOG(XLAM_LIMBDARK)) / ALN )
     >             / DXOBS

            K = INT(XK)
            Q = XK - FLOAT(K)
            P = 1. - Q
            WEIGHT_LIMBDARK(K)   = P
            WEIGHT_LIMBDARK(K+1) = Q
               

ccc               write (0,*) 'XLAM_LIMBDARK lies inside:', XLAM_LIMBDARK
ccc               xlamtest = XLAM * EXP ((XOBS0 - K*DXOBS) * ALN)  
ccc               write (0,*) 'K, lam =', K,   XLAMtest, 1-q
ccc               xlamtest = XLAM * EXP ((XOBS0 - (K+1)*DXOBS) * ALN)  
ccc               write (0,*) 'K, lam =', K+1, XLAMtest, q

C*******************************************************************
         ELSEIF (ACTPAR .EQ. 'RANGE') THEN
C                             =====
            IF (BDEFINED) GOTO 90
            BDEFINED = .TRUE.
            IF (I+2 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, I+1, ACTPAR)
            READ(ACTPAR, '(F12.0)', ERR=100) XLAM1
            CALL SARGV(LIMB_LINE, I+2, ACTPAR)
            READ(ACTPAR, '(F12.0)', ERR=100) XLAM2

C***        Sort lam1 < lam2 
            IF (XLAM1 .GT. XLAM2) THEN
                DUMMY = XLAM1
                XLAM1 = XLAM2
                XLAM2 = DUMMY
            ENDIF   

C***        Skip RANGE if not covered:
            IF ((XLAM1-XLAMMIN)*(XLAM1-XLAMMAX) .GE. .0) GOTO 201 
            IF ((XLAM2-XLAMMIN)*(XLAM2-XLAMMAX) .GE. .0) GOTO 201 

C*          Find non-integer K-index of requested XLAM1
            XK1   = ( XOBS0 + (ALOG(XLAM)-ALOG(XLAM1)) / ALN )
     >             / DXOBS

            K1 = INT(XK1)
            Q = XK1 - FLOAT(K1)
            P = 1. - Q
            WEIGHT_LIMBDARK(K1)   = P*P / 2.
            WEIGHT_LIMBDARK(K1+1) = P * (Q + 1.) / 2.
               
C*          Find non-integer K-index of requested XLAM2
            XK2   = ( XOBS0 + (ALOG(XLAM)-ALOG(XLAM2)) / ALN )
     >             / DXOBS

            K2 = INT(XK2)
            Q = XK2 - FLOAT(K2)
            P = 1. - Q
            WEIGHT_LIMBDARK(K2)   = Q * (P + 1.) / 2.
            WEIGHT_LIMBDARK(K2+1) = Q*Q / 2.
            DO K=K1, K2
               IF (K .EQ. K1 .OR. K .EQ. K1) THEN 
                  WEIGHT_LIMBDARK(K) = 0.5
               ELSE
                  WEIGHT_LIMBDARK(K) = 1.
               ENDIF
            ENDDO

C***        Integration is done over \nu. Unlike \Delta x, \Delta\nu is 
C***        not constant but propto \nu (cf. manpowr) 
            DO K=K1, K2+1
               DNUEDX = EXP(ALN)**(XOBS0-K*DXOBS)
               WEIGHT_LIMBDARK(K)   = WEIGHT_LIMBDARK(K) * DNUEDX
            ENDDO

C*******************************************************************
         ELSE IF (ACTPAR .EQ. 'BAND') THEN
C                              ====
            IF (BDEFINED) GOTO 90
            BDEFINED = .TRUE.
            IF (I+1 .GT. NPAR) GOTO 100
            CALL SARGV(LIMB_LINE, I+1, SELECTEDBAND)
            
            CALL FILTERFUNCTIONS (SELECTEDBAND, NFILT, FLAM, FILT)

C***        Skip RANGE if filter not covered:
            XLAM1 = FLAM(1)
            XLAM2 = FLAM(NFILT)
            IF ((XLAM1-XLAMMIN)*(XLAM1-XLAMMAX) .GE. .0) GOTO 202 
            IF ((XLAM2-XLAMMIN)*(XLAM2-XLAMMAX) .GE. .0) GOTO 202 

C*          Find non-integer K-index of first filter-wavelength:
            XK1   = ( XOBS0 + (ALOG(XLAM)-ALOG(XLAM1)) / ALN )
     >             / DXOBS

            K1 = INT(XK1)
            Q = XK1 - FLOAT(K1)
            P = 1. - Q
            WEIGHT_LIMBDARK(K1)   = P*P / 2.
            WEIGHT_LIMBDARK(K1+1) = P * (Q + 1.) / 2.
               
C*          Find non-integer K-index of last filter-wavelength:
            XK2   = ( XOBS0 + (ALOG(XLAM)-ALOG(XLAM2)) / ALN )
     >             / DXOBS

            K2 = INT(XK2)
            Q = XK2 - FLOAT(K2)
            P = 1. - Q
            WEIGHT_LIMBDARK(K2)   = Q * (P + 1.) / 2.
            WEIGHT_LIMBDARK(K2+1) = Q*Q / 2.
            DO K=K1, K2
               IF (K .EQ. K1 .OR. K .EQ. K1) THEN 
                  WEIGHT_LIMBDARK(K) = 0.5
               ELSE
                  WEIGHT_LIMBDARK(K) = 1.
               ENDIF
            ENDDO

C***        Integration is done over \nu. Unlike \Delta x, \Delta\nu is 
C***        not constant but propto \nu (cf. manpowr) 
ccc            DO K=K1, K2+1
            DO K=K1+1, K2
               XLAMTEST = XLAM * EXP ((XOBS0 - K*DXOBS) * ALN)  
               CALL SPLINPO (FILTVAL, XLAMTEST, FILT, FLAM, NFILT)
               DNUEDX = EXP(ALN)**(XOBS0-K*DXOBS)
               WEIGHT_LIMBDARK(K) = WEIGHT_LIMBDARK(K) * DNUEDX *FILTVAL
            ENDDO

         ENDIF
      ENDDO

      IF (.NOT. BDEFINED) GOTO 95

C***  Normalization of frequency integral
      WEIGHTSUM = .0
      DO K=1, NFOBS
         WEIGHTSUM = WEIGHTSUM + WEIGHT_LIMBDARK(K)
      ENDDO

      DO K=1, NFOBS
         WEIGHT_LIMBDARK(K) = WEIGHT_LIMBDARK(K) / WEIGHTSUM
      ENDDO

C***  Test output
      if (.false.) then 
         do k=1, nfobs
            if (WEIGHT_LIMBDARK(K) .EQ. .0) CYCLE
            xlamtest = XLAM * EXP ((XOBS0 - K*DXOBS) * ALN)  
            write (0,*) 'K, lam =', K,   XLAMtest, WEIGHT_LIMBDARK(K)
         enddo
      endif


      RETURN

C******************************************************************
C***  Skipping LIMBDARK if LAMBDA, RANGE or FILTER is not covered 
C***  by the present wavelength range

  200 WRITE (0,'(A, F8.0,A)') '*** WARNING: ' // 
     >      'LIMBDARK: LAM=', XLAM_LIMBDARK,  
     >      ' lies outside the current RANGE'
      GOTO 210

  201 WRITE (0,'(A, 2F8.0,A)') '*** WARNING: ' // 
     >      'LIMBDARK: interval LAM=', XLAM1, XLAM2,  
     >      ' not covered by current RANGE'
      GOTO 210

  202 WRITE (0,'(A, F8.0,A, F8.0, A)') '*** WARNING: ' // 
     >      'LIMBDARK: FILTER ' // SELECTEDBAND(:IDX(SELECTEDBAND)) //
     >      ' (', XLAM1, ' - ', XLAM2, 
     >             ') not covered by current RANGE'
      GOTO 210

  210 CONTINUE
      BLIMB = .FALSE.
      WRITE (0,'(A, F10.0,A,F10.0)') '*** --> LIMBDARKENING skipped ' //
     >            '- current range is ', XLAMMIN, ' to ', XLAMMAX
      RETURN 

C***  ERROR BRANCHES ************************************************

   90 WRITE (0,'(A)') 
     >  '*** ERROR: only one LAMBDA, RANGE or BAND can be selected'
      GOTO 110

   95 WRITE (0,'(A)') 
     >  '*** ERROR: either LAMBDA, RANGE or BAND must be selected'
      GOTO 110

  100 WRITE(0,'(A)') "*** ERROR: Wrong syntax"
      GOTO 110

  110 CONTINUE
      WRITE(0,'(A)') "The error occured in the following line"
      WRITE(0,'(A)') LIMB_LINE(:IDX(LIMB_LINE))
      WRITE(0,'(A)') "Possible syntax:"
      WRITE(0,'(A)') 'LIMB LAM x.x [plotoptions]'
      WRITE(0,'(A)') 'LIMB RANGE x.x y.y [plotoptions]'
      WRITE(0,'(A)') 'LIMB BAND U|B|V|u|b|v|y [plotoptions]'

      WRITE(0,'(A)')  'Fatal error in subroutine LIMBDARK_PREP'
      WRITE(0,'(A)')  '--> LIMB DARKENING not calculaated (skipped)'
      
      BLIMB = .FALSE.
      RETURN

      END
