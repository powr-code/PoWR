      SUBROUTINE ZINCONE(ZINTER, THETA, ALPHA, LPSTA, LPEND, NPHI, 
     >                   P, NP)
C******************************************************************
C***  CALCULATES THE POINTS OF INTERSECTION FOR THE CONE MODEL
C***  THETA IS THE OPENING ANGLE OF THE CONE
C***  ALPHA IS THE ANGLE OF THE INCLINATION OF THE CONE
C******************************************************************

      DIMENSION ZINTER(2, NP, NPHI), P(NP)

      DATA PI / 3.14159265358979 /

      SA  = SIN(ALPHA)
      SA2 = SA * SA
      CA  = COS(ALPHA)
      CA2 = CA * CA
      TT  = TAN(THETA)
      TT2 = TT * TT
      XNENN = (CA2 - (SA2*TT2))
c      write (*,*) 'xnenn=',xnenn
c      stop 'test in zin'
      IF (ABS(XNENN) .LT. 1.E-2) THEN
        WRITE (0,*) 'DENOMINATOR TOO SMALL :',XNENN
        WRITE (0,*) 'ALPHA, THETA=', ALPHA, THETA
        WRITE (0,*) 'PLEASE USE A COMBINATION OF ALPHA AND THETA ',
     >              'WITH ALPHA + THETA .NE. 90DEG',
     >              'THE DISTANCE SHOULD BE AT LEAST 0.5DEG'
        STOP 'ERROR IN SUBR. ZIN'
      ENDIF

      DO LPHI=LPSTA, LPEND
        IF (NPHI .GT. 1) THEN
          PHI=PI * (LPHI-1) / (NPHI-1)
        ENDIF
        SPHI  = SIN(PHI)
        SPHI2 = SPHI * SPHI
        CPHI  = COS(PHI)
        CPHI2 = CPHI * CPHI
C***  QUADP AND QUADQ ARE THE TWO ARGUMENTS OF THE SOLUTION FORMULA
C***    OF AN QUADRATIC EQUATION
        QUADP = (-2. * CPHI * CA * SA * (1.-TT2)) /
     >          XNENN
        QUADQ = (CPHI2 * (SA2 - (CA2*TT2)) + SPHI2) /
     >          XNENN
        XSQR2 = (QUADP*QUADP / 4.) - QUADQ
        IF (XSQR2 .GE. 0.) THEN
          XSQR = SQRT(XSQR2)
          Y1 = -QUADP/2. + XSQR
          Y2 = -QUADP/2. - XSQR
          DO JP=1, NP
            XP  = P(JP)
            ZINTER(1, JP, LPHI) = Y1 * XP
            ZINTER(2, JP, LPHI) = Y2 * XP
          ENDDO
        ELSE
          DO JP=1, NP
            ZINTER(1, JP, LPHI) = 0.
            ZINTER(2, JP, LPHI) = 0.
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
