      SUBROUTINE FEOP_STEAL (XLAMK, N, INDFEACT, MAXFEACT, LASTFE,
     >                    SIGMAFE, OPAFE, ETAFE, INDRB,
     >                    IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >                    VDOPFE, DXFE, XLAM0FE,
     >                    ELEVEL, WEIGHT, RSTAR, POPNUM, ENTOT, T,
     >                    BNUEFE)

C***********************************************************************
C***  NON-LTE iron opacity at given frequency for only one depth point.
C***  This subroutine corresponds exactly to CMFFEOP, but without depth loop.
C***  Called from SETXJFINE (main program STEALCL)
C***********************************************************************
 

      DIMENSION SIGMAFE(1), INDRB(LASTFE)
      DIMENSION IFRBSTA(LASTFE), IFRBEND(LASTFE)
      DIMENSION IFENUP(LASTFE), IFELOW(LASTFE)
      DIMENSION POPNUM(N), WEIGHT(N) 
      DIMENSION INDFEACT(LASTFE) 
      DIMENSION ELEVEL(N)

      LOGICAL BNUEFE

C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C     ( CGS UNITS )
      DATA C2 / 3.9724E-16 /
C***  CLIGHT = VELOCITY OF LIGHT IN KM/S
      DATA CLIGHT /2.99792458E5/

C***  CALCULATE FE-FREQUECY-INDEX
      XINDF = - ALOG10(XLAM0FE/XLAMK) /
     >          ALOG10(1. + VDOPFE*DXFE/CLIGHT)
      NINDF = NINT(XINDF-0.5)
      Q = XINDF - FLOAT(NINDF)
      P = 1. - Q

C***  WAVENUMBERS IN KAYSER ( = CM**-1)
      WAVNUMK = 1.E8 / XLAMK
      WAVNUMK3 = WAVNUMK * WAVNUMK * WAVNUMK 

C***  PRESET OPAFE, ETAFE
         OPAFE = 0.
         ETAFE = 0.
      
C***  LOOP OVER ACTIVE FE-LINES
      DO 10 INDACT=1, MAXFEACT

         IND = INDFEACT(INDACT)
         LOW = IFELOW(IND)
         NUP = IFENUP(IND)
         WEIGHTLU = WEIGHT(LOW) / WEIGHT(NUP)
         WAVNUM0 = ELEVEL(NUP) - ELEVEL(LOW)
         WAVNUM03 = WAVNUM0 * WAVNUM0 * WAVNUM0 

C***    FIND INDICES IN ARRAY 'SIGMAFE'
         IND1 = INDRB(IND) + NINDF - IFRBSTA(IND)
C***    READ NEIGHBOURING VALUES
         SIGMA1 = SIGMAFE(IND1)
         SIGMA2 = SIGMAFE(IND1 + 1)
C***    INTERPOLATION
         SIGMA = P * SIGMA1 + Q * SIGMA2

C***       Opacity and Emissivity
            IF (BNUEFE) THEN
               G = WEIGHTLU * EXP(C1*(WAVNUM0-WAVNUMK)/T)
            ELSE
               G = WEIGHTLU 
            ENDIF
            EMINDU = G * POPNUM(NUP) * SIGMA
            SUM = POPNUM(LOW) * SIGMA - EMINDU
            IF (SUM .GT. 0.) THEN
               OPAFE = OPAFE + SUM
               ETAFE = ETAFE + EMINDU 
            ENDIF
 10   CONTINUE

      OPAFE = OPAFE * ENTOT * RSTAR
      IF (BNUEFE) THEN
         ETAFE = ETAFE * ENTOT * RSTAR * C2 * WAVNUMK3
      ELSE
         ETAFE = ETAFE * ENTOT * RSTAR * C2 * WAVNUM03
      ENDIF

      RETURN
      END
