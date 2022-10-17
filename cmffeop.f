      SUBROUTINE CMFFEOP (XLAMK, ND, N, INDFEACT, MAXFEACT, LASTFE,
     >                    SIGMAFE, OPAFE, ETAFE, INDRB,
     >                    IFRBSTA, IFRBEND, IFENUP, IFELOW,
     >                    CLIGHT, VDOPFE, DXFE, XLAM0FE,
     >                    ELEVEL, WEIGHT, RSTAR, POPNUM, ENTOT,
     >                    SIGMAACT, OPAFEI, ETAFEI, T, IVERS_FE_EXPFAC, 
     >                    TEFF)

C***********************************************************************
C***  NON-LTE IRON OPACITY AT GIVEN FREQUENCY FOR ALL DEPTH POINTS
C***  CALLED FROM: COLI and FORMAL -> FORMCMF
C***********************************************************************
 

      DIMENSION SIGMAFE(1), INDRB(LASTFE)
      DIMENSION IFRBSTA(LASTFE), IFRBEND(LASTFE)
      DIMENSION IFENUP(LASTFE), IFELOW(LASTFE)
      DIMENSION OPAFE(ND), ETAFE(ND)
      DIMENSION POPNUM(ND,N), WEIGHT(N), ENTOT(ND)
      DIMENSION INDFEACT(LASTFE), SIGMAACT(LASTFE)
      DIMENSION OPAFEI(ND,LASTFE), ETAFEI(ND,LASTFE)
      DIMENSION ELEVEL(N), T(ND)

C***  C1 = H * C / K    ( CM * ANGSTROEM )
      DATA C1 / 1.4388 /
C***  C2 = 2 * H * C     ( CGS UNITS )
      DATA C2 / 3.9724E-16 /

C***  CALCULATE FE-FREQUECY-INDEX
      XINDF = - ALOG10(XLAM0FE/XLAMK) /
     >          ALOG10(1. + VDOPFE*DXFE/CLIGHT)

C***  WAVENUMBERS IN KAYSER ( = CM**-1)
      WAVNUMK = 1.E8 / XLAMK
      WAVNUMK3 = WAVNUMK * WAVNUMK * WAVNUMK 

C***  PRESET OPAFE, ETAFE
      DO L=1, ND
         OPAFE(L) = 0.
         ETAFE(L) = 0.
      ENDDO
      
C***  LOOP OVER ACTIVE FE-LINES
      DO 10 INDACT=1, MAXFEACT

         IND = INDFEACT(INDACT)
         LOW = IFELOW(IND)
         NUP = IFENUP(IND)
         WEIGHTLU = WEIGHT(LOW) / WEIGHT(NUP)
         WAVNUM0 = ELEVEL(NUP) - ELEVEL(LOW)
         WAVNUM03 = WAVNUM0 * WAVNUM0 * WAVNUM0 

C***    FIND INDICES IN ARRAY 'SIGMAFE'
         INDF = NINT(XINDF-0.5)-IFRBSTA(IND)
         INDFSTART = INDRB(IND)
C***    READ NEIGHBOURING VALUES
         SIGMA1 = SIGMAFE(INDFSTART + INDF)
         SIGMA2 = SIGMAFE(INDFSTART + INDF + 1)
C***    INTERPOLATION
         DELTA = XINDF - FLOAT(INDF+IFRBSTA(IND))
         SIGMA = (1.-DELTA)*SIGMA1 + DELTA*SIGMA2

C***    MULTIPLY WITH RSTAR (DIMENIONS OF OPA AND ETA)
         SIGMAR = SIGMA * RSTAR

ccc     For test use: if "nu-factors" are activated:
ccc     modify SIGMA --> used in FREQUINT for integrartion of XJFEMEAN
ccc     -- the following statement might be de-activated else
ccc         SIGMA = SIGMA * WAVNUM0 / WAVNUMK
ccc    Note!!! SIGMAACT is NOT MANIPULATED in the Version with correct RRATE 
ccc            integration!!

C***    STORE SIGMA FOR LATER USE
         SIGMAACT(INDACT) = SIGMA

C***    CALCULATE 'OPA' AND 'ETA' OVER DEPTH INDEX 'L' 
         DO 20 L=1, ND
            ENUP  = POPNUM(L,NUP) * ENTOT(L)
            ENLOW = POPNUM(L,LOW) * ENTOT(L)
C***       Opacity and Emissivity
C***     Different versions according to the IRONLINES-EXPFAC option:
C***     0: OFF  - no exp factor (recommended by Goetz) 
C***     1: TEFF - exp factor with max(Te,Teff) (recommended by wrh, default) 
C***     2: TRADITIONAL - exp factor with Te (standard till 17-Feb-2018)
            IF (IVERS_FE_EXPFAC .EQ. 0) THEN
              G = WEIGHTLU 
            ELSEIF (IVERS_FE_EXPFAC .EQ. 1) THEN
              G = WEIGHTLU * EXP(C1*(WAVNUM0-WAVNUMK)/MAX(T(L),TEFF))
            ELSEIF (IVERS_FE_EXPFAC .EQ. 2) THEN
              G = WEIGHTLU * EXP(C1*(WAVNUM0-WAVNUMK)/T(L))
            ELSE 
              STOP '*** FATAL INTERNAL ERROR 1 in subr. CMFFEOP'
            ENDIF
            EMINDU = G * ENUP * SIGMAR
            SUM = ENLOW * SIGMAR - EMINDU
            OPAFEI(L, IND) = SUM
            IF (IVERS_FE_EXPFAC .EQ. 0) THEN
              ETAFEI(L, IND) = EMINDU * C2 * WAVNUM03
            ELSE
              ETAFEI(L, IND) = EMINDU * C2 * WAVNUMK3
            ENDIF
ccc         in the preceding stament WAVNUM03 has been used by Goetz

            IF (OPAFEI(L,IND) .GT. 0.) THEN
               OPAFE(L) = OPAFE(L) + OPAFEI(L,IND)
               ETAFE(L) = ETAFE(L) + ETAFEI(L,IND)
            ENDIF
 20      CONTINUE

 10   CONTINUE

      RETURN
      END
