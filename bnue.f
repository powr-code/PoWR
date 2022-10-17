      FUNCTION BNUE(XLAMBDA,T)
C***********************************************************************
C***  PLANCK FUNCTION, LAMBDA IN ANGSTROEM, T IN KELVIN
C***  BNUE IN CGS UNITS: ERG PER CM**2, PER SEC AND PER HERTZ
C***  CONSTANTS :  C1 = H * C / K   (DIMENSION ANGSTROEM * KELVIN )
C***               C2 = 2 * H * C
C***********************************************************************
      DATA C1,C2 / 1.4388E8, 3.9724E+8 /
C***      DATA EXPMAX / 650. /
C***  geaendert wrh  5-Jul-2013 
      DATA EXPMAX / 500. /

C***  PREVENT ERROR FOR NEGATIVE TEMPERATURES
      IF (T .LE. .0) THEN
         BNUE=.0
      ELSE
         W=1./XLAMBDA
         HNUEKT=C1*W/T
C***  PREVENT ERROR (BAD SCALAR ARGUMENT TO ARLIB MATH ROUTINE)
C***  FOR EXTREME WIEN-DOMAIN
         IF (HNUEKT .GT. EXPMAX) THEN
           FAC=C2*W*W*W
           FACLOG=ALOG(FAC)
           ARG= HNUEKT - FACLOG
           IF (ARG .GT. EXPMAX) THEN
            BNUE = EXP(-EXPMAX)
           ELSE  
            BNUE = EXP(-ARG)
           ENDIF  
         ELSE IF (HNUEKT .LT. 1.E-10) THEN
            BNUE=C2/HNUEKT*W*W*W
         ELSE
            BNUE=C2/(EXP(HNUEKT)-1.)*W*W*W
         ENDIF
      ENDIF

      RETURN
      END
