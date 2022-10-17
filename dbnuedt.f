      FUNCTION DBNUEDT (XLAMBDA,T)
C**********************************************************************
C***  DERIVATIVE OF PLANCK FUNCTION WITH RESPECT TO TEMPERATURE
C***  XLAMBDA IN ANGSTROEM, T IN KELVIN  
C***  BNUE IN CGS UNITS: ERG PER CM**2, PER SEC AND PER HERTZ      
C**********************************************************************

C***  CONSTANTS :  C1 = H * C / K   (DIMENSION ANGSTROEM * KELVIN )
C***               C2 = 2 * H * C
      DATA C1,C2 / 1.4388E8, 3.9724E+8 /

C***  PREVENT ERROR FOR NEGATIVE TEMPERATURES 
      IF (T .LE. .0) THEN
         DBNUEDT=.0
         ELSE
         HNUEKT=C1/XLAMBDA/T
C*** PREVENT ERROR AR004 (BAD SCALAR ARGUMENT TO ARLIB MATH ROUTINE)
C*** FOR EXTREME WIEN-DOMAIN
C***                                AND FOR RALEIGH-JEANS DOMAIN

C***   CRAY VERSION
C!!!         IF (HNUEKT .GT. 5000. .OR. HNUEKT .LT. 1.E-6) THEN
C***   DEC-VERSION
         IF (HNUEKT .GT. 700. .OR. HNUEKT .LT. 1.E-6) THEN
            DBNUEDT=0.
         ELSE
            X2=XLAMBDA*XLAMBDA
            X4=X2*X2
            T2=T*T
            EXFAC=EXP(HNUEKT)
            DBNUEDT=C1*C2/X4/T2/(EXFAC + 1./EXFAC - 2.)
            ENDIF
         ENDIF

      RETURN
      END
