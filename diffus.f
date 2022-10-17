      SUBROUTINE DIFFUS (XLAM,T,RADIUS,ND,BCORE,DBDR,DTDR,TEFF,NOTEMP)
C***********************************************************************
C***  CALLED FROM: WRCONT, ETL
C***  GIVES THE PLANCK FUNCTION, BCORE, AND ITS RADIUS-DERIVATIVE, DBDR, AT
C***  THE INNER BOUNDARY FROM ...
C***  IF (NOTEMP) ... THE GIVEN TEMPERATURE STRATIFICATION
C***  IF (.NOT. NOTEMP) ... "DTDR", THE DERIVATIVE OF THE TEMPERATURE WITH
C***                    RESPECT TO THE RADIUS, AS CALCULATED IN SUBR. DIFDTDR
C***                    TO YIELD THE CORRECT TOTAL FLUX  H = L/(4*PI*RSTAR)**2
C***  IN DIFFUSION APPROXIMATION, THEN THE INCIDENT INTENSITY WILL BE
C***      IPLUS = BCORE + DBDR * Z / X
C***  Z = MUE, X = OPACITY
C***********************************************************************
 
      DIMENSION T(ND),RADIUS(ND)
      LOGICAL NOTEMP
 
C***  Prevent Zeros in the radiation field XJC in case of 
c     small (X-ray) wavelength together with low temperatures
      BCORE=MAX(BNUE(XLAM,T(ND)),1.E-286)
 
      IF (NOTEMP) THEN
         BNDM=BNUE(XLAM,T(ND-1))
         DBDR=(BCORE-BNDM)/(RADIUS(ND-1)-1.)
      ELSE
C***     DERIVATIVE OF THE PLANCK FUNCTION WITH RESPECT TO T
         DBDT=DBNUEDT(XLAM,T(ND))
         DBDR=DBDT*DTDR
      ENDIF
 
      RETURN
      END
