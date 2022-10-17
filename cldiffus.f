      SUBROUTINE CLDIFFUS (XLAM,T,RADIUS,ND,BCORE,DBDR,DTDR,TEFF,
     >                     OPAK,XHID)
C***********************************************************************
C***  CALLED FROM: WRCONT, ETL
C***  GIVES THE PLANCK FUNCTION, BCORE, AND ITS RADIUS-DERIVATIVE, DBDR, AT
C***  THE INNER BOUNDARY FROM ...
C***               ... "DTDR", THE DERIVATIVE OF THE TEMPERATURE WITH
C***                    RESPECT TO THE RADIUS, AS CALCULATED IN SUBR. DIFDTDR
C***                    TO YIELD THE CORRECT TOTAL FLUX  H = L/(4*PI*RSTAR)**2
C***  IN DIFFUSION APPROXIMATION, THEN THE INCIDENT INTENSITY WILL BE
C***      IPLUS = BCORE + DBDR * Z / X
C***  Z = MUE, X = OPACITY
C***********************************************************************
 
      DIMENSION T(ND),RADIUS(ND),OPAK(ND)
 
      BCORE=BNUE(XLAM,T(ND))
 
C***  DERIVATIVE OF THE PLANCK FUNCTION WITH RESPECT TO T
      DBDT=DBNUEDT(XLAM,T(ND))
      DBDR=DBDT*DTDR

C***  Calculate the corresponding H at the inner Boundary
      XHID = DBDR/OPAK(ND)/3.

      RETURN
      END
