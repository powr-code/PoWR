      SUBROUTINE HYSTRUKU(RSTART, REND, DELTAB, 
     >                    RADIUS, GEFFL, A2SUM, ND, H0, RSTAR)
C***********************************************************************
C***  simple Runge-Kutta iteration for the differential equation
C***  which has to be solved to obtain the solution of the transformed
C***  hydrostatic equation
C***  
C***  called by VELTHIN  (for each depth point!)
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
     
C***  Steps for Runge-Kutta integration     
      INTEGER, PARAMETER :: IRUKUSTEP = 10  
     
      INTEGER, INTENT(IN) :: ND
      REAL, INTENT(IN) :: H0, RSTAR, RSTART, REND
      REAL, INTENT(OUT) :: DELTAB
      REAL, DIMENSION(ND), INTENT(IN) :: RADIUS, GEFFL, A2SUM
     
      REAL :: DR, DRHALF, DELTAR, RIN, RL, X1, X2, X3, X4,
     >        H, GEFFI, A2SUMI, DBDR
     
      INTEGER :: I
     
      DELTAR = REND - RSTART 
      
      DR = DELTAR / IRUKUSTEP
      DRHALF = DR / 2.

      RL = RSTART
      DELTAB = 0.
         
      DO I = 1, IRUKUSTEP
        RIN = RL        
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        H = A2SUMI / (GEFFI / RIN / RIN) / RSTAR
        DBDR = (1./H0 - 1./H)
        X1 = DR * DBDR

        RIN = RL + DRHALF
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        H = A2SUMI / (GEFFI / RIN / RIN) / RSTAR
        DBDR = (1./H0 - 1./H)
        X2 = DR * DBDR

        RIN = RL + DR
        IF (RIN > RADIUS(1)) THEN
C***      In the rare case that this routine runs up to
C***      the outer boundary, we have to avoid that due to
C***      numerical reasons RIN is slightly larger than RMAX
C***      and thus SPLINPOX calls would fail.
          A2SUMI = A2SUM(1)
          GEFFI = GEFFL(1)
        ELSE
          CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND)
          CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        ENDIF
        H = A2SUMI / (GEFFI / RIN / RIN) / RSTAR
        DBDR = (1./H0 - 1./H)
        X3 = DR * DBDR

        RL = RL + DR
        DELTAB = DELTAB + ((X1/2.) + 2. * X2 + (X3/2.)) / 3.
      ENDDO

      RETURN
      END
