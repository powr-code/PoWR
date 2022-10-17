      SUBROUTINE HYSTHDRUKU(RSTART, REND, VSTART, VEND, 
     >                      RADIUS, GEFFL, A2SUM, ND, RSTAR)
C***********************************************************************
C***  4th-order Runge-Kutta iteration for the equation of motion
C***  which has to be solved to obtain the velocity field in the
C***  quasi-hydrostatic regime
C***  
C***  called by VELTHIN  (for each depth point!)
C***********************************************************************

      IMPLICIT NONE
      INCLUDE 'interfacebib.inc'
     
C***  Steps for Runge-Kutta integration     
      INTEGER, PARAMETER :: IRUKUSTEP = 50  
     
      INTEGER, INTENT(IN) :: ND
      REAL, INTENT(IN) :: RSTAR, RSTART, REND, VSTART
      REAL, INTENT(OUT) :: VEND
      REAL, DIMENSION(ND), INTENT(IN) :: RADIUS, GEFFL, A2SUM
     
      REAL :: DR, DRHALF, DELTAR, RIN, RL, X1, X2, X3, X4,
     >        H, GEFFI, A2SUMI, DVDR, DA2DR, VIN, VL, R
     
      INTEGER :: I
     
      DELTAR = REND - RSTART 
      
      DR = DELTAR / IRUKUSTEP
      DRHALF = DR / 2.

      RL = RSTART
      VL = VSTART * 1.E5
                 
      DO I = 1, IRUKUSTEP

        RIN = RL      
        VIN = VL
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND, DFDX=DA2DR)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        R = RIN * RSTAR
        DA2DR = DA2DR / RSTAR
        DVDR = (GEFFI/RIN/RIN - 2.*A2SUMI/R + DA2DR) / 
     >                (A2SUMI/VIN - VIN)   
        DVDR = DVDR * RSTAR
        X1 = DR * DVDR

        RIN = RL + DRHALF
        VIN = VL + X1/2.
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND, DFDX=DA2DR)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        R = RIN * RSTAR
        DA2DR = DA2DR / RSTAR
        DVDR = (GEFFI/RIN/RIN - 2.*A2SUMI/R + DA2DR) / 
     >                (A2SUMI/VIN - VIN)   
        DVDR = DVDR * RSTAR
        X2 = DR * DVDR

        RIN = RL + DRHALF
        VIN = VL + X2/2.
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND, DFDX=DA2DR)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        R = RIN * RSTAR
        DA2DR = DA2DR / RSTAR
        DVDR =  (GEFFI/RIN/RIN - 2.*A2SUMI/R + DA2DR) / 
     >                (A2SUMI/VIN - VIN)   
        DVDR = DVDR * RSTAR
        X3 = DR * DVDR

        RIN = RL + DR
        VIN = VL + X3
        CALL SPLINPOX(A2SUMI, RIN, A2SUM, RADIUS, ND, DFDX=DA2DR)
        CALL SPLINPOX(GEFFI,  RIN, GEFFL, RADIUS, ND)
        R = RIN * RSTAR
        DA2DR = DA2DR / RSTAR
        DVDR =  (GEFFI/RIN/RIN - 2.*A2SUMI/R + DA2DR) / 
     >                (A2SUMI/VIN - VIN)   
        DVDR = DVDR * RSTAR
        X4 = DR * DVDR

        RL = RL + DR
        VL = VL + ((X1/2.) + X2 + X3 + (X4/2.)) / 3.
      ENDDO
      
      VEND = VL / 1.E5

      RETURN
      END
