      FUNCTION PHIHOLTSMARK(BETA, BETADOP, bHYDROGEN)
C***********************************************************************
C***  called from STARKHOLTZMARK
C***
C***  The profile is obtained by connecting a doppler core to an
C***   asymptotic Holtsmark function
C***
C***  BETA = DELTALAM / F_0 / K_low,up           
C***  BETADOP = c * DELTANUED / F_0 / K_low,up
C***
C***  Attention: BETA and BETADOP are not truly dimensionless, but
C***             instead need to be provided in Angstroem/cm
C***
C***
C***  concept taken from Hubeny et al. 1994, A&A 282, 151 (Appendix B)
C***********************************************************************

      IMPLICIT NONE
      
      LOGICAL, INTENT(IN) :: bHYDROGEN

      REAL :: PHIHOLTSMARK
      REAL :: BETADCRIT, FHOLTSMARK, X0, X1, DX, C, BETASTAR
      REAL, INTENT(INOUT) :: BETA
      REAL, INTENT(IN) :: BETADOP
      
      REAL, PARAMETER :: XEPS = 1.E-3
      
      REAL, PARAMETER :: BETADCRITH  = 5.82         !for Hydrogen only      
c      REAL, PARAMETER :: BETADCRITHE = 3.66         !for all other hydrogenic ions
      REAL, PARAMETER :: BETADCRITHE = 3.67         !for all other hydrogenic ions
      
C***  Edmonds fitting parameters      
      REAL, PARAMETER :: A0 = 0.07209481
      REAL, PARAMETER :: A1 = 0.4796232
      REAL, PARAMETER :: A2 = -0.5758228        !Attention: The minus sign is erronously missing in Hubeny et al (1994)

C***  Stark broadening: empirical branch limits
C      Original set from Hubeny et al. (1994)
c      REAL, PARAMETER :: BL1 = 1.14
c      REAL, PARAMETER :: BL2 = 11.4
C      Revised set from Hubeny's Synspec code
C       better connection between branches
      REAL, PARAMETER :: BL1 = 1.52
      REAL, PARAMETER :: BL2 = 8.325        
      
      REAL, PARAMETER :: WPI         = 1.77245385   !sqrt(Pi)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      
C***  Ensure that beta is always a positive wavelength difference      
      BETA = ABS(BETA)      
      IF (BETADOP < 0.) STOP 'FATAL ERROR: BETADOP < 0'
      
      IF (bHYDROGEN) THEN
        BETADCRIT = BETADCRITH
        FHOLTSMARK = 1.
      ELSE
        BETADCRIT = BETADCRITHE
        FHOLTSMARK = 1./2.
      ENDIF
           
           
C**********************************************************************            
      IF (BETADOP >= BETADCRIT) THEN     
C***    Doppler broadening is dominant in the line core
C***       => join Doppler profile with Holtsmark profile in the wing
C***          switch is done at BETA > BETASTAR
C***            with BETASTAR denoting the position where both profile
C***            functions have the same value. With x = BETA/BETADOP
C***            this leads to the Equation x^2 - 2.5 * ln(x) - C = 0

C***    Calculation of the constant C
        C = 1.5 * LOG(BETADOP) - LOG(3.*FHOLTSMARK*WPI)

C***    BETASTAR is determined via a short iteration using a
C***    Newton-Raphson approach to find x where
C***       x^2 - 2.5 * ln(x) - C = 0
        
C***    N-R starting value is chosen depending on C (taken from Hubeny)        
        IF (C > 1.26) THEN
          X0 = SQRT(C) * ( 1. + 1.25 * LOG(C) / (4.*C - 5) )
        ELSE
          X0 = SQRT(C + 0.28)
        ENDIF

        DO
          X1 = X0 - (X0*X0 - 2.5*LOG(X0) - C) / (2.*X0 - 2.5/X0)
          DX = 1. - X0/X1
          IF (ABS(DX) < XEPS) THEN
C***        iteration converged          
            BETASTAR = X1 * BETADOP
            EXIT
          ELSE
            X0 = X1
            IF (X0 <= 0) THEN
              WRITE (hCPR,*) 'PROBLEM WITH LINEAR STARK BROADENING'
              WRITE (hCPR,*) 'CANNOT DETERMINE CORE-WIND CONNECTION'
              WRITE (hCPR,*) 'FATAL (X0 < 0),  X0 = ', X0
              WRITE (hCPR,*) 'BETADOP > BETACRIT', BETADOP, BETADCRIT
              WRITE (hCPR,*) 'C, bHYDROGEN = ', C, bHYDROGEN
              STOP 'FATAL ERROR IN PHIHOLTSMARK'
            ENDIF 
          ENDIF
        ENDDO
C***    -- end of N-R iteration --        
        
        
        IF (BETA <= BETASTAR) THEN
C***      Doppler profile in the core
          PHIHOLTSMARK = 1./WPI/BETADOP * EXP( -1.* (BETA/BETADOP)**2 )
        ELSE
C***      Holtsmark profile in the wing
          PHIHOLTSMARK = 3. * FHOLTSMARK * BETA**(-5./2.)
        ENDIF

        
C**********************************************************************      
      ELSE
C***    BETADOP < BETADCRIT => Stark broadening is dominant everywhere
C                               (Doppler broadening is neglected)
c        WRITE (hCPR,*) 'Dominant Stark broadening'

C***    Fitting according to calculations of Edmonds et al. (1967)
        IF (BETA < BL1) THEN
          PHIHOLTSMARK = 0.08
        ELSEIF (BETA < BL2) THEN
          PHIHOLTSMARK = A0 * EXP(A1 * LOG(BETA) + A2 * (LOG(BETA))**2 )
        ELSE
          PHIHOLTSMARK = 3. * BETA**(-5./2.)
        ENDIF
C***    PHI needs to be scaled with FHOLTSMARK-factor (different for Hydrogen)        
        PHIHOLTSMARK = PHIHOLTSMARK * FHOLTSMARK

      ENDIF
C**********************************************************************      

C***  Profile normalization
C***  (except for depth-dependent VDOP, similar to VOIGTH)
      PHIHOLTSMARK = BETADOP * PHIHOLTSMARK 

      
C***  Fatal error for fatal results
C***  (indicating that something must be very wrong in the calling parameters)      
      IF (PHIHOLTSMARK < 0) STOP 'FATAL ERROR: PHIHOLTSMARK < 0'
      
      RETURN
      END
