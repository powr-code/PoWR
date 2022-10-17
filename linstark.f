      SUBROUTINE LINSTARK (GRIEMPAR, XNE, NUCCHARG, 
     >                     MAINQNLOW, MAINQNNUP, XLAM0,
     >                     LINPRO, LEVELLOW, LEVELNUP)
C**********************************************************************
C***  Linear Stark effect, prepartion routine for GRIEMPAR vector
C***     motivated by the approach from Hubeny et al. (1994)
C***  
C***  GRIEMPAR             - conversion factor for Griem's beta parameter
C***                               = lamda_0/(F_0 * K_low,up)
C***  XNE                  - electron density
C***  NUCCHARG             - nuclear charge 
C***                           (often misleadingly labeled ion charge)
C***  MAINQNLOW, MAINQNNUP - main quantum numbers of lower and upper level 
C***  XLAM0                - transition wavelength in Angstroem
C***
C***  returns GRIEMPAR in unit of:  Angstroem/cm
C***    (i.e. needs to be mulitplied with _relative_ frequency difference)
C***
C***  written by ansander (Dec 2015)
C***
C***  called from STARKBROAD
C**********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MAINQNLOW, MAINQNNUP, NUCCHARG
      
      CHARACTER(8) :: LINPRO
      CHARACTER(10) :: LEVELLOW, LEVELNUP
      CHARACTER(40) :: CQNERRDETAIL

      REAL, EXTERNAL :: KHOLTSMARK
      
      REAL :: GRIEMPAR, F0, T, XNE, XLAM0, HOLTSK
            
C***  Holtsmark field strength normalization factor      
C***  1.25 = 2 * Pi * (4/15)^(2/3) * e 
      REAL, PARAMETER :: FNORM = 1.25E-9            ! in cm^(3/2) * g^(1/2) * s^(-1)
      

C***  The main quantum numbers of both levels need to be known for this method!
      IF (MAINQNLOW <= 0 .AND. MAINQNNUP <= 0) THEN
        CQNERRDETAIL = '  (both main quantum numbers not known) '
      ELSEIF (MAINQNLOW <= 0) THEN
        CQNERRDETAIL = '  (lower quantum number not known)      '
      ELSEIF (MAINQNNUP <= 0) THEN
        CQNERRDETAIL = '  (upper quantum number not known)      '
      ENDIF

C***  Fallback to quadratic Stark effect in case of unknown quantum numbers      
      IF (MAINQNLOW <= 0 .OR. MAINQNNUP <= 0) THEN
           LINPRO = 'Q-STARK '
           GRIEMPAR = .0
           WRITE (0,'(5A)') '*** WARNING: LINSTARK cannot handle ', 
     >           LEVELLOW, ' - ', LEVELNUP,  CQNERRDETAIL
           WRITE (*,'(5A)') '*** WARNING: LINSTARK cannot handle ', 
     >           LEVELLOW, ' - ', LEVELNUP,  CQNERRDETAIL
           RETURN
      ENDIF 

C***  F0 = Holtsmark normal field strength
      F0 = FNORM * XNE**(2./3.)
C***  HOLTSK is Griems K_low,up correction factor, generalized with
C***    a charge factor to approximate for all hydrogenic ions
      HOLTSK = KHOLTSMARK(MAINQNLOW, MAINQNNUP, NUCCHARG)

C***  GRIEMPAR contains the transition-dependent quantities which can
C***  be precalculated before the wavelength integration
C***  The factor XLAM0 in the nominator is neeed to convert the relative wavelengths used
C***  in the integration in ZONEINT back into absolute wavelengths in Angstroem which
C***  are required by the PHISTARK routine.
      GRIEMPAR = XLAM0 / (F0 * HOLTSK)
      
              
      RETURN
      END
