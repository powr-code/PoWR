      SUBROUTINE MGOETZ (XMSTAR, TEFF, RSTAR, XHY, WRTYPE)     
C***********************************************************************
C***  Calculates Stellar Mass from the Luminosity (mass-luminosity relation)
C***  according to Graefener et al. (2011: A&A 535, 56)
C***  Notes: 
C***  - masses are for chemically homogeneous stars, and thus upper limits
C***  - two different M-L relations: 
C***       WRTYPE='WN': He-burning stars (H-free)
C***       WRTYPE='OB': Hydrogen-burning stars
C***  - resulting mass XMSTAR is in units of MSUN  
C***********************************************************************

      IMPLICIT NONE

      REAL, INTENT(OUT) :: XMSTAR
      REAL, INTENT(IN) :: TEFF, RSTAR, XHY
      
      REAL :: F1, F2, F3, F4, F5, F6, F7, F8, F9, F1r, F2r, F3r, F4r,
     >        XLSTAR, XLSTARS, XLOGL, XLOGM, Radikant
      
      !Constants
      REAL, PARAMETER :: STEBOL = 5.6705E-5   !STEBOL = STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      REAL, PARAMETER :: PI4 = 12.5663706144  !PI4 = 4*PI
      REAL, PARAMETER :: XLSUN = 3.85E33      !Solar Luminosity (CGS-Units)
      REAL, PARAMETER :: XMSUN = 1.989E33     !Solar Mass (g)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)

      CHARACTER WRTYPE*2

C***  Fit parameters for hydrogen-burning stars, Eqs. 11 and 12, Table A.1  
      F1=4.026 
      F2=4.277
      F3=-1.0
      F4=25.48
      F5=36.93
      F6=-2.792
      F7=-3.226
      F8=-5.317
      F9=1.648

C***  Fit parameters for helium-burning stars, Eq. 13, Table A.1
      F1r=3.997
      F2r=-1.0
      F3r=25.83
      F4r=-3.268

            
      XLSTAR = PI4 * STEBOL * RSTAR*RSTAR * TEFF*TEFF*TEFF*TEFF
      XLSTARS = XLSTAR / XLSUN
      XLOGL = ALOG10(XLSTARS)

C***  Hydrogen-burning star 
      IF (WRTYPE == 'OB') THEN
        IF (XHY .LT. 0.1) GOTO 90
        Radikant = F4 + F5 * XHY + F6 * XHY**2 + (F7+ F8*XHY) * XLOGL
        XLOGM = ( F1+F2*XHY + F3*SQRT( Radikant ) ) / ( 1.+F9*XHY )

C***  Helium-burning star
      ELSEIF (WRTYPE == 'WN') THEN
        Radikant = F3r + F4r * XLOGL
        XLOGM = F1r + F2r * SQRT( Radikant )
      ELSE
        GOTO 91
      ENDIF
      
      XMSTAR = 10.**XLOGM

      RETURN

C***  ERROR BRANCHES
   90 WRITE (0,*) '*** ERROR: M-L relation only valid for X_H > 0.1' 
      GOTO 99

   91 WRITE (0,*) '*** ERROR: M-L relation (Graefener version) not '
     >             // 'valid for type ' // WRTYPE 

   99 STOP '*** INTERNAL ERROR in Subr. MGOETZ'

      END
