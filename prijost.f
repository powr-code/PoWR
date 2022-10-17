      SUBROUTINE PRIJOST (NF,XLAMBDA,EMFLUX,DISMO, FWEIGHT)
C***********************************************************************
C***  PRINTING ABSOLUTE MAGNITUDES: JOHNSON, STROEMGREN etc.
C***  REFERENCES: Heber et al. (1984: A&A 130, 119)
C***  https://ui.adsabs.harvard.edu/abs/1984A%26A...130..119H/
C***  Revised on  3-Mar-2014: integration sum now over EMFLUX index
C***                          before it was over filter-function grid
C***  Revised on  6-Jun-2019 by wrh:
C***    Calibration constants changed (see
C***       ~wrh/science.dir/colorfilters.dir/calcon-Holber+Bergeron.plot
C***  New filters added: 2MASS J, H, Ks
C***  New filters added: Gaia Gbp, Grp, G; 
C***      their calibration constants FCAL from
C***      http://svo2.cab.inta-csic.es/svo/theory/fps3/
C***
C***  Definition of filters now outsourced to subr. FILTERFUNCTIONS 
C***    wrh 26-May-2021
C***********************************************************************

C***  NUMBER OF FILTERS
      INTEGER, PARAMETER :: NFILS = 13 
      REAL, DIMENSION(NFILS) :: FCAL, FCENTER

C***  Max. number of tabulated points in any of the filter functions 
      INTEGER, PARAMETER :: MAXFILPOINTS = 200 
      REAL, DIMENSION(MAXFILPOINTS) :: FLAM, FILT

      DIMENSION XLAMBDA(NF), EMFLUX(NF), FWEIGHT(NF)
      CHARACTER STRING*10

      CHARACTER*8 IBAND(NFILS)

      REAL, PARAMETER :: PI = 3.1415926535898 

C***  NAMES OF THE FILTERS:
C***  1 - 3: JOHNSON U,B,V; 4 - 7: STROEMGREN u,v,b,y; 8-10: 2MASS J, H, K
C***  11-13: Gaia Gbp, Grp, G
      DATA IBAND / 'JOHN U  ', 'JOHN B  ', 'JOHN V  ',
     -             'STROEM U', 'STROEM V', 'STROEM B', 'STROEM Y',
     -             'J_2MASS ', 'H_2MASS ', 'Ks_2MASS', 
     -             'Gaia_Gbp', 'Gaia_Grp', 'Gaia_G  ' /

C***  CENTER WAVELENGHTS OF THE BANDS (ONLY FOR INFORMAION) IN ANG
      DATA FCENTER / 3500., 4350., 5550., 3500., 4100., 4700., 5500., 
     >               12350., 16620., 21590.0, 5270., 7760., 6290. /

C***  CALIBRATION constants
C***  Commented: calibration constants before  6-Jun-2019
cc      DATA FCAL /  49.424, 48.525, 48.624, 
cc     -             48.324, 48.316, 48.452, 48.600 /

C***  Commented: calibration constants based on VO webpage
cc      DATA FCAL / 49.514, 48.488, 48.620,
cc     >            49.742, 48.463, 48.447, 48.606,
cc     >            49.490, 49.974, 50.442 /

C***  Based on Holberg & Bergeron (2006)
      DATA FCAL /  49.282, 48.395, 48.570,
     >             48.452, 48.656, 48.393, 48.587,
     >             49.503, 49.944, 50.461, 48.655, 49.024, 48.870 /


      PRINT 1
    1 FORMAT (//, 10X, 'ABSOLUTE MAGNITUDES IN FILTERS (JOHNSON, ', 
     $      'STROEMGREN, 2MASS, GAIA_EDR3)', /, 10X, 70('='), //,
     $      10X, 'BAND       WAVELENGTH    MAGNITUDE', /)

C***  INTEGRATION OF FILTER AND FILTER*EMFLUX
C***  --> CALCULATION OF MAGNITUDES

C***  Loop over the filters
      DO 102 J=1, NFILS
         CALL FILTERFUNCTIONS (IBAND(J), NFILT, FLAM, FILT)

C***     Quadrature sum; flux and filter norm are both integrated over dnu
         FFLUX = .0
         FNORM = .0
         NPOINTS = 0

C***     Loop over intensity's frequency points
         DO 101 K=1, NF
C***        Interpolate value of filter function
C***        Filter wavelength range
            FLAM1 = FLAM(1)
            FLAM2 = FLAM(NFILT)
            XLAM = XLAMBDA(K)
            IF ( (XLAM-FLAM1)*(XLAM-FLAM2) .LT. .0) THEN
               CALL SPLINPO (FILFUN, XLAM, FILT, FLAM, NFILT) 
               FFLUX = FFLUX + FILFUN * EMFLUX(K) * FWEIGHT(K)   
               FNORM = FNORM + FILFUN * FWEIGHT(K)
               NPOINTS = NPOINTS + 1
            ENDIF
  101    CONTINUE

C***     Normalization of that filter
         IF (FNORM .GT. .0) FFLUX = FFLUX / FNORM

C***     Special branch in (unlikely) case of too few points within filter
C***     Take flux at center wavelength
         IF (NPOINTS .LE. 1) THEN
            CALL SPLINPO (FFLUX, FCENTER(J), EMFLUX, XLAMBDA, NF)
         ENDIF

         IF (FFLUX .GT. .0) THEN
            COLOR = -2.5 * ALOG10(FFLUX*PI) + DISMO - FCAL(J)
            WRITE (STRING, '(F10.3)') COLOR
         ELSE
            STRING = '  * ZERO *'
         ENDIF

         PRINT 3, IBAND(J), FCENTER(J), STRING
    3    FORMAT (10X, A8, F10.0, ' A ', 1X,  A)

  102 CONTINUE

      PRINT ('(/)')

      RETURN
      END
