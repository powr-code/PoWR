      SUBROUTINE PRIFLUX (NF,XLAMBDA,EMFLUX,TOTIN,TOTOUT,RSTAR,JOBNUM,
     >          FWEIGHT,MODHEAD,KEY,EMCOLI)
C***********************************************************************
C***   PRINTOUT OF THE EMERGENT CONTINUUM FLUX
C***********************************************************************
 
      DIMENSION XLAMBDA(NF),EMFLUX(NF),FWEIGHT(NF),KEY(NF), EMCOLI(NF)
      CHARACTER*30 MODHEAD*100
      CHARACTER*8 ITCOL, IFSUM, JUMP
 
      LOGICAL BEMCOLI

C***  CLIGHT = SPEED OF LIGHT ( CM/SEC )
      DATA CLIGHT / 2.99792458E10 /
C***  STEBOL = STEFAN-BOLTZMANN CONSTANT / PI  (ERG/CM**2/S/STERAD/KELVIN**4)
      DATA STEBOL / 1.8046E-5 /
 
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT (1X,  A,  20X,'JOB NO.',I7,
     > //,30X,
     > 'E M E R G E N T  A S T R O P H Y S I C A L  F L U X',/,30X
     > ,51('='),/,34X,'( FLUX = PI * ASTROPHYS.FLUX = 4 * PI * H )',//,
     > ' FREQUENCY      LAMBDA    LOG F-NUE   LOG F-LAMBDA  ',
     > '   T-RAD     T-COLOR    F-NUE /     LOG INTEGRATED F-NUE ',
     > '  CONT. JUMP',
     > /,'   INDEX     (ANGSTROEM)  (ERG/CM+2) (ERG/CM+2/S/A) ',
     > ' (KELVIN)    (KELVIN)    B(TEFF)    (ERG/CM2/SEC)    ',
     > '    (MAGNITUDES)',/)
 
C***  CALCULATION OF EFFECTIVE TEMPERATURE
      TEFF=(TOTOUT /STEBOL)**.25
C***  CALCULATION OF LUMINOSITY
C***  THE NUMERICAL CONSTANT DENOTES THE LOG OF 4*PI*STEBOL/LSOLAR
      ALUMI=2*ALOG10(RSTAR*TEFF*TEFF)-36.730
 
      FSUM=.0
C***  LOOP OVER ALL FREQUENCY POINTS  ----------------------------------
      DO 2 K=1,NF
      XLAM=XLAMBDA(K)
      FK=EMFLUX(K)
      W=1.E8/XLAM
      IF (FK .LE. .0) THEN
            FNUE=.0
            FLAMBDA=.0
            ELSE
            FNUE=ALOG10(FK)
            FLAMBDA=FNUE+ALOG10(CLIGHT*W/XLAM)
            ENDIF
      TRAD=TRADFUN(XLAM,FK)
C***  CALCULATION OF TCOLOR FROM THE SLOPE OF EMFLUX BETWEEN K-1 AND K
      ITCOL='        '
      IF (K.GT.1) THEN
         CALL TCOLOR (XLAMBDA(K-1),XLAM,EMFLUX(K-1),FK,ITCOL)
         ENDIF
C***  Avoid zero Planck Function, Lars 18-Mar-1997 10:28:42
      BN = BNUE(XLAM,TEFF)
      IF (BN .GT. 0) THEN
        FOVERB=FK/BN
      ELSE
        FOVERB = UNDEF
      ENDIF
C***  INTEGRATION OF F-NUE FROM ZERO TO CURRENT WAVELENGTH
      FSUM=FSUM+FK*FWEIGHT(K)
      IF (FSUM .GT. .0) THEN
         FSUMLOG=ALOG10(FSUM)
         ENCODE (8,5,IFSUM) FSUMLOG
    5    FORMAT (F8.3)
         ELSE
         IFSUM='    -INF'
         ENDIF
C***  CONTINUUM JUMPS AT EDGE FREQUENCIES
      JUMP=' '
      DECODE (5,99,KEY(K)) NEDGE
   99 FORMAT (A5)
      IF (NEDGE .EQ. 5HEDGE+ .OR. NEDGE .EQ. 5HK-ED+) THEN
         IF (EMFLUX(K-1) .LE. .0 .OR. EMFLUX(K) .LE. .0) THEN
            JUMP='   UNDEF'
            ELSE
            ENCODE (8,8,JUMP) ALOG10(EMFLUX(K)/EMFLUX(K-1))*2.5
    8       FORMAT (F8.3)
            ENDIF
         ENDIF
    2 PRINT 3, K,XLAM,FNUE,FLAMBDA,TRAD,ITCOL,FOVERB,IFSUM,JUMP
    3 FORMAT (I7,F15.2,2F13.3,F12.0,4X,A8,G14.3,6X,A8,10X,A8)
C***  ENDLOOP  ---------------------------------------------------------
 
      PRINT 4, TEFF, ALUMI 
ccccccccccc,TOTIN/TOTOUT
    4 FORMAT (/,10X,'EFFECTIVE TEMPERATURE:',F8.0,' KELVIN ',
     >          10X,'LOG OF LUMINOSITY (SOLAR UNITS):',F6.2,/)
ccc     >          10X,'RADIATIVE ENERGY INPUT/OUTPUT:',F8.3,/,/)
 
C***  CALCULATE AND PRINT No. OF UV-PHOTONS AND ZANSTRA TEMPERATURES
      PRINT *
      PRINT *, '-----------------------------------------------'
      PRINT *, 'Emergent Flux from WRCONT:'
      CALL ZANSTRA (NF, EMFLUX, XLAMBDA, FWEIGHT, TEFF, RSTAR)
      PRINT *
      PRINT *, '-----------------------------------------------'
      BEMCOLI = .FALSE.
      DO L=1, NF
         IF (EMCOLI(L).NE.0.) BEMCOLI = .TRUE.
      ENDDO 
      IF (BEMCOLI) THEN
         PRINT *, 'Emergent Flux from COLI:'
         CALL ZANSTRA (NF, EMCOLI, XLAMBDA, FWEIGHT, TEFF, RSTAR)
         CALL PRICOLR (NF,XLAMBDA,EMCOLI,RSTAR,FWEIGHT)
      ELSE
         PRINT *, 'NO FLUXES FROM COLI FOUND IN MODEL FILE !'
         CALL PRICOLR (NF,XLAMBDA,EMFLUX,RSTAR,FWEIGHT)
      ENDIF
 
 
      RETURN
      END
