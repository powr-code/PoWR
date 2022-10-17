      SUBROUTINE MLANGER (XMSTAR, TEFF, RSTAR, YHE, WRTYPE) 
C***  Calculates Stellar Mass for a H-free Star
C***  from the Luminosity (mass-luminosity relations)
C***  Note: the resulting mass XMSTAR is given in gramms !!  

      CHARACTER*2 WRTYPE
      
C***  STEBOL = STEFAN-BOLTZMANN CONSTANT (CGS-UNITS)
      DATA STEBOL / 5.6705E-5 /
C***  PI4 = 4*PI
      DATA PI4 / 12.5663706144 /
C***  XLSUN = Solar Luminosity (CGS-Units)
      DATA XLSUN / 3.85E33 /
C***  XMSUN = Solar Mass (g)
      DATA XMSUN / 1.989E33 /

      XLSTAR = PI4 * STEBOL * RSTAR*RSTAR * TEFF*TEFF*TEFF*TEFF
      XLSTARS = XLSTAR / XLSUN

C***  Mass-Luminosity relation from Langer N., A&A 210, 93 (1989)


      XLOGL = ALOG10(XLSTARS)

      IF (WRTYPE .EQ. 'WC') THEN
         A0 = -0.487870 * YHE
         A0 = A0 + 2.971463
         A1 = +0.434909 * YHE
         A1 = A1 + 2.771634
         A2 = -0.093793 * YHE
         A2 = A2 - 0.487209
         P = A1/A2
         Q = (A0 - XLOGL)/A2
         PQS = P*P/4. - Q
         IF (PQS .GE. 0.) THEN
cc            XLOGMP = -P/2. + SQRT(PQS)
            XLOGMM = -P/2. - SQRT(PQS)
            XMSTARS = 10.**XLOGMM
            XMSTAR = XMSTARS * XMSUN
         ELSE
            WRITE (0,*) '*** ERROR: Langer-Formula failed!'
            STOP '*** ERROR IN SUBROUTINE MLANGER'
         ENDIF
      ELSE IF (WRTYPE .EQ. 'WN') THEN
         XLOGM = - 0.158206 - 0.053868*XLOGL + 0.055467*XLOGL*XLOGL
         XMSTARS = 10.**XLOGM
         XMSTAR = XMSTARS * XMSUN
      ELSE
         WRITE (0,*) '*** ERROR: ILLEGAL WRTYPE: ', WRTYPE
         STOP '*** ERROR IN SUBROUTINE MLANGER'
      ENDIF

c      WRITE (0, *) '*** MASS FROM LUMINOSITY ***'
C      WRITE (0, *) XLOGMP, XLOGMM
C      WRITE (0, *) XLSTARS, 10**(A0 + A1*XLOGMM + A2*XLOGMM*XLOGMM)
c      WRITE (0, *) 'L=', XLSTARS, 'L_SUN'
c      WRITE (0, *) 'M=', XMSTARS, 'M_SUN'
      

C***  Zeta(M) correcture from Heger A., A&A 315, 421

C      ZETA = 1.8548 - 0.039711*MSTAR +6.1946E-4*MSTAR*MSTAR


C***  Mass-Radius relation from Langer N., A&A 210, 93 (1989)

C      B0 = -0.114128 * YHE
C      B0 = B0 - 0.497381
C      B1 = +0.208654 * YHE
C      B1 = B1 + 0.371859
C      B2 = -0.156504 * YHE
C      B2 = B2 + 0.156504

C      RTH = B0 + B1*XLOGM + B2*XLOGM2

      RETURN
      END
