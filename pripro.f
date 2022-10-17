      SUBROUTINE PRIPRO (XLAM,VDOP,NFOBS,PROFILE,XOBS0,DXOBS,JOBNUM,
     $        VSINI,MODHEAD,DLAM,LSPRO,IFIRST,NPHI,LPSTA,LPEND,DXMAX,
     $        TAUMAX,XMAX,TAUMINBROAD,JFIRST,JLAST,
     >        P, PWEIGHT, ELEVEL,NDIM,INDNUP,INDLOW,
     $        LASTIND,INDLAP,MAXLAP,FREMAX,FREMIN,NBLINE,RMAX,XLAMLAP,
     $        WEIGHT,LEVEL,EINST,VMAX,IVERSION,LINE, REDIS,
     $        BWESEX,BROAD,LINPRO,AVOIGT, BNOCONT, NMOD, 
     >        DENSCON, FILLFAC, 
     >        BFEWARNING, SPZ1, SPZ2, MANIPOP_OPTIONS, ND,
     >        MULTIIND, NMULTI, BAIRWAVELENGTH, MAXSPZ, RCOROT,
     >        BDD_VDOP, DD_VDOP_LINE, BMICROTURB, EL_MIN, 
     >        IVDOPSTATUS, BVSINI_AT_RCOROT)
C***********************************************************************
C***  PRINTOUT OF THE EMERGENT LINE PROFILES
C***********************************************************************

      DIMENSION PROFILE(NFOBS),DLAM(NFOBS),P(JFIRST)
      DIMENSION ELEVEL(NDIM), WEIGHT(NDIM), EINST(NDIM,NDIM)
      DIMENSION INDNUP(LASTIND), INDLOW(LASTIND), MULTIIND(LASTIND)
      DIMENSION INDLAP(NBLINE), XLAMLAP(NBLINE)
      DIMENSION AVOIGT(NBLINE)
      DIMENSION JOBNUM(NMOD)

      LOGICAL SPIKE, REDIS, BROAD, BNOCONT, BFEWARNING, BAIRWAVELENGTH
      LOGICAL BDD_VDOP, BMICROTURB, BVSINI_AT_RCOROT
      CHARACTER LEVEL(NDIM)*10, OVERLAP*30, SUBLINE*1, EL_MIN*2
      CHARACTER*100 MODHEAD(NMOD), ACTPAR
      CHARACTER*(*) DD_VDOP_LINE
      CHARACTER LINPRO(NBLINE)*8, AVOIGTSTR*8
      CHARACTER REDIMES*30, BROAMES*6
      CHARACTER TXHALF*8, TXKM*8
      CHARACTER*(*) SPZ1(MAXSPZ), SPZ2(MAXSPZ), MANIPOP_OPTIONS

      DIMENSION DENSCON(ND,NMOD),FILLFAC(ND,NMOD)

      INTEGER :: KSPIKE
      REAL :: SPIKELAM

C***  CLIGHT = VELOCITY OF LIGHT IN KM/SEC
      DATA CLIGHT /2.99792458E+5/
C***  C4 = 1 / (SIGMA-CLASSIC * 8 PI)     (IN ANGSTROEM**2)
      DATA C4 /1.499E-16/
C***  CFACS: ALL CONSTANT FACTORS USED IN THE CALCULATION OF THE ABSOLUTE
C***  LINE EMISSION (ABLIEM) IN UNITS OF THE SOLAR LUMINOSITY LSUN=3.82E33 5I7/S
C***  CFACS = 1.E8*PI*4.*PI*C/LSUN
      DATA CFACS /3.098E-14/
      DATA PI / 3.14159 /

C***  logarithmic doppler unit (here negative!
      ALN = ALOG (1. - VDOP / CLIGHT)

 
C***  OUTPUT OF HEADER: ONE TIME (VARIABLE "IFIRST" IS SET TO ZERO AFTERWARDS)
      IF (IFIRST .EQ. 1 .OR. LSPRO .GT. 0 .OR. NMOD .GT. 1) THEN
        WRITE (*,*)
        DO IMOD=1, NMOD
          PRINT 4, MODHEAD(IMOD)( :IDX(MODHEAD(IMOD))),JOBNUM(IMOD)
        ENDDO
    4   FORMAT (1X, A, 4X, 'AFTER JOB NO.', I7)
      ENDIF

      XLAMMAX = XLAM * EXP(ALN*FREMIN)
      XLAMMIN = XLAM * EXP(ALN*FREMAX)
      WRITE (*,'(/,A,2F10.2)') 
     >     ' WAVELENGTH RANGE [Angstroem]:', XLAMMIN, XLAMMAX

      IF ( BAIRWAVELENGTH ) THEN
       WRITE (*,'(A,/)') ' WAVELENGTHS REFER TO AIR'
      ELSE
       WRITE (*,'(A,/)') ' WAVELENGTHS REFER TO VACUUM'
      ENDIF


      WRITE (*,'(A,I12)')
     >      ' Total number of lines considered = ', NBLINE
      WRITE (*,'(A,F5.2,A,F5.2)') 
     >      ' Clumping parameters: DENSCON(1,1) = ', DENSCON(1,1), 
     >      '   FILLFAC(1,1) = ', FILLFAC(1,1) 
      IF (MANIPOP_OPTIONS .NE. ' ') WRITE (*,'(2A)') 
     > ' Hot LTE component added: MANIPOP_OPTIONS = ', MANIPOP_OPTIONS


cc      IF (IFIRST .EQ. 1 .OR. LSPRO .GT. 0 .OR.
cc     >  (LPSTA .NE. 1 .OR. JFIRST .NE. 1 .OR. JLAST .EQ. 1) .OR. 
cc     >  NMOD .GT. 1) THEN
        PRINT 1, DXMAX, TAUMAX, XMAX, TAUMINBROAD, NFOBS
    1   FORMAT(' NUMERICAL PARAMETERS: DXMAX=',F4.2,', TAUMAX=',F4.1,
     >        ' , XMAX=',F8.1,', TAUMINBROAD=',F6.3,' , NFOBS=',I6)

        PRINT 111, NPHI, LPSTA, LPEND, JFIRST, JLAST
  111   FORMAT (' AZIMUTH ANGLES NPHI=',I4, ', LPHI FROM',I4,' TO',I4, 
     >  /, ' IMPACT PARAMETER INDEX JP FROM',I4,' TO',I4 )

      IF (IDX(DD_VDOP_LINE) .NE. 0)
     > WRITE(*,'(A, A)') ' VDOP was specified by following line: ', 
     >                  DD_VDOP_LINE

      IF (BDD_VDOP) THEN
        WRITE(*, '(A)') ' VDOP is depth-dependent'
        IF (BMICROTURB) THEN
           WRITE(*, '(A)') ' VDOP depends on element, since VMIC was '
     >                  // 'specified: VDOP^2 = vmic^2 + vtherm^2'
        ENDIF
      ENDIF

      WRITE (*,*) 
     >   'Criterion that defined VDOP for wavelength resolution:'  
      IF (IVDOPSTATUS == 1) THEN
         WRITE (*,*) 'VDOP from MODEL file' 
      ELSEIF (IVDOPSTATUS == 2) THEN
         WRITE (*,*) 'Minimum from all Doppler line widths, element = '
     >               // EL_MIN 
      ELSEIF (IVDOPSTATUS == 3) THEN
         WRITE (*,*) 'specified as VDOP in FORMAL_CARDS'
      ELSEIF (IVDOPSTATUS == 4) THEN
         WRITE (*,*) 'VDOPFE from FEDAT_FORMAL, NO-FECONVOL requested'
      ENDIF

      IF (VSINI .GT. .0 .AND. RCOROT .GT. 1.) THEN
        IF (BVSINI_AT_RCOROT) THEN
           WRITE (*,'(A,G14.4)') 'VSINI refers to RCOROT=', RCOROT 
        ELSE
           WRITE (*,*) 'VSINI refers to RSTAR'
        ENDIF
      ENDIF

      IFIRST=0

      WRITE (*,*)
      DO ISPZ=1, MAXSPZ 
        IF (SPZ1(ISPZ) .EQ. '') EXIT
        WRITE (*,'(1X,A,A10,A,A10,A)') 
     >   'Levels starting with >', SPZ1(ISPZ), 
     >   '<, have been set to Zero, except those starting with >', 
     >   SPZ2(ISPZ), '<'
      ENDDO
      WRITE (*,*)

C***  PRINT HEADER FOR ONE BLEND BLOCK (INDEX OF REFERENCE LINE: >LINE<)
      IF (REDIS) THEN
        WRITE (REDIMES, '(F3.1)') BWESEX
        REDIMES = 'ELECTRON REDISTR. BWESEX=' // REDIMES
      ELSE
        REDIMES = ' '
      ENDIF
      IF (BROAD) THEN
        BROAMES = 'BROAD.'
      ELSE
        BROAMES = ' '
      ENDIF
      PRINT 13, LINE, XLAM, VSINI, RCOROT, VDOP, REDIMES, BROAMES
   13 FORMAT (/, 132('-'),/,'LINE', I5, '   AT', F10.2, ' A',
     $ 3X, 'VSINI/(KM/S)=', F5.0, 
     > 3X, 'RCOROT/RSTAR=', F5.1,
     $ 3X, 'VDOP/(KM/S)=', F5.1, 3X, A, 3X, A,
     $ /, 132('-'), / )

C***  LOOP OVER ALL BLEND COMPONENTS
      DO 14 NBL = NBLINE ,1 ,-1
      IND = INDLAP(NBL)
      LOW = INDLOW(IND)
      NUP = INDNUP(IND)
      XLAMI=1.E8/(ELEVEL(NUP)-ELEVEL(LOW))
      F = C4 * XLAMI * XLAMI * 
     >    EINST(NUP,LOW) * WEIGHT(NUP) / WEIGHT(LOW)
C***  MARK SUBLINES BY AN ASTERISK
      IF (IND .LE. LASTIND) THEN
         SUBLINE = ' '
      ELSE
         SUBLINE = '*'
      ENDIF
      IF (LINPRO(NBL) .EQ. 'VOIGT   ') THEN
         WRITE (AVOIGTSTR,'(1PG8.2)') AVOIGT(NBL)
      ELSE
         AVOIGTSTR = ' '
      ENDIF

      XLAMBL = XLAMLAP(NBL)
      IF ( BAIRWAVELENGTH ) THEN
         XLAM2   = XLAMBL * XLAMBL
         XLAMBL = XLAMBL - XLAMBL*(2.735182E-4 + 131.4182
     >                 / XLAM2 + 2.76249E8 / (XLAM2*XLAM2))
      ENDIF

         PRINT 6, IND, SUBLINE, XLAMBL, XLAMI,  
     $    LEVEL(LOW), LEVEL(NUP), F, LINPRO(NBL), AVOIGTSTR
    6    FORMAT ('LINE', I5, A1, '  AT', F10.2, ' A',
     $      '  (From dE:', F10.2, ' A)', 3X,
     $      'LEVELS: ', A, ' (LOW) - ', A, ' (UP)', 3X, 'F= ', 
     >      1P,G10.3, 1X,A8,1X,A)
   14 CONTINUE


C***  PRINT WARNINGS ABOUT POSSIBLE FURTHER BLENDS **********************

C***  ESTABLISH BANDWIDTH OF ONE LINE:
      FMIN = XMAX + VMAX / VDOP * SQRT(1.-1./RMAX/RMAX) 
      FMAX = XMAX + VMAX / VDOP

      NADDBL=0
      DO 30 IND=1, LASTIND

C***  CHECK WHETHER THIS LINE IS ALREADY TAKEN INTO ACCOUNT
      DO 25 NBL=1,NBLINE
      IF ((IND .EQ. INDLAP(NBL)) .OR. (IND .EQ. LINE)) GOTO 30
   25 CONTINUE

C***  CHECK WHETHER THIS LINE IS SPLIT INTO A MULTIPLET
      DO 27 IM=1, NMULTI
      IF (IND .EQ. MULTIIND(IM)) GOTO 30
   27 CONTINUE

      LOW=INDLOW(IND)
      NUP=INDNUP(IND)
      XLAMI = 1.E8 / (ELEVEL(NUP)-ELEVEL(LOW))
C***  logarithmic scale
      DELTAX = ALOG (XLAMI / XLAM) / ALN

C***  BLEND, IF DELTAX IN THE INTERVAL [ FREMIN-FMAX) , FREMAX+FMIN ]:
      IF ((DELTAX-FREMIN+FMAX)*(DELTAX-FREMAX-FMIN) .LT. .0) THEN
          F = C4 * XLAMI * XLAMI * 
     *        EINST(NUP,LOW) * WEIGHT(NUP) / WEIGHT(LOW)
C***  ATTENTION: NO WARNING ISSUED IF F-VALUE IS VERY SMALL!
          IF (F .LE. 1.E-8) GOTO 30

          NADDBL=NADDBL+1
          IF (NADDBL .EQ. 1) PRINT 28
   28     FORMAT (/,'WARNING - FURTHER LINE BLENDS EXISTING:')

          IF (XLAMI .LT. XLAMMIN) THEN
             DIFF = XLAMI - XLAMMIN
             WRITE (OVERLAP,'(A,F11.1,A)') 
     >             '(BELOW RANGE BY', DIFF, ' A)'
          ELSE IF (XLAMI .GT. XLAMMAX) THEN 
             DIFF = XLAMI - XLAMMAX
             WRITE (OVERLAP,'(A,F11.1,A)')
     >             '(ABOVE RANGE BY', DIFF, ' A)'
          ELSE
             OVERLAP = '(INSIDE RANGE                                )'
          ENDIF
          PRINT 29, IND, XLAMI, OVERLAP, LEVEL(LOW), LEVEL(NUP), F
   29     FORMAT (2X, 'LINE', I5, ' AT', F10.2,' A',  
     $        4X, A, 
     $        4X, 'LEVELS: ', A, ' (LOW) - ', A, ' (UP)',
     $        4X, 'F= ', 1PG10.3) 
      ENDIF
   30 CONTINUE


C***  PRINT LIST OF TRANSITION WHICH ARE SPLIT IN MULTIPLETS ***************

      IF (NMULTI .GT. 0) THEN
      WRITE (*,'(/,A)') 
     >       'THE FOLLOWING TRANSITIONS WERE SPLIT INTO SUBLINES:'  
      DO IND=1, LASTIND
        DO 36 IM=1, NMULTI
          IF (IND .EQ. MULTIIND(IM)) GOTO 37
   36   CONTINUE
C       Index was not found in MULTIIND:
        CYCLE

   37   CONTINUE
        LOW=INDLOW(IND)
        NUP=INDNUP(IND)
        XLAMI = 1.E8 / (ELEVEL(NUP)-ELEVEL(LOW))

        PRINT 38, IND, XLAMI, LEVEL(LOW), LEVEL(NUP)
   38   FORMAT (2X, 'LINE', I5, ' AT', F10.2,' A',  
     $        4X, 'LEVELS: ', A, ' (LOW) - ', A, ' (UP)')
      ENDDO
      ENDIF


C**************************************
C***  Print Line profile (if requested)
C**************************************

      IF (LSPRO.GT.0) PRINT 8
    8 FORMAT (/,11X,'INDEX    ',  
     $ 'DELTA-NUE           DELTA LAMBDA          LAMBDA        ',
     $ '      RELATIVE         EQUIV.WIDTH',/,
     $          11X,' (K)     ',
     $ '(DOPPLER UNITS)     (ANGSTROEM)         (ANGSTROEM)     ',
     $ '        FLUX           (ANGSTROEM)',/)
 
C***  START VALUES:
      KHALF=0
      KMAX = 1
      PMAX=1.
      PHALF=1.
      EQWI=.0
      ASYM=.0
      ABEQWI=0.
      EMEQWI=0.
      DLAM2=(DLAM(2)-DLAM(1))/2.
      REFDIFF=0.05
      SPIKE=.FALSE.
 
C***  LOOP OVER ALL FREQUENCY POINTS -----------------------------------
      X=XOBS0+DXOBS
      IF (LSPRO.GT.0) PRINT 2, 1,X,DLAM(1),DLAM(1)+XLAM,PROFILE(1),EQWI
    2 FORMAT ( 11X,I4,F13.2,F20.2,F20.2,F20.3,F20.3)
      DO 3 K=2, NFOBS
      X=XOBS0+K*DXOBS
      EQWI=EQWI+(2.-PROFILE(K-1)-PROFILE(K))*DLAM2
      ASYM=ASYM+(PROFILE(K-1)*DLAM(K-1)+PROFILE(K)*DLAM(K))*DLAM2
      IF (PROFILE(K).LT.1.) THEN
         ABEQWI=ABEQWI+(2.-PROFILE(K-1)-PROFILE(K))*DLAM2
      ELSE
         EMEQWI=EMEQWI+(2.-PROFILE(K-1)-PROFILE(K))*DLAM2
      ENDIF
      IF (PROFILE(K).GT.PMAX) THEN
         KMAX=K
         PMAX=PROFILE(K)
      ENDIF
      IF (LSPRO.GT.0) PRINT 2,K,X,DLAM(K),DLAM(K)+XLAM,PROFILE(K),EQWI
      IF (K .EQ. NFOBS) GOTO 3
C***  SPIKE DETECTOR:
      PROFILK=PROFILE(K)
      DIFF1=PROFILK-PROFILE(K-1)
      DIFF2=PROFILK-PROFILE(K+1)
      IF ((DIFF1 .GT. REFDIFF) .AND. (DIFF2 .GT. REFDIFF)
     $     .AND. (DIFF1*DIFF2 .GT. 0.0)) THEN
             SPIKE =.TRUE.
             KSPIKE=K
             SPIKELAM = DLAM(K)+XLAM
      ENDIF
    3 CONTINUE
C***  ENDLOOP ----------------------------------------------------------

C***  OUTPUT: USED METHOD OF INTEGRATION IN SUBR. ZONEINT
      IF (IVERSION .EQ. 1) PRINT 21
   21 FORMAT (/,31X, '=====  CAUTION: INTEGRATION WAS PERFORMED IN Z ',
     $          ' (IVERSION = 1)  =====')
      IF (IVERSION .EQ. 2) PRINT 22
   22 FORMAT (/,31X, '=====  CAUTION: INTEGRATION VERSION: OLDDTAU ',
     $          ' (IVERSION = 2)  =====')

C***  LOOP OUTPUT:
      IF (SPIKE) PRINT 20, SPIKELAM, KSPIKE
   20 FORMAT (/,24X, '===== SPIKE at ', F20.2,
     >          ' (K=',I3,') IN LINE PROFILE DETECTED',
     >          ' - WARNING FOR QUICK LOOK ONLY =====')

      IF (BNOCONT) PRINT 26, ' '
   26 FORMAT (/,24X,'+++++ WARNING: FORMAL INTEGRAL OF AN OLD MODEL',A1,
     >        'CALCULATED WITHOUT INTERPOLATED XJC +++++')


      IF (LSPRO.LE.0)  PRINT 10, EQWI, EMEQWI, ABEQWI
   10 FORMAT (/,5X, 'EQUIVALENT WIDTHS: TOTAL:',F10.3,' A',
     $        10X, 'EMISSION:', F10.3, ' A',
     $        10X, 'ABSORPTION:', F10.3, ' A', /, 30X, 12('=') )
 
C***  OUTPUT OF INFORMATION ABOUT THE LINE PROFILE
      PHALF=PMAX/2.+0.5
      XPEAK=DLAM(KMAX)+XLAM
      DO 23 K=KMAX+1, NFOBS
        IF (PROFILE(K).GT.PHALF) KHALF=K
   23 CONTINUE
      IF (((KHALF .EQ. 0).OR.(KHALF .EQ. NFOBS)).AND.(.NOT. SPIKE)) THEN
         PRINT 9
    9    FORMAT(/,28X,'===== ATTENTION: ABSORPTION PROFILE! =====',/)
         ENDIF

      IF (KHALF .LT. 1 .OR. KHALF .GE. NFOBS) THEN
         TXHALF = ' UNDEF. '
         TXKM   = ' UNDEF. '
         ELSE
         DENOM = PROFILE(KHALF) - PROFILE(KHALF+1)
         IF (DENOM .NE. 0.) THEN 
            XHALF = DLAM(KHALF) - (DLAM(KHALF)-DLAM(KHALF+1)) 
     >              * (PROFILE(KHALF)-PHALF) / DENOM
            XKM = XHALF / XLAM * CLIGHT
            WRITE (TXHALF, '(F8.2)') XHALF
            WRITE (TXKM  , '(F8.2)') XKM
            ENDIF
         ENDIF

         PRINT 7, PMAX, XPEAK, TXHALF, TXKM
    7    FORMAT(/,5X, 'PEAK INTENSITY:', F7.3, ' AT ', F8.2,' A',
     $          10X, 'HALF WIDTH AT HALF MAXIMUM (RED WING):', A8,
     $          ' A  =', A8, ' KM/SEC', /, 20X, '=======')


C***  CALCULATION OF DIMENSIONLESS MOMENTS (LINE PROFILE)
C***  SOURCE: CASTOR ET AL. 1981, MNRAS 194, 547
      FMOM = CLIGHT / XLAM / VMAX
      W0   = -FMOM * EQWI
      W1   = FMOM * FMOM * ASYM
C>>>>      PRINT 17, W0, W1
   17    FORMAT(/,10X,'DIMENSIONLESS MOMENTS (CASTOR ',
     $          'ET AL. 1981):   W0 = ',F9.5,/,56X,'W1 = ',F9.5)
 
      IF (JFIRST .EQ. JLAST)
     > WRITE (*,'(A,I3,A,F10.4,A,F10.4,A)') 
     $   '          INTENSITY PROFILE AT P(',JFIRST,') = ', P(JFIRST), 
     $   '    FLUX integration weight p*dp=',  PWEIGHT, ' not applied'
 
      IF (BFEWARNING) THEN
      WRITE(*,*)' '
        WRITE(*,901)'WARNING:  ACTIVE IRON-LINES DETECTED BUT NOT 
     >           ACCOUNTED IN SPECTRA'
 901  FORMAT(A)
      ENDIF

      WRITE(*,'(///)')


      RETURN
      END
 
