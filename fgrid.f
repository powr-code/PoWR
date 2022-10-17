      SUBROUTINE FGRID (NFDIM, NF, XLAMBDA, FWEIGHT, KEY, NOM, SYMBOL, 
     $                  N, NCHARG, ELEVEL, EION, EINST, NDIM,
     $                  EDGEK, KODAT, MAXATOM,
     $                  INDNUP, INDLOW, LASTIND, KONTNUP, KONTLOW, 
     >                  LASTKON, OLDFGRID, NF2, XLAMBDA2, VDOP, XLAMBLUE)
 
C***********************************************************************
C***  GENERATION OF THE FREQUENCY GRID AND INTEGRATION WEIGHTS
C***  INCLUDING PREDEFINED FREQUENCY POINTS (FROM TAPE6)
C***  AND THE CONTINUUM FREQUENCY POINTS (NO LINE FREQUENCY POINTS)
C***  XLAMBDA = CORRESPONDING WAVELENGTH POINTS IN ANGSTROEMS
C***********************************************************************
 
      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE

      REAL, DIMENSION(NDIM,NDIM) :: EINST
      REAL, DIMENSION(NFDIM) :: XLAMBDA, XLAMBDA2, FWEIGHT
      REAL, DIMENSION(N) :: ELEVEL, EION
      REAL, DIMENSION(MAXATOM,MAXATOM) :: EDGEK
      INTEGER, DIMENSION(N) :: NCHARG, NOM
      INTEGER, DIMENSION(LASTIND) :: INDNUP, INDLOW
      INTEGER, DIMENSION(LASTKON) :: KONTNUP, KONTLOW
      INTEGER, DIMENSION(MAXATOM) :: KODAT
      LOGICAL :: OLDFGRID, BADD
      CHARACTER(2), DIMENSION(MAXATOM) :: SYMBOL
      CHARACTER(8) :: NAME, NEDGE, CKEY
      CHARACTER(8), DIMENSION(NFDIM) :: KEY
      CHARACTER(100) :: CFORMAT, MODOLD
 
C***  Constants
      REAL, PARAMETER :: STEBOL = 1.8046E-5     !STEFAN-BOLTZMANN CONSTANT / PI  (ERG/CM**2/S/STERAD/KELVIN**4)
      REAL, PARAMETER :: CLIGHT = 2.99792458E18 !SPEED OF LIGHT IN ANGSTROEM / SECOND

C***  File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)
      
C***  BEGIN OF THE LISTING
      WRITE (hOUT,FMT='(A,//,20X,A,/)') '1',
     >        '2. FREQUENCY POINTS AND INTEGRATION WEIGHTS'
 
C***  TAKE THE FREQUENCY GRID (HAND-MADE WAVELENGTH POINTS) 
C***      EITHER FROM TAPE6 = FGRID, OR FROM OLD MODEL FILE
      CALL       DECFREQ (XLAMBDA, NF, NFDIM, TREF, OLDFGRID, KEY, 
     >                    MODOLD, XLAMBLUE)

C***  The FREQUENCY GRID IS TESTED FOR TEFF AND FOR 2.0 * TEFF;
C***  THE SECOND VALUE IS THE EXPECTED TEMPERATURE AT TAUROSS = 20
      TREF2 = 2.0 * TREF
 
C***  SET KEYWORD ARRAY TO BLANKS
      DO K=1,NF
        KEY(K) = '        '
      ENDDO
 
C***  ADDITION OF THE REDMOST LINE FREQUENCY POINT (MINUS EPSILON!)
 
C***  FIND THE REDMOST LINE CENTER FREQUENCY POINT
      IF (XLAMBDA(NF) .GT. XLAMBDA(1)) THEN
          REDMOST=XLAMBDA(NF)
      ELSE
          REDMOST=XLAMBDA(1)
      ENDIF
      DO IND=1,LASTIND
        J=INDNUP(IND)
        I=INDLOW(IND)
        WLENG=1.E8/(ELEVEL(J)-ELEVEL(I))
        IF (WLENG .GT. REDMOST) THEN
            REDMOST=WLENG
            INDRED=IND
            IRED=I
            JRED=J
        ENDIF
      ENDDO
      IND=INDRED
      WLENG = REDMOST * (1. + VFINAL / 100000.)
C***  Minimum redmost wavelength is 100000. A
      IF (WLENG .LT. 100000.) THEN
        WLENG = 100000.
        WRITE (UNIT=NAME, FMT='(A8)') 'REDMIN  '
      ELSE
        IF (EINST(IRED,JRED) .EQ. -2.) THEN
            WRITE (UNIT=NAME, FMT='(A4,I4)') 'RUD ', IND
        ELSE
            WRITE (UNIT=NAME, FMT='(A4,I4)') 'LINE', IND
        ENDIF
      ENDIF

      IF (NF .GT. 1) THEN 
         CALL SEQUIN(NFDIM,NF,XLAMBDA,K,WLENG)
         IF (K.GT.0) CALL SEQUINE (NFDIM,NF-1,KEY,K,NAME)
         IF (K.LT.0) KEY(-K)=NAME
      ELSE
C***  If FGRID gives only one lambda point, 
C***    the redmost lambda becomes the second point  
         NF = 2
         XLAMBDA(2) = WLENG
         KEY(2) = NAME
      ENDIF 
 
C***  INSERTION OF CONTINUUM EDGES
C***  One point each is inserted before and after the edge
C***  at minus and plus 0.3*VDOP  
      EDGEWIDTH = 0.3 * VDOP * 1.E13 / CLIGHT   

      DO KON=1,LASTKON
        J=KONTNUP(KON)
        I=KONTLOW(KON)
        WLENG=1.E8/(ELEVEL(J)+EION(I)-ELEVEL(I))
        WPLUS = (1. + EDGEWIDTH) * WLENG
        CALL SEQUIN(NFDIM,NF,XLAMBDA,K,WPLUS)
        IF (SYMBOL(NOM(I)) /= 'G ') THEN
          WRITE (UNIT=NEDGE, FMT='(A5,A2,I1)') 
     >       'EDGE+', SYMBOL(NOM(I)),NCHARG(I)+1
        ELSE 
          WRITE (UNIT=NEDGE, FMT='(A5,A1,I2)') 
     >       'EDGE+', SYMBOL(NOM(I))(1:1),NCHARG(I)+1
        ENDIF
        IF (K.GT.0) CALL SEQUINE (NFDIM,NF-1,KEY,K,NEDGE)
        WMINUS = (1. - EDGEWIDTH) * WLENG
        CALL SEQUIN(NFDIM,NF,XLAMBDA,K,WMINUS)
        IF (SYMBOL(NOM(I)) /= 'G ') THEN
          WRITE (UNIT=NEDGE, FMT='(A5,A2,I1)') 
     >      'EDGE-', SYMBOL(NOM(I)),NCHARG(I)+1
        ELSE 
          WRITE (UNIT=NEDGE, FMT='(A5,A1,I2)') 
     >      'EDGE-', SYMBOL(NOM(I)),NCHARG(I)+1
        ENDIF
        IF (K.GT.0) CALL SEQUINE (NFDIM,NF-1,KEY,K,NEDGE)
      ENDDO
 
C***  INSERTION OF K-SHELL-IONISATION EDGES

      DO NZ=1,MAXATOM
         KOJ=KODAT(NZ)
         IF  (KOJ == 0) CYCLE
         DO ISTATE=1, NZ-2
            IF (EDGEK(KOJ,ISTATE) .EQ. .0) CYCLE
            WLENG = 1.E8 / EDGEK(KOJ,ISTATE)
            WPLUS = (1. + EDGEWIDTH) * WLENG
            CALL SEQUIN (NFDIM, NF, XLAMBDA, K, WPLUS)
            WRITE (UNIT=NEDGE, FMT='(A5,A2,I1)') 
     >          'K-ED+', SYMBOL(KOJ), ISTATE
            IF (K .GT. 0) CALL SEQUINE (NFDIM,NF-1,KEY,K,NEDGE)
            WMINUS = (1. - EDGEWIDTH) * WLENG
            CALL SEQUIN (NFDIM,NF,XLAMBDA,K,WMINUS)
            WRITE (UNIT=NEDGE, FMT='(A5,A2,I1)') 
     >          'K-ED-', SYMBOL(KOJ), ISTATE
            IF (K .GT. 0) CALL SEQUINE (NFDIM,NF-1,KEY,K,NEDGE)
         ENDDO
      ENDDO

      BTOT = STEBOL * TREF**4.
      BTOT2 = STEBOL * TREF2**4.

      NF2 = NF - 2
      DO K=1, NF2
        XLAMBDA2(K) = XLAMBDA(K+1)
      ENDDO

C***  Entry-Point for additional frequencies
      IADD = 0
  100 CONTINUE

C***  Calculate RELMAX, and the Index of its Maximum (KMAX) for both
C***  Temperatures
      KMAX1 = 0
      KMAX2 = 0
      RELMAX1 = 0.
      RELMAX2 = 0.
      DO K=1,NF-1
        XLAM1 = XLAMBDA(K)
        XLAM2 = XLAMBDA(K+1)
        W = (1./XLAMBDA(K) - 1./XLAMBDA(K+1)) * CLIGHT
        RELCONT1 = 0.5 * W * 
     >             (BNUE(XLAM1,TREF) + BNUE(XLAM2,TREF)) / BTOT
        IF (RELCONT1 .GT. RELMAX1) THEN
          KMAX1 = K
          RELMAX1 = RELCONT1
        ENDIF
        RELCONT2 = 0.5 * W * 
     >             (BNUE(XLAM1,TREF2) + BNUE(XLAM2,TREF2)) / BTOT2
        IF (RELCONT2 .GT. RELMAX2) THEN
          KMAX2 = K
          RELMAX2 = RELCONT2
        ENDIF
      ENDDO

C***  Insert additional frequency points
C***  First: to ensure small relative contributions
      BADD = .FALSE.
      XNF = FLOAT(NF)
      IF (RELMAX1 .GT. RELMAX2) THEN
        IMAX = KMAX1
        RELMAX = RELMAX1
      ELSE
        IMAX = KMAX2
        RELMAX = RELMAX2
      ENDIF
      IF (NF .GT. NFDIM-1) THEN
        WRITE (0,*) 
     >        'WARNING : NFDIM insufficient for automatic Continuum', 
     >        '  Frequency-Point adjustment'
        WRITE (0,'(A, I5)') 'NFDIM=', NFDIM
        WRITE (*,*) 
     >        'WARNING : NFDIM insufficient for automatic Continuum', 
     >        '  Frequency-Point adjustment'
        WRITE (*,'(A, I5)') 'NFDIM=', NFDIM
      ELSE
C***  Criterium set to 2.5 percent
        IF (RELMAX .GT. 0.025) THEN
C***    Red: Put New Point in the Intervall IMAX and IMAX+1
          IF (IMAX .LT. NF) THEN
            WNEWR = (XLAMBDA(IMAX) + XLAMBDA(IMAX+1)) / 2.
            WRITE (UNIT=CKEY, FMT='(A8)') KEY(IMAX+1)
            CALL SEQUIN(NFDIM,NF,XLAMBDA,K,WNEWR)
            IADD = IADD + 1
            NAME = 'ADD     '
            IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)
            BADD = .TRUE.
          ENDIF
        ENDIF
      ENDIF
      IF (BADD) GOTO 100

C***  Second: To ensure small Delta-Lambda / Lambda steps
      IADD_DL = 0
C***  Define Treshhold to distiguish between intervalls
C***    Treshhold is the first Continuum frequency point redwards starting value
ccc      TRESH_START  = 2000.
      TRESH_START2 = 227.83774
      TRESH_START3 = 504.259
      THRESH_IRMID = 30000.  ! 3 micron
      CRIT1 = 0.1
      CRIT2 = 1.0
      IF (TRESH_START2 .LT. XLAMBDA(1) .OR. 
     >    TRESH_START2 .GT. XLAMBDA(NF-1)) THEN
        WRITE (0,'(2A,3(G12.5,1X))') 
     >              'Troubles in determining TRESHHOLD: ', 
     >              'TRESH_START2, XLAMBDA(1), XLAMBDA(NF-1)=', 
     >               TRESH_START2, XLAMBDA(1), XLAMBDA(NF-1)
      ENDIF
      IF (TRESH_START3 .LT. XLAMBDA(1) .OR. 
     >    TRESH_START3 .GT. XLAMBDA(NF-1)) THEN
        WRITE (0,'(2A,3(G12.5,1X))') 
     >              'Troubles in determining TRESHHOLD: ', 
     >              'TRESH_START3, XLAMBDA(1), XLAMBDA(NF-1)=', 
     >               TRESH_START3, XLAMBDA(1), XLAMBDA(NF-1)
      ENDIF

  110 CONTINUE
      BADD = .FALSE.
C***  Find first Step larger than criterion
      DO K=1, NF-1
        IF (XLAMBDA(K) >= THRESH_IRMID) THEN
          CRIT = CRIT2
        ELSEIF (XLAMBDA(K) .GE. TRESH_START3) THEN
          CRIT = 0.1 * CRIT1
        ELSE
          IF (XLAMBDA(K) .LE. TRESH_START2 .AND. 
     >        (XLAMBDA(K)-20.) .GT. 0.) THEN
            F2 = ( (XLAMBDA(K)-20.)**8 / 6.4E16 ) + 1
C!!!        write (0,*) '!!!', K, XLAMBDA(K), F2
          ELSE
            F2 = 1.
          ENDIF
          IF (XLAMBDA(K) .LE. TRESH_START3 .AND. 
     >        (XLAMBDA(K)-300.) .GT. 0.) THEN
            F3 = ( (XLAMBDA(K)-300.)**8 / 6.4E16 ) + 10
C!!!        write (0,*) '!!!', K, XLAMBDA(K), F3
          ELSE
            F3 = 1.
          ENDIF
          CRIT = CRIT1 / (F2 * F3)
        ENDIF
        XL_MID = 0.5 * (XLAMBDA(K) + XLAMBDA(K+1))
        DLL = (XLAMBDA(K+1) - XLAMBDA(K)) / XL_MID
        IF (DLL .GT. CRIT) THEN
          KINS = K
          BADD = .TRUE.
          EXIT
        ENDIF
      ENDDO

C***  Insert new Point in Interval KINS, KINS+1
      IF (BADD) THEN
        IF (NF .GT. NFDIM-1) THEN
          WRITE (0,*) 
     >          'WARNING : NFDIM insufficient for automatic Continuum', 
     >          '  Frequency-Point adjustment'
          WRITE (0,'(A, I5)') 'NFDIM=', NFDIM
          WRITE (*,*) 
     >          'WARNING : NFDIM insufficient for automatic Continuum', 
     >          '  Frequency-Point adjustment'
          WRITE (*,'(A, I5)') 'NFDIM=', NFDIM
        ELSE
          IF (IMAX .LT. NF) THEN
            WRITE (UNIT=CKEY, FMT='(A8)') KEY(IMAX+1)
            CALL SEQUIN(NFDIM,NF,XLAMBDA,K,XL_MID)
            IADD_DL = IADD_DL + 1
            NAME = 'ADD_DL  '
            IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)
          ENDIF
        ENDIF
        GOTO 110
      ENDIF

C***  New 17-Aug-2015: 
C***  insert exact wavelengths of Smith ubv monochromatic magnitudes
C***  at 3650., 4270., 5160. Ang
      XLAMSMITH = 3650.
      CALL SEQUIN(NFDIM,NF,XLAMBDA,K,XLAMSMITH)
      NAME = 'SMITH u '
      IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)

      XLAMSMITH = 4270.
      CALL SEQUIN(NFDIM,NF,XLAMBDA,K,XLAMSMITH)
      NAME = 'SMITH b '
      IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)

      XLAMSMITH = 5160.
      CALL SEQUIN(NFDIM,NF,XLAMBDA,K,XLAMSMITH)
      NAME = 'SMITH v '
      IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)

      XLAMMASSEY = 6000.
      CALL SEQUIN(NFDIM,NF,XLAMBDA,K,XLAMMASSEY)
      NAME = 'MASSEY r'
      IF (K .GT. 0) CALL SEQUINE(NFDIM, NF-1, KEY, K, NAME)

C***  All additional Frequency Points are now inserted

C***  Now OMIT frequency points which are too close!
C***  Criterion is 1.0 * VDOP
      KLAST=2
  200 DO K=KLAST, NF
        IF (XLAMBDA(K)/XLAMBDA(K-1) .LT. (1.+EDGEWIDTH) ) THEN
          WRITE (0,'(A,F10.2,2X,A)') 
     >    'FGRID: Cont. frequency point omitted: ', XLAMBDA(K), KEY(K) 
          NF = NF - 1
          KLAST=K
          DO KK=K, NF
            XLAMBDA(KK) = XLAMBDA(KK+1)
            KEY    (KK) = KEY    (KK+1)
          ENDDO
          GOTO 200
        ENDIF
      ENDDO

C***  FREQUENCY INTEGRATION WEIGHTS ACCORDING TO THE TRAPEZOIDAL RULE **********
      NFM=NF-1
      DO K=2,NFM
        FWEIGHT(K)=.5*(1./XLAMBDA(K-1) - 1./XLAMBDA(K+1))*CLIGHT
      ENDDO
      FWEIGHT(1)=.5*(1./XLAMBDA(1) - 1./XLAMBDA(2))*CLIGHT
      FWEIGHT(NF)=.5*(1./XLAMBDA(NFM) - 1./XLAMBDA(NF))*CLIGHT
 
C***  RENORMALIZATION TO RETAIN THE EXACT INTEGRAL SUM OF PLANCKS FUNCTION
C***   AT REFERENCE TEMPERATURE TREF
      SUM=.0
      DO K=1,NF
        SUM=SUM+FWEIGHT(K)*BNUE(XLAMBDA(K),TREF)
      ENDDO
      RENORM=BTOT/SUM
      DO K=1,NF
        FWEIGHT(K)=FWEIGHT(K)*RENORM
      ENDDO
 
C***  OUTPUT
C***  CONTINUATION OF THE LISTING
      IF (OLDFGRID) THEN
        WRITE (hOUT, FMT='(A,A)') ' FGRID TAKEN FROM OLD MODEL: ',MODOLD
      ENDIF

      WRITE (hOUT,*) 
      WRITE (hOUT,'(A)') 
     >      'NOTE : Continuum frequencies inserted'
      WRITE (hOUT,'(A,I4)')
     >      '       Relative Contribution Criterion: ', IADD
      WRITE (hOUT,'(A,I4)')
     >      '       Delta-Lambda / Lambda Criterion: ', IADD_DL

      PRINT 18
   18 FORMAT(//,1X,
     > '  NR   LAMBDA/ANGSTROEM    KEY     WAVENUMBER(KAYSER)  ',
     > ' FREQUENCY(HERTZ)        WEIGHT    REL.INTEGRAL CONTRIBUTION',/)

      CFORMAT = '(1X,I5,F15.2,2X,A8,1X,F15.2,2E22.6,F8.3,3X,F8.3)'
      DO K=1,NF
        XLAM=XLAMBDA(K)
        RELCO=NF*FWEIGHT(K)*BNUE(XLAM,TREF)/BTOT
        RELCO2=NF*FWEIGHT(K)*BNUE(XLAM,TREF2)/BTOT2
        WRITE (hOUT, FMT=CFORMAT)  K, XLAM, KEY(K), 1.E8/XLAM,
     >             CLIGHT/XLAM, FWEIGHT(K), RELCO, RELCO2       
      ENDDO

      WRITE (hOUT, FMT='(//,A,F10.6,A,F8.0,A)')
     >    ' RENORMALIZATION FACTOR :', RENORM,
     >    '     REFERENCE TEMPERATURE :', TREF, ' K'

      IF (NF > 9999) THEN
        WRITE (0,*)
     >     '*** MORE THAN 9999 COARSE FREQUENCY POINTS ENCOUNTERED ***'
        WRITE (0,*) 'This is not compatible with the encoding of the'
        WRITE (0,'(A)') ' frequency index in the MODEL file variables'
     >     // ' XJCnnnn, WJCnnnn and EDDInnnn.'
        STOP 'FATAL ERROR IN FGRID'
      ENDIF
      
      RETURN
      END
