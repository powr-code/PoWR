      SUBROUTINE TRAPLO (KANAL,PROFILE,DLAM,NFOBS,LINE,
     >                   FREQIN,MODHEAD,JOBNUM,
     $                   DISP,N,NBLINE,IDENT,OSMIN,
     >                   XLAM,XLAMLAP,INDLAP,INDLOW,INDNUP,LEVEL,
     >                   ELEVEL,WEIGHT,EINST,NDIM,ABSWAV,BCONT,FUNIT, 
     >                   NMOD, LPSTA, LPEND, NPHI, XUNIT, 
     >                   BCALIBRATED, BAIRWAVELENGTH, RSTAR, 
     >                   KODATIND, MAXATOM, NOM)

C**********************************************************************
C***               DIRECT TRANSFER OF PLOT DATA                     ***
C**********************************************************************
 
      DIMENSION DLAM(NFOBS),PROFILE(NFOBS),DLAMPLO(NFOBS)
      DIMENSION XLAMLAP(NBLINE),INDLAP(NBLINE)
      DIMENSION INDLOW(N),INDNUP(N),ELEVEL(N),WEIGHT(N)
      DIMENSION EINST(NDIM,NDIM)
      DIMENSION JOBNUM(NMOD)
      DIMENSION KODATIND(MAXATOM), NOM(NDIM)

      CHARACTER(*) :: FREQIN    !Name des RANGE-Bereiches
      CHARACTER(60) :: XCHAR 
      CHARACTER(65) :: YCHAR
      CHARACTER(100) :: HEADER
      CHARACTER(*) :: MODHEAD
      CHARACTER*10 LEVEL(N), FUNIT*(*)
      CHARACTER*1 BACKSLASH
      CHARACTER*8 IDSTR
      CHARACTER LAMBDASUBSCRIPT*7
      CHARACTER*(*) XUNIT
      LOGICAL AUTOMAX, IDENT, ABSWAV, BCONT, BCALIBRATED
      LOGICAL BAIRWAVELENGTH, BIDENT_ON
 
      REAL, PARAMETER :: PI = 3.14159   !Pi
      REAL, PARAMETER :: C4 = 1.499E-16 !C4 = 1 / (SIGMA-CLASSIC * 8 PI)     (IN ANGSTROEM**2)

C***  CLIGHT = SPEED OF LIGHT in Angstroem per second
      DATA CLIGHT2 / 2.99792458E18 /

C***  CALIBRATION CONSTANT TO USE FOR THE CONVERSION OF F-NUE TO MV
      DATA CONMV / 48.64 /

C***  PARSEC IN CM
      DATA PARSEC / 3.08561E18 /


      BACKSLASH = CHAR(92)

C**********************************************************************
C***  DEFINE HEADER LINE                                            ***
C**********************************************************************
 
      HEADER         = ' '
      HEADER(1:20) = FREQIN(1:20)

      IF (BCONT .AND..NOT. BCALIBRATED) HEADER(23:31) = 'Continuum'
      IF (BCONT .AND. BCALIBRATED     ) HEADER(23:31) = 'Absolut'

      HEADER (36:36) = 'M'
C***  DATE AND TIME
      WRITE (HEADER(37:),'(7A2)') MODHEAD(13:14),
     >       MODHEAD(16:17),MODHEAD(19:20),'  ',MODHEAD(25:26),
     >       MODHEAD(28:29),MODHEAD(31:32)
      HEADER (53:55) = 'J >'
      WRITE (HEADER (56:61),'(I6)') JOBNUM(1)

C***  Subscript to Lambda: Air or Vacuum
      IF (BAIRWAVELENGTH) THEN
         LAMBDASUBSCRIPT = '&Tair&M'
      ELSE
         LAMBDASUBSCRIPT = '&Tvac&M'
      ENDIF
 
C**********************************************************************
C***                SCALING OF Y-AXIS                               ***
C**********************************************************************

      IF (BCONT .OR. BCALIBRATED) THEN
                IF    (FUNIT .NE. 'MV'        
     >           .AND. FUNIT .NE. 'FLAM10PC'
     >           .AND. FUNIT .NE. 'LOGFLAM'
     >           .AND. FUNIT .NE. 'LOGFNUE'       ) THEN        
                    FUNIT = 'LOGFNUE'
                    WRITE(0,'(A)') '*** WARNING: ' //
     >                 'Flux units invalid: taking LOGFNUE instead'
                ENDIF

C***     Convert flux units if appropriate
         IF (FUNIT .EQ. 'LOGFLAM') THEN
            DO K=1, NFOBS
C***           Save for negative results
               IF (PROFILE(K) < 0.) THEN
                   PROFILE(K) = -100.
               ELSE
                   PROFILE(K) = ALOG10 (PROFILE(K) * CLIGHT2 / 
     >                     (DLAM(K)+XLAM)**2)
               ENDIF
            ENDDO
         ELSE IF (FUNIT .EQ. 'LOGFNUE') THEN
            DO K=1, NFOBS
C***           Save for negative results
               IF (PROFILE(K) < 0.) THEN
                   PROFILE(K) = -100.
               ELSE
                   PROFILE(K) =
     >               ALOG10 (PROFILE(K))
               ENDIF
            ENDDO
         ELSE IF (FUNIT .EQ. 'MV') THEN
             DISMO = -2.5 * ALOG10 (RSTAR*RSTAR /
     >                         100. / (PARSEC * PARSEC))
             DO K=1, NFOBS
C***           Save for negative results
               IF (PROFILE(K) < 0.) THEN
                   PROFILE(K) = -100.
               ELSE
                   PROFILE(K) =
     >               -2.5 * ALOG10 (PI*PROFILE(K)) + DISMO - CONMV
 
              ENDIF
            ENDDO
         ELSE IF (FUNIT .EQ. 'FLAM10PC') THEN
               FAC = PI * RSTAR*RSTAR / (PARSEC*PARSEC) / 100.
               DO K=1, NFOBS
                  PROFILE(K) = PROFILE(K) * CLIGHT2 / 
     >                         (DLAM(K)+XLAM)**2 * FAC 
               ENDDO 
          ENDIF
      ENDIF

C***  Invoke automatic scaling of Y-axis by WRplot
C***  Note: the formar YMAX option is now disabled (wrh 9-Jun-2014)
      YMIN=0.0
      YMAX=0.0

      IF (.NOT.((BCONT .OR. BCALIBRATED) .AND. FUNIT .EQ. 'MV')) 
     >   GOTO 1
C***  Manual scaling of Y axis for MV (negative direction!)
         YMIN = PROFILE(1)    
         YMAX = PROFILE(1)    
         DO K=2,NFOBS
            IF (PROFILE(K) .GT. YMIN) YMIN=PROFILE(K) 
            IF (PROFILE(K) .LT. YMAX) YMAX=PROFILE(K) 
         ENDDO
         YDIFF = YMAX - YMIN
         YMIN = YMIN - 0.05 * YDIFF 
         YMAX = YMAX + 0.05 * YDIFF 

         IF (YMAX .LT. -100.) THEN
           YMAX=-100
         ELSEIF (YMIN .GT. 100.) THEN
           YMIN=100.
         ENDIF

         YTICK = FLOAT (NINT(YDIFF*20.))/400.
         YABST = 10.*YTICK

    1 YSCALE=0.

C**********************************************************************
C***                SCALING OF X-AXIS                               ***
C**********************************************************************

C***  X-Units in micrometer (default: Angstroem)
      IF (XUNIT(:1) .EQ. 'M') THEN
         XUNITFAC = 1.E-4
      ELSE
         XUNITFAC = 1.
      ENDIF

      XLAMPLO = .0
      IF (ABSWAV) XLAMPLO = XLAM 

      DO 110 I=1,NFOBS
        DLAMPLO(I) = ( DLAM(I) + XLAMPLO ) * XUNITFAC
  110 CONTINUE

C**********************************************************************
C***  BREITE DES PLOTS IN CM                                        ***
C**********************************************************************

      BREITE= 20.
      IF (DISP .EQ. .0) THEN

C**********************************************************************
C***  DISPERSION NICHT SPEZIFIZIERT: PLOT UEBER DAS GERECHNETE      ***
C***                                 INTERVALL                      ***
C**********************************************************************

        XMIN = DLAMPLO(1)
        XMAX = DLAMPLO(NFOBS)

      ELSE

C**********************************************************************
C***  DISPERSION GEMAESS DER INPUT OPTION                           ***
C**********************************************************************

        XMAX= 0.5*BREITE*DISP + XLAMPLO
        XMIN=-0.5*BREITE*DISP + XLAMPLO

      ENDIF

      XSCALE=0.
 
C**********************************************************************
C***  BESCHRIFTUNGEN UND TICK-MARKS                                 ***
C**********************************************************************
      IF (XMAX-XMIN .LE. 1.) THEN
         XABST=.1
         XTICK=.02
      ELSE IF (XMAX-XMIN .LE. 3.) THEN
         XABST=.5
         XTICK=.1
      ELSE IF (XMAX-XMIN .LE. 10.) THEN
         XABST=1.
         XTICK=0.2
      ELSE IF (XMAX-XMIN .LE. 30.) THEN
         XABST=5.
         XTICK=1.
      ELSE IF (XMAX-XMIN .LE. 100) THEN
         XABST=10.
         XTICK=2.
      ELSE IF (XMAX-XMIN .LE. 300.) THEN
         XABST=50.
         XTICK=10.
      ELSE IF (XMAX-XMIN .LE. 1000) THEN
         XABST=100.
         XTICK=20.
      ELSE IF (XMAX-XMIN .LE. 3000.) THEN
         XABST=500.
         XTICK=100.
      ELSE IF (XMAX-XMIN .LE. 10000) THEN
         XABST=1000.
         XTICK=200.
      ELSE
         XABST= 5000.
         XTICK= 1000.
      ENDIF
 
C**********************************************************************
C     XMIN=XTICK*(IFIX(DLAM(  1  )/XTICK)-1.)                       ***
C     XMAX=XTICK*(IFIX(DLAM(NFOBS)/XTICK)+1.)                       ***
C**********************************************************************
 
C**********************************************************************
C***  WRITE MODEL HEADERS VERTICALLY
C**********************************************************************
      WRITE (KANAL,'(A)') 'PLOT: ' // HEADER
      XPOS = 1.5

C      DO IMOD=1, NMOD
        WRITE (KANAL, '(A,F3.1,2A)') 
     >    '\LUNA XMAX YMAX ', XPOS, ' 0. 0.3 -90. '
C        WRITE (KANAL, '(A)') 
C     >    '\> &E' // MODHEAD(IMOD)( :IDX(MODHEAD(IMOD)))
        WRITE (KANAL, '(A)') 
     >    '\> &E' // MODHEAD(:IDX(MODHEAD))
        XPOS = XPOS - 0.50
C      ENDDO

C**********************************************************************
C     IF (IDENT): LINE IDENTIFICATION 
C**********************************************************************

      IF (IDENT) THEN
          WRITE (KANAL, '(A,G8.2)') 
     >          '\LAB 17.5 15.3 0.2 IDENT: &Rf&N >= ',OSMIN
          WRITE (KANAL,*) '\ID_INBOX'
          WRITE (KANAL,*) '\IDSIZE=0.2'
          DO 105 NBL=NBLINE,1,-1
          IF (BAIRWAVELENGTH) THEN
             XLAMREF = XLAMLAP(NBL)
             XLAM2   = XLAMREF * XLAMREF
             XLAMREF = XLAMREF - XLAMREF*(2.735182E-4 + 131.4182
     >                 / XLAM2 + 2.76249E8 / (XLAM2*XLAM2))
          ELSE
             XLAMREF = XLAMLAP(NBL)
          ENDIF
          DXLAM = XLAMREF - XLAM + XLAMPLO
          DXLAM = DXLAM * XUNITFAC
          IND=INDLAP(NBL)
          NUP=INDNUP(IND)
          LOW=INDLOW(IND)
          XLAMI=1.E8/(ELEVEL(NUP)-ELEVEL(LOW))
          F = C4 * XLAMI * XLAMI * 
     >        EINST(NUP,LOW) * WEIGHT(NUP) / WEIGHT(LOW)
C***      IDENT shall not be commented for H, He or if f-value > OSMIN
          BIDENT_ON = (F .GT. OSMIN) .OR.( KODATIND(NOM(LOW)) .LE. 2)
          IF (BIDENT_ON) THEN
              IDSTR = ' \IDENT '
          ELSE
              IDSTR = '*\IDENT '
          ENDIF
          WRITE (KANAL,101) IDSTR, DXLAM, LEVEL(NUP), LEVEL(LOW)
  101     FORMAT (A8, G14.7,' &E',A10, ' - ', A10,'&N')

  105     CONTINUE
      ENDIF

      IF (XUNIT(:1) .EQ. 'M') THEN
       IF (ABSWAV) THEN
         XCHAR  = 
     >     BACKSLASH // 'CENTER' // BACKSLASH // '#l#'
     >     // LAMBDASUBSCRIPT // '  / #m#m'
       ELSE
         XCHAR = 
     >     BACKSLASH // 'CENTER' // BACKSLASH // '#D l#'
     >     // LAMBDASUBSCRIPT // '  / #m#m'
       ENDIF
      ELSE
       IF (ABSWAV) THEN
         XCHAR  = 
     >     BACKSLASH // 'CENTER' // BACKSLASH // '#l#' 
     >     // LAMBDASUBSCRIPT // '  / \A'
       ELSE
         XCHAR = 
     >     BACKSLASH // 'CENTER' // BACKSLASH // '#D l#'
     >     // LAMBDASUBSCRIPT // '  / \A'
       ENDIF
      ENDIF

C***  Y-AXIS descriptor
      IF (BCONT .OR. BCALIBRATED) THEN
        IF (FUNIT .EQ. 'MV  ') THEN
          YCHAR = BACKSLASH // 'CENTER' // BACKSLASH // 'M&T#l#&M [mag]'
        ELSE IF (FUNIT .EQ. 'FLAM10PC') THEN
          YCHAR = BACKSLASH // 'CENTER' // BACKSLASH // 
     >     'f&T#l#&M / (erg cm&H-2&M s&H-1&M ' //
     >     BACKSLASH // 'A&H-1&M) at 10 pc'
        ELSE IF (FUNIT .EQ. 'LOGFLAM') THEN
          YCHAR = BACKSLASH // 'CENTER' // BACKSLASH // 
     >     'log F&T#l#&M / (erg cm&H-2&M s&H-1&M ' //
     >     BACKSLASH // 'A&H-1&M)'
        ELSE IF (FUNIT .EQ. 'LOGFNUE') THEN
          YCHAR = BACKSLASH // 'CENTER' // BACKSLASH //
     >     'log F&T#n#&M / (erg cm&H-2&M s&H-1&M Hz&H-1&M)'
        ELSE
          YCHAR = 'units undefined -- internal error'
        ENDIF
      ELSE
C***  else: normalized profile
         YCHAR= BACKSLASH // 'CENTER' // BACKSLASH // 'rel. Flux'
         WRITE (KANAL, '(A)') BACKSLASH // 'LINUN XMIN 1 XMAX 1' 
      ENDIF

      CALL PLOTANFS (KANAL,HEADER,HEADER
     $ ,XCHAR,YCHAR
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,DLAMPLO,PROFILE,NFOBS,'COLOR=2')
      RETURN
      END
