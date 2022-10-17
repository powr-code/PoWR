      SUBROUTINE DECSTAR (MODHEAD, FM, RSTAR, VDOP, RMAX, TTABLE, 
     >                    NDDIM, OLDTEMP, MAXATOM, NATOM, ABXYZ, KODAT,
     >                    VPLOT, ATMASS, XDATA, MAXXDAT, OLDFGRID, 
     >                    THIN, ThinCard, GLOG, GEFFLOG, GEDD, bSaveGEFF,
     >                    XMSTAR, WRTYPE, TAUMAX, TAUACC,
     >                    bTauFix, BTWOT, TFAC, DENSCON_LINE, BLACKEDGE, 
     >                    bOLDSTART, RadiusGridParameters, XMDOT, XLOGL,
     >                    RTRANS, BTAUR, DENSCON_FIX, MASSORIGIN, 
     >                    LRTinput, ELEMENT, bOldStratification, 
     >                    bHYDROSOLVE, GEddFix, 
     >                    MLRELATION, NC, VTURB, 
     >                    bOLDRAD, RADGAMMASTART, fHYDROSTART, 
     >                    bFULLHYDROSTAT, bGAMMARADMEAN, GEFFKEY,
     >                    bGREYSTART, XLAMBLUE, MDOTINPUT, LTESTART)
C***********************************************************************
C***  DECODES INPUT CARDS, CALLED FROM WRSTART
C***********************************************************************
 
      IMPLICIT NONE

C***  COMMON/VELPAR/ TRANSFERS VELOCITY-FIELD PARAMETERS TO FUNCTION WRVEL
C***  The second line has additional parameters for the 2-beta-law 
C***     -- wrh  6-Apr-2006 17:21:57
      COMMON/VELPAR/ VFINAL,VMIN,BETA,VPAR1,VPAR2,RCON,HSCALE,
     >     BETA2, BETA2FRACTION, VPAR1_2, VPAR2_2

C***  COMMON /COMTEFF/  TRANSFERS THE EFF. TEMPERATURE TO: 
C***                    WRSTART, GREY, PRIMOD, DECFREQ, DTDR, JSTART 
      COMMON /COMTEFF/ TEFF,TMIN,TMODIFY,SPHERIC

C***  Operating system:
      COMMON / COMOS / OPSYS
      CHARACTER(8) :: OPSYS

      REAL :: FM, RSTAR, VDOP, RMAX, RTRANS, TEFF, 
     >          TMIN, TMODIFY, TFAC,
     >          VMIN, VFINAL, BETA, VPAR1, VPAR2,
     >          RCON, HSCALE, BETA2, BETA2FRACTION,
     >          VPAR1_2, VPAR2_2, fHYDROSTART,
     >          XMDOT, XMSTAR, XMSTARG, XLOGL,
     >          ABUND, ABREST, DENSUM, SUM,
     >          XHY, YHE, XC, XO, tempREAL, VTURB, XLAMBLUE, XMDOTTRANS

      INTEGER :: NDDIM, MAXATOM, NATOM, MAXXDAT,
     >           MASSORIGIN, NPAR, IDXA, IDX,
     >           KHY, KHE, KC, KO, IERR, MFORM, NC, LcalcCond
      INTEGER, INTENT(OUT) :: GEddFix, LRTinput, MDOTINPUT

      REAL, DIMENSION(MAXATOM) :: ABXYZ, ATMASS
      REAL, DIMENSION(MAXXDAT) :: XDATA
      INTEGER, DIMENSION(MAXATOM) :: KODAT
      INTEGER, DIMENSION(8) :: DTVALUES
      LOGICAL :: LTESTART, TTABLE, SPHERIC, 
     >           OLDTEMP, VPLOT, OLDFGRID,
     >           ABMASS,ABNUMB,THIN, BTWOT, RMAX_IN_RSUN, BTAUR
      LOGICAL, INTENT(OUT) :: bSaveGEFF, bOLDSTART, bGREYSTART,
     >                        bOldStratification, bTauFix, 
     >                        bHYDROSOLVE, 
     >                        bFULLHYDROSTAT, bGAMMARADMEAN, bOLDRAD
      LOGICAL bCalcGLOGfromM, bLRTcomplete
      CHARACTER(100) :: MODHEAD
      CHARACTER(8)  :: DAT, TIM, GEFFKEY
      CHARACTER(80) :: KARTE, ACTPAR, NEXTPAR, ThinCard
      CHARACTER(2)  :: WRTYPE
      CHARACTER(20) :: ACTPAR2
      CHARACTER(10) :: SYS, NODE, REL, VER, MACH
      CHARACTER(33) :: HOST
      CHARACTER*(*) :: DENSCON_LINE
      CHARACTER(10), DIMENSION(NATOM) :: ELEMENT
      CHARACTER(80), DIMENSION(3) :: RadiusGridParameters     !contains all RADIUS-GRID CARDS (for subr. RGRID)
      CHARACTER(9) :: MLRELATION

      REAL :: BLACKEDGE, DENSCON_FIX
      REAL :: GLOG, GEFFLOG, GEDD, q, QLOG, RADGAMMASTART
      REAL :: TAUMAX, TAUACC

      !Laufvariablen und co.
      INTEGER :: I, K, NA, NZ, IPAR,
     >           IFOUNDELEMENT

      !Konstanten
      REAL, PARAMETER :: RSUN = 6.96E10     !SOLAR RADIUS ( CM )
      REAL, PARAMETER :: GCONST = 6.670E-8  !GRAVITATION CONSTANT (CGS UNITS)
      REAL, PARAMETER :: TEFFSUN = 5780.    !EFFECTIVE TEMPERATURE OF THE SUN
      REAL, PARAMETER :: XMSUN = 1.989E33   !XMSUN = Solar Mass (g)

      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)


C***  DEFAULT VALUES **********
      XDATA(1)=.0
      LTESTART=.FALSE.
      DO i=1, 3, 1
        RadiusGridParameters(i) = ' '
      ENDDO
      TEFF=-1.
      RSTAR=-1.
      RMAX=-1.
      XLOGL = -99.
      RTRANS = -99.
      QLOG = -99.
      VDOP=-1.
      XMDOT= -99.
      XMDOTTRANS= -99.
      RMAX_IN_RSUN = .FALSE.
      DO NA=1,MAXATOM
        ABXYZ(NA)=0.
      ENDDO
      TTABLE=.FALSE.
      SPHERIC=.TRUE.
      OLDTEMP=.FALSE.
      OLDFGRID = .FALSE.
      TMODIFY=0.0
      TMIN=.0
      VPLOT=.FALSE.
      ABMASS=.FALSE.
      ABNUMB=.FALSE.
      THIN=.FALSE.
      TAUMAX = 0.
      TAUACC = 1.E-4
      bTauFix = .FALSE.
      BTWOT = .FALSE.
      DENSCON_FIX = 1.
      DENSCON_LINE = ' '
      BLACKEDGE = .0
      BTAUR = .FALSE.
      BETA2FRACTION = .0
      GEDD = -1.0
      XLAMBLUE = -1. 
C***  Number of core-rays
      NC = 4
C***  determines if GEDD is directly specified (=2) or implied (=1):
      GEddFix = 0          
C***  nedative default value indicates that GEFF has not been set in CARDS
      GEFFLOG = -1.0
      bSaveGEFF = .FALSE.    !Determines if GEFFLOG is fixed in the MODEL file
      bCalcGLOGfromM = .FALSE.
      bOLDSTART = .FALSE.   !true if OLDSTART CARDS option has been set

ccc   wrstart used the subroutine GREY for establishing the first tau-scale
ccc   (Rosseland mean, continuum opacities only). Andreas Sander
ccc   suggested to use the subr. TAUSCAL for calculating the Rosseland mean. 
ccc   To my opinion, this should be identical, but Andreas says results
ccc   differ. Till this is settled, the switch for
ccc   the new version is de-activated by setting the switch to TRUE
      bGREYSTART = .FALSE.
ccc      bGREYSTART = .TRUE.

      bOldStratification = .FALSE.  !true if stratification (ND, VELO, RADIUS, DENSCON) from old model should be kept (@todo: also MDOT?)
      bHYDROSOLVE = .FALSE. !true if hydro iteration will be done in STEAL
      fHYDROSTART = -99.    !starting from OLD MODEL using v infered from hydro (fraction > 0.)
      bFULLHYDROSTAT = .FALSE.
      bGAMMARADMEAN = .FALSE.
      bOLDRAD = .FALSE.
      RADGAMMASTART = -99.
      GEFFKEY = '        '
C***  Mass-Luminosity relation default: 
      MFORM = 2     ! 1 = Langer (1989); 2 = Graefener et al. (2011) 
      WRTYPE = ''   ! default: automatic choice of type for M-L relation
      LcalcCond = 0

C***  MASSORIGIN = Codenumber for source of the stellar mass (for PRIPARAM)
C***    0 = Mass-luminosity relation   
C***    1 = input   
C***    2 = log g   
C***    3 = log geff plus Eddington Gamma
      MASSORIGIN = 0
      VTURB = .0

C***  MDOTINPUT: origin of the mass-loss rate
C      MDOTINPUT = -1 : not specified --> ERROR
C      MDOTINPUT = 1 :  MDOT specified 
C      MDOTINPUT = 2 :  RTRANS specified 
C      MDOTINPUT = 3 :  MDOTTRANS specified (def. Graefener & Vink 2013) 
C      MDOTINPUT = 4 :  LOG Q specified
      MDOTINPUT = -1 


C***  CONSTRUCT MODEL HEADER  *********
      IF (OPSYS .EQ. 'DEC/UNIX') THEN
        CALL MY_DATE (DAT)
C        CALL DATE_AND_TIME(VALUES=DTVALUES)
C        WRITE(UNIT=DAT,FMT='(I4,"/",I2,"/",I2)') 
C     >          DTVALUES(1), DTVALUES(2), DTVALUES(3)
        CALL MY_CLOCK(TIM)
      ELSE IF (OPSYS .EQ. 'CRAY') THEN
        !CALL DATE (DAT)
        CALL DATE_AND_TIME(VALUES=DTVALUES)
        WRITE(UNIT=DAT,FMT='(I4,"/",I2,"/",I2)') 
     >          DTVALUES(1), DTVALUES(2), DTVALUES(3)
        !CALL DATE_AND_TIME(DATE=DAT)
        CALL CLOCK(TIM)
      ELSE
        WRITE (0,*) 'OPSYS NOT RECOGNIZED'
        WRITE (0,*) 'OPSYS=', OPSYS
        STOP 'ERROR IN SUBR. DECSTAR'
      ENDIF
      MODHEAD='MODEL START'
      MODHEAD(13:)=DAT
      MODHEAD(25:)=TIM

C***  READ, ECHO AND DECODE EACH INPUT CARDS ********
C***  AND ADD INFORMATION ABOUT HOST (SUBR. UNAME: SYSTEM CALL) ********
C!!!      CALL UNAME(SYS, NODE, REL, VER, MACH)
          sys = ' '
      HOST = OPSYS
      IF (SYS .EQ. 'craSH') THEN
         HOST(10:) = ' M92 2/64 (craSH Kiel)'
      ELSE IF (SYS .EQ. 'sn5208') THEN
         HOST(10:) = ' EL 4/1024 (craSHi Kiel)'
      ELSE IF (SYS .EQ. 'sn422') THEN 
         HOST(10:) = '/216 (crax Berlin)'
      ELSE IF (SYS .EQ. 'sn1619') THEN
         HOST(10:) = '2E/264 (cray Berlin)'
      ELSE IF (SYS .EQ. 'DEC/UNIX') THEN
         HOST      = 'DEC/UNIX at Potsdam'
      ELSE
         HOST(11:) = '(unknown)'
      ENDIF

C***  CHECK IF SYSTEM IS INSTALLED CORRECTLY
      IF (SYS .EQ. 'DEC/UNIX' .AND. SYS .NE. OPSYS) THEN
        WRITE (0,*) 'SYSTEM IS NOT INSTALLED CORRECTLY'
        WRITE (0,'(3A)') 'SYS, OPSYS=', SYS, OPSYS
      ENDIF

      PRINT 1, HOST
    1 FORMAT ('1>>> ECHO OF INPUT CARDS <<<', 53X, 'HOST: ', A33, /, 
     >        1X, 27('-'), 52X, 7('='), //)
      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
      REWIND (1)

      !-- Zeilenweises Abarbeiten der CARDS-Datei ----
   10 READ (1,11, END=100) KARTE
   11 FORMAT (A)

      PRINT 2,KARTE
    2 FORMAT (1X,A)
 

      IF (KARTE .EQ. '') GOTO 10
      CALL SARGV (KARTE, 1, ACTPAR)
      CALL SARGC (KARTE, NPAR)

      IF (KARTE(:5) .EQ. 'TABLE') THEN
C                         =====
            TTABLE=.TRUE.

      ELSE IF ((KARTE(:4) == 'THIN') .OR.
     >         (KARTE(:9) == 'HYDROSTAT')) THEN
C                              ====
            THIN=.TRUE.
            ThinCard = KARTE
            IF (NPAR > 2) THEN
              DO IPAR=3, NPAR
                CALL SARGV(ThinCard,IPAR,ACTPAR) 
                IF (ACTPAR == 'FULL') THEN
                  bFULLHYDROSTAT = .TRUE.
                ENDIF
                IF (ACTPAR == 'MEAN') THEN
                  bGAMMARADMEAN = .TRUE.
                ENDIF
              ENDDO
            ENDIF
            
      ELSE IF (KARTE(:5) .EQ. 'PLANE') THEN
C                              =====
            SPHERIC=.FALSE.

C***  Option GREY_TAUSCALE_START for compatibility with versions before
C***  23-Jan-2016
      ELSE IF (KARTE(:9) == 'GREYSTART' .OR. 
     >         KARTE(:9) == 'GREY_TAUS'     ) THEN
C                            =========
            bGREYSTART = .TRUE.
            
      ELSE IF (KARTE(:8) == 'OLDSTART') THEN
C                            ========
            bOLDSTART = .TRUE.

      ELSE IF (KARTE(:5) .EQ. 'OLD T') THEN
C                              =====
            OLDTEMP=.TRUE.
            IF (KARTE(:9) .EQ. 'OLD T TAU') BTAUR = .TRUE.

      ELSE IF (KARTE(:9) .EQ. 'OLD FGRID') THEN
C                              =========
            OLDFGRID = .TRUE.

      ELSE IF ((KARTE(:9) .EQ. 'OLD STRAT')
                               !=========
     >    .OR. (KARTE(:5) .EQ. 'OLD V')) THEN
                                !=====   
         bOldStratification = .TRUE.
            
      ELSE IF (KARTE(:7) .EQ. 'TMODIFY') THEN
C                              =======
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F20.0)',ERR=92) TMODIFY

      ELSE IF (KARTE(:4) .EQ. 'TMIN') THEN
C                              ====
         IF (KARTE(:10) /= 'TMIN-START') THEN
           WRITE (hCPR, '(A)') '*** WARNING: Deprecated TMIN card'
     >       // ' used: Please use TMIN-START card instead!'
         ENDIF
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F20.0)',ERR=92) TMIN

      ELSE IF (KARTE(:8) .EQ. 'LTESTART') THEN
C                              ========
         LTESTART=.TRUE.
         IF (NPAR >= 3) THEN
            CALL SARGV (KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'BLACKEDGE') THEN
               CALL SARGV (KARTE, 3, ACTPAR)
               READ (ACTPAR, '(F10.0)', ERR=92) BLACKEDGE
            ENDIF
         ENDIF

      ELSE IF (KARTE(:6) .EQ. 'PLOT V') THEN
C                              ======
            VPLOT=.TRUE.

      ELSE IF (KARTE(:5) .EQ. 'RGRID') THEN
C                              =====
            RadiusGridParameters(1) = KARTE

      ELSE IF (ACTPAR .EQ. 'RTRANS') THEN
C                           ======
         CALL SARGV (KARTE, 2, ACTPAR)
         IDXA = IDX(ACTPAR)
         IF (ACTPAR(IDXA-2:IDXA) .EQ. 'DEX') THEN
            READ (ACTPAR(:IDXA-3), '(F20.0)',ERR=92) RTRANS
            RTRANS = 10.**RTRANS
         ELSE
            READ (ACTPAR, '(F20.0)',ERR=92) RTRANS
         ENDIF
         MDOTINPUT = 2

      ELSE IF (KARTE(:5) .EQ. 'LOG Q') THEN
C                              =====
        CALL SARGV (KARTE, 3, ACTPAR)
        READ (ACTPAR, '(F20.0)',ERR=92) QLOG
        MDOTINPUT = 4
      
      ELSE IF (ACTPAR .EQ. 'RSTAR') THEN
C                           =====
      CALL SARGV (KARTE, 2, ACTPAR)
      READ (ACTPAR, '(F20.0)',ERR=92) RSTAR
 
      ELSE IF (KARTE(:5) .EQ. 'LOG L') THEN
C                              =====
      CALL SARGV (KARTE, 3, ACTPAR)
      READ (ACTPAR, '(F20.0)',ERR=92) XLOGL
  
      ELSE IF (KARTE(:6) .EQ. 'VELPAR') THEN
C                              ======
         CALL DECVELPAR(KARTE, VFINAL, VMIN, BETA, RMAX)
         
      ELSE IF (ACTPAR .EQ. '2BETALAW') THEN
         CALL SARGV (KARTE,3,ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=99) BETA2 
         CALL SARGV (KARTE,5,ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=99) BETA2FRACTION 
 
      ELSE IF (KARTE(:4) .EQ. 'VDOP') THEN
C                              ====
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F20.0)',ERR=92) VDOP
 
      ELSE IF (KARTE(:8) .EQ. 'HEADLINE') THEN
C                              ========
        MODHEAD(35:) = KARTE(10:)
 
      ELSE IF (ACTPAR .EQ. 'MDOT') THEN
C                           ====
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F20.0)',ERR=92) XMDOT
        MDOTINPUT = 1 

      ELSE IF (ACTPAR(1:5) == 'MDOTT' .OR. ACTPAR == 'MDTRANS') THEN
C                              =====                  =======
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F20.0)',ERR=92) XMDOTTRANS
        MDOTINPUT = 3

      ELSE IF (ACTPAR .EQ. 'TEFF') THEN
C                           ====
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F20.0)',ERR=92) TEFF

      ELSE IF (ACTPAR .EQ. 'NCORE') THEN
C                           =====
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(I3)', ERR=92) NC
        
      ELSE IF (KARTE(:5) .EQ. 'MSTAR') THEN
C                              =====
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR .NE. '?') THEN
          IF (MASSORIGIN .NE. 0) THEN
            WRITE (0, '(A)') '*** DOUBLE DEFINITION OF STELLAR MASS'
            GOTO 92
          ENDIF
          READ (ACTPAR, '(F20.0)', ERR=92) XMSTAR
          MASSORIGIN = 1
          LcalcCond = LcalcCond + 1
        ENDIF
 
      ELSE IF ((KARTE(:10) .EQ. 'LOG G_GRAV') .OR.
     >         (KARTE(:9) .EQ. 'LOG GGRAV') .OR.
     >         (KARTE(:8) .EQ. 'LOG GRAV')) THEN
C                              =====
        CALL SARGV (KARTE, 3, ACTPAR)
        IF (ACTPAR .NE. '?') THEN
          IF (MASSORIGIN .NE. 0) THEN
            WRITE (0, '(A)') '*** DOUBLE DEFINITION OF STELLAR MASS'
            GOTO 92
          ENDIF
          READ (ACTPAR, '(F20.0)', ERR=92) GLOG
          MASSORIGIN = 2
        ENDIF
 
      ELSE IF ((KARTE(:9) == 'LOG G_EFF') .OR.
C                               =========
     >         (KARTE(:8) == 'LOG GEFF') .OR.
C                               ========
     >         ((KARTE(:5) == 'LOG G') .AND. 
C                              =====
        !note: for backward compartibility log g is always interpreted as log geff
        !      To avoid confusion LOG G is considered as deprecated, instead you
        !      should always write LOG GEFF or LOG GGRAV
     >          ((KARTE(6:6) == ' ') .OR. (KARTE(6:6) == '=') .OR.
     >           (KARTE(6:6) == ',') .OR. (KARTE(6:6) == ':') ) ) 
     >        ) THEN

        CALL SARGV (KARTE, 3, ACTPAR)
        IF (ACTPAR /= '?') THEN
          READ (ACTPAR, '(F20.0)', ERR=92) GEFFLOG
          bSaveGEFF = .TRUE.
        ENDIF
        !Warning if deprecated syntax is used
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR == 'G') THEN
          WRITE (hCPR, '(A)') '*** WARNING: Deprecated syntax LOG G '
     >      // ' used: Please use LOG GEFF (or LOG GGRAV) instead!'
        ENDIF
        IF (bSaveGEFF) THEN
          !Check if interpretation keyword exists
          GEFFKEY = 'AUTO    '
          CALL SARGV (KARTE, 4, ACTPAR)
          IF (NPAR > 4 .AND. ACTPAR == 'RADFORCE') THEN
            CALL SARGV (KARTE, 5, ACTPAR)
            SELECTCASE(ACTPAR)
              CASE ('RAD', 'FULL')
                GEFFKEY = 'RAD     '
              CASE ('THOMSON', 'ELECTRON', 'THOM', 'E', 'e')
                GEFFKEY = 'THOM    '
            ENDSELECT
          ELSEIF (NPAR >= 4) THEN
            WRITE (hCPR, '(A)') '*** ERROR: Invalid or missing '
     >        // ' parameter on LOG GEFF card!'
            WRITE (hCPR, '(A)') '    Use RADFORCE=FULL or '
     >        // 'RADFORCE=ELECTRON to specify the intention'
     >        // ' of the LOG GEFF card!'
            GOTO 92          
          ENDIF
        ENDIF

      ELSE IF (ACTPAR == 'EDDINGTON-GAMMA') THEN
C                         ===============
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR /= 'AUTO') THEN
          READ (ACTPAR, '(F20.0)', ERR=92) GEDD
          bSaveGEFF = .TRUE.
          GEddFix = 2
          LcalcCond = LcalcCond + 2
        ENDIF

      ELSE IF (ACTPAR == 'RADGAMMA-START') THEN
C                         ==============
        CALL SARGV (KARTE, 2, ACTPAR)
        IF (ACTPAR /= 'OLD') THEN
          READ (ACTPAR, '(F20.0)', ERR=92) RADGAMMASTART
        ELSE
          bOLDRAD = .TRUE.
          IF (NPAR > 2) THEN
            CALL SARGV (KARTE, 3, ACTPAR)
            !Use MEAN value only? (for WRSTART)
            IF (ACTPAR == 'MEAN') bGAMMARADMEAN = .TRUE.
          ENDIF
        ENDIF  

      ELSE IF (KARTE(:7) == 'MLANGER') THEN
C                            =======
         MFORM = 1
      ELSE IF (KARTE(:6) == 'MGOETZ') THEN
C                            ======
         MFORM = 2
      ELSE IF (ACTPAR .EQ. 'TAUMAX') THEN
C                           ======
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F10.0)',ERR=92) TAUMAX
        IF (NPAR > 2) THEN
          DO i=3, NPAR 
            CALL SARGV (KARTE, i, ACTPAR)
            SELECTCASE (ACTPAR)
              CASE ('FIX')
                bTauFix = .TRUE.
                IF (NPAR >= (i+1)) THEN
                  CALL SARGV (KARTE, i+1, NEXTPAR)
                  READ (NEXTPAR, '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    TAUACC = tempREAL
                  ENDIF
                ENDIF
C              CASE ('MIN')             !not used in wrstart
C                bTauStrict = .FALSE.
              CASE ('EPS', 'ACC')
                IF (NPAR >= (i+1)) THEN
                  CALL SARGV (KARTE, i+1, NEXTPAR)
                  READ (NEXTPAR, '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    TAUACC = tempREAL
                  ENDIF                  
                ENDIF           
              CASE ('REPS', 'RELEPS', 'RELACC')
                IF (NPAR >= (i+1)) THEN
                  CALL SARGV (KARTE, i+1, NEXTPAR)
                  READ (NEXTPAR, '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    TAUACC = tempREAL * TAUMAX
                  ENDIF                  
                ENDIF           
C              CASE ('REDUCE')          !not used in wrstart
C                IF (NPAR >= (i+1)) THEN
C                  CALL SARGV (KARTE, i+1, NEXTPAR)
C                  READ (NEXTPAR, '(F10.0)', IOSTAT=IERR) tempREAL
C                  IF (IERR == 0) THEN
C                    ReduceTauCorrections = tempREAL
C                  ELSE
C                    ReduceTauCorrections = 0.5
C                  ENDIF                  
C                ENDIF
            ENDSELECT
          ENDDO
        ENDIF


      ELSE IF (ACTPAR == 'TAUFIX') THEN
C                           ======
        bTauFix = .TRUE.
        IF (NPAR > 1) THEN
          CALL SARGV (KARTE, 2, ACTPAR)
          IF (ACTPAR /= 'STRICT') THEN
            READ (ACTPAR, '(F10.0)', ERR=92) TAUACC
          ELSEIF (NPAR > 2) THEN
            CALL SARGV (KARTE, 3, ACTPAR)
            READ (ACTPAR, '(F10.0)', ERR=92) TAUACC                                     
          ENDIF
        ENDIF

      ELSEIF (ACTPAR == 'HYDRO') THEN
C                        =====
        bHYDROSOLVE = .TRUE.
        
      ELSE IF (ACTPAR .EQ. 'TWOTEMP') THEN
C                           ======
        BTWOT = .TRUE.
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F10.0)',ERR=92) TFAC

      ELSE IF (ACTPAR .EQ. 'DENSCON') THEN
C                           ======
        CALL SARGV (KARTE, 2, ACTPAR)
        READ (ACTPAR, '(F10.0)',ERR=92) DENSCON_FIX
        DENSCON_LINE = KARTE
        
      ELSE IF (ACTPAR .EQ. 'WRTYPE' .OR. ACTPAR .EQ. 'STARTYPE') THEN
        CALL SARGV (KARTE, 2, WRTYPE)
        IF (WRTYPE .NE. 'OB' .AND. WRTYPE .NE. 'WN' .AND.
     >      WRTYPE .NE. 'WC') THEN
            WRITE (0, *) '*** ERROR: Invalid choice of WRTYPE' 
            WRITE (0, *) '*** Valid types are: OB, WN, WC'
            GOTO 92
        ENDIF
C***    VTURB enters the hydrostatic equation via 
C***          P_gas = rho * (v_sound^2 + v_turb^2)
C***          and, therefore, differs from the definition of v_mic
C***          (hitherto only in FORMAL) for the line broadening:
C***          v_dop = sqrt( v_th^2 + v_mic^2) 
C***          Hence: v_mic = sqrt(2) * v_turb
      ELSE IF (ACTPAR .EQ. 'VTURB') THEN
        IF (NPAR .GT. 1) THEN
           CALL SARGV (KARTE, 2, ACTPAR2)
           READ (ACTPAR2,'(F20.0)',ERR=98) VTURB
        ELSE
           GOTO 97
        ENDIF
      ELSE IF (ACTPAR .EQ. 'VMIC') THEN
        IF (NPAR .GT. 1) THEN
           CALL SARGV (KARTE, 2, ACTPAR2)
           READ (ACTPAR2,'(F20.0)',ERR=98) VTURB
           VTURB = VTURB / SQRT(2.) 
        ELSE
           GOTO 97
        ENDIF

      ELSE IF (ACTPAR(:8) .EQ. 'BLUEMOST') THEN
C                           ========
         CALL SARGV (KARTE, 2, ACTPAR2)
         READ (ACTPAR2, '(F20.0)',ERR=92) XLAMBLUE


C*********** SPECIAL BLOCK FOR X-RAYS  **************************
      ELSE IF (ACTPAR .EQ. 'XRAY') THEN
C                           ====
        DO I=2, NPAR-1
          CALL SARGV (KARTE, I, ACTPAR)
          IF (ACTPAR .EQ. 'XFILL') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(1)
          ELSEIF (ACTPAR .EQ. 'XRAYT') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(2)
          ELSEIF (ACTPAR .EQ. 'XRMIN') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(3)
C         Differential Emission Measure for X-rays
          ELSEIF (ACTPAR .EQ. 'DIFF-EM-EXP') THEN
           CALL SARGV (KARTE, I+1, ACTPAR)
           READ (ACTPAR, '(F20.0)', ERR=92) XDATA(4)
           IF (XDATA(4) .NE. 1.5 .AND. XDATA(4) .NE. 2.5) THEN
              WRITE (0,*) '*** SORRY, only values 1.5 or 2.5 permitted'
              WRITE (0,*) '*** as exponents for diff. emmission measure'
              GOTO 92
           ENDIF
C         Or second component
          ELSEIF (ACTPAR .EQ. 'XFILL2') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(5)
          ELSEIF (ACTPAR .EQ. 'XRAYT2') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(6)
          ELSEIF (ACTPAR .EQ. 'XRMIN2') THEN
            CALL SARGV (KARTE, I+1, ACTPAR)
            READ (ACTPAR, '(F20.0)', ERR=92) XDATA(7)
          ENDIF
        ENDDO
C       Check for consistency
        IF (XDATA(4) .GT. 0. .AND. XDATA(5) .GT. 0.) THEN
           WRITE (0,*) '*** SORRY, DIFF-EM-EXP does NOT allow 2nd XFILL'
           GOTO 92         
        ENDIF

C***  BLACK EDGE IN THE RADIATION FIELD FOR NEW-START
      ELSE IF (ACTPAR .EQ. 'JSTART') THEN
         LTESTART = .FALSE.
         IF (NPAR .GE. 3) THEN 
            CALL SARGV (KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'BLACKEDGE') THEN
               CALL SARGV (KARTE, 3, ACTPAR)
               READ (ACTPAR, '(F10.0)', ERR=92) BLACKEDGE
            ENDIF
         ENDIF

      ELSE IF (ACTPAR .EQ. 'SPECIAL_OUTER_POINTS') THEN
         RadiusGridParameters(2) = KARTE

      ELSE IF (ACTPAR .EQ. 'SPECIAL_INNER_POINTS') THEN
         RadiusGridParameters(3) = KARTE

      ELSE IF (ACTPAR .EQ. 'RMAX_IN_RSUN') THEN
         RMAX_IN_RSUN = .TRUE.

C*********** END OF SPECIAL BLOCK FOR X-RAYS  **************************

      ELSE IF (ACTPAR .EQ. 'HELIUM') THEN
         WRITE (0,*) 'NOT ALLOWED TO SPECIFY THE HELIUM ABUNDANCE!'
         GOTO 92

      ELSE
C***  Check if the card refers to an element abundance
         CALL FINDCHARGE (ACTPAR, NZ)
C***     NZ > 0 means that an element of this name is known
         IF (NZ .GT. 0) THEN
            IFOUNDELEMENT = 0
            DO K=1, NATOM
               IF (ACTPAR .EQ. ELEMENT(K)) THEN
                  IFOUNDELEMENT = 1
                  CALL SARGV (KARTE, 2, ACTPAR2)
                  READ (ACTPAR2, '(G20.0)', ERR=99) ABUND
                  IF (ABUND .LT. 0.) THEN
                     WRITE (0,*) 'NEGATIVE ABUNDANCE ENCOUNTERED'
                     GOTO 92
                  ENDIF
                  ABXYZ(K)=ABUND
C*                Default: by number 
                  IF (NPAR .EQ. 2) THEN
                     IF (ABMASS) GOTO 91
                     ABNUMB=.TRUE.
                  ENDIF
C*                Check for "mass" or "number"
                  DO IPAR=3, NPAR
                     CALL SARGV (KARTE, IPAR, ACTPAR2)
                     CALL LOWERCASE (ACTPAR2)
                     IF (ACTPAR2(:4) .EQ. 'mass' .OR. 
     >                   ACTPAR2(:5) .EQ. '(mass') THEN
                        IF (ABNUMB) GOTO 91
                        ABMASS=.TRUE.
                     ELSEIF (ACTPAR2(:6) .EQ. 'number' .OR. 
     >                       ACTPAR2(:7) .EQ. '(number') THEN
                        IF (ABMASS) GOTO 91
                        ABNUMB=.TRUE.
                     ENDIF
                  ENDDO
                  EXIT
               ENDIF
            ENDDO

            IF (IFOUNDELEMENT .EQ. 0) THEN
               WRITE (0,*) 
     >             'ABUNDANCE GIVEN, BUT ELEMENT NOT FOUND IN DATOM'  
               GOTO 92
            ENDIF
         ENDIF


      ENDIF

      GOTO 10

 
  100 CONTINUE
      CLOSE (1)
C******************************************************************

C***  COMPUTATION OF THE RELATIVE ABUNDANCE (BY NUMBER OR MASS FRACTION)
C***  OF HELIUM
      ABREST=0.
      DO NA=1,MAXATOM
        IF (KODAT(NA) .GT. 0) THEN
          ABREST=ABREST+ABXYZ(KODAT(NA))
        ENDIF
      ENDDO

      IF (ABREST .GT. 1.) THEN
         WRITE(hCPR,'(A)') 'REL. ABUNDANCES add up to more than 100%'
         STOP 'ERROR detected in DECSTAR'
      ENDIF

C***  Helium abundance is ALWAYS the complement to 100% 
      ABXYZ(KODAT(2)) = 1.-ABREST

C***  CHECK OF UNNECESSARY ATOMIC DATA 
      DO NA=1,MAXATOM
        IF (KODAT(NA) .GT. 0) THEN
          IF (ABXYZ(KODAT(NA)) .EQ. 0.) THEN
            WRITE (0,*) 'ELEMENT ', ELEMENT(KODAT(NA)), 
     >                  ' has ZERO abundance'
            STOP 'ERROR detected in DECSTAR'
          ENDIF
        ENDIF
      ENDDO
                 
      
C***  CONVERTING (POSSSIBLE) MASS FRACTIONS INTO FRACTIONAL ABUNDANCES
C***  BY NUMBER
      IF (ABMASS) THEN
         DENSUM=0.0
         DO NA=1,NATOM
           DENSUM=DENSUM+ABXYZ(NA)/ATMASS(NA)
         ENDDO
         DO NA=1,NATOM
           ABXYZ(NA)=ABXYZ(NA)/ATMASS(NA)/DENSUM
         ENDDO
      ENDIF

C***  Calculate Helium and Carbon Mass fraction from the number fractions
      SUM = .0
      DO NA=1, NATOM
        SUM = SUM + ABXYZ(NA)*ATMASS(NA)
      ENDDO
C***  Hydrogen mass fraction
      KHY = KODAT(1)
      IF (KHY == 0) THEN
        XHY = .0
      ELSE
        XHY = ABXYZ(KHY)*ATMASS(KHY) / SUM
      ENDIF
C***  Helium mass fraction: 
      KHE = KODAT(2)
      IF (KHE .EQ. 0) THEN 
        YHE = .0
      ELSE
        YHE = ABXYZ(KHE)*ATMASS(KHE) / SUM
      ENDIF
C***  Carbon mass fraction; 
      KC = KODAT(6)
      IF (KC .EQ. 0) THEN
        XC = .0
      ELSE
        XC = ABXYZ(KC)*ATMASS(KC) / SUM
      ENDIF
C***  Oxygen mass fraction;
      KO = KODAT(8)
      IF (KO .EQ. 0) THEN
        XO = .0
      ELSE
        XO = ABXYZ(KO)*ATMASS(KO) / SUM
      ENDIF
      
C***  Automatic definition of type (OB, WN or WC)
C***   this type may be used for the M-L relation
      IF (WRTYPE .EQ. '') THEN
         IF (XHY .GE. 0.45) THEN 
           WRTYPE = 'OB' 
         ELSE IF (XC .GE. 0.1) THEN 
           WRTYPE = 'WC' 
         ELSE
           WRTYPE = 'WN' 
         ENDIF
      ENDIF      

C***  Start approximation for q (number of free electrons per mass unit)
      !@todo: Set Q depending on type of star or s.th. like that
      !q = n_e / (AMU n_i) aber diese Groessen sind erst spaeter bekannt
      !using the following approximations:
      IF (XHY >= 0.5) THEN
        q = 0.88     !H-rich Of/WN Star (fully ionized solar comp.)
      ELSEIF (XHY >= 0.1) THEN
        q = 0.39    !WNL
      ELSEIF (XC >= 0.2) THEN
        IF (XO > 0.1) THEN
          q = 0.37  !WO
        ELSE
          q = 0.21  !WC 
        ENDIF
      ELSE
        q = 0.25    !WNE
      ENDIF


C***  CHECK OF MISSING SPECIFICATIONS

      bLRTcomplete = .FALSE.
      IF (TEFF < .0) THEN
        IF ((RSTAR > 0.) .AND. (XLOGL /= -99)) THEN
          !TEFF has not set in the CARDS file but can be inferred from L and R
          !Note: RSTAR should be still in RSUN at this point
          TEFF = ( 10**(XLOGL) / RSTAR**2 )**(0.25) * TEFFSUN
          WRITE (hCPR, *) 'NOTE: Temperature was infered from',
     >               ' radius and given luminosity: ', TEFF
          bLRTcomplete = .TRUE.
          LRTinput = 3
        ELSEIF ((RSTAR > 0.) .AND. (LcalcCond == 3)) THEN
          !In the rare case that MSTAR and EDDINGTON-GAMMA are given, L can be calculated:
          XLOGL = LOG10( GEDD * XMSTAR / ( 10**(-4.51) * q ) )
          TEFF = ( 10**(XLOGL) / RSTAR**2 )**(0.25) * TEFFSUN
          WRITE (hCPR, *) 'NOTE: Temperature was infered from',
     >              ' radius and calculated luminosity: ', TEFF
          bLRTcomplete = .TRUE.
          LRTinput = 5
        ELSE
          !Not enough information to calculate TEFF
          WRITE(hCPR,'(A)') 'EFFECTIVE TEMPERATURE NOT SPECIFIED'
          STOP 'ERROR'
        ENDIF
      ENDIF
      
      IF (RadiusGridParameters(1) .EQ. ' ') THEN
            WRITE(hCPR,'(A)') 'RADIUS GRID NOT SPECIFIED'
            STOP 'ERROR'
      ENDIF

      IF (.NOT. bLRTcomplete) THEN
        !This is only called if TEFF has been set in the CARDS file
        IF (RSTAR < .0) THEN
          !Rstar has not been specified in the CARDS
          IF (XLOGL .NE. -99.) THEN
            RSTAR = 10.**(0.5*XLOGL) / (TEFF/TEFFSUN)**2
            LRTinput = 2
          ELSEIF (LcalcCond == 3) THEN
            !In the rare case that MSTAR and EDDINGTON-GAMMA are given, L can be calculated:
            XLOGL = LOG10( GEDD * XMSTAR / ( 10**(-4.51) * q ) )
            RSTAR = 10.**(0.5*XLOGL) / (TEFF/TEFFSUN)**2
            LRTinput = 4
          ELSE
              WRITE(hCPR,'(A)') 
     >          'NEITHER RSTAR NOR LOG L OR EDDINGTON-GAMMA SPECIFIED'
              STOP 'ERROR'
          ENDIF
        ELSE 
          !Rstar has been preset in the CARDS
          IF (XLOGL .NE. -99.) THEN
              !TEFF, Rstar und L gesetzt => System ueberbestimmt
              WRITE(hCPR,'(A)') 'RSTAR AND LOG L SPECIFIED BOTH'
              WRITE(hCPR,'(A)') 'BUT ONLY ONE OF THEM IS INDEPENDENT!'
              STOP 'ERROR'
          ELSE
              XLOGL = 2 * ALOG10(RSTAR) + 4 * ALOG10(TEFF/TEFFSUN)
              LRTinput = 1
          ENDIF
        ENDIF
      ENDIF

      IF (MDOTINPUT < 0) THEN
            WRITE(hCPR,'(A)') 
     >        'NEITHER MDOT NOR RTRANS NOR LOG Q SPECIFIED'
            STOP 'ERROR'
      ENDIF

C***  Calculate XMDOT from transformed mass-loss rate
      IF (MDOTINPUT == 3 .AND. XMDOT == -99.) THEN
         XMDOT = XMDOTTRANS + ALOG10(VFINAL / 1000.)
     >                      - 0.5 * ALOG10(DENSCON_FIX)
     >                      - 0.75 * ( XLOGL - 6. )
      ENDIF
      IF (RTRANS /= -99.  .AND. QLOG /= -99.) THEN
            WRITE(hCPR,'(A)') 'BOTH, LOG Q AND RTRANS, SPECIFIED'
            STOP 'ERROR'
      ENDIF
      IF (XMDOT /= -99.  .AND. QLOG /= -99.) THEN
            WRITE(hCPR,'(A)') 'BOTH, LOG Q AND MDOT, SPECIFIED'
            STOP 'ERROR'
      ENDIF
      IF (RTRANS /= -99.  .AND. XMDOT /= -99.) THEN
            WRITE(hCPR,'(A)') 'BOTH, MDOT AND RTRANS, SPECIFIED'
            STOP 'ERROR'
      ENDIF
      IF ( QLOG /= -99. .AND. RTRANS == -99.) THEN
         !Rt aus log Q berechnen
         RTRANS = (10**(-4.) / 2500.)**(2./3.) *
     >      (VFINAL * (10**(QLOG))**(2.) )**(-1./3.)
      ENDIF
      IF ( XMDOT == -99.) THEN
         !Mdot aus Rt-Def. berechnen
         XMDOT = ALOG10(VFINAL / 2500.) - 4.0 + (3./2.) * 
     >      ALOG10(RSTAR/RTRANS) - 0.5 * ALOG10(DENSCON_FIX)
      ENDIF

      IF (RMAX < .0) THEN
         WRITE(hCPR,'(A)') 'RMAX NOT SPECIFIED'
         STOP 'ERROR'
      ENDIF
      
      IF (RMAX_IN_RSUN) THEN
          RMAX = RMAX / RSTAR
          IF (RMAX < 1.) THEN
            WRITE(hCPR,'(A)') '*** FATAL ERROR: RMAX < RSTAR ***'
            WRITE(hCPR,'(A,G12.5)') ' RSTAR/SUN = ', RSTAR/RSUN
            WRITE(hCPR,'(A,G12.5)') ' RMAX/SUN  = ', RMAX/RSUN
            STOP 'ERROR IN DECSTAR'
          ENDIF
      ENDIF

      IF (RMAX < 1.) THEN
         WRITE(hCPR,'(A)') 'RMAX LESS THAN RSTAR'
         STOP 'ERROR'
      ENDIF
      
      IF (VDOP < .0) THEN
         WRITE(hCPR,'(A)') 'VDOP NOT SPECIFIED'
         STOP 'ERROR'
      ENDIF

      IF (TAUMAX <= .0) THEN
         WRITE(hCPR,'(A)') '**** WARNING: TAUMAX NOT SPECIFIED ****'
         WRITE(hCPR,'(A)')
     >     'Beware that PoWR might crash during JSTART->SPLINPO ',
     >     'if vmin is not high enough.'
      ENDIF
 
C***  OPTION 'OLDTEMP' OVERWRITES OTHER OPTIONS:
      IF (OLDTEMP) THEN
         TTABLE=.FALSE.
         SPHERIC=.FALSE.
         TMIN=.0
      ENDIF

C***  OPTION 'TTABLE' OVERWRITES OTHER OPTIONS:
      IF (TTABLE) THEN
         SPHERIC=.FALSE.
         TMIN=.0
      ENDIF
      

 
C***  COMPUTATION OF THE MASS FLUX
      IF (FM .NE. .0) FM= 10.**FM
      IF (XMDOT .NE. .0) FM= 10.**(XMDOT+3.02) /RSTAR/RSTAR

     
      !Warning if full a_rad integration should be performed, but no good start is given
      IF (bFULLHYDROSTAT .AND. (.NOT. bOldStratification)
     >      .AND. (.NOT. bOLDRAD) .AND. RADGAMMASTART < 0.) THEN
        WRITE (hCPR,'(A)') 'WARNING: HYDROSTATIC INTEGRATION'
     >                          // ' will have a poor start'
      ENDIF
      
      !System g/M, geff, GEDD ueberbestimmt        
      IF ((GEDD >= 0.) .AND. (GEFFLOG > 0.) .AND. 
     >    ((MASSORIGIN == 1) .OR. (MASSORIGIN == 2))) THEN
        WRITE(hCPR,'(A)') 'LOG G, GEDD AND LOG GEFF SPECIFIED'
        STOP 'ERROR detected in DECSTAR'
      ENDIF

      

      IF ((GEFFLOG > 0.) .AND. (MASSORIGIN == 0)) THEN
         !Spezialfall: Nur LOG GEFF (und ggf. EDDINGTON GAMMA) gegeben, aber LOG G oder MSTAR nicht
         MASSORIGIN = 3
         IF (GEDD >= 0.) THEN
            GLOG = GEFFLOG - ALOG10 ( 1. - GEDD )          
         ELSEIF (RADGAMMASTART >= 0.) THEN
            GLOG = GEFFLOG - ALOG10 ( 1. - RADGAMMASTART )          
         ELSE
C***        (CAUTION: This should not be used if FULL HYDROSTATIC INTEGRATION is performed)
            XMSTAR = 10.**( GEFFLOG - 4.4371 ) * RSTAR**2. 
     >                 + 10.**(-4.51) * q * (10.**XLOGL)
            bCalcGLOGfromM = .TRUE.
         ENDIF
      ELSEIF ( (GEFFLOG > 0.) .AND.
     >         ((MASSORIGIN == 1) .OR. (MASSORIGIN == 2)) ) THEN
        !System stark bestimmt, Gamma implizit fest (=> keine Berechnung ueber q im weiteren Code)
        GEddFix = 1
        WRITE(hCPR,*)
     >    "WARNING: Eddington Gamma has been implicitly fixed"
        WRITE(hCPR,*)
     >    ' (Remove either g_eff, g_grav or Mstar from CARDS',
     >    ' file if you want to avoid this.)'
      ENDIF

C***  STELLAR RADIUS IN CM
      RSTAR=RSTAR*RSUN


C***  LOG G or MSTAR might not be specified; in that case, MSTAR is 
C***     calculated from the mass-luminosity-relation
      IF (MASSORIGIN == 0) THEN
C***     M-L relations: 
C***        OB type: Goetz (mandatory) - H-burner
C***        WN type: Goetz (default) or Langer '89 (optional) - He burner
C***        WC type: Langer (mandatory)
         IF (WRTYPE == 'WC' .OR. 
     >      (WRTYPE == 'WN' .AND. MFORM == 1) .OR.
     >       WRTYPE == 'WC') THEN
            CALL MLANGER (XMSTARG, TEFF, RSTAR, YHE, WRTYPE)         
            XMSTAR = XMSTARG / XMSUN
            MLRELATION='Langer'
         ELSE
            CALL MGOETZ (XMSTAR, TEFF, RSTAR, XHY, WRTYPE)
            XMSTARG = XMSTAR * XMSUN
            MLRELATION='Graefener'
         ENDIF
      ENDIF

      IF ((MASSORIGIN >= 2) .AND. (.NOT. bCalcGLOGfromM)) THEN
         !LOG G vorgegeben oder aus LOG GEFF und GEDD berechnet
         XMSTARG = 10.**GLOG * RSTAR * RSTAR / GCONST 
         XMSTAR = XMSTARG / XMSUN
      ELSE
         !Masse direkt vorgegeben
         bCalcGLOGfromM = .TRUE.
      ENDIF

      IF (bCalcGLOGfromM) THEN         
         XMSTARG = XMSTAR * XMSUN          
         GLOG = ALOG10(GCONST * XMSTARG / RSTAR / RSTAR)
      ENDIF

C***  Calculate GEFFLOG or GEDD or both (if not specified from input)
      IF (GEDD >= 0.) THEN    
C***    EDDINGTON-GAMMA was specified:
        IF (GEDD > 1.) THEN
          WRITE(hCPR,*) "DECSTAR: EDDINGTON-GAMMA may not exceed 1"
          STOP 'FATAL ERROR in DECSTAR'
        ENDIF
        GEFFLOG = ALOG10( (10**GLOG) * (1. - GEDD) ) 
        GEDDFIX = 2         !Gamma explizit fest
      ELSEIF (GEFFLOG >= 0.) THEN
C***    EDDINGTON-GAMMA implicitely specified (from LOG GEFF and LOG GRAV)     
        GEDD = 1. - 10**( GEFFLOG - GLOG )
      ELSEIF (RADGAMMASTART >= 0.) THEN
        GEFFLOG = ALOG10( (10.**GLOG) * (1. - RADGAMMASTART) ) 
      ELSE
C***    EDDINGTON-GAMMA must be estimated (using Thomson only)
C***    CAREFUL: This branch is also called if RADGAMMASTART: OLD is used
C***    Real GEFFLOG is calculated in PREP_GAMMARAD instead
        GEDD = 10**(-4.51) * q * (10.**XLOGL) / XMSTAR
        IF (GEDD >= 1.) THEN
          WRITE(hCPR,*) 'GEDD:', GEDD
          WRITE(hCPR,*) "DECSTAR: Initial GEDD guessing failed"
          WRITE(hCPR,*) "-- Star is beyond the Eddington limit!"
          WRITE(hCPR,'(A,F7.2)') '  Suggestion: Set MSTAR > ', 
     >                   10**(-4.51) * q * (10**XLOGL)
          STOP 'ERROR in DECSTAR'
        ENDIF
        GEFFLOG = ALOG10( (10**GLOG) * (1. - GEDD) ) 
      ENDIF      

      IF (GEddFix > 0) THEN
        WRITE(hCPR,*) "Fixed Eddington Gamma value: ", GEDD
        IF (GEDD < 0.) THEN
          WRITE(hCPR,*) "DECSTAR: unphysical negative Eddington Gamma"
          STOP 'FATAL ERROR in DECSTAR'
        ENDIF
      ENDIF
  
C***  Consistency check between GEFF keyword and HYDROSTATIC INTEGRATION
      IF (THIN .AND. TRIM(GEFFKEY) /= '') THEN
        IF (TRIM(GEFFKEY) == 'AUTO') THEN
          WRITE(hCPR,*) "DECSTAR: Deprecated Syntax for LOG GEFF"
          WRITE(hCPR,*) 'When using the HYDROSTATIC INTEGRATION card, ' 
     >            // 'the intention of LOG GEFF should be specified.'
          WRITE(hCPR,*) 'Add RADFORCE=ELECTRON or RADFORCE=FULL as a '
     >            // 'keyword after LOG GEFF!'                  
          STOP 'ERROR in DECSTAR'
        ELSEIF (bFULLHYDROSTAT .AND. TRIM(GEFFKEY) /= 'RAD') THEN
          WRITE(hCPR,*) "DECSTAR: Inconsistent meaning for LOG GEFF"
          WRITE(hCPR,*) 'When using the FULL option on the HYDROSTATIC' 
     >      // ' INTEGRATION card, '
          WRITE(hCPR,*) 'LOG GEFF must also be flagged with '
     >      // ' RADFORCE=FULL.'
          STOP 'ERROR in DECSTAR'
        ELSEIF (.NOT. bFULLHYDROSTAT .AND. TRIM(GEFFKEY) /= 'THOM') THEN
          WRITE(hCPR,*) "DECSTAR: Inconsistent meaning for LOG GEFF"
          WRITE(hCPR,*) 'When using the HYDROSTATIC INTEGRATION card' 
     >      // ' without FULL option, '
          WRITE(hCPR,*) 'LOG GEFF must be flagged with '
     >      // ' RADFORCE=ELECTRON.'
          STOP 'ERROR in DECSTAR'
        ENDIF
      ENDIF
  
  
      IF (bOldStratification) THEN
        THIN = .FALSE.      !Ignore HYDROSTATIC INTEGRATION in WRSTART if OLD STRAT is used
      ENDIF
  
      RETURN


C***  ERROR EXITS ****************************************************

   91 CONTINUE
      WRITE (0,*) 'ERROR: ' //  
     >  'ALL abundances must be given EITHER by mass OR by number!'
      WRITE (0,*)   
     >   'Default if Element line has no further option: by number' 
      GOTO 92

   92 CONTINUE
      WRITE (0,*)'DECSTAR: ERROR WHILE DECODING THE FOLLOWING LINE:'
      WRITE (0,*) KARTE
      STOP 'ERROR'

   97 WRITE (0,'(A)') '*** ERROR: PARAMETER MISSING'
      WRITE (0,'(A)') 'THE ERROR OCCURED IN THE FOLLOWING LINE: '
      WRITE (0,'(A)') KARTE(:IDX(KARTE))
      STOP 'ERROR IN DECSTAR'

   98 WRITE (0,'(A)') 
     >  '*** ERROR: THE FOLLOWING STRING COULD NOT BE DECODED AS A '
     > // 'FLOATING POINT NUMBER:', ACTPAR2(:IDX(ACTPAR2)), 
     >    'THE ERROR OCCURED IN THE FOLLOWING LINE:', KARTE
      STOP 'ERROR IN DECSTAR'

   99 WRITE (0,'(A)') 
     >  '*** ERROR: THE FOLLOWING STRING COULD NOT BE DECODED AS A '
     > // 'FLOATING POINT NUMBER:', ACTPAR(:IDX(ACTPAR)), 
     >    'THE ERROR OCCURED IN THE FOLLOWING LINE:', KARTE
      STOP 'ERROR IN DECSTAR'

      END
