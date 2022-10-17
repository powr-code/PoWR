      SUBROUTINE DECFORM (KARTE, LSOPA, RANGERED, RANGEBLUE,
     $                   VSINI,DISP, REDIS, BWESEX, 
     $                   LSPRO, LSDWL, FIN, VDOP, IDTERM, 
     $                   JPFIRST_ORIG,JPLAST_ORIG,IVERSION,
     >                   TAUMAX, XMAX, DXMAX, PATH_VCSSB,
     $                   IDENT,OSMIN,BROAD, MODHEAD,JOBNUM,
     >                   STRING1,MAXSTRI,NSTRING,ABSWAV,BCONT,FUNIT, 
     >                   LINELIST,
     >                 BIRONLINES, BNOCONT, SPZ1, SPZ2, MANIPOP_OPTIONS, 
     >                 XUNIT, BCALIBRATED, FREQIN, MACROCLUMPLINE, 
     >                 PATH_LEMKE_DAT,BAIRWAVELENGTHSET,
     >                 BAIRWAVELENGTH, MAXSPZ, TRANSDWLLINE, RCOROTLINE,
     >                 DD_VDOP_LINE, BPLOTVDOP,
     >                 BMICROTURB, LIMB_LINE,
     >                 VDOPPLOT_LINE, bBIGBANDLIMIT, 
     >                 NOWIND_LINE, TAUMINBROAD, 
     >                 bDDOPAFORMCMF, bDDFECONVOL, 
     >                 LPHISTA_ORIG, LPHIEND_ORIG,
     >                 BVSINI_AT_RCOROT, DX_WINDROT, SECONDMODEL_LINE)
C***********************************************************************
C***  DECODING INPUT CARDS FOR PROGRAM "FORMAL" FROM FTO5 = $IN (DEFAULT)
C***  THIS ROUTINE IS LEFT WHEN A "LINE" OR "BLEND" OPTION IS FOUND,
C***  AND ENTERED AGAIN WHEN THE LINE HAS BEEN COMPUTED.
C***********************************************************************
 
      LOGICAL FIN, REDIS, IDENT, BMICROTURB
      LOGICAL BVSINI_AT_RCOROT
      LOGICAL BROAD, BAIRWAVELENGTHSET,BAIRWAVELENGTH
      LOGICAL ABSWAV, BCONT, LINELIST, BIRONLINES, BNOCONT, BCALIBRATED
      LOGICAL :: BTRANSVELO, BTRANSOD,  BPLOTVDOP, 
     >           bBIGBANDLIMIT, bDDOPAFORMCMF, bDDFECONVOL
      CHARACTER*132 STRING1(0:MAXSTRI), FUNIT*(*)
      CHARACTER KARTE*(*), ACTPAR*80, MODHEAD*100, MACROCLUMPLINE*(*)
      CHARACTER ACTPAR2*80
      CHARACTER TRANSDWLLINE*(*), RCOROTLINE*(*), DD_VDOP_LINE*(*)
      CHARACTER LIMB_LINE*(*), VDOPPLOT_LINE*(*), NOWIND_LINE*(*)
      CHARACTER*20 FREQIN 
      CHARACTER*256 PATH_VCSSB, PATH_LEMKE_DAT
      CHARACTER*(*) SPZ1(MAXSPZ), SPZ2(MAXSPZ), MANIPOP_OPTIONS, XUNIT
      CHARACTER*(*) SECONDMODEL_LINE

      REAL :: TAUMINBROAD
      
C***  Set Defaults
      DO I=1, MAXSPZ
         SPZ1(I) = ''
         SPZ2(I) = ''
      ENDDO
      NSPZ = 0
      RANGERED  = .0
      RANGEBLUE = .0

      JPFIRST_ORIG = 0
      JPLAST_ORIG  = 999

      LPHISTA_ORIG = 0
      LPHIEND_ORIG = 999

      SECONDMODEL_LINE = ''

C***  BEFORE ANY NEW DECODED 'LINE'-CARD OR 'BLEND'-BLOCK:
C***  FIRST, CLEAR THE STRING ARRAY FOR COMMENTS WHICH ARE TO BE 
C***  WRITTEN TO FILE 'PLOT' (KANAL = 1)
      NSTRING = 0
      DO 11 NSTRI=0,MAXSTRI
   11 STRING1(NSTRI) = ' '

    1 READ (2,2, IOSTAT=IOSTAT) KARTE
    2 FORMAT (A)

C*** END OF INPUT DATA REACHED
      IF (IOSTAT .LT. 0 .OR. KARTE(:3) .EQ. 'END') THEN
            FIN=.TRUE.
            RETURN
            ENDIF

      CALL SARGC(KARTE,NPAR)
      IF (NPAR. LE. 0) GOTO 1
      IF (KARTE(1:1) .EQ. '*') GOTO 1
      CALL SARGV (KARTE, 1, ACTPAR)

C*** NOTE: 'LINE'-CARDS OR 'BLEND'-BLOCKS ARE FURTHER DECODED IN SUBR. PREFORM
      IF ( KARTE(:4) .EQ. 'LINE' .OR. KARTE(:5) .EQ. 'BLEND') RETURN
C                          ====                       =====

      IF ( KARTE(:5) .EQ. 'RANGE') THEN
C                          =====
         IF (NPAR .LT. 3) THEN
            WRITE (0,*) '*** ERROR: RANGE option needs two parameters'
            GOTO 99
         ENDIF 
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) RANGEBLUE
         CALL SARGV (KARTE, 3, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) RANGERED
         IF (RANGEBLUE .EQ. RANGERED) THEN
            WRITE (0,*) '*** ERROR: RANGE has zero length'
            GOTO 99
         ELSEIF (RANGEBLUE .GT. RANGERED) THEN
            SCRATCH   = RANGEBLUE
            RANGEBLUE = RANGERED
            RANGERED  = SCRATCH
         ENDIF

      ELSE IF ( ACTPAR .EQ. 'PLOT' ) THEN
C                               ====
            CALL SARGV(KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'VDOP') THEN
               BPLOTVDOP = .TRUE.
               VDOPPLOT_LINE = KARTE
            ENDIF

      ELSE IF (ACTPAR(:4) .EQ. 'LIMB') THEN
C                               ====
            CALL SARGV(KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'OFF') THEN
               LIMB_LINE = "OFF"
            ELSE
               LIMB_LINE = KARTE
            ENDIF

      ELSE IF (KARTE(:9) .EQ. 'TRANS DWL') THEN
C                              =========
            IF (LSDWL.LT.0) LSDWL=2
            IF (LSDWL.EQ.1) LSDWL=4            
C** TRANSDWLLINE holds TRANS DWL x-axis options. See routine "tradwl" for description.           
            TRANSDWLLINE = KARTE
            
      ELSE IF ( KARTE(:9) .EQ. 'PRINT DWL') THEN
C                               ========
            IF (LSDWL.LT.0) LSDWL=1
            IF (LSDWL.EQ.2) LSDWL=4

      ELSE IF ( ACTPAR .EQ. 'VSINI' ) THEN
C                            =====
            CALL SARGV (KARTE,2,ACTPAR2)
            READ (ACTPAR2, '(F20.0)', ERR=98) VSINI
            DX_WINDROT = 1. ! Default
            DO IPAR=3, NPAR
              CALL SARGV (KARTE,IPAR,ACTPAR2)
              IF (ACTPAR2     .EQ. 'RSTAR') BVSINI_AT_RCOROT = .FALSE. 
              IF (ACTPAR2(:4) .EQ. 'RCOR' ) BVSINI_AT_RCOROT = .TRUE. 
              IF (ACTPAR2(:2) .EQ. 'DX') THEN
                 CALL SARGV (KARTE,IPAR+1,ACTPAR2)
                 READ (ACTPAR2, '(F20.0)', ERR=98) DX_WINDROT
                 IF (DX_WINDROT .LE. .0) THEN
                    WRITE (0,*) '*** ERROR: DX must be positive!'
                    GOTO 99
                 ENDIF
              ENDIF
            ENDDO
              
C** RCOROTLINE contains RCOROT and unit description. See routine "rotation_prep" for description.            
      ELSE IF ( ACTPAR .EQ. 'RCOROT' ) THEN
C                            ======
         RCOROTLINE = KARTE
      ELSE IF (KARTE(:4) .EQ. 'DISP') THEN
C                              ====
         CALL SARGV(KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) DISP

      ELSE IF (ACTPAR .EQ. 'SECONDMODEL') THEN
C                           ===========
         SECONDMODEL_LINE(1+IDX(SECONDMODEL_LINE):)
     >         = ' ' // KARTE(12:IDX(KARTE))

      ELSE IF ( ACTPAR .EQ. 'LPHISTA' ) THEN
C                            =======      
         CALL SARGV(KARTE, 2, ACTPAR)
         READ (ACTPAR, '(I10)', ERR=98) LPHISTA_ORIG

      ELSE IF ( ACTPAR .EQ. 'LPHIEND' ) THEN
C                            =======      
      
         CALL SARGV(KARTE, 2, ACTPAR)
         READ (ACTPAR, '(I10)', ERR=98) LPHIEND_ORIG
         
      ELSE IF ( ACTPAR .EQ. 'JPFIRST' ) THEN
C                               ======
         CALL SARGV(KARTE, 2, ACTPAR)
         READ (ACTPAR, '(I3)', ERR=98) JPFIRST_ORIG

      ELSE IF ( ACTPAR .EQ. 'JPLAST' ) THEN
C                               ======
         CALL SARGV(KARTE, 2, ACTPAR)
         READ (ACTPAR, '(I3)', ERR=98) JPLAST_ORIG

      ELSE IF ( ACTPAR == "NOWIND" ) THEN
C                          ======      
         NOWIND_LINE = KARTE 

      ELSE IF ( KARTE(:9) .EQ. 'PRINT PRO') THEN
C                               ========
            LSPRO=1

      ELSE IF ( KARTE(:10) .EQ. 'PRINT OPAL' ) THEN
C                                ==========
            CALL SARGV(KARTE, 3, ACTPAR)
            READ (ACTPAR, '(F10.0)', ERR=98) XL
            LSOPA=1
            IF (XL .GT. 1.) LSOPA=IFIX(XL)

      ELSE IF ( ACTPAR .EQ. 'VDOP' ) THEN
C                               ====
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) VDOP 
         DD_VDOP_LINE = KARTE
      ELSE IF ( ACTPAR .EQ. 'VMIC' ) THEN
C                               ====
         BMICROTURB = .TRUE.
         DD_VDOP_LINE = KARTE

      ELSE IF ( KARTE(:8) .EQ. 'IVERSION' ) THEN
C                               ========
            CALL SARGV (KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'Z'  ) THEN 
               IVERSION=1
            ELSE IF (ACTPAR .EQ. 'TAU') THEN 
               IVERSION=0
            ELSE
               IVERSION=0 
               WRITE (0, '( A)') '************************************'
               WRITE (0, '(2A)') 'WARNING: INVALID IVERSION: ', ACTPAR
               WRITE (0, '( A)') 'TAU VERSION WILL BE APPLIED INSTEAD!'
               WRITE (0, '( A)') '************************************'
            ENDIF

      ELSE IF ( KARTE(:5) .EQ. 'REDIS' ) THEN
C                               =====
            REDIS = .TRUE.

      ELSE IF ( KARTE(:7) .EQ. 'NOREDIS' ) THEN
C                               =======
            REDIS = .FALSE.

      ELSE IF ( KARTE(:5) .EQ. 'CALIB' ) THEN
C                               =====
            BCALIBRATED = .TRUE.
            
            CALL SARGC (KARTE, NPAR)
            IF (NPAR .GT. 1) CALL SARGV (KARTE, 2, FUNIT)
                
      ELSE IF ( KARTE(:7) .EQ. 'NOCALIB' ) THEN
C                               =======
            BCALIBRATED = .FALSE. 

      ELSE IF ( KARTE(:4) .EQ. 'CONT' ) THEN
C                               =====
            BCONT = .TRUE.
            
            CALL SARGC (KARTE, NPAR)
            IF (NPAR .GT. 1) CALL SARGV (KARTE, 2, FUNIT)

      ELSE IF ( KARTE(:6) .EQ. 'NOCONT' ) THEN
C                               =====
            BCONT = .FALSE.

      ELSE IF ( KARTE(:5) .EQ. 'IDENT' ) THEN
C                               =====
            IDENT = .TRUE.
            CALL SARGV (KARTE, 2, ACTPAR)
            IF (ACTPAR .EQ. 'OSMIN') THEN
                CALL SARGV (KARTE, 3, ACTPAR)
                READ (ACTPAR, '(F10.0)', ERR=98) OSMIN
            ENDIF

      ELSE IF ( KARTE(:7) .EQ. 'NOIDENT' ) THEN
C                               =======
            IDENT = .FALSE.

      ELSE IF ( KARTE(:6) .EQ. 'BWESEX' ) THEN
C                               ======
            CALL SARGV(KARTE, 2, ACTPAR)
            READ (ACTPAR, '(F10.0)', ERR=98) BWESEX

      ELSE IF ( KARTE(:5) .EQ. 'BROAD' .OR. 
C                               =====
     >          KARTE(:8) .EQ. 'ALLBROAD') THEN
C                               ========
            BROAD = .TRUE.

      ELSE IF ( KARTE(:7) .EQ. 'NOBROAD' .OR.
C                               =======
     >         KARTE(:10) .EQ. 'NOALLBROAD' ) THEN
C                               ==========
            BROAD = .FALSE.

      ELSE IF ( KARTE(:8) .EQ. 'ABS WAVE' ) THEN
C                               ========
            ABSWAV = .TRUE.

      ELSE IF ( KARTE(:8) .EQ. 'REL WAVE' ) THEN
C                               ========
            ABSWAV = .FALSE.

      ELSE IF ( KARTE(:8) .EQ. 'LISTONLY' ) THEN
C                               ========
            LINELIST = .TRUE.

      ELSE IF ( KARTE(:6) .EQ. 'STRING' ) THEN
C                               ======
         CALL SARGC (KARTE, NPAR)
         IF (NPAR .LE. 2) THEN 
            WRITE (*,*) 'STRING option has not enough parameters'
            WRITE (0,*) 'STRING option has not enough parameters'
            STOP '*** ERROR detected by Subr. DECFORM'
         ENDIF
         CALL SARGV (KARTE, 2, ACTPAR)
         IF (ACTPAR .NE. 'COMMENT') THEN
            WRITE (*,*) 'STRING COMMENT is the only supported option'
            WRITE (0,*) 'STRING COMMENT is the only supported option'
            STOP '*** ERROR detected by Subr. DECFORM'
         ENDIF
         CALL SARGREST (KARTE,NPAR, 3, ISTART, IEND)         
         IEND = MIN0(IEND, 126)
         IF (NSTRING .EQ. 0) THEN
            WRITE (STRING1(0),'(A)') 'PLOT: ' // KARTE(ISTART:IEND)
            FREQIN=KARTE(ISTART:IEND)
         ENDIF
         NSTRING = NSTRING+1
         IF (NSTRING .LE. MAXSTRI) THEN
               WRITE (STRING1(NSTRING),'(A)') KARTE(ISTART:IEND)
         ELSE
            WRITE (*,*) 'WARNING: too many string options - ignored:'
            WRITE (*,*) KARTE
            WRITE (0,*) 'WARNING: too many string options - ignored:'
            WRITE (0,*) KARTE
         ENDIF

      ELSE IF ( ACTPAR .EQ. 'TAUMAX' ) THEN
C                            ======
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) TAUMAX 

      ELSE IF ( ACTPAR .EQ. 'XMAX' ) THEN
C                            ======
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) XMAX 

      ELSE IF ( ACTPAR .EQ. 'DXMAX' ) THEN
C                            ======
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) DXMAX 

      ELSE IF ( ACTPAR == 'TAUBROAD' ) THEN
C                          ========
         CALL SARGV (KARTE, 2, ACTPAR)
         READ (ACTPAR, '(F10.0)', ERR=98) TAUMINBROAD

      ELSE IF ( ACTPAR .EQ. 'NO-IRONLINES' ) THEN
C                            ===========
              BIRONLINES =.FALSE.

      ELSE IF ( ACTPAR .EQ. 'IRONLINES' ) THEN
C                            =========
              BIRONLINES = .TRUE.

      ELSE IF ( ACTPAR .EQ. 'WAVELENGTH' ) THEN
C                            ==========
         CALL SARGV (KARTE, 2, ACTPAR)
         IF ( ACTPAR .EQ. 'AIR') THEN
              BAIRWAVELENGTHSET=.TRUE.
              BAIRWAVELENGTH=.TRUE.
         ELSE IF ( ACTPAR .EQ. 'VACUUM' ) THEN
              BAIRWAVELENGTHSET=.TRUE.
              BAIRWAVELENGTH=.FALSE.
         ELSE
              WRITE (0,*) '*** INVALID OPTION'
              WRITE (0,*) '*** Allowed are only AIR or VACUUM'
              GOTO 99
         ENDIF

      ELSE IF ( ACTPAR .EQ. 'NO-MODELCONT' ) THEN
C                            ===========
              BNOCONT=.TRUE.

      ELSE IF ( ACTPAR .EQ. 'SET_POP_ZERO' ) THEN
C                            ============
         CALL SARGV (KARTE, 2, ACTPAR)
         NSPZ = NSPZ + 1
         IF (NSPZ .GT. MAXSPZ) THEN
            WRITE (0,'(2A, I4)') 
     >       '*** WARNING: more SET_POP_ZERO commands than ',
     >       'dimensioned - MAXSPZ =', MAXSPZ
         ELSE
            SPZ1(NSPZ) = ACTPAR
            IF (NPAR .GE. 4) THEN
              CALL SARGV (KARTE, 3, ACTPAR)
              IF (ACTPAR .EQ. 'EXCEPT') THEN
                CALL SARGV (KARTE, 4, ACTPAR)
                SPZ2(NSPZ) = ACTPAR
              ELSE
                WRITE (0,'(A,A10)') 'Keyword wrong: ACTPAR = ', ACTPAR
                WRITE (0,'(2A)') 'KARTE = ', KARTE( :IDX(KARTE))
                STOP 'ERROR in Subr. DECFORM'
              ENDIF
            ENDIF
         ENDIF
      ELSE IF ( ACTPAR .EQ. 'MANIPOP_OPTIONS' ) THEN
C                            ===============
         CALL SARGP (KARTE, NPAR, 2, ISTART, IEND)
         IF (ISTART .GT. 0) MANIPOP_OPTIONS = KARTE(ISTART:)

      ELSE IF ( ACTPAR .EQ. 'XUNIT' ) THEN
C                            =====
         CALL SARGV (KARTE, 2, XUNIT)

      ELSE IF ( ACTPAR .EQ. 'MACROCLUMP') THEN
C                            ==========
         MACROCLUMPLINE = KARTE

      ELSE IF ( ACTPAR .EQ. 'PATH_VCSSB' ) THEN
C                            ==========
         CALL SARGV (KARTE, 2, PATH_VCSSB)

      ELSE IF ( ACTPAR .EQ. 'PATH_LEMKE_DAT' ) THEN
C                            ==========
         CALL SARGV (KARTE, 2, PATH_LEMKE_DAT)

      ELSE IF ( ACTPAR == 'NO-BIGBANDCUT' ) THEN
C                          =============
         bBIGBANDLIMIT = .FALSE.
         
      ELSE IF ( ACTPAR == 'BIGBANDCUT' ) THEN
C                          ==========
         bBIGBANDLIMIT = .TRUE.
         
      ELSE IF ( ACTPAR == 'NO-DDVDOPREDIS' ) THEN
C                          ==============
C***     Switches off depth-dependent opacity profile functions in FORMCMF
C***     (use for comparison calculations with models before Feb 2017)
         bDDOPAFORMCMF = .FALSE.
         
      ELSE IF ( ACTPAR == 'DDVDOPREDIS' ) THEN
C                          ===========
         bDDOPAFORMCMF = .TRUE.
         
      ELSE IF ( ACTPAR == 'NO-FECONVOL' ) THEN
C                          ===========
C***     Switches off depth-dependent iron opacity convolution
C***     (use for comparison calculations with models before Feb 2017)
         bDDFECONVOL = .FALSE.
         
      ELSE IF ( ACTPAR == 'FECONVOL' ) THEN
C                          ========
         bDDFECONVOL = .TRUE.
         
      ENDIF
 
      GOTO 1

C***  Error branches
   98 WRITE (0,*) '*** ERROR when decoding parameter as number'

   99 WRITE (0,*) '*** The error occured in the following line:'
      WRITE (0,*) KARTE(:IDX(KARTE))
      STOP ' *** FATAL ERROR DETECTED BY SUBR. DECFORM'

 
      END
 

