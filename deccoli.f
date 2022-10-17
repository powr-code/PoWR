      SUBROUTINE DECCOLI (LSOPA,LSINT,LINE,NLINE,MAXIND, MODHIST,
     >             LASTIND,LPOPAB,LPOPABD,LEVDEBUG,LPJNUE,
     >             LPJNUED,LASERV,PARALAS, LPSNUE,LPSNUED,
     $             MAXPLOT, RANGE1, RANGE2, EXLAM1, EXLAM2,MAXEXT,
     $             LAINDHI, DRNORUD, BLLIST, 
     >             NEWWRC, BCOLIRAY, 
     >             CLHLP, BITCONT, BPLOT,  
     >             IPLOT, LPLOT, ND, OPC, 
     >             IVERS, BEMIX, EMIXSTART, 
     >             BEMIXFIX, EMIXFIX, IVERS_FE_EXPFAC, BPLOTALPHA,
     >             XLAM_FINE_START, XLAM_FINE_END, LPLOT_WCHARM, 
     >             XLP1, XLP2, GAMMACOLI, GAMMAT, UNLU_TAUMAX, 
     >             UNLU_TAUMAX2, bKALPHA, bHYDROSOLVE, VELOMODFAK,
     >             bForceCOLIP, bTDIFFUS)
C*******************************************************************************
C***  DECODES INPUT CARDS FOR MAIN PROGRAM "COLI"
C*******************************************************************************

      !IMPLICIT NONE
 
      DIMENSION LINE(1),MODHIST(1)
      CHARACTER(80) :: KARTE
      LOGICAL :: DRNORUD, BLLIST, 
     >           BCOLIRAY, CLHLP, BITCONT, 
     >           BPLOT, BEMIX, BEMIXFIX, BPLOTALPHA
      INTEGER, DIMENSION(MAXPLOT) :: LPOPAB, LPOPABD, LPJNUE, LPJNUED,
     >                               LPSNUE, LPSNUED
      REAL, INTENT(INOUT) :: VELOMODFAK
      REAL :: PARALAS, tempREAL
      CHARACTER(30), DIMENSION(20) :: ACTPAR
      CHARACTER(8) :: OPC
      LOGICAL, INTENT(INOUT) :: bKALPHA, bForceCOLIP, bTDIFFUS
      LOGICAL, INTENT(OUT) :: bHYDROSOLVE
       
      INTEGER :: LSOPA, LSINT, NLINE, IPLOT, LPLOT
      REAL :: RANGE1, RANGE2, DUNLU, DUNLUR


      !File and channel handles (=KANAL)
      INTEGER, PARAMETER :: hOUT = 6        !write to wruniqX.out (stdout)
      INTEGER, PARAMETER :: hCPR = 0        !write to wruniqX.cpr (stderr)


C***  DEFAULT VALUES
      LSOPA   = -1
      LSINT   = -1
      NLINE   = 0
      LEVDEBUG= 0
      LASERV  = 1
      PARALAS = 0.01
      DRNORUD = .FALSE.
      BLLIST  = .FALSE.
      BCOLIRAY = .FALSE.
      CLHLP = .FALSE.
      BITCONT = .TRUE.
      BPLOT = .FALSE.
      IPLOT = ND / 2
      LPLOT = 10000
      OPC = 'NONE'
      IVERS = 4
      BPLOTALPHA = .FALSE.
      bHYDROSOLVE = .FALSE.
      bTDIFFUS = .TRUE.
    
C***  Version of exponential factor in Fe lines (rates and eta) 
      IVERS_FE_EXPFAC = 1

C***  NUMBER OF REPEAT ITERATIONS BETWEEN TWO COLI+ JOBS
C***   This default must agree with the definition in DECSTE (STEAL)
      NEWWRC=6

C***  Default for the LINE card: ALL lines (must agree with decste!)
      IND1 = 1
      IND2 = LASTIND

C***  Default for the DRLINES card: no lines 
      INDDR1 = 1
      INDDR2 = 0

C***  EDDIMIX enabled
      BEMIX = .TRUE.
      EMIXSTART = 1.0
      BEMIXFIX = .FALSE.
      EMIXFIX = 1.0

      DO I=1,MAXPLOT
        LPOPAB (I) = 0
        LPOPABD(I) = 0
        LPJNUE (I) = 0
        LPJNUED(I) = 0
        LPSNUE (I) = 0
        LPSNUED(I) = 0
      ENDDO

      NOPLOT = 0
      NJPLOT = 0
      NSPLOT = 0

      RANGE1 = -1.
      RANGE2 = -1.

      EXLAM1 = -1.
      EXLAM2 = -1.

C***  Default values for extension of the fine grid (for Integration in STEALCL)
      XLAM_FINE_START =  100.
      XLAM_FINE_END   = 2000.

C***  Default-Range which is plotted
      XLP1 = 100.
      XLP2 = 1000.

C***  Special PLOT of WCHARM factors; If this set to zero, no plot is prepared
      LPLOT_WCHARM = 0

c      GAMMACOLI = 0.  !this will ALWAYS lead to a crash
      GAMMACOLI = 1.

C***  Parameters for Unsoeld-Lucy method, needed by FREQUINT
      GAMMAT = 0.
      UNLU_TAUMAX = 1000.
      UNLU_TAUMAX2 = 100.

      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
      REWIND 1
    6 READ (1,FMT='(A)', END=100) KARTE
 
      CALL SARGC(KARTE,NPAR)
      IF ( NPAR < 1) GOTO 6
      IF ( NPAR > 20) NPAR = 20

      DO I=1, NPAR
        CALL SARGV(KARTE,I,ACTPAR(I))
      ENDDO

      IF (ACTPAR(1) .EQ. 'PRINT') THEN
C                         =====
        IF ( ACTPAR(2) .EQ. 'INTL') THEN
C                            ====
            READ (ACTPAR(3),8, ERR=99) XL
    8       FORMAT (F10.0)
            LSINT=IFIX(XL)
            IF (LSINT.EQ.0) LSINT=1

        ELSE IF ( ACTPAR(2) .EQ. 'OPAL') THEN
C                                 ====
            READ (ACTPAR(3),8, ERR=99) XL
            LSOPA=IFIX(XL)
            IF (LSOPA.EQ.0) LSOPA=1

        ELSE IF ( ACTPAR(2) .EQ. 'LINELIST') THEN
C                                 ========
            BLLIST = .TRUE.
        ENDIF
        GOTO 6
      ENDIF

      IF ( ACTPAR(1)(:4) .EQ. 'LINE' ) THEN
C                              ====
         IF (ACTPAR(2) .EQ. 'NONE') THEN
            IND1=1
            IND2=0
         ELSEIF (ACTPAR(2) .NE. 'ALL') THEN
            READ (ACTPAR(2),1, ERR=99) IND1
    1       FORMAT (I4)
            IF ((NPAR .GT. 3) .AND. (ACTPAR(3) .EQ. 'TO')) THEN
               READ (ACTPAR(4),1, ERR=99) IND2
            ELSE
               IND2=IND1
            ENDIF
         ENDIF
         GOTO 6
      ENDIF

      IF ( ACTPAR(1)(:6) .EQ. 'DRLINE' ) THEN
C                              ======
         DRNORUD = .TRUE.
         IF (ACTPAR(2) .EQ. 'ALL') THEN
            INDDR1 = LASTIND + 1
            INDDR2 = LAINDHI
         ELSE
            READ (ACTPAR(2),1, ERR=99) INDDR1
            IF ((NPAR .GT. 3) .AND. (ACTPAR(3) .EQ. 'TO')) THEN
               READ (ACTPAR(4),1, ERR=99) INDDR2
            ELSE
               INDDR2=INDDR1
            ENDIF
         ENDIF
         GOTO 6
      ENDIF
 
      IF (ACTPAR(1) .EQ. 'CMFDLEV') THEN
C                         =======
        READ(ACTPAR(2),'(I10)', ERR=99) LEVDEBUG
        GOTO 6
        ENDIF


      IF (ACTPAR(1) .EQ. 'LASERV') THEN
C                         ======
        READ (ACTPAR(2),'(I10)', ERR=99) LASERV
        IF (NPAR == 3) THEN
          READ (ACTPAR(3),'(G10.1)', ERR=99) PARALAS
        ELSEIF (NPAR > 3) THEN
          DO IPAR=3, NPAR-1
            SELECTCASE (ACTPAR(IPAR))
              CASE ('LINES', 'LINE', 'L')
                IF (NPAR >= (IPAR+1)) THEN
                  READ (ACTPAR(IPAR+1), '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    PARALAS = tempREAL
                  ENDIF                    
                ENDIF
            ENDSELECT
          ENDDO
        ENDIF
        GOTO 6
        ENDIF

      IF (ACTPAR(1) .EQ. 'EXTEND') THEN
C                         ======
        READ (ACTPAR(2), '(F10.0)', ERR=99) EXLAM1
        READ (ACTPAR(3), '(F10.0)', ERR=99) EXLAM2
        GOTO 6
      ENDIF

      IF (ACTPAR(1) .EQ. 'NEWWRC' ) THEN
C                         ======
         READ (ACTPAR(2),'(F10.0)',ERR=99) XL
         NEWWRC=IFIX(XL)
      ENDIF

      IF (ACTPAR(1) .EQ. 'OB-VERS') THEN
C                         =======
         IF (NPAR .EQ. 2 .OR. ACTPAR(3) .EQ. 'COLI') THEN
           READ (ACTPAR(2),'(I10.0)', ERR=99) IVERS
         ENDIF
      ENDIF

      IF (ACTPAR(1) == 'NOTDIFFUS') THEN
C                       =========
         bTDIFFUS = .FALSE.
      ENDIF
         
C********************  CMF PLOT OPTION  **********************************

      IF (ACTPAR(1) .EQ. 'CMFPLOT') THEN
C                         =======
        IF (NPAR .LT. 2) THEN
         PRINT *,'CMFPLOT: MISSING OPTION(S)'
         GOTO 6
         ENDIF
C***  CMFPLOT RANGE lambda1 lambda2
C***    - WITH THIS OPTION, ONLY ONE PLOT OF EACH QUANTITY (J, KAPPA, S) 
C***     IS POSSIBLE
        IF (ACTPAR(2) .EQ. 'RANGE') THEN
C                           =====
           MAXPLOT = 1
           READ (ACTPAR(3), '(F10.0)', ERR=99) RANGE1
           READ (ACTPAR(4), '(F10.0)', ERR=99) RANGE2
           GOTO 6
           ENDIF

        IF (ACTPAR(2) .EQ. 'OPAB') THEN
C                           ====
         IF (NOPLOT.EQ.MAXPLOT) THEN
          PRINT *,'CMF DECODE WARNING > EXTRA PLOTOPTIONS IGNORED'
          GOTO 6
         ENDIF
         NOPLOT=NOPLOT+1
         READ(ACTPAR(3),'(I10)', ERR=99) LPOPAB(NOPLOT)
         READ(ACTPAR(4),'(I10)', ERR=99) LPOPABD(NOPLOT)
         GOTO 6
         ENDIF

        IF (ACTPAR(2) .EQ. 'JNUE') THEN
C                           ====
         IF (NJPLOT.EQ.MAXPLOT) THEN
          PRINT *,'CMF DECODE WARNING > EXTRA PLOTOPTIONS IGNORED'
          GOTO 6
         ENDIF
         NJPLOT = NJPLOT + 1
         READ (ACTPAR(3),'(I10)', ERR=99) LPJNUE(NJPLOT)
         READ (ACTPAR(4),'(I10)', ERR=99) LPJNUED(NJPLOT)
         GOTO 6
         ENDIF

        IF (ACTPAR(2) .EQ. 'SNUE') THEN
C                           ====
         IF (NSPLOT.EQ.MAXPLOT) THEN
          PRINT *,'CMF DECODE WARNING > EXTRA PLOTOPTIONS IGNORED'
          GOTO 6
         ENDIF
         NSPLOT=NSPLOT+1
         READ (ACTPAR(3), '(I10)', ERR=99) LPSNUE(NSPLOT)
         READ (ACTPAR(4), '(I10)', ERR=99) LPSNUED(NSPLOT)
         GOTO 6
         ENDIF

        ENDIF

        IF (ACTPAR(1) .EQ. 'PLOT') THEN
C                           ====
           IF (ACTPAR(2) == 'ALPHA') THEN
              BPLOTALPHA = .TRUE.
              GOTO 6
           ENDIF
        ENDIF

        IF (ACTPAR(1) .EQ. 'PURE' .AND. ACTPAR(2) .EQ. 'COLIRAY') THEN
C                           ====                       =======
         BCOLIRAY = .TRUE.
         GOTO 6
         ENDIF

        IF (ACTPAR(1) .EQ. 'OPC') THEN
C                           ===
         OPC = ACTPAR(2)
         GOTO 6
         ENDIF

        IF (ACTPAR(1) .EQ. 'COLI') THEN
C                           ====
          IF (ACTPAR(2) .EQ. 'OUTPUT' .AND. ACTPAR(3) .EQ. 'ONLY') THEN
            CLHLP = .TRUE.
            GOTO 6
          ENDIF

          IF (ACTPAR(2) .EQ. 'PLOT') THEN
            BPLOT = .TRUE.
            IF (NPAR .GE. 3) THEN
              READ (ACTPAR(3),'(I10)', ERR=99) IPLOT
            ENDIF
            IF (NPAR .GE. 4) THEN
              READ (ACTPAR(4),'(I10)', ERR=99) LPLOT
            ENDIF
            GOTO 6

          ELSE IF (ACTPAR(2) .EQ. 'PLOT_WCHARM') THEN
            READ (ACTPAR(3),'(I10)', ERR=99) LPLOT_WCHARM

          ELSE IF (ACTPAR(2) .EQ. 'UNLUPAR') THEN
            DO IPAR=3, NPAR-1
               IF (ACTPAR(IPAR) .EQ. 'GAMMAT') THEN
                  READ (ACTPAR(IPAR+1),'(F10.0)', ERR=99) GAMMAT
               ELSE IF (ACTPAR(IPAR) .EQ. 'TAUMAX') THEN
                  READ (ACTPAR(IPAR+1),'(F10.0)', ERR=99) UNLU_TAUMAX
               ELSE IF (ACTPAR(IPAR) .EQ. 'TAUMAX2') THEN
                  READ (ACTPAR(IPAR+1),'(F10.0)', ERR=99) UNLU_TAUMAX2
               ENDIF
            ENDDO
          ENDIF
          IF (ACTPAR(2) .EQ. 'COLI+') THEN
C                             =====
            NEWWRC = 1
            bForceCOLIP = .TRUE.
            GOTO 6
          ENDIF
          IF (ACTPAR(2) .EQ. 'PLOTRANGE') THEN
            READ (ACTPAR(3),'(F10.0)', ERR=99) XLP1
            READ (ACTPAR(4),'(F10.0)', ERR=99) XLP2
            GOTO 6
          ENDIF
          IF (ACTPAR(2) .EQ. 'GAMMA') THEN
            READ (ACTPAR(3),'(F10.0)', ERR=99) GAMMACOLI
            GOTO 6
          ENDIF
        ENDIF
C***  End of Options beginning with COLI ...

        IF (KARTE(:22) == 'NO CONTINUUM ITERATION') THEN
C                          ======================
         BITCONT = .FALSE.
         GOTO 6
         ENDIF

        IF (ACTPAR(1) .EQ. 'NO' .AND. 
     >      ACTPAR(2) .EQ. 'EDDIMIX') THEN
C                           =======
         BEMIX = .FALSE.
         GOTO 6
         ENDIF

        IF (ACTPAR(1) .EQ. 'EDDIMIX' .AND. 
     >      ACTPAR(2) .EQ. 'START') THEN
C                           =======
         READ (ACTPAR(3),'(F10.0)', ERR=99) EMIXSTART
         GOTO 6
         ENDIF

        IF (ACTPAR(1) .EQ. 'EDDIMIX' .AND. 
     >      ACTPAR(2) .EQ. 'FIX') THEN
C                           =======
         BEMIXFIX = .TRUE.
         IF (NPAR .GT. 2) THEN
           READ (ACTPAR(3),'(F10.0)', ERR=99) EMIXFIX
         ENDIF
         GOTO 6
         ENDIF

        IF (ACTPAR(1)(:16) == 'IRONLINES-EXPFAC') THEN
C                              ================
           IF (ACTPAR(2) == 'OFF') THEN
              IVERS_FE_EXPFAC = 0
              WRITE (0,*) '*** WARNING: non-standard ' // KARTE
           ELSEIF (ACTPAR(2) == 'TEFF') THEN
              IVERS_FE_EXPFAC = 1
           ELSEIF (ACTPAR(2) == 'TRADITIONAL') THEN
              IVERS_FE_EXPFAC = 2
              WRITE (0,*) '*** WARNING: non-standard ' // KARTE
           ELSE
              WRITE (0,*) '*** WARNING: Invalid option ' // KARTE
           ENDIF
         GOTO 6
         ENDIF

        !Force calculation of force multipliers
        ! Note: Usually this is automatically forced from the HYDRO interval option
        !       By using the FORCEMULTIPLIERS card you can ensure that the calculation
        !       is done at every COLI even if the HYDRO calculations are never performed
        !Full syntax is FORCEMULTIPLIERS but ACTPAR currently only reads up to 14 chars
        IF (ACTPAR(1)(1:9) == 'FORCEMULT') THEN
C                              =========
         bKALPHA = .TRUE.
         IF (NPAR > 2) THEN
           DO I=2, NPAR
             IF ((ACTPAR(I) == 'VMOD') .AND. (NPAR >= (I+1))) THEN
               READ (ACTPAR(I+1), '(F20.0)', IOSTAT=IERR) tempREAL
               IF (IERR == 0) THEN
                 VELOMODFAK = tempREAL
               ENDIF               
             ENDIF
           ENDDO
         ENDIF
         GOTO 6
        ENDIF

        IF (ACTPAR(1) == 'HYDRO') THEN
C                         =====
         bHYDROSOLVE = .TRUE.
         IF (NPAR > 2) THEN
           DO I=2, NPAR
             IF ((ACTPAR(I) == 'VMOD') .AND. (NPAR >= (I+1))) THEN
               READ (ACTPAR(I+1), '(F20.0)', IOSTAT=IERR) tempREAL
               IF (IERR == 0) THEN
                 VELOMODFAK = tempREAL
               ENDIF               
             ENDIF
           ENDDO
         ENDIF
         GOTO 6
        ENDIF

        IF (ACTPAR(1) .EQ. 'XJLAPP' .AND.
     >         ACTPAR(2) .EQ. 'NEW') THEN
C                              ===========
         IF (NPAR .GT. 2) THEN
           READ (ACTPAR(3),'(F10.0)', ERR=99) XLAM_FINE_START
           READ (ACTPAR(4),'(F10.0)', ERR=99) XLAM_FINE_END
         ENDIF 
        ENDIF

        IF (ACTPAR(1) .EQ. 'XJCAPP' .AND.
     >         ACTPAR(2) .EQ. 'NEW') THEN
C                              ===========
         IF (NPAR .GT. 2) THEN
           READ (ACTPAR(3),'(F10.0)', ERR=99) XLAM_FINE_START
           READ (ACTPAR(4),'(F10.0)', ERR=99) XLAM_FINE_END
         ENDIF 
        ENDIF

      GOTO 6


C***  END-OF-DATA REACHED: REGULAR EXIT
  100 CONTINUE
      CLOSE (1)

C***  option LINES ... either specified, or by default ALL lines
         DO 5 IND=IND1,IND2
         NLINE=NLINE+1
         IF (NLINE .GT. MAXIND) THEN
            PRINT *,' ERROR STOP IN DECCMF: NLINE .GT. MAXIND'
            CALL REMARK ('NLINE .GT. MAXIND')
            STOP 'ERROR'
            ENDIF
         LINE(NLINE) = IND
    5    CONTINUE

C***  option DRLINES ... either specified, or by default NO drlines
         DO 11 IND=INDDR1,INDDR2
         NLINE=NLINE+1
         IF (NLINE .GT. MAXIND) THEN
            PRINT *,' ERROR STOP IN DECCMF: NLINE .GT. MAXIND'
            CALL REMARK ('NLINE .GT. MAXIND')
            STOP 'ERROR'
            ENDIF
         LINE(NLINE) = IND
   11    CONTINUE


      RETURN

C***  ERROR EXIT:  ************************************************
   99 CONTINUE
      WRITE (hCPR,*) 'ERROR WHEN DECODING THE FOLLOWING LINE:'
      WRITE (hCPR,*) KARTE
      STOP 'ERROR'

      END
