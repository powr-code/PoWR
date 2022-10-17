      SUBROUTINE PLOTT (TPLOTOPT,ND,R,TAUROSS,T,MODHEAD,JOBNUM,KANAL,
     >                  BINBOX,UNLUTECLINE)
C******************************************************************************
C***  DIRECT TRANSFER OF THE TEMPERATURE STRATIFICATION
C***  DEFAULT: T(R) VERSUS LOG(R/R*-1)
C***  TPLTAU=.TRUE.: T(R) VERSUS LOG(TAUROSS)
C******************************************************************************
 
      IMPLICIT NONE

      INTEGER, PARAMETER :: NDMAX = 210 
      INTEGER, INTENT(IN) :: ND, JOBNUM, KANAL

      REAL, DIMENSION(NDMAX) :: X, Y
      REAL, DIMENSION(ND) :: R, TAUROSS, T
      LOGICAL :: TPLTAU, BINBOX, BCOMPARTIBLE
      CHARACTER(70) :: HEAD1, HEAD2
      CHARACTER(100) :: MODHEAD
      CHARACTER(8) :: CENTER, BUFFER8
      CHARACTER(40) :: CUROPT, NEXTOPT
      CHARACTER(*) :: UNLUTECLINE
      CHARACTER(256) :: TPLOTOPT
      INTEGER, PARAMETER :: NPARMAX = 40
      CHARACTER(40), DIMENSION(NPARMAX) :: ACTPAR

      INTEGER :: I, J, L, L2, LMIN, LMAX, NPAR, NEXTPARS, IERR, IPAR, 
     >           NDL
      REAL :: TAUMAX, tempREAL, TMIN, TAUMIN, RMAX, RMIN, TMAX,
     >        YMIN, YMAX, XMIN, XMAX, XSCALE, YSCALE, XTICK, YTICK,
     >        XABST, YABST
      INTEGER, EXTERNAL :: IDX, ISRCHFGE, ISRCHFGT, ISRCHFLE, ISRCHFLT

      XMIN = 0. 
      XMAX = 0.
      YMIN = 0.
      YMAX = 0.
      ACTPAR = '0.0'
      TPLTAU = .FALSE.
      BCOMPARTIBLE = .FALSE.    !if true, the old CARDS format is also interpreted

C***  READ CARDS LINE
      IF (BCOMPARTIBLE) THEN
        DO I=1, IDX(TPLOTOPT)
          !backward compartibility - remove all "_" 
          IF (TPLOTOPT(I:I) == '_') THEN
            TPLOTOPT(I:I) = ' '
          ENDIF
        ENDDO
      ENDIF
      CALL SARGC (TPLOTOPT, NPAR)
      IF (NPAR > 2) THEN
C***  Get actual parameters
        DO I=1, NPAR
          CALL SARGV(TPLOTOPT, I, ACTPAR(I))
        ENDDO
        DO I=3, NPAR
          CUROPT = ACTPAR(I)
          SELECTCASE (CUROPT)
            CASE ('Y-AXIS')
              IF (.NOT. BCOMPARTIBLE) CYCLE
              !this is only due to backward compartibility
              NEXTPARS = MIN(NPAR, I+4)
              DO J=I+1, NEXTPARS
                CALL SARGV (TPLOTOPT, J, NEXTOPT)
                IF (NPAR >= (J+1)) THEN
                  READ (ACTPAR(J+1), '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    IF (NEXTOPT == 'MAX') THEN
                      YMAX = tempREAL
                    ELSEIF (NEXTOPT == 'MIN') THEN
                      YMIN = tempREAL
                    ENDIF                  
                  ENDIF
                ENDIF
              ENDDO
            CASE ('X-AXIS')
              IF (.NOT. BCOMPARTIBLE) CYCLE
              NEXTPARS = MIN(NPAR, I+4)
              DO J=I+1, NEXTPARS
                CALL SARGV (TPLOTOPT, J, NEXTOPT)
                IF (NPAR >= (J+1)) THEN
                  READ (ACTPAR(J+1), '(F10.0)', IOSTAT=IERR) tempREAL
                  IF (IERR == 0) THEN
                    IF (NEXTOPT == 'MAX') THEN
                      XMAX = tempREAL
                    ELSEIF (NEXTOPT == 'MIN') THEN
                      XMIN = tempREAL
                    ENDIF                  
                  ENDIF
                ENDIF
              ENDDO
            CASE ('XMIN','XMAX','YMIN','YMAX')
              IF (NPAR >= (I+1)) THEN
                READ (ACTPAR(I+1), '(F10.0)', IOSTAT=IERR) tempREAL
                IF (IERR == 0) THEN
                  IF (CUROPT == 'XMAX') THEN
                    XMAX = tempREAL
                  ELSEIF (CUROPT == 'XMIN') THEN
                    XMIN = tempREAL
                  ELSEIF (CUROPT == 'YMAX') THEN
                    YMAX = tempREAL
                  ELSEIF (CUROPT == 'YMIN') THEN
                    YMIN = tempREAL
                  ENDIF
                ENDIF
              ENDIF
            CASE ('TAU', 'TAUR')
              TPLTAU = .TRUE.
          ENDSELECT
        ENDDO
      ENDIF


C***  READ IN UNLUTEC-CARD DO OBTAIN RECENT TMIN FOR PLOT
      TMIN = 6000.                          !default value
      CALL SARGC (UNLUTECLINE,NPAR)
      DO I=2, MIN(NPAR,NPARMAX)
       CALL SARGV(UNLUTECLINE,I,ACTPAR(I))
      ENDDO

      IF (NPAR .GT. 1) THEN
       DO IPAR = 2 , NPAR
        IF (ACTPAR(IPAR) .EQ. 'TMIN') THEN
         IF (IPAR+1 .GT. NPAR) THEN
          EXIT
         ENDIF
         READ (ACTPAR(IPAR+1),'(F10.0)', ERR=1) TMIN
        ENDIF
       ENDDO
      ENDIF


      CALL JSYMSET ('G2','TRANSFER')
      CALL REMARK ('T-STRATIFICATION TO BE ROUTED')
      CENTER = CHAR(92) // 'CENTER' // CHAR(92)
 
C***  X-AXIS: ----------------------------------------------------------
C***  TPLTAU=.TRUE.: ROSSELAND OPTICAL DEPTH ( RANGE: -3. <= LOG(TAUROSS) <= 1.)
      IF (TPLTAU) THEN
         IF (XMAX == 0.) THEN
            TAUMAX=19.
            XMAX = LOG10(TAUMAX)
C            XMAX=1.
         ELSE
            TAUMAX=10.**XMAX
         ENDIF
         IF (XMIN == 0.) THEN
            TAUMIN= MIN(0.001, TAUROSS(INT(0.2*ND)))
            XMIN = LOG10(TAUMIN)
C            XMIN=-3.
         ELSE
            TAUMIN=10.**XMIN
         ENDIF
         LMIN=ISRCHFGE(ND,TAUROSS(1),1,TAUMIN)
         LMAX=ISRCHFGT(ND,TAUROSS(1),1,TAUMAX)-1
      ELSE
C***  DEFAULT: RADIUS ( RANGE: -3.5 <= LOG(R/R*-1) <= 2. )
         IF (XMAX == 0.0) THEN
            RMAX = R(1)
            XMAX = LOG10(RMAX)
C            RMAX=101.
C            XMAX=2.
         ELSE
            RMAX=1.+10.**XMAX
         ENDIF
         IF (XMIN == 0.) THEN
             RMIN=1.000316228
             RMIN = MIN(RMIN, R(ND-3))
             XMIN = MIN(-3.5, LOG10(RMIN - 1.))
C            XMIN=-3.5
         ELSE
            RMIN=1.+10.**XMIN
         ENDIF
         LMIN=ISRCHFLT(ND,R(1),1,RMAX)
         LMAX=ISRCHFLE(ND,R(1),1,RMIN)-1
      ENDIF

C***  TRANSFER ALL DATA (WITH OPTION "KASDEF INBOX")
      IF (TPLTAU) THEN
         LMIN = 2
         LMAX = ND
      ELSE
         LMIN = 1
         LMAX = ND-1
      ENDIF
      NDL=LMAX-LMIN+1
      IF (NDL .GT. NDMAX) THEN
         CALL REMARK ('T-PLOT IMPOSSIBLE')
         RETURN
         ENDIF
      IF (TPLTAU) THEN
        DO L=1, NDL
          X(L)=LOG10(TAUROSS(LMIN-1+L))
          Y(L)=T(LMIN-1+L)/1.E3
        ENDDO
      ELSE
        DO L=1, NDL
          X(L)=LOG10(R(LMAX+1-L)-1.)
          Y(L)=T(LMAX+1-L)/1.E3
        ENDDO
      ENDIF
 
C***  Y-AXIS: TEMPERATURE RANGE  ---------------------------------------
C***  DEFAULT:     10.*(AINT(TMAX/10.)+1.) > T(KK) > 0.
C***               TMAX: MAXIMUM OF T(R) IN THE SPECIFIED X-RANGE
      IF (YMAX == 0.) THEN
         TMAX=Y(1)
         DO L=2, NDL
           IF (Y(L) > TMAX) TMAX = Y(L)
         ENDDO
         YMAX=10.*(INT(TMAX/10.)+1.)
      ENDIF
 
C***  PLOT PARAMETER
      XSCALE=0.
      YSCALE=0.
      XTICK=0.5
      YTICK=10.
      IF (YMAX >= 250) YTICK = 50.
      IF (YMAX >= 750) YTICK = 100.
      XABST=1.
      YABST=10.
      WRITE (UNIT=BUFFER8, FMT='(I7)') JOBNUM
      WRITE (UNIT=HEAD2,FMT=9) MODHEAD(13:32), ADJUSTL(BUFFER8)
    9 FORMAT (A20,1X,'J',A7,1X)
      IF (TPLTAU) THEN
         HEAD1=' WR TEMPERATURE STRATIFICATION T(R) VERSUS LOG(TAUROSS)'
         HEAD2(30:65)=' TEMPERATURE STRATIFICATION T(TAU-R)'
      ELSE
         HEAD1=' WR TEMPERATURE STRATIFICATION T(R) VERSUS LOG(R/R*-1)'
         HEAD2(30:61)=' TEMPERATURE STRATIFICATION T(R)'
      ENDIF
      !HEAD2(1:20)=MODHEAD(13:32)

C***  WRITE "INBOX" OPTION
      WRITE (KANAL,*) 'PLOT: ', HEAD1
      IF (BINBOX) THEN
        WRITE (KANAL,'(A)') 'KASDEF INBOX'
      ELSE
        WRITE (KANAL,'(A)') '* KASDEF INBOX'
      ENDIF

      DO L=1, NDL
        IF (TPLTAU) THEN
          L2 = LMIN - 1 + L
        ELSE
          L2 = LMAX + 1 - L
        ENDIF
        IF ( (MOD(L2, 10) == 0)
     >       .AND. (X(L) > XMIN) .AND. (X(L) < XMAX) ) THEN
          WRITE (KANAL,41) X(L), X(L), X(L), L2
   41     FORMAT ('KASDEF LINREL ', F7.3, ' YMAX 0. -0.5', /,
     >            'KASDEF LINREL ', F7.3, ' YMIN 0.  0.5', /,
     >            'KASDEF LUN    ', F7.3, ' YMIN -0.2 0.7 0.3 ', I3)
        ENDIF
      ENDDO

        WRITE (KANAL,'(A)') 'KASDEF COLOR=2'
        IF (T(ND) .LT. T(ND-1)) THEN
          WRITE (KANAL, 43) X(1), Y(1)
   43     FORMAT('KASDEF ARR ',2(F7.3,1X), 
     >           '1.0 225. 0.4 0.3 0. 1 FILLED')
        ENDIF
        WRITE (KANAL,'(A)') 'KASDEF COLOR=1'

C***  Additional horizontal lines at 10 and 20kK
        WRITE (KANAL, '(A)') 'KASDEF COLOR=8'
        WRITE (KANAL, '(A)') 
     >    'KASDEF LINUN XMIN 10. XMAX 10.  0. 0.'
        WRITE (KANAL, '(A)')
     >    'KASDEF LINUN XMIN 20. XMAX 20.  0. 0.'
c        WRITE (KANAL, '(A)') 'KASDEF COLOR=1'
        WRITE (KANAL, '(A)') 'KASDEF COLOR=4'
        WRITE (KANAL, '(A,F10.3,A,F10.3,A)')
     >    'KASDEF LINUN XMIN ',TMIN/1.E3,' XMAX ',TMIN/1.E3,' 0. 0.'
        WRITE (KANAL, '(A)') 'KASDEF FONT=HELVET'
        WRITE (KANAL, '(A,F10.3,A,F5.1)')
     >   'KASDEF LUN XMIN ',TMIN/1.E3,' L0.5 D0.2 0.35 T\!&TMIN&M\,=\,'
     >     ,TMIN/1.E3
        WRITE (KANAL, '(A)') 'KASDEF COLOR=1'

      IF (TPLTAU) THEN
        WRITE (KANAL, '(A)') 'KASDEF COLOR=2'
        WRITE (KANAL, '(A,F9.3,A)')
     >    'KASDEF SYM XMAX ', T(ND)/1000., ' 0. 0. 0.3 8'
        IF (R(1) > RMAX) THEN
          !draw outer circle only if outside of plot
          WRITE (KANAL, '(A,F9.3,A)')
     >      'KASDEF SYM XMIN ', T(1)/1000., ' 0. 0. 0.3 8'
        ENDIF
        WRITE (KANAL, '(A)') 'KASDEF COLOR=1'
          CALL PLOTANF (KANAL,HEAD1,HEAD2
     $         ,CENTER//'log #t#&TRosseland&M'
     $         ,CENTER//'T / kK'
     $         ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $         ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $         ,X,Y,NDL,1)
      ELSE
C***  MARKIERUNG DER TIEFENPUNKTE 10, 20, 30, USW. 
        WRITE (KANAL, '(A)') 'KASDEF COLOR=2'
        WRITE (KANAL, '(A,F9.3,A)')
     >    'KASDEF SYM XMIN ', T(ND)/1000., ' 0. 0. 0.3 8'
        IF (R(1) > RMAX) THEN
          !draw outer circle only if outside of plot
          WRITE (KANAL, '(A,F9.3,A)')
     >      'KASDEF SYM XMAX ', T(1)/1000., ' 0. 0. 0.3 8'
        ENDIF
        WRITE (KANAL, '(A)') 'KASDEF COLOR=1'
C***        

         CALL PLOTANF (KANAL,HEAD1,HEAD2
     $        ,CENTER//'log (R/R&T*&M - 1)'
     $        ,CENTER//'T / kK'
     $        ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $        ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $        ,X,Y,NDL,1)
C***
      ENDIF
 
      RETURN
1     CONTINUE
      TMIN = 6000.
      END
