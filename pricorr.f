      SUBROUTINE PRICORR (POPNUM, POP1, LEVEL, N, ND, MODHEAD, LSPOP, 
     $                   CORMAX, RTCMAX, JOBNUM, REDUCE, 
     >                   GAMMAC, GAMMAL, GAMMAR, GAMMAD, 
     >                   TOLD, TNEW, EPSILON, DELTAC, SMALLPOP, BUNLU, 
     $                   DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB,
     >                   BTDIFFUS, TNDCORR, HNDCORFAC, GAHIST, MAXGAHIST, 
     >                   STHLP, IWARN_NEG_XJCAPP, IWARN_NEG_XJLAPP,
     >                   TBTAU, TAUINT)
C*******************************************************************************
C***  PRINTOUT OF CORRECTION FACTORS OF POPULATION NUMBERS
C***  RELATIVE TO THE LAST ITERATION
C***  Called from: STEAL, EXTRAP
C*******************************************************************************
 
      DIMENSION POPNUM(ND,N),  POP1(ND,N)
      DIMENSION TOLD(ND), TNEW(ND)
      DIMENSION GAHIST(26,MAXGAHIST)
      CHARACTER LEVEL(N)*10
      CHARACTER MODHEAD*100
      CHARACTER PRILINE*130,NUMBER*12
      CHARACTER*10 LEVMAX,LEVMIN,LEVMA2,LEVMI2
      LOGICAL BUNLU, BTDIFFUS, STHLP
C***  COMMON /GIIIERR/ COUNTS THE ERROR MESSAGES FROM SUBR. GAUNTFF
      COMMON / GIIIERR /  NTUP,NTLOW,NFUP,NFLOW,NHELP
C***  NEGINTL COUNTS THE ERROR MESSAGES  FROM SUBR. LINPOP
      COMMON / COMNEGI / NEGINTL,INEGMIN,INEGMAX,LNEGMIN,LNEGMAX
C***  Warning counter for negative continuum opacities <-- LINPOP
      COMMON / COMNEGC / NEGINTC,INEGMINC,INEGMAXC,LNEGMINC,LNEGMAXC
C***  ITWARN COUNTS THE NOT CONVERGED NEWTON-RAPHSON ITERATIONS FROM
C***  SUBR. LINPOP
      COMMON / COMITWA / ITWARN, ITMAX

 
C*********************************************************
C***  SET NDMIN = FIRST DEPTH INDEX WHICH IS ACCOUNTED FOR
C***               IN THE CONVERGENCE CRITERION
C**********************************************************
C***      NDMIN = 9
C***  New Version, Lars 30-Jul-1998 10:15:26
      NDMIN = ND / 8 + 1

C**********************************************************************
C***  SET SMALLPOP : SMALLER POPNUMBERS ARE NOT ACCOUNTED FOR
C***                 IN THE CONVERGENCE CRITERIUM!
C***  NOTE: THIS CONVERGENCE CRITERION SHOULD HARMONIZE WITH THE
C***        CONVERGENCE CRITERION IN SUBR. LINPOP !!
C**********************************************************************

C**********************************************************************
C***  TABLE OF RELATIVE CORRECTIONS (IF REQUESTED)
C**********************************************************************
      PRINT 1,MODHEAD,JOBNUM
    1 FORMAT ( /, 1X, A, 15X, 'JOB NO.', I7, //, 10X,
     $ 'RATIO OF POPULATION NUMBERS TO THOSE OF THE LAST ITERATION:')
      IF (LSPOP.LT.1) GOTO 8
      J1=1
    4 J2=MIN0(N,J1+9)
      PRINT 2,(LEVEL(J),J=J1,J2)
    2 FORMAT (//,' DEPTH',10(2X,A10))
      IF (N.LT.J2) PRINT 5
    5 FORMAT (1X)
      PRINT 5
 
      DO 3 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 3
      ENCODE (130,6,PRILINE) L
    6 FORMAT (I6)
 
      DO 12 J=J1,J2
      IF (POP1(L,J) .NE. .0) THEN
            ENCODE (12,11,NUMBER) POPNUM(L,J)/POP1(L,J)
   11       FORMAT (F12.4)
            ELSE
            NUMBER='    INFINITE'
            ENDIF
      I=7+(J-J1)*12
      PRILINE(I:I+11)=NUMBER
   12 CONTINUE
 
      PRINT 13,PRILINE
   13 FORMAT (A)
    3 CONTINUE
 
      IF (J2.EQ.N) GOTO 8
      J1=J1+10
      GOTO 4
    8 CONTINUE
 
C***  PRINTOUT OF THE TEMPERATURE CORRECTIONS
      IF (BUNLU .AND. LSPOP .GT. 0) THEN
      PRINT 21
   21 FORMAT (////,40X,'RELATIVE TEMPERATURE CORRECTIONS',/,40X,32('='),
     $   //,40X,'DEPTH INDEX      REL.TEMP.CORR.',/)
      DO 22 L=1,ND
      IF(((L-1)/LSPOP)*LSPOP.NE.(L-1) .AND. L.NE.ND) GOTO 22
      RTC=(TNEW(L)-TOLD(L))/TOLD(L)
      PRINT 23,L,RTC
   23 FORMAT (45X,I3,F15.4)
   22 CONTINUE
      ENDIF
 

C**********************************************************************
C***  CALCULATE MINIMUM AND MAXIMUM OF ALL CORRECTION FACTORS
C**********************************************************************
      QMAX=1.
      QMA2=1.
      QMIN=1.
      QMI2=1.
      LEVMAX='NONE'
      LEVMA2='NONE'
      LEVMIN='NONE'
      LEVMI2='NONE'
      LMAX=0
      LMA2=0
      LMIN=0
      LMI2=0
      NSMALL = 0

      DO 7 L=NDMIN,ND
      DO 7 J = 1, N
      POP1LJ=POP1(L,J)
C***  TOO SMALL OR NEGATIVE POPNUMBERS ARE CONSIDERED AS BEEING = SMALLPOP  **
      IF (POP1LJ .LT. SMALLPOP) POP1LJ = SMALLPOP
      POPLJ=POPNUM(L,J)
C***  POPNUMBER LESS THAN SMALLPOP ARE NOT ACCOUNTED FOR THE CALCULATION
C***    OF THE MAXMIMUM CORRECTION
      IF (POPLJ .LT. SMALLPOP) THEN
        POPLJ = SMALLPOP
        NSMALL = NSMALL + 1
      ELSE
        Q=POPLJ/POP1LJ
        IF (Q .GT. QMA2) THEN
           IF (Q .GT. QMAX) THEN
              QMA2=QMAX
              QMAX=Q
              LEVMA2=LEVMAX
              LEVMAX=LEVEL(J)
              LMA2=LMAX
              LMAX=L
              JMAX = J
           ELSE
              QMA2=Q
              LEVMA2=LEVEL(J)
              LMA2=L
              JMA2 = J
           ENDIF
        ENDIF
        IF (Q .LT. QMI2) THEN
           IF (Q .LT. QMIN) THEN
              QMI2=QMIN
              QMIN=Q
              LEVMI2=LEVMIN
              LEVMIN=LEVEL(J)
              LMI2=LMIN
              LMIN=L
              JMIN = J
           ELSE
              QMI2=Q
              LEVMI2=LEVEL(J)
              LMI2=L
              JMI2 = J
           ENDIF
        ENDIF
      ENDIF
    7 CONTINUE

C***********************************************************************
      CORMAX=AMAX1(QMAX-1.,1./QMIN-1.)
C***********************************************************************

      IF (CORMAX .GT. .0) CORLOG=ALOG10(CORMAX)
      PRINT 9,QMAX,LEVMAX,LMAX,QMA2,LEVMA2,LMA2,QMIN,LEVMIN,LMIN,QMI2,
     $        LEVMI2,LMI2,NDMIN,ND,CORMAX,CORLOG
    9 FORMAT (/,15X,'MAX:',F8.4,'  (',A10,'  L=',I3,')',
     $           7X,'2ND:',F8.4,'  (',A10,'  L=',I3,')',
     $        /,15X,'MIN:',F8.4,'  (',A10,'  L=',I3,')',
     $           7X,'2ND:',F8.4,'  (',A10,'  L=',I3,')',
     $        /,15X,'DEPTH POINTS CONSIDERED:',I5,'  TO',I5,
     $          15X,'CORMAX=',F7.4,5X,'LOG=',F6.2,//)

      WRITE (0, '(A,I7,A,F6.2,A)') 'JOBNUM=', JOBNUM, 
     >       '   log max correction =', CORLOG, ' <<<<<<<'
 
C***  MAXIMUM TEMPERATURE CORRECTION
      IF (BUNLU) THEN
        RTCMAX=.0
        LRTCMAX=0
        DO L=1,ND
          RTC=(TNEW(L)-TOLD(L))/TOLD(L)
          IF (ABS(RTC) .GT. ABS(RTCMAX)) THEN
            RTCMAX=RTC
            LRTCMAX=L
          ENDIF
        ENDDO

        PRINT 35, RTCMAX, LRTCMAX, 
     >            DUNLU_LOC, DUNLU_INT, DUNLU_RMAX, DUNLU_TB
   35   FORMAT (15X,'MAX.REL.TEMP.CORR.=',F8.4,'  AT L=',I3,
     >     '    (UNLU SCALINGS LOC, INT, RMAX, TB=', 4(F5.3,1X), ')')

        IF (DUNLU_TB .GT. .0) PRINT 36, TBTAU 
   36    FORMAT 
     >   (15X, 'Thermal Balance term suppressed for tau_Ross >', F6.3)

        IF (TAUINT > .0) PRINT 37, TAUINT
   37    FORMAT 
     >   (15X, 'INT and RMAX term damped depth-dependent, TAUINT=', F5.1)

        PRINT *     

      ENDIF
         
      IF (BTDIFFUS) 
     >   WRITE (*,'(15X,2A,E10.3,/)')
     >    'TEMPERATURE AT INNER BOUNDARY ADJUSTED ',
     >    'BY DIFFUSION APPROXIMATION : REL.CORR.= ', TNDCORR
 
      IF (.NOT. BTDIFFUS .AND. BUNLU) 
     >   WRITE (*, '(15X,2A,E10.3,/)') 
     >    'TEMPERATURE GRADIENT AT INNER BOUNDARY CORRECTED ',
     >    'FOR FLUX CONSISTENCY : REL.CORR.= ', HNDCORFAC-1.
     
      IF (GAMMAC.GT..0 .OR. GAMMAL.GT..0)
     $   PRINT 14, GAMMAC, GAMMAL, DELTAC, GAMMAR, GAMMAD
   14    FORMAT (10X,'SCHARMER-PARAMETERS:',  10X,
     $          'GAMMAC=',F9.1,'   GAMMAL=',F9.1,/,
     $     40X,'DELTAC=',F9.1,'   GAMMAR=',F9.1,'   GAMMAD=',F9.1,//)
 
      IF (REDUCE .NE. 1.) PRINT 10,REDUCE
   10 FORMAT (10X,'CORRECTIONS REDUCED BY FACTOR',F5.2,//)
 

C**********************************************************************
C***  PRINTOUT OF WARNINGS AND ERROR MESSAGES  ************************
C**********************************************************************
      IF (NSMALL .GT. 0) PRINT 27, NSMALL, SMALLPOP
   27 FORMAT ( 8X, I6, 
     >  ' WARNINGS: POP. NUMBERS .LT.', 1PE8.1, 
     >  ' NOT ACCOUNTED FOR IN THE MAX. CORRECTIONS')

C***  PRINTOUT OF WARNINGS FROM SUBROUTINE LINPOP  *********************

      IF (NEGINTL .GT. 0)  
     $   PRINT 28, NEGINTL,INEGMIN, INEGMAX, LNEGMIN, LNEGMAX
   28 FORMAT ( 8X,I6,' WARNINGS: NEGATIVE LINE INTENSITIES REPLACED',
     $   ' BY ZERO  --  ', I4, ' < IND <', I4, ';', I5, ' < L <', I3)

      IF (NEGINTC .GT. 0)  
     $   PRINT 29, NEGINTC,INEGMINC, INEGMAXC, LNEGMINC, LNEGMAXC
   29 FORMAT ( 8X,I6,' WARNINGS: NEGATIVE CONT INTENSITIES REPLACED',
     $   ' BY ZERO  --  ', I4, ' <  K  <', I4, ';', I5, ' < L <', I3)

      IF (ITWARN .GT. 0) PRINT 33, ITWARN, ITMAX
   33 FORMAT ( 8X,I6,' WARNINGS: NEWTON-RAPHSON ITERATION ',
     $      'NOT CONVERGED (MAX. NUMBER OF ITERATIONS: ',I3,')')
 

C***  PRINTOUT OF WARNINGS FROM SUBROUTINE LINPOP  *********************

      IF (IWARN_NEG_XJCAPP .GT. 0) THEN
        write (0,*) 'Neg XJCAPP in XJAPP:', IWARN_NEG_XJCAPP
        PRINT 40, IWARN_NEG_XJCAPP
   40   FORMAT ( 8X,I6,' WARNINGS: Negative XJCAPP appeared in XJAPP')
      ENDIF

      IF (IWARN_NEG_XJLAPP .GT. 0) THEN
        write (0,*) 'Neg XJLAPP in XJAPP:', IWARN_NEG_XJLAPP
        PRINT 41, IWARN_NEG_XJLAPP
   41   FORMAT ( 8X,I6,' WARNINGS: Negative XJLAPP appeared in XJAPP')
      ENDIF

C***  PRINTOUT OF ERROR MESSAGES FROM SUBROUTINE GAUNTFF

      IF (CORMAX .LT. EPSILON) THEN
   26 FORMAT (7X,I7,' WARNINGS: CALLS OF GAUNTFF BEYOND ',A,' BOUND')
        IF (NTUP .GT.0) PRINT 26, NTUP , 'UPPER TEMPERATURE'
        IF (NTLOW.GT.0) PRINT 26, NTLOW, 'LOWER TEMPERATURE'
        IF (NFUP .GT.0) PRINT 26, NFUP , 'UPPER FREQUENCY'
        IF (NFLOW.GT.0) PRINT 26, NFLOW, 'LOWER FREQUENCY'
        ENDIF
 
C***  Update GAMMA HISTORY
      IF (.NOT. STHLP) THEN
        GAHIST(1,1)  = FLOAT(JOBNUM)
        GAHIST(3,1)  = GAMMAC
        GAHIST(4,1)  = GAMMAL
        GAHIST(5,1)  = GAMMAR
        GAHIST(6,1)  = GAMMAD
        GAHIST(7,1)  = FLOAT(LMAX)
        GAHIST(8,1)  = FLOAT(JMAX)
        GAHIST(9,1)  = QMAX
        GAHIST(10,1) = FLOAT(LMA2)
        GAHIST(11,1) = FLOAT(JMA2)
        GAHIST(12,1) = QMA2
        GAHIST(13,1) = FLOAT(LMIN)
        GAHIST(14,1) = FLOAT(JMIN)
        GAHIST(15,1) = QMIN
        GAHIST(16,1) = FLOAT(LMI2)
        GAHIST(17,1) = FLOAT(JMI2)
        GAHIST(18,1) = QMI2
        GAHIST(19,1) = CORMAX
        GAHIST(20,1) = LRTMAX
        GAHIST(21,1) = RTCMAX
      ENDIF

      RETURN
      END
