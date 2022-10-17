      SUBROUTINE PLOTPOP (LEVELPL, NPLOT, N, ND, LEVEL, ENTOT,
     $                    POPNUM, MODHEAD, JOBNUM, KANAL, MAXSETS, 
     $                    BINBOX, POPRANG, NATOM, NFIRST, NLAST, NCHARG, 
     >                    SYMBOL)
C***********************************************************************
C***  DIRECT TRANSFER OF POPNUMBER PLOT
C***********************************************************************

      INTEGER, PARAMETER :: NDMAX = 210
      INTEGER, INTENT(IN) :: N, ND, JOBNUM, KANAL, MAXSETS
      INTEGER, INTENT(INOUT) :: NPLOT
      CHARACTER(100), INTENT(IN) :: MODHEAD
      CHARACTER LABELPOS*8
      REAL, DIMENSION(NDMAX) :: X, Y
      REAL, DIMENSION(ND) :: ENTOT
      REAL, DIMENSION(ND,N) :: POPNUM
      CHARACTER(70) :: HEADER
      CHARACTER(50) :: MORE
      CHARACTER(10) :: ACTLEV
      CHARACTER(8) :: CENTER, BUFFER8
      CHARACTER(10), DIMENSION(N) :: LEVEL
      CHARACTER(10), DIMENSION(MAXSETS,N) :: LEVELPL
      CHARACTER(2), DIMENSION(N) :: SYMBOL
      INTEGER, DIMENSION(MAXSETS) :: LI, ISYMBOL
      LOGICAL :: BINBOX, BPOPALL
      INTEGER, DIMENSION(NATOM) :: NFIRST, NLAST
      INTEGER, DIMENSION(N) :: NCHARG

      INTEGER :: L, IPLOT, NEWNPLOT
      REAL :: ALOGMIN

      IF (ND > NDMAX) THEN
        WRITE (0,*) 'WARNING : FROM PLOTPOP'
        WRITE (0,*) 'WARNING : PLOTNUMBERS COULD NOT BE PLOTTED'
        WRITE (0,'(A,2(1XI3))') 
     >         'NDMAX INSUFFICIENT : ND, NDMAX=', ND, NDMAX
        RETURN
      ENDIF

C***  DEFINE PLOTSYMBOLS
      ISYMBOL(1) = 1
      ISYMBOL(2) = 2
      ISYMBOL(3) = 3
      ISYMBOL(4) = 4
      ISYMBOL(5) = 8
      ISYMBOL(6) = 11
      ISYMBOL(7) = 15
      ISYMBOL(8) = 21
      ISYMBOL(9) = 22
      ISYMBOL(10)= 23
 
      CALL JSYMSET ('G2','TRANSFER')
      CALL REMARK ('POP PLOT -- DATA TO BE ROUTED')

      CENTER = CHAR(92) // 'CENTER' // CHAR(92)

C***  MINIMUM VALUE TO PREVENT FLOATING EXCEPTION WITH ALOG
      ALOGMIN = 1.E-300

      DO L=1,ND
        X(L)=ALOG10(ENTOT(L))
      ENDDO

      XMIN = FLOAT(INT(X(1)))      
      XMAX = 1. + FLOAT(INT(X(ND)))
      YMIN=POPRANG
      YMAX=.0
      XSCALE=0.
      YSCALE=0.
      XTICK=1.
      YTICK=1.
      XABST=3.
      YABST=5.
 
C***  Check for GROUNDSTATES option
      NEWNPLOT = NPLOT
      DO IPLOT=1, NPLOT
         IF (LEVELPL(1,IPLOT) .EQ. 'GROUNDSTAT') THEN
C***        Loop over all elements
            DO J=1, NATOM
C***           Restricted to one element? 
               IF (LEVELPL(2,IPLOT) .NE. ' ' .AND. 
     >             LEVELPL(2,IPLOT) .NE. SYMBOL(J)) CYCLE
               NEWNPLOT = NEWNPLOT + 1            
               IF (NEWNPLOT .GT. N) THEN 
                 WRITE (0,*) '*** TOO MANY PLOT POP CARDS'   
                 STOP  '*** FATAL ERROR IN Subr. PLOTPOP'
               ENDIF
               LEVELPL(1,NEWNPLOT) = LEVEL(NFIRST(J))
                  NPLEVEL = 1
                  DO K = NFIRST(J)+1, NLAST(J)
                     IF (NCHARG(K) .NE. NCHARG(K-1)) THEN
                        NPLEVEL = NPLEVEL + 1
                        IF (NPLEVEL .GT. MAXSETS) EXIT
                        LEVELPL(NPLEVEL,NEWNPLOT) = LEVEL(K)
                     ENDIF
                  ENDDO
            ENDDO
         ENDIF
      ENDDO
      NPLOT = NEWNPLOT

C***  Loop over all PLOT POP cards
      DO 10 IPLOT=1, NPLOT
         BPOPALL = .FALSE.
         ILAST = 0
C***     DECODE OR FIND LEVEL INDICES
         III = 0
         I4 = 0
C***     Loop over the levels to be plotted
         DO 11 II=1, MAXSETS
            I4 = I4 + 1
            ACTLEV = LEVELPL(I4,IPLOT)
            LI(II) = -1
            READ (ACTLEV(:IDX(ACTLEV)), '(I3)', ERR=12) ITEMP
            III = III + 1
            LI(III) = ITEMP
         CYCLE
C***        LEVEL NAMES GIVEN INSTEAD OF INDICES --> FIND INDICES
   12       CONTINUE
            DO 13 I=1, N
               IF (.NOT. BPOPALL) THEN
                  IF (ACTLEV .EQ. LEVEL(I)) THEN
                     III = III + 1
                     LI(III) = I
                     GOTO 11
                  ENDIF
                  IF (ACTLEV .EQ. 'ALL') BPOPALL = .TRUE.
               ELSE
                  ID = IDX(ACTLEV)
                  IF (ACTLEV( :ID) .EQ. LEVEL(I)( :ID) .AND. 
     >                      I .GT. ILAST) THEN
                     I4 = I4 - 1
                     III = III + 1
                     LI(III) = I
                     ILAST = I
                     GOTO 11
                  ENDIF
               ENDIF
   13       CONTINUE
   11    CONTINUE

C**   No single level found --> skip this plot
      IF (LI(1) .LE. 0 ) GOTO 10

      WRITE (UNIT=BUFFER8, FMT='(I7)') JOBNUM
      WRITE (UNIT=HEADER, FMT=4) MODHEAD(13:32), ADJUSTL(BUFFER8)
    4 FORMAT (A20,1X,'J',A8)
      !HEADER (1:20)=MODHEAD(13:32)
      IF (LI(1) .GT. 0) HEADER(30:39)=LEVEL(LI(1))
      IF (LI(2) .GT. 0) HEADER(41:50)=LEVEL(LI(2))
      IF (LI(3) .GT. 0) HEADER(52:61)=LEVEL(LI(3))
      IF (LI(4) .GT. 0 .OR. LI(5) .GT. 0 .OR. LI(6) .GT. 0 .OR. 
     >    LI(7) .GT. 0) HEADER(63:65)='etc'
      WRITE (KANAL, '(A)') 'PLOT   :'//HEADER

C***  MORE THAN THREE LEVELS TO BE PLOTTED: CONTINUE HEADER VERTICALLY
      MORE = ' '
      IF (LI(4) .GT. 0) MORE( 1:10) = LEVEL(LI(4))
      IF (LI(5) .GT. 0) MORE(12:21) = LEVEL(LI(5))
      IF (LI(6) .GT. 0) MORE(23:32) = LEVEL(LI(6))
      IF (LI(7) .GT. 0) MORE(34:43) = LEVEL(LI(7))
      IF (LI(8) .GT. 0) MORE(44:49) = ', etc.'

      IF (LI(4) .GT. 0 .OR. LI(5) .GT. 0 .OR. LI(6) .GT. 0 .OR. 
     >    LI(7) .GT. 0) WRITE (KANAL, '(2A)') 
     >  'KASDEF LUNA XMAX YMAX 0.7 0. 0.4 -90. ', MORE 

      IF (BINBOX) THEN
        WRITE (KANAL, '(A)') 'KASDEF INBOX'
      ELSE
        WRITE (KANAL, '(A)') '* KASDEF INBOX'
      ENDIF

C***  MARKIERUNG DER TIEFENPUNKTE 10, 20, 30, USW. 
      DO 40 L=10,ND,10
   40 WRITE (KANAL,41) X(L), YMAX, X(L), YMIN, X(L), YMIN, L
   41 FORMAT ('KASDEF LINREL ', F7.3, ' ', F7.3, ' 0. -0.5', /,
     >        'KASDEF LINREL ', F7.3, ' ', F7.3, ' 0.  0.5', /,
     >        'KASDEF LUN    ', F7.3, ' ', F7.3, ' -0.2 0.7 0.3 ', I3)


      DO 3 L=1,ND
    3 Y(L)=ALOG10(MAX(ALOGMIN,POPNUM(L,LI(1))))

      CALL PLOTANF (KANAL,HEADER,HEADER
     $ ,CENTER//'log(&Rn&N&Ttot&M/cm&H-3&M)'
     $ ,CENTER//'log(&Rn&N&Ti&M/&Rn&N&Ttot&M)'
     $ ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $ ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $ ,X,Y,ND, ISYMBOL(1))

      BACKSPACE (KANAL)
      IF (Y(1) .LT. -3.) THEN
         LABELPOS = ' 0 D+0.3'
      ELSE
         LABELPOS = ' 0 U-0.3'
      ENDIF
      WRITE (KANAL, '(A,F8.3,1X,F8.3,A)') 'KASDEF LUN ', X(1), Y(1), 
     >       LABELPOS // ' .4 &4' // LEVEL(LI(1))
      WRITE (KANAL, '(A)') 'END'


C***  SECOND  AND FURTHER DATASETS
      DO 22 II=2, MAXSETS
        IF (LI(II) .LE. 0 ) GOTO 22


        DO 2 L=1,ND
    2   Y(L) =ALOG10(MAX(ALOGMIN,POPNUM(L,LI(II))))

        IF (II .LE. 10) THEN
          CALL PLOTCON (KANAL,X,Y,ND, ISYMBOL(II))
        ELSE
          CALL PLOTCON (KANAL,X,Y,ND, 5)
        ENDIF

        BACKSPACE (KANAL)
        IF (Y(1) .LT. -3.) THEN
           LABELPOS = ' 0 D+0.3'
        ELSE
           LABELPOS = ' 0 U-0.3'
        ENDIF
        WRITE (KANAL, '(A,F8.3,1X,F8.3,A)') 'KASDEF LUN ', X(1), Y(1), 
     >         LABELPOS // ' .4 &4' // LEVEL(LI(II))
        WRITE (KANAL, '(A)') 'END'

   22 CONTINUE

   10 CONTINUE

      RETURN
      END
