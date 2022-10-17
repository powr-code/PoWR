      SUBROUTINE PLOTDEP (LEVELPLDEP, NPLOTDEP, N, ND, LEVEL, ENTOT,
     $                    DEPART, MODHEAD, JOBNUM, KANAL, MAXSETS, 
     $                    BINBOX, POPNUM, POPLTE, NATOM, NFIRST, NLAST,
     >                    NCHARG, SYMBOL)

C***********************************************************************
C***  Plotting of the non-LTE departure coefficients
C***********************************************************************

      INTEGER, PARAMETER :: NDMAX = 210
      REAL, DIMENSION(NDMAX) :: X, Y
      REAL, DIMENSION(ND) :: ENTOT
      REAL, DIMENSION(ND,N) :: DEPART, POPNUM, POPLTE
      CHARACTER HEADER*60,MODHEAD*100,LEVEL(N)*10,CENTER*8
      CHARACTER(5) :: MORE
      CHARACTER*10 LEVELPLDEP(MAXSETS,N), ACTLEV
      INTEGER, DIMENSION(MAXSETS) :: LI, ISYMBOL
      LOGICAL :: BINBOX, BPOPALL
      DIMENSION NFIRST(NATOM), NLAST(NATOM), NCHARG(N)
      CHARACTER(2), DIMENSION(N) :: SYMBOL

      IF (ND .GT. NDMAX) THEN
        WRITE (0,*) 'WARNING : FROM PLOTDEP'
        WRITE (0,*) 'WARNING : DEPART COULD NOT BE PLOTTED'
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
      CALL REMARK ('PLOT DEPART -- DATA TO BE ROUTED')

      CENTER = CHAR(92) // 'CENTER' // CHAR(92)

      DO L=1,ND
        X(L)=ALOG10(ENTOT(L))
      ENDDO

      XMIN = FLOAT(INT(X(1)))
      XMAX = 1. + FLOAT(INT(X(ND)))
      YMIN= -5.
      YMAX= +5.
      XSCALE=0.
      YSCALE=0.
      XTICK=1.
      YTICK=1.
      XABST=3.
      YABST=5.
 
C***  Check for GROUNDSTATES option
      NEWNPLOTDEP = NPLOTDEP
      DO IPLOT=1, NPLOTDEP
         IF (LEVELPLDEP(1,IPLOT) .EQ. 'GROUNDSTAT') THEN
C***        Loop over all elements
            DO J=1, NATOM
C***           Restricted to one element?
               IF (LEVELPLDEP(2,IPLOT) .NE. ' ' .AND.
     >             LEVELPLDEP(2,IPLOT) .NE. SYMBOL(J)) CYCLE
               NEWNPLOTDEP = NEWNPLOTDEP + 1
               IF (NEWNPLOTDEP .GT. N) THEN
                 WRITE (0,*) '*** TOO MANY PLOT DEPART CARDS'
                 STOP  '*** FATAL ERROR IN Subr. PLOTDEP'
               ENDIF
               !Add lowest state of element
               LEVELPLDEP(1,NEWNPLOTDEP) = LEVEL(NFIRST(J))
                  !Read additional ground states of the same element
                  NPLEVEL = 1
                  DO K = NFIRST(J)+1, NLAST(J)
                     IF (NCHARG(K) .NE. NCHARG(K-1)) THEN
                        NPLEVEL = NPLEVEL + 1
                        IF (NPLEVEL .GT. MAXSETS) EXIT
                        LEVELPLDEP(NPLEVEL,NEWNPLOTDEP) = LEVEL(K)
                     ENDIF
                  ENDDO
            ENDDO
         ENDIF
      ENDDO
      NPLOTDEP = NEWNPLOTDEP

C***  Loop over all PLOT POP cards
      plotloop: DO IPLOT=1, NPLOTDEP ! - - - - - - - - - - - - - - - - -
         BPOPALL = .FALSE.
         ILAST = 0
C***     DECODE OR FIND LEVEL INDICES
         III = 0
         I4 = 0
C***     Loop over the levels to be plotted
         DO 11 II=1, MAXSETS
            I4 = I4 + 1
            ACTLEV = LEVELPLDEP(I4,IPLOT)
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
      IF (LI(1) .LE. 0 ) CYCLE plotloop

      WRITE (HEADER, FMT="(19X,'J',I7)") JOBNUM
      HEADER (1:20)=MODHEAD(13:32)
      IF (LI(1) .GT. 0) HEADER(29:38)=LEVEL(LI(1))
      IF (LI(2) .GT. 0) HEADER(40:49)=LEVEL(LI(2))
      IF (LI(3) .GT. 0) HEADER(51:60)=LEVEL(LI(3))
      IF (LI(4) .GT. 0 .OR. LI(5) .GT. 0 .OR. LI(6) .GT. 0 .OR. 
     >    LI(7) .GT. 0) HEADER(62:64)='etc'
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
      DO L=10,ND,10
        WRITE (KANAL,41) X(L), YMAX, X(L), YMIN, X(L), YMIN, L
      ENDDO
   41 FORMAT ('KASDEF LINREL ', F7.3, ' ', F7.3, ' 0. -0.5', /,
     >        'KASDEF LINREL ', F7.3, ' ', F7.3, ' 0.  0.5', /,
     >        'KASDEF LUN    ', F7.3, ' ', F7.3, ' -0.2 0.7 0.3 ', I3)

C*** Loop over all levels to be plotted
      ISET=0
      DO II=1, MAXSETS
C***    GIVEN LEVEL INVALID?
        IF (LI(II) .LE.0 ) CYCLE
 
        ISET = ISET + 1
        DO L=1,ND
           DEP = DEPART(L,LI(II))
           IF (DEP .GT. 1.E-99) THEN
             Y(L) = ALOG10(DEP)
C***         dep. may noch have been calculated yet - try to calculate it here
           ELSE IF (POPLTE(L,LI(II)) .GT. 1.E-100) THEN 
             DEP = POPNUM(L,LI(II)) / POPLTE(L,LI(II))
             IF (DEP .GT. 1.E-99) THEN 
                Y(L) = ALOG10(DEP)
             ELSE 
                Y(L) = -99.
             ENDIF
           ELSE
             Y(L) = -99.
           ENDIF
        ENDDO

        IF (ISET .EQ. 1) THEN 
          CALL PLOTANF (KANAL,HEADER,HEADER
     $         ,CENTER//'log(&Rn&N&Ttot&M/cm&H-3&M)'
     $         ,CENTER//'log(&Rn&N&Ti&M/&Rn&N&Ti&HLTE&M)'
     $         ,XSCALE,XMIN,XMAX,XTICK,XABST,.0
     $         ,YSCALE,YMIN,YMAX,YTICK,YABST,.0
     $         ,X,Y,ND, ISYMBOL(1))
        ELSEIF (ISET .LE. 10) THEN
           CALL PLOTCON (KANAL,X,Y,ND, ISYMBOL(ISET))
        ELSE
           CALL PLOTCON (KANAL,X,Y,ND, 5)
        ENDIF
      ENDDO   ! Loop over the levels to be plotted 

      ENDDO plotloop ! - - - - - - - - - - - - - - - - - - - - - - - - -

      RETURN
      END
