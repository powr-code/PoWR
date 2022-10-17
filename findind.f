      SUBROUTINE FINDIND (IND, INDSTR, LEVEL, N, INDNUP,INDLOW,LASTIND, 
     >                    BFINDERR, LOW, NUP)
C***********************************************************************
C***  FIND LINE-INDEX "IND" FROM GIVEN LEVEL-NAMES
C***  - ACTION: IF INDEX-FIELD "INDSTR" CONTAINS A "?":
C***    READ NEXT INPUT LINE, WHICH MUST HAVE THE FOLLOWING FORMAT:
C***    UPPERLEVEL=.......... LOWERLEVEL=..........
C***    ELSE: READ IND FROM "INDSTR"
C***  - CALLED FROM: PREFORM < DECFORM < FORMAL
C***********************************************************************

      CHARACTER INDSTR*4, KARTE*80, ERRMES*80, LEVEL(N)*10
      DIMENSION INDNUP(LASTIND), INDLOW(LASTIND)
      LOGICAL BFINDERR

      BFINDERR = .FALSE.

      IF (INDSTR(1:1) .NE. '?' .AND. INDSTR(2:2) .NE. '?' .AND.
     $    INDSTR(3:3) .NE. '?' .AND. INDSTR(4:4) .NE. '?') THEN

C***     INDEX IS ASSUMED TO BE GIVEN AS NUMBER IN THE "INDSTR" FIELD
         READ (INDSTR, '(I4)', ERR=110) IND
         NUP = INDNUP(IND)
         LOW = INDLOW(IND)
         IF (IND .LE. 0 .OR. IND .GT. LASTIND) THEN
            ERRMES = 'LINE INDEX OUT-OF-RANGE: ' // INDSTR
            GOTO 100
         ENDIF

      ELSE

C***     READ NEXT INPUT LINE WHICH MUST GIVE THE LEVEL NAMES
         READ (2, '(A)') KARTE
         IF (KARTE( 1:10) .NE. 'UPPERLEVEL' .OR.
     $       KARTE(23:32) .NE. 'LOWERLEVEL') THEN
            ERRMES = 'INCORRECT INPUT CARD: ' // KARTE(:43)
            GOTO 100
         ENDIF

C***     FIND UPPER INDEX
         DO 2 I=2, N
            IF (KARTE(12:21) .EQ. LEVEL(I)) THEN
               NUP = I
               GOTO 3
            ENDIF
    2    CONTINUE
         ERRMES = 'UPPER LEVEL NOT FOUND: ' // KARTE(12:21)
         GOTO 100
    3    CONTINUE

C***     FIND LOWER INDEX
         DO 4 I=1, N
         IF (KARTE(34:43) .EQ. LEVEL(I)) THEN
            LOW = I
            GOTO 5
            ENDIF
    4    CONTINUE
         ERRMES = 'LOWER LEVEL NOT FOUND: ' // KARTE(34:43)
         GOTO 100
    5    CONTINUE

C***     FIND LINE INDEX
         IF (LOW .EQ. NUP) THEN
C***        lower and upper level have the same index
            IND = -1
            GOTO 7
         ELSE
            DO 6 I=1, LASTIND
            IF (INDLOW(I) .EQ. LOW .AND. INDNUP(I) .EQ. NUP) THEN
               IND = I
               GOTO 7
            ENDIF
    6       CONTINUE
         ENDIF

         WRITE (0,'(A)') 'LINE INDEX NOT FOUND FROM THE FOLLOWING CARD:'
         ERRMES = KARTE(:IDX(KARTE))
         GOTO 100

    7    CONTINUE
      ENDIF

      GOTO 200

C***  ERROR EXIT: 
  110 WRITE (0,'(A)') '*** ERROR: LINE Index entry invalid:'
      ERRMES = KARTE(:IDX(KARTE))
    
  100 CONTINUE
      BFINDERR = .TRUE.
      PRINT 101, ERRMES
  101 FORMAT ('THIS LINE or MULTIPLET SKIPPED! ', /, 
     >        ' NON-FATAL ERROR ',
     $          'DETECTED BY SUBROUTINE FINDIND:', /, A,/)

  200 RETURN

      END        
