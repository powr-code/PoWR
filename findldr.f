      SUBROUTINE FINDLDR (LOW, NUP, LEVEL, NOM, NCHARG, N, BFINDERR)
C***********************************************************************
C***  FIND LEVEL-INDICES FROM GIVEN LEVEL-NAMES
C***  - READ NEXT INPUT LINE, WHICH MUST HAVE THE FOLLOWING FORMAT:
C***    UPPERLEVEL=.......... LOWERLEVEL=..........
C***  - CALLED FROM: FORMAL > PREFORM
C***********************************************************************

      DIMENSION NOM(N), NCHARG(N)
      CHARACTER KARTE*80, ERRMES*80, LEVEL(N)*10
      LOGICAL BFINDERR

      BFINDERR = .FALSE.

C***  READ NEXT INPUT LINE WHICH MUST GIVE THE LEVEL NAMES
      READ (2, '(A)') KARTE
      IF (KARTE( 1:10) .NE. 'UPPERLEVEL' .OR.
     $    KARTE(23:32) .NE. 'LOWERLEVEL') THEN
            ERRMES = 'INCORRECT INPUT CARD: ' // KARTE(:43)
            GOTO 100
      ENDIF

C***  FIND UPPER INDEX
      DO 2 I=2, N
         IF (KARTE(12:21) .EQ. LEVEL(I)) THEN
            NUP = I
            GOTO 3
         ENDIF
    2 CONTINUE
      ERRMES = 'UPPER LEVEL NOT FOUND: ' // KARTE(12:21)
      GOTO 100

    3 CONTINUE

C***  FIND LOWER INDEX
      DO 4 I=1, N
         IF (KARTE(34:43) .EQ. LEVEL(I)) THEN
            LOW = I
            GOTO 5
         ENDIF
    4 CONTINUE
      ERRMES = 'LOWER LEVEL NOT FOUND: ' // KARTE(34:43)
      GOTO 100

    5 CONTINUE

C***  CHECK FOR CORRECT IONIZATION STAGES: CHARGE DIFFERENCE = 1
      IF ((NOM(LOW) .NE. NOM(NUP)) .OR. 
     >    (NCHARG(LOW) .NE. NCHARG(NUP)-1)) THEN
         ERRMES = '*** ERROR: INVALID COMBINATION OF LEVELS'
         GOTO 100
      ENDIF

      RETURN

C***  ERROR EXIT
  100 BFINDERR = .TRUE.

      PRINT 101, ERRMES
      WRITE (0,101) ERRMES
  101 FORMAT (/, 1X, 80('*'), /, ' NON-FATAL ERROR ON INPUT OPTIONS, ',
     $        'DETECTED BY SUBROUTINE FINDLDR:', /,
     $        1X, A, /, 1X, 80('*'), /)
      PRINT *, KARTE(:IDX(KARTE))

      RETURN

      END        
