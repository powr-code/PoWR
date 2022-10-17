      SUBROUTINE DECADP (OLDSTART, BDEPART, BTAUR, POPMIN, 
     >                LEVELCARD, NLEVELCARD, MAXLEVELCARD)
C*****************************************************************
C***  DEECODE INPUT OPTIONS FOR PROGRAM 'ADAPTER'
C*****************************************************************

      LOGICAL OLDSTART, BDEPART, BTAUR
      CHARACTER KARTE*80, ACTPAR1*20
      CHARACTER LEVELCARD*80(MAXLEVELCARD)

C***  DEFAULT VALUES: 
      OLDSTART=.FALSE.
      BDEPART = .FALSE. 
      BTAUR   = .FALSE.
C***  POPMIN: NULL popnums are set to this value
      POPMIN = 1.E-25
      NLEVELCARD = 0

      OPEN (1, FILE='CARDS', STATUS='UNKNOWN')
      REWIND 1

    5 READ (1,9,END=99) KARTE
    9 FORMAT (A)

C***  LEVEL option, modernized by wrh 29-Jul-2007 

      IF (KARTE(:5) .EQ. 'LEVEL') THEN
C                         =====
         NLEVELCARD = NLEVELCARD + 1
         IF (NLEVELCARD .GT. MAXLEVELCARD) GOTO 93
         LEVELCARD(NLEVELCARD) = KARTE

      ELSE IF (KARTE(:6) .EQ. 'POPMIN' ) THEN
C                              ======
         CALL SARGV (KARTE, 2, ACTPAR1)
         READ (ACTPAR1,'(F10.0)', ERR=95) POPMIN
 
      ELSE IF (KARTE(:8) .EQ. 'OLDSTART') THEN
C                              ========
         OLDSTART=.TRUE.
         CALL SARGC (KARTE, NPAR)
         DO I=2, NPAR         
            CALL SARGV (KARTE, I, ACTPAR1)
            IF (ACTPAR1(:6) .EQ. 'DEPART') THEN
C                                 ======
               BDEPART = .TRUE.
            ELSE IF (ACTPAR1(:3) .EQ. 'TAU') THEN
C                                      ===
               BTAUR   = .TRUE.
            ENDIF
         ENDDO
      ENDIF

      GOTO 5

C***  End of CARDS reached
   99 CLOSE (1)

      RETURN

C***  ERROR EXITS *******************************

   93 WRITE (0,*) '*** FATAL ERROR:' 
      WRITE (0,*) '*** More LEVEL cards found than dimensioned'
      WRITE (0,*) '*** MAXLEVELCARD = ', MAXLEVELCARD
      GOTO 96

   95 WRITE (0,*) '*** FATAL ERROR:' 
      WRITE (0,*) '*** Parameter not readable as floating-point number'
      GOTO 96

   96 WRITE (0,*) 'THE ERROR OCCURED IN THE FOLLOWING LINE:'
      WRITE (0,*) KARTE(:IDX(KARTE)) 
      STOP 'FATAL ERROR IN ADAPTER:DECADP'

      END 
