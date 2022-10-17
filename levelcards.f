      SUBROUTINE LEVELCARDS
C***************************************************************
C***  creates LEVEL Cards from comparing DATOM with DATOM_OLD
C***************************************************************

      PARAMETER (NDIM = 1560)
      CHARACTER*10 OLDLEVEL(NDIM), NEWLEVEL(NDIM)
      CHARACTER*80 LINE, OUTLINE(NDIM)
      CHARACTER*20 ACTPAR
      LOGICAL FEDATNEW, FEDATOLD

C***  number of superlevels in the standard FEDAT file with parity
C***  ~wrh/Blanket/Models/G2_BIG_VD100_FeXVI-K2015-parity
      DIMENSION NSUPER(16)
      DATA NSUPER /3, 3, 13, 18, 22, 29, 19, 14, 15, 28, 26, 13, 
     >             15, 14, 10, 9/
 
C***  Iron included?
      FEDATNEW = .FALSE. 
      FEDATOLD = .FALSE. 

      ISHIFT = 0
C***  read new DATOM
      OPEN (10, FILE='DATOM', STATUS='OLD', ERR=90) 
      I = 0
   10 READ (10,'(A80)', END=11) LINE
      IF (LINE(:5) .EQ. 'LEVEL') THEN
         I = I + 1
         IF (I .GT. NDIM) GOTO 95
         NEWLEVEL(I) = LINE(13:22)
      ELSEIF (LINE(:19) .EQ. 'ELEMENT     GENERIC') THEN
         FEDATNEW = .TRUE.
         CALL SARGV (LINE(31:), 1, ACTPAR)
         READ(ACTPAR,'(I20)') IONLOWNEW
      ENDIF
      GOTO 10
   11 CLOSE (10)
      NNEW=I

C***  read old DATOM
      OPEN (11, FILE='DATOM_OLD', STATUS='OLD', ERR=91) 
      I = 0
   20 READ (11,'(A80)', END=21) LINE
      IF (LINE(:5) .EQ. 'LEVEL') THEN
         I = I + 1
         IF (I .GT. NDIM) GOTO 96
         OLDLEVEL(I) = LINE(13:22)
      ELSEIF (LINE(:19) .EQ. 'ELEMENT     GENERIC') THEN
         FEDATOLD = .TRUE.
         CALL SARGV (LINE(31:), 1, ACTPAR)
         READ(ACTPAR,'(I20)') IONLOWOLD
      ENDIF
      GOTO 20
   21 CLOSE (11)
      NOLD=I

C***  Compare each new LEVEL with *all* LEVELs in old DATOM
      IOUT=0
      DO INEW = 1, NNEW
        DO IOLD = 1, NOLD
          IF (OLDLEVEL(IOLD) .EQ. NEWLEVEL(INEW)) THEN
            ISHIFT = INEW - IOLD
C***        first output line!
            IF (IOUT .EQ. 0) THEN  
                IOUT = IOUT + 1
                IOLD1 = IOLD
                IOLD2 = IOLD
                WRITE (OUTLINE(IOUT), '(A,I4,A,I4,A,I4)') 
     >            'LEVEL ', IOLD1, ' - ', IOLD2, ' SHIFT ', ISHIFT
                IOLD1LAST = IOLD1
                IOLD2LAST = IOLD2
                ISHIFTLAST = ISHIFT
            ELSE
C***          Do last LEVEL-CARDS form a connected range?
              IF (ISHIFTLAST .EQ. ISHIFT .AND.
     >           IOLD .EQ. IOLD2LAST+1) THEN           
                 IOLD2 = IOLD 
                 WRITE (OUTLINE(IOUT), '(A,I4,A,I4,A,I4)') 
     >             'LEVEL ', IOLD1, ' - ', IOLD2, ' SHIFT ', ISHIFT
                 IOLD1LAST = IOLD1
                 IOLD2LAST = IOLD2

              ELSE
                IOUT = IOUT + 1
                IOLD1 = IOLD
                IOLD2 = IOLD
                WRITE (OUTLINE(IOUT), '(A,I4,A,I4,A,I4)') 
     >            'LEVEL ', IOLD1, ' - ', IOLD2, ' SHIFT ', ISHIFT
                IOLD1LAST = IOLD1
                IOLD2LAST = IOLD2
                ISHIFTLAST = ISHIFT
              ENDIF
            ENDIF
            EXIT
          ENDIF
        ENDDO
      ENDDO

      IF (FEDATNEW .AND. FEDATOLD) THEN
         IOUT = IOUT + 1
         WRITE (OUTLINE(IOUT), '(A)') 
     >     '* The following lines are for iron levels'
         IOUT = IOUT + 1
         WRITE (OUTLINE(IOUT), '(A)') 
     >     '* and can ONLY be used if both models (OLD and current)'
         IOUT = IOUT + 1
         WRITE (OUTLINE(IOUT), '(A)') 
     >     '* have been used FEDAT with parity splitting!'

         ISHIFT_BEFOREIRON = NNEW - NOLD
C***     Lowest old Fe level (-1) is now the start of all calculations
         IOLD2 = NOLD
         
C***     lowest Fe ion unchanged:
         IF (IONLOWOLD .EQ. IONLOWNEW) THEN
           IOUT = IOUT + 1
           IOLD1 = IOLD2 + 1
           ISHIFT = ISHIFT_BEFOREIRON
           WRITE (OUTLINE(IOUT), '(A,I4,A,I4,A,I4)') 
     >          'LEVEL ', IOLD1, ' - ', NDIM, ' SHIFT ', ISHIFT

C***     lowest Fe ions removed:
         ELSEIF (IONLOWOLD .LT. IONLOWNEW) THEN
           NLEVELREMOVED = 1 + SUM (NSUPER(IONLOWOLD+1:IONLOWNEW-1))
           IOLD1 = IOLD2 + 1 + NLEVELREMOVED
           ISHIFT = -NLEVELREMOVED + ISHIFT_BEFOREIRON
           IOUT = IOUT + 1
           WRITE (OUTLINE(IOUT), '(A,I4,A,I4,A,I5)')
     >          'LEVEL ', IOLD1, ' - ', IOLD1, ' SHIFT ', ISHIFT

           IOLD1 = IOLD1 + NSUPER(IONLOWNEW)
           ISHIFT = ISHIFT - NSUPER(IONLOWNEW) + 1
           IOUT = IOUT + 1
           WRITE (OUTLINE(IOUT), '(A,I4,A,I4,A,I5)')
     >          'LEVEL ', IOLD1, ' - ', NDIM, ' SHIFT ', ISHIFT

C***     lowest Fe ions added:
         ELSEIF (IONLOWOLD .GT. IONLOWNEW) THEN
           NLEVELADDED = 1 + SUM (NSUPER(IONLOWNEW+1:IONLOWOLD-1))
           IOLD1 = IOLD2 + 1 
           ISHIFT = NLEVELADDED + ISHIFT_BEFOREIRON
           IOUT = IOUT + 1
           WRITE (OUTLINE(IOUT), '(A,I4,A,I4,A,I5)')
     >          'LEVEL ', IOLD1, ' - ', IOLD1, ' SHIFT ', ISHIFT

           IOLD1 = IOLD1 + 1
           ISHIFT = ISHIFT + NSUPER(IONLOWOLD) - 1
           IOUT = IOUT + 1
           WRITE (OUTLINE(IOUT), '(A,I4,A,I4,A,I5)')
     >          'LEVEL ', IOLD1, ' - ', NDIM, ' SHIFT ', ISHIFT

         ENDIF
      ENDIF

C***  final result:

      OPEN (13, FILE='LEVELCARDS', STATUS='UNKNOWN') 
      DO I=1, IOUT
         WRITE (13,'(A)') OUTLINE(I)(:IDX(OUTLINE(I))) 
         WRITE (*,'(A)') OUTLINE(I)(:IDX(OUTLINE(I))) 
      ENDDO
      CLOSE (13)
              
      WRITE (*,'(A)') '* File LEVELCARDS written'
      STOP '* O.K.'

C***  Error messages *********************************
   90 WRITE (0,*) '*** ERROR when opening file DATOM'
      STOP '*** FATAL ERROR in program LEVELCARDS'

   91 WRITE (0,*) '*** ERROR when opening file DATOM_OLD'
      STOP '*** FATAL ERROR in program LEVELCARDS'

   95 WRITE (0,*) '*** ERROR: file DATOM has more levels than dimensioned'
      STOP '*** FATAL ERROR in program LEVELCARDS'

   96 WRITE (0,*) '*** ERROR: file DATOM_OLD ' // 
     >          'has more levels than dimensioned'
      STOP '*** FATAL ERROR in program LEVELCARDS'

      END

