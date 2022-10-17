C**********************************************************************
C***  This subroutine copies the specified SECOND MODEL to fort.4
C***     and checks that the DATOM files are identical 
C***  Called from: FORMAL
C**********************************************************************
      SUBROUTINE COPY_SECONDMODEL 
     >           (SECONDMODEL_PATH, IGNORE_DIFF, BIRONLINES) 

      CHARACTER SECONDMODEL_PATH*(*), TESTLINE*80
      LOGICAL IGNORE_DIFF, BIRONLINES, SECONDMODEL_EXISTS

C***  Test if SECONDMODEL_PATH exists, i.e. if output of ls not empty
      CALL SYSTEM  ('echo `ls -d ' // 
     >         SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) //
     >          ' `   > secmod_path ')

      OPEN (44, FILE='secmod_path', STATUS='OLD', ERR=990)
      READ (44, '(A)') TESTLINE
      IF (TESTLINE .NE. '' ) THEN
         WRITE (0,*)
     >    SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) // ' was found'
      ELSE
         WRITE (0,*) '*** ERROR: ' //
     >    SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) // ' not found!'
         GOTO 999
      ENDIF
      CLOSE (44)

ccc   Helge's suggestion with INQUIRE (does not work!??)
ccc      INQUIRE (FILE='SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH))', 
ccc     >         EXIST=SECONDMODEL_EXISTS)
ccc      IF (SECONDMODEL_EXISTS) THEN
ccc         WRITE (0,*)
ccc     >    SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) // ' was found'
ccc      ELSE
ccc         WRITE (0,*) '*** ERROR: ' //
ccc     >    SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) // ' not found!'
ccc         GOTO 999
ccc      ENDIF

C***  Note: mass-storage files cannot be opened by name; hence copy first
      CALL SYSTEM ('cp ' // SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) 
     >             // '/MODEL fort.4')
      CALL SYSTEM ('chmod u+w fort.4')

      WRITE (0,*) 'Copying secondmodel: ', 
     >             SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH))

C******************************************************************
C***  Check if DATOM files are identical:
      CALL SYSTEM ('diff ' // SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) 
     >             // '/DATOM DATOM | wc -c > DATOM.diff')
      OPEN (44, FILE='DATOM.diff', STATUS='OLD', ERR=991)
      READ (44, '(A)') TESTLINE
      CLOSE (44)
      IF (TESTLINE .NE. '0') THEN
         WRITE (0,*) 
     >    '*** ERROR: DATOM files of first and second model differ'
         IF (IGNORE_DIFF) THEN
            WRITE (0,*) 
     >      '*** NO ABORT since option IGNORE_DIFF has been requested'
         ELSE
            WRITE (0,'(A,/,A)') 
     >         '*** If you are sure that this difference is ' 
     >         // 'not significant:', '*** use the option: IGNORE_DIFF'
            GOTO 999
         ENDIF
      ELSE
         WRITE (0,*) 'DATOM files of first and second model agree'
      ENDIF
C*******************************************************************

C******************************************************************
C***  Check if FEDAT_FORMAL files differ:
      IF (BIRONLINES) THEN
       CALL SYSTEM ('diff ' // SECONDMODEL_PATH(:IDX(SECONDMODEL_PATH)) 
     >             // '/FEDAT_FORMAL fort.21 | wc -c > FEDAT.diff')
       OPEN (45, FILE='FEDAT.diff', STATUS='OLD', ERR=992)
       READ (45, '(A)') TESTLINE
       CLOSE (45)
       IF (TESTLINE .NE. '0') THEN
         WRITE (0,*) '*** IMPORTANT WARNING:  ********************'
         WRITE (0,*) 
     >     '*** ERROR: FEDAT files of first and second model differ'
         WRITE (0,*) '*** Make sure that both FEDAT files have same ' //
     >                    'superlevel structure!'
         WRITE (0,*) '*** In particular, check for parity splitting!'
         WRITE (0,*) '********************************************'
         IF (IGNORE_DIFF) THEN
            WRITE (0,*) 
     >      '*** NO ABORT since option IGNORE_DIFF has been requested'
         ELSE
            WRITE (0,*) 
            WRITE (0,'(A,/,A)') 
     >         '*** If you are sure that this difference is ' 
     >         // 'not significant:', '*** use the option: IGNORE_DIFF'
            GOTO 999
         ENDIF
       ELSE
         WRITE (0,*) 'FEDAT_FORMAL files of first and second model agree'
       ENDIF
      ENDIF
C*******************************************************************

  100 CONTINUE
      RETURN

C**********************************************************************
C***  ERROR BRANCHES 
C**********************************************************************

  990 WRITE (0,*) '*** Internal ERROR when opening secmod_exist'
      GOTO 999

  991 WRITE (0,*) '*** Internal ERROR when opening DATOM.diff'
      GOTO 999

  992 WRITE (0,*) '*** Internal ERROR when opening FEDAT.diff'
      GOTO 999

  999 STOP '*** FATAL ERROR detected by subr. COPY_SECONDMODEL' 

      END
