      SUBROUTINE OPENMSR (KANAL,IADR,MAXADR,MODE)
C***********************************************************************
C***  SUBROUTINE TO OPAN A RANDOM MASS STORAGE FILE (NAME INDEX).
C***  ALL EXISTING RECORDS ARE READ (ONLY THE FIRST WORD) TO AVOID
C***  THE COS "DATASET SECURITY VIOLATION ERROR"
C***********************************************************************
 
      DIMENSION IADR(MAXADR)
 
      IF (MODE .NE. 1) THEN
         CALL REMARK ('WRONG MODE OF RANDOM MASS STORAGE FILE')
         STOP 'ERROR'
         ENDIF
 
      CALL OPENMS (KANAL,IADR,MAXADR,MODE, IERR)
 
      DO 1 I=1,MAXADR,2
      NAME=IADR(I)
      IF (NAME .EQ. 0) GOTO 2
      CALL READMS (KANAL,DUMMY,1,NAME, IERR)
    1 CONTINUE
 
    2 CONTINUE
      RETURN
      END
