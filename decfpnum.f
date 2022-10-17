      SUBROUTINE DECFPNUM (LINE, IARG, VALUE)
C**********************************************************************
C***  Decodes a floting-point type argument from LINE
C***  ALLOWS THE USE OF THE "SAME" KEYWORD
C***  and the prefixes LOG and DEX
C***  Note that there exists a special SUBR. DECXYUNITS for X, Y 
C**********************************************************************
      CHARACTER LINE*(*), ACTPAR*40

      CALL SARGV (LINE, IARG, ACTPAR)
      LAST = IDX(ACTPAR)

      IF (ACTPAR(:4) .EQ. 'SAME') RETURN
      
      IF (ACTPAR(1:3) .EQ. 'LOG') THEN
         READ (ACTPAR(4: ), '(F10.0)') VALUE
         VALUE = ALOG10 (VALUE)
      ELSE IF (ACTPAR(1:3) .EQ. 'DEX') THEN
         READ (ACTPAR(4: ), '(F10.0)') VALUE
         VALUE = 10. ** (VALUE)
      ELSE IF (ACTPAR(1:1) .EQ. '(' .AND. 
     >        ACTPAR(LAST:LAST) .EQ. ')' ) THEN 
C***     Arithmetic expression
         CALL ARIPAR (ACTPAR, VALUE, IERR)
         IF (IERR .NE. 0) GOTO 99
      ELSE
         READ (ACTPAR, '(F20.0)', ERR=99) VALUE
      ENDIF

      RETURN

   99 PRINT *, '*** ERROR IN KASDEF OPTION: ***'
      PRINT *, '>>> ', ACTPAR(:IDX(ACTPAR)), 
     >         ' <<<< CANNOT BE DECODED AS FLOATING POINT NUMBER!'
      PRINT *, 'THE ERROR OCCURED IN THE FOLLOWING LINE:'
      PRINT *, LINE(:IDX(LINE))
      STOP '*** ERROR IN SUBROUTINE DECFPNUM ***'

      END

